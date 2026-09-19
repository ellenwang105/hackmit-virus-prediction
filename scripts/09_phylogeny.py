"""Build a subtype-level phylogeny for the antigen picker.

The picker lists 684 chains, which is far too many to draw as a tree. Subtypes
are the useful grain: 11 tips, each standing for every chain of that subtype.

The tree is computed from real pairwise sequence identity rather than drawn by
hand, because its job is to show how far each subtype sits from the training
set. Group 1 is what the model learned; everything else is held out, and the
scores fall off with distance along these branches.
"""

import itertools
import json
import pathlib

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

PREDICTIONS_PATH = "results/epitope_predictions.csv"
FASTA_PATH = "antigen_sequences.fasta"
OUTPUT_PATH = pathlib.Path("web/public/data/phylogeny.json")

# subtypes with no coherent sequence of their own
EXCLUDED = {"unknown", "chimeric", ""}
GROUP_1 = {"H1", "H2", "H5", "H6", "H8", "H9", "H11", "H12", "H13", "H16", "H17", "H18"}


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            sequences[name] = []
        elif name:
            sequences[name].append(line)
    return {k: "".join(v) for k, v in sequences.items()}


def identity(first, second, aligner):
    """Percent identity from a global alignment, over the shorter sequence."""
    alignment = aligner.align(first, second)[0]
    matches = sum(1 for a, b in zip(alignment[0], alignment[1]) if a == b and a != "-")
    return matches / min(len(first), len(second))


def upgma(names, distances):
    """Average-linkage tree as nested dicts, each node carrying its height."""
    clusters = {
        i: {"name": names[i], "height": 0.0, "members": [i], "children": []}
        for i in range(len(names))
    }
    matrix = {
        (i, j): distances[i][j]
        for i in range(len(names)) for j in range(i + 1, len(names))
    }
    next_key = len(names)

    while len(clusters) > 1:
        (a, b), best = min(matrix.items(), key=lambda kv: kv[1])
        node = {
            "name": "",
            "height": best / 2,
            "members": clusters[a]["members"] + clusters[b]["members"],
            "children": [clusters[a], clusters[b]],
        }
        del clusters[a], clusters[b]
        for key in [k for k in matrix if a in k or b in k]:
            del matrix[key]
        for other in clusters:
            total = sum(
                distances[m][n]
                for m in node["members"] for n in clusters[other]["members"]
            )
            average = total / (len(node["members"]) * len(clusters[other]["members"]))
            matrix[(min(other, next_key), max(other, next_key))] = average
        clusters[next_key] = node
        next_key += 1

    return next(iter(clusters.values()))


def strip(node):
    """Drop the bookkeeping members list before writing the tree out."""
    return {
        "name": node["name"],
        "height": round(node["height"], 5),
        "children": [strip(c) for c in node["children"]],
    }


predictions = pd.read_csv(PREDICTIONS_PATH, low_memory=False)
predictions = predictions[predictions["ha_subtype"].notna()]
fasta = read_fasta(FASTA_PATH)

# one representative per subtype: the longest chain, which is the most complete
per_subtype = {}
for subtype, group in predictions.groupby("ha_subtype"):
    if subtype in EXCLUDED:
        continue
    lengths = {
        antigen_id: len(fasta[antigen_id])
        for antigen_id in group["antigen_id"].unique()
        if antigen_id in fasta
    }
    if not lengths:
        continue
    representative = max(lengths, key=lengths.get)
    per_subtype[subtype] = {
        "representative": representative,
        "sequence": fasta[representative],
        "chains": int(group["antigen_id"].nunique()),
        "structures": int(group["PDB"].nunique()),
        "residues": int(len(group)),
        "heldOut": bool(group["split_group"].iloc[0] in ("test_group2", "test_B")),
        "split": str(group["split_group"].mode().iloc[0]),
        "group": "group1" if subtype in GROUP_1 else ("B" if subtype == "B" else "group2"),
    }

names = sorted(per_subtype, key=lambda s: -per_subtype[s]["chains"])
print(f"{len(names)} subtypes: {', '.join(names)}")

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score = -11, -1
aligner.mode = "global"

size = len(names)
identities = np.eye(size)
for i, j in itertools.combinations(range(size), 2):
    value = identity(per_subtype[names[i]]["sequence"], per_subtype[names[j]]["sequence"], aligner)
    identities[i][j] = identities[j][i] = value
    print(f"  {names[i]:>3} vs {names[j]:<3} {value:5.1%}")

tree = upgma(names, (1 - identities).tolist())

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH.write_text(json.dumps({
    "tree": strip(tree),
    "subtypes": {
        name: {k: v for k, v in per_subtype[name].items() if k != "sequence"}
        for name in names
    },
    "identity": {
        "names": names,
        "matrix": (identities * 100).round(1).tolist(),
    },
}, indent=1))

print(f"\nwrote {OUTPUT_PATH}")
print(pd.DataFrame((identities * 100).round(0).astype(int), index=names, columns=names).to_string())
