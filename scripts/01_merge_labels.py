"""Merge per-antibody epitope labels into one row per antigen residue, then assign splits.

A residue is positive if ANY antibody in ANY complex of that antigen contacts it.
Labelling per antibody instance leaves the same residue positive in one row and
negative in another; training on that teaches the model contradictions.

The same union is then taken a second time, across structures. 586 of 684 antigen
chains have exactly one antibody solved, so a per-structure negative mostly means
"the one antibody anybody crystallised here missed this spot", not "antibodies do
not bind here". Chains sharing an identical sequence are the same physical protein
seen with different antibodies, so their labels are pooled by position into
`is_epitope_union`. No sequence spans two splits or two subtypes, so pooling
cannot move a label across the split boundary.
"""

import difflib
import itertools
import pathlib

import pandas as pd

LABELS_PATH = "epitope_labels.csv"
FASTA_PATH = "antigen_sequences.fasta"
OUTPUT_PATH = pathlib.Path("data/residues_v1.parquet")

# sequence identity above which two chains count as the same antigen
IDENTITY_THRESHOLD = 0.95

# influenza A HA splits into two phylogenetic groups; this is the real
# generalisation axis, since all influenza A HA sits above the 30% identity
# threshold that sequence clustering would use.
GROUP_1 = {"H1", "H2", "H5", "H6", "H8", "H9", "H11", "H12", "H13", "H16", "H17", "H18"}
GROUP_2 = {"H3", "H4", "H7", "H10", "H14", "H15"}

# held out of training as the validation subtype
VAL_SUBTYPE = "H5"
TEMPORAL_CUTOFF = 2020


def phylo_group(subtype):
    if subtype in GROUP_1:
        return "group1"
    if subtype in GROUP_2:
        return "group2"
    if subtype == "B":
        return "B"
    return "other"


def assign_split(row):
    if row.phylo_group == "group1":
        return "val" if row.ha_subtype == VAL_SUBTYPE else "train"
    if row.phylo_group == "group2":
        return "test_group2"
    if row.phylo_group == "B":
        return "test_B"
    return "excluded"


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


def cluster_sequences(sequences, threshold=IDENTITY_THRESHOLD):
    """Single-linkage clusters of near-identical sequences.

    The 684 HA chains are only ~96 distinct antigens: the same protein is
    re-crystallised with antibody after antibody, and strains within a subtype
    differ by a handful of residues. Counting chains therefore overstates how
    much independent data there is by roughly 9x, so `seq_cluster` exists to
    weight or group by antigen rather than by deposition.

    Single linkage chains through intermediates, so an occasional cluster is
    too greedy — one merges H1 with influenza B via an unknown-subtype
    construct. Fine for weighting; do not use it to define splits.
    """
    ordered = sorted(sequences)
    parent = list(range(len(ordered)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(range(len(ordered)), 2):
        if find(a) == find(b):
            continue
        first, second = ordered[a], ordered[b]
        if abs(len(first) - len(second)) / max(len(first), len(second)) > 0.1:
            continue
        matcher = difflib.SequenceMatcher(None, first, second, autojunk=False)
        if matcher.quick_ratio() < threshold:
            continue
        if matcher.ratio() >= threshold:
            parent[find(a)] = find(b)

    return {sequence: find(i) for i, sequence in enumerate(ordered)}


def pool_across_structures(frame):
    """Union each position's label over every chain with an identical sequence.

    Keyed on (sequence, seq_index) rather than residue_number, because numbering
    is not consistent between constructs while seq_index indexes the FASTA that
    label_epitopes.py wrote. Identical sequence implies identical subtype, so
    this never pools across a split.
    """
    fasta = read_fasta(FASTA_PATH)
    frame["sequence"] = frame["antigen_id"].map(fasta)
    missing = frame["sequence"].isna().sum()
    if missing:
        print(f"WARNING: {missing:,} rows have no FASTA sequence and cannot be pooled")

    frame["seq_cluster"] = frame["sequence"].map(cluster_sequences(set(fasta.values())))

    key = ["sequence", "seq_index"]
    pooled = frame.dropna(subset=["sequence"]).groupby(key, sort=False).agg(
        is_epitope_union=("is_epitope", "max"),
        n_chains_covering=("antigen_id", "nunique"),
        n_chains_contacting=("is_epitope", "sum"),
    )
    frame = frame.merge(pooled, left_on=key, right_index=True, how="left")

    # Chains with no sequence keep their per-structure label rather than a NaN.
    frame["is_epitope_union"] = (
        frame["is_epitope_union"].fillna(frame["is_epitope"]).astype(int)
    )
    frame["n_chains_covering"] = frame["n_chains_covering"].fillna(1).astype(int)
    frame["n_chains_contacting"] = frame["n_chains_contacting"].fillna(
        frame["is_epitope"]).astype(int)
    return frame.drop(columns=["sequence"])


# keep_default_na=False matters: the target column uses "NA" for neuraminidase,
# which pandas otherwise reads as a missing value and blanks out 184 complexes.
labels = pd.read_csv(LABELS_PATH, keep_default_na=False)
labels["target"] = labels["target"].replace("", "NA")
labels["insertion_code"] = labels["insertion_code"].fillna("").astype(str)
print(f"{len(labels):,} per-antibody rows across {labels['PDB'].nunique()} structures")

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]

merged = labels.groupby(KEY, as_index=False).agg(
    target=("target", "first"),
    ha_subtype=("ha_subtype", "first"),
    year=("year", "first"),
    resolution=("resolution", "first"),
    seq_index=("seq_index", "first"),
    residue=("residue", "first"),
    min_distance=("min_distance", "min"),
    is_epitope=("is_epitope", "max"),          # union across antibodies
    n_antibodies=("is_epitope", "size"),       # how many antibodies saw this residue
    n_contacting=("is_epitope", "sum"),        # how many actually bound it
)

merged["phylo_group"] = merged["ha_subtype"].map(phylo_group)
merged["split_group"] = merged.apply(assign_split, axis=1)
merged["split_temporal"] = merged["year"].map(
    lambda y: "train" if int(y) <= TEMPORAL_CUTOFF else "test"
)
merged["antigen_id"] = merged["PDB"] + "_" + merged["antigen_chain"]

merged = pool_across_structures(merged)

OUTPUT_PATH.parent.mkdir(exist_ok=True)
merged.to_parquet(OUTPUT_PATH, index=False)

print(f"{len(merged):,} unique antigen residues")
print(f"positives: {int(merged['is_epitope'].sum()):,} ({100 * merged['is_epitope'].mean():.1f}%)")
print(f"antigen chains: {merged['antigen_id'].nunique()}")
print()

print("--- cross-structure pooling ---")
flipped = int((merged["is_epitope_union"] > merged["is_epitope"]).sum())
shared = merged["n_chains_covering"] > 1
print(f"positions pooled with >=2 chains: {int(shared.sum()):,} of {len(merged):,}")
print(f"  chains per position: median {merged['n_chains_covering'].median():.0f}, "
      f"max {merged['n_chains_covering'].max()}")
print(f"negatives flipped to positive:    {flipped:,}")
print(f"positive rate {merged['is_epitope'].mean():.4f} -> "
      f"{merged['is_epitope_union'].mean():.4f}")
print()

print(f"--- sequence clusters at {IDENTITY_THRESHOLD:.0%} identity ---")
ha_only = merged[merged["target"] == "HA"]
per_split = ha_only.groupby("split_group").agg(
    chains=("antigen_id", "nunique"), clusters=("seq_cluster", "nunique"))
per_split["chains_per_cluster"] = (
    per_split["chains"] / per_split["clusters"]).round(1)
print(per_split.reindex(["train", "val", "test_group2", "test_B"]).to_string())
print()

print("--- target ---")
print(merged.groupby("target").agg(
    chains=("antigen_id", "nunique"), residues=("is_epitope", "size"),
    positive_rate=("is_epitope", "mean"),
).round(3))
print()

print("--- phylogenetic split (HA only) ---")
ha = merged[merged["target"] == "HA"]
print(ha.groupby("split_group").agg(
    chains=("antigen_id", "nunique"), structures=("PDB", "nunique"),
    residues=("is_epitope", "size"), positives=("is_epitope", "sum"),
    positive_rate=("is_epitope", "mean"),
).round(3))
print()

print("--- temporal split (HA only) ---")
print(ha.groupby("split_temporal").agg(
    structures=("PDB", "nunique"), residues=("is_epitope", "size"),
    positive_rate=("is_epitope", "mean"),
).round(3))

print(f"\nwrote {OUTPUT_PATH}")
