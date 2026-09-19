"""Label antigen residues as epitope / non-epitope from antibody contact distances."""

import gzip
import pathlib
import warnings

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended as THREE_TO_ONE
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

DATASET_PATH = "flu_dataset.csv"
STRUCTURE_DIR = pathlib.Path("structures")
LABELS_PATH = "epitope_labels.csv"
FASTA_PATH = "antigen_sequences.fasta"

CONTACT_ANGSTROMS = 4.5
# distances are computed in blocks so the antigen x antibody matrix stays small
CHUNK_SIZE = 512

warnings.simplefilter("ignore", PDBConstructionWarning)
parser = MMCIFParser(QUIET=True)


def to_classic_id(extended_id):
    return str(extended_id).split("_")[-1][-4:].lower()


def load_model(pdb_id):
    path = STRUCTURE_DIR / f"{pdb_id}.cif.gz"
    with gzip.open(path, "rt") as handle:
        structure = parser.get_structure(pdb_id, handle)
    return next(structure.get_models())


def amino_acids(chain):
    """Standard and modified amino acids, skipping waters, glycans, ions and ligands."""
    return [r for r in chain if r.get_resname().strip().upper() in THREE_TO_ONE]


def heavy_atom_coords(residues):
    coords, owners = [], []
    for index, residue in enumerate(residues):
        for atom in residue:
            if atom.element != "H":
                coords.append(atom.coord)
                owners.append(index)
    return np.asarray(coords, dtype=np.float32), np.asarray(owners, dtype=np.int32)


def min_distance_per_residue(antigen_coords, antigen_owners, antibody_coords, residue_count):
    """Smallest distance from each antigen residue to any antibody heavy atom."""
    per_atom = np.full(len(antigen_coords), np.inf, dtype=np.float32)
    for start in range(0, len(antigen_coords), CHUNK_SIZE):
        block = antigen_coords[start:start + CHUNK_SIZE]
        distances = np.linalg.norm(block[:, None, :] - antibody_coords[None, :, :], axis=-1)
        per_atom[start:start + CHUNK_SIZE] = distances.min(axis=1)

    per_residue = np.full(residue_count, np.inf, dtype=np.float32)
    np.minimum.at(per_residue, antigen_owners, per_atom)
    return per_residue


df = pd.read_csv(DATASET_PATH)
df["antigen_chain_list"] = df["antigen_chain"].str.split(r"\s*\|\s*")
print(f"{len(df)} complexes across {df['PDB'].nunique()} structures")

model_cache = {}
label_rows = []
sequences = {}
skipped = []

for position, row in enumerate(df.itertuples(index=False), start=1):
    pdb_id = to_classic_id(row.PDB)
    try:
        if pdb_id not in model_cache:
            model_cache[pdb_id] = load_model(pdb_id)
        model = model_cache[pdb_id]
    except Exception as exc:
        skipped.append((row.PDB, "*", f"parse failed: {exc}"))
        continue

    antibody_chain_ids = [c for c in (row.Hchain, row.Lchain) if isinstance(c, str) and c.strip()]
    antibody_residues = []
    for chain_id in antibody_chain_ids:
        if chain_id in model:
            antibody_residues.extend(amino_acids(model[chain_id]))
    if not antibody_residues:
        skipped.append((row.PDB, ",".join(antibody_chain_ids), "antibody chains absent"))
        continue
    antibody_coords, _ = heavy_atom_coords(antibody_residues)

    for antigen_chain_id in row.antigen_chain_list:
        antigen_chain_id = antigen_chain_id.strip()
        if antigen_chain_id not in model:
            skipped.append((row.PDB, antigen_chain_id, "antigen chain absent"))
            continue

        antigen_residues = amino_acids(model[antigen_chain_id])
        if not antigen_residues:
            skipped.append((row.PDB, antigen_chain_id, "no amino acids (likely glycan or ion)"))
            continue

        antigen_coords, owners = heavy_atom_coords(antigen_residues)
        distances = min_distance_per_residue(
            antigen_coords, owners, antibody_coords, len(antigen_residues)
        )

        sequence = "".join(
            THREE_TO_ONE.get(r.get_resname().strip().upper(), "X") for r in antigen_residues
        )
        sequences[(pdb_id, antigen_chain_id)] = sequence

        for index, (residue, distance) in enumerate(zip(antigen_residues, distances)):
            hetflag, resseq, icode = residue.id
            label_rows.append({
                "PDB": pdb_id,
                "instance": row.INSTANCE,
                "Hchain": row.Hchain,
                "Lchain": row.Lchain,
                "antigen_chain": antigen_chain_id,
                "target": row.target,
                "ha_subtype": row.ha_subtype,
                "year": row.year,
                "resolution": row.resolution,
                "seq_index": index,
                "residue_number": resseq,
                "insertion_code": icode.strip(),
                "residue": THREE_TO_ONE.get(residue.get_resname().strip().upper(), "X"),
                "min_distance": round(float(distance), 3) if np.isfinite(distance) else None,
                "is_epitope": int(distance <= CONTACT_ANGSTROMS),
            })

    if position % 100 == 0 or position == len(df):
        print(f"  {position}/{len(df)} complexes processed")

labels = pd.DataFrame(label_rows)
labels.to_csv(LABELS_PATH, index=False)

with open(FASTA_PATH, "w") as handle:
    for (pdb_id, chain_id), sequence in sorted(sequences.items()):
        handle.write(f">{pdb_id}_{chain_id}\n{sequence}\n")

print()
print(f"labelled residues: {len(labels):,}")
print(f"epitope residues:  {int(labels['is_epitope'].sum()):,} "
      f"({100 * labels['is_epitope'].mean():.1f}%)")
print(f"unique antigen chains: {len(sequences)}")

per_complex = labels.groupby(["instance", "antigen_chain"])["is_epitope"].sum()
print(f"epitope residues per complex: median {per_complex.median():.0f}, "
      f"mean {per_complex.mean():.1f}, max {per_complex.max()}")
print(f"complexes with zero contacts: {(per_complex == 0).sum()}")

print(f"\nwrote {LABELS_PATH} and {FASTA_PATH}")
if skipped:
    print(f"\nskipped {len(skipped)} chain entries; first 10:")
    for entry in skipped[:10]:
        print("  ", entry)
