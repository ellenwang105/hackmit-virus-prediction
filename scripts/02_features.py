"""Per-residue structural features for every antigen chain, antibody removed.

Features are computed on the antigen alone. Leaving the antibody in place makes
contact residues look buried, which leaks the label straight into the input.

SASA is computed twice: once on the isolated chain and once on all antigen
chains together. The difference exposes residues buried by the trimer interface,
which look bindable on a monomer but are inaccessible on a real virion.
"""

import concurrent.futures as futures
import gzip
import pathlib
import warnings

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended as THREE_TO_ONE
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.SASA import ShrakeRupley

RESIDUES_PATH = "data/residues_v1.parquet"
STRUCTURE_DIR = pathlib.Path("structures")
OUTPUT_PATH = pathlib.Path("data/features_v1.parquet")

MAX_WORKERS = 12
NEIGHBOUR_RADII = (8.0, 12.0)
AGGREGATION_RADIUS = 10.0

# sugar residues modelled in the structures; 254 of 290 carry them
GLYCANS = {"NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "SIA", "GLC", "XYS"}

# Tien et al. 2013 theoretical maximum accessible surface area
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.5}
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}

# aggregated over each residue's 10 A neighbourhood, giving the model a sense of
# the surface patch a residue sits in rather than the residue alone
AGGREGATED = ["rel_sasa", "hydrophobicity", "charge", "protrusion", "glycan_distance"]

warnings.simplefilter("ignore", PDBConstructionWarning)


def amino_acids(chain):
    return [r for r in chain if r.get_resname().strip().upper() in THREE_TO_ONE]


def one_letter(residue):
    return THREE_TO_ONE.get(residue.get_resname().strip().upper(), "X")


def sidechain_centre(residue):
    """CB where it exists, CA otherwise — glycine and incomplete side chains."""
    for name in ("CB", "CA"):
        if name in residue:
            return residue[name].coord
    return np.mean([a.coord for a in residue], axis=0)


def extract(job):
    pdb_id, antigen_chain_ids = job
    parser = MMCIFParser(QUIET=True)
    try:
        with gzip.open(STRUCTURE_DIR / f"{pdb_id}.cif.gz", "rt") as handle:
            model = next(parser.get_structure(pdb_id, handle).get_models())
    except Exception as exc:
        return pdb_id, None, f"parse failed: {exc}"

    glycan_coords = np.asarray(
        [a.coord for r in model.get_residues()
         if r.get_resname().strip().upper() in GLYCANS for a in r],
        dtype=np.float32,
    )

    # drop the antibody and every non-amino-acid residue before measuring surface
    for chain in list(model):
        if chain.id not in antigen_chain_ids:
            model.detach_child(chain.id)
            continue
        for residue in list(chain):
            if residue.get_resname().strip().upper() not in THREE_TO_ONE:
                chain.detach_child(residue.id)

    present = [c for c in model if amino_acids(c)]
    if not present:
        return pdb_id, None, "no antigen amino acids"

    sasa = ShrakeRupley(n_points=100)
    sasa.compute(model, level="R")
    assembly_sasa = {(c.id, r.id): r.sasa for c in model for r in amino_acids(c)}

    rows = []
    for chain in present:
        residues = amino_acids(chain)
        sasa.compute(chain, level="R")  # isolated-chain surface, overwrites r.sasa
        monomer_sasa = {r.id: r.sasa for r in residues}

        centres = np.asarray([sidechain_centre(r) for r in residues], dtype=np.float32)
        centroid = centres.mean(axis=0)
        distances = np.linalg.norm(centres[:, None, :] - centres[None, :, :], axis=-1)

        for index, residue in enumerate(residues):
            letter = one_letter(residue)
            max_asa = MAX_ASA.get(letter, 200.0)
            hetflag, resseq, icode = residue.id

            if len(glycan_coords):
                glycan_distance = float(
                    np.linalg.norm(glycan_coords - centres[index], axis=1).min()
                )
            else:
                glycan_distance = 50.0

            rows.append({
                "PDB": pdb_id,
                "antigen_chain": chain.id,
                "residue_number": resseq,
                "insertion_code": icode.strip(),
                "residue": letter,
                "rel_sasa": min(monomer_sasa[residue.id] / max_asa, 1.5),
                "rel_sasa_assembly": min(assembly_sasa[(chain.id, residue.id)] / max_asa, 1.5),
                "hydrophobicity": KYTE_DOOLITTLE.get(letter, 0.0),
                "charge": CHARGE.get(letter, 0.0),
                "volume": VOLUME.get(letter, 140.0),
                "protrusion": float(np.linalg.norm(centres[index] - centroid)),
                "glycan_distance": min(glycan_distance, 50.0),
                "bfactor": float(np.mean([a.bfactor for a in residue])),
                "neighbours_8": int((distances[index] < NEIGHBOUR_RADII[0]).sum() - 1),
                "neighbours_12": int((distances[index] < NEIGHBOUR_RADII[1]).sum() - 1),
                "chain_length": len(residues),
                "_index": index,
            })

        # neighbourhood aggregation, within the chain
        chain_rows = rows[-len(residues):]
        frame = pd.DataFrame(chain_rows)
        mask = distances < AGGREGATION_RADIUS
        for column in AGGREGATED:
            values = frame[column].to_numpy(dtype=np.float32)
            means = (mask * values[None, :]).sum(axis=1) / mask.sum(axis=1)
            for row, mean in zip(chain_rows, means):
                row[f"nbhd_{column}"] = float(mean)

    return pdb_id, rows, None


def main():
    residues = pd.read_parquet(RESIDUES_PATH)
    jobs = [
        (pdb_id, set(group["antigen_chain"]))
        for pdb_id, group in residues.groupby("PDB")
    ]
    print(f"{len(jobs)} structures, {residues['antigen_id'].nunique()} antigen chains")

    all_rows, failures = [], []
    with futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        for done, (pdb_id, rows, error) in enumerate(pool.map(extract, jobs), start=1):
            if error:
                failures.append((pdb_id, error))
            else:
                all_rows.extend(rows)
            if done % 50 == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} structures")

    features = pd.DataFrame(all_rows).drop(columns=["_index"])
    features["buried_by_assembly"] = features["rel_sasa"] - features["rel_sasa_assembly"]

    OUTPUT_PATH.parent.mkdir(exist_ok=True)
    features.to_parquet(OUTPUT_PATH, index=False)

    print()
    print(f"{len(features):,} residues with {features.shape[1]} columns")
    print(f"wrote {OUTPUT_PATH}")
    if failures:
        print(f"\n{len(failures)} structures failed:")
        for pdb_id, error in failures[:10]:
            print(f"  {pdb_id}: {error}")


if __name__ == "__main__":
    main()
