"""Per-residue features for one antigen structure, shared by training and prediction.

Everything here is computed on the antigen alone, with the antibody removed.
Leaving the antibody in makes contact residues look buried, which leaks the label
straight into the input.

SASA is computed twice: once on the isolated chain and once on all antigen
chains together. The difference exposes residues buried by the trimer interface,
which look bindable on a monomer but are inaccessible on a real virion.

scripts/02 builds the training table with these functions and epitope/predict.py
builds the same columns for an uploaded structure. Keeping one copy is the point:
a second one would drift, and the model would be scoring features it never saw.
"""

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended as THREE_TO_ONE
from Bio.PDB.SASA import ShrakeRupley

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

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# Columns the model reads, in the order it was trained on. chain_length is
# deliberately absent: it is constant within a chain, so the model would use it
# as a construct-size shortcut. protrusion is z-scored per chain for the same
# reason — raw distance from the centroid scales with how much of the trimer a
# structure happens to contain.
STRUCTURAL = [
    "rel_sasa", "rel_sasa_assembly", "buried_by_assembly", "protrusion_z",
    "glycan_distance", "hydrophobicity", "charge", "volume", "bfactor_z",
    "neighbours_8", "neighbours_12",
    "nbhd_rel_sasa", "nbhd_hydrophobicity", "nbhd_charge",
    "nbhd_protrusion_z", "nbhd_glycan_distance",
]

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]


def is_amino_acid(residue):
    return residue.get_resname().strip().upper() in THREE_TO_ONE


def amino_acids(chain):
    return [r for r in chain if is_amino_acid(r)]


def one_letter(residue):
    return THREE_TO_ONE.get(residue.get_resname().strip().upper(), "X")


def chain_sequence(chain):
    """One-letter sequence of a chain's amino acids, in the order features use."""
    return "".join(one_letter(r) for r in amino_acids(chain))


def sidechain_centre(residue):
    """CB where it exists, CA otherwise — glycine and incomplete side chains."""
    for name in ("CB", "CA"):
        if name in residue:
            return residue[name].coord
    return np.mean([a.coord for a in residue], axis=0)


def structural_features(model, antigen_chain_ids, label):
    """Geometry features for every residue of the given antigen chains.

    `model` is consumed: every other chain and every non-amino-acid residue is
    detached from it before the surface is measured, which is what removes the
    antibody. Returns one row per residue, keyed like the labels table, with
    `label` in the PDB column.
    """
    antigen_chain_ids = set(antigen_chain_ids)

    glycan_coords = np.asarray(
        [a.coord for r in model.get_residues()
         if r.get_resname().strip().upper() in GLYCANS for a in r],
        dtype=np.float32,
    )

    for chain in list(model):
        if chain.id not in antigen_chain_ids:
            model.detach_child(chain.id)
            continue
        for residue in list(chain):
            if not is_amino_acid(residue):
                chain.detach_child(residue.id)

    present = [c for c in model if amino_acids(c)]
    if not present:
        raise ValueError("no antigen amino acids")

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
                "PDB": label,
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
                "seq_index": index,
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

    features = pd.DataFrame(rows)
    features["buried_by_assembly"] = features["rel_sasa"] - features["rel_sasa_assembly"]
    return features


def add_derived_features(frame, chain_column="antigen_id"):
    """Per-chain z-scores and residue one-hots the model reads.

    B-factor scale differs between X-ray and cryo-EM, and protrusion scales with
    how much of the trimer was solved, so both are compared within a chain.
    """
    frame = frame.copy()
    for column in ("bfactor", "protrusion", "nbhd_protrusion"):
        grouped = frame.groupby(chain_column)[column]
        frame[f"{column}_z"] = (frame[column] - grouped.transform("mean")) / (
            grouped.transform("std") + 1e-6
        )

    letters = frame["residue"].where(frame["residue"].isin(AMINO_ACIDS), "X")
    for amino_acid in AMINO_ACIDS:
        frame[f"aa_{amino_acid}"] = (letters == amino_acid).astype(float)
    return frame
