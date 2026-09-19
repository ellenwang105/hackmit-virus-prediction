"""Build the constraint score C from HA surveillance sequences.

Pulls H1 and H3 HA protein sequences from the NCBI Influenza Virus Resource,
stratified by year, maps them onto the reference numbering frame, and writes one
row per (reference, piece, ha_number).

  data/alignments/{subtype}_{year_range}.fasta   cached download
  data/constraint_v1.parquet                     the score

Rerunning reuses the cached FASTA, so the download happens once.
"""

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from epitope import constraint as C  # noqa: E402
from epitope.annotations import (  # noqa: E402
    FUSION_PEPTIDE, LONG_HELIX, RECEPTOR_BINDING_SITE,
)

ALIGNMENT_DIR = pathlib.Path("data/alignments")
OUTPUT_PATH = pathlib.Path("data/constraint_v1.parquet")

YEARS = range(2010, 2026)
PER_YEAR = 300

# a column backed by too few sequences cannot support an entropy estimate
MIN_OBSERVATIONS = 50


def load_or_fetch(subtype):
    ALIGNMENT_DIR.mkdir(parents=True, exist_ok=True)
    path = ALIGNMENT_DIR / f"{subtype}_{YEARS[0]}_{YEARS[-1]}.fasta"

    if path.exists():
        sequences, header = [], None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                header = line
                sequences.append("")
            elif header is not None:
                sequences[-1] += line.strip()
        print(f"  {subtype}: {len(sequences):,} sequences from cache")
        return [s for s in sequences if s]

    print(f"  {subtype}: downloading {PER_YEAR}/year for {YEARS[0]}-{YEARS[-1]} ...")
    records = C.fetch_surveillance(subtype, YEARS, PER_YEAR)
    with path.open("w") as handle:
        for year, header, sequence in records:
            handle.write(f">{year}|{header}\n{sequence}\n")
    print(f"  {subtype}: {len(records):,} sequences downloaded -> {path}")
    return [sequence for _, _, sequence in records]


def build(subtype):
    reference, signal_length, ha2_start = C.fetch_reference(subtype)
    ha1_length = C.REFERENCES[subtype]["ha1_length"]
    print(f"  {subtype}: reference {C.REFERENCES[subtype]['accession']}, "
          f"{len(reference)} aa, HA1 1 at index {signal_length}, HA2 1 at index {ha2_start}")

    sequences = load_or_fetch(subtype)
    counts, stats = C.column_counts(sequences, reference)
    print(f"  {subtype}: {stats['direct']:,} direct, {stats['aligned']:,} aligned, "
          f"{stats['failed']:,} failed")

    entropies = [C.shannon_entropy(counter) for counter in counts]
    observations = [sum(counter.values()) for counter in counts]

    entropies = np.asarray(entropies, dtype=float)
    observations = np.asarray(observations)
    entropies[observations < MIN_OBSERVATIONS] = np.nan

    constraint = C.normalise_entropy(entropies)

    esm_component = None
    if C.esm_available():
        print(f"  {subtype}: ESM-2 masked marginals ...")
        esm_component = C.normalise_log_probs(C.esm_masked_marginal(reference))
    else:
        print(f"  {subtype}: ESM-2 unavailable (no torch/esm) -- entropy only")

    combined = C.combine(constraint, esm_component)

    mapping = C.position_map(signal_length, ha1_length, ha2_start, len(reference))
    rows = []
    for index, (piece, number) in sorted(mapping.items()):
        rows.append({
            "reference": subtype,
            "piece": piece,
            "ha_number": number,
            "ref_residue": reference[index],
            "n_observed": int(observations[index]),
            "entropy_bits": float(entropies[index]),
            "constraint_entropy": float(constraint[index]),
            "constraint_esm": float(esm_component[index]) if esm_component is not None else np.nan,
            "constraint_score": float(combined[index]),
        })

    frame = pd.DataFrame(rows)
    functional = [
        (piece == "HA1" and number in RECEPTOR_BINDING_SITE)
        or (piece == "HA2" and (number in FUSION_PEPTIDE or number in LONG_HELIX))
        for piece, number in zip(frame.piece, frame.ha_number)
    ]
    frame["is_functional_prior"] = functional
    frame["constraint_score_prior"] = C.apply_prior(frame.constraint_score, functional)
    return frame


def main():
    frames = []
    for subtype in ("H3", "H1"):
        print(f"{subtype}:")
        frames.append(build(subtype))
    frame = pd.concat(frames, ignore_index=True)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(OUTPUT_PATH, index=False)

    print(f"\nwrote {OUTPUT_PATH}  {len(frame):,} rows")
    for reference, group in frame.groupby("reference"):
        usable = group.constraint_score.notna().sum()
        print(f"  {reference}: {usable:,}/{len(group):,} positions scored, "
              f"median C = {group.constraint_score.median():.3f}")


if __name__ == "__main__":
    main()
