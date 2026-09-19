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
import sys
import warnings

import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from epitope.features import structural_features  # noqa: E402

RESIDUES_PATH = "data/residues_v1.parquet"
STRUCTURE_DIR = pathlib.Path("structures")
OUTPUT_PATH = pathlib.Path("data/features_v1.parquet")

MAX_WORKERS = 12

warnings.simplefilter("ignore", PDBConstructionWarning)


def extract(job):
    pdb_id, antigen_chain_ids = job
    parser = MMCIFParser(QUIET=True)
    try:
        with gzip.open(STRUCTURE_DIR / f"{pdb_id}.cif.gz", "rt") as handle:
            model = next(parser.get_structure(pdb_id, handle).get_models())
    except Exception as exc:
        return pdb_id, None, f"parse failed: {exc}"

    try:
        features = structural_features(model, antigen_chain_ids, label=pdb_id)
    except ValueError as exc:
        return pdb_id, None, str(exc)
    return pdb_id, features, None


def main():
    residues = pd.read_parquet(RESIDUES_PATH)
    jobs = [
        (pdb_id, set(group["antigen_chain"]))
        for pdb_id, group in residues.groupby("PDB")
    ]
    print(f"{len(jobs)} structures, {residues['antigen_id'].nunique()} antigen chains")

    frames, failures = [], []
    with futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        for done, (pdb_id, features, error) in enumerate(pool.map(extract, jobs), start=1):
            if error:
                failures.append((pdb_id, error))
            else:
                frames.append(features)
            if done % 50 == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} structures")

    features = pd.concat(frames, ignore_index=True).drop(columns=["seq_index"])

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
