"""Download mmCIF structure files for every PDB entry in the flu dataset."""

import concurrent.futures as futures
import pathlib

import pandas as pd
import requests

DATASET_PATH = "flu_dataset.csv"
STRUCTURE_DIR = pathlib.Path("structures")
# gzipped mmCIF keeps ~300MB of cryo-EM structures down to a manageable size
URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
MAX_WORKERS = 8
TIMEOUT_SECONDS = 60


def to_classic_id(extended_id):
    """SAbDab stores extended ids like 'pdb_00004fqj'; RCSB downloads use '4fqj'."""
    return str(extended_id).split("_")[-1][-4:].lower()


def download(pdb_id):
    destination = STRUCTURE_DIR / f"{pdb_id}.cif.gz"
    if destination.exists() and destination.stat().st_size > 0:
        return pdb_id, "cached", None

    try:
        response = requests.get(URL_TEMPLATE.format(pdb_id=pdb_id.upper()), timeout=TIMEOUT_SECONDS)
        response.raise_for_status()
    except requests.RequestException as exc:
        return pdb_id, "failed", str(exc)

    destination.write_bytes(response.content)
    return pdb_id, "downloaded", None


df = pd.read_csv(DATASET_PATH)
pdb_ids = sorted({to_classic_id(v) for v in df["PDB"].unique()})
STRUCTURE_DIR.mkdir(exist_ok=True)
print(f"{len(pdb_ids)} unique structures to fetch into {STRUCTURE_DIR}/")

counts = {"downloaded": 0, "cached": 0, "failed": 0}
failures = []

with futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
    for done, (pdb_id, status, error) in enumerate(pool.map(download, pdb_ids), start=1):
        counts[status] += 1
        if status == "failed":
            failures.append((pdb_id, error))
        if done % 25 == 0 or done == len(pdb_ids):
            print(f"  {done}/{len(pdb_ids)}  downloaded={counts['downloaded']} "
                  f"cached={counts['cached']} failed={counts['failed']}")

total_mb = sum(p.stat().st_size for p in STRUCTURE_DIR.glob("*.cif.gz")) / 1024 / 1024
print(f"\ndownloaded={counts['downloaded']} cached={counts['cached']} failed={counts['failed']}")
print(f"{STRUCTURE_DIR}/ holds {total_mb:.1f} MB")

if failures:
    print("\nfailed ids:")
    for pdb_id, error in failures:
        print(f"  {pdb_id}: {error}")
