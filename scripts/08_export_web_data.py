"""Turn results/epitope_predictions.csv into static files the web app can load.

Writes, under web/public/:
  data/index.json            one summary row per antigen chain (drives the picker)
  data/metrics.json          model vs surface-exposure baseline per split
  data/antigens/<id>.json    per-residue arrays for one chain, column-oriented
  structures/<pdb>.cif.gz    the mmCIF for every antigen the app can show

Per-antigen files are column-oriented (one array per field) so a 300-residue
chain is a few KB. Re-run this after scripts/06 or once constraint_score is
filled in; the app treats null constraint as "not computed yet".
"""

import json
import pathlib
import shutil

import numpy as np
import pandas as pd

PREDICTIONS_CSV = pathlib.Path("results/epitope_predictions.csv")
METRICS_CSV = pathlib.Path("results/metrics_summary.csv")
STRUCTURE_DIR = pathlib.Path("structures")

WEB_PUBLIC = pathlib.Path("web/public")
DATA_DIR = WEB_PUBLIC / "data"
ANTIGEN_DIR = DATA_DIR / "antigens"
WEB_STRUCTURE_DIR = WEB_PUBLIC / "structures"


def average_precision(labels, scores):
    """AUPRC as the mean precision at each true positive, best score first."""
    order = np.argsort(-scores, kind="stable")
    hits = np.asarray(labels)[order] == 1
    precision = np.cumsum(hits) / np.arange(1, len(hits) + 1)
    return float(precision[hits].mean())


def nullable(series, digits):
    """Round floats, turning NaN into None so JSON gets null instead of NaN."""
    return [None if pd.isna(v) else round(float(v), digits) for v in series]


def antigen_payload(chain):
    first = chain.iloc[0]
    return {
        "id": first["antigen_id"],
        "pdb": first["PDB"],
        "chain": first["antigen_chain"],
        "subtype": first["ha_subtype"],
        "year": int(first["year"]),
        "phylo": first["phylo_group"],
        "split": first["split_group"],
        "chainType": first["chain_type"],
        # a residue is identified by number + insertion code, never number alone
        "num": chain["residue_number"].astype(int).tolist(),
        "ic": chain["insertion_code"].fillna("").astype(str).str.strip().tolist(),
        "aa": chain["residue"].tolist(),
        "ha": chain["ha_number"].astype(int).tolist(),
        "region": chain["region"].tolist(),
        "site": chain["antigenic_site"].fillna("").tolist(),
        "rbs": chain["is_rbs"].astype(int).tolist(),
        "fusion": chain["is_fusion_machinery"].astype(int).tolist(),
        "sasa": nullable(chain["rel_sasa_assembly"], 3),
        "glycan": nullable(chain["glycan_distance"], 1),
        "epitope": chain["is_epitope"].astype(int).tolist(),
        "nAb": chain["n_antibodies"].astype(int).tolist(),
        "nContact": chain["n_contacting"].astype(int).tolist(),
        "score": nullable(chain["epitope_score"], 4),
        "constraint": nullable(chain["constraint_score"], 4),
        "durability": nullable(chain["durability_score"], 4),
    }


def summary_row(payload):
    epitope = np.array(payload["epitope"])
    score = np.array(payload["score"], dtype=float)
    has_positives = 0 < epitope.sum() < len(epitope)
    return {
        "id": payload["id"],
        "pdb": payload["pdb"],
        "chain": payload["chain"],
        "subtype": payload["subtype"],
        "year": payload["year"],
        "phylo": payload["phylo"],
        "split": payload["split"],
        "chainType": payload["chainType"],
        "nResidues": len(epitope),
        "nEpitope": int(epitope.sum()),
        # AUPRC on this one chain; only meaningful when it has both classes
        "auprc": round(average_precision(epitope, score), 3)
        if has_positives else None,
    }


predictions = pd.read_csv(PREDICTIONS_CSV, low_memory=False)
ANTIGEN_DIR.mkdir(parents=True, exist_ok=True)
WEB_STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)

# clear stale chains so a shrunken export doesn't leave orphans behind
for stale in ANTIGEN_DIR.glob("*.json"):
    stale.unlink()

index = []
for antigen_id, chain in predictions.groupby("antigen_id", sort=True):
    chain = chain.sort_values(["residue_number", "insertion_code"], na_position="first")
    payload = antigen_payload(chain)
    (ANTIGEN_DIR / f"{antigen_id}.json").write_text(
        json.dumps(payload, separators=(",", ":")))
    index.append(summary_row(payload))

(DATA_DIR / "index.json").write_text(json.dumps(index, separators=(",", ":")))

metrics = pd.read_csv(METRICS_CSV)
(DATA_DIR / "metrics.json").write_text(
    json.dumps(metrics.replace({np.nan: None}).to_dict(orient="records"), indent=1))

copied = 0
for pdb in sorted(predictions["PDB"].unique()):
    source = STRUCTURE_DIR / f"{pdb}.cif.gz"
    target = WEB_STRUCTURE_DIR / source.name
    if not target.exists():
        shutil.copy2(source, target)
        copied += 1

data_mb = sum(p.stat().st_size for p in DATA_DIR.rglob("*.json")) / 1e6
structure_mb = sum(p.stat().st_size for p in WEB_STRUCTURE_DIR.glob("*.gz")) / 1e6
print(f"{len(index)} antigen chains -> {ANTIGEN_DIR} ({data_mb:.1f} MB of JSON)")
print(f"{predictions['PDB'].nunique()} structures in {WEB_STRUCTURE_DIR} "
      f"({copied} newly copied, {structure_mb:.0f} MB)")
