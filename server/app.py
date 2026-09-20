"""HTTP API for scoring a structure the model has never seen.

    .venv/bin/uvicorn server.app:app --port 8000

The scoring itself is epitope/predict.py; this only moves files in and results
out, in the column-oriented shape web/src/types.ts already defines for a chain,
so the viewer needs no special case for an uploaded structure beyond where its
coordinates come from.

Everything is held in memory. Nothing a user uploads is written to disk, and the
last few structures are kept only so the viewer can fetch the coordinates back.
"""

import datetime
import os
import re
import threading
import urllib.error
import urllib.request
import uuid
from collections import OrderedDict
from contextlib import asynccontextmanager

import numpy as np
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse

# imports xgboost before torch, which predict.py needs on macOS
from epitope import predict as P

MAX_UPLOAD_BYTES = 80 * 1024 * 1024
KEEP_STRUCTURES = 16
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

# One request at a time. The ESM model and the tree model share a process that is
# only stable single-threaded (see the OpenMP note in predict.py), and a second
# concurrent request would only slow the first down.
_lock = threading.Lock()
_structures = OrderedDict()   # token -> (text, format)


@asynccontextmanager
async def lifespan(app):
    # loading the language model takes seconds; do it before the first upload
    with _lock:
        P._booster(), P._projection(), P._references(), P._constraint_table(), P._esm()
    yield


app = FastAPI(title="Epitope Explorer", lifespan=lifespan)
# Where the web app is served from. Localhost is always allowed for development;
# a deployed site's origin goes in CORS_ORIGINS, comma separated, e.g.
# CORS_ORIGINS=https://epitope-explorer.vercel.app
_origins = [o.strip().rstrip("/") for o in os.environ.get("CORS_ORIGINS", "").split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_origins,
    allow_origin_regex=r"http://(localhost|127\.0\.0\.1)(:\d+)?",
    allow_methods=["*"],
    allow_headers=["*"],
)


def _remember(text, fmt):
    token = uuid.uuid4().hex
    _structures[token] = (text, fmt)
    while len(_structures) > KEEP_STRUCTURES:
        _structures.popitem(last=False)
    return token


def _nullable(values, digits):
    """Round floats, turning NaN into None so JSON gets null instead of NaN."""
    return [None if v is None or np.isnan(v) else round(float(v), digits) for v in values]


def _phylo(subtype):
    if subtype in P.GROUP_1:
        return "group1"
    if subtype in P.GROUP_2:
        return "group2"
    return "B" if subtype == "B" else "other"


def _antigen(prediction, report, token, label, fmt):
    rows = prediction.residues[prediction.residues["antigen_chain"] == report.chain]
    zeros = [0] * len(rows)
    return {
        "id": f"upload-{token[:8]}_{report.chain}",
        "pdb": label,
        "chain": report.chain,
        "subtype": report.closest_subtype or "HA",
        "year": 0,
        "phylo": _phylo(report.closest_subtype),
        "split": "upload",
        "chainType": rows["chain_type"].iloc[0],
        "num": rows["residue_number"].astype(int).tolist(),
        "ic": rows["insertion_code"].fillna("").astype(str).str.strip().tolist(),
        "aa": rows["residue"].tolist(),
        "ha": rows["ha_number"].astype(int).tolist(),
        "region": rows["region"].fillna("").tolist(),
        "site": rows["antigenic_site"].fillna("").tolist(),
        "rbs": rows["is_rbs"].astype(int).tolist(),
        "fusion": rows["is_fusion_machinery"].astype(int).tolist(),
        "sasa": _nullable(rows["rel_sasa_assembly"], 3),
        "glycan": _nullable(rows["glycan_distance"], 1),
        # an upload has no antibody, so there is nothing observed to compare with
        "epitope": zeros,
        "nAb": zeros,
        "nContact": zeros,
        "score": _nullable(rows["epitope_score"], 4),
        "constraint": _nullable(rows["constraint_score"], 4),
        "durability": _nullable(rows["durability_score"], 4),
        "source": {"url": f"/api/structures/{token}", "format": fmt},
    }


def _score(data, filename, label):
    if len(data) > MAX_UPLOAD_BYTES:
        raise HTTPException(413, f"File exceeds the {MAX_UPLOAD_BYTES // 1024 // 1024} MB limit. "
                                 "Submit a single hemagglutinin assembly.")
    try:
        text, fmt = P.decode_structure(data, filename)
        with _lock:
            prediction = P.predict(data, filename)
    except P.UnsupportedInput as exc:
        raise HTTPException(422, str(exc)) from exc

    token = _remember(text, fmt)
    return {
        "token": token,
        "label": label,
        "format": fmt,
        "structureUrl": f"/api/structures/{token}",
        "applicability": prediction.applicability,
        "warnings": prediction.warnings,
        "chains": [
            {"chain": c.chain, "length": c.length, "isAntigen": c.is_antigen,
             "haScore": round(c.ha_score), "subtype": c.closest_subtype,
             "identity": round(c.identity, 3)}
            for c in prediction.chains
        ],
        "antigens": [_antigen(prediction, c, token, label, fmt)
                     for c in prediction.chains if c.is_antigen],
    }


def _label_from(filename):
    stem = re.sub(r"\.(gz|cif|mmcif|pdb|ent)$", "", filename or "upload", flags=re.I)
    stem = re.sub(r"\.(cif|mmcif|pdb|ent)$", "", stem, flags=re.I)
    return re.sub(r"[^A-Za-z0-9_-]", "", stem)[:16] or "upload"


@app.get("/api/health")
def health():
    return {"ok": True}


@app.post("/api/predict")
def predict_upload(file: UploadFile = File(...)):
    data = file.file.read(MAX_UPLOAD_BYTES + 1)
    return _score(data, file.filename or "upload", _label_from(file.filename))


def fetch_rcsb(pdb_id):
    """Bytes of an RCSB entry. Separate so tests can stand in for the network."""
    try:
        with urllib.request.urlopen(RCSB_URL.format(pdb_id=pdb_id.upper()), timeout=30) as reply:
            return reply.read(MAX_UPLOAD_BYTES + 1)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise HTTPException(404, f"No RCSB entry found for {pdb_id.upper()}.") from exc
        raise HTTPException(502, f"RCSB returned an error ({exc.code}). Retry, or upload the file directly.") from exc
    except (urllib.error.URLError, TimeoutError) as exc:
        raise HTTPException(502, "Unable to reach RCSB. Check the connection or upload the file directly.") from exc


@app.post("/api/predict/pdb/{pdb_id}")
def predict_pdb(pdb_id: str):
    if not PDB_ID.match(pdb_id):
        raise HTTPException(422, "A PDB ID is four characters: a digit followed by three letters or digits (e.g., 5K9K).")
    return _score(fetch_rcsb(pdb_id), f"{pdb_id}.cif.gz", pdb_id.upper())


@app.get("/api/structures/{token}", response_class=PlainTextResponse)
def structure(token: str):
    if token not in _structures:
        raise HTTPException(404, "This structure is no longer cached. Submit it again.")
    return _structures[token][0]
