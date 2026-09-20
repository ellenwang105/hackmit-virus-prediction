"""The HTTP layer: shape of the response, and every way a request can go wrong."""

import gzip
import pathlib
import sys

import numpy as np
import pytest
from fastapi.testclient import TestClient

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from server import app as server  # noqa: E402

pytestmark = pytest.mark.timeout(300)

client = TestClient(server.app)

# what web/src/types.ts declares for one chain
ARRAYS = ["num", "ic", "aa", "ha", "region", "site", "rbs", "fusion", "sasa", "glycan",
          "epitope", "nAb", "nContact", "score", "constraint", "durability"]


def cif_gz(pdb):
    return (ROOT / "structures" / f"{pdb}.cif.gz").read_bytes()


def upload(pdb, name=None):
    return client.post("/api/predict", files={"file": (name or f"{pdb}.cif.gz", cif_gz(pdb))})


@pytest.fixture(scope="module")
def result():
    reply = upload("5k9k")
    assert reply.status_code == 200, reply.text
    return reply.json()


def test_response_carries_one_antigen_per_hemagglutinin_chain(result):
    assert {a["chain"] for a in result["antigens"]} == {"F", "I"}
    assert {c["chain"] for c in result["chains"] if not c["isAntigen"]} >= {"H", "L"}
    assert result["applicability"]["level"] == "held_out"


def test_every_array_lines_up_with_the_residues(result):
    for antigen in result["antigens"]:
        n = len(antigen["num"])
        assert n > 100
        for key in ARRAYS:
            assert len(antigen[key]) == n, key
        assert antigen["split"] == "upload"
        assert all(0.0 <= s <= 1.0 for s in antigen["score"])


def test_json_has_no_nan(result):
    # NaN is not valid JSON; the parse above would have failed, but be explicit
    for antigen in result["antigens"]:
        for key in ("score", "constraint", "durability", "sasa", "glycan"):
            assert not any(isinstance(v, float) and np.isnan(v) for v in antigen[key])


def test_durability_is_present_for_a_standard_h3(result):
    antigen = next(a for a in result["antigens"] if a["chain"] == "F")
    scored = [d for d in antigen["durability"] if d is not None]
    assert len(scored) > 0.8 * len(antigen["durability"])


def test_the_viewer_can_fetch_the_coordinates_back(result):
    reply = client.get(result["structureUrl"])
    assert reply.status_code == 200
    assert reply.text.lstrip().startswith("data_")
    assert result["antigens"][0]["source"] == {"url": result["structureUrl"], "format": "cif"}


def test_unknown_structure_token_is_404():
    assert client.get("/api/structures/deadbeef").status_code == 404


def test_pdb_id_path_fetches_then_scores(monkeypatch):
    monkeypatch.setattr(server, "fetch_rcsb", lambda pdb_id: cif_gz(pdb_id.lower()))
    reply = client.post("/api/predict/pdb/5k9k")
    assert reply.status_code == 200
    assert reply.json()["label"] == "5K9K"
    assert {a["chain"] for a in reply.json()["antigens"]} == {"F", "I"}


@pytest.mark.parametrize("bad", ["abc", "12345", "ab!d", "5k9k-"])
def test_malformed_pdb_id_is_refused(bad):
    assert client.post(f"/api/predict/pdb/{bad}").status_code in (404, 422)


def test_missing_rcsb_entry_is_a_clear_404(monkeypatch):
    def gone(pdb_id):
        raise server.HTTPException(404, f"RCSB has no entry {pdb_id.upper()}.")
    monkeypatch.setattr(server, "fetch_rcsb", gone)
    reply = client.post("/api/predict/pdb/0000")
    assert reply.status_code == 404
    assert "no entry" in reply.json()["detail"]


def test_garbage_upload_is_a_422_with_a_reason():
    reply = client.post("/api/predict", files={"file": ("notes.txt", b"not a structure")})
    assert reply.status_code == 422
    assert reply.json()["detail"]


def test_oversized_upload_is_413(monkeypatch):
    monkeypatch.setattr(server, "MAX_UPLOAD_BYTES", 1000)
    reply = client.post("/api/predict", files={"file": ("big.cif", b"x" * 2000)})
    assert reply.status_code == 413


def test_only_the_last_few_structures_are_kept(monkeypatch):
    monkeypatch.setattr(server, "KEEP_STRUCTURES", 2)
    tokens = [server._remember(f"t{i}", "cif") for i in range(4)]
    assert tokens[0] not in server._structures and tokens[-1] in server._structures
