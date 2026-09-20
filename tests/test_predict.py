"""Does scoring an uploaded file reproduce what training produced?

The one that matters is the round trip: score a structure from the dataset as if
it were a fresh upload, and compare with the scores stored in results/. If they
disagree, feature construction at prediction time has drifted from feature
construction at training time, and the model is being fed inputs it never saw.

Run from the repo root:  .venv/bin/python -m pytest tests -q
"""

import gzip
import io
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import PDBIO, Select
from scipy.stats import spearmanr

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from epitope import predict as P  # noqa: E402

pytestmark = pytest.mark.timeout(300)

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]

# Structures whose antigen chains SAbDab labelled exactly as the detector finds
# them, so training and prediction saw the same assembly. Spread over every split.
ROUND_TRIP = {
    "5k9k": "H3, held out",
    "3ztj": "H3, held out",
    "6wf0": "H3, held out",
    "9jno": "influenza B, held out",
    "3gbm": "H5, validation",
    "6wex": "H6, training",
}


def structure_bytes(pdb):
    return gzip.open(ROOT / "structures" / f"{pdb}.cif.gz").read()


@pytest.fixture(scope="module")
def stored():
    table = pd.read_csv(
        ROOT / "results" / "epitope_predictions.csv", low_memory=False,
        usecols=KEY + ["epitope_score"],
    )
    table["insertion_code"] = table["insertion_code"].fillna("").astype(str)
    return table


def compare(prediction, pdb, stored):
    new = prediction.residues.assign(PDB=pdb)
    merged = new.merge(stored, on=KEY, suffixes=("_new", "_old"))
    return merged, np.abs(merged["epitope_score_new"] - merged["epitope_score_old"])


@pytest.mark.parametrize("pdb", ROUND_TRIP, ids=[f"{k} ({v})" for k, v in ROUND_TRIP.items()])
def test_round_trip_reproduces_stored_scores(pdb, stored):
    prediction = P.predict(structure_bytes(pdb), f"{pdb}.cif")
    merged, diff = compare(prediction, pdb, stored)

    expected = stored[stored["PDB"] == pdb]
    assert len(merged) == len(expected), "every stored residue should be scored"
    assert set(prediction.antigen_chains) == set(expected["antigen_chain"]), (
        "detector and SAbDab should agree on which chains are antigen for these structures")

    # ranking is what the product shows, and it should be identical
    assert spearmanr(merged["epitope_score_new"], merged["epitope_score_old"])[0] > 0.999
    # values differ only by float noise between two runs of ESM-2
    assert diff.mean() < 1e-3
    assert diff.max() < 0.06


def test_split_chain_file_still_tracks_stored_scores(stored):
    """7K39 lists A, C, E as antigen but also carries B, D, F, which are HA too.

    Training stripped those, so its assembly surface was slightly wrong; scoring
    keeps them. The values move, but the ranking should still follow.
    """
    prediction = P.predict(structure_bytes("7k39"), "7k39.cif")
    assert set(prediction.antigen_chains) == set("ABCDEF")
    merged, _ = compare(prediction, "7k39", stored)
    assert len(merged) > 900
    assert spearmanr(merged["epitope_score_new"], merged["epitope_score_old"])[0] > 0.95


def test_antibody_chains_are_set_aside():
    prediction = P.predict(structure_bytes("5k9k"), "5k9k.cif")
    assert set(prediction.antigen_chains) == {"F", "I"}
    assert set(prediction.residues["antigen_chain"]) == {"F", "I"}
    assert any("set aside" in note for note in prediction.warnings)


def test_pdb_format_gives_the_same_scores():
    """The whole file, not just the HA chains: glycans live in chains of their own
    in 5K9K, and glycan distance is one of the model's features."""
    model = P.parse_structure(structure_bytes("5k9k"), "5k9k.cif")

    buffer = io.StringIO()
    io_ = PDBIO()
    io_.set_structure(model)
    io_.save(buffer)

    from_pdb = P.predict(buffer.getvalue().encode(), "5k9k.pdb")
    from_cif = P.predict(structure_bytes("5k9k"), "5k9k.cif")
    a = from_pdb.residues.set_index(["antigen_chain", "residue_number", "insertion_code"])["epitope_score"]
    b = from_cif.residues.set_index(["antigen_chain", "residue_number", "insertion_code"])["epitope_score"]
    b = b.loc[a.index]
    assert spearmanr(a, b)[0] > 0.999


@pytest.mark.parametrize("pdb,level", [("6urm", "near"), ("5k9k", "held_out"), ("9jno", "far")])
def test_applicability_follows_distance_from_training(pdb, level):
    result = P.predict(structure_bytes(pdb), f"{pdb}.cif").applicability
    assert result["level"] == level
    assert result["expected_auprc"] is not None


def test_a_file_with_no_hemagglutinin_is_refused():
    model = P.parse_structure(structure_bytes("5k9k"), "5k9k.cif")

    class AntibodyOnly(Select):
        def accept_chain(self, chain):
            return chain.id in {"H", "L"}

    buffer = io.StringIO()
    io_ = PDBIO()
    io_.set_structure(model)
    io_.save(buffer, AntibodyOnly())

    with pytest.raises(P.UnsupportedInput, match="No influenza hemagglutinin"):
        P.predict(buffer.getvalue().encode(), "antibody_only.pdb")


@pytest.mark.parametrize("payload", [b"this is not a structure file", b"", b"\x00\x01\x02\x03"])
def test_garbage_is_refused_with_a_clear_message(payload):
    with pytest.raises(P.UnsupportedInput):
        P.predict(payload, "upload.cif")


@pytest.mark.parametrize("pdb", ["5k9k", "6urm", "3gbm"])
def test_constraint_matches_the_stored_measurement(pdb, stored):
    """Constraint is a lookup by HA number, so it must reproduce exactly what the
    constraint pipeline stored for the same residues, and durability is E x C."""
    table = pd.read_csv(
        ROOT / "results" / "epitope_predictions.csv", low_memory=False,
        usecols=KEY + ["constraint_score"],
    )
    table["insertion_code"] = table["insertion_code"].fillna("").astype(str)
    prediction = P.predict(structure_bytes(pdb), f"{pdb}.cif")
    merged = prediction.residues.assign(PDB=pdb).merge(table, on=KEY, suffixes=("_new", "_old"))

    both = merged.dropna(subset=["constraint_score_new", "constraint_score_old"])
    assert len(both) > 0.8 * len(merged)
    assert np.allclose(both["constraint_score_new"], both["constraint_score_old"], atol=1e-4)

    scored = prediction.residues.dropna(subset=["durability_score"])
    assert np.allclose(scored["durability_score"], scored["epitope_score"] * scored["constraint_score"])


def test_influenza_b_has_no_durability_rather_than_a_borrowed_one():
    prediction = P.predict(structure_bytes("9jno"), "9jno.cif")
    assert prediction.residues["constraint_score"].isna().all()
    assert any("outside the H1 and H3" in note for note in prediction.warnings)


def test_unstandard_numbering_is_not_given_a_constraint():
    """Renumber a real HA 1..N: the lookup would land on the wrong positions."""
    from Bio.PDB import MMCIFParser
    import gzip as gz

    model = P.parse_structure(structure_bytes("5k9k"), "5k9k.cif")
    for chain in model:
        for n, residue in enumerate(list(chain), start=1):
            residue.id = (residue.id[0], 10000 + n, " ")   # avoid clashes while renumbering
        for n, residue in enumerate(list(chain), start=1):
            residue.id = (residue.id[0], n, " ")

    buffer = io.StringIO()
    io_ = PDBIO()
    io_.set_structure(model)
    io_.save(buffer)
    prediction = P.predict(buffer.getvalue().encode(), "renumbered.pdb")
    assert prediction.residues["constraint_score"].isna().all()
    assert any("does not look like standard HA numbering" in note for note in prediction.warnings)


def test_split_head_and_stalk_file_is_judged_by_the_head():
    """7K39 is an H3 stored as separate HA1 and HA2 chains. The stalk is conserved
    across subtypes, so judging by the closest chain would call it a training-branch
    antigen; the head is what diverges and what accuracy depends on."""
    result = P.predict(structure_bytes("7k39"), "7k39.cif").applicability
    assert result["level"] == "held_out"
    assert result["identity_to_training"] < 0.5
