"""Score an uploaded antigen structure.

    result = predict(open("5k9k.cif", "rb").read(), "5k9k.cif")
    result.residues        # one row per residue: epitope_score plus HA annotations
    result.applicability   # how far the antigen sits from what the model learned

The model only ever sees the antigen, so the antibody is not an input. If the
file is a complex, the antibody chains are recognised and set aside rather than
scored; if it holds no influenza hemagglutinin at all, that is reported as an
error instead of returning confident numbers for something the model has never
seen.

Feature construction is epitope/features.py, the same code that built the
training table, and the ESM-2 projection uses the basis saved by
scripts/04_esm_embeddings.py. tests/test_predict.py checks that scoring a
structure here reproduces the scores stored in results/.
"""

import functools
import gzip
import io
import pathlib
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# On macOS torch and xgboost each bring their own OpenMP runtime, and the two do
# not coexist in one process. Loaded in this order xgboost is imported first (the
# other order segfaults), and both then run single-threaded (otherwise torch's
# first parallel region deadlocks at 0% CPU). Training never met this because the
# two ran in separate scripts. One structure is ~1,000 rows, so threads buy nothing.
from xgboost import XGBClassifier

from epitope import ESM_COMPONENTS
from epitope.annotations import annotate
from epitope.features import add_derived_features, chain_sequence, structural_features

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "models" / "xgb_v1.json"
PCA_PATH = ROOT / "models" / "esm_pca_v1.npz"
FASTA_PATH = ROOT / "antigen_sequences.fasta"
PREDICTIONS_PATH = ROOT / "results" / "epitope_predictions.csv"
METRICS_PATH = ROOT / "results" / "metrics_summary.csv"

# A chain counts as hemagglutinin when its best local alignment score against a
# known HA passes this. Measured on real chains: antibodies, lysozyme and random
# sequences top out at 53; HA chains start at ~680 (one odd 126-residue
# fragment scores 47). Percent identity does not separate them, because a short
# chain can gap-align to anything at 30-38%, above influenza B's ~25% against A.
HA_MIN_SCORE = 120.0
MIN_CHAIN_LENGTH = 40
# ESM-2 has a 1024-token context, two of which are start and end markers
MAX_ESM_LENGTH = 1022

# where the closest known subtype falls decides which measured accuracy applies
NEAR_IDENTITY = 0.55
FAR_IDENTITY = 0.35

warnings.simplefilter("ignore", PDBConstructionWarning)


class UnsupportedInput(ValueError):
    """The upload cannot be scored; the message says why, in words for the user."""


@dataclass
class ChainReport:
    chain: str
    length: int
    ha_score: float
    is_antigen: bool
    closest_subtype: str
    identity: float
    #: nearest subtype the model was actually trained on, and identity to it
    train_subtype: str
    train_identity: float


@dataclass
class Prediction:
    residues: pd.DataFrame
    chains: list[ChainReport]
    applicability: dict
    warnings: list[str] = field(default_factory=list)

    @property
    def antigen_chains(self):
        return [c.chain for c in self.chains if c.is_antigen]


# ---------------------------------------------------------------- reading input

def _decompress(data, filename):
    if data[:2] == b"\x1f\x8b" or filename.lower().endswith(".gz"):
        try:
            return gzip.decompress(data)
        except OSError as exc:
            raise UnsupportedInput("The file looks gzipped but could not be decompressed.") from exc
    return data


def parse_structure(data, filename="upload"):
    """First model of an mmCIF or PDB file, given as bytes."""
    text = _decompress(data, filename).decode("utf-8", errors="replace")
    head = text.lstrip()[:200]
    is_cif = head.startswith("data_") or filename.lower().replace(".gz", "").endswith((".cif", ".mmcif"))
    parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("upload", io.StringIO(text))
        return next(structure.get_models())
    except Exception as exc:
        raise UnsupportedInput(
            "This is not a readable structure file. Upload a .cif or .pdb "
            "(optionally gzipped)."
        ) from exc


# ------------------------------------------------------------ reference library

@functools.lru_cache(maxsize=1)
def _references():
    """One reference sequence per known subtype: its longest solved chain."""
    sequences, name = {}, None
    for line in FASTA_PATH.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            sequences[name] = []
        elif name:
            sequences[name].append(line)
    sequences = {k: "".join(v) for k, v in sequences.items()}

    table = pd.read_csv(
        PREDICTIONS_PATH, usecols=["antigen_id", "ha_subtype", "split_group"], low_memory=False,
    ).drop_duplicates("antigen_id")
    references, trained_on = {}, set()
    for subtype, group in table.groupby("ha_subtype"):
        if subtype in ("unknown", "chimeric"):
            continue
        known = [i for i in group["antigen_id"] if i in sequences]
        if known:
            references[subtype] = sequences[max(known, key=lambda i: len(sequences[i]))]
            if (group["split_group"] == "train").any():
                trained_on.add(subtype)
    return references, trained_on


@functools.lru_cache(maxsize=2)
def _aligner(mode):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    aligner.mode = mode
    return aligner


def _identity(first, second):
    alignment = _aligner("global").align(first, second)[0]
    matches = sum(1 for a, b in zip(alignment[0], alignment[1]) if a == b and a != "-")
    return matches / min(len(first), len(second))


def assess_chain(chain_id, sequence):
    references, trained_on = _references()
    local = _aligner("local")
    ha_score = max(local.score(sequence, ref) for ref in references.values())
    is_antigen = ha_score >= HA_MIN_SCORE

    subtype, identity, train_subtype, train_identity = "", 0.0, "", 0.0
    if is_antigen:
        identities = {name: _identity(sequence, ref) for name, ref in references.items()}
        subtype = max(identities, key=identities.get)
        identity = identities[subtype]
        # how far this is from what the model learned is what governs how far to
        # trust it, so it is measured against the training subtypes only. Held-out
        # subtypes are in the reference set (to recognise and name them) but were
        # never fitted, and matching one of them says nothing about accuracy.
        train_subtype = max(trained_on, key=identities.get)
        train_identity = identities[train_subtype]
    return ChainReport(chain_id, len(sequence), float(ha_score), is_antigen,
                       subtype, float(identity), train_subtype, float(train_identity))


# ------------------------------------------------------------------ the models

@functools.lru_cache(maxsize=1)
def _booster():
    model = XGBClassifier(n_jobs=1)
    model.load_model(str(MODEL_PATH))
    return model


@functools.lru_cache(maxsize=1)
def _projection():
    if not PCA_PATH.exists():
        raise RuntimeError(
            f"{PCA_PATH.name} is missing. Re-run scripts/04_esm_embeddings.py, which saves "
            "the projection the model was trained on."
        )
    saved = np.load(PCA_PATH)
    return saved["mean"], saved["components"][:ESM_COMPONENTS]


@functools.lru_cache(maxsize=1)
def _esm():
    import esm
    import torch

    torch.set_num_threads(1)   # see the note on OpenMP at the top of this file
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    return model.to(device).eval(), alphabet, device


def embed(sequence):
    """Per-residue ESM-2 features in the model's 128-component space."""
    import torch

    if len(sequence) > MAX_ESM_LENGTH:
        raise UnsupportedInput(
            f"A chain of {len(sequence)} residues is longer than the language model "
            f"can read ({MAX_ESM_LENGTH}). Upload a single hemagglutinin assembly."
        )
    model, alphabet, device = _esm()
    _, _, tokens = alphabet.get_batch_converter()([("chain", sequence)])
    with torch.no_grad():
        out = model(tokens.to(device), repr_layers=[33])
    raw = out["representations"][33][0, 1:len(sequence) + 1].cpu().numpy().astype(np.float32)
    mean, components = _projection()
    return (raw - mean) @ components.T


# ------------------------------------------------------------------- applicability

def _applicability(reports):
    """How far this antigen is from the training data, and the accuracy measured there."""
    antigen = [r for r in reports if r.is_antigen]
    best = max(antigen, key=lambda r: r.train_identity)
    shown = f"{best.closest_subtype}" if best.closest_subtype else "HA"

    if best.train_identity >= NEAR_IDENTITY:
        level, split = "near", "val"
        message = (f"This looks like {shown}. It is {best.train_identity:.0%} identical to {best.train_subtype}, "
                   "a subtype the model was trained on.")
        if best.train_identity >= 0.98:
            message += " It may be a training antigen itself, in which case the scores are optimistic."
    elif best.train_identity >= FAR_IDENTITY:
        level, split = "held_out", "test_group2"
        message = (f"This looks like {shown}. The nearest subtype the model was trained on is "
                   f"{best.train_subtype} at {best.train_identity:.0%}, so this is a branch it has not seen.")
    else:
        level, split = "far", "test_B"
        message = (f"This looks like {shown}. It is only {best.train_identity:.0%} identical to anything "
                   "the model was trained on, as distant as influenza B, where the ranking is low confidence.")

    expected = None
    if METRICS_PATH.exists():
        metrics = pd.read_csv(METRICS_PATH)
        row = metrics[(metrics["split"] == split) & (metrics["scorer"] == "epitope_score")]
        if len(row):
            expected = float(row["auprc"].iloc[0])

    return {
        "level": level,
        "closest_subtype": best.closest_subtype,
        "identity": round(best.identity, 3),
        "nearest_trained_subtype": best.train_subtype,
        "identity_to_training": round(best.train_identity, 3),
        "expected_auprc": expected,
        "message": message,
    }


# --------------------------------------------------------------------- the entry

def predict(data, filename="upload"):
    """Score every hemagglutinin residue in an uploaded structure file."""
    model = parse_structure(data, filename)

    reports = []
    for chain in model:
        sequence = chain_sequence(chain)
        if len(sequence) >= MIN_CHAIN_LENGTH:
            reports.append(assess_chain(chain.id, sequence))
    if not reports:
        raise UnsupportedInput("No protein chains were found in this file.")

    antigen_ids = [r.chain for r in reports if r.is_antigen]
    if not antigen_ids:
        best = max(reports, key=lambda r: r.ha_score)
        raise UnsupportedInput(
            "No influenza hemagglutinin was found in this file. The model only scores "
            f"influenza HA (closest chain {best.chain}, alignment score {best.ha_score:.0f}; "
            f"HA scores at least {HA_MIN_SCORE:.0f})."
        )

    # unique sequences are embedded once; a trimer is three identical chains
    sequences = {c.id: chain_sequence(c) for c in model if c.id in antigen_ids}
    cache = {seq: embed(seq) for seq in set(sequences.values())}

    features = structural_features(model, antigen_ids, label="upload")
    features["antigen_id"] = "upload_" + features["antigen_chain"]

    projected = np.vstack([
        cache[sequences[chain]][group["seq_index"].to_numpy()]
        for chain, group in features.groupby("antigen_chain", sort=False)
    ])
    # groupby preserves row order within each chain, and the chains were emitted
    # contiguously, so the stacked rows line up with `features`
    features = pd.concat(
        [features, pd.DataFrame(projected, columns=[f"esm_{i}" for i in range(projected.shape[1])])],
        axis=1,
    )

    frame = add_derived_features(features, chain_column="antigen_id")

    booster = _booster()
    names = list(booster.get_booster().feature_names)
    missing = [n for n in names if n not in frame.columns]
    if missing:
        raise RuntimeError(f"model expects features this build does not produce: {missing[:5]}")
    frame["epitope_score"] = booster.predict_proba(frame[names])[:, 1]

    frame = annotate(frame, chain_column="antigen_id")
    frame["constraint_score"] = np.nan   # needs the surveillance alignment; not wired up yet
    frame["durability_score"] = np.nan

    columns = [
        "antigen_chain", "residue_number", "insertion_code", "residue", "seq_index",
        "epitope_score", "constraint_score", "durability_score",
        "rel_sasa_assembly", "glycan_distance",
        "chain_type", "piece", "ha_number", "on_frame", "region", "antigenic_site",
        "is_rbs", "is_fusion_machinery",
    ]
    notes = []
    ignored = [r.chain for r in reports if not r.is_antigen]
    if ignored:
        notes.append(f"Chains {', '.join(ignored)} were set aside as not hemagglutinin (antibody or other).")
    off_frame = int((~frame["on_frame"]).sum())
    if off_frame:
        notes.append(f"{off_frame} residues fall outside the standard HA numbering and have no region label.")

    return Prediction(
        residues=frame[columns].reset_index(drop=True),
        chains=reports,
        applicability=_applicability(reports),
        warnings=notes,
    )
