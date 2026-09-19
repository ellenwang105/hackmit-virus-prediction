"""Constraint score C(r): how much the virus can afford to change a position.

C is computed, not trained. A position that never varies across two decades of
circulating strains is one the virus cannot afford to lose, regardless of whether
any antibody has been crystallised against it.

Three components, per the spec:

  1. Shannon entropy across a surveillance alignment  (primary signal)
  2. ESM-2 masked marginal                            (optional; fills thin alignments)
  3. Functional-site prior for RBS / fusion machinery (off by default, see below)

The prior is off by default on purpose. The headline validation asks whether
D = E x C ranks the receptor pocket above antigenic site B. Bumping C for
annotated RBS residues would decide that comparison by assumption rather than by
measurement -- the same circularity that kept distance-to-RBS out of the model's
features. Callers who want it get it in a separate column.

Numbering
---------
Surveillance sequences are full HA0 precursors: signal peptide, then HA1, then
HA2. The residue tables downstream are keyed by (piece, ha_number), so every
alignment column is mapped back onto that frame through a reference strain whose
offsets are asserted at load time.
"""

import collections
import math
import time

import numpy as np
import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# HA2 always opens with the fusion peptide; it is the one motif conserved well
# enough to locate the HA1/HA2 boundary in any subtype without alignment.
FUSION_PEPTIDE_MOTIF = "GLFGAIAGFIE"

# Reference strains define the numbering frame. signal_length and ha1_length are
# asserted against the fetched sequence, so a wrong constant fails loudly at load
# rather than silently shifting every position in the join.
REFERENCES = {
    "H3": {
        "accession": "AAA43239.1",        # A/Aichi/2/1968 -- the H3 numbering standard
        "signal_length": 16,
        "ha1_length": 328,
        "first_ha1_residue": "Q",
        "query": "H3N2",
    },
    "H1": {
        "accession": "ACP41105.1",        # A/California/07/2009
        "signal_length": 17,
        "ha1_length": 327,
        "first_ha1_residue": "D",
        "query": "H1N1",
    },
}

# group 1 subtypes take the H1 alignment, group 2 takes H3; anything else falls
# back to H3, which is the deeper and better-sampled of the two
GROUP_1 = {"H1", "H2", "H5", "H6", "H8", "H9", "H11", "H12", "H13", "H16", "H17", "H18"}


# --------------------------------------------------------------------------
# fetching
# --------------------------------------------------------------------------

def _entrez(endpoint, params, retries=3):
    params = dict(params, tool="hackmit-epitope", retmode=params.get("retmode", "json"))
    for attempt in range(retries):
        try:
            response = requests.get(f"{EUTILS}/{endpoint}.fcgi", params=params, timeout=120)
            response.raise_for_status()
            return response
        except requests.RequestException:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def fetch_reference(subtype):
    """Fetch the reference HA0 sequence and verify its numbering offsets."""
    spec = REFERENCES[subtype]
    text = _entrez("efetch", {"db": "protein", "id": spec["accession"],
                              "rettype": "fasta", "retmode": "text"}).text
    sequence = "".join(text.strip().split("\n")[1:])

    signal, ha1_length = spec["signal_length"], spec["ha1_length"]
    ha2_start = sequence.find(FUSION_PEPTIDE_MOTIF)
    if ha2_start < 0:
        raise ValueError(f"{subtype} reference {spec['accession']}: no fusion peptide found")

    # 0-based: HA1 residue 1 sits at index signal_length
    if sequence[signal] != spec["first_ha1_residue"]:
        raise ValueError(
            f"{subtype} reference: HA1 residue 1 is {sequence[signal]!r}, "
            f"expected {spec['first_ha1_residue']!r} -- signal_length is wrong"
        )
    if not signal + ha1_length <= ha2_start:
        raise ValueError(
            f"{subtype} reference: HA1 of {ha1_length} residues overruns HA2, "
            f"which starts at index {ha2_start}"
        )
    return sequence, signal, ha2_start


def position_map(signal_length, ha1_length, ha2_start, reference_length):
    """0-based index in the reference -> (piece, number).

    Positions outside the returned mapping are the signal peptide and the
    cleavage residues between HA1 and HA2: absent from mature structures, so
    they map nowhere and are dropped.
    """
    mapping = {}
    for number in range(1, ha1_length + 1):
        mapping[signal_length + number - 1] = ("HA1", number)
    for index in range(ha2_start, reference_length):
        mapping[index] = ("HA2", index - ha2_start + 1)
    return mapping


def fetch_surveillance(subtype, years, per_year, sleep=0.34):
    """Sequences stratified by year, so no single heavily-sequenced season dominates.

    Entropy is meant to measure drift across seasons. Pulling an unstratified
    sample would weight whichever year happened to be sequenced most.
    """
    spec = REFERENCES[subtype]
    records = []
    for year in years:
        term = (
            '"Influenza A virus"[Organism] AND hemagglutinin[Protein Name] '
            f'AND {spec["query"]}[All Fields] AND ("{year}"[PDAT] : "{year}"[PDAT])'
        )
        search = _entrez("esearch", {"db": "protein", "term": term,
                                     "retmax": per_year, "retmode": "json"})
        ids = search.json()["esearchresult"]["idlist"]
        time.sleep(sleep)
        if not ids:
            continue
        fasta = _entrez("efetch", {"db": "protein", "id": ",".join(ids),
                                   "rettype": "fasta", "retmode": "text"}).text
        for block in fasta.strip().split("\n>"):
            lines = block.lstrip(">").split("\n")
            sequence = "".join(lines[1:]).upper()
            if sequence:
                records.append((year, lines[0], sequence))
        time.sleep(sleep)
    return records


# --------------------------------------------------------------------------
# alignment and entropy
# --------------------------------------------------------------------------

def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    return aligner


def column_counts(sequences, reference):
    """Per reference position, how often each amino acid is observed.

    HA tolerates almost no indels within a subtype, so any sequence already the
    reference's length maps position-for-position. Only the ragged ones -- partial
    depositions, the occasional insertion -- pay for an alignment.
    """
    counts = [collections.Counter() for _ in reference]
    aligner = None
    aligned, direct, failed = 0, 0, 0

    for sequence in sequences:
        if len(sequence) == len(reference):
            for index, residue in enumerate(sequence):
                if residue in AMINO_ACIDS:
                    counts[index][residue] += 1
            direct += 1
            continue

        if aligner is None:
            aligner = _aligner()
        try:
            alignment = aligner.align(reference, sequence)[0]
        except (ValueError, IndexError):
            failed += 1
            continue
        for (ref_start, ref_end), (qry_start, qry_end) in zip(*alignment.aligned):
            for step in range(ref_end - ref_start):
                residue = sequence[qry_start + step]
                if residue in AMINO_ACIDS:
                    counts[ref_start + step][residue] += 1
        aligned += 1

    return counts, {"direct": direct, "aligned": aligned, "failed": failed}


def shannon_entropy(counter):
    """Entropy in bits over the observed amino acids at one position."""
    total = sum(counter.values())
    if total == 0:
        return float("nan")
    entropy = 0.0
    for observed in counter.values():
        p = observed / total
        entropy -= p * math.log2(p)
    return entropy


def normalise_entropy(entropies, percentile=99.0):
    """Map entropy onto 0-1 against the spread actually seen in this alignment.

    Dividing by log2(20) is the textbook move and the wrong one here: HA entropy
    is low nearly everywhere, so every position would land near C = 1 and the
    score would carry no ranking information. Scaling against the 99th percentile
    of observed entropy keeps the variable positions distinguishable, at the cost
    of making C relative to the subtype rather than absolute.
    """
    values = np.asarray(entropies, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return values
    ceiling = np.percentile(finite, percentile)
    if ceiling <= 0:
        return np.where(np.isfinite(values), 1.0, np.nan)
    return np.clip(1.0 - values / ceiling, 0.0, 1.0)


# --------------------------------------------------------------------------
# optional ESM-2 component
# --------------------------------------------------------------------------

def esm_available():
    try:
        import esm  # noqa: F401
        import torch  # noqa: F401
    except ImportError:
        return False
    return True


def esm_masked_marginal(sequence, batch_size=16):
    """Wild-type log-probability at each position, with that position masked.

    A residue the language model is confident about is one the protein family
    does not vary. This works where the alignment is thin, which is why it is
    worth having alongside entropy rather than instead of it.
    """
    import esm
    import torch

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    converter = alphabet.get_batch_converter()
    _, _, tokens = converter([("wt", sequence)])

    scores = np.full(len(sequence), np.nan)
    with torch.no_grad():
        for start in range(0, len(sequence), batch_size):
            positions = range(start, min(start + batch_size, len(sequence)))
            batch = tokens.repeat(len(list(positions)), 1)
            for row, position in enumerate(positions):
                batch[row, position + 1] = alphabet.mask_idx
            logits = model(batch)["logits"]
            for row, position in enumerate(positions):
                log_probs = torch.log_softmax(logits[row, position + 1], dim=-1)
                scores[position] = log_probs[alphabet.get_idx(sequence[position])].item()
    return scores


def normalise_log_probs(scores):
    """Log-probabilities onto 0-1 by rank, which avoids picking an arbitrary floor."""
    values = np.asarray(scores, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return values
    ranked = np.full(values.shape, np.nan)
    order = values[finite].argsort().argsort()
    ranked[finite] = order / max(finite.sum() - 1, 1)
    return ranked


# --------------------------------------------------------------------------
# combination
# --------------------------------------------------------------------------

PRIOR_BUMP = 0.15


def combine(entropy_component, esm_component=None, weights=(0.7, 0.3)):
    """Weighted mean of whichever components are present."""
    entropy_component = np.asarray(entropy_component, dtype=float)
    if esm_component is None:
        return entropy_component

    esm_component = np.asarray(esm_component, dtype=float)
    entropy_weight, esm_weight = weights
    stacked = np.vstack([entropy_component * entropy_weight, esm_component * esm_weight])
    present = np.vstack([np.isfinite(entropy_component) * entropy_weight,
                         np.isfinite(esm_component) * esm_weight])
    total = present.sum(axis=0)
    combined = np.nansum(stacked, axis=0) / np.where(total > 0, total, np.nan)
    return combined


def apply_prior(constraint, is_functional, bump=PRIOR_BUMP):
    """Move annotated functional residues a fixed fraction toward 1.

    Reported separately from the headline C -- see the module docstring.
    """
    constraint = np.asarray(constraint, dtype=float)
    is_functional = np.asarray(is_functional, dtype=bool)
    return np.where(is_functional, constraint + bump * (1.0 - constraint), constraint)
