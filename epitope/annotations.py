"""Functional regions of influenza HA, in H3 numbering.

These are published residue sets, not predictions. They are the yardstick the
model is measured against, so nothing here may be derived from antibody contacts.

HA is made of two pieces. HA1 carries the globular head that grabs the host cell;
HA2 is the stalk that drives membrane fusion. The head drifts every season, the
stalk cannot drift without breaking fusion, which is the whole premise of the
project.

Structures number these inconsistently: some give HA1 and HA2 separate chains
each starting at 1, others number the uncleaved precursor straight through.
classify_chain handles both.
"""

# the head domain is closed by the C52-C277 disulfide; outside it, HA1 folds
# back down into the stalk
HEAD_RANGE = (52, 277)
HA1_LENGTH = 328

# Wiley & Skehel antigenic sites: the patches that drift between seasons and
# where most head antibodies land
ANTIGENIC_SITES = {
    "A": {122, 124, 126, 128, 129, 130, 131, 132, 133, 135, 137, 138, 140,
          142, 143, 144, 145, 146},
    "B": {155, 156, 157, 158, 159, 160, 163, 164, 186, 187, 188, 189, 190,
          192, 193, 194, 196, 197, 198},
    "C": {44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280,
          294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312},
    "D": {96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177,
          179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218,
          219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248},
    "E": {57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94,
          109, 260, 261, 262, 263, 265},
}
ALL_ANTIGENIC = set().union(*ANTIGENIC_SITES.values())

# the receptor binding pocket: the walls that hold sialic acid, plus the four
# conserved residues at its base
RECEPTOR_BINDING_SITE = {
    98, 134, 135, 136, 137, 138, 152, 153, 154, 155, 156, 183,
    188, 189, 190, 191, 192, 193, 194, 195, 221, 222, 223, 224, 225, 226, 227, 228,
}

# HA2 features, in HA2 numbering
FUSION_PEPTIDE = set(range(1, 24))
LONG_HELIX = set(range(76, 131))


# a residue outside these ranges cannot be placed on the standard frame, so it
# is left unannotated rather than given a fabricated number
HA1_FRAME = (1, HA1_LENGTH)
HA2_FRAME = (1, 190)

# depositors number HA2 either straight on from HA1 or at a +1000 offset
HA2_OFFSET = 1000


def numbering_scheme(residue_numbers):
    """How a chain numbers its two pieces.

    Four conventions appear in the PDB, and telling them apart matters: reading
    an offset chain as a continuous one turns HA2 residue 1005 into "HA1 677",
    which then collects meaningless region and antigenic-site labels.
    """
    numbers = sorted(set(residue_numbers))
    highest = numbers[-1]
    if highest >= HA2_OFFSET:
        return "offset"        # HA1 numbered normally, HA2 at +1000
    if highest > HA1_LENGTH + 40:
        return "continuous"    # HA0 numbered straight through
    if highest <= HA2_FRAME[1] and len(numbers) < 220:
        return "HA2"           # the stalk piece alone
    return "HA1"               # the head piece alone


def classify_chain(residue_numbers):
    """Label a chain HA1, HA2 or HA0. Both two-piece conventions are HA0."""
    scheme = numbering_scheme(residue_numbers)
    return "HA0" if scheme in ("offset", "continuous") else scheme


def to_ha_numbering(residue_number, scheme):
    """Map a residue onto (piece, number-within-that-piece, is-on-the-frame)."""
    if scheme == "HA2":
        piece, number = "HA2", residue_number
    elif scheme == "HA1":
        piece, number = "HA1", residue_number
    elif scheme == "offset":
        if residue_number >= HA2_OFFSET:
            piece, number = "HA2", residue_number - HA2_OFFSET
        else:
            piece, number = "HA1", residue_number
    elif residue_number <= HA1_LENGTH:
        piece, number = "HA1", residue_number
    else:
        piece, number = "HA2", residue_number - HA1_LENGTH

    frame = HA1_FRAME if piece == "HA1" else HA2_FRAME
    return piece, number, frame[0] <= number <= frame[1]


def region_of(piece, number, in_frame):
    """head or stem. HA1 outside the head domain folds back into the stalk."""
    if not in_frame:
        return ""
    if piece == "HA2":
        return "stem"
    return "head" if HEAD_RANGE[0] <= number <= HEAD_RANGE[1] else "stem"


def antigenic_site_of(piece, number, in_frame):
    if piece != "HA1" or not in_frame:
        return ""
    for name, residues in ANTIGENIC_SITES.items():
        if number in residues:
            return name
    return ""


def annotate(frame, chain_column="antigen_id", residue_column="residue_number"):
    """Add piece / region / antigenic_site / is_rbs columns to a residue table."""
    frame = frame.copy()
    schemes = frame.groupby(chain_column)[residue_column].apply(
        lambda s: numbering_scheme(s.tolist())
    )
    frame["chain_type"] = frame[chain_column].map(
        schemes.map(lambda s: "HA0" if s in ("offset", "continuous") else s)
    )

    mapped = [to_ha_numbering(n, s)
              for n, s in zip(frame[residue_column], frame[chain_column].map(schemes))]
    frame["piece"] = [p for p, _, _ in mapped]
    frame["ha_number"] = [n for _, n, _ in mapped]
    frame["on_frame"] = [f for _, _, f in mapped]
    frame["region"] = [region_of(p, n, f) for p, n, f in mapped]
    frame["antigenic_site"] = [antigenic_site_of(p, n, f) for p, n, f in mapped]
    frame["is_rbs"] = [
        f and p == "HA1" and n in RECEPTOR_BINDING_SITE for p, n, f in mapped
    ]
    frame["is_fusion_machinery"] = [
        f and p == "HA2" and (n in FUSION_PEPTIDE or n in LONG_HELIX)
        for p, n, f in mapped
    ]
    return frame
