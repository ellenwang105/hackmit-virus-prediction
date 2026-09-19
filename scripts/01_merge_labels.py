"""Merge per-antibody epitope labels into one row per antigen residue, then assign splits.

A residue is positive if ANY antibody in ANY complex of that antigen contacts it.
Labelling per antibody instance leaves the same residue positive in one row and
negative in another; training on that teaches the model contradictions.
"""

import pathlib

import pandas as pd

LABELS_PATH = "epitope_labels.csv"
OUTPUT_PATH = pathlib.Path("data/residues_v1.parquet")

# influenza A HA splits into two phylogenetic groups; this is the real
# generalisation axis, since all influenza A HA sits above the 30% identity
# threshold that sequence clustering would use.
GROUP_1 = {"H1", "H2", "H5", "H6", "H8", "H9", "H11", "H12", "H13", "H16", "H17", "H18"}
GROUP_2 = {"H3", "H4", "H7", "H10", "H14", "H15"}

# held out of training as the validation subtype
VAL_SUBTYPE = "H5"
TEMPORAL_CUTOFF = 2020


def phylo_group(subtype):
    if subtype in GROUP_1:
        return "group1"
    if subtype in GROUP_2:
        return "group2"
    if subtype == "B":
        return "B"
    return "other"


def assign_split(row):
    if row.phylo_group == "group1":
        return "val" if row.ha_subtype == VAL_SUBTYPE else "train"
    if row.phylo_group == "group2":
        return "test_group2"
    if row.phylo_group == "B":
        return "test_B"
    return "excluded"


# keep_default_na=False matters: the target column uses "NA" for neuraminidase,
# which pandas otherwise reads as a missing value and blanks out 184 complexes.
labels = pd.read_csv(LABELS_PATH, keep_default_na=False)
labels["target"] = labels["target"].replace("", "NA")
labels["insertion_code"] = labels["insertion_code"].fillna("").astype(str)
print(f"{len(labels):,} per-antibody rows across {labels['PDB'].nunique()} structures")

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]

merged = labels.groupby(KEY, as_index=False).agg(
    target=("target", "first"),
    ha_subtype=("ha_subtype", "first"),
    year=("year", "first"),
    resolution=("resolution", "first"),
    seq_index=("seq_index", "first"),
    residue=("residue", "first"),
    min_distance=("min_distance", "min"),
    is_epitope=("is_epitope", "max"),          # union across antibodies
    n_antibodies=("is_epitope", "size"),       # how many antibodies saw this residue
    n_contacting=("is_epitope", "sum"),        # how many actually bound it
)

merged["phylo_group"] = merged["ha_subtype"].map(phylo_group)
merged["split_group"] = merged.apply(assign_split, axis=1)
merged["split_temporal"] = merged["year"].map(
    lambda y: "train" if int(y) <= TEMPORAL_CUTOFF else "test"
)
merged["antigen_id"] = merged["PDB"] + "_" + merged["antigen_chain"]

OUTPUT_PATH.parent.mkdir(exist_ok=True)
merged.to_parquet(OUTPUT_PATH, index=False)

print(f"{len(merged):,} unique antigen residues")
print(f"positives: {int(merged['is_epitope'].sum()):,} ({100 * merged['is_epitope'].mean():.1f}%)")
print(f"antigen chains: {merged['antigen_id'].nunique()}")
print()

print("--- target ---")
print(merged.groupby("target").agg(
    chains=("antigen_id", "nunique"), residues=("is_epitope", "size"),
    positive_rate=("is_epitope", "mean"),
).round(3))
print()

print("--- phylogenetic split (HA only) ---")
ha = merged[merged["target"] == "HA"]
print(ha.groupby("split_group").agg(
    chains=("antigen_id", "nunique"), structures=("PDB", "nunique"),
    residues=("is_epitope", "size"), positives=("is_epitope", "sum"),
    positive_rate=("is_epitope", "mean"),
).round(3))
print()

print("--- temporal split (HA only) ---")
print(ha.groupby("split_temporal").agg(
    structures=("PDB", "nunique"), residues=("is_epitope", "size"),
    positive_rate=("is_epitope", "mean"),
).round(3))

print(f"\nwrote {OUTPUT_PATH}")
