"""Score every HA residue and export a flat CSV for the rest of the team.

One row per antigen residue, carrying the model's epitope score alongside the
published region annotations. The constraint column is left empty on purpose:
it is the next person's job, and D = E * C cannot be computed without it.
"""

import pathlib
import sys

import pandas as pd
from xgboost import XGBClassifier

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from epitope.annotations import annotate  # noqa: E402
from epitope.metrics import best_threshold, comparison_table  # noqa: E402

RESIDUES_PATH = "data/residues_v1.parquet"
FEATURES_PATH = "data/features_v1.parquet"
ESM_PATH = "data/esm_v1.parquet"
MODEL_PATH = "models/xgb_v1.json"

PREDICTIONS_CSV = pathlib.Path("results/epitope_predictions.csv")
METRICS_CSV = pathlib.Path("results/metrics_summary.csv")

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
STRUCTURAL = [
    "rel_sasa", "rel_sasa_assembly", "buried_by_assembly", "protrusion_z",
    "glycan_distance", "hydrophobicity", "charge", "volume", "bfactor_z",
    "neighbours_8", "neighbours_12",
    "nbhd_rel_sasa", "nbhd_hydrophobicity", "nbhd_charge",
    "nbhd_protrusion_z", "nbhd_glycan_distance",
]

EXPORT_COLUMNS = [
    "antigen_id", "PDB", "antigen_chain", "residue_number", "insertion_code",
    "residue", "ha_subtype", "year", "phylo_group", "split_group",
    "chain_type", "piece", "ha_number", "region", "antigenic_site",
    "is_rbs", "is_fusion_machinery",
    "rel_sasa_assembly", "glycan_distance",
    "is_epitope", "n_antibodies", "n_contacting",
    "epitope_score", "constraint_score", "durability_score",
]

residues = pd.read_parquet(RESIDUES_PATH)
features = pd.read_parquet(FEATURES_PATH)
esm = pd.read_parquet(ESM_PATH)

data = residues.merge(features.drop(columns=["residue"]), on=KEY, how="inner")
data = data[data["target"] == "HA"].copy()
data = data.merge(esm, on=KEY, how="left")

for column in ("bfactor", "protrusion", "nbhd_protrusion"):
    grouped = data.groupby("antigen_id")[column]
    data[f"{column}_z"] = (data[column] - grouped.transform("mean")) / (
        grouped.transform("std") + 1e-6
    )

onehot = pd.get_dummies(data["residue"].where(data["residue"].isin(AMINO_ACIDS), "X"),
                        prefix="aa")
for amino_acid in AMINO_ACIDS:
    if f"aa_{amino_acid}" not in onehot:
        onehot[f"aa_{amino_acid}"] = False
onehot = onehot[[f"aa_{a}" for a in AMINO_ACIDS]].astype(float)
data = pd.concat([data, onehot], axis=1)

esm_columns = [c for c in esm.columns if c.startswith("esm_")]
data[esm_columns] = data[esm_columns].fillna(0.0)
feature_columns = STRUCTURAL + list(onehot.columns) + esm_columns

model = XGBClassifier()
model.load_model(MODEL_PATH)
data["epitope_score"] = model.predict_proba(data[feature_columns])[:, 1].round(5)

data = annotate(data)

# handed over empty; scripts/07 fills these in from the surveillance alignment
data["constraint_score"] = pd.NA
data["durability_score"] = pd.NA

data["rel_sasa_assembly"] = data["rel_sasa_assembly"].round(4)
data["glycan_distance"] = data["glycan_distance"].round(3)

PREDICTIONS_CSV.parent.mkdir(exist_ok=True)
data[EXPORT_COLUMNS].sort_values(["antigen_id", "residue_number"]).to_csv(
    PREDICTIONS_CSV, index=False
)

# headline numbers, same table the README quotes
validation = data[data["split_group"] == "val"]
rows = []
for split in ("val", "test_group2", "test_B"):
    subset = data[data["split_group"] == split]
    table = comparison_table(
        subset, ["epitope_score", "rel_sasa_assembly"],
        threshold_frame=validation.rename(columns={"epitope_score": "epitope_score"}),
    )
    table.insert(0, "split", split)
    table.insert(1, "scorer", table.index)
    rows.append(table)
metrics = pd.concat(rows, ignore_index=True)
metrics.to_csv(METRICS_CSV, index=False)

print(f"{len(data):,} HA residues scored across {data['antigen_id'].nunique()} chains")
print(f"  split_group: {data['split_group'].value_counts().to_dict()}")
print(f"  regions:     {data['region'].value_counts().to_dict()}")
print()
print(metrics.to_string(index=False))
print()
print(f"wrote {PREDICTIONS_CSV} ({PREDICTIONS_CSV.stat().st_size / 1e6:.1f} MB)")
print(f"wrote {METRICS_CSV}")
