"""Train the baseline epitope model and evaluate it on the held-out HA groups.

Trains on influenza A group 1 HA (H1, H2, H6, H18), validates on H5, and tests on
group 2 (H3, H4, H7, H10, H14) and on influenza B — antigens whose phylogenetic
group the model has never seen.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
from xgboost import XGBClassifier

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from epitope.metrics import best_threshold, comparison_table  # noqa: E402

RESIDUES_PATH = "data/residues_v1.parquet"
FEATURES_PATH = "data/features_v1.parquet"
ESM_PATH = pathlib.Path("data/esm_v1.parquet")  # optional; added when it exists
MODEL_PATH = pathlib.Path("models/xgb_v1.json")
PREDICTIONS_PATH = pathlib.Path("results/predictions_v1.parquet")

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# chain_length is deliberately absent: it is constant within a chain, so the
# model uses it as a construct-size shortcut rather than as biology. protrusion
# is z-scored per chain for the same reason — raw distance from the centroid
# scales with how much of the trimer a given structure happens to contain.
STRUCTURAL = [
    "rel_sasa", "rel_sasa_assembly", "buried_by_assembly", "protrusion_z",
    "glycan_distance", "hydrophobicity", "charge", "volume", "bfactor_z",
    "neighbours_8", "neighbours_12",
    "nbhd_rel_sasa", "nbhd_hydrophobicity", "nbhd_charge",
    "nbhd_protrusion_z", "nbhd_glycan_distance",
]

PARAMS = dict(
    n_estimators=1500, max_depth=6, learning_rate=0.03,
    subsample=0.8, colsample_bytree=0.8, min_child_weight=5,
    reg_lambda=2.0, eval_metric="aucpr", early_stopping_rounds=75,
    tree_method="hist", n_jobs=12, random_state=0,
)

residues = pd.read_parquet(RESIDUES_PATH)
features = pd.read_parquet(FEATURES_PATH)
data = residues.merge(features.drop(columns=["residue"]), on=KEY, how="inner")
print(f"{len(data):,} residues joined ({len(residues) - len(data)} unmatched)")

# HA only; neuraminidase is a different protein and a separate model
data = data[data["target"] == "HA"].copy()

# B-factor scale differs between X-ray and cryo-EM, and protrusion scales with
# how much of the trimer was solved, so both are compared within a chain
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

FEATURE_COLUMNS = STRUCTURAL + list(onehot.columns)

# ESM-2 principal components, when scripts/04_esm_embeddings.py has been run.
# Pass --no-esm to train on structural features alone, for the ablation.
if ESM_PATH.exists() and "--no-esm" not in sys.argv:
    esm = pd.read_parquet(ESM_PATH)
    esm_columns = [c for c in esm.columns if c.startswith("esm_")]
    data = data.merge(esm, on=KEY, how="left")
    missing = data[esm_columns[0]].isna().sum()
    data[esm_columns] = data[esm_columns].fillna(0.0)
    FEATURE_COLUMNS += esm_columns
    print(f"ESM-2: {len(esm_columns)} components ({missing:,} residues without an embedding)")
else:
    print("ESM-2: not found, training on structural features only")

print(f"{len(FEATURE_COLUMNS)} features")

splits = {name: data[data["split_group"] == name] for name in
          ("train", "val", "test_group2", "test_B")}
for name, frame in splits.items():
    print(f"  {name:<12} {frame['antigen_id'].nunique():>4} chains  "
          f"{len(frame):>7,} residues  {frame['is_epitope'].mean():.1%} positive")

train, val = splits["train"], splits["val"]
positive_rate = train["is_epitope"].mean()
model = XGBClassifier(scale_pos_weight=(1 - positive_rate) / positive_rate, **PARAMS)
model.fit(
    train[FEATURE_COLUMNS], train["is_epitope"],
    eval_set=[(val[FEATURE_COLUMNS], val["is_epitope"])],
    verbose=False,
)
print(f"\ntrained; best iteration {model.best_iteration} of {PARAMS['n_estimators']}")

for name, frame in splits.items():
    frame["score"] = model.predict_proba(frame[FEATURE_COLUMNS])[:, 1]

# the surface-exposure baseline: if the model cannot beat one column, nothing else matters
BASELINES = ["rel_sasa_assembly", "protrusion"]
threshold, val_mcc = best_threshold(val, "score")
print(f"threshold tuned on validation: {threshold:.3f} (MCC {val_mcc:.3f})\n")

for name in ("val", "test_group2", "test_B"):
    print(f"=== {name} ===")
    print(comparison_table(splits[name], ["score"] + BASELINES, threshold_frame=val))
    print()

importance = pd.Series(model.feature_importances_, index=FEATURE_COLUMNS)
print("top 12 features by gain:")
print(importance.nlargest(12).round(4).to_string())

suffix = "_structural" if "--no-esm" in sys.argv else ""
model_path = MODEL_PATH.with_stem(MODEL_PATH.stem + suffix)
predictions_path = PREDICTIONS_PATH.with_stem(PREDICTIONS_PATH.stem + suffix)

model_path.parent.mkdir(exist_ok=True)
predictions_path.parent.mkdir(exist_ok=True)
model.save_model(model_path)
keep = KEY + ["antigen_id", "ha_subtype", "split_group", "residue", "is_epitope", "score"]
pd.concat([f[keep] for f in splits.values()]).to_parquet(predictions_path, index=False)
print(f"\nwrote {model_path} and {predictions_path}")
