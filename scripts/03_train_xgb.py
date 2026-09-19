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
from epitope import ESM_COMPONENTS  # noqa: E402
from epitope.metrics import best_threshold, comparison_table  # noqa: E402

RESIDUES_PATH = "data/residues_v1.parquet"
FEATURES_PATH = "data/features_v1.parquet"
ESM_PATH = pathlib.Path("data/esm_v1.parquet")  # optional; added when it exists
MODEL_PATH = pathlib.Path("models/xgb_v1.json")
PREDICTIONS_PATH = pathlib.Path("results/predictions_v1.parquet")
METRICS_PATH = pathlib.Path("results/training_metrics.csv")
IMPORTANCE_PATH = pathlib.Path("results/importance_v1.csv")

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# override for a single run with --esm-components N
DEFAULT_COMPONENTS = ESM_COMPONENTS

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
use_esm = ESM_PATH.exists() and "--no-esm" not in sys.argv
if use_esm:
    esm = pd.read_parquet(ESM_PATH)
    esm_columns = [c for c in esm.columns if c.startswith("esm_")]
    # Principal components are ordered by variance explained and nested, so
    # keeping the leading N is equivalent to having fitted the PCA at width N,
    # and scripts/04 can stay at its widest setting.
    #
    # 128 covers 81% of ESM's variance against 71% at the 64 this started at,
    # which is worth ~3% AUPRC on the unseen phylogenetic group. 256 scores
    # higher still on the test splits but lower on validation, so it is not
    # the default: picking it would mean selecting on the test set.
    width = DEFAULT_COMPONENTS
    if "--esm-components" in sys.argv:
        width = int(sys.argv[sys.argv.index("--esm-components") + 1])
    esm_columns = sorted(esm_columns, key=lambda c: int(c.split("_")[1]))[:width]
    esm = esm[KEY + esm_columns]
    data = data.merge(esm, on=KEY, how="left")
    missing = data[esm_columns[0]].isna().sum()
    data[esm_columns] = data[esm_columns].fillna(0.0)
    FEATURE_COLUMNS += esm_columns
    print(f"ESM-2: {len(esm_columns)} components ({missing:,} residues without an embedding)")
else:
    print("ESM-2: not found, training on structural features only")

# The default run keeps the plain names, so models/xgb_v1.json stays the file
# scripts/05 and 06 load; only an explicit override gets its own filenames.
if not use_esm:
    variant, suffix = "structural", "_structural"
elif "--esm-components" in sys.argv:
    variant, suffix = f"esm{len(esm_columns)}", f"_esm{len(esm_columns)}"
else:
    variant, suffix = "esm", ""

# --union-labels trains against is_epitope_union, which pools each position's
# label over every chain with an identical sequence. Most chains have only one
# antibody solved, so the per-structure negatives are full of epitopes nobody
# has crystallised yet; pooling recovers some of them.
LABEL = "is_epitope_union" if "--union-labels" in sys.argv else "is_epitope"
if LABEL not in data.columns:
    raise SystemExit("is_epitope_union missing — re-run scripts/01_merge_labels.py")
if LABEL != "is_epitope":
    variant += "_union"
    suffix += "_union"
    print(f"labels: {LABEL} ({data[LABEL].mean():.1%} positive, "
          f"vs {data['is_epitope'].mean():.1%} per-structure)")

print(f"{len(FEATURE_COLUMNS)} features")

splits = {name: data[data["split_group"] == name] for name in
          ("train", "val", "test_group2", "test_B")}
for name, frame in splits.items():
    print(f"  {name:<12} {frame['antigen_id'].nunique():>4} chains  "
          f"{len(frame):>7,} residues  {frame[LABEL].mean():.1%} positive")

train, val = splits["train"], splits["val"]


def redundancy_weights(frame):
    """One unit of influence per antigen, not per deposited chain.

    Training has 259 chains but only 29 sequence clusters, so a popular antigen
    solved 44 times would otherwise count 44x as much as one solved once, and
    the model fits crystallographers' interests rather than biology.
    """
    chains_per_cluster = frame.groupby("seq_cluster")["antigen_id"].transform("nunique")
    weights = 1.0 / chains_per_cluster
    return weights / weights.mean()


if "--weight-redundancy" in sys.argv:
    if "seq_cluster" not in data.columns:
        raise SystemExit("seq_cluster missing — re-run scripts/01_merge_labels.py")
    variant += "_weighted"
    suffix += "_weighted"
    train_weights, val_weights = redundancy_weights(train), redundancy_weights(val)
    print(f"redundancy weighting: {train['seq_cluster'].nunique()} clusters over "
          f"{train['antigen_id'].nunique()} train chains, "
          f"weights {train_weights.min():.2f}-{train_weights.max():.2f}")
else:
    train_weights = val_weights = None

positive_rate = train[LABEL].mean()
model = XGBClassifier(scale_pos_weight=(1 - positive_rate) / positive_rate, **PARAMS)
model.fit(
    train[FEATURE_COLUMNS], train[LABEL], sample_weight=train_weights,
    eval_set=[(val[FEATURE_COLUMNS], val[LABEL])],
    sample_weight_eval_set=None if val_weights is None else [val_weights],
    verbose=False,
)
print(f"\ntrained; best iteration {model.best_iteration} of {PARAMS['n_estimators']}")

for name, frame in splits.items():
    frame["score"] = model.predict_proba(frame[FEATURE_COLUMNS])[:, 1]

# the surface-exposure baseline: if the model cannot beat one column, nothing else matters
BASELINES = ["rel_sasa_assembly", "protrusion"]
threshold, val_mcc = best_threshold(val, "score", LABEL)
print(f"threshold tuned on validation: {threshold:.3f} (MCC {val_mcc:.3f})\n")

# Scored against both label sets. "chain" is the per-structure label every
# earlier run used, so it is the only column comparable across variants;
# "union" is the pooled label and is the fairer target when it was trained on.
LABEL_SETS = {"chain": "is_epitope", "union": "is_epitope_union"}

tables = []
for name in ("val", "test_group2", "test_B"):
    for label_set, column in LABEL_SETS.items():
        table = comparison_table(splits[name], ["score"] + BASELINES,
                                 threshold_frame=val, label_column=column)
        if column == LABEL:
            print(f"=== {name}  (labels: {label_set}) ===")
            print(table)
            print()
        table.insert(0, "variant", variant)
        table.insert(1, "label_set", label_set)
        table.insert(2, "split", name)
        table.insert(3, "scorer", table.index)
        tables.append(table)

metrics = pd.concat(tables, ignore_index=True)
metrics["threshold"] = round(threshold, 4)
metrics["n_features"] = len(FEATURE_COLUMNS)
metrics["esm_components"] = len(esm_columns) if use_esm else 0
metrics["best_iteration"] = model.best_iteration

importance = pd.Series(model.feature_importances_, index=FEATURE_COLUMNS)
print("top 12 features by gain:")
print(importance.nlargest(12).round(4).to_string())

model_path = MODEL_PATH.with_stem(MODEL_PATH.stem + suffix)
predictions_path = PREDICTIONS_PATH.with_stem(PREDICTIONS_PATH.stem + suffix)
importance_path = IMPORTANCE_PATH.with_stem(IMPORTANCE_PATH.stem + suffix)

model_path.parent.mkdir(exist_ok=True)
predictions_path.parent.mkdir(exist_ok=True)
model.save_model(model_path)
keep = KEY + ["antigen_id", "ha_subtype", "split_group", "residue",
              "is_epitope", "is_epitope_union", "score"]
pd.concat([f[keep] for f in splits.values()]).to_parquet(predictions_path, index=False)

importance.round(6).rename("gain").rename_axis("feature").sort_values(
    ascending=False).to_csv(importance_path)

# One row per (variant, split, scorer). Re-running a variant replaces only its
# own rows, so the structural and ESM runs accumulate into one comparable file
# instead of each overwriting the other.
if METRICS_PATH.exists():
    previous = pd.read_csv(METRICS_PATH)
    # runs from before --union-labels existed were scored on the per-structure label
    if "label_set" not in previous.columns:
        previous.insert(1, "label_set", "chain")
    metrics = pd.concat([previous[previous["variant"] != variant], metrics],
                        ignore_index=True)
metrics.to_csv(METRICS_PATH, index=False)

print(f"\nwrote {model_path}")
print(f"      {predictions_path}")
print(f"      {importance_path}")
print(f"      {METRICS_PATH}  ({metrics['variant'].nunique()} variant(s): "
      f"{', '.join(sorted(metrics['variant'].unique()))})")
