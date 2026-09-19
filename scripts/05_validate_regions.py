"""Does the model agree with what immunologists already know about H3 HA?

Three questions, all answerable from published annotations:

  1. Do real antibody contacts concentrate on the head? (checks the labels)
  2. Does the model put its high scores there too? (checks the model)
  3. Do the model's top residues land in the known antigenic sites A-E?

A model that passes these has learned the biology rather than a dataset quirk.
The point of the project is that the answer to 1 and 2 is "yes, the head" and
that this is exactly the wrong place to aim a vaccine.
"""

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from epitope.annotations import ALL_ANTIGENIC, annotate  # noqa: E402

PREDICTIONS_PATH = "results/predictions_v1.parquet"
OUTPUT_PATH = pathlib.Path("results/region_validation.csv")
SUBTYPE = "H3"
TOP_FRACTION = 0.10

pd.set_option("display.width", 200)

predictions = pd.read_parquet(PREDICTIONS_PATH)
h3 = predictions[predictions["ha_subtype"] == SUBTYPE].copy()
h3 = annotate(h3)
print(f"{SUBTYPE}: {h3['PDB'].nunique()} structures, {h3['antigen_id'].nunique()} chains, "
      f"{len(h3):,} residues")
print(f"all held out of training (split: {sorted(h3['split_group'].unique())})")
print()

print("--- chain types found ---")
print(h3.groupby("chain_type")["antigen_id"].nunique().to_string())
print()

# 1 and 2: head versus stem
print("=== head vs stem ===")
regions = h3.groupby("region").agg(
    residues=("is_epitope", "size"),
    true_epitope_rate=("is_epitope", "mean"),
    mean_model_score=("score", "mean"),
)
regions["score_ratio_vs_stem"] = regions["mean_model_score"] / regions.loc["stem", "mean_model_score"]
print(regions.round(4).to_string())
print()

# where the model's most confident calls actually fall
cutoff = h3["score"].quantile(1 - TOP_FRACTION)
top = h3[h3["score"] >= cutoff]
share = top["region"].value_counts(normalize=True)
baseline = h3["region"].value_counts(normalize=True)
print(f"top {TOP_FRACTION:.0%} of predicted residues ({len(top):,}):")
for region in ("head", "stem"):
    print(f"  {region:<5} {share.get(region, 0):.1%} of predictions "
          f"vs {baseline.get(region, 0):.1%} of all residues "
          f"(enrichment {share.get(region, 0) / baseline.get(region, 1):.2f}x)")
print()

# 3: the classic drift sites
print("=== antigenic sites A-E (HA1 only) ===")
ha1 = h3[h3["piece"] == "HA1"].copy()
ha1["in_antigenic_site"] = ha1["ha_number"].isin(ALL_ANTIGENIC)
sites = ha1.groupby("in_antigenic_site").agg(
    residues=("is_epitope", "size"),
    true_epitope_rate=("is_epitope", "mean"),
    mean_model_score=("score", "mean"),
)
print(sites.round(4).to_string())

top_ha1 = ha1[ha1["score"] >= ha1["score"].quantile(1 - TOP_FRACTION)]
enrichment = top_ha1["in_antigenic_site"].mean() / ha1["in_antigenic_site"].mean()
print(f"\ntop {TOP_FRACTION:.0%} of HA1 predictions are {enrichment:.2f}x enriched "
      f"in antigenic-site residues")
print(f"  ({top_ha1['in_antigenic_site'].mean():.1%} of top predictions "
      f"vs {ha1['in_antigenic_site'].mean():.1%} of all HA1 residues)")
print()

print("--- per site ---")
per_site = ha1[ha1["antigenic_site"] != ""].groupby("antigenic_site").agg(
    residues=("is_epitope", "size"),
    true_epitope_rate=("is_epitope", "mean"),
    mean_model_score=("score", "mean"),
)
print(per_site.round(4).to_string())
print()

# the receptor pocket and the fusion machinery: the two constrained regions
print("=== functionally constrained regions ===")
for column, label in (("is_rbs", "receptor binding site"),
                      ("is_fusion_machinery", "fusion peptide + long helix")):
    subset = h3.groupby(column).agg(
        residues=("is_epitope", "size"),
        true_epitope_rate=("is_epitope", "mean"),
        mean_model_score=("score", "mean"),
    )
    print(f"{label}:")
    print(subset.round(4).to_string())
    print()

h3.to_csv(OUTPUT_PATH, index=False)
print(f"wrote {OUTPUT_PATH}")
