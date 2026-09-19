"""Join the constraint score onto the exported predictions and compute D = E x C.

Fills the two columns left empty by scripts/06_export_results.py and runs the
validation checks the project stands on.

  results/epitope_predictions.csv     constraint_score and durability_score filled
  results/durability_validation.csv   the checks, as a table

The headline check is the inversion: the head carries the higher epitope score,
the stem should carry the higher durability. If that does not hold, the thesis
does not hold, and this script says so rather than burying it.
"""

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from epitope.constraint import GROUP_1  # noqa: E402

PREDICTIONS_PATH = pathlib.Path("results/epitope_predictions.csv")
CONSTRAINT_PATH = pathlib.Path("data/constraint_v1.parquet")
VALIDATION_PATH = pathlib.Path("results/durability_validation.csv")

# the conserved floor of the receptor pocket: the residues that actually cannot
# move, as opposed to the pocket walls, which drift with the rest of the head
RBS_CONSERVED_BASE = {98, 153, 183, 195}


def constraint_reference(subtype):
    """Which surveillance alignment speaks for this subtype.

    Influenza B is left unscored. It sits outside both influenza A references,
    and borrowing H3 entropy for it would invent a constraint estimate rather
    than measure one.
    """
    if subtype == "B":
        return None
    if subtype in GROUP_1:
        return "H1"
    return "H3"


def attach(predictions, constraint):
    predictions = predictions.copy()
    predictions["constraint_reference"] = predictions.ha_subtype.map(constraint_reference)

    keyed = constraint.set_index(["reference", "piece", "ha_number"])
    index = pd.MultiIndex.from_arrays([
        predictions.constraint_reference, predictions.piece, predictions.ha_number,
    ])
    predictions["constraint_score"] = keyed.constraint_score.reindex(index).to_numpy()
    predictions["constraint_score_prior"] = keyed.constraint_score_prior.reindex(index).to_numpy()
    predictions["constraint_observations"] = keyed.n_observed.reindex(index).to_numpy()
    predictions["durability_score"] = predictions.epitope_score * predictions.constraint_score
    predictions["escape_risk"] = predictions.epitope_score * (1.0 - predictions.constraint_score)

    # full float repr costs ~15 MB of trailing digits across 197k rows, and the
    # scores are not meaningful past six places
    for column in ("constraint_score", "constraint_score_prior",
                   "durability_score", "escape_risk"):
        predictions[column] = predictions[column].round(6)
    return predictions


def checks(scored):
    """The validation table.

    Each row is one atomic claim, with the two numbers it compares. Claims are
    kept separate rather than conjoined: the head-vs-stem premise and the
    head-vs-stem conclusion fail and pass respectively in this dataset, and a
    combined pass/fail would hide which.
    """
    # the thesis is about H3 drift, and H3 is where the published annotations apply
    h3 = scored[(scored.ha_subtype == "H3") & scored.constraint_score.notna()]
    head, stem = h3[h3.region == "head"], h3[h3.region == "stem"]
    sites = h3[h3.antigenic_site.notna() & (h3.antigenic_site != "")]
    site_b = h3[h3.antigenic_site == "B"]
    pocket = h3[h3.is_rbs]
    base = h3[h3.is_rbs & h3.ha_number.isin(RBS_CONSERVED_BASE) & (h3.piece == "HA1")]

    def row(check, quantity, left_name, left, right_name, right, expect_greater):
        return {
            "check": check, "quantity": quantity,
            "left": left_name, "left_value": left,
            "right": right_name, "right_value": right,
            "expectation": f"{left_name} > {right_name}",
            "passed": bool(left > right),
            "note": "" if expect_greater else "premise, not conclusion",
        }

    return pd.DataFrame([
        # 1. the inversion. The premise is reported separately because SAbDab is
        # enriched for stem-binder structures, so mean E is not in fact higher on
        # the head here -- see results/README.md.
        row("1a. E premise", "mean E", "head", head.epitope_score.mean(),
            "stem", stem.epitope_score.mean(), False),
        row("1b. D inversion", "mean D", "stem", stem.durability_score.mean(),
            "head", head.durability_score.mean(), True),

        # 2. the drift sites are decoys: antibodies land there, the virus escapes
        row("2. sites are decoys", "mean C", "all H3", h3.constraint_score.mean(),
            "sites A-E", sites.constraint_score.mean(), True),

        # 3. the headline comparison
        row("3a. site B wins on E", "mean E", "site B", site_b.epitope_score.mean(),
            "receptor pocket", pocket.epitope_score.mean(), False),
        row("3b. pocket wins on D", "mean D", "receptor pocket",
            pocket.durability_score.mean(), "site B", site_b.durability_score.mean(), True),

        # 4. the pocket walls drift with the head; only the base is immovable
        row("4. conserved base wins on D", "mean D", "Y98/W153/H183/Y195",
            base.durability_score.mean(), "site B", site_b.durability_score.mean(), True),
    ])


def main():
    predictions = pd.read_csv(PREDICTIONS_PATH, low_memory=False)
    constraint = pd.read_parquet(CONSTRAINT_PATH)
    scored = attach(predictions, constraint)

    covered = scored.constraint_score.notna()
    print(f"constraint joined to {covered.sum():,}/{len(scored):,} residues "
          f"({covered.mean():.1%})")
    print("\nunscored by subtype:")
    print(scored[~covered].ha_subtype.value_counts().to_string())

    scored.to_csv(PREDICTIONS_PATH, index=False)
    table = checks(scored)
    table.to_csv(VALIDATION_PATH, index=False)

    print(f"\nwrote {PREDICTIONS_PATH} and {VALIDATION_PATH}\n")
    for _, row in table.iterrows():
        mark = "PASS" if row.passed else "FAIL"
        note = f"  ({row.note})" if row.note else ""
        print(f"{mark}  {row.check}{note}")
        print(f"      {row.quantity}: {row.left} {row.left_value:.4f} "
              f"vs {row.right} {row.right_value:.4f}")


if __name__ == "__main__":
    main()
