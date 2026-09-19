import type { AntigenSummary, MetricRow } from "../types";
import { SPLIT_LABEL, isHeldOut } from "../data";

interface Props {
  antigen: AntigenSummary;
  metrics: MetricRow[];
}

/** How much to believe the numbers on screen, said before anyone reads them. */
export function TrustBar({ antigen, metrics }: Props) {
  const model = metrics.find((m) => m.split === antigen.split && m.scorer === "epitope_score");
  const baseline = metrics.find((m) => m.split === antigen.split && m.scorer === "rel_sasa_assembly");
  const heldOut = isHeldOut(antigen.split);

  let tone = "neutral";
  let headline = SPLIT_LABEL[antigen.split] ?? antigen.split;
  let detail = "";
  if (antigen.split === "train") {
    tone = "warn";
    detail = "The model was fit on this antigen, so these scores look better than they are. Quote held-out antigens instead.";
  } else if (antigen.split === "val") {
    detail = "Used to tune the model, so slightly optimistic.";
  } else if (antigen.split === "test_group2") {
    tone = "good";
    detail = "From a phylogenetic group the model never saw in training, so these numbers are honest.";
  } else if (antigen.split === "test_B") {
    tone = "warn";
    detail = "Influenza B is far from the training data. Treat the ranking as low confidence; a structure-only model did better here.";
  } else if (antigen.split === "excluded") {
    detail = "Left out of training and evaluation.";
  }

  return (
    <div className={`trust ${tone}`}>
      <div className="trust-head">
        <span className="badge">{headline}</span>
        <span>{detail}</span>
      </div>
      <div className="trust-stats">
        <div>
          <span className="stat-value">{antigen.auprc === null ? "n/a" : antigen.auprc.toFixed(2)}</span>
          <span className="stat-label">
            this chain (AUPRC){antigen.nEpitope ? ` · ${antigen.nEpitope} of ${antigen.nResidues} residues observed` : ""}
          </span>
        </div>
        {model && (
          <div>
            <span className="stat-value">{model.auprc.toFixed(2)}</span>
            <span className="stat-label">
              {heldOut ? "held-out" : "split"} average vs {baseline ? baseline.auprc.toFixed(2) : "?"} for surface exposure alone
              (random ≈ {model.positive_rate.toFixed(2)})
            </span>
          </div>
        )}
      </div>
    </div>
  );
}
