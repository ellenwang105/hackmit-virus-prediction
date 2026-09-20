import type { Applicability, AntigenSummary, MetricRow } from "../types";
import { SPLIT_LABEL, isHeldOut } from "../data";

interface Props {
  antigen: AntigenSummary;
  metrics: MetricRow[];
  /** set for an uploaded structure, which has no split to be judged by */
  applicability?: Applicability;
  warnings?: string[];
}

/** How much to believe the numbers on screen, said before anyone reads them. */
export function TrustBar({ antigen, metrics, applicability, warnings = [] }: Props) {
  const model = metrics.find((m) => m.split === antigen.split && m.scorer === "epitope_score");
  const baseline = metrics.find((m) => m.split === antigen.split && m.scorer === "rel_sasa_assembly");
  const heldOut = isHeldOut(antigen.split);

  let tone = "neutral";
  let headline = SPLIT_LABEL[antigen.split] ?? antigen.split;
  let detail = "";
  if (antigen.split === "train") {
    tone = "warn";
    detail = "Training antigen. Scores are optimistic because the model was fit on this structure; use held-out antigens for performance estimates.";
  } else if (antigen.split === "val") {
    detail = "Validation antigen, used for model selection. Scores are slightly optimistic.";
  } else if (antigen.split === "test_group2") {
    tone = "good";
    detail = "Held-out antigen from a phylogenetic group excluded from training. Performance here is an unbiased estimate.";
  } else if (antigen.split === "test_B") {
    tone = "warn";
    detail = "Influenza B, phylogenetically distant from the training data (23–27% identity). Rankings are low-confidence; a structure-only model outperformed the full model on this set.";
  } else if (antigen.split === "excluded") {
    detail = "Excluded from training and evaluation.";
  } else if (antigen.split === "upload" && applicability) {
    tone = applicability.level === "far" ? "warn" : applicability.level === "near" ? "good" : "neutral";
    headline = SPLIT_LABEL.upload;
    detail = applicability.message;
  }

  if (antigen.split === "upload") {
    return (
      <div className={`trust ${tone}`}>
        <div className="trust-head">
          <span className="badge">{headline}</span>
          <span>{detail}</span>
        </div>
        <div className="trust-stats">
          {applicability?.expected_auprc != null && (
            <div>
              <span className="stat-value">{applicability.expected_auprc.toFixed(2)}</span>
              <span className="stat-label">
                Expected <abbr title="Area under the precision–recall curve">AUPRC</abbr> at this distance from the
                training data (random ≈ 0.05). Estimated from held-out antigens; a submitted structure has no ground truth.
              </span>
            </div>
          )}
        </div>
        {warnings.length > 0 && (
          <ul className="trust-notes">
            {warnings.map((w) => (
              <li key={w}>{w}</li>
            ))}
          </ul>
        )}
      </div>
    );
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
            <abbr title="Area under the precision–recall curve">AUPRC</abbr>, this chain
            {antigen.nEpitope ? ` · ${antigen.nEpitope} of ${antigen.nResidues} residues are observed contacts` : ""}
          </span>
        </div>
        {model && (
          <div>
            <span className="stat-value">{model.auprc.toFixed(2)}</span>
            <span className="stat-label">
              <abbr title="Area under the precision–recall curve">AUPRC</abbr>, {heldOut ? "held-out" : "split"} average · surface-exposure
              baseline {baseline ? baseline.auprc.toFixed(2) : "?"} · random {model.positive_rate.toFixed(2)}
            </span>
          </div>
        )}
      </div>
    </div>
  );
}
