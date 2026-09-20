import type { Applicability, AntigenSummary, MetricRow } from "../types";
import { isHeldOut } from "../data";

interface Props {
  antigen: AntigenSummary;
  metrics: MetricRow[];
  /** set for a submitted structure, which has no evaluation split to be judged by */
  applicability?: Applicability;
  warnings?: string[];
}

type Level = "high" | "moderate" | "low" | "seen" | "unrated";

const SENTENCE: Record<Level, string> = {
  high: "This virus type is close to the ones the model learned from.",
  moderate: "This virus type differs from the ones the model learned from, but the ranking held up well on similar viruses.",
  low: "This virus is unlike anything the model learned from. Treat the ranking as rough.",
  seen: "The model learned from this structure, so its scores look better here than they would for a new virus.",
  unrated: "This structure was left out of model building and testing.",
};
const BADGE: Record<Level, string> = {
  high: "Confidence · High",
  moderate: "Confidence · Moderate",
  low: "Confidence · Low",
  seen: "Used to build the model",
  unrated: "Not evaluated",
};
const TONE: Record<Level, string> = { high: "good", moderate: "neutral", low: "warn", seen: "warn", unrated: "neutral" };

/** Which of those five the structure on screen gets. Uploads are judged by distance from the training data. */
function levelOf(split: string, applicability?: Applicability): Level {
  if (split === "upload" && applicability) {
    return applicability.level === "near" ? "high" : applicability.level === "held_out" ? "moderate" : "low";
  }
  if (split === "train") return "seen";
  if (split === "val") return "high";
  if (split === "test_group2") return "moderate";
  if (split === "test_B") return "low";
  return "unrated";
}

/**
 * How far to trust the numbers on screen, in a sentence, before anyone reads them.
 * The measurements behind it (AUPRC, baselines) are one click away rather than in the way.
 */
export function TrustBar({ antigen, metrics, applicability, warnings = [] }: Props) {
  const level = levelOf(antigen.split, applicability);
  const upload = antigen.split === "upload";
  const model = metrics.find((m) => m.split === antigen.split && m.scorer === "epitope_score");
  const baseline = metrics.find((m) => m.split === antigen.split && m.scorer === "rel_sasa_assembly");
  const auprc = <abbr title="Area under the precision–recall curve">AUPRC</abbr>;

  return (
    <div className={`trust ${TONE[level]}`}>
      <div className="trust-head">
        <span className="badge">{BADGE[level]}</span>
        <span>{SENTENCE[level]}</span>
      </div>

      {warnings.length > 0 && (
        <ul className="trust-notes">
          {warnings.map((w) => (
            <li key={w}>{w}</li>
          ))}
        </ul>
      )}

      <details className="trust-details">
        <summary>How is this measured?</summary>
        <div className="trust-stats">
          {upload && applicability && (
            <>
              <div>
                <span className="stat-value">{applicability.expected_auprc?.toFixed(2) ?? "n/a"}</span>
                <span className="stat-label">
                  Expected {auprc} at this distance from the training data (chance ≈ 0.05), estimated from other viruses.
                  A submitted structure has no ground truth.
                </span>
              </div>
              <p className="trust-note">{applicability.message}</p>
            </>
          )}
          {!upload && (
            <>
              <div>
                <span className="stat-value">{antigen.auprc === null ? "n/a" : antigen.auprc.toFixed(2)}</span>
                <span className="stat-label">
                  {auprc} on this protein
                  {antigen.nEpitope ? ` · ${antigen.nEpitope} of ${antigen.nResidues} amino acids are known antibody contacts` : ""}
                </span>
              </div>
              {model && (
                <div>
                  <span className="stat-value">{model.auprc.toFixed(2)}</span>
                  <span className="stat-label">
                    {auprc}, {isHeldOut(antigen.split) ? "held-out" : "split"} average · surface-exposure baseline{" "}
                    {baseline ? baseline.auprc.toFixed(2) : "?"} · chance {model.positive_rate.toFixed(2)}
                  </span>
                </div>
              )}
            </>
          )}
        </div>
      </details>
    </div>
  );
}
