import type { Antigen } from "../types";
import { isUpload } from "../data";
import { SITE_COLOR, TRUTH_COLOR } from "../color";

interface Group {
  label: string;
  color: string;
  test: (a: Antigen, i: number) => boolean;
}

const GROUPS: Group[] = [
  { label: "Head", color: "#7f8ea3", test: (a, i) => a.region[i] === "head" },
  { label: "Stem", color: "#c2cad6", test: (a, i) => a.region[i] === "stem" },
  ...["A", "B", "C", "D", "E"].map((s) => ({
    label: `Site ${s}`,
    color: SITE_COLOR[s],
    test: (a: Antigen, i: number) => a.site[i] === s,
  })),
  { label: "Receptor-binding site", color: "#d1495b", test: (a, i) => a.rbs[i] === 1 },
  { label: "Fusion machinery", color: "#3b7ea1", test: (a, i) => a.fusion[i] === 1 },
];

/** Predicted mean score next to the observed contact rate, per structural region. */
export function RegionSummary({ antigen }: { antigen: Antigen }) {
  const showObserved = !isUpload(antigen);
  const rows = GROUPS.map((g) => {
    const members = antigen.score.map((_, i) => i).filter((i) => g.test(antigen, i));
    if (!members.length) return null;
    const mean = members.reduce((s, i) => s + antigen.score[i], 0) / members.length;
    const observed = members.reduce((s, i) => s + antigen.epitope[i], 0) / members.length;
    return { ...g, n: members.length, mean, observed };
  }).filter((r): r is NonNullable<typeof r> => r !== null);

  const top = Math.max(0.05, ...rows.flatMap((r) => (showObserved ? [r.mean, r.observed] : [r.mean])));

  return (
    <div className="panel-body">
      <p className="panel-note">
        {showObserved
          ? "Mean predicted probability (solid) against the observed contact rate (teal) for each structural region. "
          : "Mean predicted probability in each structural region. A submitted structure has no bound antibody, so there is no observed rate to compare. "}
        Regions overlap: a residue can belong to an antigenic site and the receptor-binding site.
      </p>
      <div className="bars">
        {rows.map((r) => (
          <div key={r.label} className="bar-row">
            <div className="bar-label">
              <span className="swatch" style={{ background: r.color }} />
              {r.label}
              <span className="muted"> · {r.n}</span>
            </div>
            <div className="bar-tracks">
              <div className="bar-track">
                <div className="bar-fill" style={{ width: `${(r.mean / top) * 100}%`, background: r.color }} />
                <span className="bar-value">{r.mean.toFixed(3)}</span>
              </div>
              {showObserved && (
                <div className="bar-track thin">
                  <div className="bar-fill" style={{ width: `${(r.observed / top) * 100}%`, background: TRUTH_COLOR }} />
                  <span className="bar-value">{r.observed.toFixed(3)}</span>
                </div>
              )}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
