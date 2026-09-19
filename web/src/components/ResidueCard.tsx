import type { Antigen } from "../types";
import { residueLabel } from "../data";
import { scoreCss } from "../color";

/** One-line detail for the residue under the cursor, so the plot and 3D view never need tooltips. */
export function ResidueCard({ antigen, index }: { antigen: Antigen; index: number | null }) {
  if (index === null) {
    return <div className="residue-card muted">Hover a residue in the structure or the plot for details.</div>;
  }
  const i = index;
  const where = [
    antigen.region[i] === "head" ? "Head" : "Stem",
    antigen.site[i] && `Site ${antigen.site[i]}`,
    antigen.rbs[i] ? "Receptor pocket" : "",
    antigen.fusion[i] ? "Fusion machinery" : "",
  ].filter(Boolean);
  const sasa = antigen.sasa[i];
  const glycan = antigen.glycan[i];

  return (
    <div className="residue-card">
      <span className="residue-name">{residueLabel(antigen, i)}</span>
      <span className="score-chip" style={{ background: scoreCss(antigen.score[i]) }}>
        {antigen.score[i].toFixed(3)}
      </span>
      <span>{where.join(" · ")}</span>
      {sasa !== null && <span className="muted">exposure {sasa.toFixed(2)}</span>}
      {glycan !== null && <span className="muted">glycan {glycan.toFixed(0)} Å away</span>}
      <span className={antigen.epitope[i] ? "truth-tag" : "muted"}>
        {antigen.epitope[i]
          ? `observed epitope (${antigen.nContact[i]} of ${antigen.nAb[i]} antibodies)`
          : `not touched by ${antigen.nAb[i]} antibod${antigen.nAb[i] === 1 ? "y" : "ies"}`}
      </span>
      <span className="muted">
        PDB {antigen.num[i]}
        {antigen.ic[i]}
      </span>
    </div>
  );
}
