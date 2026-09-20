import type { Antigen } from "../types";
import { isUpload, residueLabel } from "../data";
import { scoreCss } from "../color";

/** One-line detail for the residue under the cursor, so the plot and 3D view never need tooltips. */
export function ResidueCard({ antigen, index }: { antigen: Antigen; index: number | null }) {
  if (index === null) {
    return <div className="residue-card muted">Hover a residue in the structure or track for details.</div>;
  }
  const i = index;
  const where = [
    antigen.region[i] === "head" ? "Head" : "Stem",
    antigen.site[i] && `Site ${antigen.site[i]}`,
    antigen.rbs[i] ? "Receptor-binding site" : "",
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
      {sasa !== null && <span className="muted">rel. SASA {sasa.toFixed(2)}</span>}
      {glycan !== null && <span className="muted">nearest glycan {glycan.toFixed(0)} Å</span>}
      {!isUpload(antigen) && (
        <span className={antigen.epitope[i] ? "truth-tag" : "muted"}>
          {antigen.epitope[i]
            ? `observed contact (${antigen.nContact[i]} of ${antigen.nAb[i]} antibodies)`
            : `no contact in ${antigen.nAb[i]} antibody complex${antigen.nAb[i] === 1 ? "" : "es"}`}
        </span>
      )}
      <span className="muted">
        PDB residue {antigen.num[i]}
        {antigen.ic[i]}
      </span>
    </div>
  );
}
