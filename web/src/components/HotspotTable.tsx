import type { Patch } from "../patches";
import { scoreCss } from "../color";

interface Props {
  patches: Patch[];
  loading: boolean;
  selected: number[];
  onPick: (indices: number[]) => void;
  onClear: () => void;
  /** false for an upload, which has no antibody to have been observed touching anything */
  showObserved?: boolean;
}

const sameSet = (a: number[], b: number[]) => a.length === b.length && a.every((v, i) => v === b[i]);

function csv(patches: Patch[]) {
  const rows = patches.map((p) =>
    [p.rank, p.indices.length, p.meanScore.toFixed(3), p.peakScore.toFixed(3), `"${p.span}"`, `"${p.tags.join("; ")}"`, p.observed.toFixed(2)].join(","),
  );
  return ["rank,residues,mean_score,peak_score,ha_numbers,annotation,observed_fraction", ...rows].join("\n");
}

export function HotspotTable({ patches, loading, selected, onPick, onClear, showObserved = true }: Props) {
  const download = () => {
    const url = URL.createObjectURL(new Blob([csv(patches)], { type: "text/csv" }));
    const link = Object.assign(document.createElement("a"), { href: url, download: "hotspots.csv" });
    link.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="panel-body">
      <p className="panel-note">
        Spatially clustered high-scoring residues (Cα within 9 Å), ranked by summed score. Select a patch to focus it in the structure.
      </p>
      {loading && <p className="muted">Loading structure…</p>}
      {!loading && patches.length === 0 && <p className="muted">No compact high-scoring patch on this chain.</p>}
      {patches.length > 0 && (
        <table className="table">
          <thead>
            <tr>
              <th>#</th>
              <th>Residues (H3 numbering)</th>
              <th title="Mean predicted probability">Mean score</th>
              {showObserved && <th title="Fraction of patch residues with an observed antibody contact">Observed</th>}
            </tr>
          </thead>
          <tbody>
            {patches.map((p) => (
              <tr
                key={p.rank}
                className={sameSet(selected, p.indices) ? "row active" : "row"}
                onClick={() => onPick(p.indices)}
                tabIndex={0}
                onKeyDown={(e) => e.key === "Enter" && onPick(p.indices)}
              >
                <td>{p.rank}</td>
                <td>
                  <div className="span">{p.span}</div>
                  <div className="tags">
                    {p.tags.map((t) => (
                      <span key={t} className="tag">
                        {t}
                      </span>
                    ))}
                    <span className="muted">{p.indices.length} residues</span>
                  </div>
                </td>
                <td>
                  <span className="score-chip" style={{ background: scoreCss(p.meanScore) }}>
                    {p.meanScore.toFixed(2)}
                  </span>
                </td>
                {showObserved && <td>{Math.round(p.observed * 100)}%</td>}
              </tr>
            ))}
          </tbody>
        </table>
      )}
      <div className="panel-actions">
        <button className="ghost-button" onClick={onClear} disabled={selected.length === 0}>
          Clear selection{selected.length ? ` (${selected.length})` : ""}
        </button>
        <button className="ghost-button" onClick={download} disabled={patches.length === 0}>
          Export CSV
        </button>
      </div>
    </div>
  );
}
