import type { Antigen } from "../types";
import { residueLabel } from "../data";
import { FUNCTION_COLOR, SELECT_COLOR, SITE_COLOR } from "../color";

interface Props {
  antigen: Antigen;
  selected: number[];
  onHover: (i: number | null) => void;
  onSelect: (indices: number[], additive: boolean) => void;
}

const W = 340;
const H = 300;
const M = { l: 38, r: 12, t: 12, b: 34 };

function colourOf(a: Antigen, i: number) {
  if (a.rbs[i]) return FUNCTION_COLOR.rbs;
  if (a.site[i]) return SITE_COLOR[a.site[i]];
  return "#9aa5b4";
}

/**
 * Epitope score against evolutionary constraint. Durable targets sit top-right:
 * antibodies bind them AND the virus cannot change them. Constraint is not
 * computed yet, so until scripts/08 exports it this shows what is missing.
 */
export function DurabilityPlot({ antigen, selected, onHover, onSelect }: Props) {
  const points = antigen.score
    .map((score, i) => ({ i, score, constraint: antigen.constraint[i] }))
    .filter((p): p is { i: number; score: number; constraint: number } => p.constraint !== null);

  if (points.length === 0 && antigen.source) {
    return (
      <div className="panel-body">
        <div className="empty">
          <strong>No durability score for this structure</strong>
          <p>
            Durability requires standard HA numbering and an influenza A subtype in the H1 or H3 reference groups.
            The notes in the banner above state which condition was not met.
          </p>
        </div>
      </div>
    );
  }

  if (points.length === 0) {
    return (
      <div className="panel-body">
        <div className="empty">
          <strong>Constraint score not computed yet</strong>
          <p>
            Durability is <em>epitope score × constraint</em>: a residue is only a good vaccine target if antibodies bind it
            and the virus cannot afford to mutate it.
          </p>
          <ol>
            <li>Align HA sequences from many seasons and compute per-position entropy.</li>
            <li>Set <code>constraint_score = 1 − normalized entropy</code> and <code>durability_score = epitope × constraint</code>.</li>
            <li>Re-run <code>scripts/08_export_web_data.py</code>; this tab fills in on its own.</li>
          </ol>
          <p className="muted">
            The check to look for: antigenic site B scores highest as an epitope but drifts every season, so it should fall
            to the left of the receptor pocket.
          </p>
        </div>
      </div>
    );
  }

  const x = (v: number) => M.l + v * (W - M.l - M.r);
  const y = (v: number) => H - M.b - v * (H - M.t - M.b);
  const chosen = new Set(selected);

  return (
    <div className="panel-body">
      <p className="panel-note">
        Each point is a residue. Upper right marks a candidate durable target: high predicted epitope probability and high evolutionary constraint. Select a point to locate it.
      </p>
      <svg viewBox={`0 0 ${W} ${H}`} className="scatter" role="img" aria-label="Epitope score against constraint">
        <rect x={x(0.5)} y={y(1)} width={x(1) - x(0.5)} height={y(0.5) - y(1)} className="quadrant" />
        <text x={x(1) - 4} y={y(1) + 12} textAnchor="end" className="axis-label">
          candidate durable target
        </text>
        {[0, 0.25, 0.5, 0.75, 1].map((t) => (
          <g key={t}>
            <line x1={M.l} x2={W - M.r} y1={y(t)} y2={y(t)} className="gridline" />
            <text x={M.l - 6} y={y(t) + 3} textAnchor="end" className="axis-label">{t}</text>
            <text x={x(t)} y={H - M.b + 13} textAnchor="middle" className="axis-label">{t}</text>
          </g>
        ))}
        <text x={(M.l + W - M.r) / 2} y={H - 4} textAnchor="middle" className="axis-label">Epitope score</text>
        <text transform={`translate(10 ${(M.t + H - M.b) / 2}) rotate(-90)`} textAnchor="middle" className="axis-label">
          Constraint
        </text>
        {points.map(({ i, score, constraint }) => (
          <circle
            key={i}
            cx={x(score)}
            cy={y(constraint)}
            r={chosen.has(i) ? 5 : 3}
            fill={colourOf(antigen, i)}
            stroke={chosen.has(i) ? SELECT_COLOR : "white"}
            strokeWidth={chosen.has(i) ? 2 : 0.6}
            opacity={0.85}
            style={{ cursor: "pointer" }}
            onMouseEnter={() => onHover(i)}
            onMouseLeave={() => onHover(null)}
            onClick={(e) => onSelect([i], e.shiftKey)}
          >
            <title>{`${residueLabel(antigen, i)}: score ${score.toFixed(2)}, constraint ${constraint.toFixed(2)}`}</title>
          </circle>
        ))}
      </svg>
    </div>
  );
}
