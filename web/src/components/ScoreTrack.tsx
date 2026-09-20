import { useRef, useState } from "react";
import type { Antigen } from "../types";
import { isUpload } from "../data";
import { useElementWidth } from "../hooks";
import { FUNCTION_COLOR, REGION_COLOR, SELECT_COLOR, SITE_COLOR, TRUTH_COLOR, scoreCss } from "../color";

interface Props {
  antigen: Antigen;
  hovered: number | null;
  selected: number[];
  onHover: (i: number | null) => void;
  /** `additive` adds to / toggles within the current selection instead of replacing it */
  onSelect: (indices: number[], additive: boolean) => void;
}

const LEFT = 68;
const RIGHT = 12;
const STRIP = 9;
const GAP = 3;
const PLOT_H = 120;
const AXIS_H = 20;
const BASE_RATE = 0.05;

/** Runs of equal, non-empty values, for drawing annotation strips. */
function runs(values: string[]): { start: number; end: number; value: string }[] {
  const out: { start: number; end: number; value: string }[] = [];
  for (let i = 0; i < values.length; ) {
    let j = i;
    while (j + 1 < values.length && values[j + 1] === values[i]) j++;
    if (values[i]) out.push({ start: i, end: j, value: values[i] });
    i = j + 1;
  }
  return out;
}

export function ScoreTrack({ antigen, hovered, selected, onHover, onSelect }: Props) {
  const [ref, width] = useElementWidth<HTMLDivElement>();
  const [drag, setDrag] = useState<{ from: number; to: number } | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);

  const n = antigen.score.length;
  const plotW = Math.max(0, width - LEFT - RIGHT);
  const step = n ? plotW / n : 0;

  const stripsTop = 4;
  const rows = ["Region", "Site", "Function"];
  const plotTop = stripsTop + rows.length * (STRIP + GAP) + 6;
  const truthTop = plotTop + PLOT_H + 6;
  const height = truthTop + STRIP + AXIS_H + 4;

  const y = (score: number) => plotTop + PLOT_H * (1 - Math.min(1, score));
  const indexAt = (clientX: number) => {
    const box = svgRef.current!.getBoundingClientRect();
    return Math.max(0, Math.min(n - 1, Math.floor((clientX - box.left - LEFT) / step)));
  };

  const region = antigen.region.map((r) => r as string);
  const site = antigen.site;
  const fn = antigen.rbs.map((r, i) => (r ? "rbs" : antigen.fusion[i] ? "fusion" : ""));
  const selectedSet = new Set(selected);
  const wanted = n / Math.max(3, Math.floor(plotW / 70));
  const tickEvery = [10, 20, 25, 50, 100].find((t) => t >= wanted) ?? 100;

  const finishDrag = (event: React.PointerEvent) => {
    if (!drag) return;
    const lo = Math.min(drag.from, drag.to);
    const hi = Math.max(drag.from, drag.to);
    onSelect(
      Array.from({ length: hi - lo + 1 }, (_, k) => lo + k),
      event.shiftKey || event.ctrlKey || event.metaKey,
    );
    setDrag(null);
  };

  return (
    <div ref={ref} className="track">
      {width > 0 && (
        <svg ref={svgRef} width={width} height={height} role="img" aria-label="Predicted epitope probability along the chain">
          {/* annotation strips */}
          {rows.map((label, r) => (
            <text key={label} x={LEFT - 6} y={stripsTop + r * (STRIP + GAP) + STRIP - 1} textAnchor="end" className="axis-label">
              {label}
            </text>
          ))}
          {runs(region).map((run) => (
            <rect
              key={`r${run.start}`}
              x={LEFT + run.start * step}
              y={stripsTop}
              width={(run.end - run.start + 1) * step}
              height={STRIP}
              fill={REGION_COLOR[run.value as "head" | "stem"]}
            />
          ))}
          {runs(site).map((run) => {
            const w = (run.end - run.start + 1) * step;
            return (
              <g key={`s${run.start}`}>
                <rect x={LEFT + run.start * step} y={stripsTop + STRIP + GAP} width={w} height={STRIP} fill={SITE_COLOR[run.value]} />
                {w > 9 && (
                  <text x={LEFT + run.start * step + w / 2} y={stripsTop + 2 * STRIP + GAP - 1} textAnchor="middle" className="strip-letter">
                    {run.value}
                  </text>
                )}
              </g>
            );
          })}
          {runs(fn).map((run) => (
            <rect
              key={`f${run.start}`}
              x={LEFT + run.start * step}
              y={stripsTop + 2 * (STRIP + GAP)}
              width={(run.end - run.start + 1) * step}
              height={STRIP}
              fill={FUNCTION_COLOR[run.value as "rbs" | "fusion"]}
            />
          ))}

          {/* y axis + base-rate line */}
          {[0, 0.25, 0.5, 0.75, 1].map((t) => (
            <g key={t}>
              <line x1={LEFT} x2={LEFT + plotW} y1={y(t)} y2={y(t)} className="gridline" />
              <text x={LEFT - 6} y={y(t) + 3} textAnchor="end" className="axis-label">
                {t}
              </text>
            </g>
          ))}
          <line x1={LEFT} x2={LEFT + plotW} y1={y(BASE_RATE)} y2={y(BASE_RATE)} className="baseline" />

          {/* score bars */}
          {antigen.score.map((score, i) => (
            <rect
              key={i}
              x={LEFT + i * step}
              y={y(score)}
              width={Math.max(1, step - 0.4)}
              height={PLOT_H * Math.min(1, score)}
              fill={scoreCss(score)}
              className={selectedSet.has(i) ? "bar selected" : "bar"}
              stroke={selectedSet.has(i) ? SELECT_COLOR : "none"}
              strokeWidth={selectedSet.has(i) ? 1.2 : 0}
            />
          ))}

          {/* observed epitope ticks; an upload has no antibody, so no row */}
          {!isUpload(antigen) && (
            <text x={LEFT - 6} y={truthTop + STRIP - 1} textAnchor="end" className="axis-label">
              Observed
            </text>
          )}
          {!isUpload(antigen) && antigen.epitope.map(
            (e, i) =>
              e === 1 && (
                <rect key={i} x={LEFT + i * step} y={truthTop} width={Math.max(1, step - 0.4)} height={STRIP} fill={TRUTH_COLOR} />
              ),
          )}

          {/* x axis in standard HA numbering */}
          {antigen.ha.map((ha, i) =>
            ha % tickEvery === 0 && (i === 0 || antigen.ha[i - 1] !== ha) ? (
              <text key={`t${i}`} x={LEFT + (i + 0.5) * step} y={truthTop + STRIP + 14} textAnchor="middle" className="axis-label">
                {ha}
              </text>
            ) : null,
          )}

          {/* hover column and drag preview */}
          {hovered !== null && (
            <rect x={LEFT + hovered * step} y={stripsTop} width={Math.max(1, step)} height={height - AXIS_H - stripsTop} className="hover-column" />
          )}
          {drag && (
            <rect
              x={LEFT + Math.min(drag.from, drag.to) * step}
              y={plotTop}
              width={(Math.abs(drag.to - drag.from) + 1) * step}
              height={PLOT_H}
              fill={SELECT_COLOR}
              opacity={0.15}
            />
          )}

          <rect
            x={LEFT}
            y={0}
            width={plotW}
            height={height - AXIS_H}
            fill="transparent"
            style={{ cursor: "crosshair", touchAction: "none" }}
            onPointerDown={(e) => {
              e.currentTarget.setPointerCapture(e.pointerId);
              const i = indexAt(e.clientX);
              setDrag({ from: i, to: i });
            }}
            onPointerMove={(e) => {
              const i = indexAt(e.clientX);
              onHover(i);
              if (drag) setDrag({ ...drag, to: i });
            }}
            onPointerUp={finishDrag}
            onPointerLeave={() => !drag && onHover(null)}
          />
        </svg>
      )}
      <div className="track-legend">
        {[
          ["Head", REGION_COLOR.head],
          ["Stem", REGION_COLOR.stem],
          ...Object.entries(SITE_COLOR).map(([s, c]) => [`Site ${s}`, c]),
          ["Receptor-binding site", FUNCTION_COLOR.rbs],
          ["Fusion machinery", FUNCTION_COLOR.fusion],
          ...(isUpload(antigen) ? [] : [["Observed contact", TRUTH_COLOR]]),
        ].map(([label, color]) => (
          <span key={label}>
            <i className="key-swatch" style={{ background: color }} />
            {label}
          </span>
        ))}
      </div>
      <div className="track-caption">
        H3 numbering · dashed line: 5% base rate of epitope residues · click or drag to select; Shift to add
      </div>
    </div>
  );
}
