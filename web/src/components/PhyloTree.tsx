import { useMemo } from "react";
import type { Phylogeny, PhyloNode } from "../types";

interface Props {
  phylogeny: Phylogeny;
  /** subtype of the antigen on screen — the branch that gets highlighted */
  value: string;
  /** subtype the list is filtered to, or "all" — only the All button cares */
  filter: string;
  onChange: (subtype: string) => void;
}

const ROW = 26;
const PAD_TOP = 10;
const PAD_BOTTOM = 6;
// room for a mono subtype name beside its structure count
const LABEL_W = 96;
const LEFT = 4;

function tipsOf(node: PhyloNode, out: PhyloNode[] = []) {
  if (!node.children.length) out.push(node);
  else node.children.forEach((child) => tipsOf(child, out));
  return out;
}

/**
 * Subtypes arranged by real sequence identity.
 *
 * Picking a branch loads a representative antigen of that subtype and narrows
 * the list below to it; the highlighted branch is always whichever antigen is
 * on screen. The tree earns its place because distance from the training set is
 * what governs how much to trust a score — group 1 is what the model learned,
 * and accuracy falls off along these branches.
 */
export function PhyloTree({ phylogeny, value, filter, onChange }: Props) {
  const { tips, height, segments, xOf, yOf } = useMemo(() => {
    const tips = tipsOf(phylogeny.tree);
    const height = PAD_TOP + tips.length * ROW + PAD_BOTTOM;
    const root = phylogeny.tree.height || 1;
    const plotW = 226 - LABEL_W - LEFT;

    const xOf = (h: number) => LEFT + (1 - h / root) * plotW;
    const rows = new Map(tips.map((t, i) => [t, PAD_TOP + i * ROW + ROW / 2]));
    const segments: { d: string }[] = [];

    const walk = (node: PhyloNode): number => {
      if (!node.children.length) return rows.get(node)!;
      const ys = node.children.map(walk);
      const x = xOf(node.height);
      segments.push({ d: `M${x} ${Math.min(...ys)} L${x} ${Math.max(...ys)}` });
      node.children.forEach((child, i) => {
        const childX = child.children.length ? xOf(child.height) : xOf(0);
        segments.push({ d: `M${x} ${ys[i]} L${childX} ${ys[i]}` });
      });
      return (Math.min(...ys) + Math.max(...ys)) / 2;
    };
    walk(phylogeny.tree);

    return { tips, height, segments, xOf, yOf: (t: PhyloNode) => rows.get(t)! };
  }, [phylogeny]);

  return (
    <div className="phylo">
      <div className="phylo-head">
        <h2>Virus type</h2>
        <button
          className={filter === "all" ? "phylo-all active" : "phylo-all"}
          onClick={() => onChange("all")}
        >
          All
        </button>
      </div>

      <svg
        viewBox={`0 0 226 ${height}`}
        className="phylo-svg"
        role="group"
        aria-label="Influenza hemagglutinin subtypes arranged by sequence identity"
      >
        {segments.map((s, i) => (
          <path key={i} d={s.d} className="phylo-branch" />
        ))}

        {tips.map((tip) => {
          const meta = phylogeny.subtypes[tip.name];
          const y = yOf(tip);
          const active = value === tip.name;
          return (
            <g
              key={tip.name}
              className={active ? "phylo-tip active" : "phylo-tip"}
              role="button"
              tabIndex={0}
              aria-pressed={active}
              aria-label={`${tip.name}, ${meta?.structures ?? 0} structures`}
              onClick={() => onChange(tip.name)}
              onKeyDown={(e) => {
                if (e.key === "Enter" || e.key === " ") {
                  e.preventDefault();
                  onChange(tip.name);
                }
              }}
            >
              <rect className="phylo-row" x={xOf(0) - 6} y={y - ROW / 2 + 2} width={232 - xOf(0)} height={ROW - 4} rx={5} />
              <circle className="phylo-dot" cx={xOf(0)} cy={y} r={2.8} />
              <text className="phylo-name" x={xOf(0) + 8} y={y + 4}>
                {tip.name}
              </text>
              <text className="phylo-count" x={222} y={y + 4} textAnchor="end">
                {meta?.structures ?? 0}
              </text>
            </g>
          );
        })}
      </svg>

      <p className="phylo-note">
        Branches closer together are more similar viruses. Numbers are solved structures.
      </p>
    </div>
  );
}
