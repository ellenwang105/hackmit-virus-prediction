import { useMemo, useState } from "react";
import type { AntigenSummary, Phylogeny } from "../types";
import { SPLIT_LABEL, isHeldOut } from "../data";
import { PhyloTree } from "./PhyloTree";

interface Props {
  index: AntigenSummary[];
  phylogeny: Phylogeny | null;
  selectedId: string;
  onSelect: (id: string) => void;
}

const SPLIT_ORDER = ["test_group2", "test_B", "val", "train", "excluded"];

/** Enough observed contacts that a per-chain AUPRC means something. */
const SCOREABLE = 15;

/**
 * The chain to land on when a subtype is picked from the tree.
 *
 * Prefers held-out chains, since those are the ones worth quoting, and among
 * chains with enough observed contacts to score takes the one nearest the
 * median rather than the best — a subtype should open on a representative
 * example, not its most flattering one.
 */
function representativeOf(index: AntigenSummary[], subtype: string): AntigenSummary | null {
  const all = index.filter((a) => a.subtype === subtype);
  if (!all.length) return null;

  const heldOut = all.filter((a) => isHeldOut(a.split));
  const pool = heldOut.length ? heldOut : all;

  const scoreable = pool.filter((a) => a.nEpitope >= SCOREABLE && a.auprc !== null);
  if (!scoreable.length) {
    return pool.reduce((best, a) => (a.nEpitope > best.nEpitope ? a : best));
  }

  const sorted = [...scoreable].sort((a, b) => a.auprc! - b.auprc!);
  return sorted[Math.floor(sorted.length / 2)];
}

export function AntigenPicker({ index, phylogeny, selectedId, onSelect }: Props) {
  const [query, setQuery] = useState("");
  const [subtype, setSubtype] = useState("all");
  const [split, setSplit] = useState("all");

  // the tree covers the subtypes with a coherent reference sequence; the rest
  // (unknown, chimeric constructs) stay reachable from this list
  const strays = useMemo(
    () =>
      phylogeny
        ? [...new Set(index.map((a) => a.subtype))].filter((s) => !(s in phylogeny.subtypes)).sort()
        : [],
    [index, phylogeny],
  );

  const shown = useMemo(() => {
    const q = query.trim().toLowerCase();
    return index
      .filter((a) => (subtype === "all" || a.subtype === subtype) && (split === "all" || a.split === split))
      .filter((a) => !q || a.id.toLowerCase().includes(q))
      .sort((a, b) => SPLIT_ORDER.indexOf(a.split) - SPLIT_ORDER.indexOf(b.split) || a.id.localeCompare(b.id));
  }, [index, query, subtype, split]);

  return (
    <aside className="picker">
      {phylogeny && (
        <PhyloTree
          phylogeny={phylogeny}
          value={index.find((a) => a.id === selectedId)?.subtype ?? "all"}
          filter={subtype}
          onChange={(next) => {
            setSubtype(next);
            // picking a branch should show that virus, not just narrow the list
            if (next !== "all") {
              const representative = representativeOf(index, next);
              if (representative) onSelect(representative.id);
            }
          }}
        />
      )}

      <h2>Antigen</h2>
      <input
        className="input"
        type="search"
        placeholder="Search PDB id, e.g. 3sdy"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
        aria-label="Search antigens"
      />
      <div className="filters">
        {strays.length > 0 && (
          <select className="input" value={subtype} onChange={(e) => setSubtype(e.target.value)} aria-label="Other subtypes">
            <option value="all">Other subtypes</option>
            {strays.map((s) => (
              <option key={s} value={s}>
                {s}
              </option>
            ))}
          </select>
        )}
        <select className="input" value={split} onChange={(e) => setSplit(e.target.value)} aria-label="Data split">
          <option value="all">All splits</option>
          {SPLIT_ORDER.map((s) => (
            <option key={s} value={s}>
              {SPLIT_LABEL[s]}
            </option>
          ))}
        </select>
      </div>
      <div className="count muted">{shown.length} of {index.length} chains</div>
      <ul className="antigen-list">
        {shown.map((a) => (
          <li key={a.id}>
            <button className={a.id === selectedId ? "antigen active" : "antigen"} onClick={() => onSelect(a.id)}>
              <span className={`dot ${isHeldOut(a.split) ? "held" : ""}`} title={SPLIT_LABEL[a.split]} />
              <span className="antigen-id">{a.pdb} · {a.chain}</span>
              <span className="muted">{a.subtype} · {a.year}</span>
              <span className="antigen-score">{a.auprc === null ? "" : a.auprc.toFixed(2)}</span>
            </button>
          </li>
        ))}
      </ul>
    </aside>
  );
}
