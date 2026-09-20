import { useMemo, useState, type ReactNode } from "react";
import type { AntigenSummary, Phylogeny } from "../types";
import { isHeldOut } from "../data";
import { PhyloTree } from "./PhyloTree";

interface Props {
  index: AntigenSummary[];
  phylogeny: Phylogeny | null;
  selectedId: string;
  onSelect: (id: string) => void;
  /** rendered above the tree, for the upload panel */
  top?: ReactNode;
  /** subtype to light up in the tree when the antigen on screen is not in the index */
  highlightSubtype?: string;
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

export function AntigenPicker({ index, phylogeny, selectedId, onSelect, top, highlightSubtype }: Props) {
  const [query, setQuery] = useState("");
  const [subtype, setSubtype] = useState("all");

  // the tree covers the subtypes with a coherent reference sequence; the rest
  // (unknown, chimeric constructs) stay reachable from this list
  const strays = useMemo(
    () =>
      phylogeny
        ? [...new Set(index.map((a) => a.subtype))].filter((s) => !(s in phylogeny.subtypes)).sort()
        : [],
    [index, phylogeny],
  );

  // One row per solved structure. The same protein is often solved as several chains,
  // and a list of 684 chain letters means nothing to someone browsing; chains appear
  // under the structure that is open.
  const structures = useMemo(() => {
    const byPdb = new Map<string, AntigenSummary[]>();
    for (const a of index) byPdb.set(a.pdb, [...(byPdb.get(a.pdb) ?? []), a]);
    return [...byPdb.values()].map((chains) => {
      const lead = chains.reduce((best, c) => (c.nEpitope > best.nEpitope ? c : best));
      return { pdb: lead.pdb, subtype: lead.subtype, year: lead.year, split: lead.split, lead, chains };
    });
  }, [index]);

  const selectedPdb = index.find((a) => a.id === selectedId)?.pdb ?? null;

  // held-out structures first: they are the ones whose scores are not flattered by training
  const shown = useMemo(() => {
    const q = query.trim().toLowerCase();
    return structures
      .filter((s) => (subtype === "all" || s.subtype === subtype) && (!q || s.pdb.toLowerCase().includes(q)))
      .sort((a, b) => SPLIT_ORDER.indexOf(a.split) - SPLIT_ORDER.indexOf(b.split) || a.pdb.localeCompare(b.pdb));
  }, [structures, query, subtype]);

  return (
    <aside className="picker">
      {top}
      {phylogeny && (
        <PhyloTree
          phylogeny={phylogeny}
          value={index.find((a) => a.id === selectedId)?.subtype ?? highlightSubtype ?? "all"}
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

      <h2>Structures</h2>
      <input
        className="input"
        type="search"
        placeholder="Search by PDB ID"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
        aria-label="Search structures by PDB ID"
      />
      {strays.length > 0 && (
        <div className="filters">
          <select className="input" value={subtype} onChange={(e) => setSubtype(e.target.value)} aria-label="Other virus types">
            <option value="all">Other virus types</option>
            {strays.map((s) => (
              <option key={s} value={s}>
                {s}
              </option>
            ))}
          </select>
        </div>
      )}
      <div className="count muted">{shown.length} structures</div>
      <ul className="antigen-list">
        {shown.map((s) => {
          const active = s.pdb === selectedPdb;
          return (
            <li key={s.pdb}>
              <button
                className={active ? "antigen active" : "antigen"}
                aria-current={active}
                onClick={() => !active && onSelect(s.lead.id)}
              >
                <span className="antigen-id">{s.pdb.toUpperCase()}</span>
                <span className="muted">{s.subtype} · {s.year}</span>
              </button>
              {active && s.chains.length > 1 && (
                <div className="antigen-chains" role="group" aria-label="Chains in this structure">
                  <span className="muted">Chain</span>
                  {s.chains.map((c) => (
                    <button key={c.id} className="chip" aria-pressed={c.id === selectedId} onClick={() => onSelect(c.id)}>
                      {c.chain}
                    </button>
                  ))}
                </div>
              )}
            </li>
          );
        })}
      </ul>
    </aside>
  );
}
