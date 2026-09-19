import { useMemo, useState } from "react";
import type { AntigenSummary } from "../types";
import { SPLIT_LABEL, isHeldOut } from "../data";

interface Props {
  index: AntigenSummary[];
  selectedId: string;
  onSelect: (id: string) => void;
}

const SPLIT_ORDER = ["test_group2", "test_B", "val", "train", "excluded"];

export function AntigenPicker({ index, selectedId, onSelect }: Props) {
  const [query, setQuery] = useState("");
  const [subtype, setSubtype] = useState("all");
  const [split, setSplit] = useState("all");

  const subtypes = useMemo(() => [...new Set(index.map((a) => a.subtype))].sort(), [index]);
  const shown = useMemo(() => {
    const q = query.trim().toLowerCase();
    return index
      .filter((a) => (subtype === "all" || a.subtype === subtype) && (split === "all" || a.split === split))
      .filter((a) => !q || a.id.toLowerCase().includes(q))
      .sort((a, b) => SPLIT_ORDER.indexOf(a.split) - SPLIT_ORDER.indexOf(b.split) || a.id.localeCompare(b.id));
  }, [index, query, subtype, split]);

  return (
    <aside className="picker">
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
        <select className="input" value={subtype} onChange={(e) => setSubtype(e.target.value)} aria-label="Subtype">
          <option value="all">All subtypes</option>
          {subtypes.map((s) => (
            <option key={s} value={s}>
              {s}
            </option>
          ))}
        </select>
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
