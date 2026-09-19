import { useEffect, useMemo, useRef, useState } from "react";
import type { AntigenSummary } from "../types";
import { SPLIT_LABEL, isHeldOut, loadSequences } from "../data";
import {
  MIN_SCORE,
  exampleQuery,
  parseQuery,
  searchSequences,
  type SearchResult,
  type SequenceHit,
} from "../search";
import "../search.css";

interface Props {
  open: boolean;
  onClose: () => void;
  index: AntigenSummary[];
  selectedId: string;
  /** chain whose sequence seeds the "try an example" button */
  exampleId: string;
  onPick: (hit: SequenceHit, result: SearchResult) => void;
  /** reported so the page can show the last search's matches next to the structure */
  onResult: (result: SearchResult) => void;
}

const PAGE = 10;
const pct = (fraction: number) => `${(fraction * 100).toFixed(fraction === 1 ? 0 : 1)}%`;

function verdict(identity: number, coverage: number): { label: string; tone: "good" | "neutral" | "warn" } {
  // a perfect match to a small piece of the query is not a close match to the query
  if (coverage < 0.6) return { label: "Partial match", tone: "neutral" };
  if (identity >= 0.9) return { label: "Close match", tone: "good" };
  if (identity >= 0.6) return { label: "Related", tone: "neutral" };
  return { label: "Distant", tone: "warn" };
}

export function SequenceSearch({ open, onClose, index, selectedId, exampleId, onPick, onResult }: Props) {
  const dialog = useRef<HTMLDialogElement>(null);
  const request = useRef(0);
  const [text, setText] = useState("");
  const [sequences, setSequences] = useState<Record<string, string> | null>(null);
  const [problem, setProblem] = useState<string | null>(null);
  const [phase, setPhase] = useState<"idle" | "searching" | "done">("idle");
  const [result, setResult] = useState<SearchResult | null>(null);
  const [visible, setVisible] = useState(PAGE);

  const parsed = useMemo(() => parseQuery(text), [text]);
  const canSearch = Boolean(sequences) && parsed.sequence.length > 0 && !parsed.error && phase !== "searching";

  useEffect(() => {
    const node = dialog.current;
    if (!node) return;
    if (open && !node.open) node.showModal();
    if (!open && node.open) node.close();
  }, [open]);

  useEffect(() => {
    if (!open || sequences) return;
    loadSequences()
      .then(setSequences)
      .catch(() =>
        setProblem("Sequence data not found. Run scripts/08_export_web_data.py from the repo root, then reload."),
      );
  }, [open, sequences]);

  const search = async () => {
    if (!canSearch || !sequences) return;
    const mine = ++request.current;
    setPhase("searching");
    setProblem(null);
    try {
      const found = await searchSequences(parsed, sequences, index);
      if (mine !== request.current) return;
      setResult(found);
      setVisible(PAGE);
      setPhase("done");
      onResult(found);
    } catch (error) {
      if (mine !== request.current) return;
      setProblem(`Search failed: ${(error as Error).message}`);
      setPhase("idle");
    }
  };

  const fillExample = () => {
    if (!sequences) return;
    const example = exampleQuery(sequences, exampleId);
    if (example) setText(example.text);
  };

  const best = result?.hits[0];
  const shown = result?.hits.slice(0, visible) ?? [];

  return (
    <dialog
      ref={dialog}
      className="seq-dialog"
      aria-labelledby="seq-title"
      onClose={onClose}
      onClick={(event) => event.target === event.currentTarget && onClose()}
    >
      <div className="seq-body">
        <header className="seq-head">
          <div>
            <h2 id="seq-title">Find the closest structures</h2>
            <p className="muted">
              Paste a hemagglutinin protein sequence. It is aligned against every antigen in the data, and you can open
              the best matches.
            </p>
          </div>
          <button className="seq-close" onClick={onClose} aria-label="Close">
            ×
          </button>
        </header>

        <label className="seq-label" htmlFor="seq-input">
          Amino-acid sequence, raw or FASTA
        </label>
        <textarea
          id="seq-input"
          className="seq-input"
          spellCheck={false}
          autoComplete="off"
          placeholder={">my strain\nMKTIIALSYIFCLALG…"}
          value={text}
          onChange={(event) => setText(event.target.value)}
          onKeyDown={(event) => {
            if (event.key === "Enter" && (event.ctrlKey || event.metaKey)) search();
          }}
        />

        <div className="seq-actions">
          <button className="primary-button" onClick={search} disabled={!canSearch}>
            {phase === "searching" ? "Searching…" : "Search"}
          </button>
          <button className="ghost-button" onClick={fillExample} disabled={!sequences}>
            Try an example
          </button>
          {text && (
            <button
              className="ghost-button"
              onClick={() => {
                setText("");
                setResult(null);
                setPhase("idle");
              }}
            >
              Clear
            </button>
          )}
          <span className="seq-status" role="status">
            {parsed.error ? (
              <span className="seq-error">{parsed.error}</span>
            ) : parsed.sequence ? (
              <>
                {parsed.sequence.length} residues{parsed.name ? ` · ${parsed.name}` : ""}
              </>
            ) : (
              <span className="muted">Ctrl + Enter to search</span>
            )}
          </span>
        </div>
        {parsed.warnings.map((warning) => (
          <p key={warning} className="seq-warning">
            {warning}
          </p>
        ))}
        {problem && <p className="seq-error">{problem}</p>}

        {result && (
          <section className="seq-results" aria-live="polite">
            {best ? (
              <>
                <div className={`seq-summary ${verdict(best.identity, best.queryCoverage).tone}`}>
                  <strong>{verdict(best.identity, best.queryCoverage).label}</strong>
                  <span>
                    The best structure has {best.matches} of your {result.query.sequence.length} residues identical (
                    {pct(best.match)}).
                  </span>
                  {best.identity < 0.6 && (
                    <span>
                      That is distant, so predictions on these structures may not carry over to your sequence.
                    </span>
                  )}
                </div>

                <ol className="matches">
                  {shown.map((hit, rank) => (
                    <li key={`${hit.pdb}-${hit.id}`}>
                      <button
                        className={hit.chains.includes(selectedId) ? "match active" : "match"}
                        onClick={() => onPick(hit, result)}
                      >
                        <span className="match-rank">{rank + 1}</span>
                        <span className="match-main">
                          <span className="match-title">
                            {hit.pdb} · {hit.chain}
                            <span className="chip">{hit.chainType}</span>
                            {hit.chains.length > 1 && (
                              <span className="chip" title={hit.chains.join(", ")}>
                                +{hit.chains.length - 1} identical chain{hit.chains.length > 2 ? "s" : ""}
                              </span>
                            )}
                          </span>
                          <span className="match-sub">
                            {hit.subtype} · {hit.year} ·{" "}
                            <span className={isHeldOut(hit.split) ? "held" : ""}>{SPLIT_LABEL[hit.split]}</span>
                            {hit.nEpitope > 0 && <> · {hit.nEpitope} observed contacts</>}
                          </span>
                        </span>
                        <span className="match-stats">
                          <span className="identity" title="Share of your sequence found identically in this structure">
                            {pct(hit.match)}
                          </span>
                          <span className="identity-bar" aria-hidden>
                            <i style={{ width: `${hit.match * 100}%` }} />
                          </span>
                          <span
                            className="match-sub"
                            title={`${pct(hit.identity)} identical across the aligned stretch`}
                          >
                            {hit.matches} of your {result.query.sequence.length} residues identical · aligns{" "}
                            {Math.round(hit.queryCoverage * 100)}% of yours
                          </span>
                        </span>
                      </button>
                    </li>
                  ))}
                </ol>

                {result.hits.length > visible && (
                  <button className="ghost-button more" onClick={() => setVisible(result.hits.length)}>
                    Show {result.hits.length - visible} more
                  </button>
                )}
                <p className="seq-note">
                  Sorted by match, highest first: the share of your sequence found identically in each structure, so
                  100% means every residue of yours. Scores on a match describe that structure, not your
                  sequence: positions where yours differs are not re-predicted. Runs in your browser; nothing is
                  uploaded.
                </p>
              </>
            ) : (
              <div className="seq-summary warn">
                <strong>No significant match</strong>
                <span>
                  The best alignment scored {result.bestScore}; a real match needs at least {MIN_SCORE}. This may not be
                  influenza hemagglutinin.
                </span>
              </div>
            )}
          </section>
        )}
      </div>
    </dialog>
  );
}
