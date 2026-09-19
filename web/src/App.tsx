import { useCallback, useEffect, useMemo, useState } from "react";
import type { Antigen, AntigenSummary, MetricRow } from "./types";
import { loadAntigen, loadIndex, loadMetrics } from "./data";
import { findPatches } from "./patches";
import { AntigenPicker } from "./components/AntigenPicker";
import { StructureViewer, type ViewerOptions } from "./components/StructureViewer";
import { ScoreTrack } from "./components/ScoreTrack";
import { ResidueCard } from "./components/ResidueCard";
import { HotspotTable } from "./components/HotspotTable";
import { RegionSummary } from "./components/RegionSummary";
import { DurabilityPlot } from "./components/DurabilityPlot";
import { TrustBar } from "./components/TrustBar";

type Tab = "hotspots" | "regions" | "durability";
const TABS: [Tab, string][] = [
  ["hotspots", "Hotspots"],
  ["regions", "Regions"],
  ["durability", "Durability"],
];

const DEFAULT_ID = "3sdy_A";

export function App() {
  const [index, setIndex] = useState<AntigenSummary[] | null>(null);
  const [metrics, setMetrics] = useState<MetricRow[]>([]);
  const [error, setError] = useState<string | null>(null);

  const [selectedId, setSelectedId] = useState(() => window.location.hash.slice(1));
  const [antigen, setAntigen] = useState<Antigen | null>(null);

  const [options, setOptions] = useState<ViewerOptions>({
    showTruth: false,
    showAntibodies: true,
    showGlycans: true,
    showSurface: false,
  });
  const [selected, setSelected] = useState<number[]>([]);
  const [hovered, setHovered] = useState<number | null>(null);
  const [focus, setFocus] = useState<{ indices: number[]; nonce: number } | null>(null);
  const [ca, setCa] = useState<Float32Array | null>(null);
  const [tab, setTab] = useState<Tab>("hotspots");

  useEffect(() => {
    Promise.all([loadIndex(), loadMetrics()])
      .then(([rows, metricRows]) => {
        setIndex(rows);
        setMetrics(metricRows);
        setSelectedId((current) =>
          rows.some((r) => r.id === current)
            ? current
            : (rows.find((r) => r.id === DEFAULT_ID) ?? rows.find((r) => r.split === "test_group2") ?? rows[0]).id,
        );
      })
      .catch((e: Error) => setError(e.message));
  }, []);

  useEffect(() => {
    if (!selectedId) return;
    let cancelled = false;
    window.history.replaceState(null, "", `#${selectedId}`);
    setSelected([]);
    setHovered(null);
    setFocus(null);
    setCa(null);
    loadAntigen(selectedId)
      .then((a) => !cancelled && setAntigen(a))
      .catch((e: Error) => !cancelled && setError(e.message));
    return () => {
      cancelled = true;
    };
  }, [selectedId]);

  const summary = index?.find((a) => a.id === selectedId) ?? null;
  const current = antigen && antigen.id === selectedId ? antigen : null;

  // stable identity: the viewer reloads its structure whenever this changes
  const antigenChains = useMemo(
    () => (index && current ? index.filter((a) => a.pdb === current.pdb).map((a) => a.chain) : []),
    [index, current],
  );
  const patches = useMemo(() => (current ? findPatches(current, ca) : []), [current, ca]);

  const select = useCallback((indices: number[], additive: boolean) => {
    setSelected((prev) => {
      if (additive) {
        if (indices.length === 1) {
          return prev.includes(indices[0]) ? prev.filter((i) => i !== indices[0]) : [...prev, indices[0]];
        }
        return [...new Set([...prev, ...indices])];
      }
      // clicking the lone selected residue again deselects it
      if (indices.length === 1 && prev.length === 1 && prev[0] === indices[0]) return [];
      return indices;
    });
  }, []);

  const pickPatch = useCallback((indices: number[]) => {
    setSelected(indices);
    setFocus((f) => ({ indices, nonce: (f?.nonce ?? 0) + 1 }));
  }, []);

  if (error) {
    return (
      <div className="fatal">
        <h1>Couldn't load data</h1>
        <p>{error}</p>
        <p className="muted">
          Run <code>python scripts/08_export_web_data.py</code> from the repo root to generate <code>web/public/data</code>.
        </p>
      </div>
    );
  }

  return (
    <div className="app">
      <header className="header">
        <div>
          <h1>Epitope Explorer</h1>
          <p className="muted">Where on influenza hemagglutinin antibodies are likely to bind</p>
        </div>
        <div className="key">
          <span><i className="key-swatch" style={{ background: "var(--truth)" }} /> observed epitope</span>
          <span><i className="key-swatch" style={{ background: "var(--select)" }} /> selected</span>
          <span><i className="key-swatch" style={{ background: "var(--glycan)" }} /> glycan</span>
        </div>
      </header>

      <div className="layout">
        {index ? (
          <AntigenPicker index={index} selectedId={selectedId} onSelect={setSelectedId} />
        ) : (
          <aside className="picker muted">Loading antigens…</aside>
        )}

        <main className="main">
          {summary && <TrustBar antigen={summary} metrics={metrics} />}

          {current ? (
            <>
              <div className="workspace">
                <StructureViewer
                  antigen={current}
                  antigenChains={antigenChains}
                  options={options}
                  onOptions={setOptions}
                  selected={selected}
                  hovered={hovered}
                  focus={focus}
                  onHover={setHovered}
                  onClick={(i, additive) => select([i], additive)}
                  onCoords={setCa}
                />
                <section className="panel">
                  <div className="tabs" role="tablist">
                    {TABS.map(([id, label]) => (
                      <button
                        key={id}
                        role="tab"
                        aria-selected={tab === id}
                        className={tab === id ? "tab active" : "tab"}
                        onClick={() => setTab(id)}
                      >
                        {label}
                      </button>
                    ))}
                  </div>
                  {tab === "hotspots" && (
                    <HotspotTable
                      patches={patches}
                      loading={ca === null}
                      selected={selected}
                      onPick={pickPatch}
                      onClear={() => setSelected([])}
                    />
                  )}
                  {tab === "regions" && <RegionSummary antigen={current} />}
                  {tab === "durability" && (
                    <DurabilityPlot antigen={current} selected={selected} onHover={setHovered} onSelect={select} />
                  )}
                </section>
              </div>

              <ResidueCard antigen={current} index={hovered ?? (selected.length === 1 ? selected[0] : null)} />
              <ScoreTrack antigen={current} hovered={hovered} selected={selected} onHover={setHovered} onSelect={select} />
            </>
          ) : (
            <div className="loading-main muted">Loading antigen…</div>
          )}
        </main>
      </div>
    </div>
  );
}
