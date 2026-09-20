import { useCallback, useEffect, useMemo, useState } from "react";
import type { Antigen, AntigenSummary, MetricRow, Phylogeny, UploadResult } from "./types";
import { loadAntigen, loadIndex, loadMetrics, loadPhylogeny } from "./data";
import { findPatches } from "./patches";
import { AntigenPicker } from "./components/AntigenPicker";
import { StructureViewer, type ViewerOptions } from "./components/StructureViewer";
import { ScoreTrack } from "./components/ScoreTrack";
import { ResidueCard } from "./components/ResidueCard";
import { HotspotTable } from "./components/HotspotTable";
import { RegionSummary } from "./components/RegionSummary";
import { DurabilityPlot } from "./components/DurabilityPlot";
import { TrustBar } from "./components/TrustBar";
import { Mechanism } from "./components/Mechanism";
import { UploadPanel } from "./components/UploadPanel";
import { UPLOAD_ENABLED } from "./api";

type Tab = "hotspots" | "regions" | "durability";
const TABS: [Tab, string][] = [
  ["hotspots", "Patches"],
  ["regions", "Regions"],
  ["durability", "Durability"],
];

// A full-length H3 from the held-out branch with 26 observed contacts. Chains
// with only a handful of observed residues produce an AUPRC that reads as a
// failing grade when it is really just an undefined one, which is a poor thing
// to land on; this sits just below the median (0.46 against 0.55) for held-out
// chains with enough coverage to score, so it is representative, not flattering.
const DEFAULT_ID = "5kaq_C";

export function App() {
  const [index, setIndex] = useState<AntigenSummary[] | null>(null);
  const [metrics, setMetrics] = useState<MetricRow[]>([]);
  const [phylogeny, setPhylogeny] = useState<Phylogeny | null>(null);
  const [upload, setUpload] = useState<UploadResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  const [selectedId, setSelectedId] = useState(() => window.location.hash.slice(1));
  const [antigen, setAntigen] = useState<Antigen | null>(null);

  const [options, setOptions] = useState<ViewerOptions>({
    showSites: true,
    showTruth: false,
    showAntibodies: true,
    showOtherCopies: true,
    showGlycans: true,
    showSurface: false,
  });
  const [selected, setSelected] = useState<number[]>([]);
  const [hovered, setHovered] = useState<number | null>(null);
  const [focus, setFocus] = useState<{ indices: number[]; nonce: number } | null>(null);
  const [ca, setCa] = useState<Float32Array | null>(null);
  const [tab, setTab] = useState<Tab>("hotspots");

  useEffect(() => {
    loadPhylogeny().then(setPhylogeny).catch(() => setPhylogeny(null));
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
    // an uploaded structure lives in memory, not in the dataset
    if (!selectedId || selectedId.startsWith("upload:")) return;
    let cancelled = false;
    window.history.replaceState(null, "", `#${selectedId}`);
    loadAntigen(selectedId)
      .then((a) => {
        if (cancelled) return;
        // Swap in one go. Clearing selection up front would blank the panels
        // while the next antigen is still in flight, and the page would jump.
        setSelected([]);
        setHovered(null);
        setFocus(null);
        setCa(null);
        setAntigen(a);
      })
      .catch((e: Error) => !cancelled && setError(e.message));
    return () => {
      cancelled = true;
    };
  }, [selectedId]);

  // Keep the loaded antigen on screen while the next one fetches. Dropping to a
  // placeholder unmounts the viewer, which tears down its WebGL context and
  // collapses the page height, so every switch reads as a full reload.
  const uploadChain = selectedId.startsWith("upload:") ? selectedId.slice("upload:".length) : null;
  const uploaded = upload && uploadChain ? (upload.antigens.find((a) => a.chain === uploadChain) ?? null) : null;
  const current = uploaded ?? antigen;
  const pending = !uploaded && !uploadChain && Boolean(selectedId) && antigen?.id !== selectedId;
  const summary: AntigenSummary | null = uploaded
    ? {
        id: uploaded.id, pdb: uploaded.pdb, chain: uploaded.chain, subtype: uploaded.subtype, year: uploaded.year,
        phylo: uploaded.phylo, split: "upload", chainType: uploaded.chainType, nResidues: uploaded.num.length,
        nEpitope: 0, auprc: null,
      }
    : (index?.find((a) => a.id === (current?.id ?? selectedId)) ?? null);

  // stable identity: the viewer reloads its structure whenever this changes
  const antigenChains = useMemo(
    () =>
      upload && uploaded
        ? upload.antigens.map((a) => a.chain)
        : index && current
          ? index.filter((a) => a.pdb === current.pdb).map((a) => a.chain)
          : [],
    [index, current, upload, uploaded],
  );

  const resetSelection = useCallback(() => {
    setSelected([]);
    setHovered(null);
    setFocus(null);
    setCa(null);
  }, []);
  const showUpload = useCallback(
    (result: UploadResult) => {
      resetSelection();
      setUpload(result);
      setSelectedId(`upload:${result.antigens[0].chain}`);
    },
    [resetSelection],
  );
  const pickUploadChain = useCallback(
    (chain: string) => {
      resetSelection();
      setSelectedId(`upload:${chain}`);
    },
    [resetSelection],
  );
  const clearUpload = useCallback(() => {
    setUpload(null);
    if (uploadChain && index) {
      resetSelection();
      setSelectedId((index.find((r) => r.id === DEFAULT_ID) ?? index[0]).id);
    }
  }, [uploadChain, index, resetSelection]);
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
        <h1>Unable to load data</h1>
        <p>{error}</p>
        <p className="muted">
          Run <code>python scripts/08_export_web_data.py</code> from the repo root to generate <code>web/public/data</code>.
        </p>
      </div>
    );
  }

  return (
    <div className="app">
      <header className="masthead">
        <div>
          <p className="eyebrow">Influenza hemagglutinin · Preclinical target selection</p>
          <h1>
            Epitope prioritization for <em>durable</em> antibody targets
          </h1>
          <p className="standfirst">
            Residue-level prediction of antibody-binding sites, combined with cross-strain conservation to rank epitopes
            that are unlikely to be lost to viral escape. Trained on solved antibody–antigen structures and validated on
            held-out influenza subtypes.
          </p>
        </div>
        <div className="key">
          <span><i className="key-swatch" style={{ background: "var(--antibody)" }} /> bound antibody</span>
          <span><i className="key-swatch" style={{ background: "var(--ghost)" }} /> other protomers</span>
          {!uploaded && <span><i className="key-swatch" style={{ background: "var(--truth)" }} /> observed contact</span>}
          <span><i className="key-swatch" style={{ background: "var(--select)" }} /> selected</span>
          <span><i className="key-swatch" style={{ background: "var(--glycan)" }} /> glycan</span>
        </div>
      </header>

      <div className="layout">
        {index ? (
          <AntigenPicker
            index={index}
            phylogeny={phylogeny}
            selectedId={selectedId}
            onSelect={setSelectedId}
            highlightSubtype={uploaded?.subtype}
            top={
              UPLOAD_ENABLED ? (
                <UploadPanel
                  upload={upload}
                  activeChain={uploadChain}
                  onResult={showUpload}
                  onPickChain={pickUploadChain}
                  onClear={clearUpload}
                />
              ) : undefined
            }
          />
        ) : (
          <aside className="picker muted">Loading antigen index…</aside>
        )}

        <main className={pending ? "main is-switching" : "main"} aria-busy={pending}>
          {summary && (
            <TrustBar
              antigen={summary}
              metrics={metrics}
              applicability={uploaded ? upload?.applicability : undefined}
              warnings={uploaded ? upload?.warnings : undefined}
            />
          )}

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
                      showObserved={!uploaded}
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
            <div className="loading-main muted">Loading structure…</div>
          )}

          <Mechanism metrics={metrics} />
        </main>
      </div>

      <footer className="page-foot">
        <span>Data: PDB via SAbDab</span>
        <span>Language model: ESM-2 (650M)</span>
        <span>Classifier: gradient-boosted trees (XGBoost)</span>
        <span>Research prototype for target prioritization; not for clinical decision-making</span>
      </footer>
    </div>
  );
}
