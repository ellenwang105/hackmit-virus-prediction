import { useEffect, useRef, useState } from "react";
import * as $3Dmol from "3dmol";
import type { Antigen } from "../types";
import { loadStructure } from "../data";
import {
  GHOST_COLOR,
  GLYCAN_COLOR,
  LEGEND_TICKS,
  SELECT_COLOR,
  TRUTH_COLOR,
  hex,
  legendGradient,
  scoreInt,
} from "../color";

export interface ViewerOptions {
  showTruth: boolean;
  showAntibodies: boolean;
  showGlycans: boolean;
  showSurface: boolean;
}

interface Props {
  antigen: Antigen;
  /** every antigen chain of this PDB entry; anything else in the file is treated as antibody */
  antigenChains: string[];
  options: ViewerOptions;
  onOptions: (next: ViewerOptions) => void;
  selected: number[];
  hovered: number | null;
  /** bump `nonce` to zoom the camera onto `indices` */
  focus: { indices: number[]; nonce: number } | null;
  onHover: (i: number | null) => void;
  onClick: (i: number, additive: boolean) => void;
  onCoords: (ca: Float32Array | null) => void;
}

interface Loaded {
  /** 3Dmol atom serial -> index into the antigen's arrays */
  atomIndex: Map<number, number>;
  /** CA coordinates, 3 floats per residue, NaN where missing */
  ca: Float32Array;
  antibodyChains: string[];
  otherAntigenChains: string[];
}

const GLYCAN_RESN = ["NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "SIA", "GLC", "XYS"];

/**
 * Map 3Dmol atoms onto the antigen's residue arrays. 3Dmol keeps the PDB
 * residue number but drops the insertion code, so 52 and 52A look identical.
 * Residues are told apart by walking the chain: a backbone N more than 2 A
 * from the current residue's first N starts a new residue (a closer one is
 * just an alternate conformation of the same residue).
 */
function mapResidues(model: any, antigen: Antigen, allChains: string[], antigenChains: string[]): Loaded {
  const numToIndices = new Map<number, number[]>();
  antigen.num.forEach((n, i) => numToIndices.set(n, [...(numToIndices.get(n) ?? []), i]));

  const atomIndex = new Map<number, number>();
  const ca = new Float32Array(antigen.num.length * 3).fill(NaN);
  const seen = new Map<number, number>();
  let currentResi = NaN;
  let currentIndex = -1;
  let firstN: [number, number, number] | null = null;

  for (const atom of model.selectedAtoms({ chain: antigen.chain, hetflag: false })) {
    let startsResidue = atom.resi !== currentResi;
    if (!startsResidue && atom.atom === "N") {
      if (!firstN) firstN = [atom.x, atom.y, atom.z];
      else if (Math.hypot(atom.x - firstN[0], atom.y - firstN[1], atom.z - firstN[2]) > 2) {
        startsResidue = true;
      }
    }
    if (startsResidue) {
      currentResi = atom.resi;
      const occurrence = seen.get(currentResi) ?? 0;
      seen.set(currentResi, occurrence + 1);
      currentIndex = numToIndices.get(currentResi)?.[occurrence] ?? -1;
      firstN = atom.atom === "N" ? [atom.x, atom.y, atom.z] : null;
    }
    if (currentIndex < 0) continue;
    atomIndex.set(atom.serial, currentIndex);
    if (atom.atom === "CA" && Number.isNaN(ca[3 * currentIndex])) {
      ca.set([atom.x, atom.y, atom.z], 3 * currentIndex);
    }
  }

  const placed = Array.from({ length: antigen.num.length }, (_, i) => !Number.isNaN(ca[3 * i])).filter(Boolean).length;
  if (placed < 0.9 * antigen.num.length) {
    console.warn(
      `[epitope-explorer] ${antigen.id}: only ${placed} of ${antigen.num.length} residues matched the structure; ` +
        "colours may be misplaced (residue numbering mismatch?)",
    );
  }

  return {
    atomIndex,
    ca,
    antibodyChains: allChains.filter((c) => !antigenChains.includes(c)),
    otherAntigenChains: antigenChains.filter((c) => c !== antigen.chain),
  };
}

function applyStyles(
  viewer: any,
  loaded: Loaded,
  antigen: Antigen,
  options: ViewerOptions,
  selected: number[],
) {
  const own = { chain: antigen.chain, hetflag: false };
  const indexOf = (atom: any) => loaded.atomIndex.get(atom.serial) ?? -1;
  const colourOf = (atom: any) => {
    const i = indexOf(atom);
    return i < 0 ? hex(GHOST_COLOR) : scoreInt(antigen.score[i]);
  };

  viewer.setStyle({}, {});
  viewer.removeAllSurfaces();

  if (options.showAntibodies && loaded.antibodyChains.length) {
    viewer.setStyle({ chain: loaded.antibodyChains, hetflag: false }, {
      cartoon: { color: "#8a94a6", opacity: 0.35 },
    });
  }
  if (loaded.otherAntigenChains.length) {
    viewer.setStyle({ chain: loaded.otherAntigenChains, hetflag: false }, {
      cartoon: { color: GHOST_COLOR, opacity: 0.6 },
    });
  }
  viewer.setStyle(own, { cartoon: { colorfunc: colourOf } });

  if (options.showGlycans) {
    viewer.setStyle({ hetflag: true, resn: GLYCAN_RESN }, {
      stick: { radius: 0.14, color: hex(GLYCAN_COLOR) },
    });
  }
  if (options.showTruth) {
    const truth = new Set(antigen.epitope.flatMap((e, i) => (e ? [i] : [])));
    viewer.addStyle(
      { ...own, predicate: (a: any) => truth.has(indexOf(a)) },
      { stick: { radius: 0.15, color: hex(TRUTH_COLOR) } },
    );
  }
  if (selected.length) {
    const chosen = new Set(selected);
    viewer.addStyle(
      { ...own, predicate: (a: any) => chosen.has(indexOf(a)) },
      { stick: { radius: 0.28, color: hex(SELECT_COLOR) } },
    );
  }
  if (options.showSurface) {
    viewer.addSurface($3Dmol.SurfaceType.VDW, { colorfunc: colourOf, opacity: 0.92 }, own);
  }
  viewer.render();
}

export function StructureViewer(props: Props) {
  const { antigen, antigenChains, options, selected, hovered, focus } = props;
  const container = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const loadedRef = useRef<Loaded | null>(null);
  const markerRef = useRef<any>(null);
  const handlers = useRef(props);
  handlers.current = props;
  const [status, setStatus] = useState<"loading" | "ready" | "error">("loading");
  const [message, setMessage] = useState("");

  // one WebGL viewer for the life of the component
  useEffect(() => {
    const node = container.current!;
    const viewer = $3Dmol.createViewer(node, { backgroundAlpha: 0 });
    viewerRef.current = viewer;
    const observer = new ResizeObserver(() => viewer.resize());
    observer.observe(node);
    return () => {
      observer.disconnect();
      viewer.clear();
      node.replaceChildren();
      viewerRef.current = null;
    };
  }, []);

  // load the structure whenever the antigen changes
  useEffect(() => {
    let cancelled = false;
    const viewer = viewerRef.current;
    setStatus("loading");
    handlers.current.onCoords(null);
    loadedRef.current = null;

    loadStructure(antigen.pdb)
      .then((text) => {
        if (cancelled || !viewer) return;
        viewer.removeAllModels();
        viewer.removeAllShapes();
        markerRef.current = null;
        const model = viewer.addModel(text, "cif");
        const allChains = [...new Set<string>(model.selectedAtoms({}).map((a: any) => a.chain))];
        const loaded = mapResidues(model, antigen, allChains, antigenChains);
        loadedRef.current = loaded;

        const own = { chain: antigen.chain, hetflag: false };
        const index = (atom: any) => loaded.atomIndex.get(atom.serial);
        viewer.setHoverable(
          own,
          true,
          (atom: any) => {
            const i = index(atom);
            if (i !== undefined) handlers.current.onHover(i);
          },
          () => handlers.current.onHover(null),
        );
        viewer.setClickable(own, true, (atom: any, _v: unknown, event: MouseEvent) => {
          const i = index(atom);
          if (i !== undefined) handlers.current.onClick(i, event.shiftKey || event.ctrlKey || event.metaKey);
        });

        applyStyles(viewer, loaded, antigen, handlers.current.options, handlers.current.selected);
        viewer.zoomTo(own);
        viewer.zoom(1.3, 0);
        viewer.render();
        handlers.current.onCoords(loaded.ca);
        setStatus("ready");
      })
      .catch((error: Error) => {
        if (cancelled) return;
        setMessage(error.message);
        setStatus("error");
      });
    return () => {
      cancelled = true;
    };
  }, [antigen, antigenChains]);

  // restyle when what is shown changes
  useEffect(() => {
    const viewer = viewerRef.current;
    const loaded = loadedRef.current;
    if (status === "ready" && viewer && loaded) applyStyles(viewer, loaded, antigen, options, selected);
  }, [status, antigen, options, selected]);

  // hover marker: a translucent sphere is far cheaper than restyling the cartoon
  useEffect(() => {
    const viewer = viewerRef.current;
    const loaded = loadedRef.current;
    if (status !== "ready" || !viewer || !loaded) return;
    if (markerRef.current) viewer.removeShape(markerRef.current);
    markerRef.current = null;
    if (hovered !== null && !Number.isNaN(loaded.ca[3 * hovered])) {
      markerRef.current = viewer.addSphere({
        center: { x: loaded.ca[3 * hovered], y: loaded.ca[3 * hovered + 1], z: loaded.ca[3 * hovered + 2] },
        radius: 1.9,
        color: SELECT_COLOR,
        alpha: 0.55,
      });
    }
    viewer.render();
  }, [status, hovered]);

  useEffect(() => {
    const viewer = viewerRef.current;
    const loaded = loadedRef.current;
    if (!focus || status !== "ready" || !viewer || !loaded || !focus.indices.length) return;
    const wanted = new Set(focus.indices);
    viewer.zoomTo(
      { chain: antigen.chain, hetflag: false, predicate: (a: any) => wanted.has(loaded.atomIndex.get(a.serial) ?? -1) },
      600,
    );
    // only a new nonce should move the camera, not later status changes
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [focus?.nonce]);

  const toggle = (key: keyof ViewerOptions) => props.onOptions({ ...options, [key]: !options[key] });

  return (
    <div className="viewer">
      <div ref={container} className="viewer-canvas" />
      <div className="viewer-toolbar">
        {(
          [
            ["showTruth", "Observed epitope"],
            ["showSurface", "Surface"],
            ["showAntibodies", "Antibodies"],
            ["showGlycans", "Glycans"],
          ] as const
        ).map(([key, label]) => (
          <label key={key} className="check">
            <input type="checkbox" checked={options[key]} onChange={() => toggle(key)} />
            {label}
          </label>
        ))}
        <button
          className="ghost-button"
          onClick={() => viewerRef.current?.zoomTo({ chain: antigen.chain, hetflag: false }, 500)}
        >
          Reset view
        </button>
      </div>
      {status !== "ready" && (
        <div className="viewer-status" role="status">
          {status === "loading" ? "Loading structure…" : `Could not load structure: ${message}`}
        </div>
      )}
      <div className="viewer-legend">
        <span className="legend-title">Predicted probability</span>
        <div className="legend-bar" style={{ background: legendGradient() }} />
        <div className="legend-ticks">
          {LEGEND_TICKS.map((t) => (
            <span key={t} style={{ left: `${Math.sqrt(t) * 100}%` }}>
              {t}
            </span>
          ))}
        </div>
      </div>
    </div>
  );
}
