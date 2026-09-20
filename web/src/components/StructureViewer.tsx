import { useEffect, useRef, useState } from "react";
import * as $3Dmol from "3dmol";
import type { Antigen } from "../types";
import { isUpload, loadStructure, loadUploadedStructure } from "../data";
import {
  GHOST_COLOR,
  ANTIBODY_COLOR,
  GLYCAN_COLOR,
  LEGEND_TICKS,
  SELECT_COLOR,
  TRUTH_COLOR,
  hex,
  legendGradient,
  scoreInt,
} from "../color";

export interface ViewerOptions {
  /** filled circles on the residues the model scores highest */
  showSites: boolean;
  showTruth: boolean;
  showAntibodies: boolean;
  /** the other protomers of the same trimer, which are HA too, not antibody */
  showOtherCopies: boolean;
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
 * Which checkboxes have anything to act on in the structure that is loaded. A
 * box that silently does nothing looks broken, so the ones with nothing to draw
 * are disabled and say why: a structure with no bound antibody (an upload of an
 * apo HA, say) has no antibody to show, and 31 of the 216 structures model no glycans.
 */
type Available = Record<keyof ViewerOptions, boolean>;

/**
 * What most people want to switch is the prediction, the surface and the antibody.
 * The rest (observed contacts, trimer partners, sugars) are context, and six boxes
 * at once made the toolbar something to decode rather than something to use.
 */
const PRIMARY_LAYERS = [
  ["showSites", "Predicted epitope"],
  ["showSurface", "Surface"],
  ["showAntibodies", "Bound antibody"],
] as const;
const MORE_LAYERS = [
  ["showTruth", "Observed contacts"],
  ["showOtherCopies", "Other protomers"],
  ["showGlycans", "Glycans"],
] as const;
const ALL_AVAILABLE: Available = {
  showSites: true, showTruth: true, showSurface: true, showAntibodies: true, showOtherCopies: true, showGlycans: true,
};
const UNAVAILABLE_REASON: Partial<Record<keyof ViewerOptions, string>> = {
  showAntibodies: "This file contains no antibody chain",
  showOtherCopies: "This file contains no other protomers of this antigen",
  showGlycans: "No modeled glycans in this file",
  showTruth: "No observed antibody contacts on this chain",
};

/**
 * Camera rotation that stands the spike upright with its head at the top.
 *
 * HA is a long thin molecule, and 3Dmol's default camera lands wherever the
 * deposited coordinate frame happens to point — often straight down the long
 * axis, where the structure reads as a shapeless blob. Taking the dominant axis
 * of the alpha carbons and rotating it onto the screen's vertical gives the
 * silhouette everyone recognises, with the drifting head above the stem.
 *
 * Returned as a quaternion in 3Dmol's (x, y, z, w) order, ready for setView.
 */
function uprightQuaternion(ca: Float32Array, region: string[]): [number, number, number, number] | null {
  const points: number[][] = [];
  const isHead: boolean[] = [];
  for (let i = 0; i < region.length; i++) {
    const x = ca[3 * i];
    if (Number.isNaN(x)) continue;
    points.push([x, ca[3 * i + 1], ca[3 * i + 2]]);
    isHead.push(region[i] === "head");
  }
  if (points.length < 20) return null;

  const centre = [0, 1, 2].map((k) => points.reduce((s, p) => s + p[k], 0) / points.length);
  const centred = points.map((p) => [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]]);

  // power iteration converges on the covariance's dominant eigenvector, which
  // is cheaper than pulling in a linear-algebra dependency for one vector
  let axis = [0, 0, 1];
  for (let step = 0; step < 32; step++) {
    const next = [0, 0, 0];
    for (const p of centred) {
      const dot = p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2];
      next[0] += dot * p[0];
      next[1] += dot * p[1];
      next[2] += dot * p[2];
    }
    const norm = Math.hypot(next[0], next[1], next[2]);
    if (norm < 1e-9) return null;
    axis = next.map((v) => v / norm);
  }

  const headPoints = centred.filter((_, i) => isHead[i]);
  if (headPoints.length) {
    const mean = headPoints.reduce((s, p) => s + p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2], 0) / headPoints.length;
    if (mean < 0) axis = axis.map((v) => -v);
  }

  // shortest rotation taking the spike axis onto screen-up
  const up = [0, 1, 0];
  const dot = axis[0] * up[0] + axis[1] * up[1] + axis[2] * up[2];
  if (dot > 0.9999) return [0, 0, 0, 1];
  if (dot < -0.9999) return [0, 0, 1, 0];      // 180 degrees about z
  const cross = [
    axis[1] * up[2] - axis[2] * up[1],
    axis[2] * up[0] - axis[0] * up[2],
    axis[0] * up[1] - axis[1] * up[0],
  ];
  const w = 1 + dot;
  const length = Math.hypot(cross[0], cross[1], cross[2], w);
  return [cross[0] / length, cross[1] / length, cross[2] / length, w / length];
}

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
        "colors may be misplaced (residue numbering mismatch?)",
    );
  }

  return {
    atomIndex,
    ca,
    antibodyChains: allChains.filter((c) => !antigenChains.includes(c)),
    otherAntigenChains: antigenChains.filter((c) => c !== antigen.chain && allChains.includes(c)),
  };
}

/**
 * Predicted sites are drawn as filled circles on the alpha carbons of the
 * top-scoring residues. The cutoff is relative to the chain, not fixed: scores
 * are low on most chains (epitopes are ~5% of residues), and a fixed 0.5 would
 * mark only a handful of residues on a typical held-out antigen.
 */
const SITE_TOP_FRACTION = 0.07;
const SITE_MIN_SCORE = 0.1;
/** circle radius in angstrom for scores below 0.3, 0.3-0.5 and above 0.5 */
const SITE_RADIUS = [1.0, 1.5, 2.0] as const;
const SITE_BREAKS = [0.3, 0.5] as const;
/** ring radius around residues an antibody was actually seen touching */
const HALO_RADIUS = 2.7;

/**
 * Three orthogonal circles round a point: reads as a ring from any angle, and
 * stays sparse. A wireframe sphere did the same job but is a dense mesh, and at
 * this scale it read as a solid blob that hid the circle it was meant to frame.
 */
function ringsAround(viewer: any, centre: { x: number; y: number; z: number }, radius: number) {
  const segments = 28;
  return (["xy", "yz", "xz"] as const).map((plane) => {
    const points = Array.from({ length: segments + 1 }, (_, k) => {
      const angle = (2 * Math.PI * k) / segments;
      const point: Record<string, number> = { ...centre };
      point[plane[0]] += radius * Math.cos(angle);
      point[plane[1]] += radius * Math.sin(angle);
      return point as { x: number; y: number; z: number };
    });
    return viewer.addCurve({ points, radius: 0.13, smooth: 1, fill: false, color: hex(TRUTH_COLOR) });
  });
}

function siteCutoff(score: number[]) {
  const sorted = [...score].sort((a, b) => b - a);
  const rank = Math.max(1, Math.round(score.length * SITE_TOP_FRACTION));
  return Math.max(SITE_MIN_SCORE, sorted[Math.min(rank, sorted.length) - 1] ?? 1);
}

function applyStyles(
  viewer: any,
  loaded: Loaded,
  antigen: Antigen,
  options: ViewerOptions,
  selected: number[],
  halos: { current: any[] },
) {
  const own = { chain: antigen.chain, hetflag: false };
  const indexOf = (atom: any) => loaded.atomIndex.get(atom.serial) ?? -1;
  const colourOf = (atom: any) => {
    const i = indexOf(atom);
    return i < 0 ? hex(GHOST_COLOR) : scoreInt(antigen.score[i]);
  };

  viewer.setStyle({}, {});
  viewer.removeAllSurfaces();

  // Three roles, told apart by hue as well as weight. The other HA protomers
  // used to share a grey with the antibody, so the big pale mass in a trimer
  // read as "the antibody" and the Antibodies checkbox appeared to do nothing.
  if (options.showOtherCopies && loaded.otherAntigenChains.length) {
    viewer.setStyle({ chain: loaded.otherAntigenChains, hetflag: false }, {
      cartoon: { color: GHOST_COLOR, opacity: 0.45 },
    });
  }
  if (options.showAntibodies && loaded.antibodyChains.length) {
    viewer.setStyle({ chain: loaded.antibodyChains, hetflag: false }, {
      cartoon: { color: ANTIBODY_COLOR, opacity: 0.85 },
    });
  }
  viewer.setStyle(own, { cartoon: { colorfunc: colourOf } });

  if (options.showSites) {
    const cutoff = siteCutoff(antigen.score);
    const size = (i: number) => (antigen.score[i] >= SITE_BREAKS[1] ? 2 : antigen.score[i] >= SITE_BREAKS[0] ? 1 : 0);
    // radius is a constant per style, so one style per size bucket
    for (let bucket = 0; bucket < SITE_RADIUS.length; bucket++) {
      viewer.addStyle(
        {
          ...own,
          atom: "CA",
          predicate: (a: any) => {
            const i = indexOf(a);
            return i >= 0 && antigen.score[i] >= cutoff && size(i) === bucket;
          },
        },
        { sphere: { radius: SITE_RADIUS[bucket], colorfunc: colourOf } },
      );
    }
  }

  if (options.showGlycans) {
    viewer.setStyle({ hetflag: true, resn: GLYCAN_RESN }, {
      stick: { radius: 0.22, color: hex(GLYCAN_COLOR) },
    });
  }
  // A ring on every observed contact, whether or not the model marked it, so a
  // ring around a filled circle is a hit and a bare ring is a miss.
  //
  // Shapes rather than atom styles, because 3Dmol keeps one sphere style per atom
  // and merges later ones into it, so a halo added as a style would replace the
  // predicted-site circle instead of surrounding it. Thin rings rather than
  // shells, because alpha shapes did not render here at all and an opaque one
  // buries the circle it is meant to highlight.
  halos.current.forEach((shape) => viewer.removeShape(shape));
  halos.current = [];
  if (options.showTruth) {
    antigen.epitope.forEach((observed, i) => {
      if (!observed || Number.isNaN(loaded.ca[3 * i])) return;
      halos.current.push(
        ...ringsAround(viewer, { x: loaded.ca[3 * i], y: loaded.ca[3 * i + 1], z: loaded.ca[3 * i + 2] }, HALO_RADIUS),
      );
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

/** Breathing room so the whole model clears the frame instead of touching it. */
const FIT_MARGIN = 1.06;

/**
 * Stop the camera backing out past the point where the whole file is in view.
 *
 * `fit` is the camera distance at which 3Dmol frames everything. It sizes a
 * bounding sphere against the vertical field of view, so a canvas taller than
 * it is wide needs proportionally more distance to fit horizontally too.
 */
function limitZoomOut(viewer: any, node: HTMLElement, fit: number | null) {
  if (fit === null) return;
  const aspect = node.clientWidth / Math.max(1, node.clientHeight);
  viewer.setZoomLimits(0, fit * FIT_MARGIN * Math.max(1, 1 / aspect));
}

export function StructureViewer(props: Props) {
  const { antigen, antigenChains, options, selected, hovered, focus } = props;
  const container = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const loadedRef = useRef<Loaded | null>(null);
  const markerRef = useRef<any>(null);
  /** observed-contact shells currently in the scene, so a restyle can remove them */
  const haloRef = useRef<any[]>([]);
  /** camera distance at which the whole file fits; null until a structure loads */
  const fitRef = useRef<number | null>(null);
  const handlers = useRef(props);
  handlers.current = props;
  const [status, setStatus] = useState<"loading" | "ready" | "error">("loading");
  const [available, setAvailable] = useState<Available>(ALL_AVAILABLE);
  const [showMore, setShowMore] = useState(false);
  const [message, setMessage] = useState("");

  // one WebGL viewer for the life of the component
  useEffect(() => {
    const node = container.current!;
    const viewer = $3Dmol.createViewer(node, { backgroundAlpha: 0 });
    viewerRef.current = viewer;
    const observer = new ResizeObserver(() => {
      viewer.resize();
      limitZoomOut(viewer, node, fitRef.current);
    });
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

    (antigen.source ? loadUploadedStructure(antigen.source) : loadStructure(antigen.pdb))
      .then((text) => {
        if (cancelled || !viewer) return;
        viewer.removeAllModels();
        viewer.removeAllShapes();
        markerRef.current = null;
        haloRef.current = [];
        const model = viewer.addModel(text, antigen.source?.format ?? "cif");
        // chains that hold protein; a chain of only sugars or waters is not an antibody
        const allChains = [...new Set<string>(model.selectedAtoms({ hetflag: false }).map((a: any) => a.chain))];
        const loaded = mapResidues(model, antigen, allChains, antigenChains);
        loadedRef.current = loaded;
        setAvailable({
          showSites: true,
          showSurface: true,
          showTruth: antigen.epitope.some((e) => e === 1),
          showAntibodies: loaded.antibodyChains.length > 0,
          showOtherCopies: loaded.otherAntigenChains.length > 0,
          showGlycans: model.selectedAtoms({ hetflag: true, resn: GLYCAN_RESN }).length > 0,
        });

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

        applyStyles(viewer, loaded, antigen, handlers.current.options, handlers.current.selected, haloRef);
        // The viewer outlives any one structure, and 3Dmol clamps every zoomTo
        // against the current limit, so the previous structure's limit has to
        // go before this one is measured or a larger structure would be cut off.
        viewer.setZoomLimits(0, Infinity);
        viewer.zoomTo({});
        fitRef.current = viewer.getPerceivedDistance();
        limitZoomOut(viewer, container.current!, fitRef.current);

        viewer.zoomTo(own);
        viewer.zoom(1.3, 0);

        // keep zoomTo's centre and distance, replace only the orientation
        const upright = uprightQuaternion(loaded.ca, antigen.region);
        if (upright) {
          const view = viewer.getView();
          viewer.setView([view[0], view[1], view[2], view[3], ...upright]);
        }
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
    if (status === "ready" && viewer && loaded) applyStyles(viewer, loaded, antigen, options, selected, haloRef);
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

  const renderLayer = ([key, label]: readonly [keyof ViewerOptions, string]) => {
    const usable = available[key];
    return (
      <label key={key} className={usable ? "check" : "check unavailable"} title={usable ? undefined : UNAVAILABLE_REASON[key]}>
        <input type="checkbox" checked={options[key] && usable} disabled={!usable} onChange={() => toggle(key)} />
        {label}
      </label>
    );
  };

  const toggle = (key: keyof ViewerOptions) => props.onOptions({ ...options, [key]: !options[key] });

  return (
    <div className="viewer">
      <div ref={container} className="viewer-canvas" />
      <div className="viewer-toolbar">
        {PRIMARY_LAYERS.map(renderLayer)}
        <button className="ghost-button" aria-expanded={showMore} onClick={() => setShowMore((open) => !open)}>
          {showMore ? "Fewer layers" : "More layers"}
        </button>
        {showMore &&
          MORE_LAYERS.filter(([key]) => key !== "showTruth" || !isUpload(antigen)).map(renderLayer)}
        <button
          className="ghost-button"
          onClick={() => viewerRef.current?.zoomTo({ chain: antigen.chain, hetflag: false }, 500)}
        >
          Reset view
        </button>
      </div>
      {status !== "ready" && (
        <div className="viewer-status" role="status">
          {status === "loading" ? "Loading structure…" : `Unable to load structure: ${message}`}
        </div>
      )}
      <div className="viewer-legend">
        <span className="legend-title">Epitope probability</span>
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
