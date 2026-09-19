import type { Antigen } from "./types";

export interface Patch {
  rank: number;
  indices: number[];
  meanScore: number;
  peakScore: number;
  /** fraction of the patch's residues that are observed epitope contacts */
  observed: number;
  /** compressed HA numbering, e.g. "144-146, 155, 158" */
  span: string;
  tags: string[];
}

interface Options {
  /** consider the top fraction of residues by score */
  topFraction: number;
  /** ...but never anything scoring below this */
  minScore: number;
  /** CA-CA distance (angstrom) that joins two residues into one patch */
  radius: number;
  minSize: number;
  maxPatches: number;
}

const DEFAULTS: Options = { topFraction: 0.15, minScore: 0.05, radius: 9, minSize: 3, maxPatches: 8 };

/** Group the highest-scoring residues into spatially connected surface patches. */
export function findPatches(
  antigen: Antigen,
  ca: Float32Array | null,
  overrides: Partial<Options> = {},
): Patch[] {
  if (!ca) return [];
  const o = { ...DEFAULTS, ...overrides };
  const n = antigen.score.length;
  const ranked = [...Array(n).keys()].sort((a, b) => antigen.score[b] - antigen.score[a]);
  const cutoff = Math.max(o.minScore, antigen.score[ranked[Math.floor(n * o.topFraction)]] ?? 0);
  const candidates = ranked.filter((i) => antigen.score[i] >= cutoff && !Number.isNaN(ca[3 * i]));

  // union-find over candidates that sit close together in 3D
  const parent = new Map<number, number>(candidates.map((i) => [i, i]));
  const find = (i: number): number => {
    let root = i;
    while (parent.get(root) !== root) root = parent.get(root)!;
    parent.set(i, root);
    return root;
  };
  const r2 = o.radius * o.radius;
  for (let x = 0; x < candidates.length; x++) {
    for (let y = x + 1; y < candidates.length; y++) {
      const a = candidates[x];
      const b = candidates[y];
      const dx = ca[3 * a] - ca[3 * b];
      const dy = ca[3 * a + 1] - ca[3 * b + 1];
      const dz = ca[3 * a + 2] - ca[3 * b + 2];
      if (dx * dx + dy * dy + dz * dz <= r2) parent.set(find(a), find(b));
    }
  }

  const groups = new Map<number, number[]>();
  for (const i of candidates) {
    const root = find(i);
    groups.set(root, [...(groups.get(root) ?? []), i]);
  }

  return [...groups.values()]
    .filter((g) => g.length >= o.minSize)
    .map((indices) => ({ indices, mass: indices.reduce((s, i) => s + antigen.score[i], 0) }))
    .sort((a, b) => b.mass - a.mass)
    .slice(0, o.maxPatches)
    .map(({ indices }, k) => {
      const scores = indices.map((i) => antigen.score[i]);
      return {
        rank: k + 1,
        indices: indices.sort((a, b) => a - b),
        meanScore: scores.reduce((s, v) => s + v, 0) / scores.length,
        peakScore: Math.max(...scores),
        observed: indices.filter((i) => antigen.epitope[i]).length / indices.length,
        span: compressRanges(indices.map((i) => antigen.ha[i])),
        tags: tagsFor(antigen, indices),
      };
    });
}

function compressRanges(numbers: number[]): string {
  const sorted = [...new Set(numbers)].sort((a, b) => a - b);
  const runs: string[] = [];
  for (let i = 0; i < sorted.length; ) {
    let j = i;
    while (j + 1 < sorted.length && sorted[j + 1] === sorted[j] + 1) j++;
    runs.push(j > i ? `${sorted[i]}–${sorted[j]}` : `${sorted[i]}`);
    i = j + 1;
  }
  return runs.length > 5 ? `${runs.slice(0, 5).join(", ")}, …` : runs.join(", ");
}

/** Annotations covering at least a third of the patch, most specific first. */
function tagsFor(antigen: Antigen, indices: number[]): string[] {
  const share = (test: (i: number) => boolean) => indices.filter(test).length / indices.length;
  const tags: string[] = [];
  if (share((i) => antigen.rbs[i] === 1) >= 1 / 3) tags.push("Receptor pocket");
  if (share((i) => antigen.fusion[i] === 1) >= 1 / 3) tags.push("Fusion machinery");
  for (const site of ["A", "B", "C", "D", "E"]) {
    if (share((i) => antigen.site[i] === site) >= 1 / 3) tags.push(`Site ${site}`);
  }
  tags.push(share((i) => antigen.region[i] === "head") >= 0.5 ? "Head" : "Stem");
  return tags;
}
