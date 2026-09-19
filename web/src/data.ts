import type { Antigen, AntigenSummary, MetricRow, Phylogeny } from "./types";

const BASE = import.meta.env.BASE_URL;

async function getJson<T>(path: string): Promise<T> {
  const response = await fetch(`${BASE}${path}`);
  if (!response.ok) throw new Error(`${path}: ${response.status} ${response.statusText}`);
  return response.json() as Promise<T>;
}

export const loadIndex = () => getJson<AntigenSummary[]>("data/index.json");
export const loadMetrics = () => getJson<MetricRow[]>("data/metrics.json");
export const loadPhylogeny = () => getJson<Phylogeny>("data/phylogeny.json");

const antigenCache = new Map<string, Promise<Antigen>>();
export function loadAntigen(id: string): Promise<Antigen> {
  let cached = antigenCache.get(id);
  if (!cached) {
    cached = getJson<Antigen>(`data/antigens/${id}.json`);
    cached.catch(() => antigenCache.delete(id));
    antigenCache.set(id, cached);
  }
  return cached;
}

const structureCache = new Map<string, Promise<string>>();
/** mmCIF text for a PDB entry. Copes with hosts that gunzip on the fly and ones that don't. */
export function loadStructure(pdb: string): Promise<string> {
  let cached = structureCache.get(pdb);
  if (!cached) {
    cached = (async () => {
      const response = await fetch(`${BASE}structures/${pdb}.cif.gz`);
      if (!response.ok) throw new Error(`structure ${pdb}: ${response.status}`);
      const bytes = new Uint8Array(await response.arrayBuffer());
      const isGzip = bytes[0] === 0x1f && bytes[1] === 0x8b;
      if (!isGzip) return new TextDecoder().decode(bytes);
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      return new Response(stream).text();
    })();
    cached.catch(() => structureCache.delete(pdb));
    structureCache.set(pdb, cached);
  }
  return cached;
}

/** "N158": residue letter plus standard-H3 number, the label people recognise. */
export const residueLabel = (a: Antigen, i: number) => `${a.aa[i]}${a.ha[i]}`;

/** A residue can belong to several groups at once (e.g. site B and the receptor pocket). */
export const SITE_LETTERS = ["A", "B", "C", "D", "E"] as const;

export const SPLIT_LABEL: Record<string, string> = {
  train: "Training",
  val: "Validation",
  test_group2: "Test · group 2",
  test_B: "Test · influenza B",
  excluded: "Excluded",
};

export const isHeldOut = (split: string) => split === "test_group2" || split === "test_B";
