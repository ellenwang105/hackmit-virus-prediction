import type { AntigenSummary } from "./types";
import type { AlignRequest, AlignResponse } from "./align.worker";

/** Below this raw BLOSUM62 score a local match is indistinguishable from chance in a set this size. */
export const MIN_SCORE = 50;
const MIN_LENGTH = 20;
const MAX_LENGTH = 3000;
const AMINO_ACIDS = new Set("ACDEFGHIKLMNPQRSTVWY");

export interface ParsedQuery {
  sequence: string;
  /** FASTA header without the ">", when there was one */
  name: string | null;
  /** how many FASTA records were pasted; only the first is searched */
  records: number;
  warnings: string[];
  /** set when the input cannot be searched; `sequence` is then not to be used */
  error: string | null;
}

/** Turn whatever was pasted (raw residues or FASTA) into one clean protein sequence. */
export function parseQuery(text: string): ParsedQuery {
  const lines = text.split(/\r?\n/);
  const headers = lines.filter((l) => l.trimStart().startsWith(">"));
  const warnings: string[] = [];

  let body: string[];
  let name: string | null = null;
  if (headers.length) {
    const first = lines.findIndex((l) => l.trimStart().startsWith(">"));
    name = lines[first].trim().slice(1).trim() || null;
    const rest = lines.slice(first + 1);
    const next = rest.findIndex((l) => l.trimStart().startsWith(">"));
    body = next < 0 ? rest : rest.slice(0, next);
    if (headers.length > 1) warnings.push(`${headers.length} FASTA records pasted; searching the first only.`);
  } else {
    body = lines;
  }

  // numbers, spaces, gap dashes and a trailing stop are formatting, not residues
  const sequence = body.join("").toUpperCase().replace(/[^A-Z]/g, "");
  const base = { sequence, name, records: headers.length, warnings, error: null };

  if (!sequence) return base;
  if (sequence.length < MIN_LENGTH) {
    return { ...base, error: `Too short: ${sequence.length} residues. Paste at least ${MIN_LENGTH} amino acids.` };
  }
  if (sequence.length > MAX_LENGTH) {
    return { ...base, error: `Too long: ${sequence.length} residues. HA is about 570; the limit is ${MAX_LENGTH}.` };
  }
  if (sequence.length >= 30 && [...sequence].filter((c) => "ACGTUN".includes(c)).length / sequence.length > 0.95) {
    return { ...base, error: "This looks like DNA or RNA. Paste the translated protein (amino-acid) sequence instead." };
  }
  const ambiguous = [...sequence].filter((c) => !AMINO_ACIDS.has(c)).length;
  if (ambiguous > 0) warnings.push(`${ambiguous} ambiguous or unusual residues will be treated as unknown (X).`);
  return base;
}

export interface SequenceHit {
  /** the chain to open: one representative of everything in `chains` */
  id: string;
  pdb: string;
  chain: string;
  /** every chain of this entry carrying the identical sequence, `id` included */
  chains: string[];
  subtype: string;
  year: number;
  split: AntigenSummary["split"];
  chainType: string;
  nEpitope: number;
  /** raw BLOSUM62 alignment score */
  score: number;
  /** identical residues in the alignment */
  matches: number;
  /**
   * Share of the pasted sequence found identically in this structure, 0-1: identical residues
   * over the query's length. 1 means every residue of the query matches. This is the headline
   * number and what hits are sorted on; unlike identity below it cannot be inflated by a short
   * piece that matches perfectly.
   */
  match: number;
  /** identical residues / alignment columns, 0-1; how alike the aligned stretch is */
  identity: number;
  /** residues of the pasted sequence inside the alignment, 0-1 */
  queryCoverage: number;
}

export interface SearchResult {
  query: ParsedQuery;
  hits: SequenceHit[];
  bestScore: number;
}

// Sorting on identity can promote a chain with a modest raw score (a short piece that matches
// perfectly), so the pool that gets a full alignment has to be wide enough to include those.
const SHORTLIST = 150;
const MAX_HITS = 30;

function runWorker(request: AlignRequest): Promise<AlignResponse> {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL("./align.worker.ts", import.meta.url), { type: "module" });
    worker.onmessage = (event: MessageEvent<AlignResponse>) => {
      worker.terminate();
      resolve(event.data);
    };
    worker.onerror = (event) => {
      worker.terminate();
      reject(new Error(event.message || "alignment worker failed"));
    };
    worker.postMessage(request);
  });
}

/**
 * Rank every antigen chain by local alignment to the query.
 *
 * Many chains are identical (the copies of one trimer, or the same construct
 * in different entries), so each distinct sequence is aligned once and fanned
 * back out to its chains. Chains of one entry that share a sequence collapse
 * into a single hit; the same sequence in two different entries stays as two,
 * because they are different structures with different bound antibodies.
 */
export async function searchSequences(
  query: ParsedQuery,
  sequences: Record<string, string>,
  index: AntigenSummary[],
): Promise<SearchResult> {
  const owners = new Map<string, string[]>();
  for (const [id, sequence] of Object.entries(sequences)) {
    owners.set(sequence, [...(owners.get(sequence) ?? []), id]);
  }
  const subjects = [...owners.keys()];
  const chainsOf = subjects.map((s) => owners.get(s)!);
  const meta = new Map(index.map((a) => [a.id, a]));

  const { hits: aligned } = await runWorker({ query: query.sequence, subjects, top: SHORTLIST });
  const bestScore = aligned[0]?.score ?? 0;

  const grouped = new Map<string, SequenceHit>();
  for (const hit of aligned) {
    if (hit.score < MIN_SCORE) continue;
    for (const id of chainsOf[hit.index]) {
      const antigen = meta.get(id);
      if (!antigen) continue;
      const key = `${antigen.pdb}|${hit.index}`;
      const existing = grouped.get(key);
      if (existing) {
        existing.chains.push(id);
        continue;
      }
      grouped.set(key, {
        id,
        pdb: antigen.pdb,
        chain: antigen.chain,
        chains: [id], // sorted below, when the representative is chosen
        subtype: antigen.subtype,
        year: antigen.year,
        split: antigen.split,
        chainType: antigen.chainType,
        nEpitope: antigen.nEpitope,
        score: hit.score,
        matches: hit.matches,
        match: hit.matches / query.sequence.length,
        identity: hit.columns ? hit.matches / hit.columns : 0,
        queryCoverage: (hit.queryEnd - hit.queryStart) / query.sequence.length,
      });
    }
  }

  // Highest match first, so the list reads top to bottom in the order of the numbers shown.
  // Sorting on identity over the aligned stretch instead would put short fragments that match
  // perfectly above the full-length structures that are actually closest. Exact ties are common
  // (one construct solved many times); those go to the tighter alignment, then to the entries
  // with the most observed antibody contacts, since they have the most to look at.
  const hits = [...grouped.values()]
    .map((h) => ({ ...h, chains: [...h.chains].sort() }))
    .sort(
      (a, b) =>
        b.match - a.match ||
        b.identity - a.identity ||
        b.score - a.score ||
        b.nEpitope - a.nEpitope ||
        a.id.localeCompare(b.id),
    )
    .map((h) => ({ ...h, id: h.chains[0], chain: meta.get(h.chains[0])!.chain }))
    .slice(0, MAX_HITS);

  return { query, hits, bestScore };
}

/** A database sequence with a few substitutions, to show what a near-match looks like. */
export function exampleQuery(sequences: Record<string, string>, preferred: string): { text: string; source: string; changes: number } | null {
  const ids = Object.keys(sequences);
  if (!ids.length) return null;
  const source = sequences[preferred] ? preferred : ids.reduce((a, b) => (sequences[b].length > sequences[a].length ? b : a));
  const swap: Record<string, string> = { K: "R", R: "K", D: "E", E: "D", S: "T", T: "S", N: "D", I: "V", V: "I", L: "M", A: "S", G: "S", Q: "E", F: "Y", Y: "F" };
  const letters = [...sequences[source]];
  let changes = 0;
  for (let i = 25; i < letters.length; i += 40) {
    const replacement = swap[letters[i]];
    if (replacement) {
      letters[i] = replacement;
      changes++;
    }
  }
  return { text: `>example ${source} with ${changes} substitutions\n${letters.join("")}`, source, changes };
}
