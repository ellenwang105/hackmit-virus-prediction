/**
 * Local protein alignment for matching a pasted sequence against the antigens
 * in the data: Smith-Waterman with affine gaps and BLOSUM62, the same scoring
 * BLASTP uses by default (a gap of length k costs 11 + k).
 *
 * Deliberately free of imports so it runs unchanged in a worker, in the page
 * and under plain Node, where it is checked against Biopython.
 */

export const GAP_OPEN = 11;
export const GAP_EXTEND = 1;

/** Row and column order of BLOSUM62 below; anything else is folded into X. */
const LETTERS = "ARNDCQEGHILKMFPSTWYV";
const X = LETTERS.length;
const SIZE = LETTERS.length + 1;

const BLOSUM62 = [
  "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0", // A
  "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3", // R
  "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3", // N
  "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3", // D
  "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1", // C
  "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2", // Q
  "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2", // E
  "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3", // G
  "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3", // H
  "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3", // I
  "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1", // L
  "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2", // K
  "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1", // M
  "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1", // F
  "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2", // P
  "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2", // S
  "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0", // T
  "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3", // W
  "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1", // Y
  "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4", // V
];

/** X against each letter in LETTERS order (X against X is -1); symmetric with the X column. */
const X_SCORES = [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1];

/** Flat SIZE x SIZE table matching the published BLOSUM62, including the X row and column. */
export const MATRIX = new Int8Array(SIZE * SIZE).fill(-1);
BLOSUM62.forEach((row, i) => {
  row.split(" ").forEach((value, j) => {
    MATRIX[i * SIZE + j] = Number(value);
  });
  MATRIX[i * SIZE + X] = X_SCORES[i];
  MATRIX[X * SIZE + i] = X_SCORES[i];
});

export function encode(sequence: string): Uint8Array {
  const out = new Uint8Array(sequence.length);
  for (let i = 0; i < sequence.length; i++) {
    const index = LETTERS.indexOf(sequence[i]);
    out[i] = index < 0 ? X : index;
  }
  return out;
}

const NEG = -1_000_000_000;

/** Best local alignment score, in O(n) memory. Used to rank every candidate cheaply. */
export function scoreOnly(query: Uint8Array, subject: Uint8Array): number {
  const m = query.length;
  const n = subject.length;
  const H = new Int32Array(n + 1);
  const E = new Int32Array(n + 1).fill(NEG);
  let best = 0;

  for (let i = 1; i <= m; i++) {
    const row = query[i - 1] * SIZE;
    let diagonal = 0; // H[i-1][j-1]
    let left = 0; // H[i][j-1]
    let F = NEG;
    for (let j = 1; j <= n; j++) {
      const up = H[j]; // H[i-1][j]
      const e = Math.max(E[j] - GAP_EXTEND, up - GAP_OPEN - GAP_EXTEND);
      E[j] = e;
      F = Math.max(F - GAP_EXTEND, left - GAP_OPEN - GAP_EXTEND);
      let h = diagonal + MATRIX[row + subject[j - 1]];
      if (e > h) h = e;
      if (F > h) h = F;
      if (h < 0) h = 0;
      diagonal = up;
      H[j] = h;
      left = h;
      if (h > best) best = h;
    }
  }
  return best;
}

export interface Alignment {
  score: number;
  /** identical residues in the alignment */
  matches: number;
  /** alignment columns, counting gap columns */
  columns: number;
  gaps: number;
  /** half-open 0-based spans of the aligned region in each sequence */
  queryStart: number;
  queryEnd: number;
  subjectStart: number;
  subjectEnd: number;
}

/** Best local alignment with traceback, for the few candidates that survive ranking. */
export function align(query: Uint8Array, subject: Uint8Array): Alignment {
  const m = query.length;
  const n = subject.length;
  const w = n + 1;
  const H = new Int32Array((m + 1) * w);
  const E = new Int32Array((m + 1) * w).fill(NEG);
  const F = new Int32Array((m + 1) * w).fill(NEG);
  let best = 0;
  let bestI = 0;
  let bestJ = 0;

  for (let i = 1; i <= m; i++) {
    const row = query[i - 1] * SIZE;
    for (let j = 1; j <= n; j++) {
      const at = i * w + j;
      E[at] = Math.max(E[at - w] - GAP_EXTEND, H[at - w] - GAP_OPEN - GAP_EXTEND);
      F[at] = Math.max(F[at - 1] - GAP_EXTEND, H[at - 1] - GAP_OPEN - GAP_EXTEND);
      let h = H[at - w - 1] + MATRIX[row + subject[j - 1]];
      if (E[at] > h) h = E[at];
      if (F[at] > h) h = F[at];
      if (h < 0) h = 0;
      H[at] = h;
      if (h > best) {
        best = h;
        bestI = i;
        bestJ = j;
      }
    }
  }

  let i = bestI;
  let j = bestJ;
  let state: "H" | "E" | "F" = "H";
  let matches = 0;
  let columns = 0;
  let gaps = 0;

  while (i > 0 && j > 0) {
    const at = i * w + j;
    if (state === "H") {
      const h = H[at];
      if (h === 0) break;
      if (h === H[at - w - 1] + MATRIX[query[i - 1] * SIZE + subject[j - 1]]) {
        columns++;
        if (query[i - 1] === subject[j - 1] && query[i - 1] !== X) matches++;
        i--;
        j--;
      } else if (h === E[at]) {
        state = "E";
      } else {
        state = "F";
      }
    } else if (state === "E") {
      // a gap in the subject: this column consumes query residue i
      columns++;
      gaps++;
      if (E[at] !== E[at - w] - GAP_EXTEND) state = "H";
      i--;
    } else {
      // a gap in the query: this column consumes subject residue j
      columns++;
      gaps++;
      if (F[at] !== F[at - 1] - GAP_EXTEND) state = "H";
      j--;
    }
  }

  return {
    score: best,
    matches,
    columns,
    gaps,
    queryStart: i,
    queryEnd: bestI,
    subjectStart: j,
    subjectEnd: bestJ,
  };
}
