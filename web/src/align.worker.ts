import { align, encode, scoreOnly, type Alignment } from "./align";

export interface AlignRequest {
  query: string;
  /** distinct subject sequences; the caller maps indices back to antigen chains */
  subjects: string[];
  /** how many of the best-scoring subjects to align in full */
  top: number;
}

export interface AlignedSubject extends Alignment {
  index: number;
}

export interface AlignResponse {
  hits: AlignedSubject[];
}

const ctx = self as unknown as Worker;

ctx.onmessage = (event: MessageEvent<AlignRequest>) => {
  const { query, subjects, top } = event.data;
  const q = encode(query);
  const encoded = subjects.map(encode);

  // rank everything with the cheap pass, then trace back only the leaders
  const ranked = encoded
    .map((subject, index) => ({ index, score: scoreOnly(q, subject) }))
    .sort((a, b) => b.score - a.score || a.index - b.index)
    .slice(0, top);

  const response: AlignResponse = {
    hits: ranked.map(({ index }) => ({ index, ...align(q, encoded[index]) })),
  };
  ctx.postMessage(response);
};
