import type { UploadResult } from "./types";

/** Empty in development, where Vite proxies /api to the scoring server. */
const API = (import.meta.env.VITE_API_URL as string | undefined) ?? "";
export const apiUrl = (path: string) => `${API}${path}`;

/**
 * Scoring a structure needs the Python server. In development Vite proxies /api
 * to it; a production build has one only if VITE_API_URL says where, so a static
 * deployment simply has no upload rather than a box that always errors.
 */
export const UPLOAD_ENABLED = import.meta.env.DEV || Boolean(API);

const UNREACHABLE =
  "The scoring service is not reachable. Start it from the repository root with " +
  ".venv/bin/uvicorn server.app:app --port 8000 and retry.";

async function send(path: string, init: RequestInit): Promise<UploadResult> {
  let response: Response;
  try {
    response = await fetch(apiUrl(path), init);
  } catch {
    throw new Error(UNREACHABLE);
  }
  if (!response.ok) {
    let detail: string | null = null;
    try {
      const body = await response.json();
      if (typeof body?.detail === "string") detail = body.detail;
    } catch {
      /* not JSON: the dev proxy answers 5xx with plain text when nothing is listening */
    }
    if (detail) throw new Error(detail);
    throw new Error(response.status >= 500 ? UNREACHABLE : `Scoring failed (${response.status}).`);
  }
  return response.json() as Promise<UploadResult>;
}

export function scoreFile(file: File) {
  const body = new FormData();
  body.append("file", file);
  return send("/api/predict", { method: "POST", body });
}

export const scorePdb = (id: string) =>
  send(`/api/predict/pdb/${encodeURIComponent(id.trim())}`, { method: "POST" });
