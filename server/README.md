# Scoring server

Lets the web app score a structure that is not in the dataset: upload a `.cif` or
`.pdb` (gzipped is fine), or give a PDB id and the server fetches it from RCSB.

```
.venv/bin/uvicorn server.app:app --port 8000      # repo root; ~8 s to load ESM-2
cd web && npm run dev                              # Vite proxies /api to :8000
```

| Endpoint | Does |
|---|---|
| `POST /api/predict` | multipart `file` → scores every hemagglutinin chain |
| `POST /api/predict/pdb/{id}` | fetches the entry from RCSB, then scores it |
| `GET /api/structures/{token}` | the coordinates back, for the 3D viewer |

The response holds one object per HA chain in the same column-oriented shape as
`web/src/types.ts` `Antigen`, plus `applicability` (how far the antigen is from the
training data) and `warnings`. Errors are a 4xx with a `detail` written for the user.

Notes:

- The model needs `models/xgb_v1.json` and `models/esm_pca_v1.npz`. The second is
  written by `scripts/04_esm_embeddings.py`; without it uploads cannot be scored.
- One request at a time, on purpose. torch and xgboost only coexist in one process
  on macOS if both run single-threaded; see the note at the top of `epitope/predict.py`.
- Nothing is written to disk. The last 16 structures are held in memory so the
  viewer can fetch coordinates back.
- Set `VITE_API_URL` when building the web app if the API lives on another origin.
- Tests: `.venv/bin/python -m pytest tests -q`
