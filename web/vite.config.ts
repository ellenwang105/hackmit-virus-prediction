import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  // relative asset URLs so `vite build` output works from any static host path
  base: "./",
  build: { chunkSizeWarningLimit: 1500 },
  // the scoring server; run it with `uvicorn server.app:app --port 8000`
  server: { proxy: { "/api": "http://localhost:8000" } },
});
