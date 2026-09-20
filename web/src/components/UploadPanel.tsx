import { useRef, useState, type DragEvent } from "react";
import type { UploadResult } from "../types";
import { scoreFile, scorePdb } from "../api";

interface Props {
  upload: UploadResult | null;
  activeChain: string | null;
  onResult: (result: UploadResult) => void;
  onPickChain: (chain: string) => void;
  onClear: () => void;
}

const PDB_ID = /^[0-9][A-Za-z0-9]{3}$/;

/**
 * Score a structure the model has not precomputed: drop a file, or give a PDB id
 * and the server fetches it. The search below only browses results that were
 * computed ahead of time; this is the part that runs the model.
 */
export function UploadPanel({ upload, activeChain, onResult, onPickChain, onClear }: Props) {
  const [busy, setBusy] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [dragging, setDragging] = useState(false);
  const [pdbId, setPdbId] = useState("");
  const fileInput = useRef<HTMLInputElement>(null);

  async function run(label: string, work: () => Promise<UploadResult>) {
    setBusy(label);
    setError(null);
    try {
      onResult(await work());
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
  }

  const scoreThis = (file: File | undefined) => {
    if (file) void run(file.name, () => scoreFile(file));
  };

  const onDrop = (event: DragEvent) => {
    event.preventDefault();
    setDragging(false);
    scoreThis(event.dataTransfer.files[0]);
  };

  const submitId = (event: React.FormEvent) => {
    event.preventDefault();
    const id = pdbId.trim();
    if (!PDB_ID.test(id)) {
      setError("A PDB ID is four characters: a digit followed by three letters or digits (e.g., 5K9K).");
      return;
    }
    void run(id.toUpperCase(), () => scorePdb(id));
  };

  return (
    <section className="upload" aria-label="Score your own structure">
      <h2>Score a structure</h2>

      <button
        type="button"
        className={dragging ? "dropzone over" : "dropzone"}
        disabled={busy !== null}
        onClick={() => fileInput.current?.click()}
        onDragOver={(e) => {
          e.preventDefault();
          setDragging(true);
        }}
        onDragLeave={() => setDragging(false)}
        onDrop={onDrop}
      >
        <strong>Drop a structure file</strong>
        <span>mmCIF or PDB · .gz accepted</span>
      </button>
      <input
        ref={fileInput}
        id="structure-file"
        type="file"
        accept=".cif,.mmcif,.pdb,.ent,.gz"
        hidden
        onChange={(e) => {
          scoreThis(e.target.files?.[0]);
          e.target.value = "";
        }}
      />

      <form className="upload-id" onSubmit={submitId}>
        <input
          id="pdb-id"
          className="input"
          value={pdbId}
          maxLength={4}
          placeholder="or enter a PDB ID"
          aria-label="PDB ID"
          disabled={busy !== null}
          onChange={(e) => setPdbId(e.target.value)}
        />
        <button className="ghost-button" type="submit" disabled={busy !== null || pdbId.trim() === ""}>
          Score
        </button>
      </form>

      {busy && (
        <p className="upload-busy" role="status">
          Scoring {busy}…
        </p>
      )}
      {error && (
        <p className="upload-error" role="alert">
          {error}
        </p>
      )}

      {upload && (
        <div className="upload-result">
          <div className="upload-result-head">
            <strong>{upload.label}</strong>
            <button className="ghost-button" onClick={onClear}>
              Clear
            </button>
          </div>
          <div className="upload-chains" role="group" aria-label="Hemagglutinin chains found">
            {upload.antigens.map((a) => (
              <button
                key={a.chain}
                className="chip"
                aria-pressed={a.chain === activeChain}
                onClick={() => onPickChain(a.chain)}
              >
                chain {a.chain}
              </button>
            ))}
          </div>
          <p className="muted upload-note">
            {upload.antigens.length} hemagglutinin chain{upload.antigens.length === 1 ? "" : "s"} scored
            {upload.chains.some((c) => !c.isAntigen) ? "; non-HA chains (e.g., antibody) excluded" : ""}.
          </p>
        </div>
      )}
    </section>
  );
}
