import type { MetricRow } from "../types";
import { SPLIT_LABEL } from "../data";

interface Props {
  metrics: MetricRow[];
}

const SPLIT_ORDER = ["val", "test_group2", "test_B"];
const SPLIT_NOTE: Record<string, string> = {
  val: "H5 — same branch as training",
  test_group2: "H3, H7, H10, H14 — a branch never trained on",
  test_B: "influenza B — 23–27% identity, the far edge",
};

/**
 * Why any of this matters, for someone who has not spent a week on it.
 *
 * The tool above assumes you know what an epitope is and why the stem is worth
 * more than the head. This says it once, with the mechanism drawn rather than
 * described, and states the limits before a reader has to find them.
 */
export function Mechanism({ metrics }: Props) {
  const rows = SPLIT_ORDER.map((split) => ({
    split,
    model: metrics.find((m) => m.split === split && m.scorer === "epitope_score"),
    baseline: metrics.find((m) => m.split === split && m.scorer === "rel_sasa_assembly"),
  })).filter((r) => r.model);

  return (
    <section className="explain">
      <div className="explain-figures">
        <figure className="figure">
          <svg viewBox="0 0 420 210" role="img" aria-label="An influenza spike reaches a receptor on the host cell; an antibody clamped to the spike head blocks that contact">
            <defs>
              <marker id="mech-tip" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
                <polygon points="0,1 9,5 0,9" fill="currentColor" />
              </marker>
            </defs>

            <text x="8" y="14" className="fig-label">VIRUS</text>
            <rect x="8" y="22" width="148" height="14" rx="7" fill="currentColor" opacity=".14" />
            <line x1="54" y1="36" x2="54" y2="88" stroke="currentColor" strokeWidth="7" strokeLinecap="round" opacity=".4" />
            <circle cx="54" cy="100" r="16" fill="currentColor" opacity=".26" />
            <line x1="116" y1="36" x2="116" y2="88" stroke="currentColor" strokeWidth="7" strokeLinecap="round" opacity=".4" />
            <circle cx="116" cy="100" r="16" fill="currentColor" opacity=".26" />
            <text x="72" y="70" className="fig-small">stem</text>
            <text x="72" y="104" className="fig-small">head</text>

            <path d="M230 98 L212 90 L212 106 Z" fill="var(--accent)" />
            <rect x="230" y="82" width="28" height="32" rx="6" fill="var(--accent)" />
            <path d="M258 90 L288 72" stroke="var(--accent)" strokeWidth="7" strokeLinecap="round" />
            <path d="M258 106 L288 124" stroke="var(--accent)" strokeWidth="7" strokeLinecap="round" />
            <text x="248" y="146" className="fig-small accent" textAnchor="middle">antibody</text>
            <line x1="192" y1="98" x2="210" y2="98" stroke="var(--accent)" strokeWidth="1.4" markerEnd="url(#mech-tip)" />

            <text x="412" y="14" className="fig-label" textAnchor="end">HOST CELL</text>
            <rect x="332" y="22" width="80" height="166" rx="8" fill="currentColor" opacity=".09" />
            <circle cx="346" cy="98" r="5.5" fill="currentColor" opacity=".5" />
            <circle cx="346" cy="140" r="5.5" fill="currentColor" opacity=".5" />
            <text x="358" y="174" className="fig-small">receptor</text>

            <line x1="134" y1="140" x2="338" y2="140" stroke="currentColor" strokeWidth="1.2" strokeDasharray="4 4" markerEnd="url(#mech-tip)" />
            <text x="236" y="156" className="fig-small" textAnchor="middle">head must reach the receptor</text>

            <line x1="148" y1="98" x2="174" y2="98" stroke="currentColor" strokeWidth="1.2" strokeDasharray="3 3" />
            <line x1="174" y1="89" x2="190" y2="107" stroke="var(--accent)" strokeWidth="2.2" />
            <line x1="174" y1="107" x2="190" y2="89" stroke="var(--accent)" strokeWidth="2.2" />
            <text x="182" y="192" className="fig-small accent" textAnchor="middle">blocked</text>
          </svg>
          <figcaption>
            <b>What an epitope is.</b> The patch an antibody clamps onto. Cover the right
            patch and the spike can never reach the cell, so the virus cannot get in.
          </figcaption>
        </figure>

        <figure className="figure">
          <svg viewBox="0 0 420 210" role="img" aria-label="The spike head mutates between seasons so old antibodies stop fitting, while the stem stays unchanged">
            <text x="8" y="14" className="fig-label">2015</text>
            <line x1="60" y1="142" x2="60" y2="92" stroke="currentColor" strokeWidth="7" strokeLinecap="round" opacity=".4" />
            <circle cx="60" cy="76" r="19" fill="currentColor" opacity=".26" />
            <rect x="8" y="142" width="104" height="13" rx="6" fill="currentColor" opacity=".14" />
            <path d="M42 52 L60 36 L78 52" fill="none" stroke="var(--accent)" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" />
            <text x="60" y="26" className="fig-small accent" textAnchor="middle">fits</text>

            <line x1="146" y1="56" x2="192" y2="56" stroke="currentColor" strokeWidth="1.2" markerEnd="url(#mech-tip)" />
            <text x="169" y="48" className="fig-small" textAnchor="middle">drifts</text>
            <line x1="146" y1="124" x2="192" y2="124" stroke="currentColor" strokeWidth="1.2" markerEnd="url(#mech-tip)" />
            <text x="169" y="116" className="fig-small" textAnchor="middle">cannot drift</text>

            <text x="228" y="14" className="fig-label">2025</text>
            <line x1="282" y1="142" x2="282" y2="92" stroke="currentColor" strokeWidth="7" strokeLinecap="round" opacity=".4" />
            <path d="M262 70 q8 -12 18 -5 q4 -14 16 -6 q10 8 2 18 q6 12 -8 15 q-12 8 -20 -3 q-14 -1 -8 -19 z" fill="currentColor" opacity=".26" />
            <rect x="228" y="142" width="104" height="13" rx="6" fill="currentColor" opacity=".14" />
            <path d="M264 52 L282 36 L300 52" fill="none" stroke="currentColor" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" strokeDasharray="4 3" opacity=".45" />
            <line x1="272" y1="30" x2="292" y2="48" stroke="var(--accent)" strokeWidth="2.2" />
            <line x1="272" y1="48" x2="292" y2="30" stroke="var(--accent)" strokeWidth="2.2" />

            <text x="210" y="180" className="fig-small" textAnchor="middle">the head changes shape every season; the stem cannot</text>
            <text x="210" y="196" className="fig-small muted-text" textAnchor="middle">a stem that drifts is a spike that can no longer fuse</text>
          </svg>
          <figcaption>
            <b>Why last year's shot stops working.</b> The head tolerates change, so it drifts
            and old antibodies stop fitting. The stem drives membrane fusion and cannot drift
            without breaking the spike — which is what the <b>Durability</b> tab ranks.
          </figcaption>
        </figure>
      </div>

      <div className="explain-cols">
        <div>
          <h3>How the score is made</h3>
          <ol className="explain-steps">
            <li>
              <b>Describe each residue</b> with the antibody deleted — exposure on the chain and
              on the assembled spike, protrusion, distance to the nearest glycan, neighbourhood
              averages, and 128 components of an ESM-2 embedding.
            </li>
            <li>
              <b>Learn from solved structures.</b> A residue is an epitope when any atom sits
              within 4.5 Å of a bound antibody, pooled over every antibody solved against that
              antigen.
            </li>
            <li>
              <b>Hold out whole branches.</b> Training saw group 1 only, so the numbers below
              come from subtypes on a different branch of the tree.
            </li>
          </ol>
        </div>

        <div>
          <h3>How well it does</h3>
          <table className="explain-table">
            <thead>
              <tr>
                <th scope="col">Evaluated on</th>
                <th scope="col">AUPRC</th>
                <th scope="col">vs exposure</th>
              </tr>
            </thead>
            <tbody>
              {rows.map(({ split, model, baseline }) => (
                <tr key={split} className={split === "test_group2" ? "headline" : undefined}>
                  <th scope="row">
                    {SPLIT_LABEL[split] ?? split}
                    <span className="muted"> {SPLIT_NOTE[split]}</span>
                  </th>
                  <td>{model!.auprc.toFixed(3)}</td>
                  <td className="muted">{baseline ? baseline.auprc.toFixed(3) : "—"}</td>
                </tr>
              ))}
              <tr className="muted">
                <th scope="row">Random guessing</th>
                <td>~0.05</td>
                <td>—</td>
              </tr>
            </tbody>
          </table>
          <p className="explain-foot">
            Epitopes are ~5% of residues, so <b>0.05 is the floor, not zero</b>. On the held-out
            branch the model runs about <b>7×</b> that floor and <b>3×</b> surface exposure alone.
          </p>
        </div>
      </div>

      <div>
        <h3>What this does not tell you</h3>
        <ul className="caveats">
          <li>
            <b>These are likely antibody targets, not a vaccine.</b> A patch antibodies can reach
            is not automatically one the immune system responds to — the stem is conserved
            <em> and</em> poorly immunogenic, which is why a universal flu vaccine is still unsolved.
          </li>
          <li>
            <b>Labels follow research attention, not immunity.</b> Stem complexes are
            over-represented in the PDB because they have been the field's focus since 2009, so
            head-versus-stem comparisons are confounded. Site-level comparisons are not — they
            compare regions within the same chains.
          </li>
          <li>
            <b>Influenza B sits near the baseline.</b> At 23–27% identity the language-model
            features stop transferring. That limit is measured, not assumed.
          </li>
          <li>
            <b>Constraint is relative, not absolute.</b> It is entropy over circulating strains
            normalised within a subtype, and influenza B is deliberately unscored rather than
            borrowed from H3.
          </li>
        </ul>
      </div>
    </section>
  );
}
