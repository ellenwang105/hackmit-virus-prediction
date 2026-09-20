import type { MetricRow } from "../types";
import { SPLIT_LABEL } from "../data";

interface Props {
  metrics: MetricRow[];
}

const SPLIT_ORDER = ["val", "test_group2", "test_B"];
const SPLIT_NOTE: Record<string, string> = {
  val: "H5 · same clade as training",
  test_group2: "H3, H4, H7, H10, H14 · clade unseen in training",
  test_B: "23–27% identity to the training set",
};

/**
 * The rationale, the method and the limits, for a reader deciding whether to trust the tool.
 *
 * The tool above assumes familiarity with epitopes and with why a conserved
 * region is a better target than a variable one. This states that once, with the
 * mechanism drawn rather than described, and gives the limits before a reader
 * has to go looking for them.
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
            <text x="72" y="70" className="fig-small">stalk</text>
            <text x="72" y="104" className="fig-small">head</text>

            <path d="M230 98 L212 90 L212 106 Z" fill="var(--accent)" />
            <rect x="230" y="82" width="28" height="32" rx="6" fill="var(--accent)" />
            <path d="M258 90 L288 72" stroke="var(--accent)" strokeWidth="7" strokeLinecap="round" />
            <path d="M258 106 L288 124" stroke="var(--accent)" strokeWidth="7" strokeLinecap="round" />
            <text x="252" y="66" className="fig-small accent" textAnchor="middle">antibody</text>
            <line x1="192" y1="98" x2="210" y2="98" stroke="var(--accent)" strokeWidth="1.4" markerEnd="url(#mech-tip)" />

            <text x="412" y="14" className="fig-label" textAnchor="end">HOST CELL</text>
            <rect x="332" y="22" width="80" height="166" rx="8" fill="currentColor" opacity=".09" />
            <circle cx="346" cy="98" r="5.5" fill="currentColor" opacity=".5" />
            <circle cx="346" cy="140" r="5.5" fill="currentColor" opacity=".5" />
            <text x="358" y="174" className="fig-small">receptor</text>

            <line x1="134" y1="140" x2="338" y2="140" stroke="currentColor" strokeWidth="1.2" strokeDasharray="4 4" markerEnd="url(#mech-tip)" />
            <text x="236" y="156" className="fig-small" textAnchor="middle">head engages receptor</text>

            <line x1="148" y1="98" x2="174" y2="98" stroke="currentColor" strokeWidth="1.2" strokeDasharray="3 3" />
            <line x1="174" y1="89" x2="190" y2="107" stroke="var(--accent)" strokeWidth="2.2" />
            <line x1="174" y1="107" x2="190" y2="89" stroke="var(--accent)" strokeWidth="2.2" />
            <text x="182" y="126" className="fig-small accent" textAnchor="middle">blocked</text>
          </svg>
          <figcaption>
            <b>Neutralization.</b> An epitope is the surface patch an antibody binds. An antibody bound
            to the receptor-binding head sterically blocks attachment to host receptors and prevents entry.
          </figcaption>
        </figure>

        <figure className="figure">
          <svg viewBox="0 0 420 210" role="img" aria-label="The spike head mutates between seasons so old antibodies stop fitting, while the stem stays unchanged">
            <text x="8" y="14" className="fig-label">EARLIER STRAIN</text>
            <line x1="60" y1="142" x2="60" y2="92" stroke="currentColor" strokeWidth="7" strokeLinecap="round" opacity=".4" />
            <circle cx="60" cy="76" r="19" fill="currentColor" opacity=".26" />
            <rect x="8" y="142" width="104" height="13" rx="6" fill="currentColor" opacity=".14" />
            <path d="M42 52 L60 36 L78 52" fill="none" stroke="var(--accent)" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" />
            <text x="60" y="26" className="fig-small accent" textAnchor="middle">binds</text>

            <line x1="146" y1="56" x2="192" y2="56" stroke="currentColor" strokeWidth="1.2" markerEnd="url(#mech-tip)" />
            <text x="169" y="48" className="fig-small" textAnchor="middle">substitutions</text>
            <line x1="146" y1="124" x2="192" y2="124" stroke="currentColor" strokeWidth="1.2" markerEnd="url(#mech-tip)" />
            <text x="169" y="116" className="fig-small" textAnchor="middle">conserved</text>

            <text x="228" y="14" className="fig-label">LATER STRAIN</text>
            <line x1="282" y1="142" x2="282" y2="92" stroke="currentColor" strokeWidth="7" strokeLinecap="round" opacity=".4" />
            <path d="M262 70 q8 -12 18 -5 q4 -14 16 -6 q10 8 2 18 q6 12 -8 15 q-12 8 -20 -3 q-14 -1 -8 -19 z" fill="currentColor" opacity=".26" />
            <rect x="228" y="142" width="104" height="13" rx="6" fill="currentColor" opacity=".14" />
            <path d="M264 52 L282 36 L300 52" fill="none" stroke="currentColor" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" strokeDasharray="4 3" opacity=".45" />
            <line x1="272" y1="30" x2="292" y2="48" stroke="var(--accent)" strokeWidth="2.2" />
            <line x1="272" y1="48" x2="292" y2="30" stroke="var(--accent)" strokeWidth="2.2" />

            <text x="210" y="180" className="fig-small" textAnchor="middle">head domain: variable · stalk: conserved</text>
            <text x="210" y="196" className="fig-small muted-text" textAnchor="middle">stalk substitutions impair membrane fusion</text>
          </svg>
          <figcaption>
            <b>Antigenic drift and escape.</b> The head domain tolerates substitutions, so antibodies raised
            against earlier strains lose binding. The stalk mediates membrane fusion and is functionally
            constrained, which is what the <b>Durability</b> view ranks. Several SARS-CoV-2 monoclonal
            antibodies lost activity when variants altered their epitopes.
          </figcaption>
        </figure>
      </div>

      <div className="explain-cols">
        <div>
          <h3>Method</h3>
          <ol className="explain-steps">
            <li>
              <b>Featurization.</b> Each residue is described from the antigen alone, with the antibody
              removed: solvent exposure on the chain and on the assembled trimer, protrusion, distance
              to the nearest glycan, neighborhood averages, and 128 principal components of ESM-2
              embeddings.
            </li>
            <li>
              <b>Labels.</b> A residue is an epitope if any atom lies within 4.5 Å of a bound antibody,
              pooled across every antibody complex solved for that antigen.
            </li>
            <li>
              <b>Evaluation.</b> Training uses phylogenetic group 1 only. Performance is reported on
              held-out subtypes, so near-duplicate structures cannot inflate the results.
            </li>
          </ol>
        </div>

        <div>
          <h3>Performance</h3>
          <table className="explain-table">
            <thead>
              <tr>
                <th scope="col">Evaluation set</th>
                <th scope="col">AUPRC</th>
                <th scope="col">Exposure baseline</th>
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
                <th scope="row">Random (base rate)</th>
                <td>~0.05</td>
                <td>—</td>
              </tr>
            </tbody>
          </table>
          <p className="explain-foot">
            Epitope residues are ~5% of all residues, so the random baseline for AUPRC is <b>~0.05, not 0</b>.
            On the held-out clade the model achieves about <b>7×</b> the random baseline and <b>3×</b> the
            surface-exposure baseline.
          </p>
        </div>
      </div>

      <div>
        <h3>Limitations</h3>
        <ul className="caveats">
          <li>
            <b>Predictions are candidates, not validated targets.</b> Accessibility to an antibody does not
            establish that the immune response targets a site. The stalk is conserved yet poorly
            immunogenic, which is why a universal influenza vaccine remains unsolved.
          </li>
          <li>
            <b>Labels reflect what has been crystallized.</b> Structural databases over-represent
            well-studied complexes, such as stem-directed antibodies since 2009, so head-versus-stem
            comparisons are confounded. Site-level comparisons within the same chains are not.
          </li>
          <li>
            <b>Influenza B performance is near baseline.</b> At 23–27% identity to the training set,
            language-model features do not transfer. This limit is measured, not assumed.
          </li>
          <li>
            <b>Constraint is relative.</b> It is Shannon entropy across circulating strains, normalized
            within a subtype. Influenza B is left unscored rather than borrowing values from H3.
          </li>
          <li>
            <b>Research prototype.</b> Intended for prioritizing candidates for experimental follow-up,
            not for clinical decision-making.
          </li>
        </ul>
      </div>
    </section>
  );
}
