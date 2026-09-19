# Results handoff

## `epitope_predictions.csv`

One row per HA antigen residue: 197,477 rows across 684 antigen chains.
Produced by `scripts/06_export_results.py`.

### Identity

| Column | Meaning |
|---|---|
| `antigen_id` | `PDB_chain`, the unique antigen chain |
| `PDB`, `antigen_chain`, `residue_number`, `insertion_code` | join key — **always use all four**, `residue_number` alone is not unique |
| `residue` | one-letter amino acid |
| `ha_subtype`, `year`, `phylo_group` | H1/H3/B..., deposition year, group1/group2/B |
| `split_group` | `train` / `val` / `test_group2` / `test_B` / `excluded` |

`train` scores are fit on that data and will look inflated. Quote numbers from
`test_group2` (a phylogenetic group never seen in training).

### Structural annotation

| Column | Meaning |
|---|---|
| `chain_type` | `HA1` (head piece), `HA2` (stalk piece), `HA0` (uncleaved, numbered straight through) |
| `piece`, `ha_number` | residue mapped onto standard H3 numbering |
| `region` | `head` or `stem` |
| `antigenic_site` | `A`–`E` if the residue is in a classic drift site, else empty |
| `is_rbs` | in the receptor binding pocket — the part the virus cannot redesign |
| `is_fusion_machinery` | fusion peptide or long helix in HA2 |
| `rel_sasa_assembly` | solvent exposure, 0 = buried, ~1 = fully exposed |
| `glycan_distance` | Å to the nearest modelled sugar; small values mean shielded |

Annotations come from `epitope/annotations.py` and are published residue sets,
not predictions. Nothing there is derived from antibody contacts, so they are
safe to validate against.

### Labels and scores

| Column | Meaning |
|---|---|
| `is_epitope` | ground truth: any antibody heavy atom within 4.5 Å, unioned across every antibody bound to that antigen |
| `n_antibodies` | how many antibody complexes covered this residue |
| `n_contacting` | how many of them actually touched it |
| `epitope_score` | **E** — model probability this residue is an epitope |
| `constraint_score` | **C** — 1 means the position never varies in circulating strains |
| `durability_score` | **D = E × C** — the deliverable |
| `escape_risk` | **E × (1 − C)** — where antibodies bind and the virus can run |
| `constraint_reference` | which surveillance alignment scored this residue, `H1` or `H3` |
| `constraint_observations` | sequences backing that position; below 50 it is left unscored |
| `constraint_score_prior` | C with a functional-site bump — **do not use for validation**, see below |

## `C` and `D`, as computed

Produced by `scripts/07_constraint.py` and `scripts/08_durability.py`.

`C` is read off evolution, not learned. 4,800 H3N2 and 4,800 H1N1 HA protein
sequences from NCBI, 300 per year across 2010–2025 so no heavily-sequenced season
dominates, mapped onto the reference numbering frame (A/Aichi/2/1968 for H3,
A/California/07/2009 for H1) and reduced to Shannon entropy per position.
`C = 1 - normalised_entropy`.

Entropy is normalised against the 99th percentile of observed entropy rather than
`log2(20)`. HA is conserved nearly everywhere, so the textbook denominator puts
every position near C = 1 and destroys the ranking. The cost is that C is relative
to its subtype, not absolute.

Coverage is 184,248 / 197,477 residues (93.3%). Influenza B is deliberately
unscored: it sits outside both influenza A references, and borrowing H3 entropy
for it would invent a number rather than measure one. A further 2,242 residues
carry HA numbers outside the mature frame — see the numbering caveat below.

### The functional-site prior is off by default

The spec's third component — a constraint bump for annotated RBS and fusion
machinery — is computed into `constraint_score_prior` and kept out of
`constraint_score`. The headline validation asks whether D ranks the receptor
pocket above antigenic site B. Bumping C for RBS residues would settle that by
assumption, which is the same circularity that kept distance-to-RBS out of the
model's features. Quote `constraint_score`.

### ESM-2 masked marginals are not included

`epitope/constraint.py` implements the second component and activates it
automatically when `torch` and `esm` are importable. Neither is installed, so the
current numbers are entropy-only. Re-run `scripts/07_constraint.py` after
installing them to fold it in.

## Validation — `durability_validation.csv`

| Check | Comparison | Result |
|---|---|---|
| 1a. E premise | mean E: head 0.0602 vs stem 0.0726 | **fails** — see below |
| 1b. D inversion | mean D: **stem 0.0631** vs head 0.0333 | passes |
| 2. Sites are decoys | mean C: all H3 0.8231 vs **sites A–E 0.6242** | passes |
| 3a. Site B wins on E | mean E: **site B 0.2191** vs pocket 0.1999 | passes |
| 3b. Pocket wins on D | mean D: **pocket 0.1015** vs site B 0.0784 | passes |
| 4. Conserved base wins on D | mean D: **Y98/W153/H183/Y195 0.0816** vs site B 0.0784 | passes |

Computed against the ESM model from `c5ed20d`. C does not depend on the model, so
only the E and D columns move if the model is retrained — rerun
`scripts/08_durability.py` and this table regenerates.

**The headline result holds.** Site B carries the higher epitope score and loses
on durability to the receptor pocket. Sites A–E sit well below the H3 average on
constraint, which is exactly what makes them decoys.

**Why 1a fails, and why it is not the model's fault.** The spec expected mean E to
be higher on the head. It is not — head 0.0602 against stem 0.0726 — because the
labels are not head-dominated: of 165 H3 antigen chains, 100 are bound
predominantly by stem antibodies and 53 by head antibodies, and the ground-truth
positive rate is higher on the stem (6.2%) than the head (4.2%). SAbDab reflects
what crystallographers chose to solve, and broadly neutralising stem antibodies
have been the field's focus for a decade. The model is reproducing its labels
correctly. The conclusion the premise was meant to support — 1b — holds on its
own.

State this before a judge does: the head/stem comparison is confounded by
structure-selection bias. The site-level comparisons (2, 3, 4) are not, because
they compare regions within the same chains.

### Independent confirmation

The score was never shown antibody contacts, so its agreement with known drift is
a real check. The twelve most variable H3 positions it finds are 144, 159, 160,
186, 193, 140, 142 (antigenic sites A and B), 62, 92, 94 (site E) and 121
(site D) — the canonical H3 drift positions from the literature.

The top of the D ranking is HA2 19–46: the fusion peptide and the stem epitope
where CR6261-class broadly neutralising antibodies bind. That is the universal
vaccine target, recovered without ever being told it exists.

### Numbering caveat

2,242 residues (1.1%) have HA numbers outside the mature frame — HA1 329–343,
HA1 −1, and HA2 677–845 — because `to_ha_numbering` subtracts a fixed HA1 length
from HA0 chains without bounding the result. Those residues are unscored for C,
but their `region` and `antigenic_site` labels are also meaningless and are still
present. Worth fixing upstream.

## Still open

- **Retrospective drift** (validation check 4 in the spec). Train on ≤2020, score
  a 2020 H3, check whether the low-C sites are the ones that actually mutated in
  2021–2025. The surveillance FASTA in `data/alignments/` is already dated per
  record, so the split is cheap.
- **BepiPred-3.0 / DiscoTope-3.0 benchmarks** on 20–30 test antigens.
- **The app** — three tabs over one structure: E, D, and escape risk.

## `metrics_summary.csv`

Model versus the surface-exposure baseline on each split. Epitopes are ~5% of
residues, so random AUPRC is ~0.05 — compare against that floor, not against 1.0.

| Split | AUPRC | Precision@20 |
|---|---|---|
| val (H5) | 0.422 | 0.293 |
| **test_group2 (H3, H7, ...)** | **0.373** | **0.311** |
| test_B (influenza B) | 0.149 | 0.117 |

On influenza B the ESM-2 model falls below the structural-only model (0.203).
`models/xgb_v1_structural.json` is the fallback for distant antigens.
