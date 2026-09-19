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
| `constraint_score` | **C** — empty, to be filled in |
| `durability_score` | **D = E × C** — empty, to be filled in |

## What's next: filling in `C` and `D`

`C` asks how much the virus can afford to change a position. Read it off
evolution rather than learning it:

1. Pull H3 (and H1) HA sequences across seasons from the NCBI Influenza Virus
   Resource.
2. Align them and compute Shannon entropy per alignment column.
3. Map alignment columns onto `ha_number` so it joins to this file.
4. `C = 1 - normalised_entropy`, so 1 means never varies.
5. `D = epitope_score * C`.

The check that matters: **antigenic site B should lose to the receptor pocket.**
Site B currently has the highest epitope score on H3 (0.230) but drifts every
season. The receptor pocket scores 0.201 and cannot drift. If `D` flips their
order, the thesis holds.

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
