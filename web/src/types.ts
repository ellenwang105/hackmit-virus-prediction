/**
 * The contract with scripts/08_export_web_data.py. When results change, re-run
 * that script; as long as these shapes hold, the UI needs no changes.
 */

export type Split = "train" | "val" | "test_group2" | "test_B" | "excluded";

/** One row of data/index.json, enough to draw the antigen picker. */
export interface AntigenSummary {
  id: string;
  pdb: string;
  chain: string;
  subtype: string;
  year: number;
  phylo: string;
  split: Split;
  chainType: string;
  nResidues: number;
  nEpitope: number;
  /** AUPRC on this single chain; null when it has no observed epitope */
  auprc: number | null;
}

/** data/antigens/<id>.json: column-oriented, index i is one residue in every array. */
export interface Antigen {
  id: string;
  pdb: string;
  chain: string;
  subtype: string;
  year: number;
  phylo: string;
  split: Split;
  chainType: string;
  /** PDB (author) residue number; identifies a residue only together with `ic` */
  num: number[];
  ic: string[];
  aa: string[];
  /** residue on the standard H3 numbering */
  ha: number[];
  region: ("head" | "stem")[];
  /** antigenic site A-E, or "" */
  site: string[];
  rbs: number[];
  fusion: number[];
  sasa: (number | null)[];
  glycan: (number | null)[];
  epitope: number[];
  nAb: number[];
  nContact: number[];
  score: number[];
  constraint: (number | null)[];
  durability: (number | null)[];
}

export interface MetricRow {
  split: string;
  scorer: string;
  auprc: number;
  "precision@20": number;
  mcc: number;
  chains: number;
  residues: number;
  positive_rate: number;
}

/** data/phylogeny.json: subtypes arranged by pairwise sequence identity. */
export interface PhyloNode {
  /** subtype name at a tip; empty at an internal node */
  name: string;
  /** 1 - mean identity to the other side of the split */
  height: number;
  children: PhyloNode[];
}

export interface SubtypeMeta {
  representative: string;
  chains: number;
  structures: number;
  residues: number;
  heldOut: boolean;
  split: string;
  group: "group1" | "group2" | "B";
}

export interface Phylogeny {
  tree: PhyloNode;
  subtypes: Record<string, SubtypeMeta>;
  identity: { names: string[]; matrix: number[][] };
}
