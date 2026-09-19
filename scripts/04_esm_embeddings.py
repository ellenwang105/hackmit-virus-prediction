"""ESM-2 per-residue embeddings for every antigen chain, reduced with PCA.

The 896 antigen chains collapse to 361 unique sequences, so only the unique ones
are run through the model and the result is mapped back by sequence.

PCA is fitted on training-split residues alone. Fitting it on everything would
let the held-out phylogenetic groups shape the feature basis they are tested on.
"""

import pathlib

import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import PCA

import esm

RESIDUES_PATH = "data/residues_v1.parquet"
FASTA_PATH = "antigen_sequences.fasta"
OUTPUT_PATH = pathlib.Path("data/esm_v1.parquet")

MODEL_NAME = "esm2_t33_650M_UR50D"
REPRESENTATION_LAYER = 33
COMPONENTS = 64
BATCH_TOKENS = 8000

KEY = ["PDB", "antigen_chain", "residue_number", "insertion_code"]


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            sequences[name] = []
        elif name:
            sequences[name].append(line)
    return {k: "".join(v) for k, v in sequences.items()}


def embed(sequences, model, alphabet, device):
    """Per-residue representations for each unique sequence, batched by token count."""
    converter = alphabet.get_batch_converter()
    order = sorted(sequences, key=lambda s: -len(s))
    embeddings = {}

    batch, tokens = [], 0
    batches = []
    for sequence in order:
        if batch and tokens + len(sequence) > BATCH_TOKENS:
            batches.append(batch)
            batch, tokens = [], 0
        batch.append(sequence)
        tokens += len(sequence)
    if batch:
        batches.append(batch)

    for done, group in enumerate(batches, start=1):
        _, _, batch_tokens = converter([(str(i), s) for i, s in enumerate(group)])
        with torch.no_grad():
            out = model(batch_tokens.to(device), repr_layers=[REPRESENTATION_LAYER])
        representations = out["representations"][REPRESENTATION_LAYER].cpu().numpy()
        for index, sequence in enumerate(group):
            # strip the leading BOS token; keep one vector per residue
            embeddings[sequence] = representations[index, 1:len(sequence) + 1].astype(np.float32)
        print(f"  batch {done}/{len(batches)}  ({len(group)} sequences)")
    return embeddings


residues = pd.read_parquet(RESIDUES_PATH)
fasta = read_fasta(FASTA_PATH)
unique = sorted(set(fasta.values()))
print(f"{len(fasta)} antigen chains, {len(unique)} unique sequences, "
      f"longest {max(len(s) for s in unique)} residues")

device = "mps" if torch.backends.mps.is_available() else "cpu"
model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
model = model.to(device).eval()
print(f"{MODEL_NAME} on {device}")

embeddings = embed(unique, model, alphabet, device)

# map each chain's residues onto its embedding rows via seq_index, which
# label_epitopes.py assigned in the same order the FASTA was written
frames = []
missing_chains = []
for antigen_id, group in residues.groupby("antigen_id"):
    sequence = fasta.get(antigen_id)
    if sequence is None or antigen_id not in fasta:
        missing_chains.append(antigen_id)
        continue
    matrix = embeddings[sequence]
    group = group[group["seq_index"] < len(matrix)]
    frame = pd.DataFrame(matrix[group["seq_index"].to_numpy()])
    frame.columns = [f"raw_{i}" for i in range(matrix.shape[1])]
    for column in KEY:
        frame[column] = group[column].to_numpy()
    frame["split_group"] = group["split_group"].to_numpy()
    frames.append(frame)

stacked = pd.concat(frames, ignore_index=True)
raw_columns = [c for c in stacked.columns if c.startswith("raw_")]
print(f"\n{len(stacked):,} residues embedded ({len(missing_chains)} chains without a sequence)")

train_mask = stacked["split_group"] == "train"
pca = PCA(n_components=COMPONENTS, random_state=0)
pca.fit(stacked.loc[train_mask, raw_columns].to_numpy(dtype=np.float32))
reduced = pca.transform(stacked[raw_columns].to_numpy(dtype=np.float32))
print(f"PCA fitted on {int(train_mask.sum()):,} training residues; "
      f"{COMPONENTS} components explain {pca.explained_variance_ratio_.sum():.1%} of variance")

output = stacked[KEY].copy()
for index in range(COMPONENTS):
    output[f"esm_{index}"] = reduced[:, index].astype(np.float32)

OUTPUT_PATH.parent.mkdir(exist_ok=True)
output.to_parquet(OUTPUT_PATH, index=False)
print(f"wrote {OUTPUT_PATH}")
