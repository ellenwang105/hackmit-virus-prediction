"""Shared constants for the epitope pipeline."""

# Leading ESM-2 principal components used as features. scripts/04 fits the PCA
# wider than this and the components are nested, so training and export just
# take the first N — but they must take the SAME N, or the exported scores come
# from a feature matrix the model was never fitted on.
ESM_COMPONENTS = 128
