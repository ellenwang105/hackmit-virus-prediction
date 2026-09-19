"""Evaluation for a rare-positive, per-residue task.

Epitopes are ~5% of residues, so accuracy and ROC-AUC both look impressive while
meaning nothing. Metrics are computed per antigen chain and then averaged:
pooling lets the largest chains dominate the score.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, matthews_corrcoef

TOP_K = 20


def per_antigen_auprc(frame, score_column="score", label_column="is_epitope"):
    """Mean average-precision across antigen chains that have both classes."""
    values = []
    for _, group in frame.groupby("antigen_id"):
        labels = group[label_column].to_numpy()
        if labels.min() == labels.max():
            continue
        values.append(average_precision_score(labels, group[score_column].to_numpy()))
    return float(np.mean(values)), len(values)


def precision_at_k(frame, score_column="score", label_column="is_epitope", k=TOP_K):
    """Of the k residues we would actually show a vaccine designer, how many are real."""
    values = []
    for _, group in frame.groupby("antigen_id"):
        if group[label_column].sum() == 0:
            continue
        top = group.nlargest(min(k, len(group)), score_column)
        values.append(top[label_column].mean())
    return float(np.mean(values)), len(values)


def best_threshold(frame, score_column="score", label_column="is_epitope"):
    """Threshold maximising MCC. Tune on validation only, never on test."""
    labels = frame[label_column].to_numpy()
    scores = frame[score_column].to_numpy()
    candidates = np.quantile(scores, np.linspace(0.5, 0.995, 60))
    best, best_mcc = candidates[0], -1.0
    for threshold in candidates:
        mcc = matthews_corrcoef(labels, (scores >= threshold).astype(int))
        if mcc > best_mcc:
            best, best_mcc = threshold, mcc
    return float(best), float(best_mcc)


def evaluate(frame, threshold, score_column="score", label_column="is_epitope"):
    auprc, n_chains = per_antigen_auprc(frame, score_column, label_column)
    p_at_k, _ = precision_at_k(frame, score_column, label_column)
    predicted = (frame[score_column].to_numpy() >= threshold).astype(int)
    return {
        "auprc": round(auprc, 4),
        f"precision@{TOP_K}": round(p_at_k, 4),
        "mcc": round(float(matthews_corrcoef(frame[label_column], predicted)), 4),
        "chains": n_chains,
        "residues": len(frame),
        "positive_rate": round(float(frame[label_column].mean()), 4),
    }


def comparison_table(frame, score_columns, threshold_frame=None,
                     label_column="is_epitope"):
    """Evaluate several score columns side by side, e.g. model against baselines."""
    reference = threshold_frame if threshold_frame is not None else frame
    rows = {}
    for column in score_columns:
        threshold, _ = best_threshold(reference, column, label_column)
        rows[column] = evaluate(frame, threshold, column, label_column)
    return pd.DataFrame(rows).T
