"""Clonal recovery from single-cell mutation matrices.

Input: a binary cells x variants matrix (1 = variant present in that cell).
Task: cluster cells into clones and score against ground-truth clone labels.
This mirrors single-cell mutational profiling used to track AML clones from
diagnosis to relapse (MRD).
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score


@dataclass
class ClonalResult:
    ari: float
    n_clones_found: int


def recover_clones(mut_matrix: np.ndarray, n_clones: int,
                   true_labels=None) -> ClonalResult:
    """Cluster cells by shared mutations (Hamming distance, average linkage)."""
    X = np.asarray(mut_matrix, dtype=float)
    clust = AgglomerativeClustering(n_clusters=n_clones, metric="hamming",
                                    linkage="average")
    pred = clust.fit_predict(X)
    ari = adjusted_rand_score(true_labels, pred) if true_labels is not None else float("nan")
    return ClonalResult(float(ari), int(len(np.unique(pred))))
