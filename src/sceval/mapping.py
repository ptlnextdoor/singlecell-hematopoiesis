"""Reference mapping: label transfer from a labeled reference to query cells.

Given a labeled reference (cells x genes, with cell-type labels) and unlabeled
query cells, assign each query cell a type. This is the core operation behind
projecting leukemia cells onto a normal-hematopoiesis atlas.

Uses PCA + nearest-centroid (fast, interpretable) with an optional kNN mode.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, f1_score


@dataclass
class MappingResult:
    accuracy: float
    macro_f1: float
    n_query: int


def _embed(ref_X, query_X, n_pcs=20):
    # Drop zero-variance genes (constant across the reference); scaling them
    # divides by zero and poisons PCA with inf/nan on real data.
    ref_X = np.asarray(ref_X, dtype=float)
    query_X = np.asarray(query_X, dtype=float)
    keep = ref_X.std(axis=0) > 1e-8
    if not keep.any():
        raise ValueError("reference has no informative (non-constant) features")
    ref_X, query_X = ref_X[:, keep], query_X[:, keep]
    scaler = StandardScaler().fit(ref_X)
    ref_s = scaler.transform(ref_X)
    q_s = scaler.transform(query_X)
    n_pcs = max(1, min(n_pcs, ref_s.shape[1], ref_s.shape[0] - 1))
    pca = PCA(n_components=n_pcs, random_state=0).fit(ref_s)
    return pca.transform(ref_s), pca.transform(q_s)


def map_query(ref_X, ref_labels, query_X, query_labels=None,
              method="centroid", n_pcs=20) -> MappingResult:
    ref_e, q_e = _embed(ref_X, query_X, n_pcs)
    if method == "centroid":
        classes = np.unique(ref_labels)
        cents = np.stack([ref_e[ref_labels == c].mean(0) for c in classes])
        d = ((q_e[:, None, :] - cents[None]) ** 2).sum(-1)
        pred = classes[d.argmin(1)]
    else:  # knn
        clf = KNeighborsClassifier(n_neighbors=15).fit(ref_e, ref_labels)
        pred = clf.predict(q_e)
    if query_labels is None:
        return MappingResult(float("nan"), float("nan"), len(pred)), pred
    acc = accuracy_score(query_labels, pred)
    f1 = f1_score(query_labels, pred, average="macro")
    return MappingResult(float(acc), float(f1), len(pred)), pred
