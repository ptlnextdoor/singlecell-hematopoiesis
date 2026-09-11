"""Synthetic single-cell data with known cell types and known clones."""
from __future__ import annotations
import numpy as np

CELL_TYPES = ["HSC", "MPP", "CMP", "GMP", "MEP", "Erythroid", "Myeloid", "Lymphoid"]


def make_scrna(n_per_type=120, n_genes=200, seed=0, noise=1.0):
    """Return (X (N,genes), labels (N,)) with type-specific marker programs."""
    rng = np.random.default_rng(seed)
    Xs, ys = [], []
    # each type has a distinct set of upregulated marker genes
    markers = {t: rng.choice(n_genes, size=15, replace=False) for t in CELL_TYPES}
    for ti, t in enumerate(CELL_TYPES):
        base = rng.standard_normal((n_per_type, n_genes)) * noise
        base[:, markers[t]] += 3.0
        Xs.append(base)
        ys.append(np.full(n_per_type, ti))
    X = np.vstack(Xs)
    y = np.concatenate(ys)
    perm = rng.permutation(len(y))
    return X[perm], y[perm]


def make_clonal(n_cells=300, n_variants=12, n_clones=4, seed=0, dropout=0.1):
    """Return (mut_matrix (cells,variants) binary, clone_labels (cells,)).

    Each clone owns a nested set of variants (clonal evolution). Dropout flips
    some present variants to absent, simulating single-cell allele dropout.
    """
    rng = np.random.default_rng(seed)
    # nested variant sets: clone k has variants 0..(k+1)*step
    step = max(1, n_variants // n_clones)
    clone_variants = [set(range(0, min(n_variants, (k + 1) * step))) for k in range(n_clones)]
    labels = rng.integers(0, n_clones, size=n_cells)
    M = np.zeros((n_cells, n_variants), dtype=int)
    for i in range(n_cells):
        for v in clone_variants[labels[i]]:
            if rng.random() > dropout:
                M[i, v] = 1
    return M, labels
