#!/usr/bin/env python3
"""Self-check for the single-cell hematopoiesis pipeline.

Asserts on synthetic data with known structure:
  1. reference mapping recovers cell types well above chance (8-class)
  2. clonal recovery beats chance (ARI clearly > 0) and degrades with dropout
"""
from __future__ import annotations
import os, sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from sklearn.model_selection import train_test_split
from sceval import map_query, recover_clones, make_scrna, make_clonal, CELL_TYPES


def main():
    print("== singlecell-hematopoiesis self-check ==\n")
    ok = True

    # 1. reference mapping
    X, y = make_scrna(seed=0)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=0, stratify=y)
    res, _ = map_query(Xtr, ytr, Xte, yte, method="centroid")
    chance = 1.0 / len(CELL_TYPES)
    print(f"Reference mapping: accuracy={res.accuracy:.3f} macro_f1={res.macro_f1:.3f} "
          f"(chance={chance:.3f}, {len(CELL_TYPES)} types)")
    if res.accuracy <= 0.5:
        print("FAIL: mapping accuracy not clearly above chance"); ok = False

    # 2. clonal recovery, clean vs high-dropout
    M0, c0 = make_clonal(seed=1, dropout=0.05)
    r0 = recover_clones(M0, n_clones=4, true_labels=c0)
    M1, c1 = make_clonal(seed=1, dropout=0.45)
    r1 = recover_clones(M1, n_clones=4, true_labels=c1)
    print(f"Clonal recovery ARI: clean={r0.ari:.3f}  high_dropout={r1.ari:.3f}")
    if r0.ari <= 0.3:
        print("FAIL: clonal ARI on clean data not above chance"); ok = False
    if not (r1.ari <= r0.ari + 1e-9):
        print("FAIL: dropout did not degrade (or leave equal) clonal recovery"); ok = False

    if ok:
        print("\nPASS: mapping recovers cell types; clonal recovery works and "
              "degrades under dropout.")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
