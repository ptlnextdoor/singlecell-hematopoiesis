#!/usr/bin/env python3
"""Generate figures for singlecell-hematopoiesis from committed code.

Panel A: reference-mapping confusion (query cell types recovered), a heatmap.
Panel B: clonal recovery ARI vs allele-dropout rate, showing graceful
degradation, the failure mode single-cell mutational profiling must handle.

Run: python scripts/make_figures.py
"""
from __future__ import annotations
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)  # vendored figure_style.py lives next to this script

from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sceval import map_query, recover_clones, make_scrna, make_clonal, CELL_TYPES
import figure_style as fs


def main():
    fs.use_paper_style()

    # Panel A data: reference mapping confusion
    X, y = make_scrna(seed=0)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=0, stratify=y)
    res, pred = map_query(Xtr, ytr, Xte, yte, method="centroid")
    cm = confusion_matrix(yte, pred, normalize="true")

    # Panel B data: ARI vs dropout
    dropouts = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    aris = {d: [] for d in dropouts}
    for seed in range(20):
        for d in dropouts:
            M, c = make_clonal(seed=seed, dropout=d)
            aris[d].append(recover_clones(M, n_clones=4, true_labels=c).ari)
    means = [np.mean(aris[d]) for d in dropouts]
    stds = [np.std(aris[d]) for d in dropouts]

    fig, axes = fs.new_figure(cols=2, panels=2)
    axA, axB = axes

    im = axA.imshow(cm, cmap=fs.SEQUENTIAL, vmin=0, vmax=1, aspect="auto")
    axA.set_xticks(range(len(CELL_TYPES))); axA.set_yticks(range(len(CELL_TYPES)))
    axA.set_xticklabels(CELL_TYPES, rotation=45, ha="right", fontsize=6)
    axA.set_yticklabels(CELL_TYPES, fontsize=6)
    axA.set_xlabel("predicted cell type")
    axA.set_ylabel("true cell type")
    axA.set_title(f"Reference mapping: {res.accuracy*100:.0f}% accuracy, "
                  f"{len(CELL_TYPES)} types")
    fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04, label="row-normalized rate")

    axB.errorbar(dropouts, means, yerr=stds, color=fs.PALETTE[0], lw=2,
                 marker="o", capsize=3)
    axB.axhline(0, color=fs.PALETTE[7], lw=0.8, ls="--")
    axB.set_xlabel("allele-dropout rate")
    axB.set_ylabel("clonal recovery ARI (1 = perfect, 0 = chance)")
    axB.set_title("Clonal recovery degrades gracefully with dropout")

    fig.suptitle("Single-cell hematopoiesis: cell-type mapping and clonal "
                 "recovery on synthetic ground truth")
    fs.assert_no_clip(fig)
    out = os.path.join(HERE, "..", "figures", "singlecell_eval")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fs.save(fig, out)
    print("saved", out + ".png/.pdf")


if __name__ == "__main__":
    main()
