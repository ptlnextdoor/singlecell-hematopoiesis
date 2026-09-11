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
    import json
    fs.use_paper_style()

    # Panel A data: REAL result on the Asiri Lab's GSE262440 (3 donors), from
    # results/gse262440_mapping.json. Falls back to synthetic if it is absent.
    res_path = os.path.join(HERE, "..", "results", "gse262440_mapping.json")
    real = json.load(open(res_path)) if os.path.exists(res_path) else None

    # Panel B data: ARI vs dropout (synthetic ground truth, 20 seeds)
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

    if real:
        donors = real["per_donor"]
        labels = [f"donor {i+1}\n({d['n_cells_gated']} cells, {d['n_populations']} pops)"
                  for i, d in enumerate(donors)]
        acc = [100 * d["heldout_accuracy_rna_only"] for d in donors]
        chance = [100 * d["chance_accuracy"] for d in donors]
        x = np.arange(len(donors)); w = 0.38
        axA.bar(x - w/2, chance, w, color=fs.PALETTE[7], label="random guess")
        axA.bar(x + w/2, acc, w, color=fs.PALETTE[0], label="RNA-only mapper")
        for xi, a, c in zip(x, acc, chance):
            axA.annotate(f"{a/c:.1f}x", (xi + w/2, a), ha="center", va="bottom",
                         fontsize=8, xytext=(0, 3), textcoords="offset points")
        axA.set_xticks(x); axA.set_xticklabels(labels, fontsize=7)
        axA.set_ylabel("held-out accuracy (%)")
        axA.set_title("Real data: RNA recovers protein-defined cell types at ~2x chance")
        axA.set_ylim(0, max(acc) * 1.3)
        fs.opaque_legend(axA, loc="upper left")
    else:
        axA.text(0.5, 0.5, "run scripts/gse262440_mapping.py first", ha="center",
                 va="center", transform=axA.transAxes)

    axB.errorbar(dropouts, means, yerr=stds, color=fs.PALETTE[0], lw=2,
                 marker="o", capsize=3, label="mean ± sd, 20 runs")
    axB.axhline(0, color=fs.PALETTE[7], lw=1, ls="--", label="random grouping")
    axB.set_xlabel("fraction of true mutations missed (allele dropout)")
    axB.set_ylabel("clone recovery score (1 = perfect, 0 = random)")
    axB.set_title("Clone recovery degrades gracefully as calls go missing")
    axB.set_ylim(-0.08, 1.08)
    fs.opaque_legend(axB)

    fig.suptitle("Single-cell hematopoiesis: mapping cells to types on the Asiri Lab's "
                 "GSE262440, and recovering clones from mutations")
    fs.assert_no_clip(fig)
    out = os.path.join(HERE, "..", "figures", "singlecell_eval")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fs.save(fig, out)
    print("saved", out + ".png/.pdf")


if __name__ == "__main__":
    main()
