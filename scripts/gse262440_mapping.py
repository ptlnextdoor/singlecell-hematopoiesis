#!/usr/bin/env python3
"""Real-data run on the Asiri Lab's published reference (GEO GSE262440).

Loads one healthy bone-marrow donor from the Cell Reports 2025 dataset
(BD Rhapsody WTA + ADT surface-protein panel). Two steps:

1. Gate cells into the paper's HSPC populations using the SAME surface markers
   the paper uses (CD34, CD38, CD90, CD45RA, CD69, CLL1, CD2), from the ADT
   counts. This gives protein-defined labels, independent of gene expression.
2. Ask: can `sceval.map_query` recover those protein-defined labels from gene
   expression alone? Held-out accuracy tells us whether the RNA carries the
   population identity the paper defined by protein.

No data is shipped. Download from GEO:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262440/suppl/
Run:
  python scripts/gse262440_mapping.py --csv GSE262440_Healthy1-WTA-ADT_RSEC_MolsPerCell.csv.gz
"""
from __future__ import annotations
import os, sys, gzip, json, argparse, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
import numpy as np

MARKERS = ["CD34", "CD38", "CD90", "CD45RA", "CD69", "CLL1", "CD2"]


def load(csv_gz):
    import csv
    with gzip.open(csv_gz, "rt") as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")]
    header, data = rows[0], rows[1:]
    adt_idx = {h.split("|")[0]: i for i, h in enumerate(header) if h.endswith("|pAbO")}
    gene_idx = [i for i, h in enumerate(header) if i > 0 and not h.endswith("|pAbO")]
    M = np.array([[float(x) for x in r[1:]] for r in data], dtype=np.float32)
    adt = {m: M[:, adt_idx[m] - 1] for m in MARKERS if m in adt_idx}
    genes = M[:, [i - 1 for i in gene_idx]]
    return adt, genes, [header[i] for i in gene_idx]


def gate(adt):
    """Protein gating following the paper's marker logic (thresholds are
    per-marker medians of the positive tail; simple and transparent)."""
    def pos(m, q=0.5):
        v = adt[m]
        thr = np.quantile(v[v > 0], q) if (v > 0).any() else np.inf
        return v > thr
    cd34 = pos("CD34"); cd38lo = ~pos("CD38"); cd90 = pos("CD90"); ra = pos("CD45RA")
    cd69 = pos("CD69"); cll1 = pos("CLL1"); cd2 = pos("CD2")
    hspc = cd34 & cd38lo
    lab = np.full(len(cd34), "other", dtype=object)
    lab[hspc & cd90 & ~ra] = "HSC (CD90+CD45RA-)"
    mpp = hspc & ~cd90 & ~ra
    lab[mpp & cd69] = "MPP CD69+"
    lab[mpp & ~cd69 & cll1] = "MPP CLL1+ (myeloid-biased)"
    lab[mpp & ~cd69 & ~cll1] = "MPP CLL1-CD69- (erythroid-biased)"
    lmpp = hspc & ~cd90 & ra
    lab[lmpp & cd2] = "LMPP CD2+"
    lab[lmpp & ~cd2 & ~cll1] = "LMPP CD2-"
    lab[lmpp & ~cd2 & cll1] = "GMP CLL1+"
    lab[cd34 & ~cd38lo] = "CD34+CD38+ progenitor"
    return lab


def main():
    from sklearn.model_selection import train_test_split
    from sceval import map_query
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", default="results/gse262440_mapping.json")
    ap.add_argument("--min-cells", type=int, default=15)
    args = ap.parse_args()

    adt, genes, gene_names = load(args.csv)
    print(f"loaded {genes.shape[0]} cells x {genes.shape[1]} genes, ADT markers: {list(adt)}")
    labels = gate(adt)
    uniq, counts = np.unique(labels, return_counts=True)
    for u, c in sorted(zip(uniq, counts), key=lambda x: -x[1]):
        print(f"  {c:5d}  {u}")

    # keep populations with enough cells for a held-out split
    keep_labels = {u for u, c in zip(uniq, counts) if c >= args.min_cells and u != "other"}
    m = np.array([l in keep_labels for l in labels])
    X, y = np.log1p(genes[m]), labels[m]
    # top variable genes to keep it fast and stable
    var = X.var(0); top = np.argsort(var)[-2000:]
    X = X[:, top]
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=0, stratify=y)
    res, pred = map_query(Xtr, ytr, Xte, yte, method="centroid", n_pcs=30)
    chance = 1.0 / len(keep_labels)
    out = {
        "dataset": "GSE262440 Healthy1 (Asiri Lab, Cell Reports 2025)",
        "n_cells_total": int(genes.shape[0]),
        "n_cells_gated": int(m.sum()),
        "populations": sorted(keep_labels),
        "n_populations": len(keep_labels),
        "chance_accuracy": round(chance, 3),
        "heldout_accuracy_rna_only": round(res.accuracy, 3),
        "heldout_macro_f1": round(res.macro_f1, 3),
        "note": ("Labels come from the paper's surface-protein markers (ADT). "
                 "Mapping uses gene expression only. Accuracy > chance means RNA "
                 "carries the protein-defined population identity."),
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
