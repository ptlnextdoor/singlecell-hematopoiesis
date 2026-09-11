#!/usr/bin/env python3
"""Stub: load a real scRNA-seq .h5ad reference/query for label transfer.

Kept as a stub so the repo runs end-to-end on synthetic data (selfcheck.py)
without a large download. Fill against a real hematopoiesis atlas.

Suggested: scanpy (pip install scanpy anndata).
  import scanpy as sc
  ref = sc.read_h5ad("reference_atlas.h5ad")   # ref.X cells x genes, ref.obs['cell_type']
  q   = sc.read_h5ad("query.h5ad")
Then pass ref.X, ref.obs['cell_type'].values, q.X to sceval.map_query.
"""
import sys
print("TODO: implement with scanpy; see docstring.", file=sys.stderr)
raise SystemExit(2)
