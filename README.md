# singlecell-hematopoiesis

A small, reproducible pipeline for **single-cell analysis of hematopoiesis**:
cell-type reference mapping and clonal-structure recovery, the two computational
problems at the center of systems hematology.

Built as an independent study inspired by the Asiri Lab (systems hematology,
Stanford; systemshematology.org) line of work, e.g. the single-cell reference
atlas of human hematopoiesis and single-cell mutational profiling for AML MRD.
Not affiliated with the lab.

## What this does

Two classic single-cell tasks, on data small enough to run anywhere:

1. **Reference mapping** - project query cells onto a labeled reference of
   hematopoietic cell states with a nearest-centroid / kNN classifier, and report
   label-transfer accuracy. Mirrors mapping leukemia cells onto a normal
   hematopoiesis atlas.
2. **Clonal recovery** - from a cells x mutations binary matrix, cluster cells
   into clones and measure how well recovered clones match ground-truth clones
   (adjusted Rand index). Mirrors single-cell mutational profiling of AML.

## Quick start

```bash
python -m pip install -r requirements.txt
python scripts/selfcheck.py     # synthetic scRNA + mutation data, no download
```

The self-check builds synthetic single-cell data with known cell types and known
clones, runs both tasks, and asserts recovery well above chance.

## Using real data

- Reference mapping: any labeled scRNA-seq reference (e.g. a hematopoiesis atlas)
  as an `.h5ad`; adapter stub in `scripts/prep_scrna.py`.
- Clonal recovery: a single-cell DNA variant matrix (cells x variants).

## Repo layout

```
src/sceval/     reference mapping, clonal clustering, metrics
scripts/        selfcheck, prep stub
```

## Before you email about this

Read one Asiri Lab paper end to end (start with the Cell Reports single-cell
hematopoiesis framework, 2025) and add a note reproducing a specific idea from
it. That specificity is what makes the outreach land.

## License

MIT
