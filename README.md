# singlecell-hematopoiesis

Two questions blood-cancer researchers ask about single cells, answered with
small, testable code:

1. **What kind of cell is this?** Given a cell we've never seen, match it to
   the closest known type in a labeled reference (like matching a face to a
   photo album).
2. **Which cells came from the same clone?** Cancer grows as families of cells
   that share the same mutations. Given a table of which cell has which
   mutation, group the cells back into their families.

I built this after reading the Asiri Lab's (Stanford) single-cell map of human
blood-cell development and their work tracking leukemia clones through
treatment. It's my own independent project, not affiliated with the lab.

## What it does

![result](figures/singlecell_eval.png)

**Left:** on test data with 8 known cell types, the matcher gets every cell
right (each row lights up only on its own diagonal).

**Right:** clone recovery works well when mutation calls are clean, then
degrades as more calls go missing. In real single-cell DNA sequencing, a
mutation is often missed even when it's there (called "dropout"), so I test
how much of that the method can tolerate. Score of 1 = perfect grouping,
0 = random.

## Try it in 30 seconds (no data needed)

```bash
pip install -r requirements.txt
python scripts/selfcheck.py
```

This makes fake cells and fake clones where I *know* the right answer, then
checks the code finds it.

## Run it on real data

I ran the matcher on the Asiri Lab's own published data (GEO GSE262440, one
healthy bone-marrow donor, 3,630 cells). Their dataset measures both gene
expression and surface proteins on every cell, so I could:

1. Label each cell using the **same surface-protein markers the paper uses**
   (CD34, CD38, CD90, CD45RA, CD69, CLL1, CD2). This found all 8 of the paper's
   stem/progenitor populations, from 21 to 896 cells each.
2. Ask whether **gene expression alone** can recover those protein-defined
   labels on held-out cells.

| | held-out accuracy |
|---|---|
| Random guess (8 populations) | 12.5% |
| Gene expression only | 25.2% |

So RNA carries some of the population identity, but only modestly. These
progenitor subtypes are hard to tell apart by transcriptome alone, which is
exactly why the paper needed protein markers to define them in the first
place. Honest result, and a useful baseline.

```bash
python scripts/gse262440_mapping.py --csv GSE262440_Healthy1-WTA-ADT_RSEC_MolsPerCell.csv.gz
```

No data is shipped. The file is public at the GEO link in the script.

## What's in here

```
src/sceval/mapping.py   match new cells to known types
src/sceval/clonal.py    group cells into clones by shared mutations
src/sceval/synth.py     make fake data with known answers, for testing
scripts/selfcheck.py    quick test on fake data
scripts/make_figures.py makes the figure above
scripts/gse262440_mapping.py  real run on the Asiri Lab's published data
```

## License

MIT
