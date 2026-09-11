# Upstream: build on THEIR code

Extend the Asiri Lab's own published code and reference, not a blind rebuild.

Lab GitHub org: **systemsheme** (Asiri Lab) / personal: **ediriwas**

- **ediriwas/adult-human-hspc** (Cell Reports 2025 MPP paper support code)
  https://github.com/ediriwas/adult-human-hspc
  - R markdown scripts (Seurat/ArchR), reference metadata (RDS), supplemental tables.
  - Processed data: GEO **GSE262440**; raw under dbGaP phs003690.v1.p1.
- **systemsheme/hspc** - lab-org copy of the HSPC reference work.
- **systemsheme/hrlsc2** - code for the 2026 AML LSC preprint (hrLSC2 signature).
- **systemsheme/PrismR, PowerR** - R/Shiny plotting + power-analysis tools.

Adjacent standard tool (label transfer onto a bone-marrow reference):
- **andygxzeng/BoneMarrowMap**
  https://github.com/andygxzeng/BoneMarrowMap

The outreach extension:
> Take their published MPP/OPP reference (GSE262440) and build a reproducible
> label-transfer + benchmarking notebook that projects a query AML dataset onto
> their populations (BoneMarrowMap-style), reporting mapping confidence and where
> hrLSC2-marked cells land. Connects the 2025 atlas to the 2026 LSC preprint.

Note on stack: their code is **R (Seurat/ArchR)**. Best move is an R notebook.
The Python `sceval` package here is the language-agnostic capability proof
(mapping + clonal recovery with a synthetic self-check); frame it as the eval
layer, and add an R label-transfer notebook against their real reference.

Papers (read before emailing):
- Single-cell multi-omics purifies human AML LSCs, bioRxiv 2026.07.12.737989
- A single-cell framework for adult human MPPs, Cell Reports 2025 (S2211-1247(25)01007-1)
- Single-cell genomics in AML, Blood 2023
- Single-cell mutational profiling of AML MRD, Blood Advances 2020
