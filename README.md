# FFED-RNASeq

This repository contains the R-native analysis code for the FFED RNA-seq study.
The original analysis pipelines remain the primary runtime surface. The new LVFM
layer is an added research and validation core that sits on top of the existing
transcriptomics workflows.

## What Is Here

- `RNAseq_Analysis_Levels_I_II_III/`
  - differential gene expression and gene-set enrichment for the independent
    effect levels I to III
- `RNAseq_Analysis_Level_IV/`
  - differential expression and enrichment for interaction effects
- `DEU_Analysis/`
  - differential exon usage analysis, now wired to the shared LVFM core
- `Gene_Enrichment_Analysis/`
  - enrichment against public gene lists
- `BrainSpan_Analysis/`
  - BrainSpan-based enrichment checks for developmental expression patterns
- `Metabolomics_Analysis/`
  - metabolomics analysis for treated vs control comparisons
- `lvfm_core/`
  - the new package-style LVFM core with lattice operators, adapters, and
    legacy compatibility wrappers
- `tests/`
  - smoke and compatibility tests for the new core

## Current Architecture

The repository is now organized around two layers:

1. the existing analysis pipelines, which remain the scientific backbone
2. the `lvfm_core` layer, which provides lattice-valued interpretation and
   validation helpers for those pipelines

The LVFM core includes:

- `StateField` for `mu`, `nu`, `pi`, and provenance
- `StepOperator` for iterative updates
- `ObservableProjector` for converting tables into lattice-valued states
- `join`, `meet`, and `tensor` operators
- pipeline adapters for RNA-seq, DEU, enrichment, BrainSpan, and metabolomics
- compatibility wrappers for the legacy `*_F.R` helpers and `Sig_F`
- a candidate master-formula object based on `phi^3`, `CubNu`, `D_f_hat`,
  `I_system`, and `Adm`

## LVFM Research Note

The LVFM layer is experimental. It formalizes a candidate working formula and
an admissible lattice model, but it does not replace the existing statistical
analysis outputs.

The initial lattice used here is the intuitionistic slice:

- `L* = {(mu, nu) in [0,1]^2 : mu + nu <= 1}`
- `pi = 1 - mu - nu`
- `join(mu, nu) = (max(mu), min(nu))`
- `meet(mu, nu) = (min(mu), max(nu))`

## AlphaTensor Reference

`AlphaTensor` is treated as an external research reference, not as a runtime
dependency.

Why it matters here:

- it shows how tensor decomposition can be turned into a search and benchmark
  problem
- it is useful as a design reference if the LVFM line later needs operator
  search, tensor factorization experiments, or decomposition benchmarks
- it is not a fit for the current R-native pipeline because the repo here is
  centered on transcriptomics analysis, not matrix multiplication optimization

## Legacy Compatibility

The old factor-mapping helpers are preserved as wrappers:

- `BPA_F`
- `Pb_F`
- `FH_F`
- `Eth_F`
- `VPA_F`
- `Zn_F`
- `risk_F`
- `Sig_F`

The DEU notebooks now source the shared LVFM core instead of duplicating the
legacy helper files.

## How To Run The New Tests

From the repository root:

```bash
Rscript tests/test_lvfm_core.R
Rscript tests/test_lvfm_compat.R
```

## Scope

This repository is for research and reproducible analysis. It does not claim
clinical validation, medical performance, or life-saving capability.

