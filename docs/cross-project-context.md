# Cross-Project Context

This note provides lightweight context for the maintenance and integration work
currently being prepared around `FFED-RNASeq`.

## Primary Repository

- `SeCuReDmE-main-dev/FFED-RNASeq`
  - <https://github.com/SeCuReDmE-main-dev/FFED-RNASeq>

This repository keeps the original R-native analysis structure and adds a small
LVFM-oriented core for reusable lattice-valued interpretation, validation, and
compatibility wrappers.

## Related Research Prototype

- `SeCuReDmE-main-dev/FNP-QNN-MVP-version-desise-simulator-`
  - <https://github.com/SeCuReDmE-main-dev/FNP-QNN-MVP-version-desise-simulator->

This prototype explores a separate but related research lane focused on
crossmodal event encoding, testable QNN execution lanes, and disease-simulator
refactoring.

## Why Link Them

The two repositories address different layers of the same broader research
direction:

- `FFED-RNASeq` stays focused on transcriptomics workflows and reproducible
  analysis
- `FNP-QNN-MVP-version-desise-simulator-` explores modeling and comparative
  execution lanes outside the R-native analysis surface

The intent is not to merge these repositories into one runtime. The intent is
to keep a public reference showing how maintainable bioinformatics workflows and
separate modeling infrastructure can evolve in parallel.

## Scope of This Note

This is a documentation-only addition. It does not alter the statistical
pipeline, runtime dependencies, or scientific claims of `FFED-RNASeq`.
