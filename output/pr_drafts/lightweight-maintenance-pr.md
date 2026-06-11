# PR Draft: Lightweight LVFM Integration

## Title
Add `lvfm_core` as a lightweight LVFM integration layer

## Context

This PR is intentionally small. It is meant to show maintenance activity,
preserve the current transcriptomics workflows, and add a reusable LVFM
interpretation layer without turning the repository into a rewrite project.

The repository has been quiet for a long time, so the goal here is to land a
low-risk maintenance patch first and keep the bigger integration work strictly
adapter-first.

## What This PR Adds

- `lvfm_core/` as a package-style R core
- `StateField`, `StepOperator`, and `ObservableProjector`
- lattice operators:
  - `join`
  - `meet`
  - `tensor`
- pipeline adapters for:
  - RNA-seq levels I to III
  - RNA-seq level IV
  - DEU
  - gene enrichment
  - BrainSpan
  - metabolomics
- legacy wrappers for the existing `*_F.R` helpers and `Sig_F`
- public README refresh so the repo description matches the real structure
- a short research brief documenting the candidate Master Formula work
- smoke tests for the LVFM core and compatibility layer

## Why This Is a Good First PR

- It does not disturb the current analysis outputs.
- It keeps the old helper surface alive.
- It adds test coverage before any broader refactor.
- It gives maintainers a clear, reviewable delta instead of a long RFC.

## Related Context

- Issue #1 is vague and only signals that maintainers may still be reachable.
- This PR should be submitted as a maintenance-oriented follow-up, not as a
  proposal for a full project rewrite.

## Validation

- LVFM core tests pass locally with `Rscript`.
- Legacy compatibility tests pass locally with `Rscript`.
- The DEU notebooks now source the shared `lvfm_core` file instead of
  duplicating the factor mapping logic.

## Follow-Up After Merge

1. expand adapter coverage to the remaining pipeline notebooks if needed
2. keep the scientific/statistical outputs unchanged unless a maintainer asks
   for behavior changes
3. reserve the larger integration discussion for a later PR if there is
   interest

