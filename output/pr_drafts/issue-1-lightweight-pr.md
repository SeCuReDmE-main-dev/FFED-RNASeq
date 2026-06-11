# Draft PR for Issue #1

## Title
Introduce `lvfm_core` as a lightweight LVFM layer for FFED-RNASeq

## Summary

This PR adds a small, adapter-first LVFM layer on top of the existing R
analysis pipelines. It does not rewrite the repository and does not change the
core statistical workflows.

## What Changed

- Added a new package-style core under `lvfm_core/`
- Defined the minimal contract:
  - `StateField`
  - `StepOperator`
  - `ObservableProjector`
- Added lattice operators:
  - `join`
  - `meet`
  - `tensor`
- Kept the legacy factor helpers as compatibility wrappers:
  - `BPA_F`
  - `Pb_F`
  - `FH_F`
  - `Eth_F`
  - `VPA_F`
  - `Zn_F`
  - `risk_F`
  - `Sig_F`
- Redirected the DEU notebooks to the shared LVFM core
- Rewrote the public README so it reflects the current repository shape
- Added a research brief documenting the candidate Master Formula and the
  integration rationale
- Added smoke/compatibility tests for the new core

## Why This Is Small

- No RFC-style redesign
- No monolithic migration
- No change to the underlying RNA-seq/DEU/statistical pipeline outputs
- The LVFM layer is additive and can be ignored by users who only need the
  original analysis flow

## Issue Context

Related issue: #1

The issue body is vague and sounds like a request to confirm that the
repository is still maintained and that communication with maintainers is
possible. This PR draft keeps the scope narrow and low-risk so it can be sent as
a concise maintenance-focused follow-up.

## Validation

- Static inspection completed
- Documentation updated
- Core contract and compatibility tests added

## Follow-up If Maintainers Respond

1. run the R tests in a local runtime
2. verify one representative DEU notebook end-to-end
3. keep the PR small unless maintainers request a broader integration pass

