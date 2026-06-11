# Objective

Integrate a new LVFM research layer into the existing FFED-RNASeq repository
without breaking the current R analysis pipelines. The goal is to add a
package-style `lvfm_core` that can interpret and validate outputs from the
existing transcriptomics workflows.

This pass also evaluates whether `AlphaTensor` adds useful design guidance for
the project.

# Environment / Stack Context

- Repository: `C:\Users\jeans\Desktop\Case study\modele\FFED-RNAeq\FFED-RNASeq`
- Primary language: R
- Existing pipelines:
  - RNA-seq levels I to III
  - RNA-seq level IV interactions
  - DEU
  - gene enrichment
  - BrainSpan enrichment
  - metabolomics
- Legacy helper surface: `*_F.R` and `Sig_F`
- New core folder: `lvfm_core/`

# Research Questions

1. What is the correct role of LVFM in this repository?
2. Which parts of the legacy DEU helper surface should be preserved?
3. Can the new core stay R-native and still provide a stable contract?
4. Does `AlphaTensor` change the plan in a useful way?

# Findings

- The repository is already a valid R analysis project with multiple stable
  pipelines. The safest integration path is adapter-first, not a monolithic
  rewrite.
- The LVFM layer fits best as an interpretation and validation core:
  - `StateField` for `mu`, `nu`, `pi`
  - `StepOperator` for iterative transitions
  - `ObservableProjector` for mapping data tables to lattice-valued states
  - `join`, `meet`, `tensor` as lattice operators
- The legacy DEU helper functions are simple and should remain available as
  compatibility wrappers rather than being removed outright.
- The candidate master formula should stay an object of work, not a claim of
  truth. The current working form is a scored combination of `phi^3`,
  `CubNu`, `D_f_hat`, `I_system`, and `Adm`.
- `AlphaTensor` is useful as a research reference for tensor decomposition and
  benchmark discipline. It is not a runtime dependency and does not fit the
  current R-native transcriptomics stack directly.

# Recommended Path

1. Keep the existing FFED-RNASeq pipelines intact.
2. Use `lvfm_core` as a shared R layer for lattice projection, validation, and
   compatibility wrappers.
3. Source `lvfm_core/R/lvfm_core.R` from the DEU notebooks instead of the
   scattered `*_F.R` files.
4. Keep the candidate master formula as a documented research object with a
   reproducible scoring path.
5. Treat `AlphaTensor` as an external idea source for later tensor search or
   decomposition work, not as an implementation dependency.

# Alternatives Considered

- Full rewrite of the repository into a new monolithic package
  - rejected because it would risk the current scientific outputs
- Importing `AlphaTensor` directly into the runtime
  - rejected because it is a different problem domain and a different language
- Leaving the legacy helpers untouched
  - rejected because it preserves duplication and blocks a clean LVFM contract

# Risks / Unknowns

- The current LVFM formula is still a candidate, not a validated law.
- Some notebooks still contain old explanatory text and may need cleanup after
  the core is sourced consistently.
- External tools or future benchmarks may require a separate adapter if tensor
  search becomes a real task.

# Sources

- Local repository files:
  - `README.md`
  - `DEU_Analysis/RMarkdowns/*.Rmd`
  - `DEU_Analysis/R_codes/*.R`
  - `lvfm_core/R/lvfm_core.R`
- Google Doc:
  - `https://docs.google.com/document/d/1yLObOE6pJrCofoFUgKlgHDCsR9fYKEjWqHMPzfF_QjM/edit?usp=sharing`
- AlphaTensor:
  - `https://github.com/google-deepmind/alphatensor`
  - `https://www.nature.com/articles/s41586-022-05172-4`

