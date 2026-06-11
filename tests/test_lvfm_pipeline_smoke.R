script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grepl("^--file=", script_args)][1]
if (is.na(file_arg) || !nzchar(file_arg)) {
  stop("This test must be run with Rscript", call. = FALSE)
}

script_file <- sub("^--file=", "", file_arg)
repo_root <- normalizePath(file.path(dirname(script_file), ".."))
source(file.path(repo_root, "lvfm_core", "R", "lvfm_core.R"))

assert_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

message("Testing LVFM representative pipeline smoke path")

sample_info <- data.frame(
  sample_id = paste0("S", 1:6),
  cell_line = c("AF22", "AF22", "CTRL9II", "CTRL9II", "ASD12BI", "ASD17AII"),
  exposure = c("E1", "E5", "E2", "E7", "E3", "E8"),
  stringsAsFactors = FALSE
)

sample_info$BPA <- BPA_F(data.frame(matrix(NA, nrow = nrow(sample_info), ncol = 5), stringsAsFactors = FALSE))
sample_info$Pb <- Pb_F(data.frame(matrix(NA, nrow = nrow(sample_info), ncol = 5), stringsAsFactors = FALSE))
sample_info$risk <- risk_F(data.frame(matrix(NA, nrow = nrow(sample_info), ncol = 4), stringsAsFactors = FALSE))

counts <- data.frame(
  exon_count = c(100, 120, 80, 90, 40, 60),
  exon_signal = c(12, 18, 7, 9, 4, 6),
  exon_noise = c(1, 2, 1, 2, 3, 2)
)

rnaseq_state <- lvfm_project_rnaseq_levels(counts, numeric_cols = c("exon_count", "exon_signal", "exon_noise"))
deu_state <- lvfm_project_deu(counts, numeric_cols = c("exon_count", "exon_signal"))
enrichment_state <- lvfm_project_enrichment(counts, numeric_cols = c("exon_signal", "exon_noise"))

assert_true(inherits(rnaseq_state, "StateField"), "RNA-seq representative projection failed")
assert_true(inherits(deu_state, "StateField"), "DEU representative projection failed")
assert_true(inherits(enrichment_state, "StateField"), "enrichment representative projection failed")
assert_true(all(rnaseq_state$mu + rnaseq_state$nu <= 1), "RNA-seq representative projection not admissible")
assert_true(all(deu_state$mu + deu_state$nu <= 1), "DEU representative projection not admissible")
assert_true(all(enrichment_state$mu + enrichment_state$nu <= 1), "enrichment representative projection not admissible")

combined <- join(rnaseq_state, tensor(deu_state, enrichment_state))
assert_true(inherits(combined, "StateField"), "combined pipeline state failed")
assert_true(all(combined$mu + combined$nu <= 1), "combined pipeline state not admissible")

candidate <- lvfm_candidate_master_formula(
  phi = mean(rnaseq_state$mu),
  cubnu = mean(deu_state$nu),
  d_f_hat = mean(enrichment_state$mu),
  i_system = mean(combined$mu),
  adm = 0.9
)
assert_true(inherits(candidate$state, "StateField"), "candidate master formula pipeline linkage failed")
assert_true(all(candidate$state$mu + candidate$state$nu <= 1), "candidate master formula state not admissible")

message("LVFM representative pipeline smoke test passed")
