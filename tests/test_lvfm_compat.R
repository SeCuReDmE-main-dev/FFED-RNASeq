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

assert_equal <- function(x, y, message) {
  if (!isTRUE(all.equal(x, y))) {
    stop(message, call. = FALSE)
  }
}

message("Testing LVFM legacy compatibility")

sample_data <- data.frame(
  sample_number = paste0("S", 1:8),
  batch = rep("B1", 8),
  replicate = paste0("R", 1:8),
  cell_line = c("AF22", "CTRL9II", "ASD17AII", "ASD12BI", "AF22", "CTRL9II", "ASD17AII", "ASD12BI"),
  exposure = paste0("E", 1:8),
  stringsAsFactors = FALSE
)

assert_equal(
  BPA_F(sample_data),
  c("control", "control", "treated", "treated", "treated", "treated", "control", "control"),
  "BPA_F mapping failed"
)
assert_equal(
  Pb_F(sample_data),
  c("control", "control", "control", "control", "treated", "treated", "treated", "treated"),
  "Pb_F mapping failed"
)
assert_equal(
  FH_F(sample_data),
  c("control", "treated", "control", "treated", "treated", "control", "treated", "control"),
  "FH_F mapping failed"
)
assert_equal(
  Eth_F(sample_data),
  c("control", "treated", "control", "treated", "control", "treated", "control", "treated"),
  "Eth_F mapping failed"
)
assert_equal(
  VPA_F(sample_data),
  c("control", "control", "treated", "treated", "control", "control", "treated", "treated"),
  "VPA_F mapping failed"
)
assert_equal(
  Zn_F(sample_data),
  c("control", "treated", "treated", "control", "control", "treated", "treated", "control"),
  "Zn_F mapping failed"
)
assert_equal(
  risk_F(sample_data),
  c("control", "control", "risk", "risk", "control", "control", "risk", "risk"),
  "risk_F mapping failed"
)

dxr <- data.frame(
  padj = c(0.01, 0.2, 0.001, 0.03),
  log2fold_treated_control = c(-1.5, 2, 0.2, 1.2),
  exonBaseMean = c(20, 30, 5, 100),
  stringsAsFactors = FALSE
)

sig <- Sig_F(dxr, pval = 0.05, log2fc = 1, ebm = 10)
assert_true(nrow(sig) == 2, "Sig_F returned the wrong number of rows")
assert_equal(sig$padj, c(0.01, 0.03), "Sig_F did not sort results by padj")

message("LVFM legacy compatibility tests passed")
