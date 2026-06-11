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

message("Testing LVFM core")

sf <- StateField(c(0.2, 0.6), c(0.3, 0.1), provenance = list(source = "unit"), label = "state")
assert_true(all(sf$mu + sf$nu <= 1), "StateField admissibility failed")
assert_equal(sf$pi, c(0.5, 0.3), "StateField pi computation failed")

a <- StateField(0.2, 0.5, provenance = list(id = "a"))
b <- StateField(0.7, 0.1, provenance = list(id = "b"))
j <- join(a, b)
m <- meet(a, b)
t <- tensor(a, b)

assert_equal(j$mu, 0.7, "join mu failed")
assert_equal(j$nu, 0.1, "join nu failed")
assert_equal(m$mu, 0.2, "meet mu failed")
assert_equal(m$nu, 0.5, "meet nu failed")
assert_true(all(t$mu + t$nu <= 1), "tensor admissibility failed")
assert_equal(t$mu, 0.14, "tensor mu failed")
assert_equal(round(t$nu, 4), 0.55, "tensor nu failed")

tab <- data.frame(
  feature_a = c(1, 2, 3),
  feature_b = c(10, 20, 30),
  feature_c = c(5, NA, 15)
)
projected <- lvfm_project_table(tab, label = "test_table")
assert_true(inherits(projected, "StateField"), "table projection did not produce a StateField")
assert_true(length(projected$mu) == nrow(tab), "projected mu length mismatch")
assert_true(all(projected$mu + projected$nu <= 1 + 1e-8, na.rm = TRUE), "projected state is not admissible")

meta <- data.frame(group = c("risk", "control", "risk", "unknown"), stringsAsFactors = FALSE)
meta_state <- lvfm_project_metadata(meta, "group", positive = "risk", negative = "control", label = "meta")
assert_equal(meta_state$mu, c(1, 0, 1, 0), "metadata positive mapping failed")
assert_equal(meta_state$nu, c(0, 1, 0, 0), "metadata negative mapping failed")

candidate <- lvfm_candidate_master_formula(
  phi = 0.4,
  cubnu = 0.3,
  d_f_hat = 0.2,
  i_system = 0.1,
  adm = 0.8
)
assert_true(!is.null(candidate$state), "candidate formula did not return a state")
assert_true(all(candidate$state$mu + candidate$state$nu <= 1), "candidate state is not admissible")
assert_true(grepl("phi\\^3", candidate$formula), "candidate formula text missing phi^3")

step <- StepOperator(function(state) {
  StateField(state$mu * 0.5, pmin(state$nu + 0.1, 1 - state$mu * 0.5), provenance = state$provenance, label = "stepped")
}, name = "halving_step")
stepped <- lvfm_apply_step(a, step)
assert_true(inherits(stepped, "StateField"), "step operator failed")

projector <- ObservableProjector(function(data) {
  lvfm_project_table(data, label = "projected_via_operator")
}, name = "table_projector")
projected_via_operator <- lvfm_apply_projector(projector, tab)
assert_true(inherits(projected_via_operator, "StateField"), "projector failed")

message("LVFM core tests passed")
