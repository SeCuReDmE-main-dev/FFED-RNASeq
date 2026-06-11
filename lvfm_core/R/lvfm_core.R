lvfm_clamp01 <- function(x) {
  x <- as.numeric(x)
  x[x < 0] <- 0
  x[x > 1] <- 1
  x
}

lvfm_fill_na <- function(x, fill = 0.5) {
  if (!length(x)) {
    return(x)
  }
  x[is.na(x)] <- fill
  x
}

lvfm_scale01 <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) {
    return(numeric())
  }

  keep <- !is.na(x)
  if (!any(keep)) {
    return(rep(0.5, length(x)))
  }

  rng <- range(x[keep])
  if (isTRUE(all.equal(rng[1], rng[2]))) {
    out <- rep(0.5, length(x))
  } else {
    out <- (x - rng[1]) / (rng[2] - rng[1])
  }
  out[!keep] <- NA_real_
  lvfm_clamp01(out)
}

lvfm_normalize_provenance <- function(provenance = list()) {
  if (is.null(provenance)) {
    return(list())
  }
  if (is.list(provenance)) {
    return(provenance)
  }
  list(value = provenance)
}

lvfm_recycle_to_length <- function(x, n) {
  if (length(x) == n) {
    return(x)
  }
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  stop("length mismatch: expected length 1 or ", n, ", got ", length(x), call. = FALSE)
}

lvfm_state_length <- function(state) {
  if (inherits(state, "StateField")) {
    return(length(state$mu))
  }
  if (is.list(state) && !is.null(state$mu)) {
    return(length(state$mu))
  }
  stop("state must be a StateField or a list with mu/nu", call. = FALSE)
}

lvfm_validate_state_field <- function(mu, nu) {
  if (length(mu) != length(nu)) {
    stop("mu and nu must have the same length", call. = FALSE)
  }

  bad <- !is.na(mu) & !is.na(nu) & (mu + nu > 1 + 1e-8)
  if (any(bad)) {
    stop("state violates admissibility: mu + nu must be <= 1", call. = FALSE)
  }

  invisible(TRUE)
}

StateField <- function(mu, nu, provenance = list(), label = NULL) {
  mu <- lvfm_clamp01(mu)
  nu <- lvfm_clamp01(nu)
  lvfm_validate_state_field(mu, nu)

  structure(
    list(
      mu = mu,
      nu = nu,
      pi = lvfm_clamp01(1 - mu - nu),
      provenance = lvfm_normalize_provenance(provenance),
      label = label
    ),
    class = "StateField"
  )
}

StepOperator <- function(step, name = "step_operator", parameters = list()) {
  if (!is.function(step)) {
    stop("step must be a function", call. = FALSE)
  }

  structure(
    list(
      step = step,
      name = name,
      parameters = lvfm_normalize_provenance(parameters)
    ),
    class = "StepOperator"
  )
}

ObservableProjector <- function(project, name = "observable_projector", metadata = list()) {
  if (!is.function(project)) {
    stop("project must be a function", call. = FALSE)
  }

  structure(
    list(
      project = project,
      name = name,
      metadata = lvfm_normalize_provenance(metadata)
    ),
    class = "ObservableProjector"
  )
}

lvfm_as_state_field <- function(x, label = NULL) {
  if (inherits(x, "StateField")) {
    return(x)
  }
  if (is.list(x) && !is.null(x$mu) && !is.null(x$nu)) {
    return(StateField(x$mu, x$nu, provenance = x$provenance, label = if (is.null(label)) x$label else label))
  }
  stop("cannot coerce object to StateField", call. = FALSE)
}

lvfm_recycle_state <- function(state, n) {
  state <- lvfm_as_state_field(state)
  if (length(state$mu) == n) {
    return(state)
  }
  StateField(
    lvfm_recycle_to_length(state$mu, n),
    lvfm_recycle_to_length(state$nu, n),
    provenance = state$provenance,
    label = state$label
  )
}

lvfm_apply_step <- function(state, operator, ...) {
  state <- lvfm_as_state_field(state)
  if (!inherits(operator, "StepOperator")) {
    stop("operator must be a StepOperator", call. = FALSE)
  }
  operator$step(state, ...)
}

lvfm_apply_projector <- function(projector, data, ...) {
  if (!inherits(projector, "ObservableProjector")) {
    stop("projector must be an ObservableProjector", call. = FALSE)
  }
  projector$project(data, ...)
}

print.StateField <- function(x, ...) {
  cat("StateField <", length(x$mu), ">\n", sep = "")
  cat("  mu: ", paste(utils::head(signif(x$mu, 4), 4), collapse = ", "), "\n", sep = "")
  cat("  nu: ", paste(utils::head(signif(x$nu, 4), 4), collapse = ", "), "\n", sep = "")
  cat("  pi: ", paste(utils::head(signif(x$pi, 4), 4), collapse = ", "), "\n", sep = "")
  invisible(x)
}

print.StepOperator <- function(x, ...) {
  cat("StepOperator <", x$name, ">\n", sep = "")
  invisible(x)
}

print.ObservableProjector <- function(x, ...) {
  cat("ObservableProjector <", x$name, ">\n", sep = "")
  invisible(x)
}

lvfm_join <- function(a, b) {
  a <- lvfm_as_state_field(a)
  b <- lvfm_as_state_field(b)
  n <- max(length(a$mu), length(b$mu))
  a <- lvfm_recycle_state(a, n)
  b <- lvfm_recycle_state(b, n)

  StateField(
    pmax(a$mu, b$mu),
    pmin(a$nu, b$nu),
    provenance = list(
      operation = "join",
      left = a$provenance,
      right = b$provenance
    ),
    label = "join"
  )
}

lvfm_meet <- function(a, b) {
  a <- lvfm_as_state_field(a)
  b <- lvfm_as_state_field(b)
  n <- max(length(a$mu), length(b$mu))
  a <- lvfm_recycle_state(a, n)
  b <- lvfm_recycle_state(b, n)

  StateField(
    pmin(a$mu, b$mu),
    pmax(a$nu, b$nu),
    provenance = list(
      operation = "meet",
      left = a$provenance,
      right = b$provenance
    ),
    label = "meet"
  )
}

lvfm_tensor <- function(a, b) {
  a <- lvfm_as_state_field(a)
  b <- lvfm_as_state_field(b)
  n <- max(length(a$mu), length(b$mu))
  a <- lvfm_recycle_state(a, n)
  b <- lvfm_recycle_state(b, n)

  mu <- lvfm_clamp01(a$mu * b$mu)
  nu_candidate <- a$nu + b$nu - (a$nu * b$nu)
  nu <- pmin(lvfm_clamp01(nu_candidate), pmax(0, 1 - mu))

  StateField(
    mu,
    nu,
    provenance = list(
      operation = "tensor",
      left = a$provenance,
      right = b$provenance
    ),
    label = "tensor"
  )
}

join <- lvfm_join
meet <- lvfm_meet
tensor <- lvfm_tensor

lvfm_project_table <- function(data, numeric_cols = NULL, label = NULL, provenance = list(), fallback_mu = 0.5, fallback_nu = 0.5) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }

  n <- nrow(data)
  if (n == 0L) {
    return(StateField(numeric(), numeric(), provenance = provenance, label = label))
  }

  if (is.null(numeric_cols)) {
    numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  }
  numeric_cols <- intersect(numeric_cols, names(data))

  if (!length(numeric_cols)) {
    return(StateField(rep(fallback_mu, n), rep(fallback_nu, n), provenance = provenance, label = label))
  }

  matrix_list <- lapply(data[numeric_cols], lvfm_scale01)
  numeric_matrix <- do.call(cbind, matrix_list)
  if (is.null(dim(numeric_matrix))) {
    numeric_matrix <- matrix(numeric_matrix, ncol = 1L)
  }

  mu <- rowMeans(numeric_matrix, na.rm = TRUE)
  mu[is.nan(mu) | is.na(mu)] <- fallback_mu

  nu <- rowMeans(1 - numeric_matrix, na.rm = TRUE)
  nu[is.nan(nu) | is.na(nu)] <- fallback_nu

  missing_rate <- rowMeans(is.na(numeric_matrix))
  missing_rate[is.nan(missing_rate) | is.na(missing_rate)] <- 1

  nu <- lvfm_clamp01(nu + 0.25 * missing_rate)
  nu <- pmin(nu, pmax(0, 1 - mu))

  StateField(
    mu,
    nu,
    provenance = c(
      list(source_columns = numeric_cols, rows = n),
      lvfm_normalize_provenance(provenance)
    ),
    label = label
  )
}

lvfm_project_metadata <- function(data, column, positive = NULL, negative = NULL, label = NULL, provenance = list()) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  if (is.numeric(column)) {
    column <- names(data)[column]
  }
  if (!column %in% names(data)) {
    stop("column not found in data", call. = FALSE)
  }

  values <- as.character(data[[column]])
  n <- length(values)
  if (!length(values)) {
    return(StateField(numeric(), numeric(), provenance = provenance, label = label))
  }

  if (!is.null(positive) || !is.null(negative)) {
    mu <- if (is.null(positive)) rep(0, n) else as.numeric(values %in% positive)
    nu <- if (is.null(negative)) rep(0, n) else as.numeric(values %in% negative)
    nu <- pmin(nu, pmax(0, 1 - mu))
    return(StateField(mu, nu, provenance = c(list(column = column, positive = positive, negative = negative), lvfm_normalize_provenance(provenance)), label = label))
  }

  codes <- as.numeric(factor(values, levels = unique(values)))
  mu <- lvfm_scale01(codes)
  nu <- pmax(0, 1 - mu)
  StateField(
    mu,
    nu,
    provenance = c(list(column = column, levels = unique(values)), lvfm_normalize_provenance(provenance)),
    label = label
  )
}

lvfm_project_rnaseq_levels <- function(data, numeric_cols = NULL, provenance = list()) {
  lvfm_project_table(data, numeric_cols = numeric_cols, label = "rnaseq_levels_i_to_iii", provenance = provenance)
}

lvfm_project_rnaseq_level_iv <- function(data, numeric_cols = NULL, provenance = list()) {
  lvfm_project_table(data, numeric_cols = numeric_cols, label = "rnaseq_level_iv", provenance = provenance)
}

lvfm_project_deu <- function(data, numeric_cols = NULL, provenance = list()) {
  lvfm_project_table(data, numeric_cols = numeric_cols, label = "deu", provenance = provenance)
}

lvfm_project_enrichment <- function(data, numeric_cols = NULL, provenance = list()) {
  lvfm_project_table(data, numeric_cols = numeric_cols, label = "enrichment", provenance = provenance)
}

lvfm_project_brainspan <- function(data, numeric_cols = NULL, provenance = list()) {
  lvfm_project_table(data, numeric_cols = numeric_cols, label = "brainspan", provenance = provenance)
}

lvfm_project_metabolomics <- function(data, numeric_cols = NULL, provenance = list()) {
  lvfm_project_table(data, numeric_cols = numeric_cols, label = "metabolomics", provenance = provenance)
}

lvfm_legacy_factor_map <- function(SampleData, column, mapping, default = NA_character_) {
  if (!is.data.frame(SampleData)) {
    stop("SampleData must be a data.frame", call. = FALSE)
  }

  if (is.numeric(column)) {
    column <- column[[1L]]
  }

  if (is.numeric(column)) {
    values <- SampleData[[column]]
  } else {
    if (!column %in% names(SampleData)) {
      stop("column not found in SampleData", call. = FALSE)
    }
    values <- SampleData[[column]]
  }

  out <- rep(default, nrow(SampleData))
  mapped <- unname(mapping[as.character(values)])
  keep <- !is.na(mapped)
  out[keep] <- mapped[keep]
  out
}

BPA_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    5,
    c(
      E1 = "control",
      E2 = "control",
      E3 = "treated",
      E4 = "treated",
      E5 = "treated",
      E6 = "treated",
      E7 = "control",
      E8 = "control"
    )
  )
}

Pb_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    5,
    c(
      E1 = "control",
      E2 = "control",
      E3 = "control",
      E4 = "control",
      E5 = "treated",
      E6 = "treated",
      E7 = "treated",
      E8 = "treated"
    )
  )
}

FH_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    5,
    c(
      E1 = "control",
      E2 = "treated",
      E3 = "control",
      E4 = "treated",
      E5 = "treated",
      E6 = "control",
      E7 = "treated",
      E8 = "control"
    )
  )
}

Eth_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    5,
    c(
      E1 = "control",
      E2 = "treated",
      E3 = "control",
      E4 = "treated",
      E5 = "control",
      E6 = "treated",
      E7 = "control",
      E8 = "treated"
    )
  )
}

VPA_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    5,
    c(
      E1 = "control",
      E2 = "control",
      E3 = "treated",
      E4 = "treated",
      E5 = "control",
      E6 = "control",
      E7 = "treated",
      E8 = "treated"
    )
  )
}

Zn_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    5,
    c(
      E1 = "control",
      E2 = "treated",
      E3 = "treated",
      E4 = "control",
      E5 = "control",
      E6 = "treated",
      E7 = "treated",
      E8 = "control"
    )
  )
}

risk_F <- function(SampleData) {
  lvfm_legacy_factor_map(
    SampleData,
    4,
    c(
      AF22 = "control",
      CTRL9II = "control",
      ASD17AII = "risk",
      ASD12BI = "risk"
    )
  )
}

Sig_F <- function(data, pval, log2fc, ebm) {
  required <- c("padj", "log2fold_treated_control", "exonBaseMean")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  sig <- as.data.frame(data[!is.na(data$padj) & data$padj < pval, , drop = FALSE])
  sig <- sig[sig$log2fold_treated_control < -log2fc | sig$log2fold_treated_control > log2fc, , drop = FALSE]
  sig <- sig[sig$exonBaseMean >= ebm, , drop = FALSE]
  sig <- sig[order(sig$padj), , drop = FALSE]
  sig
}

lvfm_score_to_state <- function(score, adm = 1, label = "candidate_master_formula", provenance = list()) {
  score <- as.numeric(score)
  adm <- lvfm_clamp01(adm)
  mu <- stats::plogis(score)
  nu <- pmax(0, (1 - mu) * (1 - adm))
  StateField(
    mu,
    nu,
    provenance = c(list(score = score, adm = adm), lvfm_normalize_provenance(provenance)),
    label = label
  )
}

lvfm_candidate_master_formula <- function(phi, cubnu, d_f_hat, i_system, adm, weights = list(), offset = 0, return_state = TRUE) {
  weight_value <- function(name, default = 1) {
    if (!is.list(weights) && !is.atomic(weights)) {
      return(default)
    }
    if (!(name %in% names(weights))) {
      return(default)
    }
    as.numeric(weights[[name]])
  }

  phi3 <- phi^3
  score <- offset +
    weight_value("phi", 1) * phi3 +
    weight_value("cubnu", 1) * cubnu +
    weight_value("d_f_hat", 1) * d_f_hat +
    weight_value("i_system", 1) * i_system +
    weight_value("adm", 1) * adm

  formula_text <- paste(
    "candidate: score = w_phi * phi^3 + w_cubnu * CubNu + w_d_f_hat * D_f_hat + w_i_system * I_system + w_adm * Adm",
    sep = ""
  )

  result <- list(
    formula = formula_text,
    terms = list(
      phi3 = phi3,
      cubnu = cubnu,
      d_f_hat = d_f_hat,
      i_system = i_system,
      adm = adm
    ),
    weights = list(
      phi = weight_value("phi", 1),
      cubnu = weight_value("cubnu", 1),
      d_f_hat = weight_value("d_f_hat", 1),
      i_system = weight_value("i_system", 1),
      adm = weight_value("adm", 1)
    ),
    score = score,
    provenance = list(role = "candidate_master_formula")
  )

  if (isTRUE(return_state)) {
    result$state <- lvfm_score_to_state(score, adm = adm, provenance = result$provenance)
  }

  result
}
