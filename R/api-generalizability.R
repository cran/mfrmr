# ==============================================================================
# Generalizability theory adapter (random rater / random item effects)
# ==============================================================================
#
# `mfrm_generalizability()` formalises the generalizability-theory
# (Cronbach et al. 1972 / Brennan 2001) decomposition that
# `compute_facet_icc()` already performs internally via
# `lme4::lmer`. The MFRM frames non-person facets as fixed effects
# for measurement; for reporting a G-coefficient, the same facets
# are re-fit as crossed random effects so the variance components
# can be combined into the canonical G / Phi indices.
#
# This helper does NOT re-fit the MFRM. It treats the rating data
# as a crossed random-effects ANOVA (Person + each non-person facet
# + residual) and returns the variance decomposition + G / Phi
# coefficients.

#' Generalizability-theory variance decomposition for an MFRM design
#'
#' Re-fits the rating data underlying an `mfrm_fit` as a crossed
#' random-effects model
#' `Score ~ 1 + (1 | Person) + (1 | Facet1) + ... + Residual`
#' via `lme4::lmer`, and returns the canonical G-theory variance
#' components plus G / Phi coefficients. Useful when reviewers ask
#' for a generalizability-theory complement to the Rasch-style
#' separation / reliability statistics that `diagnose_mfrm()`
#' already emits.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param data Optional data frame. When `NULL`, the rating data
#'   stored on `fit$prep$data` is used.
#' @param object_facet Facet that plays the role of the "object of
#'   measurement" -- typically `"Person"` (default).
#' @param random_facets Character vector of non-person facets to
#'   treat as random conditions of measurement. Default uses every
#'   facet other than `object_facet`.
#' @param reml Logical, passed to [lme4::lmer()] (default `TRUE`).
#'
#' @return An object of class `mfrm_generalizability` with:
#' \describe{
#'   \item{`variance_components`}{One row per random effect plus
#'     residual, with columns `Source`, `Variance`, and
#'     `ProportionVariance`.}
#'   \item{`coefficients`}{One-row data frame with `G`
#'     (generalizability coefficient, relative decision) and
#'     `Phi` (dependability coefficient, absolute decision), using the
#'     single-observation-per-cell convention.}
#'   \item{`design`}{Description of the crossed-random model.}
#' }
#'
#' @section Interpretation:
#' - `G` is appropriate for **relative** decisions (rank-ordering
#'   persons): `G = sigma2(p) / (sigma2(p) + sigma2(Residual))`.
#' - The reported `Phi` is appropriate for **absolute** decisions (cut-score
#'   classification): `Phi = sigma2(p) / (sigma2(p) + sigma2(facet
#'   main effects) + sigma2(Residual))`, before D-study scaling.
#' - Use [mfrm_d_study()] to project `G` / `Phi` under planned numbers of
#'   raters, items, criteria, or other random measurement facets.
#' - Reporting bands follow Brennan (2001): G / Phi >= 0.8 for
#'   high-stakes decisions, >= 0.7 for routine reporting.
#'
#' @section Limitations:
#' This helper formulates the random-effects model with main effects
#' only (`Score ~ 1 + (1|Person) + (1|Facet1) + ... + Residual`); no
#' explicit `(1 | Person:Rater)`, `(1 | Person:Criterion)`, or
#' `(1 | Rater:Criterion)` interaction terms are estimated. All
#' two-way and higher interaction variance is therefore folded into
#' the `Residual` term -- the standard one-observation-per-cell
#' approximation -- which can bias `G` downward when person x facet
#' interactions are substantively large. This function reports the
#' one-observation-per-cell baseline. [mfrm_d_study()] applies D-study
#' scaling, including residual-scaling sensitivity checks, to the same
#' simplified variance-component decomposition.
#' Because person-by-facet interaction terms are not estimated separately,
#' D-study projections remain practical planning evidence rather than a
#' replacement for a fully specified G-theory design.
#'
#' @section References:
#' - Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N.
#'   (1972). *The dependability of behavioral measurements: Theory
#'   of generalizability for scores and profiles*. Wiley.
#' - Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [mfrm_d_study()], [compute_facet_icc()], [diagnose_mfrm()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   gt <- mfrm_generalizability(fit)
#'   gt$variance_components
#'   # Look for: a Person variance component well above any single
#'   #   non-person facet's variance share. Large rater or criterion
#'   #   variance shares mean those conditions add measurement error
#'   #   relative to person spread.
#'   gt$coefficients
#'   # Look for: G >= 0.7 for routine reporting, >= 0.8 for high-stakes.
#'   #   G < Phi means absolute decisions are noisier than relative
#'   #   decisions; review whether facet main effects need anchoring.
#' }
#' }
#' @export
mfrm_generalizability <- function(fit,
                                  data = NULL,
                                  object_facet = "Person",
                                  random_facets = NULL,
                                  reml = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("`mfrm_generalizability()` requires the `lme4` package ",
         "(in Suggests). Install it and retry.", call. = FALSE)
  }
  if (is.null(data)) {
    data <- as.data.frame(fit$prep$data %||% data.frame(),
                          stringsAsFactors = FALSE)
  }
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` is empty or not a data frame.", call. = FALSE)
  }

  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (is.null(random_facets)) {
    random_facets <- setdiff(facet_names, object_facet)
  }
  random_facets <- as.character(random_facets)
  if (length(random_facets) == 0L) {
    stop("At least one non-person facet is required as a random ",
         "condition of measurement.", call. = FALSE)
  }
  needed_cols <- c(object_facet, random_facets, "Score")
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Data frame is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }

  for (col in c(object_facet, random_facets)) {
    if (!is.factor(data[[col]])) {
      data[[col]] <- as.factor(as.character(data[[col]]))
    }
  }
  data$Score <- suppressWarnings(as.numeric(data$Score))
  data <- data[is.finite(data$Score), , drop = FALSE]

  random_terms <- c(object_facet, random_facets)
  formula_str <- paste0(
    "Score ~ 1 + ",
    paste0("(1 | ", random_terms, ")", collapse = " + ")
  )
  formula <- stats::as.formula(formula_str)

  lmer_warnings <- character(0)
  fit_lmer <- tryCatch(
    withCallingHandlers(
      lme4::lmer(formula, data = data, REML = isTRUE(reml)),
      warning = function(w) {
        lmer_warnings <<- c(lmer_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit_lmer, "error")) {
    stop("lme4::lmer failed: ", conditionMessage(fit_lmer), call. = FALSE)
  }

  vc <- as.data.frame(lme4::VarCorr(fit_lmer))
  vc <- vc[is.na(vc$var2), c("grp", "vcov")]
  total_var <- sum(vc$vcov, na.rm = TRUE)
  zero_tol <- sqrt(.Machine$double.eps)
  var_components <- data.frame(
    Source = as.character(vc$grp),
    Variance = round(as.numeric(vc$vcov), 6),
    ProportionVariance = if (is.finite(total_var) && total_var > zero_tol) {
      round(vc$vcov / total_var, 4)
    } else {
      rep(NA_real_, length(vc$vcov))
    },
    stringsAsFactors = FALSE
  )

  # G / Phi coefficients (single observation per cell convention).
  # Following Brennan (2001) for a crossed p x i design:
  #   G  = sigma2(p) / (sigma2(p) + sigma2(pi))
  #   Phi= sigma2(p) / (sigma2(p) + sigma2(i) + sigma2(pi))
  # Without an explicit interaction term in the formula above the
  # interaction variance is folded into the Residual term, which is
  # the standard one-observation-per-cell approximation.
  v <- stats::setNames(var_components$Variance, var_components$Source)
  sigma2_p <- as.numeric(v[object_facet] %||% NA_real_)
  sigma2_residual <- as.numeric(v["Residual"] %||% NA_real_)
  sigma2_main <- if (length(random_facets) > 0L) {
    sum(as.numeric(v[random_facets] %||% 0), na.rm = TRUE)
  } else 0
  if (!is.finite(sigma2_p) || sigma2_p <= 0) {
    G_coef <- NA_real_
    Phi_coef <- NA_real_
  } else {
    G_coef <- sigma2_p / (sigma2_p + sigma2_residual)
    Phi_coef <- sigma2_p / (sigma2_p + sigma2_main + sigma2_residual)
  }

  out <- list(
    variance_components = var_components,
    coefficients = data.frame(
      G = round(G_coef, 4),
      Phi = round(Phi_coef, 4),
      stringsAsFactors = FALSE
    ),
    design = list(
      object_facet = object_facet,
      random_facets = random_facets,
      observed_levels = stats::setNames(
        vapply(c(object_facet, random_facets), function(col) {
          dplyr::n_distinct(data[[col]])
        }, integer(1)),
        c(object_facet, random_facets)
      ),
      formula = format(formula),
      reml = isTRUE(reml),
      lmer_warnings = lmer_warnings
    )
  )
  class(out) <- c("mfrm_generalizability", "list")
  out
}

#' Project G-theory coefficients under alternative D-study designs
#'
#' @description
#' `mfrm_d_study()` applies a practical D-study projection to the
#' variance components from [mfrm_generalizability()]. It answers questions such
#' as "what happens to `G` and `Phi` if we use 2, 3, or 4 raters?" without
#' re-fitting the Rasch/MFRM model.
#'
#' @param x Output from [mfrm_generalizability()] or an `mfrm_fit`. If an
#'   `mfrm_fit` is supplied, [mfrm_generalizability()] is called first.
#' @param design_grid Data frame or named list giving planned counts for each
#'   random measurement facet. Column names may be the facet names themselves
#'   (for example `Rater`) or `n_` plus the facet name (for example
#'   `n_Rater`). When `NULL`, one row using the observed number of levels is
#'   returned.
#' @param object_facet,random_facets Passed to [mfrm_generalizability()] when
#'   `x` is an `mfrm_fit`.
#' @param residual_scaling How the collapsed residual variance should be scaled
#'   when planned facet counts increase. `"highest_order"` treats the residual
#'   as highest-order person-by-all-conditions/error variance and divides by
#'   the product of planned counts. `"single_condition"` divides by the smallest
#'   planned facet count, a conservative sensitivity check when unmodeled
#'   person-by-one-facet interactions may dominate. `"none"` leaves the residual
#'   unscaled. `"sensitivity"` returns all three assumptions for each design
#'   row.
#' @param ... Additional arguments passed to [mfrm_generalizability()] when `x`
#'   is an `mfrm_fit`.
#'
#' @details
#' The projection uses the variance decomposition already estimated by
#' [mfrm_generalizability()]. For a random measurement facet `j`, main-effect
#' variance contributes `sigma2_j / n_j` to the absolute-error denominator.
#' The residual term contains unmodeled person-by-facet and higher-order
#' interaction variance in the current simplified G-study, so the selected
#' `residual_scaling` assumption is reported explicitly. The relative-decision
#' denominator uses only this scaled residual term.
#'
#' This is a pragmatic D-study planning layer, not a full p x r x i ANOVA
#' decomposition. If person-by-rater or person-by-item interactions are a
#' primary estimand, use `residual_scaling = "sensitivity"` and treat the output
#' as planning evidence; fit a fully crossed G-theory model externally when
#' those interaction components must be estimated separately.
#'
#' The `G` and `Phi` values returned here belong to the generalizability-theory
#' metric family. They should not be interpreted as coefficient alpha, omega,
#' KR-20, or IRT marginal/separation reliability, even though all of those
#' summaries may be displayed on a 0--1 scale in broader reporting dashboards.
#'
#' @return An object of class `mfrm_d_study`, a data.frame with one row per
#'   design scenario and columns for planned facet counts, variance terms,
#'   projected `G`, projected `Phi`, and interpretation bands.
#'
#' @references
#' Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N.
#' (1972). *The dependability of behavioral measurements: Theory of
#' generalizability for scores and profiles*. Wiley.
#'
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [mfrm_generalizability()], [evaluate_mfrm_design()],
#'   [recommend_mfrm_design()], [plot_data()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   gt <- mfrm_generalizability(fit)
#'   mfrm_d_study(gt, data.frame(Rater = c(2, 3, 4), Criterion = 4))
#' }
#' }
#' @export
mfrm_d_study <- function(x,
                         design_grid = NULL,
                         object_facet = "Person",
                         random_facets = NULL,
                         residual_scaling = c("highest_order", "single_condition", "none", "sensitivity"),
                         ...) {
  residual_scaling <- match.arg(residual_scaling)
  if (inherits(x, "mfrm_fit")) {
    x <- mfrm_generalizability(
      x,
      object_facet = object_facet,
      random_facets = random_facets,
      ...
    )
  }
  if (!inherits(x, "mfrm_generalizability")) {
    stop("`x` must be output from mfrm_generalizability() or an mfrm_fit object.", call. = FALSE)
  }

  random_facets <- as.character(x$design$random_facets %||% character(0))
  if (length(random_facets) == 0L) {
    stop("No random measurement facets are available for D-study projection.", call. = FALSE)
  }
  observed_levels <- x$design$observed_levels %||% NULL
  if (is.null(design_grid)) {
    if (is.null(observed_levels)) {
      stop("`design_grid` is required because observed facet counts are unavailable.", call. = FALSE)
    }
    design_grid <- as.data.frame(as.list(observed_levels[random_facets]), stringsAsFactors = FALSE)
  } else if (is.list(design_grid) && !is.data.frame(design_grid)) {
    design_grid <- as.data.frame(design_grid, stringsAsFactors = FALSE)
  } else {
    design_grid <- as.data.frame(design_grid, stringsAsFactors = FALSE)
  }
  if (nrow(design_grid) == 0L) {
    stop("`design_grid` must contain at least one design row.", call. = FALSE)
  }

  count_cols <- vapply(random_facets, function(facet) {
    prefixed <- paste0("n_", facet)
    if (facet %in% names(design_grid)) {
      facet
    } else if (prefixed %in% names(design_grid)) {
      prefixed
    } else {
      NA_character_
    }
  }, character(1))
  if (anyNA(count_cols)) {
    missing <- random_facets[is.na(count_cols)]
    stop(
      "`design_grid` is missing count column(s) for random facet(s): ",
      paste(missing, collapse = ", "),
      ". Use either the facet name or `n_<facet>`.",
      call. = FALSE
    )
  }

  counts <- as.data.frame(
    lapply(count_cols, function(col) suppressWarnings(as.numeric(design_grid[[col]]))),
    stringsAsFactors = FALSE
  )
  names(counts) <- paste0("n_", random_facets)
  counts_matrix <- as.matrix(counts)
  if (any(!is.finite(counts_matrix) | counts_matrix <= 0)) {
    stop("All D-study facet counts must be positive finite numbers.", call. = FALSE)
  }
  scaling_levels <- if (identical(residual_scaling, "sensitivity")) {
    c("highest_order", "single_condition", "none")
  } else {
    residual_scaling
  }
  if (length(scaling_levels) > 1L) {
    row_index <- rep(seq_len(nrow(counts)), each = length(scaling_levels))
    counts <- counts[row_index, , drop = FALSE]
    scaling_col <- rep(scaling_levels, times = length(unique(row_index)))
  } else {
    scaling_col <- rep(scaling_levels, nrow(counts))
  }

  vc <- as.data.frame(x$variance_components, stringsAsFactors = FALSE)
  v <- stats::setNames(suppressWarnings(as.numeric(vc$Variance)), as.character(vc$Source))
  sigma2_p <- as.numeric(v[x$design$object_facet] %||% NA_real_)
  sigma2_residual <- as.numeric(v["Residual"] %||% NA_real_)
  if (!is.finite(sigma2_residual)) sigma2_residual <- 0
  if (!is.finite(sigma2_p) || sigma2_p <= 0) {
    projected_g <- rep(NA_real_, nrow(counts))
    projected_phi <- rep(NA_real_, nrow(counts))
    rel_error <- rep(NA_real_, nrow(counts))
    abs_error <- rep(NA_real_, nrow(counts))
    residual_divisor <- rep(NA_real_, nrow(counts))
  } else {
    residual_divisor <- vapply(seq_len(nrow(counts)), function(i) {
      row_counts <- unlist(counts[i, , drop = TRUE], use.names = FALSE)
      switch(
        scaling_col[i],
        highest_order = prod(row_counts),
        single_condition = min(row_counts),
        none = 1,
        prod(row_counts)
      )
    }, numeric(1))
    rel_error <- sigma2_residual / residual_divisor
    facet_main_error <- numeric(nrow(counts))
    for (facet in random_facets) {
      sigma2_f <- as.numeric(v[facet] %||% 0)
      if (!is.finite(sigma2_f)) sigma2_f <- 0
      facet_main_error <- facet_main_error + sigma2_f / counts[[paste0("n_", facet)]]
    }
    abs_error <- facet_main_error + rel_error
    projected_g <- sigma2_p / (sigma2_p + rel_error)
    projected_phi <- sigma2_p / (sigma2_p + abs_error)
  }

  classify_coef <- function(value) {
    dplyr::case_when(
      !is.finite(value) ~ "unavailable",
      value >= 0.80 ~ "high_stakes_candidate",
      value >= 0.70 ~ "routine_candidate",
      TRUE ~ "review"
    )
  }
  out <- cbind(
    data.frame(Scenario = seq_len(nrow(counts)), counts, stringsAsFactors = FALSE),
    data.frame(
      ResidualScaling = scaling_col,
      ResidualDivisor = residual_divisor,
      ObjectVariance = sigma2_p,
      RelativeErrorVariance = rel_error,
      AbsoluteErrorVariance = abs_error,
      G = round(projected_g, 4),
      Phi = round(projected_phi, 4),
      GStatus = classify_coef(projected_g),
      PhiStatus = classify_coef(projected_phi),
      stringsAsFactors = FALSE
    )
  )
  attr(out, "object_facet") <- x$design$object_facet
  attr(out, "random_facets") <- random_facets
  attr(out, "residual_scaling") <- residual_scaling
  attr(out, "source") <- "mfrm_generalizability"
  class(out) <- c("mfrm_d_study", "data.frame")
  out
}

#' @export
print.mfrm_d_study <- function(x, ...) {
  cat("mfrmr D-study projection\n")
  cat("  Object of measurement:", attr(x, "object_facet") %||% NA_character_, "\n")
  cat("  Random facets:", paste(attr(x, "random_facets") %||% character(0), collapse = ", "), "\n\n")
  cat("  Residual scaling:", attr(x, "residual_scaling") %||% paste(unique(x$ResidualScaling), collapse = ", "), "\n\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @export
plot.mfrm_d_study <- function(x,
                              y = NULL,
                              type = c("coefficients", "error_variance", "heatmap", "contour", "surface3d"),
                              x_var = NULL,
                              y_var = NULL,
                              group_var = NULL,
                              panel_by = NULL,
                              panel_grid = NULL,
                              metric = NULL,
                              draw = TRUE,
                              main = NULL,
                              palette = NULL,
                              preset = c("standard", "publication", "compact", "monochrome"),
                              ...) {
  type <- match.arg(type)
  tbl <- as.data.frame(x, stringsAsFactors = FALSE)
  n_cols <- grep("^n_", names(tbl), value = TRUE)
  if (length(n_cols) == 0L) {
    stop("D-study table does not contain planned-count columns.", call. = FALSE)
  }
  if (is.null(x_var)) {
    x_var <- n_cols[1L]
  }
  x_var <- as.character(x_var[1L])
  if (!x_var %in% names(tbl)) {
    stop("`x_var` must be one of: ", paste(n_cols, collapse = ", "), call. = FALSE)
  }

  is_surface <- type %in% c("heatmap", "contour", "surface3d")
  if (is_surface) {
    if (is.null(y_var)) {
      candidates <- setdiff(n_cols, x_var)
      if (length(candidates) == 0L) {
        stop("`y_var` is required for D-study heatmap/contour plots.", call. = FALSE)
      }
      y_var <- candidates[1L]
    }
    y_var <- as.character(y_var[1L])
    if (!y_var %in% names(tbl)) {
      stop("`y_var` must be one of: ", paste(setdiff(n_cols, x_var), collapse = ", "), call. = FALSE)
    }
    if (identical(y_var, x_var)) {
      stop("`y_var` must differ from `x_var`.", call. = FALSE)
    }
  }

  coefficient_cols <- c("G", "Phi")
  error_cols <- c("RelativeErrorVariance", "AbsoluteErrorVariance")
  available_metrics <- if (identical(type, "coefficients")) {
    coefficient_cols
  } else if (identical(type, "error_variance")) {
    error_cols
  } else {
    c(coefficient_cols, error_cols)
  }
  if (is_surface && is.null(metric)) {
    metric <- if ("Phi" %in% available_metrics) "Phi" else available_metrics[1L]
  }
  if (!is.null(metric)) {
    metric <- as.character(metric)
    unknown_metric <- setdiff(metric, c(coefficient_cols, error_cols))
    if (length(unknown_metric) > 0L) {
      stop("`metric` must be one of: ",
           paste(c(coefficient_cols, error_cols), collapse = ", "), call. = FALSE)
    }
    metric_cols <- intersect(metric, available_metrics)
    if (length(metric_cols) == 0L) {
      stop("`metric` is not compatible with `type = \"", type, "\"`.", call. = FALSE)
    }
  } else {
    metric_cols <- available_metrics
  }
  if (is_surface && length(metric_cols) != 1L) {
    stop("D-study surface plots require exactly one `metric`.", call. = FALSE)
  }
  missing_metrics <- setdiff(metric_cols, names(tbl))
  if (length(missing_metrics) > 0L) {
    stop("D-study table is missing metric column(s): ",
         paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  scaling <- if ("ResidualScaling" %in% names(tbl)) {
    as.character(tbl$ResidualScaling)
  } else {
    rep("projection", nrow(tbl))
  }
  series_tbl <- do.call(rbind, lapply(metric_cols, function(metric_name) {
    tmp <- tbl
    tmp$Metric <- metric_name
    tmp$MetricFamily <- "G-theory"
    tmp$MetricRole <- dplyr::case_when(
      metric_name == "G" ~ "relative_decision",
      metric_name == "Phi" ~ "absolute_decision",
      metric_name == "RelativeErrorVariance" ~ "relative_error",
      metric_name == "AbsoluteErrorVariance" ~ "absolute_error",
      TRUE ~ "projection"
    )
    tmp$ResidualScaling <- scaling
    tmp$X <- suppressWarnings(as.numeric(tmp[[x_var]]))
    tmp$Y <- if (is_surface) suppressWarnings(as.numeric(tmp[[y_var]])) else NA_real_
    tmp$Value <- suppressWarnings(as.numeric(tmp[[metric_name]]))
    tmp
  }))
  series_tbl <- series_tbl[is.finite(series_tbl$X) & is.finite(series_tbl$Value), , drop = FALSE]
  if (is_surface) {
    series_tbl <- series_tbl[is.finite(series_tbl$Y), , drop = FALSE]
  }

  plot_vars <- c(names(tbl), "Metric", "MetricFamily", "MetricRole", "ResidualScaling")
  validate_plot_var <- function(value, arg_name, allow_null = TRUE, max_len = Inf) {
    if (is.null(value)) {
      if (isTRUE(allow_null)) return(NULL)
      stop("`", arg_name, "` is required.", call. = FALSE)
    }
    value <- as.character(value)
    value <- value[nzchar(value)]
    if (length(value) == 0L) {
      if (isTRUE(allow_null)) return(NULL)
      stop("`", arg_name, "` is required.", call. = FALSE)
    }
    if (length(value) > max_len) {
      stop("`", arg_name, "` must have length <= ", max_len, ".", call. = FALSE)
    }
    missing <- setdiff(value, plot_vars)
    if (length(missing) > 0L) {
      stop("`", arg_name, "` must use column(s) from the D-study plot data: ",
           paste(plot_vars, collapse = ", "), ".", call. = FALSE)
    }
    value
  }
  if (is.null(group_var) && !is_surface) {
    group_candidates <- setdiff(n_cols, x_var)
    if (length(group_candidates) > 0L) {
      group_var <- group_candidates[1L]
    }
  }
  group_var <- validate_plot_var(group_var, "group_var", allow_null = TRUE, max_len = 1L)
  panel_by <- validate_plot_var(panel_by, "panel_by", allow_null = TRUE, max_len = 1L)
  panel_grid <- validate_plot_var(panel_grid, "panel_grid", allow_null = TRUE, max_len = 2L)
  if (!is.null(panel_by) && !is.null(panel_grid)) {
    stop("Use either `panel_by` or `panel_grid`, not both.", call. = FALSE)
  }
  if (is_surface && is.null(panel_by) && is.null(panel_grid) &&
      "ResidualScaling" %in% names(series_tbl) &&
      dplyr::n_distinct(series_tbl$ResidualScaling) > 1L) {
    panel_by <- "ResidualScaling"
  }

  style <- resolve_plot_preset(preset)
  if (nrow(series_tbl) == 0L) {
    stop("No finite D-study values are available for plotting.", call. = FALSE)
  }

  if (is_surface) {
    panel_vars <- c(panel_by, panel_grid)
    if (length(panel_vars) == 0L) {
      series_tbl$Panel <- "All designs"
      panel_levels <- "All designs"
    } else {
      series_tbl$Panel <- do.call(paste, c(series_tbl[, panel_vars, drop = FALSE], sep = " / "))
      panel_levels <- unique(series_tbl$Panel)
    }
    fill_values <- range(series_tbl$Value, na.rm = TRUE)
    fill_cols <- if (is.null(palette)) {
      if (identical(style$name, "monochrome")) {
        grDevices::gray.colors(18L, start = 0.95, end = 0.25)
      } else {
        grDevices::hcl.colors(18L, palette = "YlGnBu", rev = TRUE)
      }
    } else {
      rep(as.character(palette), length.out = 18L)
    }
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      panel_n <- length(panel_levels)
      graphics::par(mfrow = grDevices::n2mfrow(panel_n))
      for (panel in panel_levels) {
        s <- series_tbl[series_tbl$Panel == panel, , drop = FALSE]
        x_levels <- sort(unique(s$X))
        y_levels <- sort(unique(s$Y))
        z <- matrix(NA_real_, nrow = length(x_levels), ncol = length(y_levels))
        for (i in seq_len(nrow(s))) {
          xi <- match(s$X[i], x_levels)
          yi <- match(s$Y[i], y_levels)
          z[xi, yi] <- s$Value[i]
        }
        z_finite <- z[is.finite(z)]
        has_contours <- length(unique(z_finite)) > 1L
        if (identical(type, "heatmap")) {
          graphics::image(
            x_levels, y_levels, z,
            col = fill_cols,
            zlim = fill_values,
            xlab = x_var,
            ylab = y_var,
            main = main %||% paste(metric_cols[1L], panel, sep = " / ")
          )
          if (has_contours) {
            graphics::contour(x_levels, y_levels, z, add = TRUE, drawlabels = TRUE)
          }
        } else if (identical(type, "contour")) {
          if (has_contours) {
            graphics::contour(
              x_levels, y_levels, z,
              xlab = x_var,
              ylab = y_var,
              main = main %||% paste(metric_cols[1L], panel, sep = " / "),
              drawlabels = TRUE
            )
          } else {
            graphics::plot(
              range(x_levels, na.rm = TRUE),
              range(y_levels, na.rm = TRUE),
              type = "n",
              xlab = x_var,
              ylab = y_var,
              main = main %||% paste(metric_cols[1L], panel, sep = " / ")
            )
            graphics::text(mean(range(x_levels, na.rm = TRUE)), mean(range(y_levels, na.rm = TRUE)), "constant surface")
          }
        } else {
          zlim <- range(z_finite, na.rm = TRUE)
          if (!all(is.finite(zlim))) {
            zlim <- c(0, 1)
          }
          if (isTRUE(all.equal(zlim[1L], zlim[2L]))) {
            pad <- max(1e-6, abs(zlim[1L]) * 1e-6)
            zlim <- zlim + c(-pad, pad)
          }
          z_cols <- if (is.null(palette)) {
            if (identical(style$name, "monochrome")) "gray85" else style$fill_soft
          } else {
            as.character(palette)[1L]
          }
          graphics::persp(
            x_levels, y_levels, z,
            theta = 35,
            phi = 25,
            col = z_cols,
            border = grDevices::adjustcolor(style$foreground, alpha.f = 0.35),
            ticktype = "detailed",
            xlab = x_var,
            ylab = y_var,
            zlab = metric_cols[1L],
            zlim = zlim,
            main = main %||% paste(metric_cols[1L], panel, sep = " / ")
          )
        }
      }
    }
    return(invisible(new_mfrm_plot_data(
      "d_study",
      list(
        plot = type,
        table = tbl,
        series = series_tbl,
        surface = series_tbl,
        metric_family = "G-theory",
        metric = metric_cols[1L],
        x_var = x_var,
        y_var = y_var,
        group_var = group_var %||% NA_character_,
        panel_by = panel_by %||% NA_character_,
        panel_grid = panel_grid %||% character(0),
        title = main %||% paste("D-study", metric_cols[1L], type),
        subtitle = "Projection from simplified G-study variance components",
        legend = new_plot_legend(
          label = metric_cols[1L],
          role = "metric",
          aesthetic = switch(type, heatmap = "fill", contour = "contour", surface3d = "surface", "value"),
          value = paste(fill_values, collapse = " to ")
        ),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  panel_vars <- c(panel_by, panel_grid)
  group_components <- unique(c("Metric", "ResidualScaling", group_var))
  group_components <- setdiff(group_components[!is.na(group_components) & nzchar(group_components)], panel_vars)
  if (length(group_components) == 0L) {
    series_tbl$Series <- "Projection"
  } else {
    series_tbl$Series <- do.call(paste, c(series_tbl[, group_components, drop = FALSE], sep = " / "))
  }
  if (length(panel_grid) == 2L) {
    series_tbl$PanelRow <- as.character(series_tbl[[panel_grid[1L]]])
    series_tbl$PanelCol <- as.character(series_tbl[[panel_grid[2L]]])
    series_tbl$Panel <- paste(series_tbl$PanelRow, series_tbl$PanelCol, sep = " / ")
  } else if (!is.null(panel_by)) {
    series_tbl$Panel <- as.character(series_tbl[[panel_by]])
    series_tbl$PanelRow <- series_tbl$Panel
    series_tbl$PanelCol <- "panel"
  } else {
    series_tbl$Panel <- "All designs"
    series_tbl$PanelRow <- "All designs"
    series_tbl$PanelCol <- "panel"
  }
  series_levels <- unique(series_tbl$Series)
  if (is.null(palette)) {
    series_cols <- if (identical(style$name, "monochrome")) {
      stats::setNames(rep(style$foreground, length(series_levels)), series_levels)
    } else {
      stats::setNames(
        grDevices::hcl.colors(max(3L, length(series_levels)), palette = "Dark 3")[seq_along(series_levels)],
        series_levels
      )
    }
  } else {
    palette <- as.character(palette)
    series_cols <- stats::setNames(rep(palette, length.out = length(series_levels)), series_levels)
  }
  line_types <- stats::setNames(rep(seq_len(6L), length.out = length(series_levels)), series_levels)

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    if (length(panel_grid) == 2L) {
      row_levels <- unique(series_tbl$PanelRow)
      col_levels <- unique(series_tbl$PanelCol)
      graphics::par(mfrow = c(length(row_levels), length(col_levels)))
      panel_specs <- expand.grid(PanelRow = row_levels, PanelCol = col_levels, stringsAsFactors = FALSE)
    } else {
      panel_levels <- unique(series_tbl$Panel)
      graphics::par(mfrow = grDevices::n2mfrow(length(panel_levels)))
      panel_specs <- data.frame(Panel = panel_levels, stringsAsFactors = FALSE)
    }
    y_lim <- if (identical(type, "coefficients")) {
      c(0, 1)
    } else {
      range(c(0, series_tbl$Value), na.rm = TRUE)
    }
    for (i in seq_len(nrow(panel_specs))) {
      if (length(panel_grid) == 2L) {
        s_panel <- series_tbl[
          series_tbl$PanelRow == panel_specs$PanelRow[i] &
            series_tbl$PanelCol == panel_specs$PanelCol[i],
          ,
          drop = FALSE
        ]
        panel_title <- paste(panel_specs$PanelRow[i], panel_specs$PanelCol[i], sep = " / ")
      } else {
        s_panel <- series_tbl[series_tbl$Panel == panel_specs$Panel[i], , drop = FALSE]
        panel_title <- panel_specs$Panel[i]
      }
      graphics::plot(
        s_panel$X,
        s_panel$Value,
        type = "n",
        xlab = x_var,
        ylab = if (identical(type, "coefficients")) "Coefficient" else "Error variance",
        ylim = y_lim,
        main = main %||% panel_title
      )
      graphics::grid(col = style$grid)
      for (series in unique(s_panel$Series)) {
        s <- s_panel[s_panel$Series == series, , drop = FALSE]
        s <- s[order(s$X), , drop = FALSE]
        graphics::lines(s$X, s$Value, col = series_cols[series], lty = line_types[series], lwd = 2)
        graphics::points(s$X, s$Value, col = series_cols[series], pch = 16)
      }
      if (identical(type, "coefficients")) {
        graphics::abline(h = c(0.70, 0.80), col = grDevices::adjustcolor(style$neutral, alpha.f = 0.6), lty = c(3, 2))
      }
      graphics::legend(
        "bottomright",
        legend = unique(s_panel$Series),
        col = unname(series_cols[unique(s_panel$Series)]),
        lty = unname(line_types[unique(s_panel$Series)]),
        pch = 16,
        bty = "n",
        cex = 0.72
      )
    }
  }

  invisible(new_mfrm_plot_data(
    "d_study",
    list(
      plot = type,
      table = tbl,
      series = series_tbl,
      metric_family = "G-theory",
      metric = metric_cols,
      x_var = x_var,
      y_var = y_var %||% NA_character_,
      group_var = group_var %||% NA_character_,
      panel_by = panel_by %||% NA_character_,
      panel_grid = panel_grid %||% character(0),
      title = main %||% if (identical(type, "coefficients")) "D-study G/Phi projection" else "D-study error variance projection",
      subtitle = "Projection from simplified G-study variance components",
      legend = new_plot_legend(
        label = series_levels,
        role = rep("series", length(series_levels)),
        aesthetic = rep("line", length(series_levels)),
        value = unname(series_cols[series_levels])
      ),
      reference_lines = if (identical(type, "coefficients")) {
        new_reference_lines(
          axis = rep("y", 2L),
          value = c(0.70, 0.80),
          label = c("routine", "high_stakes"),
          linetype = c("dotted", "dashed"),
          role = rep("decision_band", 2L)
        )
      } else {
        new_reference_lines()
      },
      preset = style$name
    )
  ))
}

#' @export
print.mfrm_generalizability <- function(x, ...) {
  cat("Generalizability-theory decomposition\n")
  cat(sprintf("  Object of measurement: %s\n",
              x$design$object_facet))
  cat(sprintf("  Random facets: %s\n",
              paste(x$design$random_facets, collapse = ", ")))
  cat("\nVariance components\n")
  print(x$variance_components, row.names = FALSE)
  cat(sprintf("\nG (relative): %.3f | Phi (absolute): %.3f\n",
              as.numeric(x$coefficients$G),
              as.numeric(x$coefficients$Phi)))
  if (length(x$design$lmer_warnings) > 0L) {
    cat(sprintf("\n%d lme4 warning(s) suppressed; results may be unstable.\n",
                length(x$design$lmer_warnings)))
  }
  invisible(x)
}
