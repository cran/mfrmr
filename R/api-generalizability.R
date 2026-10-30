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

mfrmr_gt_boundary_status <- function(fit_lmer,
                                     lmer_warnings = character(0),
                                     lmer_messages = character(0),
                                     tolerance = 1e-4) {
  singular_fit <- tryCatch(
    isTRUE(lme4::isSingular(fit_lmer, tol = tolerance)),
    error = function(e) NA
  )
  boundary_text <- c(lmer_warnings, lmer_messages)
  boundary_message <- any(grepl(
    "boundary|singular",
    boundary_text,
    ignore.case = TRUE
  ))
  opt <- fit_lmer@optinfo
  convergence_review <- length(lmer_warnings) > 0L ||
    any(unlist(opt$conv$opt) != 0, na.rm = TRUE) ||
    length(opt$conv$lme4$messages) > 0L
  status <- if (isTRUE(singular_fit) || isTRUE(boundary_message)) {
    "boundary_or_singular_fit"
  } else if (is.na(singular_fit)) {
    "singularity_check_unavailable"
  } else if (convergence_review) {
    "convergence_review"
  } else {
    "identified"
  }
  note <- switch(
    status,
    boundary_or_singular_fit = paste(
      "The lme4 random-effects fit is singular or on a boundary;",
      "zero variance estimates are retained. This is distinct from whether",
      "the design separates variance components; inspect component_identification.",
      "Point coefficients do not establish uncertainty or decision accuracy."
    ),
    singularity_check_unavailable = paste(
      "The lme4 singularity check was unavailable;",
      "treat G/D-study coefficients as caveated until the fitted",
      "random-effects structure is reviewed."
    ),
    convergence_review = "The mixed-model fit has unresolved convergence warnings; review its numerical stability before using G/D projections.",
    "No lme4 boundary or singular fit was detected."
  )
  list(
    identification_status = status,
    identification_note = note,
    boundary_fit = identical(status, "boundary_or_singular_fit"),
    singular_fit = singular_fit,
    singular_tolerance = tolerance,
    boundary_message = boundary_message
  )
}

mfrmr_gt_classify_coef <- function(value,
                                   identification_status = "identified") {
  if (!identical(as.character(identification_status)[1L], "identified")) {
    return(ifelse(
      is.finite(value),
      "identification_warning",
      "unavailable"
    ))
  }
  dplyr::case_when(
    !is.finite(value) ~ "unavailable",
    value >= 0.80 ~ "at_or_above_0.80_reference",
    value >= 0.70 ~ "at_or_above_0.70_reference",
    TRUE ~ "below_0.70_reference"
  )
}

validate_gtheory_output <- function(x) {
  version <- if (inherits(x, "mfrm_generalizability")) x$design$calculation_version else attr(x, "calculation_version")
  if (!identical(version, 2L)) {
    stop(paste("This saved G/D result used earlier variance handling.",
      "Recreate mfrm_generalizability(fit) and any mfrm_d_study() results;",
      "the MFRM fit need not be re-estimated."), call. = FALSE)
  }
  invisible(TRUE)
}

gtheory_scaling_label <- function(x) {
  labels <- c(highest_order = "Divide residual by all facet counts",
              single_condition = "Divide residual by smallest facet count",
              none = "Keep residual unchanged", sensitivity = "Compare three residual assumptions")
  ifelse(x %in% names(labels), unname(labels[x]), x)
}

print_gtheory_data_usage <- function(usage) {
  if (is.null(usage)) {
    cat("  Input-row accounting is unavailable for this saved result.\n")
  } else {
    counts <- usage$counts
    cat("  G-study rows:", counts[["InputRows"]], "input;",
        counts[["UsedRows"]], "used;", counts[["ExcludedRows"]], "excluded.\n")
    if (identical(usage$source, "Stored fitted rows")) {
      cat("  Counts start from stored fitted rows; earlier MFRM filtering is not included.\n")
    }
    if (counts[["ExcludedRows"]] > 0L) {
      cat("  Incomplete rows were explicitly omitted; no missing values were imputed.\n")
    }
  }
}

#' Main-effects generalizability approximation for an MFRM design
#'
#' Re-fits the rating data underlying an `mfrm_fit` as a crossed
#' random-effects model
#' `Score ~ 1 + (1 | Person) + (1 | Facet1) + ... + Residual`
#' via `lme4::lmer`, and returns main-effect variance components,
#' residual variance, and G / Phi coefficients under this simplified model.
#' Interactions are not estimated separately. For a G-study including the
#' supported interactions, use [mfrm_g_study()] with one or two
#' common random facets; it also accepts a single score column.
#'
#' The decomposition is on the observed numeric `Score` scale. It does not use
#' the fitted MFRM latent scale or estimate a latent ordinal G/Phi coefficient.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param data Optional data frame. When `NULL`, the rating data
#'   stored on `fit$prep$data` is used. Required columns are the selected facets
#'   and `Score`. Scores must be numeric or numeric character/factor labels;
#'   nonnumeric labels and infinite values are refused. Use `NA` for missing
#'   values. Facet labels must be nonblank; `Score` and `Residual` are reserved
#'   and cannot be facet names.
#' @param object_facet Facet that plays the role of the "object of
#'   measurement" -- typically `"Person"` (default).
#' @param random_facets Character vector of other facets to
#'   treat as random conditions of measurement. The default uses `Person`
#'   and the configured MFRM facets, excluding `object_facet`. For example,
#'   when raters are the objects, persons remain a random condition.
#' @param reml Logical, passed to [lme4::lmer()] (default `TRUE`).
#' @param missing Either `"error"` (default) or `"omit"`. Missing scores or
#'   selected facet values stop the analysis by default. Explicit omission
#'   fits only complete rows and records the excluded row positions and missing
#'   columns. Missingness in unselected columns does not exclude a row.
#'
#' @return An object of class `mfrm_generalizability` with:
#' \describe{
#'   \item{`variance_components`}{One row per random effect plus
#'     residual, with columns `Source`, `Variance`, `ProportionVariance`,
#'     and `Identified`. When `Identified` is false, the fitted allocation
#'     and its proportion are not uniquely determined by the source design.}
#'   \item{`coefficients`}{One-row data frame with `G`
#'     (generalizability coefficient, relative decision) and
#'     `Phi` (dependability coefficient, absolute decision), coefficient
#'     status labels, and the identification status of the fitted
#'     random-effects model. `GReason` and `PhiReason` explain missing
#'     coefficients; an empty reason means the point is calculable.
#'     `ComponentIdentificationChecked` records the design check.}
#'   \item{`design`}{Description of the crossed-random model.
#'     `component_identification` retains the covariance-kernel Gram matrix,
#'     rank, scaled eigenvalues, null space and numerical tolerance.
#'     `fit_status` retains the separate lme4 boundary/convergence assessment.}
#'   \item{`data_usage`}{Input source, omission policy, named `counts`
#'     (`InputRows`, `UsedRows`, `ExcludedRows`), `excluded_rows`, and
#'     `missing_cells` (`InputRow`, `Column`). Row positions refer to the supplied
#'     data, or stored fitted rows when `data = NULL`. They cannot recover rows
#'     previously removed during MFRM fitting. Counts also accompany the
#'     coefficient table and D-study projections, including tabular exports;
#'     `GStudyDataSource` identifies the scope of those counts.}
#' }
#'
#' @section Interpretation:
#' - `G` is appropriate for **relative** decisions (rank-ordering
#'   persons): `G = sigma2(p) / (sigma2(p) + sigma2(Residual))`.
#' - The reported `Phi` describes dependability for **absolute** decisions:
#'   `Phi = sigma2(p) / (sigma2(p) + sigma2(facet
#'   main effects) + sigma2(Residual))`, before D-study scaling.
#'   It does not estimate the probability of correct classification at a
#'   particular cut score.
#' - Use [mfrm_d_study()] to project `G` / `Phi` under planned numbers of
#'   raters, items, criteria, or other random measurement facets.
#' - Values of 0.70 and 0.80 are displayed as familiar planning references,
#'   not as universal decision rules. Required dependability depends on the
#'   decision, consequences, population, and evidence beyond a single
#'   coefficient.
#' - For ordered categories, `Score` is treated as a numeric observed response
#'   in a Gaussian linear mixed model; thresholding is not modeled.
#'
#' @section Component identification and boundary estimates:
#' A zero variance estimate can be a valid constrained ML/REML solution.
#' It is distinct from a design that cannot separate two variance components.
#' For example, if every rater identifier always corresponds to the same
#' task identifier, their main-effect covariance kernels coincide. Their
#' sum may be identified while the individual allocations are not.
#' `isSingular()` alone does not test this design property.
#'
#' The design check uses linear independence of the observation covariance
#' kernels (after projecting out the intercept for REML). Only variance
#' combinations orthogonal to the numerical null space enter coefficients.
#' Thus an identified sum can still supply Phi, while a future plan assigning
#' unequal weights to its confounded components cannot. Individual fitted
#' allocations remain available with `Identified = FALSE` for diagnosis.
#' A full numerical rank does not establish correct model assumptions,
#' optimizer accuracy or finite-sample precision.
#' This is a linear design check; fitted boundary zeros are not used to
#' resolve otherwise confounded components.
#'
#' If object variance is zero and the relevant error variance is positive,
#' G or Phi equals zero and is retained with boundary diagnostics. A zero
#' total variance gives an undefined ratio and `zero_total_variance`.
#' The historical `IdentificationStatus` field summarizes design and fit
#' cautions; inspect `design$component_identification` and `design$fit_status`
#' to distinguish them. This function does not supply coefficient intervals.
#'
#' @section Limitations:
#' This helper formulates the random-effects model with main effects
#' only (`Score ~ 1 + (1|Person) + (1|Facet1) + ... + Residual`); no
#' explicit `(1 | Person:Rater)`, `(1 | Person:Criterion)`, or
#' `(1 | Rater:Criterion)` interaction terms are estimated. All
#' interactions are omitted. The residual combines unexplained variation;
#' omitted interactions may also affect the fitted main-effect components.
#' Their separate variances and correct averaging rates are not recovered. This function reports the
#' one-observation-per-cell baseline. [mfrm_d_study()] applies D-study
#' scaling, including residual-scaling sensitivity checks, to the same
#' simplified variance-component decomposition.
#' Because person-by-facet interaction terms are not estimated separately,
#' D-study projections remain practical planning evidence rather than a
#' replacement for a fully specified G-theory design. Counts held constant in
#' a D-study do not turn a random facet into a fixed-facet universe. Variance
#' estimation uncertainty is not propagated to G/Phi. Component values and G/Phi are
#' stored at full precision; rounding is only for display. Recreate older
#' G/D results from the existing MFRM fit before reuse.
#' Boundary or singular `lme4` fits are retained as diagnostic evidence but are
#' not treated as decision-ready G/D-study evidence.
#' New results record component separation explicitly. Earlier saved G-studies
#' without this record can still supply conditional D-study calculations,
#' with `ComponentIdentificationChecked = FALSE`; recreate the G-study from
#' the original ratings to obtain the check. This needs no MFRM refit.
#'
#' Omission does not correct missing-data bias or identify why ratings are
#' absent. Entirely absent assignments are not reconstructed. Review the
#' rating design and missingness assumptions before interpreting G/D results.
#' Earlier versions silently omitted incomplete rows and unparseable scores.
#' To reproduce complete-row selection, clean invalid labels explicitly and
#' choose `missing = "omit"`. Older saved results remain usable when their
#' calculation version is current, but unavailable row counts are not guessed;
#' rerun the G-study with the original data to obtain row accounting.
#'
#' @references
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software, 67*(1),
#' 1--48. \doi{10.18637/jss.v067.i01}. Gaussian ML/REML model and optimization.
#' See also \url{https://lme4.github.io/lme4/reference/isSingular.html} for
#' the meaning of a singular fitted covariance matrix.
#'
#' Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N.
#'   (1972). *The dependability of behavioral measurements: Theory
#'   of generalizability for scores and profiles*. Wiley.
#'
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   \doi{10.1007/978-1-4757-3456-0}. General G/D framework; the main-effects
#'   restriction and residual-scaling choices here are package-specific limits.
#'
#' @seealso [mfrm_g_study()], [mfrm_d_study()],
#'   [compute_facet_icc()], [diagnose_mfrm()]
#' @examples
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   gt <- mfrm_generalizability(fit)
#'   gt$variance_components
#'   # Look for: a Person variance component well above any single
#'   #   non-person facet's variance share. Large rater or criterion
#'   #   variance shares mean those conditions add measurement error
#'   #   relative to person spread.
#'   gt$coefficients
#'   # Compare G and Phi with study-specific requirements; 0.70 and 0.80
#'   #   are reference guides only. Phi < G means absolute decisions are noisier than relative
#'   #   decisions; review whether facet main effects need anchoring.
#'   # Always check IdentificationStatus before using the bands:
#'   gt$coefficients[, c("G", "Phi", "GStatus", "PhiStatus",
#'                       "IdentificationStatus")]
#'   gt$design$identification_note
#'   # If IdentificationStatus is not "identified", treat G/Phi as
#'   # design-review evidence rather than decision-ready reliability.
#' }
#' @export
mfrm_generalizability <- function(fit,
                                  data = NULL,
                                  object_facet = "Person",
                                  random_facets = NULL,
                                  reml = TRUE,
                                  missing = c("error", "omit")) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("`mfrm_generalizability()` requires the `lme4` package ",
         "(in Suggests). Install it and retry.", call. = FALSE)
  }
  missing <- match.arg(missing)
  if (!is.logical(reml) || length(reml) != 1L || is.na(reml)) {
    stop("`reml` must be TRUE or FALSE.", call. = FALSE)
  }
  data_source <- if (is.null(data)) "Stored fitted rows" else "Supplied data"
  if (is.null(data)) {
    data <- as.data.frame(fit$prep$data %||% data.frame(),
                          stringsAsFactors = FALSE)
  }
  if (!is.data.frame(data) || nrow(data) == 0L || anyNA(names(data)) ||
      anyDuplicated(names(data)) || any(!nzchar(names(data)))) {
    stop("`data` must be a nonempty data frame with unique, nonmissing column names.", call. = FALSE)
  }
  data <- as.data.frame(data)

  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (is.null(random_facets)) {
    random_facets <- setdiff(c("Person", facet_names), object_facet)
  }
  random_facets <- as.character(random_facets)
  if (!is.character(object_facet) || length(object_facet) != 1L || is.na(object_facet) ||
      anyNA(random_facets) || anyDuplicated(c(object_facet, random_facets)) ||
      any(!nzchar(c(object_facet, random_facets)))) {
    stop("The object facet and random facets must be distinct, non-missing names.", call. = FALSE)
  }
  if (length(random_facets) == 0L) {
    stop("At least one other facet is required as a random ",
         "condition of measurement.", call. = FALSE)
  }
  if (any(c(object_facet, random_facets) %in% c("Score", "Residual"))) {
    stop("`Score` and `Residual` are reserved and cannot be facet names.", call. = FALSE)
  }
  needed_cols <- c(object_facet, random_facets, "Score")
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Data frame is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }

  for (col in c(object_facet, random_facets)) {
    value <- data[[col]]
    if (!is.null(dim(value)) || !(is.factor(value) || (!is.object(value) &&
        (is.character(value) || is.logical(value) || (is.numeric(value) && !is.complex(value))))) ||
        (is.numeric(value) && any(is.infinite(value)))) {
      stop("Facet '", col, "' must contain finite numeric, character, factor, or logical labels, or NA.", call. = FALSE)
    }
    if (any(!is.na(value) & !nzchar(trimws(as.character(value))))) {
      stop("Facet '", col, "' contains blank labels; use NA for missing values.", call. = FALSE)
    }
    if (!is.factor(data[[col]])) {
      labels <- as.character(value)
      labels[is.na(value)] <- NA_character_
      data[[col]] <- as.factor(labels)
    }
  }
  score <- data$Score
  if (!is.null(dim(score)) || !(is.factor(score) || (!is.object(score) &&
      (is.character(score) || (is.numeric(score) && !is.complex(score)))))) {
    stop("`Score` must be numeric or numeric character/factor labels, with NA for missing values.", call. = FALSE)
  }
  data$Score <- suppressWarnings(as.numeric(if (is.factor(score)) as.character(score) else score))
  if (any(!is.na(score) & !is.finite(data$Score))) {
    stop("`Score` contains nonnumeric or infinite values; correct them explicitly and use NA for missing values.", call. = FALSE)
  }
  absent <- is.na(data[, needed_cols, drop = FALSE])
  keep <- rowSums(absent) == 0L
  cells <- which(absent, arr.ind = TRUE)
  data_usage <- list(source = data_source, missing = missing,
    counts = c(InputRows = nrow(data), UsedRows = sum(keep), ExcludedRows = sum(!keep)),
    excluded_rows = which(!keep),
    missing_cells = data.frame(InputRow = cells[, 1L], Column = needed_cols[cells[, 2L]],
                               row.names = NULL))
  if (any(!keep) && missing == "error") {
    stop(sum(!keep), " row(s) have missing Score or selected facet values. Review the data or choose missing = 'omit' explicitly.", call. = FALSE)
  }
  if (!any(keep)) stop("No complete rows remain for the G-study.", call. = FALSE)
  data <- data[keep, , drop = FALSE]

  random_terms <- c(object_facet, random_facets)
  formula_str <- paste0(
    "Score ~ 1 + ",
    paste0("(1 | ", vapply(random_terms, function(name) deparse1(as.name(name), backtick = TRUE),
                          character(1)), ")", collapse = " + ")
  )
  formula <- stats::as.formula(formula_str)

  lmer_warnings <- character(0)
  lmer_messages <- character(0)
  fit_lmer <- tryCatch(
    withCallingHandlers(
      lme4::lmer(formula, data = data, REML = reml, na.action = stats::na.fail),
      warning = function(w) {
        lmer_warnings <<- c(lmer_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        lmer_messages <<- c(lmer_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit_lmer, "error")) {
    stop("lme4::lmer failed: ", conditionMessage(fit_lmer), call. = FALSE)
  }
  boundary_status <- mfrmr_gt_boundary_status(
    fit_lmer,
    lmer_warnings = lmer_warnings,
    lmer_messages = lmer_messages
  )

  vc <- as.data.frame(lme4::VarCorr(fit_lmer))
  vc <- vc[is.na(vc$var2), c("grp", "vcov")]
  variance_scale <- max(vc$vcov)
  scaled_variances <- if (is.finite(variance_scale) && variance_scale > 0)
    vc$vcov / variance_scale else rep(NA_real_, nrow(vc))
  var_components <- data.frame(
    Source = as.character(vc$grp),
    Variance = as.numeric(vc$vcov),
    ProportionVariance = scaled_variances / sum(scaled_variances),
    stringsAsFactors = FALSE
  )

  # Independent covariance kernels, not optimizer convergence, identify
  # components. Some projections remain identified when components alias.
  geometry <- .mfrmr_gt_geometry(data, random_terms, reml)
  var_components$Identified <- vapply(var_components$Source, function(source) {
    w <- stats::setNames(rep(0, length(geometry$scale)), names(geometry$scale))
    w[source] <- 1
    .mfrmr_gt_estimable(w, geometry)
  }, logical(1))
  fit_status <- boundary_status$identification_status
  if (geometry$rank < geometry$components) {
    boundary_status$identification_status <- "components_not_identified"
    boundary_status$identification_note <- paste(
      "The observed assignments do not separate all variance components.",
      "Retained component allocations are fitted candidates; only identified",
      "error sums and coefficients are reported. Review component_identification.",
      boundary_status$identification_note)
  }
  v <- stats::setNames(var_components$Variance, var_components$Source)
  baseline <- .mfrmr_gt_projection(v, object_facet, random_facets,
    as.data.frame(as.list(stats::setNames(rep(1, length(random_facets)), random_facets))),
    "highest_order", geometry)

  out <- list(
    variance_components = var_components,
    coefficients = data.frame(
      G = baseline$G,
      Phi = baseline$Phi,
      GReason = baseline$GReason,
      PhiReason = baseline$PhiReason,
      GStatus = mfrmr_gt_classify_coef(
        baseline$G,
        boundary_status$identification_status
      ),
      PhiStatus = mfrmr_gt_classify_coef(
        baseline$Phi,
        boundary_status$identification_status
      ),
      IdentificationStatus = boundary_status$identification_status,
      ComponentIdentificationChecked = TRUE,
      InputRows = data_usage$counts[["InputRows"]],
      UsedRows = data_usage$counts[["UsedRows"]],
      ExcludedRows = data_usage$counts[["ExcludedRows"]],
      GStudyDataSource = data_usage$source,
      stringsAsFactors = FALSE
    ),
    data_usage = data_usage,
    design = list(
      calculation_version = 2L,
      object_facet = object_facet,
      random_facets = random_facets,
      observed_levels = stats::setNames(
        vapply(c(object_facet, random_facets), function(col) {
          dplyr::n_distinct(stats::model.frame(fit_lmer)[[col]])
        }, integer(1)),
        c(object_facet, random_facets)
      ),
      formula = format(formula),
      estimand_scale = "observed_numeric_score",
      reml = isTRUE(reml),
      lmer_warnings = lmer_warnings,
      lmer_messages = lmer_messages,
      identification_status = boundary_status$identification_status,
      identification_note = boundary_status$identification_note,
      boundary_fit = boundary_status$boundary_fit,
      singular_fit = boundary_status$singular_fit,
      singular_tolerance = boundary_status$singular_tolerance,
      boundary_message = boundary_status$boundary_message,
      fit_status = fit_status,
      component_identification = geometry
    )
  )
  class(out) <- c("mfrm_generalizability", "list")
  out
}

#' Main-effects D-study approximation
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
#'   `n_Rater`). Supply exactly one count column per facet, with unique names
#'   and an unambiguous one-to-one correspondence. Extra columns, duplicate
#'   aliases for a facet, or a column serving two facets are refused.
#'   Counts must be positive integers. When `NULL`, one row using
#'   the observed number of levels is
#'   returned.
#' @param object_facet,random_facets Passed to [mfrm_generalizability()] when
#'   `x` is an `mfrm_fit`.
#' @param residual_scaling How the collapsed residual variance should be scaled
#'   when planned facet counts increase. `"highest_order"` treats the residual
#'   as highest-order person-by-all-conditions/error variance and divides by
#'   the product of planned counts. `"single_condition"` divides by the smallest
#'   planned facet count, a sensitivity assumption for variation associated
#'   with the least-replicated condition. It is not a confidence bound. `"none"` leaves the residual
#'   unscaled. `"sensitivity"` returns all three assumptions for each design
#'   row.
#' @param ... Additional arguments passed to [mfrm_generalizability()] when `x`
#'   is an `mfrm_fit`. With an already estimated main-effects G-study or its
#'   D-study summary, unused arguments are rejected. In particular, composite
#'   `weights` are unsupported here; they require the score-component route
#'   from [mfrm_g_study()].
#'
#' @details
#' The projection uses the variance decomposition already estimated by
#' [mfrm_generalizability()]. For a random measurement facet `j`, main-effect
#' variance contributes `sigma2_j / n_j` to the absolute-error denominator.
#' Under the fitted additive model, residual errors are independent across
#' observations. For a complete future crossing with common facet levels for
#' all persons, their mean has variance `sigma2_Residual / prod(n_j)`.
#' Facet-wide shifts cancel from relative comparisons but enter absolute error
#' as `sum(sigma2_j / n_j)`. This justifies `"highest_order"` within that model.
#' If omitted interactions actually exist, the fitted residual does not
#' recover their separate variances or correlations; dividing it by the product
#' is then an assumption. `"single_condition"` and `"none"` are alternative
#' planning assumptions, not estimates of these correlations, confidence
#' bounds, or a guarantee that the true projection lies between them.
#' The relative-decision
#' denominator uses only this scaled residual term. Missing or invalid required
#' components leave the affected coefficient unavailable; they are not zero
#' variance. These are point projections conditional on estimated components,
#' not uncertainty bounds or automatic recommendations for sample size.
#' Returned G/Phi retain full precision, as do the values used for their
#' reference-status labels. Round only for presentation, after any comparisons.
#'
#' Error variances are evaluated independently of object variance: they
#' remain available when object variance is zero or missing. Zero object
#' variance with positive error gives a coefficient of zero; zero total
#' variance is undefined. For checked sources, only identified variance sums
#' enter the output. For example, two confounded facet main effects can still
#' contribute their identified sum when their planned counts are equal;
#' unequal counts require the separate, unidentified components.
#' Coefficients use scaled sums to avoid overflow. A coefficient can therefore
#' remain calculable even when its error variance exceeds floating-point range.
#' Such unrepresentable quantities are returned as `NA` with reasons.
#'
#' This is a pragmatic D-study planning layer, not a full p x r x i ANOVA
#' decomposition. If person-by-rater or person-by-item interactions are a
#' primary estimand, consider [mfrm_g_study()] and
#' [mfrm_d_study()] for one or two common random facets,
#' including Person-by-Task, Person-by-Rater, and Person-by-Rater-by-Task. Those
#' functions also accept a single score and estimate the corresponding
#' interaction components explicitly. Changing `residual_scaling` here only
#' explores assumptions; it does not estimate the omitted interactions.
#'
#' Pass the result directly to [mfrm_report()] for full projection tables and
#' their assumptions. Newly created results retain source variance components,
#' estimation settings and observed facet counts. Reports of older compatible
#' results mark missing source details `not_recorded`; report creation never
#' reconstructs components or computes confidence intervals.
#'
#' The `G` and `Phi` values returned here belong to the generalizability-theory
#' metric family. They should not be interpreted as coefficient alpha, omega,
#' KR-20, or IRT marginal/separation reliability, even though all of those
#' summaries may be displayed on a 0--1 scale in broader reporting dashboards.
#' They remain on the observed numeric score scale inherited from
#' [mfrm_generalizability()], not the fitted MFRM latent scale.
#'
#' @return An object of class `mfrm_d_study`, a data.frame with one row per
#'   design scenario and columns for planned facet counts, variance terms,
#'   projected `G`, projected `Phi`, interpretation bands, and identification
#'   status inherited from [mfrm_generalizability()].
#'   `InputRows`, `UsedRows`, and `ExcludedRows` describe the source G-study,
#'   not the planned D-study sample; `GStudyDataSource` identifies their scope.
#'   The `data_usage` attribute retains its
#'   row accounting, including after subsetting. Older source results without
#'   accounting have `NA` counts; they are not assumed to have used all rows.
#'   `summary(d)` returns the coefficient data frame with these saved attributes,
#'   rather than generic per-column descriptive statistics. No intervals or
#'   new variance-component estimates are computed.
#'   `GReason`, `PhiReason`, `ObjectVarianceReason`, `RelativeErrorReason`,
#'   `AbsoluteErrorReason` and `ResidualDivisorReason` explain unavailable
#'   quantities; empty strings mean no calculation failure. Reasons distinguish
#'   unidentified targets, missing/invalid components, zero total variance and
#'   numerical range. Fit cautions remain separate in `IdentificationStatus`
#'   and `IdentificationNote`. `ComponentIdentificationChecked` and the
#'   `component_identification` attribute retain the source design check.
#'
#' @references
#' Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N.
#' (1972). *The dependability of behavioral measurements: Theory of
#' generalizability for scores and profiles*. Wiley.
#'
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#' \doi{10.1007/978-1-4757-3456-0}. This framework does not recover the
#' interaction components omitted by the source main-effects model.
#'
#' @seealso [plot.mfrm_d_study()], [mfrm_generalizability()],
#'   [mfrm_d_study()], [evaluate_mfrm_design()],
#'   [recommend_mfrm_design()], [plot_data()]
#' @examples
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   gt <- mfrm_generalizability(fit)
#'   ds <- mfrm_d_study(gt, data.frame(Rater = c(2, 3, 4), Criterion = 4))
#'   ds[, c("n_Rater", "n_Criterion", "G", "Phi",
#'          "GStatus", "PhiStatus", "IdentificationStatus")]
#'   # If IdentificationStatus is not "identified", even large G/Phi
#'   # values remain identification warnings, not decision-ready evidence.
#' }
#' @name mfrm_d_study_main_effects
#' @method mfrm_d_study mfrm_generalizability
#' @export
mfrm_d_study.mfrm_generalizability <- function(x,
                         design_grid = NULL,
                         residual_scaling = c("highest_order", "single_condition", "none", "sensitivity"),
                         ...) {
  rlang::check_dots_empty()
  residual_scaling <- match.arg(residual_scaling)
  if (!inherits(x, "mfrm_generalizability")) {
    stop("`x` must be output from mfrm_generalizability() or an mfrm_fit object.", call. = FALSE)
  }
  validate_gtheory_output(x)

  random_facets <- as.character(x$design$random_facets %||% character(0))
  if (length(random_facets) == 0L) {
    stop("No random measurement facets are available for D-study projection.", call. = FALSE)
  }
  observed_levels <- x$design$observed_levels %||% NULL
  if (is.null(design_grid)) {
    if (is.null(observed_levels)) {
      stop("`design_grid` is required because observed facet counts are unavailable.", call. = FALSE)
    }
    design_grid <- as.list(observed_levels[random_facets])
  }
  design_grid <- mfrm_d_study_count_columns(design_grid, random_facets,
    paste0("n_", random_facets))
  count_cols <- names(design_grid)

  counts <- as.data.frame(
    lapply(count_cols, function(col) suppressWarnings(as.numeric(as.character(design_grid[[col]])))),
    stringsAsFactors = FALSE
  )
  names(counts) <- paste0("n_", random_facets)
  counts_matrix <- as.matrix(counts)
  if (any(!is.finite(counts_matrix) | counts_matrix <= 0 | counts_matrix != floor(counts_matrix))) {
    stop("All D-study facet counts must be positive finite integers.", call. = FALSE)
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
  if (anyDuplicated(vc$Source)) stop("Variance-component sources must be unique.", call. = FALSE)
  v <- stats::setNames(suppressWarnings(as.numeric(vc$Variance)), as.character(vc$Source))
  v[!is.finite(v) | v < 0] <- NA_real_
  projection <- .mfrmr_gt_projection(v, x$design$object_facet, random_facets,
    counts, scaling_col, x$design$component_identification)

  identification_status <- as.character(x$design$identification_status %||% "review_required")[1L]
  if (is.na(identification_status) || !nzchar(identification_status)) {
    identification_status <- "review_required"
  }
  identification_note <- as.character(x$design$identification_note %||% "")[1L]
  boundary_fit <- isTRUE(x$design$boundary_fit)
  classify_coef <- function(value) {
    mfrmr_gt_classify_coef(value, identification_status)
  }
  out <- cbind(
    data.frame(Scenario = seq_len(nrow(counts)), counts, stringsAsFactors = FALSE, check.names = FALSE),
    data.frame(
      ResidualScaling = scaling_col,
      projection,
      GStatus = classify_coef(projection$G),
      PhiStatus = classify_coef(projection$Phi),
      IdentificationStatus = rep(identification_status, nrow(counts)),
      BoundaryFit = rep(boundary_fit, nrow(counts)),
      IdentificationNote = rep(identification_note, nrow(counts)),
      ComponentIdentificationChecked = !is.null(x$design$component_identification),
      InputRows = x$data_usage$counts[["InputRows"]] %||% NA_integer_,
      UsedRows = x$data_usage$counts[["UsedRows"]] %||% NA_integer_,
      ExcludedRows = x$data_usage$counts[["ExcludedRows"]] %||% NA_integer_,
      GStudyDataSource = x$data_usage$source %||% NA_character_,
      stringsAsFactors = FALSE
    )
  )
  attr(out, "object_facet") <- x$design$object_facet
  attr(out, "calculation_version") <- 2L
  attr(out, "random_facets") <- random_facets
  attr(out, "residual_scaling") <- residual_scaling
  attr(out, "source") <- "mfrm_generalizability"
  attr(out, "source_design") <- x$design
  attr(out, "source_components") <- x$variance_components
  attr(out, "estimand_scale") <- x$design$estimand_scale %||%
    "observed_numeric_score"
  attr(out, "identification_status") <- identification_status
  attr(out, "identification_note") <- identification_note
  attr(out, "boundary_fit") <- boundary_fit
  attr(out, "data_usage") <- x$data_usage
  attr(out, "component_identification") <- x$design$component_identification
  class(out) <- c("mfrm_d_study", "data.frame")
  out
}

#' @export
`[.mfrm_d_study` <- function(x, ...) {
  out <- NextMethod("[")
  if (is.data.frame(out)) {
    metadata <- setdiff(names(attributes(x)), c("names", "row.names", "class"))
    for (name in metadata) attr(out, name) <- attr(x, name, exact = TRUE)
  }
  out
}

#' @export
print.mfrm_d_study <- function(x, ...) {
  validate_gtheory_output(x)
  cat("mfrmr D-study projection\n")
  cat("  Object of measurement:", attr(x, "object_facet") %||% NA_character_, "\n")
  cat("  Random facets:", paste(attr(x, "random_facets") %||% character(0), collapse = ", "), "\n\n")
  cat("  Estimand scale: observed numeric score\n")
  print_gtheory_data_usage(attr(x, "data_usage", exact = TRUE))
  cat("  Residual assumption:", gtheory_scaling_label(attr(x, "residual_scaling") %||% unique(x$ResidualScaling)), "\n\n")
  if (!identical(attr(x, "identification_status") %||% "identified", "identified")) {
    cat("  Model fit requires review.\n")
    note <- attr(x, "identification_note") %||% ""
    if (nzchar(note)) {
      cat("  Note:", note, "\n")
    }
    cat("\n")
  }
  shown <- as.data.frame(x)
  if ("ResidualScaling" %in% names(shown)) {
    shown$ResidualScaling <- gtheory_scaling_label(shown$ResidualScaling)
  }
  shown <- shown[, setdiff(names(shown), c("GStatus", "PhiStatus", "IdentificationStatus", "BoundaryFit", "IdentificationNote",
                                         "InputRows", "UsedRows", "ExcludedRows", "GStudyDataSource",
                                         "ComponentIdentificationChecked", grep("Reason$", names(shown), value = TRUE))), drop = FALSE]
  print.data.frame(shown, row.names = FALSE, ...)
  if (any(vapply(x[grep("Reason$", names(x), value = TRUE)], function(z) any(nzchar(z)), logical(1))))
    print_wrapped_line("Some quantities are unavailable or exceed numerical range; inspect the returned ...Reason columns.")
  if (!isTRUE(all(x$ComponentIdentificationChecked)))
    print_wrapped_line("Component separation was not checked for this saved source. Recreate the G-study from the original ratings to check it.")
  print_wrapped_line("Observed-score planning projections hold estimated variance components fixed. They do not establish cut-score accuracy or an adequate rating design.")
  cat("\n  Note: 0.70 and 0.80 are reference guides, not universal decision rules.\n")
  invisible(x)
}

#' Plot design comparisons from a main-effects D-study
#'
#' Compare planned facet counts using an existing [mfrm_d_study()] result.
#' G concerns relative ordering; Phi also includes shifts in absolute score
#' levels. Larger coefficients mean greater dependability under the selected
#' model and residual assumption, not proven pass/fail accuracy.
#'
#' @param x An [mfrm_d_study()] result.
#' @param y Reserved for method compatibility.
#' @param type `"coefficients"` for G/Phi curves, `"error_variance"` for error
#'   variance curves, or `"heatmap"`, `"contour"`, or `"surface3d"` for one
#'   metric over two facet counts. Error variance is in squared score units,
#'   not SEM. Lower error variance is better.
#' @param x_var,y_var Planned-count columns such as `"n_Rater"`. The default
#'   horizontal axis is the first count column; surface plots use the next
#'   column on the other axis. The two axes must differ.
#' @param group_var Optional additional column distinguishing curves. All
#'   non-horizontal facet counts and residual assumptions remain separate
#'   within each panel, including when `group_var` is supplied.
#' @param panel_by One column defining panels, or `NULL`.
#' @param panel_grid One or two columns defining panels. Use this or
#'   `panel_by`, not both. Surface plots require every other facet count and
#'   residual assumption to be constant within each panel. Subset the result
#'   or add panels when they vary; they cannot be silently averaged or overlaid.
#' @param metric Optional selection from `"G"`, `"Phi"`,
#'   `"RelativeErrorVariance"`, or `"AbsoluteErrorVariance"`, compatible with
#'   `type`. Surface plots require one metric and default to `"Phi"`.
#' @param draw Draw when `TRUE`; `FALSE` only returns plot data.
#' @param main Optional plot title.
#' @param palette Optional colors.
#' @param preset Plot style: `"standard"`, `"publication"`, `"compact"`, or
#'   `"monochrome"`.
#' @param reference Optional distinct finite numeric reference values for curves,
#'   between 0 and 1 for coefficients or nonnegative for error variance.
#'   `NULL` draws no guides. Not supported for surface plots.
#' @param ... Reserved for method compatibility.
#' @details Points represent requested scenarios; connecting lines and
#'   contours are visual guides. Missing estimates remain missing and break
#'   curves. If none are available, inspect the D-study table and source
#'   variance components. No reference lines are drawn by default. Specify
#'   `reference = c(0.70, 0.80)` to reproduce the earlier coefficient guides.
#'   User-chosen guides are not universal acceptance criteria or confidence
#'   bounds; their values and interpretation remain in the returned plot data.
#'   Heatmaps include a numeric color key; exact values remain in the table.
#'
#'   The main-effects G-study combines unmodeled interactions in its residual.
#'   Different residual-scaling curves describe assumptions, not confidence
#'   bounds. All projections hold estimated components fixed. Check source
#'   fit warnings before interpreting even large coefficients. For supported
#'   designs with separately estimated interactions, including a single score,
#'   use [mfrm_g_study()] and [mfrm_d_study()].
#' @return Invisibly, an `mfrm_plot_data` object with the scenario `table`,
#'   metric `series`, axis/group/panel settings, and labels. Use [plot_data()]
#'   for custom graphics. Automatic [as_ggplot()] conversion is not provided
#'   for this class; the base plots preserve the chosen comparisons.
#' @seealso [mfrm_d_study()], [plot.mfrm_multivariate_d_study()]
#' @examples
#' # After creating ds with mfrm_d_study():
#' # plot(ds, x_var = "n_Rater", panel_grid = c("Metric", "ResidualScaling"))
#' # For Rater x Task x Occasion scenarios, separate occasions explicitly:
#' # plot(ds, type = "heatmap", x_var = "n_Rater", y_var = "n_Task",
#' #      metric = "Phi", panel_by = "n_Occasion")
#' # With residual_scaling = "sensitivity", use
#' # panel_grid = c("n_Occasion", "ResidualScaling") instead.
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @md
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
                              reference = NULL,
                              ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  validate_gtheory_output(x)
  type <- match.arg(type)
  projection_note <- if (!identical(attr(x, "identification_status"), "identified")) {
    "Model fit requires review; conditional observed-score projections only."
  } else {
    "Observed-score point projections; variance estimation uncertainty omitted."
  }
  display_labels <- function(z) vapply(strsplit(z, " / ", fixed = TRUE), function(parts) {
    paste(gtheory_scaling_label(parts), collapse = " / ")
  }, character(1))
  tbl <- as.data.frame(x, stringsAsFactors = FALSE)
  n_cols <- grep("^n_", names(tbl), value = TRUE)
  column_labels <- function(nm, data) {
    if (nm %in% n_cols) paste(sub("^n_", "", nm), "=", data[[nm]]) else as.character(data[[nm]])
  }
  axis_label <- function(nm) paste(sub("^n_", "", nm), "count")
  if (length(n_cols) == 0L) {
    stop("D-study table does not contain planned-count columns.", call. = FALSE)
  }
  if (is.null(x_var)) {
    x_var <- n_cols[1L]
  }
  x_var <- as.character(x_var[1L])
  if (!x_var %in% n_cols) {
    stop("`x_var` must be one of: ", paste(n_cols, collapse = ", "), call. = FALSE)
  }

  is_surface <- type %in% c("heatmap", "contour", "surface3d")
  if (is_surface && !is.null(reference)) {
    stop("`reference` is supported for curves, not surface plots.", call. = FALSE)
  }
  references <- .mfrm_d_study_references(reference, type == "coefficients")
  if (is_surface) {
    if (is.null(y_var)) {
      candidates <- setdiff(n_cols, x_var)
      if (length(candidates) == 0L) {
        stop("`y_var` is required for D-study heatmap/contour plots.", call. = FALSE)
      }
      y_var <- candidates[1L]
    }
    y_var <- as.character(y_var[1L])
    if (!y_var %in% n_cols) {
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
  series_tbl <- series_tbl[is.finite(series_tbl$X), , drop = FALSE]
  series_tbl$Value[!is.finite(series_tbl$Value)] <- NA_real_
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
  if (length(panel_grid) == 1L) {
    panel_by <- panel_grid
    panel_grid <- NULL
  }
  if (is_surface && is.null(panel_by) && is.null(panel_grid) &&
      "ResidualScaling" %in% names(series_tbl) &&
      dplyr::n_distinct(series_tbl$ResidualScaling) > 1L) {
    panel_by <- "ResidualScaling"
  }

  style <- resolve_plot_preset(preset)
  if (nrow(series_tbl) == 0L || !any(is.finite(series_tbl$Value))) {
    stop("No finite D-study values are available for plotting. Inspect the D-study table and source variance components.", call. = FALSE)
  }

  if (is_surface) {
    panel_vars <- c(panel_by, panel_grid)
    conditions <- setdiff(c(n_cols, "ResidualScaling"), c(x_var, y_var, panel_vars))
    unfixed <- conditions[vapply(conditions, function(nm) {
      values <- unique(series_tbl[c(panel_vars, nm)])
      if (!length(panel_vars)) nrow(values) > 1L else anyDuplicated(values[panel_vars]) > 0L
    }, logical(1))]
    if (length(unfixed)) {
      stop("Surface plots must hold remaining conditions fixed within each panel. Varying: ",
        paste(unfixed, collapse = ", "), ". Subset to one value or use `panel_by` / `panel_grid`.", call. = FALSE)
    }
    if (length(panel_vars) == 0L) {
      series_tbl$Panel <- "All designs"
      panel_levels <- "All designs"
    } else {
      series_tbl$Panel <- do.call(paste, c(lapply(panel_vars, column_labels, data = series_tbl), sep = " / "))
      panel_levels <- unique(series_tbl$Panel)
    }
    series_tbl$Panel <- display_labels(series_tbl$Panel)
    panel_levels <- unique(series_tbl$Panel)
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
    if (fill_values[1L] == fill_values[2L]) fill_cols[] <- fill_cols[9L]
    key_values <- unique(seq(fill_values[1L], fill_values[2L], length.out = 5L))
    key_colors <- fill_cols[if (length(key_values) == 1L) 9L else round(seq(1, length(fill_cols), length.out = 5L))]
    heat_key <- new_plot_legend(format(signif(key_values, 3), trim = TRUE),
      rep("metric_value", length(key_values)), rep("fill", length(key_values)), key_colors)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      old_par <- graphics::par()[c("mfrow", "cex", "mex", "mar", "oma")]
      on.exit(graphics::par(old_par), add = TRUE)
      panel_n <- length(panel_levels)
      graphics::par(mfrow = grDevices::n2mfrow(panel_n), oma = c(1.8, 0, 0, 0),
        mar = c(4.1, 4.1, 3.1, if (type == "heatmap") 7 else 2.1))
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
        if (!length(z_finite)) {
          plot(range(x_levels), range(y_levels), type = "n",
            xlab = axis_label(x_var), ylab = axis_label(y_var), xaxt = "n", yaxt = "n",
            main = main %||% paste(metric_cols[1L], panel, sep = " / "))
          axis(1, at = x_levels)
          axis(2, at = y_levels, las = 1)
          text(mean(range(x_levels)), mean(range(y_levels)),
            "Estimates unavailable\nInspect the D-study table.", cex = 0.8)
          next
        }
        has_contours <- length(unique(z_finite)) > 1L
        if (identical(type, "heatmap")) {
          image(
            x_levels, y_levels, z,
            col = fill_cols,
            zlim = fill_values,
            xlab = axis_label(x_var),
            ylab = axis_label(y_var), xaxt = "n", yaxt = "n",
            main = main %||% paste(metric_cols[1L], panel, sep = " / ")
          )
          axis(1, at = x_levels)
          axis(2, at = y_levels, las = 1)
          usr <- graphics::par("usr")
          legend(usr[2] + 0.03 * diff(usr[1:2]), usr[4],
            legend = heat_key$label, fill = heat_key$value, title = metric_cols[1L],
            bty = "n", xpd = NA, cex = 0.75)
        } else if (identical(type, "contour")) {
          if (has_contours) {
            graphics::contour(
              x_levels, y_levels, z,
              xlab = axis_label(x_var),
              ylab = axis_label(y_var),
              main = main %||% paste(metric_cols[1L], panel, sep = " / "),
              drawlabels = TRUE
            )
          } else {
            plot(
              range(x_levels, na.rm = TRUE),
              range(y_levels, na.rm = TRUE),
              type = "n",
              xlab = axis_label(x_var),
              ylab = axis_label(y_var),
              main = main %||% paste(metric_cols[1L], panel, sep = " / ")
            )
            text(mean(range(x_levels, na.rm = TRUE)), mean(range(y_levels, na.rm = TRUE)), "constant surface")
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
            xlab = axis_label(x_var),
            ylab = axis_label(y_var),
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
        subtitle = projection_note,
        legend = if (type == "heatmap") heat_key else new_plot_legend(
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
  group_components <- unique(c("Metric", "ResidualScaling", group_var, setdiff(n_cols, x_var)))
  group_components <- setdiff(group_components[!is.na(group_components) & nzchar(group_components)], panel_vars)
  if (length(group_components) == 0L) {
    series_tbl$Series <- "Projection"
  } else {
    series_labels <- lapply(group_components, column_labels, data = series_tbl)
    series_tbl$Series <- do.call(paste, c(series_labels, sep = " / "))
  }
  if (length(panel_grid) == 2L) {
    series_tbl$PanelRow <- column_labels(panel_grid[1L], series_tbl)
    series_tbl$PanelCol <- column_labels(panel_grid[2L], series_tbl)
    series_tbl$Panel <- paste(series_tbl$PanelRow, series_tbl$PanelCol, sep = " / ")
  } else if (!is.null(panel_by)) {
    series_tbl$Panel <- column_labels(panel_by, series_tbl)
    series_tbl$PanelRow <- series_tbl$Panel
    series_tbl$PanelCol <- "panel"
  } else {
    series_tbl$Panel <- "All designs"
    series_tbl$PanelRow <- "All designs"
    series_tbl$PanelCol <- "panel"
  }
  for (col in c("Series", "Panel", "PanelRow", "PanelCol")) {
    series_tbl[[col]] <- display_labels(series_tbl[[col]])
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
    old_par <- graphics::par()[c("mfrow", "cex", "mex", "oma")]
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(oma = c(1.8, 0, 0, 0))
    if (length(panel_grid) == 2L) {
      row_levels <- unique(series_tbl$PanelRow)
      col_levels <- unique(series_tbl$PanelCol)
      graphics::par(mfrow = c(length(row_levels), length(col_levels)))
      panel_specs <- expand.grid(PanelCol = col_levels, PanelRow = row_levels, stringsAsFactors = FALSE)
    } else {
      panel_levels <- unique(series_tbl$Panel)
      graphics::par(mfrow = grDevices::n2mfrow(length(panel_levels)))
      panel_specs <- data.frame(Panel = panel_levels, stringsAsFactors = FALSE)
    }
    y_lim <- if (identical(type, "coefficients")) {
      c(0, 1)
    } else {
      range(c(0, series_tbl$Value, references$value), na.rm = TRUE)
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
      plot(
        s_panel$X,
        s_panel$Value,
        type = "n",
        xlab = axis_label(x_var), xaxt = "n",
        ylab = if (identical(type, "coefficients")) "Coefficient" else "Error variance",
        ylim = y_lim,
        main = main %||% panel_title
      )
      axis(1, at = sort(unique(s_panel$X)))
      graphics::grid(col = style$grid)
      for (series in unique(s_panel$Series)) {
        s <- s_panel[s_panel$Series == series, , drop = FALSE]
        s <- s[order(s$X), , drop = FALSE]
        lines(s$X, s$Value, col = series_cols[series], lty = line_types[series], lwd = 2)
        points(s$X, s$Value, col = series_cols[series], pch = 16)
      }
      if (nrow(references)) abline(h = references$value,
        col = style$neutral, lty = references$linetype)
      legend(
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
      subtitle = projection_note,
      legend = new_plot_legend(
        label = series_levels,
        role = rep("series", length(series_levels)),
        aesthetic = rep("line", length(series_levels)),
        value = unname(series_cols[series_levels])
      ),
      reference_lines = references,
      interpretation_note = .mfrm_d_study_reference_note(references),
      preset = style$name
    )
  ))
}

#' @export
print.mfrm_generalizability <- function(x, ...) {
  validate_gtheory_output(x)
  cat("Generalizability-theory decomposition\n")
  cat(sprintf("  Object of measurement: %s\n",
              x$design$object_facet))
  cat(sprintf("  Random facets: %s\n",
              paste(x$design$random_facets, collapse = ", ")))
  cat("  Estimand scale: observed numeric score\n")
  print_gtheory_data_usage(x$data_usage)
  print_wrapped_line("This main-effects model does not separate person-by-facet interactions. G/Phi are point summaries; uncertainty in estimated variance components is omitted.")
  cat("\nVariance components\n")
  print(x$variance_components, row.names = FALSE, digits = 4)
  if (!is.null(x$design$component_identification))
    cat("  Component kernel rank:", x$design$component_identification$rank, "of",
      x$design$component_identification$components, "\n")
  cat(sprintf("\nG (relative): %.3f | Phi (absolute): %.3f\n",
              as.numeric(x$coefficients$G),
              as.numeric(x$coefficients$Phi)))
  for (metric in c("G", "Phi")) {
    reason <- x$coefficients[[paste0(metric, "Reason")]]
    if (length(reason) && nzchar(reason)) cat(" ", metric, "unavailable:", reason, "\n")
  }
  if (!identical(x$design$identification_status %||% "identified", "identified")) {
    cat("\nModel fit requires review.\n")
    cat(x$design$identification_note, "\n")
  }
  if (length(x$design$lmer_warnings) > 0L) {
    cat(sprintf("\n%d lme4 warning(s) suppressed; results may be unstable.\n",
                length(x$design$lmer_warnings)))
  }
  if (length(x$design$lmer_messages) > 0L) {
    cat(sprintf("\n%d lme4 message(s) suppressed; review boundary/convergence diagnostics.\n",
                length(x$design$lmer_messages)))
  }
  invisible(x)
}
