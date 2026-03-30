#' Build a specification summary report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param title Optional analysis title.
#' @param data_file Optional data-file label (for reporting only).
#' @param output_file Optional output-file label (for reporting only).
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_specifications` (`type = "facet_elements"`,
#' `"anchor_constraints"`, `"convergence"`).
#'
#' @section Interpreting output:
#' - `header` / `data_spec`: run identity and model settings.
#' - `facet_labels`: facet sizes and labels.
#' - `convergence_control`: optimizer configuration and status.
#'
#' @section Typical workflow:
#' 1. Generate `specifications_report(fit)`.
#' 2. Verify model settings and convergence metadata.
#' 3. Use the output as methods and run-documentation support in reports.
#' @return A named list with specification-report components. Class:
#'   `mfrm_specifications`.
#' @seealso [fit_mfrm()], [data_quality_report()], [estimation_iteration_report()],
#'   [mfrmr_reports_and_tables], [mfrmr_compatibility_layer]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- specifications_report(fit, title = "Toy run")
#' summary(out)
#' p_spec <- plot(out, draw = FALSE)
#' class(p_spec)
#' @export
specifications_report <- function(fit,
                                  title = NULL,
                                  data_file = NULL,
                                  output_file = NULL,
                                  include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table1_specifications(
      fit = fit,
      title = title,
      data_file = data_file,
      output_file = output_file,
      include_fixed = include_fixed
    )
  )
  as_mfrm_bundle(out, "mfrm_specifications")
}

#' Build a data quality summary report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param data Optional raw data frame used for row-level audit.
#' @param person Optional person column name in `data`.
#' @param facets Optional facet column names in `data`.
#' @param score Optional score column name in `data`.
#' @param weight Optional weight column name in `data`.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_data_quality` (`type = "row_audit"`, `"category_counts"`,
#' `"missing_rows"`).
#'
#' @section Interpreting output:
#' - `summary`: retained/dropped row overview.
#' - `row_audit`: reason-level breakdown for data issues.
#' - `category_counts`: post-filter category usage.
#' - `unknown_elements`: facet levels in raw data but not in fitted design.
#'
#' @section Typical workflow:
#' 1. Run `data_quality_report(...)` with raw data.
#' 2. Check row-audit and missing/unknown element sections.
#' 3. Resolve issues before final estimation/reporting.
#' @return A named list with data-quality report components. Class:
#'   `mfrm_data_quality`.
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [specifications_report()],
#'   [mfrmr_reports_and_tables], [mfrmr_compatibility_layer]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- data_quality_report(
#'   fit, data = toy, person = "Person",
#'   facets = c("Rater", "Criterion"), score = "Score"
#' )
#' summary(out)
#' p_dq <- plot(out, draw = FALSE)
#' class(p_dq)
#' @export
data_quality_report <- function(fit,
                                data = NULL,
                                person = NULL,
                                facets = NULL,
                                score = NULL,
                                weight = NULL,
                                include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table2_data_summary(
      fit = fit,
      data = data,
      person = person,
      facets = facets,
      score = score,
      weight = weight,
      include_fixed = include_fixed
    )
  )
  as_mfrm_bundle(out, "mfrm_data_quality")
}

#' Build an estimation-iteration report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param max_iter Maximum replay iterations (excluding optional initial row).
#' @param reltol Stopping tolerance for replayed max-logit change.
#' @param include_prox If `TRUE`, include an initial pseudo-row labeled `PROX`.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_iteration_report` (`type = "residual"`, `"logit_change"`,
#' `"objective"`).
#'
#' @section Interpreting output:
#' - `iterations`: trajectory of convergence indicators by iteration.
#' - `summary`: final status and stopping diagnostics.
#' - optional `PROX` row: pseudo-initial reference point when enabled.
#'
#' @section Typical workflow:
#' 1. Run `estimation_iteration_report(fit)`.
#' 2. Inspect plateau/stability patterns in summary/plot.
#' 3. Adjust optimization settings if convergence looks weak.
#' @return A named list with iteration-report components. Class:
#'   `mfrm_iteration_report`.
#' @seealso [fit_mfrm()], [specifications_report()], [data_quality_report()],
#'   [mfrmr_reports_and_tables], [mfrmr_compatibility_layer]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- estimation_iteration_report(fit, max_iter = 5)
#' summary(out)
#' p_iter <- plot(out, draw = FALSE)
#' class(p_iter)
#' @export
estimation_iteration_report <- function(fit,
                                        max_iter = 20,
                                        reltol = NULL,
                                        include_prox = TRUE,
                                        include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table3_iteration_report(
      fit = fit,
      max_iter = max_iter,
      reltol = reltol,
      include_prox = include_prox,
      include_fixed = include_fixed
    )
  )
  as_mfrm_bundle(out, "mfrm_iteration_report")
}

#' Build a subset connectivity report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param top_n_subsets Optional maximum number of subset rows to keep.
#' @param min_observations Minimum observations required to keep a subset row.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_subset_connectivity` (`type = "subset_observations"`,
#' `"facet_levels"`, or `"linking_matrix"` / `"coverage_matrix"` /
#' `"design_matrix"`).
#'
#' @section Interpreting output:
#' - `summary`: number and size of connected subsets.
#' - subset table: whether data are fragmented into disconnected components.
#' - facet-level columns: where connectivity bottlenecks occur.
#'
#' @section Typical workflow:
#' 1. Run `subset_connectivity_report(fit)`.
#' 2. Confirm near-single-subset structure when possible.
#' 3. Use results to justify linking/anchoring strategy.
#' @return A named list with subset-connectivity components. Class:
#'   `mfrm_subset_connectivity`.
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [data_quality_report()],
#'   [mfrmr_linking_and_dff], [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- subset_connectivity_report(fit)
#' summary(out)
#' p_sub <- plot(out, draw = FALSE)
#' p_design <- plot(out, type = "design_matrix", draw = FALSE)
#' class(p_sub)
#' class(p_design)
#' out$summary[, c("Subset", "Observations", "ObservationPercent")]
#' @export
subset_connectivity_report <- function(fit,
                                       diagnostics = NULL,
                                       top_n_subsets = NULL,
                                       min_observations = 0) {
  out <- with_legacy_name_warning_suppressed(
    table6_subsets_listing(
      fit = fit,
      diagnostics = diagnostics,
      top_n_subsets = top_n_subsets,
      min_observations = min_observations
    )
  )
  as_mfrm_bundle(out, "mfrm_subset_connectivity")
}

#' Build a facet statistics report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param metrics Numeric columns in `diagnostics$measures` to summarize.
#' @param ruler_width Width of the fixed-width ruler used for `M/S/Q/X` marks.
#' @param distribution_basis Which distribution basis to keep in the appended
#'   precision summary: `"both"` (default), `"sample"`, or `"population"`.
#' @param se_mode Which standard-error mode to keep in the appended precision
#'   summary: `"both"` (default), `"model"`, or `"fit_adjusted"`.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_facet_statistics` (`type = "means"`, `"sds"`, `"ranges"`).
#'
#' @section Interpreting output:
#' - facet-level means/SD/ranges of selected metrics (`Estimate`, fit indices, `SE`).
#' - fixed-width ruler rows (`M/S/Q/X`) for compact profile scanning.
#'
#' @section Typical workflow:
#' 1. Run `facet_statistics_report(fit)`.
#' 2. Inspect summary/ranges for anomalous facets.
#' 3. Cross-check flagged facets with fit and chi-square diagnostics.
#' The returned bundle now includes:
#' - `precision_summary`: facet precision/separation indices by
#'   `DistributionBasis` and `SEMode`
#' - `variability_tests`: fixed/random variability tests by facet
#' - `se_modes`: compact list of available SE modes by facet
#'
#' @return A named list with facet-statistics components. Class:
#'   `mfrm_facet_statistics`.
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()], [plot_facets_chisq()],
#'   [mfrmr_reports_and_tables]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- facet_statistics_report(fit)
#' summary(out)
#' p_fs <- plot(out, draw = FALSE)
#' class(p_fs)
#' @export
facet_statistics_report <- function(fit,
                                    diagnostics = NULL,
                                    metrics = c("Estimate", "Infit", "Outfit", "SE"),
                                    ruler_width = 41,
                                    distribution_basis = c("both", "sample", "population"),
                                    se_mode = c("both", "model", "fit_adjusted")) {
  distribution_basis <- match.arg(distribution_basis)
  se_mode <- match.arg(se_mode)
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  out <- with_legacy_name_warning_suppressed(
    table6_2_facet_statistics(
      fit = fit,
      diagnostics = diagnostics,
      metrics = metrics,
      ruler_width = ruler_width
    )
  )
  precision_tbl <- as.data.frame(
    diagnostics$facet_precision %||% build_facet_precision_summary(diagnostics$measures, diagnostics$facets_chisq),
    stringsAsFactors = FALSE
  )
  if (nrow(precision_tbl) > 0) {
    if (!identical(distribution_basis, "both")) {
      precision_tbl <- precision_tbl[precision_tbl$DistributionBasis == distribution_basis, , drop = FALSE]
    }
    if (!identical(se_mode, "both")) {
      precision_tbl <- precision_tbl[precision_tbl$SEMode == se_mode, , drop = FALSE]
    }
  }

  variability_tbl <- as.data.frame(diagnostics$facets_chisq %||% data.frame(), stringsAsFactors = FALSE)
  se_modes_tbl <- if (nrow(precision_tbl) == 0) {
    data.frame()
  } else {
    precision_tbl |>
      dplyr::group_by(.data$Facet, .data$SEMode, .data$SEColumn) |>
      dplyr::summarize(
        DistributionBases = paste(sort(unique(.data$DistributionBasis)), collapse = ", "),
        MeanSE = mean(.data$MeanSE, na.rm = TRUE),
        MedianSE = mean(.data$MedianSE, na.rm = TRUE),
        AvailableLevels = max(.data$SEAvailable, na.rm = TRUE),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  out$precision_summary <- precision_tbl
  out$variability_tests <- variability_tbl
  out$se_modes <- se_modes_tbl
  out$settings$distribution_basis <- distribution_basis
  out$settings$se_mode <- se_mode
  as_mfrm_bundle(out, "mfrm_facet_statistics")
}

#' Build a precision audit report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#'
#' @details
#' This helper summarizes how `mfrmr` derived SE, CI, and reliability values
#' for the current run. It is package-native and is intended to help users
#' distinguish model-based precision paths from exploratory ones without
#' requiring external software conventions.
#'
#' @section What this audit means:
#' `precision_audit_report()` is a reporting gatekeeper for precision claims.
#' It tells you how the package derived uncertainty summaries for the current
#' run and how cautiously those summaries should be written up.
#'
#' @section What this audit does not justify:
#' - It does not, by itself, validate the measurement model or substantive
#'   conclusions.
#' - A favorable precision tier does not override convergence, fit, linking,
#'   or design problems elsewhere in the analysis.
#'
#' @section Interpreting output:
#' - `profile`: one-row overview of the active precision tier and recommended use.
#' - `checks`: package-native audit checks for SE ordering, reliability ordering,
#'   coverage of sample/population summaries, and SE source labels.
#' - `approximation_notes`: method notes copied from `diagnose_mfrm()`.
#'
#' @section Recommended next step:
#' Use the `profile$PrecisionTier` and `checks` table to decide whether SE, CI,
#' and reliability language can be phrased as model-based, should be qualified
#' as hybrid, or should remain exploratory in the final report.
#'
#' @section Typical workflow:
#' 1. Run `diagnose_mfrm()` for the fitted model.
#' 2. Build `precision_audit_report(fit, diagnostics = diag)`.
#' 3. Use `summary()` to see whether the run supports model-based reporting
#'    language or should remain in exploratory/screening mode.
#'
#' @return A named list with:
#' - `profile`: one-row precision overview
#' - `checks`: package-native precision audit checks
#' - `approximation_notes`: detailed method notes
#' - `settings`: resolved model and method labels
#'
#' @seealso [diagnose_mfrm()], [facet_statistics_report()], [reporting_checklist()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' out <- precision_audit_report(fit, diagnostics = diag)
#' summary(out)
#' @export
precision_audit_report <- function(fit, diagnostics = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  profile_tbl <- as.data.frame(diagnostics$precision_profile %||% data.frame(), stringsAsFactors = FALSE)
  checks_tbl <- as.data.frame(diagnostics$precision_audit %||% data.frame(), stringsAsFactors = FALSE)
  notes_tbl <- as.data.frame(diagnostics$approximation_notes %||% data.frame(), stringsAsFactors = FALSE)
  settings <- list(
    model = as.character(fit$summary$Model[1] %||% fit$config$model %||% NA_character_),
    method = as.character(fit$summary$Method[1] %||% fit$config$method %||% NA_character_),
    precision_tier = as.character(profile_tbl$PrecisionTier[1] %||% NA_character_)
  )

  out <- list(
    profile = profile_tbl,
    checks = checks_tbl,
    approximation_notes = notes_tbl,
    settings = settings
  )
  as_mfrm_bundle(out, "mfrm_precision_audit")
}

#' Build a category structure report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param theta_range Theta/logit range used to derive transition points.
#' @param theta_points Number of grid points used for transition-point search.
#' @param drop_unused If `TRUE`, remove zero-count categories from outputs.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @param fixed_max_rows Maximum rows per fixed-width section.
#'
#' @details
#' Preferred high-level API for category-structure diagnostics.
#' This wraps the legacy-compatible bar/transition export and returns a stable
#' bundle interface for reporting and plotting.
#'
#' @section Interpreting output:
#' Key components include:
#' - category usage/fit table (count, expected, infit/outfit, ZSTD)
#' - threshold ordering and adjacent threshold gaps
#' - category transition-point table on the requested theta grid
#'
#' Practical read order:
#' 1. `summary(out)` for compact warnings and threshold ordering.
#' 2. `out$category_table` for sparse/misfitting categories.
#' 3. `out$transition_points` to inspect crossing structure.
#' 4. `plot(out)` for quick visual check.
#'
#' @section Typical workflow:
#' 1. [fit_mfrm()] -> model.
#' 2. [diagnose_mfrm()] -> residual/fit diagnostics (optional argument here).
#' 3. `category_structure_report()` -> category health snapshot.
#' 4. `summary()` and `plot()` for draft-oriented review of category structure.
#' @return A named list with category-structure components. Class:
#'   `mfrm_category_structure`.
#' @seealso [rating_scale_table()], [category_curves_report()], [plot.mfrm_fit()],
#'   [mfrmr_reports_and_tables], [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- category_structure_report(fit)
#' summary(out)
#' names(out)
#' p_cs <- plot(out, draw = FALSE)
#' class(p_cs)
#' @export
category_structure_report <- function(fit,
                                      diagnostics = NULL,
                                      theta_range = c(-6, 6),
                                      theta_points = 241,
                                      drop_unused = FALSE,
                                      include_fixed = FALSE,
                                      fixed_max_rows = 200) {
  out <- with_legacy_name_warning_suppressed(
    table8_barchart_export(
      fit = fit,
      diagnostics = diagnostics,
      theta_range = theta_range,
      theta_points = theta_points,
      drop_unused = drop_unused,
      include_fixed = include_fixed,
      fixed_max_rows = fixed_max_rows
    )
  )
  as_mfrm_bundle(out, "mfrm_category_structure")
}

#' Build a category curve export bundle (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param theta_range Theta/logit range for curve coordinates.
#' @param theta_points Number of points on the theta grid.
#' @param digits Rounding digits for numeric graph output.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @param fixed_max_rows Maximum rows shown in fixed-width graph tables.
#'
#' @details
#' Preferred high-level API for category-probability curve exports.
#' Returns tidy curve coordinates and summary metadata for quick
#' plotting/report integration without calling low-level helpers directly.
#'
#' @section Interpreting output:
#' Use this report to inspect:
#' - where each category has highest probability across theta
#' - whether adjacent categories cross in expected order
#' - whether probability bands look compressed (often sparse categories)
#'
#' Recommended read order:
#' 1. `summary(out)` for compact diagnostics.
#' 2. `out$curve_points` (or equivalent curve table) for downstream graphics.
#' 3. `plot(out)` for a default visual check.
#'
#' @section Typical workflow:
#' 1. Fit model with [fit_mfrm()].
#' 2. Run `category_curves_report()` with suitable `theta_points`.
#' 3. Use `summary()` and `plot()`; export tables for manuscripts/dashboard use.
#' @return A named list with category-curve components. Class:
#'   `mfrm_category_curves`.
#' @seealso [category_structure_report()], [rating_scale_table()], [plot.mfrm_fit()],
#'   [mfrmr_reports_and_tables], [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- category_curves_report(fit, theta_points = 101)
#' summary(out)
#' names(out)
#' p_cc <- plot(out, draw = FALSE)
#' class(p_cc)
#' @export
category_curves_report <- function(fit,
                                   theta_range = c(-6, 6),
                                   theta_points = 241,
                                   digits = 4,
                                   include_fixed = FALSE,
                                   fixed_max_rows = 400) {
  out <- with_legacy_name_warning_suppressed(
    table8_curves_export(
      fit = fit,
      theta_range = theta_range,
      theta_points = theta_points,
      digits = digits,
      include_fixed = include_fixed,
      fixed_max_rows = fixed_max_rows
    )
  )
  as_mfrm_bundle(out, "mfrm_category_curves")
}

#' Build a bias-interaction plot-data bundle (preferred alias)
#'
#' @param x Output from [estimate_bias()] or [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used when `x` is fit).
#' @param facet_a First facet name (required when `x` is fit and
#'   `interaction_facets` is not supplied).
#' @param facet_b Second facet name (required when `x` is fit and
#'   `interaction_facets` is not supplied).
#' @param interaction_facets Character vector of two or more facets.
#' @param max_abs Bound for absolute bias size when estimating from fit.
#' @param omit_extreme Omit extreme-only elements when estimating from fit.
#' @param max_iter Iteration cap for bias estimation when `x` is fit.
#' @param tol Convergence tolerance for bias estimation when `x` is fit.
#' @param top_n Maximum number of ranked rows to keep.
#' @param abs_t_warn Warning cutoff for absolute t statistics.
#' @param abs_bias_warn Warning cutoff for absolute bias size.
#' @param p_max Warning cutoff for p-values.
#' @param sort_by Ranking key: `"abs_t"`, `"abs_bias"`, or `"prob"`.
#'
#' @details
#' Preferred bundle API for interaction-bias diagnostics. The function can:
#' - use a precomputed bias object from [estimate_bias()], or
#' - estimate internally from `mfrm_fit` + facet specification.
#'
#' @section Interpreting output:
#' Focus on ranked rows where multiple screening criteria converge:
#' - large absolute t statistic
#' - large absolute bias size
#' - small screening tail area
#'
#' The bundle is optimized for downstream `summary()` and
#' [plot_bias_interaction()] views.
#'
#' @section Typical workflow:
#' 1. Run [estimate_bias()] (or provide `mfrm_fit` here).
#' 2. Build `bias_interaction_report(...)`.
#' 3. Review `summary(out)` and visualize with [plot_bias_interaction()].
#' @return A named list with bias-interaction plotting/report components. Class:
#'   `mfrm_bias_interaction`.
#' @seealso [estimate_bias()], [build_fixed_reports()], [plot_bias_interaction()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' out <- bias_interaction_report(bias, top_n = 10)
#' summary(out)
#' p_bi <- plot(out, draw = FALSE)
#' class(p_bi)
#' @export
bias_interaction_report <- function(x,
                                    diagnostics = NULL,
                                    facet_a = NULL,
                                    facet_b = NULL,
                                    interaction_facets = NULL,
                                    max_abs = 10,
                                    omit_extreme = TRUE,
                                    max_iter = 4,
                                    tol = 1e-3,
                                    top_n = 50,
                                    abs_t_warn = 2,
                                    abs_bias_warn = 0.5,
                                    p_max = 0.05,
                                    sort_by = c("abs_t", "abs_bias", "prob")) {
  out <- with_legacy_name_warning_suppressed(
    table13_bias_plot_export(
      x = x,
      diagnostics = diagnostics,
      facet_a = facet_a,
      facet_b = facet_b,
      interaction_facets = interaction_facets,
      max_abs = max_abs,
      omit_extreme = omit_extreme,
      max_iter = max_iter,
      tol = tol,
      top_n = top_n,
      abs_t_warn = abs_t_warn,
      abs_bias_warn = abs_bias_warn,
      p_max = p_max,
      sort_by = sort_by
    )
  )
  as_mfrm_bundle(out, "mfrm_bias_interaction")
}

#' Build a bias-iteration report
#'
#' @inheritParams bias_interaction_report
#' @param top_n Maximum number of iteration rows to keep in preview-oriented
#'   summaries. The full iteration table is always returned.
#'
#' @details
#' This report focuses on the recalibration path used by [estimate_bias()].
#' It provides a package-native counterpart to legacy iteration printouts by
#' exposing the iteration table, convergence summary, and orientation audit in
#' one bundle.
#'
#' @return A named list with:
#' - `table`: iteration history
#' - `summary`: one-row convergence summary
#' - `orientation_audit`: interaction-facet sign audit
#' - `settings`: resolved reporting options
#'
#' @seealso [estimate_bias()], [bias_interaction_report()], [build_fixed_reports()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' out <- bias_iteration_report(fit, diagnostics = diag, facet_a = "Rater", facet_b = "Criterion")
#' summary(out)
#' @export
bias_iteration_report <- function(x,
                                  diagnostics = NULL,
                                  facet_a = NULL,
                                  facet_b = NULL,
                                  interaction_facets = NULL,
                                  max_abs = 10,
                                  omit_extreme = TRUE,
                                  max_iter = 4,
                                  tol = 1e-3,
                                  top_n = 10) {
  bias_results <- normalize_bias_plot_input(
    x = x,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )

  iter_tbl <- as.data.frame(bias_results$iteration %||% data.frame(), stringsAsFactors = FALSE)
  summary_tbl <- if (nrow(iter_tbl) == 0) {
    data.frame()
  } else {
    tail_row <- iter_tbl[nrow(iter_tbl), , drop = FALSE]
    data.frame(
      InteractionFacets = paste(as.character(bias_results$interaction_facets %||% character(0)), collapse = " x "),
      Iterations = nrow(iter_tbl),
      FinalMaxLogitChange = suppressWarnings(as.numeric(tail_row$MaxLogitChange[1])),
      FinalBiasCells = suppressWarnings(as.numeric(tail_row$BiasCells[1])),
      FinalMaxScoreResidual = suppressWarnings(as.numeric(tail_row$MaxScoreResidual[1])),
      Converged = isTRUE(abs(suppressWarnings(as.numeric(tail_row$MaxLogitChange[1]))) < tol),
      MixedSign = isTRUE(bias_results$mixed_sign),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    table = iter_tbl,
    summary = summary_tbl,
    orientation_audit = as.data.frame(bias_results$orientation_audit %||% data.frame(), stringsAsFactors = FALSE),
    settings = list(
      tol = tol,
      max_iter = max_iter,
      top_n = top_n
    ),
    direction_note = as.character(bias_results$direction_note %||% ""),
    recommended_action = as.character(bias_results$recommended_action %||% "")
  )
  as_mfrm_bundle(out, "mfrm_bias_iteration")
}

#' Build a bias pairwise-contrast report
#'
#' @inheritParams bias_interaction_report
#' @param target_facet Facet whose local contrasts should be compared across
#'   the paired context facet. Defaults to the first interaction facet.
#' @param context_facet Optional facet to condition on. Defaults to the other
#'   facet in a 2-way interaction.
#' @param p_max Flagging cutoff for pairwise p-values.
#'
#' @details
#' This helper exposes the pairwise contrast table that was previously only
#' reachable through fixed-width output generation. It is available only for
#' 2-way interactions. The pairwise contrast statistic uses a
#' Welch/Satterthwaite approximation and is labeled as a Rasch-Welch
#' comparison in the output metadata.
#'
#' @return A named list with:
#' - `table`: pairwise contrast rows
#' - `summary`: one-row contrast summary
#' - `orientation_audit`: interaction-facet sign audit
#' - `settings`: resolved reporting options
#'
#' @seealso [estimate_bias()], [bias_interaction_report()], [build_fixed_reports()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' out <- bias_pairwise_report(fit, diagnostics = diag, facet_a = "Rater", facet_b = "Criterion")
#' summary(out)
#' @export
bias_pairwise_report <- function(x,
                                 diagnostics = NULL,
                                 facet_a = NULL,
                                 facet_b = NULL,
                                 interaction_facets = NULL,
                                 max_abs = 10,
                                 omit_extreme = TRUE,
                                 max_iter = 4,
                                 tol = 1e-3,
                                 target_facet = NULL,
                                 context_facet = NULL,
                                 top_n = 50,
                                 p_max = 0.05,
                                 sort_by = c("abs_t", "abs_contrast", "prob")) {
  sort_by <- match.arg(sort_by, c("abs_t", "abs_contrast", "prob"))
  bias_results <- normalize_bias_plot_input(
    x = x,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )

  spec <- extract_bias_facet_spec(bias_results)
  if (is.null(spec) || length(spec$facets) != 2L) {
    stop("`bias_pairwise_report()` is available only for 2-way interaction runs.")
  }
  if (is.null(target_facet)) {
    target_facet <- spec$facets[1]
  }
  target_facet <- as.character(target_facet[1])
  if (!target_facet %in% spec$facets) {
    stop("`target_facet` must be one of: ", paste(spec$facets, collapse = ", "))
  }
  if (is.null(context_facet)) {
    context_facet <- setdiff(spec$facets, target_facet)
  }
  context_facet <- as.character(context_facet[1])

  pair_tbl <- as.data.frame(
    calc_bias_pairwise(bias_results$table, target_facet = target_facet, context_facet = context_facet),
    stringsAsFactors = FALSE
  )
  if (nrow(pair_tbl) > 0) {
    pair_tbl$AbsT <- abs(suppressWarnings(as.numeric(pair_tbl$t)))
    pair_tbl$AbsContrast <- abs(suppressWarnings(as.numeric(pair_tbl$Contrast)))
    pair_tbl$Flag <- with(pair_tbl, is.finite(AbsT) & AbsT >= 2 | is.finite(`Prob.`) & `Prob.` <= p_max)
    ord <- switch(
      sort_by,
      abs_t = order(pair_tbl$AbsT, decreasing = TRUE, na.last = NA),
      abs_contrast = order(pair_tbl$AbsContrast, decreasing = TRUE, na.last = NA),
      prob = order(pair_tbl$`Prob.`, decreasing = FALSE, na.last = NA)
    )
    if (length(ord) > 0) {
      pair_tbl <- pair_tbl[ord, , drop = FALSE]
    }
    if (nrow(pair_tbl) > top_n) {
      pair_tbl <- pair_tbl[seq_len(top_n), , drop = FALSE]
    }
  }

  summary_tbl <- if (nrow(pair_tbl) == 0) {
    data.frame()
  } else {
    data.frame(
      TargetFacet = target_facet,
      ContextFacet = context_facet,
      Contrasts = nrow(pair_tbl),
      Flagged = sum(pair_tbl$Flag, na.rm = TRUE),
      MeanAbsContrast = mean(pair_tbl$AbsContrast, na.rm = TRUE),
      MeanAbsT = mean(pair_tbl$AbsT, na.rm = TRUE),
      MixedSign = isTRUE(bias_results$mixed_sign),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    table = pair_tbl,
    summary = summary_tbl,
    orientation_audit = as.data.frame(bias_results$orientation_audit %||% data.frame(), stringsAsFactors = FALSE),
    settings = list(
      target_facet = target_facet,
      context_facet = context_facet,
      top_n = top_n,
      p_max = p_max,
      sort_by = sort_by
    ),
    direction_note = as.character(bias_results$direction_note %||% ""),
    recommended_action = as.character(bias_results$recommended_action %||% "")
  )
  as_mfrm_bundle(out, "mfrm_bias_pairwise")
}

#' Plot bias interaction diagnostics (preferred alias)
#'
#' @inheritParams bias_interaction_report
#' @param plot Plot type: `"scatter"`, `"ranked"`, `"abs_t_hist"`,
#'   or `"facet_profile"`.
#' @param main Optional plot title override.
#' @param palette Optional named color overrides (`normal`, `flag`, `hist`,
#'   `profile`).
#' @param label_angle Label angle hint for ranked/profile labels.
#' @param preset Visual preset (`"standard"`, `"publication"`, or `"compact"`).
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' Visualization front-end for [bias_interaction_report()] with multiple views.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"scatter"` (default)}{Scatter plot of bias size (x) vs
#'     screening t-statistic (y). Points colored by flag status. Dashed reference
#'     lines at `abs_bias_warn` and `abs_t_warn`.  Use for overall triage
#'     of interaction effects.}
#'   \item{`"ranked"`}{Ranked bar chart of top `top_n` interactions sorted
#'     by `sort_by` criterion (absolute t, absolute bias, or probability).
#'     Bars colored red for flagged cells.}
#'   \item{`"abs_t_hist"`}{Histogram of absolute screening t-statistics across all
#'     interaction cells.  Dashed reference line at `abs_t_warn`.  Use for
#'     assessing the overall distribution of interaction effect sizes.}
#'   \item{`"facet_profile"`}{Per-facet-level aggregation showing mean
#'     absolute bias and flag rate.  Useful for identifying which
#'     individual facet levels drive systematic interaction patterns.}
#' }
#'
#' @section Interpreting output:
#' Start with `"scatter"` or `"ranked"` for triage, then confirm pattern shape
#' using `"abs_t_hist"` and `"facet_profile"`.
#'
#' Consistent flags across multiple views are stronger screening signals of
#' systematic interaction bias than a single extreme row, but they do not by
#' themselves establish formal inferential evidence.
#'
#' @section Typical workflow:
#' 1. Estimate bias with [estimate_bias()] or pass `mfrm_fit` directly.
#' 2. Plot with `plot = "ranked"` for top interactions.
#' 3. Cross-check using `plot = "scatter"` and `plot = "facet_profile"`.
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [bias_interaction_report()], [estimate_bias()], [plot_displacement()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_bias_interaction(
#'   fit,
#'   diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
#'   facet_a = "Rater",
#'   facet_b = "Criterion",
#'   preset = "publication",
#'   draw = FALSE
#' )
#' @export
plot_bias_interaction <- function(x,
                                  plot = c("scatter", "ranked", "abs_t_hist", "facet_profile"),
                                  diagnostics = NULL,
                                  facet_a = NULL,
                                  facet_b = NULL,
                                  interaction_facets = NULL,
                                  top_n = 40,
                                  abs_t_warn = 2,
                                  abs_bias_warn = 0.5,
                                  p_max = 0.05,
                                  sort_by = c("abs_t", "abs_bias", "prob"),
                                  main = NULL,
                                  palette = NULL,
                                  label_angle = 45,
                                  preset = c("standard", "publication", "compact"),
                                  draw = TRUE) {
  with_legacy_name_warning_suppressed(
    plot_table13_bias(
      x = x,
      plot = plot,
      diagnostics = diagnostics,
      facet_a = facet_a,
      facet_b = facet_b,
      interaction_facets = interaction_facets,
      top_n = top_n,
      abs_t_warn = abs_t_warn,
      abs_bias_warn = abs_bias_warn,
      p_max = p_max,
      sort_by = sort_by,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset,
      draw = draw
    )
  )
}

#' Build APA text outputs from model results
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param bias_results Optional output from [estimate_bias()].
#' @param context Optional named list for report context.
#' @param whexact Use exact ZSTD transformation.
#'
#' @details
#' `context` is an optional named list for narrative customization.
#' Frequently used fields include:
#' - `assessment`, `setting`, `scale_desc`
#' - `rater_training`, `raters_per_response`
#' - `rater_facet` (used for targeted reliability note text)
#' - `line_width` (optional text wrapping width for `report_text`; default = 92)
#'
#' Output text includes residual-PCA screening commentary if PCA diagnostics are
#' available in `diagnostics`.
#'
#' By default, `report_text` includes:
#' - model/data design summary (N, facet counts, scale range)
#' - optimization/convergence metrics (`Converged`, `Iterations`, `LogLik`, `AIC`, `BIC`)
#' - anchor/constraint summary (`noncenter_facet`, anchored levels, group anchors, dummy facets)
#' - category/threshold diagnostics (including disordered-step details when present)
#' - overall fit, misfit count, and top misfit levels
#' - facet reliability/separation, residual PCA summary, and bias-screen counts
#'
#' @section Interpreting output:
#' - `report_text`: manuscript-draft narrative covering Method (model
#'   specification, estimation, convergence) and Results (global fit,
#'   facet separation/reliability, misfit triage, category diagnostics,
#'   residual-PCA screening, bias screening).  Written in third-person past tense
#'   following APA 7th edition conventions, but still intended for human review.
#' - `table_figure_notes`: reusable draft note blocks for table/figure appendices.
#' - `table_figure_captions`: draft caption candidates aligned to generated outputs.
#'
#' When bias results or PCA diagnostics are not supplied, those sections
#' are omitted from the narrative rather than producing placeholder text.
#'
#' @section Typical workflow:
#' 1. Build diagnostics (and optional bias results).
#' 2. Run `build_apa_outputs(...)`.
#' 3. Check `summary(apa)` for completeness.
#' 4. Insert `apa$report_text` and note/caption fields into manuscript drafts
#'    after checking the listed cautions.
#'
#' @section Context template:
#' A minimal `context` list can include fields such as:
#' - `assessment`: name of the assessment task
#' - `setting`: administration context
#' - `scale_desc`: short description of the score scale
#' - `rater_facet`: rater facet label used in narrative reliability text
#'
#' @return
#' An object of class `mfrm_apa_outputs` with:
#' - `report_text`: APA-style Method/Results draft prose
#' - `table_figure_notes`: consolidated draft notes for tables/visuals
#' - `table_figure_captions`: draft caption candidates without figure numbering
#' - `section_map`: package-native section table for manuscript assembly
#' - `contract`: structured APA reporting contract used for downstream checks
#'
#' @seealso [build_visual_summaries()], [estimate_bias()],
#'   [reporting_checklist()], [mfrmr_reporting_and_apa]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(
#'   fit,
#'   diag,
#'   context = list(
#'     assessment = "Toy writing task",
#'     setting = "Demonstration dataset",
#'     scale_desc = "0-2 rating scale",
#'     rater_facet = "Rater"
#'   )
#' )
#' names(apa)
#' class(apa)
#' s_apa <- summary(apa)
#' s_apa$overview
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' head(chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#' cat(apa$report_text)
#' apa$section_map[, c("SectionId", "Available")]
#' @export
build_apa_outputs <- function(fit,
                              diagnostics,
                              bias_results = NULL,
                              context = list(),
                              whexact = FALSE) {
  contract <- build_apa_reporting_contract(
    res = fit,
    diagnostics = diagnostics,
    bias_results = bias_results,
    context = context,
    whexact = whexact
  )

  out <- list(
    report_text = structure(
      as.character(contract$report_text),
      class = c("mfrm_apa_text", "character")
    ),
    table_figure_notes = as.character(contract$note_text),
    table_figure_captions = as.character(contract$caption_text),
    section_map = as.data.frame(contract$section_table %||% data.frame(), stringsAsFactors = FALSE),
    contract = contract
  )
  class(out) <- c("mfrm_apa_outputs", "list")
  out
}

normalize_apa_component_text <- function(text) {
  text <- paste(as.character(text %||% character(0)), collapse = "\n")
  gsub("\\s+", " ", trimws(text))
}

apa_text_has_fragment <- function(text, fragment) {
  frag <- normalize_apa_component_text(fragment)
  if (!nzchar(frag)) return(TRUE)
  grepl(frag, normalize_apa_component_text(text), fixed = TRUE)
}

resolve_apa_output_checks <- function(object) {
  contract <- object$contract %||% NULL
  if (!inherits(contract, "mfrm_apa_contract")) {
    return(data.frame())
  }

  report_text <- as.character(object$report_text %||% "")
  note_text <- as.character(object$table_figure_notes %||% "")
  caption_text <- as.character(object$table_figure_captions %||% "")
  note_map <- contract$note_map %||% list()
  caption_map <- contract$caption_map %||% list()
  ordered_keys <- contract$ordered_keys %||% names(caption_map)

  add_check <- function(check, passed, detail) {
    data.frame(
      Check = as.character(check),
      Passed = isTRUE(passed),
      Detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  checks <- list(
    add_check(
      "Method section heading",
      grepl("^Method\\.", report_text),
      "APA narrative should begin with a Method heading."
    ),
    add_check(
      "Results section heading",
      grepl("Results\\.", report_text),
      "APA narrative should include a Results heading."
    ),
    add_check(
      "Precision caution alignment",
      if (nzchar(contract$precision$caution %||% "")) {
        apa_text_has_fragment(report_text, contract$precision$caution) ||
          apa_text_has_fragment(note_text, contract$precision$caution)
      } else {
        TRUE
      },
      if (nzchar(contract$precision$caution %||% "")) {
        "Precision caution should appear in the report text or note blocks."
      } else {
        "No extra precision caution required for this run."
      }
    ),
    add_check(
      "Bias screening note alignment",
      if (isTRUE(contract$availability$has_bias)) {
        grepl("screening", normalize_apa_component_text(report_text), fixed = TRUE) &&
          grepl("screening", normalize_apa_component_text(note_text), fixed = TRUE)
      } else {
        TRUE
      },
      if (isTRUE(contract$availability$has_bias)) {
        "Bias outputs should be labeled as screening results in both prose and notes."
      } else {
        "No bias screening block required."
      }
    ),
    add_check(
      "Residual PCA coverage",
      if (isTRUE(contract$availability$has_pca_overall) || isTRUE(contract$availability$has_pca_by_facet)) {
        grepl("Residual PCA", report_text, fixed = TRUE) &&
          grepl("Residual PCA", note_text, fixed = TRUE) &&
          grepl("Residual PCA", caption_text, fixed = TRUE)
      } else {
        TRUE
      },
      "Residual PCA availability should be reflected in prose, notes, and captions."
    ),
    add_check(
      "Note coverage",
      all(vapply(ordered_keys[ordered_keys %in% names(note_map)], function(key) {
        apa_text_has_fragment(note_text, note_map[[key]])
      }, logical(1))),
      "All note-map entries should be represented in the consolidated note text."
    ),
    add_check(
      "Caption coverage",
      all(vapply(ordered_keys[ordered_keys %in% names(caption_map)], function(key) {
        apa_text_has_fragment(caption_text, caption_map[[key]])
      }, logical(1))),
      "All caption-map entries should be represented in the consolidated caption text."
    ),
    add_check(
      "Core section coverage",
      {
        section_tbl <- as.data.frame(contract$section_table %||% data.frame(), stringsAsFactors = FALSE)
        required_sections <- c("method_design", "method_estimation", "results_scale", "results_fit_precision")
        all(required_sections %in% section_tbl$SectionId[section_tbl$Available])
      },
      "Core package-native sections should be available in the section map."
    )
  )

  if (isTRUE(contract$availability$has_interrater) && nzchar(contract$summaries$interrater_sentence %||% "")) {
    checks <- c(
      checks,
      list(
        add_check(
          "Interrater summary alignment",
          apa_text_has_fragment(report_text, contract$summaries$interrater_sentence) ||
            apa_text_has_fragment(note_text, contract$summaries$interrater_sentence),
          "Interrater agreement wording should appear in the report text or notes."
        )
      )
    )
  }

  do.call(rbind, checks)
}

#' Print APA narrative text with preserved line breaks
#'
#' @param x Character text object from `build_apa_outputs()$report_text`.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Prints APA narrative text with preserved paragraph breaks using `cat()`.
#' This is preferred over bare `print()` when you want readable multi-line
#' report output in the console.
#'
#' @section Interpreting output:
#' The printed text is the same content stored in
#' `build_apa_outputs(...)$report_text`, but with explicit paragraph breaks.
#'
#' @section Typical workflow:
#' 1. Generate `apa <- build_apa_outputs(...)`.
#' 2. Print readable narrative with `apa$report_text`.
#' 3. Use `summary(apa)` to check completeness before manuscript use.
#'
#' @return The input object (invisibly).
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(fit, diag)
#' apa$report_text
#' @export
print.mfrm_apa_text <- function(x, ...) {
  cat(as.character(x), "\n", sep = "")
  invisible(x)
}

#' Summarize APA report-output bundles
#'
#' @param object Output from [build_apa_outputs()].
#' @param top_n Maximum non-empty lines shown in each component preview.
#' @param preview_chars Maximum characters shown in each preview cell.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is a diagnostics layer for APA text products, not a replacement
#' for the full narrative.
#'
#' It reports component completeness, line/character volume, and a compact
#' preview for quick QA before manuscript insertion.
#'
#' @section Interpreting output:
#' - `overview`: total coverage across standard text components.
#' - `components`: per-component density and mention checks
#'   (including residual-PCA mentions).
#' - `sections`: package-native section coverage table.
#' - `content_checks`: contract-based alignment checks for APA drafting readiness.
#' - `overview$DraftContractPass`: the primary contract-completeness flag for
#'   draft text components.
#' - `overview$ReadyForAPA`: a backward-compatible alias of that contract flag,
#'   not a certification of inferential adequacy.
#' - `preview`: first non-empty lines for fast visual review.
#'
#' @section Typical workflow:
#' 1. Build outputs via [build_apa_outputs()].
#' 2. Run `summary(apa)` to screen for empty/short components.
#' 3. Use `apa$report_text`, `apa$table_figure_notes`,
#'    and `apa$table_figure_captions` as draft components for final-text review.
#'
#' @return An object of class `summary.mfrm_apa_outputs`.
#' @seealso [build_apa_outputs()], [summary()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(fit, diag)
#' summary(apa)
#' @export
summary.mfrm_apa_outputs <- function(object, top_n = 3, preview_chars = 160, ...) {
  if (!inherits(object, "mfrm_apa_outputs")) {
    stop("`object` must be an mfrm_apa_outputs object from build_apa_outputs().", call. = FALSE)
  }

  top_n <- max(1L, as.integer(top_n))
  preview_chars <- max(40L, as.integer(preview_chars))

  text_line_count <- function(text) {
    if (!nzchar(text)) return(0L)
    length(strsplit(text, "\n", fixed = TRUE)[[1]])
  }
  nonempty_line_count <- function(text) {
    if (!nzchar(text)) return(0L)
    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    sum(nzchar(trimws(lines)))
  }
  text_preview <- function(text, top_n, preview_chars) {
    if (!nzchar(text)) return("")
    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) return("")
    pv <- paste(utils::head(lines, n = top_n), collapse = " | ")
    if (nchar(pv) > preview_chars) {
      pv <- paste0(substr(pv, 1, preview_chars - 3), "...")
    }
    pv
  }

  components <- c("report_text", "table_figure_notes", "table_figure_captions")
  stats_tbl <- do.call(
    rbind,
    lapply(components, function(comp) {
      text_vec <- as.character(object[[comp]] %||% character(0))
      text <- paste(text_vec, collapse = "\n")
      data.frame(
        Component = comp,
        NonEmpty = nzchar(trimws(text)),
        Characters = nchar(text),
        Lines = text_line_count(text),
        NonEmptyLines = nonempty_line_count(text),
        ResidualPCA_Mentions = stringr::str_count(
          text,
          stringr::regex("Residual\\s*PCA", ignore_case = TRUE)
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  preview_tbl <- do.call(
    rbind,
    lapply(components, function(comp) {
      text_vec <- as.character(object[[comp]] %||% character(0))
      text <- paste(text_vec, collapse = "\n")
      data.frame(
        Component = comp,
        Preview = text_preview(text, top_n = top_n, preview_chars = preview_chars),
        stringsAsFactors = FALSE
      )
    })
  )

  content_checks <- resolve_apa_output_checks(object)
  total_checks <- nrow(content_checks)
  passed_checks <- if (total_checks > 0) sum(content_checks$Passed, na.rm = TRUE) else 0L
  sections_tbl <- as.data.frame(object$section_map %||% data.frame(), stringsAsFactors = FALSE)

  overview <- data.frame(
    Components = nrow(stats_tbl),
    NonEmptyComponents = sum(stats_tbl$NonEmpty),
    TotalCharacters = sum(stats_tbl$Characters),
    TotalNonEmptyLines = sum(stats_tbl$NonEmptyLines),
    Sections = nrow(sections_tbl),
    AvailableSections = if (nrow(sections_tbl) > 0) sum(sections_tbl$Available, na.rm = TRUE) else 0L,
    ContentChecks = total_checks,
    ContentChecksPassed = passed_checks,
    DraftContractPass = if (total_checks > 0) passed_checks == total_checks else TRUE,
    ReadyForAPA = if (total_checks > 0) passed_checks == total_checks else TRUE,
    stringsAsFactors = FALSE
  )

  empty_components <- stats_tbl$Component[!stats_tbl$NonEmpty]
  failed_checks <- if (total_checks > 0) content_checks$Check[!content_checks$Passed] else character(0)
  notes <- if (length(empty_components) == 0) {
    c("All standard APA text components are populated.")
  } else {
    c(paste0("Empty components: ", paste(empty_components, collapse = ", "), "."))
  }
  if (length(failed_checks) == 0) {
    notes <- c(notes, "Contract-based content checks passed.")
  } else {
    notes <- c(notes, paste0("Content checks needing review: ", paste(failed_checks, collapse = ", "), "."))
  }
  notes <- c(
    notes,
    "In this summary, ReadyForAPA/DraftContractPass indicates contract completeness for draft text components; it does not certify formal inferential adequacy."
  )
  notes <- c(notes, "Use object fields directly for full text; summary provides compact diagnostics.")

  out <- list(
    overview = overview,
    components = stats_tbl,
    sections = sections_tbl,
    content_checks = content_checks,
    preview = preview_tbl,
    notes = notes,
    top_n = top_n,
    preview_chars = preview_chars
  )
  class(out) <- "summary.mfrm_apa_outputs"
  out
}

#' @export
print.summary.mfrm_apa_outputs <- function(x, ...) {
  cat("mfrmr APA Outputs Summary\n")

  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = 0), row.names = FALSE)
  }
  if (!is.null(x$components) && nrow(x$components) > 0) {
    cat("\nComponent stats\n")
    print(round_numeric_df(as.data.frame(x$components), digits = 0), row.names = FALSE)
  }
  if (!is.null(x$sections) && nrow(x$sections) > 0) {
    cat("\nSections\n")
    print(as.data.frame(x$sections), row.names = FALSE)
  }
  if (!is.null(x$content_checks) && nrow(x$content_checks) > 0) {
    cat("\nContent checks\n")
    print(as.data.frame(x$content_checks), row.names = FALSE)
  }
  if (!is.null(x$preview) && nrow(x$preview) > 0) {
    cat("\nPreview\n")
    print(as.data.frame(x$preview), row.names = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Build APA-style table output using base R structures
#'
#' @param x A data.frame, `mfrm_fit`, diagnostics list, or bias-result list.
#' @param which Optional table selector when `x` has multiple tables.
#' @param diagnostics Optional diagnostics from [diagnose_mfrm()] (used when
#'   `x` is `mfrm_fit` and `which` targets diagnostics tables).
#' @param digits Number of rounding digits for numeric columns.
#' @param caption Optional caption text.
#' @param note Optional note text.
#' @param bias_results Optional output from [estimate_bias()] used when
#'   auto-generating APA metadata for fit-based tables.
#' @param context Optional context list forwarded when auto-generating APA
#'   metadata for fit-based tables.
#' @param whexact Logical forwarded to APA metadata helpers.
#' @param branch Output branch:
#'   `"apa"` for manuscript-oriented labels, `"facets"` for FACETS-aligned labels.
#'
#' @details
#' This helper avoids styling dependencies and returns a reproducible base
#' `data.frame` plus metadata.
#'
#' Supported `which` values:
#' - For `mfrm_fit`: `"summary"`, `"person"`, `"facets"`, `"steps"`
#' - For diagnostics list: `"overall_fit"`, `"measures"`, `"fit"`,
#'   `"reliability"`, `"facets_chisq"`, `"bias"`, `"interactions"`,
#'   `"interrater_summary"`, `"interrater_pairs"`, `"obs"`
#' - For bias-result list: `"table"`, `"summary"`, `"chi_sq"`
#'
#' @section Interpreting output:
#' - `table`: plain data.frame ready for export or further formatting.
#' - `which`: source component that produced the table.
#' - `caption`/`note`: manuscript-oriented metadata stored with the table.
#'
#' @section Typical workflow:
#' 1. Build table object with `apa_table(...)`.
#' 2. Inspect quickly with `summary(tbl)`.
#' 3. Render base preview via `plot(tbl, ...)` or export `tbl$table`.
#'
#' @return A list of class `apa_table` with fields:
#' - `table` (`data.frame`)
#' - `which`
#' - `caption`
#' - `note`
#' - `digits`
#' - `branch`, `style`
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [build_apa_outputs()],
#'   [reporting_checklist()], [mfrmr_reporting_and_apa]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' tbl <- apa_table(fit, which = "summary", caption = "Model summary", note = "Toy example")
#' tbl_facets <- apa_table(fit, which = "summary", branch = "facets")
#' summary(tbl)
#' p <- plot(tbl, draw = FALSE)
#' p_facets <- plot(tbl_facets, type = "numeric_profile", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     tbl,
#'     type = "numeric_profile",
#'     main = "APA Table Numeric Profile (Customized)",
#'     palette = c(numeric_profile = "#2b8cbe", grid = "#d9d9d9"),
#'     label_angle = 45
#'   )
#' }
#' class(tbl)
#' tbl$note
#' @export
apa_table <- function(x,
                      which = NULL,
                      diagnostics = NULL,
                      digits = 2,
                      caption = NULL,
                      note = NULL,
                      bias_results = NULL,
                      context = list(),
                      whexact = FALSE,
                      branch = c("apa", "facets")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("apa", "facets"))
  style <- ifelse(branch == "facets", "facets_manual", "apa")
  digits <- max(0L, as.integer(digits))
  table_out <- NULL
  source_type <- "data.frame"
  resolved_which <- NULL

  if (is.data.frame(x)) {
    table_out <- x
    source_type <- "data.frame"
  } else if (inherits(x, "mfrm_fit")) {
    source_type <- "mfrm_fit"
    opts <- c("summary", "person", "facets", "steps")
    diag_opts <- c(
      "overall_fit",
      "measures",
      "fit",
      "reliability",
      "facets_chisq",
      "bias",
      "interactions",
      "interrater_summary",
      "interrater_pairs",
      "obs"
    )
    if (is.null(which)) which <- "summary"
    which <- tolower(as.character(which[1]))
    resolved_which <- which

    if (which %in% opts) {
      table_out <- switch(
        which,
        summary = x$summary,
        person = x$facets$person,
        facets = x$facets$others,
        steps = x$steps
      )
    } else if (which %in% diag_opts) {
      if (is.null(diagnostics)) {
        diagnostics <- diagnose_mfrm(x, residual_pca = "none")
      }
      if (which == "interrater_summary") {
        table_out <- diagnostics$interrater$summary
      } else if (which == "interrater_pairs") {
        table_out <- diagnostics$interrater$pairs
      } else {
        table_out <- diagnostics[[which]]
      }
    } else {
      stop("Unsupported `which` for mfrm_fit. Use one of: ", paste(c(opts, diag_opts), collapse = ", "))
    }
  } else if (is.list(x) && !is.null(names(x))) {
    source_type <- "list"
    candidate <- names(x)
    if (is.null(which)) {
      pref <- c(
        "summary", "table", "overall_fit", "measures", "fit", "reliability", "facets_chisq",
        "bias", "interactions", "interrater_summary", "interrater_pairs", "obs", "chi_sq"
      )
      hit <- pref[pref %in% candidate]
      if (length(hit) == 0) {
        stop("Could not infer `which` from list input. Please specify `which`.")
      }
      which <- hit[1]
    }
    which <- as.character(which[1])
    resolved_which <- which
    if (!which %in% names(x)) {
      stop("Requested `which` not found in list input.")
    }
    table_out <- x[[which]]
  } else {
    stop("`x` must be a data.frame, mfrm_fit, or named list.")
  }

  if (is.null(table_out)) {
    table_out <- data.frame()
  }
  table_out <- as.data.frame(table_out, stringsAsFactors = FALSE)
  if (nrow(table_out) > 0) {
    num_cols <- vapply(table_out, is.numeric, logical(1))
    table_out[num_cols] <- lapply(table_out[num_cols], round, digits = digits)
  }

  resolve_contract_key <- function(which_value) {
    which_value <- tolower(as.character(which_value %||% ""))
    switch(
      which_value,
      summary = "table1",
      person = "table1",
      facets = "table1",
      measures = "table1",
      steps = "table2",
      obs = "table2",
      overall_fit = "table3",
      fit = "table3",
      reliability = "table3",
      facets_chisq = "table3",
      interrater_summary = "table3",
      interrater_pairs = "table3",
      bias = "table4",
      interactions = "table4",
      table = "table4",
      chi_sq = "table4",
      NULL
    )
  }

  if (branch == "apa" && (is.null(caption) || is.null(note)) && inherits(x, "mfrm_fit")) {
    diag_for_contract <- diagnostics
    if (is.null(diag_for_contract)) {
      diag_for_contract <- diagnose_mfrm(x, residual_pca = "none")
    }
    contract <- build_apa_reporting_contract(
      res = x,
      diagnostics = diag_for_contract,
      bias_results = bias_results,
      context = context,
      whexact = whexact
    )
    contract_key <- resolve_contract_key(resolved_which %||% which %||% source_type)
    if (is.null(caption) && !is.null(contract_key) && contract_key %in% names(contract$caption_map)) {
      caption <- contract$caption_map[[contract_key]]
    }
    if (is.null(note) && !is.null(contract_key) && contract_key %in% names(contract$note_map)) {
      note <- extract_apa_note_body(contract$note_map[[contract_key]])
    }
  }

  out <- list(
    table = table_out,
    which = if (is.null(which)) source_type else as.character(which),
    caption = if (is.null(caption)) {
      if (branch == "facets") {
        paste0("FACETS-aligned table: ", if (is.null(which)) source_type else as.character(which))
      } else {
        ""
      }
    } else {
      as.character(caption)
    },
    note = if (is.null(note)) "" else as.character(note),
    digits = digits,
    branch = branch,
    style = style
  )
  class(out) <- c(paste0("apa_table_", branch), "apa_table", class(out))
  out
}

#' @export
print.apa_table <- function(x, ...) {
  if (!is.null(x$caption) && nzchar(x$caption)) {
    cat(x$caption, "\n", sep = "")
  }
  if (is.data.frame(x$table) && nrow(x$table) > 0) {
    print(x$table, row.names = FALSE)
  } else {
    cat("<empty table>\n")
  }
  if (!is.null(x$note) && nzchar(x$note)) {
    cat("Note. ", x$note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an APA/FACETS table object
#'
#' @param object Output from [apa_table()].
#' @param digits Number of digits used for numeric summaries.
#' @param top_n Maximum numeric columns shown in `numeric_profile`.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Compact summary helper for QA of table payloads before manuscript export.
#'
#' @section Interpreting output:
#' - `overview`: table size/composition and missingness.
#' - `numeric_profile`: quick distribution summary of numeric columns.
#' - `caption`/`note`: text metadata readiness.
#'
#' @section Typical workflow:
#' 1. Build table with [apa_table()].
#' 2. Run `summary(tbl)` and inspect `overview`.
#' 3. Use [plot.apa_table()] for quick numeric checks if needed.
#'
#' @return An object of class `summary.apa_table`.
#' @seealso [apa_table()], [plot()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' tbl <- apa_table(fit, which = "summary")
#' summary(tbl)
#' @export
summary.apa_table <- function(object, digits = 3, top_n = 8, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  tbl <- as.data.frame(object$table %||% data.frame(), stringsAsFactors = FALSE)

  num_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  numeric_profile <- data.frame()
  if (length(num_cols) > 0) {
    numeric_profile <- do.call(
      rbind,
      lapply(num_cols, function(nm) {
        vals <- suppressWarnings(as.numeric(tbl[[nm]]))
        vals <- vals[is.finite(vals)]
        data.frame(
          Column = nm,
          N = length(vals),
          Mean = if (length(vals) > 0) mean(vals) else NA_real_,
          SD = if (length(vals) > 1) stats::sd(vals) else NA_real_,
          Min = if (length(vals) > 0) min(vals) else NA_real_,
          Max = if (length(vals) > 0) max(vals) else NA_real_,
          stringsAsFactors = FALSE
        )
      })
    )
    numeric_profile <- numeric_profile |>
      dplyr::arrange(dplyr::desc(.data$SD), .data$Column) |>
      dplyr::slice_head(n = top_n)
  }

  overview <- data.frame(
    Branch = as.character(object$branch %||% "apa"),
    Style = as.character(object$style %||% "apa"),
    Which = as.character(object$which %||% ""),
    Rows = nrow(tbl),
    Columns = ncol(tbl),
    NumericColumns = length(num_cols),
    MissingValues = sum(is.na(tbl)),
    stringsAsFactors = FALSE
  )

  out <- list(
    overview = overview,
    numeric_profile = numeric_profile,
    caption = as.character(object$caption %||% ""),
    note = as.character(object$note %||% ""),
    digits = digits
  )
  class(out) <- "summary.apa_table"
  out
}

#' @export
print.summary.apa_table <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("APA Table Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$caption) && nzchar(x$caption)) {
    cat("\nCaption\n")
    cat(" - ", x$caption, "\n", sep = "")
  }
  if (!is.null(x$note) && nzchar(x$note)) {
    cat("\nNote\n")
    cat(" - ", x$note, "\n", sep = "")
  }
  if (!is.null(x$numeric_profile) && nrow(x$numeric_profile) > 0) {
    cat("\nNumeric profile\n")
    print(round_numeric_df(as.data.frame(x$numeric_profile), digits = digits), row.names = FALSE)
  }
  invisible(x)
}

#' Plot an APA/FACETS table object using base R
#'
#' @param x Output from [apa_table()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"numeric_profile"` (column means) or
#'   `"first_numeric"` (distribution of the first numeric column).
#' @param main Optional title override.
#' @param palette Optional named color overrides.
#' @param label_angle Axis-label rotation angle for bar-type plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Quick visualization helper for numeric columns in [apa_table()] output.
#' It is intended for table QA and exploratory checks, not final publication
#' graphics.
#'
#' @section Interpreting output:
#' - `"numeric_profile"`: compares column means to spot scale/centering mismatches.
#' - `"first_numeric"`: checks distribution shape of the first numeric column.
#'
#' @section Typical workflow:
#' 1. Build table with [apa_table()].
#' 2. Run `summary(tbl)` for metadata.
#' 3. Use `plot(tbl, type = "numeric_profile")` for quick numeric QC.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [apa_table()], [summary()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' tbl <- apa_table(fit, which = "summary")
#' p <- plot(tbl, draw = FALSE)
#' p2 <- plot(tbl, type = "first_numeric", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     tbl,
#'     type = "numeric_profile",
#'     main = "APA Numeric Profile (Customized)",
#'     palette = c(numeric_profile = "#2b8cbe", grid = "#d9d9d9"),
#'     label_angle = 45
#'   )
#' }
#' }
#' @export
plot.apa_table <- function(x,
                           y = NULL,
                           type = c("numeric_profile", "first_numeric"),
                           main = NULL,
                           palette = NULL,
                           label_angle = 45,
                           draw = TRUE,
                           ...) {
  type <- match.arg(tolower(as.character(type[1])), c("numeric_profile", "first_numeric"))
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) {
    stop("`x$table` is empty.")
  }
  num_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  if (length(num_cols) == 0) {
    stop("`x$table` has no numeric columns to plot.")
  }
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      numeric_profile = "#1f78b4",
      first_numeric = "#33a02c",
      grid = "#ececec"
    )
  )

  if (type == "numeric_profile") {
    vals <- vapply(num_cols, function(nm) {
      v <- suppressWarnings(as.numeric(tbl[[nm]]))
      mean(v[is.finite(v)])
    }, numeric(1))
    ord <- order(abs(vals), decreasing = TRUE, na.last = NA)
    vals <- vals[ord]
    labels <- num_cols[ord]
    plot_title <- if (is.null(main)) "APA table numeric profile (column means)" else as.character(main[1])
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["numeric_profile"],
        main = plot_title,
        ylab = "Mean",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    out <- new_mfrm_plot_data(
      "apa_table",
      list(
        plot = "numeric_profile",
        column = labels,
        mean = vals,
        title = plot_title,
        subtitle = "Column-wise numeric means for manuscript triage",
        legend = new_plot_legend("Column mean", "summary", "bar", pal["numeric_profile"]),
        reference_lines = new_reference_lines("h", 0, "Zero reference", "dashed", "reference")
      )
    )
    return(invisible(out))
  }

  nm <- num_cols[1]
  vals <- suppressWarnings(as.numeric(tbl[[nm]]))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    stop("First numeric column does not contain finite values.")
  }
  plot_title <- if (is.null(main)) paste0("Distribution of ", nm) else as.character(main[1])
  if (isTRUE(draw)) {
    graphics::hist(
      x = vals,
      breaks = "FD",
      col = pal["first_numeric"],
      border = "white",
      main = plot_title,
      xlab = nm,
      ylab = "Count"
    )
  }
  out <- new_mfrm_plot_data(
    "apa_table",
    list(
      plot = "first_numeric",
      column = nm,
      values = vals,
      title = plot_title,
      subtitle = "Distribution of the first numeric APA table column",
      legend = new_plot_legend("Histogram", "distribution", "fill", pal["first_numeric"]),
      reference_lines = new_reference_lines()
    )
  )
  invisible(out)
}

#' List literature-based warning threshold profiles
#'
#' @return An object of class `mfrm_threshold_profiles` with
#'   `profiles` (`strict`, `standard`, `lenient`) and `pca_reference_bands`.
#' @details
#' Use this function to inspect available profile presets before calling
#' [build_visual_summaries()].
#'
#' `profiles` contains thresholds used by warning logic
#' (sample size, fit ratios, PCA cutoffs, etc.).
#' `pca_reference_bands` contains literature-oriented descriptive bands used in
#' summary text.
#'
#' @section Interpreting output:
#' - `profiles`: numeric threshold presets (`strict`, `standard`, `lenient`).
#' - `pca_reference_bands`: narrative reference bands for PCA interpretation.
#'
#' @section Typical workflow:
#' 1. Review presets with `mfrm_threshold_profiles()`.
#' 2. Pick a default profile for project policy.
#' 3. Override only selected fields in [build_visual_summaries()] when needed.
#'
#' @seealso [build_visual_summaries()]
#' @examples
#' profiles <- mfrm_threshold_profiles()
#' names(profiles)
#' names(profiles$profiles)
#' class(profiles)
#' s_profiles <- summary(profiles)
#' s_profiles$overview
#' @export
mfrm_threshold_profiles <- function() {
  out <- warning_threshold_profiles()
  class(out) <- c("mfrm_threshold_profiles", "list")
  out
}

#' Summarize threshold-profile presets for visual warning logic
#'
#' @param object Output from [mfrm_threshold_profiles()].
#' @param digits Number of digits used for numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Summarizes available warning presets and their PCA reference bands used by
#' [build_visual_summaries()].
#'
#' @section Interpreting output:
#' - `thresholds`: raw preset values by profile (`strict`, `standard`, `lenient`).
#' - `threshold_ranges`: per-threshold span across profiles (sensitivity to profile choice).
#' - `pca_reference`: literature bands used for PCA narrative labeling.
#'
#' Larger `Span` in `threshold_ranges` indicates settings that most change
#' warning behavior between strict and lenient modes.
#'
#' @section Typical workflow:
#' 1. Inspect `summary(mfrm_threshold_profiles())`.
#' 2. Choose profile (`strict` / `standard` / `lenient`) for project policy.
#' 3. Override selected thresholds in [build_visual_summaries()] only when justified.
#'
#' @return An object of class `summary.mfrm_threshold_profiles`.
#' @seealso [mfrm_threshold_profiles()], [build_visual_summaries()]
#' @examples
#' profiles <- mfrm_threshold_profiles()
#' summary(profiles)
#' @export
summary.mfrm_threshold_profiles <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_threshold_profiles")) {
    stop("`object` must be an mfrm_threshold_profiles object from mfrm_threshold_profiles().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))

  profiles <- object$profiles %||% list()
  profile_names <- names(profiles)
  if (is.null(profile_names)) profile_names <- character(0)

  threshold_names <- sort(unique(unlist(lapply(profiles, names), use.names = FALSE)))
  thresholds_tbl <- if (length(threshold_names) == 0) {
    data.frame()
  } else {
    tbl <- data.frame(Threshold = threshold_names, stringsAsFactors = FALSE)
    for (nm in profile_names) {
      vals <- vapply(
        threshold_names,
        function(key) {
          val <- profiles[[nm]][[key]]
          val <- suppressWarnings(as.numeric(val))
          ifelse(length(val) == 0, NA_real_, val[1])
        },
        numeric(1)
      )
      tbl[[nm]] <- vals
    }
    tbl
  }

  thresholds_range_tbl <- data.frame()
  if (nrow(thresholds_tbl) > 0 && length(profile_names) > 0) {
    mat <- as.matrix(thresholds_tbl[, profile_names, drop = FALSE])
    suppressWarnings(storage.mode(mat) <- "numeric")
    row_stats <- t(apply(mat, 1, function(v) {
      vv <- suppressWarnings(as.numeric(v))
      vv <- vv[is.finite(vv)]
      if (length(vv) == 0) return(c(Min = NA_real_, Median = NA_real_, Max = NA_real_, Span = NA_real_))
      c(
        Min = min(vv),
        Median = stats::median(vv),
        Max = max(vv),
        Span = max(vv) - min(vv)
      )
    }))
    thresholds_range_tbl <- data.frame(
      Threshold = thresholds_tbl$Threshold,
      row_stats,
      stringsAsFactors = FALSE
    )
  }

  band_tbl <- data.frame()
  bands <- object$pca_reference_bands %||% list()
  if (length(bands) > 0) {
    band_rows <- lapply(names(bands), function(band_name) {
      vals <- bands[[band_name]]
      if (is.null(vals) || length(vals) == 0) return(NULL)
      keys <- names(vals)
      if (is.null(keys) || length(keys) != length(vals)) {
        keys <- paste0("value_", seq_along(vals))
      }
      data.frame(
        Band = band_name,
        Key = as.character(keys),
        Value = suppressWarnings(as.numeric(vals)),
        stringsAsFactors = FALSE
      )
    })
    band_rows <- Filter(Negate(is.null), band_rows)
    if (length(band_rows) > 0) {
      band_tbl <- do.call(rbind, band_rows)
    }
  }

  overview <- data.frame(
    Profiles = length(profile_names),
    ThresholdCount = nrow(thresholds_tbl),
    PCAReferenceCount = nrow(band_tbl),
    DefaultProfile = if ("standard" %in% profile_names) "standard" else ifelse(length(profile_names) > 0, profile_names[1], ""),
    stringsAsFactors = FALSE
  )

  notes <- c(
    "Profiles tune warning strictness for build_visual_summaries().",
    "Use `thresholds` in build_visual_summaries() to override selected values."
  )
  required_profiles <- c("strict", "standard", "lenient")
  missing_profiles <- setdiff(required_profiles, profile_names)
  if (length(missing_profiles) > 0) {
    notes <- c(notes, paste0("Missing presets: ", paste(missing_profiles, collapse = ", "), "."))
  }

  out <- list(
    overview = overview,
    thresholds = thresholds_tbl,
    threshold_ranges = thresholds_range_tbl,
    pca_reference = band_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_threshold_profiles"
  out
}

#' @export
print.summary.mfrm_threshold_profiles <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrmr Threshold Profile Summary\n")

  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(as.data.frame(x$overview), row.names = FALSE)
  }
  if (!is.null(x$thresholds) && nrow(x$thresholds) > 0) {
    cat("\nProfile thresholds\n")
    print(round_numeric_df(as.data.frame(x$thresholds), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$threshold_ranges) && nrow(x$threshold_ranges) > 0) {
    cat("\nThreshold ranges across profiles\n")
    print(round_numeric_df(as.data.frame(x$threshold_ranges), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$pca_reference) && nrow(x$pca_reference) > 0) {
    cat("\nPCA reference bands\n")
    print(round_numeric_df(as.data.frame(x$pca_reference), digits = digits), row.names = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Build warning and narrative summaries for visual outputs
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param threshold_profile Threshold profile name (`strict`, `standard`, `lenient`).
#' @param thresholds Optional named overrides for profile thresholds.
#' @param summary_options Summary options for `build_visual_summary_map()`.
#' @param whexact Use exact ZSTD transformation.
#' @param branch Output branch:
#'   `"facets"` adds FACETS crosswalk metadata for manual-aligned reporting;
#'   `"original"` keeps package-native summary output.
#'
#' @details
#' This function returns visual-keyed text maps
#' to support dashboard/report rendering without hard-coding narrative strings
#' in UI code.
#'
#' `thresholds` can override any profile field by name. Common overrides:
#' - `n_obs_min`, `n_person_min`
#' - `misfit_ratio_warn`, `zstd2_ratio_warn`, `zstd3_ratio_warn`
#' - `pca_first_eigen_warn`, `pca_first_prop_warn`
#'
#' `summary_options` supports:
#' - `detail`: `"standard"` or `"detailed"`
#' - `max_facet_ranges`: max facet-range snippets shown in visual summaries
#' - `top_misfit_n`: number of top misfit entries included
#'
#' @section Interpreting output:
#' - `warning_map`: rule-triggered warning text by visual key.
#' - `summary_map`: descriptive narrative text by visual key.
#' - `warning_counts` / `summary_counts`: message-count tables for QA checks.
#'
#' @section Typical workflow:
#' 1. inspect defaults with [mfrm_threshold_profiles()]
#' 2. choose `threshold_profile` (`strict` / `standard` / `lenient`)
#' 3. optionally override selected fields via `thresholds`
#' 4. pass result maps to report/dashboard rendering logic
#'
#' @return
#' An object of class `mfrm_visual_summaries` with:
#' - `warning_map`: visual-level warning text vectors
#' - `summary_map`: visual-level descriptive text vectors
#' - `warning_counts`, `summary_counts`: message counts by visual key
#' - `crosswalk`: FACETS-reference mapping for main visual keys
#' - `branch`, `style`, `threshold_profile`: branch metadata
#'
#' @seealso [mfrm_threshold_profiles()], [build_apa_outputs()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' vis <- build_visual_summaries(fit, diag, threshold_profile = "strict")
#' vis2 <- build_visual_summaries(
#'   fit,
#'   diag,
#'   threshold_profile = "standard",
#'   thresholds = c(misfit_ratio_warn = 0.20, pca_first_eigen_warn = 2.0),
#'   summary_options = list(detail = "detailed", top_misfit_n = 5)
#' )
#' vis_facets <- build_visual_summaries(fit, diag, branch = "facets")
#' vis_facets$branch
#' summary(vis)
#' p <- plot(vis, type = "comparison", draw = FALSE)
#' p2 <- plot(vis, type = "warning_counts", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     vis,
#'     type = "comparison",
#'     draw = TRUE,
#'     main = "Warning vs Summary Counts (Customized)",
#'     palette = c(warning = "#cb181d", summary = "#3182bd"),
#'     label_angle = 45
#'   )
#' }
#' }
#' @export
build_visual_summaries <- function(fit,
                                   diagnostics,
                                   threshold_profile = "standard",
                                   thresholds = NULL,
                                   summary_options = NULL,
                                   whexact = FALSE,
                                   branch = c("original", "facets")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("original", "facets"))
  style <- ifelse(branch == "facets", "facets_manual", "original")

  warning_map <- build_visual_warning_map(
    res = fit,
    diagnostics = diagnostics,
    whexact = whexact,
    thresholds = thresholds,
    threshold_profile = threshold_profile
  )
  summary_map <- build_visual_summary_map(
    res = fit,
    diagnostics = diagnostics,
    whexact = whexact,
    options = summary_options,
    thresholds = thresholds,
    threshold_profile = threshold_profile
  )

  count_map_messages <- function(x) {
    if (is.null(x) || length(x) == 0) return(0L)
    vals <- unlist(x, use.names = FALSE)
    vals <- trimws(as.character(vals))
    sum(nzchar(vals))
  }
  to_count_table <- function(x) {
    keys <- names(x)
    if (is.null(keys) || length(keys) == 0) {
      return(tibble::tibble(Visual = character(0), Messages = integer(0)))
    }
    tibble::tibble(
      Visual = keys,
      Messages = vapply(x, count_map_messages, integer(1))
    ) |>
      dplyr::arrange(dplyr::desc(.data$Messages), .data$Visual)
  }

  crosswalk <- tibble::tibble(
    Visual = c(
      "unexpected",
      "fair_average",
      "displacement",
      "interrater",
      "facets_chisq",
      "residual_pca_overall",
      "residual_pca_by_facet"
    ),
    FACETS = c(
      "Table 4 / Table 10",
      "Table 12",
      "Table 9",
      "Inter-rater outputs",
      "Facet fixed/random chi-square",
      "Residual PCA (overall)",
      "Residual PCA (by facet)"
    )
  )

  out <- list(
    warning_map = warning_map,
    summary_map = summary_map,
    warning_counts = to_count_table(warning_map),
    summary_counts = to_count_table(summary_map),
    crosswalk = crosswalk,
    branch = branch,
    style = style,
    threshold_profile = as.character(threshold_profile[1])
  )
  out <- as_mfrm_bundle(out, "mfrm_visual_summaries")
  class(out) <- unique(c(paste0("mfrm_visual_summaries_", branch), class(out)))
  out
}

resolve_facets_contract_path <- function(contract_file = NULL) {
  if (!is.null(contract_file)) {
    path <- as.character(contract_file[1])
    if (file.exists(path)) return(path)
    stop("`contract_file` does not exist: ", path)
  }

  installed <- system.file("references", "facets_column_contract.csv", package = "mfrmr")
  if (nzchar(installed) && file.exists(installed)) return(installed)

  source_path <- file.path("inst", "references", "facets_column_contract.csv")
  if (file.exists(source_path)) return(source_path)

  stop(
    "Could not locate `facets_column_contract.csv`.\n",
    "Set `contract_file` explicitly or ensure the package was installed with `inst/references`."
  )
}

read_facets_contract <- function(contract_file = NULL, branch = c("facets", "original")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("facets", "original"))
  path <- resolve_facets_contract_path(contract_file)
  contract <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("table_id", "function_name", "object_id", "component", "required_columns")
  if (!all(need %in% names(contract))) {
    stop("FACETS contract file is missing required columns: ", paste(setdiff(need, names(contract)), collapse = ", "))
  }

  # Original branch uses compact Table 11 column names.
  if (identical(branch, "original")) {
    idx <- contract$object_id == "t11" & contract$component == "table"
    contract$required_columns[idx] <- "Count|BiasSize|SE|LowCountFlag"
  }

  list(path = path, contract = contract)
}

split_contract_tokens <- function(required_columns) {
  vals <- strsplit(as.character(required_columns[1]), "|", fixed = TRUE)[[1]]
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

contract_token_present <- function(token, columns) {
  token <- as.character(token[1])
  if (!nzchar(token)) return(TRUE)
  if (endsWith(token, "*")) {
    prefix <- substr(token, 1L, nchar(token) - 1L)
    return(any(startsWith(columns, prefix)))
  }
  token %in% columns
}

make_metric_row <- function(table_id, check, pass, actual = NA_real_, expected = NA_real_, note = "") {
  data.frame(
    Table = as.character(table_id),
    Check = as.character(check),
    Pass = if (is.na(pass)) NA else as.logical(pass),
    Actual = as.character(actual),
    Expected = as.character(expected),
    Note = as.character(note),
    stringsAsFactors = FALSE
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

build_parity_metric_audit <- function(outputs, tol = 1e-8) {
  rows <- list()

  add_row <- function(table_id, check, pass, actual = NA_real_, expected = NA_real_, note = "") {
    rows[[length(rows) + 1L]] <<- make_metric_row(table_id, check, pass, actual, expected, note)
  }

  t4 <- outputs$t4
  if (!is.null(t4) && is.data.frame(t4$summary) && nrow(t4$summary) > 0) {
    s4 <- t4$summary[1, , drop = FALSE]
    total <- safe_num(s4$TotalObservations)
    unexpected_n <- safe_num(s4$UnexpectedN)
    pct <- safe_num(s4$UnexpectedPercent)
    calc <- if (is.finite(total) && total > 0) 100 * unexpected_n / total else NA_real_
    pass <- if (is.finite(calc) && is.finite(pct)) abs(calc - pct) <= 1e-6 else NA
    add_row("T4", "UnexpectedPercent consistency", pass, pct, calc)
  }

  t10 <- outputs$t10
  if (!is.null(t10) && is.data.frame(t10$summary) && nrow(t10$summary) > 0) {
    s10 <- t10$summary[1, , drop = FALSE]
    baseline <- safe_num(s10$BaselineUnexpectedN)
    after <- safe_num(s10$AfterBiasUnexpectedN)
    reduced <- safe_num(s10$ReducedBy)
    reduced_pct <- safe_num(s10$ReducedPercent)
    calc_reduced <- if (all(is.finite(c(baseline, after)))) baseline - after else NA_real_
    calc_pct <- if (is.finite(baseline) && baseline > 0 && is.finite(reduced)) 100 * reduced / baseline else NA_real_
    pass_reduced <- if (is.finite(calc_reduced) && is.finite(reduced)) abs(calc_reduced - reduced) <= tol else NA
    pass_pct <- if (is.finite(calc_pct) && is.finite(reduced_pct)) abs(calc_pct - reduced_pct) <= 1e-6 else NA
    add_row("T10", "ReducedBy consistency", pass_reduced, reduced, calc_reduced)
    add_row("T10", "ReducedPercent consistency", pass_pct, reduced_pct, calc_pct)
  }

  t11 <- outputs$t11
  if (!is.null(t11) && is.data.frame(t11$summary) && nrow(t11$summary) > 0) {
    s11 <- t11$summary[1, , drop = FALSE]
    cells <- safe_num(s11$Cells)
    low <- safe_num(s11$LowCountCells)
    low_pct <- safe_num(s11$LowCountPercent)
    calc <- if (is.finite(cells) && cells > 0 && is.finite(low)) 100 * low / cells else NA_real_
    pass <- if (is.finite(calc) && is.finite(low_pct)) abs(calc - low_pct) <= 1e-6 else NA
    add_row("T11", "LowCountPercent consistency", pass, low_pct, calc)
  }

  t7a <- outputs$t7agree
  if (!is.null(t7a) && is.data.frame(t7a$summary) && nrow(t7a$summary) > 0) {
    s <- t7a$summary[1, , drop = FALSE]
    exact <- safe_num(s$ExactAgreement)
    expected_exact <- safe_num(s$ExpectedExactAgreement)
    adjacent <- safe_num(s$AdjacentAgreement)
    in_range <- function(v) is.finite(v) && v >= -tol && v <= 1 + tol
    add_row("T7", "ExactAgreement range", in_range(exact), exact, "[0,1]")
    add_row("T7", "ExpectedExactAgreement range", in_range(expected_exact), expected_exact, "[0,1]")
    add_row("T7", "AdjacentAgreement range", in_range(adjacent), adjacent, "[0,1]")
  }

  t7c <- outputs$t7chisq
  if (!is.null(t7c) && is.data.frame(t7c$table) && nrow(t7c$table) > 0) {
    fp <- safe_num(t7c$table$FixedProb)
    rp <- safe_num(t7c$table$RandomProb)
    in_unit <- function(v) {
      vals <- v[is.finite(v)]
      if (length(vals) == 0) return(NA)
      all(vals >= -tol & vals <= 1 + tol)
    }
    add_row("T7", "FixedProb range", in_unit(fp), "all", "[0,1]")
    add_row("T7", "RandomProb range", in_unit(rp), "all", "[0,1]")
  }

  disp <- outputs$disp
  if (!is.null(disp) && is.data.frame(disp$summary) && nrow(disp$summary) > 0) {
    s <- disp$summary[1, , drop = FALSE]
    levels_n <- safe_num(s$Levels)
    anchored <- safe_num(s$AnchoredLevels)
    flagged <- safe_num(s$FlaggedLevels)
    flagged_anch <- safe_num(s$FlaggedAnchoredLevels)
    pass1 <- if (all(is.finite(c(levels_n, anchored)))) anchored <= levels_n + tol else NA
    pass2 <- if (all(is.finite(c(levels_n, flagged)))) flagged <= levels_n + tol else NA
    pass3 <- if (all(is.finite(c(anchored, flagged_anch)))) flagged_anch <= anchored + tol else NA
    add_row("T9", "AnchoredLevels <= Levels", pass1, anchored, levels_n)
    add_row("T9", "FlaggedLevels <= Levels", pass2, flagged, levels_n)
    add_row("T9", "FlaggedAnchoredLevels <= AnchoredLevels", pass3, flagged_anch, anchored)
  }

  t81 <- outputs$t81
  if (!is.null(t81) && is.data.frame(t81$summary) && nrow(t81$summary) > 0) {
    s <- t81$summary[1, , drop = FALSE]
    cats <- safe_num(s$Categories)
    used <- safe_num(s$UsedCategories)
    pass_used <- if (all(is.finite(c(cats, used)))) used <= cats + tol else NA
    add_row("T8.1", "UsedCategories <= Categories", pass_used, used, cats)

    tt <- t81$threshold_table
    if (is.data.frame(tt) && nrow(tt) > 1 && "GapFromPrev" %in% names(tt)) {
      gaps <- safe_num(tt$GapFromPrev)
      monotonic_calc <- !any(gaps[is.finite(gaps)] < -tol)
      monotonic_flag <- isTRUE(s$ThresholdMonotonic)
      add_row("T8.1", "ThresholdMonotonic consistency", monotonic_flag == monotonic_calc, monotonic_flag, monotonic_calc)
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      Table = character(0),
      Check = character(0),
      Pass = logical(0),
      Actual = character(0),
      Expected = character(0),
      Note = character(0),
      stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

#' Build a FACETS compatibility-contract audit
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. If omitted,
#'   diagnostics are computed internally with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()]. If omitted and
#'   at least two facets exist, a 2-way bias run is computed internally.
#' @param branch Contract branch. `"facets"` checks legacy-compatible columns.
#'   `"original"` adapts branch-sensitive contracts to the package's compact
#'   naming.
#' @param contract_file Optional path to a custom contract CSV.
#' @param include_metrics If `TRUE`, run additional numerical consistency checks.
#' @param top_n_missing Number of lowest-coverage contract rows to keep in
#'   `missing_preview`.
#'
#' @details
#' This function audits produced report components against a compatibility
#' contract specification (`inst/references/facets_column_contract.csv`) and
#' returns:
#' - column-level coverage per contract row
#' - table-level coverage summaries
#' - optional metric-level consistency checks
#'
#' It is intended for compatibility-layer QA and regression auditing. It does
#' not establish external validity or software equivalence beyond the specific
#' schema/metric contract encoded in the audit file.
#'
#' Coverage interpretation in `overall`:
#' - `MeanColumnCoverage` and `MinColumnCoverage` are computed across all
#'   contract rows (unavailable rows count as 0 coverage).
#' - `MeanColumnCoverageAvailable` and `MinColumnCoverageAvailable` summarize
#'   only rows whose source component is available.
#'
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_parity_report` (`type = "column_coverage"`, `"table_coverage"`,
#' `"metric_status"`, `"metric_by_table"`).
#'
#' @section Interpreting output:
#' - `overall`: high-level compatibility-contract coverage and metric-check pass
#'   rates.
#' - `column_summary` / `column_audit`: where compatibility-schema mismatches
#'   occur.
#' - `metric_summary` / `metric_audit`: numerical consistency checks tied to the
#'   current contract.
#' - `missing_preview`: quickest path to unresolved compatibility gaps.
#'
#' @section Typical workflow:
#' 1. Run `facets_parity_report(fit, branch = "facets")`.
#' 2. Inspect `summary(parity)` and `missing_preview`.
#' 3. Patch upstream table builders, then rerun the compatibility audit.
#'
#' @return
#' An object of class `mfrm_parity_report` with:
#' - `overall`: one-row compatibility-audit summary
#' - `column_summary`: coverage summary by table ID
#' - `column_audit`: row-level contract audit
#' - `missing_preview`: lowest-coverage rows
#' - `metric_summary`: one-row metric-check summary
#' - `metric_by_table`: metric-check summary by table ID
#' - `metric_audit`: row-level metric checks
#' - `settings`: branch/contract metadata
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [build_fixed_reports()],
#'   [mfrmr_compatibility_layer]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' parity <- facets_parity_report(fit, diagnostics = diag, branch = "facets")
#' summary(parity)
#' p <- plot(parity, draw = FALSE)
#' }
#' @export
facets_parity_report <- function(fit,
                                 diagnostics = NULL,
                                 bias_results = NULL,
                                 branch = c("facets", "original"),
                                 contract_file = NULL,
                                 include_metrics = TRUE,
                                 top_n_missing = 15L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  branch <- match.arg(tolower(as.character(branch[1])), c("facets", "original"))
  include_metrics <- isTRUE(include_metrics)
  top_n_missing <- max(1L, as.integer(top_n_missing))

  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (is.null(bias_results) && length(facet_names) >= 2) {
    bias_results <- estimate_bias(
      fit = fit,
      diagnostics = diagnostics,
      facet_a = facet_names[1],
      facet_b = facet_names[2],
      max_iter = 2
    )
  }

  contract_info <- read_facets_contract(contract_file = contract_file, branch = branch)
  contract <- as.data.frame(contract_info$contract, stringsAsFactors = FALSE)

  # --- Facet-aware token filtering ---
  # The contract CSV is written for the reference data (Person, Rater, Task,

  # Criterion).  When the model has fewer facets, tokens that reference
  # non-existent facets should be excluded from coverage calculations.
  # Derive the set of "reference facets" from the contract's Table 5/subsets
  # row (tokens minus standard structural columns), then subtract the model's
  # actual facets to get the excluded set.
  model_facet_set <- c("Person", as.character(facet_names))
  subsets_row <- contract[contract$object_id == "t5" &
                            contract$component == "subsets", , drop = FALSE]
  if (nrow(subsets_row) > 0) {
    subsets_tokens <- split_contract_tokens(subsets_row$required_columns[1])
    structural_cols <- c("Subset", "Observations", "ObservationPercent")
    reference_facets <- setdiff(subsets_tokens, structural_cols)
  } else {
    reference_facets <- model_facet_set
  }
  excluded_facet_tokens <- setdiff(reference_facets, model_facet_set)

  outputs <- list(
    t1 = specifications_report(fit),
    t2 = data_quality_report(
      fit = fit,
      data = fit$prep$data,
      person = fit$config$person_col,
      facets = fit$config$facet_names,
      score = fit$config$score_col,
      weight = fit$config$weight_col
    ),
    t3 = estimation_iteration_report(fit, max_iter = 5),
    t4 = unexpected_response_table(fit, diagnostics = diagnostics, top_n = 50),
    t5 = measurable_summary_table(fit, diagnostics = diagnostics),
    t6 = subset_connectivity_report(fit, diagnostics = diagnostics),
    t62 = facet_statistics_report(fit, diagnostics = diagnostics),
    t7chisq = facets_chisq_table(fit, diagnostics = diagnostics),
    t7agree = interrater_agreement_table(fit, diagnostics = diagnostics),
    t81 = rating_scale_table(fit, diagnostics = diagnostics),
    t8bar = category_structure_report(fit, diagnostics = diagnostics),
    t8curves = category_curves_report(fit, theta_points = 101),
    out = facets_output_file_bundle(fit, diagnostics = diagnostics, include = c("graph", "score"), theta_points = 81),
    t12 = fair_average_table(fit, diagnostics = diagnostics),
    disp = displacement_table(fit, diagnostics = diagnostics)
  )
  if (!is.null(bias_results) && is.data.frame(bias_results$table) && nrow(bias_results$table) > 0) {
    outputs$t10 <- unexpected_after_bias_table(fit, bias_results, diagnostics = diagnostics, top_n = 50)
    outputs$t11 <- bias_count_table(bias_results, branch = branch)
    outputs$t13 <- bias_interaction_report(bias_results)
    outputs$t14 <- build_fixed_reports(bias_results, branch = branch)
  } else {
    outputs$t10 <- NULL
    outputs$t11 <- NULL
    outputs$t13 <- NULL
    outputs$t14 <- NULL
  }

  audit_rows <- lapply(seq_len(nrow(contract)), function(i) {
    row <- contract[i, , drop = FALSE]
    tokens <- split_contract_tokens(row$required_columns)
    # Exclude tokens for facets not in the current model
    tokens <- tokens[!tokens %in% excluded_facet_tokens]
    obj <- outputs[[row$object_id]]
    if (is.null(obj)) {
      return(data.frame(
        table_id = row$table_id,
        function_name = row$function_name,
        object_id = row$object_id,
        component = row$component,
        required_n = length(tokens),
        present_n = NA_integer_,
        coverage = NA_real_,
        available = FALSE,
        full_match = FALSE,
        status = "missing_object",
        missing = paste(tokens, collapse = " | "),
        stringsAsFactors = FALSE
      ))
    }
    comp <- obj[[row$component]]
    if (!is.data.frame(comp)) {
      return(data.frame(
        table_id = row$table_id,
        function_name = row$function_name,
        object_id = row$object_id,
        component = row$component,
        required_n = length(tokens),
        present_n = NA_integer_,
        coverage = NA_real_,
        available = FALSE,
        full_match = FALSE,
        status = "missing_component",
        missing = paste(tokens, collapse = " | "),
        stringsAsFactors = FALSE
      ))
    }
    cols <- names(comp)
    present <- vapply(tokens, contract_token_present, logical(1), columns = cols)
    missing <- tokens[!present]
    cov <- if (length(tokens) == 0) 1 else sum(present) / length(tokens)
    data.frame(
      table_id = row$table_id,
      function_name = row$function_name,
      object_id = row$object_id,
      component = row$component,
      required_n = length(tokens),
      present_n = sum(present),
      coverage = cov,
      available = TRUE,
      full_match = isTRUE(all(present)),
      status = if (isTRUE(all(present))) "match" else "partial",
      missing = paste(missing, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  column_audit <- dplyr::bind_rows(audit_rows)

  summarize_coverage <- function(v, fn) {
    vals <- suppressWarnings(as.numeric(v))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(NA_real_)
    fn(vals)
  }

  # Contract-level coverage should treat unavailable rows as zero coverage.
  # This avoids reporting perfect mean/min coverage when some contract rows
  # are entirely missing from available outputs.
  contract_coverage_values <- ifelse(
    column_audit$available %in% TRUE,
    suppressWarnings(as.numeric(column_audit$coverage)),
    0
  )
  contract_coverage_values[!is.finite(contract_coverage_values)] <- 0

  column_summary <- column_audit |>
    dplyr::group_by(.data$table_id, .data$function_name) |>
    dplyr::summarize(
      Components = dplyr::n(),
      Available = sum(.data$available, na.rm = TRUE),
      FullMatch = sum(.data$full_match, na.rm = TRUE),
      MeanCoverage = summarize_coverage(.data$coverage, mean),
      MinCoverage = summarize_coverage(.data$coverage, min),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$table_id, .data$function_name)

  missing_preview <- column_audit |>
    dplyr::filter(!.data$full_match | !.data$available) |>
    dplyr::arrange(.data$coverage, .data$table_id, .data$component) |>
    dplyr::slice_head(n = top_n_missing)

  metric_audit <- if (isTRUE(include_metrics)) {
    build_parity_metric_audit(outputs = outputs)
  } else {
    data.frame(
      Table = character(0),
      Check = character(0),
      Pass = logical(0),
      Actual = character(0),
      Expected = character(0),
      Note = character(0),
      stringsAsFactors = FALSE
    )
  }

  metric_summary <- if (nrow(metric_audit) == 0) {
    data.frame(
      Checks = 0L,
      Evaluated = 0L,
      Passed = 0L,
      Failed = 0L,
      PassRate = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    ev <- metric_audit$Pass[!is.na(metric_audit$Pass)]
    data.frame(
      Checks = nrow(metric_audit),
      Evaluated = length(ev),
      Passed = sum(ev %in% TRUE),
      Failed = sum(ev %in% FALSE),
      PassRate = if (length(ev) > 0) sum(ev %in% TRUE) / length(ev) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  metric_by_table <- if (nrow(metric_audit) == 0) {
    data.frame(
      Table = character(0),
      Checks = integer(0),
      Evaluated = integer(0),
      Passed = integer(0),
      Failed = integer(0),
      PassRate = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    metric_audit |>
      dplyr::group_by(.data$Table) |>
      dplyr::summarize(
        Checks = dplyr::n(),
        Evaluated = sum(!is.na(.data$Pass)),
        Passed = sum(.data$Pass %in% TRUE, na.rm = TRUE),
        Failed = sum(.data$Pass %in% FALSE, na.rm = TRUE),
        PassRate = ifelse(sum(!is.na(.data$Pass)) > 0, sum(.data$Pass %in% TRUE, na.rm = TRUE) / sum(!is.na(.data$Pass)), NA_real_),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Table)
  }

  mean_cov_all <- summarize_coverage(contract_coverage_values, mean)
  min_cov_all <- summarize_coverage(contract_coverage_values, min)
  mean_cov_available <- summarize_coverage(column_audit$coverage, mean)
  min_cov_available <- summarize_coverage(column_audit$coverage, min)
  contract_rows <- nrow(column_audit)
  mismatches <- sum(!column_audit$full_match, na.rm = TRUE)
  overall <- data.frame(
    Branch = branch,
    ContractRows = contract_rows,
    AvailableRows = sum(column_audit$available, na.rm = TRUE),
    FullMatchRows = sum(column_audit$full_match, na.rm = TRUE),
    ColumnMismatches = mismatches,
    ColumnMismatchRate = if (contract_rows > 0) mismatches / contract_rows else NA_real_,
    MeanColumnCoverage = mean_cov_all,
    MinColumnCoverage = min_cov_all,
    MeanColumnCoverageAvailable = mean_cov_available,
    MinColumnCoverageAvailable = min_cov_available,
    MetricChecks = metric_summary$Checks[1],
    MetricEvaluated = metric_summary$Evaluated[1],
    MetricFailed = metric_summary$Failed[1],
    MetricPassRate = metric_summary$PassRate[1],
    stringsAsFactors = FALSE
  )

  out <- list(
    overall = overall,
    column_summary = as.data.frame(column_summary, stringsAsFactors = FALSE),
    column_audit = as.data.frame(column_audit, stringsAsFactors = FALSE),
    missing_preview = as.data.frame(missing_preview, stringsAsFactors = FALSE),
    metric_summary = metric_summary,
    metric_by_table = as.data.frame(metric_by_table, stringsAsFactors = FALSE),
    metric_audit = as.data.frame(metric_audit, stringsAsFactors = FALSE),
    settings = list(
      branch = branch,
      contract_path = contract_info$path,
      intended_use = "compatibility_contract_audit",
      external_validation = FALSE,
      include_metrics = include_metrics,
      top_n_missing = top_n_missing,
      bias_included = !is.null(outputs$t10)
    )
  )
  as_mfrm_bundle(out, "mfrm_parity_report")
}

#' Build a package-native reference audit for report completeness
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. If omitted,
#'   diagnostics are computed internally with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()]. If omitted and
#'   at least two facets exist, a 2-way interaction screen is computed internally.
#' @param reference_profile Audit profile. `"core"` emphasizes package-native
#'   report contracts. `"compatibility"` exposes the manual-aligned compatibility
#'   layer used by `facets_parity_report(branch = "facets")`.
#' @param include_metrics If `TRUE`, run numerical consistency checks in addition
#'   to schema coverage checks.
#' @param top_n_attention Number of lowest-coverage components to keep in
#'   `attention_items`.
#'
#' @details
#' This function repackages the internal contract audit into package-native
#' terminology so users can review output completeness without needing external
#' manual/table numbering. It reports:
#' - component-level schema coverage
#' - numerical consistency checks for derived report tables
#' - the highest-priority attention items for follow-up
#'
#' It is an internal completeness audit for package-native outputs, not an
#' external validation study.
#'
#' Use `reference_profile = "core"` for ordinary `mfrmr` workflows.
#' Use `reference_profile = "compatibility"` only when you explicitly want to
#' inspect the compatibility layer.
#'
#' @section Interpreting output:
#' - `overall`: one-row internal audit summary with schema coverage and metric
#'   pass rate.
#' - `component_summary`: per-component coverage summary.
#' - `attention_items`: quickest list of components needing review.
#' - `metric_summary` / `metric_checks`: numerical consistency status.
#'
#' @return An object of class `mfrm_reference_audit`.
#' @seealso [facets_parity_report()], [diagnose_mfrm()], [build_fixed_reports()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' audit <- reference_case_audit(fit, diagnostics = diag)
#' summary(audit)
#' }
#' @export
reference_case_audit <- function(fit,
                                 diagnostics = NULL,
                                 bias_results = NULL,
                                 reference_profile = c("core", "compatibility"),
                                 include_metrics = TRUE,
                                 top_n_attention = 15L) {
  reference_profile <- match.arg(
    tolower(as.character(reference_profile[1] %||% "core")),
    c("core", "compatibility")
  )
  branch <- if (identical(reference_profile, "compatibility")) "facets" else "original"

  parity <- facets_parity_report(
    fit = fit,
    diagnostics = diagnostics,
    bias_results = bias_results,
    branch = branch,
    include_metrics = include_metrics,
    top_n_missing = top_n_attention
  )

  overall_src <- as.data.frame(parity$overall, stringsAsFactors = FALSE)
  overall <- tibble::tibble(
    ReferenceProfile = reference_profile,
    ContractBranch = as.character(overall_src$Branch[1] %||% branch),
    SchemaCoverage = as.numeric(overall_src$MeanColumnCoverage[1] %||% NA_real_),
    AvailableSchemaCoverage = as.numeric(overall_src$MeanColumnCoverageAvailable[1] %||% NA_real_),
    MinSchemaCoverage = as.numeric(overall_src$MinColumnCoverage[1] %||% NA_real_),
    MetricPassRate = as.numeric(overall_src$MetricPassRate[1] %||% NA_real_),
    SchemaMismatches = as.integer(overall_src$ColumnMismatches[1] %||% NA_integer_),
    AttentionItems = nrow(parity$missing_preview %||% data.frame()),
    CompatibilityLayer = if (identical(reference_profile, "compatibility")) "manual-aligned" else "package-native"
  )

  component_summary <- as.data.frame(parity$column_summary, stringsAsFactors = FALSE)
  names(component_summary) <- sub("^table_id$", "ComponentID", names(component_summary))
  names(component_summary) <- sub("^function_name$", "Builder", names(component_summary))

  attention_items <- as.data.frame(parity$missing_preview, stringsAsFactors = FALSE)
  names(attention_items) <- sub("^table_id$", "ComponentID", names(attention_items))
  names(attention_items) <- sub("^function_name$", "Builder", names(attention_items))
  names(attention_items) <- sub("^component$", "Subtable", names(attention_items))
  names(attention_items) <- sub("^coverage$", "Coverage", names(attention_items))
  names(attention_items) <- sub("^missing$", "MissingColumns", names(attention_items))

  out <- list(
    overall = overall,
    component_summary = component_summary,
    attention_items = attention_items,
    metric_summary = as.data.frame(parity$metric_summary, stringsAsFactors = FALSE),
    metric_checks = as.data.frame(parity$metric_audit, stringsAsFactors = FALSE),
    settings = list(
      reference_profile = reference_profile,
      contract_branch = branch,
      intended_use = "internal_contract_audit",
      external_validation = FALSE,
      include_metrics = isTRUE(include_metrics),
      top_n_attention = max(1L, as.integer(top_n_attention))
    ),
    parity = parity
  )
  as_mfrm_bundle(out, "mfrm_reference_audit")
}

# ============================================================================
# Differential Functioning Report
# ============================================================================

collect_bias_screening_summary <- function(diagnostics = NULL, bias_results = NULL) {
  out <- list(
    available = FALSE,
    bias_pct = NA_real_,
    flagged = NA_integer_,
    total = NA_integer_,
    inference_tier = NA_character_,
    statistic_label = "screening t",
    source = NA_character_,
    error_count = 0L,
    incomplete = FALSE,
    detail = NA_character_
  )

  extract_tbl <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
    if (is.list(x) && is.data.frame(x$table)) return(as.data.frame(x$table, stringsAsFactors = FALSE))
    if (is.list(x) && is.data.frame(x$bias_table)) return(as.data.frame(x$bias_table, stringsAsFactors = FALSE))
    NULL
  }

  compute_from_tbl <- function(tbl, source_label) {
    if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(NULL)
    t_col <- intersect(c("t_Residual", "t", "t.value", "Bias t"), names(tbl))
    if (length(t_col) == 0) return(NULL)
    t_vals <- suppressWarnings(as.numeric(tbl[[t_col[1]]]))
    t_vals <- t_vals[is.finite(t_vals)]
    if (length(t_vals) == 0) return(NULL)
    tier_col <- intersect(c("InferenceTier", "inference_tier"), names(tbl))
    metric_col <- intersect(c("ProbabilityMetric", "StatisticLabel"), names(tbl))
    list(
      available = TRUE,
      bias_pct = 100 * sum(abs(t_vals) > 2) / length(t_vals),
      flagged = sum(abs(t_vals) > 2),
      total = length(t_vals),
      inference_tier = if (length(tier_col) > 0) as.character(tbl[[tier_col[1]]][1]) else "screening",
      statistic_label = if (length(metric_col) > 0) as.character(tbl[[metric_col[1]]][1]) else "screening t",
      source = source_label
    )
  }

  diag_tbl <- NULL
  if (!is.null(diagnostics) && is.list(diagnostics) && is.data.frame(diagnostics$interactions)) {
    diag_tbl <- as.data.frame(diagnostics$interactions, stringsAsFactors = FALSE)
  }
  diag_out <- compute_from_tbl(diag_tbl, "diagnostics")
  if (!is.null(diag_out)) return(diag_out)

  if (inherits(bias_results, "mfrm_bias_collection")) {
    error_tbl <- as.data.frame(bias_results$errors %||% data.frame(), stringsAsFactors = FALSE)
    tables <- lapply(bias_results$by_pair %||% list(), extract_tbl)
    tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, tables)
    if (length(tables) > 0) {
      combined <- dplyr::bind_rows(tables)
      coll_out <- compute_from_tbl(combined, "bias_results_collection")
      if (!is.null(coll_out)) {
        coll_out$error_count <- nrow(error_tbl)
        coll_out$incomplete <- nrow(error_tbl) > 0L
        coll_out$detail <- if (nrow(error_tbl) > 0L) {
          sprintf("%d requested bias pair(s) failed during collection.", nrow(error_tbl))
        } else {
          NA_character_
        }
        return(coll_out)
      }
    }
    if (nrow(error_tbl) > 0L) {
      out$source <- "bias_results_collection"
      out$error_count <- nrow(error_tbl)
      out$incomplete <- TRUE
      out$detail <- sprintf("%d requested bias pair(s) failed during collection.", nrow(error_tbl))
      return(out)
    }
  }

  bias_out <- compute_from_tbl(extract_tbl(bias_results), "bias_results")
  if (!is.null(bias_out)) return(bias_out)

  out
}

#' Generate a differential-functioning interpretation report
#'
#' Produces APA-style narrative text interpreting the results of a differential-
#' functioning analysis or interaction table. For `method = "refit"`, the
#' report summarises the number of facet levels classified as negligible (A),
#' moderate (B), and large (C). For `method = "residual"`, it summarises
#' screening-positive results, lists the specific levels and their direction,
#' and includes a caveat about the distinction between construct-relevant
#' variation and measurement bias.
#'
#' @param dif_result Output from [analyze_dff()] / [analyze_dif()]
#'   (class `mfrm_dff` with compatibility class `mfrm_dif`) or
#'   [dif_interaction_table()] (class `mfrm_dif_interaction`).
#' @param ... Currently unused; reserved for future extensions.
#'
#' @details
#' When `dif_result` is an `mfrm_dff`/`mfrm_dif` object, the report is based on
#' the pairwise differential-functioning contrasts in `$dif_table`. When it is an
#' `mfrm_dif_interaction` object, the report uses the cell-level
#' statistics and flags from `$table`.
#'
#' For `method = "refit"`, ETS-style magnitude labels are used only when
#' subgroup calibrations were successfully linked back to a common baseline
#' scale; otherwise the report labels those contrasts as unclassified because
#' the refit difference is descriptive rather than comparable on a linked
#' logit scale. For `method = "residual"`, the report describes
#' screening-positive versus screening-negative contrasts instead of applying
#' ETS labels.
#'
#' @section Interpreting output:
#' - `$narrative`: character scalar with the full narrative text.
#' - `$counts`: named integer vector of method-appropriate counts.
#' - `$large_dif`: tibble of large ETS results (`method = "refit"`) or
#'   screening-positive contrasts/cells (`method = "residual"`).
#' - `$config`: analysis configuration inherited from the input.
#'
#' @section Typical workflow:
#' 1. Run [analyze_dff()] / [analyze_dif()] or [dif_interaction_table()].
#' 2. Pass the result to `dif_report()`.
#' 3. Print the report or extract `$narrative` for inclusion in a
#'    manuscript.
#'
#' @return Object of class `mfrm_dif_report` with `narrative`,
#'   `counts`, `large_dif`, and `config`.
#'
#' @seealso [analyze_dff()], [analyze_dif()], [dif_interaction_table()],
#'   [plot_dif_heatmap()], [build_apa_outputs()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' dif <- analyze_dff(fit, diag, facet = "Rater", group = "Group", data = toy)
#' rpt <- dif_report(dif)
#' cat(rpt$narrative)
#' @export
dif_report <- function(dif_result, ...) {
  if (inherits(dif_result, "mfrm_dff") || inherits(dif_result, "mfrm_dif")) {
    .dif_report_from_dif(dif_result)
  } else if (inherits(dif_result, "mfrm_dif_interaction")) {
    .dif_report_from_interaction(dif_result)
  } else {
    stop(
      "`dif_result` must be an `mfrm_dff`, `mfrm_dif`, or `mfrm_dif_interaction` object.",
         call. = FALSE)
  }
}

# Internal: generate report from mfrm_dff / mfrm_dif
.dif_report_from_dif <- function(dif_result) {
  cfg <- dif_result$config
  dt <- dif_result$dif_table

  facet_name <- cfg$facet
  group_name <- cfg$group
  method_label <- cfg$method %||% "refit"
  functioning_label <- cfg$functioning_label %||% "DFF"

  if (identical(method_label, "refit")) {
    n_a <- sum(dt$ETS == "A", na.rm = TRUE)
    n_b <- sum(dt$ETS == "B", na.rm = TRUE)
    n_c <- sum(dt$ETS == "C", na.rm = TRUE)
    n_total <- nrow(dt)
    n_screen_only <- sum(dt$Classification == "Linked contrast (screening only)", na.rm = TRUE)
    n_unclassified <- sum(dt$Classification == "Unclassified (insufficient linking)", na.rm = TRUE)
    n_na <- sum(is.na(dt$ETS))

    counts <- c(
      A = n_a,
      B = n_b,
      C = n_c,
      Linked_screening_only = n_screen_only,
      Unclassified = n_unclassified,
      NA_count = n_na,
      Total = n_total
    )
    large_dif <- dt[!is.na(dt$ETS) & dt$ETS == "C", , drop = FALSE]

    lines <- character()
    lines <- c(lines, paste0(
      functioning_label, " analysis was conducted for the ",
      facet_name, " facet across levels of ", group_name,
      " using the ", method_label, " method. "
    ))
    lines <- c(lines, paste0(
      "A total of ", n_total, " pairwise facet-level comparisons were evaluated. "
    ))
    lines <- c(lines, paste0(
      "Using ETS-style magnitude labels on the linked logit scale, ",
      n_a, " comparison(s) were classified as A (negligible), ",
      n_b, " as B (moderate), and ",
      n_c, " as C (large). "
    ))
    if (n_screen_only > 0) {
      lines <- c(lines, paste0(
        n_screen_only, " comparison(s) remained on a linked common scale but were retained as screening-only contrasts because the subgroup precision gate for primary reporting did not pass. "
      ))
    }
    if (n_unclassified > 0) {
      lines <- c(lines, paste0(
        n_unclassified, " comparison(s) could not be classified because subgroup refits ",
        "did not retain enough common linking anchors or failed to support a common-scale comparison. "
      ))
    }

    if (n_c > 0) {
      large_levels <- unique(as.character(large_dif$Level))
      lines <- c(lines, paste0(
        "\nThe following ", facet_name, " level(s) reached the current linked Category C threshold: ",
        paste(large_levels, collapse = ", "), ". "
      ))
      for (lev in large_levels) {
        lev_rows <- large_dif[large_dif$Level == lev, , drop = FALSE]
        for (r in seq_len(nrow(lev_rows))) {
          direction <- if (is.finite(lev_rows$Contrast[r]) && lev_rows$Contrast[r] > 0) {
            "higher"
          } else if (is.finite(lev_rows$Contrast[r]) && lev_rows$Contrast[r] < 0) {
            "lower"
          } else {
            "different"
          }
          lines <- c(lines, paste0(
            "  - ", lev, ": ",
            lev_rows$Group1[r], " vs ", lev_rows$Group2[r],
            " (contrast = ", sprintf("%.3f", lev_rows$Contrast[r]),
            " logits; ", lev_rows$Group1[r], " was ", direction, "). "
          ))
        }
      }
    } else {
      lines <- c(lines,
        "\nNo linked facet levels reached the current Category C threshold under the ETS-style labeling rule. "
      )
    }
  } else {
    class_col <- dt$Classification %||% rep(NA_character_, nrow(dt))
    n_positive <- sum(class_col == "Screen positive", na.rm = TRUE)
    n_negative <- sum(class_col == "Screen negative", na.rm = TRUE)
    n_na <- sum(is.na(class_col))
    n_total <- nrow(dt)

    counts <- c(
      Screen_positive = n_positive,
      Screen_negative = n_negative,
      Unclassified = n_na,
      Total = n_total
    )
    large_dif <- dt[class_col == "Screen positive", , drop = FALSE]

    lines <- character()
    lines <- c(lines, paste0(
      functioning_label, " screening was conducted for the ",
      facet_name, " facet across levels of ", group_name,
      " using the ", method_label, " method. "
    ))
    lines <- c(lines, paste0(
      "A total of ", n_total, " pairwise facet-level comparisons were evaluated. "
    ))
    lines <- c(lines, paste0(
      n_positive, " comparison(s) were screening-positive and ",
      n_negative, " were screening-negative based on the residual-contrast test. "
    ))
    if (n_na > 0) {
      lines <- c(lines, paste0(
        n_na, " comparison(s) were unclassified because of sparse data or unavailable statistics. "
      ))
    }

    if (n_positive > 0) {
      flagged_levels <- unique(as.character(large_dif$Level))
      lines <- c(lines, paste0(
        "\nThe following ", facet_name, " level(s) showed screening-positive residual contrasts: ",
        paste(flagged_levels, collapse = ", "), ". "
      ))
      for (lev in flagged_levels) {
        lev_rows <- large_dif[large_dif$Level == lev, , drop = FALSE]
        for (r in seq_len(nrow(lev_rows))) {
          direction <- if (is.finite(lev_rows$Contrast[r]) && lev_rows$Contrast[r] > 0) {
            "higher"
          } else if (is.finite(lev_rows$Contrast[r]) && lev_rows$Contrast[r] < 0) {
            "lower"
          } else {
            "different"
          }
          lines <- c(lines, paste0(
            "  - ", lev, ": ",
            lev_rows$Group1[r], " vs ", lev_rows$Group2[r],
            " (contrast = ", sprintf("%.3f", lev_rows$Contrast[r]),
            " on the residual scale; ", lev_rows$Group1[r], " was ", direction, "). "
          ))
        }
      }
    } else {
      lines <- c(lines,
        "\nNo pairwise contrasts were screening-positive under the residual-screening method. This does not by itself establish invariance or consistent functioning across groups. "
      )
    }
  }

  lines <- c(lines, paste0(
    "\nNote: The presence of differential functioning does not necessarily indicate measurement ",
    "bias. Differential functioning may reflect construct-relevant variation ",
    "(e.g., true group differences in the attribute being measured) rather ",
    "than unwanted measurement bias. Substantive review is recommended to ",
    "distinguish between these possibilities (cf. Eckes, 2011; McNamara & ",
    "Knoch, 2012)."
  ))

  narrative <- paste(lines, collapse = "")

  out <- list(
    narrative = narrative,
    counts = counts,
    large_dif = tibble::as_tibble(large_dif),
    config = cfg
  )
  class(out) <- c("mfrm_dif_report", class(out))
  out
}

# Internal: generate report from mfrm_dif_interaction
.dif_report_from_interaction <- function(dif_result) {
  cfg <- dif_result$config
  int_tbl <- dif_result$table

  facet_name <- cfg$facet
  group_name <- cfg$group
  functioning_label <- cfg$functioning_label %||% "DFF"

  n_total <- nrow(int_tbl)
  n_sparse <- sum(int_tbl$sparse, na.rm = TRUE)
  n_flag_t <- sum(int_tbl$flag_t == TRUE, na.rm = TRUE)
  n_flag_bias <- sum(int_tbl$flag_bias == TRUE, na.rm = TRUE)

  counts <- c(
    Total = n_total, Sparse = n_sparse,
    Flag_t = n_flag_t, Flag_bias = n_flag_bias
  )

  flagged_rows <- int_tbl[
    (!is.na(int_tbl$flag_t) & int_tbl$flag_t) |
    (!is.na(int_tbl$flag_bias) & int_tbl$flag_bias), , drop = FALSE
  ]

  lines <- character()
  lines <- c(lines, paste0(
    functioning_label, " interaction screening was conducted for the ",
    facet_name, " facet across levels of ", group_name,
    " using model-based residuals. "
  ))
  lines <- c(lines, paste0(
    "A total of ", n_total, " facet-level x group cells were examined. "
  ))
  if (n_sparse > 0) {
    lines <- c(lines, paste0(
      n_sparse, " cell(s) had fewer than ", cfg$min_obs,
      " observations and were flagged as sparse. "
    ))
  }
  lines <- c(lines, paste0(
    n_flag_t, " cell(s) exceeded the |t| > ", cfg$abs_t_warn,
    " threshold, and ", n_flag_bias,
    " cell(s) exceeded the |Obs-Exp average| > ", cfg$abs_bias_warn,
    " logit threshold. "
  ))

  if (nrow(flagged_rows) > 0) {
    lines <- c(lines, "\nFlagged cells:")
    for (r in seq_len(nrow(flagged_rows))) {
      lines <- c(lines, paste0(
        "  - ", flagged_rows$Level[r], " x ", flagged_rows$GroupValue[r],
        ": Obs-Exp Avg = ", sprintf("%.3f", flagged_rows$ObsExpAvg[r]),
        ", t = ", sprintf("%.2f", flagged_rows$t[r]),
        " (N = ", flagged_rows$N[r], "). "
      ))
    }
  } else {
    lines <- c(lines,
      "\nNo cells were flagged under the current screening thresholds. This does not by itself establish consistent functioning across groups. "
    )
  }

  lines <- c(lines, paste0(
    "\nNote: The presence of differential functioning does not necessarily ",
    "indicate measurement bias. Substantive review is recommended to ",
    "distinguish between construct-relevant variation and unwanted bias ",
    "(cf. Eckes, 2011; McNamara & Knoch, 2012)."
  ))

  narrative <- paste(lines, collapse = "")

  out <- list(
    narrative = narrative,
    counts = counts,
    large_dif = tibble::as_tibble(flagged_rows),
    config = cfg
  )
  class(out) <- c("mfrm_dif_report", class(out))
  out
}

#' @export
print.mfrm_dif_report <- function(x, ...) {
  cat("--- Differential Functioning Interpretation Report ---\n\n")
  cat(x$narrative, "\n")
  invisible(x)
}

#' @export
summary.mfrm_dif_report <- function(object, ...) {
  out <- list(
    narrative = object$narrative,
    counts = object$counts,
    large_dif = object$large_dif,
    config = object$config
  )
  class(out) <- "summary.mfrm_dif_report"
  out
}

#' @export
print.summary.mfrm_dif_report <- function(x, ...) {
  cat("--- Differential Functioning Report Summary ---\n")
  cat("Facet:", x$config$facet, " | Group:", x$config$group, "\n\n")
  cat("Classification counts:\n")
  print(x$counts)
  cat("\n")
  if (nrow(x$large_dif) > 0) {
    cat("Flagged levels:\n")
    print(as.data.frame(x$large_dif), row.names = FALSE, digits = 3)
  } else {
    cat("No levels flagged.\n")
  }
  invisible(x)
}

# ---- Phase 5: QC Pipeline ------------------------------------------------

#' Run automated quality control pipeline
#'
#' Integrates convergence, model fit, reliability, separation, element misfit,
#' unexpected responses, category structure, connectivity, inter-rater agreement,
#' and DIF/bias into a single pass/warn/fail report.
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()]. Computed automatically if NULL.
#' @param threshold_profile Threshold preset: `"strict"`, `"standard"` (default),
#'   or `"lenient"`.
#' @param thresholds Named list to override individual thresholds.
#' @param rater_facet Character name of the rater facet for inter-rater check
#'   (auto-detected if NULL).
#' @param include_bias If `TRUE` and bias available in diagnostics, check DIF/bias.
#' @param bias_results Optional pre-computed bias results from [estimate_bias()].
#'
#' @details
#' The pipeline evaluates 10 quality checks and assigns a verdict
#' (Pass / Warn / Fail) to each.  The overall status is the most severe
#' verdict across all checks.  Diagnostics are computed automatically via
#' [diagnose_mfrm()] if not supplied.
#'
#' Reliability and separation are used here as QC signals. In `mfrmr`,
#' `Reliability` / `Separation` are model-based facet indices and
#' `RealReliability` / `RealSeparation` provide more conservative lower bounds.
#' For `MML`, these rely on model-based `ModelSE` values for non-person facets;
#' for `JML`, they remain exploratory approximations.
#'
#' Three threshold presets are available via `threshold_profile`:
#'
#' | Aspect            | strict  | standard | lenient |
#' | :---------------- | :------ | :------- | :------ |
#' | Global fit warn   | 1.3     | 1.5      | 1.7     |
#' | Global fit fail   | 1.5     | 2.0      | 2.5     |
#' | Reliability pass  | 0.90    | 0.80     | 0.70    |
#' | Separation pass   | 3.0     | 2.0      | 1.5     |
#' | Misfit warn (pct) | 3       | 5        | 10      |
#' | Unexpected fail   | 3       | 5        | 10      |
#' | Min cat count     | 15      | 10       | 5       |
#' | Agreement pass    | 60      | 50       | 40      |
#' | Bias fail (pct)   | 5       | 10       | 15      |
#'
#' Individual thresholds can be overridden via the `thresholds` argument
#' (a named list keyed by the internal threshold names shown above).
#'
#' @section QC checks:
#' The 10 checks are:
#' \enumerate{
#'   \item **Convergence**: Did the model converge?
#'   \item **Global fit**: Infit/Outfit MnSq within the current review band.
#'   \item **Reliability**: Minimum non-person facet model reliability index.
#'   \item **Separation**: Minimum non-person facet model separation index.
#'   \item **Element misfit**: Percentage of elements with Infit/Outfit
#'         outside the current review band.
#'   \item **Unexpected responses**: Percentage of observations with
#'         large standardized residuals.
#'   \item **Category structure**: Minimum category count and threshold
#'         ordering.
#'   \item **Connectivity**: All observations in a single connected subset.
#'   \item **Inter-rater agreement**: Exact agreement percentage for the
#'         rater facet (if applicable).
#'   \item **Functioning/Bias screen**: Percentage of interaction cells that
#'         cross the screening threshold (if interaction results are available).
#' }
#'
#' @section Interpreting output:
#' - `$overall`: character string `"Pass"`, `"Warn"`, or `"Fail"`.
#' - `$verdicts`: tibble with columns `Check`, `Verdict`, `Value`, and
#'   `Threshold` for each of the 10 checks.
#' - `$details`: character vector of human-readable detail strings.
#' - `$raw_details`: named list of per-check numeric details for
#'   programmatic access.
#' - `$recommendations`: character vector of actionable suggestions for
#'   checks that did not pass.
#' - `$config`: records the threshold profile and effective thresholds.
#'
#' @section Typical workflow:
#' 1. Fit a model: `fit <- fit_mfrm(...)`.
#' 2. Optionally compute diagnostics and bias:
#'    `diag <- diagnose_mfrm(fit)`;
#'    `bias <- estimate_bias(fit, diag, ...)`.
#' 3. Run the pipeline: `qc <- run_qc_pipeline(fit, diag, bias_results = bias)`.
#' 4. Check `qc$overall` for the headline verdict.
#' 5. Review `qc$verdicts` for per-check details.
#' 6. Follow `qc$recommendations` for remediation.
#' 7. Visualize with [plot_qc_pipeline()].
#'
#' @return Object of class `mfrm_qc_pipeline` with verdicts, overall status,
#'   details, and recommendations.
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()],
#'   [mfrm_threshold_profiles()], [plot_qc_pipeline()],
#'   [plot_qc_dashboard()], [build_visual_summaries()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("study1")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' qc <- run_qc_pipeline(fit)
#' qc
#' summary(qc)
#' qc$verdicts
#' }
#' @export
run_qc_pipeline <- function(fit,
                            diagnostics = NULL,
                            threshold_profile = "standard",
                            thresholds = NULL,
                            rater_facet = NULL,
                            include_bias = TRUE,
                            bias_results = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }

  # -- compute diagnostics if needed --
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  # -- resolve threshold profile --
  threshold_profile <- match.arg(tolower(threshold_profile),
                                 c("strict", "standard", "lenient"))

  defaults_standard <- list(
    global_fit_warn = 1.5,
    global_fit_fail = 2.0,
    global_fit_low  = 0.5,
    reliability_pass = 0.80,
    reliability_warn = 0.50,
    separation_pass = 2.0,
    separation_warn = 1.0,
    misfit_warn_pct = 5,
    misfit_fail_pct = 15,
    misfit_high = 1.5,
    misfit_low  = 0.5,
    unexpected_warn_pct = 2,
    unexpected_fail_pct = 5,
    min_cat_count = 10,
    agreement_pass_pct = 50,
    agreement_warn_pct = 30,
    bias_warn_pct = 0,
    bias_fail_pct = 10
  )

  defaults_strict <- modifyList(defaults_standard, list(
    global_fit_warn = 1.3,
    global_fit_fail = 1.5,
    reliability_pass = 0.90,
    reliability_warn = 0.70,
    separation_pass = 3.0,
    separation_warn = 2.0,
    misfit_warn_pct = 3,
    misfit_fail_pct = 10,
    unexpected_warn_pct = 1,
    unexpected_fail_pct = 3,
    min_cat_count = 15,
    agreement_pass_pct = 60,
    agreement_warn_pct = 40,
    bias_warn_pct = 0,
    bias_fail_pct = 5
  ))

  defaults_lenient <- modifyList(defaults_standard, list(
    global_fit_warn = 1.7,
    global_fit_fail = 2.5,
    global_fit_low  = 0.3,
    reliability_pass = 0.70,
    reliability_warn = 0.40,
    separation_pass = 1.5,
    separation_warn = 0.5,
    misfit_warn_pct = 10,
    misfit_fail_pct = 25,
    misfit_high = 2.0,
    misfit_low  = 0.3,
    unexpected_warn_pct = 5,
    unexpected_fail_pct = 10,
    min_cat_count = 5,
    agreement_pass_pct = 40,
    agreement_warn_pct = 20,
    bias_warn_pct = 5,
    bias_fail_pct = 15
  ))

  defaults <- switch(threshold_profile,
                     strict   = defaults_strict,
                     lenient  = defaults_lenient,
                     defaults_standard)

  effective_thresholds <- modifyList(defaults, thresholds %||% list())
  thr <- effective_thresholds

  # -- helpers --
  fmt_pct <- function(x) {
    if (is.na(x)) return("NA")
    sprintf("%.1f%%", x)
  }
  fmt_num <- function(x, digits = 2) {
    if (is.na(x)) return("NA")
    formatC(x, format = "f", digits = digits)
  }

  verdicts <- character(10)
  values   <- character(10)
  thresh   <- character(10)

  details  <- character(10)
  raw_details <- list()
  recommendations <- character(0)

  # ---- Check 1: Convergence ----
  converged <- isTRUE(fit$summary$Converged)
  verdicts[1] <- if (converged) "Pass" else "Fail"
  values[1]   <- if (converged) "TRUE" else "FALSE"
  thresh[1]   <- "Converged = TRUE"
  details[1]  <- if (converged) "Model converged" else "Model did NOT converge"
  raw_details$convergence <- list(converged = converged,
                                  iterations = fit$summary$Iterations)
  if (!converged) {
    recommendations <- c(recommendations,
                         "Model did not converge. Consider increasing maxit, simplifying the model, or checking data quality.")
  }

  # ---- Check 2: Global Fit ----
  infit_global  <- as.numeric(diagnostics$overall_fit$Infit[1])
  outfit_global <- as.numeric(diagnostics$overall_fit$Outfit[1])
  if (is.na(infit_global))  infit_global  <- 1.0
  if (is.na(outfit_global)) outfit_global <- 1.0

  gf_max <- max(infit_global, outfit_global, na.rm = TRUE)
  gf_min <- min(infit_global, outfit_global, na.rm = TRUE)

  if (gf_max > thr$global_fit_fail || gf_min < thr$global_fit_low) {
    verdicts[2] <- "Fail"
  } else if (gf_max > thr$global_fit_warn) {
    verdicts[2] <- "Warn"
  } else if (gf_min < thr$global_fit_low) {
    verdicts[2] <- "Warn"
  } else {
    verdicts[2] <- "Pass"
  }
  values[2]  <- sprintf("Infit=%.2f, Outfit=%.2f", infit_global, outfit_global)
  thresh[2]  <- sprintf("[%.2f, %.2f]", thr$global_fit_low, thr$global_fit_warn)
  details[2] <- sprintf("Global Infit=%.3f, Outfit=%.3f", infit_global, outfit_global)
  raw_details$global_fit <- list(infit = infit_global, outfit = outfit_global)
  if (verdicts[2] != "Pass") {
    recommendations <- c(recommendations,
                         "Global fit indices fall outside the current review band. Investigate element-level misfit.")
  }

  # ---- Check 3: Reliability ----
  rel_tbl <- diagnostics$reliability
  if (!is.null(rel_tbl) && nrow(rel_tbl) > 0 && "Facet" %in% names(rel_tbl)) {
    rel_non_person <- rel_tbl[rel_tbl$Facet != "Person", , drop = FALSE]
    if (nrow(rel_non_person) > 0 && "Reliability" %in% names(rel_non_person)) {
      min_rel <- min(rel_non_person$Reliability, na.rm = TRUE)
    } else {
      min_rel <- NA_real_
    }
  } else {
    min_rel <- NA_real_
  }

  if (is.na(min_rel) || !is.finite(min_rel)) {
    verdicts[3] <- "Warn"
    values[3]   <- "NA"
    details[3]  <- "Model reliability could not be computed"
  } else if (min_rel >= thr$reliability_pass) {
    verdicts[3] <- "Pass"
    values[3]   <- fmt_num(min_rel)
    details[3]  <- sprintf("Min non-person model reliability = %.3f", min_rel)
  } else if (min_rel >= thr$reliability_warn) {
    verdicts[3] <- "Warn"
    values[3]   <- fmt_num(min_rel)
    details[3]  <- sprintf("Min non-person model reliability = %.3f (below %.2f)", min_rel, thr$reliability_pass)
  } else {
    verdicts[3] <- "Fail"
    values[3]   <- fmt_num(min_rel)
    details[3]  <- sprintf("Min non-person model reliability = %.3f (below %.2f)", min_rel, thr$reliability_warn)
  }
  thresh[3] <- sprintf("Pass>=%.2f, Warn>=%.2f", thr$reliability_pass, thr$reliability_warn)
  raw_details$reliability <- list(min_reliability = min_rel, table = rel_tbl)
  if (verdicts[3] == "Fail") {
    recommendations <- c(recommendations,
                         "Low facet reliability. Consider increasing sample size or reducing measurement noise.")
  }

  # ---- Check 4: Separation ----
  if (!is.null(rel_tbl) && nrow(rel_tbl) > 0 && "Facet" %in% names(rel_tbl)) {
    sep_non_person <- rel_tbl[rel_tbl$Facet != "Person", , drop = FALSE]
    if (nrow(sep_non_person) > 0 && "Separation" %in% names(sep_non_person)) {
      min_sep <- min(sep_non_person$Separation, na.rm = TRUE)
    } else {
      min_sep <- NA_real_
    }
  } else {
    min_sep <- NA_real_
  }

  if (is.na(min_sep) || !is.finite(min_sep)) {
    verdicts[4] <- "Warn"
    values[4]   <- "NA"
    details[4]  <- "Model separation could not be computed"
  } else if (min_sep >= thr$separation_pass) {
    verdicts[4] <- "Pass"
    values[4]   <- fmt_num(min_sep)
    details[4]  <- sprintf("Min non-person model separation = %.3f", min_sep)
  } else if (min_sep >= thr$separation_warn) {
    verdicts[4] <- "Warn"
    values[4]   <- fmt_num(min_sep)
    details[4]  <- sprintf("Min non-person model separation = %.3f (below %.2f)", min_sep, thr$separation_pass)
  } else {
    verdicts[4] <- "Fail"
    values[4]   <- fmt_num(min_sep)
    details[4]  <- sprintf("Min non-person model separation = %.3f (below %.2f)", min_sep, thr$separation_warn)
  }
  thresh[4] <- sprintf("Pass>=%.2f, Warn>=%.2f", thr$separation_pass, thr$separation_warn)
  raw_details$separation <- list(min_separation = min_sep)
  if (verdicts[4] == "Fail") {
    recommendations <- c(recommendations,
                         "Low facet separation. Elements may not be distinguishable. Review facet design.")
  }

  # ---- Check 5: Element Misfit ----
  fit_tbl <- diagnostics$fit
  if (!is.null(fit_tbl) && nrow(fit_tbl) > 0 &&
      all(c("Infit", "Outfit") %in% names(fit_tbl))) {
    n_elements <- nrow(fit_tbl)
    flagged <- (fit_tbl$Infit > thr$misfit_high | fit_tbl$Outfit > thr$misfit_high |
                  fit_tbl$Infit < thr$misfit_low | fit_tbl$Outfit < thr$misfit_low)
    flagged[is.na(flagged)] <- FALSE
    n_flagged <- sum(flagged)
    misfit_pct <- 100 * n_flagged / n_elements
  } else {
    n_elements <- 0
    n_flagged  <- 0
    misfit_pct <- 0
  }

  if (misfit_pct <= thr$misfit_warn_pct) {
    verdicts[5] <- "Pass"
  } else if (misfit_pct <= thr$misfit_fail_pct) {
    verdicts[5] <- "Warn"
  } else {
    verdicts[5] <- "Fail"
  }
  values[5]  <- sprintf("%d/%d (%.1f%%)", n_flagged, n_elements, misfit_pct)
  thresh[5]  <- sprintf("Pass<=%.0f%%, Fail>%.0f%%", thr$misfit_warn_pct, thr$misfit_fail_pct)
  details[5] <- sprintf("%d of %d elements misfitting (%.1f%%)", n_flagged, n_elements, misfit_pct)
  raw_details$element_misfit <- list(n_flagged = n_flagged, n_elements = n_elements,
                                     misfit_pct = misfit_pct)
  if (verdicts[5] != "Pass") {
    recommendations <- c(recommendations,
                         "Excessive element misfit detected. Review individual element fit statistics.")
  }

  # ---- Check 6: Unexpected Responses ----
  unexp_pct <- 0
  if (!is.null(diagnostics$unexpected$summary) &&
      "UnexpectedPercent" %in% names(diagnostics$unexpected$summary)) {
    unexp_pct <- as.numeric(diagnostics$unexpected$summary$UnexpectedPercent[1])
  }
  if (is.na(unexp_pct)) unexp_pct <- 0

  if (unexp_pct <= thr$unexpected_warn_pct) {
    verdicts[6] <- "Pass"
  } else if (unexp_pct <= thr$unexpected_fail_pct) {
    verdicts[6] <- "Warn"
  } else {
    verdicts[6] <- "Fail"
  }
  values[6]  <- fmt_pct(unexp_pct)
  thresh[6]  <- sprintf("Pass<=%.0f%%, Fail>%.0f%%", thr$unexpected_warn_pct, thr$unexpected_fail_pct)
  details[6] <- sprintf("%.1f%% unexpected responses", unexp_pct)
  raw_details$unexpected <- list(unexpected_pct = unexp_pct)
  if (verdicts[6] != "Pass") {
    recommendations <- c(recommendations,
                         "High unexpected response rate. Inspect unexpected_response_table() for patterns.")
  }

  # ---- Check 7: Category Structure ----
  step_est <- suppressWarnings(as.numeric(fit$steps$Estimate))
  ordered_steps <- if (length(step_est) > 1) {
    all(diff(step_est) > -sqrt(.Machine$double.eps), na.rm = TRUE)
  } else {
    TRUE
  }

  min_cat_count <- NA_real_
  category_error <- NULL
  category_available <- FALSE
  tryCatch({
    obs_df <- diagnostics$obs
    if (!is.null(obs_df) && nrow(obs_df) > 0) {
      category_available <- TRUE
      observed <- if ("Observed" %in% names(obs_df)) {
        suppressWarnings(as.numeric(obs_df$Observed))
      } else {
        suppressWarnings(as.numeric(obs_df$Score))
      }
      weights <- get_weights(obs_df)
      all_categories <- seq(fit$prep$rating_min, fit$prep$rating_max)
      counts <- numeric(length(all_categories))
      idx <- match(observed, all_categories)
      ok <- is.finite(idx) & is.finite(weights)
      if (any(ok)) {
        grouped <- split(weights[ok], idx[ok])
        counts[as.integer(names(grouped))] <- vapply(grouped, sum, numeric(1))
      }
      min_cat_count <- min(counts, na.rm = TRUE)
    }
  }, error = function(e) {
    category_error <<- conditionMessage(e)
    NULL
  })

  cat_count_ok <- is.null(category_error) && isTRUE(category_available) &&
    (is.na(min_cat_count) || min_cat_count >= thr$min_cat_count)

  if (!is.null(category_error)) {
    verdicts[7] <- "Skip"
    details[7]  <- paste0("Category counts could not be computed: ", category_error)
  } else if (!isTRUE(category_available)) {
    verdicts[7] <- "Skip"
    details[7]  <- "Category counts were not available from diagnostics$obs."
  } else if (ordered_steps && cat_count_ok) {
    verdicts[7] <- "Pass"
    details[7]  <- "Thresholds ordered"
    if (!is.na(min_cat_count)) {
      details[7] <- sprintf("Thresholds ordered, min category count = %d", as.integer(min_cat_count))
    }
  } else if (!ordered_steps && cat_count_ok) {
    verdicts[7] <- "Warn"
    details[7]  <- "Thresholds disordered"
  } else if (ordered_steps && !cat_count_ok) {
    verdicts[7] <- "Warn"
    details[7]  <- sprintf("Thresholds ordered but min category count = %d (< %d)",
                           as.integer(min_cat_count), as.integer(thr$min_cat_count))
  } else {
    verdicts[7] <- "Fail"
    details[7]  <- sprintf("Thresholds disordered, min category count = %d (< %d)",
                           as.integer(min_cat_count), as.integer(thr$min_cat_count))
  }
  values[7] <- sprintf("Ordered=%s, MinCount=%s",
                        if (ordered_steps) "Yes" else "No",
                        if (is.na(min_cat_count)) "NA" else as.character(as.integer(min_cat_count)))
  thresh[7] <- sprintf("Ordered + count>=%d", as.integer(thr$min_cat_count))
  raw_details$category_structure <- list(ordered = ordered_steps,
                                          min_cat_count = min_cat_count,
                                          available = category_available,
                                          error = category_error)
  if (verdicts[7] != "Pass") {
    recommendations <- c(recommendations,
                         "Category structure issues. Consider collapsing rating scale categories.")
  }

  # ---- Check 8: Connectivity ----
  n_subsets <- 1L
  if (!is.null(diagnostics$subsets$summary) && nrow(diagnostics$subsets$summary) > 0) {
    n_subsets <- nrow(diagnostics$subsets$summary)
  }

  if (n_subsets == 1L) {
    verdicts[8] <- "Pass"
  } else if (n_subsets == 2L) {
    verdicts[8] <- "Warn"
  } else {
    verdicts[8] <- "Fail"
  }
  values[8]  <- as.character(n_subsets)
  thresh[8]  <- "Pass=1, Warn=2, Fail>=3"
  details[8] <- sprintf("%d disjoint subset(s)", n_subsets)
  raw_details$connectivity <- list(n_subsets = n_subsets)
  if (n_subsets > 1L) {
    recommendations <- c(recommendations,
                         sprintf("Data has %d disjoint subsets. Measures are not directly comparable across subsets.", n_subsets))
  }

  # ---- Check 9: Inter-rater Agreement ----
  detected_rater <- rater_facet
  if (is.null(detected_rater)) {
    detected_rater <- infer_default_rater_facet(fit$config$facet_names)
  }

  ira_pct <- NA_real_
  ira_available <- FALSE
  ira_error <- NULL
  ira_summary <- diagnostics$interrater$summary
  summary_rater <- if (!is.null(ira_summary) &&
                       nrow(ira_summary) > 0 &&
                       "RaterFacet" %in% names(ira_summary)) {
    as.character(ira_summary$RaterFacet[1])
  } else {
    NA_character_
  }
  tryCatch({
    if (!is.null(detected_rater) && detected_rater %in% fit$config$facet_names) {
      if (!is.null(ira_summary) &&
          nrow(ira_summary) > 0 &&
          "ExactAgreement" %in% names(ira_summary) &&
          identical(summary_rater, detected_rater)) {
        ira_pct <- as.numeric(ira_summary$ExactAgreement[1]) * 100
        ira_available <- is.finite(ira_pct)
      }
      if (!ira_available) {
        ira <- interrater_agreement_table(fit, diagnostics,
                                          rater_facet = detected_rater)
        if (!is.null(ira$summary) && nrow(ira$summary) > 0 &&
            "ExactAgreement" %in% names(ira$summary)) {
          ira_pct <- as.numeric(ira$summary$ExactAgreement[1]) * 100
          ira_available <- TRUE
        }
      }
    }
  }, error = function(e) {
    ira_error <<- conditionMessage(e)
    NULL
  })

  if (!ira_available || is.na(ira_pct)) {
    verdicts[9] <- "Skip"
    values[9]   <- "NA"
    details[9]  <- if (!is.null(ira_error)) {
      paste0("Inter-rater agreement could not be computed: ", ira_error)
    } else {
      "No rater facet available or inter-rater agreement could not be computed"
    }
    thresh[9]   <- sprintf("Pass>=%.0f%%, Warn>=%.0f%%",
                           thr$agreement_pass_pct, thr$agreement_warn_pct)
  } else {
    if (ira_pct >= thr$agreement_pass_pct) {
      verdicts[9] <- "Pass"
    } else if (ira_pct >= thr$agreement_warn_pct) {
      verdicts[9] <- "Warn"
    } else {
      verdicts[9] <- "Fail"
    }
    values[9]  <- fmt_pct(ira_pct)
    thresh[9]  <- sprintf("Pass>=%.0f%%, Warn>=%.0f%%",
                          thr$agreement_pass_pct, thr$agreement_warn_pct)
    details[9] <- sprintf("Exact agreement = %.1f%%", ira_pct)
  }
  raw_details$interrater <- list(exact_agreement_pct = ira_pct,
                                  rater_facet = detected_rater,
                                  error = ira_error)
  if (verdicts[9] == "Fail") {
    recommendations <- c(recommendations,
                         "Low inter-rater agreement. Consider rater training or calibration.")
  }

  # ---- Check 10: Functioning/Bias screen ----
  bias_screen_error <- NULL
  bias_screen <- if (isTRUE(include_bias)) {
    tryCatch(
      collect_bias_screening_summary(diagnostics = diagnostics, bias_results = bias_results),
      error = function(e) {
        bias_screen_error <<- conditionMessage(e)
        NULL
      }
    )
  } else {
    NULL
  }
  bias_available <- is.list(bias_screen) && isTRUE(bias_screen$available)
  bias_pct <- if (bias_available) as.numeric(bias_screen$bias_pct) else NA_real_
  bias_incomplete <- is.list(bias_screen) && isTRUE(bias_screen$incomplete)
  bias_detail <- if (is.list(bias_screen)) as.character(bias_screen$detail %||% "") else ""

  if (!bias_available || is.na(bias_pct)) {
    verdicts[10] <- if (bias_incomplete) "Warn" else "Skip"
    values[10]   <- "NA"
    details[10]  <- if (!is.null(bias_screen_error)) {
      paste0("Functioning/bias screen failed: ", bias_screen_error)
    } else if (bias_incomplete && nzchar(bias_detail)) {
      paste0("Functioning/bias screen was incomplete: ", bias_detail)
    } else {
      "Functioning/bias screen not available"
    }
    thresh[10]   <- sprintf("Pass<=%.0f%%, Fail>%.0f%%", thr$bias_warn_pct, thr$bias_fail_pct)
  } else {
    if (bias_pct <= thr$bias_warn_pct) {
      verdicts[10] <- if (bias_incomplete) "Warn" else "Pass"
    } else if (bias_pct <= thr$bias_fail_pct) {
      verdicts[10] <- "Warn"
    } else {
      verdicts[10] <- "Fail"
    }
    values[10]  <- fmt_pct(bias_pct)
    thresh[10]  <- sprintf("Pass<=%.0f%%, Fail>%.0f%%", thr$bias_warn_pct, thr$bias_fail_pct)
    details[10] <- sprintf(
      "%.1f%% of screened interactions crossed |%s| > 2%s",
      bias_pct,
      as.character(bias_screen$statistic_label %||% "screening t"),
      if (bias_incomplete && nzchar(bias_detail)) paste0("; ", bias_detail) else ""
    )
  }
  raw_details$bias <- list(
    bias_pct = bias_pct,
    available = bias_available,
    flagged = if (bias_available) as.integer(bias_screen$flagged) else NA_integer_,
    total = if (bias_available) as.integer(bias_screen$total) else NA_integer_,
    inference_tier = if (bias_available) as.character(bias_screen$inference_tier) else NA_character_,
    statistic_label = if (bias_available) as.character(bias_screen$statistic_label) else NA_character_,
    source = if (bias_available) as.character(bias_screen$source) else NA_character_,
    incomplete = bias_incomplete,
    error_count = if (is.list(bias_screen)) as.integer(bias_screen$error_count %||% 0L) else 0L,
    detail = bias_detail,
    error = bias_screen_error
  )
  if (isTRUE(verdicts[10] == "Fail")) {
    recommendations <- c(recommendations,
                         "Many interaction cells were screen-positive. Review estimate_bias() or analyze_dff() before making substantive bias claims.")
  }

  # -- build verdicts tibble --
  verdicts_tbl <- tibble::tibble(
    Check     = c("Convergence", "Global Fit", "Reliability", "Separation",
                  "Element Misfit", "Unexpected Responses", "Category Structure",
                  "Connectivity", "Inter-rater Agreement", "Functioning/Bias Screen"),
    Verdict   = verdicts,
    Value     = values,
    Threshold = thresh,
    Detail    = details
  )

  # -- overall verdict --
  active_verdicts <- verdicts[verdicts != "Skip"]
  if (any(active_verdicts == "Fail")) {
    overall <- "Fail"
  } else if (any(active_verdicts == "Warn")) {
    overall <- "Warn"
  } else if (any(verdicts == "Skip")) {
    overall <- "Warn"
  } else {
    overall <- "Pass"
  }

  out <- list(
    verdicts = verdicts_tbl,
    overall  = overall,
    details  = raw_details,
    recommendations = recommendations,
    config   = list(threshold_profile = threshold_profile,
                    thresholds = effective_thresholds)
  )
  class(out) <- c("mfrm_qc_pipeline", "list")
  out
}

#' @export
print.mfrm_qc_pipeline <- function(x, ...) {
  cat("--- QC Pipeline ---\n")
  cat("Overall:", x$overall, "\n\n")
  vt <- x$verdicts
  markers <- ifelse(vt$Verdict == "Pass", "[PASS]",
                    ifelse(vt$Verdict == "Warn", "[WARN]",
                           ifelse(vt$Verdict == "Fail", "[FAIL]", "[SKIP]")))
  for (i in seq_len(nrow(vt))) {
    cat(sprintf("  %s %-25s %s\n", markers[i], vt$Check[i], vt$Detail[i]))
  }
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (r in x$recommendations) cat("  -", r, "\n")
  }
  invisible(x)
}

#' @export
summary.mfrm_qc_pipeline <- function(object, ...) {
  out <- list(
    verdicts = object$verdicts,
    overall  = object$overall,
    recommendations = object$recommendations,
    pass_count = sum(object$verdicts$Verdict == "Pass"),
    warn_count = sum(object$verdicts$Verdict == "Warn"),
    fail_count = sum(object$verdicts$Verdict == "Fail"),
    skip_count = sum(object$verdicts$Verdict == "Skip")
  )
  class(out) <- "summary.mfrm_qc_pipeline"
  out
}

#' @export
print.summary.mfrm_qc_pipeline <- function(x, ...) {
  cat("--- QC Pipeline Summary ---\n")
  cat("Overall:", x$overall, "\n")
  cat(sprintf("Pass: %d | Warn: %d | Fail: %d | Skip: %d\n\n",
              x$pass_count, x$warn_count, x$fail_count, x$skip_count))
  print(as.data.frame(x$verdicts), row.names = FALSE)
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (r in x$recommendations) cat("  -", r, "\n")
  }
  invisible(x)
}
