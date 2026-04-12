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
#' p_spec$data$plot
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
#' p_dq$data$plot
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
#' p_iter$data$plot
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
#' p_sub$data$plot
#' p_design$data$plot
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
#' p_fs$data$plot
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
    method = resolve_public_mfrm_method(
      summary_method = fit$summary$Method[1] %||% NA_character_,
      method_input = fit$config$method_input %||% NA_character_,
      method_used = fit$config$method %||% NA_character_
    ),
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
#' 3. `out$median_thresholds` for adjacent-threshold caveats when zero-count
#'    categories are retained.
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
#' head(out$category_table[, c("Category", "Count", "Infit", "Outfit")])
#' p_cs <- plot(out, draw = FALSE)
#' p_cs$data$plot
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
#' head(out$probabilities[, c("CurveGroup", "Theta", "Category", "Probability")])
#' p_cc <- plot(out, draw = FALSE)
#' p_cc$data$plot
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
#' p_bi$data$plot
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
#' For bounded `GPCM`, this helper is intentionally unavailable. Use
#' [reporting_checklist()], [precision_audit_report()], and the direct
#' table/plot helpers instead, and treat [gpcm_capability_matrix()] as the
#' formal boundary statement for that branch.
#'
#' By default, `report_text` includes:
#' - model/data design summary (N, facet counts, scale range)
#' - optimization/convergence metrics (`Converged`, `Iterations`, `LogLik`, `AIC`, `BIC`)
#' - anchor/constraint summary (`noncenter_facet`, anchored levels, group anchors, dummy facets)
#' - latent-regression population-model wording when `fit` has an active
#'   `population_formula`
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
#' - active latent-regression fits add a population-model section and Table 5
#'   notes/captions that distinguish conditional-normal coefficient reporting
#'   from post hoc regression on EAP/MLE scores.
#'
#' When bias results or PCA diagnostics are not supplied, those sections
#' are omitted from the narrative rather than producing placeholder text.
#'
#' @section Typical workflow:
#' 1. Build diagnostics (and optional bias results). For `RSM` / `PCM`
#'    reporting runs, prefer an `MML` fit and
#'    `diagnose_mfrm(..., diagnostic_mode = "both")`.
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
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", maxit = 200)
#' diag <- diagnose_mfrm(fit, residual_pca = "both", diagnostic_mode = "both")
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
#' s_apa <- summary(apa)
#' s_apa$overview
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' head(chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#' cat(apa$report_text)
#' apa$section_map[, c("SectionId", "Available")]
#' }
#'
#' @section Input validation:
#' `fit` must be an `mfrm_fit` object from [fit_mfrm()].
#' `diagnostics` must be an `mfrm_diagnostics` object from [diagnose_mfrm()].
#' `context` must be a list (use `NULL` or `list()` for no extra context).
#' If supplied, `bias_results` must come from [estimate_bias()] or another
#' package-native bias helper that provides a table component.
#' @export
build_apa_outputs <- function(fit,
                              diagnostics,
                              bias_results = NULL,
                              context = list(),
                              whexact = FALSE) {
  validated <- validate_apa_builder_inputs(
    fit = fit,
    diagnostics = diagnostics,
    bias_results = bias_results,
    context = context,
    helper = "build_apa_outputs()"
  )
  fit <- validated$fit
  diagnostics <- validated$diagnostics
  bias_results <- validated$bias_results
  context <- validated$context
  stop_if_gpcm_out_of_scope(fit, "build_apa_outputs()")
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

# Internal input validator shared by APA/report helpers that build the
# package-native reporting contract.
validate_apa_builder_inputs <- function(fit,
                                        diagnostics,
                                        bias_results = NULL,
                                        context = list(),
                                        helper = "build_apa_outputs()") {
  if (!inherits(fit, "mfrm_fit")) {
    stop(
      "`", helper, "` requires `fit` to be an `mfrm_fit` object returned by `fit_mfrm()`.",
      call. = FALSE
    )
  }
  if (missing(diagnostics) || is.null(diagnostics) || !inherits(diagnostics, "mfrm_diagnostics")) {
    stop(
      "`", helper, "` requires `diagnostics` to be an `mfrm_diagnostics` object returned by `diagnose_mfrm()`.",
      call. = FALSE
    )
  }

  context <- context %||% list()
  if (!is.list(context)) {
    stop(
      "`", helper, "` requires `context` to be a list. Use `NULL` or `list()` when no extra reporting context is needed.",
      call. = FALSE
    )
  }

  if (!is.null(bias_results)) {
    has_bias_table <- FALSE
    if (is.data.frame(bias_results)) {
      has_bias_table <- TRUE
    } else if (is.list(bias_results)) {
      has_bias_table <- is.data.frame(bias_results$table) || is.data.frame(bias_results$bias_table)
    }
    if (!isTRUE(has_bias_table)) {
      stop(
        "`", helper, "` requires `bias_results` to be `NULL` or a package-native bias result with a data-frame table component, such as `estimate_bias()` output.",
        call. = FALSE
      )
    }
  }

  list(
    fit = fit,
    diagnostics = diagnostics,
    bias_results = bias_results,
    context = context
  )
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

  if (isTRUE(contract$availability$has_population_model)) {
    population_summary <- contract$summaries$population_model %||% list()
    checks <- c(
      checks,
      list(
        add_check(
          "Latent-regression wording alignment",
          apa_text_has_fragment(report_text, population_summary$caution_sentence %||% "") &&
            grepl("Latent-regression population model", note_text, fixed = TRUE) &&
            grepl("documented latent-regression MML comparison scope", normalize_apa_component_text(report_text), fixed = TRUE),
          "Active latent-regression runs should state conditional-normal population-model interpretation, avoid post hoc score-regression wording, and keep ConQuest scope wording explicit."
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

summary_table_bundle_df <- function(x) {
  if (is.null(x)) return(data.frame())
  if (inherits(x, "tbl_df")) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  data.frame()
}

summary_table_bundle_text_df <- function(x, column = "Note") {
  if (is.null(x) || length(x) == 0L) return(data.frame())
  data.frame(
    stats::setNames(list(as.character(x)), column),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

summary_table_bundle_settings_df <- function(x) {
  if (is.null(x)) return(data.frame())
  bundle_settings_table(x)
}

summary_table_bundle_supported_summary_classes <- function() {
  c(
    "summary.mfrm_fit",
    "summary.mfrm_diagnostics",
    "summary.mfrm_data_description",
    "summary.mfrm_reporting_checklist",
    "summary.mfrm_apa_outputs",
    "summary.mfrm_design_evaluation",
    "summary.mfrm_signal_detection",
    "summary.mfrm_population_prediction",
    "summary.mfrm_future_branch_active_branch",
    "summary.mfrm_facets_run",
    "summary.mfrm_bias",
    "summary.mfrm_anchor_audit",
    "summary.mfrm_linking_review",
    "summary.mfrm_misfit_casebook",
    "summary.mfrm_weighting_audit",
    "summary.mfrm_unit_prediction",
    "summary.mfrm_plausible_values"
  )
}

summary_table_bundle_is_empty <- function(x) {
  is.data.frame(x) && nrow(x) == 0L && ncol(x) == 0L
}

resolve_summary_table_bundle_input <- function(x,
                                               digits = 3,
                                               top_n = 10,
                                               preview_chars = 160) {
  summary_classes <- summary_table_bundle_supported_summary_classes()
  if (inherits(x, summary_classes)) {
    cls <- class(x)
    cls <- cls[startsWith(cls, "summary.mfrm_")][1]
    return(list(
      summary = x,
      source_class = cls,
      summary_class = cls
    ))
  }

  if (inherits(x, "mfrm_fit")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_fit",
      summary_class = "summary.mfrm_fit"
    ))
  }
  if (inherits(x, "mfrm_diagnostics")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_diagnostics",
      summary_class = "summary.mfrm_diagnostics"
    ))
  }
  if (inherits(x, "mfrm_data_description")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_data_description",
      summary_class = "summary.mfrm_data_description"
    ))
  }
  if (inherits(x, "mfrm_reporting_checklist")) {
    return(list(
      summary = summary(x, top_n = top_n),
      source_class = "mfrm_reporting_checklist",
      summary_class = "summary.mfrm_reporting_checklist"
    ))
  }
  if (inherits(x, "mfrm_apa_outputs")) {
    return(list(
      summary = summary(x, top_n = top_n, preview_chars = preview_chars),
      source_class = "mfrm_apa_outputs",
      summary_class = "summary.mfrm_apa_outputs"
    ))
  }
  if (inherits(x, "mfrm_design_evaluation")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_design_evaluation",
      summary_class = "summary.mfrm_design_evaluation"
    ))
  }
  if (inherits(x, "mfrm_signal_detection")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_signal_detection",
      summary_class = "summary.mfrm_signal_detection"
    ))
  }
  if (inherits(x, "mfrm_population_prediction")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_population_prediction",
      summary_class = "summary.mfrm_population_prediction"
    ))
  }
  if (inherits(x, "mfrm_future_branch_active_branch")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_future_branch_active_branch",
      summary_class = "summary.mfrm_future_branch_active_branch"
    ))
  }
  if (inherits(x, "mfrm_facets_run")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_facets_run",
      summary_class = "summary.mfrm_facets_run"
    ))
  }
  if (inherits(x, "mfrm_bias")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_bias",
      summary_class = "summary.mfrm_bias"
    ))
  }
  if (inherits(x, "mfrm_anchor_audit")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_anchor_audit",
      summary_class = "summary.mfrm_anchor_audit"
    ))
  }
  if (inherits(x, "mfrm_linking_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_linking_review",
      summary_class = "summary.mfrm_linking_review"
    ))
  }
  if (inherits(x, "mfrm_misfit_casebook")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_misfit_casebook",
      summary_class = "summary.mfrm_misfit_casebook"
    ))
  }
  if (inherits(x, "mfrm_weighting_audit")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_weighting_audit",
      summary_class = "summary.mfrm_weighting_audit"
    ))
  }
  if (inherits(x, "mfrm_unit_prediction")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_unit_prediction",
      summary_class = "summary.mfrm_unit_prediction"
    ))
  }
  if (inherits(x, "mfrm_plausible_values")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_plausible_values",
      summary_class = "summary.mfrm_plausible_values"
    ))
  }

  stop(
    "`x` must be an mfrm_fit, mfrm_diagnostics, mfrm_data_description, ",
    "mfrm_reporting_checklist, mfrm_apa_outputs, mfrm_design_evaluation, ",
    "mfrm_signal_detection, mfrm_population_prediction, mfrm_future_branch_active_branch, ",
    "mfrm_facets_run, mfrm_bias, mfrm_anchor_audit, mfrm_linking_review, mfrm_misfit_casebook, ",
    "mfrm_weighting_audit, mfrm_unit_prediction, or ",
    "mfrm_plausible_values object, or one of their summary() outputs.",
    call. = FALSE
  )
}

summary_table_bundle_required_components <- function(summary_class) {
  switch(
    as.character(summary_class %||% NA_character_),
    "summary.mfrm_fit" = c("overview", "reporting_map"),
    "summary.mfrm_diagnostics" = c("overview", "reporting_map", "flags"),
    "summary.mfrm_data_description" = c("overview", "score_distribution"),
    "summary.mfrm_reporting_checklist" = c("overview", "action_items"),
    "summary.mfrm_apa_outputs" = c("overview", "components", "preview"),
    "summary.mfrm_design_evaluation" = c("overview", "design_summary"),
    "summary.mfrm_signal_detection" = c("overview", "detection_summary"),
    "summary.mfrm_population_prediction" = c("overview", "design", "forecast"),
    "summary.mfrm_future_branch_active_branch" = c("overview", "profile_summary", "recommendation_table"),
    "summary.mfrm_facets_run" = c("overview", "mapping", "run_info", "fit", "diagnostics"),
    "summary.mfrm_bias" = c("overview", "top_rows"),
    "summary.mfrm_anchor_audit" = c("facet_summary", "recommendations"),
    "summary.mfrm_linking_review" = c("overview", "top_linking_risks", "group_view_index", "reporting_map"),
    "summary.mfrm_misfit_casebook" = c("overview", "top_cases", "case_rollup", "group_view_index", "reporting_map"),
    "summary.mfrm_weighting_audit" = c("overview", "top_reweighted_levels", "reporting_map"),
    "summary.mfrm_unit_prediction" = c("estimates", "settings"),
    "summary.mfrm_plausible_values" = c("draw_summary", "settings"),
    character(0)
  )
}

validate_summary_table_bundle_summary <- function(summary_obj,
                                                 summary_class,
                                                 helper = "build_summary_table_bundle()") {
  if (!is.list(summary_obj)) {
    stop(
      "`", helper, "` requires a supported package object or a package-native `summary()` output. ",
      "The supplied summary object for class `", as.character(summary_class %||% "unknown"),
      "` is not a list and does not match the package summary contract.",
      call. = FALSE
    )
  }

  required <- summary_table_bundle_required_components(summary_class)
  if (length(required) == 0L) {
    return(invisible(summary_obj))
  }

  missing_components <- required[
    !vapply(required, function(nm) {
      nm %in% names(summary_obj) && !is.null(summary_obj[[nm]])
    }, logical(1))
  ]

  if (length(missing_components) > 0L) {
    stop(
      "`", helper, "` received a malformed `", as.character(summary_class),
      "` object. Missing required component(s): ",
      paste(missing_components, collapse = ", "),
      ". Rebuild the source object with the package helper, then call `summary()` again.",
      call. = FALSE
    )
  }

  invisible(summary_obj)
}

validate_summary_table_bundle_inputs <- function(x,
                                                 which = NULL,
                                                 appendix_preset = NULL,
                                                 include_empty = FALSE,
                                                 digits = 3,
                                                 top_n = 10,
                                                 preview_chars = 160,
                                                 helper = "build_summary_table_bundle()") {
  if (missing(x) || is.null(x)) {
    stop(
      "`", helper, "` requires `x` to be a supported package object or one of its `summary()` outputs.",
      call. = FALSE
    )
  }

  if (!is.null(which)) {
    if (!is.character(which) || length(which) == 0L) {
      stop(
        "`", helper, "` requires `which` to be `NULL` or a non-empty character vector of table names.",
        call. = FALSE
      )
    }
    which <- trimws(which)
    if (anyNA(which) || any(!nzchar(which))) {
      stop(
        "`", helper, "` requires every `which` entry to be a non-empty table name.",
        call. = FALSE
      )
    }
    which <- unique(which)
  }

  if (!is.logical(include_empty) || length(include_empty) != 1L || is.na(include_empty)) {
    stop(
      "`", helper, "` requires `include_empty` to be either `TRUE` or `FALSE`.",
      call. = FALSE
    )
  }

  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) || digits < 0) {
    stop(
      "`", helper, "` requires `digits` to be a single non-negative number.",
      call. = FALSE
    )
  }
  digits <- as.integer(digits)

  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
    stop(
      "`", helper, "` requires `top_n` to be a single positive number.",
      call. = FALSE
    )
  }
  top_n <- as.integer(top_n)

  if (!is.numeric(preview_chars) || length(preview_chars) != 1L || !is.finite(preview_chars) || preview_chars < 1) {
    stop(
      "`", helper, "` requires `preview_chars` to be a single positive number.",
      call. = FALSE
    )
  }
  preview_chars <- as.integer(preview_chars)

  if (!is.null(appendix_preset)) {
    if (!is.character(appendix_preset) || length(appendix_preset) != 1L ||
        is.na(appendix_preset) || !nzchar(trimws(appendix_preset))) {
      stop(
        "`", helper, "` requires `appendix_preset` to be `NULL` or a single preset name.",
        call. = FALSE
      )
    }
    appendix_preset <- match.arg(
      tolower(trimws(as.character(appendix_preset[1]))),
      c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")
    )
  }

  resolved <- resolve_summary_table_bundle_input(
    x,
    digits = digits,
    top_n = top_n,
    preview_chars = preview_chars
  )
  validate_summary_table_bundle_summary(
    summary_obj = resolved$summary,
    summary_class = resolved$summary_class,
    helper = helper
  )

  list(
    resolved = resolved,
    which = which,
    appendix_preset = appendix_preset,
    include_empty = include_empty,
    digits = digits,
    top_n = top_n,
    preview_chars = preview_chars
  )
}

summary_table_bundle_resolve_future_branch_summary <- function(summary_obj) {
  if (inherits(summary_obj, "summary.mfrm_future_branch_active_branch")) {
    return(summary_obj)
  }
  summary_obj$future_branch_active_summary %||% NULL
}

summary_table_bundle_future_branch_spec <- function(summary_obj,
                                                    embedded = TRUE) {
  future <- summary_table_bundle_resolve_future_branch_summary(summary_obj)
  if (!inherits(future, "summary.mfrm_future_branch_active_branch")) {
    future <- NULL
  }
  overview_desc <- if (isTRUE(embedded)) {
    "Deterministic overview of the embedded future arbitrary-facet planning scaffold."
  } else {
    "Deterministic overview of the future arbitrary-facet planning active branch."
  }
  profile_desc <- if (isTRUE(embedded)) {
    "Exact-count and balanced-expectation design metrics from the embedded future-branch scaffold."
  } else {
    "Exact-count and balanced-expectation design metrics from the future arbitrary-facet planning active branch."
  }
  load_balance_desc <- if (isTRUE(embedded)) {
    "Deterministic rater-load and integer-balance diagnostics from the embedded future-branch scaffold."
  } else {
    "Deterministic rater-load and integer-balance diagnostics from the future arbitrary-facet planning active branch."
  }
  coverage_desc <- if (isTRUE(embedded)) {
    "Deterministic coverage and connectivity summaries from the embedded future-branch scaffold."
  } else {
    "Deterministic coverage and connectivity summaries from the future arbitrary-facet planning active branch."
  }
  guardrail_desc <- if (isTRUE(embedded)) {
    "Exact structural guardrail classifications from the embedded future-branch scaffold."
  } else {
    "Exact structural guardrail classifications from the future arbitrary-facet planning active branch."
  }
  readiness_desc <- if (isTRUE(embedded)) {
    "Structural readiness tiers indicating which overlap/balance conditions currently hold."
  } else {
    "Structural readiness tiers for the future arbitrary-facet planning active branch."
  }
  recommendation_desc <- if (isTRUE(embedded)) {
    "Conservative structural recommendation derived from the embedded future-branch scaffold."
  } else {
    "Conservative structural recommendation derived from the future arbitrary-facet planning active branch."
  }
  list(
    tables = list(
      future_branch_overview = summary_table_bundle_df(future$overview),
      future_branch_profile = summary_table_bundle_df(future$profile_summary),
      future_branch_load_balance = summary_table_bundle_df(future$load_balance_summary),
      future_branch_coverage = summary_table_bundle_df(future$coverage_summary),
      future_branch_guardrails = summary_table_bundle_df(future$guardrail_summary),
      future_branch_readiness = summary_table_bundle_df(future$readiness_summary),
      future_branch_recommendation = summary_table_bundle_df(future$recommendation_table),
      future_branch_appendix_presets = summary_table_bundle_df(future$appendix_presets),
      future_branch_appendix_roles = summary_table_bundle_df(future$appendix_role_summary),
      future_branch_appendix_sections = summary_table_bundle_df(future$appendix_section_summary),
      future_branch_selection_table_presets = summary_table_bundle_df(future$selection_table_preset_summary),
      future_branch_selection_handoff_tables = summary_table_bundle_df(future$selection_handoff_table_summary),
      future_branch_selection_handoff_presets = summary_table_bundle_df(future$selection_handoff_preset_summary),
      future_branch_selection_handoff = summary_table_bundle_df(future$selection_handoff_summary),
      future_branch_selection_handoff_bundles = summary_table_bundle_df(future$selection_handoff_bundle_summary),
      future_branch_selection_handoff_roles = summary_table_bundle_df(future$selection_handoff_role_summary),
      future_branch_selection_handoff_role_sections = summary_table_bundle_df(future$selection_handoff_role_section_summary),
      future_branch_selection_tables = summary_table_bundle_df(future$selection_table_summary),
      future_branch_selection_summary = summary_table_bundle_df(future$selection_summary),
      future_branch_selection_roles = summary_table_bundle_df(future$selection_role_summary),
      future_branch_selection_sections = summary_table_bundle_df(future$selection_section_summary),
      future_branch_selection_catalog = summary_table_bundle_df(future$selection_catalog),
      future_branch_reporting_map = summary_table_bundle_df(future$reporting_map)
    ),
    roles = c(
      future_branch_overview = "future_branch_overview",
      future_branch_profile = "future_branch_profile",
      future_branch_load_balance = "future_branch_load_balance",
      future_branch_coverage = "future_branch_coverage",
      future_branch_guardrails = "future_branch_guardrails",
      future_branch_readiness = "future_branch_readiness",
      future_branch_recommendation = "future_branch_recommendation",
      future_branch_appendix_presets = "future_branch_appendix_presets",
      future_branch_appendix_roles = "future_branch_appendix_roles",
      future_branch_appendix_sections = "future_branch_appendix_sections",
      future_branch_selection_table_presets = "future_branch_selection_table_presets",
      future_branch_selection_handoff_tables = "future_branch_selection_handoff_tables",
      future_branch_selection_handoff_presets = "future_branch_selection_handoff_presets",
      future_branch_selection_handoff = "future_branch_selection_handoff",
      future_branch_selection_handoff_bundles = "future_branch_selection_handoff_bundles",
      future_branch_selection_handoff_roles = "future_branch_selection_handoff_roles",
      future_branch_selection_handoff_role_sections = "future_branch_selection_handoff_role_sections",
      future_branch_selection_tables = "future_branch_selection_tables",
      future_branch_selection_summary = "future_branch_selection_summary",
      future_branch_selection_roles = "future_branch_selection_roles",
      future_branch_selection_sections = "future_branch_selection_sections",
      future_branch_selection_catalog = "future_branch_selection_catalog",
      future_branch_reporting_map = "future_branch_reporting_map"
    ),
    descriptions = c(
      future_branch_overview = overview_desc,
      future_branch_profile = profile_desc,
      future_branch_load_balance = load_balance_desc,
      future_branch_coverage = coverage_desc,
      future_branch_guardrails = guardrail_desc,
      future_branch_readiness = readiness_desc,
      future_branch_recommendation = recommendation_desc,
      future_branch_appendix_presets = "Preset-level appendix routing counts for the future arbitrary-facet planning surface.",
      future_branch_appendix_roles = "Appendix routing counts by reporting role for the future arbitrary-facet planning surface.",
      future_branch_appendix_sections = "Appendix routing counts by manuscript section for the future arbitrary-facet planning surface.",
      future_branch_selection_table_presets = "Preset-specific appendix table selections for the future arbitrary-facet planning surface.",
      future_branch_selection_handoff_tables = "Preset-specific table-level appendix handoff crosswalk for the future arbitrary-facet planning surface.",
      future_branch_selection_handoff_presets = "Preset-level appendix handoff overview for the future arbitrary-facet planning surface.",
      future_branch_selection_handoff = "Section-aware appendix handoff summary for the future arbitrary-facet planning surface.",
      future_branch_selection_handoff_bundles = "Bundle-aware appendix handoff summary for the future arbitrary-facet planning surface.",
      future_branch_selection_handoff_roles = "Role-aware appendix handoff summary for the future arbitrary-facet planning surface.",
      future_branch_selection_handoff_role_sections = "Role-by-section appendix handoff summary for the future arbitrary-facet planning surface.",
      future_branch_selection_tables = "Preset-aware appendix table selections for the future arbitrary-facet planning surface.",
      future_branch_selection_summary = "Preset-filtered appendix selection counts for the future arbitrary-facet planning surface.",
      future_branch_selection_roles = "Preset-filtered appendix selection counts by reporting role for the future arbitrary-facet planning surface.",
      future_branch_selection_sections = "Preset-filtered appendix selection counts by manuscript section for the future arbitrary-facet planning surface.",
      future_branch_selection_catalog = "Full appendix selection catalog for the future arbitrary-facet planning surface.",
      future_branch_reporting_map = "Direct reporting-map bridge for the future arbitrary-facet planning surface."
    )
  )
}

summary_table_bundle_spec <- function(summary_obj) {
  cls <- class(summary_obj)[1]
  switch(
    cls,
    "summary.mfrm_fit" = list(
      title = "Model Summary Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        population_overview = summary_table_bundle_df(summary_obj$population_overview),
        population_design = summary_table_bundle_df(summary_obj$population_design),
        population_coefficients = summary_table_bundle_df(summary_obj$population_coefficients),
        population_coding = summary_table_bundle_df(summary_obj$population_coding),
        facet_overview = summary_table_bundle_df(summary_obj$facet_overview),
        person_overview = summary_table_bundle_df(summary_obj$person_overview),
        step_overview = summary_table_bundle_df(summary_obj$step_overview),
        slope_overview = summary_table_bundle_df(summary_obj$slope_overview),
        settings_overview = summary_table_bundle_df(summary_obj$settings_overview),
        reporting_map = summary_table_bundle_df(summary_obj$reporting_map),
        caveats = summary_table_bundle_df(summary_obj$caveats),
        facet_extremes = summary_table_bundle_df(summary_obj$facet_extremes),
        person_high = summary_table_bundle_df(summary_obj$person_high),
        person_low = summary_table_bundle_df(summary_obj$person_low)
      ),
      roles = c(
        overview = "run_overview",
        population_overview = "population_basis",
        population_design = "population_design",
        population_coefficients = "population_coefficients",
        population_coding = "population_coding",
        facet_overview = "facet_distribution",
        person_overview = "person_distribution",
        step_overview = "category_structure",
        slope_overview = "gpcm_discrimination",
        settings_overview = "estimation_settings",
        reporting_map = "reporting_map",
        caveats = "analysis_caveats",
        facet_extremes = "extreme_facet_levels",
        person_high = "extreme_person_high",
        person_low = "extreme_person_low"
      ),
      descriptions = c(
        overview = "One-row model fit, convergence, and information-criteria overview.",
        population_overview = "Population-model basis, posterior basis, and omission audit.",
        population_design = "Population-model design-matrix columns and numeric audit statistics.",
        population_coefficients = "Latent-regression coefficients when the population model is active.",
        population_coding = "Latent-regression categorical covariate levels, contrasts, and encoded model-matrix columns.",
        facet_overview = "Per-facet spread, range, and level-count summary.",
        person_overview = "Distribution of person measures and posterior SD summaries.",
        step_overview = "Threshold range and monotonicity summary.",
        slope_overview = "GPCM discrimination summary under the current identification.",
        settings_overview = "Estimation settings that affect identification and interpretation.",
        reporting_map = "Companion outputs to cite for manuscript-oriented reporting.",
        caveats = "Structured fit-level caveats such as retained zero-count categories, score-category recoding, and latent-regression population-model warnings.",
        facet_extremes = "Facet levels with the largest absolute estimates.",
        person_high = "Highest person measures from the current fit.",
        person_low = "Lowest person measures from the current fit."
      )
    ),
    "summary.mfrm_diagnostics" = list(
      title = "Diagnostics Summary Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        overall_fit = summary_table_bundle_df(summary_obj$overall_fit),
        precision_profile = summary_table_bundle_df(summary_obj$precision_profile),
        precision_audit = summary_table_bundle_df(summary_obj$precision_audit),
        reliability = summary_table_bundle_df(summary_obj$reliability),
        top_fit = summary_table_bundle_df(summary_obj$top_fit),
        reporting_map = summary_table_bundle_df(summary_obj$reporting_map),
        flags = summary_table_bundle_df(summary_obj$flags)
      ),
      roles = c(
        overview = "run_overview",
        overall_fit = "overall_fit",
        precision_profile = "precision_basis",
        precision_audit = "precision_audit",
        reliability = "facet_precision",
        top_fit = "extreme_fit_rows",
        reporting_map = "reporting_map",
        flags = "flag_counts"
      ),
      descriptions = c(
        overview = "Run-level diagnostic coverage and precision tier.",
        overall_fit = "Global fit statistics from the current diagnostic run.",
        precision_profile = "Precision basis and recommended interpretation tier.",
        precision_audit = "Precision checks marked review/warn for manuscript caution.",
        reliability = "Facet-level separation, strata, and reliability summary.",
        top_fit = "Rows with the largest absolute fit Z statistics.",
        reporting_map = "Companion outputs for manuscript reporting beyond summary(diag).",
        flags = "Counts of unexpected responses, displacement, interactions, and inter-rater pairs."
      )
    ),
    "summary.mfrm_data_description" = list(
      title = "Data Description Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        missing = summary_table_bundle_df(summary_obj$missing),
        score_distribution = summary_table_bundle_df(summary_obj$score_distribution),
        facet_overview = summary_table_bundle_df(summary_obj$facet_overview),
        agreement = summary_table_bundle_df(summary_obj$agreement),
        reporting_map = summary_table_bundle_df(summary_obj$reporting_map),
        caveats = summary_table_bundle_df(summary_obj$caveats)
      ),
      roles = c(
        overview = "run_overview",
        missing = "missingness",
        score_distribution = "score_usage",
        facet_overview = "facet_coverage",
        agreement = "agreement",
        reporting_map = "reporting_map",
        caveats = "score_category_caveats"
      ),
      descriptions = c(
        overview = "One-row sample, design, and rating-span overview.",
        missing = "Missing-value counts by selected input column.",
        score_distribution = "Observed score distribution for category-usage reporting.",
        facet_overview = "Facet-level coverage and weighted counts.",
        agreement = "Observed inter-rater agreement summary when available.",
        reporting_map = "Companion outputs for fit, reliability, and residual follow-up.",
        caveats = "Structured pre-fit score-support caveats such as retained zero-count categories."
      )
    ),
    "summary.mfrm_reporting_checklist" = list(
      title = "Reporting Checklist Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        section_summary = summary_table_bundle_df(summary_obj$section_summary),
        priority_summary = summary_table_bundle_df(summary_obj$priority_summary),
        action_items = summary_table_bundle_df(summary_obj$action_items),
        settings = summary_table_bundle_df(summary_obj$settings)
      ),
      roles = c(
        overview = "checklist_overview",
        section_summary = "section_coverage",
        priority_summary = "priority_distribution",
        action_items = "draft_actions",
        settings = "checklist_settings"
      ),
      descriptions = c(
        overview = "Overall checklist coverage across sections and draft-readiness flags.",
        section_summary = "Coverage summary by reporting section.",
        priority_summary = "High/medium/low/ready counts by severity.",
        action_items = "Top unresolved manuscript-drafting actions.",
        settings = "Checklist settings used to build the reporting contract."
      )
    ),
    "summary.mfrm_apa_outputs" = list(
      title = "APA Output Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        components = summary_table_bundle_df(summary_obj$components),
        sections = summary_table_bundle_df(summary_obj$sections),
        content_checks = summary_table_bundle_df(summary_obj$content_checks),
        preview = summary_table_bundle_df(summary_obj$preview)
      ),
      roles = c(
        overview = "draft_overview",
        components = "component_stats",
        sections = "section_coverage",
        content_checks = "draft_checks",
        preview = "text_preview"
      ),
      descriptions = c(
        overview = "Overall coverage for manuscript draft text products.",
        components = "Per-component line, character, and mention counts.",
        sections = "Availability of the package-native section map.",
        content_checks = "Contract-based checks for APA drafting completeness.",
        preview = "Compact preview of the first non-empty lines in each draft component."
      )
    ),
    "summary.mfrm_design_evaluation" = {
      future_spec <- summary_table_bundle_future_branch_spec(summary_obj)
      list(
        title = "Design Evaluation Tables",
        tables = c(
          list(
            overview = summary_table_bundle_df(summary_obj$overview),
            design_summary = summary_table_bundle_df(summary_obj$design_summary)
          ),
          future_spec$tables
        ),
        roles = c(
          overview = "run_overview",
          design_summary = "design_performance",
          future_spec$roles
        ),
        descriptions = c(
          overview = "Run-level overview for the current design-evaluation study.",
          design_summary = "Aggregated Monte Carlo design summaries for the active two-role planner.",
          future_spec$descriptions
        )
      )
    },
    "summary.mfrm_signal_detection" = {
      future_spec <- summary_table_bundle_future_branch_spec(summary_obj)
      list(
        title = "Signal Detection Tables",
        tables = c(
          list(
            overview = summary_table_bundle_df(summary_obj$overview),
            detection_summary = summary_table_bundle_df(summary_obj$detection_summary)
          ),
          future_spec$tables
        ),
        roles = c(
          overview = "run_overview",
          detection_summary = "signal_detection",
          future_spec$roles
        ),
        descriptions = c(
          overview = "Run-level overview for the current signal-detection study.",
          detection_summary = "Aggregated DIF/bias screening summaries for the active two-role planner.",
          future_spec$descriptions
        )
      )
    },
    "summary.mfrm_population_prediction" = {
      future_spec <- summary_table_bundle_future_branch_spec(summary_obj)
      list(
        title = "Population Prediction Tables",
        tables = c(
          list(
            design = summary_table_bundle_df(summary_obj$design),
            overview = summary_table_bundle_df(summary_obj$overview),
            forecast = summary_table_bundle_df(summary_obj$forecast)
          ),
          future_spec$tables
        ),
        roles = c(
          design = "design_grid",
          overview = "run_overview",
          forecast = "forecast_summary",
          future_spec$roles
        ),
        descriptions = c(
          design = "Requested future design grid used for the current forecast run.",
          overview = "Run-level overview for the current population forecast.",
          forecast = "Facet-level forecast summary for the active two-role planner.",
          future_spec$descriptions
        )
      )
    },
    "summary.mfrm_future_branch_active_branch" = {
      future_spec <- summary_table_bundle_future_branch_spec(
        summary_obj,
        embedded = FALSE
      )
      list(
        title = "Future Arbitrary-Facet Planning Tables",
        tables = future_spec$tables,
        roles = future_spec$roles,
        descriptions = future_spec$descriptions
      )
    },
    "summary.mfrm_facets_run" = list(
      title = "Workflow Summary Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        mapping = summary_table_bundle_df(summary_obj$mapping),
        run_info = summary_table_bundle_df(summary_obj$run_info),
        fit_overview = summary_table_bundle_df(summary_obj$fit$overview),
        fit_reporting_map = summary_table_bundle_df(summary_obj$fit$reporting_map),
        diagnostic_overview = summary_table_bundle_df(summary_obj$diagnostics$overview),
        diagnostic_flags = summary_table_bundle_df(summary_obj$diagnostics$flags),
        diagnostic_reporting_map = summary_table_bundle_df(summary_obj$diagnostics$reporting_map)
      ),
      roles = c(
        overview = "workflow_overview",
        mapping = "column_mapping",
        run_info = "workflow_settings",
        fit_overview = "run_overview",
        fit_reporting_map = "reporting_map",
        diagnostic_overview = "run_overview",
        diagnostic_flags = "flag_counts",
        diagnostic_reporting_map = "reporting_map"
      ),
      descriptions = c(
        overview = "Legacy-compatible workflow overview with fit metadata.",
        mapping = "Resolved column mapping for the one-shot workflow run.",
        run_info = "Workflow settings and pipeline metadata recorded by run_mfrm_facets().",
        fit_overview = "Nested model-fit overview routed from summary(out$fit).",
        fit_reporting_map = "Nested reporting-map follow-up routed from summary(out$fit).",
        diagnostic_overview = "Nested diagnostic overview routed from summary(out$diagnostics).",
        diagnostic_flags = "Nested diagnostic flag counts routed from summary(out$diagnostics).",
        diagnostic_reporting_map = "Nested reporting-map follow-up routed from summary(out$diagnostics)."
      )
    ),
    "summary.mfrm_bias" = list(
      title = "Bias Summary Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        chi_sq = summary_table_bundle_df(summary_obj$chi_sq),
        final_iteration = summary_table_bundle_df(summary_obj$final_iteration),
        top_rows = summary_table_bundle_df(summary_obj$top_rows),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "bias_overview",
        chi_sq = "bias_chi_square",
        final_iteration = "bias_iteration_status",
        top_rows = "bias_screening_rows",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Interaction-order overview and screening counts for the current bias run.",
        chi_sq = "Fixed-effect chi-square block from the current bias run.",
        final_iteration = "Final bias-iteration status row for stabilization checks.",
        top_rows = "Highest-|t| interaction rows for immediate follow-up.",
        notes = "Compact interpretation notes for screening-oriented bias reporting."
      )
    ),
    "summary.mfrm_anchor_audit" = {
      issue_tbl <- summary_table_bundle_df(summary_obj$issue_counts)
      facet_tbl <- summary_table_bundle_df(summary_obj$facet_summary)
      level_tbl <- summary_table_bundle_df(summary_obj$level_observation_summary)
      category_tbl <- summary_table_bundle_df(summary_obj$category_counts)
      rec_tbl <- summary_table_bundle_text_df(summary_obj$recommendations, column = "Recommendation")
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      overview_tbl <- data.frame(
        IssueRows = nrow(issue_tbl),
        Facets = nrow(facet_tbl),
        LevelRows = nrow(level_tbl),
        CategoryRows = nrow(category_tbl),
        Recommendations = nrow(rec_tbl),
        stringsAsFactors = FALSE
      )
      list(
        title = "Anchor Audit Tables",
        tables = list(
          overview = overview_tbl,
          issue_counts = issue_tbl,
          facet_summary = facet_tbl,
          level_observation_summary = level_tbl,
          category_counts = category_tbl,
          recommendations = rec_tbl,
          notes = notes_tbl
        ),
        roles = c(
          overview = "anchor_audit_overview",
          issue_counts = "anchor_issue_counts",
          facet_summary = "facet_coverage",
          level_observation_summary = "level_observation_audit",
          category_counts = "category_usage",
          recommendations = "repair_recommendations",
          notes = "interpretation_notes"
        ),
        descriptions = c(
          overview = "Anchor-audit overview with issue, facet, and recommendation counts.",
          issue_counts = "Observed anchor-audit issue counts ranked by frequency.",
          facet_summary = "Facet-level counts and anchor-table coverage summary.",
          level_observation_summary = "Observation counts by facet level for anchor viability checks.",
          category_counts = "Observed score-category usage for anchor-audit screening.",
          recommendations = "Compact action list for anchor repair or review.",
          notes = "One-line interpretation note from the anchor audit."
        )
      )
    },
    "summary.mfrm_linking_review" = {
      overview_tbl <- summary_table_bundle_df(summary_obj$overview)
      status_tbl <- summary_table_bundle_df(summary_obj$status)
      top_tbl <- summary_table_bundle_df(summary_obj$top_linking_risks)
      group_view_index_tbl <- summary_table_bundle_df(summary_obj$group_view_index)
      prefit_tbl <- summary_table_bundle_df(summary_obj$prefit_anchor_risks)
      drift_tbl <- summary_table_bundle_df(summary_obj$drift_risks)
      chain_tbl <- summary_table_bundle_df(summary_obj$chain_risks)
      plot_tbl <- summary_table_bundle_df(summary_obj$plot_map)
      reporting_tbl <- summary_table_bundle_df(summary_obj$reporting_map)
      support_tbl <- summary_table_bundle_df(summary_obj$support_status)
      actions_tbl <- summary_table_bundle_text_df(summary_obj$next_actions, column = "Action")
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      list(
        title = "Linking Review Tables",
        tables = list(
          overview = overview_tbl,
          status = status_tbl,
          top_linking_risks = top_tbl,
          group_view_index = group_view_index_tbl,
          prefit_anchor_risks = prefit_tbl,
          drift_risks = drift_tbl,
          chain_risks = chain_tbl,
          plot_map = plot_tbl,
          reporting_map = reporting_tbl,
          support_status = support_tbl,
          next_actions = actions_tbl,
          notes = notes_tbl,
          settings = settings_tbl
        ),
        roles = c(
          overview = "linking_review_overview",
          status = "review_status",
          top_linking_risks = "linking_risk_screen",
          group_view_index = "linking_risk_group_index",
          prefit_anchor_risks = "prefit_anchor_risks",
          drift_risks = "drift_risks",
          chain_risks = "chain_risks",
          plot_map = "plot_routing",
          reporting_map = "reporting_map",
          support_status = "capability_boundary",
          next_actions = "repair_recommendations",
          notes = "interpretation_notes",
          settings = "review_settings"
        ),
        descriptions = c(
          overview = "Overview of evidence sources and current operational linking status.",
          status = "Compact front-door status block for linking review.",
          top_linking_risks = "Highest-priority linking risks across anchor, drift, and chain evidence.",
          group_view_index = "Index of stable wave/link/facet/source-family grouping views available for operational linking triage.",
          prefit_anchor_risks = "Pre-fit anchor adequacy issues and overlap-support warnings.",
          drift_risks = "Wave-level drift and retained-common-element support warnings.",
          chain_risks = "Adjacent-link instability rows from the screened equating chain.",
          plot_map = "Routing map to existing plotting helpers for operational follow-up.",
          reporting_map = "Map from operational review outputs to manuscript/reporting companions.",
          support_status = "Current support contract for RSM/PCM versus bounded GPCM use.",
          next_actions = "Top next-step actions for anchor repair or linking follow-up.",
          notes = "Compact interpretation notes for operational linking review.",
          settings = "Settings and provenance recorded by build_linking_review()."
        )
      )
    },
    "summary.mfrm_misfit_casebook" = {
      overview_tbl <- summary_table_bundle_df(summary_obj$overview)
      status_tbl <- summary_table_bundle_df(summary_obj$status)
      top_cases_tbl <- summary_table_bundle_df(summary_obj$top_cases)
      case_rollup_tbl <- summary_table_bundle_df(summary_obj$case_rollup)
      group_view_index_tbl <- summary_table_bundle_df(summary_obj$group_view_index)
      source_summary_tbl <- summary_table_bundle_df(summary_obj$source_summary)
      plot_tbl <- summary_table_bundle_df(summary_obj$plot_map)
      reporting_tbl <- summary_table_bundle_df(summary_obj$reporting_map)
      support_tbl <- summary_table_bundle_df(summary_obj$support_status)
      warning_tbl <- summary_table_bundle_text_df(summary_obj$key_warnings, column = "Warning")
      actions_tbl <- summary_table_bundle_text_df(summary_obj$next_actions, column = "Action")
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      list(
        title = "Misfit Casebook Tables",
        tables = list(
          overview = overview_tbl,
          status = status_tbl,
          top_cases = top_cases_tbl,
          case_rollup = case_rollup_tbl,
          group_view_index = group_view_index_tbl,
          source_summary = source_summary_tbl,
          plot_map = plot_tbl,
          reporting_map = reporting_tbl,
          support_status = support_tbl,
          key_warnings = warning_tbl,
          next_actions = actions_tbl,
          notes = notes_tbl,
          settings = settings_tbl
        ),
        roles = c(
          overview = "misfit_casebook_overview",
          status = "review_status",
          top_cases = "misfit_case_rows",
          case_rollup = "misfit_case_rollup",
          group_view_index = "misfit_case_rollup",
          source_summary = "misfit_case_sources",
          plot_map = "plot_routing",
          reporting_map = "reporting_map",
          support_status = "capability_boundary",
          key_warnings = "review_status",
          next_actions = "repair_recommendations",
          notes = "interpretation_notes",
          settings = "review_settings"
        ),
        descriptions = c(
          overview = "Overview of the current operational misfit case-review queue.",
          status = "Compact front-door status block for the misfit casebook.",
          top_cases = "Highest-priority case rows preserved by source family without collapsing evidence into one opaque score.",
          case_rollup = "Secondary grouping view that summarizes where flagged cases concentrate by person, facet level, pair, or source family.",
          group_view_index = "Index of stable grouping views available for operational triage on top of the raw case rows.",
          source_summary = "Counts and maximum priority by source family for the current casebook.",
          plot_map = "Routing map from casebook source families to dedicated follow-up plotting helpers.",
          reporting_map = "Map from operational case review to reporting and appendix companions.",
          support_status = "Current support contract for Rasch-family versus bounded GPCM case review.",
          key_warnings = "Top warning lines for the current casebook build.",
          next_actions = "Top next-step actions for misfit case follow-up.",
          notes = "Compact interpretation notes for the misfit casebook.",
          settings = "Casebook settings and source-family provenance."
        )
      )
    },
    "summary.mfrm_weighting_audit" = {
      overview_tbl <- summary_table_bundle_df(summary_obj$overview)
      status_tbl <- summary_table_bundle_df(summary_obj$status)
      top_shift_tbl <- summary_table_bundle_df(summary_obj$top_measure_shifts)
      top_reweighted_tbl <- summary_table_bundle_df(summary_obj$top_reweighted_levels)
      plot_tbl <- summary_table_bundle_df(summary_obj$plot_map)
      reporting_tbl <- summary_table_bundle_df(summary_obj$reporting_map)
      support_tbl <- summary_table_bundle_df(summary_obj$support_status)
      warning_tbl <- summary_table_bundle_text_df(summary_obj$key_warnings, column = "Warning")
      actions_tbl <- summary_table_bundle_text_df(summary_obj$next_actions, column = "Action")
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      list(
        title = "Weighting Audit Tables",
        tables = list(
          overview = overview_tbl,
          status = status_tbl,
          top_measure_shifts = top_shift_tbl,
          top_reweighted_levels = top_reweighted_tbl,
          plot_map = plot_tbl,
          reporting_map = reporting_tbl,
          support_status = support_tbl,
          key_warnings = warning_tbl,
          next_actions = actions_tbl,
          notes = notes_tbl,
          settings = settings_tbl
        ),
        roles = c(
          overview = "weighting_review_overview",
          status = "review_status",
          top_measure_shifts = "reweighting_measure_shift",
          top_reweighted_levels = "gpcm_discrimination",
          plot_map = "plot_routing",
          reporting_map = "reporting_map",
          support_status = "capability_boundary",
          key_warnings = "review_status",
          next_actions = "repair_recommendations",
          notes = "interpretation_notes",
          settings = "estimation_settings"
        ),
        descriptions = c(
          overview = "Overview of the equal-weighting versus bounded GPCM weighting audit.",
          status = "Compact status block for the weighting-policy review.",
          top_measure_shifts = "Largest non-person facet-measure shifts between the Rasch-family reference and bounded GPCM.",
          top_reweighted_levels = "Largest slope-facet reweighting signals under bounded GPCM.",
          plot_map = "Public plot routes for precision redistribution and comparison follow-up.",
          reporting_map = "Bundle/report handoff map for weighting-policy review outputs.",
          support_status = "Capability-boundary statement for the bounded GPCM weighting audit.",
          key_warnings = "Top warning lines for weighting-policy review.",
          next_actions = "Recommended next-step actions after weighting-policy review.",
          notes = "Interpretation notes for the weighting audit.",
          settings = "Weighting-audit settings and theta-grid parameters."
        )
      )
    },
    "summary.mfrm_unit_prediction" = {
      estimate_tbl <- summary_table_bundle_df(summary_obj$estimates)
      audit_tbl <- summary_table_bundle_df(summary_obj$audit)
      pop_audit_tbl <- summary_table_bundle_df(summary_obj$population_audit)
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      overview_tbl <- data.frame(
        Units = nrow(estimate_tbl),
        AuditRows = nrow(audit_tbl),
        PopulationAuditRows = nrow(pop_audit_tbl),
        Settings = nrow(settings_tbl),
        Notes = nrow(notes_tbl),
        stringsAsFactors = FALSE
      )
      list(
        title = "Unit Prediction Tables",
        tables = list(
          overview = overview_tbl,
          estimates = estimate_tbl,
          audit = audit_tbl,
          population_audit = pop_audit_tbl,
          settings = settings_tbl,
          notes = notes_tbl
        ),
        roles = c(
          overview = "prediction_overview",
          estimates = "unit_estimates",
          audit = "prediction_audit",
          population_audit = "population_audit",
          settings = "scoring_settings",
          notes = "interpretation_notes"
        ),
        descriptions = c(
          overview = "Posterior-scoring overview for the current unit-prediction run.",
          estimates = "Posterior summaries for the scored persons.",
          audit = "Row-level preparation audit for the supplied scoring data.",
          population_audit = "Optional person-level omission audit for latent-regression scoring.",
          settings = "Scoring settings carried into posterior unit prediction.",
          notes = "Compact interpretation notes for posterior scoring output."
        )
      )
    },
    "summary.mfrm_plausible_values" = {
      draw_tbl <- summary_table_bundle_df(summary_obj$draw_summary)
      estimate_tbl <- summary_table_bundle_df(summary_obj$estimates)
      audit_tbl <- summary_table_bundle_df(summary_obj$audit)
      pop_audit_tbl <- summary_table_bundle_df(summary_obj$population_audit)
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      total_draws <- if ("Draws" %in% names(draw_tbl)) sum(draw_tbl$Draws, na.rm = TRUE) else nrow(draw_tbl)
      overview_tbl <- data.frame(
        Persons = nrow(draw_tbl),
        TotalDraws = total_draws,
        EstimateRows = nrow(estimate_tbl),
        AuditRows = nrow(audit_tbl),
        PopulationAuditRows = nrow(pop_audit_tbl),
        Settings = nrow(settings_tbl),
        Notes = nrow(notes_tbl),
        stringsAsFactors = FALSE
      )
      list(
        title = "Plausible Value Tables",
        tables = list(
          overview = overview_tbl,
          draw_summary = draw_tbl,
          estimates = estimate_tbl,
          audit = audit_tbl,
          population_audit = pop_audit_tbl,
          settings = settings_tbl,
          notes = notes_tbl
        ),
        roles = c(
          overview = "plausible_value_overview",
          draw_summary = "plausible_value_draws",
          estimates = "unit_estimates",
          audit = "prediction_audit",
          population_audit = "population_audit",
          settings = "scoring_settings",
          notes = "interpretation_notes"
        ),
        descriptions = c(
          overview = "Approximate plausible-value overview for the current posterior scoring run.",
          draw_summary = "Empirical summaries of the sampled posterior draws by person.",
          estimates = "Companion posterior EAP summaries paired with the draw summary.",
          audit = "Row-level preparation audit for the supplied scoring data.",
          population_audit = "Optional person-level omission audit for latent-regression scoring.",
          settings = "Scoring settings used to generate the approximate plausible values.",
          notes = "Compact interpretation notes for plausible-value reporting."
        )
      )
    },
    stop("Unsupported summary class for table-bundle conversion: ", cls, call. = FALSE)
  )
}

build_summary_table_index <- function(tables, roles, descriptions) {
  do.call(
    rbind,
    lapply(names(tables), function(nm) {
      tbl <- tables[[nm]]
      data.frame(
        Table = nm,
        Rows = nrow(tbl),
        Cols = ncol(tbl),
        Role = as.character(roles[[nm]] %||% ""),
        Description = as.character(descriptions[[nm]] %||% ""),
        stringsAsFactors = FALSE
      )
    })
  )
}

#' Build a manuscript-oriented table bundle from `summary()` outputs
#'
#' @param x An `mfrm_fit`, `mfrm_diagnostics`, `mfrm_data_description`,
#'   `mfrm_reporting_checklist`, `mfrm_apa_outputs`,
#'   `mfrm_design_evaluation`, `mfrm_signal_detection`,
#'   `mfrm_population_prediction`, `mfrm_future_branch_active_branch`,
#'   `mfrm_facets_run`, `mfrm_bias`, `mfrm_anchor_audit`,
#'   `mfrm_linking_review`, `mfrm_misfit_casebook`, `mfrm_weighting_audit`,
#'   `mfrm_unit_prediction`, or `mfrm_plausible_values` object, or one of
#'   their `summary()` outputs.
#' @param which Optional character vector selecting a subset of named tables.
#' @param appendix_preset Optional appendix-oriented table preset:
#'   `"all"`, `"recommended"`, `"compact"`, `"methods"`, `"results"`,
#'   `"diagnostics"`, or `"reporting"`. Cannot be combined with `which`.
#'   Section-aware presets keep returned tables whose bundle catalog maps to
#'   the requested appendix section.
#' @param include_empty If `TRUE`, retain empty tables in the returned bundle.
#' @param digits Digits forwarded when `summary()` must be computed from a raw
#'   object.
#' @param top_n Row cap forwarded to compact `summary()` methods when `x` is a
#'   raw object.
#' @param preview_chars Character cap forwarded to
#'   `summary.mfrm_apa_outputs()` when `x` is a raw APA-output object.
#'
#' @details
#' This helper turns the package's compact summary objects into a reproducible
#' table bundle for manuscript drafting, appendix handoff, or downstream
#' formatting. It does not replace [apa_table()]; instead, it provides a
#' consistent bridge from `summary()` to named `data.frame` components that can
#' later be rendered with [apa_table()] or exported directly.
#'
#' The public entry point validates `x` and the summary-object contract up
#' front, so malformed summaries fail with a package-level message instead of
#' falling through to opaque downstream errors.
#'
#' The function first normalizes `x` through the corresponding `summary()`
#' method when needed, then records a `table_index` describing every available
#' table and returns the selected tables in `tables`. Optional appendix presets
#' can be applied at bundle-construction time when you want a conservative
#' manuscript-facing subset before plotting or export.
#'
#' @section Supported inputs:
#' - [fit_mfrm()] or `summary(fit)`
#' - [diagnose_mfrm()] or `summary(diag)`
#' - [describe_mfrm_data()] or `summary(ds)`
#' - [reporting_checklist()] or `summary(chk)`
#' - [build_apa_outputs()] or `summary(apa)`
#' - [evaluate_mfrm_design()] or `summary(sim_eval)`
#' - [evaluate_mfrm_signal_detection()] or `summary(sig_eval)`
#' - [predict_mfrm_population()] or `summary(pred)`
#' - `planning_schema$future_branch_active_branch` or `summary(...)`
#' - [run_mfrm_facets()] or `summary(out)`
#' - [estimate_bias()] or `summary(bias)`
#' - [audit_mfrm_anchors()] or `summary(audit)`
#' - [build_linking_review()] or `summary(review)`
#' - [build_misfit_casebook()] or `summary(casebook)`
#' - [build_weighting_audit()] or `summary(audit)`
#' - [predict_mfrm_units()] or `summary(pred_units)`
#' - [sample_mfrm_plausible_values()] or `summary(pv)`
#'
#' @section Interpreting output:
#' - `overview`: one-row metadata about the source summary and table counts.
#' - `table_index`: table names, dimensions, roles, and manuscript-oriented
#'   descriptions.
#' - `plot_index`: which returned tables contain numeric content and which
#'   bundle-level plot types can use them directly.
#' - `tables`: named `data.frame` objects ready for formatting or export.
#' - `appendix_preset`: active appendix subset mode (`"none"` when not used).
#' - `notes`: short guidance about omitted empty tables or source-level caveats.
#' - fit-level caveats use the `analysis_caveats` role; pre-fit data
#'   score-support caveats use the `score_category_caveats` role. Both roles are
#'   classified as diagnostics and stay in `recommended` appendix subsets.
#' - latent-regression fit summaries expose `population_coding` in the methods
#'   appendix role so categorical levels, contrasts, and encoded columns can be
#'   documented with the coefficient table.
#'
#' @section Typical workflow:
#' 1. Build a compact object with `summary(...)`.
#' 2. Convert it with `build_summary_table_bundle(...)`.
#' 3. Use `bundle$tables[[...]]` directly, or hand a selected table to
#'    [apa_table()] for formatted manuscript output.
#' 4. If you want a manuscript appendix subset up front, use a preset such as
#'    `appendix_preset = "recommended"`, `"compact"`, or `"diagnostics"`.
#'
#' @return An object of class `mfrm_summary_table_bundle` with:
#' - `overview`
#' - `table_index`
#' - `plot_index`
#' - `tables`
#' - `appendix_preset`
#' - `notes`
#' - `source_class`
#' - `summary_class`
#'
#' @seealso [summary()], [apa_table()], [reporting_checklist()],
#'   [build_apa_outputs()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' bundle <- build_summary_table_bundle(fit)
#' bundle$table_index
#' summary(bundle)$role_summary
#' }
#' @export
build_summary_table_bundle <- function(x,
                                       which = NULL,
                                       appendix_preset = NULL,
                                       include_empty = FALSE,
                                       digits = 3,
                                       top_n = 10,
                                       preview_chars = 160) {
  validated <- validate_summary_table_bundle_inputs(
    x = x,
    which = which,
    appendix_preset = appendix_preset,
    include_empty = include_empty,
    digits = digits,
    top_n = top_n,
    preview_chars = preview_chars,
    helper = "build_summary_table_bundle()"
  )
  which <- validated$which
  appendix_preset <- validated$appendix_preset
  include_empty <- validated$include_empty
  resolved <- validated$resolved

  if (!is.null(appendix_preset) && !is.null(which)) {
    stop(
      "`build_summary_table_bundle()` requires `appendix_preset` and `which` to be used separately.",
      call. = FALSE
    )
  }
  spec <- summary_table_bundle_spec(resolved$summary)
  tables <- spec$tables
  table_index <- build_summary_table_index(tables, spec$roles, spec$descriptions)

  requested <- names(tables)
  if (!is.null(which)) {
    which <- unique(as.character(which))
    unknown <- setdiff(which, names(tables))
    if (length(unknown) > 0L) {
      stop(
        "`build_summary_table_bundle()` received unknown `which` table name(s): ",
        paste(unknown, collapse = ", "),
        ". Inspect `build_summary_table_bundle(x)$table_index$Table` for supported names.",
        call. = FALSE
      )
    }
    requested <- which
  }

  if (is.null(which) && !isTRUE(include_empty)) {
    keep <- vapply(tables[requested], function(tbl) !summary_table_bundle_is_empty(tbl), logical(1))
    requested <- requested[keep]
  }

  selected_tables <- tables[requested]
  selected_index <- table_index[match(requested, table_index$Table), , drop = FALSE]
  plot_index <- summary_table_bundle_plot_index(selected_tables)
  dropped_empty <- sum(vapply(tables, summary_table_bundle_is_empty, logical(1))) -
    sum(vapply(selected_tables, summary_table_bundle_is_empty, logical(1)))
  notes <- as.character(resolved$summary$notes %||% character(0))
  if (is.null(which) && !isTRUE(include_empty) && dropped_empty > 0L) {
    notes <- c(notes, sprintf("%d empty table(s) were omitted from `tables`; use `include_empty = TRUE` to retain them.", dropped_empty))
  }
  if (!is.null(which)) {
    notes <- c(notes, sprintf("Returned %d requested table(s): %s.", length(requested), paste(requested, collapse = ", ")))
  }

  overview <- data.frame(
    Title = spec$title,
    SourceClass = resolved$source_class,
    SummaryClass = resolved$summary_class,
    TablesAvailable = nrow(table_index),
    TablesReturned = length(selected_tables),
    AppendixPreset = if (is.null(appendix_preset)) "none" else appendix_preset,
    stringsAsFactors = FALSE
  )

  out <- list(
    overview = overview,
    table_index = selected_index,
    plot_index = plot_index,
    tables = selected_tables,
    notes = unique(notes[nzchar(notes)]),
    appendix_preset = appendix_preset %||% "none",
    source_class = resolved$source_class,
    summary_class = resolved$summary_class
  )
  class(out) <- "mfrm_summary_table_bundle"
  if (!is.null(appendix_preset)) {
    out <- summary_table_bundle_select_for_appendix(out, preset = appendix_preset)
    out$appendix_preset <- appendix_preset
    if (!is.null(out$overview) && nrow(out$overview) > 0L) {
      out$overview$AppendixPreset <- appendix_preset
    }
  }
  out
}

#' @export
print.mfrm_summary_table_bundle <- function(x, ...) {
  cat("mfrmr Summary Table Bundle\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(as.data.frame(x$overview), row.names = FALSE)
  }
  if (!is.null(x$table_index) && nrow(x$table_index) > 0) {
    cat("\nTable index\n")
    print(as.data.frame(x$table_index), row.names = FALSE)
  }
  if (!is.null(x$plot_index) && nrow(x$plot_index) > 0) {
    cat("\nPlot index\n")
    print(as.data.frame(x$plot_index), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

summary_table_bundle_first_numeric_table <- function(bundle) {
  plot_idx <- as.data.frame(bundle$plot_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(plot_idx) > 0L && all(c("Table", "PlotReady") %in% names(plot_idx))) {
    ready <- plot_idx[plot_idx$PlotReady %in% TRUE, , drop = FALSE]
    if (nrow(ready) > 0L) {
      return(as.character(ready$Table[1]))
    }
  }
  tbls <- bundle$tables %||% list()
  if (length(tbls) == 0L) return(NULL)
  for (nm in names(tbls)) {
    tbl <- as.data.frame(tbls[[nm]], stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L) next
    if (any(vapply(tbl, is.numeric, logical(1)))) {
      return(nm)
    }
  }
  NULL
}

summary_table_bundle_profile <- function(bundle) {
  tbls <- bundle$tables %||% list()
  if (length(tbls) == 0L) {
    return(data.frame())
  }
  idx <- as.data.frame(bundle$table_index %||% data.frame(), stringsAsFactors = FALSE)
  do.call(
    rbind,
    lapply(names(tbls), function(nm) {
      tbl <- as.data.frame(tbls[[nm]], stringsAsFactors = FALSE)
      idx_row <- if (nrow(idx) > 0L && "Table" %in% names(idx)) {
        idx[idx$Table %in% nm, , drop = FALSE]
      } else {
        data.frame()
      }
      data.frame(
        Table = nm,
        Rows = nrow(tbl),
        Cols = ncol(tbl),
        NumericColumns = sum(vapply(tbl, is.numeric, logical(1))),
        MissingValues = sum(is.na(tbl)),
        Role = if (nrow(idx_row) > 0L && "Role" %in% names(idx_row)) as.character(idx_row$Role[1]) else "",
        Description = if (nrow(idx_row) > 0L && "Description" %in% names(idx_row)) as.character(idx_row$Description[1]) else "",
        stringsAsFactors = FALSE
      )
    })
  )
}

summary_table_bundle_plot_index <- function(tables) {
  tbls <- tables %||% list()
  if (length(tbls) == 0L) {
    return(data.frame())
  }
  do.call(
    rbind,
    lapply(names(tbls), function(nm) {
      tbl <- as.data.frame(tbls[[nm]], stringsAsFactors = FALSE)
      numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
      plot_ready <- nrow(tbl) > 0L && length(numeric_cols) > 0L
      data.frame(
        Table = nm,
        PlotReady = plot_ready,
        NumericColumns = length(numeric_cols),
        DefaultPlotTypes = if (plot_ready) "numeric_profile, first_numeric" else "",
        stringsAsFactors = FALSE
      )
    })
  )
}

summary_table_bundle_compact_labels <- function(x, max_n = 4L) {
  vals <- unique(as.character(x %||% character(0)))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0L) return("")
  max_n <- max(1L, as.integer(max_n))
  if (length(vals) <= max_n) {
    return(paste(vals, collapse = ", "))
  }
  paste(c(vals[seq_len(max_n)], "..."), collapse = ", ")
}

summary_table_bundle_appendix_role_registry <- function() {
  out <- data.frame(
    Role = c(
      "run_overview",
      "population_basis",
      "population_coefficients",
      "population_design",
      "population_coding",
      "facet_distribution",
      "person_distribution",
      "category_structure",
      "gpcm_discrimination",
      "estimation_settings",
      "overall_fit",
      "precision_basis",
      "precision_audit",
      "facet_precision",
      "flag_counts",
      "missingness",
      "score_usage",
      "facet_coverage",
      "agreement",
      "checklist_overview",
      "section_coverage",
      "priority_distribution",
      "draft_overview",
      "component_stats",
      "draft_checks",
      "text_preview",
      "reporting_map",
      "extreme_facet_levels",
      "extreme_person_high",
      "extreme_person_low",
      "extreme_fit_rows",
      "draft_actions",
      "checklist_settings",
      "future_branch_overview",
      "future_branch_profile",
      "future_branch_load_balance",
      "future_branch_coverage",
      "future_branch_guardrails",
      "future_branch_readiness",
      "future_branch_recommendation",
      "future_branch_appendix_presets",
      "future_branch_appendix_roles",
      "future_branch_appendix_sections",
      "future_branch_selection_table_presets",
      "future_branch_selection_handoff_tables",
      "future_branch_selection_handoff_presets",
      "future_branch_selection_handoff",
      "future_branch_selection_handoff_bundles",
      "future_branch_selection_handoff_roles",
      "future_branch_selection_handoff_role_sections",
      "future_branch_selection_tables",
      "future_branch_selection_summary",
      "future_branch_selection_roles",
      "future_branch_selection_sections",
      "future_branch_selection_catalog",
      "future_branch_reporting_map",
      "linking_review_overview",
      "linking_risk_screen",
      "linking_risk_group_index",
      "misfit_casebook_overview",
      "misfit_case_rows",
      "misfit_case_rollup",
      "misfit_case_sources",
      "weighting_review_overview",
      "review_status",
      "plot_routing",
      "capability_boundary",
      "repair_recommendations",
      "interpretation_notes",
      "reweighting_measure_shift",
      "score_category_caveats",
      "analysis_caveats",
      "prediction_overview",
      "unit_estimates",
      "prediction_audit",
      "population_audit",
      "scoring_settings",
      "plausible_value_overview",
      "plausible_value_draws",
      "design_performance",
      "signal_detection",
      "design_grid",
      "forecast_summary",
      "workflow_overview",
      "column_mapping",
      "workflow_settings",
      "bias_overview",
      "bias_chi_square",
      "bias_iteration_status",
      "bias_screening_rows",
      "anchor_audit_overview",
      "anchor_issue_counts",
      "level_observation_audit",
      "category_usage",
      "prefit_anchor_risks",
      "drift_risks",
      "chain_risks",
      "review_settings"
    ),
    AppendixSection = c(
      "methods",
      "methods",
      "methods",
      "methods",
      "methods",
      "results",
      "results",
      "results",
      "results",
      "methods",
      "results",
      "results",
      "diagnostics",
      "results",
      "diagnostics",
      "methods",
      "methods",
      "results",
      "reporting",
      "reporting",
      "reporting",
      "reporting",
      "reporting",
      "reporting",
      "workflow",
      "workflow",
      "workflow",
      "exploratory",
      "exploratory",
      "exploratory",
      "exploratory",
      "workflow",
      "workflow",
      "methods",
      "methods",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "methods",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "workflow",
      "diagnostics",
      "workflow",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "methods",
      "methods",
      "diagnostics",
      "workflow",
      "methods",
      "reporting",
      "reporting",
      "results",
      "diagnostics",
      "diagnostics",
      "results",
      "results",
      "diagnostics",
      "diagnostics",
      "methods",
      "results",
      "results",
      "results",
      "diagnostics",
      "methods",
      "results",
      "workflow",
      "workflow",
      "workflow",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "methods"
    ),
    RecommendedAppendix = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
      TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    ),
    CompactAppendix = c(
      TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE,
      TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE,
      TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE
    ),
    PreferredAppendixOrder = c(
      10L, 20L, 30L, 35L, 36L, 40L, 50L, 60L, 70L, 80L,
      90L, 100L, 110L, 120L, 130L, 140L, 150L, 160L, 170L,
      180L, 190L, 200L, 210L, 220L,
      990L, 991L, 900L,
      910L, 920L, 930L, 940L, 950L, 960L,
      225L, 226L, 227L, 228L, 229L, 230L, 231L,
      970L, 971L, 972L, 973L, 974L, 975L, 976L, 977L, 978L, 979L, 980L, 981L, 982L, 983L, 984L, 985L,
      232L, 233L, 234L, 235L, 236L, 237L, 238L, 239L, 240L, 986L, 987L, 988L, 989L, 125L, 131L, 132L,
      133L, 134L, 135L, 136L, 137L, 138L, 139L,
      241L, 242L, 243L, 244L, 245L, 246L, 247L, 248L, 249L, 250L, 251L,
      252L, 253L, 254L, 255L, 256L, 257L, 258L, 259L
    ),
    AppendixRationale = c(
      "Always include the main run-identification table.",
      "Include whenever population-model interpretation is part of the report.",
      "Include when latent-regression coefficients were estimated.",
      "Include to audit latent-regression design-matrix columns and variance screening.",
      "Include when categorical population covariates were encoded through the model matrix.",
      "Core facet spread and scale-location appendix table.",
      "Useful for full appendices but omitted from compact presets.",
      "Core threshold/category appendix table.",
      "Useful when GPCM discrimination is active.",
      "Methods/settings appendix table; recommended but not compact.",
      "Core fit summary table for results appendices.",
      "Core precision-basis table for cautious interpretation.",
      "Recommended when precision caveats need explicit documentation.",
      "Core reliability/separation appendix table.",
      "Recommended diagnostic count surface for QC appendices.",
      "Core missing-data appendix table.",
      "Core score-usage appendix table.",
      "Core facet-coverage appendix table.",
      "Optional agreement appendix surface.",
      "Recommended checklist overview for reporting QA appendices.",
      "Core section-coverage appendix table.",
      "Core priority distribution for reporting follow-up.",
      "Core manuscript-draft coverage overview.",
      "Core APA component inventory.",
      "Internal draft QA surface; keep out of recommended presets.",
      "Preview-only draft text; keep out of recommended presets.",
      "Bridge metadata, useful for workflow but not manuscript appendix.",
      "Exploratory extreme table; available only in full exports.",
      "Exploratory extreme table; available only in full exports.",
      "Exploratory extreme table; available only in full exports.",
      "Exploratory extreme table; available only in full exports.",
      "Internal drafting action list; keep out of recommended presets.",
      "Internal checklist settings; keep out of recommended presets.",
      "Recommended methods appendix overview for the future arbitrary-facet planning scaffold.",
      "Recommended exact-count profile for future arbitrary-facet planning methods appendices.",
      "Detailed load-balance diagnostics; retain for full exports but omit from recommended presets.",
      "Detailed coverage/connectivity diagnostics; retain for full exports but omit from recommended presets.",
      "Detailed guardrail classifications; retain for full exports but omit from recommended presets.",
      "Core structural readiness table for future arbitrary-facet planning review.",
      "Core conservative future-branch recommendation table for methods appendices.",
      "Workflow-only appendix preset counts for direct future-branch review.",
      "Workflow-only appendix role counts for direct future-branch review.",
      "Workflow-only appendix section counts for direct future-branch review.",
      "Workflow-only preset-specific appendix table selections for direct future-branch review.",
      "Workflow-only table-level appendix handoff crosswalk for direct future-branch review.",
      "Workflow-only preset-level appendix handoff overview for direct future-branch review.",
      "Workflow-only manuscript-section handoff summary for direct future-branch review.",
      "Workflow-only bundle-aware appendix handoff summary for direct future-branch review.",
      "Workflow-only role-aware appendix handoff summary for direct future-branch review.",
      "Workflow-only role-by-section appendix handoff summary for direct future-branch review.",
      "Workflow-only preset-aware appendix table selections for direct future-branch review.",
      "Workflow-only preset-filtered appendix bundle counts for direct future-branch review.",
      "Workflow-only preset-filtered appendix role counts for direct future-branch review.",
      "Workflow-only preset-filtered appendix section counts for direct future-branch review.",
      "Workflow-only preset-filtered appendix selection catalog for direct future-branch review.",
      "Workflow-only reporting bridge metadata for the direct future-branch surface.",
      "Recommended overview table for linking-review appendix handoff.",
      "Recommended top-risk table for operational linking-review follow-up appendices.",
      "Recommended grouping-view index for operational linking-review triage.",
      "Recommended overview table for operational misfit-case review appendix handoff.",
      "Recommended top-case table for operational misfit follow-up appendices.",
      "Recommended rollup table showing where flagged cases concentrate across review groupings.",
      "Recommended source-family count table for operational misfit follow-up appendices.",
      "Recommended overview table for weighting-policy review appendix handoff.",
      "Recommended compact status table for review-oriented appendix handoff.",
      "Plot-routing metadata; keep out of recommended appendix presets.",
      "Capability-boundary statement for supported-with-caveat review helpers.",
      "Recommended action-oriented table for repair or follow-up planning.",
      "Interpretation notes; retain mainly in full reporting exports.",
      "Recommended reweighting-change table for bounded GPCM comparison review.",
      "Recommended caveat table for retained zero-count score categories and related score-support warnings.",
      "Recommended fit-level caveat table for score-support, population-model, and other analysis warnings.",
      "Recommended overview table for posterior unit-scoring appendix handoff.",
      "Recommended posterior estimate table for scored-person appendix handoff.",
      "Row-level scoring audit; retain for full exports but omit from compact presets.",
      "Recommended latent-regression scoring omission audit when population-model scoring is active.",
      "Methods/settings appendix table for posterior scoring inputs; recommended but not compact.",
      "Recommended overview table for plausible-value appendix handoff.",
      "Recommended plausible-value draw summary table for posterior-scoring appendices.",
      "Recommended design-performance table for simulation design appendices.",
      "Recommended signal-detection table for simulation diagnostics appendices.",
      "Recommended design-grid table documenting requested forecast inputs.",
      "Recommended forecast-summary table for population prediction appendices.",
      "Workflow-only overview for one-shot run handoff; omit from manuscript presets.",
      "Workflow-only column mapping for replay and audit handoff.",
      "Workflow-only run settings for one-shot workflow provenance.",
      "Recommended overview table for bias-screening appendix handoff.",
      "Recommended fixed-effect chi-square table for bias-screening appendices.",
      "Bias iteration status table; retain for full exports but omit from compact presets.",
      "Recommended ranked bias-screening row table for immediate follow-up.",
      "Recommended overview table for anchor-audit appendix handoff.",
      "Recommended anchor issue-count table for pre-fit audit appendices.",
      "Recommended level-observation audit table; retain for full appendices but omit from compact presets.",
      "Recommended score-category usage table for anchor-audit support review.",
      "Recommended pre-fit anchor-risk table for linking-review appendices.",
      "Recommended drift-risk table for linking-review appendices.",
      "Recommended chain-risk table for linking-review appendices.",
      "Methods/settings table for operational review helpers; recommended but not compact."
    ),
    stringsAsFactors = FALSE
  )
  capability_boundary <- out$Role %in% "capability_boundary"
  out$CompactAppendix[capability_boundary] <- TRUE
  out$PreferredAppendixOrder[capability_boundary] <- 240.5
  out
}

summary_table_bundle_catalog <- function(bundle) {
  idx <- as.data.frame(bundle$table_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(idx) == 0L) {
    return(data.frame())
  }

  plot_idx <- as.data.frame(bundle$plot_index %||% data.frame(), stringsAsFactors = FALSE)
  plot_keep <- intersect(c("Table", "PlotReady", "NumericColumns", "DefaultPlotTypes"), names(plot_idx))
  if (length(plot_keep) > 0L) {
    plot_idx <- plot_idx[, plot_keep, drop = FALSE]
    idx <- merge(idx, plot_idx, by = "Table", all.x = TRUE, sort = FALSE)
  }

  if (!"PlotReady" %in% names(idx)) idx$PlotReady <- FALSE
  if (!"NumericColumns" %in% names(idx)) idx$NumericColumns <- 0L
  if (!"DefaultPlotTypes" %in% names(idx)) idx$DefaultPlotTypes <- ""

  idx$PlotReady[is.na(idx$PlotReady)] <- FALSE
  idx$NumericColumns[is.na(idx$NumericColumns)] <- 0L
  idx$DefaultPlotTypes[is.na(idx$DefaultPlotTypes)] <- ""
  idx$ExportReady <- TRUE
  idx$ApaTableReady <- TRUE
  idx$RecommendedBridge <- ifelse(
    idx$PlotReady %in% TRUE,
    "apa_table() / plot(bundle)",
    "apa_table() / export_summary_appendix()"
  )
  appendix_registry <- summary_table_bundle_appendix_role_registry()
  appendix_idx <- match(as.character(idx$Role), appendix_registry$Role)
  idx$RecommendedAppendix <- appendix_registry$RecommendedAppendix[appendix_idx]
  idx$CompactAppendix <- appendix_registry$CompactAppendix[appendix_idx]
  idx$PreferredAppendixOrder <- appendix_registry$PreferredAppendixOrder[appendix_idx]
  idx$AppendixRationale <- appendix_registry$AppendixRationale[appendix_idx]
  idx$AppendixSection <- appendix_registry$AppendixSection[appendix_idx]
  idx$RecommendedAppendix[is.na(idx$RecommendedAppendix)] <- FALSE
  idx$CompactAppendix[is.na(idx$CompactAppendix)] <- FALSE
  idx$PreferredAppendixOrder[is.na(idx$PreferredAppendixOrder)] <- 999L
  idx$AppendixRationale[is.na(idx$AppendixRationale)] <- "Available only through full appendix export."
  idx$AppendixSection[is.na(idx$AppendixSection)] <- "workflow"

  idx[, c(
    "Table", "Rows", "Cols", "Role", "Description",
    "PlotReady", "NumericColumns", "DefaultPlotTypes",
    "ExportReady", "ApaTableReady", "RecommendedBridge",
    "AppendixSection",
    "RecommendedAppendix", "CompactAppendix",
    "PreferredAppendixOrder", "AppendixRationale"
  ), drop = FALSE]
}

summary_table_bundle_appendix_presets <- function(catalog) {
  catalog <- as.data.frame(catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(catalog) == 0L) {
    return(data.frame())
  }

  preset_defs <- list(
    all = rep(TRUE, nrow(catalog)),
    recommended = catalog$RecommendedAppendix %in% TRUE,
    compact = catalog$CompactAppendix %in% TRUE,
    methods = catalog$AppendixSection %in% "methods",
    results = catalog$AppendixSection %in% "results",
    diagnostics = catalog$AppendixSection %in% "diagnostics",
    reporting = catalog$AppendixSection %in% "reporting"
  )
  preset_uses <- c(
    all = "Complete appendix handoff with every returned summary table.",
    recommended = "Manuscript appendix without bridge-only or preview-only surfaces.",
    compact = "Reviewer-facing compact appendix focused on core design and fit summaries.",
    methods = "Methods appendix subset focused on design, scoring basis, and settings.",
    results = "Results appendix subset focused on fit, precision, and scale summaries.",
    diagnostics = "Diagnostics appendix subset focused on caveats, flags, and precision checks.",
    reporting = "Reporting appendix subset focused on manuscript/checklist coverage surfaces."
  )

  out <- do.call(
    rbind,
    lapply(names(preset_defs), function(preset_nm) {
      part <- catalog[preset_defs[[preset_nm]], , drop = FALSE]
      data.frame(
        Preset = preset_nm,
        Tables = nrow(part),
        PlotReadyTables = sum(part$PlotReady %in% TRUE, na.rm = TRUE),
        RolesCovered = length(unique(as.character(part$Role))),
        SectionsCovered = summary_table_bundle_compact_labels(unique(as.character(part$AppendixSection)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(part$Table, max_n = 4L),
        PrimaryUse = unname(preset_uses[[preset_nm]]),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

summary_table_bundle_selection_surface <- function(bundle, surface) {
  tables <- bundle$tables %||% list()
  if (!is.list(tables) || length(tables) == 0L) {
    return(data.frame())
  }

  candidates <- switch(
    as.character(surface[1] %||% ""),
    selection_summary = c("future_branch_selection_summary", "appendix_selection_summary", "selection_summary"),
    selection_table_summary = c("future_branch_selection_tables", "appendix_selection_table_summary", "selection_table_summary"),
    selection_table_preset_summary = c("future_branch_selection_table_presets", "selection_table_preset_summary"),
    selection_handoff_table_summary = c("future_branch_selection_handoff_tables", "appendix_selection_handoff_table_summary", "selection_handoff_table_summary"),
    selection_handoff_preset_summary = c("future_branch_selection_handoff_presets", "appendix_selection_handoff_preset_summary", "selection_handoff_preset_summary"),
    selection_handoff_summary = c("future_branch_selection_handoff", "appendix_selection_handoff_summary", "selection_handoff_summary"),
    selection_handoff_bundle_summary = c("future_branch_selection_handoff_bundles", "appendix_selection_handoff_bundle_summary", "selection_handoff_bundle_summary"),
    selection_handoff_role_summary = c("future_branch_selection_handoff_roles", "appendix_selection_handoff_role_summary", "selection_handoff_role_summary"),
    selection_handoff_role_section_summary = c("future_branch_selection_handoff_role_sections", "appendix_selection_handoff_role_section_summary", "selection_handoff_role_section_summary"),
    selection_role_summary = c("future_branch_selection_roles", "appendix_selection_role_summary", "selection_role_summary"),
    selection_section_summary = c("future_branch_selection_sections", "appendix_selection_section_summary", "selection_section_summary"),
    selection_catalog = c("future_branch_selection_catalog", "appendix_selection_catalog", "selection_catalog"),
    character(0)
  )

  hit <- candidates[candidates %in% names(tables)]
  if (length(hit) == 0L) {
    return(data.frame())
  }

  as.data.frame(tables[[hit[1]]], stringsAsFactors = FALSE)
}

summary_table_bundle_appendix_role_summary <- function(catalog) {
  catalog <- as.data.frame(catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(catalog) == 0L || !"Role" %in% names(catalog)) {
    return(data.frame())
  }

  split_tbl <- split(catalog, as.character(catalog$Role %||% ""))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(role_nm) {
      part <- split_tbl[[role_nm]]
      data.frame(
        Role = as.character(role_nm),
        Tables = nrow(part),
        PlotReadyTables = sum(part$PlotReady %in% TRUE, na.rm = TRUE),
        RecommendedTables = sum(part$RecommendedAppendix %in% TRUE, na.rm = TRUE),
        CompactTables = sum(part$CompactAppendix %in% TRUE, na.rm = TRUE),
        SectionsCovered = summary_table_bundle_compact_labels(unique(as.character(part$AppendixSection)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Tables, out$Role, decreasing = TRUE), , drop = FALSE]
}

summary_table_bundle_appendix_section_summary <- function(catalog) {
  catalog <- as.data.frame(catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(catalog) == 0L || !"AppendixSection" %in% names(catalog)) {
    return(data.frame())
  }
  sections <- split(catalog, as.character(catalog$AppendixSection))
  out <- do.call(
    rbind,
    lapply(names(sections), function(section_nm) {
      part <- sections[[section_nm]]
      data.frame(
        AppendixSection = section_nm,
        Tables = nrow(part),
        PlotReadyTables = sum(part$PlotReady %in% TRUE, na.rm = TRUE),
        RecommendedTables = sum(part$RecommendedAppendix %in% TRUE, na.rm = TRUE),
        CompactTables = sum(part$CompactAppendix %in% TRUE, na.rm = TRUE),
        RolesCovered = length(unique(as.character(part$Role))),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Tables, out$AppendixSection, decreasing = TRUE), , drop = FALSE]
}

summary_table_bundle_subset <- function(bundle, which, note = NULL) {
  if (!inherits(bundle, "mfrm_summary_table_bundle")) {
    stop("`bundle` must be an mfrm_summary_table_bundle object.", call. = FALSE)
  }
  tables <- bundle$tables %||% list()
  keep <- intersect(as.character(which %||% character(0)), names(tables))
  if (length(keep) == 0L) {
    stop("No matching tables were found in the supplied summary-table bundle.", call. = FALSE)
  }

  out <- bundle
  out$tables <- tables[keep]

  idx <- as.data.frame(bundle$table_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(idx) > 0L && "Table" %in% names(idx)) {
    keep_idx <- match(keep, idx$Table)
    keep_idx <- keep_idx[!is.na(keep_idx)]
    out$table_index <- idx[keep_idx, , drop = FALSE]
  }

  plot_idx <- as.data.frame(bundle$plot_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(plot_idx) > 0L && "Table" %in% names(plot_idx)) {
    keep_plot <- match(keep, plot_idx$Table)
    keep_plot <- keep_plot[!is.na(keep_plot)]
    out$plot_index <- plot_idx[keep_plot, , drop = FALSE]
  }

  out$overview <- as.data.frame(bundle$overview %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(out$overview) > 0L && "TablesReturned" %in% names(out$overview)) {
    out$overview$TablesReturned <- length(keep)
  }

  if (!is.null(note) && nzchar(as.character(note[1] %||% ""))) {
    out$notes <- unique(c(bundle$notes %||% character(0), as.character(note[1])))
  }
  out$appendix_preset <- as.character(bundle$appendix_preset %||% "none")

  out
}

summary_table_bundle_empty_subset <- function(bundle, note = NULL) {
  if (!inherits(bundle, "mfrm_summary_table_bundle")) {
    stop("`bundle` must be an mfrm_summary_table_bundle object.", call. = FALSE)
  }
  out <- bundle
  out$tables <- list()
  out$table_index <- as.data.frame(bundle$table_index %||% data.frame(), stringsAsFactors = FALSE)[0, , drop = FALSE]
  out$plot_index <- as.data.frame(bundle$plot_index %||% data.frame(), stringsAsFactors = FALSE)[0, , drop = FALSE]
  out$overview <- as.data.frame(bundle$overview %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(out$overview) > 0L && "TablesReturned" %in% names(out$overview)) {
    out$overview$TablesReturned <- 0L
  }
  if (!is.null(note) && nzchar(as.character(note[1] %||% ""))) {
    out$notes <- unique(c(bundle$notes %||% character(0), as.character(note[1])))
  }
  out$appendix_preset <- as.character(bundle$appendix_preset %||% "none")
  out
}

summary_table_bundle_select_for_appendix <- function(bundle,
                                                     preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  if (!inherits(bundle, "mfrm_summary_table_bundle")) {
    stop("`bundle` must be an mfrm_summary_table_bundle object.", call. = FALSE)
  }
  if (identical(preset, "all")) {
    return(bundle)
  }

  catalog <- summary_table_bundle_catalog(bundle)
  keep <- switch(
    preset,
    recommended = as.character(catalog$Table[catalog$RecommendedAppendix %in% TRUE]),
    compact = as.character(catalog$Table[catalog$CompactAppendix %in% TRUE]),
    methods = as.character(catalog$Table[catalog$AppendixSection %in% "methods"]),
    results = as.character(catalog$Table[catalog$AppendixSection %in% "results"]),
    diagnostics = as.character(catalog$Table[catalog$AppendixSection %in% "diagnostics"]),
    reporting = as.character(catalog$Table[catalog$AppendixSection %in% "reporting"])
  )
  keep <- keep[nzchar(keep)]
  if (length(keep) == 0L) {
    if (preset %in% c("methods", "results", "diagnostics", "reporting")) {
      return(summary_table_bundle_empty_subset(
        bundle,
        note = sprintf("Appendix preset `%s` matched no tables in this bundle.", preset)
      ))
    }
    keep <- if ("overview" %in% names(bundle$tables)) "overview" else names(bundle$tables)[1]
  }

  note <- sprintf(
    "Appendix preset `%s` selected %d table(s): %s.",
    preset,
    length(keep),
    paste(keep, collapse = ", ")
  )
  summary_table_bundle_subset(bundle, keep, note = note)
}

summary_table_bundle_reporting_map <- function(bundle, catalog) {
  numeric_ready <- if (nrow(catalog) > 0L) sum(catalog$PlotReady %in% TRUE, na.rm = TRUE) else 0L
  data.frame(
    Area = c(
      "Coverage overview",
      "Table catalog / manuscript selection",
      "Numeric QC and quick plotting",
      "APA / appendix bridge",
      "Source-level caveats"
    ),
    CoveredHere = c("yes", "yes", "yes", "yes", "partial"),
    CompanionOutput = c(
      "summary(bundle)$overview / role_summary",
      "summary(bundle)$table_catalog / bundle$table_index",
      "summary(bundle)$plot_index / plot(bundle, ...)",
      "apa_table(bundle, which = ...) / export_summary_appendix(bundle, preset = \"recommended\")",
      "bundle$notes and the originating summary()/diagnostics output"
    ),
    stringsAsFactors = FALSE
  )
}

#' Summarize a summary-table bundle for manuscript QC
#'
#' @param object Output from [build_summary_table_bundle()].
#' @param digits Number of digits used for numeric summaries.
#' @param top_n Maximum number of table-profile rows to keep.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is designed to answer a manuscript-facing question: which
#' reporting tables are available, how large are they, which roles do they
#' serve, and which of them contain numeric content suitable for quick plotting
#' or appendix export.
#'
#' @section Interpreting output:
#' - `overview`: source class, returned-table count, note count, and whether a
#'   numeric table is available for plotting.
#' - `role_summary`: counts and total size by reporting role.
#' - `table_catalog`: complete returned-table registry with plot/export bridges.
#' - `table_profile`: table-level dimensions, numeric-column counts, and missing
#'   values for the largest returned tables.
#' - `plot_index`: which returned tables are plot-ready and which bundle-level
#'   numeric QC routes they support.
#' - `appendix_presets`: conservative `all` / `recommended` / `compact`
#'   plus section-aware `methods` / `results` / `diagnostics` / `reporting`
#'   appendix-export presets derived from table roles.
#' - `appendix_role_summary`: counts of returned tables by reporting role under
#'   the same conservative appendix routing used by the bundle catalog.
#' - `appendix_section_summary`: counts of returned tables by manuscript-facing
#'   appendix section.
#' - `selection_handoff_table_summary`: workflow-only table-level appendix
#'   handoff crosswalk when present in the bundle.
#' - `selection_handoff_preset_summary`: workflow-only appendix handoff overview
#'   aggregated at the preset level when present in the bundle.
#' - `selection_handoff_bundle_summary`: workflow-only appendix handoff
#'   overview aggregated at the bundle-by-section level when present in the
#'   bundle.
#' - `selection_handoff_role_summary`: workflow-only appendix handoff overview
#'   aggregated at the reporting-role level when present in the bundle.
#' - `selection_handoff_role_section_summary`: workflow-only appendix handoff
#'   overview aggregated at the reporting-role by appendix-section level when
#'   present in the bundle.
#' - `selection_summary`, `selection_table_summary`,
#'   `selection_table_preset_summary`, `selection_role_summary`,
#'   `selection_section_summary`, and `selection_catalog`: preset-filtered
#'   appendix selection surfaces when workflow-only handoff tables are embedded
#'   in the bundle.
#' - `reporting_map`: where to go next for plotting, APA formatting, and export.
#' - `notes`: carried forward source-level caveats from the originating summary.
#'
#' @section Typical workflow:
#' 1. Build `bundle <- build_summary_table_bundle(summary(...))`.
#' 2. Run `summary(bundle)` to see reporting coverage.
#' 3. Use `plot(bundle, type = "table_rows")` or
#'    `plot(bundle, type = "numeric_profile", which = ...)` for quick QC.
#'
#' @return An object of class `summary.mfrm_summary_table_bundle`.
#' @seealso [build_summary_table_bundle()], [apa_table()], [plot()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' bundle <- build_summary_table_bundle(fit)
#' summary(bundle)
#' }
#' @export
summary.mfrm_summary_table_bundle <- function(object, digits = 3, top_n = 8, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  tbl_index <- as.data.frame(object$table_index %||% data.frame(), stringsAsFactors = FALSE)
  profile <- summary_table_bundle_profile(object)
  numeric_available <- !is.null(summary_table_bundle_first_numeric_table(object))

  overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
  overview$Notes <- length(object$notes %||% character(0))
  overview$NumericTables <- sum(profile$NumericColumns > 0, na.rm = TRUE)
  overview$AnyNumericTable <- numeric_available

  if (nrow(profile) > 0L) {
    ord <- order(profile$Rows, profile$Cols, decreasing = TRUE, na.last = TRUE)
    profile <- profile[ord, , drop = FALSE]
    profile <- utils::head(profile, n = top_n)
  }

  role_summary <- data.frame()
  plot_index <- as.data.frame(object$plot_index %||% data.frame(), stringsAsFactors = FALSE)
  table_catalog <- summary_table_bundle_catalog(object)
  selection_summary <- summary_table_bundle_selection_surface(object, "selection_summary")
  selection_table_summary <- summary_table_bundle_selection_surface(object, "selection_table_summary")
  selection_table_preset_summary <- summary_table_bundle_selection_surface(object, "selection_table_preset_summary")
  selection_handoff_table_summary <- summary_table_bundle_selection_surface(object, "selection_handoff_table_summary")
  selection_handoff_preset_summary <- summary_table_bundle_selection_surface(object, "selection_handoff_preset_summary")
  selection_handoff_summary <- summary_table_bundle_selection_surface(object, "selection_handoff_summary")
  selection_handoff_bundle_summary <- summary_table_bundle_selection_surface(object, "selection_handoff_bundle_summary")
  selection_handoff_role_summary <- summary_table_bundle_selection_surface(object, "selection_handoff_role_summary")
  selection_handoff_role_section_summary <- summary_table_bundle_selection_surface(object, "selection_handoff_role_section_summary")
  selection_role_summary <- summary_table_bundle_selection_surface(object, "selection_role_summary")
  selection_section_summary <- summary_table_bundle_selection_surface(object, "selection_section_summary")
  selection_catalog <- summary_table_bundle_selection_surface(object, "selection_catalog")
  appendix_presets <- summary_table_bundle_appendix_presets(table_catalog)
  appendix_role_summary <- summary_table_bundle_appendix_role_summary(table_catalog)
  appendix_section_summary <- summary_table_bundle_appendix_section_summary(table_catalog)
  overview$RecommendedAppendixTables <- sum(table_catalog$RecommendedAppendix %in% TRUE, na.rm = TRUE)
  overview$CompactAppendixTables <- sum(table_catalog$CompactAppendix %in% TRUE, na.rm = TRUE)
  if (nrow(tbl_index) > 0L && "Role" %in% names(tbl_index)) {
    roles <- split(tbl_index, tbl_index$Role %||% "")
    role_summary <- do.call(
      rbind,
      lapply(names(roles), function(role_nm) {
        part <- roles[[role_nm]]
        data.frame(
          Role = as.character(role_nm),
          Tables = nrow(part),
          TotalRows = sum(suppressWarnings(as.numeric(part$Rows)), na.rm = TRUE),
          TotalCols = sum(suppressWarnings(as.numeric(part$Cols)), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
    role_summary <- role_summary[order(role_summary$Tables, role_summary$Role, decreasing = TRUE), , drop = FALSE]
  }
  reporting_map <- summary_table_bundle_reporting_map(object, table_catalog)

  out <- list(
    overview = overview,
    role_summary = role_summary,
    table_catalog = table_catalog,
    table_profile = profile,
    plot_index = plot_index,
    appendix_presets = appendix_presets,
    appendix_role_summary = appendix_role_summary,
    appendix_section_summary = appendix_section_summary,
    selection_summary = selection_summary,
    selection_table_summary = selection_table_summary,
    selection_table_preset_summary = selection_table_preset_summary,
    selection_handoff_table_summary = selection_handoff_table_summary,
    selection_handoff_preset_summary = selection_handoff_preset_summary,
    selection_handoff_summary = selection_handoff_summary,
    selection_handoff_bundle_summary = selection_handoff_bundle_summary,
    selection_handoff_role_summary = selection_handoff_role_summary,
    selection_handoff_role_section_summary = selection_handoff_role_section_summary,
    selection_role_summary = selection_role_summary,
    selection_section_summary = selection_section_summary,
    selection_catalog = selection_catalog,
    reporting_map = reporting_map,
    notes = as.character(object$notes %||% character(0)),
    digits = digits
  )
  class(out) <- "summary.mfrm_summary_table_bundle"
  out
}

#' @export
print.summary.mfrm_summary_table_bundle <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("Summary Table Bundle Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$role_summary) && nrow(x$role_summary) > 0L) {
    cat("\nRole summary\n")
    print(round_numeric_df(as.data.frame(x$role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$table_catalog) && nrow(x$table_catalog) > 0L) {
    cat("\nTable catalog\n")
    print(round_numeric_df(as.data.frame(x$table_catalog), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$table_profile) && nrow(x$table_profile) > 0L) {
    cat("\nTable profile\n")
    print(round_numeric_df(as.data.frame(x$table_profile), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$plot_index) && nrow(x$plot_index) > 0L) {
    cat("\nPlot index\n")
    print(round_numeric_df(as.data.frame(x$plot_index), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_presets) && nrow(x$appendix_presets) > 0L) {
    cat("\nAppendix presets\n")
    print(round_numeric_df(as.data.frame(x$appendix_presets), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_role_summary) && nrow(x$appendix_role_summary) > 0L) {
    cat("\nAppendix role summary\n")
    print(round_numeric_df(as.data.frame(x$appendix_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_section_summary) && nrow(x$appendix_section_summary) > 0L) {
    cat("\nAppendix section summary\n")
    print(round_numeric_df(as.data.frame(x$appendix_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_summary) && nrow(x$selection_summary) > 0L) {
    cat("\nSelection summary\n")
    print(round_numeric_df(as.data.frame(x$selection_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_table_summary) && nrow(x$selection_table_summary) > 0L) {
    cat("\nSelection table summary\n")
    print(round_numeric_df(as.data.frame(x$selection_table_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_table_preset_summary) && nrow(x$selection_table_preset_summary) > 0L) {
    cat("\nSelection table preset summary\n")
    print(round_numeric_df(as.data.frame(x$selection_table_preset_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_table_summary) && nrow(x$selection_handoff_table_summary) > 0L) {
    cat("\nSelection handoff table summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_table_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_preset_summary) && nrow(x$selection_handoff_preset_summary) > 0L) {
    cat("\nSelection handoff preset summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_preset_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_summary) && nrow(x$selection_handoff_summary) > 0L) {
    cat("\nSelection handoff summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_bundle_summary) && nrow(x$selection_handoff_bundle_summary) > 0L) {
    cat("\nSelection handoff bundle summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_bundle_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_role_summary) && nrow(x$selection_handoff_role_summary) > 0L) {
    cat("\nSelection handoff role summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_role_section_summary) && nrow(x$selection_handoff_role_section_summary) > 0L) {
    cat("\nSelection handoff role-section summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_role_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_role_summary) && nrow(x$selection_role_summary) > 0L) {
    cat("\nSelection role summary\n")
    print(round_numeric_df(as.data.frame(x$selection_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_section_summary) && nrow(x$selection_section_summary) > 0L) {
    cat("\nSelection section summary\n")
    print(round_numeric_df(as.data.frame(x$selection_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_catalog) && nrow(x$selection_catalog) > 0L) {
    cat("\nSelection catalog\n")
    print(round_numeric_df(as.data.frame(x$selection_catalog), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$reporting_map) && nrow(x$reporting_map) > 0L) {
    cat("\nReporting map\n")
    print(as.data.frame(x$reporting_map), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

resolve_selection_plot_measure <- function(tbl,
                                           type,
                                           selection_value = c("count", "fraction")) {
  selection_value <- match.arg(selection_value)
  tbl <- as.data.frame(tbl %||% data.frame(), stringsAsFactors = FALSE)

  if (identical(type, "selection_tables")) {
    if (identical(selection_value, "fraction")) {
      stop("`selection_value = \"fraction\"` is not available for `type = \"selection_tables\"`; this surface only exposes table row counts.", call. = FALSE)
    }
    return(list(
      values = suppressWarnings(as.numeric(tbl$Rows)),
      ylab = "Rows",
      legend_label = "Rows",
      selection_value = "count"
    ))
  }

  if (identical(type, "selection_bundles")) {
    if (identical(selection_value, "count")) {
      return(list(
        values = suppressWarnings(as.numeric(tbl$TablesSelected)),
        ylab = "Tables",
        legend_label = "Tables selected",
        selection_value = "count"
      ))
    }
    if (!"SelectionFraction" %in% names(tbl)) {
      stop("`selection_value = \"fraction\"` is not available because `SelectionFraction` is missing from this surface.", call. = FALSE)
    }
    return(list(
      values = suppressWarnings(as.numeric(tbl$SelectionFraction)),
      ylab = "Selection fraction",
      legend_label = "Selection fraction",
      selection_value = "fraction"
    ))
  }

  if (identical(selection_value, "count")) {
    count_col <- if (type %in% c("selection_handoff_presets", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections")) {
      "PlotReadyTables"
    } else {
      "Tables"
    }
    count_label <- if (identical(count_col, "PlotReadyTables")) "Plot-ready tables" else "Tables"
    if (!count_col %in% names(tbl)) {
      stop("`selection_value = \"count\"` is not available because `", count_col, "` is missing from this surface.", call. = FALSE)
    }
    return(list(
      values = suppressWarnings(as.numeric(tbl[[count_col]])),
      ylab = count_label,
      legend_label = count_label,
      selection_value = "count"
    ))
  }

  if (!"PlotReadyFraction" %in% names(tbl)) {
    stop("`selection_value = \"fraction\"` is not available because `PlotReadyFraction` is missing from this surface.", call. = FALSE)
  }
  list(
    values = suppressWarnings(as.numeric(tbl$PlotReadyFraction)),
    ylab = "Plot-ready fraction",
    legend_label = "Plot-ready fraction",
    selection_value = "fraction"
  )
}

summary_table_bundle_filter_selection_tables <- function(tbl, appendix_preset) {
  tbl <- as.data.frame(tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    return(tbl)
  }
  if ("Preset" %in% names(tbl)) {
    keep <- as.character(tbl$Preset %||% "") %in% appendix_preset
    return(tbl[keep, , drop = FALSE])
  }
  if (!"Presets" %in% names(tbl)) {
    return(tbl)
  }
  keep <- vapply(as.character(tbl$Presets %||% ""), function(x) {
    tokens <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    any(tokens %in% appendix_preset)
  }, logical(1))
  tbl[keep, , drop = FALSE]
}

#' Plot a summary-table bundle for manuscript QC
#'
#' @param x Output from [build_summary_table_bundle()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"table_rows"` for returned-table sizes,
#'   `"role_tables"` for returned-table counts by reporting role,
#'   `"appendix_roles"` for returned-table counts by reporting role under the
#'   bundle's appendix-routing contract,
#'   `"appendix_sections"` for returned-table counts by manuscript-facing
#'   appendix section,
#'   `"appendix_presets"` for conservative appendix-preset counts,
#'   `"selection_handoff_presets"` for workflow-only preset-level appendix
#'   handoff counts,
#'   `"selection_tables"` / `"selection_handoff"` /
#'   `"selection_handoff_bundles"` /
#'   `"selection_handoff_roles"` / `"selection_bundles"` /
#'   `"selection_roles"` / `"selection_sections"` for workflow-only appendix
#'   selection surfaces when present in the bundle,
#'   `"numeric_profile"` for column means from a selected numeric table, or
#'   `"first_numeric"` for the distribution of the first numeric column in a
#'   selected table.
#' @param selection_value For `selection_*` plot types, whether to plot exact
#'   counts (`"count"`) or the corresponding exact fraction (`"fraction"`)
#'   when that surface exposes one.
#' @param appendix_preset Appendix preset used for `selection_*` plot types.
#' @param which Optional table selector used for numeric plot types.
#' @param main Optional title override.
#' @param palette Optional named color overrides.
#' @param label_angle Axis-label rotation angle for bar-type plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This helper keeps summary-bundle plotting conservative. It either visualizes
#' the bundle's own bundle-level indexes (`"table_rows"`, `"role_tables"`,
#' `"appendix_roles"`, `"appendix_sections"`, `"appendix_presets"`) or routes a
#' selected table through [apa_table()] and [plot.apa_table()] for numeric QC.
#'
#' @section Interpreting output:
#' - `"table_rows"`: compares returned table sizes to show where reporting mass sits.
#' - `"role_tables"`: shows how many returned tables belong to each reporting role.
#' - `"appendix_roles"`: shows how returned tables contribute to conservative
#'   appendix routing by reporting role.
#' - `"appendix_sections"`: shows how returned tables are distributed across
#'   methods/results/diagnostics/reporting sections.
#' - `"appendix_presets"`: shows how many tables the current bundle contributes
#'   to the conservative appendix presets.
#' - `"selection_handoff_presets"`: shows plot-ready appendix handoff counts by
#'   preset for workflow-only appendix routing surfaces in the bundle.
#' - `"selection_tables"` / `"selection_handoff"` /
#'   `"selection_handoff_bundles"` /
#'   `"selection_handoff_roles"` / `"selection_handoff_role_sections"` /
#'   `"selection_bundles"` /
#'   `"selection_roles"` / `"selection_sections"`: show workflow-only appendix
#'   selection surfaces already materialized inside the bundle.
#' - `"numeric_profile"` / `"first_numeric"`: reuse the same numeric QC logic as
#'   [plot.apa_table()] but start from a summary-table bundle.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [build_summary_table_bundle()], [apa_table()], [plot.apa_table()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' bundle <- build_summary_table_bundle(fit)
#' plot(bundle, draw = FALSE)
#' plot(bundle, type = "numeric_profile", which = "facet_overview", draw = FALSE)
#' }
#' @export
plot.mfrm_summary_table_bundle <- function(x,
                                           y = NULL,
                                           type = c("table_rows", "role_tables", "appendix_roles", "appendix_sections", "appendix_presets", "selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections", "numeric_profile", "first_numeric"),
                                           which = NULL,
                                           selection_value = c("count", "fraction"),
                                           appendix_preset = c("recommended", "compact", "all", "methods", "results", "diagnostics", "reporting"),
                                           main = NULL,
                                           palette = NULL,
                                           label_angle = 45,
                                           draw = TRUE,
                                           ...) {
  type <- match.arg(
    tolower(as.character(type[1])),
    c("table_rows", "role_tables", "appendix_roles", "appendix_sections", "appendix_presets", "selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections", "numeric_profile", "first_numeric")
  )
  appendix_preset <- match.arg(
    tolower(as.character(appendix_preset[1])),
    c("recommended", "compact", "all", "methods", "results", "diagnostics", "reporting")
  )
  selection_value <- match.arg(selection_value)

  if (type == "table_rows") {
    idx <- as.data.frame(x$table_index %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(idx) == 0L || !"Table" %in% names(idx) || !"Rows" %in% names(idx)) {
      stop("`x$table_index` does not contain plottable table-row information.")
    }
    rows <- suppressWarnings(as.numeric(idx$Rows))
    keep <- is.finite(rows)
    if (!any(keep)) {
      stop("`x$table_index` does not contain finite row counts to plot.")
    }
    rows <- rows[keep]
    labels <- as.character(idx$Table[keep])
    ord <- order(rows, decreasing = TRUE, na.last = NA)
    rows <- rows[ord]
    labels <- labels[ord]
    pal <- resolve_palette(
      palette = palette,
      defaults = c(table_rows = "#6a4c93", grid = "#ececec")
    )
    plot_title <- if (is.null(main)) "Summary bundle table sizes" else as.character(main[1])
    if (isTRUE(draw)) {
      barplot_rot45(
        height = rows,
        labels = labels,
        col = pal["table_rows"],
        main = plot_title,
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    return(invisible(new_mfrm_plot_data(
      "summary_table_bundle",
      list(
        plot = "table_rows",
        table = labels,
        rows = rows,
        title = plot_title,
        subtitle = "Returned summary tables ranked by row count",
        legend = new_plot_legend("Table rows", "summary_table", "bar", pal["table_rows"]),
        reference_lines = new_reference_lines("h", 0, "Zero-row reference", "dashed", "reference")
      )
    )))
  }

  if (type == "role_tables") {
    idx <- as.data.frame(x$table_index %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(idx) == 0L || !"Role" %in% names(idx)) {
      stop("`x$table_index` does not contain plottable role information.")
    }
    roles <- as.character(idx$Role)
    roles <- roles[nzchar(roles)]
    if (length(roles) == 0L) {
      stop("`x$table_index` does not contain non-empty role labels to plot.")
    }
    counts <- sort(table(roles), decreasing = TRUE)
    labels <- names(counts)
    values <- as.numeric(counts)
    pal <- resolve_palette(
      palette = palette,
      defaults = c(role_tables = "#3a7ca5", grid = "#ececec")
    )
    plot_title <- if (is.null(main)) "Summary bundle role coverage" else as.character(main[1])
    if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = pal["role_tables"],
        main = plot_title,
        ylab = "Tables",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    return(invisible(new_mfrm_plot_data(
      "summary_table_bundle",
      list(
        plot = "role_tables",
        role = labels,
        tables = values,
        title = plot_title,
        subtitle = "Returned summary tables grouped by reporting role",
        legend = new_plot_legend("Role table count", "summary_table", "bar", pal["role_tables"]),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (type == "appendix_roles") {
    role_tbl <- summary_table_bundle_appendix_role_summary(summary_table_bundle_catalog(x))
    if (nrow(role_tbl) == 0L || !"Role" %in% names(role_tbl) ||
        !"RecommendedTables" %in% names(role_tbl) ||
        !"CompactTables" %in% names(role_tbl)) {
      stop("`x` does not contain plottable appendix-role information.")
    }
    labels <- as.character(role_tbl$Role)
    recommended <- suppressWarnings(as.numeric(role_tbl$RecommendedTables))
    compact <- suppressWarnings(as.numeric(role_tbl$CompactTables))
    keep <- nzchar(labels) & is.finite(recommended) & is.finite(compact)
    if (!any(keep)) {
      stop("`x` does not contain finite appendix-role table counts to plot.")
    }
    labels <- labels[keep]
    recommended <- recommended[keep]
    compact <- compact[keep]
    values <- rbind(Recommended = recommended, Compact = compact)
    pal <- resolve_palette(
      palette = palette,
      defaults = c(appendix_role_recommended = "#2a9d8f", appendix_role_compact = "#8d99ae", grid = "#ececec")
    )
    plot_title <- if (is.null(main)) "Summary bundle appendix roles" else as.character(main[1])
    if (isTRUE(draw)) {
      graphics::barplot(
        height = values,
        beside = TRUE,
        names.arg = labels,
        col = c(pal["appendix_role_recommended"], pal["appendix_role_compact"]),
        main = plot_title,
        ylab = "Tables",
        las = 2,
        cex.names = 0.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
      graphics::legend(
        "topright",
        legend = c("Recommended", "Compact"),
        fill = c(pal["appendix_role_recommended"], pal["appendix_role_compact"]),
        bty = "n"
      )
    }
    return(invisible(new_mfrm_plot_data(
      "summary_table_bundle",
      list(
        plot = "appendix_roles",
        role = labels,
        recommended_tables = recommended,
        compact_tables = compact,
        title = plot_title,
        subtitle = "Appendix-routed table counts by reporting role",
        legend = list(
          new_plot_legend("Recommended", "summary_table", "bar", pal["appendix_role_recommended"]),
          new_plot_legend("Compact", "summary_table", "bar", pal["appendix_role_compact"])
        ),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (type == "appendix_sections") {
    section_tbl <- summary_table_bundle_appendix_section_summary(summary_table_bundle_catalog(x))
    if (nrow(section_tbl) == 0L || !"AppendixSection" %in% names(section_tbl) || !"Tables" %in% names(section_tbl)) {
      stop("`x` does not contain plottable appendix-section information.")
    }
    labels <- as.character(section_tbl$AppendixSection)
    values <- suppressWarnings(as.numeric(section_tbl$Tables))
    keep <- is.finite(values) & nzchar(labels)
    if (!any(keep)) {
      stop("`x` does not contain finite appendix-section table counts to plot.")
    }
    labels <- labels[keep]
    values <- values[keep]
    pal <- resolve_palette(
      palette = palette,
      defaults = c(appendix_sections = "#457b9d", grid = "#ececec")
    )
    plot_title <- if (is.null(main)) "Summary bundle appendix sections" else as.character(main[1])
    if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = pal["appendix_sections"],
        main = plot_title,
        ylab = "Tables",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    return(invisible(new_mfrm_plot_data(
      "summary_table_bundle",
      list(
        plot = "appendix_sections",
        appendix_section = labels,
        tables = values,
        title = plot_title,
        subtitle = "Returned summary tables grouped by manuscript appendix section",
        legend = new_plot_legend("Appendix section count", "summary_table", "bar", pal["appendix_sections"]),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (type == "appendix_presets") {
    preset_tbl <- summary_table_bundle_appendix_presets(summary_table_bundle_catalog(x))
    if (nrow(preset_tbl) == 0L || !"Preset" %in% names(preset_tbl) || !"Tables" %in% names(preset_tbl)) {
      stop("`x` does not contain plottable appendix-preset information.")
    }
    labels <- as.character(preset_tbl$Preset)
    values <- suppressWarnings(as.numeric(preset_tbl$Tables))
    keep <- is.finite(values) & nzchar(labels)
    if (!any(keep)) {
      stop("`x` does not contain finite appendix-preset table counts to plot.")
    }
    labels <- labels[keep]
    values <- values[keep]
    pal <- resolve_palette(
      palette = palette,
      defaults = c(appendix_presets = "#2a9d8f", grid = "#ececec")
    )
    plot_title <- if (is.null(main)) "Summary bundle appendix presets" else as.character(main[1])
    if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = pal["appendix_presets"],
        main = plot_title,
        ylab = "Tables",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    return(invisible(new_mfrm_plot_data(
      "summary_table_bundle",
      list(
        plot = "appendix_presets",
        preset = labels,
        tables = values,
        title = plot_title,
        subtitle = "Current bundle size under conservative appendix presets",
        legend = new_plot_legend("Appendix preset count", "summary_table", "bar", pal["appendix_presets"]),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (type %in% c("selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections")) {
    sx <- summary(x)
    measure <- NULL
    if (type == "selection_handoff_presets") {
      tbl <- as.data.frame(sx$selection_handoff_preset_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Preset", "PlotReadyTables") %in% names(tbl))) {
        stop("`x` does not contain appendix handoff-preset rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$Preset)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_handoff_presets = "#f4a261", grid = "#ececec"))
      plot_name <- "selection_handoff_presets"
      plot_title <- if (is.null(main)) paste0("Summary bundle handoff presets (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Preset-level plot-ready appendix handoff for `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_tables") {
      tbl <- as.data.frame(sx$selection_table_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- summary_table_bundle_filter_selection_tables(tbl, appendix_preset = appendix_preset)
      if (nrow(tbl) == 0L || !all(c("Table", "Rows") %in% names(tbl))) {
        stop("`x` does not contain appendix table-selection rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$Table)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_tables = "#e76f51", grid = "#ececec"))
      plot_name <- "selection_tables"
      plot_title <- if (is.null(main)) paste0("Summary bundle selection tables (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Selected appendix tables for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_handoff") {
      tbl <- as.data.frame(sx$selection_handoff_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "PlotReadyTables") %in% names(tbl))) {
        stop("`x` does not contain appendix handoff rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$AppendixSection)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_handoff = "#ff9f1c", grid = "#ececec"))
      plot_name <- "selection_handoff"
      plot_title <- if (is.null(main)) paste0("Summary bundle selection handoff (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Plot-ready appendix handoff by section for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_handoff_bundles") {
      tbl <- as.data.frame(sx$selection_handoff_bundle_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Bundle", "PlotReadyTables") %in% names(tbl))) {
        stop("`x` does not contain appendix handoff-bundle rows for preset `", appendix_preset, "`.")
      }
      labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Bundle))
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_handoff_bundles = "#5c677d", grid = "#ececec"))
      plot_name <- "selection_handoff_bundles"
      plot_title <- if (is.null(main)) paste0("Summary bundle handoff bundles (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Plot-ready appendix handoff by section and bundle for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_handoff_roles") {
      tbl <- as.data.frame(sx$selection_handoff_role_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Role", "PlotReadyTables") %in% names(tbl))) {
        stop("`x` does not contain appendix handoff-role rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$Role)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_handoff_roles = "#9c6644", grid = "#ececec"))
      plot_name <- "selection_handoff_roles"
      plot_title <- if (is.null(main)) paste0("Summary bundle handoff roles (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Plot-ready appendix handoff by role for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_handoff_role_sections") {
      tbl <- as.data.frame(sx$selection_handoff_role_section_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Role", "PlotReadyTables") %in% names(tbl))) {
        stop("`x` does not contain appendix handoff role-section rows for preset `", appendix_preset, "`.")
      }
      labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Role))
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_handoff_role_sections = "#7f5539", grid = "#ececec"))
      plot_name <- "selection_handoff_role_sections"
      plot_title <- if (is.null(main)) paste0("Summary bundle handoff role-sections (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Plot-ready appendix handoff by section and role for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_bundles") {
      tbl <- as.data.frame(sx$selection_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Bundle", "TablesSelected") %in% names(tbl))) {
        stop("`x` does not contain appendix bundle-selection rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$Bundle)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_bundles = "#54a24b", grid = "#ececec"))
      plot_name <- "selection_bundles"
      plot_title <- if (is.null(main)) paste0("Summary bundle appendix bundles (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Appendix tables by source bundle for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else if (type == "selection_roles") {
      tbl <- as.data.frame(sx$selection_role_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Role", "Tables") %in% names(tbl))) {
        stop("`x` does not contain appendix role-selection rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$Role)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_roles = "#b279a2", grid = "#ececec"))
      plot_name <- "selection_roles"
      plot_title <- if (is.null(main)) paste0("Summary bundle appendix roles (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Selected appendix roles for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    } else {
      tbl <- as.data.frame(sx$selection_section_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Tables") %in% names(tbl))) {
        stop("`x` does not contain appendix section-selection rows for preset `", appendix_preset, "`.")
      }
      labels <- as.character(tbl$AppendixSection)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      pal <- resolve_palette(palette = palette, defaults = c(selection_sections = "#2a9d8f", grid = "#ececec"))
      plot_name <- "selection_sections"
      plot_title <- if (is.null(main)) paste0("Summary bundle appendix sections (", appendix_preset, ")") else as.character(main[1])
      subtitle <- paste0("Selected appendix sections for preset `", appendix_preset, "`")
      ylab <- measure$ylab
      legend_label <- measure$legend_label
    }

    keep <- is.finite(values) & nzchar(labels)
    if (!any(keep)) {
      stop("`x` does not contain finite appendix-selection values for plot type `", type, "`.")
    }
    tbl <- tbl[keep, , drop = FALSE]
    labels <- labels[keep]
    values <- values[keep]

    if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = pal[plot_name],
        main = plot_title,
        ylab = ylab,
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    return(invisible(new_mfrm_plot_data(
      "summary_table_bundle",
      list(
        plot = plot_name,
        selection_value = measure$selection_value,
        appendix_preset = appendix_preset,
        table = tbl,
        title = plot_title,
        subtitle = subtitle,
        legend = new_plot_legend(legend_label, "summary_table", "bar", pal[plot_name]),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (is.null(which) || !nzchar(as.character(which[1] %||% ""))) {
    which <- summary_table_bundle_first_numeric_table(x)
  }
  if (is.null(which)) {
    stop("No numeric summary table is available for plot type `", type, "`.", call. = FALSE)
  }

  apa_obj <- apa_table(x, which = which)
  apa_plot <- plot.apa_table(
    apa_obj,
    type = type,
    main = main,
    palette = palette,
    label_angle = label_angle,
    draw = draw,
    ...
  )
  payload <- apa_plot$data
  payload$source_table <- as.character(which[1])
  payload$source_bundle_class <- as.character(x$source_class %||% "mfrm_summary_table_bundle")
  invisible(new_mfrm_plot_data("summary_table_bundle", payload))
}

resolve_summary_bundle_table_selection <- function(bundle, which = NULL) {
  if (!inherits(bundle, "mfrm_summary_table_bundle")) {
    stop("`bundle` must be an mfrm_summary_table_bundle object.", call. = FALSE)
  }
  available <- names(bundle$tables %||% list())
  if (length(available) == 0L) {
    stop("`bundle` does not contain any tables.", call. = FALSE)
  }
  if (is.null(which) || !nzchar(as.character(which[1] %||% ""))) {
    which <- if ("overview" %in% available) "overview" else available[1]
  } else {
    which <- as.character(which[1])
  }
  if (!which %in% available) {
    stop(
      "Requested `which` not found in summary table bundle. Available tables: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  idx <- as.data.frame(bundle$table_index %||% data.frame(), stringsAsFactors = FALSE)
  idx_row <- if (nrow(idx) > 0 && "Table" %in% names(idx)) {
    idx[idx$Table %in% which, , drop = FALSE]
  } else {
    data.frame()
  }
  list(
    which = which,
    table = as.data.frame(bundle$tables[[which]], stringsAsFactors = FALSE),
    index_row = idx_row
  )
}

#' Build APA-style table output using base R structures
#'
#' @param x A data.frame, `mfrm_fit`, `summary()` output supported by
#'   [build_summary_table_bundle()], an `mfrm_summary_table_bundle`, diagnostics
#'   list, or bias-result list.
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
#' - For `summary()` outputs or `mfrm_summary_table_bundle`:
#'   names listed in `build_summary_table_bundle(x)$table_index`
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
#' fit_bundle <- build_summary_table_bundle(summary(fit))
#' tbl_from_summary <- apa_table(fit_bundle, which = "facet_overview")
#' summary(tbl)
#' p <- plot(tbl, draw = FALSE)
#' p_facets <- plot(tbl_facets, type = "numeric_profile", draw = FALSE)
#' p$data$plot
#' p_facets$data$plot
#' if (interactive()) {
#'   plot(
#'     tbl,
#'     type = "numeric_profile",
#'     main = "APA Table Numeric Profile (Customized)",
#'     palette = c(numeric_profile = "#2b8cbe", grid = "#d9d9d9"),
#'     label_angle = 45
#'   )
#' }
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

  summary_bundle_classes <- summary_table_bundle_supported_summary_classes()

  if (inherits(x, "mfrm_summary_table_bundle")) {
    source_type <- "mfrm_summary_table_bundle"
    selected <- resolve_summary_bundle_table_selection(x, which = which)
    resolved_which <- selected$which
    table_out <- selected$table
    idx_row <- selected$index_row
    if (is.null(caption) && nrow(idx_row) > 0 && "Description" %in% names(idx_row)) {
      caption <- as.character(idx_row$Description[1])
    }
    if (is.null(note) && length(x$notes %||% character(0)) > 0L) {
      note <- paste(as.character(x$notes), collapse = " ")
    }
  } else if (inherits(x, summary_bundle_classes)) {
    source_type <- class(x)[1]
    bundle <- build_summary_table_bundle(x, include_empty = TRUE)
    selected <- resolve_summary_bundle_table_selection(bundle, which = which)
    resolved_which <- selected$which
    table_out <- selected$table
    idx_row <- selected$index_row
    if (is.null(caption) && nrow(idx_row) > 0 && "Description" %in% names(idx_row)) {
      caption <- as.character(idx_row$Description[1])
    }
    if (is.null(note) && length(bundle$notes %||% character(0)) > 0L) {
      note <- paste(as.character(bundle$notes), collapse = " ")
    }
  } else if (is.data.frame(x)) {
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
    stop("`x` must be a data.frame, mfrm_fit, supported summary/table-bundle object, or named list.")
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
    validated <- validate_apa_builder_inputs(
      fit = x,
      diagnostics = diag_for_contract,
      bias_results = bias_results,
      context = context,
      helper = "apa_table()"
    )
    x <- validated$fit
    diag_for_contract <- validated$diagnostics
    bias_results <- validated$bias_results
    context <- validated$context
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
#' For bounded `GPCM`, this helper is intentionally unavailable. Use
#' [reporting_checklist()], [plot_qc_dashboard()], the residual/category table
#' helpers, and [compute_information()] / [plot_information()] instead.
#'
#' @section Interpreting output:
#' - `warning_map`: rule-triggered warning text by visual key.
#' - `summary_map`: descriptive narrative text by visual key.
#' - strict marginal keys appear when `diagnose_mfrm(..., diagnostic_mode = "both")`
#'   supplies latent-integrated first-order and pairwise screening summaries.
#' - `warning_counts` / `summary_counts`: message-count tables for QA checks.
#' - `plot_payloads`: ready-to-reuse `mfrm_plot_data` payloads for the bundle's
#'   own comparison/count plots and, when step estimates are available, the
#'   exploratory `category_probability_surface` payload from
#'   `plot(fit, type = "ccc_surface", draw = FALSE)`. The surface payload
#'   carries `category_support`, `interpretation_guide`, and `reporting_policy`
#'   tables for zero-frequency category and reporting-boundary checks.
#' - `public_plot_routes`: draw-free helper routes for the dedicated public plot
#'   functions behind each visual family.
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
#' - `plot_payloads`: reusable draw-free payloads for `comparison`,
#'   `warning_counts`, `summary_counts`, and optionally
#'   `category_probability_surface`
#' - `public_plot_routes`: public helper / draw-free route map for follow-up
#' - `crosswalk`: FACETS-reference mapping for main visual keys
#' - `branch`, `style`, `threshold_profile`: branch metadata
#'
#' @seealso [mfrm_threshold_profiles()], [build_apa_outputs()],
#'   [plot_marginal_fit()], [plot_marginal_pairwise()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", model = "RSM", maxit = 200
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "both", diagnostic_mode = "both")
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
#' vis$plot_payloads$comparison$data$plot
#' vis$public_plot_routes[, c("Visual", "PlotHelper", "DrawFreeRoute")]
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
  stop_if_gpcm_out_of_scope(fit, "build_visual_summaries()")
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
      "strict_marginal_fit",
      "strict_pairwise_local_dependence",
      "residual_pca_overall",
      "residual_pca_by_facet",
      "category_probability_surface"
    ),
    FACETS = c(
      "Table 4 / Table 10",
      "Table 12",
      "Table 9",
      "Inter-rater outputs",
      "Facet fixed/random chi-square",
      "No direct FACETS equivalent (package-native strict marginal screen)",
      "No direct FACETS equivalent (package-native strict pairwise screen)",
      "Residual PCA (overall)",
      "Residual PCA (by facet)",
      "No direct FACETS equivalent (exploratory category-probability surface payload)"
    )
  )

  out <- list(
    warning_map = warning_map,
    summary_map = summary_map,
    warning_counts = to_count_table(warning_map),
    summary_counts = to_count_table(summary_map),
    plot_payloads = NULL,
    public_plot_routes = NULL,
    crosswalk = crosswalk,
    branch = branch,
    style = style,
    threshold_profile = as.character(threshold_profile[1])
  )
  out$plot_payloads <- build_visual_plot_payloads(out, fit = fit)
  out$public_plot_routes <- build_visual_plot_route_table()
  out <- as_mfrm_bundle(out, "mfrm_visual_summaries")
  class(out) <- unique(c(paste0("mfrm_visual_summaries_", branch), class(out)))
  out
}

build_visual_plot_route_table <- function() {
  tibble::tibble(
    Visual = c(
      "comparison",
      "warning_counts",
      "summary_counts",
      "unexpected",
      "fair_average",
      "displacement",
      "interrater",
      "facets_chisq",
      "strict_marginal_fit",
      "strict_pairwise_local_dependence",
      "residual_pca_overall",
      "residual_pca_by_facet",
      "category_probability_surface"
    ),
    PlotHelper = c(
      "plot.mfrm_bundle()",
      "plot.mfrm_bundle()",
      "plot.mfrm_bundle()",
      "plot_unexpected()",
      "plot_fair_average()",
      "plot_displacement()",
      "plot_interrater_agreement()",
      "plot_facets_chisq()",
      "plot_marginal_fit()",
      "plot_marginal_pairwise()",
      "plot_residual_pca()",
      "plot_residual_pca()",
      "plot.mfrm_fit()"
    ),
    DrawFreeRoute = c(
      "plot(vis, type = \"comparison\", draw = FALSE)",
      "plot(vis, type = \"warning_counts\", draw = FALSE)",
      "plot(vis, type = \"summary_counts\", draw = FALSE)",
      "plot_unexpected(unexpected_response_table(fit, diagnostics = diagnostics), draw = FALSE)",
      "plot_fair_average(fair_average_table(fit, diagnostics = diagnostics), draw = FALSE)",
      "plot_displacement(displacement_table(fit, diagnostics = diagnostics), draw = FALSE)",
      "plot_interrater_agreement(interrater_agreement_table(fit, diagnostics = diagnostics), draw = FALSE)",
      "plot_facets_chisq(facets_chisq_table(fit, diagnostics = diagnostics), draw = FALSE)",
      "plot_marginal_fit(diagnostics, draw = FALSE)",
      "plot_marginal_pairwise(diagnostics, draw = FALSE)",
      "plot_residual_pca(analyze_residual_pca(diagnostics, mode = \"overall\"), mode = \"overall\", plot_type = \"scree\", draw = FALSE)",
      "plot_residual_pca(analyze_residual_pca(diagnostics, mode = \"both\"), mode = \"facet\", facet = \"<facet>\", plot_type = \"loadings\", draw = FALSE)",
      "plot(fit, type = \"ccc_surface\", draw = FALSE)"
    ),
    PlotReturnClass = rep("mfrm_plot_data", 13L),
    Scope = c(
      "bundle overview",
      "bundle overview",
      "bundle overview",
      "unexpected-response follow-up",
      "fair-average follow-up",
      "displacement follow-up",
      "inter-rater follow-up",
      "facet chi-square follow-up",
      "strict marginal follow-up",
      "strict pairwise follow-up",
      "overall residual-structure follow-up",
      "facet-level residual-structure follow-up",
      "exploratory category-probability surface handoff"
    )
  )
}

build_visual_plot_payloads <- function(x, fit = NULL) {
  payloads <- list(
    comparison = plot_visual_summaries_bundle(x, plot_type = "comparison", draw = FALSE),
    warning_counts = plot_visual_summaries_bundle(x, plot_type = "warning_counts", draw = FALSE),
    summary_counts = plot_visual_summaries_bundle(x, plot_type = "summary_counts", draw = FALSE)
  )
  if (inherits(fit, "mfrm_fit")) {
    surface <- tryCatch(
      plot(fit, type = "ccc_surface", draw = FALSE),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (inherits(surface, "mfrm_plot_data")) {
      payloads$category_probability_surface <- surface
    }
  }
  payloads
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
#' 2. Inspect `summary(contract_audit)` and `missing_preview`.
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
#' contract_audit <- facets_parity_report(fit, diagnostics = diag, branch = "facets")
#' summary(contract_audit)
#' p <- plot(contract_audit, draw = FALSE)
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
  stop_if_gpcm_out_of_scope(fit, "facets_parity_report()")
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
#' For bounded `GPCM`, this pipeline is intentionally unavailable because the
#' current validated route stops before bundled pass/warn/fail synthesis for
#' the free-discrimination branch.
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
  stop_if_gpcm_out_of_scope(fit, "run_qc_pipeline()")

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
