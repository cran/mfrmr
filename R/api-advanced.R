# Advanced analysis functions: differential functioning, design-weighted precision curves,
# unified Wright map.

# ============================================================================
# B. DFF Analysis
# ============================================================================

functioning_label_for_facet <- function(facet) {
  facet_lower <- tolower(as.character(facet)[1])
  if (facet_lower %in% c("rater", "raters")) {
    return("DRF")
  }
  if (facet_lower %in% c("prompt", "prompts", "task", "tasks")) {
    return("DPF")
  }
  if (facet_lower %in% c("item", "items", "criterion", "criteria")) {
    return("DIF")
  }
  "DFF"
}

welch_satterthwaite_df <- function(components, dfs) {
  components <- suppressWarnings(as.numeric(components))
  dfs <- suppressWarnings(as.numeric(dfs))
  ok <- is.finite(components) & components > 0 & is.finite(dfs) & dfs > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  num <- sum(components[ok])^2
  den <- sum((components[ok]^2) / dfs[ok])
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  num / den
}

resolve_dff_subgroup_precision <- function(sub_fit, sub_diag = NULL, diagnostics_error = NULL) {
  precision_profile <- as.data.frame(sub_diag$precision_profile %||% data.frame(), stringsAsFactors = FALSE)
  method <- as.character(sub_fit$summary$Method[1] %||% sub_fit$config$method %||% NA_character_)
  method <- ifelse(identical(method, "JMLE"), "JML", method)
  converged <- if (nrow(precision_profile) > 0 && "Converged" %in% names(precision_profile)) {
    isTRUE(precision_profile$Converged[1])
  } else {
    isTRUE(as.logical(sub_fit$summary$Converged[1] %||% FALSE))
  }

  precision_tier <- if (!is.null(diagnostics_error) && nzchar(diagnostics_error)) {
    "diagnostics_unavailable"
  } else if (nrow(precision_profile) > 0 && "PrecisionTier" %in% names(precision_profile)) {
    as.character(precision_profile$PrecisionTier[1] %||% NA_character_)
  } else if (identical(method, "MML")) {
    "hybrid"
  } else {
    "exploratory"
  }

  supports_formal <- if (nrow(precision_profile) > 0 && "SupportsFormalInference" %in% names(precision_profile)) {
    isTRUE(precision_profile$SupportsFormalInference[1])
  } else {
    FALSE
  }

  list(
    method = method,
    converged = converged,
    precision_tier = precision_tier,
    supports_formal = supports_formal
  )
}

annotate_dff_table <- function(tbl, method) {
  if (nrow(tbl) == 0) {
    tbl$EffectMetric <- character(0)
    tbl$ContrastBasis <- character(0)
    tbl$SEBasis <- character(0)
    tbl$StatisticLabel <- character(0)
    tbl$ProbabilityMetric <- character(0)
    tbl$DFBasis <- character(0)
    tbl$ClassificationSystem <- character(0)
    tbl$Classification <- character(0)
    tbl$ETS <- character(0)
    tbl$ReportingUse <- character(0)
    tbl$PrimaryReportingEligible <- logical(0)
    return(tbl)
  }

  if (identical(method, "refit")) {
    abs_diff <- abs(tbl$Contrast)
    linked <- if ("ContrastComparable" %in% names(tbl)) {
      as.logical(tbl$ContrastComparable)
    } else {
      rep(TRUE, nrow(tbl))
    }
    formal <- if ("FormalInferenceEligible" %in% names(tbl)) {
      as.logical(tbl$FormalInferenceEligible)
    } else {
      rep(TRUE, nrow(tbl))
    }
    comparable <- linked & formal
    linked_descriptive <- linked & !comparable
    ets <- ifelse(
      comparable & is.finite(abs_diff),
      ifelse(abs_diff < 0.43, "A", ifelse(abs_diff < 0.64, "B", "C")),
      NA_character_
    )
    tbl$ContrastBasis <- ifelse(
      linked,
      "linked subgroup facet-measure difference",
      "descriptive subgroup facet-measure difference"
    )
    tbl$SEBasis <- ifelse(
      comparable,
      "joint subgroup-calibration standard error",
      ifelse(
        linked_descriptive,
        "not reported without model-based subgroup precision",
        "not reported without common-scale linking"
      )
    )
    tbl$StatisticLabel <- ifelse(
      comparable,
      "Welch t",
      ifelse(linked_descriptive, "linked descriptive contrast", "descriptive contrast")
    )
    tbl$ProbabilityMetric <- ifelse(comparable, "Welch t tail area", "not reported")
    tbl$DFBasis <- ifelse(comparable, "Welch-Satterthwaite approximation", "not reported")
    tbl$EffectMetric <- ifelse(
      comparable,
      "linked_logit_difference",
      ifelse(linked_descriptive, "linked_descriptive_logit_difference", "descriptive_refit_difference")
    )
    tbl$ClassificationSystem <- ifelse(comparable, "ETS", "descriptive")
    tbl$Classification <- dplyr::case_when(
      ets == "A" ~ "A (Negligible)",
      ets == "B" ~ "B (Moderate)",
      ets == "C" ~ "C (Large)",
      linked_descriptive %in% TRUE & is.finite(tbl$Contrast) ~ "Linked contrast (screening only)",
      linked %in% FALSE & is.finite(tbl$Contrast) ~ "Unclassified (insufficient linking)",
      TRUE ~ NA_character_
    )
    tbl$ETS <- ets
    existing_primary <- tbl$PrimaryReportingEligible %||% rep(NA, nrow(tbl))
    existing_use <- tbl$ReportingUse %||% rep(NA_character_, nrow(tbl))
    tbl$PrimaryReportingEligible <- dplyr::coalesce(as.logical(existing_primary), comparable)
    tbl$ReportingUse <- dplyr::coalesce(
      as.character(existing_use),
      dplyr::case_when(
        comparable ~ "primary_reporting",
        linked_descriptive %in% TRUE ~ "review_before_reporting",
        TRUE ~ "screening_only"
      )
    )
    return(tbl)
  }

  sig <- dplyr::if_else(
    is.finite(tbl$p_adjusted),
    tbl$p_adjusted <= 0.05,
    dplyr::if_else(is.finite(tbl$p_value), tbl$p_value <= 0.05, NA)
  )
  tbl$ContrastBasis <- "group difference in mean observed-minus-expected residuals"
  tbl$SEBasis <- "Welch contrast of residual cell means"
  tbl$StatisticLabel <- "Welch screening t"
  tbl$ProbabilityMetric <- "Welch t tail area"
  tbl$DFBasis <- "Welch-Satterthwaite approximation"
  tbl$EffectMetric <- "mean_obs_minus_exp_difference"
  tbl$ClassificationSystem <- "screening"
  tbl$Classification <- dplyr::case_when(
    !is.finite(tbl$Contrast) ~ NA_character_,
    sig %in% TRUE ~ "Screen positive",
    sig %in% FALSE ~ "Screen negative",
    TRUE ~ NA_character_
  )
  tbl$ETS <- NA_character_
  existing_primary <- tbl$PrimaryReportingEligible %||% rep(NA, nrow(tbl))
  existing_use <- tbl$ReportingUse %||% rep(NA_character_, nrow(tbl))
  tbl$PrimaryReportingEligible <- dplyr::coalesce(as.logical(existing_primary), FALSE)
  tbl$ReportingUse <- dplyr::coalesce(as.character(existing_use), "screening_only")
  tbl
}

build_dff_summary <- function(tbl, method) {
  if (identical(method, "refit")) {
    return(tibble(
      Classification = c(
        "A (Negligible)",
        "B (Moderate)",
        "C (Large)",
        "Linked contrast (screening only)",
        "Unclassified (insufficient linking)"
      ),
      Count = c(
        sum(tbl$ETS == "A", na.rm = TRUE),
        sum(tbl$ETS == "B", na.rm = TRUE),
        sum(tbl$ETS == "C", na.rm = TRUE),
        sum(tbl$Classification == "Linked contrast (screening only)", na.rm = TRUE),
        sum(tbl$Classification == "Unclassified (insufficient linking)", na.rm = TRUE)
      )
    ))
  }

  tibble(
    Classification = c("Screen positive", "Screen negative", "Unclassified"),
    Count = c(
      sum(tbl$Classification == "Screen positive", na.rm = TRUE),
      sum(tbl$Classification == "Screen negative", na.rm = TRUE),
      sum(is.na(tbl$Classification), na.rm = TRUE)
    )
  )
}

resolve_dff_refit_controls <- function(fit) {
  control <- fit$config$estimation_control %||% list()
  list(
    model = fit$config$model %||% fit$summary$Model[1] %||% "RSM",
    method = fit$config$method %||% fit$summary$Method[1] %||% "JMLE",
    step_facet = fit$config$step_facet %||% NULL,
    weight = fit$config$weight_col %||% NULL,
    noncenter_facet = fit$config$noncenter_facet %||% "Person",
    dummy_facets = fit$config$dummy_facets %||% NULL,
    positive_facets = fit$config$positive_facets %||% NULL,
    quad_points = as.integer(control$quad_points %||% 15L),
    maxit = max(25L, min(as.integer(control$maxit %||% 50L), 100L)),
    reltol = as.numeric(control$reltol %||% 1e-6)
  )
}

build_dff_linking_setup <- function(fit, facet, facet_names) {
  linking_facets <- setdiff(as.character(facet_names), as.character(facet))
  anchor_tbl <- if (length(linking_facets) > 0) {
    make_anchor_table(fit, facets = linking_facets, include_person = FALSE)
  } else {
    tibble::tibble(Facet = character(0), Level = character(0), Anchor = numeric(0))
  }
  min_common_anchors <- fit$config$anchor_audit$thresholds$min_common_anchors %||% 5L
  list(
    linking_facets = linking_facets,
    anchor_tbl = tibble::as_tibble(anchor_tbl),
    min_common_anchors = as.integer(min_common_anchors)
  )
}

summarize_dff_group_linkage <- function(sub_fit, linking_setup) {
  if (length(linking_setup$linking_facets) == 0 || nrow(linking_setup$anchor_tbl) == 0) {
    return(list(
      status = "unlinked",
      ets_eligible = FALSE,
      anchored_levels = 0L,
      detail = "No non-target linking facets were available for anchored subgroup refits."
    ))
  }

  audit <- sub_fit$config$anchor_audit %||% list()
  facet_summary <- as.data.frame(audit$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(facet_summary) == 0 || !"Facet" %in% names(facet_summary)) {
    return(list(
      status = "weak_link",
      ets_eligible = FALSE,
      anchored_levels = NA_integer_,
      detail = "Linking anchors were requested, but subgroup anchor coverage could not be audited."
    ))
  }

  link_tbl <- facet_summary[facet_summary$Facet %in% linking_setup$linking_facets, , drop = FALSE]
  if (nrow(link_tbl) == 0) {
    return(list(
      status = "unlinked",
      ets_eligible = FALSE,
      anchored_levels = 0L,
      detail = "None of the requested linking facets were present in the subgroup data."
    ))
  }

  anchored_levels <- if ("AnchoredLevels" %in% names(link_tbl)) {
    suppressWarnings(as.numeric(link_tbl$AnchoredLevels))
  } else {
    rep(NA_real_, nrow(link_tbl))
  }
  anchored_levels[!is.finite(anchored_levels)] <- 0

  total_anchored <- sum(anchored_levels, na.rm = TRUE)
  strong_link <- is.finite(total_anchored) && total_anchored >= linking_setup$min_common_anchors
  weak_link <- any(anchored_levels > 0)
  status <- if (strong_link) {
    "linked"
  } else if (weak_link) {
    "weak_link"
  } else {
    "unlinked"
  }

  detail <- paste0(
    "Linking facets: ",
    paste0(link_tbl$Facet, "=", anchored_levels, collapse = ", "),
    " anchored level(s); threshold=",
    linking_setup$min_common_anchors,
    "."
  )

  list(
    status = status,
    ets_eligible = identical(status, "linked"),
    anchored_levels = as.integer(total_anchored),
    detail = detail
  )
}

extract_dff_group_estimates <- function(sub_fit, sub_diag, facet, fallback_levels, n_obs,
                                        linking_setup, linkage, diagnostics_error = NULL) {
  precision_meta <- resolve_dff_subgroup_precision(
    sub_fit,
    sub_diag = sub_diag,
    diagnostics_error = diagnostics_error
  )
  link_detail <- linkage$detail
  if (!is.null(diagnostics_error) && nzchar(diagnostics_error)) {
    diag_note <- paste0("Subgroup diagnostics failed: ", diagnostics_error)
    link_detail <- if (!is.null(link_detail) && nzchar(link_detail)) {
      paste(link_detail, diag_note)
    } else {
      diag_note
    }
  }
  if (!is.null(sub_diag) && !is.null(sub_diag$measures)) {
    sub_measures <- tibble::as_tibble(sub_diag$measures)
    sub_est <- sub_measures |>
      filter(.data$Facet == facet) |>
      select("Level", "Estimate", "SE")
  } else if (!is.null(sub_fit) && !is.null(sub_fit$facets$others)) {
    sub_others <- tibble::as_tibble(sub_fit$facets$others)
    sub_est <- sub_others |>
      filter(.data$Facet == facet) |>
      select("Level", "Estimate") |>
      mutate(SE = NA_real_)
  } else {
    sub_est <- tibble(Level = fallback_levels, Estimate = NA_real_, SE = NA_real_)
  }

  if (nrow(sub_est) == 0) {
    sub_est <- tibble(Level = fallback_levels, Estimate = NA_real_, SE = NA_real_)
  }

  sub_est |>
    mutate(
      N = as.integer(n_obs),
      LinkingFacets = if (length(linking_setup$linking_facets) > 0) {
        paste(linking_setup$linking_facets, collapse = ", ")
      } else {
        NA_character_
      },
      LinkingThreshold = as.integer(linking_setup$min_common_anchors),
      LinkingStatus = linkage$status,
      LinkingAnchoredLevels = linkage$anchored_levels,
      LinkingDetail = link_detail,
      LinkComparable = isTRUE(linkage$ets_eligible),
      Converged = isTRUE(precision_meta$converged),
      PrecisionTier = precision_meta$precision_tier,
      SupportsFormalInference = isTRUE(precision_meta$supports_formal),
      SubgroupMethod = precision_meta$method,
      DiagnosticsStatus = if (!is.null(diagnostics_error) && nzchar(diagnostics_error)) "failed" else "available",
      DiagnosticsDetail = if (!is.null(diagnostics_error) && nzchar(diagnostics_error)) diagnostics_error else NA_character_,
      ETS_Eligible = isTRUE(linkage$ets_eligible) && isTRUE(precision_meta$supports_formal)
    )
}

#' Differential facet functioning analysis
#'
#' Tests whether the difficulty of facet levels differs across a grouping
#' variable (e.g., whether rater severity differs for male vs. female
#' examinees, or whether item difficulty differs across rater subgroups).
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet Character scalar naming the facet whose elements are tested
#'   for differential functioning (for example, `"Criterion"` or `"Rater"`).
#' @param group Character scalar naming the column in the data that
#'   defines the grouping variable (e.g., `"Gender"`, `"Site"`).
#' @param data Optional data frame containing at least the group column
#'   and the same person/facet/score columns used to fit the model. If
#'   `NULL` (default), the data stored in `fit$prep$data` is used.
#' @param focal Optional character vector of group levels to treat as focal.
#'   If `NULL` (default), all pairwise group comparisons are performed.
#' @param method Analysis method: `"residual"` (default) uses the fitted
#'   model's residuals without re-estimation; `"refit"` re-estimates the
#'   model within each group subset. The residual method is faster and
#'   avoids convergence issues with small subsets.
#' @param min_obs Minimum number of observations per cell (facet-level x
#'   group). Cells below this threshold are flagged as sparse and their
#'   statistics set to `NA`. Default `10`.
#' @param p_adjust Method for multiple-comparison adjustment, passed to
#'   [stats::p.adjust()]. Default is `"holm"`.
#'
#' @details
#' **Differential facet functioning (DFF)** occurs when the
#' difficulty or severity of a facet element differs across subgroups
#' of the population, after controlling for overall ability.  In an
#' MFRM context this generalises classical DIF (which applies to
#' items) to any facet: raters, criteria, tasks, etc.
#'
#' Differential functioning is a threat to measurement fairness: if Criterion 1 is harder
#' for Group A than Group B at the same ability level, the measurement
#' scale is no longer group-invariant.
#'
#' Two methods are available:
#'
#' **Residual method** (`method = "residual"`): Uses the existing fitted
#' model's observation-level residuals.  For each facet-level \eqn{\times}
#' group cell, the observed and expected score sums are aggregated and
#' a standardized residual is computed as:
#' \deqn{z = \frac{\sum (X_{obs} - E_{exp})}{\sqrt{\sum \mathrm{Var}}}}
#' Pairwise contrasts between groups compare the mean observed-minus-expected
#' difference for each facet level, with uncertainty summarized by a
#' Welch/Satterthwaite approximation. This method is fast, stable with small
#' subsets, and does not require re-estimation. Because the resulting contrast
#' is not a logit-scale parameter difference, the residual method is treated as
#' a screening procedure rather than an ETS-style classifier.
#'
#' **Refit method** (`method = "refit"`): Subsets the data by group, refits
#' the MFRM model within each subset, anchors all non-target facets back to
#' the baseline calibration when possible, and compares the resulting
#' facet-level estimates using a Welch t-statistic:
#' \deqn{t = \frac{\hat{\delta}_1 - \hat{\delta}_2}
#'                {\sqrt{SE_1^2 + SE_2^2}}}
#' This provides group-specific parameter estimates on a common scale when
#' linking anchors are available, but is slower and may encounter convergence
#' issues with small subsets.  ETS categories are reported only for contrasts
#' whose subgroup calibrations retained enough linking anchors to support a
#' common-scale interpretation and whose subgroup precision remained on the
#' package's model-based MML path.
#'
#' When `facet` refers to an item-like facet (for example `Criterion`), this
#' recovers the familiar DIF case. When `facet` refers to raters or
#' prompts/tasks, the same machinery supports DRF/DPF-style analyses.
#'
#' For the refit method only, effect size is classified following the ETS
#' (Educational Testing Service) DIF guidelines when subgroup calibrations are
#' both linked and eligible for model-based inference:
#' - **A (Negligible)**: \eqn{|\Delta| <} 0.43 logits
#' - **B (Moderate)**: 0.43 \eqn{\le |\Delta| <} 0.64 logits
#' - **C (Large)**: \eqn{|\Delta| \ge} 0.64 logits
#'
#' Multiple comparisons are adjusted using Holm's step-down procedure by
#' default, which controls the family-wise error rate without assuming
#' independence.  Alternative methods (e.g., `"BH"` for false discovery
#' rate) can be specified via `p_adjust`.
#'
#' @section Choosing a method:
#' In most first-pass DFF screening, start with `method = "residual"`. It is
#' faster, reuses the fitted model, and is less fragile in smaller subsets.
#' Use `method = "refit"` when you specifically want group-specific parameter
#' estimates and can tolerate extra computation.  Both methods should yield
#' similar conclusions when sample sizes are adequate (\eqn{N \ge 100} per
#' group is a useful guideline for stable differential-functioning detection).
#'
#' @section Interpreting output:
#' - `$dif_table`: one row per facet-level x group-pair with contrast,
#'   SE, t-statistic, p-value, adjusted p-value, effect metric, and
#'   method-appropriate classification. Includes `Method`, `N_Group1`,
#'   `N_Group2`, `EffectMetric`, `ClassificationSystem`, `ContrastBasis`,
#'   `SEBasis`, `StatisticLabel`, `ProbabilityMetric`, `DFBasis`,
#'   `ReportingUse`, `PrimaryReportingEligible`, and `sparse` columns.
#' - `$cell_table`: (residual method only) per-cell detail with N,
#'   ObsScore, ExpScore, ObsExpAvg, StdResidual.
#' - `$summary`: counts by screening result (`method = "residual"`) or ETS
#'   category plus linked-screening and insufficient-linking rows
#'   (`method = "refit"`).
#' - `$group_fits`: (refit method only) list of per-group facet estimates and
#'   subgroup linking diagnostics.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Run `analyze_dff(fit, diagnostics, facet = "Criterion", group = "Gender", data = my_data)`.
#' 3. Inspect `$dif_table` for flagged levels and `$summary` for counts.
#' 4. Use [dif_interaction_table()] when you need cell-level diagnostics.
#' 5. Use [plot_dif_heatmap()] or [dif_report()] for communication.
#'
#' @return
#' An object of class `mfrm_dff` (with compatibility class `mfrm_dif`) with:
#' - `dif_table`: data.frame of differential-functioning contrasts.
#' - `cell_table`: (residual method) per-cell detail table.
#' - `summary`: counts by screening or ETS classification.
#' - `group_fits`: (refit method) per-group facet estimates.
#' - `config`: list with facet, group, method, min_obs, p_adjust settings.
#'
#' @seealso [fit_mfrm()], [estimate_bias()], [compare_mfrm()],
#'   [dif_interaction_table()], [plot_dif_heatmap()], [dif_report()],
#'   [subset_connectivity_report()], [mfrmr_linking_and_dff]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' dff <- analyze_dff(fit, diag, facet = "Rater", group = "Group", data = toy)
#' dff$summary
#' head(dff$dif_table[, c("Level", "Group1", "Group2", "Contrast", "Classification")])
#' sc <- subset_connectivity_report(fit, diagnostics = diag)
#' plot(sc, type = "design_matrix", draw = FALSE)
#' if ("ScaleLinkStatus" %in% names(dff$dif_table)) {
#'   unique(dff$dif_table$ScaleLinkStatus)
#' }
#' @rdname analyze_dff
#' @export
analyze_dff <- function(fit,
                        diagnostics,
                        facet,
                        group,
                        data = NULL,
                        focal = NULL,
                        method = c("residual", "refit"),
                        min_obs = 10,
                        p_adjust = "holm") {
  method <- match.arg(method)
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an `mfrm_fit` object.", call. = FALSE)
  }
  if (!is.character(facet) || length(facet) != 1) {
    stop("`facet` must be a single character string.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1) {
    stop("`group` must be a single character string naming a column in the ",
         "original data.", call. = FALSE)
  }
  if (!is.numeric(min_obs) || length(min_obs) != 1 || min_obs < 1) {
    stop("`min_obs` must be a positive integer.", call. = FALSE)
  }

  # Recover data
  orig_data <- if (!is.null(data)) data else fit$prep$data
  if (is.null(orig_data) || !is.data.frame(orig_data)) {
    stop("No data available. Pass the original data via the `data` argument.",
         call. = FALSE)
  }
  if (!group %in% names(orig_data)) {
    stop("`group` column '", group, "' not found in the data. ",
         "Available columns: ", paste(names(orig_data), collapse = ", "),
         call. = FALSE)
  }

  facet_names <- fit$config$facet_cols
  if (is.null(facet_names)) facet_names <- fit$prep$facet_names
  if (!facet %in% facet_names) {
    stop("`facet` '", facet, "' is not one of the model facets: ",
         paste(facet_names, collapse = ", "), ".", call. = FALSE)
  }

  # Group levels
  group_vals <- as.character(orig_data[[group]])
  group_levels <- sort(unique(group_vals))
  if (length(group_levels) < 2) {
    stop("Grouping variable '", group, "' must have at least 2 levels. ",
         "Found: ", length(group_levels), ".", call. = FALSE)
  }

  if (method == "residual") {
    out <- .analyze_dif_residual(
      fit = fit, facet = facet, group = group, data = data,
      orig_data = orig_data, facet_names = facet_names,
      group_levels = group_levels, focal = focal,
      min_obs = min_obs, p_adjust = p_adjust
    )
  } else {
    out <- .analyze_dif_refit(
      fit = fit, diagnostics = diagnostics,
      facet = facet, group = group,
      orig_data = orig_data, facet_names = facet_names,
      group_vals = as.character(orig_data[[group]]),
      group_levels = group_levels, focal = focal,
      min_obs = min_obs, p_adjust = p_adjust
    )
  }
  out
}

#' Backward-compatible alias for [analyze_dff()]
#'
#' `analyze_dif()` is retained for compatibility with earlier package versions.
#' In many-facet workflows, prefer [analyze_dff()] as the primary entry point.
#' @param ... Passed directly to [analyze_dff()].
#' @rdname analyze_dff
#' @export
analyze_dif <- function(...) {
  analyze_dff(...)
}

# Internal: residual-based DIF analysis
.analyze_dif_residual <- function(fit, facet, group, data, orig_data,
                                  facet_names, group_levels, focal,
                                  min_obs, p_adjust) {
  # Compute observation table from the fitted model

  obs_tbl <- compute_obs_table(fit)

  # The obs_tbl has columns from prep$data (Person, facet cols, Score,

  # Weight, score_k) plus PersonMeasure, Observed, Expected, Var, Residual,
  # StdResidual, StdSq. We need to merge with orig_data to get the group col.
  # Use row-order matching: prep$data is a cleaned subset of the original
  # data, so we need the group column from the user-supplied data.
  person_col <- fit$config$person_col %||% "Person"
  score_col <- fit$config$score_col %||% "Score"

  # Build a merge key from the obs_tbl (Person + facets)
  merge_cols <- c("Person", facet_names)

  # Add group column from orig_data by joining on person + facets + score
  orig_for_merge <- orig_data
  orig_for_merge$.group_var <- as.character(orig_data[[group]])
  # Rename to match internal names
  if (person_col != "Person") {
    orig_for_merge$Person <- as.character(orig_for_merge[[person_col]])
  } else {
    orig_for_merge$Person <- as.character(orig_for_merge$Person)
  }
  for (fn in facet_names) {
    orig_for_merge[[fn]] <- as.character(orig_for_merge[[fn]])
  }

  # Ensure obs_tbl facet columns are character for merging
  obs_tbl_chr <- obs_tbl
  obs_tbl_chr$Person <- as.character(obs_tbl_chr$Person)
  for (fn in facet_names) {
    obs_tbl_chr[[fn]] <- as.character(obs_tbl_chr[[fn]])
  }

  # Use left_join on all merge keys; include Score to handle duplicate combos
  merge_key <- c("Person", facet_names)
  join_cols <- merge_key
  names(join_cols) <- merge_key

  # Add row index for stable matching
  obs_tbl_chr$.obs_row <- seq_len(nrow(obs_tbl_chr))
  orig_for_merge$.orig_row <- seq_len(nrow(orig_for_merge))

  merged <- left_join(
    obs_tbl_chr |> select(all_of(c(".obs_row", merge_key))),
    orig_for_merge |> select(all_of(c(merge_key, ".group_var"))) |> distinct(),
    by = merge_key
  )

  # Handle duplicates from many-to-many: keep first group per obs row
  merged <- merged |>
    group_by(.data$.obs_row) |>
    slice(1L) |>
    ungroup() |>
    arrange(.data$.obs_row)

  obs_tbl_chr$.group_var <- merged$.group_var

  # Filter to rows with valid group
  obs_work <- obs_tbl_chr |>
    filter(!is.na(.data$.group_var))

  # Aggregate by facet level x group
  cell_table <- obs_work |>
    group_by(.data[[facet]], .data$.group_var) |>
    summarise(
      N = n(),
      ObsScore = sum(.data$Observed, na.rm = TRUE),
      ExpScore = sum(.data$Expected, na.rm = TRUE),
      ObsExpAvg = mean(.data$Observed - .data$Expected, na.rm = TRUE),
      Var_sum = sum(.data$Var, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      sparse = .data$N < min_obs,
      StdResidual = ifelse(
        .data$sparse | .data$Var_sum <= 0,
        NA_real_,
        (.data$ObsScore - .data$ExpScore) / sqrt(.data$Var_sum)
      ),
      t = .data$StdResidual,
      df = ifelse(.data$sparse, NA_real_, .data$N - 1),
      p_value = ifelse(
        is.finite(.data$t) & is.finite(.data$df) & .data$df > 0,
        2 * stats::pt(abs(.data$t), df = .data$df, lower.tail = FALSE),
        NA_real_
      )
    )
  names(cell_table)[names(cell_table) == facet] <- "Level"
  names(cell_table)[names(cell_table) == ".group_var"] <- "GroupValue"

  # Build pairwise contrasts
  method_label <- as.character(fit$summary$Method[1] %||% fit$config$method %||% NA_character_)
  method_label <- ifelse(identical(method_label, "JMLE"), "JML", method_label)
  if (!is.null(focal)) {
    pairs <- expand_grid(
      Group1 = setdiff(group_levels, focal),
      Group2 = focal
    )
  } else {
    pairs <- as_tibble(as.data.frame(
      t(combn(group_levels, 2)),
      stringsAsFactors = FALSE
    ))
    names(pairs) <- c("Group1", "Group2")
  }

  facet_levels <- unique(cell_table$Level)
  dif_rows <- list()
  for (i in seq_len(nrow(pairs))) {
    g1 <- pairs$Group1[i]
    g2 <- pairs$Group2[i]
    for (lev in facet_levels) {
      c1 <- cell_table |> filter(.data$Level == lev, .data$GroupValue == g1)
      c2 <- cell_table |> filter(.data$Level == lev, .data$GroupValue == g2)
      n1 <- if (nrow(c1) > 0) c1$N[1] else 0L
      n2 <- if (nrow(c2) > 0) c2$N[1] else 0L
      is_sparse <- (n1 < min_obs) || (n2 < min_obs)
      avg1 <- if (nrow(c1) > 0 && !is_sparse) c1$ObsExpAvg[1] else NA_real_
      avg2 <- if (nrow(c2) > 0 && !is_sparse) c2$ObsExpAvg[1] else NA_real_
      contrast <- if (is.finite(avg1) && is.finite(avg2)) avg1 - avg2 else NA_real_
      # Welch-style SE for a contrast of cell-level mean residuals
      var1 <- if (nrow(c1) > 0 && !is_sparse) c1$Var_sum[1] else NA_real_
      var2 <- if (nrow(c2) > 0 && !is_sparse) c2$Var_sum[1] else NA_real_
      comp1 <- if (is.finite(var1) && n1 > 0) var1 / n1^2 else NA_real_
      comp2 <- if (is.finite(var2) && n2 > 0) var2 / n2^2 else NA_real_
      se_diff <- if (is.finite(comp1) && is.finite(comp2) && (comp1 + comp2) > 0) {
        sqrt(comp1 + comp2)
      } else {
        NA_real_
      }
      t_val <- if (is.finite(contrast) && is.finite(se_diff) && se_diff > 0) {
        contrast / se_diff
      } else {
        NA_real_
      }
      df_val <- if (!is_sparse) {
        welch_satterthwaite_df(c(comp1, comp2), c(n1 - 1, n2 - 1))
      } else {
        NA_real_
      }
      p_val <- if (is.finite(t_val) && is.finite(df_val) && df_val > 0) {
        2 * stats::pt(abs(t_val), df = df_val, lower.tail = FALSE)
      } else {
        NA_real_
      }
      abs_diff <- abs(contrast)
      dif_rows[[length(dif_rows) + 1]] <- tibble(
        Level = lev,
        Group1 = g1,
        Group2 = g2,
        Contrast = contrast,
        SE = se_diff,
        t = t_val,
        df = df_val,
        p_value = p_val,
        AbsDiff = abs_diff,
        Method = "residual",
        N_Group1 = as.integer(n1),
        N_Group2 = as.integer(n2),
        sparse = is_sparse,
        ContrastComparable = FALSE,
        FormalInferenceEligible = FALSE,
        PrimaryReportingEligible = FALSE,
        InferenceTier = "screening",
        ComparisonMethod = method_label,
        ScaleLinkStatus = "not_applicable",
        ReportingUse = "screening_only"
      )
    }
  }
  dif_table <- bind_rows(dif_rows)

  # Adjust p-values
  if (nrow(dif_table) > 0 && any(is.finite(dif_table$p_value))) {
    dif_table$p_adjusted <- stats::p.adjust(dif_table$p_value, method = p_adjust)
  } else {
    dif_table$p_adjusted <- NA_real_
  }
  dif_table <- annotate_dff_table(dif_table, method = "residual")

  # Summary counts
  dif_summary <- build_dff_summary(dif_table, method = "residual")
  functioning_label <- functioning_label_for_facet(facet)

  out <- list(
    dif_table = dif_table,
    cell_table = cell_table,
    summary = dif_summary,
    group_fits = NULL,
    config = list(facet = facet, group = group, method = "residual",
                  min_obs = min_obs, p_adjust = p_adjust,
                  focal = focal, group_levels = group_levels,
                  functioning_label = functioning_label)
  )
  class(out) <- c("mfrm_dff", "mfrm_dif", class(out))
  out
}

# Internal: refit-based DIF analysis (original approach)
.analyze_dif_refit <- function(fit, diagnostics, facet, group, orig_data,
                               facet_names, group_vals, group_levels, focal,
                               min_obs, p_adjust) {
  # Get full-sample facet estimates
  measures <- tibble::as_tibble(diagnostics$measures)
  facet_estimates <- measures |>
    filter(.data$Facet == facet) |>
    select("Level", "Estimate", "SE")

  person_col <- fit$config$person_col %||% "Person"
  score_col <- fit$config$score_col %||% "Score"
  refit_controls <- resolve_dff_refit_controls(fit)
  linking_setup <- build_dff_linking_setup(fit, facet = facet, facet_names = facet_names)
  baseline_precision_meta <- resolve_dff_subgroup_precision(fit, diagnostics)

  group_fits <- list()
  for (g in group_levels) {
    idx <- group_vals == g
    sub_data <- orig_data[idx, , drop = FALSE]
    if (nrow(sub_data) < 5) {
      group_fits[[g]] <- facet_estimates |>
        mutate(
          N = 0L,
          LinkingFacets = if (length(linking_setup$linking_facets) > 0) {
            paste(linking_setup$linking_facets, collapse = ", ")
          } else {
            NA_character_
          },
          LinkingThreshold = as.integer(linking_setup$min_common_anchors),
          LinkingStatus = "insufficient_data",
          LinkingAnchoredLevels = 0L,
          LinkingDetail = "Subgroup had fewer than 5 observations; anchored refit was skipped.",
          LinkComparable = FALSE,
          Converged = FALSE,
          PrecisionTier = NA_character_,
          SupportsFormalInference = FALSE,
          SubgroupMethod = refit_controls$method,
          ETS_Eligible = FALSE
        )
      next
    }
    fit_args <- list(
      data = sub_data,
      person = person_col,
      facets = facet_names,
      score = score_col,
      weight = refit_controls$weight,
      method = refit_controls$method,
      model = refit_controls$model,
      step_facet = refit_controls$step_facet,
      anchors = if (nrow(linking_setup$anchor_tbl) > 0) linking_setup$anchor_tbl else NULL,
      noncenter_facet = refit_controls$noncenter_facet,
      dummy_facets = refit_controls$dummy_facets,
      positive_facets = refit_controls$positive_facets,
      anchor_policy = "silent",
      quad_points = refit_controls$quad_points,
      maxit = refit_controls$maxit,
      reltol = refit_controls$reltol
    )
    sub_fit <- tryCatch(
      suppressWarnings(do.call(fit_mfrm, fit_args)),
      error = function(e) structure(list(message = conditionMessage(e)), class = "mfrm_dff_fit_error")
    )
    if (inherits(sub_fit, "mfrm_dff_fit_error")) {
      group_fits[[g]] <- facet_estimates |>
        mutate(
          N = sum(idx),
          LinkingFacets = if (length(linking_setup$linking_facets) > 0) {
            paste(linking_setup$linking_facets, collapse = ", ")
          } else {
            NA_character_
          },
          LinkingThreshold = as.integer(linking_setup$min_common_anchors),
          LinkingStatus = if (nrow(linking_setup$anchor_tbl) > 0) "failed" else "unlinked",
          LinkingAnchoredLevels = NA_integer_,
          LinkingDetail = sub_fit$message %||% "Anchored subgroup refit failed.",
          LinkComparable = FALSE,
          Converged = FALSE,
          PrecisionTier = NA_character_,
          SupportsFormalInference = FALSE,
          SubgroupMethod = refit_controls$method,
          ETS_Eligible = FALSE
        )
    } else {
      linkage <- summarize_dff_group_linkage(sub_fit, linking_setup = linking_setup)
      sub_diag_error <- NULL
      sub_diag <- tryCatch(
        suppressWarnings(diagnose_mfrm(sub_fit, residual_pca = "none")),
        error = function(e) {
          sub_diag_error <<- conditionMessage(e)
          NULL
        }
      )
      group_fits[[g]] <- extract_dff_group_estimates(
        sub_fit = sub_fit,
        sub_diag = sub_diag,
        facet = facet,
        fallback_levels = facet_estimates$Level,
        n_obs = sum(idx),
        linking_setup = linking_setup,
        linkage = linkage,
        diagnostics_error = sub_diag_error
      )
    }
  }

  # Build DIF contrasts
  if (!is.null(focal)) {
    pairs <- expand_grid(
      Group1 = setdiff(group_levels, focal),
      Group2 = focal
    )
  } else {
    pairs <- as_tibble(as.data.frame(
      t(combn(group_levels, 2)),
      stringsAsFactors = FALSE
    ))
    names(pairs) <- c("Group1", "Group2")
  }

  dif_rows <- list()
  for (i in seq_len(nrow(pairs))) {
    g1 <- pairs$Group1[i]
    g2 <- pairs$Group2[i]
    est1 <- group_fits[[g1]]
    est2 <- group_fits[[g2]]
    merged <- merge(est1, est2, by = "Level", suffixes = c("_1", "_2"))
    for (j in seq_len(nrow(merged))) {
      e1 <- merged$Estimate_1[j]
      e2 <- merged$Estimate_2[j]
      se1 <- merged$SE_1[j]
      se2 <- merged$SE_2[j]
      n1 <- merged$N_1[j]
      n2 <- merged$N_2[j]
      link_comparable <- isTRUE(merged$LinkComparable_1[j]) && isTRUE(merged$LinkComparable_2[j])
      subgroup_formal <- isTRUE(merged$SupportsFormalInference_1[j]) &&
        isTRUE(merged$SupportsFormalInference_2[j])
      subgroup_converged <- isTRUE(merged$Converged_1[j]) &&
        isTRUE(merged$Converged_2[j])
      comparison_method <- dplyr::coalesce(
        merged$SubgroupMethod_1[j],
        merged$SubgroupMethod_2[j],
        baseline_precision_meta$method
      )
      inference_tier <- dplyr::case_when(
        all(c(merged$PrecisionTier_1[j], merged$PrecisionTier_2[j]) == "model_based") ~ "model_based",
        any(c(merged$PrecisionTier_1[j], merged$PrecisionTier_2[j]) == "exploratory") ~ "exploratory",
        any(c(merged$PrecisionTier_1[j], merged$PrecisionTier_2[j]) == "hybrid") ~ "hybrid",
        TRUE ~ NA_character_
      )
      scale_link_status <- dplyr::case_when(
        any(c(merged$LinkingStatus_1[j], merged$LinkingStatus_2[j]) == "failed") ~ "failed",
        any(c(merged$LinkingStatus_1[j], merged$LinkingStatus_2[j]) == "insufficient_data") ~ "insufficient_data",
        link_comparable ~ "linked",
        any(c(merged$LinkingStatus_1[j], merged$LinkingStatus_2[j]) == "weak_link") ~ "weak_link",
        any(c(merged$LinkingStatus_1[j], merged$LinkingStatus_2[j]) == "unlinked") ~ "unlinked",
        TRUE ~ "unlinked"
      )
      is_sparse <- (n1 < min_obs) || (n2 < min_obs)
      formal_eligible <- link_comparable &&
        subgroup_formal &&
        subgroup_converged &&
        isTRUE(baseline_precision_meta$supports_formal) &&
        isTRUE(baseline_precision_meta$converged) &&
        identical(comparison_method, "MML") &&
        !is_sparse
      comparable <- formal_eligible
      reporting_use <- dplyr::case_when(
        formal_eligible ~ "primary_reporting",
        link_comparable && identical(comparison_method, "MML") ~ "review_before_reporting",
        link_comparable && identical(inference_tier, "hybrid") ~ "review_before_reporting",
        TRUE ~ "screening_only"
      )
      contrast <- e1 - e2
      se_diff <- if (comparable) sqrt(se1^2 + se2^2) else NA_real_
      t_val <- if (is.finite(se_diff) && se_diff > 0) contrast / se_diff else NA_real_
      df_welch <- if (comparable && is.finite(se1) && is.finite(se2) && se1 > 0 && se2 > 0) {
        welch_satterthwaite_df(c(se1^2, se2^2), c(n1 - 1, n2 - 1))
      } else {
        NA_real_
      }
      p_val <- if (is.finite(t_val) && is.finite(df_welch) && df_welch > 0) {
        2 * stats::pt(abs(t_val), df = df_welch, lower.tail = FALSE)
      } else {
        NA_real_
      }
      abs_diff <- abs(contrast)
      dif_rows[[length(dif_rows) + 1]] <- tibble(
        Level = merged$Level[j],
        Group1 = g1,
        Group2 = g2,
        Estimate1 = e1,
        Estimate2 = e2,
        Contrast = contrast,
        SE = se_diff,
        t = t_val,
        df = df_welch,
        p_value = p_val,
        AbsDiff = abs_diff,
        Method = "refit",
        N_Group1 = as.integer(n1),
        N_Group2 = as.integer(n2),
        sparse = is_sparse,
        ContrastComparable = link_comparable,
        FormalInferenceEligible = formal_eligible,
        PrimaryReportingEligible = formal_eligible,
        InferenceTier = inference_tier,
        ComparisonMethod = comparison_method,
        ReportingUse = reporting_use,
        ETS_Eligible = comparable,
        ScaleLinkStatus = scale_link_status,
        BaselineMethod = baseline_precision_meta$method,
        BaselineConverged = isTRUE(baseline_precision_meta$converged),
        BaselinePrecisionTier = baseline_precision_meta$precision_tier,
        BaselineSupportsFormalInference = isTRUE(baseline_precision_meta$supports_formal),
        SubgroupConverged1 = isTRUE(merged$Converged_1[j]),
        SubgroupConverged2 = isTRUE(merged$Converged_2[j]),
        LinkingFacets = merged$LinkingFacets_1[j] %||% merged$LinkingFacets_2[j],
        LinkingThreshold = merged$LinkingThreshold_1[j] %||% merged$LinkingThreshold_2[j],
        LinkingStatus1 = merged$LinkingStatus_1[j],
        LinkingStatus2 = merged$LinkingStatus_2[j],
        LinkingAnchoredLevels1 = merged$LinkingAnchoredLevels_1[j],
        LinkingAnchoredLevels2 = merged$LinkingAnchoredLevels_2[j],
        LinkingDetail1 = merged$LinkingDetail_1[j],
        LinkingDetail2 = merged$LinkingDetail_2[j]
      )
    }
  }
  dif_table <- bind_rows(dif_rows)

  # Adjust p-values
  if (nrow(dif_table) > 0 && any(is.finite(dif_table$p_value))) {
    dif_table$p_adjusted <- stats::p.adjust(dif_table$p_value, method = p_adjust)
  } else {
    dif_table$p_adjusted <- NA_real_
  }
  dif_table <- annotate_dff_table(dif_table, method = "refit")

  # Summary counts
  dif_summary <- build_dff_summary(dif_table, method = "refit")
  functioning_label <- functioning_label_for_facet(facet)

  out <- list(
    dif_table = dif_table,
    cell_table = NULL,
    summary = dif_summary,
    group_fits = group_fits,
    config = list(facet = facet, group = group, method = "refit",
                  min_obs = min_obs, p_adjust = p_adjust,
                  focal = focal, group_levels = group_levels,
                  linking_facets = linking_setup$linking_facets,
                  linking_threshold = linking_setup$min_common_anchors,
                  functioning_label = functioning_label)
  )
  class(out) <- c("mfrm_dff", "mfrm_dif", class(out))
  out
}

#' @export
summary.mfrm_dif <- function(object, ...) {
  out <- list(
    dif_table = object$dif_table,
    cell_table = object$cell_table,
    summary = object$summary,
    config = object$config
  )
  class(out) <- "summary.mfrm_dif"
  out
}

#' @export
summary.mfrm_dff <- function(object, ...) {
  summary.mfrm_dif(object, ...)
}

#' @export
print.summary.mfrm_dif <- function(x, ...) {
  label <- x$config$functioning_label %||% "DFF"
  cat("--- ", label, " Analysis ---\n", sep = "")
  cat("Method:", x$config$method %||% "refit", "\n")
  cat("Facet:", x$config$facet, " | Group:", x$config$group, "\n")
  cat("Groups:", paste(x$config$group_levels, collapse = ", "), "\n")
  if (identical(x$config$method, "refit")) {
    link_txt <- if (!is.null(x$config$linking_facets) && length(x$config$linking_facets) > 0) {
      paste(x$config$linking_facets, collapse = ", ")
    } else {
      "none"
    }
    cat("Linking facets:", link_txt)
    if (!is.null(x$config$linking_threshold)) {
      cat(" | Anchor threshold:", x$config$linking_threshold)
    }
    cat("\n")
  }
  if (!is.null(x$config$min_obs)) {
    cat("Min observations per cell:", x$config$min_obs, "\n")
  }
  cat("\n")

  if (nrow(x$dif_table) > 0) {
    show_cols <- intersect(
      c("Level", "Group1", "Group2", "Contrast", "SE", "t",
        "p_adjusted", "Classification", "ETS",
        "ReportingUse", "PrimaryReportingEligible",
        "N_Group1", "N_Group2", "sparse"),
      names(x$dif_table)
    )
    print_tbl <- x$dif_table |> select(all_of(show_cols))
    print(as.data.frame(print_tbl), row.names = FALSE, digits = 3)
  } else {
    cat("No differential-functioning contrasts computed.\n")
  }

  if (identical(x$config$method, "refit")) {
    cat("\nRefit Classification Summary:\n")
  } else {
    cat("\nScreening Summary:\n")
  }
  print(as.data.frame(x$summary), row.names = FALSE)
  invisible(x)
}

#' @export
print.summary.mfrm_dff <- function(x, ...) {
  print.summary.mfrm_dif(x, ...)
}

#' @export
print.mfrm_dif <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' @export
print.mfrm_dff <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

# ============================================================================
# B2. Differential Functioning Interaction Table
# ============================================================================

#' Compute interaction table between a facet and a grouping variable
#'
#' Produces a cell-level interaction table showing Obs-Exp differences,
#' standardized residuals, and screening statistics for each
#' facet-level x group-value cell.
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet Character scalar naming the facet.
#' @param group Character scalar naming the grouping column.
#' @param data Optional data frame with the group column. If `NULL`
#'   (default), the data stored in `fit$prep$data` is used, but it
#'   must contain the `group` column.
#' @param min_obs Minimum observations per cell. Cells with fewer than
#'   this many observations are flagged as sparse and their test
#'   statistics set to `NA`. Default `10`.
#' @param p_adjust P-value adjustment method, passed to
#'   [stats::p.adjust()]. Default `"holm"`.
#' @param abs_t_warn Threshold for flagging cells by absolute t-value.
#'   Default `2`.
#' @param abs_bias_warn Threshold for flagging cells by absolute
#'   Obs-Exp average (in logits). Default `0.5`.
#'
#' @details
#' This function uses the fitted model's observation-level residuals
#' (from the internal `compute_obs_table()` function) rather than
#' re-estimating the model. For each facet-level x group-value cell,
#' it computes:
#' \itemize{
#'   \item N: number of observations in the cell
#'   \item ObsScore: sum of observed scores
#'   \item ExpScore: sum of expected scores
#'   \item ObsExpAvg: mean observed-minus-expected difference
#'   \item Var_sum: sum of model variances
#'   \item StdResidual: (ObsScore - ExpScore) / sqrt(Var_sum)
#'   \item t: approximate t-statistic (equal to StdResidual)
#'   \item df: N - 1
#'   \item p_value: two-tailed p-value from the t-distribution
#' }
#'
#' @section When to use this instead of analyze_dff():
#' Use `dif_interaction_table()` when you want cell-level screening for a
#' single facet-by-group table. Use [analyze_dff()] when you want group-pair
#' contrasts summarized into differential-functioning effect sizes and
#' method-appropriate classifications.
#'
#' @section Further guidance:
#' For plot selection and follow-up diagnostics, see
#' [mfrmr_visual_diagnostics].
#'
#' @section Interpreting output:
#' - `$table`: the full interaction table with one row per cell.
#' - `$summary`: overview counts of flagged and sparse cells.
#' - `$config`: analysis configuration parameters.
#' - Cells with `|t| > abs_t_warn` or `|ObsExpAvg| > abs_bias_warn`
#'   are flagged in the `flag_t` and `flag_bias` columns.
#' - Sparse cells (N < min_obs) have `sparse = TRUE` and NA statistics.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Run `dif_interaction_table(fit, diag, facet = "Rater", group = "Gender", data = df)`.
#' 3. Inspect `$table` for flagged cells.
#' 4. Visualize with [plot_dif_heatmap()].
#'
#' @return Object of class `mfrm_dif_interaction` with:
#' - `table`: tibble with per-cell statistics and flags.
#' - `summary`: tibble summarizing flagged and sparse cell counts.
#' - `config`: list of analysis parameters.
#'
#' @seealso [analyze_dff()], [analyze_dif()], [plot_dif_heatmap()], [dif_report()],
#'   [estimate_bias()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' int <- dif_interaction_table(fit, diag, facet = "Rater",
#'                              group = "Group", data = toy, min_obs = 2)
#' int$summary
#' head(int$table[, c("Level", "GroupValue", "ObsExpAvg", "flag_bias")])
#' @export
dif_interaction_table <- function(fit, diagnostics, facet, group, data = NULL,
                                  min_obs = 10, p_adjust = "holm",
                                  abs_t_warn = 2, abs_bias_warn = 0.5) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an `mfrm_fit` object.", call. = FALSE)
  }
  if (!is.character(facet) || length(facet) != 1) {
    stop("`facet` must be a single character string.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1) {
    stop("`group` must be a single character string.", call. = FALSE)
  }
  if (!is.numeric(min_obs) || length(min_obs) != 1 || min_obs < 1) {
    stop("`min_obs` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(abs_t_warn) || length(abs_t_warn) != 1) {
    stop("`abs_t_warn` must be a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(abs_bias_warn) || length(abs_bias_warn) != 1) {
    stop("`abs_bias_warn` must be a single numeric value.", call. = FALSE)
  }

  # Recover data
  orig_data <- if (!is.null(data)) data else fit$prep$data
  if (is.null(orig_data) || !is.data.frame(orig_data)) {
    stop("No data available. Pass the original data via the `data` argument.",
         call. = FALSE)
  }
  if (!group %in% names(orig_data)) {
    stop("`group` column '", group, "' not found in the data.",
         call. = FALSE)
  }

  facet_names <- fit$config$facet_cols
  if (is.null(facet_names)) facet_names <- fit$prep$facet_names
  if (!facet %in% facet_names) {
    stop("`facet` '", facet, "' is not one of the model facets: ",
         paste(facet_names, collapse = ", "), ".", call. = FALSE)
  }

  group_levels <- sort(unique(as.character(orig_data[[group]])))
  if (length(group_levels) < 2) {
    stop("Grouping variable '", group, "' must have at least 2 levels.",
         call. = FALSE)
  }

  # Compute observation table
  obs_tbl <- compute_obs_table(fit)

  person_col <- fit$config$person_col %||% "Person"

  # Prepare merge keys
  merge_cols <- c("Person", facet_names)
  obs_chr <- obs_tbl
  obs_chr$Person <- as.character(obs_chr$Person)
  for (fn in facet_names) {
    obs_chr[[fn]] <- as.character(obs_chr[[fn]])
  }

  orig_for_merge <- orig_data
  orig_for_merge$.group_var <- as.character(orig_data[[group]])
  if (person_col != "Person") {
    orig_for_merge$Person <- as.character(orig_for_merge[[person_col]])
  } else {
    orig_for_merge$Person <- as.character(orig_for_merge$Person)
  }
  for (fn in facet_names) {
    orig_for_merge[[fn]] <- as.character(orig_for_merge[[fn]])
  }

  obs_chr$.obs_row <- seq_len(nrow(obs_chr))
  merged <- left_join(
    obs_chr |> select(all_of(c(".obs_row", merge_cols))),
    orig_for_merge |> select(all_of(c(merge_cols, ".group_var"))) |> distinct(),
    by = merge_cols
  )
  merged <- merged |>
    group_by(.data$.obs_row) |>
    slice(1L) |>
    ungroup() |>
    arrange(.data$.obs_row)

  obs_chr$.group_var <- merged$.group_var
  obs_work <- obs_chr |> filter(!is.na(.data$.group_var))

  # Aggregate by facet level x group
  int_table <- obs_work |>
    group_by(.data[[facet]], .data$.group_var) |>
    summarise(
      N = n(),
      ObsScore = sum(.data$Observed, na.rm = TRUE),
      ExpScore = sum(.data$Expected, na.rm = TRUE),
      ObsExpAvg = mean(.data$Observed - .data$Expected, na.rm = TRUE),
      Var_sum = sum(.data$Var, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      sparse = .data$N < min_obs,
      StdResidual = ifelse(
        .data$sparse | .data$Var_sum <= 0,
        NA_real_,
        (.data$ObsScore - .data$ExpScore) / sqrt(.data$Var_sum)
      ),
      t = .data$StdResidual,
      df = ifelse(.data$sparse, NA_real_, .data$N - 1),
      p_value = ifelse(
        is.finite(.data$t) & is.finite(.data$df) & .data$df > 0,
        2 * stats::pt(abs(.data$t), df = .data$df, lower.tail = FALSE),
        NA_real_
      )
    )
  names(int_table)[names(int_table) == facet] <- "Level"
  names(int_table)[names(int_table) == ".group_var"] <- "GroupValue"

  # Adjust p-values
  if (nrow(int_table) > 0 && any(is.finite(int_table$p_value))) {
    int_table$p_adjusted <- stats::p.adjust(int_table$p_value, method = p_adjust)
  } else {
    int_table$p_adjusted <- NA_real_
  }

  # Flag cells
  int_table <- int_table |>
    mutate(
      flag_t = ifelse(.data$sparse, NA, abs(.data$t) > abs_t_warn),
      flag_bias = ifelse(.data$sparse, NA, abs(.data$ObsExpAvg) > abs_bias_warn)
    )

  # Summary
  n_total <- nrow(int_table)
  n_sparse <- sum(int_table$sparse, na.rm = TRUE)
  n_flag_t <- sum(int_table$flag_t == TRUE, na.rm = TRUE)
  n_flag_bias <- sum(int_table$flag_bias == TRUE, na.rm = TRUE)
  int_summary <- tibble(
    Metric = c("Total cells", "Sparse cells (N < min_obs)",
               "Flagged by |t|", "Flagged by |Obs-Exp Avg|"),
    Count = c(n_total, n_sparse, n_flag_t, n_flag_bias)
  )

  out <- list(
    table = int_table,
    summary = int_summary,
    config = list(facet = facet, group = group, min_obs = min_obs,
                  p_adjust = p_adjust, abs_t_warn = abs_t_warn,
                  abs_bias_warn = abs_bias_warn,
                  group_levels = group_levels,
                  functioning_label = functioning_label_for_facet(facet))
  )
  class(out) <- c("mfrm_dif_interaction", class(out))
  out
}

#' @export
summary.mfrm_dif_interaction <- function(object, ...) {
  out <- list(
    table = object$table,
    summary = object$summary,
    config = object$config
  )
  class(out) <- "summary.mfrm_dif_interaction"
  out
}

#' @export
print.summary.mfrm_dif_interaction <- function(x, ...) {
  label <- x$config$functioning_label %||% "DFF"
  cat("--- ", label, " Interaction Table ---\n", sep = "")
  cat("Facet:", x$config$facet, " | Group:", x$config$group, "\n")
  cat("Groups:", paste(x$config$group_levels, collapse = ", "), "\n")
  cat("Min obs:", x$config$min_obs, " | |t| warn:", x$config$abs_t_warn,
      " | |bias| warn:", x$config$abs_bias_warn, "\n\n")

  cat("Cell Summary:\n")
  print(as.data.frame(x$summary), row.names = FALSE)
  cat("\n")

  if (nrow(x$table) > 0) {
    show_cols <- intersect(
      c("Level", "GroupValue", "N", "ObsExpAvg", "StdResidual",
        "p_adjusted", "sparse", "flag_t", "flag_bias"),
      names(x$table)
    )
    print(as.data.frame(x$table |> select(all_of(show_cols))),
          row.names = FALSE, digits = 3)
  }
  invisible(x)
}

#' @export
print.mfrm_dif_interaction <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

# ============================================================================
# B3. Differential Functioning Heatmap
# ============================================================================

#' Plot a differential-functioning heatmap
#'
#' Visualizes the interaction between a facet and a grouping variable
#' as a heatmap. Rows represent facet levels, columns represent group
#' values, and cell color indicates the selected metric.
#'
#' @param x Output from [dif_interaction_table()], [analyze_dff()], or
#'   [analyze_dif()]. When an `mfrm_dff`/`mfrm_dif` object is passed,
#'   the `cell_table` element
#'   is used (requires `method = "residual"`).
#' @param metric Which metric to plot: `"obs_exp"` for observed-minus-expected
#'   average (default), `"t"` for the standardized residual / t-statistic,
#'   or `"contrast"` for pairwise differential-functioning contrast (only for `mfrm_dff`
#'   objects with `dif_table`).
#' @param draw If `TRUE` (default), draw the plot.
#' @param ... Additional graphical parameters passed to [graphics::image()].
#'
#' @section Interpreting output:
#' - Warm colors (red) indicate positive Obs-Exp values (the model
#'   underestimates the facet level for that group).
#' - Cool colors (blue) indicate negative Obs-Exp values (the model
#'   overestimates).
#' - White/neutral indicates no systematic difference.
#' - The `"contrast"` view is best for pairwise differential-functioning
#'   summaries, whereas
#'   `"obs_exp"` and `"t"` are best for cell-level diagnostics.
#'
#' @section Typical workflow:
#' 1. Compute interaction with [dif_interaction_table()] or differential-
#'    functioning contrasts with [analyze_dff()].
#' 2. Plot with `plot_dif_heatmap(...)`.
#' 3. Identify extreme cells or contrasts for follow-up.
#'
#' @return Invisibly, the matrix used for plotting.
#'
#' @seealso [dif_interaction_table()], [analyze_dff()], [analyze_dif()], [dif_report()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' int <- dif_interaction_table(fit, diag, facet = "Rater",
#'                              group = "Group", data = toy, min_obs = 2)
#' heat <- plot_dif_heatmap(int, metric = "obs_exp", draw = FALSE)
#' dim(heat)
#' @export
plot_dif_heatmap <- function(x, metric = c("obs_exp", "t", "contrast"),
                             draw = TRUE, ...) {
  metric <- match.arg(metric)

  # Resolve input: accept mfrm_dif_interaction or mfrm_dff/mfrm_dif
  if (inherits(x, "mfrm_dif_interaction")) {
    tbl <- x$table
    value_col <- switch(metric,
      obs_exp = "ObsExpAvg",
      t       = "StdResidual",
      contrast = {
        stop("metric = 'contrast' requires an `mfrm_dff`/`mfrm_dif` object with `dif_table`.",
             call. = FALSE)
      }
    )
    row_var <- "Level"
    col_var <- "GroupValue"
  } else if (inherits(x, "mfrm_dif")) {
    if (metric == "contrast") {
      tbl <- x$dif_table
      if (is.null(tbl) || nrow(tbl) == 0) {
        stop("No differential-functioning contrasts available.", call. = FALSE)
      }
      value_col <- "Contrast"
      row_var <- "Level"
      # For contrast, pivot: rows = Level, columns = Group pairs
      tbl$col_label <- paste0(tbl$Group1, " vs ", tbl$Group2)
      col_var <- "col_label"
    } else {
      # Use cell_table
      tbl <- x$cell_table
      if (is.null(tbl) || nrow(tbl) == 0) {
        stop("No cell_table available. Use method = 'residual' in analyze_dff().",
             call. = FALSE)
      }
      value_col <- switch(metric,
        obs_exp = "ObsExpAvg",
        t       = "StdResidual"
      )
      row_var <- "Level"
      col_var <- "GroupValue"
    }
  } else {
    stop("`x` must be an `mfrm_dif_interaction`, `mfrm_dff`, or `mfrm_dif` object.",
         call. = FALSE)
  }

  # Build matrix
  rows <- sort(unique(as.character(tbl[[row_var]])))
  cols <- sort(unique(as.character(tbl[[col_var]])))
  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(cols),
                dimnames = list(rows, cols))
  row_idx <- match(as.character(tbl[[row_var]]), rows)
  col_idx <- match(as.character(tbl[[col_var]]), cols)
  ok <- !is.na(row_idx) & !is.na(col_idx)
  if (any(ok)) {
    mat[cbind(row_idx[ok], col_idx[ok])] <- tbl[[value_col]][ok]
  }

  if (draw) {
    # Color scale: blue-white-red
    n_colors <- 64
    max_abs <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)
    blue_white_red <- grDevices::colorRampPalette(
      c("steelblue", "white", "firebrick")
    )(n_colors)

    old_par <- graphics::par(mar = c(6, 8, 4, 2), no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    metric_label <- switch(metric,
      obs_exp = "Obs - Exp Average",
      t       = "Standardized Residual (t)",
      contrast = "Differential-Functioning Contrast"
    )
    label <- NULL
    if (!is.null(x$config$functioning_label)) {
      label <- x$config$functioning_label
    } else if (inherits(x, "mfrm_dif_interaction")) {
      label <- "DFF"
    }

    graphics::image(
      x = seq_len(ncol(mat)),
      y = seq_len(nrow(mat)),
      z = t(mat),
      col = blue_white_red,
      breaks = breaks,
      axes = FALSE,
      xlab = "", ylab = "",
      main = paste(label %||% "DFF", "Heatmap:", metric_label),
      ...
    )
    graphics::axis(1, at = seq_len(ncol(mat)), labels = cols,
                   las = 2, cex.axis = 0.8)
    graphics::axis(2, at = seq_len(nrow(mat)), labels = rows,
                   las = 1, cex.axis = 0.8)
    graphics::box()

    # Add cell text
    for (ri in seq_len(nrow(mat))) {
      for (ci in seq_len(ncol(mat))) {
        val <- mat[ri, ci]
        if (is.finite(val)) {
          graphics::text(ci, ri, sprintf("%.2f", val), cex = 0.6)
        }
      }
    }
  }

  invisible(mat)
}

# ============================================================================
# C. Information Function Computation and Plotting
# ============================================================================

#' Compute design-weighted precision curves for RSM fits
#'
#' Calculates design-weighted score-variance curves across the latent
#' trait (theta) for a fitted RSM many-facet Rasch model. Returns both
#' an overall precision curve (`$tif`) and per-facet-level contribution
#' curves (`$iif`) based on the realized observation pattern.
#'
#' @param fit Output from [fit_mfrm()].
#' @param theta_range Numeric vector of length 2 giving the range of theta
#'   values. Default `c(-6, 6)`.
#' @param theta_points Integer number of points at which to evaluate
#'   information. Default `201`.
#'
#' @details
#' For a polytomous Rasch model with K+1 categories, the score variance at
#' theta for one observed design cell is:
#' \deqn{I(\theta) = \sum_{k=0}^{K} P_k(\theta) \left(k - E(\theta)\right)^2}
#' where \eqn{P_k} is the category probability and \eqn{E(\theta)} is the
#' expected score at theta. In `mfrmr`, these cell-level variances are then
#' aggregated with weights taken from the realized observation counts in
#' `fit$prep$data`.
#'
#' The resulting total curve is therefore a design-weighted precision screen
#' rather than a pure textbook test-information function for an abstract fixed
#' item set. The associated standard error summary is still
#' \eqn{SE(\theta) = 1 / \sqrt{I(\theta)}} for positive information values.
#'
#' @section What `tif` and `iif` mean here:
#' In `mfrmr`, this helper currently supports only RSM fits. The total
#' curve (`$tif`) is the sum of design-weighted cell contributions across all
#' non-person facet levels in the fitted model. The facet-level contribution
#' curves (`$iif`) keep those weighted contributions separated, so you can see
#' which observed rater levels, criteria, or other facet levels are driving
#' precision at different parts of the scale.
#'
#' @section What this quantity does not justify:
#' - It is not a textbook many-facet test-information function for an abstract
#'   fixed item set.
#' - It should not be used as if it were design-free evidence about a form's
#'   precision independent of the realized observation pattern.
#' - It does not currently extend to PCM fits; the helper stops for
#'   `model = "PCM"`.
#'
#' @section When to use this:
#' Use `compute_information()` when you want a design-weighted precision screen
#' for an RSM fit along the latent continuum. In practice:
#' - start with the total precision curve for overall targeting across the
#'   realized observation pattern
#' - inspect facet-level contribution curves when you want to see which raters,
#'   criteria, or other facet levels account for more of that design-weighted
#'   precision
#' - widen `theta_range` if you expect extreme measures and want to inspect the
#'   tails explicitly
#'
#' @section Choosing the theta grid:
#' The defaults (`theta_range = c(-6, 6)`, `theta_points = 201`) work well for
#' routine inspection. Expand the range if person or facet measures extend into
#' the tails, and increase `theta_points` only when you need a smoother grid
#' for reporting or custom graphics.
#'
#' @section Interpreting output:
#' - `$tif`: design-weighted precision curve data with theta, Information, and SE.
#' - `$iif`: design-weighted facet-level contribution curves for an RSM fit.
#' - Higher information implies more precise measurement at that theta.
#' - SE is inversely related to information.
#' - Peaks in the total curve show the trait region where the realized
#'   calibration is most informative.
#' - Facet-level curves help explain *which observed facet levels* contribute
#'   to those peaks; they are not standalone item-information curves and should
#'   be read as design contributions.
#'
#' @section How to read the main columns:
#' - `Theta`: point on the latent continuum where the curve is evaluated.
#' - `Information`: design-weighted precision value at that theta.
#' - `SE`: approximate `1 / sqrt(Information)` summary for positive values.
#' - `Exposure`: total realized observation weight contributing to a facet-level
#'   curve in `$iif`.
#'
#' @section Recommended next step:
#' Compare the precision peak with person/facet locations from a Wright map or
#' related diagnostics. If you need to decide how strongly SE/CI language can
#' be used in reporting, follow with [precision_audit_report()].
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Run `compute_information(fit)`.
#' 3. Plot with `plot_information(info, type = "tif")`.
#' 4. If needed, inspect facet contributions with
#'    `plot_information(info, type = "iif", facet = "Rater")`.
#'
#' @return
#' An object of class `mfrm_information` (named list) with:
#' - `tif`: tibble with columns `Theta`, `Information`, `SE`. The
#'   `Information` column stores the design-weighted precision value.
#' - `iif`: tibble with columns `Theta`, `Facet`, `Level`, `Information`,
#'   and `Exposure`. Here too, `Information` stores a design-weighted
#'   contribution value retained under that column name for compatibility.
#' - `theta_range`: the evaluated theta range.
#'
#' @seealso [fit_mfrm()], [plot_information()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 25)
#' info <- compute_information(fit)
#' head(info$tif)
#' info$tif$Theta[which.max(info$tif$Information)]
#' @export
compute_information <- function(fit,
                                theta_range = c(-6, 6),
                                theta_points = 201L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an `mfrm_fit` object.", call. = FALSE)
  }
  model <- as.character(fit$config$model %||% NA_character_)
  if (!identical(model, "RSM")) {
    stop(
      "`compute_information()` currently supports only `model = \"RSM\"` fits. ",
      "PCM information requires step-facet-specific thresholds and is not yet implemented.",
      call. = FALSE
    )
  }

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)

  # Extract model parameters and realized observation design.
  steps <- fit$steps
  if (is.null(steps) || nrow(steps) == 0 || !"Estimate" %in% names(steps)) {
    stop("Step/threshold estimates are required for information computation.",
         call. = FALSE)
  }
  step_est <- steps$Estimate
  step_cum <- c(0, cumsum(step_est))
  categories <- 0:length(step_est)
  category_vec <- matrix(categories, ncol = 1)
  category_sq_vec <- matrix(categories^2, ncol = 1)

  facet_tbl <- tibble::as_tibble(fit$facets$others)
  if (nrow(facet_tbl) == 0) {
    stop("Facet estimates are required for information computation.",
         call. = FALSE)
  }
  obs_df <- as.data.frame(fit$prep$data %||% NULL, stringsAsFactors = FALSE)
  if (nrow(obs_df) == 0) {
    stop("Prepared observation data are required for information computation.",
         call. = FALSE)
  }
  facet_names <- as.character(fit$config$facet_names %||% unique(as.character(facet_tbl$Facet)))
  facet_names <- facet_names[facet_names %in% names(obs_df)]
  if (length(facet_names) == 0) {
    stop("Facet columns were not found in the prepared response data.",
         call. = FALSE)
  }
  facet_signs <- fit$config$facet_signs %||% stats::setNames(rep(-1, length(facet_names)), facet_names)
  facet_signs <- facet_signs[facet_names]
  facet_signs[!is.finite(facet_signs)] <- -1

  # Design-weighted information for a single observed design cell at each theta.
  compute_cell_info <- function(offset) {
    probs <- category_prob_rsm(theta_grid + offset, step_cum)
    expected <- as.vector(probs %*% category_vec)
    second_moment <- as.vector(probs %*% category_sq_vec)
    pmax(second_moment - expected^2, 0)
  }

  facet_tbl <- facet_tbl[is.finite(facet_tbl$Estimate), , drop = FALSE]
  if (nrow(facet_tbl) == 0) {
    stop("Facet estimates are required for information computation.",
         call. = FALSE)
  }

  obs_weights <- suppressWarnings(as.numeric(obs_df$Weight %||% rep(1, nrow(obs_df))))
  obs_weights[!is.finite(obs_weights)] <- 0
  design_cells <- obs_df[, facet_names, drop = FALSE]
  design_cells$Exposure <- obs_weights
  design_cells <- design_cells |>
    dplyr::group_by(dplyr::across(dplyr::all_of(facet_names))) |>
    dplyr::summarize(Exposure = sum(.data$Exposure, na.rm = TRUE), .groups = "drop")
  design_cells <- as.data.frame(design_cells, stringsAsFactors = FALSE)

  est_key <- paste(facet_tbl$Facet, facet_tbl$Level, sep = "||")
  est_lookup <- stats::setNames(facet_tbl$Estimate, est_key)
  cell_offset <- numeric(nrow(design_cells))
  cell_ok <- design_cells$Exposure > 0
  for (facet in facet_names) {
    keys <- paste(facet, design_cells[[facet]], sep = "||")
    est_vals <- suppressWarnings(as.numeric(est_lookup[keys]))
    sign_val <- suppressWarnings(as.numeric(facet_signs[[facet]]))
    if (!is.finite(sign_val)) sign_val <- -1
    cell_ok <- cell_ok & is.finite(est_vals)
    cell_offset <- cell_offset + sign_val * est_vals
  }
  design_cells <- design_cells[cell_ok, , drop = FALSE]
  cell_offset <- cell_offset[cell_ok]
  if (nrow(design_cells) == 0) {
    stop("No valid observed design cells were available for information computation.",
         call. = FALSE)
  }

  info_mat <- vapply(cell_offset, compute_cell_info, numeric(length(theta_grid)))
  if (!is.matrix(info_mat)) {
    info_mat <- matrix(info_mat, ncol = 1)
  }
  weighted_info_mat <- sweep(info_mat, 2, design_cells$Exposure, `*`)

  total_info <- rowSums(weighted_info_mat)
  iif_rows <- lapply(seq_len(nrow(facet_tbl)), function(i) {
    facet_i <- as.character(facet_tbl$Facet[i])
    level_i <- as.character(facet_tbl$Level[i])
    mask <- as.character(design_cells[[facet_i]]) == level_i
    tibble(
      Theta = theta_grid,
      Facet = facet_i,
      Level = level_i,
      Information = if (any(mask)) rowSums(weighted_info_mat[, mask, drop = FALSE]) else 0,
      Exposure = if (any(mask)) sum(design_cells$Exposure[mask], na.rm = TRUE) else 0
    )
  })
  iif <- dplyr::bind_rows(iif_rows)
  tif <- tibble(
    Theta = theta_grid,
    Information = total_info,
    SE = ifelse(total_info > 0, 1 / sqrt(total_info), NA_real_)
  )

  out <- list(tif = tif, iif = iif, theta_range = theta_range)
  class(out) <- c("mfrm_information", class(out))
  out
}

#' Plot design-weighted precision curves
#'
#' Visualize the design-weighted precision curve and optionally
#' per-facet-level contribution curves from [compute_information()].
#'
#' @param x Output from [compute_information()].
#' @param type `"tif"` for the overall precision curve (default), `"iif"` for
#'   facet-level contribution curves, `"se"` for the approximate standard error
#'   implied by that curve, or `"both"` for precision with approximate SE on a
#'   secondary axis.
#' @param facet For `type = "iif"`, which facet to plot. If `NULL`,
#'   the first facet is used.
#' @param draw If `TRUE` (default), draw the plot. If `FALSE`, return
#'   the plot data invisibly.
#' @param ... Additional graphical parameters.
#'
#' @section Plot types:
#' - `"tif"`: overall design-weighted precision across theta.
#' - `"se"`: approximate standard error across theta.
#' - `"both"`: precision and approximate SE together, useful for presentations.
#' - `"iif"`: facet-level contribution curves for one selected facet in an
#'   RSM fit.
#'
#' @section Which type should I use?:
#' - Use `"tif"` for a quick overall read on precision.
#' - Use `"se"` when standard-error language is easier to communicate than
#'   precision.
#' - Use `"both"` when you want both views in one figure.
#' - Use `"iif"` when you want to see which facet levels are shaping the total
#'   precision curve.
#'
#' @section Interpreting output:
#' - The total curve peaks where the realized design is most precise.
#' - SE is derived as `1 / sqrt(precision)`; lower is better.
#' - Facet-level curves show which facet levels contribute most to that
#'   realized precision at each theta.
#' - If the precision peak sits far from the bulk of person measures, the
#'   realized design may be poorly targeted.
#'
#' @section Returned data when draw = FALSE:
#' For `type = "tif"`, `"se"`, or `"both"`, the returned data come from
#' `x$tif`. For `type = "iif"`, the returned data are the rows of `x$iif`
#' filtered to the requested facet.
#'
#' @section Typical workflow:
#' 1. Compute information with [compute_information()].
#' 2. Plot with `plot_information(info)` for the total precision curve.
#' 3. Use `plot_information(info, type = "iif", facet = "Rater")` for
#'    facet-level contributions.
#' 4. Use `draw = FALSE` when you want the plotting data for custom graphics.
#'
#' @return Invisibly, the plot data (tibble).
#'
#' @seealso [compute_information()], [fit_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 25)
#' info <- compute_information(fit)
#' tif_data <- plot_information(info, type = "tif", draw = FALSE)
#' head(tif_data)
#' iif_data <- plot_information(info, type = "iif", facet = "Rater", draw = FALSE)
#' head(iif_data)
#' @export
plot_information <- function(x,
                             type = c("tif", "iif", "se", "both"),
                             facet = NULL,
                             draw = TRUE,
                             ...) {
  if (!inherits(x, "mfrm_information")) {
    stop("`x` must be an `mfrm_information` object.", call. = FALSE)
  }
  type <- match.arg(type)

  if (type == "tif" || type == "both") {
    plot_data <- x$tif
    if (draw) {
      if (type == "both") {
        par_old <- graphics::par(mar = c(5, 4, 4, 4) + 0.1)
        on.exit(graphics::par(par_old), add = TRUE)
      }
      plot(plot_data$Theta, plot_data$Information,
           type = "l", lwd = 2, col = "steelblue",
           xlab = expression(theta), ylab = "Information (precision)",
           main = "Design-Weighted Precision Curve", ...)
      graphics::grid()
      if (type == "both") {
        graphics::par(new = TRUE)
        plot(plot_data$Theta, plot_data$SE,
             type = "l", lwd = 2, col = "coral", lty = 2,
             axes = FALSE, xlab = "", ylab = "")
        graphics::axis(4, col = "coral", col.axis = "coral")
        graphics::mtext("Approx. SE", side = 4, line = 2.5, col = "coral")
        graphics::legend("topright",
                         legend = c("Information (precision)", "Approx. SE"),
                         col = c("steelblue", "coral"),
                         lty = c(1, 2), lwd = 2, bty = "n")
      }
    }
    invisible(plot_data)
  } else if (type == "se") {
    plot_data <- x$tif
    if (draw) {
      plot(plot_data$Theta, plot_data$SE,
           type = "l", lwd = 2, col = "coral",
           xlab = expression(theta), ylab = "Approx. SE",
           main = "Approx. SE from Design-Weighted Precision", ...)
      graphics::grid()
    }
    invisible(plot_data)
  } else {
    # Facet-level contribution curves
    iif <- x$iif
    if (is.null(facet)) {
      facet <- unique(iif$Facet)[1]
    }
    plot_data <- iif |> filter(.data$Facet == facet)
    if (nrow(plot_data) == 0) {
      stop("No information data for facet '", facet, "'.", call. = FALSE)
    }
    if (draw) {
      levels_u <- unique(plot_data$Level)
      n_lev <- length(levels_u)
      cols <- grDevices::rainbow(n_lev, s = 0.7, v = 0.8)
      yr <- range(plot_data$Information, na.rm = TRUE)
      plot(NA, xlim = range(plot_data$Theta), ylim = yr,
           xlab = expression(theta), ylab = "Information contribution",
           main = paste("Facet-Level Precision Contributions:", facet), ...)
      for (k in seq_along(levels_u)) {
        sub <- plot_data |> filter(.data$Level == levels_u[k])
        graphics::lines(sub$Theta, sub$Information, col = cols[k], lwd = 1.5)
      }
      graphics::legend("topright", legend = levels_u, col = cols,
                       lty = 1, lwd = 1.5, bty = "n", cex = 0.8)
      graphics::grid()
    }
    invisible(plot_data)
  }
}

# ============================================================================
# D. Unified Wright Map (persons + all facets on shared logit scale)
# ============================================================================

#' Plot a unified Wright map with all facets on a shared logit scale
#'
#' Produces a shared-logit variable map showing person ability distribution
#' alongside measure estimates for every facet in side-by-side columns on
#' the same scale.
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param bins Integer number of bins for the person histogram. Default `20`.
#' @param show_thresholds Logical; if `TRUE`, display threshold/step
#'   positions on the map. Default `TRUE`.
#' @param top_n Maximum number of facet/step points retained for labeling.
#' @param show_ci Logical; if `TRUE`, draw approximate confidence intervals when
#'   standard errors are available.
#' @param ci_level Confidence level used when `show_ci = TRUE`.
#' @param draw If `TRUE` (default), draw the plot. If `FALSE`, return
#'   plot data invisibly.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`).
#' @param palette Optional named color overrides passed to the shared Wright-map
#'   drawer.
#' @param label_angle Rotation angle for group labels on the facet panel.
#' @param ... Additional graphical parameters.
#'
#' @details
#' This unified map arranges:
#' - Column 1: Person measure distribution (horizontal histogram)
#' - Shared facet/step panel: facet levels and optional threshold positions on
#'   the same vertical logit axis
#' - Range and interquartile overlays for each facet group to show spread
#'
#' This is the package's most compact targeting view when you want one display
#' that shows where persons, facet levels, and category thresholds sit
#' relative to the same latent scale.
#'
#' The logit scale on the y-axis is shared, allowing direct visual
#' comparison of all facets and persons.
#'
#' @section Interpreting output:
#' - Facet levels at the same height on the map are at similar difficulty.
#' - The person histogram shows where examinees cluster relative to the
#'   facet scale.
#' - Thresholds (if shown) indicate category boundary positions.
#' - Large gaps between the person distribution and facet locations can signal
#'   targeting problems.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Plot with `plot_wright_unified(fit)`.
#' 3. Compare person distribution with facet level locations.
#' 4. Use `show_thresholds = TRUE` when you want the category structure in the
#'    same view.
#'
#' @section When to use this instead of plot_information:
#' Use `plot_wright_unified()` when your main question is targeting or coverage
#' on the shared logit scale. Use [plot_information()] when your main question
#' is measurement precision across theta.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @return Invisibly, a list with `persons`, `facets`, and `thresholds`
#'   data used for the plot.
#'
#' @seealso [fit_mfrm()], [plot.mfrm_fit()], [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(toy_small, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 10)
#' map_data <- plot_wright_unified(fit, draw = FALSE)
#' names(map_data)
#' @export
plot_wright_unified <- function(fit,
                                diagnostics = NULL,
                                bins = 20L,
                                show_thresholds = TRUE,
                                top_n = 30L,
                                show_ci = FALSE,
                                ci_level = 0.95,
                                draw = TRUE,
                                preset = c("standard", "publication", "compact"),
                                palette = NULL,
                                label_angle = 45,
                                ...) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an `mfrm_fit` object.", call. = FALSE)
  }
  top_n <- max(1L, as.integer(top_n))
  bins <- max(5L, as.integer(bins))
  style <- resolve_plot_preset(preset)
  se_tbl_ci <- if (isTRUE(show_ci)) compute_se_for_plot(fit, ci_level = ci_level) else NULL
  plot_core <- build_wright_map_data(
    fit,
    top_n = top_n,
    se_tbl = se_tbl_ci,
    include_steps = isTRUE(show_thresholds)
  )
  plot_core$person_hist <- graphics::hist(plot_core$person$Estimate, breaks = bins, plot = FALSE)
  plot_data <- c(
    list(
      persons = plot_core$person$Estimate,
      facets = tibble::as_tibble(fit$facets$others),
      thresholds = if (isTRUE(show_thresholds) && !is.null(fit$steps)) tibble::as_tibble(fit$steps) else NULL,
      facet_names = unique(as.character(fit$facets$others$Facet)),
      y_lim = plot_core$y_range
    ),
    plot_core
  )
  if (!draw) return(invisible(plot_data))

  apply_plot_preset(style)
  draw_wright_map(
    plot_core,
    title = "Unified Wright Map",
    palette = resolve_palette(
      palette = palette,
      defaults = c(
        facet_level = style$accent_tertiary,
        step_threshold = style$accent_secondary,
        person_hist = style$fill_muted,
        grid = style$grid,
        range = style$accent_primary,
        iqr = style$foreground
      )
    ),
    label_angle = label_angle,
    show_ci = show_ci,
    ci_level = ci_level
  )

  invisible(plot_data)
}

# ============================================================================
# E. Anchoring & Equating Workflow (Phase 4)
# ============================================================================

# --- Internal helper: compute drift between anchored fit and baseline --------
measure_se_table <- function(fit, include_person = FALSE, diagnostics = NULL) {
  if (is.null(diagnostics)) {
    if (is.null(fit)) {
      stop("`fit` or `diagnostics` must be supplied to `measure_se_table()`.", call. = FALSE)
    }
    diagnostics <- diagnose_mfrm(fit)
  }
  measures <- tibble::as_tibble(diagnostics$measures)
  if (!isTRUE(include_person)) {
    measures <- measures |>
      dplyr::filter(.data$Facet != "Person")
  }
  measures |>
    dplyr::transmute(
      Facet = as.character(.data$Facet),
      Level = as.character(.data$Level),
      SE = as.numeric(.data$SE)
    ) |>
    dplyr::distinct()
}

compute_equating_offset <- function(diffs, se_from = NULL, se_to = NULL,
                                    drift_threshold = NULL) {
  diffs <- as.numeric(diffs)
  ok <- is.finite(diffs)
  if (!any(ok)) {
    return(list(
      offset_prelim = NA_real_,
      offset = NA_real_,
      residual = rep(NA_real_, length(diffs)),
      retained = rep(FALSE, length(diffs)),
      n_retained = 0L,
      weighting = "none"
    ))
  }

  se_from <- if (is.null(se_from)) rep(NA_real_, length(diffs)) else as.numeric(se_from)
  se_to <- if (is.null(se_to)) rep(NA_real_, length(diffs)) else as.numeric(se_to)
  weight_ok <- ok & is.finite(se_from) & is.finite(se_to) & se_from > 0 & se_to > 0
  weights <- ifelse(weight_ok, 1 / (se_from^2 + se_to^2), NA_real_)

  offset_prelim <- if (any(weight_ok)) {
    stats::weighted.mean(diffs[weight_ok], w = weights[weight_ok])
  } else {
    mean(diffs[ok])
  }

  residual_prelim <- diffs - offset_prelim
  retained <- ok
  if (!is.null(drift_threshold) && is.finite(drift_threshold)) {
    retained <- retained & abs(residual_prelim) <= drift_threshold
  }
  if (!any(retained)) {
    retained <- ok
  }

  weight_retained <- retained & weight_ok
  offset <- if (any(weight_retained)) {
    stats::weighted.mean(diffs[weight_retained], w = weights[weight_retained])
  } else {
    mean(diffs[retained])
  }

  list(
    offset_prelim = offset_prelim,
    offset = offset,
    residual = diffs - offset,
    retained = retained,
    n_retained = sum(retained, na.rm = TRUE),
    weighting = if (any(weight_retained)) "inverse_variance" else "unweighted"
  )
}

.summarise_link_support <- function(common_tbl,
                                    retained = NULL,
                                    guideline = 5L) {
  if (is.null(common_tbl) || nrow(common_tbl) == 0) {
    return(tibble::tibble(
      Facet = character(),
      N_Common = integer(),
      N_Retained = integer(),
      GuidelineMinCommon = integer(),
      LinkSupportAdequate = logical()
    ))
  }

  base_tbl <- tibble::as_tibble(common_tbl)[, c("Facet", "Level")]
  retained <- if (is.null(retained)) rep(TRUE, nrow(base_tbl)) else as.logical(retained)
  if (length(retained) != nrow(base_tbl)) {
    retained <- rep(FALSE, nrow(base_tbl))
  }

  common_counts <- base_tbl |>
    dplyr::count(.data$Facet, name = "N_Common")
  retained_counts <- base_tbl[retained, , drop = FALSE] |>
    dplyr::count(.data$Facet, name = "N_Retained")

  common_counts |>
    dplyr::left_join(retained_counts, by = "Facet") |>
    dplyr::mutate(
      N_Retained = dplyr::coalesce(.data$N_Retained, 0L),
      GuidelineMinCommon = as.integer(guideline),
      LinkSupportAdequate = .data$N_Retained >= .data$GuidelineMinCommon
    )
}

.compute_drift <- function(fit, anchor_tbl, diagnostics = NULL, baseline_diagnostics = NULL) {
  # Get new estimates
  new_est <- make_anchor_table(fit, include_person = FALSE)

  # Join with baseline anchors
  joined <- dplyr::inner_join(
    anchor_tbl |> dplyr::rename(Baseline = "Anchor"),
    new_est |> dplyr::rename(New = "Anchor"),
    by = c("Facet", "Level")
  )

  if (nrow(joined) == 0) {
    return(tibble::tibble(
      Facet = character(), Level = character(), Baseline = numeric(),
      New = numeric(), Drift = numeric(), SE_Baseline = numeric(),
      SE_New = numeric(), SE_Diff = numeric(),
      Drift_SE_Ratio = numeric(), Flag = logical()
    ))
  }

  baseline_se <- measure_se_table(
    fit = NULL,
    include_person = FALSE,
    diagnostics = baseline_diagnostics
  )
  new_se <- measure_se_table(fit, include_person = FALSE, diagnostics = diagnostics)

  joined <- joined |>
    dplyr::left_join(
      baseline_se |> dplyr::rename(SE_Baseline = "SE"),
      by = c("Facet", "Level")
    ) |>
    dplyr::left_join(
      new_se |> dplyr::rename(SE_New = "SE"),
      by = c("Facet", "Level")
    )

  joined |>
    dplyr::mutate(
      Drift = .data$New - .data$Baseline,
      SE_Diff = ifelse(
        is.finite(.data$SE_Baseline) & is.finite(.data$SE_New),
        sqrt(.data$SE_Baseline^2 + .data$SE_New^2),
        NA_real_
      ),
      Drift_SE_Ratio = ifelse(
        is.na(.data$SE_Diff) | .data$SE_Diff == 0,
        NA_real_,
        abs(.data$Drift) / .data$SE_Diff
      ),
      Flag = abs(.data$Drift) > 0.5 | (!is.na(.data$Drift_SE_Ratio) & .data$Drift_SE_Ratio > 2)
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$Drift)))
}

# --- anchor_to_baseline ------------------------------------------------------

#' Fit new data anchored to a baseline calibration
#'
#' Re-estimates a many-facet Rasch model on new data while holding selected
#' facet parameters fixed at the values from a previous (baseline) calibration.
#' This is the standard workflow for placing new data onto an existing scale,
#' linking test forms, or carrying a baseline calibration across
#' administration windows.
#'
#' @param new_data Data frame in long format (one row per rating).
#' @param baseline_fit An `mfrm_fit` object from a previous calibration.
#' @param person Character column name for person/examinee.
#' @param facets Character vector of facet column names.
#' @param score Character column name for the rating score.
#' @param anchor_facets Character vector of facets to anchor (default: all
#'   non-Person facets).
#' @param include_person If `TRUE`, also anchor person estimates.
#' @param weight Optional character column name for observation weights.
#' @param model Scale model override; defaults to baseline model.
#' @param method Estimation method override; defaults to baseline method.
#' @param anchor_policy How to handle anchor issues: `"warn"`, `"error"`,
#'   `"silent"`.
#' @param ... Additional arguments passed to [fit_mfrm()].
#'
#' @details
#' This function automates the baseline-anchored calibration workflow:
#'
#' 1. Extracts anchor values from the baseline fit using [make_anchor_table()].
#' 2. Re-estimates the model on `new_data` with those anchors fixed via
#'    `fit_mfrm(..., anchors = anchor_table)`.
#' 3. Runs [diagnose_mfrm()] on the anchored fit.
#' 4. Computes element-level differences (new estimate minus baseline
#'    estimate) for every common element.
#'
#' The `model` and `method` arguments default to the baseline fit's settings
#' so the calibration framework remains consistent.  Elements present in the
#' anchor table but absent from the new data are handled according to
#' `anchor_policy`: `"warn"` (default) emits a message, `"error"` stops
#' execution, and `"silent"` ignores silently.
#'
#' The returned `drift` table is best interpreted as an anchored consistency
#' check. When a facet is fixed through `anchor_facets`, those anchored levels
#' are constrained in the new run, so their reported differences are not an
#' independent drift analysis. For genuine cross-wave drift monitoring, fit the
#' waves separately and use [detect_anchor_drift()] on the resulting fits.
#'
#' Element-level differences are calculated for every element that appears in
#' both the baseline and the new calibration:
#' \deqn{\Delta_e = \hat{\delta}_{e,\text{new}} - \hat{\delta}_{e,\text{base}}}
#' An element is **flagged** when \eqn{|\Delta_e| > 0.5} logits or
#' \eqn{|\Delta_e / SE_{\Delta_e}| > 2.0}, where
#' \eqn{SE_{\Delta_e} = \sqrt{SE_{\mathrm{base}}^2 + SE_{\mathrm{new}}^2}}.
#'
#' @section Which function should I use?:
#' - Use `anchor_to_baseline()` when you have one new dataset and want to place
#'   it directly on a baseline scale.
#' - Use [detect_anchor_drift()] when you already have multiple fitted waves
#'   and want to compare their stability.
#' - Use [build_equating_chain()] when you need cumulative offsets across an
#'   ordered series of waves.
#'
#' @section Interpreting output:
#' - `$drift`: one row per common element with columns `Facet`, `Level`,
#'   `Baseline`, `New`, `Drift`, `SE_Baseline`, `SE_New`, `SE_Diff`,
#'   `Drift_SE_Ratio`, and `Flag`.
#'   Read this as an anchored consistency table. Small absolute differences
#'   indicate that the anchored re-fit stayed close to the baseline scale.
#'   Flagged rows warrant review, but they are not a substitute for a separate
#'   drift study on unanchored common elements.
#' - `$fit`: the full anchored `mfrm_fit` object, usable with
#'   [diagnose_mfrm()], [measurable_summary_table()], etc.
#' - `$diagnostics`: pre-computed diagnostics for the anchored calibration.
#' - `$baseline_anchors`: the anchor table fed to [fit_mfrm()], useful for
#'   auditing which elements were constrained.
#'
#' @section Typical workflow:
#' 1. Fit the baseline model: `fit1 <- fit_mfrm(...)`.
#' 2. Collect new data (e.g., a later administration).
#' 3. Call `res <- anchor_to_baseline(new_data, fit1, ...)`.
#' 4. Inspect `summary(res)` to confirm the anchored run remains close to the
#'    baseline scale.
#' 5. For multi-wave drift monitoring, fit waves separately and pass the fits to
#'    [detect_anchor_drift()] or [build_equating_chain()].
#'
#' @return Object of class `mfrm_anchored_fit` with components:
#'   \describe{
#'     \item{fit}{The anchored `mfrm_fit` object.}
#'     \item{diagnostics}{Output of [diagnose_mfrm()] on the anchored fit.}
#'     \item{baseline_anchors}{Anchor table extracted from the baseline.}
#'     \item{drift}{Tibble of element-level drift statistics.}
#'   }
#'
#' @seealso [fit_mfrm()], [make_anchor_table()], [detect_anchor_drift()],
#'   [diagnose_mfrm()], [build_equating_chain()], [mfrmr_linking_and_dff]
#' @export
#' @examples
#' \donttest{
#' d1 <- load_mfrmr_data("study1")
#' keep1 <- unique(d1$Person)[1:15]
#' d1 <- d1[d1$Person %in% keep1, , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 15)
#' d2 <- load_mfrmr_data("study2")
#' keep2 <- unique(d2$Person)[1:15]
#' d2 <- d2[d2$Person %in% keep2, , drop = FALSE]
#' res <- anchor_to_baseline(d2, fit1, "Person",
#'                           c("Rater", "Criterion"), "Score",
#'                           anchor_facets = "Criterion")
#' summary(res)
#' head(res$drift[, c("Facet", "Level", "Drift", "Flag")])
#' res$baseline_anchors[1:3, ]
#' }
anchor_to_baseline <- function(new_data, baseline_fit,
                               person, facets, score,
                               anchor_facets = NULL,
                               include_person = FALSE,
                               weight = NULL,
                               model = NULL, method = NULL,
                               anchor_policy = "warn",
                               ...) {
  # Validate baseline_fit
  stopifnot(inherits(baseline_fit, "mfrm_fit"))

  # Inherit model/method from baseline if not specified
  if (is.null(model))  model  <- baseline_fit$config$model
  if (is.null(method)) method <- baseline_fit$config$method

  # Extract anchor table from baseline
  anchor_tbl <- make_anchor_table(baseline_fit, facets = anchor_facets,
                                   include_person = include_person)

  if (nrow(anchor_tbl) == 0) {
    stop("No anchors could be extracted from the baseline fit.", call. = FALSE)
  }

  # Fit new data with anchors
  new_fit <- fit_mfrm(new_data, person = person, facets = facets, score = score,
                      weight = weight, model = model, method = method,
                      anchors = anchor_tbl, anchor_policy = anchor_policy, ...)

  # Compute diagnostics
  baseline_diag <- diagnose_mfrm(baseline_fit)
  new_diag <- diagnose_mfrm(new_fit)

  # Compute drift: compare new estimates to baseline anchors for common elements
  drift <- .compute_drift(
    new_fit,
    anchor_tbl,
    diagnostics = new_diag,
    baseline_diagnostics = baseline_diag
  )

  out <- list(
    fit = new_fit,
    diagnostics = new_diag,
    baseline_anchors = anchor_tbl,
    drift = drift
  )
  class(out) <- c("mfrm_anchored_fit", "list")
  out
}

#' @rdname anchor_to_baseline
#' @param x An `mfrm_anchored_fit` object.
#' @param ... Ignored.
#' @export
print.mfrm_anchored_fit <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' @rdname anchor_to_baseline
#' @param object An `mfrm_anchored_fit` object (for `summary`).
#' @export
summary.mfrm_anchored_fit <- function(object, ...) {
  drift <- object$drift
  n_anchored <- nrow(object$baseline_anchors)
  n_common <- nrow(drift)
  n_flagged <- sum(drift$Flag, na.rm = TRUE)

  out <- list(
    n_anchored = n_anchored, n_common = n_common, n_flagged = n_flagged,
    drift_summary = if (n_common > 0) {
      drift |> dplyr::group_by(.data$Facet) |>
        dplyr::summarise(N = dplyr::n(), Mean_Drift = mean(abs(.data$Drift)),
                         Max_Drift = max(abs(.data$Drift)),
                         N_Flagged = sum(.data$Flag, na.rm = TRUE),
                         .groups = "drop")
    } else {
      tibble::tibble()
    },
    flagged = drift |> dplyr::filter(.data$Flag),
    converged = object$fit$summary$Converged
  )
  class(out) <- "summary.mfrm_anchored_fit"
  out
}

#' @rdname anchor_to_baseline
#' @export
print.summary.mfrm_anchored_fit <- function(x, ...) {
  cat("--- Anchored Fit Summary ---\n")
  cat("Converged:", x$converged, "\n")
  cat("Anchors used:", x$n_anchored, "| Common elements:", x$n_common,
      "| Flagged:", x$n_flagged, "\n\n")
  if (nrow(x$drift_summary) > 0) {
    cat("Drift by facet:\n")
    print(as.data.frame(x$drift_summary), row.names = FALSE, digits = 3)
  }
  if (nrow(x$flagged) > 0) {
    cat("\nFlagged elements (|Drift| > 0.5 or |Drift|/SE > 2):\n")
    print(as.data.frame(x$flagged), row.names = FALSE, digits = 3)
  }
  invisible(x)
}

# --- detect_anchor_drift -----------------------------------------------------

#' Detect anchor drift across multiple calibrations
#'
#' Compares facet estimates across two or more calibration waves to identify
#' elements whose difficulty/severity has shifted beyond acceptable thresholds.
#' Useful for monitoring rater drift over time or checking the stability of
#' item banks.
#'
#' @param fits Named list of `mfrm_fit` objects (e.g.,
#'   `list(Year1 = fit1, Year2 = fit2)`).
#' @param facets Character vector of facets to compare (default: all
#'   non-Person facets).
#' @param drift_threshold Absolute drift threshold for flagging (logits,
#'   default 0.5).
#' @param flag_se_ratio Drift/SE ratio threshold for flagging (default 2.0).
#' @param reference Index or name of the reference fit (default: first).
#' @param include_person Include person estimates in comparison.
#'
#' @details
#' For each non-reference wave, the function extracts facet-level estimates
#' using [make_anchor_table()] and computes the element-by-element difference
#' against the reference wave.  Standard errors are obtained from
#' [diagnose_mfrm()] applied to each fit.  Only elements common to both the
#' reference and a comparison wave are included. Before reporting drift, the
#' function removes the weighted common-element link offset between the two
#' waves so that `Drift` represents residual instability rather than the
#' overall shift between calibrations. The function also records how many
#' common elements survive the screening step within each linking facet and
#' treats fewer than 5 retained common elements per facet as thin support.
#'
#' An element is **flagged** when either condition is met:
#' \deqn{|\Delta_e| > \texttt{drift\_threshold}}
#' \deqn{|\Delta_e / SE_{\Delta_e}| > \texttt{flag\_se\_ratio}}
#' The dual-criterion approach guards against flagging elements with large
#' but imprecise estimates, and against missing small but precisely estimated
#' shifts.
#'
#' When `facets` is `NULL`, all non-Person facets are compared.  Providing a
#' subset (e.g., `facets = "Criterion"`) restricts comparison to those facets
#' only.
#'
#' @section Which function should I use?:
#' - Use [anchor_to_baseline()] when your starting point is raw new data plus a
#'   single baseline fit.
#' - Use `detect_anchor_drift()` when you already have multiple fitted waves
#'   and want a reference-versus-wave comparison.
#' - Use [build_equating_chain()] when the waves form a sequence and you need
#'   cumulative linking offsets.
#'
#' @section Interpreting output:
#' - `$drift_table`: one row per element x wave combination, with columns
#'   `Facet`, `Level`, `Wave`, `Ref_Est`, `Wave_Est`, `LinkOffset`, `Drift`,
#'   `SE_Ref`, `SE_Wave`, `SE`, `Drift_SE_Ratio`, `LinkSupportAdequate`, and
#'   `Flag`.  Large drift signals instability after alignment to the
#'   common-element link.
#' - `$summary`: aggregated statistics by facet and wave: number of elements,
#'   mean/max absolute drift, and count of flagged elements.
#' - `$common_elements`: pairwise common-element counts in tidy table form.
#'   Small
#'   overlap weakens the comparison and results should be interpreted
#'   cautiously.
#' - `$common_by_facet`: retained common-element counts by linking facet for
#'   each reference-vs-wave comparison. `LinkSupportAdequate = FALSE` means the
#'   link rests on fewer than 5 retained common elements in at least one facet.
#' - `$config`: records the analysis parameters for reproducibility.
#' - A practical reading order is `summary(drift)` first, then
#'   `drift$drift_table`, then `drift$common_by_facet` if overlap looks thin.
#'
#' @section Typical workflow:
#' 1. Fit separate models for each administration wave.
#' 2. Combine into a named list: `fits <- list(Spring = fit_s, Fall = fit_f)`.
#' 3. Call `drift <- detect_anchor_drift(fits)`.
#' 4. Review `summary(drift)` and `plot_anchor_drift(drift)`.
#' 5. Flagged elements may need to be removed from anchor sets or
#'    investigated for substantive causes (e.g., rater re-training).
#'
#' @return Object of class `mfrm_anchor_drift` with components:
#'   \describe{
#'     \item{drift_table}{Tibble of element-level drift statistics.}
#'     \item{summary}{Drift summary aggregated by facet and wave.}
#'     \item{common_elements}{Tibble of pairwise common-element counts.}
#'     \item{common_by_facet}{Tibble of retained common-element counts by facet.}
#'     \item{config}{List of analysis configuration.}
#'   }
#'
#' @seealso [anchor_to_baseline()], [build_equating_chain()],
#'   [make_anchor_table()], [plot_anchor_drift()], [mfrmr_linking_and_dff]
#' @export
#' @examples
#' \donttest{
#' d1 <- load_mfrmr_data("study1")
#' d2 <- load_mfrmr_data("study2")
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 15)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 15)
#' drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
#' summary(drift)
#' head(drift$drift_table[, c("Facet", "Level", "Wave", "Drift", "Flag")])
#' drift$common_elements
#' }
detect_anchor_drift <- function(fits,
                                facets = NULL,
                                drift_threshold = 0.5,
                                flag_se_ratio = 2.0,
                                reference = 1L,
                                include_person = FALSE) {
  # Validate
  stopifnot(is.list(fits), length(fits) >= 2)
  for (f in fits) stopifnot(inherits(f, "mfrm_fit"))
  if (is.null(names(fits))) names(fits) <- paste0("Wave", seq_along(fits))

  ref_idx <- if (is.character(reference)) match(reference, names(fits)) else as.integer(reference)
  stopifnot(!is.na(ref_idx), ref_idx >= 1, ref_idx <= length(fits))

  # Extract estimates from each fit
  est_list <- lapply(fits, function(f) {
    make_anchor_table(f, facets = facets, include_person = include_person)
  })

  # Get SE from diagnostics$measures for each fit
  se_list <- lapply(fits, function(f) {
    measure_se_table(f, include_person = include_person)
  })

  ref_est <- est_list[[ref_idx]]
  wave_names <- names(fits)
  support_guideline <- 5L

  # Build drift table: for each non-reference wave, compute drift vs reference
  drift_rows <- list()
  support_rows <- list()
  for (i in seq_along(fits)) {
    if (i == ref_idx) next
    wave_est <- est_list[[i]]
    ref_se <- se_list[[ref_idx]]
    wave_se <- se_list[[i]]

    joined <- dplyr::inner_join(
      ref_est |> dplyr::rename(Ref_Est = "Anchor"),
      wave_est |> dplyr::rename(Wave_Est = "Anchor"),
      by = c("Facet", "Level")
    )

    if (nrow(joined) > 0) {
      joined <- joined |>
        dplyr::left_join(
          ref_se |> dplyr::rename(SE_Ref = "SE"),
          by = c("Facet", "Level")
        ) |>
        dplyr::left_join(
          wave_se |> dplyr::rename(SE_Wave = "SE"),
          by = c("Facet", "Level")
        )

      offset_info <- compute_equating_offset(
        diffs = joined$Wave_Est - joined$Ref_Est,
        se_from = joined$SE_Ref,
        se_to = joined$SE_Wave,
        drift_threshold = drift_threshold
      )
      support_tbl <- .summarise_link_support(
        joined[, c("Facet", "Level"), drop = FALSE],
        retained = offset_info$retained,
        guideline = support_guideline
      ) |>
        dplyr::mutate(
          Reference = wave_names[ref_idx],
          Wave = wave_names[i],
          .before = 1
        )
      support_rows <- c(support_rows, list(support_tbl))
      link_support_ok <- nrow(support_tbl) > 0 && all(support_tbl$LinkSupportAdequate)
      if (!link_support_ok) {
        weak_facets <- support_tbl$Facet[!support_tbl$LinkSupportAdequate]
        warning(
          sprintf(
            "Thin linking support between '%s' and '%s': fewer than %d retained common elements in %s.",
            wave_names[ref_idx],
            wave_names[i],
            support_guideline,
            paste(weak_facets, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      joined <- joined |>
        dplyr::mutate(
          Reference = wave_names[ref_idx],
          Wave = wave_names[i],
          LinkOffset = offset_info$offset,
          Drift = (.data$Wave_Est - .data$Ref_Est) - .data$LinkOffset,
          SE = ifelse(
            is.finite(.data$SE_Ref) & is.finite(.data$SE_Wave),
            sqrt(.data$SE_Ref^2 + .data$SE_Wave^2),
            NA_real_
          ),
          Drift_SE_Ratio = ifelse(
            is.na(.data$SE) | .data$SE == 0,
            NA_real_,
            abs(.data$Drift) / .data$SE
          ),
          LinkSupportAdequate = link_support_ok,
          Flag = abs(.data$Drift) > drift_threshold |
            (!is.na(.data$Drift_SE_Ratio) & .data$Drift_SE_Ratio > flag_se_ratio)
        )
      drift_rows <- c(drift_rows, list(joined))
    }
  }

  drift_table <- if (length(drift_rows) > 0) {
    dplyr::bind_rows(drift_rows) |>
      dplyr::select("Facet", "Level", "Reference", "Wave",
                     "Ref_Est", "Wave_Est", "LinkOffset", "Drift",
                     "SE_Ref", "SE_Wave", "SE", "Drift_SE_Ratio",
                     "LinkSupportAdequate", "Flag") |>
      dplyr::arrange(dplyr::desc(abs(.data$Drift)))
  } else {
    tibble::tibble(Facet = character(), Level = character(),
                   Reference = character(), Wave = character(),
                   Ref_Est = numeric(), Wave_Est = numeric(),
                   LinkOffset = numeric(), Drift = numeric(),
                   SE_Ref = numeric(), SE_Wave = numeric(), SE = numeric(),
                   Drift_SE_Ratio = numeric(), LinkSupportAdequate = logical(),
                   Flag = logical())
  }

  # Summary by facet
  drift_summary <- if (nrow(drift_table) > 0) {
    drift_table |>
      dplyr::group_by(.data$Facet, .data$Wave) |>
      dplyr::summarise(
        N = dplyr::n(), Mean_Drift = mean(abs(.data$Drift)),
        Max_Drift = max(abs(.data$Drift)),
        N_Flagged = sum(.data$Flag, na.rm = TRUE), .groups = "drop"
      )
  } else {
    tibble::tibble()
  }

  # Common elements count
  common_counts <- tibble::tibble(
    Wave1 = character(), Wave2 = character(), N_Common = integer()
  )
  for (i in seq_along(fits)) {
    for (j in seq_along(fits)) {
      if (j <= i) next
      n_common <- nrow(dplyr::inner_join(est_list[[i]], est_list[[j]],
                                          by = c("Facet", "Level")))
      common_counts <- dplyr::bind_rows(common_counts,
        tibble::tibble(Wave1 = wave_names[i], Wave2 = wave_names[j],
                       N_Common = as.integer(n_common)))
    }
  }

  out <- list(
    drift_table = drift_table, summary = drift_summary,
    common_elements = common_counts,
    common_by_facet = if (length(support_rows) > 0) dplyr::bind_rows(support_rows) else tibble::tibble(),
    config = list(reference = wave_names[ref_idx],
                  method = "screened_common_element_alignment",
                  intended_use = "review_screen",
                  drift_threshold = drift_threshold,
                  min_common_per_facet = support_guideline,
                  flag_se_ratio = flag_se_ratio,
                  facets = facets, waves = wave_names)
  )
  class(out) <- c("mfrm_anchor_drift", "list")
  out
}

#' @rdname detect_anchor_drift
#' @param x An `mfrm_anchor_drift` object.
#' @param ... Ignored.
#' @export
print.mfrm_anchor_drift <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' @rdname detect_anchor_drift
#' @param object An `mfrm_anchor_drift` object (for `summary`).
#' @export
summary.mfrm_anchor_drift <- function(object, ...) {
  dt <- object$drift_table
  out <- list(
    n_comparisons = nrow(dt), n_flagged = sum(dt$Flag, na.rm = TRUE),
    summary = object$summary, common_elements = object$common_elements,
    common_by_facet = object$common_by_facet,
    flagged = dt |> dplyr::filter(.data$Flag),
    config = object$config
  )
  class(out) <- "summary.mfrm_anchor_drift"
  out
}

#' @rdname detect_anchor_drift
#' @export
print.summary.mfrm_anchor_drift <- function(x, ...) {
  cat("--- Anchor Drift Screen ---\n")
  cat("Reference:", x$config$reference, "\n")
  cat("Method:", x$config$method, "| Intended use:", x$config$intended_use, "\n")
  cat("Comparisons:", x$n_comparisons, "| Flagged:", x$n_flagged, "\n\n")
  if (nrow(x$summary) > 0) {
    cat("Drift summary by facet and wave:\n")
    print(as.data.frame(x$summary), row.names = FALSE, digits = 3)
  }
  if (nrow(x$common_elements) > 0) {
    cat("\nCommon elements:\n")
    print(as.data.frame(x$common_elements), row.names = FALSE)
  }
  if (nrow(x$common_by_facet) > 0) {
    cat("\nRetained common elements by facet:\n")
    print(as.data.frame(x$common_by_facet), row.names = FALSE)
  }
  if (nrow(x$flagged) > 0) {
    cat("\nFlagged elements:\n")
    print(as.data.frame(x$flagged |> utils::head(20)), row.names = FALSE, digits = 3)
    if (nrow(x$flagged) > 20) cat("... (", nrow(x$flagged) - 20, " more)\n")
  }
  invisible(x)
}

# --- build_equating_chain ----------------------------------------------------

#' Build a screened linking chain across ordered calibrations
#'
#' Links a series of calibration waves by computing mean offsets between
#' adjacent pairs of fits. Common linking elements (e.g., raters or items
#' that appear in consecutive administrations) are used to estimate the
#' scale shift. Cumulative offsets place all waves on a common metric
#' anchored to the first wave. The procedure is intended as a practical
#' screened linking aid, not as a full general-purpose equating framework.
#'
#' @param fits Named list of `mfrm_fit` objects in chain order.
#' @param anchor_facets Character vector of facets to use as linking
#'   elements.
#' @param include_person Include person estimates in linking.
#' @param drift_threshold Threshold for flagging large residuals in links.
#'
#' @details
#' The screened linking chain uses a screened link-offset method.  For each pair of
#' adjacent waves \eqn{(A, B)}, the function:
#'
#' 1. Identifies common linking elements (facet levels present in both fits).
#' 2. Computes per-element differences:
#'    \deqn{d_e = \hat{\delta}_{e,B} - \hat{\delta}_{e,A}}
#' 3. Computes a preliminary link offset using the inverse-variance weighted
#'    mean of these differences when standard errors are available (otherwise
#'    an unweighted mean).
#' 4. Screens out elements whose residual from that preliminary offset exceeds
#'    `drift_threshold`, then recomputes the final offset on the retained set.
#' 5. Records `Offset_SD` (standard deviation of retained residuals) and
#'    `Max_Residual` (maximum absolute deviation from the mean) as
#'    indicators of link quality.
#' 6. Flags links with fewer than 5 retained common elements in any linking
#'    facet as having thin support.
#'
#' Cumulative offsets are computed by chaining link offsets from Wave 1
#' forward, placing all waves onto the metric of the first wave.
#'
#' Elements whose per-link residual exceeds `drift_threshold` are flagged
#' in `$element_detail$Flag`.  A high `Offset_SD`, many flagged elements, or a
#' thin retained anchor set signals an unstable link that may compromise the
#' resulting scale placement.
#'
#' @section Which function should I use?:
#' - Use [anchor_to_baseline()] for a single new wave anchored to a known
#'   baseline.
#' - Use [detect_anchor_drift()] when you want direct comparison against one
#'   reference wave.
#' - Use `build_equating_chain()` when no single wave should dominate and you
#'   want ordered, adjacent links across the series.
#'
#' @section Interpreting output:
#' - `$links`: one row per adjacent pair with `From`, `To`, `N_Common`,
#'   `N_Retained`, `Offset_Prelim`, `Offset`, `Offset_SD`, and
#'   `Max_Residual`. Small `Offset_SD`
#'   relative to the offset indicates a consistent shift across elements.
#'   `LinkSupportAdequate = FALSE` means at least one linking facet retained
#'   fewer than 5 common elements after screening.
#' - `$cumulative`: one row per wave with its cumulative offset from Wave 1.
#'   Wave 1 always has offset 0.
#' - `$element_detail`: per-element linking statistics (estimate in each
#'   wave, difference, residual from mean offset, and flag status).
#'   Flagged elements may indicate DIF or rater re-training effects.
#' - `$common_by_facet`: retained common-element counts by linking facet for
#'   each adjacent link.
#' - `$config`: records wave names and analysis parameters.
#' - Read `links` before `cumulative`: weak adjacent links can make later
#'   cumulative offsets less trustworthy.
#'
#' @section Typical workflow:
#' 1. Fit each administration wave separately: `fit_a <- fit_mfrm(...)`.
#' 2. Combine into an ordered named list:
#'    `fits <- list(Spring23 = fit_s, Fall23 = fit_f, Spring24 = fit_s2)`.
#' 3. Call `chain <- build_equating_chain(fits)`.
#' 4. Review `summary(chain)` for link quality.
#' 5. Visualize with `plot_anchor_drift(chain, type = "chain")`.
#' 6. For problematic links, investigate flagged elements in
#'    `chain$element_detail` and consider removing them from the anchor set.
#'
#' @return Object of class `mfrm_equating_chain` with components:
#'   \describe{
#'     \item{links}{Tibble of link-level statistics (offset, SD, etc.).}
#'     \item{cumulative}{Tibble of cumulative offsets per wave.}
#'     \item{element_detail}{Tibble of element-level linking details.}
#'     \item{common_by_facet}{Tibble of retained common-element counts by facet.}
#'     \item{config}{List of analysis configuration.}
#'   }
#'
#' @seealso [detect_anchor_drift()], [anchor_to_baseline()],
#'   [make_anchor_table()], [plot_anchor_drift()]
#' @export
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' people <- unique(toy$Person)
#' d1 <- toy[toy$Person %in% people[1:12], , drop = FALSE]
#' d2 <- toy[toy$Person %in% people[13:24], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 10)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 10)
#' chain <- build_equating_chain(list(Form1 = fit1, Form2 = fit2))
#' summary(chain)
#' chain$cumulative
#' }
build_equating_chain <- function(fits,
                                 anchor_facets = NULL,
                                 include_person = FALSE,
                                 drift_threshold = 0.5) {
  stopifnot(is.list(fits), length(fits) >= 2)
  for (f in fits) stopifnot(inherits(f, "mfrm_fit"))
  if (is.null(names(fits))) names(fits) <- paste0("Form", seq_along(fits))

  wave_names <- names(fits)
  n_waves <- length(fits)

  # Extract estimates
  est_list <- lapply(fits, function(f) {
    make_anchor_table(f, facets = anchor_facets, include_person = include_person)
  })
  se_list <- lapply(fits, function(f) {
    measure_se_table(f, include_person = include_person)
  })

  # Build links between adjacent pairs
  links <- list()
  element_details <- list()
  support_rows <- list()
  support_guideline <- 5L

  for (i in seq_len(n_waves - 1)) {
    from <- est_list[[i]]
    to <- est_list[[i + 1]]

    common <- dplyr::inner_join(
      from |> dplyr::rename(Est_From = "Anchor"),
      to |> dplyr::rename(Est_To = "Anchor"),
      by = c("Facet", "Level")
    )

    n_common <- nrow(common)

    if (n_common == 0) {
      warning(sprintf("No common elements between '%s' and '%s'.",
                       wave_names[i], wave_names[i + 1]),
              call. = FALSE)
      offset <- NA_real_
      offset_sd <- NA_real_
      max_drift <- NA_real_
    } else {
      common <- common |>
        dplyr::left_join(
          se_list[[i]] |> dplyr::rename(SE_From = "SE"),
          by = c("Facet", "Level")
        ) |>
        dplyr::left_join(
          se_list[[i + 1]] |> dplyr::rename(SE_To = "SE"),
          by = c("Facet", "Level")
        )
      diffs <- common$Est_To - common$Est_From
      offset_info <- compute_equating_offset(
        diffs = diffs,
        se_from = common$SE_From,
        se_to = common$SE_To,
        drift_threshold = drift_threshold
      )
      support_tbl <- .summarise_link_support(
        common[, c("Facet", "Level"), drop = FALSE],
        retained = offset_info$retained,
        guideline = support_guideline
      ) |>
        dplyr::mutate(
          Link = i,
          From = wave_names[i],
          To = wave_names[i + 1],
          .before = 1
        )
      support_rows <- c(support_rows, list(support_tbl))
      link_support_ok <- nrow(support_tbl) > 0 && all(support_tbl$LinkSupportAdequate)
      if (!link_support_ok) {
        weak_facets <- support_tbl$Facet[!support_tbl$LinkSupportAdequate]
        warning(
          sprintf(
            "Thin linking support between '%s' and '%s': fewer than %d retained common elements in %s.",
            wave_names[i],
            wave_names[i + 1],
            support_guideline,
            paste(weak_facets, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      offset <- offset_info$offset
      offset_sd <- if (sum(offset_info$retained, na.rm = TRUE) > 1) {
        stats::sd(offset_info$residual[offset_info$retained], na.rm = TRUE)
      } else {
        0
      }
      max_drift <- max(abs(offset_info$residual), na.rm = TRUE)

      common <- common |>
        dplyr::mutate(
          Link = paste0(wave_names[i], " -> ", wave_names[i + 1]),
          Diff = .data$Est_To - .data$Est_From,
          Offset_Prelim = offset_info$offset_prelim,
          Offset = offset,
          Residual = offset_info$residual,
          Retained = offset_info$retained,
          Flag = abs(.data$Residual) > drift_threshold
        )
      element_details <- c(element_details, list(common))
    }

    links <- c(links, list(tibble::tibble(
      Link = i,
      From = wave_names[i], To = wave_names[i + 1],
      N_Common = as.integer(n_common),
      N_Retained = if (n_common > 0) offset_info$n_retained else 0L,
      Min_Common_Per_Facet = if (n_common > 0 && nrow(support_tbl) > 0) min(support_tbl$N_Common) else 0L,
      Min_Retained_Per_Facet = if (n_common > 0 && nrow(support_tbl) > 0) min(support_tbl$N_Retained) else 0L,
      Offset_Prelim = if (n_common > 0) offset_info$offset_prelim else NA_real_,
      Offset = offset,
      Offset_SD = offset_sd,
      Max_Residual = max_drift,
      LinkSupportAdequate = if (n_common > 0) link_support_ok else FALSE,
      Offset_Method = if (n_common > 0) offset_info$weighting else "none"
    )))
  }

  links_tbl <- dplyr::bind_rows(links)

  # Cumulative offsets
  cum_offset <- cumsum(c(0, links_tbl$Offset))
  cumulative <- tibble::tibble(
    Wave = wave_names,
    Cumulative_Offset = cum_offset
  )

  element_detail <- if (length(element_details) > 0) {
    dplyr::bind_rows(element_details)
  } else {
    tibble::tibble()
  }

  out <- list(
    links = links_tbl, cumulative = cumulative,
    element_detail = element_detail,
    common_by_facet = if (length(support_rows) > 0) dplyr::bind_rows(support_rows) else tibble::tibble(),
    config = list(anchor_facets = anchor_facets,
                  method = "screened_common_element_alignment",
                  intended_use = "screened_linking_aid",
                  min_common_per_facet = support_guideline,
                  drift_threshold = drift_threshold,
                  waves = wave_names)
  )
  class(out) <- c("mfrm_equating_chain", "list")
  out
}

#' @rdname build_equating_chain
#' @param x An `mfrm_equating_chain` object.
#' @param ... Ignored.
#' @export
print.mfrm_equating_chain <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' @rdname build_equating_chain
#' @param object An `mfrm_equating_chain` object (for `summary`).
#' @export
summary.mfrm_equating_chain <- function(object, ...) {
  out <- list(links = object$links, cumulative = object$cumulative,
              common_by_facet = object$common_by_facet,
              n_flagged = sum(object$element_detail$Flag, na.rm = TRUE),
              config = object$config)
  class(out) <- "summary.mfrm_equating_chain"
  out
}

#' @rdname build_equating_chain
#' @export
print.summary.mfrm_equating_chain <- function(x, ...) {
  cat("--- Screened Linking Chain ---\n")
  cat("Method:", x$config$method, "| Intended use:", x$config$intended_use, "\n")
  cat("Links:", nrow(x$links), "| Waves:",
      paste(x$config$waves, collapse = " -> "), "\n\n")
  cat("Link details:\n")
  print(as.data.frame(x$links), row.names = FALSE, digits = 3)
  if (nrow(x$common_by_facet) > 0) {
    cat("\nRetained common elements by facet:\n")
    print(as.data.frame(x$common_by_facet), row.names = FALSE, digits = 3)
  }
  cat("\nCumulative offsets:\n")
  print(as.data.frame(x$cumulative), row.names = FALSE, digits = 3)
  if (x$n_flagged > 0) cat("\nFlagged linking elements:", x$n_flagged, "\n")
  invisible(x)
}
