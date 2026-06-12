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
  method <- resolve_public_mfrm_method(
    summary_method = sub_fit$summary$Method[1] %||% NA_character_,
    method_input = sub_fit$config$method_input %||% NA_character_,
    method_used = sub_fit$config$method %||% NA_character_
  )
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
  anchor_review_obj <- anchor_review(fit, required = FALSE) %||% list()
  min_common_anchors <- anchor_review_obj$thresholds$min_common_anchors %||% 5L
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

  review <- anchor_review(sub_fit, required = FALSE) %||% list()
  facet_summary <- as.data.frame(review$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(facet_summary) == 0 || !"Facet" %in% names(facet_summary)) {
    return(list(
      status = "weak_link",
      ets_eligible = FALSE,
      anchored_levels = NA_integer_,
      detail = "Linking anchors were requested, but subgroup anchor coverage could not be reviewed."
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
                                        linking_setup, linkage, diagnostics_error = NULL,
                                        linking_review_note = NA_character_) {
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
      LinkingReview = linking_review_note,
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
#'   `NULL` (default), mfrmr tries to recover the data from
#'   `fit$prep$data`. That slot only holds the columns that
#'   `fit_mfrm()` actually modelled; if the grouping column was not
#'   among them (common for DIF screening), pass the original data
#'   frame via `data = <df>` explicitly. The same applies when the
#'   fit object has been serialized without the prep slot.
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
#' - `$gpcm_boundary`: for bounded `GPCM` fits, a capability-boundary table
#'   marking the DFF/DIF output as caveated screening evidence.
#'
#' @section GPCM boundary:
#' For bounded `GPCM`, DFF/DIF rows are available as slope-aware screening
#' evidence over the fitted expected-score and residual scale. Keep
#' residual-method contrasts and interaction cells in screening language.
#' Refit contrasts require explicit subgroup linking and precision support
#' before stronger subgroup-comparison language is used.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()]. For `RSM` / `PCM` fairness review, prefer
#'    `method = "MML"`.
#' 2. Run [diagnose_mfrm()] and, for `RSM` / `PCM`, prefer
#'    `diagnostic_mode = "both"` so legacy and strict marginal screens remain
#'    visible together.
#' 3. Run `analyze_dff(fit, diagnostics, facet = "Criterion", group = "Gender", data = my_data)`.
#' 4. Inspect `$dif_table` for flagged levels and `$summary` for counts.
#' 5. Use [dif_interaction_table()] when you need cell-level diagnostics.
#' 6. Use [plot_dif_heatmap()] or [dif_report()] for communication.
#'
#' @return
#' An object of class `mfrm_dff` (with compatibility class `mfrm_dif`) with:
#' - `dif_table`: data.frame of differential-functioning contrasts.
#' - `cell_table`: (residual method) per-cell detail table.
#' - `summary`: counts by screening or ETS classification.
#' - `group_fits`: (refit method) per-group facet estimates.
#' - `gpcm_boundary`: for bounded `GPCM` fits, a capability-boundary table.
#' - `config`: list with facet, group, method, min_obs, p_adjust settings.
#'
#' @seealso [fit_mfrm()], [estimate_bias()], [compare_mfrm()],
#'   [dif_interaction_table()], [plot_dif_heatmap()], [dif_report()],
#'   [subset_connectivity_report()], [mfrmr_linking_and_dff]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", model = "RSM", quad_points = 7, maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#' dff <- analyze_dff(fit, diag, facet = "Rater", group = "Group", data = toy)
#' dff$summary
#' # Look for: a small `FlaggedPairs` count relative to `Pairs`. Under
#' #   method = "residual", `ClassificationSystem` is "screening", not
#' #   ETS. "Screen positive" rows are prompts for substantive review.
#' head(dff$dif_table[, c("Level", "Group1", "Group2", "Contrast",
#'                        "Classification", "ClassificationSystem")])
#' # The residual contrast is an observed-minus-expected average contrast
#' # between groups. It is useful for screening, but it is not an ETS
#' # A/B/C logit-delta classification.
#' dff_refit <- analyze_dff(fit, diag, facet = "Rater", group = "Group",
#'                          data = toy, method = "refit")
#' unique(dff_refit$dif_table$ClassificationSystem)
#' # Look for: "ETS" only when subgroup calibration, linking, and precision
#' #   checks all support a common-scale model-based contrast.
#' sc <- subset_connectivity_report(fit, diagnostics = diag)
#' plot(sc, type = "design_matrix", draw = FALSE)
#' if ("ScaleLinkStatus" %in% names(dff_refit$dif_table)) {
#'   unique(dff_refit$dif_table$ScaleLinkStatus)
#' }
#' # Look for: "linked" in `ScaleLinkStatus` confirms the focal and
#' #   reference groups share enough common elements for a comparable
#' #   contrast; "demoted_*" rows lose linking under the refit branch
#' #   and should be read as exploratory.
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
  # `method = "refit"` needs diagnostics for the anchored-fit SE path;
  # `method = "residual"` does not. Explicitly guard the refit branch
  # so the error is mfrmr-styled and locale-independent rather than
  # R's default "argument ... is missing".
  if (identical(method, "refit") && missing(diagnostics)) {
    stop("`diagnostics` is required for analyze_dff(method = \"refit\"). ",
         "Call `diagnose_mfrm(fit)` first and pass the result as the ",
         "second argument, or use `method = \"residual\"` for the ",
         "diagnostics-free screening path.",
         call. = FALSE)
  }
  if (!is.character(facet) || length(facet) != 1) {
    stop("`facet` must be a single character string.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1) {
    stop("`group` must be a single character string naming a column in the ",
         "original data.", call. = FALSE)
  }
  min_obs <- .validate_dff_count_arg(min_obs, "min_obs")
  p_adjust <- .validate_p_adjust_method(p_adjust)

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
  group_info <- .sanitize_dff_group_data(orig_data, group, "DFF analysis")
  orig_data <- group_info$data
  group_vals <- group_info$values
  group_levels <- group_info$levels
  if (length(group_levels) < 2) {
    stop("Grouping variable '", group, "' must have at least 2 levels. ",
         "Found ", length(group_levels), " after removing missing or empty ",
         "group values.", call. = FALSE)
  }
  focal <- .validate_dff_focal(focal, group_levels)

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
      group_vals = group_vals,
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

.validate_dff_count_arg <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 1 ||
      abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    stop("`", arg, "` must be a single positive integer.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_p_adjust_method <- function(p_adjust) {
  if (length(p_adjust) != 1L || is.na(p_adjust)) {
    stop("`p_adjust` must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".",
         call. = FALSE)
  }
  p_adjust <- as.character(p_adjust)
  if (!p_adjust %in% stats::p.adjust.methods) {
    stop("`p_adjust` must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".",
         call. = FALSE)
  }
  p_adjust
}

.sanitize_dff_group_data <- function(orig_data, group, context) {
  group_vals <- trimws(as.character(orig_data[[group]]))
  valid_group <- !is.na(group_vals) & nzchar(group_vals)
  if (!any(valid_group)) {
    stop("`group` column '", group,
         "' has no non-missing, non-empty values after trimming.",
         call. = FALSE)
  }
  if (any(!valid_group)) {
    message("Dropped ", sum(!valid_group), " row(s) with missing or empty `",
            group, "` values before ", context, ".")
    orig_data <- orig_data[valid_group, , drop = FALSE]
    group_vals <- group_vals[valid_group]
  }
  orig_data[[group]] <- group_vals
  list(
    data = orig_data,
    values = group_vals,
    levels = sort(unique(group_vals))
  )
}

.validate_dff_focal <- function(focal, group_levels) {
  if (is.null(focal)) return(NULL)
  if (!is.character(focal) || length(focal) < 1L || anyNA(focal)) {
    stop("`focal` must be a character vector of observed group levels.",
         call. = FALSE)
  }
  focal <- unique(trimws(focal))
  if (any(!nzchar(focal))) {
    stop("`focal` must be a character vector of observed group levels.",
         call. = FALSE)
  }
  unknown <- setdiff(focal, group_levels)
  if (length(unknown) > 0L) {
    stop("`focal` contains level(s) not found in `group`: ",
         paste(unknown, collapse = ", "), ". Observed group levels: ",
         paste(group_levels, collapse = ", "), ".",
         call. = FALSE)
  }
  if (length(setdiff(group_levels, focal)) == 0L) {
    stop("`focal` cannot include every observed group level; at least one ",
         "reference group is required.", call. = FALSE)
  }
  focal
}

.validate_dff_probability <- function(x, arg, allow_null = TRUE) {
  if (is.null(x) && isTRUE(allow_null)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must be a single probability strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(x)
}

.validate_dff_nonnegative_scalar <- function(x, arg, allow_null = TRUE) {
  if (is.null(x) && isTRUE(allow_null)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0) {
    stop("`", arg, "` must be a single non-negative finite number.",
         call. = FALSE)
  }
  as.numeric(x)
}

.validate_dff_positive_scalar <- function(x, arg, allow_null = TRUE) {
  if (is.null(x) && isTRUE(allow_null)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive finite number.",
         call. = FALSE)
  }
  as.numeric(x)
}

.validate_dff_nonnegative_integer <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 ||
      abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    stop("`", arg, "` must be a single non-negative integer.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_dff_threshold_vector <- function(x, arg) {
  if (is.null(x)) return(numeric(0))
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) ||
      any(!is.finite(x)) || any(x < 0)) {
    stop("`", arg, "` must be a numeric vector of non-negative finite values.",
         call. = FALSE)
  }
  vals <- as.numeric(x)
  labs <- names(x)
  if (is.null(labs)) labs <- rep("", length(vals))
  ord <- order(vals)
  vals <- vals[ord]
  labs <- labs[ord]
  keep <- !duplicated(vals)
  vals <- vals[keep]
  labs <- labs[keep]
  names(vals) <- labs
  vals
}

.validate_dff_logical_scalar <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(x)
}

.dff_effect_axis_label <- function(method) {
  method <- method %||% NA_character_
  if (identical(method, "residual")) {
    "Observed-minus-expected screening contrast"
  } else if (identical(method, "refit")) {
    "Differential-functioning contrast (logit)"
  } else {
    "Effect (native contrast scale)"
  }
}

.dff_interpretation_guide <- function(metric,
                                      method = NA_character_,
                                      classification_system = NA_character_,
                                      flag_threshold = NULL,
                                      effect_thresholds = numeric()) {
  rows <- list(
    data.frame(
      Item = "Zero reference",
      Meaning = "Values near zero indicate little systematic group-by-facet departure on the selected scale.",
      ReportingNote = "Use the sign convention and contrast basis reported in the source table.",
      stringsAsFactors = FALSE
    )
  )
  if (identical(metric, "obs_exp")) {
    rows[[length(rows) + 1L]] <- data.frame(
      Item = "Warm/cool colors",
      Meaning = "Positive values mean observed scores exceed model expectation; negative values mean observed scores fall below model expectation.",
      ReportingNote = "Read as residual screening evidence, not as a standalone inferential test.",
      stringsAsFactors = FALSE
    )
  } else if (identical(metric, "t")) {
    rows[[length(rows) + 1L]] <- data.frame(
      Item = "Standardized residual",
      Meaning = "Larger absolute values identify cells with stronger observed-minus-expected departure relative to model variance.",
      ReportingNote = "Use as a screening statistic; combine with sample size and substantive review.",
      stringsAsFactors = FALSE
    )
  } else if (identical(metric, "contrast") || identical(metric, "summary")) {
    rows[[length(rows) + 1L]] <- data.frame(
      Item = "Differential-functioning contrast",
      Meaning = if (identical(method, "residual")) {
        "Residual-method effects are observed-minus-expected average contrasts between groups."
      } else {
        "Refit-method effects are subgroup parameter differences when linking supports a comparable scale."
      },
      ReportingNote = if (identical(classification_system, "ETS")) {
        "ETS A/B/C labels may be reported for rows that remain model-based and linked."
      } else {
        "Treat classifications as screening labels unless the row reports ClassificationSystem == 'ETS'."
      },
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(flag_threshold)) {
    rows[[length(rows) + 1L]] <- data.frame(
      Item = "Flag threshold",
      Meaning = paste0("Cells with absolute value >= ", flag_threshold,
                       " are marked in the plot data and outlined when drawn."),
      ReportingNote = "The threshold is user-specified plotting guidance unless it comes from a documented analysis rule.",
      stringsAsFactors = FALSE
    )
  }
  if (length(effect_thresholds) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      Item = "Effect thresholds",
      Meaning = paste(effect_thresholds, collapse = ", "),
      ReportingNote = "Threshold guide lines are display aids; check the method and classification system before interpreting them as severity bands.",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
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
  method_label <- resolve_public_mfrm_method(
    summary_method = fit$summary$Method[1] %||% NA_character_,
    method_input = fit$config$method_input %||% NA_character_,
    method_used = fit$config$method %||% NA_character_
  )
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
      # ContrastDirection resolves the historical sign flip between
      # residual and refit methods for DFF screening. Under the
      # residual method `contrast = avg1 - avg2` is the mean observed-
      # minus-expected difference: contrast > 0 means Group1 scored
      # higher than expected, i.e. the facet was *easier* for Group1.
      direction <- if (is.finite(contrast)) {
        if (contrast > 0) "easier_for_group1" else "harder_for_group1"
      } else {
        NA_character_
      }
      dif_rows[[length(dif_rows) + 1]] <- tibble(
        Level = lev,
        Group1 = g1,
        Group2 = g2,
        Contrast = contrast,
        ContrastDirection = direction,
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "analyze_dff()",
      area = "Differential facet functioning screening under bounded GPCM"
    ),
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
    # Capture the anchor-review issue messages emitted while refitting
    # the subgroup so a silent anchor_policy no longer hides
    # unknown-level / invalid-anchor issues. The outer tryCatch /
    # suppressWarnings is retained so a fit failure is still handled
    # gracefully.
    linking_review_msgs <- character(0)
    sub_fit <- tryCatch(
      withCallingHandlers(
        do.call(fit_mfrm, fit_args),
        message = function(m) {
          msg <- conditionMessage(m)
          if (grepl("Anchor review", msg, fixed = TRUE) ||
              grepl("anchor_policy", msg, fixed = TRUE)) {
            linking_review_msgs <<- c(linking_review_msgs, msg)
          }
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) structure(list(message = conditionMessage(e)), class = "mfrm_dff_fit_error")
    )
    linking_review_text <- if (length(linking_review_msgs) > 0L) {
      paste(unique(linking_review_msgs), collapse = " | ")
    } else {
      NA_character_
    }
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
          LinkingReview = linking_review_text,
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
        diagnostics_error = sub_diag_error,
        linking_review_note = linking_review_text
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
      # Under the refit method `contrast = e1 - e2` is the difference
      # in difficulty parameters: contrast > 0 means Group1's facet
      # estimate is *larger* than Group2's, i.e. the facet was *harder*
      # for Group1. This is the opposite sign convention to the residual
      # method; the ContrastDirection column is provided for both so
      # downstream consumers never need to reason about the sign.
      direction <- if (is.finite(contrast)) {
        if (contrast > 0) "harder_for_group1" else "easier_for_group1"
      } else {
        NA_character_
      }
      dif_rows[[length(dif_rows) + 1]] <- tibble(
        Level = merged$Level[j],
        Group1 = g1,
        Group2 = g2,
        Estimate1 = e1,
        Estimate2 = e2,
        Contrast = contrast,
        ContrastDirection = direction,
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "analyze_dff()",
      area = "Differential facet functioning screening under bounded GPCM"
    ),
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
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
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
  .print_dff_gpcm_boundary(x$gpcm_boundary)
  invisible(x)
}

#' @export
print.summary.mfrm_dff <- function(x, ...) {
  print.summary.mfrm_dif(x, ...)
}

#' @export
print.mfrm_dif <- function(x, ...) {
  .print_dif_overview(x, label = "DIF")
  invisible(x)
}

#' @export
print.mfrm_dff <- function(x, ...) {
  .print_dif_overview(x, label = "DFF")
  invisible(x)
}

# Short scalar-ish overview for print() that complements the full table
# produced by summary(). Shows method / facet / group / flagged-count
# so the user can tell at a glance whether a detailed inspection is
# warranted; callers who want the full contrast table should call
# summary() explicitly.
.print_dif_overview <- function(x, label) {
  cfg <- x$config %||% list()
  method <- cfg$method %||% "refit"
  facet <- cfg$facet %||% NA_character_
  group <- cfg$group %||% NA_character_
  groups <- cfg$group_levels %||% character(0)
  tbl <- x$dif_table
  n_rows <- if (is.data.frame(tbl)) nrow(tbl) else 0L
  n_flag <- if (n_rows > 0L && "Classification" %in% names(tbl)) {
    sum(!is.na(tbl$Classification) & tbl$Classification != "Negligible" &
          tbl$Classification != "None")
  } else NA_integer_
  cat(sprintf("mfrm_%s (%s)\n", tolower(label), label))
  cat(sprintf("  Method: %s | Facet: %s | Group: %s\n",
              method, facet, group))
  if (length(groups) > 0L) {
    cat(sprintf("  Group levels: %s\n",
                paste(groups, collapse = ", ")))
  }
  if (n_rows > 0L) {
    cat(sprintf("  Contrasts: %d row(s)", n_rows))
    if (is.finite(n_flag)) cat(sprintf(", %d flagged", n_flag))
    cat("\n")
  } else {
    cat("  Contrasts: 0 row(s)\n")
  }
  boundary <- x$gpcm_boundary %||% data.frame()
  if (is.data.frame(boundary) && nrow(boundary) > 0L) {
    status <- as.character(boundary$Status[1] %||% "supported_with_caveat")
    cat(sprintf("  GPCM boundary: %s\n", status))
  }
  cat("  Use `summary()` for the contrast table and classification breakdown.\n")
  invisible(NULL)
}

.print_dff_gpcm_boundary <- function(boundary) {
  if (is.null(boundary) || !is.data.frame(boundary) || nrow(boundary) == 0L) {
    return(invisible(NULL))
  }
  keep <- intersect(c("Area", "Status"), names(boundary))
  if (length(keep) == 0L) {
    return(invisible(NULL))
  }
  cat("\nGPCM Boundary:\n")
  print(as.data.frame(boundary[, keep, drop = FALSE]), row.names = FALSE)
  invisible(NULL)
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
#' - `$gpcm_boundary`: for bounded `GPCM` fits, a capability-boundary table
#'   marking the table as caveated DFF screening evidence.
#' - Cells with `|t| > abs_t_warn` or `|ObsExpAvg| > abs_bias_warn`
#'   are flagged in the `flag_t` and `flag_bias` columns.
#' - Sparse cells (N < min_obs) have `sparse = TRUE` and NA statistics.
#'
#' @section GPCM boundary:
#' For bounded `GPCM`, the interaction table uses the fitted slope-aware
#' expected-score/residual scale and should be reported as screening evidence,
#' not as a standalone fairness, invariance, or operational subgroup decision.
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
#' - `gpcm_boundary`: for bounded `GPCM` fits, a capability-boundary table.
#' - `config`: list of analysis parameters.
#'
#' @seealso [analyze_dff()], [analyze_dif()], [plot_dif_heatmap()], [dif_report()],
#'   [estimate_bias()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 30)
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
  min_obs <- .validate_dff_count_arg(min_obs, "min_obs")
  p_adjust <- .validate_p_adjust_method(p_adjust)
  if (!is.numeric(abs_t_warn) || length(abs_t_warn) != 1 ||
      is.na(abs_t_warn) || !is.finite(abs_t_warn) || abs_t_warn < 0) {
    stop("`abs_t_warn` must be a single non-negative numeric value.",
         call. = FALSE)
  }
  if (!is.numeric(abs_bias_warn) || length(abs_bias_warn) != 1 ||
      is.na(abs_bias_warn) || !is.finite(abs_bias_warn) ||
      abs_bias_warn < 0) {
    stop("`abs_bias_warn` must be a single non-negative numeric value.",
         call. = FALSE)
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

  group_info <- .sanitize_dff_group_data(orig_data, group,
                                         "DFF interaction analysis")
  orig_data <- group_info$data
  group_levels <- group_info$levels
  if (length(group_levels) < 2) {
    stop("Grouping variable '", group, "' must have at least 2 levels. ",
         "Found ", length(group_levels), " after removing missing or empty ",
         "group values.", call. = FALSE)
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "dif_interaction_table()",
      area = "Differential facet functioning screening under bounded GPCM"
    ),
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
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
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
  .print_dff_gpcm_boundary(x$gpcm_boundary)
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
#' @param show_values Logical. If `TRUE` (default), print rounded cell values
#'   on top of the heatmap.
#' @param value_digits Non-negative integer number of digits after the decimal
#'   point for cell labels.
#' @param flag_threshold Optional non-negative absolute-value threshold. When
#'   supplied, cells with `abs(value) >= flag_threshold` are recorded in
#'   `$data$flag_matrix` and outlined on the drawn heatmap.
#' @param scale_limit Optional positive scalar for a symmetric color scale
#'   from `-scale_limit` to `+scale_limit`. Use this to make several heatmaps
#'   visually comparable.
#' @param flag_color Border color for cells meeting `flag_threshold`.
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
#' @return Invisibly, an `mfrm_plot_data` object whose `data` slot bundles
#'   the row x column metric matrix (`$matrix`), the source long table
#'   (`$pairs`), and the metric label. Earlier 0.1.x releases returned the
#'   bare matrix; consume `$data$matrix` to keep code forward-compatible.
#'
#' @seealso [dif_interaction_table()], [analyze_dff()], [analyze_dif()], [dif_report()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' int <- dif_interaction_table(fit, diag, facet = "Rater",
#'                              group = "Group", data = toy, min_obs = 2)
#' heat <- plot_dif_heatmap(int, metric = "obs_exp", draw = FALSE)
#' dim(heat$data$matrix)
#' # Look for (`metric = "obs_exp"`): cells near 0 are aligned with
#' #   model expectation; |Obs - Exp| > 0.5 logits is a substantive
#' #   gap. With `metric = "t"` the cell scale becomes a standardized
#' #   residual where |t| > 2 is a screening flag, not a standalone
#' #   hypothesis test. With `metric = "contrast"` the layout switches
#' #   to Level x GroupPair and reads as the pairwise differential-
#' #   functioning contrast (use `analyze_dff()`).
#' @export
plot_dif_heatmap <- function(x, metric = c("obs_exp", "t", "contrast"),
                             draw = TRUE,
                             show_values = TRUE,
                             value_digits = 2L,
                             flag_threshold = NULL,
                             scale_limit = NULL,
                             flag_color = "black",
                             ...) {
  metric <- match.arg(metric)
  show_values <- .validate_dff_logical_scalar(show_values, "show_values")
  value_digits <- .validate_dff_nonnegative_integer(value_digits, "value_digits")
  flag_threshold <- .validate_dff_nonnegative_scalar(flag_threshold,
                                                     "flag_threshold")
  scale_limit <- .validate_dff_positive_scalar(scale_limit, "scale_limit")
  if (!is.character(flag_color) || length(flag_color) != 1L ||
      is.na(flag_color) || !nzchar(flag_color)) {
    stop("`flag_color` must be a single non-empty character string.",
         call. = FALSE)
  }

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

  method <- x$config$method %||% NA_character_
  classification_system <- NA_character_
  if (inherits(x, "mfrm_dif") && !is.null(x$dif_table) &&
      "ClassificationSystem" %in% names(x$dif_table)) {
    cs <- unique(as.character(x$dif_table$ClassificationSystem))
    cs <- cs[!is.na(cs)]
    classification_system <- cs[1] %||% NA_character_
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
  flag_matrix <- if (!is.null(flag_threshold)) {
    abs(mat) >= flag_threshold
  } else {
    matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat),
           dimnames = dimnames(mat))
  }
  flag_matrix[is.na(flag_matrix)] <- FALSE

  if (draw) {
    # Color scale: blue-white-red
    n_colors <- 64
    max_abs <- scale_limit %||% max(abs(mat), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)
    blue_white_red <- grDevices::colorRampPalette(
      c("steelblue", "white", "firebrick")
    )(n_colors)
    mat_draw <- pmax(pmin(mat, max_abs), -max_abs)

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
      z = t(mat_draw),
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

    if (!is.null(flag_threshold) && any(flag_matrix, na.rm = TRUE)) {
      flagged <- which(flag_matrix, arr.ind = TRUE)
      for (ii in seq_len(nrow(flagged))) {
        ri <- flagged[ii, 1]
        ci <- flagged[ii, 2]
        graphics::rect(ci - 0.5, ri - 0.5, ci + 0.5, ri + 0.5,
                       border = flag_color, lwd = 1.2)
      }
    }

    if (isTRUE(show_values)) {
      for (ri in seq_len(nrow(mat))) {
        for (ci in seq_len(ncol(mat))) {
          val <- mat[ri, ci]
          if (is.finite(val)) {
            graphics::text(ci, ri, formatC(val, format = "f",
                                           digits = value_digits), cex = 0.6)
          }
        }
      }
    }
  }

  metric_label <- switch(metric,
    obs_exp = "Obs - Exp Average",
    t       = "Standardized Residual (t)",
    contrast = "Differential-Functioning Contrast"
  )
  out <- new_mfrm_plot_data(
    "dif_heatmap",
    list(
      matrix = mat,
      flag_matrix = flag_matrix,
      pairs = as.data.frame(tbl, stringsAsFactors = FALSE),
      metric = metric,
      value_column = value_col,
      title = paste("DFF Heatmap:", metric_label),
      subtitle = sprintf("%d row x %d column matrix; metric column = `%s`",
                         nrow(mat), ncol(mat), value_col),
      thresholds = data.frame(
        Metric = metric,
        Threshold = flag_threshold %||% NA_real_,
        Rule = if (!is.null(flag_threshold)) {
          paste0("abs(value) >= ", flag_threshold)
        } else {
          "No plotting threshold supplied"
        },
        stringsAsFactors = FALSE
      ),
      interpretation_guide = .dff_interpretation_guide(
        metric = metric,
        method = method,
        classification_system = classification_system,
        flag_threshold = flag_threshold
      ),
      gpcm_boundary = x$gpcm_boundary %||% data.frame(),
      settings = list(
        show_values = show_values,
        value_digits = value_digits,
        flag_threshold = flag_threshold,
        scale_limit = scale_limit,
        flag_color = flag_color
      )
    )
  )
  invisible(out)
}

# ============================================================================
# C. Information Function Computation and Plotting
# ============================================================================

information_build_step_structure <- function(fit, model) {
  steps <- tibble::as_tibble(fit$steps %||% tibble::tibble())
  if (nrow(steps) == 0 || !"Estimate" %in% names(steps)) {
    stop("Step/threshold estimates are required for information computation.",
         call. = FALSE)
  }

  if (identical(model, "RSM")) {
    step_est <- suppressWarnings(as.numeric(steps$Estimate))
    if (length(step_est) == 0 || any(!is.finite(step_est))) {
      stop("RSM step estimates must be finite numeric values.",
           call. = FALSE)
    }
    return(list(
      kind = "common",
      categories = 0:length(step_est),
      step_facet = NULL,
      step_levels = NULL,
      compute = function(eta) category_prob_rsm(eta, c(0, cumsum(step_est)))
    ))
  }

  if (!identical(model, "PCM") && !identical(model, "GPCM")) {
    stop(
      "`compute_information()` currently supports only ordered-category ",
      "`model = \"RSM\"`, `model = \"PCM\"`, and bounded ",
      "`model = \"GPCM\"` fits.",
      call. = FALSE
    )
  }

  step_facet <- as.character(fit$config$step_facet %||% NA_character_)
  if (!nzchar(step_facet)) {
    stop("PCM information requires a valid `step_facet` in `fit$config`.",
         call. = FALSE)
  }
  if (!"StepFacet" %in% names(steps)) {
    stop("PCM step estimates must include a `StepFacet` column.",
         call. = FALSE)
  }
  if (!"StepIndex" %in% names(steps)) {
    if (!"Step" %in% names(steps)) {
      stop("PCM step estimates must include either `StepIndex` or `Step`.",
           call. = FALSE)
    }
    steps$StepIndex <- suppressWarnings(as.integer(gsub("[^0-9]+", "", as.character(steps$Step))))
  }

  steps <- steps |>
    dplyr::transmute(
      StepFacet = as.character(.data$StepFacet),
      StepIndex = suppressWarnings(as.integer(.data$StepIndex)),
      Estimate = suppressWarnings(as.numeric(.data$Estimate))
    ) |>
    dplyr::arrange(.data$StepFacet, .data$StepIndex)

  if (any(is.na(steps$StepFacet) | !nzchar(steps$StepFacet))) {
    stop("PCM `StepFacet` labels must be non-empty strings.",
         call. = FALSE)
  }
  if (any(!is.finite(steps$StepIndex)) || any(steps$StepIndex < 1L)) {
    stop("PCM `StepIndex` values must be positive integers.",
         call. = FALSE)
  }
  if (any(!is.finite(steps$Estimate))) {
    stop("PCM step estimates must be finite numeric values.",
         call. = FALSE)
  }

  step_counts <- steps |>
    dplyr::count(.data$StepFacet, name = "n_steps")
  if (length(unique(step_counts$n_steps)) != 1L) {
    stop("Each PCM `StepFacet` level must supply the same number of thresholds.",
         call. = FALSE)
  }
  n_steps <- unique(step_counts$n_steps)
  expected_index <- seq_len(n_steps)
  bad_order <- steps |>
    dplyr::group_by(.data$StepFacet) |>
    dplyr::summarize(ok = identical(sort(unique(.data$StepIndex)), expected_index), .groups = "drop")
  if (any(!bad_order$ok)) {
    stop("PCM step estimates must provide contiguous `StepIndex` values starting at 1 for each `StepFacet` level.",
         call. = FALSE)
  }

  step_levels <- as.character(fit$prep$levels[[step_facet]] %||% character())
  observed_levels <- unique(as.character(steps$StepFacet))
  if (length(step_levels) == 0L) {
    step_levels <- observed_levels
  } else {
    step_levels <- c(step_levels, setdiff(observed_levels, step_levels))
  }

  step_mat <- matrix(NA_real_, nrow = length(step_levels), ncol = n_steps)
  rownames(step_mat) <- step_levels
  for (i in seq_len(nrow(steps))) {
    row_idx <- match(steps$StepFacet[i], step_levels)
    col_idx <- steps$StepIndex[i]
    step_mat[row_idx, col_idx] <- steps$Estimate[i]
  }
  if (any(!is.finite(step_mat))) {
    stop("PCM step estimates did not cover all expected `StepFacet` / `StepIndex` combinations.",
         call. = FALSE)
  }
  step_cum_mat <- t(apply(step_mat, 1, function(x) c(0, cumsum(x))))

  if (identical(model, "PCM")) {
    return(list(
      kind = "step_facet_specific",
      categories = 0:n_steps,
      step_facet = step_facet,
      step_levels = step_levels,
      slope_facet = NULL,
      slope_levels = NULL,
      slopes = NULL,
      compute = function(eta, step_idx) {
        category_prob_pcm(
          eta = eta,
          step_cum_mat = step_cum_mat,
          criterion_idx = step_idx
        )
      }
    ))
  }

  slope_facet <- as.character(fit$config$slope_facet %||% NA_character_)
  if (!nzchar(slope_facet)) {
    stop("GPCM information requires a valid `slope_facet` in `fit$config`.",
         call. = FALSE)
  }
  slope_tbl <- tibble::as_tibble(fit$slopes %||% tibble::tibble())
  if (nrow(slope_tbl) == 0 || !all(c("SlopeFacet", "Estimate") %in% names(slope_tbl))) {
    stop("GPCM information requires a `fit$slopes` table with `SlopeFacet` and `Estimate` columns.",
         call. = FALSE)
  }
  slope_tbl <- slope_tbl |>
    dplyr::transmute(
      SlopeFacet = as.character(.data$SlopeFacet),
      Estimate = suppressWarnings(as.numeric(.data$Estimate))
    )
  if (any(is.na(slope_tbl$SlopeFacet) | !nzchar(slope_tbl$SlopeFacet))) {
    stop("GPCM `SlopeFacet` labels must be non-empty strings.",
         call. = FALSE)
  }
  if (any(!is.finite(slope_tbl$Estimate)) || any(slope_tbl$Estimate <= 0)) {
    stop("GPCM slopes must be finite and strictly positive for information computation.",
         call. = FALSE)
  }

  slope_levels <- as.character(fit$prep$levels[[slope_facet]] %||% character())
  observed_slope_levels <- unique(as.character(slope_tbl$SlopeFacet))
  if (length(slope_levels) == 0L) {
    slope_levels <- observed_slope_levels
  } else {
    slope_levels <- c(slope_levels, setdiff(observed_slope_levels, slope_levels))
  }
  slope_vals <- rep(NA_real_, length(slope_levels))
  names(slope_vals) <- slope_levels
  slope_vals[slope_tbl$SlopeFacet] <- slope_tbl$Estimate
  if (any(!is.finite(slope_vals))) {
    stop("GPCM slopes did not cover all expected `SlopeFacet` levels.",
         call. = FALSE)
  }

  list(
    kind = "step_and_slope_specific",
    categories = 0:n_steps,
    step_facet = step_facet,
    step_levels = step_levels,
    slope_facet = slope_facet,
    slope_levels = slope_levels,
    slopes = unname(slope_vals),
    compute = function(eta, step_idx, slope_idx) {
      category_prob_gpcm(
        eta = eta,
        step_cum_mat = step_cum_mat,
        criterion_idx = step_idx,
        slopes = unname(slope_vals),
        slope_idx = slope_idx
      )
    }
  )
}

#' Compute design-weighted precision curves for ordered many-facet fits
#'
#' Calculates design-weighted score-variance curves across the latent
#' trait (theta) for a fitted ordered-category `RSM`, `PCM`, or bounded
#' `GPCM` model. Returns both an overall precision curve (`$tif`) and
#' per-facet-level contribution curves (`$iif`) based on the realized
#' observation pattern.
#'
#' @param fit Output from [fit_mfrm()].
#' @param theta_range Numeric vector of length 2 giving the range of theta
#'   values. Default `c(-6, 6)`.
#' @param theta_points Integer number of points at which to evaluate
#'   information. Default `201`.
#'
#' @details
#' For `RSM` / `PCM`, the score variance at theta for one observed design cell
#' is:
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
#' In an ordered Rasch-family model, category discrimination is fixed at 1, so
#' this score-variance representation is the natural conditional information
#' identity rather than a separate approximation. For binary data it reduces to
#' the familiar \eqn{p(\theta)\{1 - p(\theta)\}} form. For `PCM`, the package
#' evaluates each observed design cell using the threshold vector associated
#' with that cell's realized `step_facet` level. For bounded `GPCM`, the
#' same design-weighted score variance is scaled by the squared discrimination
#' attached to the realized `slope_facet` level, which is the
#' \eqn{a_j^2 \cdot \mathrm{Var}(T \mid \theta)} item-information identity that
#' Muraki (1993, Equation 10) derives by applying Samejima's (1974)
#' polytomous information formula to the GPCM kernel of Muraki (1992).
#'
#' @section What `tif` and `iif` mean here:
#' In `mfrmr`, this helper supports ordered-category `RSM`, `PCM`, and the
#' current bounded `GPCM` fit. The total curve (`$tif`) is the sum of
#' design-weighted cell contributions across all non-person facet levels in the
#' fitted model. The facet-level contribution curves (`$iif`) keep those
#' weighted contributions separated, so you can see which observed rater
#' levels, criteria, or other facet levels are driving precision at different
#' parts of the scale. For `PCM`, step-facet-specific thresholds are respected
#' when each observed design cell is evaluated. For bounded `GPCM`, those
#' same cell-level variances are additionally scaled by the squared
#' discrimination associated with the realized `slope_facet` level.
#'
#' @section What this quantity does not justify:
#' - It is not a textbook many-facet test-information function for an abstract
#'   fixed item set.
#' - It should not be used as if it were design-free evidence about a form's
#'   precision independent of the realized observation pattern.
#' - It does not currently extend beyond the ordered-category `RSM` / `PCM` /
#'   bounded `GPCM` family implemented by [fit_mfrm()].
#'
#' @section When to use this:
#' Use `compute_information()` when you want a design-weighted precision screen
#' for an `RSM`, `PCM`, or bounded `GPCM` fit along the latent
#' continuum. In practice:
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
#' @section References:
#' The ordered-category probability structures come from Andrich's `RSM`
#' formulation and Masters' `PCM`. The bounded `GPCM` information identity
#' \eqn{a_j^2 \cdot \mathrm{Var}(T \mid \theta)} is derived in Muraki
#' (1993, Equation 10) by applying Samejima's (1974) general polytomous
#' information formula \eqn{I_j(\theta) = \sum_k P_{jk}(\theta)
#' [-\partial^2 \ln P_{jk} / \partial \theta^2]} to the GPCM probability
#' kernel of Muraki (1992). For the integer scoring function
#' \eqn{T_k = k} used by `mfrmr`, this reduces to
#' \eqn{a_j^2 \cdot \mathrm{Var}(K \mid \theta)}. In `mfrmr`, those formulas
#' are applied to the realized many-facet observation design, so the output
#' should be read as a design-weighted precision summary rather than as a
#' design-free abstract test function.
#'
#' - Andrich, D. (1978). *A rating formulation for ordered response
#'   categories*. Psychometrika, 43(4), 561-573.
#' - Masters, G. N. (1982). *A Rasch model for partial credit scoring*.
#'   Psychometrika, 47(2), 149-174.
#' - Muraki, E. (1992). *A generalized partial credit model: Application
#'   of an EM algorithm*. Applied Psychological Measurement, 16(2),
#'   159-176. \doi{10.1177/014662169201600206} (See Equations 6, 10, and
#'   13 for the probability kernel and the
#'   \eqn{\partial P_k / \partial \theta = a_j P_k (k - E[K])}
#'   derivative used by all GPCM helpers in `mfrmr`.)
#' - Muraki, E. (1993). *Information functions of the generalized
#'   partial credit model*. Applied Psychological Measurement, 17(4),
#'   351-363. \doi{10.1177/014662169301700403} (Equation 10 derives the
#'   item information function for the GPCM,
#'   \eqn{I_j(\theta) = D^2 a_j^2 \mathrm{Var}(T \mid \theta)}, by
#'   applying Samejima's (1974) polytomous information formula to the
#'   GPCM kernel; this is the canonical reference for `compute_information()`
#'   under bounded `GPCM`.)
#' - Samejima, F. (1974). *Normal ogive model on the continuous
#'   response level in the multidimensional latent space*.
#'   Psychometrika, 39, 111-121. (Source for the general polytomous
#'   information formula that Muraki 1993 specializes to the GPCM.)
#'
#' @section Interpreting output:
#' - `$tif`: design-weighted precision curve data with theta, Information, and SE.
#' - `$iif`: design-weighted facet-level contribution curves for the fitted
#'   non-person facets.
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
#' be used in reporting, follow with [precision_review_report()].
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
#'                  method = "JML", model = "RSM", maxit = 30)
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

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)

  # Extract model parameters and realized observation design.
  step_structure <- information_build_step_structure(fit, model)
  categories <- step_structure$categories
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
  compute_cell_info <- function(offset, step_idx = NULL, slope_idx = NULL) {
    eta <- theta_grid + offset
    probs <- if (identical(step_structure$kind, "common")) {
      step_structure$compute(eta)
    } else if (identical(step_structure$kind, "step_facet_specific")) {
      step_structure$compute(eta, rep(step_idx, length(eta)))
    } else {
      step_structure$compute(eta, rep(step_idx, length(eta)), rep(slope_idx, length(eta)))
    }
    expected <- as.vector(probs %*% category_vec)
    second_moment <- as.vector(probs %*% category_sq_vec)
    variance <- pmax(second_moment - expected^2, 0)
    if (identical(step_structure$kind, "step_and_slope_specific")) {
      variance <- variance * (step_structure$slopes[slope_idx]^2)
    }
    variance
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

  cell_step_idx <- NULL
  cell_slope_idx <- NULL
  if (identical(step_structure$kind, "step_facet_specific") ||
      identical(step_structure$kind, "step_and_slope_specific")) {
    step_facet <- step_structure$step_facet
    if (!step_facet %in% names(design_cells)) {
      stop("Step facet column was not found in the observed design cells.",
           call. = FALSE)
    }
    cell_step_idx <- match(as.character(design_cells[[step_facet]]), step_structure$step_levels)
    cell_ok <- cell_ok & is.finite(cell_step_idx)
  }
  if (identical(step_structure$kind, "step_and_slope_specific")) {
    slope_facet <- step_structure$slope_facet
    if (!slope_facet %in% names(design_cells)) {
      stop("GPCM slope facet column was not found in the observed design cells.",
           call. = FALSE)
    }
    cell_slope_idx <- match(as.character(design_cells[[slope_facet]]), step_structure$slope_levels)
    cell_ok <- cell_ok & is.finite(cell_slope_idx)
  }

  design_cells <- design_cells[cell_ok, , drop = FALSE]
  cell_offset <- cell_offset[cell_ok]
  if (!is.null(cell_step_idx)) {
    cell_step_idx <- cell_step_idx[cell_ok]
  }
  if (!is.null(cell_slope_idx)) {
    cell_slope_idx <- cell_slope_idx[cell_ok]
  }
  if (nrow(design_cells) == 0) {
    stop("No valid observed design cells were available for information computation.",
         call. = FALSE)
  }

  if (is.null(cell_step_idx)) {
    info_mat <- vapply(cell_offset, compute_cell_info, numeric(length(theta_grid)))
  } else if (is.null(cell_slope_idx)) {
    info_mat <- vapply(
      seq_along(cell_offset),
      function(i) compute_cell_info(cell_offset[i], cell_step_idx[i]),
      numeric(length(theta_grid))
    )
  } else {
    info_mat <- vapply(
      seq_along(cell_offset),
      function(i) compute_cell_info(cell_offset[i], cell_step_idx[i], cell_slope_idx[i]),
      numeric(length(theta_grid))
    )
  }
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
  conditional_sem <- tibble(
    Theta = theta_grid,
    ConditionalSEM = tif$SE,
    Information = tif$Information
  )
  information_long <- dplyr::bind_rows(
    tibble(
      PlotType = "tif",
      Metric = "Information",
      Series = "Information",
      Theta = theta_grid,
      Value = tif$Information,
      ValueName = "Information",
      DisplayedByDefault = TRUE
    ),
    tibble(
      PlotType = "conditional_sem",
      Metric = "Conditional SEM",
      Series = "Conditional SEM",
      Theta = theta_grid,
      Value = tif$SE,
      ValueName = "ConditionalSEM",
      DisplayedByDefault = TRUE
    )
  )
  finite_info <- is.finite(tif$Information)
  finite_sem <- is.finite(tif$SE)
  summary_tbl <- tibble(
    ThetaAtMaxInformation = if (any(finite_info)) tif$Theta[which.max(ifelse(finite_info, tif$Information, -Inf))] else NA_real_,
    MaxInformation = if (any(finite_info)) max(tif$Information[finite_info], na.rm = TRUE) else NA_real_,
    ThetaAtMinConditionalSEM = if (any(finite_sem)) tif$Theta[which.min(ifelse(finite_sem, tif$SE, Inf))] else NA_real_,
    MinConditionalSEM = if (any(finite_sem)) min(tif$SE[finite_sem], na.rm = TRUE) else NA_real_,
    ThetaPoints = length(theta_grid),
    ThetaMin = min(theta_grid, na.rm = TRUE),
    ThetaMax = max(theta_grid, na.rm = TRUE)
  )

  out <- list(
    tif = tif,
    conditional_sem = conditional_sem,
    information_long = information_long,
    iif = iif,
    summary = summary_tbl,
    theta_range = theta_range
  )
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
#'   facet-level contribution curves, `"se"` / `"sem"` / `"csem"` for the
#'   conditional standard error of measurement implied by that curve, or
#'   `"both"` for precision with conditional SEM on a secondary axis.
#' @param facet For `type = "iif"`, which facet to plot. If `NULL`,
#'   the first facet is used.
#' @param draw If `TRUE` (default), draw the plot. If `FALSE`, return
#'   reusable `mfrm_plot_data` invisibly.
#' @param ... Additional graphical parameters.
#'
#' @section Plot types:
#' - `"tif"`: overall design-weighted precision across theta.
#' - `"se"` / `"sem"` / `"csem"`: conditional SEM across theta.
#' - `"both"`: precision and conditional SEM together, useful for presentations.
#' - `"iif"`: facet-level contribution curves for one selected facet in a
#'   supported `RSM`, `PCM`, or bounded `GPCM` fit.
#'
#' @section Which type should I use?:
#' - Use `"tif"` for a quick overall read on precision.
#' - Use `"sem"` or `"csem"` when standard-error language is easier to
#'   communicate than precision.
#' - Use `"both"` when you want both views in one figure.
#' - Use `"iif"` when you want to see which facet levels are shaping the total
#'   precision curve.
#'
#' @section Interpreting output:
#' - The total curve peaks where the realized design is most precise.
#' - Conditional SEM is derived as `1 / sqrt(precision)`; lower is better.
#' - Facet-level curves show which facet levels contribute most to that
#'   realized precision at each theta.
#' - For bounded `GPCM`, those contributions include the squared
#'   discrimination scaling implied by the fitted `slope_facet`.
#' - If the precision peak sits far from the bulk of person measures, the
#'   realized design may be poorly targeted.
#'
#' @section Returned data when draw = FALSE:
#' `draw = FALSE` returns an `mfrm_plot_data` object. The underlying plotting
#' data are stored in `$data$plot`. For `type = "tif"`, `"se"`, or `"both"`,
#' those rows come from `x$tif`. For `type = "iif"`, the returned rows come
#' from `x$iif` filtered to the requested facet. The plot data also include
#' `plot_long`, `information_long`, `conditional_sem`, `summary`, and
#' `settings` so ggplot2, plotly, Quarto, and table workflows can reuse the
#' information and conditional-SEM series without parsing the drawn figure.
#'
#' @section Typical workflow:
#' 1. Compute information with [compute_information()].
#' 2. Plot with `plot_information(info)` for the total precision curve.
#' 3. Use `plot_information(info, type = "iif", facet = "Rater")` for
#'    facet-level contributions.
#' 4. Use `draw = FALSE` when you want reusable plot data for custom graphics
#'    or reporting helpers.
#'
#' @return Invisibly, an `mfrm_plot_data` object.
#'
#' @seealso [compute_information()], [fit_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 30)
#' info <- compute_information(fit)
#' tif_data <- plot_information(info, type = "tif", draw = FALSE)
#' head(tif_data$data$plot)
#' iif_data <- plot_information(info, type = "iif", facet = "Rater", draw = FALSE)
#' head(iif_data$data$plot)
#' @export
plot_information <- function(x,
                             type = c("tif", "iif", "se", "sem", "csem", "both"),
                             facet = NULL,
                             draw = TRUE,
                             ...) {
  if (!inherits(x, "mfrm_information")) {
    stop("`x` must be an `mfrm_information` object.", call. = FALSE)
  }
  requested_type <- match.arg(type)
  type <- if (requested_type %in% c("sem", "csem")) "se" else requested_type
  plot_name <- if (requested_type %in% c("sem", "csem")) "information_sem" else switch(type,
    tif = "information_tif",
    se = "information_se",
    both = "information_tif_se",
    iif = "information_iif"
  )
  title <- switch(type,
    tif = "Design-weighted precision curve",
    se = "Conditional SEM curve",
    both = "Design-weighted precision and conditional SEM",
    iif = paste("Facet-level precision contributions:", facet %||% unique(x$iif$Facet)[1])
  )
  subtitle <- switch(type,
    tif = "Overall precision across the latent continuum",
    se = "Conditional standard error of measurement implied by the design-weighted precision curve",
    both = "Precision and conditional SEM on a shared theta grid",
    iif = "Contribution curves for the selected facet"
  )
  legend <- switch(type,
    tif = new_plot_legend("Information (precision)", "precision", "line", "steelblue"),
    se = new_plot_legend("Conditional SEM", "conditional_sem", "line", "coral"),
    both = new_plot_legend(
      label = c("Information (precision)", "Conditional SEM"),
      role = c("precision", "conditional_sem"),
      aesthetic = c("line", "line"),
      value = c("steelblue", "coral")
    ),
    iif = new_plot_legend("Facet level contributions", "facet_level", "line", "rainbow")
  )
  series <- switch(type,
    tif = "Information",
    se = "ConditionalSEM",
    both = c("Information", "ConditionalSEM"),
    iif = "Information"
  )
  payload <- function(data, title_override = title, subtitle_override = subtitle) {
    full_long <- as.data.frame(x$information_long %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(full_long) == 0L && !is.null(x$tif)) {
      full_long <- dplyr::bind_rows(
        tibble(
          PlotType = "tif",
          Metric = "Information",
          Series = "Information",
          Theta = x$tif$Theta,
          Value = x$tif$Information,
          ValueName = "Information",
          DisplayedByDefault = TRUE
        ),
        tibble(
          PlotType = "conditional_sem",
          Metric = "Conditional SEM",
          Series = "Conditional SEM",
          Theta = x$tif$Theta,
          Value = x$tif$SE,
          ValueName = "ConditionalSEM",
          DisplayedByDefault = TRUE
        )
      )
    }
    if (nrow(full_long) > 0L) {
      full_long$DisplayedByDefault <- switch(
        type,
        tif = full_long$ValueName == "Information",
        se = full_long$ValueName == "ConditionalSEM",
        both = full_long$ValueName %in% c("Information", "ConditionalSEM"),
        iif = FALSE
      )
    }
    iif_long <- if (identical(type, "iif")) {
      data.frame(
        PlotType = "iif",
        Metric = "Facet-level information contribution",
        Series = as.character(data$Level),
        Theta = suppressWarnings(as.numeric(data$Theta)),
        Value = suppressWarnings(as.numeric(data$Information)),
        ValueName = "InformationContribution",
        Facet = as.character(data$Facet),
        Level = as.character(data$Level),
        Exposure = suppressWarnings(as.numeric(data$Exposure)),
        DisplayedByDefault = TRUE,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame()
    }
    plot_long <- if (identical(type, "iif")) iif_long else full_long
    conditional_sem <- as.data.frame(x$conditional_sem %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(conditional_sem) == 0L && !is.null(x$tif)) {
      conditional_sem <- data.frame(
        Theta = x$tif$Theta,
        ConditionalSEM = x$tif$SE,
        Information = x$tif$Information,
        stringsAsFactors = FALSE
      )
    }
    new_mfrm_plot_data(
      plot_name,
      list(
        plot = data,
        plot_long = plot_long,
        information_long = full_long,
        conditional_sem = conditional_sem,
        summary = as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE),
        settings = data.frame(
          RequestedType = requested_type,
          CanonicalType = type,
          Facet = as.character(facet %||% NA_character_),
          Draw = isTRUE(draw),
          stringsAsFactors = FALSE
        ),
        series = series,
        title = title_override,
        subtitle = subtitle_override,
        legend = legend,
        reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
      )
    )
  }

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
        graphics::mtext("Conditional SEM", side = 4, line = 2.5, col = "coral")
        graphics::legend("topright",
                         legend = c("Information (precision)", "Conditional SEM"),
                         col = c("steelblue", "coral"),
                         lty = c(1, 2), lwd = 2, bty = "n")
      }
    }
    invisible(payload(plot_data))
  } else if (type == "se") {
    plot_data <- x$tif
    if (draw) {
      plot(plot_data$Theta, plot_data$SE,
           type = "l", lwd = 2, col = "coral",
           xlab = expression(theta), ylab = "Conditional SEM",
           main = "Conditional SEM from Design-Weighted Precision", ...)
      graphics::grid()
    }
    invisible(payload(plot_data))
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
    invisible(payload(plot_data,
      title_override = paste("Facet-level precision contributions:", facet),
      subtitle_override = "Contribution curves for the selected facet"
    ))
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
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
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
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept Wright maps
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(toy_small, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 30)
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
                                preset = c("standard", "publication", "compact", "monochrome"),
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
# E. Anchoring & Equating Workflow
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
#' Re-estimates a fitted many-facet model on new data while holding selected
#' facet parameters fixed at the values from a previous (baseline) calibration.
#' This is the standard workflow for placing new data onto an existing scale,
#' linking test forms, or carrying a baseline calibration across
#' administration windows.
#' For bounded `GPCM`, treat this as direct exploratory anchor/drift support
#' rather than as the package's formal linking-synthesis route.
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
#'   reviewing which elements were constrained.
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
#'                  method = "JML", maxit = 30)
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
#'     \item{common_vs_reference}{Tibble of common-element counts
#'       between each wave and the reference wave (i.e., which
#'       elements remain comparable across the entire chain).}
#'     \item{n_common_all_waves}{Integer count of elements that are
#'       common across every wave; used by `summary()` to gauge how
#'       robust the chain is to chained linking error.}
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
#'                  method = "JML", maxit = 30)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
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
        weak_rows <- support_tbl[!support_tbl$LinkSupportAdequate, , drop = FALSE]
        weak_detail <- if ("CommonElements" %in% names(weak_rows)) {
          paste0(weak_rows$Facet, " (", weak_rows$CommonElements, "/", support_guideline, ")")
        } else {
          as.character(weak_rows$Facet)
        }
        warning(
          sprintf(
            "Thin linking support between '%s' and '%s': fewer than %d retained common elements in %s.",
            wave_names[ref_idx],
            wave_names[i],
            support_guideline,
            paste(weak_detail, collapse = ", ")
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

  # Common elements count (pairwise, every combination; existing).
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

  # Reference-vs-each-wave view for 3+ wave reviews. The existing
  # `common_elements` table enumerates all pairs, which is fine for 2
  # waves but noisy for 3+. Add a dedicated table focused on the
  # reference wave, and record the count common to ALL waves.
  common_vs_reference <- if (length(fits) >= 2L) {
    ref_table <- do.call(rbind, lapply(seq_along(fits), function(i) {
      if (i == ref_idx) return(NULL)
      n_common <- nrow(dplyr::inner_join(
        est_list[[ref_idx]], est_list[[i]], by = c("Facet", "Level")
      ))
      data.frame(
        Reference = wave_names[ref_idx],
        Wave = wave_names[i],
        N_Common = as.integer(n_common),
        stringsAsFactors = FALSE
      )
    }))
    if (is.null(ref_table)) {
      data.frame(Reference = character(0), Wave = character(0),
                 N_Common = integer(0), stringsAsFactors = FALSE)
    } else {
      ref_table
    }
  } else {
    data.frame(Reference = character(0), Wave = character(0),
               N_Common = integer(0), stringsAsFactors = FALSE)
  }
  n_common_all_waves <- if (length(est_list) >= 2L) {
    shared <- Reduce(
      function(a, b) dplyr::inner_join(a, b, by = c("Facet", "Level")),
      lapply(est_list, function(e) e[, c("Facet", "Level")])
    )
    as.integer(nrow(shared))
  } else NA_integer_

  out <- list(
    drift_table = drift_table, summary = drift_summary,
    common_elements = common_counts,
    common_vs_reference = common_vs_reference,
    n_common_all_waves = n_common_all_waves,
    common_by_facet = if (length(support_rows) > 0) dplyr::bind_rows(support_rows) else tibble::tibble(),
    config = list(reference = wave_names[ref_idx],
                  method = "screened_common_element_alignment",
                  intended_use = "review_screen",
                  models = unique(stats::na.omit(vapply(
                    fits,
                    function(f) {
                      as.character(f$config$model %||% f$summary$Model[1] %||% NA_character_)
                    },
                    character(1)
                  ))),
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
#'                  method = "JML", maxit = 30)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
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
        weak_rows <- support_tbl[!support_tbl$LinkSupportAdequate, , drop = FALSE]
        weak_detail <- if ("CommonElements" %in% names(weak_rows)) {
          paste0(weak_rows$Facet, " (", weak_rows$CommonElements, "/", support_guideline, ")")
        } else {
          as.character(weak_rows$Facet)
        }
        warning(
          sprintf(
            "Thin linking support between '%s' and '%s': fewer than %d retained common elements in %s.",
            wave_names[i],
            wave_names[i + 1],
            support_guideline,
            paste(weak_detail, collapse = ", ")
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
                  models = unique(stats::na.omit(vapply(
                    fits,
                    function(f) {
                      as.character(f$config$model %||% f$summary$Model[1] %||% NA_character_)
                    },
                    character(1)
                  ))),
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
#' @param y Unused (S3 plot signature requirement).
#' @param type One of `"graph"` (bipartite Wave x anchor-element graph;
#'   requires the `igraph` package), `"common_anchors"` (default; bar
#'   chart of common-anchor counts per wave pair), or `"chain"`.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw the plot with base graphics.
#' @export
plot.mfrm_equating_chain <- function(x, y = NULL,
                                     type = c("common_anchors", "graph", "chain"),
                                     preset = c("standard", "publication", "compact", "monochrome"),
                                     draw = TRUE, ...) {
  if (!inherits(x, "mfrm_equating_chain")) {
    stop("`x` must be an mfrm_equating_chain object.", call. = FALSE)
  }
  type <- match.arg(type)
  style <- resolve_plot_preset(preset)
  if (identical(type, "chain")) {
    if (isTRUE(draw)) apply_plot_preset(style)
    return(.plot_equating_chain(x, draw = draw, style = style, ...))
  }
  if (identical(type, "graph")) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      message("`plot(..., type = \"graph\")` requires the `igraph` package ",
              "(in Suggests). Falling back to type = \"common_anchors\".")
      type <- "common_anchors"
    } else {
      detail <- as.data.frame(x$element_detail %||% data.frame(),
                              stringsAsFactors = FALSE)
      if (nrow(detail) == 0L ||
          !all(c("Facet", "Level") %in% names(detail))) {
        stop("Graph view requires Facet/Level columns in element_detail.",
             call. = FALSE)
      }
      detail$AnchorId <- paste0(detail$Facet, ":", detail$Level)
      if ("Wave" %in% names(detail)) {
        edges <- unique(detail[, c("Wave", "AnchorId")])
      } else if ("Link" %in% names(detail)) {
        # Parse "<from> -> <to>" / "<from> | <to>" / "<from> <-> <to>"
        # link strings into both endpoints so each wave-anchor pair
        # becomes an edge.
        link_str <- as.character(detail$Link)
        parts <- strsplit(link_str, "\\s*(->|<->|<-|\\|)\\s*", perl = TRUE)
        endpoints <- do.call(rbind, lapply(seq_along(parts), function(i) {
          p <- parts[[i]]
          if (length(p) >= 2L) {
            data.frame(
              Wave = c(p[1], p[2]),
              AnchorId = rep(detail$AnchorId[i], 2L),
              stringsAsFactors = FALSE
            )
          } else if (length(p) == 1L) {
            data.frame(Wave = p[1], AnchorId = detail$AnchorId[i],
                       stringsAsFactors = FALSE)
          } else NULL
        }))
        edges <- if (!is.null(endpoints)) unique(endpoints) else
          data.frame(Wave = character(0), AnchorId = character(0),
                     stringsAsFactors = FALSE)
      } else {
        stop("Graph view requires either a Wave column or a Link column ",
             "(parseable as 'WaveA -> WaveB').", call. = FALSE)
      }
      g <- igraph::graph_from_data_frame(
        d = edges, directed = FALSE,
        vertices = data.frame(
          name = c(unique(edges$Wave), unique(edges$AnchorId)),
          type = c(rep(TRUE, length(unique(edges$Wave))),
                   rep(FALSE, length(unique(edges$AnchorId))))
        )
      )
      out <- new_mfrm_plot_data(
        "equating_chain_graph",
        list(
          data = list(edges = edges, n_waves = length(unique(edges$Wave)),
                       n_anchors = length(unique(edges$AnchorId))),
          title = "Equating chain (bipartite graph)",
          subtitle = sprintf("%d wave(s), %d anchor element(s)",
                              length(unique(edges$Wave)),
                              length(unique(edges$AnchorId))),
          preset = style$name
        )
      )
      if (isTRUE(draw)) {
        apply_plot_preset(style)
        igraph::V(g)$color <- ifelse(igraph::V(g)$type,
                                      style$accent_primary,
                                      style$accent_tertiary)
        igraph::V(g)$shape <- ifelse(igraph::V(g)$type, "square", "circle")
        igraph::V(g)$label.cex <- 0.7
        plot(g, layout = igraph::layout_as_bipartite(g),
             vertex.size = 14, vertex.label.color = "black",
             edge.color = style$grid, edge.width = 1,
             main = "Equating chain (bipartite graph)")
      }
      return(invisible(out))
    }
  }
  # type == "common_anchors": bar chart of pairwise common anchor counts.
  links <- as.data.frame(x$links %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(links) == 0L) {
    stop("`x$links` is empty; nothing to plot.", call. = FALSE)
  }
  # Tolerate both From/To and Wave1/Wave2 column names.
  if (all(c("From", "To") %in% names(links))) {
    names(links)[names(links) == "From"] <- "Wave1"
    names(links)[names(links) == "To"] <- "Wave2"
  }
  if (!all(c("Wave1", "Wave2", "N_Common") %in% names(links))) {
    stop("Common-anchor view requires links table with Wave1/Wave2/N_Common ",
         "(or From/To/N_Common).", call. = FALSE)
  }
  links$Pair <- paste(links$Wave1, links$Wave2, sep = " <-> ")
  out <- new_mfrm_plot_data(
    "equating_chain_common_anchors",
    list(
      data = links,
      title = "Common anchors per wave pair",
      subtitle = sprintf("%d pairwise comparison(s)", nrow(links)),
      preset = style$name
    )
  )
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(8, 5, 3, 1))
    graphics::barplot(
      height = links$N_Common,
      names.arg = links$Pair,
      las = 2,
      col = style$accent_primary,
      main = "Common anchors per wave pair",
      ylab = "N common"
    )
  }
  invisible(out)
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

# --- build_linking_review ----------------------------------------------------

.validate_linking_review_input <- function(x, arg, classes) {
  if (is.null(x)) return(NULL)
  if (!inherits(x, classes)) {
    stop(
      "`", arg, "` must be one of: ", paste(classes, collapse = ", "), ".",
      call. = FALSE
    )
  }
  x
}

.linking_review_support_status <- function(source_models) {
  gpcm_detected <- any(identical(source_models, "GPCM") | source_models == "GPCM")
  tibble::tibble(
    Scope = c("RSM / PCM", "bounded GPCM"),
    Status = c("supported", "supported_with_caveat"),
    Note = c(
      "Supported as a synthesis layer over validated anchor-review, drift, and equating-chain objects.",
      if (gpcm_detected) {
        "Supported with caveat: bounded GPCM source objects are summarized as exploratory anchor/drift/chain evidence, not an operational linking decision."
      } else {
        "Supported with caveat when bounded GPCM source objects are supplied; not active for this RSM/PCM review."
      }
    )
  )
}

.linking_review_empty <- function() {
  tibble::tibble(
    RiskID = character(),
    Area = character(),
    SourceFamily = character(),
    SourceTable = character(),
    SourceRowKey = character(),
    AdministrationID = character(),
    WaveID = character(),
    LinkKey = character(),
    Facet = character(),
    Level = character(),
    Wave = character(),
    Link = character(),
    Signal = character(),
    Magnitude = numeric(),
    SeverityGroup = character(),
    ReviewPriority = numeric(),
    Guidance = character(),
    PrimaryPlotRoute = character(),
    SupportStatus = character()
  )
}

.linking_review_standardize <- function(tbl) {
  template <- .linking_review_empty()
  if (is.null(tbl) || nrow(tbl) == 0L) {
    return(template)
  }
  for (nm in names(template)) {
    if (!nm %in% names(tbl)) {
      proto <- template[[nm]]
      tbl[[nm]] <- if (is.character(proto)) {
        rep(NA_character_, nrow(tbl))
      } else if (is.numeric(proto)) {
        rep(NA_real_, nrow(tbl))
      } else {
        rep(NA_character_, nrow(tbl))
      }
    }
  }
  tibble::as_tibble(tbl[, names(template), drop = FALSE])
}

.linking_review_group_view_empty <- function() {
  tibble::tibble(
    GroupType = character(),
    GroupKey = character(),
    GroupLabel = character(),
    Facet = character(),
    Cases = integer(),
    MaxPriority = numeric(),
    MeanPriority = numeric(),
    DominantSourceFamily = character(),
    TopRiskID = character(),
    PlotRoutes = character()
  )
}

.linking_review_group_view_index_empty <- function() {
  tibble::tibble(
    View = character(),
    Rows = integer(),
    Description = character()
  )
}

.linking_review_group_summary <- function(tbl, group_type, group_key, group_label, facet = NULL) {
  tbl <- .linking_review_standardize(tbl)
  if (nrow(tbl) == 0L) {
    return(.linking_review_group_view_empty())
  }
  group_key_quo <- rlang::enquo(group_key)
  group_label_quo <- rlang::enquo(group_label)
  facet_quo <- rlang::enquo(facet)

  out <- tbl |>
    dplyr::arrange(.data$SeverityGroup, dplyr::desc(.data$ReviewPriority), dplyr::desc(.data$Magnitude), .data$Area, .data$SourceFamily) |>
    dplyr::mutate(
      .GroupType = as.character(group_type),
      .GroupKey = as.character(!!group_key_quo),
      .GroupLabel = as.character(!!group_label_quo),
      .FacetGroup = if (rlang::quo_is_null(facet_quo)) NA_character_ else as.character(!!facet_quo)
    ) |>
    dplyr::filter(!is.na(.data$.GroupKey), nzchar(.data$.GroupKey)) |>
    dplyr::group_by(.data$.GroupType, .data$.GroupKey, .data$.GroupLabel, .data$.FacetGroup) |>
    dplyr::summarize(
      Cases = dplyr::n(),
      MaxPriority = max(.data$ReviewPriority, na.rm = TRUE),
      MeanPriority = mean(.data$ReviewPriority, na.rm = TRUE),
      DominantSourceFamily = dplyr::first(.data$SourceFamily),
      TopRiskID = dplyr::first(.data$RiskID),
      PlotRoutes = paste(sort(unique(.data$PrimaryPlotRoute)), collapse = " | "),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      GroupType = .data$.GroupType,
      GroupKey = .data$.GroupKey,
      GroupLabel = .data$.GroupLabel,
      Facet = .data$.FacetGroup,
      Cases = as.integer(.data$Cases),
      MaxPriority = .data$MaxPriority,
      MeanPriority = .data$MeanPriority,
      DominantSourceFamily = .data$DominantSourceFamily,
      TopRiskID = .data$TopRiskID,
      PlotRoutes = .data$PlotRoutes
    ) |>
    dplyr::arrange(dplyr::desc(.data$MaxPriority), dplyr::desc(.data$Cases), .data$GroupLabel)

  .linking_review_group_view_empty() |>
    dplyr::slice(0) |>
    dplyr::bind_rows(out)
}

.linking_review_build_group_views <- function(all_risks) {
  all_risks <- .linking_review_standardize(all_risks)
  by_wave <- .linking_review_group_summary(
    all_risks |>
      dplyr::filter(!is.na(.data$WaveID), nzchar(.data$WaveID)),
    group_type = "wave",
    group_key = .data$WaveID,
    group_label = paste0("Wave: ", .data$WaveID)
  )
  by_link <- .linking_review_group_summary(
    all_risks |>
      dplyr::filter(!is.na(.data$LinkKey), nzchar(.data$LinkKey)),
    group_type = "link",
    group_key = .data$LinkKey,
    group_label = paste0("Link: ", .data$LinkKey)
  )
  by_facet <- .linking_review_group_summary(
    all_risks |>
      dplyr::filter(!is.na(.data$Facet), nzchar(.data$Facet)),
    group_type = "facet",
    group_key = .data$Facet,
    group_label = paste0("Facet: ", .data$Facet),
    facet = .data$Facet
  )
  by_source_family <- .linking_review_group_summary(
    all_risks |>
      dplyr::filter(!is.na(.data$SourceFamily), nzchar(.data$SourceFamily)),
    group_type = "source_family",
    group_key = .data$SourceFamily,
    group_label = paste0("Source: ", .data$SourceFamily)
  )
  list(
    by_wave = by_wave,
    by_link = by_link,
    by_facet = by_facet,
    by_source_family = by_source_family
  )
}

.linking_review_group_view_index <- function(group_views) {
  descriptions <- c(
    by_wave = "Concentrated linking risks by fitted wave.",
    by_link = "Concentrated linking risks by adjacent screened link.",
    by_facet = "Concentrated linking risks by facet.",
    by_source_family = "Volume and priority by evidence source family."
  )
  out <- lapply(names(descriptions), function(nm) {
    tbl <- tibble::as_tibble(group_views[[nm]] %||% .linking_review_group_view_empty())
    tibble::tibble(
      View = nm,
      Rows = nrow(tbl),
      Description = descriptions[[nm]]
    )
  }) |>
    dplyr::bind_rows()
  .linking_review_group_view_index_empty() |>
    dplyr::slice(0) |>
    dplyr::bind_rows(out)
}

.review_plot_routes <- function(plot_map) {
  tbl <- tibble::as_tibble(plot_map %||% tibble::tibble())
  if (nrow(tbl) == 0L || !"Available" %in% names(tbl)) {
    return(tbl)
  }
  tbl |>
    dplyr::filter(.data$Available %in% TRUE)
}

.linking_review_anchor_risks <- function(anchor_review) {
  if (is.null(anchor_review)) {
    return(tibble::tibble())
  }

  issue_rows <- tibble::tibble()
  issue_tbl <- tibble::as_tibble(anchor_review$issue_counts %||% tibble::tibble())
  if (nrow(issue_tbl) > 0 && all(c("Issue", "N") %in% names(issue_tbl))) {
    issue_rows <- issue_tbl |>
      dplyr::filter(.data$N > 0) |>
      dplyr::transmute(
        RiskID = paste0("anchor_issue:", .data$Issue),
        Area = "pre_fit_anchor_adequacy",
        SourceFamily = "anchor_issue_count",
        SourceTable = "anchor_review$issue_counts",
        SourceRowKey = as.character(.data$Issue),
        AdministrationID = NA_character_,
        WaveID = NA_character_,
        LinkKey = NA_character_,
        Facet = NA_character_,
        Level = as.character(.data$Issue),
        Wave = NA_character_,
        Link = NA_character_,
        Signal = paste0("Anchor issue count: ", .data$Issue),
        Magnitude = as.numeric(.data$N),
        SeverityGroup = "review",
        ReviewPriority = as.numeric(.data$N),
        Guidance = "Revise anchor/group-anchor tables and rerun review_mfrm_anchors().",
        PrimaryPlotRoute = "plot(anchor_review, type = \"issue_counts\")",
        SupportStatus = "supported"
      )
  }

  overlap_rows <- tibble::tibble()
  facet_tbl <- tibble::as_tibble(anchor_review$facet_summary %||% tibble::tibble())
  min_common <- suppressWarnings(as.integer(anchor_review$thresholds$min_common_anchors %||% NA_integer_))
  # Only flag overlap inadequacy when the user actually supplied anchor
  # or group-anchor constraints. A single-wave fit with no linking
  # design would otherwise show all facets as "high severity" risks
  # because `OverlapLevels == 0` everywhere.
  has_linking_constraint <- nrow(facet_tbl) > 0 &&
    any(c("AnchoredLevels", "GroupedLevels", "ConstrainedLevels") %in% names(facet_tbl)) &&
    sum(
      suppressWarnings(as.numeric(facet_tbl$AnchoredLevels %||% 0)),
      suppressWarnings(as.numeric(facet_tbl$GroupedLevels %||% 0)),
      suppressWarnings(as.numeric(facet_tbl$ConstrainedLevels %||% 0)),
      na.rm = TRUE
    ) > 0
  if (has_linking_constraint &&
      nrow(facet_tbl) > 0 &&
      is.finite(min_common) &&
      all(c("Facet", "OverlapLevels") %in% names(facet_tbl))) {
    overlap_rows <- facet_tbl |>
      dplyr::filter(.data$OverlapLevels < min_common) |>
      dplyr::transmute(
        RiskID = paste0("anchor_overlap:", .data$Facet),
        Area = "pre_fit_anchor_adequacy",
        SourceFamily = "anchor_overlap_support",
        SourceTable = "anchor_review$facet_summary",
        SourceRowKey = as.character(.data$Facet),
        AdministrationID = NA_character_,
        WaveID = NA_character_,
        LinkKey = NA_character_,
        Facet = as.character(.data$Facet),
        Level = NA_character_,
        Wave = NA_character_,
        Link = NA_character_,
        Signal = paste0("Overlap levels below linking guideline (< ", min_common, ")"),
        Magnitude = as.numeric(min_common - .data$OverlapLevels),
        SeverityGroup = "high",
        ReviewPriority = as.numeric(min_common - .data$OverlapLevels),
        Guidance = "Increase common anchor coverage before relying on drift or chain review.",
        PrimaryPlotRoute = "plot(anchor_review, type = \"facet_constraints\")",
        SupportStatus = "supported"
      )
  }

  .linking_review_standardize(dplyr::bind_rows(issue_rows, overlap_rows))
}

.linking_review_drift_risks <- function(drift) {
  if (is.null(drift)) {
    return(tibble::tibble())
  }

  flagged_rows <- tibble::tibble()
  drift_summary <- tibble::as_tibble(drift$summary %||% tibble::tibble())
  drift_threshold <- as.numeric(drift$config$drift_threshold %||% NA_real_)
  if (nrow(drift_summary) > 0 &&
      all(c("Facet", "Wave", "N_Flagged", "Max_Drift") %in% names(drift_summary))) {
    flagged_rows <- drift_summary |>
      dplyr::filter(.data$N_Flagged > 0 | .data$Max_Drift > drift_threshold) |>
      dplyr::transmute(
        RiskID = paste0("drift:", .data$Facet, "::", .data$Wave),
        Area = "post_fit_element_drift",
        SourceFamily = "anchor_drift",
        SourceTable = "drift$summary",
        SourceRowKey = paste(.data$Facet, .data$Wave, sep = "::"),
        AdministrationID = NA_character_,
        WaveID = as.character(.data$Wave),
        LinkKey = NA_character_,
        Facet = as.character(.data$Facet),
        Level = NA_character_,
        Wave = as.character(.data$Wave),
        Link = NA_character_,
        Signal = "Flagged element drift in fitted-wave comparison",
        Magnitude = as.numeric(.data$Max_Drift),
        SeverityGroup = "review",
        ReviewPriority = pmax(as.numeric(.data$N_Flagged), as.numeric(.data$Max_Drift), na.rm = TRUE),
        Guidance = "Inspect flagged rows in detect_anchor_drift() and review wave-specific drift plots.",
        PrimaryPlotRoute = "plot_anchor_drift(drift, type = \"drift\")",
        SupportStatus = "supported"
      )
  }

  support_rows <- tibble::tibble()
  common_tbl <- tibble::as_tibble(drift$common_by_facet %||% tibble::tibble())
  if (nrow(common_tbl) > 0 &&
      all(c("Facet", "Wave", "N_Retained", "GuidelineMinCommon", "LinkSupportAdequate") %in% names(common_tbl))) {
    support_rows <- common_tbl |>
      dplyr::filter(!.data$LinkSupportAdequate) |>
      dplyr::transmute(
        RiskID = paste0("thin_link:", .data$Facet, "::", .data$Wave),
        Area = "post_fit_element_drift",
        SourceFamily = "thin_link_support",
        SourceTable = "drift$common_by_facet",
        SourceRowKey = paste(.data$Facet, .data$Wave, sep = "::"),
        AdministrationID = NA_character_,
        WaveID = as.character(.data$Wave),
        LinkKey = NA_character_,
        Facet = as.character(.data$Facet),
        Level = NA_character_,
        Wave = as.character(.data$Wave),
        Link = NA_character_,
        Signal = "Retained common-element support is below the package guideline",
        Magnitude = as.numeric(.data$GuidelineMinCommon - .data$N_Retained),
        SeverityGroup = "high",
        ReviewPriority = as.numeric(.data$GuidelineMinCommon - .data$N_Retained),
        Guidance = "Treat drift flags as low-support until more retained common elements are available.",
        PrimaryPlotRoute = "plot_anchor_drift(drift, type = \"drift\")",
        SupportStatus = "supported"
      )
  }

  .linking_review_standardize(dplyr::bind_rows(flagged_rows, support_rows))
}

.linking_review_chain_risks <- function(chain) {
  if (is.null(chain)) {
    return(tibble::tibble())
  }

  links_tbl <- tibble::as_tibble(chain$links %||% tibble::tibble())
  if (nrow(links_tbl) == 0 ||
      !all(c("From", "To", "Max_Residual", "LinkSupportAdequate") %in% names(links_tbl))) {
    return(tibble::tibble())
  }

  drift_threshold <- as.numeric(chain$config$drift_threshold %||% NA_real_)
  links_tbl |>
    dplyr::mutate(
      LinkLabel = paste0(.data$From, " -> ", .data$To)
    ) |>
    dplyr::filter(!.data$LinkSupportAdequate | .data$Max_Residual > drift_threshold) |>
    dplyr::transmute(
      RiskID = paste0(ifelse(.data$LinkSupportAdequate, "chain_link:", "chain_support:"), .data$LinkLabel),
      Area = "chain_level_stability",
      SourceFamily = ifelse(.data$LinkSupportAdequate, "equating_chain_link", "equating_chain_support"),
      SourceTable = "chain$links",
      SourceRowKey = as.character(.data$LinkLabel),
      AdministrationID = NA_character_,
      WaveID = NA_character_,
      LinkKey = as.character(.data$LinkLabel),
      Facet = NA_character_,
      Level = NA_character_,
      Wave = NA_character_,
      Link = as.character(.data$LinkLabel),
      Signal = ifelse(
        .data$LinkSupportAdequate,
        "Large residual spread in an adjacent screened link",
        "Thin retained support in an adjacent screened link"
      ),
      Magnitude = as.numeric(.data$Max_Residual),
      SeverityGroup = ifelse(.data$LinkSupportAdequate, "review", "high"),
      ReviewPriority = ifelse(
        .data$LinkSupportAdequate,
        as.numeric(.data$Max_Residual),
        as.numeric(pmax(.data$Min_Common_Per_Facet - .data$Min_Retained_Per_Facet, 1))
      ),
      Guidance = "Inspect the adjacent link and cumulative offsets before using the chain for operational review.",
      PrimaryPlotRoute = "plot_anchor_drift(chain, type = \"chain\")",
      SupportStatus = "supported"
    ) |>
    .linking_review_standardize()
}

#' Build a linking-review synthesis object
#'
#' @param anchor_review Optional output from [review_mfrm_anchors()].
#' @param drift Optional output from [detect_anchor_drift()].
#' @param chain Optional output from [build_equating_chain()].
#' @param top_n Maximum number of linking-risk rows to highlight in summary
#'   outputs. The full object keeps the full risk tables.
#'
#' @details
#' `build_linking_review()` does not recompute anchor, drift, or chain
#' statistics. It is a synthesis layer that organizes package-native evidence
#' into one operational review surface with:
#'
#' - a front-door status block,
#' - ranked linking risks,
#' - explicit next actions,
#' - plot routing metadata,
#' - a reporting/export handoff map.
#'
#' The helper keeps the current conservative interpretation policy:
#' anchor drift and screened links are operational review tools, not automatic
#' proofs of scale equivalence or score comparability.
#'
#' @section Recommended input route:
#' Use existing package-native outputs in this order:
#' 1. [review_mfrm_anchors()] for pre-fit anchor adequacy.
#' 2. [detect_anchor_drift()] for direct wave-to-reference drift screening.
#' 3. [build_equating_chain()] for adjacent screened-link review across waves.
#'
#' @section Interpreting output:
#' - `overview`: which evidence sources were supplied and the current review status.
#' - `top_linking_risks`: primary operational triage table.
#' - `group_view_index`: stable wave/link/facet/source-family grouping routes.
#' - `plot_map`: which existing plotting helper should be used next.
#' - `reporting_map`: what is covered here versus which manuscript-oriented
#'   helper should be used separately.
#'
#' @section GPCM boundary:
#' This helper is currently intended for the validated `RSM` / `PCM` linking
#' workflow. If the supplied drift/chain sources resolve to bounded `GPCM`,
#' the helper stops with a package-level message rather than silently implying
#' support.
#'
#' @return An object of class `mfrm_linking_review`.
#' @seealso [review_mfrm_anchors()], [detect_anchor_drift()],
#'   [build_equating_chain()], [plot_anchor_drift()], [mfrmr_linking_and_dff]
#' @examples
#' \donttest{
#' d1 <- load_mfrmr_data("study1")
#' d2 <- load_mfrmr_data("study2")
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' anchor_review_obj <- review_mfrm_anchors(d1, "Person", c("Rater", "Criterion"), "Score")
#' drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
#' chain <- build_equating_chain(list(Wave1 = fit1, Wave2 = fit2))
#' review <- build_linking_review(anchor_review = anchor_review_obj, drift = drift, chain = chain)
#' summary(review)
#' review$top_linking_risks
#' review$group_view_index
#' }
#' @export
build_linking_review <- function(anchor_review = NULL,
                                 drift = NULL,
                                 chain = NULL,
                                 top_n = 10) {
  anchor_review <- .validate_linking_review_input(
    anchor_review,
    "anchor_review",
    "mfrm_anchor_review"
  )
  drift <- .validate_linking_review_input(drift, "drift", "mfrm_anchor_drift")
  chain <- .validate_linking_review_input(chain, "chain", "mfrm_equating_chain")
  top_n <- max(1L, as.integer(top_n %||% 10L))

  if (is.null(anchor_review) && is.null(drift) && is.null(chain)) {
    stop(
      "build_linking_review() requires at least one of `anchor_review`, `drift`, or `chain`.",
      call. = FALSE
    )
  }

  source_models <- unique(stats::na.omit(c(
    as.character(drift$config$models %||% character(0)),
    as.character(chain$config$models %||% character(0))
  )))
  gpcm_detected <- any(source_models == "GPCM")

  anchor_risks <- .linking_review_anchor_risks(anchor_review)
  drift_risks <- .linking_review_drift_risks(drift)
  chain_risks <- .linking_review_chain_risks(chain)
  all_risks <- .linking_review_standardize(dplyr::bind_rows(anchor_risks, drift_risks, chain_risks))
  if (nrow(all_risks) > 0) {
    severity_levels <- c("high", "review")
    all_risks <- all_risks |>
      dplyr::mutate(
        SeverityGroup = factor(.data$SeverityGroup, levels = severity_levels, ordered = TRUE)
      ) |>
      dplyr::arrange(.data$SeverityGroup, dplyr::desc(.data$ReviewPriority), dplyr::desc(.data$Magnitude), .data$Area, .data$SourceFamily) |>
      dplyr::mutate(
        RiskRank = dplyr::row_number()
      ) |>
      dplyr::mutate(
        SeverityGroup = as.character(.data$SeverityGroup)
      )
  }
  group_views <- .linking_review_build_group_views(all_risks)
  group_view_index <- .linking_review_group_view_index(group_views)

  insufficient_support <- any(all_risks$SourceFamily %in% c(
    "anchor_overlap_support",
    "thin_link_support",
    "equating_chain_support"
  ))
  has_review_risks <- nrow(all_risks) > 0

  review_status <- dplyr::case_when(
    insufficient_support ~ "insufficient_anchor_evidence",
    has_review_risks ~ "review_required",
    gpcm_detected ~ "exploratory_gpcm_review",
    TRUE ~ "stable_for_linking_review"
  )

  support_status <- .linking_review_support_status(source_models)
  overview <- tibble::tibble(
    AnchorReviewAvailable = !is.null(anchor_review),
    DriftAvailable = !is.null(drift),
    ChainAvailable = !is.null(chain),
    ReviewStatus = review_status,
    TopRiskRows = nrow(all_risks),
    GroupViews = sum(group_view_index$Rows > 0L),
    SourceModels = if (length(source_models) > 0) paste(source_models, collapse = ", ") else NA_character_,
    GPCMSupport = as.character(support_status$Status[support_status$Scope == "bounded GPCM"][1] %||% NA_character_)
  )

  key_warnings <- character(0)
  if (nrow(anchor_risks) > 0) {
    key_warnings <- c(
      key_warnings,
      paste0("Anchor review flagged ", nrow(anchor_risks), " pre-fit anchor risk rows.")
    )
  }
  if (nrow(drift_risks) > 0) {
    key_warnings <- c(
      key_warnings,
      paste0("Drift review flagged ", nrow(drift_risks), " wave/facet support or drift rows.")
    )
  }
  if (nrow(chain_risks) > 0) {
    key_warnings <- c(
      key_warnings,
      paste0("Chain review flagged ", nrow(chain_risks), " adjacent-link instability rows.")
    )
  }
  key_warnings <- clean_summary_lines(key_warnings, max_n = 4L)
  if (length(key_warnings) == 0) {
    key_warnings <- "No immediate linking-review warnings."
  }

  next_actions <- character(0)
  if (nrow(anchor_risks) > 0) {
    next_actions <- c(
      next_actions,
      "Revise anchors with make_anchor_table() / review_mfrm_anchors() before relying on drift or chain review."
    )
  }
  if (nrow(drift_risks) > 0) {
    next_actions <- c(
      next_actions,
      "Inspect detect_anchor_drift() and plot_anchor_drift(drift, type = \"drift\") for wave-level follow-up."
    )
  }
  if (nrow(chain_risks) > 0) {
    next_actions <- c(
      next_actions,
      "Inspect build_equating_chain() and plot_anchor_drift(chain, type = \"chain\") before using cumulative offsets operationally."
    )
  }
  if (length(next_actions) == 0) {
    next_actions <- c(
      next_actions,
      "After confirming connectivity and anchor support, continue with anchor_to_baseline() or DFF follow-up as needed."
    )
  }
  next_actions <- clean_summary_lines(next_actions, max_n = 4L)

  status <- make_summary_block(
    "Overall status" = review_status,
    "Evidence sources" = paste(c(
      if (!is.null(anchor_review)) "anchor_review",
      if (!is.null(drift)) "drift",
      if (!is.null(chain)) "chain"
    ), collapse = ", "),
    "Bounded GPCM" = as.character(support_status$Status[support_status$Scope == "bounded GPCM"][1] %||% NA_character_)
  )

  plot_map <- tibble::tibble(
    ReviewArea = c("Anchor adequacy", "Wave-level drift", "Screened chain"),
    Available = c(!is.null(anchor_review), !is.null(drift), !is.null(chain)),
    PlotHelper = c(
      "plot(anchor_review, type = \"issue_counts\")",
      "plot_anchor_drift(drift, type = \"drift\")",
      "plot_anchor_drift(chain, type = \"chain\")"
    ),
    Trigger = c(
      "Use when anchor issues or overlap warnings are present.",
      "Use when fitted waves show flagged drift or thin retained common-element support.",
      "Use when adjacent links show thin support or large residual spread."
    )
  )

  reporting_map <- tibble::tibble(
    Area = c(
      "Operational linking review",
      "Pre-fit anchor repair",
      "Wave-level drift detail",
      "Chain-level cumulative offsets",
      "Manuscript reporting / appendix"
    ),
    CoveredHere = c("yes", "partial", "partial", "partial", "partial"),
    CompanionOutput = c(
      "summary(build_linking_review(...))",
      "review_mfrm_anchors() / make_anchor_table()",
      "detect_anchor_drift() / plot_anchor_drift()",
      "build_equating_chain() / plot_anchor_drift(type = \"chain\")",
      "build_summary_table_bundle(summary(linking_review)) / reporting_checklist()"
    )
  )

  notes <- clean_summary_lines(c(
    if (gpcm_detected) {
      "Linking review is an exploratory synthesis layer over existing package-native bounded GPCM anchor, drift, and chain evidence."
    } else {
      "Linking review is an operational synthesis layer over existing package-native anchor, drift, and chain evidence."
    },
    "Drift or thin-support warnings do not prove scale breakdown by themselves; they indicate where review is needed.",
    "Repeated signals across anchor, drift, and chain evidence deserve priority, but this helper does not collapse them into one opaque composite score.",
    if (gpcm_detected) {
      paste(
        "Bounded GPCM linking review is exploratory: it indexes direct",
        "anchor/drift/chain evidence and should not be reported as an",
        "operational GPCM linking decision or as evidence that drift is absent."
      )
    } else {
      ""
    }
  ))

  out <- list(
    overview = overview,
    status = status,
    key_warnings = key_warnings,
    next_actions = next_actions,
    top_linking_risks = tibble::as_tibble(all_risks),
    group_view_index = group_view_index,
    group_views = group_views,
    prefit_anchor_risks = tibble::as_tibble(anchor_risks),
    drift_risks = tibble::as_tibble(drift_risks),
    chain_risks = tibble::as_tibble(chain_risks),
    plot_map = plot_map,
    reporting_map = reporting_map,
    support_status = support_status,
    gpcm_boundary = gpcm_capability_boundary_table(
      fit = if (gpcm_detected) {
        structure(
          list(config = list(model = "GPCM"), summary = data.frame(Model = "GPCM")),
          class = "mfrm_fit"
        )
      } else {
        NULL
      },
      helper = "build_linking_review()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review"
      )
    ),
    notes = notes,
    settings = list(
      top_n = top_n,
      intended_use = if (gpcm_detected) "exploratory_gpcm_linking_review" else "operational_linking_review",
      source_models = source_models,
      source_profile = data.frame(
        Source = c("anchor_review", "drift", "chain"),
        Available = c(!is.null(anchor_review), !is.null(drift), !is.null(chain)),
        Class = c(
          class(anchor_review)[1] %||% NA_character_,
          class(drift)[1] %||% NA_character_,
          class(chain)[1] %||% NA_character_
        ),
        stringsAsFactors = FALSE
      )
    )
  )
  as_mfrm_bundle(out, "mfrm_linking_review")
}

#' @export
print.mfrm_linking_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Summarize a linking-review object
#'
#' @param object Output from [build_linking_review()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of top linking-risk rows to keep in the compact summary.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_linking_review`.
#' @seealso [build_linking_review()]
#' @export
summary.mfrm_linking_review <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_linking_review")) {
    stop("`object` must be output from build_linking_review().", call. = FALSE)
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  top_risks <- tibble::as_tibble(object$top_linking_risks %||% tibble::tibble())
  if (nrow(top_risks) > 0) {
    top_risks <- top_risks |>
      dplyr::slice_head(n = top_n)
  }

  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    status = tibble::as_tibble(object$status %||% tibble::tibble()),
    key_warnings = clean_summary_lines(object$key_warnings %||% character(0), max_n = 4L),
    next_actions = clean_summary_lines(object$next_actions %||% character(0), max_n = 4L),
    top_linking_risks = top_risks,
    group_view_index = tibble::as_tibble(object$group_view_index %||% tibble::tibble()),
    group_views = object$group_views %||% list(),
    prefit_anchor_risks = tibble::as_tibble(object$prefit_anchor_risks %||% tibble::tibble()),
    drift_risks = tibble::as_tibble(object$drift_risks %||% tibble::tibble()),
    chain_risks = tibble::as_tibble(object$chain_risks %||% tibble::tibble()),
    plot_routes = .review_plot_routes(object$plot_map),
    plot_map = tibble::as_tibble(object$plot_map %||% tibble::tibble()),
    reporting_map = tibble::as_tibble(object$reporting_map %||% tibble::tibble()),
    support_status = tibble::as_tibble(object$support_status %||% tibble::tibble()),
    gpcm_boundary = tibble::as_tibble(object$gpcm_boundary %||% tibble::tibble()),
    notes = clean_summary_lines(object$notes %||% character(0)),
    settings = object$settings %||% list(),
    digits = digits
  )
  class(out) <- "summary.mfrm_linking_review"
  out
}

#' @export
print.summary.mfrm_linking_review <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Linking Review Summary\n")
  if (nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$status) > 0) {
    cat("\nStatus\n")
    print(as.data.frame(x$status), row.names = FALSE)
  }
  print_bullet_section("Key Warnings", x$key_warnings)
  print_bullet_section("Next Actions", x$next_actions)
  if (nrow(x$top_linking_risks) > 0) {
    cat("\nTop Linking Risks\n")
    print(round_numeric_df(as.data.frame(x$top_linking_risks), digits = digits), row.names = FALSE)
  }
  if (nrow(x$group_view_index) > 0) {
    cat("\nGrouping Views\n")
    print(as.data.frame(x$group_view_index), row.names = FALSE)
  }
  if (nrow(x$plot_routes) > 0) {
    cat("\nPlot Follow-up\n")
    print(as.data.frame(x$plot_routes), row.names = FALSE)
  }
  if (nrow(x$support_status) > 0) {
    cat("\nSupport Status\n")
    print(as.data.frame(x$support_status), row.names = FALSE)
  }
  if (nrow(x$gpcm_boundary) > 0) {
    cat("\nGPCM Boundary\n")
    print(as.data.frame(x$gpcm_boundary)[, c("Area", "Status"), drop = FALSE], row.names = FALSE)
  }
  print_bullet_section("Notes", x$notes)
  invisible(x)
}

# --- build_misfit_casebook ---------------------------------------------------

.validate_misfit_casebook_input <- function(x, arg, classes) {
  if (is.null(x)) return(NULL)
  if (!inherits(x, classes)) {
    stop(
      "`", arg, "` must be one of: ", paste(classes, collapse = ", "), ".",
      call. = FALSE
    )
  }
  x
}

.misfit_casebook_support_status <- function(model) {
  model <- as.character(model %||% NA_character_)[1]
  gpcm <- identical(model, "GPCM")
  tibble::tibble(
    Scope = c("RSM / PCM", "bounded GPCM"),
    Status = c(
      "supported",
      if (gpcm) "supported_with_caveat" else "deferred"
    ),
    Note = c(
      "Supported as a synthesis layer over package-native screening outputs.",
      if (gpcm) {
        paste(
          "Supported with caveat: bounded GPCM casebook rows inherit",
          "exploratory screening semantics from residual and strict marginal sources."
        )
      } else {
        "Deferred unless a bounded GPCM source fit is supplied."
      }
    )
  )
}

.misfit_casebook_source_support <- function(model, diagnostics, unexpected, displacement) {
  model <- as.character(model %||% NA_character_)[1]
  gpcm <- identical(model, "GPCM")
  marginal_fit <- diagnostics$marginal_fit %||% list()
  pairwise_fit <- marginal_fit$pairwise %||% list()

  rows <- list(
    list(
      SourceFamily = "marginal_cell",
      Available = isTRUE(marginal_fit$available),
      SupportBasis = "marginal_fit",
      Status = if (isTRUE(marginal_fit$available)) {
        if (gpcm) "supported_with_caveat" else "supported"
      } else {
        "deferred"
      },
      Note = if (isTRUE(marginal_fit$available)) {
        if (gpcm) {
          paste(
            "Strict marginal cell screening is available, but bounded GPCM rows remain",
            "exploratory and should not be treated as formal item-fit tests."
          )
        } else {
          "Strict marginal cell screening is available for operational follow-up."
        }
      } else {
        as.character(
          marginal_fit$summary$Reason[1] %||%
            marginal_fit$notes[1] %||%
            "Strict marginal cell screening is unavailable for the current run."
        )
      }
    ),
    list(
      SourceFamily = "marginal_pair",
      Available = isTRUE(pairwise_fit$available),
      SupportBasis = "marginal_fit",
      Status = if (isTRUE(pairwise_fit$available)) {
        if (gpcm) "supported_with_caveat" else "supported"
      } else {
        "deferred"
      },
      Note = if (isTRUE(pairwise_fit$available)) {
        if (gpcm) {
          paste(
            "Strict pairwise screening is available, but bounded GPCM rows remain",
            "exploratory and should not be treated as formal local-dependence tests."
          )
        } else {
          "Strict pairwise screening is available for operational follow-up."
        }
      } else {
        as.character(
          pairwise_fit$summary$Reason[1] %||%
            pairwise_fit$notes[1] %||%
            "Strict pairwise screening is unavailable for the current run."
        )
      }
    ),
    list(
      SourceFamily = "unexpected",
      Available = inherits(unexpected, "mfrm_unexpected"),
      SupportBasis = "legacy",
      Status = if (inherits(unexpected, "mfrm_unexpected")) {
        if (gpcm) "supported_with_caveat" else "supported"
      } else {
        "deferred"
      },
      Note = if (inherits(unexpected, "mfrm_unexpected")) {
        if (gpcm) {
          "Unexpected-response rows are available, but bounded GPCM interpretation remains operational rather than inferential."
        } else {
          "Unexpected-response rows are available for operational follow-up."
        }
      } else {
        "Unexpected-response screening was not available for the current run."
      }
    ),
    list(
      SourceFamily = "displacement",
      Available = inherits(displacement, "mfrm_displacement"),
      SupportBasis = "legacy",
      Status = if (inherits(displacement, "mfrm_displacement")) {
        if (gpcm) "supported_with_caveat" else "supported"
      } else {
        "deferred"
      },
      Note = if (inherits(displacement, "mfrm_displacement")) {
        if (gpcm) {
          "Displacement rows are available, but bounded GPCM interpretation remains operational rather than inferential."
        } else {
          "Displacement rows are available for operational follow-up."
        }
      } else {
        "Displacement screening was not available for the current run."
      }
    )
  )

  tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))) |>
    dplyr::mutate(
      Available = as.logical(.data$Available),
      SourceFamily = as.character(.data$SourceFamily),
      SupportBasis = as.character(.data$SupportBasis),
      Status = as.character(.data$Status),
      Note = as.character(.data$Note)
    )
}

.misfit_casebook_context_key <- function(tbl, cols) {
  if (is.null(tbl) || nrow(tbl) == 0L || length(cols) == 0L) {
    return(rep(NA_character_, nrow(tbl %||% data.frame())))
  }
  cols <- cols[cols %in% names(tbl)]
  if (length(cols) == 0L) {
    return(rep(NA_character_, nrow(tbl)))
  }
  apply(as.data.frame(tbl[, cols, drop = FALSE]), 1, function(x) {
    paste(as.character(x), collapse = " | ")
  })
}

.misfit_casebook_empty <- function() {
  tibble::tibble(
    CaseID = character(),
    CaseType = character(),
    SourceFamily = character(),
    SourceTable = character(),
    SourceRowKey = character(),
    AdministrationID = character(),
    WaveID = character(),
    PrimaryUnit = character(),
    PrimaryUnitType = character(),
    Person = character(),
    Facet = character(),
    Level = character(),
    Category = integer(),
    PairKey = character(),
    ContextKey = character(),
    Wave = character(),
    Signal = character(),
    Direction = character(),
    Magnitude = numeric(),
    ReviewPriority = numeric(),
    WithinSourceRank = integer(),
    EvidenceN = integer(),
    SupportBasis = character(),
    InterpretationTier = character(),
    PrimaryPlotRoute = character(),
    SupportStatus = character()
  )
}

.misfit_casebook_standardize <- function(tbl) {
  template <- .misfit_casebook_empty()
  if (is.null(tbl) || nrow(tbl) == 0L) {
    return(template)
  }
  for (nm in names(template)) {
    if (!nm %in% names(tbl)) {
      proto <- template[[nm]]
      tbl[[nm]] <- if (is.character(proto)) {
        rep(NA_character_, nrow(tbl))
      } else if (is.integer(proto)) {
        rep(NA_integer_, nrow(tbl))
      } else if (is.numeric(proto)) {
        rep(NA_real_, nrow(tbl))
      } else if (is.logical(proto)) {
        rep(NA, nrow(tbl))
      } else {
        rep(NA_character_, nrow(tbl))
      }
    }
  }
  tbl <- tbl[, names(template), drop = FALSE]
  tibble::as_tibble(tbl)
}

.misfit_casebook_rollup_empty <- function() {
  tibble::tibble(
    AdministrationID = character(),
    WaveID = character(),
    RollupType = character(),
    RollupKey = character(),
    RollupLabel = character(),
    SourceFamily = character(),
    Facet = character(),
    SupportBasis = character(),
    InterpretationTier = character(),
    PrimaryPlotRoute = character(),
    SupportStatus = character(),
    Cases = integer(),
    DistinctSourceRows = integer(),
    PersonsFlagged = integer(),
    MaxPriority = numeric(),
    MeanPriority = numeric(),
    EvidenceN = integer(),
    TopCaseID = character()
  )
}

.misfit_casebook_group_view_empty <- function() {
  tibble::tibble(
    AdministrationID = character(),
    WaveID = character(),
    GroupType = character(),
    GroupKey = character(),
    GroupLabel = character(),
    Facet = character(),
    Cases = integer(),
    DistinctSourceRows = integer(),
    PersonsFlagged = integer(),
    MaxPriority = numeric(),
    MeanPriority = numeric(),
    EvidenceN = integer(),
    DominantSourceFamily = character(),
    TopCaseID = character(),
    PlotRoutes = character()
  )
}

.misfit_casebook_group_view_index_empty <- function() {
  tibble::tibble(
    View = character(),
    Rows = integer(),
    Description = character()
  )
}

.normalize_misfit_casebook_id <- function(x, arg) {
  if (is.null(x)) {
    return(NA_character_)
  }
  if (is.list(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be NULL or a single non-missing scalar identifier.", call. = FALSE)
  }
  value <- trimws(as.character(x)[1])
  if (!nzchar(value)) {
    stop("`", arg, "` must be NULL or a non-empty scalar identifier.", call. = FALSE)
  }
  value
}

.misfit_casebook_apply_provenance <- function(tbl,
                                              administration_id = NA_character_,
                                              wave_id = NA_character_) {
  tbl <- .misfit_casebook_standardize(tbl)
  if (nrow(tbl) == 0L) {
    return(tbl)
  }
  tbl$AdministrationID <- rep(as.character(administration_id %||% NA_character_), nrow(tbl))
  tbl$WaveID <- rep(as.character(wave_id %||% NA_character_), nrow(tbl))
  tbl
}

.misfit_casebook_build_rollup <- function(tbl) {
  tbl <- .misfit_casebook_standardize(tbl)
  if (nrow(tbl) == 0L) {
    return(.misfit_casebook_rollup_empty())
  }

  rollup_tbl <- tbl |>
    dplyr::arrange(dplyr::desc(.data$ReviewPriority), .data$SourceFamily, .data$WithinSourceRank) |>
    dplyr::mutate(
      RollupType = dplyr::case_when(
        !is.na(.data$Person) & nzchar(.data$Person) ~ "person",
        !is.na(.data$PairKey) & nzchar(.data$PairKey) ~ "facet_pair",
        !is.na(.data$Facet) & nzchar(.data$Facet) & !is.na(.data$Level) & nzchar(.data$Level) ~ "facet_level",
        TRUE ~ "source_family"
      ),
      RollupKey = dplyr::case_when(
        .data$RollupType == "person" ~ .data$Person,
        .data$RollupType == "facet_pair" ~ paste(.data$Facet, .data$PairKey, sep = "::"),
        .data$RollupType == "facet_level" ~ paste(.data$Facet, .data$Level, sep = "::"),
        TRUE ~ .data$SourceFamily
      ),
      RollupLabel = dplyr::case_when(
        .data$RollupType == "person" ~ paste0("Person: ", .data$Person),
        .data$RollupType == "facet_pair" ~ paste0(.data$Facet, " pair: ", .data$PairKey),
        .data$RollupType == "facet_level" ~ paste0(.data$Facet, ": ", .data$Level),
        TRUE ~ paste0("Source: ", .data$SourceFamily)
      )
    ) |>
    dplyr::group_by(
      .data$AdministrationID,
      .data$WaveID,
      .data$RollupType,
      .data$RollupKey,
      .data$RollupLabel,
      .data$SourceFamily,
      .data$Facet,
      .data$SupportBasis,
      .data$InterpretationTier,
      .data$PrimaryPlotRoute,
      .data$SupportStatus
    ) |>
    dplyr::summarize(
      Cases = dplyr::n(),
      DistinctSourceRows = dplyr::n_distinct(.data$SourceRowKey),
      PersonsFlagged = dplyr::n_distinct(.data$Person[!is.na(.data$Person) & nzchar(.data$Person)]),
      MaxPriority = max(.data$ReviewPriority, na.rm = TRUE),
      MeanPriority = mean(.data$ReviewPriority, na.rm = TRUE),
      EvidenceN = sum(.data$EvidenceN, na.rm = TRUE),
      TopCaseID = dplyr::first(.data$CaseID),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$MaxPriority), dplyr::desc(.data$Cases), .data$SourceFamily, .data$RollupLabel)

  .misfit_casebook_rollup_empty() |>
    dplyr::slice(0) |>
    dplyr::bind_rows(rollup_tbl)
}

.misfit_casebook_group_summary <- function(tbl, group_type, group_key, group_label, facet = NULL) {
  tbl <- .misfit_casebook_standardize(tbl)
  if (nrow(tbl) == 0L) {
    return(.misfit_casebook_group_view_empty())
  }
  group_key_quo <- rlang::enquo(group_key)
  group_label_quo <- rlang::enquo(group_label)
  facet_quo <- rlang::enquo(facet)

  tbl <- tbl |>
    dplyr::arrange(dplyr::desc(.data$ReviewPriority), .data$SourceFamily, .data$WithinSourceRank) |>
    dplyr::mutate(
      .GroupType = as.character(group_type),
      .GroupKey = as.character(!!group_key_quo),
      .GroupLabel = as.character(!!group_label_quo),
      .FacetGroup = if (rlang::quo_is_null(facet_quo)) NA_character_ else as.character(!!facet_quo)
    ) |>
    dplyr::filter(!is.na(.data$.GroupKey), nzchar(.data$.GroupKey))

  if (nrow(tbl) == 0L) {
    return(.misfit_casebook_group_view_empty())
  }

  out <- tbl |>
    dplyr::group_by(
      .data$AdministrationID,
      .data$WaveID,
      .data$.GroupType,
      .data$.GroupKey,
      .data$.GroupLabel,
      .data$.FacetGroup
    ) |>
    dplyr::summarize(
      Cases = dplyr::n(),
      DistinctSourceRows = dplyr::n_distinct(.data$SourceRowKey),
      PersonsFlagged = dplyr::n_distinct(.data$Person[!is.na(.data$Person) & nzchar(.data$Person)]),
      MaxPriority = max(.data$ReviewPriority, na.rm = TRUE),
      MeanPriority = mean(.data$ReviewPriority, na.rm = TRUE),
      EvidenceN = sum(.data$EvidenceN, na.rm = TRUE),
      DominantSourceFamily = dplyr::first(.data$SourceFamily),
      TopCaseID = dplyr::first(.data$CaseID),
      PlotRoutes = paste(sort(unique(.data$PrimaryPlotRoute)), collapse = " | "),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      AdministrationID = .data$AdministrationID,
      WaveID = .data$WaveID,
      GroupType = .data$.GroupType,
      GroupKey = .data$.GroupKey,
      GroupLabel = .data$.GroupLabel,
      Facet = .data$.FacetGroup,
      Cases = as.integer(.data$Cases),
      DistinctSourceRows = as.integer(.data$DistinctSourceRows),
      PersonsFlagged = as.integer(.data$PersonsFlagged),
      MaxPriority = .data$MaxPriority,
      MeanPriority = .data$MeanPriority,
      EvidenceN = as.integer(.data$EvidenceN),
      DominantSourceFamily = .data$DominantSourceFamily,
      TopCaseID = .data$TopCaseID,
      PlotRoutes = .data$PlotRoutes
    ) |>
    dplyr::arrange(dplyr::desc(.data$MaxPriority), dplyr::desc(.data$Cases), .data$GroupLabel)

  .misfit_casebook_group_view_empty() |>
    dplyr::slice(0) |>
    dplyr::bind_rows(out)
}

.misfit_casebook_build_group_views <- function(all_cases, case_rollup) {
  all_cases <- .misfit_casebook_standardize(all_cases)
  case_rollup <- tibble::as_tibble(case_rollup %||% .misfit_casebook_rollup_empty())

  by_person <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$Person), nzchar(.data$Person)),
    group_type = "person",
    group_key = .data$Person,
    group_label = paste0("Person: ", .data$Person)
  )

  by_facet_level <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$Facet), nzchar(.data$Facet), !is.na(.data$Level), nzchar(.data$Level)),
    group_type = "facet_level",
    group_key = paste(.data$Facet, .data$Level, sep = "::"),
    group_label = paste0(.data$Facet, ": ", .data$Level),
    facet = .data$Facet
  )

  by_facet_pair <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$Facet), nzchar(.data$Facet), !is.na(.data$PairKey), nzchar(.data$PairKey)),
    group_type = "facet_pair",
    group_key = paste(.data$Facet, .data$PairKey, sep = "::"),
    group_label = paste0(.data$Facet, " pair: ", .data$PairKey),
    facet = .data$Facet
  )

  by_source_family <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$SourceFamily), nzchar(.data$SourceFamily)),
    group_type = "source_family",
    group_key = .data$SourceFamily,
    group_label = paste0("Source: ", .data$SourceFamily)
  )

  by_facet <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$Facet), nzchar(.data$Facet)),
    group_type = "facet",
    group_key = .data$Facet,
    group_label = paste0("Facet: ", .data$Facet),
    facet = .data$Facet
  )

  by_administration <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$AdministrationID), nzchar(.data$AdministrationID)),
    group_type = "administration",
    group_key = .data$AdministrationID,
    group_label = paste0("Administration: ", .data$AdministrationID)
  )

  by_wave <- .misfit_casebook_group_summary(
    all_cases |>
      dplyr::filter(!is.na(.data$WaveID), nzchar(.data$WaveID)),
    group_type = "wave",
    group_key = .data$WaveID,
    group_label = paste0("Wave: ", .data$WaveID)
  )

  facet_rollup <- case_rollup |>
    dplyr::filter(!is.na(.data$Facet), nzchar(.data$Facet))
  facet_views <- split(facet_rollup, facet_rollup$Facet)
  facet_views <- lapply(facet_views, tibble::as_tibble)

  list(
    by_person = by_person,
    by_facet_level = by_facet_level,
    by_facet_pair = by_facet_pair,
    by_source_family = by_source_family,
    by_facet = by_facet,
    by_administration = by_administration,
    by_wave = by_wave,
    facet_views = facet_views
  )
}

.misfit_casebook_group_view_index <- function(group_views) {
  fixed_names <- c(
    "by_person",
    "by_facet_level",
    "by_facet_pair",
    "by_source_family",
    "by_facet",
    "by_administration",
    "by_wave"
  )
  descriptions <- c(
    by_person = "Repeated signals concentrated on the same person across evidence families.",
    by_facet_level = "Repeated signals concentrated on the same facet level.",
    by_facet_pair = "Repeated pairwise signals within the same facet.",
    by_source_family = "Volume and priority by evidence source family.",
    by_facet = "All flagged evidence grouped by facet.",
    by_administration = "Operational concentration by administration/form when provided.",
    by_wave = "Operational concentration by wave/occasion when provided."
  )

  index_tbl <- lapply(fixed_names, function(nm) {
    tbl <- tibble::as_tibble(group_views[[nm]] %||% .misfit_casebook_group_view_empty())
    tibble::tibble(
      View = nm,
      Rows = nrow(tbl),
      Description = descriptions[[nm]]
    )
  }) |>
    dplyr::bind_rows()

  facet_views <- group_views$facet_views %||% list()
  if (length(facet_views) > 0L) {
    facet_index <- lapply(names(facet_views), function(nm) {
      tibble::tibble(
        View = paste0("facet_views$", nm),
        Rows = nrow(tibble::as_tibble(facet_views[[nm]])),
        Description = paste0("Case-rollup rows restricted to facet `", nm, "`.")
      )
    }) |>
      dplyr::bind_rows()
    index_tbl <- dplyr::bind_rows(index_tbl, facet_index)
  }

  .misfit_casebook_group_view_index_empty() |>
    dplyr::slice(0) |>
    dplyr::bind_rows(index_tbl)
}

.misfit_casebook_from_marginal_cells <- function(diagnostics, support_status = "supported") {
  src <- tibble::as_tibble(diagnostics$marginal_fit$top_cells %||% tibble::tibble())
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  if ("FlaggedAbsZ" %in% names(src)) {
    src <- src |>
      dplyr::filter(isTRUE(.data$FlaggedAbsZ) | (.data$FlaggedAbsZ %in% TRUE))
  }
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  src <- src |>
    dplyr::arrange(dplyr::desc(.data$AbsStdResidual), dplyr::desc(abs(.data$PropDiff)))

  out <- src |>
    dplyr::mutate(
      SourceRowKey = paste(
        .data$CellType,
        dplyr::coalesce(.data$Facet, "<none>"),
        dplyr::coalesce(.data$Level, "<none>"),
        .data$Category,
        dplyr::coalesce(.data$StepFacet, "<none>"),
        sep = "::"
      ),
      PrimaryUnit = dplyr::if_else(
        !is.na(.data$Level) & nzchar(.data$Level),
        .data$Level,
        paste0("Category ", .data$Category)
      ),
      Direction = dplyr::if_else(
        .data$StdResidual >= 0,
        "Higher than expected",
        "Lower than expected"
      ),
      Signal = "Strict marginal cell screen",
      Magnitude = as.numeric(.data$AbsStdResidual),
      ReviewPriority = as.numeric(.data$AbsStdResidual),
      WithinSourceRank = dplyr::row_number(),
      CaseID = paste0("marginal_cell:", .data$SourceRowKey)
    ) |>
    dplyr::transmute(
      CaseID = .data$CaseID,
      CaseType = "marginal_cell_case",
      SourceFamily = "marginal_cell",
      SourceTable = "diagnostics$marginal_fit$top_cells",
      SourceRowKey = .data$SourceRowKey,
      PrimaryUnit = as.character(.data$PrimaryUnit),
      PrimaryUnitType = dplyr::if_else(!is.na(.data$Level) & nzchar(.data$Level), "facet_level", "score_category"),
      Person = NA_character_,
      Facet = as.character(.data$Facet),
      Level = as.character(.data$Level),
      Category = suppressWarnings(as.integer(.data$Category)),
      PairKey = NA_character_,
      ContextKey = dplyr::coalesce(as.character(.data$StepFacet), as.character(.data$CellType)),
      Wave = NA_character_,
      Signal = .data$Signal,
      Direction = .data$Direction,
      Magnitude = .data$Magnitude,
      ReviewPriority = .data$ReviewPriority,
      WithinSourceRank = as.integer(.data$WithinSourceRank),
      EvidenceN = 1L,
      SupportBasis = "marginal_fit",
      InterpretationTier = "screening_only",
      PrimaryPlotRoute = "plot_marginal_fit(diagnostics, draw = FALSE)",
      SupportStatus = support_status
    )
  .misfit_casebook_standardize(out)
}

.misfit_casebook_from_pairwise <- function(diagnostics, support_status = "supported") {
  src <- tibble::as_tibble(diagnostics$marginal_fit$pairwise$top_pairs %||% tibble::tibble())
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  if ("Flagged" %in% names(src)) {
    src <- src |>
      dplyr::filter(isTRUE(.data$Flagged) | (.data$Flagged %in% TRUE))
  }
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  src <- src |>
    dplyr::mutate(
      PairKey = paste(.data$Level1, .data$Level2, sep = "::"),
      AbsResidual = pmax(
        as.numeric(.data$AbsExactStdResidual %||% 0),
        as.numeric(.data$AbsAdjacentStdResidual %||% 0),
        na.rm = TRUE
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$AbsResidual), dplyr::desc(.data$LevelPairCount))

  out <- src |>
    dplyr::mutate(
      SourceRowKey = paste(.data$Facet, .data$PairKey, sep = "::"),
      CaseID = paste0("pairwise:", .data$SourceRowKey),
      Direction = dplyr::case_when(
        abs(.data$ExactStdResidual) >= abs(.data$AdjacentStdResidual) & .data$ExactStdResidual >= 0 ~ "Higher than expected exact agreement",
        abs(.data$ExactStdResidual) >= abs(.data$AdjacentStdResidual) & .data$ExactStdResidual < 0 ~ "Lower than expected exact agreement",
        .data$AdjacentStdResidual >= 0 ~ "Higher than expected adjacent agreement",
        TRUE ~ "Lower than expected adjacent agreement"
      )
    ) |>
    dplyr::transmute(
      CaseID = .data$CaseID,
      CaseType = "pairwise_local_dependence_case",
      SourceFamily = "marginal_pair",
      SourceTable = "diagnostics$marginal_fit$pairwise$top_pairs",
      SourceRowKey = .data$SourceRowKey,
      PrimaryUnit = paste(.data$Level1, "vs", .data$Level2),
      PrimaryUnitType = "facet_pair",
      Person = NA_character_,
      Facet = as.character(.data$Facet),
      Level = NA_character_,
      Category = NA_integer_,
      PairKey = as.character(.data$PairKey),
      ContextKey = as.character(.data$Facet),
      Wave = NA_character_,
      Signal = "Strict pairwise local-dependence screen",
      Direction = .data$Direction,
      Magnitude = as.numeric(.data$AbsResidual),
      ReviewPriority = as.numeric(.data$AbsResidual),
      WithinSourceRank = as.integer(dplyr::row_number()),
      EvidenceN = 1L,
      SupportBasis = "marginal_fit",
      InterpretationTier = as.character(.data$InferenceTier %||% "screening_only"),
      PrimaryPlotRoute = "plot_marginal_pairwise(diagnostics, draw = FALSE)",
      SupportStatus = support_status
    )
  .misfit_casebook_standardize(out)
}

.misfit_casebook_from_unexpected <- function(unexpected, facet_names, support_status = "supported") {
  src <- tibble::as_tibble(unexpected$table %||% tibble::tibble())
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  src <- src |>
    dplyr::arrange(dplyr::desc(.data$Severity), dplyr::desc(abs(.data$StdResidual)), .data$ObsProb)
  context_key <- .misfit_casebook_context_key(src, facet_names)

  out <- src |>
    dplyr::mutate(
      SourceRowKey = as.character(.data$Row),
      CaseID = paste0("unexpected:", .data$SourceRowKey),
      ContextKey = context_key
    ) |>
    dplyr::transmute(
      CaseID = .data$CaseID,
      CaseType = "unexpected_response_case",
      SourceFamily = "unexpected",
      SourceTable = "unexpected_response_table()",
      SourceRowKey = .data$SourceRowKey,
      PrimaryUnit = as.character(.data$Person %||% .data$Row),
      PrimaryUnitType = "person_observation",
      Person = as.character(.data$Person %||% NA_character_),
      Facet = NA_character_,
      Level = NA_character_,
      Category = suppressWarnings(as.integer(.data$Score)),
      PairKey = NA_character_,
      ContextKey = as.character(.data$ContextKey),
      Wave = NA_character_,
      Signal = "Unexpected response screen",
      Direction = as.character(.data$Direction),
      Magnitude = as.numeric(.data$Severity),
      ReviewPriority = as.numeric(.data$Severity),
      WithinSourceRank = as.integer(dplyr::row_number()),
      EvidenceN = 1L,
      SupportBasis = "legacy",
      InterpretationTier = "operational_review",
      PrimaryPlotRoute = "plot_unexpected(unexpected, draw = FALSE)",
      SupportStatus = support_status
    )
  .misfit_casebook_standardize(out)
}

.misfit_casebook_from_element_fit <- function(diagnostics,
                                              support_status = "supported",
                                              lower = NULL,
                                              upper = NULL) {
  band <- mfrm_misfit_thresholds(lower = lower, upper = upper)
  lower <- as.numeric(band["lower"])
  upper <- as.numeric(band["upper"])
  fit_tbl <- tibble::as_tibble(diagnostics$fit %||% tibble::tibble())
  if (nrow(fit_tbl) == 0L ||
      !all(c("Facet", "Level", "Infit", "Outfit") %in% names(fit_tbl))) {
    return(.misfit_casebook_empty())
  }
  flagged <- fit_tbl |>
    dplyr::filter(
      (is.finite(.data$Infit) & (.data$Infit < lower | .data$Infit > upper)) |
        (is.finite(.data$Outfit) & (.data$Outfit < lower | .data$Outfit > upper))
    )
  if (nrow(flagged) == 0L) {
    return(.misfit_casebook_empty())
  }

  flagged <- flagged |>
    dplyr::mutate(
      Magnitude = pmax(
        ifelse(is.finite(.data$Outfit), abs(log(pmax(.data$Outfit, 1e-6))), 0),
        ifelse(is.finite(.data$Infit), abs(log(pmax(.data$Infit, 1e-6))), 0),
        na.rm = TRUE
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$Magnitude))

  out <- flagged |>
    dplyr::mutate(
      SourceRowKey = paste(.data$Facet, .data$Level, sep = "::"),
      CaseID = paste0("element_fit:", .data$SourceRowKey),
      Direction = dplyr::case_when(
        is.finite(.data$Outfit) & .data$Outfit > upper ~ "Outfit MnSq above band",
        is.finite(.data$Infit) & .data$Infit > upper ~ "Infit MnSq above band",
        is.finite(.data$Outfit) & .data$Outfit < lower ~ "Outfit MnSq below band",
        is.finite(.data$Infit) & .data$Infit < lower ~ "Infit MnSq below band",
        TRUE ~ "Mixed"
      )
    ) |>
    dplyr::transmute(
      CaseID = .data$CaseID,
      CaseType = "element_fit_case",
      SourceFamily = "element_fit",
      SourceTable = "diagnostics$fit",
      SourceRowKey = .data$SourceRowKey,
      PrimaryUnit = paste(.data$Facet, .data$Level, sep = ": "),
      PrimaryUnitType = "facet_level",
      Person = NA_character_,
      Facet = as.character(.data$Facet),
      Level = as.character(.data$Level),
      Category = NA_integer_,
      PairKey = NA_character_,
      ContextKey = sprintf("Infit=%.2f, Outfit=%.2f", .data$Infit, .data$Outfit),
      Wave = NA_character_,
      Signal = sprintf("MnSq misfit (band %.1f-%.1f)", lower, upper),
      Direction = .data$Direction,
      Magnitude = as.numeric(.data$Magnitude),
      ReviewPriority = as.numeric(.data$Magnitude),
      WithinSourceRank = as.integer(dplyr::row_number()),
      EvidenceN = 1L,
      SupportBasis = "legacy",
      InterpretationTier = "operational_review",
      PrimaryPlotRoute = "plot_qc_dashboard(fit, diagnostics = diagnostics, draw = FALSE)",
      SupportStatus = support_status
    )
  .misfit_casebook_standardize(out)
}

.misfit_casebook_from_displacement <- function(displacement, support_status = "supported") {
  src <- tibble::as_tibble(displacement$table %||% tibble::tibble())
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  if ("Flag" %in% names(src)) {
    src <- src |>
      dplyr::filter(isTRUE(.data$Flag) | (.data$Flag %in% TRUE))
  }
  if (nrow(src) == 0L) {
    return(.misfit_casebook_empty())
  }
  src <- src |>
    dplyr::mutate(
      AbsMagnitude = pmax(abs(.data$Displacement), abs(.data$DisplacementT), na.rm = TRUE)
    ) |>
    dplyr::arrange(dplyr::desc(.data$AbsMagnitude), dplyr::desc(abs(.data$AnchorGap)))

  out <- src |>
    dplyr::mutate(
      SourceRowKey = paste(.data$Facet, .data$Level, sep = "::"),
      CaseID = paste0("displacement:", .data$SourceRowKey),
      Direction = dplyr::if_else(
        .data$Displacement >= 0,
        "Higher than anchored",
        "Lower than anchored"
      )
    ) |>
    dplyr::transmute(
      CaseID = .data$CaseID,
      CaseType = "displacement_case",
      SourceFamily = "displacement",
      SourceTable = "displacement_table()",
      SourceRowKey = .data$SourceRowKey,
      PrimaryUnit = paste(.data$Facet, .data$Level, sep = ": "),
      PrimaryUnitType = "facet_level",
      Person = NA_character_,
      Facet = as.character(.data$Facet),
      Level = as.character(.data$Level),
      Category = NA_integer_,
      PairKey = NA_character_,
      ContextKey = as.character(.data$AnchorType),
      Wave = NA_character_,
      Signal = "Displacement screen",
      Direction = .data$Direction,
      Magnitude = as.numeric(.data$AbsMagnitude),
      ReviewPriority = as.numeric(.data$AbsMagnitude),
      WithinSourceRank = as.integer(dplyr::row_number()),
      EvidenceN = 1L,
      SupportBasis = "legacy",
      InterpretationTier = "operational_review",
      PrimaryPlotRoute = "plot_displacement(displacement, draw = FALSE)",
      SupportStatus = support_status
    )
  .misfit_casebook_standardize(out)
}

#' Build a case-level misfit review bundle
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param unexpected Optional output from [unexpected_response_table()].
#' @param displacement Optional output from [displacement_table()].
#' @param administration_id Optional scalar identifier describing the current
#'   administration or form. It is stored in row-level provenance and summary
#'   outputs when supplied.
#' @param wave_id Optional scalar identifier for the current wave or occasion.
#'   It is stored in row-level provenance and summary outputs when supplied.
#' @param top_n Maximum number of rows to keep in compact summary outputs.
#'
#' @details
#' `build_misfit_casebook()` is a synthesis layer over package-native screening
#' outputs. It does not invent a new misfit statistic. Instead, it organizes
#' existing evidence families into one case-level review surface:
#'
#' - element-level Infit / Outfit MnSq misfit from `diagnostics$fit`
#'   (rows whose Infit or Outfit MnSq falls outside the 0.5-1.5 Linacre
#'   acceptance band)
#' - strict marginal cell screens from `diagnostics$marginal_fit$top_cells`
#' - strict pairwise screens from `diagnostics$marginal_fit$pairwise$top_pairs`
#' - unexpected responses from [unexpected_response_table()]
#' - displacement flags from [displacement_table()]
#'
#' The result is an operational review bundle. It is not a formal adjudication
#' system, and repeated signals across evidence families should be prioritized
#' over any single isolated case row. In addition to raw case rows, the object
#' includes stable grouping views such as `by_person`, `by_facet_level`,
#' `by_source_family`, and `by_wave` to support operational triage. The
#' `source_support` component records which evidence families are currently
#' supported, caveated, or deferred under the active model.
#'
#' @section Recommended input route:
#' 1. Fit with [fit_mfrm()].
#' 2. Build diagnostics with [diagnose_mfrm()].
#' 3. Optionally build [unexpected_response_table()] and [displacement_table()]
#'    yourself when you want custom thresholds before synthesizing the casebook.
#'
#' @section GPCM boundary:
#' For bounded `GPCM`, the helper is available with caveat. The casebook inherits
#' exploratory screening semantics from the underlying residual and strict
#' marginal sources; it should not be read as a formal inferential case test.
#'
#' @return An object of class `mfrm_misfit_casebook`.
#' @seealso [diagnose_mfrm()], [unexpected_response_table()],
#'   [displacement_table()], [plot_unexpected()], [plot_displacement()],
#'   [plot_marginal_fit()], [plot_marginal_pairwise()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", model = "RSM", quad_points = 5)
#' diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
#' casebook <- build_misfit_casebook(fit, diagnostics = diag, top_n = 10)
#' summary(casebook)
#' casebook$top_cases
#' }
#' @export
build_misfit_casebook <- function(fit,
                                  diagnostics = NULL,
                                  unexpected = NULL,
                                  displacement = NULL,
                                  administration_id = NULL,
                                  wave_id = NULL,
                                  top_n = 25) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  diagnostics <- .validate_misfit_casebook_input(diagnostics, "diagnostics", "mfrm_diagnostics")
  unexpected <- .validate_misfit_casebook_input(unexpected, "unexpected", "mfrm_unexpected")
  displacement <- .validate_misfit_casebook_input(displacement, "displacement", "mfrm_displacement")
  administration_id <- .normalize_misfit_casebook_id(administration_id, "administration_id")
  wave_id <- .normalize_misfit_casebook_id(wave_id, "wave_id")
  top_n <- max(1L, as.integer(top_n %||% 25L))

  if (is.null(diagnostics)) {
    diag_mode <- if (identical(as.character(fit$config$model %||% NA_character_), "GPCM")) {
      "both"
    } else {
      "both"
    }
    diagnostics <- diagnose_mfrm(fit, diagnostic_mode = diag_mode, residual_pca = "none")
  }
  if (is.null(unexpected)) {
    unexpected <- unexpected_response_table(fit, diagnostics = diagnostics, top_n = max(top_n, 50L))
  }
  if (is.null(displacement)) {
    displacement <- displacement_table(fit, diagnostics = diagnostics, anchored_only = FALSE, top_n = max(top_n, 50L))
  }

  model <- as.character(fit$config$model %||% NA_character_)[1]
  gpcm_status <- if (identical(model, "GPCM")) "supported_with_caveat" else "supported"
  support_status <- .misfit_casebook_support_status(model)
  source_support <- .misfit_casebook_source_support(
    model = model,
    diagnostics = diagnostics,
    unexpected = unexpected,
    displacement = displacement
  )

  marginal_cell_cases <- .misfit_casebook_from_marginal_cells(diagnostics, support_status = gpcm_status)
  pairwise_cases <- .misfit_casebook_from_pairwise(diagnostics, support_status = gpcm_status)
  unexpected_cases <- .misfit_casebook_from_unexpected(
    unexpected = unexpected,
    facet_names = as.character(fit$config$facet_names %||% character(0)),
    support_status = gpcm_status
  )
  displacement_cases <- .misfit_casebook_from_displacement(displacement, support_status = gpcm_status)
  element_fit_cases <- .misfit_casebook_from_element_fit(diagnostics, support_status = gpcm_status)

  marginal_cell_cases <- .misfit_casebook_apply_provenance(marginal_cell_cases, administration_id = administration_id, wave_id = wave_id)
  pairwise_cases <- .misfit_casebook_apply_provenance(pairwise_cases, administration_id = administration_id, wave_id = wave_id)
  unexpected_cases <- .misfit_casebook_apply_provenance(unexpected_cases, administration_id = administration_id, wave_id = wave_id)
  displacement_cases <- .misfit_casebook_apply_provenance(displacement_cases, administration_id = administration_id, wave_id = wave_id)
  element_fit_cases <- .misfit_casebook_apply_provenance(element_fit_cases, administration_id = administration_id, wave_id = wave_id)

  all_cases <- dplyr::bind_rows(
    marginal_cell_cases,
    pairwise_cases,
    unexpected_cases,
    displacement_cases,
    element_fit_cases
  )
  if (nrow(all_cases) > 0) {
    all_cases <- all_cases |>
      dplyr::arrange(dplyr::desc(.data$ReviewPriority), .data$SourceFamily, .data$WithinSourceRank)
  } else {
    all_cases <- .misfit_casebook_empty()
  }

  review_status <- if (nrow(all_cases) > 0L) "review_required" else "no_flagged_cases"
  case_rollup <- .misfit_casebook_build_rollup(all_cases)
  group_views <- .misfit_casebook_build_group_views(all_cases, case_rollup)
  group_view_index <- .misfit_casebook_group_view_index(group_views)
  source_summary <- if (nrow(all_cases) == 0L) {
    tibble::tibble()
  } else {
    all_cases |>
      dplyr::group_by(.data$SourceFamily, .data$SupportBasis, .data$InterpretationTier, .data$PrimaryPlotRoute) |>
      dplyr::summarize(
        Cases = dplyr::n(),
        MaxPriority = max(.data$ReviewPriority, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data$Cases), dplyr::desc(.data$MaxPriority))
  }

  top_cases <- all_cases
  if (nrow(top_cases) > 0L) {
    top_cases <- top_cases |>
      dplyr::slice_head(n = top_n)
  }

  key_warnings <- character(0)
  if (nrow(marginal_cell_cases) > 0L) {
    key_warnings <- c(key_warnings, paste0("Strict marginal cell screening contributed ", nrow(marginal_cell_cases), " flagged case rows."))
  }
  if (nrow(pairwise_cases) > 0L) {
    key_warnings <- c(key_warnings, paste0("Strict pairwise screening contributed ", nrow(pairwise_cases), " flagged pair rows."))
  }
  if (nrow(unexpected_cases) > 0L) {
    key_warnings <- c(key_warnings, paste0("Unexpected-response screening contributed ", nrow(unexpected_cases), " case rows."))
  }
  if (nrow(displacement_cases) > 0L) {
    key_warnings <- c(key_warnings, paste0("Displacement screening contributed ", nrow(displacement_cases), " flagged facet-level rows."))
  }
  key_warnings <- clean_summary_lines(key_warnings, max_n = 4L)
  if (length(key_warnings) == 0L) {
    key_warnings <- "No flagged casebook rows met the current source-specific rules."
  }

  next_actions <- character(0)
  if (nrow(marginal_cell_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot_marginal_fit(diagnostics, draw = FALSE) to inspect the largest first-order strict marginal cells.")
  }
  if (nrow(pairwise_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot_marginal_pairwise(diagnostics, draw = FALSE) to inspect the strongest pairwise local-dependence signals.")
  }
  if (nrow(unexpected_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot_unexpected(unexpected, draw = FALSE) to review the most surprising person-level observations.")
  }
  if (nrow(displacement_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot_displacement(displacement, draw = FALSE) when flagged facet levels suggest anchor or stability review.")
  }
  if (length(next_actions) == 0L) {
    next_actions <- "If no case rows are flagged, continue with reporting_checklist() or linking review as needed."
  }
  next_actions <- clean_summary_lines(next_actions, max_n = 4L)

  status <- make_summary_block(
    "Overall status" = review_status,
    "Model" = model,
    "Administration ID" = administration_id,
    "Wave ID" = wave_id,
    "Bounded GPCM" = as.character(support_status$Status[support_status$Scope == "bounded GPCM"][1] %||% NA_character_)
  )

  plot_map <- tibble::tibble(
    SourceFamily = c("marginal_cell", "marginal_pair", "unexpected", "displacement"),
    Available = c(
      nrow(marginal_cell_cases) > 0L,
      nrow(pairwise_cases) > 0L,
      nrow(unexpected_cases) > 0L,
      nrow(displacement_cases) > 0L
    ),
    PlotHelper = c(
      "plot_marginal_fit(diagnostics, draw = FALSE)",
      "plot_marginal_pairwise(diagnostics, draw = FALSE)",
      "plot_unexpected(unexpected, draw = FALSE)",
      "plot_displacement(displacement, draw = FALSE)"
    ),
    Trigger = c(
      "Use when strict first-order category cells are flagged.",
      "Use when strict pairwise local-dependence rows are flagged.",
      "Use when person-level unexpected responses dominate review.",
      "Use when anchor or facet-level displacement rows are flagged."
    )
  )

  reporting_map <- tibble::tibble(
    Area = c(
      "Operational misfit review",
      "Strict marginal screening follow-up",
      "Observation-level appendix",
      "Manuscript/reporting companion"
    ),
    CoveredHere = c("yes", "partial", "partial", "partial"),
    CompanionOutput = c(
      "summary(build_misfit_casebook(...))",
      "plot_marginal_fit() / plot_marginal_pairwise()",
      "unexpected_response_table() / displacement_table()",
      "reporting_checklist() / build_summary_table_bundle() / build_apa_outputs()"
    )
  )

  overview <- tibble::tibble(
    Model = model,
    DiagnosticMode = as.character(diagnostics$diagnostic_mode %||% NA_character_),
    ReviewStatus = review_status,
    AdministrationID = administration_id,
    WaveID = wave_id,
    TotalCases = nrow(all_cases),
    RollupRows = nrow(case_rollup),
    GroupViews = sum(group_view_index$Rows > 0L),
    TopCaseRows = nrow(top_cases),
    SourcesAvailable = sum(c(
      nrow(marginal_cell_cases) > 0L,
      nrow(pairwise_cases) > 0L,
      nrow(unexpected_cases) > 0L,
      nrow(displacement_cases) > 0L
    )),
    GPCMSupport = as.character(support_status$Status[support_status$Scope == "bounded GPCM"][1] %||% NA_character_)
  )

  notes <- clean_summary_lines(c(
    "Misfit casebook rows are operational review units, not formal case decisions.",
    "The helper preserves source-family-specific screening logic rather than collapsing all evidence into one opaque score.",
    "Repeated signals across strict marginal, unexpected-response, and displacement sources deserve priority."
  ))

  out <- list(
    overview = overview,
    status = status,
    key_warnings = key_warnings,
    next_actions = next_actions,
    top_cases = top_cases,
    case_rollup = case_rollup,
    group_view_index = group_view_index,
    group_views = group_views,
    source_summary = source_summary,
    source_support = source_support,
    marginal_cell_cases = marginal_cell_cases,
    pairwise_cases = pairwise_cases,
    unexpected_cases = unexpected_cases,
    displacement_cases = displacement_cases,
    plot_map = plot_map,
    reporting_map = reporting_map,
    support_status = support_status,
    notes = notes,
    settings = list(
      top_n = top_n,
      model = model,
      administration_id = administration_id,
      wave_id = wave_id,
      intended_use = "operational_misfit_review",
      source_profile = data.frame(
        Source = c("marginal_cell", "marginal_pair", "unexpected", "displacement"),
        Available = c(
          nrow(marginal_cell_cases) > 0L,
          nrow(pairwise_cases) > 0L,
          nrow(unexpected_cases) > 0L,
          nrow(displacement_cases) > 0L
        ),
        stringsAsFactors = FALSE
      )
    )
  )
  as_mfrm_bundle(out, "mfrm_misfit_casebook")
}

#' @export
print.mfrm_misfit_casebook <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Summarize a misfit-casebook object
#'
#' @param object Output from [build_misfit_casebook()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of top case rows to keep in the compact summary.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_misfit_casebook`.
#' @seealso [build_misfit_casebook()]
#' @export
summary.mfrm_misfit_casebook <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_misfit_casebook")) {
    stop("`object` must be output from build_misfit_casebook().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  top_cases <- tibble::as_tibble(object$top_cases %||% tibble::tibble())
  if (nrow(top_cases) > 0L) {
    top_cases <- top_cases |>
      dplyr::slice_head(n = top_n)
  }

  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    status = tibble::as_tibble(object$status %||% tibble::tibble()),
    key_warnings = clean_summary_lines(object$key_warnings %||% character(0), max_n = 4L),
    next_actions = clean_summary_lines(object$next_actions %||% character(0), max_n = 4L),
    top_cases = top_cases,
    case_rollup = tibble::as_tibble(object$case_rollup %||% tibble::tibble()),
    group_view_index = tibble::as_tibble(object$group_view_index %||% tibble::tibble()),
    group_views = object$group_views %||% list(),
    source_summary = tibble::as_tibble(object$source_summary %||% tibble::tibble()),
    source_support = tibble::as_tibble(object$source_support %||% tibble::tibble()),
    plot_routes = .review_plot_routes(object$plot_map),
    plot_map = tibble::as_tibble(object$plot_map %||% tibble::tibble()),
    reporting_map = tibble::as_tibble(object$reporting_map %||% tibble::tibble()),
    support_status = tibble::as_tibble(object$support_status %||% tibble::tibble()),
    notes = clean_summary_lines(object$notes %||% character(0)),
    settings = object$settings %||% list(),
    digits = digits
  )
  class(out) <- "summary.mfrm_misfit_casebook"
  out
}

#' @export
print.summary.mfrm_misfit_casebook <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Misfit Casebook Summary\n")
  if (nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$status) > 0) {
    cat("\nStatus\n")
    print(as.data.frame(x$status), row.names = FALSE)
  }
  print_bullet_section("Key Warnings", x$key_warnings)
  print_bullet_section("Next Actions", x$next_actions)
  if (nrow(x$case_rollup) > 0) {
    cat("\nCase Rollup\n")
    print(round_numeric_df(as.data.frame(x$case_rollup), digits = digits), row.names = FALSE)
  }
  if (nrow(x$group_view_index) > 0) {
    cat("\nGrouping Views\n")
    print(as.data.frame(x$group_view_index), row.names = FALSE)
  }
  if (nrow(x$plot_routes) > 0) {
    cat("\nPlot Follow-up\n")
    print(as.data.frame(x$plot_routes), row.names = FALSE)
  }
  if (nrow(x$source_summary) > 0) {
    cat("\nSource Summary\n")
    print(round_numeric_df(as.data.frame(x$source_summary), digits = digits), row.names = FALSE)
  }
  if (nrow(x$source_support) > 0) {
    cat("\nSource Support\n")
    print(as.data.frame(x$source_support), row.names = FALSE)
  }
  if (nrow(x$top_cases) > 0) {
    cat("\nTop Cases\n")
    print(round_numeric_df(as.data.frame(x$top_cases), digits = digits), row.names = FALSE)
  }
  if (nrow(x$support_status) > 0) {
    cat("\nSupport Status\n")
    print(as.data.frame(x$support_status), row.names = FALSE)
  }
  print_bullet_section("Notes", x$notes)
  invisible(x)
}

# --- build_weighting_review --------------------------------------------------

.validate_weighting_review_fit <- function(x, arg, models) {
  if (!inherits(x, "mfrm_fit")) {
    stop("`", arg, "` must be an `mfrm_fit` object from fit_mfrm().", call. = FALSE)
  }
  model <- as.character(x$config$model %||% x$summary$Model[1] %||% NA_character_)[1]
  if (!model %in% models) {
    stop(
      "`", arg, "` must use one of: ", paste(models, collapse = ", "),
      ". Got: ", model, ".",
      call. = FALSE
    )
  }
  x
}

.weighting_review_support_status <- function() {
  tibble::tibble(
    Scope = c("RSM / PCM reference", "bounded GPCM comparison"),
    Status = c("supported", "supported_with_caveat"),
    Note = c(
      "Supported as the equal-weighting reference side of the review.",
      paste(
        "Supported with caveat as a slope-aware comparison model.",
        "Use it to inspect discrimination-based reweighting, not as an automatic replacement",
        "for the Rasch-family route."
      )
    )
  )
}

.weighting_review_facet_shift <- function(rasch_fit, gpcm_fit) {
  ref_tbl <- tibble::as_tibble(rasch_fit$facets$others %||% tibble::tibble())
  gpcm_tbl <- tibble::as_tibble(gpcm_fit$facets$others %||% tibble::tibble())
  if (nrow(ref_tbl) == 0L || nrow(gpcm_tbl) == 0L) {
    return(tibble::tibble())
  }

  ref_tbl <- ref_tbl |>
    dplyr::rename(
      ReferenceEstimate = "Estimate"
    ) |>
    dplyr::group_by(.data$Facet) |>
    dplyr::arrange(.data$ReferenceEstimate, .by_group = TRUE) |>
    dplyr::mutate(
      ReferenceRank = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  gpcm_tbl <- gpcm_tbl |>
    dplyr::rename(
      ComparisonEstimate = "Estimate"
    ) |>
    dplyr::group_by(.data$Facet) |>
    dplyr::arrange(.data$ComparisonEstimate, .by_group = TRUE) |>
    dplyr::mutate(
      ComparisonRank = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  ref_tbl |>
    dplyr::inner_join(gpcm_tbl, by = c("Facet", "Level")) |>
    dplyr::mutate(
      DeltaEstimate = .data$ComparisonEstimate - .data$ReferenceEstimate,
      AbsDeltaEstimate = abs(.data$DeltaEstimate),
      RankShift = .data$ComparisonRank - .data$ReferenceRank,
      Direction = dplyr::case_when(
        .data$DeltaEstimate > 0 ~ "Higher in bounded GPCM",
        .data$DeltaEstimate < 0 ~ "Lower in bounded GPCM",
        TRUE ~ "No change"
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$AbsDeltaEstimate), dplyr::desc(abs(.data$RankShift)), .data$Facet, .data$Level)
}

.weighting_review_slope_profile <- function(gpcm_fit) {
  slope_tbl <- tibble::as_tibble(gpcm_fit$slopes %||% tibble::tibble())
  if (nrow(slope_tbl) == 0L) {
    return(tibble::tibble())
  }

  slope_facet <- as.character(gpcm_fit$config$slope_facet %||% NA_character_)[1]
  obs_df <- as.data.frame(gpcm_fit$prep$data %||% NULL, stringsAsFactors = FALSE)
  exposure_tbl <- tibble::tibble(SlopeFacet = character(), Exposure = numeric(), ExposureShare = numeric())
  if (nrow(obs_df) > 0L && slope_facet %in% names(obs_df)) {
    w <- suppressWarnings(as.numeric(obs_df$Weight %||% rep(1, nrow(obs_df))))
    w[!is.finite(w)] <- 0
    obs_df$..Exposure <- w
    exposure_tbl <- obs_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(slope_facet))) |>
      dplyr::summarize(Exposure = sum(.data$..Exposure, na.rm = TRUE), .groups = "drop")
    names(exposure_tbl)[1] <- "SlopeFacet"
    exposure_tbl <- tibble::as_tibble(exposure_tbl) |>
      dplyr::mutate(
        ExposureShare = ifelse(sum(.data$Exposure, na.rm = TRUE) > 0,
                               .data$Exposure / sum(.data$Exposure, na.rm = TRUE),
                               NA_real_)
      )
  }

  slope_tbl |>
    dplyr::left_join(exposure_tbl, by = "SlopeFacet") |>
    dplyr::mutate(
      RelativeWeight = .data$Estimate,
      AbsLogDeviation = abs(.data$LogEstimate),
      WeightingDirection = dplyr::case_when(
        .data$Estimate > 1.05 ~ "Upweighted",
        .data$Estimate < 0.95 ~ "Downweighted",
        TRUE ~ "Near unit"
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$AbsLogDeviation), dplyr::desc(dplyr::coalesce(.data$Exposure, 0)))
}

.weighting_review_information_profile <- function(fit,
                                                 theta_range,
                                                 theta_points) {
  info <- compute_information(fit, theta_range = theta_range, theta_points = theta_points)
  tibble::as_tibble(info$iif) |>
    dplyr::group_by(.data$Facet, .data$Level) |>
    dplyr::summarize(
      IntegratedInfo = sum(.data$Information, na.rm = TRUE),
      Exposure = max(.data$Exposure, na.rm = TRUE),
      .groups = "drop_last"
    ) |>
    dplyr::mutate(
      InfoShare = ifelse(sum(.data$IntegratedInfo, na.rm = TRUE) > 0,
                         .data$IntegratedInfo / sum(.data$IntegratedInfo, na.rm = TRUE),
                         NA_real_),
      ExposureShare = ifelse(sum(.data$Exposure, na.rm = TRUE) > 0,
                             .data$Exposure / sum(.data$Exposure, na.rm = TRUE),
                             NA_real_)
    ) |>
    dplyr::ungroup()
}

#' Build a weighting-policy review between Rasch-family and bounded GPCM fits
#'
#' @param rasch_fit Output from [fit_mfrm()] using `model = "RSM"` or `"PCM"`.
#' @param gpcm_fit Output from [fit_mfrm()] using bounded `model = "GPCM"`.
#' @param theta_range Numeric vector of length 2 passed to [compute_information()]
#'   for the information-redistribution comparison.
#' @param theta_points Integer number of theta grid points passed to
#'   [compute_information()].
#' @param top_n Maximum number of rows to keep in compact summary outputs.
#'
#' @details
#' `build_weighting_review()` is an operational model-choice review helper. It
#' is designed for the common question:
#'
#' - what changes when a Rasch-family equal-weighting model is replaced with a
#'   bounded `GPCM` that allows discrimination-based reweighting?
#'
#' The helper does not estimate a new model. Instead, it synthesizes four
#' package-native evidence sources:
#'
#' - [compare_mfrm()] for same-data model comparison
#' - the non-person facet measures from each fit
#' - the bounded `GPCM` slope table
#' - [compute_information()] for design-weighted information redistribution
#'
#' The result is intended for substantive review, not for automatic model
#' selection. In particular, a better-fitting `GPCM` should not by itself be
#' interpreted as a reason to discard an equal-weighting Rasch-family route.
#'
#' @section Recommended input route:
#' 1. Fit an equal-weighting reference model with `model = "RSM"` or `"PCM"`.
#' 2. Fit a bounded `GPCM` on the same prepared response data.
#' 3. Run `build_weighting_review(rasch_fit, gpcm_fit)`.
#' 4. Read `summary(review)` before deciding whether the discrimination-based
#'    reweighting is substantively acceptable.
#'
#' @section What the returned tables mean:
#' - `model_comparison`: same-data model-comparison bundle from [compare_mfrm()].
#' - `facet_shift`: how non-person facet estimates move under bounded `GPCM`.
#' - `slope_profile`: which `slope_facet` levels are upweighted or downweighted.
#' - `information_redistribution`: within-facet information-share changes
#'   between the Rasch-family fit and bounded `GPCM`.
#' - `top_reweighted_levels`: compact triage table for the strongest
#'   slope-facet-level redistribution signals.
#'
#' @section GPCM boundary:
#' This helper is available only for the current bounded `GPCM` branch. It
#' requires the package's existing `slope_facet == step_facet` contract and
#' should be read as an operational weighting-policy review, not as a formal
#' validity adjudication.
#'
#' @return An object of class `mfrm_weighting_review`.
#' @seealso [compare_mfrm()], [compute_information()], [gpcm_capability_matrix()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' rasch_fit <- fit_mfrm(
#'   toy,
#'   "Person",
#'   c("Rater", "Criterion"),
#'   "Score",
#'   method = "MML",
#'   model = "RSM",
#'   quad_points = 9
#' )
#' gpcm_fit <- fit_mfrm(
#'   toy,
#'   "Person",
#'   c("Rater", "Criterion"),
#'   "Score",
#'   method = "MML",
#'   model = "GPCM",
#'   step_facet = "Criterion",
#'   slope_facet = "Criterion",
#'   quad_points = 9
#' )
#' review <- build_weighting_review(rasch_fit, gpcm_fit, theta_points = 41)
#' summary(review)
#' review$top_reweighted_levels
#' }
#' @name build_weighting_review
#' @export
build_weighting_review <- function(rasch_fit,
                                   gpcm_fit,
                                   theta_range = c(-6, 6),
                                   theta_points = 101L,
                                   top_n = 10L) {
  rasch_fit <- .validate_weighting_review_fit(rasch_fit, "rasch_fit", c("RSM", "PCM"))
  gpcm_fit <- .validate_weighting_review_fit(gpcm_fit, "gpcm_fit", "GPCM")
  top_n <- max(1L, as.integer(top_n %||% 10L))
  theta_points <- max(11L, as.integer(theta_points %||% 101L))

  slope_facet <- as.character(gpcm_fit$config$slope_facet %||% NA_character_)[1]
  step_facet <- as.character(gpcm_fit$config$step_facet %||% NA_character_)[1]
  if (!identical(slope_facet, step_facet)) {
    stop(
      "build_weighting_review() currently requires the bounded `GPCM` branch with `slope_facet == step_facet`.",
      call. = FALSE
    )
  }

  comparison <- suppressWarnings(compare_mfrm(
    rasch_fit,
    gpcm_fit,
    labels = c(
      paste0(as.character(rasch_fit$config$model %||% "RSM")[1], " reference"),
      "bounded GPCM"
    ),
    warn_constraints = FALSE,
    nested = FALSE
  ))

  basis <- comparison$comparison_basis %||% list()
  if (!isTRUE(basis$same_data)) {
    stop(
      "build_weighting_review() requires the two fits to share the same prepared response data.",
      call. = FALSE
    )
  }
  if (!identical(as.character(rasch_fit$config$facet_names %||% character(0)),
                 as.character(gpcm_fit$config$facet_names %||% character(0)))) {
    stop(
      "The two fits must use the same facet columns in the same order.",
      call. = FALSE
    )
  }

  facet_shift <- .weighting_review_facet_shift(rasch_fit, gpcm_fit)
  slope_profile <- .weighting_review_slope_profile(gpcm_fit)
  rasch_info <- .weighting_review_information_profile(rasch_fit, theta_range = theta_range, theta_points = theta_points)
  gpcm_info <- .weighting_review_information_profile(gpcm_fit, theta_range = theta_range, theta_points = theta_points)

  information_redistribution <- rasch_info |>
    dplyr::rename(
      ReferenceIntegratedInfo = "IntegratedInfo",
      ReferenceExposure = "Exposure",
      ReferenceInfoShare = "InfoShare",
      ReferenceExposureShare = "ExposureShare"
    ) |>
    dplyr::inner_join(
      gpcm_info |>
        dplyr::rename(
          ComparisonIntegratedInfo = "IntegratedInfo",
          ComparisonExposure = "Exposure",
          ComparisonInfoShare = "InfoShare",
          ComparisonExposureShare = "ExposureShare"
        ),
      by = c("Facet", "Level")
    ) |>
    dplyr::mutate(
      InfoShareDelta = .data$ComparisonInfoShare - .data$ReferenceInfoShare,
      ExposureShareDelta = .data$ComparisonExposureShare - .data$ReferenceExposureShare,
      IntegratedInfoRatio = dplyr::if_else(
        is.finite(.data$ReferenceIntegratedInfo) & .data$ReferenceIntegratedInfo > 0,
        .data$ComparisonIntegratedInfo / .data$ReferenceIntegratedInfo,
        NA_real_
      ),
      AbsInfoShareDelta = abs(.data$InfoShareDelta)
    ) |>
    dplyr::mutate(
      AbsLogInfoRatio = abs(log(dplyr::coalesce(.data$IntegratedInfoRatio, 1)))
    ) |>
    dplyr::arrange(dplyr::desc(.data$AbsInfoShareDelta), dplyr::desc(.data$AbsLogInfoRatio), .data$Facet, .data$Level)

  top_reweighted_levels <- information_redistribution |>
    dplyr::filter(.data$Facet == slope_facet) |>
    dplyr::left_join(
      slope_profile |>
        dplyr::select("SlopeFacet", "Estimate", "LogEstimate", "WeightingDirection", "Exposure", "ExposureShare"),
      by = c("Level" = "SlopeFacet"),
      suffix = c("", "_Slope")
    ) |>
    dplyr::rename(
      SlopeEstimate = "Estimate",
      SlopeLogEstimate = "LogEstimate",
      SlopeDirection = "WeightingDirection",
      SlopeExposure = "Exposure",
      SlopeExposureShare = "ExposureShare"
    ) |>
    dplyr::slice_head(n = top_n)

  top_measure_shifts <- facet_shift |>
    dplyr::slice_head(n = top_n)

  max_abs_log_slope <- if (nrow(slope_profile) > 0L) max(slope_profile$AbsLogDeviation, na.rm = TRUE) else 0
  max_abs_info_delta <- if (nrow(information_redistribution) > 0L) max(information_redistribution$AbsInfoShareDelta, na.rm = TRUE) else 0
  review_status <- dplyr::case_when(
    max_abs_log_slope <= log(1.05) && max_abs_info_delta <= 0.02 ~ "minimal_reweighting_detected",
    TRUE ~ "reweighting_review_required"
  )

  support_status <- .weighting_review_support_status()
  comparison_mode <- if (isTRUE(basis$ic_comparable)) "same_basis_fit_comparison" else "descriptive_model_contrast_only"
  overview <- tibble::tibble(
    ReferenceModel = as.character(rasch_fit$config$model %||% NA_character_)[1],
    ComparisonModel = as.character(gpcm_fit$config$model %||% NA_character_)[1],
    ReferenceMethod = public_mfrm_method_label(as.character(rasch_fit$config$method %||% NA_character_)[1]),
    ComparisonMethod = public_mfrm_method_label(as.character(gpcm_fit$config$method %||% NA_character_)[1]),
    SlopeFacet = slope_facet,
    ReviewStatus = review_status,
    ComparisonMode = comparison_mode,
    MaxAbsLogSlope = max_abs_log_slope,
    MaxAbsInfoShareDelta = max_abs_info_delta
  )

  status <- make_summary_block(
    "Overall status" = review_status,
    "Weighting principle" = "Rasch-family equal weighting vs bounded GPCM discrimination-based reweighting",
    "Comparison basis" = comparison_mode
  )

  key_warnings <- character(0)
  if (!isTRUE(basis$ic_comparable)) {
    key_warnings <- c(
      key_warnings,
      "Model-comparison weights are descriptive only because the two fits do not share a fully comparable formal MML basis."
    )
  }
  if (nrow(slope_profile) > 0L) {
    lead_slope <- slope_profile[1, , drop = FALSE]
    key_warnings <- c(
      key_warnings,
      paste0(
        "Largest bounded GPCM slope deviation is at ",
        slope_facet, " = ", lead_slope$SlopeFacet[[1]],
        " (Estimate = ", format(round(lead_slope$Estimate[[1]], 3), nsmall = 3), ")."
      )
    )
  }
  if (nrow(top_reweighted_levels) > 0L) {
    lead_delta <- top_reweighted_levels[1, , drop = FALSE]
    key_warnings <- c(
      key_warnings,
      paste0(
        "Largest within-facet information-share shift is ",
        format(round(lead_delta$InfoShareDelta[[1]], 3), nsmall = 3),
        " for ", slope_facet, " = ", lead_delta$Level[[1]], "."
      )
    )
  }
  if (nrow(top_measure_shifts) > 0L) {
    lead_shift <- top_measure_shifts[1, , drop = FALSE]
    key_warnings <- c(
      key_warnings,
      paste0(
        "Largest facet-measure shift is ",
        format(round(lead_shift$DeltaEstimate[[1]], 3), nsmall = 3),
        " for ", lead_shift$Facet[[1]], " = ", lead_shift$Level[[1]], "."
      )
    )
  }
  key_warnings <- clean_summary_lines(key_warnings, max_n = 4L)

  next_actions <- clean_summary_lines(c(
    "Read summary(model_comparison) before interpreting any fit advantage as a scoring recommendation.",
    paste0("Use slope_profile and top_reweighted_levels to inspect whether ", slope_facet, " levels are being upweighted or downweighted in substantively acceptable ways."),
    paste0("Use plot_information(compute_information(rasch_fit), type = \"iif\", facet = \"", slope_facet, "\", draw = FALSE) and the bounded GPCM analogue to inspect precision redistribution visually."),
    "If equal contributions of items and raters are part of the score interpretation, retain the Rasch-family fit as the operational reference even when bounded GPCM fits better."
  ), max_n = 4L)

  plot_map <- tibble::tibble(
    ReviewArea = c("Reference precision", "bounded GPCM precision", "Fit comparison"),
    Available = c(TRUE, TRUE, TRUE),
    PlotHelper = c(
      paste0("plot_information(compute_information(rasch_fit), type = \"iif\", facet = \"", slope_facet, "\", draw = FALSE)"),
      paste0("plot_information(compute_information(gpcm_fit), type = \"iif\", facet = \"", slope_facet, "\", draw = FALSE)"),
      "summary(compare_mfrm(rasch_fit, gpcm_fit))"
    ),
    Trigger = c(
      "Use to inspect the equal-weighting reference precision split across slope-facet levels.",
      "Use to inspect how bounded GPCM redistributes precision across the same levels.",
      "Use to review AIC/BIC and evidence ratios before making a model-choice argument."
    )
  )

  reporting_map <- tibble::tibble(
    Area = c(
      "Operational weighting review",
      "Model-comparison companion",
      "Score-semantics follow-up"
    ),
    CoveredHere = c("yes", "partial", "partial"),
    CompanionOutput = c(
      "summary(build_weighting_review(...))",
      "compare_mfrm() / summary(compare_mfrm(...))",
      "fair_average_table() for the Rasch-family route; keep bounded GPCM score semantics separate"
    )
  )

  notes <- clean_summary_lines(c(
    "Observation weights and discrimination-based reweighting are separate concepts in this package.",
    "The review is intended to make reweighting visible; it does not decide by itself whether bounded GPCM should replace the Rasch-family operational model.",
    "Information-share changes are computed within each facet because the same total information is partitioned separately by facet."
  ))

  out <- list(
    overview = overview,
    status = status,
    key_warnings = key_warnings,
    next_actions = next_actions,
    model_comparison = comparison,
    facet_shift = facet_shift,
    top_measure_shifts = top_measure_shifts,
    slope_profile = slope_profile,
    information_redistribution = information_redistribution,
    top_reweighted_levels = top_reweighted_levels,
    plot_map = plot_map,
    reporting_map = reporting_map,
    support_status = support_status,
    notes = notes,
    settings = list(
      theta_range = theta_range,
      theta_points = theta_points,
      top_n = top_n,
      slope_facet = slope_facet,
      intended_use = "weighting_policy_review"
    )
  )
  as_mfrm_bundle(out, "mfrm_weighting_review")
}

#' @export
print.mfrm_weighting_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Summarize a weighting-review object
#'
#' @param object Output from [build_weighting_review()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of top rows to retain in compact summary tables.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_weighting_review`.
#' @seealso [build_weighting_review()]
#' @export
summary.mfrm_weighting_review <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_weighting_review")) {
    stop("`object` must be output from build_weighting_review().", call. = FALSE)
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    status = tibble::as_tibble(object$status %||% tibble::tibble()),
    key_warnings = clean_summary_lines(object$key_warnings %||% character(0), max_n = 4L),
    next_actions = clean_summary_lines(object$next_actions %||% character(0), max_n = 4L),
    top_measure_shifts = tibble::as_tibble(object$top_measure_shifts %||% tibble::tibble()) |>
      dplyr::slice_head(n = top_n),
    top_reweighted_levels = tibble::as_tibble(object$top_reweighted_levels %||% tibble::tibble()) |>
      dplyr::slice_head(n = top_n),
    plot_map = tibble::as_tibble(object$plot_map %||% tibble::tibble()),
    reporting_map = tibble::as_tibble(object$reporting_map %||% tibble::tibble()),
    support_status = tibble::as_tibble(object$support_status %||% tibble::tibble()),
    notes = clean_summary_lines(object$notes %||% character(0)),
    settings = object$settings %||% list(),
    digits = digits
  )
  class(out) <- "summary.mfrm_weighting_review"
  out
}

#' @export
print.summary.mfrm_weighting_review <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Weighting Review Summary\n")
  if (nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$status) > 0) {
    cat("\nStatus\n")
    print(as.data.frame(x$status), row.names = FALSE)
  }
  print_bullet_section("Key Warnings", x$key_warnings)
  print_bullet_section("Next Actions", x$next_actions)
  if (nrow(x$top_measure_shifts) > 0) {
    cat("\nTop Measure Shifts\n")
    print(round_numeric_df(as.data.frame(x$top_measure_shifts), digits = digits), row.names = FALSE)
  }
  if (nrow(x$top_reweighted_levels) > 0) {
    cat("\nTop Reweighted Levels\n")
    print(round_numeric_df(as.data.frame(x$top_reweighted_levels), digits = digits), row.names = FALSE)
  }
  if (nrow(x$support_status) > 0) {
    cat("\nSupport Status\n")
    print(as.data.frame(x$support_status), row.names = FALSE)
  }
  print_bullet_section("Notes", x$notes)
  invisible(x)
}

# --- build_model_choice_review ----------------------------------------------

.model_choice_fit_labels <- function(fits, labels = NULL) {
  if (!is.null(labels)) {
    labels <- as.character(labels)
    if (length(labels) != length(fits)) {
      stop("`labels` must have the same length as the number of fits.",
           call. = FALSE)
    }
    return(labels)
  }
  labels <- names(fits)
  if (is.null(labels) || any(!nzchar(labels))) {
    labels <- vapply(fits, function(fit) {
      model <- as.character(fit$config$model %||% "?")[1]
      method <- public_mfrm_method_label(as.character(fit$config$method %||% "?")[1])
      paste0(model, "/", method)
    }, character(1))
  }
  if (anyDuplicated(labels)) labels <- paste0(labels, " (", seq_along(labels), ")")
  labels
}

.model_choice_role <- function(model) {
  model <- toupper(as.character(model %||% NA_character_)[1])
  switch(
    model,
    RSM = "equal_weighting_reference",
    PCM = "equal_weighting_reference",
    GPCM = "slope_aware_sensitivity",
    "ordered_response_candidate"
  )
}

.model_choice_score_contract <- function(model) {
  model <- toupper(as.character(model %||% NA_character_)[1])
  switch(
    model,
    RSM = "Common threshold structure; equal discrimination fixed at 1.",
    PCM = "Step thresholds vary by the designated step facet; equal discrimination fixed at 1.",
    GPCM = paste(
      "Positive slopes are estimated for the designated slope facet;",
      "the current bounded route requires slope_facet == step_facet and",
      "identifies slopes with geometric mean 1."
    ),
    "Ordered-response model; inspect the fitted configuration before reporting."
  )
}

.model_choice_report_template <- function(model) {
  model <- toupper(as.character(model %||% NA_character_)[1])
  switch(
    model,
    RSM = "We fit a many-facet rating-scale Rasch model, treating category thresholds as common across the step facet.",
    PCM = "We fit a many-facet partial-credit Rasch model, allowing step thresholds to vary by the designated step facet while retaining equal discrimination.",
    GPCM = "We fit a bounded generalized partial-credit many-facet model as a slope-aware sensitivity analysis; interpretation focused on whether discrimination-based reweighting changed the substantive conclusions.",
    "We fit a many-facet ordered-response model; report the fitted response-model contract explicitly."
  )
}

.model_choice_downstream_route <- function(model) {
  model <- toupper(as.character(model %||% NA_character_)[1])
  if (model %in% c("RSM", "PCM")) {
    return(tibble::tibble(
      FullAPARoute = "supported",
      ScoreSideExport = "supported",
      LinkingSynthesis = "supported",
      RecoveryChecks = "supported",
      FairAverage = "supported",
      BiasScreening = "supported",
      SummaryAppendix = "supported",
      PrimaryHelpers = paste(
        "diagnose_mfrm(..., diagnostic_mode = \"both\");",
        "compare_mfrm(); build_apa_outputs(); build_linking_review();",
        "fair_average_table(); estimate_bias()"
      )
    ))
  }
  if (identical(model, "GPCM")) {
    return(tibble::tibble(
      FullAPARoute = "blocked",
      ScoreSideExport = "blocked",
      LinkingSynthesis = "deferred",
      RecoveryChecks = "supported_with_caveat",
      FairAverage = "supported_with_caveat",
      BiasScreening = "supported_with_caveat",
      SummaryAppendix = "supported_with_caveat",
      PrimaryHelpers = paste(
        "gpcm_capability_matrix(); compare_mfrm(); build_weighting_review();",
        "compute_information(); evaluate_mfrm_recovery(); fair_average_table();",
        "estimate_bias(); export_summary_appendix()"
      )
    ))
  }
  tibble::tibble(
    FullAPARoute = "review_required",
    ScoreSideExport = "review_required",
    LinkingSynthesis = "review_required",
    RecoveryChecks = "review_required",
    FairAverage = "review_required",
    BiasScreening = "review_required",
    SummaryAppendix = "review_required",
    PrimaryHelpers = "summary(); diagnose_mfrm(); inspect the fitted configuration"
  )
}

.model_choice_role_table <- function(fits, labels) {
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    model <- toupper(as.character(fit$config$model %||% NA_character_)[1])
    method <- public_mfrm_method_label(as.character(fit$config$method %||% NA_character_)[1])
    step_facet <- as.character(fit$config$step_facet %||% NA_character_)[1]
    slope_facet <- as.character(fit$config$slope_facet %||% NA_character_)[1]
    tibble::tibble(
      Label = labels[[i]],
      Model = model,
      Method = method,
      StepFacet = step_facet,
      SlopeFacet = slope_facet,
      RecommendedRole = .model_choice_role(model),
      ScoreContract = .model_choice_score_contract(model),
      ReportTemplate = .model_choice_report_template(model)
    )
  })
  dplyr::bind_rows(rows)
}

.model_choice_downstream_table <- function(role_table) {
  rows <- lapply(seq_len(nrow(role_table)), function(i) {
    dplyr::bind_cols(
      role_table[i, c("Label", "Model", "RecommendedRole"), drop = FALSE],
      .model_choice_downstream_route(role_table$Model[[i]])
    )
  })
  dplyr::bind_rows(rows)
}

.model_choice_report_templates <- function(role_table) {
  role_table |>
    dplyr::transmute(
      Label = .data$Label,
      Model = .data$Model,
      UseFor = dplyr::case_when(
        .data$Model %in% c("RSM", "PCM") ~ "Operational equal-weighting report route",
        .data$Model == "GPCM" ~ "Slope-aware sensitivity report route",
        TRUE ~ "Model-specific report route"
      ),
      Template = .data$ReportTemplate,
      Avoid = dplyr::case_when(
        .data$Model == "GPCM" ~ "Do not write that better fit alone makes bounded GPCM the operational scoring model.",
        .data$Model %in% c("RSM", "PCM") ~ "Do not describe the fit as free-discrimination or slope-weighted.",
        TRUE ~ "Do not omit the model's response-kernel and score-contract assumptions."
      )
    )
}

.model_choice_route_map <- function(has_gpcm) {
  base <- tibble::tibble(
    Question = c(
      "Which model fits better on the same likelihood basis?",
      "Does bounded GPCM change the score interpretation?",
      "Where did precision move under bounded GPCM?",
      "Can I use manuscript APA/report bundles?",
      "Can I use fair averages?",
      "Can I screen bias/interactions?",
      "Can I export a summary appendix?"
    ),
    PrimaryHelper = c(
      "compare_mfrm()",
      "build_weighting_review()",
      "compute_information(); plot_information()",
      "build_apa_outputs(); build_visual_summaries()",
      "fair_average_table()",
      "estimate_bias()",
      "build_summary_table_bundle(); export_summary_appendix()"
    ),
    Interpretation = c(
      "Use AIC/BIC/logLik as fit evidence, not as a standalone scoring decision.",
      "Use only when comparing an RSM/PCM reference to bounded GPCM.",
      "Review item/rater/criterion information redistribution before changing the operational model.",
      "Supported for RSM/PCM; blocked for bounded GPCM in this release.",
      "Supported for RSM/PCM; supported with explicit SE caveat for bounded GPCM.",
      "Supported for RSM/PCM; conditional screening with profile-likelihood follow-up for bounded GPCM.",
      "Available for direct supported outputs; fit-based bundles remain RSM/PCM only."
    )
  )
  if (isTRUE(has_gpcm)) {
    base <- dplyr::bind_rows(
      tibble::tibble(
        Question = "What is the exact bounded GPCM support boundary?",
        PrimaryHelper = "gpcm_capability_matrix()",
        Interpretation = "Treat supported_with_caveat rows as review outputs and blocked/deferred rows as out of scope."
      ),
      base
    )
  }
  base
}

.model_choice_reference_labels <- function(role_table) {
  refs <- role_table$Label[role_table$Model %in% c("RSM", "PCM")]
  gpcm <- role_table$Label[role_table$Model == "GPCM"]
  list(
    reference = if (length(refs) > 0L) refs[[1]] else NA_character_,
    sensitivity = if (length(gpcm) > 0L) gpcm[[1]] else NA_character_
  )
}

.model_choice_find_weighting_pair <- function(fits, role_table) {
  ref_idx <- which(role_table$Model %in% c("RSM", "PCM"))[1]
  gpcm_idx <- which(role_table$Model == "GPCM")[1]
  if (is.na(ref_idx) || is.na(gpcm_idx)) return(NULL)
  list(reference = fits[[ref_idx]], gpcm = fits[[gpcm_idx]])
}

#' Build a model-choice review across RSM, PCM, and bounded GPCM fits
#'
#' @param ... Two or more fitted `mfrm_fit` objects from [fit_mfrm()].
#' @param labels Optional labels for the supplied fits. If omitted, names from
#'   `...` are used when available; otherwise labels are generated from
#'   model/method combinations.
#' @param run_weighting_review Logical. If `TRUE` and the supplied fits include
#'   at least one `RSM`/`PCM` reference plus one bounded `GPCM` fit, also run
#'   [build_weighting_review()] for the first such pair.
#' @param theta_range,theta_points,top_n Passed to [build_weighting_review()] when
#'   `run_weighting_review = TRUE`.
#' @param warn_constraints Passed to [compare_mfrm()].
#'
#' @details
#' `build_model_choice_review()` is a user-facing synthesis helper. It does not
#' estimate new models. It bundles:
#'
#' - [compare_mfrm()] for AIC/BIC/log-likelihood comparison;
#' - model-role guidance for `RSM`, `PCM`, and bounded `GPCM`;
#' - downstream-route availability for APA output, score-side export, linking,
#'   recovery, fair averages, bias screening, and summary-appendix handoff;
#' - report wording templates that avoid treating better bounded-`GPCM` fit as
#'   an automatic operational-scoring decision;
#' - [gpcm_capability_matrix()] when bounded `GPCM` is present;
#' - optionally, [build_weighting_review()] for the first Rasch-family reference
#'   versus bounded-`GPCM` pair.
#'
#' The word "bounded" is intentional: the package implements a bounded GPCM
#' route, not every possible generalized partial-credit many-facet extension.
#' The current route uses positive slopes, requires `slope_facet == step_facet`,
#' identifies slopes on the log scale with geometric mean 1, and keeps several
#' downstream score-side/reporting helpers outside the validated boundary.
#'
#' @return An object of class `mfrm_model_choice_review`.
#' @seealso [compare_mfrm()], [build_weighting_review()],
#'   [gpcm_capability_matrix()], [compute_information()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit_rsm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                     method = "MML", model = "RSM", quad_points = 7)
#' fit_pcm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                     method = "MML", model = "PCM", step_facet = "Criterion",
#'                     quad_points = 7)
#' review <- build_model_choice_review(RSM = fit_rsm, PCM = fit_pcm)
#' summary(review)
#' }
#' @export
build_model_choice_review <- function(...,
                                      labels = NULL,
                                      run_weighting_review = NULL,
                                      theta_range = c(-6, 6),
                                      theta_points = 61L,
                                      top_n = 10L,
                                      warn_constraints = TRUE) {
  run_weighting_review <- isTRUE(run_weighting_review)

  fits <- list(...)
  if (length(fits) < 2L) {
    stop("`build_model_choice_review()` requires at least two `mfrm_fit` objects.",
         call. = FALSE)
  }
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "mfrm_fit")) {
      stop("Argument ", i, " is not an `mfrm_fit` object. Got: ",
           paste(class(fits[[i]]), collapse = "/"), ".", call. = FALSE)
    }
  }

  labels <- .model_choice_fit_labels(fits, labels = labels)
  comparison <- suppressWarnings(do.call(
    compare_mfrm,
    c(fits, list(labels = labels, warn_constraints = warn_constraints, nested = FALSE))
  ))
  role_table <- .model_choice_role_table(fits, labels)
  downstream_routes <- .model_choice_downstream_table(role_table)
  report_templates <- .model_choice_report_templates(role_table)
  has_gpcm <- any(role_table$Model == "GPCM")
  ref_labels <- .model_choice_reference_labels(role_table)
  basis <- comparison$comparison_basis %||% list()

  overview <- tibble::tibble(
    FitCount = length(fits),
    Models = paste(role_table$Model, collapse = ", "),
    HasBoundedGPCM = has_gpcm,
    OperationalReference = ref_labels$reference,
    SensitivityModel = ref_labels$sensitivity,
    ICComparable = isTRUE(basis$ic_comparable),
    ReviewStatus = dplyr::case_when(
      has_gpcm && !is.na(ref_labels$reference) ~ "reference_plus_sensitivity_review",
      has_gpcm ~ "gpcm_without_rasch_family_reference",
      TRUE ~ "rasch_family_model_choice_review"
    )
  )

  key_warnings <- character(0)
  if (isTRUE(has_gpcm)) {
    key_warnings <- c(
      key_warnings,
      "Bounded GPCM is slope-aware: better fit is sensitivity evidence, not an automatic operational-scoring decision.",
      "Bounded means the current route is constrained to positive slopes, slope_facet == step_facet, and documented downstream support."
    )
  }
  if (!isTRUE(basis$ic_comparable)) {
    key_warnings <- c(
      key_warnings,
      "Information-criterion ranking is descriptive because the compared fits do not all share a comparable formal MML basis, observation set, and convergence status."
    )
  }
  if (isTRUE(has_gpcm) && is.na(ref_labels$reference)) {
    key_warnings <- c(
      key_warnings,
      "No RSM/PCM equal-weighting reference was supplied; add one before making an operational model-choice argument."
    )
  }
  key_warnings <- clean_summary_lines(key_warnings, max_n = 5L)

  next_actions <- clean_summary_lines(c(
    "Read comparison_table for fit evidence, but decide operational use from the score interpretation.",
    if (isTRUE(has_gpcm) && !is.na(ref_labels$reference)) {
      "Run with `run_weighting_review = TRUE` or call `build_weighting_review()` to inspect discrimination-based reweighting."
    },
    "Use downstream_routes before calling APA, score-side export, linking, recovery, fair-average, or bias-screening helpers.",
    "Use report_templates when drafting methods text so the fitted model and score contract are described accurately."
  ), max_n = 5L)

  support_status <- if (isTRUE(has_gpcm)) {
    gpcm_capability_matrix()
  } else {
    tibble::tibble()
  }

  weighting_review <- NULL
  weighting_review_error <- NA_character_
  if (isTRUE(run_weighting_review)) {
    pair <- .model_choice_find_weighting_pair(fits, role_table)
    if (is.null(pair)) {
      weighting_review_error <- "No RSM/PCM reference plus bounded GPCM pair was available."
    } else {
      weighting_review <- tryCatch(
        build_weighting_review(
          pair$reference,
          pair$gpcm,
          theta_range = theta_range,
          theta_points = theta_points,
          top_n = top_n
        ),
        error = function(e) {
          weighting_review_error <<- conditionMessage(e)
          NULL
        }
      )
    }
  }

  weighting_review_status <- tibble::tibble(
    Requested = isTRUE(run_weighting_review),
    Available = inherits(weighting_review, "mfrm_weighting_review"),
    Message = if (inherits(weighting_review, "mfrm_weighting_review")) {
      "Available in `weighting_review`."
    } else if (isTRUE(run_weighting_review)) {
      weighting_review_error
    } else {
      "Not requested; set `run_weighting_review = TRUE` for the first RSM/PCM versus bounded GPCM pair."
    }
  )

  notes <- clean_summary_lines(c(
    "This review is a decision aid; it does not refit models or choose an operational model automatically.",
    "Observation weights and GPCM discrimination-based reweighting are separate concepts.",
    "Use bounded GPCM wording only for the current constrained implementation, not for an unrestricted GPCM family claim."
  ))

  out <- list(
    overview = overview,
    key_warnings = key_warnings,
    next_actions = next_actions,
    comparison = comparison,
    comparison_table = comparison$table,
    model_roles = role_table,
    downstream_routes = downstream_routes,
    report_templates = report_templates,
    route_map = .model_choice_route_map(has_gpcm),
    support_status = support_status,
    weighting_review_status = weighting_review_status,
    weighting_review = weighting_review,
    notes = notes,
    settings = list(
      run_weighting_review = isTRUE(run_weighting_review),
      theta_range = theta_range,
      theta_points = theta_points,
      top_n = top_n,
      intended_use = "model_choice_review"
    )
  )
  as_mfrm_bundle(out, "mfrm_model_choice_review")
}

#' @export
print.mfrm_model_choice_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Summarize a model-choice review
#'
#' @param object Output from [build_model_choice_review()].
#' @param digits Number of digits for printed numeric values.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_model_choice_review`.
#' @seealso [build_model_choice_review()]
#' @export
summary.mfrm_model_choice_review <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_model_choice_review")) {
    stop("`object` must be output from build_model_choice_review().",
         call. = FALSE)
  }
  digits <- max(0L, as.integer(digits %||% 3L))
  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    key_warnings = clean_summary_lines(object$key_warnings %||% character(0)),
    next_actions = clean_summary_lines(object$next_actions %||% character(0)),
    comparison_table = tibble::as_tibble(object$comparison_table %||% tibble::tibble()),
    model_roles = tibble::as_tibble(object$model_roles %||% tibble::tibble()),
    downstream_routes = tibble::as_tibble(object$downstream_routes %||% tibble::tibble()),
    report_templates = tibble::as_tibble(object$report_templates %||% tibble::tibble()),
    route_map = tibble::as_tibble(object$route_map %||% tibble::tibble()),
    weighting_review_status = tibble::as_tibble(object$weighting_review_status %||% tibble::tibble()),
    support_status = tibble::as_tibble(object$support_status %||% tibble::tibble()),
    notes = clean_summary_lines(object$notes %||% character(0)),
    settings = object$settings %||% list(),
    digits = digits
  )
  class(out) <- "summary.mfrm_model_choice_review"
  out
}

#' @export
print.summary.mfrm_model_choice_review <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("mfrm Model Choice Review\n")
  if (nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  print_bullet_section("Key Warnings", x$key_warnings)
  print_bullet_section("Next Actions", x$next_actions)
  if (nrow(x$comparison_table) > 0L) {
    cat("\nComparison Table\n")
    print(round_numeric_df(as.data.frame(x$comparison_table), digits = digits), row.names = FALSE)
  }
  if (nrow(x$model_roles) > 0L) {
    cat("\nModel Roles\n")
    print(as.data.frame(x$model_roles[, c("Label", "Model", "RecommendedRole", "ScoreContract"), drop = FALSE]),
          row.names = FALSE)
  }
  if (nrow(x$downstream_routes) > 0L) {
    cat("\nDownstream Routes\n")
    keep <- intersect(
      c("Label", "Model", "FullAPARoute", "ScoreSideExport", "LinkingSynthesis",
        "RecoveryChecks", "FairAverage", "BiasScreening", "SummaryAppendix"),
      names(x$downstream_routes)
    )
    print(as.data.frame(x$downstream_routes[, keep, drop = FALSE]), row.names = FALSE)
  }
  if (nrow(x$weighting_review_status) > 0L) {
    cat("\nWeighting Review Status\n")
    print(as.data.frame(x$weighting_review_status), row.names = FALSE)
  }
  print_bullet_section("Notes", x$notes)
  invisible(x)
}
