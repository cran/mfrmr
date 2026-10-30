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
  if (!length(components) || length(components) != length(dfs)) return(NA_real_)
  ok <- is.finite(components) & components > 0 & is.finite(dfs) & dfs > 0
  if (!all(ok)) {
    return(NA_real_)
  }
  components <- components / max(components)
  num <- sum(components)^2
  den <- sum(components^2 / dfs)
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
  converged <- if (nrow(precision_profile) > 0 && "InferenceReady" %in% names(precision_profile)) {
    isTRUE(precision_profile$InferenceReady[1])
  } else {
    mfrm_inference_ready(sub_fit)
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
    tbl$Interpretation <- character(0)
    tbl$EffectMetric <- character(0)
    tbl$ContrastBasis <- character(0)
    tbl$SEBasis <- character(0)
    tbl$StatisticLabel <- character(0)
    tbl$ProbabilityMetric <- character(0)
    tbl$DFBasis <- character(0)
    tbl$ClassificationSystem <- character(0)
    tbl$Classification <- character(0)
    tbl$ETS <- character(0)
    tbl$ETS_Eligible <- logical(0)
    tbl$ConditionalRefitScreenEligible <- logical(0)
    tbl$FormalInferenceEligible <- logical(0)
    tbl$SupportsFormalInference <- logical(0)
    tbl$ReportingUse <- character(0)
    tbl$PrimaryReportingEligible <- logical(0)
    return(tbl)
  }

  if (identical(method, "refit")) {
    linked <- if ("ContrastComparable" %in% names(tbl)) {
      as.logical(tbl$ContrastComparable)
    } else {
      rep(TRUE, nrow(tbl))
    }
    conditional <- if ("ConditionalRefitScreenEligible" %in% names(tbl)) {
      as.logical(tbl$ConditionalRefitScreenEligible)
    } else {
      rep(FALSE, nrow(tbl))
    }
    tbl$ContrastBasis <- ifelse(
      linked,
      "linked subgroup facet-measure difference",
      "descriptive subgroup facet-measure difference"
    )
    tbl$SEBasis <- ifelse(
      conditional,
      paste0(
        "conditional independent-subgroup plug-in SE; baseline-anchor ",
        "uncertainty and cross-refit covariance omitted"
      ),
      ifelse(
        linked,
        "not reported without complete conditional-screen inputs",
        "not reported without common-scale linking"
      )
    )
    tbl$StatisticLabel <- ifelse(
      conditional,
      "conditional plug-in Welch t screen",
      ifelse(linked, "linked descriptive contrast", "descriptive contrast")
    )
    tbl$ProbabilityMetric <- ifelse(
      conditional,
      "conditional plug-in Welch tail area (screening only)",
      "not reported"
    )
    tbl$DFBasis <- ifelse(
      conditional,
      "Welch-Satterthwaite plug-in approximation; not a validated joint-refit df",
      "not reported"
    )
    tbl$EffectMetric <- ifelse(
      linked,
      "linked_descriptive_logit_difference",
      "descriptive_refit_difference"
    )
    tbl$ClassificationSystem <- "descriptive"
    tbl$Classification <- dplyr::case_when(
      linked %in% TRUE & is.finite(tbl$Contrast) ~ "Linked contrast (screening only)",
      linked %in% FALSE & is.finite(tbl$Contrast) ~ "Unclassified (insufficient linking)",
      TRUE ~ "Unavailable"
    )
    tbl$ETS <- NA_character_
    tbl$ETS_Eligible <- FALSE
    tbl$FormalInferenceEligible <- FALSE
    tbl$SupportsFormalInference <- FALSE
    tbl$PrimaryReportingEligible <- FALSE
    tbl$ReportingUse <- dplyr::case_when(
      !is.finite(tbl$Contrast) ~ "unavailable",
      conditional ~ "screening_only_conditional_plugin",
      linked %in% TRUE ~ "screening_only_linked_descriptive",
      TRUE ~ "screening_only"
    )
    if (all(c("Fixed1", "Fixed2") %in% names(tbl))) {
      fixed <- tbl$Fixed1 %in% TRUE | tbl$Fixed2 %in% TRUE
      if (all(c("Estimate1", "Estimate2") %in% names(tbl))) {
        fixed <- fixed & is.finite(tbl$Estimate1) & is.finite(tbl$Estimate2)
      }
      tbl$Classification[fixed] <- "Fixed target; no comparison"
      tbl$ContrastBasis[fixed] <- "Target fixed by constraints in at least one subgroup"
      tbl$SEBasis[fixed] <- "Not applicable to fixed targets"
      tbl$StatisticLabel[fixed] <- "No differential-functioning comparison"
      tbl$EffectMetric[fixed] <- "not_applicable"
      tbl$ReportingUse[fixed] <- "not_applicable_fixed_target"
    }
    return(tbl)
  }

  for (field in c("SE", "t", "df", "p_value", "p_adjusted")) {
    tbl[[field]] <- rep(NA_real_, nrow(tbl))
  }
  tbl$ContrastBasis <- "group difference in mean observed-minus-expected residuals"
  tbl$SEBasis <- "not reported: residual-contrast uncertainty is not established"
  tbl$StatisticLabel <- "residual mean difference"
  tbl$ProbabilityMetric <- "not reported"
  tbl$DFBasis <- "not applicable"
  tbl$EffectMetric <- "mean_obs_minus_exp_difference"
  tbl$ClassificationSystem <- "descriptive"
  tbl$Classification <- ifelse(is.finite(tbl$Contrast), "Residual contrast", "Unavailable")
  tbl$Interpretation <- "Residual differences do not isolate differential functioning; no test or confidence interval is provided."
  tbl$ETS <- NA_character_
  tbl$ETS_Eligible <- FALSE
  tbl$FormalInferenceEligible <- FALSE
  tbl$SupportsFormalInference <- FALSE
  tbl$PrimaryReportingEligible <- FALSE
  tbl$ReportingUse <- "descriptive"
  tbl
}

validate_dff_residual_output <- function(x) {
  cfg <- x$config %||% list()
  residual <- identical(cfg$method, "residual") || "abs_t_warn" %in% names(cfg)
  if (residual && !identical(cfg$residual_output_version, 3L)) {
    stop("This saved result uses earlier residual inference or observation matching. Recompute it with analyze_dff() or dif_interaction_table() using the fitted model; no model refit is needed.", call. = FALSE)
  }
  if (identical(cfg$method, "refit") && !identical(cfg$refit_output_version, 1L)) {
    stop("This saved refit result predates subgroup availability and cell-count checks. Rerun analyze_dff(method = 'refit') from the baseline fit and original data to recompute subgroup estimates; regenerate exported tables and figures.", call. = FALSE)
  }
  invisible(x)
}

build_dff_summary <- function(tbl, method) {
  if (identical(method, "refit")) {
    out <- tibble(
      Classification = c(
        "Linked contrast (screening only)",
        "Unclassified (insufficient linking)"
      ),
      Count = c(
        sum(tbl$Classification == "Linked contrast (screening only)", na.rm = TRUE),
        sum(tbl$Classification == "Unclassified (insufficient linking)", na.rm = TRUE)
      )
    )
    if (any(tbl$Classification == "Fixed target; no comparison", na.rm = TRUE)) {
      out <- dplyr::bind_rows(out, tibble(
        Classification = "Fixed target; no comparison",
        Count = sum(tbl$Classification == "Fixed target; no comparison", na.rm = TRUE)
      ))
    }
    if (any(tbl$Classification == "Unavailable", na.rm = TRUE)) {
      out <- dplyr::bind_rows(out, tibble(
        Classification = "Unavailable",
        Count = sum(tbl$Classification == "Unavailable", na.rm = TRUE)
      ))
    }
    return(out)
  }

  tibble(
    Classification = c("Residual contrast", "Unavailable"),
    Count = c(sum(is.finite(tbl$Contrast)), sum(!is.finite(tbl$Contrast)))
  )
}

resolve_dff_refit_controls <- function(fit) {
  replay <- fit$config$replay_inputs %||% list()
  control <- fit$config$estimation_control %||% list()
  if (is.list(replay$mml_starts))
    stop("Named starting coordinates cannot be reused after subgroup linking changes the free coordinates. Use numeric `mml_starts` offsets or fit the subgroup models explicitly.", call. = FALSE)

  population_active <- isTRUE(fit$config$population_active) ||
    isTRUE(fit$config$population_spec$active) ||
    !is.null(replay$population_formula)
  auto_gpcm_population <- population_active &&
    identical(as.character(fit$population$source %||% ""),
              "gpcm_mml_default_identification") &&
    identical(as.character(fit$config$model %||% ""), "GPCM") &&
    identical(as.character(fit$config$method %||% ""), "MML")
  if (population_active && !auto_gpcm_population) {
    stop(
      "`analyze_dff(method = \"refit\")` cannot currently reproduce an active ",
      "latent-regression population model within subgroup fits without changing ",
      "the fitted population contract. Use `method = \"residual\"` for screening, ",
      "or fit and link a prespecified subgroup population model outside ",
      "`analyze_dff()`.",
      call. = FALSE
    )
  }

  interaction_terms <- unique(c(
    as.character(fit$config$facet_interactions %||% character(0)),
    names(fit$config$interaction_specs %||% list())
  ))
  interaction_terms <- interaction_terms[!is.na(interaction_terms) & nzchar(interaction_terms)]
  if (length(interaction_terms) > 0L) {
    stop(
      "`analyze_dff(method = \"refit\")` cannot currently reproduce fitted ",
      "facet interactions and anchor their subgroup interaction coordinates on ",
      "a common scale. Use `method = \"residual\"` for screening, or refit and ",
      "link the interaction model explicitly outside `analyze_dff()`.",
      call. = FALSE
    )
  }

  anchor_tables <- tryCatch(
    extract_anchor_tables(fit$config),
    error = function(e) list(groups = data.frame())
  )
  group_anchors <- as.data.frame(
    replay$group_anchors %||% anchor_tables$groups %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(group_anchors) > 0L) {
    stop(
      "`analyze_dff(method = \"refit\")` cannot currently preserve a baseline ",
      "group-anchor constraint while constructing its separate subgroup-linking ",
      "anchors. Use `method = \"residual\"` for screening, or specify and audit ",
      "the subgroup linking model outside `analyze_dff()`.",
      call. = FALSE
    )
  }

  rating_min <- suppressWarnings(as.numeric(
    fit$config$rating_min %||% fit$prep$rating_min %||% NA_real_
  )[1])
  rating_max <- suppressWarnings(as.numeric(
    fit$config$rating_max %||% fit$prep$rating_max %||% NA_real_
  )[1])
  if (!is.finite(rating_min) || !is.finite(rating_max) || rating_min >= rating_max) {
    stop(
      "`analyze_dff(method = \"refit\")` requires the baseline fit's resolved ",
      "finite `rating_min` and `rating_max` so subgroup fits cannot redefine the ",
      "observed score scale. Refit the baseline model with the current mfrmr ",
      "version or use `method = \"residual\"`.",
      call. = FALSE
    )
  }

  list(
    model = fit$config$model %||% fit$summary$Model[1] %||% "RSM",
    method = fit$config$method %||% fit$summary$Method[1] %||% "JML",
    step_facet = fit$config$step_facet %||% NULL,
    slope_facet = fit$config$slope_facet %||% NULL,
    gpcm_mml_identification = as.character(
      replay$gpcm_mml_identification %||%
        fit$config$gpcm_mml_identification %||%
        if (identical(as.character(fit$config$model %||% ""), "GPCM") &&
            identical(as.character(fit$config$method %||% ""), "MML") &&
            !population_active) {
          "fixed_standard_normal"
        } else {
          "free_population"
        }
    )[1],
    rating_min = rating_min,
    rating_max = rating_max,
    weight = fit$config$weight_col %||% NULL,
    keep_original = isTRUE(replay$keep_original %||% fit$config$keep_original),
    repeated_ratings = replay$repeated_ratings %||% fit$config$repeated_ratings %||% "review",
    missing_codes = replay$missing_codes %||% NULL,
    noncenter_facet = fit$config$noncenter_facet %||% "Person",
    centering_scope = fit$config$centering_scope %||% "remainder",
    dummy_facets = fit$config$dummy_facets %||% NULL,
    positive_facets = fit$config$positive_facets %||% NULL,
    min_obs_per_element = as.numeric(replay$min_obs_per_element %||% 30),
    min_obs_per_category = as.numeric(replay$min_obs_per_category %||% 10),
    quad_points = as.integer(replay$quad_points %||% control$quad_points %||% 31L),
    mml_integration = as.character(control$mml_integration %||% "fixed"),
    mml_starts = replay$mml_starts,
    maxit = as.integer(replay$maxit %||% control$maxit %||% 400L),
    reltol = as.numeric(replay$reltol %||% control$reltol %||% 1e-9),
    optimizer = as.character(
      replay$optimizer %||% control$optimizer_requested %||% "auto"
    )[1],
    mml_engine = as.character(
      replay$mml_engine %||% control$mml_engine_requested %||% "direct"
    )[1],
    identity_contract = "baseline_response_and_estimation_controls_v1"
  )
}

build_dff_linking_setup <- function(fit, facet, facet_names) {
  linking_facets <- setdiff(as.character(facet_names), as.character(facet))
  anchor_tbl <- if (length(linking_facets) > 0) {
    make_anchor_table(
      fit,
      facets = linking_facets,
      include_person = FALSE,
      readiness_policy = "review"
    )
  } else {
    tibble::tibble(Facet = character(0), Level = character(0), Anchor = numeric(0))
  }
  # Linking fixes non-target estimates; it must not release supplied target
  # or JML Person anchors from the baseline model.
  original <- extract_anchor_tables(fit$config)$anchors
  if (nrow(original) > 0L) {
    retained <- original[!original$Facet %in% linking_facets,
                         c("Facet", "Level", "Anchor"), drop = FALSE]
    anchor_tbl <- dplyr::bind_rows(anchor_tbl, retained)
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
  # Keep the baseline target roster, including absent subgroup levels.
  sub_est <- dplyr::left_join(tibble(Level = as.character(fallback_levels)),
    sub_est, by = "Level")
  used <- sub_fit$prep$data
  counts <- used |>
    group_by(.data[[facet]]) |>
    summarise(N = n(), NPersons = n_distinct(.data$Person), .groups = "drop")
  names(counts)[1] <- "Level"
  counts$Level <- as.character(counts$Level)
  sub_est <- dplyr::left_join(sub_est, counts, by = "Level")
  sub_est$N <- dplyr::coalesce(sub_est$N, 0L)
  sub_est$NPersons <- dplyr::coalesce(sub_est$NPersons, 0L)
  absent <- sub_est$N == 0L
  sub_est$Estimate[absent] <- sub_est$SE[absent] <- NA_real_
  sub_est$Facet <- facet
  sub_est <- apply_fixed_measure_precision(sub_est, sub_fit$config)
  sub_est$Facet <- NULL
  sub_est$Fixed[absent] <- FALSE

  sub_est |>
    mutate(
      N_Subgroup = as.integer(nrow(used)),
      N_Input = as.integer(n_obs),
      EstimateStatus = ifelse(is.finite(.data$Estimate), "available", "unavailable"),
      EstimateReason = ifelse(absent, "Target level absent from retained subgroup responses.",
        ifelse(is.finite(.data$Estimate), NA_character_, "Subgroup estimate unavailable.")),
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
      ETS_Eligible = FALSE
    )
}

#' Differential facet functioning analysis
#'
#' Compares observed-minus-expected scores between groups, or describes linked
#' subgroup facet estimates after refitting. Residual differences do not isolate
#' differential functioning and are returned without tests or classifications.
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()]. With `method = "refit"`,
#'   it must match `fit`; see [mfrmr_compatibility_layer] for saved-diagnostic
#'   migration. The residual method does not use this argument.
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
#'   Residual analyses also require the fitted weight column when one was used.
#'   See the response-matching contract below for subsets and repeated rows.
#' @param focal Optional character vector of group levels to treat as focal.
#'   If `NULL` (default), all pairwise group comparisons are performed.
#' @param method Analysis method: `"residual"` (default) uses the fitted
#'   model's residuals without re-estimation; `"refit"` re-estimates the
#'   model within each group subset. The residual method is faster and
#'   avoids additional subgroup fits; neither method establishes adequate sample size.
#' @param min_obs Minimum number of observations per cell (facet-level x
#'   group). Cells below this threshold are flagged as sparse. Residual
#'   contrasts/scaled residuals and refit screening statistics are withheld;
#'   actual refit point estimates and contrasts remain available. Default `10`.
#' @param p_adjust Adjustment for refit screening tail areas; default `"holm"`.
#'   Retained but unused for residual comparisons, which return no p-values.
#' @param comparison_review Optional result of
#'   [describe_mfrm_data()] with `comparisons` specified, or its
#'   `$comparison_review` component. The retained response rows, origin groups,
#'   facet and score coding must match this request. If `focal` is supplied it
#'   must equal the review's first group. The prospective interaction states
#'   are retained separately; they do not classify the residual/refit contrasts
#'   as identified, change estimates, or permit an otherwise refused fit.
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
#' Differences between group ability distributions are not themselves DFF.
#' Residual screens inherit the fitted model's ability and population
#' assumptions. With fixed-standard-normal RSM/PCM MML, subgroup refits retain
#' that population assumption; linking anchors do not estimate subgroup ability
#' distributions. A residual difference can therefore reflect an inadequately
#' represented group difference as well as differential facet functioning.
#' For MML, observation expectations are evaluated at each Person's EAP
#' ability, rather than integrated over the conditional ability distribution.
#' This nonlinear substitution can produce different mean residuals across
#' groups even under a correctly specified no-DFF model. A zero mean-residual
#' contrast and absence of differential functioning are therefore distinct
#' null hypotheses; correcting the SE alone does not make them equivalent.
#'
#' Two methods are available:
#'
#' **Residual method** (`method = "residual"`): Uses the existing fitted
#' model's observation-level residuals.  For each facet-level \eqn{\times}
#' group cell, the observed and expected score sums are aggregated and
#' a standardized residual is computed as:
#' \deqn{z = \frac{\sum (X_{obs} - E_{exp})}{\sqrt{\sum \mathrm{Var}}}}
#' Pairwise contrasts compare the mean observed-minus-expected score in Group1
#' minus that in Group2. Positive values mean higher residual scores in Group1,
#' not greater rater leniency for that group. Differences use score units.
#' `StdResidual` is a descriptive scaling by model response variance, not a
#' t-statistic. No residual-contrast SE, confidence interval, p-value or
#' positive/negative classification is provided. Compatibility columns `SE`,
#' `t`, `df`, `p_value` and `p_adjusted` contain `NA`.
#'
#' @section Response matching for residual analyses:
#' `data` selects fitted response rows using Person, all fitted facets, the
#' observed score and fitted weight. Source column names, identifier trimming,
#' the saved original/internal category map and the saved missing-code policy
#' are respected. A changed score does not inherit a residual from another
#' response with the same Person and facets. Rows with missing groups, invalid
#' responses or no fitted match are excluded with reasons in `$cases`; unmatched
#' complete rows also produce a message. Requesting more copies of a response
#' than were fitted is an error. No requested match is an error.
#'
#' Repeated identical response rows are matched by their multiplicities, with
#' a canonical group order. They must have identical fitted expectations and
#' variances. `FitRow` then identifies an exchangeable fitted copy for this
#' calculation, not a recovered occasion/event identity. This permits descriptive
#' repeated-rating summaries without first-row group assignment or replication
#' of unselected rows. It does not supply cluster-adjusted inference.
#'
#' `$cases` contains `InputRow` (row of supplied data, or `fit$prep$data` when
#' `data = NULL`), `FitRow` (row of `fit$prep$data`), `GroupValue`, `Included`
#' and `Reason`. Preserve the data and fit alongside this audit. `N` and
#' `N_Group1`/`N_Group2` count retained response rows; `NPersons` and
#' `NPersons_Group1`/`NPersons_Group2` count distinct Persons. The same Person
#' may occur in multiple groups. Neither count establishes independence.
#' Aggregates remain unweighted row summaries even when fitting used weights.
#' Subsetting data retains the original model's fitted expectations; it does
#' not recalibrate the model on the subset.
#'
#' @details
#' **Refit method** (`method = "refit"`): Subsets the data by group, refits
#' the MFRM model within each subset, anchors all non-target facets back to
#' the baseline calibration when possible, and compares the resulting
#' facet-level estimates. When linking, convergence, model-based MML precision,
#' finite positive SEs, disjoint Person rosters and cell-count screens all pass,
#' a conditional plug-in Welch statistic is also
#' shown:
#' \deqn{t = \frac{\hat{\delta}_1 - \hat{\delta}_2}
#'                {\sqrt{SE_1^2 + SE_2^2}}}
#' This provides group-specific point estimates on a common scale when linking
#' anchors are available. However, the plug-in standard error conditions on
#' the baseline anchors and omits baseline-anchor uncertainty and cross-refit
#' covariance. Refit statistics therefore remain screening evidence and set
#' `FormalInferenceEligible`, `SupportsFormalInference`,
#' `PrimaryReportingEligible`, and `ETS_Eligible` to `FALSE`.
#' The Welch degrees of freedom use a response-count approximation, not
#' independently estimated degrees of freedom for model-based SEs. Disjoint
#' Person IDs do not establish independence or valid matching.
#'
#' `N_Group1` / `N_Group2` count retained responses for the particular target
#' level in each group; `N_Subgroup1` / `N_Subgroup2` count all retained subgroup
#' responses. Failed or skipped fits have unknown retained counts (`NA`), not
#' zero; `$group_fits[[group]]$N_Input` retains the input count. An absent target
#' in a successfully prepared subgroup has a known count of zero. All baseline
#' target levels remain in the result. Failed/skipped fits and absent levels
#' have missing estimates and contrasts, with `ComparisonStatus` and
#' `ComparisonReason`; they never inherit the baseline estimate.
#' Finite estimates from sparse cells or incomplete diagnostics are retained,
#' while `ConditionalScreenReason` explains withheld screening uncertainty.
#' `NPersons_Group1` / `NPersons_Group2` count Persons in each target cell,
#' and `NPersons_Overlap` counts shared Persons across the entire two subgroup
#' fits. Any overlap withholds the independent-subgroup screen.
#'
#' `ContrastDirection` describes the location parameter only:
#' `larger_location_in_group1`, `smaller_location_in_group1`, or
#' `equal_locations`. Translating this into higher/lower expected scores
#' requires the facet orientation and the remaining response parameters.
#' Subgroup step ladders or slopes can differ; a location difference alone
#' does not order their entire response curves.
#' The subgroup fits replay the baseline response family, resolved score range,
#' step/slope facets, weighting, repeated-rating independence declaration,
#' optimizer, MML engine, and numerical controls;
#' the non-target anchors are the intentional linking change. Direct anchors
#' on the target facet and JML Person anchors from the baseline are preserved.
#' A target fixed by constraints in either subgroup keeps its point estimates,
#' with `Fixed1` / `Fixed2` identifying the fixed values; its contrast, SE and
#' screening tests are unavailable. Equal supplied anchors do not demonstrate
#' absence of differential functioning. Use residual DFF to describe response
#' departures, or explicitly fit a baseline with different target constraints
#' before requesting free-parameter comparisons. Saved refit results created
#' before this preservation rule require rerunning `analyze_dff()` from the
#' baseline fit; loading an old result does not repair its subgroup estimates.
#' Refit analysis
#' fails closed when the baseline uses a user-specified active latent-regression
#' population model, facet interactions, or group-anchor constraints, because
#' those structures do not yet have a complete subgroup-replay and linking
#' contract. The automatically generated intercept-only population model used
#' for default GPCM-MML scale identification is replayed and is not treated as
#' a substantive covariate model.
#'
#' When `facet` refers to an item-like facet (for example `Criterion`), this
#' recovers the familiar DIF case. When `facet` refers to raters or
#' prompts/tasks, the same machinery supports DRF/DPF-style analyses.
#'
#' Refit contrasts are not assigned ETS A/B/C labels. Even when subgroup point
#' estimates share a linked logit scale, a validated joint, bootstrap, or
#' replicate covariance contract is required before formal refit inference.
#'
#' Refit screening tail areas are adjusted using Holm's step-down procedure by
#' default, jointly across the returned facet-level/group-pair rows in this
#' call. Holm's family-wise error control requires valid unadjusted p-values;
#' the approximate screening tail areas here have not been shown to meet that
#' requirement. Adjustment therefore does not establish an error-rate guarantee
#' or formal inference eligibility. Alternative adjustments can be specified
#' via `p_adjust`; see [stats::p.adjust()] for their assumptions.
#'
#' @section Choosing a method:
#' Use `method = "residual"` to describe group differences in model residuals
#' without fitting subgroup models. This does not test differential functioning.
#' Use `method = "refit"` when you specifically want group-specific parameter
#' estimates and can tolerate extra computation. Agreement between methods is
#' not guaranteed by a universal per-group sample-size threshold: stability and
#' detection depend jointly on effect size, category support, response-pattern
#' overlap, linking strength, slope heterogeneity, and estimator behavior.
#' `min_obs` is only a cell-computability/sparsity guard; it is not evidence of
#' power, parameter stability, or sample-size adequacy.
#'
#' @section Comparing with MH or ordinal logistic DIF:
#' [mantel_haenszel_table()] supplies a separate binary conditional comparison
#' with explicit matching strata and one response per Person-target cell.
#' [ordinal_logistic_table()] supplies a separate cumulative-logit comparison
#' with uniform/nonuniform terms and a proportional-odds diagnostic. Neither
#' is a `method` of this function. A comparison needs the same persons, target
#' response, group orientation and defensible matching variable. With a common
#' binary Rasch item set and invariant other items, the total including the
#' studied item removes the ability distribution from the conditional item
#' response; a rest score generally does not. This result does not justify
#' pooling different rater/task plans at the same raw total.
#'
#' Ordinary MH targets a common conditional odds ratio and can miss opposing
#' effects across matching strata. Ordinal logistic models distinguish group
#' and matching-by-group terms, but their cumulative proportional-odds
#' assumption is not generally implied by RSM/PCM. Repeated ratings from one
#' Person cannot be counted as independent persons for either ordinary method.
#' Agreement or disagreement with these procedures does not by itself validate
#' a DFF conclusion. See the linking-and-DFF vignette for the comparison scope
#' and primary references.
#'
#' @section Interpreting output:
#' - `$dif_table`: one row per facet-level x group-pair with contrast,
#'   a residual mean difference or linked subgroup difference. Residual rows have
#'   no SE, test or binary classification. Includes `Method`, `N_Group1`,
#'   `N_Group2`, `EffectMetric`, `ClassificationSystem`, `ContrastBasis`,
#'   `SEBasis`, `StatisticLabel`, `ProbabilityMetric`, `DFBasis`,
#'   `ReportingUse`, `PrimaryReportingEligible`, and `sparse` columns.
#' - `$cell_table`: (residual method only) per-cell detail with N, NPersons,
#'   ObsScore, ExpScore, ObsExpAvg, StdResidual and an interpretation note.
#' - `$cases`: (residual method only) requested-row matching and exclusion audit,
#'   also retained by `summary()`.
#' - `$summary`: available/unavailable residual comparisons or linked-
#'   screening, insufficient-linking, fixed-target and unavailable rows
#'   (`method = "refit"`).
#' - `$group_fits`: (refit method only) list of per-group facet estimates and
#'   subgroup linking diagnostics, with `EstimateStatus` / `EstimateReason`
#'   and cell, subgroup and input counts.
#' - `$gpcm_boundary`: for `GPCM` fits, a capability-boundary table
#'   marking the DFF/DIF output as caveated screening evidence.
#' - `$comparison_review`: when supplied, the separate prospective interaction
#'   question, weights, counts, assumptions and identification states. Preserved
#'   by `summary()`, [dif_report()] and both DFF plot views. Its model can differ
#'   from the fitted screening model; read its `fit_arguments` and `notes`.
#'
#' @section GPCM boundary:
#' For `GPCM`, DFF/DIF rows are available as slope-aware screening
#' evidence over the fitted expected-score and residual scale. Keep
#' residual-method contrasts and interaction cells in screening language.
#' Refit contrasts require explicit subgroup linking and precision support for
#' conditional screening, but remain in screening language.
#'
#' @section Saved results:
#' Residual results, summaries and reports saved by earlier versions must be
#' recomputed with this function using the fitted model and original data.
#' This does not require refitting the model. Previously exported tables and
#' figures should also be regenerated.
#' Refit results created before the subgroup availability/cell-count checks
#' must instead be recomputed with `method = "refit"`. This reruns the subgroup
#' fits from the saved baseline fit and original data. Reading an old result
#' cannot repair borrowed estimates, dropped target levels or incorrect counts;
#' summary, report and plotting methods reject these earlier refit objects.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()]. For `RSM` / `PCM` fairness review, prefer
#'    `method = "MML"`.
#' 2. Run [diagnose_mfrm()] and, for `RSM` / `PCM`, prefer
#'    `diagnostic_mode = "both"` so legacy and strict marginal screens remain
#'    visible together.
#' 3. Run `analyze_dff(fit, diagnostics, facet = "Criterion", group = "Gender", data = my_data)`.
#' 4. Inspect `$dif_table` for contrasts and their reporting status, and
#'    `$summary` for method-specific counts. Residual rows have no flags.
#' 5. Use [dif_interaction_table()] when you need cell-level diagnostics.
#' 6. Use [plot(x, type = "heatmap")][plot.mfrm_dif] or [dif_report()] for communication.
#'
#' @return
#' An object of class `mfrm_dff` (with compatibility class `mfrm_dif`) with:
#' - `dif_table`: data.frame of differential-functioning contrasts. Refit rows
#'   include `Fixed1` and `Fixed2`; fixed-target rows have
#'   `ReportingUse = "not_applicable_fixed_target"` and no contrast or test.
#' - `cell_table`: (residual method) per-cell detail table.
#' - `summary`: counts by method-appropriate screening classification.
#' - `group_fits`: (refit method) per-group facet estimates.
#' - `gpcm_boundary`: for `GPCM` fits, a capability-boundary table.
#' - `config`: list with facet, group, method, min_obs, p_adjust settings.
#'
#' @seealso [fit_mfrm()], [estimate_bias()], [compare_mfrm()],
#'   [dif_interaction_table()], [plot(x, type = "heatmap")][plot.mfrm_dif], [dif_report()],
#'   [subset_connectivity_report()], [mfrmr_linking_and_dff]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("analyze_dff", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", model = "RSM", quad_points = 7, maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#' dff <- analyze_dff(fit, diag, facet = "Rater", group = "Group", data = toy)
#' dff$summary
#' # Read the residual mean difference and the number of observations per group.
#' head(dff$dif_table[, c("Level", "Group1", "Group2", "Contrast",
#'                        "N_Group1", "N_Group2")])
#' # The residual contrast is an observed-minus-expected average contrast
#' # between groups. It is useful for screening, but it is not an ETS
#' # A/B/C logit-delta classification.
#' dff_refit <- analyze_dff(fit, diag, facet = "Rater", group = "Group",
#'                          data = toy, method = "refit")
#' unique(dff_refit$dif_table$ClassificationSystem)
#' # Look for: "descriptive". Linked refit contrasts remain screening-only
#' #   because their plug-in uncertainty omits anchor uncertainty and
#' #   cross-refit covariance.
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
                        p_adjust = "holm",
                        comparison_review = NULL) {
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

  comparison_review <- mfrmr_bind_comparison_review(comparison_review, fit,
    group_info$input_data, facet, group, focal)

  if (method == "residual") {
    out <- .analyze_dif_residual(
      fit = fit, facet = facet, group = group, data = data,
      orig_data = group_info$input_data, facet_names = facet_names,
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
  if (!is.null(comparison_review)) out$comparison_review <- comparison_review
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
  input_data <- orig_data
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
    input_data = input_data,
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
    "Group residual mean difference (score units)"
  } else if (identical(method, "refit")) {
    "Subgroup location difference (logits)"
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
      Meaning = "Zero denotes no difference on the displayed scale.",
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
      ReportingNote = "This scaling is descriptive; it does not supply a t test or a differential-functioning decision.",
      stringsAsFactors = FALSE
    )
  } else if (identical(metric, "contrast") || identical(metric, "summary")) {
    rows[[length(rows) + 1L]] <- data.frame(
      Item = "Group comparison",
      Meaning = if (identical(method, "residual")) {
        "Group1 minus Group2 in mean observed-minus-expected scores (score units): a positive contrast means higher residual scores in Group1, not established leniency or favoritism."
      } else {
        "Group1 minus Group2 in facet locations (logits), when linking supports a common scale: a positive contrast means a larger Group1 location. Holding other parameters fixed, the usual difficulty/severity orientation lowers expected scores; subgroup step or slope differences can prevent ordering the whole response curves."
      },
      ReportingNote = if (identical(method, "residual")) {
        "Residual differences do not isolate differential functioning; no test or confidence interval is provided."
      } else {
        "Finite values do not certify a common scale; check ScaleLinkStatus. Linked refit differences condition on the estimated anchors; intervals omit their uncertainty and cross-refit covariance. Formal DFF inference remains unavailable."
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
  matched <- .dff_match_observations(fit, orig_data, group, normalized = is.null(data))
  obs_work <- matched$data

  # Aggregate by facet level x group
  cell_table <- obs_work |>
    group_by(.data[[facet]], .data$.group_var) |>
    summarise(
      N = n(),
      NPersons = n_distinct(.data$Person),
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
      t = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      FormalInferenceEligible = FALSE,
      ReportingUse = "descriptive",
      Interpretation = "Observed-minus-expected scores; scaled residuals are descriptive and are not t statistics."
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
      abs_diff <- abs(contrast)
      direction <- if (!is.finite(contrast)) NA_character_ else if (contrast > 0) {
        "higher_residual_in_group1"
      } else if (contrast < 0) "lower_residual_in_group1" else "equal_residuals"
      dif_rows[[length(dif_rows) + 1]] <- tibble(
        Level = lev,
        Group1 = g1,
        Group2 = g2,
        Contrast = contrast,
        ContrastDirection = direction,
        SE = NA_real_,
        t = NA_real_,
        df = NA_real_,
        p_value = NA_real_,
        AbsDiff = abs_diff,
        Method = "residual",
        N_Group1 = as.integer(n1),
        N_Group2 = as.integer(n2),
        NPersons_Group1 = if (nrow(c1)) c1$NPersons[1] else 0L,
        NPersons_Group2 = if (nrow(c2)) c2$NPersons[1] else 0L,
        sparse = is_sparse,
        ContrastComparable = FALSE,
        FormalInferenceEligible = FALSE,
        PrimaryReportingEligible = FALSE,
        InferenceTier = "descriptive",
        ComparisonMethod = method_label,
        ScaleLinkStatus = "not_applicable",
        ReportingUse = "descriptive"
      )
    }
  }
  dif_table <- bind_rows(dif_rows)

  dif_table <- annotate_dff_table(dif_table, method = "residual")

  # Summary counts
  dif_summary <- build_dff_summary(dif_table, method = "residual")
  functioning_label <- functioning_label_for_facet(facet)

  out <- list(
    dif_table = dif_table,
    cell_table = cell_table,
    cases = matched$cases,
    summary = dif_summary,
    group_fits = NULL,
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "analyze_dff()",
      area = "Differential facet functioning screening under GPCM"
    ),
    config = list(facet = facet, group = group, method = "residual",
                  residual_output_version = 3L,
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
  mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "analyze_dff(method = 'refit')")
  # Get full-sample facet estimates
  measures <- tibble::as_tibble(diagnostics$measures)
  facet_estimates <- measures |>
    filter(.data$Facet == facet) |>
    select("Facet", "Level", "Estimate", "SE")
  facet_estimates <- apply_fixed_measure_precision(facet_estimates, fit$config)
  facet_estimates$Facet <- NULL

  person_col <- fit$config$person_col %||% "Person"
  score_col <- fit$config$score_col %||% "Score"
  refit_controls <- resolve_dff_refit_controls(fit)
  linking_setup <- build_dff_linking_setup(fit, facet = facet, facet_names = facet_names)
  baseline_precision_meta <- resolve_dff_subgroup_precision(fit, diagnostics)

  group_fits <- list()
  group_persons <- list()
  for (g in group_levels) {
    idx <- group_vals == g
    sub_data <- orig_data[idx, , drop = FALSE]
    if (nrow(sub_data) < 5) {
      group_fits[[g]] <- facet_estimates |>
        mutate(
          Estimate = NA_real_, SE = NA_real_, Fixed = FALSE,
          N = NA_integer_, NPersons = NA_integer_, N_Subgroup = NA_integer_,
          N_Input = as.integer(nrow(sub_data)),
          EstimateStatus = "unavailable",
          EstimateReason = "Subgroup had fewer than 5 input observations; refit was skipped.",
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
      rating_min = refit_controls$rating_min,
      rating_max = refit_controls$rating_max,
      weight = refit_controls$weight,
      keep_original = refit_controls$keep_original,
      repeated_ratings = refit_controls$repeated_ratings,
      missing_codes = refit_controls$missing_codes,
      method = refit_controls$method,
      model = refit_controls$model,
      step_facet = refit_controls$step_facet,
      slope_facet = refit_controls$slope_facet,
      gpcm_mml_identification = refit_controls$gpcm_mml_identification,
      anchors = if (nrow(linking_setup$anchor_tbl) > 0) linking_setup$anchor_tbl else NULL,
      noncenter_facet = refit_controls$noncenter_facet,
      centering_scope = refit_controls$centering_scope,
      dummy_facets = refit_controls$dummy_facets,
      positive_facets = refit_controls$positive_facets,
      anchor_policy = "silent",
      min_common_anchors = linking_setup$min_common_anchors,
      min_obs_per_element = refit_controls$min_obs_per_element,
      min_obs_per_category = refit_controls$min_obs_per_category,
      quad_points = refit_controls$quad_points,
      maxit = refit_controls$maxit,
      reltol = refit_controls$reltol,
      optimizer = refit_controls$optimizer,
      mml_engine = refit_controls$mml_engine,
      mml_starts = refit_controls$mml_starts,
      mml_integration = refit_controls$mml_integration
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
          Estimate = NA_real_, SE = NA_real_, Fixed = FALSE,
          N = NA_integer_, NPersons = NA_integer_, N_Subgroup = NA_integer_,
          N_Input = as.integer(nrow(sub_data)),
          EstimateStatus = "unavailable",
          EstimateReason = paste0("Subgroup refit failed: ", sub_fit$message),
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
      group_persons[[g]] <- unique(as.character(sub_fit$prep$data$Person))
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
    overlap <- if (g1 %in% names(group_persons) && g2 %in% names(group_persons)) {
      length(intersect(group_persons[[g1]], group_persons[[g2]]))
    } else NA_integer_
    merged <- merge(est1, est2, by = "Level", suffixes = c("_1", "_2"))
    for (j in seq_len(nrow(merged))) {
      e1 <- merged$Estimate_1[j]
      e2 <- merged$Estimate_2[j]
      se1 <- merged$SE_1[j]
      se2 <- merged$SE_2[j]
      n1 <- merged$N_1[j]
      n2 <- merged$N_2[j]
      fixed1 <- isTRUE(merged$Fixed_1[j])
      fixed2 <- isTRUE(merged$Fixed_2[j])
      fixed_comparison <- (fixed1 || fixed2) && is.finite(e1) && is.finite(e2)
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
      is_sparse <- if (is.finite(n1) && is.finite(n2)) {
        (n1 < min_obs) || (n2 < min_obs)
      } else NA
      points_available <- is.finite(e1) && is.finite(e2) && is.finite(e1 - e2)
      se_scale <- max(se1, se2)
      se_diff <- if (is.finite(se1) && is.finite(se2) && se1 > 0 && se2 > 0) {
        se_scale * sqrt((se1 / se_scale)^2 + (se2 / se_scale)^2)
      } else NA_real_
      complete_precision <- points_available && is.finite(se_diff) && se_diff > 0 &&
        is.finite((e1 - e2) / se_diff) &&
        is.finite(n1) && is.finite(n2) && n1 > 1 && n2 > 1
      conditional_screen_eligible <- link_comparable &&
        subgroup_formal &&
        subgroup_converged &&
        isTRUE(baseline_precision_meta$supports_formal) &&
        isTRUE(baseline_precision_meta$converged) &&
        identical(comparison_method, "MML") &&
        isFALSE(is_sparse) && !fixed_comparison && complete_precision &&
        isTRUE(overlap == 0L)
      comparison_reason <- if (fixed_comparison) {
        "Target fixed by constraints in at least one subgroup."
      } else if (!points_available) {
        reasons <- c(merged$EstimateReason_1[j], merged$EstimateReason_2[j])
        reasons <- unique(reasons[!is.na(reasons) & nzchar(reasons)])
        if (length(reasons)) paste(reasons, collapse = " | ") else "Nonfinite subgroup location or difference."
      } else NA_character_
      screen_reason <- if (fixed_comparison || !points_available) {
        comparison_reason
      } else if (!complete_precision) {
        "Finite positive subgroup SEs, a representable standardized difference and at least two retained responses per cell are required."
      } else if (isTRUE(overlap > 0L)) {
        "Subgroups contain shared Persons; cross-refit covariance is unavailable."
      } else if (!isFALSE(is_sparse)) {
        "Target-by-group cell is sparse or its retained count is unavailable."
      } else if (!conditional_screen_eligible) {
        "Linking, baseline/subgroup readiness or model-based MML precision requirements are not met."
      } else NA_character_
      formal_eligible <- FALSE
      reporting_use <- dplyr::case_when(
        conditional_screen_eligible ~ "screening_only_conditional_plugin",
        link_comparable ~ "screening_only_linked_descriptive",
        TRUE ~ "screening_only"
      )
      contrast <- if (fixed_comparison || !points_available) NA_real_ else e1 - e2
      if (!conditional_screen_eligible) se_diff <- NA_real_
      t_val <- if (is.finite(se_diff) && se_diff > 0) contrast / se_diff else NA_real_
      df_welch <- if (conditional_screen_eligible && is.finite(se1) && is.finite(se2) && se1 > 0 && se2 > 0) {
        welch_satterthwaite_df(c((se1 / se_scale)^2, (se2 / se_scale)^2),
          c(n1 - 1, n2 - 1))
      } else {
        NA_real_
      }
      p_val <- if (is.finite(t_val) && is.finite(df_welch) && df_welch > 0) {
        2 * stats::pt(abs(t_val), df = df_welch, lower.tail = FALSE)
      } else {
        NA_real_
      }
      abs_diff <- abs(contrast)
      # Location alone does not order response curves with different steps/slopes.
      direction <- if (is.finite(contrast)) {
        if (contrast > 0) "larger_location_in_group1" else if (contrast < 0) {
          "smaller_location_in_group1"
        } else "equal_locations"
      } else {
        NA_character_
      }
      dif_rows[[length(dif_rows) + 1]] <- tibble(
        Level = merged$Level[j],
        Group1 = g1,
        Group2 = g2,
        Estimate1 = e1,
        Estimate2 = e2,
        Fixed1 = fixed1,
        Fixed2 = fixed2,
        Contrast = contrast,
        ContrastDirection = direction,
        ComparisonStatus = if (fixed_comparison) "fixed_target" else if (points_available) "available" else "unavailable",
        ComparisonReason = comparison_reason,
        ConditionalScreenReason = screen_reason,
        SE = se_diff,
        t = t_val,
        df = df_welch,
        p_value = p_val,
        AbsDiff = abs_diff,
        Method = "refit",
        N_Group1 = as.integer(n1),
        N_Group2 = as.integer(n2),
        NPersons_Group1 = merged$NPersons_1[j],
        NPersons_Group2 = merged$NPersons_2[j],
        N_Subgroup1 = merged$N_Subgroup_1[j],
        N_Subgroup2 = merged$N_Subgroup_2[j],
        NPersons_Overlap = as.integer(overlap),
        sparse = is_sparse,
        ContrastComparable = link_comparable,
        ConditionalRefitScreenEligible = conditional_screen_eligible,
        FormalInferenceEligible = formal_eligible,
        PrimaryReportingEligible = formal_eligible,
        SupportsFormalInference = formal_eligible,
        InferenceTier = inference_tier,
        ComparisonMethod = comparison_method,
        ReportingUse = reporting_use,
        ETS_Eligible = FALSE,
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
    dif_table$p_adjusted <- stats::p.adjust(dif_table$p_value, method = p_adjust,
      n = nrow(dif_table))
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
      area = "Differential facet functioning screening under GPCM"
    ),
    config = list(facet = facet, group = group, method = "refit",
                  refit_output_version = 1L,
                  min_obs = min_obs, p_adjust = p_adjust,
                  focal = focal, group_levels = group_levels,
                  linking_facets = linking_setup$linking_facets,
                  linking_threshold = linking_setup$min_common_anchors,
                  refit_identity_contract = refit_controls$identity_contract,
                  refit_anchor_policy = "preserve_baseline_direct_anchors",
                  refit_rating_range = c(refit_controls$rating_min,
                                         refit_controls$rating_max),
                  refit_slope_facet = refit_controls$slope_facet,
                  refit_optimizer = refit_controls$optimizer,
                  refit_mml_engine = refit_controls$mml_engine,
                  refit_mml_integration = refit_controls$mml_integration,
                  functioning_label = functioning_label)
  )
  class(out) <- c("mfrm_dff", "mfrm_dif", class(out))
  out
}

#' Review differential-functioning results
#'
#' `summary()` returns the contrast, cell, summary, and configuration tables
#' without refitting. A summary from [analyze_dff()] keeps its DFF label when
#' printed; a summary from [analyze_dif()] keeps its DIF label.
#'
#' @param object,x An `mfrm_dff` or `mfrm_dif` result, or its summary.
#' @param ... Reserved for generic compatibility.
#' @return `summary()` returns a structured differential-functioning summary.
#'   `print()` returns its input invisibly.
#' @seealso [analyze_dff()], [analyze_dif()]
#' @name mfrm_dff_methods
NULL

#' @rdname mfrm_dff_methods
#' @method summary mfrm_dif
#' @export
summary.mfrm_dif <- function(object, ...) {
  validate_dff_residual_output(object)
  out <- list(
    dif_table = object$dif_table,
    cell_table = object$cell_table,
    cases = object$cases,
    comparison_review = object$comparison_review,
    summary = object$summary,
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
    config = object$config
  )
  class(out) <- "summary.mfrm_dif"
  out
}

#' @rdname mfrm_dff_methods
#' @method summary mfrm_dff
#' @export
summary.mfrm_dff <- function(object, ...) {
  out <- summary.mfrm_dif(object, ...)
  class(out) <- c("summary.mfrm_dff", class(out))
  out
}

#' @rdname mfrm_dff_methods
#' @method print summary.mfrm_dif
#' @export
print.summary.mfrm_dif <- function(x, ...) {
  validate_dff_residual_output(x)
  label <- x$config$functioning_label %||% "DFF"
  cat(if (identical(x$config$method, "residual")) "--- Group residual comparisons ---\n" else paste0("--- ", label, " Analysis ---\n"))
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
        "N_Group1", "N_Group2", "sparse"),
      names(x$dif_table)
    )
    if (identical(x$config$method, "residual")) {
      show_cols <- intersect(c("Level", "Group1", "Group2", "Contrast", "N_Group1", "N_Group2", "sparse"), names(x$dif_table))
    }
    print_tbl <- x$dif_table |> select(all_of(show_cols))
    print(as.data.frame(print_tbl), row.names = FALSE, digits = 3)
  } else {
    cat("No differential-functioning contrasts computed.\n")
  }

  if (identical(x$config$method, "refit")) {
    cat("\nRefit Classification Summary:\n")
  } else {
    cat("\nResidual differences do not isolate differential functioning; no test or confidence interval is provided.\n\nAvailable comparisons:\n")
  }
  print(as.data.frame(x$summary), row.names = FALSE)
  mfrmr_print_comparison_review(x$comparison_review)
  .print_dff_gpcm_boundary(x$gpcm_boundary)
  invisible(x)
}

#' @rdname mfrm_dff_methods
#' @method print summary.mfrm_dff
#' @export
print.summary.mfrm_dff <- function(x, ...) {
  print.summary.mfrm_dif(x, ...)
  invisible(x)
}

#' @rdname mfrm_dff_methods
#' @method print mfrm_dif
#' @export
print.mfrm_dif <- function(x, ...) {
  .print_dif_overview(x, label = "DIF")
  invisible(x)
}

#' @rdname mfrm_dff_methods
#' @method print mfrm_dff
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
  validate_dff_residual_output(x)
  cfg <- x$config %||% list()
  method <- cfg$method %||% "refit"
  facet <- cfg$facet %||% NA_character_
  group <- cfg$group %||% NA_character_
  groups <- cfg$group_levels %||% character(0)
  tbl <- x$dif_table
  n_rows <- if (is.data.frame(tbl)) nrow(tbl) else 0L
  cat(if (identical(method, "residual")) "Group residual comparisons\n" else sprintf("mfrm_%s (%s)\n", tolower(label), label))
  cat(sprintf("  Method: %s | Facet: %s | Group: %s\n",
              method, facet, group))
  if (length(groups) > 0L) {
    cat(sprintf("  Group levels: %s\n",
                paste(groups, collapse = ", ")))
  }
  if (n_rows > 0L) {
    cat(sprintf("  Contrasts: %d row(s)", n_rows))
    cat("\n")
  } else {
    cat("  Contrasts: 0 row(s)\n")
  }
  boundary <- x$gpcm_boundary %||% data.frame()
  if (is.data.frame(boundary) && nrow(boundary) > 0L) {
    cat("  GPCM: differences use the fitted category slopes and expected scores.\n")
  }
  if (identical(method, "residual")) cat("  Residual differences do not isolate differential functioning; no test or confidence interval is provided.\n")
  cat("  Use `summary()` for the contrast table.\n")
  invisible(NULL)
}

.print_dff_gpcm_boundary <- function(boundary) {
  if (is.null(boundary) || !is.data.frame(boundary) || nrow(boundary) == 0L) {
    return(invisible(NULL))
  }
  cat("\nGPCM: these summaries use the fitted category slopes and expected scores.\n")
  invisible(NULL)
}

# ============================================================================
# B2. Differential Functioning Interaction Table
# ============================================================================

#' Compute interaction table between a facet and a grouping variable
#'
#' Produces a cell-level interaction table showing Obs-Exp differences,
#' scaled residuals, and absolute residual mean comparisons for each
#' facet-level x group-value cell.
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet Character scalar naming the facet.
#' @param group Character scalar naming the grouping column.
#' @param data Optional data frame with the group column. If `NULL`
#'   (default), the data stored in `fit$prep$data` is used, but it
#'   must contain the `group` column.
#'   All fitted response columns, including the weight column if used, are
#'   required. Subsets and repeated responses follow [analyze_dff()]'s
#'   response-matching contract; unmatched rows remain visible in `$cases`.
#' @param min_obs Minimum observations per cell. Cells with fewer than
#'   this many observations are marked sparse; scaled residuals and magnitude
#'   flags are `NA`. Observed and expected score summaries remain available.
#' @param p_adjust Retained for compatibility; unused because no p-values are
#'   reported. Default `"holm"`.
#' @param abs_t_warn Retained for compatibility; unused because scaled residuals
#'   are not t statistics. `flag_t` is `NA`.
#' @param abs_bias_warn Threshold for marking the absolute observed-minus-expected
#'   average, in score units. Default `0.5`. This is a descriptive magnitude rule.
#'
#' @details
#' This function uses observation-level residuals computed from the fitted
#' model rather than re-estimating it. For each facet-level x group-value cell,
#' it computes:
#' \itemize{
#'   \item N: number of observations in the cell
#'   \item NPersons: number of distinct Persons in the cell, not an independence certificate
#'   \item ObsScore: sum of observed scores
#'   \item ExpScore: sum of expected scores
#'   \item ObsExpAvg: mean observed-minus-expected difference
#'   \item Var_sum: sum of model variances
#'   \item StdResidual: (ObsScore - ExpScore) / sqrt(Var_sum)
#'   \item t, df, p_value, p_adjusted, flag_t: `NA`, retained for compatibility
#' }
#'
#' @section When to use this instead of analyze_dff():
#' Use `dif_interaction_table()` when you want cell-level screening for a
#' single facet-by-group table. Use [analyze_dff()] when you want group-pair
#' comparisons. Neither residual output tests differential functioning.
#'
#' @section Further guidance:
#' For plot selection and follow-up diagnostics, see
#' [mfrmr_visual_diagnostics].
#'
#' @section Interpreting output:
#' - `$table`: the full interaction table with one row per cell.
#' - `$cases`: requested-row matching and exclusion audit, also kept by `summary()`.
#' - `$summary`: overview counts of flagged and sparse cells.
#' - `$config`: analysis configuration parameters.
#' - `$gpcm_boundary`: for `GPCM` fits, a capability-boundary table
#'   marking the table as caveated DFF screening evidence.
#' - `flag_bias` records `|ObsExpAvg| > abs_bias_warn` in score units.
#'   It does not establish differential functioning. `flag_t` is unavailable.
#' - Sparse cells (N < min_obs) have `sparse = TRUE` and unavailable scaled
#'   residuals and magnitude flags. The score means and counts are retained.
#'
#' @section GPCM boundary:
#' For `GPCM`, the interaction table uses the fitted slope-aware
#' expected-score/residual scale and should be reported as screening evidence,
#' not as a standalone fairness, invariance, or operational subgroup decision.
#'
#' Residual means can differ even when response parameters are the same between
#' groups. These summaries provide no p-values or t-based decisions. Recompute
#' older saved tables with this function using the fitted model; no refit is needed.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Run `dif_interaction_table(fit, diag, facet = "Rater", group = "Gender", data = df)`.
#' 3. Inspect `$table` for flagged cells.
#' 4. Visualize with [plot(x, type = "heatmap")][plot.mfrm_dif].
#'
#' @return Object of class `mfrm_dif_interaction` with:
#' - `table`: tibble with per-cell statistics and flags.
#' - `cases`: requested-row matching and exclusion audit.
#' - `summary`: tibble summarizing flagged and sparse cell counts.
#' - `gpcm_boundary`: for `GPCM` fits, a capability-boundary table.
#' - `config`: list of analysis parameters.
#'
#' @seealso [analyze_dff()], [analyze_dif()], [plot(x, type = "heatmap")][plot.mfrm_dif], [dif_report()],
#'   [estimate_bias()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("dif_interaction_table", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 300)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' int <- dif_interaction_table(fit, diag, facet = "Rater",
#'                              group = "Group", data = toy, min_obs = 2)
#' int$summary
#' head(int$table[, c("Level", "GroupValue", "ObsExpAvg", "flag_bias")])
#' }
#' @export
#' @order 1
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

  matched <- .dff_match_observations(fit, group_info$input_data, group, normalized = is.null(data))
  obs_work <- matched$data

  # Aggregate by facet level x group
  int_table <- obs_work |>
    group_by(.data[[facet]], .data$.group_var) |>
    summarise(
      N = n(),
      NPersons = n_distinct(.data$Person),
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
      t = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      FormalInferenceEligible = FALSE,
      ReportingUse = "descriptive",
      Interpretation = "Observed-minus-expected scores; scaled residuals are descriptive and are not t statistics."
    )
  names(int_table)[names(int_table) == facet] <- "Level"
  names(int_table)[names(int_table) == ".group_var"] <- "GroupValue"

  int_table$p_adjusted <- NA_real_
  int_table$flag_t <- NA
  int_table$flag_bias <- ifelse(int_table$sparse, NA, abs(int_table$ObsExpAvg) > abs_bias_warn)

  # Summary
  n_total <- nrow(int_table)
  n_sparse <- sum(int_table$sparse, na.rm = TRUE)
  n_flag_bias <- sum(int_table$flag_bias == TRUE, na.rm = TRUE)
  int_summary <- tibble(
    Metric = c("Total cells", "Sparse cells (N < min_obs)",
               "Above absolute residual mean threshold"),
    Count = c(n_total, n_sparse, n_flag_bias)
  )

  out <- list(
    table = int_table,
    cases = matched$cases,
    summary = int_summary,
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "dif_interaction_table()",
      area = "Differential facet functioning screening under GPCM"
    ),
    config = list(facet = facet, group = group, min_obs = min_obs,
                  method = "residual", residual_output_version = 3L,
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
  validate_dff_residual_output(object)
  out <- list(
    table = object$table,
    cases = object$cases,
    summary = object$summary,
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
    config = object$config
  )
  class(out) <- "summary.mfrm_dif_interaction"
  out
}

#' @export
print.summary.mfrm_dif_interaction <- function(x, ...) {
  validate_dff_residual_output(x)
  label <- x$config$functioning_label %||% "DFF"
  cat("--- Residuals by group and facet ---\n")
  cat("Facet:", x$config$facet, " | Group:", x$config$group, "\n")
  cat("Groups:", paste(x$config$group_levels, collapse = ", "), "\n")
  cat("Min obs:", x$config$min_obs, " | Absolute residual mean threshold:",
      x$config$abs_bias_warn, "score units\n")
  cat("Residual differences do not isolate differential functioning; no test is provided.\n\n")

  cat("Cell Summary:\n")
  print(as.data.frame(x$summary), row.names = FALSE)
  cat("\n")

  if (nrow(x$table) > 0) {
    show_cols <- intersect(
      c("Level", "GroupValue", "N", "ObsExpAvg", "StdResidual",
        "sparse", "flag_bias"),
      names(x$table)
    )
    print(as.data.frame(x$table |> select(all_of(show_cols))),
          row.names = FALSE, digits = 3)
  }
  .print_dff_gpcm_boundary(x$gpcm_boundary)
  invisible(x)
}

#' @rawNamespace S3method(print,mfrm_dif_interaction,print.mfrm_comparison)
NULL

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
#'   average (default), `"t"` for the descriptive scaled residual (not a t-test),
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
#' 2. Plot with `plot(..., type = "heatmap")`.
#' 3. Identify extreme cells or contrasts for follow-up.
#'
#' @return Invisibly, an `mfrm_plot_data` object whose `data` slot bundles
#'   the row x column metric matrix (`$matrix`), the source long table
#'   (`$pairs`), and the metric label. Earlier 0.1.x releases returned the
#'   bare matrix; consume `$data$matrix` to keep code forward-compatible.
#'
#' @seealso [dif_interaction_table()], [analyze_dff()], [analyze_dif()], [dif_report()]
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
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

  validate_dff_residual_output(x)

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
    max_abs <- scale_limit %||% max(c(0, abs(mat[is.finite(mat)])))
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
      t       = "Scaled residual",
      contrast = "Group contrast"
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
      main = "",
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
    t       = "Scaled residual",
    contrast = "Group contrast"
  )
  out <- new_mfrm_plot_data(
    "dif_heatmap",
    list(
      matrix = mat,
      flag_matrix = flag_matrix,
      pairs = as.data.frame(tbl, stringsAsFactors = FALSE),
      cases = x$cases,
      analysis_config = x$config,
      comparison_review = x$comparison_review,
      metric = metric,
      value_column = value_col,
      title = paste("Group comparisons:", metric_label),
      subtitle = if (identical(method, "refit")) {
        "Linked differences condition on the estimated anchors; formal DFF inference remains unavailable."
      } else "Residual differences do not isolate differential functioning.",
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
  guide <- out$data$interpretation_guide
  out$data$notes <- rbind(.mfrm_plot_annotation_notes(out$data),
    data.frame(Type = guide$Item, Text = paste(guide$Meaning, guide$ReportingNote)))
  if (!is.null(x$comparison_review)) out$data$notes <- rbind(out$data$notes,
    data.frame(Type = "Comparison review", Text = "A separate prospective interaction review is retained in comparison_review; its states do not classify the displayed residual/refit effects."))
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
#' trait (theta) for a fitted ordered-category `RSM`, `PCM`, or
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
#' with that cell's realized `step_facet` level. For `GPCM`, the
#' same design-weighted score variance is scaled by the squared discrimination
#' attached to the realized `slope_facet` level, which is the
#' \eqn{a_j^2 \cdot \mathrm{Var}(T \mid \theta)} item-information identity that
#' Muraki (1993, Equation 10) derives by applying Samejima's (1974)
#' polytomous information formula to the GPCM kernel of Muraki (1992).
#'
#' @section What `tif` and `iif` mean here:
#' In `mfrmr`, this helper supports ordered-category `RSM`, `PCM`, and the
#' current `GPCM` fit. The total curve (`$tif`) is the sum of
#' design-weighted cell contributions across all non-person facet levels in the
#' fitted model. The facet-level contribution curves (`$iif`) keep those
#' weighted contributions separated, so you can see which observed rater
#' levels, criteria, or other facet levels are driving precision at different
#' parts of the scale. For `PCM`, step-facet-specific thresholds are respected
#' when each observed design cell is evaluated. For `GPCM`, those
#' same cell-level variances are additionally scaled by the squared
#' discrimination associated with the realized `slope_facet` level.
#'
#' @section What this quantity does not justify:
#' - It is not a textbook many-facet test-information function for an abstract
#'   fixed item set.
#' - It should not be used as if it were design-free evidence about a form's
#'   precision independent of the realized observation pattern.
#' - It does not currently extend beyond the ordered-category `RSM` / `PCM` /
#'   `GPCM` family implemented by [fit_mfrm()].
#'
#' @section When to use this:
#' Use `compute_information()` when you want a design-weighted precision screen
#' for an `RSM`, `PCM`, or `GPCM` fit along the latent
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
#' formulation and Masters' `PCM`. The `GPCM` information identity
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
#'   under `GPCM`.)
#' - Samejima, F. (1974). *Normal ogive model on the continuous
#'   response level in the multidimensional latent space*.
#'   Psychometrika, 39, 111-121. (Source for the general polytomous
#'   information formula that Muraki 1993 specializes to the GPCM.)
#'
#' @section Interpreting output:
#' Category variance is computed from squared deviations about its mean to
#' retain small positive information near an endpoint category. This avoids
#' cancellation in raw second moments; it does not eliminate floating-point
#' underflow or account for calibration uncertainty.
#'
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
#' 3. Plot with `plot(info, type = "tif")`.
#' 4. If needed, inspect facet contributions with
#'    `plot(info, type = "iif", facet = "Rater")`.
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
#' @seealso [fit_mfrm()], [plot(x)][plot.mfrm_information]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("compute_information", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", model = "RSM",
#'                  quad_points = 7, maxit = 30)
#' info <- compute_information(fit)
#' head(info$tif)
#' info$tif$Theta[which.max(info$tif$Information)]
#' }
#' @export
#' @order 1
compute_information <- function(fit,
                                theta_range = c(-6, 6),
                                theta_points = 201L) {
  stop_if_product_slopes(fit, "compute_information()")
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an `mfrm_fit` object.", call. = FALSE)
  }
  model <- as.character(fit$config$model %||% NA_character_)

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)

  # Extract model parameters and realized observation design.
  step_structure <- information_build_step_structure(fit, model)
  categories <- step_structure$categories

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
    variance <- mfrm_category_variance(probs, categories)
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
#'   supported `RSM`, `PCM`, or `GPCM` fit.
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
#' - For `GPCM`, those contributions include the squared
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
#' 2. Plot with `plot(info)` for the total precision curve.
#' 3. Use `plot(info, type = "iif", facet = "Rater")` for
#'    facet-level contributions.
#' 4. Use `draw = FALSE` when you want reusable plot data for custom graphics
#'    or reporting helpers.
#'
#' @return Invisibly, an `mfrm_plot_data` object.
#'
#' @seealso [compute_information()], [fit_mfrm()]
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
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
           main = "", ...)
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
           main = "", ...)
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
           main = "", ...)
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
#' @param diagnostics Optional output from [diagnose_mfrm()]. When supplied,
#'   the map uses matching standard errors and precision metadata while keeping
#'   coordinates from `fit`.
#' @param bins Integer number of bins for the person histogram. Default `20`.
#' @param show_thresholds Logical; if `TRUE`, display threshold/step
#'   positions on the map. Default `TRUE`.
#' @param top_n Maximum number of facet/step locations retained by the native
#'   renderer for a compact display. Step transitions are always retained and
#'   omitted facet locations are reported in `retention`; use `Inf` for the
#'   complete final map. Every retained native location is labelled using
#'   collision-aware displaced text and leader lines; step transitions are
#'   shown as one vertical ladder with fitted logits in their labels. The
#'   FACETS-style payload retains and labels every fitted location, grouping
#'   coincident labels when needed.
#' @param show_ci Logical or `NULL`. `NULL` (the default) draws available
#'   uncertainty intervals for the native renderer and omits them from the
#'   FACETS-style renderer. Explicit `TRUE` with `renderer = "facets"` creates
#'   a hybrid FACETS-style ruler with mfrmr uncertainty intervals.
#' @param ci_level Confidence level used when `show_ci = TRUE`.
#' @param draw If `TRUE` (default), draw the plot. If `FALSE`, return
#'   plot data invisibly.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
#' @param palette Optional named color overrides passed to the shared Wright-map
#'   drawer.
#' @param label_angle Rotation angle for group labels on the facet panel.
#' @param renderer Canonical Wright-map renderer selector: `"native"`
#'   (default) or `"facets"` for the FACETS Table 6-style visual layout.
#' @param wright_style Wright-map renderer: `"native"` preserves the histogram,
#'   point, range, and facet-SE display; `"facets_style"` adds a FACETS Table
#'   6-style text ruler. The latter is a visual layout, not a claim of numerical
#'   equivalence with FACETS, and is equivalent to `renderer = "facets"`.
#' @param category_labels Optional score-rubric labels for
#'   `wright_style = "facets_style"`. Supply a named character vector keyed by
#'   every retained original score, an unnamed vector with one label per
#'   retained category, or a data frame with `Score` and `Label` columns.
#' @param rows_per_logit Number of rows per logit on the FACETS-style ruler.
#' @param wright_range Optional finite increasing length-2 logit range.
#' @param extreme_placement Place extreme-score persons at ruler `"ends"` or at
#'   their fitted `"estimate"` in the FACETS-style renderer.
#' @param persons_per_star Number of persons represented by one `*`; `NULL`
#'   selects a compact value automatically.
#' @inheritParams plot.mfrm_fit
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
#' If the fit records boundary-separated facet levels and `wright_range` is
#' `NULL`, both renderers derive the display range from supported locations and
#' place separated levels at ruler ends. The returned tables retain exact
#' `OriginalEstimate` and `CI_Lower` / `CI_Upper` values alongside display and
#' clipping metadata; endpoint triangles and footers disclose the adjustment.
#'
#' With `renderer = "facets"` (or `wright_style = "facets_style"`), the
#' draw-free result additionally
#' contains tidy ruler rows, person star frequencies, signed facet headers,
#' all facet levels, step lines, original-score transitions, mean half-score
#' boundaries, category labels, and display settings under `facets_style`.
#' The payload keeps both the nearest line-printer `RulerValue` and the exact
#' step/midpoint `DrawValue`; the current renderer draws threshold lines at the
#' exact value and prints fitted logits in step labels.
#' These tables support custom ggplot2/plotly rendering without parsing the
#' base plot.
#' Use `show_ci = FALSE` for the closest FACETS-style presentation. A
#' FACETS-style ruler drawn with `show_ci = TRUE` is intentionally labelled as
#' a hybrid because its uncertainty intervals are supplied by `mfrmr`.
#'
#' The returned list separates plotting availability from interpretation. It
#' includes `fit_readiness`, `interpretation_status`, and
#' `interpretation_note`. If numerical, data, connectivity, or stability
#' review is unresolved, the map is returned for diagnosis but the call warns
#' and prefixes the returned subtitle and drawn title with `REVIEW ONLY`.
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
#' 2. Plot with `plot(fit, type = "wright", view = "unified")`.
#' 3. Compare person distribution with facet level locations.
#' 4. Use `show_thresholds = TRUE` when you want the category structure in the
#'    same view.
#'
#' @section When to use this instead of plot_information:
#' Use `plot()` when your main question is targeting or coverage
#' on the shared logit scale. Use [plot(x)][plot.mfrm_information] when your main question
#' is measurement precision across theta.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @return Invisibly, a list with `persons`, `facets`, `thresholds`, and the
#'   underlying Wright-map tables used for the plot. Native output includes a
#'   `retention` table and `retention_note` documenting compact-display
#'   omissions. All renderers include fit-readiness and interpretation-status
#'   metadata.
#'
#' @seealso [fit_mfrm()], [plot.mfrm_fit()], [mfrmr_visual_diagnostics]
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept Wright maps
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_wright_unified", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(toy_small, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 30)
#' map_data <- plot(fit, type = "wright", view = "unified", draw = FALSE)
#' names(map_data)
#' facets_map <- plot(fit, type = "wright", view = "unified", renderer = "facets", category_labels = c(
#'     `1` = "Beginning", `2` = "Developing", `3` = "Secure", `4` = "Advanced"
#'   ), draw = FALSE)
#' facets_map$facets_style$settings
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_wright_unified <- function(fit,
                                diagnostics = NULL,
                                bins = 20L,
                                show_thresholds = TRUE,
                                top_n = 30L,
                                show_ci = NULL,
                                ci_level = 0.95,
                                draw = TRUE,
                                preset = c("standard", "publication", "compact", "monochrome"),
                                palette = NULL,
                                label_angle = 45,
                                renderer = NULL,
                                wright_style = c("native", "facets_style"),
                                category_labels = NULL,
                                rows_per_logit = 2L,
                                wright_range = NULL,
                                extreme_placement = c("ends", "estimate"),
                                persons_per_star = NULL,
                                highlight = NULL,
                                ...) {
  stop_if_product_slopes(fit, "plot(x, type = \"wright\", view = \"unified\")")
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an `mfrm_fit` object.", call. = FALSE)
  }
  show_ci_auto <- missing(show_ci) || is.null(show_ci)
  top_n <- suppressWarnings(as.numeric(top_n[1]))
  if (length(top_n) != 1L || is.na(top_n) || top_n <= 0) {
    stop("`top_n` must be a positive number or `Inf`.", call. = FALSE)
  }
  if (is.finite(top_n)) {
    top_n <- as.integer(floor(top_n))
    if (top_n < 1L) stop("`top_n` must be at least 1 or `Inf`.", call. = FALSE)
  }
  bins <- max(5L, as.integer(bins))
  style <- resolve_plot_preset(preset)
  wright_style <- match_wright_style(wright_style, renderer = renderer)
  if (isTRUE(show_ci_auto)) {
    show_ci <- identical(wright_style, "native")
  } else if (!is.logical(show_ci) || length(show_ci) != 1L || is.na(show_ci)) {
    stop("`show_ci` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  se_tbl_ci <- if (isTRUE(show_ci) || !is.null(diagnostics)) {
    compute_se_for_plot(fit, ci_level = ci_level, diagnostics = diagnostics)
  } else {
    NULL
  }
  plot_core <- build_wright_map_data(
    fit,
    top_n = top_n,
    highlight = highlight,
    se_tbl = se_tbl_ci,
    include_steps = isTRUE(show_thresholds),
    wright_style = wright_style,
    category_labels = category_labels,
    rows_per_logit = rows_per_logit,
    wright_range = wright_range,
    extreme_placement = extreme_placement,
    persons_per_star = persons_per_star
  )
  person_hist_values <- plot_core$person$Estimate[
    is.finite(plot_core$person$Estimate)
  ]
  plot_core$person_hist <- if (length(person_hist_values) > 0L) {
    graphics::hist(person_hist_values, breaks = bins, plot = FALSE)
  } else {
    structure(
      list(
        breaks = c(-0.5, 0.5), counts = 0L, density = 0,
        intensities = 0, mids = 0, xname = "finite Person estimates",
        equidist = TRUE
      ),
      class = "histogram"
    )
  }
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
  plot_data$show_ci <- isTRUE(show_ci)
  plot_data$uncertainty_display <- if (identical(wright_style, "facets_style")) {
    if (isTRUE(show_ci)) "hybrid_mfrmr_ci" else "facets_style_no_ci"
  } else if (isTRUE(show_ci)) {
    "native_mfrmr_ci"
  } else {
    "native_no_ci"
  }
  plot_data$legend <- .wright_map_legend(
    plot_core,
    style = style,
    show_ci = show_ci,
    ci_level = ci_level
  )
  plot_data$title <- if (identical(plot_data$uncertainty_display, "hybrid_mfrmr_ci")) {
    "Unified Wright Map (FACETS-style + mfrmr CI hybrid)"
  } else {
    "Unified Wright Map"
  }
  plot_data$subtitle <- if (identical(wright_style, "facets_style")) {
    if (isTRUE(show_ci)) {
      "Hybrid display: FACETS-style ruler with mfrmr uncertainty intervals"
    } else {
      "FACETS Table 6-style visual layout (mfrmr estimates)"
    }
  } else {
    "Shared logit scale for persons, facets, and thresholds"
  }
  plot_readiness <- .mfrm_fit_plot_readiness(fit)
  plot_data$fit_readiness <- plot_readiness$table
  plot_data$interpretation_status <- plot_readiness$status
  plot_data$interpretation_note <- plot_readiness$detail
  if (!isTRUE(plot_readiness$ready)) {
    plot_data$subtitle <- paste0("REVIEW ONLY - ", plot_data$subtitle)
    warning(plot_readiness$detail, call. = FALSE)
  }
  if (!draw) return(invisible(plot_data))

  apply_plot_preset(style)
  draw_wright_map(
    plot_data,
    title = plot_data$title,
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

.validate_linking_threshold <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0) {
    stop("`", arg, "` must be one non-negative number (Inf disables the threshold).",
      call. = FALSE)
  }
  as.numeric(x)
}

.linking_difference_se <- function(a, b) {
  scale <- pmax(a, b)
  valid <- is.finite(a) & is.finite(b) & a >= 0 & b >= 0
  out <- rep(NA_real_, length(scale))
  out[valid & scale == 0] <- 0
  use <- valid & scale > 0
  out[use] <- scale[use] * sqrt((a[use] / scale[use])^2 + (b[use] / scale[use])^2)
  out[!is.finite(out)] <- NA_real_
  out
}

.linking_cumulative_offsets <- function(offsets) {
  out <- cumsum(c(0, offsets))
  out[!is.finite(out)] <- NA_real_
  out
}

compute_equating_offset <- function(diffs, se_from = NULL, se_to = NULL,
                                    drift_threshold = NULL) {
  diffs <- as.numeric(diffs)
  if (!is.null(drift_threshold)) {
    drift_threshold <- .validate_linking_threshold(drift_threshold, "drift_threshold")
  }
  se_from <- if (is.null(se_from)) rep(NA_real_, length(diffs)) else as.numeric(se_from)
  se_to <- if (is.null(se_to)) rep(NA_real_, length(diffs)) else as.numeric(se_to)
  if (length(se_from) != length(diffs) || length(se_to) != length(diffs)) {
    stop("Differences and paired SEs must have the same length.", call. = FALSE)
  }
  ok <- is.finite(diffs)
  if (!any(ok)) {
    return(list(
      offset_prelim = NA_real_,
      offset = NA_real_,
      residual = rep(NA_real_, length(diffs)),
      retained = rep(FALSE, length(diffs)),
      n_retained = 0L,
      weighting = "none",
      contributing = rep(FALSE, length(diffs)),
      screening_fallback = FALSE,
      status = if (length(diffs)) "no_finite_differences" else "no_common",
      reason = "No finite common-element differences are available."
    ))
  }

  weight_ok <- ok & is.finite(se_from) & is.finite(se_to) & se_from > 0 & se_to > 0
  log_weights <- rep(NA_real_, length(diffs))
  high <- pmax(se_from[weight_ok], se_to[weight_ok])
  low <- pmin(se_from[weight_ok], se_to[weight_ok])
  log_weights[weight_ok] <- -2 * log(high) - log1p((low / high)^2)
  average <- function(retained) {
    weighted <- retained & weight_ok
    if (!any(weighted)) return(mean(diffs[retained]))
    d <- diffs[weighted]
    if (all(d == d[1L])) return(d[1L])
    lw <- log_weights[weighted]
    lw <- lw - max(lw)
    weights <- exp(lw)
    if (all(weights > 0)) {
      scale <- max(abs(d))
      center <- stats::weighted.mean(d / scale, weights)
      value <- center * scale
      if (value == 0 && center != 0) return(NA_real_)
      return(max(min(d), min(max(d), value)))
    }
    lw <- lw - log(sum(weights))
    nonzero <- d != 0
    if (!any(nonzero)) return(0)
    # Form weighted terms before exponentiating: a tiny weight can multiply
    # a large difference to give a representable contribution.
    terms <- log(abs(d[nonzero])) + lw[nonzero]
    shift <- max(terms)
    total <- sum(sign(d[nonzero]) * exp(terms - shift))
    if (total == 0) return(0)
    value <- sign(total) * exp(shift + log(abs(total)))
    if (value == 0) return(NA_real_)
    # The exact positive-weight mean is bounded by the input range.
    max(min(d), min(max(d), value))
  }
  offset_prelim <- average(ok)
  if (!is.finite(offset_prelim)) {
    return(list(offset_prelim = offset_prelim, offset = NA_real_,
      residual = rep(NA_real_, length(diffs)), retained = rep(FALSE, length(diffs)),
      n_retained = 0L, weighting = "none", contributing = rep(FALSE, length(diffs)),
      screening_fallback = FALSE, status = "numerical_failure",
      reason = "The preliminary offset is not numerically representable."))
  }

  residual_prelim <- diffs - offset_prelim
  retained <- ok
  if (!is.null(drift_threshold) && is.finite(drift_threshold)) {
    retained <- retained & is.finite(residual_prelim) & abs(residual_prelim) <= drift_threshold
  }
  screening_fallback <- !any(retained)
  if (screening_fallback) {
    retained <- ok
  }

  weight_retained <- retained & weight_ok
  offset <- average(retained)
  residual <- diffs - offset
  residual[!is.finite(residual)] <- NA_real_

  list(
    offset_prelim = offset_prelim,
    offset = offset,
    residual = residual,
    retained = retained,
    n_retained = sum(retained, na.rm = TRUE),
    weighting = if (any(weight_retained)) "inverse_variance" else "unweighted",
    contributing = if (any(weight_retained)) {
      weight_retained
    } else retained,
    screening_fallback = screening_fallback,
    status = if (is.finite(offset)) "computed" else "numerical_failure",
    reason = if (!is.finite(offset)) {
      "The final offset is not numerically representable."
    } else if (screening_fallback) {
      "No row passed preliminary screening; all finite differences were used."
    } else NA_character_
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
  new_est <- make_anchor_table(
    fit,
    include_person = FALSE,
    readiness_policy = "review"
  )

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
      SE_Diff = .linking_difference_se(.data$SE_Baseline, .data$SE_New),
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
#' This transfers baseline coordinates into a new fit through direct equality
#' constraints. Common-scale interpretation is conditional on a defensible
#' baseline and invariant cross-run element identity.
#' For `GPCM`, treat this as direct exploratory anchor/drift support
#' rather than as the package's formal linking-synthesis route.
#'
#' @param new_data Data frame in long format (one row per rating).
#' @param baseline_fit An inference-ready `mfrm_fit` object from a previous
#'   calibration under the current readiness contract.
#' @param person Character column name for person/examinee.
#' @param facets Character vector of facet column names.
#' @param score Character column name for the rating score.
#' @param anchor_facets Character vector of facets to anchor (default: all
#'   non-Person facets).
#' @param include_person If `TRUE`, also export Person values for a receiving
#'   JML fit. MML rejects Person anchors. Use only when constant ability for the
#'   same Persons is substantively intended; shared Persons need not be fixed.
#' @param weight Optional character column name for observation weights.
#' @param model Scale model override; defaults to baseline model.
#' @param method Estimation method override; defaults to baseline method.
#' @param anchor_policy How to handle anchor issues: `"warn"`, `"error"`,
#'   `"silent"`.
#' @param ... For `anchor_to_baseline()`, additional arguments passed to
#'   [fit_mfrm()]. Ignored by the S3 print and summary methods on this page.
#'
#' @details
#' This function automates the baseline-anchored calibration workflow:
#'
#' 1. Extracts anchor values from the baseline fit using [make_anchor_table()].
#' 2. Re-estimates the model on `new_data` with those anchors fixed via
#'    `fit_mfrm(..., anchors = anchor_table)`.
#' 3. Runs [diagnose_mfrm()] on the anchored fit.
#' 4. Compares the fixed non-Person anchor values in the new fit with their
#'    supplied values. This checks constraint consistency, not invariance.
#'
#' The helper refuses a baseline unless its current readiness record has
#' `InferenceReady = TRUE`. Also confirm compatible model, score coding,
#' orientation, and population conventions, and document that anchor labels
#' denote the same elements with invariant meaning. The function checks object
#' shape and receiving-data anchor compatibility, but does not establish those
#' substantive conditions. Fixing parameters can identify coordinates; it does
#' not create observed common-rating links.
#'
#' Only `model`, `method`, the centering convention and MML integration mode inherit baseline
#' settings automatically. Supply other settings explicitly through `...`,
#' including rubric bounds, category policy, facet orientation, PCM/GPCM facet
#' roles and any population model. RSM/PCM MML without a population model uses
#' a fixed standard-normal ability distribution; anchors do not automatically
#' estimate a new cohort's mean or variance.
#'
#' The table fixes locations only: thresholds and GPCM slopes are re-estimated.
#' Use [extract_mfrm_calibration()] and [score_mfrm_calibration()] when the question is scoring
#' against a complete fixed calibration. Partial anchors retain the centering
#' restrictions described in [fit_mfrm()]; this is not an unrestricted
#' fixed-common-parameter linking implementation for arbitrary new levels.
#' Elements present in the
#' anchor table but absent from the new data are handled according to
#' `anchor_policy`: `"warn"` (default) emits a message, `"error"` stops
#' execution, and `"silent"` ignores silently.
#'
#' The returned `drift` table contains only non-Person anchor rows that match
#' the new data. Unanchored common levels and Person levels are absent. The
#' compared values are fixed by construction, so zero differences do not show
#' that raters or tasks remained stable. For separate cross-wave screening,
#' fit the waves independently and use [detect_anchor_drift()], retaining its
#' offset and uncertainty limitations.
#'
#' Differences compare each retained fixed value to the supplied anchor:
#' \deqn{\Delta_e = \hat{\delta}_{e,\text{new}} - \hat{\delta}_{e,\text{base}}}
#' An element is **flagged** when \eqn{|\Delta_e| > 0.5} logits or
#' \eqn{|\Delta_e / SE_{\Delta_e}| > 2.0}, where
#' \eqn{SE_{\Delta_e} = \sqrt{SE_{\mathrm{base}}^2 + SE_{\mathrm{new}}^2}}.
#' `SE_Diff` is a plug-in independence calculation. It does not propagate
#' baseline-anchor uncertainty or cross-fit covariance, and anchored estimates
#' in the new fit are constrained rather than independently re-estimated.
#' Consequently, `Drift_SE_Ratio` and `Flag` are descriptive consistency
#' screens, not formal z tests or calibrated error-rate decisions.
#'
#' @section Which function should I use?:
#' - Use `anchor_to_baseline()` when you have one new dataset and want to place
#'   it on baseline coordinates within the documented constraint system.
#' - Use [detect_anchor_drift()] when you already have multiple fitted waves
#'   and want to compare their stability.
#' - Use [build_equating_chain()] when you need cumulative offsets across an
#'   ordered series of waves.
#'
#' @section Interpreting output:
#' - `$drift`: one row per retained non-Person anchor with columns `Facet`, `Level`,
#'   `Baseline`, `New`, `Drift`, `SE_Baseline`, `SE_New`, `SE_Diff`,
#'   `Drift_SE_Ratio`, and `Flag`.
#'   Read this as a constraint-consistency table. Zero differences are expected
#'   by construction and provide no independent evidence of stability.
#' - `$fit`: the full anchored `mfrm_fit` object, usable with
#'   [diagnose_mfrm()], [measurable_summary_table()], etc.
#' - `$diagnostics`: pre-computed diagnostics for the anchored calibration.
#' - `$baseline_anchors`: the candidate table fed to [fit_mfrm()]. Some rows may
#'   be excluded when `anchor_policy` permits unmatched IDs. Inspect
#'   `$fit$config$anchor_review$anchors` for the retained direct constraints.
#'
#' @section Typical workflow:
#' 1. Fit the baseline model: `fit1 <- fit_mfrm(...)`.
#' 2. Collect new data (e.g., a later administration).
#' 3. Call `res <- anchor_to_baseline(new_data, fit1, ...)`.
#' 4. Inspect `summary(res)` for numerical readiness and fixed-value consistency;
#'    neither establishes anchor invariance or a successful substantive link.
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
#' # Long-running workflow; run explicitly when needed.
#' # example("anchor_to_baseline", run.dontrun = TRUE)
#' \dontrun{
#' # Packaged synthetic teaching data; no data-generation code is needed.
#' ratings <- load_mfrmr_data("example_core")
#' people <- unique(ratings$Person)
#' # Two balanced 24-Person waves retain every Rater and Criterion.
#' d1 <- ratings[ratings$Person %in% people[1:24], , drop = FALSE]
#' d2 <- ratings[ratings$Person %in% people[25:48], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", rating_min = 1, rating_max = 4,
#'                  category_policy = "preserve")
#' res <- anchor_to_baseline(d2, fit1, "Person",
#'                           c("Rater", "Criterion"), "Score",
#'                           anchor_facets = "Criterion", anchor_policy = "error",
#'                           rating_min = 1, rating_max = 4,
#'                           category_policy = "preserve")
#' summary(res)
#' # These fixed-value differences are not a test of Criterion stability.
#' head(res$drift[, c("Facet", "Level", "Drift", "Flag")])
#' res$baseline_anchors[1:3, ]
#' }
#' @param mml_integration MML integration mode passed to [fit_mfrm()]. `NULL`
#'   inherits the baseline fit's mode for MML; JML uses `"fixed"`.
#' @param centering_scope `NULL` inherits the baseline's convention. Older fits
#'   without this setting retain `"remainder"`. Supply `"facet"` or `"remainder"`
#'   explicitly only when that constraint system is intended; see [fit_mfrm()].
anchor_to_baseline <- function(new_data, baseline_fit,
                               person, facets, score,
                               anchor_facets = NULL,
                               include_person = FALSE,
                               weight = NULL,
                               model = NULL, method = NULL,
                               anchor_policy = "warn",
                               mml_integration = NULL,
                               centering_scope = NULL,
                               ...) {
  # Validate baseline_fit
  stopifnot(inherits(baseline_fit, "mfrm_fit"))

  # Inherit model/method from baseline if not specified
  if (is.null(model))  model  <- baseline_fit$config$model
  if (is.null(method)) method <- baseline_fit$config$method
  if (is.null(centering_scope)) {
    centering_scope <- baseline_fit$config$centering_scope %||% "remainder"
  }
  if (is.null(mml_integration)) {
    mml_integration <- if (identical(toupper(method), "MML")) {
      baseline_fit$config$estimation_control$mml_integration %||% "fixed"
    } else "fixed"
  }

  # Extract anchor table from baseline
  anchor_tbl <- make_anchor_table(baseline_fit, facets = anchor_facets,
                                   include_person = include_person)

  if (nrow(anchor_tbl) == 0) {
    stop("No anchors could be extracted from the baseline fit.", call. = FALSE)
  }

  # Fit new data with anchors
  new_fit <- fit_mfrm(new_data, person = person, facets = facets, score = score,
                      weight = weight, model = model, method = method,
                      anchors = anchor_tbl, anchor_policy = anchor_policy,
                      mml_integration = mml_integration,
                      centering_scope = centering_scope, ...)

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
#' @name print.mfrm_anchored_fit
#' @aliases print.mfrm_anchored_fit
#' @usage \method{print}{mfrm_anchored_fit}(x, ...)
#' @rawNamespace S3method(print,mfrm_anchored_fit,print.mfrm_comparison)
NULL

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
    converged = mfrm_inference_ready(object$fit),
    optimizer_code_zero = mfrm_convergence_state(object$fit)$code_converged,
    convergence_status = mfrm_convergence_state(object$fit)$status
  )
  class(out) <- "summary.mfrm_anchored_fit"
  out
}

#' @rdname anchor_to_baseline
#' @export
print.summary.mfrm_anchored_fit <- function(x, ...) {
  cat("--- Anchored Fit Summary ---\n")
  cat("Inference-ready:", x$converged, "\n")
  cat("Optimizer returned code 0:", x$optimizer_code_zero,
      "| Status:", x$convergence_status, "\n")
  cat("Anchor candidates:", x$n_anchored, "| Matched non-Person anchors:", x$n_common,
      "| Flagged:", x$n_flagged, "\n\n")
  cat("Fixed-anchor differences check constraints, not stability across administrations.\n")
  if (nrow(x$drift_summary) > 0) {
    cat("Fixed-anchor differences by facet:\n")
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
#' elements whose difficulty/severity has shifted beyond supplied screening thresholds.
#' Useful for monitoring rater drift over time or checking the stability of
#' item banks.
#'
#' @param fits Named list of `mfrm_fit` objects (e.g.,
#'   `list(Year1 = fit1, Year2 = fit2)`).
#' @param facets Character vector of facets to compare (default: all
#'   non-Person facets).
#' @param drift_threshold Absolute screening threshold for flagging (logits,
#'   default 0.5). A non-negative scalar; `Inf` disables magnitude screening.
#' @param flag_se_ratio Drift/SE ratio screening threshold for flagging
#'   (default 2.0). A non-negative scalar; `Inf` disables SE-ratio flagging.
#' @param reference Index or name of the reference fit (default: first).
#' @param include_person Include person estimates in comparison. Use only when
#'   the person identifiers refer to the same persons across fitted waves.
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
#' This count is a package screening convention, not a universal adequacy
#' threshold.
#'
#' One `LinkOffset` is pooled across all selected common elements and facets.
#' If at least one pair has finite positive SEs, only pairs with such SEs
#' contribute to the inverse-variance weighted offset; rows without usable SEs
#' can still appear in retained/support counts. If preliminary residual
#' screening would retain no row, the current implementation falls back to all
#' finite differences. Use a single substantively coherent facet when a common
#' shift across facet blocks is not defensible, and inspect `drift_table` rather
#' than treating `LinkSupportAdequate` as sufficient evidence.
#' The fallback is recorded as `ScreeningFallback`, and `Contributing` /
#' `N_Contributing` identify the elements actually used in the final offset.
#' Retention and contribution are different: a retained row with a missing SE
#' does not enter a weighted offset when other rows have usable SEs.
#'
#' `OffsetStatus` and `OffsetReason` distinguish computed offsets, absent or
#' nonfinite common differences, and numerical failure. Missing residuals stay
#' missing. `SERatioStatus` explains unavailable SE ratios; `FlagBasis` identifies
#' whether the flag used both criteria, only residual magnitude, or neither.
#' In the magnitude-only case, `Flag = FALSE` is not a negative SE-ratio screen.
#' Summaries report unavailable counts, and withhold the complete flagged count
#' when any residual flag is unavailable.
#'
#' Matching `Facet` and `Level` labels are treated as common elements; the
#' helper cannot verify that they have invariant meaning across waves. The
#' input fits must therefore use compatible model, score, orientation, and
#' population conventions and must be independently reviewed for fit readiness.
#' Offset removal aligns coordinates conditionally on those assumptions; it is
#' not proof of scale equivalence.
#'
#' This is a relative drift screen. For example, if the selected judges'
#' estimates all differ from the reference by +0.40 logits, `LinkOffset` is
#' +0.40 and each residual `Drift` is zero. This does not establish panel
#' stability. The offset itself also cannot distinguish real panel change
#' from a change of coordinates. Common recorded performances can provide
#' a separate reference, but a complete recording-based annual analysis and
#' its formal difference intervals are not yet verified here. See the annual
#' assessment example in [mfrmr_linking_and_dff].
#'
#' An element is **flagged** when either condition is met:
#' \deqn{|\Delta_e| > \texttt{drift\_threshold}}
#' \deqn{|\Delta_e / SE_{\Delta_e}| > \texttt{flag\_se\_ratio}}
#' The two criteria are joined by `OR`: a large absolute residual is flagged
#' even when imprecise, while a smaller residual can be flagged by its SE ratio.
#' The ratio uses only the two element SEs and omits uncertainty in the
#' estimated `LinkOffset` and any cross-fit covariance. It is therefore a
#' descriptive screen, not a formal z test with a calibrated error rate.
#' Residuals share an estimated offset, and contributing elements also enter
#' that offset. Offset covariance can either reduce or increase residual
#' uncertainty relative to the displayed denominator. The ratio is not
#' guaranteed to be conservative.
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
#' 4. Review `summary(drift)` and `plot(drift)`.
#' 5. Flagged elements may need to be removed from anchor sets or
#'    investigated for substantive causes (e.g., rater re-training).
#'
#' @return Object of class `mfrm_anchor_drift` with components:
#'   \describe{
#'     \item{drift_table}{Tibble of element-level drift statistics.}
#'     \item{summary}{Drift summary aggregated by facet and wave.}
#'     \item{common_elements}{Tibble of pairwise common-element counts.}
#'     \item{common_vs_reference}{Tibble of common-element counts
#'       between each wave and the reference wave (i.e., which labels remain
#'       available as candidate common elements).}
#'     \item{n_common_all_waves}{Integer count of elements that are
#'       common across every wave; used by `summary()` to gauge how
#'       robust the chain is to chained linking error.}
#'     \item{common_by_facet}{Tibble of retained common-element counts by facet.}
#'     \item{config}{List of analysis configuration.}
#'   }
#'
#' @seealso [anchor_to_baseline()], [build_equating_chain()],
#'   [make_anchor_table()], [plot(x)][plot.mfrm_anchor_drift], [mfrmr_linking_and_dff]
#' @export
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("detect_anchor_drift", run.dontrun = TRUE)
#' \dontrun{
#' # Deliberately linked teaching waves: common labels below represent the
#' # same rater and criterion identities by construction.
#' toy <- load_mfrmr_data("example_core")
#' people <- unique(toy$Person)
#' # Two balanced 24-Person waves retain every Rater and Criterion.
#' d1 <- toy[toy$Person %in% people[1:24], , drop = FALSE]
#' d2 <- toy[toy$Person %in% people[25:48], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", quad_points = 7, maxit = 30)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", quad_points = 7, maxit = 30)
#' drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
#' summary(drift)
#' head(drift$drift_table[, c("Facet", "Level", "Wave", "Drift", "Flag")])
#' drift$common_elements
#' }
#' @order 1
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
  drift_threshold <- .validate_linking_threshold(drift_threshold, "drift_threshold")
  flag_se_ratio <- .validate_linking_threshold(flag_se_ratio, "flag_se_ratio")
  valid_reference <- length(reference) == 1L && !is.na(reference) &&
    ((is.character(reference) && sum(names(fits) == reference) == 1L) ||
      (is.numeric(reference) && is.finite(reference) && reference == floor(reference) &&
        reference >= 1 && reference <= length(fits)))
  if (!valid_reference) stop("`reference` must identify one fit by a unique name or integer index.",
    call. = FALSE)
  ref_idx <- if (is.character(reference)) match(reference, names(fits)) else as.integer(reference)

  # Extract estimates from each fit
  est_list <- lapply(fits, function(f) {
    make_anchor_table(
      f,
      facets = facets,
      include_person = include_person,
      readiness_policy = "review"
    )
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
          OffsetStatus = offset_info$status,
          OffsetReason = offset_info$reason,
          ScreeningFallback = offset_info$screening_fallback,
          N_Contributing = sum(offset_info$contributing),
          Contributing = offset_info$contributing,
          Retained = offset_info$retained,
          Drift = offset_info$residual,
          SE = .linking_difference_se(.data$SE_Ref, .data$SE_Wave),
          Drift_SE_Ratio = ifelse(
            is.na(.data$SE) | .data$SE == 0,
            NA_real_,
            abs(.data$Drift) / .data$SE
          ),
          SERatioStatus = dplyr::case_when(
            !is.finite(.data$Drift) ~ "drift_unavailable",
            !is.finite(.data$SE) | .data$SE <= 0 ~ "paired_se_unavailable",
            !is.finite(.data$Drift_SE_Ratio) ~ "ratio_unrepresentable",
            TRUE ~ "available"
          ),
          Drift_SE_Ratio = ifelse(is.finite(.data$Drift_SE_Ratio), .data$Drift_SE_Ratio, NA_real_),
          FlagBasis = ifelse(!is.finite(.data$Drift), "unavailable",
            ifelse(is.finite(.data$Drift_SE_Ratio), "magnitude_and_se_ratio", "magnitude_only")),
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
                     "LinkSupportAdequate", "Flag", "OffsetStatus", "OffsetReason",
                     "ScreeningFallback", "N_Contributing", "Contributing", "Retained",
                     "SERatioStatus", "FlagBasis") |>
      dplyr::arrange(dplyr::desc(abs(.data$Drift)))
  } else {
    tibble::tibble(Facet = character(), Level = character(),
                   Reference = character(), Wave = character(),
                   Ref_Est = numeric(), Wave_Est = numeric(),
                   LinkOffset = numeric(), Drift = numeric(),
                   SE_Ref = numeric(), SE_Wave = numeric(), SE = numeric(),
                   Drift_SE_Ratio = numeric(), LinkSupportAdequate = logical(),
                   Flag = logical(), OffsetStatus = character(), OffsetReason = character(),
                   ScreeningFallback = logical(), N_Contributing = integer(),
                   Contributing = logical(), Retained = logical(),
                   SERatioStatus = character(), FlagBasis = character())
  }

  # Summary by facet
  drift_summary <- if (nrow(drift_table) > 0) {
    drift_table |>
      dplyr::group_by(.data$Facet, .data$Wave) |>
      dplyr::summarise(
        N = dplyr::n(), Mean_Drift = mean(abs(.data$Drift)),
        Max_Drift = max(abs(.data$Drift)),
        N_Flagged = if (anyNA(.data$Flag)) NA_integer_ else sum(.data$Flag),
        N_FlaggedAvailable = sum(.data$Flag, na.rm = TRUE),
        N_Unavailable = sum(is.na(.data$Flag)), .groups = "drop"
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
                  calculation_version = 2L,
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
#' @param x For `print.mfrm_anchor_drift()`: An `mfrm_anchor_drift` object.
#'
#' For `plot()`: A saved result from [detect_anchor_drift()].
#' @param ... For `print.mfrm_anchor_drift()`: Ignored.
#'
#' For `plot()`: `facet` optionally selects facets; `ci_level = 0.95` controls
#'   forest whiskers only; `preset` selects `"standard"`, `"publication"`,
#'   `"compact"`, or `"monochrome"`; `draw = FALSE` returns data without drawing.
#'   Additional base-graphics arguments depend on the selected view.
#' @name print.mfrm_anchor_drift
#' @aliases print.mfrm_anchor_drift
#' @usage \method{print}{mfrm_anchor_drift}(x, ...)
#' @rawNamespace S3method(print,mfrm_anchor_drift,print.mfrm_comparison)
NULL

#' @rdname detect_anchor_drift
#' @param object An `mfrm_anchor_drift` object (for `summary`).
#' @export
summary.mfrm_anchor_drift <- function(object, ...) {
  dt <- object$drift_table
  out <- list(
    n_comparisons = nrow(dt),
    n_flagged = if (anyNA(dt$Flag)) NA_integer_ else sum(dt$Flag),
    n_flagged_available = sum(dt$Flag, na.rm = TRUE),
    n_unavailable = sum(is.na(dt$Flag)),
    summary = object$summary, common_elements = object$common_elements,
    common_by_facet = object$common_by_facet,
    flagged = dt |> dplyr::filter(.data$Flag),
    config = object$config
  )
  class(out) <- "summary.mfrm_anchor_drift"
  out
}

linking_review_notes <- function(gpcm_detected = FALSE) {
  clean_summary_lines(c(
    "Review anchor overlap, fitted-wave drift and adjacent-link residuals together before comparing results across administrations.",
    "Drift or limited-overlap flags identify comparisons to investigate; they do not establish that the measurement scale has changed.",
    if (gpcm_detected) {
      "GPCM linking results are exploratory and do not establish that scores can be used interchangeably across administrations or that drift is absent."
    }
  ))
}

print_linking_uncertainty_note <- function() {
  print_wrapped_line(paste(
    "Linking flags are review screens, not tests of anchor invariance.",
    "Source-parameter, estimated-offset and cross-fit covariance are not fully propagated.",
    "A sufficient element count alone does not establish a common scale."
  ))
}

#' @rdname detect_anchor_drift
#' @export
print.summary.mfrm_anchor_drift <- function(x, ...) {
  cat("--- Anchor Drift Screen ---\n")
  cat("Reference:", x$config$reference, "\n")
  print_linking_uncertainty_note()
  cat("Comparisons:", x$n_comparisons, "| Flagged:", x$n_flagged, "\n\n")
  if (isTRUE(x$n_unavailable > 0)) cat("Unavailable:", x$n_unavailable,
    "| Flagged among available:", x$n_flagged_available, "\n\n")
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
#' scale shift. Cumulative offsets express all waves relative to the first
#' wave, conditional on the common-element assumptions. The procedure is a
#' practical screened linking aid, not a full general-purpose equating
#' framework or proof of common-scale comparability.
#'
#' @param fits Named list of `mfrm_fit` objects in chain order.
#' @param anchor_facets Character vector of facets to use as linking
#'   elements.
#' @param include_person Include person estimates in linking.
#' @param drift_threshold Non-negative scalar threshold for screening and
#'   flagging large residuals in links. `Inf` retains all finite differences
#'   and disables magnitude flagging.
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
#' The five-element rule is a package screening convention, not a universal
#' adequacy threshold. Matching labels are assumed to identify the same
#' invariant elements, and input fits are assumed to use compatible model,
#' score, orientation, and population conventions. The helper does not verify
#' those assumptions or source-fit readiness. In particular, adjacent label
#' overlap and a finite offset do not by themselves establish a common scale.
#'
#' Each adjacent `Offset` is one value pooled across all selected facets. When
#' any common rows have finite positive SEs, only those rows contribute to the
#' inverse-variance weighted offset; retained/support counts can still include
#' rows without usable SEs. If preliminary screening would remove every row,
#' the implementation falls back to all finite differences. Prefer one
#' substantively coherent `anchor_facets` block unless a common shift across
#' facets is justified.
#'
#' Cumulative offsets are computed by chaining link offsets from Wave 1
#' forward, placing all waves onto the metric of the first wave.
#' `Offset_SD` is retained-element residual spread, not the standard error of
#' `Offset`. No offset SE, confidence interval, cross-fit covariance, or
#' cumulative uncertainty propagation is currently returned, so uncertainty
#' can compound along the chain without appearing in `cumulative`.
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
#' - `OffsetStatus` / `OffsetReason` explain unavailable link offsets.
#'   `ScreeningFallback` records use of all finite differences when no row
#'   passed preliminary screening. `N_Contributing` and the element-level
#'   `Contributing` flag identify the rows actually used in the final average;
#'   retained rows without usable SEs may not contribute.
#' - `OffsetSDStatus` distinguishes available residual spread from fewer than
#'   two retained elements, unavailable offsets, or a numerical failure.
#'   A singleton has an undefined sample SD, not an estimated SD of zero.
#'   Unavailable `Max_Residual` and flagged counts stay missing rather than
#'   representing no drift. Summaries also give the unavailable count and
#'   the flagged count among available elements.
#' - `cumulative$Status` / `Reason` distinguish unavailable preceding links
#'   from an unrepresentable cumulative sum. Later offsets do not become zero
#'   when a link is unavailable.
#' - Read `links` before `cumulative`: weak adjacent links can make later
#'   cumulative offsets less trustworthy.
#'
#' @section Common-element graph:
#' `plot(chain, type = "graph")` displays waves as squares and common elements
#' as circles, including waves without retained connections. Each element node
#' belongs to one reviewed link: `[L1]` refers to link 1 in the returned links
#' table. The same element can have different screening results in different
#' comparisons, so it is repeated per link. Solid lines denote retained,
#' unflagged elements; dot-dash denotes retained elements with a drift flag or
#' an unavailable flag; dashed denotes excluded elements; dotted denotes
#' unknown retention. These distinctions also apply in monochrome.
#'
#' This is a graph of screened common-element comparisons, not an inventory
#' of fixed anchors or an all-pairs assessment. Retention does not guarantee
#' positive weight in the offset, invariance, adequate support or precision.
#' New chains store positional `FromID`/`ToID` in `links` and `LinkID` in
#' `element_detail`; display names may repeat or contain punctuation. Older
#' chains are matched by their complete link labels; ambiguous references
#' require rebuilding the chain.
#'
#' With `draw = FALSE`, the plot returns an `mfrm_plot_data` object. Its
#' `$data$data` contains `nodes`, `edges`, `links`, `elements`, `n_waves`,
#' `n_anchors` (unique facet/level identities), and `n_element_nodes`
#' (link-specific nodes). Edge `from`/`to` refer to node `NodeId`; `AnchorId`
#' identifies the same facet/level across links. `Retained`, `Flag`, `Status`
#' and `Reason` preserve screening information. `RetainedDegree` counts
#' retained incidences, not statistical information. Interpretation notes
#' remain in `$data$notes` and in `print()` even when `show_notes = FALSE`.
#' Full labels remain in the node table if the device is too narrow. Base
#' graphics are sufficient; no optional graph package is needed.
#'
#' @section Wave links and common-element removal:
#' `plot(chain, type = "links")` summarizes the recorded wave-to-wave
#' comparisons. Numbers count common elements with `Retained = TRUE`, including
#' flagged retained elements. A dotted comparison has no recorded retained
#' element; this is not proof that the original rating design is disconnected.
#' Missing retention is counted separately. Dot-dash connections have flagged
#' retained elements or inadequate/unknown recorded link support.
#'
#' `plot(chain, type = "anchor_removal")` removes each facet/level identity,
#' one at a time, from all recorded comparisons. Other retention decisions stay
#' fixed. The horizontal axis counts newly disconnected unordered wave pairs,
#' including indirect paths. Pairs already disconnected at baseline are not
#' counted. Losing a direct link need not disconnect its endpoints if another
#' path remains. Zero means no new graph disconnection, not negligible impact
#' on estimates or adequate linking support. No screening, offsets, estimates,
#' SEs or confidence intervals are recomputed by either view.
#'
#' Both views return `$data$data$nodes`, `links`, `elements` and `summary`.
#' The removal view additionally returns `removal` (one row per common-element
#' identity), `lost_links` (recorded comparisons losing their last retained
#' element), `lost_pairs` (newly disconnected wave pairs), and `components_after`
#' (baseline and post-removal membership for every wave). Use `AnchorId` to
#' join the removal table to these details; `FromID`/`ToID` and `WaveID` identify
#' waves independently of display labels. Component numbers are local to each
#' partition: compare membership relations, not numeric labels across Before,
#' After or different removal scenarios. Removal counts
#' include `RetainedOccurrences`, `UnknownOccurrences`, `LostLinks`,
#' `ComponentsAfter` and `NewlyDisconnectedPairs`.
#'
#' These calculations reuse the recorded common-element graph; they do not
#' search for new links or assess statistical influence. The full lost-pair
#' table can be large for chains with many waves and common elements. Compact
#' plots abbreviate long labels and warn when IDs or a larger device are needed;
#' full identities and interpretation notes remain in the returned tables.
#'
#' @section Conditional offset sensitivity:
#' `plot(chain, type = "offset_sensitivity")` removes each common `(Facet,
#' Level)` identity from every adjacent link, then reruns the existing offset
#' calculation with its recorded drift threshold. Preliminary offsets,
#' screening, final weighted/unweighted offsets and cumulative offsets are
#' recomputed. The source wave estimates and SEs remain fixed: this is a
#' conditional linking sensitivity calculation, not a refit of item, rater or
#' person parameters or a new uncertainty estimate.
#'
#' The view requires an ordered adjacent chain with the recorded method,
#' threshold, support guideline, source estimate/SE columns and offsets.
#' Recomputed baseline offsets must agree with recorded offsets within
#' `1e-8 * max(1, abs(recorded_offset))`; inconsistent or incomplete source
#' records require rebuilding the chain. Numerical overflow/underflow that
#' prevents a finite preliminary offset is reported as a numerical failure.
#'
#' The figure shows each deletion's maximum absolute cumulative offset change
#' over finite comparisons, excluding the fixed-zero first wave. Shapes
#' distinguish complete comparisons, partial comparisons and no finite
#' comparison. An unavailable comparison remains `NA`, not zero. If one link
#' becomes unavailable, later cumulative offsets also become unavailable.
#' Read the detailed tables when a figure marks an incomplete comparison.
#'
#' The `$data$data` payload contains `settings`, `baseline_links`, `baseline_cumulative`,
#' `baseline_elements`, `removal`, `links`, `cumulative`, `retention_changes`
#' and `common_by_facet`. Join `RemovedAnchorId` in scenario tables to
#' `removal$AnchorId`. Signed `Change` is recalculated minus baseline offset.
#' `Status` distinguishes `computed`, `no_common`, `no_finite_differences`
#' and `numerical_failure`. `N_Contributing` and element `Contributing` identify
#' positive finite weight contributors (or retained rows in an unweighted
#' offset). `retention_changes` records other elements whose retention,
#' contribution or residual flag changes; the deliberately removed element is
#' omitted. Facets whose last element is removed remain in support tables with
#' zero counts. No support failure is converted into a statement of readiness.
#'
#' `ScreeningFallback` identifies the existing rule that restores all finite
#' differences if preliminary screening would remove every one. SEs may be
#' missing, in which case the existing unweighted fallback is used. No offset
#' SE, confidence interval, cross-fit covariance, or uncertainty propagation
#' is supplied. A small conditional change is not evidence of anchor invariance
#' or negligible change under a full model refit. Notes remain available when
#' figure titles and footers are omitted.
#'
#' @section Typical workflow:
#' 1. Fit each administration wave separately: `fit_a <- fit_mfrm(...)`.
#' 2. Combine into an ordered named list:
#'    `fits <- list(Spring23 = fit_s, Fall23 = fit_f, Spring24 = fit_s2)`.
#' 3. Call `chain <- build_equating_chain(fits)`.
#' 4. Review `summary(chain)` for link quality.
#' 5. Visualize with `plot(chain, type = "chain")`.
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
#'   [make_anchor_table()], [plot(x)][plot.mfrm_anchor_drift]
#' @export
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("build_equating_chain", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' people <- unique(toy$Person)
#' d1 <- toy[toy$Person %in% people[1:24], , drop = FALSE]
#' d2 <- toy[toy$Person %in% people[25:48], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML")
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML")
#' # Check convergence before using the estimates for linking.
#' summary(fit1)
#' summary(fit2)
#' # The four shared levels per facet trigger a thin-support warning.
#' # Inspect that warning; matching labels alone do not establish a common scale.
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
  drift_threshold <- .validate_linking_threshold(drift_threshold, "drift_threshold")

  wave_names <- names(fits)
  n_waves <- length(fits)

  # Extract estimates
  est_list <- lapply(fits, function(f) {
    make_anchor_table(
      f,
      facets = anchor_facets,
      include_person = include_person,
      readiness_policy = "review"
    )
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
      offset_sd_status <- "no_common"
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
      retained_residual <- offset_info$residual[offset_info$retained]
      offset_sd <- if (length(retained_residual) >= 2L && all(is.finite(retained_residual))) {
        scale <- max(abs(retained_residual))
        if (scale == 0) 0 else stats::sd(retained_residual / scale) * scale
      } else {
        NA_real_
      }
      offset_sd_status <- if (is.finite(offset_sd)) "available" else if (!is.finite(offset)) {
        offset_info$status
      } else if (length(retained_residual) < 2L) "insufficient_elements" else "numerical_failure"
      if (!is.finite(offset_sd)) offset_sd <- NA_real_
      max_drift <- if (all(is.finite(offset_info$residual))) max(abs(offset_info$residual)) else NA_real_

      common <- common |>
        dplyr::mutate(
          Link = paste0(wave_names[i], " -> ", wave_names[i + 1]),
          LinkID = i,
          Diff = .data$Est_To - .data$Est_From,
          Offset_Prelim = offset_info$offset_prelim,
          Offset = offset,
          Residual = offset_info$residual,
          Retained = offset_info$retained,
          Contributing = offset_info$contributing,
          Flag = abs(.data$Residual) > drift_threshold
        )
      element_details <- c(element_details, list(common))
    }

    links <- c(links, list(tibble::tibble(
      Link = i,
      FromID = i, ToID = i + 1L,
      From = wave_names[i], To = wave_names[i + 1],
      N_Common = as.integer(n_common),
      N_Retained = if (n_common > 0) offset_info$n_retained else 0L,
      Min_Common_Per_Facet = if (n_common > 0 && nrow(support_tbl) > 0) min(support_tbl$N_Common) else 0L,
      Min_Retained_Per_Facet = if (n_common > 0 && nrow(support_tbl) > 0) min(support_tbl$N_Retained) else 0L,
      Offset_Prelim = if (n_common > 0) offset_info$offset_prelim else NA_real_,
      Offset = offset,
      Offset_SD = offset_sd,
      OffsetSDStatus = offset_sd_status,
      OffsetStatus = if (n_common > 0) offset_info$status else "no_common",
      OffsetReason = if (n_common > 0) offset_info$reason else "No common elements.",
      N_Contributing = if (n_common > 0) sum(offset_info$contributing) else 0L,
      ScreeningFallback = if (n_common > 0) offset_info$screening_fallback else FALSE,
      Max_Residual = max_drift,
      LinkSupportAdequate = if (n_common > 0) link_support_ok else FALSE,
      Offset_Method = if (n_common > 0) offset_info$weighting else "none"
    )))
  }

  links_tbl <- dplyr::bind_rows(links)

  # Cumulative offsets
  cum_offset <- .linking_cumulative_offsets(links_tbl$Offset)
  cumulative <- tibble::tibble(
    Wave = wave_names,
    Cumulative_Offset = cum_offset,
    Status = ifelse(is.finite(cum_offset), "computed", "unavailable"),
    Reason = ifelse(is.finite(cum_offset), NA_character_,
      ifelse(cumsum(c(0L, !is.finite(links_tbl$Offset))) > 0L,
        "An earlier or current link offset is unavailable.", "Cumulative sum is not numerically representable."))
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
                  calculation_version = 2L,
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
#' @name print.mfrm_equating_chain
#' @aliases print.mfrm_equating_chain
#' @usage \method{print}{mfrm_equating_chain}(x, ...)
#' @rawNamespace S3method(print,mfrm_equating_chain,print.mfrm_comparison)
NULL

#' @rdname build_equating_chain
#' @param y Unused (S3 plot signature requirement).
#' @param type One of `"graph"` (bipartite wave x link-specific common-element
#'   graph), `"common_anchors"` (default; bar
#'   chart of common-anchor counts per wave pair), `"chain"` (cumulative
#'   offsets), `"links"` (recorded wave comparisons), or `"anchor_removal"`
#'   (newly disconnected wave pairs after deleting each common element), or
#'   `"offset_sensitivity"` (conditional rescreened linking-offset changes).
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw the plot with base graphics.
#' @param show_title Logical; display the main title.
#' @param show_notes Logical; display graph/topology interpretation footers.
#'   Returned notes and ordinary R warnings are unaffected. Axes, legends and
#'   data labels remain visible; other chain views have no interpretation footer.
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_equating_chain <- function(x, y = NULL,
                                     type = c("common_anchors", "graph", "chain", "links", "anchor_removal", "offset_sensitivity"),
                                     preset = c("standard", "publication", "compact", "monochrome"),
                                     draw = TRUE, show_title = FALSE, show_notes = FALSE, ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(x, "mfrm_equating_chain")) {
    stop("`x` must be an mfrm_equating_chain object.", call. = FALSE)
  }
  for (arg in c("show_title", "show_notes")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
    }
  }
  type <- match.arg(type)
  style <- resolve_plot_preset(preset)
  if (identical(type, "chain")) {
    if (isTRUE(draw)) apply_plot_preset(style)
    return(.plot_equating_chain(x, draw = draw, style = style,
      show_title = show_title, show_notes = show_notes, ...))
  }
  if (identical(type, "offset_sensitivity")) {
    sensitivity <- .equating_offset_sensitivity(x)
    key <- data.frame(Label = c("All cumulative comparisons available", "Some comparisons unavailable", "No finite cumulative comparison"),
      Shape = c(19, 17, 4), Color = unname(.plot_series_colors(1:3, style$name)))
    out <- new_mfrm_plot_data("equating_chain_offset_sensitivity", list(
      data = sensitivity$data, notes = sensitivity$notes, styles = key,
      title = "Common-element offset sensitivity",
      subtitle = "Rescreened linking offsets; source estimates and SEs held fixed",
      display = list(show_title = show_title, show_notes = show_notes),
      legend = new_plot_legend(key$Label, rep("availability", 3), rep("shape", 3), key$Shape),
      preset = style$name
    ))
    if (isTRUE(draw)) .draw_equating_offset_sensitivity(out, style)
    return(invisible(out))
  }
  if (type %in% c("links", "anchor_removal")) {
    removal <- identical(type, "anchor_removal")
    topology <- .equating_link_topology(.equating_graph_data(x), removal = removal)
    if (removal) {
      key <- data.frame(Status = c("disconnected", "unchanged", "not_retained"),
        Label = c("New disconnection", "No new disconnection", "No retained occurrence"),
        Shape = c(19, 1, 4), stringsAsFactors = FALSE)
    } else {
      key <- data.frame(Status = c("retained", "retained_review", "no_retained"),
        Label = c("Retained connection", "Retained; review support", "No recorded retained connection"),
        Linetype = c("solid", "dotdash", "dotted"), stringsAsFactors = FALSE)
    }
    key$Color <- unname(.plot_series_colors(key$Status, style$name))
    out <- new_mfrm_plot_data(paste0("equating_chain_", type), list(
      data = topology$data, notes = topology$notes, styles = key,
      title = if (removal) "Common-element removal" else "Recorded wave links",
      subtitle = "Connectivity of recorded retained elements; no re-estimation",
      display = list(show_title = show_title, show_notes = show_notes),
      legend = new_plot_legend(key$Label, rep("topology", 3),
        rep(if (removal) "shape" else "linetype", 3),
        if (removal) key$Shape else key$Linetype),
      preset = style$name
    ))
    if (isTRUE(draw)) .draw_equating_topology(out, style, removal = removal)
    return(invisible(out))
  }
  if (identical(type, "graph")) {
    graph_data <- .equating_graph_data(x)
    status <- c("retained", "retained_review", "excluded", "candidate")
    colors <- .plot_series_colors(status, style$name)
    key <- data.frame(
      Status = status,
      Label = c("Retained", "Retained; review drift", "Excluded", "Retention unknown"),
      Linetype = c("solid", "dotdash", "dashed", "dotted"),
      Color = unname(colors), stringsAsFactors = FALSE
    )
    out <- new_mfrm_plot_data("equating_chain_graph", list(
      data = graph_data$data, notes = graph_data$notes,
      styles = key,
      title = "Screened common-element links",
      subtitle = sprintf("%d waves; %d unique elements; %d link-specific element nodes",
        graph_data$data$n_waves, graph_data$data$n_anchors, graph_data$data$n_element_nodes),
      display = list(show_title = show_title, show_notes = show_notes),
      legend = new_plot_legend(key$Label, rep("screening", 4), rep("linetype", 4), key$Linetype),
      preset = style$name
    ))
    if (isTRUE(draw)) .draw_equating_graph(out, style)
    return(invisible(out))
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
      display = list(show_title = show_title, show_notes = show_notes),
      subtitle = sprintf("%d pairwise comparison(s)", nrow(links)),
      preset = style$name
    )
  )
  out$data$notes <- .mfrm_plot_notes(out$data)
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par()["mar"]
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(8, 5, 3, 1))
    graphics::barplot(
      height = links$N_Common,
      names.arg = links$Pair,
      las = 2,
      col = style$accent_primary,
      main = .mfrm_plot_display_title(out$data),
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
              n_flagged = if (anyNA(object$element_detail$Flag)) NA_integer_ else sum(object$element_detail$Flag),
              n_flagged_available = sum(object$element_detail$Flag, na.rm = TRUE),
              n_unavailable = sum(is.na(object$element_detail$Flag)),
              config = object$config)
  class(out) <- "summary.mfrm_equating_chain"
  out
}

#' @rdname build_equating_chain
#' @export
print.summary.mfrm_equating_chain <- function(x, ...) {
  cat("--- Screened Linking Chain ---\n")
  print_linking_uncertainty_note()
  print_wrapped_line("Offset_SD is residual spread, not the SE of the estimated offset; cumulative offsets omit uncertainty across links.")
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
  if (isTRUE(x$n_flagged > 0)) cat("\nFlagged linking elements:", x$n_flagged, "\n")
  if (isTRUE(x$n_unavailable > 0)) cat("\nUnavailable linking-element flags:", x$n_unavailable,
    "| Flagged among available:", x$n_flagged_available, "\n")
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
    Scope = c("RSM / PCM", "GPCM"),
    Status = c("supported", "supported_with_caveat"),
    Note = c(
      "Supported as a synthesis layer over documented anchor-review, drift, and equating-chain objects.",
      if (gpcm_detected) {
        "Supported with caveat: GPCM source objects are summarized as exploratory anchor/drift/chain evidence, not an operational linking decision."
      } else {
        "Supported with caveat when GPCM source objects are supplied; not active for this RSM/PCM review."
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
  drift_summary$N_Unavailable <- drift_summary[["N_Unavailable"]] %||% rep(0L, nrow(drift_summary))
  drift_threshold <- as.numeric(drift$config$drift_threshold %||% NA_real_)
  if (nrow(drift_summary) > 0 &&
      all(c("Facet", "Wave", "N_Flagged", "Max_Drift") %in% names(drift_summary))) {
    flagged_rows <- drift_summary |>
      dplyr::filter(.data$N_Flagged > 0 | .data$Max_Drift > drift_threshold | .data$N_Unavailable > 0) |>
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
        Signal = dplyr::if_else(.data$N_Unavailable > 0,
          "Unavailable drift comparisons; inspect source estimates and offset status",
          "Flagged element drift in fitted-wave comparison"),
        Magnitude = as.numeric(.data$Max_Drift),
        SeverityGroup = "review",
        ReviewPriority = pmax(as.numeric(.data$N_Flagged), as.numeric(.data$Max_Drift),
          as.numeric(.data$N_Unavailable), na.rm = TRUE),
        Guidance = "Inspect flagged rows in detect_anchor_drift() and review wave-specific drift plots.",
        PrimaryPlotRoute = "plot(drift, type = \"drift\")",
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
        PrimaryPlotRoute = "plot(drift, type = \"drift\")",
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
      PrimaryPlotRoute = "plot(chain, type = \"chain\")",
      SupportStatus = "supported"
    ) |>
    .linking_review_standardize()
}

#' Build a linking-review synthesis object
#' @param anchor_review Optional output from [review_mfrm_anchors()].
#' @param drift Optional output from [detect_anchor_drift()].
#' @param chain Optional output from [build_equating_chain()].
#' @param top_n For `build_linking_review()`: Maximum number of linking-risk rows to highlight in summary
#'   outputs. The full object keeps the full risk tables.
#'
#' For `summary()`: Number of top linking-risk rows to keep in the compact summary.
#' @details `build_linking_review()` does not recompute anchor, drift, or chain
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
#' It also does not verify source-fit readiness, cross-wave element identity,
#' invariance, or the external assumptions behind group anchors; these must be
#' established before promoting the synthesized review.
#' @section Recommended input route:
#' Use existing package-native outputs in this order:
#' 1. [review_mfrm_anchors()] for pre-fit anchor adequacy.
#' 2. [detect_anchor_drift()] for direct wave-to-reference drift screening.
#' 3. [build_equating_chain()] for adjacent screened-link review across waves.
#' @section Interpreting output:
#' - `overview`: which evidence sources were supplied and the current review status.
#' - `top_linking_risks`: primary operational triage table.
#' - `group_view_index`: stable wave/link/facet/source-family grouping routes.
#' - `plot_map`: which plotting function to use for each result.
#' - `reporting_map`: what is covered here versus which manuscript-oriented
#'   helper should be used separately.
#' @section GPCM boundary:
#' This helper summarizes the documented `RSM` / `PCM` linking workflow.
#' When supplied drift/chain sources resolve to `GPCM`, it retains an explicitly
#' exploratory review and marks that scope as `supported_with_caveat`.
#' This does not establish score interchangeability, invariance or a common
#' scale across administrations.
#' @return An object of class `mfrm_linking_review`.
#' @seealso [review_mfrm_anchors()], [detect_anchor_drift()],
#'   [build_equating_chain()], [plot(x)][plot.mfrm_anchor_drift], [mfrmr_linking_and_dff]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("build_linking_review", run.dontrun = TRUE)
#' \dontrun{
#' # Deliberately linked teaching waves: common labels below represent the
#' # same rater and criterion identities by construction.
#' toy <- load_mfrmr_data("example_core")
#' people <- unique(toy$Person)
#' # Two balanced 24-Person waves retain every Rater and Criterion.
#' d1 <- toy[toy$Person %in% people[1:24], , drop = FALSE]
#' d2 <- toy[toy$Person %in% people[25:48], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", quad_points = 7, maxit = 30)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", quad_points = 7, maxit = 30)
#' anchor_review_obj <- review_mfrm_anchors(d1, "Person", c("Rater", "Criterion"), "Score")
#' drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
#' chain <- build_equating_chain(list(Wave1 = fit1, Wave2 = fit2))
#' review <- build_linking_review(anchor_review = anchor_review_obj, drift = drift, chain = chain)
#' summary(review)
#' review$top_linking_risks
#' review$group_view_index
#' }
#' @export
#' @order 1
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
    GPCMSupport = as.character(support_status$Status[support_status$Scope == "GPCM"][1] %||% NA_character_)
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
      "Inspect detect_anchor_drift() and plot(drift, type = \"drift\") for wave-level follow-up."
    )
  }
  if (nrow(chain_risks) > 0) {
    next_actions <- c(
      next_actions,
      "Inspect build_equating_chain() and plot(chain, type = \"chain\") before using cumulative offsets operationally."
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
    "GPCM" = as.character(support_status$Status[support_status$Scope == "GPCM"][1] %||% NA_character_)
  )

  plot_map <- tibble::tibble(
    ReviewArea = c("Anchor adequacy", "Wave-level drift", "Screened chain"),
    Available = c(!is.null(anchor_review), !is.null(drift), !is.null(chain)),
    PlotHelper = c(
      "plot(anchor_review, type = \"issue_counts\")",
      "plot(drift, type = \"drift\")",
      "plot(chain, type = \"chain\")"
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
      "detect_anchor_drift() / plot(x)",
      "build_equating_chain() / plot(type = \"chain\")",
      "build_summary_table_bundle(summary(linking_review)) / reporting_checklist()"
    )
  )

  notes <- linking_review_notes(gpcm_detected)

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
        "Score-side scorefile export under GPCM",
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

#' @rawNamespace S3method(print,mfrm_linking_review,print.mfrm_bundle)
NULL

#' @rdname build_linking_review
#' @details `summary()` method: Summarize a linking-review object
#' @param object Output from [build_linking_review()].
#' @param digits Number of digits for printed numeric values.
#' @param ... Reserved for generic compatibility.
#' @return For `summary()`: An object of class `summary.mfrm_linking_review`.
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
  print_linking_uncertainty_note()
  if (nrow(x$overview) > 0) {
    cat("\nOverview\n")
    keep <- intersect(c("AnchorReviewAvailable", "DriftAvailable", "ChainAvailable", "SourceModels"), names(x$overview))
    print(as.data.frame(x$overview[, keep, drop = FALSE]), row.names = FALSE)
  }
  print_bullet_section("Key Warnings", x$key_warnings)
  print_bullet_section("Next Actions", x$next_actions)
  if (nrow(x$top_linking_risks) > 0) {
    cat("\nTop Linking Risks\n")
    keep <- intersect(c("Facet", "Level", "Wave", "Link", "Signal", "Magnitude", "Guidance"), names(x$top_linking_risks))
    print(round_numeric_df(as.data.frame(x$top_linking_risks[, keep, drop = FALSE]), digits = digits), row.names = FALSE)
  }
  if (nrow(x$plot_routes) > 0) {
    cat("\nPlot Follow-up\n")
    print(as.data.frame(x$plot_routes), row.names = FALSE)
  }
  print_wrapped_line("Detailed risk values, grouped views and support information remain in the returned review tables.")
  print_bullet_section("Notes", linking_review_notes(any(x$settings$source_models == "GPCM")))
  invisible(x)
}

# --- build_misfit_casebook ---------------------------------------------------


.misfit_casebook_support_status <- function(model) {
  model <- as.character(model %||% NA_character_)[1]
  gpcm <- identical(model, "GPCM")
  tibble::tibble(
    Scope = c("RSM / PCM", "GPCM"),
    Status = c(
      "supported",
      if (gpcm) "supported_with_caveat" else "deferred"
    ),
    Note = c(
      "Supported as a synthesis layer over package-native screening outputs.",
      if (gpcm) {
        paste(
          "Supported with caveat: GPCM casebook rows inherit",
          "exploratory screening semantics from residual and strict marginal sources."
        )
      } else {
        "Deferred unless a GPCM source fit is supplied."
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
            "Strict marginal cell screening is available, but GPCM rows remain",
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
            "Strict pairwise screening is available, but GPCM rows remain",
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
          "Unexpected-response rows are available, but GPCM interpretation remains operational rather than inferential."
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
          "Displacement rows are available, but GPCM interpretation remains operational rather than inferential."
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

  .linking_review_group_view_index_empty() |>
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
      PrimaryPlotRoute = "plot(diagnostics, type = \"marginal_fit\", draw = FALSE)",
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
      PrimaryPlotRoute = "plot(diagnostics, type = \"marginal_pairwise\", draw = FALSE)",
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
      PrimaryPlotRoute = "plot(unexpected, type = \"unexpected\", draw = FALSE)",
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
      PrimaryPlotRoute = "plot(fit, type = \"qc\", diagnostics = diagnostics, draw = FALSE)",
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
      PrimaryPlotRoute = "plot(displacement, type = \"displacement\", draw = FALSE)",
      SupportStatus = support_status
    )
  .misfit_casebook_standardize(out)
}

#' Build a case-level misfit review bundle
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#'   Must match this fitted analysis; see [mfrmr_compatibility_layer] for
#'   migration of saved diagnostics.
#' @param unexpected Optional output from [unexpected_response_table()].
#' @param displacement Optional output from [displacement_table()].
#' @param administration_id Optional scalar identifier describing the current
#'   administration or form. It is stored in row-level provenance and summary
#'   outputs when supplied.
#' @param wave_id Optional scalar identifier for the current wave or occasion.
#'   It is stored in row-level provenance and summary outputs when supplied.
#' @param top_n For `build_misfit_casebook()`: Maximum number of rows to keep in compact summary outputs.
#'
#' For `summary()`: Number of top case rows to keep in the compact summary.
#' @details `build_misfit_casebook()` is a synthesis layer over package-native screening
#' outputs. It does not invent a new misfit statistic. Instead, it organizes
#' existing evidence families into one case-level review surface:
#'
#' - element-level Infit / Outfit MnSq misfit from `diagnostics$fit`
#'   (rows whose Infit or Outfit MnSq falls outside the configured Linacre
#'   heuristic review band, 0.5-1.5 by default)
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
#' @section Recommended input route:
#' 1. Fit with [fit_mfrm()].
#' 2. Build diagnostics with [diagnose_mfrm()].
#' 3. Optionally build [unexpected_response_table()] and [displacement_table()]
#'    yourself when you want custom thresholds before synthesizing the casebook.
#' @section GPCM boundary:
#' For `GPCM`, the helper is available with caveat. The casebook inherits
#' exploratory screening semantics from the underlying residual and strict
#' marginal sources; it should not be read as a formal inferential case test.
#' @return An object of class `mfrm_misfit_casebook`.
#' @seealso [diagnose_mfrm()], [unexpected_response_table()],
#'   [displacement_table()], [plot(x, type = "unexpected")][plot.mfrm_fit], [plot(x, type = "displacement")][plot.mfrm_fit],
#'   [plot(x, type = "marginal_fit")][plot.mfrm_fit], [plot(x, type = "marginal_pairwise")][plot.mfrm_fit]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("build_misfit_casebook", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' # A balanced slice retains every Rater and Criterion while running quickly.
#' toy <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", model = "RSM", quad_points = 5)
#' diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
#' casebook <- build_misfit_casebook(fit, diagnostics = diag, top_n = 10)
#' summary(casebook)
#' casebook$top_cases
#' }
#' @section Saved-result reuse:
#' Supplied screening results must describe this fitted analysis. Results from
#' another fit or older results without identity metadata are refused. See
#' [mfrmr_compatibility_layer] to recreate them from the saved fit and matching
#' diagnostics without refitting the model.
#' @export
#' @order 1
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
  diagnostics <- .validate_linking_review_input(diagnostics, "diagnostics", "mfrm_diagnostics")
  unexpected <- .validate_linking_review_input(unexpected, "unexpected", "mfrm_unexpected")
  displacement <- .validate_linking_review_input(displacement, "displacement", "mfrm_displacement")
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
  mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "build_misfit_casebook()")
  identity <- mfrm_result_fit_identity(fit)
  if (!is.null(unexpected)) {
    mfrm_validate_result_identity(unexpected, identity, "unexpected",
                                 "unexpected_response_table()", "build_misfit_casebook()")
  }
  if (!is.null(displacement)) {
    mfrm_validate_result_identity(displacement, identity, "displacement",
                                 "displacement_table()", "build_misfit_casebook()")
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
    next_actions <- c(next_actions, "Use plot(diagnostics, type = \"marginal_fit\", draw = FALSE) to inspect the largest first-order strict marginal cells.")
  }
  if (nrow(pairwise_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot(diagnostics, type = \"marginal_pairwise\", draw = FALSE) to inspect the strongest pairwise local-dependence signals.")
  }
  if (nrow(unexpected_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot(unexpected, type = \"unexpected\", draw = FALSE) to review the most surprising person-level observations.")
  }
  if (nrow(displacement_cases) > 0L) {
    next_actions <- c(next_actions, "Use plot(displacement, type = \"displacement\", draw = FALSE) when flagged facet levels suggest anchor or stability review.")
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
    "GPCM" = as.character(support_status$Status[support_status$Scope == "GPCM"][1] %||% NA_character_)
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
      "plot(diagnostics, type = \"marginal_fit\", draw = FALSE)",
      "plot(diagnostics, type = \"marginal_pairwise\", draw = FALSE)",
      "plot(unexpected, type = \"unexpected\", draw = FALSE)",
      "plot(displacement, type = \"displacement\", draw = FALSE)"
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
      "plot(x, type = \"marginal_fit\") / plot(x, type = \"marginal_pairwise\")",
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
    GPCMSupport = as.character(support_status$Status[support_status$Scope == "GPCM"][1] %||% NA_character_)
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

#' @rawNamespace S3method(print,mfrm_misfit_casebook,print.mfrm_bundle)
NULL

#' @rdname build_misfit_casebook
#' @details `summary()` method: Summarize a misfit-casebook object
#' @param object Output from [build_misfit_casebook()].
#' @param digits Number of digits for printed numeric values.
#' @param ... Reserved for generic compatibility.
#' @return For `summary()`: An object of class `summary.mfrm_misfit_casebook`.
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
    Scope = c("RSM / PCM reference", "GPCM comparison"),
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
        .data$DeltaEstimate > 0 ~ "Higher in GPCM",
        .data$DeltaEstimate < 0 ~ "Lower in GPCM",
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

.weighting_review_comparison_contract <- function(rasch_fit,
                                                  gpcm_fit,
                                                  comparison,
                                                  reference_model,
                                                  aligned_pcm_owner,
                                                  pcm_gpcm_lrt) {
  basis <- comparison$comparison_basis %||% list()
  comparison_tbl <- as.data.frame(
    comparison$table %||% data.frame(),
    stringsAsFactors = FALSE
  )
  reference_method <- public_mfrm_method_label(
    as.character(rasch_fit$config$method %||% NA_character_)[1]
  )
  comparison_method <- public_mfrm_method_label(
    as.character(gpcm_fit$config$method %||% NA_character_)[1]
  )
  same_method <- identical(reference_method, comparison_method)
  both_jml <- same_method && identical(reference_method, "JML")
  both_mml <- same_method && identical(reference_method, "MML")
  same_data <- isTRUE(basis$same_data)
  both_ready <- isTRUE(basis$all_inference_ready)
  both_numerical <- all(vapply(list(rasch_fit, gpcm_fit), function(fit) {
    identical(as.character(mfrmr_get_readiness_record(fit)$fit$NumericalState[1]), "ready")
  }, logical(1)))
  ic_selectable <- isTRUE(basis$all_ic_selectable)
  ic_available <- both_mml && same_data && isTRUE(basis$ic_comparable) &&
    nrow(comparison_tbl) == 2L && "ICFitEligible" %in% names(comparison_tbl) &&
    all(!is.na(comparison_tbl$ICFitEligible) & comparison_tbl$ICFitEligible)

  loglik_difference <- NA_real_
  if (nrow(comparison_tbl) >= 2L && "LogLik" %in% names(comparison_tbl)) {
    loglik <- suppressWarnings(as.numeric(comparison_tbl$LogLik[seq_len(2L)]))
    if (length(loglik) == 2L && all(is.finite(loglik))) {
      loglik_difference <- loglik[2] - loglik[1]
    }
  }

  evidence_tier <- if (ic_available) {
    "mml_information_criterion_comparison"
  } else if (both_jml && same_data && !both_numerical) {
    "jml_numerical_review_only"
  } else if (both_jml && same_data) {
    "jml_descriptive_reweighting_only"
  } else if (both_mml && same_data && !both_numerical) {
    "mml_numerical_review_only"
  } else if (both_mml && same_data && !both_ready) {
    "mml_descriptive_not_inference_ready"
  } else if (both_mml && same_data && !ic_selectable) {
    "mml_screening_or_review_grid_only"
  } else if (both_mml && same_data) {
    "mml_noncomparable_descriptive_only"
  } else {
    "cross_method_or_data_not_comparable"
  }

  loglik_status <- if (!is.finite(loglik_difference)) {
    "unavailable"
  } else if (!both_numerical) {
    "numerical_review_required"
  } else if (ic_available) {
    "comparable_mml_likelihood"
  } else if (both_jml) {
    "descriptive_unpenalized_gain_not_selection"
  } else {
    "descriptive_noncomparable"
  }

  selection_route <- if (ic_available) {
    "compare_mfrm_information_criteria"
  } else if (both_jml) {
    "withheld_JML_has_no_automatic_PCM_GPCM_selection"
  } else {
    "withheld_comparison_basis_not_selectable"
  }

  facets_role <- if (identical(reference_model, "PCM") && both_jml &&
                      aligned_pcm_owner) {
    "PCM_JML_side_only_no_FACETS_free_slope_GPCM_counterpart"
  } else if (identical(reference_model, "PCM") && both_jml) {
    "PCM_JML_side_only_owner_mismatch_blocks_reduction"
  } else {
    "not_an_internal_free_slope_comparator"
  }

  recommended_use <- paste(
    "Inspect changes in facet measures, relative slopes and information shares.",
    if (ic_available) {
      "Information criteria compare the supplied MML candidates; they do not select an operational scoring policy."
    } else {
      "The supplied fits do not satisfy the information-criterion comparison checks."
    },
    if (identical(pcm_gpcm_lrt, "computed")) {
      "The PCM/GPCM LRT tests equal relative slopes; it does not decide the scoring policy."
    } else if (isTRUE(basis$nested_requested)) {
      paste("The requested likelihood-ratio test was not computed:", basis$lrt_reason)
    } else if (!both_mml || !identical(reference_model, "PCM") || !aligned_pcm_owner) {
      "An equal-slope PCM/GPCM test requires two MML fits with the same step facet and population model."
    } else {
      "For an equal-slope test, set nested = TRUE and inspect the separate LRT status and reason."
    }
  )
  numerical_note <- if (both_numerical) {
    "Numerical convergence checks passed for both fits; inference eligibility is assessed separately."
  } else {
    "Numerical convergence is incomplete or requires review for at least one fit; interpret differences cautiously."
  }

  tibble::tibble(
    ReferenceModel = reference_model,
    ComparisonModel = as.character(gpcm_fit$config$model %||% NA_character_)[1],
    ReferenceMethod = reference_method,
    ComparisonMethod = comparison_method,
    SameMethod = same_method,
    SamePreparedData = same_data,
    BothInferenceReady = both_ready,
    BothNumericallyReady = both_numerical,
    NumericalReview = numerical_note,
    EvidenceTier = evidence_tier,
    FormalModelSelectionAvailable = ic_available,
    SelectionRoute = selection_route,
    ObservedLogLikDifference = loglik_difference,
    LogLikDifferenceStatus = loglik_status,
    AICPreferred = if (ic_available) comparison$preferred$AIC %||% NA_character_ else NA_character_,
    PersonBICPreferred = if (ic_available) comparison$preferred$BIC %||% NA_character_ else NA_character_,
    SABICPreferred = if (ic_available) comparison$preferred$SABIC %||% NA_character_ else NA_character_,
    PCMvsGPCMLRT = pcm_gpcm_lrt,
    FACETSComparisonRole = facets_role,
    RecommendedUse = recommended_use
  )
}

#' Build a weighting-policy review between Rasch-family and GPCM fits
#' @param rasch_fit Output from [fit_mfrm()] using `model = "RSM"` or `"PCM"`.
#' @param gpcm_fit Output from [fit_mfrm()] using `model = "GPCM"`.
#' @param theta_range Numeric vector of length 2 passed to [compute_information()]
#'   for the information-redistribution comparison.
#' @param theta_points Integer number of theta grid points passed to
#'   [compute_information()].
#' @param top_n For `build_weighting_review()`: Maximum number of rows to keep in compact summary outputs.
#'
#' For `summary()`: Number of top rows to retain in compact summary tables.
#' @param nested Request the PCM/GPCM equal-slope likelihood-ratio test. Default
#'   `FALSE`. Requires matched MML fits and the checks in [compare_mfrm()].
#' @details `build_weighting_review()` is an operational model-choice review helper. It
#' is designed for the common question:
#'
#' - what changes when a Rasch-family equal-weighting model is replaced with a
#'   `GPCM` that allows discrimination-based reweighting?
#'
#' The helper does not estimate a new model. Instead, it synthesizes four
#' package-native evidence sources:
#'
#' - [compare_mfrm()] for same-data model comparison
#' - the non-person facet measures from each fit
#' - the `GPCM` slope table
#' - [compute_information()] for design-weighted information redistribution
#'
#' The result is intended for substantive review, not for automatic model
#' selection. In particular, a better-fitting `GPCM` should not by itself be
#' interpreted as a reason to discard an equal-weighting Rasch-family route.
#' The fitted GPCM contains one slope for every level of one designated facet,
#' not one common slope and not simultaneous criterion-by-rater slope blocks.
#' The overview records the slope owner, step owner, level count, free relative
#' slope contrasts, and whether the supplied reference is the exact unit-slope
#' PCM response-kernel reduction. MML information-criterion ranking requires
#' the likelihood and local-solution checks in [compare_mfrm()]. With
#' `nested = TRUE`, a PCM/GPCM asymptotic chi-square LRT additionally requires
#' matching population, step and facet settings and G-1 free slope contrasts.
#' The returned `comparison_contract` records the comparison and test status;
#' observed changes in scores and information need substantive interpretation. A JML
#' log-likelihood increase is not promoted to automatic PCM-versus-GPCM model
#' selection because it is unpenalized and the GPCM contains additional slope
#' parameters. FACETS may serve as a direct comparator for the PCM/JML side
#' only; its post-fit discrimination statistic is not a jointly estimated
#' free-slope GPCM counterpart.
#' @section Recommended input route:
#' 1. Fit an equal-weighting reference model with `model = "RSM"` or `"PCM"`.
#' 2. Fit a `GPCM` on the same prepared response data.
#' 3. Run `build_weighting_review(rasch_fit, gpcm_fit)`.
#' 4. Read `summary(review)` before deciding whether the discrimination-based
#'    reweighting is substantively acceptable.
#' @section What the returned tables mean:
#' - `model_comparison`: same-data model-comparison bundle from [compare_mfrm()].
#'   AIC/Person-BIC/SABIC ranking is available only when `ICComparable` is true.
#'   Inspect `$lrt` and `$comparison_basis$lrt_reason` for a requested test.
#' - `comparison_contract`: one-row evidence-tier table stating whether formal
#'   model selection is available, how any observed log-likelihood difference
#'   may be read, and the bounded role of FACETS in a JML review.
#' - `facet_shift`: how non-person facet estimates move under `GPCM`.
#' - `slope_profile`: which `slope_facet` levels are upweighted or downweighted.
#' - `information_redistribution`: within-facet information-share changes
#'   between the Rasch-family fit and `GPCM`.
#' - `top_reweighted_levels`: compact triage table for the strongest
#'   slope-facet-level redistribution signals.
#' @section GPCM boundary:
#' This helper is available only for the current `GPCM` branch. It
#' requires `slope_facet == step_facet` even though MML fitting permits separate owners, and
#' should be read as an operational weighting-policy review, not as a formal
#' validity adjudication.
#' @return An object of class `mfrm_weighting_review`.
#' @seealso [compare_mfrm()], [compute_information()], [gpcm_capability_matrix()]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("build_weighting_review", run.dontrun = TRUE)
#' \dontrun{
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
#' @order 1
build_weighting_review <- function(rasch_fit,
                                   gpcm_fit,
                                   theta_range = c(-6, 6),
                                   theta_points = 101L,
                                   top_n = 10L,
                                   nested = FALSE) {
  stop_if_product_slopes(gpcm_fit, "build_weighting_review()")
  rasch_fit <- .validate_weighting_review_fit(rasch_fit, "rasch_fit", c("RSM", "PCM"))
  gpcm_fit <- .validate_weighting_review_fit(gpcm_fit, "gpcm_fit", "GPCM")
  top_n <- max(1L, as.integer(top_n %||% 10L))
  theta_points <- max(11L, as.integer(theta_points %||% 101L))

  slope_facet <- as.character(gpcm_fit$config$slope_facet %||% NA_character_)[1]
  step_facet <- as.character(gpcm_fit$config$step_facet %||% NA_character_)[1]
  if (!identical(slope_facet, step_facet)) {
    stop(
      "build_weighting_review() currently requires the `GPCM` branch with `slope_facet == step_facet`.",
      call. = FALSE
    )
  }

  comparison <- suppressWarnings(compare_mfrm(
    rasch_fit,
    gpcm_fit,
    labels = c(
      paste0(as.character(rasch_fit$config$model %||% "RSM")[1], " reference"),
      "GPCM"
    ),
    warn_constraints = FALSE,
    nested = nested
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
  reference_model <- as.character(
    rasch_fit$config$model %||% NA_character_
  )[1]
  reference_step_facet <- as.character(
    rasch_fit$config$step_facet %||% NA_character_
  )[1]
  aligned_pcm_owner <- identical(reference_model, "PCM") &&
    !is.na(reference_step_facet) && nzchar(reference_step_facet) &&
    identical(reference_step_facet, step_facet) &&
    identical(step_facet, slope_facet)
  pcm_gpcm_lrt <- if (!identical(reference_model, "PCM")) {
    "not_applicable_reference_is_not_pcm"
  } else if (aligned_pcm_owner) {
    as.character(basis$lrt_status %||% "not_requested")
  } else {
    "not_available_owner_mismatch"
  }
  comparison_contract <- .weighting_review_comparison_contract(
    rasch_fit = rasch_fit,
    gpcm_fit = gpcm_fit,
    comparison = comparison,
    reference_model = reference_model,
    aligned_pcm_owner = aligned_pcm_owner,
    pcm_gpcm_lrt = pcm_gpcm_lrt
  )
  slope_level_count <- nrow(slope_profile)
  overview <- tibble::tibble(
    ReferenceModel = reference_model,
    ComparisonModel = as.character(gpcm_fit$config$model %||% NA_character_)[1],
    ReferenceMethod = public_mfrm_method_label(as.character(rasch_fit$config$method %||% NA_character_)[1]),
    ComparisonMethod = public_mfrm_method_label(as.character(gpcm_fit$config$method %||% NA_character_)[1]),
    SlopeFacet = slope_facet,
    StepFacet = step_facet,
    ReferenceStepFacet = reference_step_facet,
    AlignedStepSlopeOwner = aligned_pcm_owner,
    SlopeLevelCount = as.integer(slope_level_count),
    FreeRelativeSlopeContrasts = as.integer(max(slope_level_count - 1L, 0L)),
    SlopeComposition = "single_aligned_facet_level_specific_relative_slopes",
    SimultaneousCriterionRaterSlopeBlocks = FALSE,
    UnitSlopePCMReduction = aligned_pcm_owner,
    PCMvsGPCMLRT = pcm_gpcm_lrt,
    EvidenceTier = comparison_contract$EvidenceTier,
    FormalModelSelectionAvailable = comparison_contract$FormalModelSelectionAvailable,
    ObservedLogLikDifference = comparison_contract$ObservedLogLikDifference,
    LogLikDifferenceStatus = comparison_contract$LogLikDifferenceStatus,
    ReviewStatus = review_status,
    ComparisonMode = comparison_mode,
    MaxAbsLogSlope = max_abs_log_slope,
    MaxAbsInfoShareDelta = max_abs_info_delta
  )

  status <- make_summary_block(
    "Overall status" = review_status,
    "Weighting principle" = "Rasch-family equal weighting vs GPCM discrimination-based reweighting",
    "Slope ownership" = paste0(
      slope_facet, " levels; one aligned slope/step owner only"
    ),
    "Comparison basis" = comparison_mode,
    "Evidence tier" = comparison_contract$EvidenceTier,
    "Formal model selection" = if (isTRUE(
      comparison_contract$FormalModelSelectionAvailable
    )) "available" else "withheld",
    "PCM vs GPCM LRT" = pcm_gpcm_lrt
  )

  key_warnings <- character(0)
  if (!isTRUE(basis$ic_comparable)) {
    key_warnings <- c(
      key_warnings,
      "Model-comparison weights are descriptive only because the two fits do not share a fully comparable formal MML basis."
    )
  }
  if (identical(comparison_contract$ReferenceMethod, "JML") &&
      identical(comparison_contract$ComparisonMethod, "JML")) {
    key_warnings <- c(
      key_warnings,
      paste(
        "The JML log-likelihood difference is unpenalized and cannot select",
        "PCM versus GPCM. Review numerical convergence separately from",
        "the conditions needed for statistical inference."
      )
    )
  }
  if (!aligned_pcm_owner && identical(reference_model, "PCM")) {
    key_warnings <- c(
      key_warnings,
      paste(
        "PCM and GPCM do not share the same explicit aligned step/slope owner;",
        "the unit-slope response-kernel reduction and a PCM-versus-GPCM LRT",
        "are unavailable for this pair."
      )
    )
  }
  if (nrow(slope_profile) > 0L) {
    lead_slope <- slope_profile[1, , drop = FALSE]
    key_warnings <- c(
      key_warnings,
      paste0(
        "Largest GPCM slope deviation is at ",
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
    "Use IC ranking only when ICComparable is true; weighting consequences require a separate substantive decision.",
    paste0("Use slope_profile and top_reweighted_levels to inspect whether ", slope_facet, " levels are being upweighted or downweighted in substantively acceptable ways."),
    paste0("Use plot_information(compute_information(rasch_fit), type = \"iif\", facet = \"", slope_facet, "\", draw = FALSE) and the GPCM analogue to inspect precision redistribution visually."),
    "If equal contributions of items and raters are part of the score interpretation, retain the Rasch-family fit as the operational reference even when GPCM fits better."
  ), max_n = 4L)

  plot_map <- tibble::tibble(
    ReviewArea = c("Reference precision", "GPCM precision", "Fit comparison"),
    Available = c(TRUE, TRUE, TRUE),
    PlotHelper = c(
      paste0("plot_information(compute_information(rasch_fit), type = \"iif\", facet = \"", slope_facet, "\", draw = FALSE)"),
      paste0("plot_information(compute_information(gpcm_fit), type = \"iif\", facet = \"", slope_facet, "\", draw = FALSE)"),
      "summary(compare_mfrm(rasch_fit, gpcm_fit))"
    ),
    Trigger = c(
      "Use to inspect the equal-weighting reference precision split across slope-facet levels.",
      "Use to inspect how GPCM redistributes precision across the same levels.",
      "Use to review q>=31 comparable AIC/Person-BIC/SABIC and relative candidate-set weights before making a model-choice argument; close decisions need a denser common-grid check."
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
      "fair_average_table() for the Rasch-family route; keep GPCM score semantics separate"
    )
  )

  notes <- clean_summary_lines(c(
    "Observation weights and discrimination-based reweighting are separate concepts in this package.",
    paste0(
      "The fitted slopes vary across levels of `", slope_facet,
      "`; other facets have no separate slope block."
    ),
    "Criterion-owned and rater-owned GPCM fits are separate restricted models; both blocks cannot be estimated together by the current GPCM interface.",
    "FACETS is a direct JML comparator only for the aligned equal-discrimination PCM side; its reported discrimination is a post-fit diagnostic, not the fitted free-slope GPCM parameter.",
    "The review is intended to make reweighting visible; it does not decide by itself whether GPCM should replace the Rasch-family operational model.",
    "Information-share changes are computed within each facet because the same total information is partitioned separately by facet."
  ))

  out <- list(
    overview = overview,
    status = status,
    key_warnings = key_warnings,
    next_actions = next_actions,
    model_comparison = comparison,
    comparison_contract = comparison_contract,
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
      intended_use = "weighting_policy_review",
      nested = isTRUE(nested)
    )
  )
  as_mfrm_bundle(out, "mfrm_weighting_review")
}

validate_model_review_convergence <- function(x) {
  weighting <- inherits(x, c("mfrm_weighting_review", "summary.mfrm_weighting_review"))
  table <- if (weighting) x$comparison_contract else x$model_roles
  field <- if (weighting) "BothNumericallyReady" else "NumericalConvergence"
  if (!is.data.frame(table) || !field %in% names(table)) {
    stop(paste0(
      "This saved review does not distinguish numerical convergence from inference eligibility. ",
      "Recreate it with ", if (weighting) "build_weighting_review()" else "build_model_choice_review()",
      " using the existing fits; no model refit is needed."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' @rawNamespace S3method(print,mfrm_weighting_review,print.mfrm_bundle)
NULL

#' @rdname build_weighting_review
#' @details `summary()` method: Summarize a weighting-review object
#' @param object Output from [build_weighting_review()].
#' @param digits Number of digits for printed numeric values.
#' @param ... Reserved for generic compatibility.
#' @return For `summary()`: An object of class `summary.mfrm_weighting_review`, including the
#'   evidence-tier `comparison_contract` table and any computed `lrt`.
#' @seealso [build_weighting_review()]
#' @export
summary.mfrm_weighting_review <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_weighting_review")) {
    stop("`object` must be output from build_weighting_review().", call. = FALSE)
  }

  validate_model_review_convergence(object)
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  out <- list(
    lrt = object$model_comparison$lrt,
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    status = tibble::as_tibble(object$status %||% tibble::tibble()),
    comparison_contract = tibble::as_tibble(
      object$comparison_contract %||% tibble::tibble()
    ),
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
  validate_model_review_convergence(x)
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Weighting Review Summary\n")
  if (nrow(x$overview) > 0) {
    cat("\nOverview\n")
    keep <- intersect(c("ReferenceModel", "ComparisonModel", "ReferenceMethod",
                        "ComparisonMethod", "SlopeFacet", "StepFacet",
                        "SlopeLevelCount", "FreeRelativeSlopeContrasts"), names(x$overview))
    print(as.data.frame(x$overview[, keep, drop = FALSE]), row.names = FALSE)
  }
  if (nrow(x$comparison_contract) > 0) {
    contract <- x$comparison_contract
    cat("\nComparison interpretation\n")
    if ("NumericalReview" %in% names(contract)) print_wrapped_line(contract$NumericalReview[1])
    cat("  Log-likelihood difference:",
        format(round(contract$ObservedLogLikDifference[1], digits)), "\n")
    print_wrapped_line(contract$RecommendedUse[1])
  }
  if (!is.null(x$lrt)) print_mfrm_lrt(x$lrt)
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
      "One positive relative slope is estimated for every level of one",
      "designated slope facet; MML allows a separate step owner, whereas JML requires",
      "slope_facet == step_facet. The model identifies slopes with geometric mean 1",
      "and does not combine criterion and rater slope blocks."
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
    GPCM = "We fit a generalized partial-credit many-facet model as a slope-aware sensitivity analysis; interpretation focused on whether discrimination-based reweighting changed the substantive conclusions.",
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
    capability <- .gpcm_capability_registry()
    status_for <- function(capability_id) {
      idx <- which(capability$CapabilityID == capability_id)
      if (length(idx) != 1L) {
        stop(
          "Expected one GPCM capability for ID `",
          capability_id,
          "`; found ",
          length(idx),
          ".",
          call. = FALSE
        )
      }
      capability$Status[[idx]]
    }
    return(tibble::tibble(
      FullAPARoute = status_for("apa_export_bundles"),
      ScoreSideExport = status_for("scorefile_export"),
      LinkingSynthesis = status_for("linking_synthesis"),
      RecoveryChecks = status_for("simulation_recovery"),
      FairAverage = status_for("fair_average"),
      BiasScreening = status_for("residual_bias_screening"),
      SummaryAppendix = status_for("summary_appendix"),
      PrimaryHelpers = paste(
        "gpcm_capability_matrix(); build_apa_outputs();",
        "facets_output_file_bundle(include = \"score\"); build_linking_review();",
        "simulate_mfrm_recovery(); fair_average_table(); estimate_bias();",
        "export_summary_appendix()"
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
    sizes <- build_param_sizes(fit$config)
    readiness <- mfrmr_get_readiness_record(fit)$fit
    decision <- mfrm_fit_decision_summary(readiness)
    tibble::tibble(
      Label = labels[[i]],
      Model = model,
      Method = method,
      StepFacet = step_facet,
      StepStructure = if (identical(model, "RSM")) {
        "shared_ladder"
      } else {
        "step_facet_specific_ladders"
      },
      StepCoordinates = nrow(as.data.frame(fit$steps %||% data.frame())),
      FreeStepParameters = as.integer(sizes$steps %||% 0L),
      SlopeFacet = slope_facet,
      SlopeStructure = if (identical(model, "GPCM")) {
        "relative_by_slope_facet_geometric_mean_1"
      } else {
        "fixed_equal_at_1"
      },
      SlopeCoordinates = nrow(as.data.frame(fit$slopes %||% data.frame())),
      FreeSlopeParameters = as.integer(sizes[["log_slopes"]] %||% 0L),
      NumericalConvergence = if (identical(as.character(readiness$NumericalState[1]), "ready")) {
        "Passed"
      } else {
        "Incomplete or requires review"
      },
      FitReadiness = decision$FitReadiness[1L],
      FormalInference = decision$FormalInference[1L],
      Interpretation = decision$Interpretation[1L],
      ReadinessReason = decision$Why[1L],
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
        .data$Model == "GPCM" ~ "Do not write that better fit alone makes GPCM the operational scoring model.",
        .data$Model %in% c("RSM", "PCM") ~ "Do not describe the fit as free-discrimination or slope-weighted.",
        TRUE ~ "Do not omit the model's response-kernel and score-contract assumptions."
      )
    )
}

.model_choice_route_map <- function(has_gpcm) {
  base <- tibble::tibble(
    Question = c(
      "Which model fits better on the same likelihood basis?",
      "Does GPCM change the score interpretation?",
      "Where did precision move under GPCM?",
      "Can I use manuscript APA/report bundles?",
      "Can I use fair averages?",
      "Can I screen bias/interactions?",
      "Can I export a summary appendix?"
    ),
    PrimaryHelper = c(
      "compare_mfrm()",
      "build_weighting_review()",
      "compute_information(); plot(x)",
      "build_apa_outputs(); build_visual_summaries()",
      "fair_average_table()",
      "estimate_bias()",
      "build_summary_table_bundle(); export_summary_appendix()"
    ),
    Interpretation = c(
      "Use q>=31 comparable AIC/Person-BIC/SABIC/logLik as candidate-set fit evidence, not as a standalone scoring decision; close decisions need a denser common-grid check.",
      "Use only when comparing an RSM/PCM reference to GPCM.",
      "Review item/rater/criterion information redistribution before changing the operational model.",
      "Supported for RSM/PCM; available with explicit sensitivity-reporting caveats for GPCM.",
      "Supported for RSM/PCM; supported with explicit SE caveat for GPCM.",
      "Supported for RSM/PCM; conditional screening with profile-likelihood follow-up for GPCM.",
      "Available for direct supported outputs; GPCM fit bundles remain sensitivity-reporting outputs with explicit caveats."
    )
  )
  if (isTRUE(has_gpcm)) {
    base <- dplyr::bind_rows(
      tibble::tibble(
        Question = "What is the exact GPCM support boundary?",
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

#' Build a model-choice review across RSM, PCM, and GPCM fits
#' @param ... For `build_model_choice_review()`: Two or more fitted `mfrm_fit` objects from [fit_mfrm()].
#'
#' For `summary()`: Reserved for generic compatibility.
#' @param labels Optional labels for the supplied fits. If omitted, names from
#'   `...` are used when available; otherwise labels are generated from
#'   model/method combinations.
#' @param run_weighting_review Logical. If `TRUE` and the supplied fits include
#'   at least one `RSM`/`PCM` reference plus one `GPCM` fit, also run
#'   [build_weighting_review()] for the first such pair.
#' @param theta_range,theta_points,top_n Passed to [build_weighting_review()] when
#'   `run_weighting_review = TRUE`.
#' @param warn_constraints Passed to [compare_mfrm()].
#' @details `build_model_choice_review()` is a user-facing synthesis helper. It does not
#' estimate new models. It bundles:
#'
#' - [compare_mfrm()] for verified AIC/Person-BIC/SABIC/log-likelihood
#'   comparison;
#' - comparison-boundary warnings captured from [compare_mfrm()] and retained
#'   in `comparison_warnings` for printing and appendix export;
#' - model-role guidance for `RSM`, `PCM`, and `GPCM`;
#' - reported step/slope coordinate counts, identified free-parameter counts,
#'   and the stored fit-readiness decision for every supplied model;
#' - downstream-route availability for APA output, score-side export, linking,
#'   recovery, fair averages, bias screening, and summary-appendix handoff;
#' - report wording templates that avoid treating better `GPCM` fit as
#'   an automatic operational-scoring decision;
#' - [gpcm_capability_matrix()] when `GPCM` is present;
#' - optionally, [build_weighting_review()] for the first Rasch-family reference
#'   versus `GPCM` pair.
#'
#' The word "bounded" describes the documented model and workflow scope: the
#' package does not implement every possible generalized partial-credit
#' many-facet extension. It does **not** mean that finite optimizer box bounds
#' define the estimator. The current route uses positive slopes, allows separate
#' slope/step owners in MML (JML requires `slope_facet == step_facet`), identifies
#' slopes on the log scale with geometric mean 1, and keeps several downstream score-side/reporting helpers
#' outside the documented boundary.
#' Model-choice ranking also requires the current selectable IC contract:
#' q<31 fits retain raw criteria but produce a screening/review-only bundle.
#' @return An object of class `mfrm_model_choice_review`.
#' @seealso [compare_mfrm()], [build_weighting_review()],
#'   [gpcm_capability_matrix()], [compute_information()]
#' @examples toy <- load_mfrmr_data("example_core")
#' fit_rsm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                     method = "MML", model = "RSM", quad_points = 31)
#' fit_pcm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                     method = "MML", model = "PCM", step_facet = "Criterion",
#'                     quad_points = 31)
#' review <- build_model_choice_review(RSM = fit_rsm, PCM = fit_pcm)
#' summary(review)
#' @export
#' @order 1
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
  comparison_warnings <- character(0)
  comparison <- withCallingHandlers(
    do.call(
      compare_mfrm,
      c(fits, list(labels = labels, warn_constraints = warn_constraints, nested = FALSE))
    ),
    warning = function(w) {
      comparison_warnings <<- c(comparison_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  comparison_warnings <- clean_summary_lines(comparison_warnings)
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
      "GPCM is slope-aware: better fit is sensitivity evidence, not an automatic operational-scoring decision.",
      "GPCM uses one positive slope family; MML permits separate step owners and JML requires slope_facet == step_facet. Weighting reviews and simulation/design workflows require the same owner."
    )
  }
  if (!isTRUE(basis$ic_comparable)) {
    key_warnings <- c(
      key_warnings,
      if (isTRUE(has_gpcm)) {
        "GPCM IC comparison requires its likelihood and solution checks to pass; slope intervals and PCM/GPCM LRT have separate eligibility."
      } else {
        "Information-criterion ranking is unavailable. Review the comparison warnings for data, estimation, integration and inference requirements."
      }
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
      weighting_review_error <- "No RSM/PCM reference plus GPCM pair was available."
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
      "Not requested; set `run_weighting_review = TRUE` for the first RSM/PCM versus GPCM pair."
    }
  )

  notes <- clean_summary_lines(c(
    "This review is a decision aid; it does not refit models or choose an operational model automatically.",
    "Observation weights and GPCM discrimination-based reweighting are separate concepts.",
    "Use GPCM wording only for the current constrained implementation, not for an unrestricted GPCM family claim."
  ))

  out <- list(
    overview = overview,
    key_warnings = key_warnings,
    comparison_warnings = comparison_warnings,
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

#' @rawNamespace S3method(print,mfrm_model_choice_review,print.mfrm_bundle)
NULL

#' @rdname build_model_choice_review
#' @details `summary()` method: Summarize a model-choice review
#' @param object Output from [build_model_choice_review()].
#' @param digits Number of digits for printed numeric values.
#' @return For `summary()`: An object of class `summary.mfrm_model_choice_review`.
#' @seealso [build_model_choice_review()]
#' @export
summary.mfrm_model_choice_review <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_model_choice_review")) {
    stop("`object` must be output from build_model_choice_review().",
         call. = FALSE)
  }
  validate_model_review_convergence(object)
  digits <- max(0L, as.integer(digits %||% 3L))
  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    key_warnings = clean_summary_lines(object$key_warnings %||% character(0)),
    comparison_warnings = clean_summary_lines(object$comparison_warnings %||% character(0)),
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
  validate_model_review_convergence(x)
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("mfrm Model Choice Review\n")
  if (nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    keep <- intersect(c("FitCount", "Models", "OperationalReference", "SensitivityModel"),
                      names(x$overview))
    print(as.data.frame(x$overview[, keep, drop = FALSE]), row.names = FALSE)
  }
  print_bullet_section("Key Warnings", x$key_warnings)
  print_bullet_section("Comparison Warnings", x$comparison_warnings)
  print_bullet_section("Next Actions", x$next_actions)
  if (nrow(x$comparison_table) > 0L) {
    cat("\nComparison Table\n")
    comparison <- round_numeric_df(
      as.data.frame(x$comparison_table),
      digits = digits
    )
    keep <- intersect(
      c(
        "Label", "Model", "Method", "Persons", "Npar", "LogLik",
        "AIC", "BIC", "SABIC", "Delta_AIC", "Delta_BIC",
        "Delta_SABIC", "InferenceReady"
      ),
      names(comparison)
    )
    print(comparison[, keep, drop = FALSE], row.names = FALSE)
    cat("  Full IC audit fields remain available in `$comparison_table`.\n")
  }
  if (nrow(x$model_roles) > 0L) {
    cat("\nModel Contracts and Readiness\n")
    roles <- as.data.frame(x$model_roles, stringsAsFactors = FALSE)
    roles$StepStructure <- dplyr::recode(
      roles$StepStructure,
      shared_ladder = "shared ladder",
      step_facet_specific_ladders = "facet-specific ladders"
    )
    roles$SlopeStructure <- dplyr::recode(
      roles$SlopeStructure,
      fixed_equal_at_1 = "fixed at 1",
      relative_by_slope_facet_geometric_mean_1 =
        "relative by slope facet (GM = 1)"
    )
    cat("\nStep parameters (reported values versus independent parameters)\n")
    steps <- roles[, c(
      "Label", "Model", "StepStructure", "StepCoordinates",
      "FreeStepParameters"
    ), drop = FALSE]
    names(steps)[4:5] <- c("Reported", "Free")
    print(steps, row.names = FALSE)
    cat("\nSlope parameters (reported values versus independent parameters)\n")
    slopes <- roles[, c(
      "Label", "Model", "SlopeStructure", "SlopeCoordinates",
      "FreeSlopeParameters"
    ), drop = FALSE]
    names(slopes)[4:5] <- c("Reported", "Free")
    print(slopes, row.names = FALSE)
    cat("\nFit readiness\n")
    if ("NumericalConvergence" %in% names(roles)) {
      for (i in seq_len(nrow(roles))) {
        print_wrapped_line(paste0(roles$Label[i], ": numerical convergence: ",
                                  roles$NumericalConvergence[i], "."))
      }
    }
    print(roles[, c(
      "Label", "Model", "FormalInference", "Interpretation"
    ), drop = FALSE], row.names = FALSE)
    cat("  Full score contracts and readiness reasons remain in `$model_roles`.\n")
  }
  if (nrow(x$downstream_routes) > 0L) {
    cat("\nDownstream Routes\n")
    keep <- intersect(
      c("Label", "Model", "FullAPARoute", "ScoreSideExport", "LinkingSynthesis",
        "RecoveryChecks", "FairAverage", "BiasScreening", "SummaryAppendix"),
      names(x$downstream_routes)
    )
    routes <- as.data.frame(x$downstream_routes[, keep, drop = FALSE])
    for (field in setdiff(keep, c("Label", "Model"))) {
      routes[[field]] <- dplyr::recode(routes[[field]],
        supported = "Available", supported_with_caveat = "Available with limitations",
        review_required = "Review required", blocked = "Unavailable",
        unsupported = "Unavailable", not_supported = "Unavailable")
    }
    print(routes, row.names = FALSE)
  }
  if (nrow(x$weighting_review_status) > 0L) {
    cat("\nWeighting Review Status\n")
    print(as.data.frame(x$weighting_review_status), row.names = FALSE)
  }
  print_bullet_section("Notes", x$notes)
  invisible(x)
}
