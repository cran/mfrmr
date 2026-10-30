category_model_context <- function(fit) {
  common <- identical(fit$config$model, "RSM")
  data.frame(Model = fit$config$model, Method = fit$config$method,
    StepOwner = if (common) "Common" else fit$config$step_facet,
    StepSelection = fit$config$step_facet_source %||% "not_recorded",
    SlopeOwner = paste(fit$config$slope_facet %||% character(), collapse = ", "),
    StepScale = if (inherits(fit, "mfrm_imported_fit")) "imported_source_threshold" else "centered_step_offset",
    stringsAsFactors = FALSE)
}

category_original_labels <- function(categories, fit) {
  labels <- as.character(categories)
  map <- fit$prep$score_map
  if (is.data.frame(map) && all(c("InternalScore", "OriginalScore") %in% names(map))) {
    idx <- match(categories, map$InternalScore)
    labels[!is.na(idx)] <- as.character(map$OriginalScore[idx[!is.na(idx)]])
  }
  labels
}

category_owner_usage <- function(fit, obs) {
  context <- category_model_context(fit)
  owner <- context$StepOwner
  common <- identical(fit$config$model, "RSM")
  levels <- if (common) "Common" else fit$config$facet_levels[[owner]]
  categories <- seq(fit$prep$rating_min, fit$prep$rating_max)
  owner_known <- common || (owner %in% names(obs) &&
    all(!is.na(obs[[owner]]) & as.character(obs[[owner]]) %in% levels))
  weights <- obs$Weight %||% rep(1, nrow(obs))
  do.call(rbind, lapply(as.character(levels), function(level) {
    rows <- if (common) rep(TRUE, nrow(obs)) else if (owner_known) as.character(obs[[owner]]) == level else rep(FALSE, nrow(obs))
    scores <- obs$Observed[rows]
    score_ok <- owner_known && length(scores) == sum(rows) &&
      all(is.finite(scores) & scores %in% categories)
    weight_ok <- score_ok && all(is.finite(weights[rows]) & weights[rows] >= 0)
    raw <- if (score_ok) vapply(categories, function(k) sum(scores == k), integer(1)) else rep(NA_integer_, length(categories))
    count <- if (weight_ok) vapply(categories, function(k) sum(weights[rows][scores == k]), numeric(1)) else rep(NA_real_, length(categories))
    total <- if (weight_ok) sum(weights[rows]) else NA_real_
    data.frame(StepOwner = owner, StepFacet = level, Category = categories,
      CategoryLabel = category_original_labels(categories, fit), RawCount = raw,
      Count = count, OwnerRatings = if (owner_known) sum(rows) else NA_integer_,
      OwnerWeight = total, Percent = if (is.finite(total) && total > 0) 100 * count / total else NA_real_,
      ZeroCount = count == 0,
      UsageStatus = if (!weight_ok) "unavailable" else if (!sum(rows)) "no_exposure" else
        ifelse(raw == 0, "not_observed", ifelse(count == 0, "zero_weight", "observed")),
      stringsAsFactors = FALSE)
  }))
}

category_threshold_table <- function(fit, diagnostics, usage) {
  steps <- calc_step_order(fit$steps, expected_steps = fit$config$n_cat - 1L,
    expected_facets = fit$config$facet_levels[[fit$config$step_facet %||% ""]])
  coverage <- summarize_threshold_order(steps)
  tab <- as.data.frame(steps)
  if (!nrow(tab)) return(list(table = tab, coverage = coverage))
  tab$GapFromPrev <- tab$Spacing
  tab <- as.data.frame(annotate_threshold_caveats(tab, prep = fit$prep))
  tab$StepOwner <- category_model_context(fit)$StepOwner
  tab$LowerCategoryLabel <- category_original_labels(tab$LowerCategory, fit)
  tab$UpperCategoryLabel <- category_original_labels(tab$UpperCategory, fit)
  if (is.data.frame(usage) && nrow(usage)) {
    for (side in c("Lower", "Upper")) {
      counts <- usage[c("StepFacet", "Category", "Count")]
      names(counts)[2:3] <- c(paste0(side, "Category"), paste0(side, "OwnerCount"))
      tab <- dplyr::left_join(tab, counts, by = c("StepFacet", paste0(side, "Category")))
    }
    tab$OwnerSupportStatus <- ifelse(is.na(tab$LowerOwnerCount) | is.na(tab$UpperOwnerCount),
      "unavailable", ifelse(tab$LowerOwnerCount == 0 | tab$UpperOwnerCount == 0,
        "zero_count_adjacent_category", "observed_adjacent_categories"))
    weak <- tab$OwnerSupportStatus == "zero_count_adjacent_category"
    tab$WeaklyIdentified[weak] <- TRUE
    tab$ThresholdCaveat[weak] <- "Adjacent category has zero weighted count within this step owner; review exposure and category support."
  }

  tab$SE <- tab$CI_Lower <- tab$CI_Upper <- rep(NA_real_, nrow(tab))
  tab$SE_Status <- "not_supplied"
  tab$CIEligible <- FALSE
  tab$CIUse <- "unavailable"
  tab$UncertaintySourceStatus <- "not_supplied"
  uncertainty <- diagnostics$parameter_uncertainty
  se <- as.data.frame(uncertainty$steps %||% data.frame())
  if (nrow(se)) {
    source_verified <- !is.null(uncertainty$source)
    if (source_verified && !identical(uncertainty$source, mfrm_native_location_identity(fit))) {
      stop("fit/diagnostics mismatch: stored step uncertainty belongs to a different calibration.", call. = FALSE)
    }
    if (!"StepFacet" %in% names(se) && identical(fit$config$model, "RSM")) se$StepFacet <- "Common"
    keys <- c("StepFacet", "Step")
    if (!all(c(keys, "Estimate") %in% names(se)) || anyNA(se[keys]) || anyDuplicated(se[keys]) ||
        anyDuplicated(tab[keys]) || nrow(dplyr::anti_join(se, tab, by = keys))) {
      stop("Stored step uncertainty must have unique matching owner/level/step keys.", call. = FALSE)
    }
    idx <- dplyr::left_join(tab[keys], cbind(se[keys], .Index = seq_len(nrow(se))), by = keys)$.Index
    matched <- !is.na(idx)
    if (!isTRUE(all.equal(as.numeric(tab$Estimate[matched]), as.numeric(se$Estimate[idx[matched]]),
                         tolerance = 1e-10, check.attributes = FALSE))) {
      stop("fit/diagnostics mismatch: stored step estimates differ from the fitted thresholds.", call. = FALSE)
    }
    fields <- intersect(c("SE", "CI_Lower", "CI_Upper", "CI_Level", "SE_Method", "SE_Status",
      "SE_Detail", "InferenceReady", "SupportsFormalInference", "CIEligible", "CIUse"), names(se))
    for (field in fields) {
      if (!field %in% names(tab)) tab[[field]] <- se[[field]][rep(NA_integer_, nrow(tab))]
      tab[[field]][matched] <- se[[field]][idx[matched]]
    }
    tab$UncertaintySourceStatus[matched] <- if (source_verified) "source_verified" else "legacy_structural_match"
    if (!source_verified) {
      # Old diagnostics retain their numbers, but cannot certify the entire calibration.
      tab$SourceCIEligible <- tab$CIEligible
      tab$SourceCIUse <- tab$CIUse
      tab$CIEligible[matched] <- FALSE
      tab$CIUse[matched & (is.finite(tab$CI_Lower) | is.finite(tab$CI_Upper))] <- "review_before_reporting"
      if ("SupportsFormalInference" %in% names(tab)) tab$SupportsFormalInference[matched] <- FALSE
    }
  }
  if (identical(as.integer(fit$config$n_cat), 2L) && !inherits(fit, "mfrm_imported_fit")) {
    tab$SE_Status <- "fixed_by_step_centering"
    tab$CIEligible <- FALSE
    tab$CIUse <- "not_applicable"
  }
  list(table = tab, coverage = coverage)
}

category_owner_notes <- function(context) {
  c(paste0("Step owner: ", context$StepOwner, "; selection: ", context$StepSelection, "."),
    .mfrm_step_interpretation_note(),
    "Category counts are observed ratings, not independent Persons. Owner percentages use each owner's total weight; no-exposure percentages are unavailable.",
    "An unobserved category does not establish structural absence. Review the rating plan and rubric before merging categories.",
    "Step SEs describe the stored step offsets, not full ability-axis crossings, adjacent-step gaps or decision cuts. No additional covariance calculation is requested.",
    "Legacy step uncertainty without full source identity is retained for review; inspect SE_Status, CIUse and UncertaintySourceStatus.")
}
