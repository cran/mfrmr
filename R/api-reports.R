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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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

fit_measure_status_label <- function(values, lower, upper, zstd_cut, kind = c("mnsq", "zstd")) {
  kind <- match.arg(kind)
  vals <- suppressWarnings(as.numeric(values))
  out <- rep("not_available", length(vals))
  ok <- is.finite(vals)
  if (identical(kind, "mnsq")) {
    out[ok & vals < lower] <- "overfit"
    out[ok & vals > upper] <- "underfit"
    out[ok & vals >= lower & vals <= upper] <- "within_band"
  } else {
    out[ok & vals <= -abs(zstd_cut)] <- "overfit"
    out[ok & vals >= abs(zstd_cut)] <- "underfit"
    out[ok & abs(vals) < abs(zstd_cut)] <- "within_band"
  }
  out
}

first_existing_fit_measure_column <- function(tbl, columns) {
  hit <- intersect(columns, names(tbl))
  if (length(hit) == 0L) return(rep(NA_real_, nrow(tbl)))
  suppressWarnings(as.numeric(tbl[[hit[1L]]]))
}

fit_measure_reason <- function(infit_band, outfit_band, infit_z_band, outfit_z_band) {
  reasons <- character(0)
  if (identical(infit_band, "underfit")) reasons <- c(reasons, "Infit MnSq high")
  if (identical(outfit_band, "underfit")) reasons <- c(reasons, "Outfit MnSq high")
  if (identical(infit_z_band, "underfit")) reasons <- c(reasons, "Infit ZSTD high")
  if (identical(outfit_z_band, "underfit")) reasons <- c(reasons, "Outfit ZSTD high")
  if (identical(infit_band, "overfit")) reasons <- c(reasons, "Infit MnSq low")
  if (identical(outfit_band, "overfit")) reasons <- c(reasons, "Outfit MnSq low")
  if (identical(infit_z_band, "overfit")) reasons <- c(reasons, "Infit ZSTD low")
  if (identical(outfit_z_band, "overfit")) reasons <- c(reasons, "Outfit ZSTD low")
  if (length(reasons) == 0L) return("")
  paste(reasons, collapse = "; ")
}

make_facets_fit_measure_labels <- function(tbl) {
  keep <- intersect(
    c("Facet", "Level", "Measure", "S.E.", "Lower CI", "Upper CI",
      "CI Level", "Obs", "Infit MnSq", "Infit ZStd", "Outfit MnSq",
      "Outfit ZStd", "Infit df", "Outfit df", "Fit df method",
      "FACETS Infit df", "FACETS Outfit df", "FACETS Infit ZStd",
      "FACETS Outfit ZStd", "Max ZStd shift", "Flag changed by df",
      "Max df rel shift", "df review", "Fit Status", "Review Reason"),
    names(tbl)
  )
  tbl[, keep, drop = FALSE]
}

fit_measure_validate_ci_level <- function(ci_level) {
  ci_level <- suppressWarnings(as.numeric(ci_level[1]))
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  ci_level
}

fit_measure_validate_nonnegative_finite <- function(x, arg) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x) || x < 0) {
    stop("`", arg, "` must be a non-negative finite number.", call. = FALSE)
  }
  x
}

#' Guide FACETS-style fit df and ZSTD standardization
#'
#' @description
#' `facets_fit_df_guide()` gives a compact user-facing guide to the degrees of
#' freedom and ZSTD standardization choices used when comparing mfrmr fit
#' output with FACETS-style fit tables.
#'
#' @param include_references If `TRUE`, include source-reference rows for the
#'   FACETS/Winsteps documentation and Rasch measurement texts that motivate the
#'   guide.
#'
#' @details
#' The guide separates mean-square size from ZSTD standardization. Infit and
#' outfit MnSq values answer how large the residual noise or predictability
#' signal is. ZSTD values standardize those MnSq values using a degrees-of-
#' freedom convention and a Wilson-Hilferty-style transformation, so ZSTD can
#' differ even when the underlying MnSq values are nearly identical.
#'
#' Two boundaries sit upstream of any df comparison. First, the residual
#' basis: `method = "MML"` fits evaluate residuals at shrunken EAP person
#' measures, whereas FACETS evaluates them at JMLE estimates, so MnSq values
#' themselves can differ before any standardization is applied; refit with
#' `method = "JML"` when the comparison requires a JMLE-style residual basis.
#' Second, small df: `mfrmr` returns `NA` ZSTD when `df < 1` because the
#' Wilson-Hilferty transformation is numerically unstable there, while
#' FACETS/Winsteps under `WHEXACT` can continue with a linear approximation,
#' so sparse cells can show `NA` against a finite external value without
#' indicating a fit difference.
#'
#' @return A bundle of class `mfrm_facets_fit_df_guide` with:
#' - `summary`: one-row scope summary
#' - `formula_guide`: formulas and package columns
#' - `column_guide`: where engine and FACETS-style columns appear
#' - `decision_guide`: recommended comparison steps
#' - `interpretation_guide`: how to read common difference patterns
#' - `references`: optional source-reference rows
#' - `settings`: guide metadata
#'
#' @seealso [diagnose_mfrm()], [fit_measures_table()],
#'   [facets_fit_review()]
#' @examples
#' facets_fit_df_guide()
#' facets_fit_df_guide()$decision_guide
#' @export
facets_fit_df_guide <- function(include_references = TRUE) {
  if (!is.logical(include_references) || length(include_references) != 1L ||
      is.na(include_references)) {
    stop("`include_references` must be TRUE or FALSE.", call. = FALSE)
  }

  summary_tbl <- data.frame(
    Scope = "FACETS-style fit df and ZSTD comparison",
    PrimaryRule = "Compare MnSq first; compare ZSTD only after checking df and transformation settings.",
    RecommendedRoute = "diagnose_mfrm(fit_df_method = \"both\") -> facets_fit_review()",
    DefaultMfrmrPrimary = "engine df unless fit_df_method = \"facets\" is requested",
    stringsAsFactors = FALSE
  )

  formula_guide <- data.frame(
    Quantity = c(
      "Engine infit df",
      "Engine outfit df",
      "FACETS-style df",
      "Wilson-Hilferty ZSTD",
      "Linear normal approximation"
    ),
    Formula = c(
      "sum(Var * Weight)",
      "sum(Weight)",
      "2 / q^2, implemented as 2 * numerator^2 / fourth-moment denominator",
      "(MnSq^(1/3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))",
      "(MnSq - 1) * sqrt(df / 2)"
    ),
    MfrmrColumns = c(
      "DF_Infit or DF_Infit_ENGINE",
      "DF_Outfit or DF_Outfit_ENGINE",
      "DF_Infit_FACETS / DF_Outfit_FACETS",
      "InfitZSTD / OutfitZSTD, or *_FACETS companion columns",
      "Used when whexact = TRUE"
    ),
    Use = c(
      "Package-native standardization and ordinary diagnostics.",
      "Package-native standardization and ordinary diagnostics.",
      "FACETS comparison layer; affects ZSTD but not MnSq.",
      "Default ZSTD transformation used for fit-review output.",
      "Optional comparison mode for WHEXACT-style review."
    ),
    stringsAsFactors = FALSE
  )

  column_guide <- data.frame(
    Route = c(
      "diagnose_mfrm(fit_df_method = \"engine\")",
      "diagnose_mfrm(fit_df_method = \"facets\")",
      "diagnose_mfrm(fit_df_method = \"both\")",
      "fit_measures_table(fit_df_method = \"both\")",
      "facets_fit_review()"
    ),
    PrimaryColumns = c(
      "DF_Infit, DF_Outfit, InfitZSTD, OutfitZSTD use engine df.",
      "DF_Infit, DF_Outfit, InfitZSTD, OutfitZSTD use FACETS-style df.",
      "DF_Infit, DF_Outfit, InfitZSTD, OutfitZSTD use engine df.",
      "FACETS-style table keeps primary df and companion FACETS df/ZSTD columns.",
      "Internal comparison places engine and FACETS-style df/ZSTD side by side."
    ),
    CompanionColumns = c(
      "No FACETS companion columns are retained.",
      "Engine companion columns are retained where available.",
      "DF_*_ENGINE, DF_*_FACETS, *ZSTD_ENGINE, and *ZSTD_FACETS are retained.",
      "FACETS Infit df, FACETS Outfit df, FACETS Infit ZStd, FACETS Outfit ZStd.",
      "DFRatio, ZSTDDiff, and FlagChangedByDf columns."
    ),
    UseWhen = c(
      "Routine package-native diagnostics.",
      "You want the primary ZSTD columns to mimic the FACETS-style df convention.",
      "You need to explain why ZSTD flags change under a different df convention.",
      "You want a FACETS-readable table while preserving R-friendly columns.",
      "You are comparing imported FACETS output or preparing a methods note."
    ),
    stringsAsFactors = FALSE
  )

  decision_guide <- data.frame(
    Step = seq_len(5L),
    Question = c(
      "Are MnSq values close?",
      "Are df values close under the same convention?",
      "Do ZSTD values differ after MnSq and df agree?",
      "Does |ZSTD| > 2 status change only after changing df convention?",
      "Is an external FACETS table supplied?"
    ),
    RecommendedAction = c(
      "If MnSq differs materially, treat this as a fit-statistic or estimation difference before discussing ZSTD.",
      "If df differs, classify the ZSTD gap as a df-convention issue unless MnSq also differs.",
      "Check WHEXACT/normalization settings and rounding/truncation before making a substantive claim.",
      "Report the flag as convention-sensitive; inspect MnSq and substantive context before acting on it.",
      "Use read_facets_fit_table() or normalize_facets_fit_frame(), then facets_fit_review()."
    ),
    stringsAsFactors = FALSE
  )

  interpretation_guide <- data.frame(
    Pattern = c(
      "MnSq same, df different, ZSTD different",
      "MnSq different",
      "MnSq different and the mfrmr fit used method = \"MML\"",
      "Small df with counterintuitive ZSTD sign",
      "mfrmr ZSTD is NA but the external table reports a value",
      "FACETS-style flag but engine flag absent",
      "Engine flag but FACETS-style flag absent"
    ),
    Interpretation = c(
      "Usually a standardization-convention difference, not a different residual fit signal.",
      "Potential estimation, weighting, missing-data, or table-matching difference.",
      "MML residuals are evaluated at shrunken EAP person measures, while FACETS uses JMLE estimates; the residual basis itself differs, most visibly for extreme-scoring persons.",
      "Known small-df behavior of Wilson-Hilferty-style standardization; prioritize MnSq and context.",
      "mfrmr withholds ZSTD when df < 1 (Wilson-Hilferty instability); FACETS/Winsteps WHEXACT can continue with a linear approximation on the same cell.",
      "The FACETS-style df makes the same MnSq more statistically extreme.",
      "The engine df makes the same MnSq more statistically extreme."
    ),
    ReportingAction = c(
      "State the df convention and compare MnSq separately from ZSTD.",
      "Do not explain the difference as only a df issue until MnSq matching is resolved.",
      "Refit with method = \"JML\" before attributing the gap to fit computation; report the residual basis in methods notes.",
      "Avoid strong claims from ZSTD alone; show MnSq and df together.",
      "Treat as a small-df availability difference, not a fit difference; compare MnSq for that row instead.",
      "Label as convention-sensitive and review the actual MnSq band.",
      "Label as convention-sensitive and review the actual MnSq band."
    ),
    stringsAsFactors = FALSE
  )

  references <- if (isTRUE(include_references)) {
    data.frame(
      Source = c(
        "Winsteps WHEXACT documentation",
        "FACETS Diagnosing Misfit documentation",
        "Facets Tutorial 2",
        "Wright & Masters (1982)"
      ),
      Supports = c(
        "Wilson-Hilferty/WHEXACT normalization and small-df caveats.",
        "Separating MnSq size from ZSTD significance-style interpretation.",
        "Mean-square as chi-square divided by df; ZStd as Wilson-Hilferty standardization.",
        "Rasch rating-scale fit-statistic interpretation used by FACETS-style reporting."
      ),
      URL = c(
        "https://www.winsteps.com/winman/whexact.htm",
        "https://winsteps.com/facetman/diagnosingmisfit.htm",
        "https://www.winsteps.com/a/ftutorial2.pdf",
        NA_character_
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  out <- list(
    summary = summary_tbl,
    formula_guide = formula_guide,
    column_guide = column_guide,
    decision_guide = decision_guide,
    interpretation_guide = interpretation_guide,
    references = references,
    settings = list(include_references = isTRUE(include_references))
  )
  as_mfrm_bundle(out, "mfrm_facets_fit_df_guide")
}

fit_measure_threshold_profile_table <- function(lower,
                                                upper,
                                                threshold_profiles = c("literature", "active", "all", "none")) {
  threshold_profiles <- match.arg(
    tolower(as.character(threshold_profiles[1])),
    c("literature", "active", "all", "none")
  )
  active <- data.frame(
    Profile = "active",
    ProfileLabel = "Active review band",
    Lower = lower,
    Upper = upper,
    Source = "Current call/options",
    SuggestedUse = "The band used for the main fit-measures table",
    stringsAsFactors = FALSE
  )
  literature <- data.frame(
    Profile = c(
      "linacre_productive",
      "wright_linacre_high_stakes_mcq",
      "wright_linacre_routine_mcq",
      "wright_linacre_rating_scale",
      "wright_linacre_clinical_observation",
      "wright_linacre_judged_performance"
    ),
    ProfileLabel = c(
      "Productive measurement",
      "High-stakes multiple-choice",
      "Routine multiple-choice",
      "Rating-scale surveys",
      "Clinical observation",
      "Judged performance"
    ),
    Lower = c(0.5, 0.8, 0.7, 0.6, 0.5, 0.4),
    Upper = c(1.5, 1.2, 1.3, 1.4, 1.7, 1.2),
    Source = c(
      "Linacre (2002); Bond & Fox (2015)",
      rep("Wright & Linacre (1994)", 5)
    ),
    SuggestedUse = c(
      "Broad screening band for productive measurement",
      "High-stakes selected-response tests",
      "Routine selected-response tests",
      "Rating-scale surveys and questionnaires",
      "Clinical observation ratings",
      "Judged performance ratings"
    ),
    stringsAsFactors = FALSE
  )
  if (identical(threshold_profiles, "none")) {
    return(literature[0, , drop = FALSE])
  }
  if (identical(threshold_profiles, "active")) return(active)
  if (identical(threshold_profiles, "all")) {
    return(rbind(active, literature))
  }
  literature
}

fit_measure_status_for_band <- function(tbl, lower, upper, zstd_cut) {
  infit_band <- fit_measure_status_label(tbl$Infit, lower, upper, zstd_cut, "mnsq")
  outfit_band <- fit_measure_status_label(tbl$Outfit, lower, upper, zstd_cut, "mnsq")
  infit_z_band <- fit_measure_status_label(tbl$InfitZSTD, lower, upper, zstd_cut, "zstd")
  outfit_z_band <- fit_measure_status_label(tbl$OutfitZSTD, lower, upper, zstd_cut, "zstd")
  underfit <- infit_band == "underfit" | outfit_band == "underfit" |
    infit_z_band == "underfit" | outfit_z_band == "underfit"
  overfit <- infit_band == "overfit" | outfit_band == "overfit" |
    infit_z_band == "overfit" | outfit_z_band == "overfit"
  available <- is.finite(tbl$Infit) | is.finite(tbl$Outfit) |
    is.finite(tbl$InfitZSTD) | is.finite(tbl$OutfitZSTD)
  ifelse(
    !available, "not_available",
    ifelse(underfit & overfit, "mixed",
           ifelse(underfit, "underfit",
                  ifelse(overfit, "overfit", "within_band")))
  )
}

fit_measure_profile_counts <- function(tbl, facet_value) {
  available <- tbl$FitStatus != "not_available"
  denom <- sum(available, na.rm = TRUE)
  safe_rate <- function(n) if (denom > 0L) n / denom else NA_real_
  under_n <- sum(tbl$FitStatus == "underfit", na.rm = TRUE)
  over_n <- sum(tbl$FitStatus == "overfit", na.rm = TRUE)
  mixed_n <- sum(tbl$FitStatus == "mixed", na.rm = TRUE)
  within_n <- sum(tbl$FitStatus == "within_band", na.rm = TRUE)
  not_avail_n <- sum(tbl$FitStatus == "not_available", na.rm = TRUE)
  data.frame(
    Facet = facet_value,
    Rows = nrow(tbl),
    AvailableRows = denom,
    UnderfitRows = under_n,
    OverfitRows = over_n,
    MixedRows = mixed_n,
    WithinBandRows = within_n,
    NotAvailableRows = not_avail_n,
    UnderfitRate = safe_rate(under_n),
    OverfitRate = safe_rate(over_n),
    MixedRate = safe_rate(mixed_n),
    AnyFlagRate = safe_rate(under_n + over_n + mixed_n),
    stringsAsFactors = FALSE
  )
}

summarize_fit_measure_profiles <- function(tbl, profile_tbl, zstd_cut) {
  if (nrow(profile_tbl) == 0L || nrow(tbl) == 0L) {
    out <- profile_tbl[0, , drop = FALSE]
    out$ZSTDCut <- numeric(0)
    out$Facet <- character(0)
    out$Rows <- integer(0)
    out$AvailableRows <- integer(0)
    out$UnderfitRows <- integer(0)
    out$OverfitRows <- integer(0)
    out$MixedRows <- integer(0)
    out$WithinBandRows <- integer(0)
    out$NotAvailableRows <- integer(0)
    out$UnderfitRate <- numeric(0)
    out$OverfitRate <- numeric(0)
    out$MixedRate <- numeric(0)
    out$AnyFlagRate <- numeric(0)
    return(out)
  }
  chunks <- vector("list", nrow(profile_tbl))
  for (i in seq_len(nrow(profile_tbl))) {
    profile <- profile_tbl[i, , drop = FALSE]
    tmp <- tbl
    tmp$FitStatus <- fit_measure_status_for_band(
      tmp,
      lower = profile$Lower,
      upper = profile$Upper,
      zstd_cut = zstd_cut
    )
    facets <- sort(unique(as.character(tmp$Facet)))
    facet_counts <- do.call(rbind, lapply(facets, function(facet) {
      fit_measure_profile_counts(tmp[tmp$Facet == facet, , drop = FALSE], facet)
    }))
    overall_counts <- fit_measure_profile_counts(tmp, "All facets")
    counts <- rbind(overall_counts, facet_counts)
    chunks[[i]] <- cbind(
      profile[rep(1L, nrow(counts)), , drop = FALSE],
      ZSTDCut = zstd_cut,
      counts,
      row.names = NULL
    )
  }
  out <- do.call(rbind, chunks)
  row.names(out) <- NULL
  out
}

build_fit_measure_df_sensitivity <- function(tbl,
                                             zstd_cut = 2,
                                             df_zstd_tolerance = 0.05,
                                             df_zstd_large_shift = 0.5,
                                             df_ratio_tolerance = 0.05) {
  if (nrow(tbl) == 0L) {
    return(data.frame())
  }
  required <- c(
    "Facet", "Level", "Infit", "Outfit",
    "DF_Infit_ENGINE", "DF_Outfit_ENGINE",
    "DF_Infit_FACETS", "DF_Outfit_FACETS",
    "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
    "InfitZSTD_FACETS", "OutfitZSTD_FACETS"
  )
  missing <- setdiff(required, names(tbl))
  if (length(missing) > 0L) {
    return(data.frame(
      Facet = as.character(tbl$Facet %||% character(0)),
      Level = as.character(tbl$Level %||% character(0)),
      DfSensitivityStatus = rep("not_available", nrow(tbl)),
      Interpretation = rep("FACETS-style df/ZSTD companion columns are not available.", nrow(tbl)),
      stringsAsFactors = FALSE
    ))
  }

  infit_df_ratio <- suppressWarnings(tbl$DF_Infit_ENGINE / tbl$DF_Infit_FACETS)
  outfit_df_ratio <- suppressWarnings(tbl$DF_Outfit_ENGINE / tbl$DF_Outfit_FACETS)
  infit_df_relative_diff <- ifelse(
    is.finite(tbl$DF_Infit_ENGINE) & is.finite(tbl$DF_Infit_FACETS) & tbl$DF_Infit_FACETS != 0,
    abs(tbl$DF_Infit_ENGINE - tbl$DF_Infit_FACETS) / abs(tbl$DF_Infit_FACETS),
    NA_real_
  )
  outfit_df_relative_diff <- ifelse(
    is.finite(tbl$DF_Outfit_ENGINE) & is.finite(tbl$DF_Outfit_FACETS) & tbl$DF_Outfit_FACETS != 0,
    abs(tbl$DF_Outfit_ENGINE - tbl$DF_Outfit_FACETS) / abs(tbl$DF_Outfit_FACETS),
    NA_real_
  )
  max_df_relative_diff <- fit_review_pmax_na(infit_df_relative_diff, outfit_df_relative_diff)
  infit_z_diff <- suppressWarnings(tbl$InfitZSTD_FACETS - tbl$InfitZSTD_ENGINE)
  outfit_z_diff <- suppressWarnings(tbl$OutfitZSTD_FACETS - tbl$OutfitZSTD_ENGINE)
  max_abs_z_diff <- fit_review_pmax_na(abs(infit_z_diff), abs(outfit_z_diff))
  max_abs_log_df_ratio <- fit_review_pmax_na(abs(log(infit_df_ratio)), abs(log(outfit_df_ratio)))
  engine_flag <- abs(tbl$InfitZSTD_ENGINE) >= zstd_cut | abs(tbl$OutfitZSTD_ENGINE) >= zstd_cut
  facets_flag <- abs(tbl$InfitZSTD_FACETS) >= zstd_cut | abs(tbl$OutfitZSTD_FACETS) >= zstd_cut
  engine_flag[is.na(engine_flag)] <- FALSE
  facets_flag[is.na(facets_flag)] <- FALSE
  flag_changed <- engine_flag != facets_flag

  available <- is.finite(max_abs_z_diff) | is.finite(max_abs_log_df_ratio) | is.finite(max_df_relative_diff)
  status <- ifelse(
    !available, "not_available",
    ifelse(flag_changed, "flag_changed_by_df",
           ifelse(is.finite(max_abs_z_diff) & max_abs_z_diff >= df_zstd_large_shift,
                  "large_zstd_shift",
                  ifelse(is.finite(max_df_relative_diff) & max_df_relative_diff > df_ratio_tolerance,
                         "df_convention_difference",
                         ifelse(is.finite(max_abs_z_diff) & max_abs_z_diff > df_zstd_tolerance,
                                "small_zstd_shift",
                                "same_or_rounding"))))
  )
  interpretation <- dplyr::case_when(
    status == "flag_changed_by_df" ~ "The same MnSq values cross the ZSTD flag threshold under one df convention but not the other.",
    status == "large_zstd_shift" ~ "The FACETS-style df changes ZSTD substantially; interpret ZSTD only with the df convention stated.",
    status == "df_convention_difference" ~ "The df convention differs enough to affect ZSTD interpretation even if the flag status is unchanged.",
    status == "small_zstd_shift" ~ "ZSTD differs slightly after changing df convention; usually a convention/rounding note.",
    status == "same_or_rounding" ~ "Engine and FACETS-style standardization are practically the same for this row.",
    TRUE ~ "FACETS-style df/ZSTD companion columns are not available."
  )

  out <- data.frame(
    Facet = as.character(tbl$Facet),
    Level = as.character(tbl$Level),
    Infit = suppressWarnings(as.numeric(tbl$Infit)),
    Outfit = suppressWarnings(as.numeric(tbl$Outfit)),
    DF_Infit_ENGINE = suppressWarnings(as.numeric(tbl$DF_Infit_ENGINE)),
    DF_Infit_FACETS = suppressWarnings(as.numeric(tbl$DF_Infit_FACETS)),
    DF_Outfit_ENGINE = suppressWarnings(as.numeric(tbl$DF_Outfit_ENGINE)),
    DF_Outfit_FACETS = suppressWarnings(as.numeric(tbl$DF_Outfit_FACETS)),
    InfitZSTD_ENGINE = suppressWarnings(as.numeric(tbl$InfitZSTD_ENGINE)),
    InfitZSTD_FACETS = suppressWarnings(as.numeric(tbl$InfitZSTD_FACETS)),
    OutfitZSTD_ENGINE = suppressWarnings(as.numeric(tbl$OutfitZSTD_ENGINE)),
    OutfitZSTD_FACETS = suppressWarnings(as.numeric(tbl$OutfitZSTD_FACETS)),
    InfitDFRatio_ENGINE_over_FACETS = infit_df_ratio,
    OutfitDFRatio_ENGINE_over_FACETS = outfit_df_ratio,
    InfitDFRelativeDifference_ENGINE_vs_FACETS = infit_df_relative_diff,
    OutfitDFRelativeDifference_ENGINE_vs_FACETS = outfit_df_relative_diff,
    InfitZSTDDiff_FACETS_minus_ENGINE = infit_z_diff,
    OutfitZSTDDiff_FACETS_minus_ENGINE = outfit_z_diff,
    MaxAbsZSTDDiff_FACETS_vs_ENGINE = max_abs_z_diff,
    MaxAbsLogDFRatio_ENGINE_over_FACETS = max_abs_log_df_ratio,
    MaxDFRelativeDifference_ENGINE_vs_FACETS = max_df_relative_diff,
    EngineFlagAbsZ = engine_flag,
    FacetsStyleFlagAbsZ = facets_flag,
    FlagChangedByDf = flag_changed,
    DfSensitivityStatus = status,
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  )
  out |>
    dplyr::arrange(
      dplyr::desc(.data$FlagChangedByDf),
      dplyr::desc(.data$MaxAbsZSTDDiff_FACETS_vs_ENGINE),
      .data$Facet,
      .data$Level
    ) |>
    as.data.frame(stringsAsFactors = FALSE)
}

summarize_fit_measure_df_sensitivity <- function(df_sensitivity) {
  if (nrow(df_sensitivity) == 0L) {
    return(data.frame(
      Rows = 0L,
      ComparedRows = 0L,
      FlagChangedByDfRows = 0L,
      LargeZSTDShiftRows = 0L,
      DfConventionDifferenceRows = 0L,
      SmallZSTDShiftRows = 0L,
      SameOrRoundingRows = 0L,
      stringsAsFactors = FALSE
    ))
  }
  status <- as.character(df_sensitivity$DfSensitivityStatus %||% "not_available")
  data.frame(
    Rows = nrow(df_sensitivity),
    ComparedRows = sum(status != "not_available", na.rm = TRUE),
    FlagChangedByDfRows = sum(status == "flag_changed_by_df", na.rm = TRUE),
    LargeZSTDShiftRows = sum(status == "large_zstd_shift", na.rm = TRUE),
    DfConventionDifferenceRows = sum(status == "df_convention_difference", na.rm = TRUE),
    SmallZSTDShiftRows = sum(status == "small_zstd_shift", na.rm = TRUE),
    SameOrRoundingRows = sum(status == "same_or_rounding", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Build a FACETS-style fit-measures review table
#'
#' @param x Output from [fit_mfrm()] or [diagnose_mfrm()].
#' @param diagnostics Optional diagnostics object. If supplied, `x` may be the
#'   fitted object used only for provenance.
#' @param facet Optional facet-name filter, for example `"Rater"`.
#' @param include_person Logical; if `FALSE` (default), excludes the `Person`
#'   facet so operational facet elements are shown first.
#' @param lower,upper Optional mean-square review band. Defaults to
#'   [mfrm_misfit_thresholds()].
#' @param zstd_cut Absolute ZSTD cutoff used for directional underfit/overfit
#'   flags. Default `2`.
#' @param ci_level Confidence level used to add approximate Wald intervals for
#'   facet measures. Default `0.95`.
#' @param threshold_profiles Which mean-square threshold profiles to summarize
#'   in addition to the active table band. `"literature"` (default) returns
#'   commonly cited bands from Linacre, Bond & Fox, and Wright & Linacre;
#'   `"active"` returns only the active band; `"all"` returns both; `"none"`
#'   suppresses profile summaries.
#' @param fit_df_method Degrees-of-freedom convention used when `diagnostics`
#'   is computed inside the helper. `"engine"` keeps the package-native fit df,
#'   `"facets"` makes primary ZSTD columns use the FACETS/Wright-Masters
#'   fourth-moment df convention, and `"both"` keeps engine columns primary
#'   while adding FACETS-style companion df/ZSTD columns for comparison.
#' @param df_zstd_tolerance Smallest absolute engine-vs-FACETS-style ZSTD
#'   difference treated as interpretively visible rather than rounding noise
#'   in `df_sensitivity`. Default `0.05`.
#' @param df_zstd_large_shift Absolute engine-vs-FACETS-style ZSTD difference
#'   labeled `large_zstd_shift` when the `zstd_cut` flag status is unchanged.
#'   Default `0.5`.
#' @param df_ratio_tolerance Relative df-difference threshold used to label
#'   `df_convention_difference`; for example, `0.05` means a 5 percent
#'   engine-vs-FACETS-style df difference. Default `0.05`.
#' @param sort_by Sorting rule: `"status"` prioritizes underfit/overfit rows,
#'   `"abs_zstd"` sorts by largest absolute ZSTD, and `"facet"` / `"level"`
#'   sort alphabetically.
#' @param top_n Optional maximum number of rows in the returned main table.
#'
#' @details
#' This helper gives users a direct table route for the common FACETS-style
#' question: which raters, criteria, or other facet elements show underfit or
#' overfit? It uses the fit statistics already computed by [diagnose_mfrm()].
#'
#' Directional labels are based on both mean-square and ZSTD evidence:
#' high MnSq or positive large ZSTD is labeled `underfit`; low MnSq or negative
#' large ZSTD is labeled `overfit`. Rows with conflicting directions are labeled
#' `mixed`. Treat the table as a review screen and inspect substantive context
#' before removing raters or changing an instrument.
#'
#' FACETS-style ZSTD comparison is controlled by `fit_df_method`. MnSq values
#' should be compared first; df and ZSTD columns explain how the same MnSq values
#' are standardized. Use `fit_df_method = "both"` when preparing a table for
#' FACETS users or when explaining why |ZSTD| flags change across df
#' conventions. The `df_zstd_tolerance`, `df_zstd_large_shift`, and
#' `df_ratio_tolerance` arguments make the df-sensitivity screen explicit so
#' the same table can be reproduced under stricter or more permissive review
#' rules.
#'
#' @return A bundle of class `mfrm_fit_measures` with:
#' - `table`: R-friendly fit-measure table with status columns
#' - `facets_table`: FACETS-style column labels for reporting/review
#' - `status_summary`: counts by facet and fit status
#' - `profile_summary_by_facet`: underfit/overfit rates for each threshold
#'   profile and facet
#' - `profile_summary_overall`: threshold-profile rates pooled over facets
#' - `df_sensitivity`: row-level engine-vs-FACETS-style df/ZSTD comparison
#' - `df_sensitive`: subset of rows where df convention changes the ZSTD flag
#'   or materially changes ZSTD interpretation
#' - `df_sensitivity_summary`: counts of df-sensitive rows
#' - `underfit`, `overfit`, `mixed`: filtered row subsets
#' - `df_conversion_guide`: FACETS-style df/ZSTD comparison guide
#' - `settings`: thresholds and filters used
#'
#' @seealso [diagnose_mfrm()], [facets_fit_review()], [plot_bubble()],
#'   [mfrm_misfit_thresholds()]
#' @concept confidence intervals
#' @concept fit statistics
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' fm <- fit_measures_table(fit, facet = "Rater")
#' fm$facets_table
#' fm$underfit
#'
#' # Include FACETS-style df/ZSTD companion columns for comparison.
#' fm_facets <- fit_measures_table(fit, facet = "Rater", fit_df_method = "both")
#' fm_facets$df_conversion_guide$decision_guide
#' }
#' @export
fit_measures_table <- function(x,
                               diagnostics = NULL,
                               facet = NULL,
                               include_person = FALSE,
                               lower = NULL,
                               upper = NULL,
                               zstd_cut = 2,
                               ci_level = 0.95,
                               threshold_profiles = c("literature", "active", "all", "none"),
                               fit_df_method = c("engine", "facets", "both"),
                               df_zstd_tolerance = 0.05,
                               df_zstd_large_shift = 0.5,
                               df_ratio_tolerance = 0.05,
                               sort_by = c("status", "abs_zstd", "facet", "level"),
                               top_n = Inf) {
  sort_by <- match.arg(tolower(as.character(sort_by[1])), c("status", "abs_zstd", "facet", "level"))
  threshold_profiles <- match.arg(
    tolower(as.character(threshold_profiles[1])),
    c("literature", "active", "all", "none")
  )
  fit_df_method <- match_fit_df_method(fit_df_method)
  if (!is.logical(include_person) || length(include_person) != 1L || is.na(include_person)) {
    stop("`include_person` must be TRUE or FALSE.", call. = FALSE)
  }
  zstd_cut <- suppressWarnings(as.numeric(zstd_cut[1]))
  if (!is.finite(zstd_cut) || zstd_cut <= 0) {
    stop("`zstd_cut` must be a positive finite number.", call. = FALSE)
  }
  ci_level <- fit_measure_validate_ci_level(ci_level)
  df_zstd_tolerance <- fit_measure_validate_nonnegative_finite(df_zstd_tolerance, "df_zstd_tolerance")
  df_zstd_large_shift <- fit_measure_validate_nonnegative_finite(df_zstd_large_shift, "df_zstd_large_shift")
  df_ratio_tolerance <- fit_measure_validate_nonnegative_finite(df_ratio_tolerance, "df_ratio_tolerance")
  if (df_zstd_large_shift < df_zstd_tolerance) {
    stop("`df_zstd_large_shift` must be greater than or equal to `df_zstd_tolerance`.", call. = FALSE)
  }

  diagnostics_supplied <- !is.null(diagnostics)
  if (is.null(diagnostics)) {
    diagnostics <- if (inherits(x, "mfrm_fit")) {
      diagnose_mfrm(x, residual_pca = "none", fit_df_method = fit_df_method)
    } else if (is.list(x) && !is.null(x$measures)) {
      x
    } else {
      stop("`x` must be output from fit_mfrm() or diagnose_mfrm().", call. = FALSE)
    }
  }
  if (!is.list(diagnostics) || is.null(diagnostics$measures)) {
    stop("`diagnostics` must be output from diagnose_mfrm() with a `measures` table.", call. = FALSE)
  }
  if (isTRUE(diagnostics_supplied) && inherits(x, "mfrm_fit") &&
      !identical(fit_df_method, "engine")) {
    measure_names <- names(as.data.frame(diagnostics$measures, stringsAsFactors = FALSE))
    needs_facets_df <- !all(c(
      "DF_Infit_FACETS", "DF_Outfit_FACETS",
      "InfitZSTD_FACETS", "OutfitZSTD_FACETS"
    ) %in% measure_names)
    if (needs_facets_df) {
      diagnostics <- diagnose_mfrm(x, residual_pca = "none", fit_df_method = fit_df_method)
    }
  }

  band <- mfrm_misfit_thresholds(lower = lower, upper = upper)
  lower <- as.numeric(band["lower"])
  upper <- as.numeric(band["upper"])

  measures <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  if (nrow(measures) == 0L || !all(c("Facet", "Level") %in% names(measures))) {
    stop("diagnostics$measures must contain `Facet` and `Level` rows.", call. = FALSE)
  }
  if (!isTRUE(include_person)) {
    measures <- measures[as.character(measures$Facet) != "Person", , drop = FALSE]
  }
  if (!is.null(facet)) {
    facet <- unique(as.character(facet))
    measures <- measures[as.character(measures$Facet) %in% facet, , drop = FALSE]
    if (nrow(measures) == 0L) {
      stop("No fit-measure rows matched `facet`.", call. = FALSE)
    }
  }

  n_obs <- first_existing_fit_measure_column(measures, c("N", "Obs", "Count", "N.x", "N.y"))
  estimate <- first_existing_fit_measure_column(measures, c("Estimate", "Measure"))
  se <- first_existing_fit_measure_column(measures, c("SE", "ModelSE", "S.E."))
  infit <- first_existing_fit_measure_column(measures, "Infit")
  outfit <- first_existing_fit_measure_column(measures, "Outfit")
  infit_z <- first_existing_fit_measure_column(measures, c("InfitZSTD", "InfitZStd", "Infit ZStd"))
  outfit_z <- first_existing_fit_measure_column(measures, c("OutfitZSTD", "OutfitZStd", "Outfit ZStd"))
  df_infit <- first_existing_fit_measure_column(measures, "DF_Infit")
  df_outfit <- first_existing_fit_measure_column(measures, "DF_Outfit")
  df_infit_engine <- first_existing_fit_measure_column(measures, "DF_Infit_ENGINE")
  df_outfit_engine <- first_existing_fit_measure_column(measures, "DF_Outfit_ENGINE")
  df_infit_facets <- first_existing_fit_measure_column(measures, "DF_Infit_FACETS")
  df_outfit_facets <- first_existing_fit_measure_column(measures, "DF_Outfit_FACETS")
  infit_z_engine <- first_existing_fit_measure_column(measures, "InfitZSTD_ENGINE")
  outfit_z_engine <- first_existing_fit_measure_column(measures, "OutfitZSTD_ENGINE")
  infit_z_facets <- first_existing_fit_measure_column(measures, "InfitZSTD_FACETS")
  outfit_z_facets <- first_existing_fit_measure_column(measures, "OutfitZSTD_FACETS")
  fit_df_method_col <- if ("FitDfMethod" %in% names(measures)) {
    as.character(measures$FitDfMethod)
  } else {
    std <- as.data.frame(diagnostics$fit_standardization %||% data.frame(), stringsAsFactors = FALSE)
    primary <- if (nrow(std) > 0 && "PrimaryFitDfMethod" %in% names(std)) {
      as.character(std$PrimaryFitDfMethod[1])
    } else {
      fit_df_method
    }
    rep(primary, nrow(measures))
  }
  fit_zstd_transform <- if ("FitZSTDTransform" %in% names(measures)) {
    as.character(measures$FitZSTDTransform)
  } else {
    std <- as.data.frame(diagnostics$fit_standardization %||% data.frame(), stringsAsFactors = FALSE)
    transform <- if (nrow(std) > 0 && "ZSTDTransform" %in% names(std)) {
      as.character(std$ZSTDTransform[1])
    } else {
      NA_character_
    }
    rep(transform, nrow(measures))
  }
  if (!any(is.finite(df_infit_engine)) && !identical(fit_df_method, "facets")) {
    df_infit_engine <- df_infit
  }
  if (!any(is.finite(df_outfit_engine)) && !identical(fit_df_method, "facets")) {
    df_outfit_engine <- df_outfit
  }
  if (!any(is.finite(infit_z_engine)) && !identical(fit_df_method, "facets")) {
    infit_z_engine <- infit_z
  }
  if (!any(is.finite(outfit_z_engine)) && !identical(fit_df_method, "facets")) {
    outfit_z_engine <- outfit_z
  }
  if (!any(is.finite(df_infit_facets)) && identical(fit_df_method, "facets")) {
    df_infit_facets <- df_infit
  }
  if (!any(is.finite(df_outfit_facets)) && identical(fit_df_method, "facets")) {
    df_outfit_facets <- df_outfit
  }
  if (!any(is.finite(infit_z_facets)) && identical(fit_df_method, "facets")) {
    infit_z_facets <- infit_z
  }
  if (!any(is.finite(outfit_z_facets)) && identical(fit_df_method, "facets")) {
    outfit_z_facets <- outfit_z
  }

  infit_band <- fit_measure_status_label(infit, lower, upper, zstd_cut, "mnsq")
  outfit_band <- fit_measure_status_label(outfit, lower, upper, zstd_cut, "mnsq")
  infit_z_band <- fit_measure_status_label(infit_z, lower, upper, zstd_cut, "zstd")
  outfit_z_band <- fit_measure_status_label(outfit_z, lower, upper, zstd_cut, "zstd")
  underfit <- infit_band == "underfit" | outfit_band == "underfit" |
    infit_z_band == "underfit" | outfit_z_band == "underfit"
  overfit <- infit_band == "overfit" | outfit_band == "overfit" |
    infit_z_band == "overfit" | outfit_z_band == "overfit"
  available <- is.finite(infit) | is.finite(outfit) | is.finite(infit_z) | is.finite(outfit_z)
  status <- ifelse(
    !available, "not_available",
    ifelse(underfit & overfit, "mixed",
           ifelse(underfit, "underfit",
                  ifelse(overfit, "overfit", "within_band")))
  )
  max_abs_z <- apply(cbind(abs(infit_z), abs(outfit_z)), 1L, function(v) {
    if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
  })
  max_mnsq_distance <- apply(cbind(abs(infit - 1), abs(outfit - 1)), 1L, function(v) {
    if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
  })
  reasons <- mapply(
    fit_measure_reason,
    infit_band,
    outfit_band,
    infit_z_band,
    outfit_z_band,
    USE.NAMES = FALSE
  )
  reason_out <- ifelse(nzchar(reasons), reasons, "Within selected review band")
  reason_out[!available] <- "Fit statistics unavailable"
  z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
  ci_ok <- is.finite(estimate) & is.finite(se) & se >= 0
  ci_lower <- ifelse(ci_ok, estimate - z_ci * se, NA_real_)
  ci_upper <- ifelse(ci_ok, estimate + z_ci * se, NA_real_)

  out <- data.frame(
    Facet = as.character(measures$Facet),
    Level = as.character(measures$Level),
    Measure = estimate,
    SE = se,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    CI_Level = ci_level,
    N = n_obs,
    Infit = infit,
    Outfit = outfit,
    InfitZSTD = infit_z,
    OutfitZSTD = outfit_z,
    DF_Infit = df_infit,
    DF_Outfit = df_outfit,
    DF_Infit_ENGINE = df_infit_engine,
    DF_Outfit_ENGINE = df_outfit_engine,
    DF_Infit_FACETS = df_infit_facets,
    DF_Outfit_FACETS = df_outfit_facets,
    InfitZSTD_ENGINE = infit_z_engine,
    OutfitZSTD_ENGINE = outfit_z_engine,
    InfitZSTD_FACETS = infit_z_facets,
    OutfitZSTD_FACETS = outfit_z_facets,
    FitDfMethod = fit_df_method_col,
    FitZSTDTransform = fit_zstd_transform,
    InfitBand = infit_band,
    OutfitBand = outfit_band,
    InfitZSTDBand = infit_z_band,
    OutfitZSTDBand = outfit_z_band,
    Underfit = underfit,
    Overfit = overfit,
    FitStatus = status,
    ReviewReason = reason_out,
    MaxAbsZSTD = max_abs_z,
    MaxMnSqDistance = max_mnsq_distance,
    stringsAsFactors = FALSE
  )
  df_sensitivity_all <- build_fit_measure_df_sensitivity(
    out,
    zstd_cut = zstd_cut,
    df_zstd_tolerance = df_zstd_tolerance,
    df_zstd_large_shift = df_zstd_large_shift,
    df_ratio_tolerance = df_ratio_tolerance
  )
  if (nrow(df_sensitivity_all) > 0L) {
    out_key <- paste(out$Facet, out$Level, sep = "\r")
    sens_key <- paste(df_sensitivity_all$Facet, df_sensitivity_all$Level, sep = "\r")
    sens_idx <- match(out_key, sens_key)
    for (nm in c(
      "InfitZSTDDiff_FACETS_minus_ENGINE",
      "OutfitZSTDDiff_FACETS_minus_ENGINE",
      "MaxAbsZSTDDiff_FACETS_vs_ENGINE",
      "MaxAbsLogDFRatio_ENGINE_over_FACETS",
      "MaxDFRelativeDifference_ENGINE_vs_FACETS",
      "EngineFlagAbsZ",
      "FacetsStyleFlagAbsZ",
      "FlagChangedByDf",
      "DfSensitivityStatus"
    )) {
      out[[nm]] <- df_sensitivity_all[[nm]][sens_idx]
    }
  } else {
    out$InfitZSTDDiff_FACETS_minus_ENGINE <- NA_real_
    out$OutfitZSTDDiff_FACETS_minus_ENGINE <- NA_real_
    out$MaxAbsZSTDDiff_FACETS_vs_ENGINE <- NA_real_
    out$MaxAbsLogDFRatio_ENGINE_over_FACETS <- NA_real_
    out$MaxDFRelativeDifference_ENGINE_vs_FACETS <- NA_real_
    out$EngineFlagAbsZ <- NA
    out$FacetsStyleFlagAbsZ <- NA
    out$FlagChangedByDf <- NA
    out$DfSensitivityStatus <- "not_available"
  }
  status_rank <- c(mixed = 5, underfit = 4, overfit = 3, within_band = 2, not_available = 1)
  out$FitStatusRank <- unname(status_rank[out$FitStatus])
  out$FitStatusRank[is.na(out$FitStatusRank)] <- 0
  out_full <- out
  df_sensitivity_summary <- summarize_fit_measure_df_sensitivity(df_sensitivity_all)
  df_sensitive <- df_sensitivity_all[
    !as.character(df_sensitivity_all$DfSensitivityStatus %||% "not_available") %in%
      c("same_or_rounding", "not_available"),
    ,
    drop = FALSE
  ]

  ord <- switch(
    sort_by,
    status = order(-out$FitStatusRank, -ifelse(is.finite(out$MaxAbsZSTD), out$MaxAbsZSTD, -Inf),
                   -ifelse(is.finite(out$MaxMnSqDistance), out$MaxMnSqDistance, -Inf),
                   out$Facet, out$Level),
    abs_zstd = order(-ifelse(is.finite(out$MaxAbsZSTD), out$MaxAbsZSTD, -Inf), out$Facet, out$Level),
    facet = order(out$Facet, out$Level),
    level = order(out$Level, out$Facet)
  )
  out_sorted <- out[ord, , drop = FALSE]
  out_display <- out_sorted

  if (!is.null(top_n)) {
    top_n_num <- suppressWarnings(as.numeric(top_n[1]))
    if (is.finite(top_n_num)) {
      out_display <- utils::head(out_display, n = max(1L, as.integer(top_n_num)))
    }
  }

  facets_table <- data.frame(
    Facet = out_display$Facet,
    Level = out_display$Level,
    Measure = out_display$Measure,
    `S.E.` = out_display$SE,
    `Lower CI` = out_display$CI_Lower,
    `Upper CI` = out_display$CI_Upper,
    `CI Level` = out_display$CI_Level,
    Obs = out_display$N,
    `Infit MnSq` = out_display$Infit,
    `Infit ZStd` = out_display$InfitZSTD,
    `Outfit MnSq` = out_display$Outfit,
    `Outfit ZStd` = out_display$OutfitZSTD,
    `Infit df` = out_display$DF_Infit,
    `Outfit df` = out_display$DF_Outfit,
    `Fit df method` = out_display$FitDfMethod,
    `Max ZStd shift` = out_display$MaxAbsZSTDDiff_FACETS_vs_ENGINE,
    `Flag changed by df` = out_display$FlagChangedByDf,
    `Max df rel shift` = out_display$MaxDFRelativeDifference_ENGINE_vs_FACETS,
    `df review` = out_display$DfSensitivityStatus,
    `Fit Status` = out_display$FitStatus,
    `Review Reason` = out_display$ReviewReason,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (any(is.finite(out_display$DF_Infit_FACETS)) ||
      any(is.finite(out_display$DF_Outfit_FACETS)) ||
      any(is.finite(out_display$InfitZSTD_FACETS)) ||
      any(is.finite(out_display$OutfitZSTD_FACETS))) {
    facets_table$`FACETS Infit df` <- out_display$DF_Infit_FACETS
    facets_table$`FACETS Outfit df` <- out_display$DF_Outfit_FACETS
    facets_table$`FACETS Infit ZStd` <- out_display$InfitZSTD_FACETS
    facets_table$`FACETS Outfit ZStd` <- out_display$OutfitZSTD_FACETS
  }
  facets_table <- make_facets_fit_measure_labels(facets_table)
  status_summary <- out_full |>
    dplyr::count(.data$Facet, .data$FitStatus, name = "Rows") |>
    dplyr::arrange(.data$Facet, dplyr::desc(.data$Rows), .data$FitStatus) |>
    as.data.frame(stringsAsFactors = FALSE)
  overall_summary <- data.frame(
    Rows = nrow(out_full),
    DisplayedRows = nrow(out_display),
    UnderfitRows = sum(out_full$FitStatus == "underfit", na.rm = TRUE),
    OverfitRows = sum(out_full$FitStatus == "overfit", na.rm = TRUE),
    MixedRows = sum(out_full$FitStatus == "mixed", na.rm = TRUE),
    WithinBandRows = sum(out_full$FitStatus == "within_band", na.rm = TRUE),
    NotAvailableRows = sum(out_full$FitStatus == "not_available", na.rm = TRUE),
    DfComparedRows = df_sensitivity_summary$ComparedRows[1],
    DfSensitiveRows = nrow(df_sensitive),
    FlagChangedByDfRows = df_sensitivity_summary$FlagChangedByDfRows[1],
    LargeZSTDShiftRows = df_sensitivity_summary$LargeZSTDShiftRows[1],
    DfConventionDifferenceRows = df_sensitivity_summary$DfConventionDifferenceRows[1],
    stringsAsFactors = FALSE
  )
  profile_tbl <- fit_measure_threshold_profile_table(
    lower = lower,
    upper = upper,
    threshold_profiles = threshold_profiles
  )
  profile_summary <- summarize_fit_measure_profiles(
    out_full[, setdiff(names(out_full), "FitStatusRank"), drop = FALSE],
    profile_tbl,
    zstd_cut = zstd_cut
  )
  profile_summary_overall <- profile_summary[
    profile_summary$Facet == "All facets",
    ,
    drop = FALSE
  ]
  profile_summary_by_facet <- profile_summary[
    profile_summary$Facet != "All facets",
    ,
    drop = FALSE
  ]

  bundle <- list(
    table = out_display[, setdiff(names(out_display), "FitStatusRank"), drop = FALSE],
    facets_table = facets_table,
    status_summary = status_summary,
    summary = overall_summary,
    threshold_profiles = profile_tbl,
    profile_summary = profile_summary,
    profile_summary_by_facet = profile_summary_by_facet,
    profile_summary_overall = profile_summary_overall,
    df_sensitivity = df_sensitivity_all,
    df_sensitive = df_sensitive,
    df_sensitivity_summary = df_sensitivity_summary,
    underfit = out_sorted[out_sorted$FitStatus == "underfit", setdiff(names(out_sorted), "FitStatusRank"), drop = FALSE],
    overfit = out_sorted[out_sorted$FitStatus == "overfit", setdiff(names(out_sorted), "FitStatusRank"), drop = FALSE],
    mixed = out_sorted[out_sorted$FitStatus == "mixed", setdiff(names(out_sorted), "FitStatusRank"), drop = FALSE],
    df_conversion_guide = facets_fit_df_guide(include_references = TRUE),
    settings = list(
      facet = facet %||% NA_character_,
      include_person = isTRUE(include_person),
      lower = lower,
      upper = upper,
      zstd_cut = zstd_cut,
      ci_level = ci_level,
      df_zstd_tolerance = df_zstd_tolerance,
      df_zstd_large_shift = df_zstd_large_shift,
      df_ratio_tolerance = df_ratio_tolerance,
      threshold_profiles = threshold_profiles,
      fit_df_method = fit_df_method,
      sort_by = sort_by,
      top_n = top_n
    )
  )
  as_mfrm_bundle(bundle, "mfrm_fit_measures")
}

#' Build a data quality summary report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param data Optional raw data frame used for row-level review.
#' @param person Optional person column name in `data`.
#' @param facets Optional facet column names in `data`.
#' @param score Optional score column name in `data`.
#' @param weight Optional weight column name in `data`.
#' @param min_category_count Minimum raw or weighted count used to label a
#'   non-zero facet-level score category as sparse. Default `10`.
#' @param dominant_category_cutoff Proportion in `(0, 1]` used to flag a
#'   facet level whose responses are dominated by one score category. Default
#'   `0.95`.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_data_quality` (`type = "dashboard"`, `"quality_flags"`,
#' `"row_review"`, `"category_counts"`, `"score_support"`,
#' `"facet_category_usage"`, `"facet_response_patterns"`, `"score_map"`,
#' `"missing_rows"`).
#'
#' @section Interpreting output:
#' - `summary`: retained/dropped row overview.
#' - `quality_overview`: area-level QC status for rows, score support,
#'   facet-category use, and design matching.
#' - `quality_flags`: prioritized QC flags with counts and recommended next
#'   actions. This is not an item/person/rater table.
#' - `row_review`: reason-level breakdown for data issues.
#' - `category_counts`: post-filter category usage, including retained
#'   zero-count score-support categories.
#' - `score_support_review`: quick view of zero-count boundary/intermediate
#'   categories and their threshold-functioning caveats.
#' - `category_usage_by_facet`: facet-level category counts over the retained
#'   score support.
#' - `category_usage_summary`: per-facet-level zero/sparse category summary.
#' - `facet_response_patterns`: facet-level response-pattern summaries,
#'   including single-category and dominant-category use.
#' - `caveats`: user-facing score-support warnings, including cases where
#'   non-consecutive original labels such as `1, 2, 4, 5` were recoded because
#'   `keep_original = FALSE`.
#' - `score_map`: original-to-internal score mapping used when labels are
#'   recoded.
#' - `unknown_elements`: facet levels in raw data but not in fitted design.
#'
#' @section Typical workflow:
#' 1. Run `data_quality_report(...)` with raw data.
#' 2. Check `summary(out)` and `plot(out, type = "dashboard")`, then inspect
#'    `quality_flags`, score-support, score-map, facet-response-pattern, and missing/unknown element
#'    sections as needed.
#' 3. Resolve missing values, score-support gaps, and sparse categories before
#'    final estimation/reporting.
#' @return A named list with data-quality report components. Class:
#'   `mfrm_data_quality`.
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [specifications_report()],
#'   [mfrmr_reports_and_tables], [mfrmr_compatibility_layer]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
                                min_category_count = 10,
                                dominant_category_cutoff = 0.95,
                                include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table2_data_summary(
      fit = fit,
      data = data,
      person = person,
      facets = facets,
      score = score,
      weight = weight,
      min_category_count = min_category_count,
      dominant_category_cutoff = dominant_category_cutoff,
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
#' `"design_matrix"` / `"network"`). The network route returns reusable node
#' and edge tables with `draw = FALSE`; drawing uses `igraph` when available.
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
#' @seealso [diagnose_mfrm()], [mfrm_network_analysis()],
#'   [measurable_summary_table()], [data_quality_report()], [mfrmr_linking_and_dff],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' out <- subset_connectivity_report(fit)
#' summary(out)
#' p_sub <- plot(out, draw = FALSE)
#' p_design <- plot(out, type = "design_matrix", draw = FALSE)
#' p_net <- plot(out, type = "network", draw = FALSE)
#' p_sub$data$plot
#' p_design$data$plot
#' p_net$data$edges
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

#' Analyze the MFRM design network
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param top_n_subsets Optional maximum number of connected-subset rows to
#'   retain before constructing the graph.
#' @param min_observations Minimum observations required to keep a subset row.
#' @param include_graph Logical; if `TRUE`, include the underlying `igraph`
#'   object in the returned bundle. Defaults to `FALSE` so outputs remain easy
#'   to serialize.
#'
#' @details
#' `mfrm_network_analysis()` treats the person/facet-level observation design as
#' an undirected weighted graph. Nodes are person or facet levels; edges connect
#' levels that co-occur in at least one observed rating; edge weights are
#' co-observation counts. The resulting network metrics are design diagnostics,
#' not psychometric measures of person ability or rater quality.
#' `plot(net, type = "centrality")`, `plot(net, type = "facet_summary")`, and
#' `plot(net, type = "network")` provide immediate visual checks; use
#' `draw = FALSE` to extract reusable plot data.
#'
#' The most useful review columns are:
#' - `Components`: more than one component means the design has disconnected
#'   measurement subsets.
#' - `IsArticulationPoint`: a node whose removal would increase disconnectedness.
#' - `IsBridge`: an edge whose removal would increase disconnectedness.
#' - `Betweenness`: a routing-dependence indicator; high values identify levels
#'   that carry many shortest paths through the design graph.
#'
#' In incomplete rater-mediated designs, these graph summaries help identify
#' fragile linking structures before interpreting facet measures or planning
#' additional data collection.
#'
#' @section References:
#' - McEwen, M. R. (2015). *Development of a Software Prototype for Generating
#'   and Classifying Incomplete Many-Facet-Rasch Model Rating Designs*.
#'   Brigham Young University.
#' - Csardi, G., Nepusz, T., Traag, V., Horvat, S., Zanini, F., Noom, D., &
#'   Muller, K. (2026). *igraph: Network Analysis and Visualization*.
#'
#' @return A bundle of class `mfrm_network_analysis` containing:
#' - `summary`: graph-level connectedness and vulnerability metrics
#' - `node_metrics`: node-level degree, strength, centrality, and cutpoint flags
#' - `edge_metrics`: edge-level weights, betweenness, and bridge flags
#' - `facet_summary`: facet-level aggregation of node/bridge indicators
#' - `cut_nodes`: articulation-point rows from `node_metrics`
#' - `bridge_edges`: bridge rows from `edge_metrics`
#'
#' @seealso [subset_connectivity_report()], [diagnose_mfrm()],
#'   [mfrmr_linking_and_dff], [mfrmr_visual_diagnostics]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   net <- mfrm_network_analysis(fit)
#'   net$summary
#'   head(net$node_metrics)
#'   net$cut_nodes
#'   plot(net, type = "centrality", draw = FALSE)
#' }
#' }
#' @export
mfrm_network_analysis <- function(fit,
                                  diagnostics = NULL,
                                  top_n_subsets = NULL,
                                  min_observations = 0,
                                  include_graph = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("`mfrm_network_analysis()` requires the `igraph` package ",
         "(in Suggests). Install it and retry.", call. = FALSE)
  }

  sc <- subset_connectivity_report(
    fit = fit,
    diagnostics = diagnostics,
    top_n_subsets = top_n_subsets,
    min_observations = min_observations
  )
  nodes_tbl <- as.data.frame(sc$nodes %||% data.frame(), stringsAsFactors = FALSE)
  edges_tbl <- as.data.frame(sc$edges %||% data.frame(), stringsAsFactors = FALSE)

  empty_summary <- data.frame(
    Nodes = 0L,
    Edges = 0L,
    Components = NA_integer_,
    LargestComponentNodes = NA_integer_,
    LargestComponentShare = NA_real_,
    Density = NA_real_,
    MeanDegree = NA_real_,
    MeanStrength = NA_real_,
    ArticulationPoints = NA_integer_,
    Bridges = NA_integer_,
    Connected = NA,
    Diameter = NA_real_,
    MeanDistance = NA_real_,
    stringsAsFactors = FALSE
  )
  empty_out <- list(
    summary = empty_summary,
    node_metrics = data.frame(),
    edge_metrics = data.frame(),
    facet_summary = data.frame(),
    cut_nodes = data.frame(),
    bridge_edges = data.frame(),
    source_connectivity = sc,
    caveats = data.frame(
      Area = "network",
      Severity = "high",
      Message = "No node/edge graph could be constructed from the fitted design.",
      stringsAsFactors = FALSE
    ),
    settings = list(
      top_n_subsets = top_n_subsets %||% NA_integer_,
      min_observations = min_observations,
      include_graph = isTRUE(include_graph),
      graph_definition = "undirected weighted co-observation graph"
    )
  )
  if (nrow(nodes_tbl) == 0L || nrow(edges_tbl) == 0L ||
      !all(c("Node", "Facet", "Level", "Subset") %in% names(nodes_tbl)) ||
      !all(c("From", "To", "Weight") %in% names(edges_tbl))) {
    return(as_mfrm_bundle(empty_out, "mfrm_network_analysis"))
  }

  vertices <- nodes_tbl |>
    dplyr::mutate(
      Node = as.character(.data$Node),
      Facet = as.character(.data$Facet),
      Level = as.character(.data$Level),
      Subset = suppressWarnings(as.integer(.data$Subset))
    ) |>
    dplyr::distinct(.data$Node, .keep_all = TRUE) |>
    dplyr::arrange(.data$Subset, .data$Facet, .data$Level)

  edges <- edges_tbl |>
    dplyr::mutate(
      From = as.character(.data$From),
      To = as.character(.data$To),
      Weight = suppressWarnings(as.numeric(.data$Weight)),
      DistanceWeight = 1 / pmax(suppressWarnings(as.numeric(.data$Weight)), 1)
    ) |>
    dplyr::filter(.data$From %in% vertices$Node, .data$To %in% vertices$Node) |>
    dplyr::arrange(.data$Subset, dplyr::desc(.data$Weight), .data$From, .data$To)
  if (nrow(edges) == 0L) {
    return(as_mfrm_bundle(empty_out, "mfrm_network_analysis"))
  }

  graph_edges <- edges |>
    dplyr::select("From", "To", dplyr::everything())
  graph <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = FALSE,
    vertices = vertices
  )
  comp <- igraph::components(graph)
  comp_membership <- as.integer(comp$membership)
  comp_size <- as.integer(comp$csize[comp_membership])
  node_names <- igraph::V(graph)$name
  edge_weights <- suppressWarnings(as.numeric(igraph::E(graph)$Weight))
  edge_dist <- suppressWarnings(as.numeric(igraph::E(graph)$DistanceWeight))
  edge_dist[!is.finite(edge_dist) | edge_dist <= 0] <- 1

  articulation_names <- igraph::as_ids(igraph::articulation_points(graph))
  bridge_ids <- as.integer(igraph::bridges(graph))
  degree <- igraph::degree(graph, mode = "all", loops = FALSE)
  strength <- igraph::strength(graph, mode = "all", weights = edge_weights)
  betweenness <- igraph::betweenness(
    graph,
    directed = FALSE,
    weights = edge_dist,
    normalized = igraph::vcount(graph) > 2L
  )
  closeness <- suppressWarnings(igraph::closeness(
    graph,
    mode = "all",
    weights = edge_dist,
    normalized = TRUE
  ))

  node_metrics <- data.frame(
    Node = node_names,
    Facet = as.character(igraph::V(graph)$Facet),
    Level = as.character(igraph::V(graph)$Level),
    Subset = suppressWarnings(as.integer(igraph::V(graph)$Subset)),
    Component = comp_membership,
    ComponentSize = comp_size,
    Degree = as.numeric(degree),
    Strength = as.numeric(strength),
    Betweenness = as.numeric(betweenness),
    Closeness = as.numeric(closeness),
    IsArticulationPoint = node_names %in% articulation_names,
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(
      dplyr::desc(.data$IsArticulationPoint),
      dplyr::desc(.data$Betweenness),
      dplyr::desc(.data$Strength),
      .data$Facet,
      .data$Level
    )

  edge_metrics <- igraph::as_data_frame(graph, what = "edges")
  edge_metrics$EdgeId <- seq_len(nrow(edge_metrics))
  names(edge_metrics)[names(edge_metrics) == "from"] <- "From"
  names(edge_metrics)[names(edge_metrics) == "to"] <- "To"
  edge_metrics$Weight <- suppressWarnings(as.numeric(edge_metrics$Weight))
  edge_metrics$DistanceWeight <- suppressWarnings(as.numeric(edge_metrics$DistanceWeight))
  edge_metrics$EdgeBetweenness <- as.numeric(igraph::edge_betweenness(
    graph,
    directed = FALSE,
    weights = edge_dist
  ))
  edge_metrics$IsBridge <- edge_metrics$EdgeId %in% bridge_ids
  edge_metrics <- edge_metrics |>
    dplyr::arrange(
      dplyr::desc(.data$IsBridge),
      dplyr::desc(.data$EdgeBetweenness),
      dplyr::desc(.data$Weight),
      .data$From,
      .data$To
    ) |>
    as.data.frame(stringsAsFactors = FALSE)

  bridge_incident <- edge_metrics[edge_metrics$IsBridge, , drop = FALSE]
  bridge_facet_counts <- if (nrow(bridge_incident) > 0L &&
                             all(c("FromFacet", "ToFacet") %in% names(bridge_incident))) {
    all_facets <- c(as.character(bridge_incident$FromFacet), as.character(bridge_incident$ToFacet))
    as.data.frame(table(Facet = all_facets), stringsAsFactors = FALSE)
  } else {
    data.frame(Facet = character(), Freq = integer(), stringsAsFactors = FALSE)
  }

  facet_summary <- node_metrics |>
    dplyr::group_by(.data$Facet) |>
    dplyr::summarise(
      Levels = dplyr::n(),
      MeanDegree = mean(.data$Degree, na.rm = TRUE),
      MeanStrength = mean(.data$Strength, na.rm = TRUE),
      MaxBetweenness = max(.data$Betweenness, na.rm = TRUE),
      MeanCloseness = mean(.data$Closeness, na.rm = TRUE),
      ArticulationPoints = sum(.data$IsArticulationPoint, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      bridge_facet_counts |>
        dplyr::rename(BridgeIncidentEdges = "Freq"),
      by = "Facet"
    ) |>
    dplyr::mutate(
      BridgeIncidentEdges = dplyr::if_else(
        is.na(.data$BridgeIncidentEdges),
        0L,
        as.integer(.data$BridgeIncidentEdges)
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$ArticulationPoints),
      dplyr::desc(.data$BridgeIncidentEdges),
      dplyr::desc(.data$MaxBetweenness)
    ) |>
    as.data.frame(stringsAsFactors = FALSE)

  connected <- igraph::is_connected(graph)
  diameter <- if (isTRUE(connected)) {
    as.numeric(igraph::diameter(graph, directed = FALSE, weights = edge_dist))
  } else {
    NA_real_
  }
  mean_distance <- if (isTRUE(connected)) {
    as.numeric(igraph::mean_distance(graph, directed = FALSE, weights = edge_dist))
  } else {
    NA_real_
  }
  summary_tbl <- data.frame(
    Nodes = igraph::vcount(graph),
    Edges = igraph::ecount(graph),
    Components = as.integer(comp$no),
    LargestComponentNodes = max(comp$csize),
    LargestComponentShare = max(comp$csize) / igraph::vcount(graph),
    Density = igraph::edge_density(graph, loops = FALSE),
    MeanDegree = mean(node_metrics$Degree, na.rm = TRUE),
    MeanStrength = mean(node_metrics$Strength, na.rm = TRUE),
    ArticulationPoints = length(articulation_names),
    Bridges = length(bridge_ids),
    Connected = isTRUE(connected),
    Diameter = diameter,
    MeanDistance = mean_distance,
    stringsAsFactors = FALSE
  )

  caveats <- data.frame()
  if (!isTRUE(connected)) {
    caveats <- rbind(caveats, data.frame(
      Area = "connectedness",
      Severity = "high",
      Message = "The design graph has more than one connected component; measures may require explicit linking or anchoring before being interpreted on one scale.",
      stringsAsFactors = FALSE
    ))
  }
  if (length(articulation_names) > 0L) {
    caveats <- rbind(caveats, data.frame(
      Area = "node_vulnerability",
      Severity = "review",
      Message = paste0(length(articulation_names), " articulation point(s) indicate levels whose removal would fragment the design graph."),
      stringsAsFactors = FALSE
    ))
  }
  if (length(bridge_ids) > 0L) {
    caveats <- rbind(caveats, data.frame(
      Area = "edge_vulnerability",
      Severity = "review",
      Message = paste0(length(bridge_ids), " bridge edge(s) indicate one-link dependencies between graph regions."),
      stringsAsFactors = FALSE
    ))
  }

  out <- list(
    summary = summary_tbl,
    node_metrics = as.data.frame(node_metrics, stringsAsFactors = FALSE),
    edge_metrics = as.data.frame(edge_metrics, stringsAsFactors = FALSE),
    facet_summary = facet_summary,
    cut_nodes = as.data.frame(
      node_metrics[node_metrics$IsArticulationPoint, , drop = FALSE],
      stringsAsFactors = FALSE
    ),
    bridge_edges = as.data.frame(
      edge_metrics[edge_metrics$IsBridge, , drop = FALSE],
      stringsAsFactors = FALSE
    ),
    source_connectivity = sc,
    caveats = caveats,
    settings = list(
      top_n_subsets = top_n_subsets %||% NA_integer_,
      min_observations = min_observations,
      include_graph = isTRUE(include_graph),
      graph_definition = "undirected weighted co-observation graph",
      weight_interpretation = "Weight is the number of observations in which the two levels co-occur; DistanceWeight = 1 / max(Weight, 1)."
    )
  )
  if (isTRUE(include_graph)) {
    out$graph <- graph
  }
  as_mfrm_bundle(out, "mfrm_network_analysis")
}

network_review_top_rows <- function(x, top_n = 10) {
  x <- as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(x) == 0L) return(x)
  top_n <- max(1L, as.integer(top_n))
  x[seq_len(min(nrow(x), top_n)), , drop = FALSE]
}

network_review_sparse_design_table <- function(sparse_design = NULL) {
  if (is.null(sparse_design)) {
    return(data.frame())
  }
  if (is.data.frame(sparse_design)) {
    sparse_attr <- attr(sparse_design, "mfrm_sparse_design")
    if (is.list(sparse_attr) && is.data.frame(sparse_attr$overview)) {
      tbl <- as.data.frame(sparse_attr$overview, stringsAsFactors = FALSE)
    } else {
      tbl <- as.data.frame(sparse_design, stringsAsFactors = FALSE)
    }
  } else if (is.list(sparse_design) && is.data.frame(sparse_design$overview)) {
    tbl <- as.data.frame(sparse_design$overview, stringsAsFactors = FALSE)
  } else {
    stop(
      "`sparse_design` must be `NULL`, a sparse-design overview data frame, ",
      "or the `mfrm_sparse_design` attribute from simulate_mfrm_data().",
      call. = FALSE
    )
  }
  if (nrow(tbl) == 0L) {
    return(tbl)
  }
  if (!"SparseDesignActive" %in% names(tbl)) {
    tbl$SparseDesignActive <- if ("Active" %in% names(tbl)) {
      simulation_sparse_design_active(tbl$Active)
    } else {
      TRUE
    }
  }
  tbl
}

network_review_peer_review_table <- function(peer_review_design = NULL) {
  if (is.null(peer_review_design)) {
    return(data.frame())
  }
  if (is.data.frame(peer_review_design)) {
    peer_attr <- attr(peer_review_design, "mfrm_peer_review_design")
    if (is.list(peer_attr) && is.data.frame(peer_attr$overview)) {
      tbl <- as.data.frame(peer_attr$overview, stringsAsFactors = FALSE)
    } else {
      tbl <- as.data.frame(peer_review_design, stringsAsFactors = FALSE)
    }
  } else if (is.list(peer_review_design) && is.data.frame(peer_review_design$overview)) {
    tbl <- as.data.frame(peer_review_design$overview, stringsAsFactors = FALSE)
  } else {
    stop(
      "`peer_review_design` must be `NULL`, a peer-review overview data frame, ",
      "or the `mfrm_peer_review_design` attribute from simulate_mfrm_data().",
      call. = FALSE
    )
  }
  if (nrow(tbl) == 0L) {
    return(tbl)
  }
  if (!"Active" %in% names(tbl)) {
    tbl$Active <- TRUE
  }
  if (!"Scenario" %in% names(tbl)) {
    tbl$Scenario <- "peer_review"
  }
  if (!"ReviewUse" %in% names(tbl)) {
    tbl$ReviewUse <- "design_diagnostic_not_measurement_gate"
  }
  tbl
}

peer_review_design_bundle <- function(peer_review_design) {
  if (is.null(peer_review_design)) {
    stop("`peer_review_design` must not be `NULL`.", call. = FALSE)
  }
  if (is.data.frame(peer_review_design)) {
    peer_attr <- attr(peer_review_design, "mfrm_peer_review_design")
    if (is.list(peer_attr) && is.data.frame(peer_attr$overview)) {
      src <- peer_attr
    } else {
      src <- list(overview = as.data.frame(peer_review_design, stringsAsFactors = FALSE))
    }
  } else if (is.list(peer_review_design) && is.data.frame(peer_review_design$overview)) {
    src <- peer_review_design
  } else {
    stop(
      "`peer_review_design` must be a generated data frame carrying ",
      "`mfrm_peer_review_design`, the attribute itself, or a peer-review ",
      "overview data frame.",
      call. = FALSE
    )
  }
  list(
    overview = as.data.frame(src$overview %||% data.frame(), stringsAsFactors = FALSE),
    submission_load = as.data.frame(src$submission_load %||% data.frame(), stringsAsFactors = FALSE),
    reviewer_load = as.data.frame(src$reviewer_load %||% data.frame(), stringsAsFactors = FALSE),
    reviewer_pair_common_submissions = as.data.frame(
      src$reviewer_pair_common_submissions %||% data.frame(),
      stringsAsFactors = FALSE
    ),
    review_pairs = as.data.frame(src$review_pairs %||% data.frame(), stringsAsFactors = FALSE),
    notes = clean_summary_lines(src$notes %||% character(0))
  )
}

peer_review_design_status <- function(overview) {
  overview <- as.data.frame(overview %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(overview) == 0L) {
    return(data.frame(
      PeerReviewStatus = "insufficient_data",
      PeerReviewReason = "No peer-review design overview is available.",
      ReviewUse = "design_diagnostic_not_measurement_gate",
      stringsAsFactors = FALSE
    ))
  }
  active <- if ("Active" %in% names(overview)) {
    simulation_sparse_design_active(overview$Active[1])
  } else {
    TRUE
  }
  self_reviews <- suppressWarnings(as.integer(overview$SelfReviews[1] %||% NA_integer_))
  avoid_self <- if ("AvoidSelfReview" %in% names(overview)) {
    simulation_sparse_design_active(overview$AvoidSelfReview[1])
  } else {
    NA
  }
  min_reviewer_load <- suppressWarnings(as.numeric(overview$MinSubmissionsPerReviewer[1] %||% NA_real_))
  min_common <- suppressWarnings(as.numeric(overview$MinCommonSubmissionsPerReviewerPair[1] %||% NA_real_))
  zero_common <- suppressWarnings(as.integer(overview$ZeroCommonReviewerPairs[1] %||% NA_integer_))

  status <- "ok"
  reason <- "Recorded peer-review assignment checks are satisfied."
  if (!isTRUE(active)) {
    status <- "insufficient_data"
    reason <- "The supplied metadata does not describe an active peer-review design."
  } else if (isTRUE(avoid_self) && is.finite(self_reviews) && self_reviews > 0L) {
    status <- "warning"
    reason <- "Self-review was requested to be excluded, but at least one self-review pair is present."
  } else if (is.finite(min_reviewer_load) && min_reviewer_load <= 0) {
    status <- "review"
    reason <- "At least one reviewer has no assigned submissions."
  } else if (is.finite(zero_common) && zero_common > 0L) {
    status <- "review"
    reason <- "At least one reviewer pair has no common submissions."
  } else if (!is.finite(min_common)) {
    status <- "review"
    reason <- "Reviewer-pair common-submission counts are unavailable."
  }

  data.frame(
    PeerReviewStatus = status,
    PeerReviewReason = reason,
    ReviewUse = "design_diagnostic_not_measurement_gate",
    stringsAsFactors = FALSE
  )
}

peer_review_load_summary <- function(submission_load, reviewer_load) {
  submission_load <- as.data.frame(submission_load %||% data.frame(), stringsAsFactors = FALSE)
  reviewer_load <- as.data.frame(reviewer_load %||% data.frame(), stringsAsFactors = FALSE)
  sub_counts <- if ("ReviewersAssigned" %in% names(submission_load)) {
    suppressWarnings(as.numeric(submission_load$ReviewersAssigned))
  } else {
    numeric(0)
  }
  rev_counts <- if ("SubmissionsReviewed" %in% names(reviewer_load)) {
    suppressWarnings(as.numeric(reviewer_load$SubmissionsReviewed))
  } else {
    numeric(0)
  }
  safe_min <- function(x) if (length(x) == 0L || all(!is.finite(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_mean <- function(x) if (length(x) == 0L || all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
  safe_max <- function(x) if (length(x) == 0L || all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE)
  data.frame(
    Submissions = nrow(submission_load),
    Reviewers = nrow(reviewer_load),
    MinReviewersPerSubmission = safe_min(sub_counts),
    MeanReviewersPerSubmission = safe_mean(sub_counts),
    MaxReviewersPerSubmission = safe_max(sub_counts),
    MinSubmissionsPerReviewer = safe_min(rev_counts),
    MeanSubmissionsPerReviewer = safe_mean(rev_counts),
    MaxSubmissionsPerReviewer = safe_max(rev_counts),
    LoadReviewUse = "assignment_diagnostic_not_quality_gate",
    stringsAsFactors = FALSE
  )
}

peer_review_reciprocal_pair_rows <- function(review_pairs) {
  review_pairs <- as.data.frame(review_pairs %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(review_pairs) == 0L || !all(c("Person", "Reviewer") %in% names(review_pairs))) {
    return(data.frame())
  }
  pairs <- unique(review_pairs[, c("Person", "Reviewer"), drop = FALSE])
  pairs$Person <- as.character(pairs$Person)
  pairs$Reviewer <- as.character(pairs$Reviewer)
  key <- paste(pairs$Person, pairs$Reviewer, sep = "\r")
  reciprocal <- paste(pairs$Reviewer, pairs$Person, sep = "\r")
  idx <- key %in% reciprocal & pairs$Person != pairs$Reviewer
  if (!any(idx)) {
    return(data.frame())
  }
  p1 <- pmin(pairs$Person[idx], pairs$Reviewer[idx])
  p2 <- pmax(pairs$Person[idx], pairs$Reviewer[idx])
  unique(data.frame(
    Participant1 = p1,
    Participant2 = p2,
    ReciprocalReviewPair = TRUE,
    ReviewUse = "design_diagnostic_not_measurement_gate",
    stringsAsFactors = FALSE
  ))
}

peer_review_reporting_map <- function() {
  data.frame(
    Area = c(
      "Peer-review assignment",
      "Reviewer load",
      "Common-submission links",
      "Reciprocal review pairs",
      "MFRM measurement model"
    ),
    PrimaryTable = c(
      "overview",
      "load_summary; reviewer_load; submission_load",
      "reviewer_pair_common_submissions; low_common_pairs",
      "reciprocal_pairs",
      "fit_mfrm(); diagnose_mfrm()"
    ),
    Use = c(
      "Check self-review exclusion, design density, anchor counts, and front-door assignment status.",
      "Inspect whether submissions and reviewers received the intended assignment load.",
      "Inspect reviewer linkage through shared reviewed submissions.",
      "Flag participant pairs that reviewed each other's submissions for design transparency.",
      "Estimate and diagnose person/submission measures, reviewer severity, criterion difficulty, fit, and precision."
    ),
    Boundary = c(
      "Design diagnostic; not evidence of peer quality or rating fairness.",
      "Assignment diagnostic; not a reviewer-quality score.",
      "Design-link diagnostic; not a universal adequacy threshold.",
      "Design-transparency diagnostic; not automatically a bias finding.",
      "Measurement model; keep separate from assignment-design diagnostics."
    ),
    stringsAsFactors = FALSE
  )
}

peer_review_top_rows <- function(x, top_n = 10) {
  x <- as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(x) == 0L) return(x)
  top_n <- max(1L, as.integer(top_n))
  x[seq_len(min(nrow(x), top_n)), , drop = FALSE]
}

#' Build a peer-review design review
#'
#' @param peer_review_design A generated data frame carrying the
#'   `mfrm_peer_review_design` attribute, the attribute itself, or its
#'   `overview` data frame.
#' @param top_n Number of reviewer-load, submission-load, common-link, and
#'   reciprocal-pair rows to keep in compact summary tables.
#'
#' @details
#' `build_peer_review_design_review()` converts peer-review simulation metadata
#' into a reportable design-review object. The review summarizes self-review
#' checks, reviewer and submission load, common submissions per reviewer pair,
#' and reciprocal review pairs. These rows are assignment-design diagnostics:
#' they do not replace MFRM estimates, fit, separation, reliability, or
#' substantive review-quality evidence.
#'
#' Peer-review use of MFRM follows studies that model peer/self/teacher rater
#' severity and leniency. Common-link anchor interpretation follows sparse
#' rater-mediated design work; the review status is therefore descriptive and
#' conservative rather than a literature-derived universal adequacy cutoff.
#'
#' @section References:
#' - Farrokhi, F., Esfandiari, R., & Schaefer, E. (2012). A many-facet Rasch
#'   measurement of differential rater severity/leniency in three types of
#'   assessment. *JALT Journal*, 34(1), 79-102.
#'   doi:10.37546/JALTJJ34.1-3.
#' - Uto, M., & Ueno, M. (2020). A generalized many-facet Rasch model and its
#'   Bayesian estimation using Hamiltonian Monte Carlo. *Behaviormetrika*,
#'   47, 469-496. doi:10.1007/s41237-020-00115-7.
#' - DeMars, C. E., Shapovalov, Y. A., & Hathcoat, J. D. (2023).
#'   *Many-Facet Rasch Designs: How Should Raters be Assigned to Examinees?*
#'   NCME presentation.
#'
#' @return A bundle of class `mfrm_peer_review_design_review`.
#' @seealso [build_peer_review_sim_spec()], [simulate_mfrm_data()],
#'   [build_mfrm_network_review()], [build_summary_table_bundle()]
#' @examples
#' peer_spec <- build_peer_review_sim_spec(
#'   n_submission = 12,
#'   n_criterion = 3,
#'   reviewers_per_submission = 2,
#'   anchor_submissions = 2
#' )
#' peer_sim <- simulate_mfrm_data(sim_spec = peer_spec, seed = 123)
#' review <- build_peer_review_design_review(peer_sim)
#' summary(review)$overview
#' @export
build_peer_review_design_review <- function(peer_review_design, top_n = 10) {
  top_n <- max(1L, as.integer(top_n))
  src <- peer_review_design_bundle(peer_review_design)
  overview <- src$overview
  status <- peer_review_design_status(overview)
  if ("ReviewUse" %in% names(overview)) {
    overview$ReviewUse <- NULL
  }
  overview <- data.frame(overview, status, check.names = FALSE, stringsAsFactors = FALSE)

  load_summary <- peer_review_load_summary(src$submission_load, src$reviewer_load)
  reviewer_load <- src$reviewer_load
  if (nrow(reviewer_load) > 0L && "SubmissionsReviewed" %in% names(reviewer_load)) {
    reviewer_load <- reviewer_load |>
      dplyr::arrange(dplyr::desc(.data$SubmissionsReviewed), .data$Reviewer)
  }
  submission_load <- src$submission_load
  if (nrow(submission_load) > 0L && "ReviewersAssigned" %in% names(submission_load)) {
    submission_load <- submission_load |>
      dplyr::arrange(dplyr::desc(.data$ReviewersAssigned), .data$Person)
  }
  common_pairs <- src$reviewer_pair_common_submissions
  low_common <- common_pairs
  if (nrow(low_common) > 0L && "CommonSubmissions" %in% names(low_common)) {
    low_common <- low_common |>
      dplyr::arrange(.data$CommonSubmissions, .data$Reviewer1, .data$Reviewer2)
  }
  reciprocal <- peer_review_reciprocal_pair_rows(src$review_pairs)

  caveats <- data.frame(
    Area = "interpretation_boundary",
    Severity = "info",
    Message = "Peer-review design diagnostics summarize assignment structure and linkage; they are not MFRM fit statistics, reviewer-quality estimates, or fairness evidence.",
    stringsAsFactors = FALSE
  )
  out <- list(
    overview = overview,
    load_summary = load_summary,
    submission_load = peer_review_top_rows(submission_load, top_n = top_n),
    reviewer_load = peer_review_top_rows(reviewer_load, top_n = top_n),
    reviewer_pair_common_submissions = peer_review_top_rows(common_pairs, top_n = top_n),
    low_common_pairs = peer_review_top_rows(low_common, top_n = top_n),
    reciprocal_pairs = peer_review_top_rows(reciprocal, top_n = top_n),
    reporting_map = peer_review_reporting_map(),
    caveats = caveats,
    settings = list(
      top_n = top_n,
      review_use = "design_diagnostic_not_measurement_gate"
    ),
    notes = c(
      src$notes,
      "Peer-review design status is a routing aid for assignment review, not a measurement adequacy decision."
    )
  )
  as_mfrm_bundle(out, "mfrm_peer_review_design_review")
}

#' @export
print.mfrm_peer_review_design_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Summarize a peer-review design review
#'
#' @param object Output from [build_peer_review_design_review()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of rows to keep in compact follow-up tables.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_peer_review_design_review`.
#' @seealso [build_peer_review_design_review()]
#' @export
summary.mfrm_peer_review_design_review <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_peer_review_design_review")) {
    stop("`object` must be output from build_peer_review_design_review().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    load_summary = tibble::as_tibble(object$load_summary %||% tibble::tibble()),
    submission_load = tibble::as_tibble(peer_review_top_rows(object$submission_load, top_n = top_n)),
    reviewer_load = tibble::as_tibble(peer_review_top_rows(object$reviewer_load, top_n = top_n)),
    reviewer_pair_common_submissions = tibble::as_tibble(
      peer_review_top_rows(object$reviewer_pair_common_submissions, top_n = top_n)
    ),
    low_common_pairs = tibble::as_tibble(peer_review_top_rows(object$low_common_pairs, top_n = top_n)),
    reciprocal_pairs = tibble::as_tibble(peer_review_top_rows(object$reciprocal_pairs, top_n = top_n)),
    reporting_map = tibble::as_tibble(object$reporting_map %||% tibble::tibble()),
    caveats = tibble::as_tibble(object$caveats %||% tibble::tibble()),
    notes = clean_summary_lines(object$notes %||% character(0)),
    settings = object$settings %||% list(),
    digits = digits
  )
  class(out) <- "summary.mfrm_peer_review_design_review"
  out
}

#' @export
print.summary.mfrm_peer_review_design_review <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Peer-Review Design Review Summary\n")
  if (is.data.frame(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$load_summary) && nrow(x$load_summary) > 0L) {
    cat("\nLoad summary\n")
    print(round_numeric_df(as.data.frame(x$load_summary), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$low_common_pairs) && nrow(x$low_common_pairs) > 0L) {
    cat("\nLowest common-submission reviewer pairs\n")
    print(round_numeric_df(as.data.frame(x$low_common_pairs), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$reciprocal_pairs) && nrow(x$reciprocal_pairs) > 0L) {
    cat("\nReciprocal review pairs\n")
    print(as.data.frame(x$reciprocal_pairs), row.names = FALSE)
  }
  if (is.data.frame(x$caveats) && nrow(x$caveats) > 0L) {
    cat("\nCaveats\n")
    print(as.data.frame(x$caveats), row.names = FALSE)
  }
  if (length(x$notes) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

network_review_status <- function(summary_tbl) {
  summary_tbl <- as.data.frame(summary_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_tbl) == 0L) {
    return(data.frame(
      NetworkReviewStatus = "insufficient_data",
      NetworkReviewReason = "No design-network summary is available.",
      ReviewUse = "design_diagnostic_not_measurement_gate",
      stringsAsFactors = FALSE
    ))
  }

  nodes <- suppressWarnings(as.integer(summary_tbl$Nodes[1] %||% NA_integer_))
  edges <- suppressWarnings(as.integer(summary_tbl$Edges[1] %||% NA_integer_))
  components <- suppressWarnings(as.integer(summary_tbl$Components[1] %||% NA_integer_))
  connected <- if ("Connected" %in% names(summary_tbl)) {
    as.logical(summary_tbl$Connected[1])
  } else {
    NA
  }
  articulation <- suppressWarnings(as.integer(summary_tbl$ArticulationPoints[1] %||% NA_integer_))
  bridges <- suppressWarnings(as.integer(summary_tbl$Bridges[1] %||% NA_integer_))

  status <- "ok"
  reason <- "The observed design graph is connected and has no recorded articulation points or bridge edges."
  if (!is.finite(nodes) || !is.finite(edges) || nodes == 0L || edges == 0L) {
    status <- "insufficient_data"
    reason <- "No node/edge design graph could be constructed from the fitted data."
  } else if (identical(connected, FALSE) || (is.finite(components) && components > 1L)) {
    status <- "warning"
    reason <- "The design graph has more than one connected component; interpret common-scale claims only with explicit linking or anchoring support."
  } else if ((is.finite(articulation) && articulation > 0L) ||
             (is.finite(bridges) && bridges > 0L)) {
    status <- "review"
    reason <- "The design graph is connected but contains articulation points or bridge edges that indicate linking vulnerability."
  }

  data.frame(
    NetworkReviewStatus = status,
    NetworkReviewReason = reason,
    ReviewUse = "design_diagnostic_not_measurement_gate",
    stringsAsFactors = FALSE
  )
}

network_review_overview <- function(summary_tbl) {
  summary_tbl <- as.data.frame(summary_tbl %||% data.frame(), stringsAsFactors = FALSE)
  status <- network_review_status(summary_tbl)
  if (nrow(summary_tbl) == 0L) {
    base <- data.frame(
      Nodes = NA_integer_,
      Edges = NA_integer_,
      Components = NA_integer_,
      Connected = NA,
      ArticulationPoints = NA_integer_,
      Bridges = NA_integer_,
      LargestComponentShare = NA_real_,
      Density = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    keep <- intersect(
      c(
        "Nodes", "Edges", "Components", "Connected", "ArticulationPoints",
        "Bridges", "LargestComponentShare", "Density"
      ),
      names(summary_tbl)
    )
    base <- summary_tbl[1, keep, drop = FALSE]
  }
  data.frame(base, status, check.names = FALSE, stringsAsFactors = FALSE)
}

network_review_reporting_map <- function() {
  data.frame(
    Area = c(
      "MFRM measurement model",
      "Design network",
      "Sparse linked design",
      "Peer-review design",
      "Rater-effect network"
    ),
    PrimaryHelper = c(
      "fit_mfrm(); diagnose_mfrm()",
      "mfrm_network_analysis(); build_mfrm_network_review()",
      "simulate_mfrm_data(..., assignment = \"sparse_linked\")",
      "build_peer_review_sim_spec(); simulate_mfrm_data()",
      "rater_network_analysis(); rater_halo_network_analysis()"
    ),
    Use = c(
      "Estimate and diagnose person/facet measures, fit, precision, and bias.",
      "Inspect observed co-observation connectedness, articulation points, bridge edges, and design-level linking vulnerability.",
      "Inspect planned missingness, rater coverage, and common-person links in sparse simulation designs.",
      "Simulate peer-review or peer-assessment assignments with shared submission/reviewer IDs, no-self-review checks, reviewer-load diagnostics, and common-submission reviewer links.",
      "Screen observed rater relationship or halo patterns as descriptive diagnostics."
    ),
    Boundary = c(
      "Primary measurement model.",
      "Design diagnostic; not person ability, rater quality, or formal fit.",
      "Design diagnostic; not a recovery, fit, or separation gate.",
      "Design diagnostic; not peer quality, reviewer fairness, or a universal common-link adequacy threshold.",
      "Descriptive network diagnostic; not a Rasch logit estimate or causal halo conclusion."
    ),
    stringsAsFactors = FALSE
  )
}

#' Build an MFRM network review
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param sparse_design Optional sparse-design metadata. Supply either the
#'   generated data frame that carries the `mfrm_sparse_design` attribute, the
#'   attribute itself, or a data frame with sparse design columns such as
#'   `SparseDesignActive`, `DesignDensity`, `MinCommonPersonsPerRaterPair`,
#'   and `ZeroCommonRaterPairs`.
#' @param peer_review_design Optional peer-review design metadata. Supply
#'   either the generated data frame that carries the `mfrm_peer_review_design`
#'   attribute, the attribute itself, or its `overview` data frame.
#' @param top_n_subsets Optional maximum number of connected-subset rows to
#'   retain before constructing the graph; passed to [mfrm_network_analysis()].
#' @param min_observations Minimum observations required to keep a subset row;
#'   passed to [mfrm_network_analysis()].
#' @param top_n Number of central/cut/bridge rows to retain in the review.
#' @param include_graph Logical; if `TRUE`, keep the underlying `igraph` object
#'   in the nested `source_network` bundle.
#'
#' @details
#' `build_mfrm_network_review()` is a synthesis layer over
#' [mfrm_network_analysis()]. It keeps the measurement model and graph view in
#' separate lanes: MFRM estimates remain the measurement results, while the
#' network review summarizes co-observation connectedness and linking
#' vulnerability in the observed design. This is especially useful for sparse
#' or incomplete rater-mediated designs, where common-person links, connected
#' subsets, articulation points, and bridge edges can explain why an otherwise
#' estimable model depends on fragile design links.
#'
#' The review status is deliberately conservative and descriptive. It is not a
#' literature-derived adequacy cut point for fit, separation, recovery, or
#' rater quality. Use it to decide which design links, anchors, or additional
#' observations need inspection before making common-scale claims.
#'
#' @section References:
#' - Wind, S. A., & Jones, E. (2018). The stabilizing influences of linking set
#'   size and model-data fit in sparse rater-mediated assessment networks.
#'   *Educational and Psychological Measurement*. doi:10.1177/0013164417703733.
#' - Wind, S. A., Jones, E., & Grajeda, S. (2023). Does sparseness matter?
#'   Examining the use of generalizability theory and many-facet Rasch
#'   measurement in sparse rating designs. *Applied Psychological
#'   Measurement*, 47(5-6), 351-364. doi:10.1177/01466216231182148.
#' - DeMars, C. E., Shapovalov, Y. A., & Hathcoat, J. D. (2023).
#'   *Many-Facet Rasch Designs: How Should Raters be Assigned to Examinees?*
#'   NCME presentation.
#'
#' @return A bundle of class `mfrm_network_review` containing:
#' - `overview`: connectedness and front-door review status
#' - `network_summary`: graph-level metrics from [mfrm_network_analysis()]
#' - `facet_summary`: facet-level vulnerability summaries
#' - `top_central_nodes`, `top_cut_nodes`, `top_bridge_edges`: follow-up rows
#' - `sparse_review`: optional sparse-design linking review
#' - `peer_review`: optional peer-review assignment and linkage diagnostics
#' - `reporting_map`: boundary between MFRM, design network, sparse design,
#'   peer-review design, and rater-effect network routes
#'
#' @seealso [mfrm_network_analysis()], [subset_connectivity_report()],
#'   [build_summary_table_bundle()], [rater_network_analysis()],
#'   [rater_halo_network_analysis()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   review <- build_mfrm_network_review(fit)
#'   summary(review)
#'   build_summary_table_bundle(review)
#' }
#' }
#' @export
build_mfrm_network_review <- function(fit,
                                      diagnostics = NULL,
                                      sparse_design = NULL,
                                      peer_review_design = NULL,
                                      top_n_subsets = NULL,
                                      min_observations = 0,
                                      top_n = 10,
                                      include_graph = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  top_n <- max(1L, as.integer(top_n))

  net <- mfrm_network_analysis(
    fit = fit,
    diagnostics = diagnostics,
    top_n_subsets = top_n_subsets,
    min_observations = min_observations,
    include_graph = include_graph
  )
  net_summary <- as.data.frame(net$summary %||% data.frame(), stringsAsFactors = FALSE)
  sparse_tbl <- network_review_sparse_design_table(sparse_design)
  sparse_review <- simulation_sparse_design_review_summary(sparse_tbl)
  peer_tbl <- network_review_peer_review_table(peer_review_design)

  caveats <- as.data.frame(net$caveats %||% data.frame(), stringsAsFactors = FALSE)
  boundary_caveat <- data.frame(
    Area = "interpretation_boundary",
    Severity = "info",
    Message = "Network-review metrics summarize observation-design connectedness and linking vulnerability; they are not Rasch fit statistics, person measures, or rater-quality estimates.",
    stringsAsFactors = FALSE
  )
  caveats <- rbind(caveats, boundary_caveat)

  central <- as.data.frame(net$node_metrics %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(central) > 0L && all(c("Betweenness", "Strength") %in% names(central))) {
    central <- central |>
      dplyr::arrange(
        dplyr::desc(.data$Betweenness),
        dplyr::desc(.data$Strength),
        .data$Facet,
        .data$Level
      )
  }

  out <- list(
    overview = network_review_overview(net_summary),
    network_summary = net_summary,
    facet_summary = as.data.frame(net$facet_summary %||% data.frame(), stringsAsFactors = FALSE),
    top_central_nodes = network_review_top_rows(central, top_n = top_n),
    top_cut_nodes = network_review_top_rows(net$cut_nodes, top_n = top_n),
    top_bridge_edges = network_review_top_rows(net$bridge_edges, top_n = top_n),
    sparse_review = as.data.frame(sparse_review, stringsAsFactors = FALSE),
    peer_review = as.data.frame(peer_tbl, stringsAsFactors = FALSE),
    reporting_map = network_review_reporting_map(),
    caveats = caveats,
    source_network = net,
    settings = list(
      top_n = top_n,
      top_n_subsets = top_n_subsets %||% NA_integer_,
      min_observations = min_observations,
      include_graph = isTRUE(include_graph),
      review_use = "design_diagnostic_not_measurement_gate"
    ),
    notes = c(
      "MFRM estimates remain the measurement-model results; network rows summarize observed design links.",
      "Articulation points and bridge edges identify levels or links whose removal would fragment the co-observation graph.",
      "Sparse-design review rows, when supplied, report planned-missingness and rater-link diagnostics rather than recovery or fit gates.",
      "Peer-review design rows, when supplied, report assignment structure and reviewer linkage rather than reviewer quality or fairness."
    )
  )
  as_mfrm_bundle(out, "mfrm_network_review")
}

#' @export
print.mfrm_network_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Summarize an MFRM network review
#'
#' @param object Output from [build_mfrm_network_review()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of central/cut/bridge rows to keep in the compact
#'   summary.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_network_review`.
#' @seealso [build_mfrm_network_review()]
#' @export
summary.mfrm_network_review <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_network_review")) {
    stop("`object` must be output from build_mfrm_network_review().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    network_summary = tibble::as_tibble(object$network_summary %||% tibble::tibble()),
    facet_summary = tibble::as_tibble(object$facet_summary %||% tibble::tibble()),
    top_central_nodes = tibble::as_tibble(network_review_top_rows(object$top_central_nodes, top_n = top_n)),
    top_cut_nodes = tibble::as_tibble(network_review_top_rows(object$top_cut_nodes, top_n = top_n)),
    top_bridge_edges = tibble::as_tibble(network_review_top_rows(object$top_bridge_edges, top_n = top_n)),
    sparse_review = tibble::as_tibble(object$sparse_review %||% tibble::tibble()),
    peer_review = tibble::as_tibble(object$peer_review %||% tibble::tibble()),
    reporting_map = tibble::as_tibble(object$reporting_map %||% tibble::tibble()),
    caveats = tibble::as_tibble(object$caveats %||% tibble::tibble()),
    notes = clean_summary_lines(object$notes %||% character(0)),
    settings = object$settings %||% list(),
    digits = digits
  )
  class(out) <- "summary.mfrm_network_review"
  out
}

#' @export
print.summary.mfrm_network_review <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Network Review Summary\n")
  if (is.data.frame(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$network_summary) && nrow(x$network_summary) > 0L) {
    cat("\nNetwork summary\n")
    print(round_numeric_df(as.data.frame(x$network_summary), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$facet_summary) && nrow(x$facet_summary) > 0L) {
    cat("\nFacet vulnerability summary\n")
    print(round_numeric_df(as.data.frame(x$facet_summary), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$top_central_nodes) && nrow(x$top_central_nodes) > 0L) {
    cat("\nTop central nodes\n")
    print(round_numeric_df(as.data.frame(x$top_central_nodes), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$top_cut_nodes) && nrow(x$top_cut_nodes) > 0L) {
    cat("\nArticulation nodes\n")
    print(round_numeric_df(as.data.frame(x$top_cut_nodes), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$top_bridge_edges) && nrow(x$top_bridge_edges) > 0L) {
    cat("\nBridge edges\n")
    print(round_numeric_df(as.data.frame(x$top_bridge_edges), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$sparse_review) && nrow(x$sparse_review) > 0L) {
    cat("\nSparse design review\n")
    print(round_numeric_df(as.data.frame(x$sparse_review), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$peer_review) && nrow(x$peer_review) > 0L) {
    cat("\nPeer-review design\n")
    print(round_numeric_df(as.data.frame(x$peer_review), digits = digits), row.names = FALSE)
  }
  if (is.data.frame(x$caveats) && nrow(x$caveats) > 0L) {
    cat("\nCaveats\n")
    print(as.data.frame(x$caveats), row.names = FALSE)
  }
  if (length(x$notes) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

rater_network_score_wide <- function(obs_df, facet_cols, rater_facet) {
  if (is.null(obs_df) || nrow(obs_df) == 0L) {
    return(list(wide = data.frame(), raters = character(), context_cols = character()))
  }
  context_cols <- setdiff(facet_cols, rater_facet)
  if (length(context_cols) == 0L) {
    return(list(wide = data.frame(), raters = character(), context_cols = character()))
  }
  df <- obs_df |>
    dplyr::mutate(dplyr::across(dplyr::all_of(context_cols), as.character)) |>
    tidyr::unite(".context", dplyr::all_of(context_cols), sep = "|", remove = FALSE) |>
    dplyr::select(".context", dplyr::all_of(rater_facet), "Observed", dplyr::any_of("Weight"))
  df$.Weight <- get_weights(df)
  df <- df |>
    dplyr::group_by(.data$.context, .data[[rater_facet]]) |>
    dplyr::summarise(Score = weighted_mean(.data$Observed, .data$.Weight), .groups = "drop")
  if (nrow(df) == 0L) {
    return(list(wide = data.frame(), raters = character(), context_cols = context_cols))
  }
  wide <- tryCatch(
    tidyr::pivot_wider(
      df,
      id_cols = ".context",
      names_from = !!rlang::sym(rater_facet),
      values_from = "Score"
    ),
    error = function(e) NULL
  )
  if (is.null(wide)) {
    return(list(wide = data.frame(), raters = character(), context_cols = context_cols))
  }
  raters <- setdiff(names(wide), ".context")
  list(
    wide = as.data.frame(wide, stringsAsFactors = FALSE),
    raters = raters,
    context_cols = context_cols
  )
}

rater_network_direction_pairs <- function(wide, raters, score_diff_tolerance = 0) {
  if (is.null(wide) || nrow(wide) == 0L || length(raters) < 2L) {
    return(data.frame())
  }
  score_diff_tolerance <- max(0, suppressWarnings(as.numeric(score_diff_tolerance[1])))
  pairs <- utils::combn(raters, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    v1 <- suppressWarnings(as.numeric(wide[[pair[1]]]))
    v2 <- suppressWarnings(as.numeric(wide[[pair[2]]]))
    ok <- is.finite(v1) & is.finite(v2)
    n_ok <- sum(ok)
    if (n_ok == 0L) {
      return(data.frame(
        Rater1 = pair[1],
        Rater2 = pair[2],
        DirectionN = 0L,
        Rater1HigherCount = 0L,
        Rater2HigherCount = 0L,
        TiedOrWithinToleranceCount = 0L,
        Rater1HigherProp = NA_real_,
        Rater2HigherProp = NA_real_,
        TiedOrWithinToleranceProp = NA_real_,
        NetLeniencyProp = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    diff <- v1[ok] - v2[ok]
    r1_higher <- sum(diff > score_diff_tolerance, na.rm = TRUE)
    r2_higher <- sum(diff < -score_diff_tolerance, na.rm = TRUE)
    tied <- sum(abs(diff) <= score_diff_tolerance, na.rm = TRUE)
    data.frame(
      Rater1 = pair[1],
      Rater2 = pair[2],
      DirectionN = n_ok,
      Rater1HigherCount = r1_higher,
      Rater2HigherCount = r2_higher,
      TiedOrWithinToleranceCount = tied,
      Rater1HigherProp = r1_higher / n_ok,
      Rater2HigherProp = r2_higher / n_ok,
      TiedOrWithinToleranceProp = tied / n_ok,
      NetLeniencyProp = (r1_higher - r2_higher) / n_ok,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

empty_rater_network_bundle <- function(settings, source_interrater = NULL, message = NULL) {
  summary_tbl <- data.frame(
    RaterFacet = as.character(settings$rater_facet %||% NA_character_),
    Mode = as.character(settings$mode %||% NA_character_),
    Raters = 0L,
    PairRows = 0L,
    Edges = 0L,
    Directed = isTRUE(settings$directed),
    WeightMetric = as.character(settings$weight_metric %||% NA_character_),
    Density = NA_real_,
    MeanWeight = NA_real_,
    MeanDegree = NA_real_,
    MeanStrength = NA_real_,
    stringsAsFactors = FALSE
  )
  caveats <- data.frame(
    Area = "rater_network",
    Severity = "high",
    Message = message %||% "No rater network could be constructed from shared scoring contexts.",
    stringsAsFactors = FALSE
  )
  as_mfrm_bundle(list(
    summary = summary_tbl,
    node_metrics = data.frame(),
    edge_metrics = data.frame(),
    pair_metrics = data.frame(),
    caveats = caveats,
    source_interrater = source_interrater,
    settings = settings
  ), "mfrm_rater_network")
}

#' Analyze rater agreement, disagreement, and severity-direction networks
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param rater_facet Name of the rater-like facet. If omitted, mfrmr uses the
#'   same heuristic as [interrater_agreement_table()].
#' @param context_facets Facets defining shared scoring contexts. By default,
#'   the person facet and all non-rater facets are used.
#' @param mode Network definition. `"agreement"` builds an undirected network
#'   whose edge weights represent observed agreement. `"disagreement"` builds
#'   an undirected network whose edge weights represent observed disagreement.
#'   `"severity_direction"` builds a directed network: an edge from rater A to
#'   rater B means A assigned higher scores than B in shared contexts and is
#'   therefore relatively more lenient under the usual higher-score-is-better
#'   rating convention.
#' @param weight_metric Pair-level weight used for `"agreement"` or
#'   `"disagreement"` networks. Defaults to `Exact` for agreement and `MAD`
#'   for disagreement. Available pair columns include `Exact`, `Adjacent`,
#'   `Corr`, `MAD`, `OneMinusExact`, and `AbsMeanDiff`.
#' @param min_pair_n Minimum number of shared contexts required for a rater
#'   pair to contribute an edge.
#' @param min_weight Minimum edge weight retained in the graph.
#' @param score_diff_tolerance Score-difference tolerance for directed
#'   severity networks. With the default `0`, any higher score contributes to
#'   the outgoing leniency edge. Larger values reproduce thresholded
#'   disagreement displays such as "only differences greater than 3 marks".
#' @param severity_continuity Continuity constant added to incoming and
#'   outgoing strengths before computing the finite severity index
#'   `-log((OutStrength + c) / (InStrength + c))`.
#' @param exact_warn,corr_warn Passed to [interrater_agreement_table()] to keep
#'   pair flags consistent with the tabular agreement view.
#' @param include_graph If `TRUE`, include the underlying `igraph` object in the
#'   returned bundle.
#'
#' @details
#' This function implements a package-native rater-effect network view
#' complementary to MFRM output. It follows the pairwise-network logic used in
#' Lamprianou's rater-effect network work: nodes are raters, edges summarize
#' pairwise relationships among raters in shared scoring contexts, and directed
#' disagreement edges can be interpreted as relative leniency/severity
#' indicators. These network summaries are descriptive diagnostics, not Rasch
#' logit estimates and not formal fit statistics.
#'
#' For `mode = "severity_direction"`, outgoing strength means the rater more
#' often assigned higher scores than comparison raters; incoming strength means
#' comparison raters more often assigned higher scores than this rater. The
#' reported `SeverityIndex` is positive for relatively severe raters and
#' negative for relatively lenient raters, but it is on a network-analysis scale
#' and should not be read as an MFRM severity logit.
#'
#' @return A bundle of class `mfrm_rater_network` containing:
#' \describe{
#'   \item{`summary`}{One-row graph summary.}
#'   \item{`node_metrics`}{Rater-level degree, strength, centrality, and
#'     severity-direction summaries.}
#'   \item{`edge_metrics`}{Retained rater-pair network edges.}
#'   \item{`pair_metrics`}{All eligible pairwise agreement and directional
#'     comparison metrics before edge thresholding.}
#'   \item{`caveats`}{Interpretation notes and sparse-design warnings.}
#'   \item{`source_interrater`}{The underlying [interrater_agreement_table()]
#'     output used for agreement statistics.}
#' }
#'
#' @references
#' Lamprianou, I. (2018). Investigation of rater effects using Social Network
#' Analysis and Exponential Random Graph Models. *Educational and Psychological
#' Measurement, 78*(3), 430-459.
#'
#' Lamprianou, I. (2025). Network Analysis for the investigation of rater
#' effects in language assessment: A comparison of ChatGPT vs human raters.
#' *Research Methods in Applied Linguistics, 4*, 100205.
#'
#' @seealso [interrater_agreement_table()], [plot_interrater_agreement()],
#'   [mfrm_network_analysis()], [plot.mfrm_bundle()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   rn <- rater_network_analysis(fit, mode = "severity_direction")
#'   rn$summary
#'   head(rn$node_metrics)
#'   plot(rn, type = "severity", draw = FALSE)
#' }
#' }
#' @export
rater_network_analysis <- function(fit,
                                   diagnostics = NULL,
                                   rater_facet = NULL,
                                   context_facets = NULL,
                                   mode = c("agreement", "disagreement", "severity_direction"),
                                   weight_metric = NULL,
                                   min_pair_n = 1,
                                   min_weight = 0,
                                   score_diff_tolerance = 0,
                                   severity_continuity = 0.5,
                                   exact_warn = 0.50,
                                   corr_warn = 0.30,
                                   include_graph = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("`rater_network_analysis()` requires the `igraph` package ",
         "(in Suggests). Install it and retry.", call. = FALSE)
  }
  mode <- match.arg(tolower(as.character(mode[1])),
                    c("agreement", "disagreement", "severity_direction"))
  min_pair_n <- max(1L, as.integer(min_pair_n[1]))
  min_weight <- max(0, suppressWarnings(as.numeric(min_weight[1])))
  score_diff_tolerance <- max(0, suppressWarnings(as.numeric(score_diff_tolerance[1])))
  severity_continuity <- max(0, suppressWarnings(as.numeric(severity_continuity[1])))

  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0L) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.", call. = FALSE)
  }

  known_facets <- c("Person", fit$config$facet_names)
  if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
    if (!is.null(diagnostics$interrater$summary) &&
        nrow(diagnostics$interrater$summary) > 0L &&
        "RaterFacet" %in% names(diagnostics$interrater$summary)) {
      rater_facet <- as.character(diagnostics$interrater$summary$RaterFacet[1])
    } else {
      rater_facet <- infer_default_rater_facet(fit$config$facet_names)
    }
  } else {
    rater_facet <- as.character(rater_facet[1])
  }
  if (is.null(rater_facet) || !rater_facet %in% known_facets) {
    stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "),
         call. = FALSE)
  }
  if (identical(rater_facet, "Person")) {
    stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.",
         call. = FALSE)
  }

  if (is.null(context_facets)) {
    facet_cols <- known_facets
  } else {
    context_facets <- unique(as.character(context_facets))
    unknown <- setdiff(context_facets, known_facets)
    if (length(unknown) > 0L) {
      stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    context_facets <- setdiff(context_facets, rater_facet)
    if (length(context_facets) == 0L) {
      stop("`context_facets` must include at least one facet different from `rater_facet`.",
           call. = FALSE)
    }
    facet_cols <- c(rater_facet, context_facets)
  }

  source_interrater <- interrater_agreement_table(
    fit = fit,
    diagnostics = diagnostics,
    rater_facet = rater_facet,
    context_facets = setdiff(facet_cols, rater_facet),
    exact_warn = exact_warn,
    corr_warn = corr_warn,
    top_n = NULL
  )
  pair_metrics <- as.data.frame(source_interrater$pairs %||% data.frame(),
                                stringsAsFactors = FALSE)
  wide_info <- rater_network_score_wide(diagnostics$obs, facet_cols, rater_facet)
  direction_pairs <- rater_network_direction_pairs(
    wide = wide_info$wide,
    raters = wide_info$raters,
    score_diff_tolerance = score_diff_tolerance
  )
  settings <- list(
    rater_facet = rater_facet,
    context_facets = setdiff(facet_cols, rater_facet),
    mode = mode,
    weight_metric = weight_metric,
    min_pair_n = min_pair_n,
    min_weight = min_weight,
    score_diff_tolerance = score_diff_tolerance,
    severity_continuity = severity_continuity,
    exact_warn = exact_warn,
    corr_warn = corr_warn,
    include_graph = isTRUE(include_graph),
    directed = identical(mode, "severity_direction"),
    edge_definition = if (identical(mode, "severity_direction")) {
      "directed edge from relatively higher-scoring/lenient rater to lower-scoring/severe rater"
    } else {
      paste0(mode, " edge between rater pairs in shared scoring contexts")
    }
  )
  if (nrow(pair_metrics) == 0L || length(wide_info$raters) < 2L) {
    return(empty_rater_network_bundle(
      settings = settings,
      source_interrater = source_interrater,
      message = "Fewer than two raters share scorable contexts."
    ))
  }

  pair_metrics <- pair_metrics |>
    dplyr::left_join(direction_pairs, by = c("Rater1", "Rater2")) |>
    dplyr::mutate(
      OneMinusExact = ifelse(is.finite(.data$Exact), 1 - .data$Exact, NA_real_),
      AbsMeanDiff = abs(.data$MeanDiff),
      EligiblePair = is.finite(.data$N) & .data$N >= min_pair_n
    ) |>
    dplyr::arrange(dplyr::desc(.data$EligiblePair), dplyr::desc(.data$MAD),
                   .data$Rater1, .data$Rater2) |>
    as.data.frame(stringsAsFactors = FALSE)

  if (is.null(weight_metric) || !nzchar(as.character(weight_metric[1]))) {
    weight_metric <- if (identical(mode, "agreement")) "Exact" else if (identical(mode, "disagreement")) "MAD" else "DirectionalHigherProp"
  } else {
    weight_metric <- as.character(weight_metric[1])
  }
  settings$weight_metric <- weight_metric

  eligible_pairs <- pair_metrics[pair_metrics$EligiblePair, , drop = FALSE]
  if (identical(mode, "severity_direction")) {
    fwd <- data.frame(
      From = as.character(eligible_pairs$Rater1),
      To = as.character(eligible_pairs$Rater2),
      Pair = paste(eligible_pairs$Rater1, eligible_pairs$Rater2, sep = " | "),
      Weight = suppressWarnings(as.numeric(eligible_pairs$Rater1HigherProp)),
      Count = suppressWarnings(as.numeric(eligible_pairs$Rater1HigherCount)),
      OpportunityCount = suppressWarnings(as.numeric(eligible_pairs$DirectionN)),
      Direction = "Rater1Higher",
      WeightMetric = "DirectionalHigherProp",
      stringsAsFactors = FALSE
    )
    rev <- data.frame(
      From = as.character(eligible_pairs$Rater2),
      To = as.character(eligible_pairs$Rater1),
      Pair = paste(eligible_pairs$Rater1, eligible_pairs$Rater2, sep = " | "),
      Weight = suppressWarnings(as.numeric(eligible_pairs$Rater2HigherProp)),
      Count = suppressWarnings(as.numeric(eligible_pairs$Rater2HigherCount)),
      OpportunityCount = suppressWarnings(as.numeric(eligible_pairs$DirectionN)),
      Direction = "Rater2Higher",
      WeightMetric = "DirectionalHigherProp",
      stringsAsFactors = FALSE
    )
    edges <- dplyr::bind_rows(fwd, rev) |>
      dplyr::filter(is.finite(.data$Weight), .data$Weight >= min_weight,
                    is.finite(.data$Count), .data$Count > 0)
  } else {
    if (!weight_metric %in% names(eligible_pairs)) {
      valid_cols <- names(eligible_pairs)[vapply(eligible_pairs, is.numeric, logical(1))]
      stop("`weight_metric` must be a numeric pair_metrics column: ",
           paste(valid_cols, collapse = ", "), call. = FALSE)
    }
    signed_weight <- suppressWarnings(as.numeric(eligible_pairs[[weight_metric]]))
    graph_weight <- signed_weight
    if (identical(mode, "agreement")) {
      graph_weight <- pmax(graph_weight, 0)
    }
    edges <- data.frame(
      From = as.character(eligible_pairs$Rater1),
      To = as.character(eligible_pairs$Rater2),
      Pair = paste(eligible_pairs$Rater1, eligible_pairs$Rater2, sep = " | "),
      Weight = graph_weight,
      SignedWeight = signed_weight,
      Count = suppressWarnings(as.numeric(eligible_pairs$N)),
      OpportunityCount = suppressWarnings(as.numeric(eligible_pairs$N)),
      Direction = "undirected",
      WeightMetric = weight_metric,
      stringsAsFactors = FALSE
    ) |>
      dplyr::filter(is.finite(.data$Weight), .data$Weight >= min_weight)
  }
  edges <- edges |>
    dplyr::mutate(
      DistanceWeight = 1 / pmax(.data$Weight, .Machine$double.eps),
      EdgeId = dplyr::row_number()
    ) |>
    dplyr::arrange(dplyr::desc(.data$Weight), .data$From, .data$To) |>
    as.data.frame(stringsAsFactors = FALSE)

  vertices <- data.frame(
    name = sort(unique(c(as.character(wide_info$raters),
                         as.character(pair_metrics$Rater1),
                         as.character(pair_metrics$Rater2)))),
    Rater = sort(unique(c(as.character(wide_info$raters),
                          as.character(pair_metrics$Rater1),
                          as.character(pair_metrics$Rater2)))),
    stringsAsFactors = FALSE
  )
  graph_edges <- edges |>
    dplyr::select("From", "To", dplyr::everything())
  directed <- identical(mode, "severity_direction")
  graph <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = directed,
    vertices = vertices
  )
  edge_weights <- suppressWarnings(as.numeric(igraph::E(graph)$Weight))
  if (length(edge_weights) == 0L) edge_weights <- NULL
  edge_dist <- suppressWarnings(as.numeric(igraph::E(graph)$DistanceWeight))
  if (length(edge_dist) == 0L) edge_dist <- NULL
  if (length(edge_dist) > 0L) {
    edge_dist[!is.finite(edge_dist) | edge_dist <= 0] <- 1
  }

  degree_all <- igraph::degree(graph, mode = "all", loops = FALSE)
  degree_in <- igraph::degree(graph, mode = "in", loops = FALSE)
  degree_out <- igraph::degree(graph, mode = "out", loops = FALSE)
  strength_all <- igraph::strength(graph, mode = "all", weights = edge_weights)
  strength_in <- igraph::strength(graph, mode = "in", weights = edge_weights)
  strength_out <- igraph::strength(graph, mode = "out", weights = edge_weights)
  betweenness <- igraph::betweenness(
    graph,
    directed = directed,
    weights = edge_dist,
    normalized = igraph::vcount(graph) > 2L
  )
  closeness <- suppressWarnings(igraph::closeness(
    graph,
    mode = if (directed) "out" else "all",
    weights = edge_dist,
    normalized = TRUE
  ))
  severity_ratio_raw <- suppressWarnings(as.numeric(strength_out) / as.numeric(strength_in))
  severity_ratio <- (as.numeric(strength_out) + severity_continuity) /
    (as.numeric(strength_in) + severity_continuity)
  severity_index <- if (directed) -log(severity_ratio) else rep(NA_real_, length(severity_ratio))

  node_metrics <- data.frame(
    Rater = igraph::V(graph)$name,
    Degree = as.numeric(degree_all),
    InDegree = as.numeric(degree_in),
    OutDegree = as.numeric(degree_out),
    Strength = as.numeric(strength_all),
    InStrength = as.numeric(strength_in),
    OutStrength = as.numeric(strength_out),
    Betweenness = as.numeric(betweenness),
    Closeness = as.numeric(closeness),
    SeverityRatioRaw = if (directed) severity_ratio_raw else NA_real_,
    SeverityRatio = if (directed) severity_ratio else NA_real_,
    SeverityIndex = severity_index,
    RelativePattern = if (directed) {
      dplyr::case_when(
        is.finite(severity_index) & severity_index > 0 ~ "more_severe",
        is.finite(severity_index) & severity_index < 0 ~ "more_lenient",
        is.finite(severity_index) ~ "balanced",
        TRUE ~ "insufficient_directional_edges"
      )
    } else {
      rep(NA_character_, length(severity_index))
    },
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(
      if (directed) dplyr::desc(.data$SeverityIndex) else dplyr::desc(.data$Strength),
      dplyr::desc(.data$Betweenness),
      .data$Rater
    ) |>
    as.data.frame(stringsAsFactors = FALSE)

  edge_metrics <- igraph::as_data_frame(graph, what = "edges")
  if (nrow(edge_metrics) > 0L) {
    names(edge_metrics)[names(edge_metrics) == "from"] <- "From"
    names(edge_metrics)[names(edge_metrics) == "to"] <- "To"
    edge_metrics$Weight <- suppressWarnings(as.numeric(edge_metrics$Weight))
    edge_metrics$DistanceWeight <- suppressWarnings(as.numeric(edge_metrics$DistanceWeight))
    edge_metrics$EdgeBetweenness <- if (igraph::ecount(graph) > 0L) {
      as.numeric(igraph::edge_betweenness(graph, directed = directed, weights = edge_dist))
    } else {
      numeric(0)
    }
    edge_metrics <- edge_metrics |>
      dplyr::arrange(dplyr::desc(.data$Weight), .data$From, .data$To) |>
      as.data.frame(stringsAsFactors = FALSE)
  } else {
    edge_metrics <- data.frame(
      From = character(),
      To = character(),
      Weight = numeric(),
      DistanceWeight = numeric(),
      EdgeBetweenness = numeric(),
      stringsAsFactors = FALSE
    )
  }

  comp <- igraph::components(graph, mode = if (directed) "weak" else "strong")
  density <- igraph::edge_density(graph, loops = FALSE)
  mean_dist <- tryCatch(
    igraph::mean_distance(graph, directed = directed, weights = edge_dist, unconnected = TRUE),
    error = function(e) NA_real_
  )
  diameter <- tryCatch(
    igraph::diameter(graph, directed = directed, weights = edge_dist, unconnected = TRUE),
    error = function(e) NA_real_
  )
  if (!is.finite(mean_dist)) mean_dist <- NA_real_
  if (!is.finite(diameter)) diameter <- NA_real_
  summary_tbl <- data.frame(
    RaterFacet = rater_facet,
    Mode = mode,
    Raters = igraph::vcount(graph),
    PairRows = nrow(eligible_pairs),
    Edges = igraph::ecount(graph),
    Directed = directed,
    WeightMetric = weight_metric,
    Density = density,
    MeanWeight = if (nrow(edge_metrics) > 0L) mean(edge_metrics$Weight, na.rm = TRUE) else NA_real_,
    MeanDegree = mean(node_metrics$Degree, na.rm = TRUE),
    MeanStrength = mean(node_metrics$Strength, na.rm = TRUE),
    Components = as.integer(comp$no),
    Diameter = as.numeric(diameter),
    MeanDistance = as.numeric(mean_dist),
    MeanSeverityIndex = if (directed) mean(node_metrics$SeverityIndex, na.rm = TRUE) else NA_real_,
    SeverityContinuity = severity_continuity,
    ScoreDiffTolerance = score_diff_tolerance,
    MinPairN = min_pair_n,
    MinWeight = min_weight,
    stringsAsFactors = FALSE
  )

  caveats <- data.frame(
    Area = c("scale", "evidence"),
    Severity = c("high", "review"),
    Message = c(
      "Network indices are not MFRM logit estimates and should be compared only as descriptive network diagnostics.",
      "Edges summarize observed score comparisons within shared contexts; they do not replace MFRM fit, bias, or fair-average inference."
    ),
    stringsAsFactors = FALSE
  )
  dropped_pairs <- sum(!pair_metrics$EligiblePair, na.rm = TRUE)
  if (dropped_pairs > 0L) {
    caveats <- rbind(caveats, data.frame(
      Area = "sparse_pairs",
      Severity = "review",
      Message = paste0(dropped_pairs, " rater pair(s) had fewer than min_pair_n shared contexts and were excluded from graph edges."),
      stringsAsFactors = FALSE
    ))
  }
  if (nrow(edge_metrics) == 0L) {
    caveats <- rbind(caveats, data.frame(
      Area = "empty_edges",
      Severity = "high",
      Message = "No edges remained after pair-count and edge-weight thresholds.",
      stringsAsFactors = FALSE
    ))
  }
  if (identical(mode, "agreement") && identical(weight_metric, "Corr") &&
      any(is.finite(pair_metrics$Corr) & pair_metrics$Corr < 0, na.rm = TRUE)) {
    caveats <- rbind(caveats, data.frame(
      Area = "signed_weights",
      Severity = "review",
      Message = "Negative correlations are retained in SignedWeight but truncated to zero for graph-weight centrality.",
      stringsAsFactors = FALSE
    ))
  }

  out <- list(
    summary = summary_tbl,
    node_metrics = node_metrics,
    edge_metrics = edge_metrics,
    pair_metrics = pair_metrics,
    caveats = caveats,
    source_interrater = source_interrater,
    settings = settings
  )
  if (isTRUE(include_graph)) {
    out$graph <- graph
  }
  as_mfrm_bundle(out, "mfrm_rater_network")
}

infer_default_criterion_facet <- function(facet_names, rater_facet = NULL) {
  candidates <- setdiff(as.character(facet_names), as.character(rater_facet %||% character()))
  if (length(candidates) == 0L) return(NULL)
  lower <- tolower(candidates)
  preferred <- candidates[grepl("criterion|criteria|rubric|domain|dimension", lower)]
  if (length(preferred) > 0L) return(preferred[1])
  item_like <- candidates[grepl("item|task|prompt|occasion|category", lower)]
  if (length(item_like) > 0L) return(item_like[1])
  candidates[1]
}

halo_network_wide_scores <- function(obs_df, context_cols, rater_facet, criterion_facet) {
  if (is.null(obs_df) || nrow(obs_df) == 0L || length(context_cols) == 0L) {
    return(list(wide = data.frame(), nodes = data.frame()))
  }
  df <- obs_df |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c(context_cols, rater_facet, criterion_facet)), as.character)
    ) |>
    tidyr::unite(".context", dplyr::all_of(context_cols), sep = "|", remove = FALSE) |>
    dplyr::mutate(
      .node = paste(.data[[rater_facet]], .data[[criterion_facet]], sep = "::")
    ) |>
    dplyr::select(".context", ".node", dplyr::all_of(rater_facet),
                  dplyr::all_of(criterion_facet), "Observed", dplyr::any_of("Weight"))
  df$.Weight <- get_weights(df)
  node_tbl <- df |>
    dplyr::distinct(.data$.node, .data[[rater_facet]], .data[[criterion_facet]]) |>
    dplyr::rename(
      Node = ".node",
      Rater = dplyr::all_of(rater_facet),
      Criterion = dplyr::all_of(criterion_facet)
    ) |>
    dplyr::arrange(.data$Rater, .data$Criterion) |>
    as.data.frame(stringsAsFactors = FALSE)
  scores <- df |>
    dplyr::group_by(.data$.context, .data$.node) |>
    dplyr::summarise(Score = weighted_mean(.data$Observed, .data$.Weight), .groups = "drop")
  wide <- tryCatch(
    tidyr::pivot_wider(
      scores,
      id_cols = ".context",
      names_from = ".node",
      values_from = "Score"
    ),
    error = function(e) NULL
  )
  if (is.null(wide)) wide <- data.frame()
  list(
    wide = as.data.frame(wide, stringsAsFactors = FALSE),
    nodes = node_tbl
  )
}

halo_pair_correlations <- function(wide, node_tbl, method = "spearman",
                                   min_pair_n = 5, p_adjust = "bonferroni") {
  if (is.null(wide) || nrow(wide) == 0L || nrow(node_tbl) < 2L) {
    return(data.frame())
  }
  node_names <- intersect(as.character(node_tbl$Node), names(wide))
  if (length(node_names) < 2L) return(data.frame())
  method <- match.arg(method, c("spearman", "pearson", "kendall"))
  if (!p_adjust %in% stats::p.adjust.methods) {
    stop("`p_adjust` must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  }
  min_pair_n <- max(2L, as.integer(min_pair_n[1]))
  meta <- node_tbl[match(node_names, node_tbl$Node), , drop = FALSE]
  pairs <- utils::combn(seq_along(node_names), 2, simplify = FALSE)
  rows <- lapply(pairs, function(idx) {
    n1 <- node_names[idx[1]]
    n2 <- node_names[idx[2]]
    v1 <- suppressWarnings(as.numeric(wide[[n1]]))
    v2 <- suppressWarnings(as.numeric(wide[[n2]]))
    ok <- is.finite(v1) & is.finite(v2)
    n_ok <- sum(ok)
    r1 <- meta$Rater[idx[1]]
    r2 <- meta$Rater[idx[2]]
    c1 <- meta$Criterion[idx[1]]
    c2 <- meta$Criterion[idx[2]]
    if (n_ok < min_pair_n || length(unique(v1[ok])) < 2L || length(unique(v2[ok])) < 2L) {
      return(data.frame(
        From = n1,
        To = n2,
        Rater1 = r1,
        Criterion1 = c1,
        Rater2 = r2,
        Criterion2 = c2,
        N = n_ok,
        Estimate = NA_real_,
        AbsEstimate = NA_real_,
        PValue = NA_real_,
        EdgeType = if (identical(r1, r2)) "halo" else "non_halo",
        stringsAsFactors = FALSE
      ))
    }
    ct <- suppressWarnings(tryCatch(
      stats::cor.test(v1[ok], v2[ok], method = method, exact = FALSE),
      error = function(e) NULL
    ))
    estimate <- if (is.null(ct)) {
      suppressWarnings(stats::cor(v1[ok], v2[ok], method = method, use = "complete.obs"))
    } else {
      suppressWarnings(as.numeric(ct$estimate[1]))
    }
    p_val <- if (is.null(ct)) NA_real_ else suppressWarnings(as.numeric(ct$p.value[1]))
    data.frame(
      From = n1,
      To = n2,
      Rater1 = r1,
      Criterion1 = c1,
      Rater2 = r2,
      Criterion2 = c2,
      N = n_ok,
      Estimate = estimate,
      AbsEstimate = abs(estimate),
      PValue = p_val,
      EdgeType = if (identical(r1, r2)) "halo" else "non_halo",
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) > 0L) {
    ok <- is.finite(out$PValue)
    out$PAdjusted <- NA_real_
    out$PAdjusted[ok] <- stats::p.adjust(out$PValue[ok], method = p_adjust)
  }
  out
}

empty_halo_network_bundle <- function(settings, message = NULL) {
  summary_tbl <- data.frame(
    RaterFacet = as.character(settings$rater_facet %||% NA_character_),
    CriterionFacet = as.character(settings$criterion_facet %||% NA_character_),
    Nodes = 0L,
    PairRows = 0L,
    Edges = 0L,
    HaloEdges = 0L,
    NonHaloEdges = 0L,
    MeanHaloWeight = NA_real_,
    MeanNonHaloWeight = NA_real_,
    MeanRetainedHaloWeight = NA_real_,
    MeanRetainedNonHaloWeight = NA_real_,
    HaloMinusNonHalo = NA_real_,
    HaloRatio = NA_real_,
    RatersWarning = 0L,
    RatersReview = 0L,
    RatersOk = 0L,
    stringsAsFactors = FALSE
  )
  caveats <- data.frame(
    Area = "halo_network",
    Severity = "high",
    Message = message %||% "No rater-by-criterion halo network could be constructed.",
    stringsAsFactors = FALSE
  )
  as_mfrm_bundle(list(
    summary = summary_tbl,
    node_metrics = data.frame(),
    edge_metrics = data.frame(),
    pair_metrics = data.frame(),
    halo_summary_by_rater = data.frame(),
    caveats = caveats,
    settings = settings
  ), "mfrm_halo_network")
}

#' Analyze rater-by-criterion halo-effect networks
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param rater_facet Name of the rater-like facet.
#' @param criterion_facet Name of the criterion, rubric, task, or item-like
#'   facet used to form rater-by-criterion nodes.
#' @param context_facets Facets defining rows in the reshaped wide matrix.
#'   Defaults to the person facet plus any facets other than the rater and
#'   criterion facets.
#' @param method Correlation method used for rater-by-criterion node pairs.
#' @param min_pair_n Minimum shared contexts required to estimate a node-pair
#'   relationship.
#' @param alpha Adjusted p-value threshold for retaining edges. Set to `1` to
#'   retain all finite correlations after `min_abs_weight` filtering.
#' @param p_adjust Multiple-comparison adjustment passed to [stats::p.adjust()].
#'   The default `"bonferroni"` follows the conservative screening used in
#'   Lamprianou's halo-network example.
#' @param min_abs_weight Minimum absolute correlation retained as a graph edge.
#' @param halo_weight_review Same-rater cross-criterion mean absolute
#'   correlation at or above which a rater is marked for review.
#' @param halo_contrast_review Minimum difference between a rater's mean halo
#'   edge weight and incident non-halo edge weight for a stronger review flag.
#' @param min_retained_halo_edges Minimum retained halo edges required before a
#'   strong `"warning"` status is assigned.
#' @param positive_only If `TRUE`, negative correlations are kept in
#'   `pair_metrics` but excluded from the graph edge table.
#' @param include_graph If `TRUE`, include the underlying `igraph` object.
#'
#' @details
#' `rater_halo_network_analysis()` reshapes rating data so that each
#' rater-by-criterion combination is a node. Edges are correlations between
#' those node score vectors across shared contexts. Edges connecting two nodes
#' from the same rater but different criteria are labelled `"halo"`; all other
#' retained edges are labelled `"non_halo"`.
#'
#' Per-rater `ReviewStatus` combines same-rater cross-criterion mean weight,
#' incident non-halo comparison weight, and the number of retained halo edges.
#' A `"warning"` means these criteria converge strongly enough to prioritize
#' follow-up; `"review"` means at least one screening criterion is elevated.
#' Neither label is a causal halo diagnosis.
#'
#' The key descriptive comparison is the distribution of halo-edge weights
#' versus non-halo-edge weights. A larger halo-edge distribution is consistent
#' with a halo pattern, but this function deliberately reports it as a
#' screening diagnostic. The included Welch test is descriptive only because
#' edge weights are clustered by rater and node.
#'
#' @return A bundle of class `mfrm_halo_network` containing:
#' \describe{
#'   \item{`summary`}{One-row halo-network summary and halo/non-halo contrast.}
#'   \item{`node_metrics`}{Rater-by-criterion node strength and centrality.}
#'   \item{`edge_metrics`}{Retained graph edges.}
#'   \item{`pair_metrics`}{All estimated node-pair correlations before edge
#'     filtering.}
#'   \item{`halo_summary_by_rater`}{Per-rater summaries of same-rater
#'     criterion-pair edges, including `ReviewStatus` and `ReviewReason`.}
#'   \item{`caveats`}{Interpretation notes.}
#' }
#'
#' @references
#' Lai, E. R., Wolfe, E. W., & Vickers, D. (2015). Differentiation of
#' illusory and true halo in writing scores. *Educational and Psychological
#' Measurement, 75*(1), 102-125.
#'
#' Lamprianou, I. (2025). Network Analysis for the investigation of rater
#' effects in language assessment: A comparison of ChatGPT vs human raters.
#' *Research Methods in Applied Linguistics, 4*, 100205.
#'
#' @seealso [rater_network_analysis()], [interrater_agreement_table()],
#'   [plot.mfrm_bundle()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   halo <- rater_halo_network_analysis(fit)
#'   halo$summary
#'   head(halo$halo_summary_by_rater)
#'   plot(halo, type = "edge_distribution", draw = FALSE)
#' }
#' }
#' @export
rater_halo_network_analysis <- function(fit,
                                        diagnostics = NULL,
                                        rater_facet = NULL,
                                        criterion_facet = NULL,
                                        context_facets = NULL,
                                        method = c("spearman", "pearson", "kendall"),
                                        min_pair_n = 5,
                                        alpha = 0.05,
                                        p_adjust = "bonferroni",
                                        min_abs_weight = 0,
                                        halo_weight_review = 0.50,
                                        halo_contrast_review = 0.10,
                                        min_retained_halo_edges = 1,
                                        positive_only = TRUE,
                                        include_graph = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("`rater_halo_network_analysis()` requires the `igraph` package ",
         "(in Suggests). Install it and retry.", call. = FALSE)
  }
  method <- match.arg(tolower(as.character(method[1])), c("spearman", "pearson", "kendall"))
  min_pair_n <- suppressWarnings(as.integer(min_pair_n[1]))
  if (!is.finite(min_pair_n)) min_pair_n <- 5L
  min_pair_n <- max(2L, min_pair_n)
  alpha <- suppressWarnings(as.numeric(alpha[1]))
  if (!is.finite(alpha)) alpha <- 0.05
  alpha <- max(0, min(1, alpha))
  min_abs_weight <- suppressWarnings(as.numeric(min_abs_weight[1]))
  if (!is.finite(min_abs_weight)) min_abs_weight <- 0
  min_abs_weight <- max(0, min_abs_weight)
  halo_weight_review <- suppressWarnings(as.numeric(halo_weight_review[1]))
  if (!is.finite(halo_weight_review)) halo_weight_review <- 0.50
  halo_weight_review <- max(0, halo_weight_review)
  halo_contrast_review <- suppressWarnings(as.numeric(halo_contrast_review[1]))
  if (!is.finite(halo_contrast_review)) halo_contrast_review <- 0.10
  halo_contrast_review <- max(0, halo_contrast_review)
  min_retained_halo_edges <- suppressWarnings(as.integer(min_retained_halo_edges[1]))
  if (!is.finite(min_retained_halo_edges)) min_retained_halo_edges <- 1L
  min_retained_halo_edges <- max(1L, min_retained_halo_edges)
  if (!p_adjust %in% stats::p.adjust.methods) {
    stop("`p_adjust` must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0L) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.", call. = FALSE)
  }

  known_facets <- c("Person", fit$config$facet_names)
  if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
    rater_facet <- infer_default_rater_facet(fit$config$facet_names)
  } else {
    rater_facet <- as.character(rater_facet[1])
  }
  if (is.null(rater_facet) || !rater_facet %in% known_facets || identical(rater_facet, "Person")) {
    stop("`rater_facet` must match a non-person facet: ",
         paste(setdiff(known_facets, "Person"), collapse = ", "), call. = FALSE)
  }
  if (is.null(criterion_facet) || !nzchar(as.character(criterion_facet[1]))) {
    criterion_facet <- infer_default_criterion_facet(fit$config$facet_names, rater_facet)
  } else {
    criterion_facet <- as.character(criterion_facet[1])
  }
  if (is.null(criterion_facet) || !criterion_facet %in% known_facets ||
      identical(criterion_facet, "Person") || identical(criterion_facet, rater_facet)) {
    stop("`criterion_facet` must match a non-person facet different from `rater_facet`.",
         call. = FALSE)
  }
  if (is.null(context_facets)) {
    context_facets <- setdiff(known_facets, c(rater_facet, criterion_facet))
  } else {
    context_facets <- unique(as.character(context_facets))
    unknown <- setdiff(context_facets, known_facets)
    if (length(unknown) > 0L) {
      stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    context_facets <- setdiff(context_facets, c(rater_facet, criterion_facet))
  }
  if (length(context_facets) == 0L) {
    stop("`context_facets` must include at least one facet different from rater and criterion facets.",
         call. = FALSE)
  }

  settings <- list(
    rater_facet = rater_facet,
    criterion_facet = criterion_facet,
    context_facets = context_facets,
    method = method,
    min_pair_n = min_pair_n,
    alpha = alpha,
    p_adjust = p_adjust,
    min_abs_weight = min_abs_weight,
    halo_weight_review = halo_weight_review,
    halo_contrast_review = halo_contrast_review,
    min_retained_halo_edges = min_retained_halo_edges,
    positive_only = isTRUE(positive_only),
    include_graph = isTRUE(include_graph),
    node_definition = "rater-by-criterion score profile",
    halo_edge_definition = "edge connecting two criteria scored by the same rater"
  )
  wide_info <- halo_network_wide_scores(
    obs_df = diagnostics$obs,
    context_cols = context_facets,
    rater_facet = rater_facet,
    criterion_facet = criterion_facet
  )
  node_tbl <- as.data.frame(wide_info$nodes %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(node_tbl) < 2L) {
    return(empty_halo_network_bundle(settings, "Fewer than two rater-by-criterion nodes are available."))
  }
  pair_metrics <- halo_pair_correlations(
    wide = wide_info$wide,
    node_tbl = node_tbl,
    method = method,
    min_pair_n = min_pair_n,
    p_adjust = p_adjust
  )
  if (nrow(pair_metrics) == 0L) {
    return(empty_halo_network_bundle(settings, "No rater-by-criterion node pairs could be estimated."))
  }
  pair_metrics <- pair_metrics |>
    dplyr::mutate(
      RetainedByN = is.finite(.data$N) & .data$N >= min_pair_n,
      RetainedByP = is.finite(.data$PAdjusted) & .data$PAdjusted <= alpha,
      RetainedByWeight = is.finite(.data$AbsEstimate) & .data$AbsEstimate >= min_abs_weight,
      RetainedBySign = if (isTRUE(positive_only)) is.finite(.data$Estimate) & .data$Estimate > 0 else is.finite(.data$Estimate),
      RetainedEdge = .data$RetainedByN & .data$RetainedByP & .data$RetainedByWeight & .data$RetainedBySign
    ) |>
    dplyr::arrange(dplyr::desc(.data$RetainedEdge), .data$EdgeType,
                   dplyr::desc(.data$AbsEstimate), .data$From, .data$To) |>
    as.data.frame(stringsAsFactors = FALSE)

  edges <- pair_metrics[pair_metrics$RetainedEdge, , drop = FALSE] |>
    dplyr::transmute(
      From = .data$From,
      To = .data$To,
      Rater1 = .data$Rater1,
      Criterion1 = .data$Criterion1,
      Rater2 = .data$Rater2,
      Criterion2 = .data$Criterion2,
      EdgeType = .data$EdgeType,
      Weight = .data$AbsEstimate,
      SignedWeight = .data$Estimate,
      N = .data$N,
      PValue = .data$PValue,
      PAdjusted = .data$PAdjusted,
      DistanceWeight = 1 / pmax(.data$AbsEstimate, .Machine$double.eps)
    ) |>
    dplyr::arrange(.data$EdgeType, dplyr::desc(.data$Weight), .data$From, .data$To) |>
    as.data.frame(stringsAsFactors = FALSE)

  vertices <- node_tbl |>
    dplyr::rename(name = "Node") |>
    dplyr::arrange(.data$Rater, .data$Criterion) |>
    as.data.frame(stringsAsFactors = FALSE)
  graph_edges <- edges |>
    dplyr::select("From", "To", dplyr::everything())
  graph <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = FALSE,
    vertices = vertices
  )
  edge_weights <- suppressWarnings(as.numeric(igraph::E(graph)$Weight))
  if (length(edge_weights) == 0L) edge_weights <- NULL
  edge_dist <- suppressWarnings(as.numeric(igraph::E(graph)$DistanceWeight))
  if (length(edge_dist) == 0L) edge_dist <- NULL
  if (length(edge_dist) > 0L) {
    edge_dist[!is.finite(edge_dist) | edge_dist <= 0] <- 1
  }
  degree <- igraph::degree(graph, mode = "all", loops = FALSE)
  strength <- igraph::strength(graph, mode = "all", weights = edge_weights)
  betweenness <- igraph::betweenness(
    graph,
    directed = FALSE,
    weights = edge_dist,
    normalized = igraph::vcount(graph) > 2L
  )
  closeness <- suppressWarnings(igraph::closeness(
    graph,
    mode = "all",
    weights = edge_dist,
    normalized = TRUE
  ))
  node_names <- igraph::V(graph)$name
  halo_strength <- rep(0, length(node_names))
  non_halo_strength <- rep(0, length(node_names))
  if (nrow(edges) > 0L) {
    for (i in seq_len(nrow(edges))) {
      w <- suppressWarnings(as.numeric(edges$Weight[i]))
      if (!is.finite(w)) next
      idx <- match(c(edges$From[i], edges$To[i]), node_names)
      idx <- idx[is.finite(idx)]
      if (identical(as.character(edges$EdgeType[i]), "halo")) {
        halo_strength[idx] <- halo_strength[idx] + w
      } else {
        non_halo_strength[idx] <- non_halo_strength[idx] + w
      }
    }
  }
  node_metrics <- data.frame(
    Node = node_names,
    Rater = as.character(igraph::V(graph)$Rater),
    Criterion = as.character(igraph::V(graph)$Criterion),
    Degree = as.numeric(degree),
    Strength = as.numeric(strength),
    HaloStrength = halo_strength,
    NonHaloStrength = non_halo_strength,
    HaloStrengthShare = ifelse((halo_strength + non_halo_strength) > 0,
                               halo_strength / (halo_strength + non_halo_strength),
                               NA_real_),
    Betweenness = as.numeric(betweenness),
    Closeness = as.numeric(closeness),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(.data$HaloStrengthShare),
                   dplyr::desc(.data$HaloStrength), .data$Rater, .data$Criterion) |>
    as.data.frame(stringsAsFactors = FALSE)

  edge_metrics <- if (nrow(edges) > 0L) {
    ed <- igraph::as_data_frame(graph, what = "edges")
    names(ed)[names(ed) == "from"] <- "From"
    names(ed)[names(ed) == "to"] <- "To"
    ed$Weight <- suppressWarnings(as.numeric(ed$Weight))
    ed$SignedWeight <- suppressWarnings(as.numeric(ed$SignedWeight))
    ed$DistanceWeight <- suppressWarnings(as.numeric(ed$DistanceWeight))
    ed$EdgeBetweenness <- as.numeric(igraph::edge_betweenness(
      graph,
      directed = FALSE,
      weights = edge_dist
    ))
    ed |>
      dplyr::arrange(.data$EdgeType, dplyr::desc(.data$Weight), .data$From, .data$To) |>
      as.data.frame(stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  retained_weights <- pair_metrics[pair_metrics$RetainedByN & is.finite(pair_metrics$AbsEstimate), , drop = FALSE]
  halo_vals <- retained_weights$AbsEstimate[retained_weights$EdgeType == "halo"]
  non_halo_vals <- retained_weights$AbsEstimate[retained_weights$EdgeType == "non_halo"]
  wt <- if (length(halo_vals) >= 2L && length(non_halo_vals) >= 2L) {
    suppressWarnings(tryCatch(
      stats::t.test(halo_vals, non_halo_vals),
      error = function(e) NULL
    ))
  } else {
    NULL
  }
  review_pair_metrics <- pair_metrics[pair_metrics$RetainedByN &
                                        is.finite(pair_metrics$AbsEstimate), ,
                                      drop = FALSE]
  rater_levels <- sort(unique(as.character(node_tbl$Rater)))
  halo_summary_by_rater <- dplyr::bind_rows(lapply(rater_levels, function(rater) {
    halo_rows <- review_pair_metrics[
      review_pair_metrics$EdgeType == "halo" &
        as.character(review_pair_metrics$Rater1) == rater,
      ,
      drop = FALSE
    ]
    non_halo_rows <- review_pair_metrics[
      review_pair_metrics$EdgeType == "non_halo" &
        (as.character(review_pair_metrics$Rater1) == rater |
           as.character(review_pair_metrics$Rater2) == rater),
      ,
      drop = FALSE
    ]
    retained_halo <- halo_rows[halo_rows$RetainedEdge, , drop = FALSE]
    mean_halo <- if (nrow(halo_rows) > 0L) mean(halo_rows$AbsEstimate, na.rm = TRUE) else NA_real_
    mean_non_halo <- if (nrow(non_halo_rows) > 0L) mean(non_halo_rows$AbsEstimate, na.rm = TRUE) else NA_real_
    mean_retained_halo <- if (nrow(retained_halo) > 0L) mean(retained_halo$AbsEstimate, na.rm = TRUE) else NA_real_
    halo_contrast <- if (is.finite(mean_halo) && is.finite(mean_non_halo)) {
      mean_halo - mean_non_halo
    } else {
      NA_real_
    }
    halo_ratio <- if (is.finite(mean_halo) && is.finite(mean_non_halo) && mean_non_halo > 0) {
      mean_halo / mean_non_halo
    } else {
      NA_real_
    }
    retained_n <- nrow(retained_halo)
    retained_share <- if (nrow(halo_rows) > 0L) retained_n / nrow(halo_rows) else NA_real_
    status <- dplyr::case_when(
      !is.finite(mean_halo) ~ "insufficient_data",
      retained_n >= min_retained_halo_edges &&
        mean_halo >= halo_weight_review &&
        is.finite(halo_contrast) &&
        halo_contrast >= halo_contrast_review ~ "warning",
      retained_n > 0L ||
        mean_halo >= halo_weight_review ||
        (is.finite(halo_contrast) && halo_contrast >= halo_contrast_review) ~ "review",
      TRUE ~ "ok"
    )
    reason <- if (identical(status, "insufficient_data")) {
      "Too few estimable same-rater cross-criterion pairs for halo review."
    } else {
      sprintf(
        "Mean halo weight %.3f; incident non-halo mean %.3f; halo contrast %.3f; retained halo edges %d/%d.",
        mean_halo,
        mean_non_halo,
        halo_contrast,
        retained_n,
        nrow(halo_rows)
      )
    }
    data.frame(
      Rater = rater,
      HaloPairs = nrow(halo_rows),
      RetainedHaloEdges = retained_n,
      RetainedHaloShare = retained_share,
      MeanHaloWeight = mean_halo,
      MeanRetainedHaloWeight = mean_retained_halo,
      MaxHaloWeight = if (nrow(halo_rows) > 0L) max(halo_rows$AbsEstimate, na.rm = TRUE) else NA_real_,
      MeanSignedHaloWeight = if (nrow(halo_rows) > 0L) mean(halo_rows$Estimate, na.rm = TRUE) else NA_real_,
      IncidentNonHaloPairs = nrow(non_halo_rows),
      MeanIncidentNonHaloWeight = mean_non_halo,
      HaloMinusIncidentNonHalo = halo_contrast,
      HaloRatioIncident = halo_ratio,
      ReviewStatus = status,
      ReviewReason = reason,
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::arrange(
      factor(.data$ReviewStatus, levels = c("warning", "review", "ok", "insufficient_data")),
      dplyr::desc(.data$HaloMinusIncidentNonHalo),
      dplyr::desc(.data$MeanHaloWeight),
      .data$Rater
    ) |>
    as.data.frame(stringsAsFactors = FALSE)

  retained_halo_vals <- edge_metrics$Weight[edge_metrics$EdgeType == "halo"]
  retained_non_halo_vals <- edge_metrics$Weight[edge_metrics$EdgeType == "non_halo"]

  summary_tbl <- data.frame(
    RaterFacet = rater_facet,
    CriterionFacet = criterion_facet,
    Nodes = igraph::vcount(graph),
    PairRows = nrow(pair_metrics),
    Edges = igraph::ecount(graph),
    HaloEdges = sum(edge_metrics$EdgeType == "halo", na.rm = TRUE),
    NonHaloEdges = sum(edge_metrics$EdgeType == "non_halo", na.rm = TRUE),
    MeanHaloWeight = if (length(halo_vals) > 0L) mean(halo_vals, na.rm = TRUE) else NA_real_,
    MeanNonHaloWeight = if (length(non_halo_vals) > 0L) mean(non_halo_vals, na.rm = TRUE) else NA_real_,
    MeanRetainedHaloWeight = if (length(retained_halo_vals) > 0L) mean(retained_halo_vals, na.rm = TRUE) else NA_real_,
    MeanRetainedNonHaloWeight = if (length(retained_non_halo_vals) > 0L) mean(retained_non_halo_vals, na.rm = TRUE) else NA_real_,
    HaloMinusNonHalo = if (length(halo_vals) > 0L && length(non_halo_vals) > 0L) {
      mean(halo_vals, na.rm = TRUE) - mean(non_halo_vals, na.rm = TRUE)
    } else {
      NA_real_
    },
    HaloRatio = if (length(halo_vals) > 0L && length(non_halo_vals) > 0L &&
                    is.finite(mean(non_halo_vals, na.rm = TRUE)) &&
                    mean(non_halo_vals, na.rm = TRUE) > 0) {
      mean(halo_vals, na.rm = TRUE) / mean(non_halo_vals, na.rm = TRUE)
    } else {
      NA_real_
    },
    WelchT = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
    WelchDF = if (!is.null(wt)) unname(wt$parameter) else NA_real_,
    WelchP = if (!is.null(wt)) wt$p.value else NA_real_,
    RatersWarning = sum(halo_summary_by_rater$ReviewStatus == "warning", na.rm = TRUE),
    RatersReview = sum(halo_summary_by_rater$ReviewStatus == "review", na.rm = TRUE),
    RatersOk = sum(halo_summary_by_rater$ReviewStatus == "ok", na.rm = TRUE),
    Method = method,
    PAdjust = p_adjust,
    Alpha = alpha,
    MinPairN = min_pair_n,
    MinAbsWeight = min_abs_weight,
    PositiveOnly = isTRUE(positive_only),
    stringsAsFactors = FALSE
  )

  caveats <- data.frame(
    Area = c("construct", "inference"),
    Severity = c("review", "high"),
    Message = c(
      "Halo edges are same-rater cross-criterion correlations; they screen for halo-like score-profile similarity, not causal halo by themselves.",
      "Welch halo/non-halo comparisons are descriptive because network edges are clustered and statistically dependent."
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(edge_metrics) == 0L) {
    caveats <- rbind(caveats, data.frame(
      Area = "empty_edges",
      Severity = "review",
      Message = "No edges remained after adjusted-p and weight filtering; inspect pair_metrics or relax alpha/min_abs_weight for exploratory visualization.",
      stringsAsFactors = FALSE
    ))
  }
  if (any(!pair_metrics$RetainedByN, na.rm = TRUE)) {
    caveats <- rbind(caveats, data.frame(
      Area = "sparse_pairs",
      Severity = "review",
      Message = paste0(sum(!pair_metrics$RetainedByN, na.rm = TRUE),
                       " rater-by-criterion pair(s) had fewer than min_pair_n shared contexts."),
      stringsAsFactors = FALSE
    ))
  }

  out <- list(
    summary = summary_tbl,
    node_metrics = node_metrics,
    edge_metrics = edge_metrics,
    pair_metrics = pair_metrics,
    halo_summary_by_rater = halo_summary_by_rater,
    caveats = caveats,
    settings = settings
  )
  if (isTRUE(include_graph)) {
    out$graph <- graph
  }
  as_mfrm_bundle(out, "mfrm_halo_network")
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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

build_fit_separation_reporting_basis <- function(fit, diagnostics) {
  measures_tbl <- as.data.frame(
    diagnostics$measures %||% data.frame(),
    stringsAsFactors = FALSE
  )
  reliability_tbl <- as.data.frame(
    diagnostics$reliability %||% data.frame(),
    stringsAsFactors = FALSE
  )
  facet_precision_tbl <- as.data.frame(
    diagnostics$facet_precision %||% data.frame(),
    stringsAsFactors = FALSE
  )
  model <- toupper(as.character(fit$summary$Model[1] %||% fit$config$model %||% NA_character_))

  has_mnsq <- nrow(measures_tbl) > 0L &&
    any(c("Infit", "Outfit", "InfitMnSq", "OutfitMnSq") %in% names(measures_tbl))
  has_zstd <- nrow(measures_tbl) > 0L &&
    any(c(
      "InfitZSTD", "OutfitZSTD",
      "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
      "InfitZSTD_FACETS", "OutfitZSTD_FACETS"
    ) %in% names(measures_tbl))
  has_df_review <- nrow(measures_tbl) > 0L &&
    any(c(
      "DF_Infit_ENGINE", "DF_Outfit_ENGINE",
      "DF_Infit_FACETS", "DF_Outfit_FACETS"
    ) %in% names(measures_tbl))
  has_separation <- nrow(reliability_tbl) > 0L ||
    (nrow(facet_precision_tbl) > 0L &&
       any(c(
         "Separation", "RealSeparation", "Reliability",
         "RealReliability", "Strata", "RealStrata"
       )
           %in% names(facet_precision_tbl)))

  availability <- c(
    if (has_mnsq) "available_in_diagnostics" else "not_available_in_diagnostics",
    if (has_zstd && has_df_review) {
      "available_with_df_review"
    } else if (has_zstd) {
      "available_without_df_review"
    } else {
      "not_available_in_diagnostics"
    },
    if (has_separation) "available_in_precision_tables" else "not_available_in_precision_tables",
    if (identical(model, "GPCM")) {
      "restricted_for_gpcm_bundles"
    } else {
      "available_where_qc_pipeline_supports_model"
    }
  )

  data.frame(
    Topic = c(
      "Fit MnSq",
      "Fit ZSTD",
      "Separation reliability and strata",
      "Operational QC thresholds"
    ),
    SourceBasis = c(
      "Wright & Linacre (1994); Linacre (2002)",
      "Linacre (2002); FACETS WHEXACT documentation",
      "Wright & Masters (1982); Wright & Masters (2002); FACETS manual",
      "Package QC policy layered on the fit and separation conventions above"
    ),
    PackageSurface = c(
      "diagnose_mfrm(); fit_measures_table()",
      "diagnose_mfrm(fit_df_method = \"both\"); facets_fit_review()",
      "diagnostics$reliability; facet_statistics_report(); precision_review_report()",
      "run_qc_pipeline(); evaluate_mfrm_design()"
    ),
    Interpretation = c(
      "Mean-square fit is the primary size diagnostic; values are read relative to the expected value of 1.",
      "ZSTD standardizes mean-square fit and is sensitive to df, transformation, and sample-size conventions.",
      "Separation/reliability/strata summarize spread relative to average measurement error; they are not inter-rater agreement.",
      "Pass/warn/fail cutoffs are reporting policy overlays and should remain separate from formula validation."
    ),
    ValidationUse = c(
      "Use as diagnostic evidence and external comparison input; not a standalone validation success criterion.",
      "Compare MnSq first; label df-driven ZSTD changes convention-sensitive when validating against FACETS-style output.",
      "Report with precision tier and model/real basis; do not use separation alone as measurement-quality proof.",
      "Use for operational triage after formula/source checks; calibrate thresholds with simulations or external reference cases."
    ),
    Availability = availability,
    stringsAsFactors = FALSE
  )
}

#' Build a precision review report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#'
#' @details
#' This helper summarizes how `mfrmr` derived SE, CI, and reliability values
#' for the current run. It also includes a source-grounded fit/separation
#' basis table so users can keep mean-square fit, ZSTD standardization,
#' Rasch/FACETS-style separation, and package QC thresholds in separate
#' reporting lanes.
#'
#' @section What this review means:
#' `precision_review_report()` is a reporting gatekeeper for precision claims.
#' It tells you how the package derived uncertainty summaries for the current
#' run and how cautiously those summaries should be written up.
#'
#' @section What this review does not justify:
#' - It does not, by itself, validate the measurement model or substantive
#'   conclusions.
#' - A favorable precision tier does not override convergence, fit, linking,
#'   or design problems elsewhere in the analysis.
#' - Fit and separation rows in this report are reporting/validation
#'   boundaries, not standalone success criteria.
#'
#' @section Interpreting output:
#' - `profile`: one-row overview of the active precision tier and recommended use.
#' - `checks`: package-native review checks for SE ordering, reliability ordering,
#'   coverage of sample/population summaries, and SE source labels.
#' - `fit_separation_basis`: source-grounded boundary table for fit and
#'   separation reporting.
#' - `approximation_notes`: method notes copied from `diagnose_mfrm()`.
#'
#' @section Recommended next step:
#' Use the `profile$PrecisionTier` and `checks` table to decide whether SE, CI,
#' and reliability language can be phrased as model-based, should be qualified
#' as hybrid, or should remain exploratory in the final report.
#'
#' @section Typical workflow:
#' 1. Run `diagnose_mfrm()` for the fitted model.
#' 2. Build `precision_review_report(fit, diagnostics = diag)`.
#' 3. Use `summary()` to see whether the run supports model-based reporting
#'    language or should remain in exploratory/screening mode.
#'
#' @return A named list with:
#' - `profile`: one-row precision overview
#' - `checks`: package-native precision review checks
#' - `fit_separation_basis`: source-grounded fit/separation reporting boundary
#' - `approximation_notes`: detailed method notes
#' - `settings`: resolved model and method labels
#'
#' @seealso [diagnose_mfrm()], [facet_statistics_report()], [reporting_checklist()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' out <- precision_review_report(fit, diagnostics = diag)
#' summary(out)
#' @name precision_review_report
#' @export
precision_review_report <- function(fit, diagnostics = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  profile_tbl <- as.data.frame(diagnostics$precision_profile %||% data.frame(), stringsAsFactors = FALSE)
  checks_tbl <- as.data.frame(precision_review(diagnostics, required = FALSE) %||% data.frame(), stringsAsFactors = FALSE)
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
    fit_separation_basis = build_fit_separation_reporting_basis(fit, diagnostics),
    approximation_notes = notes_tbl,
    settings = settings
  )
  as_mfrm_bundle(out, "mfrm_precision_review")
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
#' The expected-score table also carries the per-curve score variance and
#' information function. For `GPCM`, the information column follows the
#' Muraki/Samejima identity \eqn{a^2 \mathrm{Var}(X \mid \theta)};
#' for `RSM` / `PCM`, this reduces to the usual score variance because
#' discrimination is fixed at one. The `category_information` table decomposes
#' that total into category-level contributions,
#' \eqn{a^2 P_k(\theta)(k - E[X \mid \theta])^2}, whose sum equals the
#' reported information at the same theta value. The
#' `cumulative_probabilities` table follows the FACETS / Winsteps graph
#' convention of accumulating modeled probabilities across ordered categories
#' (`P(X <= k)` by default, with `P(X >= k)` also returned for flipped curves).
#' `cumulative_boundaries` reports approximate theta values where
#' `P(X <= k) = .5`, with `BoundaryStatus` and `CrossingCount` to avoid
#' over-interpreting boundaries outside the requested theta range or with
#' multiple crossings.
#'
#' @section Interpreting output:
#' Use this report to inspect:
#' - where each category has highest probability across theta
#' - where cumulative category probabilities cross .5
#' - whether adjacent categories cross in expected order
#' - whether probability bands look compressed (often sparse categories)
#'
#' Recommended read order:
#' 1. `summary(out)` for compact diagnostics.
#' 2. `out$probabilities`, `out$expected_ogive`, and
#'    `out$category_information` for custom graphics.
#' 3. `plot(out)` for a default visual check, or
#'    `plot(out, type = "cumulative")` to inspect cumulative probabilities.
#'    `plot(out, type = "information")` to inspect curve-level information.
#'    Use `plot(out, type = "category_information")` when category-level
#'    contributions are needed.
#'
#' @section References:
#' Category response curves follow Andrich's rating-scale formulation,
#' Masters' partial-credit model, and Muraki's generalized partial-credit
#' model. The `Information` column for bounded `GPCM` uses Muraki's
#' item-information result obtained from Samejima's general polytomous
#' information formula.
#'
#' - Andrich, D. (1978). *A rating formulation for ordered response
#'   categories*. Psychometrika, 43(4), 561-573.
#' - Masters, G. N. (1982). *A Rasch model for partial credit scoring*.
#'   Psychometrika, 47(2), 149-174.
#' - Muraki, E. (1992). *A generalized partial credit model: Application
#'   of an EM algorithm*. Applied Psychological Measurement, 16(2),
#'   159-176. \doi{10.1177/014662169201600206}
#' - Muraki, E. (1993). *Information functions of the generalized
#'   partial credit model*. Applied Psychological Measurement, 17(4),
#'   351-363. \doi{10.1177/014662169301700403}
#'
#' @section Typical workflow:
#' 1. Fit model with [fit_mfrm()].
#' 2. Run `category_curves_report()` with suitable `theta_points`.
#' 3. Use `summary()` and `plot()`; export tables for manuscripts/dashboard use.
#'    `plot(out)` gives a four-panel overview. Use
#'    `preset = "monochrome"` for grayscale/line-type output and
#'    `boundary_status = "none"` when cumulative `.5` boundary lines should
#'    be suppressed. `plot(out, type = "category_probability")` and
#'    `plot(out, type = "conditional_probability")` are explicit aliases for
#'    the same category-probability curves as `type = "ccc"`. Use
#'    `plot_data(out, component = "plot_long")` when rebuilding the curves with
#'    ggplot2, plotly, or another R graphics system.
#' @return A named list with category-curve components. Class:
#'   `mfrm_category_curves`.
#' @seealso [category_structure_report()], [rating_scale_table()], [plot.mfrm_fit()],
#'   [mfrmr_reports_and_tables], [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' out <- category_curves_report(fit, theta_points = 101)
#' summary(out)
#' head(out$probabilities[, c("CurveGroup", "Theta", "Category", "Probability")])
#' p_overview <- plot(out, draw = FALSE)
#' p_overview$data$plot
#' p_cum <- plot(out, type = "cumulative", draw = FALSE)
#' head(p_cum$data$cumulative_boundaries)
#' p_info <- plot(out, type = "category_information", draw = FALSE)
#' head(p_info$data$category_information)
#' curve_long <- plot_data(out, component = "plot_long")
#' head(curve_long[, c("PlotType", "Theta", "Series", "Value")])
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

#' Build a bias-interaction plot-data bundle (FACETS Table 13: ranked bias list)
#'
#' Bundles the **ranked flagged-cells** view of a bias-interaction run for
#' downstream printing and plotting. The three sibling reports in this
#' family are intentionally distinct:
#' - [bias_interaction_report()] (this one) = FACETS Table 13: a ranked
#'   list of interaction cells with `t`, `bias size`, and screening tail
#'   area -- use when reviewing which `(facet_a, facet_b)` cells deserve
#'   follow-up.
#' - [bias_iteration_report()] = iteration history / convergence trace
#'   for the bias recalibration (FACETS Table 9 territory) -- use when
#'   diagnosing whether the bias run itself stabilised.
#' - [bias_pairwise_report()] = pairwise contrast table for a target
#'   facet (FACETS Table 14 territory) -- use when comparing levels
#'   within a facet while controlling for the other.
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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

#' Build a bias-iteration report (FACETS Table 9: iteration / convergence trace)
#'
#' This report is NOT an alias of [bias_interaction_report()] despite the
#' similar name. It focuses on the **recalibration path** of a bias run:
#' iteration table, convergence summary, and orientation review. Use this
#' to confirm that the bias recalibration itself converged; use
#' [bias_interaction_report()] to review the ranked flagged cells from
#' the converged run.
#'
#' @inheritParams bias_interaction_report
#' @param top_n Maximum number of iteration rows to keep in preview-oriented
#'   summaries. The full iteration table is always returned.
#'
#' @details
#' This report focuses on the recalibration path used by [estimate_bias()].
#' It provides a package-native counterpart to legacy iteration printouts by
#' exposing the iteration table, convergence summary, and orientation review in
#' one bundle.
#'
#' @return A named list with:
#' - `table`: iteration history
#' - `summary`: one-row convergence summary
#' - `orientation_review`: interaction-facet sign review
#' - `settings`: resolved reporting options
#' - `direction_note`: one-line interpretive note describing which
#'   direction the iteration moved (carried from the bias estimator;
#'   empty string when the underlying estimator does not emit one)
#' - `recommended_action`: one-line recommended action label
#'   (e.g. `"converged"`, `"increase max_iter"`); empty string when
#'   the underlying estimator does not emit one
#'
#' @seealso [estimate_bias()], [bias_interaction_report()], [build_fixed_reports()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
    orientation_review = as.data.frame(bias_results$orientation_review %||% data.frame(), stringsAsFactors = FALSE),
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

#' Build a bias pairwise-contrast report (FACETS Table 14: pairwise contrasts)
#'
#' Build a pairwise contrast table that, for a chosen target facet
#' (e.g. raters), compares each pair of target-facet levels while
#' holding a context facet (e.g. items / criteria) constant. This is
#' the FACETS Table 14 view: it answers "is rater A consistently
#' more severe than rater B on the same items?" rather than "which
#' (rater, item) cell has the largest local bias?" -- the latter is
#' covered by [bias_interaction_report()].
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
#' @section Interpreting output:
#' - `table`: one row per ordered (target_level_1, target_level_2)
#'   pair, with `Bias_diff`, `SE_diff`, `t_diff`, `df_diff`,
#'   `p_diff`, and the underlying per-level bias rows. Rows are
#'   sorted so that the largest-magnitude `|t_diff|` rises to the
#'   top.
#' - `summary`: one-row screening summary with `MaxAbsBiasDiff`,
#'   `MaxAbsT`, `Significant` (count of flagged pairs at `p_max`),
#'   `BonferroniSignificant`, and `HolmSignificant`.
#' - `orientation_review` carries the same facet-orientation sign
#'   review as the parent `estimate_bias()` run.
#' - The SE caveat below applies: read `Significant` /
#'   `BonferroniSignificant` as a screening triage, not as formal
#'   inferential tests.
#'
#' @section Typical workflow:
#' 1. Fit and diagnose the model.
#' 2. Run `estimate_bias()` to get the underlying interaction effects.
#' 3. Pass that result to `bias_pairwise_report()` for the rater-pair
#'    contrast table.
#' 4. Use `summary(out)$MaxAbsT` and the top rows of `out$table` to
#'    flag rater-pair systematic differences for follow-up review.
#' 5. For the ranked flagged-cells view (which (rater, item) pairs
#'    have the largest local bias), use `bias_interaction_report()`
#'    on the same `estimate_bias()` output.
#'
#' @section Standard-error caveat:
#' The contrast standard error is computed as
#' `SE(b_i - b_j) = sqrt(SE_i^2 + SE_j^2)` -- the independence
#' approximation. For same-facet bias values that share a sum-to-zero
#' identification, `Cov(b_i, b_j) < 0`, so the true contrast variance
#' is `SE_i^2 + SE_j^2 - 2 * Cov(b_i, b_j)`, which is **smaller**
#' than the reported value. The reported t-statistics and p-values
#' are therefore conservative for same-facet contrasts (the true
#' significance is higher than reported). For across-facet contrasts
#' the covariance term is approximately zero and the approximation
#' is appropriate. Use the report as a screening / triage table; for
#' inferential claims that hinge on a marginally-significant
#' same-facet contrast, follow up with a contrast that uses the full
#' parameter covariance.
#'
#' @return A named list with:
#' - `table`: pairwise contrast rows
#' - `summary`: one-row contrast summary
#' - `orientation_review`: interaction-facet sign review
#' - `settings`: resolved reporting options
#' - `direction_note`: one-line interpretive note describing the
#'   dominant pairwise-contrast direction (carried from the
#'   underlying bias estimator; empty string when not applicable)
#' - `recommended_action`: one-line recommended-action label
#'   (e.g. routing the user to follow-up review of the largest
#'   flagged pairs); empty string when the underlying estimator
#'   does not emit one
#'
#' @section References:
#' - Linacre, J. M. (1989). *Many-Facet Rasch Measurement*. MESA Press.
#' - Eckes, T. (2005). Examining rater effects in TestDaF writing and
#'   speaking performance assessments: A many-facet Rasch analysis.
#'   *Language Assessment Quarterly, 2*(3), 197-221.
#' - Myford, C. M., & Wolfe, E. W. (2003). Detecting and measuring
#'   rater effects using many-facet Rasch measurement: Part I.
#'   *Journal of Applied Measurement, 4*(4), 386-422.
#' - Myford, C. M., & Wolfe, E. W. (2004). Detecting and measuring
#'   rater effects using many-facet Rasch measurement: Part II.
#'   *Journal of Applied Measurement, 5*(2), 189-227.
#'
#' @seealso [estimate_bias()], [bias_interaction_report()], [build_fixed_reports()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' out <- bias_pairwise_report(fit, diagnostics = diag, facet_a = "Rater", facet_b = "Criterion")
#' s <- summary(out)
#' s$summary
#' # Look for: `MaxAbsBiasDiff` < ~0.5 logits and `Significant = 0` mean
#' #   no rater pair contrasts above the screen. The `BonferroniSignificant`
#' #   / `HolmSignificant` columns count pairs that survive multiple-
#' #   testing correction; both being 0 is a stronger "no rater-pair
#' #   inconsistency" signal than the raw screen-positive count alone.
#' head(out$table)
#' # Look for: top rows with `|t_diff|` > 2 and |Bias_diff| > 0.5 logits
#' #   warrant content-review of the two raters' scoring conventions on
#' #   the conditioning context facet (e.g. compare their item-level
#' #   marks for systematic strictness/leniency patterns).
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
    orientation_review = as.data.frame(bias_results$orientation_review %||% data.frame(), stringsAsFactors = FALSE),
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
#' @param plot Plot type: `"scatter"`, `"ranked"`, `"heatmap"`,
#'   `"abs_t_hist"`, or `"facet_profile"`.
#' @param show_ci Logical. When `TRUE` and `plot` is `"scatter"` or
#'   `"ranked"`, draw confidence-interval whiskers for `Bias Size`.
#'   Bounded `GPCM` rows use the conditional profile-likelihood limits
#'   returned by [estimate_bias()] when available; otherwise the interval
#'   uses the per-cell standard error from [estimate_bias()]. Ignored for
#'   `"heatmap"`, `"abs_t_hist"`, and `"facet_profile"`.
#' @param ci_level Confidence level used when `show_ci = TRUE`; default
#'   `0.95`. The returned plot-data object gains `CI_Lower` / `CI_Upper`
#'   / `CI_Level` columns on the `ranked_table` and `scatter_data`
#'   elements for downstream reuse.
#' @param main Optional plot title override.
#' @param palette Optional named color overrides (`normal`, `flag`, `hist`,
#'   `profile`).
#' @param label_angle Label angle hint for ranked/profile labels.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' Visualization front-end for [bias_interaction_report()] with multiple views.
#' With `draw = FALSE`, the returned plot data include `plot_long`,
#' `plot_annotations`, `flag_summary`, and `plot_settings` in addition to the
#' view-specific `ranked_table`, `scatter_data`, `facet_profile`, and heatmap
#' components. Use these fields when rebuilding the same screening view in
#' ggplot2, plotly, Quarto, or a dashboard.
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
#'   \item{`"heatmap"`}{Facet A by facet B matrix of signed bias size.
#'     Cells retain reusable matrix and flag tables for dashboards. This is
#'     a Table 13 follow-up display: it supports pattern recognition but does
#'     not turn screening rows into confirmatory tests.}
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
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept bias screening
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
                                  plot = c("scatter", "ranked", "heatmap", "abs_t_hist", "facet_profile"),
                                  diagnostics = NULL,
                                  facet_a = NULL,
                                  facet_b = NULL,
                                  interaction_facets = NULL,
                                  top_n = 40,
                                  abs_t_warn = 2,
                                  abs_bias_warn = 0.5,
                                  p_max = 0.05,
                                  sort_by = c("abs_t", "abs_bias", "prob"),
                                  show_ci = FALSE,
                                  ci_level = 0.95,
                                  main = NULL,
                                  palette = NULL,
                                  label_angle = 45,
                                  preset = c("standard", "publication", "compact", "monochrome"),
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
      show_ci = show_ci,
      ci_level = ci_level,
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
#' For bounded `GPCM`, this helper returns a caveated partial reporting bundle
#' over supported diagnostics, direct tables, and plots. It also includes a
#' `gpcm_boundary` table. Treat the output as slope-aware sensitivity-reporting
#' text, not FACETS score-side equivalence, automatic operational scoring, or
#' design-forecasting evidence.
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
#' # Fast smoke run: a JML fit and a legacy diagnostic let us build the
#' # APA bundle and confirm `report_text` is non-empty in well under
#' # a second.
#' toy <- load_mfrmr_data("example_core")
#' fit_quick <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                       method = "JML", maxit = 30)
#' diag_quick <- diagnose_mfrm(fit_quick, residual_pca = "none",
#'                              diagnostic_mode = "legacy")
#' apa_quick <- build_apa_outputs(fit_quick, diag_quick)
#' nchar(apa_quick$report_text) > 0
#'
#' \donttest{
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 7, maxit = 30)
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
#' # Look for: `SentenceCount` non-zero in every section that the run
#' #   should support (Method / Results / fit / reliability / bias).
#' #   Zero counts mean that section's prose is empty and the
#' #   manuscript will need to fill it manually.
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' head(chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#' # Look for: rows with `DraftReady = "yes"` are ready to paste into
#' #   the manuscript. `"no"` rows tell you which helper / setting
#' #   needs to run before that paragraph can be drafted, via
#' #   `NextAction`. Aim for every Visual Displays / Reliability /
#' #   Diagnostics row to be `"yes"` before submitting.
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "build_apa_outputs()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review",
        "Design planning and forecasting"
      )
    ),
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
        # Match "residual PCA" or "Residual PCA" and also the longer
        # "Exploratory residual PCA" wording that the APA contract uses.
        pat <- "[Rr]esidual PCA"
        grepl(pat, report_text) &&
          grepl(pat, note_text) &&
          grepl(pat, caption_text)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
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
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
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
  if (!is.null(x$gpcm_boundary) && nrow(as.data.frame(x$gpcm_boundary)) > 0) {
    cat("\nGPCM Boundary\n")
    print(as.data.frame(x$gpcm_boundary)[, c("Area", "Status"), drop = FALSE], row.names = FALSE)
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

summary_table_bundle_collapse_value <- function(x) {
  if (is.null(x)) return("")
  if (is.data.frame(x)) return(paste0("<table ", nrow(x), "x", ncol(x), ">"))
  if (is.list(x)) {
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }
  if (length(x) == 0L) return("")
  paste(as.character(x), collapse = "; ")
}

summary_table_bundle_recovery_ademp_df <- function(ademp) {
  if (is.null(ademp) || !is.list(ademp)) return(data.frame())
  dgm <- ademp$data_generating_mechanism %||% list()
  methods <- ademp$methods %||% list()
  out <- data.frame(
    Section = c(
      "Aim",
      "Data-generating mechanism",
      "Data-generating mechanism",
      "Data-generating mechanism",
      "Analysis method",
      "Analysis method",
      "Estimands",
      "Performance measures"
    ),
    Item = c(
      "Aim",
      "Model",
      "Assignment",
      "Step facet",
      "Fitting method",
      "Fitted model",
      "Recovered quantities",
      "Reported metrics"
    ),
    Value = c(
      summary_table_bundle_collapse_value(ademp$aims),
      summary_table_bundle_collapse_value(dgm$model),
      summary_table_bundle_collapse_value(dgm$assignment),
      summary_table_bundle_collapse_value(dgm$step_facet),
      summary_table_bundle_collapse_value(methods$fit_method),
      summary_table_bundle_collapse_value(methods$fitted_model),
      summary_table_bundle_collapse_value(ademp$estimands),
      summary_table_bundle_collapse_value(ademp$performance_measures)
    ),
    stringsAsFactors = FALSE
  )
  out[nzchar(out$Value), , drop = FALSE]
}

summary_table_bundle_sparse_active <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  if (is.numeric(x)) return(is.finite(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "yes", "1")
}

summary_table_bundle_sparse_design_df <- function(x) {
  tbl <- summary_table_bundle_df(x)
  if (nrow(tbl) == 0L || !"SparseDesignActive" %in% names(tbl)) {
    return(data.frame())
  }
  active <- summary_table_bundle_sparse_active(tbl$SparseDesignActive)
  active[is.na(active)] <- FALSE
  if (!any(active)) {
    return(data.frame())
  }
  tbl <- tbl[active, , drop = FALSE]
  review_tbl <- simulation_sparse_design_review_fields(tbl)
  tbl <- cbind(tbl, review_tbl)
  id_cols <- intersect(
    c("Facet", "design_id", "rep", "Seed", "RunOK", "Converged",
      "Observations", "RecoveryRows", "n_person", "n_rater", "n_criterion",
      "raters_per_person"),
    names(tbl)
  )
  alias_cols <- names(tbl)[grepl("^(n_|[A-Za-z0-9_.]+_per_person$)", names(tbl))]
  sparse_cols <- intersect(
    c(
      "SparseDesignActive",
      "DesignDensity",
      "PlannedMissingRate",
      "LinkPersons",
      "LinkFractionActual",
      "LinkRatersPerPerson",
      "MinCommonPersonsPerRaterPair",
      "ZeroCommonRaterPairs",
      "RaterPairsBelowTarget",
      "TargetCommonPersonsPerRaterPair",
      "MeanDesignDensity",
      "MeanPlannedMissingRate",
      "MeanLinkPersons",
      "MeanLinkFractionActual",
      "MeanLinkRatersPerPerson",
      "MeanMinCommonPersonsPerRaterPair",
      "MaxZeroCommonRaterPairs",
      "MaxRaterPairsBelowTarget"
    ),
    names(tbl)
  )
  review_cols <- intersect(c("LinkReviewStatus", "LinkReviewReason", "ReviewUse"), names(tbl))
  keep <- unique(c(id_cols, setdiff(alias_cols, c(sparse_cols, review_cols)), review_cols, sparse_cols))
  tbl[, keep, drop = FALSE]
}

summary_table_bundle_supported_summary_classes <- function() {
  c(
    "summary.mfrm_fit",
    "summary.mfrm_diagnostics",
    "summary.mfrm_precision_review",
    "summary.mfrm_fit_measures",
    "summary.mfrm_facets_fit_review",
    "summary.mfrm_person_fit_indices",
    "summary.mfrm_data_description",
    "summary.mfrm_reporting_checklist",
    "summary.mfrm_apa_outputs",
    "summary.mfrm_design_evaluation",
    "summary.mfrm_signal_detection",
    "summary.mfrm_diagnostic_screening",
    "summary.mfrm_recovery_simulation",
    "summary.mfrm_recovery_assessment",
    "summary.mfrmr_recovery_validation",
    "summary.mfrm_population_prediction",
    "summary.mfrm_future_branch_active_branch",
    "summary.mfrm_facets_run",
    "summary.mfrm_results",
    "summary.mfrm_report",
    "summary.mfrm_bias",
    "summary.mfrm_anchor_review",
    "summary.mfrm_peer_review_design_review",
    "summary.mfrm_network_review",
    "summary.mfrm_linking_review",
    "summary.mfrm_misfit_casebook",
    "summary.mfrm_weighting_review",
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
    cls <- intersect(class(x), summary_classes)[1]
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
  if (inherits(x, "mfrm_precision_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_precision_review",
      summary_class = "summary.mfrm_precision_review"
    ))
  }
  if (inherits(x, "mfrm_fit_measures")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_fit_measures",
      summary_class = "summary.mfrm_fit_measures"
    ))
  }
  if (inherits(x, "mfrm_facets_fit_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_facets_fit_review",
      summary_class = "summary.mfrm_facets_fit_review"
    ))
  }
  if (inherits(x, "mfrm_person_fit_indices")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_person_fit_indices",
      summary_class = "summary.mfrm_person_fit_indices"
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
  if (inherits(x, "mfrm_diagnostic_screening")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_diagnostic_screening",
      summary_class = "summary.mfrm_diagnostic_screening"
    ))
  }
  if (inherits(x, "mfrm_recovery_simulation")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_recovery_simulation",
      summary_class = "summary.mfrm_recovery_simulation"
    ))
  }
  if (inherits(x, "mfrm_recovery_assessment")) {
    return(list(
      summary = summary(x, digits = digits),
      source_class = "mfrm_recovery_assessment",
      summary_class = "summary.mfrm_recovery_assessment"
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
  if (inherits(x, "mfrm_results")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_results",
      summary_class = "summary.mfrm_results"
    ))
  }
  if (inherits(x, "mfrm_report")) {
    return(list(
      summary = summary(x, top_n = top_n),
      source_class = "mfrm_report",
      summary_class = "summary.mfrm_report"
    ))
  }
  if (inherits(x, "mfrm_bias")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_bias",
      summary_class = "summary.mfrm_bias"
    ))
  }
  if (inherits(x, "mfrm_anchor_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_anchor_review",
      summary_class = "summary.mfrm_anchor_review"
    ))
  }
  if (inherits(x, "mfrm_peer_review_design_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_peer_review_design_review",
      summary_class = "summary.mfrm_peer_review_design_review"
    ))
  }
  if (inherits(x, "mfrm_network_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_network_review",
      summary_class = "summary.mfrm_network_review"
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
  if (inherits(x, "mfrm_weighting_review")) {
    return(list(
      summary = summary(x, digits = digits, top_n = top_n),
      source_class = "mfrm_weighting_review",
      summary_class = "summary.mfrm_weighting_review"
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
    "`x` must be an mfrm_fit, mfrm_diagnostics, mfrm_precision_review, ",
    "mfrm_fit_measures, mfrm_facets_fit_review, mfrm_person_fit_indices, ",
    "mfrm_data_description, mfrm_reporting_checklist, mfrm_apa_outputs, ",
    "mfrm_design_evaluation, ",
    "mfrm_signal_detection, mfrm_diagnostic_screening, ",
    "mfrm_recovery_simulation, mfrm_recovery_assessment, ",
    "mfrm_population_prediction, mfrm_future_branch_active_branch, ",
    "mfrm_facets_run, mfrm_results, mfrm_report, mfrm_bias, mfrm_anchor_review, ",
    "mfrm_peer_review_design_review, mfrm_network_review, ",
    "mfrm_linking_review, mfrm_misfit_casebook, ",
    "mfrm_weighting_review, mfrm_unit_prediction, or ",
    "mfrm_plausible_values object, or one of their summary() outputs.",
    call. = FALSE
  )
}

summary_table_bundle_required_components <- function(summary_class) {
  switch(
    as.character(summary_class %||% NA_character_),
    "summary.mfrm_fit" = c("overview", "reporting_map"),
    "summary.mfrm_diagnostics" = c("overview", "reporting_map", "flags"),
    "summary.mfrm_precision_review" = c("overview", "summary", "profile", "checks", "fit_separation_basis"),
    "summary.mfrm_fit_measures" = c("overview", "summary", "status_summary", "table"),
    "summary.mfrm_facets_fit_review" = c("overview", "summary", "df_sensitivity", "guidance"),
    "summary.mfrm_person_fit_indices" = c("overview", "status_summary", "top_review"),
    "summary.mfrm_data_description" = c("overview", "score_distribution"),
    "summary.mfrm_reporting_checklist" = c("overview", "action_items"),
    "summary.mfrm_apa_outputs" = c("overview", "components", "preview"),
    "summary.mfrm_design_evaluation" = c("overview", "design_summary"),
    "summary.mfrm_signal_detection" = c("overview", "detection_summary"),
    "summary.mfrm_diagnostic_screening" = c("overview", "reading_order", "next_actions", "reporting_notes", "figure_recipes", "scenario_summary", "performance_summary", "plot_overview_rate"),
    "summary.mfrm_recovery_simulation" = c("overview", "recovery_summary", "rep_overview"),
    "summary.mfrm_recovery_assessment" = c("overview", "reading_order", "checklist", "condition_reporting_notes", "condition_review", "diagnostic_reporting_notes", "diagnostic_review", "metric_review", "uncertainty_review"),
    "summary.mfrmr_recovery_validation" = c("topline_release_decision", "reading_order", "release_decision_table", "case_summary", "condition_reporting_notes", "condition_summary", "diagnostic_reporting_notes", "diagnostic_oc_summary", "domain_decision_table"),
    "summary.mfrm_population_prediction" = c("overview", "design", "forecast"),
    "summary.mfrm_future_branch_active_branch" = c("overview", "profile_summary", "recommendation_table"),
    "summary.mfrm_facets_run" = c("overview", "mapping", "run_info", "fit", "diagnostics"),
    "summary.mfrm_results" = c("overview", "triage", "status", "component_index", "table_index", "plot_map", "next_actions"),
    "summary.mfrm_report" = c("overview", "first_screen", "status_counts", "immediate_actions", "optional_sections", "claim_readiness", "report_gaps", "boundary_index", "routes"),
    "summary.mfrm_bias" = c("overview", "top_rows"),
    "summary.mfrm_anchor_review" = c("facet_summary", "recommendations"),
    "summary.mfrm_peer_review_design_review" = c("overview", "load_summary", "low_common_pairs", "reporting_map"),
    "summary.mfrm_network_review" = c("overview", "network_summary", "reporting_map"),
    "summary.mfrm_linking_review" = c("overview", "top_linking_risks", "group_view_index", "reporting_map"),
    "summary.mfrm_misfit_casebook" = c("overview", "top_cases", "case_rollup", "group_view_index", "reporting_map"),
    "summary.mfrm_weighting_review" = c("overview", "top_reweighted_levels", "reporting_map"),
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
        population_overview = "Population-model basis, posterior basis, and omission review.",
        population_design = "Population-model design-matrix columns and numeric check statistics.",
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
        precision_review = summary_table_bundle_df(summary_obj$precision_review),
        reliability = summary_table_bundle_df(summary_obj$reliability),
        top_fit = summary_table_bundle_df(summary_obj$top_fit),
        reporting_map = summary_table_bundle_df(summary_obj$reporting_map),
        flags = summary_table_bundle_df(summary_obj$flags)
      ),
      roles = c(
        overview = "run_overview",
        overall_fit = "overall_fit",
        precision_profile = "precision_basis",
        precision_review = "precision_review",
        reliability = "facet_precision",
        top_fit = "extreme_fit_rows",
        reporting_map = "reporting_map",
        flags = "flag_counts"
      ),
      descriptions = c(
        overview = "Run-level diagnostic coverage and precision tier.",
        overall_fit = "Global fit statistics from the current diagnostic run.",
        precision_profile = "Precision basis and recommended interpretation tier.",
        precision_review = "Precision checks marked review/warn for manuscript caution.",
        reliability = "Facet-level separation, strata, and reliability summary.",
        top_fit = "Rows with the largest absolute fit Z statistics.",
        reporting_map = "Companion outputs for manuscript reporting beyond summary(diag).",
        flags = "Counts of unexpected responses, displacement, interactions, and inter-rater pairs."
      )
    ),
    "summary.mfrm_precision_review" = list(
      title = "Precision Review Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        summary = summary_table_bundle_df(summary_obj$summary),
        profile = summary_table_bundle_df(summary_obj$profile),
        checks = summary_table_bundle_df(summary_obj$checks),
        fit_separation_basis = summary_table_bundle_df(summary_obj$fit_separation_basis),
        approximation_notes = summary_table_bundle_df(summary_obj$approximation_notes),
        settings = summary_table_bundle_df(summary_obj$settings),
        caveats = summary_table_bundle_df(summary_obj$caveats),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "run_overview",
        summary = "precision_basis",
        profile = "precision_basis",
        checks = "precision_review",
        fit_separation_basis = "precision_review",
        approximation_notes = "interpretation_notes",
        settings = "review_settings",
        caveats = "analysis_caveats",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level precision-review bundle metadata.",
        summary = "One-row precision tier, formal-inference support, and fit/separation row counts.",
        profile = "Precision basis, method tier, and formal-inference support for the current run.",
        checks = "Precision checks that should be reviewed before reporting SE, CI, or reliability claims.",
        fit_separation_basis = "Source-grounded boundary table for fit MnSq, ZSTD, separation/reliability/strata, and QC-threshold interpretation.",
        approximation_notes = "Approximation notes carried from diagnostics for uncertainty-reporting caveats.",
        settings = "Precision-review settings used for the current run.",
        caveats = "Structured fit-level caveats carried with the precision review.",
        notes = "Compact interpretation notes for precision, fit, and separation reporting."
      )
    ),
    "summary.mfrm_fit_measures" = list(
      title = "Fit-Measure Review Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        summary = summary_table_bundle_df(summary_obj$summary),
        table = summary_table_bundle_df(summary_obj$table),
        facets_table = summary_table_bundle_df(summary_obj$facets_table),
        status_summary = summary_table_bundle_df(summary_obj$status_summary),
        threshold_profiles = summary_table_bundle_df(summary_obj$threshold_profiles),
        profile_summary = summary_table_bundle_df(summary_obj$profile_summary),
        profile_summary_by_facet = summary_table_bundle_df(summary_obj$profile_summary_by_facet),
        profile_summary_overall = summary_table_bundle_df(summary_obj$profile_summary_overall),
        df_sensitivity = summary_table_bundle_df(summary_obj$df_sensitivity),
        df_sensitive = summary_table_bundle_df(summary_obj$df_sensitive),
        df_sensitivity_summary = summary_table_bundle_df(summary_obj$df_sensitivity_summary),
        underfit = summary_table_bundle_df(summary_obj$underfit),
        overfit = summary_table_bundle_df(summary_obj$overfit),
        mixed = summary_table_bundle_df(summary_obj$mixed),
        settings = summary_table_bundle_df(summary_obj$settings),
        caveats = summary_table_bundle_df(summary_obj$caveats),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "run_overview",
        summary = "overall_fit",
        table = "extreme_fit_rows",
        facets_table = "extreme_fit_rows",
        status_summary = "review_status",
        threshold_profiles = "review_settings",
        profile_summary = "overall_fit",
        profile_summary_by_facet = "facet_distribution",
        profile_summary_overall = "overall_fit",
        df_sensitivity = "precision_review",
        df_sensitive = "precision_review",
        df_sensitivity_summary = "review_status",
        underfit = "extreme_fit_rows",
        overfit = "extreme_fit_rows",
        mixed = "extreme_fit_rows",
        settings = "review_settings",
        caveats = "analysis_caveats",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level fit-measure bundle metadata.",
        summary = "Overall fit-status counts, df-sensitivity counts, and displayed-row coverage.",
        table = "Selected fit-measure rows with MnSq, ZSTD, confidence intervals, status labels, and df-sensitivity columns.",
        facets_table = "FACETS-style fit table columns for manuscript or external-output review.",
        status_summary = "Counts by facet and fit-review status.",
        threshold_profiles = "Fit-threshold profiles used for status and sensitivity summaries.",
        profile_summary = "Fit-flag rates under the requested threshold profiles.",
        profile_summary_by_facet = "Facet-level fit-flag rates under the requested threshold profiles.",
        profile_summary_overall = "Overall fit-flag rates under the requested threshold profiles.",
        df_sensitivity = "Engine-vs-FACETS-style df/ZSTD comparison rows.",
        df_sensitive = "Subset of df-sensitivity rows where df convention changes flag status or materially shifts ZSTD.",
        df_sensitivity_summary = "Counts by df-sensitivity status for the selected fit-measure surface.",
        underfit = "Rows labelled underfit under the selected fit bands.",
        overfit = "Rows labelled overfit under the selected fit bands.",
        mixed = "Rows with mixed underfit and overfit evidence across fit columns.",
        settings = "Fit-measure review settings and thresholds.",
        caveats = "Structured fit-level caveats carried with the fit-measure review.",
        notes = "Compact interpretation notes for fit-measure reporting."
      )
    ),
    "summary.mfrm_facets_fit_review" = list(
      title = "FACETS Fit-Review Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        summary = summary_table_bundle_df(summary_obj$summary),
        standardization = summary_table_bundle_df(summary_obj$standardization),
        df_sensitivity = summary_table_bundle_df(summary_obj$df_sensitivity),
        df_sensitive = summary_table_bundle_df(summary_obj$df_sensitive),
        df_sensitivity_summary = summary_table_bundle_df(summary_obj$df_sensitivity_summary),
        external_table_quality = summary_table_bundle_df(summary_obj$external_table_quality),
        external_comparison = summary_table_bundle_df(summary_obj$external_comparison),
        guidance = summary_table_bundle_df(summary_obj$guidance),
        settings = summary_table_bundle_df(summary_obj$settings),
        caveats = summary_table_bundle_df(summary_obj$caveats),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "run_overview",
        summary = "overall_fit",
        standardization = "review_settings",
        df_sensitivity = "precision_review",
        df_sensitive = "precision_review",
        df_sensitivity_summary = "review_status",
        external_table_quality = "review_status",
        external_comparison = "precision_review",
        guidance = "interpretation_notes",
        settings = "review_settings",
        caveats = "analysis_caveats",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level FACETS fit-review bundle metadata.",
        summary = "Overview of internal df/ZSTD sensitivity and optional external FACETS comparison coverage.",
        standardization = "Primary fit-df method, companion columns, and ZSTD transform metadata from diagnostics.",
        df_sensitivity = "Engine-vs-FACETS-style df/ZSTD comparison rows.",
        df_sensitive = "Subset of df-sensitivity rows where df convention changes flag status or materially shifts ZSTD.",
        df_sensitivity_summary = "Counts by df-sensitivity status.",
        external_table_quality = "Completeness and duplicate-key review for supplied external FACETS fit rows.",
        external_comparison = "Matched external FACETS-vs-mfrmr fit comparison rows when supplied.",
        guidance = "Interpretation notes separating MnSq, ZSTD, df convention, external matching, and GPCM scope.",
        settings = "FACETS fit-review tolerances and metadata.",
        caveats = "Structured fit-level caveats carried with the FACETS fit review.",
        notes = "Compact interpretation notes for FACETS fit-review reporting."
      )
    ),
    "summary.mfrm_person_fit_indices" = list(
      title = "Person-Fit Summary Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        status_summary = summary_table_bundle_df(summary_obj$status_summary),
        report_index_summary = summary_table_bundle_df(summary_obj$report_index_summary),
        lz_star_status_summary = summary_table_bundle_df(summary_obj$lz_star_status_summary),
        top_review = summary_table_bundle_df(summary_obj$top_review),
        caveats = summary_table_bundle_df(summary_obj$caveats),
        thresholds = summary_table_bundle_df(summary_obj$thresholds),
        reporting_map = summary_table_bundle_df(summary_obj$reporting_map),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "overall_fit",
        status_summary = "review_status",
        report_index_summary = "review_status",
        lz_star_status_summary = "review_status",
        top_review = "extreme_fit_rows",
        caveats = "analysis_caveats",
        thresholds = "review_settings",
        reporting_map = "reporting_map",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level person-fit counts, reportable rows, Snijders-corrected rows, and flag rates.",
        status_summary = "Counts by report-level person-fit review status.",
        report_index_summary = "Counts showing whether the report index used lz_star, lz, or no finite statistic.",
        lz_star_status_summary = "Counts by Snijders-correction availability/status.",
        top_review = "Highest-priority person rows for response-level follow-up.",
        caveats = "Visible caveats explaining whether lz_star or uncorrected lz underlies each reporting route.",
        thresholds = "Practical two-sided z thresholds used for person-fit flags.",
        reporting_map = "Companion outputs for response-level follow-up and draw-free plot-data handoff.",
        notes = "Compact interpretation notes for person-fit reporting."
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
        facets_positioning = summary_table_bundle_df(summary_obj$facets_positioning),
        priority_summary = summary_table_bundle_df(summary_obj$priority_summary),
        action_items = summary_table_bundle_df(summary_obj$action_items),
        settings = summary_table_bundle_df(summary_obj$settings)
      ),
      roles = c(
        overview = "checklist_overview",
        section_summary = "section_coverage",
        facets_positioning = "facets_relationship_wording",
        priority_summary = "priority_distribution",
        action_items = "draft_actions",
        settings = "checklist_settings"
      ),
      descriptions = c(
        overview = "Overall checklist coverage across sections and draft-readiness flags.",
        section_summary = "Coverage summary by reporting section.",
        facets_positioning = "Report-ready wording that separates mfrmr estimation from FACETS-style handoff or external-table review.",
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
            design_summary = summary_table_bundle_df(summary_obj$design_summary),
            sparse_review = summary_table_bundle_df(summary_obj$sparse_review),
            sparse_design = summary_table_bundle_sparse_design_df(summary_obj$design_summary)
          ),
          future_spec$tables
        ),
        roles = c(
          overview = "run_overview",
          design_summary = "design_performance",
          sparse_review = "sparse_design_diagnostics",
          sparse_design = "sparse_design_diagnostics",
          future_spec$roles
        ),
        descriptions = c(
          overview = "Run-level overview for the current design-evaluation study.",
          design_summary = "Aggregated Monte Carlo design summaries for the active two-role planner.",
          sparse_review = "Compact sparse linked design-review counts for planned missingness and rater-pair linkage.",
          sparse_design = "Sparse linked planned-missingness and rater-link diagnostics for design-evaluation rows.",
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
    "summary.mfrm_diagnostic_screening" = list(
      title = "Diagnostic Screening Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        reading_order = summary_table_bundle_df(summary_obj$reading_order),
        next_actions = summary_table_bundle_df(summary_obj$next_actions),
        reporting_notes = summary_table_bundle_df(summary_obj$reporting_notes),
        figure_recipes = summary_table_bundle_df(summary_obj$figure_recipes),
        scenario_summary = summary_table_bundle_df(summary_obj$scenario_summary),
        performance_summary = summary_table_bundle_df(summary_obj$performance_summary),
        report_signal_summary = summary_table_bundle_df(summary_obj$report_signal_summary),
        scenario_contrast = summary_table_bundle_df(summary_obj$scenario_contrast),
        plot_overview_rate = summary_table_bundle_df(summary_obj$plot_overview_rate),
        plot_overview_count = summary_table_bundle_df(summary_obj$plot_overview_count),
        plot_report_rate = summary_table_bundle_df(summary_obj$plot_report_rate),
        plot_contrast_count = summary_table_bundle_df(summary_obj$plot_contrast_count),
        plot_runtime = summary_table_bundle_df(summary_obj$plot_runtime),
        ademp = summary_table_bundle_recovery_ademp_df(summary_obj$ademp),
        settings = summary_table_bundle_settings_df(summary_obj$settings),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "diagnostic_screening_overview",
        reading_order = "diagnostic_screening_reading_order",
        next_actions = "diagnostic_screening_next_actions",
        reporting_notes = "diagnostic_screening_reporting_notes",
        figure_recipes = "diagnostic_screening_figure_recipes",
        scenario_summary = "diagnostic_screening_scenario_summary",
        performance_summary = "diagnostic_screening_performance",
        report_signal_summary = "diagnostic_screening_report_signals",
        scenario_contrast = "diagnostic_screening_contrast",
        plot_overview_rate = "diagnostic_screening_plot_data",
        plot_overview_count = "diagnostic_screening_plot_data",
        plot_report_rate = "diagnostic_screening_plot_data",
        plot_contrast_count = "diagnostic_screening_plot_data",
        plot_runtime = "diagnostic_screening_runtime",
        ademp = "diagnostic_screening_design_basis",
        settings = "diagnostic_screening_settings",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level overview for the current diagnostic-screening simulation study.",
        reading_order = "Recommended reading order for diagnostic-screening summaries and appendix tables.",
        next_actions = "Action-oriented triage for replication count, run completion, contrasts, report signals, and export.",
        reporting_notes = "Report-facing boundaries and recommended wording safeguards for diagnostic-screening output.",
        figure_recipes = "Figure/display recipes linking plot calls, plot_data extraction, caption focus, and interpretation boundaries.",
        scenario_summary = "Aggregated scenario-by-design legacy and strict marginal/pairwise screening summaries.",
        performance_summary = "Scenario-by-design any-flag rates, agreement rates, elapsed time, and Type I/sensitivity proxies.",
        report_signal_summary = "Optional report-index availability/readiness and review-signal summaries when include_report = TRUE.",
        scenario_contrast = "Misspecification-minus-well-specified contrasts for legacy and strict screening signals.",
        plot_overview_rate = "Long-form overview-rate plot data for legacy, strict, and optional report-review signals.",
        plot_overview_count = "Long-form overview-count plot data for flagged counts and report-signal counts.",
        plot_report_rate = "Long-form report-rate plot data from report_signal_summary.",
        plot_contrast_count = "Long-form scenario-contrast count plot data.",
        plot_runtime = "Long-form runtime plot data.",
        ademp = "ADEMP-style description of the diagnostic-screening simulation aim, data-generating mechanism, estimands, methods, and performance measures.",
        settings = "Diagnostic-screening simulation settings used to generate, fit, diagnose, and optionally report repeated datasets.",
        notes = "Compact interpretation notes for diagnostic-screening reporting."
      )
    ),
    "summary.mfrm_recovery_simulation" = list(
      title = "Recovery Simulation Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        recovery_summary = summary_table_bundle_df(summary_obj$recovery_summary),
        rep_overview = summary_table_bundle_df(summary_obj$rep_overview),
        sparse_review = summary_table_bundle_df(summary_obj$sparse_review),
        sparse_design = summary_table_bundle_sparse_design_df(summary_obj$rep_overview),
        diagnostic_oc = summary_table_bundle_df(summary_obj$diagnostic_oc),
        diagnostic_oc_summary = summary_table_bundle_df(summary_obj$diagnostic_oc_summary),
        ademp = summary_table_bundle_recovery_ademp_df(summary_obj$ademp),
        settings = summary_table_bundle_settings_df(summary_obj$settings),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "recovery_overview",
        recovery_summary = "recovery_performance",
        rep_overview = "recovery_replications",
        sparse_review = "recovery_sparse_design_diagnostics",
        sparse_design = "recovery_sparse_design_diagnostics",
        diagnostic_oc = "recovery_diagnostic_operating_characteristics",
        diagnostic_oc_summary = "recovery_diagnostic_oc_summary",
        ademp = "recovery_design_basis",
        settings = "recovery_settings",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level overview for the current parameter-recovery simulation.",
        recovery_summary = "Parameter-group recovery metrics, including bias, RMSE, coverage, and Monte Carlo SE.",
        rep_overview = "Replication-level fit status, convergence status, recovery-row counts, and elapsed time.",
        sparse_review = "Compact replication-level sparse linked design-review counts kept separate from recovery metrics.",
        sparse_design = "Replication-level sparse linked planned-missingness and rater-link diagnostics retained separately from recovery metrics.",
        diagnostic_oc = "Optional replication-by-facet fit/separation operating characteristics retained when include_diagnostics = TRUE.",
        diagnostic_oc_summary = "Optional facet-level fit/separation operating-characteristic summary for diagnostic context.",
        ademp = "ADEMP-style description of the recovery simulation aim, data-generating mechanism, estimands, methods, and performance measures.",
        settings = "Recovery simulation settings used to generate and refit repeated datasets.",
        notes = "Compact interpretation notes for recovery simulation reporting."
      )
    ),
    "summary.mfrm_recovery_assessment" = list(
      title = "Recovery Assessment Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        reading_order = summary_table_bundle_df(summary_obj$reading_order),
        checklist = summary_table_bundle_df(summary_obj$checklist),
        condition_reporting_notes = summary_table_bundle_df(summary_obj$condition_reporting_notes),
        condition_review = summary_table_bundle_df(summary_obj$condition_review),
        diagnostic_reporting_notes = summary_table_bundle_df(summary_obj$diagnostic_reporting_notes),
        diagnostic_review = summary_table_bundle_df(summary_obj$diagnostic_review),
        metric_review = summary_table_bundle_df(summary_obj$metric_review),
        uncertainty_review = summary_table_bundle_df(summary_obj$uncertainty_review),
        next_actions = summary_table_bundle_text_df(summary_obj$next_actions, column = "Action"),
        thresholds = summary_table_bundle_settings_df(summary_obj$thresholds),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "recovery_assessment_overview",
        reading_order = "recovery_assessment_reading_order",
        checklist = "recovery_assessment_checklist",
        condition_reporting_notes = "recovery_condition_reporting_notes",
        condition_review = "recovery_condition_review",
        diagnostic_reporting_notes = "recovery_diagnostic_reporting_notes",
        diagnostic_review = "recovery_diagnostic_review",
        metric_review = "recovery_metric_review",
        uncertainty_review = "recovery_uncertainty_review",
        next_actions = "repair_recommendations",
        thresholds = "review_settings",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "Run-level recovery adequacy status for the current assessment.",
        reading_order = "Recommended first-read order for recovery assessment summary, condition, plot, and row-level outputs.",
        checklist = "Reviewer-facing adequacy checklist for replication count, convergence, uncertainty, Monte Carlo precision, and practical thresholds.",
        condition_reporting_notes = "Reporter-facing generator-condition notes for bounded-GPCM slope stress and generated score-category support.",
        condition_review = "Generator-condition metadata for interpreting recovery evidence, including bounded-GPCM slope-regime labels and generated score-category support when available.",
        diagnostic_reporting_notes = "Reporter-facing fit/separation diagnostic notes that flag caveats without treating them as recovery gates.",
        diagnostic_review = "Fit/separation operating-characteristic review retained as diagnostic context rather than a release-recovery gate.",
        metric_review = "Parameter-group recovery review with threshold status and next-action guidance.",
        uncertainty_review = "Parameter-group coverage and standard-error availability interpretation.",
        next_actions = "Prioritized follow-up actions for strengthening or documenting the recovery evidence.",
        thresholds = "Assessment thresholds used to classify recovery adequacy.",
        notes = "Compact interpretation notes for recovery assessment reporting."
      )
    ),
    "summary.mfrmr_recovery_validation" = list(
      title = "Recovery Validation Tables",
      tables = list(
        topline_release_decision = summary_table_bundle_df(summary_obj$topline_release_decision),
        reading_order = summary_table_bundle_df(summary_obj$reading_order),
        release_decision_table = summary_table_bundle_df(summary_obj$release_decision_table),
        case_summary = summary_table_bundle_df(summary_obj$case_summary),
        condition_reporting_notes = summary_table_bundle_df(summary_obj$condition_reporting_notes),
        condition_summary = summary_table_bundle_df(summary_obj$condition_summary),
        diagnostic_reporting_notes = summary_table_bundle_df(summary_obj$diagnostic_reporting_notes),
        diagnostic_oc_summary = summary_table_bundle_df(summary_obj$diagnostic_oc_summary),
        domain_decision_table = summary_table_bundle_df(summary_obj$domain_decision_table)
      ),
      roles = c(
        topline_release_decision = "recovery_validation_topline",
        reading_order = "recovery_validation_reading_order",
        release_decision_table = "recovery_validation_release_decisions",
        case_summary = "recovery_validation_case_summary",
        condition_reporting_notes = "recovery_validation_condition_reporting_notes",
        condition_summary = "recovery_validation_condition_summary",
        diagnostic_reporting_notes = "recovery_validation_diagnostic_reporting_notes",
        diagnostic_oc_summary = "recovery_validation_diagnostic_oc_summary",
        domain_decision_table = "recovery_validation_domain_decisions"
      ),
      descriptions = c(
        topline_release_decision = "Top-line release-recovery decision across the summarized validation cases.",
        reading_order = "Recommended first-read order for recovery-validation summary outputs.",
        release_decision_table = "Case-level release-recovery decision table with recovery, uncertainty, and Monte Carlo status.",
        case_summary = "Case-level validation summary including recovery status and generator-condition fields.",
        condition_reporting_notes = "Reporter-facing generator-condition notes for slope-regime and score-support caveats.",
        condition_summary = "Generator-condition summary for slope-regime and score-support evidence across validation cases.",
        diagnostic_reporting_notes = "Reporter-facing fit/separation diagnostic notes that flag caveats without treating them as release gates.",
        diagnostic_oc_summary = "Fit/separation operating-characteristic summary across validation cases, retained as diagnostic-only context.",
        domain_decision_table = "Long-form validation-domain status table for recovery metrics, uncertainty, Monte Carlo precision, score support, and overall status."
      )
    ),
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
    "summary.mfrm_results" = list(
      title = "Comprehensive Results Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        triage = summary_table_bundle_df(summary_obj$triage),
        status = summary_table_bundle_df(summary_obj$status),
        component_index = summary_table_bundle_df(summary_obj$component_index),
        table_index = summary_table_bundle_df(summary_obj$table_index),
        plot_map = summary_table_bundle_df(summary_obj$plot_map),
        next_actions = summary_table_bundle_df(summary_obj$next_actions),
        mapping = summary_table_bundle_df(summary_obj$mapping),
        reproducible_code = summary_table_bundle_df(summary_obj$reproducible_code),
        notes = summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      ),
      roles = c(
        overview = "workflow_overview",
        triage = "review_status",
        status = "review_status",
        component_index = "reporting_map",
        table_index = "reporting_map",
        plot_map = "plot_routing",
        next_actions = "draft_actions",
        mapping = "column_mapping",
        reproducible_code = "workflow_settings",
        notes = "interpretation_notes"
      ),
      descriptions = c(
        overview = "High-level mfrm_results run overview with component, table, and plot-route counts.",
        triage = "First-screen reading order for unavailable, review, information, and OK signals.",
        status = "Section-level availability table for automatically assembled result components.",
        component_index = "Classes of the fitted, diagnostic, report, review, and table components retained in the object.",
        table_index = "Available table registry with size and numeric plotting readiness.",
        plot_map = "User-facing plot routes exposed by plot.mfrm_results().",
        next_actions = "Prioritized next-action routes after the comprehensive first screen.",
        mapping = "Column mapping used when mfrm_results() started from a data.frame or run_mfrm_facets() object.",
        reproducible_code = "Line-by-line replay scaffold for the mfrm_results() route.",
        notes = "Compact interpretation notes carried by mfrm_results()."
      )
    ),
    "summary.mfrm_report" = list(
      title = "Report Summary Tables",
      tables = list(
        overview = summary_table_bundle_df(summary_obj$overview),
        first_screen = summary_table_bundle_df(summary_obj$first_screen),
        status_counts = summary_table_bundle_df(summary_obj$status_counts),
        immediate_actions = summary_table_bundle_df(summary_obj$immediate_actions),
        optional_sections = summary_table_bundle_df(summary_obj$optional_sections),
        claim_readiness = summary_table_bundle_df(summary_obj$claim_readiness),
        report_gaps = summary_table_bundle_df(summary_obj$report_gaps),
        boundary_index = summary_table_bundle_df(summary_obj$boundary_index),
        routes = summary_table_bundle_df(summary_obj$routes)
      ),
      roles = c(
        overview = "report_overview",
        first_screen = "first_screen",
        status_counts = "review_status",
        immediate_actions = "draft_actions",
        optional_sections = "requested_followup",
        claim_readiness = "claim_readiness",
        report_gaps = "draft_actions",
        boundary_index = "claim_boundary",
        routes = "reporting_map"
      ),
      descriptions = c(
        overview = "Reader-facing report status row derived from mfrm_report() first_screen.",
        first_screen = "Compact FACETS-like first screen with status, main issue, next action, and route.",
        status_counts = "Counts of evidence areas by first-screen status.",
        immediate_actions = "Unavailable, review, and caveated areas to inspect before drafting.",
        optional_sections = "Not-requested evidence areas to request only when the claim is needed.",
        claim_readiness = "Claim-readiness counts summarized from the report object.",
        report_gaps = "Highest-priority report gaps and recommended actions.",
        boundary_index = "Template-boundary rows to review before using report wording.",
        routes = "Standard report-reading and export routes."
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
    "summary.mfrm_anchor_review" = {
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
        title = "Anchor Review Tables",
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
          overview = "anchor_review_overview",
          issue_counts = "anchor_issue_counts",
          facet_summary = "facet_coverage",
          level_observation_summary = "level_observation_review",
          category_counts = "category_usage",
          recommendations = "repair_recommendations",
          notes = "interpretation_notes"
        ),
        descriptions = c(
          overview = "Anchor-review overview with issue, facet, and recommendation counts.",
          issue_counts = "Observed anchor-review issue counts ranked by frequency.",
          facet_summary = "Facet-level counts and anchor-table coverage summary.",
          level_observation_summary = "Observation counts by facet level for anchor viability checks.",
          category_counts = "Observed score-category usage for anchor-review screening.",
          recommendations = "Compact action list for anchor repair or review.",
          notes = "One-line interpretation note from the anchor review."
        )
      )
    },
    "summary.mfrm_peer_review_design_review" = {
      overview_tbl <- summary_table_bundle_df(summary_obj$overview)
      load_summary_tbl <- summary_table_bundle_df(summary_obj$load_summary)
      submission_load_tbl <- summary_table_bundle_df(summary_obj$submission_load)
      reviewer_load_tbl <- summary_table_bundle_df(summary_obj$reviewer_load)
      common_tbl <- summary_table_bundle_df(summary_obj$reviewer_pair_common_submissions)
      low_common_tbl <- summary_table_bundle_df(summary_obj$low_common_pairs)
      reciprocal_tbl <- summary_table_bundle_df(summary_obj$reciprocal_pairs)
      reporting_tbl <- summary_table_bundle_df(summary_obj$reporting_map)
      caveat_tbl <- summary_table_bundle_df(summary_obj$caveats)
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      list(
        title = "Peer-Review Design Tables",
        tables = list(
          overview = overview_tbl,
          load_summary = load_summary_tbl,
          submission_load = submission_load_tbl,
          reviewer_load = reviewer_load_tbl,
          reviewer_pair_common_submissions = common_tbl,
          low_common_pairs = low_common_tbl,
          reciprocal_pairs = reciprocal_tbl,
          reporting_map = reporting_tbl,
          caveats = caveat_tbl,
          notes = notes_tbl,
          settings = settings_tbl
        ),
        roles = c(
          overview = "peer_review_design_diagnostics",
          load_summary = "peer_review_load_summary",
          submission_load = "peer_review_submission_load",
          reviewer_load = "peer_review_reviewer_load",
          reviewer_pair_common_submissions = "peer_review_common_links",
          low_common_pairs = "peer_review_low_common_links",
          reciprocal_pairs = "peer_review_reciprocal_pairs",
          reporting_map = "reporting_map",
          caveats = "analysis_caveats",
          notes = "interpretation_notes",
          settings = "review_settings"
        ),
        descriptions = c(
          overview = "Front-door peer-review assignment status and design-density diagnostics.",
          load_summary = "Submission and reviewer assignment-load summary.",
          submission_load = "Highest-load submission rows for assignment follow-up.",
          reviewer_load = "Highest-load reviewer rows for assignment follow-up.",
          reviewer_pair_common_submissions = "Reviewer-pair common-submission links retained for design review.",
          low_common_pairs = "Reviewer pairs with the fewest common submissions.",
          reciprocal_pairs = "Participant pairs that reviewed each other's submissions.",
          reporting_map = "Map separating peer-review assignment diagnostics from MFRM measurement evidence.",
          caveats = "Interpretation caveats for peer-review design diagnostics.",
          notes = "Compact interpretation notes for peer-review design review.",
          settings = "Settings recorded by build_peer_review_design_review()."
        )
      )
    },
    "summary.mfrm_network_review" = {
      overview_tbl <- summary_table_bundle_df(summary_obj$overview)
      network_tbl <- summary_table_bundle_df(summary_obj$network_summary)
      facet_tbl <- summary_table_bundle_df(summary_obj$facet_summary)
      central_tbl <- summary_table_bundle_df(summary_obj$top_central_nodes)
      cut_tbl <- summary_table_bundle_df(summary_obj$top_cut_nodes)
      bridge_tbl <- summary_table_bundle_df(summary_obj$top_bridge_edges)
      sparse_tbl <- summary_table_bundle_df(summary_obj$sparse_review)
      peer_tbl <- summary_table_bundle_df(summary_obj$peer_review)
      reporting_tbl <- summary_table_bundle_df(summary_obj$reporting_map)
      caveat_tbl <- summary_table_bundle_df(summary_obj$caveats)
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      list(
        title = "Network Review Tables",
        tables = list(
          overview = overview_tbl,
          network_summary = network_tbl,
          facet_summary = facet_tbl,
          top_central_nodes = central_tbl,
          top_cut_nodes = cut_tbl,
          top_bridge_edges = bridge_tbl,
          sparse_review = sparse_tbl,
          peer_review = peer_tbl,
          reporting_map = reporting_tbl,
          caveats = caveat_tbl,
          notes = notes_tbl,
          settings = settings_tbl
        ),
        roles = c(
          overview = "network_review_overview",
          network_summary = "network_design_summary",
          facet_summary = "network_facet_vulnerability",
          top_central_nodes = "network_central_nodes",
          top_cut_nodes = "network_articulation_nodes",
          top_bridge_edges = "network_bridge_edges",
          sparse_review = "sparse_design_diagnostics",
          peer_review = "peer_review_design_diagnostics",
          reporting_map = "reporting_map",
          caveats = "analysis_caveats",
          notes = "interpretation_notes",
          settings = "review_settings"
        ),
        descriptions = c(
          overview = "Front-door connectedness and design-network review status.",
          network_summary = "Graph-level connectedness, density, articulation-point, and bridge-edge metrics.",
          facet_summary = "Facet-level aggregation of node centrality and design-link vulnerability.",
          top_central_nodes = "Highest-betweenness design-network nodes for follow-up inspection.",
          top_cut_nodes = "Articulation-point rows whose removal would fragment the design graph.",
          top_bridge_edges = "Bridge-edge rows indicating one-link dependencies between graph regions.",
          sparse_review = "Optional sparse linked design-review counts for planned missingness and rater-pair linkage.",
          peer_review = "Optional peer-review assignment, load, and common-submission linkage diagnostics.",
          reporting_map = "Map separating MFRM measurement, design-network review, sparse-design review, peer-review design, and rater-effect network routes.",
          caveats = "Interpretation caveats for network-design diagnostics.",
          notes = "Compact interpretation notes for network-design review.",
          settings = "Settings and provenance recorded by build_mfrm_network_review()."
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
    "summary.mfrm_weighting_review" = {
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
        title = "Weighting Review Tables",
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
          overview = "Overview of the equal-weighting versus bounded GPCM weighting review.",
          status = "Compact status block for the weighting-policy review.",
          top_measure_shifts = "Largest non-person facet-measure shifts between the Rasch-family reference and bounded GPCM.",
          top_reweighted_levels = "Largest slope-facet reweighting signals under bounded GPCM.",
          plot_map = "Public plot routes for precision redistribution and comparison follow-up.",
          reporting_map = "Bundle/report handoff map for weighting-policy review outputs.",
          support_status = "Capability-boundary statement for the bounded GPCM weighting review.",
          key_warnings = "Top warning lines for weighting-policy review.",
          next_actions = "Recommended next-step actions after weighting-policy review.",
          notes = "Interpretation notes for the weighting review.",
          settings = "Weighting-review settings and theta-grid parameters."
        )
      )
    },
    "summary.mfrm_unit_prediction" = {
      estimate_tbl <- summary_table_bundle_df(summary_obj$estimates)
      row_review_tbl <- summary_table_bundle_df(summary_obj$row_review)
      population_review_tbl <- summary_table_bundle_df(summary_obj$population_review)
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      overview_tbl <- data.frame(
        Units = nrow(estimate_tbl),
        RowReviewRows = nrow(row_review_tbl),
        PopulationReviewRows = nrow(population_review_tbl),
        Settings = nrow(settings_tbl),
        Notes = nrow(notes_tbl),
        stringsAsFactors = FALSE
      )
      list(
        title = "Unit Prediction Tables",
        tables = list(
          overview = overview_tbl,
          estimates = estimate_tbl,
          row_review = row_review_tbl,
          population_review = population_review_tbl,
          settings = settings_tbl,
          notes = notes_tbl
        ),
        roles = c(
          overview = "prediction_overview",
          estimates = "unit_estimates",
          row_review = "prediction_row_review",
          population_review = "prediction_population_review",
          settings = "scoring_settings",
          notes = "interpretation_notes"
        ),
        descriptions = c(
          overview = "Posterior-scoring overview for the current unit-prediction run.",
          estimates = "Posterior summaries for the scored persons.",
          row_review = "Row-level preparation review for the supplied scoring data.",
          population_review = "Optional person-level omission review for latent-regression scoring.",
          settings = "Scoring settings carried into posterior unit prediction.",
          notes = "Compact interpretation notes for posterior scoring output."
        )
      )
    },
    "summary.mfrm_plausible_values" = {
      draw_tbl <- summary_table_bundle_df(summary_obj$draw_summary)
      estimate_tbl <- summary_table_bundle_df(summary_obj$estimates)
      row_review_tbl <- summary_table_bundle_df(summary_obj$row_review)
      population_review_tbl <- summary_table_bundle_df(summary_obj$population_review)
      settings_tbl <- summary_table_bundle_settings_df(summary_obj$settings)
      notes_tbl <- summary_table_bundle_text_df(summary_obj$notes, column = "Note")
      total_draws <- if ("Draws" %in% names(draw_tbl)) sum(draw_tbl$Draws, na.rm = TRUE) else nrow(draw_tbl)
      overview_tbl <- data.frame(
        Persons = nrow(draw_tbl),
        TotalDraws = total_draws,
        EstimateRows = nrow(estimate_tbl),
        RowReviewRows = nrow(row_review_tbl),
        PopulationReviewRows = nrow(population_review_tbl),
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
          row_review = row_review_tbl,
          population_review = population_review_tbl,
          settings = settings_tbl,
          notes = notes_tbl
        ),
        roles = c(
          overview = "plausible_value_overview",
          draw_summary = "plausible_value_draws",
          estimates = "unit_estimates",
          row_review = "prediction_row_review",
          population_review = "prediction_population_review",
          settings = "scoring_settings",
          notes = "interpretation_notes"
        ),
        descriptions = c(
          overview = "Approximate plausible-value overview for the current posterior scoring run.",
          draw_summary = "Empirical summaries of the sampled posterior draws by person.",
          estimates = "Companion posterior EAP summaries paired with the draw summary.",
          row_review = "Row-level preparation review for the supplied scoring data.",
          population_review = "Optional person-level omission review for latent-regression scoring.",
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
#' @param x An `mfrm_fit`, `mfrm_diagnostics`, `mfrm_precision_review`,
#'   `mfrm_fit_measures`, `mfrm_facets_fit_review`,
#'   `mfrm_person_fit_indices`, `mfrm_data_description`,
#'   `mfrm_reporting_checklist`, `mfrm_apa_outputs`,
#'   `mfrm_design_evaluation`, `mfrm_signal_detection`,
#'   `mfrm_recovery_simulation`, `mfrm_recovery_assessment`,
#'   `mfrm_population_prediction`, `mfrm_future_branch_active_branch`,
#'   `mfrm_facets_run`, `mfrm_bias`, `mfrm_anchor_review`,
#'   `mfrm_linking_review`, `mfrm_misfit_casebook`, `mfrm_weighting_review`,
#'   `mfrm_unit_prediction`, or `mfrm_plausible_values` object, one of their
#'   `summary()` outputs, or a `summary.mfrmr_recovery_validation` object from
#'   the packaged validation protocol.
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
#' - [precision_review_report()] or `summary(precision_review)`
#' - [fit_measures_table()] or `summary(fit_measures)`
#' - [facets_fit_review()] or `summary(facets_fit_review)`
#' - [compute_person_fit_indices()] or `summary(person_fit)`
#' - [describe_mfrm_data()] or `summary(ds)`
#' - [reporting_checklist()] or `summary(chk)`
#' - [build_apa_outputs()] or `summary(apa)`
#' - [evaluate_mfrm_design()] or `summary(sim_eval)`
#' - [evaluate_mfrm_signal_detection()] or `summary(sig_eval)`
#' - [evaluate_mfrm_recovery()] or `summary(rec)`
#' - [assess_mfrm_recovery()] or `summary(rec_assessment)`
#' - `summary(validation)` from `recovery-validation.R`
#' - [predict_mfrm_population()] or `summary(pred)`
#' - `planning_schema$future_branch_active_branch` or `summary(...)`
#' - [run_mfrm_facets()] or `summary(out)`
#' - [estimate_bias()] or `summary(bias)`
#' - [review_mfrm_anchors()] or `summary(review)`
#' - [build_linking_review()] or `summary(review)`
#' - [build_misfit_casebook()] or `summary(casebook)`
#' - [build_weighting_review()] or `summary(review)`
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
#' - recovery-assessment and recovery-validation summaries expose
#'   `diagnostic_reporting_notes` before `diagnostic_review` or
#'   `diagnostic_oc_summary` so fit/separation caveats can be reported without
#'   treating them as recovery or release gates.
#' - recovery-validation summaries expose `condition_reporting_notes` before
#'   `condition_summary` so GPCM generator stress and sparse score support are
#'   not mistaken for recovery-metric failures.
#' - precision-review summaries expose `fit_separation_basis` so fit,
#'   ZSTD, separation/reliability/strata, and QC thresholds remain separate
#'   reporting surfaces rather than implicit validation gates.
#' - fit-measure and FACETS fit-review summaries expose df/ZSTD sensitivity
#'   tables under precision-review roles, keeping MnSq status, ZSTD
#'   standardization, and external FACETS matching distinct in appendix
#'   handoffs.
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
#' 5. For recovery-assessment or recovery-validation summaries, inspect
#'    `bundle$tables$reading_order` first when it is available.
#' 6. For recovery-assessment or recovery-validation summaries with retained
#'    diagnostics, read `diagnostic_reporting_notes` before the raw
#'    `diagnostic_review` or `diagnostic_oc_summary`. Read
#'    `condition_reporting_notes` before `condition_review` or
#'    `condition_summary` when bounded `GPCM` generator stress is part of the
#'    plan.
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
#'                 method = "JML", maxit = 30)
#' bundle <- build_summary_table_bundle(fit)
#' bundle$table_index
#' summary(bundle)$role_summary
#' }
#'
#' # Recovery-validation output can be converted to appendix-ready tables.
#' \dontrun{
#' source(system.file("validation", "recovery-validation.R", package = "mfrmr"))
#' validation <- mfrmr_run_recovery_validation(
#'   case_ids = c("gpcm_slope_profile", "gpcm_high_dispersion_sparse"),
#'   quick = TRUE,
#'   seed = 20260525
#' )
#' validation_bundle <- build_summary_table_bundle(summary(validation))
#' validation_bundle$tables$reading_order
#' validation_bundle$tables$topline_release_decision
#' validation_bundle$tables$condition_reporting_notes
#' validation_bundle$tables$condition_summary
#' validation_bundle$tables$diagnostic_reporting_notes
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
      "precision_review",
      "facet_precision",
      "flag_counts",
      "missingness",
      "score_usage",
      "facet_coverage",
      "agreement",
      "checklist_overview",
      "section_coverage",
      "facets_relationship_wording",
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
      "prediction_row_review",
      "prediction_population_review",
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
      "anchor_review_overview",
      "anchor_issue_counts",
      "level_observation_review",
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
      "methods",
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
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
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
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
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
      180L, 190L, 195L, 200L, 210L, 220L,
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
      "Include to review latent-regression design-matrix columns and variance screening.",
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
      "Core FACETS-positioning wording for methods or migration appendices.",
      "Core priority distribution for reporting follow-up.",
      "Core manuscript-draft coverage overview.",
      "Core APA component inventory.",
      "Draft QA surface; keep out of recommended presets.",
      "Preview-only draft text; keep out of recommended presets.",
      "Bridge metadata, useful for workflow but not manuscript appendix.",
      "Exploratory extreme table; available only in full exports.",
      "Exploratory extreme table; available only in full exports.",
      "Exploratory extreme table; available only in full exports.",
      "Exploratory extreme table; available only in full exports.",
      "Drafting action list; keep out of recommended presets.",
      "Checklist settings; keep out of recommended presets.",
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
      "Row-level scoring review; retain for full exports but omit from compact presets.",
      "Recommended latent-regression scoring omission review when population-model scoring is active.",
      "Methods/settings appendix table for posterior scoring inputs; recommended but not compact.",
      "Recommended overview table for plausible-value appendix handoff.",
      "Recommended plausible-value draw summary table for posterior-scoring appendices.",
      "Recommended design-performance table for simulation design appendices.",
      "Recommended signal-detection table for simulation diagnostics appendices.",
      "Recommended design-grid table documenting requested forecast inputs.",
      "Recommended forecast-summary table for population prediction appendices.",
      "Workflow-only overview for one-shot run handoff; omit from manuscript presets.",
      "Workflow-only column mapping for replay and traceability handoff.",
      "Workflow-only run settings for one-shot workflow provenance.",
      "Recommended overview table for bias-screening appendix handoff.",
      "Recommended fixed-effect chi-square table for bias-screening appendices.",
      "Bias iteration status table; retain for full exports but omit from compact presets.",
      "Recommended ranked bias-screening row table for immediate follow-up.",
      "Recommended overview table for anchor-review appendix handoff.",
      "Recommended anchor issue-count table for pre-fit review appendices.",
      "Recommended level-observation review table; retain for full appendices but omit from compact presets.",
      "Recommended score-category usage table for anchor-review support review.",
      "Recommended pre-fit anchor-risk table for linking-review appendices.",
      "Recommended drift-risk table for linking-review appendices.",
      "Recommended chain-risk table for linking-review appendices.",
      "Methods/settings table for operational review helpers; recommended but not compact."
    ),
    stringsAsFactors = FALSE
  )
  recovery_roles <- data.frame(
    Role = c(
      "recovery_overview",
      "recovery_design_basis",
      "recovery_settings",
      "recovery_performance",
      "recovery_replications",
      "recovery_diagnostic_operating_characteristics",
      "recovery_diagnostic_oc_summary",
      "recovery_assessment_overview",
      "recovery_assessment_reading_order",
      "recovery_assessment_checklist",
      "recovery_condition_reporting_notes",
      "recovery_condition_review",
      "recovery_diagnostic_reporting_notes",
      "recovery_diagnostic_review",
      "recovery_metric_review",
      "recovery_uncertainty_review"
    ),
    AppendixSection = c(
      "methods",
      "methods",
      "methods",
      "results",
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
      "diagnostics"
    ),
    RecommendedAppendix = rep(TRUE, 16),
    CompactAppendix = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    PreferredAppendixOrder = 260:275,
    AppendixRationale = c(
      "Recommended overview table for parameter-recovery simulation appendix handoff.",
      "Recommended ADEMP-style methods table documenting the recovery simulation basis.",
      "Methods/settings table for recovery simulation provenance; recommended but not compact.",
      "Recommended parameter-group recovery-performance table for simulation appendices.",
      "Recommended replication-status table for recovery simulation diagnostics.",
      "Replication-by-facet fit/separation operating characteristics; retain for diagnostics exports.",
      "Recommended facet-level fit/separation operating-characteristic summary for diagnostic context.",
      "Recommended overview table for recovery adequacy assessment appendix handoff.",
      "Recommended reading-order table for recovery assessment handoff.",
      "Recommended checklist table for reviewer-facing recovery adequacy decisions.",
      "Recommended reporter-facing table for generator-condition caveats kept separate from recovery metrics.",
      "Recommended generator-condition table for interpreting GPCM slope-regime recovery evidence.",
      "Recommended reporter-facing table for fit/separation diagnostic caveats kept separate from recovery gates.",
      "Recommended diagnostic-only fit/separation review table for recovery assessment handoff.",
      "Recommended parameter-group review table for recovery adequacy follow-up.",
      "Recommended uncertainty-evidence table separating coverage availability, SE availability, and coverage decision status."
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, recovery_roles)
  recovery_validation_roles <- data.frame(
    Role = c(
      "recovery_validation_topline",
      "recovery_validation_reading_order",
      "recovery_validation_release_decisions",
      "recovery_validation_case_summary",
      "recovery_validation_condition_reporting_notes",
      "recovery_validation_condition_summary",
      "recovery_validation_diagnostic_reporting_notes",
      "recovery_validation_diagnostic_oc_summary",
      "recovery_validation_domain_decisions"
    ),
    AppendixSection = c(
      "results",
      "results",
      "results",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics"
    ),
    RecommendedAppendix = rep(TRUE, 9),
    CompactAppendix = rep(TRUE, 9),
    PreferredAppendixOrder = 276:284,
    AppendixRationale = c(
      "Recommended top-line recovery-validation decision table for release-review appendices.",
      "Recommended reading-order table for recovery-validation handoff.",
      "Recommended case-level release-decision table for validation handoff.",
      "Recommended validation case summary table for release-review traceability.",
      "Recommended reporter-facing table for generator-condition caveats kept out of release gates.",
      "Recommended generator-condition summary table separating GPCM slope-regime and score-support stress.",
      "Recommended reporter-facing table for fit/separation diagnostic caveats kept out of release gates.",
      "Recommended fit/separation operating-characteristic summary kept separate from release-recovery gates.",
      "Recommended long-form domain-decision table for validation diagnostics."
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, recovery_validation_roles)
  sparse_roles <- data.frame(
    Role = c(
      "sparse_design_diagnostics",
      "recovery_sparse_design_diagnostics",
      "peer_review_design_diagnostics",
      "peer_review_load_summary",
      "peer_review_submission_load",
      "peer_review_reviewer_load",
      "peer_review_common_links",
      "peer_review_low_common_links",
      "peer_review_reciprocal_pairs"
    ),
    AppendixSection = c(
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics",
      "diagnostics"
    ),
    RecommendedAppendix = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    CompactAppendix = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
    PreferredAppendixOrder = c(259.5, 263.5, 264.5, 264.6, 264.7, 264.8, 264.9, 265.1, 265.2),
    AppendixRationale = c(
      "Recommended planned-missingness and rater-link diagnostics for sparse linked design appendices.",
      "Recommended replication-level sparse linked design diagnostics kept separate from recovery metrics.",
      "Recommended peer-review assignment and reviewer-link diagnostics kept separate from measurement estimates.",
      "Recommended compact load-balance summary for peer-review assignment appendices.",
      "Full submission-load table; retain for detailed diagnostics exports.",
      "Recommended reviewer-load table for assignment follow-up.",
      "Full reviewer-pair common-submission table; retain for detailed diagnostics exports.",
      "Recommended reviewer-pair low-common-link table for peer-review assignment follow-up.",
      "Recommended reciprocal review-pair table for peer-review assignment transparency."
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, sparse_roles)
  diagnostic_screening_roles <- data.frame(
    Role = c(
      "diagnostic_screening_design_basis",
      "diagnostic_screening_settings",
      "diagnostic_screening_overview",
      "diagnostic_screening_reading_order",
      "diagnostic_screening_next_actions",
      "diagnostic_screening_reporting_notes",
      "diagnostic_screening_figure_recipes",
      "diagnostic_screening_scenario_summary",
      "diagnostic_screening_performance",
      "diagnostic_screening_report_signals",
      "diagnostic_screening_contrast",
      "diagnostic_screening_plot_data",
      "diagnostic_screening_runtime"
    ),
    AppendixSection = c(
      "methods",
      "methods",
      "methods",
      "reporting",
      "workflow",
      "reporting",
      "workflow",
      "diagnostics",
      "diagnostics",
      "reporting",
      "diagnostics",
      "diagnostics",
      "methods"
    ),
    RecommendedAppendix = rep(TRUE, 13),
    CompactAppendix = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    PreferredAppendixOrder = 291:303,
    AppendixRationale = c(
      "Recommended ADEMP-style methods table documenting the diagnostic-screening simulation basis.",
      "Methods/settings table for diagnostic-screening provenance; recommended but not compact.",
      "Recommended overview table for diagnostic-screening appendix handoff.",
      "Recommended reading-order table for diagnostic-screening summaries and appendices.",
      "Recommended action table for diagnostic-screening follow-up and export routing.",
      "Recommended reporting-boundary table for cautious diagnostic-screening interpretation.",
      "Recommended figure/display recipe table so custom graphics keep plot calls, captions, and boundaries explicit.",
      "Recommended scenario-by-design screening summary for diagnostic appendices.",
      "Recommended operating-characteristic table for diagnostic-screening interpretation.",
      "Recommended report-signal table when mfrm_report() review signals were retained.",
      "Recommended misspecification-minus-baseline contrast table for diagnostic follow-up.",
      "Long-form draw-free plot data for interactive or custom visualization handoff; recommended but not compact.",
      "Runtime operating-characteristic table; recommended for methods appendices but not compact."
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, diagnostic_screening_roles)
  network_review_roles <- data.frame(
    Role = c(
      "network_review_overview",
      "network_design_summary",
      "network_facet_vulnerability",
      "network_central_nodes",
      "network_articulation_nodes",
      "network_bridge_edges"
    ),
    AppendixSection = rep("diagnostics", 6),
    RecommendedAppendix = rep(TRUE, 6),
    CompactAppendix = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    PreferredAppendixOrder = 285:290,
    AppendixRationale = c(
      "Recommended overview table for design-network connectedness and vulnerability review.",
      "Recommended graph-level design-network summary for connectivity diagnostics.",
      "Recommended facet-level design-network vulnerability table for follow-up appendices.",
      "High-betweenness node table; retain for full diagnostics exports but omit from compact presets.",
      "Recommended articulation-point table for fragile-link follow-up.",
      "Recommended bridge-edge table for one-link dependency follow-up."
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, network_review_roles)
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
#'                 method = "JML", maxit = 30)
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
#'                 method = "JML", maxit = 30)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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

#' Convert an `apa_table` to a `knitr::kable()` object
#'
#' Renders the table data for direct inclusion in RMarkdown,
#' Quarto, or HTML reports, wiring the `caption` and `note` slots
#' into the standard APA placement (caption above, note below).
#' When `kableExtra` is installed the note is attached as a footer;
#' otherwise the note is appended as a `knitr::asis_output()` block.
#'
#' @param x An `apa_table` object from [apa_table()].
#' @param format One of `"pipe"` (default, Markdown), `"html"`, or
#'   `"latex"`, passed through to `knitr::kable()`.
#' @param digits Numeric; passed to `knitr::kable()`.
#' @param ... Additional arguments forwarded to `knitr::kable()`.
#'
#' @return A `knitr_kable` object ready to be printed inline in a
#'   report, or a message when `knitr` is unavailable.
#' @seealso [as_flextable.apa_table()], [apa_table()].
#' @export
as_kable.apa_table <- function(x, format = c("pipe", "html", "latex"),
                               digits = 3L, ...) {
  format <- match.arg(format)
  if (!requireNamespace("knitr", quietly = TRUE)) {
    message("`as_kable.apa_table()` requires the `knitr` package (in Suggests).")
    return(invisible(NULL))
  }
  tbl <- if (is.data.frame(x$table)) x$table else as.data.frame(x$table %||% list())
  caption <- as.character(x$caption %||% "")
  note <- as.character(x$note %||% "")
  k <- knitr::kable(tbl, format = format, digits = digits,
                    caption = if (nzchar(caption)) caption else NULL, ...)
  if (nzchar(note)) {
    # `kableExtra::footnote()` internally converts the kable to HTML, so
    # only route through it when the user actually wants HTML or LaTeX.
    # For Markdown / "pipe" output we fall back to the safe append path;
    # otherwise a user asking for "pipe" would silently get an HTML
    # table, which then breaks Quarto / RMarkdown paragraph-mode paste.
    use_kableextra <- format %in% c("html", "latex") &&
      requireNamespace("kableExtra", quietly = TRUE)
    if (use_kableextra) {
      k <- kableExtra::footnote(k, general = note,
                                general_title = "Note.",
                                footnote_as_chunk = TRUE)
    } else {
      k <- paste0(k, "\n\nNote. ", note)
      class(k) <- c("knitr_kable", class(k))
    }
  }
  k
}

#' Convert an `apa_table` to a `flextable`
#'
#' Produces a Word / PowerPoint-friendly `flextable` with the
#' caption and note wired in. Requires `flextable` (in Suggests).
#'
#' @param x An `apa_table` object from [apa_table()].
#' @param ... Additional arguments reserved for future use.
#'
#' @return A `flextable` object, or a message when `flextable` is
#'   unavailable.
#' @seealso [as_kable.apa_table()], [apa_table()].
#' @export
as_flextable.apa_table <- function(x, ...) {
  if (!requireNamespace("flextable", quietly = TRUE)) {
    message("`as_flextable.apa_table()` requires the `flextable` package (in Suggests).")
    return(invisible(NULL))
  }
  tbl <- if (is.data.frame(x$table)) x$table else as.data.frame(x$table %||% list())
  caption <- as.character(x$caption %||% "")
  note <- as.character(x$note %||% "")
  ft <- flextable::flextable(tbl)
  if (nzchar(caption)) {
    ft <- flextable::set_caption(ft, caption)
  }
  if (nzchar(note)) {
    ft <- flextable::add_footer_lines(ft, values = paste0("Note. ", note))
  }
  ft
}

#' Generic for converting objects to a `knitr::kable`
#'
#' @param x Object to convert.
#' @param ... Passed to methods.
#'
#' @return A `knitr::kable` object (concrete return type from the
#'   underlying method, e.g. `[as_kable.apa_table()]` returns a
#'   `kableExtra` object when the package is installed).
#'
#' @seealso [as_kable.apa_table()] for the `apa_table` method;
#'   [as_flextable()] for a `flextable`-targeted alternative;
#'   [apa_table()] for constructing an `apa_table` in the first place.
#' @export
as_kable <- function(x, ...) UseMethod("as_kable")

#' Generic for converting objects to a `flextable`
#'
#' @param x Object to convert.
#' @param ... Passed to methods.
#'
#' @return A `flextable` object (concrete return type from the
#'   underlying method, e.g. `[as_flextable.apa_table()]` returns a
#'   `flextable` ready for `flextable::save_as_docx()`).
#'
#' @seealso [as_flextable.apa_table()] for the `apa_table` method;
#'   [as_kable()] for a `knitr::kable`-targeted alternative;
#'   [apa_table()] for constructing an `apa_table` in the first place.
#' @export
as_flextable <- function(x, ...) UseMethod("as_flextable")

#' Summarize an APA/FACETS table object
#'
#' @param object Output from [apa_table()].
#' @param digits Number of digits used for numeric summaries.
#' @param top_n Maximum numeric columns shown in `numeric_profile`.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Compact summary helper for QA of table data before manuscript export.
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
#' For bounded `GPCM`, this helper returns caveated warning/summary maps over
#' supported diagnostics, direct tables, and plots. The returned object includes
#' `gpcm_boundary` so score-side, design-forecasting, DFF, and linking routes
#' remain visibly separate capability rows.
#'
#' @section Interpreting output:
#' - `warning_map`: rule-triggered warning text by visual key.
#' - `summary_map`: descriptive narrative text by visual key.
#' - strict marginal keys appear when `diagnose_mfrm(..., diagnostic_mode = "both")`
#'   supplies latent-integrated first-order and pairwise screening summaries.
#' - `warning_counts` / `summary_counts`: message-count tables for QA checks.
#' - `plot_payloads`: ready-to-reuse `mfrm_plot_data` objects for the bundle's
#'   own comparison/count plots and, when step estimates are available, the
#'   exploratory `category_probability_surface` data from
#'   `plot(fit, type = "ccc_surface", draw = FALSE)`. The surface data carry
#'   `category_support`, `interpretation_guide`, and `reporting_policy` tables
#'   for zero-frequency category and reporting-boundary checks.
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
#' - `plot_payloads`: reusable draw-free `mfrm_plot_data` objects for
#'   `comparison`, `warning_counts`, `summary_counts`, and optionally
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
#'   method = "MML", model = "RSM", quad_points = 7, maxit = 30
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
      "No direct FACETS equivalent (exploratory category-probability surface data)"
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "build_visual_summaries()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review",
        "Design planning and forecasting"
      )
    ),
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

build_contract_metric_review <- function(outputs, tol = 1e-8) {
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

fit_review_normalize_name <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", as.character(x)))
}

fit_review_find_col <- function(df, candidates, explicit = NULL, required = FALSE,
                               label = NULL, data_label = "facets_fit") {
  if (!is.null(explicit)) {
    explicit <- as.character(explicit[1])
    if (explicit %in% names(df)) return(explicit)
    stop("Column `", explicit, "` was not found in `", data_label, "`.", call. = FALSE)
  }
  nm <- names(df)
  norm <- fit_review_normalize_name(nm)
  hit <- match(fit_review_normalize_name(candidates), norm)
  hit <- hit[is.finite(hit) & !is.na(hit)]
  if (length(hit) > 0L) return(nm[hit[1]])
  if (isTRUE(required)) {
    stop(
      "Could not infer the ", label %||% "required", " column in `", data_label, "`. ",
      "Available columns: ", paste(nm, collapse = ", "), ".",
      call. = FALSE
    )
  }
  NA_character_
}

fit_review_clean_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub("[<>]", "", x)
  x[x %in% c("", ".", "*", "NA", "N/A", "NaN", "Inf", "-Inf")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

fit_review_numeric_col <- function(df, col) {
  if (is.na(col) || !nzchar(col) || !col %in% names(df)) {
    return(rep(NA_real_, nrow(df)))
  }
  fit_review_clean_numeric(df[[col]])
}

fit_review_pmax_na <- function(...) {
  out <- pmax(..., na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

fit_review_prepare_facet_map <- function(facet_map = NULL) {
  if (is.null(facet_map)) {
    return(c("1" = "Person"))
  }
  facet_map <- as.character(facet_map)
  if (is.null(names(facet_map)) || any(!nzchar(names(facet_map)))) {
    names(facet_map) <- as.character(seq_along(facet_map))
  }
  facet_map
}

fit_review_scorefile_number <- function(path) {
  base <- tolower(basename(path))
  hit <- regmatches(base, regexec("^score[._-]([0-9]+)\\.txt$", base))[[1]]
  if (length(hit) >= 2L) hit[2] else NA_character_
}

fit_review_resolve_scorefile_facet <- function(path, facet = NULL, facet_map = NULL) {
  if (!is.null(facet)) {
    facet <- as.character(facet[1])
    if (nzchar(facet)) return(facet)
  }
  facet_num <- fit_review_scorefile_number(path)
  fmap <- fit_review_prepare_facet_map(facet_map)
  if (!is.na(facet_num) && facet_num %in% names(fmap)) {
    return(unname(fmap[[facet_num]]))
  }
  if (!is.na(facet_num)) {
    return(paste0("Facet", facet_num))
  }
  stop(
    "Could not infer the facet name for FACETS score file `", path, "`. ",
    "Supply `facet` or a named `facet_map`.",
    call. = FALSE
  )
}

fit_review_read_delimited_file <- function(path, delimiter = NULL, encoding = "UTF-8") {
  lines <- readLines(path, warn = FALSE, encoding = encoding)
  first <- lines[nzchar(trimws(lines))][1] %||% ""
  if (is.null(delimiter)) {
    count_delim <- function(pattern) {
      hit <- gregexpr(pattern, first, fixed = TRUE)[[1]]
      if (length(hit) == 1L && hit[1] == -1L) 0L else length(hit)
    }
    comma_n <- count_delim(",")
    tab_n <- count_delim("\t")
    semicolon_n <- count_delim(";")
    delimiter <- if (tab_n > comma_n) {
      "\t"
    } else if (semicolon_n > comma_n) {
      ";"
    } else {
      ","
    }
  }
  utils::read.table(
    file = path,
    sep = delimiter,
    header = TRUE,
    quote = "\"",
    comment.char = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", ".", "*")
  )
}

fit_review_repair_csv_lines <- function(lines) {
  vapply(lines, function(line) {
    if ((stringr::str_count(line, "\"") %% 2L) == 1L) {
      paste0(line, "\"")
    } else {
      line
    }
  }, character(1))
}

fit_review_scorefile_fixed_data <- function(path, raw = NULL, encoding = "UTF-8") {
  if (is.null(raw)) {
    raw <- readLines(path, warn = FALSE, encoding = encoding)
  }
  raw <- raw[nzchar(trimws(raw))]
  if (length(raw) == 0L) {
    stop("Could not find fixed-field FACETS score-file rows in `", path, "`.", call. = FALSE)
  }
  field_num <- function(lines, start, end) {
    vals <- vapply(lines, function(line) {
      if (nchar(line, type = "chars") < start) {
        ""
      } else {
        substr(line, start, min(end, nchar(line, type = "chars")))
      }
    }, character(1))
    fit_review_clean_numeric(vals)
  }
  field_chr <- function(lines, start, end = NA_integer_) {
    vapply(lines, function(line) {
      n <- nchar(line, type = "chars")
      if (n < start) {
        ""
      } else {
        end_i <- if (is.na(end)) n else min(end, n)
        trimws(substr(line, start, end_i))
      }
    }, character(1))
  }

  raw_score <- field_num(raw, 1, 10)
  t_count <- field_num(raw, 11, 20)
  measure <- field_num(raw, 41, 50)
  se <- field_num(raw, 51, 60)
  infit <- field_num(raw, 61, 70)
  infit_z <- field_num(raw, 71, 80)
  outfit <- field_num(raw, 81, 90)
  outfit_z <- field_num(raw, 91, 100)
  infit_df <- field_num(raw, 191, 200)
  outfit_df <- field_num(raw, 221, 230)
  element_number <- field_chr(raw, 241, 250)
  element_label <- field_chr(raw, 251, NA_integer_)
  element_label[!nzchar(element_label)] <- element_number[!nzchar(element_label)]
  keep <- nzchar(element_label) & (
    is.finite(t_count) | is.finite(measure) | is.finite(infit) |
      is.finite(infit_z) | is.finite(outfit) | is.finite(outfit_z)
  )
  if (!any(keep)) {
    stop("Could not find fixed-field FACETS score-file rows in `", path, "`.", call. = FALSE)
  }
  data.frame(
    Level = element_label[keep],
    T.Score = raw_score[keep],
    T.Count = t_count[keep],
    Measure = measure[keep],
    S.E. = se[keep],
    InfitMS = infit[keep],
    InfitZ = infit_z[keep],
    OutfitMS = outfit[keep],
    OutfitZ = outfit_z[keep],
    InfitDF = infit_df[keep],
    OutfitDF = outfit_df[keep],
    ElementNumber = element_number[keep],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

fit_review_scorefile_data <- function(path, encoding = "UTF-8") {
  raw <- readLines(path, warn = FALSE, encoding = encoding)
  raw <- raw[nzchar(trimws(raw))]
  header_idx <- which(grepl("Measure", raw, fixed = TRUE) & grepl(",", raw, fixed = TRUE))[1]
  if (is.na(header_idx)) {
    return(fit_review_scorefile_fixed_data(path, raw = raw, encoding = encoding))
  }
  csv_lines <- raw[seq.int(header_idx, length(raw))]
  csv_lines <- csv_lines[grepl(",", csv_lines, fixed = TRUE)]
  csv_lines <- fit_review_repair_csv_lines(csv_lines)
  utils::read.table(
    text = paste(csv_lines, collapse = "\n"),
    sep = ",",
    header = TRUE,
    quote = "\"",
    comment.char = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", ".", "*")
  )
}

fit_review_standardize_frame <- function(df,
                                        facet = NULL,
                                        facet_col = NULL,
                                        level_col = NULL,
                                        source = "facets_fit_table",
                                        data_label = "FACETS fit table") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) {
    return(tibble::tibble(
      Facet = character(0), Level = character(0), Estimate = numeric(0),
      SE = numeric(0), N = numeric(0), Infit = numeric(0), Outfit = numeric(0),
      InfitZSTD = numeric(0), OutfitZSTD = numeric(0),
      DF_Infit = numeric(0), DF_Outfit = numeric(0), Source = character(0)
    ))
  }

  person_col <- fit_review_find_col(
    df, c("Person", "PersonID", "PersonId"),
    required = FALSE,
    data_label = data_label
  )
  facet_col <- fit_review_find_col(
    df,
    c("Facet", "FacetName", "Facets"),
    explicit = facet_col,
    required = is.null(facet) && is.na(person_col),
    label = "facet",
    data_label = data_label
  )
  level_col <- fit_review_find_col(
    df,
    c("Level", "Element", "ElementName", "Name", "Label", "Person", "Rater",
      "Item", "Task", "Criterion", "Criteria"),
    explicit = level_col,
    required = is.na(person_col),
    label = "level",
    data_label = data_label
  )

  if (!is.null(facet)) {
    facet_values <- rep(as.character(facet[1]), nrow(df))
    level_values <- as.character(df[[level_col]])
  } else if (!is.na(person_col) && (is.na(facet_col) || is.na(level_col))) {
    facet_values <- rep("Person", nrow(df))
    level_values <- as.character(df[[person_col]])
  } else {
    facet_values <- as.character(df[[facet_col]])
    level_values <- as.character(df[[level_col]])
  }

  estimate_col <- fit_review_find_col(
    df, c("Estimate", "Measure", "Logit", "FACETSMeasure", "FACETS_Measure"),
    required = FALSE, data_label = data_label
  )
  se_col <- fit_review_find_col(
    df, c("SE", "S.E.", "StdError", "StandardError", "ModelSE"),
    required = FALSE, data_label = data_label
  )
  infit_col <- fit_review_find_col(
    df, c("Infit", "InfitMS", "InfitMnSq", "InfitMNSQ", "InfitMeanSquare",
          "InfitMSQ", "FACETS_Infit"),
    required = FALSE, data_label = data_label
  )
  outfit_col <- fit_review_find_col(
    df, c("Outfit", "OutfitMS", "OutfitMnSq", "OutfitMNSQ", "OutfitMeanSquare",
          "OutfitMSQ", "FACETS_Outfit"),
    required = FALSE, data_label = data_label
  )
  infit_z_col <- fit_review_find_col(
    df, c("InfitZSTD", "InfitZstd", "InfitZ", "ZSTDInfit", "ZstdInfit",
          "FACETS_InfitZSTD"),
    required = FALSE, data_label = data_label
  )
  outfit_z_col <- fit_review_find_col(
    df, c("OutfitZSTD", "OutfitZstd", "OutfitZ", "ZSTDOutfit", "ZstdOutfit",
          "FACETS_OutfitZSTD"),
    required = FALSE, data_label = data_label
  )
  df_infit_col <- fit_review_find_col(
    df, c("DF_Infit", "DFInfit", "InfitDF", "InfitDf", "Infitdf",
          "InfitDegreesFreedom", "FACETS_DF_Infit"),
    required = FALSE, data_label = data_label
  )
  df_outfit_col <- fit_review_find_col(
    df, c("DF_Outfit", "DFOutfit", "OutfitDF", "OutfitDf", "Outfitdf",
          "OutfitDegreesFreedom", "FACETS_DF_Outfit"),
    required = FALSE, data_label = data_label
  )
  n_col <- fit_review_find_col(
    df, c("N", "Count", "TCount", "T.Count", "TotalCount", "Observations",
          "FACETS_N"),
    required = FALSE, data_label = data_label
  )

  tibble::tibble(
    Facet = facet_values,
    Level = level_values,
    Estimate = fit_review_numeric_col(df, estimate_col),
    SE = fit_review_numeric_col(df, se_col),
    N = fit_review_numeric_col(df, n_col),
    Infit = fit_review_numeric_col(df, infit_col),
    Outfit = fit_review_numeric_col(df, outfit_col),
    InfitZSTD = fit_review_numeric_col(df, infit_z_col),
    OutfitZSTD = fit_review_numeric_col(df, outfit_z_col),
    DF_Infit = fit_review_numeric_col(df, df_infit_col),
    DF_Outfit = fit_review_numeric_col(df, df_outfit_col),
    Source = as.character(source[1])
  ) |>
    dplyr::filter(!is.na(.data$Facet), !is.na(.data$Level),
                  nzchar(.data$Facet), nzchar(.data$Level))
}

fit_review_read_score_file <- function(path,
                                      facet = NULL,
                                      facet_map = NULL,
                                      level_col = NULL,
                                      encoding = "UTF-8") {
  df <- fit_review_scorefile_data(path, encoding = encoding)
  fnum_col <- fit_review_find_col(
    df, c("FNumber", "F-Number", "FacetNumber"),
    required = FALSE,
    data_label = "FACETS score file"
  )
  if (is.null(level_col) && !is.na(fnum_col)) {
    idx <- match(fnum_col, names(df))
    if (is.finite(idx) && idx > 1L) {
      level_col <- names(df)[idx - 1L]
    }
  }
  facet_name <- fit_review_resolve_scorefile_facet(
    path,
    facet = facet,
    facet_map = facet_map
  )
  out <- fit_review_standardize_frame(
    df,
    facet = facet_name,
    level_col = level_col,
    source = basename(path),
    data_label = "FACETS score file"
  )
  if (!is.na(fit_review_scorefile_number(path))) {
    out$RawFacetNumber <- rep(fit_review_scorefile_number(path), nrow(out))
  }
  out
}

fit_review_read_table_file <- function(path,
                                      facet = NULL,
                                      facet_col = NULL,
                                      level_col = NULL,
                                      delimiter = NULL,
                                      encoding = "UTF-8") {
  df <- fit_review_read_delimited_file(path, delimiter = delimiter, encoding = encoding)
  fit_review_standardize_frame(
    df,
    facet = facet,
    facet_col = facet_col,
    level_col = level_col,
    source = basename(path),
    data_label = "FACETS fit table"
  )
}

fit_review_read_one_path <- function(path,
                                    facet = NULL,
                                    facet_map = NULL,
                                    format = c("auto", "delimited", "scorefile"),
                                    facet_col = NULL,
                                    level_col = NULL,
                                    delimiter = NULL,
                                    encoding = "UTF-8") {
  format <- match.arg(format)
  if (!file.exists(path)) {
    stop("FACETS fit table file was not found: `", path, "`.", call. = FALSE)
  }
  if (dir.exists(path)) {
    score_files <- list.files(path, pattern = "^score[._-][0-9]+\\.txt$",
                              full.names = TRUE, ignore.case = TRUE)
    if (length(score_files) == 0L) {
      stop("Directory `", path, "` does not contain FACETS score.N.txt files.", call. = FALSE)
    }
    return(dplyr::bind_rows(lapply(score_files, fit_review_read_one_path,
                                   facet = NULL, facet_map = facet_map,
                                   format = "scorefile", facet_col = facet_col,
                                   level_col = level_col, delimiter = delimiter,
                                   encoding = encoding)))
  }

  if (identical(format, "auto")) {
    lines <- readLines(path, warn = FALSE, encoding = encoding, n = 50L)
    has_score_header <- any(grepl("Measure", lines, fixed = TRUE) &
                              grepl(",", lines, fixed = TRUE) &
                              grepl("F[- .]?Number|FNumber|FacetNumber", lines,
                                    ignore.case = TRUE))
    is_score_name <- !is.na(fit_review_scorefile_number(path))
    format <- if (is_score_name || has_score_header) "scorefile" else "delimited"
  }
  if (identical(format, "scorefile")) {
    fit_review_read_score_file(
      path,
      facet = facet,
      facet_map = facet_map,
      level_col = level_col,
      encoding = encoding
    )
  } else {
    fit_review_read_table_file(
      path,
      facet = facet,
      facet_col = facet_col,
      level_col = level_col,
      delimiter = delimiter,
      encoding = encoding
    )
  }
}

#' Read a FACETS fit table for fit review
#'
#' @param file Path to a FACETS-derived fit table. A character vector of files is
#'   accepted. A directory containing `score.N.txt` files is also accepted.
#' @param facet Optional facet name to assign when the file does not contain a
#'   facet column. Use this for one-facet CSV exports.
#' @param facet_map Optional character vector mapping FACETS score-file numbers
#'   to facet names, for example `c("1" = "Person", "2" = "Rater")`. If unnamed,
#'   positions are used as score-file numbers.
#' @param format File format. `"auto"` detects FACETS `score.N.txt` files from
#'   their name/header; `"delimited"` reads a CSV/TSV/semicolon table; and
#'   `"scorefile"` reads a FACETS score-file table.
#' @param facet_col,level_col Optional explicit column names for delimited
#'   tables when automatic detection is not sufficient. For score files,
#'   `level_col` can override the column immediately before `F-Number`.
#' @param delimiter Optional delimiter for delimited tables. If omitted, comma,
#'   tab, and semicolon are detected from the header line.
#' @param encoding File encoding passed to [readLines()].
#'
#' @details
#' This helper does not run FACETS. It reads FACETS output that already exists
#' on disk and normalizes it to columns that [facets_fit_review()] can consume:
#' `Facet`, `Level`, `Estimate`, `SE`, `N`, `Infit`, `Outfit`, `InfitZSTD`,
#' `OutfitZSTD`, `DF_Infit`, and `DF_Outfit`.
#'
#' Two common workflows are supported:
#' - a FACETS score file such as `score.2.txt`, where the facet name is supplied
#'   by `facet_map` or inferred as `Facet2`. Both comma-delimited score files
#'   with field names and fixed-field score files using the FACETS manual
#'   column positions are supported;
#' - a CSV/TSV table already exported from FACETS or a harmonization script,
#'   with FACETS-style column names such as `Infit MnSq`, `Outfit ZStd`, and
#'   `T.Count`.
#'
#' After import, pass the table to [facets_fit_review()]. Inspect
#' `review$external_table_quality` first when the FACETS export is partial,
#' duplicated, or missing MnSq/df columns. Then inspect
#' `review$external_comparison` for supplied FACETS-vs-mfrmr differences and
#' `review$df_sensitivity` / `review$df_sensitive` for engine-vs-FACETS-style
#' df/ZSTD convention sensitivity. Use `plot(review, type = "df_sensitivity")`
#' for a quick visual check of the largest ZSTD shifts caused by df convention.
#'
#' @return A tibble with standardized fit-table columns suitable for
#'   `facets_fit_review(fit, facets_fit = read_facets_fit_table(...))`.
#'
#' @seealso [facets_fit_review()], [diagnose_mfrm()]
#' @examples
#' path <- tempfile(fileext = ".csv")
#' write.csv(
#'   data.frame(
#'     Facet = "Rater", Level = "R1", Infit = 1.02, Outfit = 0.98,
#'     InfitZSTD = 0.3, OutfitZSTD = -0.2, DF_Infit = 12, DF_Outfit = 13
#'   ),
#'   path,
#'   row.names = FALSE
#' )
#' read_facets_fit_table(path)
#' @export
read_facets_fit_table <- function(file,
                                  facet = NULL,
                                  facet_map = NULL,
                                  format = c("auto", "delimited", "scorefile"),
                                  facet_col = NULL,
                                  level_col = NULL,
                                  delimiter = NULL,
                                  encoding = "UTF-8") {
  format <- match.arg(format)
  if (is.data.frame(file)) {
    return(fit_review_standardize_frame(
      file,
      facet = facet,
      facet_col = facet_col,
      level_col = level_col,
      source = "data_frame",
      data_label = "FACETS fit table"
    ))
  }
  paths <- as.character(file)
  if (length(paths) == 0L) {
    stop("`file` must contain at least one path.", call. = FALSE)
  }
  pieces <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    facet_i <- if (!is.null(facet) && length(facet) == length(paths)) facet[i] else facet
    pieces[[i]] <- fit_review_read_one_path(
      paths[i],
      facet = facet_i,
      facet_map = facet_map,
      format = format,
      facet_col = facet_col,
      level_col = level_col,
      delimiter = delimiter,
      encoding = encoding
    )
  }
  dplyr::bind_rows(pieces)
}

#' @rdname read_facets_fit_table
#' @export
import_facets_fit_table <- read_facets_fit_table

normalize_facets_fit_frame <- function(x,
                                       facet_col = NULL,
                                       level_col = NULL,
                                       source = "facets_fit") {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) {
    return(tibble::tibble(
      Source = character(0), Facet = character(0), Level = character(0),
      FACETS_Infit = numeric(0), FACETS_Outfit = numeric(0),
      FACETS_InfitZSTD = numeric(0), FACETS_OutfitZSTD = numeric(0),
      FACETS_DF_Infit = numeric(0), FACETS_DF_Outfit = numeric(0),
      FACETS_N = numeric(0)
    ))
  }

  person_col <- fit_review_find_col(df, c("Person", "PersonID", "PersonId"), required = FALSE)
  facet_col <- fit_review_find_col(
    df,
    c("Facet", "FacetName", "Facets"),
    explicit = facet_col,
    required = is.na(person_col),
    label = "facet"
  )
  level_col <- fit_review_find_col(
    df,
    c("Level", "Element", "ElementName", "Name", "Label"),
    explicit = level_col,
    required = is.na(person_col),
    label = "level"
  )

  if (!is.na(person_col) && (is.na(facet_col) || is.na(level_col))) {
    facet <- rep("Person", nrow(df))
    level <- as.character(df[[person_col]])
  } else {
    facet <- as.character(df[[facet_col]])
    level <- as.character(df[[level_col]])
  }

  infit_col <- fit_review_find_col(
    df,
    c("Infit", "InfitMS", "InfitMnSq", "InfitMNSQ", "InfitMeanSquare",
      "InfitMSQ", "FACETS_Infit")
  )
  outfit_col <- fit_review_find_col(
    df,
    c("Outfit", "OutfitMS", "OutfitMnSq", "OutfitMNSQ", "OutfitMeanSquare",
      "OutfitMSQ", "FACETS_Outfit")
  )
  infit_z_col <- fit_review_find_col(
    df,
    c("InfitZSTD", "InfitZstd", "InfitZ", "ZSTDInfit", "ZstdInfit",
      "FACETS_InfitZSTD")
  )
  outfit_z_col <- fit_review_find_col(
    df,
    c("OutfitZSTD", "OutfitZstd", "OutfitZ", "ZSTDOutfit", "ZstdOutfit",
      "FACETS_OutfitZSTD")
  )
  df_infit_col <- fit_review_find_col(
    df,
    c("DF_Infit", "DFInfit", "InfitDF", "InfitDf", "Infitdf",
      "InfitDegreesFreedom", "FACETS_DF_Infit")
  )
  df_outfit_col <- fit_review_find_col(
    df,
    c("DF_Outfit", "DFOutfit", "OutfitDF", "OutfitDf", "Outfitdf",
      "OutfitDegreesFreedom", "FACETS_DF_Outfit")
  )
  n_col <- fit_review_find_col(
    df,
    c("N", "Count", "TCount", "T.Count", "TotalCount", "Observations",
      "FACETS_N")
  )

  tibble::tibble(
    Source = as.character(source[1]),
    Facet = facet,
    Level = level,
    FACETS_Infit = fit_review_numeric_col(df, infit_col),
    FACETS_Outfit = fit_review_numeric_col(df, outfit_col),
    FACETS_InfitZSTD = fit_review_numeric_col(df, infit_z_col),
    FACETS_OutfitZSTD = fit_review_numeric_col(df, outfit_z_col),
    FACETS_DF_Infit = fit_review_numeric_col(df, df_infit_col),
    FACETS_DF_Outfit = fit_review_numeric_col(df, df_outfit_col),
    FACETS_N = fit_review_numeric_col(df, n_col)
  ) |>
    dplyr::filter(!is.na(.data$Facet), !is.na(.data$Level),
                  nzchar(.data$Facet), nzchar(.data$Level))
}

normalize_facets_fit_input <- function(facets_fit, facet_col = NULL, level_col = NULL) {
  if (is.null(facets_fit)) {
    return(tibble::tibble(
      Source = character(0), Facet = character(0), Level = character(0),
      FACETS_Infit = numeric(0), FACETS_Outfit = numeric(0),
      FACETS_InfitZSTD = numeric(0), FACETS_OutfitZSTD = numeric(0),
      FACETS_DF_Infit = numeric(0), FACETS_DF_Outfit = numeric(0),
      FACETS_N = numeric(0)
    ))
  }
  if (is.data.frame(facets_fit)) {
    return(normalize_facets_fit_frame(
      facets_fit,
      facet_col = facet_col,
      level_col = level_col,
      source = "facets_fit"
    ))
  }
  if (is.list(facets_fit)) {
    pieces <- list()
    nms <- names(facets_fit)
    if (is.null(nms)) nms <- paste0("facets_fit_", seq_along(facets_fit))
    for (i in seq_along(facets_fit)) {
      if (!is.data.frame(facets_fit[[i]])) next
      pieces[[length(pieces) + 1L]] <- normalize_facets_fit_frame(
        facets_fit[[i]],
        facet_col = facet_col,
        level_col = level_col,
        source = nms[i]
      )
    }
    if (length(pieces) > 0L) {
      return(dplyr::bind_rows(pieces))
    }
  }
  stop("`facets_fit` must be a data frame or a list of data frames.", call. = FALSE)
}

summarize_external_facets_fit_quality <- function(external_tbl) {
  external_tbl <- as.data.frame(external_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(external_tbl) == 0L) {
    return(data.frame(
      Rows = 0L,
      UniqueFacetLevelRows = 0L,
      DuplicateFacetLevelRows = 0L,
      CompleteMnSqRows = 0L,
      CompleteZSTDRows = 0L,
      CompleteDFRows = 0L,
      CompleteNRows = 0L,
      CompleteExternalFitRows = 0L,
      stringsAsFactors = FALSE
    ))
  }
  key <- paste(external_tbl$Facet, external_tbl$Level, sep = "\r")
  duplicated_key <- duplicated(key) | duplicated(key, fromLast = TRUE)
  finite_col <- function(nm) {
    if (!nm %in% names(external_tbl)) {
      return(rep(FALSE, nrow(external_tbl)))
    }
    is.finite(suppressWarnings(as.numeric(external_tbl[[nm]])))
  }
  has_infit <- finite_col("FACETS_Infit")
  has_outfit <- finite_col("FACETS_Outfit")
  has_infit_z <- finite_col("FACETS_InfitZSTD")
  has_outfit_z <- finite_col("FACETS_OutfitZSTD")
  has_df_infit <- finite_col("FACETS_DF_Infit")
  has_df_outfit <- finite_col("FACETS_DF_Outfit")
  has_n <- finite_col("FACETS_N")
  complete_mnsq <- has_infit & has_outfit
  complete_zstd <- has_infit_z & has_outfit_z
  complete_df <- has_df_infit & has_df_outfit
  data.frame(
    Rows = nrow(external_tbl),
    UniqueFacetLevelRows = length(unique(key)),
    DuplicateFacetLevelRows = sum(duplicated_key, na.rm = TRUE),
    CompleteMnSqRows = sum(complete_mnsq, na.rm = TRUE),
    CompleteZSTDRows = sum(complete_zstd, na.rm = TRUE),
    CompleteDFRows = sum(complete_df, na.rm = TRUE),
    CompleteNRows = sum(has_n, na.rm = TRUE),
    CompleteExternalFitRows = sum(complete_mnsq & complete_zstd & complete_df, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

facets_fit_review_prepare_diagnostics <- function(fit, diagnostics = NULL) {
  needs_diagnostics <- is.null(diagnostics) ||
    !is.list(diagnostics) ||
    is.null(diagnostics$fit) ||
    !all(c("InfitZSTD_FACETS", "OutfitZSTD_FACETS",
           "DF_Infit_FACETS", "DF_Outfit_FACETS") %in% names(diagnostics$fit))

  if (needs_diagnostics) {
    mode <- if (is.list(diagnostics)) {
      as.character(diagnostics$diagnostic_mode %||% "legacy")
    } else {
      "legacy"
    }
    if (!mode %in% c("legacy", "marginal_fit", "both")) mode <- "legacy"
    diagnostics <- diagnose_mfrm(
      fit,
      residual_pca = "none",
      diagnostic_mode = mode,
      fit_df_method = "both"
    )
  }
  diagnostics
}

build_internal_fit_standardization_review <- function(fit_tbl,
                                                     df_zstd_tolerance = 0.05,
                                                     df_zstd_large_shift = 0.5,
                                                     df_ratio_tolerance = 0.05) {
  fit_tbl <- as.data.frame(fit_tbl, stringsAsFactors = FALSE)
  required <- c("Facet", "Level", "Infit", "Outfit",
                "InfitZSTD", "OutfitZSTD",
                "DF_Infit", "DF_Outfit",
                "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
                "DF_Infit_FACETS", "DF_Outfit_FACETS")
  missing <- setdiff(required, names(fit_tbl))
  if (length(missing) > 0L) {
    stop("Diagnostics fit table is missing required FACETS comparison columns: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  df_infit_engine <- if ("DF_Infit_ENGINE" %in% names(fit_tbl)) fit_tbl$DF_Infit_ENGINE else fit_tbl$DF_Infit
  df_outfit_engine <- if ("DF_Outfit_ENGINE" %in% names(fit_tbl)) fit_tbl$DF_Outfit_ENGINE else fit_tbl$DF_Outfit
  z_infit_engine <- if ("InfitZSTD_ENGINE" %in% names(fit_tbl)) fit_tbl$InfitZSTD_ENGINE else fit_tbl$InfitZSTD
  z_outfit_engine <- if ("OutfitZSTD_ENGINE" %in% names(fit_tbl)) fit_tbl$OutfitZSTD_ENGINE else fit_tbl$OutfitZSTD

  source_tbl <- data.frame(
    Facet = as.character(fit_tbl$Facet),
    Level = as.character(fit_tbl$Level),
    N = suppressWarnings(as.numeric(fit_tbl$N %||% NA_real_)),
    Infit = suppressWarnings(as.numeric(fit_tbl$Infit)),
    Outfit = suppressWarnings(as.numeric(fit_tbl$Outfit)),
    DF_Infit_ENGINE = suppressWarnings(as.numeric(df_infit_engine)),
    DF_Infit_FACETS = suppressWarnings(as.numeric(fit_tbl$DF_Infit_FACETS)),
    DF_Outfit_ENGINE = suppressWarnings(as.numeric(df_outfit_engine)),
    DF_Outfit_FACETS = suppressWarnings(as.numeric(fit_tbl$DF_Outfit_FACETS)),
    InfitZSTD_ENGINE = suppressWarnings(as.numeric(z_infit_engine)),
    InfitZSTD_FACETS = suppressWarnings(as.numeric(fit_tbl$InfitZSTD_FACETS)),
    OutfitZSTD_ENGINE = suppressWarnings(as.numeric(z_outfit_engine)),
    OutfitZSTD_FACETS = suppressWarnings(as.numeric(fit_tbl$OutfitZSTD_FACETS)),
    stringsAsFactors = FALSE
  )
  out <- build_fit_measure_df_sensitivity(
    source_tbl,
    zstd_cut = 2,
    df_zstd_tolerance = df_zstd_tolerance,
    df_zstd_large_shift = df_zstd_large_shift,
    df_ratio_tolerance = df_ratio_tolerance
  )
  if (nrow(out) > 0L && "N" %in% names(source_tbl)) {
    source_key <- paste(source_tbl$Facet, source_tbl$Level, sep = "\r")
    out_key <- paste(out$Facet, out$Level, sep = "\r")
    out$N <- source_tbl$N[match(out_key, source_key)]
    out <- out[, c("Facet", "Level", "N", setdiff(names(out), c("Facet", "Level", "N"))), drop = FALSE]
  }
  out
}

external_fit_status <- function(max_mnsq_delta, max_zstd_delta, max_df_delta,
                                mnsq_tolerance, external_zstd_tolerance,
                                df_tolerance) {
  has_any <- is.finite(max_mnsq_delta) || is.finite(max_zstd_delta) || is.finite(max_df_delta)
  if (!has_any) return("insufficient_external_columns")
  mnsq_ok <- !is.finite(max_mnsq_delta) || max_mnsq_delta <= mnsq_tolerance
  zstd_ok <- !is.finite(max_zstd_delta) || max_zstd_delta <= external_zstd_tolerance
  df_ok <- !is.finite(max_df_delta) || max_df_delta <= df_tolerance
  if (mnsq_ok && zstd_ok && df_ok) return("same")
  if (mnsq_ok &&
      ((!zstd_ok && is.finite(max_zstd_delta) && max_zstd_delta <= 2 * external_zstd_tolerance) ||
       (!df_ok && is.finite(max_df_delta) && max_df_delta <= 2 * df_tolerance))) {
    return("rounding")
  }
  if (mnsq_ok && (!zstd_ok || !df_ok)) return("df_or_whexact_difference")
  if (!mnsq_ok) return("mnsq_or_measure_difference")
  "needs_review"
}

build_external_facets_fit_comparison <- function(internal_tbl,
                                                 external_tbl,
                                                 mnsq_tolerance = 0.01,
                                                 external_zstd_tolerance = 0.05,
                                                 df_tolerance = 0.5) {
  if (nrow(external_tbl) == 0L) {
    return(tibble::tibble())
  }
  external_tbl <- external_tbl |>
    dplyr::group_by(.data$Facet, .data$Level) |>
    dplyr::slice(1L) |>
    dplyr::ungroup()

  joined <- internal_tbl |>
    dplyr::transmute(
      Facet = .data$Facet,
      Level = .data$Level,
      MFRMR_Infit = .data$Infit,
      MFRMR_Outfit = .data$Outfit,
      MFRMR_InfitZSTD_FACETS = .data$InfitZSTD_FACETS,
      MFRMR_OutfitZSTD_FACETS = .data$OutfitZSTD_FACETS,
      MFRMR_DF_Infit_FACETS = .data$DF_Infit_FACETS,
      MFRMR_DF_Outfit_FACETS = .data$DF_Outfit_FACETS
    ) |>
    dplyr::left_join(external_tbl, by = c("Facet", "Level")) |>
    dplyr::mutate(
      ExternalMatched = !is.na(.data$Source),
      InfitDelta_FACETS_minus_mfrmr = .data$FACETS_Infit - .data$MFRMR_Infit,
      OutfitDelta_FACETS_minus_mfrmr = .data$FACETS_Outfit - .data$MFRMR_Outfit,
      InfitZSTDDelta_FACETS_minus_mfrmr = .data$FACETS_InfitZSTD - .data$MFRMR_InfitZSTD_FACETS,
      OutfitZSTDDelta_FACETS_minus_mfrmr = .data$FACETS_OutfitZSTD - .data$MFRMR_OutfitZSTD_FACETS,
      DFInfitDelta_FACETS_minus_mfrmr = .data$FACETS_DF_Infit - .data$MFRMR_DF_Infit_FACETS,
      DFOutfitDelta_FACETS_minus_mfrmr = .data$FACETS_DF_Outfit - .data$MFRMR_DF_Outfit_FACETS,
      MaxAbsMnSqDelta = fit_review_pmax_na(
        abs(.data$InfitDelta_FACETS_minus_mfrmr),
        abs(.data$OutfitDelta_FACETS_minus_mfrmr)
      ),
      MaxAbsZSTDDelta = fit_review_pmax_na(
        abs(.data$InfitZSTDDelta_FACETS_minus_mfrmr),
        abs(.data$OutfitZSTDDelta_FACETS_minus_mfrmr)
      ),
      MaxAbsDFDelta = fit_review_pmax_na(
        abs(.data$DFInfitDelta_FACETS_minus_mfrmr),
        abs(.data$DFOutfitDelta_FACETS_minus_mfrmr)
      )
    )

  joined$ExternalStatus <- ifelse(!joined$ExternalMatched, "no_external_match", vapply(
    seq_len(nrow(joined)),
    function(i) external_fit_status(
      joined$MaxAbsMnSqDelta[i],
      joined$MaxAbsZSTDDelta[i],
      joined$MaxAbsDFDelta[i],
      mnsq_tolerance = mnsq_tolerance,
      external_zstd_tolerance = external_zstd_tolerance,
      df_tolerance = df_tolerance
    ),
    character(1)
  ))

  joined |>
    dplyr::arrange(.data$ExternalStatus != "same",
                   dplyr::desc(.data$MaxAbsZSTDDelta),
                   dplyr::desc(.data$MaxAbsMnSqDelta),
                   .data$Facet, .data$Level)
}

facets_fit_review_guidance <- function(model, external_supplied) {
  tibble::tibble(
    Topic = c(
      "Primary comparison",
      "Residual basis",
      "ZSTD convention",
      "Small df",
      "External FACETS fit",
      "Bounded GPCM"
    ),
    Guidance = c(
      "Compare MnSq values separately from ZSTD values; MnSq differences indicate fit-statistic or estimation differences.",
      "MML fits evaluate residuals at shrunken EAP person measures while FACETS uses JMLE estimates, so MnSq itself can differ across the two bases; refit with method = \"JML\" before attributing MnSq gaps to fit computation.",
      "Use the FACETS-style companion columns for FACETS ZSTD comparison; engine ZSTD columns retain the package-native df convention.",
      "FACETS permits chi-square df below 1 under WHEXACT=Y; mfrmr withholds ZSTD as NA when df < 1, so an NA-vs-finite ZSTD pair on a sparse cell is an availability difference (compare MnSq for that row), not a fit difference.",
      if (isTRUE(external_supplied)) {
        "External rows are matched by Facet and Level. Rows without a match are marked no_external_match."
      } else {
        "No external FACETS table was supplied; the review reports internal engine-vs-FACETS-style standardization only."
      },
      if (identical(model, "GPCM")) {
        "Bounded GPCM has no direct FACETS free-slope counterpart; read this as an internal standardization review, not external FACETS equivalence."
      } else {
        "For RSM/PCM this review supports FACETS comparison, but it still does not prove full software equivalence."
      }
    )
  )
}

#' Review fit standardization against FACETS-style ZSTD conventions
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. If it does not
#'   contain FACETS-style fit columns, diagnostics are recomputed with
#'   `fit_df_method = "both"` and `residual_pca = "none"`.
#' @param facets_fit Optional external FACETS fit table, or a list of such
#'   tables. The helper matches rows by `Facet` and `Level`; a person-only table
#'   with a `Person` column is also accepted.
#' @param facet_col,level_col Optional explicit column names for the external
#'   FACETS table when automatic detection is not sufficient.
#' @param mnsq_tolerance,external_zstd_tolerance,df_tolerance Numeric
#'   tolerances used to classify external FACETS-vs-mfrmr differences.
#' @param df_zstd_tolerance Smallest absolute engine-vs-FACETS-style ZSTD
#'   difference treated as interpretively visible rather than rounding noise
#'   in `df_sensitivity`. Default `0.05`.
#' @param df_zstd_large_shift Absolute engine-vs-FACETS-style ZSTD difference
#'   labeled `large_zstd_shift` when the |ZSTD| flag status is unchanged.
#'   Default `0.5`.
#' @param df_ratio_tolerance Relative df-difference tolerance used to classify
#'   the internal engine-vs-FACETS-style df difference; for example, `0.05`
#'   means a 5 percent df difference.
#'
#' @details
#' This helper separates two questions that are often conflated when comparing
#' mfrmr output with FACETS:
#' - how much the package-native `engine` ZSTD changes when the same MnSq values
#'   are standardized with the FACETS/Wright-Masters fourth-moment df convention;
#' - when an external FACETS table is supplied, whether the FACETS-reported rows
#'   match mfrmr's FACETS-style companion columns closely enough for practical
#'   reporting.
#'
#' The review is row-matched by `Facet` and `Level`. It treats MnSq, ZSTD, and df
#' differences separately because FACETS documentation makes the df convention
#' and Wilson-Hilferty/WHEXACT handling central to ZSTD interpretation.
#'
#' Two upstream boundaries also apply. For `method = "MML"` fits, residuals
#' are evaluated at shrunken EAP person measures while FACETS uses JMLE
#' estimates, so MnSq itself can differ before standardization; refit with
#' `method = "JML"` for a JMLE-style residual basis. And mfrmr withholds
#' ZSTD as `NA` when the applicable df falls below 1 (Wilson-Hilferty
#' instability), while FACETS under `WHEXACT` can report a value on the same
#' sparse cell; such NA-vs-finite pairs are availability differences, not
#' fit differences. Both notes are repeated in the returned `guidance` table.
#'
#' @return An `mfrm_facets_fit_review` bundle with:
#' - `summary`: one-row overview of internal and external comparison counts
#' - `standardization`: the fit-standardization guide from diagnostics
#' - `df_sensitivity`: engine-vs-FACETS-style df/ZSTD comparison using
#'   the same row-level status taxonomy as `fit_measures_table()$df_sensitivity`
#' - `df_sensitive`: subset of `df_sensitivity` whose df convention changes
#'   the |ZSTD| flag or materially changes ZSTD interpretation
#' - `df_sensitivity_summary`: counts by df-sensitivity status
#' - `external_table_quality`: completeness and duplicate-key review for the
#'   supplied FACETS fit table
#' - `external_comparison`: optional external FACETS-vs-mfrmr comparison
#' - `df_conversion_guide`: formulas, column map, and comparison decisions for
#'   FACETS-style df/ZSTD review
#' - `guidance`: interpretation notes
#' - `settings`: tolerances and review metadata
#'
#' @seealso [diagnose_mfrm()], [facets_output_contract_review()],
#'   [mfrmr_compatibility_layer]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' review <- facets_fit_review(fit)
#' summary(review)
#' @export
facets_fit_review <- function(fit,
                             diagnostics = NULL,
                             facets_fit = NULL,
                             facet_col = NULL,
                             level_col = NULL,
                             mnsq_tolerance = 0.01,
                             external_zstd_tolerance = 0.05,
                             df_tolerance = 0.5,
                             df_zstd_tolerance = 0.05,
                             df_zstd_large_shift = 0.5,
                             df_ratio_tolerance = 0.05) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  mnsq_tolerance <- fit_measure_validate_nonnegative_finite(mnsq_tolerance, "mnsq_tolerance")
  external_zstd_tolerance <- fit_measure_validate_nonnegative_finite(
    external_zstd_tolerance,
    "external_zstd_tolerance"
  )
  df_tolerance <- fit_measure_validate_nonnegative_finite(df_tolerance, "df_tolerance")
  df_zstd_tolerance <- fit_measure_validate_nonnegative_finite(df_zstd_tolerance, "df_zstd_tolerance")
  df_zstd_large_shift <- fit_measure_validate_nonnegative_finite(df_zstd_large_shift, "df_zstd_large_shift")
  df_ratio_tolerance <- fit_measure_validate_nonnegative_finite(df_ratio_tolerance, "df_ratio_tolerance")
  if (df_zstd_large_shift < df_zstd_tolerance) {
    stop("`df_zstd_large_shift` must be greater than or equal to `df_zstd_tolerance`.", call. = FALSE)
  }
  model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  external_supplied <- !is.null(facets_fit)
  if (identical(model, "GPCM") && isTRUE(external_supplied)) {
    stop(
      "External FACETS fit comparison is not defined for bounded GPCM, ",
      "because FACETS does not estimate the package's free-slope GPCM route. ",
      "Run without `facets_fit` for an internal standardization review.",
      call. = FALSE
    )
  }

  diagnostics <- facets_fit_review_prepare_diagnostics(fit, diagnostics = diagnostics)
  fit_tbl <- tibble::as_tibble(diagnostics$fit %||% tibble::tibble())
  internal_tbl <- build_internal_fit_standardization_review(
    fit_tbl,
    df_zstd_tolerance = df_zstd_tolerance,
    df_zstd_large_shift = df_zstd_large_shift,
    df_ratio_tolerance = df_ratio_tolerance
  )
  df_sensitivity_summary <- summarize_fit_measure_df_sensitivity(internal_tbl)
  df_sensitive <- internal_tbl[
    !as.character(internal_tbl$DfSensitivityStatus %||% "not_available") %in%
      c("same_or_rounding", "not_available"),
    ,
    drop = FALSE
  ]
  external_tbl <- normalize_facets_fit_input(
    facets_fit,
    facet_col = facet_col,
    level_col = level_col
  )
  external_table_quality <- summarize_external_facets_fit_quality(external_tbl)
  external_comparison <- build_external_facets_fit_comparison(
    internal_tbl = internal_tbl,
    external_tbl = external_tbl,
    mnsq_tolerance = mnsq_tolerance,
    external_zstd_tolerance = external_zstd_tolerance,
    df_tolerance = df_tolerance
  )

  flag_changed <- sum(internal_tbl$FlagChangedByDf %in% TRUE, na.rm = TRUE)
  external_matched <- if (nrow(external_comparison) > 0L) {
    sum(external_comparison$ExternalMatched %in% TRUE, na.rm = TRUE)
  } else {
    0L
  }
  external_review <- if (nrow(external_comparison) > 0L) {
    sum(!external_comparison$ExternalStatus %in% c("same", "rounding") &
          external_comparison$ExternalMatched %in% TRUE, na.rm = TRUE)
  } else {
    0L
  }

  summary_tbl <- tibble::tibble(
    Model = model,
    Elements = nrow(internal_tbl),
    DfComparedRows = df_sensitivity_summary$ComparedRows[1],
    DfSensitiveRows = nrow(df_sensitive),
    DfSameOrRoundingRows = df_sensitivity_summary$SameOrRoundingRows[1],
    LargeZSTDShiftRows = df_sensitivity_summary$LargeZSTDShiftRows[1],
    DfConventionDifferenceRows = df_sensitivity_summary$DfConventionDifferenceRows[1],
    FlagChangedByDf = flag_changed,
    ExternalRows = nrow(external_tbl),
    ExternalDuplicateKeyRows = external_table_quality$DuplicateFacetLevelRows[1],
    ExternalCompleteMnSqRows = external_table_quality$CompleteMnSqRows[1],
    ExternalCompleteZSTDRows = external_table_quality$CompleteZSTDRows[1],
    ExternalCompleteDFRows = external_table_quality$CompleteDFRows[1],
    ExternalMatched = external_matched,
    ExternalNeedsReview = external_review,
    ExternalComparison = if (external_supplied) "supplied" else "not_supplied"
  )

  out <- list(
    summary = summary_tbl,
    standardization = tibble::as_tibble(diagnostics$fit_standardization %||% tibble::tibble()),
    df_sensitivity = internal_tbl,
    df_sensitive = df_sensitive,
    df_sensitivity_summary = df_sensitivity_summary,
    external_table_quality = external_table_quality,
    external_comparison = external_comparison,
    df_conversion_guide = facets_fit_df_guide(include_references = TRUE),
    guidance = facets_fit_review_guidance(model, external_supplied = external_supplied),
    settings = list(
      intended_use = "fit_standardization_review",
      external_validation = isTRUE(external_supplied) && !identical(model, "GPCM"),
      fit_df_method = "both",
      mnsq_tolerance = mnsq_tolerance,
      external_zstd_tolerance = external_zstd_tolerance,
      df_tolerance = df_tolerance,
      df_zstd_tolerance = df_zstd_tolerance,
      df_zstd_large_shift = df_zstd_large_shift,
      df_ratio_tolerance = df_ratio_tolerance
    )
  )
  as_mfrm_bundle(out, "mfrm_facets_fit_review")
}

#' Build a FACETS output-contract review
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
#' This function checks produced report components against a FACETS-style
#' output-contract specification (`inst/references/facets_column_contract.csv`) and
#' returns:
#' - column-level coverage per contract row
#' - table-level coverage summaries
#' - optional metric-level consistency checks
#'
#' It is intended for output-contract QA and regression review. It does
#' not establish external validity or software equivalence beyond the specific
#' schema/metric contract encoded in the contract file.
#'
#' @section Bounded GPCM boundary:
#' This helper remains blocked for bounded `GPCM` fits in 0.2.1. The FACETS
#' output contract includes score-side rows whose measure-to-score and
#' uncertainty semantics are validated for the current Rasch-family route, not
#' for free-discrimination bounded `GPCM`. Use [gpcm_capability_matrix()] before
#' routing a bounded `GPCM` fit into score-side compatibility-output helpers.
#'
#' Coverage interpretation in `overall`:
#' - `MeanColumnCoverage` and `MinColumnCoverage` are computed across all
#'   contract rows (unavailable rows count as 0 coverage).
#' - `MeanColumnCoverageAvailable` and `MinColumnCoverageAvailable` summarize
#'   only rows whose source component is available.
#'
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_facets_contract_review` (`type = "column_coverage"`, `"table_coverage"`,
#' `"metric_status"`, `"metric_by_table"`).
#'
#' @section Interpreting output:
#' - `overall`: high-level output-contract coverage and metric-check pass
#'   rates.
#' - `column_summary` / `column_review`: where output-schema mismatches
#'   occur.
#' - `metric_summary` / `metric_checks`: numerical consistency checks tied to the
#'   current contract.
#' - `missing_preview`: direct path to unresolved output-contract gaps.
#'
#' @section Typical workflow:
#' 1. Run `facets_output_contract_review(fit, branch = "facets")`.
#' 2. Inspect `summary(contract_review)` and `missing_preview`.
#' 3. Patch upstream table builders, then rerun the output-contract review.
#'
#' @return
#' An object of class `mfrm_facets_contract_review` with:
#' - `overall`: one-row output-contract review summary
#' - `column_summary`: coverage summary by table ID
#' - `column_review`: row-level output-contract review
#' - `missing_preview`: lowest-coverage rows
#' - `metric_summary`: one-row metric-check summary
#' - `metric_by_table`: metric-check summary by table ID
#' - `metric_checks`: row-level metric checks
#' - `settings`: branch/contract metadata
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [build_fixed_reports()],
#'   [mfrmr_compatibility_layer]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' contract_review <- facets_output_contract_review(fit, diagnostics = diag, branch = "facets")
#' summary(contract_review)
#' p <- plot(contract_review, draw = FALSE)
#' }
#' @export
facets_output_contract_review <- function(fit,
                                 diagnostics = NULL,
                                 bias_results = NULL,
                                 branch = c("facets", "original"),
                                 contract_file = NULL,
                                 include_metrics = TRUE,
                                 top_n_missing = 15L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  stop_if_gpcm_out_of_scope(fit, "facets_output_contract_review()")
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

  column_review_rows <- lapply(seq_len(nrow(contract)), function(i) {
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
  column_review <- dplyr::bind_rows(column_review_rows)

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
    column_review$available %in% TRUE,
    suppressWarnings(as.numeric(column_review$coverage)),
    0
  )
  contract_coverage_values[!is.finite(contract_coverage_values)] <- 0

  column_summary <- column_review |>
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

  missing_preview <- column_review |>
    dplyr::filter(!.data$full_match | !.data$available) |>
    dplyr::arrange(.data$coverage, .data$table_id, .data$component) |>
    dplyr::slice_head(n = top_n_missing)

  metric_checks <- if (isTRUE(include_metrics)) {
    build_contract_metric_review(outputs = outputs)
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

  metric_summary <- if (nrow(metric_checks) == 0) {
    data.frame(
      Checks = 0L,
      Evaluated = 0L,
      Passed = 0L,
      Failed = 0L,
      PassRate = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    ev <- metric_checks$Pass[!is.na(metric_checks$Pass)]
    data.frame(
      Checks = nrow(metric_checks),
      Evaluated = length(ev),
      Passed = sum(ev %in% TRUE),
      Failed = sum(ev %in% FALSE),
      PassRate = if (length(ev) > 0) sum(ev %in% TRUE) / length(ev) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  metric_by_table <- if (nrow(metric_checks) == 0) {
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
    metric_checks |>
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
  mean_cov_available <- summarize_coverage(column_review$coverage, mean)
  min_cov_available <- summarize_coverage(column_review$coverage, min)
  contract_rows <- nrow(column_review)
  mismatches <- sum(!column_review$full_match, na.rm = TRUE)
  overall <- data.frame(
    Branch = branch,
    ContractRows = contract_rows,
    AvailableRows = sum(column_review$available, na.rm = TRUE),
    FullMatchRows = sum(column_review$full_match, na.rm = TRUE),
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
    column_review = as.data.frame(column_review, stringsAsFactors = FALSE),
    missing_preview = as.data.frame(missing_preview, stringsAsFactors = FALSE),
    metric_summary = metric_summary,
    metric_by_table = as.data.frame(metric_by_table, stringsAsFactors = FALSE),
    metric_checks = as.data.frame(metric_checks, stringsAsFactors = FALSE),
    settings = list(
      branch = branch,
      contract_path = contract_info$path,
      intended_use = "facets_output_contract_review",
      external_validation = FALSE,
      include_metrics = include_metrics,
      top_n_missing = top_n_missing,
      bias_included = !is.null(outputs$t10)
    )
  )
  as_mfrm_bundle(out, "mfrm_facets_contract_review")
}

#' Build a package-native reference review for report completeness
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. If omitted,
#'   diagnostics are computed internally with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()]. If omitted and
#'   at least two facets exist, a 2-way interaction screen is computed internally.
#' @param reference_profile Review profile. `"core"` emphasizes package-native
#'   report contracts. `"compatibility"` exposes the manual-aligned compatibility
#'   layer used by `facets_output_contract_review(branch = "facets")`.
#' @param include_metrics If `TRUE`, run numerical consistency checks in addition
#'   to schema coverage checks.
#' @param top_n_attention Number of lowest-coverage components to keep in
#'   `attention_items`.
#'
#' @details
#' This function repackages the output-contract review into package-native
#' terminology so users can review output completeness without needing external
#' manual/table numbering. It reports:
#' - component-level schema coverage
#' - numerical consistency checks for derived report tables
#' - the highest-priority attention items for follow-up
#'
#' It is a package-output completeness review, not an external validation
#' study.
#'
#' Use `reference_profile = "core"` for ordinary `mfrmr` workflows.
#' Use `reference_profile = "compatibility"` only when you explicitly want to
#' inspect the compatibility layer.
#'
#' @section Interpreting output:
#' - `overall`: one-row review summary with schema coverage and metric
#'   pass rate.
#' - `component_summary`: per-component coverage summary.
#' - `attention_items`: direct list of components needing review.
#' - `metric_summary` / `metric_checks`: numerical consistency status.
#'
#' @return An object of class `mfrm_reference_review`.
#' @seealso [facets_output_contract_review()], [diagnose_mfrm()], [build_fixed_reports()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' review <- reference_case_review(fit, diagnostics = diag)
#' summary(review)
#' }
#' @name reference_case_review
#' @export
reference_case_review <- function(fit,
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

  contract_review <- facets_output_contract_review(
    fit = fit,
    diagnostics = diagnostics,
    bias_results = bias_results,
    branch = branch,
    include_metrics = include_metrics,
    top_n_missing = top_n_attention
  )

  overall_src <- as.data.frame(contract_review$overall, stringsAsFactors = FALSE)
  overall <- tibble::tibble(
    ReferenceProfile = reference_profile,
    ContractBranch = as.character(overall_src$Branch[1] %||% branch),
    SchemaCoverage = as.numeric(overall_src$MeanColumnCoverage[1] %||% NA_real_),
    AvailableSchemaCoverage = as.numeric(overall_src$MeanColumnCoverageAvailable[1] %||% NA_real_),
    MinSchemaCoverage = as.numeric(overall_src$MinColumnCoverage[1] %||% NA_real_),
    MetricPassRate = as.numeric(overall_src$MetricPassRate[1] %||% NA_real_),
    SchemaMismatches = as.integer(overall_src$ColumnMismatches[1] %||% NA_integer_),
    AttentionItems = nrow(contract_review$missing_preview %||% data.frame()),
    CompatibilityLayer = if (identical(reference_profile, "compatibility")) "manual-aligned" else "package-native"
  )

  component_summary <- as.data.frame(contract_review$column_summary, stringsAsFactors = FALSE)
  names(component_summary) <- sub("^table_id$", "ComponentID", names(component_summary))
  names(component_summary) <- sub("^function_name$", "Builder", names(component_summary))

  attention_items <- as.data.frame(contract_review$missing_preview, stringsAsFactors = FALSE)
  names(attention_items) <- sub("^table_id$", "ComponentID", names(attention_items))
  names(attention_items) <- sub("^function_name$", "Builder", names(attention_items))
  names(attention_items) <- sub("^component$", "Subtable", names(attention_items))
  names(attention_items) <- sub("^coverage$", "Coverage", names(attention_items))
  names(attention_items) <- sub("^missing$", "MissingColumns", names(attention_items))

  out <- list(
    overall = overall,
    component_summary = component_summary,
    attention_items = attention_items,
    metric_summary = as.data.frame(contract_review$metric_summary, stringsAsFactors = FALSE),
    metric_checks = as.data.frame(contract_review$metric_checks, stringsAsFactors = FALSE),
    settings = list(
      reference_profile = reference_profile,
      contract_branch = branch,
      intended_use = "reference_contract_review",
      external_validation = FALSE,
      include_metrics = isTRUE(include_metrics),
      top_n_attention = max(1L, as.integer(top_n_attention))
    ),
    contract_review = contract_review
  )
  as_mfrm_bundle(out, "mfrm_reference_review")
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
#' - `$gpcm_boundary`: for bounded `GPCM` inputs, a capability-boundary table
#'   marking the narrative as caveated DFF screening output.
#' - `$config`: analysis configuration inherited from the input.
#'
#' @section GPCM boundary:
#' If the input comes from a bounded `GPCM` fit, the narrative includes a
#' bounded-`GPCM` note and the returned report carries `gpcm_boundary`.
#' Treat the text as slope-aware screening/reporting support, not as a
#' standalone fairness, invariance, or operational subgroup decision.
#'
#' @section Typical workflow:
#' 1. Run [analyze_dff()] / [analyze_dif()] or [dif_interaction_table()].
#' 2. Pass the result to `dif_report()`.
#' 3. Print the report or extract `$narrative` for inclusion in a
#'    manuscript.
#'
#' @return Object of class `mfrm_dif_report` with `narrative`,
#'   `counts`, `large_dif`, `gpcm_boundary`, and `config`.
#'
#' @section References:
#' The narrative caveat about distinguishing construct-relevant variation
#' from unwanted measurement bias is grounded in:
#'
#' - Eckes, T. (2011). *Introduction to Many-Facet Rasch Measurement:
#'   Analyzing and Evaluating Rater-Mediated Assessments*. Frankfurt am
#'   Main: Peter Lang. ISBN 978-3-631-61350-4.
#' - McNamara, T., & Knoch, U. (2012). The Rasch wars: The emergence of
#'   Rasch measurement in language testing. *Language Testing*, 29(4),
#'   555--576. \doi{10.1177/0265532211430367}
#'
#' @seealso [analyze_dff()], [analyze_dif()], [dif_interaction_table()],
#'   [plot_dif_heatmap()], [build_apa_outputs()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 30)
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
  gpcm_boundary <- dif_result$gpcm_boundary %||% data.frame()
  if (is.data.frame(gpcm_boundary) && nrow(gpcm_boundary) > 0L) {
    lines <- c(lines, paste0(
      "\nBounded GPCM note: Treat these differential-functioning rows as ",
      "slope-aware screening evidence under the current bounded-GPCM fit. ",
      "They do not by themselves establish fairness, invariance, or an ",
      "operational subgroup decision."
    ))
  }

  narrative <- paste(lines, collapse = "")

  out <- list(
    narrative = narrative,
    counts = counts,
    large_dif = tibble::as_tibble(large_dif),
    gpcm_boundary = gpcm_boundary,
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
  gpcm_boundary <- dif_result$gpcm_boundary %||% data.frame()
  if (is.data.frame(gpcm_boundary) && nrow(gpcm_boundary) > 0L) {
    lines <- c(lines, paste0(
      "\nBounded GPCM note: Treat these interaction-screening rows as ",
      "slope-aware residual evidence under the current bounded-GPCM fit. ",
      "They do not by themselves establish fairness, invariance, or an ",
      "operational subgroup decision."
    ))
  }

  narrative <- paste(lines, collapse = "")

  out <- list(
    narrative = narrative,
    counts = counts,
    large_dif = tibble::as_tibble(flagged_rows),
    gpcm_boundary = gpcm_boundary,
    config = cfg
  )
  class(out) <- c("mfrm_dif_report", class(out))
  out
}

#' @export
print.mfrm_dif_report <- function(x, ...) {
  cat("--- Differential Functioning Interpretation Report ---\n\n")
  cat(x$narrative, "\n")
  .print_dff_gpcm_boundary(x$gpcm_boundary)
  invisible(x)
}

#' @export
summary.mfrm_dif_report <- function(object, ...) {
  out <- list(
    narrative = object$narrative,
    counts = object$counts,
    large_dif = object$large_dif,
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
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
  .print_dff_gpcm_boundary(x$gpcm_boundary)
  invisible(x)
}

# ---- QC Pipeline ---------------------------------------------------------

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
#' For bounded `GPCM`, this pipeline is available as caveated operational
#' triage over supported diagnostics. Its pass/warn/fail labels remain package
#' QC policy overlays; they are not FACETS score-side equivalence, operational
#' scoring decisions, design-forecasting evidence, or automatic fairness /
#' validity decisions.
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
#'                 method = "JML", maxit = 30)
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "run_qc_pipeline()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review",
        "Design planning and forecasting"
      )
    ),
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
  if (!is.null(x$gpcm_boundary) && nrow(as.data.frame(x$gpcm_boundary)) > 0) {
    cat("\nGPCM Boundary:\n")
    print(as.data.frame(x$gpcm_boundary)[, c("Area", "Status"), drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.mfrm_qc_pipeline <- function(object, ...) {
  out <- list(
    verdicts = object$verdicts,
    overall  = object$overall,
    recommendations = object$recommendations,
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
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
