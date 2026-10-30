#' Facet-quality dashboard for facet-level screening
#'
#' Build a compact dashboard for one facet at a time, combining facet
#' severity, misfit, central-tendency screening, and optional bias counts.
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] for the same `fit`.
#'   Mismatched or outdated diagnostics are rejected. Recompute diagnostics
#'   after refitting; the model itself does not need to be fitted again.
#' @param facet Optional facet name. When `NULL`, the function tries to infer
#'   a rater-like facet and otherwise falls back to the first modeled facet.
#' @param bias_results Optional output from [estimate_bias()] or a named list
#'   of such outputs. Non-matching bundles are skipped quietly.
#' @param severity_warn Absolute estimate cutoff used to flag severity
#'   outliers.
#' @param misfit_warn `NULL` (default) uses both bounds returned by
#'   [mfrm_misfit_thresholds()], including any session options. A numeric value
#'   instead uses that upper cutoff and its reciprocal as the lower cutoff;
#'   it does not preserve the current lower bound. For example, `1.5` selects
#'   about 0.67--1.5, whereas the unmodified package defaults are 0.5--1.5.
#' @param central_tendency_max Legacy opt-in absolute estimate cutoff for
#'   marking facet estimates near the fitted origin. The default `NULL`
#'   disables this flag because origin proximity is not evidence that a rater
#'   avoids extreme score categories. Use [data_quality_report()] for observed
#'   category-use screening and [rater_centrality_table()] for model-based
#'   indirect centrality indices.
#' @param bias_count_warn Minimum flagged-bias row count required to flag a
#'   level.
#' @param bias_abs_t_warn Absolute `t` cutoff used when deriving bias-row
#'   flags from a raw bias bundle.
#' @param bias_abs_size_warn Absolute bias-size cutoff used when deriving
#'   bias-row flags from a raw bias bundle.
#' @param bias_p_max Probability cutoff used when deriving bias-row flags
#'   from a raw bias bundle.
#'
#' @details
#' The dashboard screens individual facet elements across four
#' complementary criteria:
#'
#' - **Severity**: elements with \eqn{|\mathrm{Estimate}| >}
#'   `severity_warn` logits are flagged as unusually harsh or lenient.
#' - **Misfit**: elements with Infit or Outfit MnSq outside the
#'   configured heuristic review band are flagged. The band defaults to the package
#'   pair returned by [mfrm_misfit_thresholds()] (Linacre 0.5-1.5);
#'   pass `misfit_warn = 1.5` to keep the older symmetric
#'   \eqn{[1/}\code{misfit_warn}\eqn{,\;}\code{misfit_warn}\eqn{]}
#'   form (0.67-1.5).
#' - **Reference proximity (legacy `CentralTendencyFlag` label)**: when
#'   `central_tendency_max` is supplied, elements with
#'   \eqn{|\mathrm{Estimate}| <} `central_tendency_max` logits are marked.
#'   This only describes proximity to the fitted origin and must not be
#'   interpreted as central-category use or restriction of range. Those are
#'   response-pattern questions; see [data_quality_report()] and
#'   [rater_centrality_table()] for their distinct descriptive indices.
#' - **Bias**: elements involved in \eqn{\ge} `bias_count_warn`
#'   screen-positive interaction cells (from [estimate_bias()]) are flagged.
#'
#' A **flag density** score counts how many enabled criteria each element
#' triggers. Elements flagged on multiple criteria warrant priority
#' review and may motivate training or a documented data-quality review; the
#' dashboard does not justify automatic row, person, or rater exclusion.
#'
#' Default thresholds are screening heuristics. Prespecify and justify any
#' application-specific alternatives rather than treating them as universal
#' validity or acceptance criteria.
#' `MissingMetrics` lists unavailable estimate, SE, Infit or Outfit values.
#' A severity or misfit flag is `NA` when it cannot be evaluated (an observed
#' misfit exceedance still flags even if the other fit index is unavailable).
#' `FlagCount` counts observed flags only; zero is not a complete pass.
#' `IncompleteLevels` counts levels with any missing diagnostic. Bias counts
#' describe supplied results only, not tests of absence of bias.
#' Fit-readiness restrictions are retained in `fit_readiness`,
#' `interpretation_status` and `notes`. Review overlap with
#' [subset_connectivity_report()] and category use with [data_quality_report()]
#' before interpreting between-rater differences.
#'
#' @return An object of class `mfrm_facet_dashboard` (also inheriting from
#'   `mfrm_bundle` and `list`). The object summarizes one target facet:
#'   `overview` reports the facet-level screening totals, `summary` provides
#'   aggregate estimates and flag counts, `detail` contains one row per facet
#'   level with the computed screening indicators, `ranked` orders levels by
#'   review priority, `flagged` keeps only levels requiring follow-up,
#'   `bias_sources` records which bias-result bundles contributed to the
#'   counts, `settings` stores the resolved thresholds, and `notes` gives short
#'   interpretation messages about how to read the dashboard.
#'
#' @section Output:
#' The returned object is a bundle-like list with class
#' `mfrm_facet_dashboard` and components:
#' - `facet`: character scalar naming the dashboard's target facet
#' - `facet_source`: character scalar describing whether the target
#'   facet was inferred from the fit configuration or supplied
#'   explicitly
#' - `overview`: one-row structural overview
#' - `summary`: one-row screening summary
#' - `detail`: level-level detail table
#' - `ranked`: detail ordered by flag density / severity
#' - `flagged`: flagged levels only
#' - `bias_sources`: per-bundle bias aggregation metadata
#' - `settings`: resolved threshold settings
#' - `notes`: short interpretation notes
#' - `diagnostics`: the `mfrm_diagnostics` bundle the dashboard was
#'   built from (echoed for downstream helpers that need to traverse
#'   the same diagnostics object)
#' - `bias_results`: the `mfrm_bias` bundle (or list of bundles)
#'   when `bias_results` was supplied; `NULL` otherwise
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()], [plot(x, type = "qc")][plot.mfrm_fit]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("facet_quality_dashboard", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' toy <- toy[toy$Person %in% unique(toy$Person)[1:8], ]
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' dash <- facet_quality_dashboard(fit, diagnostics = diag)
#' summary(dash)
#' }
#' @section Saved-result reuse:
#' Supplied screening results must describe this fitted analysis. Results from
#' another fit or older results without identity metadata are refused. See
#' [mfrmr_compatibility_layer] to recreate them from the saved fit and matching
#' diagnostics without refitting the model.
#' @export
#' @order 1
facet_quality_dashboard <- function(fit,
                                    diagnostics = NULL,
                                    facet = NULL,
                                    bias_results = NULL,
                                    severity_warn = 1.0,
                                    misfit_warn = NULL,
                                    central_tendency_max = NULL,
                                    bias_count_warn = 1L,
                                    bias_abs_t_warn = 2,
                                    bias_abs_size_warn = 0.5,
                                    bias_p_max = 0.05) {
  # When `misfit_warn` is NULL, defer to the package-level threshold
  # pair so all helpers (summary.mfrm_diagnostics, build_apa_outputs,
  # build_misfit_casebook, this dashboard) can be steered together via
  # `mfrm_misfit_thresholds()`.
  thresholds <- mfrm_misfit_thresholds()
  if (is.null(misfit_warn)) {
    misfit_warn <- as.numeric(thresholds["upper"])
    misfit_lower_band <- as.numeric(thresholds["lower"])
  } else {
    misfit_lower_band <- 1 / abs(as.numeric(misfit_warn[1]))
  }
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }

  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (!is.list(diagnostics) || is.null(diagnostics$measures)) {
    stop("`diagnostics` must be output from diagnose_mfrm().", call. = FALSE)
  }
  mfrm_results_validate_diagnostics_identity(fit, diagnostics,
                                             helper = "facet_quality_dashboard()")
  mfrm_validate_bias_identity(fit, bias_results, helper = "facet_quality_dashboard()")

  facet_names <- fit$config$facet_names
  if (is.null(facet_names)) facet_names <- fit$prep$facet_names
  facet_names <- as.character(facet_names)
  facet_names <- facet_names[!is.na(facet_names) & nzchar(facet_names)]

  if (length(facet_names) == 0) {
    stop("`fit` does not expose any facet names.", call. = FALSE)
  }

  facet_input <- facet
  if (is.null(facet_input) || !nzchar(as.character(facet_input[1]))) {
    facet <- infer_default_rater_facet(facet_names)
  } else {
    facet <- as.character(facet_input[1])
  }
  if (!facet %in% facet_names) {
    stop("`facet` must be one of the modeled facets: ", paste(facet_names, collapse = ", "), ".",
         call. = FALSE)
  }

  measures <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  needed <- c("Facet", "Level", "Estimate")
  if (!all(needed %in% names(measures))) {
    stop("`diagnostics$measures` does not include the required facet columns.", call. = FALSE)
  }

  detail <- measures[as.character(measures$Facet) == facet, , drop = FALSE]
  if (nrow(detail) == 0) {
    stop("The selected facet was not found in `diagnostics$measures`.", call. = FALSE)
  }

  num_or_na <- function(x) suppressWarnings(as.numeric(x))
  first_existing <- function(tbl, candidates) {
    for (nm in candidates) {
      if (nm %in% names(tbl)) return(num_or_na(tbl[[nm]]))
    }
    rep(NA_real_, nrow(tbl))
  }

  detail$Level <- as.character(detail$Level)
  detail$Estimate <- num_or_na(detail$Estimate)
  detail$SE <- first_existing(detail, c("SE", "ModelSE", "RealSE"))
  detail$Infit <- first_existing(detail, c("Infit", "InfitMnSq", "Infit Mnsq"))
  detail$Outfit <- first_existing(detail, c("Outfit", "OutfitMnSq", "Outfit Mnsq"))
  detail$AbsEstimate <- abs(detail$Estimate)
  detail$MissingMetrics <- vapply(seq_len(nrow(detail)), function(i) {
    metrics <- c("Estimate", "SE", "Infit", "Outfit")
    paste(metrics[!is.finite(as.numeric(detail[i, metrics]))], collapse = ", ")
  }, character(1))
  detail$SeverityFlag <- ifelse(is.finite(detail$AbsEstimate),
                                detail$AbsEstimate >= abs(severity_warn), NA)
  fit_indices <- as.matrix(detail[c("Infit", "Outfit")])
  fit_indices[!is.finite(fit_indices)] <- NA_real_
  detail$MisfitFlag <- rowSums(
    fit_indices >= abs(misfit_warn) | fit_indices <= misfit_lower_band,
    na.rm = TRUE
  ) > 0L
  detail$MisfitFlag[!detail$MisfitFlag &
    (!is.finite(detail$Infit) | !is.finite(detail$Outfit))] <- NA
  central_tendency_value <- suppressWarnings(as.numeric(central_tendency_max[1]))
  central_tendency_enabled <- length(central_tendency_value) == 1L &&
    is.finite(central_tendency_value)
  central_tendency_value <- if (isTRUE(central_tendency_enabled)) {
    abs(central_tendency_value)
  } else {
    NA_real_
  }
  detail$CentralTendencyFlag <- isTRUE(central_tendency_enabled) &
    is.finite(detail$AbsEstimate) & detail$AbsEstimate <= central_tendency_value
  detail$BiasCount <- 0L
  detail$BiasSources <- 0L

  bias_meta <- dashboard_bias_level_counts(
    bias_results = bias_results,
    target_facet = facet,
    bias_abs_t_warn = abs(bias_abs_t_warn),
    bias_abs_size_warn = abs(bias_abs_size_warn),
    bias_p_max = max(0, min(1, as.numeric(bias_p_max[1])))
  )
  if (nrow(bias_meta$levels) > 0) {
    idx <- match(detail$Level, bias_meta$levels$Level)
    hit <- !is.na(idx)
    detail$BiasCount[hit] <- as.integer(bias_meta$levels$BiasCount[idx[hit]])
    detail$BiasSources[hit] <- as.integer(bias_meta$levels$BiasSources[idx[hit]])
  }
  detail$BiasFlag <- is.finite(detail$BiasCount) & detail$BiasCount >= as.integer(bias_count_warn)
  detail$FlagCount <- rowSums(cbind(
    detail$SeverityFlag,
    detail$MisfitFlag,
    detail$CentralTendencyFlag,
    detail$BiasFlag
  ), na.rm = TRUE)
  detail$AnyFlag <- detail$FlagCount > 0
  detail$FlagLabel <- vapply(seq_len(nrow(detail)), function(i) {
    labs <- character(0)
    if (isTRUE(detail$SeverityFlag[i])) labs <- c(labs, "severity")
    if (isTRUE(detail$MisfitFlag[i])) labs <- c(labs, "misfit")
    if (isTRUE(detail$CentralTendencyFlag[i])) labs <- c(labs, "reference proximity")
    if (isTRUE(detail$BiasFlag[i])) labs <- c(labs, "bias")
    if (length(labs) == 0) "" else paste(labs, collapse = ", ")
  }, character(1))

  ranked <- detail |>
    dplyr::mutate(
      .AnyFlagRank = as.integer(.data$AnyFlag),
      .FlagScore = .data$FlagCount,
      .AbsEstimate = abs(.data$Estimate)
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$.AnyFlagRank),
      dplyr::desc(.data$.FlagScore),
      dplyr::desc(.data$.AbsEstimate),
      .data$Level
    ) |>
    dplyr::select(-".AnyFlagRank", -".FlagScore")

  flagged <- ranked[ranked$AnyFlag %in% TRUE, , drop = FALSE]

  overview <- data.frame(
    Facet = facet,
    FacetSource = if (is.null(facet_input) || !nzchar(as.character(facet_input[1]))) "inferred" else "user",
    Levels = nrow(detail),
    FlaggedLevels = sum(detail$AnyFlag, na.rm = TRUE),
    IncompleteLevels = sum(nzchar(detail$MissingMetrics)),
    BiasSourceBundles = nrow(bias_meta$sources[bias_meta$sources$Used %in% TRUE, , drop = FALSE]),
    stringsAsFactors = FALSE
  )

  summary_tbl <- data.frame(
    Facet = facet,
    Levels = nrow(detail),
    MeanEstimate = mean(detail$Estimate, na.rm = TRUE),
    SD = stats::sd(detail$Estimate, na.rm = TRUE),
    MinEstimate = min(detail$Estimate, na.rm = TRUE),
    MaxEstimate = max(detail$Estimate, na.rm = TRUE),
    MeanInfit = mean(detail$Infit, na.rm = TRUE),
    MeanOutfit = mean(detail$Outfit, na.rm = TRUE),
    SeverityFlagged = sum(detail$SeverityFlag, na.rm = TRUE),
    MisfitFlagged = sum(detail$MisfitFlag, na.rm = TRUE),
    CentralTendencyFlagged = sum(detail$CentralTendencyFlag, na.rm = TRUE),
    BiasFlagged = sum(detail$BiasFlag, na.rm = TRUE),
    AnyFlagged = sum(detail$AnyFlag, na.rm = TRUE),
    BiasRows = sum(detail$BiasCount, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  settings <- bundle_settings_table(list(
    facet = facet,
    facet_source = overview$FacetSource[1],
    severity_warn = abs(severity_warn),
    misfit_warn = abs(misfit_warn),
    misfit_lower = misfit_lower_band,
    central_tendency_max = central_tendency_value,
    bias_count_warn = as.integer(bias_count_warn),
    bias_abs_t_warn = abs(bias_abs_t_warn),
    bias_abs_size_warn = abs(bias_abs_size_warn),
    bias_p_max = max(0, min(1, as.numeric(bias_p_max[1]))),
    bias_source_bundles = nrow(bias_meta$sources[bias_meta$sources$Used %in% TRUE, , drop = FALSE])
  ))

  readiness <- .mfrm_fit_plot_readiness(fit)
  notes <- c(
    "Flags are screening prompts, not evidence of invalid ratings or grounds for automatic exclusion.",
    "Severity is relative to the fitted reference; inspect workload, category use and common ratings before comparing levels.",
    "FlagCount counts observed flags only. MissingMetrics identifies unavailable diagnostics; zero flags does not mean all checks passed.",
    "BiasCount counts flagged cells in supplied bias results only; zero does not establish absence of bias.",
    readiness$detail
  )
  notes <- c(
    notes,
    if (isTRUE(central_tendency_enabled)) {
      paste(
        "Legacy CentralTendencyFlag is an origin-proximity marker only;",
        "use data_quality_report() for observed category avoidance or range restriction."
      )
    } else {
      paste(
        "Legacy CentralTendencyFlag is disabled by default because origin proximity",
        "does not diagnose observed category avoidance or range restriction."
      )
    }
  )
  if (sum(detail$AnyFlag, na.rm = TRUE) == 0L) {
    notes <- c(notes, "No level-level flags were triggered under the current thresholds.")
  }
  if (nrow(bias_meta$sources) > 0 && any(grepl("^pair error:", bias_meta$sources$Reason))) {
    notes <- c(notes, "Some requested bias bundles failed and were excluded from the dashboard counts.")
  }
  if (nrow(bias_meta$sources) > 0 &&
      any(bias_meta$sources$Used %in% FALSE & !grepl("^pair error:", bias_meta$sources$Reason))) {
    notes <- c(notes, "Some bias bundles were skipped because they did not involve the target facet.")
  }
  if (length(notes) == 0) {
    notes <- "Dashboard constructed successfully."
  }

  out <- list(
    facet = facet,
    facet_source = overview$FacetSource[1],
    overview = overview,
    summary = summary_tbl,
    detail = detail,
    ranked = ranked,
    flagged = flagged,
    bias_sources = bias_meta$sources,
    settings = settings,
    fit_readiness = readiness$table,
    interpretation_status = readiness$status,
    notes = notes,
    diagnostics = diagnostics,
    bias_results = bias_results
  )
  class(out) <- c("mfrm_facet_dashboard", "mfrm_bundle", "list")
  out
}


dashboard_round_numeric_df <- function(df, digits = 3L) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  out <- df
  is_num <- vapply(out, is.numeric, logical(1))
  out[is_num] <- lapply(out[is_num], function(x) round(x, digits = digits))
  out
}

dashboard_normalize_bias_inputs <- function(bias_results) {
  if (is.null(bias_results)) return(list())
  error_tbl <- data.frame()
  if (inherits(bias_results, "mfrm_bias_collection")) {
    error_tbl <- as.data.frame(bias_results$errors %||% data.frame(), stringsAsFactors = FALSE)
    bias_results <- bias_results$by_pair %||% list()
  }
  if (is.list(bias_results) && !is.null(bias_results$table) && nrow(bias_results$table) > 0) {
    out <- stats::setNames(list(bias_results), "bias_1")
    attr(out, "errors") <- error_tbl
    return(out)
  }
  if (!is.list(bias_results) || length(bias_results) == 0) {
    out <- list()
    attr(out, "errors") <- error_tbl
    return(out)
  }

  out <- list()
  nm <- names(bias_results)
  if (is.null(nm)) nm <- rep("", length(bias_results))
  for (i in seq_along(bias_results)) {
    item <- bias_results[[i]]
    if (is.list(item) && !is.null(item$table) && nrow(item$table) > 0) {
      key <- nm[i]
      if (!nzchar(key)) key <- paste0("bias_", i)
      out[[key]] <- item
    }
  }
  attr(out, "errors") <- error_tbl
  out
}

dashboard_bias_level_counts <- function(bias_results,
                                        target_facet,
                                        bias_abs_t_warn = 2,
                                        bias_abs_size_warn = 0.5,
                                        bias_p_max = 0.05) {
  sources <- dashboard_normalize_bias_inputs(bias_results)
  error_tbl <- attr(sources, "errors", exact = TRUE)
  if (length(sources) == 0) {
    source_tbl <- if (is.data.frame(error_tbl) && nrow(error_tbl) > 0) {
      data.frame(
        Source = as.character(error_tbl$Interaction %||% paste0("bias_error_", seq_len(nrow(error_tbl)))),
        Used = FALSE,
        Reason = paste0("pair error: ", as.character(error_tbl$Error)),
        Facets = as.character(error_tbl$Facets %||% NA_character_),
        Rows = 0L,
        FlaggedRows = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Source = character(0),
        Used = logical(0),
        Reason = character(0),
        Facets = character(0),
        Rows = integer(0),
        FlaggedRows = integer(0),
        stringsAsFactors = FALSE
      )
    }
    return(list(
      levels = data.frame(Level = character(0), BiasCount = integer(0), BiasSources = integer(0), stringsAsFactors = FALSE),
      sources = source_tbl
    ))
  }

  level_counts <- list()
  source_rows <- list()

  for (i in seq_along(sources)) {
    source_name <- names(sources)[i]
    bundle <- sources[[i]]
    spec <- extract_bias_facet_spec(bundle)
    reason <- ""
    used <- FALSE
    rows_n <- 0L
    flagged_n <- 0L
    facets_label <- ""

    if (is.null(spec) || length(spec$facets) < 2) {
      reason <- "unrecognized bundle"
    } else if (!target_facet %in% spec$facets) {
      reason <- "target facet not involved"
      facets_label <- paste(spec$facets, collapse = " x ")
    } else {
      tbl <- as.data.frame(bundle$table, stringsAsFactors = FALSE)
      level_col <- spec$level_cols[match(target_facet, spec$facets)]
      if (is.na(level_col) || !nzchar(level_col) || !level_col %in% names(tbl)) {
        reason <- "target level column missing"
        facets_label <- paste(spec$facets, collapse = " x ")
      } else {
        used <- TRUE
        facets_label <- paste(spec$facets, collapse = " x ")
        tbl[[level_col]] <- as.character(tbl[[level_col]])
        rows_n <- nrow(tbl)
        bias_size <- numeric_or_na(tbl, c("Bias Size", "BiasSize"))
        t_val <- numeric_or_na(tbl, c("t", "AbsT"))
        p_val <- numeric_or_na(tbl, c("Prob.", "Prob"))
        if ("Flag" %in% names(tbl)) {
          row_flag <- as.logical(tbl$Flag)
        } else {
          row_flag <- (is.finite(t_val) & abs(t_val) >= abs(bias_abs_t_warn)) |
            (is.finite(bias_size) & abs(bias_size) >= abs(bias_abs_size_warn)) |
            (is.finite(p_val) & p_val <= bias_p_max)
        }
        flagged_n <- sum(row_flag, na.rm = TRUE)
        lvl <- tibble::as_tibble(tbl[, c(level_col), drop = FALSE])
        names(lvl) <- "Level"
        agg <- lvl |>
          dplyr::mutate(Flag = row_flag) |>
          dplyr::group_by(.data$Level) |>
          dplyr::summarise(
            BiasCount = sum(.data$Flag, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::mutate(BiasSources = as.integer(.data$BiasCount > 0))
        level_counts[[length(level_counts) + 1L]] <- as.data.frame(agg, stringsAsFactors = FALSE)
      }
    }

    source_rows[[length(source_rows) + 1L]] <- data.frame(
      Source = source_name,
      Used = used,
      Reason = reason,
      Facets = facets_label,
      Rows = rows_n,
      FlaggedRows = flagged_n,
      stringsAsFactors = FALSE
    )
  }

  level_tbl <- if (length(level_counts) == 0) {
    data.frame(Level = character(0), BiasCount = integer(0), BiasSources = integer(0), stringsAsFactors = FALSE)
  } else {
    merged <- dplyr::bind_rows(level_counts)
    merged <- merged |>
      dplyr::group_by(.data$Level) |>
      dplyr::summarise(
        BiasCount = sum(.data$BiasCount, na.rm = TRUE),
        BiasSources = sum(.data$BiasSources, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data$BiasCount), .data$Level)
    as.data.frame(merged, stringsAsFactors = FALSE)
  }

  source_tbl <- dplyr::bind_rows(source_rows)
  if (is.data.frame(error_tbl) && nrow(error_tbl) > 0) {
    source_tbl <- dplyr::bind_rows(
      source_tbl,
      data.frame(
        Source = as.character(error_tbl$Interaction %||% paste0("bias_error_", seq_len(nrow(error_tbl)))),
        Used = FALSE,
        Reason = paste0("pair error: ", as.character(error_tbl$Error)),
        Facets = as.character(error_tbl$Facets %||% NA_character_),
        Rows = 0L,
        FlaggedRows = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }
  list(levels = level_tbl, sources = source_tbl)
}

numeric_or_na <- function(tbl, candidates) {
  for (nm in candidates) {
    if (nm %in% names(tbl)) return(suppressWarnings(as.numeric(tbl[[nm]])))
  }
  rep(NA_real_, nrow(tbl))
}

dashboard_draw_plot <- function(tbl,
                                plot_type,
                                facet,
                                thresholds,
                                main = NULL,
                                palette = NULL,
                                label_angle = 45) {
  pal <- stats::setNames(
    c("#2b8cbe", "#cb181d", "#756bb1", "#bdbdbd", "#238b45"),
    c("neutral", "flag", "bias", "central", "misfit")
  )
  if (!is.null(palette) && length(palette) > 0) {
    nm <- intersect(names(palette), names(pal))
    pal[nm] <- as.character(palette[nm])
  }

  labels <- utils::head(as.character(tbl$Level), n = nrow(tbl))
  if (plot_type == "severity") {
    ord <- order(abs(tbl$Estimate), decreasing = TRUE, na.last = NA)
    tbl <- tbl[ord, , drop = FALSE]
    labels <- utils::head(as.character(tbl$Level), n = nrow(tbl))
    cols <- rep(unname(pal["neutral"]), nrow(tbl))
    cols[tbl$CentralTendencyFlag %in% TRUE] <- unname(pal["central"])
    cols[tbl$MisfitFlag %in% TRUE] <- unname(pal["misfit"])
    cols[tbl$SeverityFlag %in% TRUE] <- unname(pal["flag"])
    y <- rev(seq_len(nrow(tbl)))
    plot(
      x = rev(tbl$Estimate),
      y = y,
      type = "n",
      yaxt = "n",
      xlab = "Estimate (logits)",
      ylab = "",
      main = ""
    )
    segments(0, y, rev(tbl$Estimate), y, col = "gray70")
    points(rev(tbl$Estimate), y, pch = 16, col = rev(cols))
    axis(side = 2, at = y, labels = rev(labels), las = 2, cex.axis = 0.75)
    abline(v = c(-abs(thresholds$severity_warn[1]), 0, abs(thresholds$severity_warn[1])),
                     lty = c(2, 1, 2), col = c("gray45", "gray30", "gray45"))
    return(invisible(NULL))
  }

  tbl <- tbl[tbl$AnyFlag %in% TRUE, , drop = FALSE]
  if (nrow(tbl) == 0) {
    plot.new()
    text(0.5, 0.5, "No observed flags; inspect missing diagnostics and notes")
    return(invisible(NULL))
  }
  ord <- order(tbl$FlagCount, abs(tbl$Estimate), decreasing = TRUE, na.last = TRUE)
  tbl <- tbl[ord, , drop = FALSE]
  labels <- utils::head(as.character(tbl$Level), n = nrow(tbl))
  cols <- ifelse(tbl$BiasFlag %in% TRUE, pal["bias"], pal["flag"])
  bp <- barplot(
    height = tbl$FlagCount,
    names.arg = labels,
    col = cols,
    las = if (isTRUE(label_angle >= 45)) 2 else 1,
    main = "",
    ylab = "Flag count"
  )
  abline(h = 0, col = "gray60")
  invisible(bp)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' @rdname facet_quality_dashboard
#' @details `summary()` method: Summarize a facet-quality dashboard
#' @param object Output from [facet_quality_dashboard()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of flagged levels to preview.
#' @return For `summary()`: An object of class `summary.mfrm_facet_dashboard`.
#' @seealso [facet_quality_dashboard()], [plot(x, type = "facet_dashboard")][plot.mfrm_fit]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("summary.mfrm_facet_dashboard", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' diag <- example_fit$diagnostics
#' summary(facet_quality_dashboard(fit, diagnostics = diag))
#' }
#' @export
summary.mfrm_facet_dashboard <- function(object, digits = 3, top_n = 10, ...) {
  if (!is.list(object) || is.null(object$detail)) {
    stop("`object` must be output from facet_quality_dashboard().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  flagged <- as.data.frame(object$flagged %||% data.frame(), stringsAsFactors = FALSE)
  preview <- utils::head(flagged, n = top_n)
  settings <- as.data.frame(object$settings %||% data.frame(), stringsAsFactors = FALSE)
  bias_sources <- as.data.frame(object$bias_sources %||% data.frame(), stringsAsFactors = FALSE)

  out <- list(
    summary_kind = "facet_dashboard",
    overview = overview,
    summary = summary_tbl,
    preview_name = "flagged",
    preview = preview,
    settings = settings,
    bias_sources = bias_sources,
    notes = object$notes %||% character(0),
    digits = digits
  )
  class(out) <- "summary.mfrm_facet_dashboard"
  out
}

#' @export
print.summary.mfrm_facet_dashboard <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrmr Facet Quality Dashboard Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(dashboard_round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$summary) && nrow(x$summary) > 0) {
    cat("\nSummary\n")
    print(dashboard_round_numeric_df(as.data.frame(x$summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$preview) && nrow(x$preview) > 0) {
    cat("\nFlagged levels\n")
    print(dashboard_round_numeric_df(as.data.frame(x$preview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$bias_sources) && nrow(x$bias_sources) > 0) {
    cat("\nBias sources\n")
    print(dashboard_round_numeric_df(as.data.frame(x$bias_sources), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$settings) && nrow(x$settings) > 0) {
    cat("\nSettings\n")
    print(as.data.frame(x$settings), row.names = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Plot a facet-quality dashboard
#'
#' @param x Output from [facet_quality_dashboard()] or [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is a fit.
#' @param facet Optional facet name.
#' @param bias_results Optional bias bundle or list of bundles.
#' @param severity_warn Absolute estimate cutoff used to flag severity
#'   outliers.
#' @param misfit_warn Mean-square cutoff used to flag misfit.
#' @param central_tendency_max Legacy opt-in absolute estimate cutoff for
#'   marking facet estimates near the fitted origin. Default `NULL`; use
#'   [data_quality_report()] for observed category-use screening or
#'   [rater_centrality_table()] for indirect centrality indices.
#' @param bias_count_warn Minimum flagged-bias row count required to flag a
#'   level.
#' @param bias_abs_t_warn Absolute `t` cutoff used when deriving bias-row
#'   flags from a raw bias bundle.
#' @param bias_abs_size_warn Absolute bias-size cutoff used when deriving
#'   bias-row flags from a raw bias bundle.
#' @param bias_p_max Probability cutoff used when deriving bias-row flags
#'   from a raw bias bundle.
#' @param plot_type Plot type, `"severity"` or `"flags"`.
#' @param top_n Number of rows to keep in the plot data.
#' @inheritParams plot_marginal_fit
#' @param palette Optional named color overrides.
#' @param label_angle Label angle hint for the `"flags"` plot.
#' @param draw If `TRUE`, draw with base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details When `x` is a dashboard, its stored screening settings define both
#' flags and plot guides; threshold arguments apply only when `x` is a fit.
#' Returned plot data retain those settings, interpretation notes and fit
#' readiness. Restricted or older dashboards without readiness are labeled
#' `REVIEW ONLY`. Zero observed flags is not a complete diagnostic pass.
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [facet_quality_dashboard()], [summary.mfrm_facet_dashboard()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_facet_quality_dashboard", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' diag <- example_fit$diagnostics
#' p <- plot(fit, type = "facet_dashboard", diagnostics = diag, draw = FALSE)
#' p$data$plot
#' }
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_facet_quality_dashboard <- function(x,
                                         diagnostics = NULL,
                                         facet = NULL,
                                         bias_results = NULL,
                                         severity_warn = 1.0,
                                         misfit_warn = NULL,
                                         central_tendency_max = NULL,
                                         bias_count_warn = 1L,
                                         bias_abs_t_warn = 2,
                                         bias_abs_size_warn = 0.5,
                                         bias_p_max = 0.05,
                                         plot_type = c("severity", "flags"),
                                         top_n = 20,
                                         main = NULL,
                                         palette = NULL,
                                         label_angle = 45,
                                         draw = TRUE,
                                         ..., title = NULL) {
  main <- .mfrm_plot_title_alias(main, title, missing(main), missing(title))
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("severity", "flags"))
  top_n <- max(1L, as.integer(top_n))

  bundle <- if (inherits(x, "mfrm_facet_dashboard")) {
    x
  } else if (inherits(x, "mfrm_fit")) {
    facet_quality_dashboard(
      fit = x,
      diagnostics = diagnostics,
      facet = facet,
      bias_results = bias_results,
      severity_warn = severity_warn,
      misfit_warn = misfit_warn,
      central_tendency_max = central_tendency_max,
      bias_count_warn = bias_count_warn,
      bias_abs_t_warn = bias_abs_t_warn,
      bias_abs_size_warn = bias_abs_size_warn,
      bias_p_max = bias_p_max
    )
  } else {
    stop("`x` must be an mfrm_fit object or a facet dashboard bundle.", call. = FALSE)
  }

  # A saved dashboard's settings generated its flags and must also define its guides.
  threshold_names <- c("severity_warn", "misfit_warn", "misfit_lower",
                       "central_tendency_max", "bias_count_warn", "bias_abs_t_warn",
                       "bias_abs_size_warn", "bias_p_max")
  thresholds <- lapply(threshold_names, function(key) {
    suppressWarnings(as.numeric(bundle$settings$Value[match(key, bundle$settings$Setting)]))
  })
  names(thresholds) <- threshold_names
  plot_title <- main %||% paste0("Facet quality: ", bundle$facet)
  if ((missing(title) || !identical(main, "")) &&
      !identical(bundle$interpretation_status, "ready_for_diagnostic_interpretation")) {
    plot_title <- paste("REVIEW ONLY -", plot_title)
  }

  tbl <- as.data.frame(bundle$detail, stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) {
    stop("Facet dashboard does not contain any level rows.", call. = FALSE)
  }

  if (plot_type == "severity") {
    tbl <- tbl[is.finite(tbl$Estimate), , drop = FALSE]
    if (nrow(tbl) == 0) stop("No finite severity estimates are available.", call. = FALSE)
    ord <- order(abs(tbl$Estimate), decreasing = TRUE, na.last = NA)
    tbl <- tbl[ord, , drop = FALSE]
    if (nrow(tbl) > top_n) tbl <- tbl[seq_len(top_n), , drop = FALSE]
  } else {
    tbl <- tbl[tbl$AnyFlag %in% TRUE, , drop = FALSE]
    if (nrow(tbl) > 0) {
      ord <- order(tbl$FlagCount, abs(tbl$Estimate), decreasing = TRUE, na.last = TRUE)
      tbl <- tbl[ord, , drop = FALSE]
      if (nrow(tbl) > top_n) tbl <- tbl[seq_len(top_n), , drop = FALSE]
    }
  }

  if (isTRUE(draw)) {
    opar <- graphics::par()["mar"]
    on.exit(graphics::par(opar), add = TRUE)
    graphics::par(mar = c(4.5, 8.5, 3.5, 1.5))
    dashboard_draw_plot(
      tbl = tbl,
      plot_type = plot_type,
      facet = bundle$facet %||% bundle$overview$Facet[1],
      thresholds = thresholds,
      main = plot_title,
      palette = palette,
      label_angle = label_angle
    )
  }

  out <- new_mfrm_plot_data(
    "facet_quality_dashboard",
    list(
      plot = plot_type,
      facet = bundle$facet %||% bundle$overview$Facet[1],
      table = tbl,
      summary = bundle$summary,
      overview = bundle$overview,
      settings = bundle$settings,
      flagged = bundle$flagged,
      ranked = bundle$ranked,
      thresholds = thresholds,
      title = plot_title,
      notes = data.frame(Type = "Interpretation", Text = bundle$notes),
      fit_readiness = bundle$fit_readiness,
      interpretation_status = bundle$interpretation_status
    )
  )
  invisible(out)
}

#' @rdname facet_quality_dashboard
#' @details `plot()` method: Plot a saved facet-quality dashboard
#' @param x A result from [facet_quality_dashboard()].
#' @param y Ignored; present for generic compatibility.
#' @param type `"severity"` or `"flags"`.
#' @param ... For `plot()`: Display controls from [plot_facet_quality_dashboard()].
#'
#' For `summary()`: Reserved for generic compatibility.
#' @return For `plot()`: The existing plot payload, including screening cautions.
#' @seealso [mfrmr_plot_migration]
#' @export
plot.mfrm_facet_dashboard <- function(x, y = NULL, type = c("severity", "flags"), ...) {
  args <- list(x = x, ...)
  if (identical(type, "facet_dashboard")) {
    type <- NULL
    args <- .mfrm_plot_rename(args, "view", "plot_type")
  }
  if (!is.null(type)) args$plot_type <- type
  do.call(plot_facet_quality_dashboard, args)
}
