# Reporting and narrative helpers (package-native / APA-style)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

py_style_format <- function(fmt, value) {
  if (is.na(value)) return("")
  if (is.null(fmt)) return(as.character(value))
  if (is.function(fmt)) return(as.character(fmt(value)))
  if (!is.character(fmt) || length(fmt) != 1) return(as.character(value))
  if (identical(fmt, "{}")) return(as.character(value))

  m <- regexec("^\\{:\\.(\\d+)f\\}$", fmt)
  g <- regmatches(fmt, m)[[1]]
  if (length(g) == 2) {
    digits <- as.integer(g[2])
    return(sprintf(paste0("%.", digits, "f"), as.numeric(value)))
  }
  as.character(value)
}

format_fixed_width_table <- function(df,
                                     columns,
                                     formats = list(),
                                     right_align = NULL,
                                     max_col_width = 16,
                                     min_col_width = 6) {
  if (is.null(df) || nrow(df) == 0) return("No data")

  if (is.null(right_align)) {
    right_align <- columns[columns %in% names(df)][vapply(columns[columns %in% names(df)], function(col) {
      is.numeric(df[[col]])
    }, logical(1))]
  }

  str_cols <- list()
  widths <- list()

  for (col in columns) {
    if (!col %in% names(df)) {
      vals <- rep("", nrow(df))
    } else {
      vals <- vapply(df[[col]], function(v) {
        if (is.na(v)) return("")
        py_style_format(formats[[col]], v)
      }, character(1))
    }

    str_cols[[col]] <- vals
    max_len <- max(c(nchar(col), nchar(vals)), na.rm = TRUE)
    widths[[col]] <- max(min_col_width, min(max_len, max_col_width))
  }

  pad <- function(col, text) {
    text <- substr(text, 1, widths[[col]])
    if (col %in% right_align) {
      stringr::str_pad(text, widths[[col]], side = "left")
    } else {
      stringr::str_pad(text, widths[[col]], side = "right")
    }
  }

  header <- paste(vapply(columns, function(col) pad(col, col), character(1)), collapse = " ")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    paste(vapply(columns, function(col) pad(col, str_cols[[col]][i]), character(1)), collapse = " ")
  }, character(1))

  paste(c(header, rows), collapse = "\n")
}

build_sectioned_fixed_report <- function(title = NULL,
                                         sections = list(),
                                         max_col_width = 18,
                                         min_col_width = 6) {
  lines <- character(0)
  if (!is.null(title) && nzchar(as.character(title[1]))) {
    lines <- c(lines, as.character(title[1]))
  }

  for (i in seq_along(sections)) {
    section <- sections[[i]]
    section_title <- section$title %||% names(sections)[i] %||% paste0("Section ", i)
    section_data <- section$data %||% section$table
    section_columns <- section$columns %||% if (is.data.frame(section_data)) names(section_data) else NULL
    section_formats <- section$formats %||% list()
    section_max_col_width <- as.integer(section$max_col_width %||% max_col_width)
    section_min_col_width <- as.integer(section$min_col_width %||% min_col_width)
    section_max_rows <- suppressWarnings(as.integer(section$max_rows %||% NA_integer_))
    total_rows <- if (is.data.frame(section_data)) nrow(section_data) else NA_integer_

    if (length(lines) > 0) {
      lines <- c(lines, "")
    }
    if (!is.null(section_title) && nzchar(as.character(section_title[1]))) {
      lines <- c(lines, as.character(section_title[1]))
    }

    if (is.null(section_data)) {
      lines <- c(lines, "No data")
      next
    }

    if (is.character(section_data) && length(section_data) > 0) {
      lines <- c(lines, as.character(section_data))
      next
    }

    if (!is.data.frame(section_data)) {
      lines <- c(lines, as.character(section_data))
      next
    }

    if (nrow(section_data) == 0) {
      lines <- c(lines, "No data")
      next
    }

    table_for_print <- section_data
    if (is.finite(section_max_rows) && section_max_rows > 0 && nrow(table_for_print) > section_max_rows) {
      table_for_print <- table_for_print[seq_len(section_max_rows), , drop = FALSE]
    }

    lines <- c(
      lines,
      format_fixed_width_table(
        df = table_for_print,
        columns = section_columns,
        formats = section_formats,
        max_col_width = section_max_col_width,
        min_col_width = section_min_col_width
      )
    )
    if (is.finite(section_max_rows) && section_max_rows > 0 &&
        is.finite(total_rows) && total_rows > section_max_rows) {
      lines <- c(lines, paste0("Showing first ", section_max_rows, " rows of ", total_rows, "."))
    }
  }

  paste(lines, collapse = "\n")
}

build_bias_fixed_text <- function(table_df,
                                  summary_df,
                                  chi_df,
                                  facet_a,
                                  facet_b,
                                  interaction_label = NULL,
                                  columns,
                                  formats) {
  if (is.null(table_df) || nrow(table_df) == 0) return("No bias data")

  fixed_table <- format_fixed_width_table(table_df, columns, formats = formats, max_col_width = 18)
  label <- if (!is.null(interaction_label) && nzchar(interaction_label)) {
    interaction_label
  } else if (!is.null(facet_b) && nzchar(as.character(facet_b))) {
    paste0(facet_a, " x ", facet_b)
  } else {
    as.character(facet_a)
  }
  lines <- c(paste0("Bias/Interaction: ", label), "", fixed_table)

  if (!is.null(summary_df) && nrow(summary_df) > 0) {
    lines <- c(lines, "", "Summary")
    for (i in seq_len(nrow(summary_df))) {
      row <- summary_df[i, , drop = FALSE]
      line <- paste0(
        as.character(row$Statistic),
        ": Bias Size=", fmt_num(row$`Bias Size`, 2),
        ", Obs-Exp Avg=", fmt_num(row$`Obs-Exp Average`, 2),
        ", Model S.E.=", fmt_num(row$`S.E.`, 2)
      )
      lines <- c(lines, line)
    }
  }

  if (!is.null(chi_df) && nrow(chi_df) > 0) {
    chi <- chi_df[1, , drop = FALSE]
    line <- if (!is.na(chi$FixedChiSq)) {
      paste0(
        "Fixed (all = 0) chi-squared: ", sprintf("%.2f", as.numeric(chi$FixedChiSq)),
        "  d.f.: ", ifelse(is.na(chi$FixedDF), "", as.character(as.integer(round(chi$FixedDF)))),
        "  significance (probability): ", sprintf("%.4f", as.numeric(chi$FixedProb))
      )
    } else {
      "Fixed (all = 0) chi-squared: N/A"
    }
    lines <- c(lines, "", line)
  }

  paste(lines, collapse = "\n")
}

build_pairwise_fixed_text <- function(pair_df, target_facet, context_facet, columns, formats) {
  if (is.null(pair_df) || nrow(pair_df) == 0) return("No pairwise data")
  fixed_table <- format_fixed_width_table(pair_df, columns, formats = formats, max_col_width = 18)
  paste(
    paste0("Bias/Interaction Pairwise Report: Target=", target_facet, "  Context=", context_facet),
    "",
    fixed_table,
    sep = "\n"
  )
}

to_float <- function(value) {
  out <- suppressWarnings(as.numeric(value))
  ifelse(length(out) == 0, NA_real_, out)
}

fmt_count <- function(value) {
  val <- to_float(value)
  if (!is.finite(val)) return("NA")
  if (abs(val - round(val)) < 1e-6) return(as.character(as.integer(round(val))))
  sprintf("%.0f", val)
}

fmt_num <- function(value, decimals = 2) {
  val <- to_float(value)
  if (!is.finite(val)) return("NA")
  sprintf(paste0("%.", decimals, "f"), val)
}

fmt_pvalue <- function(value) {
  val <- to_float(value)
  if (!is.finite(val)) return("NA")
  if (val < 0.001) return("< .001")
  paste0("= ", sprintf("%.3f", val))
}

describe_series <- function(series) {
  if (is.null(series)) return(NULL)
  arr <- suppressWarnings(as.numeric(series))
  arr <- arr[is.finite(arr)]
  if (length(arr) == 0) return(NULL)
  list(
    min = min(arr),
    max = max(arr),
    mean = mean(arr),
    sd = if (length(arr) > 1) stats::sd(arr) else NA_real_
  )
}

safe_residual_pca <- function(diagnostics, mode = "both", pca_max_factors = 10L) {
  if (is.null(diagnostics)) return(NULL)
  tryCatch(
    {
      # Reporting summarizes stored diagnostics; it must not add an analysis.
      available_mode <- infer_export_residual_pca_mode(diagnostics)
      if (available_mode == "none" ||
          (mode != "both" && available_mode != "both" && mode != available_mode)) {
        return(NULL)
      }
      if (mode == "both") mode <- available_mode
      stored <- c(
        if (mode %in% c("overall", "both")) list(diagnostics$residual_pca_overall),
        if (mode %in% c("facet", "both")) diagnostics$residual_pca_by_facet
      )
      stored <- Filter(Negate(is.null), stored)
      if (!all(vapply(stored, function(x) identical(x$calculation_version, 2L), logical(1)))) {
        stop("Recreate the stored residual PCA with diagnose_mfrm(fit, residual_pca = ...) before reporting it.")
      }
      pca_max_factors <- stored[[1]]$max_factors
      if (!all(vapply(stored, function(x) identical(x$max_factors, pca_max_factors), logical(1)))) {
        stop("The stored residual PCA analyses have different component limits. Recreate them together before reporting.")
      }
      analyze_residual_pca(
        diagnostics = diagnostics,
        mode = mode,
        pca_max_factors = pca_max_factors
      )
    },
    error = function(e) {
      structure(
        list(
          overall_table = data.frame(),
          by_facet_table = data.frame(),
          overall = NULL,
          by_facet = list(),
          errors = list(
            overall = conditionMessage(e),
            by_facet = data.frame(Facet = character(0), Error = character(0), stringsAsFactors = FALSE)
          )
        ),
        class = c("mfrm_bundle", "mfrm_residual_pca", "list")
      )
    }
  )
}

extract_diagnostic_basis_status <- function(diagnostics, diagnostic_path) {
  basis <- as.data.frame(diagnostics$diagnostic_basis %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(basis) == 0L || !all(c("DiagnosticPath", "Status") %in% names(basis))) {
    return(NA_character_)
  }
  idx <- which(as.character(basis$DiagnosticPath) == diagnostic_path)
  if (length(idx) == 0L) return(NA_character_)
  as.character(basis$Status[idx[1]])
}

marginal_coverage_notes <- function(coverage) {
  if (is.null(coverage) || !nrow(coverage)) return(character(0))
  vapply(seq_len(nrow(coverage)), function(i) {
    r <- coverage[i, ]
    if (!r$Classified) {
      sprintf("%s: no classifications available; %d of %d unavailable.", r$Screen, r$Unclassified, r$Total)
    } else {
      sprintf("%s: %d flagged among %d classified; %d of %d unavailable.",
              r$Screen, r$Flagged, r$Classified, r$Unclassified, r$Total)
    }
  }, character(1))
}

extract_strict_marginal_visual_state <- function(diagnostics) {
  marginal_fit <- diagnostics$marginal_fit %||% NULL
  validate_marginal_coverage(marginal_fit)
  marginal_summary <- as.data.frame(marginal_fit$summary %||% data.frame(), stringsAsFactors = FALSE)
  marginal_available <- is.list(marginal_fit) && isTRUE(marginal_fit$available)
  pairwise_bundle <- if (is.list(marginal_fit)) marginal_fit$pairwise %||% NULL else NULL
  pairwise_available <- is.list(pairwise_bundle) && isTRUE(pairwise_bundle$available)

  marginal_status <- extract_diagnostic_basis_status(diagnostics, "strict_marginal_fit")
  pairwise_status <- extract_diagnostic_basis_status(diagnostics, "strict_pairwise_local_dependence")
  if (!is.finite(match(marginal_status, c("computed", "requested_not_available", "not_requested", "available_but_not_requested", "not_available_for_run")))) {
    marginal_status <- if (marginal_available) "computed" else "not_requested"
  }
  if (!is.finite(match(pairwise_status, c("computed", "requested_not_available", "not_requested", "available_but_not_requested", "not_available_for_run")))) {
    pairwise_status <- if (pairwise_available) "computed" else "not_requested"
  }

  top_cell <- {
    top_cells <- as.data.frame(marginal_fit$top_cells %||% data.frame(), stringsAsFactors = FALSE)
    top_cells <- top_cells[is.finite(top_cells$StdResidual), , drop = FALSE]
    if (nrow(top_cells) > 0L) top_cells[1, , drop = FALSE] else data.frame()
  }
  top_pair <- {
    top_pairs <- as.data.frame(pairwise_bundle$top_pairs %||% data.frame(), stringsAsFactors = FALSE)
    top_pairs <- top_pairs[is.finite(top_pairs$ExactStdResidual) | is.finite(top_pairs$AdjacentStdResidual), , drop = FALSE]
    if (nrow(top_pairs) > 0L) top_pairs[1, , drop = FALSE] else data.frame()
  }

  marginal_reason <- if (marginal_available) {
    ""
  } else if (nrow(marginal_summary) > 0L && "Reason" %in% names(marginal_summary)) {
    as.character(marginal_summary$Reason[1] %||% "")
  } else if (identical(marginal_status, "not_requested") || identical(marginal_status, "available_but_not_requested")) {
    "Rerun diagnose_mfrm(..., diagnostic_mode = \"both\") to compute strict marginal diagnostics."
  } else {
    "Strict marginal diagnostics were not available for this run."
  }

  pairwise_reason <- if (pairwise_available) {
    ""
  } else if (is.list(pairwise_bundle) && nzchar(as.character(pairwise_bundle$reason %||% ""))) {
    as.character(pairwise_bundle$reason)
  } else if (identical(pairwise_status, "not_requested") || identical(pairwise_status, "available_but_not_requested")) {
    "Rerun diagnose_mfrm(..., diagnostic_mode = \"both\") to compute strict pairwise local-dependence diagnostics."
  } else {
    "Strict pairwise local-dependence diagnostics were not available for this run."
  }

  list(
    coverage = marginal_fit$coverage %||% data.frame(),
    coverage_notes = marginal_coverage_notes(marginal_fit$coverage),
    marginal_status = marginal_status,
    marginal_available = marginal_available,
    marginal_reason = marginal_reason,
    overall_rmsd = if (nrow(marginal_summary) > 0L) to_float(marginal_summary$OverallRMSD[1]) else NA_real_,
    overall_max_abs_std_residual = if (nrow(marginal_summary) > 0L) to_float(marginal_summary$OverallMaxAbsStdResidual[1]) else NA_real_,
    step_groups_flagged = if (nrow(marginal_summary) > 0L) to_float(marginal_summary$StepGroupsFlagged[1]) else NA_real_,
    facet_levels_flagged = if (nrow(marginal_summary) > 0L) to_float(marginal_summary$FacetLevelsFlagged[1]) else NA_real_,
    top_cell = top_cell,
    marginal_thresholds = marginal_fit$thresholds %||% list(),
    pairwise_status = pairwise_status,
    pairwise_available = pairwise_available,
    pairwise_reason = pairwise_reason,
    pairwise_flagged_pairs = if (nrow(marginal_summary) > 0L) to_float(marginal_summary$PairwiseFlaggedLevelPairs[1]) else NA_real_,
    top_pair = top_pair,
    pairwise_thresholds = pairwise_bundle$thresholds %||% list()
  )
}

format_reporting_marginal_cell_label <- function(cell_row) {
  cell_row <- as.data.frame(cell_row %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(cell_row) == 0L) return("")
  format_marginal_cell_label(
    cell_type = cell_row$CellType[1] %||% "",
    step_facet = cell_row$StepFacet[1] %||% "",
    facet = cell_row$Facet[1] %||% "",
    level = cell_row$Level[1] %||% "",
    category = cell_row$Category[1] %||% ""
  )
}

format_reporting_marginal_pair_label <- function(pair_row) {
  pair_row <- as.data.frame(pair_row %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(pair_row) == 0L) return("")
  format_marginal_pair_label(
    facet = pair_row$Facet[1] %||% "",
    level1 = pair_row$Level1[1] %||% "",
    level2 = pair_row$Level2[1] %||% ""
  )
}

#' MnSq misfit threshold pair used across mfrmr screening helpers
#'
#' Returns the lower / upper bounds that mfrmr screens treat as the
#' heuristic mean-square (Infit / Outfit MnSq) review band when flagging
#' element-level misfit. Defaults use the published 0.5-1.5 interval as a
#' screening convention; both ends can be overridden via R options. The
#' interval is not a universal acceptance rule or an automatic exclusion rule.
#' Low mean squares describe low residual variability rather than poor rater
#' quality. [fit_measures_table()] uses this mean-square screen by default;
#' ZSTD-only evidence is separate unless its combined rule is requested.
#' Compare declared bands against known simulation truth with
#' [mfrm_screening_sensitivity()].
#'
#' Helpers that consume the band include
#' [summary.mfrm_diagnostics()] (`misfit_flagged` block and
#' `key_warnings` auto-flag), [build_misfit_casebook()] (the new
#' `element_fit` source family), the bias / misfit narrative inside
#' [build_apa_outputs()], and [facet_quality_dashboard()] when
#' `misfit_warn = NULL`. Setting the options once at the top of an
#' analysis script therefore changes every downstream screen at once.
#'
#' @section Configuration:
#' Two scalar R options drive the band:
#' \describe{
#'   \item{`mfrmr.misfit_lower`}{Lower review bound. Default `0.5`.}
#'   \item{`mfrmr.misfit_upper`}{Upper review bound. Default `1.5`.}
#' }
#'
#' Pass scalar arguments to override the options for a single call,
#' e.g. `mfrm_misfit_thresholds(lower = 0.7, upper = 1.3)` for the
#' tighter Bond & Fox (2015) reporting band.
#'
#' @param lower Optional lower bound. When `NULL` (default), the
#'   active option / package default is used.
#' @param upper Optional upper bound.
#'
#' @return A named numeric vector `c(lower = ..., upper = ...)` with
#'   `lower < upper`.
#' @seealso [summary.mfrm_diagnostics()], [build_misfit_casebook()],
#'   [facet_quality_dashboard()]
#' @examples
#' mfrm_misfit_thresholds()
#' old <- options(mfrmr.misfit_lower = 0.7, mfrmr.misfit_upper = 1.3)
#' mfrm_misfit_thresholds()
#' options(old)
#' @export
mfrm_misfit_thresholds <- function(lower = NULL, upper = NULL) {
  lo <- if (is.null(lower)) {
    suppressWarnings(as.numeric(getOption("mfrmr.misfit_lower", 0.5)))
  } else {
    suppressWarnings(as.numeric(lower[1]))
  }
  up <- if (is.null(upper)) {
    suppressWarnings(as.numeric(getOption("mfrmr.misfit_upper", 1.5)))
  } else {
    suppressWarnings(as.numeric(upper[1]))
  }
  if (!is.finite(lo) || !is.finite(up) || lo <= 0 || up <= 0 || lo >= up) {
    stop("`mfrm_misfit_thresholds` requires `0 < lower < upper`. ",
         "Got lower = ", lo, ", upper = ", up, ".",
         call. = FALSE)
  }
  c(lower = lo, upper = up)
}

fit_screening_summary <- function(tbl, band = mfrm_misfit_thresholds()) {
  n <- nrow(tbl) %||% 0L
  column <- function(name) suppressWarnings(as.numeric(tbl[[name]] %||% rep(NA_real_, n)))
  infit <- column("Infit"); outfit <- column("Outfit")
  inz <- column("InfitZSTD"); outz <- column("OutfitZSTD")
  infit[!is.finite(infit) | infit < 0] <- NA_real_
  outfit[!is.finite(outfit) | outfit < 0] <- NA_real_
  inz[!is.finite(inz)] <- NA_real_; outz[!is.finite(outz)] <- NA_real_
  complete_mnsq <- !is.na(infit) & !is.na(outfit)
  complete_z <- !is.na(inz) & !is.na(outz)
  flags <- data.frame(
    MeanSquare = infit < band["lower"] | infit > band["upper"] |
      outfit < band["lower"] | outfit > band["upper"],
    ZSTD2 = abs(inz) >= 2 | abs(outz) >= 2,
    ZSTD3 = abs(inz) >= 3 | abs(outz) >= 3
  )
  counts <- do.call(rbind, lapply(seq_along(flags), function(i) {
    flag <- flags[[i]]
    available <- sum(!is.na(flag))
    data.frame(
      Screen = c(sprintf("MnSq outside [%g, %g]", band["lower"], band["upper"]), "|ZSTD| >= 2", "|ZSTD| >= 3")[i],
      Elements = n, Classified = available, Unclassified = sum(is.na(flag)),
      IncompleteStatistics = sum(!(if (i == 1L) complete_mnsq else complete_z)),
      Flagged = if (available) sum(flag, na.rm = TRUE) else NA_integer_,
      FlagRate = if (n > 0L && available == n) mean(flag) else NA_real_
    )
  }))
  list(flags = flags, counts = counts,
       max_abs_z = pmax(abs(inz), abs(outz)),
       max_mnsq_deviation = pmax(abs(infit - 1), abs(outfit - 1)))
}

fit_screening_note <- function(row) {
  if (row$Classified == 0L) {
    return(sprintf("%s: no elements could be classified (%d elements total).", row$Screen, row$Elements))
  }
  text <- sprintf("%s: %d of %d classified elements flagged; %d of %d elements unclassified.",
                  row$Screen, row$Flagged, row$Classified, row$Unclassified, row$Elements)
  if (row$IncompleteStatistics > 0L) {
    text <- paste(text, sprintf("%d elements have an unavailable statistic.", row$IncompleteStatistics))
  }
  text
}

# Warning threshold profiles for MFRM quality control.
# Sources:
#   - n_obs_min, n_person_min: Linacre (1994), sample size guidelines for stable estimates.
#   - low_cat_min: Linacre (2002), minimum 10 observations per category for stable thresholds.
#   - misfit_ratio_warn: Bond & Fox (2015), MnSq 0.5-1.5 acceptable range; >10% flagged.
#   - zstd thresholds: |ZSTD| > 2 and >3 are heuristic normal-reference bands,
#     not calibrated 5% and 1% tests (Wright & Linacre, 1994).
#   - pca_first_eigen_warn: Linacre-style residual-PCA heuristic band; use only as exploratory screening, not direct proof of multidimensionality.
#   - pca_first_prop_warn: Smith (2002), unexplained variance > 5-10% merits investigation.
#   - pca_reference_bands: uncalibrated descriptive reference values only.
warning_threshold_profiles <- function() {
  list(
    profiles = list(
      strict = list(
        n_obs_min = 200,
        n_person_min = 50,
        low_cat_min = 15,
        min_facet_levels = 4,
        misfit_ratio_warn = 0.08,
        missing_fit_ratio_warn = 0.15,
        zstd2_ratio_warn = 0.08,
        zstd3_ratio_warn = 0.03,
        expected_var_min = 0.30,
        pca_first_eigen_warn = 1.5,
        pca_first_prop_warn = 0.10
      ),
      standard = list(
        n_obs_min = 100,
        n_person_min = 30,
        low_cat_min = 10,
        min_facet_levels = 3,
        misfit_ratio_warn = 0.10,
        missing_fit_ratio_warn = 0.20,
        zstd2_ratio_warn = 0.10,
        zstd3_ratio_warn = 0.05,
        expected_var_min = 0.20,
        pca_first_eigen_warn = 2.0,
        pca_first_prop_warn = 0.10
      ),
      lenient = list(
        n_obs_min = 60,
        n_person_min = 20,
        low_cat_min = 5,
        min_facet_levels = 2,
        misfit_ratio_warn = 0.15,
        missing_fit_ratio_warn = 0.30,
        zstd2_ratio_warn = 0.15,
        zstd3_ratio_warn = 0.08,
        expected_var_min = 0.10,
        pca_first_eigen_warn = 3.0,
        pca_first_prop_warn = 0.20
      )
    ),
    pca_reference_bands = list(
      eigenvalue = c(
        critical_minimum = 1.4,
        caution = 1.5,
        common = 2.0,
        strong = 3.0
      ),
      proportion = c(
        minor = 0.05,
        caution = 0.10,
        strong = 0.20
      )
    )
  )
}

resolve_warning_thresholds <- function(thresholds = NULL, threshold_profile = "standard") {
  refs <- warning_threshold_profiles()
  profile_name <- tolower(as.character(threshold_profile %||% "standard"))
  if (!profile_name %in% names(refs$profiles)) profile_name <- "standard"
  active <- refs$profiles[[profile_name]]
  if (!is.null(thresholds) && length(thresholds) > 0) {
    active[names(thresholds)] <- thresholds
  }
  list(
    profile_name = profile_name,
    thresholds = active,
    pca_reference_bands = refs$pca_reference_bands
  )
}

build_pca_reference_text <- function(reference_bands) {
  eigen <- reference_bands$eigenvalue
  prop <- reference_bands$proportion
  paste0(
    "Heuristic reference bands (uncalibrated): EV >= ",
    paste(vapply(unname(eigen), fmt_num, character(1), decimals = 1), collapse = ", "),
    "; variance >= ", paste(vapply(100 * unname(prop), fmt_num, character(1), decimals = 0), collapse = "%, "),
    "%. These values are descriptive references, not critical values or evidence for a number of dimensions."
  )
}

build_pca_check_text <- function(eigenvalue, proportion, reference_bands) {
  eigen <- to_float(eigenvalue)
  prop <- to_float(proportion)
  if (!is.finite(eigen) && !is.finite(prop)) {
    return("Current PC1 threshold checks are unavailable.")
  }

  e <- reference_bands$eigenvalue
  p <- reference_bands$proportion
  checks <- c(
    paste0("EV>=", fmt_num(e[["caution"]], 1), ":", ifelse(is.finite(eigen) && eigen >= e[["caution"]], "Y", "N")),
    paste0("EV>=", fmt_num(e[["common"]], 1), ":", ifelse(is.finite(eigen) && eigen >= e[["common"]], "Y", "N")),
    paste0("EV>=", fmt_num(e[["strong"]], 1), ":", ifelse(is.finite(eigen) && eigen >= e[["strong"]], "Y", "N")),
    paste0("Var>=", fmt_num(100 * p[["caution"]], 0), "%:", ifelse(is.finite(prop) && prop >= p[["caution"]], "Y", "N")),
    paste0("Var>=", fmt_num(100 * p[["strong"]], 0), "%:", ifelse(is.finite(prop) && prop >= p[["strong"]], "Y", "N"))
  )
  paste0("Current exploratory PC1 checks: ", paste(checks, collapse = ", "), ".")
}

extract_overall_pca_first <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$overall_table) || nrow(pca_obj$overall_table) == 0) return(NULL)
  tbl <- pca_obj$overall_table
  idx <- which.min(suppressWarnings(as.numeric(tbl$Component)))
  if (!length(idx) || !is.finite(idx)) return(NULL)
  tbl[idx[1], , drop = FALSE]
}

extract_overall_pca_second <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$overall_table) || nrow(pca_obj$overall_table) == 0) return(NULL)
  tbl <- pca_obj$overall_table
  comp <- suppressWarnings(as.numeric(tbl$Component))
  idx <- which(comp == 2)
  if (length(idx) == 0) return(NULL)
  tbl[idx[1], , drop = FALSE]
}

extract_facet_pca_first <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$by_facet_table) || nrow(pca_obj$by_facet_table) == 0) {
    return(data.frame())
  }

  tbl <- pca_obj$by_facet_table
  if (!"Facet" %in% names(tbl) || !"Component" %in% names(tbl)) return(data.frame())

  split_tbl <- split(tbl, as.character(tbl$Facet))
  out <- lapply(split_tbl, function(df) {
    comp <- suppressWarnings(as.numeric(df$Component))
    idx <- which.min(comp)
    if (!length(idx) || !is.finite(idx)) return(NULL)
    df[idx[1], , drop = FALSE]
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.frame())
  out <- dplyr::bind_rows(out)
  out[order(suppressWarnings(as.numeric(out$Eigenvalue)), decreasing = TRUE), , drop = FALSE]
}

extract_overall_pca_error <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$errors) || is.null(pca_obj$errors$overall)) return("")
  err <- as.character(pca_obj$errors$overall[1] %||% "")
  if (is.na(err)) "" else err
}

extract_facet_pca_errors <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$errors) || is.null(pca_obj$errors$by_facet)) {
    return(data.frame(Facet = character(0), Error = character(0), stringsAsFactors = FALSE))
  }
  as.data.frame(pca_obj$errors$by_facet, stringsAsFactors = FALSE)
}

collapse_apa_paragraph <- function(sentences, width = 92L,
                                    output_mode = c("wrapped", "reflow")) {
  lines <- trimws(as.character(sentences %||% character(0)))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return("")

  output_mode <- match.arg(output_mode)
  width <- suppressWarnings(as.integer(width))
  if (!is.finite(width) || width < 40L) width <- 92L

  txt <- paste(lines, collapse = " ")
  if (identical(output_mode, "reflow")) {
    # Reflow mode returns a single line without hard breaks, which is
    # what Word / RMarkdown / Quarto want when pasting into a running
    # paragraph. The default "wrapped" mode retains the classic 92-
    # char console readability.
    return(txt)
  }
  wrapped <- strwrap(txt, width = width)
  paste(wrapped, collapse = "\n")
}

mfrm_centering_scope <- function(config) {
  config$centering_scope %||%
    if (length(config$facet_specs)) "remainder" else "not_recorded"
}

summarize_anchor_constraints <- function(config) {
  anchor_tbl <- as.data.frame(config$anchor_summary %||% data.frame(), stringsAsFactors = FALSE)
  noncenter <- as.character(config$noncenter_facet %||% "Person")
  dummy_facets <- as.character(config$dummy_facets %||% character(0))
  centering <- switch(mfrm_centering_scope(config),
    facet = " Centered facets sum to zero across all levels, including anchors and group totals.",
    remainder = " Legacy centering sets the ungrouped, unanchored remainder to sum to zero.",
    anchor_origin = " Direct or group anchors set their facet's origin without zero-sum centering; other facets are centered except the declared unanchored noncenter facet.",
    " The centering convention was not recorded.")
  if (nrow(anchor_tbl) == 0) {
    dummy_txt <- if (length(dummy_facets) > 0) paste(dummy_facets, collapse = ", ") else "none"
    return(
      paste0(
        "Constraint settings: noncenter facet = ", noncenter,
        "; anchored levels = 0; group anchors = 0; dummy facets = ", dummy_txt, ".", centering
      )
    )
  }

  anchored <- suppressWarnings(as.numeric(anchor_tbl$AnchoredLevels))
  grouped <- suppressWarnings(as.numeric(anchor_tbl$GroupAnchors))
  anchored_total <- sum(ifelse(is.finite(anchored), anchored, 0), na.rm = TRUE)
  grouped_total <- sum(ifelse(is.finite(grouped), grouped, 0), na.rm = TRUE)

  anchor_facets <- if ("Facet" %in% names(anchor_tbl) && length(anchored) == nrow(anchor_tbl)) {
    as.character(anchor_tbl$Facet[is.finite(anchored) & anchored > 0])
  } else {
    character(0)
  }
  group_facets <- if ("Facet" %in% names(anchor_tbl) && length(grouped) == nrow(anchor_tbl)) {
    as.character(anchor_tbl$Facet[is.finite(grouped) & grouped > 0])
  } else {
    character(0)
  }
  dummy_from_tbl <- if ("Facet" %in% names(anchor_tbl) && "DummyFacet" %in% names(anchor_tbl)) {
    as.character(anchor_tbl$Facet[isTRUE(anchor_tbl$DummyFacet) | (is.logical(anchor_tbl$DummyFacet) & anchor_tbl$DummyFacet)])
  } else {
    character(0)
  }
  dummy_all <- unique(c(dummy_facets, dummy_from_tbl))

  anchor_txt <- if (length(anchor_facets) > 0) paste(anchor_facets, collapse = ", ") else "none"
  group_txt <- if (length(group_facets) > 0) paste(group_facets, collapse = ", ") else "none"
  dummy_txt <- if (length(dummy_all) > 0) paste(dummy_all, collapse = ", ") else "none"

  paste0(
    "Constraint settings: noncenter facet = ", noncenter,
    "; anchored levels = ", fmt_count(anchored_total), " (facets: ", anchor_txt, ")",
    "; group anchors = ", fmt_count(grouped_total), " (facets: ", group_txt, ")",
    "; dummy facets = ", dummy_txt, ".", centering
  )
}

summarize_convergence_metrics <- function(summary_row) {
  if (is.null(summary_row) || nrow(summary_row) == 0) {
    return("Optimization summary was not available.")
  }

  convergence <- mfrm_convergence_state(summary_row)
  converged <- convergence$code_converged
  iter <- if ("Iterations" %in% names(summary_row)) to_float(summary_row$Iterations[1]) else NA_real_
  fn_eval <- if ("FunctionEvaluations" %in% names(summary_row)) to_float(summary_row$FunctionEvaluations[1]) else iter
  gr_eval <- if ("GradientEvaluations" %in% names(summary_row)) to_float(summary_row$GradientEvaluations[1]) else NA_real_
  loglik <- if ("LogLik" %in% names(summary_row)) to_float(summary_row$LogLik[1]) else NA_real_
  aic <- if ("AIC" %in% names(summary_row)) to_float(summary_row$AIC[1]) else NA_real_
  bic <- if ("BIC" %in% names(summary_row)) to_float(summary_row$BIC[1]) else NA_real_
  sabic <- if ("SABIC" %in% names(summary_row)) to_float(summary_row$SABIC[1]) else NA_real_
  contract_current <- "ICContractVersion" %in% names(summary_row) &&
    identical(
      as.character(summary_row$ICContractVersion[1]),
      mfrm_ic_contract_version()
    )
  ic_eligible <- contract_current && "ICEligible" %in% names(summary_row) &&
    isTRUE(as.logical(summary_row$ICEligible[1]))
  ic_selectable <- ic_eligible && "ICSelectable" %in% names(summary_row) &&
    isTRUE(as.logical(summary_row$ICSelectable[1]))
  integration_q <- if ("ICQuadraturePoints" %in% names(summary_row)) {
    suppressWarnings(as.integer(summary_row$ICQuadraturePoints[1]))
  } else {
    NA_integer_
  }
  status <- if ("ConvergenceStatus" %in% names(summary_row)) as.character(summary_row$ConvergenceStatus[1]) else NA_character_
  detail <- if ("ConvergenceDetail" %in% names(summary_row)) as.character(summary_row$ConvergenceDetail[1]) else NA_character_
  grad_sup <- if ("TerminalGradientSupNorm" %in% names(summary_row)) to_float(summary_row$TerminalGradientSupNorm[1]) else NA_real_
  grad_tol <- if ("GradientReviewTolerance" %in% names(summary_row)) to_float(summary_row$GradientReviewTolerance[1]) else NA_real_

  conv_txt <- if (isTRUE(converged) && identical(convergence$severity, "pass")) {
    "met the numerical convergence checks"
  } else if (isTRUE(converged) && identical(status, "converged_gradient_review")) {
    "finished, but the terminal gradient requires review"
  } else if (identical(status, "reviewable_warning")) {
    "ended with a reviewable optimizer warning"
  } else if (identical(converged, FALSE)) {
    "did not meet the package convergence checks"
  } else {
    "had unknown convergence status"
  }
  iter_txt <- if (is.finite(fn_eval)) fmt_count(fn_eval) else if (is.finite(iter)) fmt_count(iter) else "NA"
  gr_txt <- if (is.finite(gr_eval)) fmt_count(gr_eval) else "NA"
  ll_txt <- if (is.finite(loglik)) fmt_num(loglik, 3) else "NA"
  aic_txt <- if (is.finite(aic)) fmt_num(aic, 3) else "NA"
  bic_txt <- if (is.finite(bic)) fmt_num(bic, 3) else "NA"
  sabic_txt <- if (is.finite(sabic)) fmt_num(sabic, 3) else "NA"
  criterion_txt <- if (ic_eligible) {
    paste0(
      ", canonical MML AIC = ", aic_txt,
      ", Person-BIC = ", bic_txt,
      ", Sclove SABIC = ", sabic_txt
    )
  } else {
    ""
  }
  out <- paste0(
    "Optimization ", conv_txt, " after ", iter_txt,
    " function evaluations",
    if (is.finite(gr_eval)) paste0(" and ", gr_txt, " gradient evaluations") else "",
    " (LogLik = ", ll_txt,
    criterion_txt, ")."
  )
  if (!ic_eligible) {
    out <- paste0(out, " MML model-comparison criteria are unavailable for this fit; numerical completion alone does not establish comparability.")
  }
  if (ic_eligible && !ic_selectable) {
    out <- paste0(
      out,
      " These criteria are screening/review values only at q=",
      if (is.finite(integration_q)) integration_q else "unknown",
      ". Automatic model comparisons remain unavailable; review the model and integration assumptions."
    )
  }
  if (identical(as.character(summary_row$Method[1]), "MML")) {
    out <- paste0(out, " MML integration used ",
      as.character(summary_row$MMLIntegration[1] %||% "fixed"),
      " Gauss-Hermite quadrature (q=",
      if (is.finite(integration_q)) integration_q else "unknown", ").")
  }
  persons <- if ("Persons" %in% names(summary_row)) {
    suppressWarnings(as.integer(summary_row$Persons[1]))
  } else {
    NA_integer_
  }
  if (ic_eligible && is.finite(persons) && persons <= 22L) {
    out <- paste0(
      out,
      " SABIC was retained for sensitivity only; automatic selection is ",
      "disabled at 22 or fewer Persons."
    )
  } else if (!ic_eligible) {
    legacy_aic <- if ("LegacyAIC" %in% names(summary_row)) {
      to_float(summary_row$LegacyAIC[1])
    } else if (!contract_current && "AIC" %in% names(summary_row)) {
      to_float(summary_row$AIC[1])
    } else {
      NA_real_
    }
    legacy_bic <- if ("LegacyBIC" %in% names(summary_row)) {
      to_float(summary_row$LegacyBIC[1])
    } else if (!contract_current && "BIC" %in% names(summary_row)) {
      to_float(summary_row$BIC[1])
    } else {
      NA_real_
    }
    if (is.finite(legacy_aic) || is.finite(legacy_bic)) {
      out <- paste0(
        out,
        " Legacy descriptive AIC = ",
        if (is.finite(legacy_aic)) fmt_num(legacy_aic, 3) else "NA",
        "; legacy descriptive BIC = ",
        if (is.finite(legacy_bic)) fmt_num(legacy_bic, 3) else "NA",
        "; neither enters the common MML ranking panel."
      )
    }
  }
  if (is.finite(grad_sup)) {
    out <- paste0(
      out,
      " Terminal gradient sup-norm = ", fmt_num(grad_sup, 4),
      if (is.finite(grad_tol)) paste0(" (review threshold = ", fmt_num(grad_tol, 4), ")") else "",
      "."
    )
  }
  if (!is.na(detail) && nzchar(detail) && !grepl("^Optimizer returned convergence code", detail)) {
    out <- paste0(out, " ", detail)
  }
  out
}

category_usage_note <- function(usage, low_count = 10) {
  if (!usage$Categories) return("Category counts were not available.")
  if (usage$UnavailableCounts > 0L) {
    return(sprintf("Category counts were available for %d of %d categories; usage totals remain unavailable.",
                   usage$AvailableCounts, usage$Categories))
  }
  sprintf("Category counts were available for all %d categories: %d unused and %d below %g. Counts alone do not establish category adequacy.",
          usage$Categories, usage$UnusedCategories, usage$LowCountCategories, low_count)
}

threshold_order_note <- function(coverage) {
  if (isTRUE(coverage$NotApplicable)) {
    return("Threshold ordering is not applicable: each supplied ladder has a single threshold.")
  }
  if (coverage$Available == 0L) {
    return("Threshold-order comparisons were unavailable; no ordering conclusion was drawn.")
  }
  if (!is.finite(coverage$Comparisons)) {
    return(sprintf("Adjacent threshold comparisons: %d decreasing among %d available; the total comparison count was unavailable.",
                   coverage$Decreasing, coverage$Available))
  }
  sprintf("Adjacent threshold comparisons: %d decreasing among %d available; %d of %d comparisons unavailable.",
          coverage$Decreasing, coverage$Available, coverage$Unavailable, coverage$Comparisons)
}

summarize_step_estimates <- function(step_tbl, expected_steps = NULL, expected_facets = NULL) {
  step_order <- calc_step_order(step_tbl, expected_steps = expected_steps, expected_facets = expected_facets)
  if (!nrow(step_order)) return("Step/threshold estimates and ordering comparisons were not available.")
  coverage <- summarize_threshold_order(step_order)
  est <- suppressWarnings(as.numeric(step_order$Estimate))
  est <- est[is.finite(est)]
  range_text <- if (length(est)) sprintf("Available estimates range from %s to %s logits.",
    fmt_num(min(est)), fmt_num(max(est))) else "Step/threshold estimates were not available."
  bad <- step_order[step_order$Ordered %in% FALSE, , drop = FALSE]
  detail <- if (nrow(bad)) paste0("Decreasing pairs end at ",
    paste(head(paste0(bad$StepFacet, ":", bad$Step), 3L), collapse = ", "), ".") else character(0)
  paste(c(range_text, threshold_order_note(coverage), detail), collapse = " ")
}

summarize_top_misfit_levels <- function(fit_tbl, top_n = 3L) {
  if (is.null(fit_tbl) || nrow(fit_tbl) == 0) {
    return("Top misfit levels were not available.")
  }

  tbl <- as.data.frame(fit_tbl, stringsAsFactors = FALSE)
  screening <- fit_screening_summary(tbl)
  absz <- screening$max_abs_z
  metric_label <- "|ZSTD|"
  if (!any(is.finite(absz))) {
    absz <- screening$max_mnsq_deviation
    metric_label <- "|MnSq - 1|"
  }

  tbl$AbsMetric <- absz
  tbl <- tbl[is.finite(tbl$AbsMetric), , drop = FALSE]
  if (nrow(tbl) == 0) {
    return("Top misfit levels were not available.")
  }

  tbl <- tbl[order(tbl$AbsMetric, decreasing = TRUE), , drop = FALSE]
  show <- utils::head(tbl, n = max(1L, as.integer(top_n)))
  labels <- vapply(seq_len(nrow(show)), function(i) {
    facet <- if ("Facet" %in% names(show)) as.character(show$Facet[i]) else "Facet"
    level <- if ("Level" %in% names(show)) as.character(show$Level[i]) else as.character(i)
    paste0(facet, ":", level,
           " (", metric_label, " = ", fmt_num(show$AbsMetric[i]), ")")
  }, character(1))

  paste0("Largest misfit signals among ", nrow(tbl), " elements with complete paired statistics: ",
         paste(labels, collapse = "; "), ".")
}

summarize_bias_counts <- function(bias_results) {
  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return("Bias analysis was not estimated in this run.")
  }

  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  pvals <- suppressWarnings(as.numeric(tbl$`Prob.`))
  sig_n <- sum(is.finite(pvals) & pvals < 0.05, na.rm = TRUE)
  total_n <- nrow(tbl)

  eff <- suppressWarnings(as.numeric(tbl$`Bias Size`))
  large_n <- sum(is.finite(eff) & abs(eff) >= 0.5, na.rm = TRUE)

  paste0(
    "Bias screening evaluated ", fmt_count(total_n),
    " interaction cell(s); ", fmt_count(sig_n),
    " met a screening tail area below .05 and ", fmt_count(large_n),
    " had |Bias Size| >= 0.50 logits."
  )
}

summarize_population_model_for_apa <- function(res) {
  population <- res$population %||% list()
  active <- isTRUE(population$active)
  empty <- list(
    active = FALSE,
    identification_only = FALSE,
    formula_label = "",
    coefficient_count = 0L,
    residual_variance = NA_real_,
    design_columns = character(0),
    coding = data.frame(),
    coding_sentence = "",
    method_sentence = "",
    result_sentence = "",
    omission_sentence = "",
    caution_sentence = "",
    conquest_sentence = ""
  )
  if (!active) {
    return(empty)
  }
  gpcm_identification_only <- identical(
    as.character(population$source %||% ""),
    "gpcm_mml_default_identification"
  )

  formula_label <- if (!is.null(population$formula)) {
    paste(deparse(population$formula), collapse = " ")
  } else {
    as.character(res$config$population_formula %||% "<unspecified>")
  }
  coefficients <- suppressWarnings(as.numeric(population$coefficients %||% numeric(0)))
  design_columns <- as.character(population$design_columns %||% names(population$coefficients) %||% character(0))
  if (length(design_columns) == 0L && length(coefficients) > 0L) {
    design_columns <- paste0("b", seq_along(coefficients))
  }
  coefficient_count <- length(coefficients)
  finite_coefficient_count <- sum(is.finite(coefficients))
  sigma2 <- suppressWarnings(as.numeric(population$sigma2 %||% NA_real_))
  policy <- as.character(population$policy %||% NA_character_)
  omitted_persons <- length(population$omitted_persons %||% character(0))
  omitted_rows <- suppressWarnings(as.integer(population$response_rows_omitted %||% 0L))
  if (!is.finite(omitted_rows)) {
    omitted_rows <- 0L
  }

  coding <- population_coding_summary_table(population)
  coding_sentence <- if (nrow(coding) > 0L) {
    coding_labels <- vapply(seq_len(nrow(coding)), function(i) {
      row <- coding[i, , drop = FALSE]
      encoded <- as.character(row$EncodedColumns[1] %||% "")
      if (!nzchar(encoded)) {
        encoded <- "none"
      }
      paste0(
        as.character(row$Variable[1] %||% "covariate"),
        " levels [", as.character(row$Levels[1] %||% ""), "]",
        ", contrast ", as.character(row$Contrast[1] %||% ""),
        ", encoded columns [", encoded, "]"
      )
    }, character(1))
    paste0(
    "Categorical population covariate coding was recorded: ",
      paste(coding_labels, collapse = "; "),
      "."
    )
  } else if (length(design_columns) > 0L) {
    "No categorical covariate coding was recorded for the population model; covariates appear intercept-only, numeric, or logical."
  } else {
    "Population-model covariate coding details were not recorded."
  }

  result_sentence <- paste0(
    if (gpcm_identification_only) {
      "The intercept-only GPCM population-scale model returned "
    } else {
      "The latent-regression population model returned "
    },
    fmt_count(finite_coefficient_count), " finite coefficient(s) out of ",
    fmt_count(coefficient_count),
    if (length(design_columns) > 0L) {
      paste0(" design column(s) [", paste(design_columns, collapse = ", "), "]")
    } else {
      " design column(s)"
    },
    if (is.finite(sigma2)) {
      paste0(", with residual variance = ", fmt_num(sigma2, 3))
    } else {
      ", with residual variance unavailable"
    },
    "."
  )
  omission_sentence <- paste0(
    "The population-model covariate policy was ", policy,
    "; omitted persons = ", fmt_count(omitted_persons),
    " and omitted response rows = ", fmt_count(omitted_rows), "."
  )
  method_sentence <- if (gpcm_identification_only) {
    paste0(
      "An intercept-only conditional-normal population model was activated via ",
      formula_label,
      " to estimate the GPCM MML location and common discrimination scale; it is not a substantive covariate regression."
    )
  } else {
    paste0(
      "A conditional-normal latent-regression population model was included via ",
      formula_label,
      "; coefficients were estimated jointly in the MML model, not as post hoc regression on EAP or MLE person scores."
    )
  }
  caution_sentence <- if (gpcm_identification_only) {
    "Interpret the population intercept and variance as GPCM identification coordinates; fixed-latent-SD discrimination equals population SD times the relative slope."
  } else {
    "Latent-regression coefficients should be interpreted as conditional-normal population-model parameters, not as post hoc regression on estimated person scores."
  }
  conquest_sentence <- if (gpcm_identification_only) {
    "The exact ConQuest scoresfree map is item-only; multifacet equivalence and external acceptance remain unestablished."
  } else {
    "ConQuest overlap is limited to the documented latent-regression MML comparison scope."
  }

  list(
    active = TRUE,
    identification_only = gpcm_identification_only,
    formula_label = formula_label,
    coefficient_count = coefficient_count,
    residual_variance = sigma2,
    design_columns = design_columns,
    coding = coding,
    coding_sentence = coding_sentence,
    method_sentence = method_sentence,
    result_sentence = result_sentence,
    omission_sentence = omission_sentence,
    caution_sentence = caution_sentence,
    conquest_sentence = conquest_sentence
  )
}

build_apa_table_figure_key_order <- function(include_population = FALSE) {
  keys <- c(
    "table1", "table2", "table3", "table4",
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds", "category_curves",
    "observed_expected", "fit_diagnostics", "fit_zstd_distribution", "misfit_levels",
    "strict_marginal_fit", "strict_pairwise_local_dependence",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  if (isTRUE(include_population)) {
    keys <- append(keys, "population_model", after = match("table4", keys))
  }
  keys
}

extract_apa_note_body <- function(note_text) {
  lines <- strsplit(as.character(note_text %||% ""), "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return("")
  note_idx <- grep("^Note\\.", lines)
  if (length(note_idx) == 0) return(paste(lines, collapse = " "))
  body <- sub("^Note\\.\\s*", "", lines[note_idx[1]])
  tail <- if (note_idx[1] < length(lines)) lines[(note_idx[1] + 1L):length(lines)] else character(0)
  tail <- tail[nzchar(tail)]
  paste(c(body, tail), collapse = " ")
}

build_apa_note_map_from_contract <- function(contract) {
  meta <- contract$metadata
  precision <- contract$precision
  summaries <- contract$summaries
  availability <- contract$availability

  positive_facets <- as.character(meta$positive_facets %||% character(0))
  orientation_clause <- if (length(positive_facets) > 0) {
    paste0(
      "higher person values indicate higher ability; for non-person facets, higher values indicate greater severity/difficulty under the default negative orientation, except ",
      paste(positive_facets, collapse = ", "),
      ", which were fit with positive orientation (higher values raise expected scores)."
    )
  } else {
    "higher person values indicate higher ability, and higher non-person facet values indicate greater severity/difficulty (all non-person facets used the default negative orientation in this fit)."
  }

  note_map <- list()
  note_map$table1 <- paste0(
    "Table 1. Facet summary\n",
    "Note. Measures are reported in logits; ", orientation_clause, " ",
    precision$se_note,
    summaries$ci_table_note %||% "",
    "Model = ", meta$model, "; estimation = ", meta$method,
    "; N = ", fmt_count(meta$n_obs), " observations from ", fmt_count(meta$n_person),
    " persons on a ", fmt_count(meta$n_cat), "-category scale (",
    fmt_count(meta$rating_min), "-", fmt_count(meta$rating_max), ")."
  )

  note_map$table2 <- paste0(
    "Table 2. Rating scale diagnostics\n",
    "Note. ", category_usage_note(summaries$category_usage %||% summarize_category_usage(NULL)), " ",
    threshold_order_note(summaries$threshold_coverage %||% summarize_threshold_order(NULL))
  )

  fit_sentence <- if (is.finite(summaries$overall_fit_infit) && is.finite(summaries$overall_fit_outfit)) {
    paste0(
      "Overall fit: infit MnSq = ", fmt_num(summaries$overall_fit_infit),
      ", outfit MnSq = ", fmt_num(summaries$overall_fit_outfit), "."
    )
  } else {
    "Overall fit indices are reported as mean infit and outfit MnSq."
  }

  note_map$table3 <- paste0(
    "Table 3. Fit and reliability summary\n",
    if (precision$supports_formal_inference) {
      "Note. Separation and reliability are based on observed variance, measurement error, and adjusted true variance. "
    } else if (identical(precision$tier, "hybrid")) {
      "Note. Separation and reliability combine model-based and fallback precision summaries and should be interpreted cautiously. "
    } else {
      "Note. Separation and reliability are exploratory summaries based on observed variance, measurement error, and adjusted true variance. "
    },
    fit_sentence
  )

  if (availability$has_rater_reliability) {
    note_map$table3 <- paste0(
      note_map$table3,
      " Rater facet (", summaries$rater_facet, ") ",
      precision$reliability_label,
      " = ", fmt_num(summaries$rater_reliability$Reliability),
      ", separation = ", fmt_num(summaries$rater_reliability$Separation), "."
    )
  }
  if (nzchar(summaries$interrater_sentence %||% "")) {
    note_map$table3 <- paste0(note_map$table3, " ", summaries$interrater_sentence)
  }

  if (availability$has_bias) {
    note_map$table4 <- paste0(
      "Table 4. Bias/interaction effects\n",
      "Note. Bias contrasts are in logits and represent observed minus expected scores with main effects held fixed. ",
      "The reported t and probability columns are screening metrics based on conditional plug-in information, not formal hypothesis tests. ",
      "Rows with screening tail area below .05 = ", fmt_count(summaries$bias_sig_n), ". ",
      mfrm_bias_uncertainty_note()
    )
  } else {
    note_map$table4 <- paste0(
      "Table 4. Bias/interaction effects\n",
      "Note. Bias contrasts are in logits and represent observed minus expected scores with main effects held fixed. ",
      "Any t or probability columns should be read as screening metrics rather than formal hypothesis tests."
    )
  }

  if (isTRUE(availability$has_population_model)) {
    population_summary <- summaries$population_model
    note_map$population_model <- paste0(
      if (isTRUE(population_summary$identification_only)) {
        "Table 5. GPCM population-scale identification model\n"
      } else {
        "Table 5. Latent-regression population model\n"
      },
      "Note. ", population_summary$caution_sentence, " ",
      "Formula = ", population_summary$formula_label,
      "; residual variance = ",
      if (is.finite(population_summary$residual_variance)) fmt_num(population_summary$residual_variance, 3) else "NA",
      ". ",
      population_summary$coding_sentence,
      " ",
      population_summary$omission_sentence
    )
  }

  note_map$wright_map <- paste0(
    "Wright map\nNote. Persons and facet elements are located on a common logit scale; ",
    orientation_clause
  )
  note_map$pathway_map <- "Pathway map\nNote. Curves show expected score across theta/logit levels from estimated thresholds. This expected-score display is distinct from the Bond-and-Fox-style measure-versus-fit pathway bubble chart, which is available via plot(x, type = \"bubble\")."
  note_map$facet_distribution <- "Facet estimate distribution\nNote. Distributions summarize severity/difficulty spread within each facet."
  note_map$step_thresholds <- "Step/threshold estimates\nNote. Adjacent threshold estimates are compared within each family; unavailable comparisons and decreasing pairs require separate review."
  note_map$category_curves <- "Category characteristic curves\nNote. Curves show category response probability across theta/logit levels; peaks describe the fitted model and do not by themselves establish category adequacy."
  note_map$observed_expected <- "Observed vs expected scores\nNote. Points summarize mean observed and expected scores by bin; deviations from the diagonal suggest local misfit."
  note_map$fit_diagnostics <- "Fit diagnostics (Infit vs Outfit)\nNote. Each point represents an element within a facet. Values near 1.0 indicate expected fit; values substantially above 1.0 suggest misfit."
  note_map$fit_zstd_distribution <- "Fit ZSTD distribution\nNote. ZSTD standardizes mean-square fit. Cutoffs are descriptive review references, without calibrated individual or multiple-element error rates; incomplete classifications remain visible."
  note_map$misfit_levels <- "Misfit levels\nNote. Narrative rankings use complete paired ZSTD values, or complete paired mean-square deviations when no paired ZSTD is available. The two metrics are not mixed in one ranking."
  if (isTRUE(availability$has_strict_marginal)) {
    marginal_tail <- paste0(
      " Overall RMSD = ", fmt_num(summaries$strict_marginal_overall_rmsd),
      ", overall max |standardized residual| = ", fmt_num(summaries$strict_marginal_max_abs_std_residual), "."
    )
    if (nzchar(summaries$strict_marginal_top_cell_label %||% "")) {
      marginal_tail <- paste0(
        marginal_tail,
        " Largest screening cell: ", summaries$strict_marginal_top_cell_label,
        " (standardized residual = ", fmt_num(summaries$strict_marginal_top_cell_std_residual),
        ", proportion difference = ", fmt_num(summaries$strict_marginal_top_cell_prop_diff), ")."
      )
    }
    note_map$strict_marginal_fit <- paste0(
      "Strict marginal fit screen\n",
      "Note. This latent-integrated screen summarizes category-level residuals for step/scale groups and facet levels. ",
      "Treat flagged cells as exploratory screening evidence rather than formal inferential tests.",
      marginal_tail, " ", paste(marginal_coverage_notes(summaries$strict_marginal_coverage), collapse = " ")
    )
  } else {
    note_map$strict_marginal_fit <- paste0(
      "Strict marginal fit screen\n",
      "Note. Strict marginal diagnostics were not available for this run. ",
      summaries$strict_marginal_reason %||% ""
    )
  }

  if (isTRUE(availability$has_strict_pairwise)) {
    pair_tail <- paste0(
      " Flagged level pairs = ", fmt_count(summaries$strict_pairwise_flagged_pairs), "."
    )
    if (nzchar(summaries$strict_pairwise_top_pair_label %||% "")) {
      pair_tail <- paste0(
        pair_tail,
        " Largest screening pair: ", summaries$strict_pairwise_top_pair_label,
        " (exact-agreement standardized residual = ", fmt_num(summaries$strict_pairwise_top_pair_exact_std_residual),
        ", adjacent-agreement standardized residual = ", fmt_num(summaries$strict_pairwise_top_pair_adjacent_std_residual), ")."
      )
    }
    note_map$strict_pairwise_local_dependence <- paste0(
      "Strict pairwise local dependence\n",
      "Note. This latent-integrated follow-up compares observed and expected exact/adjacent agreement within Person x remaining-facets contexts. ",
      "Treat flagged pairs as exploratory local-dependence evidence rather than standalone inferential tests.",
      pair_tail
    )
  } else {
    note_map$strict_pairwise_local_dependence <- paste0(
      "Strict pairwise local dependence\n",
      "Note. Strict pairwise local-dependence diagnostics were not available for this run. ",
      summaries$strict_pairwise_reason %||% ""
    )
  }

  if (availability$has_pca_overall) {
    overall_tail <- paste0(
      " PC1 eigenvalue = ", fmt_num(summaries$pca_overall_1$Eigenvalue),
      " (", fmt_num(100 * to_float(summaries$pca_overall_1$Proportion), 1), "% variance)."
    )
    if (!is.null(summaries$pca_overall_2)) {
      overall_tail <- paste0(
        overall_tail,
        " PC2 eigenvalue = ", fmt_num(summaries$pca_overall_2$Eigenvalue), "."
      )
    }
    note_map$residual_pca_overall <- paste0(
      "Residual PCA scree (overall)\n",
      "Note. Eigenvalues are from PCA of the person x facet-combination standardized residual correlation matrix.",
      overall_tail,
      " ",
      summaries$pca_reference_text
    )
  } else {
    note_map$residual_pca_overall <- paste0(
      "Residual PCA scree (overall)\n",
      "Note. Overall residual PCA was not available for this run. ",
      summaries$pca_reference_text
    )
  }

  if (availability$has_pca_by_facet) {
    top <- utils::head(summaries$pca_by_facet_first, 3)
    labels <- vapply(seq_len(nrow(top)), function(i) {
      paste0(
        top$Facet[i], ": ", fmt_num(top$Eigenvalue[i]), " (",
        fmt_num(100 * to_float(top$Proportion[i]), 1), "%)"
      )
    }, character(1))

    note_map$residual_pca_by_facet <- paste0(
      "Residual PCA by facet\n",
      "Note. Each facet is analyzed using a person x facet-level standardized residual matrix. ",
      "Largest PC1 signals: ", paste(labels, collapse = "; "), ". ",
      summaries$pca_reference_text
    )
  } else {
    note_map$residual_pca_by_facet <- paste0(
      "Residual PCA by facet\n",
      "Note. Facet-specific residual PCA was not available for this run. ",
      summaries$pca_reference_text
    )
  }

  note_map
}

build_apa_caption_map_from_contract <- function(contract) {
  assessment <- trimws(as.character(contract$context$assessment %||% ""))
  facet_pair <- contract$summaries$bias_facet_pair %||% ""
  assessment_phrase <- if (nzchar(assessment)) paste0(" for ", assessment) else ""

  caption_map <- list(
    table1 = paste0("Table 1\nFacet Summary (Measures, Precision, Fit, Reliability)", assessment_phrase),
    table2 = "Table 2\nRating Scale Diagnostics (Category Counts and Thresholds)",
    table3 = "Table 3\nFit and Reliability Summary",
    table4 = if (nzchar(facet_pair)) {
      paste0("Table 4\nBias/Interaction Effects for ", facet_pair)
    } else {
      "Table 4\nBias/Interaction Effects"
    },
    wright_map = paste0("Wright Map\nPerson and Facet Measures", assessment_phrase),
    pathway_map = "Pathway Map\nExpected Score by Theta",
    facet_distribution = "Facet Estimate Distribution",
    step_thresholds = "Step/Threshold Estimates",
    category_curves = "Category Characteristic Curves",
    observed_expected = "Observed vs. Expected Scores",
    fit_diagnostics = "Fit Diagnostics (Infit vs Outfit)",
    fit_zstd_distribution = "Fit ZSTD Distribution",
    misfit_levels = "Misfit Levels (Max |ZSTD|)",
    strict_marginal_fit = "Strict Marginal Fit Screen\nLatent-Integrated Category Residuals",
    strict_pairwise_local_dependence = "Strict Pairwise Local Dependence\nLatent-Integrated Agreement Residuals",
    residual_pca_overall = "Residual PCA Scree (Overall)",
    residual_pca_by_facet = "Residual PCA by Facet"
  )
  if (isTRUE(contract$availability$has_population_model)) {
    caption_map$population_model <- "Table 5\nLatent-Regression Population Model Coefficients and Coding"
  }
  caption_map
}

build_apa_section_entry <- function(parent, heading, sentences, width = 92L,
                                    output_mode = c("wrapped", "reflow")) {
  output_mode <- match.arg(output_mode)
  sentences <- trimws(as.character(sentences %||% character(0)))
  sentences <- sentences[nzchar(sentences)]
  list(
    Parent = as.character(parent),
    Heading = as.character(heading),
    Sentences = sentences,
    SentenceCount = length(sentences),
    Text = collapse_apa_paragraph(sentences, width = width,
                                  output_mode = output_mode),
    Available = length(sentences) > 0
  )
}

build_apa_section_map_from_contract <- function(contract) {
  width <- contract$metadata$line_width %||% 92L
  output_mode <- contract$metadata$output_mode %||% "wrapped"
  if (!output_mode %in% c("wrapped", "reflow")) output_mode <- "wrapped"
  entry <- function(parent, heading, sentences) {
    build_apa_section_entry(parent, heading, sentences,
                            width = width, output_mode = output_mode)
  }
  sections <- list(
    method_design = entry("Method", "Design and data", contract$method_design_sentences),
    method_estimation = entry("Method", "Estimation settings", contract$method_estimation_sentences),
    results_scale = entry("Results", "Scale functioning", contract$results_scale_sentences),
    results_measures = entry("Results", "Facet measures", contract$results_measure_sentences),
    results_population_model = entry("Results", "Latent-regression population model", contract$results_population_sentences),
    results_fit_precision = entry("Results", "Fit and precision", contract$results_fit_precision_sentences),
    results_residual_structure = entry("Results", "Residual structure", contract$results_residual_sentences),
    results_bias_screening = entry("Results", "Bias screening", contract$results_bias_sentences),
    results_cautions = entry("Results", "Reporting cautions", contract$caution_sentences)
  )
  sections
}

flatten_apa_section_map <- function(section_map, order = names(section_map)) {
  keys <- order[order %in% names(section_map)]
  if (length(keys) == 0) return(data.frame())
  do.call(
    rbind,
    lapply(keys, function(key) {
      entry <- section_map[[key]]
      data.frame(
        SectionId = as.character(key),
        Parent = as.character(entry$Parent %||% ""),
        Heading = as.character(entry$Heading %||% ""),
        Available = isTRUE(entry$Available),
        SentenceCount = as.integer(entry$SentenceCount %||% 0L),
        Text = as.character(entry$Text %||% ""),
        stringsAsFactors = FALSE
      )
    })
  )
}

build_apa_report_text_from_contract <- function(contract) {
  section_map <- contract$section_map %||% build_apa_section_map_from_contract(contract)
  method_keys <- c("method_design", "method_estimation")
  results_keys <- c(
    "results_scale", "results_measures", "results_population_model", "results_fit_precision",
    "results_residual_structure", "results_bias_screening", "results_cautions"
  )

  render_sections <- function(keys) {
    blocks <- vapply(keys[keys %in% names(section_map)], function(key) {
      entry <- section_map[[key]]
      if (!isTRUE(entry$Available) || !nzchar(entry$Text)) return("")
      paste0(entry$Heading, ".\n", entry$Text)
    }, character(1))
    blocks <- blocks[nzchar(blocks)]
    paste(blocks, collapse = "\n\n")
  }

  method_text <- paste0("Method.\n\n", render_sections(method_keys))
  results_text <- paste0("Results.\n\n", render_sections(results_keys))
  paste0(method_text, "\n\n", results_text)
}

mfrm_separation_interpretation <- function() {
  paste(
    "Separation reliability describes differences among facet locations relative to their estimation errors.",
    "High Person values indicate distinguishable abilities; high rater values indicate distinguishable severities, not inter-rater agreement.",
    "Low rater values can also reflect imprecision and do not establish agreement.",
    "The index does not test every rater pair or quantify the practical size of severity differences.",
    "Reliability, separation and Strata are approximate compatibility indices using marginal SEs; joint-information SEs do not validate this variance decomposition.",
    "They do not account for covariance among estimates. EAP-based separation is not posterior-variance EAP reliability.",
    "Model and fit-adjusted values are companion indices, not confidence limits. Strata describes hypothetical distinguishable strata, not discovered groups.",
    "Facet heterogeneity chi-square summaries use diagonal SE weights, not full-covariance joint tests."
  )
}

reporting_interrater_text <- function(diagnostics) {
  if (is.null(diagnostics)) return("")
  s <- as.data.frame(diagnostics$interrater$summary %||% data.frame())
  if (!nrow(s)) return(paste(
    "No inter-rater agreement summary was stored in the supplied diagnostics.",
    "Its absence does not establish agreement; check the rater facet and shared ratings."
  ))
  value <- function(name) to_float(s[[name]][1L])
  facet <- as.character(s$RaterFacet[1L] %||% "")
  if (is.na(facet) || !nzchar(facet)) facet <- "the stored rater facet"
  opportunities <- value("OpportunityCount")
  if (!is.finite(opportunities)) opportunities <- value("TotalPairs")
  text <- paste0("Inter-rater agreement for ", facet,
    " is reported separately from separation reliability.")
  text <- c(text, paste(
    "Diagnostics match scores within Person and all remaining non-rater facets;",
    "each compared pair contributes one opportunity per matched context."
  ))
  if (is.finite(opportunities)) {
    text <- c(text, paste0("The denominator is ", fmt_count(opportunities),
      " matched-context opportunities across rater pairs, not a count of independent people."))
    if (opportunities == 0) text <- c(text,
      "No matched opportunities were available; agreement is unavailable, not zero.")
  } else {
    text <- c(text, "The matched-opportunity denominator was not stored.")
  }
  if (all(is.finite(c(value("AvailablePairs"), value("Pairs"))))) {
    text <- c(text, paste0(fmt_count(value("AvailablePairs")), " of ",
      fmt_count(value("Pairs")), " candidate rater pairs have matched scores."))
  }
  metric <- function(name) if (is.finite(value(name))) fmt_num(value(name)) else "unavailable"
  text <- c(text, paste0("Observed exact agreement = ", metric("ExactAgreement"),
    "; observed adjacent agreement (absolute score difference at most 1, including exact matches) = ",
    metric("AdjacentAgreement"), "; model-expected exact agreement = ",
    metric("ExpectedExactAgreement"), "."))
  if (is.finite(value("RepeatedCells")) && value("RepeatedCells") > 0) {
    text <- c(text, paste0(fmt_count(value("RepeatedCells")),
      " rater/context cells contain repeated ratings; their scores are averaged using observation weights,",
      " and model-expected agreement is unavailable for comparisons involving these cells."))
  }
  if (is.finite(value("ExpectedUnavailableContexts")) && value("ExpectedUnavailableContexts") > 0) {
    text <- c(text, paste0("Valid fitted probabilities are unavailable for ",
      fmt_count(value("ExpectedUnavailableContexts")),
      " matched opportunities; the overall model-expected value is withheld rather than computed on a different denominator."))
  }
  if (any(grepl("primary_estimate_unavailable", s$ExpectedReason, fixed = TRUE))) {
    text <- c(text, "Some primary Person, facet, step or slope estimates are unavailable; finite optimizer proxies are not used for expected agreement.")
  }
  if (!"ExpectedReason" %in% names(s)) {
    text <- c(text, "The stored summary lacks the primary-estimate availability audit. Regenerate diagnose_mfrm() from the saved fit before interpreting its model-expected agreement.")
  }
  text <- c(text, paste(
    "The model-expected baseline uses fitted category probabilities with calibration fixed;",
    "it is not chance agreement or an inferential test.",
    "It assumes conditional independence at the saved Person point score (EAP for MML), rather than integrating over the Person posterior.",
    "These descriptive summaries do not classify rater quality or establish fairness."
  ))
  paste(text, collapse = " ")
}

build_apa_reporting_contract <- function(res, diagnostics, bias_results = NULL, context = list(), whexact = FALSE) {
  mfrm_results_validate_diagnostics_identity(res, diagnostics,
                                             helper = "APA reporting")
  mfrm_validate_bias_identity(res, bias_results, helper = "APA reporting")
  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  prep <- res$prep
  config <- res$config
  model <- toupper(as.character(config$model %||% "RSM"))
  fit_readiness_record <- mfrmr_get_readiness_record(res)
  fit_inference_ready <- isTRUE(
    fit_readiness_record$fit$InferenceReady[1]
  )

  n_obs <- to_float(prep$n_obs %||% summary$N)
  weight_sum <- to_float(prep$weighted_n %||% summary$N)
  n_person <- if (!is.null(summary)) to_float(summary$Persons) else nrow(res$facets$person)
  n_cat <- if (!is.null(summary)) to_float(summary$Categories) else to_float(config$n_cat)
  rating_min <- to_float(prep$rating_min)
  rating_max <- to_float(prep$rating_max)

  facet_names <- as.character(config$facet_names %||% character(0))
  facet_levels <- config$facet_levels %||% list()
  facet_counts <- if (length(facet_names) > 0) {
    vapply(facet_names, function(f) length(facet_levels[[f]] %||% character(0)), numeric(1))
  } else {
    numeric(0)
  }
  facets_text <- if (length(facet_counts) > 0) {
    paste(paste0(names(facet_counts), " (n = ", vapply(facet_counts, fmt_count, character(1)), ")"), collapse = ", ")
  } else {
    "no additional facets"
  }

  assessment <- trimws(as.character(context$assessment %||% ""))
  setting <- trimws(as.character(context$setting %||% ""))
  rater_training <- trimws(as.character(context$rater_training %||% ""))
  raters_per_response <- trimws(as.character(context$raters_per_response %||% ""))
  scale_desc <- trimws(as.character(context$scale_desc %||% ""))
  line_width <- suppressWarnings(as.integer(context$line_width %||% 92L))
  if (!is.finite(line_width) || line_width < 40L) line_width <- 92L

  precision_profile <- as.data.frame(diagnostics$precision_profile %||% data.frame(), stringsAsFactors = FALSE)
  precision_tier <- trimws(as.character(precision_profile$PrecisionTier[1] %||% NA_character_))
  if (!nzchar(precision_tier)) {
    precision_tier <- if (identical(config$method, "MML")) "hybrid" else "exploratory"
  }
  supports_formal_inference <- nrow(precision_profile) > 0 &&
    isTRUE(precision_profile$SupportsFormalInference[1]) &&
    fit_inference_ready
  precision_label <- if (supports_formal_inference) {
    "model-based"
  } else if (identical(precision_tier, "hybrid")) {
    "hybrid"
  } else {
    "exploratory"
  }
  recommended_use <- trimws(as.character(precision_profile$RecommendedUse[1] %||% ""))
  reliability_label <- if (supports_formal_inference) {
    "reliability"
  } else if (identical(precision_tier, "hybrid")) {
    "hybrid reliability summary"
  } else {
    "exploratory reliability summary"
  }
  se_note <- if (supports_formal_inference) {
    "Model S.E. = model-based standard error; Real S.E. = fit-adjusted standard error; MnSq = mean-square fit. "
  } else if (identical(precision_tier, "hybrid")) {
    "Model S.E. = primarily model-based standard error with fallback approximations when required; Real S.E. = fit-adjusted precision summary; MnSq = mean-square fit. "
  } else {
    "Model S.E. = exploratory standard error; Real S.E. = fit-adjusted exploratory standard error; MnSq = mean-square fit. "
  }
  precision_sentence <- if (supports_formal_inference) {
    "Model-based precision summaries were available for this run."
  } else if (identical(precision_tier, "hybrid")) {
    "Precision summaries combined model-based quantities with fallback approximations when needed."
  } else {
    "Precision summaries were exploratory in this run."
  }
  precision_caution <- if (supports_formal_inference) {
    ""
  } else if (identical(precision_tier, "hybrid")) {
    "Precision note: this run mixed model-based and fallback approximations, so confidence intervals and reliability summaries should be interpreted cautiously."
  } else {
    "Precision note: this run relies on exploratory precision summaries, so confidence intervals and reliability summaries should not be treated as formal inferential quantities."
  }

  method_design_sentences <- character(0)
  method_estimation_sentences <- character(0)
  method_sentences <- character(0)
  if (nzchar(assessment)) {
    assessment_sentence <- if (nzchar(setting)) {
      paste0("The analysis focused on ", assessment, " in ", setting, ".")
    } else {
      paste0("The analysis focused on ", assessment, ".")
    }
    method_design_sentences <- c(method_design_sentences, assessment_sentence)
    method_sentences <- c(method_sentences, assessment_sentence)
  }

  model_overview_label <- switch(
    model,
    RSM = "A many-facet rating-scale Rasch model",
    PCM = "A many-facet partial-credit Rasch model",
    GPCM = "A generalized partial-credit many-facet model",
    "A many-facet ordered-response model"
  )
  design_overview_sentence <- paste0(
    model_overview_label, " was fit to ", fmt_count(n_obs),
    " observations from ", fmt_count(n_person),
    " persons scored on a ", fmt_count(n_cat),
    "-category scale (", fmt_count(rating_min), "-", fmt_count(rating_max), ")."
  )
  if (any(abs(prep$data$Weight - 1) > 1e-12, na.rm = TRUE)) {
    design_overview_sentence <- paste0(
      design_overview_sentence, " The sum of observation weights was ",
      format(weight_sum, digits = 10, trim = TRUE),
      "; the fit summary's N field denotes this weight sum, not the number of rating rows."
    )
  }
  design_facets_sentence <- if (length(facet_names) > 0) {
    paste0("The design included facets for ", facets_text, ".")
  } else {
    "No additional facets beyond Person were modeled."
  }
  method_design_sentences <- c(method_design_sentences, design_overview_sentence, design_facets_sentence)
  method_sentences <- c(method_sentences, design_overview_sentence, design_facets_sentence)

  # Describe retained rating counts without equating them with independent
  # examinees, effective sample size, or qualified measurement precision.
  facet_flag <- as.character(summary$FacetSampleSizeFlag %||% NA_character_)
  facet_min_n <- suppressWarnings(as.integer(summary$FacetMinLevelN %||% NA_integer_))
  if (!is.na(facet_flag) && facet_flag %in% c("sparse", "marginal", "standard", "strong")) {
    adequacy_sentence <- paste0(
      "The smallest non-person facet-level count was ",
      if (is.na(facet_min_n)) "NA" else facet_min_n,
      " retained rating rows, placing it in the package's '", facet_flag,
      "' screening band. These counts are not counts of independent examinees ",
      "or effective sample sizes and do not establish precision or interval ",
      "coverage. Review targeting, connectedness and estimator-specific ",
      "uncertainty alongside `facet_small_sample_review()`."
    )
    method_design_sentences <- c(method_design_sentences, adequacy_sentence)
    method_sentences <- c(method_sentences, adequacy_sentence)
  }

  # Empirical-Bayes / Laplace shrinkage note. When the caller opted
  # into empirical-Bayes shrinkage, record the fact in Methods so reviewers see
  # both the fixed-effects fit and the post-hoc partial-pooling layer.
  shrinkage_mode <- as.character(config$facet_shrinkage %||% "none")
  if (!identical(shrinkage_mode, "none")) {
    validate_shrinkage_output(res)
    shrink_report <- res$shrinkage_report
    if (is.data.frame(shrink_report)) shrink_report <- shrink_report[shrink_report$Facet != "Person", , drop = FALSE]
    tau_mean <- if (!is.null(shrink_report) && nrow(shrink_report) > 0L) {
      mean(as.numeric(shrink_report$Tau2), na.rm = TRUE)
    } else NA_real_
    mean_shrink <- if (!is.null(shrink_report) && nrow(shrink_report) > 0L) {
      mean(as.numeric(shrink_report$MeanShrinkage), na.rm = TRUE)
    } else NA_real_
    shrinkage_sentence <- paste0(
      "Zero-centered normal-model shrinkage was applied after fitting. ",
      "Across non-Person facets, the mean selected prior variance was ",
      if (is.finite(tau_mean)) sprintf("tau^2 = %.3f", tau_mean) else "not identifiable",
      ", and the mean shrinkage factor was ",
      if (is.finite(mean_shrink)) sprintf("%.2f", mean_shrink) else "NA",
      ". Original fixed-effects estimates were retained. ",
      "The adjustment and plug-in SEs are descriptive; prior-variance uncertainty ",
      "and cross-level covariance are omitted, and zero SE after full pooling ",
      "does not establish perfect precision."
    )
    method_estimation_sentences <- c(method_estimation_sentences, shrinkage_sentence)
    method_sentences <- c(method_sentences, shrinkage_sentence)
  }

  # Extreme-score person warning. Under JMLE the
  # theta for such persons diverges; under MML the EAP is finite but
  # the information is small. Flag either tail in the Methods section
  # so reviewers understand why those persons may have been dropped or
  # reported at the truncation limit.
  ext_hi <- suppressWarnings(as.integer(summary$ExtremeHighN %||% NA_integer_))
  ext_lo <- suppressWarnings(as.integer(summary$ExtremeLowN %||% NA_integer_))
  if (isTRUE(is.finite(ext_hi) && is.finite(ext_lo) &&
             (ext_hi + ext_lo) > 0)) {
    extreme_sentence <- paste0(
      "Extreme-score persons were observed at the ceiling (n = ", ext_hi,
      ") and floor (n = ", ext_lo, "). Under JML the corresponding ",
      "theta estimates diverge toward infinity (Wright, 1998); under ",
      "MML they remain finite but carry little information. ",
      "Fit$facets$person$Extreme records the per-person flag."
    )
    method_design_sentences <- c(method_design_sentences, extreme_sentence)
    method_sentences <- c(method_sentences, extreme_sentence)
  }

  if (nzchar(scale_desc)) {
    scale_sentence <- paste0("The rating scale was described as ", scale_desc, ".")
    method_design_sentences <- c(method_design_sentences, scale_sentence)
    method_sentences <- c(method_sentences, scale_sentence)
  }
  if (nzchar(rater_training)) {
    training_sentence <- paste0("Raters received ", rater_training, ".")
    method_design_sentences <- c(method_design_sentences, training_sentence)
    method_sentences <- c(method_sentences, training_sentence)
  }
  if (nzchar(raters_per_response)) {
    rater_load_sentence <- paste0("Each response was scored by ", raters_per_response, " raters on average.")
    method_design_sentences <- c(method_design_sentences, rater_load_sentence)
    method_sentences <- c(method_sentences, rater_load_sentence)
  }

  population_summary <- summarize_population_model_for_apa(res)
  method <- config$method
  model_sentence <- paste0("The ", model, " specification was estimated using ", method, " with mfrmr.")
  if (identical(model, "PCM") && !is.null(config$step_facet) && nzchar(config$step_facet)) {
    model_sentence <- paste0(model_sentence, " The step structure varied by ", config$step_facet, ".")
  }
  estimation_basis_sentence <- if (identical(method, "MML")) {
    paste0(
      "Person measures are expected a posteriori (EAP) estimates under the ",
      "marginal person distribution, and residual-based fit statistics are ",
      "evaluated at these EAP measures rather than at joint maximum ",
      "likelihood (JMLE) estimates."
    )
  } else {
    ""
  }
  method_estimation_sentences <- c(method_estimation_sentences, model_sentence, precision_sentence)
  method_sentences <- c(method_sentences, model_sentence, precision_sentence)
  if (nzchar(estimation_basis_sentence)) {
    method_estimation_sentences <- c(method_estimation_sentences, estimation_basis_sentence)
    method_sentences <- c(method_sentences, estimation_basis_sentence)
  }
  if (isTRUE(population_summary$active)) {
    method_estimation_sentences <- c(
      method_estimation_sentences,
      population_summary$method_sentence,
      population_summary$coding_sentence
    )
    method_sentences <- c(
      method_sentences,
      population_summary$method_sentence,
      population_summary$coding_sentence
    )
  }
  if (nzchar(recommended_use)) {
    recommended_use_sentence <- paste0("Recommended use for this precision profile: ", recommended_use, ".")
    method_estimation_sentences <- c(method_estimation_sentences, recommended_use_sentence)
    method_sentences <- c(method_sentences, recommended_use_sentence)
  }

  if (!is.null(config$weight_col) && nzchar(config$weight_col)) {
    weight_sentence <- paste0(
      "Positive observation weights from `", as.character(config$weight_col[1]),
      "` multiplied row-level ordered-category log-likelihood contributions.",
      if (identical(method, "MML")) {
        paste0(
          " The Person distribution was integrated once per Person; the row weights ",
          "were not interpreted as replicated independent Person response patterns."
        )
      } else {
        ""
      }
    )
    method_estimation_sentences <- c(method_estimation_sentences, weight_sentence)
    method_sentences <- c(method_sentences, weight_sentence)
  }
  convergence_sentence <- summarize_convergence_metrics(summary)
  anchor_sentence <- summarize_anchor_constraints(config)
  method_estimation_sentences <- c(method_estimation_sentences, convergence_sentence, anchor_sentence)
  method_sentences <- c(method_sentences, convergence_sentence, anchor_sentence)

  cat_tbl <- if (!is.null(diagnostics)) calc_category_stats(diagnostics$obs, res = res, whexact = whexact) else tibble::tibble()
  step_order <- calc_step_order(res$steps, expected_steps = res$config$n_cat - 1L,
    expected_facets = res$config$facet_levels[[res$config$step_facet %||% ""]])
  category_usage <- summarize_category_usage(cat_tbl)
  threshold_coverage <- summarize_threshold_order(step_order)
  unused <- category_usage$UnusedCategories
  low_count <- category_usage$LowCountCategories
  threshold_text <- threshold_order_note(threshold_coverage)

  results_scale_sentences <- character(0)
  results_measure_sentences <- character(0)
  results_fit_precision_sentences <- character(0)
  results_population_sentences <- character(0)
  results_residual_sentences <- character(0)
  results_bias_sentences <- character(0)
  results_sentences <- character(0)
  category_sentence <- paste(category_usage_note(category_usage), threshold_text)

  results_scale_sentences <- c(results_scale_sentences, category_sentence)
  results_sentences <- c(results_sentences, category_sentence)
  step_summary_sentence <- summarize_step_estimates(res$steps, expected_steps = res$config$n_cat - 1L,
    expected_facets = res$config$facet_levels[[res$config$step_facet %||% ""]])
  results_scale_sentences <- c(results_scale_sentences, step_summary_sentence)
  results_sentences <- c(results_sentences, step_summary_sentence)

  person_stats <- describe_series(res$facets$person$Estimate)
  if (!is.null(person_stats)) {
    person_sentence <- paste0(
      "Person measures ranged from ", fmt_num(person_stats$min), " to ", fmt_num(person_stats$max),
      " logits (M = ", fmt_num(person_stats$mean), ", SD = ", fmt_num(person_stats$sd), ")."
    )
    results_measure_sentences <- c(results_measure_sentences, person_sentence)
    results_sentences <- c(results_sentences, person_sentence)
  }

  if (!is.null(res$facets$others) && nrow(res$facets$others) > 0) {
    for (facet in facet_names) {
      df_f <- res$facets$others |> dplyr::filter(Facet == facet)
      stats_f <- describe_series(df_f$Estimate)
      if (!is.null(stats_f)) {
        facet_sentence <- paste0(
          facet, " measures ranged from ", fmt_num(stats_f$min), " to ", fmt_num(stats_f$max),
          " logits (M = ", fmt_num(stats_f$mean), ", SD = ", fmt_num(stats_f$sd), ")."
        )
        results_measure_sentences <- c(results_measure_sentences, facet_sentence)
        results_sentences <- c(results_sentences, facet_sentence)
      }
    }
  }

  if (isTRUE(population_summary$active)) {
    results_population_sentences <- c(
      results_population_sentences,
      population_summary$result_sentence,
      population_summary$omission_sentence
    )
    results_sentences <- c(results_sentences, results_population_sentences)
  }

  band <- mfrm_misfit_thresholds()
  band_lower <- as.numeric(band["lower"])
  band_upper <- as.numeric(band["upper"])
  band_text <- sprintf("%.1f-%.1f", band_lower, band_upper)

  overall_fit <- if (!is.null(diagnostics$overall_fit) && nrow(diagnostics$overall_fit) > 0) diagnostics$overall_fit[1, , drop = FALSE] else NULL
  infit <- outfit <- NA_real_
  if (!is.null(overall_fit)) {
    infit <- to_float(overall_fit$Infit)
    outfit <- to_float(overall_fit$Outfit)
    overall_flag <- fit_screening_summary(overall_fit, band)$flags$MeanSquare[1]
    fit_label <- if (is.na(overall_flag)) "unavailable" else if (overall_flag) "outside" else "within"
    fit_sentence <- paste0(
      if (is.na(overall_flag)) "Overall mean-square screening was unavailable" else
        paste0("Overall mean-square fit was ", fit_label, " the ", band_text, " screening band"),
      " (infit MnSq = ", fmt_num(infit),
      ", outfit MnSq = ", fmt_num(outfit), "). ",
      "This band is the package's review convention; published mean-square ",
      "guidelines differ, and band position is screening evidence rather ",
      "than a model-validity decision."
    )
    results_fit_precision_sentences <- c(results_fit_precision_sentences, fit_sentence)
    results_sentences <- c(results_sentences, fit_sentence)
  }

  fit_tbl <- diagnostics$fit
  fit_screening <- fit_screening_summary(fit_tbl, band)$counts
  misfit_n <- fit_screening$Flagged[1]
  misfit_total <- fit_screening$Elements[1]
  top_misfit_sentence <- "Top misfit levels were not available."
  misfit_sentence <- fit_screening_note(fit_screening[1, ])
  results_fit_precision_sentences <- c(results_fit_precision_sentences, misfit_sentence)
  results_sentences <- c(results_sentences, misfit_sentence)
  if (!is.null(fit_tbl) && nrow(fit_tbl) > 0) {
    top_misfit_sentence <- summarize_top_misfit_levels(fit_tbl, top_n = 3L)
    results_fit_precision_sentences <- c(results_fit_precision_sentences, top_misfit_sentence)
    results_sentences <- c(results_sentences, top_misfit_sentence)
  }

  rel_tbl <- diagnostics$reliability
  rater_facet <- trimws(as.character(context$rater_facet %||% ""))
  interrater_summary <- as.data.frame(diagnostics$interrater$summary %||% data.frame(), stringsAsFactors = FALSE)
  if (!nzchar(rater_facet) && nrow(interrater_summary) > 0 && "RaterFacet" %in% names(interrater_summary)) {
    rater_facet <- trimws(as.character(interrater_summary$RaterFacet[1] %||% ""))
  }
  rater_rel <- NULL
  if (nzchar(rater_facet) && !is.null(rel_tbl) && nrow(rel_tbl) > 0) {
    match <- rel_tbl |> dplyr::filter(.data$Facet == rater_facet)
    if (nrow(match) > 0) rater_rel <- match[1, , drop = FALSE]
  }
  if (!is.null(rel_tbl) && nrow(rel_tbl) > 0) {
    rel_lines <- vapply(seq_len(nrow(rel_tbl)), function(i) {
      row <- rel_tbl[i, , drop = FALSE]
      paste0(row$Facet, " ", reliability_label, " = ", fmt_num(row$Reliability), " (separation = ", fmt_num(row$Separation), ").")
    }, character(1))
    reliability_basis_sentence <- mfrm_separation_interpretation()
    person_reliability_note <- if (identical(method, "MML") &&
                                   any(rel_tbl$Facet == "Person", na.rm = TRUE)) {
      paste0(
        "The Person row uses EAP measures with posterior SDs, which yields ",
        "a conservative summary that is not numerically comparable to ",
        "JMLE-based person reliability from FACETS."
      )
    } else {
      ""
    }
    reliability_sentence <- paste(
      Filter(nzchar, c(rel_lines, reliability_basis_sentence, person_reliability_note)),
      collapse = " "
    )
    results_fit_precision_sentences <- c(results_fit_precision_sentences, reliability_sentence)
    results_sentences <- c(results_sentences, reliability_sentence)
  }
  interrater_sentence <- reporting_interrater_text(diagnostics)
  if (nzchar(interrater_sentence)) {
    results_fit_precision_sentences <- c(results_fit_precision_sentences, interrater_sentence)
    results_sentences <- c(results_sentences, interrater_sentence)
  }

  measures_ci_tbl <- as.data.frame(diagnostics$measures %||% data.frame(), stringsAsFactors = FALSE)
  ci_available_n <- 0L
  ci_sentence <- ""
  ci_table_note <- ""
  if (nrow(measures_ci_tbl) > 0 && all(c("CI_Lower", "CI_Upper") %in% names(measures_ci_tbl))) {
    ci_lower <- suppressWarnings(as.numeric(measures_ci_tbl$CI_Lower))
    ci_upper <- suppressWarnings(as.numeric(measures_ci_tbl$CI_Upper))
    ci_available_n <- sum(is.finite(ci_lower) & is.finite(ci_upper))
    if (ci_available_n > 0) {
      ci_level <- suppressWarnings(as.numeric(measures_ci_tbl$CI_Level[1] %||% NA_real_))
      ci_level_text <- if (is.finite(ci_level)) paste0(format(round(100 * ci_level)), "%") else "95%"
      ci_method_text <- trimws(as.character(measures_ci_tbl$CI_Method[1] %||% ""))
      if (!nzchar(ci_method_text)) ci_method_text <- "normal approximation"
      ci_eligible_n <- if ("CIEligible" %in% names(measures_ci_tbl)) {
        sum(is.finite(ci_lower) & is.finite(ci_upper) &
              vapply(measures_ci_tbl$CIEligible, isTRUE, logical(1)))
      } else {
        NA_integer_
      }
      ci_sentence <- paste0(
        "Element-level ", ci_level_text, " approximate intervals (",
        ci_method_text, ") accompany ", fmt_count(ci_available_n), " of ",
        fmt_count(nrow(measures_ci_tbl)), " estimates",
        if (is.finite(ci_eligible_n)) {
          paste0("; ", fmt_count(ci_eligible_n), " of ", fmt_count(nrow(measures_ci_tbl)),
                 " estimates have intervals eligible for primary reporting")
        } else {
          ""
        },
        "."
      )
      ci_table_note <- paste0(
        "Report eligible interval limits (", ci_level_text, ", ", ci_method_text,
        ") alongside measures; intervals that do not meet the reporting requirements remain descriptive. "
      )
      results_fit_precision_sentences <- c(results_fit_precision_sentences, ci_sentence)
      results_sentences <- c(results_sentences, ci_sentence)
    }
  }

  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_overall_2 <- extract_overall_pca_second(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  pca_overall_error <- extract_overall_pca_error(pca_obj)
  pca_facet_errors <- extract_facet_pca_errors(pca_obj)
  pca_reference_text <- build_pca_reference_text(warning_threshold_profiles()$pca_reference_bands)
  marginal_state <- extract_strict_marginal_visual_state(diagnostics)

  if (!is.null(pca_overall_1)) {
    ev1 <- to_float(pca_overall_1$Eigenvalue)
    pr1 <- to_float(pca_overall_1$Proportion) * 100
    if (!is.null(pca_overall_2)) {
      ev2 <- to_float(pca_overall_2$Eigenvalue)
      pca_overall_sentence <- paste0(
        "Exploratory residual PCA (overall standardized residual matrix) showed PC1 eigenvalue = ",
        fmt_num(ev1), " (", fmt_num(pr1, 1), "% variance), with PC2 eigenvalue = ", fmt_num(ev2), "."
      )
    } else {
      pca_overall_sentence <- paste0(
        "Exploratory residual PCA (overall standardized residual matrix) showed PC1 eigenvalue = ",
        fmt_num(ev1), " (", fmt_num(pr1, 1), "% variance)."
      )
    }
    results_residual_sentences <- c(results_residual_sentences, pca_overall_sentence)
    results_sentences <- c(results_sentences, pca_overall_sentence)
  } else if (nzchar(pca_overall_error)) {
    unavailable_msg <- paste0("Overall residual PCA was not available for this run: ", pca_overall_error)
    results_residual_sentences <- c(results_residual_sentences, unavailable_msg)
    results_sentences <- c(results_sentences, unavailable_msg)
  }

  if (nrow(pca_facet_1) > 0) {
    top <- pca_facet_1[1, , drop = FALSE]
    facet_pca_sentence <- paste0(
      "Facet-specific exploratory residual PCA showed the largest first-component signal in ",
      as.character(top$Facet), " (eigenvalue = ", fmt_num(top$Eigenvalue),
      ", ", fmt_num(100 * to_float(top$Proportion), 1), "% variance)."
    )
    results_residual_sentences <- c(results_residual_sentences, facet_pca_sentence)
    results_sentences <- c(results_sentences, facet_pca_sentence)
  }
  if (nrow(pca_facet_errors) > 0) {
    unavailable_msg <- paste0(
      "Facet-specific residual PCA was unavailable: ",
      paste(paste0(pca_facet_errors$Facet, " (", pca_facet_errors$Error, ")"), collapse = "; "), "."
    )
    results_residual_sentences <- c(results_residual_sentences, unavailable_msg)
    results_sentences <- c(results_sentences, unavailable_msg)
  }
  if (!is.null(pca_overall_1) || nrow(pca_facet_1) > 0) {
    results_residual_sentences <- c(results_residual_sentences, pca_reference_text)
    results_sentences <- c(results_sentences, pca_reference_text)
  }

  if (isTRUE(marginal_state$marginal_available)) {
    strict_marginal_sentence <- paste0(
      paste(marginal_state$coverage_notes, collapse = " "), " ",
      "Expected counts condition on the same responses through Person posteriors, with fitted calibration held fixed. ",
      "Residual scales omit cross-response covariance and calibration-parameter uncertainty; thresholds are descriptive, not calibrated tests. ",
      "Strict marginal screening gives an overall RMSD of ",
      fmt_num(marginal_state$overall_rmsd),
      ", overall max |standardized residual| = ",
      fmt_num(marginal_state$overall_max_abs_std_residual), "."
    )
    results_residual_sentences <- c(results_residual_sentences, strict_marginal_sentence)
    results_sentences <- c(results_sentences, strict_marginal_sentence)
    if (nrow(marginal_state$top_cell) > 0L) {
      strict_marginal_top_sentence <- paste0(
        "The largest strict marginal cell involved ",
        format_reporting_marginal_cell_label(marginal_state$top_cell),
        " (standardized residual = ", fmt_num(marginal_state$top_cell$StdResidual[1]),
        ", proportion difference = ", fmt_num(marginal_state$top_cell$PropDiff[1]), ")."
      )
      results_residual_sentences <- c(results_residual_sentences, strict_marginal_top_sentence)
      results_sentences <- c(results_sentences, strict_marginal_top_sentence)
    }
  }

  if (isTRUE(marginal_state$pairwise_available)) {
    strict_pairwise_sentence <- paste0(
      "Strict pairwise local-dependence follow-up flagged ",
      fmt_count(marginal_state$pairwise_flagged_pairs),
      " level pair(s) under the latent-integrated agreement screen."
    )
    results_residual_sentences <- c(results_residual_sentences, strict_pairwise_sentence)
    results_sentences <- c(results_sentences, strict_pairwise_sentence)
    if (nrow(marginal_state$top_pair) > 0L) {
      strict_pairwise_top_sentence <- paste0(
        "The largest strict pairwise signal involved ",
        format_reporting_marginal_pair_label(marginal_state$top_pair),
        " (exact-agreement standardized residual = ", fmt_num(marginal_state$top_pair$ExactStdResidual[1]),
        ", adjacent-agreement standardized residual = ", fmt_num(marginal_state$top_pair$AdjacentStdResidual[1]), ")."
      )
      results_residual_sentences <- c(results_residual_sentences, strict_pairwise_top_sentence)
      results_sentences <- c(results_sentences, strict_pairwise_top_sentence)
    }
  }

  if (!is.null(bias_results) && !is.null(bias_results$table) && nrow(bias_results$table) > 0) {
    bias_summary_sentence <- summarize_bias_counts(bias_results)
    results_bias_sentences <- c(results_bias_sentences, bias_summary_sentence)
    results_sentences <- c(results_sentences, bias_summary_sentence)
    bias_tbl <- bias_results$table |> dplyr::filter(is.finite(t))
    bias_sig_n <- sum(is.finite(suppressWarnings(as.numeric(bias_results$table$`Prob.`))) &
      suppressWarnings(as.numeric(bias_results$table$`Prob.`)) < 0.05, na.rm = TRUE)
    bias_facet_pair <- ""
    if (nrow(bias_tbl) > 0) {
      idx <- which.max(abs(bias_tbl$t))
      row <- bias_tbl[idx, , drop = FALSE]
      bias_spec <- extract_bias_facet_spec(bias_results)
      interaction_label <- if (!is.null(bias_spec) && length(bias_spec$facets) > 0) {
        paste(bias_spec$facets, collapse = " x ")
      } else {
        paste0(bias_results$facet_a, " x ", bias_results$facet_b)
      }
      bias_facet_pair <- interaction_label
      bias_detail_sentence <- paste0(
        "Bias analysis for ", interaction_label,
        " showed a largest contrast of ", fmt_num(row$`Bias Size`),
        " logits (screening t = ", fmt_num(row$t),
        ", screening tail area ", fmt_pvalue(row$`Prob.`), ")."
      )
      results_bias_sentences <- c(results_bias_sentences, bias_detail_sentence)
      results_sentences <- c(results_sentences, bias_detail_sentence)
    } else {
      bias_facet_pair <- if (!is.null(bias_results$facet_a) && !is.null(bias_results$facet_b)) {
        paste0(bias_results$facet_a, " x ", bias_results$facet_b)
      } else {
        ""
      }
    }
  } else {
    bias_sig_n <- 0L
    bias_facet_pair <- ""
  }

  bias_caution <- if (!is.null(bias_results) && !is.null(bias_results$table) && nrow(bias_results$table) > 0) {
    "Bias note: bias contrasts and screening tail areas are screening metrics based on plug-in information rather than formal hypothesis tests."
  } else {
    ""
  }
  gpcm_caution <- if (identical(model, "GPCM")) {
    paste(
      "GPCM note: this APA/reporting bundle is a slope-aware",
      "sensitivity-reporting surface over supported diagnostics, direct",
      "tables, and plots. It is not FACETS score-side equivalence, an",
      "automatic operational-scoring decision, or design-forecasting evidence."
    )
  } else {
    ""
  }
  eap_fit_caution <- if (identical(method, "MML")) {
    paste(
      "Fit-basis note: MnSq/ZSTD fit statistics in this run were computed at",
      "EAP person measures, which are shrunken toward the population mean;",
      "they are therefore not numerically interchangeable with JMLE-based",
      "engines such as FACETS. Refit with method = \"JML\" when a JMLE-style",
      "residual basis is required for external comparison."
    )
  } else {
    ""
  }

  contract <- list(
    readiness = list(
      fit = as.data.frame(
        fit_readiness_record$fit, stringsAsFactors = FALSE
      ),
      components = as.data.frame(
        fit_readiness_record$components %||% data.frame(),
        stringsAsFactors = FALSE
      ),
      parameters = as.data.frame(
        fit_readiness_record$parameters %||%
          mfrmr_readiness_empty_parameter_table(),
        stringsAsFactors = FALSE
      )
    ),
    metadata = list(
      model = model,
      method = method,
      n_obs = n_obs,
      weight_sum = weight_sum,
      n_person = n_person,
      n_cat = n_cat,
      rating_min = rating_min,
      rating_max = rating_max,
      facet_names = facet_names,
      facet_counts = facet_counts,
      facets_text = facets_text,
      positive_facets = as.character(config$positive_facets %||% character(0)),
      line_width = line_width,
      # Reflow vs wrapped controls whether section text uses hard line
      # breaks (wrapped, the current 92-col default) or returns one
      # long line per sentence-joined paragraph (reflow, suitable for
      # pasting into Word / RMarkdown / Quarto paragraphs without
      # breaking re-flow).
      output_mode = {
        mode <- as.character(context$output_mode %||% "wrapped")
        if (!mode %in% c("wrapped", "reflow")) "wrapped" else mode
      }
    ),
    context = list(
      assessment = assessment,
      setting = setting,
      rater_training = rater_training,
      raters_per_response = raters_per_response,
      scale_desc = scale_desc,
      rater_facet = rater_facet
    ),
    precision = list(
      tier = precision_tier,
      label = precision_label,
      supports_formal_inference = supports_formal_inference,
      recommended_use = recommended_use,
      reliability_label = reliability_label,
      se_note = se_note,
      caution = precision_caution
    ),
    availability = list(
      has_bias = !is.null(bias_results) && !is.null(bias_results$table) && nrow(bias_results$table) > 0,
      has_pca_overall = !is.null(pca_overall_1),
      has_pca_by_facet = nrow(pca_facet_1) > 0,
      has_strict_marginal = isTRUE(marginal_state$marginal_available),
      has_strict_pairwise = isTRUE(marginal_state$pairwise_available),
      has_rater_reliability = !is.null(rater_rel),
      has_interrater = nrow(interrater_summary) > 0,
      has_measure_ci = ci_available_n > 0,
      has_population_model = isTRUE(population_summary$active),
      has_population_coding = isTRUE(population_summary$active) &&
        nrow(population_summary$coding) > 0L
    ),
    summaries = list(
      threshold_text = threshold_text,
      category_usage = category_usage,
      threshold_coverage = threshold_coverage,
      unused_categories = unused,
      low_count_categories = low_count,
      step_summary = step_summary_sentence,
      person_stats = person_stats,
      overall_fit_infit = infit,
      overall_fit_outfit = outfit,
      misfit_n = misfit_n,
      misfit_total = misfit_total,
      fit_screening = fit_screening,
      top_misfit_sentence = top_misfit_sentence,
      reliability = as.data.frame(rel_tbl %||% data.frame(), stringsAsFactors = FALSE),
      rater_facet = rater_facet,
      rater_reliability = rater_rel,
      interrater_summary = interrater_summary,
      interrater_sentence = interrater_sentence,
      strict_marginal_coverage = marginal_state$coverage,
      strict_marginal_status = marginal_state$marginal_status,
      strict_marginal_reason = marginal_state$marginal_reason,
      strict_marginal_overall_rmsd = marginal_state$overall_rmsd,
      strict_marginal_max_abs_std_residual = marginal_state$overall_max_abs_std_residual,
      strict_marginal_top_cell_label = if (nrow(marginal_state$top_cell) > 0L) format_reporting_marginal_cell_label(marginal_state$top_cell) else "",
      strict_marginal_top_cell_std_residual = if (nrow(marginal_state$top_cell) > 0L) to_float(marginal_state$top_cell$StdResidual[1]) else NA_real_,
      strict_marginal_top_cell_prop_diff = if (nrow(marginal_state$top_cell) > 0L) to_float(marginal_state$top_cell$PropDiff[1]) else NA_real_,
      strict_pairwise_status = marginal_state$pairwise_status,
      strict_pairwise_reason = marginal_state$pairwise_reason,
      strict_pairwise_flagged_pairs = marginal_state$pairwise_flagged_pairs,
      strict_pairwise_top_pair_label = if (nrow(marginal_state$top_pair) > 0L) format_reporting_marginal_pair_label(marginal_state$top_pair) else "",
      strict_pairwise_top_pair_exact_std_residual = if (nrow(marginal_state$top_pair) > 0L) to_float(marginal_state$top_pair$ExactStdResidual[1]) else NA_real_,
      strict_pairwise_top_pair_adjacent_std_residual = if (nrow(marginal_state$top_pair) > 0L) to_float(marginal_state$top_pair$AdjacentStdResidual[1]) else NA_real_,
      pca_overall_1 = pca_overall_1,
      pca_overall_2 = pca_overall_2,
      pca_by_facet_first = pca_facet_1,
      pca_reference_text = pca_reference_text,
      ci_availability_sentence = ci_sentence,
      ci_table_note = ci_table_note,
      ci_available_n = ci_available_n,
      bias_summary = summarize_bias_counts(bias_results),
      bias_sig_n = bias_sig_n,
      bias_facet_pair = bias_facet_pair,
      population_model = population_summary
    ),
    method_design_sentences = method_design_sentences,
    method_estimation_sentences = method_estimation_sentences,
    results_scale_sentences = results_scale_sentences,
    results_measure_sentences = results_measure_sentences,
    results_population_sentences = results_population_sentences,
    results_fit_precision_sentences = results_fit_precision_sentences,
    results_residual_sentences = results_residual_sentences,
    results_bias_sentences = results_bias_sentences,
    method_sentences = method_sentences,
    results_sentences = results_sentences,
    caution_sentences = Filter(nzchar, c(
      precision_caution,
      eap_fit_caution,
      bias_caution,
      gpcm_caution,
      population_summary$caution_sentence,
      population_summary$conquest_sentence
    )),
    section_order = c(
      "method_design", "method_estimation",
      "results_scale", "results_measures", "results_population_model", "results_fit_precision",
      "results_residual_structure", "results_bias_screening", "results_cautions"
    ),
    ordered_keys = build_apa_table_figure_key_order(include_population = population_summary$active)
  )

  contract$section_map <- build_apa_section_map_from_contract(contract)
  contract$section_table <- flatten_apa_section_map(contract$section_map, order = contract$section_order)
  contract$note_map <- build_apa_note_map_from_contract(contract)
  contract$caption_map <- build_apa_caption_map_from_contract(contract)
  contract$note_text <- paste(
    vapply(contract$ordered_keys[contract$ordered_keys %in% names(contract$note_map)], function(k) contract$note_map[[k]], character(1)),
    collapse = "\n\n"
  )
  contract$caption_text <- paste(
    vapply(contract$ordered_keys[contract$ordered_keys %in% names(contract$caption_map)], function(k) contract$caption_map[[k]], character(1)),
    collapse = "\n\n"
  )
  contract$report_text <- build_apa_report_text_from_contract(contract)
  class(contract) <- c("mfrm_apa_contract", "list")
  contract
}





build_visual_warning_map <- function(res,
                                     diagnostics,
                                     whexact = FALSE,
                                     thresholds = NULL,
                                     threshold_profile = "standard") {
  # Plot-level warning text is accumulated independently for each visual output.
  visual_keys <- c(
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds", "category_curves",
    "observed_expected", "fit_diagnostics", "fit_zstd_distribution", "misfit_levels",
    "strict_marginal_fit", "strict_pairwise_local_dependence",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  warnings <- stats::setNames(replicate(length(visual_keys), character(0), simplify = FALSE), visual_keys)
  if (is.null(res) || is.null(diagnostics)) return(warnings)

  # Stage 1: Resolve active threshold profile (strict/standard/lenient + overrides).
  resolved <- resolve_warning_thresholds(thresholds = thresholds, threshold_profile = threshold_profile)
  active <- resolved$thresholds
  profile_name <- resolved$profile_name
  pca_reference_text <- build_pca_reference_text(resolved$pca_reference_bands)

  n_obs_min <- active$n_obs_min %||% 100
  n_person_min <- active$n_person_min %||% 30
  low_cat_min <- active$low_cat_min %||% 10
  min_facet_levels <- active$min_facet_levels %||% 3
  misfit_ratio_warn <- active$misfit_ratio_warn %||% 0.10
  missing_fit_ratio_warn <- active$missing_fit_ratio_warn %||% 0.20
  zstd2_ratio_warn <- active$zstd2_ratio_warn %||% 0.10
  zstd3_ratio_warn <- active$zstd3_ratio_warn %||% 0.05
  expected_var_min <- active$expected_var_min %||% 0.20
  pca_first_eigen_warn <- active$pca_first_eigen_warn %||% 2.0
  pca_first_prop_warn <- active$pca_first_prop_warn %||% 0.10

  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  n_obs <- to_float(res$prep$n_obs %||% summary$N)
  n_person <- if (!is.null(res$facets$person)) nrow(res$facets$person) else 0

  # Stage 2: Sample-size and design warnings (Wright/pathway/observed-expected plots).
  if (is.finite(n_obs) && n_obs < n_obs_min) {
    warnings$wright_map <- c(warnings$wright_map, paste0("Small number of observations (N = ", fmt_count(n_obs), " < ", fmt_count(n_obs_min), ")."))
    warnings$pathway_map <- c(warnings$pathway_map, paste0("Small number of observations (N = ", fmt_count(n_obs), " < ", fmt_count(n_obs_min), "); pathway curves may be unstable."))
    warnings$observed_expected <- c(warnings$observed_expected, paste0("Small number of observations (N = ", fmt_count(n_obs), " < ", fmt_count(n_obs_min), "); bin averages may be noisy."))
  }

  if (n_person < n_person_min) {
    warnings$wright_map <- c(warnings$wright_map, paste0("Small person sample (n = ", fmt_count(n_person), " < ", fmt_count(n_person_min), "); interpret spread cautiously."))
  }

  facet_levels <- res$config$facet_levels
  small_facets <- names(facet_levels)[vapply(facet_levels, length, integer(1)) < min_facet_levels]
  if (length(small_facets) > 0) {
    warnings$wright_map <- c(warnings$wright_map, paste0("Facets with very few levels: ", paste(small_facets, collapse = ", "), "."))
    warnings$facet_distribution <- c(warnings$facet_distribution, paste0("Facet distributions are based on few levels: ", paste(small_facets, collapse = ", "), "."))
  }

  # Stage 3: Rating-scale warnings.
  cat_tbl <- calc_category_stats(diagnostics$obs, res = res, whexact = whexact)
  usage <- summarize_category_usage(cat_tbl, low_count = low_cat_min)
  if (!usage$Categories || usage$UnavailableCounts > 0L || isTRUE(usage$LowCountCategories > 0L)) {
    warnings$category_curves <- c(warnings$category_curves, category_usage_note(usage, low_cat_min))
  }
  step_order <- calc_step_order(res$steps, expected_steps = res$config$n_cat - 1L,
    expected_facets = res$config$facet_levels[[res$config$step_facet %||% ""]])
  coverage <- summarize_threshold_order(step_order)
  if (!isTRUE(coverage$NotApplicable) &&
      (!isTRUE(coverage$Unavailable == 0L) || isTRUE(coverage$Decreasing > 0L))) {
    warnings$step_thresholds <- c(warnings$step_thresholds, threshold_order_note(coverage))
  }

  # Stage 4: Fit-based warnings.
  measures <- diagnostics$measures
  if (is.null(measures) || nrow(measures) == 0) {
    warnings$pathway_map <- c(warnings$pathway_map, "Fit statistics are not available for this run.")
    warnings$fit_diagnostics <- c(warnings$fit_diagnostics, "Fit statistics are not available for this run.")
    warnings$fit_zstd_distribution <- c(warnings$fit_zstd_distribution, "ZSTD distributions are not available for this run.")
    warnings$misfit_levels <- c(warnings$misfit_levels, "Misfit ranking requires fit statistics.")
    return(warnings)
  }

  screening <- fit_screening_summary(measures)$counts
  limits <- c(misfit_ratio_warn, zstd2_ratio_warn, zstd3_ratio_warn)
  for (i in seq_len(nrow(screening))) {
    row <- screening[i, ]
    key <- if (i == 1L) "fit_diagnostics" else "fit_zstd_distribution"
    if (row$IncompleteStatistics > 0L ||
        (is.finite(row$FlagRate) && row$FlagRate > limits[i])) {
      warnings[[key]] <- c(warnings[[key]], fit_screening_note(row))
    }
  }
  if (isTRUE(screening$IncompleteStatistics[1] / screening$Elements[1] >= missing_fit_ratio_warn)) {
    warnings$fit_diagnostics <- c(warnings$fit_diagnostics,
      "The share of incomplete Infit/Outfit pairs reaches the configured missing-fit warning threshold.")
  }

  obs <- diagnostics$obs
  if (!is.null(obs) && nrow(obs) > 0 && "Expected" %in% names(obs)) {
    exp_var <- stats::var(suppressWarnings(as.numeric(obs$Expected)), na.rm = TRUE)
    if (is.finite(exp_var) && exp_var < expected_var_min) {
      warnings$observed_expected <- c(warnings$observed_expected, "Expected scores have limited spread; trends may be muted.")
    }
  }

  # Stage 4b: Strict marginal fit warnings.
  marginal_state <- extract_strict_marginal_visual_state(diagnostics)
  if (!isTRUE(marginal_state$marginal_available)) {
    warnings$strict_marginal_fit <- c(
      warnings$strict_marginal_fit,
      paste0("Strict marginal diagnostics are not available: ", marginal_state$marginal_reason)
    )
  } else {
    warnings$strict_marginal_fit <- c(
      warnings$strict_marginal_fit,
      "Strict marginal diagnostics are exploratory posterior-expected screens, not formal inferential tests.",
      marginal_state$coverage_notes
    )
    flagged_groups <- sum(
      c(marginal_state$step_groups_flagged, marginal_state$facet_levels_flagged),
      na.rm = TRUE
    )
    if (is.finite(flagged_groups) && flagged_groups > 0) {
      warnings$strict_marginal_fit <- c(
        warnings$strict_marginal_fit,
        paste0(
          "Flagged strict marginal groups were detected (step/scale = ",
          fmt_count(marginal_state$step_groups_flagged),
          ", facet levels = ", fmt_count(marginal_state$facet_levels_flagged), ")."
        )
      )
    }
    rmsd_warn <- to_float(marginal_state$marginal_thresholds$rmsd_warn %||% 0.05)
    if (is.finite(marginal_state$overall_rmsd) && is.finite(rmsd_warn) && marginal_state$overall_rmsd >= rmsd_warn) {
      warnings$strict_marginal_fit <- c(
        warnings$strict_marginal_fit,
        paste0(
          "Overall strict marginal RMSD = ", fmt_num(marginal_state$overall_rmsd),
          " exceeds the current screening band (", fmt_num(rmsd_warn), ")."
        )
      )
    }
    if (nrow(marginal_state$top_cell) > 0L) {
      warnings$strict_marginal_fit <- c(
        warnings$strict_marginal_fit,
        paste0(
          "Inspect the largest cell with plot(x, type = \"marginal_fit\"): ",
          format_reporting_marginal_cell_label(marginal_state$top_cell),
          " (|standardized residual| = ", fmt_num(abs(marginal_state$top_cell$StdResidual[1])), ")."
        )
      )
    }
  }

  if (!isTRUE(marginal_state$pairwise_available)) {
    warnings$strict_pairwise_local_dependence <- c(
      warnings$strict_pairwise_local_dependence,
      paste0("Strict pairwise local-dependence diagnostics are not available: ", marginal_state$pairwise_reason)
    )
  } else {
    warnings$strict_pairwise_local_dependence <- c(
      warnings$strict_pairwise_local_dependence,
      "Strict pairwise local-dependence diagnostics are exploratory follow-ups to first-order strict marginal flags."
    )
    if (is.finite(marginal_state$pairwise_flagged_pairs) && marginal_state$pairwise_flagged_pairs > 0) {
      warnings$strict_pairwise_local_dependence <- c(
        warnings$strict_pairwise_local_dependence,
        paste0(
          "Flagged strict pairwise level pairs were detected (n = ",
          fmt_count(marginal_state$pairwise_flagged_pairs), ")."
        )
      )
    }
    if (nrow(marginal_state$top_pair) > 0L) {
      top_pair_max_abs <- max(
        c(
          abs(to_float(marginal_state$top_pair$ExactStdResidual[1])),
          abs(to_float(marginal_state$top_pair$AdjacentStdResidual[1]))
        ),
        na.rm = TRUE
      )
      if (!is.finite(top_pair_max_abs)) top_pair_max_abs <- NA_real_
      warnings$strict_pairwise_local_dependence <- c(
        warnings$strict_pairwise_local_dependence,
        paste0(
          "Inspect the largest pair with plot(x, type = \"marginal_pairwise\"): ",
          format_reporting_marginal_pair_label(marginal_state$top_pair),
          " (maximum available |standardized residual| = ",
          fmt_num(top_pair_max_abs),
          ")."
        )
      )
    }
  }

  # Stage 5: Residual PCA warnings (overall and by facet).
  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  pca_overall_error <- extract_overall_pca_error(pca_obj)
  pca_facet_errors <- extract_facet_pca_errors(pca_obj)
  warnings$residual_pca_overall <- c(
    warnings$residual_pca_overall,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )
  warnings$residual_pca_by_facet <- c(
    warnings$residual_pca_by_facet,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )

  if (is.null(pca_overall_1)) {
    warnings$residual_pca_overall <- c(
      warnings$residual_pca_overall,
      if (nzchar(pca_overall_error)) paste0("Overall residual PCA is not available: ", pca_overall_error) else "Overall residual PCA is not available."
    )
  } else {
    warnings$residual_pca_overall <- c(
      warnings$residual_pca_overall,
      build_pca_check_text(
        eigenvalue = pca_overall_1$Eigenvalue,
        proportion = pca_overall_1$Proportion,
        reference_bands = resolved$pca_reference_bands
      )
    )
    if (to_float(pca_overall_1$Eigenvalue) > pca_first_eigen_warn) {
      warnings$residual_pca_overall <- c(
        warnings$residual_pca_overall,
        paste0("Overall residual PCA PC1 exceeds the current heuristic eigenvalue band (", fmt_num(pca_overall_1$Eigenvalue), ").")
      )
    }
    if (to_float(pca_overall_1$Proportion) > pca_first_prop_warn) {
      warnings$residual_pca_overall <- c(
        warnings$residual_pca_overall,
        paste0("Overall residual PCA PC1 explains ", fmt_num(100 * to_float(pca_overall_1$Proportion), 1), "% variance.")
      )
    }
  }

  if (nrow(pca_facet_1) == 0) {
    warnings$residual_pca_by_facet <- c(
      warnings$residual_pca_by_facet,
      if (nrow(pca_facet_errors) > 0) {
        paste0(
          "Facet-specific residual PCA is not available: ",
          paste(paste0(pca_facet_errors$Facet, " (", pca_facet_errors$Error, ")"), collapse = "; ")
        )
      } else {
        "Facet-specific residual PCA is not available."
      }
    )
  } else {
    top <- pca_facet_1[1, , drop = FALSE]
    warnings$residual_pca_by_facet <- c(
      warnings$residual_pca_by_facet,
      paste0(
        "Top facet PC1 (", as.character(top$Facet), "): ",
        build_pca_check_text(
          eigenvalue = top$Eigenvalue,
          proportion = top$Proportion,
          reference_bands = resolved$pca_reference_bands
        )
      )
    )
    flagged <- pca_facet_1[
      (suppressWarnings(as.numeric(pca_facet_1$Eigenvalue)) > pca_first_eigen_warn) |
        (suppressWarnings(as.numeric(pca_facet_1$Proportion)) > pca_first_prop_warn),
      ,
      drop = FALSE
    ]
    if (nrow(flagged) > 0) {
      warnings$residual_pca_by_facet <- c(
        warnings$residual_pca_by_facet,
        paste0("Facet residual PCA shows stronger exploratory PC1 signal in: ", paste(flagged$Facet, collapse = ", "), ".")
      )
    }
  }

  warnings
}

build_visual_summary_map <- function(res,
                                     diagnostics,
                                     whexact = FALSE,
                                     options = NULL,
                                     thresholds = NULL,
                                     threshold_profile = "standard") {
  # Summary map mirrors the warning map but provides descriptive, non-binary text.
  options <- options %||% list()
  detail <- tolower(as.character(options$detail %||% "standard"))
  max_facet_ranges <- as.integer(options$max_facet_ranges %||% 4)
  top_misfit_n <- as.integer(options$top_misfit_n %||% 3)
  include_top_misfit <- top_misfit_n > 0
  resolved <- resolve_warning_thresholds(thresholds = thresholds, threshold_profile = threshold_profile)
  active <- resolved$thresholds
  profile_name <- resolved$profile_name
  pca_reference_text <- build_pca_reference_text(resolved$pca_reference_bands)
  pca_first_eigen_warn <- active$pca_first_eigen_warn %||% 2.0
  pca_first_prop_warn <- active$pca_first_prop_warn %||% 0.10

  visual_keys <- c(
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds", "category_curves",
    "observed_expected", "fit_diagnostics", "fit_zstd_distribution", "misfit_levels",
    "strict_marginal_fit", "strict_pairwise_local_dependence",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  summaries <- stats::setNames(replicate(length(visual_keys), character(0), simplify = FALSE), visual_keys)
  if (is.null(res) || is.null(diagnostics)) return(summaries)

  # Stage 1: Global design summary.
  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  n_obs <- to_float(res$prep$n_obs %||% summary$N)
  n_person <- if (!is.null(res$facets$person)) nrow(res$facets$person) else 0

  if (is.finite(n_obs)) summaries$wright_map <- c(summaries$wright_map, paste0("Observations: N = ", fmt_count(n_obs), "."))
  summaries$wright_map <- c(summaries$wright_map, paste0("Persons: n = ", fmt_count(n_person), "."))

  person_stats <- describe_series(res$facets$person$Estimate)
  if (!is.null(person_stats)) {
    summaries$wright_map <- c(summaries$wright_map,
      paste0("Person range ", fmt_num(person_stats$min), " to ", fmt_num(person_stats$max),
             " (M = ", fmt_num(person_stats$mean), ", SD = ", fmt_num(person_stats$sd), ").")
    )
    if (detail == "detailed" && is.finite(person_stats$sd)) {
      summaries$wright_map <- c(summaries$wright_map, paste0("Person spread (SD) = ", fmt_num(person_stats$sd), " logits."))
    }
  }

  if (!is.null(res$facets$others) && nrow(res$facets$others) > 0) {
    facet_stats <- c()
    for (facet in unique(res$facets$others$Facet)) {
      df <- res$facets$others |> dplyr::filter(Facet == facet)
      st <- describe_series(df$Estimate)
      if (!is.null(st)) {
        facet_stats <- c(facet_stats, paste0(facet, ": n = ", fmt_count(nrow(df)), ", range ", fmt_num(st$min), " to ", fmt_num(st$max)))
      }
    }
    if (length(facet_stats) > 0) {
      summaries$wright_map <- c(summaries$wright_map, paste0("Facet ranges: ", paste(head(facet_stats, max_facet_ranges), collapse = "; "), "."))
      if (length(facet_stats) > max_facet_ranges) summaries$wright_map <- c(summaries$wright_map, "Additional facets omitted for brevity.")
    }
  }

  # Stage 2: Fit and category summaries.
  measures <- diagnostics$measures
  step_tbl <- res$steps
  step_order <- calc_step_order(step_tbl, expected_steps = res$config$n_cat - 1L,
    expected_facets = res$config$facet_levels[[res$config$step_facet %||% ""]])
  order_note <- threshold_order_note(summarize_threshold_order(step_order))
  summaries$pathway_map <- c(summaries$pathway_map, order_note)
  summaries$step_thresholds <- c(summaries$step_thresholds, order_note)
  if (!is.null(res$facets$others) && nrow(res$facets$others) > 0) {
    summaries$facet_distribution <- c(summaries$facet_distribution, "Distributions show the spread of severity/difficulty within each facet.")
  }
  cat_tbl <- calc_category_stats(diagnostics$obs, res = res, whexact = whexact)
  summaries$category_curves <- c(summaries$category_curves,
    category_usage_note(summarize_category_usage(cat_tbl)))

  obs <- diagnostics$obs
  if (!is.null(obs) && nrow(obs) > 0 && all(c("Observed", "Expected") %in% names(obs))) {
    resid <- suppressWarnings(as.numeric(obs$Observed) - as.numeric(obs$Expected))
    if ("Weight" %in% names(obs)) {
      w <- suppressWarnings(as.numeric(obs$Weight))
      w <- ifelse(is.finite(w) & w > 0, w, 0)
      mean_resid <- if (sum(w) > 0) sum(resid * w, na.rm = TRUE) / sum(w) else NA_real_
      mae <- if (sum(w) > 0) sum(abs(resid) * w, na.rm = TRUE) / sum(w) else NA_real_
    } else {
      mean_resid <- mean(resid, na.rm = TRUE)
      mae <- mean(abs(resid), na.rm = TRUE)
    }
    summaries$observed_expected <- c(summaries$observed_expected, paste0("Mean residual: ", fmt_num(mean_resid), "."))
    summaries$observed_expected <- c(summaries$observed_expected, paste0("Mean absolute residual: ", fmt_num(mae), "."))
  }

  screening <- fit_screening_summary(measures)
  summaries$fit_diagnostics <- c(summaries$fit_diagnostics,
    fit_screening_note(screening$counts[1, ]))
  summaries$fit_zstd_distribution <- c(summaries$fit_zstd_distribution,
    fit_screening_note(screening$counts[2, ]), fit_screening_note(screening$counts[3, ]))
  if (!is.null(measures) && nrow(measures) > 0) {
    if (detail == "detailed") {
      for (name in c("Infit", "Outfit")) {
        values <- suppressWarnings(as.numeric(measures[[name]]))
        valid <- is.finite(values) & values >= 0
        summaries$fit_diagnostics <- c(summaries$fit_diagnostics,
          sprintf("Mean %s = %s (%d of %d elements available).", name,
                  fmt_num(if (any(valid)) mean(values[valid]) else NA_real_), sum(valid), nrow(measures)))
      }
    }
    if (include_top_misfit) {
      summaries$misfit_levels <- c(summaries$misfit_levels,
        summarize_top_misfit_levels(measures, top_n = top_misfit_n))
    }
  }

  # Stage 2b: Strict marginal fit summaries.
  marginal_state <- extract_strict_marginal_visual_state(diagnostics)
  if (!isTRUE(marginal_state$marginal_available)) {
    summaries$strict_marginal_fit <- c(
      summaries$strict_marginal_fit,
      paste0("Strict marginal diagnostics unavailable: ", marginal_state$marginal_reason)
    )
  } else {
    summaries$strict_marginal_fit <- c(
      summaries$strict_marginal_fit,
      "Expected counts condition on the same responses through Person posteriors, with fitted calibration held fixed.",
      "Residual scales omit cross-response covariance and calibration uncertainty; cutoffs are descriptive review rules.",
      marginal_state$coverage_notes
    )
    summaries$strict_marginal_fit <- c(
      summaries$strict_marginal_fit,
      paste0(
        "Overall RMSD = ", fmt_num(marginal_state$overall_rmsd),
        "; overall max |standardized residual| = ",
        fmt_num(marginal_state$overall_max_abs_std_residual), "."
      )
    )
    summaries$strict_marginal_fit <- c(
      summaries$strict_marginal_fit,
      paste0(
        "Flagged step/scale groups = ", fmt_count(marginal_state$step_groups_flagged),
        "; flagged facet levels = ", fmt_count(marginal_state$facet_levels_flagged), "."
      )
    )
    if (nrow(marginal_state$top_cell) > 0L) {
      summaries$strict_marginal_fit <- c(
        summaries$strict_marginal_fit,
        paste0(
          "Top cell: ", format_reporting_marginal_cell_label(marginal_state$top_cell),
          " (standardized residual = ", fmt_num(marginal_state$top_cell$StdResidual[1]),
          ", proportion difference = ", fmt_num(marginal_state$top_cell$PropDiff[1]), ")."
        )
      )
    }
    summaries$strict_marginal_fit <- c(
      summaries$strict_marginal_fit,
      "Use plot(x, type = \"marginal_fit\") to inspect the highest-residual category cells."
    )
  }

  if (!isTRUE(marginal_state$pairwise_available)) {
    summaries$strict_pairwise_local_dependence <- c(
      summaries$strict_pairwise_local_dependence,
      paste0("Strict pairwise local-dependence diagnostics unavailable: ", marginal_state$pairwise_reason)
    )
  } else {
    summaries$strict_pairwise_local_dependence <- c(
      summaries$strict_pairwise_local_dependence,
      "Strict pairwise local dependence is an exploratory second-order follow-up."
    )
    summaries$strict_pairwise_local_dependence <- c(
      summaries$strict_pairwise_local_dependence,
      paste0("Flagged level pairs = ", fmt_count(marginal_state$pairwise_flagged_pairs), ".")
    )
    if (nrow(marginal_state$top_pair) > 0L) {
      summaries$strict_pairwise_local_dependence <- c(
        summaries$strict_pairwise_local_dependence,
        paste0(
          "Top pair: ", format_reporting_marginal_pair_label(marginal_state$top_pair),
          " (exact-agreement standardized residual = ", fmt_num(marginal_state$top_pair$ExactStdResidual[1]),
          ", adjacent-agreement standardized residual = ", fmt_num(marginal_state$top_pair$AdjacentStdResidual[1]), ")."
        )
      )
    }
    summaries$strict_pairwise_local_dependence <- c(
      summaries$strict_pairwise_local_dependence,
      "Use plot(x, type = \"marginal_pairwise\") to inspect exact and adjacent agreement gaps."
    )
  }

  # Stage 3: Residual PCA summaries with threshold profile context.
  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_overall_2 <- extract_overall_pca_second(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  pca_overall_error <- extract_overall_pca_error(pca_obj)
  pca_facet_errors <- extract_facet_pca_errors(pca_obj)
  summaries$residual_pca_overall <- c(
    summaries$residual_pca_overall,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )
  summaries$residual_pca_by_facet <- c(
    summaries$residual_pca_by_facet,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )

  if (!is.null(pca_overall_1)) {
    summaries$residual_pca_overall <- c(
      summaries$residual_pca_overall,
      paste0(
        "Overall residual PCA PC1: eigenvalue = ", fmt_num(pca_overall_1$Eigenvalue),
        ", variance = ", fmt_num(100 * to_float(pca_overall_1$Proportion), 1), "%."
      )
    )
    summaries$residual_pca_overall <- c(
      summaries$residual_pca_overall,
      build_pca_check_text(
        eigenvalue = pca_overall_1$Eigenvalue,
        proportion = pca_overall_1$Proportion,
        reference_bands = resolved$pca_reference_bands
      )
    )
    if (!is.null(pca_overall_2)) {
      summaries$residual_pca_overall <- c(
        summaries$residual_pca_overall,
        paste0("Overall residual PCA PC2: eigenvalue = ", fmt_num(pca_overall_2$Eigenvalue), ".")
      )
    }
  } else {
    summaries$residual_pca_overall <- c(
      summaries$residual_pca_overall,
      if (nzchar(pca_overall_error)) paste0("Overall residual PCA unavailable: ", pca_overall_error) else "Overall residual PCA unavailable."
    )
  }

  if (nrow(pca_facet_1) > 0) {
    show_n <- min(nrow(pca_facet_1), ifelse(detail == "detailed", 5L, 3L))
    top <- pca_facet_1[seq_len(show_n), , drop = FALSE]
    summaries$residual_pca_by_facet <- c(
      summaries$residual_pca_by_facet,
      paste0(
        "Top facet PC1 (", as.character(top$Facet[1]), "): ",
        build_pca_check_text(
          eigenvalue = top$Eigenvalue[1],
          proportion = top$Proportion[1],
          reference_bands = resolved$pca_reference_bands
        )
      )
    )
    labels <- vapply(seq_len(nrow(top)), function(i) {
      paste0(
        top$Facet[i], "=", fmt_num(top$Eigenvalue[i]),
        " (", fmt_num(100 * to_float(top$Proportion[i]), 1), "%)"
      )
    }, character(1))
    summaries$residual_pca_by_facet <- c(
      summaries$residual_pca_by_facet,
      paste0("Facet residual PCA PC1 signals: ", paste(labels, collapse = "; "), ".")
    )
    if (nrow(pca_facet_1) > show_n) {
      summaries$residual_pca_by_facet <- c(summaries$residual_pca_by_facet, "Additional facets omitted for brevity.")
    }
  } else {
    summaries$residual_pca_by_facet <- c(
      summaries$residual_pca_by_facet,
      if (nrow(pca_facet_errors) > 0) {
        paste0(
          "Facet-specific residual PCA unavailable: ",
          paste(paste0(pca_facet_errors$Facet, " (", pca_facet_errors$Error, ")"), collapse = "; ")
        )
      } else {
        "Facet-specific residual PCA unavailable."
      }
    )
  }

  summaries
}
