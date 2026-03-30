#' Build an auto-filled MFRM reporting checklist
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are computed with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()] or a named list of
#'   such outputs.
#' @param include_references If `TRUE`, include a compact reference table in the
#'   returned bundle.
#'
#' @details
#' This helper ports the app-level reporting checklist into a package-native
#' bundle. It does not try to judge substantive reporting quality; instead, it
#' checks whether the fitted object and related diagnostics contain the evidence
#' typically reported in MFRM write-ups.
#'
#' Checklist items are grouped into seven sections:
#' - Method section
#' - Global fit
#' - Facet-level statistics
#' - Element-level statistics
#' - Rating scale diagnostics
#' - Bias/interaction analysis
#' - Visual displays
#'
#' The output is designed for manuscript preparation, audit trails, and
#' reproducible reporting workflows.
#'
#' @section What this checklist means:
#' `reporting_checklist()` is a content-availability contract. It tells you
#' which reporting elements are already present in the current analysis
#' objects and which still need to be generated or documented. The primary
#' draft-status column is `DraftReady`; `ReadyForAPA` is retained as a
#' backward-compatible alias.
#'
#' @section What this checklist does not justify:
#' - It is not a single run-level pass/fail decision for publication.
#' - `DraftReady = TRUE` / `ReadyForAPA = TRUE` does not certify formal
#'   inferential adequacy.
#' - Missing bias rows may simply mean `bias_results` were not supplied.
#'
#' @section Interpreting output:
#' - `checklist`: one row per reporting item with `Available = TRUE/FALSE`.
#'   `DraftReady = TRUE` means the item can be drafted into a report with the
#'   package's documented caveats. `ReadyForAPA` is a backward-compatible alias
#'   of the same flag; neither field certifies formal inferential adequacy.
#' - `section_summary`: available items by section.
#' - `references`: core background references when requested.
#'
#' @section Recommended next step:
#' Review the rows with `Available = FALSE` or `DraftReady = FALSE`, then add
#' the missing diagnostics, bias results, or narrative context before calling
#' [build_apa_outputs()] for draft text generation.
#'
#' @section Typical workflow:
#' 1. Fit with [fit_mfrm()].
#' 2. Compute diagnostics with [diagnose_mfrm()].
#' 3. Run `reporting_checklist()` to see which reporting elements are already
#'    available from the current analysis objects.
#'
#' @return A named list with checklist tables. Class:
#'   `mfrm_reporting_checklist`.
#' @seealso [build_apa_outputs()], [build_visual_summaries()],
#'   [specifications_report()], [data_quality_report()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' summary(chk)
#' apa <- build_apa_outputs(fit, diag)
#' head(chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#' nchar(apa$report_text)
#' }
#' @export
reporting_checklist <- function(fit,
                                diagnostics = NULL,
                                bias_results = NULL,
                                include_references = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (!is.list(diagnostics) || is.null(diagnostics$obs)) {
    stop("`diagnostics` must be output from diagnose_mfrm().", call. = FALSE)
  }

  config <- fit$config %||% list()
  prep <- fit$prep %||% list()
  measures <- as.data.frame(diagnostics$measures %||% data.frame(), stringsAsFactors = FALSE)
  obs_df <- as.data.frame(diagnostics$obs %||% data.frame(), stringsAsFactors = FALSE)
  rel_df <- as.data.frame(diagnostics$reliability %||% data.frame(), stringsAsFactors = FALSE)
  precision_profile_df <- as.data.frame(diagnostics$precision_profile %||% data.frame(), stringsAsFactors = FALSE)
  steps_df <- as.data.frame(fit$steps %||% data.frame(), stringsAsFactors = FALSE)
  pca_obj <- diagnostics$pca %||%
    diagnostics$residual_pca_overall %||%
    diagnostics$residual_pca_by_facet %||%
    NULL
  bias_tbls <- collect_bias_tables_for_checklist(bias_results)
  bias_error_tbl <- attr(bias_tbls, "errors", exact = TRUE)
  n_bias_errors <- if (is.data.frame(bias_error_tbl)) nrow(bias_error_tbl) else 0L

  converged <- isTRUE(fit$summary$Converged %||% FALSE)
  summary_msg <- if (is.data.frame(fit$summary) && "Message" %in% names(fit$summary)) {
    as.character(fit$summary$Message[1] %||% "")
  } else {
    ""
  }
  conv_msg <- as.character(fit$opt$message %||% summary_msg %||% "")
  n_bias_pairs <- length(bias_tbls)
  has_bias_sig <- FALSE
  bias_stat_parse_issue <- FALSE
  if (n_bias_pairs > 0) {
    has_bias_sig <- any(vapply(bias_tbls, function(tbl) {
      if (!is.data.frame(tbl) || nrow(tbl) == 0 || !"t" %in% names(tbl)) return(FALSE)
      raw_t <- as.character(tbl$t)
      t_vals <- suppressWarnings(as.numeric(raw_t))
      bad_t <- is.na(t_vals) & !is.na(raw_t) & nzchar(trimws(raw_t))
      if (any(bad_t)) {
        bias_stat_parse_issue <<- TRUE
      }
      any(is.finite(t_vals) & abs(t_vals) >= 2)
    }, logical(1)))
  }

  has_resid <- nrow(obs_df) > 0 && "StdResidual" %in% names(obs_df)
  has_meas <- nrow(measures) > 0 && "Estimate" %in% names(measures)
  has_fit <- has_meas && all(c("Infit", "Outfit") %in% names(measures))
  has_ci <- has_meas && "SE" %in% names(measures)
  has_rel <- nrow(rel_df) > 0
  precision_tier <- as.character(precision_profile_df$PrecisionTier[1] %||% NA_character_)
  formal_precision <- isTRUE(precision_profile_df$SupportsFormalInference[1] %||% FALSE) && converged
  measure_formal <- if (has_meas && "SupportsFormalInference" %in% names(measures)) {
    any(as.logical(measures$SupportsFormalInference), na.rm = TRUE)
  } else {
    formal_precision
  }
  ci_ready <- if (has_meas && "CIEligible" %in% names(measures)) {
    any(as.logical(measures$CIEligible), na.rm = TRUE)
  } else {
    has_ci && formal_precision
  }
  has_steps <- nrow(steps_df) > 0
  has_pca <- !is.null(pca_obj) && length(pca_obj) > 0
  has_counts <- has_resid && "Observed" %in% names(obs_df)
  has_person_measure <- has_resid && "PersonMeasure" %in% names(obs_df)

  add_item <- function(section,
                       item,
                       available,
                       detail,
                       source_component,
                       severity = c("required", "recommended", "optional"),
                       ready_for_apa = available,
                       missing_action = "Compute or document this component before manuscript export.",
                       available_action = NULL) {
    severity <- match.arg(severity)
    ready_for_apa <- isTRUE(ready_for_apa)
    priority <- if (ready_for_apa) {
      "ready"
    } else {
      switch(
        severity,
        required = "high",
        recommended = "medium",
        optional = "low",
        "medium"
      )
    }
    if (is.null(available_action) || !nzchar(as.character(available_action))) {
      available_action <- if (ready_for_apa) {
        "Available; adapt this evidence into the manuscript draft after methodological review."
      } else {
        "Available, but keep the documented cautionary language when drafting."
      }
    }
    data.frame(
      Section = as.character(section),
      Item = as.character(item),
      Available = isTRUE(available),
      DraftReady = ready_for_apa,
      ReadyForAPA = ready_for_apa,
      Severity = as.character(severity),
      Priority = as.character(priority),
      SourceComponent = as.character(source_component),
      Detail = as.character(detail),
      NextAction = as.character(if (isTRUE(available)) available_action else missing_action),
      stringsAsFactors = FALSE
    )
  }

  checklist <- do.call(
    rbind,
    list(
      add_item(
        "Method Section",
        "Model specification",
        !is.null(config$model) && !is.null(config$method),
        detail = sprintf(
          "Model=%s; Method=%s",
          as.character(config$model %||% "NA"),
          as.character(config$method %||% "NA")
        ),
        source_component = "fit$config",
        missing_action = "Fit a model first so the APA report can name the model and estimation method."
      ),
      add_item(
        "Method Section",
        "Data description",
        is.finite(as.numeric(prep$n_obs %||% NA)),
        detail = sprintf(
          "%s observations; %s persons; %s categories (%s-%s)",
          format(prep$n_obs %||% NA, big.mark = ","),
          format(prep$n_person %||% NA, big.mark = ","),
          as.character(config$n_cat %||% NA),
          as.character(prep$rating_min %||% NA),
          as.character(prep$rating_max %||% NA)
        ),
        source_component = "fit$prep + fit$config",
        missing_action = "Populate the basic design counts so the manuscript can describe the sample and scale."
      ),
      add_item(
        "Method Section",
        "Precision basis",
        nrow(precision_profile_df) > 0,
        detail = if (nrow(precision_profile_df) > 0) {
          paste0(
            "Precision tier = ", precision_tier,
            "; ", as.character(precision_profile_df$RecommendedUse[1] %||% "")
          )
        } else {
          "No precision-profile summary"
        },
        source_component = "diagnostics$precision_profile",
        ready_for_apa = nrow(precision_profile_df) > 0,
        missing_action = "Run diagnostics so the report can explain whether precision is model-based, hybrid, or exploratory.",
        available_action = if (formal_precision) {
          "Report the precision tier as model-based in the APA narrative."
        } else {
          "Report the precision tier explicitly and keep the exploratory/hybrid caution in the APA narrative."
        }
      ),
      add_item(
        "Method Section",
        "Convergence",
        converged,
        detail = if (nzchar(conv_msg)) conv_msg else if (converged) "Converged" else "Convergence status unavailable",
        source_component = "fit$summary + fit$opt",
        missing_action = "Resolve convergence before reporting model results."
      ),
      add_item(
        "Method Section",
        "Connectivity assessed",
        !is.null(diagnostics$subsets),
        detail = if (!is.null(diagnostics$subsets)) "Connectivity/subset output available" else "No subset output",
        source_component = "diagnostics$subsets",
        severity = "recommended",
        missing_action = "Run the subset/connectivity diagnostics and summarize whether the design is connected.",
        available_action = "Document the connectivity result before making common-scale or linking claims."
      ),
      add_item(
        "Global Fit",
        "Standardized residuals",
        has_resid,
        detail = if (has_resid) "Observation-level standardized residuals available" else "No standardized residuals",
        source_component = "diagnostics$obs",
        missing_action = "Compute diagnostics so global fit and local residual screening can be reported.",
        available_action = "Use standardized residuals as screening diagnostics, not as standalone proof of model adequacy."
      ),
      add_item(
        "Global Fit",
        "PCA of residuals",
        has_pca,
        detail = if (has_pca) "Residual PCA output available" else "Residual PCA not computed",
        source_component = "diagnostics$pca",
        severity = "recommended",
        missing_action = "Run residual PCA if you want to comment on unexplained residual structure.",
        available_action = "Report residual PCA as exploratory residual-structure follow-up, not as a standalone dimensionality test."
      ),
      add_item(
        "Facet-Level Statistics",
        "Separation / strata / reliability",
        has_rel,
        detail = if (has_rel) {
          if (formal_precision) {
            "Facet separation/reliability table available with model-based precision"
          } else {
            "Facet separation/reliability table available as an exploratory precision summary"
          }
        } else {
          "No reliability table"
        },
        source_component = "diagnostics$reliability",
        ready_for_apa = has_rel && formal_precision,
        missing_action = "Compute facet reliability/separation before describing facet spread and precision.",
        available_action = if (formal_precision) {
          "Report facet reliability/separation directly in the APA results section."
        } else {
          "Report facet reliability/separation as exploratory or hybrid, not as formal inferential evidence."
        }
      ),
      add_item(
        "Facet-Level Statistics",
        "Fixed/random variability summary",
        has_rel,
        detail = if (has_rel) {
          if (formal_precision) {
            "Facet-level variability summary available with fixed/random reference statistics"
          } else {
            "Facet-level variability summary available, but precision is not fully model-based"
          }
        } else {
          "No variability summary"
        },
        source_component = "diagnostics$reliability + diagnostics$facets_chisq",
        ready_for_apa = has_rel && formal_precision,
        missing_action = "Compute facet variability summaries before discussing fixed/random spread across facets.",
        available_action = if (formal_precision) {
          "Use the fixed/random variability summary in the results text or table notes."
        } else {
          "Describe the variability summary as exploratory or screening-oriented."
        }
      ),
      add_item(
        "Facet-Level Statistics",
        "RMSE and true SD",
        has_rel && any(c("RMSE", "TrueSD") %in% names(rel_df)),
        detail = if (has_rel && any(c("RMSE", "TrueSD") %in% names(rel_df))) "RMSE/True SD columns present" else "RMSE/True SD not found",
        source_component = "diagnostics$reliability",
        severity = "recommended",
        missing_action = "Expose RMSE/TrueSD columns if you want the manuscript to summarize facet spread in detail."
      ),
      add_item(
        "Element-Level Statistics",
        "Measures with SE",
        has_ci,
        detail = if (has_ci) {
          if (formal_precision) {
            "Element estimates with model-based SE / ModelSE are available"
          } else {
            "Element estimates with exploratory SE / ModelSE are available"
          }
        } else {
          "Estimate/SE columns not complete"
        },
        source_component = "diagnostics$measures",
        ready_for_apa = has_ci && measure_formal,
        missing_action = "Compute element-level measures and SE before reporting facet-level results.",
        available_action = if (formal_precision) {
          "Use ModelSE/RealSE language consistently in tables and notes."
        } else {
          "Use the available SE columns, but label them as exploratory or fit-adjusted summaries."
        }
      ),
      add_item(
        "Element-Level Statistics",
        "95% confidence intervals",
        has_ci,
        detail = if (has_ci) {
          if (formal_precision) {
            "Approximate normal CIs can be derived from Estimate +/- 1.96 * SE"
          } else {
            "Exploratory normal CIs can be derived from Estimate +/- 1.96 * SE"
          }
        } else {
          "SE not available"
        },
        source_component = "diagnostics$measures",
        severity = "recommended",
        ready_for_apa = ci_ready,
        missing_action = "Add SE first if you plan to report approximate confidence intervals.",
        available_action = if (formal_precision) {
          "Report approximate normal intervals if they are substantively useful."
        } else {
          "If reported, label the intervals as approximate and exploratory."
        }
      ),
      add_item(
        "Element-Level Statistics",
        "Infit and Outfit statistics",
        has_fit,
        detail = if (has_fit) "Infit/Outfit columns available" else "Fit statistics not complete",
        source_component = "diagnostics$measures or diagnostics$fit",
        missing_action = "Compute fit statistics before reporting misfit or element-level screening."
      ),
      add_item(
        "Element-Level Statistics",
        "Misfit counts",
        has_fit,
        detail = if (has_fit) "Misfit counts can be summarized from element fit" else "No fit table for misfit counting",
        source_component = "diagnostics$fit",
        severity = "recommended",
        missing_action = "Summarize element-level misfit counts if the manuscript reports quality-control findings."
      ),
      add_item(
        "Rating Scale Diagnostics",
        "Category counts",
        has_counts,
        detail = if (has_counts) "Observed score distribution available" else "Observed score counts unavailable",
        source_component = "diagnostics$obs",
        missing_action = "Expose observed score counts before commenting on category use.",
        available_action = "Use category counts as descriptive support for category use, not as a standalone quality verdict."
      ),
      add_item(
        "Rating Scale Diagnostics",
        "Average measures by category",
        has_person_measure,
        detail = if (has_person_measure) "Person measures can be summarized by observed category" else "Person measures unavailable",
        source_component = "diagnostics$obs",
        severity = "recommended",
        missing_action = "Retain person-measure information if the manuscript discusses average measures by category.",
        available_action = "Report average measures by category as descriptive scale-functioning evidence."
      ),
      add_item(
        "Rating Scale Diagnostics",
        "Threshold ordering",
        has_steps,
        detail = if (has_steps) "Step/threshold table available" else "No step table",
        source_component = "fit$steps",
        missing_action = "Fit step/threshold estimates before reporting category structure.",
        available_action = "Describe threshold ordering as category-structure evidence under the fitted model, not as a standalone proof of scale validity."
      ),
      add_item(
        "Rating Scale Diagnostics",
        "Category probability curves",
        has_steps,
        detail = if (has_steps) "Curve inputs available from threshold table" else "No step table for curves",
        source_component = "fit$steps",
        severity = "recommended",
        missing_action = "Retain the threshold table if you plan to include category probability curves.",
        available_action = "Use category probability curves as descriptive follow-up for scale structure."
      ),
      add_item(
        "Bias / Interaction Analysis",
        "Facet pairs tested",
        n_bias_pairs > 0 || n_bias_errors > 0,
        detail = if (n_bias_errors > 0) {
          sprintf("%d bias result bundle(s); %d requested pair(s) failed", n_bias_pairs, n_bias_errors)
        } else {
          sprintf("%d bias result bundle(s)", n_bias_pairs)
        },
        source_component = "bias_results",
        severity = "recommended",
        missing_action = "Run bias screening if the manuscript needs interaction-level follow-up."
      ),
      add_item(
        "Bias / Interaction Analysis",
        "Screen-positive interactions",
        n_bias_pairs > 0 || n_bias_errors > 0,
        detail = if (n_bias_errors > 0) {
          sprintf(
            "Bias screening completed for %d bundle(s), but %d requested pair(s) failed; review the bias collection errors before interpreting screen-positive counts.",
            n_bias_pairs,
            n_bias_errors
          )
        } else if (bias_stat_parse_issue) {
          "Bias screening output included non-numeric screening statistics; verify the bias tables before interpreting screen-positive counts."
        } else if (has_bias_sig) {
          "At least one interaction row exceeded the screening threshold."
        } else if (n_bias_pairs > 0) {
          "Bias screening was run and no interaction rows crossed the current screening threshold."
        } else {
          "No bias screening output available."
        },
        source_component = "bias_results$table",
        severity = "recommended",
        ready_for_apa = FALSE,
        missing_action = "Run bias screening before discussing interaction-level anomalies.",
        available_action = "Report these findings as screening results, not as formal hypothesis tests."
      ),
      add_item(
        "Visual Displays",
        "Wright map",
        TRUE,
        detail = "Supported by plot.mfrm_fit() / plot_wright_unified()",
        source_component = "plot.mfrm_fit",
        severity = "recommended",
        available_action = "Include a Wright map if the manuscript benefits from a scale-location display."
      ),
      add_item(
        "Visual Displays",
        "Fit/category visuals",
        has_fit || has_steps,
        detail = if (has_fit || has_steps) "Plotting inputs available for fit/category visuals" else "No fit or step visuals available",
        source_component = "diagnostics$fit + fit$steps",
        severity = "optional",
        missing_action = "Add fit or threshold inputs if you want figure-ready diagnostics."
      )
    )
  )

  section_summary <- checklist |>
    dplyr::group_by(.data$Section) |>
    dplyr::summarise(
      Items = dplyr::n(),
      Available = sum(.data$Available, na.rm = TRUE),
      DraftReady = sum(.data$DraftReady, na.rm = TRUE),
      ReadyForAPA = sum(.data$ReadyForAPA, na.rm = TRUE),
      Missing = .data$Items - .data$Available,
      NeedsDraftWork = .data$Items - .data$DraftReady,
      NeedsAction = .data$Items - .data$ReadyForAPA,
      .groups = "drop"
    )

  references <- if (isTRUE(include_references)) {
    data.frame(
      Citation = c(
        "Eckes (2005)",
        "Koizumi et al. (2019)",
        "Myford & Wolfe (2003, 2004)",
        "Linacre (1989, 2004)",
        "Wright & Masters (1982)"
      ),
      Topic = c(
        "Rater effects in MFRM",
        "Validity / MFRM task reporting",
        "Bias and interaction analysis",
        "MFRM and rating scale guidance",
        "Rating scale analysis"
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  settings <- list(
    include_references = isTRUE(include_references),
    diagnostics_supplied = !missing(diagnostics) && !is.null(diagnostics),
    bias_result_count = n_bias_pairs,
    bias_error_count = n_bias_errors,
    precision_tier = precision_tier
  )

  out <- list(
    checklist = checklist,
    summary = as.data.frame(section_summary, stringsAsFactors = FALSE),
    section_summary = as.data.frame(section_summary, stringsAsFactors = FALSE),
    references = references,
    settings = settings
  )
  as_mfrm_bundle(out, "mfrm_reporting_checklist")
}

collect_bias_tables_for_checklist <- function(bias_results) {
  if (is.null(bias_results)) return(list())
  error_tbl <- data.frame()
  if (inherits(bias_results, "mfrm_bias_collection")) {
    error_tbl <- as.data.frame(bias_results$errors %||% data.frame(), stringsAsFactors = FALSE)
    bias_results <- bias_results$by_pair %||% list()
  }
  if (inherits(bias_results, "mfrm_bias")) {
    out <- list(bias = as.data.frame(bias_results$table %||% data.frame(), stringsAsFactors = FALSE))
    attr(out, "errors") <- error_tbl
    return(out)
  }
  if (is.list(bias_results) && !is.data.frame(bias_results)) {
    out <- list()
    nms <- names(bias_results)
    if (is.null(nms)) {
      nms <- paste0("bias_", seq_along(bias_results))
    }
    for (i in seq_along(bias_results)) {
      obj <- bias_results[[i]]
      if (inherits(obj, "mfrm_bias")) {
        out[[nms[i]]] <- as.data.frame(obj$table %||% data.frame(), stringsAsFactors = FALSE)
      }
    }
    attr(out, "errors") <- error_tbl
    return(out)
  }
  out <- list()
  attr(out, "errors") <- error_tbl
  out
}
