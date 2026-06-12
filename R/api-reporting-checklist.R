#' Build an auto-filled MFRM reporting checklist
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are computed with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()] or a named list of
#'   such outputs.
#' @param hierarchical_structure Optional output from
#'   [analyze_hierarchical_structure()]. When supplied, the
#'   "Hierarchical structure review" checklist item is flipped to
#'   `DraftReady = TRUE` and its `Detail` column surfaces the number
#'   of nested / crossed facet pairs and whether the ICC table is
#'   available.
#' @param include_references If `TRUE`, include a compact reference table in the
#'   returned bundle.
#'
#' @details
#' This helper ports the app-level reporting checklist into a package-native
#' bundle. It does not try to judge substantive reporting quality; instead, it
#' checks whether the fitted object and related diagnostics contain the evidence
#' typically reported in MFRM write-ups.
#'
#' Checklist items are grouped into seven core sections:
#' - Method section
#' - Global fit
#' - Facet-level statistics
#' - Element-level statistics
#' - Rating scale diagnostics
#' - Bias/interaction analysis
#' - Visual displays
#'
#' When a fit uses the latent-regression population-model branch, the checklist
#' also adds a `Population Model` section covering coefficient reporting,
#' categorical model-matrix coding, complete-case omissions, posterior-basis
#' wording, and ConQuest scope wording.
#'
#' The output is designed for manuscript preparation, reproducibility records, and
#' reproducible reporting workflows.
#'
#' @section What this checklist means:
#' `reporting_checklist()` is a manuscript-preparation guide. It tells you
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
#' - The Global Fit section includes a "Fit/separation reporting boundary"
#'   row that points to [precision_review_report()], [fit_measures_table()],
#'   and [facets_fit_review()] before users phrase fit, ZSTD, separation, or
#'   reliability claims.
#' - `software_scope`: external-software relationship summary for `mfrmr`,
#'   FACETS, ConQuest, and SPSS-style tabular handoffs.
#' - `facets_positioning`: report-ready wording that states `mfrmr` is not a
#'   FACETS numerical clone and separates native estimation from FACETS-style
#'   handoff or external-table review.
#' - `visual_scope`: plotting-route summary that separates report-default
#'   2D figures from exploratory surface/3D-ready data handoffs, including a
#'   short `InterpretationCheck` for the main user-facing caveat.
#' - `references`: core background references when requested.
#'
#' @section Recommended next step:
#' Review the rows with `Available = FALSE` or `DraftReady = FALSE`, then add
#' the missing diagnostics, bias results, or narrative context before calling
#' [build_apa_outputs()] for draft text generation. For `RSM` / `PCM`
#' reporting runs, the preferred route is an `MML` fit plus
#' `diagnose_mfrm(..., diagnostic_mode = "both")` so the checklist can see the
#' legacy and strict marginal screens together.
#'
#' @section How this differs from operational review:
#' `reporting_checklist()` is the manuscript/reporting branch of the package.
#' Use it when the question is "what is still missing from the report?" rather
#' than "which observations or links need follow-up?" For operational review:
#' - Use [build_misfit_casebook()] after [diagnose_mfrm()] when you need ranked
#'   misfit cases and grouping views for local follow-up.
#' - Use [build_linking_review()] after anchor/drift/chain helpers when you
#'   need operational linking triage rather than manuscript-oriented reporting
#'   tables.
#'
#' @section Typical workflow:
#' 1. Fit with [fit_mfrm()]. For `RSM` / `PCM` reporting runs, prefer
#'    `method = "MML"`.
#' 2. Compute diagnostics with [diagnose_mfrm()]. For `RSM` / `PCM`, prefer
#'    `diagnostic_mode = "both"`.
#' 3. Run `reporting_checklist()` to see which reporting elements are already
#'    available from the current analysis objects.
#' 4. If the issue is operational rather than manuscript-facing, branch to
#'    [build_misfit_casebook()] or [build_linking_review()] instead of treating
#'    `reporting_checklist()` as the single review hub.
#'
#' @return A named list with checklist tables. Class:
#'   `mfrm_reporting_checklist`.
#' @seealso [build_apa_outputs()], [build_visual_summaries()],
#'   [specifications_report()], [data_quality_report()],
#'   [build_misfit_casebook()], [build_linking_review()]
#' @examples
#' # Fast smoke run: a JML fit + legacy-only diagnostic produces a
#' # populated checklist in well under a second.
#' toy <- load_mfrmr_data("example_core")
#' fit_quick <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                       method = "JML", maxit = 30)
#' diag_quick <- diagnose_mfrm(fit_quick, residual_pca = "none",
#'                              diagnostic_mode = "legacy")
#' chk_quick <- reporting_checklist(fit_quick, diagnostics = diag_quick)
#' head(chk_quick$checklist[, c("Section", "Item", "DraftReady")])
#'
#' \donttest{
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 7, maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "both", diagnostic_mode = "both")
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' summary(chk)
#' # Look for: a high `Ready` / `Total` ratio in the summary block.
#' #   Sections with `Ready = 0` need follow-up before submitting
#' #   (typically diagnostic_mode = "both" or a residual-PCA pass).
#' apa <- build_apa_outputs(fit, diag)
#' head(chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#' # Look for: every row where `DraftReady = "yes"` is ready to paste
#' #   into the manuscript. `"no"` rows include a concrete `NextAction`
#' #   step (e.g. "run plot_qc_dashboard()") so the gap can be closed
#' #   without re-reading the methodology guide.
#' nchar(apa$report_text)
#' }
#' @export
reporting_checklist <- function(fit,
                                diagnostics = NULL,
                                bias_results = NULL,
                                hierarchical_structure = NULL,
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
  if (!is.null(hierarchical_structure) &&
      !inherits(hierarchical_structure, "mfrm_hierarchical_structure")) {
    stop("`hierarchical_structure` must come from ",
         "analyze_hierarchical_structure().", call. = FALSE)
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
  has_subsets <- !is.null(diagnostics$subsets)
  marginal_fit_bundle <- diagnostics$marginal_fit %||% list()
  strict_marginal_available <- isTRUE(marginal_fit_bundle$available)
  strict_pairwise_available <- isTRUE(marginal_fit_bundle$pairwise$available)
  strict_marginal_reason <- if (strict_marginal_available) {
    if (strict_pairwise_available) {
      "Strict marginal first-order and pairwise bundles are available for plotting."
    } else {
      "Strict marginal first-order bundle is available, but pairwise follow-up is unavailable."
    }
  } else {
    as.character(
      marginal_fit_bundle$summary$Reason[1] %||%
        marginal_fit_bundle$notes[1] %||%
        "Strict marginal diagnostics are unavailable for the current run."
    )
  }
  fit_df_review_available <- nrow(measures) > 0 &&
    all(c("DF_Infit_FACETS", "DF_Outfit_FACETS",
          "InfitZSTD_FACETS", "OutfitZSTD_FACETS") %in% names(measures))
  fit_separation_boundary_available <- has_fit || has_rel || fit_df_review_available
  fit_separation_boundary_detail <- paste0(
    "MnSq fit=", if (has_fit) "available" else "missing",
    "; separation/reliability=", if (has_rel) "available" else "missing",
    "; df/ZSTD review=", if (fit_df_review_available) "available" else "not requested"
  )
  population <- fit$population %||% list()
  population_active <- isTRUE(population$active)
  population_coefficients <- as.numeric(population$coefficients %||% numeric(0))
  population_coefficients_finite <- length(population_coefficients) > 0L &&
    all(is.finite(population_coefficients))
  population_sigma2 <- suppressWarnings(as.numeric(population$sigma2 %||% NA_real_))
  population_formula <- if (!is.null(population$formula)) {
    paste(deparse(population$formula), collapse = " ")
  } else {
    ""
  }
  population_posterior_basis <- as.character(population$posterior_basis %||% "legacy_mml")
  population_coding <- population_coding_summary_table(population)
  population_design_columns <- as.character(population$design_columns %||% character(0))
  population_omitted_persons <- length(population$omitted_persons %||% character(0))
  population_omitted_rows <- suppressWarnings(as.integer(population$response_rows_omitted %||% 0L))
  if (!is.finite(population_omitted_rows)) population_omitted_rows <- 0L

  add_item <- function(section,
                       item,
                       available,
                       detail,
                       source_component,
                       severity = c("required", "recommended", "optional"),
                       ready_for_apa = available,
                       missing_action = "Compute or document this component before manuscript export.",
                       available_action = NULL,
                       plot_helper = NA_character_,
                       draw_free_route = NA_character_,
                       plot_return_class = NA_character_) {
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
      PlotHelper = as.character(plot_helper),
      DrawFreeRoute = as.character(draw_free_route),
      PlotReturnClass = as.character(plot_return_class),
      NextAction = as.character(if (isTRUE(available)) available_action else missing_action),
      stringsAsFactors = FALSE
    )
  }

  population_items <- list()
  if (population_active) {
    coding_detail <- if (nrow(population_coding) > 0L) {
      paste0(
        "Stored coding for ",
        paste(population_coding$Variable, collapse = ", "),
        "; encoded columns: ",
        paste(population_coding$EncodedColumns[nzchar(population_coding$EncodedColumns)], collapse = ", ")
      )
    } else {
      "No categorical coding was recorded; formula appears intercept-only, numeric, or logical."
    }
    if (!nzchar(coding_detail)) {
      coding_detail <- "No encoded categorical columns recorded."
    }
    omission_detail <- paste0(
      "Omitted persons = ", population_omitted_persons,
      "; omitted response rows = ", population_omitted_rows,
      "; policy = ", as.character(population$policy %||% NA_character_)
    )
    population_items <- list(
      add_item(
        "Population Model",
        "Latent-regression basis",
        population_active && identical(population_posterior_basis, "population_model"),
        detail = paste0(
          "Formula=", population_formula,
          "; posterior basis=", population_posterior_basis,
          "; model=", as.character(config$model %||% NA_character_),
          "; method=", as.character(config$method %||% NA_character_)
        ),
        source_component = "fit$population",
        severity = "required",
        ready_for_apa = population_active && identical(population_posterior_basis, "population_model"),
        missing_action = "Fit with `method = \"MML\"`, `population_formula`, and one-row-per-person `person_data` before reporting latent regression.",
        available_action = "Describe the fit as a first-version conditional-normal latent-regression MML model, not as a post hoc regression on EAP/MLE scores."
      ),
      add_item(
        "Population Model",
        "Population coefficients and residual variance",
        population_coefficients_finite && is.finite(population_sigma2),
        detail = paste0(
          length(population_coefficients), " coefficient(s); residual variance = ",
          if (is.finite(population_sigma2)) signif(population_sigma2, 4) else NA_real_
        ),
        source_component = "summary(fit)$population_coefficients + summary(fit)$population_overview",
        severity = "required",
        ready_for_apa = population_coefficients_finite && is.finite(population_sigma2) && population_sigma2 > 0,
        missing_action = "Inspect `summary(fit)$population_coefficients` and `summary(fit)$population_overview` before reporting latent-regression effects.",
        available_action = "Report coefficients and residual variance as conditional-normal population-model parameters, with scale/coding notes."
      ),
      add_item(
        "Population Model",
        "Model-matrix covariate coding",
        length(population_design_columns) > 0L,
        detail = coding_detail,
        source_component = "summary(fit)$population_coding + fit$population$xlevels + fit$population$contrasts",
        severity = "required",
        ready_for_apa = length(population_design_columns) > 0L,
        missing_action = "Inspect `summary(fit)$population_coding` and document categorical levels/contrasts before scoring or reporting the population model.",
        available_action = "Use `summary(fit)$population_coding` to document categorical levels, contrasts, and encoded columns used by scoring/replay."
      ),
      add_item(
        "Population Model",
        "Complete-case omission review",
        TRUE,
        detail = omission_detail,
        source_component = "summary(fit)$population_overview + summary(fit)$caveats",
        severity = "required",
        ready_for_apa = population_omitted_persons == 0L && population_omitted_rows == 0L,
        missing_action = "Review omitted persons/rows and the population-data policy before reporting the population model.",
        available_action = if (population_omitted_persons == 0L && population_omitted_rows == 0L) {
          "State that no persons were omitted by the population-model covariate policy."
        } else {
          "Document the complete-case policy and omitted-person/row counts before reporting latent-regression results."
        }
      ),
      add_item(
        "Population Model",
        "Population-model posterior scoring wording",
        TRUE,
        detail = "Active latent-regression scoring should condition on the fitted population model; new-person scoring requires matching `person_data` when covariates are present.",
        source_component = "predict_mfrm_units() / sample_mfrm_plausible_values()",
        severity = "recommended",
        ready_for_apa = TRUE,
        available_action = "When reporting EAP/PV outputs, state whether they use the fitted population-model posterior and document the required `person_data` input."
      ),
      add_item(
        "Population Model",
        "ConQuest overlap wording",
        TRUE,
        detail = "Current overlap is narrow RSM/PCM unidimensional conditional-normal latent regression; ConQuest comparison is scoped to the documented external-table workflow.",
        source_component = "README latent-regression status + review_conquest_overlap()",
        severity = "recommended",
        ready_for_apa = FALSE,
        available_action = "Use conservative wording: ConQuest overlap is limited to the documented latent-regression MML comparison scope."
      )
    )
  }

  checklist <- do.call(
    rbind,
    c(
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
        "Method Section",
        "Empirical-Bayes shrinkage when small-N facets are present",
        {
          shrink_mode <- as.character(config$facet_shrinkage %||% "none")
          sparse_n <- suppressWarnings(as.integer(
            fit$summary$FacetSparseCount %||% NA_integer_
          ))
          # Ready if either (a) shrinkage was applied, or (b) there are no
          # sparse facets so shrinkage isn't needed.
          (!identical(shrink_mode, "none")) ||
            (is.finite(sparse_n) && sparse_n == 0L)
        },
        detail = {
          shrink_mode <- as.character(config$facet_shrinkage %||% "none")
          sparse_n <- suppressWarnings(as.integer(
            fit$summary$FacetSparseCount %||% NA_integer_
          ))
          if (!identical(shrink_mode, "none")) {
            paste0(
              "Shrinkage active: ", shrink_mode, ". See `fit$shrinkage_report`."
            )
          } else if (is.finite(sparse_n) && sparse_n == 0L) {
            "No sparse facets detected; fixed-effects estimates are stable without shrinkage."
          } else {
            paste0(
              "Sparse facet(s) detected (count = ", sparse_n,
              ") but no shrinkage was applied. Consider ",
              "`fit_mfrm(..., facet_shrinkage = 'empirical_bayes')`."
            )
          }
        },
        source_component = "fit$config$facet_shrinkage + fit$shrinkage_report",
        severity = "recommended",
        missing_action = paste0(
          "Re-run with `facet_shrinkage = 'empirical_bayes'` or apply ",
          "`apply_empirical_bayes_shrinkage(fit)` post-hoc when small-N ",
          "facets are present."
        ),
        available_action = paste0(
          "Report both the fixed-effects and shrunk estimates; cite ",
          "Efron & Morris (1973) for the empirical-Bayes rationale."
        )
      ),
      add_item(
        "Method Section",
        "Facet sample-size adequacy",
        {
          flag <- suppressWarnings(as.character(fit$summary$FacetSampleSizeFlag %||% NA_character_))
          !is.na(flag) && flag %in% c("standard", "strong")
        },
        detail = {
          flag <- suppressWarnings(as.character(fit$summary$FacetSampleSizeFlag %||% NA_character_))
          min_n <- suppressWarnings(as.integer(fit$summary$FacetMinLevelN %||% NA_integer_))
          sparse_n <- suppressWarnings(as.integer(fit$summary$FacetSparseCount %||% NA_integer_))
          if (is.na(flag)) {
            "Sample-size flag unavailable in fit summary."
          } else {
            sprintf(
              "Worst facet band: %s (min level N = %s; sparse facets = %s).",
              flag,
              if (is.na(min_n)) NA_character_ else as.character(min_n),
              if (is.na(sparse_n)) NA_character_ else as.character(sparse_n)
            )
          }
        },
        source_component = "fit$summary$FacetSampleSizeFlag + facet_small_sample_review()",
        severity = "recommended",
        missing_action = paste0(
          "Run `facet_small_sample_review(fit)` and report the bands. ",
          "mfrmr treats facets as fixed effects with no shrinkage, so small-N ",
          "levels keep wide SEs."
        ),
        available_action = paste0(
          "Report the per-facet adequacy bands and discuss any sparse/marginal ",
          "levels; cite Linacre (1994) sample-size guidance where relevant."
        )
      ),
      add_item(
        "Method Section",
        "Hierarchical structure review",
        # Ready when the user actually ran analyze_hierarchical_structure()
        # and passed the result in. Previously this item was hard-coded to
        # FALSE so there was no way to mark it draft-ready; now callers can
        # supply `hierarchical_structure = hs`.
        !is.null(hierarchical_structure),
        detail = if (!is.null(hierarchical_structure)) {
          hs_sum <- hierarchical_structure$summary
          paste0(
            "Hierarchical review complete: ",
            hs_sum$NFacets %||% NA, " facets, ",
            hs_sum$NestedPairs %||% 0L, " nested pair(s), ",
            hs_sum$CrossedPairs %||% 0L, " crossed pair(s)",
            if (isTRUE(hs_sum$ICCAvailable)) "; ICC available" else "",
            "."
          )
        } else {
          paste0(
            "Structural review (nesting, ICC, design effect) is optional but ",
            "recommended when raters, criteria, or persons span strata ",
            "(regions, schools, cohorts) that additive fixed-effects MFRM cannot ",
            "partition out."
          )
        },
        source_component = "analyze_hierarchical_structure(data, facets)",
        severity = "recommended",
        missing_action = paste0(
          "Run `analyze_hierarchical_structure(fit)` once per design and pass ",
          "the result to `reporting_checklist(..., hierarchical_structure = hs)`."
        ),
        available_action = paste0(
          "Report nesting classifications and (where applicable) the Kish design ",
          "effect before generalizing rater severity beyond the sampled raters."
        )
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
        "Global Fit",
        "Fit/separation reporting boundary",
        fit_separation_boundary_available,
        detail = fit_separation_boundary_detail,
        source_component = "precision_review_report()$fit_separation_basis + fit_measures_table() + facets_fit_review()",
        severity = "recommended",
        ready_for_apa = fit_separation_boundary_available,
        missing_action = paste0(
          "Run `diagnose_mfrm()` before writing fit, separation, reliability, ",
          "or ZSTD language. Use `diagnose_mfrm(..., fit_df_method = \"both\")` ",
          "when FACETS-style df/ZSTD review is needed."
        ),
        available_action = paste0(
          "Use `precision_review_report(fit, diagnostics)$fit_separation_basis` ",
          "before drafting fit/separation wording; export `fit_measures_table()` ",
          "or `facets_fit_review()` when df/ZSTD sensitivity or external FACETS ",
          "matching is part of the report."
        )
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
        plot_helper = "plot.mfrm_fit() / plot_wright_unified()",
        draw_free_route = "plot(fit, type = \"wright\", draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        available_action = "Include a Wright map when the manuscript benefits from a shared-scale targeting display."
      ),
      add_item(
        "Visual Displays",
        "QC / facet dashboard",
        nrow(obs_df) > 0,
        detail = if (nrow(obs_df) > 0) {
          "plot_qc_dashboard() / plot_facet_quality_dashboard() can use the current diagnostics bundle"
        } else {
          "No observation-level diagnostics for dashboard plotting"
        },
        source_component = "diagnostics$obs + diagnostics$fit",
        severity = "recommended",
        plot_helper = "plot_qc_dashboard() / plot_facet_quality_dashboard()",
        draw_free_route = "plot_qc_dashboard(fit, diagnostics = diagnostics, draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Run diagnose_mfrm() so the QC and facet dashboards can be rendered from the current analysis bundle.",
        available_action = "Use the dashboard as a first-pass triage view, then move to the specific follow-up plot behind each flag."
      ),
      add_item(
        "Visual Displays",
        "Residual PCA visuals",
        has_pca,
        detail = if (has_pca) {
          "plot_residual_pca() can render scree/loadings from the current residual PCA bundle"
        } else {
          "Residual PCA not computed"
        },
        source_component = "diagnostics$pca",
        severity = "recommended",
        plot_helper = "plot_residual_pca()",
        draw_free_route = "plot_residual_pca(analyze_residual_pca(diagnostics, mode = \"overall\"), mode = \"overall\", plot_type = \"scree\", draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Run residual PCA if you want scree/loadings visuals for residual-structure follow-up.",
        available_action = "Use residual PCA plots as exploratory structure follow-up, not as standalone dimensionality proof."
      ),
      add_item(
        "Visual Displays",
        "Connectivity / design-matrix visual",
        has_subsets,
        detail = if (has_subsets) {
          "subset_connectivity_report() and plot(..., type = \"design_matrix\") can use the current subset bundle"
        } else {
          "No subset/connectivity bundle available"
        },
        source_component = "diagnostics$subsets",
        severity = "recommended",
        plot_helper = "plot.mfrm_subset_connectivity()",
        draw_free_route = "plot(subset_connectivity_report(fit, diagnostics = diagnostics), type = \"design_matrix\", draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Run subset/connectivity diagnostics before showing design-matrix or linkage visuals.",
        available_action = "Use the design-matrix view to support linkage and comparability claims."
      ),
      add_item(
        "Visual Displays",
        "Inter-rater / displacement visuals",
        has_fit || has_resid,
        detail = if (has_fit || has_resid) {
          "plot_displacement() is available; plot_interrater_agreement() is available when a rater facet is present"
        } else {
          "No fit/residual inputs for displacement or inter-rater visuals"
        },
        source_component = "diagnostics$obs + diagnostics$measures",
        severity = "recommended",
        plot_helper = "plot_displacement() / plot_interrater_agreement()",
        draw_free_route = "plot_displacement(displacement_table(fit, diagnostics = diagnostics), draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Retain fit and residual outputs if you want displacement or inter-rater follow-up figures.",
        available_action = "Use displacement and inter-rater views to localize QC issues after dashboard screening."
      ),
      add_item(
        "Visual Displays",
        "Strict marginal visuals",
        strict_marginal_available,
        detail = strict_marginal_reason,
        source_component = "diagnostics$marginal_fit",
        severity = "recommended",
        ready_for_apa = FALSE,
        plot_helper = "plot_marginal_fit() / plot_marginal_pairwise()",
        draw_free_route = "plot_marginal_fit(diagnostics, draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "For MML reporting runs, call diagnose_mfrm(..., diagnostic_mode = \"both\") to enable strict marginal follow-up visuals where supported.",
        available_action = "Treat strict marginal plots as exploratory corroboration screens, then corroborate with design review and legacy diagnostics."
      ),
      add_item(
        "Visual Displays",
        "Bias / DIF visuals",
        n_bias_pairs > 0 || n_bias_errors > 0,
        detail = if (n_bias_errors > 0) {
          "plot_bias_interaction() can use the current bias bundle, but failed pair requests should be reviewed before interpretation"
        } else if (n_bias_pairs > 0) {
          "plot_bias_interaction() can use the current bias bundle; use plot_dif_heatmap() when a DIF result is available"
        } else {
          "No bias screening bundle available"
        },
        source_component = "bias_results",
        severity = "recommended",
        ready_for_apa = FALSE,
        plot_helper = "plot_bias_interaction() / plot_dif_heatmap()",
        draw_free_route = "plot_bias_interaction(bias_results[[1]], draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Run bias or DIF screening before discussing interaction-level visuals.",
        available_action = "Use bias/DIF plots as screening follow-up, not as formal hypothesis tests."
      ),
      add_item(
        "Visual Displays",
        "Precision / information curves",
        converged,
        detail = if (converged) {
          "compute_information() / plot_information() are available for the current fitted model"
        } else {
          "Information curves are not recommended until convergence is resolved"
        },
        source_component = "fit + compute_information()",
        severity = "recommended",
        ready_for_apa = converged && formal_precision,
        plot_helper = "plot_information()",
        draw_free_route = "plot_information(compute_information(fit), draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Resolve convergence before using information or precision curves in reporting.",
        available_action = if (formal_precision) {
          "Use information curves to describe precision across theta when that is the reporting question."
        } else {
          "Use information curves descriptively and keep the current precision-tier caveat in the narrative."
        }
      ),
      add_item(
        "Visual Displays",
        "Fit/category visuals",
        has_fit || has_steps,
        detail = if (has_fit || has_steps) "Plotting inputs available for fit/category visuals" else "No fit or step visuals available",
        source_component = "diagnostics$fit + fit$steps",
        severity = "optional",
        plot_helper = "plot.mfrm_fit() / plot.mfrm_category_curves()",
        draw_free_route = "plot(fit, type = \"ccc\", draw = FALSE)",
        plot_return_class = "mfrm_plot_data",
        missing_action = "Add fit or threshold inputs if you want figure-ready diagnostics.",
        available_action = "Use category curves and fit visuals as local descriptive follow-up after QC screening."
      )
      ),
      population_items
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
        "Linacre (1989, 2002)",
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
    software_scope = external_software_scope_table(fit),
    facets_positioning = facets_positioning_guide(),
    visual_scope = visual_scope_table(fit, checklist),
    references = references,
    settings = settings
  )
  as_mfrm_bundle(out, "mfrm_reporting_checklist")
}

external_software_scope_table <- function(fit) {
  cfg <- fit$config %||% list()
  population <- fit$population %||% list()
  model <- toupper(as.character(cfg$model %||% NA_character_))[1]
  method <- toupper(as.character(cfg$method %||% cfg$method_input %||% NA_character_))
  population_active <- isTRUE(population$active)
  conquest_candidate <- population_active &&
    identical(method, "MML") &&
    model %in% c("RSM", "PCM")

  data.frame(
    Software = c("mfrmr native", "FACETS", "ConQuest", "SPSS"),
    Relationship = c(
      "primary estimation/reporting surface",
      "FACETS-style reporting and handoff surface",
      "scoped external-table comparison for latent-regression overlap",
      "downstream table/report handoff only"
    ),
    CurrentSupport = c(
      "active",
      "active compatibility layer",
      if (conquest_candidate) "candidate fit; additional exact-overlap restrictions apply" else "not active for this fit",
      "CSV/data.frame outputs only"
    ),
    PrimaryHelpers = c(
      "fit_mfrm() -> diagnose_mfrm() -> reporting_checklist() -> build_apa_outputs()",
      "facets_positioning_guide(), facets_feature_coverage(), run_mfrm_facets(), facets_output_file_bundle(), facets_output_contract_review()",
      "build_conquest_overlap_bundle() -> normalize_conquest_overlap_*() -> review_conquest_overlap()",
      "export_mfrm_bundle(), export_summary_appendix(), as.data.frame()"
    ),
    Boundary = c(
      "Package-native results are the authoritative analysis objects.",
      "Results remain mfrmr estimates; use an external FACETS run plus review helpers only when numerical comparison is needed.",
      "Requires an external ConQuest run and extracted output tables for the documented overlap case.",
      "CSV/data-frame outputs support reporting handoff; native SPSS integration is not implemented."
    ),
    RecommendedWording = c(
      "Estimated with mfrmr under the stated model/method settings.",
      "Estimated with mfrmr; FACETS-style outputs were used for handoff or report organization unless external FACETS output is explicitly compared.",
      "ConQuest overlap is limited to the documented latent-regression MML comparison scope.",
      "Tables were exported for possible SPSS/reporting use; analysis was not performed in SPSS."
    ),
    stringsAsFactors = FALSE
  )
}

visual_scope_table <- function(fit, checklist) {
  cfg <- fit$config %||% list()
  model <- toupper(as.character(cfg$model %||% NA_character_))[1]
  steps <- as.data.frame(fit$steps %||% data.frame(), stringsAsFactors = FALSE)
  has_steps <- nrow(steps) > 0 && "Estimate" %in% names(steps)

  visual_status <- function(item) {
    row <- checklist[checklist$Section == "Visual Displays" & checklist$Item == item, , drop = FALSE]
    if (nrow(row) == 0) return("not listed in reporting checklist")
    if (isTRUE(row$Available[1])) return("available for current fit")
    paste0("not ready for current fit: ", as.character(row$NextAction[1]))
  }

  surface_status <- if (has_steps && model %in% c("RSM", "PCM", "GPCM")) {
    "active plot-data route for current fit"
  } else if (!(model %in% c("RSM", "PCM", "GPCM"))) {
    paste0("not active for current model: ", model)
  } else {
    "not ready for current fit: step/threshold estimates are required"
  }

  data.frame(
    Visualization = c(
      "QC dashboard",
      "Wright map",
      "Pathway / CCC",
      "Category probability surface",
      "Information curves",
      "Strict marginal visuals",
      "Bias / DIF visuals",
      "Residual PCA visuals"
    ),
    Role = c(
      "first-pass diagnostic triage",
      "shared-scale targeting and spread",
      "ordered-category structure and category dominance",
      "exploratory theta x category x probability review",
      "precision and targeting across theta",
      "strict marginal fit follow-up",
      "interaction and group-functioning screening",
      "exploratory residual-structure follow-up"
    ),
    CurrentSupport = c(
      visual_status("QC / facet dashboard"),
      visual_status("Wright map"),
      visual_status("Fit/category visuals"),
      surface_status,
      visual_status("Precision / information curves"),
      visual_status("Strict marginal visuals"),
      visual_status("Bias / DIF visuals"),
      visual_status("Residual PCA visuals")
    ),
    SupportedModels = c(
      "RSM/PCM; GPCM residual stack with documented fair-average boundary",
      "RSM/PCM/GPCM when person, facet, and step locations are available",
      "RSM/PCM/GPCM when step/threshold estimates are available",
      "RSM/PCM/GPCM when step/threshold estimates are available",
      "RSM/PCM/GPCM where compute_information() supports the fitted object",
      "RSM/PCM MML diagnostic_mode = \"both\" where strict marginal bundles are available",
      "Model-dependent; requires supplied bias/DIF screening results",
      "Model-dependent; requires residual PCA computation"
    ),
    PrimaryHelper = c(
      "plot_qc_dashboard() / plot_facet_quality_dashboard()",
      "plot(fit, type = \"wright\") / plot_wright_unified()",
      "plot(fit, type = \"pathway\") / plot(fit, type = \"ccc\")",
      "plot(fit, type = \"ccc_surface\")",
      "compute_information() -> plot_information()",
      "plot_marginal_fit() / plot_marginal_pairwise()",
      "plot_bias_interaction() / plot_dif_heatmap()",
      "analyze_residual_pca() -> plot_residual_pca()"
    ),
    DrawFreeRoute = c(
      "plot_qc_dashboard(fit, diagnostics = diagnostics, draw = FALSE)",
      "plot(fit, type = \"wright\", draw = FALSE)",
      "plot(fit, type = \"ccc\", draw = FALSE)",
      "plot(fit, type = \"ccc_surface\", draw = FALSE)",
      "plot_information(compute_information(fit), draw = FALSE)",
      "plot_marginal_fit(diagnostics, draw = FALSE)",
      "plot_bias_interaction(bias_results[[1]], draw = FALSE)",
      "plot_residual_pca(analyze_residual_pca(diagnostics, mode = \"overall\"), draw = FALSE)"
    ),
    ReportUse = c(
      "triage figure; usually not the final manuscript figure by itself",
      "good candidate for reports when targeting/spread is the point",
      "good candidate for reports as descriptive category-functioning evidence",
      "exploratory or teaching appendix only; not a default APA figure",
      "report when precision/targeting is a substantive question and precision tier supports the wording",
      "exploratory follow-up, not a standalone inferential test",
      "screening follow-up, not a formal hypothesis-test figure by itself",
      "exploratory follow-up, not a standalone dimensionality test"
    ),
    ThreeDStatus = c(
      "2D dashboard only; 3D not recommended",
      "2D recommended; 3D Wright maps are discouraged",
      "2D report default; surface plot data available through the category route",
      "advanced surface data only; no package-native interactive renderer",
      "2D curve route active; 3D information surface is a future data handoff candidate",
      "2D heatmap/bar style preferred; 3D not recommended",
      "2D heatmap/profile preferred; 3D not recommended",
      "2D scree/loadings preferred; 3D not recommended"
    ),
    Boundary = c(
      "Summarizes flags; inspect the component plots before making claims.",
      "Shows common-scale locations but does not prove model fit.",
      "Descriptive category-functioning evidence, not proof of rating-scale validity.",
      "Returned data are exploratory mfrmr output for advanced visualization.",
      "Depends on the fitted model and precision tier; interpret with checklist caveats.",
      "Requires strict marginal diagnostics; use as corroborating evidence.",
      "Requires a bias/DIF result object; treat as screening unless separately justified.",
      "Residual PCA is exploratory residual-structure review, not a formal dimensionality decision."
    ),
    InterpretationCheck = c(
      "Review component warnings before citing the dashboard.",
      "Report targeting/spread, then check fit diagnostics separately.",
      "Check ordered peaks and dominance bands in the 2D CCC/pathway views.",
      "Inspect category_support before using the surface in an appendix or downstream renderer.",
      "State whether the curve is model-based, design-weighted, or approximate.",
      "Confirm diagnostic_mode = \"both\" and describe this as follow-up evidence.",
      "Confirm the tested pair, low-count cells, and screening threshold.",
      "Use scree/loadings as exploratory residual-structure evidence only."
    ),
    stringsAsFactors = FALSE
  )
}

#' Summarize a reporting-checklist bundle for manuscript work
#'
#' @param object Output from [reporting_checklist()].
#' @param top_n Maximum number of draft-action rows shown in the compact action
#'   table.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_reporting_checklist` with:
#' - `overview`: run-level counts of available and draft-ready items
#' - `section_summary`: section-level checklist coverage
#' - `software_scope`: external-software relationship summary
#' - `facets_positioning`: report-ready FACETS relationship wording
#' - `visual_scope`: plotting-route and 3D-ready data-handoff summary, including
#'   the main `InterpretationCheck` caveat for each visual family
#' - `priority_summary`: counts by priority/severity
#' - `action_items`: highest-priority rows that still need draft work
#' - `settings`: checklist settings rendered as a compact table
#' - `notes`: interpretation notes
#' @seealso [reporting_checklist()], [summary.mfrm_apa_outputs]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 7, maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "both", diagnostic_mode = "both")
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' summary(chk)
#' }
#' @method summary mfrm_reporting_checklist
#' @export
summary.mfrm_reporting_checklist <- function(object, top_n = 10, ...) {
  if (!inherits(object, "mfrm_reporting_checklist")) {
    stop("`object` must be output from reporting_checklist().", call. = FALSE)
  }

  top_n <- max(1L, as.integer(top_n[1]))
  checklist <- as.data.frame(object$checklist %||% data.frame(), stringsAsFactors = FALSE)
  section_summary <- as.data.frame(object$section_summary %||% object$summary %||% data.frame(), stringsAsFactors = FALSE)
  software_scope <- as.data.frame(object$software_scope %||% data.frame(), stringsAsFactors = FALSE)
  facets_positioning <- as.data.frame(object$facets_positioning %||% facets_positioning_guide(), stringsAsFactors = FALSE)
  visual_scope <- as.data.frame(object$visual_scope %||% data.frame(), stringsAsFactors = FALSE)

  overview <- data.frame(
    Sections = if (nrow(checklist) > 0) length(unique(checklist$Section)) else 0L,
    Items = nrow(checklist),
    Available = if (nrow(checklist) > 0) sum(checklist$Available %in% TRUE, na.rm = TRUE) else 0L,
    DraftReady = if (nrow(checklist) > 0) sum(checklist$DraftReady %in% TRUE, na.rm = TRUE) else 0L,
    Missing = if (nrow(checklist) > 0) sum(checklist$Available %in% FALSE, na.rm = TRUE) else 0L,
    NeedsDraftWork = if (nrow(checklist) > 0) sum(checklist$DraftReady %in% FALSE, na.rm = TRUE) else 0L,
    stringsAsFactors = FALSE
  )

  priority_summary <- data.frame()
  if (nrow(checklist) > 0 && all(c("Priority", "Severity") %in% names(checklist))) {
    priority_summary <- checklist |>
      dplyr::count(.data$Priority, .data$Severity, name = "Items") |>
      dplyr::arrange(
        factor(.data$Priority, levels = c("high", "medium", "low", "ready")),
        factor(.data$Severity, levels = c("required", "recommended", "optional"))
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  action_items <- data.frame()
  if (nrow(checklist) > 0) {
    action_items <- checklist |>
      dplyr::filter(.data$DraftReady %in% FALSE | .data$Available %in% FALSE) |>
      dplyr::arrange(
        factor(.data$Priority, levels = c("high", "medium", "low", "ready")),
        .data$Section,
        .data$Item
      ) |>
      dplyr::select(
        "Section", "Item", "Available", "DraftReady", "Severity",
        "Priority", "NextAction"
      ) |>
      utils::head(n = top_n) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  settings_tbl <- bundle_settings_table(object$settings)
  notes <- c(
    "This summary is a manuscript-preparation guide.",
    "DraftReady indicates that the corresponding reporting element can be drafted with the package's documented caveats; it does not certify inferential adequacy.",
    "Detailed FACETS positioning, software scope, and visual scope tables are available in `$facets_positioning`, `$software_scope`, and `$visual_scope`."
  )
  if (nrow(action_items) == 0) {
    notes <- c(notes, "No remaining draft-action rows were detected in the current checklist.")
  }

  out <- list(
    overview = overview,
    section_summary = section_summary,
    software_scope = software_scope,
    facets_positioning = facets_positioning,
    visual_scope = visual_scope,
    priority_summary = priority_summary,
    action_items = action_items,
    settings = settings_tbl,
    notes = notes,
    top_n = top_n
  )
  class(out) <- "summary.mfrm_reporting_checklist"
  out
}

#' @export
print.summary.mfrm_reporting_checklist <- function(x, ...) {
  cat("mfrmr Reporting Checklist Summary\n")

  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(as.data.frame(x$overview), row.names = FALSE)
  }
  if (!is.null(x$section_summary) && nrow(x$section_summary) > 0) {
    cat("\nSection summary\n")
    print(as.data.frame(x$section_summary), row.names = FALSE)
  }
  if (!is.null(x$priority_summary) && nrow(x$priority_summary) > 0) {
    cat("\nPriority summary\n")
    print(as.data.frame(x$priority_summary), row.names = FALSE)
  }
  if (!is.null(x$action_items) && nrow(x$action_items) > 0) {
    cat("\nAction items (preview)\n")
    print(as.data.frame(x$action_items), row.names = FALSE)
  }
  if (!is.null(x$facets_positioning) && nrow(x$facets_positioning) > 0) {
    cat("\nFACETS positioning\n")
    pos_cols <- intersect(c("Topic", "RecommendedWording"), names(x$facets_positioning))
    print(as.data.frame(x$facets_positioning[, pos_cols, drop = FALSE]), row.names = FALSE)
  }
  if (!is.null(x$settings) && nrow(x$settings) > 0) {
    cat("\nSettings\n")
    print(as.data.frame(x$settings), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

collect_bias_tables_for_checklist <- function(bias_results) {
  bias_results <- validate_bias_results_input(
    bias_results,
    helper = "reporting_checklist()"
  )
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
