#' Portable fixed-calibration capabilities
#'
#' Returns the model and estimator combinations supported by the portable
#' calibration workflow. This matrix concerns saved calibration artifacts;
#' it does not replace the wider fitted-object capabilities of [fit_mfrm()] or
#' [predict_mfrm_units()].
#' For JML fits, fitted-object scoring is post-hoc EAP with a standard-normal
#' reference prior by default, or an explicit `scoring_prior`. It is not
#' ML/WLE scoring or a population distribution estimated by JML. Portable
#' extraction supports RSM/PCM and shared-owner GPCM JML within their distinct
#' source-check scopes. No artifact stores training Person estimates.
#'
#' @return A data frame with one row per model, estimator, and scoring-basis
#'   combination. `PortableCalibration` is either `"available"` or
#'   `"unavailable"`; the listed scope and source checks still apply.
#' @export
#'
#' @examples
#' mfrm_calibration_capabilities()
mfrm_calibration_capabilities <- function() {
  .mfrmr_table_mfrm_calibration_capabilities_data
}

mfrmr_public_calibration_extraction_error <- function(error) {
  detail <- switch(
    as.character(error$code),
    MODEL_FAMILY_UNSUPPORTED = paste(
      "portable calibration is available only for RSM, PCM or GPCM fits;",
      "use fitted-object scoring for other model families"
    ),
    MODEL_ESTIMATOR_UNSUPPORTED = paste(
      "portable calibration supports the documented MML and JML models;",
      "check model and method compatibility before portable or fitted-object scoring"
    ),
    SCORING_BASIS_UNSUPPORTED = paste(
      "portable calibration supports fixed-normal RSM/PCM and two-family GPCM,",
      "or intercept-only estimated-normal one-family GPCM; use fitted-object scoring for other population models"
    ),
    NULL
  )
  if (is.null(detail)) stop(error)
  mfrmr_calibration_abort(error$code, error$field_path, detail)
}

mfrmr_validate_calibration_quadrature_review <- function(fit, review) {
  if (is.null(review)) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_REQUIRED", "quadrature_review",
      paste(
        "run `mml_quadrature_sensitivity()` on the source data, inspect its",
        "continuous differences, and pass the review with its highest-grid fit"
      )
    )
  }
  if (!inherits(review, "mfrm_quadrature_sensitivity")) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_INVALID", "quadrature_review",
      "the review must come from `mml_quadrature_sensitivity()`"
    )
  }
  nodes <- as.integer(review$settings$quad_points %||% integer())
  fits <- review$fits %||% list()
  runs <- as.data.frame(review$runs %||% data.frame(),
                        stringsAsFactors = FALSE)
  expected_names <- paste0("q", nodes)
  valid_shape <- length(nodes) >= 2L && !anyNA(nodes) &&
    identical(nodes, sort(unique(nodes))) &&
    length(fits) == length(nodes) && !anyDuplicated(names(fits)) &&
    setequal(names(fits), expected_names) &&
    nrow(runs) == length(nodes) &&
    all(c("Nodes", "EstimationConverged") %in% names(runs)) &&
    identical(as.integer(runs$Nodes), nodes) &&
    all(vapply(fits, inherits, logical(1L), what = "mfrm_fit"))
  if (!valid_shape) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_INVALID", "quadrature_review",
      "the review must retain one fitted model and run record for every evaluated quadrature order"
    )
  }
  fits <- fits[expected_names]
  model <- toupper(as.character(fit$config$model %||% "")[1L])
  if (!identical(as.character(review$settings$model %||% ""), model)) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_INVALID", "quadrature_review.settings.model",
      "the review model does not match the selected source fit"
    )
  }
  if (!all(vapply(fits, function(candidate) {
    identical(toupper(as.character(candidate$config$model %||% "")[1L]), model) &&
      identical(public_mfrm_method_label(candidate$config$method %||% ""),
                "MML") &&
      isTRUE(tryCatch(
        mfrmr_gqs_same_prepared_data(fits[[1L]], candidate),
        error = function(condition) FALSE
      ))
  }, logical(1L)))) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_INVALID", "quadrature_review.fits",
      "the reviewed fits do not share one model and prepared response dataset"
    )
  }
  same_settings <- tryCatch({
    comparison_model <- function(candidate) {
      model <- mfrm_checkpoint_objective_components(list(), candidate$config)$model
      if (identical(candidate$config$model, "GPCM")) {
        # Population estimates vary across grids by design. Compare their
        # specification, while convergence and the exact selected fit are
        # checked separately below.
        model$population_spec[c("coefficients", "sigma2", "source", "notes")] <- NULL
      }
      model
    }
    reference_arguments <- mfrmr_gqs_refit_arguments(fits[[1L]], NULL, 2L)
    reference_model <- comparison_model(fits[[1L]])
    all(vapply(seq_along(fits), function(i) {
      candidate <- fits[[i]]
      candidate_nodes <- candidate$config$estimation_control$quad_points
      identical(as.numeric(candidate_nodes), as.numeric(nodes[i])) &&
        identical(mfrmr_gqs_refit_arguments(candidate, NULL, 2L),
                  reference_arguments) &&
        identical(comparison_model(candidate), reference_model) &&
        identical(candidate$prep$score_map, fits[[1L]]$prep$score_map)
    }, logical(1L)))
  }, error = function(condition) FALSE)
  if (!isTRUE(same_settings)) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_INVALID", "quadrature_review.fits",
      paste(
        "each reviewed fit must use its recorded quadrature order and the same",
        "score map, model, anchors, integration mode and fitting settings"
      )
    )
  }
  if (!isTRUE(all(as.logical(runs$EstimationConverged))) ||
      !all(vapply(fits, function(candidate) {
        identical(as.integer(candidate$opt$convergence), 0L)
      }, logical(1L)))) {
    mfrmr_calibration_abort(
      "QUADRATURE_REVIEW_INCOMPLETE", "quadrature_review.runs",
      "every evaluated quadrature fit must complete estimation"
    )
  }
  source_nodes <- as.integer(
    fit$config$estimation_control$quad_points %||%
      fit$config$quad_points %||% NA_integer_
  )[1L]
  source_name <- paste0("q", source_nodes)
  if (!identical(source_nodes, max(nodes)) ||
      !source_name %in% names(fits) || !identical(fit, fits[[source_name]])) {
    mfrmr_calibration_abort(
      "QUADRATURE_SOURCE_NOT_HIGHEST", "fit",
      paste(
        "select the exact highest-grid fit stored in the reviewed object;",
        "the package does not choose a numerical tolerance for you"
      )
    )
  }
  invisible(TRUE)
}

#' Create and use a portable fixed calibration
#'
#' These functions implement a strict lifecycle for a saved, versioned
#' calibration. [extract_mfrm_calibration()] creates a draft from an eligible
#' `RSM` or `PCM` MML fit under the fixed standard-normal scoring basis, or
#' a one-family GPCM MML fit with an estimated intercept-only normal population that passes
#' the conditional source checks in [predict_mfrm_units()]. Experimental
#' two-family GPCM MML retains a fixed N(0,1) prior and its own source checks.
#' RSM/PCM and shared-owner GPCM JML support portable post-hoc EAP with a
#' standard-normal reference prior; explicit corrected GPCM JML preserves its
#' correction order and uses adjusted-equation checks, as described below.
#' [validate_mfrm_calibration()] and [freeze_mfrm_calibration()] are separate,
#' fail-closed transitions. Only a frozen artifact can be passed to
#' [score_mfrm_calibration()].
#'
#' The portable workflow supports one observed score scale, one latent
#' dimension and known non-Person facet levels. RSM/PCM MML retain stored anchors
#' and two-way facet interactions in file format 1. One-family GPCM MML uses file format 2,
#' retaining both step and slope owners, positive geometric-mean-one relative
#' slopes and the estimated population mean and SD. The slope multiplies the
#' entire adjacent-category predictor. Portable GPCM requires unit weights
#' and no anchors or interactions. One-family MML uses an estimated intercept-only
#' normal population. Shared-owner GPCM JML uses file format 4 and a post-hoc
#' N(0,1) reference prior; its distinct checks are described below.
#' Experimental two-family GPCM MML uses file format 6, with two ordered slope
#' owners, second-owner steps and the fixed N(0,1) prior. The first owner's
#' locations and log slopes are centered; the second owner's locations and
#' slopes remain free. Both fitting engines retain their fixed or adaptive
#' scoring algorithm. Source identity/convergence and finer-integration checks
#' must pass; review-only sources cannot be exported for frozen scoring.
#' No separate `quadrature_review` is accepted for this route. The file contains
#' component slopes and observed facet combinations, without training responses,
#' Person estimates or an executable fit. Known levels may form new combinations;
#' `row_dispositions$ObservedContext` identifies whether each pair occurred in
#' calibration. New combinations are model-based, not empirically validated.
#' Prior overrides and `event_id` extensions are unsupported; repeated Person-
#' facet cells are rejected. Facet names must not collide with scoring/disposition
#' columns (`Person`, `Score`, `Weight`, `InputRow`, `EventId`, `Disposition`,
#' `ReasonCode`, `CalibrationId`, `ObservedContext`).
#'
#' For format 6, use the same extract, validate, freeze, save, load and score
#' functions as below. Each new batch must pass its own integration checks;
#' validation and freezing do not establish sampling coverage, global
#' identification, boundary absence or population transport. Posterior
#' intervals exclude calibration uncertainty. Older readers refuse format 6;
#' existing formats 1--5 keep their original interpretation.
#' RSM/PCM JML uses file format 3, with unit weights, no anchors or interactions,
#' a finite identified source and a post-hoc N(0,1) reference prior. This prior
#' is not estimated by JML. Estimated-population RSM/PCM and latent regression
#' use fitted-object routes; other structures follow their documented capabilities.
#' See [mfrm_calibration_capabilities()].
#'
#' @section Anchors, centering and score uncertainty:
#' New exports distinguish whole-facet centering, legacy remainder centering
#' and uncentered facets in `constraints$identification`. Direct values and
#' group means are retained in `constraints$anchors`; validation checks the
#' explicit facet maps against the stored coordinates and free dimensions.
#' For example, if four raters are centered as a whole and the first three
#' have mean 0.10, the fourth is fixed at -0.30. Freezing retains all four
#' calibrated values, including the estimates of the free group members.
#'
#' Scoring uses those expanded values without centering or estimating them
#' again. A returned `SD` is a Person posterior SD conditional on the entire
#' frozen calibration and scoring prior. It is not a rater SE and excludes
#' uncertainty in both the estimated calibration and supplied anchors.
#' Use [mfrm_facet_intervals()] for supported conditional facet intervals.
#' With an anchored fitted model, [predict_mfrm_units()] accepts supplied
#' response rows; the narrower [score_mfrm_persons()] entry excludes anchors.
#'
#' Loading an older artifact preserves its values and original constraint
#' labels. Older generic labels do not fully identify the centering convention;
#' consult the source fit or re-extract it to obtain explicit labels. Changing
#' a label cannot convert a calibration to another scale.
#'
#' @details
#' One-family GPCM MML extraction freshly evaluates the native fit's local likelihood,
#' gradient, unregularized information and source-integration stability using
#' the same conditional source checks as [predict_mfrm_units()]. A failed or
#' unresolved source cannot be frozen by selecting a review policy. The
#' artifact retains the passing extraction decision, numerical comparison,
#' and original global boundary/identification states. Reading it validates
#' these records against its stored identity; it does not refit the source or
#' claim a global maximum, valid population transport or calibration intervals.
#' Every GPCM scoring call also checks EAP/SD against two adaptive reference
#' orders under the actual scoring prior. Failed score integration stops with
#' `SCORING_INTEGRATION_FAILED`; extract with more `scoring_quad_points` before
#' retrying. Endpoint and sparse-pattern review labels remain separate from
#' numerical accuracy. File format 1 retains its original default semantics.
#'
#' For RSM/PCM MML extraction, run [mml_quadrature_sensitivity()] on user-selected grids
#' and inspect its continuous differences. Pass the exact highest-grid fit in
#' that object together with `quadrature_review`. The retained fits must use
#' their recorded grid sizes and the same data, score map and fitting settings,
#' including anchors and integration mode, and must all have converged.
#' This procedural requirement
#' does not declare the fit numerically stable: the package does not choose the
#' application-specific tolerance or decide whether more grids are needed.
#' Archive the review separately when it is part of the audit trail, because
#' response-linked fits are deliberately not embedded in the portable artifact.
#' GPCM uses the fresh source checks above and does not require an additional
#' refitting study. If supplied, `quadrature_review` is still validated.
#'
#' RSM/PCM JML extraction requires current passing input, identification, category,
#' boundary and numerical checks from the fitted-object scoring policy. It
#' freshly evaluates the joint likelihood and gradient: the stored objective
#' must agree within 1e-6 and the maximum absolute gradient must be <= 1e-4.
#' These are numerical source criteria, not interval-coverage or bias guarantees.
#' Boundary or unresolved sources are refused; there is no review override.
#' No MML quadrature review is used because JML does not integrate over Persons
#' during calibration. Every new JML score batch receives adaptive-reference
#' EAP/SD checks under the actual scoring prior. Increase `scoring_quad_points`
#' on extraction if the scoring grid is inadequate. Source checks and the prior
#' origin persist through saving, loading and score summaries; replay checks
#' the stored evidence without re-estimating the training model.
#'
#' Experimental corrected GPCM JML uses file format 5 and its own source
#' checks. Extraction reconstructs the adjusted equation from the observed
#' ratings and checks the saved root, full equation Jacobian and parameter
#' tables. The maximum mean-equation residual must be <= 1e-7 and the
#' remaining Newton step <= 1e-5. These are numerical checks, not a global
#' uniqueness or bias guarantee. An unavailable covariance does not erase a
#' valid calibration. The artifact records correction order, assignment-sampling
#' assumption, local-root status and numerical checks without training responses
#' or Person estimates. It is never labelled an ordinary JML likelihood maximum.
#' Its EAP and continuous posterior intervals use a separate normal scoring
#' prior; they exclude calibration uncertainty and do not supply corrected
#' Person ML/WLE or structural confidence intervals. Older file formats retain
#' their original checks. See [fit_mfrm()] for corrected estimation limits.
#'
#' For shared-owner uncorrected GPCM JML, extraction checks the joint likelihood and
#' gradient, and positive-definite, unregularized curvature over all free Person
#' and structural parameters (reciprocal condition number > 1e-10). Recorded
#' Person, additive or slope boundary certificates prevent automatic scoring;
#' finite optimizer traces cannot override them. Detected identification failure,
#' inadequate categories, failed numerical checks and unsupported structures
#' also prevent extraction. Otherwise, global identification/boundary audits
#' labelled as incomplete stay incomplete in the artifact. This qualifies
#' conditional EAP only, not global optimality or formal JML slope inference.
#' Source curvature currently uses a dense joint matrix, so its memory cost
#' grows quadratically with the number of free parameters, including Persons.
#' Frozen-artifact scoring does not repeat that calculation or use training data.
#' The objective, gradient and condition-number cutoffs are package numerical
#' criteria, not statistical thresholds established by the scoring literature.
#' A refused source does not imply that all its structural parameters are
#' inestimable: portable scoring does not currently qualify boundary-profile
#' calibrations with infinite training-Person estimates.
#'
#' Posterior EAP estimates, posterior standard deviations, and intervals are
#' conditional on the frozen point calibration and the scoring prior.
#' For GPCM MML, the prior mean and SD estimated during calibration are frozen;
#' a new cohort does not reestimate or replace them. Use `scoring_prior` to
#' examine another shared normal prior without changing the artifact. The
#' original and scoring prior identities remain separate in estimates and
#' settings, including summaries, plot data and exports. A supplied prior is
#' an analyst assumption, not a population estimate. Intervals exclude both
#' population-parameter and calibration-parameter
#' uncertainty. Loading validates structure and semantic consistency, but does
#' not authenticate an artifact from an untrusted source.
#' New v2 scoring algorithms invert the continuous posterior CDF for equal-tail
#' intervals. EAP and SD retain the stored quadrature rule. Saved v1 artifacts
#' preserve their discrete grid interval endpoints and carry a note that the
#' continuous posterior mass can differ from the requested interval level.
#'
#' @section Statistical basis and external comparisons:
#' Bock and Mislevy (1982, pp. 432-433) describe posterior-mean and
#' posterior-SD scoring from fixed response functions and a specified prior,
#' including multiple-category responses. This supports the scoring layer;
#' it does not validate JML calibration bias, the default reference prior for
#' a new cohort, or the numerical tolerances used to check a fitted calibration. Their adaptive
#' testing simulations do not establish coverage for sparse many-facet designs.
#' Muraki (1992) supplies the GPCM model basis and an EM estimation method,
#' not evidence that many-facet JML intervals have accurate coverage.
#'
#' External checks of fixed probabilities, likelihoods and EAP/SD must align
#' the response functions and prior. They do not compare free JML estimation.
#' ConQuest documents that JML cannot estimate item scores (the free scores
#' used for GPCM discrimination). TAM's `tam.jml()` uses a supplied loading
#' array and defaults to extreme-score adjustment and item-bias correction;
#' its default output is not the same estimator as unadjusted JML. Do not
#' infer agreement of estimators from agreement of fixed-calibration scores.
#'
#' @section Choose what to save:
#' Save a fitted analysis or [mfrm_results()] with `saveRDS()` when you need
#' to revisit its diagnostics and reports. A portable `mfrm_calibration` has a
#' different purpose: scoring new response rows without the training fit or
#' training responses. Use `save_mfrm_calibration()` and
#' `load_mfrm_calibration()` for that artifact. Save an already computed
#' `mfrm_calibration_score` with `saveRDS()` to reproduce its tables and figures
#' without scoring the responses again. A CSV is a table, not a replacement
#' for any of these complete objects.
#'
#' Validation checks the saved parameters, scoring settings and model checks;
#' freezing marks that calibration as ready for scoring without changing it.
#' These steps do not establish model fit, lack of bias or suitability for another
#' population. New Persons may use known calibrated facet levels. An unknown
#' rater, task or score category requires an appropriate calibration workflow;
#' assigning an existing level's name is not a way to extend the calibration.
#' Missing responses are not zero scores. With `missing_response = "omit"`,
#' inspect row dispositions and Persons who could not be scored.
#'
#' Older saved objects are checked against their recorded format and evidence.
#' Inspect `review_mfrm_calibration()` if a file has missing or inconsistent
#' records. Do not change saved check results to force acceptance. If checks
#' of the original fitted model are missing, re-extract from that saved model
#' if it passes the current checks; when only the
#' score records are incomplete, score again from a valid frozen artifact.
#' Neither action means that model refitting is automatically necessary.
#'
#' @references
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
#' in a microcomputer environment. *Applied Psychological Measurement*,
#' 6(4), 431-444. \doi{10.1177/014662168200600405}.
#'
#' Muraki, E. (1992). A generalized partial credit model: Application of an
#' EM algorithm. *Applied Psychological Measurement*, 16(2), 159-176.
#' \doi{10.1177/014662169201600206}.
#'
#' [TAM joint maximum likelihood documentation](https://alexanderrobitzsch.github.io/TAM/reference/tam.jml.html)
#' and [ACER ConQuest command reference](https://conquestmanual.acer.org/s4-00.html).
#'
#' @param fit An eligible `mfrm_fit` produced by [fit_mfrm()].
#' @param calibration_id Optional nonempty calibration identifier.
#' @param source_fit_id Optional nonempty source-fit identifier.
#' @param created_at_utc Optional RFC3339 UTC timestamp. Omit it to use the
#'   current time.
#' @param scoring_quad_points Integer quadrature order of at least 2 used for
#'   later artifact scoring. It is independent of the fit-time quadrature and
#'   defaults to 31. The fixed or adaptive integration mode is inherited from
#'   the source fit and stored in the artifact's scoring algorithm identity.
#'   With adaptive scoring, stored nodes and weights define the base Hermite
#'   rule; their person-specific locations are recomputed for each new batch.
#' @param quadrature_review Required for RSM/PCM MML, optional for GPCM MML
#'   (which checks the fitted calibration directly); not used when checking JML fits.
#'   An `mfrm_quadrature_sensitivity` from
#'   [mml_quadrature_sensitivity()] for the same data and model. `fit` must be
#'   the exact highest-grid fit stored in this object. The package checks this
#'   procedural evidence but leaves acceptable numerical movement to the user.
#' @param calibration An `mfrm_calibration` object.
#' @param validated_at_utc,frozen_at_utc Optional RFC3339 UTC timestamps for
#'   reproducible lifecycle records. Omit them to use the current time.
#' @param record_id A nonempty identifier for a terminal superseded or retired
#'   record. It must differ from the frozen parent identifier.
#' @param superseded_at_utc,retired_at_utc Optional RFC3339 UTC timestamps.
#' @param file A path ending in `.rds`.
#' @param overwrite Whether an existing persistence target may be replaced.
#' @param new_data A data frame of response rows for new Persons.
#' @param person,score,weight,event_id Optional input-column names. Stored
#'   source-column names are used when applicable. `event_id` distinguishes
#'   otherwise duplicate response events.
#' @param facets Optional named character vector mapping stored facet names to
#'   input columns.
#' @param interval_level Central posterior interval level, strictly between 0
#'   and 1.
#' @param missing_response Either `"error"` or `"omit"`.
#' @param adaptive_quad_points Optional vector of at least two distinct integer
#'   orders >= 3, for example `c(31, 61)`, used only by
#'   `score_mfrm_calibration()`. Adds an unrounded `quadrature_review` comparing
#'   stored-grid results with mode/curvature-adapted integration for each Person.
#'   Inspect both fixed/adaptive differences and changes between adaptive orders.
#'   This does not change the artifact or replace the reported scores/intervals;
#'   adaptive results condition on the same point calibration.
#'   For GPCM or an explicitly supplied scoring prior, reference orders are
#'   added automatically and the comparison governs numerical score acceptance.
#'   Only Persons with scored responses have numerical review rows. `Status =
#'   "computed"` means the calculation finished, not that accuracy is certified;
#'   an unavailable row retains the reason in `Detail`.
#'
#' @param scoring_prior Optional `list(mean = ..., sd = ...)` for
#'   `score_mfrm_calibration()`. `NULL` retains the artifact's normal prior.
#'   Values use the unchanged calibration scale and must be finite with positive
#'   SD and representable positive variance. A supplied prior automatically
#'   receives numerical scoring checks and does not modify the frozen artifact.
#'   Legacy artifacts with discrete grid-endpoint intervals require re-extraction
#'   before using this option. The prior is common to all scored persons;
#'   covariate-dependent overrides are not supported. Two-family GPCM format 6
#'   requires `scoring_prior = NULL` and retains the fixed N(0,1) prior.
#' @param cut_scores Optional named, strictly increasing finite ability cuts
#'   for `score_mfrm_calibration()`. Ordinary RSM/PCM with unit weights and
#'   continuous posterior intervals only. Returns a `decisions` table of
#'   probabilities above cuts and in successive bands on the calibration scale.
#'   Calibration, prior and cuts are held fixed; this excludes their estimation
#'   uncertainty. Unscored Persons remain unavailable. Cut names do not establish
#'   CEFR alignment. See [predict_mfrm_units()] for interpretation and reporting.
#'
#' @inheritSection predict_mfrm_units Conditional ability decisions
#' @return `extract_mfrm_calibration()`, `validate_mfrm_calibration()`,
#'   `freeze_mfrm_calibration()`, `supersede_mfrm_calibration()`,
#'   `retire_mfrm_calibration()`, and `load_mfrm_calibration()` return an
#'   `mfrm_calibration`. `review_mfrm_calibration()` returns a data frame of
#'   structured refusals, with zero rows for an acceptable object.
#'   `save_mfrm_calibration()` invisibly returns the normalized path.
#'   `score_mfrm_calibration()` returns an `mfrm_calibration_score` containing
#'   estimates plus row and Person dispositions and scoring identities.
#'   With `cut_scores`, it also retains `decisions`: conditional event
#'   probabilities with boundaries, labels, actual priors and availability.
#'   When requested and scored rows exist, `quadrature_review` contains the
#'   fixed/adaptive comparison; `summary()` preserves it and adds a compact
#'   `quadrature_overview`.
#'
#' @seealso [mml_quadrature_sensitivity()], [mfrm_calibration_score_methods]
#'   for an executable scoring example, concise review and visualization of
#'   returned score batches; [predict_mfrm_units()], [score_mfrm_persons()],
#'   [apa_table()] and [export_summary_appendix()] for saved-score tables.
#'
#' @name mfrm_calibration_workflow
NULL

#' @rdname mfrm_calibration_workflow
#' @export
extract_mfrm_calibration <- function(fit, calibration_id = NULL,
                                     source_fit_id = NULL,
                                     created_at_utc = NULL,
                                     scoring_quad_points = 31L,
                                     quadrature_review = NULL) {
  draft <- tryCatch(
    mfrmr_extract_calibration_draft(
      fit = fit,
      calibration_id = calibration_id,
      source_fit_id = source_fit_id,
      created_at_utc = created_at_utc,
      scoring_quad_points = scoring_quad_points
    ),
    mfrm_calibration_error = mfrmr_public_calibration_extraction_error
  )
  if (mfrm_has_product_slopes(fit)) {
    if (!is.null(quadrature_review)) stop("Two-family extraction uses its own source identity, convergence and finer-integration checks; quadrature_review is not accepted.", call. = FALSE)
  } else if (identical(fit$config$method, "JML")) {
    if (!is.null(quadrature_review)) stop("JML extraction does not use an MML quadrature review; scoring integration is checked for each new batch.", call. = FALSE)
  } else if (!identical(fit$config$model, "GPCM") || !is.null(quadrature_review)) {
    mfrmr_validate_calibration_quadrature_review(fit, quadrature_review)
  }
  draft
}

#' @rdname mfrm_calibration_workflow
#' @export
review_mfrm_calibration <- function(calibration) {
  mfrmr_review_calibration(calibration)
}

#' @rdname mfrm_calibration_workflow
#' @export
validate_mfrm_calibration <- function(calibration, validated_at_utc = NULL) {
  mfrmr_validate_calibration_draft(calibration, validated_at_utc)
}

#' @rdname mfrm_calibration_workflow
#' @export
freeze_mfrm_calibration <- function(calibration, frozen_at_utc = NULL) {
  mfrmr_freeze_calibration(calibration, frozen_at_utc)
}

#' @rdname mfrm_calibration_workflow
#' @export
supersede_mfrm_calibration <- function(calibration, record_id,
                                       superseded_at_utc = NULL) {
  mfrmr_supersede_calibration(calibration, record_id, superseded_at_utc)
}

#' @rdname mfrm_calibration_workflow
#' @export
retire_mfrm_calibration <- function(calibration, record_id,
                                    retired_at_utc = NULL) {
  mfrmr_retire_calibration(calibration, record_id, retired_at_utc)
}

#' @rdname mfrm_calibration_workflow
#' @export
save_mfrm_calibration <- function(calibration, file, overwrite = FALSE) {
  mfrmr_save_calibration(calibration, file, overwrite)
}

#' @rdname mfrm_calibration_workflow
#' @export
load_mfrm_calibration <- function(file) {
  mfrmr_load_calibration(file)
}

#' @rdname mfrm_calibration_workflow
#' @export
score_mfrm_calibration <- function(calibration,
                                   new_data,
                                   person = NULL,
                                   facets = NULL,
                                   score = NULL,
                                   weight = NULL,
                                   interval_level = 0.95,
                                   missing_response = "error",
                                   event_id = NULL,
                                   adaptive_quad_points = NULL,
                                   scoring_prior = NULL,
                                   cut_scores = NULL) {
  mfrmr_score_calibration(
    calibration = calibration,
    new_data = new_data,
    person = person,
    facets = facets,
    score = score,
    weight = weight,
    interval_level = interval_level,
    missing_response = missing_response,
    event_id = event_id,
    adaptive_quad_points = adaptive_quad_points,
    scoring_prior = scoring_prior,
    cut_scores = cut_scores
  )
}
