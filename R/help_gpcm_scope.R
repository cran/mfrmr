#' GPCM workflow and score availability
#'
#' @description
#' Check which `GPCM` workflows can be used in `mfrmr`, the limits that
#' apply to each workflow, and the recommended alternative when a route is not
#' available.
#'
#' The table is intended for route selection before or after fitting and is
#' limited to workflow availability, interpretive constraints, and the route
#' to use next.
#' Except for the multiple-slope and EM rows, the table describes one facet
#' supplying relative discriminations. The provisional two-family route has
#' a separate scope below. The one-family model uses one slope facet. MML permits
#' a different facet for category steps; JML requires `slope_facet == step_facet`.
#' It has one substantive ability dimension. Separate-owner MML supports
#' fitted-object scoring, information, slope/curve uncertainty, matched PCM
#' comparison and saved inference. Weighting reviews and simulation/design
#' workflows still require the same owner; same-design [bootstrap_mfrm_gpcm()]
#' supports separate owners. Numerical eligibility is not a coverage guarantee. These structural choices are stated separately from the
#' availability of each output. An available probability or descriptive comparison does
#' not establish eligibility for a confidence interval or model-selection rule.
#'
#' JML fits estimate relative slopes but do not support the MML slope interval,
#' bootstrap, curve-uncertainty or inferential comparison routes described here.
#' Their available facet/location SEs are exploratory observation-table
#' approximations, not slope SEs or nuisance-adjusted joint-information SEs.
#' Local curvature checks do not establish bias control or interval coverage.
#' Fitted-object
#' JML scoring uses post-hoc EAP with a default standard-normal reference prior
#' or an explicit `scoring_prior`; it is not ML/WLE scoring. Portable GPCM extraction
#' supports MML and shared-owner JML within their distinct source-check scopes.
#' JML local checks leave incomplete global audits unchanged. See [fit_mfrm()]
#' for JML estimation and boundary
#' conventions and [predict_mfrm_units()] for conditional scoring.
#'
#' @param status Which rows to return. For `view = "workflows"` or
#'   `"unavailable"`, use `"all"` (default), `"supported"`,
#'   `"supported_with_caveat"`, `"blocked"`, or `"deferred"`.
#'   For `view = "scores"`, use `"all"`, `"available_with_caveat"`,
#'   `"supporting_route"`, or `"unavailable"`.
#' @param view The table to read: `"workflows"` (default) gives the general
#'   map; `"scores"` describes package-native versus FACETS-style score
#'   quantities; `"unavailable"` lists blocked/deferred routes and alternatives.
#'
#' @details
#' `Status` has the following user-facing meanings:
#'
#' - `supported`: the helper is available within the stated boundary;
#' - `supported_with_caveat`: the helper runs, but its interpretation is
#'   restricted as described in `Boundary`;
#' - `blocked`: the helper intentionally stops for a `GPCM` fit;
#' - `deferred`: no public `mfrmr` route is currently available.
#'
#' Read `Boundary` before interpreting a caveated result. For a blocked or
#' deferred row, use `RecommendedRoute` to choose a supported analysis or a
#' Rasch-family alternative.
#'
#' In the score view, `available_with_caveat` denotes an available quantity
#' within its `Limitation`, `supporting_route` denotes related evidence or a
#' reference model rather than a replacement for that quantity, and
#' `unavailable` denotes a quantity not supplied. Package-native GPCM scores
#' and uncertainty are not FACETS score-side equivalents. These statuses are
#' kept separate from the general workflow statuses to preserve their meaning.
#'
#' In 0.2.4, this is the single public entry for all three tables. Replace
#' `gpcm_score_side_contract(status)` with
#' `gpcm_capability_matrix(status, view = "scores")`, and
#' `gpcm_runtime_guard_coverage()` with
#' `gpcm_capability_matrix(view = "unavailable")`. Each view retains its former
#' column schema and status meanings. These tables describe supported routes;
#' they do not inspect a particular fit or certify its inferential eligibility.
#'
#' @section Connected outputs and unavailable extensions:
#' For the current one-slope-family model, use `summary(fit)` to review the
#' model and numerical status, `plot(fit, type = ...)` for available location,
#' fit and category views, and [stats::confint()] / [mfrm_curve_intervals()] for
#' separately checked MML uncertainty. An interval result has its own `print()`
#' and `plot()` methods; attach selected results to [mfrm_results()] for
#' [mfrm_report()] and [export_mfrm_results()]. A plot of locations or fit
#' statistics is not a plot of slope uncertainty.
#' [mfrm_facet_intervals()] separately supplies native-scale model-information
#' Wald intervals for one-family MML locations and within-facet contrasts.
#' Direct fixed/adaptive MML, centered additive facets/steps and unit weights
#' are required, with an estimated intercept-only normal population or the
#' explicit fixed-N(0,1) restriction. Shared/separate slope and step owners are
#' supported. Full-information and finer-grid checks preserve refusal reasons;
#' regular asymptotic validity assumes a correctly specified, identified
#' interior marginal model with nonsingular limiting information.
#' Finite-sample coverage is not guaranteed. Explicit zero-sum contrasts can
#' be standardized by the same fit's population SD using the full delta method;
#' individual standardized-location intervals are not supplied.
#'
#' Portable calibration has its own row and [mfrm_calibration_capabilities()]
#' gives estimator-specific source restrictions. Two slope families are available
#' provisionally through [fit_mfrm()] with fixed-standard-normal MML, exactly
#' two facets, no anchors and unit weights. Choose fixed-grid EM or adaptive
#' direct MML explicitly; see [fit_mfrm()]. Use `summary(fit)`, then
#' inspect adaptive initialization in `fit$opt$mml_initialization`: the default
#' compares neutral and EM-derived starts with the same adaptive likelihood,
#' retaining failed starts and the selected candidate's convergence status.
#' The saved initialization policy also follows quadrature-order refits. Use
#' `curves <- mfrm_curve_intervals(fit, newdata)` and `plot(curves)`. Despite the
#' function name, this route returns fitted values with unavailable intervals.
#' Attach them with `mfrm_results(fit, intervals = curves)` for reports and
#' saved exports. Without attachments, `mfrm_results(fit)` collects the saved
#' fit and its numerical status; neither call calculates new diagnostics.
#' Separately use `ci <- confint(fit, parm = "slopes")` for experimental component log-Wald intervals
#' and attach `intervals = list(slopes = ci, curves = curves)`. Local numerical
#' checks do not establish global identification or sampling coverage; failed
#' checks retain missing bounds.
#' An explicit `confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t1"))`
#' profiles one two-family component, using the actual owner/level names.
#' It preserves the fitted fixed or adaptive integration method; adaptive nodes,
#' scales and Jacobians are recomputed throughout the constrained searches.
#' It retains nuisance reoptimization, numerical checks, unavailable endpoints
#' and a same-target Wald comparison. Neither method has qualified coverage.
#' Use [mml_quadrature_sensitivity()] to compare
#' refits with the same engine and integration method. Both routes retain
#' each order's experimental log-Wald checks. Optional `adaptive_quad_points` also adds fixed-calibration
#' integration checks for two families; the posterior moments in that review
#' are numerical diagnostics, not a Person-scoring workflow. Person-score
#' comparisons remain unavailable for two families.
#' For both engines, [mfrm_response_diagnostics()] supplies same-data posterior
#' predictive residuals with fixed calibration, including descriptive Infit/Outfit
#' without reference cutoffs. It retains the fitted integration method and
#' complete Person conditioning record, including when selecting output rows.
#' Attach the saved object through `response_diagnostics`
#' to the results call above for plots, reports and exports.
#' [predict_mfrm_units()] separately supplies experimental conditional new-Person
#' EAP, posterior SD and continuous posterior intervals under the retained
#' N(0,1) prior, with source/batch numerical checks and no calibration uncertainty.
#' [extract_mfrm_calibration()] provides portable two-family format 6 with the
#' same fixed prior and separate source/batch checks.
#' [mfrm_facet_intervals()] supplies separately checked experimental normal
#' intervals for either owner's locations or prespecified within-facet contrasts.
#' It uses the constrained location block of the inverse full marginal information;
#' failed numerical checks retain point estimates with missing intervals.
#' Attach the result as `intervals = list(locations = ci)` for saved reports and
#' plots. Coverage is unqualified; location differences are not uniform rating
#' differences when slopes or steps vary. Ordinary fit diagnostics, step/curve
#' intervals and model ranking
#' remain unavailable. Shared-owner corrected
#' JML has an explicit experimental point-estimation and reporting route through
#' `jml_correction_order`; formal structural intervals remain unavailable.
#' Its descriptive response diagnostics and conditional EAP scoring use their
#' own saved identities and source checks. Portable corrected calibration uses
#' file format 5; no ordinary JML likelihood checks or structural intervals
#' are inherited. Separate slope and step owners alone still
#' define one slope family.
#'
#' @section Model, estimation and algorithm:
#' RSM, PCM and GPCM describe response probabilities, not a particular fitting
#' algorithm. MML integrates over an ability distribution; JML estimates the
#' training persons' abilities jointly with the other parameters. EM and direct
#' optimization are numerical ways to fit an MML model. RSM/PCM can therefore
#' also use MML--EM. In this package, MML defaults to direct optimization;
#' RSM/PCM EM and hybrid support additive, fixed-population models with fixed
#' integration (`population = NULL`). One-family GPCM EM/hybrid requests fall back to direct
#' and record the engine in `fit$summary`.
#'
#' The provisional two-slope-family GMFRM route uses fixed-standard-normal MML
#' with fixed-grid generalized EM (`mml_engine = "em"`) or adaptive direct
#' optimization (`mml_engine = "direct", mml_integration = "adaptive"`).
#' `mml_integration = "auto"` preserves the selected engine and increases
#' integration points up to `mml_accuracy_control = list(max_points = ...)`,
#' including reference rules, while retaining the numerical and information
#' checks. Fixed integration remains the package default; two-family direct
#' fitting requires explicit adaptive or automatic integration. Adaptive EM
#' and two-family hybrid fitting are unavailable. EM is not a defining property
#' of GMFRM. See the GPCM
#' scope vignette for the distinction between the model, estimation method
#' and algorithm, with links to Muraki, TAM and sirt.
#'
#' @section Applications of two slope families:
#' In a speaking assessment, task difficulty and assessor severity describe
#' shifts in ratings, whereas task and assessor slopes describe how responses
#' change with ability. Their product is the effective slope for a rating.
#' For example, 0.8 times 1.2 gives 0.96; it is a log-odds multiplier, not an
#' accuracy percentage or a score weight. The same roles can use piece/judge
#' or station/examiner labels, provided the single-ability assumption is
#' appropriate. Two slope families do not create two latent abilities or a
#' separate free slope for every task-by-assessor pair.
#'
#' Potential uses include task review and assessor feedback. A high assessor
#' slope is not a competence threshold, and a low task slope is not a rule for
#' deleting content. Read category and information curves over the relevant
#' ability range, inspect uncertainty and assignment overlap, and bring in
#' reference ratings or substantive evidence for accuracy claims. The vignette
#' explains these uses and a provisional two-family fitting example. It does
#' not yet provide a two-family model comparison or qualified feedback decision.
#' It also shows how to import the empirical writing table
#' `sirt::data.ratings1` from that separately installed package, preserve its
#' categories, and review unequal assignment. The example distinguishes new
#' fitted response curves from evidence of predictive accuracy. The original
#' fixed-grid fits had integration-sensitive curves and incomplete residual
#' summaries. A subsequent adaptive fit resolved the retained numerical example;
#' that agreement does not validate model fit or feedback decisions. All score
#' datasets bundled with mfrmr are synthetic.
#'
#' @section Rankings and consequential decisions:
#' An official competition result, highest latent ability and a future winner
#' are different targets. The scope vignette reviews actual figure-skating
#' score protocols, changing judge panels and advancement to a final.
#' Individual score intervals, high rank correlation or a high G coefficient
#' do not establish the probability of selecting the correct champion.
#' Differences need covariance and selection-aware uncertainty; independently
#' drawing from printed SEs omits shared calibration uncertainty.
#' There is no winner-probability or simultaneous Person-rank confidence-set
#' API. Available one-family and experimental two-family new-Person
#' EAP intervals condition on the calibration and prior; they do not provide
#' a validated winner-selection procedure. Preserve competition-specific
#' aggregation, rounding, tie rules and advancement separately from modelled
#' ability. Do not interpret sensitivity to a judge's marks as proof of bias.
#'
#' @section Local independence, testlets and random effects:
#' Two fixed slope families change response sensitivity; they do not by
#' themselves remove dependence between ratings from the same performance.
#' A testlet effect is itself a random effect, defined by its sharing unit.
#' A testlet model assumes independence conditional on ability and its local
#' effect, with dependence remaining after the effect is integrated out.
#' One substantive ability can therefore coexist with additional latent
#' dependence variables.
#'
#' [fit_mfrm_testlet()] supplies a Person-local RSM block effect; reusing a
#' block label for another Person creates a different effect.
#' [fit_mfrm_random_rater()] instead supplies a rater severity effect shared
#' across persons. These separate RSM routes neither add effects to GPCM nor
#' jointly estimate shared-rater and testlet effects. Random discrimination is
#' also different from estimating two fixed slope families.
#' The GPCM scope vignette explains effect-sharing units, observed versus new
#' block/rater prediction and the distinction from multivariate G-theory.
#' A testlet variance does not diagnose halo or enforce equal task weights.
#'
#' @section Interpreting the former bounded GPCM label:
#' Use GPCM as the model name and state the structure and output restrictions
#' explicitly. The older label described limited implementation scope, not a
#' distinct unidimensional response model. Positive log-parameterized relative
#' slopes retain their geometric-mean-one identification; there is no extra
#' user-facing finite slope box. Rejection of numerical overflow/underflow is
#' not clipping to a valid estimate. Separate owners, MML inference/comparison
#' and portable calibration are available within their respective scopes;
#' two-family fitting provides provisional estimates and conditional curves,
#' with separately checked experimental component-slope intervals. JML slope
#' intervals remain unavailable.
#' A data-specific boundary or unbounded interval is a separate statistical
#' issue, and retirement of the label does not certify interval coverage.
#'
#' @section Estimates, intervals and comparisons:
#' For free slopes, these are different questions:
#' - **What did the numerical fit return?** `fit$slopes$OptimizerEstimate`
#'   retains the fitted relative slopes for descriptive sensitivity analysis.
#'   The compatibility column `Estimate` contains the same numerical values;
#'   it does not override `ParameterStatus` or `PrimaryEstimate`.
#' - **How uncertain is a slope?** `confint(fit, parm = "slopes")` returns
#'   approximate pointwise intervals for eligible GPCM MML fits.
#'   `diagnose_mfrm(fit)$parameter_uncertainty$slopes` supplies the same 95%
#'   calculation with `CIEligible` and `InferenceReview`. Ineligible solutions
#'   retain missing ordinary bounds and explicitly labelled `Optimizer*`
#'   diagnostic quantities. Old eligibility flags do not authorize an interval.
#' - **Which model should be selected?** MML information criteria may be
#'   retained numerically. `compare_mfrm()` ranks GPCM MML candidates when its
#'   separate solution and comparison checks pass. `ICSelectable` describes
#'   likelihood and integration requirements; `ICFitEligible` describes the
#'   fit's solution check; `ICComparable` is the final comparison decision.
#'   The weighting review preserves this decision without making an
#'   operational-scoring recommendation. With `nested = TRUE`, it can also
#'   request the separately checked PCM/GPCM equal-slope test.
#'
#' The inference restrictions above concern the current package implementation;
#' they are not a claim that GPCM inference is impossible in general. A local
#'   rank or curvature check does not by itself assess competing solutions,
#'   numerical integration error or the performance of an interval procedure.
#' Relative-slope intervals have their own checks. A PCM/GPCM LRT requires the matched
#' comparison described below; more iterations or quadrature points alone do
#' not establish its structural assumptions.
#'
#' @section Different requirements for intervals and model comparison:
#' These decisions are separate, and do not follow from unidimensionality:
#' - **Slope intervals:** [confint.mfrm_fit()] with `parm = "slopes"` uses the inverse joint observed
#'   information, including estimated population parameters, and the sum-zero
#'   log-slope transformation. It checks likelihood consistency, convergence,
#'   positive unregularized information, unit weights and a grid of at least
#'   31 points. The exponentiated log-Wald limits are pointwise model-based
#'   approximations for geometric-mean-one relative slopes by default.
#'   Explicit options add population-SD-standardized slopes, named ratios or
#'   differences, Bonferroni adjustment and independent-cluster sandwich
#'   covariance. An experimental `method = "profile"` with one named `slope`
#'   reoptimizes all nuisance parameters, including the normal population, and
#'   retains endpoint failures and a saved likelihood plot. Its initial scope
#'   excludes anchors, interactions, covariates and standardized targets.
#'   Profile coverage and superiority over Wald have not been established.
#'   Small samples and misspecification can still affect coverage.
#'   Failed checks retain missing limits and a reason. [bootstrap_mfrm_gpcm()]
#'   provides fitted-model bootstrap intervals or a matched PCM/GPCM test;
#'   [mfrm_curve_intervals()] propagates calibration uncertainty to curves.
#'   Probability-curve intervals showed undercoverage in a saved-fit study of
#'   small incomplete designs, including finite-grid Bonferroni families.
#'   Numerical availability and multiplicity adjustment do not certify nominal
#'   coverage; bootstrap coverage requires separate evidence too.
#'   Saved results support [apa_table()], [plot_data()] and [as_ggplot()];
#'   attach selected intervals to [mfrm_results()] for reports and exports.
#'   Their targets remain separate from Wright/Pathway location and fit displays.
#' - **Information criteria:** AIC/BIC compare maximized likelihoods on the same
#'   response data, with appropriate free-parameter counts and numerical
#'   accuracy. They do not require a slope confidence interval or nested
#'   models. The GPCM MML solution check reevaluates the retained likelihood,
#'   terminal gradient and positive unregularized local information without
#'   refitting. It uses the existing numerical-gradient tolerance (at most
#'   \eqn{10^{-4}}) and information-inversion eigenvalue tolerance. The existing
#'   joint information calculation is shared with slope intervals and governed
#'   by `options(mfrmr.max_information_bytes = 256 * 1024^2)` rather than an
#'   80-coordinate cutoff. The budget estimates dense matrix workspace, not
#'   total process memory. An unavailable check gives a reason, not permission
#'   to rank. Local checks do
#'   not prove global optimality or integration accuracy; inspect different
#'   starts and [mml_quadrature_sensitivity()] when the decision is close.
#' - **PCM/GPCM likelihood-ratio test:** with the same population model, step
#'   structure and other constraints, setting all relative slopes to one gives
#'   PCM. One is an interior positive slope value, not a variance-zero boundary.
#'   With \eqn{G} slope levels and no other differing free parameters, the null
#'   imposes \eqn{G-1} independent log-slope restrictions. A chi-square reference
#'   additionally requires identified, regular solutions and adequate sample
#'   information. [compare_mfrm()] with `nested = TRUE` checks the matched
#'   model settings, G-1 free dimensions and regular local MML solutions before
#'   reporting an asymptotic chi-square p-value. Default PCM and GPCM calls can
#'   use different population models: supply `population_formula = ~1` and the
#'   same person data to both fits to compare estimated-normal models.
#'   Small or sparse samples can give inaccurate asymptotic p-values; examine
#'   starting-value and quadrature sensitivity and report the test's assumptions.
#'
#' The official \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/AIC.html}{R AIC documentation}
#' describes likelihood comparability. The \href{https://philchalmers.github.io/mirt/reference/mirt.html}{mirt model documentation}
#' documents GPCM and information-matrix SEs, and its
#' \href{https://philchalmers.github.io/mirt/reference/anova-method.html}{model-comparison documentation}
#' describes likelihood-ratio and information-criterion comparisons. These are
#' examples of supported statistical methods, not validation of mfrmr's
#' many-facet implementation.
#'
#' @section Comparing models with ConQuest and TAM:
#' Match the response formula before comparing estimates. In mfrmr, the
#' selected positive slope multiplies ability, facet locations and the category
#' step together. A slope multiplying ability alone, with separately additive
#' rater severity, generally specifies a different many-facet model.
#'
#' TAM's `tam.mml.mfr()` does not estimate slopes itself, but its documented
#' Example 14, Model 14c combines a facet intercept design with grouped slopes
#' in `tam.mml.2pl(irtmodel = "GPCM.design")`. ConQuest estimates GPCM scores
#' with `scoresfree`; its default slopes belong to combinations of facets
#' (generalized items), with further grouping available through a scoring
#' design. Neither construction automatically reproduces mfrmr's single
#' slope family and complete-predictor multiplication. See the
#' \href{https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html}{TAM fitting documentation}
#' and \href{https://www.acer.org/files/Note_8--The_ConQuest_4_Model.pdf}{ConQuest Note 8}.
#'
#' A matched item-only, positive-slope GPCM with an estimated normal population
#' can be expressed on either scale. mfrmr fixes the geometric mean of relative
#' slopes \eqn{\alpha_i} to one and estimates the population SD \eqn{\sigma}
#' conditional on any population covariates. On the unit-variance scale the
#' slopes become \eqn{a_i=\sigma\alpha_i}.
#' Locations, steps and any population regression also need transformation.
#' This equivalence does not include imposing both unit variance and
#' geometric-mean-one slopes: together they impose an additional restriction.
#'
#' Transformed intervals require the joint parameter covariance. Multiplying
#' relative-slope interval endpoints by an estimated population SD omits its
#' uncertainty and covariance with the slopes. Use [confint.mfrm_fit()] with
#' `parm = "slopes", scale = "standardized"` for the full transformation; with covariates the SD
#' is the residual population SD. TAM's documented `tam.se()` omits
#' parameter covariances; ConQuest distinguishes the covariance of parameter
#' estimates from the latent-population covariance. Neither marginal SEs nor
#' a latent-population covariance matrix replace the required joint matrix.
#' Numerical replication requires matching data, model and identification.
#' IC selection can compare different, nonnested models on compatible
#' likelihoods; verify each model's free-parameter count and integration
#' accuracy. A likelihood-ratio test additionally requires a nested null.
#' Consult `vignette("mfrmr-gpcm-scope")` for examples, sources and the scope
#' of existing numerical comparisons. These do not supply a general GPCM
#' import or comparison API.
#'
#' @section Conditional Person uncertainty:
#' Person posterior SDs and intervals, where returned, condition on the fitted
#' calibration; they are not slope intervals or calibration-aware confidence
#' intervals. An unstable calibration also limits their interpretation.
#'
#' @return For `view = "workflows"`, a data.frame of class
#' `mfrmr_gpcm_capabilities` with one row per workflow family and columns:
#' - `Area`
#' - `Helpers`
#' - `Status`
#' - `Boundary`
#' - `RecommendedRoute`
#'
#' For `view = "scores"`, an ordinary data.frame with `Capability`, `Status`,
#' `Limitation` and `Alternative`. For `view = "unavailable"`, an ordinary
#' data.frame with `Area`, `Helper`, `Status`, `Boundary` and `RecommendedRoute`.
#'
#' @section Typical workflow:
#' 1. Call `gpcm_capability_matrix()` before using `GPCM` in a new workflow.
#' 2. For `supported_with_caveat`, read `Boundary` before interpreting output.
#' 3. For `blocked` or `deferred`, follow `RecommendedRoute` instead.
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [compute_information()],
#'   [predict_mfrm_units()], [sample_mfrm_plausible_values()],
#'   [reporting_checklist()], [mfrmr_workflow_methods], [mfrmr-package]
#' @examples
#' gpcm_capability_matrix()
#' gpcm_capability_matrix("supported")
#' gpcm_capability_matrix("blocked")
#' gpcm_capability_matrix(view = "scores")
#' gpcm_capability_matrix("unavailable", view = "scores")
#' gpcm_capability_matrix(view = "unavailable")
#' @aliases gpcm_score_side_contract gpcm_runtime_guard_coverage
#' @concept GPCM boundaries
#' @concept route selection
#' @export
gpcm_capability_matrix <- function(status = "all", view = c("workflows", "scores", "unavailable")) {
  view <- match.arg(view)
  if (view == "scores") return(.gpcm_score_capabilities(status))
  status <- match.arg(status, c("all", "supported", "supported_with_caveat", "blocked", "deferred"))
  if (view == "unavailable") {
    out <- .gpcm_unavailable_routes()
    if (status != "all") out <- out[out$Status == status, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  out <- .gpcm_capability_registry(status)
  out <- out[, c(
    "Area", "Helpers", "Status", "Boundary", "RecommendedRoute"
  ), drop = FALSE]
  out$Helpers[out$Status == "deferred"] <- NA_character_
  rownames(out) <- NULL
  class(out) <- c("mfrmr_gpcm_capabilities", "data.frame")
  out
}

#' @export
#' @method print mfrmr_gpcm_capabilities
#' @noRd
print.mfrmr_gpcm_capabilities <- function(x, ...) {
  status_levels <- c(
    "supported", "supported_with_caveat", "blocked", "deferred"
  )
  counts <- table(factor(x$Status, levels = status_levels))
  status_summary <- data.frame(
    Status = names(counts),
    Routes = as.integer(counts),
    stringsAsFactors = FALSE
  )
  status_summary <- status_summary[status_summary$Routes > 0L, , drop = FALSE]

  cat("mfrmr GPCM workflow availability\n")
  cat("Most rows describe one slope family; two-family MML has a separate provisional scope.\n")
  cat("MML IC comparison and PCM/GPCM tests have separate checks; relative-slope intervals use separate MML checks.\n\n")
  print.data.frame(status_summary, row.names = FALSE)

  max_rows <- min(nrow(x), 8L)
  if (max_rows > 0L) {
    cat("\nRoute preview\n")
    preview <- x[seq_len(max_rows), c("Area", "Status"), drop = FALSE]
    print.data.frame(preview, row.names = FALSE)
  }
  if (nrow(x) > max_rows) {
    cat("\n... ", nrow(x) - max_rows, " more route(s).\n", sep = "")
  }
  cat(
    "\nFilter by status, for example ",
    "gpcm_capability_matrix(\"supported_with_caveat\").\n",
    "Read Boundary and RecommendedRoute before interpreting a caveated or ",
    "unavailable route.\n",
    sep = ""
  )
  invisible(x)
}

.gpcm_capability_registry <- function(status = c("all", "supported", "supported_with_caveat", "blocked", "deferred")) {
  status <- match.arg(status)

  out <- .mfrmr_table_.gpcm_capability_registry_out

  out <- rbind(out,.mfrmr_table_.gpcm_capability_registry_data)
  if (!identical(status, "all")) {
    out <- out[out$Status == status, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

.gpcm_score_capabilities <- function(status = "all") {
  status <- as.character(status[1])
  public_status <- c(
    "all", "available_with_caveat", "supporting_route", "unavailable"
  )
  legacy_status <- c(
    "implemented_with_caveat", "required_for_full_facets_review",
    "validated_dependency"
  )
  if (is.na(status) || !nzchar(status) ||
      !status %in% c(public_status, legacy_status)) {
    stop(
      "`status` must be one of: ",
      paste(sprintf('"%s"', public_status), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  registry <- .gpcm_score_side_contract_registry()
  if (status %in% legacy_status) {
    registry <- registry[
      registry$CompatibilityStatus == status, , drop = FALSE
    ]
  } else if (!identical(status, "all")) {
    registry <- registry[registry$PublicStatus == status, , drop = FALSE]
  }
  out <- registry[, c(
    "Capability", "PublicStatus", "Limitation", "Alternative"
  ), drop = FALSE]
  names(out)[names(out) == "PublicStatus"] <- "Status"
  rownames(out) <- NULL
  out
}

.gpcm_score_side_contract_registry <- function() {
  .mfrmr_table_.gpcm_score_side_contract_registry_data
}

.gpcm_unavailable_routes <- function() {
  matrix <- .gpcm_capability_registry()
  registry <- .gpcm_runtime_guard_registry()
  idx <- match(registry$CapabilityID, matrix$CapabilityID)
  out <- data.frame(
    Area = registry$Area,
    Helper = registry$Helper,
    Status = matrix$Status[idx],
    Boundary = matrix$Boundary[idx],
    RecommendedRoute = matrix$RecommendedRoute[idx],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

.gpcm_runtime_guard_registry <- function() {
  matrix <- .gpcm_capability_registry()
  guard <- data.frame(
    CapabilityID = c(
      "facets_score_review",
      "mcmc_backends",
      "corrected_jml_inference"
    ),
    Helper = c(
      "facets_output_contract_review()",
      NA_character_,
      NA_character_
    ),
    AvailabilityMode = c(
      "structured_error",
      "no_public_helper",
      "no_public_helper"
    ),
    ConditionClass = c(
      "mfrmr_gpcm_scope_error",
      NA_character_,
      NA_character_
    ),
    stringsAsFactors = FALSE
  )

  idx <- match(guard$CapabilityID, matrix$CapabilityID)
  guard$Area <- matrix$Area[idx]
  guard$Status <- matrix$Status[idx]
  guard$RecommendedRoute <- matrix$RecommendedRoute[idx]
  guard <- guard[, c(
    "CapabilityID", "Area", "Helper", "Status", "AvailabilityMode",
    "ConditionClass", "RecommendedRoute"
  )]
  rownames(guard) <- NULL
  guard
}
