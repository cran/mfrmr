#' Fit a many-facet Rasch model with a flexible number of facets
#'
#' This is the package entry point. It wraps `mfrm_estimate()` and defaults to
#' `method = "MML"`. Any number of facet columns can be supplied via `facets`.
#'
#' @param data A data.frame in long format with one row per observed rating
#'   event.
#' @param person Column name for the person (character scalar).
#' @param facets Character vector of facet column names.
#' @param score Column name for the observed ordered category score. Values
#'   must be coercible to numeric integer category codes. Fractional values are
#'   rejected. Binary `0/1` or `1/2` responses are supported as the ordered
#'   two-category special case. When `keep_original = FALSE`, unused
#'   intermediate categories are collapsed to a contiguous internal scale and
#'   the mapping is recorded in `fit$prep$score_map`. If `rating_min` /
#'   `rating_max` are supplied and the observed scores are a contiguous
#'   subset of that range (for example a 1-5 scale with only 2-5 observed),
#'   the supplied full range is retained so zero-count boundary categories
#'   remain part of the fitted score support.
#' @param rating_min Optional minimum category value. Supply this with
#'   `rating_max` when the intended score scale includes unobserved boundary
#'   categories.
#' @param rating_max Optional maximum category value. Supply this with
#'   `rating_min` when the intended score scale includes unobserved boundary
#'   categories.
#' @param weight Optional weight column name.
#' @param keep_original Keep original category values.
#' @param model `"RSM"`, `"PCM"`, or bounded `"GPCM"`.
#' @param method `"MML"` (default) or `"JML"`. `"JMLE"` is accepted as a
#'   backward-compatible alias for the same joint-maximum-likelihood path.
#' @param step_facet Step facet for `PCM` and the bounded `GPCM`
#'   branch. For `GPCM`, this should be supplied explicitly rather than
#'   relying on an implicit default.
#' @param slope_facet Slope facet for the bounded `GPCM` branch. The
#'   current release requires `slope_facet == step_facet` and uses a
#'   positive-slope identification convention on the log scale with geometric
#'   mean discrimination fixed to 1.
#' @param anchors Optional anchor table.
#' @param group_anchors Optional group-anchor table.
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#' @param positive_facets Facets with positive orientation.
#' @param anchor_policy How to handle anchor-audit issues: `"warn"` (default),
#'   `"error"`, or `"silent"`.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   anchor-audit recommendations.
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in anchor-audit recommendations.
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in anchor-audit recommendations.
#' @param quad_points Quadrature points for MML.
#' @param maxit Maximum optimizer iterations.
#' @param reltol Optimization tolerance.
#' @param mml_engine MML optimization engine for `method = "MML"`:
#'   `"direct"` (default) uses direct BFGS on the marginal log-likelihood,
#'   `"em"` uses an EM loop for `RSM` / `PCM` with `population = NULL`, and
#'   `"hybrid"` uses EM as a warm start before the direct optimizer. Unsupported
#'   combinations currently fall back to `"direct"` and record that fallback in
#'   `fit$summary`.
#' @param population_formula Optional one-sided formula for a person-level
#'   latent-regression population model, for example `~ grade + ses`. In the
#'   current release, latent regression is implemented only for
#'   `method = "MML"` with a unidimensional conditional-normal population
#'   model.
#' @param person_data Optional one-row-per-person data.frame holding background
#'   variables for `population_formula`. Numeric, logical, factor, ordered
#'   factor, and character predictors are expanded through `stats::model.matrix()`;
#'   categorical xlevels and contrasts are stored for replay and scoring.
#'   Required when `population_formula` is supplied.
#' @param person_id Optional person-ID column in `person_data`. Defaults to
#'   `person` when that column exists in `person_data`.
#' @param population_policy How missing background data are handled for a
#'   latent-regression fit. `"error"` (default) requires complete person-level
#'   covariates; `"omit"` fits the model on the complete-case subset and records
#'   omitted persons / omitted response rows in the returned `population`
#'   metadata while retaining the observed-person-aligned pre-omit table for
#'   replay/export provenance.
#'
#' @details
#' Data must be in **long format** (one row per observed rating event).
#'
#' @section Model:
#' `fit_mfrm()` estimates the many-facet Rasch model (Linacre, 1989).
#' For a two-facet design (rater \eqn{j}, criterion \eqn{i}) the model is:
#'
#' \deqn{\ln\frac{P(X_{nij} = k)}{P(X_{nij} = k-1)} =
#'   \theta_n - \delta_j - \beta_i - \tau_k}
#'
#' where \eqn{\theta_n} is person ability, \eqn{\delta_j} rater severity,
#' \eqn{\beta_i} criterion difficulty, and \eqn{\tau_k} the \eqn{k}-th
#' Rasch-Andrich threshold.  Any number of facets may be specified via the
#' `facets` argument; each enters as an additive term in the linear
#' predictor \eqn{\eta}.
#'
#' With `model = "RSM"`, thresholds \eqn{\tau_k} are shared across all
#' levels of all facets.
#' With `model = "PCM"`, each level of `step_facet` receives its own
#' threshold vector \eqn{\tau_{i,k}} on the package's shared observed
#' score scale.
#'
#' With only two ordered categories (\eqn{K = 1}), the same adjacent-category
#' formulation reduces to the usual binary Rasch logit for the single category
#' boundary:
#'
#' \deqn{\ln\frac{P(X_{n\cdot} = 1)}{P(X_{n\cdot} = 0)} = \eta - \tau_1}
#'
#' With `method = "MML"`, person parameters are integrated out using
#' Gauss-Hermite quadrature and EAP estimates are computed post-hoc.
#' With `method = "JML"`, all parameters are estimated jointly as fixed
#' effects. `"JMLE"` remains an accepted compatibility alias, but package
#' output now uses `"JML"` as the public label. See the "Estimation methods"
#' section of [mfrmr-package] for details.
#'
#' @section Weighting policy:
#' `mfrmr` treats `RSM` / `PCM` as the equal-weighting reference route for
#' operational many-facet measurement. In that Rasch-family branch,
#' discrimination is fixed, so the scoring model does not differentially
#' reweight item-facet combinations through estimated slopes.
#'
#' bounded `GPCM` is supported as an alternative when users explicitly accept
#' discrimination-based reweighting. This often improves model fit, but the
#' package does not treat better fit alone as a sufficient reason to replace an
#' equal-weighting Rasch-family model.
#'
#' The `weight` argument is separate from that modeling choice. It supplies an
#' observation-weight column; it does not create a free-form facet-weighting
#' scheme and does not change the fixed-discrimination contract of `RSM` /
#' `PCM`.
#'
#' @section Input requirements:
#' Minimum required columns are:
#' - person identifier (`person`)
#' - one or more facet identifiers (`facets`)
#' - observed score (`score`)
#'
#' Scores are treated as ordered categories.
#' Non-numeric score labels are dropped with a warning after coercion, whereas
#' fractional numeric scores are rejected with an error instead of being
#' silently truncated.
#'
#' Binary responses are therefore supported as ordered two-category scores
#' (for example `0/1` or `1/2`) under the same `RSM` / `PCM` interface.
#' If your observed categories do not start at 0, set `rating_min`/`rating_max`
#' explicitly to avoid unintended recoding assumptions. For example, if the
#' intended instrument is a 1-5 scale but the current sample only uses 2-5,
#' set `rating_min = 1, rating_max = 5` to retain the zero-count category 1
#' in the score support.
#'
#' When `keep_original = FALSE`, observed gaps such as `1, 3, 5` are recoded
#' internally to a contiguous scale (`1, 2, 3`) and the mapping is stored in
#' `fit$prep$score_map`. To retain zero-count intermediate categories as part
#' of the original scale, set `keep_original = TRUE` in addition to supplying
#' the full `rating_min` / `rating_max` range.
#'
#' This is ordered binary support, not a separate nominal-response model.
#' In `PCM`, a binary fit still uses one threshold per `step_facet` level on
#' the shared observed-score scale.
#'
#' Supported model/estimation combinations in the current release:
#' - `model = "RSM"` with `method = "MML"` or `"JML"/"JMLE"`
#' - `model = "PCM"` with a designated `step_facet` (defaults to first facet)
#' - `model = "GPCM"` is currently implemented only for the narrow bounded
#'   branch with `slope_facet == step_facet`; `MML` and `JML` fitting, core
#'   summaries, fixed-calibration posterior scoring, [compute_information()],
#'   Wright/pathway/CCC fit plots, [diagnose_mfrm()], residual-PCA follow-up,
#'   [interrater_agreement_table()], [unexpected_response_table()],
#'   [displacement_table()], [measurable_summary_table()],
#'   [rating_scale_table()], [facet_quality_dashboard()],
#'   [reporting_checklist()], [category_structure_report()],
#'   [category_curves_report()], and graph-only
#'   [facets_output_file_bundle()] are available. Direct simulation
#'   specifications and data generation are also supported through
#'   [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()], and
#'   [simulate_mfrm_data()] when the slope-aware generator contract is stored
#'   explicitly. Fair-average reporting, planning/forecasting, scorefile
#'   exports, and broader APA/QC pipelines should still be treated as
#'   unsupported unless documented otherwise. Use [gpcm_capability_matrix()] as
#'   the formal boundary statement for the current `GPCM` scope.
#'
#' Latent-regression status:
#' - `population_formula = NULL` keeps the legacy unconditional `MML` / `JML`
#'   behavior.
#' - Supplying `population_formula` activates a first-version latent-regression
#'   branch for `method = "MML"` only.
#' - The current branch assumes a one-dimensional conditional-normal population
#'   model with person-specific quadrature nodes
#'   \eqn{\theta_{nq} = x_n^\top \beta + \sigma z_q}.
#' - Background variables must be supplied in `person_data`; numeric/logical
#'   columns and categorical factor/character columns are expanded through
#'   `stats::model.matrix()`.
#' - Current overlap with the ConQuest latent-regression documentation is
#'   limited to direct estimation from response data under a unidimensional
#'   `MML` population model with package-built model-matrix covariates. It
#'   should not be described as parity for arbitrary imported design matrices,
#'   multidimensional models, or the full ConQuest plausible-values workflow.
#' - `predict_mfrm_units()` and `sample_mfrm_plausible_values()` can score
#'   latent-regression fits under the fitted population model, but they require
#'   one-row-per-person background data for scored units when the fitted
#'   population model includes covariates. Intercept-only latent-regression
#'   fits (`population_formula = ~ 1`) can reconstruct that minimal person
#'   table internally during scoring.
#'
#' @section Latent-regression quick start:
#' For a first latent-regression run, keep the setup explicit:
#' 1. Put response data in `data`, with one row per rating event.
#' 2. Put background variables in `person_data`, with exactly one row per
#'    person. The ID column must match `person`, or be supplied through
#'    `person_id`.
#' 3. Use `method = "MML"` and a one-sided formula such as
#'    `population_formula = ~ Grade + Group`.
#' 4. Numeric/logical and factor/character predictors are expanded with
#'    `stats::model.matrix()`. After fitting, inspect
#'    `summary(fit)$population_coding` to see the fitted levels, contrasts, and
#'    encoded design columns that will be reused for scoring/replay.
#' 5. Start with `population_policy = "error"` while preparing data. Use
#'    `"omit"` only when complete-case removal is intended, and then inspect
#'    `summary(fit)$population_overview` and `summary(fit)$caveats` before
#'    reporting results.
#' 6. Report `summary(fit)$population_coefficients` as coefficients of the
#'    conditional-normal latent population model, not as a post hoc regression
#'    on EAP or MLE scores.
#'
#' Anchor inputs are optional:
#' - `anchors` should contain facet/level/fixed-value information.
#' - `group_anchors` should contain facet/level/group/group-value information.
#' Both are normalized internally, so column names can be flexible
#' (`facet`, `level`, `anchor`, `group`, `groupvalue`, etc.).
#'
#' Anchor audit behavior:
#' - `fit_mfrm()` runs an internal anchor audit.
#' - invalid rows are removed before estimation.
#' - duplicate rows keep the last occurrence for each key.
#' - `anchor_policy` controls whether detected issues are warned, treated as
#'   errors, or kept silent.
#'
#' Facet sign orientation:
#' - facets listed in `positive_facets` are treated as `+1`
#' - all other facets are treated as `-1`
#' This affects interpretation of reported facet measures.
#'
#' @section Performance tips:
#' For exploratory work, `method = "JML"` is usually faster than `method = "MML"`,
#' but it may require a larger `maxit` to converge on larger datasets.
#'
#' For MML runs, `quad_points` is the main accuracy/speed trade-off:
#' - `quad_points = 7` is a good lightweight default for quick iteration.
#' - `quad_points = 15` gives a more stable approximation for final reporting.
#' - `mml_engine = "direct"` remains the most stable general-purpose path.
#' - `mml_engine = "em"` or `"hybrid"` currently target `RSM` / `PCM` fits
#'   without a latent-regression population model.
#' - Benchmark your own workload before using `mml_engine = "em"` or
#'   `"hybrid"` for final reporting; `direct` remains the safer default when
#'   you have not compared engines for your data.
#'
#' Downstream diagnostics can also be staged:
#' - use `diagnose_mfrm(fit, residual_pca = "none")` for a quick first pass
#' - add residual PCA only when you need exploratory residual-structure evidence
#'
#' Downstream diagnostics report `ModelSE` / `RealSE` columns and related
#' reliability indices. For `MML`, non-person facet `ModelSE` values are based
#' on the observed information of the marginal log-likelihood and person rows
#' use posterior SDs from EAP scoring. For `JML`, these quantities remain
#' exploratory approximations and should not be treated as equally formal.
#'
#' For bounded `GPCM`, residual-based mean-square fit screens are also
#' best treated as exploratory diagnostics rather than strict Rasch-style
#' invariance tests, because the discrimination parameter is free.
#'
#' @section Interpreting output:
#' A typical first-pass read is:
#' 1. `fit$summary` for convergence and global fit indicators.
#' 2. `summary(fit)` for human-readable overviews.
#' 3. for `RSM` / `PCM`, `diagnose_mfrm(fit)` for element-level fit,
#'    approximate separation/reliability, and warning tables.
#' 4. for bounded `GPCM`, use [diagnose_mfrm()] and the residual-based
#'    table helpers as exploratory screens, together with posterior scoring /
#'    [compute_information()] where documented.
#'
#' @section Typical workflow:
#' 1. Fit the model with `fit_mfrm(...)`.
#' 2. Validate convergence and scale structure with `summary(fit)`.
#' 3. For `RSM` / `PCM`, run [diagnose_mfrm()] and proceed to reporting with
#'    [build_apa_outputs()].
#' 4. For bounded `GPCM`, use the fitted object, slope summary,
#'    [diagnose_mfrm()], residual-based table helpers, posterior scoring
#'    helpers, and [compute_information()] while broader downstream
#'    validation is still being completed. Use [gpcm_capability_matrix()] to
#'    confirm which helper families are currently supported, caveated, blocked,
#'    or deferred.
#'
#' @section References:
#' The ordered-category many-facet formulation follows Linacre (1989), with
#' the `RSM` and `PCM` branches grounded in Andrich (1978) and Masters (1982).
#' The bounded `GPCM` branch follows the generalized partial credit
#' formulation of Muraki (1992) under a package-specific positive
#' log-slope identification convention. The `MML` route follows the
#' quadrature-based marginal-likelihood framework of Bock and Aitkin (1981).
#'
#' - Andrich, D. (1978). *A rating formulation for ordered response
#'   categories*. Psychometrika, 43(4), 561-573.
#' - Bock, R. D., & Aitkin, M. (1981). *Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm*. Psychometrika, 46(4),
#'   443-459.
#' - Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#' - Masters, G. N. (1982). *A Rasch model for partial credit scoring*.
#'   Psychometrika, 47(2), 149-174.
#' - Muraki, E. (1992). *A generalized partial credit model: Application of an
#'   EM algorithm*. Applied Psychological Measurement, 16(2), 159-176.
#' - Robitzsch, A., & Steinfeld, J. (2018). *Modeling rater effects in
#'   achievements tests by item response models: Facets, generalized linear
#'   mixed models, or signal detection models?* Journal of Educational and
#'   Behavioral Statistics, 43(2), 218-244.
#'
#' @return
#' An object of class `mfrm_fit` (named list) with:
#' - `summary`: one-row model summary (`LogLik`, `AIC`, `BIC`, convergence)
#'   including public `Method`, internal `MethodUsed`, and
#'   `MMLEngineRequested`, `MMLEngineUsed`, and `EMIterations` for MML fits
#' - `facets$person`: person estimates (`Estimate`; plus `SD` for MML)
#' - `facets$others`: facet-level estimates for each facet
#' - `steps`: estimated threshold/step parameters
#' - `slopes`: estimated discrimination parameters for `GPCM` fits
#' - `population`: population-model metadata. Ordinary fits keep an inactive
#'   scaffold (`active = FALSE`, `posterior_basis = "legacy_mml"`). Active
#'   latent-regression fits store the fitted design matrix, regression
#'   coefficients, residual variance, omission audit, the complete-case
#'   estimation table (`person_table`), and the observed-person-aligned
#'   replay/export provenance table retained before complete-case omission
#'   (`person_table_replay`), plus stored categorical `xlevels` / `contrasts`
#'   for model-matrix replay and scoring, together with
#'   `posterior_basis = "population_model"`.
#' - `config`: resolved model configuration used for estimation
#'   (includes `config$anchor_audit`)
#' - `prep`: preprocessed data/level metadata
#' - `opt`: raw optimizer result from [stats::optim()]
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()], [build_apa_outputs()],
#'   [gpcm_capability_matrix], [mfrmr_workflow_methods],
#'   [mfrmr_reporting_and_apa]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#'
#' fit <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   model = "RSM",
#'   maxit = 25
#' )
#' fit$summary
#' s_fit <- summary(fit)
#' s_fit$overview[, c("Model", "Method", "Converged")]
#' p_fit <- plot(fit, draw = FALSE)
#' p_fit$wright_map$data$plot
#'
#' # MML is the default:
#' fit_mml <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   model = "RSM",
#'   quad_points = 7,
#'   maxit = 25
#' )
#' summary(fit_mml)
#'
#' # Latent regression (MML only) uses person-level background variables:
#' person_tbl <- unique(toy[c("Person")])
#' person_tbl$Grade <- seq_len(nrow(person_tbl))
#' person_tbl$Group <- rep(c("A", "B"), length.out = nrow(person_tbl))
#' \dontrun{
#' fit_pop <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   population_formula = ~ Grade + Group,
#'   person_data = person_tbl
#' )
#' summary(fit_pop)$population_overview
#' summary(fit_pop)$population_coding
#' }
#'
#' # Binary responses are supported as ordered two-category scores:
#' set.seed(1)
#' binary_toy <- expand.grid(
#'   Person = paste0("P", 1:30),
#'   Item = paste0("I", 1:4),
#'   stringsAsFactors = FALSE
#' )
#' theta <- stats::rnorm(length(unique(binary_toy$Person)))
#' beta <- seq(-0.8, 0.8, length.out = length(unique(binary_toy$Item)))
#' eta <- theta[match(binary_toy$Person, unique(binary_toy$Person))] -
#'   beta[match(binary_toy$Item, unique(binary_toy$Item))]
#' binary_toy$Score <- stats::rbinom(nrow(binary_toy), 1, stats::plogis(eta))
#' fit_binary <- fit_mfrm(
#'   data = binary_toy,
#'   person = "Person",
#'   facets = "Item",
#'   score = "Score",
#'   model = "RSM",
#'   method = "JML",
#'   maxit = 50
#' )
#' fit_binary$summary[, c("Model", "Categories", "Converged")]
#'
#' # Next steps after fitting:
#' diag_mml <- diagnose_mfrm(fit_mml, residual_pca = "none")
#' chk <- reporting_checklist(fit_mml, diagnostics = diag_mml)
#' head(chk$checklist[, c("Section", "Item", "DraftReady")])
#' }
#' @export
fit_mfrm <- function(data,
                     person,
                     facets,
                     score,
                     rating_min = NULL,
                     rating_max = NULL,
                     weight = NULL,
                     keep_original = FALSE,
                     model = c("RSM", "PCM", "GPCM"),
                     method = c("MML", "JML", "JMLE"),
                     step_facet = NULL,
                     slope_facet = NULL,
                     anchors = NULL,
                     group_anchors = NULL,
                     noncenter_facet = "Person",
                     dummy_facets = NULL,
                     positive_facets = NULL,
                     anchor_policy = c("warn", "error", "silent"),
                     min_common_anchors = 5L,
                     min_obs_per_element = 30,
                     min_obs_per_category = 10,
                     quad_points = 15,
                     maxit = 400,
                     reltol = 1e-6,
                     mml_engine = c("direct", "em", "hybrid"),
                     population_formula = NULL,
                     person_data = NULL,
                     person_id = NULL,
                     population_policy = c("error", "omit")) {
  # -- input validation --
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame. Got: ", class(data)[1], ". ",
         "Convert with as.data.frame() if needed.", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("`data` has zero rows. ",
         "Supply a data.frame with at least one observation.", call. = FALSE)
  }
  if (!is.character(person) || length(person) != 1 || !nzchar(person)) {
    stop("`person` must be a single non-empty character string ",
         "naming the person column.", call. = FALSE)
  }
  if (!is.character(facets) || length(facets) == 0) {
    stop("`facets` must be a character vector of one or more facet column names.",
         call. = FALSE)
  }
  if (!is.character(score) || length(score) != 1 || !nzchar(score)) {
    stop("`score` must be a single non-empty character string ",
         "naming the score column.", call. = FALSE)
  }
  if (!is.null(weight) && (!is.character(weight) || length(weight) != 1)) {
    stop("`weight` must be NULL or a single character string ",
         "naming the weight column.", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1) {
    stop("`maxit` must be a positive integer. Got: ", deparse(maxit), ".",
         call. = FALSE)
  }
  if (!is.numeric(reltol) || length(reltol) != 1 || reltol <= 0) {
    stop("`reltol` must be a positive number. Got: ", deparse(reltol), ".",
         call. = FALSE)
  }
  if (!is.numeric(quad_points) || length(quad_points) != 1 || quad_points < 1) {
    stop("`quad_points` must be a positive integer. Got: ", deparse(quad_points), ".",
         call. = FALSE)
  }
  if (!is.null(person_id) && (!is.character(person_id) || length(person_id) != 1 || !nzchar(person_id))) {
    stop("`person_id` must be NULL or a single non-empty character string naming the person column in `person_data`.",
         call. = FALSE)
  }
  if (!is.null(slope_facet) && (!is.character(slope_facet) || length(slope_facet) != 1 || !nzchar(slope_facet))) {
    stop("`slope_facet` must be NULL or a single non-empty character string naming a facet column.",
         call. = FALSE)
  }

  model <- toupper(match.arg(model))
  method_input <- toupper(match.arg(method))
  method <- ifelse(method_input == "JML", "JMLE", method_input)
  mml_engine <- tolower(match.arg(mml_engine))
  anchor_policy <- tolower(match.arg(anchor_policy))
  population_policy <- tolower(match.arg(population_policy))

  population <- prepare_mfrm_population_scaffold(
    data = data,
    person = person,
    population_formula = population_formula,
    person_data = person_data,
    person_id = person_id,
    population_policy = population_policy
  )
  if (isTRUE(population$active)) {
    if (!identical(method_input, "MML")) {
      stop("Latent-regression scaffolding currently requires `method = 'MML'`. ",
           "The requested population model can currently be estimated only in the MML branch.",
           call. = FALSE)
    }
  }
  estimation_data <- data
  if (isTRUE(population$active) && length(population$included_persons) > 0) {
    estimation_mask <- as.character(data[[person]]) %in% population$included_persons
    estimation_data <- data[estimation_mask, , drop = FALSE]
  }

  anchor_audit <- audit_mfrm_anchors(
    data = estimation_data,
    person = person,
    facets = facets,
    score = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight = weight,
    keep_original = keep_original,
    anchors = anchors,
    group_anchors = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  anchors <- anchor_audit$anchors
  group_anchors <- anchor_audit$group_anchors

  issue_counts <- anchor_audit$issue_counts
  issue_total <- if (is.null(issue_counts) || nrow(issue_counts) == 0) 0L else sum(issue_counts$N, na.rm = TRUE)
  if (issue_total > 0) {
    msg <- format_anchor_audit_message(anchor_audit)
    if (anchor_policy == "error") {
      stop(msg, call. = FALSE)
    } else if (anchor_policy == "warn") {
      warning(msg, call. = FALSE)
    }
  }

  fit <- mfrm_estimate(
    data = estimation_data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    model = model,
    method = method,
    step_facet = step_facet,
    slope_facet = slope_facet,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    population = population,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol,
    mml_engine = mml_engine
  )

  fit$config$anchor_audit <- anchor_audit
  fit$config$method_input <- method_input
  fit$population <- finalize_mfrm_population_fit(fit, population)
  fit$config$population_spec <- compact_population_spec(fit$population, fit$prep$levels$Person)
  fit$config$population_active <- isTRUE(population$active)
  fit$config$posterior_basis <- as.character(fit$population$posterior_basis %||% "legacy_mml")
  fit$config$population_policy <- fit$population$policy %||% NULL
  fit$config$population_formula <- if (!is.null(fit$population$formula)) {
    paste(deparse(fit$population$formula), collapse = " ")
  } else {
    NULL
  }

  class(fit) <- c("mfrm_fit", class(fit))
  fit
}

finalize_mfrm_population_fit <- function(fit, population) {
  pop <- population %||% list()
  if (!isTRUE(pop$active)) {
    return(pop)
  }

  sizes <- build_param_sizes(fit$config)
  params <- expand_params(fit$opt$par, sizes, fit$config)
  coeff <- as.numeric(params$population$coefficients %||% numeric(0))
  coeff_names <- names(params$population$coefficients %||% NULL)
  if (!is.null(coeff_names) && length(coeff_names) == length(coeff)) {
    names(coeff) <- coeff_names
  }

  notes <- as.character(pop$notes %||% character(0))
  notes <- c(
    notes,
    "Latent-regression parameters were estimated under the MML population-model branch."
  )

  pop$coefficients <- coeff
  pop$sigma2 <- as.numeric(params$population$sigma2[1] %||% NA_real_)
  pop$converged <- isTRUE(fit$summary$Converged[1])
  pop$logLik_component <- as.numeric(fit$summary$LogLik[1] %||% NA_real_)
  pop$posterior_basis <- "population_model"
  pop$design_columns <- pop$design_columns %||% names(coeff)
  pop$notes <- unique(notes)
  pop
}

prepare_mfrm_population_scaffold <- function(data,
                                             person,
                                             population_formula = NULL,
                                             person_data = NULL,
                                             person_id = NULL,
                                             population_policy = c("error", "omit"),
                                             population_xlevels = NULL,
                                             population_contrasts = NULL,
                                             require_full_rank = TRUE) {
  if (is.null(population_formula)) {
    return(list(
      active = FALSE,
      formula = NULL,
      person_id = if (is.null(person_id)) person else person_id,
      person_table = NULL,
      person_table_replay = NULL,
      person_table_replay_scope = NULL,
      design_matrix = NULL,
      design_columns = NULL,
      xlevels = NULL,
      contrasts = NULL,
      coefficients = NULL,
      sigma2 = NULL,
      converged = FALSE,
      posterior_basis = "legacy_mml",
      policy = NULL,
      included_persons = character(0),
      omitted_persons = character(0),
      response_rows_retained = nrow(data),
      response_rows_omitted = 0L,
      notes = "No population model was requested; this fit uses the package's legacy unconditional estimation path."
    ))
  }

  population_policy <- match.arg(population_policy)
  population_formula <- tryCatch(
    stats::as.formula(population_formula),
    error = function(e) {
      stop("`population_formula` must be coercible to a one-sided formula such as `~ grade + ses`. ",
           "Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (length(population_formula) != 2L) {
    stop("`population_formula` must be a one-sided formula such as `~ grade + ses`.",
         call. = FALSE)
  }
  if (!is.data.frame(person_data)) {
    stop("`person_data` must be a data.frame with one row per person when `population_formula` is supplied.",
         call. = FALSE)
  }
  if (is.null(person_id)) {
    if (!person %in% names(person_data)) {
      stop("`person_data` must contain the person column `", person,
           "` or you must supply `person_id` explicitly when `population_formula` is used.",
           call. = FALSE)
    }
    person_id <- person
  }
  if (!person_id %in% names(person_data)) {
    stop("`person_id = '", person_id, "'` is not a column in `person_data`.",
         call. = FALSE)
  }

  observed_persons <- unique(as.character(data[[person]]))
  person_tbl <- as.data.frame(person_data, stringsAsFactors = FALSE)
  person_tbl[[person_id]] <- as.character(person_tbl[[person_id]])
  if (anyNA(person_tbl[[person_id]]) || any(!nzchar(person_tbl[[person_id]]))) {
    stop("`person_data` contains missing or empty person IDs in column `", person_id, "`.",
         call. = FALSE)
  }
  dup_ids <- unique(person_tbl[[person_id]][duplicated(person_tbl[[person_id]])])
  if (length(dup_ids) > 0) {
    preview <- paste(utils::head(dup_ids, 5), collapse = ", ")
    stop("`person_data` must contain one unique row per person. Duplicate IDs detected in `",
         person_id, "`: ", preview, if (length(dup_ids) > 5) ", ..." else ".", call. = FALSE)
  }
  missing_persons <- setdiff(observed_persons, person_tbl[[person_id]])
  if (length(missing_persons) > 0) {
    preview <- paste(utils::head(missing_persons, 5), collapse = ", ")
    stop("`person_data` is missing ", length(missing_persons), " person(s) observed in `data`: ",
         preview, if (length(missing_persons) > 5) ", ..." else ".", call. = FALSE)
  }

  keep_idx <- match(observed_persons, person_tbl[[person_id]])
  person_tbl <- person_tbl[keep_idx, , drop = FALSE]
  rownames(person_tbl) <- NULL
  person_tbl_replay <- person_tbl

  model_frame_args <- list(
    formula = population_formula,
    data = person_tbl,
    na.action = stats::na.pass
  )
  if (!is.null(population_xlevels)) {
    model_frame_args$xlev <- population_xlevels
  }
  mf <- tryCatch(
    do.call(stats::model.frame, model_frame_args),
    error = function(e) {
      stop("Could not build the latent-regression model frame from `person_data`. ",
           "Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  mf_names <- names(mf)
  if (length(mf_names) > 0) {
    unsupported <- mf_names[vapply(mf[mf_names], function(x) is.list(x) || is.data.frame(x), logical(1))]
    if (length(unsupported) > 0) {
      stop("Variables referenced in `population_formula` must be atomic columns usable by stats::model.matrix(). ",
           "Numeric, logical, factor, ordered factor, and character predictors are supported; unsupported columns: ",
           paste(unsupported, collapse = ", "), ".", call. = FALSE)
    }
  }

  complete_mask <- if (ncol(mf) == 0) rep(TRUE, nrow(person_tbl)) else stats::complete.cases(mf)
  omitted_persons <- person_tbl[[person_id]][!complete_mask]
  if (length(omitted_persons) > 0 && identical(population_policy, "error")) {
    preview <- paste(utils::head(omitted_persons, 5), collapse = ", ")
    stop("`person_data` contains missing covariate values for persons referenced by `population_formula`. ",
         "Use `population_policy = 'omit'` to build a complete-case scaffold instead. Affected IDs: ",
         preview, if (length(omitted_persons) > 5) ", ..." else ".", call. = FALSE)
  }
  if (length(omitted_persons) > 0) {
    person_tbl <- person_tbl[complete_mask, , drop = FALSE]
    mf <- mf[complete_mask, , drop = FALSE]
    rownames(person_tbl) <- NULL
    rownames(mf) <- NULL
  }
  if (nrow(person_tbl) == 0) {
    stop("No persons remain in `person_data` after applying the latent-regression covariate policy.",
         call. = FALSE)
  }

  included_persons <- as.character(person_tbl[[person_id]])
  response_keep <- as.character(data[[person]]) %in% included_persons
  response_rows_retained <- sum(response_keep)
  response_rows_omitted <- sum(!response_keep)

  terms_obj <- attr(mf, "terms") %||% stats::terms(population_formula)
  mm <- tryCatch(
    stats::model.matrix(terms_obj, data = mf, contrasts.arg = population_contrasts),
    error = function(e) {
      stop("Could not build the latent-regression design matrix from `person_data`. ",
           "Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  xlevels <- tryCatch(stats::.getXlevels(terms_obj, mf), error = function(e) NULL)
  mm_contrasts <- attr(mm, "contrasts", exact = TRUE)
  qr_mm <- qr(mm)
  if (isTRUE(require_full_rank) && qr_mm$rank < ncol(mm)) {
    stop("The latent-regression design matrix is rank-deficient. ",
         "Adjust `population_formula` or the coding of `person_data` before fitting.",
         call. = FALSE)
  }

  notes <- c(
    "Population-model covariate scaffolding is active.",
    "This scaffold validates the person-level covariates and design-matrix contract used by the latent-regression branch."
  )
  if (length(omitted_persons) > 0) {
    notes <- c(
      notes,
      paste0("Complete-case scaffolding retained ", nrow(person_tbl), " person(s) and omitted ",
             length(omitted_persons), " person(s) under `population_policy = 'omit'`.")
    )
  }
  if (response_rows_omitted > 0) {
    notes <- c(
      notes,
      paste0("If this scaffold is activated in estimation, ", response_rows_omitted,
             " response row(s) would be excluded because their persons lack complete background data.")
    )
  }

  list(
    active = TRUE,
    formula = population_formula,
    person_id = person_id,
    person_table = person_tbl,
    person_table_replay = person_tbl_replay,
    person_table_replay_scope = "observed_person_subset_pre_omit",
    design_matrix = mm,
    design_columns = colnames(mm),
    xlevels = xlevels,
    contrasts = mm_contrasts,
    coefficients = NULL,
    sigma2 = NULL,
    converged = FALSE,
    posterior_basis = "population_model",
    policy = population_policy,
    included_persons = included_persons,
    omitted_persons = as.character(omitted_persons),
    response_rows_retained = as.integer(response_rows_retained),
    response_rows_omitted = as.integer(response_rows_omitted),
    notes = notes
  )
}

format_anchor_audit_message <- function(anchor_audit) {
  if (is.null(anchor_audit$issue_counts) || nrow(anchor_audit$issue_counts) == 0) {
    return("Anchor audit detected no issues.")
  }
  nonzero <- anchor_audit$issue_counts |>
    dplyr::filter(.data$N > 0)
  if (nrow(nonzero) == 0) {
    return("Anchor audit detected no issues.")
  }
  labels <- paste0(nonzero$Issue, "=", nonzero$N)
  paste0(
    "Anchor audit detected ", sum(nonzero$N), " issue row(s): ",
    paste(labels, collapse = "; "),
    ". Invalid rows were removed; duplicate keys keep the last row."
  )
}

summarize_linkage_by_facet <- function(df, facet) {
  by_level <- df |>
    dplyr::group_by(.data[[facet]]) |>
    dplyr::summarize(
      PersonsPerLevel = dplyr::n_distinct(.data$Person),
      Observations = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    )

  by_person <- df |>
    dplyr::group_by(.data$Person) |>
    dplyr::summarize(
      LevelsPerPerson = dplyr::n_distinct(.data[[facet]]),
      .groups = "drop"
    )

  tibble::tibble(
    Facet = facet,
    Levels = nrow(by_level),
    MinPersonsPerLevel = min(by_level$PersonsPerLevel, na.rm = TRUE),
    MedianPersonsPerLevel = stats::median(by_level$PersonsPerLevel, na.rm = TRUE),
    MinLevelsPerPerson = min(by_person$LevelsPerPerson, na.rm = TRUE),
    MedianLevelsPerPerson = stats::median(by_person$LevelsPerPerson, na.rm = TRUE)
  )
}

normalize_compare_signature <- function(fit) {
  cfg <- fit$config %||% list()
  anchor_tables <- extract_anchor_tables(cfg)

  anchor_tbl <- as.data.frame(anchor_tables$anchors %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(anchor_tbl) > 0) {
    anchor_tbl <- anchor_tbl[, intersect(c("Facet", "Level", "Anchor"), names(anchor_tbl)), drop = FALSE]
    anchor_tbl <- anchor_tbl[do.call(order, c(anchor_tbl[intersect(c("Facet", "Level"), names(anchor_tbl))], na.last = TRUE)), , drop = FALSE]
  }

  group_tbl <- as.data.frame(anchor_tables$groups %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(group_tbl) > 0) {
    group_tbl <- group_tbl[, intersect(c("Facet", "Level", "Group", "GroupValue"), names(group_tbl)), drop = FALSE]
    group_tbl <- group_tbl[do.call(order, c(group_tbl[intersect(c("Facet", "Level", "Group"), names(group_tbl))], na.last = TRUE)), , drop = FALSE]
  }

  list(
    model = as.character(cfg$model %||% NA_character_),
    method = as.character(cfg$method %||% NA_character_),
    person = as.character(cfg$person_col %||% NA_character_),
    facets = sort(as.character(cfg$facet_cols %||% fit$prep$facet_names %||% character(0))),
    score = as.character(cfg$score_col %||% NA_character_),
    weight = as.character(cfg$weight_col %||% NA_character_),
    step_facet = as.character(cfg$step_facet %||% NA_character_),
    slope_facet = as.character(cfg$slope_facet %||% NA_character_),
    noncenter_facet = as.character(cfg$noncenter_facet %||% "Person"),
    dummy_facets = sort(as.character(cfg$dummy_facets %||% character(0))),
    positive_facets = sort(as.character(cfg$positive_facets %||% character(0))),
    anchors = anchor_tbl,
    group_anchors = group_tbl
  )
}

same_signature_component <- function(x, y) {
  if (is.data.frame(x) || is.data.frame(y)) {
    return(identical(
      as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE),
      as.data.frame(y %||% data.frame(), stringsAsFactors = FALSE)
    ))
  }
  identical(x, y)
}

audit_compare_mfrm_nesting <- function(fits, labels) {
  if (length(fits) != 2L) {
    return(list(
      eligible = FALSE,
      reason = "Likelihood-ratio tests are only audited for two-model comparisons.",
      simpler = NA_character_,
      complex = NA_character_,
      relation = "unsupported"
    ))
  }

  sigs <- lapply(fits, normalize_compare_signature)
  lbls <- as.character(labels)

  same_components <- c(
    person = same_signature_component(sigs[[1]]$person, sigs[[2]]$person),
    facets = same_signature_component(sigs[[1]]$facets, sigs[[2]]$facets),
    score = same_signature_component(sigs[[1]]$score, sigs[[2]]$score),
    weight = same_signature_component(sigs[[1]]$weight, sigs[[2]]$weight),
    noncenter_facet = same_signature_component(sigs[[1]]$noncenter_facet, sigs[[2]]$noncenter_facet),
    dummy_facets = same_signature_component(sigs[[1]]$dummy_facets, sigs[[2]]$dummy_facets),
    positive_facets = same_signature_component(sigs[[1]]$positive_facets, sigs[[2]]$positive_facets),
    anchors = same_signature_component(sigs[[1]]$anchors, sigs[[2]]$anchors),
    group_anchors = same_signature_component(sigs[[1]]$group_anchors, sigs[[2]]$group_anchors)
  )

  if (!all(same_components)) {
    mismatch <- names(same_components)[!same_components]
    return(list(
      eligible = FALSE,
      reason = paste0(
        "Models differ in structural comparison settings: ",
        paste(mismatch, collapse = ", "),
        "."
      ),
      simpler = NA_character_,
      complex = NA_character_,
      relation = "unsupported"
    ))
  }

  model_pair <- toupper(c(sigs[[1]]$model, sigs[[2]]$model))
  if (identical(model_pair[1], model_pair[2])) {
    return(list(
      eligible = FALSE,
      reason = "Both fits use the same model family, so there is no supported nested restriction to test.",
      simpler = lbls[1],
      complex = lbls[2],
      relation = "same_model"
    ))
  }

  if (setequal(model_pair, c("RSM", "PCM"))) {
    idx_simple <- which(model_pair == "RSM")[1]
    idx_complex <- which(model_pair == "PCM")[1]
    step_facet <- sigs[[idx_complex]]$step_facet
    if (!is.na(step_facet) && nzchar(step_facet)) {
      return(list(
        eligible = TRUE,
        reason = paste0(
          "Supported nesting audit passed: shared design/constraints with RSM nested inside PCM on step facet '",
          step_facet,
          "'."
        ),
        simpler = lbls[idx_simple],
        complex = lbls[idx_complex],
        relation = "RSM_in_PCM"
      ))
    }
  }

  list(
    eligible = FALSE,
    reason = "Automatic nesting audit currently supports only RSM nested inside PCM under shared design and constraints.",
    simpler = NA_character_,
    complex = NA_character_,
    relation = "unsupported"
  )
}

#' Summarize MFRM input data (TAM-style descriptive snapshot)
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value. Supply with
#'   `rating_max` to retain unused boundary categories in the intended score
#'   support.
#' @param rating_max Optional maximum category value. Supply with
#'   `rating_min` to retain unused boundary categories in the intended score
#'   support.
#' @param keep_original Keep original category values. Use this with
#'   `rating_min` / `rating_max` when the intended scale has unused
#'   intermediate categories such as `1, 2, 4, 5` on a 1-5 scale.
#' @param include_person_facet If `TRUE`, include person-level rows in
#'   `facet_level_summary`.
#' @param include_agreement If `TRUE`, include an observed-score inter-rater
#'   agreement bundle (summary/pairs/settings) in the output.
#' @param rater_facet Optional rater facet name used for agreement summaries.
#'   If `NULL`, inferred from facet names.
#' @param context_facets Optional facets used to define matched contexts for
#'   agreement. If `NULL`, all remaining facets (including `Person`) are used.
#' @param agreement_top_n Optional maximum number of agreement pair rows.
#'
#' @details
#' This function provides a compact descriptive bundle similar to the
#' pre-fit summaries commonly checked in TAM workflows:
#' sample size, score distribution, per-facet coverage, and linkage counts.
#' `psych::describe()` is used for numeric descriptives of score and weight.
#'
#' **Key data-quality checks to perform before fitting:**
#' - *Sparse categories*: any score category with fewer than 10 weighted
#'   observations may produce unstable threshold estimates
#'   (Linacre, 2002).  Consider collapsing adjacent categories.
#' - *Unlinked elements*: if a facet level has zero overlap with one or
#'   more levels of another facet, the design is disconnected and
#'   parameters cannot be placed on a common scale.  Check
#'   `linkage_summary` for low connectivity.
#' - *Extreme scores*: persons or facet levels with all-minimum or
#'   all-maximum scores yield infinite logit estimates under JML;
#'   they are handled via Bayesian shrinkage under MML.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `overview`: confirms sample size, facet count, and category span.
#'   The `MinWeightedN` column shows the smallest weighted observation
#'   count across all facet levels; values below 30 may lead to
#'   unstable parameter estimates.
#' - `missing_by_column`: identifies immediate data-quality risks.
#'   Any non-zero count warrants investigation before fitting.
#' - `score_distribution`: checks sparse/unused score categories.
#'   Balanced usage across categories is ideal; heavily skewed
#'   distributions may compress the measurement range.
#' - `facet_level_summary` and `linkage_summary`: checks per-level
#'   support and person-facet connectivity.  Low linkage ratios
#'   indicate sparse or disconnected design blocks.
#' - `agreement`: optional observed inter-rater consistency summary
#'   (exact agreement, correlation, mean differences per rater pair).
#'
#' @section Typical workflow:
#' 1. Run `describe_mfrm_data()` on long-format input.
#' 2. Review `summary(ds)` and `plot(ds, ...)`.
#' 3. Resolve missingness/sparsity issues before [fit_mfrm()].
#'
#' @return A list of class `mfrm_data_description` with:
#' - `overview`: one-row run-level summary
#' - `missing_by_column`: missing counts in selected input columns
#' - `score_descriptives`: output from [psych::describe()] for score
#' - `weight_descriptives`: output from [psych::describe()] for weight
#' - `score_distribution`: weighted and raw score frequencies over the prepared
#'   score support. Unused boundary categories are retained when the rating
#'   range was supplied explicitly; unused intermediate categories require
#'   `keep_original = TRUE`.
#' - `facet_level_summary`: per-level usage and score summaries
#' - `linkage_summary`: person-facet connectivity diagnostics
#' - `agreement`: observed-score inter-rater agreement bundle
#' - `score_support`: minimal prepared score-support metadata used by
#'   `summary(ds)$caveats`
#'
#' @seealso [fit_mfrm()], [audit_mfrm_anchors()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score"
#' )
#' s_ds <- summary(ds)
#' s_ds$overview
#' p_ds <- plot(ds, draw = FALSE)
#' p_ds$data$plot
#' @export
describe_mfrm_data <- function(data,
                               person,
                               facets,
                               score,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               include_person_facet = FALSE,
                               include_agreement = TRUE,
                               rater_facet = NULL,
                               context_facets = NULL,
                               agreement_top_n = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original
  )

  df <- prep$data |>
    dplyr::mutate(
      Person = as.character(.data$Person),
      dplyr::across(dplyr::all_of(prep$facet_names), as.character)
    )

  selected_cols <- unique(c(person, facets, score, if (!is.null(weight)) weight))
  missing_by_column <- tibble::tibble(
    Column = selected_cols,
    Missing = vapply(selected_cols, function(col) sum(is.na(data[[col]])), integer(1))
  )

  score_desc <- psych::describe(df$Score, fast = TRUE)
  weight_desc <- psych::describe(df$Weight, fast = TRUE)

  total_weight <- sum(df$Weight, na.rm = TRUE)
  observed_score_distribution <- df |>
    dplyr::group_by(.data$Score) |>
    dplyr::summarize(
      RawN = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      Percent = ifelse(total_weight > 0, 100 * .data$WeightedN / total_weight, NA_real_),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Score)
  score_distribution <- tibble::tibble(Score = seq(prep$rating_min, prep$rating_max)) |>
    dplyr::left_join(observed_score_distribution, by = "Score") |>
    dplyr::mutate(
      RawN = tidyr::replace_na(.data$RawN, 0L),
      WeightedN = tidyr::replace_na(.data$WeightedN, 0),
      Percent = tidyr::replace_na(.data$Percent, 0)
    )

  report_facets <- prep$facet_names
  if (isTRUE(include_person_facet)) {
    report_facets <- c("Person", report_facets)
  }

  facet_level_summary <- purrr::map_dfr(report_facets, function(facet) {
    df |>
      dplyr::group_by(.data[[facet]]) |>
      dplyr::summarize(
        RawN = dplyr::n(),
        WeightedN = sum(.data$Weight, na.rm = TRUE),
        MeanScore = weighted_mean(.data$Score, .data$Weight),
        SDScore = stats::sd(.data$Score, na.rm = TRUE),
        MinScore = min(.data$Score, na.rm = TRUE),
        MaxScore = max(.data$Score, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::rename(Level = dplyr::all_of(facet)) |>
      dplyr::mutate(
        Facet = facet,
        Level = as.character(.data$Level),
        .before = 1
      )
  })

  linkage_summary <- if (length(prep$facet_names) == 0) {
    tibble::tibble()
  } else {
    purrr::map_dfr(prep$facet_names, function(facet) summarize_linkage_by_facet(df, facet))
  }

  agreement_bundle <- list(
    summary = data.frame(),
    pairs = data.frame(),
    settings = list(
      included = FALSE,
      rater_facet = NA_character_,
      context_facets = character(0),
      expected_exact_from_model = FALSE,
      top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n))
    )
  )

  if (isTRUE(include_agreement) && length(prep$facet_names) > 0) {
    known_facets <- c("Person", prep$facet_names)
    if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
      rater_facet <- infer_default_rater_facet(prep$facet_names)
    } else {
      rater_facet <- as.character(rater_facet[1])
    }
    if (is.null(rater_facet) || !rater_facet %in% known_facets) {
      stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
    }
    if (identical(rater_facet, "Person")) {
      stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
    }

    if (is.null(context_facets)) {
      facet_cols <- known_facets
      resolved_context <- setdiff(facet_cols, rater_facet)
    } else {
      context_facets <- unique(as.character(context_facets))
      unknown <- setdiff(context_facets, known_facets)
      if (length(unknown) > 0) {
        stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "))
      }
      resolved_context <- setdiff(context_facets, rater_facet)
      if (length(resolved_context) == 0) {
        stop("`context_facets` must include at least one facet different from `rater_facet`.")
      }
      facet_cols <- c(rater_facet, resolved_context)
    }

    obs_agreement <- df |>
      dplyr::select(dplyr::all_of(unique(c("Person", prep$facet_names, "Score", "Weight")))) |>
      dplyr::rename(Observed = "Score")

    agreement <- calc_interrater_agreement(
      obs_df = obs_agreement,
      facet_cols = facet_cols,
      rater_facet = rater_facet,
      res = NULL
    )
    agreement_pairs <- as.data.frame(agreement$pairs, stringsAsFactors = FALSE)
    if (!is.null(agreement_top_n) && nrow(agreement_pairs) > 0) {
      agreement_pairs <- agreement_pairs |>
        dplyr::slice_head(n = max(1L, as.integer(agreement_top_n)))
    }

    agreement_bundle <- list(
      summary = as.data.frame(agreement$summary, stringsAsFactors = FALSE),
      pairs = agreement_pairs,
      settings = list(
        included = TRUE,
        rater_facet = rater_facet,
        context_facets = resolved_context,
        expected_exact_from_model = FALSE,
        top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n))
      )
    )
  }

  overview <- tibble::tibble(
    Observations = nrow(df),
    TotalWeight = total_weight,
    Persons = length(prep$levels$Person),
    Facets = length(prep$facet_names),
    Categories = prep$rating_max - prep$rating_min + 1,
    RatingMin = prep$rating_min,
    RatingMax = prep$rating_max
  )

  out <- list(
    overview = overview,
    missing_by_column = missing_by_column,
    score_descriptives = score_desc,
    weight_descriptives = weight_desc,
    score_distribution = score_distribution,
    facet_level_summary = facet_level_summary,
    linkage_summary = linkage_summary,
    agreement = agreement_bundle,
    score_support = list(
      data = data.frame(Score = sort(unique(df$Score))),
      rating_min = prep$rating_min,
      rating_max = prep$rating_max,
      unused_score_categories = prep$unused_score_categories,
      score_map = prep$score_map
    )
  )
  class(out) <- c("mfrm_data_description", class(out))
  out
}

#' @export
print.mfrm_data_description <- function(x, ...) {
  cat("mfrm data description\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    print(x$overview, row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(x$score_distribution, row.names = FALSE)
  }
  if (!is.null(x$agreement$summary) && nrow(x$agreement$summary) > 0) {
    cat("\nInter-rater agreement (observed)\n")
    print(x$agreement$summary, row.names = FALSE)
  }
  invisible(x)
}

#' Summarize a data-description object
#'
#' @param object Output from [describe_mfrm_data()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in preview blocks.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is intended as a compact pre-fit quality snapshot for
#' manuscripts and analysis logs.
#'
#' @section Interpreting output:
#' Recommended read order:
#' - `overview`: sample size, persons/facets/categories.
#' - `missing`: missingness hotspots by selected input columns.
#' - `score_distribution`: category usage balance.
#' - `notes` / printed `Caveats`: retained zero-count score categories and
#'   related score-support caveats; intermediate unused categories should be
#'   treated as threshold-functioning warnings before model fitting.
#' - `facet_overview`: coverage per facet (minimum/maximum weighted counts).
#' - `agreement`: observed-score inter-rater agreement (when available).
#'
#' Very low `MinWeightedN` in `facet_overview` is a practical warning for
#' unstable downstream facet estimates.
#'
#' @section Typical workflow:
#' 1. Run [describe_mfrm_data()] on raw long-format data.
#' 2. Inspect `summary(ds)` before model fitting.
#' 3. Resolve sparse/missing issues, then run [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_data_description`.
#' - `overview`: design/sample counts
#' - `missing`: top columns by missingness
#' - `score_distribution`: compact score-usage table, including zero-count
#'   categories retained by the prepared score support
#' - `facet_overview`: facet-level coverage summary
#' - `agreement`: inter-rater agreement summary when available
#' - `reporting_map`: manuscript-oriented guide to what is covered here versus
#'   which companion outputs should be consulted
#' - `caveats`: structured warning/review rows for score-support issues;
#'   `print(summary(ds))` shows a compact `Caveats` block when rows are present
#' @seealso [describe_mfrm_data()], [summary.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(ds)
#' @export
summary.mfrm_data_description <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
  missing_tbl <- as.data.frame(object$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(missing_tbl) > 0 && all(c("Column", "Missing") %in% names(missing_tbl))) {
    missing_tbl <- missing_tbl |>
      dplyr::arrange(dplyr::desc(.data$Missing), .data$Column) |>
      dplyr::slice_head(n = top_n)
  }

  score_dist_full <- as.data.frame(object$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
  score_dist <- score_dist_full
  if (nrow(score_dist) > 0) {
    score_dist <- utils::head(score_dist, n = top_n)
  }

  facet_tbl <- as.data.frame(object$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
  facet_overview <- data.frame()
  if (nrow(facet_tbl) > 0 && all(c("Facet", "Level", "WeightedN") %in% names(facet_tbl))) {
    facet_overview <- facet_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n_distinct(.data$Level),
        TotalWeightedN = sum(.data$WeightedN, na.rm = TRUE),
        MeanWeightedN = mean(.data$WeightedN, na.rm = TRUE),
        MinWeightedN = min(.data$WeightedN, na.rm = TRUE),
        MaxWeightedN = max(.data$WeightedN, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Facet) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  agreement_tbl <- as.data.frame(object$agreement$summary %||% data.frame(), stringsAsFactors = FALSE)
  reporting_map <- data.frame(
    Area = c(
      "Sample / design counts",
      "Missingness audit",
      "Score usage / category distribution",
      "Facet coverage",
      "Inter-rater agreement",
      "Fit / reliability / residual PCA"
    ),
    CoveredHere = c("yes", "yes", "yes", "yes", if (nrow(agreement_tbl) > 0) "yes" else "partial", "no"),
    CompanionOutput = c(
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...)) / plot_interrater_agreement()",
      "summary(diagnose_mfrm(fit))"
    ),
    stringsAsFactors = FALSE
  )
  notes <- character(0)
  if (nrow(missing_tbl) > 0 && any(suppressWarnings(as.numeric(missing_tbl$Missing)) > 0, na.rm = TRUE)) {
    notes <- c(notes, "Missing values were detected in one or more input columns.")
  } else {
    notes <- c(notes, "No missing values were detected in selected input columns.")
  }
  caveat_prep <- object$score_support %||% object$prep %||% NULL
  caveats <- collect_mfrm_caveats(
    prep = caveat_prep,
    score_distribution = score_dist_full,
    include_recode = TRUE,
    context = "data"
  )
  if (nrow(caveats) > 0 && "Message" %in% names(caveats)) {
    notes <- c(notes, as.character(caveats$Message))
  }

  out <- list(
    overview = overview,
    missing = missing_tbl,
    score_distribution = score_dist,
    facet_overview = facet_overview,
    agreement = agreement_tbl,
    reporting_map = reporting_map,
    caveats = caveats,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_data_description"
  out
}

#' @export
print.summary.mfrm_data_description <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("mfrm Data Description Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$missing) && nrow(x$missing) > 0) {
    cat("\nMissing by column\n")
    print(round_numeric_df(as.data.frame(x$missing), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(round_numeric_df(as.data.frame(x$score_distribution), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_overview) && nrow(x$facet_overview) > 0) {
    cat("\nFacet coverage\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$agreement) && nrow(x$agreement) > 0) {
    cat("\nInter-rater agreement\n")
    print(round_numeric_df(as.data.frame(x$agreement), digits = digits), row.names = FALSE)
  }
  print_caveat_section(x$caveats)
  if (!is.null(x$reporting_map) && nrow(x$reporting_map) > 0) {
    cat("\nPaper reporting map\n")
    print(as.data.frame(x$reporting_map), row.names = FALSE)
  }
  note_lines <- as.character(x$notes %||% character(0))
  note_lines <- note_lines[nzchar(note_lines)]
  if (length(note_lines) > 0L) {
    cat("\nNotes\n")
    cat(paste0(" - ", note_lines, collapse = "\n"), "\n", sep = "")
  }
  invisible(x)
}

#' Plot a data-description object
#'
#' @param x Output from [describe_mfrm_data()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"score_distribution"`, `"facet_levels"`, or `"missing"`.
#' @param main Optional title override.
#' @param palette Optional named colors (`score`, `facet`, `missing`).
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method draws quick pre-fit quality views from [describe_mfrm_data()]:
#' - score distribution balance
#' - facet-level structure size
#' - missingness by selected columns
#'
#' @section Interpreting output:
#' - `"score_distribution"`: bar chart of weighted observation counts per
#'   score category.  Y-axis is `WeightedN` (sum of weights for each
#'   category).  Categories with very few observations (< 10) may produce
#'   unstable threshold estimates.  A roughly uniform or unimodal
#'   distribution is ideal; heavy floor/ceiling effects compress the
#'   measurement range.
#' - `"facet_levels"`: bar chart showing the number of distinct levels
#'   per facet.  Useful for verifying that the design structure matches
#'   expectations (e.g., expected number of raters or criteria).  Very
#'   large numbers of levels increase computation time and may require
#'   higher `maxit` in [fit_mfrm()].
#' - `"missing"`: bar chart of missing-value counts per input column.
#'   Columns with non-zero counts should be investigated before
#'   fitting---rows with missing scores, persons, or facet IDs are
#'   dropped during estimation.
#'
#' @section Typical workflow:
#' 1. Run [describe_mfrm_data()] before fitting.
#' 2. Inspect `summary(ds)` and `plot(ds, type = "missing")`.
#' 3. Check category/facet balance with other plot types.
#' 4. Fit model after resolving obvious data issues.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [describe_mfrm_data()], `plot()`
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(ds, draw = FALSE)
#' @export
plot.mfrm_data_description <- function(x,
                                       y = NULL,
                                       type = c("score_distribution", "facet_levels", "missing"),
                                       main = NULL,
                                       palette = NULL,
                                       label_angle = 45,
                                       draw = TRUE,
                                       ...) {
  type <- match.arg(tolower(as.character(type[1])), c("score_distribution", "facet_levels", "missing"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      score = "#2b8cbe",
      facet = "#31a354",
      missing = "#756bb1"
    )
  )

  if (type == "score_distribution") {
    tbl <- as.data.frame(x$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Score", "WeightedN") %in% names(tbl))) {
      stop("Score distribution is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$WeightedN)),
        labels = as.character(tbl$Score),
        col = pal["score"],
        main = if (is.null(main)) "Score distribution" else as.character(main[1]),
        ylab = "Weighted N",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_description",
      list(plot = "score_distribution", table = tbl)
    )))
  }

  if (type == "facet_levels") {
    tbl <- as.data.frame(x$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "Level") %in% names(tbl))) {
      stop("Facet level summary is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
    }
    agg <- tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(Levels = dplyr::n_distinct(.data$Level), .groups = "drop") |>
      dplyr::arrange(.data$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(agg$Levels)),
        labels = as.character(agg$Facet),
        col = pal["facet"],
        main = if (is.null(main)) "Facet levels" else as.character(main[1]),
        ylab = "Levels",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_description",
      list(plot = "facet_levels", table = as.data.frame(agg, stringsAsFactors = FALSE))
    )))
  }

  tbl <- as.data.frame(x$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Column", "Missing") %in% names(tbl))) {
    stop("Missing-by-column table is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
  }
  if (isTRUE(draw)) {
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl$Missing)),
      labels = as.character(tbl$Column),
      col = pal["missing"],
      main = if (is.null(main)) "Missing values by column" else as.character(main[1]),
      ylab = "Missing",
      label_angle = label_angle,
      mar_bottom = 8.0
    )
  }
  invisible(new_mfrm_plot_data(
    "data_description",
    list(plot = "missing", table = tbl)
  ))
}

#' Audit and normalize anchor/group-anchor tables
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param anchors Optional anchor table (Facet, Level, Anchor).
#' @param group_anchors Optional group-anchor table
#'   (Facet, Level, Group, GroupValue).
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param keep_original Keep original category values.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   recommendations (default `5`).
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in recommendations (default `30`).
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in recommendations (default `10`).
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#'
#' @details
#' **Anchoring** (also called "fixing" or scale linking) constrains selected
#' parameter estimates to pre-specified values, placing the current
#' analysis on a previously established scale.  This is essential when
#' comparing results across administrations, linking test forms, or
#' monitoring rater drift over time.
#'
#' This function applies the same preprocessing and key-resolution rules
#' as `fit_mfrm()`, but returns an audit object so constraints can be
#' checked *before* estimation.  Running the audit first helps avoid
#' estimation failures caused by misspecified or data-incompatible
#' anchors.
#'
#' **Anchor types:**
#' - *Direct anchors* fix individual element measures to specific logit
#'   values (e.g., Rater R1 anchored at 0.35 logits).
#' - *Group anchors* constrain the mean of a set of elements to a
#'   target value, allowing individual elements to vary freely around
#'   that mean.
#' - When both types overlap for the same element, the direct anchor
#'   takes precedence.
#'
#' **Design checks** verify that each anchored element has at least
#' `min_obs_per_element` weighted observations (default 30) and each
#' score category has at least `min_obs_per_category` (default 10).
#' These thresholds follow standard Rasch sample-size recommendations
#' (Linacre, 1994).
#'
#' @section Interpreting output:
#' - `issue_counts`/`issues`: concrete data or specification problems.
#' - `facet_summary`: constraint coverage by facet.
#' - `design_checks`: whether anchor targets have enough observations.
#' - `recommendations`: action items before estimation.
#'
#' @section Typical workflow:
#' 1. Build candidate anchors (e.g., with [make_anchor_table()]).
#' 2. Run `audit_mfrm_anchors(...)`.
#' 3. Resolve issues, then fit with [fit_mfrm()].
#'
#' @return A list of class `mfrm_anchor_audit` with:
#' - `anchors`: cleaned anchor table used by estimation
#' - `group_anchors`: cleaned group-anchor table used by estimation
#' - `facet_summary`: counts of levels, constrained levels, and free levels
#' - `design_checks`: observation-count checks by level/category
#' - `thresholds`: active threshold settings used for recommendations
#' - `issue_counts`: issue-type counts
#' - `issues`: list of issue tables
#' - `recommendations`: package-native anchor guidance strings
#'
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [make_anchor_table()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#'
#' anchors <- data.frame(
#'   Facet = c("Rater", "Rater"),
#'   Level = c("R1", "R1"),
#'   Anchor = c(0, 0.1),
#'   stringsAsFactors = FALSE
#' )
#' aud <- audit_mfrm_anchors(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   anchors = anchors
#' )
#' aud$issue_counts
#' summary(aud)
#' p_aud <- plot(aud, draw = FALSE)
#' p_aud$data$plot
#' @export
audit_mfrm_anchors <- function(data,
                               person,
                               facets,
                               score,
                               anchors = NULL,
                               group_anchors = NULL,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               min_common_anchors = 5L,
                               min_obs_per_element = 30,
                               min_obs_per_category = 10,
                               noncenter_facet = "Person",
                               dummy_facets = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original
  )

  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)

  audit <- audit_anchor_tables(
    prep = prep,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )
  class(audit) <- c("mfrm_anchor_audit", class(audit))
  audit
}

#' @export
print.mfrm_anchor_audit <- function(x, ...) {
  issue_total <- if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) sum(x$issue_counts$N) else 0
  cat("mfrm anchor audit\n")
  cat("  issue rows: ", issue_total, "\n", sep = "")

  if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
    nonzero <- x$issue_counts |>
      dplyr::filter(.data$N > 0)
    if (nrow(nonzero) > 0) {
      cat("\nIssue counts\n")
      print(nonzero, row.names = FALSE)
    }
  }

  if (!is.null(x$facet_summary) && nrow(x$facet_summary) > 0) {
    cat("\nFacet summary\n")
    print(x$facet_summary, row.names = FALSE)
  }

  if (!is.null(x$design_checks) &&
      !is.null(x$design_checks$level_observation_summary) &&
      nrow(x$design_checks$level_observation_summary) > 0) {
    cat("\nLevel observation summary\n")
    print(x$design_checks$level_observation_summary, row.names = FALSE)
  }

  invisible(x)
}

#' Summarize an anchor-audit object
#'
#' @param object Output from [audit_mfrm_anchors()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in issue previews.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary provides a compact pre-estimation audit of anchor and
#' group-anchor specifications.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `issue_counts`: primary triage table (non-zero issues first).
#' - `facet_summary`: anchored/grouped/free-level balance by facet.
#' - `level_observation_summary` and `category_counts`: sparse-cell diagnostics.
#' - `recommendations`: concrete remediation suggestions.
#'
#' If `issue_counts` is non-empty, treat anchor constraints as provisional and
#' resolve issues before final estimation.
#'
#' @section Typical workflow:
#' 1. Run [audit_mfrm_anchors()] with intended anchors/group anchors.
#' 2. Review `summary(aud)` and recommendations.
#' 3. Revise anchor tables, then call [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_anchor_audit`.
#' @seealso [audit_mfrm_anchors()], [fit_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' aud <- audit_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(aud)
#' @export
summary.mfrm_anchor_audit <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  issue_counts <- as.data.frame(object$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(issue_counts) > 0 && all(c("Issue", "N") %in% names(issue_counts))) {
    issue_counts <- issue_counts |>
      dplyr::filter(.data$N > 0) |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue) |>
      dplyr::slice_head(n = top_n)
  }

  facet_summary <- as.data.frame(object$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  level_summary <- as.data.frame(object$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  category_summary <- as.data.frame(object$design_checks$category_counts %||% data.frame(), stringsAsFactors = FALSE)

  recommendations <- as.character(object$recommendations %||% character(0))
  if (length(recommendations) > top_n) {
    recommendations <- recommendations[seq_len(top_n)]
  }

  notes <- if (nrow(issue_counts) > 0) {
    "Anchor-audit issues were detected. Review issue counts and recommendations."
  } else {
    "No anchor-table issue rows were detected."
  }

  out <- list(
    issue_counts = issue_counts,
    facet_summary = facet_summary,
    level_observation_summary = level_summary,
    category_counts = category_summary,
    recommendations = recommendations,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_anchor_audit"
  out
}

#' @export
print.summary.mfrm_anchor_audit <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Anchor Audit Summary\n")
  if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
    cat("\nIssue counts\n")
    print(round_numeric_df(as.data.frame(x$issue_counts), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_summary) && nrow(x$facet_summary) > 0) {
    cat("\nFacet summary\n")
    print(round_numeric_df(as.data.frame(x$facet_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$level_observation_summary) && nrow(x$level_observation_summary) > 0) {
    cat("\nLevel observation summary\n")
    print(round_numeric_df(as.data.frame(x$level_observation_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$category_counts) && nrow(x$category_counts) > 0) {
    cat("\nCategory counts\n")
    print(round_numeric_df(as.data.frame(x$category_counts), digits = digits), row.names = FALSE)
  }
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations\n")
    for (line in x$recommendations) cat(" - ", line, "\n", sep = "")
  }
  if (!is.null(x$notes) && nzchar(x$notes)) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Plot an anchor-audit object
#'
#' @param x Output from [audit_mfrm_anchors()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"issue_counts"`, `"facet_constraints"`,
#'   or `"level_observations"`.
#' @param main Optional title override.
#' @param palette Optional named colors.
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Base-R visualization helper for anchor audit outputs.
#'
#' @section Interpreting output:
#' - `"issue_counts"`: volume of each issue class.
#' - `"facet_constraints"`: anchored/grouped/free mix by facet.
#' - `"level_observations"`: observation support across levels.
#'
#' @section Typical workflow:
#' 1. Run [audit_mfrm_anchors()].
#' 2. Start with `plot(aud, type = "issue_counts")`.
#' 3. Inspect constraint and support plots before fitting.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [audit_mfrm_anchors()], [make_anchor_table()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' aud <- audit_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(aud, draw = FALSE)
#' @export
plot.mfrm_anchor_audit <- function(x,
                                   y = NULL,
                                   type = c("issue_counts", "facet_constraints", "level_observations"),
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45,
                                   draw = TRUE,
                                   ...) {
  type <- match.arg(tolower(as.character(type[1])), c("issue_counts", "facet_constraints", "level_observations"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      issues = "#cb181d",
      anchored = "#756bb1",
      grouped = "#9ecae1",
      levels = "#2b8cbe"
    )
  )

  if (type == "issue_counts") {
    tbl <- as.data.frame(x$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Issue", "N") %in% names(tbl))) {
      stop("Issue-count table is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    tbl <- tbl |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$N)),
        labels = as.character(tbl$Issue),
        col = pal["issues"],
        main = if (is.null(main)) "Anchor-audit issue counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_audit",
      list(plot = "issue_counts", table = tbl)
    )))
  }

  if (type == "facet_constraints") {
    tbl <- as.data.frame(x$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "AnchoredLevels", "GroupedLevels", "FreeLevels") %in% names(tbl))) {
      stop("Facet summary with constraint columns is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar), add = TRUE)
      mar <- old_mar
      mar[1] <- max(mar[1], 8.8)
      graphics::par(mar = mar)
      mat <- rbind(
        Anchored = suppressWarnings(as.numeric(tbl$AnchoredLevels)),
        Grouped = suppressWarnings(as.numeric(tbl$GroupedLevels)),
        Free = suppressWarnings(as.numeric(tbl$FreeLevels))
      )
      mids <- graphics::barplot(
        height = mat,
        beside = FALSE,
        names.arg = FALSE,
        col = c(pal["anchored"], pal["grouped"], "#d9d9d9"),
        border = "white",
        ylab = "Levels",
        main = if (is.null(main)) "Constraint profile by facet" else as.character(main[1])
      )
      draw_rotated_x_labels(
        at = mids,
        labels = as.character(tbl$Facet),
        srt = label_angle,
        cex = 0.82,
        line_offset = 0.085
      )
      graphics::legend(
        "topright",
        legend = c("Anchored", "Grouped", "Free"),
        fill = c(pal["anchored"], pal["grouped"], "#d9d9d9"),
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_audit",
      list(plot = "facet_constraints", table = tbl)
    )))
  }

  tbl <- as.data.frame(x$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Facet", "MinObsPerLevel") %in% names(tbl))) {
    stop("Level observation summary is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
  }
  if (isTRUE(draw)) {
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl$MinObsPerLevel)),
      labels = as.character(tbl$Facet),
      col = pal["levels"],
      main = if (is.null(main)) "Minimum observations per level" else as.character(main[1]),
      ylab = "Min observations",
      label_angle = label_angle,
      mar_bottom = 8.0
    )
    if ("RecommendedMinObs" %in% names(tbl)) {
      r <- suppressWarnings(as.numeric(tbl$RecommendedMinObs))
      r <- r[is.finite(r)]
      if (length(r) > 0) graphics::abline(h = unique(r)[1], lty = 2, col = "gray45")
    }
  }
  invisible(new_mfrm_plot_data(
    "anchor_audit",
    list(plot = "level_observations", table = tbl)
  ))
}

#' Build an anchor table from fitted estimates
#'
#' @param fit Output from [fit_mfrm()].
#' @param facets Optional subset of facets to include.
#' @param include_person Include person estimates as anchors.
#' @param digits Rounding digits for anchor values.
#'
#' @details
#' This function exports estimated facet parameters as an anchor table
#' for use in subsequent calibrations.  This is the standard approach
#' for **linking** across administrations: a reference
#' run establishes the measurement scale, and anchored re-analyses
#' place new data on that same scale.
#'
#' Anchor values should be exported from a well-fitting reference run
#' with adequate sample size.  If the reference model has convergence
#' issues or large misfit, the exported anchors may propagate
#' instability.  Re-run [audit_mfrm_anchors()] on the receiving data
#' to verify compatibility before estimation.
#'
#' The `digits` parameter controls rounding precision.  Use at least 4
#' digits for research applications; excessive rounding (e.g., 1 digit)
#' can introduce avoidable calibration error.
#'
#' @section Interpreting output:
#' - `Facet`: facet name to be anchored in later runs.
#' - `Level`: specific element/level name inside that facet.
#' - `Anchor`: fixed logit value (rounded by `digits`).
#'
#' @section Typical workflow:
#' 1. Fit a reference run with [fit_mfrm()].
#' 2. Export anchors with `make_anchor_table(fit)`.
#' 3. Pass selected rows back into `fit_mfrm(..., anchors = ...)`.
#'
#' @return A data.frame with `Facet`, `Level`, and `Anchor`.
#' @seealso [fit_mfrm()], [audit_mfrm_anchors()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' anchors_tbl <- make_anchor_table(fit)
#' head(anchors_tbl)
#' summary(anchors_tbl$Anchor)
#' @export
make_anchor_table <- function(fit,
                              facets = NULL,
                              include_person = FALSE,
                              digits = 6) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }

  digits <- max(0L, as.integer(digits))
  out <- tibble::tibble()

  if (isTRUE(include_person) && !is.null(fit$facets$person) && nrow(fit$facets$person) > 0) {
    per <- tibble::as_tibble(fit$facets$person)
    if ("Person" %in% names(per) && "Estimate" %in% names(per)) {
      out <- dplyr::bind_rows(
        out,
        per |>
          dplyr::transmute(
            Facet = "Person",
            Level = as.character(.data$Person),
            Anchor = round(as.numeric(.data$Estimate), digits = digits)
          )
      )
    }
  }

  others <- tibble::as_tibble(fit$facets$others)
  if (nrow(others) > 0 && all(c("Facet", "Level", "Estimate") %in% names(others))) {
    out <- dplyr::bind_rows(
      out,
      others |>
        dplyr::transmute(
          Facet = as.character(.data$Facet),
          Level = as.character(.data$Level),
          Anchor = round(as.numeric(.data$Estimate), digits = digits)
        )
    )
  }

  if (!is.null(facets)) {
    keep <- as.character(facets)
    out <- out |>
      dplyr::filter(.data$Facet %in% keep)
  }

  out |>
    dplyr::arrange(.data$Facet, .data$Level)
}

#' Compute diagnostics for an `mfrm_fit` object
#'
#' @param fit Output from [fit_mfrm()].
#' @param interaction_pairs Optional list of facet pairs.
#' @param top_n_interactions Number of top interactions.
#' @param whexact Use exact ZSTD transformation.
#' @param diagnostic_mode Diagnostic basis to compute: `"legacy"` keeps the
#'   residual/EAP-based stack only, `"marginal_fit"` adds the strict
#'   latent-integrated first-order marginal-fit companion, and `"both"`
#'   computes both paths.
#' @param residual_pca Residual PCA mode: `"none"`, `"overall"`, `"facet"`, or `"both"`.
#' @param pca_max_factors Maximum number of PCA factors to retain per matrix.
#'
#' @details
#' This function computes a diagnostic bundle used by downstream reporting.
#' It calculates element-level fit statistics, approximate facet
#' separation/reliability summaries, residual-based QC diagnostics, and
#' optionally residual PCA for
#' exploratory residual-structure screening.
#'
#' `diagnostic_mode` keeps the legacy residual fit path explicit rather than
#' silently replacing it. The legacy path is a compatibility-oriented
#' residual/EAP stack, whereas the strict marginal path targets
#' latent-integrated first-order category counts. When `diagnostic_mode =
#' "both"`, the output includes a `diagnostic_basis` guide so downstream
#' tables and summaries can distinguish these targets.
#'
#' Choosing `diagnostic_mode`:
#' - `"legacy"`: use when continuity with historical residual-based workflows is
#'   the priority.
#' - `"marginal_fit"`: use when you want the strict latent-integrated screen
#'   without the extra legacy bundle.
#' - `"both"`: recommended when you want continuity with the legacy residual
#'   stack while making the strict marginal path explicit for `RSM`, `PCM`,
#'   and bounded `GPCM` fits.
#'
#' For bounded `GPCM`, the same generalized partial credit kernel now
#' drives both the residual/probability tables and the strict marginal
#' category-fit companion. Residual-based MnSq summaries should still be read
#' as exploratory screening tools rather than strict Rasch-style invariance
#' tests because discrimination is free, and the strict marginal companion
#' should likewise be treated as a slope-aware screen rather than a finalized
#' inferential test family.
#'
#' **Key fit statistics computed for each element:**
#' - **Infit MnSq**: information-weighted mean-square residual; sensitive
#'   to on-target misfitting patterns.  Expected value = 1.0.
#' - **Outfit MnSq**: unweighted mean-square residual; sensitive to
#'   off-target outliers.  Expected value = 1.0.
#' - **ZSTD**: Wilson-Hilferty cube-root transformation of MnSq to an
#'   approximate standard normal deviate.
#' - **PTMEA**: point-measure correlation (item-rest correlation in MFRM
#'   context); positive values confirm alignment with the latent trait.
#'
#' **Misfit flagging guidelines (Bond & Fox, 2015):**
#' - MnSq < 0.5: overfit (too predictable; may inflate reliability)
#' - MnSq 0.5--1.5: productive for measurement
#' - MnSq > 1.5: underfit (noise degrades measurement)
#' - \eqn{|\mathrm{ZSTD}| > 2}: statistically significant misfit (5\%)
#'
#' When Infit and Outfit disagree, Infit is generally more informative
#' because it downweights extreme observations.  Large Outfit with
#' acceptable Infit typically indicates a few outlying responses rather
#' than systematic misfit.
#'
#' `interaction_pairs` controls which facet interactions are summarized.
#' Each element can be:
#' - a length-2 character vector such as `c("Rater", "Criterion")`, or
#' - omitted (`NULL`) to let the function select top interactions automatically.
#'
#' Residual PCA behavior:
#' - `"none"`: skip PCA (fastest; recommended for initial exploration)
#' - `"overall"`: compute overall residual PCA across all facets
#' - `"facet"`: compute facet-specific residual PCA for each facet
#' - `"both"`: compute both overall and facet-specific PCA
#'
#' Overall PCA examines the person \eqn{\times} combined-facet residual
#' matrix; facet-specific PCA examines person \eqn{\times} facet-level
#' matrices. These summaries are exploratory screens for residual
#' structure, not standalone proofs for or against unidimensionality.
#' Facet-specific PCA can help localise where a stronger residual signal
#' is concentrated.
#'
#' @section Reading key components:
#' Practical interpretation often starts with:
#' - `overall_fit`: global infit/outfit and degrees of freedom.
#' - `reliability`: facet-level model/real separation and reliability. `MML`
#'   uses model-based `ModelSE` values where available; `JML` keeps these
#'   quantities as exploratory approximations.
#' - `fit`: element-level misfit scan (`Infit`, `Outfit`, `ZSTD`).
#' - `unexpected`, `fair_average`, `displacement`: targeted QC bundles.
#'   For bounded `GPCM`, `fair_average` is retained as an unavailable
#'   placeholder because that compatibility calculation has not yet been
#'   validated for the generalized model.
#' - `approximation_notes`: method notes for SE/CI/reliability summaries.
#'
#' @section Interpreting output:
#' Start with `overall_fit` and `reliability`, then move to element-level
#' diagnostics (`fit`) and targeted bundles (`unexpected`, `displacement`,
#' `interrater`, `facets_chisq`). Treat `fair_average` as available only for
#' the `RSM` / `PCM` branch.
#'
#' Consistent signals across multiple components are typically more robust than
#' a single isolated warning.  For example, an element flagged for both high
#' Outfit and high displacement is more concerning than one flagged on a
#' single criterion.
#'
#' `SE` is kept as a compatibility alias for `ModelSE`. `RealSE` is a
#' fit-adjusted companion defined as `ModelSE * sqrt(max(Infit, 1))`.
#' Reliability tables report model and fit-adjusted bounds from observed
#' variance, error variance, and true variance; `JML` entries should still be
#' treated as exploratory.
#'
#' @section Typical workflow:
#' 1. Start with `diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")`.
#' 2. Inspect `summary(diag)` and use `diagnostic_basis` to separate legacy residual evidence from strict marginal evidence.
#' 3. If needed, rerun with residual PCA (`"overall"` or `"both"`).
#'
#' @return
#' An object of class `mfrm_diagnostics` including:
#' - `obs`: observed/expected/residual-level table
#' - `measures`: facet/person fit table (`Infit`, `Outfit`, `ZSTD`, `PTMEA`)
#' - `overall_fit`: overall fit summary
#' - `fit`: element-level fit diagnostics
#' - `reliability`: facet-level model/real separation and reliability
#' - `precision_profile`: one-row summary of the active precision tier and its
#'   recommended use
#' - `precision_audit`: package-native checks for SE, CI, and reliability
#' - `facet_precision`: facet-level precision summary by distribution basis and
#'   SE mode
#' - `facets_chisq`: fixed/random facet variability summary
#' - `interactions`: top interaction diagnostics
#' - `interrater`: inter-rater agreement bundle (`summary`, `pairs`) including
#'   agreement and rater-severity spread indices
#' - `unexpected`: unexpected-response bundle
#' - `fair_average`: adjusted-score reference bundle (placeholder only for
#'   bounded `GPCM`)
#' - `displacement`: displacement diagnostics bundle
#' - `approximation_notes`: method notes for SE/CI/reliability summaries
#' - `diagnostic_basis`: guide to the statistical target of each diagnostic path
#' - `marginal_fit`: optional strict marginal-fit companion based on
#'   posterior-expected first-order category counts
#' - `residual_pca_overall`: optional overall PCA object
#' - `residual_pca_by_facet`: optional facet PCA objects
#'
#' @seealso [fit_mfrm()], [analyze_residual_pca()], [build_visual_summaries()],
#'   [mfrmr_visual_diagnostics], [mfrmr_reporting_and_apa]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
#' s_diag <- summary(diag)
#' s_diag$overview[, c("Observations", "Facets", "Categories")]
#' s_diag$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]
#' p_qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE)
#' p_qc$data$plot
#'
#' # Optional: include residual PCA in the diagnostic bundle
#' diag_pca <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag_pca, mode = "overall")
#' head(pca$overall_table)
#'
#' # Reporting route:
#' prec <- precision_audit_report(fit, diagnostics = diag)
#' summary(prec)
#' }
#' @export
diagnose_mfrm <- function(fit,
                          interaction_pairs = NULL,
                          top_n_interactions = 20,
                          whexact = FALSE,
                          diagnostic_mode = c("legacy", "marginal_fit", "both"),
                          residual_pca = c("none", "overall", "facet", "both"),
                          pca_max_factors = 10L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  diagnostic_mode <- match.arg(diagnostic_mode)
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))

  out <- mfrm_diagnostics(
    fit,
    interaction_pairs = interaction_pairs,
    top_n_interactions = top_n_interactions,
    whexact = whexact,
    diagnostic_mode = diagnostic_mode,
    residual_pca = residual_pca,
    pca_max_factors = pca_max_factors
  )
  class(out) <- c("mfrm_diagnostics", class(out))
  out
}

#' Compare two or more fitted MFRM models
#'
#' Produce a side-by-side comparison of multiple [fit_mfrm()] results using
#' information criteria, log-likelihood, and parameter counts. When exactly
#' two models are supplied and the current conservative nesting audit passes,
#' a likelihood-ratio test is included.
#'
#' @param ... Two or more `mfrm_fit` objects to compare.
#' @param labels Optional character vector of labels for each model.
#'   If `NULL`, labels are generated from model/method combinations.
#' @param warn_constraints Logical. If `TRUE` (the default), emit a warning
#'   when models use different centering constraints (`noncenter_facet` or
#'   `dummy_facets`), which can make information-criterion comparisons
#'   misleading.
#' @param nested Logical. Set to `TRUE` only when the supplied models are
#'   known to be nested and fitted with the same likelihood basis on the same
#'   observations. The default is `FALSE`, in which case no likelihood-ratio
#'   test is reported. When `TRUE`, the function still runs a conservative
#'   structural audit and computes the LRT only for supported nesting patterns.
#'
#' @details
#' Models should be fit to the **same data** (same rows, same person/facet
#' columns) for the comparison to be meaningful. The function checks that
#' observation counts match and warns otherwise.
#'
#' Information-criterion ranking is reported only when all candidate models
#' use the package's `MML` estimation path, analyze the same observations, and
#' converge successfully. Raw `AIC` and `BIC` values are still shown for each
#' model, but `Delta_*`, weights, and preferred-model summaries are suppressed
#' when the likelihood basis is not comparable enough for primary reporting.
#'
#' **Nesting**: Two models are *nested* when one is a special case of the
#' other obtained by imposing equality constraints.  The most common
#' nesting in MFRM is RSM (shared thresholds) inside PCM
#' (item-specific thresholds).  Models that differ only in estimation
#' method (MML vs JML) on the same specification are not nested in the
#' usual sense---use information criteria rather than LRT for that
#' comparison.
#'
#' In the **current `mfrmr` model space**, the automatic nesting audit is
#' intentionally conservative: it treats `RSM` nested inside `PCM` under shared
#' data and shared constraints as the only supported automatic relation.
#' Same-family comparisons, cross-method comparisons, or comparisons that
#' change anchors/dummying/centering are not automatically promoted to LRT
#' claims.
#'
#' The **likelihood-ratio test (LRT)** is reported only when exactly two
#' models are supplied, `nested = TRUE`, the structural audit passes, and the
#' difference in the number of parameters is positive:
#'
#' \deqn{\Lambda = -2 (\ell_{\mathrm{restricted}} - \ell_{\mathrm{full}})
#'   \sim \chi^2_{\Delta p}}
#'
#' The LRT is asymptotically valid when models are nested and the data
#' are independent.  With small samples or boundary conditions (e.g.,
#' variance components near zero), treat p-values as approximate.
#'
#' @section Information-criterion diagnostics:
#' In addition to raw AIC and BIC values, the function computes:
#' - **Delta_AIC / Delta_BIC**: difference from the best (minimum) value.
#'   A Delta < 2 is typically considered negligible; 4--7 suggests
#'   moderate evidence; > 10 indicates strong evidence against the
#'   higher-scoring model (Burnham & Anderson, 2002).
#' - **AkaikeWeight / BICWeight**: model probabilities derived from
#'   `exp(-0.5 * Delta)`, normalised across the candidate set.  An
#'   Akaike weight of 0.90 means the model has a 90\% probability of
#'   being the best in the candidate set.
#' - **Evidence ratios**: pairwise ratios of Akaike weights, quantifying
#'   the relative evidence for one model over another (e.g., an
#'   evidence ratio of 5 means the preferred model is 5 times more
#'   likely).
#'
#' AIC penalises complexity less than BIC; when they disagree, AIC
#' favours the more complex model and BIC the simpler one.
#'
#' @section What this comparison means:
#' `compare_mfrm()` is a same-basis model-comparison helper. Its strongest
#' claims apply only when the models were fit to the same response data,
#' under a compatible likelihood basis, and with compatible constraint
#' structure.
#'
#' @section What this comparison does not justify:
#' - Do not treat AIC/BIC differences as primary evidence when
#'   `table$ICComparable` is `FALSE`.
#' - Do not interpret the LRT unless `nested = TRUE` and the structural audit
#'   in `comparison_basis$nesting_audit` passes.
#' - Do not assume that `nested = TRUE` overrides the package's conservative
#'   nesting boundary; unsupported relations remain unsupported.
#' - Do not compare models fit to different datasets, different score codings,
#'   or materially different constraint systems as if they were commensurate.
#'
#' @section Interpreting output:
#' - Lower AIC/BIC values indicate better parsimony-accuracy trade-off only
#'   when `table$ICComparable` is `TRUE`.
#' - A significant LRT p-value suggests the more complex model provides a
#'   meaningfully better fit only when the nesting assumption truly holds.
#' - `preferred` indicates the model preferred by each criterion.
#' - `evidence_ratios` gives pairwise Akaike-weight ratios (returned only
#'   when Akaike weights can be computed for at least two models).
#' - When comparing more than two models, interpret evidence ratios
#'   cautiously---they do not adjust for multiple comparisons.
#'
#' @section How to read the main outputs:
#' - `table`: first-pass comparison table; start with `ICComparable`,
#'   `Model`, `Method`, `AIC`, and `BIC`.
#' - `comparison_basis`: records whether IC and LRT claims are defensible for
#'   the supplied models. Inspect `comparison_basis$nesting_audit$relation`
#'   and `reason` before reading any LRT output.
#' - `lrt`: nested-model test summary, present only when the requested and
#'   audited conditions are met.
#' - `preferred`: candidate preferred by each criterion when those summaries
#'   are available.
#'
#' @section Recommended next step:
#' Inspect `comparison_basis` before writing conclusions. If comparability is
#' weak, treat the result as descriptive and revise the model setup (for
#' example, explicit `step_facet`, common data, or common constraints) before
#' using IC or LRT results in reporting.
#'
#' @section Typical workflow:
#' 1. Fit two models with [fit_mfrm()] (e.g., RSM and PCM).
#' 2. Compare with `compare_mfrm(fit_rsm, fit_pcm)`.
#' 3. Inspect `summary(comparison)` for AIC/BIC diagnostics and, when
#'    appropriate, an LRT.
#'
#' @return
#' An object of class `mfrm_comparison` (named list) with:
#' - `table`: data.frame of model-level statistics (LogLik, AIC, BIC,
#'   Delta_AIC, AkaikeWeight, Delta_BIC, BICWeight, npar, nobs, Model,
#'   Method, Converged, ICComparable).
#' - `lrt`: data.frame with likelihood-ratio test result (only when two models
#'   are supplied and `nested = TRUE`). Contains `ChiSq`, `df`, `p_value`.
#' - `evidence_ratios`: data.frame of pairwise Akaike-weight ratios (Model1,
#'   Model2, EvidenceRatio). `NULL` when weights cannot be computed.
#' - `preferred`: named list with the preferred model label by each criterion.
#' - `comparison_basis`: list describing whether IC and LRT comparisons were
#'   considered comparable. Includes a conservative `nesting_audit`.
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#'
#' fit_rsm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "MML", model = "RSM", maxit = 25)
#' fit_pcm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "MML", model = "PCM",
#'                      step_facet = "Criterion", maxit = 25)
#' comp <- compare_mfrm(fit_rsm, fit_pcm, labels = c("RSM", "PCM"))
#' comp$table
#' comp$evidence_ratios
#' }
#' @export
compare_mfrm <- function(..., labels = NULL, warn_constraints = TRUE, nested = FALSE) {
  fits <- list(...)
  if (length(fits) < 2) {
    stop("`compare_mfrm()` requires at least two `mfrm_fit` objects.",
         call. = FALSE)
  }
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "mfrm_fit")) {
      stop("Argument ", i, " is not an `mfrm_fit` object. Got: ",
           class(fits[[i]])[1], ".", call. = FALSE)
    }
  }

  # -- build labels --
  if (is.null(labels)) {
    labels <- vapply(fits, function(f) {
      m <- if (!is.null(f$config$model)) toupper(f$config$model[1]) else "?"
      e <- if (!is.null(f$config$method)) toupper(f$config$method[1]) else "?"
      paste0(m, "/", e)
    }, character(1))
    if (anyDuplicated(labels)) {
      labels <- paste0(labels, " (", seq_along(labels), ")")
    }
  }
  if (length(labels) != length(fits)) {
    stop("`labels` must have the same length as the number of models (",
         length(fits), ").", call. = FALSE)
  }

  normalize_compare_col <- function(col) {
    if (is.factor(col)) return(as.character(col))
    if (inherits(col, "POSIXt")) return(format(col, tz = "UTC", usetz = TRUE))
    if (is.numeric(col)) return(format(col, digits = 17, trim = TRUE, scientific = FALSE))
    if (is.logical(col)) return(ifelse(is.na(col), "NA", ifelse(col, "TRUE", "FALSE")))
    if (is.character(col)) return(col)
    as.character(col)
  }

  canonicalize_compare_data <- function(dat) {
    if (is.null(dat) || !is.data.frame(dat)) return(NULL)
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
    dat <- dat[, sort(names(dat)), drop = FALSE]
    dat[] <- lapply(dat, normalize_compare_col)
    if (nrow(dat) > 0 && ncol(dat) > 0) {
      ord <- do.call(order, c(dat, list(na.last = TRUE, method = "radix")))
      dat <- dat[ord, , drop = FALSE]
    }
    rownames(dat) <- NULL
    dat
  }

  compare_data_equal <- function(fit_a, fit_b) {
    dat_a <- canonicalize_compare_data(fit_a$prep$data %||% NULL)
    dat_b <- canonicalize_compare_data(fit_b$prep$data %||% NULL)
    if (is.null(dat_a) || is.null(dat_b)) return(FALSE)
    identical(dat_a, dat_b)
  }

  # -- constraint compatibility check --
  if (isTRUE(warn_constraints) && length(fits) >= 2) {
    noncenter <- vapply(fits, function(f) {
      nc <- f$config$noncenter_facet
      if (is.null(nc)) "Person" else as.character(nc[1])
    }, character(1))
    if (length(unique(noncenter)) > 1) {
      warning("Models use different centering constraints (",
              paste(unique(noncenter), collapse = ", "),
              "). IC comparisons may be misleading.", call. = FALSE)
    }
    dummy_sets <- lapply(fits, function(f) {
      d <- f$config$dummy_facets
      if (is.null(d)) character(0) else sort(as.character(d))
    })
    dummy_sigs <- vapply(dummy_sets, paste, character(1), collapse = ",")
    if (length(unique(dummy_sigs)) > 1) {
      warning("Models use different dummy-facet constraints (",
              paste(unique(dummy_sigs), collapse = " vs "),
              "). IC comparisons may be misleading.", call. = FALSE)
    }
  }

  # -- extract summary statistics --
  rows <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    s <- f$summary
    nobs <- if (!is.null(f$prep$data)) nrow(f$prep$data) else NA_integer_
    npar <- if (!is.null(f$opt$par)) length(f$opt$par) else NA_integer_
    method <- if (!is.null(f$config$method)) toupper(f$config$method[1]) else NA_character_
    method <- ifelse(identical(method, "JMLE"), "JML", method)
    tibble(
      Label     = labels[i],
      Model     = if (!is.null(f$config$model)) toupper(f$config$model[1]) else NA_character_,
      Method    = method,
      nobs      = nobs,
      npar      = npar,
      LogLik    = if ("LogLik" %in% names(s)) s$LogLik[1] else NA_real_,
      AIC       = if ("AIC" %in% names(s)) s$AIC[1] else NA_real_,
      BIC       = if ("BIC" %in% names(s)) s$BIC[1] else NA_real_,
      Converged = if ("Converged" %in% names(s)) s$Converged[1] else NA
    )
  })
  tbl <- bind_rows(rows)

  method_vals <- tbl$Method[!is.na(tbl$Method)]
  same_method <- length(unique(method_vals)) <= 1
  all_mml <- length(method_vals) > 0 && all(method_vals == "MML")
  obs_vals <- tbl$nobs[is.finite(tbl$nobs)]
  same_nobs <- length(unique(obs_vals)) <= 1
  same_data <- if (length(fits) >= 2) {
    all(vapply(fits[-1], function(f) compare_data_equal(fits[[1]], f), logical(1)))
  } else {
    TRUE
  }
  conv_vals <- tbl$Converged
  all_converged <- length(conv_vals) > 0 && all(!is.na(conv_vals) & as.logical(conv_vals))
  ic_comparable <- same_method && same_nobs && same_data && all_converged && all_mml

  if (!isTRUE(nested)) {
    if (!same_method) {
      warning(
        "Models use different estimation methods (",
        paste(unique(method_vals), collapse = ", "),
        "). Raw AIC/BIC values are shown, but cross-method deltas, weights, ",
        "and automatic preferences are suppressed.",
        call. = FALSE
      )
    } else if (!all_mml && all_converged) {
      warning(
        "Information-criterion ranking is limited to converged MML fits in this package. ",
        "Raw AIC/BIC values are shown for ",
        paste(unique(method_vals), collapse = ", "),
        " models, but deltas, weights, automatic preferences, and LRT were suppressed.",
        call. = FALSE
      )
    }

    # -- warn if observation counts differ --
    if (!same_nobs) {
      warning("Models were fit to different numbers of observations (",
              paste(obs_vals, collapse = ", "),
              "). Raw AIC/BIC values are shown, but cross-sample deltas, weights, ",
              "and automatic preferences are suppressed.", call. = FALSE)
    }
    if (!same_data) {
      warning(
        "Models were not fit to the same prepared response data. Raw AIC/BIC values are shown, ",
        "but deltas, weights, automatic preferences, and likelihood-ratio testing were suppressed.",
        call. = FALSE
      )
    }
    if (!all_converged) {
      warning(
        "At least one compared model did not converge. Raw AIC/BIC values are shown, ",
        "but IC ranking, weights, and likelihood-ratio testing were suppressed.",
        call. = FALSE
      )
    }
  }

  tbl$ICComparable <- ic_comparable
  preferred <- list()

  # -- Delta AIC and Akaike Weights --
  if (ic_comparable && any(is.finite(tbl$AIC))) {
    min_aic <- min(tbl$AIC, na.rm = TRUE)
    tbl$Delta_AIC <- tbl$AIC - min_aic
    raw_w <- exp(-0.5 * tbl$Delta_AIC)
    tbl$AkaikeWeight <- ifelse(is.finite(raw_w),
                               raw_w / sum(raw_w, na.rm = TRUE),
                               NA_real_)
    preferred$AIC <- tbl$Label[which.min(tbl$AIC)]
  } else {
    tbl$Delta_AIC <- NA_real_
    tbl$AkaikeWeight <- NA_real_
  }

  # -- Delta BIC and BIC Weights --
  if (ic_comparable && any(is.finite(tbl$BIC))) {
    min_bic <- min(tbl$BIC, na.rm = TRUE)
    tbl$Delta_BIC <- tbl$BIC - min_bic
    raw_bw <- exp(-0.5 * tbl$Delta_BIC)
    tbl$BICWeight <- ifelse(is.finite(raw_bw),
                            raw_bw / sum(raw_bw, na.rm = TRUE),
                            NA_real_)
    preferred$BIC <- tbl$Label[which.min(tbl$BIC)]
  } else {
    tbl$Delta_BIC <- NA_real_
    tbl$BICWeight <- NA_real_
  }

  # -- likelihood-ratio test (two models only) --
  lrt <- NULL
  nesting_audit <- audit_compare_mfrm_nesting(fits, labels = labels)
  if (isTRUE(nested) && length(fits) == 2) {
    if (!same_method || !same_nobs || !same_data || !all_converged || !all_mml) {
      warning(
        "`nested = TRUE` was requested, but the models do not share the same ",
        "formal MML likelihood basis, observation set, and convergence status. LRT was not computed.",
        call. = FALSE
      )
    } else if (!isTRUE(nesting_audit$eligible)) {
      warning(
        "`nested = TRUE` was requested, but the structural nesting audit did not pass. ",
        nesting_audit$reason,
        " LRT was not computed.",
        call. = FALSE
      )
    } else {
      ll <- tbl$LogLik
      np <- tbl$npar
      if (all(is.finite(ll)) && all(is.finite(np)) && np[1] != np[2]) {
        idx_simple <- which.min(np)
        idx_complex <- which.max(np)
        chi_sq <- 2 * (ll[idx_complex] - ll[idx_simple])
        df_diff <- np[idx_complex] - np[idx_simple]
        if (chi_sq >= 0 && df_diff > 0) {
          p_val <- stats::pchisq(chi_sq, df = df_diff, lower.tail = FALSE)
          lrt <- tibble(
            Simple   = tbl$Label[idx_simple],
            Complex  = tbl$Label[idx_complex],
            ChiSq    = chi_sq,
            df       = df_diff,
            p_value  = p_val
          )
          preferred$LRT <- if (p_val < 0.05) tbl$Label[idx_complex] else tbl$Label[idx_simple]
        }
      }
    }
  }

  # -- evidence ratios (pairwise Akaike-weight ratios) --
  evidence_ratios <- NULL
  if (ic_comparable && "AkaikeWeight" %in% names(tbl) && nrow(tbl) >= 2) {
    er_rows <- list()
    for (i in 1:(nrow(tbl) - 1)) {
      for (j in (i + 1):nrow(tbl)) {
        w_i <- tbl$AkaikeWeight[i]
        w_j <- tbl$AkaikeWeight[j]
        er <- if (is.finite(w_i) && is.finite(w_j) && w_j > 0) {
          w_i / w_j
        } else {
          NA_real_
        }
        er_rows[[length(er_rows) + 1]] <- tibble(
          Model1        = tbl$Label[i],
          Model2        = tbl$Label[j],
          EvidenceRatio = er
        )
      }
    }
    evidence_ratios <- bind_rows(er_rows)
  }

  comparison_basis <- list(
    same_method = same_method,
    all_mml = all_mml,
    same_nobs = same_nobs,
    same_data = same_data,
    all_converged = all_converged,
    ic_comparable = ic_comparable,
    nested_requested = isTRUE(nested),
    nesting_audit = nesting_audit
  )

  out <- list(
    table           = tbl,
    lrt             = lrt,
    evidence_ratios = evidence_ratios,
    preferred       = preferred,
    comparison_basis = comparison_basis
  )
  class(out) <- c("mfrm_comparison", class(out))
  out
}

#' @export
summary.mfrm_comparison <- function(object, ...) {
  out <- list(
    table           = object$table,
    lrt             = object$lrt,
    evidence_ratios = object$evidence_ratios,
    preferred       = object$preferred,
    comparison_basis = object$comparison_basis
  )
  class(out) <- "summary.mfrm_comparison"
  out
}

#' @export
print.summary.mfrm_comparison <- function(x, ...) {
  cat("--- MFRM Model Comparison ---\n\n")

  # -- main comparison table --
  tbl <- x$table
  # Format weight columns with 4 decimal places for readability
  fmt_tbl <- as.data.frame(tbl)
  weight_cols <- intersect(c("AkaikeWeight", "BICWeight"), names(fmt_tbl))
  for (wc in weight_cols) {
    fmt_tbl[[wc]] <- ifelse(is.na(fmt_tbl[[wc]]), NA_character_,
                            sprintf("%.4f", fmt_tbl[[wc]]))
  }
  delta_cols <- intersect(c("Delta_AIC", "Delta_BIC"), names(fmt_tbl))
  for (dc in delta_cols) {
    fmt_tbl[[dc]] <- ifelse(is.na(fmt_tbl[[dc]]), NA_character_,
                            sprintf("%.2f", fmt_tbl[[dc]]))
  }
  print(fmt_tbl, row.names = FALSE)

  if (!isTRUE(x$comparison_basis$ic_comparable)) {
    cat("\nInformation-criterion ranking was suppressed because the models do not share\n")
    cat("a comparable formal MML likelihood basis, observation set, and convergence status.\n")
  }

  if (!is.null(x$lrt)) {
    cat("\nLikelihood-ratio test:\n")
    cat(sprintf("  Chi-sq = %.3f, df = %d, p = %.4f\n",
                x$lrt$ChiSq[1], x$lrt$df[1], x$lrt$p_value[1]))
    cat(sprintf("  %s vs %s\n", x$lrt$Simple[1], x$lrt$Complex[1]))
  } else if (isTRUE(x$comparison_basis$nested_requested)) {
    cat("\nLikelihood-ratio test was not reported.\n")
    audit <- x$comparison_basis$nesting_audit %||% list()
    if (!is.null(audit$reason) && nzchar(audit$reason)) {
      cat("  Nesting audit:", audit$reason, "\n")
    }
  }

  # -- evidence ratios --
  if (!is.null(x$evidence_ratios) && nrow(x$evidence_ratios) > 0) {
    cat("\nEvidence ratios (Akaike weights):\n")
    er <- x$evidence_ratios
    for (k in seq_len(nrow(er))) {
      ratio_str <- if (is.finite(er$EvidenceRatio[k])) {
        sprintf("%.2f", er$EvidenceRatio[k])
      } else {
        "NA"
      }
      cat(sprintf("  %s / %s = %s\n",
                  er$Model1[k], er$Model2[k], ratio_str))
    }
  }

  if (length(x$preferred) > 0) {
    cat("\nPreferred model:\n")
    for (nm in names(x$preferred)) {
      cat(sprintf("  By %s: %s\n", nm, x$preferred[[nm]]))
    }
  }
  invisible(x)
}

#' @export
print.mfrm_comparison <- function(x, ...) {
  print(summary(x))
  invisible(x)
}
