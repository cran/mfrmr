# Validate the shared native observation route without computing unused diagnostics.
.mfrm_validate_diagnostic_fit <- function(fit) {
  stop_if_product_slopes(fit, "diagnose_mfrm()")
  if (inherits(fit, c("mfrm_testlet", "mfrm_random_rater"))) {
    stop("Use mfrm_response_diagnostics() for descriptive testlet/shared-rater posterior predictive residuals. Ordinary diagnose_mfrm() indices are unsupported; numerical checks are not model-fit diagnostics.", call. = FALSE)
  }
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  if (inherits(fit, "mfrm_imported_fit")) {
    stop("Native response-level diagnostics require a fit_mfrm() result. For imported fit statistics, re-import the source fit with compute_fit = TRUE.", call. = FALSE)
  }
  fit$config$method <- public_mfrm_method_label(fit$config$method)
  fit
}

#' Fit many-facet ordered-response models with a flexible number of facets
#'
#' Estimate effects on ordered ratings, such as rater severity and criterion
#' difficulty, and Person abilities where the estimator supplies them. A facet is one such
#' source of variation; its levels are the individual raters or criteria.
#' Each data row is one rating. Supply its column names with `person`,
#' `facets`, and `score`, as in the complete example below.
#' The default is `method = "MML"` (marginal maximum likelihood).
#' The `RSM` / `PCM` branches are the package's many-facet Rasch-family
#' reference route. `GPCM` adds positive, level-specific discriminations to
#' one facet, or provisionally to two facets under MML. With one family,
#' MML permits a different facet to supply category steps; JML requires the
#' same facet for both roles. Start with "Choose a model and its population"
#' below; the model-specific sections explain uncertainty and comparisons.
#' In the example, `ratings` stores the data and `fit` stores the fitted model.
#' The loading line reads an existing packaged synthetic dataset; you do not
#' need to generate ratings or create an empty object first.
#' Quoted column names such as `"Person"` must match the data, including case.
#' For your own CSV, see the "Use your own CSV" section of
#' `vignette("mfrmr-workflow", package = "mfrmr")`. If the vignette is not
#' installed, [mfrmr_workflow_methods] and [describe_mfrm_data()] explain the
#' input checks. Pass the reviewed rating data frame to `data`, not the object
#' returned by `describe_mfrm_data()`.
#' Before adapting the example, read "Check defaults before adapting an example"
#' in [mfrmr_workflow_methods]. The score ladder, missing-row handling and
#' population distribution are analysis choices even when arguments are omitted.
#'
#' @param data A data.frame in long format with one row per observed rating
#'   event. Grouped or rowwise tibbles are analyzed as one table; dplyr grouping
#'   metadata is ignored. The named model columns define the design.
#' @param person Column name for the person (character scalar).
#' @param facets Character vector of facet column names. Facets cannot be named
#'   `Person`, `Score`, `Weight`, `score_k`, `theta`, `steps`, `interactions`,
#'   `log_slopes`, `beta` or `log_sigma2`, which are reserved data/parameter
#'   names. Rename such facet columns before analysis. Similar names such as
#'   `beta_group` are supported and matched exactly. Data-column names remain
#'   available for their intended roles, such as `person = "Person"`,
#'   `score = "Score"` and `weight = "Weight"`.
#' @param score Column name for the observed ordered category score. Values
#'   must be coercible to real integer category codes within R's integer range.
#'   Shared preparation rounds values within `1e-6` of an integer to the
#'   nearest category, recording normalization in `fit$prep$preparation_notes`;
#'   other fractional values are rejected. Corrected JML and two-family GPCM
#'   require exact integer inputs; normalize roundoff explicitly before those
#'   specialized fits. Binary `0/1` or `1/2` responses are supported as the ordered
#'   two-category special case. When `keep_original = FALSE`, unused
#'   intermediate categories are collapsed to a contiguous internal scale and
#'   the mapping is recorded in `fit$prep$score_map`. If `rating_min` /
#'   `rating_max` are supplied and the observed scores are a contiguous
#'   subset of that range (for example a 1-5 scale with only 2-5 observed),
#'   the supplied full range is retained so zero-count boundary categories
#'   remain visible in the data-support review. A zero-count boundary category
#'   is retained as review evidence. With `keep_original = TRUE`, however, an
#'   unobserved internal category in a polytomous fitted ladder creates an
#'   unsupported adjacent-step contrast and fitting stops before optimization
#'   with a structured category-support error.
#'   The explicit Jeffreys route below retains empty categories under a
#'   different estimating equation; it does not certify ordinary finite MLEs.
#' @param rating_min Optional minimum category value. Supply this with
#'   `rating_max` when the intended score scale includes unobserved boundary
#'   categories.
#' @param rating_max Optional maximum category value. Supply this with
#'   `rating_min` when the intended score scale includes unobserved boundary
#'   categories.
#' @param weight Optional weight column name. Numeric weights retain their
#'   supplied precision; character/factor weights are parsed from their labels.
#'   Retained weights must have a finite sum and are not automatically rescaled.
#' @param category_policy Optional explicit category choice: `"collapse"`
#'   maps gaps in the observed categories to consecutive scores;
#'   `"preserve"` keeps the intended ladder, declared with `rating_min` and
#'   `rating_max`. This changes the fitted category steps, not just labels.
#'   `NULL` (default) uses `keep_original`, whose default is `FALSE`
#'   (`"collapse"`). Supplying both choices is allowed only when they agree.
#'   Preservation does not estimate unsupported steps: fitting stops if a
#'   retained internal category has no observations. Use the same policy in
#'   [describe_mfrm_data()] and [review_mfrm_anchors()].
#' @param keep_original Logical. `FALSE` (the current default) collapses
#'   non-consecutive observed categories to a contiguous internal scale and
#'   records the mapping in `fit$prep$score_map` (the downstream Count = 0
#'   rows are consequently absent). `TRUE` preserves the declared scale so
#'   unused intermediate categories remain visible in
#'   [rating_scale_table()] and APA outputs, which is recommended for
#'   preserving the rubric. This changes which category steps are fitted, not
#'   just their displayed labels. Fitting stops if a retained internal category
#'   has no observations; reviewing or revising that ladder is a substantive
#'   decision, not a formatting option.
#'   Retained for compatibility; `category_policy = "preserve"` or
#'   `"collapse"` makes the choice explicit in new code.
#' @param missing_codes Optional pre-processing step that converts sentinel
#'   missing-code values to `NA` before any downstream logic. One of:
#'   \itemize{
#'     \item `NULL` (default): no recoding; strictly backward-compatible.
#'     \item `TRUE` or `"default"`: FACETS / SPSS / SAS convention set
#'       (`"99"`, `"999"`, `"-1"`, `"N"`, `"NA"`, `"n/a"`, `"."`, `""`)
#'       on the score column only. Person and facet identifiers are preserved
#'       because short codes such as `"N"` can be legitimate labels.
#'     \item Character vector: an explicit code set, e.g.
#'       `c("99", "999", ".a")`, applied across the person, facet, and
#'       score columns.
#'   }
#'   Replacement counts are recorded in `fit$prep$missing_recoding` and
#'   surfaced by [build_mfrm_manifest()]. Equivalent to calling
#'   [recode_missing_codes()] manually before the fit.
#' @param model `"RSM"` (default), `"PCM"`, or `"GPCM"`.
#' @param method `"MML"` (default), `"JML"`, or `"CML"`. `"JMLE"` is accepted as a
#'   backward-compatible alias for the same joint-maximum-likelihood path.
#'   Native CML supports RSM/PCM with preserved categories and unit weights;
#'   see "Conditional maximum likelihood" for its distinct result class.
#' @param step_facet Facet whose levels receive separate step parameters in
#'   `PCM` and `GPCM`. Supply exactly one facet name for a final analysis. If it
#'   is omitted for `PCM`, mfrmr uses a unique item-like facet name (for
#'   example, `Item`, `Task`, or `Criterion`) when available; otherwise it
#'   retains the first-facet fallback with a warning. `GPCM` always requires an
#'   explicit value. CML also requires an explicit PCM step owner.
#'   This argument is not used by `RSM`, which has one shared
#'   set of rating-scale thresholds.
#' @param slope_facet Slope facet for the `GPCM` branch. mfrmr
#'   estimates one positive slope for every level of this designated facet.
#'   Thus `slope_facet = "Criterion"` gives criterion-specific slopes, whereas
#'   `slope_facet = "Rater"` gives rater-specific slopes. For one slope family,
#'   MML allows a distinct `step_facet`;
#'   JML requires `slope_facet == step_facet`. Criterion and rater slope
#'   blocks can also be estimated together with an ordered pair of column
#'   names; see "Two slope families" below for the explicit MML contract.
#'   With one family, slopes are identified on the log scale with their
#'   geometric mean fixed to 1, so the table reports relative discrimination
#'   across the selected facet's levels rather than unrelated absolute weights.
#' @param gpcm_mml_identification Scale-identification convention for
#'   `model = "GPCM"` with `method = "MML"`. `"free_population"` (the
#'   default) estimates an intercept-only person distribution
#'   \eqn{N(\beta_0,\sigma^2)} when `population_formula` is omitted, while
#'   retaining geometric-mean-one relative slopes. This restores the common
#'   discrimination degree of freedom used by a conventional fixed-latent-
#'   variance GPCM; in the documented item-only overlap it is a one-to-one
#'   reparameterization of ConQuest `scoresfree` GPCM. An explicitly supplied
#'   `population_formula` is retained under this convention.
#'   With one slope family, `"fixed_standard_normal"` is the legacy restricted
#'   branch: it requires
#'   `population_formula = NULL`, fixes the person distribution to
#'   \eqn{N(0,1)}, and also fixes the slope geometric mean to one. The latter
#'   is then a substantive relative-discrimination restriction rather than an
#'   identification requirement. With two families, fixed N(0,1) instead
#'   accompanies a geometric-mean-one first family and free second-family slopes.
#'   This argument does not change JML, whose
#'   geometric-mean-one slope constraint is required to identify its freely
#'   estimated person coordinates.
#' @param facet_interactions Optional confirmatory two-way interaction terms
#'   between non-person facets, supplied as explicit character terms such as
#'   `"Rater:Criterion"` or as a list of length-two character vectors. These
#'   interactions are estimated simultaneously as fixed effects in `RSM` and
#'   `PCM` fits. Person-involving interactions, higher-order interactions, and
#'   random-effect interaction terms are outside the current scope.
#' @param min_obs_per_interaction Minimum weighted observations recommended for
#'   each interaction cell. Cells below this value are flagged in
#'   `interaction_effect_table()` and handled according to `interaction_policy`.
#' @param interaction_policy How to handle sparse interaction cells:
#'   `"warn"` (default), `"error"`, or `"silent"`.
#' @param anchors Optional table with `Facet`, `Level`, and `Anchor` columns.
#'   Each retained row fixes a location in internal logit coordinates; it is
#'   not merely a declaration that an element is common. Person anchors require
#'   JML; MML rejects them. See **Anchors and scale constraints** below.
#' @param group_anchors Optional group-mean constraint table with facet, level,
#'   group, and target-value columns. Its use requires a defensible external
#'   assumption about the group target; it does not create observed overlap.
#'   Use columns `Facet`, `Level`, `Group`, and `GroupValue`. Person group
#'   constraints require JML.
#' @param noncenter_facet One facet to leave non-centered (default `"Person"`).
#'   Other facets use the convention selected by `centering_scope`.
#'   To estimate a cohort mean while leaving an anchored RSM/PCM facet
#'   uncentered, also set `population_formula = ~ 1` and `population_sd`;
#'   see **Anchor origin and cohort means**.
#' @param population_sd Optional positive residual population SD, fixed during
#'   estimation. `NULL` preserves the existing behavior: an active population
#'   model estimates its residual variance; RSM/PCM without a population
#'   formula uses N(0,1).
#'   A numeric value requires `population_formula`, full direct RSM/PCM MML,
#'   and no interactions. It does not activate a population model implicitly.
#'   Supplied SD and anchor values are treated as fixed, without source-error
#'   propagation. GPCM and other estimation engines do not support this option.
#' @param centering_scope `"facet"` (default) centers each declared facet across
#'   all its levels, including direct anchors and group totals. `"remainder"`
#'   reproduces the earlier convention that centers only levels outside direct
#'   and group constraints. For ordinary additive RSM/PCM JML,
#'   `"anchor_origin"` lets direct or group anchors replace centering on their facet;
#'   select an unanchored, active `noncenter_facet`. Other facets remain centered.
#'   This option excludes comparison-only fitting. None of
#'   these conventions changes the population model or category steps.
#'   Unsupported specialized routes reject an explicit setting.
#'   See **Anchors and scale constraints**.
#' @param fit_scope `"full"` (default) fits the ordinary model and retains its
#'   identification guards. `"comparisons"` explicitly fits only the targets in
#'   `comparison_review`; see **Comparison-only estimation**.
#' @param comparison_review For `fit_scope = "comparisons"`, the result of
#'   [describe_mfrm_data()] with `comparisons`, or its `$comparison_review`.
#'   Pass the review's `fit_arguments` and the same original ratings.
#' @param dummy_facets Facets to fix at zero.
#' @param repeated_ratings Treatment of repeated Person-by-facet cells.
#'   `"review"` (default) retains the duplicate-cell warning and ordinary
#'   MML/JML input-review state. Use `"independent"` only to declare that the
#'   rows are distinct, conditionally independent rating events under the
#'   specified model. Every row and its weight are retained; this argument
#'   does not aggregate scores, add a facet, change the likelihood or establish
#'   independence from the data. Other input, numerical and inference checks
#'   still apply. The declaration is saved in preparation notes, configuration
#'   and replay settings. Corrected JML already requires conditionally
#'   independent repeated ratings; the argument records that assumption
#'   explicitly there. It is not supported by the two-family GPCM route.
#' @param positive_facets Facets with positive orientation.
#' @param anchor_policy How to handle anchor-review issues: `"warn"` (default),
#'   `"error"`, or `"silent"`.
#' @param min_common_anchors Minimum directly anchored levels per non-Person
#'   facet used in the package's local count recommendation. This does not
#'   verify cross-run element identity, invariance, or empirical connectedness.
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in anchor-review recommendations.
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in anchor-review recommendations.
#' @param quad_points Integer number of Gauss-Hermite quadrature points
#'   used for MML integration over the person distribution. The default is
#'   `31`. Useful accuracy/runtime settings are:
#'   \tabular{ll}{
#'     `7`  \tab lightweight screening run; information-criterion deltas,
#'                 weights, preferences, and LRT are disabled. Helpers such as
#'                 [predict_mfrm_population()] and
#'                 [reference_case_benchmark()] use this value. \cr
#'     `15` \tab intermediate review run when runtime matters; automatic
#'                 model ranking remains disabled. \cr
#'     `31` \tab package default and a starting grid, not evidence by itself
#'                 that numerical integration is adequate. \cr
#'     `61+` \tab user-selected denser grids for same-data sensitivity review;
#'                 no single order is adequate for every response pattern.
#'   }
#'   Quadrature adequacy depends on the fitted distribution and score support.
#'   Orders whose weights cannot all be represented as positive finite doubles
#'   produce an error; increasing the order indefinitely is not supported.
#'   Use [mml_quadrature_sensitivity()] to compare the same model and data on
#'   user-selected grids before portable calibration and whenever numerical
#'   movement could affect a consequential result. The helper reports
#'   continuous differences without choosing a cutoff. Raw AIC/BIC/SABIC
#'   remain visible below 31 points for diagnosis, but [compare_mfrm()] fails
#'   closed rather than turning a screening/review grid into automatic
#'   selection.
#' @param maxit Computational ceiling on optimizer iterations. The default is
#'   `400`; for two-family EM this counts outer EM iterations.
#'   Two-family adaptive initialization also allows this many EM seed iterations;
#'   it then applies separately to every direct optimizer stage at each start.
#'   For corrected JML it applies separately to each recorded root-solving stage.
#'   This is not a convergence criterion or a model-selection control:
#'   a fit that reaches the ceiling remains non-ready until the common
#'   convergence and terminal-gradient checks pass. Smaller values used in
#'   brief examples are for demonstration only and should not be copied into a
#'   substantive analysis without an explicit computational protocol.
#' @param reltol Portable tolerance setting for the initial optimizer stage.
#'   For two-family fixed-grid EM, omit this argument and use `em_score_tol`.
#'   Adaptive two-family direct MML uses `reltol` and the direct optimizer's
#'   gradient check. New two-family adaptive fits cap the ordinary total-gradient
#'   review tolerance at `N * 1e-6`, where N is the number of fitted Persons;
#'   a stricter ordinary tolerance is retained. This aligns optimizer finishing
#'   with the subsequent per-Person score check. The rule is recorded in the
#'   fitted controls; older saved fits retain their recorded review convention.
#'   Curvature, derivative and integration checks remain separate.
#'   The default is `1e-9`. For BFGS this is passed as `reltol`; for L-BFGS-B
#'   it is mapped to `factr` and `pgtol`, whose actual values are recorded in
#'   the fit. When this setting is at least as strict as the public default
#'   (`reltol <= 1e-9`), optimizer code zero followed by a failed common
#'   terminal-gradient review triggers a bounded warm-started polish ladder.
#'   The best non-worsening stage under the recorded selection rule is
#'   retained. Requested and selected-stage settings remain in `fit$summary`,
#'   and the complete stage history remains in `fit$opt$optimizer_polish`.
#'   Ordinary additive RSM/PCM JML without shrinkage can finish a remaining
#'   stall with at most three sparse Newton steps from the exact likelihood
#'   information. Each step must reduce the gradient without likelihood loss
#'   beyond roundoff, pass its linear-solve check and move no coordinate by more
#'   than one logit. Extreme Persons, nonpositive information and failed
#'   proposals retain their diagnostic state. The final stage is labelled
#'   `jml_information_polish`, with method `Newton`, `InformationSteps` and
#'   `InformationSolveResidual`. This does not correct JML bias or qualify
#'   interval coverage; the original terminal-gradient criterion still applies.
#'   If ordinary polishing still stalls, fixed-grid RSM/PCM/GPCM MML fits
#'   (including estimated population parameters) or adaptive GPCM MML fits
#'   with a fixed population can use one local curvature step to restart the
#'   selected optimizer, provided there are at most 64 free parameters. The step requires
#'   positive-definite, well-conditioned curvature, a smaller gradient and
#'   an objective that does not worsen beyond floating-point roundoff.
#'   The original convergence and terminal-gradient criteria still apply;
#'   failed proposals retain their reasons in the stage history.
#'   Fixed-grid GPCM MML fits with at most 64 free parameters also check
#'   numerical curvature after optimizer code zero. Negative curvature can
#'   trigger up to three BFGS restarts in rescaled search coordinates, even
#'   when the raw gradient is small. Each restart uses the requested `maxit`
#'   ceiling. A replacement must pass the original gradient/convergence checks,
#'   have no detected negative curvature and not worsen the objective beyond
#'   roundoff. Failed recovery retains the estimate with a numerical warning.
#'   For fixed-population, fixed-grid RSM/PCM interaction fits, a remaining
#'   gradient stall also permits one L-BFGS-B polish after subtracting a
#'   constant from the objective (`centered_objective_polish`). This removes
#'   the large likelihood constant from relative stopping, without changing
#'   derivatives or precision requirements. Stored likelihoods retain their
#'   original scale; a resolved likelihood loss is never accepted.
#'   This changes only the search coordinates, not the model or information
#'   matrix used for inference. It does not establish a global optimum,
#'   adequate quadrature, or valid confidence intervals. Inspect the
#'   `SmallestCurvature`, `CurvatureScale` and `CurvatureReviewError` fields
#'   in the stage history alongside the objective and terminal gradient.
#' @param gpcm_mml_start Initial values for two-family adaptive direct MML,
#'   including its `mml_integration="auto"` refinement stages.
#'   `NULL` (default) selects `"neutral_em"`: optimize from both the neutral vector
#'   and a same-data fixed-grid EM vector, and select the lowest finite adaptive
#'   negative log likelihood, retaining its convergence status. `"neutral"`
#'   retains the earlier neutral-only initialization. Numerical optimizer repairs
#'   still apply; exact historical reproduction requires the original source.
#'   An EM vector is only a start;
#'   its fixed-grid likelihood is never compared with an adaptive likelihood.
#'   The seed uses the requested quadrature order, `maxit` outer iterations,
#'   100 M-step iterations and per-Person score tolerance `1e-6`. A finite seed
#'   need not have converged. Each direct start uses the requested optimizer,
#'   `reltol` and its existing polishing sequence. All starts, errors, warnings,
#'   stage histories, seed trace and elapsed costs are stored in
#'   `fit$opt$mml_initialization`; results and reports include its comparison table.
#'   This comparison does not prove a global optimum or interval validity.
#'   Alternatively, supply a plain named list of at least two finite native
#'   free-parameter vectors. For example, reuse `p <- previous_fit$opt$par`
#'   from the same model coding with `list(previous = p, neutral = 0 * p)`.
#'   The order is first-owner free locations, second-owner locations,
#'   second-owner free steps (within each level), first-owner free log slopes,
#'   and second-owner log slopes, using the fitted level order. Vector lengths
#'   must match; optional parameter names must exactly match the model's native
#'   names. Reported positive slopes or expanded constrained parameters cannot
#'   be supplied in place of native coordinates. Explicit starts replace the
#'   default start roster and do not run an EM initialization. Their values are
#'   reevaluated and optimized under the requested adaptive objective; no
#'   likelihood, convergence flag, covariance or scoring eligibility is imported
#'   from a previous calculation. The roster and vectors remain in the saved
#'   initialization record and replay inputs, including quadrature refits.
#' @param em_score_tol Stopping tolerance for the two-family GPCM MML-EM route
#'   only. `NULL` selects `1e-6` for that route. Stops when the largest absolute
#'   derivative of the negative marginal log likelihood, divided by the number
#'   of Persons, is at most this value. It does not control interval eligibility.
#'   The ascent-checked BFGS M step allows 100 iterations with `reltol = 1e-12`;
#'   both the marginal likelihood and EM auxiliary objective must not decrease
#'   beyond numerical roundoff. If BFGS stops before the requested score
#'   accuracy, usable positive curvature can refine the auxiliary objective
#'   with the same E-step counts. The ascent checks still apply. Dense
#'   curvature refinement is limited to 1--64 free parameters; otherwise
#'   the ordinary M-step proposal is retained. The saved `fit$opt$em_trace`
#'   records refinement attempts and failures;
#'   `mstep_score` is the proposal's auxiliary-objective score before any
#'   step halving. Other model routes require `NULL`.
#' @param optimizer Direct-optimization method. `"auto"` (default) uses the
#'   BFGS M step for the two-family EM route. For other routes it uses the
#'   limited-memory `"L-BFGS-B"` method for MML and for larger JML parameter
#'   vectors (at least 200 free parameters), while retaining BFGS for smaller
#'   JML fits. Use `"BFGS"` or `"L-BFGS-B"` to request one method explicitly.
#'   The method actually used is recorded in `fit$summary$OptimizerMethod`.
#'   For L-BFGS-B, inspect `OptimizerFactr` and `OptimizerPgtol` rather than
#'   interpreting `EffectiveReltol` as a native [stats::optim()] control.
#' @param mml_integration Numerical integration for MML: `"fixed"` (default)
#'   uses a common Gauss-Hermite grid relative to each person's normal prior;
#'   `"adaptive"` centers and scales the grid at each person's posterior mode
#'   on every objective evaluation. Adaptive fitting uses the gradient of the
#'   moving-node objective, including changes in its center and width. It
#'   requires `method = "MML"`, `mml_engine = "direct"`, and no checkpoint.
#'   The choice also controls the fitted likelihood, Hessian, posterior person
#'   summaries, fitted-object scoring, and supported portable calibration.
#'   For these two modes `quad_points` remains the order; neither guarantees
#'   adequate accuracy at a particular order. Compare orders before reporting.
#'   During adaptive optimization, trial values that cause numerical overflow
#'   are rejected so the search can shorten its step. Starting values must still
#'   be evaluable, and final estimates must pass the usual convergence checks;
#'   rejecting a trial does not establish a parameter boundary or impose a
#'   discrimination upper limit.
#'   Adaptive fits with nonlinear coordinates have no completed probability-
#'   map identification audit; GPCM also lacks an adaptive slope-boundary audit.
#'   Their overall `InferenceReady` flag remains false. Supported native-location
#'   and RSM/PCM population intervals have separate source and numerical checks;
#'   see the MML integration workflow section below and [mfrm_facet_intervals()].
#'   Available intervals retain their asymptotic assumptions and coverage cautions.
#'   ConQuest exports currently require fixed integration.
#'   `"auto"` first fits the fixed grid at `quad_points`, then checks every
#'   Person's log marginal likelihood, EAP and posterior SD against two adaptive
#'   orders. It also checks the change in the full marginal gradient and its
#'   absolute maximum at the retained and finer rules. If any criterion fails,
#'   it warm-starts adaptive refits at the orders in `mml_accuracy_control`.
#'   This fixed-to-adaptive sequence applies to ordinary full direct MML
#'   (RSM, PCM or one-slope-family GPCM), with no checkpoint. Population regression, anchors and interactions retain
#'   their existing model restrictions. On success the resolved fixed/adaptive
#'   rule and order control all subsequent fit calculations; inspect
#'   `fit$config$mml_accuracy` for controls, stage history and final Person
#'   checks. These checks concern the fitted response patterns; scoring new
#'   patterns reuses the resolved rule and needs its own integration review.
#'   The initial fixed rule is retained if its checks pass. Requested
#'   settings are preserved for replay; explicit sensitivity refits use the
#'   resolved rule at each requested order. If optimization fails or all
#'   refinement orders are exhausted, an error of class
#'   `mfrmr_mml_accuracy_error` carries `history`, `control` and the last
#'   parameter vector instead of returning a fit.
#'   In current errors (`provenance_version = 2`), `last_candidate` holds the
#'   last returned optimizer vector, its parameter-block names and sizes, its
#'   matching `config`, and its stage diagnostics. These are free optimization
#'   coordinates, not a table of facet estimates or a precision-qualified fit.
#'   `last_attempt` records the final attempted stage and failure phase; it may
#'   describe a later failed optimization. `last_integration` and
#'   `last_quad_points` now describe `last_parameters`, not that later attempt.
#'   If no optimizer returned, `last_candidate` is NULL, `parameter_source` is
#'   `"initial_values"`, and the integration aliases are NA. Initial values are
#'   retained separately as `initial_parameters`. Earlier saved errors lack
#'   this provenance version; their integration aliases can describe the failed
#'   attempt rather than the retained vector. Inspect their stage history.
#'   Capture the condition with `tryCatch(fit_mfrm(...),
#'   mfrmr_mml_accuracy_error = function(e) e)` and inspect both records.
#'   Saving the condition preserves the values and reasons without converting
#'   unsuccessful accuracy checks into a successful model fit.
#'
#'   For two-family GPCM, `"auto"` instead preserves the requested engine:
#'   `mml_engine="em"` uses fixed integration and `mml_engine="direct"`
#'   uses adaptive integration at every stage. Each refit retains the existing
#'   neutral EM start or the specified adaptive `gpcm_mml_start` policy;
#'   no cross-order warm start or engine switch is introduced. The default
#'   hard ceiling is 241 points, including the `2*q-1` reference checks.
#'   Person log marginal likelihood, EAP and posterior SD differences must
#'   pass, together with the existing two-family local-information and
#'   stationarity checks described in [confint.mfrm_fit()]. This conservative
#'   automatic route therefore requires more than delivery of a point estimate.
#'   A sole quadrature-sensitivity failure or inaccurate Person moments causes
#'   refinement; another failed information, convergence or model check stops
#'   with its own reason. It does not automatically diagnose a boundary.
#'   In two-family errors, `last_candidate$accuracy` retains available checks
#'   and Person comparisons. If fitting failed before a candidate was returned,
#'   `parameter_source="unavailable"` and `last_parameters=NULL`; there is no
#'   single initial vector to report for a multi-start policy. No two-family
#'   global optimum or interval-coverage guarantee follows from auto success.
#'   The latest successful automatic fit can reuse its complete qualification
#'   for subsequent two-family intervals in the same R session, provided the
#'   entire fit, information-memory budget and numerical backend are unchanged.
#'   Source and result retention is bounded by that budget. Warning-producing
#'   or refused qualification is not reused.
#'   This transient cache is not saved in the fit; a new session recomputes.
#' @param mml_accuracy_control Named list used only with
#'   `mml_integration = "auto"`. To specify only a ceiling, use
#'   `list(max_points=121)`. It bounds every quadrature rule used by automatic
#'   fitting and its acceptance checks, including reference rules; it does not
#'   bound optimizer iterations, total runtime or later separately requested
#'   inference/scoring. Orders increase automatically and stop at the first
#'   accepted stage. `max_points` must be an integer of at least 7 for ordinary
#'   MML, or 3 for two-family GPCM; requested rules must be representable.
#'   Unused finer rules are not constructed. A ceiling and an explicit
#'   `adaptive_points` vector cannot both be supplied. An explicitly supplied
#'   `quad_points` must fit within the ceiling (including its finer reference
#'   for two-family fitting); an omitted starting order is reduced if needed.
#'   For example, the two-family ceiling 121 permits fitting/checking orders
#'   31/61 and 61/121. A ceiling of 100 permits 31/61 and 50/99, preserving the
#'   same reference-order formula. Reaching the ceiling without passing all
#'   criteria raises `mfrmr_mml_accuracy_error`, retaining the attempted rules
#'   and source-bound candidate; it does not relax tolerances or return a
#'   precision-qualified fit.
#'
#'   For ordinary MML without a ceiling, the retained defaults are
#'   `adaptive_points = c(31, 61, 121)`
#'   (strictly increasing integer refit orders, each checked against `2*q+1`),
#'   `loglik_tol = 1e-7` (maximum absolute log-marginal difference per Person),
#'   `posterior_tol = 1e-6` (maximum absolute EAP or posterior-SD difference,
#'   on the ability scale), and `gradient_tol = 1e-4` (maximum absolute gradient
#'   or gradient difference, in free optimization coordinates for the summed
#'   negative log likelihood). At the initial fixed fit the first adaptive
#'   order and its `2*q+1` refinement are used. All criteria and optimizer
#'   convergence must pass. These are empirical numerical tolerances, not
#'   certified error bounds, relative tolerances, a bound on parameter movement,
#'   or a guarantee of a global optimum or interval coverage. Weak curvature
#'   can amplify small gradient errors. Hessian precision and formal inference
#'   require their separate checks; automatic integration does not qualify
#'   previously unavailable identification or inference. `maxit` applies to
#'   each fit, so refinement increases total computation.
#'   Two-family controls instead default to `max_points=241`,
#'   `loglik_tol=1e-7` and `posterior_tol=1e-6`, with `2*q-1` reference rules.
#'   They do not accept `adaptive_points` or `gradient_tol`: their established
#'   engine-specific stationarity and information criteria are retained.
#'   Passing component-interval checks at an order does not ensure that its
#'   Person likelihood and posterior moments meet these additional tolerances.
#'   A prerequisite check that stops two-family refinement records
#'   `history$Status = "check_refused"`. Inspect the failed rows of
#'   `last_candidate$accuracy$checks`; this status does not identify the cause.
#'   A failed higher-order information check also stops refinement. Earlier
#'   saved histories use `"nonintegration_refusal"` for the same stopping rule;
#'   neither label establishes a boundary or whether more points would help.
#' @param progress Logical; `TRUE` displays a progress indicator with elapsed
#'   time. Two-family GPCM reports EM iterations, adaptive starts, quadrature
#'   stages and precision checks; direct optimization reports objective
#'   evaluations. Other routes may provide only a coarse phase display.
#'   Iteration/stage limits are budgets, not percent completion or a remaining-
#'   time forecast: convergence and required refinement are not known in advance.
#'   Display follows `cli` console settings and is cleared on return, error or
#'   interruption. The default is silent; it does not change numerical controls
#'   or become part of saved calibration/replay settings.
#' @param mml_starts Optional starting-value review for ordinary full RSM/PCM
#'   or single-family GPCM MML. `NULL` keeps the usual single start. A numeric
#'   vector such as `c(0, 0.5, -0.5)` adds each nonzero offset times the
#'   alternating vector `(1, -1, 1, ...)` to the default free coordinates.
#'   These deterministic offsets explore one direction in the fitted free-
#'   coordinate order, not random or exhaustive searches. The default start
#'   is always included. To explore other directions, supply a uniquely named
#'   list of named numeric coordinate values; obtain names from
#'   `fit$opt$start_review$parameter_map$Coordinate` after `mml_starts = 0`.
#'   Unlisted coordinates retain their default initial values. Anchored and
#'   dependent coordinates are expanded from the existing constraints.
#'   All attempts share the same fixed or adaptive integration rule and order;
#'   automatic integration, checkpoints and two-family GPCM are unsupported.
#'   Each attempt uses the requested engine, with existing fallbacks recorded.
#'   The smallest reevaluated negative log likelihood is retained even when
#'   that attempt is unfinished; a worse converged attempt is not substituted.
#'   Inspect `summary(fit)$starting_values` and `fit$opt$start_review` for all
#'   starts, optimizer outcomes, warnings, errors and coordinate differences.
#'   If all attempts fail, the error condition retains `$start_review`.
#'   Missing curvature means it was not evaluated. This review increases fitting
#'   time and does not establish a global optimum, adequate integration or valid
#'   inference. Use [mml_quadrature_sensitivity()] for integration sensitivity.
#'   Same-model replay and resampling retain this policy. Subgroup DFF refits
#'   accept numeric offsets; named starts require an explicit subgroup fit
#'   because linking changes the free-coordinate basis.
#' @param mml_engine MML optimization engine for `method = "MML"`:
#'   `"direct"` (default) uses the selected direct optimizer on the marginal
#'   log-likelihood,
#'   `"em"` uses an EM loop for `RSM` / `PCM` with `population = NULL`, including
#'   model-estimated `facet_interactions`, or the explicitly scoped two-family
#'   GPCM route described below. `"hybrid"` uses EM as a warm start before the
#'   direct optimizer. Ordinary EM/hybrid requires fixed quadrature. Numerical
#'   M steps implement generalized EM, checking both the auxiliary and marginal
#'   log likelihood for non-decrease beyond roundoff before accepting an update;
#'   failed checks stop without replacing a valid checkpoint. This is not a
#'   block ECM algorithm or a guarantee of global convergence or faster fitting.
#'   Use `mml_starts` to review specified starts and keep integration sensitivity
#'   separate. For a raw EM fit with a large terminal gradient, use `"hybrid"`
#'   or `"direct"` and review the resulting readiness. Unsupported
#'   one-family combinations fall back to `"direct"` and record that fallback in
#'   `fit$summary`. Direct, hybrid, and EM engines all require the common
#'   terminal-gradient check for the Numerical component of fit readiness; EM
#'   relative log-likelihood convergence alone does not establish numerical
#'   readiness. `InferenceReady` is `TRUE` only when every stored fit-readiness
#'   component passes.
#' @param population_formula Optional one-sided formula for a person-level
#'   latent-regression population model, for example `~ grade + ses`. Latent
#'   regression is implemented only for
#'   `method = "MML"` with a unidimensional conditional-normal population
#'   model. With `NULL`, RSM/PCM MML uses a fixed \eqn{N(0,1)} basis;
#'   default GPCM MML instead estimates an intercept-only normal population.
#'   Version 0.2.4 does not accept arbitrary fixed normal means or standard
#'   deviations: portable RSM/PCM calibration and its anchors are defined on
#'   the standard-normal basis, and a silent change of basis would change the
#'   meaning of those stored values.
#'   A formula-compatible `terms` object returned in
#'   `attr(simulated_data, "mfrm_population_data")$population_formula` is also
#'   accepted, preserving the simulation template's transformation basis,
#'   factor levels and contrasts when refitting generated data.
#' @param person_data Optional one-row-per-person data.frame holding background
#'   variables for `population_formula`. Numeric, logical, factor, ordered
#'   factor, and character predictors are expanded through `stats::model.matrix()`;
#'   categorical xlevels and contrasts are stored for replay and scoring.
#'   Prediction-aware transformations such as `scale()`, `poly()` and
#'   `splines::ns()` retain their training basis when scoring new persons.
#'   Required when `population_formula` is supplied.
#' @param person_id Optional person-ID column in `person_data`. Defaults to
#'   `person` when that column exists in `person_data`.
#' @param population_policy How missing background data are handled for a
#'   latent-regression fit. `"error"` (default) requires complete person-level
#'   covariates; `"omit"` fits the model on the complete-case subset and records
#'   omitted persons / omitted response rows in the returned `population`
#'   metadata while retaining the observed-person-aligned pre-omit table for
#'   replay/export provenance.
#' @param facet_shrinkage Character. `"none"` (default) keeps the unshrunk
#'   fixed-effects estimates. `"empirical_bayes"` applies a post-hoc
#'   zero-centered normal-model adjustment to each non-person facet
#'   (Efron & Morris, 1973); `fit$facets$others` gains `ShrunkEstimate`,
#'   `ShrunkSE`, and `ShrinkageFactor` columns, and
#'   `fit$shrinkage_report` records the per-facet prior variance and
#'   sum of retained weights (not model degrees of freedom). This descriptive
#'   approximation has no general James-Stein risk guarantee; see
#'   [apply_empirical_bayes_shrinkage()]. `"laplace"` is retained as a compatibility
#'   alias for `"empirical_bayes"`; it does not fit a penalized likelihood.
#' @param facet_prior_sd Optional numeric scalar. When supplied, the
#'   shrinkage prior variance is fixed at `facet_prior_sd^2` instead
#'   of being estimated by method of moments. Useful for eliciting a
#'   prior from domain knowledge or a previous fit.
#' @param shrink_person Logical. When `TRUE` and `facet_shrinkage` is
#'   active, the same descriptive adjustment is applied to ordinary-JML
#'   Person estimates with available sampling SEs. Default `FALSE`.
#'   MML EAP/posterior SDs cannot be used as observations/sampling SEs
#'   for a second normal-model update; use `FALSE` with other estimators.
#' @param attach_diagnostics Logical. When `TRUE`, [diagnose_mfrm()] is
#'   run once after the fit with `residual_pca = "none"`, and the
#'   per-level `SE`, `Infit`, `Outfit`, `InfitZSTD`, `OutfitZSTD`, and
#'   `PtMeaCorr` columns from `diagnostics$measures` are merged onto
#'   `fit$facets$others` (non-person facets) and `fit$facets$person`
#'   (Person rows). This is convenient when downstream code expects a
#'   FACETS Table 7 style facet table with fit statistics in one place,
#'   and lets `summary(fit)` show per-person fit columns alongside the
#'   measure. For person rows, an existing posterior `SE` (typical for
#'   `method = "MML"`) is preserved and the diagnostic `SE` is only
#'   attached when the existing column is empty. Adds diagnostic
#'   runtime (typically +1-2 s on moderate designs) and sets
#'   `fit$config$attached_diagnostics = TRUE`. Default `FALSE`
#'   preserves the minimal `Facet` / `Level` / `Estimate` layout.
#' @param checkpoint Optional `list(file = ..., every_iter = ...)`.
#'   When supplied, the MML EM engine writes its state to `file`
#'   every `every_iter` outer EM iterations using checked same-directory
#'   replacement. If the file already exists when the fit starts, the engine
#'   resumes only when its versioned identity exactly matches the current data,
#'   model, parameter layout, constraints, quadrature, package version, and
#'   engine stage. Legacy, corrupt, or incompatible files fail closed. A
#'   non-converged pure-EM run may resume with a larger `maxit`; a completed
#'   checkpoint cannot re-enter the same iteration boundary. Only the EM engine
#'   (`mml_engine = "em"` or the EM warm-start step of
#'   `mml_engine = "hybrid"`) honours the checkpoint; the direct
#'   `optim()` engine ignores it. Use this to make long MML EM
#'   fits crash-resilient on shared compute environments.
#'
#' @param jml_penalty `"none"` (default) preserves the ordinary estimator.
#'   Explicit `"jeffreys"` with additive RSM/PCM `method = "JML"` selects
#'   canonical information-penalized joint estimation. It returns a separate
#'   `mfrm_jeffreys_fit` with point estimates and numerical checks, without
#'   ordinary SEs or FACETS diagnostics. See **Jeffreys-penalized JML**.
#' @param jml_correction_order `NULL` (default) leaves the usual estimator
#'   unchanged. A positive integer explicitly selects an experimental
#'   profile-score adjustment for shared-owner `model = "GPCM", method = "JML"`.
#'   There is no automatic choice of order. See **Corrected JML** for the
#'   input restrictions, uncertainty meaning and connected outputs.
#' @param jml_correction_sampling Sampling assumption for corrected JML's
#'   Person covariance: `"fixed_rosters"` (default) centers contributions
#'   within each observed assignment pattern; `"random_rosters"` centers
#'   globally when assignment patterns are sampled. This changes covariance,
#'   not the point equation. The fixed-roster calculation needs at least two
#'   Persons per pattern. Both describe new independent Persons, including
#'   variation in their ability composition; neither estimates the variance
#'   from reassessing the same Persons at fixed abilities. Omit this argument
#'   when no correction is requested.
#'
#' @section Choose a model and its population:
#' For the `Rater` and `Criterion` columns in the Examples, the model choices
#' are below. A slope family gives one discrimination per level of its owner;
#' "two families" does not mean two ability dimensions or just two slopes.
#' No Q matrix is required or accepted by this unidimensional API.
#'
#' | Model | Step owner | Slope owner | MML population |
#' | --- | --- | --- | --- |
#' | `RSM` | Shared | Fixed at one | Fixed N(0,1) |
#' | `PCM` | `Criterion` | Fixed at one | Fixed N(0,1) |
#' | One-family `GPCM` | `Criterion` | `Criterion` | Estimated mean and SD |
#' | Two-family `GPCM` | `Rater` | `Criterion`, then `Rater` | Explicit fixed N(0,1) |
#'
#' Supply these owner names through `step_facet` and `slope_facet`; for two
#' families use `slope_facet = c("Criterion", "Rater")`. The first three rows
#' show their default MML populations. Two-family settings must be chosen
#' explicitly; see "Two slope families" below. One-family slopes are relative.
#'
#' Changing only `model` therefore need not hold the population assumptions
#' constant. For RSM/PCM, estimating its mean and SD requires
#' `population_formula = ~ 1` and `person_data = unique(ratings["Person"])`.
#' This is a model choice, not a universal recommendation. Fixing both the
#' population SD and the one-family GPCM slope geometric mean to one is a
#' narrower model, not just an alternative reporting scale.
#'
#' Start with the complete RSM Examples. For PCM, add `step_facet` and change
#' `model`; for one-family GPCM, also explicitly choose `slope_facet` and
#' read the population choice above. Keep the data-column and rubric arguments.
#' The complete two-family call and its narrower output workflow are in
#' `vignette("mfrmr-gpcm-scope", package = "mfrmr")`, under "A provisional
#' two-family workflow". Its zero-based scores and ordered owners require
#' the explicit preparation described there; changing `slope_facet` alone
#' is insufficient.
#'
#' @section Choose an estimator and the result you need:
#' The ordinary worked example uses RSM/MML. Choose the model and estimation
#' assumptions for your question before choosing a table or figure.
#' - **Abilities, a population model or GPCM inference:** use `method = "MML"`
#'   and inspect `summary(fit)$decision`. Its interpretation depends on the
#'   declared ability distribution. Ordinary RSM/PCM fixes N(0,1) by default;
#'   `population_formula = ~ 1` instead estimates its mean and SD. GPCM
#'   estimates them by default. See the population and estimator sections.
#' - **RSM/PCM structural effects without an ability distribution:** use
#'   `method = "CML", category_policy = "preserve"`. Read
#'   `summary(fit)$locations` and `$steps`, then [mfrm_facet_intervals()].
#'   CML returns no Person abilities or population SD; ordinary diagnostics
#'   and consolidated reports do not apply. The conditional-likelihood
#'   section below and the facet-interval vignette give its assumptions.
#' - **Joint point estimates or FACETS comparisons:** explicitly select JML
#'   when that estimand is intended. With few ratings per Person, adding
#'   Persons does not remove incidental-parameter bias. Stored structural
#'   SEs or corrected-JML root uncertainty are not automatically truth intervals.
#'
#' Facet names are your column names, not fixed roles imposed by the package.
#' For example, `facets = c("Judge", "Task")` uses those two columns; it does
#' not require columns called Rater or Criterion. Supply the same names to
#' follow-up functions and use the fitted level labels for contrasts or curves.
#' Begin with the Examples section, then read only the specialist sections
#' for your selected model. [mfrmr_methods_reference] separates measurement,
#' G-theory assessment planning and observed-score agreement questions.
#'
#' @section Choose the arguments by their purpose:
#' - **Identify the ratings:** `data` is the rating table; `person`, `facets`
#'   and `score` are quoted column names, not the values in those columns.
#' - **Declare the rubric:** set `rating_min`, `rating_max` and
#'   `category_policy` consistently with [describe_mfrm_data()]. The default
#'   can collapse gaps; preservation can reveal an unsupported category step.
#' - **Choose the statistical model:** `model`, `step_facet`, `slope_facet`,
#'   anchors and population arguments determine what is estimated. Ordinary
#'   RSM/PCM MML with `population_formula = NULL` fixes N(0,1); the default
#'   GPCM MML instead estimates the normal population mean and variance.
#' - **Control computation:** `quad_points`, `maxit`, `reltol`, `optimizer`
#'   and `mml_engine` govern numerical fitting. Increasing them does not change
#'   the model's support or automatically justify statistical inference.
#' - **Choose follow-up output:** `attach_diagnostics = TRUE` computes and
#'   attaches diagnostics; the default leaves that separate. `summary(fit)`
#'   explains the result's status; [diagnose_mfrm()] reviews response fit.
#'
#' @section Corrected JML:
#' When each Person has few ratings, ordinary JML can retain structural bias
#' even with many Persons. The explicit `jml_correction_order` option adjusts
#' the profile-score equation to address this source of bias. It changes the
#' estimator, not the response model. A higher order is not necessarily more
#' accurate: bias, variance and numerical stability can change in different
#' directions. Choose the order as part of the analysis plan and examine
#' sensitivity; neither AIC nor the smallest RootSE chooses a justified order.
#'
#' This experimental route requires the optional `nleqslv` package, observed
#' integer scores, declared `rating_min` and `rating_max`, additive fixed facets,
#' and a single explicit `step_facet` that also owns the slopes. Facet names are
#' unrestricted. Use `category_policy = "preserve"` to retain the declared
#' score ladder. All facet locations and each step ladder sum to zero; slopes
#' have geometric mean one. Persons are independent sampling units, with
#' conditionally independent ratings within Person. Only observed assignments
#' contribute; unassigned or missing scores are not imputed. Repeated cell rows
#' are treated as independent responses, not correlated repeated measurements.
#' Fixed anchors, nonunit weights, interactions, shrinkage, separate slope/step
#' owners, corrected RSM/PCM, and ordinary optimizer controls are unavailable.
#' Explicit unsupported arguments produce an error rather than being ignored.
#'
#' For an order \eqn{k}, the equation uses \eqn{U_k=(I-K_\beta)^k U_0}, where
#' \eqn{U_0} is the negative profile-likelihood score and \eqn{K_\beta} averages
#' over responses at their profiled Person MLE. This is the MLE plug-in
#' construction discussed by Dhaene and Weidner (2023, section 8.1), applied
#' here to the declared GPCM structure. Exact owner-total enumeration is
#' limited to 5,000 joint states per assignment pattern; larger calculations
#' stop explicitly. No approximate pruning is substituted. The full equation
#' Jacobian and empirical Person contributions give a sandwich covariance,
#' with the selected assignment-sampling assumption. This is not an inverse
#' likelihood Hessian or a general proof of bias removal.
#'
#' For example, `"fixed_rosters"` describes new cohorts with the same number
#' of Persons assigned to each rater/task combination. It does not hold each
#' Person's ability fixed across repeated cohorts. Ability distributions are
#' unspecified and may differ between assignment patterns. With
#' `"random_rosters"`, both Persons and their assignments are sampled from
#' the same joint population; assignment need not be independent of ability.
#' Neither option models shared random raters or dependent Persons.
#' The empirical sandwich has no small-sample degrees-of-freedom correction.
#' Its fixed-roster interpretation requires enough independent Persons within
#' each pattern; two Persons is a computational minimum, not an assurance of
#' accurate uncertainty. Taking `"random_rosters"` solely to obtain a RootSE
#' changes the sampling assumption and is not a repair for sparse information.
#'
#' Use `summary(fit)$tables` for locations, steps, relative slopes and the
#' uncertainty explanation; `include_person = TRUE` adds Person profiles.
#' `RootSE` describes local variation around the adjusted-equation solution,
#' which may remain biased for the true structural parameter. `LogRootSE`
#' describes log-slope variation. Neither is supplied as a standard error
#' with established structural coverage; no confidence limits are constructed.
#' Persons are reprofiled at the adjusted calibration, with infinite limits
#' for extreme scores and no invented Person SE. Numerical attempts remain
#' available in `summary(fit)$attempts` and `fit$jml_adjustment`.
#'
#' A valid point estimate is retained if covariance is unavailable, with its
#' reason. Different accepted roots yield missing primary estimates rather
#' than selection by the ordinary likelihood. Agreement between starting
#' values is a local check, not proof of a unique global solution. A failure
#' remains a reportable fitted object with missing estimates and its attempts.
#'
#' `plot(fit, type = "slopes")`, `"locations"` (select one `facet`) and
#' `"steps"` show point estimates. Choose `style = "distribution"` for an
#' empirical cumulative view; `as_ggplot(plot(fit, draw = FALSE))` is supported.
#' Use `mfrm_results(fit)` then `mfrm_report()` or `export_mfrm_results()` to
#' retain these meanings through reporting and saved replay. For conditional
#' probabilities and descriptive residuals on the observed fitted rows, use
#' [mfrm_response_diagnostics()] and attach the saved result with
#' `mfrm_results(fit, response_diagnostics = result)`. This holds corrected
#' calibration and Person profiles fixed; no posterior averaging, uncertainty
#' intervals or fit cutoffs are supplied. Extreme-score probabilities remain
#' available; zero variance prevents standardization, with separate availability
#' for Infit and Outfit. See that helper's help before interpreting a summary.
#' For new Persons, [predict_mfrm_units()] and [extract_mfrm_calibration()]
#' provide conditional EAP scoring after their corrected-equation source checks.
#' These scores add a separate normal reference prior. Their posterior intervals
#' condition on calibration and do not establish structural coverage or remove
#' residual bias. The portable artifact omits training responses and Person
#' estimates; its own summary and score plots retain the correction order.
#' Ordinary fit tests, Wright/Pathway maps, structural confidence intervals,
#' likelihood ranking and corrected Person ML/WLE remain unavailable.
#' The default MML and uncorrected JML workflows are unchanged.
#'
#' @section Two slope families:
#' To estimate, for example, criterion and assessor discrimination together,
#' use `slope_facet = c("Criterion", "Assessor")`. The order declares roles:
#' the first facet has centered locations and slopes with geometric mean one;
#' the second has free locations/slopes and owns the centered category steps.
#' The product of their slopes multiplies the complete adjacent-category
#' predictor. There is still one ability dimension. Column names and data
#' column order do not select these roles.
#'
#' This route requires exactly those two non-Person facets, `model = "GPCM"`,
#' `method = "MML"`, `gpcm_mml_identification = "fixed_standard_normal"`,
#' either fixed integration with `mml_engine = "em"` or
#' `mml_integration = "adaptive"` with `mml_engine = "direct"`,
#' and both `step_facet` and `noncenter_facet` set to the second slope facet.
#' These options must be chosen explicitly; the one-family defaults are not
#' silently replaced. Use observed, unweighted integer scores on a scale from
#' zero to `rating_max`, without missing IDs, duplicate person-facet rows or
#' surrounding ID whitespace. Categories are not collapsed; if any are absent,
#' declare `rating_max` and `category_policy = "preserve"` so category-support
#' checks can assess that scale. Unobserved assignments are not zero scores.
#' Anchors, population covariates, shrinkage, interactions, positive/dummy
#' facets, checkpoints and automatic diagnostics are not supported in this route;
#' explicitly supplying their arguments, including unused policy controls,
#' produces an error. Fixed-grid EM uses `em_score_tol`; adaptive direct MML
#' uses `reltol` and the direct optimizer's gradient check. Passing the other
#' engine's tolerance is an error. Adaptive integration preserves the same
#' N(0,1) population and response equation while moving each Person's grid.
#' It is direct maximization of the marginal likelihood, not adaptive EM.
#' By default it compares neutral and EM-derived starts through that same
#' adaptive objective; see `gpcm_mml_start`. A better unfinished solution is
#' retained as unfinished, rather than replaced by a worse converged candidate.
#' If a retained starting point is better than all terminal candidates beyond
#' roundoff, the returned fit remains numerically unresolved. A failed alternative
#' is disclosed; if every direct attempt fails, the error carries an
#' `initialization` record. Neither fixed-grid likelihoods nor interval outcomes
#' select a candidate. Extreme slopes and weak information can still prevent
#' inference even after a successful start comparison.
#' Use [mml_quadrature_sensitivity()] to compare refits at different orders;
#' it preserves the chosen integration method, engine and initialization policy.
#'
#' `summary(fit)` and `print(fit)` retain numerical status and the two slope
#' references. [mfrm_curve_intervals()] evaluates provisional category or
#' per-rating information curves at supplied abilities; all their intervals
#' are unavailable. Use `plot(curves)` and
#' `mfrm_results(fit, include = c("fit", "plots"), compute = "never",
#' intervals = list(curves = curves))` to report saved curves. Their values can
#' be inspected even after nonconvergence, but are not qualified estimates.
#' Separately request `confint(fit, parm = "slopes")` for experimental component log-Wald intervals
#' from the full observed marginal information. Their numerical checks do not
#' qualify global identification or sampling coverage; failed checks retain
#' missing bounds. Attach the result alongside curves for saved plots/reports.
#' To profile one component with the fitted integration method, explicitly request
#' `confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t1"))`, replacing the
#' named owner and level with your fitted identifiers. This reoptimizes other
#' coefficients; it is an experimental local interval without established coverage.
#' [mfrm_response_diagnostics()] integrates each Person's ability posterior at
#' the saved calibration for descriptive residuals. Both slopes are retained;
#' unavailable integration stays explicit. There are no reference fit cutoffs.
#' Attach this result through `response_diagnostics` in the results call above.
#' Adaptive two-family log-Wald checks use the moving-node marginal objective
#' and compare adaptive quadrature orders. Posterior residuals preserve the
#' fitted integration method and complete conditioning record, with separate
#' row-wise integration checks. Adaptive profile searches recalculate the nodes,
#' scales and Jacobians at every constrained trial and keep separate checks.
#' [predict_mfrm_units()] separately supplies experimental conditional new-Person
#' EAP and posterior intervals with the retained N(0,1) prior, known levels,
#' unit weights and separate source/batch checks; calibration uncertainty is excluded.
#' [mfrm_facet_intervals()] separately supplies experimental model-based
#' intervals for either owner's locations and within-facet contrasts, retaining
#' slope/step nuisance uncertainty. Failed numerical checks leave missing bounds;
#' location differences need not imply uniform rating differences. Attach these
#' intervals to [mfrm_results()] for plots, reports and exports.
#' Step/curve intervals, sandwich inference, model ranking/LRT, ordinary fit/bias
#' diagnostics and Wright/Pathway plots remain unavailable.
#' [extract_mfrm_calibration()] provides a separately checked portable two-family
#' route. Matching native or portable scores can be attached with
#' `mfrm_results(fit, scores = scores)` without rescoring.
#' Numerical agreement of the fitting implementation is not a
#' general identification, convergence or coverage guarantee. See
#' `vignette("mfrmr-gpcm-scope")` for a complete example and interpretation.
#'
#' @section GPCM model and inference:
#' One selected facet supplies level-specific positive discriminations. With
#' MML, a different facet may supply category steps. For example,
#' `slope_facet = "Criterion", step_facet = "Rater"` estimates a relative
#' discrimination for each criterion and category-step contrasts for each
#' rater. It does not simultaneously estimate rater discriminations. JML
#' requires the same owner for both blocks. The model has one
#' substantive ability dimension; its structural choices and currently
#' unavailable inferential outputs are separate considerations.
#'
#' Free-slope fits retain numerical estimates for review. [confint.mfrm_fit()]
#' with `parm = "slopes"` supplies approximate relative-slope intervals for eligible MML solutions.
#' MML information-criterion ranking uses separate
#' likelihood and local-solution checks in [compare_mfrm()].
#' [compare_mfrm()] and [build_weighting_review()] also accept `nested = TRUE`
#' for an equal-slope PCM/GPCM test after verifying matching population and
#' constraint settings. Fitted-object scoring and information have their own
#' scope and do not imply a portable GPCM calibration artifact. Consult
#' [gpcm_capability_matrix()] for each operation before using its output.
#'
#' @details
#' Data must be in **long format** (one row per observed rating event).
#' Exact duplicate Person-by-facet combinations are retained, warned once, and
#' propagated as a Data review state by default. Repeated cell keys alone do
#' not establish dependence. With `repeated_ratings = "independent"`, the
#' analyst declares distinct rating events that are independent conditional
#' on Person ability and the model's effects. Keep event IDs in the source
#' data for audit. Event IDs alone do not demonstrate independence; correlated
#' re-ratings require a suitable dependence model. Adding an estimated occasion
#' facet changes the response model and is appropriate only when intended.
#' Repeated ratings do not increase the number of independent Persons.
#'
#' @section Model:
#' `fit_mfrm()` estimates many-facet ordered-response models. The `RSM` and
#' `PCM` branches follow the many-facet Rasch-family tradition (Linacre, 1989);
#' the `GPCM` branch extends the partial-credit kernel with estimated
#' positive slopes under the package's documented identification constraints.
#' For the equal-slope `RSM`/`PCM` branch, a two-facet design
#' (rater \eqn{j}, criterion \eqn{i}) is:
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
#' One response-model family is used per `fit_mfrm()` call. The current public
#' interface does not combine binary, RSM, PCM, or GPCM observations in one
#' fit, define multiple independent rating scales, or accept general
#' threshold/scale anchors and fixed-calibration starting values.
#'
#' With `model = "GPCM"`, the adjacent-category kernel is multiplied by
#' a positive slope for the designated slope-facet level:
#'
#' \deqn{\ln\frac{P(X_{nij} = k)}{P(X_{nij} = k-1)} =
#'   \alpha_g(\eta - \tau_{h,k}),\quad \alpha_g > 0.}
#'
#' Here \eqn{g} indexes the slope facet and \eqn{h} the step facet; MML permits
#' these to differ. Slopes have a sum-to-zero constraint on log slopes, so
#' their geometric mean is 1. Step contrasts are centered within each step
#' owner, separately from facet location effects. Crossing the facets and
#' observing enough categories is necessary to distinguish their effects;
#' allowing separate owners does not ensure identification in a given design.
#' A selected facet owns a vector rather than one common number: if there are
#' \eqn{G} levels, the fit returns \eqn{G} positive slopes with \eqn{G-1} free
#' log-slope contrasts. Every other facet remains additive inside \eqn{\eta}
#' and receives no separate slope. Selecting a rater facet is therefore a
#' different restricted model from selecting a criterion or task facet.
#' The placement of the slope is part of the model identity: it multiplies the
#' complete adjacent-category predictor, including the person coordinate, all
#' additive facet locations and fitted facet interactions inside \eqn{\eta},
#' and the owned step. It is not a loading-only formulation in which the slope
#' multiplies ability while rater severity and other intercept terms remain
#' unscaled. Such a formulation, including TAM multifacet `GPCM.design`
#' constructions with separate linear intercept and slope designs, is a
#' different model unless an algebraic reduction establishes equivalence.
#' In this one-family many-facet GPCM, exactly one facet supplies slopes. MML permits
#' a different step owner; JML requires a shared owner. It is not the broader Uto--Ueno generalized MFRM,
#' whose task and rater slopes enter multiplicatively and whose step owner must
#' be stated separately. The provisional two-family route above fits that product
#' equation with its own scale and output restrictions. Setting every slope to
#' one recovers the equal-discrimination PCM kernel; neither route supplies
#' multidimensional traits or response-style parameters.
#' Under the default one-family `gpcm_mml_identification = "free_population"` branch,
#' the population standard deviation carries the common discrimination scale
#' while the geometric-mean-one slopes describe relative discrimination.
#' Equivalently, on a standardized latent variable the absolute slopes are
#' \eqn{\sigma\alpha_g}. Under
#' `gpcm_mml_identification = "fixed_standard_normal"`, both the population
#' standard deviation and slope geometric mean are fixed to one; that legacy
#' branch is a narrower relative-discrimination model. Under JML, the
#' geometric-mean-one constraint is required to resolve the ability/slope
#' scale because person coordinates are estimated jointly.
#'
#' The model name does not imply finite parameter bounds. The JML branch
#' maximizes the identified joint log-likelihood
#' without a statistical penalty or finite bounds on person, location, step, or
#' slope coordinates. Numerical line-search rejection of non-representable
#' slope proposals is not regularization. When a recession direction is
#' certified, the finite optimizer iterate remains a numerical trace and the
#' primary result uses the appropriate extended-real or typed boundary status;
#' it is not relabelled as a finite maximizer of the original JML objective.
#'
#' With only two ordered categories (\eqn{K = 1}), the `RSM`/`PCM`
#' branch reduces to the usual binary Rasch logit for the single category
#' boundary:
#'
#' \deqn{\ln\frac{P(X_{n\cdot} = 1)}{P(X_{n\cdot} = 0)} = \eta - \tau_1}
#'
#' `GPCM` uses the slope-scaled counterpart
#' \eqn{\alpha_g(\eta - \tau_{g,1})}.
#'
#' With `method = "MML"`, person parameters are integrated out using
#' Gauss-Hermite quadrature and EAP estimates are computed post-hoc.
#' With `method = "JML"`, all parameters are estimated jointly as fixed
#' effects. `"JMLE"` remains an accepted compatibility alias, but package
#' output now uses `"JML"` as the public label. See the "Estimation methods"
#' section of [mfrmr-package] for details.
#'
#' @section Logistic metric (D = 1):
#' Native `RSM`, `PCM` and `GPCM` models, including two-family GMFRM, use the
#' logistic metric (L-metric), with scaling constant \eqn{D = 1}. No additional
#' normal-ogive approximation factor such as 1.7 or 1.702 multiplies the
#' adjacent-category predictor. This response-scale convention is separate
#' from the ability distribution: fixing a population to \eqn{N(0,1)} does
#' not change the link to a normal-ogive model or introduce a 1.7 multiplier.
#'
#' In `RSM`/`PCM`, increasing ability by one unit, holding everything else
#' fixed, increases each adjacent-category log odds by one (one logit).
#' In one-family `GPCM`, the change is \eqn{\alpha_g} logits; in two-family
#' GMFRM it is the product of the two component slopes. Thus a one-unit
#' change in ability is not universally a one-logit change when slopes vary.
#'
#' Before comparing or importing parameter values from other software,
#' reconcile its scaling constant, slope placement and identification
#' constraints with this model. Changing only the multiplier while retaining
#' all parameter values changes response probabilities. An equivalent
#' reparameterization must also transform any affected population parameters,
#' anchors and, for Bayesian estimation, priors. See [export_mfrm_model()]
#' for the same L-metric in the supported Stan/brms exports.
#'
#' @section Weighting policy:
#' `mfrmr` treats `RSM` / `PCM` as the equal-weighting reference route for
#' operational many-facet measurement. In that Rasch-family branch,
#' discrimination is fixed, so the scoring model does not differentially
#' reweight item-facet combinations through estimated slopes.
#'
#' `GPCM` is supported as an alternative when users explicitly accept
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
#' Scores are treated as ordered categories. Although the fitted category
#' probabilities form a multinomial probability vector, category order is part
#' of the likelihood: unordered nominal/multinomial-logit responses are not
#' supported. Poisson, negative-binomial, and grouped binomial-trial counts are
#' also not response families in `fit_mfrm()`. Integer counts supplied as
#' `score` are interpreted only as ordered category codes.
#' Non-numeric score labels are dropped with a warning after coercion, whereas
#' fractional numeric scores outside the rounding tolerance described above
#' are rejected with an error instead of being silently truncated.
#'
#' A positive numeric `weight` may encode a replication/likelihood weight for
#' an ordered-rating row when weighting that conditional contribution is the
#' intended estimand. It is not a general collapsed-person frequency-table
#' interface: under MML, powering responses inside one Person's conditional
#' pattern is not the same as replicating a complete Person pattern after
#' marginalization. A weight also does not turn the score into a count outcome
#' or model dependence among repeated ratings. Non-positive finite weights are
#' excluded during preparation, and non-unit observation-weight fits are not
#' eligible for ordinary inference, facet equivalence, or the common MML
#' information-criterion panel under the current package contract. Normalizing
#' weights to mean one does not remove this restriction. Point estimates and
#' computed curvature/posterior precision remain available for diagnostic
#' review; they do not establish sampling SEs or confidence coverage for the
#' weighted objective. Omitted weights and explicitly all-unit weights use
#' the same eligibility rules. Earlier saved fits require refitting or a
#' current readiness audit; their old inference flags are not carried forward.
#'
#' The fitted many-facet ordered-response model assumes conditional
#' independence of observations given the person and facet parameters
#' (Linacre, 1989). Repeated ratings of the same
#' person-criterion combination by the same rater violate this assumption.
#' When such structures may be present, follow fitting with
#' `diagnose_mfrm(fit, diagnostic_mode = "both")`; its
#' `strict_pairwise_local_dependence` screen is an exploratory check for
#' residual dependence beyond what the additive linear predictor absorbs.
#'
#' Binary responses are therefore supported as ordered two-category scores
#' (for example `0/1` or `1/2`) under the same ordered-response interface.
#' If your observed categories do not start at 0, set `rating_min`/`rating_max`
#' explicitly to avoid unintended recoding assumptions. For example, if the
#' intended instrument is a 1-5 scale but the current sample only uses 2-5,
#' set `rating_min = 1, rating_max = 5` to retain the zero-count category 1
#' in the data-support review. That boundary absence is a review condition for
#' the separate element-boundary contract, not by itself an unsupported free
#' step contrast. By contrast, retaining an unobserved internal category in a
#' polytomous fitted ladder creates an adjacent-step recession direction, so
#' [fit_mfrm()] stops before optimization rather than reporting finite step
#' estimates for that ladder.
#' If these bounds are omitted, the observed score range is used and the
#' provenance is stored in `fit$prep` and `summary(fit)$settings_overview`.
#' Set `options(mfrmr.show_inferred_rating_range = TRUE)` when you want an
#' interactive reminder whenever a bound is inferred.
#' Data-preparation events such as row drops, ID trimming, duplicate
#' person-by-facet cells, and single-level facets are stored in
#' `fit$prep$row_retention` and `fit$prep$preparation_notes`. Routine
#' row-drop/trim/single-level messages are quiet by default; set
#' `options(mfrmr.show_preparation_messages = TRUE)` to show them during
#' interactive checks.
#'
#' When `keep_original = FALSE`, observed gaps such as `1, 3, 5` are recoded
#' internally to a contiguous scale (`1, 2, 3`) and the mapping is stored in
#' `fit$prep$score_map`. To retain zero-count intermediate categories as part
#' of the original scale, set `keep_original = TRUE` in addition to supplying
#' the full `rating_min` / `rating_max` range.
#'
#' @section Facet effects and rating counts:
#' In MML, person ability is integrated out under a `N(0, 1)` distribution (or under
#' the `N(X\beta, \sigma^2)` population model when `population_formula` is
#' supplied). GPCM MML instead activates an intercept-only
#' `N(\beta_0, \sigma^2)` population model by default so its common
#' discrimination scale is estimable. JML instead estimates person parameters;
#' native CML conditions on person totals. Non-person facet parameters
#' (`Rater`, `Criterion`, `Task`, ...) are fitted as fixed effects,
#' identified by the configured centering or anchor constraints. By default,
#' these facets have no hierarchical prior or shrinkage. Optional
#' `facet_shrinkage` applies a separate post-estimation adjustment; see its
#' argument documentation for scope and uncertainty limitations.
#'
#' Few ratings can limit information about an individual facet level; few
#' distinct raters also limit what can be learned about a wider rater population.
#' These are different concerns, and neither is resolved by a count threshold.
#' Precision depends on targeting, connectedness, category use, dependence and
#' the estimator. Jones and Wind (2018) examine repeated ratings in incomplete
#' designs, finding particular sensitivity in examinee estimates. Useful reviews
#' include:
#'
#' - [`facet_small_sample_review()`] for per-level retained rating counts,
#'   count-screening labels and available SEs.
#' - [`detect_facet_nesting()`] and
#'   [`analyze_hierarchical_structure()`] when raters are nested in
#'   regions, schools, or other strata that the additive fixed-effects
#'   MFRM cannot partition out.
#' - [`compute_facet_icc()`] and
#'   [`compute_facet_design_effect()`] for descriptive variance-
#'   component summaries based on `lme4` (optional).
#'
#' `fit$summary$FacetSampleSizeFlag` summarizes the lowest rating-count band
#' across non-person facet levels (`"sparse"` < 10, `"marginal"` < 30,
#' `"standard"` < 50, `"strong"` >= 50). These package screening bands are
#' not guarantees of precision, interval coverage or sample adequacy.
#'
#' @section Jeffreys-penalized JML:
#' `jml_penalty = "jeffreys"` minimizes
#' \eqn{-\ell(\beta)-\tfrac12\log|I(\beta)|} in the actual affine free
#' canonical coordinates, with anchor offsets retained. Here \eqn{I} is the
#' full joint expected response information, including free Person parameters.
#' The optimizer starts at the free-coordinate origin, uses BFGS with
#' `reltol = min(reltol, 1e-14)` and the requested `maxit`, and requires maximum
#' penalized gradient at most `1e-4`, positive information with checked solves,
#' and numerically interior category probabilities. It uses no probability
#' clipping, parameter box, ridge or ordinary extreme-score replacement.
#' If BFGS stops on objective change before satisfying the gradient criterion,
#' at most three local Newton steps use the analytic penalized Jacobian. Steps
#' must be numerically identified, at most one in every free coordinate, not
#' increase the same objective and reduce the full penalized gradient; backtracking is
#' bounded. `checks$polish` records an attempted finish. The original gradient
#' criterion and failure status remain in force if finishing fails.
#'
#' The supported scope is an identified additive RSM/PCM design, unit weights,
#' explicitly declared `rating_min` and `rating_max`, and unchanged full score
#' coding (`category_policy = "preserve"`). Direct finite anchors and ordinary
#' centering conventions are supported; group anchors, interactions, GPCM,
#' population models, shrinkage, corrected-JML arguments and attached ordinary
#' diagnostics are not. `optimizer` must be `"auto"` or `"BFGS"`.
#' Repeated cells require the explicit independent-occasion declaration when
#' substantively justified. Missing categories and constant outcomes are
#' retained on the declared rubric; rank deficiency is still refused.
#'
#' [summary.mfrm_jeffreys_fit()] describes Person, facet and step points,
#' category probabilities and numerical diagnostics. For common-plan Person
#' expected scores use [fair_average_table()] with `reference_plan`, followed
#' by [confint.mfrm_fair_average()] for experimental refit or delta intervals. Calibration,
#' scoring and generating fits must all select this estimator; every replicate
#' re-estimates its calibration and anchors using the same equation. The ordinary
#' delta derivative and inverse information SEs do not describe this estimator;
#' its delta route uses analytic penalized score/anchor derivatives and shared
#' categorical response covariance. It is not an inverse-curvature interval.
#' Parameter covariance, Wald intervals, AIC, likelihood-ratio inference,
#' ordinary FairM/FairZ diagnostics and FACETS export are not supplied.
#'
#' This is an opt-in alternative to ordinary FACETS-oriented JML, separate from
#' both the generator-only Jeffreys option and corrected GPCM. Firth's (1993)
#' fixed-dimensional canonical bias result does not establish unbiased nonlinear
#' expected scores, coverage, or incidental-parameter consistency when the number
#' of Person parameters grows. The binomial finiteness theorem of Kosmidis and
#' Firth (2021) is not asserted as a general many-facet PCM theorem. A qualified
#' local stationary point does not establish a global optimum. Small-sample,
#' sparse and weak-link bias/coverage remain to be qualified independently.
#' References: Firth (1993), \doi{10.1093/biomet/80.1.27}; Kosmidis and Firth
#' (2021), \doi{10.1093/biomet/asaa052}.
#'
#' @section Estimator choice and the JML incidental-parameter caveat:
#' Joint maximum likelihood (`method = "JML"` / `"JMLE"`) estimates
#' both the structural parameters (facets, thresholds, slopes) and
#' every person measure as fixed parameters in one optimization. This
#' is the **incidental-parameter problem** of Neyman & Scott (1948):
#' structural-parameter bias can persist as the number of persons grows
#' with the number of items per person held fixed. In classical Rasch
#' settings this bias can be of order \eqn{1/L} (where \eqn{L} is the number of
#' items per person) and therefore need not vanish by adding persons alone. Wright &
#' Stone (1979) and Wright & Masters (1982, ch. 5) document an
#' empirical \eqn{(L-1)/L} correction that approximately removes the
#' bias for the dichotomous Rasch model; mfrmr does **not** apply
#' that correction (no `bias_correction` argument exists). A common positive
#' multiplier is not a correction for relative GPCM slopes: multiplying all
#' slopes by the same factor leaves their ratios unchanged, and rescaling
#' them to geometric mean one returns the original slopes. A correction for
#' relative slopes therefore needs its own justification.
#'
#' Uncorrected JML facet/location SEs use observation-table information:
#' \deqn{\mathrm{SE}_{g} \approx
#'   \left[\sum_{r \in g} w_r a_r^2
#'   \mathrm{Var}(X_r \mid \widehat\eta_r)\right]^{-1/2}.}
#' Here \eqn{g} is a facet level, \eqn{r} an observed response row,
#' \eqn{w_r} its weight, and \eqn{a_r} its GPCM slope (one for RSM/PCM).
#' Only finite information contributions are used. These exploratory SEs
#' treat the other fitted parameters as fixed; they are not slope SEs or a
#' joint covariance adjusted for estimating Person and structural parameters
#' together. The local joint-curvature check used for portable GPCM JML
#' scoring does not supply such an inferential covariance.
#' Values fixed by anchors or identification constraints are not estimated:
#' diagnostic tables mark them `Fixed = TRUE` and leave their sampling SEs
#' and intervals missing, rather than assigning an observation-information SE.
#'
#' Nuisance-parameter adjustment and estimation bias are separate issues.
#' Haberman (2004, Sections 1.4-1.5) shows for a binary Rasch setting that,
#' with fixed test length, JML can concentrate around a biased limit as
#' the number of persons grows. Even a variance appropriate to that limit
#' does not establish coverage of the true parameter. This result motivates
#' checking test length separately from sample size; it does not validate
#' a correction or an interval for sparse many-facet GPCM. Formal JML
#' slope intervals are not currently available in mfrmr.
#'
#' When comparing software, report response-score adjustment and post-fit
#' bias correction separately, including how the correction defines exposure
#' when responses are missing or unequal across Persons. Matching the label
#' "JML" alone does not establish matching estimates or uncertainty.
#'
#' Practical recommendation:
#'
#' - For manuscript or operational reporting, choose the estimator from the
#'   inferential target and assumptions, and report the choice. MML integrates
#'   person measures under a specified population model and provides marginal
#'   observed-information SEs; consistency of its structural estimates is
#'   conditional on an adequate response model, population distribution, and
#'   regularity conditions.
#' - JML remains useful for a JMLE-oriented FACETS comparison, descriptive or
#'   exploratory work, and designs with substantial information per person.
#'   Report its incidental-parameter limitation and the exploratory basis of
#'   this package's JML structural SEs rather than treating estimator choice as
#'   a universal reporting rule.
#' - Native RSM/PCM CML conditions out Person parameters, requiring no ability
#'   population distribution. It supplies structural-parameter intervals under
#'   conditional-likelihood regularity; see "Conditional maximum likelihood".
#'   General GPCM inference uses MML. JML remains a point-estimation/FACETS
#'   comparison route, not a fixed-exposure truth-centered Wald procedure.
#'   Third-party CML fits can still be imported with [import_erm_fit()].
#'
#' @section MML integration workflow:
#' Choose integration together with the result you need. The default remains
#' `mml_integration = "fixed"`, `quad_points = 31`; this order does not certify
#' accuracy. [mml_quadrature_sensitivity()] compares refits at specified orders.
#' Review changes in estimates and uncertainty before reporting a fixed-grid fit.
#'
#' For ordinary full direct RSM/PCM or one-family GPCM fitting, `"auto"` checks
#' the fitted response patterns and refines integration when needed. Read
#' `fit$config$mml_accuracy` and the resolved rule in
#' `fit$config$estimation_control$mml_integration`. Downstream support depends
#' on this resolved rule, which can be adaptive even when the initial rule was
#' fixed. Check the intended output's model and integration scope before fitting:
#' - Ordinary fixed-N(0,1) RSM/PCM model and sandwich facet intervals require
#'   fixed integration, as do joint cell intervals requested with `interaction`
#'   and independent anchor-source propagation. See [mfrm_facet_intervals()].
#' - Native-location intervals for estimated-normal RSM/PCM and one-family
#'   GPCM accept fixed or adaptive direct MML within their stated scope.
#'   RSM/PCM population-coefficient and residual-variance intervals use the same
#'   qualification; see [mfrm_population_intervals()]. GPCM latent-regression
#'   location intervals are unavailable.
#' - Automatic integration does not support EM/hybrid engines, checkpoints,
#'   comparison-only fitting or two slope families. ConQuest overlap export
#'   currently requires a resolved fixed rule.
#'
#' Overall fit readiness, interval eligibility and coverage are different
#' statements. A supported target-specific interval can have `CIEligible = TRUE`
#' while its unchanged source fit has `InferenceReady = FALSE`. Inspect the
#' interval's `Status`, `checks`, `numerical_checks` and cautions. Eligibility
#' means the procedure's checks passed; it does not establish global
#' identification, nominal finite-sample coverage or availability of other
#' inference procedures. Automatic likelihood/gradient/posterior checks also
#' do not replace the interval's full-information and finer-order review.
#'
#' @section MML inference conditions:
#' There is no ten-response-row minimum. Full RSM/PCM JML, native CML and
#' fixed-population MML can accept one Person when their model-specific checks
#' pass. One-Person JML/MML also requires `lpSolve` for a pre-fit additive
#' recession check; a detected or unresolved boundary returns no finite estimate.
#' Separate rating events must satisfy the assumed conditional independence;
#' explicitly declared repeated ratings are not evidence that independence holds.
#' The returned preparation note, interval cautions and summary explain the
#' single-Person limitation. Any available model-information interval remains
#' an unqualified normal approximation; one Person supplies neither independent-
#' Person sandwich uncertainty nor empirical between-Person variation.
#' Estimated-population MML, GPCM, corrected JML and comparison-only fitting
#' retain their at-least-two-Person software scope. Data and anchor review can
#' inspect one Person. The calibration interface still requires two observed
#' score categories; these input checks do not establish identification or
#' precision. Scoring one Person against saved calibration is a separate task;
#' see [score_mfrm_calibration()]. Small supported datasets can
#' return point estimates and asymptotic uncertainty when the model-specific
#' checks pass. A singleton category count alone need not prevent RSM/PCM
#' model-information intervals: [mfrm_facet_intervals()] performs a separate
#' joint-information and quadrature review and retains the small-sample caution.
#'
#' MML integrates the joint response vector for each Person under the declared
#' ability distribution. Model-based structural SEs use the inverse full
#' marginal observed information, including estimated population and other
#' nuisance coordinates. Holding these coordinates fixed gives a different
#' uncertainty calculation. EAP Person scores are predictions under the fitted
#' model, not additional free parameters in this information matrix.
#'
#' Regular likelihood inference requires the response and population models
#' to be correct for the actual assignments and covariates, independent Persons,
#' a fixed-dimensional identifiable model, a consistent interior solution,
#' smooth derivatives and growing information in every free direction.
#' This permits fixed exposure per Person under suitable designs; it does not
#' establish consistency merely by increasing Person count. For an estimated
#' normal population, zero variance is a boundary outside this regular argument.
#' Numerical convergence or a positive local Hessian does not prove global
#' identification, exclude a boundary solution or verify the population model.
#'
#' Finite quadrature approximates the continuous-normal likelihood. Optimizer
#' and integration errors must be small relative to sampling uncertainty;
#' a fixed quadrature order does not automatically maintain this property as
#' sample size grows. Comparing quadrature orders diagnoses sensitivity but
#' does not provide a continuous-integral error bound. See
#' [mfrm_facet_intervals()] for the assumptions, available checks and the joint
#' delta method for contrasts divided by the same fit's estimated population SD.
#'
#' @section All-minimum and all-maximum Persons:
#' `fit_mfrm()` has no public option to remove Persons before fitting because
#' all their responses are at the minimum or maximum. Both JML and MML retain
#' those observed responses. Extreme-response flags use the usable rows after
#' data preparation; missing responses do not count as intermediate scores.
#' Declare `rating_min` and `rating_max` from the rubric so that an observed
#' maximum is not mistaken for the intended scale maximum. Inspect
#' `fit$facets$person$Extreme` for `"low"`, `"high"`, or `"none"`.
#'
#' **Ordinary JML.** `fit_mfrm()` does not replace extreme response scores
#' before JML fitting. For an independently free Person with all-minimum or
#' all-maximum responses, the primary `Estimate` is `-Inf` or `Inf`;
#' `OptimizerEstimate` retains the finite computational trace, not a finite
#' Person MLE. Fixed Person anchors retain their supplied values, and coupled
#' constraints require their own boundary review. The Person audit's
#' `BoundaryState = "has_exclusions"` identifies nonfinite parameters; it
#' does not mean that these Persons or their input rows were deleted.
#' Finite structural estimates alone do not establish a finite joint maximum.
#'
#' **Corrected JML.** In the supported shared-owner GPCM route, extreme
#' Persons also remain in the data and assignment-pattern accounting. Their
#' profiled abilities have infinite limits and their structural estimating-
#' equation contributions are zero at those limits. This boundary calculation
#' is not an input filter or a general proof of bias removal.
#'
#' **MML.** Extreme Persons contribute to the marginal likelihood through
#' integration over the specified or estimated population distribution. Their
#' reported abilities are posterior EAPs, with posterior SDs rather than
#' frequentist Person-MLE standard errors. A proper normal population model
#' with finite parameters gives finite EAPs, subject to valid numerical
#' integration; an extreme response pattern alone does not require deletion.
#' This does not guarantee convergence, identification or valid structural
#' intervals for the fitted model. Removing these Persons would change the
#' observed sample and marginal likelihood.
#'
#' New-Person scoring with [predict_mfrm_units()] is a separate operation.
#' With an eligible calibration and an admitted scoring prior, it can return
#' finite EAPs for extreme new Persons, including after JML calibration. Those
#' scores do not replace the original JML Person MLEs. Calibration-source
#' checks and scoring-batch integration checks must both pass; increasing
#' scoring nodes does not resolve a boundary in the source calibration.
#' A finite display from `fair_average_table(..., xtreme = ...)`, or placement
#' at the end of a Wright map, likewise does not change the fitted model or
#' correct JML bias.
#'
#' @section Model-estimated facet interactions:
#' `facet_interactions` adds confirmatory fixed-effect interaction terms to the
#' linear predictor. For example, `facet_interactions = "Rater:Criterion"`
#' estimates a rater-by-criterion deviation matrix in the same likelihood as
#' the main MFRM fit. The additive reference is
#'
#' \deqn{\eta_{nij} = \theta_n - \delta_j - \beta_i}
#'
#' and the interaction extension is
#'
#' \deqn{\eta_{nij} = \theta_n - \delta_j - \beta_i + \gamma_{ji}}
#'
#' where the interaction block is identified by zero marginal sums:
#'
#' \deqn{\sum_j \gamma_{ji} = 0,\quad \sum_i \gamma_{ji} = 0.}
#'
#' With \eqn{J} levels of the first facet and \eqn{I} levels of the second
#' facet, this contributes \eqn{(J - 1)(I - 1)} free parameters. Positive
#' interaction estimates indicate scores higher than expected under the
#' additive main-effects model for that facet-level combination; negative
#' estimates indicate lower-than-expected scores.
#'
#' This is a model-estimated interaction term, not the residual screening
#' reported by [estimate_bias()] or [estimate_all_bias()]. In line with the
#' MFRM bias-interaction literature, the facet pair should be named explicitly
#' before fitting. Exploratory use is possible, but should be reported as
#' screening, with sparse-cell and multiplicity caveats. The current
#' implementation is intentionally narrow: two-way non-person facet
#' interactions for `RSM` and `PCM` only, estimated as fixed effects. GPCM
#' interactions, person interactions, higher-order interactions, and
#' random-effect facet interactions are deferred.

#' This is ordered binary support, not a separate nominal-response model.
#' In `PCM`, a binary fit still uses one threshold per `step_facet` level on
#' the shared observed-score scale.
#'
#' Supported model/estimation combinations:
#' - `model = "RSM"` with `method = "MML"`, `"JML"`/`"JMLE"`, or `"CML"`
#' - `model = "PCM"` with a designated `step_facet`; CML requires it explicitly
#' - Native CML returns a separate structural-parameter result with
#'   conditional-information intervals; ordinary fit/scoring/diagnostic methods
#'   below apply to MML/JML results
#' - `facet_interactions` with `model = "RSM"` or `"PCM"` for explicit
#'   two-way non-person facet interactions
#' - `model = "GPCM"` uses one positive slope family. MML permits separate
#'   slope/step owners; JML requires `slope_facet == step_facet`. Core
#'   summaries, fitted-object posterior scoring, [compute_information()],
#'   Wright/pathway/CCC fit plots, [diagnose_mfrm()], residual-PCA follow-up,
#'   [interrater_agreement_table()], [unexpected_response_table()],
#'   [displacement_table()], [measurable_summary_table()],
#'   [rating_scale_table()], [facet_quality_dashboard()],
#'   [reporting_checklist()], [category_structure_report()],
#'   [category_curves_report()], and graph/scorefile
#'   [facets_output_file_bundle()] routes are available with score-side
#'   caveats. Separate-owner MML supports slope and curve intervals,
#'   same-design [bootstrap_mfrm_gpcm()], and matched PCM comparison.
#'   [build_weighting_review()] and the following simulation/design workflows
#'   still require the same slope/step owner. Direct simulation specifications
#'   and data generation are also supported through
#'   [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()], and
#'   [simulate_mfrm_data()] when the slope-aware generator contract is stored
#'   explicitly; direct recovery checks are available through
#'   [simulate_mfrm_recovery()] and [assess_mfrm_recovery()]. Slope-aware
#'   [fair_average_table()] and [estimate_bias()] are available with their
#'   documented caveats. Role-based design evaluation, population forecasting,
#'   diagnostic-screening, and signal-detection helpers are available as
#'   caveated sensitivity evidence. Full FACETS-style score-side contract
#'   review, posterior predictive checks, and MCMC estimation are not available
#'   for `GPCM`. Use
#'   [gpcm_capability_matrix()] as the formal boundary statement for the
#'   current `GPCM` scope.
#'
#' Latent-regression status:
#' - `population_formula = NULL` keeps the standard unconditional behavior for
#'   RSM/PCM and JML. For GPCM MML, the default
#'   `gpcm_mml_identification = "free_population"` constructs an intercept-only
#'   population model internally; use `"fixed_standard_normal"` only to
#'   reproduce the legacy restricted likelihood.
#' - Supplying `population_formula` activates latent regression for
#'   `method = "MML"` only.
#' - This implementation assumes a one-dimensional conditional-normal population
#'   model with person-specific quadrature nodes
#'   \eqn{\theta_{nq} = x_n^\top \beta + \sigma z_q}.
#' - Background variables must be supplied in `person_data`; numeric/logical
#'   columns and categorical factor/character columns are expanded through
#'   `stats::model.matrix()`.
#' - Documented overlap with the ConQuest latent-regression model is
#'   limited to direct estimation from response data under a unidimensional
#'   `MML` population model with package-built model-matrix covariates. It
#'   should not be described as numerical equivalence for arbitrary imported design matrices,
#'   multidimensional models, or the full ConQuest plausible-values workflow.
#' - `predict_mfrm_units()` and `sample_mfrm_plausible_values()` can score
#'   latent-regression fits under the fitted population model, but they require
#'   one-row-per-person background data for scored units when the fitted
#'   population model includes covariates. Intercept-only latent-regression
#'   fits (`population_formula = ~ 1`) can reconstruct that minimal person
#'   table internally during scoring.
#'
#' @section Latent-regression workflow:
#' For an initial latent-regression run, keep the setup explicit:
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
#' With covariates, the estimated `sigma2` is the residual variance of ability
#' conditional on those covariates, not the marginal population variance.
#' Marginal variance also depends on the distribution of the covariates.
#' Standardization and polynomial/spline bases learned during fitting are
#' reused for new persons; changing the scoring cohort does not redefine them.
#' For custom transformations, supply a prediction-aware R transformation or
#' precompute predictors with fixed training constants. Older saved fits can
#' reconstruct transformed terms only from retained person data that reproduce
#' the training design; otherwise scoring requests a new calibration fit.
#'
#' For an intercept-only model, `population_formula = ~ 1` estimates a single
#' population mean and variance by default. In supported RSM/PCM direct MML,
#' `population_sd = 1` instead fixes the SD and estimates only the mean.
#' Training still requires a one-row-per-person
#' ID table in `person_data`. Adjusting these two moments retains a normal
#' population shape; it does not learn skewness or establish that the training
#' population represents the people who will receive scores. New-Person
#' intervals condition on the estimated calibration and population parameters
#' and exclude uncertainty from estimating them. Current population-model
#' scoring remains a review workflow; see [predict_mfrm_units()].
#'
#' @section Latent-regression standard-error caveat:
#' `summary(fit)$population_coefficients` reports point estimates of
#' \eqn{\hat{\boldsymbol{\beta}}}; `population_overview` reports the estimated
#' population variance. These tables do **not** provide standard errors,
#' confidence intervals, or asymptotic z / Wald statistics for the population
#' parameters. For an estimated normal RSM/PCM population,
#' [mfrm_population_intervals()] separately requests experimental pointwise
#' coefficient, prespecified coefficient-contrast and log-variance intervals
#' after full joint numerical checks.
#' Its result has a `vcov()` method on coefficient/log-variance coordinates;
#' the ordinary fit summary does not acquire automatic intervals.
#' The internal MML observed-information calculation includes
#' \eqn{(\boldsymbol{\beta}, \log\sigma^2)} when computing joint covariance
#' for structural-parameter SEs. With explicit `population_sd`, variance is
#' fixed and its log-variance coordinate is excluded from that covariance;
#' coefficient and structural uncertainty are conditional on the supplied SD.
#' Population intervals for GPCM remain
#' unavailable through this new helper. Sampling coverage is not established
#' by numerical checks. Treat the population summary tables as descriptive
#' point estimates; use the separately saved interval result when reporting
#' the experimental approximation and its limitations.
#'
#' The default latent-regression identification uses
#' `noncenter_facet = "Person"`. Without anchors, this sum-to-zero-centers
#' all non-Person facets. With anchors, `centering_scope = "facet"` includes
#' their contributions in each facet's zero sum; `"remainder"` retains the
#' earlier restriction on the unconstrained remainder. The population model
#' is never changed implicitly to accommodate anchors. The restricted
#' anchor-origin alternative below uses an explicitly fixed population SD.
#'
#' @section Anchor origin and cohort means:
#' When an old rater panel is anchored, leaving Rater uncentered allows new
#' raters to have a nonzero mean severity. Fixing the Person mean to zero can
#' nevertheless distort estimates when the receiving cohort differs. For
#' RSM/PCM direct MML, `population_formula = ~ 1` with `population_sd = 1`
#' estimates the receiving mean on an explicitly fixed-SD normal population.
#' Set `noncenter_facet = "Rater"` and supply a retained direct Rater anchor.
#' An anchor on another facet does not resolve the mean/Rater alias.
#' The complete constrained design must still have full rank; one anchor is
#' sufficient only in an otherwise identified design. Group anchors, dummy
#' facets, population covariates and interactions are outside this alternative.
#'
#' For example, begin with the bundled ratings and an externally supplied
#' anchor value (0.3 below is illustrative, not a validated calibration):
#' ```r
#' ratings <- load_mfrmr_data("example_core")
#' fixed_rater <- data.frame(Facet = "Rater", Level = "R01", Anchor = 0.3)
#' cohort_fit <- fit_mfrm(
#'   ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#'   model = "RSM", method = "MML", anchors = fixed_rater,
#'   noncenter_facet = "Rater", population_formula = ~ 1,
#'   person_data = unique(ratings["Person"]), population_sd = 1
#' )
#' summary(cohort_fit)$population_coefficients # Estimated receiving-cohort mean
#' summary(cohort_fit)$population_overview     # Fixed versus estimated variance
#' ```
#' Use [mfrm_population_intervals()] for the mean and [mfrm_facet_intervals()]
#' for rater uncertainty after their independent numerical checks. The fixed
#' variance remains visible without an inferential CI. Covariance includes
#' all free calibration coordinates and the population mean, but excludes
#' uncertainty in the supplied anchors and SD. No general coverage guarantee
#' follows from agreement with another program.
#'
#' Other facets retain their declared centering. Rater anchors do not establish
#' item invariance: a changed item pool or its difficulty distribution can
#' change the mean's substantive interpretation. Report the anchored panel,
#' item/criterion reference, fixed SD, and assumed stability across cohorts.
#'
#' @section Anchors and scale constraints:
#' Anchor inputs are optional and have distinct roles:
#' - `anchors` contains facet/level/fixed-value information and imposes direct
#'   equality constraints on selected parameters.
#' - `group_anchors` contains facet/level/group/group-value information and
#'   constrains each declared group mean. Its interpretation is conditional on
#'   the externally justified target or equal-mean assumption.
#' - Common Persons, raters, items, or rating events are properties of the
#'   observed design. Neither constraint type creates empirical overlap or
#'   proves that disconnected subsets are substantively comparable.
#' Both are normalized internally, so column names can be flexible
#' (`facet`, `level`, `anchor`, `group`, `groupvalue`, etc.).
#'
#' Direct anchors fix facet locations, not thresholds or GPCM slopes. Values
#' must use the source model's internal metric and orientation, not transformed
#' display units. [make_anchor_table()] exports this metric. For reuse of a
#' complete fixed calibration, see [extract_mfrm_calibration()] and [score_mfrm_calibration()].
#'
#' With the default `noncenter_facet = "Person"`, a partial Rater anchor such
#' as R01 = 0.7 makes the other Rater levels sum to -0.7 under the default
#' `centering_scope = "facet"`. With one remaining level, it is fixed at -0.7.
#' Group means contribute their target multiplied by the group's size, including
#' directly fixed members. If anchors/groups cover all levels, their total must
#' agree with the facet's zero sum. Contradictions stop review and fitting even
#' with `anchor_policy = "silent"`.
#'
#' Use `centering_scope = "remainder"` to reproduce earlier partial anchoring:
#' the ungrouped, unanchored levels sum to zero and a singleton is fixed at zero.
#' Saved results retain their original constraints; reloading does not migrate
#' estimates. Changing conventions requires refitting.
#'
#' For ordinary additive RSM/PCM JML, `centering_scope = "anchor_origin"`
#' implements the FACETS convention that a directly or group-anchored facet is not
#' additionally centered. With a Person anchor, choose an unanchored facet such
#' as `noncenter_facet = "Rater"`; other unanchored facets remain centered.
#' A Rater anchor can instead retain the default noncenter Person facet.
#' The noncenter facet must be active and have neither direct nor group anchors.
#' A single direct anchor or group mean supplies
#' its facet's origin; multiple fixed values can impose additional substantive
#' constraints. Without Person group constraints, if all retained responses of
#' the directly anchored Persons are
#' at the same declared minimum or maximum, fitting stops before optimization,
#' without requiring `lpSolve`: free Person and noncenter-facet locations can
#' diverge together while preserving every fixed value. FACETS 4.5.0 can report
#' finite values for such inputs, but these are not an ordinary finite JML
#' solution and are not reproduced by this route. An extreme anchor does not
#' alone trigger this rule when another Person anchor supplies interior or
#' opposite-boundary responses. Other detected joint additive recession
#' directions affecting calibration also stop this route; passing the checks
#' does not certify every boundary configuration. General extreme-score
#' exclusion and every FACETS option are not equivalent. See the FACETS
#' [centering specification](https://www.winsteps.com/facetman64/centerfacet.htm).
#' Match category support, signs, the response design and population assumptions
#' for TAM/ConQuest comparisons; an anchor value alone does not establish
#' equivalence. All these routes condition on supplied anchor values and do not
#' propagate error from an external estimated calibration.
#'
#' If an old zero-sum panel is fully anchored and new raters are added,
#' whole-facet centering also forces the new raters to sum to zero. This is
#' an additional constraint, not evidence that the new panel has the same
#' average severity. In supported fixed-population MML models, choosing
#' `noncenter_facet = "Rater"` can leave that facet uncentered. This is not
#' a general remedy: only one facet can be noncentered. The RSM/PCM
#' intercept-only, fixed-SD alternative is described above. Other latent
#' regression, including the usual one-family GPCM population model, requires
#' `noncenter_facet = "Person"`. Review the complete identification and
#' population assumptions before using anchors to compare cohorts.
#'
#' A group target is an unweighted arithmetic mean of its declared levels,
#' including any directly anchored members. It is not weighted by rating
#' counts. FACETS `G` supplies per-element contributions to a group total;
#' convert their arithmetic mean to one repeated `GroupValue`. For example,
#' contributions 0.2 and 1.0 become `GroupValue = 0.6` for both members.
#' Unequal `GroupValue` entries are conflicting targets, not contributions;
#' use `anchor_policy = "error"` to stop on review issues. FACETS group 0
#' denotes direct anchors and must be translated to `anchors`; mfrmr treats
#' `Group = "0"` as an ordinary group identifier. The FACETS `X` convention
#' that removes extreme members is not applied automatically. See
#' [FACETS element specifications](https://www.winsteps.com/facetman64/elements.htm).
#'
#' Group-origin support retains constrained maximum likelihood, not every
#' native FACETS update. Matched interior RSM/PCM examples agree when a single
#' group only sets the origin. With additional restrictions, FACETS 4.5.0 can
#' satisfy the group means while leaving a nonzero likelihood derivative along
#' a constraint-preserving direction. mfrmr optimizes that likelihood rather
#' than reproducing those numerical differences. Matching means alone does
#' not establish parameter, SE or interval equivalence; group-anchor extreme
#' extrapolation remains unsupported. In MML, individual abilities are
#' integrated out: Person direct and
#' group anchors are unsupported and rejected even with `anchor_policy =
#' "silent"`. Sharing Persons across raters requires no fixed Person values.
#'
#' Ordinary fit summaries condition on supplied anchor values. Independent
#' source-anchor covariance can be propagated explicitly in the bounded
#' fixed-SD RSM/PCM route through [mfrm_facet_intervals()] and
#' [mfrm_population_intervals()]; it is not added automatically.
#' By default, downstream SEs and intervals condition on those values and retain the
#' selected estimator's other limitations. A missing SE for a fixed value
#' does not mean its source calibration was error-free. See
#' [mfrmr_linking_and_dff] for the theoretical background and a loading-first
#' example, and [review_mfrm_anchors()] for receiving-data checks.
#'
#' Anchor review behavior:
#' - `fit_mfrm()` automatically runs an anchor review.
#' - invalid rows are removed before estimation.
#' - duplicate rows keep the last occurrence for each key.
#' - `anchor_policy` controls whether detected issues are warned, treated as
#'   errors, or kept silent.
#' - the review checks table/data compatibility and local support counts; it
#'   does not establish source-fit readiness, cross-run identity, parameter
#'   invariance, or the validity of a group-mean assumption.
#'
#' Facet sign orientation:
#' - facets listed in `positive_facets` are treated as `+1`
#' - all other facets are treated as `-1`
#' This affects interpretation of reported facet measures.
#'
#' @section Estimator-specific checks before fitting:
#' Before optimization, mfrmr builds a sparse adjacent-category-logit design
#' in the same constrained free coordinates used by the optimizer. The check
#' includes Person coordinates for JML, integrates them out for MML, and
#' includes facet anchors, group constraints, signs, supported two-way
#' interactions, and RSM/PCM step coordinates.
#'
#' An exactly rank-deficient design stops with a structured
#' `mfrmr_estimability_error`; its `estimability` field records rank, nullity,
#' parameter blocks, tolerance checks, and a bounded null-direction
#' explanation. Optimization is not run. A full-rank MML fixed-effect design
#' whose corresponding free-Person JML design is rank deficient returns a
#' fit with an `mfrmr_estimability_warning`: its cross-panel contrasts rely on
#' the common latent-population assumption and remain review-only.
#'
#' Inspect `fit$data_review$estimability`. RSM and PCM use the full linear
#' free-coordinate check. For GPCM and an active latent-regression
#' residual variance, the additive block is audited before fitting. A retained
#' vector also records the analytic free-to-expanded log/natural-scale
#' transformation Jacobians and a central-difference check in
#' `fit$data_review$estimability$nonlinear_transformation`. This verifies the
#' parameterization only; it is not a response-likelihood Jacobian or a
#' structural-identification result. A stationary retained solution of modest
#' free dimension also receives a local observed-information Hessian and a
#' recorded eigenvalue-tolerance ladder in
#' `fit$data_review$estimability$fitted_information`. Nonstationary or larger
#' fits retain an explicit not-evaluated status. This fitted-information layer
#' is diagnostic only: it does not yet classify weak information, make the
#' nonlinear check complete, or turn full additive rank into a full-model
#' estimability claim. Eligible nonlinear MML fits also receive bounded
#' observed-pattern and all-response-pattern score checks. The latter operates
#' on each Person's retained observation design under unit row weights and
#' records probability-normalization, zero-expected-score, expected-information,
#' and selected numerical-derivative summaries. Missing rows are not imputed;
#' nonunit weights and excessive pattern grids retain explicit not-evaluated
#' states. Mathematically identical Person observation designs are evaluated
#' once and reconstructed by exact multiplicity; active latent-regression
#' covariate rows are part of this identity. Only conceptual and evaluated
#' workload summaries are retained. These retained-point diagnostics do not
#' by themselves establish global structural identification, weak-information
#' status, or readiness.
#'
#' `fit$data_review$estimability$nonlinear_local_estimability` interprets only
#' the first-order local rank that these maps support. For JML GPCM, full
#' column rank of the complete conditional adjacent-logit Jacobian is a
#' sufficient retained-point local certificate. For fixed-quadrature MML with
#' unit row weights and finite parameters, the Person-specific observed-
#' pattern score vectors are part of the positive finite response-pattern
#' support. If those vectors span every optimizer free coordinate, the full
#' expected score information is positive definite; exhaustive enumeration is
#' unnecessary for this sufficient direction. A rank-deficient observed subset
#' is inconclusive and is classified only when the all-pattern enumeration is
#' available. The record explicitly leaves continuous-integral and global
#' identification, boundary status, weak information, and inference readiness
#' unclassified.
#'
#' @section Choosing maxit without result-driven tuning:
#' Treat `maxit` as a predeclared computational budget, not as a value to tune
#' until preferred estimates appear.
#'
#' 1. Choose the model, estimation method, optimizer, tolerance, quadrature
#'    rule, and initial `maxit` before examining coefficient or fit results.
#'    The default `maxit = 400` is the package starting point for an analysis.
#' 2. Inspect `FitReadiness` and `InferenceReady`, together with the eligibility
#'    and reasons for the specific calculation. RSM/PCM model intervals can
#'    pass their own singleton-category and information checks while the fit
#'    retains its category review; see [mfrm_facet_intervals()]. This exception
#'    does not waive numerical convergence. Also inspect purpose-specific Design,
#'    Stability, Diagnostics, and Reporting workflow rows. Optimizer code zero
#'    alone is insufficient. A missing, nonfinite or incomplete terminal gradient
#'    produces `ConvergenceReason = "code_zero_unavailable_gradient"` and review
#'    status, even when the optimizer reports code zero. Resolve the numerical
#'    evaluation before using the fit for inference; more iterations alone do
#'    not establish a valid gradient.
#' 3. If `ConvergenceStatus == "iteration_limit"`, keep that fit review-only.
#'    Refit the same data, model, method, anchors, optimizer, tolerance, and
#'    quadrature rule with the next ceiling in a prespecified sequence, such as
#'    400, 800, then 1600. Do not choose among runs by coefficient size,
#'    statistical significance, fit statistics, or agreement with an expected
#'    answer.
#' 4. The first run in that sequence that satisfies the numerical-readiness
#'    criteria becomes
#'    eligible for interpretation. If separately ready runs differ materially,
#'    treat the difference as numerical instability and review the model,
#'    identification, data support, and optimizer rather than selecting the
#'    preferred result.
#' 5. Report the requested `maxit`, actual iteration/evaluation counts,
#'    convergence status and reason, optimizer, terminal gradient, and any
#'    polishing stages. These are retained in `fit$summary` and
#'    `fit$opt$optimizer_polish`.
#'
#' This rule applies to both JML and MML. JML can require a larger computation
#' budget because it estimates one fixed effect per person; increasing `maxit`
#' does not make JML equivalent to MML and must not be used to switch the
#' estimand after seeing results.
#'
#' @section Performance tips:
#' When JML is the prespecified estimand, it is often faster than MML but may
#' require a larger `maxit` on larger datasets. Do not switch from MML to JML
#' only to shorten runtime: integrating over a population distribution and
#' treating person parameters as fixed effects are different analysis choices.
#'
#' For MML runs, `quad_points` is the main accuracy/speed trade-off.
#' The `@param quad_points` tier table is the authoritative reference;
#' in short:
#' - `quad_points = 7` is a lightweight screening setting; do not use its IC
#'   values for automatic model selection.
#' - `quad_points = 15` is an intermediate review option when runtime matters;
#'   automatic IC ranking remains disabled.
#' - `quad_points = 31` is the package default and a starting point, not a
#'   guarantee that integration error is negligible.
#' - `quad_points = 61` (or higher) supplies candidate denser grids at
#'   additional computational cost; no fixed order is sufficient for every
#'   response pattern.
#' - Use [mml_quadrature_sensitivity()] to inspect same-data movement without
#'   an automatic stable/unstable decision.
#' - `mml_engine = "direct"` remains the most stable general-purpose path.
#' - `mml_engine = "em"` or `"hybrid"` currently target `RSM` / `PCM` fits
#'   without a latent-regression population model.
#' - Benchmark your own workload before using `mml_engine = "em"` or
#'   `"hybrid"` for final reporting; `direct` remains the documented default
#'   when you have not compared engines for your data.
#' - When a direct code-zero stage stops ahead of the terminal-gradient check,
#'   bounded polishing is automatic. Inspect `fit$opt$optimizer_polish$Stages`
#'   rather than repeatedly lowering `reltol` without reviewing the retained
#'   objective, gradient, and parameter changes.
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
#' For `GPCM`, residual-based mean-square fit screens are also
#' best treated as exploratory diagnostics rather than strict Rasch-style
#' invariance tests, because the discrimination parameter is free.
#'
#' @section Interpreting output:
#' Start with `summary(fit)` and the estimator selected above. For ordinary
#' MML/JML results, read its `$decision` before interpreting the Person and
#' facet overviews. For CML, read `$locations` and `$steps`, or `$comparisons`
#' for a comparison-only fit; those results do not contain Person abilities.
#' A reported estimate, numerical convergence and an eligible confidence
#' interval are different claims. Read the associated checks and limitations.
#'
#' @section Typical workflow:
#' 1. Fit the prespecified model and inspect `summary(fit)`.
#' 2. For ordinary MML/JML RSM/PCM, collect `res <- mfrm_results(fit)` and
#'    inspect `summary(res)$triage` and `$next_actions`. Pass already computed
#'    diagnostics through `diagnostics = diagnostics` to reuse them.
#' 3. Choose `plot(res, ...)`, [mfrm_report()], [export_mfrm_results()] or
#'    `saveRDS(res, file)` according to the output you need. See
#'    [mfrmr_reports_and_tables] for the complete example.
#' 4. For CML, use [mfrm_facet_intervals()] and its own summary, plot,
#'    [apa_table()] and RDS saving. For GPCM or extended models, use
#'    [gpcm_capability_matrix()] and [mfrm_results()] to select supported
#'    outputs. Shared function names do not make ordinary diagnostics or
#'    scoring available for every estimator or result class.
#'
#' @section Information-criterion contract:
#' For an eligible fixed-facet `MML` fit, let `D = -2 * LogLik`, let `k` be
#' `Npar` (the retained free optimization-vector dimension after constraints),
#' and let `N_person` be the number of unique prepared Persons. The canonical
#' panel is `AIC = D + 2 * k`, `BIC = D + log(N_person) * k`, and
#' `SABIC = D + log((N_person + 2) / 24) * k`.
#'
#' `ResponseRows`, `WeightedResponseTotal`, `Persons`, and `ICSampleSize` are
#' separate fields. The compatibility field `N` retains its earlier
#' response-row or summed-observation-weight meaning and is not the BIC sample
#' size. Explicit all-unit weights remain eligible; every non-unit
#' observation-weight fit, JML fit, and object without the current contract
#' identity is excluded from the common MML panel. Its canonical
#' `AIC`/`BIC`/`SABIC` fields are `NA`, while any retained raw values are
#' explicitly named `LegacyAIC` and `LegacyBIC`. At 22 or fewer Persons,
#' SABIC is displayed only as sensitivity evidence and
#' `SABICSelectable = FALSE`.
#'
#' Integration adequacy is recorded separately from formula eligibility.
#' `ICIntegrationTier` is `"coarse_screening"` below 15 points,
#' `"intermediate_review"` at 15--30, `"standard_start"` at 31--60, and
#' `"dense_sensitivity"` at 61 or more. Raw canonical criteria remain visible
#' in every eligible MML tier, but `ICSelectable = FALSE` below 31 points;
#' automatic deltas, criterion weights, preferences, and LRT are suppressed.
#' A close or consequential q>=31 comparison should still be reevaluated on a
#' denser common grid.
#'
#' @section References:
#' The ordered-category many-facet formulation follows Linacre (1989), with
#' the `RSM` and `PCM` branches grounded in Andrich (1978) and Masters (1982).
#' The `GPCM` branch follows the generalized partial credit
#' formulation of Muraki (1992) under a package-specific positive
#' log-slope identification convention. The `MML` route follows the
#' quadrature-based marginal-likelihood framework of Bock and Aitkin (1981).
#'
#' - Akaike, H. (1974). *A new look at the statistical model identification*.
#'   IEEE Transactions on Automatic Control, 19(6), 716-723.
#' - Andrich, D. (1978). *A rating formulation for ordered response
#'   categories*. Psychometrika, 43(4), 561-573.
#' - Bock, R. D., & Aitkin, M. (1981). *Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm*. Psychometrika, 46(4),
#'   443-459.
#' - Dhaene, G., & Weidner, M. (2023). *Approximate functional differencing*.
#'   SERIEs, 14, 379-416. \doi{10.1007/s13209-023-00283-1}.
#' - Haberman, S. J. (2004). *Joint and conditional maximum likelihood
#'   estimation for the Rasch model for binary responses*. ETS Research
#'   Report RR-04-20. \doi{10.1002/j.2333-8504.2004.tb01947.x}.
#' - Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#' - Masters, G. N. (1982). *A Rasch model for partial credit scoring*.
#'   Psychometrika, 47(2), 149-174.
#' - Myford, C. M., & Wolfe, E. W. (2003). Detecting and measuring rater
#'   effects using many-facet Rasch measurement: Part I. *Journal of Applied
#'   Measurement*, 4(4), 386-422.
#' - Myford, C. M., & Wolfe, E. W. (2004). Detecting and measuring rater
#'   effects using many-facet Rasch measurement: Part II. *Journal of Applied
#'   Measurement*, 5(2), 189-227.
#' - Muraki, E. (1992). *A generalized partial credit model: Application of an
#'   EM algorithm*. Applied Psychological Measurement, 16(2), 159-176.
#' - Uto, M., & Ueno, M. (2020). *A generalized many-facet Rasch model and its
#'   Bayesian estimation using Hamiltonian Monte Carlo*. Behaviormetrika, 47,
#'   469-496.
#' - Robitzsch, A., & Steinfeld, J. (2018). *Item response models for human
#'   ratings: Overview, estimation methods, and implementation in R*.
#'   Psychological Test and Assessment Modeling, 60(1), 101-139.
#' - Schwarz, G. (1978). *Estimating the dimension of a model*.
#'   Annals of Statistics, 6(2), 461-464.
#' - Sclove, S. L. (1987). *Application of model-selection criteria to some
#'   problems in multivariate analysis*. Psychometrika, 52(3), 333-343.
#'
#' @section Conditional maximum likelihood:
#' Set `method = "CML"` and `category_policy = "preserve"` for native many-facet
#' RSM/PCM estimation. Declare the rubric range with `rating_min` and `rating_max`;
#' PCM requires `step_facet`. Non-person facets and category steps have the
#' ordinary sum-to-zero constraints. Supported two-way facet interactions keep
#' zero marginal sums. Unit observation weights are required. The default
#' `fit_scope = "full"` estimates the full identified structural model; the
#' comparison-only route below reports selected contrasts in deficient models.
#' The suggested package `lpSolve` supplies the finite-MLE check.
#' Repeated Person-by-facet cells require correction or an explicit justified
#' `repeated_ratings = "independent"` declaration before CML intervals are supplied.
#'
#' For Person \eqn{p}, let \eqn{s_p=\sum_j y_{pj}} and let \eqn{x_{pjk}}
#' be the cumulative structural design vector for category \eqn{k}, with
#' \eqn{x_{pj0}=0}. Conditioning on the observed roster and \eqn{s_p} gives
#' \deqn{\ell_c(\beta)=\sum_{p,j}x_{pj,y_{pj}}^T\beta-
#'   \sum_p\log [t^{s_p}]\prod_j\sum_{k=0}^K
#'   \exp(x_{pjk}^T\beta)t^k.}
#' A generalized ESF recursion computes the coefficients in log scale,
#' conditional sufficient-statistic means and covariances. The score is observed
#' minus conditional expected statistics; information is their conditional
#' covariance. It contains cross-information for all estimated structural
#' coordinates. No ability quadrature or ability population model is used.
#'
#' Person-constant effects disappear on conditioning. The information at zero
#' checks identification for the actual rosters and totals; an unidentified or
#' numerically rank-sensitive full model is refused. A dynamic-program graph
#' represents conditional support for a linear-program recession check without
#' enumerating response patterns. A detected or inconclusive recession check,
#' or failed stationarity check prevents delivery of a full fit. If those checks
#' pass but the fitted information cannot be evaluated or safely inverted,
#' the finite numerical point estimates and conditional likelihood are retained
#' with a warning. `checks$information_status` is `"unavailable"` and
#' `checks$information_reason` explains why; the joint covariance is an `NA`
#' matrix, not a regularized substitute. Non-fixed targets have missing SE and
#' bounds, `CIEligible = FALSE` and an `InferenceReview` reason in summaries and
#' facet intervals. Small score residuals alone do not certify accurate points
#' in weak directions; examine the information warning before interpreting them.
#' Inverse information, when available, is unregularized. Extreme-total Persons have
#' zero conditional likelihood contribution and remain in `person_totals`.
#'
#' The full-scope result is an `mfrm_cml_fit`. `summary(fit)` contains `locations`, `steps`,
#' `interactions` and checks; these tables use pointwise 95 percent normal
#' intervals. [mfrm_facet_intervals()] accepts locations or named contrasts at
#' another confidence level, with `method = "model"`. Plot the interval object
#' with `plot(ci)` and retain it with `saveRDS()`. `vcov(fit)` is the covariance
#' in free structural coordinates, or an `NA` matrix with an `unavailable_reason`
#' attribute when information fails. Default `confint(fit)` retains missing
#' bounds and reasons in its `diagnostics` attribute. `logLik(fit)` is the conditional log likelihood.
#' Conditional and marginal/joint likelihood values are different objectives.
#'
#' The conditional model must hold for the actual rating roster: conditional
#' on ability and the roster, ratings follow the specified independent RSM/PCM
#' response law. A prespecified incomplete assignment is compatible with this
#' assumption. Response-dependent assignment, stopping or missingness requires
#' a separate justification that the conditional law is preserved. Conditioning
#' alone does not supply that justification or repair local dependence.
#'
#' One sufficient regular regime has independent Persons, fixed structural
#' dimension and category range, bounded ratings per Person, a finite true
#' parameter, and conditional information \eqn{I_N(\beta_0)/N} converging to a
#' positive-definite matrix. For example, the proportions of roster/total strata
#' can stabilize while their combined support identifies every free direction.
#' Then CML is consistent as Person count \eqn{N} grows, even with fixed ratings
#' per Person, and inverse-information Wald intervals are asymptotically valid.
#' Increasing the total number of Persons alone is insufficient: extreme totals
#' add no information, and a direction supported by only a few Persons can
#' remain imprecise. Growing numbers of facet levels/parameters require a
#' separate asymptotic argument. Numerical checks on one fit do not establish
#' these sampling conditions or guarantee finite-sample coverage.
#' It does not estimate abilities, population means or a population SD. In
#' particular, it cannot supply a contrast divided by that fit's estimated
#' population SD; the MML delta-method route remains separate.
#'
#' Anchors, population regression, shrinkage, corrected JML and GPCM are outside
#' this CML route. Full-scope CML also excludes dummy/uncentered facets; the
#' comparison route uses the review's zero origin main facet.
#' Unsupported explicit arguments are rejected even when equal to a default.
#' Ordinary `mfrm_fit` scoring, diagnostics and consolidated reports are not
#' attached to this separate class. These restrictions do not affect imported
#' eRm fits. See Mair and Hatzinger (2007),
#' \doi{10.18637/jss.v020.i09}, for the linear Rasch-family
#' CML framework.
#'
#' @section Comparison-only estimation:
#' An incomplete assignment can identify a contrast even when the full model
#' has aliased parameters. Set `fit_scope = "comparisons"` and supply an explicit
#' `comparison_review` to estimate such contrasts. This route supports the
#' review's unadjusted RSM/PCM models, unit-weight retained ratings, preserved
#' categories, one facet-by-origin interaction, zero origin main effect and
#' centered other facets/interaction margins. JML has free Person locations;
#' MML uses the declared fixed normal or normal population regression model.
#' CML conditions out Person abilities and requires no population distribution.
#' PCM requires an explicit step owner. Repeated ratings retain their declared
#' treatment and every response; a repeated row is not another independent Person.
#'
#' The retained responses, score map and actual model constraints must match
#' the review. Its map is rebuilt. For JML/MML, a column-scaled SVD removes only
#' null directions of the conditional adjacent-logit map. Every observed probability
#' and invariant contrast is preserved; no substantive effect is fixed to zero
#' to identify a previously unidentified target. For MML the transformation
#' preserves the conditional law at every latent node and leaves the variance
#' coordinate free. It therefore preserves the marginal likelihood, but does
#' not prove global marginal identification. The bounded review requires
#' `nrow(A) * ncol(A)` and `ncol(A)^2` no greater than two million, where `A`
#' is the conditional adjacent-logit design matrix.
#' Numerically sensitive or unevaluated maps yield unavailable estimates.
#'
#' JML first checks the observed-category recession cone using the optional
#' `lpSolve` package. A detected recession or an inconclusive/unavailable check
#' withholds all estimates. This route does not deliver extended-MLE contrasts
#' when some parameters diverge. Unsupported category information also blocks
#' estimates. Numerical acceptance requires optimizer convergence, full native
#' gradient sup norm at most `max(1e-4, 10 * reltol)` and reduced observed
#' information eigenvalues with minimum/maximum ratio greater than `1e-8`
#' (the maximum is floored at one). No covariance is inferred from this check.
#' MML population fits use three deterministic starting variances, requiring
#' all starts to converge, target agreement within `1e-4` logits and objective
#' agreement per Person within `1e-6`. Fixed-population MML has one start.
#' MML also checks `2 * quad_points + 1` nodes at the retained point, requiring
#' the same gradient tolerance and objective change per Person at most `1e-6`.
#' These numerical criteria are not sampling-precision guarantees or a search
#' proof of the global maximum. There is no automatic finer-grid refit.
#'
#' `interaction_effect_table(fit)` and `summary(fit)` return the declared
#' comparisons, retaining unidentified and failed rows with `Estimate = NA`
#' and `EstimateReason`. Available JML/MML estimates are `descriptive_only`, in logits.
#' Positive MML conditional ranks keep primary `State = "not_evaluated"`.
#' JML/MML SEs and interval endpoints remain missing; `FormalInferenceEligible`
#' is false. Save the full object or table as RDS to retain weights, empty cells,
#' assumptions and checks. CSV retains visible columns only.
#'
#' For CML, declare `method = "CML"` in the prospective review. Identification
#' uses differences between every possible response pattern at the observed
#' Person totals and rosters, constructed by a dynamic program. A target is
#' identifiable only if its coefficient vector is in that row space. The same
#' dense size bound applies to the support-difference matrix and intermediate
#' state vectors. Extreme totals add no information. The zero origin main facet
#' is coding for an effect eliminated by conditioning, not an equal-ability
#' population assumption. `FreePersonState` remains a separate response-map audit.
#'
#' CML optimizes the conditional likelihood in the identified subspace. This
#' preserves probabilities given the Person total, not arbitrary unconditional
#' category probabilities. If \eqn{\beta=T\gamma+v}, where \eqn{v} is a
#' conditional null direction, an identified target \eqn{C\beta} satisfies
#' \eqn{Cv=0}. Its covariance is
#' \deqn{CT(T^T I_N T)^{-1}T^T C^T.}
#' This expression is invariant to the choice of basis \eqn{T} for the
#' identified coordinates and includes covariance with all estimable nuisance
#' directions. The asymptotic argument above applies in that subspace when it
#' stabilizes and its normalized information has a positive-definite limit.
#' Numerical identification in one sample does not establish those conditions.
#' It checks finite-MLE existence on the complete
#' conditional support and refuses all estimates if any recession is detected
#' or cannot be checked; extended-MLE target inference is outside scope.
#' Stationarity and unregularized information checks match full-scope CML.
#' If all checks pass, identified targets have `EstimateStatus = "conditional_wald"`,
#' `FormalInferenceEligible = TRUE`, `SE`, and pointwise 95 percent `Lower`/`Upper`
#' limits. If only information fails, identified finite numerical points instead
#' have `EstimateStatus = "conditional_point"`, `FormalInferenceEligible = FALSE`,
#' missing SE/bounds and an `EstimateReason` explaining the failure. Unidentified,
#' boundary or unqualified-stationarity targets do not acquire point estimates.
#' `$covariance` retains available joint covariance; unidentified or
#' failed covariance entries remain `NA`. The result also retains
#' `person_totals`. These are asymptotic intervals under the CML assumptions above,
#' not a finite-sample coverage guarantee or a calibrated DFF classification.
#' CML comparison fits estimate no population SD and cannot standardize by one.
#' Use `mfrm_facet_intervals(fit, level = .95)` to reuse their saved targets
#' through the existing interval plots and APA tables, preserving unavailable
#' rows; changing the confidence level does not refit the model.
#' Population and quadrature arguments are rejected for CML.
#'
#' The result is an `mfrm_comparison_fit`, not an ordinary `mfrm_fit`.
#' Individual cell estimates, Person scoring, ordinary fit diagnostics and
#' prediction are not exposed through this class. `representation` and
#' `optimization` retain numerical coordinates for auditing only; their
#' representatives are not identified parameter estimates. Additional anchors,
#' shrinkage, corrected JML, GPCM, alternate MML engines/integration, checkpoints
#' and attached diagnostics are outside this route. Explicit arguments outside
#' data preparation, the reviewed model and basic optimizer/quadrature controls
#' are rejected, including unsupported arguments set to their default values.
#'
#' The estimable-function criterion is discussed by Lenth (2015),
#' <https://journal.r-project.org/articles/RJ-2015-016/>. Recession directions
#' and finite-MLE existence in exponential families are distinguished by
#' Geyer (2009), \doi{10.1214/08-EJS349}. These principles do not
#' establish global MML identification or finite-sample inference here.
#'
#' @return
#' With `method = "JML", jml_penalty = "jeffreys"`, a separate
#' `mfrm_jeffreys_fit`; see [summary.mfrm_jeffreys_fit()] for its point tables,
#' probabilities, numerical checks and unsupported inference.
#' With `method = "CML", fit_scope = "full"`, an `mfrm_cml_fit` with structural estimates,
#' joint `covariance` and `information` (missing where unavailable), `person_totals`, `checks`, preparation,
#' configuration and optimizer controls; see "Conditional maximum likelihood".
#' With `fit_scope = "comparisons"`, an `mfrm_comparison_fit` containing
#' `comparisons`, the prospective `comparison_review`, numerical/category/boundary
#' `checks`, preparation/configuration, controls and audit-only parameter
#' `representation`/`optimization`. CML adds contrast `covariance` and
#' `person_totals`. Its `summary()` returns the comparison table.
#' Otherwise, an object of class `mfrm_fit` (named list) with:
#' - `summary`: one-row model summary including `LogLik`, `Deviance`, canonical
#'   free dimension `Npar`, response-row/weight/Person counts, the versioned
#'   information-criterion contract, Person-based MML `AIC`/`BIC`/`SABIC`,
#'   explicitly descriptive legacy fields when the common panel is ineligible,
#'   and convergence; it also includes user-facing `Method`, engine-facing
#'   `MethodUsed`, MML-engine
#'   fields, terminal-gradient readiness, requested/selected-stage tolerance
#'   settings, and the actual L-BFGS-B `OptimizerFactr` / `OptimizerPgtol`
#'   controls when applicable
#' - `facets$person`: person estimates (`Estimate`; plus `SD` for MML), with
#'   `SourceFitReadiness`, `SourceInferenceReady`, and `EstimateUse` separating
#'   a defined Person summary from the source fit's permission to interpret it.
#'   In particular, a finite prior-regularized MML EAP does not become
#'   reportable when the source fit is blocked.
#' - `facets$others`: facet-level estimates for each facet
#' - `steps`: estimated threshold/step parameters as a one-row-per-step
#'   `tibble` with `Estimate`. Bare fits keep this table as point estimates.
#'   `diagnose_mfrm()` exposes MML observed-information step uncertainty in
#'   `diagnostics$parameter_uncertainty$steps`. Check `CIEligible` and `CIUse`:
#'   finite curvature-based bands alone do not establish ordinary inference. When
#'   `attach_diagnostics = TRUE`, those `SE`, confidence-limit, and status
#'   columns are attached to `fit$steps` when the Hessian is available.
#'   For step-structure quality, also use the step-collapse and disordering
#'   warnings from [diagnose_mfrm()] and [category_structure_report()].
#' - `slopes`: discrimination parameters for `GPCM` fits as a
#'   one-row-per-slope-element `tibble`. `LogEstimate` and `Estimate` retain
#'   the finite optimizer values for compatibility and numerical diagnosis;
#'   `OptimizerLogEstimate` / `OptimizerEstimate` name that role explicitly.
#'   Read `ParameterStatus`, `PrimaryLogEstimate`, `PrimaryEstimate`,
#'   `SEEligible`, `CIEligible`, and `ReasonCodes` before interpretation. A
#'   certified JML slope-only path receives a typed extended-real primary
#'   boundary. `config$boundary_audit` retains the supporting fixed-objective,
#'   joint-path, and terminal-gradient records. A certified path can establish
#'   that a finite JML maximum is unattained for the evaluated case. The
#'   converse is deliberately not used: failure to find a path in the
#'   evaluated families, or retention of a finite optimizer point, does not
#'   establish existence of a finite global maximum for the non-concave GPCM
#'   likelihood. These JML diagnostic records do not supply uncertainty
#'   estimates, MML results or evidence of agreement with other software.
#'   `config$boundary_audit$gpcm_terminal_gradient_stability` reconstructs the
#'   same fixed JML objective and analytic terminal gradient, checks stored
#'   optimizer/polish summaries and deterministic central-difference probes,
#'   and reports gradient norms by free-parameter block. Positive boundary
#'   certificates take precedence over a finite-point zero or small gradient;
#'   otherwise a coherent small gradient is retained-point first-order
#'   evidence only. The numerical tolerance is not a statistical significance
#'   threshold. These checks do not establish a finite global maximum,
#'   boundary absence, uncertainty, external comparability, or readiness.
#'   The conditional JML boundary checks are not reused for MML.
#'   `confint(fit, parm = "slopes")` and `diagnose_mfrm()` separately check
#'   the current local MML solution before supplying approximate pointwise
#'   relative-slope intervals. `CIEligible` and `InferenceReview` record this
#'   output-specific decision, without overriding the global boundary audit.
#'   Ineligible local calculations remain in `Optimizer*SE` / `Optimizer*CI`.
#'   Refresh saved fits through either function; no refitting is necessary.
#'   The identification convention pins the geometric mean of finite optimizer
#'   slopes at 1.
#' - `readiness`: the versioned fit record, five component rows, and current
#'   parameter-level rows. The current parameter slice includes GPCM slopes;
#'   other non-Person parameter classes remain scheduled for later propagation.
#' - `interactions`: model-estimated facet interaction effects and metadata
#'   when `facet_interactions` is supplied
#' - `population`: population-model metadata. Ordinary RSM/PCM and JML fits
#'   keep an inactive record (`active = FALSE`,
#'   `posterior_basis = "legacy_mml"`). Default GPCM MML and active
#'   latent-regression fits store the fitted design matrix, regression
#'   coefficients, residual variance, omission review, the complete-case
#'   estimation table (`person_table`), and the observed-person-aligned
#'   replay/export provenance table retained before complete-case omission
#'   (`person_table_replay`), plus stored categorical `xlevels` / `contrasts`
#'   for model-matrix replay and scoring, together with
#'   `posterior_basis = "population_model"`. `estimation_converged` records
#'   optimizer convergence and `inference_ready` records the separate formal
#'   readiness decision. The older `population$converged` field is retained as
#'   a compatibility alias of `inference_ready`; its basis is recorded in
#'   `population$converged_basis`.
#' - `data_review`: pre-fit Data, Design, Stability, and Reporting readiness
#'   evidence propagated into summaries and plot-interpretation checks
#' - `config`: resolved model configuration used for estimation, including
#'   `config$anchor_review` and the recorded estimation controls
#' - `prep`: preprocessed data/level metadata
#' - `opt`: optimizer result augmented with `optimizer_diagnostics`, the
#'   complete `optimizer_polish` stage history, method-selection metadata, and
#'   evaluation-cache counters. For direct fitting, its core fields originate
#'   from [stats::optim()]; EM additionally records engine-specific diagnostics.
#'
#' @seealso [mfrm_fit_extractors], [mfrm_fit_responses], [diagnose_mfrm()], [estimate_bias()], [build_apa_outputs()],
#'   [gpcm_capability_matrix], [mfrmr_workflow_methods],
#'   [mfrmr_reporting_and_apa]
#' @examples
#' # Load the package
#' library(mfrmr)
#'
#' # Load example ratings and look at the first six rows
#' ratings <- load_mfrmr_data("example_operational")
#' head(ratings)
#'
#' # Fit the model
#' fit <- fit_mfrm(
#'   data = ratings,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   rating_min = 1, rating_max = 4, category_policy = "preserve",
#'   method = "MML",
#'   model = "RSM"
#' )
#'
#' # Plot the results (Wright map)
#' plot(fit)
#'
#' # Save the summary, then display its tables
#' results <- summary(fit)
#' results$person_overview # One row summarizing person ability estimates
#' results$facet_overview  # One row per facet: number of levels, mean, SD, range
#'
#' # Check the interpretation status and recommended next step
#' results$decision
#'
#' \dontrun{
#' # Optional repeated fitting: compare specified starting values.
#' reviewed <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   mml_starts = c(0, 0.5, -0.5))
#' summary(reviewed)$starting_values # Every attempt, including warnings/failures
#' reviewed$opt$start_review$parameter_map # Meaning and constraints of coordinates
#' write.csv(summary(reviewed)$starting_values, "starting-values.csv", row.names = FALSE)
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
                     missing_codes = NULL,
                     model = c("RSM", "PCM", "GPCM"),
                     method = c("MML", "JML", "JMLE", "CML"),
                     step_facet = NULL,
                     slope_facet = NULL,
                     facet_interactions = NULL,
                     min_obs_per_interaction = 10,
                     interaction_policy = c("warn", "error", "silent"),
                     anchors = NULL,
                     group_anchors = NULL,
                     noncenter_facet = "Person",
                     dummy_facets = NULL,
                     positive_facets = NULL,
                     anchor_policy = c("warn", "error", "silent"),
                     min_common_anchors = 5L,
                     min_obs_per_element = 30,
                     min_obs_per_category = 10,
                     quad_points = 31,
                     maxit = 400,
                     reltol = 1e-9,
                     optimizer = c("auto", "BFGS", "L-BFGS-B"),
                     mml_engine = c("direct", "em", "hybrid"),
                     population_formula = NULL,
                     person_data = NULL,
                     person_id = NULL,
                     population_policy = c("error", "omit"),
                     facet_shrinkage = c("none", "empirical_bayes", "laplace"),
                     facet_prior_sd = NULL,
                     shrink_person = FALSE,
                     attach_diagnostics = FALSE,
                     checkpoint = NULL,
                     gpcm_mml_identification = c("free_population", "fixed_standard_normal"),
                     mml_integration = c("fixed", "adaptive", "auto"),
                     category_policy = NULL,
                     em_score_tol = NULL,
                     jml_correction_order = NULL,
                     jml_correction_sampling = c("fixed_rosters", "random_rosters"),
                     gpcm_mml_start = NULL,
                     repeated_ratings = c("review", "independent"),
                     centering_scope = c("facet", "remainder", "anchor_origin"),
                     fit_scope = c("full", "comparisons"),
                     comparison_review = NULL,
                     mml_accuracy_control = NULL,
                     population_sd = NULL,
                     jml_penalty = c("none", "jeffreys"),
                     progress = FALSE,
                     mml_starts = NULL) {
  close_progress <- mfrmr_start_fit_progress(progress)
  on.exit(close_progress(), add = TRUE)
  supplied_arguments <- setdiff(names(as.list(match.call())[-1L]), "progress")
  if (is.null(mml_starts)) supplied_arguments <- setdiff(supplied_arguments, "mml_starts")
  fit_scope <- match.arg(fit_scope)
  if (fit_scope == "full" && !is.null(comparison_review))
    stop("Use fit_scope = 'comparisons' to fit a comparison review, or extract comparisons from a full fit with interaction_effect_table().", call. = FALSE)
  centering_scope <- match.arg(centering_scope)
  repeated_ratings <- match.arg(repeated_ratings)
  keep_original <- resolve_mfrm_category_policy(
    keep_original, category_policy, !missing(keep_original))
  # If users opt in to preparation messages, suppress duplicates that would
  # otherwise fire once in review_mfrm_anchors() and again in mfrm_estimate()
  # -> prepare_mfrm_data(). By default, preparation provenance is stored in the
  # returned object and displayed by summary/data-description helpers instead
  # of being emitted as routine messages.
  prior_announce_opt <- getOption("mfrmr._rating_range_announced")
  options(mfrmr._rating_range_announced = FALSE)
  on.exit(options(mfrmr._rating_range_announced = prior_announce_opt),
          add = TRUE)
  prior_prep_message_opt <- getOption("mfrmr._preparation_messages_announced")
  options(mfrmr._preparation_messages_announced = FALSE)
  on.exit(options(mfrmr._preparation_messages_announced = prior_prep_message_opt),
          add = TRUE)
  prior_score_recode_opt <- getOption("mfrmr._score_recode_announced")
  options(mfrmr._score_recode_announced = FALSE)
  on.exit(options(mfrmr._score_recode_announced = prior_score_recode_opt),
          add = TRUE)
  prior_duplicate_warning_opt <- getOption("mfrmr._duplicate_warning_announced")
  options(mfrmr._duplicate_warning_announced = FALSE)
  on.exit(options(mfrmr._duplicate_warning_announced = prior_duplicate_warning_opt),
          add = TRUE)
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
  integer_tolerance <- sqrt(.Machine$double.eps)
  if (!is.numeric(maxit) || length(maxit) != 1 ||
      !is.finite(maxit) || maxit < 1 ||
      abs(maxit - round(maxit)) > integer_tolerance) {
    stop("`maxit` must be a finite positive integer. Got: ",
         deparse(maxit), ".", call. = FALSE)
  }
  if (!is.numeric(reltol) || length(reltol) != 1 ||
      !is.finite(reltol) || reltol <= 0) {
    stop("`reltol` must be a finite positive number. Got: ",
         deparse(reltol), ".", call. = FALSE)
  }
  if (!is.numeric(quad_points) || length(quad_points) != 1 ||
      !is.finite(quad_points) || quad_points < 1 ||
      abs(quad_points - round(quad_points)) > integer_tolerance) {
    stop("`quad_points` must be a finite positive integer. Got: ",
         deparse(quad_points), ".", call. = FALSE)
  }
  if (!is.logical(attach_diagnostics) || length(attach_diagnostics) != 1L ||
      is.na(attach_diagnostics)) {
    stop("`attach_diagnostics` must be a single logical value (TRUE or FALSE).",
         call. = FALSE)
  }
  if (!is.null(person_id) && (!is.character(person_id) || length(person_id) != 1 || !nzchar(person_id))) {
    stop("`person_id` must be NULL or a single non-empty character string naming the person column in `person_data`.",
         call. = FALSE)
  }
  if (!is.null(slope_facet) && (!is.character(slope_facet) || !length(slope_facet) ||
                              anyNA(slope_facet) || any(!nzchar(slope_facet)))) {
    stop("`slope_facet` must be NULL or non-empty character strings naming facet columns.",
         call. = FALSE)
  }

  model <- toupper(match.arg(model))
  method_input <- public_mfrm_method_label(toupper(match.arg(method)))
  method <- method_input
  if (identical(centering_scope, "anchor_origin") &&
      (!identical(method, "JML") || !model %in% c("RSM", "PCM") ||
       length(facet_interactions) > 0L || fit_scope != "full")) {
    stop("`centering_scope = 'anchor_origin'` requires ordinary additive RSM/PCM JML with fit_scope = 'full'.",
         call. = FALSE)
  }
  optimizer <- normalize_mfrm_optimizer(match.arg(optimizer))
  mml_engine <- tolower(match.arg(mml_engine))
  mml_integration <- match.arg(mml_integration)
  jml_penalty <- match.arg(jml_penalty)
  if (!is.null(population_sd)) {
    if (!is.numeric(population_sd) || is.complex(population_sd) ||
        length(population_sd) != 1L || !is.finite(population_sd) ||
        population_sd <= 0 || !is.finite(population_sd^2) || population_sd^2 <= 0)
      stop("`population_sd` must be NULL or one positive finite number with a representable positive variance.", call. = FALSE)
    if (!model %in% c("RSM", "PCM") || method != "MML" ||
        mml_engine != "direct" || is.null(population_formula) ||
        fit_scope != "full" || length(facet_interactions))
      stop("`population_sd` requires a population_formula with full direct RSM/PCM MML and no interactions.", call. = FALSE)
    population_sd <- as.numeric(population_sd)
  }
  if (identical(mml_integration, "auto")) {
    product_auto <- identical(model, "GPCM") && length(slope_facet) == 2L &&
      is.null(jml_correction_order) && identical(jml_penalty, "none")
    if (!identical(method, "MML") || !is.null(checkpoint) || fit_scope != "full" ||
        (!product_auto && (!identical(mml_engine, "direct") || length(slope_facet) > 1L)) ||
        (product_auto && !mml_engine %in% c("em", "direct"))) {
      stop("`mml_integration = 'auto'` requires ordinary full MML with the direct engine, or two-family GPCM with the em/direct engine, and no checkpoint.",
        call. = FALSE)
    }
    mml_accuracy_control <- if (product_auto) mfrmr_gmfrm_accuracy_control(mml_accuracy_control) else
      mfrmr_mml_accuracy_control(mml_accuracy_control)
    if (!is.null(mml_accuracy_control$max_points)) {
      maximum_initial <- if (product_auto) floor((mml_accuracy_control$max_points + 1) / 2) else
        mml_accuracy_control$max_points
      if ("quad_points" %in% supplied_arguments && quad_points > maximum_initial)
        stop("`quad_points` exceeds the initial order allowed by `max_points`, including reference checks.", call. = FALSE)
      quad_points <- min(quad_points, maximum_initial)
    }
  } else if (!is.null(mml_accuracy_control)) {
    stop("`mml_accuracy_control` requires `mml_integration = 'auto'`.", call. = FALSE)
  }
  if (identical(mml_integration, "adaptive") &&
      (!identical(method, "MML") || !identical(mml_engine, "direct") || !is.null(checkpoint))) {
    stop("`mml_integration = 'adaptive'` requires MML with the direct engine and no checkpoint.",
         call. = FALSE)
  }
  if (!is.null(mml_starts) &&
      (method != "MML" || fit_scope != "full" || length(slope_facet) > 1L ||
       !is.null(jml_correction_order) || jml_penalty != "none" ||
       !is.null(checkpoint) || mml_integration == "auto"))
    stop("`mml_starts` requires ordinary full MML, fixed/adaptive integration and no checkpoint; two-family GPCM uses gpcm_mml_start.", call. = FALSE)
  gpcm_mml_identification <- match.arg(gpcm_mml_identification)
  interaction_policy <- tolower(match.arg(interaction_policy))
  anchor_policy <- tolower(match.arg(anchor_policy))
  population_policy <- tolower(match.arg(population_policy))
  facet_shrinkage <- match.arg(facet_shrinkage)
  if (jml_penalty != "none") {
    arguments <- mget(names(formals(fit_mfrm)), envir = environment(), inherits = FALSE)
    return(mfrm_fit_jeffreys(arguments, supplied_arguments))
  }
  if (!is.null(facet_prior_sd)) {
    if (!is.numeric(facet_prior_sd) || length(facet_prior_sd) != 1L ||
        !is.finite(facet_prior_sd) || facet_prior_sd < 0 ||
        !is.finite(facet_prior_sd^2) ||
        (facet_prior_sd > 0 && facet_prior_sd^2 == 0)) {
      stop("`facet_prior_sd` must be NULL or a non-negative finite scalar whose square is representable (positive when the SD is positive).",
           call. = FALSE)
    }
  }
  if (!is.logical(shrink_person) || length(shrink_person) != 1L ||
      is.na(shrink_person)) {
    stop("`shrink_person` must be a single logical value.", call. = FALSE)
  }
  if (!identical(facet_shrinkage, "none") && isTRUE(shrink_person) &&
      (!identical(method, "JML") || !is.null(jml_correction_order))) {
    stop("`shrink_person = TRUE` requires ordinary JML estimates with sampling SEs. EAP posterior SDs are not sampling SEs. Use shrink_person = FALSE to retain the original Person estimates.", call. = FALSE)
  }

  if (method == "CML" && fit_scope == "full") {
    arguments <- mget(names(formals(fit_mfrm)), envir = environment(), inherits = FALSE)
    return(mfrmr_fit_cml(arguments, supplied_arguments))
  }

  if (fit_scope == "comparisons") {
    arguments <- mget(names(formals(fit_mfrm)), envir = environment(), inherits = FALSE)
    return(mfrmr_fit_comparisons(arguments, supplied_arguments))
  }

  supplied_arguments <- setdiff(supplied_arguments, c("fit_scope", "comparison_review"))

  if (!is.null(gpcm_mml_start) &&
      !(length(slope_facet) == 2L && identical(model, "GPCM") &&
        identical(method_input, "MML") && identical(mml_engine, "direct") &&
        mml_integration %in% c("adaptive", "auto") && is.null(jml_correction_order)))
    stop("`gpcm_mml_start` applies only to two-family adaptive direct MML.", call. = FALSE)
  if (!is.null(jml_correction_order)) {
    arguments <- mget(names(formals(fit_mfrm)), envir = environment(), inherits = FALSE)
    return(mfrm_fit_adjusted_jml(arguments, supplied_arguments))
  }
  if (!missing(jml_correction_sampling))
    stop("`jml_correction_sampling` requires an explicit `jml_correction_order`.", call. = FALSE)
  if (length(slope_facet) > 1L) {
    arguments <- mget(names(formals(fit_mfrm)), envir = environment(), inherits = FALSE)
    if (identical(mml_integration, "auto"))
      return(mfrm_fit_product_slopes_auto(arguments, supplied_arguments))
    return(mfrm_fit_product_slopes(arguments, supplied_arguments))
  }
  if (!is.null(em_score_tol)) {
    stop("`em_score_tol` applies only to the two-family GPCM MML-EM route.", call. = FALSE)
  }

  gpcm_mml_active <- identical(model, "GPCM") && identical(method_input, "MML")
  if (gpcm_mml_active &&
      identical(gpcm_mml_identification, "fixed_standard_normal") &&
      !is.null(population_formula)) {
    stop(
      "`gpcm_mml_identification = \"fixed_standard_normal\"` requires ",
      "`population_formula = NULL`. Use `\"free_population\"` when estimating ",
      "a GPCM population model.",
      call. = FALSE
    )
  }

  auto_gpcm_population <- gpcm_mml_active &&
    identical(gpcm_mml_identification, "free_population") &&
    is.null(population_formula)
  effective_population_formula <- population_formula
  effective_person_data <- person_data
  effective_person_id <- person_id
  population_scaffold_data <- data
  estimation_person_ids <- as.character(data[[person]])
  if (auto_gpcm_population) {
    estimation_person_ids <- trimws(estimation_person_ids)
    default_missing_codes <- isTRUE(missing_codes) ||
      (is.character(missing_codes) && length(missing_codes) == 1L &&
       identical(tolower(missing_codes), "default"))
    if (!is.null(missing_codes) && !isFALSE(missing_codes) &&
        !default_missing_codes) {
      estimation_person_ids[
        estimation_person_ids %in% trimws(as.character(missing_codes))
      ] <- NA_character_
    }
    population_scaffold_data[[person]] <- estimation_person_ids
    population_scaffold_data <- population_scaffold_data[
      !is.na(estimation_person_ids), , drop = FALSE
    ]
    effective_population_formula <- stats::as.formula("~ 1")
    effective_person_data <- data.frame(
      .mfrmr_person_id = unique(estimation_person_ids[!is.na(estimation_person_ids)]),
      stringsAsFactors = FALSE
    )
    names(effective_person_data) <- person
    effective_person_id <- person
  }

  population <- prepare_mfrm_population_scaffold(
    data = population_scaffold_data,
    person = person,
    population_formula = effective_population_formula,
    person_data = effective_person_data,
    person_id = effective_person_id,
    population_policy = population_policy
  )
  if (!is.null(population_sd)) population$fixed_sd <- population_sd
  population$source <- if (auto_gpcm_population) {
    "gpcm_mml_default_identification"
  } else if (isTRUE(population$active)) {
    "user_supplied"
  } else {
    "none"
  }
  population$identification_role <- if (gpcm_mml_active &&
                                         isTRUE(population$active)) {
    "free_population_scale_with_geometric_mean_one_relative_slopes"
  } else if (gpcm_mml_active) {
    "fixed_standard_normal_with_geometric_mean_one_relative_slopes"
  } else {
    "not_applicable"
  }
  if (auto_gpcm_population) {
    population$response_rows_retained <- as.integer(nrow(population_scaffold_data))
    population$response_rows_omitted <- as.integer(
      nrow(data) - nrow(population_scaffold_data)
    )
    population$notes <- unique(c(
      as.character(population$notes %||% character(0)),
      paste(
        "An intercept-only person population model was activated automatically",
        "to identify GPCM MML without fixing both latent variance and",
        "the slope geometric mean."
      ),
      if (population$response_rows_omitted > 0L) {
        paste0(
          population$response_rows_omitted,
          " response row(s) with missing person IDs were excluded from the ",
          "automatic population scaffold; ordinary response preparation ",
          "applies the same exclusion."
        )
      } else {
        character(0)
      }
    ))
  }
  if (isTRUE(population$active)) {
    if (!identical(method_input, "MML")) {
      stop("Latent-regression estimation requires `method = 'MML'`. ",
           "The requested population model is not available with JML.",
           call. = FALSE)
    }
    if (!identical(as.character(noncenter_facet[1]), "Person")) {
      if (auto_gpcm_population) {
        stop(
          "Default GPCM MML identification estimates an intercept-only ",
          "population model and therefore requires `noncenter_facet = \"Person\"`. ",
          "Use the default centering, or request ",
          "`gpcm_mml_identification = \"fixed_standard_normal\"` only when ",
          "the legacy restricted likelihood is intended.",
          call. = FALSE
        )
      }
      if (!model %in% c("RSM", "PCM") || is.null(population_sd) ||
          !identical(population$design_columns, "(Intercept)") ||
          any(population$design_matrix != 1) || length(dummy_facets) ||
          (!is.null(group_anchors) && NROW(group_anchors) > 0L))
        stop("Latent-regression identification requires `noncenter_facet = \"Person\"`, or an intercept-only RSM/PCM population with explicit `population_sd` and a direct anchor on the noncentered facet, without groups or dummy facets.", call. = FALSE)
      population$identification_role <- "direct_anchor_origin_fixed_population_sd"
    }
  }
  estimation_data <- data
  if (isTRUE(population$active) && !auto_gpcm_population &&
      length(population$included_persons) > 0) {
    estimation_mask <- estimation_person_ids %in% population$included_persons
    estimation_data <- data[estimation_mask, , drop = FALSE]
  }

  anchor_review <- review_mfrm_anchors(
    data = estimation_data,
    person = person,
    facets = facets,
    score = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight = weight,
    keep_original = keep_original,
    missing_codes = missing_codes,
    anchors = anchors,
    group_anchors = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    repeated_ratings = repeated_ratings,
    centering_scope = centering_scope
  )

  schema_issues <- anchor_review$issues[c("anchor_schema_mismatch", "group_anchor_schema_mismatch")]
  if (any(vapply(schema_issues, nrow, integer(1)) > 0L))
    stop("Malformed anchor schema: nonempty constraints cannot be discarded. Direct anchors require facet, level and anchor/value/measure; group anchors also require group/subset and groupvalue/value/anchor.", call. = FALSE)
  anchors <- anchor_review$anchors
  group_anchors <- anchor_review$group_anchors

  issue_counts <- anchor_review$issue_counts
  issue_total <- if (is.null(issue_counts) || nrow(issue_counts) == 0) 0L else sum(issue_counts$N, na.rm = TRUE)
  if (issue_total > 0) {
    msg <- format_anchor_review_message(anchor_review)
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
    missing_codes = missing_codes,
    model = model,
    method = method,
    step_facet = step_facet,
    slope_facet = slope_facet,
    facet_interactions = facet_interactions,
    min_obs_per_interaction = min_obs_per_interaction,
    interaction_policy = interaction_policy,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    population = population,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol,
    optimizer = optimizer,
    mml_engine = mml_engine,
    mml_integration = mml_integration,
    mml_accuracy_control = mml_accuracy_control,
    mml_starts = mml_starts,
    checkpoint = checkpoint,
    repeated_ratings = repeated_ratings,
    centering_scope = centering_scope
  )

  fit$config$anchor_review <- anchor_review
  fit$config$method_input <- method_input
  fit$population <- finalize_mfrm_population_fit(fit, population)
  fit$config$population_spec <- compact_population_spec(fit$population, fit$prep$levels$Person)
  fit$config$population_active <- isTRUE(population$active)
  fit$config$posterior_basis <- as.character(fit$population$posterior_basis %||% "legacy_mml")
  fit$config$population_policy <- fit$population$policy %||% NULL
  fit$config$gpcm_mml_identification_requested <-
    as.character(gpcm_mml_identification)
  fit$config$gpcm_mml_identification <- if (gpcm_mml_active) {
    as.character(gpcm_mml_identification)
  } else if (identical(model, "GPCM") && identical(method_input, "JML")) {
    "not_applicable_jml"
  } else {
    "not_applicable"
  }
  fit$config$gpcm_common_discrimination <- if (
    gpcm_mml_active && identical(gpcm_mml_identification, "free_population")
  ) {
    "estimated_via_population_sd"
  } else if (gpcm_mml_active) {
    "fixed_to_one_legacy_restriction"
  } else {
    "not_applicable"
  }
  fit$config$gpcm_estimator_family <- if (!identical(model, "GPCM")) {
    "not_applicable"
  } else if (identical(method_input, "JML")) {
    "unpenalized_identified_jml"
  } else {
    "marginal_maximum_likelihood"
  }
  gpcm_identity <- mfrmr_gpcm_model_identity(model)
  fit$config$gpcm_model_family <- gpcm_identity$model_family
  fit$config$gpcm_slope_action <- gpcm_identity$slope_action
  fit$config$gpcm_slope_composition <- gpcm_identity$slope_composition
  fit$config$gpcm_latent_dimension_count <-
    gpcm_identity$latent_dimension_count
  # These fields describe the statistical objective, not numerical line-search
  # safeguards. The current GPCM routes use neither a regularization penalty
  # nor finite box constraints on person, location, step, or slope coordinates.
  fit$config$gpcm_statistical_penalty <- if (identical(model, "GPCM")) {
    "none"
  } else {
    "not_applicable"
  }
  fit$config$gpcm_finite_parameter_box <- if (identical(model, "GPCM")) {
    FALSE
  } else {
    NA
  }
  fit$config$gpcm_extreme_person_policy <- if (!identical(model, "GPCM")) {
    "not_applicable"
  } else if (identical(method_input, "JML")) {
    "extended_real_primary_when_certified"
  } else {
    "posterior_eap_under_population_model"
  }
  fit$config$population_formula <- if (!is.null(fit$population$formula)) {
    paste(deparse(fit$population$formula), collapse = " ")
  } else {
    NULL
  }

  if (identical(model, "GPCM") && identical(method_input, "MML") &&
      nrow(fit$slopes %||% data.frame()) > 0L) {
    fitted_population_sd <- if (isTRUE(fit$population$active)) {
      sqrt(as.numeric(fit$population$sigma2 %||% NA_real_)[1])
    } else {
      1
    }
    optimizer_relative <- as.numeric(
      fit$slopes$OptimizerEstimate %||% fit$slopes$Estimate
    )
    primary_relative <- as.numeric(
      fit$slopes$PrimaryEstimate %||% rep(NA_real_, nrow(fit$slopes))
    )
    fit$slopes$PopulationSD <- rep(fitted_population_sd, nrow(fit$slopes))
    fit$slopes$FixedLatentSDOptimizerEstimate <-
      fitted_population_sd * optimizer_relative
    fit$slopes$FixedLatentSDPrimaryEstimate <-
      fitted_population_sd * primary_relative
    fit$slopes$FixedLatentSDBasis <- if (isTRUE(fit$population$active)) {
      "population_sd_times_relative_slope"
    } else {
      "legacy_unit_population_sd_times_relative_slope"
    }
  }

  class(fit) <- c("mfrm_fit", class(fit))

  # Optional empirical-Bayes / Laplace shrinkage applied post-fit
  # post-hoc on the fitted facet estimates. The `"none"` default
  # preserves the no-shrinkage behavior; other choices add
  # `ShrunkEstimate` / `ShrunkSE` / `ShrinkageFactor` columns and a
  # `shrinkage_report` table. See `?apply_empirical_bayes_shrinkage`
  # and the "Empirical-Bayes shrinkage" section of `?fit_mfrm`.
  if (!identical(facet_shrinkage, "none")) {
    fit <- .apply_shrinkage_to_fit(
      fit = fit,
      method = facet_shrinkage,
      facet_prior_sd = facet_prior_sd,
      shrink_person = isTRUE(shrink_person)
    )
  } else {
    fit$config$facet_shrinkage <- "none"
    fit$config$facet_prior_sd <- NULL
  }

  # Capture every `fit_mfrm()` argument that affects the fit so that
  # `export_mfrm_bundle()` can write a complete replay script. We
  # store inputs as supplied (post `match.arg`) rather than rederiving
  # them from `fit$config`, because some arguments (e.g.
  # `missing_codes`, `min_obs_per_*`, `anchor_policy`) only take
  # effect during preparation and are not otherwise echoed back.
  fit$config$replay_inputs <- list(
    person = as.character(person),
    facets = as.character(facets),
    score = as.character(score),
    weight = if (is.null(weight)) NULL else as.character(weight),
    rating_min = rating_min,
    rating_max = rating_max,
    keep_original = isTRUE(keep_original),
    missing_codes = missing_codes,
    model = as.character(model),
    method = as.character(method_input),
    step_facet = if (is.null(step_facet)) NULL else as.character(step_facet),
    slope_facet = if (is.null(slope_facet)) NULL else as.character(slope_facet),
    gpcm_mml_identification = as.character(gpcm_mml_identification),
    facet_interactions = facet_interactions,
    min_obs_per_interaction = as.numeric(min_obs_per_interaction),
    interaction_policy = as.character(interaction_policy),
    anchors = anchors,
    group_anchors = group_anchors,
    noncenter_facet = as.character(noncenter_facet),
    centering_scope = centering_scope,
    dummy_facets = if (length(dummy_facets) > 0L) as.character(dummy_facets) else NULL,
    positive_facets = if (length(positive_facets) > 0L) as.character(positive_facets) else NULL,
    anchor_policy = as.character(anchor_policy),
    min_common_anchors = as.integer(min_common_anchors),
    min_obs_per_element = as.numeric(min_obs_per_element),
    min_obs_per_category = as.numeric(min_obs_per_category),
    quad_points = as.integer(quad_points),
    maxit = as.integer(maxit),
    reltol = as.numeric(reltol),
    optimizer = as.character(optimizer),
    mml_engine = as.character(mml_engine),
    mml_integration = mml_integration,
    mml_accuracy_control = mml_accuracy_control,
    population_formula = population_formula,
    population_sd = population_sd,
    person_id = if (is.null(person_id)) NULL else as.character(person_id),
    population_policy = as.character(population_policy),
    facet_shrinkage = as.character(facet_shrinkage),
    facet_prior_sd = facet_prior_sd,
    shrink_person = isTRUE(shrink_person),
    attach_diagnostics = isTRUE(attach_diagnostics),
    repeated_ratings = repeated_ratings,
    package_version = as.character(utils::packageVersion("mfrmr"))
  )

  if (!is.null(mml_starts)) fit$config$replay_inputs$mml_starts <- mml_starts

  if (isTRUE(attach_diagnostics)) {
    fit <- attach_diagnostics_to_fit(fit)
  } else {
    fit$config$attached_diagnostics <- FALSE
  }

  fit
}

#' Extract model-estimated facet interaction effects
#'
#' `interaction_effect_table()` returns the fixed-effect interaction block
#' estimated by [fit_mfrm()] when `facet_interactions` is supplied. These are
#' model-estimated deviations from the additive main-effects MFRM, not the
#' residual screening statistics returned by [estimate_bias()].
#'
#' @param fit An `mfrm_fit` or `mfrm_comparison_fit` returned by [fit_mfrm()].
#' @param comparison_review Optional result of [describe_mfrm_data()] with
#'   `comparisons`, or its `$comparison_review` component. Returns the declared
#'   interaction contrasts from a matching joint fit instead of individual cells.
#'
#' @details
#' mfrmr supports two-way interactions between non-person facets,
#' for example `facet_interactions = "Rater:Criterion"`. Each interaction matrix
#' is identified by zero marginal sums across both participating facets, so the
#' interaction estimates are separable from the two main effects. Positive values
#' indicate higher-than-expected scores for the facet-level combination under the
#' additive model; negative values indicate lower-than-expected scores.
#'
#' Use this table for confirmatory model review after specifying the facet pair
#' of substantive interest. For exploratory screening without adding parameters
#' to the fitted model, use [estimate_bias()] or [estimate_all_bias()].
#'
#' With `comparison_review`, first fit the declared model using the review's
#' `fit_arguments` and the original ratings. The fit must match the retained
#' responses, origins, score coding, model, constraints, facet orientation,
#' step owner and population design. Only native, unadjusted RSM/PCM fits
#' without shrinkage are supported. The contrast is the declared weighted
#' sum of Group1-minus-Group2 interaction effects, in logits. Positive means
#' a relative shift toward higher coded scores; it is not an absolute school
#' effect or evidence of unfairness. Location and step effects are not added.
#'
#' This route reevaluates the retained likelihood and full gradient without
#' optimization. Objective agreement must be within `1e-6` and the gradient
#' sup norm within `max(1e-4, 10 * reltol)`, matching the ordinary numerical
#' stopping tolerance. For MML it also evaluates `2 * q + 1` quadrature nodes,
#' requiring the finer-grid gradient to meet that tolerance and the absolute
#' negative-log-likelihood change per Person to be at most `1e-6`.
#' These are numerical checks, not statistical precision criteria. A failed
#' check retains the row with `Estimate = NA` and an explanation; the separate
#' `OptimizerEstimate` is a diagnostic trace. No refit or partial fit is run.
#'
#' `EstimateStatus = "descriptive_only"` requires full conditional-map rank,
#' numerical convergence, current readiness information and no recorded
#' blocking condition or boundary exclusion. It does not certify finite-MLE
#' existence, global marginal identification or sampling performance. MML
#' retains `State = "not_evaluated"` for positive conditional ranks even when
#' a descriptive fitted contrast is available. `BoundaryState`,
#' `FitReadiness`, `NumericalState` and `QuadratureState` remain separate.
#' `SE`, `Lower`, `Upper` remain missing and `FormalInferenceEligible` remains
#' false. A target-only identified contrast in a rejected full model remains
#' unavailable through this route. Save the returned table with `saveRDS()`
#' to preserve its review and checks without recalculation; CSV retains the
#' visible columns but not the attributes.
#'
#' For a `mfrm_comparison_fit` created with `fit_scope = "comparisons"`, omit
#' `comparison_review`: the table already contains the declared targets and
#' their partial-fit checks. Extraction does not fit again. See [fit_mfrm()]
#' for its likelihood-preserving parameterization, availability requirements
#' and limits. For full RSM/PCM MML fits, [mfrm_facet_intervals()] with the same
#' `comparison_review` separately supplies model-information intervals and
#' supported population-SD standardization; extraction here stays descriptive.
#' CML comparison fits add conditional-information SEs and pointwise
#' 95 percent Wald limits for identified targets, with `conditional_wald` status;
#' unidentified rows remain missing. The ordinary full-fit route above remains unchanged.
#'
#' @return A tibble with one row per interaction cell. Returns an empty tibble
#'   when the fit has no model-estimated facet interactions. With
#'   `comparison_review`, one row per requested contrast with availability,
#'   numerical and readiness fields; attributes `comparison_review`,
#'   `fit_checks` and `notes` preserve the prospective question, source
#'   fingerprint and numerical evidence, and interpretation limits.
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [summary.mfrm_fit()],
#'   [estimate_bias()], [compare_mfrm()], [mfrmr_linking_and_dff]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("interaction_effect_table", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(
#'   toy, person = "Person", facets = c("Rater", "Criterion"),
#'   score = "Score", method = "MML", model = "RSM",
#'   facet_interactions = "Rater:Criterion",
#'   quad_points = 7, maxit = 30
#' )
#' head(interaction_effect_table(fit))
#' }
#' @seealso [mfrm_facet_intervals()] with `interaction` for supported
#'   fixed-population RSM/PCM MML cell intervals from full joint information.
#' @export
interaction_effect_table <- function(fit, comparison_review = NULL) {
  if (inherits(fit, "mfrm_comparison_fit")) {
    if (!is.null(comparison_review))
      stop("A comparison fit retains its declared targets; omit comparison_review when extracting it.", call. = FALSE)
    return(fit$comparisons)
  }
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object returned by fit_mfrm().",
         call. = FALSE)
  }
  if (!is.null(comparison_review)) {
    return(mfrmr_fitted_comparison_table(fit, comparison_review))
  }
  tbl <- fit$interactions$effects
  if (is.null(tbl)) {
    return(tibble::tibble())
  }
  tibble::as_tibble(tbl)
}

# Internal: merge per-level fit statistics from diagnose_mfrm() onto
# fit$facets$others so downstream code can read SE / Infit / Outfit
# alongside Estimate without re-running the diagnostics pass.
attach_diagnostics_to_fit <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit.", call. = FALSE)
  }
  others <- fit$facets$others
  person <- fit$facets$person
  if ((is.null(others) || nrow(others) == 0L) &&
      (is.null(person) || nrow(person) == 0L)) {
    fit$config$attached_diagnostics <- TRUE
    return(fit)
  }
  diag <- tryCatch(
    suppressMessages(suppressWarnings(
      diagnose_mfrm(fit, residual_pca = "none")
    )),
    error = function(e) NULL
  )
  if (is.null(diag) || is.null(diag$measures)) {
    fit$config$attached_diagnostics <- FALSE
    fit$config$attached_diagnostics_note <-
      "diagnose_mfrm() did not return a measures table; attach skipped."
    return(fit)
  }
  structural_uncertainty <- diag$parameter_uncertainty %||% list()
  m <- as.data.frame(diag$measures, stringsAsFactors = FALSE)
  # diagnose_mfrm() exposes the point-measure correlation under the key
  # `PTMEA` (FACETS fixed-width label); the attach layer renames it to
  # `PtMeaCorr` so fit$facets$others matches the FACETS Table 7 column
  # naming the user most often expects.
  merge_cols <- intersect(
    c("ModelSE", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD",
      "PTMEA", "PtMeaCorr", "InferenceReady", "SupportsFormalInference",
      "SEUse", "CIEligible", "CIUse", "CIBasis", "CILabel"),
    names(m)
  )
  if (length(merge_cols) == 0L) {
    fit$config$attached_diagnostics <- FALSE
    return(fit)
  }
  keep <- c("Facet", "Level", merge_cols)
  m <- m[, keep, drop = FALSE]
  m$Facet <- as.character(m$Facet)
  m$Level <- as.character(m$Level)
  # Rename to FACETS Table 7 conventions: ModelSE -> SE, PTMEA -> PtMeaCorr.
  if ("ModelSE" %in% names(m)) {
    names(m)[names(m) == "ModelSE"] <- "SE"
  }
  if ("PTMEA" %in% names(m)) {
    names(m)[names(m) == "PTMEA"] <- "PtMeaCorr"
  }
  attached_cols <- setdiff(names(m), c("Facet", "Level"))

  # Non-person facets table.
  if (!is.null(others) && nrow(others) > 0L) {
    others$Facet <- as.character(others$Facet)
    others$Level <- as.character(others$Level)
    others_m <- m[m$Facet != "Person", , drop = FALSE]
    # Drop any columns that are being newly attached to avoid duplicates
    # when the helper is called a second time (idempotent merge).
    others[, intersect(attached_cols, names(others))] <- NULL
    others <- merge(others, others_m, by = c("Facet", "Level"),
                    all.x = TRUE, sort = FALSE)
    fit$facets$others <- others
  }

  # Person facet table (added so summary(fit) can show per-person fit).
  if (!is.null(person) && nrow(person) > 0L) {
    person_df <- as.data.frame(person, stringsAsFactors = FALSE)
    person_df$Person <- as.character(person_df$Person)
    person_m <- m[m$Facet == "Person", , drop = FALSE]
    if (nrow(person_m) > 0L) {
      person_m$Person <- as.character(person_m$Level)
      person_m <- person_m[, c("Person", attached_cols), drop = FALSE]
      # When the existing person table already carries a model SE column,
      # prefer it (e.g. MML posterior SE) and only attach the diagnostic
      # SE if the existing column is empty / NA.
      person_attach_cols <- attached_cols
      if ("SE" %in% person_attach_cols && "SE" %in% names(person_df) &&
          any(is.finite(suppressWarnings(as.numeric(person_df$SE))))) {
        person_attach_cols <- setdiff(person_attach_cols, "SE")
        person_m$SE <- NULL
      }
      person_df[, intersect(person_attach_cols, names(person_df))] <- NULL
      person_df <- merge(person_df, person_m, by = "Person",
                         all.x = TRUE, sort = FALSE)
      fit$facets$person <- tibble::as_tibble(person_df)
    }
  }

  step_unc <- as.data.frame(structural_uncertainty$steps %||% data.frame(), stringsAsFactors = FALSE)
  if (!is.null(fit$steps) && nrow(fit$steps) > 0L && nrow(step_unc) > 0L) {
    step_keys <- intersect(c("StepFacet", "Step"), names(fit$steps))
    if (!"StepFacet" %in% names(fit$steps) && "Step" %in% names(fit$steps)) {
      step_keys <- "Step"
    }
    if (length(step_keys) > 0L && all(step_keys %in% names(step_unc))) {
      attach_step_cols <- setdiff(names(step_unc), c(step_keys, "Estimate"))
      if (length(attach_step_cols) > 0L) {
        step_df <- as.data.frame(fit$steps, stringsAsFactors = FALSE)
        step_df[, intersect(attach_step_cols, names(step_df))] <- NULL
        fit$steps <- merge(
          step_df,
          step_unc[, c(step_keys, attach_step_cols), drop = FALSE],
          by = step_keys,
          all.x = TRUE,
          sort = FALSE
        )
      }
    }
  }

  slope_unc <- as.data.frame(structural_uncertainty$slopes %||% data.frame(), stringsAsFactors = FALSE)
  if (!is.null(fit$slopes) && nrow(fit$slopes) > 0L && nrow(slope_unc) > 0L &&
      "SlopeFacet" %in% names(fit$slopes) && "SlopeFacet" %in% names(slope_unc)) {
    attach_slope_cols <- setdiff(names(slope_unc), c("SlopeFacet", "LogEstimate", "Estimate"))
    if (length(attach_slope_cols) > 0L) {
      slope_df <- as.data.frame(fit$slopes, stringsAsFactors = FALSE)
      slope_df[, intersect(attach_slope_cols, names(slope_df))] <- NULL
      fit$slopes <- merge(
        slope_df,
        slope_unc[, c("SlopeFacet", attach_slope_cols), drop = FALSE],
        by = "SlopeFacet",
        all.x = TRUE,
        sort = FALSE
      )
    }
  }

  if (!is.null(fit$readiness$parameters)) {
    fit$readiness$parameters <- mfrm_update_slope_readiness_parameters(
      fit$readiness$parameters, slope_unc)
  }
  fit$config$attached_diagnostics <- TRUE
  fit$config$attached_diagnostics_cols <- attached_cols
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
  convergence <- mfrm_convergence_state(fit)
  pop$estimation_converged <- isTRUE(convergence$code_converged)
  pop$inference_ready <- isTRUE(convergence$inference_ready)
  # Retain the historical field without allowing its ambiguous name to hide
  # the distinction between optimizer convergence and inference readiness.
  pop$converged <- pop$inference_ready
  pop$converged_basis <- "legacy_alias_of_inference_ready"
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
                                             require_full_rank = TRUE,
                                             population_terms = NULL) {
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
      notes = "No population model was requested; this fit uses an unconditional normal person distribution."
    ))
  }

  population_policy <- match.arg(population_policy)
  if (inherits(population_formula, "terms")) {
    population_contrasts <- population_contrasts %||%
      attr(population_formula, "mfrm_population_contrasts", exact = TRUE)
    population_xlevels <- population_xlevels %||%
      attr(population_formula, "mfrm_population_xlevels", exact = TRUE)
  }
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
    formula = population_terms %||% population_formula,
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
         "Use `population_policy = 'omit'` to fit the complete-case analysis set instead. Affected IDs: ",
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
  # Rank-deficiency alone will not catch near-singular designs (e.g. two
  # covariates with correlation 0.9999). Warn when the reciprocal condition
  # number is very small so users are not silently handed a fit whose
  # coefficients are dominated by numerical noise.
  mm_rcond <- tryCatch(suppressWarnings(rcond(mm)), error = function(e) NA_real_)
  if (isTRUE(require_full_rank) && is.finite(mm_rcond) && mm_rcond < 1e-8) {
    warning(sprintf(
      paste0("Latent-regression design matrix is near-singular (rcond = %.1e). ",
             "Estimated coefficients may be unstable; consider dropping ",
             "redundant covariates or rescaling `person_data`."),
      mm_rcond
    ), call. = FALSE)
  }

  notes <- c(
    "Latent-regression covariates are active.",
    "Person-level covariates and their design matrix were prepared for MML estimation."
  )
  if (length(omitted_persons) > 0) {
    notes <- c(
      notes,
      paste0("Complete-case preparation retained ", nrow(person_tbl), " person(s) and omitted ",
             length(omitted_persons), " person(s) under `population_policy = 'omit'`.")
    )
  }
  if (response_rows_omitted > 0) {
    notes <- c(
      notes,
      paste0("For latent-regression estimation, ", response_rows_omitted,
             " response row(s) would be excluded because their persons lack complete background data.")
    )
  }

  list(
    active = TRUE,
    formula = population_formula,
    terms = terms_obj,
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

format_anchor_review_message <- function(anchor_review) {
  if (is.null(anchor_review$issue_counts) || nrow(anchor_review$issue_counts) == 0) {
    return("Anchor review detected no issues.")
  }
  nonzero <- anchor_review$issue_counts |>
    dplyr::filter(.data$N > 0)
  if (nrow(nonzero) == 0) {
    return("Anchor review detected no issues.")
  }
  labels <- paste0(nonzero$Issue, "=", nonzero$N)
  paste0(
    "Anchor review detected ", sum(nonzero$N), " issue row(s): ",
    paste(labels, collapse = "; "),
    ". Invalid rows, if any, were removed; duplicate keys keep the last row.",
    if (any(nonzero$Issue == "overlap_anchor_group"))
      " Overlapping direct and group anchors retain both constraints." else ""
  )
}

normalize_expected_mfrm_design <- function(expected_design, person, facets) {
  if (is.null(expected_design)) return(NULL)
  if (!is.data.frame(expected_design)) {
    stop("`expected_design` must be a data.frame.", call. = FALSE)
  }
  required <- c(person, facets)
  missing_cols <- setdiff(required, names(expected_design))
  if (length(missing_cols) > 0L) {
    stop(
      "`expected_design` is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (nrow(expected_design) == 0L) {
    stop("`expected_design` must contain at least one planned rating cell.", call. = FALSE)
  }

  out <- as.data.frame(expected_design[, required, drop = FALSE], stringsAsFactors = FALSE)
  names(out)[1L] <- "Person"
  out[] <- lapply(out, function(x) trimws(as.character(x)))
  invalid <- vapply(out, function(x) any(is.na(x) | !nzchar(x)), logical(1))
  if (any(invalid)) {
    stop(
      "`expected_design` contains missing or blank IDs in: ",
      paste(names(out)[invalid], collapse = ", "), ".",
      call. = FALSE
    )
  }
  duplicate_rows <- duplicated(out) | duplicated(out, fromLast = TRUE)
  if (any(duplicate_rows)) {
    stop(
      "`expected_design` must contain one row per planned Person x ",
      paste(facets, collapse = " x "),
      " cell; duplicate planned cells were found.",
      call. = FALSE
    )
  }
  out
}

audit_mfrm_expected_design <- function(df, expected_design, facets) {
  key_cols <- c("Person", facets)
  observed <- unique(as.data.frame(df[, key_cols, drop = FALSE], stringsAsFactors = FALSE))
  empty_cells <- observed[0, , drop = FALSE]

  if (is.null(expected_design)) {
    level_coverage <- purrr::map_dfr(facets, function(facet) {
      tibble::tibble(
        Facet = facet,
        ExpectedLevels = NA_integer_,
        ObservedLevels = dplyr::n_distinct(observed[[facet]]),
        ExpectedOnlyLevels = NA_integer_,
        UnexpectedLevels = NA_integer_
      )
    })
    return(list(
      summary = data.frame(
        Status = "not_declared",
        ExpectedCells = NA_integer_,
        ObservedCells = nrow(observed),
        MatchedCells = NA_integer_,
        MissingExpectedCells = NA_integer_,
        UnexpectedObservedCells = NA_integer_,
        CoverageRate = NA_real_,
        ExpectedOnlyPersons = NA_integer_,
        UnexpectedPersons = NA_integer_,
        stringsAsFactors = FALSE
      ),
      missing_expected_cells = empty_cells,
      unexpected_observed_cells = empty_cells,
      level_coverage = as.data.frame(level_coverage, stringsAsFactors = FALSE),
      settings = list(
        declared = FALSE,
        key = key_cols,
        note = paste0(
          "Structural missingness was not assessed because `expected_design` was not supplied. ",
          "Absent rows cannot be distinguished from cells that were never assigned."
        )
      )
    ))
  }

  matched <- dplyr::semi_join(expected_design, observed, by = key_cols)
  missing_expected <- dplyr::anti_join(expected_design, observed, by = key_cols)
  unexpected_observed <- dplyr::anti_join(observed, expected_design, by = key_cols)
  level_coverage <- purrr::map_dfr(facets, function(facet) {
    expected_levels <- unique(expected_design[[facet]])
    observed_levels <- unique(observed[[facet]])
    tibble::tibble(
      Facet = facet,
      ExpectedLevels = length(expected_levels),
      ObservedLevels = length(observed_levels),
      ExpectedOnlyLevels = length(setdiff(expected_levels, observed_levels)),
      UnexpectedLevels = length(setdiff(observed_levels, expected_levels))
    )
  })

  list(
    summary = data.frame(
      Status = "declared",
      ExpectedCells = nrow(expected_design),
      ObservedCells = nrow(observed),
      MatchedCells = nrow(matched),
      MissingExpectedCells = nrow(missing_expected),
      UnexpectedObservedCells = nrow(unexpected_observed),
      CoverageRate = nrow(matched) / nrow(expected_design),
      ExpectedOnlyPersons = length(setdiff(expected_design$Person, observed$Person)),
      UnexpectedPersons = length(setdiff(observed$Person, expected_design$Person)),
      stringsAsFactors = FALSE
    ),
    missing_expected_cells = as.data.frame(missing_expected, stringsAsFactors = FALSE),
    unexpected_observed_cells = as.data.frame(unexpected_observed, stringsAsFactors = FALSE),
    level_coverage = as.data.frame(level_coverage, stringsAsFactors = FALSE),
    settings = list(
      declared = TRUE,
      key = key_cols,
      note = "Observed unique cells were compared with the explicitly declared assignment roster."
    )
  )
}

summarize_duplicate_mfrm_cells <- function(df, facets) {
  key_cols <- c("Person", facets)
  detail <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) |>
    dplyr::summarise(
      Observations = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$Observations > 1L) |>
    dplyr::arrange(dplyr::desc(.data$Observations)) |>
    as.data.frame(stringsAsFactors = FALSE)

  list(
    summary = data.frame(
      DuplicateCells = nrow(detail),
      RowsInDuplicateCells = if (nrow(detail) > 0L) sum(detail$Observations) else 0L,
      ExtraRows = if (nrow(detail) > 0L) sum(detail$Observations - 1L) else 0L,
      HasDuplicates = nrow(detail) > 0L,
      stringsAsFactors = FALSE
    ),
    detail = detail
  )
}

mfrm_bipartite_components <- function(cells, facet, basis, include_persons = FALSE) {
  edges <- unique(as.data.frame(cells[, c("Person", facet), drop = FALSE], stringsAsFactors = FALSE))
  edges$Person <- as.character(edges$Person)
  edges[[facet]] <- as.character(edges[[facet]])
  persons_to_levels <- split(edges[[facet]], edges$Person)
  levels_to_persons <- split(edges$Person, edges[[facet]])
  unseen <- unique(edges$Person)
  components <- list()

  while (length(unseen) > 0L) {
    component_persons <- unseen[1L]
    component_levels <- character(0)
    repeat {
      next_levels <- unique(unname(unlist(
        persons_to_levels[component_persons], use.names = FALSE
      )))
      next_persons <- unique(unname(unlist(
        levels_to_persons[next_levels], use.names = FALSE
      )))
      expanded_persons <- union(component_persons, next_persons)
      expanded_levels <- union(component_levels, next_levels)
      if (setequal(expanded_persons, component_persons) &&
          setequal(expanded_levels, component_levels)) break
      component_persons <- expanded_persons
      component_levels <- expanded_levels
    }
    component_edges <- sum(
      edges$Person %in% component_persons & edges[[facet]] %in% component_levels
    )
    components[[length(components) + 1L]] <- data.frame(
      Basis = basis,
      Facet = facet,
      Component = length(components) + 1L,
      Persons = length(component_persons),
      FacetLevels = length(component_levels),
      Edges = component_edges,
      LevelLabels = paste(sort(component_levels), collapse = ", "),
      PersonLabels = if (isTRUE(include_persons)) {
        paste(sort(component_persons), collapse = ", ")
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
    unseen <- setdiff(unseen, component_persons)
  }

  component_table <- do.call(rbind, components)
  largest_index <- which.max(component_table$Persons)
  summary <- data.frame(
    Basis = basis,
    Facet = facet,
    PersonNodes = length(unique(edges$Person)),
    FacetLevelNodes = length(unique(edges[[facet]])),
    Edges = nrow(edges),
    Components = nrow(component_table),
    LargestComponentPersons = component_table$Persons[largest_index],
    LargestComponentLevels = component_table$FacetLevels[largest_index],
    LargestComponentPercent = 100 * component_table$Persons[largest_index] /
      length(unique(edges$Person)),
    Connected = nrow(component_table) == 1L,
    stringsAsFactors = FALSE
  )
  list(summary = summary, components = component_table)
}

audit_mfrm_connectivity <- function(df, expected_design, facets, include_persons = FALSE) {
  bases <- list(observed = df)
  if (!is.null(expected_design)) bases$declared_expected <- expected_design
  audits <- unlist(lapply(names(bases), function(basis) {
    lapply(facets, function(facet) {
      mfrm_bipartite_components(
        bases[[basis]], facet, basis, include_persons = include_persons
      )
    })
  }), recursive = FALSE)

  list(
    summary = do.call(rbind, lapply(audits, `[[`, "summary")),
    components = do.call(rbind, lapply(audits, `[[`, "components"))
  )
}

summarize_linkage_by_facet <- function(df, facet, min_linking_persons = 2L) {
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
    Persons = nrow(by_person),
    MinPersonsPerLevel = min(by_level$PersonsPerLevel, na.rm = TRUE),
    MedianPersonsPerLevel = stats::median(by_level$PersonsPerLevel, na.rm = TRUE),
    MinLevelsPerPerson = min(by_person$LevelsPerPerson, na.rm = TRUE),
    MedianLevelsPerPerson = stats::median(by_person$LevelsPerPerson, na.rm = TRUE),
    LinkingPersons = sum(by_person$LevelsPerPerson > 1L),
    LinkingPersonRate = mean(by_person$LevelsPerPerson > 1L),
    SingleLevelPersons = sum(by_person$LevelsPerPerson == 1L),
    SparseLevels = sum(by_level$PersonsPerLevel < min_linking_persons),
    SparseLevelThreshold = min_linking_persons
  )
}

canonical_compare_interaction_terms <- function(x) {
  x <- as.character(x %||% character(0))
  if (length(x) == 0L) return(character(0))
  out <- vapply(strsplit(x, ":", fixed = TRUE), function(parts) {
    paste(sort(parts), collapse = ":")
  }, character(1))
  sort(unique(out))
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

  signature <- list(
    model = as.character(cfg$model %||% NA_character_),
    method = as.character(cfg$method %||% NA_character_),
    person = as.character(cfg$person_col %||% NA_character_),
    facets = sort(as.character(cfg$facet_cols %||% fit$prep$facet_names %||% character(0))),
    score = as.character(cfg$score_col %||% NA_character_),
    rating_min = suppressWarnings(as.numeric(
      cfg$rating_min %||% fit$prep$rating_min %||% NA_real_
    )),
    rating_max = suppressWarnings(as.numeric(
      cfg$rating_max %||% fit$prep$rating_max %||% NA_real_
    )),
    score_map = cfg$score_map %||% fit$prep$score_map %||% NULL,
    weight = as.character(cfg$weight_col %||% NA_character_),
    step_facet = as.character(cfg$step_facet %||% NA_character_),
    slope_facet = as.character(cfg$slope_facet %||% NA_character_),
    facet_interactions = canonical_compare_interaction_terms(
      cfg$facet_interactions
    ),
    noncenter_facet = as.character(cfg$noncenter_facet %||% "Person"),
    centering_scope = cfg$centering_scope %||% "remainder",
    dummy_facets = sort(as.character(cfg$dummy_facets %||% character(0))),
    positive_facets = sort(as.character(cfg$positive_facets %||% character(0))),
    anchors = anchor_tbl,
    group_anchors = group_tbl
  )
  # Preserve saved legacy identities when the constraint equations are unchanged.
  if (nrow(anchor_tbl) + nrow(group_tbl) == 0L ||
      identical(signature$centering_scope, "remainder")) {
    signature$centering_scope <- NULL
  }
  signature
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

mfrm_population_design <- function(fit) {
    pop <- fit$config$population_spec %||% list()
    if (!isTRUE(pop$active)) return(list(active = FALSE))
    persons <- as.character(fit$prep$levels$Person %||% character(0))
    design <- pop$design_matrix
    columns <- pop$design_columns
    lookup <- pop$person_lookup
    if (!is.matrix(design) || !is.numeric(design) || ncol(design) == 0L ||
        length(persons) == 0L || anyNA(persons) || anyDuplicated(persons) ||
        !is.character(columns) || length(columns) != ncol(design) ||
        anyNA(columns) || anyDuplicated(columns) ||
        !is.numeric(lookup) || length(lookup) != length(persons) ||
        any(!is.finite(lookup)) || any(lookup != floor(lookup)) ||
        any(lookup < 1L | lookup > nrow(design))) return(NULL)
    design <- design[lookup[order(persons, method = "radix")],
                     order(columns, method = "radix"), drop = FALSE]
    if (any(!is.finite(design))) return(NULL)
    out <- list(active = TRUE, persons = sort(persons, method = "radix"),
         columns = sort(columns, method = "radix"),
         design = matrix(as.numeric(design), nrow = length(persons)))
    if (!is.null(pop$fixed_sd)) out$fixed_sd <- pop$fixed_sd
    out
}

audit_compare_mfrm_nesting <- function(fits, labels) {
  if (length(fits) != 2L) {
    return(list(
      eligible = FALSE,
      reason = "Likelihood-ratio tests are only reviewed for two-model comparisons.",
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
    rating_min = same_signature_component(sigs[[1]]$rating_min, sigs[[2]]$rating_min),
    rating_max = same_signature_component(sigs[[1]]$rating_max, sigs[[2]]$rating_max),
    score_map = same_signature_component(sigs[[1]]$score_map, sigs[[2]]$score_map),
    weight = same_signature_component(sigs[[1]]$weight, sigs[[2]]$weight),
    noncenter_facet = same_signature_component(sigs[[1]]$noncenter_facet, sigs[[2]]$noncenter_facet),
    centering_scope = same_signature_component(sigs[[1]]$centering_scope, sigs[[2]]$centering_scope),
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

  # For the supported response-model restrictions, require a shared population
  # model. Compare the declared design at the same Persons; coefficients and
  # variance are reestimated under each hypothesis.
  population_designs <- lapply(fits, mfrm_population_design)
  if (any(vapply(population_designs, is.null, logical(1))) ||
      !identical(population_designs[[1]], population_designs[[2]])) {
    return(list(
      eligible = FALSE,
      reason = paste(
        "A shared population model could not be verified.",
        "Automatic nesting requires the same fixed population or the same",
        "estimated-normal population design (declared columns and values",
        "aligned by Person). Changed or unavailable designs require separate review."
      ),
      simpler = NA_character_,
      complex = NA_character_,
      relation = "population_model_unverified"
    ))
  }

  model_pair <- toupper(c(sigs[[1]]$model, sigs[[2]]$model))
  interaction_sets <- lapply(sigs, function(sig) {
    sort(unique(as.character(sig$facet_interactions %||% character(0))))
  })
  if (identical(model_pair[1], model_pair[2])) {
    same_family_components <- c(
      step_facet = same_signature_component(sigs[[1]]$step_facet, sigs[[2]]$step_facet),
      slope_facet = same_signature_component(sigs[[1]]$slope_facet, sigs[[2]]$slope_facet)
    )
    if (!all(same_family_components)) {
      mismatch <- names(same_family_components)[!same_family_components]
      return(list(
        eligible = FALSE,
        reason = paste0(
          "Same-family fits differ in structural model settings: ",
          paste(mismatch, collapse = ", "),
          "."
        ),
        simpler = NA_character_,
        complex = NA_character_,
        relation = "same_model"
      ))
    }

    int_1 <- interaction_sets[[1]]
    int_2 <- interaction_sets[[2]]
    first_in_second <- all(int_1 %in% int_2)
    second_in_first <- all(int_2 %in% int_1)
    if (length(int_1) < length(int_2) && first_in_second) {
      added <- setdiff(int_2, int_1)
      return(list(
        eligible = TRUE,
        reason = paste0(
          "Supported nesting review passed: shared model family and constraints; ",
          "the complex fit adds fixed facet interaction term(s): ",
          paste(added, collapse = ", "),
          "."
        ),
        simpler = lbls[1],
        complex = lbls[2],
        relation = "facet_interaction_extension"
      ))
    }
    if (length(int_2) < length(int_1) && second_in_first) {
      added <- setdiff(int_1, int_2)
      return(list(
        eligible = TRUE,
        reason = paste0(
          "Supported nesting review passed: shared model family and constraints; ",
          "the complex fit adds fixed facet interaction term(s): ",
          paste(added, collapse = ", "),
          "."
        ),
        simpler = lbls[2],
        complex = lbls[1],
        relation = "facet_interaction_extension"
      ))
    }
    if (!identical(int_1, int_2)) {
      return(list(
        eligible = FALSE,
        reason = "Both fits use the same model family, but their fixed facet-interaction sets are not nested.",
        simpler = NA_character_,
        complex = NA_character_,
        relation = "same_model"
      ))
    }

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
    missing_interactions <- setdiff(interaction_sets[[idx_simple]], interaction_sets[[idx_complex]])
    if (length(missing_interactions)) {
      return(list(
        eligible = FALSE,
        reason = paste0("RSM interaction terms are absent from the PCM: ",
          paste(missing_interactions, collapse = ", "),
          ". Freeing category steps does not preserve these interaction effects."),
        simpler = NA_character_, complex = NA_character_,
        relation = "RSM_PCM_interactions_not_nested"
      ))
    }
    if (!is.na(step_facet) && nzchar(step_facet)) {
      return(list(
        eligible = TRUE,
        reason = paste0(
          "Supported nesting review passed: shared design/constraints with RSM nested inside PCM on step facet '",
          step_facet,
          "'."
        ),
        simpler = lbls[idx_simple],
        complex = lbls[idx_complex],
        relation = "RSM_in_PCM"
      ))
    }
  }

  if (setequal(model_pair, c("PCM", "GPCM"))) {
    idx_pcm <- which(model_pair == "PCM")[1]
    idx_gpcm <- which(model_pair == "GPCM")[1]
    pcm_step <- sigs[[idx_pcm]]$step_facet
    gpcm_step <- sigs[[idx_gpcm]]$step_facet
    gpcm_slope <- sigs[[idx_gpcm]]$slope_facet
    aligned_owner <- !is.na(pcm_step) && nzchar(pcm_step) &&
      identical(pcm_step, gpcm_step) && !is.na(gpcm_slope) && nzchar(gpcm_slope)
    review <- list(
      eligible = FALSE,
      reason = paste(
        "PCM and GPCM must share an explicit step owner and GPCM must name its slope owner,",
        "so even the unit-slope response-kernel reduction is not established."
      ),
      simpler = if (aligned_owner) lbls[idx_pcm] else NA_character_,
      complex = if (aligned_owner) lbls[idx_gpcm] else NA_character_,
      relation = "PCM_GPCM_owner_mismatch"
    )
    if (!aligned_owner) return(review)
    review$relation <- "PCM_GPCM_restriction_unverified"
    pcm <- fits[[idx_pcm]]$config
    gpcm <- fits[[idx_gpcm]]$config
    if (!identical(pcm$method, "MML") || !identical(gpcm$method, "MML")) {
      review$reason <- "The PCM/GPCM chi-square approximation requires MML fits; JML is not supported."
      return(review)
    }
    # Only slopes may differ: equal parameter counts alone cannot establish
    # matching step anchors, interactions or facet constraints.
    shared <- c("facet_names", "facet_levels", "facet_signs", "facet_specs", "step_specs")
    if (!identical(interaction_sets[[1]], interaction_sets[[2]]) ||
        !all(vapply(shared, function(nm) {
          !is.null(pcm[[nm]]) && identical(pcm[[nm]], gpcm[[nm]])
        }, logical(1))) ||
        !identical(pcm$interaction_specs, gpcm$interaction_specs)) {
      review$reason <- "PCM/GPCM testing requires identical facet constraints, step constraints and fixed interactions."
      return(review)
    }
    spec <- gpcm$gpcm_spec
    levels <- gpcm$facet_levels[[gpcm_slope]]
    sizes <- tryCatch(lapply(list(pcm, gpcm), build_param_sizes), error = function(e) NULL)
    restrictions <- length(levels) - 1L
    if (is.null(sizes) || restrictions < 1L ||
        !isTRUE(spec$active) || !identical(spec$levels, levels) ||
        !identical(spec$slope_facet, gpcm_slope) || !identical(spec$step_facet, gpcm_step) ||
        !identical(spec$identification, "sum_to_zero_log_slopes") ||
        !identical(spec$scale_reference, "geometric_mean_one") ||
        !identical(as.integer(sizes[[2]][["log_slopes"]]), restrictions) ||
        !identical(sizes[[1]], sizes[[2]][names(sizes[[2]]) != "log_slopes"])) {
      review$reason <- "The G-1 free relative-slope restrictions and otherwise identical free-parameter blocks could not be verified."
      return(review)
    }
    review$eligible <- TRUE
    review$relation <- "PCM_in_GPCM"
    review$df <- restrictions
    review$reason <- paste0(
      "PCM is the unit-slope response-kernel reduction on facet '", pcm_step,
      "' with slopes on facet '", gpcm_slope, "', a shared population model and ", restrictions,
      " independent log-slope restrictions. Unit slopes are interior; the",
      " reference is an ordinary asymptotic chi-square distribution. Numerical",
      " solution and likelihood-comparison checks are required separately."
    )
    return(review)
  }

  list(
    eligible = FALSE,
    reason = paste0(
      "Automatic nesting review supports RSM nested inside PCM, aligned PCM inside GPCM, ",
      "or same-family fixed facet-interaction extensions under shared design ",
      "and constraints."
    ),
    simpler = NA_character_,
    complex = NA_character_,
    relation = "unsupported"
  )
}

resolve_mfrm_category_policy <- function(keep_original, category_policy,
                                         keep_original_supplied) {
  if (!is.logical(keep_original) || length(keep_original) != 1L || is.na(keep_original)) {
    stop("`keep_original` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  if (is.null(category_policy)) return(keep_original)
  if (!is.character(category_policy) || length(category_policy) != 1L || is.na(category_policy)) {
    stop("`category_policy` must be NULL, 'collapse' or 'preserve'.", call. = FALSE)
  }
  category_policy <- match.arg(category_policy, c("collapse", "preserve"))
  preserve <- identical(category_policy, "preserve")
  if (keep_original_supplied && keep_original != preserve) {
    stop("`category_policy` and `keep_original` conflict; supply one choice or matching values.",
         call. = FALSE)
  }
  preserve
}

#' Check rating data before fitting an MFRM
#'
#' Inspect how many rating rows can be used, how often each score category
#' occurs, and whether facet levels connect through shared persons. This
#' function prepares descriptive checks; it does not fit a model or change
#' the data object supplied by the caller. Each row should represent one
#' rating, and `person`, `facets`, and `score` name its columns.
#'
#' @param data A data.frame in long format (one row per rating event).
#'   Grouped or rowwise tibbles are treated as one table; only the selected
#'   model columns enter row filtering and the design review.
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names. The reserved data
#'   and parameter names listed in [fit_mfrm()] cannot name facets; rename
#'   such columns before analysis.
#' @param score Column name for observed score. The category-code and
#'   near-integer normalization rules are the same as in [fit_mfrm()].
#' @param weight Optional weight/frequency column name, with the numeric
#'   precision and finite-sum rules of [fit_mfrm()].
#' @param rating_min Optional minimum category value. Supply with
#'   `rating_max` to retain unused boundary categories in the intended score
#'   support.
#' @param rating_max Optional maximum category value. Supply with
#'   `rating_min` to retain unused boundary categories in the intended score
#'   support.
#' @param keep_original Keep original category values. Use this with
#'   `rating_min` / `rating_max` when the intended scale has unused
#'   intermediate categories such as `1, 2, 4, 5` on a 1-5 scale.
#'   New code can instead use `category_policy = "preserve"`.
#' @inheritParams fit_mfrm
#' @param missing_codes Optional. `NULL` (default) is a no-op;
#'   `TRUE` or `"default"` activates the FACETS / SPSS / SAS
#'   convention (`c("99", "999", "-1", "N", "NA", "n/a", ".", "")`)
#'   for the score column while preserving person/facet identifiers;
#'   supply a character vector to apply a custom code set across all model
#'   columns. Replacement
#'   counts are returned in the `missing_recoding` component when
#'   supported by the calling helper. See [recode_missing_codes()]
#'   for the standalone version.
#' @param include_person_facet If `TRUE`, include person-level rows in
#'   `facet_level_summary`.
#' @param include_agreement If `TRUE`, include an observed-score agreement
#'   bundle (summary/pairs/settings) for a selected non-person facet.
#' @param rater_facet Optional facet name used to identify repeated scorers for
#'   agreement summaries. If `NULL`, a rater-like name such as `Rater`, `Judge`,
#'   or `Scorer` is inferred. No agreement analysis is run when such a name is
#'   absent; set this argument explicitly only when another facet genuinely
#'   represents repeated scorers.
#' @param context_facets Optional facets used to define matched contexts for
#'   agreement. If `NULL`, all remaining facets (including `Person`) are used.
#' @param agreement_top_n Optional maximum number of agreement pair rows.
#' @param expected_design Optional data frame declaring the planned assignment
#'   roster. It must contain the columns named by `person` and `facets`, with
#'   one row per planned Person x facet cell. Extra columns are ignored. When
#'   supplied, observed cells are compared with the roster so planned omissions
#'   can be distinguished from cells that were never assigned.
#' @param min_linking_persons Positive integer used as a descriptive sparse-link
#'   flag. A facet level observed for fewer than this many distinct persons is
#'   counted in `linkage_summary$SparseLevels`. This is a review threshold, not
#'   a model-acceptance rule.
#' @param comparisons Optional named specification for a prospective
#'   facet-by-origin interaction review. See "Comparison-specific review".
#'   Default `NULL` performs no model/target review or fitting.
#'
#' @section Comparison-specific review:
#' Use `comparisons = list(facet = ..., group = ..., groups = c("A", "B"),
#' contrasts = W, model = "RSM", method = "JML")` to review explicitly
#' declared interaction comparisons before estimation. Both `facet` and
#' `group` must occur in `facets`; each Person/recording must keep one origin
#' group. `groups` gives the subtraction order, including when additional
#' neutral origins are present. Unit weights and preserved categories are
#' required. `W` is a numeric matrix with unique comparison row names and
#' every reviewed facet level as a named column, including zero weights.
#' Each nonzero row must sum to zero. For four raters, weights
#' `c(.75, -.25, -.25, -.25)` compare the first rater's A-minus-B
#' interaction with the equal-weight panel; `c(1, -1, 0, 0)` compares the
#' first two raters' A-minus-B interactions. These are different targets.
#'
#' The prospective model has centered measurement facets, a zero-valued
#' group main facet, and the named facet-by-group interaction with zero
#' marginal sums. JML has free Person locations. Choose `model = "PCM"`
#' with an explicit `step_facet` for facet-owned steps. Optional
#' `positive_facets` records the location orientation. Anchors, extra
#' interactions and nonunit weights are outside this review. GPCM interaction
#' requests return `unsupported`, not an additive-model substitute.
#'
#' `method = "MML"` uses a fixed common standard-normal population unless
#' `population_formula` is supplied (an intercept or the origin group only).
#' `ConditionalState` refers to the linear response map; `FreePersonState`
#' reviews the same response structure with unrestricted Persons. A same-law
#' alias remains unidentified after mixing. A positive conditional rank alone
#' does not certify the marginal model: those MML rows have
#' `State = "not_evaluated"`. No fit or numerical marginal-rank calculation
#' is performed. JML rows distinguish `all_parameters_identified`,
#' `target_identified_only`, and `target_unidentified`. Numerically sensitive
#' or oversized dense maps remain `not_evaluated`.
#'
#' `IdentificationBasis` interprets these checks for each requested target:
#' `free_person_response_model` means the declared response map identifies
#' it even with unrestricted Person locations. `population_restriction` means
#' an unrestricted-Person alias is removed only after imposing the declared
#' population restrictions in the conditional map. This does **not** certify
#' marginal identification, and the corresponding MML `State` stays
#' `not_evaluated`. `conditional_support` identifies a CML target from response
#' patterns at the observed totals. `unidentified`, `not_evaluated` and
#' `unsupported` remain separate. `IdentificationDetail` explains the scope;
#' `Basis` retains the broader estimator/model basis.
#'
#' A fixed common N(0,1) population or estimated common normal population
#' (`~1`) restricts origin groups to the same ability distribution. This is
#' substantive, not arbitrary scale coding. Origin regression relaxes only
#' the means allowed by its design and retains a common residual SD. Interpret
#' the zero origin main effect jointly with those restrictions. A contrast
#' identified in a free-Person response map does not make its MML estimator
#' distribution-free. All conclusions still depend on the declared response
#' model, including equal slopes and shared/facet-owned step restrictions.
#'
#' `method = "CML"` instead reviews conditional support at the actual Person
#' totals and rating rosters. A dynamic program spans differences between all
#' possible patterns at those totals, so extreme scores and unavailable category
#' combinations do not manufacture information. `State` and `ConditionalState`
#' distinguish full identification, target-only identification and unidentified
#' targets; `FreePersonState` is the separate unrestricted-Person response-map
#' review. Person-constant main effects disappear, and no population distribution
#' or SD is estimated. Population formulas are not accepted for CML.
#'
#' Read `$comparison_review$table`, `$weights`, `$cells` and `$notes` together.
#' The table preserves group direction, model, step owner and counts; counts
#' cover nonzero-weight facet levels, while `cells` includes empty cells.
#' `$fit_arguments` records the prospective structural settings; use matching
#' input preparation when fitting. An identified target does not permit fitting
#' a rejected full parameterization. `Estimate` stays missing and formal
#' inference is unavailable. Save the review with its retained response rows.
#' Pass it as `comparison_review` to [analyze_dff()] to preserve this separate
#' prospective question alongside residual/refit effects; it does not certify
#' those different effects. See `vignette("mfrmr-linking-and-dff")` for a
#' school-assignment example and [mfrmr_linking_and_dff] for interpretation.
#' A matching, fully supported joint fit can return descriptive contrasts via
#' [interaction_effect_table()] or `summary(fit, comparison_review = review)`.
#' These explicit calls add retained-point numerical checks, not partial
#' estimation or formal uncertainty. To fit selected CML contrasts and obtain
#' conditional-information Wald intervals, use `fit_mfrm(fit_scope = "comparisons",
#' comparison_review = review)` with the CML review's `$fit_arguments`; see
#' the comparison-only section of [fit_mfrm()].
#'
#' @details
#' Set `rating_min` and `rating_max` from the rubric, including categories
#' nobody received. Use `keep_original = TRUE` to preserve its category
#' structure in the review. Numeric descriptives of score and weight use
#' `psych::describe()`.
#'
#' **Key data-quality checks to perform before fitting:**
#' - *Sparse categories*: review categories with little weighted support because
#'   their threshold estimates may be imprecise. Do not collapse categories
#'   solely from a package warning; also consider the rubric, intended score
#'   interpretation, and category diagnostics after fitting.
#' - *Unlinked elements*: inspect `design_connectivity` for the observed
#'   Person-facet graph. More than one component means that the levels of that
#'   facet are not connected through shared persons. This facet-specific check
#'   is conservative and does not by itself prove full model identification.
#' - *Extreme scores*: MML uses a person distribution to obtain posterior
#'   person scores, including for persons with all-minimum or all-maximum
#'   scores. Non-person facets remain fixed effects: an extreme rater or
#'   criterion is not given a prior or automatically shrunk by choosing MML.
#'   Review the fitted boundary and precision evidence before interpretation.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `overview`: confirms retained ratings (`Observations`), persons, facets,
#'   and category span. Use `row_retention` to compare input and retained
#'   `Rows`; `DroppedRows` counts exclusions during preparation.
#' - `missing_by_column`: counts `NA` values in the input columns. A missing
#'   score or required identifier excludes that rating row, not automatically
#'   the person's other ratings. The package does not fill missing ratings.
#'   When `missing_codes` is supplied, these counts still describe the original
#'   input; inspect `missing_recoding` and `preparation_notes` as well.
#' - `structural_missingness`: compares observed rating cells with
#'   `expected_design`, when supplied. Without a declared roster, structural
#'   missingness is reported as not assessed rather than assumed to be zero.
#' - `score_distribution`: checks sparse/unused score categories.
#'   Skew can be substantively expected, but weakly supported or unused
#'   categories need explicit interpretation.
#' - `facet_level_summary` and `linkage_summary`: checks per-level support,
#'   shared-person counts, and sparse levels. Use `design_connectivity` for the
#'   separate graph-component result.
#' - `agreement`: optional observed agreement summary for the selected scorer
#'   facet (exact agreement, correlation, and mean differences per pair).
#'
#' `data_review <- describe_mfrm_data(...)` saves all these checks.
#' `review <- summary(data_review)` provides a compact view; its missingness
#' table is named `review$missing`, while the original full table is
#' `data_review$missing_by_column`. Summary previews use `top_n = 10` by
#' default. Use the original tables to inspect all rows or categories.
#'
#' @section If the input needs attention:
#' - **Column name not found:** run `names(ratings)` and match spelling,
#'   spaces, and case in `person`, `facets`, and `score`.
#' - **Unexpected row loss:** inspect `row_retention`, `missing_by_column`,
#'   and `preparation_notes`. Resolve unintended missing IDs and invalid score
#'   text in the input data. For documented score markers such as `99` or `.`,
#'   use [recode_missing_codes()] with explicit `columns` and `codes`.
#' - **Repeated person-by-facet cells:** inspect `duplicate_cell_detail`.
#'   Correct accidental duplicates; declare `repeated_ratings = "independent"`
#'   for distinct conditionally independent events, or model dependence.
#'   Include a task or occasion facet when
#'   ratings represent distinct events in the design.
#' - **Unused category or disconnected design:** inspect
#'   `score_distribution` and `design_connectivity`. Review the rubric and
#'   assignments before changing the model. A retained internal zero-count
#'   category stops fitting; extra optimizer iterations cannot supply the
#'   missing category information.
#'
#' After editing or recoding ratings, rerun this review on the corrected data
#' and pass that same data to [fit_mfrm()], using the same columns, score
#' bounds, and `keep_original` setting. If you use `missing_codes` within the
#' review instead of recoding first, supply the same option to the fit:
#' reviewing does not modify the original data. For CSV import, column mapping,
#' and a complete worked example, see
#' `vignette("mfrmr-workflow", package = "mfrmr")`.
#'
#' @section Typical workflow:
#' 1. Run `data_review <- describe_mfrm_data(...)` on the rating data.
#' 2. Inspect row retention, category counts, and design connectivity.
#' 3. Correct input issues, repeat the review, and fit the reviewed data with
#'    [fit_mfrm()].
#'
#' @return A list of class `mfrm_data_description` with:
#' - `overview`: one-row run-level summary
#'   including `CategoryPolicy` and `ScoreRecoded`. The former records the
#'   selected category handling; the latter indicates whether original score
#'   values actually changed. A `"collapse"` policy can leave a contiguous
#'   scale unchanged. Inspect `score_support$score_map` for the mapping.
#' - `missing_by_column`: missing counts in selected input columns
#' - `missing_rate_summary`: per-column missingness rate summary
#'   (one row per input column, with raw and proportion-of-N columns)
#' - `score_descriptives`: output from [psych::describe()] for score
#' - `weight_descriptives`: output from [psych::describe()] for weight
#' - `score_distribution`: weighted and raw score frequencies over the prepared
#'   score support. Unused boundary categories are retained when the rating
#'   range was supplied explicitly; unused intermediate categories require
#'   `keep_original = TRUE`.
#' - `facet_level_summary`: per-level usage and score summaries
#' - `facet_crosstabs`: pairwise observation-count crosstabs between
#'   non-person facets (named list keyed `"facetA__facetB"`) for optional
#'   downstream coverage displays
#' - `linkage_summary`: person-facet connectivity diagnostics
#' - `structural_missingness`: declared-design comparison bundle containing a
#'   one-row summary, missing expected cells, unexpected observed cells,
#'   per-facet level coverage, and settings
#' - `design_connectivity`: component counts for each observed Person-facet
#'   graph and, when declared, each expected Person-facet graph
#' - `design_components`: component-level counts and facet-level labels;
#'   person labels are included only when `include_person_facet = TRUE`
#' - `duplicate_cell_summary`: counts of duplicate Person x facet cells
#' - `duplicate_cell_detail`: duplicate-cell keys and row counts
#' - `agreement`: observed-score agreement bundle for the selected scorer facet
#' - `row_retention`: row counts before and after preparation filters
#' - `preparation_notes`: structured notes for row drops, ID trimming, and
#'   design conditions detected during preparation
#' - `missing_recoding`: per-column counts of declared missing-code values
#'   replaced with `NA` before row filtering
#' - `score_support`: minimal prepared score-support metadata used by
#'   `summary(ds)$caveats`
#' - `comparison_review`: optional prospective interaction review, including
#'   all requested rows and their explicit weights; also retained by `summary()`.
#'
#' @seealso [fit_mfrm()], [summary.mfrm_data_description()],
#'   [recode_missing_codes()], [mfrmr_workflow_methods], [review_mfrm_anchors()]
#' @examples
#' library(mfrmr)
#' ratings <- load_mfrmr_data("example_operational")
#' head(ratings)
#'
#' # Check the data before fitting; the intended score categories are 1 to 4
#' data_review <- describe_mfrm_data(
#'   data = ratings,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   rating_min = 1,
#'   rating_max = 4,
#'   category_policy = "preserve"
#' )
#' data_review$row_retention       # Input and retained rows; check DroppedRows
#' data_review$missing_by_column   # Missing input values in each model column
#' data_review$score_distribution  # RawN is the number of ratings per category
#' data_review$design_connectivity # Components = 1 means connected for that facet
#' # Here all 282 rows are retained, and all four categories have observations
#'
#' # Save a compact summary when you want the overview and review notes
#' review <- summary(data_review)
#' review$overview
#' review$notes
#' # For the next fit, use data = ratings; data_review is a set of checks, not ratings
#' @export
#' @order 1
describe_mfrm_data <- function(data,
                               person,
                               facets,
                               score,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               missing_codes = NULL,
                               include_person_facet = FALSE,
                               include_agreement = TRUE,
                               rater_facet = NULL,
                               context_facets = NULL,
                               agreement_top_n = NULL,
                               expected_design = NULL,
                               min_linking_persons = 2L,
                               category_policy = NULL,
                               repeated_ratings = c("review", "independent"),
                               comparisons = NULL) {
  keep_original <- resolve_mfrm_category_policy(
    keep_original, category_policy, !missing(keep_original))
  if (!is.numeric(min_linking_persons) || length(min_linking_persons) != 1L ||
      !is.finite(min_linking_persons) || min_linking_persons < 1L ||
      abs(min_linking_persons - round(min_linking_persons)) > 1e-8) {
    stop("`min_linking_persons` must be a positive integer.", call. = FALSE)
  }
  min_linking_persons <- as.integer(round(min_linking_persons))
  expected_design_normalized <- normalize_expected_mfrm_design(
    expected_design = expected_design,
    person = person,
    facets = facets
  )
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    missing_codes = missing_codes,
    repeated_ratings = repeated_ratings,
    allow_single_person = TRUE
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
    purrr::map_dfr(prep$facet_names, function(facet) {
      summarize_linkage_by_facet(
        df,
        facet,
        min_linking_persons = min_linking_persons
      )
    })
  }
  structural_missingness <- audit_mfrm_expected_design(
    df = df,
    expected_design = expected_design_normalized,
    facets = prep$facet_names
  )
  connectivity <- audit_mfrm_connectivity(
    df = df,
    expected_design = expected_design_normalized,
    facets = prep$facet_names,
    include_persons = include_person_facet
  )
  duplicate_cells <- summarize_duplicate_mfrm_cells(df, prep$facet_names)

  agreement_bundle <- list(
    summary = data.frame(),
    pairs = data.frame(),
    settings = list(
      included = FALSE,
      status = if (isTRUE(include_agreement)) "not_computed" else "not_requested",
      rater_facet = NA_character_,
      context_facets = character(0),
      expected_exact_from_model = FALSE,
      top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n)),
      note = if (isTRUE(include_agreement)) {
        "No rater-like facet was detected. Set `rater_facet` explicitly to request facet-agreement summaries."
      } else {
        "Agreement summaries were not requested."
      }
    )
  )

  if (isTRUE(include_agreement) && length(prep$facet_names) > 0) {
    known_facets <- c("Person", prep$facet_names)
    if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
      rater_facet <- infer_default_rater_facet(
        prep$facet_names,
        fallback_first = FALSE
      )
    } else {
      rater_facet <- as.character(rater_facet[1])
    }
    if (is.null(rater_facet)) {
      rater_facet <- NULL
    } else if (!rater_facet %in% known_facets) {
      stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
    }
    if (!is.null(rater_facet) && identical(rater_facet, "Person")) {
      stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
    }

    if (!is.null(rater_facet) && is.null(context_facets)) {
      facet_cols <- known_facets
      resolved_context <- setdiff(facet_cols, rater_facet)
    } else if (!is.null(rater_facet)) {
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

    if (!is.null(rater_facet)) {
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
          status = "computed",
          rater_facet = rater_facet,
          context_facets = resolved_context,
          expected_exact_from_model = FALSE,
          top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n)),
          note = "Observed-score agreement was computed for the selected scorer facet."
        )
      )
    }
  }

  overview <- tibble::tibble(
    Observations = nrow(df),
    TotalWeight = total_weight,
    Persons = length(prep$levels$Person),
    Facets = length(prep$facet_names),
    Categories = prep$rating_max - prep$rating_min + 1,
    RatingMin = prep$rating_min,
    RatingMax = prep$rating_max,
    RatingRangeSource = prep$rating_range_source %||% "unknown",
    RatingMinSource = prep$rating_min_source %||% "unknown",
    RatingMaxSource = prep$rating_max_source %||% "unknown",
    !!!mfrm_category_settings(prep)
  )

  # Cross-tabulations between facet pairs give the
  # raw per-cell observation count for each pair of facets, which is
  # a heatmap-ready structure for downstream dashboards.
  facet_crosstabs <- list()
  if (length(prep$facet_names) >= 2L) {
    for (i in seq_len(length(prep$facet_names) - 1L)) {
      for (j in seq(i + 1L, length(prep$facet_names))) {
        a <- prep$facet_names[i]; b <- prep$facet_names[j]
        ctab <- df |>
          dplyr::count(
            LevelA = as.character(.data[[a]]),
            LevelB = as.character(.data[[b]]),
            name = "N"
          ) |>
          dplyr::mutate(FacetA = a, FacetB = b)
        facet_crosstabs[[paste(a, b, sep = "__")]] <-
          as.data.frame(ctab, stringsAsFactors = FALSE)
      }
    }
  }

  # Missing-rate summary at both column and facet-cell level.
  n_total <- nrow(data)
  missing_rate_summary <- data.frame(
    Column = missing_by_column$Column,
    Missing = as.integer(missing_by_column$Missing),
    NonMissing = n_total - as.integer(missing_by_column$Missing),
    MissingRate = ifelse(
      rep(n_total > 0, nrow(missing_by_column)),
      missing_by_column$Missing / n_total,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )

  out <- list(
    overview = overview,
    comparison_review = mfrmr_build_comparison_review(prep, comparisons),
    missing_by_column = missing_by_column,
    missing_rate_summary = missing_rate_summary,
    score_descriptives = score_desc,
    weight_descriptives = weight_desc,
    score_distribution = score_distribution,
    facet_level_summary = facet_level_summary,
    facet_crosstabs = facet_crosstabs,
    linkage_summary = linkage_summary,
    structural_missingness = structural_missingness,
    design_connectivity = as.data.frame(connectivity$summary, stringsAsFactors = FALSE),
    design_components = as.data.frame(connectivity$components, stringsAsFactors = FALSE),
    repeated_ratings = prep$repeated_ratings,
    duplicate_cell_summary = duplicate_cells$summary,
    duplicate_cell_detail = duplicate_cells$detail,
    agreement = agreement_bundle,
    row_retention = as.data.frame(prep$row_retention %||% data.frame(), stringsAsFactors = FALSE),
    preparation_notes = as.data.frame(prep$preparation_notes %||% data.frame(), stringsAsFactors = FALSE),
    missing_recoding = if (is.data.frame(prep$missing_recoding)) {
      as.data.frame(prep$missing_recoding, stringsAsFactors = FALSE)
    } else {
      prep$missing_recoding
    },
    score_support = list(
      keep_original = prep$keep_original,
      data = data.frame(Score = sort(unique(df$Score))),
      rating_min = prep$rating_min,
      rating_max = prep$rating_max,
      rating_range_source = prep$rating_range_source %||% "unknown",
      rating_min_source = prep$rating_min_source %||% "unknown",
      rating_max_source = prep$rating_max_source %||% "unknown",
      unused_score_categories = prep$unused_score_categories,
      score_map = prep$score_map
    )
  )
  class(out) <- c("mfrm_data_description", class(out))
  out
}

#' Recode common missing-value sentinels to `NA`
#'
#' Convenience helper that replaces the standard non-`NA` missing-code
#' sentinels used in SPSS / SAS / FACETS exports (`99`, `999`, `-1`,
#' `"N"`, `"NA"`, `"n/a"`, `"."`, `""`) with `NA` across the columns
#' you select. It is useful before calling [fit_mfrm()] on data exported with
#' those conventions.
#' A sentinel is a value used to mean "missing" instead of an actual score.
#' Only recode values that your data documentation defines as missing. For
#' custom score markers, specify both `columns = "Score"` and `codes`; this
#' preserves a person or rater identifier such as `99`.
#'
#' @param data A data frame.
#' @param columns Character vector of column names to recode. Defaults
#'   to `NULL`, in which case all columns are scanned.
#' @param codes Character vector of code values to convert to `NA`.
#'   Defaults to the FACETS / SPSS / SAS conventions; override when
#'   your instrument uses different sentinels.
#' @param numeric_codes Logical; if `TRUE` (default), numeric columns
#'   are also compared against the numeric conversion of `codes`.
#' @param verbose Logical; if `TRUE`, emits a `message()` summary of
#'   per-column replacement counts.
#'
#' @return The input `data` with the specified missing sentinels
#'   replaced by `NA`. A `mfrm_missing_recoding` attribute records the
#'   per-column replacement counts for traceability logs.
#'
#' @details
#' Save the returned data, for example `cleaned <- recode_missing_codes(...)`.
#' The original object is unchanged. This helper replaces cells and retains
#' every row; later, [describe_mfrm_data()] and [fit_mfrm()] exclude rating rows
#' with missing scores or required identifiers. Inspect the replacement counts
#' with `attr(cleaned, "mfrm_missing_recoding")`, then review and fit `cleaned`.
#'
#' The defaults differ across entry points: this helper scans all columns when
#' `columns` is omitted. In `describe_mfrm_data()` and `fit_mfrm()`,
#' `missing_codes = TRUE` uses the conventional code set on the score column
#' only; an explicit `missing_codes` vector applies to the person, facet, and
#' score columns. Use this helper with an explicit score column when your custom code
#' could also be a legitimate identifier.
#'
#' @seealso [describe_mfrm_data()], [fit_mfrm()].
#'
#' @examples
#' library(mfrmr)
#'
#' # A small input example: 99 and . mean missing only in the Score column
#' ratings <- data.frame(
#'   Person = c("001", "001", "002", "002"),
#'   Rater = c("R1", "99", "R1", "99"),
#'   Score = c("3", "99", ".", "2")
#' )
#' cleaned <- recode_missing_codes(
#'   ratings,
#'   columns = "Score",
#'   codes = c("99", ".")
#' )
#' cleaned # Two scores become NA; rater ID 99 and all four rows remain
#' attr(cleaned, "mfrm_missing_recoding") # Score: Replaced = 2
#' # Use the returned cleaned data for subsequent data review and fitting
#' @export
recode_missing_codes <- function(data,
                                 columns = NULL,
                                 codes = c("99", "999", "-1", "N", "NA",
                                           "n/a", ".", ""),
                                 numeric_codes = TRUE,
                                 verbose = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (is.null(columns)) columns <- names(data)
  bad_cols <- setdiff(columns, names(data))
  if (length(bad_cols) > 0L) {
    stop("Unknown column(s): ", paste(bad_cols, collapse = ", "), call. = FALSE)
  }
  codes <- as.character(codes)
  num_codes <- if (isTRUE(numeric_codes)) {
    suppressWarnings(as.numeric(codes))
  } else {
    numeric(0)
  }
  num_codes <- num_codes[is.finite(num_codes)]

  counts <- integer(length(columns))
  names(counts) <- columns
  for (col in columns) {
    x <- data[[col]]
    if (is.character(x)) {
      hit <- trimws(x) %in% codes
    } else if (is.numeric(x)) {
      hit <- x %in% num_codes
    } else {
      # factor or other: compare as character
      hit <- trimws(as.character(x)) %in% codes
    }
    n_hit <- sum(hit, na.rm = TRUE)
    counts[col] <- as.integer(n_hit)
    if (n_hit > 0L) {
      data[[col]][hit] <- NA
    }
  }
  attr(data, "mfrm_missing_recoding") <- data.frame(
    Column = columns,
    Replaced = as.integer(counts),
    stringsAsFactors = FALSE
  )
  if (isTRUE(verbose)) {
    total <- sum(counts, na.rm = TRUE)
    message(sprintf(
      "recode_missing_codes(): replaced %d cell(s) across %d column(s).",
      total, sum(counts > 0L, na.rm = TRUE)
    ))
  }
  data
}

#' @export
print.mfrm_data_description <- function(x, ...) {
  cat("mfrm data description\n")
  mfrmr_print_comparison_review(x$comparison_review)
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    print(x$overview, row.names = FALSE)
  }
  structural <- as.data.frame(
    x$structural_missingness$summary %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(structural) > 0L) {
    cat("\nPlanned assignment coverage\n")
    print(structural, row.names = FALSE)
  }
  if (!is.null(x$design_connectivity) && nrow(x$design_connectivity) > 0L) {
    cat("\nPerson-facet connectivity\n")
    print(as.data.frame(x$design_connectivity), row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(x$score_distribution, row.names = FALSE)
  }
  if (!is.null(x$agreement$summary) && nrow(x$agreement$summary) > 0) {
    agreement_facet <- as.character(x$agreement$settings$rater_facet %||% "selected facet")
    cat("\nObserved agreement by ", agreement_facet, "\n", sep = "")
    print(x$agreement$summary, row.names = FALSE)
  }
  invisible(x)
}

collect_mfrm_design_caveats <- function(object) {
  out <- empty_mfrm_caveats()
  add_caveat <- function(area, severity, condition, facets, message, action, details = "") {
    out <<- rbind(
      out,
      data.frame(
        Area = area,
        Severity = severity,
        Condition = condition,
        Categories = paste(as.character(facets), collapse = ", "),
        CategoryType = "design",
        Message = message,
        RecommendedAction = action,
        Details = details,
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }

  structural <- as.data.frame(
    object$structural_missingness$summary %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(structural) > 0L && identical(as.character(structural$Status[1L]), "declared")) {
    missing_n <- suppressWarnings(as.integer(structural$MissingExpectedCells[1L]))
    unexpected_n <- suppressWarnings(as.integer(structural$UnexpectedObservedCells[1L]))
    if (is.finite(missing_n) && missing_n > 0L) {
      add_caveat(
        area = "structural_missingness",
        severity = "review",
        condition = "planned_cells_not_observed",
        facets = object$structural_missingness$settings$key %||% character(0),
        message = paste0(
          missing_n,
          " planned rating cell(s) in `expected_design` were not observed."
        ),
        action = "Inspect the returned object's `structural_missingness$missing_expected_cells` table and document whether each absence was planned, unavailable, or lost."
      )
    }
    if (is.finite(unexpected_n) && unexpected_n > 0L) {
      add_caveat(
        area = "structural_missingness",
        severity = "warning",
        condition = "observed_cells_not_in_design",
        facets = object$structural_missingness$settings$key %||% character(0),
        message = paste0(
          unexpected_n,
          " observed rating cell(s) were not present in `expected_design`."
        ),
        action = "Reconcile the assignment roster and observed data before fitting."
      )
    }
  }

  connectivity <- as.data.frame(object$design_connectivity %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(connectivity) > 0L && all(c("Connected", "Basis", "Facet", "Components") %in% names(connectivity))) {
    disconnected <- connectivity[!as.logical(connectivity$Connected), , drop = FALSE]
    for (i in seq_len(nrow(disconnected))) {
      basis <- as.character(disconnected$Basis[i])
      facet <- as.character(disconnected$Facet[i])
      components <- suppressWarnings(as.integer(disconnected$Components[i]))
      add_caveat(
        area = "design_connectivity",
        severity = "warning",
        condition = paste0("disconnected_", basis, "_person_facet_graph"),
        facets = facet,
        message = paste0(
          "The ", gsub("_", " ", basis), " Person-", facet,
          " graph has ", components, " disconnected components."
        ),
        action = paste0(
          "Inspect the returned object's `design_components` table for ", facet,
          " and add defensible linking observations or analyze components separately."
        )
      )
    }
  }

  duplicates <- as.data.frame(object$duplicate_cell_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(duplicates) > 0L && isTRUE(duplicates$HasDuplicates[1L])) {
    add_caveat(
      area = "duplicate_cells",
      severity = if (identical(object$repeated_ratings, "independent")) "info" else "warning",
      condition = if (identical(object$repeated_ratings, "independent"))
        "independent_repeated_ratings_declared" else "duplicate_person_facet_cells",
      facets = object$structural_missingness$settings$key %||% character(0),
      message = paste0(
        duplicates$DuplicateCells[1L],
        " Person x facet cell(s) contain multiple observed rows."
      ),
      action = if (identical(object$repeated_ratings, "independent"))
        "Conditional independence was declared by the analyst; the data check does not verify it." else
        "Review repeated events: correct accidental copies, model dependence, or explicitly declare conditional independence."
    )
  }

  linkage <- as.data.frame(object$linkage_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(linkage) > 0L && all(c("Facet", "SparseLevels", "SparseLevelThreshold") %in% names(linkage))) {
    sparse <- linkage[suppressWarnings(as.integer(linkage$SparseLevels)) > 0L, , drop = FALSE]
    for (i in seq_len(nrow(sparse))) {
      add_caveat(
        area = "design_linkage",
        severity = "review",
        condition = "sparse_facet_levels",
        facets = sparse$Facet[i],
        message = paste0(
          sparse$SparseLevels[i], " level(s) of ", sparse$Facet[i],
          " were observed for fewer than ", sparse$SparseLevelThreshold[i],
          " distinct persons."
        ),
        action = "Review the affected level support and intended inference; the threshold is descriptive, not an automatic exclusion rule."
      )
    }
  }
  out
}

#' @rdname describe_mfrm_data
#' @details `summary()` method: Summarize a data-description object
#'
#' Read a compact summary of the checks from [describe_mfrm_data()] before
#' fitting a model. Save it with `review <- summary(data_review)` and select
#' the tables you need with `$`, as in the example.
#' @param object Output from [describe_mfrm_data()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in preview blocks.
#' @details For `summary()`: `data_review` holds the complete data checks; `review` holds summary tables
#' and notes. Neither object contains model estimates. The default `top_n = 10`
#' limits the missing-column and score-distribution previews; use
#' `data_review$missing_by_column` and `data_review$score_distribution` to see
#' the complete tables.
#' @section Interpreting output:
#' For `summary()`:
#'
#' Recommended read order:
#' - `overview`: retained ratings (`Observations`), persons, facets, and
#'   categories. Compare input and retained `Rows` in `row_retention` and
#'   investigate unexpected `DroppedRows`.
#'   `CategoryPolicy` and `ScoreRecoded` distinguish the selected policy from
#'   actual changes to score values. Older results without the policy record
#'   report `"not_recorded"`; an absent map gives `NA` for recoding.
#' - `missing`: input `NA` counts by column. This table is named
#'   `missing_by_column` in the original `data_review` object. Declared
#'   missing-code replacements and invalid score text can cause additional
#'   row loss; inspect `data_review$missing_recoding` and `preparation_notes`.
#' - `score_distribution`: category usage balance.
#' - `notes` / printed `Caveats`: retained zero-count score categories and
#'   related score-support caveats. With `keep_original = TRUE`, a retained
#'   unused internal category stops fitting; review the data and rubric first.
#' - `facet_overview`: coverage per facet (minimum/maximum weighted counts).
#' - `agreement`: observed-score agreement for the selected scorer facet (when
#'   available).
#' - `design_connectivity`: check for more than one observed component before
#'   comparing facet levels. `structural_missingness` reports planned omissions
#'   only when an assignment roster was supplied; `not_declared` does not mean
#'   that no ratings are missing.
#'
#' Very low `MinWeightedN` in `facet_overview` is a practical warning for
#' unstable downstream facet estimates.
#' @section Typical workflow:
#' For `summary()`:
#'
#' 1. Run [describe_mfrm_data()] on raw long-format data.
#' 2. Inspect `review <- summary(data_review)` before model fitting.
#' 3. Correct input issues and repeat the review, then pass the corrected
#'    rating data to [fit_mfrm()] with the same preparation settings.
#' @return For `summary()`: An object of class `summary.mfrm_data_description`.
#' - `overview`: design/sample counts
#' - `missing`: top columns by missingness
#' - `score_distribution`: compact score-usage table, including zero-count
#'   categories retained by the prepared score support
#' - `facet_overview`: facet-level coverage summary
#' - `structural_missingness`: declared assignment coverage summary; status is
#'   `"not_declared"` when no assignment roster was supplied
#' - `structural_level_coverage`: expected versus observed level counts
#' - `design_connectivity`: Person-facet component counts for observed and
#'   declared-expected designs
#' - `design_components`: component sizes and facet-level labels; person labels
#'   are suppressed unless explicitly requested in [describe_mfrm_data()]
#' - `linkage_summary`: sparse-support and shared-person counts by facet
#' - `duplicate_cell_summary`: aggregate duplicate-cell counts
#' - `agreement`: selected-facet agreement summary when available
#' - `agreement_settings`: selected scorer facet, matching context, and status
#' - `row_retention`: row counts before and after preparation filters
#' - `preparation_notes`: structured preparation notes retained from
#'   [describe_mfrm_data()]
#' - `reporting_map`: manuscript-oriented guide to what is covered here versus
#'   which companion outputs should be consulted
#' - `caveats`: structured warning/review rows for score-support issues;
#'   `print(summary(ds))` shows a compact `Caveats` block when rows are present
#' - `notes`: plain-language explanations of missingness, preparation, score
#'   support, and design-review findings
#' @seealso [describe_mfrm_data()], [summary.mfrm_fit()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("summary.mfrm_data_description", run.dontrun = TRUE)
#' \dontrun{
#' library(mfrmr)
#' ratings <- load_mfrmr_data("example_operational")
#'
#' # Demonstrate two missing scores in a copy of the example data
#' ratings$Score[1:2] <- NA
#' data_review <- describe_mfrm_data(
#'   data = ratings,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   rating_min = 1,
#'   rating_max = 4,
#'   category_policy = "preserve"
#' )
#' review <- summary(data_review)
#' review$row_retention # 282 input rows, 280 retained rows
#' review$missing       # Score has 2 missing input values
#' review$overview      # Counts describe the retained ratings
#' review$notes         # Explanations to read before fitting
#'
#' # The original description retains the full missingness table
#' data_review$missing_by_column
#' # Investigate missingness before using ratings in fit_mfrm()
#' }
#' @export
summary.mfrm_data_description <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(overview) > 0L) {
    category_settings <- mfrm_category_settings(object$score_support %||% object$prep)
    for (field in setdiff(names(category_settings), names(overview))) {
      overview[[field]] <- category_settings[[field]]
    }
  }
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
  row_retention <- as.data.frame(object$row_retention %||% data.frame(), stringsAsFactors = FALSE)
  preparation_notes <- as.data.frame(object$preparation_notes %||% data.frame(), stringsAsFactors = FALSE)
  structural_missingness <- as.data.frame(
    object$structural_missingness$summary %||% data.frame(),
    stringsAsFactors = FALSE
  )
  structural_level_coverage <- as.data.frame(
    object$structural_missingness$level_coverage %||% data.frame(),
    stringsAsFactors = FALSE
  )
  design_connectivity <- as.data.frame(object$design_connectivity %||% data.frame(), stringsAsFactors = FALSE)
  design_components <- as.data.frame(object$design_components %||% data.frame(), stringsAsFactors = FALSE)
  linkage_summary <- as.data.frame(object$linkage_summary %||% data.frame(), stringsAsFactors = FALSE)
  duplicate_cell_summary <- as.data.frame(object$duplicate_cell_summary %||% data.frame(), stringsAsFactors = FALSE)
  reporting_map <- data.frame(
    Area = c(
      "Sample / design counts",
      "Missingness review",
      "Planned assignment coverage",
      "Person-facet connectivity",
      "Score usage / category distribution",
      "Facet coverage",
      "Rater-facet agreement",
      "Fit / reliability / residual PCA"
    ),
    CoveredHere = c(
      "yes", "yes",
      if (nrow(structural_missingness) > 0L &&
          identical(as.character(structural_missingness$Status[1L]), "declared")) "yes" else "not declared",
      "yes", "yes", "yes",
      if (nrow(agreement_tbl) > 0) "yes" else "no", "no"
    ),
    CompanionOutput = c(
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...))",
      "data_review$structural_missingness$missing_expected_cells",
      "data_review$design_components",
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...))",
      "summary(describe_mfrm_data(...)) / plot(x, type = \"agreement\")",
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
  design_caveats <- collect_mfrm_design_caveats(object)
  caveats <- rbind(caveats, design_caveats)
  if (nrow(design_caveats) > 0L && "Message" %in% names(design_caveats)) {
    notes <- c(notes, as.character(design_caveats$Message))
  }
  structural_note <- as.character(object$structural_missingness$settings$note %||% "")
  if (nzchar(structural_note)) notes <- c(notes, structural_note)
  prep_review_messages <- if (nrow(preparation_notes) > 0L &&
                              all(c("Severity", "Message") %in% names(preparation_notes))) {
    preparation_notes$Message[
      tolower(as.character(preparation_notes$Severity)) %in% c("review", "warning", "error")
    ]
  } else {
    character(0)
  }
  notes <- c(notes, prep_review_messages)

  out <- list(
    overview = overview,
    comparison_review = object$comparison_review,
    missing = missing_tbl,
    score_distribution = score_dist,
    facet_overview = facet_overview,
    structural_missingness = structural_missingness,
    structural_level_coverage = structural_level_coverage,
    design_connectivity = design_connectivity,
    design_components = design_components,
    linkage_summary = linkage_summary,
    duplicate_cell_summary = duplicate_cell_summary,
    agreement = agreement_tbl,
    agreement_settings = object$agreement$settings %||% list(),
    row_retention = row_retention,
    preparation_notes = preparation_notes,
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
  mfrmr_print_comparison_review(x$comparison_review)
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$missing) && nrow(x$missing) > 0) {
    cat("\nMissing by column\n")
    print(round_numeric_df(as.data.frame(x$missing), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$structural_missingness) && nrow(x$structural_missingness) > 0) {
    cat("\nPlanned assignment coverage\n")
    print(round_numeric_df(as.data.frame(x$structural_missingness), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$design_connectivity) && nrow(x$design_connectivity) > 0) {
    cat("\nPerson-facet connectivity\n")
    print(round_numeric_df(as.data.frame(x$design_connectivity), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$linkage_summary) && nrow(x$linkage_summary) > 0) {
    cat("\nFacet linkage support\n")
    print(round_numeric_df(as.data.frame(x$linkage_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(round_numeric_df(as.data.frame(x$score_distribution), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_overview) && nrow(x$facet_overview) > 0) {
    cat("\nFacet coverage\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$row_retention) && nrow(x$row_retention) > 0 &&
      "DroppedRows" %in% names(x$row_retention) &&
      any(suppressWarnings(as.numeric(x$row_retention$DroppedRows)) > 0, na.rm = TRUE)) {
    cat("\nRow retention\n")
    print(round_numeric_df(as.data.frame(x$row_retention), digits = digits), row.names = FALSE)
  }
  print_preparation_section(x$preparation_notes)
  if (!is.null(x$agreement) && nrow(x$agreement) > 0) {
    agreement_facet <- as.character(x$agreement_settings$rater_facet %||% "selected facet")
    cat("\nObserved agreement by ", agreement_facet, "\n", sep = "")
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

#' @rdname describe_mfrm_data
#' @details `plot()` method: Plot a data-description object
#' @param x Output from [describe_mfrm_data()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"score_distribution"`, `"facet_levels"`, or `"missing"`.
#' @param main Optional title override.
#' @param palette Optional named colors (`score`, `facet`, `missing`).
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#' @details For `plot()`: This method draws quick pre-fit quality views from [describe_mfrm_data()]:
#' - score distribution balance
#' - facet-level structure size
#' - missingness by selected columns
#' @section Interpreting output:
#' For `plot()`:
#'
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
#' @section Typical workflow:
#' For `plot()`:
#'
#' 1. Run [describe_mfrm_data()] before fitting.
#' 2. Inspect `summary(ds)` and `plot(ds, type = "missing")`.
#' 3. Check category/facet balance with other plot types.
#' 4. Fit model after resolving obvious data issues.
#' @return For `plot()`: A plotting-data object of class `mfrm_plot_data`.
#' @seealso [describe_mfrm_data()], `plot()`
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_data_description", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(ds, draw = FALSE)
#' }
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

#' Review and normalize anchor/group-anchor tables
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param anchors Optional direct-anchor table (Facet, Level, Anchor). Retained
#'   values are fixed during estimation.
#' @param group_anchors Optional group-anchor table
#'   (Facet, Level, Group, GroupValue) defining group-mean constraints.
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param keep_original Keep original category values.
#'   New code can instead use `category_policy = "preserve"`.
#' @inheritParams fit_mfrm
#' @param missing_codes Optional. `NULL` (default) is a no-op;
#'   `TRUE` or `"default"` converts the FACETS / SPSS / SAS sentinel
#'   set to `NA` on the score column before review while preserving person and
#'   facet identifiers. Supply a character vector to apply a custom code set
#'   across the person, facet, and score columns.
#' @param min_common_anchors Minimum directly anchored levels per non-Person
#'   facet used in the package count recommendation (default `5`). The function
#'   cannot verify that those levels have invariant cross-run identity.
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in recommendations (default `30`).
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in recommendations (default `10`).
#' @param noncenter_facet One facet to leave non-centered. Other facets use
#'   `centering_scope`, as in [fit_mfrm()].
#' @inheritParams fit_mfrm
#' @param dummy_facets Facets to fix at zero.
#'
#' @details
#' This helper reviews computational constraints. A *direct anchor* fixes an
#' individual parameter to a supplied value. A *group anchor* constrains the
#' mean of declared elements to a supplied target. A *common-element link* is
#' different: it is observed overlap among administrations or subsets. Direct
#' or group constraints can transfer coordinates from a defensible reference,
#' but they do not manufacture empirical overlap.
#'
#' This function applies the same preprocessing and key-resolution rules
#' as `fit_mfrm()`, but returns a review object so constraints can be
#' checked *before* estimation. Running the review first helps avoid
#' estimation failures caused by misspecified or data-incompatible
#' anchors. It does not inspect source-fit readiness or establish that labels
#' denote the same invariant elements across runs. Those substantive checks
#' remain the caller's responsibility.
#' This review has no estimation-method argument. A clean table can still
#' request an unsupported model combination: in particular, `fit_mfrm()`
#' rejects Person direct/group anchors with MML. See the scale-constraint
#' section of [fit_mfrm()] before using partial anchors.
#'
#' **Anchor types:**
#' - *Direct anchors* fix individual element measures to specific logit
#'   values (e.g., Rater R01 anchored at 0.35 logits). Use internal coordinates,
#'   not values transformed for display with `Umean` or `Uscale`.
#' - *Group anchors* constrain the mean of a set of elements to a
#'   target value, allowing individual elements to vary freely around
#'   that mean. The target requires an external justification such as a
#'   defensible equal-mean or known-scale assumption.
#'   The mean is unweighted over declared levels, not over rating rows.
#'   `GroupValue` repeats one target per group. It is not a FACETS per-element
#'   contribution: values 0.2 and 1.0 in a native two-member group translate
#'   to 0.6 for both rows. See [fit_mfrm()] for group-origin centering and
#'   FACETS group 0 / extreme-exclusion differences.
#' - When both types include the same element, both constraints are retained:
#'   the direct anchor fixes that element and the group constraint still fixes
#'   the declared group mean. Incompatible combinations are rejected by the
#'   constraint/estimability checks.
#'
#' **Design checks** report whether each observed facet level has at least
#' `min_obs_per_element` weighted observations (default 30) and each
#' score category has at least `min_obs_per_category` (default 10).
#' These user-configurable counts are package screening recommendations, not
#' universal adequacy thresholds and not tests of connectedness or invariance.
#'
#' @section Interpreting output:
#' - `issue_counts`/`issues`: concrete data or specification problems.
#' - `facet_summary`: constraint coverage by facet. `FreeLevels` counts levels
#'   outside direct and group constraints, not free parameters. A group member
#'   may still vary, whereas centering removes a degree of freedom outside groups.
#' - `design_checks`: observation counts against the requested screening values.
#' - `recommendations`: action items before estimation.
#'
#' @section Typical workflow:
#' 1. Build candidate anchors (e.g., with [make_anchor_table()]).
#' 2. Run `review_mfrm_anchors(...)`.
#' 3. Resolve issues, then fit with [fit_mfrm()].
#'
#' @return A list of class `mfrm_anchor_review` with:
#' - `anchors`: cleaned direct constraints used by estimation
#' - `group_anchors`: cleaned group-mean constraints used by estimation
#' - `centering_scope`: the reviewed `"facet"`, `"remainder"` or
#'   `"anchor_origin"` convention
#' - `facet_summary`: counts of levels inside/outside declared constraints,
#'   `Centered` indicating the declared centering constraint, and `OriginBasis`
#'   describing how the facet's origin is set
#' - `design_checks`: observation-count checks by level/category
#' - `thresholds`: active threshold settings used for recommendations
#' - `issue_counts`: issue-type counts
#' - `issues`: list of issue tables
#' - `recommendations`: package-native anchor guidance strings
#'
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [make_anchor_table()],
#'   [mfrmr_linking_and_dff]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("review_mfrm_anchors", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_operational")
#' # Illustrative external value, not an estimate supplied with this dataset.
#' anchors <- data.frame(Facet = "Rater", Level = "R01", Anchor = 0.35)
#' review <- review_mfrm_anchors(
#'   data = ratings,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   anchors = anchors,
#'   rating_min = 1, rating_max = 4, category_policy = "preserve"
#' )
#' review$anchors                 # Check the retained value and matching ID.
#' review$issue_counts
#' review$recommendations         # Explains the selected centering convention.
#'
#' # A duplicate is a separate input problem: resolve it before fitting.
#' duplicate <- rbind(anchors, transform(anchors, Anchor = 0.10))
#' problem <- review_mfrm_anchors(
#'   ratings, "Person", c("Rater", "Criterion"), "Score",
#'   anchors = duplicate
#' )
#' problem$issues$duplicate_anchors
#' # Normalization retains the last value; it does not decide which is correct.
#' }
#' @name review_mfrm_anchors
#' @export
#' @order 1
review_mfrm_anchors <- function(data,
                                person,
                                facets,
                                score,
                                anchors = NULL,
                                group_anchors = NULL,
                                weight = NULL,
                                rating_min = NULL,
                                rating_max = NULL,
                                keep_original = FALSE,
                                missing_codes = NULL,
                                min_common_anchors = 5L,
                                min_obs_per_element = 30,
                                min_obs_per_category = 10,
                                noncenter_facet = "Person",
                                dummy_facets = NULL,
                                category_policy = NULL,
                                repeated_ratings = c("review", "independent"),
                                centering_scope = c("facet", "remainder", "anchor_origin")) {
  centering_scope <- match.arg(centering_scope)
  keep_original <- resolve_mfrm_category_policy(
    keep_original, category_policy, !missing(keep_original))
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    missing_codes = missing_codes,
    repeated_ratings = repeated_ratings,
    allow_single_person = TRUE
  )

  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)

  review <- audit_anchor_tables(
    prep = prep,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    centering_scope = centering_scope
  )
  # Check the same equalities as fitting, reusing the cleaned tables.
  prepare_constraint_specs(prep, noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets, centering_scope = centering_scope,
    anchor_review = review)
  class(review) <- unique(c("mfrm_anchor_review", class(review)))
  review
}

#' @export
print.mfrm_anchor_review <- function(x, ...) {
  issue_total <- if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) sum(x$issue_counts$N) else 0
  cat("mfrm anchor review\n")
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

#' @rdname review_mfrm_anchors
#' @details `summary()` method: Summarize an anchor-review object
#' @param object Output from [review_mfrm_anchors()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in issue previews.
#' @details For `summary()`: This summary provides a compact pre-estimation review of anchor and
#' group-anchor specifications.
#' @section Interpreting output:
#' For `summary()`:
#'
#' Recommended order:
#' - `issue_counts`: primary triage table (non-zero issues first).
#' - `facet_summary`: anchored/grouped/free-level balance by facet.
#' - `level_observation_summary` and `category_counts`: sparse-cell diagnostics.
#' - `recommendations`: concrete remediation suggestions.
#'
#' If `issue_counts` is non-empty, treat anchor constraints as provisional and
#' resolve issues before final estimation.
#' @section Typical workflow:
#' For `summary()`:
#'
#' 1. Run [review_mfrm_anchors()] with intended anchors/group anchors.
#' 2. Review `summary(review)` and recommendations.
#' 3. Revise anchor tables, then call [fit_mfrm()].
#' @return For `summary()`: An object of class `summary.mfrm_anchor_review`.
#' @seealso [review_mfrm_anchors()], [fit_mfrm()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("summary.mfrm_anchor_review", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' review <- review_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(review)
#' }
#' @export
summary.mfrm_anchor_review <- function(object, digits = 3, top_n = 10, ...) {
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
    "Anchor-review issues were detected. Review issue counts and recommendations."
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
  class(out) <- "summary.mfrm_anchor_review"
  out
}

#' @export
print.summary.mfrm_anchor_review <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Anchor Review Summary\n")
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

#' @rdname review_mfrm_anchors
#' @details `plot()` method: Plot an anchor-review object
#' @param x Output from [review_mfrm_anchors()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"issue_counts"`, `"facet_constraints"`,
#'   or `"level_observations"`.
#' @param main Optional title override.
#' @param palette Optional named colors.
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#' @details For `plot()`: Base-R visualization helper for anchor-review outputs.
#' @section Interpreting output:
#' For `plot()`:
#'
#' - `"issue_counts"`: volume of each issue class.
#' - `"facet_constraints"`: anchored/grouped/free mix by facet.
#' - `"level_observations"`: observation support across levels.
#' @section Typical workflow:
#' For `plot()`:
#'
#' 1. Run [review_mfrm_anchors()].
#' 2. Start with `plot(review, type = "issue_counts")`.
#' 3. Inspect constraint and support plots before fitting.
#' @return For `plot()`: A plotting-data object of class `mfrm_plot_data`.
#' @seealso [review_mfrm_anchors()], [make_anchor_table()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_anchor_review", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' review <- review_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(review, draw = FALSE)
#' }
#' @export
plot.mfrm_anchor_review <- function(x,
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
      stop("Issue-count table is not available. Ensure review_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    tbl <- tbl |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$N)),
        labels = as.character(tbl$Issue),
        col = pal["issues"],
        main = if (is.null(main)) "Anchor-review issue counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_review",
      list(plot = "issue_counts", table = tbl)
    )))
  }

  if (type == "facet_constraints") {
    tbl <- as.data.frame(x$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "AnchoredLevels", "GroupedLevels", "FreeLevels") %in% names(tbl))) {
      stop("Facet summary with constraint columns is not available. Ensure review_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
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
        main = ""
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
      "anchor_review",
      list(plot = "facet_constraints", table = tbl)
    )))
  }

  tbl <- as.data.frame(x$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Facet", "MinObsPerLevel") %in% names(tbl))) {
    stop("Level observation summary is not available. Ensure review_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
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
    "anchor_review",
    list(plot = "level_observations", table = tbl)
  ))
}

#' Build a candidate direct-anchor table from fitted estimates
#'
#' @param fit Output from [fit_mfrm()].
#' @param facets Optional subset of facets to include.
#' @param include_person Include person estimates as candidate anchors. Use only
#'   when cross-run person identity and the intended longitudinal constraint
#'   are substantively justified. Fixed Person values are supported only by
#'   a receiving JML fit; MML rejects Person anchors.
#' @param digits `NULL` (default) preserves the fitted numeric values. A
#'   nonnegative integer explicitly rounds them for display or external use.
#' @param readiness_policy How a source fit that is not inference-ready is
#'   handled. `"error"` (default) refuses anchor export or reuse. `"review"`
#'   permits explicit review-only extraction; those values must not be used as
#'   operational anchors.
#'
#' @details
#' This function performs a mechanical conversion from fitted estimates to the
#' `Facet`/`Level`/`Anchor` schema accepted by `fit_mfrm()`. The returned rows
#' are candidate direct constraints, not an approved anchor set.
#' `Anchor` uses the internal, untransformed location metric. The table does
#' not include category thresholds, GPCM slopes, model settings or source SEs;
#' it is not a complete saved calibration. See [extract_mfrm_calibration()] for that
#' separate workflow. Fixing an estimated value does not propagate its source
#' uncertainty into the new fit.
#'
#' By default, the function refuses export when the source fit is not
#' inference-ready under the current readiness contract. Set
#' `readiness_policy = "review"` only to inspect candidate values; this does not
#' make them eligible for reuse. The function cannot verify cross-run element
#' identity or invariance. Before reuse, document why selected elements retain
#' the same meaning, check the observed design's connectedness, and run
#' [review_mfrm_anchors()] on the receiving data. That review checks schema and
#' receiving-data support; it does not validate the source fit or the
#' substantive invariance assumption.
#'
#' Keep `digits = NULL` for recalibration. Rounding individual anchors can
#' break their joint sum or group-mean constraint, even at six decimal places.
#' Explicitly rounded values are different constraints and must be reviewed
#' before reuse; the fitter does not silently recenter supplied anchors.
#'
#' @section Interpreting output:
#' - `Facet`: facet name to be anchored in later runs.
#' - `Level`: specific element/level name inside that facet.
#' - `Anchor`: fixed logit value (rounded only when `digits` is supplied).
#'
#' @section Typical workflow:
#' 1. Fit and diagnose a defensible reference run with [fit_mfrm()].
#' 2. Confirm source readiness, element identity, and the intended link.
#' 3. Export candidates with `make_anchor_table(fit)` and review them.
#' 4. Pass selected rows back into `fit_mfrm(..., anchors = ...)`.
#'
#' @return A data.frame with `Facet`, `Level`, and `Anchor`.
#' @seealso [fit_mfrm()], [review_mfrm_anchors()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("make_anchor_table", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(
#'   ratings, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", rating_min = 1, rating_max = 4,
#'   category_policy = "preserve"
#' )
#' # Select a facet deliberately; matching labels alone do not prove invariance.
#' anchors_tbl <- make_anchor_table(fit, facets = "Rater")
#' anchors_tbl
#' }
#' @export
make_anchor_table <- function(fit,
                              facets = NULL,
                              include_person = FALSE,
                              digits = NULL,
                              readiness_policy = c("error", "review")) {
  stop_if_product_slopes(fit, "make_anchor_table()")
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }

  readiness_policy <- match.arg(readiness_policy)
  source_readiness <- mfrm_convergence_state(fit)
  fit_readiness <- as.character(
    source_readiness$fit_readiness %||% "legacy_unknown"
  )[1]
  if (is.na(fit_readiness) || !nzchar(fit_readiness)) {
    fit_readiness <- "legacy_unknown"
  }
  if (!isTRUE(source_readiness$inference_ready) &&
      identical(readiness_policy, "error")) {
    stop(
      "`fit` is not eligible for anchor export or reuse (FitReadiness = '",
      fit_readiness,
      "', InferenceReady = FALSE). Resolve the source fit or use ",
      "`readiness_policy = \"review\"` only to inspect review-only ",
      "candidate values.",
      call. = FALSE
    )
  }

  if (!is.null(digits) && (!is.numeric(digits) || is.complex(digits) ||
      length(digits) != 1L || !is.finite(digits) || digits < 0 || digits != floor(digits))) {
    stop("`digits` must be NULL or one nonnegative integer.", call. = FALSE)
  }
  anchor_values <- function(x) {
    if (is.null(digits)) as.numeric(x) else round(as.numeric(x), digits = digits)
  }
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
            Anchor = anchor_values(.data$Estimate)
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
          Anchor = anchor_values(.data$Estimate)
        )
    )
  }

  if (!is.null(facets)) {
    keep <- as.character(facets)
    out <- out |>
      dplyr::filter(.data$Facet %in% keep)
  }

  out <- out |>
    dplyr::arrange(.data$Facet, .data$Level)
  attr(out, "source_fit_readiness") <- fit_readiness
  attr(out, "source_inference_ready") <-
    isTRUE(source_readiness$inference_ready)
  attr(out, "readiness_policy") <- readiness_policy
  out
}

#' Compute diagnostics for an `mfrm_fit` object
#'
#' @param fit Output from [fit_mfrm()].
#' @param interaction_pairs Optional list of facet pairs.
#' @param top_n_interactions Number of top interactions.
#' @param whexact ZSTD transformation. `FALSE` (default) uses
#'   Wilson-Hilferty with the q/3 term (FACETS `WHEXACT=YES`);
#'   `TRUE` omits that term (FACETS `WHEXACT=NO`), where
#'   q = sqrt(2/df). The historical logical names are retained.
#'   `"legacy_linear"` explicitly requests the former TRUE calculation,
#'   (MnSq - 1) * sqrt(df/2), which is not FACETS WHEXACT=NO.
#' @param fit_df_method Degrees-of-freedom convention used for fit ZSTD.
#'   `"facets"` (default) uses the Wright-Masters fourth-moment approximation
#'   `df = 2 / q^2` and caps reported ZSTD at +/-9. `"engine"` retains
#'   `DF_Infit = sum(Var * Weight)` and `DF_Outfit = sum(Weight)`.
#'   `"both"` uses engine primary columns and adds `*_FACETS` companions.
#'   Engine transformations require df >= 1; the FACETS route admits df > 0.
#' @param diagnostic_mode Diagnostic basis to compute: `"both"` (the
#'   current default) computes both the residual/EAP-based stack and
#'   the strict latent-integrated first-order marginal-fit companion;
#'   `"legacy"` keeps the residual/EAP-based stack only;
#'   `"marginal_fit"` returns only the marginal-fit companion. The
#'   `"both"` path adds a posterior-integrated pass that typically
#'   doubles to quintuples wall-clock time relative to `"legacy"`;
#'   pass `"legacy"` explicitly when iterating on large designs and
#'   only the residual stack is needed. Use `"both"` for RSM/PCM
#'   reporting fits because it enables [plot(x, type = "marginal_fit")][plot.mfrm_fit] and
#'   [plot(x, type = "marginal_pairwise")][plot.mfrm_fit] follow-up.
#' @section Saved diagnostic definitions:
#' New diagnostics carry `calculation_version = 2L`. Recreate older diagnostics
#' with `diagnose_mfrm(fit)` before using current reports; no refit is needed.
#' For historical ZSTD use `fit_df_method = "engine"` and, if applicable,
#' `whexact = "legacy_linear"`. PTMEA, displacement, eligible chi-square df,
#' fair scores and unexpected-response defaults use the corrected definitions.
#'
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
#' **Legacy response calculations.** Saved MML EAPs are matched to responses
#' by Person identifier. JML retains optimizer-coordinate probabilities for
#' compatibility diagnostics. In `obs`, `PersonMeasure` records the ability
#' used, `PrimaryPersonMeasure` records the primary estimate (possibly
#' infinite), and `PersonMeasureBasis` identifies the calculation. A finite
#' JML optimizer value is not a finite primary estimate at a boundary.
#' For supported ordinary fits, [fitted()][mfrm_fit_responses] instead uses
#' primary JML limits. The posterior/corrected-JML targets have their own
#' documented routes in [mfrm_response_diagnostics()].
#'
#' `obs$ConditionalVariance` retains the actual conditional variance.
#' Legacy `Var`, `StdResidual` and `StdSq` are withheld when that variance is
#' at or below `1e-10`. A positive value below this cutoff is a diagnostic
#' screening convention, not a mathematically undefined Pearson residual;
#' [residuals(fit, "pearson")][mfrm_fit_responses] retains representable positive
#' variances within its supported scope. A zero variance makes the Pearson
#' residual undefined. Expected ratings and raw residuals can remain available
#' in either case. `PearsonAvailable`, `ResponseStatus` and `ResponseReason`
#' explain each row, and are retained in score/residual file tables.
#'
#' Overall and element-level mean squares require every contributing row's
#' standardized residual, variance and positive weight. If any is unavailable,
#' the aggregate and its degrees of freedom remain `NA`; missing residuals
#' are not counted as zero misfit or silently excluded. `FitResponseRows`,
#' `FitAvailableRows` and `FitUnavailableRows` count rows (not weight totals),
#' while `FitStatus` and `FitReason` explain availability. These fields accompany
#' `overall_fit`, `fit` and `measures`. FACETS-style degrees of freedom also
#' require complete fourth moments. Regenerate older diagnostics and exports
#' from the saved fit to obtain these explanations; model refitting is not
#' required.
#'
#' Marginal expected counts integrate over each Person's posterior conditioned
#' on the same observed responses, holding fitted calibration fixed. They are
#' not expectations from an independent replication or a prior-only population
#' margin. First-order residual scales use `sum(w^2 * p * (1-p))`; pairwise
#' scales use the analogous formula with products of row weights. These scales
#' omit cross-response/opportunity covariance and calibration-parameter
#' uncertainty. They are descriptive screens, not calibrated residual tests.
#'
#' Missing/invalid contributing probabilities, scores or weights withhold the
#' affected complete-scope aggregate rather than selecting usable rows. Missing
#' standardized residuals and flags remain `NA`. A known cutoff crossing stays
#' flagged when a companion rule is unavailable. `marginal_fit$coverage` records
#' classified/unclassified cells, groups and level pairs. Row/opportunity counts
#' remain visible; zero-weight opportunities contribute no information.
#' Available maxima are accompanied by classification/residual counts; they
#' are not maxima over unavailable values. RMSDs require the complete scope.
#' In PCM and GPCM, step-group summaries retain the declared `step_facet`.
#' Regenerate older marginal diagnostics and downstream summaries/plots/exports
#' from the existing fit and original diagnostic settings; no model refit is
#' required.
#'
#' Choosing `diagnostic_mode`:
#' - `"legacy"`: use when continuity with historical residual-based workflows is
#'   the priority.
#' - `"marginal_fit"`: use when you want the strict latent-integrated screen
#'   without the extra legacy bundle.
#' - `"both"`: recommended when you want continuity with the legacy residual
#'   stack while making the strict marginal path explicit for `RSM`, `PCM`,
#'   and `GPCM` fits.
#'
#' For `GPCM`, the same generalized partial credit kernel now
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
#' - **PTMEA**: within-element point-measure correlation between observed
#'   scores and the oriented sum of model measures (`obs$MeasureSum`),
#'   including Person rows. A positive value is directionally
#'   consistent with the fitted orientation; it is not a confirmatory test.
#'
#' The MnSq values and the ZSTD values should be read separately. `mfrmr` keeps
#' the fourth-moment df convention by default. The former engine convention
#' remains explicit. This changes some ZSTD screening decisions, not the
#' underlying MnSq. MML residual diagnostics remain EAP plug-in quantities;
#' this standardization does not make them posterior-integrated fit tests. Use
#' `fit_df_method = "both"` to review these two standardization conventions
#' side by side, using engine primary columns and FACETS companion columns.
#'
#' **Residual basis under MML.** For `method = "MML"` fits, residuals,
#' MnSq, and ZSTD are computed at the EAP person measures from the
#' marginal model. EAP measures are shrunken toward the population mean,
#' so expected scores -- and therefore fit statistics -- differ
#' systematically from JMLE-based engines such as FACETS, especially for
#' persons with extreme raw scores. The df conventions above do not remove
#' this difference: it is a residual-basis difference, not a
#' standardization difference. Refit with `method = "JML"` when an
#' external FACETS fit comparison requires a JMLE-style residual basis
#' (see [facets_fit_review()]).
#' ConQuest MML parameter fit also has a different target: residuals are
#' aggregated through the response design for each parameter, then averaged
#' over posterior ability distributions. Its Infit averages a ratio, rather
#' than dividing separately averaged residual sums and variances. Posterior
#' integration of row residuals alone does not reproduce that statistic. Native
#' marginal category-count diagnostics are another target. See the
#' [ConQuest technical manual](https://conquestmanual.acer.org/s3-00.html).
#'
#' **Heuristic misfit-screening guidelines (Bond & Fox, 2015):**
#' - MnSq < 0.5: overfit (too predictable; may inflate reliability)
#' - MnSq 0.5--1.5: productive for measurement
#' - MnSq > 1.5: underfit (noise degrades measurement)
#' - \eqn{|\mathrm{ZSTD}| \ge 2}: package convention for the approximate-normal review flag;
#'   not a calibrated 5\% hypothesis test, especially after parameter estimation
#'   and repeated screening across elements
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
#' These residual-PCA summaries are not a DIMTEST/UNIDIM implementation.
#' DIMTEST-style essential-unidimensionality tests work at an item-response
#' layer and require an explicit decision about how many-facet rating data are
#' collapsed, conditioned, or adjusted for rater/task/facet effects. For
#' manuscripts, combine global/element fit, residual PCA, and local-dependence
#' screens, and use limited wording such as "evidence consistent with essential
#' unidimensionality under the specified facet structure" rather than
#' "unidimensionality was established."
#'
#' @section Reading key components:
#' Practical interpretation often starts with:
#' - `overall_fit`: global infit/outfit and degrees of freedom.
#' - `reliability`: facet-level model/real separation and reliability. `MML`
#'   uses model-based `ModelSE` values where available; `JML` keeps these
#'   quantities as exploratory approximations.
#' - `fit`: element-level misfit scan (`Infit`, `Outfit`, `ZSTD`).
#' - `unexpected`, `fair_average`, `displacement`: targeted QC bundles.
#'   For `GPCM`, `fair_average` is retained with an unavailable
#'   status because that compatibility calculation is outside the documented
#'   generalized-model contract.
#' - `approximation_notes`: method notes for SE/CI/reliability summaries.
#'
#' @section Interpreting output:
#' Testlet and shared-rater fits use [mfrm_response_diagnostics()] instead of
#' this ordinary-model route. Their same-data posterior predictive residual
#' summaries do not inherit ordinary fit cutoffs, ZSTD or p-values.
#'
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
#' Reliability tables report model and fit-adjusted indices from observed
#' variance minus mean squared SE, truncated at zero. Fit-adjusted values are
#' not confidence bounds; `JML` entries remain exploratory. Separation, strata, and reliability follow the
#' Wright & Masters (1982) conventions:
#' \eqn{G = \mathrm{TrueSD}/\mathrm{RMSE}},
#' \eqn{R = G^2 / (1 + G^2)}, and \eqn{H = (4G + 1) / 3}.
#'
#' Tables record finite-estimate counts and the available SEs on those same
#' levels. Non-finite estimates are excluded from the spread and SE summaries.
#' If any finite estimate lacks a valid SE, reliability, separation, strata and
#' error-adjusted spread are unavailable, rather than combining different sets
#' of levels. Excluded levels or incomplete uncertainty prevent a facet summary
#' from supporting formal reporting. A facet containing fixed values also has
#' no aggregate reliability, separation or strata: a fixed anchor's zero
#' conditional variance is not evidence of perfect measurement. Fixing a group
#' mean alone does not fix each group member. For EAP Persons, this separation-based
#' index is distinct from posterior-variance EAP reliability. High rater
#' separation means distinguishable rater measures, not high rater agreement.
#'
#' Facet SEs from a regularized information matrix or an observation-table
#' fallback remain diagnostic approximations. They are labelled explicitly and
#' cannot authorize ordinary confidence-interval reporting. Numerical convergence
#' is reviewed separately from support for inference; switching from JML to MML
#' does not by itself establish valid SEs or intervals. Recompute older diagnostic
#' objects with `diagnose_mfrm(fit)` before reporting; the existing fit can be used
#' without refitting the model.
#'
#' @section Typical workflow:
#' 1. Start with `diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")`.
#' 2. Inspect `summary(diag)` and use `diagnostic_basis` to separate legacy residual evidence from strict marginal evidence.
#' 3. If needed, rerun with residual PCA (`"overall"` or `"both"`).
#'
#' @return
#' An object of class `mfrm_diagnostics` including:
#' - `obs`: observed/expected/residual-level table
#' - `measures`: facet/person fit table (`Infit`, `Outfit`, `ZSTD`,
#'   `PTMEA`, `ModelSE`, `RealSE`, `CI_Lower`, `CI_Upper`, `CI_Level`,
#'   `CI_Method`). `Fixed = TRUE` identifies location values fixed by anchors
#'   or identification constraints. Their sampling SEs and normal intervals
#'   are `NA` (not estimated, not a numerical failure). Free levels retain
#'   their estimator-specific SEs. Aggregate reliability/separation and facet
#'   variability tests are unavailable for facets containing fixed values;
#'   response-fit diagnostics remain available. This convention differs from
#'   [mfrm_facet_intervals()], which represents a fixed target's conditional
#'   covariance as zero, also without an interval. Neither route propagates
#'   uncertainty in supplied anchor values.
#' - `overall_fit`: overall fit summary
#' - `fit`: element-level fit diagnostics
#' - `reliability`: facet-level model/real separation and reliability
#' - `precision_profile`: one-row summary of the active precision tier and its
#'   recommended use
#' - `precision_review`: package-native checks for SE, CI, and reliability
#' - `parameter_uncertainty`: MML observed-information uncertainty for
#'   structural parameters when available (`steps`, and `GPCM`
#'   `slopes` on both log and positive scales), plus covariance status metadata.
#'   Step `CIEligible` and `CIUse` retain the source fit's restrictions;
#'   non-unit observation-weight bands are diagnostic only. GPCM MML relative
#'   slopes use the separate [confint.mfrm_fit()] checks and pointwise log-Wald
#'   approximation; inspect `CIEligible`, `CIUse` and `InferenceReview`.
#' - `facet_precision`: facet-level precision summary by distribution basis and
#'   SE mode
#' - `facets_chisq`: fixed/random facet variability summary
#' - `interactions`: top interaction diagnostics
#' - `interrater`: inter-rater agreement bundle (`summary`, `pairs`) including
#'   agreement and rater-severity spread indices
#' - `unexpected`: unexpected-response bundle
#' - `fair_average`: adjusted-score reference bundle (reported as unavailable
#'   for `GPCM`)
#' - `displacement`: displacement diagnostics bundle
#' - `approximation_notes`: method notes for SE/CI/reliability summaries
#' - `diagnostic_basis`: guide to the statistical target of each diagnostic path
#' - `fit_standardization`: guide to the df convention behind fit ZSTD values
#' - `fit_readiness`, `fit_readiness_components`, and
#'   `fit_readiness_parameters`: the source fit's versioned readiness decision;
#'   diagnostic computation does not override global fit readiness. GPCM slope
#'   parameter rows additionally retain the output-specific interval decision
#' - `marginal_fit`: optional strict marginal-fit companion based on
#'   posterior-expected first-order category counts, with classification coverage
#' - `residual_pca_overall`: optional overall PCA object
#' - `residual_pca_by_facet`: optional facet PCA objects
#' - `replay_inputs`: diagnostic settings retained for reproducible export,
#'   including fit standardization, interaction selection, and PCA limits
#'
#' @seealso [fit_mfrm()], [analyze_residual_pca()], [build_visual_summaries()],
#'   [mfrmr_visual_diagnostics], [mfrmr_reporting_and_apa]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("diagnose_mfrm", run.dontrun = TRUE)
#' \dontrun{
#' # Load the package and example ratings
#' library(mfrmr)
#' ratings <- load_mfrmr_data("example_operational")
#'
#' # Fit the model
#' fit <- fit_mfrm(
#'   data = ratings,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   rating_min = 1, rating_max = 4, category_policy = "preserve",
#'   method = "MML",
#'   model = "RSM"
#' )
#'
#' # Check model fit and the support for standard errors and intervals
#' diagnostics <- diagnose_mfrm(fit)
#' # Check each residual's basis and the coverage of aggregate fit statistics
#' head(diagnostics$obs[, c("PersonMeasureBasis", "ConditionalVariance",
#'   "PearsonAvailable", "ResponseReason")])
#' diagnostics$overall_fit[, c("Infit", "Outfit", "FitUnavailableRows", "FitReason")]
#' diagnostic_summary <- summary(diagnostics)
#' diagnostic_summary$decision
#'
#' diagnostic_summary$key_warnings # Issues to investigate, if present
#' diagnostic_summary$top_fit      # Most unusual residual-based fit statistics
#'
#' # Distinguish residual-based checks from marginal model checks
#' diagnostic_summary$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]
#' }
#' @section References:
#' - Wright, B. D., & Masters, G. N. (1982). *Rating scale analysis*.
#'   MESA Press. (G/R/H separation, reliability, and strata
#'   formulas summarized in `s_diag$reliability` follow this
#'   convention.)
#' - Wright, B. D., & Linacre, J. M. (1994). Reasonable mean-square
#'   fit values. *Rasch Measurement Transactions, 8*(3), 370.
#'   (Source for the 0.5-1.5 Infit / Outfit heuristic review interval that
#'   `s_diag$key_warnings` and `misfit_thresholds` apply.)
#' - Linacre, J. M. (1989). *Many-Facet Rasch Measurement*. MESA
#'   Press. (FACETS Tables 6 + 7 correspond to the per-facet
#'   element measures, fit, and chi-square heterogeneity screen
#'   exposed via `s_diag$reliability` and `s_diag$facets_chisq`.)
#' - Bond, T. G., & Fox, C. M. (2015). *Applying the Rasch model:
#'   Fundamental measurement in the human sciences* (3rd ed.).
#'   Routledge. (Reference text for the Rasch-family fit
#'   conventions exposed by this helper.)
#' - Linacre, J. M. (2002). What do Infit and Outfit, Mean-square
#'   and Standardized mean? *Rasch Measurement Transactions,
#'   16*(2), 878.
#' - Linacre, J. M. (2026). *A user's guide to Facets Rasch-model computer
#'   programs*. Winsteps.com. (WHEXACT / FACETS standardized fit df notes.)
#' @export
#' @order 1
diagnose_mfrm <- function(fit,
                          interaction_pairs = NULL,
                          top_n_interactions = 20,
                          whexact = FALSE,
                          fit_df_method = "facets",
                          diagnostic_mode = c("both", "legacy", "marginal_fit"),
                          residual_pca = c("none", "overall", "facet", "both"),
                          pca_max_factors = 10L) {
  fit <- .mfrm_validate_diagnostic_fit(fit)
  fit_df_method <- match_fit_df_method(fit_df_method)
  diagnostic_mode <- match.arg(diagnostic_mode)
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))

  out <- mfrm_diagnostics(
    fit,
    interaction_pairs = interaction_pairs,
    top_n_interactions = top_n_interactions,
    whexact = whexact,
    fit_df_method = fit_df_method,
    diagnostic_mode = diagnostic_mode,
    residual_pca = residual_pca,
    pca_max_factors = pca_max_factors
  )
  out$parameter_uncertainty$source <- mfrm_native_location_identity(fit)
  class(out) <- c("mfrm_diagnostics", class(out))
  out
}

mfrm_ic_contract_required_fields <- function() {
  c(
    "ICContractVersion", "N", "ResponseRows", "WeightedResponseTotal",
    "Persons", "Npar", "LogLik", "Deviance", "WeightPolicy",
    "ICEligible", "ICSelectable", "ICStatus", "ICSampleSize",
    "ICSampleSizeBasis",
    "AIC", "BIC", "SABIC", "SABICSelectable", "AICFormula",
    "BICFormula", "SABICFormula", "IntegrationEvaluationId",
    "ICQuadraturePoints", "ICIntegrationTier", "ICIntegrationStatus",
    "ICIntegrationSelectable"
  )
}

mfrm_ic_stored_value_matches <- function(actual,
                                         expected,
                                         tolerance = 1e-10) {
  if (length(actual) == 0L || length(expected) == 0L) return(FALSE)
  actual <- actual[1]
  expected <- expected[1]
  if (is.numeric(expected) || is.integer(expected)) {
    actual <- suppressWarnings(as.numeric(actual))
    expected <- suppressWarnings(as.numeric(expected))
    if (is.na(actual) && is.na(expected)) return(TRUE)
    if (!is.finite(actual) || !is.finite(expected)) {
      return(identical(actual, expected))
    }
    return(abs(actual - expected) <=
             tolerance * max(1, abs(expected)))
  }
  if (is.logical(expected)) {
    actual <- suppressWarnings(as.logical(actual))
    if (is.na(actual) && is.na(expected)) return(TRUE)
    return(identical(actual, expected))
  }
  actual <- as.character(actual)
  expected <- as.character(expected)
  if (is.na(actual) && is.na(expected)) return(TRUE)
  identical(actual, expected)
}

mfrm_extract_fit_ic_contract <- function(fit, tolerance = 1e-10) {
  summary_row <- as.data.frame(fit$summary %||% data.frame(),
                               stringsAsFactors = FALSE)
  if (nrow(summary_row) > 0L) summary_row <- summary_row[1, , drop = FALSE]
  config <- fit$config %||% list()
  prep <- fit$prep %||% list()
  method <- public_mfrm_method_label(toupper(as.character(
    config$method %||%
      if ("Method" %in% names(summary_row)) summary_row$Method[1] else NA_character_
  )[1]))
  model <- toupper(as.character(
    config$model %||%
      if ("Model" %in% names(summary_row)) summary_row$Model[1] else NA_character_
  )[1])
  current_contract <- nrow(summary_row) > 0L &&
    "ICContractVersion" %in% names(summary_row) &&
    identical(
      as.character(summary_row$ICContractVersion[1]),
      mfrm_ic_contract_version()
    )
  can_rebuild <- !is.null(fit$opt$value) && length(fit$opt$value) > 0L &&
    !is.null(fit$opt$par) && is.list(prep) && is.data.frame(prep$data)

  expected <- if (can_rebuild) {
    tryCatch(
      build_mfrm_ic_contract(
        loglik = -as.numeric(fit$opt$value[1]),
        npar = length(fit$opt$par),
        prep = prep,
        config = config,
        method = method
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  if (current_contract && !is.null(expected)) {
    required <- mfrm_ic_contract_required_fields()
    missing_fields <- setdiff(required, names(summary_row))
    comparable_fields <- intersect(required, names(expected))
    field_checks <- vapply(comparable_fields, function(field) {
      field %in% names(summary_row) && mfrm_ic_stored_value_matches(
        summary_row[[field]][1], expected[[field]], tolerance = tolerance
      )
    }, logical(1))
    stored_consistent <- length(missing_fields) == 0L && all(field_checks)
    mismatch <- c(missing_fields, names(field_checks)[!field_checks])
    contract_state <- if (stored_consistent) "current" else "current_inconsistent"
    values <- expected
  } else {
    stored_consistent <- NA
    mismatch <- if (current_contract) "contract_rebuild_unavailable" else
      "missing_current_contract"
    contract_state <- "legacy_or_unknown"
    response_rows <- if (is.data.frame(prep$data)) nrow(prep$data) else NA_integer_
    weighted_total <- if (is.data.frame(prep$data) && "Weight" %in% names(prep$data)) {
      legacy_weights <- suppressWarnings(as.numeric(prep$data$Weight))
      if (length(legacy_weights) > 0L && all(is.finite(legacy_weights))) {
        sum(legacy_weights)
      } else {
        NA_real_
      }
    } else {
      as.numeric(response_rows)
    }
    persons <- if (is.data.frame(prep$data) && "Person" %in% names(prep$data)) {
      person_values <- as.character(prep$data$Person)
      length(unique(person_values[!is.na(person_values)]))
    } else {
      suppressWarnings(as.integer(config$n_person %||% NA_integer_))
    }
    stored_loglik <- if ("LogLik" %in% names(summary_row)) {
      suppressWarnings(as.numeric(summary_row$LogLik[1]))
    } else if (!is.null(fit$opt$value)) {
      -suppressWarnings(as.numeric(fit$opt$value[1]))
    } else {
      NA_real_
    }
    stored_aic <- if ("LegacyAIC" %in% names(summary_row)) {
      suppressWarnings(as.numeric(summary_row$LegacyAIC[1]))
    } else if ("AIC" %in% names(summary_row)) {
      suppressWarnings(as.numeric(summary_row$AIC[1]))
    } else {
      NA_real_
    }
    stored_bic <- if ("LegacyBIC" %in% names(summary_row)) {
      suppressWarnings(as.numeric(summary_row$LegacyBIC[1]))
    } else if ("BIC" %in% names(summary_row)) {
      suppressWarnings(as.numeric(summary_row$BIC[1]))
    } else {
      NA_real_
    }
    legacy_n <- if ("N" %in% names(summary_row)) {
      suppressWarnings(as.numeric(summary_row$N[1]))
    } else {
      NA_real_
    }
    values <- list(
      ICContractVersion = "legacy_or_unknown",
      N = legacy_n,
      ResponseRows = as.integer(response_rows),
      WeightedResponseTotal = as.numeric(weighted_total),
      Persons = as.integer(persons),
      Npar = if (!is.null(fit$opt$par)) length(fit$opt$par) else NA_integer_,
      LogLik = stored_loglik,
      Deviance = if (is.finite(stored_loglik)) -2 * stored_loglik else NA_real_,
      WeightPolicy = if (is.data.frame(prep$data)) {
        mfrm_ic_weight_policy(prep, config)
      } else {
        "unknown"
      },
      ICEligible = FALSE,
      ICSelectable = FALSE,
      ICStatus = "suppressed_legacy_contract",
      ICSampleSize = NA_real_,
      ICSampleSizeBasis = "legacy_or_unknown",
      AIC = NA_real_,
      BIC = NA_real_,
      SABIC = NA_real_,
      SABICSelectable = FALSE,
      AICFormula = NA_character_,
      BICFormula = NA_character_,
      SABICFormula = NA_character_,
      IntegrationEvaluationId = mfrm_ic_integration_evaluation_id(method, config),
      ICQuadraturePoints = NA_integer_,
      ICIntegrationTier = "legacy_or_unknown",
      ICIntegrationStatus = "legacy_or_unknown",
      ICIntegrationSelectable = FALSE,
      LegacyAIC = stored_aic,
      LegacyBIC = stored_bic,
      LegacyICSampleSize = legacy_n,
      LegacyICSampleSizeBasis = if (!is.null(config$weight_col)) {
        "sum_weights"
      } else {
        "row_count"
      }
    )
  }

  tibble::tibble(
    Model = model,
    Method = method,
    nobs = values$ResponseRows,
    ResponseRows = values$ResponseRows,
    WeightedN = values$WeightedResponseTotal,
    WeightedResponseTotal = values$WeightedResponseTotal,
    Persons = values$Persons,
    ICSampleSize = values$ICSampleSize,
    ICSampleSizeBasis = values$ICSampleSizeBasis,
    Npar = values$Npar,
    npar = values$Npar,
    LogLik = values$LogLik,
    Deviance = values$Deviance,
    AIC = values$AIC,
    BIC = values$BIC,
    SABIC = values$SABIC,
    SABICSelectable = values$SABICSelectable,
    WeightPolicy = values$WeightPolicy,
    ICEligible = values$ICEligible,
    ICSelectable = values$ICSelectable,
    ICStatus = values$ICStatus,
    ICContractVersion = values$ICContractVersion,
    AICFormula = values$AICFormula,
    BICFormula = values$BICFormula,
    SABICFormula = values$SABICFormula,
    IntegrationEvaluationId = values$IntegrationEvaluationId,
    ICQuadraturePoints = values$ICQuadraturePoints,
    ICIntegrationTier = values$ICIntegrationTier,
    ICIntegrationStatus = values$ICIntegrationStatus,
    ICIntegrationSelectable = values$ICIntegrationSelectable,
    LegacyAIC = values$LegacyAIC,
    LegacyBIC = values$LegacyBIC,
    LegacyICSampleSize = values$LegacyICSampleSize,
    LegacyICSampleSizeBasis = values$LegacyICSampleSizeBasis,
    ICContractState = contract_state,
    StoredICConsistent = stored_consistent,
    StoredICMismatch = paste(unique(mismatch), collapse = ";")
  )
}

mfrm_ic_fit_check <- function(fit, contract, information = NULL, allow_singleton = FALSE) {
  readiness <- mfrmr_get_readiness_record(fit)$fit
  ordinary_ready <- mfrm_inference_ready(fit)
  out <- list(
    eligible = ordinary_ready,
    basis = "fit_inference_readiness",
    review = if (ordinary_ready) "" else as.character(mfrm_fit_decision_summary(readiness)$Why[1])
  )
  rasch_singleton <- isTRUE(allow_singleton) && fit$config$model %in% c("RSM", "PCM") &&
    identical(as.character(readiness$CategoryState[1]), "weak_information")
  if ((!identical(fit$config$model, "GPCM") &&
       !isTRUE(fit$config$population_spec$active) && !rasch_singleton) ||
      !identical(fit$config$method, "MML")) return(out)

  out$eligible <- FALSE
  out$basis <- "mml_local_solution_information"
  refuse <- function(reason) {
    out$review <- reason
    out
  }
  if (!isTRUE(contract$StoredICConsistent[1]) ||
      !isTRUE(contract$ICSelectable[1])) {
    return(refuse("Current, consistent MML likelihood metadata and a sufficiently resolved quadrature grid are required."))
  }
  identification <- mfrm_mml_inference_identification(fit)
  if (isTRUE(identification$excluded)) return(refuse(identification$detail))
  category_ready <- identical(as.character(readiness$CategoryState[1]), "adequate")
  singleton_review <- FALSE
  if (!category_ready && isTRUE(allow_singleton) &&
      fit$config$model %in% c("RSM", "PCM", "GPCM") &&
      identical(as.character(readiness$CategoryState[1]), "weak_information")) {
    # Count-based warnings do not replace a likelihood/information check.
    # GPCM retains its bootstrap-only exception; Rasch locations admit Wald
    # approximations after the separate full-information and quadrature review.
    audit <- tryCatch(audit_mfrm_category_support(fit$prep, fit$config,
      build_param_sizes(fit$config)), error = function(e) NULL)
    counts <- audit$category_table$WithinScopeCount
    local <- audit$local_support
    local_singleton <- rasch_singleton && is.data.frame(local) && nrow(local) > 0L &&
      all(local$MissingCategories == "") && any(local$MinimumObservedCategoryCount == 1L)
    singleton_review <- isTRUE(audit$readiness$Complete) &&
      identical(audit$readiness$CategoryState, "weak_information") &&
      identical(audit$readiness$ReasonCodes, "weak_category_information") &&
      isTRUE(audit$readiness$UnsupportedStepCoordinates == 0L) &&
      isTRUE(audit$readiness$UnsupportedCategoryContrasts == 0L) &&
      length(counts) > 0L && all(is.finite(counts)) &&
      all(counts > 0L) && (any(counts == 1L) || local_singleton)
    if (rasch_singleton) singleton_review <- singleton_review &&
      is.data.frame(local) && nrow(local) > 0L && all(local$MissingCategories == "")
  }
  if (!identical(as.character(readiness$InputState[1]), "pass") ||
      !(category_ready || singleton_review)) {
    return(refuse("Input or category support requires review before MML inference."))
  }
  if (!identical(as.character(readiness$NumericalState[1]), "ready") ||
      !identical(fit$opt$optimizer_diagnostics$ConvergenceSeverity, "pass")) {
    return(refuse("The MML optimizer has not passed the numerical convergence checks."))
  }
  gradient_tolerance <- as.numeric(fit$opt$optimizer_diagnostics$GradientReviewTolerance)
  if (length(gradient_tolerance) != 1L || !is.finite(gradient_tolerance) ||
      gradient_tolerance <= 0 || gradient_tolerance > 1e-4) gradient_tolerance <- 1e-4
  # Reevaluate the retained solution; a saved positive-curvature label must
  # not authorize comparison after the data or parameter vector has changed.
  information <- information %||% tryCatch({
    covariance <- compute_mml_parameter_covariance(fit,
      solution_tolerance = gradient_tolerance)
    covariance$solution_information %||% list(status = "unavailable", detail = covariance$detail)
  }, error = function(e) list(status = "unavailable", detail = conditionMessage(e)))
  if (identical(information$status, "rejected_local_solution")) {
    return(refuse(information$detail))
  }
  if (!identical(information$status, "evaluated_diagnostic_only")) {
    return(refuse(paste("MML solution information is unavailable:", information$detail)))
  }
  evaluation <- information$evaluation_summary
  eigenvalues <- information$eigenvalue_summary
  objective <- fit$opt$value
  reevaluated <- evaluation$ReevaluatedObjective
  gradient <- evaluation$GradientMaxAbs
  smallest <- eigenvalues$Smallest
  scale <- eigenvalues$AbsoluteScale
  scalars <- list(objective, reevaluated, gradient, smallest, scale, information$free_dimension)
  if (!all(vapply(scalars, function(x) is.numeric(x) && !is.complex(x) &&
      length(x) == 1L && is.finite(x), logical(1))) ||
      length(fit$opt$par) != information$free_dimension || gradient < 0 || smallest > scale) {
    return(refuse("Finite, dimensionally consistent likelihood, gradient and information are required."))
  }
  if (abs(reevaluated - objective) > 1e-10 * max(1, abs(objective))) {
    return(refuse("The retained objective does not match reevaluation at the supplied data and parameters."))
  }
  if (gradient > gradient_tolerance) {
    return(refuse("The reevaluated MML gradient exceeds the numerical convergence tolerance."))
  }
  # Weak raw-coordinate conditioning alone need not reject a freshly verified,
  # unregularized inverse. Indefinite/zero curvature still cannot be admitted.
  inverse_review <- information$inverse_review
  refined <- isTRUE(inverse_review$Verified) &&
    all(c("RelativeChange", "InverseResidual", "CurvatureScaledGradient") %in% names(inverse_review)) &&
    all(vapply(inverse_review[c("RelativeChange", "InverseResidual", "CurvatureScaledGradient")],
               length, integer(1)) == 1L) &&
    all(is.finite(c(inverse_review$RelativeChange, inverse_review$InverseResidual,
                    inverse_review$CurvatureScaledGradient))) &&
    min(inverse_review$RelativeChange, inverse_review$InverseResidual,
        inverse_review$CurvatureScaledGradient) >= 0 &&
    inverse_review$RelativeChange <= 1e-3 && inverse_review$InverseResidual <= 1e-6 &&
    inverse_review$CurvatureScaledGradient <= 1e-4
  if (!is.null(inverse_review) && !isTRUE(refined)) {
    return(refuse(inverse_review$Detail %||% "Information refinement did not pass."))
  }
  if (scale <= 0 || smallest <= 0 ||
      (smallest <= scale * sqrt(.Machine$double.eps) && !refined)) {
    return(refuse("MML information is singular, indefinite or too ill-conditioned for ordinary inference."))
  }
  out$eligible <- TRUE
  if (refined) out$caution <- inverse_review$Detail
  if (singleton_review) out$caution <- trimws(paste(out$caution %||% "",
    if (rasch_singleton)
      "Some score categories occur only once within a fitted or local scope. Model-information intervals remain asymptotic approximations; small samples can affect bias and coverage. The source fit retains its category-support review."
    else paste("Singleton score categories: the point estimate is retained for the basic bootstrap",
      "after numerical and information checks. Sparse category support may affect accuracy;",
      "this does not establish Wald-interval eligibility.")))
  out$review <- paste(
    "The retained MML solution has matching likelihood, a small gradient and",
    "positive unregularized local information. Check starting-value and",
    "integration sensitivity when comparing close criteria; this is not a",
    "global-maximum certificate or permission for slope intervals or an LRT.", out$caution %||% ""
  )
  out
}

#' Compare two or more fitted MFRM models
#'
#' Produce a side-by-side comparison of multiple [fit_mfrm()] results using
#' AIC, Person-based BIC, Sclove SABIC, log-likelihood, and free-parameter
#' counts. An ordinary RSM paired with a testlet or shared-rater fit instead
#' follows a descriptive facet-effect route, with no automatic ranking.
#' For ordinary models, when exactly
#' two models are supplied and the current conservative nesting review passes,
#' a likelihood-ratio test is included.
#'
#' GPCM MML and RSM/PCM MML with an estimated normal population use a separate
#' local-solution check for information criteria, without refitting. It checks
#' the retained likelihood, terminal gradient and positive unregularized
#' observed information. `ICFitEligible`, `ICFitBasis` and `ICFitReview` record
#' that decision; `ICComparable` also requires the common data/likelihood and
#' integration checks. This can allow IC comparison while `InferenceReady`
#' remains false for intervals or tests. It shares the joint-information
#' calculation and workspace budget with [confint.mfrm_fit()]. Unknown, unstable or unavailable results
#' remain ineligible with a reason. A local check does not prove a global
#' maximum or negligible integration error; examine starting-value and
#' [mml_quadrature_sensitivity()] results when the decision is close.
#' Known marginal response-dimension deficits and complete-pattern first-order
#' rank failures also prevent inferential model comparison, even with a positive
#' numerical Hessian. More Persons with the same response designs improve
#' precision but cannot repair a structural dimension deficit.
#' Positive but ill-conditioned information may pass with a caution after
#' numerical refinement, unregularized inversion and curvature-scaled gradient
#' checks. A warning is emitted and retained in `ICFitCaution`, `ICFitReview`
#' and the requested LRT's interpretation; it does not establish
#' finite-sample accuracy for model ranking or tests.
#' @param ... For `compare_mfrm()`: Two or more `mfrm_fit` objects, or exactly one ordinary
#'   `mfrm_fit` and one [fit_mfrm_testlet()] or [fit_mfrm_random_rater()] result.
#'   The latter route is descriptive; see the extended-model section.
#'
#' For `plot()`: Unused.
#' @param object An extended-model comparison returned by `compare_mfrm()`.
#' @param labels Optional character vector of labels for each model.
#'   If `NULL`, labels are generated from model/method combinations.
#' @param warn_constraints Logical. If `TRUE` (the default), emit a warning
#'   when models use different score coding or constraint settings. Ranking is
#'   suppressed whenever this basis differs, regardless of whether the warning
#'   is displayed.
#' @param nested Logical. Set to `TRUE` only when the supplied models are
#'   known to be nested and fitted with the same likelihood basis on the same
#'   observations. The default is `FALSE`, in which case no likelihood-ratio
#'   test is reported. When `TRUE`, the function still runs a conservative
#'   structural nesting review and computes the LRT only for supported nesting patterns.
#' @param response_diagnostics Optional list of two saved
#'   [mfrm_response_diagnostics()] outputs in the same order as the fits.
#'   Supported for an ordinary RSM paired with a testlet or shared-rater fit.
#'   This adds descriptive predictive comparisons without computing integrals.
#' @param person_scores Optional list of two saved [score_mfrm_persons()]
#'   results in fit order for an ordinary versus extended RSM comparison.
#'   Both must score the same complete source roster, return the same Person
#'   IDs and use the same conditional interval level. No scoring is run here.
#' @details Models should be fit to the **same data** (same rows, same person/facet
#' columns) for the comparison to be meaningful. The function checks that
#' observation counts match and warns otherwise.
#'
#' Information-criterion ranking is reported only when all candidates are
#' eligible solutions under the package's current `MML` contract, use the same
#' prepared observations, score coding, constraints, formula contract, and
#' integration-evaluation identity, are eligible under the weighting policy,
#' and have a selectable integration tier. For fixed-facet MML,
#' `ICSampleSize` is the number of independent
#' Persons and `ICSampleSizeBasis` is `"person_count"`; response rows and their
#' weighted total are retained separately. Explicit all-unit weights are
#' eligible, whereas every non-unit observation-weight fit fails closed under
#' the current package contract.
#'
#' Canonical `AIC`, `BIC`, and `SABIC` are recomputed from the retained
#' objective, optimizer-vector dimension, and Person count. A stale stored
#' value, a legacy object without the current contract identity, a JML fit, or
#' an ineligible weighted fit cannot enter ranking. Earlier raw values may be
#' shown only as explicitly labelled `LegacyAIC` and `LegacyBIC` fields.
#'
#' **Integration selection guard**: raw canonical criteria are retained at
#' every valid MML quadrature count, but automatic comparison is screening-only
#' below 15 points and review-only at 15--30 points. `ICSelectable` and
#' `ICIntegrationSelectable` become `TRUE` at q>=31. q=31 is a starting grid,
#' not a universal guarantee: close or consequential comparisons should be
#' reevaluated at a denser shared grid (normally q>=61). The guard also applies
#' to `nested = TRUE`, so a coarse-grid LRT cannot bypass it.
#'
#' **Nesting**: Two models are *nested* when one is a special case of the
#' other obtained by imposing equality constraints.  The most common
#' nesting in MFRM is RSM (shared thresholds) inside PCM
#' (item-specific thresholds).  Models that differ only in estimation
#' method (MML vs JML) on the same specification are not nested in the
#' usual sense, and their information criteria do not share the common MML
#' contract. Do not use either the LRT or IC ranking as a direct cross-method
#' comparison.
#'
#' In the **current `mfrmr` model space**, the automatic nesting review is
#' intentionally conservative. It currently supports the following
#' restrictions under shared data and shared constraints:
#' - `RSM` nested inside `PCM` when the `PCM` fit has an explicit
#'   `step_facet` and retains every `facet_interactions` term in the RSM.
#'   The PCM may add interaction terms; omitting an RSM interaction is not
#'   repaired by freeing category steps;
#' - `PCM` nested inside `GPCM` with the same step facet, population
#'   design, other facet/step constraints and interactions. The only additional
#'   parameters must be G-1 relative log-slope contrasts for G slope levels.
#'   The GPCM slope owner may differ from the shared step owner;
#' - same-family additive-vs-interaction comparisons when the smaller fit's
#'   `facet_interactions` set is a subset of the larger fit's set.
#'
#' Cross-method comparisons, comparisons that change anchors/dummying/centering,
#' and same-family comparisons that do not add fixed interaction terms are not
#' automatically promoted to LRT claims.
#'
#' Automatic nesting also requires a shared population specification: either
#' both fits use the fixed standard-normal population or both use the same
#' estimated-normal design, with the same declared columns and design values
#' aligned by Person. Coefficient and variance estimates may differ. Row and
#' column permutations are aligned; other recodings, changed designs, and
#' unavailable design metadata require separate review. Passing this check
#' does not grant inference readiness to an estimated-population fit. For a
#' matched PCM/GPCM pair, the local-solution checks described above supply
#' the numerical requirement independently of slope-interval availability.
#'
#' The **likelihood-ratio test (LRT)** is reported only when exactly two
#' models are supplied, `nested = TRUE`, the structural nesting review passes, and the
#' difference in the number of parameters is positive:
#'
#' \deqn{\Lambda = -2 (\ell_{\mathrm{restricted}} - \ell_{\mathrm{full}})
#'   \sim \chi^2_{\Delta p}}
#'
#' The LRT is asymptotically valid only under its regularity assumptions.
#' With small samples or boundary/singular conditions, the reference
#' chi-square p-value can be incorrect. At large Person counts, a very small
#' practical improvement can also become statistically significant. Read the
#' p-value as formal nested-fit evidence, not as an automatic practical model
#' preference or evidence that a subscore is useful.
#' @section Ordinary versus extended RSMs:
#' Exactly one ordinary RSM MML fit and one testlet/shared-rater fit return an
#' `mfrm_extended_comparison`, with `models`, `checks`, `effects` and `notes`.
#' The information-criterion and LRT sections describe ordinary `mfrm_fit`
#' comparisons only. This route does not rank models, compute difference
#' intervals, or provide an AIC/BIC preference. `nested = TRUE` is refused.
#'
#' Both fits must use identical observed rating events, including repeated
#' events, category coding, unit weights and fixed facets except for the
#' deliberately changed rater treatment. Use unanchored additive severity
#' facets, `noncenter_facet = "Person"`, no facet shrinkage, and matching
#' person/score column roles. Category recoding must preserve the same
#' consecutive adjacent-category scale. These checks cannot be disabled by
#' `warn_constraints = FALSE`.
#'
#' Match the population assumption: ordinary `population_formula = ~1` with
#' both extension defaults estimates a common normal ability SD; the ordinary
#' default and extension `person_sd = 1` impose known N(0,1). Covariate-dependent
#' populations and fixed nonunit SDs have no matching ordinary route here.
#' Population intercepts and free step locations may use different origins;
#' they are recorded, not compared as substantive mean differences.
#'
#' The source roster and omitted counts must agree. New ordinary fits retain
#' `prep$omitted_data` and `prep$omitted_input_rows`. When scores were omitted,
#' their identities must match, including repeated events. Older fits without
#' this provenance require refitting from the complete assigned-score roster.
#' Missing IDs, excluded weights and population-data omissions are unsupported.
#'
#' Effects are centered at the unweighted mean of the same complete set of
#' levels within each facet. This preserves pairwise contrasts while removing
#' arbitrary facet origins. `SourceReference`/`SourceComparison` retain raw
#' values and `CenterReference`/`CenterComparison` retain the subtracted means.
#' `Difference` is comparison minus reference. Conditional rater modes and
#' fixed coefficients are labeled separately; shrinkage is not evidence of
#' greater accuracy or rater quality. Failed numerical checks or a missing
#' level estimate withhold affected differences without dropping source rows.
#' Native inference readiness remains explicit even when descriptive numerical
#' checks pass. No SE or confidence interval for model differences is implied.
#'
#' With `response_diagnostics`, `responses` retains matched category
#' probabilities, means, full mixture variances and grouped descriptive
#' Infit/Outfit. Both saved diagnostics must use the same probability target,
#' selected event multiset and shared `group_by` identifier columns. Matching
#' uses event contents, including repeated-event multiplicities and missing
#' scores, rather than assuming row numbers agree. Identical repeated events
#' match in source occurrence order. `responses$events` and `responses$rows`
#' retain the matching and original row numbers. Missing or unavailable rows
#' remain explicit and differences involving them are withheld.
#'
#' All observed source events condition each model; selecting a Person's rows
#' does not remove other Persons from shared-rater inference. These are
#' same-data descriptive comparisons with calibration held fixed, not
#' held-out predictive performance. A smaller Infit/Outfit is not evidence
#' of improvement. Do not substitute ordinary plug-in diagnostics or their
#' reference cutoffs. See [plot.mfrm_extended_comparison()] for predictive
#' paired and difference displays.
#'
#' With `person_scores`, `persons$table` compares conditional EAPs after
#' subtracting each fitted population mean, retaining original estimates,
#' origins, posterior SDs and aligned conditional endpoints. Unit Rasch slopes
#' retain the logit unit; centering removes the arbitrary origin, not shrinkage
#' or changes in the assumed population variance. Prior-only and unavailable
#' differences are withheld. Separate conditional intervals are not intervals
#' for model differences or tests of Person differences. Use `metric = "person"`
#' in [plot.mfrm_extended_comparison()] for a paired/difference display.
#'
#' The default facet-only output does not compare Person scores, steps,
#' predictive Infit/Outfit or replacement-rater predictions. For testlets, zero local
#' variance retains the fixed facets; zero shared-rater variance removes rater
#' differences and is not the fixed-rater model. Numerical zero-reduction
#' agreement does not establish model adequacy or variance-boundary inference.
#'
#' Use `plot(comparison, style = "paired")` for matched effects, or
#' `style = "difference"` for differences versus means. See
#' [plot.mfrm_extended_comparison()] for display and accessibility controls.
#' Attach saved results with `mfrm_results(extended_fit, comparison = comparison)`
#' for static reports, tables, figures and replay. No comparison/report step
#' fits, scores or resamples either model.
#' @section Information-criterion diagnostics:
#' In addition to the canonical criteria, the function computes:
#' - **Delta_AIC / Delta_BIC / Delta_SABIC**: difference from the minimum
#'   value in the supplied candidate set.
#'   For AIC, the often-used Delta ranges below 2, 4--7, and above 10 are
#'   rough descriptions of decreasing relative support (Burnham & Anderson,
#'   2002), not significance thresholds, adequacy checks, or universal rules
#'   for BIC and SABIC.
#' - **AkaikeWeight / BICWeight / SABICWeight**: relative candidate-set
#'   weights derived from `exp(-0.5 * Delta)` and normalized only across the
#'   supplied models. They are not posterior probabilities, probabilities of
#'   model truth, or guarantees that the candidate set is adequate.
#' - **Evidence ratios**: pairwise ratios of Akaike weights, quantifying
#'   relative support within this candidate set. A ratio of 5 means only that
#'   one normalized Akaike weight is five times the other; it does not mean
#'   that one model is five times more likely to be true.
#'   Ratios are computed from AIC differences, so they remain available when
#'   both normalized weights round to zero. `LogEvidenceRatio` retains the
#'   natural log ratio; `RatioReason` explains `NA` ratios outside numeric range
#'   (`"overflow"`, `"underflow"`, or `"nonfinite_difference"`, otherwise
#'   `"available"`).
#'
#' AIC, BIC, and SABIC answer different approximation/penalty questions. SABIC
#' is a sensitivity criterion, not a universal tie-breaker. At 22 or fewer
#' Persons its Sclove penalty is non-positive, so `SABICSelectable` and
#' `SABICComparable` are `FALSE` and no automatic SABIC preference is returned.
#' @section What this comparison means:
#' `compare_mfrm()` is a same-basis model-comparison helper. Its strongest
#' claims apply only when the models were fit to the same response data,
#' under a compatible likelihood basis, and with compatible constraint
#' structure.
#' @section What this comparison does not justify:
#' - Do not treat AIC/BIC/SABIC differences as primary evidence when
#'   `table$ICComparable` is `FALSE`.
#' - Do not infer numerical adequacy merely because all models share the same
#'   coarse quadrature identity. Below q=31, raw criteria are diagnostic only
#'   and automatic deltas, weights, preferences, and LRT are suppressed.
#' - Do not interpret the LRT unless `nested = TRUE` and the structural nesting review
#'   in `comparison_basis$nesting_review` passes.
#' - Same-family additive-vs-interaction fits are considered nested only when
#'   all other structural settings match and the smaller model's
#'   `facet_interactions` set is a subset of the larger model's set.
#' - Do not assume that `nested = TRUE` overrides the package's conservative
#'   nesting boundary; unsupported relations remain unsupported.
#' - PCM is the unit-slope reduction of an aligned GPCM. `nested = TRUE` tests
#'   equal relative slopes only after population, constraint, free-dimension
#'   and numerical checks pass. The relation is recorded as `PCM_in_GPCM`.
#'   Default PCM uses a fixed standard-normal population, whereas default GPCM
#'   estimates its mean and variance: changing only `model` does not isolate
#'   slope differences. Supply `population_formula = ~1` and the same
#'   `person_data` to both fits for an estimated-normal comparison.
#'   Unit slopes are interior positive values; the reference is ordinary
#'   asymptotic chi-square with G-1 degrees of freedom, not a boundary mixture.
#'   IC eligibility alone does not establish nesting or qualify slope intervals.
#' - Do not compare models fit to different datasets, different score codings,
#'   or materially different constraint systems as if they were commensurate.
#' - At large Person counts, a small systematic likelihood gain can dominate an
#'   IC difference; boundary or singular fits can also invalidate routine
#'   asymptotics. Evaluate effect size, stability, interpretability, and score
#'   consequences separately.
#' @section Interpreting output:
#' - Lower AIC/BIC values indicate relative support only when
#'   `table$ICComparable` is `TRUE`; use SABIC ranking only when
#'   `table$SABICComparable` is also `TRUE`.
#' - A significant LRT p-value is formal evidence against the restricted
#'   model under the stated assumptions; practical gain, score utility, and
#'   dimensional interpretation require separate checks.
#' - `preferred` identifies the minimum criterion within the supplied
#'   candidate set; it is not an unconditional model verdict.
#' - `evidence_ratios` gives pairwise Akaike-weight ratios (returned only
#'   when Akaike weights can be computed for at least two models).
#' - When comparing more than two models, interpret evidence ratios
#'   cautiously---they do not adjust for multiple comparisons.
#' @section How to read the main outputs:
#' - `table`: first-pass comparison table; start with `ICComparable`,
#'   `ICSelectable`, `ICIntegrationTier`, `Model`, `Method`, `ICStatus`, `AIC`,
#'   `BIC`, and `SABIC`.
#' - `comparison_basis`: records whether IC and LRT claims are defensible for
#'   the supplied models. Inspect `comparison_basis$nesting_review$relation`
#'   and `reason` before reading any LRT output.
#' - `lrt`: nested-model test summary, present only when the requested and
#'   reviewed conditions are met.
#' - `preferred`: candidate preferred by each criterion when those summaries
#'   are available.
#' @section Recommended next step:
#' Inspect `comparison_basis` before writing conclusions. If comparability is
#' weak, treat the result as descriptive and revise the model setup (for
#' example, explicit `step_facet`, common data, or common constraints) before
#' using IC or LRT results in reporting.
#' @section Typical workflow:
#' 1. Fit two models with [fit_mfrm()] (e.g., PCM and GPCM) on the
#'    same prepared rows, explicit step owner, constraints, and MML quadrature
#'    setting. Use at least 31 common quadrature points for selectable ICs.
#' 2. Compare with `compare_mfrm(fit_pcm, fit_gpcm)` and start by checking
#'    `comparison$table$ICComparable`.
#' 3. Inspect `summary(comparison)` for AIC/BIC/SABIC diagnostics and the reasons
#'    for withholding ranking or tests. GPCM IC ranking has separate solution
#'    checks; request the matched PCM/GPCM test with `nested = TRUE`.
#' 4. Use `build_weighting_review(fit_pcm, fit_gpcm)` to inspect which selected
#'    facet levels and information shares were reweighted by the slopes.
#' @return With one extended model, an `mfrm_extended_comparison` containing `models`,
#' `checks`, centered `effects`, `notes`, and source/omission metadata.
#' With ordinary models only, an object of class `mfrm_comparison` (named list) with:
#' - `table`: data.frame of model-level statistics including `LogLik`,
#'   `Deviance`, `Npar` (and equal compatibility alias `npar`), `AIC`, `BIC`,
#'   `SABIC`, criterion deltas and candidate-set weights, `ResponseRows`,
#'   `WeightedResponseTotal`, `Persons`, `ICSampleSize`, formula and integration
#'   identities, integration tier/selectability, stored-value consistency
#'   fields, convergence fields,
#'   `ICComparable`, `SABICComparable`, and `ConstraintComparable`.
#' - `lrt`: data.frame with likelihood-ratio test result (only when two models
#'   are supplied and `nested = TRUE`). Contains `ChiSq`, `df`, `p_value`.
#' - `evidence_ratios`: data.frame of pairwise Akaike-weight ratios (Model1,
#'   Model2, EvidenceRatio, LogEvidenceRatio, RatioReason). `NULL` when
#'   comparison is ineligible. Extremely small normalized weights can round to zero.
#' - `preferred`: named list with the preferred model label by each criterion.
#' - `comparison_basis`: list describing whether IC and LRT comparisons were
#'   considered comparable. Includes a conservative `nesting_review` plus
#'   `lrt_status` / `lrt_reason` so withheld LRTs are explicit rather than
#'   silently absent.
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("compare_mfrm", run.dontrun = TRUE)
#' \dontrun{
#' # Load one dataset for both models
#' library(mfrmr)
#' toy <- load_mfrmr_data("example_operational")
#'
#' # RSM: shared category thresholds
#' fit_rsm <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   model = "RSM"
#' )
#'
#' # PCM: separate category thresholds for each criterion
#' fit_pcm <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   model = "PCM",
#'   step_facet = "Criterion"
#' )
#'
#' # Check that the fitted models can be compared
#' comparison <- compare_mfrm(fit_rsm, fit_pcm, labels = c("RSM", "PCM"))
#' comparison$table[, c("Label", "Converged", "ICComparable", "ICSelectable")]
#'
#' # Smaller AIC/BIC indicates better relative support among these models
#' comparison$table[, c("Label", "AIC", "Delta_AIC", "BIC", "Delta_BIC")]
#' # Delta is the difference from the lowest value of that criterion
#' # For close or consequential comparisons, check a denser shared grid with
#' # mml_quadrature_sensitivity() before choosing a model
#' # For a PCM/GPCM equal-slope test, see the matched-population example in
#' # vignette("mfrmr-gpcm-scope"); request compare_mfrm(..., nested = TRUE).
#' }
#' @section References:
#' - Burnham, K. P., & Anderson, D. R. (2002). *Model selection and
#'   multimodel inference: A practical information-theoretic
#'   approach* (2nd ed.). Springer.
#' - Akaike, H. (1974). A new look at the statistical model
#'   identification. *IEEE Transactions on Automatic Control,
#'   19*(6), 716-723.
#' - Schwarz, G. (1978). Estimating the dimension of a model.
#'   *Annals of Statistics, 6*(2), 461-464.
#' - Sclove, S. L. (1987). Application of model-selection criteria to some
#'   problems in multivariate analysis. *Psychometrika, 52*(3), 333-343.
#' @export
#' @order 1
compare_mfrm <- function(..., labels = NULL, warn_constraints = TRUE, nested = FALSE,
    response_diagnostics = NULL, person_scores = NULL) {
  if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
    stop("`nested` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  fits <- list(...)
  for (fit in fits) stop_if_product_slopes(fit, "compare_mfrm()")
  if (any(vapply(fits, mfrm_extended_fit, logical(1)))) return(mfrm_compare_extended(fits, labels, nested, response_diagnostics, person_scores))
  if (!is.null(response_diagnostics) || !is.null(person_scores)) stop("Predictive or Person comparison currently requires one ordinary RSM MML and one testlet/shared-rater fit.", call. = FALSE)
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

  # -- constraint and score-coding compatibility check --
  comparison_signatures <- lapply(fits, normalize_compare_signature)
  constraint_component_names <- c(
    "person", "facets", "score", "rating_min", "rating_max", "score_map",
    "noncenter_facet", "dummy_facets", "positive_facets", "anchors",
    "group_anchors", "centering_scope"
  )
  constraint_components <- vapply(constraint_component_names, function(nm) {
    all(vapply(comparison_signatures[-1], function(sig) {
      same_signature_component(comparison_signatures[[1]][[nm]], sig[[nm]])
    }, logical(1)))
  }, logical(1))
  same_constraint_basis <- all(constraint_components)

  if (isTRUE(warn_constraints) && !same_constraint_basis) {
    mismatch <- names(constraint_components)[!constraint_components]
    if ("noncenter_facet" %in% mismatch) {
      noncenter <- vapply(comparison_signatures, function(sig) {
        as.character(sig$noncenter_facet[1])
      }, character(1))
      warning(
        "Models use different centering constraints (",
        paste(unique(noncenter), collapse = ", "),
        "). Information-criterion ranking is suppressed.",
        call. = FALSE
      )
    }
    if ("dummy_facets" %in% mismatch) {
      dummy_sigs <- vapply(comparison_signatures, function(sig) {
        paste(sig$dummy_facets, collapse = ",")
      }, character(1))
      warning(
        "Models use different dummy-facet constraints (",
        paste(unique(dummy_sigs), collapse = " vs "),
        "). Information-criterion ranking is suppressed.",
        call. = FALSE
      )
    }
    other_mismatch <- setdiff(mismatch, c("noncenter_facet", "dummy_facets"))
    if (length(other_mismatch) > 0L) {
      warning(
        "Models differ in likelihood-basis coding or constraint settings: ",
        paste(other_mismatch, collapse = ", "),
        ". Information-criterion ranking is suppressed.",
        call. = FALSE
      )
    }
  }

  # -- extract summary statistics --
  rows <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    convergence <- mfrm_convergence_state(f)
    contract <- mfrm_extract_fit_ic_contract(f)
    ic_fit <- mfrm_ic_fit_check(f, contract)
    readiness <- mfrmr_get_readiness_record(f)$fit
    tibble::add_column(
      contract,
      Label = labels[i],
      Converged = convergence$code_converged,
      InferenceReady = convergence$inference_ready,
      ICFitEligible = ic_fit$eligible,
      ICFitBasis = ic_fit$basis,
      ICFitReview = ic_fit$review,
      ICFitCaution = ic_fit$caution %||% "",
      InferenceReview = if (convergence$inference_ready) "" else
        as.character(mfrm_fit_decision_summary(readiness)$Why[1]),
      ReadinessReasonCodes = as.character(readiness$ReasonCodes[1]),
      ConvergenceSeverity = convergence$severity,
      .before = 1L
    )
  })
  tbl <- bind_rows(rows)
  ic_cautions <- unique(tbl$ICFitCaution[nzchar(tbl$ICFitCaution)])
  if (length(ic_cautions)) warning(paste(ic_cautions, collapse = " "), call. = FALSE)

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
  conv_vals <- tbl$InferenceReady
  all_converged <- length(conv_vals) > 0 && all(!is.na(conv_vals) & as.logical(conv_vals))
  all_ic_fit_eligible <- nrow(tbl) > 0L && all(!is.na(tbl$ICFitEligible) & tbl$ICFitEligible)
  ic_fit_review_reason <- if (!all_ic_fit_eligible) {
    pending <- which(is.na(tbl$ICFitEligible) | !tbl$ICFitEligible)
    paste0("IC solution checks are not satisfied: ",
           paste(paste0(labels[pending], ": ", tbl$ICFitReview[pending]), collapse = "; "))
  } else ""
  inference_review_reason <- if (!all_converged) {
    pending <- which(is.na(conv_vals) | !as.logical(conv_vals))
    details <- vapply(pending, function(i) {
      paste0(labels[i], ": ", sub("[.]$", "", tbl$InferenceReview[i]))
    }, character(1))
    paste0("Inference readiness is not satisfied: ",
           paste(details, collapse = "; "), ".")
  } else {
    ""
  }
  all_current_contract <- nrow(tbl) > 0L &&
    all(!is.na(tbl$ICContractVersion) &
          tbl$ICContractVersion == mfrm_ic_contract_version())
  all_stored_consistent <- nrow(tbl) > 0L &&
    all(!is.na(tbl$StoredICConsistent) & tbl$StoredICConsistent)
  all_ic_eligible <- nrow(tbl) > 0L &&
    all(!is.na(tbl$ICEligible) & tbl$ICEligible)
  all_ic_selectable <- nrow(tbl) > 0L &&
    all(!is.na(tbl$ICSelectable) & tbl$ICSelectable) &&
    all(!is.na(tbl$ICIntegrationSelectable) & tbl$ICIntegrationSelectable)
  contract_vals <- tbl$ICContractVersion[!is.na(tbl$ICContractVersion)]
  same_contract <- length(contract_vals) > 0L &&
    length(unique(contract_vals)) == 1L &&
    identical(unique(contract_vals), mfrm_ic_contract_version())
  integration_vals <- tbl$IntegrationEvaluationId[
    !is.na(tbl$IntegrationEvaluationId) & nzchar(tbl$IntegrationEvaluationId)
  ]
  same_integration <- length(integration_vals) == nrow(tbl) &&
    length(unique(integration_vals)) == 1L
  same_formula <- length(unique(tbl$AICFormula)) == 1L &&
    length(unique(tbl$BICFormula)) == 1L &&
    length(unique(tbl$SABICFormula)) == 1L
  ic_comparable <- same_method && same_nobs && same_data && all_ic_fit_eligible &&
    all_mml && all_current_contract && all_stored_consistent &&
    all_ic_eligible && all_ic_selectable && same_contract && same_integration &&
    same_formula && same_constraint_basis
  sabic_comparable <- ic_comparable &&
    all(!is.na(tbl$SABICSelectable) & tbl$SABICSelectable)

  if (!isTRUE(nested)) {
    if (!same_method) {
      warning(
        "Models use different estimation methods (",
        paste(unique(method_vals), collapse = ", "),
        "). Canonical IC ranking, deltas, weights, ",
        "and automatic preferences are suppressed.",
        call. = FALSE
      )
    } else if (!all_mml) {
      warning(
        "Information-criterion ranking is limited to converged MML fits in this package. ",
        "Legacy descriptive values remain separately labelled for ",
        paste(unique(method_vals), collapse = ", "),
        " models, but deltas, weights, automatic preferences, and LRT were suppressed.",
        call. = FALSE
      )
    }

    # -- warn if observation counts differ --
    if (!same_nobs) {
      warning("Models were fit to different numbers of observations (",
              paste(obs_vals, collapse = ", "),
              "). Cross-sample deltas, weights, ",
              "and automatic preferences are suppressed.", call. = FALSE)
    }
    if (!same_data) {
      warning(
        "Models were not fit to the same prepared response data. ",
        "Deltas, weights, automatic preferences, and likelihood-ratio testing were suppressed.",
        call. = FALSE
      )
    }
    if (!all_ic_fit_eligible) {
      warning(
        ic_fit_review_reason, " ",
        "IC ranking, weights, and likelihood-ratio testing were suppressed.",
        call. = FALSE
      )
    }
    if (!all_current_contract) {
      warning(
        "At least one fit lacks current, complete information for comparing its information criteria. ",
        "LegacyAIC/LegacyBIC remain descriptive, but canonical ranking was suppressed; ",
        "refit under the current package before comparison.",
        call. = FALSE
      )
    } else if (!all_stored_consistent) {
      warning(
        "At least one stored information-criterion value or identity does not match ",
        "the retained objective, parameter vector, prepared data, or integration setting. ",
        "Canonical ranking was suppressed.",
        call. = FALSE
      )
    }
    if (all_mml && all_current_contract && !all_ic_eligible) {
      warning(
        "At least one MML fit cannot be ranked by information criteria. ",
        "Check for non-unit observation weights or incomplete likelihood information.",
        call. = FALSE
      )
    }
    if (all_mml && all_current_contract && all_ic_eligible &&
        !all_ic_selectable) {
      warning(
        "Information-criterion ranking is unavailable at the supplied ",
        "quadrature resolution. Raw criteria remain visible, but deltas, weights, preferences, and ",
        "LRT are suppressed below q = 31. Refit every candidate with a common ",
        "q >= 31 grid and check a denser grid when the decision is close.",
        call. = FALSE
      )
    }
    if (all_current_contract && all_mml && !same_integration) {
      warning(
        "Models do not share one integration-evaluation identity. Refit or reevaluate ",
        "all candidates with the same locked quadrature setting before IC ranking.",
        call. = FALSE
      )
    }
  }

  tbl$ICComparable <- ic_comparable
  tbl$SABICComparable <- sabic_comparable
  tbl$ConstraintComparable <- same_constraint_basis
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

  # -- Delta SABIC and adjusted-BIC candidate-set weights --
  if (sabic_comparable && any(is.finite(tbl$SABIC))) {
    min_sabic <- min(tbl$SABIC, na.rm = TRUE)
    tbl$Delta_SABIC <- tbl$SABIC - min_sabic
    raw_sw <- exp(-0.5 * tbl$Delta_SABIC)
    tbl$SABICWeight <- ifelse(
      is.finite(raw_sw),
      raw_sw / sum(raw_sw, na.rm = TRUE),
      NA_real_
    )
    preferred$SABIC <- tbl$Label[which.min(tbl$SABIC)]
  } else {
    tbl$Delta_SABIC <- NA_real_
    tbl$SABICWeight <- NA_real_
  }

  # -- likelihood-ratio test (two models only) --
  lrt <- NULL
  lrt_status <- if (isTRUE(nested)) "requested" else "not_requested"
  lrt_reason <- if (isTRUE(nested)) {
    "Likelihood-ratio test was requested but not yet evaluated."
  } else {
    "Likelihood-ratio test was not requested; set `nested = TRUE` only for supported nested comparisons."
  }
  nesting_review <- audit_compare_mfrm_nesting(fits, labels = labels)
  pcm_gpcm_regular <- isTRUE(nesting_review$eligible) &&
    identical(nesting_review$relation, "PCM_in_GPCM") && all_ic_fit_eligible
  if (isTRUE(nested)) {
    if (length(fits) != 2L) {
      lrt_status <- "not_computed"
      lrt_reason <- "Likelihood-ratio testing requires exactly two fitted models."
      warning("`nested = TRUE` was requested, but LRT requires exactly two fitted models. LRT was not computed.",
              call. = FALSE)
    } else if (!ic_comparable || (!all_converged && !pcm_gpcm_regular)) {
      lrt_status <- "not_computed"
      lrt_reason <- paste(
        "Not all comparison requirements are satisfied: formal MML likelihood basis,",
        "verified information criteria, common response data, constraints and",
        "integration settings, and adequate support for inference."
      )
      if (!all_ic_fit_eligible) lrt_reason <- paste(lrt_reason, ic_fit_review_reason)
      if (!all_converged && !pcm_gpcm_regular) lrt_reason <- paste(lrt_reason, inference_review_reason)
      warning(
        "`nested = TRUE` was requested. ", lrt_reason, " LRT was not computed.",
        call. = FALSE
      )
    } else if (!isTRUE(nesting_review$eligible)) {
      lrt_status <- "not_computed"
      lrt_reason <- paste0("Structural nesting review did not pass: ", nesting_review$reason)
      warning(
        "`nested = TRUE` was requested, but the structural nesting review did not pass. ",
        nesting_review$reason,
        " LRT was not computed.",
        call. = FALSE
      )
    } else {
      ll <- tbl$LogLik
      np <- tbl$npar
      if (!all(is.finite(ll)) || !all(is.finite(np))) {
        lrt_status <- "not_computed"
        lrt_reason <- "Log-likelihoods or parameter counts were non-finite."
        warning(
          "`nested = TRUE` was requested, but at least one log-likelihood or parameter count was non-finite. LRT was not computed.",
          call. = FALSE
        )
      } else if (np[1] == np[2]) {
        lrt_status <- "not_computed"
        lrt_reason <- "Compared models have the same number of free parameters."
        warning(
          "`nested = TRUE` was requested, but the compared models have the same number of free parameters. LRT was not computed.",
          call. = FALSE
        )
      } else {
        # User labels need not be unique. Identify the PCM null by model, not
        # by matching its display label or guessing from a possibly bad count.
        idx_simple <- if (pcm_gpcm_regular) which(tbl$Model == "PCM")[1L] else which.min(np)
        idx_complex <- if (pcm_gpcm_regular) which(tbl$Model == "GPCM")[1L] else which.max(np)
        chi_sq <- 2 * (ll[idx_complex] - ll[idx_simple])
        df_diff <- np[idx_complex] - np[idx_simple]
        df_verified <- is.null(nesting_review$df) || identical(as.numeric(df_diff), as.numeric(nesting_review$df))
        if (is.finite(chi_sq) && chi_sq >= 0 && df_diff > 0 && df_verified) {
          p_val <- stats::pchisq(chi_sq, df = df_diff, lower.tail = FALSE)
          lrt <- tibble(
            Simple   = tbl$Label[idx_simple],
            Complex  = tbl$Label[idx_complex],
            ChiSq    = chi_sq,
            df       = df_diff,
            p_value  = p_val,
            ReferenceDistribution = "asymptotic_chi_square",
            Interpretation = "nested_fit_evidence_not_practical_preference"
          )
          lrt_status <- "computed"
          lrt_reason <- paste(
            "Likelihood-ratio test computed after the structural nesting review passed;",
            "the chi-square reference is asymptotic and assumes independent Persons,",
            "a correctly specified population/response model and regular interior solutions.",
            "Check starting-value and integration sensitivity; the p-value is not an automatic practical model preference."
          )
        } else {
          lrt_status <- "not_computed"
          lrt_reason <- if (!is.finite(chi_sq) || !is.finite(df_diff)) {
            "Likelihood-ratio statistic or parameter-count difference was non-finite."
          } else if (!df_verified) {
            "The parameter-count difference does not equal the verified number of null restrictions."
          } else if (df_diff <= 0) {
            "Parameter-count difference was not positive."
          } else {
            "The more complex model had a lower log-likelihood, producing a negative likelihood-ratio statistic."
          }
          warning(
            "`nested = TRUE` was requested, but the LRT statistic was not interpretable. ",
            lrt_reason,
            " LRT was not computed.",
            call. = FALSE
          )
        }
      }
    }
  }

  if (isTRUE(nested) && length(ic_cautions)) {
    lrt_reason <- paste(lrt_reason, paste(ic_cautions, collapse = " "))
  }

  # -- evidence ratios (pairwise Akaike-weight ratios) --
  evidence_ratios <- NULL
  if (ic_comparable && "AkaikeWeight" %in% names(tbl) && nrow(tbl) >= 2) {
    er_rows <- list()
    for (i in 1:(nrow(tbl) - 1)) {
      for (j in (i + 1):nrow(tbl)) {
        log_er <- (tbl$AIC[j] - tbl$AIC[i]) / 2
        if (!is.finite(log_er)) log_er <- tbl$AIC[j] / 2 - tbl$AIC[i] / 2
        er <- exp(log_er)
        ratio_reason <- if (!is.finite(log_er)) "nonfinite_difference" else
          if (!is.finite(er)) "overflow" else
          if (er == 0) "underflow" else "available"
        if (ratio_reason != "available") er <- NA_real_
        er_rows[[length(er_rows) + 1]] <- tibble(
          Model1        = tbl$Label[i],
          Model2        = tbl$Label[j],
          EvidenceRatio = er,
          LogEvidenceRatio = log_er,
          RatioReason = ratio_reason
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
    all_inference_ready = all_converged,
    all_ic_fit_eligible = all_ic_fit_eligible,
    all_converged = all_converged,
    all_current_contract = all_current_contract,
    all_stored_ic_consistent = all_stored_consistent,
    all_ic_eligible = all_ic_eligible,
    all_ic_selectable = all_ic_selectable,
    same_ic_contract = same_contract,
    same_integration_evaluation = same_integration,
    same_formula_contract = same_formula,
    same_constraint_basis = same_constraint_basis,
    constraint_components = constraint_components,
    ic_comparable = ic_comparable,
    sabic_comparable = sabic_comparable,
    nested_requested = isTRUE(nested),
    nesting_review = nesting_review,
    lrt_status = lrt_status,
    lrt_reason = lrt_reason
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

print_mfrm_lrt <- function(lrt) {
  cat("\nLikelihood-ratio test:\n")
  p_text <- trimws(format.pval(lrt$p_value[1], digits = 4, eps = 1e-4))
  cat(sprintf("  Chi-sq = %.3f, df = %d, p %s%s\n",
              lrt$ChiSq[1], lrt$df[1],
              if (startsWith(p_text, "<")) "" else "= ", p_text))
  cat(sprintf("  %s vs %s\n", lrt$Simple[1], lrt$Complex[1]))
  cat("  Reference: asymptotic chi-square; small or sparse samples may need further calibration.\n")
  cat("  Interpretation: formal nested-fit evidence; practical gain and score utility require separate checks.\n")
  invisible(lrt)
}

#' @export
print.summary.mfrm_comparison <- function(x, ...) {
  cat("--- MFRM Model Comparison ---\n\n")

  # -- main comparison table --
  tbl <- x$table
  # Format weight columns with 4 decimal places for readability
  fmt_tbl <- as.data.frame(tbl)
  weight_cols <- intersect(
    c("AkaikeWeight", "BICWeight", "SABICWeight"),
    names(fmt_tbl)
  )
  for (wc in weight_cols) {
    fmt_tbl[[wc]] <- ifelse(is.na(fmt_tbl[[wc]]), NA_character_,
                            sprintf("%.4f", fmt_tbl[[wc]]))
  }
  delta_cols <- intersect(
    c("Delta_AIC", "Delta_BIC", "Delta_SABIC"),
    names(fmt_tbl)
  )
  for (dc in delta_cols) {
    fmt_tbl[[dc]] <- ifelse(is.na(fmt_tbl[[dc]]), NA_character_,
                            sprintf("%.2f", fmt_tbl[[dc]]))
  }
  display_cols <- c(
    "Label", "Model", "Method", "Persons", "Npar", "LogLik",
    "AIC", "BIC", "SABIC"
  )
  if (isTRUE(x$comparison_basis$ic_comparable)) {
    display_cols <- c(
      display_cols, "Delta_AIC", "AkaikeWeight", "Delta_BIC", "BICWeight",
      "Delta_SABIC", "SABICWeight"
    )
  }
  display_cols <- intersect(display_cols, names(fmt_tbl))
  print(fmt_tbl[, display_cols, drop = FALSE], row.names = FALSE)
  print_wrapped_line("Full statistics, eligibility and numerical settings remain in `$table`.")

  if (!isTRUE(x$comparison_basis$ic_comparable)) {
    cat("\nInformation-criterion ranking is unavailable. Displayed criteria are descriptive only.\n")
    reasons <- c(
      all_mml = "Ranking requires MML fits; JML criteria are descriptive only.",
      same_data = "Models must use the same response data.",
      same_constraint_basis = "Models must use compatible score coding and constraints.",
      all_current_contract = "Stored comparison information is incomplete or outdated; refit before comparing.",
      all_stored_ic_consistent = "Stored criteria do not match the fitted values or numerical settings.",
      all_ic_eligible = "One or more likelihoods are ineligible, for example because of non-unit observation weights.",
      all_ic_selectable = "One or more fits lack eligible information criteria or sufficient quadrature resolution.",
      same_integration_evaluation = "Models must use the same integration settings.",
      all_ic_fit_eligible = "One or more models lack adequate numerical or identification checks for IC comparison."
    )
    for (field in names(reasons)) {
      if (field %in% c("all_ic_eligible", "all_ic_selectable") &&
          !isTRUE(x$comparison_basis$all_mml)) next
      if (identical(x$comparison_basis[[field]], FALSE)) print_wrapped_line(reasons[[field]])
    }
    review_column <- if ("ICFitReview" %in% names(tbl)) "ICFitReview" else "InferenceReview"
    if (review_column %in% names(tbl)) {
      for (i in which(!is.na(tbl[[review_column]]) & nzchar(tbl[[review_column]]))) {
        print_wrapped_line(paste0(tbl$Label[i], ": ", tbl[[review_column]][i]))
      }
    }
  } else if (!isTRUE(x$comparison_basis$sabic_comparable)) {
    cat("\nSABIC is displayed for sensitivity only; automatic SABIC selection is disabled\n")
    cat("at 22 or fewer Persons.\n")
  }

  if (!is.null(x$lrt)) {
    print_mfrm_lrt(x$lrt)
  } else if (isTRUE(x$comparison_basis$nested_requested)) {
    cat("\nLikelihood-ratio test was not reported.\n")
    if (!is.null(x$comparison_basis$lrt_status) &&
        !is.null(x$comparison_basis$lrt_reason)) {
      if (nrow(tbl) != 2L || isTRUE(x$comparison_basis$ic_comparable)) {
        print_wrapped_line(paste0("Reason: ", x$comparison_basis$lrt_reason))
      }
    }
    review <- x$comparison_basis$nesting_review %||% list()
    if (!is.null(review$reason) && nzchar(review$reason)) {
      cat("  Nesting review:", review$reason, "\n")
    }
  }

  # -- evidence ratios --
  for (caution in unique(tbl$ICFitCaution[nzchar(tbl$ICFitCaution)])) {
    print_wrapped_line(paste0("Caution: ", caution))
  }
  if (!is.null(x$evidence_ratios) && nrow(x$evidence_ratios) > 0) {
    cat("\nRelative candidate-set ratios (Akaike weights; not model probabilities):\n")
    er <- x$evidence_ratios
    for (k in seq_len(nrow(er))) {
      ratio_str <- if (is.finite(er$EvidenceRatio[k])) {
        format(signif(er$EvidenceRatio[k], 3), trim = TRUE)
      } else if (all(c("RatioReason", "LogEvidenceRatio") %in% names(er))) {
        sprintf("NA (%s; log ratio = %s)", er$RatioReason[k],
                format(signif(er$LogEvidenceRatio[k], 4), trim = TRUE))
      } else {
        "NA"
      }
      cat(sprintf("  %s / %s = %s\n",
                  er$Model1[k], er$Model2[k], ratio_str))
    }
  }

  if (length(x$preferred) > 0) {
    cat("\nMinimum criterion within this candidate set:\n")
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
