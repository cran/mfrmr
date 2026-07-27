#' Fit many-facet ordered-response models with a flexible number of facets
#'
#' This is the package entry point. It wraps `mfrm_estimate()` and defaults to
#' `method = "MML"`. Any number of facet columns can be supplied via `facets`.
#' The `RSM` / `PCM` branches are the package's many-facet Rasch-family
#' reference route; the bounded `GPCM` branch is available where explicitly
#' documented.
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
#' @param keep_original Logical. `FALSE` (the current default) collapses
#'   non-consecutive observed categories to a contiguous internal scale and
#'   records the mapping in `fit$prep$score_map` (the downstream Count = 0
#'   rows are consequently absent). `TRUE` preserves the declared scale so
#'   unused intermediate categories remain visible in
#'   [rating_scale_table()] and APA outputs, which is recommended for
#'   publication reporting.
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
#' @param model `"RSM"`, `"PCM"`, or bounded `"GPCM"`.
#' @param method `"MML"` (default) or `"JML"`. `"JMLE"` is accepted as a
#'   backward-compatible alias for the same joint-maximum-likelihood path.
#' @param step_facet Facet whose levels receive separate step parameters in
#'   `PCM` and bounded `GPCM`. Supply it explicitly for a final analysis. If it
#'   is omitted for `PCM`, mfrmr uses a unique item-like facet name (for
#'   example, `Item`, `Task`, or `Criterion`) when available; otherwise it
#'   retains the first-facet fallback with a warning. `GPCM` always requires an
#'   explicit value. This argument is not used by `RSM`, which has one shared
#'   set of rating-scale thresholds.
#' @param slope_facet Slope facet for the bounded `GPCM` branch. mfrmr
#'   requires `slope_facet == step_facet` and uses a
#'   positive-slope identification convention on the log scale with geometric
#'   mean discrimination fixed to 1.
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
#' @param anchors Optional anchor table.
#' @param group_anchors Optional group-anchor table.
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#' @param positive_facets Facets with positive orientation.
#' @param anchor_policy How to handle anchor-review issues: `"warn"` (default),
#'   `"error"`, or `"silent"`.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   anchor-review recommendations.
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in anchor-review recommendations.
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in anchor-review recommendations.
#' @param quad_points Integer number of Gauss-Hermite quadrature points
#'   used for MML integration over the person distribution. The default is
#'   `31`. Useful accuracy/runtime settings are:
#'   \tabular{ll}{
#'     `7`  \tab lightweight exploratory run; helpers such as
#'                 [predict_mfrm_population()] and
#'                 [reference_case_benchmark()] use this value. \cr
#'     `15` \tab intermediate analysis when runtime matters. \cr
#'     `31` \tab package default and a starting point for final analysis. \cr
#'     `61+` \tab sensitivity analysis for narrow score distributions or
#'                 demanding numerical comparisons.
#'   }
#'   Quadrature adequacy depends on the fitted distribution and score support.
#'   When substantive conclusions are sensitive, compare results under a
#'   denser rule and report the setting used.
#' @param maxit Computational ceiling on optimizer iterations. The default is
#'   `400`. This is not a convergence criterion or a model-selection control:
#'   a fit that reaches the ceiling remains non-ready until the common
#'   convergence and terminal-gradient checks pass. Smaller values used in
#'   executable examples shorten package checks and should not be copied into
#'   a final analysis without an explicit computational protocol.
#' @param reltol Portable tolerance setting for the initial optimizer stage.
#'   The default is `1e-9`. For BFGS this is passed as `reltol`; for L-BFGS-B
#'   it is mapped to `factr` and `pgtol`, whose actual values are recorded in
#'   the fit. When this setting is at least as strict as the public default
#'   (`reltol <= 1e-9`), optimizer code zero followed by a failed common
#'   terminal-gradient review triggers a bounded warm-started polish ladder.
#'   The best non-worsening stage under the recorded selection rule is
#'   retained. Requested and selected-stage settings remain in `fit$summary`,
#'   and the complete stage history remains in `fit$opt$optimizer_polish`.
#' @param optimizer Direct-optimization method. `"auto"` (default) uses the
#'   limited-memory `"L-BFGS-B"` method for MML and for larger JML parameter
#'   vectors (at least 200 free parameters), while retaining BFGS for smaller
#'   JML fits. Use `"BFGS"` or `"L-BFGS-B"` to request one method explicitly.
#'   The method actually used is recorded in `fit$summary$OptimizerMethod`.
#'   For L-BFGS-B, inspect `OptimizerFactr` and `OptimizerPgtol` rather than
#'   interpreting `EffectiveReltol` as a native [stats::optim()] control.
#' @param mml_engine MML optimization engine for `method = "MML"`:
#'   `"direct"` (default) uses the selected direct optimizer on the marginal
#'   log-likelihood,
#'   `"em"` uses an EM loop for `RSM` / `PCM` with `population = NULL`, and
#'   `"hybrid"` uses EM as a warm start before the direct optimizer. Unsupported
#'   combinations currently fall back to `"direct"` and record that fallback in
#'   `fit$summary`. Direct, hybrid, and EM engines all require the common
#'   terminal-gradient gate for `InferenceReady`; EM relative log-likelihood
#'   convergence alone does not establish numerical readiness.
#' @param population_formula Optional one-sided formula for a person-level
#'   latent-regression population model, for example `~ grade + ses`. Latent
#'   regression is implemented only for
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
#' @param facet_shrinkage Character. `"none"` (default) keeps the unshrunk
#'   fixed-effects estimates. `"empirical_bayes"` applies a post-hoc
#'   James-Stein / empirical-Bayes shrinkage to each non-person facet
#'   (Efron & Morris, 1973); `fit$facets$others` gains `ShrunkEstimate`,
#'   `ShrunkSE`, and `ShrinkageFactor` columns, and
#'   `fit$shrinkage_report` records the per-facet prior variance and
#'   effective degrees of freedom. `"laplace"` is retained as a compatibility
#'   alias for `"empirical_bayes"`; it does not fit a penalized likelihood.
#' @param facet_prior_sd Optional numeric scalar. When supplied, the
#'   shrinkage prior variance is fixed at `facet_prior_sd^2` instead
#'   of being estimated by method of moments. Useful for eliciting a
#'   prior from domain knowledge or a previous fit.
#' @param shrink_person Logical. When `TRUE` and `facet_shrinkage` is
#'   active, the same empirical-Bayes shrinkage is applied to
#'   `fit$facets$person`. Default `FALSE`, since MML already integrates
#'   over an N(0, 1) prior on theta; the option mainly benefits JML.
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
#'   every `every_iter` outer EM iterations using `saveRDS()`.
#'   If the file already exists when the fit starts, the engine
#'   resumes from the recorded iteration. Only the EM engine
#'   (`mml_engine = "em"` or the EM warm-start step of
#'   `mml_engine = "hybrid"`) honours the checkpoint; the direct
#'   `optim()` engine ignores it. Use this to make long MML EM
#'   fits crash-resilient on shared compute environments.
#'
#' @details
#' Data must be in **long format** (one row per observed rating event).
#' Exact duplicate Person-by-facet combinations are retained, warned once, and
#' propagated as a Data review state. They are not treated as independent
#' replication evidence. A legitimate re-rating or replicated scoring event
#' should be represented by an event, occasion, or other distinguishing facet
#' before fitting.
#'
#' @section Model:
#' `fit_mfrm()` estimates many-facet ordered-response models. The `RSM` and
#' `PCM` branches follow the many-facet Rasch-family tradition (Linacre, 1989);
#' the bounded `GPCM` branch extends the partial-credit kernel with estimated
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
#' With bounded `model = "GPCM"`, the adjacent-category kernel is multiplied by
#' a positive slope for the designated slope-facet level:
#'
#' \deqn{\ln\frac{P(X_{nij} = k)}{P(X_{nij} = k-1)} =
#'   \alpha_g(\eta - \tau_{g,k}),\quad \alpha_g > 0.}
#'
#' The current implementation requires `slope_facet == step_facet` and
#' identifies slopes by a sum-to-zero constraint on log slopes, so their
#' geometric mean is 1.
#'
#' With only two ordered categories (\eqn{K = 1}), the `RSM`/`PCM`
#' branch reduces to the usual binary Rasch logit for the single category
#' boundary:
#'
#' \deqn{\ln\frac{P(X_{n\cdot} = 1)}{P(X_{n\cdot} = 0)} = \eta - \tau_1}
#'
#' Bounded `GPCM` uses the slope-scaled counterpart
#' \eqn{\alpha_g(\eta - \tau_{g,1})}.
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
#' Bounded `GPCM` is supported as an alternative when users explicitly accept
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
#' in the score support.
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
#' @section Fixed effects assumption (facets have no prior):
#' `fit_mfrm()` follows the Linacre (1989) many-facet Rasch specification:
#' person ability is integrated out under a `N(0, 1)` prior (or under the
#' `N(X\beta, \sigma^2)` latent-regression population model when
#' `population_formula` is supplied), but every facet parameter
#' (`Rater`, `Criterion`, `Task`, ...) is estimated as a fixed effect
#' identified by a sum-to-zero constraint. There is no hierarchical
#' prior, no shrinkage, and no variance component for the facets.
#'
#' Practical implication: when a facet has very few observed levels
#' (for example 3 raters) or some of its levels have very few ratings
#' (for example 5 ratings per rater), the fixed-effect estimates retain
#' wide SEs, and extreme estimates are not pulled toward the facet
#' mean. Jones and Wind (2018) note that rater estimates in particular
#' are "more sensitive to link reductions" than examinee or task
#' estimates. For a publication-workflow review of this, use:
#'
#' - [`facet_small_sample_review()`] for per-level N and SE bands against
#'   Linacre (1994) sample-size guidelines.
#' - [`detect_facet_nesting()`] and
#'   [`analyze_hierarchical_structure()`] when raters are nested in
#'   regions, schools, or other strata that the additive fixed-effects
#'   MFRM cannot partition out.
#' - [`compute_facet_icc()`] and
#'   [`compute_facet_design_effect()`] for descriptive variance-
#'   component summaries based on `lme4` (optional).
#'
#' `fit$summary$FacetSampleSizeFlag` summarizes the worst Linacre band
#' across non-person facet levels (`"sparse"` < 10, `"marginal"` < 30,
#' `"standard"` < 50, `"strong"` >= 50).
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
#' that correction (no `bias_correction` argument exists). The JML
#' branch also does not produce a profile-likelihood Hessian for the
#' structural parameters: SEs reported under JML are observation-table
#' approximations (\eqn{1/\sqrt{\sum \mathrm{Var}(X_{pi})}}) and are
#' marked as exploratory in the diagnostics output.
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
#' - For supported Rasch-family formulations, conditional maximum likelihood
#'   is a distribution-free alternative that conditions out person parameters.
#'   A third-party CML fit can be imported from `eRm` with
#'   [`import_erm_fit()`].
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
#'

#' This is ordered binary support, not a separate nominal-response model.
#' In `PCM`, a binary fit still uses one threshold per `step_facet` level on
#' the shared observed-score scale.
#'
#' Supported model/estimation combinations:
#' - `model = "RSM"` with `method = "MML"` or `"JML"/"JMLE"`
#' - `model = "PCM"` with a designated `step_facet` (defaults to first facet)
#' - `facet_interactions` with `model = "RSM"` or `"PCM"` for explicit
#'   two-way non-person facet interactions
#' - `model = "GPCM"` is currently implemented only for the narrow bounded
#'   branch with `slope_facet == step_facet`; `MML` and `JML` fitting, core
#'   summaries, fixed-calibration posterior scoring, [compute_information()],
#'   Wright/pathway/CCC fit plots, [diagnose_mfrm()], residual-PCA follow-up,
#'   [interrater_agreement_table()], [unexpected_response_table()],
#'   [displacement_table()], [measurable_summary_table()],
#'   [rating_scale_table()], [facet_quality_dashboard()],
#'   [reporting_checklist()], [category_structure_report()],
#'   [category_curves_report()], and graph/scorefile
#'   [facets_output_file_bundle()] routes are available with score-side
#'   caveats. Direct simulation
#'   specifications and data generation are also supported through
#'   [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()], and
#'   [simulate_mfrm_data()] when the slope-aware generator contract is stored
#'   explicitly; direct recovery checks are available through
#'   [evaluate_mfrm_recovery()] and [assess_mfrm_recovery()]. Slope-aware
#'   [fair_average_table()] and [estimate_bias()] are available with their
#'   documented caveats. Role-based design evaluation, population forecasting,
#'   diagnostic-screening, and signal-detection helpers are available as
#'   caveated sensitivity evidence. Full FACETS-style score-side contract
#'   review, posterior predictive checks, and MCMC estimation are not available
#'   for bounded `GPCM`. Use
#'   [gpcm_capability_matrix()] as the formal boundary statement for the
#'   current `GPCM` scope.
#'
#' Latent-regression status:
#' - `population_formula = NULL` keeps the standard unconditional `MML` / `JML`
#'   behavior.
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
#' @section Latent-regression standard-error caveat:
#' `summary(fit)$population_coefficients` reports point estimates of
#' \eqn{\hat{\boldsymbol{\beta}}} and \eqn{\hat{\sigma}^2} only. mfrmr does
#' **not** currently compute standard errors, confidence intervals, or
#' asymptotic z / Wald statistics for the population-model parameters: no
#' Hessian on \eqn{(\boldsymbol{\beta}, \log\sigma^2)} is extracted from the
#' marginal log-likelihood, and no `vcov()` method is exposed for these
#' coefficients. Treat the coefficient table as point estimates suitable
#' for descriptive reporting; **do not** quote \eqn{\hat{\beta}_j \pm 1.96
#' \cdot \mathrm{SE}} bounds because the SE column is not provided. A
#' marginal-Hessian-based SE for \eqn{(\boldsymbol{\beta}, \sigma^2)} is not
#' available from this function.
#'
#' Identification: the latent-regression intercept is identifiable only
#' under the default `noncenter_facet = "Person"` (which sum-to-zero-
#' centers all non-Person facets). `fit_mfrm()` therefore rejects an active
#' latent-regression model with a different `noncenter_facet` rather than
#' returning a confounded intercept.
#'
#' Anchor inputs are optional:
#' - `anchors` should contain facet/level/fixed-value information.
#' - `group_anchors` should contain facet/level/group/group-value information.
#' Both are normalized internally, so column names can be flexible
#' (`facet`, `level`, `anchor`, `group`, `groupvalue`, etc.).
#'
#' Anchor review behavior:
#' - `fit_mfrm()` automatically runs an anchor review.
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
#' @section Choosing maxit without result-driven tuning:
#' Treat `maxit` as a predeclared computational budget, not as a value to tune
#' until preferred estimates appear.
#'
#' 1. Choose the model, estimation method, optimizer, tolerance, quadrature
#'    rule, and initial `maxit` before examining coefficient or fit results.
#'    The default `maxit = 400` is the package starting point for an analysis.
#' 2. Use estimates substantively only when `Converged` and `InferenceReady`
#'    are both `TRUE` and the Numerical row of `summary(fit)$readiness` is
#'    `pass`. Optimizer code zero alone is insufficient.
#' 3. If `ConvergenceStatus == "iteration_limit"`, keep that fit review-only.
#'    Refit the same data, model, method, anchors, optimizer, tolerance, and
#'    quadrature rule with the next ceiling in a prespecified sequence, such as
#'    400, 800, then 1600. Do not choose among runs by coefficient size,
#'    statistical significance, fit statistics, or agreement with an expected
#'    answer.
#' 4. The first run in that sequence that clears the numerical gate becomes
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
#' - `quad_points = 7` is a lightweight setting for quick iteration.
#' - `quad_points = 15` is an intermediate option when runtime matters.
#' - `quad_points = 31` is the package default and a suitable starting point
#'   for a final analysis; always review convergence and, when conclusions are
#'   sensitive, compare a denser quadrature rule.
#' - `quad_points = 61` (or higher) supports sensitivity checks on narrow score
#'   distributions at additional computational cost.
#' - `mml_engine = "direct"` remains the most stable general-purpose path.
#' - `mml_engine = "em"` or `"hybrid"` currently target `RSM` / `PCM` fits
#'   without a latent-regression population model.
#' - Benchmark your own workload before using `mml_engine = "em"` or
#'   `"hybrid"` for final reporting; `direct` remains the documented default
#'   when you have not compared engines for your data.
#' - When a direct code-zero stage stops ahead of the terminal-gradient gate,
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
#'    helpers, [compute_information()], direct simulation/recovery helpers,
#'    [fair_average_table()], and [estimate_bias()] with their documented
#'    caveats. Use [gpcm_capability_matrix()] to confirm which helper families
#'    are currently supported, caveated, blocked, or deferred.
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
#' - Myford, C. M., & Wolfe, E. W. (2003). Detecting and measuring rater
#'   effects using many-facet Rasch measurement: Part I. *Journal of Applied
#'   Measurement*, 4(4), 386-422.
#' - Myford, C. M., & Wolfe, E. W. (2004). Detecting and measuring rater
#'   effects using many-facet Rasch measurement: Part II. *Journal of Applied
#'   Measurement*, 5(2), 189-227.
#' - Muraki, E. (1992). *A generalized partial credit model: Application of an
#'   EM algorithm*. Applied Psychological Measurement, 16(2), 159-176.
#' - Robitzsch, A., & Steinfeld, J. (2018). *Item response models for human
#'   ratings: Overview, estimation methods, and implementation in R*.
#'   Psychological Test and Assessment Modeling, 60(1), 101-139.
#'
#' @return
#' An object of class `mfrm_fit` (named list) with:
#' - `summary`: one-row model summary (`LogLik`, `AIC`, `BIC`, convergence),
#'   including user-facing `Method`, engine-facing `MethodUsed`, MML-engine
#'   fields, terminal-gradient readiness, requested/selected-stage tolerance
#'   settings, and the actual L-BFGS-B `OptimizerFactr` / `OptimizerPgtol`
#'   controls when applicable
#' - `facets$person`: person estimates (`Estimate`; plus `SD` for MML)
#' - `facets$others`: facet-level estimates for each facet
#' - `steps`: estimated threshold/step parameters as a one-row-per-step
#'   `tibble` with `Estimate`. Bare fits keep this table as point estimates.
#'   `diagnose_mfrm()` exposes MML observed-information step uncertainty in
#'   `diagnostics$parameter_uncertainty$steps`; when
#'   `attach_diagnostics = TRUE`, those `SE`, confidence-limit, and status
#'   columns are attached to `fit$steps` when the Hessian is available.
#'   For step-structure quality, also use the step-collapse and disordering
#'   warnings from [diagnose_mfrm()] and [category_structure_report()].
#' - `slopes`: estimated discrimination parameters for `GPCM` fits as
#'   a one-row-per-slope-element `tibble` with `LogEstimate` and
#'   `Estimate`. Bare fits keep this table as point estimates. For MML
#'   bounded-`GPCM` fits, `diagnose_mfrm()` exposes log-slope SEs plus
#'   positive-scale delta-method SEs and confidence limits in
#'   `diagnostics$parameter_uncertainty$slopes`; when
#'   `attach_diagnostics = TRUE`, those columns are attached to
#'   `fit$slopes` when the Hessian is available. The identification
#'   convention pins the geometric mean of slopes at 1.
#' - `interactions`: model-estimated facet interaction effects and metadata
#'   when `facet_interactions` is supplied
#' - `population`: population-model metadata. Ordinary fits keep an inactive
#'   record (`active = FALSE`, `posterior_basis = "legacy_mml"`). Active
#'   latent-regression fits store the fitted design matrix, regression
#'   coefficients, residual variance, omission review, the complete-case
#'   estimation table (`person_table`), and the observed-person-aligned
#'   replay/export provenance table retained before complete-case omission
#'   (`person_table_replay`), plus stored categorical `xlevels` / `contrasts`
#'   for model-matrix replay and scoring, together with
#'   `posterior_basis = "population_model"`.
#' - `data_review`: pre-fit Data, Design, Stability, and Reporting readiness
#'   evidence propagated into summaries and plot-interpretation gates
#' - `config`: resolved model configuration used for estimation, including
#'   `config$anchor_review` and the recorded estimation controls
#' - `prep`: preprocessed data/level metadata
#' - `opt`: optimizer result augmented with `optimizer_diagnostics`, the
#'   complete `optimizer_polish` stage history, method-selection metadata, and
#'   evaluation-cache counters. For direct fitting, its core fields originate
#'   from [stats::optim()]; EM additionally records engine-specific diagnostics.
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()], [build_apa_outputs()],
#'   [gpcm_capability_matrix], [mfrmr_workflow_methods],
#'   [mfrmr_reporting_and_apa]
#' @examples
#' # Lightweight executable mechanics example on the connected teaching data.
#' # The small quadrature grid keeps CRAN example time short; the tighter
#' # portable tolerance setting keeps this reduced example numerically stable.
#' # Use the documented default grid and a sensitivity check for final work.
#' toy <- load_mfrmr_data("example_operational")
#' fit_quick <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", model = "RSM", quad_points = 7, maxit = 30,
#'   reltol = 1e-11
#' )
#' fit_quick$summary[, c(
#'   "Model", "Method", "N", "Converged", "InferenceReady",
#'   "ConvergenceSeverity"
#' )]
#'
#' \donttest{
#' # Full run with the package default MML estimator. This route integrates
#' # person parameters under an N(0, 1) population model, so its reporting
#' # value depends on the response-model and population assumptions. The
#' # default `quad_points = 31` is a practical starting value; compare a
#' # larger grid when quadrature sensitivity matters.
#' fit <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   model = "RSM",
#'   quad_points = 31
#' )
#' fit$summary
#' s_fit <- summary(fit)
#' s_fit$overview[, c("Model", "Method", "Converged", "InferenceReady",
#'                    "ConvergenceSeverity")]
#' # `InferenceReady = FALSE` is a numerical stop signal. A TRUE value only
#' # clears the package's optimizer review; model specification, design,
#' # identification, and inferential assumptions still require review.
#' s_fit$person_overview
#' # Compare the person distribution with the facet and step locations. The
#' # scale identification does not create universal targeting thresholds.
#' s_fit$targeting
#' # Interpret targeting magnitude against the intended population and score
#' # use rather than a universal pass/fail cutoff.
#' p_fit <- plot(fit, draw = FALSE)
#' p_fit$name
#' head(p_fit$data$locations)
#' # The bare plot route is the native Wright map and includes available
#' # facet uncertainty. Use plot(fit, type = "bundle") for the three-plot
#' # Wright/pathway/category overview.
#'
#' # JML is a distinct fixed-person-effects route, not a drop-in speed setting:
#' fit_jml <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   model = "RSM"
#' )
#' summary(fit_jml)$overview[, c(
#'   "Model", "Method", "Converged", "InferenceReady",
#'   "ConvergenceSeverity"
#' )]
#'
#' # Latent regression (MML only) uses person-level background variables:
#' person_tbl <- unique(toy[c("Person")])
#' person_tbl$Grade <- seq_len(nrow(person_tbl))
#' person_tbl$Group <- rep(c("A", "B"), length.out = nrow(person_tbl))
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
#'   maxit = 30
#' )
#' fit_binary$summary[, c("Model", "Categories", "Converged")]
#'
#' # Next steps after fitting:
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' chk <- reporting_checklist(fit, diagnostics = diag)
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
                     missing_codes = NULL,
                     model = c("RSM", "PCM", "GPCM"),
                     method = c("MML", "JML", "JMLE"),
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
                     checkpoint = NULL) {
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
  if (!is.logical(keep_original) || length(keep_original) != 1L ||
      is.na(keep_original)) {
    stop("`keep_original` must be a single logical value (TRUE or FALSE).",
         call. = FALSE)
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
  if (!is.null(slope_facet) && (!is.character(slope_facet) || length(slope_facet) != 1 || !nzchar(slope_facet))) {
    stop("`slope_facet` must be NULL or a single non-empty character string naming a facet column.",
         call. = FALSE)
  }

  model <- toupper(match.arg(model))
  method_input <- public_mfrm_method_label(toupper(match.arg(method)))
  method <- method_input
  optimizer <- normalize_mfrm_optimizer(match.arg(optimizer))
  mml_engine <- tolower(match.arg(mml_engine))
  interaction_policy <- tolower(match.arg(interaction_policy))
  anchor_policy <- tolower(match.arg(anchor_policy))
  population_policy <- tolower(match.arg(population_policy))
  facet_shrinkage <- match.arg(facet_shrinkage)
  if (!is.null(facet_prior_sd)) {
    if (!is.numeric(facet_prior_sd) || length(facet_prior_sd) != 1L ||
        !is.finite(facet_prior_sd) || facet_prior_sd < 0) {
      stop("`facet_prior_sd` must be NULL or a single non-negative finite number.",
           call. = FALSE)
    }
  }
  if (!is.logical(shrink_person) || length(shrink_person) != 1L ||
      is.na(shrink_person)) {
    stop("`shrink_person` must be a single logical value.", call. = FALSE)
  }

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
      stop("Latent-regression estimation requires `method = 'MML'`. ",
           "The requested population model is not available with JML.",
           call. = FALSE)
    }
    if (!identical(as.character(noncenter_facet[1]), "Person")) {
      stop(
        "Latent-regression identification requires `noncenter_facet = \"Person\"` ",
        "so all non-person facets remain centered and the population-model ",
        "intercept is identified.",
        call. = FALSE
      )
    }
  }
  estimation_data <- data
  if (isTRUE(population$active) && length(population$included_persons) > 0) {
    estimation_mask <- as.character(data[[person]]) %in% population$included_persons
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
    dummy_facets = dummy_facets
  )

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
    checkpoint = checkpoint
  )

  fit$config$anchor_review <- anchor_review
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
    facet_interactions = facet_interactions,
    min_obs_per_interaction = as.numeric(min_obs_per_interaction),
    interaction_policy = as.character(interaction_policy),
    anchors = anchors,
    group_anchors = group_anchors,
    noncenter_facet = as.character(noncenter_facet),
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
    population_formula = population_formula,
    person_id = if (is.null(person_id)) NULL else as.character(person_id),
    population_policy = as.character(population_policy),
    facet_shrinkage = as.character(facet_shrinkage),
    facet_prior_sd = facet_prior_sd,
    shrink_person = isTRUE(shrink_person),
    attach_diagnostics = isTRUE(attach_diagnostics),
    package_version = as.character(utils::packageVersion("mfrmr"))
  )

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
#' @param fit An `mfrm_fit` object returned by [fit_mfrm()].
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
#' @return A tibble with one row per interaction cell. Returns an empty tibble
#'   when the fit has no model-estimated facet interactions.
#' @seealso [fit_mfrm()], [estimate_bias()], [compare_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(
#'   toy, person = "Person", facets = c("Rater", "Criterion"),
#'   score = "Score", method = "MML", model = "RSM",
#'   facet_interactions = "Rater:Criterion",
#'   quad_points = 7, maxit = 30
#' )
#' head(interaction_effect_table(fit))
#' @export
interaction_effect_table <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object returned by fit_mfrm().",
         call. = FALSE)
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
      "PTMEA", "PtMeaCorr"),
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
  pop$converged <- mfrm_inference_ready(fit)
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
      notes = "No population model was requested; this fit uses an unconditional normal person distribution."
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
  if (is.finite(mm_rcond) && mm_rcond < 1e-8) {
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
    ". Invalid rows were removed; duplicate keys keep the last row."
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

  list(
    model = as.character(cfg$model %||% NA_character_),
    method = as.character(cfg$method %||% NA_character_),
    person = as.character(cfg$person_col %||% NA_character_),
    facets = sort(as.character(cfg$facet_cols %||% fit$prep$facet_names %||% character(0))),
    score = as.character(cfg$score_col %||% NA_character_),
    weight = as.character(cfg$weight_col %||% NA_character_),
    step_facet = as.character(cfg$step_facet %||% NA_character_),
    slope_facet = as.character(cfg$slope_facet %||% NA_character_),
    facet_interactions = canonical_compare_interaction_terms(
      cfg$facet_interactions
    ),
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

  list(
    eligible = FALSE,
    reason = paste0(
      "Automatic nesting review currently supports only RSM nested inside PCM ",
      "or same-family fixed facet-interaction extensions under shared design ",
      "and constraints."
    ),
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
#'
#' @details
#' This function provides a compact descriptive bundle similar to the
#' pre-fit summaries commonly checked in TAM workflows:
#' sample size, score distribution, per-facet coverage, and linkage counts.
#' `psych::describe()` is used for numeric descriptives of score and weight.
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
#' - *Extreme scores*: persons or facet levels with all-minimum or
#'   all-maximum scores yield infinite logit estimates under JML;
#'   they are handled via Bayesian shrinkage under MML.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `overview`: confirms sample size, facet count, and category span.
#' - `missing_by_column`: identifies immediate data-quality risks.
#'   Understand why values are missing and whether the fitted missing-data
#'   handling matches the study design.
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
#' @section Typical workflow:
#' 1. Run `describe_mfrm_data()` on long-format input.
#' 2. Review `summary(ds)` and `plot(ds, ...)`.
#' 3. Resolve missingness/sparsity issues before [fit_mfrm()].
#'
#' @return A list of class `mfrm_data_description` with:
#' - `overview`: one-row run-level summary
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
#'
#' @seealso [fit_mfrm()], [review_mfrm_anchors()]
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
                               missing_codes = NULL,
                               include_person_facet = FALSE,
                               include_agreement = TRUE,
                               rater_facet = NULL,
                               context_facets = NULL,
                               agreement_top_n = NULL,
                               expected_design = NULL,
                               min_linking_persons = 2L) {
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
    missing_codes = missing_codes
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
    RatingMaxSource = prep$rating_max_source %||% "unknown"
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
#' @seealso [describe_mfrm_data()], [fit_mfrm()].
#'
#' @examples
#' dat <- data.frame(
#'   Person = paste0("P", 1:5),
#'   Rater = c("R1", "R1", "R2", "R2", "R2"),
#'   Score = c(1, 99, 2, -1, 3)
#' )
#' cleaned <- recode_missing_codes(dat, columns = "Score")
#' cleaned$Score
#' attr(cleaned, "mfrm_missing_recoding")
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
      severity = "warning",
      condition = "duplicate_person_facet_cells",
      facets = object$structural_missingness$settings$key %||% character(0),
      message = paste0(
        duplicates$DuplicateCells[1L],
        " Person x facet cell(s) contain multiple observed rows."
      ),
      action = "Aggregate, deduplicate, or add a distinguishing facet before fitting."
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
#' - `agreement`: observed-score agreement for the selected scorer facet (when
#'   available).
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

#' Review and normalize anchor/group-anchor tables
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
#' @param missing_codes Optional. `NULL` (default) is a no-op;
#'   `TRUE` or `"default"` converts the FACETS / SPSS / SAS sentinel
#'   set to `NA` on the score column before review while preserving person and
#'   facet identifiers. Supply a character vector to apply a custom code set
#'   across the person, facet, and score columns.
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
#' as `fit_mfrm()`, but returns a review object so constraints can be
#' checked *before* estimation.  Running the review first helps avoid
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
#' 2. Run `review_mfrm_anchors(...)`.
#' 3. Resolve issues, then fit with [fit_mfrm()].
#'
#' @return A list of class `mfrm_anchor_review` with:
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
#' review <- review_mfrm_anchors(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   anchors = anchors
#' )
#' review$issue_counts
#' summary(review)
#' p_review <- plot(review, draw = FALSE)
#' p_review$data$plot
#' @name review_mfrm_anchors
#' @export
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
                                dummy_facets = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    missing_codes = missing_codes
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
    dummy_facets = dummy_facets
  )
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

#' Summarize an anchor-review object
#'
#' @param object Output from [review_mfrm_anchors()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in issue previews.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary provides a compact pre-estimation review of anchor and
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
#' 1. Run [review_mfrm_anchors()] with intended anchors/group anchors.
#' 2. Review `summary(review)` and recommendations.
#' 3. Revise anchor tables, then call [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_anchor_review`.
#' @seealso [review_mfrm_anchors()], [fit_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' review <- review_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(review)
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

#' Plot an anchor-review object
#'
#' @param x Output from [review_mfrm_anchors()].
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
#' Base-R visualization helper for anchor-review outputs.
#'
#' @section Interpreting output:
#' - `"issue_counts"`: volume of each issue class.
#' - `"facet_constraints"`: anchored/grouped/free mix by facet.
#' - `"level_observations"`: observation support across levels.
#'
#' @section Typical workflow:
#' 1. Run [review_mfrm_anchors()].
#' 2. Start with `plot(review, type = "issue_counts")`.
#' 3. Inspect constraint and support plots before fitting.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [review_mfrm_anchors()], [make_anchor_table()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' review <- review_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(review, draw = FALSE)
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
#' instability.  Re-run [review_mfrm_anchors()] on the receiving data
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
#' @seealso [fit_mfrm()], [review_mfrm_anchors()]
#' @examples
#' toy <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", quad_points = 7, maxit = 30
#' )
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
#' @param whexact Logical controlling the ZSTD standardisation of
#'   mean-square fit statistics. `FALSE` (default) applies the
#'   Wilson-Hilferty cube-root transformation
#'   \eqn{(\mathrm{MnSq}^{1/3} - (1 - 2/(9\,\mathit{df})))/\sqrt{2/(9\,\mathit{df})}}
#'   (recommended; the Winsteps/FACETS convention for `WHEXACT=Y`).
#'   `TRUE` uses the simpler linear-normal standardisation
#'   \eqn{(\mathrm{MnSq} - 1)\sqrt{\mathit{df}/2}}, which is kept for
#'   backward compatibility with earlier mfrmr summaries and with
#'   FACETS' `WHEXACT=N` mode.
#' @param fit_df_method Degrees-of-freedom convention used for fit ZSTD.
#'   `"engine"` (default) keeps the package-native convention
#'   `DF_Infit = sum(Var * Weight)` and `DF_Outfit = sum(Weight)`.
#'   `"facets"` uses the FACETS/Wright-Masters fourth-moment approximation
#'   `df = 2 / q^2` as the primary `InfitZSTD` / `OutfitZSTD` basis and caps
#'   reported ZSTD values at +/-9. `"both"` keeps the engine convention as the
#'   primary columns and adds `*_FACETS` companion columns for comparison.
#' @param diagnostic_mode Diagnostic basis to compute: `"both"` (the
#'   current default) computes both the residual/EAP-based stack and
#'   the strict latent-integrated first-order marginal-fit companion;
#'   `"legacy"` keeps the residual/EAP-based stack only;
#'   `"marginal_fit"` returns only the marginal-fit companion. The
#'   `"both"` path adds a posterior-integrated pass that typically
#'   doubles to quintuples wall-clock time relative to `"legacy"`;
#'   pass `"legacy"` explicitly when iterating on large designs and
#'   only the residual stack is needed. Use `"both"` for RSM/PCM
#'   reporting fits because it enables [plot_marginal_fit()] and
#'   [plot_marginal_pairwise()] follow-up.
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
#' - **PTMEA**: within-element point-measure correlation between observed
#'   scores and fitted person measures. A positive value is directionally
#'   consistent with the fitted orientation; it is not a confirmatory test.
#'
#' The MnSq values and the ZSTD values should be read separately. `mfrmr` keeps
#' the package-native engine df convention by default because it is the basis
#' used by the fitted observation-level diagnostics. FACETS reports closely related
#' MnSq values but standardizes them with a Wright-Masters fourth-moment df
#' approximation (`df = 2 / q^2`) and caps reported ZSTD values. Use
#' `fit_df_method = "both"` to review these two standardization conventions
#' side by side without changing the primary `InfitZSTD` / `OutfitZSTD`
#' columns.
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
#'
#' **Heuristic misfit-screening guidelines (Bond & Fox, 2015):**
#' - MnSq < 0.5: overfit (too predictable; may inflate reliability)
#' - MnSq 0.5--1.5: productive for measurement
#' - MnSq > 1.5: underfit (noise degrades measurement)
#' - \eqn{|\mathrm{ZSTD}| > 2}: conventional approximate-normal review flag;
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
#'   For bounded `GPCM`, `fair_average` is retained with an unavailable
#'   status because that compatibility calculation is outside the documented
#'   generalized-model contract.
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
#' treated as exploratory. Separation, strata, and reliability follow the
#' Wright & Masters (1982) conventions:
#' \eqn{G = \mathrm{TrueSD}/\mathrm{RMSE}},
#' \eqn{R = G^2 / (1 + G^2)}, and \eqn{H = (4G + 1) / 3}.
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
#'   `CI_Method`)
#' - `overall_fit`: overall fit summary
#' - `fit`: element-level fit diagnostics
#' - `reliability`: facet-level model/real separation and reliability
#' - `precision_profile`: one-row summary of the active precision tier and its
#'   recommended use
#' - `precision_review`: package-native checks for SE, CI, and reliability
#' - `parameter_uncertainty`: MML observed-information uncertainty for
#'   structural parameters when available (`steps`, and bounded-`GPCM`
#'   `slopes` on both log and positive scales), plus covariance status metadata
#' - `facet_precision`: facet-level precision summary by distribution basis and
#'   SE mode
#' - `facets_chisq`: fixed/random facet variability summary
#' - `interactions`: top interaction diagnostics
#' - `interrater`: inter-rater agreement bundle (`summary`, `pairs`) including
#'   agreement and rater-severity spread indices
#' - `unexpected`: unexpected-response bundle
#' - `fair_average`: adjusted-score reference bundle (reported as unavailable
#'   for bounded `GPCM`)
#' - `displacement`: displacement diagnostics bundle
#' - `approximation_notes`: method notes for SE/CI/reliability summaries
#' - `diagnostic_basis`: guide to the statistical target of each diagnostic path
#' - `fit_standardization`: guide to the df convention behind fit ZSTD values
#' - `marginal_fit`: optional strict marginal-fit companion based on
#'   posterior-expected first-order category counts
#' - `residual_pca_overall`: optional overall PCA object
#' - `residual_pca_by_facet`: optional facet PCA objects
#'
#' @seealso [fit_mfrm()], [analyze_residual_pca()], [build_visual_summaries()],
#'   [mfrmr_visual_diagnostics], [mfrmr_reporting_and_apa]
#' @examples
#' \donttest{
#' # Minimal diagnostic example without residual PCA.
#' toy <- load_mfrmr_data("example_operational")
#' fit_quick <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", model = "RSM", quad_points = 7, maxit = 30
#' )
#' diag_quick <- diagnose_mfrm(fit_quick, diagnostic_mode = "both",
#'                             residual_pca = "none")
#' summary(diag_quick)$overview[, c("Observations", "Facets", "Categories")]
#'
#' fit <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", model = "RSM", quad_points = 7, maxit = 30
#' )
#' diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
#' s_diag <- summary(diag)
#' s_diag$overview[, c("Observations", "Facets", "Categories")]
#' s_diag$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]
#' s_diag$key_warnings
#' # Look for: lines starting with "MnSq misfit:" name the element +
#' #   Infit / Outfit values outside the configured heuristic review band.
#' #   Review those signals in context; an empty warning list is not an
#' #   automatic all-clear decision.
#' s_diag$facets_chisq
#' # Look for: `FixedProb` < 0.05 is evidence against the fixed-effect
#' #   "all elements equal" null under the reported chi-square approximation.
#' #   Interpret the magnitude and precision as well; a non-significant result
#' #   does not demonstrate homogeneous elements or negligible facet spread.
#' s_diag$interrater
#' # Look for: ExactAgreement >= ExpectedExactAgreement and
#' #   AgreementMinusExpected >= 0 indicate raters agree at least as
#' #   often as the model expects. Negative values warrant a closer
#' #   look at `diag$interrater$pairs`.
#' p_qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE)
#' p_qc$data$plot
#'
#' # Optional: include residual PCA in the diagnostic bundle
#' diag_pca <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag_pca, mode = "overall")
#' head(pca$overall_table)
#'
#' # Reporting route:
#' prec <- precision_review_report(fit, diagnostics = diag)
#' summary(prec)
#'
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
diagnose_mfrm <- function(fit,
                          interaction_pairs = NULL,
                          top_n_interactions = 20,
                          whexact = FALSE,
                          fit_df_method = c("engine", "facets", "both"),
                          diagnostic_mode = c("both", "legacy", "marginal_fit"),
                          residual_pca = c("none", "overall", "facet", "both"),
                          pca_max_factors = 10L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
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
  class(out) <- c("mfrm_diagnostics", class(out))
  out
}

#' Compare two or more fitted MFRM models
#'
#' Produce a side-by-side comparison of multiple [fit_mfrm()] results using
#' information criteria, log-likelihood, and parameter counts. When exactly
#' two models are supplied and the current conservative nesting review passes,
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
#'   structural nesting review and computes the LRT only for supported nesting patterns.
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
#' The comparison table records both row count (`nobs`) and the sample-size
#' basis used for the BIC penalty (`ICSampleSize`, `ICSampleSizeBasis`); for
#' weighted fits this is the sum of weights rather than the number of rows.
#'
#' **Nesting**: Two models are *nested* when one is a special case of the
#' other obtained by imposing equality constraints.  The most common
#' nesting in MFRM is RSM (shared thresholds) inside PCM
#' (item-specific thresholds).  Models that differ only in estimation
#' method (MML vs JML) on the same specification are not nested in the
#' usual sense---use information criteria rather than LRT for that
#' comparison.
#'
#' In the **current `mfrmr` model space**, the automatic nesting review is
#' intentionally conservative. It currently supports two fixed-effect
#' restrictions under shared data and shared constraints:
#' - `RSM` nested inside `PCM` when the `PCM` fit has an explicit
#'   `step_facet`;
#' - same-family additive-vs-interaction comparisons when the smaller fit's
#'   `facet_interactions` set is a subset of the larger fit's set.
#'
#' Cross-method comparisons, comparisons that change anchors/dummying/centering,
#' and same-family comparisons that do not add fixed interaction terms are not
#' automatically promoted to LRT claims.
#'
#' The **likelihood-ratio test (LRT)** is reported only when exactly two
#' models are supplied, `nested = TRUE`, the structural nesting review passes, and the
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
#' - Do not interpret the LRT unless `nested = TRUE` and the structural nesting review
#'   in `comparison_basis$nesting_review` passes.
#' - Same-family additive-vs-interaction fits are considered nested only when
#'   all other structural settings match and the smaller model's
#'   `facet_interactions` set is a subset of the larger model's set.
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
#'   the supplied models. Inspect `comparison_basis$nesting_review$relation`
#'   and `reason` before reading any LRT output.
#' - `lrt`: nested-model test summary, present only when the requested and
#'   reviewed conditions are met.
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
#'   Delta_AIC, AkaikeWeight, Delta_BIC, BICWeight, npar, nobs, WeightedN,
#'   ICSampleSize, ICSampleSizeBasis, Model, Method, Converged,
#'   InferenceReady, ConvergenceSeverity, ICComparable).
#' - `lrt`: data.frame with likelihood-ratio test result (only when two models
#'   are supplied and `nested = TRUE`). Contains `ChiSq`, `df`, `p_value`.
#' - `evidence_ratios`: data.frame of pairwise Akaike-weight ratios (Model1,
#'   Model2, EvidenceRatio). `NULL` when weights cannot be computed.
#' - `preferred`: named list with the preferred model label by each criterion.
#' - `comparison_basis`: list describing whether IC and LRT comparisons were
#'   considered comparable. Includes a conservative `nesting_review` plus
#'   `lrt_status` / `lrt_reason` so withheld LRTs are explicit rather than
#'   silently absent.
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#'
#' fit_rsm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "MML", model = "RSM", quad_points = 7, maxit = 30)
#' fit_pcm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "MML", model = "PCM",
#'                      step_facet = "Criterion", quad_points = 7, maxit = 30)
#' comp <- compare_mfrm(fit_rsm, fit_pcm, labels = c("RSM", "PCM"))
#' comp$table
#' comp$evidence_ratios
#' }
#'
#' @section References:
#' - Burnham, K. P., & Anderson, D. R. (2002). *Model selection and
#'   multimodel inference: A practical information-theoretic
#'   approach* (2nd ed.). Springer.
#' - Akaike, H. (1974). A new look at the statistical model
#'   identification. *IEEE Transactions on Automatic Control,
#'   19*(6), 716-723.
#' - Schwarz, G. (1978). Estimating the dimension of a model.
#'   *Annals of Statistics, 6*(2), 461-464.
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
    weighted_n <- if (!is.null(f$config$weight_col) &&
        !is.null(f$prep$data) && "Weight" %in% names(f$prep$data)) {
      sum(f$prep$data$Weight, na.rm = TRUE)
    } else if (is.finite(nobs)) {
      as.numeric(nobs)
    } else {
      NA_real_
    }
    summary_n <- if ("N" %in% names(s)) {
      suppressWarnings(as.numeric(s$N[1]))
    } else {
      NA_real_
    }
    ic_sample_size <- if (is.finite(summary_n)) summary_n else weighted_n
    ic_sample_size_basis <- if (!is.null(f$config$weight_col) &&
        !is.null(f$prep$data) && "Weight" %in% names(f$prep$data)) {
      "sum_weights"
    } else {
      "row_count"
    }
    npar <- if (!is.null(f$opt$par)) length(f$opt$par) else NA_integer_
    method <- if (!is.null(f$config$method)) toupper(f$config$method[1]) else NA_character_
    method <- ifelse(identical(method, "JMLE"), "JML", method)
    convergence <- mfrm_convergence_state(f)
    tibble(
      Label     = labels[i],
      Model     = if (!is.null(f$config$model)) toupper(f$config$model[1]) else NA_character_,
      Method    = method,
      nobs      = nobs,
      WeightedN = weighted_n,
      ICSampleSize = ic_sample_size,
      ICSampleSizeBasis = ic_sample_size_basis,
      npar      = npar,
      LogLik    = if ("LogLik" %in% names(s)) s$LogLik[1] else NA_real_,
      AIC       = if ("AIC" %in% names(s)) s$AIC[1] else NA_real_,
      BIC       = if ("BIC" %in% names(s)) s$BIC[1] else NA_real_,
      Converged = convergence$code_converged,
      InferenceReady = convergence$inference_ready,
      ConvergenceSeverity = convergence$severity
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
  conv_vals <- tbl$InferenceReady
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
        "At least one compared model requires optimizer or convergence review. Raw AIC/BIC values are shown, ",
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
  lrt_status <- if (isTRUE(nested)) "requested" else "not_requested"
  lrt_reason <- if (isTRUE(nested)) {
    "Likelihood-ratio test was requested but not yet evaluated."
  } else {
    "Likelihood-ratio test was not requested; set `nested = TRUE` only for supported nested comparisons."
  }
  nesting_review <- audit_compare_mfrm_nesting(fits, labels = labels)
  if (isTRUE(nested)) {
    if (length(fits) != 2L) {
      lrt_status <- "not_computed"
      lrt_reason <- "Likelihood-ratio testing requires exactly two fitted models."
      warning("`nested = TRUE` was requested, but LRT requires exactly two fitted models. LRT was not computed.",
              call. = FALSE)
    } else if (!same_method || !same_nobs || !same_data || !all_converged || !all_mml) {
      lrt_status <- "not_computed"
      lrt_reason <- paste(
        "Models do not share the same formal MML likelihood basis,",
        "observation set, and convergence status."
      )
      warning(
        "`nested = TRUE` was requested, but the models do not share the same ",
        "formal MML likelihood basis, observation set, and convergence status. LRT was not computed.",
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
          lrt_status <- "computed"
          lrt_reason <- "Likelihood-ratio test computed after the structural nesting review passed."
        } else {
          lrt_status <- "not_computed"
          lrt_reason <- if (!is.finite(chi_sq) || !is.finite(df_diff)) {
            "Likelihood-ratio statistic or parameter-count difference was non-finite."
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
    all_inference_ready = all_converged,
    all_converged = all_converged,
    ic_comparable = ic_comparable,
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
    cat("a comparable formal MML likelihood basis, observation set, and inference-ready status.\n")
  }

  if (!is.null(x$lrt)) {
    cat("\nLikelihood-ratio test:\n")
    cat(sprintf("  Chi-sq = %.3f, df = %d, p = %.4f\n",
                x$lrt$ChiSq[1], x$lrt$df[1], x$lrt$p_value[1]))
    cat(sprintf("  %s vs %s\n", x$lrt$Simple[1], x$lrt$Complex[1]))
  } else if (isTRUE(x$comparison_basis$nested_requested)) {
    cat("\nLikelihood-ratio test was not reported.\n")
    if (!is.null(x$comparison_basis$lrt_status) &&
        !is.null(x$comparison_basis$lrt_reason)) {
      cat("  LRT status:", x$comparison_basis$lrt_status, "\n")
      cat("  Reason:", x$comparison_basis$lrt_reason, "\n")
    }
    review <- x$comparison_basis$nesting_review %||% list()
    if (!is.null(review$reason) && nzchar(review$reason)) {
      cat("  Nesting review:", review$reason, "\n")
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
