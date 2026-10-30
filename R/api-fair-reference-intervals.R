#' Intervals for common-plan Person expected scores
#'
#' Fixed-Person parametric refit intervals or a local delta approximation for a
#' declared common rating plan. Both preserve the calibration-to-anchor-to-scoring
#' estimator and account for explicitly shared response events. A separate
#' conservative bounded-response method uses the scoring observations directly.
#'
#' @param object For `confint()`, output from [fair_average_table()] with
#'   `reference_plan`. For `summary()`, a saved `mfrm_fair_intervals` result.
#' @param parm Character Person identifiers. Default selects all Persons.
#' @param level Pointwise confidence level in (0, 1).
#' @param nsim For `method = "refit"`, required integer number of replications,
#'   at least two. Small
#'   values are useful for checking the procedure, not precise tail estimation.
#'   Omit for `method = "delta"` or `"bounded"`.
#' @param mc_level For `method = "refit"`, confidence level in (0, 1) for
#'   Monte Carlo enclosures of both interval endpoints together, per target.
#'   Default `.95`. This measures finite-replication error, separately from
#'   `level`, and does not qualify sampling coverage. Omit for other methods.
#' @param seed For `method = "refit"`, required nonnegative integer seed.
#'   The caller's RNG state is restored and the RNG kind is recorded. Omit
#'   for `method = "delta"` or `"bounded"`, which do not simulate or change the RNG state.
#' @param method `"refit"` (default) for empirical basic limits from whole-pipeline
#'   refits; `"delta"` for pointwise normal-reference limits using the local
#'   derivative of the same estimator. These two methods are experimental.
#'   `"bounded"` gives conservative pointwise Hoeffding bounds for the same
#'   true expected-score target from independent bounded scoring responses,
#'   retaining the full score range for unobserved reference-plan cells.
#' @param generator_method For `method = "refit"` only, `"ordinary"` (default)
#'   uses the supplied ordinary-JML generating estimates. Explicit `"jeffreys"`
#'   re-estimates that generating model with a canonical Jeffreys information
#'   penalty. It retains ordinary JML for the original point and every refit.
#'   This experimental option can provide a finite generator when Person MLEs
#'   are infinite; it does not qualify coverage or correct the reported points.
#'   For a score result from `fit_mfrm(jml_penalty = "jeffreys")`, the default
#'   is instead `"jeffreys"`: use that already fitted generating model, and
#'   retain Jeffreys estimation for original points and every refit. Explicit
#'   `"ordinary"` is refused for these results.
#' @param contrasts Optional numeric matrix with named rows and columns equal
#'   to `parm` (order may differ). Rows must sum to zero and define differences
#'   of common-plan expected scores. Without it, return individual scores.
#' @param anchor_source Native ordinary or explicit Jeffreys JML fit that supplied every direct
#'   anchor in `object`. Values must match the unrounded primary estimates.
#'   Its own anchors are not supported in this initial two-stage route.
#'   Omit calibration/generator inputs, `response_map`, `anchor_uncertainty`
#'   and `gradient_tol` for `method = "bounded"`: its bounds do not use fitted
#'   anchors or an estimated generating model.
#' @param generating_fit Fit on the union of calibration and scoring response
#'   events, supplying the fitted response model. Defaults to `anchor_source`
#'   for estimated anchors and the scoring fit otherwise. A separate union fit
#'   is needed when scoring includes events outside the calibration data.
#' @param response_map For estimated anchors, a required list with numeric
#'   vectors `calibration` and `scoring`. Each maps rows of that fit's
#'   `prep$data` to rows of `generating_fit$prep$data`. Indices identify actual
#'   response events, not just equal Person/facet labels. Shared indices mean
#'   the same observation; distinct occasions require distinct generating rows.
#'   Identifiers and observed scores must agree, and the maps must cover the
#'   union. Providing these maps asserts the event correspondence.
#' @param anchor_uncertainty `"estimated"` (default) requires `anchor_source`
#'   when direct anchors are present. `"fixed"` explicitly conditions on them.
#'   It does not propagate estimated-anchor uncertainty.
#' @param gradient_tol Positive tolerance for the maximum absolute estimating
#'   gradient (ordinary or penalized according to the fitted estimator), default
#'   `1e-4`, applied to sources and every refit. Optimizer
#'   convergence, identification and boundary checks are also required.
#' @param ... Must be empty.
#'
#' @details
#' For refit and delta intervals, each Person is fixed at the generating fit's
#' primary estimate. Responses are
#' locally independent under its RSM/PCM model, with the retained assignment and
#' missingness pattern. The procedure does not draw a new Person population or
#' resample rows. Calibration and scoring must share the same model, category
#' coding, facet signs, estimator and parameter meanings. Only unit-weight additive
#' RSM/PCM JML without shrinkage or group anchors is supported. Score recoding
#' must be the identity; transformed category labels are not accepted here.
#' Ordinary points and Jeffreys points cannot be mixed across stages.
#' With explicit Jeffreys fits, `method = "refit"` uses the same penalty for
#' original calibration, anchor extraction, scoring and every replication.
#' Numerical checks evaluate the penalized equations; ordinary likelihood
#' recession does not turn these penalized points into infinities. Rank
#' deficiency and numerical failures still stop or remain unresolved in the
#' planned denominator. No failed replicate is redrawn or replaced by an
#' ordinary estimate. `method = "delta"` uses the analytic penalized score
#' Jacobian, including the penalty's dependence on estimated anchors. It does
#' not use ordinary information or inverse penalized curvature as covariance.
#' `method = "bounded"`
#' retains its separate response-based assumptions and fitted-point display.
#' The penalty changes the points; it does not qualify score bias or coverage.
#'
#' With `method = "refit"`, the original two-stage estimator is preserved: first refit
#' calibration, extract the specified anchors, then refit the scoring data under
#' those anchors and compute scores on the saved common plan. The generator is
#' a plug-in response model, not a replacement joint estimator for the data.
#' The bootstrap error is replicate score minus the common-plan score under the
#' **generating** model; this need not equal the original two-stage estimate.
#' Basic limits subtract its empirical quantiles from the original estimate.
#' Contrasts use paired scores from each whole replication, retaining covariance.
#'
#' Unresolved refits remain in the planned denominator. Replacing their unknown
#' errors by minus/plus infinity gives outer empirical-quantile bounds; the
#' function does not silently condition on successful refits. Nonfinite bounds
#' are returned as missing with `Status = "unresolved_tail"`. Finite bounds
#' outside the score range are retained, not clipped to disguise approximation
#' error. Type-1 empirical quantiles and expected observations per tail are
#' recorded. These are pointwise intervals, without multiplicity adjustment.
#' Inspect `LowerMCLower`, `LowerMCUpper`, `UpperMCLower` and `UpperMCUpper`
#' for Monte Carlo enclosures of the lower and upper basic endpoints. Unlike
#' the reported empirical interval, these enclosures retain infinite limits
#' when replication counts or unresolved refits cannot bound an endpoint.
#' `MCFinite` means all four limits are finite, not sufficiently precise.
#' The `monte_carlo` attribute retains quantile probabilities and order ranks.
#' For each error quantile p, the ranks are `qbinom(a, nsim, p)` and
#' `qbinom(1-a, nsim, p)+1`, where `a = (1-mc_level)/4`. Ranks zero and
#' `nsim+1` denote infinities. The four one-sided allocations jointly cover
#' both infinite-replication basic endpoints for one target at least at
#' `mc_level`, conditional on a fixed generating model and independent draws.
#' They are conservative for discrete errors/ties. There is no simultaneous
#' guarantee across targets, reruns or adaptively chosen stopping times.
#' Unknown errors retain their minus/plus-infinity envelopes; the enclosures
#' include every completion of failed draws. If the estimator is not defined
#' on some outcomes, this does not create a sampling distribution for it.
#' Increasing `nsim` cannot repair structural nonidentification or qualify bias.
#' With finite qualified structural calibration, a replicate Person at minus
#' or plus infinity has a well-defined limiting expected score (the scale
#' endpoint). That limit is retained and counted in `trials$BoundaryPersons`,
#' without using a finite optimizer proxy. An infinite estimated anchor cannot
#' be passed into the next fitting stage and remains unresolved. The generating
#' fit and originally selected Persons must have finite primary estimates for
#' the default ordinary generator. The explicit Jeffreys option below also
#' accepts original selected limiting expected scores under finite structural
#' calibration; its own generating parameters must be finite and qualified.
#'
#' With `method = "delta"`, the derivative is evaluated at the observed
#' calibration and scoring fits. Each stage's information uses that stage's
#' fitted category probabilities; the variance and cross-covariance of response
#' events use the declared generating fit. Shared events contribute to both
#' stages before their categorical covariance is taken. Distinct occasions with
#' identical labels are independent under this fixed-Person model unless their
#' event indices are shared. The target covariance retains cross-Person and
#' cross-stage terms, which can increase or decrease a contrast's variance.
#' A direct anchor changes both its own location and any dependent centering
#' coordinate; the derivative includes that affine change. This is model-based
#' propagation, not an empirical robust/cluster sandwich estimate.
#'
#' Delta propagation requires finite primary Person estimates throughout the
#' generating, calibration and scoring fits and numerically positive-definite
#' score Jacobians. Selected-target solves use diagonal scaling without
#' regularization. For ordinary JML no full parameter inverse is formed.
#' The Jeffreys derivative additionally solves for centered category designs
#' and processes one derivative direction at a time, without event-squared or
#' parameter-cubed arrays. Its scaled Jacobian reciprocal condition must exceed
#' `sqrt(.Machine$double.eps)`; numerical availability can fail in weak designs.
#' Information, conditioning or boundary failures
#' preserve the point estimates, return missing SEs/limits and a `Reason`, and
#' set `Status = "unavailable_information"`. Normal limits are not clipped to the
#' score support; large widths or impossible endpoints can expose a poor local
#' approximation. No simulations, resolved-replication counts or tail counts are
#' implied by this method (those count columns are `NA`).
#' Even an unselected extreme Person in the generating fit can make uncertainty
#' unavailable for otherwise finite selected targets; their point estimates and
#' the reason are retained. Do not remove that Person merely to obtain intervals.
#' Delta limits are centered on the observed estimate and do not perform a
#' bootstrap bias correction or correct a discrepancy between stage models.
#'
#' Refit and delta are experimental fitted-model approximations. Neither removes JML
#' incidental-parameter bias, certify nominal coverage, include uncertainty in
#' the chosen plan, predict a future observed score, or establish fairness.
#' Model misspecification, dependence beyond the model and changing Persons or
#' calibrations across occasions require different procedures. `CIEligible`
#' remains false pending finite-sample qualification; `Available` describes
#' numerical availability only. No external-software coverage guarantee is implied.
#' Correct local covariance does not guarantee adequate normal limits. In the
#' supported binary special case with six independent responses, the Jeffreys
#' expected-score estimate is \eqn{(S+0.5)/7}. At true success probability 0.5,
#' the plug-in delta 95% interval covers only for \eqn{S=2,3,4}: exact coverage
#' is 0.78125, despite an unbiased point at that probability. At probability
#' 0.25, exact coverage is 0.784423828125. These are exhaustive binomial
#' calculations, not Monte Carlo estimates. They also hold when the estimated
#' Person anchor is passed to fully fixed scoring on the same responses.
#' The normal-reference approximation, including plug-in standard errors and
#' discrete outcomes, therefore has a demonstrated small-sample limitation.
#' More Persons do not necessarily increase the responses for each focal Person.
#' The basic refit route is not a general remedy: exact, infinite-replication
#' binomial refit quantiles for the same estimator give coverage 0.78125 at
#' probability 0.5 and 0.452709 at probability 0.1. For an all-zero sample,
#' its upper basic limit equals the fitted score \eqn{0.5/7}; the bootstrap
#' estimator cannot fall below its own score floor. This is an interval-method
#' limitation, not inner Monte Carlo error, and more refits do not remove it.
#'
#' @section Explicit Jeffreys generating model:
#' `generator_method = "jeffreys"` minimizes
#' \eqn{-\ell(b)-\tfrac12\log|I(b)|} in the generating fit's actual affine free
#' coordinates, where \eqn{I} is the ordinary categorical model information.
#' Direct anchor values enter as fixed offsets in this generating model.
#' The calibrated anchors used by the scoring estimator are still re-estimated
#' on every replication, with shared response events generated only once.
#' This option requires declared rubric endpoints and unchanged full score
#' coding; it is unavailable for delta and bounded-response intervals.
#'
#' Computation starts at the free-coordinate origin and uses BFGS, the generating
#' fit's `maxit`, and the stricter of its `reltol` and `1e-14`. The analytic
#' penalized gradient must pass `gradient_tol`. Information must be numerically
#' positive definite; scaled sparse solves must pass a `1e-8` residual check.
#' There is no ridge, finite parameter box or probability clipping. Optimization,
#' rank or numerically noninterior probability failures stop with a reason;
#' success is not a certificate of a global optimum. Ordinary source fits and
#' every refit retain their existing stationarity and structural checks.
#'
#' The original ordinary-JML point and its estimator are unchanged. Basic limits
#' use replicate ordinary-JML errors relative to the Jeffreys generating
#' expected score, not the original fitted generating score. Finite expected-score
#' limits for original extreme Persons are retained; unresolved structural
#' calibration or an infinite anchor passed between stages remains unsupported.
#' Inspect `generating_model` and `generating_target` attributes for the distinct
#' generating parameters, probabilities, target and numerical diagnostics.
#'
#' Firth's first-order bias reduction concerns canonical parameters in regular
#' fixed-dimensional asymptotics. It does not establish bias reduction for
#' nonlinear expected scores, fixed-short-exposure JML as Person count increases,
#' or this two-stage interval procedure. For example, with \eqn{S\sim Bin(n,p)},
#' the induced estimate \eqn{(S+1/2)/(n+1)} has bias \eqn{(1/2-p)/(n+1)}, whereas
#' \eqn{S/n} is unbiased for \eqn{p}. Thus a finite generator can worsen a score
#' bias. Sparse/weak-link coverage and the bias/variance tradeoff require separate
#' qualification; `CIEligible` remains false. This is distinct from FACETS Xtreme
#' and from the corrected-GPCM `jml_correction_order` estimator.
#'
#' @section Conservative bounded-response intervals:
#' `method = "bounded"` targets the same true common-plan expected score or
#' contrast, using the original scoring rows. It does not propagate a fitted
#' covariance or generate responses from extreme Person estimates. A finite-sample
#' probability inequality applies under independent responses, a fixed assignment,
#' plan and contrast chosen independently of the observed scores, known score
#' bounds, and stable expected scores within each Person-by-plan cell. Repeated
#' rows must represent independent observations, not copied ratings. Dependence,
#' adaptive allocation and response-dependent missingness require other methods.
#' Both `rating_min` and `rating_max` must have been explicitly declared from
#' the rubric when fitting, with an unchanged full `score_map`. Use
#' `category_policy = "preserve"` to retain that scoring rule. Limits inferred
#' from observed minima/maxima, missing range provenance and internally collapsed
#' categories are rejected: a small sample can omit a valid category, and its
#' observed range is not known population support. Declaring the sample extrema
#' after seeing the data does not satisfy this assumption. If a retained empty
#' category prevents fitting, that failure remains in an end-to-end denominator;
#' do not collapse the category merely to obtain an interval.
#'
#' Each observed cell contributes its sample mean times its normalized plan
#' weight; cells outside the plan receive zero weight. An unobserved cell retains
#' the full declared score range. For a contrast, signed event weights preserve
#' its definition. Bounds are intersected with the target's mathematically known
#' range. This preserves the coverage inequality and is not clipping a normal
#' approximation. Sparse plans can produce very wide or uninformative bounds;
#' the full range has `Informative = FALSE`, even though it is a valid bound.
#' `Informative = TRUE` means only narrower than the full target range, not
#' necessarily precise enough for a substantive decision.
#'
#' `Estimate` retains the original JML score, including an endpoint limit or
#' missing estimate. The interval need not be centered on or contain that point.
#' `DesignEstimate` is a separate unbiased weighted-mean estimate when every
#' relevant Person has all reference-plan cells observed; otherwise it is missing.
#' It is not a corrected JML estimate. The confidence statement concerns the true
#' expected score, not the accuracy or bias of the retained JML point.
#'
#' Estimated anchors, even from the same sample, affect the retained JML point
#' but are not inputs to these bounds. Their uncertainty is therefore neither
#' conditioned away nor approximated by a fixed-anchor SE. The price is foregoing
#' model-based extrapolation to unobserved plan cells. `CIEligible = TRUE` here
#' denotes the finite-sample inequality under the stated assumptions, not empirical
#' qualification of JML or guaranteed usefulness. Bounds are pointwise; inspecting
#' many Persons or selecting a contrast after seeing results is not covered.
#' The inequality is unconditional over the prespecified response experiment.
#' It does not imply nominal coverage conditional on successful model fitting,
#' an extreme-score subgroup or another response-dependent selection. This API
#' requires an existing common-plan result: if fitting cannot create that result,
#' retain the failure in the planned denominator of an end-to-end evaluation.
#'
#' If \eqn{a_e} is a signed event weight and scores lie in \eqn{[l,u]}, write
#' \eqn{T=\sum_e a_eY_e} for the observed contribution. With
#' \eqn{\alpha=1-\mathrm{level}}, the sampling radius is
#' \deqn{r=\sqrt{\tfrac12\log(2/\alpha)\sum_e a_e^2(u-l)^2}.}
#' Each missing cell with signed plan weight \eqn{b} contributes between
#' \eqn{\min(bl,bu)} and \eqn{\max(bl,bu)}. Add these lower/upper sums to
#' \eqn{T-r,T+r}. Hoeffding's two tails bound the probability of missing the
#' true expected-score target by at most \eqn{\alpha}, uniformly over those
#' missing means. Extreme observed totals do not invalidate this argument.
#'
#' @section Local two-stage derivative:
#' Let \eqn{b_c, b_s} be the actual free calibration/scoring coordinates,
#' \eqn{a(b_c)} the extracted direct anchors, and \eqn{t(b_s,a)} the selected
#' expected scores or contrasts. Category-logit coefficients in the scoring
#' stage are \eqn{X_s} for free coordinates and \eqn{O_s} for anchors, including
#' centering changes. Write \eqn{H_c,H_s} for stage information,
#' \eqn{B = \sum_e \operatorname{Cov}_s(X_s,O_s)},
#' \eqn{A = \partial a/\partial b_c}, and \eqn{G_s,G_a} for target derivatives.
#' For explicit Jeffreys estimation replace each \eqn{H} by the Hessian of
#' \eqn{-\ell-\tfrac12\log|I|} in its stage's free coordinates, and replace
#' \eqn{B} by the mixed free-coordinate/anchor derivative of that same
#' objective. This includes the anchor dependence of \eqn{I}; holding the
#' penalty fixed during anchor propagation would give a different derivative.
#' The penalty does not depend on response outcomes at fixed design/parameters,
#' so the response-event coefficients \eqn{X_s,X_c} remain unchanged.
#' Differentiating the two score equations gives event-category coefficients
#' \deqn{U_s = X_s H_s^{-1}G_s^T,\qquad
#' U_c = X_c H_c^{-1} A^T (G_a-G_s H_s^{-1}B)^T.}
#' Put both coefficients on their explicit union-event indices (zero when an
#' event is absent from a stage), then compute
#' \eqn{V = \sum_e \operatorname{Cov}_{p_{g,e}}(U_s+U_c)} over categories.
#' Centering under the generating probabilities matters when stage and generator
#' fits differ. The saved decomposition is direct scoring, propagation through
#' calibration, and the signed shared-event cross term. Their sum is the target
#' covariance; the cross term alone need not be positive semidefinite.
#' This local derivation does not establish fixed-item JML asymptotics or
#' finite-sample normality. Estimated plans and unmodelled dependence are excluded.
#'
#' @return An `mfrm_fair_intervals` data frame with estimates, numerical bounds,
#'   availability and interpretation columns. Refit attributes retain `draws`,
#'   `person_draws`, `anchor_draws`, `trials`, `settings`, `reference_plan`,
#'   `generating_target`, `response_map`, `source_checks` and source identities.
#'   Refit results also retain `MCLevel`, four endpoint-enclosure columns,
#'   `MCFinite` and the `monte_carlo` rank table, without changing `Available`
#'   or promoting `CIEligible` on the strength of numerical precision.
#'   The explicit Jeffreys route additionally retains `generating_model`, with
#'   free parameters, probabilities, numerical checks and optimization controls.
#'   Delta results additionally have `SE` and `Reason` columns and attributes
#'   `covariance` (ordered by Target), `covariance_components` and `linearization`
#'   (target/anchor derivatives, selected solve loadings, parameter coordinates
#'   and solve residuals). The source/plan/settings attributes are shared by both
#'   methods. Unavailable delta covariance is `NA`, not zero.
#'   Bounded results instead retain `Informative`, `DesignEstimate`,
#'   `ObservedContribution`, `SamplingRadius`, unobserved-cell and target-support
#'   limits, a per-Person `support` attribute and sparse `person_event_weights`.
#'   Their `settings$score_support` retains the declared limits and coding basis.
#'   Their SE and simulation-count columns are missing; the radius is not an SE.
#'   `summary()` reports the same table without recalculation; `as.data.frame()`
#'   supports export and `saveRDS()` preserves the full evidence. See
#'   [plot.mfrm_fair_intervals()] for display and interpretation.
#' @references Murphy, K. M., and Topel, R. H. (1985). Estimation and inference
#'   in two-step econometric models. Journal of Business & Economic Statistics,
#'   3(4), 370--379. \doi{10.1080/07350015.1985.10509471}.
#'   Related estimating-equation framework: Hardin, J. W. (2002). The robust
#'   variance estimator for two-stage models. The Stata Journal, 2(3), 253--266.
#'   \doi{10.1177/1536867X0200200302}. The present method uses a model-based event
#'   covariance and does not claim the robust estimator described there.
#'   Hoeffding, W. (1963). Probability inequalities for sums of bounded random
#'   variables. Journal of the American Statistical Association, 58(301), 13--30,
#'   Theorem 2. \doi{10.1080/01621459.1963.10500830}.
#'   Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#'   Biometrika, 80(1), 27--38. \doi{10.1093/biomet/80.1.27}.
#'   Burke, S. (2016, revised 2018). Confidence Intervals for the Median and
#'   Other Percentiles. STAT COE-Report-08-2016. Binomial order-statistic
#'   construction; the present Monte Carlo diagnostic additionally allocates
#'   errors across both endpoints and retains unresolved-draw envelopes.
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("confint.mfrm_fair_average", run.dontrun = TRUE)
#' \dontrun{
#' d <- load_mfrmr_data("example_core")
#' # Calibrate the observed ratings with ordinary joint maximum likelihood.
#' f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
#'   method = "JML", rating_min = 1, rating_max = 4, category_policy = "preserve",
#'   optimizer = "BFGS", maxit = 500, reltol = 1e-11)
#' # The bundled example's rubric is 1--4; declare it independently of sample extrema.
#' # Every Person is evaluated against these same Rater-Criterion cells.
#' plan <- unique(d[c("Rater", "Criterion")])
#' fair <- fair_average_table(f, reference_plan = plan)
#' head(fair$reference_scores) # Expected average, in internal score units.
#' # Two replications exercise the API only; they cannot qualify coverage.
#' # Explicit experimental generator; points and every refit remain ordinary JML.
#' ci <- confint(fair, parm = c("P001", "P002"), nsim = 2, seed = 71,
#'   generator_method = "jeffreys",
#'   contrasts = matrix(c(1, -1), 1,
#'     dimnames = list("P001 minus P002", c("P001", "P002"))))
#' summary(ci) # Read Status and Available alongside the limits.
#' attr(ci, "trials") # Includes unresolved refits and their reasons.
#' attr(ci, "monte_carlo") # Infinite enclosures expose the two-draw tail uncertainty.
#' plot(ci, draw = FALSE) # Saved plotting data, including the cautions.
#' # Deterministic local approximation to the same expected-score target:
#' delta <- confint(fair, parm = c("P001", "P002"), method = "delta")
#' delta[, c("Target", "Estimate", "SE", "Available", "CIEligible")]
#' attr(delta, "covariance") # Includes the covariance between Persons.
#' # Conservative alternative with explicit observed/unobserved plan support:
#' bounded <- confint(fair, parm = c("P001", "P002"), method = "bounded")
#' bounded[, c("Target", "Estimate", "DesignEstimate", "CI_Lower", "CI_Upper", "Informative")]
#' attr(bounded, "support") # Check how much of the plan was actually observed.
#' }
#' @export
#' @order 1
confint.mfrm_fair_average <- function(object, parm = NULL, level = .95,
    nsim, seed, contrasts = NULL, anchor_source = NULL, generating_fit = NULL,
    response_map = NULL, anchor_uncertainty = c("estimated", "fixed"),
    gradient_tol = 1e-4, method = c("refit", "delta", "bounded"),
    generator_method = c("ordinary", "jeffreys"), mc_level = .95, ...) {
  rlang::check_dots_empty()
  anchor_uncertainty_supplied <- !missing(anchor_uncertainty)
  anchor_uncertainty <- match.arg(anchor_uncertainty)
  method <- match.arg(method)
  generator_method_supplied <- !missing(generator_method)
  if (method != "refit" && !missing(generator_method))
    stop("`generator_method` applies only to method = 'refit'.", call. = FALSE)
  generator_method <- match.arg(generator_method)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(level) || level <= 0 || level >= 1) stop("`level` must be in (0, 1).", call. = FALSE)
  if (method != "refit" && !missing(mc_level))
    stop("`mc_level` applies only to method = 'refit'.", call. = FALSE)
  if (!scalar(mc_level) || mc_level <= 0 || mc_level >= 1)
    stop("`mc_level` must be in (0, 1).", call. = FALSE)
  if (method == "refit" && (missing(nsim) || !scalar(nsim) || nsim < 2 || nsim != floor(nsim) || nsim > .Machine$integer.max))
    stop("Supply an integer `nsim` of at least two.", call. = FALSE)
  if (method == "refit" && (missing(seed) || !scalar(seed) || seed < 0 || seed != floor(seed) || seed > .Machine$integer.max))
    stop("Supply a nonnegative integer `seed`.", call. = FALSE)
  if (method != "refit" && (!missing(nsim) || !missing(seed)))
    stop("`nsim` and `seed` apply only to method = 'refit'; omit them for deterministic intervals.", call. = FALSE)
  if (!scalar(gradient_tol) || gradient_tol <= 0) stop("`gradient_tol` must be positive.", call. = FALSE)
  fit <- object$reference_fit; plan <- object$reference_plan
  if (is.null(fit) || is.null(plan))
    stop("Recreate fair_average_table() with an explicit reference_plan.", call. = FALSE)
  mfrm_fair_reference_fit(fit)
  penalized <- mfrm_is_jeffreys_fit(fit)
  if (penalized && method == "refit") {
    if (generator_method_supplied && generator_method != "jeffreys")
      stop("Jeffreys point/refit estimation requires generator_method = 'jeffreys'.", call. = FALSE)
    generator_method <- "jeffreys"
  }
  if (!identical(object$reference_identity, mfrm_fair_fit_identity(fit)))
    stop("The retained reference fit changed; recreate fair_average_table().", call. = FALSE)
  plan <- mfrm_fair_reference_plan(fit, plan[setdiff(names(plan), "NormalizedWeight")])
  scores <- mfrm_fair_reference_scores(fit, plan)
  if (!isTRUE(all.equal(scores, object$reference_scores, tolerance = 1e-12)))
    stop("The retained reference plan or scores changed; recreate fair_average_table().", call. = FALSE)
  if (is.null(parm)) parm <- scores$Person
  if (!is.character(parm) || !length(parm) || anyNA(parm) || anyDuplicated(parm) ||
      any(!parm %in% scores$Person)) stop("`parm` must select distinct known Person identifiers.", call. = FALSE)
  if (is.null(contrasts)) {
    L <- diag(length(parm)); dimnames(L) <- list(parm, parm)
  } else {
    if (!is.matrix(contrasts) || !is.numeric(contrasts) || !nrow(contrasts) ||
        any(!is.finite(contrasts)) || is.null(rownames(contrasts)) ||
        anyNA(rownames(contrasts)) || any(!nzchar(rownames(contrasts))) ||
        anyDuplicated(rownames(contrasts)) || anyDuplicated(colnames(contrasts)) ||
        !setequal(colnames(contrasts), parm) ||
        any(abs(rowSums(contrasts)) > 1e-12 * pmax(1, rowSums(abs(contrasts)))) ||
        any(rowSums(abs(contrasts)) == 0))
      stop("`contrasts` needs finite nonzero, zero-sum rows and unique row/Person column names.", call. = FALSE)
    L <- contrasts[, parm, drop = FALSE]
  }
  selected <- scores$Estimate[match(parm, scores$Person)]
  estimate <- vapply(seq_len(nrow(L)), function(i) {
    use <- L[i, ] != 0
    sum(L[i, use] * selected[use])
  }, 0)
  if (method == "bounded") {
    if (!is.null(anchor_source) || !is.null(generating_fit) || !is.null(response_map) ||
        anchor_uncertainty_supplied || !missing(gradient_tol))
      stop("Bounded-response intervals use scoring responses directly; omit calibration/generator inputs, anchor_uncertainty and gradient_tol.", call. = FALSE)
    return(mfrm_fair_bounded_result(fit, plan, parm, L, estimate, level,
      has_contrasts = !is.null(contrasts)))
  }
  anchors <- fit$config$replay_inputs$anchors
  two_stage <- nrow(anchors) > 0L && anchor_uncertainty == "estimated"
  if (two_stage && is.null(anchor_source))
    stop("Estimated-anchor intervals are unavailable without the calibration fit and response-event maps. Use anchor_uncertainty = 'fixed' only for explicitly conditional inference.", call. = FALSE)
  if (!two_stage && (!is.null(anchor_source) || !is.null(response_map) || !is.null(generating_fit)))
    stop("Separate calibration/generator inputs require estimated anchors.", call. = FALSE)
  generator <- generating_fit %||% (if (two_stage) anchor_source else fit)
  sources <- if (two_stage) list(scoring = fit, calibration = anchor_source, generating = generator) else list(scoring = fit)
  source_checks <- lapply(sources, function(f) {
    mfrm_fair_reference_fit(f)
    if (!identical(mfrm_is_jeffreys_fit(f), penalized))
      stop("Calibration, scoring and generating sources must use the same estimator; ordinary and Jeffreys points cannot be mixed.", call. = FALSE)
    replay <- f$config$replay_inputs
    if (is.null(replay) || any(f$prep$data$Weight != 1) || nrow(replay$group_anchors %||% data.frame()) ||
        !identical(as.numeric(f$prep$score_map$OriginalScore), as.numeric(f$prep$score_map$InternalScore)))
      stop("Common-plan intervals require replay inputs, unit weights, identity score coding and no group anchors.", call. = FALSE)
    keys <- c("model", "n_cat", "rating_min", "rating_max", "facet_names", "facet_signs", "step_facet", "dummy_facets")
    if (!isTRUE(all.equal(f$config[keys], fit$config[keys], check.attributes = FALSE)))
      stop("Calibration, scoring and generating models/coding must agree.", call. = FALSE)
    mfrm_fair_refit_check(f, gradient_tol)
  })
  if (two_stage) {
    if (nrow(anchor_source$config$replay_inputs$anchors))
      stop("The initial two-stage route requires an unanchored calibration source.", call. = FALSE)
    values <- mfrm_fair_anchor_values(anchor_source, anchors)
    if (any(abs(values - anchors$Anchor) > 1e-10 * pmax(1, abs(values))))
      stop("Scoring anchors do not match the unrounded calibration estimates.", call. = FALSE)
    if (!is.list(response_map) || !setequal(names(response_map), c("calibration", "scoring")))
      stop("Provide explicit calibration/scoring response-event row maps.", call. = FALSE)
    maps <- list(calibration = mfrm_fair_response_map(anchor_source, generator,
      response_map$calibration, "response_map$calibration"),
      scoring = mfrm_fair_response_map(fit, generator, response_map$scoring, "response_map$scoring"))
    if (!setequal(unlist(maps), seq_len(nrow(generator$prep$data))))
      stop("The response maps must cover the generating fit's union of events.", call. = FALSE)
  } else maps <- list(scoring = seq_len(nrow(fit$prep$data)))
  generator_theta <- generator$facets$person$PrimaryEstimate %||% generator$facets$person$Estimate
  # A nonselected extreme Person makes this local approximation unavailable,
  # but must not discard finite selected points and their explicit diagnostic.
  generating_model <- if (penalized) c(generator$checks,
    list(par = generator$opt$par, probabilities = mfrm_fair_generator_probs(generator),
      parameter_map = generator$parameter_map)) else if (generator_method == "jeffreys")
    mfrm_fair_jeffreys_generator(generator, gradient_tol) else NULL
  probs <- if (!is.null(generating_model)) generating_model$probabilities else
    if (method == "delta" && any(!is.finite(generator_theta))) NULL else
      mfrm_fair_generator_probs(generator)
  if (!is.null(probs) && (any(!is.finite(probs)) || any(probs < 0) || any(abs(rowSums(probs) - 1) > 1e-10)))
    stop("Invalid generating probabilities.", call. = FALSE)
  truth_scores <- if (!is.null(generating_model))
    mfrm_fair_jeffreys_target(generator, generating_model$par, plan, parm) else {
      tab <- mfrm_fair_reference_scores(generator, plan)
      tab$Estimate[match(parm, tab$Person)]
    }
  truth <- drop(L %*% truth_scores)
  allowed_status <- if (generator_method == "jeffreys") c("point_estimate", "boundary_limit") else "point_estimate"
  if (any(!is.finite(truth)) || any(!is.finite(estimate)) ||
      any(!scores$Status[match(parm, scores$Person)] %in% allowed_status))
    stop("Selected targets require finite primary point estimates and generating targets.", call. = FALSE)
  if (method == "delta") return(mfrm_fair_delta_result(fit,
    if (two_stage) anchor_source else NULL, generator, plan, parm, L, estimate, truth,
    anchors, maps, probs, level, gradient_tol, anchor_uncertainty, sources, source_checks,
    has_contrasts = !is.null(contrasts)))
  draws <- matrix(NA_real_, nsim, nrow(L), dimnames = list(NULL, rownames(L)))
  person_draws <- matrix(NA_real_, nsim, length(parm), dimnames = list(NULL, parm))
  anchor_draws <- matrix(NA_real_, nsim, nrow(anchors),
    dimnames = list(NULL, if (nrow(anchors)) paste(anchors$Facet, anchors$Level, sep = ":") else NULL))
  trials <- data.frame(Replication = seq_len(nsim), CalibrationGradient = NA_real_,
    ScoringGradient = NA_real_, BoundaryPersons = NA_integer_, Available = FALSE,
    Reason = "", Warnings = "", Seconds = NA_real_)
  rng_kind <- RNGkind()
  with_preserved_rng_seed(seed, {
    cumulative <- t(apply(probs, 1L, cumsum))
    for (b in seq_len(nsim)) {
      started <- proc.time()[[3]]; messages <- character()
      result <- tryCatch(withCallingHandlers({
        # One variate per union EVENT, shared by both mapped stages.
        u <- stats::runif(nrow(probs))
        y <- pmin(rowSums(u > cumulative), ncol(probs) - 1L) + generator$prep$rating_min
        a <- anchors
        if (two_stage) {
          calibration <- mfrm_fair_refit(anchor_source, y[maps$calibration])
          check <- mfrm_fair_refit_check(calibration, gradient_tol)
          trials$CalibrationGradient[b] <- check$max_gradient
          a$Anchor <- mfrm_fair_anchor_values(calibration, anchors)
        }
        if (nrow(a)) anchor_draws[b, ] <- a$Anchor
        scored <- mfrm_fair_refit(fit, y[maps$scoring], a)
        check <- mfrm_fair_refit_check(scored, gradient_tol)
        trials$ScoringGradient[b] <- check$max_gradient
        s <- mfrm_fair_reference_scores(scored, plan)
        v <- s$Estimate[match(parm, s$Person)]
        status <- s$Status[match(parm, s$Person)]
        trials$BoundaryPersons[b] <- sum(status == "boundary_limit")
        v[!status %in% c("point_estimate", "boundary_limit")] <- NA_real_
        person_draws[b, ] <- v
        # Zero coefficients must not propagate irrelevant missing Persons.
        draws[b, ] <- vapply(seq_len(nrow(L)), function(j) {
          use <- L[j, ] != 0
          sum(L[j, use] * v[use])
        }, 0)
        if (all(is.finite(draws[b, ]))) "" else "Selected Person expected score unavailable."
      }, warning = function(w) {
        messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) conditionMessage(e))
      trials$Available[b] <- all(is.finite(draws[b, ]))
      trials$Reason[b] <- result
      trials$Warnings[b] <- paste(unique(messages), collapse = " | ")
      trials$Seconds[b] <- proc.time()[[3]] - started
    }
  })
  bounds <- mfrm_fair_bootstrap_bounds(estimate, truth, draws, level)
  mc <- mfrm_fair_bootstrap_mc(estimate, bounds$roots, level, mc_level)
  mc_lower <- mc[mc$Endpoint == "lower", ]; mc_upper <- mc[mc$Endpoint == "upper", ]
  available <- is.finite(bounds$lower) & is.finite(bounds$upper)
  out <- data.frame(Target = rownames(L), Estimate = unname(estimate),
    CI_Lower = ifelse(available, bounds$lower, NA_real_),
    CI_Upper = ifelse(available, bounds$upper, NA_real_), Level = level,
    Available = available, Status = ifelse(available, "experimental", "unresolved_tail"),
    Planned = nsim, Resolved = bounds$available, Unresolved = nsim - bounds$available,
    ExpectedTailDraws = nsim * (1 - level) / 2, CIEligible = FALSE,
    LowerMCLower = mc_lower$Lower, LowerMCUpper = mc_lower$Upper,
    UpperMCLower = mc_upper$Lower, UpperMCUpper = mc_upper$Upper,
    MCLevel = mc_level, MCFinite = is.finite(mc_lower$Lower) & is.finite(mc_lower$Upper) &
      is.finite(mc_upper$Lower) & is.finite(mc_upper$Upper),
    Uncertainty = if (two_stage) "calibration_and_scoring_refit" else if (nrow(anchors))
      "scoring_refit_conditional_on_fixed_anchors" else "joint_calibration_and_person_refit",
    InferenceReview = "Experimental pointwise fitted-model approximation; finite-sample coverage and JML bias remain unqualified.",
    row.names = NULL)
  if (generator_method == "jeffreys") out$InferenceReview <- paste(out$InferenceReview,
    if (penalized) "Jeffreys estimation throughout original calibration, scoring and refits; this is not ordinary FACETS JML." else
      "Jeffreys generating model; the reported point and refits remain ordinary JML.")
  if (penalized) out$InferenceReview <- paste(out$InferenceReview,
    "Basic refit limits can substantially under-cover in small samples even with exact refit quantiles; see help.")
  attr(out, "draws") <- draws; attr(out, "person_draws") <- person_draws
  attr(out, "monte_carlo") <- mc
  attr(out, "anchor_draws") <- anchor_draws; attr(out, "trials") <- trials
  attr(out, "generating_target") <- setNames(truth, rownames(L))
  if (!is.null(generating_model)) attr(out, "generating_model") <- generating_model
  attr(out, "reference_plan") <- plan; attr(out, "response_map") <- maps
  attr(out, "source_checks") <- source_checks
  attr(out, "source_identity") <- lapply(sources, mfrm_fair_fit_identity)
  attr(out, "settings") <- list(parm = parm, contrasts = L, level = level, nsim = nsim,
    seed = seed, rng_kind = rng_kind, gradient_tol = gradient_tol, mc_level = mc_level,
    anchor_uncertainty = anchor_uncertainty, anchors = anchors,
    method = "basic_fixed_person_parametric_refit", generator_method = generator_method,
    estimator = if (penalized) "jeffreys_canonical_joint" else "ordinary_jml", quantile_type = 1L,
    target = if (is.null(contrasts)) "common_plan_expected_average" else "common_plan_expected_average_contrast",
    shared_responses = if (two_stage) length(intersect(maps$calibration, maps$scoring)) else 0L)
  class(out) <- c("mfrm_fair_intervals", "data.frame")
  out
}

#' @rdname confint.mfrm_fair_average
#' @export
summary.mfrm_fair_intervals <- function(object, ...) {
  rlang::check_dots_empty()
  object
}

#' @rdname confint.mfrm_fair_average
#' @param x For `print.mfrm_fair_intervals()`: A saved `mfrm_fair_intervals` result.
#'
#' For `plot()`: Saved output of [confint.mfrm_fair_average()].
#' @export
print.mfrm_fair_intervals <- function(x, ...) {
  bounded <- identical(attr(x, "settings")$method, "hoeffding_bounded_responses")
  cat(if (bounded) "Common-plan expected-score bounds (independent bounded responses)\n" else
    "Common-plan expected-score intervals (experimental)\n")
  print.data.frame(x, ...)
  cat(if (bounded) "Coverage requires the declared independence and fixed-plan assumptions; inspect Informative and plan support.\n" else
    "Numerical availability is separate from coverage qualification.\n")
  if (identical(attr(x, "settings")$generator_method, "jeffreys"))
    cat(if (identical(attr(x, "settings")$estimator, "jeffreys_canonical_joint"))
      "Jeffreys estimation throughout original calibration, scoring and refits.\n" else
      "Jeffreys generating model; displayed points and refits retain ordinary JML.\n")
  if (identical(attr(x, "settings")$method, "local_fixed_person_delta") &&
      identical(attr(x, "settings")$estimator, "jeffreys_canonical_joint"))
    cat("Analytic Jeffreys two-stage propagation; score bias and coverage remain unqualified.\n")
  if (!is.null(attr(x, "monte_carlo")))
    cat("Endpoint Monte Carlo enclosures measure simulation precision, not sampling coverage; MCFinite alone does not establish precision.\n")
  invisible(x)
}

#' @rdname confint.mfrm_fair_average
#' @details `plot()` method: Plot common-plan expected-score intervals
#' @param draw Logical; draw with base graphics when true.
#' @param show_title,show_notes Logical; display title/cautions, off by default.
#' @return For `plot()`: Invisibly an `mfrm_plot_data` object. Its `data$data` retains all
#'   result rows and uncertainty labels; `data$notes` retains interpretation.
#' @details For `plot()`: Points are common-plan expected scores or declared paired differences.
#'   Bars show only available pointwise bounds under the saved method. Missing
#'   bars are not zero uncertainty. Plotting never recalculates the fit or removes
#'   unresolved rows from the saved result. For refit results the bars show the
#'   empirical interval; endpoint Monte Carlo enclosures are retained in the
#'   returned table and `data$monte_carlo`, not drawn as sampling intervals. See
#'   [confint.mfrm_fair_average()] for
#'   the estimand, assumptions and finite-sample limitations.
#' @examples # Continue the confint.mfrm_fair_average() example:
#' # plot(ci)
#' # p <- plot(ci, draw = FALSE)
#' # p$data$notes
#' @export
plot.mfrm_fair_intervals <- function(x, draw = TRUE, show_title = FALSE, show_notes = FALSE, ...) {
  rlang::check_dots_empty()
  d <- as.data.frame(x); y <- rev(seq_len(nrow(d)))
  note <- if (identical(attr(x, "settings")$method, "hoeffding_bounded_responses"))
    paste("Conservative pointwise bounds require independent bounded responses and a fixed plan.",
      "Points retain fitted estimates; inspect Informative and observed plan support. These bounds do not correct point bias.") else
    if (identical(attr(x, "settings")$method, "local_fixed_person_delta"))
    "Experimental pointwise delta intervals; absent bars mean unavailable information. Coverage is not qualified." else
    "Experimental pointwise expected-score intervals; absent bars mean unresolved tails. Coverage is not qualified."
  if (identical(attr(x, "settings")$generator_method, "jeffreys")) note <- paste(note,
    if (identical(attr(x, "settings")$estimator, "jeffreys_canonical_joint"))
      "Jeffreys estimation throughout original calibration, scoring and refits." else
      "Jeffreys generating model; points and refits retain ordinary JML.")
  if (!is.null(attr(x, "monte_carlo"))) note <- paste(note,
    "Bars show empirical limits; inspect endpoint Monte Carlo enclosures for simulation precision.")
  if (identical(attr(x, "settings")$method, "local_fixed_person_delta") &&
      identical(attr(x, "settings")$estimator, "jeffreys_canonical_joint")) note <- paste(note,
    "Analytic Jeffreys score and anchor derivatives; this is not inverse penalized curvature.")
  xlab <- if (identical(attr(x, "settings")$target, "common_plan_expected_average_contrast"))
    "Expected score contrast (internal units)" else "Expected average (internal score units)"
  out <- new_mfrm_plot_data("fair_reference_intervals", list(data = d,
    settings = attr(x, "settings"), reference_plan = attr(x, "reference_plan"),
    monte_carlo = attr(x, "monte_carlo"),
    plot_data = data.frame(Label = d$Target, X = d$Estimate, Y = y,
      Lower = d$CI_Lower, Upper = d$CI_Upper, Available = d$Available),
    notes = data.frame(Type = "uncertainty", Text = note),
    title = "Common-plan expected-score intervals", xlab = xlab,
    display = list(show_title = show_title, show_notes = show_notes)))
  if (isTRUE(draw)) {
    label_margin <- max(graphics::strwidth(d$Target, units = "inches",
      cex = graphics::par("cex.axis"))) /
      (graphics::par("csi") * graphics::par("mex")) + 2
    old <- graphics::par(mar = c(if (show_notes) 6 else 4.5,
      max(7, label_margin), if (show_title) 3 else 1, 1))
    on.exit(graphics::par(old), add = TRUE)
    graphics::plot(d$Estimate, y, xlim = range(c(d$Estimate, d$CI_Lower, d$CI_Upper), finite = TRUE),
      ylim = c(.5, nrow(d) + .5), yaxt = "n", ylab = "", xlab = out$data$xlab,
      main = if (show_title) out$data$title else "", pch = 19)
    graphics::axis(2, at = y, labels = d$Target, las = 2)
    graphics::segments(d$CI_Lower, y, d$CI_Upper, y)
    if (show_notes) graphics::mtext(note, side = 1, line = 4, cex = .7)
  }
  invisible(out)
}
