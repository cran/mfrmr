#' Parametric bootstrap for GPCM slopes or a matched PCM comparison
#'
#' Simulate Persons and scores on the analyzed rating assignment and refit the
#' model. With no `null_fit`, saved draws support slope intervals. With a matched
#' PCM `null_fit`, simulate under PCM and refit both models for a bootstrap LRT.
#'
#' @param fit A native GPCM MML fit with eligible local information.
#' @param nsim Planned datasets (default 499, minimum 2). Each needs one or two
#'   refits. Small values demonstrate mechanics, not accurate tail probabilities.
#' @param seed Required nonnegative integer. The caller's RNG state is restored.
#' @param null_fit Optional PCM MML fit accepted by [compare_mfrm()] with
#'   `nested = TRUE`. This changes the target to an equal-slope LRT; its draws
#'   must not be used as confidence intervals around the alternative GPCM fit.
#' @param object,x A saved `mfrm_gpcm_bootstrap` result.
#' @param parm Must be `"slopes"`.
#' @param level Nominal confidence level.
#' @param scale,contrasts,contrast_scale,simultaneous As in [confint.mfrm_fit()].
#' @param ... Unused.
#'
#' @details One ability is drawn per Person from the fitted conditional-normal
#'   population, shared by all that Person's rows. Facets and covariates are fixed.
#'   The observed assignment is preserved; omitted ratings remain omitted. No
#'   assignment or missingness mechanism is simulated. Population coefficients
#'   and variance are reestimated if estimated in the source model. Retained
#'   formula, factor coding and person data must reproduce the population design.
#'
#'   `confint()` uses basic bootstrap errors on the log scale for slopes/ratios,
#'   and the identity scale for differences. The point estimate minus reversed
#'   empirical error quantiles (type 1) supplies the limits. Failed refits and
#'   rejected returned estimates are never removed or replaced: their unknown errors are placed at both extremes
#'   to enclose the empirical limits for every completion. Limits can be zero,
#'   infinite or unbounded. `availability` and `expected_tail_draws` describe this
#'   uncertainty. Bonferroni covers only the requested finite family, assuming
#'   adequate marginal bootstrap approximations; it is not a coverage guarantee.
#'
#'   For an LRT, `test` reports `(1 + exceedances)/(nsim + 1)`. If any replicate
#'   test is unresolved, `PValue` is missing and `PValueLower`/`PValueUpper` bound
#'   all completions. Monte Carlo binomial bounds and resolution are also retained.
#'   This is fitted-model calibration, not an exact finite-sample test or a remedy
#'   for misspecification, dependence between Persons, or informative assignment.
#'   Increasing `nsim` improves simulation precision, not the underlying model.
#'   A profile-likelihood method is not used by this function.
#'
#' @section Eligibility and unresolved replicates:
#'   A returned fit is not the same as an accepted replicate. For slope
#'   bootstrapping, a singleton score category (observed once) is a warning,
#'   not by itself a reason to exclude the point estimate. This applies to the
#'   source fit and refits only when a fresh category audit confirms all score
#'   categories are observed in every step scope, with no unsupported step
#'   coordinates or category contrasts. Model identity, numerical convergence,
#'   reevaluated likelihood/gradient and unregularized joint-information checks
#'   must still pass. The saved fit's readiness is not changed.
#'
#'   Basic bootstrap quantiles do not themselves require a variance estimate
#'   for every replicate. The retained information check is a conservative
#'   solution-quality restriction, not a requirement of the basic formula.
#'   Singular information or unstable numerical refinement remains unresolved,
#'   with returned estimates saved for diagnosis. Positive ill-conditioned
#'   information can be used with a caution after refinement, unregularized
#'   inversion and scaled-gradient checks pass. This is distinct from a random-effect
#'   variance estimated at zero; neither proves optimization failed by itself.
#'   Wald intervals, information-criterion comparisons and LRT eligibility
#'   share the numerical information review but retain their other checks.
#'
#'   `checks` separates `BootstrapEligible` from `WaldEligible` and retains
#'   `BootstrapCaution`. Accepted singleton or weak-information cases produce an aggregate warning;
#'   cautions also follow `confint()`, its printed output, default plot subtitle,
#'   APA tables and reports. A custom subtitle (including NULL) overrides the
#'   plot text, not the saved diagnostics. `source_checks` records the source
#'   decision. Use `apa_table(result, which = "checks")` to inspect refits.
#'   `InformationRefinementVerified`, `InformationRelativeChange`,
#'   `InformationInverseResidual` and `InformationScaledGradient` record a
#'   numerical refinement when needed; missing values mean it was not run.
#'   Do not drop unresolved replicates or substitute diagnostic `refit_draws`
#'   for accepted `draws`. A larger `nsim` does not remove an unresolved
#'   estimation problem or guarantee coverage.
#'
#'   `PopulationSD`, `MinimumStandardizedSlope` and `MaximumStandardizedSlope`
#'   describe the returned optimizer solution. Standardized slopes multiply
#'   relative slopes by the fitted population SD (one for a fixed standard-normal
#'   population). With covariates, this is the residual population SD. These
#'   diagnostics help distinguish scale changes from a slope approaching zero;
#'   they do not certify a boundary solution or label rater quality. A missing
#'   field means the required estimate was not retained. A near-zero slope,
#'   an unused score category and an ill-conditioned information matrix can
#'   have different consequences for different parameters: stable slopes do not
#'   establish finite thresholds or valid Wald intervals. Use `checks` to identify
#'   cases needing further numerical review before interpreting an unresolved case.
#'   APA check tables display scale, slope, gradient and information diagnostics
#'   with significant digits so that small positive values are not rounded to
#'   zero. The original `checks` fields remain unrounded numeric values.
#'
#' @return A serializable result containing `source`, optional `null_fit`, all
#'   `trials` with errors/warnings and seeds, `draws` on the free parameter scale,
#'   and `settings`. `refit_draws` additionally retains returned alternative-model
#'   parameter vectors before eligibility checks; rejected rows are diagnostic
#'   values, never inputs to `confint()`. `trials` identifies the last stage and
#'   whether each model refit returned. `checks` and `source_checks` retain category counts/states,
#'   numerical status and the information diagnostics already computed; missing
#'   fields mean not recorded, not a successful check. No extra information
#'   calculation is performed for recording. Older saved results may lack these
#'   fields. LRT results also have `comparison` and `test`. Printing or
#'   changing the interval level never refits. Failed replicates remain present.
#'   `summary()` returns planned/available/unresolved counts, settings and either
#'   the saved likelihood-ratio test or default slope intervals, plus retained
#'   numerical cautions. It does not refit or discard unresolved trials.
#'   If a saved analysis records selected refit updates in `settings$recheck`,
#'   printing and slope intervals retain a caution that it is not a complete
#'   rerun under one procedure. Sampling tables preserve that record and any
#'   `settings$diagnostic_checks_scope`. Missing history means not recorded;
#'   it does not certify that all draws used the currently installed estimator.
#'   After changing estimation or acceptance rules, run a separate complete
#'   bootstrap to evaluate the changed procedure; retain the original result.
#' @seealso [confint.mfrm_fit()], [compare_mfrm()], [mml_quadrature_sensitivity()]
#' @references Chalmers, R. P. (2012). mirt: A Multidimensional Item Response
#'   Theory Package for the R Environment. Journal of Statistical Software,
#'   48(6), 1--29. \doi{10.18637/jss.v048.i06}.
#'
#'   Davison, A. C., and Hinkley, D. V. (1997). Bootstrap Methods and Their
#'   Application. Cambridge University Press, Chapter 5.
#' @examples
#' # After fitting compatible MML models:
#' # boot <- bootstrap_mfrm_gpcm(gpcm_fit, nsim = 499, seed = 92401)
#' # confint(boot, scale = "standardized")
#' # test <- bootstrap_mfrm_gpcm(gpcm_fit, nsim = 499, seed = 92402,
#' #                             null_fit = pcm_fit)
#' # test$test
#' @export
bootstrap_mfrm_gpcm <- function(fit, nsim = 499L, seed, null_fit = NULL) {
  stop_if_product_slopes(fit, "bootstrap_mfrm_gpcm()")
  info <- mfrm_gpcm_inference(fit)
  source_check <- if (is.null(null_fit)) mfrm_gpcm_bootstrap_check(fit, info) else info$check
  if (!isTRUE(source_check$eligible)) stop(source_check$review, call. = FALSE)
  if (!is.numeric(nsim) || is.complex(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim != floor(nsim) || nsim < 2 || nsim > 100000) stop("`nsim` must be an integer from 2 to 100000.", call. = FALSE)
  if (missing(seed) || !is.numeric(seed) || is.complex(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("Supply a nonnegative integer `seed`.", call. = FALSE)
  }
  comparison <- NULL
  if (!is.null(null_fit)) {
    if (!inherits(null_fit, "mfrm_fit") || !identical(null_fit$config$model, "PCM")) stop("`null_fit` must be a native PCM fit.", call. = FALSE)
    comparison <- compare_mfrm(null_fit, fit, nested = TRUE)
    if (!identical(comparison$comparison_basis$lrt_status, "computed")) {
      stop("The source PCM/GPCM comparison is not eligible: ", comparison$comparison_basis$lrt_reason, call. = FALSE)
    }
  }
  generator <- null_fit %||% fit
  # Validate replay availability before planning any random draws.
  mfrmr_gqs_population_arguments(fit)
  if (!is.null(null_fit)) mfrmr_gqs_population_arguments(null_fit)
  if (!is.null(source_check$caution)) warning(source_check$caution, call. = FALSE)
  with_preserved_rng_seed(seed, {
    seeds <- sample.int(.Machine$integer.max, nsim)
    draws <- matrix(NA_real_, nsim, length(fit$opt$par))
    refit_draws <- draws
    checks <- vector("list", nsim)
    trials <- lapply(seq_len(nsim), function(i) {
      warnings <- character(); reason <- ""; value <- NA_real_; ready <- FALSE
      stage <- "generation"; alternative_returned <- null_returned <- FALSE
      trial_checks <- list()
      tryCatch(withCallingHandlers({
        dat <- mfrm_gpcm_bootstrap_generate(generator, seeds[i])
        stage <- "alternative_refit"
        alternative <- mfrm_gpcm_bootstrap_refit(fit, dat)
        alternative_returned <- TRUE
        refit_draws[i, ] <<- alternative$opt$par
        trial_checks[[1L]] <- mfrm_gpcm_bootstrap_fit_checks(alternative, i, "alternative")
        if (is.null(null_fit)) {
          stage <- "slope_eligibility"
          inference <- mfrm_gpcm_inference(alternative)
          check <- mfrm_gpcm_bootstrap_check(alternative, inference)
          trial_checks[[1L]] <- mfrm_gpcm_bootstrap_fit_checks(alternative, i,
            "alternative", inference, check)
          if (!isTRUE(check$eligible)) stop(check$review, call. = FALSE)
          if (!is.null(check$caution)) warning(check$caution, call. = FALSE)
          draws[i, ] <<- alternative$opt$par
        } else {
          stage <- "null_refit"
          restricted <- mfrm_gpcm_bootstrap_refit(null_fit, dat)
          null_returned <- TRUE
          trial_checks[[2L]] <- mfrm_gpcm_bootstrap_fit_checks(restricted, i, "null")
          stage <- "comparison"
          result <- compare_mfrm(restricted, alternative, nested = TRUE)
          if (!identical(result$comparison_basis$lrt_status, "computed")) {
            stop(result$comparison_basis$lrt_reason, call. = FALSE)
          }
          value <- result$lrt$ChiSq
        }
        ready <- TRUE; stage <- "complete"
      }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) reason <<- conditionMessage(e))
      checks[i] <<- list(if (length(trial_checks)) do.call(rbind, trial_checks) else NULL)
      data.frame(Replicate = i, Seed = seeds[i], Available = ready, LR = value,
        Stage = stage, AlternativeReturned = alternative_returned, NullReturned = null_returned,
        Reason = reason, Warnings = paste(unique(warnings), collapse = " | "))
    })
    trials <- do.call(rbind, trials)
    colnames(draws) <- colnames(refit_draws) <- mfrm_checkpoint_parameter_names(info$information$sizes)
    checks <- do.call(rbind, checks)
    if (is.null(checks)) checks <- data.frame()
    out <- list(source = fit, null_fit = null_fit, comparison = comparison,
      draws = draws, refit_draws = refit_draws, trials = trials, checks = checks,
      source_checks = mfrm_gpcm_bootstrap_fit_checks(fit, NA_integer_, "source", info,
        if (is.null(null_fit)) source_check else NULL),
      settings = list(nsim = nsim, seed = seed,
        replicate_seeds = seeds, rng_kind = RNGkind(),
        purpose = if (is.null(null_fit)) "slope_intervals" else "PCM_GPCM_test",
        sampling = "Conditional-normal independent Persons, fixed facets/covariates and analyzed assignment"))
    if (!is.null(null_fit)) out$test <- mfrm_gpcm_bootstrap_test(trials$LR, comparison$lrt$ChiSq)
    class(out) <- "mfrm_gpcm_bootstrap"
    cautioned <- sum(checks$BootstrapEligible %in% TRUE & nzchar(checks$BootstrapCaution))
    if (cautioned > 0L) warning(cautioned,
      " bootstrap refit(s) accepted with numerical or category cautions; inspect `checks` and interval cautions.",
      call. = FALSE)
    out
  })
}

# Reuse the same evaluated information; no second Hessian is calculated.
mfrm_gpcm_bootstrap_check <- function(fit, inference) {
  check <- mfrm_gpcm_slope_inference_check(fit, inference$information, allow_singleton = TRUE)
  if (isTRUE(check$eligible)) check$review <- trimws(paste(
    "Point estimate accepted for the basic slope bootstrap after identity, category,",
    "convergence and unregularized information checks.", check$caution %||% ""))
  check
}

# Read retained diagnostics only; recording these must not add numerical work.
mfrm_gpcm_bootstrap_fit_checks <- function(fit, replicate, role, inference = NULL, bootstrap = NULL) {
  readiness <- mfrmr_get_readiness_record(fit)$fit
  category <- fit$config$category_support_audit
  counts <- category$category_table$WithinScopeCount
  information <- inference$information
  evaluation <- information$solution_information$evaluation_summary
  eigenvalues <- information$solution_information$eigenvalue_summary
  population_sd <- if (isTRUE(fit$config$population_spec$active)) {
    variance <- fit$population$sigma2 %||% NA_real_
    if (length(variance) == 1L && is.finite(variance) && variance > 0) sqrt(variance) else NA_real_
  } else 1
  slopes <- fit$slopes$OptimizerEstimate %||% fit$slopes$Estimate
  standardized <- slopes * population_sd
  data.frame(Replicate = replicate, Role = role, Model = fit$config$model,
    InputState = readiness$InputState[1] %||% NA_character_,
    CategoryState = readiness$CategoryState[1] %||% NA_character_,
    CategoryReason = category$readiness$ReasonCodes[1] %||% NA_character_,
    MinimumCategoryCount = if (length(counts)) min(counts) else NA_integer_,
    EmptyCategories = if (length(counts)) sum(counts == 0L) else NA_integer_,
    SingletonCategories = if (length(counts)) sum(counts == 1L) else NA_integer_,
    NumericalState = readiness$NumericalState[1] %||% NA_character_,
    OptimizerConvergence = fit$opt$convergence %||% NA_integer_,
    Objective = fit$opt$value %||% NA_real_,
    PopulationSD = population_sd,
    MinimumStandardizedSlope = if (length(standardized)) min(standardized) else NA_real_,
    MaximumStandardizedSlope = if (length(standardized)) max(standardized) else NA_real_,
    InformationStatus = information$status %||% NA_character_,
    GradientMaxAbs = evaluation$GradientMaxAbs[1] %||% NA_real_,
    SmallestEigenvalue = eigenvalues$Smallest[1] %||% NA_real_,
    InformationScale = eigenvalues$AbsoluteScale[1] %||% NA_real_,
    InformationRefinementVerified = information$solution_information$inverse_review$Verified %||% NA,
    InformationRelativeChange = information$solution_information$inverse_review$RelativeChange %||% NA_real_,
    InformationInverseResidual = information$solution_information$inverse_review$InverseResidual %||% NA_real_,
    InformationScaledGradient = information$solution_information$inverse_review$CurvatureScaledGradient %||% NA_real_,
    WaldEligible = inference$check$eligible %||% NA,
    InferenceReview = inference$check$review %||% NA_character_,
    BootstrapEligible = bootstrap$eligible %||% NA,
    BootstrapReview = bootstrap$review %||% NA_character_,
    BootstrapCaution = bootstrap$caution %||% "")
}

mfrm_gpcm_bootstrap_generate <- function(fit, seed) {
  with_preserved_rng_seed(seed, {
    config <- fit$config; sizes <- build_param_sizes(config)
    params <- expand_params(fit$opt$par, sizes, config)
    idx <- build_indices(fit$prep, config$step_facet, config$slope_facet, config$interaction_specs)
    pop <- materialize_population_spec(config, params)
    mu <- if (isTRUE(pop$active)) drop(pop$design_matrix[pop$person_lookup, , drop = FALSE] %*% pop$coefficients) else rep(0, config$n_person)
    sd <- if (isTRUE(pop$active)) sqrt(pop$sigma2) else 1
    theta <- stats::rnorm(config$n_person, mu, sd)
    eta <- theta[idx$person] + compute_base_eta(idx, params, config)
    probs <- if (config$model == "RSM") {
      category_prob_rsm(eta, c(0, cumsum(params$steps)))
    } else {
      cumulative <- t(apply(params$steps_mat, 1L, function(x) c(0, cumsum(x))))
      if (config$model == "GPCM") category_prob_gpcm(eta, cumulative, idx$step_idx,
        params$slopes, idx$slope_idx) else category_prob_pcm(eta, cumulative, idx$step_idx)
    }
    u <- stats::runif(nrow(probs)); category <- rowSums(u > t(apply(probs, 1L, cumsum)))
    category <- pmin(category, ncol(probs)-1L)
    map <- fit$prep$score_map
    score <- map$OriginalScore[match(category + fit$prep$rating_min, map$InternalScore)]
    if (anyNA(score)) stop("The source category map cannot reproduce every category.", call. = FALSE)
    replay <- config$replay_inputs
    dat <- as.data.frame(fit$prep$data[c("Person", config$facet_names)])
    dat[] <- lapply(dat, as.character)
    names(dat)[1] <- replay$person
    dat[[replay$score]] <- score
    if (!is.null(replay$weight)) dat[[replay$weight]] <- 1
    dat
  })
}

mfrm_gpcm_bootstrap_refit <- function(source, data) {
  args <- mfrmr_gqs_refit_arguments(source, data, source$config$estimation_control$quad_points)
  fit <- do.call(fit_mfrm, args)
  if (!identical(source$config$step_facet, fit$config$step_facet) ||
      !identical(source$config$slope_facet, fit$config$slope_facet) ||
      !identical(mfrm_population_design(source), mfrm_population_design(fit)) ||
      !identical(build_param_sizes(source$config), build_param_sizes(fit$config)) ||
      !identical(source$prep$levels, fit$prep$levels) ||
      !isTRUE(all.equal(source$prep$score_map, fit$prep$score_map, check.attributes = FALSE))) {
    stop("A bootstrap refit changed the population design, category map or parameter identities.", call. = FALSE)
  }
  fit
}

mfrm_gpcm_bootstrap_test <- function(draws, observed) {
  known <- is.finite(draws); n <- length(draws); exceed <- sum(draws[known] >= observed)
  missing <- sum(!known)
  lower <- (1 + exceed)/(n+1); upper <- (1 + exceed + missing)/(n+1)
  mc_lower <- stats::binom.test(exceed, n)$conf.int[1]
  mc_upper <- stats::binom.test(exceed + missing, n)$conf.int[2]
  data.frame(LR = observed, Planned = n, Available = sum(known), Unresolved = missing,
    Exceedances = exceed, PValue = if (missing) NA_real_ else lower,
    PValueLower = lower, PValueUpper = upper, MonteCarloLower = mc_lower,
    MonteCarloUpper = mc_upper, Resolution = 1/(n+1))
}

#' @rdname bootstrap_mfrm_gpcm
#' @export
confint.mfrm_gpcm_bootstrap <- function(object, parm = "slopes", level = .95,
    scale = c("relative", "standardized"), contrasts = NULL,
    contrast_scale = c("ratio", "difference"), simultaneous = c("none", "bonferroni"), ...) {
  rlang::check_dots_empty(); mfrm_random_rater_interval_level(level)
  scale <- match.arg(scale); contrast_scale <- match.arg(contrast_scale); simultaneous <- match.arg(simultaneous)
  if (!identical(parm, "slopes")) stop("Use parm = 'slopes'.", call. = FALSE)
  if (!identical(object$settings$purpose, "slope_intervals")) stop("Null-model bootstrap draws are for the LRT, not slope intervals.", call. = FALSE)
  target <- mfrm_gpcm_slope_target(object$source, scale, contrasts, contrast_scale)
  roots <- matrix(NA_real_, nrow(object$draws), length(target$value))
  for (i in which(object$trials$Available)) {
    f <- object$source; f$opt$par <- object$draws[i, ]
    roots[i, ] <- mfrm_gpcm_slope_target(f, scale, contrasts, contrast_scale)$value - target$value
  }
  missing <- !is.finite(roots); low <- high <- roots
  low[missing] <- -Inf; high[missing] <- Inf
  tail <- (1-level)/(2 * if (simultaneous == "bonferroni") ncol(roots) else 1)
  lower <- target$value - apply(high, 2, stats::quantile, probs = 1-tail, type = 1, names = FALSE)
  upper <- target$value - apply(low, 2, stats::quantile, probs = tail, type = 1, names = FALSE)
  if (target$log_scale) { lower <- exp(lower); upper <- exp(upper) }
  cautions <- unique(c(if (!is.null(object$settings$recheck))
      "Saved result includes selected refit updates; it is not a complete rerun under one procedure.",
    object$source_checks$BootstrapCaution,
    object$checks$BootstrapCaution[object$checks$BootstrapEligible %in% TRUE]))
  cautions <- cautions[!is.na(cautions) & nzchar(cautions)]
  tab <- data.frame(SlopeFacet = target$labels, Estimate = target$estimate, CI_Lower = lower,
    CI_Upper = upper, CIEligible = is.finite(lower) & is.finite(upper) & (!target$log_scale | lower > 0),
    InferenceReview = paste(c("Basic fitted-model bootstrap; unresolved draws widen limits, possibly without finite bounds.",
      cautions), collapse = " "))
  out <- mfrm_gpcm_interval_result(tab, level, "Basic parametric bootstrap", target$target, simultaneous)
  attr(out, "cautions") <- cautions
  attr(out, "availability") <- data.frame(Target = target$labels, Planned = nrow(roots),
    Available = colSums(!missing), Unresolved = colSums(missing))
  attr(out, "expected_tail_draws") <- nrow(roots)*tail
  attr(out, "settings") <- list(scale = scale, contrasts = target$contrasts,
    contrast_scale = contrast_scale, seed = object$settings$seed,
    recheck = object$settings$recheck,
    diagnostic_checks_scope = object$settings$diagnostic_checks_scope)
  attr(out, "source") <- mfrm_gpcm_inference_source(object$source)
  out
}

#' @rdname bootstrap_mfrm_gpcm
#' @export
summary.mfrm_gpcm_bootstrap <- function(object, ...) {
  rlang::check_dots_empty()
  available <- sum(object$trials$Available %in% TRUE)
  cautions <- unique(object$checks$BootstrapCaution[object$checks$BootstrapEligible %in% TRUE])
  list(overview = data.frame(Planned = nrow(object$trials), Available = available,
      Unresolved = nrow(object$trials) - available),
    settings = object$settings, test = object$test,
    intervals = if (is.null(object$test)) confint(object) else NULL,
    cautions = cautions[!is.na(cautions) & nzchar(cautions)])
}

#' @rdname bootstrap_mfrm_gpcm
#' @export
print.mfrm_gpcm_bootstrap <- function(x, ...) {
  purpose <- if (identical(x$settings$purpose, "slope_intervals"))
    "slope intervals" else "PCM versus GPCM likelihood-ratio test"
  cat("GPCM parametric bootstrap:", purpose, "\n")
  cat("Planned:", nrow(x$trials), " Available:", sum(x$trials$Available), "\n")
  if (!is.null(x$test) && !is.null(x$settings$recheck))
    cat("Saved result includes selected refit updates; it is not a complete rerun under one procedure.\n")
  if ("AlternativeReturned" %in% names(x$trials)) {
    cat("Alternative refits returned:", sum(x$trials$AlternativeReturned),
      " | Unresolved trials with a returned alternative:",
      sum(x$trials$AlternativeReturned & !x$trials$Available), "\n")
  }
  if (!is.null(x$test)) print(x$test, row.names = FALSE) else print(confint(x))
  cat("Fitted-model approximation; unresolved replicates retained and finite-sample accuracy not guaranteed.\n")
  invisible(x)
}
