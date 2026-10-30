#' Fit a rating-scale model with dependence within Person-specific testlets
#'
#' Account for extra dependence among ratings of the same person, for example
#' several rubric criteria from one performance. A testlet identifies the
#' ratings sharing this local effect within that person. The rating-scale
#' model (RSM) estimates fixed facets, normal ability variance and one common
#' local variance. For a rater effect shared across people instead, see
#' [fit_mfrm_random_rater()].
#' @param data Long-format ratings; each row is one assigned rating.
#' @param person,score,testlet Column names. `testlet` identifies a group within
#'   a Person and may also occur in `facets`, for example a fixed Rater effect.
#' @param facets Fixed-facet column names; default none. A column used as
#'   `testlet` is not automatically a fixed facet. For example,
#'   `testlet = "Task"` models local dependence; also specify `facets = "Task"`
#'   if the model should estimate fixed task difficulties.
#' @param score_levels Consecutive integer categories in increasing order.
#' @param testlet_variance `NULL` estimates the common local variance. A
#'   nonnegative number fixes it as known; zero removes local dependence.
#' @param person_sd `NULL` (default) estimates the SD of a mean-zero normal
#'   ability population. A positive number fixes the SD as known. Use `1` to
#'   reproduce the earlier fixed N(0,1) model; this is a population assumption,
#'   not merely a choice of units.
#' @param person_variance_max Upper search bound for estimated ability
#'   variance (not SD); default 16. Reaching it withholds numerical readiness.
#' @param quad_points Standard-normal quadrature order for both integrals,
#'   from 7 to 121; default 31. Results are checked at `2 * quad_points + 1`.
#' @param maxit Maximum iterations per optimization start; default 300.
#' @param variance_max Upper search bound for an estimated variance; default
#'   16. A fit reaching this bound is not numerically ready.
#' @param missing `"fail"` (default) refuses missing assigned scores;
#'   `"omit"` analyzes observed rows and records the omissions.
#' @details The adjacent-category logit is
#'   \deqn{\log\{P(Y_{pbi}=k)/P(Y_{pbi}=k-1)\}
#'      = \theta_p + \gamma_{pb} - x_{pbi}'\beta - \delta_k.}
#'   Abilities are independent N(0,\eqn{\sigma_p^2}). Local effects are independent N(0,v),
#'   independent of abilities and assignment. A label reused by another Person
#'   represents a new local effect, not a shared random rater. Every row has one
#'   non-overlapping membership. Fixed-facet effects sum to zero; free steps
#'   determine overall location. There are no calibration priors.
#'
#'   Frequentist MML uses nested one-dimensional Gaussian quadrature: integrate
#'   each local effect conditional on ability, multiply its block likelihoods,
#'   then integrate ability. Unequal block sizes and incomplete assignments are
#'   supported. At least two Persons must have multiple observed testlets and
#'   at least two Persons must have a testlet with repeated ratings. These
#'   initial design requirements do not prove estimability. Unit weights, RSM,
#'   additive fixed facets and one common variance are the supported scope;
#'   overlapping memberships, PCM, correlated effects, random rater effects
#'   shared across Persons, covariate-dependent ability means and nonnormal
#'   ability populations are not included. Estimating one ability variance
#'   does not establish homogeneous populations across assignment groups.
#'   After integrating both random effects, different Persons are independent
#'   under the model; ratings within a Person are not. Regular likelihood
#'   inference requires an identified interior model and growing information
#'   across independent Persons. The full calibration information includes
#'   every estimated variance; treating an estimated variance as known would
#'   omit nuisance-parameter uncertainty. A positive numerical Hessian alone
#'   does not prove identification, and fixed-order quadrature agreement does
#'   not establish an asymptotic approximation rate.
#'
#'   Zero is evaluated explicitly for both estimated variances using their
#'   one-sided variance scores, including the joint-zero submodel. All
#'   optimization starts are retained. Fixed coordinates are searched within
#'   `[-20,20]`; reaching a bound withholds readiness. Integration differences,
#'   the projected gradient and information are separate checks.
#'   A converged start is preferred when its objective differs from the best
#'   failed start only within floating-point rounding precision.
#'   Increasing quadrature order may be necessary; no order is universally
#'   sufficient.
#'   Calibration tables retain estimates and approximate SEs but omit bounds
#'   by default. `confint(fit, parm = "calibration", level = 0.95)` or
#'   `summary(fit, calibration_intervals = "normal", level = 0.95)` explicitly
#'   requests observed-information normal approximations. Their finite-sample
#'   coverage is not established. Numerical/information failures and estimated
#'   variance boundaries retain missing bounds. A regular variance interval is
#'   not supplied. Rebuilding summaries or [mfrm_results()] from older saved
#'   fits applies this policy without refitting or changing the source fit.
#'   If estimated ability variance is zero, Person scoring returns unavailable
#'   rows rather than degenerate zero-width intervals. Positive information for
#'   the remaining fixed coordinates does not resolve variance-boundary inference.
#'   [predict.mfrm_testlet()] supplies conditional Person scores and continuous
#'   equal-tail intervals; these exclude calibration-estimation uncertainty.
#'   Numerical agreement does not establish coverage or model fit.
#'
#' @section Model assumptions and related research:
#' Wang and Wilson's Rasch testlet model (2005) motivates dependence within
#' a Person-specific block. With `testlet = "Rater"` and a fixed Rater facet,
#' the more direct reference is their random-effects facet model (2005,
#' equations 14--15): usual rater severity is distinct from an effect shared
#' only within a Person/rater pair. Here the local variance is common to all
#' testlets; it cannot rank raters by individual inconsistency.
#'
#' Both papers allow ability variance to be estimated. This function estimates
#' one normal ability variance by default while fixing its mean at zero and
#' keeping free step location. This avoids adding a second location parameter.
#' Unlike models with separate local variances, it retains a common testlet
#' variance. The cited results do not guarantee coverage for a new assessment
#' design or for this implementation.
#'
#' @section Relation to discrimination and shared rater effects:
#' Discrimination changes response sensitivity; it is not a substitute for
#' dependence within a performance. This function fits an RSM, not a GPCM
#' with additional slopes. Its local effect is shared only within a
#' Person/block pair. A Person/task block shared across assessors differs
#' from a Person/task/assessor rating-occasion block. Select the identifier
#' according to the assessment process; only one nonoverlapping membership
#' is supported. A nonzero local variance does not by itself establish halo.
#' For one rater effect shared across persons, see [fit_mfrm_random_rater()].
#' The two models cannot currently be combined. See [gpcm_capability_matrix()]
#' for the relationship to fixed discrimination and local independence.
#'
#' @section Comparison with ordinary MFRM:
#' Hold the observed events, categories, fixed facets, omissions and ability
#' population constant. A known `testlet_variance = 0` removes local
#' dependence while retaining the fixed Rater facet, if supplied. The
#' corresponding ordinary model is an RSM MML with a matching normal
#' population and location convention. For `person_sd = 1` this is the
#' fixed N(0,1) model; estimated-population comparisons must also align location
#' constraints rather than compare raw coefficients. An estimated zero
#' variance lies on a boundary; ordinary
#' chi-squared likelihood-ratio calibration cannot be assumed.
#' [compare_mfrm()] checks matched events and compares centered facet effects
#' with an ordinary RSM MML fit; it does not provide automatic ranking. Descriptive
#' score or interval changes do not establish which model is better.
#' @section Assessment applications:
#' With `testlet = "Task"`, the local effect is shared by a Person's ratings
#' on that task. `vignette("mfrmr-testlet-applications")` provides runnable
#' examples of allocating a fixed rating budget across tasks and comparing
#' score sensitivity when tasks have unequal numbers of criteria. Its small
#' enumeration reports model-conditional posterior precision and marginal EAP
#' reliability while holding fitted calibration fixed; it is not a general
#' design optimizer, a G/D study or a coverage guarantee. Ordinary MFRM does
#' not generally assign equal information to arbitrary added criteria or tasks.
#'
#' A common local variance cannot identify halo, establish that criteria are
#' indistinguishable, or rank tasks by dependence. Task-specific variances are
#' not estimated by this API. Fitting tasks separately does not remedy this:
#' with one local block per Person, ability and local variances are confounded.
#' Unequal block sizes are handled through the likelihood, not by imposing equal
#' task weights or removing every source of bias. Assignment, content and the
#' origin of dependence require substantive review.
#' @section Earlier saved results:
#' Saved fits without an ability-variance parameter retain their original
#' N(0,1) interpretation when scored or reported. Refitting now estimates
#' ability variance unless `person_sd = 1` is supplied. Reprinting or
#' rescoring does not refit a population. Save a new fit and regenerate its
#' scores together when changing the population assumption.
#' @return An `mfrm_testlet` object with calibration, covariance, calibration
#'   table, checks, all optimization runs, input/omission accounting and settings.
#'   It contains no native pointers; use `saveRDS()` and `readRDS()`.
#'   Calibration tables retain `IntervalStatus`, `CIEligible` and `InferenceReview`
#'   for missing or explicitly requested normal bounds; see
#'   [confint.mfrm_testlet()] for the status meanings. A retained point is not
#'   certified by interval availability. Inspect `checks` before interpreting it.
#' @param object An `mfrm_testlet` result.
#' @param ... Unused.
#' @references Wang, W.-C. and Wilson, M. (2005). The Rasch testlet model.
#'   *Applied Psychological Measurement*, 29, 126--149.
#'   \doi{10.1177/0146621604271053}.
#'
#'   Wang, W.-C. and Wilson, M. (2005). Exploring local item dependence using
#'   a random-effects facet model. *Applied Psychological Measurement*, 29,
#'   296--318. \doi{10.1177/0146621605276281}.
#' @seealso [predict.mfrm_testlet()], [plot.mfrm_testlet()], [mfrm_results()],
#'   [mfrm_response_diagnostics()] for descriptive posterior predictive residuals.
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("fit_mfrm_testlet", run.dontrun = TRUE)
#' \dontrun{
#' # Saved fit for the synthetic example_core ratings.
#' example <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr"))
#' fit <- example$testlet$fit
#' # To refit instead (this takes longer than inspecting the saved result):
#' # ratings <- load_mfrmr_data("example_core")
#' # fit <- fit_mfrm_testlet(ratings, "Person", "Score", "Rater",
#' #   facets = c("Rater", "Criterion"), score_levels = 1:4, quad_points = 121)
#' fit$checks
#' plot(fit, facet = "Rater")
#' # The complete regeneration recipe is included with the package:
#' system.file("examples", "extended-models.R", package = "mfrmr")
#' }
#' @export
#' @order 1
fit_mfrm_testlet <- function(data, person, score, testlet, facets = character(),
    score_levels, testlet_variance = NULL, quad_points = 31L, maxit = 300L,
    variance_max = 16, missing = c("fail", "omit"), person_sd = NULL,
    person_variance_max = 16) {
  if (!is.null(testlet_variance) && (!is.numeric(testlet_variance) || is.complex(testlet_variance) ||
      length(testlet_variance) != 1L || !is.finite(testlet_variance) || testlet_variance < 0)) {
    stop("`testlet_variance` must be NULL or a finite nonnegative known variance.", call. = FALSE)
  }
  if (!is.null(person_sd) && (!is.numeric(person_sd) || is.complex(person_sd) ||
      length(person_sd) != 1L || !is.finite(person_sd) || person_sd <= 0 ||
      !is.finite(person_sd^2) || person_sd^2 == 0)) {
    stop("`person_sd` must be NULL or a finite positive known SD with a finite squared value.", call. = FALSE)
  }
  for (arg in c("quad_points", "maxit", "variance_max", "person_variance_max")) {
    z <- get(arg)
    if (!is.numeric(z) || is.complex(z) || length(z) != 1L || !is.finite(z) || z <= 0 ||
        (arg %in% c("quad_points", "maxit") && (z != floor(z) || z < if (arg == "quad_points") 7 else 1)) ||
        (arg == "quad_points" && z > 121)) stop("Use integer quad_points from 7 to 121, positive integer maxit and positive variance bounds.", call. = FALSE)
  }
  input <- mfrm_testlet_data(data, person, score, testlet, facets, score_levels, match.arg(missing))
  fit <- mfrm_testlet_fit(input, quad_points, maxit, testlet_variance, variance_max,
    if (is.null(person_sd)) NULL else person_sd^2, person_variance_max)
  fit$input <- input
  fit$settings <- list(model = "RSM", method = "MML: nested Gaussian quadrature",
    person_distribution = "N(0, person_variance)", fixed_person_sd = person_sd,
    person_variance_max = person_variance_max, local_distribution = "N(0, testlet_variance)",
    fixed_variance = testlet_variance, quad_points = quad_points, check_points = 2 * quad_points + 1,
    maxit = maxit, variance_max = variance_max, fixed_coordinate_bounds = c(-20, 20),
    missing = match.arg(missing))
  nb <- ncol(input$X); ns <- length(score_levels) - 1L
  fit$calibration <- list(beta = fit$parameters[seq_len(nb)],
    steps = fit$parameters[nb + seq_len(ns)], variance = fit$parameters[nb + ns + 1L],
    person_variance = fit$parameters[nb + ns + 2L], person_sd = sqrt(fit$parameters[nb + ns + 2L]))
  tables <- list(); at <- 0L
  regular <- fit$checks$NumericalReady && fit$checks$InformationPositive &&
    !fit$checks$EstimatedVarianceBoundary && !fit$checks$EstimatedPersonVarianceBoundary
  for (f in facets) {
    b <- input$basis[[f]]; index <- at + seq_len(ncol(b)); at <- max(index)
    se <- if (regular) sqrt(pmax(0, diag(b %*% fit$covariance[index, index, drop = FALSE] %*% t(b)))) else rep(NA_real_, nrow(b))
    tables[[length(tables) + 1L]] <- data.frame(Parameter = "Fixed facet", Facet = f,
      Level = rownames(b), Estimate = as.vector(b %*% fit$parameters[index]), SE = se)
  }
  tables[[length(tables) + 1L]] <- data.frame(Parameter = "Step", Facet = score,
    Level = as.character(score_levels[-1]), Estimate = fit$calibration$steps,
    SE = if (regular) sqrt(diag(fit$covariance)[nb + seq_len(ns)]) else NA_real_)
  fit$calibration_table <- do.call(rbind, tables); rownames(fit$calibration_table) <- NULL
  fit$calibration_table <- mfrm_extended_calibration_table(fit)
  fit$call <- match.call(); class(fit) <- "mfrm_testlet"
  if (!fit$checks$NumericalReady || !fit$checks$InformationPositive) warning(
    "Numerical or information checks require review; intervals and scoring are unavailable. Inspect $checks.", call. = FALSE)
  fit
}

#' @rdname fit_mfrm_testlet
#' @param calibration_intervals For summaries, `"none"` (default) omits
#'   calibration bounds; `"normal"` explicitly requests their pointwise normal
#'   approximation. Person scoring remains a separate conditional output.
#' @param level For `summary.mfrm_testlet()`: Nominal level for explicitly requested summary calibration
#'   intervals; default 0.95. Numerical and variance-boundary restrictions remain.
#'
#' For `plot()`: Nominal level of explicitly requested calibration intervals;
#'   default 0.95. This does not establish finite-sample coverage.
#' @export
summary.mfrm_testlet <- function(object, ..., calibration_intervals = c("none", "normal"), level = .95) {
  rlang::check_dots_empty()
  calibration_intervals <- match.arg(calibration_intervals)
  list(calibration = mfrm_extended_calibration_table(object, calibration_intervals, level), testlet_variance = object$calibration$variance,
    calibration_intervals = list(method = calibration_intervals, level = level),
    person_variance = object$calibration$person_variance %||% 1,
    checks = object$checks, data_usage = c(Input = object$input$input_rows,
      Observed = length(object$input$y), Omitted = length(object$input$omitted_rows)),
    settings = object$settings,
    notes = c(if (calibration_intervals == "none") "Calibration bounds are omitted by default; SEs are observed-information approximations. Request calibration_intervals = 'normal' explicitly to inspect pointwise bounds." else
      "Explicit calibration intervals use an observed-information normal approximation, not established finite-sample coverage.",
      "Estimated variance boundaries withhold regular calibration intervals; no regular variance interval is supplied."))
}

#' @rdname fit_mfrm_testlet
#' @export
print.mfrm_testlet <- function(x, ...) {
  cat("Person-specific testlet RSM (MML)\n")
  cat("Persons:", length(x$input$persons), " Observed blocks:", nrow(x$input$blocks),
    " Local variance:", format(x$calibration$variance, digits = 4),
    " Ability variance:", format(x$calibration$person_variance %||% 1, digits = 4), "\n")
  print(mfrm_extended_calibration_table(x), row.names = FALSE)
  cat("Calibration bounds omitted; SEs are observed-information approximations.\n")
  if (x$checks$EstimatedVarianceBoundary) cat("Estimated zero variance: regular calibration intervals withheld.\n")
  if (isTRUE(x$checks$EstimatedPersonVarianceBoundary)) cat("Estimated zero ability variance: regular intervals and Person scores withheld.\n")
  if (!x$checks$NumericalReady || !x$checks$InformationPositive) cat("Checks require review; scoring is unavailable.\n")
  invisible(x)
}

mfrm_testlet_data <- function(data, person, score, testlet, facets, score_levels, missing, reference = NULL) {
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) stop("Supply a nonempty data frame with distinct column names.", call. = FALSE)
  roles <- list(person, score, testlet)
  if (any(!vapply(roles, function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x), logical(1))) ||
      length(unique(c(person, score, testlet))) != 3L || !is.character(facets) || anyNA(facets) ||
      anyDuplicated(facets) || any(facets %in% c(person, score))) stop("Declare distinct person/score/testlet columns and distinct fixed facets; testlet may also be a fixed facet.", call. = FALSE)
  columns <- unique(c(person, score, testlet, facets))
  if (!all(columns %in% names(data))) stop("Required rating columns are missing.", call. = FALSE)
  if (!is.numeric(score_levels) || is.complex(score_levels) || length(score_levels) < 2L ||
      any(!is.finite(score_levels)) || any(score_levels != floor(score_levels)) || any(diff(score_levels) != 1)) {
    stop("Declare consecutive integer score categories in increasing order.", call. = FALSE)
  }
  d <- data[columns]
  for (id in unique(c(person, testlet, facets))) {
    v <- d[[id]]
    if (!is.atomic(v) || !is.null(dim(v)) || anyNA(v) || any(!nzchar(trimws(as.character(v))))) stop("Identifiers must be complete and nonempty, including on missing-score rows.", call. = FALSE)
    d[[id]] <- as.character(v)
  }
  persons <- sort(unique(d[[person]])); labels <- sort(unique(d[[testlet]]))
  y <- d[[score]]
  if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) || any(!is.na(y) & !y %in% score_levels)) stop("Scores must be numeric declared categories or NA.", call. = FALSE)
  absent <- is.na(y)
  if (any(absent) && missing == "fail") stop("Assigned scores are missing; review them or explicitly use missing = 'omit'.", call. = FALSE)
  levels <- if (is.null(reference)) lapply(d[!absent, facets, drop = FALSE], function(x) sort(unique(x))) else reference$levels
  for (f in facets) if (any(!d[[f]] %in% levels[[f]])) stop("A fixed facet contains an unknown or entirely unobserved level.", call. = FALSE)
  basis <- if (!is.null(reference)) reference$basis else lapply(levels, function(x) {
    if (length(x) < 2L) stop("Each fixed facet needs at least two observed levels.", call. = FALSE)
    b <- qr.Q(qr(stats::contr.sum(length(x)))); dimnames(b) <- list(x, paste0("C", seq_len(ncol(b)))); b
  })
  assigned <- d
  d <- d[!absent, , drop = FALSE]; y <- match(d[[score]], score_levels) - 1L
  X <- if (length(facets)) do.call(cbind, lapply(facets, function(f) basis[[f]][match(d[[f]], levels[[f]]), , drop = FALSE])) else matrix(numeric(), nrow(d), 0L)
  p <- match(d[[person]], persons); t <- match(d[[testlet]], labels)
  groups <- lapply(seq_along(persons), function(i) split(which(p == i), t[p == i]))
  blocks <- do.call(rbind, lapply(seq_along(groups), function(i) {
    if (!length(groups[[i]])) return(NULL)
    data.frame(Person = persons[i], Testlet = labels[as.integer(names(groups[[i]]))],
      Observed = lengths(groups[[i]]), row.names = NULL)
  }))
  if (is.null(blocks)) blocks <- data.frame(Person = character(), Testlet = character(), Observed = integer())
  if (is.null(reference)) {
    if (length(unique(p)) < 2 || !all(score_levels %in% d[[score]])) stop("Fitting requires multiple observed Persons and every declared category.", call. = FALSE)
    if (sum(lengths(groups) >= 2) < 2 || sum(vapply(groups, function(g) any(lengths(g) >= 2), logical(1))) < 2) stop("Fitting requires multiple testlets and repeated within-testlet ratings in at least two Persons.", call. = FALSE)
    if (qr(cbind(1, X))$rank < ncol(X) + 1) stop("Fixed-facet effects are aliased in the observed design.", call. = FALSE)
  }
  # Integer level codes avoid collisions between user labels containing delimiters.
  key <- if (length(facets)) do.call(paste, c(lapply(facets, function(f) match(d[[f]], levels[[f]])), sep = ":")) else rep("1", nrow(d))
  designs <- unique(key); design <- match(key, designs); design_X <- X[match(designs, key), , drop = FALSE]
  nd <- length(designs); nc <- length(score_levels)
  cells <- lapply(groups, function(g) lapply(g, function(rows) {
    counts <- table(factor(design[rows] + y[rows] * nd, levels = seq_len(nd * nc)))
    data.frame(Cell = which(counts > 0), Count = as.integer(counts[counts > 0]))
  }))
  list(data = d, assigned_data = assigned, input_rows = nrow(data), omitted_rows = which(absent), persons = persons,
    columns = list(person = person, score = score, testlet = testlet, facets = facets),
    levels = levels, basis = basis, X = X, y = y, groups = groups, cells = cells,
    design_X = design_X, blocks = blocks, score_levels = score_levels)
}

mfrm_testlet_logsum_rows <- function(x) {
  largest <- do.call(pmax, lapply(seq_len(ncol(x)), function(j) x[, j]))
  largest + log(rowSums(exp(x - largest)))
}

# Factorization generalizes the independently checked two-rater calculation.
mfrm_testlet_kernel <- function(input, par, theta, rule, gradient = FALSE, persons = seq_along(input$persons)) {
  nb <- ncol(input$X); ns <- length(input$score_levels) - 1L; v <- par[nb + ns + 1L]
  if (!is.finite(v) || v < 0) stop("Local variance must be nonnegative.", call. = FALSE)
  gamma <- if (v == 0) 0 else sqrt(v) * rule$nodes
  logweights <- if (v == 0) 0 else log(rule$weights)
  nq <- length(theta); ng <- length(gamma); np <- length(persons); npar <- nb + ns + 1L
  latent <- as.vector(outer(theta, gamma, "+")); steps <- par[nb + seq_len(ns)]
  nd <- nrow(input$design_X); cache <- vector("list", nd * (ns + 1L))
  designs <- seq_len(nd)
  if (nq > 0L && length(persons) < length(input$persons)) {
    cells <- unlist(lapply(persons, function(p)
      unlist(lapply(input$cells[[p]], `[[`, "Cell"), use.names = FALSE)), use.names = FALSE)
    # Cell = design + category * nd; retain the original design order.
    designs <- sort(unique((cells - 1L) %% nd + 1L))
  }
  for (j in designs) {
    eta <- latent - sum(input$design_X[j, ] * par[seq_len(nb)])
    logw <- outer(eta, 0:ns) - matrix(c(0, cumsum(steps)), length(eta), ns + 1L, byrow = TRUE)
    logw <- logw - row_max_fast(logw)
    den <- log(rowSums(exp(logw)))
    if (gradient) {
      prob <- exp(logw - den)
      mean <- as.vector(prob %*% (0:ns))
      variance <- mfrm_category_variance(prob, 0:ns, mean)
    }
    for (y in 0:ns) {
      residual <- if (gradient) y - mean else NULL
      g <- if (gradient) cbind(-outer(residual, input$design_X[j, ]),
        matrix(vapply(seq_len(ns), function(k) rowSums(prob[, (k + 1L):(ns + 1L), drop = FALSE]) - as.numeric(y >= k), numeric(length(eta))), length(eta), ns),
        if (v > 0) residual * rep(gamma / (2 * v), each = nq) else rep(0, length(eta))) else NULL
      cache[[j + y * nd]] <- list(logp = matrix(logw[, y + 1L] - den, nq, ng),
        g = g, residual = if (gradient) matrix(residual, nq, ng),
        variance = if (gradient) matrix(variance, nq, ng))
    }
  }
  loglik <- matrix(0, nq, np); scores <- if (gradient) array(0, c(nq, np, npar)) else NULL
  theta_score <- theta_curvature <- if (gradient) matrix(0, nq, np) else NULL
  for (pi in seq_along(persons)) for (block in input$cells[[persons[pi]]]) {
    lc <- S <- V <- matrix(0, nq, ng)
    g <- if (gradient) matrix(0, nq * ng, npar) else NULL
    for (j in seq_len(nrow(block))) {
      cell <- cache[[block$Cell[j]]]; count <- block$Count[j]
      lc <- lc + count * cell$logp
      if (gradient) {
        g <- g + count * cell$g
        S <- S + count * cell$residual; V <- V + count * cell$variance
      }
    }
    offset <- row_max_fast(lc)
    lw <- sweep(lc - offset, 2L, logweights, "+")
    marginal <- mfrm_testlet_logsum_rows(lw)
    loglik[, pi] <- loglik[, pi] + (offset + marginal)
    if (gradient) {
      if (v == 0) g[, npar] <- as.vector(.5 * (S^2 - V))
      weights <- exp(lw - marginal)
      for (a in seq_len(npar)) scores[, pi, a] <- scores[, pi, a] + rowSums(weights * matrix(g[, a], nq, ng))
      block_score <- rowSums(weights * S)
      theta_score[, pi] <- theta_score[, pi] + block_score
      theta_curvature[, pi] <- theta_curvature[, pi] + rowSums(weights * (S^2 - V)) - block_score^2
    }
  }
  list(loglik = loglik, score = scores, theta_score = theta_score, theta_curvature = theta_curvature)
}

mfrm_testlet_person_variance <- function(input, par) {
  at <- ncol(input$X) + length(input$score_levels) + 1L
  # Older saved fits contain only fixed coordinates and the local variance.
  if (length(par) == at - 1L) return(1)
  if (length(par) != at || !is.finite(par[at]) || par[at] < 0) stop("Invalid testlet ability variance.", call. = FALSE)
  par[at]
}

mfrm_testlet_evaluate <- function(input, par, rule) {
  person_variance <- mfrm_testlet_person_variance(input, par)
  theta <- sqrt(person_variance) * rule$nodes
  kernel <- mfrm_testlet_kernel(input, par, theta, rule, gradient = TRUE)
  offset <- row_max_fast(t(kernel$loglik))
  lw <- sweep(sweep(kernel$loglik, 2L, offset, "-"), 1L, log(rule$weights), "+")
  normalizer <- mfrm_testlet_logsum_rows(t(lw)); posterior <- exp(sweep(lw, 2L, normalizer, "-"))
  gradient <- vapply(seq_len(dim(kernel$score)[3]), function(a) sum(posterior * kernel$score[, , a]), numeric(1))
  if (length(par) > length(gradient)) {
    score <- if (person_variance > 0) kernel$theta_score * theta / (2 * person_variance) else
      .5 * (kernel$theta_score^2 + kernel$theta_curvature)
    gradient <- c(gradient, sum(posterior * score))
  }
  means <- colSums(theta * posterior)
  sds <- sqrt(colSums(sweep(matrix(theta, length(theta), ncol(posterior)), 2L, means, "-")^2 * posterior))
  normalizer <- offset + normalizer
  list(loglik = sum(normalizer), gradient = gradient, person_loglik = normalizer,
    moments = cbind(Estimate = means, ConditionalSD = sds))
}

mfrm_testlet_select_run <- function(runs) {
  values <- vapply(runs, `[[`, numeric(1), "value")
  at <- which.min(values)
  if (runs[[at]]$convergence != 0L && is.finite(values[at])) {
    slack <- 64 * .Machine$double.eps * max(1, abs(values[at]))
    eligible <- which(vapply(runs, function(x) isTRUE(x$convergence == 0L), logical(1)) &
      is.finite(values) & values <= values[at] + slack)
    if (length(eligible)) at <- eligible[which.min(values[eligible])]
  }
  at
}

mfrm_testlet_fit <- function(input, order, maxit, fixed, variance_max,
    fixed_person = 1, person_variance_max = 16) {
  rule <- gauss_hermite_normal(order); high <- gauss_hermite_normal(2 * order + 1L)
  nb <- ncol(input$X); ns <- length(input$score_levels) - 1L; nv <- nb + ns + 1L
  np <- nv + 1L
  start <- c(rep(0, nb), seq(-.5, .5, length.out = ns), 0,
    if (is.null(fixed_person)) min(1, person_variance_max / 2) else fixed_person)
  optimize <- function(start, fixed_v, fixed_p) {
    fixed_at <- integer()
    if (!is.null(fixed_v)) { start[nv] <- fixed_v; fixed_at <- c(fixed_at, nv) }
    if (!is.null(fixed_p)) { start[np] <- fixed_p; fixed_at <- c(fixed_at, np) }
    free <- setdiff(seq_len(np), fixed_at)
    expand <- function(p) {
      z <- start; z[free] <- p
      # L-BFGS-B can evaluate a boundary a few floating-point units below zero.
      if (any(z[c(nv, np)] < -1e-12)) stop("Optimizer evaluated a negative variance.")
      z[c(nv, np)] <- pmax(0, z[c(nv, np)])
      z
    }
    previous <- value <- NULL
    evaluate <- function(p) {
      if (!identical(previous, p)) { value <<- mfrm_testlet_evaluate(input, expand(p), rule); previous <<- p }
      value
    }
    fit <- stats::optim(start[free], function(p) -evaluate(p)$loglik,
      function(p) -evaluate(p)$gradient[free], method = "L-BFGS-B",
      lower = c(rep(-20, nv - 1), 0, 0)[free],
      upper = c(rep(20, nv - 1), variance_max, person_variance_max)[free],
      control = list(maxit = maxit, factr = 10, pgtol = 1e-7))
    fit$parameters <- expand(fit$par); fit$fixed_variance <- fixed_v
    fit$fixed_person_variance <- fixed_p
    fit
  }
  local_runs <- function(fixed_p) {
    if (!is.null(fixed)) return(list(known = optimize(start, fixed, fixed_p)))
    zero <- optimize(start, 0, fixed_p)
    positive <- larger <- zero$parameters
    positive[nv] <- min(.5, variance_max / 2); larger[nv] <- min(2, variance_max * .8)
    if (is.null(fixed_p)) {
      positive[np] <- min(1, person_variance_max / 2)
      larger[np] <- min(4, person_variance_max * .8)
    }
    list(zero = zero, positive = optimize(positive, NULL, fixed_p),
      larger = optimize(larger, NULL, fixed_p))
  }
  runs <- if (!is.null(fixed_person)) local_runs(fixed_person) else
    c(setNames(local_runs(0), paste0("person_zero_", if (is.null(fixed)) c("zero", "positive", "larger") else "known")),
      local_runs(NULL))
  if (is.null(fixed_person) && !is.null(fixed)) {
    larger <- start; larger[np] <- min(4, person_variance_max * .8)
    runs$larger_person <- optimize(larger, fixed, NULL)
  }
  best <- runs[[mfrm_testlet_select_run(runs)]]; par <- best$parameters
  low <- mfrm_testlet_evaluate(input, par, rule); higher <- mfrm_testlet_evaluate(input, par, high)
  boundary <- is.null(fixed) && par[nv] == 0
  person_boundary <- is.null(fixed_person) && par[np] == 0
  free <- c(seq_len(nv - 1L), if (is.null(fixed) && !boundary) nv,
    if (is.null(fixed_person) && !person_boundary) np)
  score <- higher$gradient[free]
  projected <- max(abs(score), if (boundary) max(0, higher$gradient[nv]) else 0,
    if (person_boundary) max(0, higher$gradient[np]) else 0)
  bound <- any(abs(par[seq_len(nv - 1L)]) >= 20 - 1e-6) ||
    (is.null(fixed) && par[nv] >= variance_max - 1e-6) ||
    (is.null(fixed_person) && par[np] >= person_variance_max - 1e-6)
  gradient_at <- if (is.null(fixed_person)) seq_len(np) else seq_len(nv)
  differences <- c(LogLikDifference = abs(higher$loglik - low$loglik),
    GradientDifference = max(abs(higher$gradient[gradient_at] - low$gradient[gradient_at])),
    MomentDifference = max(abs(higher$moments - low$moments)))
  ready <- best$convergence == 0L && !bound && projected < 1e-5 &&
    all(is.finite(differences)) && all(differences < c(1e-6, 1e-5, 1e-6))
  covariance <- matrix(NA_real_, np, np); H <- matrix(NA_real_, length(free), length(free))
  if (ready) for (j in seq_along(free)) {
    h <- 1e-4 * (1 + abs(par[free[j]])); if (free[j] >= nv) h <- min(h, par[free[j]] / 4)
    step <- numeric(np); step[free[j]] <- h
    H[, j] <- -(mfrm_testlet_evaluate(input, par + step, high)$gradient[free] -
      mfrm_testlet_evaluate(input, par - step, high)$gradient[free]) / (2 * h)
  }
  H <- (H + t(H)) / 2
  eigenvalues <- if (all(is.finite(H))) eigen(H, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  information <- all(is.finite(eigenvalues)) && min(eigenvalues) > 1e-7 && min(eigenvalues) / max(eigenvalues) > 1e-9
  if (information) covariance[free, free] <- solve(H)
  checks <- c(as.list(differences), list(MaxProjectedGradient = projected,
    VarianceScore = higher$gradient[nv], PersonVarianceScore = higher$gradient[np], Convergence = best$convergence,
    SearchBoundary = bound, EstimatedVarianceBoundary = boundary,
    EstimatedPersonVarianceBoundary = person_boundary,
    NumericalReady = ready, InformationPositive = information,
    MinInformationEigenvalue = min(eigenvalues)))
  list(parameters = par, covariance = covariance, loglik = higher$loglik, moments = higher$moments,
    person_loglik = higher$person_loglik, checks = checks, runs = runs)
}
