#' Describe fitted response probabilities and residuals
#'
#' Compute response means, variances and descriptive Infit/Outfit summaries.
#' Ordinary and extended RSMs and two-family GPCMs integrate uncertain latent effects; corrected
#' GPCM JML conditions on saved calibration and reprofiled Person estimates.
#' @param fit A numerically ready native RSM MML [fit_mfrm()],
#'   [fit_mfrm_testlet()] or [fit_mfrm_random_rater()] result, a converged
#'   two-family GPCM fixed-grid EM or adaptive direct MML fit, or a shared-owner
#'   GPCM fit with an explicit `jml_correction_order` and an available point solution.
#' @param rows Distinct original input row numbers to return, or `NULL` for
#'   all assigned rows. This selects outputs only: every observed source rating
#'   still informs the joint posterior or saved JML Person profile. Groups
#'   summarize the selected rows; selection does not refit the model.
#' @param group_by Fitted identifier columns to summarize separately. `NULL`
#'   uses all Person, testlet/rater and fixed-facet columns. No grouping by score.
#'   Corrected JML lists the fixed facets first, then Person, so the default
#'   plot shows the first declared facet.
#' @param quad_points Normal quadrature order, 7 to 121 for testlets or 7 to
#'   241 for ordinary RSMs, two-family GPCMs or shared raters. Defaults to the saved fit's order.
#'   Results are checked at `2 * quad_points + 1`;
#'   each category probability must change by no more than `1e-7`.
#'   Omit this argument for corrected JML, which performs no integration.
#' @details Calibration is held fixed. For a posterior replicate rating sharing the
#'   original row's latent effects, the target is
#'   \deqn{p_{ik}=E[P(Y_i^{rep}=k\mid b,\widehat\psi)\mid
#'     Y_{obs},\widehat\psi],}
#'   where b includes ability, and for extensions the local effect or jointly
#'   uncertain shared-rater vector. This uses the same observed data for
#'   estimation and checking; it is not held-out prediction. Define
#'   \deqn{\mu_i=\sum_k k p_{ik},\qquad
#'     V_i=\sum_k (k-\mu_i)^2p_{ik}.}
#'   For these posterior calculations, the mixture variance includes conditional rating variance and
#'   variance of conditional means. Substituting an ability/rater estimate or
#'   averaging conditional variances alone gives a different quantity.
#'
#'   For selected observed unit-weight rows in a group, descriptive Outfit is
#'   the mean of \eqn{(y_i-\mu_i)^2/V_i}; descriptive Infit is
#'   \eqn{\sum_i(y_i-\mu_i)^2/\sum_i V_i}. These summaries have no established
#'   expectation of one. There are no default cutoffs, ZSTD, p-values, automatic
#'   flags or rater-quality classifications. Ordinary plug-in Infit/Outfit
#'   values do not have the same probability definition and must not be
#'   compared directly. No response covariance or calibration uncertainty is
#'   supplied, so these results do not calibrate a formal goodness-of-fit test.
#'
#'   Ordinary RSMs integrate normal ability by quadrature, retaining every
#'   observed rating of the relevant Person. Supported ordinary fits have
#'   unit weights, additive unanchored severity facets, no interactions or
#'   shrinkage, Person as the noncentered facet, consecutive integer source
#'   categories, and either known N(0,1) or `population_formula = ~1`.
#'   The estimated population mean and variance are both retained. This
#'   separate route does not alter the plug-in indices from [diagnose_mfrm()].
#'   Old fits that omitted scores without retaining original row positions
#'   need to be refitted from the complete assigned-score roster.
#'
#'   Testlets use nested normal quadrature. Shared raters use category-specific
#'   Laplace integrals of the complete likelihood times the replicate category
#'   probability, after integrating abilities by quadrature. Each numerator
#'   reoptimizes the joint latent rater mode; calibration is never refitted.
#'   The positive numerator integrals are normalized across categories.
#'   `NormalizationError` records the unnormalized probability sum minus one.
#'   This defect and agreement between quadrature orders are numerical
#'   diagnostics, not bounds on the Laplace approximation error. In particular,
#'   small defects do not certify accuracy with few raters or sparse ratings.
#'   Shared-rater work grows with requested rows and categories; use `rows`
#'   to inspect a subset while retaining the complete conditioning data.
#'
#'   For these posterior routes, missing scores remain `missing_score`; they
#'   are not imputed. Nonfinite
#'   or zero predictive variances and unresolved integration remain unavailable
#'   with a reason. A group containing any unavailable observed row has no
#'   summary, rather than silently dropping that row. Entirely missing groups
#'   have no summary. A fitted zero latent variance defines the corresponding
#'   submodel's predictions; it does not establish absence of heterogeneity.
#'
#' @section Two-family GPCM posterior residuals:
#' Use this route to describe how observed scores differ from the fitted
#' model's posterior predictions. For example, summarize by a Judge column
#' with `mfrm_response_diagnostics(fit, group_by = "Judge")`. Both slope
#' families multiply the complete adjacent-category predictor, as in fitting.
#' Ability is integrated under the fitted fixed N(0,1) population and all of
#' the Person's observed ratings; probabilities are not evaluated at the EAP.
#' Only fitted observed rows are returned, preserving the original identifiers
#' and incomplete assignment pattern. No unassigned cell is filled.
#'
#' The saved solution and point calibration must agree with their fitting
#' specification. A slope covariance or slope interval is not required.
#' Fixed-grid EM and adaptive direct MML retain their respective integration
#' methods. Adaptive grids use the complete Person posterior's mode and curvature,
#' preserving the original N(0,1) prior through density-ratio and Jacobian weights.
#' Source likelihood and engine-specific gradient checks are recomputed before
#' prediction. Stored convergence flags alone do not establish a usable source.
#' This has the same scope as the two-family fit: two fixed slope facets,
#' unit weights, no anchors or covariates, and observed zero-based categories.
#' `quad_points` checks response integration at the retained parameters; it
#' does not refit the calibration. If integration is unresolved, increase it
#' explicitly. The returned probabilities use `2 * quad_points + 1` points;
#' the lower order supplies the row-wise check. Selecting fewer output rows
#' does not shorten the posterior's conditioning record or change its grid.
#' Settings retain the integration method through plots, reports and exports.
#' Separately use [mml_quadrature_sensitivity()] to examine whether
#' the fitted calibration changes across quadrature grids.
#'
#' Inspect `rows` for unavailable predictions before reading `measures`.
#' A larger descriptive Infit or Outfit is not a calibrated misfit decision;
#' there is no expectation-one reference. These same-data summaries neither
#' establish rater quality nor test unidimensionality or local independence.
#' [diagnose_mfrm()], Q3/PCA, Wright/Pathway maps and model ranking retain
#' their separate support restrictions. Attach the saved result with
#' `mfrm_results(fit, include = c("fit", "plots"), compute = "never",
#' response_diagnostics = result)` for plots, reports and exports. Collection
#' and replay do not refit or reintegrate the saved probabilities.
#'
#' @section Corrected JML conditional residuals:
#' For an explicit profile-score-adjusted GPCM fit, probabilities are evaluated
#' at the saved corrected facet/step/slope estimates and the reprofiled Person
#' abilities: \eqn{p_{ik}=P(Y_i=k\mid\widehat\theta_i,\widehat\beta_{adj})}.
#' These are conditional plug-in probabilities, not posterior predictions or
#' held-out predictions. There is no normal ability prior or quadrature.
#' Calibration and Person estimation uncertainty are excluded. Residual bias
#' may remain. Unavailable RootSE does not prevent these point calculations.
#' The response formula is the same shared-owner GPCM used for fitting.
#'
#' Only observed fitted rows are available; no unassigned cells are filled or
#' missing scores imputed. Original row numbers, repeated-event multiplicity,
#' category labels and facet names are retained. `rows` selects displayed and
#' summarized events without changing their Person profiles. Save the result
#' and attach it through `mfrm_results(fit, response_diagnostics = result)`;
#' the saved calibration, slopes, Person profiles and source rows must match.
#'
#' An extreme Person profile has a limiting probability of one at the minimum
#' or maximum score. Its expected score, zero conditional variance and raw
#' residual remain available. `StandardizedResidual` is missing with
#' `Status = "zero_variance"`: the ratio at the saved probabilities is
#' undefined. No automatic convention replacing `0/0` with zero is applied.
#' `ProbabilityAvailable` distinguishes this case from a failed probability
#' calculation. If all selected raw residuals and variances are available and
#' their total variance is positive, Infit remains available, including those
#' zero-variance rows. Outfit requires every selected standardized residual;
#' otherwise it is missing, and the group is `partially_available` with a
#' reason. The grouped `Available` column counts rows with finite standardized
#' residuals; Infit can use all `Observed` rows despite that smaller count.
#' No row is silently discarded. An all-zero-variance group has neither
#' index. Explicitly selecting other rows changes the summarized set and must
#' be reported as such.
#'
#' These same-data residual summaries have no calibrated expectation-one
#' reference, uncertainty intervals, p-values, automatic flags or rater-quality
#' classifications. They do not establish correction of structural bias or
#' replace [diagnose_mfrm()], which is unavailable for this estimator.
#' `plot()` offers paired and scatter views, including [as_ggplot()] conversion.
#' A paired view retains available Infit when Outfit is missing; a scatter view
#' needs both. Reports and exports preserve the probability definition and
#' unavailable outcomes. Corrected-JML model comparison remains unavailable;
#' portable scoring has its own calibration checks and conditional EAP target.
#' @return An `mfrm_response_diagnostics` object with `rows`, category
#'   `probabilities`, grouped `measures`, settings and exact source metadata.
#'   Save with `saveRDS()`; attach via
#'   `mfrm_results(fit, response_diagnostics = result, compute = "never")`.
#'   The resulting `mfrm_results` object contains `tables$response_overview`,
#'   distinguishing available standardized residuals, unresolved calculations,
#'   missing scores, zero-variance rows and source rows not selected. Partial
#'   results require review; availability does not establish model adequacy.
#'   No probabilities or residuals are recomputed for this table.
#'   Older saved `mfrm_results` objects without this table recover the overview
#'   and status from retained diagnostic rows when summarized, reported, exported
#'   or opened in a supported results viewer. The original object is not modified;
#'   newly exported results include the recovered overview and status.
#'   Supply two saved outputs to [compare_mfrm()] for aligned ordinary versus
#'   extended RSM comparisons, using the same selected events and group_by.
#'   Corrected-JML and two-family GPCM comparisons are not supported by that route.
#'   Collecting, plotting and exporting saved results never recomputes integrals.
#' @references Tierney, L., & Kadane, J. B. (1986). Accurate approximations
#'   for posterior moments and marginal densities. Journal of the American
#'   Statistical Association, 81(393), 82--86.
#'   \doi{10.1080/01621459.1986.10478240}. This supports the ratio-of-integrals
#'   approximation, not universal accuracy for this rating design.
#' @seealso [plot.mfrm_response_diagnostics()], [mfrm_results()],
#'   [score_mfrm_random_rater()], [predict.mfrm_testlet()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("mfrm_response_diagnostics", run.dontrun = TRUE)
#' \dontrun{
#' example <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr"))
#' fit <- example$testlet$fit
#' residuals <- example$testlet$diagnostics
#' # To recompute the posterior predictive residuals:
#' # residuals <- mfrm_response_diagnostics(fit, group_by = "Rater")
#' residuals$measures
#' plot(residuals)
#' mfrm_results(fit, response_diagnostics = residuals, compute = "never")
#' }
#' @export
#' @order 1
mfrm_response_diagnostics <- function(fit, rows = NULL, group_by = NULL,
    quad_points = fit$settings$quad_points %||% fit$config$estimation_control$quad_points) {
  if (mfrm_has_jml_adjustment(fit)) {
    if (!missing(quad_points)) stop("Corrected-JML response diagnostics use saved conditional point estimates; omit quad_points.",call.=FALSE)
    return(mfrm_jml_response_diagnostics(fit,rows,group_by))
  }
  product <- mfrm_has_product_slopes(fit)
  ordinary <- inherits(fit, "mfrm_fit")
  if (!ordinary && (!mfrm_extended_fit(fit) || !isTRUE(fit$checks$NumericalReady) ||
      !isTRUE(fit$checks$InformationPositive))) stop("Supply a numerically ready ordinary RSM MML, testlet or shared-rater fit; review its numerical and information checks.", call. = FALSE)
  testlet <- inherits(fit, "mfrm_testlet"); maximum <- if (testlet) 121 else 241
  if (!is.numeric(quad_points) || is.complex(quad_points) || length(quad_points) != 1L ||
      !is.finite(quad_points) || quad_points != floor(quad_points) || quad_points < 7 || quad_points > maximum) stop("Invalid quad_points for this model.", call. = FALSE)
  input <- if (product) mfrm_gmfrm_response_input(fit) else
    if (ordinary) mfrm_ordinary_response_input(fit) else fit$input
  columns <- input$columns
  if (is.null(input$assigned_data) && length(input$omitted_rows)) stop("The saved fit lacks its original missing-score roster; recover the complete fitted object.", call. = FALSE)
  source_data <- input$assigned_data %||% input$data
  n <- nrow(source_data); rows <- rows %||% seq_len(n)
  if (!is.numeric(rows) || is.complex(rows) || !is.null(dim(rows)) || !length(rows) || anyNA(rows) ||
      any(rows != floor(rows)) || any(rows < 1 | rows > n) || anyDuplicated(rows)) stop("`rows` must contain distinct original input row numbers.", call. = FALSE)
  identifiers <- setdiff(names(source_data), columns$score)
  group_by <- group_by %||% identifiers
  if (!is.character(group_by) || !length(group_by) || anyNA(group_by) ||
      anyDuplicated(group_by) || !all(group_by %in% identifiers)) stop("Choose distinct fitted identifier columns for group_by.", call. = FALSE)
  table <- data.frame(InputRow = rows, Score = source_data[[columns$score]][rows])
  observed <- which(!is.na(source_data[[columns$score]]))
  index <- match(rows, observed); valid <- which(!is.na(index))
  nc <- length(input$score_levels)
  probabilities <- matrix(NA_real_, length(rows), nc, dimnames = list(as.character(rows), as.character(input$score_levels)))
  table$ExpectedScore <- table$PredictiveVariance <- table$SquaredResidual <-
    table$StandardizedResidual <- table$IntegrationDifference <- table$NormalizationError <- NA_real_
  table$Status <- ifelse(is.na(index), "missing_score", "unavailable")
  table$Reason <- ifelse(is.na(index), "Assigned score is missing; no residual was computed.", "")
  # One testlet calculation per Person; shared raters require row-specific
  # numerator integrals, with all Persons retained in every integral.
  batches <- if (testlet || ordinary) split(valid, source_data[[columns$person]][rows[valid]]) else as.list(valid)
  normalizers <- list()
  for (batch in batches) {
    answer <- tryCatch(withCallingHandlers({
      calculate <- function(q) {
        if (product) return(mfrm_gmfrm_response_probabilities(input, index[batch], q))
        if (ordinary) return(mfrm_ordinary_response_probabilities(input, index[batch], q))
        if (testlet) return(mfrm_testlet_response_probabilities(input, fit$parameters, index[batch], q))
        value <- mfrm_random_rater_response_probabilities(input, fit$calibration, index[batch], q,
          logden = normalizers[[as.character(q)]])
        normalizers[[as.character(q)]] <<- value$logden
        value
      }
      low <- calculate(quad_points); high <- calculate(2L * quad_points + 1L)
      delta <- apply(abs(high$probabilities - low$probabilities), 1L, max)
      list(high = high, delta = delta)
    }, warning = function(w) stop(conditionMessage(w), call. = FALSE)),
    error = function(e) { table$Reason[batch] <<- conditionMessage(e); NULL })
    if (is.null(answer)) next
    table$IntegrationDifference[batch] <- answer$delta
    table$NormalizationError[batch] <- answer$high$normalization_error
    for (j in seq_along(batch)) {
      at <- batch[j]; p <- answer$high$probabilities[j, ]
      if (any(!is.finite(p)) || any(p < 0) || abs(sum(p) - 1) > 1e-10 ||
          !is.finite(answer$high$normalization_error[j]) ||
          !is.finite(answer$delta[j]) || answer$delta[j] > 1e-7) {
        table$Reason[at] <- "Response integration is unresolved; increase quad_points and inspect the rating design."
        next
      }
      mu <- sum(input$score_levels * p); variance <- sum((input$score_levels - mu)^2 * p)
      if (!is.finite(variance) || variance <= 0) {
        table$Reason[at] <- "Predictive variance is nonfinite or zero."
        next
      }
      probabilities[at, ] <- p
      table$ExpectedScore[at] <- mu; table$PredictiveVariance[at] <- variance
      residual <- table$Score[at] - mu
      table$SquaredResidual[at] <- residual^2
      table$StandardizedResidual[at] <- residual / sqrt(variance)
      table$Status[at] <- "available_conditional"
    }
  }
  measures <- mfrm_response_group_measures(table,source_data,rows,group_by)
  structure(list(rows = table, probabilities = probabilities, measures = measures,
    source = mfrm_response_source(fit), source_data = source_data, source_observed = input$data,
    settings = list(model = if (product) "Two-family GPCM" else if (ordinary) "Ordinary fixed-facet RSM" else if (testlet) "Testlet RSM" else "Shared-rater RSM",
      quad_points = quad_points, check_points = 2L * quad_points + 1L,
      probability_tolerance = 1e-7,
      group_by = group_by, calibration_uncertainty = FALSE,
      target = "Same-data posterior predictive replicate; original latent effects shared; calibration fixed",
      integration = if (product) paste0(if (mfrmr_adaptive_integration(fit$config))
        "Adaptive" else "Fixed-grid", " quadrature under N(0,1); both slope families, locations and steps held fixed") else if (ordinary) "Normal ability quadrature; fixed facets and calibration held fixed" else if (testlet) "Nested normal quadrature" else "Normalized category-specific joint-rater Laplace integrals with ability quadrature",
      limitation = "Descriptive only; no expectation-one reference, cutoffs, ZSTD or p-values. Numerical agreement does not certify approximation accuracy.",
      selection = "Selected original rows summarized; all observed source ratings condition every prediction")),
    class = "mfrm_response_diagnostics")
}

mfrm_response_group_measures <- function(table, source_data, rows, group_by) {
  do.call(rbind, lapply(group_by, function(facet) do.call(rbind,
    lapply(unique(as.character(source_data[[facet]][rows])), function(id) {
      z <- table[as.character(source_data[[facet]][rows]) == id, , drop = FALSE]
      available <- z$Status == "available_conditional"; obs <- !is.na(z$Score)
      complete <- all(available[obs]) && any(obs)
      infit_available <- any(obs) && all(z$Status[obs] %in% c("available_conditional","zero_variance")) &&
        all(is.finite(z$SquaredResidual[obs])) && all(is.finite(z$PredictiveVariance[obs])) &&
        sum(z$PredictiveVariance[obs])>0
      data.frame(Facet = facet, Level = id, Selected = nrow(z), Observed = sum(obs),
        Missing = sum(!obs), Available = sum(available),
        Infit = if (infit_available) sum(z$SquaredResidual[obs]) / sum(z$PredictiveVariance[obs]) else NA_real_,
        Outfit = if (complete) mean(z$StandardizedResidual[obs]^2) else NA_real_,
        Status = if (complete) "descriptive_only" else if (infit_available) "partially_available" else "unavailable",
        Reason = if (complete) "" else if (infit_available)
          "Infit retains all selected rows; Outfit is undefined because at least one conditional variance is zero." else if (!any(obs)) "No observed scores in selected rows." else
          "At least one selected observed row is unavailable.", row.names = NULL)
    }))))
}

mfrm_testlet_response_probabilities <- function(input, par, rows, order) {
  rule <- gauss_hermite_normal(order)
  theta <- sqrt(mfrm_testlet_person_variance(input, par)) * rule$nodes
  nb <- ncol(input$X); ns <- length(input$score_levels) - 1L
  variance <- par[nb + ns + 1L]; steps <- par[nb + seq_len(ns)]
  gamma <- if (variance == 0) 0 else sqrt(variance) * rule$nodes
  gw <- if (variance == 0) 1 else rule$weights
  person <- match(input$data[[input$columns$person]][rows[1L]], input$persons)
  kernel <- mfrm_testlet_kernel(input, par, theta, rule, persons = person)
  lw <- (kernel$loglik[, 1L] - max(kernel$loglik[, 1L])) + log(rule$weights)
  pw <- exp(lw - max(lw)); pw <- pw / sum(pw)
  result <- matrix(NA_real_, length(rows), ns + 1L)
  for (block in input$groups[[person]]) {
    selected <- which(rows %in% block)
    if (!length(selected)) next
    cache <- lapply(block, function(i) {
      eta <- as.vector(outer(theta, gamma, "+")) - sum(input$X[i, ] * par[seq_len(nb)])
      logw <- outer(eta, 0:ns) - matrix(c(0, cumsum(steps)), length(eta), ns + 1L, byrow = TRUE)
      logw <- logw - row_max_fast(logw)
      logw - mfrm_testlet_logsum_rows(logw)
    })
    loglik <- Reduce(`+`, lapply(seq_along(block), function(j)
      matrix(cache[[j]][, input$y[block[j]] + 1L], length(theta), length(gamma))))
    local <- sweep(loglik - row_max_fast(loglik), 2L, log(gw), "+")
    local <- exp(local - mfrm_testlet_logsum_rows(local))
    joint <- as.vector(local * pw)
    for (j in selected) result[j, ] <- colSums(exp(cache[[match(rows[j], block)]]) * joint)
  }
  list(probabilities = result, normalization_error = rowSums(result) - 1)
}

mfrm_random_rater_response_probabilities <- function(input, calibration, rows, order, logden = NULL) {
  row <- rows[1L]; steps <- calibration$steps; par <- c(calibration$beta, steps)
  make <- function(replicate = NULL) mfrm_random_rater_objective(input, order,
    fixed_sd = calibration$rater_sd, fixed_person_sd = calibration$person_sd %||% 1,
    replicate_row = replicate)
  if (is.null(logden)) {
    baseline <- make(); logden <- -as.numeric(baseline$fn(par))
  }
  numerator <- make(row)
  lognum <- vapply(0:length(steps), function(k)
    -as.numeric(numerator$fn(c(par, k, as.numeric(seq_along(steps) <= k)))), numeric(1))
  if (!is.finite(logden) || any(!is.finite(lognum))) stop("Joint-rater response integration is nonfinite.")
  maximum <- max(lognum)
  shifted <- lognum - maximum
  normalizer <- log(sum(exp(shifted)))
  list(probabilities = matrix(exp(shifted - normalizer), 1L),
    normalization_error = expm1((maximum - logden) + normalizer), logden = logden)
}

#' @rdname mfrm_response_diagnostics
#' @param object A saved response-diagnostics result.
#' @param ... For `summary.mfrm_response_diagnostics()`: Unused.
#'
#' For `plot()`: Unused; ordinary-model threshold arguments are not accepted.
#' @export
summary.mfrm_response_diagnostics <- function(object, ...) {
  list(measures = object$measures, status_counts = table(object$rows$Status), settings = object$settings)
}

#' @rdname mfrm_response_diagnostics
#' @export
print.mfrm_response_diagnostics <- function(x, ...) {
  cat(if (identical(x$settings$probability_method,"corrected_jml_plugin"))
    "Corrected JML conditional residual summaries (descriptive only)\n" else
    "Posterior predictive residual summaries (descriptive only)\n")
  cat(x$settings$target,"\n")
  print(x$measures, row.names = FALSE)
  cat("Same data used for fitting and checking; no reference cutoffs or calibrated fit test.\n")
  invisible(x)
}
