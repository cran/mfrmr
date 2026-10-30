#' Compare ordered responses using cumulative logistic models
#'
#' Compare two groups at an explicitly supplied numeric matching variable,
#' using three nested proportional-odds models. This independent-Person
#' sensitivity analysis does not estimate ability or fit an MFRM.
#'
#' @param data Long-format data with at most one row per Person and target.
#' @param person Name of the independent Person identifier column.
#' @param target Name of the item or other target identifier column.
#' @param score Name of the ordered response column. Numeric, character or
#'   factor scores are accepted, with order declared through `categories`.
#' @param group Name of the group column. Each Person belongs to one group.
#' @param matching Name of one finite numeric matching column. Values may
#'   differ across targets. Missing scores or matching values are excluded
#'   together from all models for the affected target.
#' @param groups Two distinct labels in reference, focal order.
#' @param categories Distinct numeric or character category labels, in increasing
#'   response order. Factor labels are treated as character. Category spacing
#'   has no effect. Every declared category must occur in each target's complete
#'   data; otherwise that target is unavailable. Categories are never collapsed.
#' @param p_adjust A [stats::p.adjust()] method, default `"holm"`. The DIF
#'   comparison family comprises all three tests for every observed target,
#'   including unavailable tests. Proportional-odds diagnostics form a separate
#'   family of one test per target.
#'
#' @details
#' Requires the optional package \pkg{ordinal}. Maximum likelihood is computed
#' by [ordinal::clm()] with a logit link and flexible thresholds. For each target,
#' matching is centered at its pooled midrange and divided by its half-range:
#' `z = (matching - Center) / Scale`. Group is coded 0 for reference and 1 for
#' focal. `M1` uses `z`, `M2` uses `z + g`, and `M3` uses `z * g`. All three use
#' identical complete responses. The model is
#' `logit Pr(Y <= k) = threshold[k] - linear_predictor`.
#'
#' `Uniform` compares M2 to M1 (1 df), `Nonuniform` compares M3 to M2 (1 df),
#' and `Joint` compares M3 to M1 (2 df). These are asymptotic likelihood-ratio
#' tests. A uniform comparison assumes no interaction; its coefficient must not
#' be interpreted as a constant group difference when an interaction matters.
#' No automatic selection between these models or DIF/bias flag is made.
#' `NagelkerkeR2` is `(1-exp(2*(LL0-LL)/N))/(1-exp(2*LL0/N))`, where LL0 is the
#' intercept-only multinomial log likelihood. `DeltaNagelkerkeR2` subtracts
#' the reduced model's value from the full model's value. This is not explained
#' score variance or McKelvey--Zavoina R-squared; no universal effect-size cutoff
#' is applied.
#'
#' `contrasts` gives M3's reference/focal upper-tail odds ratio at the lower
#' endpoint, midpoint and upper endpoint of the groups' overlapping matching
#' ranges. Its log is `-(beta_g + beta_z:g * z)`; a value above 1 favors higher
#' reference scores. This applies to each category cut under proportional odds.
#' These are point estimates, not simultaneous intervals. Fits use all complete
#' responses, including those outside that overlap. Range overlap alone does not
#' establish comparable distributions or adequate local support. Each group
#' must have varying matching values and their ranges must overlap with positive
#' width. The `matching` table retains the ranges and transformation.
#'
#' For three or more categories, `proportional_odds` compares M3 with a model
#' in which all three slopes may differ across cumulative category cuts
#' (`clm(y ~ 1, nominal = ~ z * g)`). This global asymptotic LR diagnostic has
#' `3 * (K - 2)` df. It concerns equality across category cuts, whereas the
#' nonuniform comparison concerns differences across matching values. They are
#' different questions. The diagnostic is not applicable to binary responses.
#' An unavailable or nonsignificant diagnostic does not certify proportional
#' odds. Its p-value does not determine availability or change the primary comparisons.
#'
#' `models` retains convergence codes, warnings, gradient, Hessian and Newton
#' step diagnostics. A fit is unavailable if the backend reports nonconvergence
#' or nonunique parameters, any coefficient is aliased/nonfinite, the Hessian's
#' minimum eigenvalue is at most `sqrt(.Machine$double.eps) * max(1, max_eigen)`,
#' the maximum absolute gradient exceeds `1e-6`, or the maximum absolute Newton
#' step exceeds `1e-4` on the scaled parameterization. Backend warnings also make
#' the fit unavailable. These conservative numerical checks detect common
#' separation failures; they are not a general proof of finite ML existence.
#' Likelihood decreases beyond `64 * .Machine$double.eps * max(1, abs(LL))`
#' invalidate a comparison; smaller rounding errors are set to zero. No penalized
#' estimator, collapsed categories or alternative inference is substituted.
#'
#' Inference requires independent Persons within each target, a suitable linear
#' matching relationship, finite identifiable parameters and adequate data for
#' chi-squared approximations. Sparse categories and small samples can invalidate
#' those approximations even when numerical checks pass. Duplicate Person-target
#' rows, including missing-score duplicates, are rejected before exclusions.
#' Renaming repeated ratings as Persons does not establish independence.
#' Dependence across targets is compatible with Holm adjustment of valid
#' individual p-values, but no dependence within targets is handled here.
#'
#' Interpreting conditional group association as DIF additionally requires
#' defensible matching, invariant anchors and comparable rater/task exposure.
#' Use a justified common rating plan; this interface does not adjust for plan
#' differences, matching error, clustered ratings or informative missingness.
#' Total and rest scores need not define the same comparison. RSM/PCM's adjacent
#' category logits do not generally imply cumulative proportional odds. Therefore
#' agreement with [analyze_dff()] or [mantel_haenszel_table()] does not by itself
#' establish robustness, and disagreement does not by itself establish DIF.
#'
#' @return A plain list containing `table` (three comparisons per target),
#'   `models` (three PO fits and, for ordinal responses, a diagnostic fit),
#'   `coefficients` (available models only), `contrasts`, `proportional_odds`,
#'   `matching`, `cases`, `settings` and `notes`. `cases$Included` denotes complete
#'   rows supplied to eligible models, not successful estimation or valid DIF
#'   inference. Precluded targets retain their counts and reasons. Save the full
#'   result with [base::saveRDS()] and individual tables with [utils::write.csv()].
#'   No backend model objects or executable environments are retained.
#'
#' @references
#' Rogers, H. J., & Swaminathan, H. (1993). A comparison of logistic regression
#' and Mantel-Haenszel procedures for detecting differential item functioning.
#' *Applied Psychological Measurement, 17*, 105--116.
#' \doi{10.1177/014662169301700201}.
#'
#' Zumbo, B. D. (1999). *A handbook on the theory and methods of differential
#' item functioning (DIF): Logistic regression modeling as a unitary framework
#' for binary and Likert-type (ordinal) item scores*. Department of National
#' Defense, Ottawa, pp. 23--32. Pseudo-R-squared definitions and thresholds in
#' that handbook are not interchangeable.
#' @seealso [analyze_dff()], [mantel_haenszel_table()], [mfrmr_linking_and_dff]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("ordinal_logistic_table", run.dontrun = TRUE)
#' \dontrun{
#' if (requireNamespace("ordinal", quietly = TRUE)) {
#'   d <- expand.grid(Match = c(-1, 0, 1), Group = c("R", "F"),
#'                    Score = 0:2, Replicate = 1:5)
#'   d$Person <- paste0("P", seq_len(nrow(d)))
#'   d$Item <- "A"
#'   ol <- ordinal_logistic_table(d, "Person", "Item", "Score", "Group",
#'                                "Match", c("R", "F"), 0:2)
#'   ol$table
#'   ol$proportional_odds
#' }
#' }
#' @export
ordinal_logistic_table <- function(data, person, target, score, group, matching,
                                   groups, categories, p_adjust = "holm") {
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  if (!is.data.frame(data) || !nrow(data) || anyNA(names(data)) ||
      anyDuplicated(names(data)) || any(!nzchar(names(data)))) stop(
    "`data` must be a nonempty data frame with unique nonblank column names.", call. = FALSE)
  data <- as.data.frame(data)
  for (nm in c("person", "target", "score", "group", "matching")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v)) stop(
      "`", nm, "` must name one column.", call. = FALSE)
  }
  selected <- c(person, target, score, group, matching)
  if (anyDuplicated(selected) || !all(selected %in% names(data))) stop(
    "Person, target, score, group and matching must name distinct existing columns.", call. = FALSE)
  if (!is.character(groups) || length(groups) != 2L || anyNA(groups) ||
      any(!nzchar(trimws(groups))) || anyDuplicated(groups)) stop(
    "`groups` must contain two distinct nonblank labels in reference, focal order.", call. = FALSE)
  for (nm in c(person, target, group)) {
    v <- data[[nm]]
    if (!is.null(dim(v)) || !(is.factor(v) || (!is.object(v) &&
        (is.character(v) || is.logical(v) || (is.numeric(v) && !is.complex(v))))) ||
        anyNA(v) || any(!nzchar(trimws(as.character(v)))) ||
        (is.numeric(v) && any(!is.finite(v))) ||
        length(unique(v)) != length(unique(as.character(v)))) stop(
      "Identifiers must have distinct finite, nonblank atomic labels: ", nm, ".", call. = FALSE)
  }
  x <- data[[matching]]
  if (!is.numeric(x) || is.complex(x) || is.object(x) || !is.null(dim(x)) ||
      any(!is.na(x) & !is.finite(x))) stop(
    "Matching must be finite numeric values or NA.", call. = FALSE)
  y <- data[[score]]
  if (is.factor(y)) y <- as.character(y)
  if (is.factor(categories)) categories <- as.character(categories)
  valid_labels <- function(v) !is.object(v) && is.null(dim(v)) &&
    (is.character(v) || (is.numeric(v) && !is.complex(v) && all(is.na(v) | is.finite(v))))
  if (!valid_labels(y) || !valid_labels(categories) || length(categories) < 2L ||
      anyNA(categories) || anyDuplicated(categories) ||
      is.numeric(y) != is.numeric(categories)) stop(
    "Scores and `categories` must have matching numeric or character types; declare at least two distinct nonmissing categories.",
    call. = FALSE)
  code <- match(y, categories)
  if (any(!is.na(y) & is.na(code))) stop("Scores outside declared `categories`.", call. = FALSE)
  if (anyDuplicated(data[c(person, target)])) stop(
    "Duplicate Person-target cells are not supported; repeated ratings are not independent Persons.", call. = FALSE)
  pid <- as.character(data[[person]]); tid <- as.character(data[[target]])
  g <- as.character(data[[group]])
  if (!all(g %in% groups)) stop("Data contain groups outside `groups`; subset explicitly.", call. = FALSE)
  if (anyDuplicated(pid[!duplicated(data[c(person, group)])])) stop(
    "Each Person must belong to exactly one group across all targets.", call. = FALSE)
  if (!requireNamespace("ordinal", quietly = TRUE)) stop(
    "Install the optional package `ordinal` to use `ordinal_logistic_table()`.", call. = FALSE)

  k <- length(categories)
  complete <- !is.na(code) & !is.na(x)
  cases <- data.frame(Row = seq_len(nrow(data)), Person = pid, Target = tid,
    Group = g, Score = y, Matching = x, Included = FALSE,
    Reason = ifelse(is.na(code) & is.na(x), "missing_score_and_matching",
      ifelse(is.na(x), "missing_matching", ifelse(is.na(code), "missing_score", "pending"))))
  targets <- unique(tid)
  comparisons <- models <- coefficients <- contrasts <- diagnostics <- matches <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    t <- targets[i]; rows <- which(tid == t); used <- rows[complete[rows]]
    counts <- tabulate(code[used], nbins = k)
    xr <- x[used[g[used] == groups[1]]]; xf <- x[used[g[used] == groups[2]]]
    rr <- if (length(xr)) range(xr) else c(NA_real_, NA_real_)
    rf <- if (length(xf)) range(xf) else c(NA_real_, NA_real_)
    lo <- max(rr[1], rf[1]); hi <- min(rr[2], rf[2])
    center <- scale <- NA_real_
    reason <- if (!length(used)) "no_complete_responses" else
      if (!length(xr) || !length(xf)) "missing_group" else
        if (any(counts == 0)) "unused_category" else
          if (rr[1] == rr[2] || rf[1] == rf[2]) "constant_matching_within_group" else
            if (lo >= hi) "no_matching_range_overlap" else "eligible"
    if (reason == "eligible") {
      center <- min(x[used]) / 2 + max(x[used]) / 2
      scale <- max(abs(x[used] - center))
      if (!is.finite(scale) || scale <= 0) reason <- "invalid_matching_scale"
    }
    matches[[i]] <- data.frame(Target = t, Reference = groups[1], Focal = groups[2],
      InputRows = length(rows), CompleteRows = length(used), NReference = length(xr), NFocal = length(xf),
      ReferenceMin = rr[1], ReferenceMax = rr[2], FocalMin = rf[1], FocalMax = rf[2],
      OverlapMin = lo, OverlapMax = hi, Center = center, Scale = scale, Status = reason)
    cases$Reason[used] <- if (reason == "eligible") "model_input" else reason
    cases$Included[used] <- reason == "eligible"
    fits <- setNames(vector("list", 3L + (k > 2L)), c("M1", "M2", "M3", if (k > 2L) "NonPO"))
    ll0 <- if (length(used)) sum(counts[counts > 0] * log(counts[counts > 0] / length(used))) else NA_real_
    forms <- list(y ~ z, y ~ z + g, y ~ z * g, y ~ 1)
    for (j in seq_along(fits)) {
      f <- mfrm_ordinal_fit(if (reason == "eligible") data.frame(
        y = ordered(code[used], levels = seq_len(k)),
        z = (x[used] - center) / scale, g = as.integer(g[used] == groups[2])) else NULL,
        forms[[j]], nominal = if (j == 4L) ~ z * g else NULL, unavailable = reason)
      f$row$Target <- t; f$row$Model <- names(fits)[j]
      f$row$N <- length(used)
      f$row$NagelkerkeR2 <- if (f$row$Status == "available")
        (-expm1(2 * (ll0 - f$row$LogLik) / length(used))) / (-expm1(2 * ll0 / length(used))) else NA_real_
      fits[[j]] <- f
    }
    models[[i]] <- do.call(rbind, lapply(fits, `[[`, "row"))
    coefficients[[i]] <- do.call(rbind, lapply(names(fits), function(name) {
      f <- fits[[name]]
      if (f$row$Status != "available") return(NULL)
      data.frame(Target = t, Model = name, Parameter = names(f$coefficients), Estimate = unname(f$coefficients))
    }))
    pairs <- list(Uniform = c("M1", "M2"), Nonuniform = c("M2", "M3"), Joint = c("M1", "M3"))
    comparisons[[i]] <- do.call(rbind, lapply(names(pairs), function(name) {
      p <- pairs[[name]]
      out <- mfrm_ordinal_lr(fits[[p[1]]]$row, fits[[p[2]]]$row)
      data.frame(Target = t, Comparison = name, Reduced = p[1], Full = p[2], out,
        DeltaNagelkerkeR2 = if (out$Status == "available")
          fits[[p[2]]]$row$NagelkerkeR2 - fits[[p[1]]]$row$NagelkerkeR2 else NA_real_)
    }))
    diagnostics[[i]] <- data.frame(Target = t,
      if (k > 2L) mfrm_ordinal_lr(fits$M3$row, fits$NonPO$row) else
        data.frame(ChiSquared = NA_real_, DF = NA_real_, PValue = NA_real_, Status = "not_applicable_binary"))
    points <- if (reason == "eligible") c(lo, lo / 2 + hi / 2, hi) else rep(NA_real_, 3)
    logor <- rep(NA_real_, 3)
    if (fits$M3$row$Status == "available") {
      b <- fits$M3$coefficients
      logor <- unname(-(b["g"] + b["z:g"] * (points - center) / scale))
    }
    contrasts[[i]] <- data.frame(Target = t, Reference = groups[1], Focal = groups[2],
      Position = c("OverlapMin", "OverlapMid", "OverlapMax"), Matching = points,
      LogOddsRatio = logor, OddsRatio = exp(logor), Status = fits$M3$row$Status)
  }
  out <- lapply(list(table = comparisons, models = models, coefficients = coefficients,
    contrasts = contrasts, proportional_odds = diagnostics, matching = matches), function(x) {
      ans <- do.call(rbind, x); rownames(ans) <- NULL; ans
    })
  if (is.null(out$coefficients)) out$coefficients <- data.frame(Target = character(),
    Model = character(), Parameter = character(), Estimate = double())
  out$table$PAdjusted <- stats::p.adjust(out$table$PValue, method = p_adjust, n = nrow(out$table))
  out$proportional_odds$PAdjusted <- stats::p.adjust(out$proportional_odds$PValue,
    method = p_adjust, n = length(targets))
  out$cases <- cases
  out$settings <- list(person = person, target = target, score = score, group = group,
    matching = matching, groups = groups, categories = categories, p_adjust = p_adjust,
    adjustment_family_size = nrow(out$table), diagnostic_family_size = length(targets),
    link = "logit", group_coding = c(reference = 0L, focal = 1L),
    backend = "ordinal::clm", backend_version = as.character(utils::packageVersion("ordinal")),
    gradient_tolerance = 1e-6, newton_step_tolerance = 1e-4,
    hessian_relative_tolerance = sqrt(.Machine$double.eps))
  out$notes <- c(
    "Likelihood-ratio p-values are asymptotic; numerical availability does not certify DIF inference.",
    "Contrasts are reference/focal upper-tail odds ratios from M3 within overlapping matching ranges.",
    "Proportional-odds diagnostics and matching-by-group interactions address different assumptions.",
    "Matching validity, common exposure, independence and sparse-sample inference require substantive justification.",
    "Repeated Person-target ratings and estimation of matching/calibration uncertainty are not supported.")
  out
}

# Keep backend failures and numerical qualification separate from test results.
mfrm_ordinal_fit <- function(data, formula, nominal, unavailable) {
  row <- data.frame(LogLik = NA_real_, Parameters = NA_real_, Convergence = NA_integer_,
    MaxGradient = NA_real_, MinEigenvalue = NA_real_, ConditionNumber = NA_real_,
    MaxNewtonStep = NA_real_, Status = unavailable, Message = "")
  if (is.null(data)) return(list(row = row, coefficients = NULL))
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(ordinal::clm(formula, data = data, nominal = nominal,
    link = "logit", threshold = "flexible", control = ordinal::clm.control(convergence = "silent")),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL })
  row$Message <- paste(unique(warnings), collapse = "; ")
  row$Status <- "computation_failed"
  if (is.null(fit)) return(list(row = row, coefficients = NULL))
  row$LogLik <- fit$logLik; row$Parameters <- fit$edf
  row$Convergence <- fit$convergence$code; row$MaxGradient <- fit$maxGradient
  row$Message <- paste(c(warnings, fit$convergence$messages, fit$convergence$alg.message), collapse = "; ")
  finite <- all(is.finite(c(fit$coefficients, fit$Hessian, fit$gradient, fit$logLik)))
  eigenvalues <- if (finite) eigen(fit$Hessian, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  row$MinEigenvalue <- min(eigenvalues)
  row$ConditionNumber <- max(eigenvalues) / min(eigenvalues)
  row$MaxNewtonStep <- if (finite) tryCatch(max(abs(solve(fit$Hessian, fit$gradient))),
    error = function(e) Inf) else NA_real_
  row$Status <- if (any(unlist(fit$aliased))) "aliased_parameters" else
    if (!finite) "nonfinite_calculation" else
      if (fit$convergence$code != 0L) "nonconvergence_or_nonunique" else
        if (length(warnings)) "backend_warning" else
          if (row$MinEigenvalue <= sqrt(.Machine$double.eps) * max(1, max(eigenvalues))) "ill_conditioned_hessian" else
            if (row$MaxGradient > 1e-6 || !is.finite(row$MaxNewtonStep) || row$MaxNewtonStep > 1e-4)
              "stationarity_not_established" else "available"
  list(row = row, coefficients = fit$coefficients)
}

mfrm_ordinal_lr <- function(reduced, full) {
  out <- data.frame(ChiSquared = NA_real_, DF = NA_real_, PValue = NA_real_, Status = "model_unavailable")
  if (reduced$Status != "available" || full$Status != "available") return(out)
  delta <- full$LogLik - reduced$LogLik
  tol <- 64 * .Machine$double.eps * max(1, abs(c(reduced$LogLik, full$LogLik)))
  df <- full$Parameters - reduced$Parameters
  if (delta < -tol || df <= 0) { out$Status <- "invalid_nested_comparison"; return(out) }
  out$ChiSquared <- 2 * max(0, delta); out$DF <- df
  out$PValue <- stats::pchisq(out$ChiSquared, df, lower.tail = FALSE)
  out$Status <- "available"
  out
}
