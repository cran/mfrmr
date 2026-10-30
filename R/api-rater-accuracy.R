#' Estimate criterion-referenced rater accuracy from exact matches
#'
#' Fit a binary many-facet Rasch model to agreement with supplied reference
#' ratings. Raters are the objects of measurement; higher accuracy logits
#' indicate a greater probability of matching the reference, after adjusting
#' for the difficulty facets.
#'
#' @param data A data frame with one row per rater and benchmark combination.
#' @param rater Name of the rater identifier column.
#' @param facets Names of one or more difficulty-facet columns. Their joint
#'   combination identifies a benchmark, for example `c("Essay", "Domain")`.
#'   For a separate difficulty for every performance-domain combination,
#'   supply a single preconstructed benchmark identifier instead. These are
#'   different models; choose the structure substantively.
#' @param score,reference Names of aligned observed and reference rating
#'   columns. Both must be numeric, or both character/factor category labels.
#'   Numeric equality is exact; factors are compared by labels. Scores must
#'   already be categorized. No rounding, averaging or partial credit occurs.
#' @param method Ordinary `"JML"` (default) or `"MML"`; see Details.
#' @param maxit,reltol,quad_points Optimization and quadrature controls passed
#'   to [fit_mfrm()]. Quadrature applies only to MML.
#'
#' @details
#' Following the dichotomous example in Engelhard, Wang and Wind (2018,
#' pp. 42-43), each complete pair becomes `Match = 1` for equality and `0`
#' otherwise. The additive model is
#' `logit Pr(Match = 1) = accuracy[rater] - sum(difficulty[facet level])`.
#' Each difficulty facet is sum-to-zero centered; accuracy locations are free.
#' There is no free binary threshold. This function reuses the native RSM
#' estimator on the declared 0/1 scale, with unit weights, no interactions,
#' anchors, shrinkage or JML correction. It does not reproduce the paper's
#' interaction tests or its FACETS results.
#'
#' References are treated as fixed supplied benchmarks. Their quality,
#' selection and adjudication require substantive justification; their
#' uncertainty is not estimated. A benchmark must have at most one distinct
#' nonmissing reference across raters, even on rows with missing observed
#' scores. Repeated rater-by-benchmark rows are rejected, including repetitions
#' with a missing score. Choose independent observations and an appropriate
#' facet structure before fitting. Identifiers must be present and nonblank.
#' Distinct numeric identifiers that become the same text label are rejected;
#' provide unambiguous character identifiers instead.
#' Missing observed/reference ratings are excluded with separate indicators
#' in `cases`; missing references are never filled from other rows. Both
#' matches and mismatches must remain. Raters with no retained pairs still
#' appear in `table` with unavailable estimates.
#'
#' JML estimates a fixed accuracy location per rater and has the usual
#' incidental-parameter and sparse-design limitations. All-match/all-mismatch
#' raters retain primary `Inf`/`-Inf` estimates, not the finite optimizer
#' proxies. MML instead estimates a normal population of rater accuracies
#' and returns rater EAPs with `PosteriorSD`. Its effective population sample
#' size is the number of raters, not the number of benchmark ratings.
#' Posterior SD conditions on fitted calibration/population parameters; it
#' is not a frequentist SE or a confidence interval. Small rater pools and
#' population misspecification need particular review. No JML SE, interval,
#' significance test or automatic accuracy cutoff is computed here.
#' MML explicitly uses `population_formula = ~ 1` to estimate the rater
#' population mean and variance; it does not use the RSM default fixed N(0,1).
#'
#' Exact agreement proportions describe the retained benchmark mix. Accuracy
#' logits adjust for the specified difficulty structure under the response
#' model, local independence and a connected identifiable design. They are
#' neither ordinary rater severity nor examinee proficiency. Binary coding
#' discards the direction and magnitude of disagreements. The native fit's
#' readiness and estimate-use labels remain visible; finite numbers do not
#' establish reportability, model fit or transport to other benchmarks.
#' High separation reliability for these accuracy locations would describe
#' differences in accuracy, not agreement or differences in severity.
#'
#' @return A plain list with `table` (rater counts, exact agreement, accuracy,
#'   posterior SD and source status), `difficulty` (facet locations), `cases`
#'   (all input rows, binary matches and missingness), `fit` (native
#'   `mfrm_fit`), `settings` and `notes`. Internal `Person` identifiers in
#'   `fit` refer to raters; `Facet1`, `Facet2`, etc. refer to the requested
#'   difficulty facets. `settings$column_map` records this mapping, also used
#'   by `cases`. Interpret native fit reports by these roles, not their default
#'   examinee/rater labels. `difficulty$Facet` restores the input facet names.
#'   Save the complete result with [base::saveRDS()]; the data frames can be
#'   exported with [utils::write.csv()]. Reading saved output does not refit.
#'
#' @references
#' Engelhard, G., Jr., Wang, J., & Wind, S. A. (2018). A tale of two models:
#' Psychometric and cognitive perspectives on rater-mediated assessments using
#' accuracy ratings. *Psychological Test and Assessment Modeling, 60*(1), 33-52.
#' Equations and binary coding: pp. 42-43.
#' @seealso [fit_mfrm()], [score_comparison_table()], [rater_centrality_table()]
#' @examples
#' # A small illustrative benchmark exercise, not a validation sample.
#' d <- expand.grid(Judge = paste0("J", 1:4), Benchmark = paste0("B", 1:6))
#' d$Reference <- rep(c(1, 2, 3, 1, 2, 3), each = 4)
#' matched <- c(1,1,0,0, 1,0,1,0, 0,1,0,1,
#'              1,1,0,0, 0,1,1,0, 0,0,1,1)
#' d$Rating <- ifelse(matched == 1, d$Reference, (d$Reference %% 3) + 1)
#' accuracy <- fit_rater_accuracy(d, "Judge", "Benchmark", "Rating", "Reference")
#' accuracy$table
#' accuracy$difficulty
#' # Review SourceFitReadiness and EstimateUse, especially for small samples.
#' @export
fit_rater_accuracy <- function(data, rater, facets, score, reference,
                               method = c("JML", "MML"), maxit = 400,
                               reltol = 1e-9, quad_points = 31) {
  method <- match.arg(method)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) stop(
    "`data` must be a nonempty data frame with unique column names.", call. = FALSE)
  columns <- list(rater = rater, facets = facets, score = score, reference = reference)
  for (nm in names(columns)) {
    x <- columns[[nm]]
    if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) ||
        (nm != "facets" && length(x) != 1L)) stop(
      "`", nm, "` must name ", if (nm == "facets") "one or more columns." else "one column.",
      call. = FALSE)
  }
  selected <- unname(unlist(columns))
  if (anyDuplicated(selected) || !all(selected %in% names(data))) stop(
    "Rater, facets, score and reference must name distinct existing columns.", call. = FALSE)
  ids <- data[c(rater, facets)]
  for (nm in names(ids)) {
    v <- ids[[nm]]
    if (!(is.numeric(v) || is.character(v) || is.factor(v)) || !is.null(dim(v)) ||
        anyNA(v) || any(!nzchar(trimws(as.character(v)))) ||
        (is.numeric(v) && any(!is.finite(v)))) stop(
      "Identifiers must be finite, nonmissing and nonblank: ", nm, ".", call. = FALSE)
    labels <- as.character(v)
    if (length(unique(labels)) != length(unique(v))) stop(
      "Distinct identifiers in `", nm, "` have identical text labels; supply unique character identifiers.",
      call. = FALSE)
    ids[[nm]] <- labels
  }
  observed <- data[[score]]; target <- data[[reference]]
  numeric_pair <- is.numeric(observed) && is.numeric(target)
  label_pair <- (is.character(observed) || is.factor(observed)) &&
    (is.character(target) || is.factor(target))
  if ((!numeric_pair && !label_pair) || !is.null(dim(observed)) || !is.null(dim(target))) stop(
    "Score and reference must both be numeric or both be character/factor categories.", call. = FALSE)
  if (label_pair) {
    observed <- as.character(observed); target <- as.character(target)
    if (any(!is.na(observed) & !nzchar(trimws(observed))) ||
        any(!is.na(target) & !nzchar(trimws(target)))) stop(
      "Blank rating labels are not categories; use NA for missing ratings.", call. = FALSE)
  } else if (any(!is.na(observed) & !is.finite(observed)) ||
             any(!is.na(target) & !is.finite(target))) stop(
    "Nonmissing numeric ratings must be finite.", call. = FALSE)
  if (anyDuplicated(ids)) stop(
    "Repeated rater-by-benchmark rows are not supported; do not silently average or duplicate ratings.",
    call. = FALSE)
  facet_names <- paste0("Facet", seq_along(facets))
  names(ids) <- c("Rater", facet_names)
  known <- !is.na(target)
  references <- unique(data.frame(ids[known, facet_names, drop = FALSE],
                                  Reference = target[known], check.names = FALSE))
  if (anyDuplicated(references[facet_names])) stop(
    "Conflicting reference ratings within a benchmark combination.", call. = FALSE)
  missing_score <- is.na(observed); missing_reference <- is.na(target)
  included <- !missing_score & !missing_reference
  binary <- as.integer(observed == target)
  if (length(unique(binary[included])) != 2L) stop(
    "Retained ratings must contain both matches and mismatches to estimate accuracy.", call. = FALSE)
  cases <- data.frame(Row = seq_len(nrow(data)), ids, Score = observed,
    Reference = target, Match = binary, MissingScore = missing_score,
    MissingReference = missing_reference, Included = included, check.names = FALSE)
  model_data <- data.frame(ids[included, , drop = FALSE], Match = binary[included])
  fit <- fit_mfrm(model_data, person = "Rater", facets = facet_names, score = "Match",
    rating_min = 0, rating_max = 1, keep_original = TRUE, model = "RSM", method = method,
    maxit = maxit, reltol = reltol, quad_points = quad_points,
    population_formula = if (method == "MML") stats::as.formula("~ 1", env = baseenv()) else NULL,
    person_data = if (method == "MML") unique(model_data["Rater"]) else NULL,
    person_id = if (method == "MML") "Rater" else NULL)
  estimates <- fit$facets$person
  raters <- unique(ids$Rater)
  rows <- lapply(raters, function(id) {
    at <- ids$Rater == id
    retained <- sum(at & included)
    j <- match(id, as.character(estimates$Person))
    data.frame(Rater = id, InputRows = sum(at), RatingRows = retained,
      ExcludedRows = sum(at & !included), Matches = sum(binary[at & included]),
      ExactAgreement = if (retained) mean(binary[at & included]) else NA_real_,
      Accuracy = estimates$Estimate[j],
      PosteriorSD = if (method == "MML") estimates$SD[j] else NA_real_,
      Extreme = if (retained) as.character(estimates$Extreme[j]) else "unavailable",
      SourceFitReadiness = as.character(estimates$SourceFitReadiness[j]),
      SourceInferenceReady = if (retained) estimates$SourceInferenceReady[j] else FALSE,
      EstimateUse = if (retained) as.character(estimates$EstimateUse[j]) else "no_retained_ratings")
  })
  difficulty <- as.data.frame(fit$facets$others)
  difficulty$Facet <- facets[match(difficulty$Facet, facet_names)]
  names(difficulty)[names(difficulty) == "Estimate"] <- "Difficulty"
  list(table = do.call(rbind, rows), difficulty = difficulty, cases = cases, fit = fit,
    settings = list(method = method, model = "binary_additive_Rasch_accuracy",
      estimate = if (method == "MML") "rater_EAP" else "rater_JML",
      population = if (method == "MML") list(
        mean = unname(fit$population$coefficients[1]),
        sd = sqrt(fit$population$sigma2), raters = nrow(estimates),
        basis = "estimated_normal_rater_population") else NULL,
      reference = "fixed_supplied_ratings", uncertainty = if (method == "MML")
        "posterior_SD_given_fitted_parameters" else "not_computed",
      column_map = data.frame(Role = c("rater", rep("difficulty", length(facets))),
        Input = c(rater, facets), Model = c("Person", facet_names),
        Cases = c("Rater", facet_names))),
    notes = c(
      "Accuracy is reference-specific exact matching, not severity or examinee proficiency; disagreement direction and size are discarded.",
      "References are fixed supplied benchmarks; their uncertainty and validity are not estimated.",
      "Review source readiness, design connectivity, local independence and benchmark selection before interpretation; no accuracy cutoff or interval is supplied.",
      "Native fit Person fields refer to raters; each difficulty facet is sum-to-zero centered. Separation concerns differences in accuracy.",
      if (method == "MML") paste(
        "MML assumes a normal rater-accuracy population. EAP and posterior SD condition on fitted parameters;",
        "the population sample size is the number of raters, and posterior SD is not a frequentist SE.") else
        "JML accuracy is exploratory; incidental-parameter and sparse-design limitations remain. Infinite primary estimates are retained."))
}
