#' Compute pairwise weighted Cohen's kappa
#'
#' Compare observed scores on matched units for every pair of raters, using
#' quadratic (QWK), linear or unweighted agreement. No model is fitted.
#'
#' @param data Long-format ratings with at most one row per rater and unit.
#' @param rater Name of the rater identifier column.
#' @param unit Names of columns jointly identifying a rated unit, such as
#'   `c("Person", "Task")`.
#' @param score Name of the score column. Values may be numeric, character,
#'   factor or logical, with `NA` for missing scores. Factors use their labels.
#' @param categories Required distinct permitted category values, including
#'   unused categories, in their substantive order. At least two are required.
#'   They must have the same numeric/character/logical representation as scores.
#'   With linear/quadratic weights, numeric values use their actual spacing;
#'   character/logical values use equally spaced positions in this vector.
#'   Factor levels never silently determine the order. To treat numeric codes
#'   as equally spaced ranks, explicitly convert scores and categories to
#'   character. Choose spacing on substantive grounds.
#' @param weights One of `"quadratic"` (default), `"linear"` or `"unweighted"`.
#'   Agreement weights are 1 minus squared normalized distance, 1 minus
#'   normalized absolute distance, or the identity matrix, respectively.
#'
#' @details
#' For each pair, both observed and expected agreement use only common complete
#' units. Expected cell proportions are products of the two rater-specific
#' margins on this same sample. Kappa is `(observed - expected)/(1 - expected)`.
#' This complete-case convention differs from the all-available rater margins
#' in Gwet's incomplete-data Example 4.1. Different pairs may retain different
#' cohorts; comparing or averaging their kappas does not establish a panel
#' coefficient. At least two complete units are required, as in the existing
#' QWK score-comparison interface.
#'
#' Computation uses the equivalent disagreement ratio, with distances rescaled
#' to the categories actually present in each retained pair. The common scale
#' cancels, avoiding loss of the coefficient when an unused numeric category
#' is very distant. `ScaledObservedDisagreement` and
#' `ScaledExpectedDisagreement` retain those ingredients; agreement columns
#' use the declared full-scale weights. Agreement columns can round to 1 at
#' extreme relative spacings even when kappa is defined. Negative estimates
#' remain visible; zero expected disagreement returns `NA` with a status.
#'
#' Duplicate rater-unit cells are rejected before dropping missing scores.
#' The case audit distinguishes an absent assignment from an explicit missing
#' score and retains the union of units assigned to either member of the pair.
#' No imputation, confidence intervals, tests or quality cutoffs are supplied.
#' Missingness may affect representativeness. Agreement does not establish
#' reference accuracy or MFRM rater-severity homogeneity. Storage of the weight
#' matrix requires space quadratic in the number of declared categories;
#' these are categorical-score interfaces, not continuous-score ICCs.
#'
#' @return A plain list with `table` (one row per pair), `counts` (nonzero
#'   retained contingency cells, with category indices), `categories`
#'   (indices, values and weighting positions), `weights` (agreement matrix),
#'   `cases`, `units`, `settings` and `notes`. `Unit` indexes `units`.
#'   Save with [base::saveRDS()] or export data frames with [utils::write.csv()].
#' @references
#' Cohen, J. (1968). Weighted kappa: nominal scale agreement with provision
#' for scaled disagreement or partial credit. *Psychological Bulletin*, 70,
#' 213-220. \doi{10.1037/h0026256}.
#'
#' Gwet, K. L. *Handbook of Inter-Rater Reliability*, fifth edition, volume 1,
#' chapter 4, pp. 100-110, equations 4.2.4-4.2.6.
#' \url{https://www.agreestat.com/assets/cac5-chapter4-excerpt.pdf}
#' @seealso [gwet_ac2_table()], [krippendorff_alpha_table()],
#'   [pairwise_icc_table()], [score_comparison_table()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("pairwise_kappa_table", run.dontrun = TRUE)
#' \dontrun{
#' data("mfrmr_example_agreement", package = "mfrmr")
#' ratings <- mfrmr_example_agreement
#' head(ratings) # One rater's score for one Person-Task response.
#' pairs <- pairwise_kappa_table(ratings, rater = "Rater",
#'   unit = c("Person", "Task"), score = "Score", categories = 0:3)
#' pairs$table
#' pairs$cases # Included marks the shared nonmissing ratings for each pair.
#' # R2 gives every response one category more than R1: equal rankings do not
#' # imply agreement. Do not treat two Tasks by the same Person as one response.
#' }
#' @export
pairwise_kappa_table <- function(data, rater, unit, score, categories,
                                 weights = c("quadratic", "linear", "unweighted")) {
  weights <- match.arg(weights)
  roster <- mfrm_rating_roster(data, rater, unit, score)
  scale <- mfrm_agreement_scale(data[[score]], categories, weights)
  y <- scale$scores; code <- scale$code; q <- length(scale$categories)
  u <- roster$unit_id
  pairs <- utils::combn(roster$raters, 2, simplify = FALSE)
  tables <- cases <- counts <- vector("list", length(pairs))
  for (i in seq_along(pairs)) {
    p <- pairs[[i]]
    a <- which(as.character(roster$ids[[rater]]) == p[1])
    b <- which(as.character(roster$ids[[rater]]) == p[2])
    scope <- sort(unique(c(u[a], u[b])))
    a <- a[match(scope, u[a])]; b <- b[match(scope, u[b])]
    complete <- !is.na(a) & !is.na(b) & !is.na(y[a]) & !is.na(y[b])
    n <- sum(complete)
    cases[[i]] <- data.frame(Pair = i, Rater1 = p[1], Rater2 = p[2], Unit = scope,
      Row1 = a, Row2 = b, Score1 = y[a], Score2 = y[b],
      Assigned1 = !is.na(a), Assigned2 = !is.na(b),
      Missing1 = !is.na(a) & is.na(y[a]), Missing2 = !is.na(b) & is.na(y[b]),
      Included = complete)
    x <- code[a[complete]]; z <- code[b[complete]]
    cells <- table(factor(x, levels = seq_len(q)), factor(z, levels = seq_len(q)))
    nonzero <- which(cells > 0, arr.ind = TRUE)
    counts[[i]] <- data.frame(Pair = rep(i, nrow(nonzero)), Category1 = nonzero[, 1],
      Category2 = nonzero[, 2], N = as.integer(cells[nonzero]))
    pa <- pe <- observed <- expected <- value <- NA_real_
    status <- "fewer_than_two_complete_units"
    if (n >= 2L) {
      margin1 <- rowSums(cells) / n; margin2 <- colSums(cells) / n
      pa <- 1 - sum(scale$distance * (cells / n))
      pe <- 1 - sum(scale$distance * outer(margin1, margin2))
      active <- sort(unique(c(x, z)))
      distance <- mfrm_agreement_distance(scale$positions[active], weights)
      observed <- sum(distance * (cells[active, active, drop = FALSE] / n))
      expected <- sum(distance * outer(margin1[active], margin2[active]))
      status <- if (expected > 0) "available" else "zero_expected_disagreement"
      if (status == "available") value <- 1 - observed / expected
    }
    tables[[i]] <- data.frame(Pair = i, Rater1 = p[1], Rater2 = p[2], N = n,
      CandidateUnits = length(scope), ExcludedUnits = length(scope) - n,
      Kappa = value, Status = status, Weighting = weights,
      ObservedAgreement = pa, ExpectedAgreement = pe,
      ScaledObservedDisagreement = observed, ScaledExpectedDisagreement = expected)
  }
  list(table = do.call(rbind, tables), counts = do.call(rbind, counts),
    categories = scale$table, weights = scale$weights, cases = do.call(rbind, cases),
    units = roster$units,
    settings = list(rater = rater, unit = unit, score = score, weights = weights,
      spacing = scale$spacing, missing = "pairwise_complete_units",
      uncertainty = "not_computed"),
    notes = c("Both margins use the same common complete units as observed agreement.",
      "Different pair cohorts affect comparability; pairwise kappas are not a panel coefficient.",
      "Agreement is distinct from reference accuracy, MFRM fit and separation reliability.",
      "No confidence intervals, tests, quality cutoffs or correction for selective missingness are supplied."))
}

#' Compute Gwet's AC2 or unweighted AC1
#'
#' Compute observed-score agreement across two or more raters, with a declared
#' category universe and quadratic, linear or unweighted agreement.
#'
#' @inheritParams pairwise_kappa_table
#' @details
#' Observed agreement is the equally weighted mean of within-unit weighted
#' agreement over distinct ordered rating pairs, using only units with at least
#' two nonmissing ratings. Category proportions are the equally weighted mean
#' of within-unit category proportions over all units with at least one rating.
#' Thus singletons contribute to expected agreement but not observed agreement;
#' units with no ratings contribute to neither. This is Gwet's raw-rating
#' convention, not Krippendorff's pairable-rating margins or Cohen's separate
#' rater margins. Unequal rating counts do not give units greater weight.
#'
#' With q declared categories, agreement matrix W and mean category proportions
#' pi, expected agreement is `sum(W) * sum(pi * (1 - pi)) / (q * (q - 1))`.
#' AC is `(observed - expected)/(1 - expected)`. Quadratic/linear weights give
#' AC2; identity weights give AC1. At two categories the supplied linear and
#' quadratic weights equal identity weights, so the numerical results coincide.
#' The complete category universe is required: adding unused categories can
#' change expected agreement and AC. It must not be selected to improve a result.
#'
#' Negative estimates are retained. No pairable units or a nonpositive chance
#' denominator gives `NA` with a status. With all raters always using the same
#' category, AC can be 1 when the declared universe has at least two categories;
#' this does not establish discrimination, validity or reference accuracy.
#' Unlike kappa, its chance model does not become undefined in that case.
#'
#' Repeated rater-unit cells are rejected, even when scores are missing. Audits
#' distinguish explicit missing scores from absent assignments relative to the
#' observed rater roster. Units with no input rows are unknown. Missingness and
#' assignment may still bias interpretation. Results are descriptive point
#' estimates: no confidence intervals, tests, quality bands or invariance to
#' prevalence are claimed. Different chance models answer different questions;
#' AC2 is not an automatic replacement or correction for kappa. Dense weight
#' storage is quadratic in the declared category count. Custom weights and
#' circular, bipolar, radical and ratio weights are not implemented.
#'
#' @return A plain list with `table` (coefficient, observed/expected agreement,
#'   status and counts), `counts` (category totals and unit-mean proportions),
#'   `categories`, `weights`, `cases` (input rows and their two inclusion flags),
#'   `unit_summary`, `units`, `settings` and `notes`. `Unit` indexes `units`.
#'   Results can be saved with [base::saveRDS()] without refitting on reload.
#' @references
#' Gwet, K. L. Author implementation of AC1/AC2, `gwet.ac1.raw`,
#' and linear/quadratic weights, in *irrCAC*.
#' \url{https://github.com/kgwet/irrCAC/blob/master/R/agree.coeff3.raw.r}
#' \url{https://github.com/kgwet/irrCAC/blob/master/R/weights.gen.r}
#' @seealso [pairwise_kappa_table()], [krippendorff_alpha_table()],
#'   [pairwise_icc_table()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("gwet_ac2_table", run.dontrun = TRUE)
#' \dontrun{
#' data("mfrmr_example_agreement", package = "mfrmr")
#' ratings <- mfrmr_example_agreement
#' ac <- gwet_ac2_table(ratings, rater = "Rater",
#'   unit = c("Person", "Task"), score = "Score", categories = 0:3)
#' ac$table
#' ac$unit_summary # Number of observed ratings and pairability for each response.
#' ac$notes # Read the panel coefficient together with its assumptions.
#' }
#' @export
gwet_ac2_table <- function(data, rater, unit, score, categories,
                            weights = c("quadratic", "linear", "unweighted")) {
  weights <- match.arg(weights)
  roster <- mfrm_rating_roster(data, rater, unit, score)
  scale <- mfrm_agreement_scale(data[[score]], categories, weights)
  y <- scale$scores; code <- scale$code; q <- length(scale$categories)
  u <- roster$unit_id; nu <- nrow(roster$units)
  observed <- tabulate(u[!is.na(y)], nbins = nu)
  assigned <- tabulate(u, nbins = nu)
  pairable <- observed >= 2L; rated <- observed > 0L; nrated <- sum(rated)
  pi <- numeric(q); unit_agreement <- rep(NA_real_, nu)
  groups <- split(code[!is.na(y)], u[!is.na(y)])
  for (name in names(groups)) {
    index <- as.integer(name); z <- groups[[name]]; m <- length(z)
    f <- tabulate(z, nbins = q)
    pi <- pi + f / m / nrated
    if (m >= 2L) unit_agreement[index] <- 1 -
      sum(scale$distance * outer(f / m, f / (m - 1)))
  }
  pa <- if (any(pairable)) mean(unit_agreement[pairable]) else NA_real_
  pe <- if (any(rated)) sum(scale$weights) / q / (q - 1) * sum(pi * (1 - pi)) else NA_real_
  status <- if (!any(pairable)) "no_pairable_ratings" else
    if (!is.finite(pa) || !is.finite(pe)) "nonfinite_calculation" else
      if (pe >= 1) "nonpositive_chance_denominator" else "available"
  value <- if (status == "available") (pa - pe) / (1 - pe) else NA_real_
  list(table = data.frame(Coefficient = if (weights == "unweighted") "AC1" else "AC2",
      AC = value, Status = status, Weighting = weights, ObservedAgreement = pa,
      ExpectedAgreement = pe, Raters = length(roster$raters), Units = nu,
      PairableUnits = sum(pairable), MarginalUnits = sum(rated),
      NonmissingRatings = sum(observed), MissingScores = sum(is.na(y)),
      SingletonUnits = sum(observed == 1L), UnratedUnits = sum(!rated)),
    counts = data.frame(Category = scale$categories,
      NonmissingRatings = tabulate(code, nbins = q),
      UnitMeanProportion = if (any(rated)) pi else rep(NA_real_, q)),
    categories = scale$table, weights = scale$weights,
    cases = data.frame(Row = seq_along(y), Rater = as.character(roster$ids[[rater]]),
      Unit = u, Score = y, MissingScore = is.na(y),
      IncludedObserved = !is.na(y) & pairable[u], IncludedMargins = !is.na(y)),
    unit_summary = data.frame(Unit = seq_len(nu), AssignedRows = assigned,
      NonmissingRatings = observed, MissingScores = assigned - observed,
      AbsentAssignments = length(roster$raters) - assigned,
      IncludedObserved = pairable, IncludedMargins = rated,
      ObservedAgreement = unit_agreement),
    units = roster$units,
    settings = list(rater = rater, unit = unit, score = score, weights = weights,
      spacing = scale$spacing, observed = "equal_weight_pairable_units",
      margins = "equal_weight_nonempty_units_including_singletons",
      uncertainty = "not_computed"),
    notes = c("Singletons enter category proportions but not observed agreement; empty units enter neither.",
      "Unused declared categories affect the chance model. Choose the category universe before inspecting agreement.",
      "AC1/AC2, kappa and alpha use different chance models and are not interchangeable.",
      "A high AC does not establish discrimination, reference accuracy or MFRM severity homogeneity.",
      "No confidence intervals, tests, quality bands or correction for selective missingness are supplied."))
}

mfrm_agreement_scale <- function(y, categories, weights) {
  if (is.factor(y)) y <- as.character(y)
  if (is.factor(categories)) categories <- as.character(categories)
  valid <- function(v) !is.object(v) && is.null(dim(v)) &&
    (is.character(v) || is.logical(v) || (is.numeric(v) && !is.complex(v))) &&
    (!is.numeric(v) || all(is.na(v) | is.finite(v))) &&
    (!is.character(v) || all(is.na(v) | nzchar(trimws(v))))
  if (!valid(y)) stop("Scores must be finite numbers or nonblank category labels, with NA for missing values.", call. = FALSE)
  if (!valid(categories) || length(categories) < 2L || anyNA(categories) ||
      anyDuplicated(categories) || is.numeric(categories) != is.numeric(y) ||
      is.character(categories) != is.character(y)) stop(
    "Supply at least two distinct nonmissing categories of the same representation as scores.", call. = FALSE)
  code <- match(y, categories)
  if (any(!is.na(y) & is.na(code))) stop("A score is outside the supplied categories.", call. = FALSE)
  positions <- if (is.numeric(categories)) as.numeric(categories) else as.numeric(seq_along(categories))
  distance <- mfrm_agreement_distance(positions, weights)
  agreement <- 1 - distance
  dimnames(agreement) <- list(as.character(seq_along(categories)), as.character(seq_along(categories)))
  list(scores = y, code = code, categories = categories, positions = positions,
    distance = distance, weights = agreement,
    spacing = if (weights == "unweighted") "identity" else
      if (is.numeric(categories)) "numeric_category_distance" else "declared_category_positions",
    table = data.frame(CategoryIndex = seq_along(categories), Category = categories,
      Position = if (weights == "unweighted") NA_real_ else positions))
}

mfrm_agreement_distance <- function(positions, weights) {
  q <- length(positions)
  if (weights == "unweighted") return(1 - diag(q))
  if (q == 1L) return(matrix(0, 1, 1))
  span <- diff(range(positions))
  # Half-scale first only when the raw span overflows; preserve tiny gaps otherwise.
  distance <- if (is.finite(span)) abs(outer(positions, positions, "-")) / span else
    abs(outer(positions / 2, positions / 2, "-")) / (max(positions) / 2 - min(positions) / 2)
  if (weights == "quadratic") distance <- distance^2
  distance
}
