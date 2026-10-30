# Known conditional calibrations exercise cell arithmetic without fitting.
bias_numeric_reference <- function(eta, y, model = "RSM", slope = rep(1, length(y)),
                                   rater = rep("r", length(y)),
                                   criterion = rep("c", length(y)), max_abs = 10,
                                   weight = rep(1, length(y))) {
  config <- list(model = model, method = "JML", facet_names = c("Rater", "Criterion"),
    facet_signs = c(Rater = -1, Criterion = -1))
  prep <- list(rating_min = 0, rating_max = 1,
    levels = list(Person = "p", Rater = unique(rater), Criterion = unique(criterion)))
  obs <- data.frame(Person = "p", Rater = rater, Criterion = criterion,
    Observed = y, Weight = weight)
  measures <- data.frame(Facet = c(rep("Rater", length(unique(rater))),
    rep("Criterion", length(unique(criterion)))),
    Level = c(unique(rater), unique(criterion)), Estimate = 0, SE = .1)
  local_mocked_bindings(
    build_param_sizes = function(...) list(),
    expand_params = function(...) list(theta = 0, steps = 0,
      steps_mat = matrix(0, 1, 1), slopes = slope),
    build_indices = function(...) list(score_k = y, weight = weight,
      step_idx = rep(1L, length(y)), slope_idx = seq_along(y)),
    compute_eta = function(...) eta, .package = "mfrmr")
  out <- estimate_bias_interaction(list(config = config, prep = prep, opt = list(par = 0)),
    list(obs = obs, measures = measures), interaction_facets = c("Rater", "Criterion"),
    omit_extreme = FALSE, max_abs = max_abs)
  class(out) <- c("mfrm_bias", class(out))
  out
}

test_that("bounded bias points do not imply available interior screening statistics", {
  for (model in c("RSM", "PCM", "GPCM")) for (y in list(c(1, 1), c(0, 0))) {
    b <- bias_numeric_reference(c(0, .2), y, model = model, max_abs = 2)
    tab <- b$table
    expect_equal(tab$`Bias Size`, if (y[1] == 1) 2 else -2)
    expect_identical(tab$OptimizationStatus, "ok")
    expect_match(tab$BiasStatus, "search_boundary")
    expect_match(tab$ScreenReason, "search bound")
    expect_true(all(is.na(unlist(tab[c("S.E.", "t", "Prob.")]))))
    if (model == "GPCM") {
      expect_true(is.na(tab$`LR Prob.`))
      expect_true(is.na(tab$`Profile CI Lower`))
      expect_match(tab$`Profile CI Status`, "search boundary")
    }
    overview <- summary(b)$overview
    expect_equal(overview$Points, 1)
    expect_equal(overview$ScreenAvailable, 0)
    expect_true(is.na(overview$ScreenPositive))
    expect_true(is.na(overview$HolmScreenPositive))
  }
})

test_that("bias information includes small positive variances without partial-row fit", {
  b <- bias_numeric_reference(c(-30, 0, 0), c(0, 0, 1))
  tab <- b$table
  p <- plogis(c(-30, 0, 0) + tab$`Bias Size`)
  expect_equal(tab$ConditionalInformation, sum(p * (1 - p)), tolerance = 1e-14)
  expect_equal(tab$`S.E.`, 1 / sqrt(sum(p * (1 - p))), tolerance = 1e-14)
  expect_identical(tab$BiasStatus, "interior")
  expect_identical(tab$ScreenStatus, "available_conditional")
  expect_true(is.na(tab$Infit))
  expect_true(is.na(tab$Outfit))
  expect_equal(tab$FitUnavailableRows, 1)
  expect_match(tab$FitReason, "no partial-row")
})

test_that("interior bias points and curvature agree with an independent binary score", {
  for (model in c("RSM", "PCM", "GPCM")) {
    a <- if (model == "GPCM") c(.7, 1.4, 2) else rep(1, 3)
    eta <- c(-.2, .3, -.4); y <- c(0, 1, 1)
    score <- function(b) sum(a * (y - plogis(a * (eta + b))))
    expected <- uniroot(score, c(-10, 10), tol = 1e-12)$root
    b <- bias_numeric_reference(eta, y, model, a)$table
    p <- plogis(a * (eta + b$`Bias Size`))
    information <- sum(a^2 * p * (1 - p))
    expect_equal(b$`Bias Size`, expected, tolerance = 1e-4)
    expect_equal(b$ConditionalInformation, information, tolerance = 1e-13)
    expect_equal(b$`S.E.`, 1 / sqrt(information), tolerance = 1e-13)
    expect_identical(b$BiasStatus, "interior")
    expect_identical(b$FitStatus, "available_descriptive")
  }
})

test_that("batched bias fits keep weighted, failed and incomplete cells separate", {
  calls <- 0L
  local_mocked_bindings(optimize = function(...) {
    calls <<- calls + 1L
    if (calls %% 3L == 2L) stop("Injected cell failure")
    stats::optimize(...)
  })
  eta <- c(-.2, .3, -.4, .1, .2, .3, -30, 0, 0)
  y <- c(0, 1, 1, 1, 0, 1, 0, 0, 1)
  w <- c(1, .5, 2, 2, 1, 1, .5, 1, 3)
  b <- bias_numeric_reference(eta, y, rater = rep(c("a", "b", "c"), each = 3), weight = w)
  a <- b$table[b$table$FacetA_Level == "a", ]
  p <- plogis(eta[1:3] + a$`Bias Size`)
  variance <- p * (1 - p)
  residual <- y[1:3] - p
  expect_equal(a$Infit, sum(w[1:3] * residual^2) / sum(w[1:3] * variance), tolerance = 1e-13)
  expect_equal(a$Outfit, sum(w[1:3] * residual^2 / variance) / sum(w[1:3]), tolerance = 1e-13)
  expect_identical(a$FitStatus, "available_descriptive")
  failed <- b$table[b$table$FacetA_Level == "b", ]
  expect_identical(failed$OptimizationStatus, "failed")
  expect_match(failed$FitReason, "cell estimate is unavailable")
  partial <- b$table[b$table$FacetA_Level == "c", ]
  expect_identical(partial$OptimizationStatus, "ok")
  expect_equal(partial$FitAvailableRows, 2)
  expect_equal(partial$FitUnavailableRows, 1)
  expect_match(partial$FitReason, "no partial-row")
  expect_true(all(is.na(c(failed$Infit, failed$Outfit, partial$Infit, partial$Outfit))))
  expect_equal(b$table$FitResponseRows, rep(3, 3))
})

test_that("all failed bias cells retain fit availability without empty aggregation", {
  local_mocked_bindings(optimize = function(...) stop("Injected cell failure"))
  b <- bias_numeric_reference(rep(0, 4), c(0, 1, 0, 1), rater = c("a", "a", "b", "b"))
  expect_identical(b$table$OptimizationStatus, rep("failed", 2))
  expect_equal(b$table$FitUnavailableRows, c(2, 2))
  expect_equal(b$table$FitAvailableRows, c(0, 0))
  expect_true(all(is.na(b$table$Infit)))
})

test_that("bias cells preserve facet labels containing the former separator", {
  b <- bias_numeric_reference(rep(0, 4), c(0, 1, 1, 0),
    rater = c("a||b", "a", "a||b", "a"),
    criterion = c("c", "b||c", "c", "b||c"))
  expect_equal(nrow(b$table), 2)
  expect_equal(b$table$ObsN, c(2, 2))
  actual <- b$table[, c("FacetA_Level", "FacetB_Level")]
  expect_setequal(paste(actual$FacetA_Level, actual$FacetB_Level, sep = "/"),
    c("a||b/c", "a/b||c"))
  expect_equal(b$table$FacetA_Measure, c(0, 0))
  expect_equal(b$table$FacetB_Measure, c(0, 0))
  b$table$`Bias Size` <- c(.2, -.3)
  d <- data.frame(Rater = c("a||b", "a", "a||b", "a"),
    Criterion = c("c", "b||c", "c", "b||c"))
  adjustment <- compute_bias_adjustment_vector(list(prep = list(data = d)), b)
  expect_equal(adjustment, b$table$`Bias Size`[match(d$Rater, b$table$FacetA_Level)])
})

test_that("unavailable bias information retains the bounded point and explains missing screens", {
  local_mocked_bindings(mfrm_category_variance = function(...) c(NA_real_, .25),
    .package = "mfrmr")
  b <- bias_numeric_reference(c(0, 0), c(0, 1))
  expect_true(is.finite(b$table$`Bias Size`))
  expect_true(is.na(b$table$`S.E.`))
  expect_true(is.na(b$table$Outfit))
  expect_identical(b$table$ScreenStatus, "unavailable_information")
  expect_match(b$table$ScreenReason, "complete")
})

test_that("heterogeneity screens require the full retained cell family", {
  b <- bias_numeric_reference(rep(0, 4), c(0, 1, 1, 1),
    rater = c("a", "a", "b", "b"))
  expect_equal(sum(is.finite(b$table$`S.E.`)), 1)
  expect_true(is.na(b$chi_sq$FixedChiSq))
  expect_true(is.na(b$chi_sq$FixedDF))
  expect_identical(b$chi_sq$Status, "unavailable")
  expect_equal(b$chi_sq$AvailableCells, 1)
  expect_match(b$chi_sq$Reason, "no selected-subset")
})

test_that("bias multiplicity summaries retain unavailable cells in the family size", {
  b <- bias_numeric_reference(rep(0, 6), rep(c(0, 1), 3),
    rater = rep(c("a", "b", "c"), each = 2))
  b$table$`Prob.` <- c(.02, NA_real_, .8)
  overview <- summary(b)$overview
  expect_equal(overview$ScreenAvailable, 2)
  expect_equal(overview$ScreenUnavailable, 1)
  expect_equal(overview$ScreenPositive, 1)
  expect_equal(overview$BonferroniScreenPositive, 0)
  expect_equal(overview$HolmScreenPositive, 0)
  expect_match(paste(summary(b)$notes, collapse = " "), "do not establish formal")
  b$table$BiasStatus <- NULL
  expect_match(paste(summary(b)$notes, collapse = " "), "predates")
})
