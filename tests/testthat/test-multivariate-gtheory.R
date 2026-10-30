# Construct balanced data with exactly specified ANOVA covariance components.
# Orthonormal tensor contrasts provide independent mean-product targets.
mvgt_fixture <- function(indefinite = FALSE, rater_task = NULL) {
  n <- c(4L, 3L, 4L)
  sources <- c("Person", "Rater", "Task", "Person:Rater", "Person:Task", "Rater:Task", "Residual")
  subsets <- list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
  values <- list(c(1.2, .3, .3, .8), c(.2, .04, .04, .3),
    c(.3, -.02, -.02, .15), c(.4, .1, .1, .5),
    c(.3, -.06, -.06, .4), c(.15, .03, .03, .1), c(.6, .2, .2, .7))
  if (indefinite) values[[1L]] <- c(.01, .03, .03, .01)
  if (!is.null(rater_task)) values[[6L]] <- rater_task
  gamma <- setNames(lapply(values, matrix, nrow = 2L,
    dimnames = list(c("Content", "Organization"), c("Content", "Organization"))), sources)
  q <- lapply(n, function(k) qr.Q(qr(stats::contr.helmert(k))))
  average <- lapply(n, function(k) matrix(1 / sqrt(k), k, 1))
  y <- matrix(0, prod(n), 2)
  for (i in seq_along(subsets)) {
    s <- subsets[[i]]
    mp <- matrix(0, 2, 2)
    for (j in seq_along(subsets)) {
      if (all(s %in% subsets[[j]])) mp <- mp + prod(n[setdiff(1:3, subsets[[j]])]) * gamma[[j]]
    }
    b <- lapply(1:3, function(j) if (j %in% s) q[[j]] else average[[j]])
    basis <- kronecker(b[[3]], kronecker(b[[2]], b[[1]]))
    y <- y + basis[, 1:2, drop = FALSE] %*% chol(prod(n[s] - 1) * mp)
  }
  data <- expand.grid(Person = paste0("P", seq_len(n[1])),
    Rater = paste0("R", seq_len(n[2])), Task = paste0("T", seq_len(n[3])))
  data$Content <- y[, 1] + 10
  data$Organization <- y[, 2] + 20
  list(data = data, components = gamma, counts = n)
}

mvgt_common_task_data <- function() {
  path <- testthat::test_path("..", "..", "inst", "extdata", "mgenova-table12.csv")
  if (!file.exists(path)) path <- system.file("extdata", "mgenova-table12.csv", package = "mfrmr")
  read.csv(path)
}

test_that("D-study printing keeps score weights and source counts distinct from plans", {
  g <- mfrm_g_study(mvgt_common_task_data(), c("V", "W"), rater = NULL)
  weights <- cbind(Equal = c(V = .5, W = .5), Difference = c(V = -1, W = 1))
  d <- mfrm_d_study(g, data.frame(Tasks = c(3, 12)), weights)
  before <- serialize(d, NULL)
  output <- capture.output(print(d))
  expect_match(paste(output, collapse = "\n"), "Person = 10, Task = 6", fixed = TRUE)
  expect_match(paste(output, collapse = "\n"), "Source rows: 60 used of 60 ; 0 explicitly omitted.", fixed = TRUE)
  expect_true(any(grepl("Future complete balanced plans", output, fixed = TRUE)))
  expect_true(any(grepl("Composite weights (not normalized)", output, fixed = TRUE)))
  expect_true(all(capture.output(print(weights)) %in% output))
  expect_identical(summary(d), d$coefficients)
  expect_identical(serialize(d, NULL), before)
  vector <- mfrm_d_study(g, weights = weights[, "Difference"])
  expected <- matrix(c(-1, 1), 2, 1, dimnames = list(c("V", "W"), "Composite"))
  expect_true(all(capture.output(print(expected)) %in% capture.output(print(vector))))
  expect_false(any(grepl("Composite weights", capture.output(print(mfrm_d_study(g))), fixed = TRUE)))
})

test_that("common-task raw data reproduce mGENOVA Appendix E G-study and D-study values", {
  # Brennan (2001), Manual for mGENOVA 2.1, Table 12 (p. 32), Appendix E
  # (pp. 74-77): synthetic data from Generalizability Theory, Table 9.3.
  # All persons and items are shared across V and W. Here the item facet is Task.
  # This estimates matrices from the original scores, without substituting them.
  data <- mvgt_common_task_data()
  expect_equal(dim(data), c(60L, 4L))
  expect_equal(unname(colMeans(data[c("V", "W")])), c(4.516666666666667, 5.083333333333333))
  g <- mfrm_g_study(data, c("V", "W"), rater = NULL)
  expect_identical(g$design$counts, c(Person = 10L, Task = 6L))
  expect_null(g$design$rater)
  expect_equal(g$degrees_of_freedom$DF, c(9, 5, 45))
  expected_mp <- list(Person = c(3.46111, 2.62037, 2.62037, 3.75000),
    Task = c(4.69667, 1.62333, 1.62333, 4.73667),
    Residual = c(1.25222, .70481, .70481, 1.53667))
  expected_vc <- list(Person = c(.36815, .31926, .31926, .36889),
    Task = c(.34444, .09185, .09185, .32000),
    Residual = c(1.25222, .70481, .70481, 1.53667))
  # The manual prints only five decimals; upper-triangle component entries
  # can be correlations. Expected symmetric matrices use the covariances.
  expect_equal(lapply(g$mean_products, function(m) round(as.vector(m), 5)), expected_mp)
  expect_equal(lapply(g$components, function(m) round(as.vector(m), 5)), expected_vc)
  expect_true(all(g$component_diagnostics$PositiveSemidefinite))
  expect_output(print(g), "10 persons, 6 tasks")
  d <- mfrm_d_study(g, weights = c(V = -1, W = 1))
  expect_identical(d$design_grid, data.frame(Tasks = 6L))
  expect_false("Raters" %in% names(d$coefficients))
  expected <- rbind(
    c(.36815, .20870, .26611, .63820, .58044, .45684, .51586),
    c(.36889, .25611, .30944, .59022, .54382, .50607, .55628),
    c(.09852, .22988, .31000, .30000, .24116, .47945, .55678))
  measures <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
                "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  expect_equal(unname(round(as.matrix(d$coefficients[measures]), 5)), expected)
  expect_true(all(d$coefficients$Status == "Available"))
  expect_equal(unname(round(d$covariances[[1]]$RelativeError, 5)),
    matrix(c(.20870, .11747, .11747, .25611), 2))
  expect_equal(unname(round(d$covariances[[1]]$AbsoluteError, 5)),
    matrix(c(.26611, .13278, .13278, .30944), 2))
  expect_output(print(d), "Means over common random tasks")
})

test_that("one-facet projections preserve averaging, score identity and contrast meaning", {
  data <- mvgt_common_task_data()
  g <- mfrm_g_study(data, c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(3L, 6L, 12L)), c(V = -1, W = 1))
  composite <- d$coefficients[d$coefficients$Kind == "Composite", ]
  expect_equal(composite$UniverseVariance, rep(composite$UniverseVariance[2], 3))
  expect_equal(composite$RelativeErrorVariance, composite$RelativeErrorVariance[2] * c(2, 1, .5))
  expect_equal(composite$AbsoluteErrorVariance, composite$AbsoluteErrorVariance[2] * c(2, 1, .5))
  data$Difference <- data$W - data$V
  direct <- mfrm_d_study(mfrm_g_study(data, "Difference", rater = NULL),
    d$design_grid)$coefficients
  expect_equal(direct$G, composite$G)
  expect_equal(direct$AbsoluteSEM, composite$AbsoluteSEM)
  one <- mfrm_d_study(mfrm_g_study(data, "V", rater = NULL),
    d$design_grid)$coefficients
  expect_equal(one$Phi, d$coefficients$Phi[d$coefficients$Score == "V"])
  shuffled <- data[rev(seq_len(nrow(data))), ]
  names(shuffled)[1:2] <- c("Candidate", "Item")
  again <- mfrm_g_study(shuffled, c("W", "V"),
    person = "Candidate", rater = NULL, task = "Item")
  expect_equal(again$components, lapply(g$components, function(m) m[2:1, 2:1]))
  # Adding an unselected column is explicit and does not introduce a facet.
  data$Rater <- NA_character_
  expect_identical(mfrm_g_study(data, c("V", "W"), rater = NULL), g)
  # Missing rater identifiers must not silently switch the default model.
  expect_error(mfrm_g_study(data, c("V", "W")), "nonmissing, nonblank")
})

test_that("one-facet designs refuse incomplete cells and unavailable rater projections", {
  data <- mvgt_common_task_data()
  scores <- c("V", "W")
  expect_error(mfrm_g_study(data, scores), "column names")
  expect_error(mfrm_g_study(data[-1, ], scores, rater = NULL), "complete, balanced")
  expect_error(mfrm_g_study(rbind(data, data[1, ]), scores, rater = NULL), "Duplicate")
  expect_error(mfrm_g_study(subset(data, Task == 1), scores, rater = NULL), "at least two")
  expect_error(mfrm_g_study(data, scores, rater = NULL, task = NULL), "distinct")
  missing <- data
  missing$V[1] <- NA_real_
  expect_error(mfrm_g_study(missing, scores, rater = NULL), "finite numeric")
  # Multiple raters cannot be dropped or averaged implicitly.
  crossed <- mvgt_fixture()$data
  expect_error(mfrm_g_study(crossed, c("Content", "Organization"), rater = NULL), "Duplicate")
  g <- mfrm_g_study(data, scores, rater = NULL)
  for (grid in list(data.frame(Raters = 1, Tasks = 6), data.frame(Raters = 6),
      data.frame(Tasks = 0), data.frame(Tasks = 1.5), data.frame(Tasks = NA_real_),
      data.frame(Tasks = numeric()))) {
    expect_error(mfrm_d_study(g, grid), "positive integer Tasks|exactly one count|missing count|at least one design row")
  }
  # Negative universe variance affects only the corresponding score's coefficients.
  data$V <- data$V - ave(data$V, data$Person, FUN = mean)
  negative <- mfrm_g_study(data, scores, rater = NULL)
  expect_lt(negative$components$Person[1, 1], 0)
  projected <- mfrm_d_study(negative)$coefficients
  expect_true(is.na(projected$G[1]))
  expect_identical(projected$GStatus[1], "Negative universe variance")
  expect_true(is.finite(projected$G[2]))
})

test_that("balanced multivariate ANOVA recovers cross-covariances and matches a QR reference", {
  f <- mvgt_fixture()
  g <- mfrm_g_study(f$data, c("Content", "Organization"))
  expect_s3_class(g, "mfrm_multivariate_gstudy")
  expect_equal(g$components, f$components, tolerance = 1e-12)
  expect_true(all(g$component_diagnostics$PositiveSemidefinite))
  expect_equal(sum(g$degrees_of_freedom$DF), nrow(f$data) - 1)
  # Independent saturated linear-model QR effects plus EMS inversion.
  mm <- model.matrix(~ Person * Rater * Task, f$data)
  fit <- lm.fit(mm, as.matrix(f$data[c("Content", "Organization")]))
  assignment <- attr(mm, "assign")
  mp <- lapply(1:7, function(i) crossprod(fit$effects[assignment == i, , drop = FALSE]) / sum(assignment == i))
  subsets <- list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
  ems <- matrix(0, 7, 7)
  for (i in 1:7) for (j in 1:7) {
    if (all(subsets[[i]] %in% subsets[[j]])) ems[i, j] <- prod(f$counts[setdiff(1:3, subsets[[j]])])
  }
  reference <- solve(ems, t(vapply(mp, as.vector, numeric(4))))
  expect_equal(unname(t(vapply(g$components, as.vector, numeric(4)))), unname(reference), tolerance = 1e-12)
  expect_identical(summary(g), g$component_diagnostics)
  expect_output(print(g), "highest-order interaction and residual are combined")
  # Row order and literal identifier column names do not change estimates.
  shuffled <- f$data[rev(seq_len(nrow(f$data))), ]
  names(shuffled)[1:3] <- c("Participant id", "Assessor id", "Task id")
  again <- mfrm_g_study(shuffled, c("Content", "Organization"),
    "Participant id", "Assessor id", "Task id")
  expect_equal(again$components, g$components, tolerance = 1e-12)
})

test_that("large crossed plan counts give the same projection with integer storage", {
  f <- mvgt_fixture()
  g <- mfrm_g_study(f$data, c("Content", "Organization"))
  grid <- data.frame(Raters = c(30000L, 50000L), Tasks = c(30000L, 50000L))
  d <- expect_warning(mfrm_d_study(g, grid), NA)
  numeric_grid <- as.data.frame(lapply(grid, as.double))
  reference <- mfrm_d_study(g, numeric_grid)
  expect_identical(d$design_grid, grid)
  expect_identical(d$covariances, reference$covariances)
  expect_equal(d$coefficients, reference$coefficients)
  # Independently formed divisors avoid relying on agreement of two API calls.
  expected <- .4 / 50000 + .3 / 50000 + .6 / 2500000000
  expect_equal(d$coefficients$RelativeErrorVariance[3L], expected, tolerance = 1e-14)
  expect_true(all(is.finite(d$coefficients$G)))
})

test_that("D-studies retain interaction-specific averaging and composite covariances", {
  f <- mvgt_fixture()
  g <- mfrm_g_study(f$data, c("Content", "Organization"))
  grid <- data.frame(Raters = c(2L, 4L), Tasks = c(3L, 3L))
  w <- c(Organization = .4, Content = .6)
  d <- mfrm_d_study(g, grid, w)
  expect_s3_class(d, "mfrm_multivariate_d_study")
  expect_equal(nrow(d$coefficients), 6)
  expect_true(all(d$coefficients$Status == "Available"))
  expect_true(all(d$coefficients$Phi <= d$coefficients$G))
  # Hand-computed variances at 2 raters x 3 tasks.
  # Composite quadratic variances include BOTH symmetric covariance entries.
  u <- .6^2 * 1.2 + .4^2 * .8 + 2 * .6 * .4 * .3
  pr <- .6^2 * .4 + .4^2 * .5 + 2 * .6 * .4 * .1
  pt <- .6^2 * .3 + .4^2 * .4 - 2 * .6 * .4 * .06
  re <- .6^2 * .6 + .4^2 * .7 + 2 * .6 * .4 * .2
  r <- .6^2 * .2 + .4^2 * .3 + 2 * .6 * .4 * .04
  t <- .6^2 * .3 + .4^2 * .15 - 2 * .6 * .4 * .02
  rt <- .6^2 * .15 + .4^2 * .1 + 2 * .6 * .4 * .03
  relative <- pr / 2 + pt / 3 + re / 6
  absolute <- relative + r / 2 + t / 3 + rt / 6
  first <- d$coefficients[3L, ]
  expect_equal(first$UniverseVariance, u)
  expect_equal(first$RelativeErrorVariance, relative)
  expect_equal(first$AbsoluteErrorVariance, absolute)
  expect_equal(first$G, u / (u + relative))
  expect_equal(first$Phi, u / (u + absolute))
  expect_equal(first$RelativeSEM, sqrt(relative))
  expect_equal(first$AbsoluteSEM, sqrt(absolute))
  expect_equal(d$covariances[[1L]]$RelativeError,
    f$components[["Person:Rater"]] / 2 + f$components[["Person:Task"]] / 3 + f$components$Residual / 6)
  expect_identical(d$weights, w[g$design$scores])
  expect_identical(d$gstudy, g)
  expect_identical(summary(d), d$coefficients)
  expect_output(print(d), "no confidence intervals")
  baseline <- mfrm_d_study(g)
  expect_equal(baseline$design_grid, data.frame(Raters = 3L, Tasks = 4L))
  expect_true(all(baseline$coefficients$Kind == "Score"))
  # Shared complete-case observations, no separate univariate averaging.
  opposite <- f$data
  opposite$Organization <- -opposite$Organization
  other <- mfrm_d_study(mfrm_g_study(opposite,
    c("Content", "Organization")), grid, w)$coefficients
  expect_equal(other$G[other$Kind == "Score"], d$coefficients$G[d$coefficients$Kind == "Score"])
  expect_true(all(abs(other$G[other$Kind == "Composite"] - d$coefficients$G[d$coefficients$Kind == "Composite"]) > .001))
  # A single column reduces to the full-interaction univariate model.
  one <- mfrm_g_study(f$data, "Content")
  one_d <- mfrm_d_study(one, grid)
  expect_equal(one_d$coefficients$G, d$coefficients$G[d$coefficients$Score == "Content"])
  expect_equal(one_d$coefficients$Phi, d$coefficients$Phi[d$coefficients$Score == "Content"])
  expect_equal(one_d$coefficients$RelativeErrorVariance[1], .4 / 2 + .3 / 3 + .6 / 6)
  # Perfectly duplicated scores give a rank-deficient but valid point projection.
  duplicated <- f$data
  duplicated$Organization <- duplicated$Content
  duplicate_g <- mfrm_g_study(duplicated, c("Content", "Organization"))
  expect_true(all(duplicate_g$component_diagnostics$Rank == 1L))
  duplicate_d <- mfrm_d_study(duplicate_g, grid,
    c(Content = .6, Organization = .4))$coefficients
  expect_equal(duplicate_d$G[duplicate_d$Kind == "Composite"], one_d$coefficients$G)
  expect_equal(duplicate_d$RelativeSEM[duplicate_d$Kind == "Composite"], one_d$coefficients$RelativeSEM)
})

test_that("score order, score units and composite units are respected", {
  f <- mvgt_fixture()
  scores <- c("Content", "Organization")
  weights <- c(Content = .6, Organization = .4)
  g <- mfrm_g_study(f$data, scores)
  d <- mfrm_d_study(g, weights = weights)$coefficients
  reversed <- mfrm_d_study(mfrm_g_study(f$data, rev(scores)), weights = weights)$coefficients
  expect_equal(reversed[c(2, 1, 3), c("G", "Phi", "RelativeSEM", "AbsoluteSEM")],
    d[c("G", "Phi", "RelativeSEM", "AbsoluteSEM")], ignore_attr = TRUE)
  factors <- c(Content = 1e-4, Organization = 1e3)
  changed <- f$data
  changed[scores] <- sweep(as.matrix(changed[scores]), 2, factors, "*")
  scaled_g <- mfrm_g_study(changed, scores)
  scaled <- mfrm_d_study(scaled_g, weights = weights / factors)$coefficients
  expect_identical(scaled_g$component_diagnostics$PositiveSemidefinite, g$component_diagnostics$PositiveSemidefinite)
  expect_equal(scaled$G, d$G)
  expect_equal(scaled$Phi, d$Phi)
  expect_equal(scaled$RelativeSEM[3], d$RelativeSEM[3])
  expect_equal(scaled$AbsoluteSEM[1:2], d$AbsoluteSEM[1:2] * unname(factors))
  double <- mfrm_d_study(g, weights = 2 * weights)$coefficients[3, ]
  expect_equal(double$G, d$G[3])
  expect_equal(double$Phi, d$Phi[3])
  expect_equal(double$AbsoluteSEM, 2 * d$AbsoluteSEM[3])
  zero_weight <- mfrm_d_study(g, weights = c(Content = 1, Organization = 0))$coefficients
  expect_equal(zero_weight$G[3], zero_weight$G[1])
  tiny <- mfrm_d_study(g, weights = 1e-100 * weights)$coefficients[3, ]
  expect_equal(tiny$G, d$G[3])
  expect_equal(tiny$AbsoluteSEM / 1e-100, d$AbsoluteSEM[3])
  # A third, linearly dependent score preserves the same composite meaning.
  third <- f$data
  third$Language <- .5 * third$Content + third$Organization
  three <- mfrm_d_study(mfrm_g_study(third,
    c(scores, "Language")), weights = c(Content = .2, Organization = .3, Language = .5))
  equivalent <- mfrm_d_study(g, weights = c(Content = .45, Organization = .8))
  expect_equal(three$coefficients[4, c("G", "Phi", "AbsoluteSEM")],
    equivalent$coefficients[3, c("G", "Phi", "AbsoluteSEM")], ignore_attr = TRUE)
  expect_true(all(three$component_diagnostics$Dimension == 3L))
})

test_that("signed weights agree with direct difference scores and respect sign reversal", {
  data <- mvgt_fixture()$data
  g <- mfrm_g_study(data, c("Content", "Organization"))
  w <- c(Content = 1, Organization = -1)
  difference <- mfrm_d_study(g, weights = w)$coefficients[3, ]
  data$Difference <- data$Content - data$Organization
  direct <- mfrm_d_study(mfrm_g_study(data, "Difference"))$coefficients
  measures <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
                "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  expect_equal(difference[measures], direct[measures], ignore_attr = TRUE)
  expect_equal(difference$UniverseVariance, 1.2 + .8 - 2 * .3)
  reversed <- mfrm_d_study(g, weights = -w)$coefficients[3, ]
  expect_equal(reversed[measures], difference[measures])
  scaled <- mfrm_d_study(g, weights = -2 * w)$coefficients[3, ]
  expect_equal(scaled[c("G", "Phi")], difference[c("G", "Phi")])
  expect_equal(scaled$AbsoluteSEM, 2 * difference$AbsoluteSEM)
  negative <- mfrm_d_study(g,
    weights = c(Content = -1, Organization = 0))$coefficients
  expect_equal(negative[3, measures], negative[1, measures], ignore_attr = TRUE)
  # A contrast of identical scores has no universe variance, not perfect reliability.
  data$Organization <- data$Content
  zero <- mfrm_d_study(mfrm_g_study(data,
    c("Content", "Organization")), weights = w)$coefficients[3, ]
  expect_equal(zero$UniverseVariance, 0)
  expect_true(is.na(zero$G) && is.na(zero$Phi))
  expect_identical(zero$GStatus, "Zero total variance")
  expect_identical(zero$PhiStatus, "Zero total variance")
  expect_equal(zero$RelativeSEM, 0)
  expect_equal(zero$AbsoluteSEM, 0)
})

test_that("supplied mGENOVA matrices reproduce Appendix F D-study values", {
  # Brennan (2001), Manual for mGENOVA Version 2.1, Table 14 and Appendix F,
  # pp. 79-81. This is a D-study-only fixture, NOT a G-study estimate.
  # The original design has common tasks but local raters. Its zero
  # rater-related cross-covariances give the same arithmetic for this balanced
  # projection. This does not validate raw-data estimation or local-rater support.
  scores <- c("Listening", "Writing")
  values <- list(Person = c(.324, .356, .356, .691),
    Rater = c(.012, 0, 0, .010), Task = c(.127, .039, .039, .025),
    `Person:Rater` = c(.014, 0, 0, .047),
    `Person:Task` = c(.393, .030, .030, .159),
    `Rater:Task` = c(.022, 0, 0, .008), Residual = c(.317, 0, 0, .218))
  g <- structure(list(
    components = lapply(values, matrix, nrow = 2, dimnames = list(scores, scores)),
    score_scale = setNames(c(1, 1), scores),
    design = list(scores = scores, calculation_version = 1L,
      counts = c(Person = 50L, Rater = 2L, Task = 6L))),
    class = "mfrm_multivariate_gstudy")
  d <- mfrm_d_study(g, weights = c(Listening = 1, Writing = -1))
  expect_true(all(d$coefficients$Status == "Available"))
  # Printed values have five decimal places. Upper triangles of the manual's
  # output matrices contain correlations; lower triangles contain covariances.
  expected <- rbind(
    c(.32400, .09892, .12792, .76611, .71695, .31451, .35765),
    c(.69100, .06817, .07800, .91021, .89857, .26109, .27928),
    c(.30300, .15708, .18292, .65858, .62356, .39634, .42769))
  measures <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
                "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  expect_equal(unname(round(as.matrix(d$coefficients[measures]), 5)), expected)
  expect_equal(unname(d$covariances[[1]]$Universe), matrix(values$Person, 2))
  expect_equal(unname(round(d$covariances[[1]]$RelativeError, 5)),
    matrix(c(.09892, .00500, .00500, .06817), 2))
  expect_equal(unname(round(d$covariances[[1]]$AbsoluteError, 5)),
    matrix(c(.12792, .01150, .01150, .07800), 2))
})

test_that("negative RT components match executed jGENOVA with explicit D-study conventions", {
  # Official jGENOVA 1.0 source, run unmodified on 2026-09-21: P=4, R=3, T=4.
  # Raw Algorithm columns for each score, their sum and their difference;
  # external output uses seven decimals. No Java dependency for package tests.
  f <- mvgt_fixture(rater_task = c(-.1, .03, .03, .1))
  g <- mfrm_g_study(f$data, c("Content", "Organization"))
  expected <- cbind(Content = c(1.2, .2, .3, .4, .3, -.1, .6),
    Organization = c(.8, .3, .15, .5, .4, .1, .7),
    Sum = c(2.6, .58, .41, 1.1, .58, .06, 1.7),
    Difference = c(1.4, .42, .49, .7, .82, -.06, .9))
  vectors <- cbind(c(1, 0), c(0, 1), c(1, 1), c(1, -1))
  actual <- vapply(seq_len(ncol(vectors)), function(j)
    vapply(g$components, function(a) drop(crossprod(vectors[, j], a %*% vectors[, j])),
      numeric(1)), numeric(7))
  expect_equal(unname(actual), unname(expected), tolerance = 5e-8)
  covariance <- (expected[, "Sum"] - expected[, "Content"] - expected[, "Organization"]) / 2
  expect_equal(unname(vapply(g$components, function(a) a[1, 2], numeric(1))), covariance,
    tolerance = 5e-8)
  # GENOVA 3.1 (executed 2026-10-23) independently matches these jGENOVA
  # references: both zero negative variances under either convention.
  # The EMS option also changes R and T; mfrmr itself retains all raw estimates.
  h <- g
  h$design$scores <- "Content"
  h$score_scale <- g$score_scale["Content"]
  h$components <- lapply(g$components, function(a) a["Content", "Content", drop = FALSE])
  grid <- data.frame(Raters = 2, Tasks = 6)
  raw <- mfrm_d_study(h, grid)$coefficients
  expect_equal(raw$Phi, .730964467005076, tolerance = 1e-12)
  h$components$`Rater:Task`[,] <- 0
  algorithm <- mfrm_d_study(h, grid)$coefficients
  h$components$Rater[,] <- .175
  h$components$Task[,] <- 4 / 15
  ems <- mfrm_d_study(h, grid)$coefficients
  metrics <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
    "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  # These rows are printed external D-study output, not mfrmr reference fits.
  expect_equal(unname(unlist(algorithm[metrics])),
    c(1.2, .3, .45, .8, .72727, .54772, .67082), tolerance = 5.1e-6)
  expect_equal(unname(unlist(ems[metrics])),
    c(1.2, .3, .43194, .8, .73532, .54772, .65722), tolerance = 5.1e-6)
  expect_equal(g$components, f$components, tolerance = 1e-12)
})

test_that("negative-component D-studies match executed multivariate mGENOVA", {
  # Supplied mGENOVA 2.1 executable, run on 2026-09-21 with the same raw
  # 48-row data as the jGENOVA benchmark. Values below are printed output;
  # neither Wine nor mGENOVA is a package/test dependency.
  g <- mfrm_g_study(
    mvgt_fixture(rater_task = c(-.1, .03, .03, .1))$data,
    c("Content", "Organization"))
  weights <- cbind(Equal = c(Content = .5, Organization = .5),
    Difference = c(Content = 1, Organization = -1))
  grid <- data.frame(Raters = 2, Tasks = 6)
  metrics <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
    "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  # DOPTIONS NEGATIVE retains the raw components, as mfrmr does.
  expected <- rbind(
    c(1.2, .3, .44167, .8, .73096, .54772, .66458),
    c(.8, .375, .55833, .68085, .58896, .61237, .74722),
    c(.65, .19708, .28792, .76734, .69303, .44394, .53658),
    c(1.4, .56167, .84833, .71368, .62268, .74944, .92105))
  raw <- mfrm_d_study(g, grid, weights)
  expect_equal(unname(as.matrix(raw$coefficients[metrics])), expected, tolerance = 5.1e-6)
  expect_false(any(raw$coefficients$ComponentPSD))
  # Match mGENOVA's default explicitly: zero negative diagonals before
  # forming composites; preserve covariances. This is not a mfrmr option.
  g$components$`Rater:Task`[1, 1] <- 0
  expected[c(1, 3, 4), c(3, 5, 7)] <- rbind(
    c(.45, .72727, .67082), c(.29, .69149, .53852),
    c(.85667, .62038, .92556))
  adjusted <- mfrm_d_study(g, grid, weights)
  expect_equal(unname(as.matrix(adjusted$coefficients[metrics])), expected,
    tolerance = 5.1e-6)
})

test_that("component diagnostics do not suppress calculable score projections", {
  f <- mvgt_fixture(indefinite = TRUE)
  g <- mfrm_g_study(f$data, c("Content", "Organization"))
  expect_equal(g$components, f$components, tolerance = 1e-12)
  expect_false(g$component_diagnostics$PositiveSemidefinite[1])
  expect_true(all(diag(g$components$Person) > 0))
  d <- mfrm_d_study(g, weights = c(Content = .5, Organization = .5))
  expect_true(all(is.finite(as.matrix(d$coefficients[c("G", "Phi", "RelativeSEM", "AbsoluteSEM")]))))
  expect_true(all(d$coefficients$Status == "Available"))
  expect_false(any(d$coefficients$ComponentPSD))
  expect_output(print(d), "Warning: non-PSD covariance components")
  expect_equal(d$coefficients$UniverseVariance[1], .01)
  # Negative diagonal estimates are also kept, not replaced with zero.
  f <- mvgt_fixture()
  f$data$Content <- f$data$Content - ave(f$data$Content, f$data$Person, FUN = mean)
  negative <- mfrm_g_study(f$data, c("Content", "Organization"))
  expect_lt(negative$components$Person[1, 1], 0)
  neg_d <- mfrm_d_study(negative)$coefficients
  expect_true(is.na(neg_d$G[1]))
  expect_identical(neg_d$GStatus[1], "Negative universe variance")
  expect_true(is.finite(neg_d$G[2]))
  expect_true(all(is.finite(neg_d$RelativeSEM)))
  constant <- f$data
  constant$Content <- constant$Organization <- 7
  zero <- mfrm_g_study(constant, c("Content", "Organization"))
  expect_true(all(unlist(zero$components) == 0))
  expect_true(all(is.na(mfrm_d_study(zero)$coefficients$G)))
  zero_d <- mfrm_d_study(zero)$coefficients
  expect_true(all(zero_d$GStatus == "Zero total variance"))
  expect_true(all(zero_d$RelativeSEM == 0 & zero_d$AbsoluteSEM == 0))
})

test_that("each D-study metric uses only its own projected variances", {
  g <- mfrm_g_study(mvgt_fixture()$data, "Content")
  grid <- data.frame(Raters = 2, Tasks = 6)
  base <- mfrm_d_study(g, grid)$coefficients
  # Rater:Task affects absolute error only, even when its estimate is negative.
  g$components$`Rater:Task`[,] <- -12
  d <- mfrm_d_study(g, grid)$coefficients
  expect_equal(d$G, base$G)
  expect_equal(d$RelativeSEM, base$RelativeSEM)
  expect_lt(d$AbsoluteErrorVariance, 0)
  expect_true(is.na(d$Phi) && is.na(d$AbsoluteSEM))
  expect_identical(d$GStatus, "Available")
  expect_identical(d$PhiStatus, "Negative projected error variance")
  expect_identical(d$Status, "Partially available")
  expect_false(d$ComponentPSD)
  # Unconstrained estimates can reverse the usual G/Phi ordering while both
  # projected errors remain positive. Retain the values and the component flag.
  g$components$`Rater:Task`[,] <- -2
  d <- mfrm_d_study(g, grid)$coefficients
  expect_equal(d$G, .8)
  expect_equal(d$Phi, 72 / 89)
  expect_gt(d$Phi, d$G)
  expect_identical(d$Status, "Available")
  expect_false(d$ComponentPSD)
  # Absolute error can be positive when relative error is negative.
  g$components$`Rater:Task`[,] <- .15
  g$components$`Person:Rater`[,] <- -.4
  g$components$Rater[,] <- .4
  d <- mfrm_d_study(g, grid)$coefficients
  expect_equal(d$RelativeErrorVariance, -.1)
  expect_equal(d$AbsoluteErrorVariance, .1625)
  expect_true(is.na(d$G) && is.na(d$RelativeSEM))
  expect_equal(d$Phi, 1.2 / (1.2 + .1625))
  expect_equal(d$AbsoluteSEM, sqrt(.1625))
  expect_identical(d$GStatus, "Negative projected error variance")
  expect_identical(d$PhiStatus, "Available")
  # Zero universe variance gives zero reliability when total variance is positive.
  g$components$`Person:Rater`[,] <- .4
  g$components$Person[,] <- 0
  d <- mfrm_d_study(g, grid)$coefficients
  expect_equal(c(d$G, d$Phi), c(0, 0))
  expect_identical(d$Status, "Available")
})

test_that("incomplete, repeated or malformed observed designs are refused", {
  data <- mvgt_fixture()$data
  scores <- c("Content", "Organization")
  expect_error(mfrm_g_study(data[-1, ], scores), "complete, balanced")
  expect_error(mfrm_g_study(rbind(data, data[1, ]), scores), "Duplicate")
  expect_error(mfrm_g_study(subset(data, Rater == "R1"), scores), "at least two")
  for (value in list(NA_real_, NaN, Inf, "unknown")) {
    bad <- data
    bad$Content[1] <- value
    expect_error(mfrm_g_study(bad, scores), "finite numeric")
  }
  for (value in list(NA_character_, " ")) {
    bad <- data
    bad$Rater <- as.character(bad$Rater)
    bad$Rater[1] <- value
    expect_error(mfrm_g_study(bad, scores), "nonmissing, nonblank")
  }
  bad <- data
  bad$Content <- factor(bad$Content)
  expect_error(mfrm_g_study(bad, scores), "finite numeric")
  expect_error(mfrm_g_study(data, c("Content", "Content")), "distinct")
  expect_error(mfrm_g_study(data, "Rater"), "distinct")
  expect_error(mfrm_g_study(data, "Absent"), "distinct")
  expect_error(mfrm_g_study(data, scores, person = c("Person", "Task")), "distinct")
  expect_error(mfrm_g_study(data, scores, rater = "Person"), "distinct")
  extra <- data
  extra$Unused <- NA_real_
  expect_identical(mfrm_g_study(extra, scores), mfrm_g_study(data, scores))
  bad <- data
  bad$Content <- bad$Content * 1e200
  expect_error(mfrm_g_study(bad, scores), "rescale")
  bad$Content <- data$Content * 1e-200
  expect_error(mfrm_g_study(bad, scores), "rescale")
})

test_that("named composites match separate projections and direct weighted scores", {
  data <- mvgt_fixture()$data
  scores <- c("Content", "Organization")
  facets <- c(Rater = "Rater", Occasion = "Task")
  g <- mfrm_g_study(data, scores, facets = facets)
  grid <- expand.grid(Rater = c(2, 4), Occasion = c(3, 6))
  weights <- cbind(Equal = c(Organization = .5, Content = .5),
    ContentFocus = c(Organization = .2, Content = .8),
    Difference = c(Organization = -1, Content = 1),
    Scaled = c(Organization = -2, Content = 2))
  d <- mfrm_d_study(g, grid, weights)
  expect_identical(d$weights, weights[scores, , drop = FALSE])
  expect_equal(nrow(d$coefficients), nrow(grid) * (length(scores) + ncol(weights)))
  expect_identical(d$coefficients$Score, rep(c(scores, colnames(weights)), nrow(grid)))
  baseline <- mfrm_d_study(g, grid)
  expect_equal(d$covariances, baseline$covariances)
  expect_equal(d$coefficients[d$coefficients$Kind == "Score", ],
    baseline$coefficients, ignore_attr = TRUE)
  measures <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
                "G", "Phi", "RelativeSEM", "AbsoluteSEM", "Status")
  for (name in colnames(weights)) {
    separate <- mfrm_d_study(g, grid, weights[, name])$coefficients
    actual <- d$coefficients[d$coefficients$Kind == "Composite" & d$coefficients$Score == name, ]
    expect_equal(actual[measures], separate[separate$Kind == "Composite", measures],
      ignore_attr = TRUE)
  }
  # Estimate one weighted score directly from ratings, independently of the projection.
  data$Weighted <- .8 * data$Content + .2 * data$Organization
  direct <- mfrm_d_study(
    mfrm_g_study(data, "Weighted", facets = facets), grid)$coefficients
  expect_equal(d$coefficients[d$coefficients$Score == "ContentFocus", measures],
    direct[measures], ignore_attr = TRUE)
  single <- mfrm_g_study(data, "Content", facets = facets)
  one <- mfrm_d_study(single, grid,
    matrix(2, 1, 1, dimnames = list("Content", "Double")))$coefficients
  expect_equal(one$G[one$Kind == "Composite"], one$G[one$Kind == "Score"])
  expect_equal(one$AbsoluteSEM[one$Kind == "Composite"], 2 * one$AbsoluteSEM[one$Kind == "Score"])

  invalid <- mfrm_g_study(mvgt_fixture(indefinite = TRUE)$data, scores)
  projected <- mfrm_d_study(invalid, weights = weights)$coefficients
  difference <- projected$Score %in% c("Difference", "Scaled")
  expect_false(any(projected$ComponentPSD))
  expect_true(all(is.na(projected$G[difference])))
  expect_true(all(projected$GStatus[difference] == "Negative universe variance"))
  expect_true(all(is.finite(projected$G[!difference])))
  expect_true(all(is.finite(projected$RelativeSEM)))
})

test_that("weight matrices require explicit identities and usable columns", {
  g <- mfrm_g_study(mvgt_fixture()$data, c("Content", "Organization"))
  w <- cbind(Equal = c(Content = .5, Organization = .5), Difference = c(Content = 1, Organization = -1))
  bad_names <- list(list(NULL, colnames(w)), list(rownames(w), NULL),
    list(c("Content", "Absent"), colnames(w)),
    list(c("Content", "Content"), colnames(w)),
    list(c("Content", NA_character_), colnames(w)),
    list(rownames(w), c("Equal", "Equal")), list(rownames(w), c("Equal", " ")),
    list(rownames(w), c("Equal", NA_character_)))
  for (nm in bad_names) {
    bad <- w
    dimnames(bad) <- nm
    expect_error(mfrm_d_study(g, weights = bad), "weights.*matrix must")
  }
  for (bad in list(w[1, , drop = FALSE], w[, FALSE, drop = FALSE],
      cbind(w, Zero = 0), w * Inf, w * NA_real_, w + 1i, matrix(as.character(w), 2))) {
    expect_error(mfrm_d_study(g, weights = bad), "weights.*matrix must")
  }
  for (scale in c(1e200, 1e-200)) {
    bad <- w
    bad[, 2] <- scale * bad[, 2]
    expect_error(mfrm_d_study(g, weights = bad), "rescale")
  }
})

test_that("D-study input checks preserve the declared covariance identities", {
  g <- mfrm_g_study(mvgt_fixture()$data, c("Content", "Organization"))
  for (grid in list(data.frame(Raters = 0, Tasks = 3), data.frame(Raters = 2.5, Tasks = 3),
      data.frame(Raters = NA_real_, Tasks = 3), data.frame(Raters = Inf, Tasks = 3),
      data.frame(Raters = "2", Tasks = 3), data.frame(Raters = 2),
      data.frame(Raters = 2, Tasks = 3, Persons = 100), data.frame(Raters = numeric(), Tasks = numeric()))) {
    expect_error(mfrm_d_study(g, grid), "positive integer|exactly one count|missing count|at least one design row")
  }
  for (w in list(c(.6, .4), c(Content = .6), c(Content = .6, Absent = .4),
      c(Content = .6, Content = .4),
      c(Content = 0, Organization = 0), c(Content = Inf, Organization = 1),
      c(Content = NA_real_, Organization = 1))) {
    expect_error(mfrm_d_study(g, weights = w), "name every score")
  }
  bad <- g
  bad$components$Person[1, 2] <- 9
  expect_error(mfrm_d_study(bad), "covariance matrices")
  bad <- g
  bad$components$Person <- bad$components$Person[2:1, 2:1]
  expect_error(mfrm_d_study(bad), "score identities")
  bad <- g
  bad$components$Task <- NULL
  expect_error(mfrm_d_study(bad), "incomplete")
  expect_error(mfrm_d_study(list()), "must be a G-study")
  expect_error(mfrm_d_study(g, weights = c(Content = 1e308, Organization = 1e308)), "rescale")
  expect_error(mfrm_d_study(g, weights = c(Content = 1e-200, Organization = 1e-200)), "rescale")
})

test_that("MINQUE0 reduces to balanced ANOVA for one and two facets", {
  for (rater in list(NULL, "Rater")) {
    data <- if (is.null(rater)) mvgt_common_task_data() else mvgt_fixture()$data
    scores <- names(data)[(ncol(data) - 1):ncol(data)]
    anova <- mfrm_g_study(data, scores, rater = rater)
    minque <- mfrm_g_study(data, scores, rater = rater, method = "minque0")
    expect_equal(minque$components, anova$components, tolerance = 1e-11)
    expect_equal(mfrm_d_study(minque)$coefficients,
      mfrm_d_study(anova)$coefficients, tolerance = 1e-11)
    expect_null(minque$mean_products)
    expect_null(minque$degrees_of_freedom)
    expect_identical(minque$design$calculation_version, 2L)
  }
})

test_that("single-score MINQUE0 preserves complete and incomplete G/D-study calculations", {
  measures <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
    "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  for (rater in list(NULL, "Rater")) {
    full <- if (is.null(rater)) mvgt_common_task_data() else mvgt_fixture()$data
    scores <- tail(names(full), 2)
    for (data in list(full, full[seq_len(nrow(full)) %% 7 != 0, ])) {
      joint <- mfrm_g_study(data, scores, rater = rater, method = "minque0")
      one <- mfrm_g_study(data, scores[1L], rater = rater, method = "minque0")
      expect_equal(one$components,
        lapply(joint$components, function(a) a[1L, 1L, drop = FALSE]), tolerance = 1e-11)
      if (nrow(data) == nrow(full)) {
        anova <- mfrm_g_study(data, scores[1L], rater = rater)
        expect_equal(one$components, anova$components, tolerance = 1e-11)
      }
      grid <- if (is.null(rater)) data.frame(Tasks = 5L) else data.frame(Raters = 2L, Tasks = 3L)
      individual <- mfrm_d_study(one, grid)$coefficients
      combined <- mfrm_d_study(joint, grid)$coefficients
      expect_equal(individual[measures], combined[combined$Score == scores[1L], measures],
        tolerance = 1e-11, ignore_attr = TRUE)
    }
  }
})

test_that("incomplete multivariate MINQUE0 matches independent dense kernel calculations", {
  for (rater in list(NULL, "Rater")) {
    data <- if (is.null(rater)) mvgt_common_task_data() else mvgt_fixture()$data
    data <- data[seq_len(nrow(data)) %% 7 != 0, ]
    scores <- tail(names(data), 2)
    g <- mfrm_g_study(data, scores, rater = rater, method = "minque0")
    ids <- if (is.null(rater)) c("Person", "Task") else c("Person", "Rater", "Task")
    subsets <- if (is.null(rater)) list(1, 2, 1:2) else list(1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
    n <- nrow(data)
    h <- diag(n) - matrix(1/n, n, n)
    # Dense pairwise equality kernels do not use the production group-count trace code.
    kernels <- lapply(subsets, function(s) {
      Reduce(`*`, lapply(data[ids[s]], function(x) outer(x, x, `==`) * 1))
    })
    b <- lapply(kernels, function(k) h %*% k %*% h)
    dense_gram <- outer(seq_along(b), seq_along(b), Vectorize(function(i, j) sum(b[[i]] * b[[j]])))
    y <- as.matrix(data[scores])
    dense_q <- t(vapply(b, function(k) as.vector(t(y) %*% k %*% y), numeric(4)))
    reference <- solve(dense_gram, dense_q)
    expect_equal(unname(g$estimation$kernel_gram), dense_gram, tolerance = 1e-11)
    expect_equal(unname(t(vapply(g$components, as.vector, numeric(4)))), reference, tolerance = 1e-10)
    # Expected quadratic products under a specified covariance model must
    # recover its components, including off-diagonals, without simulation.
    truth <- seq_along(b) / 10
    covariance <- Reduce(`+`, Map(`*`, kernels, truth))
    expected_q <- vapply(b, function(k) sum(k * covariance), numeric(1))
    expect_equal(as.vector(solve(g$estimation$kernel_gram, expected_q)), truth, tolerance = 1e-11)
    expect_false(g$design$complete)
    expect_equal(g$design$observed_fraction, nrow(data) / g$design$potential_cells)
    expect_error(mfrm_d_study(g), "explicit.*design_grid")
    grid <- if (is.null(rater)) data.frame(Tasks = 5L) else data.frame(Raters = 2L, Tasks = 3L)
    expect_s3_class(mfrm_d_study(g, grid), "mfrm_multivariate_d_study")
    reversed <- mfrm_g_study(data[rev(seq_len(nrow(data))), ], rev(scores),
      rater = rater, method = "minque0")
    expect_equal(reversed$components, lapply(g$components, function(m) m[2:1, 2:1]), tolerance = 1e-10)
    # Covariances must also transform linearly when score units change.
    rescaled <- data
    rescaled[scores] <- sweep(as.matrix(data[scores]), 2, c(.01, 100), `*`)
    other <- mfrm_g_study(rescaled, scores, rater = rater, method = "minque0")
    expect_equal(other$components, lapply(g$components, function(m) m * outer(c(.01, 100), c(.01, 100))),
      tolerance = 1e-10)
  }
})

test_that("sparse covariance identification is separate from graph connectedness", {
  # A completely aliased pattern (one task unique to each person).
  diagonal <- data.frame(Person = 1:8, Task = 1:8, V = 1:8)
  expect_error(mfrm_g_study(diagonal, "V", rater = NULL, method = "minque0"),
    "Cannot separate covariance components")
  # In a connected two-facet design, one rater assigned per Person/Task
  # still aliases Person:Task with Residual; connectedness cannot fix this.
  crossed <- expand.grid(Person = 1:12, Task = 1:6)
  crossed$Rater <- (crossed$Person + crossed$Task) %% 4
  crossed$V <- sin(seq_len(nrow(crossed)))
  expect_error(mfrm_g_study(crossed, "V", method = "minque0"),
    "Cannot separate covariance components")
  # Negative estimates are retained, even with full-rank sparse moment equations.
  data <- mvgt_common_task_data()
  thin <- data[(data$Person + data$Task) %% 3 != 0, ]
  g <- mfrm_g_study(thin, c("V", "W"), rater = NULL, method = "minque0")
  expect_identical(g$estimation$rank, 3L)
  expect_false(all(g$component_diagnostics$PositiveSemidefinite))
  d <- mfrm_d_study(g, data.frame(Tasks = 6), c(V = -1, W = 1))
  expect_false(any(d$coefficients$ComponentPSD))
  expect_true(all(is.finite(d$coefficients$G[1:2])))
  expect_equal(d$coefficients$G[1:2], with(d$coefficients[1:2, ],
    UniverseVariance / (UniverseVariance + RelativeErrorVariance)))
  expect_lt(d$coefficients$UniverseVariance[3], 0)
  expect_true(is.na(d$coefficients$G[3]))
  expect_identical(d$coefficients$GStatus[3], "Negative universe variance")
})

test_that("missing-score omission retains one shared multivariate sample and exclusions", {
  data <- mvgt_common_task_data()
  data$V[2] <- NA_real_
  data$W[7] <- NA_real_
  data$Person[13] <- NA_integer_
  expect_error(mfrm_g_study(data, c("V", "W"), rater = NULL, method = "minque0"),
    "missing = 'omit'")
  g <- mfrm_g_study(data, c("V", "W"), rater = NULL, method = "minque0", missing = "omit")
  explicit <- mfrm_g_study(data[-c(2, 7, 13), ], c("V", "W"), rater = NULL, method = "minque0")
  expect_equal(g$components, explicit$components)
  expect_equal(g$data_usage$counts, c(InputRows = 60, UsedRows = 57, ExcludedRows = 3))
  expect_identical(g$data_usage$excluded_rows, c(2L, 7L, 13L))
  expect_equal(g$data_usage$missing_cells,
    data.frame(InputRow = c(13L, 2L, 7L), Column = c("Person", "V", "W")))
  expect_error(mfrm_g_study(data, c("V", "W"), rater = NULL, missing = "omit"),
    "ANOVA requires")
  bad <- data
  bad$V[1] <- Inf
  expect_error(mfrm_g_study(bad, c("V", "W"), rater = NULL,
    method = "minque0", missing = "omit"), "finite numeric")
  duplicate <- rbind(data, data[2, ])
  expect_error(mfrm_g_study(duplicate, c("V", "W"), rater = NULL,
    method = "minque0", missing = "omit"), "Duplicate")
  data$V[] <- NA_real_
  expect_error(mfrm_g_study(data, c("V", "W"), rater = NULL,
    method = "minque0", missing = "omit"), "No complete rows")
})

test_that("literal facet labels cannot merge different interaction cells", {
  data <- mvgt_fixture()$data
  renamed <- data
  renamed$Person <- c("a.b", "a", "x", "z")[match(data$Person, unique(data$Person))]
  renamed$Rater <- c("c", "b.c", "r")[match(data$Rater, unique(data$Rater))]
  for (method in c("anova", "minque0")) {
    original <- mfrm_g_study(data, c("Content", "Organization"), method = method)
    other <- mfrm_g_study(renamed, c("Content", "Organization"), method = method)
    expect_equal(other$components, original$components, tolerance = 1e-11)
  }
})

test_that("rater-only and named one-facet models preserve the reference calculations", {
  data <- mvgt_common_task_data()
  data$Assessor <- data$Task
  for (method in c("anova", "minque0")) {
    reference <- mfrm_g_study(data, c("V", "W"), rater = NULL, method = method)
    g <- mfrm_g_study(data, c("V", "W"), rater = "Assessor", task = NULL,
      method = method)
    expect_identical(g$design$counts, c(Person = 10L, Rater = 6L))
    expect_identical(names(g$components), c("Person", "Rater", "Residual"))
    expect_equal(unname(g$components), unname(reference$components))
    expect_output(print(g), "10 persons, 6 raters")
    d <- mfrm_d_study(g, weights = c(V = -1, W = 1))
    expect_identical(d$design_grid, data.frame(Raters = 6L))
    expect_output(print(d), "Means over common random raters;")
    expect_equal(d$coefficients[-2L],
      mfrm_d_study(reference, weights = c(V = -1, W = 1))$coefficients[-2L])
    # Display labels need not be syntactic R identifiers.
    named <- mfrm_g_study(data, c("V", "W"),
      facets = c("Test occasion" = "Assessor"), method = method)
    expect_identical(named$design$facets, c("Test occasion" = "Assessor"))
    expect_equal(unname(named$components), unname(reference$components))
    named_d <- mfrm_d_study(named)
    expect_identical(names(named_d$design_grid), "Test occasion")
    expect_true("Test occasion" %in% names(named_d$coefficients))
    expect_identical(named_d$coefficients$G, d$coefficients$G[1:2])
  }
})

test_that("named facets preserve two-facet estimates and projections after reordering", {
  data <- mvgt_fixture()$data
  data$Assessor <- data$Rater
  data$Session <- data$Task
  scores <- c("Content", "Organization")
  for (method in c("anova", "minque0")) {
    if (method == "minque0") data <- data[-c(1, 12, 29), ]
    reference <- mfrm_g_study(data, scores, method = method)
    g <- mfrm_g_study(data, scores,
      facets = c(Occasion = "Session", Rater = "Assessor"), method = method)
    expect_equal(unname(g$components), unname(reference$components[c(1, 3, 2, 5, 4, 6, 7)]),
      tolerance = 1e-10)
    d <- mfrm_d_study(g, data.frame(Rater = c(2, 4), Occasion = c(3, 6)),
      c(Content = .6, Organization = .4))
    expected <- mfrm_d_study(reference,
      data.frame(Raters = c(2, 4), Tasks = c(3, 6)), c(Content = .6, Organization = .4))
    expect_identical(names(d$design_grid), c("Occasion", "Rater"))
    expect_equal(d$coefficients[-c(2, 3)], expected$coefficients[-c(2, 3)], tolerance = 1e-10)
    expect_equal(d$covariances, expected$covariances, tolerance = 1e-10)
    if (method == "minque0") expect_error(mfrm_d_study(g), "explicit.*design_grid")
  }
  # Existing saved results without the new count mapping remain readable.
  old <- mfrm_g_study(mvgt_fixture()$data, scores)
  expected <- mfrm_d_study(old)$coefficients
  old$design$count_columns <- old$design$facets <- NULL
  expect_equal(mfrm_d_study(old)$coefficients, expected)
})

test_that("facet specifications refuse unsupported models and conflicting identities", {
  data <- mvgt_fixture()$data
  scores <- c("Content", "Organization")
  for (facets in list(character(), c("Rater", "Task", "Person"), c("Rater", "Rater"),
      c(Occasion = "Absent"), c(Occasion = NA_character_), setNames("Task", ""),
      c(Score = "Task"), c(Residual = "Task"), c(Panel = "Task"), c("A:B" = "Task"),
      c(Occasion = "Rater", Occasion = "Task"))) {
    expect_error(mfrm_g_study(data, scores, facets = facets))
  }
  expect_error(mfrm_g_study(data, scores, facets = "Task", rater = NULL), "either")
  expect_error(mfrm_g_study(data, scores, facets = "Rater", task = NULL), "either")
  # Selecting a facet never averages away repetitions from an omitted facet.
  expect_error(mfrm_g_study(data, scores, facets = "Rater"), "Duplicate")
  g <- mfrm_g_study(data, scores, facets = c("Rater", "Task"))
  expect_identical(names(g$design$counts), c("Person", "Rater", "Task"))
  expect_error(mfrm_d_study(g, data.frame(Raters = 2, Tasks = 3)), "missing count.*Rater, Task")
  g$design$count_columns[1] <- "Score"
  expect_error(mfrm_d_study(g), "count identities")
})
