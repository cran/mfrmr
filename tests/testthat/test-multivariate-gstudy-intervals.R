mvgt_urgenova_interval_data <- function() {
  # Brennan, urGENOVA 2.1 manual, pp. 12 and 26 (Appendix C), cc.manual.
  # Ten persons, three raters, four items nested within each rater.
  y <- matrix(c(
    5,6,5,5,5,3,4,5,6,7,3,3,
    9,3,7,7,7,5,5,5,7,7,5,2,
    3,4,3,3,5,3,3,5,6,5,1,6,
    7,5,5,3,3,1,4,3,5,3,3,5,
    9,2,9,7,7,7,3,7,2,7,5,3,
    3,4,3,5,3,3,6,3,4,5,1,2,
    7,3,7,7,7,5,5,7,5,5,5,4,
    5,8,5,7,7,5,5,4,3,2,1,1,
    9,9,8,8,6,6,6,5,5,8,1,1,
    4,4,4,3,3,5,6,5,5,7,1,1), nrow = 10, byrow = TRUE)
  d <- expand.grid(Person = 1:10, Item = 1:4, Rater = 1:3)
  d$Score <- as.vector(y)
  d
}

test_that("narrow component intervals retain estimates without false zero-width certainty", {
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed, "Content")
  for (method in c("mls", "residual_contrast")) {
    ordinary <- confint(g, parm = "Residual", method = method, assumption = "normal")
    narrow <- confint(g, parm = "Residual", method = method, level = 1e-18, assumption = "normal")
    expect_identical(narrow$Estimate, ordinary$Estimate)
    expect_identical(narrow$SE, ordinary$SE)
    expect_true(is.na(narrow$Lower) && is.na(narrow$Upper))
    expect_match(narrow$Status, "below numeric resolution")
    if (method == "mls") {
      expect_true(is.finite(narrow$RawLower))
      expect_equal(narrow$RawLower, narrow$RawUpper)
    } else {
      expect_false(narrow$Available)
      expect_identical(narrow$IntervalCenter, ordinary$IntervalCenter)
    }
  }
})

test_that("raw urGENOVA data reproduce published SEVC and Ting intervals", {
  g <- mfrm_g_study(mvgt_urgenova_interval_data(), "Score",
    facets = c(Rater = "Rater", Item = "Item"), nesting = c(Item = "Rater"))
  ci <- confint(g, level = .80, assumption = "normal")
  expect_identical(ci$Source, c("Person", "Rater", "Item(Rater)", "Person:Rater", "Residual"))
  expect_equal(round(ci$Estimate, 5), c(.47315, .32515, .64753, .55957, 2.38025))
  expect_equal(round(ci$SE, 5), c(.38558, .43799, .37941, .37663, .36949))
  expect_equal(round(ci$Lower, 4), c(.0606, 0, .2987, .1746, 1.9738))
  expect_equal(round(ci$Upper, 4), c(1.4657, 5.4144, 1.6724, 1.3202, 2.9581))
  expect_lt(ci$RawLower[2], 0)
  expect_identical(ci$Status[2], "Includes zero")
  expect_identical(ci$Method, c(rep("Ting MLS", 4), "Chi-square"))
  expect_identical(attr(ci, "source_design")$nesting, c(Item = "Rater"))
  expect_identical(unname(attr(ci, "source_rows")), c(120L, 120L, 0L))
  # The opposite input order must describe the same nested model.
  reverse <- mfrm_g_study(mvgt_urgenova_interval_data(), "Score",
    facets = c(Item = "Item", Rater = "Rater"), nesting = c(Item = "Rater"))
  expect_equal(confint(reverse, level = .80, assumption = "normal")$SE, ci$SE)
})

test_that("crossed component SEs follow independent mean-square sampling variances", {
  fixture <- new.env()
  # First expression defines the independent tensor-contrast fixture.
  eval(parse(test_path("test-multivariate-gtheory.R"))[[1L]], fixture)
  f <- fixture$mvgt_fixture()
  g <- mfrm_g_study(f$data, c("Content", "Organization"))
  ci <- confint(g, score = "Content", assumption = "normal")
  expect_equal(ci$Estimate, unname(vapply(f$components, function(a) a[1, 1], numeric(1))))
  # From the specified population components, not the fitted ANOVA map.
  ms <- c(12*1.2 + 4*.4 + 3*.3 + .6, 16*.2 + 4*.4 + 4*.15 + .6,
    12*.3 + 3*.3 + 4*.15 + .6, 4*.4 + .6, 3*.3 + .6, 4*.15 + .6, .6)
  v <- 2 * ms^2 / (c(3, 2, 3, 6, 9, 6, 18) + 2)
  expected <- c(sum(v[c(1,4,5,7)])/12^2, sum(v[c(2,4,6,7)])/16^2,
    sum(v[c(3,5,6,7)])/12^2, sum(v[c(4,7)])/4^2,
    sum(v[c(5,7)])/3^2, sum(v[c(6,7)])/4^2, v[7])
  expect_equal(ci$SE^2, expected, tolerance = 1e-12)
  expect_true(all(is.finite(ci$RawLower) & is.finite(ci$RawUpper)))
  expect_equal(ci$Lower[7], .6 * 18/qchisq(.975, 18))
  expect_equal(ci$Upper[7], .6 * 18/qchisq(.025, 18))
  more <- confint(g, assumption = "normal", level = .99)
  expect_true(all(more$RawLower <= ci$RawLower))
  expect_true(all(more$RawUpper >= ci$RawUpper))
})

test_that("one-facet composites match directly weighted scores and preserve units", {
  path <- test_path("..", "..", "inst", "extdata", "mgenova-table12.csv")
  if (!file.exists(path)) path <- system.file("extdata", "mgenova-table12.csv", package = "mfrmr")
  d <- read.csv(path)
  g <- mfrm_g_study(d, c("V", "W"), rater = NULL)
  numeric_columns <- c("Estimate", "SE", "Lower", "Upper", "RawLower", "RawUpper")
  for (w in list(c(V = .5, W = .5), c(W = -.3, V = 2))) {
    d$Weighted <- d$V * w["V"] + d$W * w["W"]
    direct <- mfrm_g_study(d, "Weighted", rater = NULL)
    ci <- confint(g, weights = w, assumption = "normal")
    expected <- confint(direct, assumption = "normal")
    for (column in numeric_columns) expect_equal(ci[[column]], expected[[column]], tolerance = 1e-10)
    expect_equal(attr(ci, "weights"), w[c("V", "W")])
    doubled <- confint(g, weights = 2*w, assumption = "normal")
    for (column in numeric_columns) expect_equal(doubled[[column]], 4*ci[[column]], tolerance = 1e-10)
    ms <- vapply(direct$mean_products, function(a) a[1,1], numeric(1))
    v <- 2 * ms^2 / (c(9,5,45) + 2)
    expect_equal(ci$SE^2, c((v[1]+v[3])/36, (v[2]+v[3])/100, v[3]), ignore_attr = TRUE)
  }
  d[c("V", "W")] <- d[c("V", "W")] * 1e-80
  tiny <- confint(mfrm_g_study(d, c("V", "W"), rater = NULL), assumption = "normal")
  original <- confint(g, assumption = "normal")
  for (column in numeric_columns) expect_equal(tiny[[column]] / 1e-160, original[[column]], tolerance = 1e-10)
  expect_identical(tiny$Status, original$Status)
  expect_equal(confint(g, parm = c("Residual", "Person"), assumption = "normal")$SE, original$SE[c(3,1)])
  expect_equal(confint(g, parm = c(3L, 1L), assumption = "normal")$SE, original$SE[c(3,1)])
  # A selected score's precision must not disappear because an unused score
  # is measured on a vastly larger scale.
  d <- read.csv(path)
  d$V <- d$V * 1e80; d$W <- d$W * 1e-80
  unequal_units <- confint(mfrm_g_study(d, c("V", "W"), rater = NULL),
    score = "W", assumption = "normal")
  reference <- confint(g, score = "W", assumption = "normal")
  for (column in numeric_columns)
    expect_equal(unequal_units[[column]] / 1e-160, reference[[column]], tolerance = 1e-10)
  expect_identical(unequal_units$Status, reference$Status)
})

test_that("negative and degenerate estimates do not become false zero-width intervals", {
  d <- expand.grid(Person = 1:3, Task = 1:3)
  d$Score <- as.vector(outer(c(-1, 0, 1), c(-1, 0, 1)))
  g <- mfrm_g_study(d, "Score", rater = NULL)
  ci <- confint(g, assumption = "normal")
  expect_true(all(ci$Estimate[1:2] < 0))
  expect_true(all(ci$RawUpper[1:2] < 0))
  expect_true(all(is.na(ci$Lower[1:2]) & is.na(ci$Upper[1:2])))
  expect_identical(ci$Status[1:2], rep("Empty nonnegative interval", 2))
  expect_equal(ci$Estimate, unname(vapply(g$components, function(a) a[1,1], numeric(1))))
  d$Score <- 5
  zero <- confint(mfrm_g_study(d, "Score", rater = NULL), assumption = "normal")
  expect_equal(zero$SE, rep(0, 3))
  expect_true(all(is.na(zero$Lower) & is.na(zero$Upper)))
  expect_identical(zero$Status, rep("Degenerate mean squares", 3))
})

test_that("saved ANOVA summaries suffice and unsupported inputs are refused", {
  d <- mvgt_urgenova_interval_data()
  g <- mfrm_g_study(d, "Score", facets = c(Rater = "Rater", Item = "Item"),
    nesting = c(Item = "Rater"))
  original <- g
  expected <- confint(g, assumption = "normal")
  expect_identical(g, original)
  g$data <- NULL
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(g, file)
  expect_identical(confint(readRDS(file), assumption = "normal"), expected)
  saveRDS(expected, file)
  expect_identical(readRDS(file), expected)
  expect_error(confint(g), "assumption")
  expect_error(confint(g, assumption = "robust"), "assumption")
  expect_error(confint(g, assumption = "normal", level = 1), "level")
  expect_error(confint(g, assumption = "normal", score = "Missing"), "score")
  expect_error(confint(g, assumption = "normal", parm = "Missing"), "parm")
  expect_error(confint(g, assumption = "normal", parm = c(1,1)), "parm")
  expect_error(confint(g, assumption = "normal", weights = c(Other = 1)), "weights")
  expect_error(confint(g, assumption = "normal", score = "Score", weights = c(Score = 1)), "either")
  expect_error(confint(g, assumption = "normal", unknown = TRUE), "must be empty")
  bad <- g; bad$design$method <- "Multivariate MINQUE(0)"
  expect_error(confint(bad, assumption = "normal"), "ANOVA")
  bad <- g; bad$design$balanced <- FALSE
  expect_error(confint(bad, assumption = "normal"), "ANOVA")
  bad <- g; bad$mean_products <- NULL
  expect_error(confint(bad, assumption = "normal"), "matrices")
  bad <- g; bad$degrees_of_freedom$DF[1] <- 8
  expect_error(confint(bad, assumption = "normal"), "degrees of freedom")
  bad <- g; bad$components[[1]][1,1] <- bad$components[[1]][1,1] + 1
  expect_error(confint(bad, assumption = "normal"), "disagree")
  bad <- g; bad$design$child_counts[1] <- 5
  expect_error(confint(bad, assumption = "normal"), "child counts")
})
