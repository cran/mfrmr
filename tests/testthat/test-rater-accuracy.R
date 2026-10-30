accuracy_data <- function() {
  withr::local_seed(813)
  d <- expand.grid(Judge = paste0("J", 1:12), Essay = paste0("E", 1:18),
    Domain = c("A", "B"), stringsAsFactors = FALSE)
  a <- seq(-1, 1.5, length.out = 12)
  b <- seq(-.9, .9, length.out = 18)
  p <- plogis(a[match(d$Judge, unique(d$Judge))] -
    b[match(d$Essay, unique(d$Essay))] - ifelse(d$Domain == "A", -.35, .35))
  d$Reference <- (match(d$Essay, unique(d$Essay)) + as.integer(d$Domain == "B")) %% 3
  y <- rbinom(nrow(d), 1, p)
  d$Rating <- ifelse(y == 1, d$Reference, (d$Reference + 1) %% 3)
  d
}

accuracy_fixture <- function(method = "JML") {
  key <- paste0("accuracy_", method)
  if (!exists(key, envir = .mfrmr_test_cache, inherits = FALSE)) {
    assign(key, fit_rater_accuracy(accuracy_data(), "Judge", c("Essay", "Domain"),
      "Rating", "Reference", method = method, maxit = 500, quad_points = 61),
      envir = .mfrmr_test_cache)
  }
  get(key, envir = .mfrmr_test_cache, inherits = FALSE)
}

accuracy_glm <- function(d) {
  d$Match <- as.integer(d$Rating == d$Reference)
  for (nm in c("Judge", "Essay", "Domain")) d[[nm]] <- factor(d[[nm]])
  stats::glm(Match ~ Judge + Essay + Domain, data = d, family = stats::binomial(),
    contrasts = list(Judge = "contr.sum", Essay = "contr.sum", Domain = "contr.sum"),
    control = stats::glm.control(epsilon = 1e-12))
}

test_that("accuracy JML agrees with an independently fitted additive logistic model", {
  x <- accuracy_fixture()
  d <- accuracy_data(); g <- accuracy_glm(d)
  effects <- stats::predict(g, type = "terms")
  # Sum-coded GLM terms plus its adjusted constant reconstruct the same scale.
  expected_accuracy <- tapply(effects[, "Judge"] + attr(effects, "constant"),
    d$Judge, mean)
  expect_equal(x$table$Accuracy, as.numeric(expected_accuracy[x$table$Rater]), tolerance = 2e-5)
  expect_equal(as.numeric(x$fit$summary$LogLik), as.numeric(stats::logLik(g)), tolerance = 1e-7)
  for (nm in c("Essay", "Domain")) {
    expected_difficulty <- tapply(-effects[, nm], d[[nm]], mean)
    z <- x$difficulty[x$difficulty$Facet == nm, ]
    expect_equal(z$Difficulty, as.numeric(expected_difficulty[z$Level]), tolerance = 2e-5)
    expect_equal(sum(z$Difficulty), 0, tolerance = 1e-10)
  }
  expect_equal(x$cases$Match, as.integer(d$Rating == d$Reference))
  expect_equal(x$table$ExactAgreement, x$table$Matches / x$table$RatingRows)
  expect_true(all(is.na(x$table$PosteriorSD)))
  expect_identical(x$settings$column_map$Input, c("Judge", "Essay", "Domain"))
  expect_identical(x$settings$column_map$Model, c("Person", "Facet1", "Facet2"))
  expect_identical(x$settings$estimate, "rater_JML")
  expect_identical(x$table$EstimateUse,
    x$fit$facets$person$EstimateUse[match(x$table$Rater, x$fit$facets$person$Person)])
})

test_that("an incomplete connected assignment keeps its likelihood and benchmark meanings", {
  d <- accuracy_data()
  d <- d[(match(d$Judge, unique(d$Judge)) + match(d$Essay, unique(d$Essay))) %% 3 != 0, ]
  d <- d[rev(seq_len(nrow(d))), ]
  # Original Person is a benchmark facet, never the model's accuracy object.
  names(d)[names(d) == "Essay"] <- "Person"
  x <- fit_rater_accuracy(d, "Judge", c("Person", "Domain"), "Rating", "Reference", maxit = 500)
  names(d)[names(d) == "Person"] <- "Essay"
  g <- accuracy_glm(d)
  expect_equal(as.numeric(x$fit$summary$LogLik), as.numeric(stats::logLik(g)), tolerance = 1e-7)
  expect_equal(sum(x$table$RatingRows), nrow(d))
  expect_setequal(x$difficulty$Facet, c("Person", "Domain"))
  expect_identical(x$table$Rater, unique(d$Judge))
})

test_that("missing pairs are audited and never filled or assigned zero accuracy", {
  d <- accuracy_data()
  d$Rating[d$Judge == "J1"] <- NA
  d$Reference[2] <- NA; d$Rating[3] <- NA; d$Reference[3] <- NA
  x <- fit_rater_accuracy(d, "Judge", c("Essay", "Domain"), "Rating", "Reference", maxit = 500)
  j <- match("J1", x$table$Rater)
  expect_identical(x$table$EstimateUse[j], "no_retained_ratings")
  expect_true(is.na(x$table$Accuracy[j]))
  expect_true(is.na(x$table$ExactAgreement[j]))
  expect_equal(x$table$RatingRows[j], 0)
  expect_true(x$cases$MissingReference[2])
  expect_true(x$cases$MissingReference[3] && x$cases$MissingScore[3])
  expect_true(all(is.na(x$cases$Match[!x$cases$Included])))
  expect_equal(sum(x$table$RatingRows), sum(complete.cases(d[c("Rating", "Reference")])))
  expect_equal(sum(x$table$ExcludedRows), sum(!x$cases$Included))
  expect_equal(nrow(x$cases), nrow(d))
})

test_that("reference conflicts and repetitions fail before any fitting", {
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"))
  d <- accuracy_data()
  bad <- d; bad$Reference[1] <- 99; bad$Rating[1] <- NA
  expect_error(fit_rater_accuracy(bad, "Judge", c("Essay", "Domain"), "Rating", "Reference"),
    "Conflicting reference")
  bad <- rbind(d, d[1, ]); bad$Rating[nrow(bad)] <- NA
  expect_error(fit_rater_accuracy(bad, "Judge", c("Essay", "Domain"), "Rating", "Reference"),
    "Repeated rater-by-benchmark")
  bad <- d; bad$Judge[1] <- NA
  expect_error(fit_rater_accuracy(bad, "Judge", c("Essay", "Domain"), "Rating", "Reference"),
    "Identifiers must")
  bad <- d; bad$Rating <- bad$Reference
  expect_error(fit_rater_accuracy(bad, "Judge", c("Essay", "Domain"), "Rating", "Reference"),
    "both matches and mismatches")
  bad <- d; bad$Reference <- factor(bad$Reference)
  expect_error(fit_rater_accuracy(bad, "Judge", c("Essay", "Domain"), "Rating", "Reference"),
    "both be numeric")
  expect_error(fit_rater_accuracy(d, "Judge", "Essay", "Rating", "Rating"), "distinct")
  expect_error(fit_rater_accuracy(d, "Judge", NA_character_, "Rating", "Reference"), "must name")
  bad <- d; bad$Rating[1] <- Inf
  expect_error(fit_rater_accuracy(bad, "Judge", c("Essay", "Domain"), "Rating", "Reference"),
    "must be finite")
})

test_that("numeric identifier labels cannot silently merge raters or benchmarks", {
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"))
  for (column in c("Judge", "Essay", "Domain")) {
    d <- accuracy_data()
    labels <- unique(d[[column]])
    numbers <- seq_along(labels)
    numbers[1:2] <- c(1, 1 + 2e-15)
    d[[column]] <- numbers[match(d[[column]], labels)]
    expect_equal(as.character(numbers[1]), as.character(numbers[2]))
    expect_error(fit_rater_accuracy(d, "Judge", c("Essay", "Domain"),
      "Rating", "Reference"), "identical text labels")
  }
})

test_that("category labels compare by value and numeric comparisons never round", {
  original <- accuracy_fixture()
  local_mocked_bindings(fit_mfrm = function(...) original$fit)
  d <- accuracy_data()
  d$Rating <- factor(d$Rating, levels = 2:0)
  d$Reference <- factor(d$Reference, levels = 0:2)
  x <- fit_rater_accuracy(d, "Judge", c("Essay", "Domain"), "Rating", "Reference")
  expect_identical(x$cases$Match, original$cases$Match)
  d <- accuracy_data(); d$Rating[1] <- d$Reference[1] + 1e-12
  x <- fit_rater_accuracy(d, "Judge", c("Essay", "Domain"), "Rating", "Reference")
  expect_identical(x$cases$Match[1], 0L)
  d$Rating <- as.character(d$Rating); d$Reference <- as.character(d$Reference)
  d$Rating[1] <- " "
  expect_error(fit_rater_accuracy(d, "Judge", c("Essay", "Domain"), "Rating", "Reference"), "Blank")
})

test_that("JML extreme raters retain infinite primary estimates and source review", {
  d <- accuracy_data()
  d$Rating[d$Judge == "J1"] <- d$Reference[d$Judge == "J1"]
  d$Rating[d$Judge == "J2"] <- (d$Reference[d$Judge == "J2"] + 1) %% 3
  x <- suppressWarnings(fit_rater_accuracy(d, "Judge", c("Essay", "Domain"),
    "Rating", "Reference", maxit = 500))
  j <- match(c("J1", "J2"), x$table$Rater)
  expect_identical(x$table$Accuracy[j], c(Inf, -Inf))
  expect_identical(x$table$Extreme[j], c("high", "low"))
  expect_true(all(!x$table$SourceInferenceReady))
  expect_equal(x$table$ExactAgreement[j], c(1, 0))
  expect_equal(sum(x$table$RatingRows), nrow(d))
})

test_that("MML rater EAPs agree with independent continuous posterior integration", {
  x <- accuracy_fixture("MML")
  pop <- x$settings$population
  # Read fitted coordinates, then integrate the binary likelihood independently.
  mu <- pop$mean; sigma <- pop$sd
  expect_true(is.finite(mu) && is.finite(sigma) && sigma > 0)
  d <- accuracy_data()
  difficulty <- rep(0, nrow(d))
  for (nm in c("Essay", "Domain")) {
    z <- x$difficulty[x$difficulty$Facet == nm, ]
    difficulty <- difficulty + z$Difficulty[match(d[[nm]], z$Level)]
  }
  moments <- t(vapply(x$table$Rater, function(id) {
    take <- d$Judge == id
    y <- as.integer(d$Rating[take] == d$Reference[take]); b <- difficulty[take]
    logdensity <- function(theta) vapply(theta, function(t) {
      sum(stats::dbinom(y, 1, plogis(t - b), log = TRUE)) + stats::dnorm(t, mu, sigma, log = TRUE)
    }, numeric(1))
    shift <- max(logdensity(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = 101)))
    moment <- function(k) stats::integrate(function(t) t^k * exp(logdensity(t) - shift),
      -Inf, Inf, rel.tol = 1e-9)$value
    denom <- moment(0); mean <- moment(1) / denom
    c(mean, sqrt(moment(2) / denom - mean^2))
  }, numeric(2)))
  expect_equal(x$table$Accuracy, moments[, 1], ignore_attr = TRUE, tolerance = 2e-4)
  expect_equal(x$table$PosteriorSD, moments[, 2], ignore_attr = TRUE, tolerance = 2e-4)
  expect_identical(x$settings$uncertainty, "posterior_SD_given_fitted_parameters")
  expect_identical(x$settings$estimate, "rater_EAP")
  expect_identical(environment(x$fit$population$formula), baseenv())
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(x, path)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"))
  expect_identical(readRDS(path), x)
})
