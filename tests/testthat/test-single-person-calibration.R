single_person_ratings <- function() data.frame(
  Person = "P1", Rater = rep(c("A", "B"), each = 4),
  Score = c(1, 1, 1, 0, 1, 0, 0, 0))

single_person_args <- function(model = "RSM", method = "MML") {
  a <- list(data = single_person_ratings(), person = "Person", facets = "Rater",
    score = "Score", rating_min = 0, rating_max = 1, category_policy = "preserve",
    repeated_ratings = "independent", model = model, method = method,
    maxit = 300, reltol = 1e-10)
  if (model == "PCM") a$step_facet <- "Rater"
  if (method == "MML") a$quad_points <- 61L
  a
}

# Independent integral for ONE joint eight-response vector, not eight Persons
# or four independently integrated two-response patterns.
single_person_mml_reference <- function(b) {
  integrand <- function(theta, order) {
    a <- plogis(theta - b); z <- plogis(theta + b)
    probability <- a^3 * (1-a) * z * (1-z)^3
    score <- -2 + 4*a - 4*z
    multiplier <- switch(order, 1, score, score^2 - 4*a*(1-a) - 4*z*(1-z))
    probability * multiplier * dnorm(theta)
  }
  value <- vapply(1:3, function(k) integrate(integrand, -Inf, Inf,
    order = k, rel.tol = 1e-11, abs.tol = 1e-12)$value, numeric(1))
  c(NLL = -log(value[1]), Gradient = -value[2]/value[1],
    Information = (value[2]/value[1])^2 - value[3]/value[1])
}

# Given total 4, the count of A successes ranges from 0 to 4; there are
# choose(4,s)^2 response vectors for each count. No package ESF is used.
single_person_cml_reference <- function(b) {
  statistic <- 4 - 2*(0:4)
  eta <- 2*lchoose(4, 0:4) + b*statistic
  weights <- exp(eta - max(eta)); probability <- weights/sum(weights)
  expected <- sum(statistic*probability)
  c(NLL = 2*b + max(eta) + log(sum(weights)), Gradient = expected + 2,
    Information = sum((statistic-expected)^2*probability))
}

test_that("one-Person full RSM/PCM fits match independent likelihoods", {
  skip_if_not_installed("lpSolve")
  for (method in c("MML", "CML", "JML")) {
    reference <- switch(method, MML = single_person_mml_reference,
      CML = single_person_cml_reference, JML = NULL)
    root <- if (method == "JML") -log(3) else
      uniroot(function(b) reference(b)["Gradient"], c(-5, 5), tol = 1e-11)$root
    if (!is.null(reference)) {
      expected <- reference(root)
      expect_equal(unname((reference(root+1e-4)["NLL"]-reference(root-1e-4)["NLL"])/2e-4),
        unname(expected["Gradient"]), tolerance = 1e-7)
      expect_equal(unname((reference(root+1e-4)["Gradient"]-reference(root-1e-4)["Gradient"])/2e-4),
        unname(expected["Information"]), tolerance = 1e-7)
    }
    for (model in c("RSM", "PCM")) {
      fit <- suppressWarnings(do.call(fit_mfrm, single_person_args(model, method)))
      expect_equal(fit$prep$n_person, 1L)
      expect_equal(nobs(fit, "ratings"), 8L)
      expect_equal(unname(fit$facets$others$Estimate), c(root, -root), tolerance = 1e-7)
      expect_true("single_person_calibration" %in% fit$prep$preparation_notes$Condition)
      if (method == "JML") {
        expect_equal(unname(coef(fit))[1], 0, tolerance = 1e-7)
        expect_error(vcov(fit), "JML point estimates")
        expect_error(confint(fit), "JML point estimates")
        expect_true(any(grepl("Only one Person", summary(fit)$preparation_notes$Message)))
      } else {
        expect_equal(fit$opt$value, unname(expected["NLL"]), tolerance = 1e-9)
        V <- suppressWarnings(vcov(fit))
        expect_equal(as.numeric(V), 1/unname(expected["Information"]), tolerance = 1e-7)
        expect_true(any(grepl("Only one Person", attr(V, "cautions"))))
        ci <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
        expect_true(all(ci$table$CIEligible))
        expect_false(any(ci$table$CoverageQualified))
        expect_true(all(grepl("Only one Person", ci$table$InferenceCaution)))
        expect_equal(ci$table$SE, rep(sqrt(as.numeric(V)), 2), tolerance = 1e-7)
        bounds <- suppressWarnings(confint(fit))
        expect_equal(unname(bounds[1, ]), root + c(-1, 1)*qnorm(.975)*sqrt(as.numeric(V)),
          tolerance = 1e-7)
        expect_true(any(grepl("Only one Person", attr(bounds, "cautions"))))
        expect_identical(unserialize(serialize(ci, NULL)), ci)
        if (method == "MML") {
          drawing <- build_wright_map_data(fit)
          expect_equal(diff(drawing$person_hist$breaks), 1)
          expect_equal(drawing$person_hist$counts, 1L)
          expect_equal(mean(drawing$person_hist$breaks), unname(fit$facets$person$Estimate))
          expect_match(drawing$retention_note, "not an uncertainty interval")
          facets_drawing <- build_wright_map_data(fit, wright_style = "facets_style")
          expect_false(grepl("histogram bin", facets_drawing$retention_note, fixed = TRUE))
          result <- mfrm_results(fit, intervals = list(raters = ci), compute = "never")
          expect_identical(result$tables$facet_raters_intervals$InferenceCaution, ci$table$InferenceCaution)
          expect_identical(mfrm_report(result)$tables$facet_raters_intervals$InferenceCaution,
            ci$table$InferenceCaution)
        } else expect_true(any(grepl("Only one Person", summary(fit)$notes)))
      }
    }
  }
})

test_that("polytomous one-Person CML agrees with independent response enumeration", {
  skip_if_not_installed("lpSolve")
  patterns <- as.matrix(expand.grid(rep(list(0:2), 6)))
  patterns <- patterns[rowSums(patterns) == 6, , drop = FALSE]
  for (model in c("RSM", "PCM")) {
    a <- single_person_args(model, "CML")
    a$data <- data.frame(Person = "P1", Rater = rep(c("A", "B"), each = 3),
      Score = rep(0:2, 2))
    a$rating_max <- 2L
    fit <- do.call(fit_mfrm, a)
    # Centered severity (b,-b); each centered step ladder is (t,-t).
    X <- cbind(rowSums(patterns[, 4:6])-rowSums(patterns[, 1:3]),
      -rowSums(patterns[, 1:3] == 1), -rowSums(patterns[, 4:6] == 1))
    observed <- c(0, -1, -1)
    if (model == "RSM") { X <- cbind(X[, 1], rowSums(X[, 2:3])); observed <- c(0, -2) }
    eta <- drop(X %*% coef(fit)); w <- exp(eta - max(eta)); p <- w/sum(w)
    mean <- drop(crossprod(p, X)); centered <- sweep(X, 2, mean)
    information <- crossprod(centered, centered * p)
    expect_equal(fit$opt$value, max(eta)+log(sum(w))-sum(observed*coef(fit)), tolerance = 1e-10)
    expect_lt(max(abs(mean-observed)), 1e-6)
    expect_equal(unname(as.matrix(vcov(fit))), unname(solve(information)),
      tolerance = 1e-7, ignore_attr = TRUE)
  }
})

test_that("one-Person permission does not waive model, repetition or inference checks", {
  skip_if_not_installed("lpSolve")
  a <- single_person_args(method = "CML")
  unresolved <- a; unresolved$repeated_ratings <- "review"
  expect_error(suppressWarnings(do.call(fit_mfrm, unresolved)), "resolution of repeated rating cells")
  separated <- a
  separated$data$Score <- rep(c(1, 0), each = 4)
  expect_error(do.call(fit_mfrm, separated), "recession_detected")
  for (method in c("JML", "MML")) {
    rejected <- single_person_args(method = method)
    rejected$data <- separated$data
    expect_error(suppressWarnings(do.call(fit_mfrm, rejected)), "recession_detected")
  }
  confounded <- a; confounded$data$Task <- confounded$data$Rater
  confounded$facets <- c("Rater", "Task")
  expect_error(do.call(fit_mfrm, confounded), "unidentified")
  a <- single_person_args()
  a$data <- a$data[rep(seq_len(nrow(a$data)), each = 2), ]
  fit <- do.call(fit_mfrm, a) # No singleton category counts; isolate cluster requirement.
  expect_error(mfrm_facet_intervals(fit, "Rater", method = "sandwich"), "at least two independent clusters")
  pop <- a; pop$population_formula <- ~ 1
  pop$person_data <- data.frame(Person = "P1")
  expect_error(do.call(fit_mfrm, pop), "selected route retains a software support limit")
  for (method in c("MML", "JML")) {
    gpcm <- single_person_args("GPCM", method)
    gpcm$step_facet <- gpcm$slope_facet <- "Rater"
    expect_error(do.call(fit_mfrm, gpcm), "selected route retains a software support limit")
  }
  gpcm$jml_correction_order <- 1L
  gpcm$reltol <- NULL
  expect_error(do.call(fit_mfrm, gpcm), "selected route retains a software support limit")
  comparison <- single_person_args(method = "CML"); comparison$fit_scope <- "comparisons"
  expect_error(do.call(fit_mfrm, comparison), "selected route retains a software support limit")
  product_data <- expand.grid(Person = "P1", Task = c("T1", "T2", "T3"),
    Rater = c("A", "B", "C"), stringsAsFactors = FALSE)
  product_data$Score <- rep(0:2, 3)
  expect_error(fit_mfrm(product_data, "Person", c("Task", "Rater"), "Score",
    model = "GPCM", method = "MML", mml_engine = "em", slope_facet = c("Task", "Rater"),
    step_facet = "Rater", noncenter_facet = "Rater",
    gpcm_mml_identification = "fixed_standard_normal", rating_min = 0, rating_max = 2,
    category_policy = "preserve"), "selected route retains a software support limit")
  review <- describe_mfrm_data(a$data, "Person", "Rater", "Score",
    rating_min = 0, rating_max = 1, repeated_ratings = "independent", include_agreement = FALSE)
  expect_true(any(grepl("Only one Person", review$preparation_notes$Message)))
})

test_that("an unresolved one-Person recession check cannot release optimizer points", {
  local_mocked_bindings(mfrmr_additive_response_boundary = function(...)
    list(state = "not_evaluated", reason = "Independent unresolved-check sentinel"), .package = "mfrmr")
  for (method in c("JML", "MML"))
    expect_error(suppressWarnings(do.call(fit_mfrm, single_person_args(method = method))),
      "not_evaluated.*Independent unresolved-check sentinel.*No finite point")
})
