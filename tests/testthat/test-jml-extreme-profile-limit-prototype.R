jml_profile_limit_prototype_path <- function() {
  testthat::test_path(
    "..", "..", "inst", "validation",
    "jml-extreme-profile-limit-prototype-0.2.3.R"
  )
}

load_jml_profile_limit_prototype <- function() {
  prototype <- jml_profile_limit_prototype_path()
  skip_if_not(file.exists(prototype),
              "repository-internal validation artifacts are excluded")
  env <- new.env(parent = globalenv())
  sys.source(prototype, envir = env)
  env
}

make_jml_profile_limit_data <- function() {
  data <- expand.grid(
    Person = sprintf("P%02d", seq_len(16L)),
    Rater = paste0("R", seq_len(3L)),
    Criterion = paste0("C", seq_len(3L)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  person <- as.integer(sub("P", "", data$Person))
  rater <- as.integer(sub("R", "", data$Rater))
  criterion <- as.integer(sub("C", "", data$Criterion))
  data$Score <- ((person + 2L * rater + criterion) %% 4L) + 1L
  data
}

fit_jml_profile_limit_fixture <- function(model = "RSM") {
  data <- if (identical(model, "RSM")) {
    make_jml_profile_limit_data()
  } else {
    load_mfrmr_data("example_core")
  }
  persons <- unique(as.character(data$Person))[1:2]
  data$Score[as.character(data$Person) == persons[1]] <- max(data$Score)
  data$Score[as.character(data$Person) == persons[2]] <- min(data$Score)
  args <- list(
    data = data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = min(data$Score),
    rating_max = max(data$Score),
    keep_original = TRUE,
    method = "JML",
    model = model,
    maxit = 240L,
    reltol = 1e-10
  )
  if (!identical(model, "RSM")) args$step_facet <- "Criterion"
  if (identical(model, "GPCM")) args$slope_facet <- "Criterion"
  fit <- suppressMessages(suppressWarnings(do.call(fit_mfrm, args)))
  list(fit = fit, persons = persons)
}

test_that("profile-limit refit recovers the extended JML supremum", {
  env <- load_jml_profile_limit_prototype()
  for (model in c("RSM", "PCM", "GPCM")) {
    fixture <- fit_jml_profile_limit_fixture(model)
    result <- suppressWarnings(env$mfrmr_jml_profile_limit_refit(
      fixture$fit,
      maxit = 240L,
      reltol = 1e-10,
      optimizer = "BFGS"
    ))

    expect_identical(result$State, "profile_limit_refit_verified")
    expect_true(isTRUE(result$Complete))
    expect_false(result$FiniteOriginalJMLMaximum)
    expect_false(result$OriginalLikelihoodMaximumAttained)
    expect_setequal(result$ExcludedPersons, fixture$persons)
    expect_identical(
      unname(result$ExcludedDirections[fixture$persons]),
      c("high", "low"),
      info = model
    )
    expect_gte(result$ProfileLimitGain, -1e-8)
    expect_true(result$LimitPathMonotone)
    expect_true(result$LimitGapNonnegative)
    expect_true(result$LimitVerified)
    expect_lte(result$TerminalLimitGap, 1e-8)
    expect_true(all(is.finite(result$parameters$RawFiniteEstimate)))
    expect_true(all(is.finite(result$parameters$ProfileLimitEstimate)))
    expect_true(all(diff(result$path$OriginalJMLLogLikelihood) >= -1e-8))
    expect_identical(result$ReadinessEffect, "none_prototype_only")
  }
})

test_that("profile-limit prototype does not relabel anchored extremes", {
  env <- load_jml_profile_limit_prototype()
  data <- make_jml_profile_limit_data()
  data$Score[data$Person == "P01"] <- 4L
  anchors <- data.frame(
    Facet = "Person", Level = "P01", Anchor = 0.5,
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    keep_original = TRUE,
    method = "JML",
    model = "RSM",
    anchors = anchors,
    maxit = 160L,
    reltol = 1e-10
  ))
  result <- env$mfrmr_jml_profile_limit_refit(fit)
  p01 <- fit$facets$person[fit$facets$person$Person == "P01", , drop = FALSE]

  expect_identical(p01$ParameterStatus, "fixed")
  expect_equal(p01$Estimate, 0.5)
  expect_identical(result$State, "no_free_extreme_persons")
  expect_length(result$ExcludedPersons, 0L)
  expect_true(is.na(result$FiniteOriginalJMLMaximum))
})

test_that("profile-limit prototype fails closed for coupled extremes", {
  env <- load_jml_profile_limit_prototype()
  data <- make_jml_profile_limit_data()
  data$Score[data$Person == "P01"] <- 4L
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    keep_original = TRUE,
    method = "JML",
    model = "RSM",
    noncenter_facet = "Rater",
    maxit = 160L,
    reltol = 1e-10
  ))
  result <- env$mfrmr_jml_profile_limit_refit(fit)

  expect_identical(
    result$State, "constraint_coupled_extreme_not_profiled"
  )
  expect_false(result$Complete)
  expect_identical(result$ConstraintCoupledPersons, "P01")
  expect_identical(result$ReadinessEffect, "none_prototype_only")
})

test_that("profile-limit prototype rejects non-JML fits and invalid caps", {
  env <- load_jml_profile_limit_prototype()
  data <- make_jml_profile_limit_data()
  mml <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    keep_original = TRUE,
    method = "MML",
    model = "RSM",
    maxit = 40L
  ))
  expect_error(
    env$mfrmr_jml_profile_limit_refit(mml),
    "requires a fitted JML"
  )

  fixture <- fit_jml_profile_limit_fixture("RSM")
  expect_error(
    env$mfrmr_jml_profile_limit_refit(fixture$fit, caps = c(4, -1)),
    "positive finite"
  )
})
