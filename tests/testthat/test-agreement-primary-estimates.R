test_that("unavailable Person points cannot leak optimizer proxies into expected agreement", {
  fit <- make_toy_fit()
  obs <- compute_obs_table(fit)
  facets <- c("Person", fit$config$facet_names)
  original <- calc_interrater_agreement(obs, facets, "Rater", fit)
  # Presentation-state injection isolates the guard from the response backend.
  fit$facets$person$Estimate[1] <- Inf
  x <- calc_interrater_agreement(obs, facets, "Rater", fit)
  observed <- c("N", "Exact", "Adjacent", "Corr", "MeanDiff", "MAD", "ExactCount")
  expect_identical(x$pairs[observed], original$pairs[observed])
  expect_true(all(is.na(x$pairs$ExpectedExact)))
  expect_true(all(x$pairs$ExpectedReason == "primary_estimate_unavailable"))
  expect_equal(x$pairs$ExpectedAvailableContexts + x$pairs$ExpectedUnavailableContexts, x$pairs$N)
  expect_true(all(x$pairs$ExpectedUnavailableContexts > 0))
  expect_true(is.na(x$summary$ExpectedExactAgreement))
  expect_match(reporting_interrater_text(list(interrater = x)), "optimizer proxies")
  expect_match(reporting_interrater_text(list(interrater = x)), "EAP for MML")
  expect_identical(serialize(fit$opt$par, NULL), serialize(make_toy_fit()$opt$par, NULL))
})

test_that("nonfinite facets and PCM steps affect only their incident comparisons", {
  fit <- make_toy_fit(model = "PCM")
  obs <- compute_obs_table(fit)
  facets <- c("Person", fit$config$facet_names)
  original <- calc_interrater_agreement(obs, facets, "Rater", fit)
  changed <- fit
  at <- which(changed$facets$others$Facet == "Rater")[1]
  rater <- as.character(changed$facets$others$Level[at])
  changed$facets$others$Estimate[at] <- NA_real_
  x <- calc_interrater_agreement(obs, facets, "Rater", changed)
  affected <- x$pairs$Rater1 == rater | x$pairs$Rater2 == rater
  expect_true(all(is.na(x$pairs$ExpectedExact[affected])))
  expect_identical(x$pairs[!affected, ], original$pairs[!affected, ])
  changed <- fit
  changed$steps$Estimate[1] <- Inf
  x <- calc_interrater_agreement(obs, facets, "Rater", changed)
  expect_true(all(x$pairs$ExpectedUnavailableContexts > 0))
  expect_true(all(x$pairs$ExpectedAvailableContexts > 0))
  expect_true(all(is.na(x$pairs$ExpectedExact)))
  expect_true(all(x$pairs$ExpectedReason == "primary_estimate_unavailable"))
})

test_that("one-family slope availability and valid plug-in probabilities stay distinct", {
  fit <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
  obs <- compute_obs_table(fit)
  facets <- c("Person", fit$config$facet_names)
  original <- calc_interrater_agreement(obs, facets, "Rater", fit)
  expect_true(all(is.finite(original$pairs$ExpectedExact)))
  for (value in c(NA_real_, 0, Inf)) {
    changed <- fit; changed$slopes$Estimate[1] <- value
    x <- calc_interrater_agreement(obs, facets, "Rater", changed)
    expect_identical(x$pairs$Exact, original$pairs$Exact)
    expect_true(any(x$pairs$ExpectedReason == "primary_estimate_unavailable"))
  }
  changed <- fit; changed$steps$Estimate[1] <- NA_real_
  x <- calc_interrater_agreement(obs, facets, "Rater", changed)
  expect_true(all(x$pairs$ExpectedAvailableContexts > 0))
  expect_true(all(x$pairs$ExpectedUnavailableContexts > 0))
})

test_that("a native extreme-Person fit keeps observed agreement and explains missing expectations", {
  data <- load_mfrmr_data("example_core")
  id <- as.character(data$Person[1])
  data$Score[as.character(data$Person) == id] <- min(data$Score)
  fit <- suppressWarnings(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", maxit = 100))
  expect_identical(fit$facets$person$Estimate[match(id, fit$facets$person$Person)], -Inf)
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  x <- interrater_agreement_table(fit, diagnostics)
  expect_true(is.finite(x$summary$ExactAgreement))
  expect_true(is.na(x$summary$ExpectedExactAgreement))
  expect_match(x$summary$ExpectedReason, "primary_estimate_unavailable")
  expect_identical(x$pairs$ExpectedReason,
    diagnostics$interrater$pairs$ExpectedReason[match(x$pairs$Pair,
      paste(diagnostics$interrater$pairs$Rater1, diagnostics$interrater$pairs$Rater2, sep = " | "))])
  expect_match(reporting_interrater_text(diagnostics), "optimizer proxies")
  assign("agreement_native_extreme", list(fit = fit, diagnostics = diagnostics,
    result = x), envir = .mfrmr_test_cache)
})
