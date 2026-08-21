gpcm_latent_distribution_stress_runner <- function() {
  validation_path <- function(file) {
    testthat::test_path("..", "..", "inst", "validation", file)
  }
  support <- validation_path("gpcm-estimator-asymptotics-0.2.3.R")
  runner <- validation_path("gpcm-latent-distribution-stress-0.2.3.R")
  skip_if_not(file.exists(support),
              "repository-internal validation artifacts are excluded")
  skip_if_not(file.exists(runner),
              "repository-internal validation artifacts are excluded")
  env <- new.env(parent = globalenv())
  sys.source(support, envir = env)
  sys.source(runner, envir = env)
  env
}

test_that("distribution stress fails clearly without its support runner", {
  path <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-latent-distribution-stress-0.2.3.R"
  )
  skip_if_not(file.exists(path),
              "repository-internal validation artifacts are excluded")
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  expect_error(
    env$mfrmr_glds_require_support(),
    "Source `gpcm-estimator-asymptotics-0.2.3.R`",
    fixed = TRUE
  )
})

test_that("latent supports isolate shape while preserving first moments", {
  env <- gpcm_latent_distribution_stress_runner()
  audit <- env$mfrmr_glds_support_audit()
  expect_setequal(
    audit$Distribution,
    c("normal", "skewed_gamma2", "symmetric_mixture")
  )
  expect_true(all(audit$MomentContractPassed))
  expect_equal(audit$Mean, rep(0, 3), tolerance = 1e-12)
  expect_equal(audit$Variance, rep(1, 3), tolerance = 1e-12)
  expect_lt(abs(audit$Skewness[audit$Distribution == "normal"]), 0.01)
  expect_gt(audit$Skewness[audit$Distribution == "skewed_gamma2"], 1)
  expect_lt(
    abs(audit$Skewness[audit$Distribution == "symmetric_mixture"]),
    0.01
  )
  expect_lt(
    audit$ExcessKurtosis[audit$Distribution == "symmetric_mixture"],
    -1
  )
})

test_that("stress manifest crosses distributions, exposure, and estimators", {
  env <- gpcm_latent_distribution_stress_runner()
  manifest <- env$mfrmr_glds_manifest()
  plan <- env$mfrmr_run_gpcm_latent_distribution_stress()
  expect_equal(nrow(manifest), 144L)
  expect_equal(length(unique(manifest$CellId)), 6L)
  expect_true(all(table(manifest$Replicate, manifest$CellId) == 2L))
  expect_setequal(manifest$Method, c("JML", "MML"))
  expect_setequal(manifest$ObservationsPerPerson, c(8L, 24L))
  expect_true(all(!manifest$EstimatorSelectionAuthorized))
  expect_true(all(!manifest$BayesianComparatorAuthorized))
  expect_true(all(!manifest$ConfirmationAuthorized))
  expect_true(all(!manifest$CapabilityPromotionAuthorized))
  expect_s3_class(plan, "mfrmr_gpcm_latent_distribution_stress_plan")
  expect_equal(nrow(plan$manifest), 144L)
})

test_that("matched distributions preserve structural truth and nested rows", {
  env <- gpcm_latent_distribution_stress_runner()
  manifest <- env$mfrmr_glds_manifest(replicates = 1L)
  bases <- env$mfrmr_glds_generate_replicate(manifest, replicate = 1L)
  audit <- env$mfrmr_glds_coupling_audit(
    bases, replicate = 1L, seed = unique(manifest$Seed)
  )
  expect_true(all(audit$CouplingContractPassed))
  expect_true(all(audit$DesignMatched))
  expect_true(all(audit$FacetTruthMatched))
  expect_true(all(audit$StepTruthMatched))
  expect_true(all(audit$SlopeTruthMatched))
  expect_gt(
    audit$RealizedAbilitySkewness[
      audit$Distribution == "skewed_gamma2"
    ],
    0.5
  )

  cells <- env$mfrmr_glds_cells()
  for (distribution in unique(cells$Distribution)) {
    low_cell <- cells[
      cells$Distribution == distribution &
        cells$ObservationsPerPerson == 8L, , drop = FALSE
    ]
    high_cell <- cells[
      cells$Distribution == distribution &
        cells$ObservationsPerPerson == 24L, , drop = FALSE
    ]
    low <- env$mfrmr_gas_subset_cell(bases[[distribution]], low_cell)
    high <- env$mfrmr_gas_subset_cell(bases[[distribution]], high_cell)
    key <- function(data) {
      paste(data$Person, data$Rater, data$Criterion, data$Score, sep = "\r")
    }
    expect_true(all(key(low) %in% key(high)))
    expect_equal(nrow(low), low_cell$ExpectedRows)
    expect_equal(nrow(high), high_cell$ExpectedRows)
  }
})

test_that("normal contrasts use paired replicate differences", {
  env <- gpcm_latent_distribution_stress_runner()
  metrics <- data.frame(
    Replicate = rep(1:2, 2),
    CellId = rep(c("normal-L08", "skewed_gamma2-L08"), each = 2),
    NPersons = 120L,
    ObservationsPerPerson = 8L,
    Method = "MML",
    Component = "step",
    Coordinates = 12L,
    RMSE = c(0.2, 0.3, 0.4, 0.5),
    MAE = c(0.1, 0.2, 0.3, 0.4),
    stringsAsFactors = FALSE
  )
  manifest <- data.frame(
    CellId = c("normal-L08", "skewed_gamma2-L08"),
    Distribution = c("normal", "skewed_gamma2"),
    ObservationsPerPerson = 8L,
    stringsAsFactors = FALSE
  )
  contrast <- env$mfrmr_glds_normal_contrasts(metrics, manifest)
  expect_equal(contrast$Replicates, 2L)
  expect_equal(contrast$MeanRMSEDifference, 0.2)
  expect_equal(contrast$StressHigherRMSERate, 1)
  expect_identical(contrast$ReferenceDistribution, "normal")
})

test_that("population target labels do not leak onto structural rows", {
  env <- gpcm_latent_distribution_stress_runner()
  row <- data.frame(Distribution = "skewed_gamma2")
  fixture <- data.frame(
    ParameterType = c("facet", "population"),
    stringsAsFactors = FALSE
  )
  decorated <- env$mfrmr_glds_decorate(fixture, row)
  expect_identical(
    decorated$PopulationTargetUse,
    c(
      "not_population_parameter",
      "generating_moment_difference_under_normal_misspecification"
    )
  )
  status <- env$mfrmr_glds_decorate(data.frame(FitReturned = TRUE), row)
  expect_identical(status$PopulationTargetUse, "not_applicable")
})

test_that("distribution stress has no machine-identity gate", {
  runner <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-latent-distribution-stress-0.2.3.R"
  )
  contract <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-latent-distribution-stress-contract-0.2.3.md"
  )
  skip_if_not(file.exists(runner),
              "repository-internal validation artifacts are excluded")
  skip_if_not(file.exists(contract),
              "repository-internal validation artifacts are excluded")
  text <- paste(
    c(readLines(runner, warn = FALSE), readLines(contract, warn = FALSE)),
    collapse = "\n"
  )
  expect_false(grepl("sha-?256|digest::|serialize\\(", text,
                     ignore.case = TRUE))
  expect_match(
    text,
    "DistributionRobustnessDecision = pilot_completed_no_prespecified_decision_rule",
    fixed = TRUE
  )
  expect_match(text, "BayesianComparatorAuthorized = FALSE", fixed = TRUE)
})

test_that("distribution stress record preserves its bounded conclusion", {
  record <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-latent-distribution-stress-record-0.2.3.md"
  )
  skip_if_not(file.exists(record),
              "repository-internal validation artifacts are excluded")
  text <- paste(readLines(record, warn = FALSE), collapse = "\n")
  expect_match(text, "all 144 fit calls returned objects", fixed = TRUE)
  expect_match(text, "sparse-exposure MML advantage survived", fixed = TRUE)
  expect_match(text, "fitted normal population variance was shape-", fixed = TRUE)
  expect_match(
    text,
    "DistributionRobustnessDecision = pilot_completed_no_prespecified_decision_rule",
    fixed = TRUE
  )
  expect_match(text, "BayesianComparatorAuthorized = FALSE", fixed = TRUE)
  expect_match(text, "CapabilityPromotionAuthorized = FALSE", fixed = TRUE)
  expect_false(grepl("SHA[^\n]*=|MD5[^\n]*=", text))
})
