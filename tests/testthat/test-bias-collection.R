bias_collection_fixture <- local({
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  dat <- mfrmr:::sample_mfrm_data(seed = 654)
  fit <- suppressWarnings(mfrmr::fit_mfrm(
    data = dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  diagnostics <- suppressWarnings(mfrmr::diagnose_mfrm(fit, residual_pca = "none"))
  bias_all <- suppressWarnings(mfrmr::estimate_all_bias(
    fit,
    diagnostics = diagnostics,
    max_iter = 2
  ))

  list(
    fit = fit,
    diagnostics = diagnostics,
    bias_all = bias_all
  )
})

test_that("estimate_all_bias batches all modeled facet pairs", {
  bias_all <- bias_collection_fixture$bias_all

  expect_s3_class(bias_all, "mfrm_bias_collection")
  expect_true(is.data.frame(bias_all$summary))
  expect_true(is.list(bias_all$by_pair))
  expect_true(all(c(
    "Interaction", "Rows", "ScreenPositive", "Significant",
    "FormalInferenceEligible", "PrimaryReportingEligible", "Kept"
  ) %in% names(bias_all$summary)))
  expect_identical(bias_all$summary$ScreenPositive, bias_all$summary$Significant)
  expect_true(all(!bias_all$summary$FormalInferenceEligible))
  expect_true(all(!bias_all$summary$PrimaryReportingEligible))
  expect_true(all(bias_all$summary$ReportingUse == "screening_only"))
  expect_identical(
    sort(bias_all$summary$Interaction),
    sort(c("Rater x Task", "Rater x Criterion", "Task x Criterion"))
  )
  expect_true(length(bias_all$by_pair) >= 1)

  one_summary <- summary(bias_all$by_pair[[1]])$overview
  expect_true(all(c(
    "ScreenPositive", "BonferroniScreenPositive", "HolmScreenPositive",
    "Significant", "BonferroniSignificant", "HolmSignificant",
    "FormalInferenceEligible", "PrimaryReportingEligible"
  ) %in% names(one_summary)))
  expect_identical(one_summary$ScreenPositive, one_summary$Significant)
  expect_identical(
    one_summary$BonferroniScreenPositive,
    one_summary$BonferroniSignificant
  )
  expect_identical(one_summary$HolmScreenPositive, one_summary$HolmSignificant)
  expect_false(one_summary$FormalInferenceEligible)
  expect_false(one_summary$PrimaryReportingEligible)
})

test_that("estimate_all_bias accepts explicit pair specifications", {
  fit <- bias_collection_fixture$fit
  diagnostics <- bias_collection_fixture$diagnostics

  bias_subset <- suppressWarnings(mfrmr::estimate_all_bias(
    fit,
    diagnostics = diagnostics,
    pairs = list(c("Rater", "Criterion")),
    max_iter = 2
  ))

  expect_s3_class(bias_subset, "mfrm_bias_collection")
  expect_identical(bias_subset$summary$Interaction, "Rater x Criterion")
  expect_true("Rater x Criterion" %in% names(bias_subset$by_pair) || !bias_subset$summary$Kept[1])
})

test_that("explicit Person pairs do not require automatic Person enumeration", {
  one_facet_fit <- list(config = list(facet_names = "Rater"))

  resolved <- mfrmr:::resolve_bias_collection_pairs(
    one_facet_fit,
    pairs = list(c("Person", "Rater")),
    include_person = FALSE
  )

  expect_identical(resolved[[1]]$facets, c("Person", "Rater"))
  expect_identical(resolved[[1]]$label, "Person x Rater")
})

test_that("estimate_bias supports explicit Person screens with provenance", {
  fit <- bias_collection_fixture$fit
  diagnostics <- bias_collection_fixture$diagnostics

  expect_warning(
    person_bias <- mfrmr::estimate_bias(
      fit,
      diagnostics = diagnostics,
      facet_a = "Person",
      facet_b = "Rater",
      max_iter = 1
    ),
    regexp = "conditional plug-in likelihood screen",
    class = "mfrmr_person_bias_screen_warning"
  )

  expect_s3_class(person_bias, "mfrm_bias")
  expect_identical(person_bias$interaction_facets, c("Person", "Rater"))
  expect_true(isTRUE(person_bias$person_interaction))
  expect_identical(person_bias$person_source_column, "Person")
  expect_match(person_bias$person_interaction_note, "not a jointly fitted Person x facet")
  expect_gt(nrow(person_bias$table), 0L)
  expect_true(all(person_bias$table$FacetA == "Person"))
  expect_true(all(!person_bias$table$SupportsFormalInference))
  expect_true(all(!person_bias$table$FormalInferenceEligible))

  orientation <- as.data.frame(person_bias$orientation_review)
  expect_identical(orientation$Sign[orientation$Facet == "Person"], 1)
  expect_identical(
    orientation$Orientation[orientation$Facet == "Person"],
    "positive"
  )
  expect_true(isTRUE(person_bias$mixed_sign))

  aliased_note <- mfrmr:::bias_person_screen_note(list(
    config = list(source_columns = list(person = "Candidate"))
  ))
  expect_match(aliased_note, "source person column `Candidate`", fixed = TRUE)
})

test_that("estimate_all_bias warns once and runs Person pairs in fit order", {
  fit <- bias_collection_fixture$fit
  diagnostics <- bias_collection_fixture$diagnostics

  person_warnings <- list()
  bias_with_person <- withCallingHandlers(
    mfrmr::estimate_all_bias(
      fit,
      diagnostics = diagnostics,
      include_person = TRUE,
      max_iter = 1
    ),
    mfrmr_person_bias_screen_warning = function(cnd) {
      person_warnings[[length(person_warnings) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_length(person_warnings, 1L)
  expect_match(
    conditionMessage(person_warnings[[1]]),
    "estimate_all_bias.*included.*Person"
  )

  expected <- c(
    "Person x Rater", "Person x Task", "Person x Criterion",
    "Rater x Task", "Rater x Criterion", "Task x Criterion"
  )
  expect_identical(bias_with_person$summary$Interaction, expected)
  expect_true(all(bias_with_person$summary$Rows > 0L))
  expect_true(all(bias_with_person$summary$Kept))
  expect_true(all(bias_with_person$summary$Error == ""))
  expect_identical(names(bias_with_person$by_pair), expected)
  expect_equal(nrow(bias_with_person$errors), 0L)
  expect_true(any(
    bias_with_person$settings$Setting == "person_pairs_requested" &
      bias_with_person$settings$Value == "TRUE"
  ))
})

test_that("MML Person screens use fitted EAP measures and posterior SEs", {
  dat <- mfrmr:::sample_mfrm_data(seed = 655)
  keep <- unique(dat$Person)[seq_len(12)]
  dat <- dat[dat$Person %in% keep, , drop = FALSE]
  fit <- suppressWarnings(mfrmr::fit_mfrm(
    dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    quad_points = 7,
    maxit = 20
  ))
  diagnostics <- suppressWarnings(mfrmr::diagnose_mfrm(
    fit,
    residual_pca = "none"
  ))

  expect_warning(
    person_bias <- mfrmr::estimate_bias(
      fit,
      diagnostics,
      interaction_facets = c("Person", "Rater"),
      max_iter = 1
    ),
    class = "mfrmr_person_bias_screen_warning"
  )

  person_measures <- diagnostics$measures[
    diagnostics$measures$Facet == "Person",
    c("Level", "Estimate", "SE", "SE_Method"),
    drop = FALSE
  ]
  matched <- match(person_bias$table$FacetA_Level, person_measures$Level)
  expect_false(anyNA(matched))
  expect_equal(
    person_bias$table$FacetA_Measure,
    person_measures$Estimate[matched],
    tolerance = 1e-10
  )
  expect_equal(
    person_bias$table$FacetA_SE,
    person_measures$SE[matched],
    tolerance = 1e-10
  )
  expect_true(all(grepl(
    "Posterior SD",
    person_measures$SE_Method[matched],
    fixed = TRUE
  )))
})
