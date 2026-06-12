test_that("sparse linked simulation preserves linking metadata and planned missingness", {
  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 5,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "sparse_linked",
    sparse_controls = list(
      link_persons = 3,
      link_raters_per_person = 5,
      assignment_mode = "balanced",
      min_common_persons_per_rater_pair = 3
    )
  )

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$assignment, "sparse_linked")
  expect_true(isTRUE(spec$sparse_controls$active))
  expect_equal(spec$sparse_controls$link_persons, 3L)
  expect_equal(spec$sparse_controls$link_raters_per_person, 5L)

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 20260526)
  sparse <- attr(sim, "mfrm_sparse_design")
  sim_meta <- attr(sim, "mfrm_simulation_spec")
  truth <- attr(sim, "mfrm_truth")

  expect_true(is.list(sparse))
  expect_true(isTRUE(sparse$active))
  expect_true(is.data.frame(sparse$overview))
  expect_equal(sparse$overview$Rows[1], nrow(sim))
  expect_equal(nrow(sim), ((3L * 5L) + (17L * 2L)) * 2L)
  expect_lt(sparse$overview$DesignDensity[1], 1)
  expect_gt(sparse$overview$PlannedMissingRate[1], 0)
  expect_equal(sparse$overview$LinkPersons[1], 3L)
  expect_equal(sparse$overview$MinCommonPersonsPerRaterPair[1], 3L)
  expect_equal(sparse$overview$RaterPairsBelowTarget[1], 0L)
  expect_equal(nrow(sparse$rater_pair_links), choose(5, 2))
  expect_true(all(sparse$rater_pair_links$CommonPersons >= 3L))
  expect_true(all(sparse$rater_coverage$Persons > 0L))
  expect_true(all(sparse$person_assignment$RatersAssigned[sparse$person_assignment$LinkPerson] == 5L))
  expect_true(all(sparse$person_assignment$RatersAssigned[!sparse$person_assignment$LinkPerson] == 2L))
  expect_identical(sim_meta$assignment, "sparse_linked")
  expect_true(isTRUE(sim_meta$sparse_controls$active))
  expect_true(isTRUE(truth$design$sparse$active))
})

test_that("sparse linked simulation supports custom facet names through a simulation spec", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 1,
    facet_names = c(rater = "Judge", criterion = "Task"),
    assignment = "sparse_linked",
    sparse_controls = list(link_fraction = 0.25)
  )

  sim <- simulate_mfrm_data(
    sim_spec = spec,
    seed = 99
  )

  sparse <- attr(sim, "mfrm_sparse_design")
  expect_true(isTRUE(sparse$active))
  expect_equal(sparse$overview$LinkPersons[1], 3L)
  expect_equal(sparse$overview$LinkRatersPerPerson[1], 4L)
  expect_lt(sparse$overview$DesignDensity[1], 1)
  expect_true(all(c("Study", "Person", "Judge", "Task", "Score") %in% names(sim)))
})

test_that("sparse linked controls validate infeasible linking settings", {
  expect_error(
    build_mfrm_sim_spec(
      n_person = 10,
      n_rater = 4,
      n_criterion = 2,
      raters_per_person = 2,
      assignment = "sparse_linked",
      sparse_controls = list(link_persons = 11)
    ),
    "between 0 and `n_person`"
  )

  expect_error(
    build_mfrm_sim_spec(
      n_person = 10,
      n_rater = 4,
      n_criterion = 2,
      raters_per_person = 3,
      assignment = "sparse_linked",
      sparse_controls = list(link_persons = 2, link_raters_per_person = 2)
    ),
    "must be >= `raters_per_person`"
  )
})

test_that("sparse linked controls can scale with design overrides", {
  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 4,
    n_criterion = 2,
    raters_per_person = 1,
    assignment = "sparse_linked",
    sparse_controls = list(link_fraction = 0.25)
  )
  override <- getFromNamespace("simulation_override_spec_design", "mfrmr")

  smaller <- override(spec, n_person = 12, n_rater = 5)
  expect_equal(smaller$sparse_controls$link_persons, 3L)
  expect_equal(smaller$sparse_controls$link_raters_per_person, 5L)
  sim <- simulate_mfrm_data(sim_spec = smaller, seed = 20260526)
  sparse <- attr(sim, "mfrm_sparse_design")
  expect_equal(sparse$overview$LinkPersons[1], 3L)
  expect_equal(sparse$overview$LinkRatersPerPerson[1], 5L)
})

test_that("design and recovery evaluation retain sparse linked diagnostics", {
  sim_eval <- suppressWarnings(evaluate_mfrm_design(
    n_person = 12,
    n_rater = 4,
    n_criterion = 2,
    raters_per_person = 1,
    assignment = "sparse_linked",
    sparse_controls = list(link_persons = 2, link_raters_per_person = 4),
    reps = 1,
    maxit = 10,
    seed = 20260526,
    progress = FALSE
  ))

  expect_true(sim_eval$rep_overview$SparseDesignActive[1])
  expect_equal(sim_eval$rep_overview$LinkPersons[1], 2L)
  expect_gt(sim_eval$rep_overview$PlannedMissingRate[1], 0)
  s_eval <- summary(sim_eval)
  expect_true(all(s_eval$design_summary$SparseDesignActive))
  expect_true(all(is.finite(s_eval$design_summary$MeanDesignDensity)))
  expect_s3_class(s_eval$sparse_review, "data.frame")
  expect_equal(s_eval$sparse_review$ReviewUse[1], "design_diagnostic_not_recovery_gate")
  expect_equal(s_eval$sparse_review$ReviewRows[1], 0L)
  expect_true(isTRUE(sim_eval$ademp$data_generating_mechanism$sparse_controls$active))
  p_sparse <- plot(
    sim_eval,
    facet = "Rater",
    metric = "plannedmissingrate",
    x_var = "n_person",
    draw = FALSE
  )
  expect_equal(p_sparse$metric_col, "MeanPlannedMissingRate")
  expect_true(all(is.finite(p_sparse$data$y)))

  design_bundle <- build_summary_table_bundle(s_eval)
  expect_true("sparse_design" %in% names(design_bundle$tables))
  expect_true("sparse_review" %in% names(design_bundle$tables))
  expect_equal(
    design_bundle$table_index$Role[design_bundle$table_index$Table == "sparse_design"],
    "sparse_design_diagnostics"
  )
  expect_true(all(c("MeanPlannedMissingRate", "MaxZeroCommonRaterPairs") %in%
                    names(design_bundle$tables$sparse_design)))
  expect_true(all(c("LinkReviewStatus", "LinkReviewReason", "ReviewUse") %in%
                    names(design_bundle$tables$sparse_design)))
  expect_true(all(design_bundle$tables$sparse_design$LinkReviewStatus == "ok"))

  rec <- suppressWarnings(evaluate_mfrm_recovery(
    n_person = 12,
    n_rater = 4,
    n_criterion = 2,
    raters_per_person = 1,
    assignment = "sparse_linked",
    sparse_controls = list(link_persons = 2, link_raters_per_person = 4),
    reps = 1,
    maxit = 10,
    seed = 20260527
  ))

  expect_true(rec$rep_overview$SparseDesignActive[1])
  expect_equal(rec$rep_overview$LinkPersons[1], 2L)
  expect_gt(rec$rep_overview$PlannedMissingRate[1], 0)
  expect_true(isTRUE(rec$settings$sparse_controls$active))
  s_rec <- summary(rec)
  expect_s3_class(s_rec$sparse_review, "data.frame")
  expect_equal(s_rec$sparse_review$SparseDesignRows[1], 1L)
  recovery_bundle <- build_summary_table_bundle(s_rec)
  expect_true("sparse_design" %in% names(recovery_bundle$tables))
  expect_true("sparse_review" %in% names(recovery_bundle$tables))
  expect_equal(
    recovery_bundle$table_index$Role[recovery_bundle$table_index$Table == "sparse_design"],
    "recovery_sparse_design_diagnostics"
  )
  expect_true(all(c("PlannedMissingRate", "ZeroCommonRaterPairs") %in%
                    names(recovery_bundle$tables$sparse_design)))
  expect_identical(recovery_bundle$tables$sparse_design$ReviewUse[1],
                   "design_diagnostic_not_recovery_gate")
  recovery_compact <- build_summary_table_bundle(summary(rec), appendix_preset = "compact")
  expect_true("sparse_design" %in% names(recovery_compact$tables))
})

test_that("sparse linked bundle diagnostics flag weak rater-pair linkage without changing metric status", {
  sparse_tbl <- data.frame(
    Facet = c("Rater", "Rater"),
    SparseDesignActive = TRUE,
    MeanDesignDensity = c(0.42, 0.45),
    MeanPlannedMissingRate = c(0.58, 0.55),
    MeanMinCommonPersonsPerRaterPair = c(0, 2),
    MaxZeroCommonRaterPairs = c(1, 0),
    MaxRaterPairsBelowTarget = c(1, 1),
    TargetCommonPersonsPerRaterPair = c(2, 3),
    stringsAsFactors = FALSE
  )
  out <- mfrmr:::summary_table_bundle_sparse_design_df(sparse_tbl)

  expect_equal(out$LinkReviewStatus, c("review", "review"))
  expect_match(out$LinkReviewReason[1], "no common persons", fixed = TRUE)
  expect_match(out$LinkReviewReason[2], "below the requested common-person target", fixed = TRUE)
  expect_true(all(out$ReviewUse == "design_diagnostic_not_recovery_gate"))
})
