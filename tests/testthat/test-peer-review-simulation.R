test_that("peer-review simulation specs create fixed no-self-review skeletons", {
  spec <- build_peer_review_sim_spec(
    n_submission = 12,
    n_criterion = 3,
    reviewers_per_submission = 2,
    anchor_submissions = 2
  )

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$source, "peer_review")
  expect_equal(spec$assignment, "skeleton")
  expect_true(isTRUE(spec$peer_review$active))
  expect_identical(unname(spec$facet_names), c("Reviewer", "Criterion"))
  expect_true(all(c("TemplatePerson", "Rater", "Criterion", "TemplatePersonReuse") %in%
                    names(spec$design_skeleton)))
  expect_true(all(spec$design_skeleton$TemplatePersonReuse))
  expect_false(any(spec$design_skeleton$TemplatePerson == spec$design_skeleton$Rater))
  expect_equal(spec$peer_review$overview$AnchorSubmissions[1], 2L)
  expect_equal(spec$peer_review$overview$OrdinaryReviewersPerSubmission[1], 2L)
  expect_equal(spec$peer_review$overview$AnchorReviewersPerSubmission[1], 11L)

  no_anchor <- build_peer_review_sim_spec(
    n_submission = 8,
    n_criterion = 2,
    reviewers_per_submission = 2,
    anchor_fraction = 0
  )
  expect_equal(no_anchor$peer_review$overview$AnchorSubmissions[1], 0L)
  expect_equal(no_anchor$peer_review$overview$AnchorReviewersPerSubmission[1], 0L)
})

test_that("simulate_mfrm_data preserves peer-review assignment metadata", {
  spec <- build_peer_review_sim_spec(
    n_submission = 12,
    n_criterion = 2,
    reviewers_per_submission = 2,
    anchor_submissions = 2
  )

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 20260526)
  peer <- attr(sim, "mfrm_peer_review_design")
  sim_meta <- attr(sim, "mfrm_simulation_spec")
  truth <- attr(sim, "mfrm_truth")

  expect_true(is.list(peer))
  expect_true(isTRUE(peer$active))
  expect_true(is.data.frame(peer$overview))
  expect_equal(peer$overview$Scenario[1], "peer_review")
  expect_equal(peer$overview$Submissions[1], 12L)
  expect_equal(peer$overview$Reviewers[1], 12L)
  expect_equal(peer$overview$Criteria[1], 2L)
  expect_equal(peer$overview$SelfReviews[1], 0L)
  expect_true(peer$overview$DesignDensity[1] < 1)
  expect_true(nrow(peer$reviewer_pair_common_submissions) > 0L)
  expect_true(all(c("Reviewer1", "Reviewer2", "CommonSubmissions") %in%
                    names(peer$reviewer_pair_common_submissions)))
  expect_false(any(sim$Person == sim$Reviewer))
  expect_equal(sim_meta$assignment, "skeleton")
  expect_true(isTRUE(sim_meta$peer_review$active))
  expect_true(isTRUE(truth$design$peer_review$active))
})

test_that("build_peer_review_design_review exposes assignment diagnostics for reporting", {
  spec <- build_peer_review_sim_spec(
    n_submission = 12,
    n_criterion = 2,
    reviewers_per_submission = 2,
    anchor_submissions = 2
  )
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 20260526)

  review <- build_peer_review_design_review(sim, top_n = 4)
  expect_s3_class(review, "mfrm_bundle")
  expect_s3_class(review, "mfrm_peer_review_design_review")
  expect_true(all(c(
    "overview", "load_summary", "submission_load", "reviewer_load",
    "reviewer_pair_common_submissions", "low_common_pairs",
    "reciprocal_pairs", "reporting_map"
  ) %in% names(review)))
  expect_equal(review$overview$PeerReviewStatus[1], "ok")
  expect_equal(review$overview$ReviewUse[1], "design_diagnostic_not_measurement_gate")
  expect_equal(review$overview$SelfReviews[1], 0L)
  expect_true(nrow(review$low_common_pairs) <= 4L)
  expect_true(any(review$reporting_map$Area == "MFRM measurement model"))

  s <- summary(review, top_n = 2)
  expect_s3_class(s, "summary.mfrm_peer_review_design_review")
  expect_true(nrow(s$low_common_pairs) <= 2L)

  bundle <- build_summary_table_bundle(review)
  expect_identical(bundle$source_class, "mfrm_peer_review_design_review")
  expect_identical(bundle$summary_class, "summary.mfrm_peer_review_design_review")
  expect_true(all(c("overview", "load_summary", "low_common_pairs",
                    "reporting_map") %in% names(bundle$tables)))
  expect_identical(
    as.character(bundle$table_index$Role[bundle$table_index$Table == "overview"]),
    "peer_review_design_diagnostics"
  )
  expect_true("peer_review_low_common_links" %in% bundle$table_index$Role)

  attr_review <- build_peer_review_design_review(attr(sim, "mfrm_peer_review_design"))
  expect_s3_class(attr_review, "mfrm_peer_review_design_review")
})

test_that("peer-review simulation validates infeasible reviewer assignment settings", {
  expect_error(
    build_peer_review_sim_spec(
      n_submission = 6,
      reviewers_per_submission = 6
    ),
    "eligible reviewers"
  )

  expect_error(
    build_peer_review_sim_spec(
      n_submission = 6,
      reviewers_per_submission = 3,
      anchor_submissions = 2,
      anchor_reviewers_per_submission = 2
    ),
    "must be >= `reviewers_per_submission`"
  )
})
