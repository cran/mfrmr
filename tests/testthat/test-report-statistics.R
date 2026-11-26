test_that("agreement, common-plan scores and the focused card reuse fitted results", {
  ratings <- load_mfrmr_data("example_core")
  facets <- c("Rater", "Criterion")
  wf <- make_core_workflow()
  plan <- data.frame(Rater = c("R01", "R02"), Criterion = c("Content", "Organization"))
  out <- summary(wf, reference_plan = plan)
  ag <- out$report_tables[["Expected agreement and model-based kappa"]]

  rows <- split(seq_len(nrow(ratings)), interaction(ratings$Person, ratings$Criterion))
  pairs <- do.call(cbind, lapply(rows, combn, m = 2L))
  probs <- predict(wf, type = "probs")
  observed <- mean(ratings$Score[pairs[1L, ]] == ratings$Score[pairs[2L, ]])
  expected <- mean(rowSums(probs[pairs[1L, ], ] * probs[pairs[2L, ], ]))
  expect_equal(ag$Pairs, ncol(pairs))
  expect_equal(ag$ExpectedExact, expected)
  expect_equal(ag$ModelBasedKappa, (observed - expected) / (1 - expected))
  expect_equal(out$metadata$calibration$PersonScoring, "EAP")
  fair <- out$report_tables[["Common-plan Fair Average"]]
  ci <- confint(wf$fit, "person")
  plan$theta <- ci[1L, 1L]
  expect_equal(fair$Lower[1L], mean(predict(wf$fit, plan, type = "response")))
  expect_true(all(fair$Lower <= fair$Estimate & fair$Estimate <= fair$Upper))

  md <- render_mfrm_report(wf, rater = "R03")
  expect_match(md, "Rater: R03", fixed = TRUE)
  expect_match(md, "Criterion halo screen", fixed = TRUE)
  expect_false(grepl("Rater: R01", md, fixed = TRUE))
  expect_error(render_mfrm_report(wf, rater = "Absent"), "known rater ID")
})

test_that("drift comparison requires a common reference and detects no change in identical data", {
  d <- load_mfrmr_data("example_core")
  facets <- c("Rater", "Criterion")
  group <- data.frame(Facet = "Rater", Level = unique(d$Rater), Group = "all", GroupValue = .75)
  fit <- fit_mfrm(d, "Person", facets, "Score", method = "JML", group_anchors = group,
                  spec = list(centering_scope = "remainder"))
  unlinked <- anchor_to_baseline(d, fit, "Person", facets, "Score", anchor_facets = "Rater")
  expect_true(all(is.na(unlinked$drift$Drift) & is.na(unlinked$drift$Flag)))
  expect_true(is.na(summary(unlinked)$n_flagged))
  reference <- data.frame(Facet = "Rater", Level = "R01")
  linked <- anchor_to_baseline(d, fit, "Person", facets, "Score", anchor_facets = "Rater",
                               drift_reference = reference)
  comparable <- linked$drift$Status == "conditional_on_baseline"
  expect_equal(sum(comparable), 3L)
  expect_equal(linked$drift$Drift[comparable], rep(0, 3), tolerance = 1e-6)
  expect_true(all(!linked$drift$Flag[comparable]))
  expect_error(anchor_to_baseline(subset(d, Rater != "R01"), fit, "Person", facets, "Score",
                                  anchor_facets = "Rater", anchor_policy = "error"), "absent")
})

test_that("assignment robustness measures preserve disconnected and parallel paths", {
  panels <- c(R1 = "A", R2 = "A", R3 = "B", R4 = "B")
  result <- lapply(0:2, function(n) {
    pair <- rbind(matrix(rep(c(1, 2), 3), ncol = 2, byrow = TRUE),
                  matrix(rep(c(3, 4), 3), ncol = 2, byrow = TRUE))
    if (n >= 1) pair[3, 2] <- 3
    if (n >= 2) pair[6, 2] <- 1
    d <- data.frame(Person = rep(1:6, each = 2), Rater = paste0("R", c(t(pair))), Score = rep(1:2, 6))
    check_mfrm_connectivity(d, "Person", "Rater", "Score",
                            rater_panels = panels)$rater_linkage$robustness
  })
  expect_equal(vapply(result, function(x) x$overview$FiedlerValue, 0), c(0, 0.1744269, 0.2907246),
               tolerance = 1e-6)
  expect_equal(vapply(result, function(x) x$panel_resistance$EffectiveResistance, 0), c(Inf, 29 / 12, 1.5))
})

test_that("multivariate report names the covariance and the reason for unavailable weights", {
  ratings <- load_mfrmr_data("example_core")
  d <- subset(ratings, Criterion %in% c("Content", "Organization"),
              select = c("Person", "Rater", "Criterion", "Score"))
  wide <- reshape(d, idvar = c("Person", "Rater"), timevar = "Criterion", direction = "wide")
  wf <- mfrm_workflow(ratings, "Person", c("Rater", "Criterion"), "Score") |>
    step_mfrm_gtheory(g_data = wide, g_facets = "Rater", g_weights = c(.6, .4),
                      g_scores = c("Score.Content", "Score.Organization"))
  out <- summary(wf)$report_tables
  eff <- out[["Effective weights at observed rating counts"]]
  expect_true(all(eff$Status[eff$Target == "Universe score"] == "non_psd_target_covariance"))
  expect_true(all(is.na(eff$EffectiveWeight[eff$Target == "Universe score"])))
  expect_equal(sum(eff$EffectiveWeight[eff$Target == "Absolute error"]), 1)
  expect_match(render_mfrm_report(wf, "html"), "supplied G-Study data", fixed = TRUE)
})

test_that("the sparse Fiedler solver accepts symmetric Matrix input", {
  skip_if_not_installed("RSpectra")
  d <- expand.grid(Person = 1:129, Rater = c("R1", "R2"))
  d$Score <- rep(1:2, length.out = nrow(d))
  conn <- check_mfrm_connectivity(d, "Person", "Rater", "Score", network = TRUE)
  expect_identical(conn$rater_linkage$robustness$overview$Status, "available")
  expect_equal(conn$rater_linkage$robustness$overview$FiedlerValue, 2, tolerance = 1e-8)
})

test_that("constant observed ratings do not create halo or central tendency screens", {
  d <- load_mfrmr_data("example_core")
  d$Score[d$Rater == "R04"] <- 4
  wf <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score") |>
    step_mfrm_measure(method = "JML_jeffreys") |> step_mfrm_diagnose()
  screen <- .mfrm_rater_card(wf, "R04", "Criterion", 3L)$tables
  expect_true(all(screen[["Central tendency screen"]]$Status == "no_observed_score_variation"))
  expect_true(all(screen[["Criterion halo screen"]]$Status == "no_observed_score_variation"))
  expect_true(all(is.na(screen[["Criterion halo screen"]]$ResidualCorrelation)))
})
