test_that("JML preserves explicit zero-count score support when requested", {
  set.seed(20260515)
  dat <- expand.grid(
    Person = paste0("P", sprintf("%02d", seq_len(10))),
    Rater = c("R1", "R2"),
    Criterion = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  dat$Score <- sample(
    c(0L, 1L, 2L, 4L),
    size = nrow(dat),
    replace = TRUE,
    prob = c(0.15, 0.30, 0.35, 0.20)
  )
  dat$Score[seq_len(4)] <- c(0L, 1L, 2L, 4L)

  expect_false(3L %in% dat$Score)

  fit <- suppressWarnings(
    fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      rating_min = 0,
      rating_max = 4,
      keep_original = TRUE,
      method = "JML",
      model = "RSM",
      maxit = 30,
      reltol = 1e-4
    )
  )

  expect_s3_class(fit, "mfrm_fit")
  expect_equal(as.integer(fit$summary$Categories[1]), 5L)
  expect_true(isTRUE(fit$summary$Converged[1]))
  expect_equal(fit$prep$score_map$OriginalScore, 0:4)
  expect_equal(fit$prep$score_map$InternalScore, 0:4)
  expect_true(all(is.finite(fit$steps$Estimate)))

  diagnostics <- suppressWarnings(
    diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy")
  )
  expect_s3_class(diagnostics, "mfrm_diagnostics")
  expect_equal(nrow(diagnostics$obs), nrow(dat))

  dq <- data_quality_report(
    fit,
    data = dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score"
  )
  expect_s3_class(dq, "mfrm_data_quality")
  expect_true("score_support_review" %in% names(dq))

  score_3 <- dq$score_support_review[
    dq$score_support_review$Score == 3L,
    ,
    drop = FALSE
  ]
  expect_equal(nrow(score_3), 1L)
  expect_true(isTRUE(score_3$ZeroCount[1]))
  expect_true(isTRUE(score_3$WeaklyIdentified[1]))
  expect_true(any(dq$quality_flags$Area == "Score support"))
})
