test_that("JML rejects an explicit unsupported internal category contrast", {
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

  condition <- tryCatch(
    suppressWarnings(fit_mfrm(
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
    )),
    error = identity
  )

  expect_s3_class(condition, "mfrmr_category_readiness_error")
  expect_identical(condition$readiness$CategoryState,
                   "unsupported_coordinate")
  expect_identical(
    condition$category_support$support_table$UnsupportedCategory,
    "3"
  )
  expect_identical(
    condition$category_support$unsupported_contrasts$Direction,
    "Step3 - Step4"
  )
  expect_match(condition$readiness$ReasonCodes,
               "step_scope_category_unobserved", fixed = TRUE)
})
