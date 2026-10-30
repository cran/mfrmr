# MML EM checkpoint schema, identity, and exact-resume boundary tests.

checkpoint_fit <- function(data, checkpoint, maxit = 3L,
                           engine = "em", quad_points = 5L) {
  suppressWarnings(
    fit_mfrm(
      data, "Person", c("Rater", "Criterion"), "Score",
      method = "MML", quad_points = quad_points, maxit = maxit,
      mml_engine = engine,
      checkpoint = list(file = checkpoint, every_iter = 1L)
    )
  )
}

test_that("MML EM checkpoint writes a file when supplied", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)
  fit <- checkpoint_fit(toy, ckpt, maxit = 5L)
  expect_s3_class(fit, "mfrm_fit")
  expect_true(file.exists(ckpt))
  saved <- readRDS(ckpt)
  expect_identical(as.character(saved$.mfrm_checkpoint_kind), "mml_em")
  expect_identical(saved$schema_id, "mfrmr_mml_em_checkpoint")
  expect_identical(saved$schema_version, 2L)
  expect_identical(saved$identity$engine_stage, "pure_em")
  expect_identical(saved$identity$quadrature_identity$rule,
                   "gauss_hermite_standard_normal_recurrence_v2")
  legacy_identity <- saved$identity
  legacy_identity$quadrature_identity$rule <- "gauss_hermite_standard_normal_golub_welsch_v1"
  expect_error(mfrmr:::mfrm_validate_checkpoint_payload(saved, legacy_identity),
               "checkpoint identity does not match", fixed = TRUE)
  expect_true(nzchar(saved$identity$data_objective_fingerprint))
  expect_true(is.numeric(saved$par))
  expect_identical(
    names(saved$par), saved$identity$parameter_layout$parameter_names
  )
  expect_true(is.integer(saved$next_iter))
  expect_true(is.integer(saved$last_completed_iter))
  expect_true(saved$completed)
})

test_that("MML EM resumes from an existing checkpoint", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)
  # First short run -- write checkpoint.
  fit_a <- checkpoint_fit(toy, ckpt, maxit = 3L)
  expect_true(file.exists(ckpt))
  saved <- readRDS(ckpt)
  initial_next <- saved$next_iter

  # Second run -- should consume the checkpoint and continue.
  fit_b <- checkpoint_fit(toy, ckpt, maxit = 6L)
  expect_s3_class(fit_b, "mfrm_fit")
  saved2 <- readRDS(ckpt)
  expect_gte(saved2$next_iter, initial_next)
})

test_that("checkpoint rejects a changed objective with the same parameter count", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)
  checkpoint_fit(toy, ckpt, maxit = 3L)

  changed <- toy
  row <- which(changed$Score > min(changed$Score))[1]
  changed$Score[row] <- changed$Score[row] - 1
  expect_error(
    checkpoint_fit(changed, ckpt, maxit = 6L),
    "checkpoint identity does not match",
    fixed = TRUE
  )

  expect_error(
    checkpoint_fit(toy, ckpt, maxit = 6L, quad_points = 7L),
    "checkpoint identity does not match",
    fixed = TRUE
  )
})

test_that("hybrid checkpoint writes and reuses its completed warm start", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)

  first <- checkpoint_fit(toy, ckpt, maxit = 5L, engine = "hybrid")
  expect_s3_class(first, "mfrm_fit")
  expect_true(file.exists(ckpt))
  saved <- readRDS(ckpt)
  expect_identical(saved$identity$engine_stage, "hybrid_em_warm_start")
  expect_true(saved$completed)

  messages <- capture.output(
    second <- checkpoint_fit(toy, ckpt, maxit = 5L, engine = "hybrid"),
    type = "message"
  )
  expect_s3_class(second, "mfrm_fit")
  expect_true(any(grepl(
    "Loaded completed hybrid EM warm-start checkpoint", messages,
    fixed = TRUE
  )))
})

test_that("completed pure-EM checkpoint cannot re-enter the same maxit", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)
  checkpoint_fit(toy, ckpt, maxit = 3L)
  before <- readRDS(ckpt)

  expect_error(
    checkpoint_fit(toy, ckpt, maxit = 3L),
    "already reached the requested `maxit`",
    fixed = TRUE
  )
  after <- readRDS(ckpt)
  expect_identical(after$last_completed_iter, before$last_completed_iter)
  expect_identical(after$ll_trace, before$ll_trace)
})

test_that("checkpoint controls and payloads fail closed", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)

  for (invalid in list(NA_integer_, Inf, 0, 1.5)) {
    expect_error(
      suppressWarnings(fit_mfrm(
        toy, "Person", c("Rater", "Criterion"), "Score",
        method = "MML", quad_points = 5, maxit = 3,
        mml_engine = "em",
        checkpoint = list(file = ckpt, every_iter = invalid)
      )),
      "finite positive integer",
      fixed = TRUE
    )
  }

  saveRDS(list(.mfrm_checkpoint_kind = "mml_em"), ckpt)
  expect_error(
    checkpoint_fit(toy, ckpt, maxit = 6L),
    "unsupported or legacy MML EM schema",
    fixed = TRUE
  )
})
