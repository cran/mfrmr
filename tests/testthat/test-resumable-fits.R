# Tests for the MML EM checkpoint / resume scaffolding added in 0.1.6.
# We run a short EM fit with a checkpoint, force the file to be
# re-read by spawning a second `fit_mfrm()` call, and verify the
# resumed fit converges from the saved iteration rather than from
# scratch.

test_that("MML EM checkpoint writes a file when supplied", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", quad_points = 5, maxit = 5,
             mml_engine = "em",
             checkpoint = list(file = ckpt, every_iter = 1L))
  ))
  expect_s3_class(fit, "mfrm_fit")
  expect_true(file.exists(ckpt))
  saved <- readRDS(ckpt)
  expect_identical(as.character(saved$.mfrm_checkpoint_kind), "mml_em")
  expect_true(is.numeric(saved$par))
  expect_true(is.integer(saved$next_iter))
})

test_that("MML EM resumes from an existing checkpoint", {
  toy <- load_mfrmr_data("example_core")
  ckpt <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpt), add = TRUE)
  # First short run -- write checkpoint.
  fit_a <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", quad_points = 5, maxit = 3,
             mml_engine = "em",
             checkpoint = list(file = ckpt, every_iter = 1L))
  ))
  expect_true(file.exists(ckpt))
  saved <- readRDS(ckpt)
  initial_next <- saved$next_iter

  # Second run -- should consume the checkpoint and continue.
  fit_b <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", quad_points = 5, maxit = 6,
             mml_engine = "em",
             checkpoint = list(file = ckpt, every_iter = 1L))
  ))
  expect_s3_class(fit_b, "mfrm_fit")
  saved2 <- readRDS(ckpt)
  expect_gte(saved2$next_iter, initial_next)
})
