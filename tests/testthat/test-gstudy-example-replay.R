test_that("the saved incomplete G-study example retains its data and targets", {
  path <- system.file("examples", "gstudy-sampling.rds", package = "mfrmr", mustWork = TRUE)
  saved <- readRDS(path)
  observed <- mfrmr_example_gtheory_crossed[-c(1, 8, 23), ]
  expected <- mfrm_g_study(observed, c("Content", "Organization"),
    method = "minque0")
  # BLAS/platform rounding may differ from the saved example's environment.
  # Match the independent sampling-covariance tests' numerical tolerance.
  expect_equal(saved$study, expected, tolerance = 1e-9)
  expect_false(saved$study$design$complete)
  expect_identical(attr(saved$variance, "score"), "Content")
  expect_identical(attr(saved$covariance, "score"), c("Content", "Organization"))
  for (result in saved[c("variance", "covariance")]) {
    expect_identical(attr(result, "assumption"), "normal")
    expect_identical(rownames(result), names(expected$components))
    expect_equal(attr(result, "summary")$SamplingVariance, unname(diag(result)))
    expect_true(all(is.finite(result)))
  }
  # Reading the saved matrices and the regeneration recipe must not recompute P.
  local_mocked_bindings(.mfrm_mvgt_moment_polynomial = function(...) stop("Unexpected moment calculation"))
  recipe <- new.env(parent = baseenv())
  sys.source(system.file("examples", "gstudy-sampling.R", package = "mfrmr",
    mustWork = TRUE), envir = recipe)
  expect_true(is.function(recipe$gstudy_sampling_example))
  expect_identical(readRDS(path), saved)
})
