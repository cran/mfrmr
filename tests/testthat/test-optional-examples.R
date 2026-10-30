test_that("the posterior help example displays its saved results", {
  env <- new.env(parent = globalenv())
  shown <- capture.output(utils::example("read_mfrm_posterior", package = "mfrmr",
    local = env, echo = TRUE, ask = FALSE))
  expect_s3_class(env$posterior, "mfrm_posterior_summary")
  expect_match(paste(shown, collapse = "\n"),
    "Sampling checks do not establish model validity", fixed = TRUE)
})

test_that("supplementary help examples remain executable through example", {
  skip_on_cran()
  candidates <- normalizePath(c(file.path(test_path(), "..", ".."),
    file.path(test_path(), "..", "..", "00_pkg_src", "mfrmr")),
    winslash = "/", mustWork = FALSE)
  pkg_root <- candidates[file.exists(file.path(candidates, "inst", "examples", "execution-policy.csv"))][1L]
  skip_if(is.na(pkg_root), "source examples are not available")
  policy <- read.csv(file.path(pkg_root, "inst", "examples", "execution-policy.csv"))
  topics <- unique(policy$help_topic[as.logical(policy$replay_in_full_tests)])
  pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf)
  withr::defer({grDevices::dev.off(); unlink(pdf)})
  for (topic in topics) {
    env <- new.env(parent = globalenv())
    run_example <- function() utils::example(topic, package = "mfrmr",
      character.only = TRUE, local = env, echo = FALSE,
      ask = FALSE, run.dontrun = TRUE)
    set.seed(1)
    if (identical(topic, "estimate_bias")) {
      expect_warning(run_example(),
        "^Category support is retained but requires review:")
    } else expect_no_error(run_example())
  }
})
