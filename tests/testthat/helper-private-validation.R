# These checks need local evidence; missing records are not numerical successes.
skip_if_private_validation_missing <- function(files) {
  paths <- testthat::test_path("..", "..", "inst", "validation", files)
  testthat::skip_if_not(
    all(file.exists(paths)),
    "Local-only maintainer evidence is unavailable in this checkout."
  )
}
