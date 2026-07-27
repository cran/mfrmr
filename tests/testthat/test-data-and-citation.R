test_that("packaged data aliases and loaders are available", {
  aliases <- mfrmr::list_mfrmr_data()
  expect_true(is.character(aliases))
  expect_true(all(c("study1", "study2", "combined") %in% aliases))

  d <- mfrmr::load_mfrmr_data("study1")
  expect_s3_class(d, "data.frame")
  expect_true(all(c("Person", "Rater", "Criterion", "Score") %in% names(d)))
  expect_gt(nrow(d), 0)
})

test_that("citation metadata is available", {
  # `test_local()` loads the source tree without installing it, so base R may
  # warn that no installed package was found even though pkgload exposes the
  # package metadata and inst/CITATION correctly. Installed-package checks do
  # not take this branch.
  cit <- suppressWarnings(utils::citation("mfrmr"))
  expect_true(length(cit) >= 1)

  cit_file <- system.file("CITATION", package = "mfrmr")
  expect_true(nzchar(cit_file))
  expect_true(file.exists(cit_file))
})
