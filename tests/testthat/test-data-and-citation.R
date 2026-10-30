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

test_that("all dataset keys preserve data() contents and copy-on-modify", {
  keys <- list_mfrmr_data()
  objects <- c("mfrmr_example_core", "mfrmr_example_bias",
    "mfrmr_example_operational", "ej2021_study1", "ej2021_study2",
    "ej2021_combined", "ej2021_study1_itercal", "ej2021_study2_itercal",
    "ej2021_combined_itercal")
  expect_length(keys, length(objects))
  for (i in seq_along(keys)) {
    reference <- new.env(parent = emptyenv())
    utils::data(list = objects[i], package = "mfrmr", envir = reference)
    expected <- reference[[objects[i]]]
    ratings <- load_mfrmr_data(keys[i])
    expect_identical(ratings, expected, info = keys[i])
    expect_identical(load_mfrmr_data(toupper(keys[i])), expected)
    ratings$Score[1] <- -99
    names(ratings)[1] <- "changed"
    attr(ratings, "user_note") <- "local edit"
    expect_identical(load_mfrmr_data(keys[i]), expected, info = keys[i])
  }
})
