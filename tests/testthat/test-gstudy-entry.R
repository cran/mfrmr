main_effects_entry_fixture <- function() {
  structure(list(variance_components = data.frame(
    Source = c("Person", "Rater", "Task", "Residual"), Variance = c(1, .3, .2, .8)),
    design = list(object_facet = "Person", random_facets = c("Rater", "Task"),
      calculation_version = 2L, identification_status = "identified")),
    class = c("mfrm_generalizability", "list"))
}

test_that("one and several scores use the same G-to-D entry and count aliases", {
  for (scores in list("Content", c("Content", "Organization"))) {
    g <- mfrm_g_study(mfrmr_example_gtheory_crossed, scores)
    before <- serialize(g, NULL)
    expected <- mfrm_d_study(g, data.frame(Raters = c(2, 4), Tasks = c(3, 6)))
    expect_identical(mfrm_d_study(g, list(n_Task = c(3, 6), n_Rater = c(2, 4))), expected)
    expect_identical(mfrm_d_study(g, data.frame(Task = c(3, 6), Rater = c(2, 4))), expected)
    expect_identical(summary(expected), expected$coefficients)
    expect_identical(serialize(g, NULL), before)
    expect_error(mfrm_d_study(g, residual_scaling = "none"), "must be empty")
    expect_error(mfrm_d_study(g, typo = TRUE), "must be empty")
  }
  exports <- getNamespaceExports("mfrmr")
  expect_true(all(c("mfrm_g_study", "mfrm_d_study") %in% exports))
  expect_false(any(c("mfrm_multivariate_gstudy", "mfrm_multivariate_d_study") %in% exports))
})

test_that("nested counts and generic labels retain their source meanings", {
  nested <- mfrm_g_study(mfrmr_example_gtheory_nested, "Content", nesting = c(Rater = "Task"))
  d <- mfrm_d_study(nested, list(n_Rater = c(1, 2), n_Task = 4))
  expect_identical(d, mfrm_d_study(nested, data.frame(Raters = c(1, 2), Tasks = 4)))
  expect_identical(d$gstudy$design$nesting, c(Rater = "Task"))
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed, "Content",
    facets = c(Assessor = "Rater", Occasion = "Task"))
  expect_identical(mfrm_d_study(g, list(n_Occasion = 3, n_Assessor = 2)),
    mfrm_d_study(g, data.frame(Assessor = 2, Occasion = 3)))
  collision <- mfrm_g_study(mfrmr_example_gtheory_crossed, "Content",
    facets = c(Rater = "Rater", n_Rater = "Task"))
  expect_identical(mfrm_d_study(collision, list(n_Rater = 2, n_n_Rater = 3)),
    mfrm_d_study(collision, data.frame(Rater = 2, n_Rater = 3)))
  for (bad in list(data.frame(Rater = 2, n_Rater = 2, Task = 3),
      data.frame(n_Rater = 2), data.frame(n_Rater = 2, n_Task = 3, Person = 50)))
    expect_error(mfrm_d_study(nested, bad), "exactly one count|missing count")
})

test_that("main-effects dispatch preserves its formula and refuses unused requests", {
  g <- main_effects_entry_fixture()
  d <- mfrm_d_study(g, list(n_Rater = 2, n_Task = 4))
  expect_equal(d$G, 1 / (1 + .8 / 8))
  expect_equal(d$Phi, 1 / (1 + .8 / 8 + .3 / 2 + .2 / 4))
  expect_identical(summary(d), as.data.frame(d))
  for (args in list(list(weights = 1), list(typo = TRUE),
      list(object_facet = "Person"), list(random_facets = "Rater")))
    expect_error(do.call(mfrm_d_study, c(list(g, list(Rater = 2, Task = 4)), args)), "must be empty")
  expect_error(summary(d, intervals = TRUE), "must be empty")
  calls <- list()
  local_mocked_bindings(mfrm_generalizability = function(...) {
    calls[[length(calls) + 1L]] <<- list(...); g
  })
  fit <- structure(list(), class = "mfrm_fit")
  actual <- mfrm_d_study(fit, list(Rater = 2, Task = 4),
    object_facet = "Person", random_facets = c("Rater", "Task"),
    residual_scaling = "none", missing = "omit")
  expect_identical(actual, mfrm_d_study(g, list(Rater = 2, Task = 4), residual_scaling = "none"))
  expect_identical(calls[[1]]$missing, "omit")
  expect_identical(calls[[1]]$random_facets, c("Rater", "Task"))
})

test_that("saved G and D objects keep their classes and can be used without refitting", {
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed, c("Content", "Organization"))
  d <- mfrm_d_study(g, list(n_Rater = c(2, 4), n_Task = 3), weights = c(Content = 1, Organization = -1))
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(list(g = g, d = d), path)
  saved <- readRDS(path)
  local_mocked_bindings(mfrm_g_study = function(...) stop("Unexpected G-study"),
    mfrm_generalizability = function(...) stop("Unexpected main-effects fit"))
  expect_identical(mfrm_d_study(saved$g, saved$d$design_grid, weights = saved$d$weights), saved$d)
  expect_identical(summary(saved$d), saved$d$coefficients)
  expect_s3_class(saved$g, "mfrm_multivariate_gstudy")
  expect_s3_class(saved$d, "mfrm_multivariate_d_study")
})
