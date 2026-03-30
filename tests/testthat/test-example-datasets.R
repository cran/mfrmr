test_that("packaged example datasets are available and well-formed", {
  keys <- list_mfrmr_data()
  expect_true(all(c("example_core", "example_bias") %in% keys))

  core <- load_mfrmr_data("example_core")
  bias <- load_mfrmr_data("example_bias")

  for (dat in list(core, bias)) {
    expect_true(is.data.frame(dat))
    expect_true(all(c("Study", "Person", "Rater", "Criterion", "Score", "Group") %in%
                      names(dat)))
    expect_gte(length(unique(dat$Person)), 30)
    expect_true(all(table(dat$Score) >= 10))
  }
})

test_that("example datasets are available through utils::data", {
  env <- new.env(parent = emptyenv())
  utils::data(list = "mfrmr_example_core", package = "mfrmr", envir = env)
  utils::data(list = "mfrmr_example_bias", package = "mfrmr", envir = env)

  expect_true(exists("mfrmr_example_core", envir = env, inherits = FALSE))
  expect_true(exists("mfrmr_example_bias", envir = env, inherits = FALSE))
  expect_s3_class(env$mfrmr_example_core, "data.frame")
  expect_s3_class(env$mfrmr_example_bias, "data.frame")
})

test_that("example_bias supports non-null DIF and bias help examples", {
  dat <- load_mfrmr_data("example_bias")
  fit <- fit_mfrm(dat, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 50)
  diag <- diagnose_mfrm(fit, residual_pca = "none")

  dif <- analyze_dif(fit, diag, facet = "Criterion",
                     group = "Group", data = dat, method = "residual")
  expect_true(max(abs(dif$dif_table$Contrast), na.rm = TRUE) > 0.3)

  bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion",
                        max_iter = 2)
  expect_true(max(abs(suppressWarnings(as.numeric(bias$table$`Bias Size`))), na.rm = TRUE) > 0.5)
  expect_true(all(as.character(bias$table$InferenceTier) == "screening"))
})
