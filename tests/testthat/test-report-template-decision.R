test_that("GPCM report introduction uses the assembled report decision", {
  skip_if_not_installed("knitr")
  template <- readLines(system.file("templates", "gpcm-report.Rmd", package = "mfrmr"))
  introduction <- grep("^\\*\\*(Current analysis|Formal-inference status):", template, value = TRUE)
  expect_length(introduction, 2L)
  env <- new.env(parent = baseenv())
  env$fit_summary <- list(decision = data.frame(
    Interpretation = "Fit-only decision must not appear", FormalInference = "Yes"))
  env$report <- list(decision = data.frame(
    Interpretation = "Collected diagnostics require review", FormalInference = "No"))
  text <- knitr::knit(text = introduction, envir = env, quiet = TRUE)
  expect_match(text, "Collected diagnostics require review", fixed = TRUE)
  expect_match(text, "**Formal-inference status:** No", fixed = TRUE)
  expect_false(grepl("Fit-only decision|Yes", text))
})
