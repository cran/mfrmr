test_that("normal rater intervals require selection and reuse saved estimates", {
  fit <- structure(list(raters = data.frame(Rater = c("B", "A"), Estimate = c(.4, -.2),
    PredictionSE = c(.3, .1), Lower = c(-99, -99), Upper = c(99, 99)),
    settings = list(fixed_rater_sd = .7),
    checks = list(NumericalReady = TRUE, InformationPositive = TRUE)), class = "mfrm_random_rater")
  before <- fit
  local_mocked_bindings(mfrm_random_rater_fit = function(...) stop("Unexpected refit"),
    mfrm_random_rater_objective = function(...) stop("Unexpected objective"))
  ci <- confint(fit, parm = "raters", level = .8)
  expect_equal(unname(ci[, "Lower"]), c(.4, -.2) - qnorm(.9) * c(.3, .1))
  expect_equal(unname(ci[, "Upper"]), c(.4, -.2) + qnorm(.9) * c(.3, .1))
  expect_identical(rownames(ci), c("B", "A"))
  expect_match(attr(ci, "target"), "uncentered")
  expect_match(attr(ci, "note"), "not established")
  level <- 1 - .Machine$double.eps / 2
  high <- confint(fit, parm = "raters", level = level)
  z <- sqrt(qchisq(level, 1))
  expect_equal(unname(high[, "Lower"]), fit$raters$Estimate-z*fit$raters$PredictionSE)
  expect_equal(unname(high[, "Upper"]), fit$raters$Estimate+z*fit$raters$PredictionSE)
  expect_error(confint(fit, parm = "raters", level = Inf), "level")
  default <- plot(fit, draw = FALSE)$data
  expect_true(all(is.na(default$table$Lower) & is.na(default$table$Upper)))
  expect_equal(default$table$Estimate, fit$raters$Estimate)
  expect_match(default$alt_text, "Point estimates only")
  expect_error(plot(fit, draw = FALSE, style = "precision"), "explicit")
  expect_error(plot(fit, draw = FALSE, sort = "uncertainty"), "explicit")
  chosen <- plot(fit, draw = FALSE, intervals = "normal", show_notes = FALSE)$data
  expect_equal(chosen$table$Lower, unname(confint(fit, parm = "raters")[, "Lower"]))
  expect_match(chosen$notes$Text, "not established")
  expect_identical(chosen$settings$rater_intervals, "normal")
  expect_match(plot(fit, draw = FALSE, style = "distribution",
    intervals = "normal")$data$caption, "no intervals displayed")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- as_ggplot(plot(fit, draw = FALSE))
    expect_equal(nrow(p$layers[[2]]$data), 0L)
    expect_silent(ggplot2::ggplotGrob(p))
    expect_silent(ggplot2::ggplotGrob(as_ggplot(plot(fit, draw = FALSE,
      intervals = "normal", style = "precision"))))
  }
  for (key in c("NumericalReady", "InformationPositive", "EstimatedVarianceBoundary", "EstimatedPersonVarianceBoundary")) {
    changed <- fit; changed$checks[[key]] <- grepl("Boundary", key)
    expect_true(all(is.na(confint(changed, parm = "raters"))))
  }
  changed <- fit; changed$raters$PredictionSE <- c(NA_real_, .1)
  expect_true(all(is.na(confint(changed, parm = "raters")[1, ])))
  expect_true(all(is.finite(confint(changed, parm = "raters")[2, ])))
  path <- tempfile(); saveRDS(fit, path)
  expect_identical(confint(readRDS(path), parm = "raters"), confint(fit, parm = "raters"))
  expect_identical(fit, before)
})
