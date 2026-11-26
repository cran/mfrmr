test_that("shared and fixed steps survive inference, serialization and prediction", {
  ratings <- subset(load_mfrmr_data("example_core"), Criterion %in% c("Accuracy", "Content"))
  declaration <- data.frame(Level = c("Accuracy", "Content"), From = 1, To = 2, Group = "first")
  fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score", model = "PCM",
    step_facet = "Criterion", step_constraints = declaration)
  shared <- fit$steps[fit$steps$StepIndex == 1L, ]
  expect_equal(shared$Estimate[1], shared$Estimate[2], tolerance = 1e-10)
  expect_equal(shared$ModelSE[1], shared$ModelSE[2], tolerance = 1e-10)
  expect_true(all(is.finite(shared$ModelSE)))
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(predict(restored, type = "probs"), predict(fit, type = "probs"))
  expect_equal(confint(restored, component = "steps"), confint(fit, component = "steps"))
  declaration$Value <- c(-.7, NA)
  fixed <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score", model = "PCM",
    step_facet = "Criterion", step_constraints = declaration)
  shared <- fixed$steps[fixed$steps$StepIndex == 1L, ]
  expect_equal(shared$Estimate, rep(-.7, 2), tolerance = 1e-10)
  ci <- confint(fixed, component = "steps")
  expect_true(all(is.na(ci[fixed$steps$StepIndex == 1L, ])))
  free <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score", model = "PCM", step_facet = "Criterion")
  comparison <- anova(fit, free, nested = TRUE)
  expect_true(nrow(comparison$lrt) == 1L)
})

test_that("contradictory step declarations fail before fitting", {
  ratings <- load_mfrmr_data("example_core")
  fixed <- data.frame(Level = "shared", From = 1:3, To = 2:4, Value = c(-1, 0, 2))
  for (method in c("MML", "JML")) expect_error(fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
    method = method, step_constraints = fixed), "sum.to.zero|inconsistent|incompatible")
  expect_error(fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
    step_constraints = fixed[c(1, 1), ]), "once|duplicate")
})
