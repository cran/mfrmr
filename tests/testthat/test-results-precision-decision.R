test_that("result and report decisions retain saved precision without overriding fit readiness", {
  fit <- make_toy_fit(method = "MML", maxit = 200)
  dx <- make_toy_diagnostics(fit)
  expect_true(isTRUE(dx$precision_profile$SupportsFormalInference[1L]))
  res <- mfrm_results(fit, diagnostics = dx, include = "fit", compute = "never")
  # Reading and reporting saved decisions must not run another diagnostic pass.
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected diagnostics"),
                        .package = "mfrmr")
  decision <- summary(res)$decision
  expect_identical(decision$FormalInference, summary(dx)$decision$FormalInference)
  expect_identical(decision$FormalInference, "Yes")
  expect_false(grepl("not been evaluated", decision$Why, fixed = TRUE))
  expect_identical(mfrm_report(res)$decision, decision)
  expect_true(mfrm_report(res)$precision_evidence_summary$SupportsFormalInference)

  missing <- res
  missing$diagnostics <- NULL
  expect_identical(summary(missing)$decision$FormalInference, "No")
  expect_match(summary(missing)$decision$Why, "not been evaluated", fixed = TRUE)
  missing$diagnostics <- dx
  missing$diagnostics$precision_profile <- NULL
  expect_match(summary(missing)$decision$Why, "not been evaluated", fixed = TRUE)

  unsupported <- res
  unsupported$diagnostics$precision_profile$SupportsFormalInference <- FALSE
  unsupported$diagnostics$precision_profile$PrecisionTier <- "exploratory"
  expect_identical(summary(unsupported)$decision$FormalInference, "No")
  expect_match(summary(unsupported)$decision$Why, "tier: exploratory", fixed = TRUE)
  expect_identical(mfrm_report(unsupported)$decision, summary(unsupported)$decision)
  expect_false(mfrm_report(unsupported)$precision_evidence_summary$SupportsFormalInference)

  restricted <- res
  restricted$fit_readiness$InferenceReady <- FALSE
  restricted$fit_readiness$FitReadiness <- "review"
  restricted$fit_readiness$NumericalState <- "review"
  expect_identical(summary(restricted)$decision$FormalInference, "No")
  expect_match(summary(restricted)$decision$Why, "Numerical convergence requires review",
               fixed = TRUE)
  expect_identical(mfrm_report(restricted)$decision, summary(restricted)$decision)
})
