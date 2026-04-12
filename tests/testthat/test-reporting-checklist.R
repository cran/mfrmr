test_that("reporting_checklist returns a bundle with checklist coverage tables", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))

  chk <- reporting_checklist(fit, diagnostics = diag, bias_results = list(rater_criterion = bias))

  expect_s3_class(chk, "mfrm_reporting_checklist")
  expect_true(is.data.frame(chk$checklist))
  expect_true(is.data.frame(chk$summary))
  expect_true(is.data.frame(chk$section_summary))
  expect_true(is.data.frame(chk$software_scope))
  expect_true(is.data.frame(chk$visual_scope))
  expect_true("InterpretationCheck" %in% names(chk$visual_scope))
  expect_true(all(c("mfrmr native", "FACETS", "ConQuest", "SPSS") %in% chk$software_scope$Software))
  expect_true("Category probability surface" %in% chk$visual_scope$Visualization)
  expect_match(
    chk$visual_scope$ThreeDStatus[chk$visual_scope$Visualization == "Category probability surface"][1],
    "advanced surface data only",
    fixed = TRUE
  )
  expect_match(
    chk$visual_scope$InterpretationCheck[
      chk$visual_scope$Visualization == "Category probability surface"
    ][1],
    "category_support",
    fixed = TRUE
  )
  expect_match(
    chk$software_scope$Boundary[chk$software_scope$Software == "SPSS"][1],
    "native SPSS integration is not implemented",
    fixed = TRUE
  )
  expect_match(
    chk$software_scope$Boundary[chk$software_scope$Software == "FACETS"][1],
    "Results remain mfrmr estimates",
    fixed = TRUE
  )
  expect_true(all(c(
    "Section", "Item", "Available", "DraftReady", "ReadyForAPA", "Severity",
    "Priority", "SourceComponent", "Detail", "PlotHelper",
    "DrawFreeRoute", "PlotReturnClass", "NextAction"
  ) %in% names(chk$checklist)))
  expect_identical(chk$checklist$DraftReady, chk$checklist$ReadyForAPA)
  expect_true(any(chk$checklist$Item == "PCA of residuals"))
  expect_true(any(chk$checklist$Item == "Facet pairs tested"))
  expect_true(any(chk$checklist$Item == "QC / facet dashboard"))
  expect_true(any(chk$checklist$Item == "Residual PCA visuals"))
  expect_true(any(chk$checklist$Item == "Connectivity / design-matrix visual"))
  expect_true(any(chk$checklist$Item == "Inter-rater / displacement visuals"))
  expect_true(any(chk$checklist$Item == "Strict marginal visuals"))
  expect_true(any(chk$checklist$Item == "Bias / DIF visuals"))
  expect_true(any(chk$checklist$Item == "Precision / information curves"))
  expect_true(chk$checklist$Available[chk$checklist$Item == "PCA of residuals"][1])
  expect_true(chk$checklist$Available[chk$checklist$Item == "QC / facet dashboard"][1])
  expect_true(chk$checklist$Available[chk$checklist$Item == "Residual PCA visuals"][1])
  expect_false(chk$checklist$Available[chk$checklist$Item == "Strict marginal visuals"][1])
  expect_false(chk$checklist$ReadyForAPA[chk$checklist$Item == "95% confidence intervals"][1])
  expect_false(chk$checklist$ReadyForAPA[chk$checklist$Item == "Separation / strata / reliability"][1])
  expect_false(chk$checklist$ReadyForAPA[chk$checklist$Item == "Strict marginal visuals"][1])
  expect_true(any(nzchar(chk$checklist$NextAction)))
  visual_rows <- chk$checklist[chk$checklist$Section == "Visual Displays", , drop = FALSE]
  expect_true(all(nzchar(visual_rows$PlotHelper)))
  expect_true(all(nzchar(visual_rows$DrawFreeRoute)))
  expect_true(all(visual_rows$PlotReturnClass == "mfrm_plot_data"))
  expect_true(all(c("DraftReady", "ReadyForAPA", "NeedsDraftWork", "NeedsAction") %in% names(chk$summary)))

  s_chk <- summary(chk)
  expect_s3_class(s_chk, "summary.mfrm_reporting_checklist")
  expect_true(is.data.frame(s_chk$overview))
  expect_true(is.data.frame(s_chk$section_summary))
  expect_true(is.data.frame(s_chk$software_scope))
  expect_true(is.data.frame(s_chk$visual_scope))
  expect_true(is.data.frame(s_chk$priority_summary))
  expect_true(is.data.frame(s_chk$action_items))
  expect_true(is.data.frame(s_chk$settings))
  expect_gt(nrow(s_chk$section_summary), 0)
  printed <- capture.output(print(s_chk))
  expect_true(any(grepl("mfrmr Reporting Checklist Summary", printed, fixed = TRUE)))
  expect_false(any(grepl("External software scope", printed, fixed = TRUE)))
  expect_false(any(grepl("Visual scope", printed, fixed = TRUE)))
  expect_true(any(grepl("Detailed software and visual scope tables", printed, fixed = TRUE)))
})

test_that("reporting_checklist surfaces latent-regression reporting tasks", {
  set.seed(141)
  persons <- paste0("P", sprintf("%02d", seq_len(24)))
  items <- paste0("I", seq_len(4))
  group <- rep(c("high", "low"), length.out = length(persons))
  theta <- ifelse(group == "high", 0.7, -0.4) + stats::rnorm(length(persons), sd = 0.5)
  item_beta <- seq(-0.6, 0.6, length.out = length(items))
  dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
  eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
  dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))
  person_data <- data.frame(
    Person = persons,
    Group = group,
    stringsAsFactors = FALSE
  )

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    "Item",
    "Score",
    method = "MML",
    model = "RSM",
    population_formula = ~ Group,
    person_data = person_data,
    quad_points = 5,
    maxit = 35
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))

  chk <- reporting_checklist(fit, diagnostics = diag)
  pop_rows <- chk$checklist[chk$checklist$Section == "Population Model", , drop = FALSE]
  conquest_scope <- chk$software_scope[chk$software_scope$Software == "ConQuest", , drop = FALSE]

  expect_equal(nrow(pop_rows), 6L)
  expect_match(conquest_scope$CurrentSupport[1], "candidate fit", fixed = TRUE)
  expect_true(all(c(
    "Latent-regression basis",
    "Population coefficients and residual variance",
    "Model-matrix covariate coding",
    "Complete-case omission audit",
    "Population-model posterior scoring wording",
    "ConQuest overlap wording"
  ) %in% pop_rows$Item))
  expect_true(pop_rows$Available[pop_rows$Item == "Model-matrix covariate coding"][1])
  expect_match(pop_rows$Detail[pop_rows$Item == "Model-matrix covariate coding"][1], "Group", fixed = TRUE)
  expect_false(pop_rows$DraftReady[pop_rows$Item == "ConQuest overlap wording"][1])
  expect_match(pop_rows$NextAction[pop_rows$Item == "ConQuest overlap wording"][1], "documented latent-regression MML comparison scope", fixed = TRUE)

  s_chk <- summary(chk)
  expect_true("Population Model" %in% s_chk$section_summary$Section)
})

test_that("reporting_checklist surfaces non-numeric bias screening statistics", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))
  bias$table$t <- "not-a-number"

  chk <- reporting_checklist(fit, diagnostics = diag, bias_results = list(rater_criterion = bias))
  row <- chk$checklist[chk$checklist$Item == "Screen-positive interactions", , drop = FALSE]

  expect_match(row$Detail[1], "non-numeric screening statistics", fixed = TRUE)
  expect_false(row$ReadyForAPA[1])
})

test_that("reporting_checklist surfaces failed bias-collection pairs", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))
  bias_collection <- structure(
    list(
      by_pair = list(rater_criterion = bias),
      errors = data.frame(
        Interaction = "Task x Criterion",
        Facets = "Task x Criterion",
        Error = "forced pair failure",
        stringsAsFactors = FALSE
      )
    ),
    class = c("mfrm_bias_collection", "mfrm_bundle", "list")
  )

  chk <- reporting_checklist(fit, diagnostics = diag, bias_results = bias_collection)
  row <- chk$checklist[chk$checklist$Item == "Screen-positive interactions", , drop = FALSE]

  expect_match(row$Detail[1], "failed", fixed = TRUE)
  expect_false(row$ReadyForAPA[1])
  expect_identical(chk$settings$bias_error_count, 1L)
})

test_that("reporting_checklist rejects malformed bias_results inputs early", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))

  expect_error(
    reporting_checklist(fit, diagnostics = diag, bias_results = list(bad = data.frame(x = 1))),
    "`bias_results` in reporting_checklist\\(\\) must be NULL, output from estimate_bias\\(\\), an `mfrm_bias_collection`, or a list of `mfrm_bias` objects."
  )
})
