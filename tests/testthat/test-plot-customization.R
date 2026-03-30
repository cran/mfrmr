test_that("plot APIs accept title/palette/label customization", {
  d <- mfrmr:::sample_mfrm_data(seed = 321)

  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20,
    quad_points = 7
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  expect_no_error(
    plot(
      fit,
      type = "wright",
      draw = FALSE,
      title = "Custom Wright",
      palette = c(facet_level = "#1f78b4", step_threshold = "#d95f02"),
      label_angle = 45
    )
  )

  expect_no_error(
    mfrmr::plot_unexpected(
      fit,
      diagnostics = diag,
      plot_type = "severity",
      draw = FALSE,
      main = "Custom Unexpected",
      palette = c(higher = "#d95f02", lower = "#1b9e77", bar = "#2b8cbe"),
      label_angle = 45,
      preset = "publication"
    )
  )

  expect_identical(
    as.character(
      mfrmr::plot_unexpected(
        fit,
        diagnostics = diag,
        draw = FALSE,
        preset = "publication"
      )$data$preset
    ),
    "publication"
  )
  p_unexpected <- mfrmr::plot_unexpected(
    fit,
    diagnostics = diag,
    draw = FALSE,
    preset = "publication"
  )
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_unexpected$data)))
  expect_true(is.data.frame(p_unexpected$data$legend))
  expect_true(is.data.frame(p_unexpected$data$reference_lines))

  p_displacement <- mfrmr::plot_displacement(
    fit,
    diagnostics = diag,
    draw = FALSE,
    preset = "publication"
  )
  expect_identical(as.character(p_displacement$data$preset), "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_displacement$data)))

  p_fchi <- mfrmr::plot_facets_chisq(
    fit,
    diagnostics = diag,
    draw = FALSE,
    preset = "publication"
  )
  expect_identical(as.character(p_fchi$data$preset), "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_fchi$data)))

  expect_no_error(
    {
      p_ir <- mfrmr::plot_interrater_agreement(
        fit,
        diagnostics = diag,
        rater_facet = "Rater",
        plot_type = "exact",
        draw = FALSE,
        main = "Custom Inter-rater",
        palette = c(ok = "#2b8cbe", flag = "#cb181d", expected = "#08519c"),
        label_angle = 45,
        preset = "publication"
      )
      stopifnot(identical(as.character(p_ir$data$preset), "publication"))
      stopifnot(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_ir$data)))
    }
  )

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    mfrmr::plot_interrater_agreement(
      fit,
      diagnostics = diag,
      rater_facet = "Rater",
      plot_type = "exact",
      draw = TRUE,
      main = "Custom Inter-rater",
      palette = c(ok = "#2b8cbe", flag = "#cb181d", expected = "#08519c"),
      label_angle = 45,
      preset = "publication"
    )
  )

  bias3 <- mfrmr::estimate_bias(
    fit,
    diag,
    interaction_facets = c("Rater", "Task", "Criterion"),
    max_iter = 2
  )
  t13_3 <- mfrmr::bias_interaction_report(bias3, top_n = 20)

  expect_no_error(
    mfrmr::plot_bias_interaction(
      t13_3,
      plot = "facet_profile",
      draw = FALSE,
      main = "Custom Higher-Order Bias Profile",
      palette = c(normal = "#2b8cbe", flag = "#cb181d", profile = "#756bb1"),
      label_angle = 45,
      preset = "publication"
    )
  )

  p_pca <- mfrmr::plot_residual_pca(
    mfrmr::analyze_residual_pca(diag, mode = "overall"),
    mode = "overall",
    plot_type = "scree",
    draw = FALSE,
    preset = "publication"
  )
  expect_identical(as.character(p_pca$data$preset), "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_pca$data)))

  p_wright <- plot(fit, type = "wright", draw = FALSE, preset = "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_wright$data)))
})
