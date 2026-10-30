test_that("plot APIs accept title/palette/label customization", {
  d <- mfrmr:::sample_mfrm_data(seed = 321)

  fit <- suppressWarnings(mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20,
    quad_points = 7
  ))
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  expect_no_error(
    .mfrmr_muffle_expected_warnings(
      plot(
        fit,
        type = "wright",
        draw = FALSE,
        title = "Custom Wright",
        palette = c(facet_level = "#1f78b4", step_threshold = "#d95f02"),
        label_angle = 45
      ),
      "^Review-only display:"
    )
  )

  expect_no_error(
    plot(
      fit, type = "unexpected",
      diagnostics = diag,
      view = "severity",
      draw = FALSE,
      main = "Custom Unexpected",
      palette = c(higher = "#d95f02", lower = "#1b9e77", bar = "#2b8cbe"),
      label_angle = 45,
      preset = "publication"
    )
  )

  expect_identical(
    as.character(
      plot(
        fit, type = "unexpected",
        diagnostics = diag,
        draw = FALSE,
        preset = "publication"
      )$data$preset
    ),
    "publication"
  )
  p_unexpected <- plot(
    fit, type = "unexpected",
    diagnostics = diag,
    draw = FALSE,
    preset = "publication"
  )
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_unexpected$data)))
  expect_true(is.data.frame(p_unexpected$data$legend))
  expect_true(is.data.frame(p_unexpected$data$reference_lines))

  p_displacement <- plot(
    fit, type = "displacement",
    diagnostics = diag,
    draw = FALSE,
    preset = "publication"
  )
  expect_identical(as.character(p_displacement$data$preset), "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_displacement$data)))

  p_fchi <- plot(
    fit, type = "facets_chisq",
    diagnostics = diag,
    draw = FALSE,
    preset = "publication"
  )
  expect_identical(as.character(p_fchi$data$preset), "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_fchi$data)))

  expect_no_error(
    {
      p_ir <- plot(
        fit, type = "agreement",
        diagnostics = diag,
        rater_facet = "Rater",
        metric = "exact",
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
    plot(
      fit, type = "agreement",
      diagnostics = diag,
      rater_facet = "Rater",
      metric = "exact",
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
    plot(
      t13_3, type = "bias",
      view = "facet_profile",
      draw = FALSE,
      main = "Custom Higher-Order Bias Profile",
      palette = c(normal = "#2b8cbe", flag = "#cb181d", profile = "#756bb1"),
      label_angle = 45,
      preset = "publication"
    )
  )

  p_pca <- plot(
    mfrmr::analyze_residual_pca(diag, mode = "overall"), type = "residual_pca",
    mode = "overall",
    view = "scree",
    draw = FALSE,
    preset = "publication"
  )
  expect_identical(as.character(p_pca$data$preset), "publication")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_pca$data)))

  p_wright <- .mfrmr_muffle_expected_warnings(
    plot(fit, type = "wright", draw = FALSE, preset = "publication"),
    "^Review-only display:"
  )
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_wright$data)))

  p_bubble_mono <- plot(fit, type = "bubble", draw = FALSE, preset = "monochrome")
  expect_identical(as.character(p_bubble_mono$data$preset), "monochrome")
})

test_that("fit plot presentation flags retain notes and numerical payloads", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  fit <- make_toy_fit(model = "PCM", maxit = 20)
  grDevices::pdf(NULL, width = 12, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  text_drawn <- character()
  original_mtext <- graphics::mtext
  original_title <- graphics::title
  testthat::local_mocked_bindings(
    mtext = function(text, ...) {
      text_drawn <<- c(text_drawn, as.character(text))
      original_mtext(text, ...)
    },
    title = function(main = NULL, ...) {
      text_drawn <<- c(text_drawn, as.character(main))
      original_title(main = main, ...)
    }, .package = "mfrmr")
  local_graphics_headings(function(main) text_drawn <<- c(text_drawn, as.character(main)))
  for (type in c("wright", "facets", "pathway", "ccc", "fit_pathway", "person", "step", "facet")) {
    args <- list(x = fit, type = if (type == "facets") "wright" else type,
                 title = "Custom heading", top_n = 2)
    if (type == "facets") {
      args$renderer <- "facets"
      args$persons_per_star <- 10
    }
    baseline <- .mfrmr_muffle_expected_warnings(
      do.call(plot, c(args, list(draw = FALSE))), "^Review-only display:")
    text_drawn <- character()
    p <- .mfrmr_muffle_expected_warnings(
      do.call(plot, c(args, list(show_title = FALSE, show_notes = FALSE))),
      "^Review-only display:")
    expect_false(any(grepl("Custom heading|REVIEW ONLY|fixed at zero|Compact native view|rows/logit", text_drawn)))
    expect_identical(p$data$notes, baseline$data$notes)
    expect_identical(p$data$fit_readiness, baseline$data$fit_readiness)
    expect_identical(p$data$interpretation_status, baseline$data$interpretation_status)
    expect_identical(p$data[setdiff(names(p$data), "display")],
                     baseline$data[setdiff(names(baseline$data), "display")])
    g <- as_ggplot(p)
    expect_null(g$labels$title)
    expect_null(g$labels$subtitle)
    expect_null(g$labels$caption)
    expect_identical(utils::head(plot_data(g, "notes"), nrow(p$data$notes)), p$data$notes)
    expect_true(any(grepl("Review-only display", capture.output(print(p)), fixed = TRUE)))
    if (type == "wright") expect_true(any(grepl("omitted", p$data$notes$Text, fixed = TRUE)))
    if (type %in% c("pathway", "ccc")) expect_true(any(p$data$notes$Type == "reference_profile"))
    if (type == "facets") expect_true(all(p$data$facets_style$headers$Header %in%
      gsub("\n", "", text_drawn, fixed = TRUE)))
  }
  for (flags in list(c(TRUE, FALSE), c(FALSE, TRUE))) {
    text_drawn <- character()
    p <- .mfrmr_muffle_expected_warnings(
      plot(fit, type = "pathway", title = "Independent heading",
           show_title = flags[1], show_notes = flags[2]), "^Review-only display:")
    expect_identical(any(grepl("Independent heading", text_drawn, fixed = TRUE)), flags[1])
    expect_identical(any(grepl("fixed at zero", text_drawn, fixed = TRUE)), flags[2])
    g <- as_ggplot(p)
    expect_identical(is.null(g$labels$title), !flags[1])
    expect_identical(is.null(g$labels$subtitle), !flags[2])
  }
  text_drawn <- character()
  .mfrmr_muffle_expected_warnings(
    plot(fit, renderer = "facets", title = "Optional FACETS heading",
      show_title = TRUE, show_notes = FALSE), "^Review-only display:")
  expect_true(any(grepl("Optional FACETS heading", text_drawn, fixed = TRUE)))
  direct <- .mfrmr_muffle_expected_warnings(
    as_ggplot(fit, type = "wright", show_title = FALSE, show_notes = FALSE),
    "^Review-only display:")
  expect_null(direct$labels$title)
  expect_true(is.data.frame(attr(direct, "mfrmr_notes")))
  expect_warning(.mfrmr_muffle_expected_warnings(
    plot(fit, type = "wright", renderer = "facets", persons_per_star = 0.1,
         show_title = FALSE, show_notes = FALSE), "^Review-only display:"),
    "FACETS-style frequency stars exceed their column", fixed = TRUE)
  expect_warning(plot(fit, draw = FALSE, show_title = FALSE, show_notes = FALSE),
                 "^Review-only display:")
  for (flag in c("show_title", "show_notes")) {
    for (bad in list(NA, NULL, "no", 1, c(TRUE, FALSE))) {
      expect_error(do.call(plot, c(list(fit, draw = FALSE), setNames(list(bad), flag))),
                   paste0("`", flag, "` must be TRUE or FALSE"), fixed = TRUE)
    }
  }
  bundle <- .mfrmr_muffle_expected_warnings(
    plot(fit, type = "bundle", draw = FALSE, show_title = FALSE, show_notes = FALSE),
    "^Review-only display:")
  expect_true(all(vapply(bundle, function(p) identical(p$data$display,
    list(show_title = FALSE, show_notes = FALSE)), logical(1))))
})
