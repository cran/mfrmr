test_that("all direct common-preset entries consult the option only when omitted", {
  ns <- asNamespace("mfrmr")
  public <- getNamespaceExports("mfrmr")
  entries <- Filter(function(name) {
    f <- get(name, ns)
    is.function(f) && (name %in% public || startsWith(name, "plot.")) &&
      "preset" %in% names(formals(f)) &&
      grepl('"standard"', paste(deparse(formals(f)$preset), collapse = " "), fixed = TRUE)
  }, ls(ns, all.names = TRUE))
  expect_true(all(c("plot.mfrm_fit", "plot.mfrm_multivariate_d_study") %in% entries))
  withr::local_options(mfrmr.plot_preset = "not-a-preset")
  for (name in entries) {
    # NULL data allow model-routing guards to inspect the argument. This
    # checks preset handling independently of missing-argument error order.
    fun <- get(name, ns)
    required <- names(formals(fun))[vapply(formals(fun), identical, logical(1), quote(expr = ))]
    required <- setdiff(required, "...")
    args <- setNames(rep(list(NULL), length(required)), required)
    expect_error(do.call(fun, args), "Option `mfrmr.plot_preset`", info = name)
    error <- tryCatch(do.call(fun, list(preset = "standard")), error = identity)
    expect_false(inherits(error, "error") &&
      grepl("mfrmr.plot_preset", conditionMessage(error), fixed = TRUE), info = name)
  }
  for (bad in list(NA_character_, character(), c("standard", "monochrome"), 1, "pub")) {
    options(mfrmr.plot_preset = bad)
    expect_error(plot_bubble(NULL, draw = FALSE), "Option `mfrmr.plot_preset`")
  }
})

local({
  fit <- make_toy_fit(method = "MML", maxit = 500)
  diagnostics <- make_toy_diagnostics(fit, diagnostic_mode = "both")
  features <- mfrm_features(data.frame(ID = letters[1:6], Years = c(1, 2, 5, 7, 10, 12),
    Training = c(4, 8, 6, 10, 16, 14)), "ID", c("Years", "Training"))
  cases <- list(
    bubble = list(fun = plot_bubble, args = list(x = fit, diagnostics = diagnostics)),
    wright = list(fun = plot, args = list(x = fit, type = "wright", show_ci = FALSE)),
    curves = list(fun = plot, args = list(x = fit, type = "ccc")),
    paired = list(fun = plot_compare_mfrm, args = list(reference = fit, comparison = fit, type = "ccc")),
    pca = list(fun = plot, args = list(x = mfrm_pca(features), type = "scores")),
    tree = list(fun = plot, args = list(x = mfrm_cluster_hierarchical(features, k = 2))),
    quality = list(fun = plot, args = list(x = data_quality_report(fit), type = "category_counts")),
    coverage = list(fun = plot, args = list(x = subset_connectivity_report(fit, diagnostics = diagnostics), type = "design_matrix")))
  run <- function(case, ...) do.call(case$fun, c(case$args, list(draw = FALSE), list(...)))

  test_that("session defaults match explicit presets, preserve overrides and restore defaults", {
    withr::local_options(mfrmr.plot_preset = NULL)
    standard <- lapply(cases, run)
    for (preset in c("standard", "publication", "compact", "monochrome")) {
      options(mfrmr.plot_preset = preset)
      for (name in names(cases)) {
        actual <- run(cases[[name]])
        expect_identical(actual, run(cases[[name]], preset = preset), info = name)
        expect_identical(actual$data$preset, preset, info = name)
        expect_identical(run(cases[[name]], preset = "standard"), standard[[name]], info = name)
        expect_identical(run(cases[[name]], preset = NULL), standard[[name]], info = name)
      }
    }
    # One explicit parent choice must reach both source plots, even if the
    # session option is invalid. An ignored argument must not read the option.
    options(mfrmr.plot_preset = "not-a-preset")
    paired <- run(cases$paired, preset = "publication")
    expect_true(all(vapply(paired$data$source_plots, function(x)
      identical(x$data$preset, "publication"), logical(1))))
    for (name in names(cases)) expect_identical(run(cases[[name]], preset = "standard"), standard[[name]])
    options(mfrmr.plot_preset = NULL)
    expect_identical(lapply(cases, run), standard)
  })

  test_that("saved payload rendering does not read a changed session option or alter global theme", {
    skip_if_not_installed("ggplot2")
    withr::local_options(mfrmr.plot_preset = "monochrome")
    theme <- ggplot2::theme_get()
    p <- run(cases$bubble)
    g <- as_ggplot(p)
    built <- ggplot2::ggplot_build(g)$data
    file <- tempfile(fileext = ".rds")
    on.exit(unlink(file), add = TRUE)
    saveRDS(p, file)
    for (preset in c("publication", "not-a-preset")) {
      options(mfrmr.plot_preset = preset)
      expect_identical(ggplot2::ggplot_build(as_ggplot(readRDS(file)))$data, built)
      expect_identical(ggplot2::theme_get(), theme)
    }
    options(mfrmr.plot_preset = "publication")
    newer <- run(cases$bubble)
    expect_identical(newer$data$preset, "publication")
    expect_identical(newer$data$table, p$data$table)
    expect_identical(newer$data$radius, p$data$radius)
    expect_identical(newer$data$reference_lines, p$data$reference_lines)
    # Legacy renderer defaults are stable when a saved payload lacks a preset.
    old <- run(cases$curves, preset = "standard")
    old$data$preset <- NULL
    options(mfrmr.plot_preset = "not-a-preset")
    expect_no_error(ggplot2::ggplot_build(as_ggplot(old)))
  })

  test_that("results and plot-data forwarding resolve the same default without estimation", {
    withr::local_options(mfrmr.plot_preset = "monochrome")
    res <- mfrm_results(fit, diagnostics = diagnostics, include = c("fit", "plots"), compute = "never")
    local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"),
      diagnose_mfrm = function(...) stop("Unexpected diagnostic calculation"))
    p <- plot(res, type = "wright", show_ci = FALSE, draw = FALSE)
    expect_identical(p$data$preset, "monochrome")
    expect_identical(p, plot(res, type = "wright", preset = "monochrome", show_ci = FALSE, draw = FALSE))
    tables <- plot_data(fit, type = "ccc")
    expected <- plot_data(fit, type = "ccc", preset = "monochrome")
    expect_identical(tables, expected)
  })
})
