plot_title_cases <- function(fit, diagnostics) {
  bias <- bias_interaction_report(fit, diagnostics = diagnostics,
    facet_a = "Rater", facet_b = "Criterion")
  dashboard <- facet_quality_dashboard(fit, diagnostics = diagnostics, facet = "Rater")
  definitions <- list(
    plot_marginal_fit = list(x = fit, arg = "plot_type", views = c("std_residual", "prop_diff")),
    plot_marginal_pairwise = list(x = fit, arg = "metric", views = c("exact", "adjacent")),
    plot_unexpected = list(x = fit, arg = "plot_type", views = c("scatter", "severity")),
    plot_interrater_agreement = list(x = fit, arg = "plot_type", views = c("exact", "corr", "difference")),
    plot_facets_chisq = list(x = fit, arg = "plot_type", views = c("fixed", "random", "variance")),
    plot_bubble = list(x = fit, arg = "view", views = c("measure", "infit_outfit")),
    plot_bias_interaction = list(x = bias, arg = "plot", views = c("scatter", "ranked", "heatmap", "abs_t_hist", "facet_profile")),
    plot_facet_quality_dashboard = list(x = dashboard, arg = "plot_type", views = c("severity", "flags")))
  out <- list()
  for (fun in names(definitions)) {
    d <- definitions[[fun]]
    for (view in d$views) {
      args <- list(x = d$x, draw = FALSE)
      if (inherits(d$x, "mfrm_fit")) args$diagnostics <- diagnostics
      args[[d$arg]] <- view
      out[[paste(fun, view, sep = "/")]] <- list(fun = fun, args = args)
    }
  }
  out
}

local({
  fit <- make_toy_fit(method = "MML", maxit = 500)
  diagnostics <- make_toy_diagnostics(fit, diagnostic_mode = "both")
  cases <- plot_title_cases(fit, diagnostics)

  test_that("title aliases preserve all 21 numerical payloads and legacy defaults", {
    for (case in cases) {
      old <- do.call(case$fun, c(case$args, list(main = "Custom plot title")))
      new <- do.call(case$fun, c(case$args, list(title = "Custom plot title")))
      expect_identical(new, old, info = case$fun)
      original <- do.call(case$fun, case$args)
      expect_identical(do.call(case$fun, c(case$args, list(main = NULL))), original)
      blank <- do.call(case$fun, c(case$args, list(title = NULL)))
      expect_identical(blank$data$title, "", info = case$fun)
      expect_identical(do.call(case$fun, c(case$args, list(title = ""))), blank)
      blank$data$title <- original$data$title
      expect_identical(blank, original, info = case$fun)
    }
  })

  test_that("conflicting and invalid titles fail before consuming data or drawing", {
    for (fun in unique(vapply(cases, `[[`, "", "fun"))) {
      for (title in list(NULL, "same")) {
        expect_error(do.call(fun, list(x = NULL, main = title, title = title)), "only one")
      }
      for (bad in list(NA_character_, c("A", "B"), 1, character(0))) {
        expect_error(do.call(fun, list(x = NULL, title = bad)), "one character string")
      }
    }
  })

  test_that("base rendering omits headings while preserving selected title metadata", {
    seen <- list()
    real_title <- graphics::title
    local_mocked_bindings(title = function(main = NULL, ...) {
      seen[[length(seen) + 1L]] <<- list(main = main)
      real_title(main = main, ...)
    }, .package = "mfrmr")
    local_graphics_headings(function(main) seen[[length(seen) + 1L]] <<- list(main = main))
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (case in cases) {
      args <- case$args; args$draw <- TRUE
      for (text in list("Heading sent to renderer", NULL)) {
        seen <- list()
        p <- do.call(case$fun, c(args, list(title = text)))
        expect_true(length(seen) > 0L, info = paste(case$fun, args[[2]]))
        expect_false(any(vapply(seen, function(s) any(nzchar(as.character(s$main))), logical(1))),
          info = paste(case$fun, args[[2]]))
      }
    }
  })

  test_that("S3 dashboard forwarding and supported ggplot conversion retain suppression", {
    dashboard <- cases[["plot_facet_quality_dashboard/severity"]]$args$x
    dashboard$interpretation_status <- "review_only"
    hidden <- plot(dashboard, title = NULL, draw = FALSE)
    original <- plot(dashboard, draw = FALSE)
    expect_identical(hidden$data$title, "")
    expect_match(original$data$title, "REVIEW ONLY")
    expect_identical(hidden$data$notes, original$data$notes)
    expect_identical(hidden$data$fit_readiness, original$data$fit_readiness)
    expect_identical(hidden$data$interpretation_status, "review_only")
    skip_if_not_installed("ggplot2")
    args <- cases[["plot_bubble/measure"]]$args
    p <- do.call(plot_bubble, c(args, list(title = NULL)))
    g <- as_ggplot(p)
    expect_null(g$labels$title)
    expect_no_error(ggplot2::ggplot_build(g))
    custom <- do.call(plot_bubble, c(args, list(title = "Custom ggplot heading")))
    expect_identical(plot_data(as_ggplot(custom), "title"), "Custom ggplot heading")
    file <- tempfile(fileext = ".rds")
    on.exit(unlink(file), add = TRUE)
    saveRDS(p, file)
    expect_null(as_ggplot(readRDS(file))$labels$title)
  })
})

test_that("bubble conversion preserves saved radii, facet colours, order and reference lines", {
  skip_if_not_installed("ggplot2")
  # Unequal N and SE make a constant-size or area/radius swap detectable.
  d <- structure(list(measures = data.frame(Facet = c("Rater", "Criterion", "Rater"),
    Level = c("A", "B", "C"), Estimate = c(-.4, .2, .7),
    Infit = c(.8, 1.1, 1.4), Outfit = c(.9, 1.3, 1.2),
    N = c(9, 36, 81), SE = c(.1, .2, .4))), class = "mfrm_diagnostics")
  for (size in c("N", "SE", "equal")) {
    for (view in c("measure", "infit_outfit")) {
      p <- plot_bubble(d, view = view, bubble_size = size, preset = "monochrome", draw = FALSE)
      g <- as_ggplot(p)
      built <- ggplot2::ggplot_build(g)$data
      point <- built[[1]]
      expect_equal(point$size, 40 * p$data$radius)
      expect_equal(point$size / point$size[1], p$data$radius / p$data$radius[1])
      expected <- switch(size, N = sqrt(p$data$table$N), SE = 1 / p$data$table$SE,
        equal = rep(1, nrow(p$data$table)))
      expect_equal(point$size / point$size[1], expected / expected[1])
      expect_equal(point$x, if(view == "measure") p$data$table$Estimate else p$data$table$Infit)
      expect_equal(point$y, if(view == "measure") p$data$table$Infit else p$data$table$Outfit)
      colours <- setNames(p$data$legend$value, p$data$legend$label)
      expect_identical(point$colour, unname(colours[as.character(p$data$table$Facet)]))
      rgb <- grDevices::col2rgb(point$colour)
      expect_true(all(rgb[1, ] == rgb[2, ] & rgb[2, ] == rgb[3, ]))
      expect_identical(levels(g$data$Facet), unique(as.character(p$data$table$Facet)))
      reference <- p$data$reference_lines
      rendered_h <- unlist(lapply(built[-1], function(x) x$yintercept))
      rendered_v <- unlist(lapply(built[-1], function(x) x$xintercept))
      expect_equal(sort(as.numeric(rendered_h)), sort(reference$value[reference$axis == "h"]))
      expect_equal(sort(as.numeric(rendered_v)), sort(reference$value[reference$axis == "v"]))
    }
  }
  p <- plot_bubble(d, palette = c(Rater = "red", Criterion = "blue"), draw = FALSE)
  expect_equal(ggplot2::ggplot_build(as_ggplot(p))$data[[1]]$colour, c("red", "blue", "red"))
  # Physical device units differ, but save/load must not change the mapping.
  expect_identical(as_ggplot(unserialize(serialize(p, NULL)))$data, as_ggplot(p)$data)
  p$data$radius <- NULL
  expect_error(as_ggplot(p), "Saved bubble radii")
})
