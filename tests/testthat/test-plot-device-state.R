test_that("single-panel plots advance through caller layouts and restore style", {
  fit <- make_toy_fit(maxit = 20)
  grDevices::pdf(NULL, width = 12, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), col = "purple", fg = "navy")
  style_names <- c("bg", "fg", "col", "col.axis", "col.lab", "col.main",
                   "col.sub", "cex.axis", "cex.lab", "cex.main", "lend", "ljoin", "mar")
  original <- graphics::par()[style_names]
  types <- c("pathway", "person", "step", "fit_pathway")
  positions <- matrix(c(1, 1, 1, 2, 2, 1, 2, 2), ncol = 2, byrow = TRUE)
  for (i in seq_along(types)) {
    .mfrmr_muffle_expected_warnings(
      plot(fit, type = types[i], preset = "publication"), "^Review-only display:"
    )
    expect_equal(graphics::par("mfg")[1:2], positions[i, ])
    expect_identical(graphics::par()[style_names], original)
    expect_false(graphics::par("new"))
    expect_false(identical(graphics::par("usr"), c(0, 1, 0, 1)))
  }

  graphics::layout(matrix(1:2, nrow = 1), widths = c(1, 3))
  .mfrmr_muffle_expected_warnings(plot(fit, type = "person"), "^Review-only display:")
  expect_equal(graphics::par("fig"), c(0, 0.25, 0, 1))
  .mfrmr_muffle_expected_warnings(
    plot(fit, type = "wright", renderer = "facets", show_ci = FALSE),
    "^Review-only display:"
  )
  expect_equal(graphics::par("fig"), c(0.25, 1, 0, 1))
})

test_that("Wright and CCC layouts do not shrink subsequent standalone plots", {
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (model in c("RSM", "PCM")) {
    fit <- make_toy_fit(model = model, maxit = 20)
    for (type in c("wright", "pathway", "ccc", "pathway")) {
      .mfrmr_muffle_expected_warnings(plot(fit, type = type), "^Review-only display:")
      expect_equal(graphics::par("mfrow"), c(1, 1))
      expect_false(graphics::par("new"))
      if (type == "pathway") expect_equal(graphics::par("fig"), c(0, 1, 0, 1))
    }
  }
})

test_that("APA composite draws four panels on one page", {
  fit <- make_toy_fit(maxit = 20)
  grDevices::pdf(NULL, width = 12, height = 9)
  on.exit(grDevices::dev.off(), add = TRUE)
  previous_hook <- getHook("plot.new")
  on.exit(setHook("plot.new", previous_hook, action = "replace"), add = TRUE)
  positions <- list()
  setHook("plot.new", function(...) {
    positions[[length(positions) + 1L]] <<- graphics::par("mfg")
  }, action = "append")
  result <- .mfrmr_muffle_expected_warnings(
    plot_apa_figure_one(fit, diagnostics = make_toy_diagnostics(fit)),
    "^Review-only display:"
  )
  expect_equal(do.call(rbind, positions),
               cbind(matrix(c(1, 1, 1, 2, 2, 1, 2, 2), ncol = 2, byrow = TRUE), 2, 2))
  expect_equal(graphics::par("mfrow"), c(1, 1))
  expect_identical(result$data$data$wright$data$wright_style, "facets_style")
})

test_that("native labels remain distinct without moving or discarding fitted data", {
  place_labels <- mfrmr:::.place_plot_labels
  layouts <- list()
  testthat::local_mocked_bindings(.place_plot_labels = function(...) {
    args <- list(...)
    placed <- place_labels(...)
    layouts[[length(layouts) + 1L]] <<- list(
      placed = placed, preferred_y = args[[2]], points = cbind(args$point_x, args$point_y),
      usr = graphics::par("usr")
    )
    placed
  }, .package = "mfrmr")
  for (size in list(c(7, 5), c(5, 4))) {
    grDevices::pdf(NULL, width = size[1], height = size[2])
    tryCatch({
      for (model in c("RSM", "PCM")) {
        fit <- make_toy_fit(model = model, maxit = 20)
        for (type in c("wright", "pathway")) {
          before <- .mfrmr_muffle_expected_warnings(
            plot(fit, type = type, draw = FALSE), "^Review-only display:"
          )
          drawn <- .mfrmr_muffle_expected_warnings(
            plot(fit, type = type), "^Review-only display:"
          )
          expect_identical(drawn, before)
          expected_n <- if (type == "wright") nrow(before$data$label_points) else nrow(before$data$steps)
          expect_equal(nrow(tail(layouts, 1)[[1]]$placed), expected_n)
        }
      }
      fit <- make_toy_fit(model = "PCM", maxit = 20)
      narrowed <- .mfrmr_muffle_expected_warnings(
        plot(fit, type = "pathway", theta_range = c(-0.1, 0.1)), "^Review-only display:"
      )
      expect_gt(nrow(narrowed$data$steps), nrow(tail(layouts, 1)[[1]]$placed))
      expect_true(all(abs(tail(layouts, 1)[[1]]$points[, 1]) <= 0.1))
      n_layouts <- length(layouts)
      .mfrmr_muffle_expected_warnings(
        plot(fit, type = "pathway", theta_range = c(5, 6)), "^Review-only display:"
      )
      expect_length(layouts, n_layouts)
    }, finally = grDevices::dev.off())
  }
  expect_length(layouts, 10L)
  for (layout in layouts) {
    box <- layout$placed
    overlap <- pmin(outer(box$Right, box$Right, pmin) - outer(box$Left, box$Left, pmax),
                    outer(box$Top, box$Top, pmin) - outer(box$Bottom, box$Bottom, pmax))
    expect_true(all(overlap[upper.tri(overlap)] < 1e-10))
    expect_true(all(box$Left >= layout$usr[1] & box$Right <= layout$usr[2]))
    expect_true(all(box$Bottom >= layout$usr[3] & box$Top <= layout$usr[4]))
    expect_lt(max(abs(box$Y - layout$preferred_y)), diff(layout$usr[3:4]) * 0.3)
    covered <- vapply(seq_len(nrow(layout$points)), function(j) {
      p <- layout$points[j, ]
      any(p[1] > box$Left & p[1] < box$Right & p[2] > box$Bottom & p[2] < box$Top)
    }, logical(1))
    expect_false(any(covered))
  }
})

test_that("long and Japanese names stay distinguishable when abbreviated", {
  for (prefix in c("Analytical judgement of evidence number ",
                   "論理構成と根拠の説明に関する評価観点番号")) {
    original <- paste0(prefix, sprintf("C%02d", 1:8))
    labels <- mfrmr:::truncate_axis_label(c(original, original[1], NA), width = 14)
    expect_equal(anyDuplicated(labels[1:8]), 0L)
    expect_true(all(nchar(labels[1:8], type = "width") <= 14))
    expect_equal(labels[1], labels[9])
    expect_true(is.na(labels[10]))
    expect_true(all(endsWith(labels[1:8], sprintf("C%02d", 1:8))))
  }
  original <- paste0("abcdefghijk", c("A", "B"), "lmnopqrstuv")
  expect_identical(mfrmr:::truncate_axis_label(original, width = 8), original)
  expect_equal(mfrmr:::truncate_axis_label(c("abcde...", "abcdefghijkl"), width = 8)[1], "abcde...")
})

test_that("numerous PCM thresholds do not erase the default facet display", {
  # A layout fixture must retain all nine categories in each PCM ladder.
  dat <- expand.grid(Person = paste0("P", 1:30), Rater = paste0("R", 1:2),
                     Criterion = paste0("C", 1:4), stringsAsFactors = FALSE)
  dat$Score <- rep(1:9, length.out = nrow(dat))
  expect_true(all(table(dat$Criterion, dat$Score) > 0))
  fit <- suppressWarnings(fit_mfrm(dat, person = "Person",
    facets = c("Rater", "Criterion"), score = "Score", model = "PCM",
    method = "JML", maxit = 15))
  full <- .mfrmr_muffle_expected_warnings(plot(fit, draw = FALSE), "^Review-only display:")
  compact <- .mfrmr_muffle_expected_warnings(plot(fit, top_n = 2, draw = FALSE), "^Review-only display:")
  expect_equal(sum(full$data$locations$PlotType == "Facet level"), 6)
  expect_equal(sum(compact$data$locations$PlotType == "Facet level"), 2)
  expect_equal(sum(full$data$locations$PlotType == "Step threshold"), 32)
  expect_equal(sum(compact$data$locations$PlotType == "Step threshold"), 32)
  original_legend <- graphics::legend
  legends <- list()
  testthat::local_mocked_bindings(legend = function(...) {
    args <- list(...)
    result <- original_legend(...)
    box <- result$rect
    legends[[length(legends) + 1L]] <<- list(labels = args$legend,
      x = graphics::grconvertX(c(box$left, box$left + box$w), from = "user", to = "ndc"),
      y = graphics::grconvertY(c(box$top - box$h, box$top), from = "user", to = "ndc"))
    result
  }, .package = "mfrmr")
  grDevices::pdf(NULL, width = 5, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  .mfrmr_muffle_expected_warnings(plot(fit, type = "ccc"), "^Review-only display:")
  expect_length(legends, 1L)
  expect_equal(legends[[1]]$labels, paste("Category", 1:9))
  expect_true(all(legends[[1]]$x >= 0 & legends[[1]]$x <= 1))
  expect_true(all(legends[[1]]$y >= 0 & legends[[1]]$y <= 1))
  expect_equal(length(unique(mfrmr:::.plot_series_colors(as.character(1:10)))), 10)
  rsm <- suppressWarnings(mfrmr::fit_mfrm(dat, person = "Person",
    facets = c("Rater", "Criterion"), score = "Score", model = "RSM",
    method = "JML", maxit = 15))
  grDevices::dev.off()
  grDevices::pdf(NULL, width = 10, height = 5)
  graphics::par(mfrow = c(1, 2))
  .mfrmr_muffle_expected_warnings(plot(rsm, type = "ccc"), "^Review-only display:")
  expect_equal(graphics::par("mfg")[1:2], c(1, 1))
  expect_true(all(legends[[2]]$x >= 0 & legends[[2]]$x <= 0.5))
  graphics::plot(1:2)
  expect_equal(graphics::par("mfg")[1:2], c(1, 2))
  legends <- list()
  small_pcm <- make_toy_fit(model = "PCM", maxit = 15)
  .mfrmr_muffle_expected_warnings(plot(small_pcm, type = "ccc"), "^Review-only display:")
  expect_length(legends, 1L)
  expect_lte(length(legends[[1]]$labels), 5L)
})

test_that("an overcrowded canvas reports its limit without deleting labels", {
  grDevices::pdf(NULL, width = 4, height = 3)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(0:1, 0:1, type = "n")
  expect_warning(placed <- mfrmr:::.place_plot_labels(
    rep(0.5, 40), rep(0.5, 40), sprintf("Label %02d", 1:40), cex = 0.7,
    point_x = 0.5, point_y = 0.5), "^Plot labels need more space")
  expect_equal(nrow(placed), 40)
  expect_true(all(is.finite(as.matrix(placed))))
})

test_that("pathway and CCC retain the reference profile outside the default figure", {
  fit <- make_toy_fit(maxit = 20)
  bundle <- .mfrmr_muffle_expected_warnings(
    plot(fit, type = "bundle", draw = FALSE), "^Review-only display:")
  expected_basis <- bundle$category_characteristic_curves$data$curve_basis
  expect_identical(bundle$pathway_map$data$curve_basis, expected_basis)
  labels <- character()
  original_mtext <- graphics::mtext
  testthat::local_mocked_bindings(mtext = function(text, ...) {
    labels <<- c(labels, as.character(text))
    original_mtext(text, ...)
  }, .package = "mfrmr")
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (type in c("pathway", "ccc")) {
    labels <- character()
    p <- .mfrmr_muffle_expected_warnings(plot(fit, type = type), "^Review-only display:")
    expect_identical(p$data$curve_basis, expected_basis)
    expect_false(any(grepl("additive effects and interactions fixed at zero", labels, fixed = TRUE)))
    expect_match(paste(plot_data(p, "notes")$Text, collapse = " "),
      "additive effects and interactions fixed at zero", fixed = TRUE)
    expect_true(all(p$data$curve_basis$PredictorOffset == 0))
  }
})

test_that("category-count axes include expected counts above observed bars", {
  grDevices::pdf(NULL, width = 12, height = 9)
  on.exit(grDevices::dev.off(), add = TRUE)
  axes <- list()
  original <- mfrmr:::barplot_rot45
  testthat::local_mocked_bindings(barplot_rot45 = function(..., main = NULL) {
    result <- original(..., main = main)
    if (main %in% c("Category counts", "Rating-scale category counts", "QC: Category counts")) {
      axes[[main]] <<- graphics::par("usr")[3:4]
    }
    result
  }, .package = "mfrmr")

  categories <- list(calculation_version = 2L, category_table = data.frame(
    Category = 1:3, Count = c(2, 4, 3), ExpectedCount = c(NA, 30, 1)
  ))
  categories$category_usage <- mfrmr:::summarize_category_usage(categories$category_table)
  categories$threshold_coverage <- mfrmr:::summarize_threshold_order(data.frame())
  empty <- categories
  empty$category_table$Count <- 0
  empty$category_table$ExpectedCount <- NULL
  empty$category_usage <- mfrmr:::summarize_category_usage(empty$category_table)
  for (draw in list(mfrmr:::draw_category_structure_bundle, mfrmr:::draw_rating_scale_bundle)) {
    draw(categories)
    expect_gt(tail(axes, 1)[[1]][2], 30)
    draw(empty)
    expect_gt(tail(axes, 1)[[1]][2], 0)
  }

  toy <- load_mfrmr_data("example_operational")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML")
  dashboard <- plot_qc_dashboard(fit, diagnostics = diagnose_mfrm(fit))
  counts <- dashboard$data$category_stats
  expect_gt(max(counts$ExpectedCount), max(counts$Count))
  expect_gt(axes[["QC: Category counts"]][2], max(counts$ExpectedCount))
})

test_that("rater profiles reserve room for complete labels or explain the size limit", {
  fit <- make_toy_fit(maxit = 20)
  dx <- make_toy_diagnostics(fit)
  # This layout test deliberately changes only display labels. Identity is
  # checked separately; bypass it here rather than fitting a new labelled model.
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity = function(...) "layout_fixture")
  rater <- dx$measures$Facet == "Rater"
  dx$measures$Level[rater] <- paste("Long assessment team", dx$measures$Level[rater])
  grDevices::pdf(NULL, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_mar <- graphics::par("mar")
  axis_calls <- 0L
  original_axis <- graphics::axis
  testthat::local_mocked_bindings(axis = function(side, at = NULL, labels = TRUE, ...) {
    if (side == 2 && is.character(labels)) {
      axis_calls <<- axis_calls + 1L
      expect_setequal(labels, dx$measures$Level[rater])
      expect_lt(max(graphics::strwidth(labels, units = "inches", cex = 0.85)),
                graphics::par("mai")[2] - 0.1)
    }
    original_axis(side, at = at, labels = labels, ...)
  }, .package = "mfrmr")
  plot_rater_severity_profile(fit, dx)
  expect_equal(axis_calls, 1L)
  expect_identical(graphics::par("mar"), old_mar)
  dx$measures$Level[rater] <- strrep("Very long label ", 30)
  expect_error(plot_rater_severity_profile(fit, dx), "wider graphics device")
  expect_identical(graphics::par("mar"), old_mar)
})

test_that("FACETS maps warn when adjacent multi-line cells cannot fit", {
  fit <- make_toy_fit(maxit = 20)
  p <- .mfrmr_muffle_expected_warnings(
    plot(fit, renderer = "facets", show_ci = FALSE, draw = FALSE), "^Review-only display:")
  grDevices::pdf(NULL, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  cells <- p$data$facets_style$facet_cells
  # Keep the observed locations; make the labels too large for those cells.
  crowded <- cells$Facet == cells$Facet[1]
  p$data$facets_style$facet_cells$CellLabel[crowded] <- paste(cells$CellLabel[crowded],
    strrep("long label ", 30))
  expect_warning(mfrmr:::draw_wright_facets_style(p$data, show_ci = FALSE),
                 "FACETS-style facet labels need more space")
})


test_that("facet bars reserve label space and include confidence limits", {
  grDevices::png(tempfile(fileext = ".png"), width = 1200, height = 740, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  before <- graphics::par("mar")
  seen <- NULL
  previous_hook <- getHook("plot.new")
  on.exit(setHook("plot.new", previous_hook, action = "replace"), add = TRUE)
  setHook("plot.new", function(...) {
    seen <<- graphics::par("mai")[2]
  }, action = "append")
  tab <- data.frame(Facet = "Rater", Level = c("R01", "R06"),
    Estimate = c(-.61, .41), SE = c(.3, .2))
  draw_facet_plot(tab, title = "", show_ci = TRUE)
  expect_gt(seen, max(graphics::strwidth(paste0(tab$Facet, ": ", tab$Level),
    units = "inches")))
  bounds <- c(tab$Estimate - qnorm(.975)*tab$SE, tab$Estimate + qnorm(.975)*tab$SE)
  expect_lte(graphics::par("usr")[1], min(bounds))
  expect_gte(graphics::par("usr")[2], max(bounds))
  expect_identical(graphics::par("mar"), before)
})
