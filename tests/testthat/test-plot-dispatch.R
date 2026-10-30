plot_dispatch_fixture <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      fit <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[[1L]]
      d <- diagnose_mfrm(fit, residual_pca = "both", diagnostic_mode = "both")
      saved <<- list(fit = fit, diagnostics = d)
    }
    saved
  }
})

test_that("only data accessors retain the public plot_ prefix", {
  expect_setequal(grep("^plot_", getNamespaceExports("mfrmr"), value = TRUE),
    c("plot_data", "plot_data_components"))
  expect_error(getExportedValue("mfrmr", "plot_reliability_snapshot"), "not an exported")
})

test_that("consolidated fit views preserve renderer defaults and saved notes", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit; d <- fixture$diagnostics
  routes <- list(
    reliability = "plot_reliability_snapshot", rater_severity = "plot_rater_severity_profile",
    thresholds = "plot_threshold_ladder", person_fit = "plot_person_fit",
    residual_matrix = "plot_residual_matrix", residual_qq = "plot_residual_qq",
    guttman = "plot_guttman_scalogram", marginal_fit = "plot_marginal_fit",
    marginal_pairwise = "plot_marginal_pairwise", bubble = "plot_bubble",
    displacement = "plot_displacement", facets_chisq = "plot_facets_chisq",
    fair_average = "plot_fair_average", unexpected = "plot_unexpected",
    local_dependence = "plot_local_dependence_heatmap", qc = "plot_qc_dashboard",
    apa_summary = "plot_apa_figure_one", facet_dashboard = "plot_facet_quality_dashboard",
    facet_equivalence = "plot_facet_equivalence")
  for (type in names(routes)) {
    old <- get(routes[[type]], asNamespace("mfrmr"))
    args <- list(draw = FALSE)
    if ("diagnostics" %in% names(formals(old))) args$diagnostics <- d
    expected <- suppressWarnings(do.call(old, c(list(fit), args)))
    actual <- suppressWarnings(do.call(plot, c(list(fit, type = type), args)))
    expect_identical(actual, expected, info = type)
    if (inherits(actual, "mfrm_plot_data"))
      expect_identical(plot_data(actual, "notes"), plot_data(expected, "notes"), info = type)
  }
  expect_identical(plot(fit, type = "wright", view = "unified", draw = FALSE),
    plot_wright_unified(fit, draw = FALSE))
  shrunk <- apply_empirical_bayes_shrinkage(fit)
  expect_identical(plot(shrunk, type = "shrinkage", view = "funnel", draw = FALSE),
    plot_shrinkage_funnel(shrunk, draw = FALSE))
  expect_identical(plot(fit, type = "agreement", view = "heatmap", metric = "correlation",
    diagnostics = d, draw = FALSE),
    plot_rater_agreement_heatmap(fit, metric = "correlation", diagnostics = d, draw = FALSE))
  expect_identical(plot(fit, type = "agreement", metric = "correlation", diagnostics = d, draw = FALSE),
    plot_interrater_agreement(fit, plot_type = "corr", diagnostics = d, draw = FALSE))
  expect_identical(plot(fit, type = "reliability", metric = "strata", diagnostics = d, draw = FALSE),
    plot_reliability_snapshot(fit, metric = "strata", diagnostics = d, draw = FALSE))
})

test_that("Wright dispatch retains named levels without recalculation", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit
  level <- tail(as.character(fit$facets$others$Level[fit$facets$others$Facet == "Rater"]), 1)
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    diagnose_mfrm = function(...) stop("No new diagnostics"), .package = "mfrmr")
  for (style in c("native", "facets_style")) {
    p <- plot(fit, top_n = 1, highlight = list(Rater = level),
      wright_style = style, draw = FALSE)
    marked <- subset(p$data$locations, Highlighted)
    expect_identical(as.character(marked$Label), level)
    expect_identical(as.character(marked$Group), "Rater")
    expect_match(paste(plot_data(p, "notes")$Text, collapse = " "), "not a statistical flag")
  }
})

test_that("saved result dispatch avoids refitting and preserves result semantics", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit; d <- fixture$diagnostics
  local_mocked_bindings(fit_mfrm = function(...) stop("No model fitting while plotting"), .package = "mfrmr")
  info <- compute_information(fit)
  expect_identical(plot(info, type = "iif", facet = "Rater", draw = FALSE),
    plot_information(info, type = "iif", facet = "Rater", draw = FALSE))
  unavailable_pca <- analyze_residual_pca(d)
  expect_error(plot(unavailable_pca, type = "scree", draw = FALSE), "No eigenvalues")
  obs <- expand.grid(Person = 1:20, Criterion = letters[1:4])
  obs$StdResidual <- c(sin(1:20), cos(1:20), sin(2 * (1:20)), c(-10:-1, 10:1) / 10)
  pca <- analyze_residual_pca(list(obs = obs, facet_names = "Criterion"), mode = "both")
  expect_identical(plot(pca, type = "loadings", component = 1L, draw = FALSE),
    plot_residual_pca(pca, plot_type = "loadings", component = 1L, draw = FALSE))
  expect_identical(plot(pca, type = "overall_scree", draw = FALSE),
    plot(pca, type = "scree", draw = FALSE))
  qc <- run_qc_pipeline(fit, diagnostics = d, include_bias = FALSE)
  expect_identical(plot(qc, draw = FALSE), plot_qc_pipeline(qc, draw = FALSE))
  expect_warning(drift <- detect_anchor_drift(list(A = fit, B = fit)), "Thin linking support")
  expect_identical(plot(drift, type = "forest", draw = FALSE),
    plot_anchor_drift(drift, type = "forest", draw = FALSE))
  bias <- estimate_bias(fit, d, facet_a = "Rater", facet_b = "Criterion")
  expect_identical(plot(bias, type = "heatmap", draw = FALSE),
    plot_bias_interaction(bias, plot = "heatmap", draw = FALSE))
  expect_identical(plot(d, type = "marginal_fit", view = "prop_diff", draw = FALSE),
    plot_marginal_fit(d, plot_type = "prop_diff", draw = FALSE))
})

test_that("DIF and interaction result methods keep their distinct metrics", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit; d <- fixture$diagnostics
  data <- fit$prep$data
  people <- unique(data$Person)
  data$Group <- ifelse(match(data$Person, people) <= length(people)/2, "A", "B")
  dif <- analyze_dff(fit, d, facet = "Criterion", group = "Group", data = data, min_obs = 2)
  expect_identical(plot(dif, draw = FALSE), plot_dif_summary(dif, draw = FALSE))
  expect_identical(plot(dif, type = "heatmap", metric = "contrast", draw = FALSE),
    plot_dif_heatmap(dif, metric = "contrast", draw = FALSE))
  interaction <- dif_interaction_table(fit, d, facet = "Criterion", group = "Group", data = data, min_obs = 2)
  expect_identical(plot(interaction, draw = FALSE), plot_dif_heatmap(interaction, draw = FALSE))
  expect_error(plot(interaction, metric = "contrast", draw = FALSE), "requires")
})

test_that("multi-fit plots retain order and do not turn comparisons into inference", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit
  expect_identical(suppressWarnings(plot(fit, comparison = fit, type = "ccc", view = "difference", draw = FALSE)),
    suppressWarnings(plot_compare_mfrm(fit, fit, type = "ccc", view = "difference", draw = FALSE)))
  expect_identical(plot(fit, type = "trajectory", comparison = list(T2 = fit), labels = c("T1", "T2"), draw = FALSE),
    plot_rater_trajectory(list(T1 = fit, T2 = fit), draw = FALSE))
  expect_error(plot(fit, type = "trajectory", draw = FALSE), "require")
  expect_error(plot(fit, type = "trajectory", comparison = list(fit), labels = c("T1", "T1")), "distinct")
  expect_error(plot(fit, type = "reliability", comparison = fit), "only")
})

test_that("view-specific controls cannot disappear into dots", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit; d <- fixture$diagnostics
  expect_error(plot(fit, type = "reliability", diagnostics = d, typo = 1), "not supported")
  expect_error(plot(fit, type = "agreement", view = "unknown"), "pairs.*heatmap")
  expect_error(plot(fit, type = "thresholds", metric = "strata"), "not supported")
  expect_error(plot(fit, type = "marginal_fit", view = "prop_diff", plot_type = "std_residual"), "only one")
  expect_error(plot(fit, type = "rater_severity", level = .9, ci_level = .95), "only one")
  expect_error(plot(fit, type = "ccc", view = "heatmap"), "not supported")
  expect_error(plot(fit, type = "fair_average", diagnostics = d, typo = 1), "not supported")
  expect_error(plot(fit, type = "fair_average", diagnostics = d, fit = fit),
    "matches multiple formal arguments|not supported")
  expect_error(plot(fit, type = "reliability", diagnostics = d, totalscore = FALSE), "not supported")
  expect_error(do.call(plot, c(list(fit, type = "fair_average", diagnostics = d),
    list(totalscore = TRUE, totalscore = FALSE))), "Name each plot argument once")
  product <- fit
  product$config$gpcm_spec <- list(component_levels = list(Rater = "1", Criterion = "1"))
  expect_error(plot(product, type = "reliability", diagnostics = d), "two-family|product|two slope")
})

test_that("saved table bundles and dedicated results accept their documented family types", {
  fixture <- plot_dispatch_fixture(); fit <- fixture$fit; d <- fixture$diagnostics
  cases <- list(
    unexpected = list(unexpected_response_table(fit, d), plot_unexpected),
    fair_average = list(fair_average_table(fit, d), plot_fair_average),
    displacement = list(displacement_table(fit, d), plot_displacement),
    agreement = list(interrater_agreement_table(fit, d), plot_interrater_agreement),
    facets_chisq = list(facets_chisq_table(fit, d), plot_facets_chisq),
    facet_equivalence = list(suppressWarnings(analyze_facet_equivalence(fit)), plot_facet_equivalence),
    facet_dashboard = list(facet_quality_dashboard(fit, d), plot_facet_quality_dashboard))
  for (type in names(cases)) {
    x <- cases[[type]][[1]]; old <- cases[[type]][[2]]
    expect_identical(suppressWarnings(plot(x, type = type, draw = FALSE)),
      suppressWarnings(old(x, draw = FALSE)), info = type)
  }
  expect_identical(plot(cases$agreement[[1]], type = "agreement", metric = "correlation", draw = FALSE),
    plot_interrater_agreement(cases$agreement[[1]], plot_type = "corr", draw = FALSE))
})

test_that("residual heatmaps retain their facet, correlation scale and unavailable cells", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  mat <- matrix(c(1, -.4, NA, -.4, 1, 0, NA, 0, 1), 3, 3,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  saved <- new_mfrm_plot_data("local_dependence_heatmap", list(
    matrix = mat, facet = "Criterion", title = "Residual associations",
    caption = "Unavailable pairs are not zero correlations."))
  g <- as_ggplot(saved)
  expect_identical(g$labels$x, "Criterion")
  expect_identical(g$labels$y, "Criterion")
  expect_identical(g$scales$get_scales("fill")$name, "Residual correlation")
  expect_identical(g$scales$get_scales("fill")$limits, c(-1, 1))
  expect_equal(g$data$Value, as.vector(mat))
  built <- ggplot2::ggplot_build(g)$data[[1]]
  expect_equal(nrow(built), 9L)
  expect_true(all(built$fill[is.na(g$data$Value)] == "grey90"))
  expect_false(any(built$fill[!is.na(g$data$Value) & g$data$Value == 0] == "grey90"))
  expect_identical(plot_data(g, "matrix"), mat)
  expect_match(paste(plot_data(g, "notes")$Text, collapse = " "), "not zero")
  expect_null(g$labels$title)
  expect_null(g$labels$caption)
  expect_identical(ggplot2::ggplot_build(as_ggplot(saved, component = "matrix"))$data,
    ggplot2::ggplot_build(g)$data)
  dif <- new_mfrm_plot_data("dif_heatmap", list(matrix = 2 * mat))
  other <- as_ggplot(dif)
  expect_identical(other$labels$x, "Group or contrast")
  expect_identical(other$labels$y, "Facet level")
  expect_null(other$scales$get_scales("fill")$limits)
  expect_equal(other$data$Value, as.vector(2 * mat))
})


agreement_plot_fixture <- function() {
  as_mfrm_bundle(list(
    summary = data.frame(ClassifiedPairs = 3L),
    pairs = data.frame(Rater1 = "A", Rater2 = c("B", "C", "D", "E"),
      N = c(10, 100, 0, 20), Exact = c(.5, .9, NA, .75),
      ExpectedExact = c(.45, NA, NA, .6), Corr = c(-.7, .8, NA, 0),
      MeanDiff = c(-.5, .1, NA, .1), MAD = c(.8, .2, NA, .4),
      Flag = c(TRUE, FALSE, NA, TRUE)),
    settings = list(rater_facet = "Judge", context_facets = c("Person", "Task"),
      exact_warn = .8, corr_warn = .6),
    notes = "Review thresholds do not establish rater quality."
  ), "mfrm_interrater")
}

test_that("saved agreement plots use their own review thresholds", {
  a <- agreement_plot_fixture()
  expect_equal(plot(a, draw = FALSE)$data$reference_lines$value, .8)
  expect_equal(plot(a, type = "agreement", metric = "correlation", draw = FALSE)$data$reference_lines$value, .6)
  expect_error(plot(a, exact_warn = .5, draw = FALSE), "saved thresholds")
  expect_error(plot(a, corr_warn = .3, draw = FALSE), "saved thresholds")
  expect_no_error(plot(a, exact_warn = .8, corr_warn = .6, draw = FALSE))
  expect_equal(nrow(plot(a, top_n = 1, draw = FALSE)$data$table), 1L)
  expect_equal(nrow(plot(a, type = "agreement", metric = "difference", top_n = 1, draw = FALSE)$data$table), 4L)
})

test_that("agreement ggplots preserve selected metrics, missing pairs and expected values", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  a <- agreement_plot_fixture()
  for (metric in c("exact", "correlation", "difference")) {
    p <- plot(a, type = "agreement", metric = metric, draw = FALSE,
      palette = c(ok = "navy", flag = "orange", expected = "purple"))
    g <- as_ggplot(p); b <- ggplot2::ggplot_build(g)
    tab <- p$data$table
    if (metric == "difference") {
      expect_equal(b$data[[1]]$x, tab$MeanDiff)
      expect_equal(b$data[[1]]$y, tab$MAD)
      expect_identical(g$labels$y, "Mean absolute difference")
    } else {
      expect_equal(b$data[[1]]$y, if (metric == "exact") tab$Exact else tab$Corr)
      expect_equal(as.numeric(b$data[[1]]$x), seq_len(nrow(tab)))
      expect_equal(g$scales$get_scales("y")$limits, if (metric == "exact") c(0, 1) else c(-1, 1))
      expect_identical(g$labels$x, "Rater pair")
      expect_identical(b$data[[1]]$fill, ifelse(is.na(tab$Flag),
        resolve_plot_preset("standard")$neutral, ifelse(tab$Flag, "orange", "navy")))
    }
    if (metric == "exact") {
      expect_equal(b$data[[2]]$y, tab$ExpectedExact)
      expect_equal(b$data[[3]]$y, tab$ExpectedExact)
      expect_true(all(b$data[[3]]$colour == "purple"))
    }
    ref <- b$data[[length(b$data)]]
    expect_equal(if (metric == "difference") ref$xintercept else ref$yintercept,
      if (metric == "difference") 0 else if (metric == "exact") .8 else .6)
    expect_identical(plot_data(g), p$data)
    expect_null(g$labels$title)
    expect_null(g$labels$caption)
    expect_silent(ggplot2::ggplotGrob(g))
    expect_equal(ggplot2::ggplot_build(as_ggplot(p, component = "table"))$data, b$data)
    expect_error(as_ggplot(p, component = "summary"), "plot_data")
    expect_error(as_ggplot(p, typo = TRUE), "empty")
  }
  expect_silent(ggplot2::ggplotGrob(as_ggplot(plot(a, top_n = 1, draw = FALSE))))
  a$pairs[, c("Exact", "ExpectedExact", "Corr", "MeanDiff", "MAD")] <- NA_real_
  for (metric in c("exact", "correlation", "difference"))
    expect_silent(ggplot2::ggplotGrob(as_ggplot(plot(a, type = "agreement", metric = metric, draw = FALSE))))
})

test_that("agreement heatmaps retain facet labels, metric scales and monochrome", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  a <- agreement_plot_fixture()
  local_mocked_bindings(interrater_agreement_table = function(...) a, .package = "mfrmr")
  fit <- structure(list(), class = "mfrm_fit")
  for (preset in c("standard", "monochrome")) for (metric in c("exact", "correlation")) {
    p <- plot(fit, type = "agreement", view = "heatmap", rater_facet = "Judge",
      metric = metric, preset = preset, draw = FALSE)
    g <- as_ggplot(p); b <- ggplot2::ggplot_build(g)
    expect_identical(g$labels$x, "Judge")
    expect_identical(g$labels$y, "Judge")
    expect_equal(g$scales$get_scales("fill")$limits, if (metric == "exact") c(0, 1) else c(-1, 1))
    expect_true(all(is.na(diag(p$data$matrix))))
    values <- as.numeric(p$data$matrix)
    expect_identical(g$data$Value, values)
    expect_true(all(b$data[[1]]$fill[is.na(values)] == "grey90"))
    expect_true(all(b$data[[2]]$label[is.na(values)] == ""))
    expect_identical(plot_data(g), p$data)
    expect_silent(ggplot2::ggplotGrob(g))
    if (preset == "monochrome") {
      rgb <- grDevices::col2rgb(p$data$colors)
      expect_equal(rgb[1, ], rgb[2, ])
      expect_equal(rgb[2, ], rgb[3, ])
    }
    expect_error(as_ggplot(p, component = "pairs"), "plot_data")
  }
})
