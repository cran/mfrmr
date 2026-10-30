rater_effect_bias_fixture <- function() {
  tab <- data.frame(Pair = c("r1 | c1", "r2 | c1"),
    FacetA = "Rater", FacetB = "Criterion", FacetA_Level = c("r1", "r2"), FacetB_Level = "c1",
    ObsExpAverage = c(.1, -.1), BiasSize = c(.2, -.2), SE = .1,
    t = c(2, -2), Prob = .05, AbsT = 2, Flag = TRUE,
    TFlag = TRUE, BiasFlag = FALSE, PFlag = TRUE, ObservedCount = 10,
    ProfileCILower = c(0, -5), ProfileCIUpper = c(.4, .1), ProfileCILevel = .95,
    ProfileCIStatus = c("ok", "limited by search range"))
  as_mfrm_bundle(list(ranked_table = tab, scatter_data = tab,
    summary = data.frame(InteractionFacets = "Rater x Criterion"),
    thresholds = list(abs_t_warn = 3, abs_bias_warn = .8, p_max = .01)), "mfrm_bias_interaction")
}

rater_effect_halo_fixture <- function() {
  as_mfrm_bundle(list(summary = data.frame(AvailablePairs = 2, PValueAvailablePairs = 2),
    settings = list(rater_facet = "Judge", alpha = .05, min_pair_n = 5),
    caveats = data.frame(Message = "No causal halo diagnosis; correlations are dependent."),
    node_metrics = data.frame(Node = c("r1::c1", "r1::c2", "r2::c1"),
      Rater = c("r1", "r1", "r2"), Criterion = c("c1", "c2", "c1")),
    edge_metrics = data.frame(From = "r1::c1", To = "r1::c2", Weight = .8, SignedWeight = -.8, EdgeType = "halo"),
    pair_metrics = data.frame(N = c(10, 12, 3), AbsEstimate = c(.8, .3, .5),
      RetainedByN = c(TRUE, TRUE, FALSE), EdgeType = c("halo", "non_halo", "halo")),
    halo_summary_by_rater = data.frame(Rater = c("r1", "r2"), ExpectedHaloPairs = c(8, 4),
      MeanHaloWeight = c(.8, NA), ReviewStatus = c("warning", "insufficient_data"))), "mfrm_halo_network")
}

test_that("bias conversion preserves the selected geometry and conditional bounds", {
  skip_if_not_installed("ggplot2")
  x <- rater_effect_bias_fixture()
  local_mocked_bindings(estimate_bias = function(...) stop("No bias estimation"),
    fit_mfrm = function(...) stop("No refit"))
  p <- plot(x, type = "scatter", show_ci = TRUE, draw = FALSE, palette = c(flag = "purple"))
  expect_match(p$data$subtitle, "0.8", fixed = TRUE)
  expect_match(p$data$subtitle, "0.01", fixed = TRUE)
  expect_error(plot(x, abs_bias_warn = .5, draw = FALSE), "saved thresholds")
  expect_no_error(plot(x, abs_bias_warn = .8, draw = FALSE))
  g <- as_ggplot(p); b <- ggplot2::ggplot_build(g)
  expect_equal(b$data[[1]]$y, c(0, -5))
  expect_equal(b$data[[1]]$yend, c(.4, .1))
  expect_identical(b$data[[1]]$linetype, c("solid", "dashed"))
  expect_equal(b$data[[2]]$x, c(.1, -.1))
  expect_equal(b$data[[2]]$y, c(.2, -.2))
  expect_true(all(b$data[[2]]$colour == "purple"))
  expect_equal(unlist(lapply(b$data, `[[`, "yintercept")), c(-.8, 0, .8))
  expect_identical(plot_data(g), p$data)
  notes <- paste(plot_data(g, "notes")$Text, collapse = "\n")
  expect_match(notes, "not overall rater severity, estimator bias", fixed = TRUE)
  expect_match(notes, "Positive Bias Size increases expected coded scores", fixed = TRUE)
  expect_match(notes, "discards direction", fixed = TRUE)
  expect_match(notes, "calibration uncertainty is not propagated", fixed = TRUE)
  expect_null(g$labels$title); expect_null(g$labels$caption)
  p <- plot(x, type = "ranked", show_ci = TRUE, top_n = 1, draw = FALSE)
  b <- ggplot2::ggplot_build(as_ggplot(p))
  expect_equal(b$data[[1]]$xend, .4)
  expect_equal(b$data[[2]]$x, .2)
  expect_equal(nrow(b$data[[2]]), 1)
  expect_error(as_ggplot(p, component = "scatter_data"), "selected rater-effect view")
  x$ranked_table$ProfileCILower <- x$ranked_table$ProfileCIUpper <- NA_real_
  p <- plot(x, type = "ranked", show_ci = TRUE, draw = FALSE)
  g <- as_ggplot(p)
  expect_true(all(is.na(plot_data(g)$ranked_table$CI_Lower)))
  expect_equal(nrow(ggplot2::ggplot_build(g)$data[[1]]), 0)
  for(view in c("heatmap", "facet_profile", "abs_t_hist")) {
    p <- new_mfrm_plot_data("table13_bias", list(plot = view))
    expect_error(as_ggplot(p), "unavailable")
    expect_error(as_ggplot(p, component = "plot_long"), "unavailable")
  }
})

test_that("halo plots keep available-pair metrics, missing states and cautions", {
  skip_if_not_installed("ggplot2")
  x <- rater_effect_halo_fixture()
  before <- x
  p <- plot(x, type = "edge_distribution", draw = FALSE)
  g <- as_ggplot(p)
  expect_equal(g$data$AbsEstimate, c(.8, .3))
  expect_identical(g$labels$y, "Absolute correlation")
  expect_match(paste(plot_data(g, "notes")$Text, collapse = " "), "No causal halo")
  expect_identical(p$data$settings, x$settings)
  p <- plot(x, type = "halo_summary", metric = "MeanHaloWeight", draw = FALSE,
    palette = c(warning = "purple", neutral = "grey50"))
  g <- as_ggplot(p);b <- ggplot2::ggplot_build(g)
  expect_equal(b$data[[1]]$y, .8)
  expect_identical(b$data[[1]]$fill, "purple")
  expect_equal(b$data[[2]]$label, "Unavailable")
  expect_identical(g$labels$x, "Judge")
  expect_true("Insufficient information" %in% p$data$legend$label)
  expect_false("OK" %in% p$data$legend$label)
  expect_identical(plot_data(g), p$data)
  p <- plot(x, type = "halo_summary", metric = "ExpectedHaloPairs", draw = FALSE)
  expect_equal(ggplot2::ggplot_build(as_ggplot(p))$data[[1]]$y, c(8, 4))
  p <- plot(x, type = "matrix", draw = FALSE);g <- as_ggplot(p)
  expect_equal(g$scales$get_scales("fill")$limits, c(-1, 1))
  expect_identical(g$labels$x, "Rater x criterion node")
  expect_equal(sort(g$data$Value[is.finite(g$data$Value)]), c(-.8, -.8))
  expect_true(anyNA(g$data$Value))
  expect_identical(g$scales$get_scales("fill")$na.value, "grey90")
  expect_null(g$labels$title); expect_null(g$labels$caption)
  missing <- x
  missing$halo_summary_by_rater$MeanHaloWeight[] <- NA_real_
  grDevices::pdf(NULL, width = 5, height = 3); on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(missing, type = "halo_summary"))
  expect_no_error(ggplot2::ggplot_build(as_ggplot(plot(missing, type = "halo_summary", draw = FALSE))))
  expect_identical(x, before)
  expect_error(as_ggplot(new_mfrm_plot_data("halo_network", list(plot = "network"))), "unavailable")
})

test_that("halo network draws without automatic titles and retains its notes", {
  skip_if_not_installed("igraph")
  x <- rater_effect_halo_fixture()
  seen <- NULL
  local_mocked_bindings(.mfrm_plot_igraph = function(..., main = NULL) { seen <<- main }, .package = "mfrmr")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(x, type = "network")
  expect_identical(seen, "")
  expect_match(paste(plot_data(p, "notes")$Text, collapse = " "), "not rater centrality")
  expect_identical(p$data$caveats, x$caveats)
  filtered <- plot(x, type = "network", display_min_weight = 1, draw = FALSE)
  expect_identical(filtered$data$caveats, x$caveats)
  expect_identical(filtered$data$settings, x$settings)
  expect_identical(filtered$data$source_edges, x$edge_metrics)
  expect_equal(nrow(filtered$data$edges), 0L)
  expect_match(paste(plot_data(filtered, "notes")$Text, collapse = " "),
    "No causal halo diagnosis; correlations are dependent.", fixed = TRUE)
  expect_match(paste(plot_data(filtered, "notes")$Text, collapse = " "),
    "Missing or filtered edges do not mean zero correlation or absence of halo.", fixed = TRUE)
})

severity_plot_fixture <- function() {
  d <- data.frame(Facet = 'Rater', Level = c('regular', 'fixed', 'point only', 'missing', 'invalid'),
    Estimate = c(-.4, 0, 4, NA, .8), ModelSE = c(.2, 0, NA, .1, -.1), SE = 99,
    Fixed = c(FALSE, TRUE, FALSE, FALSE, FALSE), SE_Method = 'Observation-table information',
    PrecisionTier = 'exploratory', SupportsFormalInference = FALSE,
    SEUse = 'screening_only', CIBasis = 'Exploratory normal band', CIUse = 'screening_only',
    CIEligible = FALSE, CILabel = 'Screening only')
  list(fit = structure(list(config = list(facet_signs = c(Rater = -1))), class = 'mfrm_fit'),
    diagnostics = structure(list(measures = d), class = 'mfrm_diagnostics'))
}

mock_severity_source <- function(env = parent.frame()) {
  local_mocked_bindings(
    .mfrm_fit_plot_readiness = function(...) list(ready = FALSE, status = 'review_only',
      detail = 'Review source fit.', table = data.frame(Domain = 'Fit', Status = 'review')),
    mfrm_results_validate_diagnostics_identity = function(...) 'geometry_fixture',
    diagnose_mfrm = function(...) stop('No new diagnostics'),
    fit_mfrm = function(...) stop('No fit'), .env = env)
}

test_that('severity profiles retain SE basis, missing rows and fixed points', {
  mock_severity_source()
  x <- severity_plot_fixture()
  p <- plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics, ci_level = .9,
    draw = FALSE)
  d <- p$data$data
  expect_false(p$data$show_bands)
  expect_false(any(p$data$legend$role == 'band'))
  expect_identical(d$Band, c('|location| <= 0.5', '|location| <= 0.5',
    '0.5 < |location| <= 1.0', '|location| > 1.0', NA_character_))
  expect_equal(d$Level, c('regular', 'fixed', 'invalid', 'point only', 'missing'))
  expect_equal(d$SE, c(.2, 0, -.1, NA, .1))
  expect_equal(d$CI_Lower[1], -.4 - qnorm(.95)*.2)
  expect_equal(d$CI_Upper[1], -.4 + qnorm(.95)*.2)
  expect_true(all(is.na(d$CI_Lower[-1])))
  expect_equal(d$IntervalStatus, c('normal_approximation', 'fixed', 'se_nonpositive', 'se_unavailable', 'estimate_unavailable'))
  expect_identical(d$CIUse, c('screening_only', 'not_applicable', rep('unavailable', 3)))
  expect_false(any(d$CIEligible | d$SupportsFormalInference))
  expect_true(all(d$SEColumn == 'ModelSE'))
  expect_identical(p$data$interpretation_status, 'review_only')
  expect_identical(p$data$source_status, 'geometry_fixture')
  expect_identical(p$data$axis_label, 'Rater location (logit)')
  expect_match(paste(p$data$notes$Text, collapse = ' '), 'not likelihood-profile')
  expect_match(paste(p$data$notes$Text, collapse = ' '), 'not homogeneity, practical equivalence')
  expect_equal(.mfrm_severity_geometry(d)$xrange[2], 4)
  x$diagnostics$measures$ModelSE <- NULL
  expect_true(all(plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics,
    draw = FALSE)$data$data$SEColumn == 'SE'))
  x$diagnostics$measures$SE_Method <- NULL
  unknown <- plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics, draw = FALSE)$data$data
  expect_true(all(unknown$PrecisionTier == 'unknown'))
  expect_false(any(unknown$CIEligible | unknown$SupportsFormalInference))
  x$diagnostics$measures$SE <- NULL
  expect_true(all(is.na(plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics,
    draw = FALSE)$data$data$CI_Lower)))
  x$fit$config$facet_signs[] <- 1
  expect_identical(plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics,
    draw = FALSE)$data$axis_label, 'Rater location (logit)')
  expect_error(plot(x$fit, type = 'rater_severity', facet = 'Person', draw = FALSE), 'non-Person')
  expect_error(plot(x$fit, type = 'rater_severity', show_bands = NA, draw = FALSE), 'show_bands')
})

test_that('severity ggplot uses stored endpoints and complete level rows', {
  skip_if_not_installed('ggplot2')
  mock_severity_source()
  x <- severity_plot_fixture()
  p <- plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics, show_bands = FALSE, draw = FALSE)
  # Conversion must use the stored bounds, not recompute from SE.
  p$data$data$CI_Lower[1] <- -2; p$data$data$CI_Upper[1] <- .25
  g <- as_ggplot(p); b <- ggplot2::ggplot_build(g)$data
  expect_equal(b[[2]]$x, -2); expect_equal(b[[2]]$xend, .25)
  expect_equal(b[[3]]$x, c(-.4, 0, .8, 4))
  expect_equal(b[[3]]$shape, c(19, 5, 4, 4))
  expect_equal(b[[4]]$label, 'Unavailable')
  expect_equal(g$scales$get_scales('y')$labels, p$data$data$Level)
  expect_identical(plot_data(g), p$data)
  expect_null(g$labels$title); expect_null(g$labels$caption)
  expect_equal(ggplot2::ggplot_build(as_ggplot(p, component = 'data'))$data, b)
  expect_error(as_ggplot(p, component = 'notes'), 'severity profile')
  expect_error(as_ggplot(p, ci_level = .8), 'must be empty')
  p$data$show_bands <- TRUE
  bands <- ggplot2::ggplot_build(as_ggplot(p))$data[[1]]
  expect_equal(bands$xmin, c(-.5, -1, .5)); expect_equal(bands$xmax, c(.5, -.5, 1))
  expect_equal(bands$alpha, c(.1, .08, .08))
  expect_identical(bands$fill, p$data$band_colors)
  expect_length(unique(bands$fill), 1L)
  # Older saved payloads retain geometry without acquiring eligibility fields.
  p$data$data <- p$data$data[,c('Level', 'Estimate', 'SE', 'CI_Lower', 'CI_Upper', 'Band')]
  p$data$band_colors <- NULL
  p$data$data$Band <- c('gentle', 'gentle', 'moderate', 'strict', NA_character_)
  legacy <- serialize(p, NULL)
  legacy_bands <- ggplot2::ggplot_build(as_ggplot(p))$data[[1]]
  style <- resolve_plot_preset(p$data$preset)
  expect_identical(legacy_bands$fill, c(style$success, style$warn, style$warn))
  expect_identical(plot_data(as_ggplot(p)), p$data)
  expect_identical(serialize(p, NULL), legacy)
})

test_that('severity guides are opt-in and use the same neutral colors in both renderers', {
  mock_severity_source()
  x <- severity_plot_fixture()
  # Equal distance from zero does not imply a small pairwise difference.
  x$diagnostics$measures$Estimate[1:2] <- c(-.5, .5)
  rectangles <- character()
  local_mocked_bindings(rect = function(..., col = NULL) {
    rectangles <<- c(rectangles, col)
  }, .package = "mfrmr")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics)
  expect_length(rectangles, 0L)
  p <- plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics, show_bands = TRUE)
  expect_identical(p$data$data$Band[1:2], rep('|location| <= 0.5', 2L))
  expect_equal(diff(p$data$data$Estimate[1:2]), 1)
  expected <- vapply(seq_len(3L), function(i)
    grDevices::adjustcolor(p$data$band_colors[i], alpha.f = c(.10, .08, .08)[i]), '')
  expect_identical(rectangles, unname(expected))
  expect_identical(p$data$legend$value[p$data$legend$role == 'band'], p$data$band_colors[1:2])
})

test_that('severity base graphics retain unavailable rows and restore margins', {
  mock_severity_source()
  x <- severity_plot_fixture()
  grDevices::pdf(NULL, width = 5, height = 3); on.exit(grDevices::dev.off(), add = TRUE)
  old_mar <- graphics::par('mar')
  expect_no_error(plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics))
  expect_identical(graphics::par('mar'), old_mar)
  x$diagnostics$measures$Estimate[] <- NA_real_
  expect_no_error(plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics))
  if (requireNamespace('ggplot2', quietly = TRUE)) {
    p <- plot(x$fit, type = 'rater_severity', diagnostics = x$diagnostics, draw = FALSE)
    expect_no_error(ggplot2::ggplot_build(as_ggplot(p)))
  }
})

test_that('severity source and model restrictions cannot be bypassed by diagnostics', {
  x <- list(fit = make_toy_fit())
  x$diagnostics <- make_toy_diagnostics(x$fit)
  local_mocked_bindings(.mfrm_fit_plot_readiness = function(...) list(),
    diagnose_mfrm = function(...) stop('Unexpected diagnostics'))
  expect_identical(plot_rater_severity_profile(x$fit, x$diagnostics, draw = FALSE)$data$source_status,
    'structural_identity_match')
  x$diagnostics$measures$Estimate[1] <- 10
  expect_error(plot_rater_severity_profile(x$fit, x$diagnostics, draw = FALSE), 'fit/diagnostics mismatch')
  x$fit$config$jml_adjustment <- TRUE
  expect_error(plot_rater_severity_profile(x$fit, x$diagnostics, draw = FALSE), 'corrected JML')
  x$fit$config$jml_adjustment <- NULL
  x$fit$config$gpcm_spec <- list(component_levels = list(Task = c('a','b'), Rater = c('x','y')))
  expect_error(plot_rater_severity_profile(x$fit, x$diagnostics, draw = FALSE), 'two slope families')
})
