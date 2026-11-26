test_that("workflow comparison delegates nested likelihood checks", {
  d <- mfrmr::load_mfrmr_data("example_core")
  spec <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score")
  rsm <- step_mfrm_measure(spec)
  pcm <- step_mfrm_measure(spec, model = "PCM", step_facet = "Criterion")
  gpcm <- step_mfrm_measure(spec, model = "GPCM", step_facet = "Criterion", slope_facet = "Rater",
                           gpcm_mml_identification = "fixed_standard_normal")
  cmp <- anova(rsm, pcm, gpcm)
  expect_identical(cmp$table$Family, c("RSM", "PCM", "GPCM"))
  expect_equal(cmp$lrt$ChiSq, 2 * diff(cmp$table$LogLik), tolerance = 1e-12)
  expect_equal(cmp$lrt$df, diff(cmp$table$FreeParams))
  unmatched <- step_mfrm_measure(spec, model = "PCM", step_facet = "Criterion", quad_points = 41L)
  expect_warning(ineligible <- anova(rsm, unmatched), "Nested LRT requires")
  expect_null(ineligible$lrt)
})

test_that("context curves reflect fitted discrimination for one and two families", {
  d <- mfrmr::load_mfrmr_data("example_core")
  spec <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score")
  contexts <- data.frame(Rater = c("R01", "R02"), Criterion = "Content")
  for (sf in list("Rater", c("Rater", "Criterion"))) {
    formula <- if (length(sf) == 1L) "full_predictor" else "ability_only"
    wf <- step_mfrm_measure(spec, model = "GPCM", step_facet = "Criterion", slope_facet = sf, slope_type = formula)
    p <- plot(wf, "category_curves", draw = FALSE, newdata = contexts)
    expect_equal(length(unique(p$data$RatingContext)), 2L)
    f <- wf$fit
    for (i in 1:2) {
      # Evaluate the ordered-category logits independently at the same context.
      th <- 0.4
      sub <- p$data[p$data$ContextID == i & abs(p$data$theta - th) < 1e-10, ]
      locations <- f$facets[f$facets$Level %in% as.character(contexts[i, ]), ]
      b <- sum(locations$Estimate * locations$Sign)
      slopes <- f$slopes[f$slopes$Level %in% as.character(contexts[i, ]), ]
      a <- prod(slopes$Estimate)
      cum <- f$step_cumulative[match("Content", f$encoding$step_owners), ]
      lp <- if (formula == "full_predictor") a * (0:3 * (th + b) - cum) else 0:3 * (a * th + b) - cum
      expected <- exp(lp) / sum(exp(lp))
      expect_equal(sub$Probability, expected, tolerance = 1e-12)
    }
  }
})

test_that("the frontier uses the eligible grid and the two targets separately", {
  d <- mfrmr::load_mfrmr_data("example_core")
  wf <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score") |>
    step_mfrm_gtheory(d_study_grid = expand.grid(Rater = 1:5, Criterion = 1:5), target_g = .80, target_phi = .85)
  p <- plot(wf, "dstudy_frontier", draw = FALSE)
  expect_equal(nrow(p$data), 2L * nrow(wf$gtheory$d_study))
  step <- p$layers[[which(vapply(p$layers, function(x) inherits(x$geom, "GeomStep"), logical(1)))]]$data
  for (metric in unique(step$Metric)) {
    draws <- p$data[p$data$Metric == metric, ]
    boundary <- step[step$Metric == metric, ]
    expected <- vapply(boundary$Ratings, function(n) max(draws$Coefficient[draws$Ratings <= n]), numeric(1))
    expect_equal(boundary$Coefficient, expected)
    expect_true(all(diff(boundary$Coefficient) > 0))
  }
  labels <- p$layers[[which(vapply(p$layers, function(x) inherits(x$geom, "GeomLabel"), logical(1)))]]$data
  expect_equal(setNames(labels$Ratings, labels$Metric)[c("G (relative)", "Phi (absolute)")],
               c("G (relative)" = 8, "Phi (absolute)" = 15))
  standard <- ggplot2::ggplot_build(plot(wf, "dstudy", draw = FALSE))$data
  lines <- unlist(lapply(standard, `[[`, "yintercept"))
  expect_equal(sort(lines), c(.80, .85))

  d$Score <- 2
  undefined <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score", rating_min = 1, rating_max = 4) |>
    step_mfrm_gtheory()
  q <- plot(undefined, "dstudy_frontier", draw = FALSE)
  expect_true(all(is.na(q$data$Coefficient)))
  expect_true(any(grepl("Zero total variance", attr(q, "mfrm_notes"))))
  expect_silent(ggplot2::ggplot_build(q))
})

test_that("subscore and composite points use the observed assignment and weights", {
  d <- mfrmr::load_mfrmr_data("example_core")
  d <- d[d$Criterion %in% c("Content", "Organization"), c("Person", "Rater", "Criterion", "Score")]
  d <- reshape(d, idvar = c("Person", "Rater"), timevar = "Criterion", direction = "wide")
  wf <- mfrm_workflow(d, "Person", "Rater", "Score.Content") |>
    step_mfrm_gtheory(g_scores = c("Score.Content", "Score.Organization"), g_weights = c(.6, .4))
  gt <- mfrmr::mfrm_g_study(d, person = "Person", facets = "Rater", scores = c("Score.Content", "Score.Organization"),
                          weights = c(.6, .4), method = "urgenova")
  p <- plot(wf, "gstudy", draw = FALSE)
  expect_equal(p$data$Coefficient, c(gt$g_study_summary$G, gt$univariate$G, gt$g_study_summary$Phi, gt$univariate$Phi))
  expect_true(any(grepl("Score.Content 0.6", attr(p, "mfrm_notes"), fixed = TRUE)))
  q <- plot(wf, "dstudy_frontier", draw = FALSE)
  ds <- mfrmr::mfrm_d_study(gt)
  expect_equal(q$data$Coefficient, c(ds$G, ds$Phi))
})
