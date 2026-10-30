mvgt_plot_data <- function() {
  path <- testthat::test_path("..", "..", "inst", "extdata", "mgenova-table12.csv")
  if (!file.exists(path)) path <- system.file("extdata", "mgenova-table12.csv", package = "mfrmr")
  read.csv(path)
}

test_that("error components retain covariance and sum to planned error variances", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(1, 6)),
    weights = cbind(Equal = c(V = .5, W = .5), Difference = c(V = -1, W = 1)))
  before <- serialize(d, NULL)
  local_mocked_bindings(mfrm_g_study = function(...) stop("refitted"),
    confint.mfrm_multivariate_d_study = function(...) stop("intervals recomputed"))
  p <- plot(d, type = "error_components", composite = c("Difference", "Equal"), draw = FALSE)
  z <- plot_data(p, "table")
  for (s in c("Difference", "Equal")) for (i in 1:2) for (metric in unique(z$Metric)) {
    a <- z[z$Score == s & z$Scenario == i & z$Metric == metric, ]
    w <- d$weights[, s]
    for (k in seq_len(nrow(a))) {
      m <- g$components[[a$Source[k]]]
      expect_equal(a$Contribution[k], sum(outer(w, w)*m)/d$design_grid$Tasks[i])
      expect_equal(a$CovarianceContribution[k], unname(2*w[1]*w[2]*m[1, 2]/d$design_grid$Tasks[i]))
    }
    expect_equal(sum(a$Contribution), d$coefficients[
      d$coefficients$Score == s & d$coefficients$Scenario == i, metric])
    expect_equal(a$DiagonalContribution+a$CovarianceContribution, a$Contribution)
  }
  expect_false(any(z$Source == "Person"))
  expect_equal(unique(z$Source[z$Metric == "RelativeErrorVariance"]), "Residual")
  expect_identical(serialize(d, NULL), before)
  expect_identical(plot_data(p)$weights, d$weights[, c("Difference", "Equal")])
  expect_true(any(grepl("not SEM", plot_data(p)$notes$Text)))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(p, path)
  expect_identical(plot_data(readRDS(path)), plot_data(p))
})

test_that("error contributions use crossed and nested planned counts in facet order", {
  data <- mfrmr_example_gtheory_crossed
  for (nested in c(FALSE, TRUE)) for (reverse in c(FALSE, TRUE)) {
    facets <- if (reverse) c(Task = "Task", Rater = "Rater") else c(Rater = "Rater", Task = "Task")
    g <- mfrm_g_study(data, c("Content", "Organization"), facets = facets,
      nesting = if (nested) c(Rater = "Task") else NULL)
    # Inject known signed components after fitting: the plot tests the D-study
    # projection, independently of variance-component estimation.
    g$components <- lapply(seq_along(g$components), function(i)
      matrix(c(i, -.2*i, -.2*i, 2*i), 2,
        dimnames = list(c("Content", "Organization"), c("Content", "Organization")))) |>
      stats::setNames(names(g$components))
    d <- mfrm_d_study(g, data.frame(Rater = c(2, 50000L), Task = c(3, 50000L)),
      weights = c(Content = -1, Organization = 2))
    z <- plot_data(plot(d, type = "error_components", draw = FALSE), "table")
    for (i in 1:2) for (metric in unique(z$Metric)) {
      a <- z[z$Scenario == i & z$Metric == metric, ]
      nr <- d$design_grid$Rater[i]; nt <- d$design_grid$Task[i]
      expected <- if (nested) c(Task = 1/nt, "Rater(Task)" = 1/(nr*nt),
        "Person:Task" = 1/nt, Residual = 1/(nr*nt)) else
        c(Rater = 1/nr, Task = 1/nt, "Person:Rater" = 1/nr,
          "Person:Task" = 1/nt, "Rater:Task" = 1/(nr*nt),
          "Task:Rater" = 1/(nr*nt), Residual = 1/(nr*nt))
      expect_equal(a$Multiplier, unname(expected[a$Source]))
      expect_equal(a$Contribution, a$Multiplier*vapply(a$Source, function(s)
        sum(outer(c(-1, 2), c(-1, 2))*g$components[[s]]), numeric(1)),
        ignore_attr = TRUE)
      expect_equal(sum(a$Contribution), a$TotalErrorVariance[1])
    }
    if (nested) expect_true(all(grepl("Rater per Task", z$Plan, fixed = TRUE)))
  }
})

test_that("signed error components preserve unavailable totals and reject altered identities", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  g$components$Residual <- -abs(g$components$Residual)
  d <- mfrm_d_study(g, data.frame(Tasks = 1:2))
  p <- plot(d, type = "error_components", score = "V", draw = FALSE)
  expect_true(any(plot_data(p, "table")$Contribution < 0))
  expect_true(all(plot_data(p)$totals$RelativeSEMStatus != "Available"))
  expect_true(all(is.na(plot_data(p)$totals$RelativeSEM)))
  draw <- function(x = d, ...) plot(x, type = "error_components", score = "V", draw = FALSE, ...)
  for (alter in c("component", "weights", "grid", "coefficient")) {
    bad <- d
    if (alter == "component") bad$gstudy$components$Residual <- bad$gstudy$components$Residual*2
    if (alter == "weights") bad$weights <- c(V = .5, W = .5)
    if (alter == "grid") bad$design_grid$Tasks <- c(3, 4)
    if (alter == "coefficient") bad$coefficients$RelativeErrorVariance[1] <- 100
    expect_error(draw(bad), "do not match")
  }
  expect_error(draw(metric = "G"), "metric")
  expect_error(draw(ncol = 0), "ncol")
  expect_error(draw(x_var = "Tasks"), "complete plans")
  expect_error(draw(reference = .8), "complete plans")
  expect_error(draw(intervals = TRUE), "complete plans")
  expect_error(draw(method = "delta"), "intervals")
})

test_that("error component renderers preserve values, style, labels and graphics settings", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(1, 6)))
  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  before <- graphics::par()[c("mfrow", "mar", "oma", "cex", "mex", "fig", "mfg")]
  p <- expect_no_warning(plot(d, type = "error_components", preset = "monochrome"))
  expect_identical(graphics::par()[names(before)], before)
  expect_equal(length(unique(plot_data(p)$legend$value)), 1L)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p)
    expect_no_warning(ggplot2::ggplotGrob(gg))
    expect_identical(gg$data$Contribution, plot_data(p, "table")$Contribution)
    expect_identical(gg$labels$y, plot_data(p)$x_label)
    expect_null(gg$labels$title)
    expect_identical(plot_data(gg, "totals"), plot_data(p, "totals"))
  }
  one <- mfrm_d_study(g, data.frame(Tasks = 1))
  expect_no_warning(plot(one, type = "error_components", metric = "RelativeErrorVariance"))
})

test_that("D-study ribbons reuse matching saved intervals without recomputation", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(12, 1, 6)))
  ci <- lapply(c("V", "W"), function(s)
    confint(d, score = s, parm = "Phi", assumption = "normal", draws = 4000))
  before <- serialize(d, NULL)
  # The plot must not invoke either inference or model estimation.
  local_mocked_bindings(confint.mfrm_multivariate_d_study = function(...) stop("recomputed"),
    mfrm_g_study = function(...) stop("refitted"))
  p <- plot(d, score = c("V", "W"), metric = "Phi", intervals = rev(ci), draw = FALSE)
  z <- plot_data(p)
  for (interval in ci) {
    rows <- z$series[z$series$Score == interval$Score[1], ]
    expect_identical(rows$Lower, interval$Lower)
    expect_identical(rows$Upper, interval$Upper)
    expect_identical(rows$Value, interval$Estimate)
  }
  expect_equal(length(unique(z$series$IntervalRun)), 2)
  expect_equal(nrow(z$unavailable_intervals), 0)
  expect_identical(z$intervals, rev(ci))
  expect_identical(serialize(d, NULL), before)
  expect_true(any(grepl("simultaneous coverage", z$notes$Text)))
  expect_true(any(grepl("frequentist coverage", z$notes$Text)))
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, score = c("V", "W"), metric = "Phi", intervals = ci))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p)
    expect_identical(gg$data$Lower, z$series$Lower)
    expect_identical(plot_data(gg, "intervals"), rev(ci))
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})

test_that("curve intervals reject wrong sources, targets, plans and weights", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = 1:4))
  ci <- confint(d, assumption = "normal")
  draw <- function(x = d, intervals = ci, ...) plot(x, metric = "G", intervals = intervals, draw = FALSE, ...)
  expect_error(plot(d, intervals = ci, draw = FALSE), "select one coefficient")
  expect_error(plot(d, type = "sem", metric = "AbsoluteSEM", intervals = ci), "plans, metric and score")
  expect_error(draw(intervals = as.data.frame(ci)), "confint")
  expect_error(draw(intervals = ci[-2, ]), "one interval row per plan")
  expect_error(draw(intervals = list(ci, ci)), "exactly one")
  expect_error(draw(score = c("V", "W")), "every selected")
  expect_error(draw(score = "W"), "plans, metric and score")
  other <- d; other$design_grid$Tasks <- other$design_grid$Tasks + 1
  expect_error(draw(x = other), "plans, metric and score")
  # Identical point ratios do not make a different source acceptable.
  other <- d; other$gstudy$mean_products <- lapply(g$mean_products, function(m) m*4)
  expect_error(draw(x = other), "source summaries")
  old <- ci; attr(old, "source_summary") <- NULL
  expect_error(draw(intervals = old), "source summaries")
  expect_no_error(plot(old, draw = FALSE))
  difference <- confint(d, assumption = "normal", reference = 1)
  expect_error(draw(intervals = difference), "plans, metric and score")
  other <- d; other$coefficients$G[1] <- other$coefficients$G[1] + .01
  expect_error(draw(x = other), "point estimates")
  wide <- confint(d, score = "W", level = .99, assumption = "normal")
  expect_error(draw(intervals = list(ci, wide), score = c("V", "W")), "same method and confidence")
  weighted <- mfrm_d_study(g, d$design_grid, c(V = .5, W = .5))
  wc <- confint(weighted, assumption = "normal")
  weighted$weights <- c(V = 1, W = -1)
  expect_error(draw(x = weighted, intervals = wc), "weights")
})

test_that("the ordinary ribbon workflow computes selected scores in one plot call", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = 1:4))
  expected <- lapply(c("W", "V"), function(s) confint(d, score = s, assumption = "normal", level = .9))
  # The convenience route computes intervals but never refits the G-study.
  local_mocked_bindings(mfrm_g_study = function(...) stop("refitted"))
  p <- plot(d, score = c("W", "V"), metric = "G", intervals = TRUE,
    assumption = "normal", level = .9, draw = FALSE)
  expect_identical(plot_data(p, "intervals"), expected)
  expect_identical(plot_data(p), plot_data(plot(d, score = c("W", "V"), metric = "G",
    intervals = expected, draw = FALSE)))
  expect_error(plot(d, metric = "G", intervals = TRUE), "assumption")
  expect_error(plot(d, metric = "G", intervals = TRUE, assumption = "normal",
    method = c("exact_f", "gci")), "one interval method")
  expect_error(plot(d, metric = "G", assumption = "normal"), "only when")
  expect_identical(plot_data(plot(d, draw = FALSE)), plot_data(plot(d, intervals = NULL, draw = FALSE)))
})

test_that("missing ribbons break at gaps and isolated intervals remain visible", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(5, 1, 4, 2, 3)))
  ci <- confint(d, assumption = "normal")
  # A saved failed interval in the middle must not be joined across.
  gap <- ci$Tasks == 3
  ci$Lower[gap] <- ci$Upper[gap] <- NA_real_
  ci$Available[gap] <- FALSE; ci$Status[gap] <- "Numerical interval failure"
  p <- plot(d, metric = "G", intervals = ci, draw = FALSE)
  z <- plot_data(p)
  expect_equal(z$unavailable_intervals$X, 3)
  expect_equal(z$unavailable_intervals$IntervalStatus, "Numerical interval failure")
  expect_equal(length(unique(na.omit(z$series$IntervalRun))), 2)
  expect_equal(z$series$Value, ci$Estimate)
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, metric = "G", intervals = ci, preset = "monochrome"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p); built <- ggplot2::ggplot_build(gg)
    runs <- split(built$data[[1]]$x, built$data[[1]]$group)
    expect_equal(length(runs), 2)
    expect_true(all(vapply(runs, function(x) max(x) < 3 || min(x) > 3, logical(1))))
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
  # One interval amid unavailable intervals needs a visible segment, not a zero-area ribbon.
  ci$Lower[-1] <- ci$Upper[-1] <- NA_real_; ci$Available[-1] <- FALSE
  p <- plot(d, metric = "G", intervals = ci)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p)
    expect_equal(nrow(ggplot2::ggplot_build(gg)$data[[2]]), 1)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
  ci$Lower[] <- ci$Upper[] <- NA_real_; ci$Available[] <- FALSE
  expect_no_warning(p <- plot(d, metric = "G", intervals = ci))
  expect_equal(nrow(plot_data(p, "unavailable_intervals")), 5)
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_no_warning(ggplot2::ggplotGrob(as_ggplot(p)))
})

test_that("ribbons retain crossed panels, composite identity and method limits", {
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed, c("Content", "Organization"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 4), Tasks = c(3, 6, 12)),
    cbind(Equal = c(Content = .5, Organization = .5), Difference = c(Content = 1, Organization = -1)))
  ci <- confint(d, composite = "Equal", parm = "Phi", assumption = "normal", method = "approximate_f")
  p <- plot(d, composite = "Equal", metric = "Phi", intervals = ci, panel_by = "Raters", draw = FALSE)
  z <- plot_data(p)
  expect_equal(length(unique(na.omit(z$series$IntervalRun))), 2)
  expect_identical(z$intervals[[1]]$Validation, ci$Validation)
  expect_identical(attr(z$intervals[[1]], "diagnostics"), attr(ci, "diagnostics"))
  mls <- plot(d, composite = "Equal", metric = "Phi", intervals = TRUE,
    assumption = "normal", draw = FALSE)
  expect_true(all(plot_data(mls)$intervals[[1]]$MethodCode == "mls"))
  auto <- plot(d, composite = "Equal", metric = "Phi", intervals = TRUE,
    assumption = "normal", method = "approximate_f", panel_by = "Raters", draw = FALSE)
  expect_identical(plot_data(auto), z)
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, composite = "Equal", metric = "Phi", intervals = ci, panel_by = "Raters"))
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_no_warning(ggplot2::ggplotGrob(as_ggplot(p)))
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    both <- confint(d, composite = "Equal", parm = "Phi", assumption = "normal",
      method = c("approximate_f", "mean_square_region"))
    expect_error(plot(d, composite = "Equal", metric = "Phi", intervals = both), "one interval method")
  }
})

test_that("an available interval remains visible without an admissible point", {
  n <- 10; k <- 6
  p <- qr.Q(qr(contr.helmert(n))); t <- qr.Q(qr(contr.helmert(k)))
  y <- outer(p[, 1]*sqrt(.5*(n-1)/k), rep(1, k)) + sqrt((n-1)*(k-1))*outer(p[, 2], t[, 1])
  ratings <- expand.grid(Person = 1:n, Task = 1:k); ratings$Score <- as.vector(y)
  d <- mfrm_d_study(mfrm_g_study(ratings, "Score", rater = NULL),
    data.frame(Tasks = c(6, 60)))
  ci <- confint(d, assumption = "normal")
  expect_true(all(is.na(ci$Estimate)))
  expect_true(all(ci$Available))
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(p <- plot(d, metric = "G", intervals = ci))
  expect_identical(plot_data(p, "series")$Upper, ci$Upper)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p)
    expect_equal(nrow(ggplot2::ggplot_build(gg)$data[[1]]), 2)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})

test_that("D-study plots preserve stored values, units, weights and score identities", {
  tasks <- mvgt_plot_data()
  g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(12, 3, 6)), c(V = -1, W = 1))
  before <- serialize(d, NULL)
  device <- grDevices::dev.cur()
  p <- plot_data(plot(d, draw = FALSE))
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(p$kind, "Composite")
  expect_equal(p$table, d$coefficients[d$coefficients$Kind == "Composite", ])
  expect_equal(p$series$Value[p$series$Metric == "G"], p$table$G)
  expect_equal(p$series$Value[p$series$Metric == "Phi"], p$table$Phi)
  expect_identical(p$weights, c(V = -1, W = 1))
  expect_identical(p$x_var, "Tasks")
  expect_length(p$group_var, 0)
  expect_equal(nrow(p$reference_lines), 0)
  sem <- plot_data(plot(d, type = "sem", draw = FALSE))
  expect_equal(sem$series$Value[sem$series$Metric == "RelativeSEM"], p$table$RelativeSEM)
  expect_equal(sem$series$Value[sem$series$Metric == "AbsoluteSEM"], p$table$AbsoluteSEM)
  v <- plot_data(plot(d, score = "V", draw = FALSE))
  expect_null(v$weights)
  expect_equal(v$table, d$coefficients[d$coefficients$Kind == "Score" & d$coefficients$Score == "V", ])
  expect_identical(serialize(d, NULL), before)
  # Plotting only needs stored coefficients and identities, not a refittable G-study.
  d$gstudy$components <- NULL
  expect_identical(plot_data(plot(d, draw = FALSE)), p)

  names(tasks)[names(tasks) == "V"] <- "Composite"
  g <- mfrm_g_study(tasks, c("Composite", "W"), rater = NULL)
  d <- mfrm_d_study(g, weights = c(Composite = -1, W = 1))
  expect_identical(plot_data(plot(d, draw = FALSE))$kind, "Composite")
  expect_identical(plot_data(plot(d, score = "Composite", draw = FALSE))$kind, "Score")
  d <- mfrm_d_study(g)
  expect_identical(plot_data(plot(d, draw = FALSE))$score, "Composite")
  expect_identical(plot_data(plot(d, draw = FALSE))$kind, "Score")
})

test_that("named composite plots require an unambiguous selection and keep its weights", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  # A composite may share an original score's name; Kind distinguishes them.
  w <- cbind(V = c(V = .5, W = .5), Difference = c(V = -1, W = 1))
  d <- mfrm_d_study(g, data.frame(Tasks = c(3, 6, 12)), w)
  expect_error(plot(d, draw = FALSE), "Several composites.*select one")
  expect_error(plot(d, score = "V", composite = "V", draw = FALSE), "not both")
  for (bad in list("absent", c("V", "V"), NA_character_, 1, character())) {
    expect_error(plot(d, composite = bad, draw = FALSE), "name one composite")
  }
  original <- plot_data(plot(d, score = "V", draw = FALSE))
  expect_identical(original$kind, "Score")
  expect_null(original$weights)
  for (name in colnames(w)) {
    p <- plot_data(plot(d, composite = name, type = "sem", draw = FALSE))
    expect_identical(p$kind, "Composite")
    expect_identical(p$score, name)
    expect_identical(p$weights, w[, name])
    expect_identical(p$table, d$coefficients[d$coefficients$Kind == "Composite" & d$coefficients$Score == name, ])
    expect_equal(p$series$Value[p$series$Metric == "AbsoluteSEM"], p$table$AbsoluteSEM)
    expect_match(p$title, paste0("Composite ", name, ":"), fixed = TRUE)
    expect_match(p$title, "V =", fixed = TRUE)
    expect_match(p$title, "W =", fixed = TRUE)
  }
  single <- mfrm_d_study(g, weights = w[, "Difference", drop = FALSE])
  expect_identical(plot_data(plot(single, draw = FALSE))$score, "Difference")
  expect_identical(plot_data(plot(single, draw = FALSE))$weights, w[, "Difference"])
  vector <- mfrm_d_study(g, weights = w[, "Difference"])
  expect_identical(plot_data(plot(vector, composite = "Composite", draw = FALSE)),
    plot_data(plot(vector, draw = FALSE)))
  one_score <- mfrm_d_study(
    mfrm_g_study(mvgt_plot_data(), "V", rater = NULL),
    weights = matrix(2, 1, 1, dimnames = list("V", "Double")))
  expect_identical(plot_data(plot(one_score, draw = FALSE))$weights, c(V = 2))

  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, composite = "Difference", type = "sem"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_error(as_ggplot(d), "Several composites.*select one")
    p <- as_ggplot(d, composite = "Difference", type = "sem")
    payload <- plot_data(plot(d, composite = "Difference", type = "sem", draw = FALSE))
    expect_identical(p$data$Value, payload$series$Value)
    expect_identical(plot_data(p, "title"), payload$title)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
})

test_that("D-study plots hold the other facet fixed and keep nonrectangular scenarios", {
  ratings <- merge(mvgt_plot_data(), data.frame(Rater = 1:3))
  ratings$V <- ratings$V + ratings$Rater
  ratings$W <- ratings$W + 2 * ratings$Rater
  g <- mfrm_g_study(ratings, c("V", "W"))
  grid <- data.frame(Raters = c(2, 4, 2, 4), Tasks = c(3, 3, 9, 6))
  d <- mfrm_d_study(g, grid)
  expect_true(all(d$coefficients$Status == "Available"))
  p <- plot_data(plot(d, draw = FALSE))
  expect_identical(p$design_grid, grid)
  expect_equal(p$series$X, rep(grid$Tasks, 2))
  expect_equal(p$series$Group, rep(paste("Raters =", grid$Raters), 2))
  expect_equal(p$series$Scenario, rep(seq_len(nrow(grid)), 2))
  expect_equal(p$series$Value[p$series$Metric == "G"], d$coefficients$G[d$coefficients$Score == "V"])
  flipped <- plot_data(plot(d, x_var = "Raters", draw = FALSE))
  expect_equal(flipped$series$X, rep(grid$Raters, 2))
  expect_identical(flipped$group_var, "Tasks")
  only_raters <- mfrm_d_study(g, grid[1:2, ])
  expect_identical(plot_data(plot(only_raters, draw = FALSE))$x_var, "Raters")

  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mfrow", "mar", "oma", "cex", "mex", "bg", "fg", "cex.axis"))
  expect_no_warning(plot(d))
  expect_no_warning(plot(d, type = "sem", x_var = "Raters", preset = "monochrome"))
  expect_equal(graphics::par(names(old)), old)
})

test_that("D-study plots preserve independent metric availability and component notes", {
  tasks <- mvgt_plot_data()
  g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
  g$components$Task[,] <- diag(-10, 2)
  d <- mfrm_d_study(g, data.frame(Tasks = c(6, 12)))
  p <- plot_data(plot(d, draw = FALSE))
  expect_true(all(is.finite(p$series$Value[p$series$Metric == "G"])))
  expect_true(all(is.na(p$series$Value[p$series$Metric == "Phi"])))
  expect_true(all(p$unavailable$Metric == "Phi"))
  expect_true(all(p$unavailable$Status == "Negative projected error variance"))
  expect_match(p$component_note, "Non-PSD covariance components")
  sem <- plot_data(plot(d, type = "sem", draw = FALSE))
  expect_true(all(is.finite(sem$series$Value[sem$series$Metric == "RelativeSEM"])))
  expect_true(all(sem$unavailable$Metric == "AbsoluteSEM"))
  expect_match(p$subtitle, "Future complete crossed designs")
  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d))
  expect_no_warning(plot(d, type = "sem"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(d)
    expect_equal(gg$data$Value, p$series$Value)
    expect_match(paste(plot_data(gg, "notes")$Text, collapse = "\n"), "Non-PSD covariance components")
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }

  g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(3, 6, 12)))
  # Exercise a gap in a stored projection without using zero as a placeholder.
  gap <- d$coefficients$Scenario == 2
  d$coefficients[gap, c("G", "Phi", "RelativeSEM", "AbsoluteSEM")] <- NA_real_
  d$coefficients$Status[gap] <- "Negative projected error variance"
  d$coefficients[gap, c("GStatus", "PhiStatus", "RelativeSEMStatus", "AbsoluteSEMStatus")] <-
    "Negative projected error variance"
  p <- plot_data(plot(d, draw = FALSE))
  expect_equal(p$unavailable$Scenario, c(2L, 2L))
  expect_true(all(is.na(p$series$Value[p$series$Scenario == 2])))
  expect_no_warning(plot(d))
  single <- mfrm_d_study(g, data.frame(Tasks = 6))
  expect_no_warning(plot(single))
  expect_error(plot(single, x_var = "Raters"), "facet count present")
  expect_error(plot(single, score = "unknown"), "one original score")
  expect_error(plot(single, score = c("V", "V")), "one original score")
  expect_error(plot(single, draw = NA), "TRUE or FALSE")
  expect_error(plot(single, unused = TRUE), "must be empty")
})

test_that("ggplot conversion retains D-study metrics, units, groups and unavailable gaps", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  ratings <- merge(mvgt_plot_data(), data.frame(Rater = 1:3))
  ratings$V <- ratings$V + ratings$Rater
  ratings$W <- ratings$W + 2 * ratings$Rater
  g <- mfrm_g_study(ratings, c("V", "W"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 4), Tasks = c(3, 6, 12)),
    c(V = .4, W = .6))
  before <- grDevices::dev.cur()
  as_ggplot(d)
  expect_identical(grDevices::dev.cur(), before)
  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (type in c("coefficients", "sem")) {
    before <- grDevices::dev.cur()
    payload <- plot(d, type = type, preset = "monochrome", draw = FALSE)
    p <- as_ggplot(payload)
    expect_identical(grDevices::dev.cur(), before)
    expect_equal(p$data$Value, payload$data$series$Value)
    expect_equal(p$data$X, payload$data$series$X)
    expect_identical(as.character(p$data$Group), payload$data$series$Group)
    expect_identical(levels(p$data$Panel), payload$data$panel_labels)
    expect_identical(p$labels$x, "Number of tasks")
    expect_match(p$labels$y, if (type == "sem") "score units" else "Dependability")
    expect_identical(plot_data(p, "notes"), plot_data(payload, "notes"))
    built <- ggplot2::ggplot_build(p)
    expect_identical(unname(built$plot$scales$get_scales("colour")$map(payload$data$legend$label)),
      payload$data$legend$value)
    expect_equal(sort(built$data[[2]]$y), sort(payload$data$series$Value))
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  direct <- as_ggplot(d, type = "sem", score = "V", x_var = "Raters")
  expect_identical(plot_data(direct, "title"), "Score: V")
  expect_identical(direct$labels$x, "Number of raters")
  expect_true(all(direct$data$Kind == "Score"))

  gaps <- d$coefficients$Tasks == 6
  d$coefficients[gaps, c("G", "Phi", "RelativeSEM", "AbsoluteSEM")] <- NA_real_
  d$coefficients$Status[gaps] <- "Negative projected error variance"
  d$coefficients[gaps, c("GStatus", "PhiStatus", "RelativeSEMStatus", "AbsoluteSEMStatus")] <-
    "Negative projected error variance"
  p <- as_ggplot(d)
  expect_true(all(is.na(p$data$Value[p$data$X == 6])))
  expect_true(all(is.na(p$layers[[1]]$data$Value[p$layers[[1]]$data$X == 6])))
  expect_match(gsub("\n", " ", paste(plot_data(p, "notes")$Text, collapse = "\n")), "Unavailable estimates in 2 scenario")
  expect_no_warning(ggplot2::ggplotGrob(p))

  d$coefficients[c("G", "Phi", "RelativeSEM", "AbsoluteSEM")] <- NA_real_
  d$coefficients$Status <- "Non-PSD component estimates"
  # Saved results using the old global policy keep their recorded gaps.
  d$coefficients[c("GStatus", "PhiStatus", "RelativeSEMStatus", "AbsoluteSEMStatus")] <- NULL
  p <- as_ggplot(d)
  expect_true(all(is.na(p$data$Value)))
  expect_match(p$layers[[3]]$data$Label, "Covariance components failed validity checks")
  expect_match(p$layers[[4]]$data$Label, "Covariance components failed validity checks")
  expect_no_warning(ggplot2::ggplotGrob(p))
  single <- mfrm_d_study(g)
  expect_no_warning(ggplot2::ggplotGrob(as_ggplot(single)))
})

test_that("ggplot preserves separate component and unavailable notes for sparse D studies", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  data <- mvgt_plot_data()
  data <- data[(data$Person + data$Task) %% 3 != 0, ]
  data$V[1] <- NA_real_
  g <- mfrm_g_study(data, c("V", "W"), rater = NULL,
    method = "minque0", missing = "omit")
  d <- mfrm_d_study(g, data.frame(Tasks = c(6, 12)), c(V = -1, W = 1))
  expect_true(any(!g$component_diagnostics$PositiveSemidefinite))
  for (type in c("coefficients", "sem")) {
    p <- plot(d, type = type, draw = FALSE)
    before <- serialize(p, NULL)
    gg <- as_ggplot(p)
    expect_identical(plot_data(gg, "notes"), plot_data(p, "notes"))
    expect_identical(serialize(p, NULL), before)
    expect_identical(gg$data$Value, plot_data(p, "series")$Value)
    expect_null(gg$labels$caption)
    expect_no_warning(ggplot2::ggplotGrob(gg))
    # A saved payload from before returned captions still retains each caution.
    legacy <- p
    legacy$data$caption <- legacy$data$notes <- NULL
    expect_identical(trimws(plot_data(as_ggplot(legacy), "notes")$Text),
      trimws(plot_data(p, "notes")$Text))
  }
})

test_that("rater-only and named-facet plots retain their declared count identities", {
  data <- mvgt_plot_data()
  data$Assessor <- data$Task
  g <- mfrm_g_study(data, c("V", "W"), rater = "Assessor", task = NULL)
  d <- mfrm_d_study(g, data.frame(Raters = c(3, 6, 12)))
  p <- plot_data(plot(d, draw = FALSE))
  expect_identical(p$x_var, "Raters")
  expect_identical(p$x_label, "Number of raters")
  expect_identical(p$legend$label, "Raters only")
  expect_equal(p$series$X, rep(c(3, 6, 12), 2))
  named <- mfrm_g_study(data, c("V", "W"), facets = c("Test occasion" = "Assessor"))
  named_d <- mfrm_d_study(named)
  one <- plot_data(plot(named_d, draw = FALSE))
  expect_identical(one$x_label, "Test occasion count")
  expect_true("Test occasion" %in% names(one$series))

  data <- merge(data, data.frame(Session = 1:3))
  data$V <- data$V + data$Session
  data$W <- data$W + 2 * data$Session
  g <- mfrm_g_study(data, c("V", "W"),
    facets = c(Occasion = "Session", Rater = "Assessor"))
  grid <- expand.grid(Occasion = c(2, 4), Rater = c(3, 6))
  d <- mfrm_d_study(g, grid)
  p <- plot_data(plot(d, draw = FALSE))
  expect_identical(p$x_var, "Rater")
  expect_identical(p$group_var, "Occasion")
  expect_equal(p$series$X, rep(grid$Rater, 2))
  expect_equal(p$series$Group, rep(paste("Occasion =", grid$Occasion), 2))
  flipped <- plot_data(plot(d, x_var = "Occasion", type = "sem", draw = FALSE))
  expect_identical(flipped$x_label, "Occasion count")
  expect_equal(flipped$series$X, rep(grid$Occasion, 2))
  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, x_var = "Occasion", type = "sem"))
  expect_no_warning(plot(named_d))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(d, x_var = "Occasion", type = "sem")
    expect_identical(gg$labels$x, flipped$x_label)
    expect_equal(gg$data$X, flipped$series$X)
    expect_identical(as.character(gg$data$Group), flipped$series$Group)
    expect_no_warning(ggplot2::ggplotGrob(gg))
    expect_identical(as_ggplot(named_d)$labels$x, "Test occasion count")
  }
})


test_that("multiple score curves retain identity, fixed counts, weights and gaps", {
  ratings <- merge(mvgt_plot_data(), data.frame(Rater = 1:3))
  ratings$V <- ratings$V + ratings$Rater
  ratings$W <- ratings$W + 2 * ratings$Rater
  g <- mfrm_g_study(ratings, c("V", "W"))
  w <- cbind(V = c(V = .5, W = .5), Difference = c(V = -1, W = 1))
  d <- mfrm_d_study(g, expand.grid(Raters = 1:4, Tasks = c(3, 6, 12)), w)
  before <- serialize(d, NULL)
  p <- plot_data(plot(d, score = c("W", "V"), metric = "Phi", draw = FALSE))
  expect_identical(p$score, c("W", "V"))
  expect_identical(p$metric, "Phi")
  expect_true(p$show_legend)
  expect_null(p$weights)
  expect_identical(p$table, d$coefficients[d$coefficients$Kind == "Score", ])
  expect_identical(p$series$Value, p$table$Phi)
  expect_identical(p$legend$label, c(paste("W | Raters =", 1:4), paste("V | Raters =", 1:4)))
  curves <- split(p$series, p$series$Group)
  expect_length(curves, 8L)
  expect_true(all(vapply(curves, function(z) nrow(unique(z[c("Score", "Raters", "Metric")])) == 1L,
    logical(1))))
  expect_equal(anyDuplicated(data.frame(p$line_types, p$point_shapes)), 0L)
  expect_match(paste(p$notes$Text, collapse = " "), "equivalent score meaning")
  comp <- plot_data(plot(d, composite = c("Difference", "V"), metric = "G", draw = FALSE))
  expect_identical(comp$weights, w[, c("Difference", "V")])
  expect_identical(comp$table, d$coefficients[d$coefficients$Kind == "Composite", ])
  expect_true(all(comp$series$Kind == "Composite"))
  expect_identical(serialize(d, NULL), before)

  # An unavailable score remains in the table/legend; internal gaps break lines.
  missing <- d$coefficients$Kind == "Score" &
    (d$coefficients$Score == "W" | d$coefficients$Tasks == 6)
  d$coefficients$Phi[missing] <- NA_real_
  d$coefficients$PhiStatus[missing] <- "Negative projected error variance"
  payload <- plot(d, score = c("W", "V"), metric = "Phi", reference = .8, draw = FALSE)
  expect_equal(nrow(payload$data$unavailable), sum(missing))
  expect_equal(nrow(payload$data$series), sum(d$coefficients$Kind == "Score"))
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, score = c("W", "V"), metric = "Phi", reference = .8, preset = "monochrome"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(payload)
    expect_identical(plot_data(gg, "notes"), payload$data$notes)
    expect_identical(as.character(gg$data$Group), payload$data$series$Group)
    expect_true(all(is.na(gg$layers[[1]]$data$Value[gg$layers[[1]]$data$Tasks == 6])))
    expect_length(unique(ggplot2::ggplot_build(gg)$layout$layout$PANEL), 1L)
    expect_identical(gg$theme$legend.position, "bottom")
    # Long display labels must fit the article's seven-inch drawing device.
    long <- payload
    labels <- paste("Detailed criterion", long$data$legend$label)
    long$data$series$Group <- labels[match(long$data$series$Group, long$data$legend$label)]
    long$data$legend$label <- labels
    grob <- ggplot2::ggplotGrob(as_ggplot(long))
    legend <- which(startsWith(grob$layout$name, "guide-box") &
      vapply(grob$grobs, inherits, logical(1), "gtable"))
    expect_length(legend, 1L)
    expect_lte(grid::convertWidth(sum(grob$grobs[[legend]]$widths), "inches", valueOnly = TRUE), 7)
    expect_null(gg$labels$title)
    expect_null(gg$labels$caption)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})

test_that("single metrics and optional references use the requested scale", {
  g <- mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(3, 6, 12)))
  p <- plot_data(plot(d, score = c("W", "V"), metric = "Phi", reference = c(.7, .8), draw = FALSE))
  expect_identical(p$legend$label, c("W", "V"))
  expect_true(p$show_legend)
  expect_identical(p$reference_lines$value, c(.7, .8))
  expect_true(all(p$reference_lines$role == "user_reference"))
  expect_match(paste(p$notes$Text, collapse = " "), "not universal acceptance rules")
  expect_equal(nrow(plot_data(plot(d, draw = FALSE))$reference_lines), 0L)
  for (bad in list(character(), c("G", "Phi"), NA_character_, "AbsoluteSEM", "unknown", 1)) {
    expect_error(plot(d, metric = bad, draw = FALSE), "compatible with `type`")
  }
  for (bad in list(numeric(), NA_real_, Inf, -0.1, 1.1, c(.8, .8), "0.8", TRUE, .8+0i, matrix(.8))) {
    expect_error(plot(d, reference = bad, draw = FALSE), "distinct finite numeric")
  }
  expect_error(plot(d, type = "sem", metric = "Phi", draw = FALSE), "compatible")
  sem <- plot(d, score = c("V", "W"), type = "sem", metric = "AbsoluteSEM", reference = 10, draw = FALSE)
  expect_match(sem$data$y_label, "score units")
  expect_identical(sem$data$series$Value, d$coefficients$AbsoluteSEM)
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  h <- NULL
  original_abline <- graphics::abline
  local_mocked_bindings(abline = function(..., h = NULL) {
    if (!is.null(h)) h <<- h
    original_abline(..., h = h)
  }, .package = "mfrmr")
  expect_no_warning(plot(d, type = "sem", metric = "AbsoluteSEM", reference = 10))
  expect_equal(h, 10)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(sem)
    expect_equal(ggplot2::ggplot_build(gg)$layout$panel_params[[1]]$y$limits, c(0, 10))
    expect_equal(ggplot2::ggplot_build(gg)$data[[3]]$yintercept, 10)
    expect_identical(plot_data(gg, "reference_lines"), sem$data$reference_lines)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})


test_that("facet-count panels share score styles and keep only supplied scenarios", {
  ratings <- merge(mvgt_plot_data(), data.frame(Rater = 1:3))
  ratings$V <- ratings$V + ratings$Rater
  ratings$W <- ratings$W + 2 * ratings$Rater
  g <- mfrm_g_study(ratings, c("V", "W"))
  grid <- data.frame(Raters = c(2, 4, 2, 4, 2), Tasks = c(3, 3, 6, 6, 12))
  weights <- cbind(Equal = c(V = .5, W = .5), Difference = c(V = -1, W = 1))
  d <- mfrm_d_study(g, grid, weights)
  before <- serialize(d, NULL)
  p <- plot(d, score = c("W", "V"), x_var = "Tasks", panel_by = "Raters",
    metric = "Phi", reference = .8, draw = FALSE)
  expect_identical(p$data$panel_by, "Raters")
  expect_equal(p$data$ncol, 2)
  expect_equal(p$data$panels$Count, c(2, 4))
  expect_identical(p$data$panels$Label, c("Raters = 2", "Raters = 4"))
  expect_identical(p$data$table, d$coefficients[d$coefficients$Kind == "Score", ])
  expect_identical(p$data$series$Value, p$data$table$Phi)
  expect_identical(p$data$series$Group, p$data$series$Score)
  expect_identical(p$data$legend$label, c("W", "V"))
  expect_length(unique(p$data$legend$value), 2L)
  expect_equal(nrow(p$data$series), 2 * nrow(grid))
  groups <- split(p$data$series, p$data$series$Panel)
  expect_true(all(vapply(groups, function(z) length(unique(z$Raters)) == 1L, logical(1))))
  expect_equal(sort(unique(groups[[2]]$X)), c(3, 6))
  both <- plot(d, score = "V", panel_by = "Raters", draw = FALSE)
  expect_identical(both$data$panels$Metric, c("G", "G", "Phi", "Phi"))
  expect_false(both$data$show_legend)
  expect_identical(both$data$legend$label, "V")
  single_column <- plot(d, score = "V", panel_by = "Raters", ncol = 1, draw = FALSE)
  expect_equal(single_column$data$ncol, 1)
  comp <- plot(d, composite = c("Equal", "Difference"), metric = "Phi", panel_by = "Raters", draw = FALSE)
  expect_identical(comp$data$weights, weights)
  expect_identical(comp$data$series$Value, comp$data$table$Phi)
  flipped <- plot(d, score = "V", x_var = "Raters", panel_by = "Tasks", draw = FALSE)
  expect_equal(flipped$data$panels$Count, rep(c(3, 6, 12), 2))
  expect_identical(serialize(d, NULL), before)
  for (bad in list("Tasks", "Score", "missing", c("Raters", "Tasks"), NA_character_, 1, character())) {
    expect_error(plot(d, score = "V", x_var = "Tasks", panel_by = bad, draw = FALSE), "non-horizontal facet count")
  }
  for (bad in list(0, -1, 1.5, NA_real_, Inf, "2", TRUE, c(1, 2), numeric(), 1+0i)) {
    expect_error(plot(d, score = "V", panel_by = "Raters", ncol = bad, draw = FALSE), "positive integer")
  }
  expect_error(plot(d, score = "V", ncol = 2, draw = FALSE), "requires `panel_by`")
  one <- mfrm_d_study(mfrm_g_study(mvgt_plot_data(), c("V", "W"), rater = NULL))
  expect_error(plot(one, panel_by = "Tasks", draw = FALSE), "non-horizontal facet count")

  grDevices::pdf(NULL, width = 7, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mfrow", "mar", "oma", "cex", "mex", "fig", "bg", "fg"))
  expect_no_warning(plot(d, score = c("W", "V"), metric = "Phi", panel_by = "Raters", reference = .8))
  expect_no_warning(plot(d, composite = c("Equal", "Difference"), type = "sem", panel_by = "Raters"))
  expect_equal(graphics::par(names(old)), old)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p)
    built <- ggplot2::ggplot_build(gg)
    expect_length(unique(built$layout$layout$PANEL), 2L)
    expect_identical(levels(gg$data$Panel), p$data$panels$Label)
    expect_identical(plot_data(gg, "notes"), p$data$notes)
    expect_identical(gg$data$Value, p$data$series$Value)
    expect_equal(built$data[[3]]$yintercept, rep(.8, 2))
    expect_equal(built$data[[2]]$colour,
      p$data$legend$value[match(p$data$series$Score, p$data$legend$label)])
    expect_length(built$layout$panel_scales_x, 1L)
    expect_length(built$layout$panel_scales_y, 1L)
    expect_equal(built$layout$panel_params[[1]]$x$limits, c(3, 12))
    expect_equal(built$layout$panel_params[[2]]$y$limits, c(0, 1))
    expect_identical(as_ggplot(both)$theme$legend.position, "none")
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})

test_that("empty panels and missing points stay visible with metric-specific reasons", {
  ratings <- merge(mvgt_plot_data(), data.frame(Rater = 1:3))
  g <- mfrm_g_study(ratings, c("V", "W"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 4), Tasks = c(3, 6, 12)))
  missing <- d$coefficients$Raters == 4 | d$coefficients$Tasks == 6
  d$coefficients$Phi[missing] <- NA_real_
  d$coefficients$PhiStatus[missing] <- "Negative projected error variance"
  p <- plot(d, score = c("V", "W"), panel_by = "Raters", draw = FALSE)
  expect_equal(nrow(p$data$panels), 4L)
  expect_true(all(is.finite(p$data$series$Value[p$data$series$Metric == "G"])))
  expect_true(all(p$data$unavailable$Metric == "Phi"))
  expect_true(all(is.na(p$data$series$Value[p$data$series$Panel == 4L])))
  grDevices::pdf(NULL, width = 7, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, score = c("V", "W"), panel_by = "Raters"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(p)
    expect_length(levels(gg$data$Panel), 4L)
    line <- gg$layers[[1]]$data
    gap <- line$Metric == "Phi" & line$Raters == 2 & line$X == 6
    expect_equal(sum(gap), 2L)
    expect_true(all(is.na(line$Value[gap])))
    expect_false(any(line$Metric == "Phi" & line$Raters == 4))
    expect_match(gg$layers[[3]]$data$Label, "Negative projected error")
    expect_identical(as.character(gg$layers[[3]]$data$Panel), p$data$panels$Label[4])
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})

test_that("nested and named-facet panel labels retain their count meaning", {
  env <- new.env(parent = environment())
  eval(parse(test_path("test-multivariate-nested-gtheory.R"), keep.source = FALSE)[[1L]], env)
  ratings <- env$mvgt_nested_fixture()$data
  g <- mfrm_g_study(ratings, c("Content", "Organization"), nesting = c(Rater = "Task"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 3), Tasks = c(4, 8)))
  p <- plot(d, metric = "Phi", panel_by = "Raters", draw = FALSE)
  expect_identical(p$data$panels$Label, c("Raters per Task = 2", "Raters per Task = 3"))
  flipped <- plot(d, x_var = "Raters", panel_by = "Tasks", metric = "Phi", draw = FALSE)
  expect_identical(flipped$data$x_label, "Number of raters per Task")
  named <- mfrm_g_study(ratings, c("Content", "Organization"),
    facets = c("Test occasion" = "Task", Label = "Rater"))
  named_d <- mfrm_d_study(named, expand.grid(`Test occasion` = c(4, 8), Label = c(2, 3)))
  named_p <- plot(named_d, x_var = "Test occasion", panel_by = "Label", metric = "Phi", draw = FALSE)
  expect_identical(named_p$data$panels$Label, c("Label = 2", "Label = 3"))
  expect_equal(named_p$data$panels$Count, c(2, 3))
  expect_identical(named_p$data$series$Label, named_p$data$table$Label)
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, metric = "Phi", panel_by = "Raters"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_identical(levels(as_ggplot(p)$data$Panel), p$data$panels$Label)
    expect_no_warning(ggplot2::ggplotGrob(as_ggplot(named_p)))
  }
})
