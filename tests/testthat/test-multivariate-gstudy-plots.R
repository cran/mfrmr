mvgs_plot_fixture <- function(nested = FALSE) {
  e <- new.env(parent = environment())
  file <- if (nested) 'test-multivariate-nested-gtheory.R' else 'test-multivariate-gtheory.R'
  eval(parse(test_path(file), keep.source = FALSE)[[1L]], e)
  if (nested) e$mvgt_nested_fixture()$data else e$mvgt_fixture()$data
}

test_that('component plots agree with direct weighted-score estimation across designs', {
  for (nested in c(FALSE, TRUE)) for (incomplete in c(FALSE, TRUE)) {
    data <- mvgs_plot_fixture(nested)
    if (incomplete) data <- data[-c(1, 10, 35), ]
    nesting <- if (nested) c(Rater = 'Task') else NULL
    method <- if (incomplete) 'minque0' else 'anova'
    g <- mfrm_g_study(data, c('Content', 'Organization'), method = method, nesting = nesting)
    before <- serialize(g, NULL)
    weights <- c(Organization = -.7, Content = .3)
    p <- plot(g, weights = weights, draw = FALSE)
    data$Total <- .3 * data$Content - .7 * data$Organization
    direct <- mfrm_g_study(data, 'Total', method = method, nesting = nesting)
    expected <- vapply(direct$components, function(a) a[1, 1], numeric(1))
    expect_equal(p$data$table$Estimate, unname(expected), tolerance = 1e-11)
    expect_identical(p$data$table$Source, names(g$components))
    expect_identical(p$data$table$ComponentPSD, g$component_diagnostics$PositiveSemidefinite)
    expect_identical(p$data$weights, weights[g$design$scores])
    expect_identical(p$data$source_design$nesting, nesting)
    expect_identical(p$data$source_rows, g$data_usage$counts)
    expect_false('levels' %in% names(p$data$source_design))
    scaled <- plot(g, weights = -3 * weights, draw = FALSE)
    expect_equal(scaled$data$table$Estimate, 9 * p$data$table$Estimate, tolerance = 1e-12)
    expect_identical(serialize(g, NULL), before)
    one <- plot(g, score = 'Organization', draw = FALSE)
    expect_equal(one$data$table$Estimate, unname(vapply(g$components, function(a) a[2, 2], numeric(1))))
  }
})

test_that('component plots preserve negative values, zero and covariance warnings', {
  tasks <- read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr', mustWork = TRUE))
  tasks <- tasks[(tasks$Person + tasks$Task) %% 3 != 0, ]
  tasks$V[1] <- NA_real_
  g <- mfrm_g_study(tasks, c('V', 'W'), rater = NULL, method = 'minque0', missing = 'omit')
  p <- plot(g, weights = c(V = -1, W = 1), draw = FALSE)
  expect_lt(p$data$table$Estimate[1], 0)
  expect_identical(p$data$table$Status[1], 'Negative estimate')
  expect_match(paste(p$data$notes$Text, collapse = ' '), 'Negative raw estimates')
  expect_match(paste(p$data$notes$Text, collapse = ' '), 'Non-PSD')
  # Positive selected projections must not hide a non-PSD full covariance matrix.
  positive <- plot(g, score = 'V', draw = FALSE)
  expect_true(any(!positive$data$table$ComponentPSD))
  expect_match(paste(positive$data$notes$Text, collapse = ' '), 'positive selected variance')
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  before <- graphics::par(c('mfrow', 'mar', 'oma', 'cex', 'mex'))
  expect_no_warning(plot(g, weights = c(V = -1, W = 1), preset = 'monochrome'))
  expect_equal(graphics::par(names(before)), before)
  g$components <- lapply(g$components, function(a) {a[] <- 0; a})
  expect_no_warning(zero <- plot(g))
  expect_identical(zero$data$table$Estimate, rep(0, 3))
})

test_that('G-study plot selection and malformed saved values have actionable errors', {
  g <- mfrm_g_study(mvgs_plot_fixture(), c('Content', 'Organization'),
    facets = c(Assessor = 'Rater', Occasion = 'Task'))
  expect_identical(plot(g, draw = FALSE)$data$score, 'Content')
  expect_true('Assessor:Occasion' %in% plot(g, draw = FALSE)$data$table$Source)
  expect_error(plot(g, score = 'absent'), 'original score')
  expect_error(plot(g, score = 'Content', weights = c(Content = 1, Organization = 1)), 'either')
  for (w in list(c(1, 1), c(Content = 1), c(Content = NA, Organization = 1),
    c(Content = 0, Organization = 0), c(Content = 1, Content = 1),
    matrix(1, 2, 1), c(Content = 1i, Organization = 1))) {
    expect_error(plot(g, weights = w), 'weights.*name every score')
  }
  expect_error(plot(g, weights = c(Content = 1e200, Organization = 1)), 'rescale')
  expect_error(plot(g, draw = NA), 'TRUE or FALSE')
  expect_error(plot(g, unused = TRUE), 'must be empty')
  altered <- g; altered$components$Person[1, 2] <- NA_real_
  expect_error(plot(altered), 'incomplete or altered')
  altered <- g; altered$component_diagnostics <- NULL
  expect_error(plot(altered), 'diagnostics.*incomplete')
  legacy <- g; legacy$design$count_columns <- legacy$design$facets <- NULL
  expect_identical(plot(legacy, draw = FALSE)$data$table, plot(g, draw = FALSE)$data$table)
})

test_that('G-study ggplot and saved payloads retain signed values and external notes', {
  skip_if_not_installed('ggplot2', minimum_version = '3.4.0')
  g <- mfrm_g_study(mvgs_plot_fixture(), c('Content', 'Organization'))
  g$components$Person[1, 1] <- -.25
  g$component_diagnostics$PositiveSemidefinite[1] <- FALSE
  p <- plot(g, score = 'Content', preset = 'monochrome', draw = FALSE)
  path <- tempfile(fileext = '.rds'); on.exit(unlink(path))
  saveRDS(p, path)
  local_mocked_bindings(mfrm_g_study = function(...) stop('No re-estimation'),
    mfrm_d_study = function(...) stop('No projection'), .package = 'mfrmr')
  before <- grDevices::dev.cur()
  gg <- as_ggplot(readRDS(path))
  expect_identical(grDevices::dev.cur(), before)
  expect_identical(gg$data$Estimate, p$data$table$Estimate)
  expect_identical(plot_data(gg, 'notes'), p$data$notes)
  expect_identical(plot_data(gg, 'table'), p$data$table)
  expect_identical(gg$scales$get_scales('y')$labels, p$data$labels)
  expect_null(gg$labels$title); expect_null(gg$labels$caption); expect_null(gg$labels$subtitle)
  built <- ggplot2::ggplot_build(gg)
  expect_equal(built$data[[2]]$xend, p$data$table$Estimate)
  expect_equal(built$data[[3]]$shape, ifelse(p$data$table$Estimate < 0, 1, 16))
  expect_no_warning(ggplot2::ggplotGrob(gg))
  expect_identical(as_ggplot(g, score = 'Content')$data$Estimate, gg$data$Estimate)
  expect_error(as_ggplot(p, component = 'weights'), 'table')
  expect_error(as_ggplot(p, score = 'Organization'), 'must be empty')
})

test_that('interaction labels expose the residual confounding for each supported design', {
  crossed <- mfrm_g_study(mvgs_plot_fixture(), c('Content', 'Organization'),
    facets = c(Assessor = 'Rater', Occasion = 'Task'))
  nested <- mfrm_g_study(mvgs_plot_fixture(TRUE), c('Content', 'Organization'),
    nesting = c(Rater = 'Task'))
  tasks <- read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr', mustWork = TRUE))
  one <- mfrm_g_study(tasks, c('V', 'W'), rater = NULL)
  rater <- mfrm_g_study(tasks, c('V', 'W'), rater = 'Task', task = NULL)
  expected <- list(
    c('Person', 'Assessor', 'Occasion', 'Person x Assessor', 'Person x Occasion',
      'Assessor x Occasion', 'Person x Assessor x Occasion + residual'),
    c('Person', 'Task', 'Rater(Task)', 'Person x Task', 'Person x Rater(Task) + residual'),
    c('Person', 'Task', 'Person x Task + residual'),
    c('Person', 'Rater', 'Person x Rater + residual'))
  for (i in seq_along(expected)) {
    g <- list(crossed, nested, one, rater)[[i]]
    p <- plot(g, draw = FALSE)
    expect_identical(p$data$table$Label, expected[[i]])
    expect_identical(p$data$table$Source, names(g$components))
    expect_equal(p$data$table$Estimate, unname(vapply(g$components, function(a) a[1, 1], numeric(1))))
    expect_match(paste(p$data$notes$Text, collapse = ' '), 'cannot be separated with one observation')
    expect_match(paste(p$data$notes$Text, collapse = ' '), 'Relative error.*absolute error')
  }
  altered <- crossed; altered$design$counts <- NULL
  expect_error(plot(altered, draw = FALSE), 'factor identities')
  altered <- nested; altered$design$nesting <- NULL
  expect_error(plot(altered, draw = FALSE), 'component identities.*design')
})

test_that('score panels and covariance cells retain exact source estimates across designs', {
  for (nested in c(FALSE, TRUE)) for (incomplete in c(FALSE, TRUE)) {
    data <- mvgs_plot_fixture(nested)
    if (incomplete) data <- data[-c(1, 10, 35), ]
    g <- mfrm_g_study(data, c('Content', 'Organization'),
      method = if (incomplete) 'minque0' else 'anova',
      nesting = if (nested) c(Rater = 'Task') else NULL)
    before <- serialize(g, NULL)
    scores <- c('Organization', 'Content')
    points <- plot(g, score = scores, ncol = 1, draw = FALSE)
    expect_identical(points$data$panels$Label, scores)
    expect_equal(points$data$ncol, 1)
    for (s in scores) expect_equal(
      subset(points$data$table, Score == s)$Estimate,
      unname(vapply(g$components, function(a) a[s, s], numeric(1))))
    heat <- plot(g, score = scores, type = 'covariance', draw = FALSE)
    tab <- heat$data$table
    expect_equal(nrow(tab), 4L * length(g$components))
    expect_identical(unique(tab$Source), names(g$components))
    expected <- vapply(seq_len(nrow(tab)), function(i)
      g$components[[tab$Source[i]]][tab$Row[i], tab$Column[i]], numeric(1))
    expect_identical(tab$Estimate, expected)
    expect_identical(tab$Element == 'Variance', tab$Row == tab$Column)
    expect_identical(tab$ComponentPSD, rep(g$component_diagnostics$PositiveSemidefinite, each = 4))
    expect_identical(heat$data$source_rows, g$data_usage$counts)
    expect_identical(heat$data$source_design$nesting, g$design$nesting)
    expect_match(tail(heat$data$panels$Label, 1), 'residual')
    expect_null(heat$data$weights)
    expect_identical(serialize(g, NULL), before)
  }
})

test_that('covariance views distinguish negative covariance and full matrix diagnostics', {
  g <- mfrm_g_study(mvgs_plot_fixture(), c('Content', 'Organization'))
  g$components$Person[,] <- matrix(c(1, -2, -2, 1), 2)
  g$component_diagnostics$PositiveSemidefinite[1] <- FALSE
  p <- plot(g, type = 'covariance', draw = FALSE)
  person <- subset(p$data$table, Source == 'Person')
  expect_identical(person$Estimate, c(1, -2, -2, 1))
  expect_false(any(person$ComponentPSD))
  expect_match(paste(p$data$notes$Text, collapse = ' '), 'Negative covariance can be legitimate')
  expect_match(paste(p$data$notes$Text, collapse = ' '), 'Non-PSD')
  subset <- plot(g, score = 'Content', type = 'covariance', draw = FALSE)
  expect_equal(nrow(subset$data$table), length(g$components))
  expect_false(subset$data$table$ComponentPSD[1])
  expect_match(paste(subset$data$notes$Text, collapse = ' '), 'scores not selected')
  expect_error(plot(g, type = 'covariance', weights = c(Content = .5, Organization = .5)), 'omit.*weights')
  for (s in list(character(), c('Content', 'Content'), c('Content', NA), matrix('Content')))
    expect_error(plot(g, score = s, draw = FALSE), 'distinct original scores')
  for (nc in list(0, -1, 1.5, NA, Inf, c(1, 2), '2', 1i, matrix(2)))
    expect_error(plot(g, ncol = nc, draw = FALSE), 'positive integer')
  expect_error(plot(g, type = 'correlation', draw = FALSE), 'arg')
  g$components <- lapply(g$components, function(a) {a[] <- 0; a})
  grDevices::pdf(NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  before <- graphics::par(c('mfrow', 'mar', 'oma', 'cex', 'mex'))
  expect_no_warning(zero <- plot(g, type = 'covariance', preset = 'monochrome'))
  expect_equal(zero$data$fill_limits, c(-1, 1))
  expect_true(all(zero$data$table$Estimate == 0))
  expect_equal(graphics::par(names(before)), before)
  expect_no_warning(plot(g, score = g$design$scores))
})

test_that('new G-study views reopen without estimation and render the stored panels', {
  skip_if_not_installed('ggplot2', minimum_version = '3.4.0')
  g <- mfrm_g_study(mvgs_plot_fixture(), c('Content', 'Organization'))
  plots <- list(plot(g, score = c('Organization', 'Content'), draw = FALSE),
    plot(g, type = 'covariance', preset = 'monochrome', draw = FALSE))
  local_mocked_bindings(mfrm_g_study = function(...) stop('No re-estimation'),
    mfrm_d_study = function(...) stop('No projection'), .package = 'mfrmr')
  for (p in plots) {
    path <- tempfile(fileext = '.rds')
    saveRDS(p, path); saved <- readRDS(path); unlink(path)
    gg <- as_ggplot(saved)
    expect_identical(plot_data(gg, 'table'), p$data$table)
    expect_identical(plot_data(gg, 'notes'), p$data$notes)
    expect_identical(plot_data(gg, 'panels'), p$data$panels)
    expect_null(gg$labels$title); expect_null(gg$labels$caption)
    built <- ggplot2::ggplot_build(gg)
    expect_equal(nrow(built$layout$layout), nrow(p$data$panels))
    expect_length(unique(built$layout$layout$SCALE_X), 1L)
    expect_length(unique(built$layout$layout$SCALE_Y), 1L)
    expect_no_warning(ggplot2::ggplotGrob(gg))
    if (p$data$plot == 'covariance') {
      expect_identical(as.character(gg$data$Row), p$data$table$Row)
      expect_identical(as.character(gg$data$Column), p$data$table$Column)
      expect_identical(gg$data$Estimate, p$data$table$Estimate)
      expect_identical(gg$scales$get_scales('fill')$limits, p$data$fill_limits)
      expect_identical(gg$data$CellLabel, trimws(formatC(p$data$table$Estimate, format = 'g', digits = 3)))
    } else expect_equal(built$data[[2]]$xend, p$data$table$Estimate)
  }
  legacy <- plot(g, draw = FALSE)
  legacy$data[c('plot', 'panels', 'ncol')] <- NULL
  legacy$data$table[c('Panel', 'Score')] <- NULL
  gg <- as_ggplot(legacy)
  expect_equal(ggplot2::ggplot_build(gg)$data[[2]]$xend, legacy$data$table$Estimate)
})

test_that('component bars reuse interval tables for each score and supported design', {
  tasks <- read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr', mustWork = TRUE))
  studies <- list(mfrm_g_study(tasks, c('V', 'W'), rater = NULL),
    mfrm_g_study(mvgs_plot_fixture(), c('Content', 'Organization')),
    mfrm_g_study(mvgs_plot_fixture(TRUE), c('Content', 'Organization'),
      nesting = c(Rater = 'Task')))
  grDevices::pdf(NULL, width = 8, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  local_mocked_bindings(mfrm_g_study = function(...) stop('No refitting'),
    .package = 'mfrmr')
  for (g in studies) {
    before <- serialize(g, NULL)
    scores <- rev(g$design$scores)
    expect_no_warning(p <- plot(g, score = scores, intervals = TRUE,
      assumption = 'normal', level = .90))
    points <- plot(g, score = scores, draw = FALSE)
    expect_identical(p$data$table[names(points$data$table)], points$data$table)
    expect_identical(names(p$data$intervals), scores)
    expect_identical(p$data$table$Source, rep(names(g$components), 2))
    for (s in scores) {
      ci <- confint(g, score = s, assumption = 'normal', level = .90)
      expect_identical(p$data$intervals[[s]], ci)
      rows <- p$data$table[p$data$table$Score == s, ]
      for (column in c('SE', 'Lower', 'Upper', 'RawLower', 'RawUpper', 'Method'))
        expect_identical(rows[[column]], ci[[column]])
      expect_identical(rows$IntervalStatus, ci$Status)
    }
    weights <- stats::setNames(c(-.7, .3), scores)
    composite <- plot(g, weights = weights, intervals = TRUE,
      assumption = 'normal', draw = FALSE)
    expect_identical(composite$data$intervals$Composite,
      confint(g, weights = weights, assumption = 'normal'))
    expect_identical(serialize(g, NULL), before)
  }
})

test_that('component interval views keep empty and degenerate intervals visible', {
  d <- expand.grid(Person = 1:3, Task = 1:3)
  d$Score <- as.vector(outer(c(-1, 0, 1), c(-1, 0, 1)))
  grDevices::pdf(NULL, width = 7, height = 4.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (zero in c(FALSE, TRUE)) {
    if (zero) d$Score <- 5
    g <- mfrm_g_study(d, 'Score', rater = NULL)
    old <- graphics::par(c('mfrow', 'mar', 'oma', 'cex', 'mex'))
    expect_no_warning(p <- plot(g, intervals = TRUE, assumption = 'normal'))
    expect_equal(graphics::par(names(old)), old)
    missing <- p$data$unavailable_intervals
    expect_equal(nrow(missing), if (zero) 3 else 2)
    expect_true(all(is.na(missing$Lower) & is.na(missing$Upper)))
    expect_identical(missing$IntervalStatus, rep(if (zero) 'Degenerate mean squares' else
      'Empty nonnegative interval', nrow(missing)))
    expect_true(all(if (zero) missing$Estimate == 0 else missing$Estimate < 0))
    if (!zero) expect_true(all(missing$RawUpper < 0))
    expect_match(paste(p$data$notes$Text, collapse = ' '), 'Square outlines.*without intervals')
    expect_match(paste(p$data$notes$Text, collapse = ' '), 'not simultaneous')
    if (requireNamespace('ggplot2', quietly = TRUE)) {
      gg <- as_ggplot(p)
      built <- ggplot2::ggplot_build(gg)
      expect_equal(nrow(built$data[[2]]), sum(p$data$table$IntervalAvailable))
      expect_equal(built$data[[3]]$x, missing$Estimate)
      expect_true(all(built$data[[3]]$shape == 0))
      expect_equal(built$data[[4]]$shape, ifelse(p$data$table$Estimate < 0, 1, 16))
      expect_no_warning(ggplot2::ggplotGrob(gg))
    }
  }
})

test_that('saved component interval panels render without recalculation', {
  skip_if_not_installed('ggplot2', minimum_version = '3.4.0')
  tasks <- read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr', mustWork = TRUE))
  g <- mfrm_g_study(tasks, c('V', 'W'), rater = NULL)
  p <- plot(g, score = c('W', 'V'), intervals = TRUE, assumption = 'normal', draw = FALSE)
  path <- tempfile(fileext = '.rds'); on.exit(unlink(path), add = TRUE)
  saveRDS(p, path)
  local_mocked_bindings(confint.mfrm_multivariate_gstudy = function(...) stop('No interval recalculation'),
    mfrm_g_study = function(...) stop('No refitting'), .package = 'mfrmr')
  gg <- as_ggplot(readRDS(path))
  for (component in c('table', 'notes', 'intervals', 'unavailable_intervals', 'panels'))
    expect_identical(plot_data(gg, component), p$data[[component]])
  expect_null(gg$labels$title); expect_null(gg$labels$caption)
  built <- ggplot2::ggplot_build(gg)
  expect_equal(nrow(built$layout$layout), 2)
  expect_identical(as.character(built$layout$layout$Panel), c('W', 'V'))
  expect_equal(built$data[[2]]$xmin, p$data$table$Lower[p$data$table$IntervalAvailable])
  expect_equal(built$data[[2]]$xmax, p$data$table$Upper[p$data$table$IntervalAvailable])
  expect_length(unique(built$layout$layout$SCALE_X), 1L)
  expect_length(unique(built$layout$layout$SCALE_Y), 1L)
  expect_no_warning(ggplot2::ggplotGrob(gg))
})

test_that('component interval requests reject unsupported scope and unused settings', {
  tasks <- read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr', mustWork = TRUE))
  g <- mfrm_g_study(tasks, c('V', 'W'), rater = NULL)
  for (value in list(NULL, NA, 1, c(TRUE, FALSE), list(TRUE)))
    expect_error(plot(g, intervals = value), 'intervals.*TRUE or FALSE')
  expect_error(plot(g, intervals = TRUE), 'assumption')
  expect_error(plot(g, intervals = TRUE, assumption = 'normal', level = 1), 'level')
  expect_error(plot(g, assumption = 'normal'), 'only when')
  expect_error(plot(g, level = .90), 'only when')
  expect_error(plot(g, type = 'covariance', intervals = TRUE, assumption = 'normal'), 'explicitly select')
  incomplete <- mfrm_g_study(tasks[-1, ], c('V', 'W'), rater = NULL, method = 'minque0')
  expect_error(plot(incomplete, intervals = TRUE, assumption = 'normal'), 'MINQUE')
  altered <- g; altered$mean_products <- NULL
  expect_error(plot(altered, intervals = TRUE, assumption = 'normal'), 'matrices')
  expect_no_error(plot(altered, draw = FALSE))
})
