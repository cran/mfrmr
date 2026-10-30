facet_reporting_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      ratings <- load_mfrmr_data("example_core")
      cached <<- lapply(c("RSM", "PCM"), function(model) {
        fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
          model = model, step_facet = if (model == "PCM") "Criterion" else NULL)
        ci <- mfrm_facet_intervals(fit, "Rater", method = "sandwich", level = .995)
        contrast <- matrix(c(1, -1, 0, 0), 1,
          dimnames = list("R01 minus R02", ci$table$Target))
        pair <- mfrm_facet_intervals(fit, "Rater", contrasts = contrast, level = .9)
        list(fit = fit, ci = ci, pair = pair)
      })
    }
    cached
  }
})

test_that("RSM and PCM saved intervals reach tables and reports without recomputation", {
  fixtures <- facet_reporting_fixture()
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
    compute_mml_parameter_covariance = function(...) stop("unexpected covariance"))
  for (x in fixtures) {
    res <- mfrm_results(x$fit, intervals = list(Raters = x$ci, difference = x$pair),
      include = c("fit", "plots"), compute = "never")
    expect_identical(res$facet_intervals$raters, x$ci)
    code <- paste(summary(res)$reproducible_code$Code, collapse = "\n")
    expect_identical(code, res$input$reproducible_code)
    expect_match(code, "saveRDS", fixed = TRUE)
    expect_match(code, "readRDS", fixed = TRUE)
    expect_false(grepl("mfrm_results(", code, fixed = TRUE))
    historical <- res; historical$input$reproducible_code <- "stop('outdated reconstruction')"
    expect_identical(paste(summary(historical)$reproducible_code$Code, collapse = "\n"), code)
    payload <- mfrm_results_viewer_payload(historical)
    expect_identical(payload$replay_code, code)
    expect_true("facet_raters" %in% payload$plot_choices)
    expect_true("facet_raters_intervals" %in% payload$table_names)
    html <- mfrm_results_html(historical)$html
    expect_match(html, "saveRDS", fixed = TRUE)
    expect_false(grepl("outdated reconstruction", html, fixed = TRUE))
    folder <- tempfile("displayed-replay-"); dir.create(folder)
    env <- list2env(list(res = res), parent = globalenv())
    withr::with_dir(folder, capture.output(eval(parse(text = code), envir = env)))
    expect_identical(env$res, res)
    unlink(folder, recursive = TRUE)
    expect_identical(res$tables$facet_raters_intervals[names(x$ci$table)], x$ci$table)
    expect_true(all(res$tables$facet_raters_intervals$ConfidenceLevel == "99.5%"))
    expect_true(all(res$tables$facet_difference_intervals$ConfidenceLevel == "90%"))
    expect_equal(res$tables$facet_raters_clusters, x$ci$clusters)
    expect_equal(as.matrix(res$tables$facet_difference_contrasts[-1]), x$pair$contrasts,
      ignore_attr = TRUE)
    expect_equal(apa_table(x$ci, digits = 10)$table$Lower, x$ci$table$Lower, tolerance = 1e-9)
    expect_identical(apa_table(x$ci)$table$ConfidenceLevel, rep("99.5%", 4))
    expect_identical(apa_table(x$ci, which = "settings")$table$ConfidenceLevel, "99.5%")
    expect_identical(apa_table(x$ci, which = "clusters")$table, x$ci$clusters)
    expect_error(apa_table(x$ci, which = "typo"), "Choose a fixed-facet table")
    report <- mfrm_report(res)
    expect_identical(report$tables$facet_raters_intervals, res$tables$facet_raters_intervals)
    expect_match(report$markdown, "Fixed-facet uncertainty")
    expect_match(report$markdown, "99.5%", fixed = TRUE)
    expect_identical(plot_data(res, type = "facet_difference", component = "table"), x$pair$table)
    expect_s3_class(as_ggplot(res, type = "facet_raters"), "ggplot")
    single <- mfrm_results(x$fit, intervals = x$ci, include = "fit", compute = "never")
    expect_identical(single$facet_intervals$inference, x$ci)
    expect_error(as_ggplot(single, type = "facet_inference"), "not available")
  }
})

test_that("saved facet intervals reject mismatched analyses and ambiguous names", {
  x <- facet_reporting_fixture()[[1]]
  for (change in list(
    function(f) { f$opt$par[1] <- f$opt$par[1] + .01; f },
    function(f) { f$prep$data$Score[1] <- 99; f },
    function(f) { f$config$estimation_control$quad_points <- 81; f },
    function(f) { f$config$noncenter_facet <- "Rater"; f },
    function(f) { f$config$model <- "PCM"; f })) {
    expect_error(mfrm_results(change(x$fit), intervals = x$ci, compute = "never"), "must match")
  }
  expect_error(mfrm_results(x$fit, intervals = list(A = x$ci, a = x$ci)), "distinct")
  expect_error(mfrm_results(x$fit, intervals = list(x$ci)), "named list")
  expect_error(mfrm_results(x$fit, intervals = list()), "named list")
  expect_error(mfrm_results(x$fit, intervals = list(wrong = x$fit)), "named list")
})

test_that("both plotting routes retain selected and unavailable outcomes with optional text", {
  x <- facet_reporting_fixture()[[1]]
  groups <- data.frame(Person = x$fit$prep$levels$Person, Cluster = rep(1:2, each = 24))
  missing <- mfrm_facet_intervals(x$fit, "Rater", method = "sandwich", clusters = groups)
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
    compute_mml_parameter_covariance = function(...) stop("unexpected covariance"))
  p <- as_ggplot(missing, title = NULL, subtitle = NULL, caption = NULL,
    reference = NULL, show_legend = FALSE, preset = "monochrome")
  expect_null(p$labels$title); expect_null(p$labels$subtitle); expect_null(p$labels$caption)
  expect_identical(p$theme$legend.position, "none")
  expect_identical(plot_data(p, "table"), missing$table)
  expect_true(all(is.na(plot_data(p, "table")$Lower)))
  expect_equal(sum(p$data$Display == "Unavailable"), 4)
  expect_equal(sum(p$data$Display == "Available"), 4)
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(as_ggplot(missing, component = "table"), "plot_data")
  expect_error(as_ggplot(plot(missing, draw = FALSE), component = "table"), "plot_data")
  expect_match(paste(plot_data(as_ggplot(x$ci), "notes")$Text, collapse = "\n"), "99.5%", fixed = TRUE)
  expect_false(plot_data(as_ggplot(x$ci, comparison = FALSE), "comparison"))
  expect_error(plot(x$ci, reference = Inf, draw = FALSE), "finite")
  expect_error(plot(x$ci, title = NA_character_, draw = FALSE), "Plot text")
  expect_error(plot(x$ci, show_legend = NA, draw = FALSE), "TRUE or FALSE")
  path <- tempfile(fileext = ".pdf"); grDevices::pdf(path)
  plot(missing, title = NULL, subtitle = NULL, caption = NULL, reference = NULL, show_legend = FALSE)
  grDevices::dev.off()
  expect_gt(file.info(path)$size, 0)
})

test_that("export replay preserves RSM and PCM selections instead of refitting", {
  fixtures <- facet_reporting_fixture()
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
    compute_mml_parameter_covariance = function(...) stop("unexpected covariance"))
  for (x in fixtures) {
    res <- mfrm_results(x$fit, intervals = list(raters = x$ci, difference = x$pair),
      include = c("fit", "plots"), compute = "never")
    folder <- tempfile("facet-export-")
    warnings <- character()
    exported <- withCallingHandlers(export_mfrm_results(res, output_dir = folder,
      include = c("tables", "plots", "report", "replay"), preset = "starter", acknowledge_sensitive = TRUE),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    expect_true(all(startsWith(warnings, "Review-only display:")))
    files <- list.files(folder, recursive = TRUE, full.names = TRUE)
    saved <- readRDS(files[grepl("\\.rds$", files)][1])
    expect_identical(saved$facet_intervals, res$facet_intervals)
    expect_true(any(grepl("facet_raters.*\\.png$", files)))
    index <- paste(readLines(file.path(folder, "index.html")), collapse = "\n")
    expect_match(index, "Saved inference", fixed = TRUE)
    expect_match(index, 'src="mfrmr_results_plot_facet_raters.png"', fixed = TRUE)
    expect_match(index, "99.5% pointwise", fixed = TRUE)
    expect_match(index, "does not replace", fixed = TRUE)
    expect_true(any(grepl("99.5% pointwise", exported$written_files$Note, fixed = TRUE)))
    expect_false(any(grepl("facet_", exported$plot_errors$Plot)))
    replay <- files[grepl("replay.*\\.R$", files)]
    expect_length(replay, 1)
    expect_match(paste(readLines(replay), collapse = "\n"), "readRDS")
    expect_no_error(source(replay, local = new.env(), chdir = TRUE))
    unlink(folder, recursive = TRUE)
  }
})
