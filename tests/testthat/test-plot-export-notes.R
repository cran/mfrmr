test_that("fit bundle exports preserve all three images and their own notes", {
  fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  diagnostic_calls <- 0L
  fail <- function(...) stop("No fitting or information during plot export")
  local_mocked_bindings(fit_mfrm = fail, compute_mml_parameter_covariance = fail,
    mfrm_native_location_inference = fail, simulate_mfrm_data = fail,
    diagnose_mfrm = function(...) { diagnostic_calls <<- diagnostic_calls + 1L; fail() },
    .package = "mfrmr")
  res <- mfrm_results(fit, include = c("fit", "plots"), compute = "never")
  res$plot_map <- subset(res$plot_map, Type == "fit")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  device <- grDevices::dev.cur()
  grDevices::png(file.path(folder, "reference-%03d.png"), width = 1008, height = 720, res = 144)
  bundle <- suppressWarnings(plot(res, type = "fit"))
  grDevices::dev.off()
  reference <- list.files(folder, pattern = "^reference-", full.names = TRUE)
  expect_length(reference, 3L)
  route <- getS3method("plot", "mfrm_results"); calls <- 0L
  local_mocked_bindings(plot.mfrm_results = function(...) {
    calls <<- calls + 1L
    route(...)
  }, .package = "mfrmr")
  exported <- suppressWarnings(export_mfrm_results(res, file.path(folder, "export"),
    include = c("plots", "manifest", "rds"), acknowledge_sensitive = TRUE,
    plot_width = 1008, plot_height = 720, plot_res = 144))
  expect_identical(calls, 1L)
  expect_equal(nrow(exported$plot_errors), 0L)
  expect_identical(grDevices::dev.cur(), device)
  paths <- file.path(folder, "export", paste0("mfrmr_results_plot_fit_", names(bundle), ".png"))
  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), unname(tools::md5sum(reference)))
  expect_false(file.exists(file.path(folder, "export", "mfrmr_results_plot_fit.png")))
  html <- paste(mfrm_results_export_index_html("mfrmr_results", exported$written_files,
    exported$plot_errors), collapse = "\n")
  for (i in seq_along(bundle)) {
    note_file <- sub("\\.png$", "_notes.csv", paths[i])
    notes <- read.csv(note_file, stringsAsFactors = FALSE)
    expect_identical(notes$Text[notes$Type == "title"], plot_data(bundle[[i]], "title"))
    expect_equal(notes[notes$Type != "title", ], plot_data(bundle[[i]], "notes"), ignore_attr = TRUE)
    expect_match(html, basename(paths[i]), fixed = TRUE)
    expect_match(html, basename(note_file), fixed = TRUE)
  }
  manifest <- read.csv(file.path(folder, "export", "mfrmr_results_written_files.csv"))
  expect_true(all(normalizePath(paths, winslash = "/") %in% manifest$Path))
  restored <- readRDS(file.path(folder, "export", "mfrmr_results_results.rds"))
  again <- suppressWarnings(export_mfrm_results(restored, file.path(folder, "again"),
    include = "plots", acknowledge_sensitive = TRUE,
    plot_width = 1008, plot_height = 720, plot_res = 144))
  expect_identical(calls, 2L)
  expect_equal(nrow(again$plot_errors), 0L)
  expect_identical(unname(tools::md5sum(paths)),
    unname(tools::md5sum(file.path(folder, "again", basename(paths)))))
  expect_identical(diagnostic_calls, 0L)
})

test_that("result pathway plots respect disabled diagnostic computation", {
  .mfrmr_muffle_expected_warnings({
  fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  res <- mfrm_results(fit, include = c("fit", "plots"), compute = "never")
  calls <- 0L
  local_mocked_bindings(diagnose_mfrm = function(...) {
    calls <<- calls + 1L
    stop("Unexpected diagnostic calculation")
  }, .package = "mfrmr")
  reference <- plot(fit, type = "pathway", draw = FALSE, include_fit_measures = FALSE)
  for (type in c("pathway", "fit")) {
    p <- plot(res, type = type, draw = FALSE)
    if (type == "fit") p <- p$pathway_map
    expect_identical(p$data$expected, reference$data$expected)
    expect_identical(p$data$fit_measure_status$Status, "not_requested")
    expect_equal(nrow(p$data$fit_measures), 0L)
    expect_error(plot(res, type = type, draw = FALSE, include_fit_measures = TRUE),
      "compute = 'never'.*matching diagnostics")
  }
  expect_identical(calls, 0L)
  }, "^Review-only display:")
})

test_that("result bundles and pathways reuse matching stored diagnostics", {
  .mfrmr_muffle_expected_warnings({
  fit <- make_toy_fit()
  diagnostics <- make_toy_diagnostics(fit)
  reference <- plot(fit, type = "pathway", diagnostics = diagnostics, draw = FALSE)
  calls <- 0L
  local_mocked_bindings(diagnose_mfrm = function(...) {
    calls <<- calls + 1L
    stop("Unexpected diagnostic recalculation")
  }, .package = "mfrmr")
  for (compute in c("auto", "never")) {
    res <- mfrm_results(fit, diagnostics = diagnostics,
      include = c("fit", "plots"), compute = compute)
    for (type in c("pathway", "fit")) {
      p <- plot(res, type = type, draw = FALSE)
      if (type == "fit") p <- p$pathway_map
      for (field in c("fit_measures", "fit_status", "curve_fit_status", "fit_measure_status"))
        expect_identical(p$data[[field]], reference$data[[field]])
      bad <- res; bad$diagnostics$measures$Estimate[1] <- bad$diagnostics$measures$Estimate[1] + 1
      expect_error(plot(bad, type = type, draw = FALSE), "fit/diagnostics mismatch")
    }
  }
  res <- mfrm_results(fit, include = c("fit", "plots"), compute = "never")
  supplied <- plot(res, type = "pathway", diagnostics = diagnostics,
    include_fit_measures = TRUE, draw = FALSE)
  expect_identical(supplied$data$fit_measures, reference$data$fit_measures)
  expect_identical(calls, 0L)
  }, "^Review-only display:")
})

test_that("full Wright axis labels stay on the device when abbreviations collide", {
  fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  draw_labels <- draw_rotated_x_labels
  bounds <- list()
  local_mocked_bindings(draw_rotated_x_labels = function(at, labels, srt, cex, line_offset) {
    usr <- graphics::par("usr")
    anchor <- graphics::grconvertY(usr[3] - line_offset * diff(usr[3:4]), "user", "inches")
    height <- abs(sin(srt * pi / 180)) * graphics::strwidth(labels, units = "inches", cex = cex) +
      abs(cos(srt * pi / 180)) * graphics::strheight(labels, units = "inches", cex = cex)
    bounds[[length(bounds) + 1L]] <<- list(labels = labels, bottom = anchor - height)
    draw_labels(at, labels, srt, cex, line_offset)
  }, .package = "mfrmr")
  for (size in list(c(7, 5), c(5, 4))) {
    grDevices::pdf(NULL, width = size[1], height = size[2])
    tryCatch({
      before <- suppressWarnings(plot(fit, type = "wright", draw = FALSE))
      drawn <- suppressWarnings(plot(fit, type = "wright"))
      expect_identical(drawn, before)
    }, finally = grDevices::dev.off())
  }
  expect_length(bounds, 2L)
  for (bound in bounds) {
    expect_true(any(nchar(bound$labels) > 16L))
    expect_false(anyDuplicated(bound$labels) > 0)
    expect_true(all(bound$bottom >= 0))
  }
})

test_that("bundle conflicts protect every member before publishing files", {
  figure <- new_mfrm_plot_data("test", list(title = "Saved title", note = "Saved note"))
  res <- structure(list(plot_map = data.frame(Type = "fit", Available = TRUE)), class = "mfrm_results")
  local_mocked_bindings(plot.mfrm_results = function(...) {
    graphics::plot(1:2); graphics::plot(2:1)
    structure(list(first = figure, second = figure), class = "mfrm_plot_bundle")
  }, .package = "mfrmr")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  existing <- file.path(folder, "mfrmr_results_plot_fit_second_notes.csv")
  writeLines("Keep this file", existing)
  device <- grDevices::dev.cur()
  expect_error(export_mfrm_results(res, folder, include = "plots", acknowledge_sensitive = TRUE),
    "File already exists")
  expect_identical(readLines(existing), "Keep this file")
  expect_identical(list.files(folder), basename(existing))
  expect_identical(grDevices::dev.cur(), device)
  directory <- file.path(folder, "mfrmr_results_plot_fit_first.png")
  dir.create(directory)
  expect_error(export_mfrm_results(res, folder, include = "plots", overwrite = TRUE,
    acknowledge_sensitive = TRUE), "Output path is a directory")
  expect_length(list.files(directory), 0L)
  expect_identical(readLines(existing), "Keep this file")
})

test_that("plot failures cannot silently lose pages or overwrite saved figures", {
  figure <- new_mfrm_plot_data("test", list())
  res <- structure(list(plot_map = data.frame(Type = c("fit", "pathway"), Available = TRUE)), class = "mfrm_results")
  mode <- "extra_page"
  local_mocked_bindings(plot.mfrm_results = function(x, type, ...) {
    graphics::plot(1:2)
    if (type == "pathway") return(figure)
    if (mode == "error") stop("Injected plot failure")
    if (mode == "extra_page") graphics::plot(2:1)
    figure
  }, .package = "mfrmr")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  existing <- file.path(folder, "mfrmr_results_plot_fit.png")
  writeLines("Keep existing image bytes", existing)
  device <- grDevices::dev.cur()
  for (current in c("extra_page", "error")) {
    mode <- current
    exported <- export_mfrm_results(res, folder, include = "plots", overwrite = TRUE,
      acknowledge_sensitive = TRUE)
    expect_identical(exported$plot_errors$Plot, "fit")
    expect_match(exported$plot_errors$Error, if (mode == "error") "Injected plot failure" else "page count")
    expect_identical(readLines(existing), "Keep existing image bytes")
    expect_false(file.exists(sub("\\.png$", "_notes.csv", existing)))
    expect_true(file.exists(file.path(folder, "mfrmr_results_plot_pathway.png")))
    expect_false(any(exported$written_files$Component == "plot_fit"))
    expect_identical(grDevices::dev.cur(), device)
  }
})

test_that("saved profile exports carry figure text outside the PNG without recomputation", {
  skip_if_not_installed("ggplot2")
  source <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))
  fail <- function(...) stop("Unexpected numerical computation")
  local_mocked_bindings(fit_mfrm = fail, mfrm_gpcm_profile_search = fail,
    mfrm_gpcm_product_inference = fail, .package = "mfrmr")
  res <- mfrm_results(source$fit, include = c("fit", "plots"), compute = "never",
    intervals = list(component = source$profiles$Task))
  path <- tempfile(); on.exit(unlink(path, recursive = TRUE), add = TRUE)
  exported <- export_mfrm_results(res, path, preset = "starter", acknowledge_sensitive = TRUE)
  expect_equal(nrow(exported$plot_errors), 0L)
  note_file <- exported$written_files$Path[exported$written_files$Component == "plot_gpcm_component_notes"]
  expect_length(note_file, 1L)
  notes <- read.csv(note_file, stringsAsFactors = FALSE)
  figure <- plot(res, type = "gpcm_component", draw = FALSE)
  expect_identical(notes$Text[notes$Type == "title"], plot_data(figure, "title"))
  expect_equal(notes[notes$Type != "title", ], plot_data(figure, "notes"), ignore_attr = TRUE)
  expect_true(all(attr(source$profiles$Task, "cautions") %in% notes$Text))
  expect_true(any(plot_data(figure, "table")$Display == "Unresolved"))
  manifest <- read.csv(file.path(path, "mfrmr_results_written_files.csv"))
  expect_true(note_file %in% manifest$Path)
  index <- paste(readLines(file.path(path, "index.html")), collapse = "\n")
  expect_match(index, basename(note_file), fixed = TRUE)
  expect_match(index, "Figure titles and notes", fixed = TRUE)
  restored <- readRDS(file.path(path, "mfrmr_results_results.rds"))
  again <- tempfile(); on.exit(unlink(again, recursive = TRUE), add = TRUE)
  replay <- export_mfrm_results(restored, again, include = "plots", acknowledge_sensitive = TRUE)
  expect_equal(nrow(replay$plot_errors), 0L)
  expect_identical(readLines(file.path(again, basename(note_file))), readLines(note_file))
})

test_that("base plot exports preserve source columns and protect existing note files", {
  notes <- data.frame(Type = "caution", Text = c("Shared caution", "Shared caution"),
    Fit = c("Reference", "Comparison"))
  figure <- mfrmr:::new_mfrm_plot_data("test", list(notes = notes))
  res <- structure(list(plot_map = data.frame(Type = "fit", Available = TRUE)), class = "mfrm_results")
  draws <- 0L
  local_mocked_bindings(plot.mfrm_results = function(...) {
    draws <<- draws + 1L
    graphics::plot(1:2)
    figure
  }, .package = "mfrmr")
  path <- tempfile(); on.exit(unlink(path, recursive = TRUE), add = TRUE)
  exported <- export_mfrm_results(res, path, include = "plots", acknowledge_sensitive = TRUE)
  expect_equal(draws, 1L)
  expect_equal(nrow(exported$plot_errors), 0L)
  note_file <- file.path(path, "mfrmr_results_plot_fit_notes.csv")
  expect_identical(read.csv(note_file, stringsAsFactors = FALSE), notes)
  unlink(file.path(path, "mfrmr_results_plot_fit.png"))
  before <- readLines(note_file)
  expect_error(export_mfrm_results(res, path, include = "plots", acknowledge_sensitive = TRUE),
    "File already exists")
  expect_identical(readLines(note_file), before)
  expect_equal(draws, 1L)
  figure <- mfrmr:::new_mfrm_plot_data("test", list())
  empty <- export_mfrm_results(res, path, include = "plots", overwrite = TRUE,
    acknowledge_sensitive = TRUE)
  expect_equal(nrow(empty$plot_errors), 0L)
  expect_identical(names(read.csv(note_file)), c("Type", "Text"))
  expect_equal(nrow(read.csv(note_file)), 0L)
})
