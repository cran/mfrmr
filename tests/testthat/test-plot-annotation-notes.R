test_that("figure text is returned without changing data or explicit display controls", {
  skip_if_not_installed("ggplot2")
  source <- mfrmr:::new_mfrm_plot_data("test", list(
    table = data.frame(x = 1:3, y = c(2, 1, 3)),
    title = "Suggested figure title", subtitle = "Method and level",
    caption = "Conditional uncertainty only.",
    notes = data.frame(Note = "Inspect unavailable values.")))
  p <- as_ggplot(source)
  expect_null(p$labels$title)
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  expect_identical(plot_data(p, "table"), source$data$table)
  expect_true(all(c("Method and level", "Conditional uncertainty only.",
    "Inspect unavailable values.") %in% plot_data(p, "notes")$Text))
  expect_match(paste(capture.output(print(source)), collapse = "\n"), "Conditional uncertainty only")
  source$data$display <- list(show_title = TRUE, show_notes = TRUE)
  annotated <- as_ggplot(source)
  expect_identical(annotated$labels$title, source$data$title)
  expect_identical(annotated$labels$caption, source$data$caption)
  expect_identical(plot_data(annotated, "notes"), plot_data(p, "notes"))
})

test_that("plot printing counts legend entries and names public reuse methods", {
  source <- new_mfrm_plot_data("test", list(
    notes = "Uncertainty is conditional on calibration.",
    legend = new_plot_legend(label = c("Estimate", "No interval"),
      role = c("point", "status"), aesthetic = "shape", value = c("16", "4"))))
  printed <- paste(capture.output(returned <- print(source)), collapse = "\n")
  expect_identical(returned, source)
  expect_match(printed, "legend   : 2 entries", fixed = TRUE)
  expect_match(printed, "Uncertainty is conditional on calibration", fixed = TRUE)
  expect_match(printed, "plot_data()", fixed = TRUE)
  expect_match(printed, "For supported views, `as_ggplot()`", fixed = TRUE)
  expect_false(grepl("originating|inverse|plotly", printed))
})

test_that("score-direction notes use fitted signs without guessing rubric meaning", {
  fit <- structure(list(config = list(facet_names = c("Judge", "Task", "Unknown"),
    facet_signs = c(Judge = -1, Task = 1))), class = "mfrm_fit")
  notes <- .mfrm_score_direction_notes(fit)
  text <- paste(notes$Text, collapse = "\n")
  expect_match(text, "higher coded scores need not mean better", fixed = TRUE)
  expect_match(text, "Judge: higher locations reduce expected coded scores", fixed = TRUE)
  expect_match(text, "Task: higher locations increase expected coded scores", fixed = TRUE)
  expect_match(text, "Unknown: higher locations have no recorded score direction", fixed = TRUE)
  expect_false(any(grepl("Person locations", .mfrm_score_direction_notes(
    fit, facets = "Judge", include_person = FALSE)$Text)))
  class(fit) <- c("mfrm_imported_fit", class(fit))
  text <- paste(.mfrm_score_direction_notes(fit)$Text, collapse = "\n")
  expect_match(text, "source model's orientation", fixed = TRUE)
  expect_false(grepl("Judge:|Higher Person", text))
})

test_that("annotation merging ignores outer whitespace without changing saved prose", {
  stored <- data.frame(Type = "caption", Text = "Points are requested scenarios. ",
    Origin = "saved", stringsAsFactors = FALSE)
  payload <- list(notes = stored, caption = "Points are requested scenarios.",
    cautions = c(" Review uncertainty. ", "Review uncertainty.", " ", NA_character_))
  result <- mfrmr:::.mfrm_plot_annotation_notes(payload)
  expect_identical(result$Text, c(stored$Text, " Review uncertainty. "))
  expect_identical(result$Type, c("caption", "caution"))
  expect_identical(result$Origin, c("saved", NA_character_))
  expect_identical(payload$notes, stored)
  payload$notes <- result
  expect_identical(mfrmr:::.mfrm_plot_annotation_notes(payload), result)
})

test_that("base reliability graphics keep the interpretation in returned notes", {
  fit <- structure(list(), class = "mfrm_fit")
  diagnostics <- list(reliability = data.frame(Facet = c("Person", "Rater"),
    Reliability = c(.8, .9)))
  text <- character()
  original_title <- graphics::title
  local_mocked_bindings(
    title = function(main = NULL, sub = NULL, ...) {
      text <<- c(text, as.character(main), as.character(sub))
      original_title(main = main, sub = sub, ...)
    },
    mtext = function(text, ...) {
      stop("No margin notes are expected for this figure")
    }, .package = "mfrmr")
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  local_graphics_headings(function(main) text <<- c(text, as.character(main)))
  p <- plot_reliability_snapshot(fit, diagnostics = diagnostics)
  expect_false(any(nzchar(text)))
  notes <- plot_data(p, "notes")
  expect_true(any(grepl("severity differences, not agreement", notes$Text, fixed = TRUE)))
  expect_true(any(grepl("Low rater values can also reflect imprecision", notes$Text, fixed = TRUE)))
  expect_identical(p$data$table$Value, c(.9, .8))
  expect_match(paste(capture.output(p), collapse = "\n"), "not agreement")
})

test_that("saved profile plots retain cautions and unresolved points outside figures", {
  skip_if_not_installed("ggplot2")
  source <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))
  ci <- source$profiles$Task
  local_mocked_bindings(
    fit_mfrm = function(...) stop("No fitting"),
    mfrm_gpcm_profile_search = function(...) stop("No profile search"),
    mfrm_gpcm_product_inference = function(...) stop("No information calculation"),
    .package = "mfrmr")
  p <- plot(ci, draw = FALSE)
  expect_null(p$labels$title)
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  expect_match(p$labels$x, "discrimination")
  expect_identical(p$labels$shape, "Numerical check")
  expect_true(any(plot_data(p, "table")$Display == "Unresolved"))
  expect_identical(plot_data(p, "endpoints"), attr(ci, "profile")$endpoints)
  expect_true(any(grepl("Dotted: Wald limits", plot_data(p, "notes")$Text, fixed = TRUE)))
  expect_true(all(attr(ci, "cautions") %in% plot_data(p, "notes")$Text))
  custom <- plot(ci, title = NULL, subtitle = "Custom description", caption = NULL, draw = FALSE)
  expect_true(all(attr(ci, "cautions") %in% plot_data(custom, "notes")$Text))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(p, path)
  restored <- readRDS(path)
  expect_identical(plot_data(restored, "notes"), plot_data(p, "notes"))
  expect_identical(ggplot2::ggplot_build(restored)$data, ggplot2::ggplot_build(p)$data)
})

test_that("shared bars distinguish structural panel labels from figure headings", {
  seen <- character()
  local_graphics_headings(function(main) seen <<- c(seen, as.character(main)))
  original <- graphics::title
  local_mocked_bindings(title = function(main = NULL, ...) {
    seen <<- c(seen, as.character(main))
    original(main = main, ...)
  }, .package = "mfrmr")
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  mfrmr:::barplot_rot45(1:2, c("A", "B"), "grey", main = "Figure title")
  expect_false(any(nzchar(seen)))
  mfrmr:::barplot_rot45(1:2, c("A", "B"), "grey", main = "Panel A", show_title = TRUE)
  expect_true("Panel A" %in% seen)
})
