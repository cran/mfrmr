test_that("appendix export retains analysis tables with metadata names", {
  tables <- list(
    overview = data.frame(Persons = 3L, NotScored = 1L),
    reporting_map = data.frame(Source = "Stored scoring results"),
    table_overview = data.frame(Other = 7L),
    estimates = data.frame(Person = c("001", "NA"), Estimate = c(.25, NA_real_)),
    table_index = data.frame(AnalysisIndex = "Retain this source table")
  )
  roles <- setNames(rep("unit_estimates", length(tables)), names(tables))
  bundle <- structure(list(
    overview = data.frame(Title = "Saved scores", SourceClass = "fixture",
      SummaryClass = "summary.fixture", TablesAvailable = length(tables),
      TablesReturned = length(tables), AppendixPreset = "none"),
    table_index = build_summary_table_index(tables, roles, roles),
    plot_index = summary_table_bundle_plot_index(tables),
    tables = tables, notes = "Conditional on the saved calibration.",
    source_class = "fixture", summary_class = "summary.fixture",
    appendix_preset = "none"
  ), class = "mfrm_summary_table_bundle")
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    score_mfrm_calibration = function(...) stop("No scoring"))
  path <- tempfile("appendix-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  result <- export_summary_appendix(list(saved = bundle), path, include_html = TRUE)
  files <- result$written_files
  expect_false(anyDuplicated(files$Path) > 0L)
  tags <- c("table_overview", "table_reporting_map", "table_overview_1",
    "estimates", "table_table_index")
  for(i in seq_along(tables)) {
    f <- files$Path[files$Component == paste0("summary_saved_", tags[[i]]) & files$Format == "csv"]
    expect_length(f, 1L)
    expected <- tables[[i]]
    actual <- read.csv(f, colClasses = vapply(expected, function(x) class(x)[1], ""),
      na.strings = "", check.names = FALSE)
    expect_identical(actual, expected)
  }
  metadata <- read.csv(files$Path[files$Component == "summary_saved_overview"],
    check.names = FALSE)
  expect_identical(metadata$Title, "Saved scores")
  html <- paste(readLines(files$Path[files$Format == "html"], warn = FALSE), collapse = "\n")
  expect_match(html, "NotScored", fixed = TRUE)
  expect_match(html, "Stored scoring results", fixed = TRUE)
  expect_match(html, "Retain this source table", fixed = TRUE)
  expect_error(export_summary_appendix(list(saved = bundle), path), "already exists")
})

test_that("appendix summaries count only the requested preset", {
  catalog <- data.frame(
    Bundle = c("A", "A", "B"), Table = c("estimates", "fit", "estimates"),
    Rows = c(2L, 3L, 4L), Role = c("estimates", "diagnostics", "estimates"),
    AppendixSection = c("results", "diagnostics", "results"),
    PreferredAppendixOrder = c(1L, 2L, 1L),
    PlotReady = c(TRUE, FALSE, FALSE), NumericColumns = c(1L, 0L, 2L),
    ExportReady = TRUE, ApaTableReady = TRUE
  )
  all <- transform(catalog, Preset = "all", Selected = TRUE)
  compact <- transform(catalog, Preset = "compact", Selected = c(TRUE, FALSE, FALSE))
  empty <- transform(catalog, Preset = "reporting", Selected = FALSE)
  mixed <- rbind(all, compact, empty)
  summaries <- list(
    export_summary_table_selection_summary,
    export_summary_table_selection_table_summary,
    export_summary_table_selection_handoff_summary,
    export_summary_table_selection_handoff_bundle_summary,
    export_summary_table_selection_handoff_preset_summary,
    export_summary_table_selection_handoff_role_summary,
    export_summary_table_selection_handoff_role_section_summary,
    export_summary_table_selection_role_summary,
    export_summary_table_selection_section_summary
  )
  for (fun in summaries) {
    for (preset in c("all", "compact", "reporting")) {
      selected_input <- mixed[mixed$Preset == preset, , drop = FALSE]
      expect_identical(fun(mixed, preset), fun(selected_input, preset))
    }
    expect_equal(nrow(fun(mixed, "methods")), 0L)
    expect_identical(fun(NULL, "compact"), data.frame())
  }
  counts <- export_summary_table_selection_summary(mixed, "compact")
  expect_identical(counts$TablesAvailable, c(2L, 1L))
  expect_identical(counts$TablesSelected, c(1L, 0L))
  expect_equal(counts$SelectionFraction, c(.5, 0))
  expect_equal(counts$PlotReadyFraction, c(1, NA_real_))
  tables <- export_summary_table_selection_table_summary(mixed, "all")
  expect_identical(tables$Table, c("estimates", "fit"))
  expect_identical(tables$PresetsSelected, c(1L, 1L))
  roles <- export_summary_table_selection_role_summary(mixed, "compact")
  expect_identical(roles$Role, "estimates")
  expect_identical(roles$Tables, 1L)
  expect_identical(roles$BundlesCovered, 1L)
  handoff <- export_summary_table_selection_handoff_preset_summary(mixed, "compact")
  expect_identical(handoff$Tables, 1L)
  expect_equal(handoff$NumericFraction, 1)
})
