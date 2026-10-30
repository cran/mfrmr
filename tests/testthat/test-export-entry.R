export_entry_fixture <- local({
  fit <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))$RSM$fit
  list(fit = fit, diagnostics = suppressWarnings(diagnose_mfrm(fit, residual_pca = "none")),
    result = mfrm_results(fit, include = "fit", compute = "never"))
})

test_that("the CSV-only export has one descriptive public name", {
  expect_true("export_mfrm_tables" %in% getNamespaceExports("mfrmr"))
  expect_false("export_mfrm" %in% getNamespaceExports("mfrmr"))
  expect_false(exists("export_mfrm", envir = asNamespace("mfrmr"), inherits = FALSE))
  directory <- tempfile("export-tables-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  x <- export_entry_fixture
  original <- serialize(x, NULL)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"), .package = "mfrmr")
  written <- export_mfrm_tables(x$fit, diagnostics = x$diagnostics,
    output_dir = directory, acknowledge_sensitive = TRUE)
  expect_equal(sort(written$Table), sort(c("person", "facets", "summary", "steps", "measures")))
  expect_true(all(file.exists(written$Path)))
  expect_equal(read.csv(written$Path[written$Table == "person"])$Estimate,
    x$fit$facets$person$Estimate, tolerance = 1e-12)
  expect_identical(serialize(x, NULL), original)
})

test_that("fit bundles list and zip only core tables written by this call", {
  skip_if(!nzchar(Sys.which(Sys.getenv("R_ZIPCMD", "zip"))), "ZIP command unavailable")
  x <- export_entry_fixture
  directory <- tempfile("export-core-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  stale <- file.path(directory, "probe_measures.csv")
  writeLines("retained older measures", stale)
  original <- tools::md5sum(stale)
  diagnostics <- x$diagnostics
  diagnostics$measures <- NULL
  out <- export_mfrm_bundle(x$fit, diagnostics = diagnostics, output_dir = directory,
    prefix = "probe", include = "core_tables", zip_bundle = TRUE, acknowledge_sensitive = TRUE)
  expect_false("core_measures" %in% out$written_files$Component)
  expect_identical(tools::md5sum(stale), original)
  archive <- out$written_files$Path[out$written_files$Format == "zip"]
  expect_length(archive, 1L)
  expect_false(basename(stale) %in% utils::unzip(archive, list = TRUE)$Name)
  expect_equal(out$summary$FilesWritten, nrow(out$written_files))
})

test_that("facet CSVs retain qualification fields and source row order", {
  x <- export_entry_fixture
  directory <- tempfile("export-qualification-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  # Diagnostic row order is independent of the fit's facet-table order.
  diagnostics <- x$diagnostics
  diagnostics$measures <- diagnostics$measures[nrow(diagnostics$measures):1, ]
  out <- export_mfrm_tables(x$fit, diagnostics = diagnostics, output_dir = directory,
    tables = "facets", acknowledge_sensitive = TRUE)
  table <- read.csv(out$Path, na.strings = "", check.names = FALSE)
  key <- function(x) paste(x$Facet, x$Level, sep = "\r")
  expect_identical(key(table), key(x$fit$facets$others))
  idx <- match(key(table), key(diagnostics$measures))
  for (column in c("Estimate", "SE", "Infit", "Outfit", "SE_Method", "SEUse",
    "CIUse", "InferenceReady", "SupportsFormalInference", "FitStatus", "FitReason")) {
    expect_equal(table[[column]], unname(diagnostics$measures[[column]][idx]), tolerance = 1e-12,
      info = column)
  }
  expect_true(any(!table$SupportsFormalInference))
  expect_true(all(nzchar(table$FitReason)))
  expect_false(any(grepl("[.][xy]$", names(table))))
})

test_that("all archive routes replace ZIP members and report failed commands", {
  skip_if(!nzchar(Sys.which(Sys.getenv("R_ZIPCMD", "zip"))), "ZIP command unavailable")
  x <- export_entry_fixture
  routes <- list(
    results = function(path, first) export_mfrm_results(x$result, path,
      include = if (first) "tables" else "rds", zip_bundle = TRUE,
      overwrite = TRUE, acknowledge_sensitive = TRUE),
    bundle = function(path, first) export_mfrm_bundle(x$fit, diagnostics = x$diagnostics,
      output_dir = path, include = if (first) c("core_tables", "html") else "core_tables",
      zip_bundle = TRUE, overwrite = TRUE, acknowledge_sensitive = TRUE),
    appendix = function(path, first) export_summary_appendix(x$fit,
      output_dir = path, include_html = first, zip_bundle = TRUE, overwrite = TRUE)
  )
  directory <- tempfile("export-zip-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  for (name in names(routes)) {
    path <- file.path(directory, name)
    first <- routes[[name]](path, TRUE)
    second <- routes[[name]](path, FALSE)
    archive <- second$written_files$Path[second$written_files$Format == "zip"]
    expect_length(archive, 1L)
    expected <- basename(second$written_files$Path[second$written_files$Format != "zip"])
    expect_equal(sort(utils::unzip(archive, list = TRUE)$Name), sort(expected), info = name)
    removed <- setdiff(first$written_files$Path, second$written_files$Path)
    expect_gt(length(removed), 0L)
    expect_true(all(file.exists(removed)), info = "Unrequested loose files are preserved")
    old_zip <- tools::md5sum(archive)
    failed <- with_mocked_bindings(routes[[name]](path, FALSE),
      zip = function(...) 12L, .package = "mfrmr")
    expect_false(any(failed$written_files$Format == "zip"), info = name)
    expect_false(as.logical(failed$summary$ZipWritten), info = name)
    expect_match(paste(c(failed$zip_note, failed$notes), collapse = " "), "exit status: 12")
    expect_identical(tools::md5sum(archive), old_zip)
    expect_length(list.files(path, pattern = "^[.]mfrmr-export-", all.files = TRUE), 0L)
  }
})

test_that("ZIP validation preserves old archives on missing or incomplete output", {
  skip_if(!nzchar(Sys.which(Sys.getenv("R_ZIPCMD", "zip"))), "ZIP command unavailable")
  directory <- tempfile("zip-validation-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- file.path(directory, "current.txt"); writeLines("current", input)
  path <- file.path(directory, "archive.zip")
  expect_null(mfrm_write_export_zip(path, input)$note)
  before <- tools::md5sum(path)
  original_zip <- utils::zip
  other <- file.path(directory, "other.txt"); writeLines("old", other)
  scenarios <- list(
    missing = function(...) 0L,
    mismatch = function(zipfile, files, ...) original_zip(zipfile, other, ...),
    invalid = function(zipfile, ...) { writeLines("not a zip", zipfile); 0L },
    error = function(...) stop("Injected ZIP error")
  )
  for (name in names(scenarios)) {
    out <- with_mocked_bindings(mfrm_write_export_zip(path, input, overwrite = TRUE),
      zip = scenarios[[name]], .package = "mfrmr")
    expect_null(out$path)
    expect_true(nzchar(out$note), info = name)
    expect_identical(tools::md5sum(path), before)
  }
  expect_error(mfrm_write_export_zip(path, input), "already exists")
  expect_error(mfrm_write_export_zip(input, input, overwrite = TRUE), "must differ")
  expect_error(mfrm_write_export_zip(directory, input), "is a directory")
})

test_that("on-disk results manifests include plot errors and consistent counts", {
  x <- export_entry_fixture$result
  x$plot_map <- data.frame(Type = "missing_route", Available = TRUE)
  directory <- tempfile("export-manifest-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  out <- export_mfrm_results(x, directory, include = c("plots", "manifest"),
    acknowledge_sensitive = TRUE)
  expect_equal(nrow(out$plot_errors), 1L)
  manifest <- read.csv(file.path(directory, "mfrmr_results_written_files.csv"))
  overview <- read.csv(file.path(directory, "mfrmr_results_export_summary.csv"))
  expect_true("plot_errors" %in% manifest$Component)
  expect_identical(manifest$Component, out$written_files$Component)
  expect_equal(overview$FilesWritten, nrow(manifest))
  expect_equal(overview$CsvWritten, sum(manifest$Format == "csv"))
  expect_equal(overview$FilesWritten, out$summary$FilesWritten)
})
