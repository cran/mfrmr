documentation_source_root <- function() {
  test_root <- normalizePath(testthat::test_path(), mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    file.path(test_root, "..", ".."),
    getwd(),
    file.path(getwd(), "..")
  ), mustWork = FALSE))

  matches <- candidates[
    file.exists(file.path(candidates, "DESCRIPTION")) &
      file.exists(file.path(candidates, "README.md")) &
      file.exists(file.path(candidates, "R", "api-estimation.R")) &
      dir.exists(file.path(candidates, "man"))
  ]
  if (length(matches) == 0L) NA_character_ else matches[1]
}

read_public_text <- function(pkg_root, paths) {
  paths <- paths[file.exists(paths)]
  stats::setNames(
    lapply(paths, readLines, warn = FALSE, encoding = "UTF-8"),
    sub(paste0("^", pkg_root, "/?"), "", paths)
  )
}

public_documentation_files <- function(pkg_root) {
  unique(c(
    file.path(pkg_root, "DESCRIPTION"),
    file.path(pkg_root, "README.md"),
    file.path(pkg_root, "NEWS.md"),
    list.files(file.path(pkg_root, "vignettes"), pattern = "[.]Rmd$",
               full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "inst", "references"),
               pattern = "[.](md|txt|csv)$", full.names = TRUE,
               recursive = TRUE),
    list.files(file.path(pkg_root, "inst", "extdata"),
               pattern = "[.](md|txt)$", full.names = TRUE,
               recursive = TRUE),
    list.files(file.path(pkg_root, "man"), pattern = "[.]Rd$",
               full.names = TRUE, recursive = TRUE)
  ))
}

test_that("the public guide starts from data, fit, summary, and required Wright output", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  readme <- paste(readLines(file.path(pkg_root, "README.md"), warn = FALSE),
                  collapse = "\n")

  data_pos <- regexpr("dat <- load_mfrmr_data", readme, fixed = TRUE)[1]
  fit_pos <- regexpr("fit <- fit_mfrm", readme, fixed = TRUE)[1]
  summary_pos <- regexpr("fit_summary <- summary", readme, fixed = TRUE)[1]
  wright_pos <- regexpr("required native wright map", tolower(readme),
                        fixed = TRUE)[1]

  expect_true(all(c(data_pos, fit_pos, summary_pos, wright_pos) > 0L))
  expect_lt(data_pos, fit_pos)
  expect_lt(fit_pos, summary_pos)
  expect_lt(summary_pos, wright_pos)
  expect_match(readme, "method = \"MML\"", fixed = TRUE)
})

test_that("FACETS positioning distinguishes native estimates from visual emulation", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  readme <- paste(readLines(file.path(pkg_root, "README.md"), warn = FALSE),
                  collapse = "\n")

  readme_flat <- gsub("[[:space:]]+", " ", readme)
  expect_match(readme_flat, "keeps uncertainty visible", fixed = TRUE)
  expect_match(readme, "FACETS-style", fixed = TRUE)
  expect_true(grepl("asterisk", tolower(readme), fixed = TRUE))
  expect_true(grepl("rubric", tolower(readme), fixed = TRUE))
  expect_true(grepl("does not reproduce a facets run", tolower(readme),
                    fixed = TRUE))
  expect_false(grepl("FACETS produced the estimates", readme, fixed = TRUE))
  expect_false(grepl("complete facets reproduction", tolower(readme),
                     fixed = TRUE))
})

test_that("ConQuest guidance states the external MML comparison boundary", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  readme <- paste(readLines(file.path(pkg_root, "README.md"), warn = FALSE),
                  collapse = "\n")
  guide_source <- paste(readLines(
    file.path(pkg_root, "R", "help_reports_and_tables.R"), warn = FALSE
  ), collapse = "\n")
  combined <- paste(readme, guide_source, sep = "\n")

  expect_match(combined, "ConQuest", fixed = TRUE)
  expect_match(combined, "MML", fixed = TRUE)
  expect_true(grepl("does not execute conquest", tolower(combined), fixed = TRUE))
  expect_true(grepl("parse raw conquest output", tolower(combined), fixed = TRUE))
})

test_that("public maxit guidance prevents result-driven tuning", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  paths <- c(
    file.path(pkg_root, "README.md"),
    file.path(pkg_root, "man", "fit_mfrm.Rd"),
    file.path(pkg_root, "man", "mfrmr_workflow_methods.Rd"),
    file.path(pkg_root, "vignettes", "mfrmr-workflow.Rmd"),
    file.path(pkg_root, "vignettes", "mfrmr-mml-and-marginal-fit.Rmd")
  )
  docs <- paste(unlist(read_public_text(pkg_root, paths), use.names = FALSE),
                collapse = "\n")
  docs_flat <- gsub("[[:space:]]+", " ", docs)

  expect_match(docs_flat, "computational ceiling", fixed = TRUE)
  expect_match(docs_flat, "not a convergence criterion", fixed = TRUE)
  expect_match(docs_flat, "prespecified", fixed = TRUE)
  expect_match(docs_flat, "same data, model, method", fixed = TRUE)
  expect_match(docs_flat, "Do not select", fixed = TRUE)
  expect_match(docs_flat, "ConvergenceStatus", fixed = TRUE)
  expect_match(docs_flat, "InferenceReady", fixed = TRUE)
  expect_match(docs_flat, "Numerical", fixed = TRUE)
  expect_match(docs_flat, "Do not select it merely as a faster substitute", fixed = TRUE)
  expect_false(grepl("fast exploratory JML pass", docs_flat, fixed = TRUE))
})

test_that("CRAN-facing documentation excludes development-process language", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  docs <- read_public_text(pkg_root, public_documentation_files(pkg_root))
  expect_gt(length(docs), 0L)

  blocked <- c(
    "release evidence map",
    "release-readiness",
    "pre-release status",
    "checked locally",
    "regression fixture",
    "test fixture",
    "golden output",
    "golden comparison",
    "fast smoke run",
    "smoke check",
    "current public branch",
    "repository-only",
    "for efficient development",
    "public first-contact route",
    "review every artifact",
    "lower-level component",
    "full per-helper support contract",
    "main package architecture",
    "blocked explicitly",
    "binding contract",
    "manuscript surface",
    "structured drafting and review surface",
    "reportable surface",
    "fitted-scale artifact",
    "planner-schema contract",
    "second-wave visual layer",
    "package-test coverage",
    "roadmap_only",
    "schema-only future branch",
    "future-branch active planning scaffold"
  )

  hits <- character(0)
  for (path in names(docs)) {
    lines <- docs[[path]]
    for (phrase in blocked) {
      idx <- grep(tolower(phrase), tolower(lines), fixed = TRUE)
      if (length(idx) > 0L) {
        hits <- c(hits, paste0(path, ":", idx, ": ", trimws(lines[idx])))
      }
    }
  }
  expect_identical(hits, character(0))
})

test_that("internal validation records are excluded from the source package", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  ignore_path <- file.path(pkg_root, ".Rbuildignore")
  testthat::skip_if_not(
    file.exists(ignore_path),
    ".Rbuildignore is available only in a source checkout"
  )
  ignore <- readLines(ignore_path, warn = FALSE)
  expect_true("^inst/validation$" %in% ignore)
})

test_that("the shipped source excludes development-stage contract names", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source package files are not available")

  source_files <- unique(c(
    file.path(pkg_root, "NEWS.md"),
    file.path(pkg_root, "NAMESPACE"),
    list.files(file.path(pkg_root, "R"), pattern = "[.]R$",
               full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "tests", "testthat"), pattern = "[.]R$",
               full.names = TRUE, recursive = TRUE)
  ))
  source_files <- source_files[
    basename(source_files) != "test-documentation-terminology.R"
  ]
  bad_file_names <- grep(
    "coverage-(push|boost|gaps)|final-coverage|remaining-coverage",
    basename(source_files),
    value = TRUE,
    ignore.case = TRUE
  )
  expect_identical(bad_file_names, character(0))
  source_text <- read_public_text(pkg_root, source_files)

  blocked <- c(
    "future_branch",
    "future branch",
    "future-branch",
    "schema_only",
    "schema-only",
    "planning scaffold",
    "roadmap_only",
    "created with 0.2.2",
    "coverage-push",
    "coverage-boost",
    "coverage-gaps",
    "final-coverage",
    "remaining-coverage"
  )
  hits <- character(0)
  for (path in names(source_text)) {
    lines <- source_text[[path]]
    for (phrase in blocked) {
      idx <- grep(tolower(phrase), tolower(lines), fixed = TRUE)
      if (length(idx) > 0L) {
        hits <- c(hits, paste0(path, ":", idx, ": ", trimws(lines[idx])))
      }
    }
  }
  expect_identical(hits, character(0))
})

test_that("high-level help pages retain searchable user concepts", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  expected <- list(
    "mfrmr_output_guide.Rd" = c("reporting workflow", "route selection"),
    "gpcm_capability_matrix.Rd" = c("GPCM boundaries", "route selection"),
    "mfrmr_visual_diagnostics.Rd" = c("visual diagnostics"),
    "plot_wright_unified.Rd" = c("Wright maps")
  )
  missing <- character(0)
  for (rd in names(expected)) {
    path <- file.path(pkg_root, "man", rd)
    expect_true(file.exists(path), info = rd)
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    for (concept in expected[[rd]]) {
      if (!grepl(paste0("\\concept{", concept, "}"), text, fixed = TRUE)) {
        missing <- c(missing, paste(rd, concept, sep = ": "))
      }
    }
  }
  expect_identical(missing, character(0))
})
