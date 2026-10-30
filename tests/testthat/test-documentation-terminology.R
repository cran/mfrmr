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

test_that("the public guide starts from data, fit, plot, and summary", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  readme <- paste(readLines(file.path(pkg_root, "README.md"), warn = FALSE),
                  collapse = "\n")

  data_pos <- regexpr("ratings <- load_mfrmr_data", readme, fixed = TRUE)[1]
  fit_pos <- regexpr("fit <- fit_mfrm", readme, fixed = TRUE)[1]
  plot_pos <- regexpr("plot(fit)", readme, fixed = TRUE)[1]
  summary_pos <- regexpr("results <- summary(fit)", readme, fixed = TRUE)[1]

  expect_true(all(c(data_pos, fit_pos, plot_pos, summary_pos) > 0L))
  expect_lt(data_pos, fit_pos)
  expect_lt(fit_pos, plot_pos)
  expect_lt(plot_pos, summary_pos)
  expect_match(readme, "method = \"MML\"", fixed = TRUE)
})

test_that("FACETS positioning distinguishes native estimates from visual emulation", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  readme <- paste(readLines(file.path(pkg_root, "README.md"), warn = FALSE),
                  collapse = "\n")

  expect_match(readme, "mfrmr-facets-migration.html", fixed = TRUE)
  guide <- paste(readLines(file.path(pkg_root, "vignettes", "mfrmr-facets-migration.Rmd"),
    warn = FALSE), collapse = "\n")
  text <- gsub("[[:space:]]+", " ", paste(readme, guide))
  expect_match(guide, "mfrmr uncertainty estimates", fixed = TRUE)
  expect_match(guide, "FACETS-style", fixed = TRUE)
  expect_true(grepl("asterisk", tolower(guide), fixed = TRUE))
  expect_true(grepl("rubric", tolower(guide), fixed = TRUE))
  expect_match(text, paste("Neither renderer implies that FACETS performed the estimation",
    "or that the two programs are numerically equivalent."), fixed = TRUE)
  expect_false(grepl("FACETS produced the estimates", readme, fixed = TRUE))
  expect_false(grepl("complete facets reproduction", tolower(readme),
                     fixed = TRUE))
})

test_that("public documentation does not advertise later-release routes as current", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  docs <- read_public_text(pkg_root, public_documentation_files(pkg_root))
  text <- tolower(paste(unlist(docs, use.names = FALSE), collapse = "\n"))
  flat <- gsub("[[:space:]]+", " ", text)

  prohibited_claims <- c(
    "supports unrestricted gpcm",
    "supports an unrestricted gpcm",
    "provides unrestricted gpcm",
    "supports multiple observed scales",
    "supports multiple independent rating scales",
    "supports general threshold anchors",
    "supports threshold anchoring",
    "supports native multidimensional estimation",
    "provides native multidimensional estimation"
  )
  hits <- prohibited_claims[vapply(
    prohibited_claims,
    grepl,
    logical(1),
    x = flat,
    fixed = TRUE
  )]

  expect_identical(hits, character(0))
  expect_match(flat, "one observed score scale", fixed = TRUE)
  portable <- tolower(gsub("[[:space:]]+", " ", paste(readLines(file.path(
    pkg_root, "vignettes", "mfrmr-portable-calibration.Rmd"), warn = FALSE), collapse = "\n")))
  expect_match(portable, "portable workflow is deliberately narrower than fitted-object scoring", fixed = TRUE)
  expect_true(grepl(
    "estimated-population rsm/pcm and latent regression remain fitted-object-only",
    portable, fixed = TRUE
  ))
  expect_match(portable, "portable calibration are different objects", fixed = TRUE)
})

test_that("public estimator vocabulary keeps JMLE as an input alias only", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  docs <- read_public_text(pkg_root, public_documentation_files(pkg_root))
  public_text <- paste(unlist(docs, use.names = FALSE), collapse = "\n")

  fit_help <- paste(readLines(
    file.path(pkg_root, "man", "fit_mfrm.Rd"),
    warn = FALSE,
    encoding = "UTF-8"
  ), collapse = "\n")
  package_help <- paste(readLines(
    file.path(pkg_root, "man", "mfrmr-package.Rd"),
    warn = FALSE,
    encoding = "UTF-8"
  ), collapse = "\n")

  expect_match(fit_help, 'method = c("MML", "JML", "JMLE", "CML")', fixed = TRUE)
  expect_match(fit_help, '"JMLE"} is accepted as a', fixed = TRUE)
  expect_match(package_help, '"JMLE"} as a backward-compatible alias', fixed = TRUE)
  expect_false(grepl('method = "MMLE"', public_text, fixed = TRUE))
  expect_false(grepl("\\bMMLE\\b", public_text, perl = TRUE))
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

test_that("first-use guides share one decision-first readiness route", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  paths <- c(
    file.path(pkg_root, "README.md"),
    file.path(pkg_root, "vignettes", "mfrmr-workflow.Rmd"),
    file.path(pkg_root, "vignettes", "mfrmr-gpcm-scope.Rmd"),
    file.path(pkg_root, "vignettes", "mfrmr-facets-migration.Rmd")
  )
  docs <- read_public_text(pkg_root, paths)

  for (path in names(docs)) {
    text <- paste(docs[[path]], collapse = "\n")
    expect_match(text, "$decision", fixed = TRUE, info = path)
    expect_match(text, "FormalInference", fixed = TRUE, info = path)
  }
})

test_that("CRAN-facing documentation excludes development-process language", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  docs <- read_public_text(pkg_root, public_documentation_files(pkg_root))
  expect_gt(length(docs), 0L)

  blocked <- c(
    "release evidence map",
    "release-readiness",
    "repository roadmap",
    "repository-level",
    "validation record",
    "CRAN checks",
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
    "repository theorem-only",
    "internal analytic prototype",
    "retrospective 40-dataset owner panel",
    "predeclared deterministic challenge",
    "roadmap_only",
    "schema-only future branch",
    "future-branch active planning scaffold",
    "regression-tested default",
    "focused regression test",
    "test code mechanics",
    "0.2.3 evidence map",
    "repository audit",
    "candidate-linked evidence contract",
    "optional local execution",
    "for final evidence",
    "executable examples shorten package checks",
    "keeps CRAN example time short",
    "reproducible analysis, package checks",
    "documents, tests, and reproducible analyses"
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
    idx <- grep("\\bP[0-9]+[a-z]\\b", lines, perl = TRUE)
    if (length(idx) > 0L) {
      hits <- c(hits, paste0(path, ":", idx, ": ", trimws(lines[idx])))
    }
    idx <- grep(
      "\\b(lane|lanes|gate|gates|preflight)\\b",
      lines, ignore.case = TRUE, perl = TRUE
    )
    if (length(idx) > 0L) {
      hits <- c(hits, paste0(path, ":", idx, ": ", trimws(lines[idx])))
    }
    idx <- grep(
      "(CORE-[0-9]+|OPT-[0-9]+|Draft[.][0-9]+|(?<![[:alnum:]-])G[0-9](?![[:alnum:]]))",
      lines, perl = TRUE
    )
    if (length(idx) > 0L) {
      hits <- c(hits, paste0(path, ":", idx, ": ", trimws(lines[idx])))
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
  ignore_path <- file.path(pkg_root, ".Rbuildignore")
  if (file.exists(ignore_path)) {
    ignore <- readLines(ignore_path, warn = FALSE)
    ignore <- ignore[nzchar(ignore) & !startsWith(trimws(ignore), "#")]
    relative <- substring(
      normalizePath(source_files, winslash = "/", mustWork = FALSE),
      nchar(normalizePath(pkg_root, winslash = "/", mustWork = TRUE)) + 2L
    )
    ignored <- vapply(relative, function(path) {
      any(vapply(ignore, grepl, logical(1), x = path, perl = TRUE))
    }, logical(1))
    source_files <- source_files[!ignored]
  }
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
    "remaining-coverage",
    "P1z compactification",
    "P2b lexicographic-limit",
    "P2c parameter",
    "P2d finite-sequence",
    "P2e contract",
    "P2e oracle",
    "P2e flag",
    "P2f global-existence",
    "P2g source contract",
    "P2g scope",
    "P2h binary",
    "P2i fixture",
    "P2i envelope",
    "P2j face chart"
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
