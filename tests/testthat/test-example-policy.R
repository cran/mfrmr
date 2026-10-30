example_policy_source_root <- function() {
  test_root <- normalizePath(testthat::test_path(), winslash = "/", mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    file.path(test_root, "..", ".."),
    file.path(test_root, "..", "..", "00_pkg_src", "mfrmr"),
    file.path(test_root, "..", "..", "..", "00_pkg_src", "mfrmr"),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "00_pkg_src", "mfrmr"),
    file.path(getwd(), "..", "..", "00_pkg_src", "mfrmr")
  ), winslash = "/", mustWork = FALSE))

  candidates[file.exists(file.path(candidates, "DESCRIPTION")) &
               file.exists(file.path(candidates, "R", "api-estimation.R")) &
               file.exists(file.path(candidates, "tests", "testthat.R"))][1]
}

example_policy_roxygen_examples <- function(pkg_root) {
  files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                      recursive = TRUE, full.names = TRUE)
  rows <- list()
  for (path in files) {
    lines <- readLines(path, warn = FALSE)
    starts <- grep("^#' @examples", lines)
    for (start in starts) {
      block <- character()
      i <- start + 1L
      while (i <= length(lines) && grepl("^#'", lines[[i]])) {
        if (grepl("^#' @", lines[[i]])) break
        block <- c(block, lines[[i]])
        i <- i + 1L
      }
      raw <- sub("^#' ?", "", block)
      active <- example_policy_active_lines(raw)
      text <- paste(raw, collapse = "\n")
      active_text <- paste(active, collapse = "\n")
      remaining <- lines[seq.int(start + 1L, length(lines))]
      assignment <- grep(
        "^[.A-Za-z][A-Za-z0-9._]*\\s*<-\\s*function\\s*\\(",
        remaining,
        perl = TRUE
      )
      target <- if (length(assignment) == 0L) {
        NA_character_
      } else {
        sub(
          "^([.A-Za-z][A-Za-z0-9._]*)\\s*<-.*$",
          "\\1",
          remaining[[assignment[[1]]]],
          perl = TRUE
        )
      }
      rows[[length(rows) + 1L]] <- data.frame(
        file = sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path),
        line = start,
        target = target,
        active_fit_calls = sum(grepl("fit_mfrm\\s*\\(", active, perl = TRUE)),
        has_mml = grepl('method\\s*=\\s*["\']MML["\']', text, perl = TRUE),
        active_has_mml = grepl('method\\s*=\\s*["\']MML["\']', active_text, perl = TRUE),
        has_quad_points = grepl("quad_points\\s*=", text),
        active_has_quad_points = grepl("quad_points\\s*=", active_text),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

example_policy_active_lines <- function(lines) {
  active <- character()
  guard_depth <- 0L
  for (line in lines) {
    if (grepl("^\\s*#", line, perl = TRUE)) next
    if (guard_depth > 0L) {
      guard_depth <- max(0L, guard_depth + example_policy_count(line, "{") -
                           example_policy_count(line, "}"))
      next
    }

    if (grepl("\\\\donttest\\s*\\{|\\\\dontrun\\s*\\{|if\\s*\\(\\s*interactive\\s*\\(\\s*\\)\\s*\\)\\s*\\{",
              line, perl = TRUE)) {
      guard_depth <- max(0L, example_policy_count(line, "{") -
                           example_policy_count(line, "}"))
      next
    }

    active <- c(active, line)
  }
  active
}

example_policy_count <- function(x, pattern) {
  matches <- gregexpr(pattern, x, fixed = TRUE)[[1]]
  if (length(matches) == 1L && matches[[1]] == -1L) 0L else length(matches)
}

test_that("example policy excludes comments without hiding executable calls", {
  active <- example_policy_active_lines(c(
    '# After fitting with method = "MML":',
    '# fit_mfrm(data, method = "MML", parallel = TRUE)',
    '# \\donttest{',
    'fit_mfrm(data, method = "MML")',
    '\\donttest{',
    '# }',
    'fit_mfrm(other_data, method = "MML")',
    '}',
    'plot(fit)'
  ))
  expect_identical(active, c('fit_mfrm(data, method = "MML")', 'plot(fit)'))
})

example_policy_hits <- function(rows) {
  if (nrow(rows) == 0L) return(character(0))
  paste0(rows$file, ":", rows$line)
}

example_policy_rd_pages_with <- function(pkg_root, pattern) {
  files <- list.files(
    file.path(pkg_root, "man"),
    pattern = "\\.Rd$",
    full.names = TRUE
  )
  hits <- vapply(files, function(path) {
    any(grepl(pattern, readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1))
  sort(basename(files[hits]))
}

test_that("Rd execution guards have a narrow semantic allowlist", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  expect_identical(
    example_policy_rd_pages_with(pkg_root, "\\dontrun{"),
    sort(unique(paste0(subset(read.csv(file.path(pkg_root, "inst", "examples",
                                  "execution-policy.csv")), execution %in% c("dontrun", "mixed"))$help_topic, ".Rd"))),
    info = paste(
      "Optional examples must have an explicit execution-policy entry.",
      "Supplementary help examples are replayed in the full local suite."
    )
  )
  expect_identical(
    example_policy_rd_pages_with(pkg_root, "# examplesIf"),
    "launch_mfrmr_viewer.Rd",
    info = paste(
      "examplesIf interactive() is reserved for the local Shiny viewer;",
      "ordinary reports, tables, exports, and plots must remain checkable."
    )
  )
})

test_that("CRAN testthat surface is an explicit representative whitelist", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  path <- file.path(pkg_root, "tests", "testthat.R")
  lines <- readLines(path, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  expect_true(grepl("cran_light_tests", text, fixed = TRUE))
  expect_true(grepl("cran_light_filter", text, fixed = TRUE))
  expect_true(grepl('"(^|/)(test-)?("', text, fixed = TRUE))
  expect_true(grepl('")$"', text, fixed = TRUE))
  expect_false(grepl("invert\\s*=\\s*TRUE", text, perl = TRUE))

  expected <- c(
    "cran-smoke",
    "data-and-citation",
    "mml-cpp11-backend",
    "bundle-summary-privacy",
    "likelihood-numerics",
    "gpcm-inference-extensions",
    "gmfrm-em",
    "gmfrm-public-workflow",
    "random-rater",
    "testlet",
    "cml",
    "namespace-contract",
    "vignette-artifacts"
  )
  # Check membership, not a textual occurrence in the extended-suite list.
  selection <- new.env(parent = baseenv())
  for (expr in parse(text = text)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        as.character(expr[[2L]]) %in% c("cran_light_tests", "extended_regression_tests"))
      eval(expr, selection)
  }
  for (slug in expected) {
    expect_true(slug %in% selection$cran_light_tests, info = slug)
  }
  expect_true("gmfrm-portable" %in% selection$extended_regression_tests)
  expect_length(intersect(selection$cran_light_tests, selection$extended_regression_tests), 0L)
  smoke_text <- paste(
    readLines(
      file.path(pkg_root, "tests", "testthat", "test-cran-smoke.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_true(grepl('load_mfrmr_data\\("example_operational"\\)', smoke_text))
  for (route in c(
    'expected_design =',
    'method = "MML"',
    'type = "wright"',
    'fit_stat = "Infit"',
    'include_person = TRUE',
    'export_mfrm_tables('
  )) {
    expect_true(grepl(route, smoke_text, fixed = TRUE))
  }
})

test_that("front-door examples use one explicit beginner workflow", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  rd_text <- function(page) {
    rd <- tools::parse_Rd(file.path(pkg_root, "man", page))
    examples <- rd[vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\examples"),
                          logical(1))]
    paste(unlist(examples, use.names = FALSE), collapse = "")
  }

  load_text <- rd_text("load_mfrmr_data.Rd")
  expect_true(grepl('load_mfrmr_data("example_operational")', load_text, fixed = TRUE))

  data_text <- rd_text("describe_mfrm_data.Rd")
  for (snippet in c(
    'load_mfrmr_data("example_operational")',
    'person = "Person"', 'facets = c("Rater", "Criterion")',
    'score = "Score"', 'rating_min = 1', 'rating_max = 4'
  )) {
    expect_true(grepl(snippet, data_text, fixed = TRUE), info = snippet)
  }

  fit_pages <- c(
    "fit_mfrm.Rd", "summary.mfrm_fit.Rd", "mfrm_results.Rd",
    "mfrm_report.Rd", "export_mfrm_results.Rd",
    "mfrmr_workflow_methods.Rd", "mfrmr_reports_and_tables.Rd",
    "mfrmr-package.Rd", "diagnose_mfrm.Rd",
    "summary.mfrm_diagnostics.Rd", "plot.mfrm_fit.Rd",
    "data_quality_report.Rd", "rating_scale_table.Rd",
    "precision_review_report.Rd", "as.data.frame.mfrm_fit.Rd",
    "mfrmr_visual_diagnostics.Rd", "mfrmr_reporting_and_apa.Rd",
    "compare_mfrm.Rd", "plot_compare_mfrm.Rd", "plot_marginal_fit.Rd",
    "plot_marginal_pairwise.Rd", "plot_qc_dashboard.Rd", "plot_person_fit.Rd",
    "compute_person_fit_indices.Rd", "plot_rater_severity_profile.Rd",
    "plot_bubble.Rd", "fair_average_table.Rd", "plot_fair_average.Rd",
    "interrater_agreement_table.Rd", "plot_interrater_agreement.Rd",
    "apa_table.Rd", "export_mfrm_tables.Rd", "reporting_checklist.Rd",
    "build_apa_outputs.Rd", "plot_data.Rd", "plot_data_components.Rd"
  )
  for (page in fit_pages) {
    text <- rd_text(page)
    if (grepl('readRDS\\(system.file\\("examples", "(core-jml|bias-jml|operational-mml)[.]rds"', text)) {
      expect_true(grepl("example_fit$fit", text, fixed = TRUE), info = page)
      next
    }
    for (snippet in c(
      'load_mfrmr_data("example_operational")',
      'person = "Person"', 'facets = c("Rater", "Criterion")',
      'score = "Score"', 'method = "MML"', 'model = "RSM"'
    )) {
      expect_true(grepl(snippet, text, fixed = TRUE), info = paste(page, snippet))
    }
    expect_false(grepl("quad_points\\s*=|maxit\\s*=|reltol\\s*=", text), info = page)
    expect_false(grepl("JML keeps the help example fast", text, fixed = TRUE), info = page)
  }
  for (page in c("reporting_checklist.Rd", "build_apa_outputs.Rd")) {
    text <- rd_text(page)
    expect_false(grepl('DraftReady = "yes"|ready to paste', text), info = page)
  }
  expect_false(grepl("panels$Status", rd_text("plot_qc_dashboard.Rd"), fixed = TRUE))
})

test_that("repeated fitting in examples has an explicit comparison purpose", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  examples <- example_policy_roxygen_examples(pkg_root)
  testthat::skip_if(is.null(examples) || nrow(examples) == 0L,
                    "roxygen examples are not available")

  multi_fit <- examples[examples$active_fit_calls > 1L, ]
  # Linking/drift/trajectory need two occasions; weighting, model choice and
  # centrality need two fits. Intervals compare fixed-SD MML, free-SD MML and CML.
  expected <- c(detect_anchor_drift = 2L, build_equating_chain = 2L,
    build_linking_review = 2L, build_weighting_review = 2L,
    build_model_choice_review = 2L, compare_mfrm = 2L,
    mfrm_facet_intervals = 3L, plot_compare_mfrm = 2L,
    plot.mfrm_anchor_drift = 2L, plot_rater_trajectory = 2L,
    rater_centrality_table = 2L)
  observed <- setNames(multi_fit$active_fit_calls, multi_fit$target)
  expect_true(all(names(observed) %in% names(expected)))
  expect_identical(observed[order(names(observed))], expected[sort(names(observed))])
  # residual PCA's parallel=TRUE requests parallel analysis, not CPU workers.
})

test_that("vignettes keep heavy executable chunks off during CRAN checks", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  files <- list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$",
                      full.names = TRUE)
  testthat::skip_if(length(files) == 0L, "vignettes are not available")

  missing_guard <- character(0)
  unsafe_cache_chunk <- character(0)
  # These short observed-score tutorials execute in full without MFRM calibration.
  fully_executed <- c("mfrmr-gtheory.Rmd", "mfrmr-rating-agreement.Rmd")
  for (path in files) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    if (basename(path) == "mfrmr-model-choice.Rmd") {
      # This reading guide has no executable R; fail if that contract changes.
      expect_false(grepl("(?m)^```\\{r|`r[[:space:]]", text, perl = TRUE))
      next
    }
    inline_guard <- grepl('eval = isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))', text, fixed = TRUE)
    if (!basename(path) %in% fully_executed && !inline_guard &&
        (!grepl("is_cran_check", text, fixed = TRUE) ||
         !grepl("eval = !is_cran_check", text, fixed = TRUE))) {
      missing_guard <- c(missing_guard, basename(path))
    }
    if (grepl("eval = is_cran_check", text, fixed = TRUE) &&
        !grepl("vignette_artifact", text, fixed = TRUE)) {
      unsafe_cache_chunk <- c(unsafe_cache_chunk, basename(path))
    }
  }

  expect_identical(missing_guard, character(0))
  expect_identical(unsafe_cache_chunk, character(0))
})
