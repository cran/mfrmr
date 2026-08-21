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
        outer_donttest = length(raw) > 0L &&
          identical(trimws(raw[[1]]), "\\donttest{"),
        donttest_open_count = sum(grepl(
          "^\\s*\\\\donttest\\s*\\{",
          raw,
          perl = TRUE
        )),
        active_fit_calls = sum(grepl("fit_mfrm\\s*\\(", active, perl = TRUE)),
        has_mml = grepl('method\\s*=\\s*["\']MML["\']', text, perl = TRUE),
        active_has_mml = grepl('method\\s*=\\s*["\']MML["\']', active_text, perl = TRUE),
        has_quad_points = grepl("quad_points\\s*=", text),
        active_has_quad_points = grepl("quad_points\\s*=", active_text),
        active_has_high_maxit = grepl("maxit\\s*=\\s*([3-9][1-9]|[4-9][0-9]|[1-9][0-9]{2,})", active_text),
        active_has_parallel = grepl("parallel\\s*=\\s*TRUE", active_text),
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
    sort(c(
      "normalize_conquest_overlap_exports.Rd",
      "review_conquest_overlap.Rd"
    )),
    info = paste(
      "dontrun is reserved for examples that require separately generated",
      "external ConQuest output. Executable lengthy examples belong in donttest."
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

test_that("comparatively expensive demonstrations run in the donttest pass", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  expected <- sort(c(
    "gpcm_mml_quadrature_sensitivity",
    "dif_interaction_table", "plot_dif_heatmap", "facet_quality_dashboard",
    "summary.mfrm_facet_dashboard", "plot_facet_quality_dashboard",
    "diagnose_mfrm", "build_mfrm_manifest", "build_mfrm_replay_script",
    "export_mfrm_bundle", "analyze_facet_equivalence",
    "plot_facet_equivalence", "summary.mfrm_diagnostics",
    "summary.mfrm_bias", "compute_person_fit_indices", "plot_person_fit",
    "plot_rater_severity_profile", "plot_dif_summary",
    "plot_apa_figure_one", "plot_guttman_scalogram", "plot_residual_qq",
    "plot_rater_agreement_heatmap", "plot_local_dependence_heatmap",
    "plot_reliability_snapshot", "plot_residual_matrix",
    "plot_shrinkage_funnel", "plot_unexpected", "plot_fair_average",
    "plot_displacement", "plot_facets_chisq", "plot_qc_dashboard",
    "plot_bubble", "export_mfrm", "q3_statistic", "reporting_checklist",
    "subset_connectivity_report", "build_peer_review_design_review",
    "facet_statistics_report", "precision_review_report",
    "category_structure_report", "category_curves_report",
    "bias_interaction_report", "bias_iteration_report",
    "bias_pairwise_report", "plot_bias_interaction", "build_apa_outputs",
    "print.mfrm_apa_text", "summary.mfrm_apa_outputs", "summary.apa_table",
    "facets_fit_review", "dif_report", "apply_empirical_bayes_shrinkage",
    "shrinkage_report", "build_peer_review_sim_spec",
    "interrater_agreement_table", "facets_chisq_table",
    "unexpected_response_table", "displacement_table",
    "measurable_summary_table", "rating_scale_table", "bias_count_table",
    "unexpected_after_bias_table", "facets_output_file_bundle",
    "plot_residual_pca", "estimate_bias", "summary.mfrm_facets_run"
  ))
  examples <- example_policy_roxygen_examples(pkg_root)
  observed <- sort(unique(examples$target[examples$outer_donttest]))
  missing <- setdiff(expected, observed)
  nested <- examples[
    examples$target %in% expected & examples$donttest_open_count != 1L,
    c("file", "line", "target", "donttest_open_count"),
    drop = FALSE
  ]

  expect_identical(
    missing,
    character(0),
    info = paste(
      "These noninteractive but comparatively expensive examples must remain",
      "executable in the explicit --run-donttest pass."
    )
  )
  expect_equal(
    nrow(nested),
    0L,
    info = "Former interactive examples should have one outer donttest guard."
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
    "compatibility-aliases",
    "data-and-citation",
    "example-datasets",
    "mml-cpp11-backend",
    "missing-codes-integration",
    "bundle-summary-privacy",
    "gpcm-capability-matrix",
    "namespace-contract",
    "vignette-artifacts"
  )
  for (slug in expected) {
    expect_true(grepl(paste0('"', slug, '"'), text, fixed = TRUE))
  }
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
    'export_mfrm('
  )) {
    expect_true(grepl(route, smoke_text, fixed = TRUE))
  }
})

test_that("roxygen examples keep expensive demonstrations conditional", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  examples <- example_policy_roxygen_examples(pkg_root)
  testthat::skip_if(is.null(examples) || nrow(examples) == 0L,
                    "roxygen examples are not available")

  multi_fit <- examples[examples$active_fit_calls > 1L, ]
  expect_identical(
    example_policy_hits(multi_fit),
    character(0),
    info = "Standard Rd examples should not run multiple fit_mfrm() calls."
  )

  mml_without_quadrature <- examples[
    examples$has_mml & !examples$has_quad_points, ]
  expect_identical(
    example_policy_hits(mml_without_quadrature),
    character(0),
    info = "MML examples should set quad_points, including donttest examples."
  )

  active_mml_without_quadrature <- examples[
    examples$active_has_mml & !examples$active_has_quad_points, ]
  expect_identical(
    example_policy_hits(active_mml_without_quadrature),
    character(0),
    info = "Standard MML examples should set quad_points."
  )

  high_maxit <- examples[examples$active_has_high_maxit, ]
  expect_identical(
    example_policy_hits(high_maxit),
    character(0),
    info = "Standard Rd examples should keep maxit at 30 or below."
  )

  active_parallel <- examples[examples$active_has_parallel, ]
  expect_identical(
    example_policy_hits(active_parallel),
    character(0),
    info = "Parallel diagnostics should be in donttest or vignettes."
  )
})

test_that("vignettes keep heavy executable chunks off during CRAN checks", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  files <- list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$",
                      full.names = TRUE)
  testthat::skip_if(length(files) == 0L, "vignettes are not available")

  missing_guard <- character(0)
  unsafe_cache_chunk <- character(0)
  for (path in files) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    if (!grepl("is_cran_check", text, fixed = TRUE) ||
        !grepl("eval = !is_cran_check", text, fixed = TRUE)) {
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
