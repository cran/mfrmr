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
               dir.exists(file.path(candidates, "R"))][1]
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
      rows[[length(rows) + 1L]] <- data.frame(
        file = sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path),
        line = start,
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

test_that("CRAN testthat surface is an explicit lightweight whitelist", {
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
    "compatibility-aliases",
    "data-and-citation",
    "gpcm-capability-matrix",
    "namespace-contract"
  )
  for (slug in expected) {
    expect_true(grepl(paste0('"', slug, '"'), text, fixed = TRUE))
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

test_that("vignettes keep executable chunks off during CRAN checks", {
  pkg_root <- example_policy_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  files <- list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$",
                      full.names = TRUE)
  testthat::skip_if(length(files) == 0L, "vignettes are not available")

  missing_guard <- character(0)
  for (path in files) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    if (!grepl("is_cran_check", text, fixed = TRUE) ||
        !grepl("eval = !is_cran_check", text, fixed = TRUE)) {
      missing_guard <- c(missing_guard, basename(path))
    }
  }

  expect_identical(missing_guard, character(0))
})
