compiled_header_source_root <- function() {
  test_root <- normalizePath(
    testthat::test_path(),
    winslash = "/",
    mustWork = TRUE
  )
  candidates <- unique(normalizePath(c(
    file.path(test_root, "..", ".."),
    file.path(test_root, "..", "..", "00_pkg_src", "mfrmr"),
    file.path(test_root, "..", "..", "..", "00_pkg_src", "mfrmr"),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "00_pkg_src", "mfrmr"),
    file.path(getwd(), "..", "..", "00_pkg_src", "mfrmr")
  ), winslash = "/", mustWork = FALSE))

  roots <- candidates[
    file.exists(file.path(candidates, "DESCRIPTION")) &
      dir.exists(file.path(candidates, "src"))
  ]
  if (length(roots) == 0L) NA_character_ else roots[[1L]]
}

compiled_header_has_override <- function(lines) {
  any(grepl("HAVE_ENUM_BASE_TYPE", lines, fixed = TRUE))
}

test_that("compiled-header override detector covers source and build flags", {
  expect_true(all(vapply(c(
    "#define HAVE_ENUM_BASE_TYPE 1",
    "# undef HAVE_ENUM_BASE_TYPE",
    "PKG_CPPFLAGS=-DHAVE_ENUM_BASE_TYPE=1",
    "PKG_CPPFLAGS = -U HAVE_ENUM_BASE_TYPE"
  ), compiled_header_has_override, logical(1))))
  expect_false(compiled_header_has_override(c(
    "#include <Rconfig.h>",
    "PKG_CPPFLAGS = -DNDEBUG"
  )))
})

test_that("compiled sources preserve R's configured header types", {
  pkg_root <- compiled_header_source_root()
  testthat::skip_if(is.na(pkg_root), "compiled sources are not available")

  sources <- list.files(
    file.path(pkg_root, "src"),
    pattern = "[.](c|cc|cpp|cxx|h|hh|hpp)$",
    recursive = TRUE,
    full.names = TRUE
  )
  expect_gt(length(sources), 0L)
  makevars <- list.files(
    file.path(pkg_root, "src"),
    pattern = "^Makevars",
    full.names = TRUE
  )
  scanned <- c(sources, makevars)
  overrides <- vapply(scanned, function(path) {
    compiled_header_has_override(readLines(path, warn = FALSE))
  }, logical(1))

  expect_identical(
    basename(scanned[overrides]),
    character(0),
    info = paste(
      "All compiled translation units must use the Rconfig.h definition of",
      "HAVE_ENUM_BASE_TYPE so Rboolean and R_UnwindProtect remain ODR-safe",
      "under link-time optimization."
    )
  )
})
