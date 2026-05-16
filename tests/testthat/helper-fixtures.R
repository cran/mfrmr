# Shared toy-fit fixtures for testthat.
#
# testthat sources every `helper-*.R` file in this directory before
# running tests, so anything defined here is visible to all
# `test-*.R` files. The goal is to keep test setup focused on what
# is actually being tested rather than re-typing the same
# `load_mfrmr_data() + fit_mfrm() + diagnose_mfrm()` chain.
#
# Each helper is deterministic and fast (the JML toy fit takes
# well under a second on example_core / example_bias). Caching
# inside the helpers avoids paying that cost more than once when
# multiple tests in the same session ask for the same fixture.

# Internal cache, scoped to a list so concurrent test files do not
# collide. We deliberately do not export this; tests should call
# the public helpers below.
.mfrmr_test_cache <- new.env(parent = emptyenv())

.cache_key <- function(...) paste(vapply(list(...), function(x) {
  if (is.null(x)) "NULL" else paste(deparse(x), collapse = "|")
}, character(1)), collapse = ":")

#' Toy MFRM fit on `example_core` for tests
#'
#' Returns a deterministic JML fit suitable for smoke / contract
#' tests. The result is cached for repeated calls with the same
#' arguments.
make_toy_fit <- function(method = "JML",
                         maxit = 25,
                         model = "RSM",
                         attach_diagnostics = FALSE,
                         dataset = "example_core") {
  key <- .cache_key("fit", dataset, method, maxit, model, attach_diagnostics)
  if (exists(key, envir = .mfrmr_test_cache, inherits = FALSE)) {
    return(get(key, envir = .mfrmr_test_cache, inherits = FALSE))
  }
  toy <- load_mfrmr_data(dataset)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy,
             person = "Person",
             facets = c("Rater", "Criterion"),
             score = "Score",
             method = method,
             model = model,
             maxit = maxit,
             attach_diagnostics = attach_diagnostics)
  ))
  assign(key, fit, envir = .mfrmr_test_cache)
  fit
}

#' Toy diagnostics bundle for the same `example_core` fit
make_toy_diagnostics <- function(fit = NULL,
                                 diagnostic_mode = "legacy",
                                 residual_pca = "none") {
  if (is.null(fit)) fit <- make_toy_fit()
  key <- .cache_key("diag", attr(fit, "config")$method %||% "JML",
                    diagnostic_mode, residual_pca)
  if (exists(key, envir = .mfrmr_test_cache, inherits = FALSE)) {
    return(get(key, envir = .mfrmr_test_cache, inherits = FALSE))
  }
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(fit,
                  diagnostic_mode = diagnostic_mode,
                  residual_pca = residual_pca)
  ))
  assign(key, diag, envir = .mfrmr_test_cache)
  diag
}

#' Loads the toy data + fit + diagnostics into the calling test_that
#' block under the conventional `.toy` / `.fit` / `.diag` names.
#'
#' Tests that follow the `.toy / .fit / .diag` convention can replace
#' an opening `local({ ... })` block with a single call to
#' `local_toy_fit()` inside `test_that()`.
local_toy_fit <- function(envir = parent.frame(),
                          method = "JML",
                          maxit = 25,
                          model = "RSM",
                          dataset = "example_core",
                          diagnostic_mode = "legacy") {
  envir$.toy <- load_mfrmr_data(dataset)
  envir$.fit <- make_toy_fit(method = method, maxit = maxit,
                              model = model, dataset = dataset)
  envir$.diag <- make_toy_diagnostics(envir$.fit,
                                       diagnostic_mode = diagnostic_mode)
  invisible(envir)
}
