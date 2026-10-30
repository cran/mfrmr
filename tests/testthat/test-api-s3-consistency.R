api_consistency_source_root <- function() {
  test_root <- normalizePath(testthat::test_path(), mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    file.path(test_root, "..", ".."),
    getwd(),
    file.path(getwd(), "..")
  ), mustWork = FALSE))
  matches <- candidates[
    file.exists(file.path(candidates, "DESCRIPTION")) &
      file.exists(file.path(candidates, "NAMESPACE")) &
      dir.exists(file.path(candidates, "man"))
  ]
  if (length(matches) == 0L) NA_character_ else matches[1L]
}

test_that("registered print summary and plot methods use generic signatures", {
  registry <- as.data.frame(
    getNamespaceInfo(asNamespace("mfrmr"), "S3methods"),
    stringsAsFactors = FALSE
  )
  names(registry)[seq_len(3L)] <- c("generic", "class", "method")
  registry <- registry[
    registry$generic %in% c("print", "summary", "plot"), , drop = FALSE
  ]

  expect_gt(nrow(registry), 0L)
  for (i in seq_len(nrow(registry))) {
    method <- get(registry$method[i], envir = asNamespace("mfrmr"))
    arguments <- names(formals(method))
    expected_first <- if (registry$generic[i] == "summary") "object" else "x"
    expect_identical(
      arguments[1L], expected_first,
      info = paste(registry$generic[i], registry$class[i], sep = ".")
    )
    expect_true(
      "..." %in% arguments,
      info = paste(registry$generic[i], registry$class[i], sep = ".")
    )
  }
})

test_that("portable calibration methods have direct help and semantic scope", {
  root <- api_consistency_source_root()
  skip_if(is.na(root), "source help files are not available")

  expected <- list(
    mfrm_calibration_methods.Rd = c(
      "summary.mfrm_calibration", "print.mfrm_calibration",
      "print.summary.mfrm_calibration"
    ),
    mfrm_calibration_score_methods.Rd = c(
      "summary.mfrm_calibration_score", "print.mfrm_calibration_score",
      "print.summary.mfrm_calibration_score", "plot.mfrm_calibration_score"
    ),
    mfrm_dff_methods.Rd = c(
      "summary.mfrm_dif", "summary.mfrm_dff", "print.mfrm_dif",
      "print.mfrm_dff", "print.summary.mfrm_dif",
      "print.summary.mfrm_dff"
    )
  )
  for (file in names(expected)) {
    text <- paste(readLines(file.path(root, "man", file), warn = FALSE),
                  collapse = "\n")
    for (method in expected[[file]]) {
      expect_match(
        text, paste0("\\alias{", method, "}"), fixed = TRUE,
        info = paste(file, method)
      )
    }
  }

  artifact_help <- paste(readLines(
    file.path(root, "man", "mfrm_calibration_methods.Rd"), warn = FALSE
  ), collapse = "\n")
  expect_match(
    artifact_help, "intentionally has no \\code{plot()} method", fixed = TRUE
  )
  expect_null(getS3method("plot", "mfrm_calibration", optional = TRUE))
})
