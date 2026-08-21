test_that("precomputed vignette artifacts match their semantic manifest", {
  artifact_dir <- system.file(
    "extdata", "vignette-artifacts",
    package = "mfrmr",
    mustWork = TRUE
  )
  manifest <- utils::read.csv(
    file.path(artifact_dir, "manifest.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_true(nrow(manifest) > 0L)
  expect_true(all(c(
    "Artifact", "Rows", "Columns", "Schema", "Source",
    "DataKey", "GeneratedWith"
  ) %in% names(manifest)))
  expect_false(any(grepl("sha|md5|hash", names(manifest), ignore.case = TRUE)))
  generator_path <- testthat::test_path(
    "..", "..", "inst", "validation", "generate-vignette-artifacts.R"
  )
  if (file.exists(generator_path)) {
    generator_text <- paste(
      readLines(generator_path, warn = FALSE), collapse = "\n"
    )
    hash_terms <- "\\b(?:md5|sha(?:1|224|256|384|512)?|hash(?:es|ed|ing)?)\\b"
    expect_false(grepl(
      hash_terms, generator_text, ignore.case = TRUE, perl = TRUE
    ))
  }
  expect_true(all(manifest$DataKey == "example_operational"))
  expect_true(all(manifest$GeneratedWith == as.character(utils::packageVersion("mfrmr"))))

  for (i in seq_len(nrow(manifest))) {
    path <- file.path(artifact_dir, manifest$Artifact[i])
    expect_true(file.exists(path), info = manifest$Artifact[i])
    artifact <- utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = character()
    )
    expect_equal(nrow(artifact), manifest$Rows[i], info = manifest$Artifact[i])
    expect_equal(ncol(artifact), manifest$Columns[i], info = manifest$Artifact[i])
    expect_identical(
      paste(names(artifact), collapse = "|"),
      manifest$Schema[i],
      info = manifest$Artifact[i]
    )
  }
})

test_that("precomputed workflow artifacts represent the canonical successful fit", {
  artifact_dir <- system.file(
    "extdata", "vignette-artifacts",
    package = "mfrmr",
    mustWork = TRUE
  )
  fit_overview <- utils::read.csv(
    file.path(artifact_dir, "workflow_fit_overview.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  diagnostic_overview <- utils::read.csv(
    file.path(artifact_dir, "workflow_diagnostic_overview.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fit_decision <- utils::read.csv(
    file.path(artifact_dir, "workflow_fit_decision.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_identical(fit_overview$Model, "RSM")
  expect_identical(fit_overview$Method, "MML")
  expect_true(isTRUE(fit_overview$Converged))
  expect_true(isTRUE(fit_overview$InferenceReady))
  expect_identical(fit_overview$ConvergenceSeverity, "pass")
  expect_identical(fit_decision$Interpretation,
                   "Fit gates passed; formal precision review required")
  expect_identical(fit_decision$FormalInference, "No")
  expect_identical(fit_decision$FitReadiness, "ready")
  expect_identical(fit_decision$Why,
                   "Formal precision support has not been evaluated.")
  expect_match(fit_decision$NextAction, "diagnose_mfrm()", fixed = TRUE)
  expect_true(all(c(
    "OptimizerInitialMethod", "OptimizerMethod", "OptimizerPolished",
    "RequestedReltol", "EffectiveReltol", "OptimizerFactr",
    "OptimizerPgtol"
  ) %in% names(fit_overview)))
  expect_identical(fit_overview$OptimizerInitialMethod, "L-BFGS-B")
  expect_true(isTRUE(fit_overview$OptimizerPolished))
  expect_equal(fit_overview$RequestedReltol, 1e-9)
  expect_lt(fit_overview$EffectiveReltol, fit_overview$RequestedReltol)
  expect_true(is.finite(fit_overview$OptimizerFactr))
  expect_true(is.finite(fit_overview$OptimizerPgtol))
  expect_identical(diagnostic_overview$Method, "MML")
  expect_identical(diagnostic_overview$PrecisionTier, "model_based")
})

test_that("public workflow artifacts do not expose build-machine paths", {
  artifact_dir <- system.file(
    "extdata", "vignette-artifacts",
    package = "mfrmr",
    mustWork = TRUE
  )
  export_files <- utils::read.csv(
    file.path(artifact_dir, "workflow_export_files.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_identical(export_files$Path, basename(export_files$Path))
  expect_false(any(grepl(
    "([[:alpha:]]:|/)([/\\\\]|[^[:space:]])*(Users|AppData|Rtmp|Temp)([/\\\\]|$)",
    export_files$Path,
    ignore.case = TRUE,
    perl = TRUE
  )))

  workflow_file <- system.file(
    "doc", "mfrmr-workflow.html",
    package = "mfrmr",
    mustWork = FALSE
  )
  if (!nzchar(workflow_file)) {
    workflow_file <- testthat::test_path(
      "..", "..", "vignettes", "mfrmr-workflow.Rmd"
    )
  }
  skip_if_not(
    file.exists(workflow_file),
    "Built vignette and source vignette are unavailable in this check mode."
  )
  expect_true(file.exists(workflow_file))
  html <- paste(readLines(workflow_file, warn = FALSE), collapse = "\n")
  expect_false(grepl(
    "[[:alpha:]]:[/\\\\](Users|Documents and Settings)[/\\\\]",
    html,
    ignore.case = TRUE,
    perl = TRUE
  ))
  expect_false(grepl(
    "[[:alpha:]]:[/\\\\].*(AppData|Rtmp)[/\\\\]",
    html,
    ignore.case = TRUE,
    perl = TRUE
  ))
})
