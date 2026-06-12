documentation_source_root <- function() {
  test_root <- normalizePath(testthat::test_path(), mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    file.path(test_root, "..", ".."),
    file.path(test_root, "..", "..", "00_pkg_src", "mfrmr"),
    file.path(test_root, "..", "..", "..", "00_pkg_src", "mfrmr"),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "00_pkg_src", "mfrmr"),
    file.path(getwd(), "..", "..", "00_pkg_src", "mfrmr")
  ), mustWork = FALSE))

  candidates[file.exists(file.path(candidates, "README.md")) &
               dir.exists(file.path(candidates, "vignettes"))][1]
}

test_that("high-level help pages expose searchable concepts", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  expected <- list(
    "mfrmr_interval_guide.Rd" = c("confidence intervals", "uncertainty displays"),
    "mfrmr_visual_diagnostics.Rd" = c("visual diagnostics", "confidence intervals"),
    "visual_reporting_template.Rd" = c("visual diagnostics", "figure captions"),
    "mfrmr_output_guide.Rd" = c("reporting workflow", "route selection"),
    "gpcm_capability_matrix.Rd" = c("GPCM boundaries", "route selection"),
    "fit_measures_table.Rd" = c("confidence intervals"),
    "plot_fair_average.Rd" = c("confidence intervals"),
    "plot_bias_interaction.Rd" = c("confidence intervals"),
    "plot_wright_unified.Rd" = c("Wright maps"),
    "plot.mfrm_fit.Rd" = c("shrinkage")
  )
  missing <- character(0)
  for (rd in names(expected)) {
    text <- paste(readLines(file.path(pkg_root, "man", rd), warn = FALSE),
                  collapse = "\n")
    for (concept in expected[[rd]]) {
      tag <- paste0("\\concept{", concept, "}")
      if (!grepl(tag, text, fixed = TRUE)) {
        missing <- c(missing, paste(rd, concept, sep = ": "))
      }
    }
  }
  expect_identical(missing, character(0))
})

test_that("documentation keeps RSM/PCM and bounded GPCM wording separated", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  doc_files <- c(
    file.path(pkg_root, "README.md"),
    list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE)
  )
  doc_files <- doc_files[file.exists(doc_files)]
  doc_files <- doc_files[basename(doc_files) != "test-documentation-terminology.R"]
  expect_true(length(doc_files) > 0)

  stale_phrases <- c(
    "Native R package for many-facet Rasch model",
    "A many-facet Rasch model (MFRM) was fit",
    "This function generates synthetic MFRM data from the Rasch model",
    "under the same `RSM` / `PCM` interface",
    "Rasch-MFRM",
    "still-unvalidated score-side",
    "still-blocked score-side",
    "should still be treated as unsupported",
    "current public response-model scope",
    "fit-based bundle/export contract still depends",
    "current public RSM/PCM release",
    "fair-average, APA writer, and broader planning semantics remain generalized only"
  )

  hits <- character(0)
  for (path in doc_files) {
    lines <- readLines(path, warn = FALSE)
    for (phrase in stale_phrases) {
      idx <- grep(phrase, lines, fixed = TRUE)
      if (length(idx) > 0L) {
        rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path)
        hits <- c(hits, paste0(rel, ":", idx, ": ", phrase))
      }
    }
  }

  expect_identical(hits, character(0))
})

test_that("user-facing model choice guide is present", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  read_file <- function(path) paste(readLines(file.path(pkg_root, path), warn = FALSE), collapse = "\n")

  readme <- read_file("README.md")
  expect_true(grepl("Model selection guide and report wording", readme, fixed = TRUE))
  expect_true(grepl("bounded `GPCM` to inspect whether allowing discrimination-based reweighting", readme, fixed = TRUE))

  gpcm_vignette <- read_file("vignettes/mfrmr-gpcm-scope.Rmd")
  expect_true(grepl("Before fitting: model-choice triage", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("Report wording templates", gpcm_vignette, fixed = TRUE))

  mml_vignette <- read_file("vignettes/mfrmr-mml-and-marginal-fit.Rmd")
  expect_true(grepl("Practical Reading Order", mml_vignette, fixed = TRUE))
  expect_true(grepl("fit improvement alone decide the operational model", mml_vignette, fixed = TRUE))
})

test_that("FACETS positioning docs avoid complete-reproduction claims", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  read_file <- function(path) paste(readLines(file.path(pkg_root, path), warn = FALSE), collapse = "\n")

  readme <- read_file("README.md")
  mapping <- read_file(file.path("inst", "references", "FACETS_manual_mapping.md"))
  positioning <- read_file(file.path("man", "facets_positioning_guide.Rd"))
  combined <- paste(readme, mapping, positioning, sep = "\n")

  expect_true(grepl("results remain\n  `mfrmr` estimates", readme, fixed = TRUE))
  expect_true(grepl("source of truth unless external FACETS output is\nexplicitly supplied for comparison", mapping, fixed = TRUE))
  expect_true(grepl("not a FACETS numerical\nclone", positioning, fixed = TRUE))
  expect_true(grepl("not evidence that FACETS was executed", combined, fixed = TRUE))
  stale_claims <- c(
    "complete FACETS reproduction",
    "fully reproduce FACETS",
    "FACETS produced the estimates"
  )
  hits <- stale_claims[vapply(stale_claims, grepl, logical(1), x = combined, fixed = TRUE)]
  expect_identical(hits, character(0))
})

test_that("FACETS output-contract route avoids old equivalence terminology", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  doc_files <- c(
    file.path(pkg_root, "README.md"),
    file.path(pkg_root, "NEWS.md"),
    list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "man"), pattern = "\\.Rd$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "inst", "references"), pattern = "\\.(md|csv)$", full.names = TRUE, recursive = TRUE)
  )
  doc_files <- doc_files[file.exists(doc_files)]
  doc_files <- doc_files[basename(doc_files) != "test-documentation-terminology.R"]

  hits <- character(0)
  blocked <- paste(
    c(
      paste0("pa", "rity"),
      paste0("facets_", "pa", "rity", "_report"),
      paste0("mfrm_", "pa", "rity", "_report"),
      paste0("pa", "rity", "_report")
    ),
    collapse = "|"
  )
  for (path in doc_files) {
    lines <- readLines(path, warn = FALSE)
    idx <- grep(blocked, lines, ignore.case = TRUE, perl = TRUE)
    if (length(idx) > 0L) {
      rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path)
      hits <- c(hits, paste0(rel, ":", idx, ": ", lines[idx]))
    }
  }

  expect_identical(hits, character(0))
})

test_that("release validation docs separate recovery decision from uncertainty status", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  read_file <- function(path) paste(readLines(file.path(pkg_root, path), warn = FALSE), collapse = "\n")

  readme <- read_file("README.md")
  workflow <- read_file("vignettes/mfrmr-workflow.Rmd")
  gpcm_vignette <- read_file("vignettes/mfrmr-gpcm-scope.Rmd")
  readme_flat <- gsub("\\s+", " ", readme)
  workflow_flat <- gsub("\\s+", " ", workflow)

  expect_true(grepl("topline_release_decision", readme, fixed = TRUE))
  expect_true(grepl("ReleaseRecoveryStatus", readme, fixed = TRUE))
  expect_true(grepl("ExtendedSensitivityStatus", readme, fixed = TRUE))
  expect_true(grepl("summary(validation)", readme, fixed = TRUE))
  expect_true(grepl("condition_review", readme, fixed = TRUE))
  expect_true(grepl("condition_summary", readme, fixed = TRUE))
  expect_true(grepl("do not treat `OverallStatus = \"review\"` as a release-level", readme_flat, fixed = TRUE))
  expect_true(grepl("SE/coverage evidence is intentionally reported as a separate limitation", readme_flat, fixed = TRUE))
  expect_true(grepl("The recommended reading order is: `summary(recovery_review)`", readme_flat, fixed = TRUE))
  expect_true(grepl("summary(recovery_review)$reading_order", readme, fixed = TRUE))
  expect_true(grepl("`reading_order` and `guidance` fields", readme_flat, fixed = TRUE))

  expect_true(grepl("summary(recovery_review)$condition_review", workflow, fixed = TRUE))
  expect_true(grepl("summary(recovery_review)$diagnostic_reporting_notes", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$condition_reporting_notes", workflow, fixed = TRUE))
  expect_true(grepl("summary(recovery_review)$reading_order", workflow, fixed = TRUE))
  expect_true(grepl("# summary(validation)", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$reading_order", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$topline_release_decision", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$release_decision_table", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$condition_summary", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$diagnostic_reporting_notes", workflow, fixed = TRUE))
  expect_true(grepl("validation_summary$domain_decision_table", workflow, fixed = TRUE))
  expect_true(grepl("validation_bundle$tables$condition_summary", workflow, fixed = TRUE))
  expect_true(grepl("validation_bundle$tables$condition_reporting_notes", workflow, fixed = TRUE))
  expect_true(grepl("validation_bundle$tables$diagnostic_reporting_notes", workflow, fixed = TRUE))
  expect_true(grepl("validation_bundle$tables$reading_order", workflow, fixed = TRUE))
  expect_true(grepl("ExtendedSensitivityStatus reports non-core stress", workflow, fixed = TRUE))
  expect_true(grepl("so OverallStatus = \"review\" is not", workflow, fixed = TRUE))
  expect_true(grepl("status_plot$data$section_status", workflow, fixed = TRUE))
  expect_true(grepl("metric_plot$data$guidance", workflow, fixed = TRUE))

  expect_true(grepl("Source-grounded recovery interpretation", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("slope_regime", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("not literature-derived adequacy cut points", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("recovery_review$condition_review", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("recovery_review$condition_reporting_notes", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("summary(review_gpcm)$reading_order", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("ExtendedSensitivityStatus", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("diagnostic_reporting_notes", gpcm_vignette, fixed = TRUE))
  expect_true(grepl("build_summary_table_bundle(validation_summary)", gpcm_vignette, fixed = TRUE))
})

test_that("GPCM source-grounding docs keep operational labels separate from literature", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  read_file <- function(path) paste(readLines(file.path(pkg_root, path), warn = FALSE), collapse = "\n")

  readme <- read_file("README.md")
  gpcm_vignette <- read_file("vignettes/mfrmr-gpcm-scope.Rmd")
  sim_spec_help <- read_file("man/build_mfrm_sim_spec.Rd")
  evidence_021 <- read_file("inst/validation/release-evidence-map-0.2.1.md")

  combined <- paste(readme, gpcm_vignette, sim_spec_help, evidence_021, sep = "\n")
  expect_true(grepl("Muraki (1992", combined, fixed = TRUE))
  expect_true(grepl("Muraki (1993", combined, fixed = TRUE))
  expect_true(grepl("Morris, White, and Crowther (2019", combined, fixed = TRUE))
  expect_true(grepl("slope_regime", combined, fixed = TRUE))
  expect_true(grepl("not psychometric fit or\nadequacy cut points", readme, fixed = TRUE) ||
                grepl("not psychometric fit or adequacy cut points", gsub("\\s+", " ", readme), fixed = TRUE))
  expect_true(grepl("not model-fit tests and\nthey are not literature-derived adequacy cut points", gpcm_vignette, fixed = TRUE) ||
                grepl("not model-fit tests and they are not literature-derived adequacy cut points",
                      gsub("\\s+", " ", gpcm_vignette), fixed = TRUE))
  expect_true(grepl("Samejima (1974)", evidence_021, fixed = TRUE))
  expect_true(grepl("do not add a\nSamejima normal-ogive or graded-response-model implementation", evidence_021, fixed = TRUE) ||
                grepl("do not add a Samejima normal-ogive or graded-response-model implementation",
                      gsub("\\s+", " ", evidence_021), fixed = TRUE))
})

test_that("fit and separation reporting docs keep diagnostics separate from validation gates", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")
  read_file <- function(path) paste(readLines(file.path(pkg_root, path), warn = FALSE), collapse = "\n")

  readme <- read_file("README.md")
  reporting_vignette <- read_file("vignettes/mfrmr-reporting-and-apa.Rmd")
  precision_help <- read_file("man/precision_review_report.Rd")
  evidence_021 <- read_file("inst/validation/release-evidence-map-0.2.1.md")
  combined <- paste(readme, reporting_vignette, precision_help, evidence_021, sep = "\n")

  expect_true(grepl("fit_separation_basis", combined, fixed = TRUE))
  expect_true(grepl("diagnostic operating characteristics", combined, fixed = TRUE))
  expect_true(grepl("diagnostic_only_not_release_gate", combined, fixed = TRUE))
  expect_true(grepl("Wright and Linacre (1994)", evidence_021, fixed = TRUE))
  expect_true(grepl("Linacre (2002)", evidence_021, fixed = TRUE))
  expect_true(grepl("Wright and Masters (1982)", evidence_021, fixed = TRUE))
  expect_true(grepl("not inter-rater agreement", combined, fixed = TRUE))
  expect_true(grepl("should not be treated as automatic\nvalidation success criteria", readme, fixed = TRUE) ||
                grepl("should not be treated as automatic validation success criteria",
                      gsub("\\s+", " ", readme), fixed = TRUE))
  expect_true(grepl("not standalone success criteria", precision_help, fixed = TRUE))
  expect_true(grepl("ZSTD standardization", gsub("\\s+", " ", reporting_vignette), fixed = TRUE))
})

test_that("0.2.1 NEWS keeps public changes separate from validation internals", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")

  news_lines <- readLines(file.path(pkg_root, "NEWS.md"), warn = FALSE)
  starts <- grep("^# mfrmr ", news_lines)
  current <- starts[1]
  next_heading <- starts[starts > current][1]
  if (is.na(next_heading)) next_heading <- length(news_lines) + 1L
  news_021 <- paste(news_lines[current:(next_heading - 1L)], collapse = "\n")

  expect_true(grepl("condition_review", news_021, fixed = TRUE))
  expect_true(grepl("slope_regime", news_021, fixed = TRUE))
  expect_true(grepl("extended sensitivity evidence", news_021, fixed = TRUE))
  expect_true(grepl("fit/separation operating characteristics", news_021, fixed = TRUE))
  expect_true(grepl("`reading_order`", news_021, fixed = TRUE))
  expect_false(grepl("inst/validation/release-evidence-map-0.2.1.md", news_021, fixed = TRUE))

  internal_terms <- c(
    "gpcm_high_dispersion_sparse",
    "summary.mfrmr_recovery_validation",
    "registered recovery-validation appendix roles",
    "recovery_condition_review",
    "recovery_uncertainty_review",
    "ScoreLevelsObserved",
    "ZeroScoreLevels",
    "MinScoreCount",
    "ADEMP method metadata",
    "Regression tests cover"
  )
  hits <- internal_terms[vapply(internal_terms, grepl, logical(1),
                                x = news_021, fixed = TRUE)]
  expect_identical(
    hits,
    character(0),
    info = "0.2.1 NEWS should summarize public-facing changes; implementation and fixture details belong in validation artifacts or tests."
  )
})

test_that("release evidence map is source-grounded and user-facing", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")

  path <- file.path(pkg_root, "inst", "validation", "release-evidence-map-0.2.1.md")
  checklist_path <- file.path(pkg_root, "inst", "validation", "release-evidence-checklist-0.2.1.csv")
  historical_path <- file.path(pkg_root, "inst", "validation", "release-evidence-map-0.2.0.md")
  external_path <- file.path(pkg_root, "inst", "validation", "external-parameter-recovery-simulation-0.2.0.md")
  external_helper_path <- file.path(pkg_root, "inst", "validation", "external-recovery-audit.R")
  expect_true(file.exists(path))
  expect_true(file.exists(checklist_path))
  expect_true(file.exists(historical_path))
  expect_true(file.exists(external_path))
  expect_true(file.exists(external_helper_path))
  doc <- paste(readLines(path, warn = FALSE), collapse = "\n")
  historical_doc <- paste(readLines(historical_path, warn = FALSE), collapse = "\n")
  external_doc <- paste(readLines(external_path, warn = FALSE), collapse = "\n")
  external_helper <- paste(readLines(external_helper_path, warn = FALSE), collapse = "\n")
  readme <- paste(readLines(file.path(pkg_root, "README.md"), warn = FALSE), collapse = "\n")
  checklist <- utils::read.csv(checklist_path, stringsAsFactors = FALSE)

  expect_true(grepl("release-evidence-map-0.2.1.md", readme, fixed = TRUE))
  expect_true(grepl("release-evidence-checklist-0.2.1.csv", readme, fixed = TRUE))
  expect_true(grepl("external-parameter-recovery-simulation-0.2.0.md", readme, fixed = TRUE))
  expect_true(grepl("Andrich (1978)", doc, fixed = TRUE))
  expect_true(grepl("Wright and Masters (1982)", doc, fixed = TRUE))
  expect_true(grepl("Muraki (1992)", doc, fixed = TRUE))
  expect_true(grepl("Morris, White, and Crowther (2019)", doc, fixed = TRUE))
  expect_true(grepl("Implementation boundary", doc, fixed = TRUE))
  expect_true(grepl("0.2.1 checks added", doc, fixed = TRUE))
  expect_true(grepl("Parameter_Recovery_Simulation", external_doc, fixed = TRUE))
  expect_true(grepl("mfrmr_review_external_recovery_simulation", external_helper, fixed = TRUE))
  expect_true(grepl("Review limits", external_doc, fixed = TRUE))
  expect_true(grepl("Decision rule", historical_doc, fixed = TRUE))
  expect_true(grepl("Scorecard template", historical_doc, fixed = TRUE))
  expect_true(grepl("Pre-release action plan for 0.2.0", historical_doc, fixed = TRUE))
  expect_true(grepl("Post-release roadmap", historical_doc, fixed = TRUE))
  expect_true(all(c("Domain", "Item", "SourceBasis", "PackageSurface",
                    "RequiredEvidence", "ReleaseDecision", "UserImplication",
                    "FollowUp") %in% names(checklist)))
  expect_true(nrow(checklist) >= 10)
  expect_true(any(checklist$ReleaseDecision == "blocker_if_failed"))
  expect_true(any(checklist$ReleaseDecision == "caveat_if_incomplete"))
  expect_true(any(checklist$ReleaseDecision == "roadmap_if_missing"))
  expect_true(all(checklist$ReleaseDecision %in% c(
    "blocker_if_failed", "caveat_if_incomplete", "roadmap_if_missing"
  )))
  expect_false(grepl("delete alias|can remove|removal decision", doc, ignore.case = TRUE))
})

test_that("current user guides use review spellings for migrated helpers", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")

  guide_files <- c(
    file.path(pkg_root, "README.md"),
    list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "inst", "cheatsheet"), pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE)
  )
  guide_files <- guide_files[file.exists(guide_files)]
  legacy_review_names <- c(
    "audit_mfrm_anchors(",
    "precision_audit_report(",
    "audit_conquest_overlap(",
    "facet_small_sample_audit(",
    "build_weighting_audit(",
    "reference_case_audit("
  )

  hits <- character(0)
  for (path in guide_files) {
    lines <- readLines(path, warn = FALSE)
    for (legacy_name in legacy_review_names) {
      idx <- grep(legacy_name, lines, fixed = TRUE)
      if (length(idx) > 0L) {
        rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path)
        hits <- c(hits, paste0(rel, ":", idx, ": ", legacy_name))
      }
    }
  }

  expect_identical(hits, character(0))
})

test_that("generated help pages are named after review helpers", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")

  old_pages <- file.path(pkg_root, "man", c(
    "audit_mfrm_anchors.Rd",
    "precision_audit_report.Rd",
    "audit_conquest_overlap.Rd",
    "facet_small_sample_audit.Rd",
    "build_weighting_audit.Rd",
    "reference_case_audit.Rd",
    "plot.mfrm_anchor_audit.Rd",
    "plot.mfrm_facet_sample_audit.Rd",
    "summary.mfrm_anchor_audit.Rd",
    "summary.mfrm_weighting_audit.Rd"
  ))
  new_pages <- file.path(pkg_root, "man", c(
    "review_mfrm_anchors.Rd",
    "precision_review_report.Rd",
    "review_conquest_overlap.Rd",
    "facet_small_sample_review.Rd",
    "build_weighting_review.Rd",
    "reference_case_review.Rd",
    "plot.mfrm_anchor_review.Rd",
    "plot.mfrm_facet_sample_review.Rd",
    "summary.mfrm_anchor_review.Rd",
    "summary.mfrm_weighting_review.Rd"
  ))

  expect_false(any(file.exists(old_pages)))
  expect_true(all(file.exists(new_pages)))
})

test_that("migrated review internals use canonical review labels in public routes", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source files are not available")

  source_files <- c(
    file.path(pkg_root, "R", "api-advanced.R"),
    file.path(pkg_root, "R", "api-methods.R"),
    file.path(pkg_root, "R", "api-reporting-checklist.R"),
    file.path(pkg_root, "R", "api-reports.R")
  )
  source_files <- source_files[file.exists(source_files)]
  stale_patterns <- c(
    "SourceTable = \"anchor_audit",
    "plot(anchor_audit, type =",
    ".weighting_audit_",
    "summarize_precision_audit_bundle",
    "summarize_conquest_overlap_audit_bundle",
    "precision_audit = summary_table_bundle_df(summary_obj$precision_audit)",
    "precision_audit = \"precision_audit\"",
    "Reference Audit Summary",
    "Hierarchical structure audit",
    "Complete-case omission audit"
  )

  hits <- character(0)
  for (path in source_files) {
    lines <- readLines(path, warn = FALSE)
    for (pattern in stale_patterns) {
      idx <- grep(pattern, lines, fixed = TRUE)
      if (length(idx) > 0L) {
        rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path)
        hits <- c(hits, paste0(rel, ":", idx, ": ", pattern))
      }
    }
  }

  expect_identical(hits, character(0))
})

test_that("remaining audit wording in public docs is limited to source-path headers", {
  pkg_root <- documentation_source_root()
  testthat::skip_if(is.na(pkg_root), "source documentation files are not available")

  doc_files <- c(
    file.path(pkg_root, "README.md"),
    list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "inst", "cheatsheet"), pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "man"), pattern = "\\.Rd$", full.names = TRUE, recursive = TRUE)
  )
  doc_files <- doc_files[file.exists(doc_files)]
  expect_true(length(doc_files) > 0)

  allowlist <- data.frame(
    reason = c(
      "roxygen source path header"
    ),
    pattern = c(
      "^% Please edit documentation in R/.*audit.*\\.R$"
    ),
    stringsAsFactors = FALSE
  )

  rel_path <- function(path) {
    sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path)
  }
  is_allowed <- function(line) {
    any(vapply(allowlist$pattern, grepl, logical(1), x = line, perl = TRUE, ignore.case = TRUE))
  }

  hits <- character(0)
  for (path in doc_files) {
    lines <- readLines(path, warn = FALSE)
    idx <- grep("\\baudit\\b|\\bAudit\\b|_audit|audit_", lines, perl = TRUE)
    if (length(idx) == 0L) next
    for (i in idx) {
      line <- lines[[i]]
      if (!is_allowed(line)) {
        hits <- c(hits, paste0(rel_path(path), ":", i, ": ", line))
      }
    }
  }

  expect_identical(
    hits,
    character(0),
    info = paste(
      "Any remaining audit wording in current public docs must be limited to",
      "generated source-path headers covered by the allowlist."
    )
  )
})

test_that("FACETS-facing docs use output-contract wording consistently", {
  pkg_root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  doc_files <- c(
    file.path(pkg_root, "README.md"),
    list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "man"), pattern = "\\.Rd$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "inst", "references"), pattern = "\\.(md|csv)$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(pkg_root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE, recursive = TRUE)
  )
  doc_files <- doc_files[file.exists(doc_files)]
  doc_files <- doc_files[basename(doc_files) != "test-documentation-terminology.R"]

  old_contract <- paste0("compatibility", "[- ]contract")
  old_spec <- paste0("compatibility", " specification")
  old_label <- paste0("Compatibility", " Output")

  rel_path <- function(path) {
    sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_root), "/?"), "", path)
  }
  hits <- character(0)
  for (path in doc_files) {
    lines <- readLines(path, warn = FALSE)
    idx <- grep(
      paste(c(old_contract, old_spec, old_label), collapse = "|"),
      lines,
      ignore.case = FALSE,
      perl = TRUE
    )
    if (length(idx) > 0L) {
      hits <- c(hits, paste0(rel_path(path), ":", idx, ": ", lines[idx]))
    }
  }

  expect_identical(
    hits,
    character(0),
    info = paste0(
      "FACETS-facing public docs should use output-contract wording, not ",
      "compatibility", "-contract wording."
    )
  )
})
