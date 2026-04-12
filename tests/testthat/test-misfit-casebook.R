test_that("build_misfit_casebook returns a structured review bundle", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 9,
    maxit = 20
  )
  diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")

  casebook <- build_misfit_casebook(fit, diagnostics = diag, top_n = 8)

  expect_s3_class(casebook, "mfrm_misfit_casebook")
  expect_s3_class(casebook, "mfrm_bundle")
  expect_true(all(c(
    "overview", "status", "top_cases", "source_summary",
    "case_rollup", "group_view_index", "group_views",
    "plot_map", "reporting_map", "support_status", "source_support", "notes", "settings"
  ) %in% names(casebook)))
  expect_true(all(c(
    "CaseID", "CaseType", "SourceFamily", "SourceTable", "SourceRowKey",
    "AdministrationID", "WaveID", "PrimaryUnit", "PrimaryUnitType", "Magnitude", "ReviewPriority",
    "WithinSourceRank", "SupportBasis", "InterpretationTier",
    "PrimaryPlotRoute", "SupportStatus"
  ) %in% names(casebook$top_cases)))
  expect_true(all(c(
    "AdministrationID", "WaveID", "RollupType", "RollupKey", "Cases",
    "MaxPriority", "TopCaseID"
  ) %in% names(casebook$case_rollup)))
  expect_true(all(c("View", "Rows", "Description") %in% names(casebook$group_view_index)))
  expect_true(is.list(casebook$group_views))
  expect_true(all(c(
    "by_person", "by_facet_level", "by_facet_pair",
    "by_source_family", "by_facet", "by_administration", "by_wave", "facet_views"
  ) %in% names(casebook$group_views)))
  expect_true(any(casebook$support_status$Scope == "RSM / PCM"))
  expect_true(all(c("SourceFamily", "Available", "SupportBasis", "Status", "Note") %in%
                    names(casebook$source_support)))
})

test_that("summary methods for build_misfit_casebook expose a front-door summary", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 9,
    maxit = 20
  )
  diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
  unexpected <- unexpected_response_table(fit, diagnostics = diag, abs_z_min = 1.5, prob_max = 0.4, top_n = 8)
  displacement <- displacement_table(fit, diagnostics = diag, anchored_only = FALSE, top_n = 8)

  casebook <- build_misfit_casebook(
    fit,
    diagnostics = diag,
    unexpected = unexpected,
    displacement = displacement,
    top_n = 6
  )
  sx <- summary(casebook, top_n = 4)

  expect_s3_class(sx, "summary.mfrm_misfit_casebook")
  expect_true(all(c(
    "overview", "status", "key_warnings", "next_actions",
    "top_cases", "case_rollup", "group_view_index", "group_views",
    "source_summary", "source_support", "plot_routes", "plot_map", "reporting_map", "support_status"
  ) %in% names(sx)))
  expect_lte(nrow(sx$top_cases), 4)
  expect_output(print(sx), "Plot Follow-up")
})

test_that("build_misfit_casebook records administration and wave provenance", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 9,
    maxit = 20
  )
  diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")

  casebook <- build_misfit_casebook(
    fit,
    diagnostics = diag,
    administration_id = "FormA",
    wave_id = "Wave1",
    top_n = 6
  )

  expect_equal(casebook$overview$AdministrationID[[1]], "FormA")
  expect_equal(casebook$overview$WaveID[[1]], "Wave1")
  expect_equal(casebook$settings$administration_id, "FormA")
  expect_equal(casebook$settings$wave_id, "Wave1")
  if (nrow(casebook$top_cases) > 0) {
    expect_true(all(casebook$top_cases$AdministrationID == "FormA"))
    expect_true(all(casebook$top_cases$WaveID == "Wave1"))
  }
  if (nrow(casebook$case_rollup) > 0) {
    expect_true(all(casebook$case_rollup$AdministrationID == "FormA"))
    expect_true(all(casebook$case_rollup$WaveID == "Wave1"))
  }
  if (nrow(casebook$group_views$by_administration) > 0) {
    expect_true(all(casebook$group_views$by_administration$AdministrationID == "FormA"))
  }
  if (nrow(casebook$group_views$by_wave) > 0) {
    expect_true(all(casebook$group_views$by_wave$WaveID == "Wave1"))
  }
})

test_that("build_misfit_casebook can return a no-flagged-cases status", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 9,
    maxit = 20
  )
  diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
  unexpected <- unexpected_response_table(fit, diagnostics = diag, abs_z_min = 99, prob_max = 1e-12, top_n = 5)
  displacement <- displacement_table(fit, diagnostics = diag, abs_displacement_warn = 99, abs_t_warn = 99, top_n = 5)

  casebook <- build_misfit_casebook(
    fit,
    diagnostics = diag,
    unexpected = unexpected,
    displacement = displacement,
    top_n = 5
  )

  expect_equal(casebook$overview$ReviewStatus[[1]], "no_flagged_cases")
  expect_equal(nrow(casebook$top_cases), 0)
})

test_that("build_misfit_casebook marks bounded GPCM as supported with caveat", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:10]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    slope_facet = "Criterion",
    step_facet = "Criterion",
    quad_points = 5,
    maxit = 20
  ))

  casebook <- build_misfit_casebook(fit, top_n = 6)
  gpcm_row <- casebook$support_status[casebook$support_status$Scope == "bounded GPCM", , drop = FALSE]
  gpcm_sources <- casebook$source_support

  expect_equal(gpcm_row$Status[[1]], "supported_with_caveat")
  expect_true(any(casebook$top_cases$SupportStatus == "supported_with_caveat") || nrow(casebook$top_cases) == 0)
  expect_true(all(gpcm_sources$Status %in% c("supported_with_caveat", "deferred")))
})
