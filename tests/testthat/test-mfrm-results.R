test_that("mfrm_results builds a comprehensive object from a fitted model", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:8], , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 30
    )
  )

  res <- suppressWarnings(
    mfrm_results(
      fit,
      include = c("fit", "diagnostics", "tables", "precision", "reporting")
    )
  )
  expect_s3_class(res, "mfrm_results")
  expect_s3_class(res$diagnostics, "mfrm_diagnostics")
  expect_true(any(res$status$Section == "diagnostics" & res$status$Status == "ok"))
  expect_true("fit" %in% names(res$components))
  expect_true(length(res$tables) > 0)

  sx <- summary(res)
  expect_s3_class(sx, "summary.mfrm_results")
  expect_true(all(c("overview", "triage", "status", "plot_map", "next_actions") %in% names(sx)))
  expect_true(nrow(sx$triage) > 0)
  expect_true(all(c("Area", "Severity", "Signal", "Route", "Detail") %in% names(sx$triage)))
  expect_true(any(sx$triage$Area == "Diagnostics"))
  expect_true(nrow(sx$plot_map) > 0)
  expect_true(nrow(sx$next_actions) > 0)
  expect_true(any(sx$next_actions$Area %in% "Triage") || any(sx$triage$Severity %in% "ok"))
  expect_true(any(grepl("mfrm_results", sx$reproducible_code$Code)))

  bundle <- build_summary_table_bundle(res)
  expect_s3_class(bundle, "mfrm_summary_table_bundle")
  expect_identical(bundle$summary_class, "summary.mfrm_results")
  expect_true("triage" %in% names(bundle$tables))
  expect_true("next_actions" %in% names(bundle$tables))

  plt <- plot(res, type = "tables", draw = FALSE)
  expect_s3_class(plt, "mfrm_plot_data")

  report <- mfrm_report(res, style = "qc")
  expect_s3_class(report, "mfrm_report")
  printed_report <- capture.output(print(report))
  expect_true(any(grepl("Read order: summary(report)", printed_report, fixed = TRUE)))
  expect_true(any(grepl("Detailed tables are available in report$tables.", printed_report, fixed = TRUE)))
  expect_false(any(grepl("^Fit evidence summary$", printed_report)))
  expect_false(any(grepl("^Fit reporting templates$", printed_report)))
  report_summary <- summary(report)
  expect_s3_class(report_summary, "summary.mfrm_report")
  expect_true(all(c(
    "overview", "first_screen", "status_counts", "immediate_actions",
    "optional_sections", "claim_readiness", "report_gaps", "boundary_index",
    "routes"
  ) %in% names(report_summary)))
  expect_true(report_summary$overview$OverallStatus[1] %in% c(
    "ok", "review", "caveat", "request_if_needed", "unavailable"
  ))
  expect_true(any(report_summary$first_screen$Area == "Overall"))
  expect_true(any(report_summary$immediate_actions$Status %in% c("review", "caveat", "unavailable")))
  expect_true(any(report_summary$optional_sections$Status == "request_if_needed"))
  expect_true(any(report_summary$routes$Route == "report$first_screen"))
  report_summary_bundle <- build_summary_table_bundle(report)
  expect_s3_class(report_summary_bundle, "mfrm_summary_table_bundle")
  expect_identical(report_summary_bundle$summary_class, "summary.mfrm_report")
  expect_true(all(c("overview", "first_screen", "immediate_actions") %in% names(report_summary_bundle$tables)))
  expect_true(all(c(
    "first_screen", "report_index", "section_plan", "claim_readiness", "report_gaps",
    "template_index",
    "fit_criteria", "fit_evidence_summary", "fit_threshold_sensitivity",
    "fit_reporting_templates", "zstd_conventions", "fit_decision_policy",
    "precision_evidence_summary", "precision_basis", "precision_reporting_templates",
    "bias_evidence_summary", "bias_reporting_templates",
    "misfit_evidence_summary", "misfit_reporting_templates",
    "linking_evidence_summary", "linking_reporting_templates",
    "fit_df_sensitivity_summary", "fit_df_sensitive_rows",
    "evidence_boundary", "action_items"
  ) %in% names(report$tables)))
  expect_true(all(c(
    "Area", "Status", "Readiness", "MainIssue", "NextAction", "PrimaryRoute"
  ) %in% names(report$first_screen)))
  expect_true(report$first_screen$Area[1] == "Overall")
  expect_true(any(report$first_screen$Area == "Bias / DFF" &
                    report$first_screen$Status == "request_if_needed"))
  expect_true(any(grepl("report\\$report_index", report$first_screen$PrimaryRoute)))
  expect_true(all(c(
    "Area", "EvidenceStatus", "Readiness", "ReviewSignalCount",
    "PrimaryTable", "TemplateTable", "EvidenceRoute", "TemplateRoute",
    "PlotRoute", "ExportRoute", "IncludePreset", "Boundary"
  ) %in% names(report$report_index)))
  expect_true(all(c(
    "Fit", "Precision", "Bias / DFF", "Misfit / pathway", "Linking / anchors"
  ) %in% report$report_index$Area))
  expect_true(any(report$report_index$Area == "Fit" &
                    grepl("type = 'qc'|type = \"qc\"", report$report_index$PlotRoute)))
  expect_true(all(grepl("export_mfrm_results", report$report_index$ExportRoute, fixed = TRUE)))
  expect_true(any(grepl("include = \"bias\"", report$report_index$IncludePreset, fixed = TRUE)))
  expect_true(any(report$report_index$Readiness %in% c("ready", "review", "request_if_needed")))
  expect_true(all(c(
    "Area", "TemplateTable", "TemplateRow", "Topic", "BoundaryType",
    "ClaimStrength", "RecommendedUse", "EvidenceTable", "EvidenceRoute"
  ) %in% names(report$template_index)))
  expect_true(any(report$template_index$TemplateTable == "fit_reporting_templates"))
  expect_true(any(report$template_index$Area == "Bias / DFF" &
                    report$template_index$ClaimStrength == "not_supported_without_followup"))
  expect_true(any(report$template_index$BoundaryType == "fit_not_validity"))
  expect_true(any(report$sections$Section == "Fit, separation, and precision"))
  expect_true(any(report$claim_readiness$Claim == "Fit and precision evidence"))
  expect_true(any(report$claim_readiness$Readiness %in% c("ready", "write_with_caveat", "needs_requested_section")))
  expect_true(any(report$report_gaps$Section %in% "Bias screening"))
  expect_true(all(c("Priority", "GapType", "RecommendedAction") %in% names(report$report_gaps)))
  expect_true(any(report$fit_criteria$Profile == "active"))
  expect_true(any(report$fit_criteria$Profile == "linacre_productive"))
  expect_true(all(c(
    "Rows", "UnderfitRows", "FitDfMethod", "FacetsCompanionAvailable"
  ) %in% names(report$fit_evidence_summary)))
  expect_true(report$fit_evidence_summary$Status[1] %in% "available")
  expect_true(report$fit_evidence_summary$FacetsCompanionAvailable[1])
  expect_true(all(c("Profile", "AnyFlagRate", "ReportBoundary") %in% names(report$fit_threshold_sensitivity)))
  expect_true(any(report$fit_threshold_sensitivity$Profile == "active"))
  expect_true(all(c("Topic", "Template", "Caveat", "Route") %in% names(report$fit_reporting_templates)))
  expect_true(any(report$fit_reporting_templates$Topic == "DF/ZSTD sensitivity wording"))
  expect_true(any(grepl("MnSq", report$fit_reporting_templates$Caveat, fixed = TRUE)))
  expect_true(all(c(
    "PrecisionTier", "SupportsFormalInference", "ReliabilityRows",
    "MinSeparation", "MaxReliability", "Boundary"
  ) %in% names(report$precision_evidence_summary)))
  expect_true(report$precision_evidence_summary$Status[1] %in% "available")
  expect_true(any(grepl("Separation", report$precision_basis$Topic, fixed = TRUE)))
  expect_true(all(c("Topic", "Template", "Caveat", "Route") %in% names(report$precision_reporting_templates)))
  expect_true(any(report$precision_reporting_templates$Topic == "Reliability wording"))
  expect_true(any(grepl("inter-rater agreement", report$precision_reporting_templates$Caveat, fixed = TRUE)))
  expect_true(all(c("Status", "Rows", "InteractionStatus", "Boundary") %in% names(report$bias_evidence_summary)))
  expect_true(report$bias_evidence_summary$Status[1] %in% "not_requested")
  expect_true(any(report$bias_reporting_templates$Topic == "Bias/DFF evidence not requested"))
  expect_true(all(c(
    "Status", "UnexpectedRows", "DisplacementRows", "PathwayAvailable", "Boundary"
  ) %in% names(report$misfit_evidence_summary)))
  expect_true(report$misfit_evidence_summary$Status[1] %in% "not_requested")
  expect_true(any(report$misfit_reporting_templates$Topic == "Misfit/pathway evidence not requested"))
  expect_true(all(c(
    "Status", "ReviewStatus", "AnchorReviewAvailable", "DriftAvailable",
    "ChainAvailable", "Boundary"
  ) %in% names(report$linking_evidence_summary)))
  expect_true(report$linking_evidence_summary$Status[1] %in% "not_requested")
  expect_true(any(report$linking_reporting_templates$Topic == "Linking evidence not requested"))
  expect_true(any(report$zstd_conventions$Convention == "FACETS-style df"))
  expect_true(any(grepl("df < 1", report$zstd_conventions$PackageConstraint, fixed = TRUE)))
  expect_true(any(report$fit_decision_policy$Rule == "Read MnSq before ZSTD"))
  expect_true(all(c(
    "ComparedRows", "DfSensitiveRows", "FitDfMethod", "ReportBoundary"
  ) %in% names(report$fit_df_sensitivity_summary)))
  expect_true("DfSensitivityStatus" %in% names(report$fit_df_sensitive_rows))
  expect_true(any(grepl("Fit, separation", report$evidence_boundary$EvidenceSource, fixed = TRUE)))
  template_contract <- c(
    "Topic", "Template", "EvidenceUsed", "EvidenceTable", "EvidenceRoute",
    "BoundaryType", "ClaimStrength", "RecommendedUse", "Caveat", "Route"
  )
  template_tables <- list(
    fit = report$fit_reporting_templates,
    precision = report$precision_reporting_templates,
    bias = report$bias_reporting_templates,
    misfit = report$misfit_reporting_templates,
    linking = report$linking_reporting_templates
  )
  expect_true(all(vapply(
    template_tables,
    function(tbl) all(template_contract %in% names(tbl)),
    logical(1)
  )))
  expect_true(any(report$fit_reporting_templates$BoundaryType == "fit_not_validity"))
  expect_true(any(report$precision_reporting_templates$BoundaryType == "precision_not_agreement"))
  expect_true(any(report$bias_reporting_templates$BoundaryType == "screen_not_fairness_decision"))
  expect_true(any(report$misfit_reporting_templates$BoundaryType == "misfit_not_exclusion_rule"))
  expect_true(any(report$linking_reporting_templates$BoundaryType == "anchor_not_drift_absence"))
  expect_true(any(report$bias_reporting_templates$ClaimStrength == "not_supported_without_followup"))
  expect_true(any(report$fit_reporting_templates$RecommendedUse == "methods_or_appendix_caveat"))

  markdown <- mfrm_report(res, output = "markdown")
  expect_type(markdown, "character")
  expect_length(markdown, 1L)
  expect_match(markdown, "Evidence Boundary", fixed = TRUE)
  expect_match(markdown, "First Screen", fixed = TRUE)
  expect_match(markdown, "Report Index", fixed = TRUE)
  expect_match(markdown, "Template Index", fixed = TRUE)
  expect_match(markdown, "Claim Readiness", fixed = TRUE)
  expect_match(markdown, "Report Gaps", fixed = TRUE)
  expect_match(markdown, "Fit Evidence Summary", fixed = TRUE)
  expect_match(markdown, "Fit Reporting Templates", fixed = TRUE)
  expect_match(markdown, "Precision Evidence Summary", fixed = TRUE)
  expect_match(markdown, "Precision Reporting Templates", fixed = TRUE)
  expect_match(markdown, "Bias Evidence Summary", fixed = TRUE)
  expect_match(markdown, "Bias Reporting Templates", fixed = TRUE)
  expect_match(markdown, "Misfit Evidence Summary", fixed = TRUE)
  expect_match(markdown, "Misfit Reporting Templates", fixed = TRUE)
  expect_match(markdown, "Linking Evidence Summary", fixed = TRUE)
  expect_match(markdown, "Linking Reporting Templates", fixed = TRUE)
  expect_match(markdown, "ZSTD Conventions", fixed = TRUE)
  expect_match(markdown, "Fit DF Sensitivity", fixed = TRUE)

  html <- mfrm_report(res, style = "validation", output = "html")
  expect_s3_class(html, "mfrm_report_html")
  expect_true(file.exists(html$path))
  expect_match(html$html, "Reader guidance", fixed = TRUE)
  expect_match(html$html, "report_summary_overview", fixed = TRUE)
  expect_match(html$html, "report_summary_first_screen", fixed = TRUE)
  expect_match(html$html, "Report Markdown", fixed = TRUE)
  expect_lt(
    regexpr("report_summary_overview", html$html, fixed = TRUE)[1],
    regexpr("Report Markdown", html$html, fixed = TRUE)[1]
  )
  unlink(html$path)

  tables <- mfrm_report(res, style = "reviewer", output = "tables")
  expect_type(tables, "list")
  expect_true("first_screen" %in% names(tables))
  expect_true("report_index" %in% names(tables))
  expect_true("template_index" %in% names(tables))
  apa_report <- mfrm_report(res, style = "apa")
  expect_s3_class(apa_report, "mfrm_report")
  expect_true(any(apa_report$fit_reporting_templates$Audience == "APA manuscript"))
  expect_true(any(apa_report$precision_reporting_templates$Audience == "APA manuscript"))
  expect_true(any(apa_report$misfit_reporting_templates$Audience == "APA manuscript"))
  expect_true(any(apa_report$linking_reporting_templates$Audience == "APA manuscript"))
  expect_error(
    mfrm_report(fit),
    "Call `mfrm_results\\(\\)` first"
  )
})

test_that("launch_mfrmr_viewer is a viewer over mfrm_results", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:6], , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 30
    )
  )
  res <- suppressWarnings(
    mfrm_results(fit, include = c("fit", "diagnostics", "tables"))
  )

  payload <- mfrmr:::mfrm_results_viewer_payload(res)
  expect_true(all(c(
    "summary", "tables", "table_names", "qc_tables", "report_tables",
    "report_text", "bias_table", "bias_guidance", "unexpected_table",
    "tab_status", "plot_choices", "replay_code"
  ) %in% names(payload)))
  expect_s3_class(payload$summary, "summary.mfrm_results")
  expect_true(length(payload$table_names) > 0L)
  expect_true(length(payload$qc_table_names) > 0L)
  expect_true(nrow(payload$bias_guidance) > 0L)
  expect_true(nrow(payload$unexpected_table) > 0L)
  expect_true(all(c("qc", "report", "bias", "misfit") %in% names(payload$tab_status)))
  expect_true(any(payload$tab_status$report$Name == "apa_outputs" &
                    payload$tab_status$report$Status == "not_requested"))
  expect_true(any(payload$tab_status$bias$Name == "bias_screen" &
                    payload$tab_status$bias$Status == "not_requested"))
  expect_true(any(payload$tab_status$bias$Name == "bias_screen" &
                    grepl("include = \"bias\"", payload$tab_status$bias$Detail, fixed = TRUE)))
  expect_true(any(payload$tab_status$misfit$Name == "misfit_review" &
                    payload$tab_status$misfit$Status == "not_requested"))
  expect_true(any(payload$tab_status$misfit$Name == "misfit_review" &
                    grepl("include = \"misfit_review\"", payload$tab_status$misfit$Detail, fixed = TRUE)))
  expect_true(length(payload$plot_choices) > 0L)
  expect_match(payload$replay_code, "mfrm_results", fixed = TRUE)

  expect_error(
    launch_mfrmr_viewer(fit, return_app = TRUE),
    "Call `mfrm_results\\(\\)` first"
  )

  if (requireNamespace("shiny", quietly = TRUE)) {
    app <- launch_mfrmr_viewer(res, return_app = TRUE)
    expect_true(inherits(app, "shiny.appobj"))
  } else {
    expect_error(
      launch_mfrmr_viewer(res, return_app = TRUE),
      "requires the optional package `shiny`",
      fixed = TRUE
    )
  }
})

test_that("export_mfrm_results writes lightweight result downloads", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:6], , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 30
    )
  )
  res <- suppressWarnings(
    mfrm_results(fit, include = c("fit", "diagnostics", "tables"))
  )

  out_dir <- file.path(tempdir(), paste0("mfrmr_results_export_", sample.int(1e6, 1)))
  exported <- export_mfrm_results(
    res,
    output_dir = out_dir,
    prefix = "results export",
    overwrite = TRUE
  )

  expect_s3_class(exported, "mfrm_results_export")
  expect_true(dir.exists(out_dir))
  expect_true(nrow(exported$written_files) > 0L)
  expect_true(all(file.exists(exported$written_files$Path)))
  expect_true(any(exported$written_files$Component == "summary_overview"))
  expect_true(any(exported$written_files$Component == "results_html"))
  expect_true(any(exported$written_files$Component == "results_rds"))
  expect_true(any(exported$written_files$Component == "replay_code"))
  expect_true(any(exported$written_files$Component == "written_files"))
  expect_true(any(grepl("_table_", basename(exported$written_files$Path), fixed = TRUE)))
  expect_true(exported$summary$CsvWritten[1] >= 3L)

  manifest_path <- exported$written_files$Path[exported$written_files$Component == "written_files"][1]
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  expect_true(all(c("Component", "Format", "Path", "Note") %in% names(manifest)))
  expect_true(any(manifest$Component == "results_rds"))

  rds_path <- exported$written_files$Path[exported$written_files$Component == "results_rds"][1]
  expect_s3_class(readRDS(rds_path), "mfrm_results")

  expect_error(
    export_mfrm_results(res, output_dir = out_dir, prefix = "results export"),
    "File already exists"
  )

  plot_dir <- file.path(tempdir(), paste0("mfrmr_results_export_plots_", sample.int(1e6, 1)))
  plot_export <- export_mfrm_results(
    res,
    output_dir = plot_dir,
    prefix = "plots",
    include = c("summary", "plots", "manifest"),
    overwrite = TRUE
  )
  expect_true(any(plot_export$written_files$Format == "png") || nrow(plot_export$plot_errors) > 0L)

  report_dir <- file.path(tempdir(), paste0("mfrmr_results_export_report_", sample.int(1e6, 1)))
  report_export <- export_mfrm_results(
    res,
    output_dir = report_dir,
    prefix = "report",
    include = c("report", "manifest"),
    overwrite = TRUE
  )
  expect_s3_class(report_export, "mfrm_results_export")
  expect_true(any(report_export$written_files$Component == "report_first_screen"))
  expect_true(any(report_export$written_files$Component == "report_report_index"))
  expect_true(any(report_export$written_files$Component == "report_template_index"))
  expect_true(any(report_export$written_files$Component == "report_markdown"))
  expect_true(any(report_export$written_files$Component == "report_html"))
  expect_true(all(file.exists(report_export$written_files$Path)))
  first_screen_path <- report_export$written_files$Path[
    report_export$written_files$Component == "report_first_screen"
  ][1]
  first_screen <- utils::read.csv(first_screen_path, stringsAsFactors = FALSE)
  expect_true(all(c("Area", "Status", "Readiness", "PrimaryRoute") %in% names(first_screen)))
  expect_true(first_screen$Area[1] == "Overall")
  report_index_path <- report_export$written_files$Path[
    report_export$written_files$Component == "report_report_index"
  ][1]
  report_index <- utils::read.csv(report_index_path, stringsAsFactors = FALSE)
  expect_true(all(c(
    "Area", "EvidenceStatus", "Readiness", "PrimaryTable",
    "EvidenceRoute", "TemplateRoute", "PlotRoute", "ExportRoute", "IncludePreset"
  ) %in% names(report_index)))
  expect_true(any(report_index$Area == "Fit"))
  expect_true(any(grepl("export_mfrm_results", report_index$ExportRoute, fixed = TRUE)))
  template_index_path <- report_export$written_files$Path[
    report_export$written_files$Component == "report_template_index"
  ][1]
  template_index <- utils::read.csv(template_index_path, stringsAsFactors = FALSE)
  expect_true(all(c("TemplateTable", "BoundaryType", "ClaimStrength") %in% names(template_index)))
  expect_true(any(template_index$TemplateTable == "fit_reporting_templates"))
  report_md_path <- report_export$written_files$Path[
    report_export$written_files$Component == "report_markdown"
  ][1]
  expect_true(any(grepl("Report Index", readLines(report_md_path, warn = FALSE), fixed = TRUE)))
})

test_that("mfrm_results output modes and standard data-frame route work", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:6], , drop = FALSE]

  summary_out <- suppressWarnings(
    mfrm_results(toy_small, include = c("fit", "diagnostics"), output = "summary")
  )
  expect_s3_class(summary_out, "summary.mfrm_results")
  expect_identical(summary_out$overview$InputMode[1], "data.frame")
  expect_true(nrow(summary_out$mapping) > 0)
  expect_true(any(grepl("fit_mfrm", summary_out$reproducible_code$Code)))

  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 30
    )
  )
  tables <- suppressWarnings(
    mfrm_results(fit, include = c("fit", "diagnostics", "tables"), output = "tables")
  )
  expect_type(tables, "list")
  expect_true(length(tables) > 0)

  html <- suppressWarnings(mfrm_results(fit, include = "fit", output = "html"))
  expect_s3_class(html, "mfrm_results_html")
  expect_true(file.exists(html$path))
  unlink(html$path)
})

test_that("mfrm_results purpose presets resolve to documented sections", {
  presets <- mfrmr:::mfrm_results_include_preset_table()
  expect_s3_class(presets, "data.frame")
  expect_true(all(c("Preset", "Sections") %in% names(presets)))
  expect_true(all(c("standard", "publication", "validation", "facets", "linking", "network", "gpcm_review", "all") %in% presets$Preset))
  expect_true(all(c("bias", "bias_review", "misfit", "misfit_review") %in% presets$Preset))

  publication <- mfrmr:::mfrm_results_resolve_include("publication")
  expect_true(all(c("fit", "diagnostics", "reporting", "apa") %in% publication))
  expect_identical(attr(publication, "requested"), "publication")
  expect_identical(attr(publication, "presets"), "publication")

  validation <- mfrmr:::mfrm_results_resolve_include("validation")
  expect_true(all(c("fit", "diagnostics", "precision", "facets_fit") %in% validation))

  network <- mfrmr:::mfrm_results_resolve_include("network")
  expect_true(all(c("fit", "diagnostics", "network") %in% network))

  bias <- mfrmr:::mfrm_results_resolve_include("bias")
  expect_true(all(c("fit", "diagnostics", "bias") %in% bias))

  misfit <- mfrmr:::mfrm_results_resolve_include("misfit_review")
  expect_true(all(c("fit", "diagnostics", "misfit") %in% misfit))

  linking <- mfrmr:::mfrm_results_resolve_include("linking")
  expect_true(all(c("fit", "diagnostics", "linking") %in% linking))

  gpcm <- mfrmr:::mfrm_results_resolve_include("gpcm")
  expect_true(all(c("fit", "diagnostics", "precision", "reporting") %in% gpcm))
  expect_identical(attr(gpcm, "presets"), "gpcm_review")

  rt <- mfrmr:::mfrm_results_resolve_include("timing")
  expect_true("response_time" %in% rt)

  all_sections <- mfrmr:::mfrm_results_resolve_include("all")
  expect_true("response_time" %in% all_sections)

  expect_identical(
    mfrmr:::mfrm_results_maybe_response_time_column(data.frame(ResponseTimeMs = 1)),
    "ResponseTimeMs"
  )
})

test_that("mfrm_results can attach explicit response-time review", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:6], , drop = FALSE]
  toy_small$ResponseTime <- 8 + seq_len(nrow(toy_small)) %% 5 +
    as.numeric(toy_small$Score)
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 20
    )
  )

  res <- suppressWarnings(
    mfrm_results(
      fit,
      include = c("fit", "response_time"),
      response_time = "ResponseTime",
      response_time_data = toy_small
    )
  )
  expect_s3_class(res, "mfrm_results")
  expect_s3_class(res$components$response_time_review, "mfrm_response_time_review")
  expect_true(any(res$status$Section == "response_time_review" &
                    res$status$Status == "ok"))
  expect_true(any(grepl("^response_time_review_", names(res$tables))))
  expect_true(any(res$plot_map$Type == "response_time" &
                    res$plot_map$Available))
  sx <- summary(res)
  expect_true(any(sx$triage$Area == "Response-time QC"))
  expect_true(any(sx$next_actions$Area == "Response-time QC"))
  expect_true(any(grepl("response_time_data", sx$reproducible_code$Code,
                        fixed = TRUE)))

  p <- plot(res, type = "response_time", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(is.data.frame(plot_data(p, component = "table")))

  report <- mfrm_report(res, style = "qc")
  expect_true(any(report$sections$Section == "Response-time QC"))
  expect_true(any(grepl("speed parameter", report$evidence_boundary$DoNotUseAs,
                        fixed = TRUE)))

  payload <- mfrmr:::mfrm_results_viewer_payload(res)
  expect_true("response_time" %in% names(payload$tab_status))
  expect_true(length(payload$response_time_table_names) > 0L)
  expect_true(any(payload$tab_status$response_time$Status == "ok"))
  expect_true(any(grepl("joint speed-accuracy model",
                        payload$response_time_boundary$Detail,
                        fixed = TRUE)))
  expect_true(any(grepl("modified logit estimate",
                        payload$response_time_boundary$Detail,
                        fixed = TRUE)))

  out_dir <- file.path(tempdir(), paste0("mfrmr_rt_export_", sample.int(1e6, 1)))
  exported <- export_mfrm_results(
    res,
    output_dir = out_dir,
    prefix = "rt",
    include = c("summary", "tables", "manifest"),
    overwrite = TRUE
  )
  expect_true(any(grepl("response_time_review", exported$written_files$Component,
                        fixed = TRUE)))

  res_auto <- suppressWarnings(
    mfrm_results(
      fit,
      include = c("fit", "response_time"),
      response_time_data = toy_small
    )
  )
  expect_s3_class(res_auto$components$response_time_review, "mfrm_response_time_review")
  expect_identical(
    res_auto$components$response_time_review$config$results_route$TimeColumn,
    "ResponseTime"
  )
})

test_that("mfrm_results response-time route is unavailable without timing data", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:5], , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 10
    )
  )

  res <- suppressWarnings(mfrm_results(fit, include = "response_time"))
  expect_false("response_time_review" %in% names(res$components))
  expect_true(any(res$status$Section == "response_time_review" &
                    res$status$Status == "not_available"))
  expect_false(res$plot_map$Available[res$plot_map$Type == "response_time"])
  expect_error(
    plot(res, type = "response_time", draw = FALSE),
    "not available"
  )
})

test_that("mfrm_results data-frame route treats ResponseTime as metadata", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:5], , drop = FALSE]
  toy_small$ResponseTime <- 10 + seq_len(nrow(toy_small)) %% 4

  sx <- suppressWarnings(
    mfrm_results(
      toy_small,
      include = c("fit", "response_time"),
      output = "summary"
    )
  )
  expect_s3_class(sx, "summary.mfrm_results")
  facet_row <- sx$mapping$Value[sx$mapping$Key == "Facets"]
  expect_false(grepl("ResponseTime", facet_row, fixed = TRUE))
  expect_true(any(sx$status$Section == "response_time_review" &
                    sx$status$Status == "ok"))
  expect_true(any(grepl("response_time = \"ResponseTime\"",
                        sx$reproducible_code$Code,
                        fixed = TRUE)))
})

test_that("mfrm_results bias and misfit presets expose review surfaces", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:7], , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 30
    )
  )

  res_bias <- suppressWarnings(mfrm_results(fit, include = "bias"))
  expect_s3_class(res_bias, "mfrm_results")
  expect_true("bias_screen" %in% names(res_bias$components))
  expect_true("bias_screen_guidance" %in% names(res_bias$tables))
  expect_true(any(summary(res_bias)$triage$Area == "Bias screening"))
  bias_report <- mfrm_report(res_bias, style = "reviewer")
  expect_s3_class(bias_report, "mfrm_report")
  expect_true(bias_report$bias_evidence_summary$Status[1] %in% "available")
  expect_true(all(c(
    "Rows", "ResidualTScreenPositiveRows", "ChiSqScreenPositiveRows",
    "ExplicitInteractionSelected", "Boundary"
  ) %in% names(bias_report$bias_evidence_summary)))
  expect_true(any(bias_report$bias_reporting_templates$Topic == "Interaction contrast wording"))
  expect_true(any(grepl("fairness", bias_report$bias_reporting_templates$Caveat, ignore.case = TRUE)))
  expect_true(any(bias_report$bias_reporting_templates$BoundaryType == "screen_not_fairness_decision"))
  expect_true(any(bias_report$bias_reporting_templates$ClaimStrength == "not_supported_without_followup"))

  res_misfit <- suppressWarnings(mfrm_results(fit, include = "misfit_review"))
  expect_s3_class(res_misfit, "mfrm_results")
  expect_true("misfit_review" %in% names(res_misfit$components))
  expect_true(any(grepl("^misfit_review_", names(res_misfit$tables))))
  expect_true(any(summary(res_misfit)$triage$Area == "Pathway / misfit"))
  expect_true(any(res_misfit$plot_map$Type == "pathway" & res_misfit$plot_map$Available))
  expect_s3_class(plot(res_misfit, type = "pathway", draw = FALSE), "mfrm_plot_data")
  misfit_report <- mfrm_report(res_misfit, style = "reviewer")
  expect_s3_class(misfit_report, "mfrm_report")
  expect_true(any(misfit_report$report_index$Area == "Misfit / pathway" &
                    grepl("pathway", misfit_report$report_index$PlotRoute, fixed = TRUE)))
  expect_true(misfit_report$misfit_evidence_summary$Status[1] %in% "available")
  expect_true(all(c(
    "UnexpectedRows", "DisplacementRows", "PathwayAvailable", "Boundary"
  ) %in% names(misfit_report$misfit_evidence_summary)))
  expect_true(any(misfit_report$misfit_reporting_templates$Topic == "Case-review wording"))
  expect_true(any(grepl("exclusion", misfit_report$misfit_reporting_templates$Caveat, ignore.case = TRUE)))
  expect_true(any(misfit_report$misfit_reporting_templates$BoundaryType == "misfit_not_exclusion_rule"))
  expect_true(any(misfit_report$misfit_reporting_templates$EvidenceRoute == "report$misfit_evidence_summary"))
})

test_that("mfrm_results linking preset exposes anchor readiness", {
  toy <- load_mfrmr_data("example_core")
  toy_small <- toy[toy$Person %in% unique(toy$Person)[1:7], , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy_small,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      maxit = 30
    )
  )

  res <- suppressWarnings(mfrm_results(fit, include = "linking"))
  expect_s3_class(res, "mfrm_results")
  expect_s3_class(res$components$linking_review, "mfrm_linking_review")
  expect_true(any(grepl("^linking_review_", names(res$tables))))
  expect_true("linking_review_first_screen_guidance" %in% names(res$tables))
  expect_true(any(summary(res)$triage$Area == "Linking / anchors"))
  expect_true(any(summary(res)$next_actions$Area == "Linking / anchors"))
  expect_true(any(res$plot_map$Type == "anchors" & res$plot_map$Available))
  expect_s3_class(plot(res, type = "anchors", draw = FALSE), "mfrm_plot_data")
  linking_report <- mfrm_report(res, style = "reviewer")
  expect_s3_class(linking_report, "mfrm_report")
  expect_true(any(linking_report$report_index$Area == "Linking / anchors" &
                    grepl("anchors", linking_report$report_index$PlotRoute, fixed = TRUE)))
  expect_true(linking_report$linking_evidence_summary$Status[1] %in% "available")
  expect_true(all(c(
    "ReviewStatus", "TopRiskRows", "AnchoredLevels", "DriftReviewStatus",
    "EquatingChainStatus", "Boundary"
  ) %in% names(linking_report$linking_evidence_summary)))
  expect_true(any(linking_report$linking_reporting_templates$Topic == "Drift wording"))
  expect_true(any(grepl("single fitted object", linking_report$linking_reporting_templates$Caveat, ignore.case = TRUE)))
  expect_true(any(linking_report$linking_reporting_templates$BoundaryType == "anchor_not_drift_absence"))
  expect_true(any(linking_report$linking_reporting_templates$ClaimStrength == "not_supported_without_followup"))
})

test_that("mfrm_results_interactive is opt-in only", {
  testthat::skip_if(interactive(), "non-interactive guard is tested only under automated runs")
  expect_error(
    mfrm_results_interactive(data.frame(Person = "P1", Score = 1, Rater = "R1")),
    "interactive session"
  )
})
