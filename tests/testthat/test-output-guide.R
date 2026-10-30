test_that("rater-feedback guidance separates uncertainty, misfit and known-truth accuracy", {
  feedback <- mfrmr_output_guide("feedback")
  expect_true(all(feedback$Scope == "feedback"))
  expect_true(all(feedback$Question %in% mfrmr_output_guide()$Question))
  expect_identical(names(feedback), names(mfrmr_output_guide("public")))
  expect_identical(feedback$DecisionBoundary, feedback$Notes)
  # Recommended calls must resolve; the guide must not invent a wrapper/API.
  functions <- unique(unlist(regmatches(feedback$MainFunction,
    gregexpr("[A-Za-z][A-Za-z0-9_.]*(?=\\()", feedback$MainFunction, perl = TRUE))))
  for (name in functions) expect_true(is.function(get(name, mode = "function")))
  fixed <- feedback[feedback$MainFunction == "mfrm_facet_intervals()", ]
  expect_match(fixed$TypicalInput, "fixed-standard-normal RSM/PCM MML", fixed = TRUE)
  expect_match(fixed$DecisionBoundary, "do not correct a biased estimate", fixed = TRUE)
  ordinary <- feedback[feedback$MainFunction == "fit_measures_table()", ]
  expect_match(ordinary$NextStep, "profile_summary_by_facet", fixed = TRUE)
  expect_match(ordinary$DecisionBoundary, "does not estimate", fixed = TRUE)
  extended <- feedback[feedback$MainFunction == "mfrm_response_diagnostics()", ]
  expect_match(extended$DecisionBoundary, "No classic cutoffs", fixed = TRUE)
  expect_match(extended$GPCMStatus, "corrected JML", fixed = TRUE)
  expect_match(extended$DecisionBoundary, "Infit defined but Outfit unavailable", fixed = TRUE)
  uncertainty <- feedback[grepl("mfrm_random_rater_intervals", feedback$MainFunction), ]
  expect_match(uncertainty$NextStep, "parm = 'raters'", fixed = TRUE)
  expect_match(uncertainty$DecisionBoundary, "unqualified for general coverage", fixed = TRUE)
  accuracy <- feedback[grepl("mfrm_screening_performance", feedback$MainFunction), ]
  expect_match(accuracy$TypicalInput, "known Affected", fixed = TRUE)
  expect_match(accuracy$DecisionBoundary, "real-data flags alone cannot", fixed = TRUE)
})

test_that("mfrmr_output_guide returns a stable purpose-to-helper map", {
  guide <- mfrmr_output_guide()

  expect_s3_class(guide, "data.frame")
  expect_true(all(c(
    "Scope",
    "Question",
    "OutputFamily",
    "Lifecycle",
    "UserLevel",
    "APILayer",
    "ObjectRole",
    "DecisionBoundary",
    "RecommendedEntry",
    "MainFunction",
    "UseWhen",
    "TypicalInput",
    "NextStep",
    "GPCMStatus",
    "Notes"
  ) %in% names(guide)))
  expect_true(nrow(guide) >= 10L)
  expect_true(any(grepl("mfrm_results", guide$MainFunction, fixed = TRUE)))
  expect_true(any(guide$RecommendedEntry))
  expect_true(any(guide$APILayer == "top_level_public_surface"))
  expect_true(any(guide$APILayer == "specialist_followup"))
  expect_true(any(guide$ObjectRole == "comprehensive result object"))
  expect_true(any(grepl("not a new estimator", guide$DecisionBoundary, fixed = TRUE)))
  expect_true(any(guide$Lifecycle == "advanced"))
  expect_true(any(guide$Lifecycle == "compatibility"))
  expect_true(any(grepl("precision_review_report", guide$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("build_summary_table_bundle", guide$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("facets_output_file_bundle", guide$MainFunction, fixed = TRUE)))
  bias <- guide[grepl("estimate_all_bias()", guide$MainFunction, fixed = TRUE), ]
  expect_true(nrow(bias) > 0L)
  expect_true(all(grepl("collection-wide", bias$DecisionBoundary, fixed = TRUE)))
  expect_true(all(grepl("unavailable", bias$DecisionBoundary, fixed = TRUE)))
  expect_true(all(grepl("plot(screens)", bias$NextStep, fixed = TRUE)))
  plans <- mfrmr_output_guide("gtheory")
  expect_match(plans$NextStep, "error_components", fixed = TRUE)
  expect_match(plans$Notes, "confint.mfrm_multivariate_d_study", fixed = TRUE)
  expect_match(plans$Notes, "IntervalCenter", fixed = TRUE)
})

test_that("mfrmr_output_guide supports focused scopes", {
  guide_columns <- c(
    "Scope", "Question", "OutputFamily", "Lifecycle", "UserLevel",
    "APILayer", "ObjectRole", "DecisionBoundary", "RecommendedEntry",
    "MainFunction", "UseWhen", "TypicalInput", "NextStep", "GPCMStatus",
    "Notes"
  )

  public <- mfrmr_output_guide("public")
  expect_identical(names(public), guide_columns)
  expect_equal(nrow(public), 6L)
  expect_true(all(public$Scope == "public"))
  expect_true(all(public$RecommendedEntry))
  expect_true(all(public$UserLevel == "beginner"))
  expect_true(all(public$APILayer == "top_level_public_surface"))
  expect_true(any(public$ObjectRole == "model estimation and result-object entry"))
  expect_identical(public$ObjectRole[2], "fit summary surface")
  expect_identical(public$ObjectRole[3], "comprehensive result object")
  expect_match(
    public$DecisionBoundary[2],
    "without computing comprehensive diagnostics",
    fixed = TRUE
  )
  expect_identical(
    public$Question,
    c(
      "1. Check your rating data and fit the model",
      "2. Read the fit summary and its recommended next step",
      "3. Review estimates, fit and precision",
      "4. Plot abilities and facet estimates with uncertainty",
      "5. Add other maps when they answer your question",
      "6. Read the report and save the reviewed results"
    )
  )
  expect_match(public$MainFunction[1], "describe_mfrm_data", fixed = TRUE)
  expect_match(public$MainFunction[1], "rating_min", fixed = TRUE)
  expect_match(public$MainFunction[1], "rating_max", fixed = TRUE)
  expect_match(public$MainFunction[1], 'method = "MML"', fixed = TRUE)
  expect_match(public$MainFunction[2], 'profile = "fit"', fixed = TRUE)
  expect_match(public$MainFunction[3], 'mfrm_results(fit)', fixed = TRUE)
  expect_match(public$MainFunction[3], "summary(res)", fixed = TRUE)
  expect_match(public$MainFunction[4], 'renderer = "native"', fixed = TRUE)
  expect_match(public$MainFunction[4], "show_ci = TRUE", fixed = TRUE)
  expect_match(public$MainFunction[4], "top_n = Inf", fixed = TRUE)
  expect_match(public$MainFunction[5], 'renderer = "facets"', fixed = TRUE)
  expect_match(public$MainFunction[5], "category_labels", fixed = TRUE)
  expect_match(public$MainFunction[5], "show_ci = FALSE", fixed = TRUE)
  expect_match(public$MainFunction[5], 'type = "fit_pathway"', fixed = TRUE)
  expect_match(public$MainFunction[5], 'fit_stat = "Infit"', fixed = TRUE)
  expect_match(public$MainFunction[5], "include_person = TRUE", fixed = TRUE)
  expect_true(any(grepl("mfrm_report", public$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("export_mfrm_results", public$MainFunction, fixed = TRUE)))
  expect_false(any(grepl("launch_mfrmr_viewer", public$MainFunction, fixed = TRUE)))
  expect_false(any(grepl("mfrm_results_interactive", public$MainFunction, fixed = TRUE)))

  beginner <- mfrmr_output_guide("beginner")
  expect_identical(names(beginner), guide_columns)
  expect_equal(beginner, public)
  expect_true(all(beginner$Scope == "public"))
  expect_true(all(beginner$UserLevel == "beginner"))
  expect_false(any(beginner$Scope %in% c("entry", "viewer", "binary")))

  psychometric <- mfrmr_output_guide("psychometric")
  expect_true(nrow(psychometric) > 0L)
  expect_true(all(nzchar(psychometric$DecisionBoundary)))
  expect_true(any(psychometric$ObjectRole == "focused evidence table"))

  entry <- mfrmr_output_guide("entry")
  expect_true(nrow(entry) > 0L)
  expect_true(all(entry$Scope == "entry"))
  expect_true(any(grepl("mfrm_results", entry$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("fit_mfrm", entry$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("mfrm_results_interactive", entry$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("launch_mfrmr_viewer", entry$MainFunction, fixed = TRUE)))
  expect_true(all(entry$RecommendedEntry))
  expect_true(all(entry$UserLevel == "beginner"))

  binary <- mfrmr_output_guide("binary")
  expect_true(nrow(binary) > 0L)
  expect_true(all(binary$Scope == "binary"))
  expect_true(all(binary$RecommendedEntry))
  expect_true(all(binary$UserLevel == "beginner"))
  expect_true(any(grepl("fit_mfrm", binary$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("Categories is 2", binary$NextStep, fixed = TRUE)))
  expect_true(any(grepl("Rasch", binary$Notes, fixed = TRUE)))

  viewer <- mfrmr_output_guide("viewer")
  expect_true(nrow(viewer) > 0L)
  expect_true(all(viewer$Scope == "viewer"))
  expect_true(all(viewer$OutputFamily == "viewer"))
  expect_true(all(viewer$RecommendedEntry))
  expect_true(all(viewer$UserLevel == "beginner"))
  expect_true(any(grepl("launch_mfrmr_viewer", viewer$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("include = \"publication\"", viewer$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("include = \"bias\"", viewer$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("include = \"misfit_review\"", viewer$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("include = \"linking\"", viewer$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("bias_interaction_report", viewer$NextStep, fixed = TRUE)))
  expect_true(any(grepl("does not estimate", viewer$Notes, fixed = TRUE)))

  simulation <- mfrmr_output_guide("simulation")
  expect_true(nrow(simulation) > 0L)
  expect_true(all(simulation$Scope == "simulation"))
  expect_true(all(simulation$Lifecycle == "advanced"))
  expect_true(any(grepl("simulate_mfrm_data", simulation$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("evaluate_mfrm_diagnostic_screening", simulation$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("export_summary_appendix", simulation$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_overview_rate", simulation$NextStep, fixed = TRUE)))
  expect_true(any(grepl("operating characteristics", simulation$Notes, fixed = TRUE)))

  linking <- mfrmr_output_guide("linking")
  expect_true(nrow(linking) > 0L)
  expect_true(all(linking$Scope == "linking"))
  expect_true(any(grepl("mfrm_results(fit, include = \"linking\")", linking$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("review_mfrm_anchors", linking$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("detect_anchor_drift", linking$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("build_equating_chain", linking$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("one fit", linking$Notes, fixed = TRUE)))

  network <- mfrmr_output_guide("network")
  expect_true(nrow(network) > 0L)
  expect_true(all(network$Scope == "network"))
  expect_true(all(network$Lifecycle == "advanced"))
  expect_true(any(grepl("build_mfrm_network_review", network$MainFunction, fixed = TRUE)))

  reviews <- mfrmr_output_guide("reviews")
  expect_true(nrow(reviews) > 0L)
  expect_true(all(reviews$Scope == "reviews"))
  expect_true(all(reviews$OutputFamily == "review"))
  expect_true(any(grepl("response_time_review", reviews$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("descriptive QC context",
                        reviews$DecisionBoundary,
                        fixed = TRUE)))

  response_time <- mfrmr_output_guide("response_time")
  expect_true(nrow(response_time) >= 2L)
  expect_true(any(response_time$Scope == "reviews"))
  expect_true(any(response_time$Scope == "r"))
  expect_true(any(grepl("plot(timing",
                        response_time$MainFunction,
                        fixed = TRUE)))
  expect_true(any(grepl("speed-accuracy",
                        response_time$DecisionBoundary,
                        fixed = TRUE)))
  expect_true(any(grepl("summary(res)$next_actions",
                        response_time$NextStep,
                        fixed = TRUE)))
  expect_true(any(grepl("do not alter MFRM estimates",
                        response_time$Notes,
                        fixed = TRUE)))
  expect_true(any(grepl("automatic exclusion rules",
                        response_time$DecisionBoundary,
                        fixed = TRUE)))

  reports <- mfrmr_output_guide("reports")
  expect_true(nrow(reports) > 0L)
  expect_true(all(reports$Scope == "reports"))
  expect_true(any(grepl("mfrm_report", reports$MainFunction, fixed = TRUE)))

  exports <- mfrmr_output_guide("exports")
  expect_true(nrow(exports) > 0L)
  expect_true(all(exports$Scope == "exports"))
  expect_true(any(grepl("export_mfrm_results", exports$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("export_summary_appendix", exports$MainFunction, fixed = TRUE)))

  gpcm <- mfrmr_output_guide("gpcm")
  expect_true(nrow(gpcm) > 0L)
  expect_false(any(gpcm$GPCMStatus == "supported"))
  expect_false(any(grepl("fit_bundle_blocked", gpcm$GPCMStatus, fixed = TRUE)))
  expect_false(any(grepl("rsm_pcm_only_for_strict_screening_currently", gpcm$GPCMStatus, fixed = TRUE)))
  expect_false(any(grepl("rsm_pcm_linking_route", gpcm$GPCMStatus, fixed = TRUE)))
  expect_true(any(
    grepl("evaluate_mfrm_diagnostic_screening", gpcm$MainFunction, fixed = TRUE) &
      gpcm$GPCMStatus == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("detect_anchor_drift", gpcm$MainFunction, fixed = TRUE) &
      grepl("supported_with_caveat", gpcm$GPCMStatus, fixed = TRUE)
  ))
  expect_true(any(grepl("GPCM", gpcm$NextStep, ignore.case = TRUE)))
  expect_true(any(gpcm$Scope == "gpcm" &
                    grepl("gpcm_capability_matrix", gpcm$MainFunction, fixed = TRUE)))
  expect_true(any(gpcm$Scope == "gpcm" &
                    grepl("gpcm_capability_matrix(view = 'unavailable')", gpcm$MainFunction, fixed = TRUE)))
  expect_true(any(gpcm$ObjectRole == "unavailable-route alternatives"))
  expect_true(any(grepl("inspect the source object and help page",
                        gpcm$DecisionBoundary,
                        fixed = TRUE)))
})

test_that("facets_feature_coverage separates implemented and unsupported FACETS surfaces", {
  coverage <- facets_feature_coverage()

  expect_s3_class(coverage, "data.frame")
  expect_identical(names(coverage), c(
    "FACETSArea",
    "FACETSFeature",
    "FACETSReference",
    "mfrmrRoute",
    "Status",
    "SurfaceCoverage",
    "StatisticalContract",
    "ValidationEvidence",
    "OperationalStatus",
    "Capability",
    "Limitation",
    "Alternative"
  ))
  expect_false("Priority" %in% names(coverage))
  expect_true(all(nzchar(coverage$Capability)))
  expect_true(all(nzchar(coverage$Limitation)))
  expect_true(all(nzchar(coverage$Alternative)))
  expect_true(all(nzchar(coverage$SurfaceCoverage)))
  expect_true(all(nzchar(coverage$StatisticalContract)))
  expect_true(all(nzchar(coverage$ValidationEvidence)))
  expect_true(all(nzchar(coverage$OperationalStatus)))
  expect_true(all(
    coverage$ValidationEvidence[coverage$Status %in% c(
      "implemented", "supported_with_caveat", "partial"
    )] == "external_match_not_established"
  ))
  expect_true(any(grepl(
    "documented FACETS version",
    coverage$Limitation,
    fixed = TRUE
  )))
  expect_true(nrow(coverage) >= 40L)
  expect_true(any(grepl("Table 14", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented"))
  expect_true(any(grepl("Table 7: agreement", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "supported_with_caveat" &
                    grepl("model-probability", coverage$Capability, fixed = TRUE) &
                    grepl("external Table 7 parity", coverage$Limitation, fixed = TRUE)))
  expect_true(any(grepl("Table 8.1: polytomous", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "supported_with_caveat" &
                    grepl("one RSM scale", coverage$Capability, fixed = TRUE) &
                    grepl("Multiple independent scales", coverage$Limitation, fixed = TRUE)))
  expect_true(any(grepl("Wright map", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    coverage$SurfaceCoverage ==
                      "familiar_visual_grammar_available" &
                    coverage$StatisticalContract ==
                      "package_native_not_facets_equivalent" &
                    coverage$ValidationEvidence ==
                      "external_match_not_established" &
                    coverage$OperationalStatus ==
                      "package_route_available"))
  expect_true(any(grepl("Generalizability Theory", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "supported_with_caveat" &
                    grepl("mfrm_d_study", coverage$mfrmrRoute, fixed = TRUE) &
                    grepl("IdentificationStatus", coverage$Capability, fixed = TRUE)))
  gtheory <- subset(coverage, FACETSReference == "gtheory.htm")
  expect_equal(nrow(gtheory), 1L)
  expect_match(gtheory$mfrmrRoute, "mfrm_g_study()", fixed = TRUE)
  expect_match(gtheory$mfrmrRoute, "mfrm_d_study()", fixed = TRUE)
  expect_match(gtheory$Limitation, "Component SEs/intervals use balanced ANOVA", fixed = TRUE)
  expect_match(gtheory$Limitation, "two-score covariance intervals have explicit generalized and conservative choices", fixed = TRUE)
  expect_match(gtheory$Limitation, "Joint reml_profile covariance intervals also accept identifiable incomplete sources", fixed = TRUE)
  expect_match(gtheory$Limitation, "experimental reml_profile G/Phi and SEM intervals", fixed = TRUE)
  expect_match(gtheory$Limitation, "sparse future plans", fixed = TRUE)
  expect_identical(gtheory$ValidationEvidence, "external_match_not_established")
  expect_true(any(grepl("Connectivity network graph", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("mfrm_network_analysis", coverage$mfrmrRoute, fixed = TRUE) &
                    grepl("type = \"network\"", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Residuals output file", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("write_mfrm_residual_file", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Category information function", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("category_information", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Cumulative probability curves", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("type = \"cumulative\"", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Winsteps", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "not_implemented"))
  expect_true(any(grepl("command-file parser", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "not_targeted"))

  partial <- facets_feature_coverage("partial")
  expect_true(nrow(partial) > 0L)
  expect_true(all(partial$Status == "partial"))

  caveated <- facets_feature_coverage("supported_with_caveat")
  expect_true(nrow(caveated) > 0L)
  expect_true(all(caveated$Status == "supported_with_caveat"))
  expect_true(any(grepl("Generalizability Theory", caveated$FACETSFeature, fixed = TRUE)))

  missing <- facets_feature_coverage("not_implemented")
  expect_true(nrow(missing) > 0L)
  expect_true(all(missing$Status == "not_implemented"))
})

test_that("facets_feature_coverage separates native calibration from unsupported routes", {
  coverage <- facets_feature_coverage()
  future_features <- c(
    "FACETS or third-party frozen-calibration import",
    "General threshold or step anchors and starting-value import",
    "Multiple observed scales and scale-specific PCM",
    "Nominal/multinomial response models",
    "Binomial-trial and Poisson/count response models",
    "Native multidimensional estimation and dimension-specific scores",
    "Unrestricted GPCM",
    "FACETS free-slope polytomous GPCM comparison"
  )
  future <- coverage[
    coverage$FACETSArea == "Current scope boundary" &
      coverage$FACETSFeature %in% future_features,
    ,
    drop = FALSE
  ]

  expect_setequal(future$FACETSFeature, future_features)
  expect_true(all(future$Status == "not_implemented"))
  expect_true(all(future$SurfaceCoverage == "unavailable"))
  expect_true(all(future$StatisticalContract == "not_available"))
  expect_true(all(future$ValidationEvidence == "not_applicable"))
  expect_true(all(future$OperationalStatus == "blocked"))
  expect_true(any(
    future$FACETSFeature == "FACETS or third-party frozen-calibration import" &
      grepl("Native mfrmr portable", future$Capability, fixed = TRUE) &
      grepl("FACETS or another program", future$Limitation, fixed = TRUE)
  ))
  expect_false(any(vapply(
    coverage,
    function(x) any(grepl("0.2.3", x, fixed = TRUE)),
    logical(1)
  )))
  expect_true(any(
    future$FACETSFeature == "Unrestricted GPCM" &
      grepl("requires slope_facet to equal step_facet", future$Capability, fixed = TRUE)
  ))
  expect_true(any(
    future$FACETSFeature == "Nominal/multinomial response models" &
      grepl("ordered category", future$Capability, fixed = TRUE) &
      grepl("not an unordered", future$Limitation, fixed = TRUE)
  ))
  expect_true(any(
    future$FACETSFeature == "Binomial-trial and Poisson/count response models" &
      grepl("integer scores", future$Limitation, fixed = TRUE)
  ))
  expect_true(any(
    future$FACETSFeature == "FACETS free-slope polytomous GPCM comparison" &
      grepl("post-fit diagnostic", future$Limitation, fixed = TRUE) &
      grepl("PCM/JML comparison", future$Alternative, fixed = TRUE)
  ))
})

test_that("facets_feature_coverage separates frequency weights from count outcomes", {
  coverage <- facets_feature_coverage()
  frequency <- coverage[
    coverage$FACETSFeature == "Row-frequency weights for ordered ratings",
    , drop = FALSE
  ]

  expect_equal(nrow(frequency), 1L)
  expect_identical(frequency$Status, "supported_with_caveat")
  expect_match(frequency$mfrmrRoute, "fit_mfrm(weight = ...)", fixed = TRUE)
  expect_match(frequency$Capability, "positive numeric observation weight", fixed = TRUE)
  expect_match(frequency$Capability, "row-replication weight", fixed = TRUE)
  expect_match(frequency$Limitation, "not a general collapsed-person frequency table", fixed = TRUE)
  expect_match(frequency$Limitation, "after marginalization", fixed = TRUE)
  expect_match(frequency$Limitation, "does not create a count-response family", fixed = TRUE)
  expect_match(frequency$Limitation, "information-criterion panel", fixed = TRUE)
})

test_that("fit_mfrm signature does not expose later-release calibration routes", {
  current_args <- names(formals(fit_mfrm))
  later_release_args <- c(
    "calibration_bundle",
    "frozen_calibration",
    "threshold_anchors",
    "step_anchors",
    "scale_id",
    "dimensions"
  )

  expect_false("..." %in% current_args)
  expect_false(any(later_release_args %in% current_args))
})

test_that("facets_positioning_guide prevents FACETS numerical-clone wording", {
  guide <- facets_positioning_guide()

  expect_s3_class(guide, "data.frame")
  expect_true(all(c(
    "Topic",
    "Position",
    "RecommendedWording",
    "PrimaryRoute"
  ) %in% names(guide)))
  expect_true(any(guide$Topic == "Estimation authority"))
  expect_true(any(guide$Topic == "External FACETS comparison"))
  expect_true(any(guide$Topic == "Current model and calibration boundary"))
  expect_true(any(grepl("package-native", guide$Position, fixed = TRUE)))
  expect_true(any(grepl("not evidence of FACETS numerical equivalence",
                        guide$RecommendedWording,
                        fixed = TRUE)))
  expect_true(any(grepl("read_facets_fit_table", guide$PrimaryRoute, fixed = TRUE)))
  expect_true(any(grepl("not as a general FACETS operational-calibration replacement",
                        guide$RecommendedWording,
                        fixed = TRUE)))
  expect_true(any(grepl("importing FACETS or third-party frozen calibrations",
                        guide$Position,
                        fixed = TRUE)))
  expect_true(any(grepl("existing fitted object",
                        guide$Position,
                        fixed = TRUE)))
  expect_false(any(grepl("\\baudit\\b", unlist(guide), ignore.case = TRUE)))
})

test_that("FACETS term and visual contracts separate visual and numerical correspondence", {
  terms <- facets_term_crosswalk()
  expect_true(all(c(
    "FACETSTerm", "mfrmrTerm", "mfrmrRoute", "Relationship", "Boundary"
  ) %in% names(terms)))
  expect_true(any(terms$FACETSTerm == "Measure"))
  expect_true(any(terms$FACETSTerm == "ZSTD" &
                    grepl("convention", terms$Relationship, fixed = TRUE)))

  visuals <- facets_visual_contract()
  expect_true(all(c(
    "FACETSVisualSurface", "Status", "FirstMfrmrRoute",
    "EditableDataRoute", "ClaimBoundary"
  ) %in% names(visuals)))
  expect_true(any(grepl("asterisk Wright", visuals$FACETSVisualSurface, fixed = TRUE) &
                    grepl('renderer = "facets"', visuals$FirstMfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("uncertainty", visuals$FACETSVisualSurface, fixed = TRUE) &
                    visuals$Status == "mfrmr_extension"))
  expect_true(any(grepl("ruler_rows", visuals$EditableDataRoute, fixed = TRUE)))
  expect_false(any(grepl('component = "ruler"', visuals$EditableDataRoute, fixed = TRUE)))
  expect_true(any(grepl(
    "documented FACETS version",
    visuals$ClaimBoundary,
    fixed = TRUE
  )))
})

test_that("mfrmr_output_guide gives FACETS, ConQuest, and R user pathways", {
  facets <- mfrmr_output_guide("facets")
  expect_true(nrow(facets) >= 14L)
  expect_true(all(facets$Scope == "facets"))
  expect_match(facets$MainFunction[1], 'method = "MML"', fixed = TRUE)
  expect_match(facets$MainFunction[1], 'profile = "facets"', fixed = TRUE)
  expect_match(facets$MainFunction[2], 'renderer = "native"', fixed = TRUE)
  expect_match(facets$MainFunction[2], "show_ci = TRUE", fixed = TRUE)
  expect_match(facets$MainFunction[2], "top_n = Inf", fixed = TRUE)
  expect_match(facets$MainFunction[3], 'renderer = "facets"', fixed = TRUE)
  expect_match(facets$MainFunction[3], "category_labels", fixed = TRUE)
  expect_match(facets$MainFunction[3], "show_ci = FALSE", fixed = TRUE)
  expect_match(facets$MainFunction[4], 'type = "fit_pathway"', fixed = TRUE)
  expect_match(facets$MainFunction[4], 'fit_stat = "Infit"', fixed = TRUE)
  expect_match(facets$MainFunction[4], "include_person = TRUE", fixed = TRUE)
  expect_match(facets$Notes[1], "FACETS software was run", fixed = TRUE)
  expect_match(facets$Notes[3], "labelled step lines", fixed = TRUE)
  expect_match(facets$Notes[3], "uncertainty hybrid", fixed = TRUE)
  expect_true(any(grepl("facets_positioning_guide", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("not a FACETS numerical clone", facets$Notes, fixed = TRUE)))
  expect_true(any(grepl("facets_feature_coverage", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("review_mfrm_anchors", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("make_anchor_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("group_anchors", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("detect_anchor_drift", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot(drift", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("fit_measures_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("facets_fit_df_guide", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("df_sensitivity", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("df_sensitive", facets$NextStep, fixed = TRUE)))
  expect_true(any(grepl("rating_scale_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("fair_average_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("estimate_bias", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("bias_interaction_report", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("bias_pairwise_report", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot(bias", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("view = \"unified\"", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("type = \"wright\"", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("data_quality_report", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("write_mfrm_residual_file", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("write_mfrm_subset_file", facets$MainFunction, fixed = TRUE)))

  conquest <- mfrmr_output_guide("conquest")
  expect_true(nrow(conquest) >= 5L)
  expect_true(all(conquest$Scope == "conquest"))
  expect_match(conquest$MainFunction[1], 'method = "MML"', fixed = TRUE)
  expect_match(conquest$Notes[1], "does not execute ConQuest", fixed = TRUE)
  expect_match(conquest$MainFunction[2], "build_conquest_overlap_bundle", fixed = TRUE)
  expect_match(conquest$Notes[2], "does not execute ConQuest", fixed = TRUE)
  expect_match(conquest$Notes[2], "parse raw report text", fixed = TRUE)
  expect_true(any(grepl("build_conquest_overlap_bundle", conquest$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("normalize_conquest_overlap_exports", conquest$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("review_conquest_overlap", conquest$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("less flexible than ConQuest", conquest$Question, fixed = TRUE)))

  r_path <- mfrmr_output_guide("r")
  expect_true(nrow(r_path) >= 3L)
  expect_true(all(r_path$Scope == "r"))
  expect_true(any(grepl("plot_data", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_data_components", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("mfrmr_interval_guide", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_long", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("type = \"pathway\"", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot(bias", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot(info", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot(timing", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_data_components", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("mfrmr_interval_guide", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("pathway_long", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("annotations/settings", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("table, thresholds, overview, and notes",
                        r_path$NextStep,
                        fixed = TRUE)))
  expect_true(any(grepl("ggplot2", r_path$UseWhen, fixed = TRUE)))
})

test_that("plot_data extracts reusable payloads and selected components", {
  fit <- make_toy_fit(maxit = 8)
  dq <- data_quality_report(fit)

  dashboard_payload <- plot_data(dq, type = "dashboard")
  expect_true(is.list(dashboard_payload))
  expect_true(all(c(
    "quality_flags",
    "score_support",
    "facet_response_patterns"
  ) %in% names(dashboard_payload)))

  flags <- plot_data(dq, type = "dashboard", component = "quality_flags")
  expect_s3_class(flags, "data.frame")
  dashboard_components <- plot_data_components(dq, type = "dashboard")
  expect_s3_class(dashboard_components, "data.frame")
  expect_true(all(c(
    "PlotName", "Component", "Role", "ObjectType", "Rows",
    "Columns", "Accessor", "Notes"
  ) %in% names(dashboard_components)))
  expect_true(any(dashboard_components$Component == "quality_flags"))

  plotted <- plot(dq, type = "dashboard", draw = FALSE)
  expect_s3_class(plotted, "mfrm_plot_data")
  score_support <- plot_data(plotted, component = "score_support")
  expect_s3_class(score_support, "data.frame")
  plotted_components <- plot_data_components(plotted)
  expect_true(any(plotted_components$Component == "score_support"))
  expect_error(plot_data(plotted, component = "not_a_component"), "component")

  curves <- category_curves_report(fit, theta_points = 21)
  cumulative_payload <- plot_data(curves, type = "cumulative")
  expect_true(all(c(
    "cumulative_probabilities",
    "cumulative_boundaries",
    "cumulative_direction"
  ) %in% names(cumulative_payload)))
  cumulative_table <- plot_data(curves, type = "cumulative", component = "cumulative_probabilities")
  expect_s3_class(cumulative_table, "data.frame")
  category_probability_payload <- plot_data(curves, type = "category_probability")
  expect_identical(category_probability_payload$plot, "ccc")
  expect_identical(category_probability_payload$plot_settings$RequestedType[1], "category_probability")
  expect_true(any(category_probability_payload$plot_long$PlotType == "ccc" &
                    category_probability_payload$plot_long$DisplayedByDefault))
  expect_true(all(c("plot_annotations", "curve_summary") %in%
                    names(category_probability_payload)))
  category_probability_components <- plot_data_components(curves, type = "category_probability")
  expect_true(any(category_probability_components$Component == "plot_long" &
                    category_probability_components$Role == "primary_data"))
  expect_true(any(category_probability_components$Component == "plot_annotations" &
                    category_probability_components$Role == "annotation"))
  category_information <- plot_data(
    curves,
    type = "category_information",
    component = "category_information"
  )
  expect_s3_class(category_information, "data.frame")
  expect_true(all(c("CategoryInformation", "CategoryInformationShare") %in% names(category_information)))

  curve_long <- plot_data(curves, component = "plot_long")
  expect_s3_class(curve_long, "data.frame")
  expect_true(all(c(
    "PlotType", "Panel", "CurveGroup", "Theta", "Series",
    "ValueName", "Value", "DisplayedByDefault"
  ) %in% names(curve_long)))
  expect_true(all(c(
    "ogive", "ccc", "cumulative", "information", "category_information"
  ) %in% unique(curve_long$PlotType)))
  expect_true(any(curve_long$PlotType == "cumulative" &
                    curve_long$Direction == "at_or_below" &
                    curve_long$DisplayedByDefault))
  curve_style <- plot_data(curves, preset = "monochrome", component = "curve_style")
  expect_s3_class(curve_style, "data.frame")
  expect_true(all(c("Series", "Colour", "LineType", "Preset") %in% names(curve_style)))
  expect_true(all(curve_style$Preset == "monochrome"))
  expect_gt(length(unique(curve_style$LineType)), 1L)

  diag <- diagnose_mfrm(fit, residual_pca = "none")
  pathway_long <- .mfrmr_muffle_expected_warnings(
    plot_data(
      fit,
      type = "pathway",
      diagnostics = diag,
      component = "pathway_long"
    ),
    "^Review-only display:"
  )
  expect_s3_class(pathway_long, "data.frame")
  expect_true(all(c("Layer", "CurveGroup", "Theta", "Value") %in% names(pathway_long)))
  expect_true(any(pathway_long$Layer == "expected_score"))

  pathway_fit <- .mfrmr_muffle_expected_warnings(
    plot_data(
      fit,
      type = "pathway",
      diagnostics = diag,
      component = "fit_measures"
    ),
    "^Review-only display:"
  )
  expect_s3_class(pathway_fit, "data.frame")
  expect_true(all(c("Facet", "Level", "Infit", "Outfit", "FitStatus") %in% names(pathway_fit)))
  pathway_components <- .mfrmr_muffle_expected_warnings(
    plot_data_components(fit, type = "pathway", diagnostics = diag),
    "^Review-only display:"
  )
  expect_true(any(pathway_components$Component == "pathway_long" &
                    pathway_components$Role == "primary_data"))
  expect_true(any(pathway_components$Component == "fit_measures" &
                    pathway_components$Role == "fit_review"))

  pathway_without_fit <- .mfrmr_muffle_expected_warnings(
    plot_data(
      fit,
      type = "pathway",
      include_fit_measures = FALSE,
      component = "fit_measure_status"
    ),
    "^Review-only display:"
  )
  expect_s3_class(pathway_without_fit, "data.frame")
  expect_false(isTRUE(pathway_without_fit$Available[1]))
  expect_identical(pathway_without_fit$Status[1], "not_requested")

  info <- compute_information(fit, theta_points = 21)
  sem_long <- plot_data(
    plot(info, type = "sem", draw = FALSE),
    component = "plot_long"
  )
  expect_s3_class(sem_long, "data.frame")
  expect_true(any(sem_long$ValueName == "ConditionalSEM"))
  sem_components <- plot_data_components(plot(info, type = "sem", draw = FALSE))
  expect_true(any(sem_components$Component == "conditional_sem"))

  bias <- suppressWarnings(suppressMessages(
    estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 1)
  ))
  bias_long <- plot_data(
    plot(bias, type = "heatmap", draw = FALSE),
    component = "plot_long"
  )
  expect_s3_class(bias_long, "data.frame")
  expect_true(any(bias_long$Layer == "heatmap_cell"))
  bias_components <- plot_data_components(
    plot(bias, type = "heatmap", draw = FALSE)
  )
  expect_true(any(bias_components$Component == "flag_summary" &
                    bias_components$Role == "summary_or_guidance"))
})
