#' FACETS Positioning Guide
#'
#' @description
#' `facets_positioning_guide()` gives user-facing wording for the relationship
#' between `mfrmr` and FACETS. Use it when a report, migration note, or
#' methods appendix must make clear that `mfrmr` is not a FACETS numerical
#' clone.
#'
#' @details
#' The guide separates five ideas that are easy to conflate:
#'
#' - estimation authority: fitted values come from `mfrmr` unless external
#'   FACETS output is explicitly supplied;
#' - compatibility purpose: FACETS-style names and files are transition,
#'   handoff, and report-organization surfaces;
#' - external comparison: FACETS comparisons require a supplied external table
#'   and should separate MnSq differences from df/ZSTD convention differences;
#' - current model boundary: one ordered-categorical response-model family and
#'   one observed score scale are used per fit; nominal and count-response
#'   families, mixed families, multiple independent scales,
#'   general threshold anchoring, and scoring from an imported versioned frozen
#'   calibration are not part of the current public estimator; posterior
#'   scoring from an existing fitted object is a separate current route;
#' - extension surface: native R tables, plot data, GPCM diagnostics,
#'   network views, and G/D-study helpers are package extensions, not promises
#'   of FACETS menu-level reproduction.
#'
#' @return A data.frame with columns:
#' - `Topic`
#' - `Position`
#' - `RecommendedWording`
#' - `PrimaryRoute`
#'
#' @seealso [facets_feature_coverage()], [mfrmr_output_guide()],
#'   [read_facets_fit_table()], [facets_fit_review()]
#' @examples
#' facets_positioning_guide()
#' @export
facets_positioning_guide <- function() {
  data.frame(
    Topic = c(
      "Estimation authority",
      "Compatibility purpose",
      "External FACETS comparison",
      "Current model and calibration boundary",
      "Reporting source of truth",
      "Extension beyond FACETS"
    ),
    Position = c(
      "mfrmr estimates are package-native; FACETS-style names do not mean that FACETS estimated the model.",
      "FACETS-style wrappers, table labels, and files support transition, handoff, and report organization, not optimizer-level reproduction.",
      "Numerical comparison requires an explicit external FACETS output table supplied by the user.",
      "Each fit uses one ordered-categorical response-model family and one observed score scale; nominal and count-response families, mixed families, multiple independent scales, general threshold anchors, and scoring from an imported versioned frozen calibration are not current public capabilities. Positive observation weights act on row-level conditional ordered-rating contributions; they are not a general collapsed-person frequency-table interface and do not change the response family. Posterior scoring from an existing fitted object is a separate route.",
      "Inference and reporting should be based on native fit, diagnostics, review, table, and plot-data objects.",
      "GPCM, D-study, network, and reusable visualization data are extension routes rather than FACETS menu clones."
    ),
    RecommendedWording = c(
      "The model was estimated with mfrmr; FACETS-style output names are used only to organize the report.",
      "FACETS-style outputs were generated for handoff or reader familiarity; they are not evidence of FACETS numerical equivalence.",
      "When external FACETS output is supplied, compare MnSq first and report df/ZSTD convention sensitivity separately.",
      "Describe mfrmr as a native R RSM/PCM analysis, diagnostic, and reporting environment, not as a general FACETS operational-calibration replacement.",
      "Report estimates, standard errors, fit summaries, and plots from documented mfrmr objects.",
      "Use package-native extensions as additional evidence and label them as mfrmr analyses."
    ),
    PrimaryRoute = c(
      "fit_mfrm(); diagnose_mfrm(); reporting_checklist()",
      "facets_feature_coverage(); run_mfrm_facets(); facets_output_file_bundle()",
      "read_facets_fit_table(); facets_fit_review(); fit_measures_table(df_sensitivity = TRUE)",
      "fit_mfrm(); facets_feature_coverage(); mfrmr_output_guide()",
      "build_summary_table_bundle(); build_visual_summaries(); plot_data()",
      "gpcm_capability_matrix(); mfrm_d_study(); mfrm_network_analysis(); plot_data_components()"
    ),
    stringsAsFactors = FALSE
  )
}

#' FACETS-to-mfrmr term crosswalk
#'
#' @description
#' `facets_term_crosswalk()` maps common FACETS report and specification terms
#' to their closest `mfrmr` objects or routes. The relationship column makes
#' explicit whether a row is a substantive counterpart, a presentation
#' convention, or only a migration aid.
#'
#' @return A data.frame with `FACETSTerm`, `mfrmrTerm`, `mfrmrRoute`,
#'   `Relationship`, and `Boundary` columns.
#' @seealso [facets_positioning_guide()], [facets_feature_coverage()],
#'   [facets_visual_contract()]
#' @examples
#' facets_term_crosswalk()
#' @export
facets_term_crosswalk <- function() {
  data.frame(
    FACETSTerm = c(
      "Measure", "S.E.", "Infit MnSq", "Outfit MnSq", "ZSTD",
      "Rating Scale=", "Table 6 Wright/variable map", "Graphfile=",
      "Anchorfile=", "Positive="
    ),
    mfrmrTerm = c(
      "Estimate (logits)", "SE", "Infit", "Outfit", "InfitZSTD / OutfitZSTD",
      "model = \"RSM\" or \"PCM\"", "Wright map", "graph output table",
      "anchor table", "facet orientation"
    ),
    mfrmrRoute = c(
      "fit_measures_table(); summary(fit)",
      "fit_measures_table(); plot(fit, type = \"wright\", show_ci = TRUE)",
      "diagnose_mfrm(); fit_measures_table()",
      "diagnose_mfrm(); fit_measures_table()",
      "facets_fit_df_guide(); fit_measures_table(df_sensitivity = TRUE)",
      "fit_mfrm(model = \"RSM\" or \"PCM\")",
      "plot_wright_unified(); plot(fit, type = \"wright\")",
      "facets_output_file_bundle(include = \"graph\")",
      "make_anchor_table(); fit_mfrm(anchors = ...)",
      "fit_mfrm(positive_facets = ...); plot_wright_unified()"
    ),
    Relationship = c(
      "substantive counterpart", "substantive counterpart", "substantive counterpart",
      "substantive counterpart", "convention-sensitive counterpart",
      "model-setting crosswalk", "visual counterpart", "handoff-file counterpart",
      "input-table counterpart", "orientation convention"
    ),
    Boundary = c(
      "Numerical equality requires aligned model, estimator, identification, and supplied external FACETS output.",
      "SE bases differ by estimator and element type; inspect method metadata.",
      "Compare MnSq before standardized fit and document residual basis.",
      "Compare MnSq before standardized fit and document residual basis.",
      "df and Wilson-Hilferty conventions can change ZSTD without changing MnSq.",
      "The package exposes documented RSM/PCM and bounded-GPCM routes, not a FACETS command parser.",
      "FACETS-style rendering reproduces the ruler grammar, not optimizer-level numerical identity.",
      "CSV/TSV handoff is package-native rather than fixed-field FACETS syntax.",
      "The table is an R-native anchor contract, not a complete FACETS specification file.",
      "Always report which facets use the positive orientation."
    ),
    stringsAsFactors = FALSE
  )
}

#' FACETS-facing visual contract
#'
#' @description
#' `facets_visual_contract()` identifies the closest package route for common
#' FACETS visual surfaces and states what can and cannot be claimed from that
#' visual correspondence. It distinguishes the FACETS-style asterisk ruler
#' from the native uncertainty-aware Wright map.
#'
#' @return A data.frame with visual surface, status, first route, editable data
#'   route, and claim-boundary columns.
#' @seealso [plot_wright_unified()], [plot_data()], [facets_term_crosswalk()],
#'   [facets_feature_coverage()]
#' @examples
#' facets_visual_contract()
#' @export
facets_visual_contract <- function() {
  data.frame(
    FACETSVisualSurface = c(
      "Table 6 asterisk Wright ruler",
      "Native Wright map with uncertainty",
      "Table 8 rating-scale structure",
      "Graphs: expected-score ICC/IRF",
      "Bond-Fox fit pathway",
      "DIF/bias visual review",
      "Graphfile and R/Web handoff"
    ),
    Status = c(
      "implemented_visual", "mfrmr_extension", "implemented",
      "implemented", "implemented_extension", "partial", "partial"
    ),
    FirstMfrmrRoute = c(
      "plot_wright_unified(fit, renderer = \"facets\")",
      "plot_wright_unified(fit, renderer = \"native\", show_ci = TRUE)",
      "rating_scale_table(); category_structure_report()",
      "plot(fit, type = \"pathway\")",
      "plot(fit, type = \"fit_pathway\", fit_stat = \"Infit\", include_person = TRUE, top_n_person = 12)",
      "plot_bias_interaction()",
      "facets_output_file_bundle(); plot_data(); as_ggplot()"
    ),
    EditableDataRoute = c(
      "plot(..., draw = FALSE)$data$facets_style$ruler_rows",
      "plot(..., draw = FALSE); plot_data(component = \"locations\")",
      "rating_scale_table(); plot_data()",
      "plot_data(type = \"pathway\")",
      "plot_data(type = \"fit_pathway\")",
      "plot_data(); bias_interaction_report()",
      "facets_output_file_bundle(); plot_data_components()"
    ),
    ClaimBoundary = c(
      "FACETS-style visual grammar; numerical comparison requires output from a documented FACETS version under aligned settings.",
      "The SE/CI display is an mfrmr extension and should be retained for uncertainty interpretation.",
      "Structured R output replaces FACETS line-printer artwork.",
      "This is expected score over theta, not a measure-versus-fit pathway map.",
      "This is a Bond-Fox-style mfrmr extension; it is not a standard FACETS Table 6 output.",
      "Screening display only; it is not a final fairness conclusion.",
      "Editable handoff is supported, but FACETS UI, Excel, and webpage behavior is not cloned."
    ),
    stringsAsFactors = FALSE
  )
}

#' FACETS Feature Coverage Matrix
#'
#' @description
#' `facets_feature_coverage()` summarizes how `mfrmr` maps
#' the main FACETS output-table, output-file, and graph-menu surface to package
#' functions. It is a surface-coverage guide, not a statement that estimands,
#' conditioning, extreme-score handling, degrees of freedom, or numerical
#' results are equivalent.
#'
#' Use this helper before migration work when you need a public, user-facing
#' answer to three questions:
#'
#' - which FACETS outputs have a close `mfrmr` route,
#' - which outputs are only partially covered by structured R objects,
#' - which FACETS-specific outputs are not implemented or intentionally outside
#'   the current package scope.
#'
#' @param status Which rows to return. `"all"` returns the full matrix.
#'   Other values filter by the `Status` column.
#'
#' @details
#' The matrix is based on the FACETS 64-bit output index, which lists output
#' Tables 1--14, DIF/bias plots, R/Web plots, output files, and graph-menu
#' curves. `mfrmr` intentionally prioritizes structured R tables and reusable
#' plot data over exact FACETS line-printer output.
#'
#' The current software reference target is FACETS 64-bit 4.5.1 (July 2026).
#' Bibliographic references retain the title and edition of the consulted
#' manual rather than silently relabelling a 4.5.0 manual as 4.5.1. External
#' numerical validation is a separate evidence contract.
#'
#' Status meanings:
#'
#' - `implemented`: a package-native route covers the substantive output
#'   surface; this status alone does not claim an externally matched
#'   statistical contract.
#' - `supported_with_caveat`: a package-native route exists, but the output
#'   must be read with explicit identification, validation, or scope caveats.
#' - `partial`: the concept is covered, but not the full FACETS formatting,
#'   option surface, file type, or external integration.
#' - `not_implemented`: a FACETS feature has no direct package-native route in
#'   the documented package scope.
#' - `not_targeted`: the feature is tied to FACETS UI, Web/Excel handoff, or
#'   another external program format and is outside the package scope.
#'
#' The four contract axes are deliberately independent:
#'
#' - `SurfaceCoverage` records whether a corresponding package surface exists;
#'   familiar visual grammar belongs only to this axis.
#' - `StatisticalContract` records the package-native statistical scope and
#'   never implies a FACETS-matched estimand.
#' - `ValidationEvidence` records whether this table establishes matched
#'   external numerical evidence. It currently does not; validation belongs to
#'   a separate candidate-linked evidence contract.
#' - `OperationalStatus` records route availability. A package route being
#'   available does not make `mfrmr` operationally interchangeable with
#'   FACETS.
#'
#' @return A data.frame with columns:
#' - `FACETSArea`
#' - `FACETSFeature`
#' - `FACETSReference`
#' - `mfrmrRoute`
#' - `Status`
#' - `SurfaceCoverage`
#' - `StatisticalContract`
#' - `ValidationEvidence`
#' - `OperationalStatus`
#' - `Capability`
#' - `Limitation`
#' - `Alternative`
#'
#' @references
#' Linacre, J. M. (2026). *A user's guide to FACETS, version 4.5.0*.
#' Winsteps.com. See the guide's output-table index for the documented
#' FACETS files, plots, and graphs.
#'
#' @seealso [facets_positioning_guide()], [mfrmr_output_guide()],
#'   [facets_fit_df_guide()], [read_facets_fit_table()], [facets_fit_review()],
#'   [gpcm_capability_matrix()]
#' @examples
#' facets_feature_coverage()
#' facets_feature_coverage("partial")
#' facets_feature_coverage("supported_with_caveat")
#' facets_feature_coverage("not_implemented")
#' @export
facets_feature_coverage <- function(status = c("all", "implemented",
                                               "supported_with_caveat",
                                               "partial", "not_implemented",
                                               "not_targeted")) {
  status <- match.arg(status)

  package_native_alternative <- paste(
    "Use the documented mfrmr route;",
    "use FACETS externally only when its exact layout or option surface is required."
  )
  external_format_alternative <- paste(
    "Use the closest documented mfrmr table or plot for analysis;",
    "use FACETS externally when the omitted statistic or format is required."
  )
  custom_r_alternative <- paste(
    "Build the display from package plot data or another R graphics system;",
    "use FACETS externally for its menu-specific view."
  )
  external_program_alternative <- paste(
    "Use package-native R output where suitable, or run the relevant external",
    "program for its program-specific file, interface, or report."
  )

  row <- function(area, feature, reference, route, status,
                  capability, limitation, alternative) {
    data.frame(
      FACETSArea = area,
      FACETSFeature = feature,
      FACETSReference = reference,
      mfrmrRoute = route,
      Status = status,
      Capability = capability,
      Limitation = limitation,
      Alternative = alternative,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, list(
    row("Output table", "Table 1: specification summary", "table1.htm",
        "specifications_report()", "implemented",
        "Structured run settings and reproducibility context.",
        "Not an exact FACETS line-printer layout.", package_native_alternative),
    row("Output table", "Table 2: data summary report", "table2.htm",
        "data_quality_report(); describe_mfrm_data()", "implemented",
        "Rows, exclusions, missingness, score support, and response-pattern QC.",
        "Structured QC replaces FACETS text layout.", package_native_alternative),
    row("Output table", "Table 3: main iteration report", "table3.htm",
        "estimation_iteration_report()", "partial",
        "Convergence and replayed iteration evidence.",
        "Does not reproduce the complete FACETS iteration trace.", package_native_alternative),
    row("Output table", "Table 4: unexpected responses", "table4.htm",
        "unexpected_response_table(); plot_unexpected()", "implemented",
        "Case-level unexpected-response screening.",
        "Structured table and plots, not printer-identical FACETS output.", package_native_alternative),
    row("Output table", "Table 5: measurable data summary", "table5.htm",
        "measurable_summary_table(); describe_mfrm_data()", "implemented",
        "Facet coverage, category counts, and subset/connectivity checks.",
        "Column order and text layout differ from FACETS.", package_native_alternative),
    row("Output table", "Table 6.0: all-facet Wright map rulers", "table6.htm",
        "plot_wright_unified(renderer = \"facets\"); plot(fit, type = \"wright\")", "implemented",
        "Common-logit person/facet/threshold display with FACETS-style asterisk ruler or native SE/CI rendering.",
        "Visual correspondence does not establish numerical equivalence; compare output from a documented FACETS version under aligned settings.",
        package_native_alternative),
    row("Output table", "Table 6.0.0: disjoint element listing", "table6_0_0.htm",
        "subset_connectivity_report()", "implemented",
        "Disconnected subsets and facet-by-subset coverage.",
        "Network-style graph is not the default display.", package_native_alternative),
    row("Output table", "Table 6.2: graphical facet statistics", "table6_2.htm",
        "facet_statistics_report(); plot(...)", "partial",
        "Facet statistics and visual summaries.",
        "FACETS M/S/Q/X printer-graph formatting is not reproduced exactly.", package_native_alternative),
    row("Output table", "Table 7: facet measurement report", "table7.htm",
        "fit_measures_table(); diagnose_mfrm(); summary(fit)", "supported_with_caveat",
        "Measures, SEs, fit, anchoring status, and review flags.",
        "Estimator and person-score basis, extreme handling, df/ZSTD conventions, and FACETS options can differ; external numerical equivalence is not established.", package_native_alternative),
    row("Output table", "Table 7: reliability and chi-square", "table7summarystatistics.htm",
        "facets_chisq_table(); diagnose_mfrm()$reliability", "supported_with_caveat",
        "Rasch/FACETS-style separation, reliability, and chi-square summaries.",
        "Uses package-native structured output and documented df conventions; external numerical equivalence is not established.", package_native_alternative),
    row("Output table", "Table 7: agreement statistics", "table7agreementstatistics.htm",
        "interrater_agreement_table(); rater_network_analysis(); rater_halo_network_analysis(); plot_interrater_agreement()", "supported_with_caveat",
        "Raw-category observed agreement and model-probability expected agreement, plus pairwise rater-network, rater-by-criterion halo network, and rater-agreement views.",
        "The package does not translate category positions across multiple independent scales, apply FACETS agreement-based SE inflation, or establish external Table 7 parity.", package_native_alternative),
    row("Output table", "Table 8.1: dichotomous/binomial/Poisson statistics",
        "table8_1dichotomous.htm",
        "rating_scale_table() for two-category ordered scores", "partial",
        "Two-category Rasch-category summaries are available.",
        "FACETS binomial-trial and Poisson-specific reports are not implemented.", external_format_alternative),
    row("Output table", "Table 8.1: polytomous rating-scale/partial-credit statistics",
        "table8_1ratingscale.htm",
        "rating_scale_table(); category_structure_report()", "supported_with_caveat",
        "Core category diagnostics and thresholds for one RSM scale or one PCM step facet on the package's shared observed score scale.",
        "Multiple independent scales, scale-specific anchors or starting values, Thurstone thresholds, and FACETS text layout are not reproduced.", package_native_alternative),
    row("Output table", "Table 8: scale-structure bar chart", "table8barchart.htm",
        "category_structure_report()", "partial",
        "Category structure and transition summaries.",
        "FACETS line-printer artwork is not reproduced exactly.", package_native_alternative),
    row("Output table", "Table 8: scale-structure probability curves", "table8curves.htm",
        "category_curves_report(); plot(fit, type = \"ccc\")", "implemented",
        "Category probability and expected-score curve data.",
        "Uses R-native plot data rather than FACETS graph text.", package_native_alternative),
    row("Output table", "Table 9: bias-estimation iteration report", "table9.htm",
        "estimate_bias(); bias_iteration_report()", "implemented",
        "Bias recalibration path and final iteration status.",
        "Conditional screening semantics are documented separately.", package_native_alternative),
    row("Output table", "Table 10: unexpected after allowing for bias", "table10.htm",
        "unexpected_after_bias_table()", "implemented",
        "Unexpected rows after the current bias-screening layer.",
        "Structured table replaces FACETS text layout.", package_native_alternative),
    row("Output table", "Table 11: bias-calculation counts", "table11.htm",
        "bias_count_table()", "implemented",
        "Response counts behind bias estimates.",
        "Structured output replaces FACETS text layout.", package_native_alternative),
    row("Output table", "Table 12: bias summary report", "table12.htm",
        "summary(estimate_bias(...)); plot_bias_interaction()", "partial",
        "Distributional and visual bias summaries are available.",
        "FACETS vertical frequency bar-chart is not reproduced exactly.", package_native_alternative),
    row("Output table", "Table 13: DIF/bias detail report", "table13.htm",
        "estimate_bias(); bias_interaction_report()", "implemented",
        "Ranked cell-level bias/interactions with screening statistics.",
        "Reported as screening evidence, not final fairness inference.", package_native_alternative),
    row("Output table", "Table 14: pairwise bias report", "table14.htm",
        "bias_pairwise_report(); build_fixed_reports()", "implemented",
        "Pairwise contrasts for two-way bias runs.",
        "Higher-order runs omit pairwise sections by design.", package_native_alternative),
    row("Output table", "DIF/bias Excel plot", "difbiasplot.htm",
        "plot_bias_interaction(plot = ...)", "partial",
        "R-native scatter, heatmap, and facet-profile bias displays.",
        "Excel-specific output is not implemented.", external_format_alternative),
    row("R/Web plots", "Scatterplots and histograms from FACETS menus",
        "outputtableindex.htm",
        "plot_data(); package plot helpers", "partial",
        "Reusable plot data supports custom R graphics.",
        "FACETS arbitrary R/Web plotting menus are not mirrored.", custom_r_alternative),
    row("R/Web plots", "X-Y plot: R Statistics", "xyplotr.htm",
        "plot_data(); user-defined R plotting", "partial",
        "Users can build X-Y plots from returned data frames.",
        "No dedicated FACETS-style arbitrary X-Y plot wrapper.", custom_r_alternative),
    row("R/Web plots", "X-Y plot: Webpage", "xyplotwebpage.htm",
        "none", "not_targeted",
        "No package-native Webpage plot generator.",
        "Webpage menu output is a FACETS UI feature.", external_program_alternative),
    row("R/Web plots", "X-Y-Z plot: R Statistics", "xyzplotr.htm",
        "plot(fit, type = \"ccc_surface\"); plot_data()", "partial",
        "Selected 3D/surface-ready plot data are available.",
        "No arbitrary FACETS X-Y-Z plot wrapper.", custom_r_alternative),
    row("R/Web plots", "Histogram: R Statistics", "histogramr.htm",
        "plot_data(); plot(fit, type = \"wright\"); plot_qc_dashboard()", "partial",
        "Several package outputs include histogram-like summaries.",
        "No general FACETS histogram menu clone.", custom_r_alternative),
    row("R/Web plots", "Generalizability Theory via R package gtheory", "gtheory.htm",
        "mfrm_generalizability(); mfrm_d_study(); compute_facet_icc()", "supported_with_caveat",
        "Observed univariate G-study variance components plus D-study projections with residual-scaling sensitivity and `IdentificationStatus` columns.",
        "Package-native caveated G/D-study route; not a FACETS/gtheory UI clone, not multivariate/profile G-theory, and not suitable for high-stakes use when boundary or singular fits are reported.",
        package_native_alternative),
    row("R/Web plots", "Connectivity network graph via igraph", "networkgraph.htm",
        "subset_connectivity_report(); mfrm_network_analysis(); rater_network_analysis(); rater_halo_network_analysis(); plot(..., type = \"network\")", "implemented",
        "Facet-level co-observation network plus rater agreement/disagreement/severity-direction and halo networks with reusable node/edge tables.",
        "R-native igraph analysis and display rather than FACETS menu output.", package_native_alternative),
    row("Output file", "Specification settings file", "specificationfile.htm",
        "build_mfrm_manifest(); build_mfrm_replay_script()", "partial",
        "R-native reproducibility manifest and replay script.",
        "Does not write a FACETS command specification file.", package_native_alternative),
    row("Output file", "Anchor output file", "anchorfile.htm",
        "make_anchor_table(); export_mfrm_bundle(include = \"anchors\")", "implemented",
        "Reusable anchor tables from fitted estimates.",
        "Uses R/CSV tables rather than FACETS fixed syntax.", package_native_alternative),
    row("Output file", "Graph plotting file", "graphoutputfile.htm",
        "facets_output_file_bundle(include = \"graph\")", "implemented",
        "Graphfile-style category curve output.",
        "Command-level FACETS graph options are not fully mirrored.", package_native_alternative),
    row("Output file", "Output report file", "outputfile.htm",
        "export_summary_appendix(); build_fixed_reports()", "partial",
        "Structured appendix/report artifacts.",
        "Full FACETS report-file emulation is not implemented.", external_format_alternative),
    row("Output file", "Residuals output file", "residualfile.htm",
        "write_mfrm_residual_file(); diagnose_mfrm(); unexpected_response_table(); residual plot helpers", "implemented",
        "Standalone observation-level residual CSV/TSV output, residual tables, and residual visualizations are available.",
        "Uses package-native residual columns rather than exact FACETS fixed-field residual syntax.", package_native_alternative),
    row("Output file", "Score output file", "scorefile.htm",
        "facets_output_file_bundle(include = \"score\"); read_facets_fit_table()",
        "partial",
        "Score-side export/import is available for documented Rasch-family routes covered by package tests.",
        "Bounded GPCM score-side equivalence is outside the documented scope.", package_native_alternative),
    row("Output file", "Simulated data file", "simulatedfile.htm",
        "simulate_mfrm_data(); build_mfrm_sim_spec()", "partial",
        "Simulation data and explicit simulation specifications.",
        "Not a FACETS simulated-data file clone.", package_native_alternative),
    row("Output file", "Subset group-anchor file", "subsetfile.htm",
        "write_mfrm_subset_file(); group_anchors; review_mfrm_anchors(); make_anchor_table()", "partial",
        "Connected-subset summary/node files and group-anchor inputs/checks are available.",
        "The standalone subset writer exports connectivity review tables, not a full FACETS UI-compatible subset command file.", package_native_alternative),
    row("Output file", "Winsteps control and data file", "winstepsfile.htm",
        "none", "not_implemented",
        "No Winsteps control/data export route.",
        "Would require a separate Winsteps output contract.", external_program_alternative),
    row("Graph menu", "Category probability curves", "graphs.htm",
        "category_curves_report(); plot(fit, type = \"ccc\")", "implemented",
        "Category probability curve data and plots.",
        "R-native plots replace FACETS graph menu output.", package_native_alternative),
    row("Graph menu", "Expected score ICC/IRF", "graphs.htm",
        "plot(fit, type = \"pathway\"); category_curves_report()", "implemented",
        "Expected-score curves over theta.",
        "Not labeled as FACETS ICC/IRF menu output.", package_native_alternative),
    row("Graph menu", "Cumulative probability curves", "graphs.htm",
        "category_curves_report(); plot(..., type = \"cumulative\")", "implemented",
        "Cumulative category-probability curve data, flipped direction data, and approximate .5 boundaries are available.",
        "R-native plot data replace FACETS graph-menu output.", package_native_alternative),
    row("Graph menu", "Test information function", "graphs.htm",
        "compute_information(); plot_information(type = \"tif\")", "implemented",
        "Design-weighted test/scale information curves.",
        "R-native information definition and plot data.", package_native_alternative),
    row("Graph menu", "Category information function", "graphs.htm",
        "category_curves_report(); plot(..., type = \"category_information\"); compute_information(); plot_information(type = \"iif\")", "implemented",
        "Category-specific information contributions, total information curves, and facet/level contribution curves are available.",
        "R-native plot data replace FACETS graph-menu output.", package_native_alternative),
    row("Graph menu", "Conditional probability curves", "graphs.htm",
        "category_curves_report()", "partial",
        "Category probability curves conditional on theta are available.",
        "FACETS conditional-probability menu semantics are not mirrored exactly.", external_format_alternative),
    row("Specification/workflow", "Full FACETS command-file parser and UI option surface",
        "index.htm",
        "run_mfrm_facets(); fit_mfrm()", "not_targeted",
        "R function arguments are the package interface.",
        "Parsing arbitrary FACETS command files is outside the package scope.", external_program_alternative),
    row("Specification/workflow", "Exact FACETS line-printer report emulation",
        "outputtableindex.htm",
        "build_fixed_reports() for selected tables", "not_targeted",
        "Selected fixed-width handoff is available.",
        "Exact full report emulation is outside the package scope.", external_program_alternative),
    row("Specification/workflow", "Raw FACETS report-text import",
        "outputtableindex.htm",
        "read_facets_fit_table() for delimited/fixed-field score extracts", "partial",
        "Fit/score table import is supported.",
        "General raw FACETS report parsing is not implemented.", external_format_alternative),
    row("Current scope boundary",
        "Versioned frozen-calibration import and operational scoring",
        "mfrmr 0.2.3 public contract",
        "none in 0.2.3", "not_implemented",
        "No current public route imports a reusable versioned calibration bundle for operational scoring.",
        "Posterior scoring from an existing fitted object is supported separately and is not a reusable frozen-calibration contract.",
        "Use fitted-object posterior scoring for current analyses; retain frozen-calibration workflows for a later release."),
    row("Current scope boundary",
        "General threshold or step anchors and starting-value import",
        "mfrmr 0.2.3 public contract",
        "none in 0.2.3", "not_implemented",
        "No current public route accepts general threshold or step anchors or a threshold starting-value contract.",
        "Element and group anchors do not make threshold ladders fixed or supply a general calibration-import schema.",
        "Use current element/group anchor routes only for their documented scope; retain threshold anchoring for a later release."),
    row("Current scope boundary",
        "Multiple observed scales and scale-specific PCM",
        "mfrmr 0.2.3 public contract",
        "none in 0.2.3", "not_implemented",
        "Each fit uses one observed score scale and one homogeneous response-model family.",
        "There is no per-observation ScaleId contract, scale-specific category map, or ragged scale-specific PCM threshold block.",
        "Fit supported single-scale designs separately; retain multi-scale and mixed-family workflows for a later release."),
    row("Current scope boundary",
        "Nominal/multinomial response models",
        "mfrmr 0.2.3 public contract",
        "none in 0.2.3", "not_implemented",
        "The current RSM, PCM, and bounded-GPCM routes model ordered category probabilities only.",
        "A category-probability vector that sums to one is not an unordered nominal-response or multinomial-logit model; category order enters every current likelihood.",
        "Use a nominal-response or multinomial-regression implementation externally when category order is not substantively defined."),
    row("Current scope boundary",
        "Binomial-trial and Poisson/count response models",
        "models.htm",
        "none in 0.2.3", "not_implemented",
        "Binary ordered scores are available as the two-category special case of the current ordered-response kernel.",
        "Grouped binomial trials, Poisson counts, negative-binomial counts, and other count likelihoods are not implemented; integer scores are interpreted as ordered category codes.",
        "Use FACETS or another count-model implementation for an appropriate binomial-trial or Poisson estimand; do not relabel an ordered-category fit as a count model."),
    row("Observation weighting",
        "Row-frequency weights for ordered ratings",
        "mfrmr 0.2.3 weight contract",
        "fit_mfrm(weight = ...)", "supported_with_caveat",
        "A positive numeric observation weight multiplies that row's conditional ordered-category likelihood contribution and can represent a defensible row-replication weight.",
        "It is not a general collapsed-person frequency table: for MML, powering responses inside one Person pattern is not equivalent to replicating a complete Person pattern after marginalization. It also does not create a count-response family, model within-cell dependence, or make non-unit-weight fits eligible for the common information-criterion panel.",
        "Retain one distinguishable event per row when possible; preserve distinct Person response patterns, and report the exact likelihood-weight interpretation."),
    row("Current scope boundary",
        "Native multidimensional estimation and dimension-specific scores",
        "mfrmr 0.2.3 public contract",
        "none in 0.2.3", "not_implemented",
        "The current public estimator and score routes are unidimensional.",
        "Residual PCA is exploratory dimensionality evidence, not native multidimensional estimation or dimension-specific scoring.",
        "Use exploratory dimensionality diagnostics and external multidimensional software when a multidimensional estimator is required."),
    row("Current scope boundary",
        "Unrestricted GPCM",
        "mfrmr 0.2.3 public contract",
        "none for unrestricted GPCM in 0.2.3", "not_implemented",
        "The current public GPCM route is bounded and requires slope_facet to equal step_facet.",
        "Bounded GPCM support does not establish an unrestricted free-discrimination model family.",
        "Use gpcm_capability_matrix() and the documented bounded-GPCM route; retain unrestricted GPCM for a later release."),
    row("Current scope boundary",
        "FACETS free-slope polytomous GPCM comparison",
        "t7menu.htm",
        "none as a direct common-estimand lane", "not_implemented",
        "FACETS PCM/JMLE can serve as the direct equal-discrimination comparison after the full estimation contract is aligned.",
        "FACETS Table 7 Estimated Discrimination is a post-fit diagnostic that does not update other Rasch estimates, so it is not the jointly estimated bounded-GPCM slope from mfrmr.",
        "Use FACETS for the matched PCM/JML lane or as a deliberately misspecified equal-discrimination control; use a genuinely slope-estimating program only after the GPCM kernel and identification are matched.")
  ))

  row.names(out) <- NULL
  surface_by_status <- c(
    implemented = "available",
    supported_with_caveat = "available_with_caveat",
    partial = "partial",
    not_implemented = "unavailable",
    not_targeted = "out_of_scope"
  )
  contract_by_status <- c(
    implemented = "package_native_not_facets_equivalent",
    supported_with_caveat = "package_native_caveated",
    partial = "partial_package_native_contract",
    not_implemented = "not_available",
    not_targeted = "not_applicable"
  )
  validation_by_status <- c(
    implemented = "external_match_not_established",
    supported_with_caveat = "external_match_not_established",
    partial = "external_match_not_established",
    not_implemented = "not_applicable",
    not_targeted = "not_applicable"
  )
  operation_by_status <- c(
    implemented = "package_route_available",
    supported_with_caveat = "available_with_mandatory_caveat",
    partial = "partial_workflow_only",
    not_implemented = "blocked",
    not_targeted = "external_only"
  )
  out$SurfaceCoverage <- unname(surface_by_status[out$Status])
  out$StatisticalContract <- unname(contract_by_status[out$Status])
  out$ValidationEvidence <- unname(validation_by_status[out$Status])
  out$OperationalStatus <- unname(operation_by_status[out$Status])
  wright_visual <- grepl(
    "Table 6.0: all-facet Wright map rulers",
    out$FACETSFeature,
    fixed = TRUE
  )
  out$SurfaceCoverage[
    wright_visual & out$Status == "implemented"
  ] <- "familiar_visual_grammar_available"
  out <- out[, c(
    "FACETSArea", "FACETSFeature", "FACETSReference", "mfrmrRoute", "Status",
    "SurfaceCoverage", "StatisticalContract", "ValidationEvidence",
    "OperationalStatus", "Capability", "Limitation", "Alternative"
  ), drop = FALSE]
  if (identical(status, "all")) {
    return(out)
  }
  out[out$Status == status, , drop = FALSE]
}
