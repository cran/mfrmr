#' FACETS Positioning Guide
#'
#' @description
#' `facets_positioning_guide()` gives user-facing wording for the relationship
#' between `mfrmr` and FACETS. Use it when a report, migration note, or
#' methods appendix must make clear that `mfrmr` is not a FACETS numerical
#' clone.
#'
#' @details
#' The guide separates four ideas that are easy to conflate:
#'
#' - estimation authority: fitted values come from `mfrmr` unless external
#'   FACETS output is explicitly supplied;
#' - compatibility purpose: FACETS-style names and files are transition,
#'   handoff, and report-organization surfaces;
#' - external comparison: FACETS comparisons require a supplied external table
#'   and should separate MnSq differences from df/ZSTD convention differences;
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
      "Reporting source of truth",
      "Extension beyond FACETS"
    ),
    Position = c(
      "mfrmr estimates are package-native; FACETS-style names do not mean that FACETS estimated the model.",
      "FACETS-style wrappers, table labels, and files support transition, handoff, and report organization, not optimizer-level reproduction.",
      "Numerical comparison requires an explicit external FACETS output table supplied by the user.",
      "Inference and reporting should be based on native fit, diagnostics, review, table, and plot-data objects.",
      "GPCM, D-study, network, and reusable visualization data are extension routes rather than FACETS menu clones."
    ),
    RecommendedWording = c(
      "The model was estimated with mfrmr; FACETS-style output names are used only to organize the report.",
      "FACETS-style outputs were generated for handoff or reader familiarity; they are not evidence of FACETS numerical equivalence.",
      "When external FACETS output is supplied, compare MnSq first and report df/ZSTD convention sensitivity separately.",
      "Report estimates, standard errors, fit summaries, and plots from documented mfrmr objects.",
      "Use package-native extensions as additional evidence and label them as mfrmr analyses."
    ),
    PrimaryRoute = c(
      "fit_mfrm(); diagnose_mfrm(); reporting_checklist()",
      "facets_feature_coverage(); run_mfrm_facets(); facets_output_file_bundle()",
      "read_facets_fit_table(); facets_fit_review(); fit_measures_table(df_sensitivity = TRUE)",
      "build_summary_table_bundle(); build_visual_summaries(); plot_data()",
      "gpcm_capability_matrix(); mfrm_d_study(); mfrm_network_analysis(); plot_data_components()"
    ),
    stringsAsFactors = FALSE
  )
}

#' FACETS Feature Coverage Matrix
#'
#' @description
#' `facets_feature_coverage()` summarizes how the current `mfrmr` release maps
#' the main FACETS output-table, output-file, and graph-menu surface to package
#' functions.
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
#' Status meanings:
#'
#' - `implemented`: a package-native route covers the substantive output.
#' - `partial`: the concept is covered, but not the full FACETS formatting,
#'   option surface, file type, or external integration.
#' - `not_implemented`: a FACETS feature has no direct package-native route in
#'   the current release.
#' - `not_targeted`: the feature is tied to FACETS UI, Web/Excel handoff, or
#'   another external program format and is not a release goal.
#'
#' @return A data.frame with columns:
#' - `FACETSArea`
#' - `FACETSFeature`
#' - `FACETSReference`
#' - `mfrmrRoute`
#' - `Status`
#' - `Scope`
#' - `GapOrBoundary`
#' - `Priority`
#'
#' @references
#' Linacre, J. M. (2026). *A user's guide to FACETS, version 4.5.0*.
#' Output tables - files - plots - graphs:
#' <https://www.winsteps.com/facetman64/outputtableindex.htm>.
#'
#' @seealso [facets_positioning_guide()], [mfrmr_output_guide()],
#'   [facets_fit_df_guide()], [read_facets_fit_table()], [facets_fit_review()],
#'   [gpcm_capability_matrix()]
#' @examples
#' facets_feature_coverage()
#' facets_feature_coverage("partial")
#' facets_feature_coverage("not_implemented")
#' @export
facets_feature_coverage <- function(status = c("all", "implemented", "partial",
                                               "not_implemented", "not_targeted")) {
  status <- match.arg(status)

  row <- function(area, feature, reference, route, status, scope, gap, priority) {
    data.frame(
      FACETSArea = area,
      FACETSFeature = feature,
      FACETSReference = reference,
      mfrmrRoute = route,
      Status = status,
      Scope = scope,
      GapOrBoundary = gap,
      Priority = priority,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, list(
    row("Output table", "Table 1: specification summary", "table1.htm",
        "specifications_report()", "implemented",
        "Structured run settings and reproducibility context.",
        "Not an exact FACETS line-printer layout.", "release_core"),
    row("Output table", "Table 2: data summary report", "table2.htm",
        "data_quality_report(); describe_mfrm_data()", "implemented",
        "Rows, exclusions, missingness, score support, and response-pattern QC.",
        "Structured QC replaces FACETS text layout.", "release_core"),
    row("Output table", "Table 3: main iteration report", "table3.htm",
        "estimation_iteration_report()", "partial",
        "Convergence and replayed iteration evidence.",
        "Does not reproduce every FACETS optimizer-internal line.", "release_core"),
    row("Output table", "Table 4: unexpected responses", "table4.htm",
        "unexpected_response_table(); plot_unexpected()", "implemented",
        "Case-level unexpected-response screening.",
        "Structured table and plots, not printer-identical FACETS output.", "release_core"),
    row("Output table", "Table 5: measurable data summary", "table5.htm",
        "measurable_summary_table(); describe_mfrm_data()", "implemented",
        "Facet coverage, category counts, and subset/connectivity checks.",
        "Column order and text layout differ from FACETS.", "release_core"),
    row("Output table", "Table 6.0: all-facet Wright map rulers", "table6.htm",
        "plot(fit, type = \"wright\"); plot_wright_unified()", "implemented",
        "Common-logit person/facet/threshold display.",
        "R-native graphics replace FACETS ruler text.", "release_core"),
    row("Output table", "Table 6.0.0: disjoint element listing", "table6_0_0.htm",
        "subset_connectivity_report()", "implemented",
        "Disconnected subsets and facet-by-subset coverage.",
        "Network-style graph is not the default display.", "release_core"),
    row("Output table", "Table 6.2: graphical facet statistics", "table6_2.htm",
        "facet_statistics_report(); plot(...)", "partial",
        "Facet statistics and visual summaries.",
        "FACETS M/S/Q/X printer-graph formatting is not reproduced exactly.", "release_core"),
    row("Output table", "Table 7: facet measurement report", "table7.htm",
        "fit_measures_table(); diagnose_mfrm(); summary(fit)", "implemented",
        "Measures, SEs, fit, anchoring status, and review flags.",
        "FACETS column order/options are broader than the default table.", "release_core"),
    row("Output table", "Table 7: reliability and chi-square", "table7summarystatistics.htm",
        "facets_chisq_table(); diagnose_mfrm()$reliability", "implemented",
        "Rasch/FACETS-style separation, reliability, and chi-square summaries.",
        "Uses package-native structured output.", "release_core"),
    row("Output table", "Table 7: agreement statistics", "table7agreementstatistics.htm",
        "interrater_agreement_table(); rater_network_analysis(); rater_halo_network_analysis(); plot_interrater_agreement()", "implemented",
        "Observed/expected agreement, pairwise rater-network, rater-by-criterion halo network, and rater-agreement views.",
        "Structured output replaces FACETS text blocks.", "release_core"),
    row("Output table", "Table 8.1: dichotomous/binomial/Poisson statistics",
        "table8_1dichotomous.htm",
        "rating_scale_table() for two-category ordered scores", "partial",
        "Two-category Rasch-category summaries are available.",
        "FACETS binomial-trial and Poisson-specific reports are not implemented.", "defer"),
    row("Output table", "Table 8.1: polytomous rating-scale/partial-credit statistics",
        "table8_1ratingscale.htm",
        "rating_scale_table(); category_structure_report()", "implemented",
        "Rating-scale/partial-credit category diagnostics and thresholds.",
        "Exact FACETS text layout is not reproduced.", "release_core"),
    row("Output table", "Table 8: scale-structure bar chart", "table8barchart.htm",
        "category_structure_report()", "partial",
        "Category structure and transition summaries.",
        "FACETS line-printer artwork is not reproduced exactly.", "release_core"),
    row("Output table", "Table 8: scale-structure probability curves", "table8curves.htm",
        "category_curves_report(); plot(fit, type = \"ccc\")", "implemented",
        "Category probability and expected-score curve data.",
        "Uses R-native plot data rather than FACETS graph text.", "release_core"),
    row("Output table", "Table 9: bias-estimation iteration report", "table9.htm",
        "estimate_bias(); bias_iteration_report()", "implemented",
        "Bias recalibration path and final iteration status.",
        "Conditional screening semantics are documented separately.", "release_core"),
    row("Output table", "Table 10: unexpected after allowing for bias", "table10.htm",
        "unexpected_after_bias_table()", "implemented",
        "Unexpected rows after the current bias-screening layer.",
        "Structured table replaces FACETS text layout.", "release_core"),
    row("Output table", "Table 11: bias-calculation counts", "table11.htm",
        "bias_count_table()", "implemented",
        "Response counts behind bias estimates.",
        "Structured output replaces FACETS text layout.", "release_core"),
    row("Output table", "Table 12: bias summary report", "table12.htm",
        "summary(estimate_bias(...)); plot_bias_interaction()", "partial",
        "Distributional and visual bias summaries are available.",
        "FACETS vertical frequency bar-chart is not reproduced exactly.", "release_core"),
    row("Output table", "Table 13: DIF/bias detail report", "table13.htm",
        "estimate_bias(); bias_interaction_report()", "implemented",
        "Ranked cell-level bias/interactions with screening statistics.",
        "Reported as screening evidence, not final fairness inference.", "release_core"),
    row("Output table", "Table 14: pairwise bias report", "table14.htm",
        "bias_pairwise_report(); build_fixed_reports()", "implemented",
        "Pairwise contrasts for two-way bias runs.",
        "Higher-order runs omit pairwise sections by design.", "release_core"),
    row("Output table", "DIF/bias Excel plot", "difbiasplot.htm",
        "plot_bias_interaction(plot = ...)", "partial",
        "R-native scatter, heatmap, and facet-profile bias displays.",
        "Excel-specific output is not implemented.", "defer"),
    row("R/Web plots", "Scatterplots and histograms from FACETS menus",
        "outputtableindex.htm",
        "plot_data(); package plot helpers", "partial",
        "Reusable plot data supports custom R graphics.",
        "FACETS arbitrary R/Web plotting menus are not mirrored.", "low"),
    row("R/Web plots", "X-Y plot: R Statistics", "xyplotr.htm",
        "plot_data(); user-defined R plotting", "partial",
        "Users can build X-Y plots from returned data frames.",
        "No dedicated FACETS-style arbitrary X-Y plot wrapper.", "low"),
    row("R/Web plots", "X-Y plot: Webpage", "xyplotwebpage.htm",
        "none", "not_targeted",
        "No package-native Webpage plot generator.",
        "Webpage menu output is a FACETS UI feature.", "not_planned"),
    row("R/Web plots", "X-Y-Z plot: R Statistics", "xyzplotr.htm",
        "plot(fit, type = \"ccc_surface\"); plot_data()", "partial",
        "Selected 3D/surface-ready plot data are available.",
        "No arbitrary FACETS X-Y-Z plot wrapper.", "low"),
    row("R/Web plots", "Histogram: R Statistics", "histogramr.htm",
        "plot_data(); plot(fit, type = \"wright\"); plot_qc_dashboard()", "partial",
        "Several package outputs include histogram-like summaries.",
        "No general FACETS histogram menu clone.", "low"),
    row("R/Web plots", "Generalizability Theory via R package gtheory", "gtheory.htm",
        "mfrm_generalizability(); mfrm_d_study(); compute_facet_icc()", "implemented",
        "Observed G-study variance components plus D-study projections with residual-scaling sensitivity.",
        "Package-native G/D-study route; not a FACETS/gtheory UI clone.", "release_core"),
    row("R/Web plots", "Connectivity network graph via igraph", "networkgraph.htm",
        "subset_connectivity_report(); mfrm_network_analysis(); rater_network_analysis(); rater_halo_network_analysis(); plot(..., type = \"network\")", "implemented",
        "Facet-level co-observation network plus rater agreement/disagreement/severity-direction and halo networks with reusable node/edge tables.",
        "R-native igraph analysis and display rather than FACETS menu output.", "release_core"),
    row("Output file", "Specification settings file", "specificationfile.htm",
        "build_mfrm_manifest(); build_mfrm_replay_script()", "partial",
        "R-native reproducibility manifest and replay script.",
        "Does not write a FACETS command specification file.", "release_core"),
    row("Output file", "Anchor output file", "anchorfile.htm",
        "make_anchor_table(); export_mfrm_bundle(include = \"anchors\")", "implemented",
        "Reusable anchor tables from fitted estimates.",
        "Uses R/CSV tables rather than FACETS fixed syntax.", "release_core"),
    row("Output file", "Graph plotting file", "graphoutputfile.htm",
        "facets_output_file_bundle(include = \"graph\")", "implemented",
        "Graphfile-style category curve output.",
        "Command-level FACETS graph options are not fully mirrored.", "release_core"),
    row("Output file", "Output report file", "outputfile.htm",
        "export_summary_appendix(); build_fixed_reports()", "partial",
        "Structured appendix/report artifacts.",
        "Full FACETS report-file emulation is not implemented.", "defer"),
    row("Output file", "Residuals output file", "residualfile.htm",
        "write_mfrm_residual_file(); diagnose_mfrm(); unexpected_response_table(); residual plot helpers", "implemented",
        "Standalone observation-level residual CSV/TSV output, residual tables, and residual visualizations are available.",
        "Uses package-native residual columns rather than exact FACETS fixed-field residual syntax.", "release_core"),
    row("Output file", "Score output file", "scorefile.htm",
        "facets_output_file_bundle(include = \"score\"); read_facets_fit_table()",
        "partial",
        "Score-side export/import is available for validated Rasch-family routes.",
        "Bounded GPCM score-side equivalence is outside the current boundary.", "release_core"),
    row("Output file", "Simulated data file", "simulatedfile.htm",
        "simulate_mfrm_data(); build_mfrm_sim_spec()", "partial",
        "Simulation data and explicit simulation specifications.",
        "Not a FACETS simulated-data file clone.", "release_core"),
    row("Output file", "Subset group-anchor file", "subsetfile.htm",
        "write_mfrm_subset_file(); group_anchors; review_mfrm_anchors(); make_anchor_table()", "partial",
        "Connected-subset summary/node files and group-anchor inputs/checks are available.",
        "The standalone subset writer exports connectivity review tables, not a full FACETS UI-compatible subset command file.", "release_core"),
    row("Output file", "Winsteps control and data file", "winstepsfile.htm",
        "none", "not_implemented",
        "No Winsteps control/data export route.",
        "Would require a separate Winsteps output contract.", "not_planned"),
    row("Graph menu", "Category probability curves", "graphs.htm",
        "category_curves_report(); plot(fit, type = \"ccc\")", "implemented",
        "Category probability curve data and plots.",
        "R-native plots replace FACETS graph menu output.", "release_core"),
    row("Graph menu", "Expected score ICC/IRF", "graphs.htm",
        "plot(fit, type = \"pathway\"); category_curves_report()", "implemented",
        "Expected-score curves over theta.",
        "Not labeled as FACETS ICC/IRF menu output.", "release_core"),
    row("Graph menu", "Cumulative probability curves", "graphs.htm",
        "category_curves_report(); plot(..., type = \"cumulative\")", "implemented",
        "Cumulative category-probability curve data, flipped direction data, and approximate .5 boundaries are available.",
        "R-native plot data replace FACETS graph-menu output.", "release_core"),
    row("Graph menu", "Test information function", "graphs.htm",
        "compute_information(); plot_information(type = \"tif\")", "implemented",
        "Design-weighted test/scale information curves.",
        "R-native information definition and plot data.", "release_core"),
    row("Graph menu", "Category information function", "graphs.htm",
        "category_curves_report(); plot(..., type = \"category_information\"); compute_information(); plot_information(type = \"iif\")", "implemented",
        "Category-specific information contributions, total information curves, and facet/level contribution curves are available.",
        "R-native plot data replace FACETS graph-menu output.", "release_core"),
    row("Graph menu", "Conditional probability curves", "graphs.htm",
        "category_curves_report()", "partial",
        "Category probability curves conditional on theta are available.",
        "FACETS conditional-probability menu semantics are not mirrored exactly.", "defer"),
    row("Specification/workflow", "Full FACETS command-file parser and UI option surface",
        "index.htm",
        "run_mfrm_facets(); fit_mfrm()", "not_targeted",
        "R function arguments are the package interface.",
        "Parsing arbitrary FACETS command files is outside the release scope.", "not_planned"),
    row("Specification/workflow", "Exact FACETS line-printer report emulation",
        "outputtableindex.htm",
        "build_fixed_reports() for selected tables", "not_targeted",
        "Selected fixed-width handoff is available.",
        "Exact full report emulation is intentionally not a package goal.", "not_planned"),
    row("Specification/workflow", "Raw FACETS report-text import",
        "outputtableindex.htm",
        "read_facets_fit_table() for delimited/fixed-field score extracts", "partial",
        "Fit/score table import is supported.",
        "General raw FACETS report parsing is not implemented.", "defer")
  ))

  row.names(out) <- NULL
  if (identical(status, "all")) {
    return(out)
  }
  out[out$Status == status, , drop = FALSE]
}
