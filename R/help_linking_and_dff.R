#' mfrmr Linking and DFF Guide
#'
#' @description
#' Package-native guide to checking connectedness, building anchor-based links,
#' monitoring drift, and screening differential facet functioning (DFF) in
#' `mfrmr`.
#'
#' @section Start with the linking question:
#' - "Is the design connected enough to support a common scale?"
#'   Use [subset_connectivity_report()] and `plot(..., type = "design_matrix")`.
#' - "Which elements can I export as anchors from an existing fit?"
#'   Use [make_anchor_table()] and [audit_mfrm_anchors()].
#' - "How do I anchor a new administration to a baseline?"
#'   Use [anchor_to_baseline()].
#' - "Have common elements drifted across separately fitted waves?"
#'   Use [detect_anchor_drift()] and [plot_anchor_drift()].
#' - "Do specific facet levels function differently across groups?"
#'   Use [analyze_dff()] and [plot_dif_heatmap()].
#'
#' @section Recommended linking route:
#' 1. Fit with [fit_mfrm()] and diagnose with [diagnose_mfrm()].
#' 2. Check connectedness with [subset_connectivity_report()].
#' 3. Build or audit anchors with [make_anchor_table()] and
#'    [audit_mfrm_anchors()].
#' 4. Use [anchor_to_baseline()] when you need to place raw new data onto a
#'    baseline scale.
#' 5. Use [build_equating_chain()] only as a screened linking aid across
#'    already fitted waves.
#' 6. Use [detect_anchor_drift()] for stability monitoring on separately fitted
#'    waves.
#' 7. Run [analyze_dff()] only after checking connectivity and common-scale
#'    evidence.
#'
#' @section Which helper answers which task:
#' \describe{
#'   \item{[subset_connectivity_report()]}{Summarizes connected subsets,
#'   bottleneck facets, and design-matrix coverage.}
#'   \item{[make_anchor_table()]}{Extracts reusable anchor candidates from a fit.}
#'   \item{[anchor_to_baseline()]}{Anchors new raw data to a baseline fit and
#'   returns anchored diagnostics plus a consistency check against the baseline
#'   scale.}
#'   \item{[detect_anchor_drift()]}{Compares fitted waves directly to flag
#'   unstable anchor elements.}
#'   \item{[build_equating_chain()]}{Accumulates screened pairwise links across
#'   a series of administrations or forms.}
#'   \item{[analyze_dff()]}{Screens differential facet functioning with residual
#'   or refit methods, using screening-only language unless linking and
#'   precision support stronger interpretation.}
#' }
#'
#' @section Practical linking rules:
#' - Check connectedness before interpreting subgroup or wave differences.
#' - Use DFF outputs as screening results when common-scale linking is weak.
#' - Treat drift flags as prompts for review, not automatic evidence that an
#'   anchor must be removed.
#' - Treat `LinkSupportAdequate = FALSE` as a weak-link warning: at least one
#'   linking facet retained fewer than 5 common elements after screening.
#' - Rebuild anchors from a defensible baseline rather than chaining unstable
#'   links by hand.
#'
#' @section Typical workflow:
#' - Cross-sectional linkage review:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [subset_connectivity_report()] ->
#'   `plot(..., type = "design_matrix")`.
#' - Baseline placement review:
#'   [make_anchor_table()] -> [anchor_to_baseline()] -> [diagnose_mfrm()].
#' - Multi-wave drift review:
#'   fit each wave separately -> [detect_anchor_drift()] -> [plot_anchor_drift()].
#' - Group comparison route:
#'   [subset_connectivity_report()] -> [analyze_dff()] ->
#'   [dif_report()] -> [plot_dif_heatmap()].
#'
#' @section Companion guides:
#' - For visual follow-up, see [mfrmr_visual_diagnostics].
#' - For report/table selection, see [mfrmr_reports_and_tables].
#' - For end-to-end routes, see [mfrmr_workflow_methods].
#' - For a longer walkthrough, see
#'   `vignette("mfrmr-linking-and-dff", package = "mfrmr")`.
#'
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 10
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#'
#' subsets <- subset_connectivity_report(fit, diagnostics = diag)
#' subsets$summary[, c("Subset", "Observations", "ObservationPercent")]
#'
#' dff <- analyze_dff(fit, diag, facet = "Rater", group = "Group", data = toy)
#' head(dff$dif_table[, c("Level", "Group1", "Group2", "Classification")])
#'
#' @name mfrmr_linking_and_dff
NULL
