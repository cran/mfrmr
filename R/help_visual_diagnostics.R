#' mfrmr Visual Diagnostics Map
#'
#' @description
#' Quick guide to choosing the right base-R diagnostic plot in `mfrmr`.
#' Use this page when you know the analysis question but do not yet know
#' which plotting helper or `plot()` method to call.
#'
#' @section Start with the question:
#' - "Do persons and facet levels overlap on the same logit scale?"
#'   Use `plot(fit, type = "wright")` or [plot_wright_unified()].
#' - "Where do score categories transition across theta?"
#'   Use `plot(fit, type = "pathway")` and `plot(fit, type = "ccc")`.
#' - "Is the design linked well enough across subsets or administrations?"
#'   Use `plot(subset_connectivity_report(...), type = "design_matrix")` and
#'   [plot_anchor_drift()].
#' - "Which responses or levels look locally problematic?"
#'   Use [plot_unexpected()] and [plot_displacement()].
#' - "Do raters agree and do facets separate meaningfully?"
#'   Use [plot_interrater_agreement()] and [plot_facets_chisq()].
#' - "Is there notable residual structure after the main Rasch dimension?"
#'   Use [plot_residual_pca()].
#' - "Which interaction cells or facet levels drive bias screening results?"
#'   Use [plot_bias_interaction()].
#' - "I need one compact triage screen first."
#'   Use [plot_qc_dashboard()].
#'
#' @section Recommended visual route:
#' 1. Start with [plot_qc_dashboard()] for one-page triage.
#' 2. Move to [plot_unexpected()], [plot_displacement()], and
#'    [plot_interrater_agreement()] for flagged local issues.
#' 3. Use `plot(fit, type = "wright")`, `plot(fit, type = "pathway")`,
#'    and `plot_residual_pca()` for structural interpretation.
#' 4. Use `plot(..., draw = FALSE)` when you want reusable plotting payloads
#'    instead of immediate graphics.
#' 5. Use `preset = "publication"` when you want the package's cleaner
#'    manuscript-oriented styling.
#'
#' @section Which plot answers which question:
#' \describe{
#'   \item{`plot(fit, type = "wright")`}{Shared logit map of persons, facet
#'   levels, and step thresholds. Best for targeting and spread.}
#'   \item{`plot(fit, type = "pathway")`}{Expected score by theta, with
#'   dominant-category strips. Best for scale progression.}
#'   \item{`plot(fit, type = "ccc")`}{Category probability curves. Best for
#'   checking whether categories peak in sequence.}
#'   \item{`plot_unexpected()`}{Observation-level surprises. Best for case
#'   review and local misfit triage.}
#'   \item{`plot_displacement()`}{Level-wise anchor movement. Best for anchor
#'   robustness and residual calibration tension.}
#'   \item{`plot_interrater_agreement()`}{Exact agreement, expected agreement,
#'   pairwise correlation, and agreement gaps. Best for rater consistency.}
#'   \item{`plot_facets_chisq()`}{Facet variability and chi-square summaries.
#'   Best for checking whether a facet contributes meaningful spread.}
#'   \item{`plot_residual_pca()`}{Residual structure after the Rasch dimension
#'   is removed. Best for exploratory residual-structure review, not as a
#'   standalone unidimensionality test.}
#'   \item{`plot_bias_interaction()`}{Interaction-bias screening views for
#'   cells and facet profiles. Best for systematic departure from the
#'   additive main-effects model.}
#'   \item{`plot_anchor_drift()`}{Anchor drift and screened linking-chain visuals.
#'   Best for multi-form or multi-wave linking review after checking retained
#'   common-element support.}
#' }
#'
#' @section Practical interpretation rules:
#' - Wright map: look for gaps between person density and facet/step
#'   locations; large gaps indicate weaker targeting.
#' - Pathway / CCC: look for monotone progression and clear category
#'   dominance bands; flat or overlapping curves suggest weak category
#'   separation.
#' - Unexpected / displacement: use as screening tools, not final evidence
#'   by themselves.
#' - Inter-rater agreement and facet variability address different questions:
#'   agreement concerns scoring consistency, whereas variability concerns
#'   whether facet elements are statistically distinguishable.
#' - Residual PCA and bias plots should be interpreted as follow-up layers
#'   after the main fit screen, not as first-pass diagnostics.
#'
#' @section Typical workflow:
#' - Quick screening:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [plot_qc_dashboard()].
#' - Scale and targeting review:
#'   `plot(fit, type = "wright")` -> `plot(fit, type = "pathway")` ->
#'   `plot(fit, type = "ccc")`.
#' - Linking review:
#'   [subset_connectivity_report()] -> `plot(..., type = "design_matrix")` ->
#'   [plot_anchor_drift()].
#' - Interaction review:
#'   [estimate_bias()] -> [plot_bias_interaction()] ->
#'   [reporting_checklist()].
#'
#' @section Companion vignette:
#' For a longer, plot-first walkthrough, run
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @seealso [mfrmr_workflow_methods], [mfrmr_reports_and_tables],
#'   [mfrmr_reporting_and_apa], [mfrmr_linking_and_dff],
#'   [plot.mfrm_fit()], [plot_qc_dashboard()],
#'   [plot_unexpected()], [plot_displacement()], [plot_interrater_agreement()],
#'   [plot_facets_chisq()], [plot_residual_pca()], [plot_bias_interaction()],
#'   [plot_anchor_drift()]
#'
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 25
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#'
#' qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE, preset = "publication")
#' names(qc$data)
#'
#' wright <- plot(fit, type = "wright", draw = FALSE, preset = "publication")
#' wright$data$preset
#'
#' pca <- analyze_residual_pca(diag, mode = "overall")
#' scree <- plot_residual_pca(pca, plot_type = "scree", draw = FALSE, preset = "publication")
#' scree$data$preset
#' }
#'
#' @name mfrmr_visual_diagnostics
NULL
