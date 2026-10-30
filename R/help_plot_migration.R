#' Plotting after the 0.2.4 API consolidation
#'
#' @description
#' Use `plot()` on a fitted model or an analysis result. The 33 former drawing
#' helpers are no longer exported in 0.2.4. This is a breaking API change:
#' replace old `plot_*()` calls using the table below. Their numerical/drawing
#' implementations remain internal; do not depend on `mfrmr:::` access.
#' The public data accessors [plot_data()] and [plot_data_components()] remain.
#'
#' @section Choosing an input:
#' Use a fitted ordinary model for fit-level views. Use saved G/D studies,
#' information, DIF/DFF, bias, QC or PCA results for their own `plot()` methods.
#' Their public help pages explain the view controls, returned tables and
#' interpretation: [plot.mfrm_information()], [plot.mfrm_dif()],
#' [plot.mfrm_dif_interaction()], [plot.mfrm_anchor_drift()],
#' [plot.mfrm_qc_pipeline()], [plot.mfrm_bias()], [plot.mfrm_residual_pca()]
#' and [plot.mfrm_diagnostics()]. [plot.mfrm_population_intervals()] adds
#' saved latent-regression coefficients and contrasts. Start with these pages
#' for saved results.
#' A common function name does not make every view available for every model.
#' Two-family GPCM and corrected JML retain their separate supported result
#' routes. Plotting does not add intervals, repair unavailable values or equate
#' scales. Use original result objects rather than unclassed lists or CSV tables.
#'
#' @section Single-fit views:
#' `plot(fit)` still draws the ordinary Wright map. New `type` choices include
#' `reliability`, `agreement`, `person_fit`, `rater_severity`, `thresholds`,
#' `residual_qq`, `residual_matrix`, `residual_pca`, `local_dependence`,
#' `marginal_fit`, `marginal_pairwise`, `fair_average`, `displacement`,
#' `facets_chisq`, `bubble`, `guttman`, `unexpected`, `qc`, `apa_summary`,
#' `bias`, `facet_equivalence` and `facet_dashboard`.
#'
#' Use `metric` for the quantity and `view` for a display variant. Reliability
#' accepts `reliability`, `separation` or `strata`. Agreement accepts `exact`
#' or `correlation`, with `view = "pairs"` (default) or `"heatmap"`;
#' `metric = "difference"` applies only to pairs. Severity separation is
#' not agreement. The former `plot_type` argument becomes `view` for marginal
#' fit, fair averages, displacement, facet chi-square, unexpected ratings and
#' facet dashboards. Keep other named controls from the legacy parameter help.
#'
#' Explicit controls are passed only to the chosen renderer; omitted controls
#' retain that renderer's defaults. Inapplicable and conflicting controls are
#' rejected. Supply matching saved `diagnostics` for repeated fit-level views.
#' Some legacy fit-level convenience routes compute diagnostics when omitted;
#' this does not imply model refitting. Saved-result plots reuse their results.
#'
#' @section Replacement table:
#' In the table, `x` is the former first argument, including its original
#' result class. For paired models use `reference` and `comparison`. For an
#' ordered trajectory, `fits` is the existing named list of fitted models.
#'
#' | Former exported helper | Replacement |
#' | --- | --- |
#' | `plot_anchor_drift()` | `plot(x)` |
#' | `plot_apa_figure_one()` | `plot(x, type = "apa_summary")` |
#' | `plot_bias_interaction()` | `plot(x, type = "bias")` |
#' | `plot_bubble()` | `plot(x, type = "bubble")` |
#' | `plot_compare_mfrm()` | `plot(reference, comparison = comparison)` |
#' | `plot_dif_heatmap()` | `plot(x, type = "heatmap")` |
#' | `plot_dif_summary()` | `plot(x, type = "summary")` |
#' | `plot_displacement()` | `plot(x, type = "displacement")` |
#' | `plot_facet_equivalence()` | `plot(x, type = "facet_equivalence")` |
#' | `plot_facet_quality_dashboard()` | `plot(x, type = "facet_dashboard")` |
#' | `plot_facets_chisq()` | `plot(x, type = "facets_chisq")` |
#' | `plot_fair_average()` | `plot(x, type = "fair_average")` |
#' | `plot_guttman_scalogram()` | `plot(x, type = "guttman")` |
#' | `plot_information()` | `plot(x)` |
#' | `plot_interrater_agreement()` | `plot(x, type = "agreement")` |
#' | `plot_local_dependence_heatmap()` | `plot(x, type = "local_dependence")` |
#' | `plot_marginal_fit()` | `plot(x, type = "marginal_fit")` |
#' | `plot_marginal_pairwise()` | `plot(x, type = "marginal_pairwise")` |
#' | `plot_person_fit()` | `plot(x, type = "person_fit")` |
#' | `plot_qc_dashboard()` | `plot(x, type = "qc")` |
#' | `plot_qc_pipeline()` | `plot(x)` |
#' | `plot_rater_agreement_heatmap()` | `plot(x, type = "agreement", view = "heatmap")` |
#' | `plot_rater_severity_profile()` | `plot(x, type = "rater_severity")` |
#' | `plot_rater_trajectory()` | `plot((fits)[[1L]], type = "trajectory", comparison = (fits)[-1L], labels = names(fits))` |
#' | `plot_reliability_snapshot()` | `plot(x, type = "reliability")` |
#' | `plot_residual_matrix()` | `plot(x, type = "residual_matrix")` |
#' | `plot_residual_pca()` | `plot(x, type = "residual_pca")` |
#' | `plot_residual_qq()` | `plot(x, type = "residual_qq")` |
#' | `plot_response_time_review()` | `plot(x)` |
#' | `plot_shrinkage_funnel()` | `plot(x, type = "shrinkage", view = "funnel")` |
#' | `plot_threshold_ladder()` | `plot(x, type = "thresholds")` |
#' | `plot_unexpected()` | `plot(x, type = "unexpected")` |
#' | `plot_wright_unified()` | `plot(x, type = "wright", view = "unified")` |
#'
#' For result objects the object determines the target: for example
#' `plot(info, type = "iif", facet = "Rater")` or
#' `plot(dif, type = "heatmap", metric = "contrast")`.
#' For the fit convenience routes, rename the old bias `plot` argument to
#' `view`; rename the old equivalence `type` to `view`. For residual PCA
#' rename `plot_type` to `view` with `type = "residual_pca"`, or use
#' `plot(pca, type = "loadings")` directly on a saved PCA result.
#'
#' @section Multiple models and occasions:
#' Use `plot(reference, comparison = other, type = "ccc")` or
#' `type = "wright"` for paired displays; `view = "difference"` shows
#' differences without adding difference SEs or model-selection tests.
#' For trajectories use `plot(fits[[1]], type = "trajectory",
#' comparison = fits[-1], labels = names(fits))`. Labels include the first
#' model, and the supplied order is retained. Models must already be placed
#' on an interpretable common scale; plotting performs no linking.
#'
#' @section Returned data and appearance:
#' Use `draw = FALSE` to obtain the existing payload without drawing. Most
#' routes return `mfrm_plot_data`; the QC pipeline retains its verdict-table
#' return, and native interval routes retain their existing ggplot results.
#' Return classes have not been homogenized. Use [plot_data()] for supported
#' payloads and [as_ggplot()] only for supported conversions. Unsupported
#' conversions do not become available through this API change.
#' Titles and cautions stay outside figures by default; keep returned notes
#' when exporting. G/D plots retain their existing color, monochrome, score
#' and facet-panel controls. Save full result objects for later plotting.
#'
#' @seealso [plot.mfrm_fit()], [mfrmr_visual_diagnostics], [mfrmr_output_guide()]
#' @name mfrmr_plot_migration
NULL
