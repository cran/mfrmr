#' mfrmr Visual Diagnostics Map
#'
#' @description
#' Quick guide to choosing the right base-R diagnostic plot in `mfrmr`.
#' Use this page when you know the analysis question but do not yet know
#' which plotting helper or `plot()` method to call.
#'
#' If you are preparing figures for a report, start with
#' [reporting_checklist()] and inspect the `"Visual Displays"` rows first.
#' Those rows now map directly onto the public plotting family covered on this
#' page, so the checklist can act as a plot-readiness router rather than just a
#' manuscript checklist.
#'
#' This guide is primarily for diagnostics-based `RSM` / `PCM` workflows.
#' First-release `GPCM` fits now support the residual-based diagnostics stack
#' through [diagnose_mfrm()], [plot_unexpected()], [plot_displacement()],
#' [plot_interrater_agreement()], [plot_facets_chisq()],
#' [plot_residual_pca()], and [plot_qc_dashboard()] with an explicit
#' fair-average placeholder, in addition to the core summary,
#' posterior-scoring, design-weighted-information path via
#' [compute_information()] / [plot_information()], and Wright/pathway/CCC fit
#' plots. For `GPCM`, treat residual-based mean-square screens as exploratory
#' rather than strict Rasch-style invariance tests because the discrimination
#' parameter is free. FACETS-style fair averages are Rasch-family
#' measure-to-score transformations, so fair-average visuals themselves and
#' broader compatibility exports are still outside the validated `GPCM`
#' boundary. Use [gpcm_capability_matrix()] when you need the formal helper
#' boundary before choosing a `GPCM` follow-up plot route.
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
#' - "Which facet/category cells drive strict marginal misfit?"
#'   Use [plot_marginal_fit()].
#' - "Which level pairs drive strict local-dependence follow-up?"
#'   Use [plot_marginal_pairwise()].
#' - "Do raters agree and do facets separate meaningfully?"
#'   Use [plot_interrater_agreement()] and [plot_facets_chisq()].
#' - "Is there notable residual structure after the main Rasch dimension?"
#'   Use [plot_residual_pca()].
#' - "Which interaction cells or facet levels drive bias screening results?"
#'   Use [plot_bias_interaction()].
#' - "I need one compact triage screen first."
#'   Use [plot_qc_dashboard()] for `RSM` / `PCM`. First-release `GPCM`
#'   can also use [plot_qc_dashboard()], but its fair-average panel remains an
#'   explicit unavailable placeholder because that panel's score-metric
#'   semantics have not yet been generalized beyond the Rasch-family branch.
#' - "Which figures are already supported by my current run?"
#'   Use [reporting_checklist()] and review the `"Visual Displays"` rows before
#'   choosing the next plot.
#' - "Where should this figure go in a paper or appendix?"
#'   Use [visual_reporting_template()] for a static reporting-use table, then
#'   cross-check run-specific availability with `reporting_checklist()$visual_scope`.
#' - "Do I need a 3D-style category probability surface?"
#'   Use `plot(fit, type = "ccc_surface", draw = FALSE)` to get a
#'   theta-by-category-by-probability payload for exploratory teaching or
#'   downstream interactive rendering. Keep 2D pathway/CCC plots as the
#'   default reporting figures.
#'
#' @section Recommended visual route:
#' 1. If you are drafting a report, run [reporting_checklist()] first and read
#'    the `"Visual Displays"` rows as the plot-readiness layer.
#' 2. Start with [plot_qc_dashboard()] for one-page triage.
#' 3. Move to [plot_unexpected()], [plot_displacement()],
#'    [plot_marginal_fit()], [plot_marginal_pairwise()], and
#'    [plot_interrater_agreement()] for flagged local issues.
#' 4. Use `plot(fit, type = "wright")`, `plot(fit, type = "pathway")`,
#'    and `plot_residual_pca()` for structural interpretation.
#' 5. Use [plot_bias_interaction()], [plot_anchor_drift()], and
#'    [plot_information()] when the checklist or dashboard points to
#'    interaction, linking, or precision follow-up.
#' 6. Use `plot(..., draw = FALSE)` when you want reusable plotting payloads
#'    instead of immediate graphics.
#' 7. Use `plot(fit, type = "ccc_surface", draw = FALSE)` only when you need
#'    a 3D-ready category-probability payload; `mfrmr` intentionally does not
#'    add a package-native plotly/rgl renderer for this route.
#' 8. Use `preset = "publication"` when you want the package's cleaner
#'    manuscript-oriented styling.
#'
#' @section Visual coverage for this release:
#' This release treats the plotting layer as sufficient when the current run
#' supports all of the following follow-up roles through public helpers:
#' - First-pass triage:
#'   [plot_qc_dashboard()] or the `"Visual Displays"` rows from
#'   [reporting_checklist()].
#' - Structural interpretation:
#'   `plot(fit, type = "wright")`, `plot(fit, type = "pathway")`,
#'   `plot(fit, type = "ccc")`, and [plot_residual_pca()].
#' - Local issue follow-up:
#'   [plot_unexpected()], [plot_displacement()],
#'   [plot_interrater_agreement()], and [plot_bias_interaction()].
#' - Strict marginal follow-up:
#'   [plot_marginal_fit()] and [plot_marginal_pairwise()] for
#'   `diagnostic_mode = "both"`.
#' - Reporting/export handoff:
#'   [build_visual_summaries()] and `draw = FALSE` routes that return reusable
#'   `mfrm_plot_data` payloads for downstream review and export. When step
#'   estimates are available, `build_visual_summaries()` also exposes
#'   `$plot_payloads$category_probability_surface`.
#' - 3D-ready exploratory handoff:
#'   `plot(fit, type = "ccc_surface", draw = FALSE)` returns a
#'   theta-by-category-by-probability `mfrm_plot_data` payload. This is not a
#'   default APA/reporting figure and does not load plotly/rgl.
#'
#' @section 3D and surface payloads:
#' The package currently treats 3D as an exploratory data handoff, not as a
#' default plotting layer. The supported route is
#' `plot(fit, type = "ccc_surface", draw = FALSE)`, which returns
#' `surface`, `categories`, `category_support`, `groups`, `axis_contract`,
#' `renderer_contract`, `interpretation_guide`, and `reporting_policy` tables
#' inside an `mfrm_plot_data` object. These columns can be passed to an
#' external renderer if needed, while `category_support` and
#' `interpretation_guide` should be checked before interpreting retained
#' zero-frequency categories or adjacent threshold ridges.
#'
#' Do not replace the standard 2D Wright map, pathway map, CCC plot,
#' heatmap/profile diagnostics, or information curves with 3D figures in
#' routine reports. In particular, 3D Wright maps are discouraged because
#' perspective and occlusion obscure the shared-scale comparison that the
#' Wright map is meant to support.
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
#'   \item{`plot_marginal_fit()`}{Posterior-integrated first-order category
#'   residuals. Best for seeing which facet/category cells drive strict
#'   marginal flags.}
#'   \item{`plot_marginal_pairwise()`}{Posterior-integrated exact/adjacent
#'   agreement residuals. Best for exploratory local-dependence follow-up after
#'   strict marginal flags.}
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
#' - 3D-ready category surface: use as an exploratory view of the same
#'   category-probability information, not as a replacement for the 2D
#'   pathway/CCC figures in reports. Read `category_support` first when a
#'   retained category has zero observed responses.
#' - Unexpected / displacement: use as screening tools, not final evidence
#'   by themselves.
#' - Strict marginal and pairwise local-dependence plots are exploratory
#'   follow-up layers for `diagnostic_mode = "both"`, not standalone
#'   inferential tests.
#' - Inter-rater agreement and facet variability address different questions:
#'   agreement concerns scoring consistency, whereas variability concerns
#'   whether facet elements are statistically distinguishable.
#' - Residual PCA and bias plots should be interpreted as follow-up layers
#'   after the main fit screen, not as first-pass diagnostics.
#'
#' @section Typical workflow:
#' - Figure-readiness route:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()] ->
#'   inspect `"Visual Displays"` rows -> chosen public plot helper.
#' - Quick screening:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [plot_qc_dashboard()].
#' - Strict marginal follow-up:
#'   [diagnose_mfrm()] with `diagnostic_mode = "both"` ->
#'   [plot_marginal_fit()] ->
#'   [plot_marginal_pairwise()].
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
#'   [gpcm_capability_matrix], [visual_reporting_template()],
#'   [plot.mfrm_fit()], [plot_qc_dashboard()],
#'   [plot_unexpected()], [plot_displacement()], [plot_marginal_fit()],
#'   [plot_marginal_pairwise()], [plot_interrater_agreement()],
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
#'   method = "MML",
#'   maxit = 200
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#' checklist <- reporting_checklist(fit, diagnostics = diag)
#' visual_reporting_template("manuscript")
#' subset(
#'   checklist$checklist,
#'   Section == "Visual Displays" & Item %in% c("QC / facet dashboard", "Strict marginal visuals"),
#'   c("Item", "Available", "NextAction")
#' )
#'
#' qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE, preset = "publication")
#' qc$data$plot
#'
#' p_marg <- plot_marginal_fit(diag, draw = FALSE, preset = "publication")
#' p_marg$data$preset
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

#' Figure-reporting template for visual diagnostics
#'
#' @description
#' Return a compact, beginner-oriented template that explains where each
#' visual family normally belongs in a report, which helper to call, what to
#' say, and what not to claim. Use this static table together with the dynamic
#' `reporting_checklist(fit, diagnostics)$visual_scope` table: the template
#' answers "how should I use this figure?", while the checklist answers "is
#' this figure ready for the current run?".
#'
#' @param scope Which part of the template to return: `"all"` (default),
#'   `"manuscript"`, `"appendix"`, `"diagnostic"`, or `"surface"`.
#'
#' @return A data.frame with columns:
#' - `FigureFamily`: short visual family label.
#' - `Scope`: broad reporting role used for filtering.
#' - `PrimaryHelper`: public helper or plot route.
#' - `DefaultPlacement`: recommended location in a report.
#' - `WhatToReport`: wording focus for results sections or captions.
#' - `CaptionSkeleton`: caption starter that must be tailored to the study.
#' - `ResultsWording`: results-sentence starter that must be checked against
#'   the fitted object and diagnostics.
#' - `WhatNotToClaim`: common overclaim to avoid.
#' - `BeginnerCheck`: first thing a new user should inspect.
#' - `ThreeDPolicy`: whether 3D is recommended, discouraged, or payload-only.
#'
#' @details
#' This helper is intentionally conservative. It does not inspect a fitted
#' object and does not certify that a plot is available. Run
#' [reporting_checklist()] for run-specific readiness, then use this table to
#' decide how to describe the resulting figure.
#'
#' @examples
#' visual_reporting_template()
#' visual_reporting_template("manuscript")
#' visual_reporting_template("surface")
#' @export
visual_reporting_template <- function(scope = c("all", "manuscript", "appendix", "diagnostic", "surface")) {
  scope <- match.arg(scope)

  out <- data.frame(
    FigureFamily = c(
      "Wright map",
      "Pathway map",
      "Category characteristic curves",
      "Category probability surface",
      "Information curves",
      "QC dashboard",
      "Unexpected / displacement",
      "Strict marginal visuals",
      "Bias / DIF visuals",
      "Residual PCA"
    ),
    Scope = c(
      "manuscript",
      "manuscript",
      "manuscript",
      "surface",
      "manuscript",
      "diagnostic",
      "appendix",
      "appendix",
      "diagnostic",
      "appendix"
    ),
    PrimaryHelper = c(
      "plot(fit, type = \"wright\", preset = \"publication\")",
      "plot(fit, type = \"pathway\", preset = \"publication\")",
      "plot(fit, type = \"ccc\", preset = \"publication\")",
      "plot(fit, type = \"ccc_surface\", draw = FALSE)",
      "compute_information(fit) -> plot_information(..., preset = \"publication\")",
      "plot_qc_dashboard(fit, diagnostics = diagnostics, preset = \"publication\")",
      "plot_unexpected(); plot_displacement()",
      "plot_marginal_fit(); plot_marginal_pairwise()",
      "plot_bias_interaction(); plot_dif_heatmap()",
      "analyze_residual_pca() -> plot_residual_pca()"
    ),
    DefaultPlacement = c(
      "Main text when targeting, spread, or shared-logit interpretation is central.",
      "Main text or category-functioning subsection for ordered-category interpretation.",
      "Main text or appendix; pair with pathway when category behavior is central.",
      "Appendix, teaching, audit, or downstream interactive rendering only.",
      "Main text when precision or targeting across theta is a substantive claim.",
      "Screening dashboard; usually methods appendix or internal triage rather than the final main figure.",
      "Case-review appendix or quality-control supplement.",
      "Diagnostic appendix after diagnostic_mode = \"both\".",
      "Main text only if interaction/DIF is a study question; otherwise diagnostic appendix.",
      "Diagnostic appendix or sensitivity discussion."
    ),
    WhatToReport = c(
      "Describe whether persons, facet levels, and thresholds overlap on the same logit scale.",
      "Describe expected-score progression and the theta regions where categories dominate.",
      "Describe whether categories peak in the intended order and whether adjacent curves separate.",
      "Describe it as exploratory category-probability support, not as a default manuscript figure.",
      "Describe where measurement information is highest or weakest across theta.",
      "Describe which components triggered follow-up, not a single pass/fail publication verdict.",
      "Describe which responses or levels need local review.",
      "Describe which facet/category cells or pairwise structures need follow-up.",
      "Describe screened interaction patterns with low-count and threshold caveats.",
      "Describe residual structure as exploratory follow-up after the main fit screen."
    ),
    CaptionSkeleton = c(
      "Figure X. Wright map showing person measures, facet-level locations, and step thresholds on the shared logit scale.",
      "Figure X. Expected score pathway across theta, with dominant-category regions for the fitted rating scale.",
      "Figure X. Category characteristic curves showing fitted category probabilities across theta.",
      "Appendix Figure X. Exploratory category-probability surface showing theta, retained category index, and fitted probability.",
      "Figure X. Test information curve showing where the fitted model provides relatively stronger or weaker measurement precision.",
      "Appendix Figure X. Quality-control dashboard summarizing diagnostic components that require follow-up.",
      "Appendix Figure X. Local response or level-review display for unexpected responses and displacement diagnostics.",
      "Appendix Figure X. Strict marginal diagnostic display for retained facet/category or pairwise follow-up evidence.",
      "Figure/Appendix Figure X. Bias or differential-functioning screening display for the specified facet pair or group contrast.",
      "Appendix Figure X. Residual PCA scree or loading display used for exploratory residual-structure review."
    ),
    ResultsWording = c(
      "The Wright map was inspected to evaluate targeting and shared-scale overlap among persons, facet levels, and thresholds.",
      "The pathway plot was inspected to evaluate whether expected scores and dominant-category regions progressed in the intended order.",
      "The category characteristic curves were inspected to evaluate the ordering and separation of fitted response categories.",
      "The category-probability surface was used as exploratory support for understanding the fitted category-probability structure.",
      "The information curve was inspected to identify theta regions with relatively stronger or weaker measurement precision.",
      "The QC dashboard was used as a triage screen to identify components requiring more specific diagnostic follow-up.",
      "Unexpected-response and displacement displays were used to identify local cases or levels requiring review.",
      "Strict marginal displays were used as follow-up evidence for facet/category and pairwise local-dependence patterns.",
      "Bias/DIF displays were used to screen interaction or group-functioning patterns under the documented screening threshold.",
      "Residual PCA displays were used as exploratory follow-up for residual structure after the main model dimension."
    ),
    WhatNotToClaim = c(
      "Do not present targeting as proof of global model fit.",
      "Do not treat smooth category progression as proof that the rating scale is valid.",
      "Do not overstate overlapping curves as definitive category failure without category counts and context.",
      "Use the surface as exploratory mfrmr output or downstream renderer input; prefer 2D CCC/pathway plots for reports.",
      "Do not ignore the precision tier or approximation caveats used to compute the curve.",
      "Do not cite the dashboard alone as inferential evidence.",
      "Do not interpret a single flagged case as final evidence by itself.",
      "Do not treat strict marginal visuals as standalone hypothesis tests.",
      "Do not claim formal DIF unless the design and inferential route support that wording.",
      "Do not treat residual PCA as a standalone dimensionality test."
    ),
    BeginnerCheck = c(
      "Check gaps between person density and thresholds/facet levels.",
      "Check whether the dominant-category bands progress in the expected order.",
      "Check whether every retained category has a visible peak or clear role.",
      "Read surface$data$category_support and surface$data$interpretation_guide before rendering.",
      "Check whether the information peak covers the theta region of interest.",
      "Open the component rows or plots behind any dashboard warning.",
      "Sort by magnitude and inspect repeated patterns, not isolated extremes only.",
      "Confirm diagnostic_mode = \"both\" and inspect low-count or sparse-cell caveats.",
      "Confirm the tested facet pair, low-count cells, and screening threshold.",
      "Start with the scree plot, then inspect loadings only for targeted follow-up."
    ),
    ThreeDPolicy = c(
      "2D recommended; 3D Wright maps are discouraged.",
      "2D report default.",
      "2D report default.",
      "advanced surface data only; no package-native interactive renderer.",
      "2D curve route active; 3D information surface is deferred.",
      "2D dashboard only; 3D not recommended.",
      "2D point/profile views preferred.",
      "2D heatmap/bar views preferred.",
      "2D heatmap/profile views preferred.",
      "2D scree/loadings preferred."
    ),
    stringsAsFactors = FALSE
  )

  if (identical(scope, "all")) {
    return(out)
  }
  out[out$Scope == scope, , drop = FALSE]
}
