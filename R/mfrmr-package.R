#' mfrmr: Many-Facet Rasch Modeling in R
#'
#' @description
#' `mfrmr` provides estimation, diagnostics, and reporting utilities for
#' many-facet Rasch models (MFRM) using a native R implementation.
#'
#' @useDynLib mfrmr, .registration = TRUE
#'
#' @details
#' If you are new to the package, read the next four steps first and ignore the
#' longer `GPCM`, simulation, and planning notes until the basic route works:
#'
#' 1. Fit with [fit_mfrm()] using `method = "MML"`
#' 2. For `RSM` / `PCM`, run [diagnose_mfrm()] with
#'    `diagnostic_mode = "both"`
#' 3. Read `summary(fit)` and `summary(diag)` before branching
#' 4. Use [plot_qc_dashboard()] and [reporting_checklist()] as the first visual
#'    and reporting screens
#'
#' Recommended workflow:
#'
#' 1. Fit model with [fit_mfrm()]
#' 2. For `RSM` / `PCM`, compute diagnostics with
#'    [diagnose_mfrm()] and prefer `diagnostic_mode = "both"` when you want
#'    legacy residual continuity plus the newer strict marginal-fit screen
#' 3. For `RSM` / `PCM`, run residual PCA with [analyze_residual_pca()] if needed
#' 4. For `RSM` / `PCM`, estimate interactions with [estimate_bias()]
#' 5. For `RSM` / `PCM`, choose a downstream branch:
#'    [reporting_checklist()] for manuscript/report preparation, or
#'    [build_misfit_casebook()] / [build_linking_review()] for operational
#'    misfit or anchor/drift review. After
#'    [build_misfit_casebook()], inspect `casebook$group_view_index` before
#'    moving to source-specific plots.
#' 6. For `RSM` / `PCM`, build narrative/report outputs with
#'    [build_apa_outputs()] and [build_visual_summaries()]
#' 7. Treat `GPCM`, prediction, and planning helpers as advanced scope after
#'    the basic `RSM` / `PCM` route is working cleanly.
#'
#' Guide pages:
#' - [mfrmr_workflow_methods]
#' - [mfrmr_visual_diagnostics]
#' - [mfrmr_reports_and_tables]
#' - [mfrmr_reporting_and_apa]
#' - [mfrmr_linking_and_dff]
#' - [gpcm_capability_matrix]
#' - [mfrmr_compatibility_layer]
#'
#' Companion vignettes:
#' - `vignette("mfrmr-workflow", package = "mfrmr")`
#' - `vignette("mfrmr-mml-and-marginal-fit", package = "mfrmr")`
#' - `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`
#' - `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`
#' - `vignette("mfrmr-linking-and-dff", package = "mfrmr")`
#'
#' @section First 5-minute route:
#' Use this order before exploring the broader feature surface:
#' 1. [fit_mfrm()] with `method = "MML"`
#' 2. [diagnose_mfrm()] with `diagnostic_mode = "both"` for `RSM` / `PCM`
#' 3. `summary(fit)` and `summary(diag)`
#' 4. [plot_qc_dashboard()] for first-pass triage
#' 5. Choose the next branch:
#'    [reporting_checklist()] for reporting,
#'    [build_weighting_audit()] for Rasch-versus-`GPCM` weighting review,
#'    [build_misfit_casebook()] for operational case review, or
#'    [build_linking_review()] for operational linking review
#'
#' @section Advanced scope:
#' After the basic route above:
#' - the package now includes a first-version latent-regression `MML` branch
#'   for ordered-response `RSM` / `PCM` models with a one-dimensional
#'   conditional-normal population model and explicit one-row-per-person
#'   covariates expanded through `stats::model.matrix()`
#' - bounded `GPCM` support is summarized by [gpcm_capability_matrix()]
#' - bounded `GPCM` supports the core fit/summary/scoring/information
#'   path, direct Wright/pathway/CCC plots, residual-PCA follow-up, and the
#'   residual-based diagnostics tables/plots as exploratory tools
#' - posterior-predictive computation, `MCMC` engines, and Docker-based
#'   advanced runtimes are future extensions rather than requirements for the
#'   current bounded `GPCM` route
#' - direct `GPCM` data generation through [build_mfrm_sim_spec()],
#'   [extract_mfrm_sim_spec()], and [simulate_mfrm_data()] is available when
#'   the specification carries both thresholds and slopes
#' - fair-average, APA writer, and broader planning semantics remain
#'   generalized only for `RSM` / `PCM`
#' - `predict_mfrm_population()` remains a scenario-level forecast helper and
#'   should not be described as the latent-regression estimator itself
#' - the current simulation/planning layer still remains role-based for two
#'   non-person facets rather than fully arbitrary-facet, with boundaries
#'   exposed through planner metadata such as `planning_scope`,
#'   `planning_constraints`, and `planning_schema`
#'
#' @section Equal weighting versus bounded GPCM:
#' The package's operational reference route is still the Rasch-family
#' `RSM` / `PCM` branch. That route enforces fixed discrimination and therefore
#' preserves an equal-weighting scoring interpretation across observed ratings.
#'
#' bounded `GPCM` is supported because some users want a slope-aware model-
#' comparison or sensitivity layer inside the same many-facet workflow. However,
#' the package does not treat bounded `GPCM` as a universal replacement for the
#' Rasch-family route. A better fit under `GPCM` should be read as evidence
#' about discrimination-based reweighting, not as an automatic reason to
#' discard the equal-weighting model.
#'
#' Observation weights are a different concept again. Optional `Weight`
#' columns change how observed rating events enter estimation and summaries, but
#' they do not create a free-form facet-weighting scheme and do not alter the
#' fixed-discrimination meaning of `RSM` / `PCM`.
#'
#' Function families:
#' - Model fitting: [fit_mfrm()], [summary.mfrm_fit()], [plot.mfrm_fit()]
#' - Legacy-compatible workflow wrapper: [run_mfrm_facets()], [mfrmRFacets()]
#' - Diagnostics: [diagnose_mfrm()], `summary(diag)`,
#'   [analyze_residual_pca()], [plot_residual_pca()]
#' - Bias and interaction: [estimate_bias()], [estimate_all_bias()],
#'   `summary(bias)`, [bias_interaction_report()], [plot_bias_interaction()]
#' - Differential functioning: [analyze_dff()], [analyze_dif()],
#'   [dif_interaction_table()], [plot_dif_heatmap()], [dif_report()]
#' - Design simulation: [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()],
#'   [simulate_mfrm_data()], [evaluate_mfrm_design()],
#'   [evaluate_mfrm_signal_detection()], [predict_mfrm_population()],
#'   [predict_mfrm_units()], [sample_mfrm_plausible_values()] (including
#'   fit-derived empirical / resampled / skeleton-based simulation
#'   specifications; fixed-calibration unit scoring supports `MML` fits
#'   directly, latent-regression `MML` fits through the fitted population
#'   model when scored units also provide one-row-per-person background data,
#'   and `JML` fits through a post hoc reference-prior EAP layer;
#'   fit-derived simulation specifications also support direct bounded
#'   `GPCM` data generation, while planning/forecasting helpers remain
#'   restricted to `RSM` / `PCM`; curve reports and graph-only exports are
#'   also available for bounded `GPCM`)
#' - Reporting: [build_apa_outputs()], [build_visual_summaries()],
#'   [reporting_checklist()], [apa_table()] for the full `RSM` / `PCM` route;
#'   bounded `GPCM` currently stays on the checklist / direct-table / direct-
#'   plot side instead of the narrative/QC layer
#' - Weighting review: [compare_mfrm()], [build_weighting_audit()],
#'   [compute_information()], [plot_information()]
#' - Case review: [build_misfit_casebook()], [plot_unexpected()],
#'   [plot_displacement()], [plot_marginal_fit()], [plot_marginal_pairwise()]
#' - Linking and scale maintenance: [audit_mfrm_anchors()],
#'   [detect_anchor_drift()], [build_equating_chain()],
#'   [build_linking_review()], [plot_anchor_drift()]
#' - Dashboards: [facet_quality_dashboard()], [plot_facet_quality_dashboard()]
#' - Export / reproducibility: [build_mfrm_manifest()], [build_mfrm_replay_script()],
#'   [build_conquest_overlap_bundle()], [normalize_conquest_overlap_files()],
#'   [normalize_conquest_overlap_tables()],
#'   [audit_conquest_overlap()],
#'   [export_mfrm_bundle()] for the diagnostics-compatible Rasch-family route;
#'   bounded `GPCM` remains outside the current manifest/replay/bundle layer
#' - Equivalence: [analyze_facet_equivalence()], [plot_facet_equivalence()]
#' - Data and anchors: [describe_mfrm_data()], [audit_mfrm_anchors()],
#'   [make_anchor_table()], [load_mfrmr_data()]
#'
#' Data interface:
#' - Input analysis data is long format (one row per observed rating).
#' - Required columns are one person column, one ordered score column, and one
#'   or more non-person facet columns named in `facets = c(...)`.
#' - Score values should be ordered integer categories. Binary `0/1` or `1/2`
#'   input is supported as the two-category Rasch-family special case; by
#'   contrast, fractional score values should be recoded before fitting rather
#'   than relying on automatic coercion.
#' - If `keep_original = FALSE`, unused intermediate categories are collapsed
#'   to a contiguous internal scale and the mapping is stored in
#'   `fit$prep$score_map`.
#' - If the intended scale has unused boundary categories, such as a 1-5 scale
#'   with only 2-5 observed, set `rating_min = 1, rating_max = 5` so the
#'   zero-count boundary category remains in the fitted support. If unused
#'   intermediate categories should also remain in the original scale, set
#'   `keep_original = TRUE`.
#' - `summary(describe_mfrm_data(...))` reports retained zero-count categories
#'   in `Notes`, printed `Caveats`, and `$caveats`; `summary(fit)` carries full
#'   structured rows into printed `Caveats` and `$caveats`, with `Key warnings`
#'   as a short triage subset. Summary-table exports route those rows through
#'   `score_category_caveats` or `analysis_caveats`. Treat adjacent thresholds
#'   as weakly identified when an intermediate category is unobserved.
#' - Optional columns such as `Subset`, `Weight`, and `Group` support linking,
#'   weighted analysis, and fairness-focused follow-up workflows.
#' - Packaged simulation data is available via [load_mfrmr_data()] or `data()`.
#'
#' @section Interpreting output:
#' Core object classes are:
#' - `mfrm_fit`: fitted model parameters and metadata.
#' - `mfrm_diagnostics`: fit, facet-level reliability, and flag diagnostics,
#'   plus inter-rater agreement when one facet is treated as a rater facet.
#' - `mfrm_bias`: interaction bias estimates.
#' - `mfrm_dff` / `mfrm_dif`: differential-functioning contrasts and screening summaries.
#' - `mfrm_population_prediction`: scenario-level forecast summaries for one
#'   future design.
#' - `mfrm_unit_prediction`: posterior summaries for future or partially
#'   observed persons under the fitted scoring basis.
#' - `mfrm_plausible_values`: posterior draws for future or partially observed
#'   persons under the fitted scoring basis.
#' - `mfrm_bundle` families: summary/report bundles and plotting payloads.
#'
#' @section Typical workflow:
#' 1. Prepare long-format data.
#' 2. Fit with [fit_mfrm()].
#' 3. For `RSM` / `PCM`, diagnose with [diagnose_mfrm()] and prefer
#'    `diagnostic_mode = "both"` for final `MML` runs.
#' 4. For `RSM` / `PCM`, run [analyze_dff()] or [estimate_bias()] when
#'    fairness or interaction questions matter.
#' 5. For `RSM` / `PCM`, report with [build_apa_outputs()] and
#'    [build_visual_summaries()].
#' 6. For design planning, move to [build_mfrm_sim_spec()],
#'    [evaluate_mfrm_design()], and [predict_mfrm_population()]. bounded
#'    `GPCM` also supports direct simulation via
#'    [extract_mfrm_sim_spec()] / [simulate_mfrm_data()], but not the broader
#'    planning helpers. Those helpers still assume two non-person facet roles
#'    even though the estimation core supports arbitrary facet counts.
#'    `predict_mfrm_population()` remains the scenario-level forecast helper,
#'    not the latent-regression estimator.
#' 7. For future-unit scoring, retain an `MML` calibration when you want the
#'    fitted marginal model directly, use an active latent-regression `MML`
#'    fit when scored units also provide one-row-per-person background data, or
#'    use a `JML` calibration when a post hoc fixed-calibration EAP layer is
#'    acceptable; then score with
#'    [predict_mfrm_units()] or [sample_mfrm_plausible_values()].
#' 8. For bounded `GPCM`, use [summary.mfrm_fit()],
#'    [diagnose_mfrm()], [analyze_residual_pca()],
#'    [predict_mfrm_units()], [sample_mfrm_plausible_values()],
#'    [compute_information()], [plot_qc_dashboard()], [plot.mfrm_fit()],
#'    [category_structure_report()], [category_curves_report()],
#'    graph-only [facets_output_file_bundle()], direct simulation-spec
#'    generation/data generation, and the residual-based table helpers while
#'    fair-average, APA writer, fit-based export/replay, and planning semantics
#'    are still being generalized. In particular, FACETS-style fair averages
#'    are Rasch-family measure-to-score transformations, so `mfrmr` still keeps
#'    those score-side semantics blocked for bounded `GPCM`. Use
#'    [gpcm_capability_matrix()] as the formal boundary statement.
#'
#' @section Model formulation:
#' The many-facet Rasch model (MFRM; Linacre, 1989) extends the basic Rasch
#' model by incorporating multiple measurement facets into a single linear
#' model on the log-odds scale.
#'
#' **General MFRM equation**
#'
#' For an observation where person \eqn{n} with ability \eqn{\theta_n} is
#' rated by rater \eqn{j} with severity \eqn{\delta_j} on criterion \eqn{i}
#' with difficulty \eqn{\beta_i}, the probability of observing category
#' \eqn{k} (out of \eqn{K} ordered categories) is:
#'
#' \deqn{P(X_{nij} = k \mid \theta_n, \delta_j, \beta_i, \tau) =
#'   \frac{\exp\bigl[\sum_{s=1}^{k}(\theta_n - \delta_j - \beta_i - \tau_s)\bigr]}
#'        {\sum_{c=0}^{K}\exp\bigl[\sum_{s=1}^{c}(\theta_n - \delta_j - \beta_i - \tau_s)\bigr]}}
#'
#' where \eqn{\tau_s} are the Rasch-Andrich threshold (step) parameters and
#' \eqn{\sum_{s=1}^{0}(\cdot) \equiv 0} by convention.  Additional facets
#' enter as additive terms in the linear predictor
#' \eqn{\eta = \theta_n - \delta_j - \beta_i - \ldots}.
#'
#' This formulation generalises to any number of facets; the
#' `facets` argument to [fit_mfrm()] accepts an arbitrary-length
#' character vector.
#'
#' **Rating Scale Model (RSM)**
#'
#' Under the RSM (Andrich, 1978), all levels of the step facet share a
#' single set of threshold parameters \eqn{\tau_1, \ldots, \tau_K}.
#'
#' **Partial Credit Model (PCM)**
#'
#' Under the PCM (Masters, 1982), each level of the designated `step_facet`
#' has its own threshold vector on the package's common observed score scale.
#' In the current implementation, threshold locations may vary by step-facet
#' level, but the fitted score range is still defined by one global category
#' set taken from the observed data.
#'
#' **Ordered-response scope**
#'
#' The current public response-model scope is ordered categorical only.
#' Binary responses are the \eqn{K = 1} special case of the same formulation,
#' so they are handled through the ordinary ordered-score interface. This means
#' `mfrmr` supports ordered binary and ordered polytomous data under `RSM` and
#' `PCM`, plus a narrow bounded `GPCM` branch with one designated
#' `slope_facet` that currently must equal `step_facet`. Unordered
#' nominal/multinomial response models are not yet implemented.
#'
#' @section Estimation methods:
#' **Marginal Maximum Likelihood (MML)**
#'
#' MML integrates over the person ability distribution using Gauss-Hermite
#' quadrature (Bock & Aitkin, 1981):
#'
#' \deqn{L = \prod_{n} \int P(\mathbf{X}_n \mid \theta, \boldsymbol{\delta})
#'   \, \phi(\theta) \, d\theta
#'   \approx \prod_{n} \sum_{q=1}^{Q} w_q \,
#'   P(\mathbf{X}_n \mid \theta_q, \boldsymbol{\delta})}
#'
#' where \eqn{\phi(\theta)} is the assumed normal prior and
#' \eqn{(\theta_q, w_q)} are quadrature nodes and weights.  Person
#' estimates are obtained post-hoc via Expected A Posteriori (EAP):
#'
#' \deqn{\hat{\theta}_n^{\mathrm{EAP}} =
#'   \frac{\sum_q \theta_q \, w_q \, L(\mathbf{X}_n \mid \theta_q)}
#'        {\sum_q w_q \, L(\mathbf{X}_n \mid \theta_q)}}
#'
#' MML avoids the incidental-parameter problem and is generally preferred
#' for smaller samples.
#'
#' **Joint Maximum Likelihood (JML)**
#'
#' JML estimates all person and facet parameters simultaneously as fixed
#' effects by maximising the joint log-likelihood
#' \eqn{\ell(\boldsymbol{\theta}, \boldsymbol{\delta} \mid \mathbf{X})}
#' directly.  It does not assume a parametric person distribution, which
#' can be advantageous when the population shape is strongly non-normal,
#' but parameter estimates are known to be biased when the number of
#' persons is small relative to the number of items (Neyman & Scott, 1948).
#' The package still accepts `"JMLE"` as a backward-compatible alias, but
#' user-facing summaries and documentation use `"JML"` as the public label.
#'
#' See [fit_mfrm()] for practical guidance on choosing between the two.
#'
#' @section Strict marginal diagnostics and literature positioning:
#' For `RSM` / `PCM`, `diagnose_mfrm(..., diagnostic_mode = "both")`
#' separates two targets:
#'
#' - `legacy`: residual/EAP-oriented diagnostics used for continuity with the
#'   earlier package surface
#' - `marginal_fit`: latent-integrated expected counts and pairwise summaries
#'   derived from the fitted `MML` posterior bundle
#'
#' Write the posterior weight for person \eqn{n} at quadrature node
#' \eqn{q} as
#'
#' \deqn{\omega_{nq} =
#'   \frac{w_q \, P(\mathbf{X}_n \mid \theta_q, \hat{\boldsymbol{\delta}})}
#'        {\sum_{r=1}^{Q} w_r \, P(\mathbf{X}_n \mid \theta_r, \hat{\boldsymbol{\delta}})}}
#'
#' and let \eqn{g} denote a grouped cell, facet combination, or pairwise
#' comparison target. Then the package's strict first-order expected counts are
#' of the form
#'
#' \deqn{E_{\hat{\delta}}(N_{gc}) =
#'   \sum_{n=1}^{N}\sum_{q=1}^{Q}
#'   \omega_{nq} \, I(n \in g)\,
#'   P(X_n = c \mid \theta_q, \hat{\boldsymbol{\delta}}).}
#'
#' Pairwise local-dependence screens use the same posterior bundle but replace
#' the one-category event \eqn{X_n = c} with agreement or adjacency events for
#' the relevant pair of facet levels.
#'
#' This places the current package closest to limited-information item-fit and
#' generalized-residual traditions rather than to a single definitive omnibus
#' test. In the current release, these ideas are adapted to a many-facet
#' screening layer rather than implemented as literal `S-X2` or formal
#' generalized-residual tests. Orlando and Thissen (2000, 2003) motivate the
#' limited-information item-fit family, Haberman and Sinharay (2013) motivate
#' generalized residual reasoning, Sinharay et al. (2006) motivate posterior
#' predictive follow-up as a separate checking family, and Sinharay and Monroe
#' (2025) argue that practitioners should match fit procedures to intended
#' uses, examine practical significance, and avoid relying on any one statistic
#' in isolation. `mfrmr` therefore reports strict marginal diagnostics as
#' structured screening evidence, not as a completed universal accept/reject
#' test battery.
#'
#' In many-facet practice, this strict screening layer complements rather than
#' replaces the usual MFRM tools for fit, severity/leniency review, and
#' agreement. Facet-level separation/reliability summarizes how distinctly a
#' facet is measured, whereas inter-rater agreement summarizes observed
#' agreement across matched contexts; they should not be treated as
#' interchangeable quantities.
#'
#' @section Statistical background:
#' Key statistics reported throughout the package:
#'
#' **Infit (Information-Weighted Mean Square)**
#'
#' Weighted average of squared standardized residuals, where weights are the
#' model-based variance of each observation:
#'
#' \deqn{\mathrm{Infit}_j = \frac{\sum_i Z_{ij}^2 \, \mathrm{Var}_i \, w_i}
#'                               {\sum_i \mathrm{Var}_i \, w_i}}
#'
#' Expected value is 1.0 under model fit.  Values below 0.5 suggest overfit
#' (Mead-style responses); values above 1.5 suggest underfit (noise or
#' misfit).  Infit is most sensitive to unexpected patterns among on-target
#' observations (Wright & Masters, 1982).
#'
#' Note: The 0.5--1.5 range is a widely used rule of thumb (Bond & Fox,
#' 2015).  Acceptable ranges may differ by context: 0.6--1.4 for high-stakes
#' testing; 0.7--1.3 for clinical instruments; up to 0.5--1.7 for surveys
#' and exploratory work (Linacre, 2002).
#'
#' **Outfit (Unweighted Mean Square)**
#'
#' Simple average of squared standardized residuals:
#'
#' \deqn{\mathrm{Outfit}_j = \frac{\sum_i Z_{ij}^2 \, w_i}{\sum_i w_i}}
#'
#' Same expected value and flagging thresholds as Infit, but more sensitive
#' to extreme off-target outliers (e.g., a high-ability person scoring the
#' lowest category).
#'
#' **ZSTD (Standardized Fit Statistic)**
#'
#' Wilson-Hilferty cube-root transformation that converts the mean-square
#' chi-square ratio to an approximate standard normal deviate:
#'
#' \deqn{\mathrm{ZSTD} = \frac{\mathrm{MnSq}^{1/3} - (1 - 2/(9\,\mathit{df}))}
#'                            {\sqrt{2/(9\,\mathit{df})}}}
#'
#' Values near 0 indicate expected fit; \eqn{|\mathrm{ZSTD}| > 2} flags
#' potential misfit at the 5\% level, and \eqn{|\mathrm{ZSTD}| > 3} at the
#' 1\% level (Wright & Linacre, 1994).  ZSTD is reported alongside every
#' Infit and Outfit value.
#'
#' **PTMEA (Point-Measure Correlation)**
#'
#' Pearson correlation between observed scores and estimated person measures
#' within each facet level.  Positive values indicate that scoring aligns
#' with the latent trait dimension; negative values suggest reversed
#' orientation or scoring errors.
#'
#' **Separation**
#'
#' Package-reported separation is the ratio of adjusted true standard deviation
#' to root-mean-square measurement error:
#'
#' \deqn{G = \frac{\mathrm{SD}_{\mathrm{adj}}}{\mathrm{RMSE}}}
#'
#' where \eqn{\mathrm{SD}_{\mathrm{adj}} =
#' \sqrt{\mathrm{ObservedVariance} - \mathrm{ErrorVariance}}}. Higher values
#' indicate the facet discriminates more statistically distinct levels along the
#' measured variable. In `mfrmr`, `Separation` is the model-based value and
#' `RealSeparation` provides a more conservative companion based on `RealSE`.
#'
#' **Reliability**
#'
#' \deqn{R = \frac{G^2}{1 + G^2}}
#'
#' Analogous to Cronbach's alpha or KR-20 for the reproducibility of element
#' ordering. In `mfrmr`, `Reliability` is the model-based value and
#' `RealReliability` gives the conservative companion based on `RealSE`. For
#' `MML`, these are anchored to observed-information `ModelSE`
#' estimates for non-person facets; `JML` keeps them as exploratory summaries.
#'
#' **Strata**
#'
#' Number of statistically distinguishable groups of elements:
#'
#' \deqn{H = \frac{4G + 1}{3}}
#'
#' Three or more strata are commonly used as a practical target
#' (Wright & Masters, 1982), but in this package the estimate inherits the
#' same approximation limits as the separation index.
#'
#' @section Key references:
#' - Andrich, D. (1978). A rating formulation for ordered response
#'   categories. *Psychometrika*, 43, 561--573.
#' - Bond, T. G., & Fox, C. M. (2015). *Applying the Rasch model* (3rd
#'   ed.). Routledge.
#' - Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm. *Psychometrika*, 46,
#'   443--459.
#' - Haberman, S. J., & Sinharay, S. (2013). Generalized residuals for general
#'   models for contingency tables with application to item response theory.
#'   *Journal of the American Statistical Association*, 108, 1435--1444.
#' - Eckes, T. (2005). Examining rater effects in TestDaF writing and speaking
#'   performance assessments: A many-facet Rasch analysis. *Language Assessment
#'   Quarterly*, 2, 197--221.
#' - Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#' - Linacre, J. M. (2002). What do Infit and Outfit, mean-square and
#'   standardized mean? *Rasch Measurement Transactions*, 16(2), 878.
#' - Masters, G. N. (1982). A Rasch model for partial credit scoring.
#'   *Psychometrika*, 47, 149--174.
#' - Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for
#'   dichotomous item response theory models. *Applied Psychological
#'   Measurement*, 24, 50--64.
#' - Orlando, M., & Thissen, D. (2003). Further investigation of the
#'   performance of S-X2: An item fit index for use with dichotomous item
#'   response theory models. *Applied Psychological Measurement*, 27, 289--298.
#' - Sinharay, S., Johnson, M. S., & Stern, H. S. (2006). Posterior predictive
#'   assessment of item response theory models. *Applied Psychological
#'   Measurement*, 30, 298--321.
#' - Sinharay, S., & Monroe, S. (2025). Assessment of fit of item response
#'   theory models: A critical review of the status quo and some future
#'   directions. *British Journal of Mathematical and Statistical Psychology*,
#'   78, 711--733.
#' - Wright, B. D., & Masters, G. N. (1982). *Rating scale analysis*.
#'   MESA Press.
#' - Wright, B. D., & Linacre, J. M. (1994). Reasonable mean-square fit
#'   values. *Rasch Measurement Transactions*, 8(3), 370.
#'
#' @section Model selection:
#' **RSM vs PCM**
#'
#' The Rating Scale Model (RSM; Andrich, 1978) assumes all levels of the
#' step facet share identical threshold parameters.  The Partial Credit
#' Model (PCM; Masters, 1982) allows each level of the `step_facet` to have
#' its own set of thresholds on the package's shared observed score scale.
#' Use RSM when the rating rubric is identical across all items/criteria;
#' use PCM when category boundaries are expected to vary by item or criterion.
#' In the current implementation, PCM still assumes one common observed score
#' support across the fitted data, so it should not be described as a fully
#' mixed-category model with arbitrary item-specific category counts.
#'
#' **MML vs JML**
#'
#' Marginal Maximum Likelihood (MML) integrates over the person ability
#' distribution using Gauss-Hermite quadrature and does not directly estimate
#' person parameters; person estimates are computed post-hoc via Expected A
#' Posteriori (EAP).  Joint Maximum Likelihood (JML) estimates all person
#' and facet parameters simultaneously as fixed effects; `"JMLE"` remains a
#' backward-compatible alias.
#'
#' MML is generally preferred for smaller samples because it avoids the
#' incidental-parameter problem of JML.  JML does not assume a normal person
#' distribution and can be lighter computationally in some settings, which
#' may be an advantage when the population shape is strongly non-normal.
#'
#' See [fit_mfrm()] for usage.
#'
#' **Fixed-calibration scoring after fitting**
#'
#' [predict_mfrm_units()] and [sample_mfrm_plausible_values()] score future or
#' partially observed persons on a quadrature grid under the fitted scoring
#' basis. For ordinary `MML` fits, these summaries inherit the fitted marginal
#' calibration directly. For latent-regression `MML` fits, they use the fitted
#' one-dimensional conditional normal population model and therefore require
#' one-row-per-person background data for the scored units when the fitted
#' population model includes covariates. Intercept-only latent-regression fits
#' (`population_formula = ~ 1`) can reconstruct that minimal person table from
#' the scored person IDs. For `JML` fits, `mfrmr` uses the fitted facet and
#' step parameters together with a standard normal reference prior introduced
#' only for the post hoc scoring layer. This is useful for practical
#' fixed-scale scoring, but it should still be described as a limited
#' approximation rather than as full ConQuest-style population modeling.
#'
#' **Current ConQuest overlap**
#'
#' The package now includes a first-version latent-regression `MML` branch, but
#' the overlap with ConQuest should still be described conservatively. The
#' defensible shared ground is:
#' ordered-response `RSM` / `PCM`, one latent dimension, a conditional-normal
#' person population model, and person covariates supplied through an explicit
#' one-row-per-person table and expanded through the package-built model
#' matrix. Categorical person covariates carry fitted levels and contrasts into
#' scoring. This is a scoped overlap, not a claim of broad ConQuest numerical
#' equivalence for arbitrary imported design matrices, multidimensional models,
#' imported design specifications, or the full plausible-values workflow.
#'
#' @examples
#' mfrm_threshold_profiles()
#' list_mfrmr_data()
#'
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   model = "RSM",
#'   quad_points = 7
#' )
#' diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
#' summary(diag)
#' }
#'
#' @importFrom dplyr across all_of any_of arrange bind_cols bind_rows
#'   case_when coalesce count desc distinct everything filter group_by
#'   if_all inner_join lag last left_join mutate n n_distinct pull rename
#'   row_number rowwise select slice slice_head slice_tail summarise
#'   summarize transmute ungroup
#' @importFrom tidyr drop_na expand_grid pivot_wider replace_na unite
#' @importFrom tibble as_tibble column_to_rownames tibble
#' @importFrom purrr map map_dbl map_dfr pmap_chr
#' @importFrom stringr regex str_count str_detect str_extract str_pad
#'   str_replace str_replace_all str_trunc
#' @importFrom psych cor.smooth describe principal
#' @importFrom utils combn data globalVariables head modifyList packageVersion
#'   read.csv tail write.csv
#' @importFrom lifecycle deprecate_soft
#' @importFrom rlang .data
#' @importFrom stats na.omit optim p.adjust pchisq pt rnorm sd setNames uniroot
#'
#' @name mfrmr-package
"_PACKAGE"
