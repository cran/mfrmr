#' mfrmr: Many-Facet Rasch Modeling in R
#'
#' @description
#' `mfrmr` provides estimation, diagnostics, and reporting utilities for
#' many-facet Rasch models (MFRM) using a native R implementation.
#'
#' @details
#' Recommended workflow:
#'
#' 1. Fit model with [fit_mfrm()]
#' 2. Compute diagnostics with [diagnose_mfrm()]
#' 3. Run residual PCA with [analyze_residual_pca()] if needed
#' 4. Estimate interactions with [estimate_bias()]
#' 5. Build narrative/report outputs with [build_apa_outputs()] and [build_visual_summaries()]
#'
#' Guide pages:
#' - [mfrmr_workflow_methods]
#' - [mfrmr_visual_diagnostics]
#' - [mfrmr_reports_and_tables]
#' - [mfrmr_reporting_and_apa]
#' - [mfrmr_linking_and_dff]
#' - [mfrmr_compatibility_layer]
#'
#' Companion vignettes:
#' - `vignette("mfrmr-workflow", package = "mfrmr")`
#' - `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`
#' - `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`
#' - `vignette("mfrmr-linking-and-dff", package = "mfrmr")`
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
#'   specifications; fixed-calibration unit scoring currently requires
#'   `method = "MML"`)
#' - Reporting: [build_apa_outputs()], [build_visual_summaries()],
#'   [reporting_checklist()], [apa_table()]
#' - Dashboards: [facet_quality_dashboard()], [plot_facet_quality_dashboard()]
#' - Export / reproducibility: [build_mfrm_manifest()], [build_mfrm_replay_script()],
#'   [export_mfrm_bundle()]
#' - Equivalence: [analyze_facet_equivalence()], [plot_facet_equivalence()]
#' - Data and anchors: [describe_mfrm_data()], [audit_mfrm_anchors()],
#'   [make_anchor_table()], [load_mfrmr_data()]
#'
#' Data interface:
#' - Input analysis data is long format (one row per observed rating).
#' - Packaged simulation data is available via [load_mfrmr_data()] or `data()`.
#'
#' @section Interpreting output:
#' Core object classes are:
#' - `mfrm_fit`: fitted model parameters and metadata.
#' - `mfrm_diagnostics`: fit/reliability/flag diagnostics.
#' - `mfrm_bias`: interaction bias estimates.
#' - `mfrm_dff` / `mfrm_dif`: differential-functioning contrasts and screening summaries.
#' - `mfrm_population_prediction`: scenario-level forecast summaries for one
#'   future design.
#' - `mfrm_unit_prediction`: fixed-calibration posterior summaries for future
#'   or partially observed persons.
#' - `mfrm_plausible_values`: approximate fixed-calibration posterior draws for
#'   future or partially observed persons.
#' - `mfrm_bundle` families: summary/report bundles and plotting payloads.
#'
#' @section Typical workflow:
#' 1. Prepare long-format data.
#' 2. Fit with [fit_mfrm()].
#' 3. Diagnose with [diagnose_mfrm()].
#' 4. Run [analyze_dff()] or [estimate_bias()] when fairness or interaction
#'    questions matter.
#' 5. Report with [build_apa_outputs()] and [build_visual_summaries()].
#' 6. For design planning, move to [build_mfrm_sim_spec()],
#'    [evaluate_mfrm_design()], and [predict_mfrm_population()].
#' 7. For future-unit scoring, refit or retain an `MML` calibration and then
#'    use [predict_mfrm_units()] or [sample_mfrm_plausible_values()].
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
#' **Joint Maximum Likelihood (JML/JMLE)**
#'
#' JMLE estimates all person and facet parameters simultaneously as fixed
#' effects by maximising the joint log-likelihood
#' \eqn{\ell(\boldsymbol{\theta}, \boldsymbol{\delta} \mid \mathbf{X})}
#' directly.  It does not assume a parametric person distribution, which
#' can be advantageous when the population shape is strongly non-normal,
#' but parameter estimates are known to be biased when the number of
#' persons is small relative to the number of items (Neyman & Scott, 1948).
#'
#' See [fit_mfrm()] for practical guidance on choosing between the two.
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
#' - Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#' - Linacre, J. M. (2002). What do Infit and Outfit, mean-square and
#'   standardized mean? *Rasch Measurement Transactions*, 16(2), 878.
#' - Masters, G. N. (1982). A Rasch model for partial credit scoring.
#'   *Psychometrika*, 47, 149--174.
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
#' Posteriori (EAP).  Joint Maximum Likelihood (JML/JMLE) estimates all
#' person and facet parameters simultaneously as fixed effects.
#'
#' MML is generally preferred for smaller samples because it avoids the
#' incidental-parameter problem of JML.  JML is faster and does not assume
#' a normal person distribution, which may be an advantage when the
#' population shape is strongly non-normal.
#'
#' See [fit_mfrm()] for usage.
#'
#' @examples
#' mfrm_threshold_profiles()
#' list_mfrmr_data()
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
