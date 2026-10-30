#' What mfrmr Calculates and Where to Read More
#'
#' @description
#' Find an analysis by its question, understand what its numbers mean, and
#' follow links to the calculation details and original references. For a
#' first runnable analysis, use [mfrmr_workflow_methods] or
#' `vignette("mfrmr-workflow", package = "mfrmr")` first. This is a help topic,
#' opened with `help("mfrmr_methods_reference")`, not a function to run.
#' For an integrated guide to **when to use each model and estimator**, read
#' `vignette("mfrmr-model-choice", package = "mfrmr")` or its
#' [online version](https://ryuya-dot-com.github.io/mfrmr/articles/mfrmr-model-choice.html).
#' It compares practical uses and limitations, including why a more complex
#' model is not automatically preferable. The sections below provide the
#' calculation references after that choice.
#'
#' @section Start with the question:
#' Suppose several judges score performances on several criteria. A low score
#' might reflect the performance, a stricter judge, a harder criterion, or a
#' combination the chosen model does not explain. MFRM estimates these effects
#' under a stated model. Raw agreement asks a different question: did judges
#' give the same ratings to the same performances?
#'
#' A **facet** is a modeled role, such as Rater or Task; a **level** is one
#' rater or task. Person is specified separately in `fit_mfrm()`, so
#' `facets = c("Rater", "Criterion")` gives three total facets including
#' Person. More facets do not mean more ability dimensions. The rating
#' assignments must contain information that distinguishes the requested
#' effects. Arbitrary names and flexible facet counts do not make every
#' crossed, nested or sparse design identifiable.
#' Before interpreting a diagram, read the score-direction and diagnostic-term
#' section of [mfrmr_visual_diagnostics]. It explains reverse-coded
#' comprehensibility, calibrated judges and fit, conditional interaction bias,
#' and adjacent-category Steps. A higher coded score is not inherently better.
#'
#' @section Check the rating design before fitting:
#' Start with [describe_mfrm_data()] and the intended comparison. A successful
#' data review is not a certificate of identification or adequate precision.
#'
#' 1. **Count the information, not only Persons.** Inspect ratings per Person,
#'    observations per facet level and shared ratings. Adding raters while each
#'    Person still receives two ratings spreads the same workload across more
#'    rater effects. More Persons do not repair an assignment that cannot
#'    separate the effects. There is no package-wide minimum sample size or
#'    facet-count rule that guarantees reliable estimates or intervals.
#' 2. **Separate the effects you want to compare.** For example, if Task A always
#'    uses Criterion A and Task B always uses Criterion B, their separate
#'    difficulties are confounded without additional information or constraints.
#'    Each Person-facet graph can be connected while this joint ambiguity
#'    remains. Review `design_connectivity`, then the estimator-specific
#'    identification and readiness results from [fit_mfrm()]. A local nonlinear
#'    rank check is not a proof of global identification.
#' 3. **Declare the assignments and rubric.** Supply `expected_design` to the
#'    data review when the planned roster is known. Compare planned cells,
#'    observed cells and retained rating rows; repeated events need their own
#'    interpretation. An unassigned cell is not a zero score. Set the intended
#'    category bounds and inspect category use within the selected step owner.
#'    A category observed elsewhere may still be unsupported for that owner.
#'    See [rating_scale_table()] and [fit_mfrm()] for the model-specific checks.
#' 4. **State what missingness and anchors mean.** Using available ratings does
#'    not correct ability-dependent assignment or score-dependent missingness.
#'    Missing population predictors have a separate `population_policy`.
#'    Shared Persons provide linking observations; they are not fixed anchors.
#'    [review_mfrm_anchors()] checks declared anchors and constraints. Treating
#'    an external estimate as an exact anchor does not propagate its uncertainty.
#'    Person anchors require JML; MML rejects them. By default,
#'    `centering_scope = "facet"` centers the whole declared facet, including
#'    fixed values. Use `"remainder"` only to reproduce legacy centering of
#'    remaining ungrouped levels. Read [mfrmr_linking_and_dff]
#'    for a complete example, the distinction from FACETS, and the references
#'    supporting each linking claim.
#' 5. **Check the requested output separately.** A fitted point, covariance,
#'    interval, residual analysis and G/D projection have different requirements.
#'    Flexible facet names/counts in one route do not extend another route.
#'    For example, two-family GPCM currently has exactly two non-Person facets;
#'    incomplete G-study observations do not authorize arbitrary unbalanced
#'    future D-study plans. Follow [gpcm_capability_matrix()],
#'    [mfrmr_interval_guide()] and [mfrm_d_study()] for those limits.
#'
#' @section Find the appropriate calculation:
#' Follow the help link for arguments, equations, returned fields, examples
#' and references. These are complementary analyses; running every row is not
#' a requirement for every assessment.
#'
#' - **Can my ratings support the intended comparison?** Retained ratings, category use,
#'   exposure and design/identification checks. Connectedness alone does not establish
#'   precision. [describe_mfrm_data()], [subset_connectivity_report()] and [fit_mfrm()];
#'   Linacre (1989).
#'
#' - **How do ability, severity and category use explain scores?** Ordered-category
#'   response probabilities with specified facet effects. RSM shares steps; PCM assigns
#'   steps to one selected facet; GPCM additionally allows positive slopes. [fit_mfrm()]
#'   and [gpcm_capability_matrix()]; Andrich (1978), Masters (1982), Muraki (1992).
#'
#' - **How are the parameters estimated?** MML integrates over a declared ability
#'   distribution; JML estimates Person parameters jointly with structural parameters.
#'   Native RSM/PCM CML conditions on each Person's total, eliminating ability
#'   without estimating a population distribution. It supplies structural
#'   conditional-information intervals and no ability population SD. Explicit
#'   comparison-only CML also supports identified interaction contrasts when
#'   other structural directions are unidentified; unavailable targets stay missing.
#'   These routes have different assumptions and uncertainty targets. The estimator and population
#'   sections of [fit_mfrm()]; Bock and Aitkin (1981), and the JML references on that
#'   page.
#'
#' - **How uncertain is a measure or difference?** A model-specific SE or interval for a
#'   named target, conditional on specified assumptions. A difference may need the
#'   covariance between estimates. [mfrmr_interval_guide()], [mfrm_facet_intervals()] and
#'   [mfrm_population_intervals()]; each page describes its calculation and limitations.
#'
#' - **Which ratings are unexpected?** Infit/Outfit and residual summaries compare
#'   ratings with model predictions; ZSTD transforms fit using the selected reference
#'   convention. [diagnose_mfrm()] and [fit_measures_table()]; Wright and Masters (1982),
#'   Linacre (2002). Extended models use [mfrm_response_diagnostics()].
#'
#' - **Are categories functioning as intended?** Usage counts, owner-specific step
#'   estimates and probability curves describe category use and transitions.
#'   [rating_scale_table()], [category_structure_report()] and
#'   [category_curves_report()]; Andrich (1978), Masters (1982).
#'   RSM has a common step set; PCM/GPCM has one declared step-owning facet,
#'   which can differ from the single slope-owning facet in GPCM. Review
#'   `owner_category_table` before interpreting pooled counts. Stored step
#'   SEs describe centered offsets, not full crossing, gap or curve uncertainty.
#'   Disorder does not automatically invalidate the model or justify collapsing
#'   categories. Two-family curves use their separate fitted-curve route with
#'   unavailable curve/step intervals; these ordinary helpers do not extend to it.
#'
#' - **Is a rater unusually strict in a particular context?** Conditional bias screens
#'   describe departures from the additive model for a combination of facet levels.
#'   [estimate_bias()] and [estimate_all_bias()]; Linacre (1989) and the FACETS sources
#'   in their help. These screens hold calibration fixed.
#'   Use [bias_interaction_report()] to keep the chosen review rules with saved
#'   tables and plots. Its scatter compares observed-minus-expected score
#'   averages with conditional logit shifts; neither its normal/profile bounds
#'   nor repeated views provide a calibrated fairness test.
#'
#' - **Are group differences specific to a facet?** DFF screens and separately matched
#'   MH/ordinal-logistic analyses use different adjustment models and sampling
#'   assumptions. [analyze_dff()], [mantel_haenszel_table()] and
#'   [ordinal_logistic_table()]; see their References and matching rules.
#'
#' - **How well are measures separated relative to error?** Separation, reliability and
#'   strata summarize estimated spread and measurement error for a facet. The
#'   interpretation and References sections of [diagnose_mfrm()] and
#'   [plot(x, type = "reliability")][plot.mfrm_fit]; Wright and Masters (1982).
#'
#' - **Do raters or scoring procedures agree?** Matched observed scores supply agreement
#'   coefficients, mean differences or ICCs. These do not estimate MFRM severity.
#'   [pairwise_kappa_table()], [gwet_ac2_table()], [krippendorff_alpha_table()],
#'   [pairwise_icc_table()] and [score_comparison_table()]; each has method-specific
#'   References.
#'
#' - **Do unexplained scores vary together?** Residual association concerns
#'   departures from the fitted model, not agreement of the original scores.
#'   [q3_statistic()] uses standardized residuals averaged within Person-by-level
#'   cells; it differs from Yen's raw-residual Q3. Its help identifies the
#'   original definition and the limits of the package's screening cutoffs.
#'   For organized residual patterns, [analyze_residual_pca()] distinguishes
#'   overall combined-facet columns from separate facet-level matrices.
#'   Eigenvalues across these different matrices do not rank facet importance.
#'   Its reference comparisons are exploratory; neither component counts nor
#'   fixed cutoffs establish a number of substantive abilities. Save the PCA
#'   result before repeated plots so the reference is not generated again.
#'
#' - **Do raters match an expert standard or show centrality?** Agreement with supplied
#'   references and concentration of category use are separate rater properties.
#'   [fit_rater_accuracy()] and [rater_centrality_table()]; their References identify the
#'   accuracy and centrality methods and current restrictions.
#'   Centrality's RSM residual correlation, supporting fit and rater-owned PCM
#'   threshold SD are indirect indices with separate availability and source
#'   status, not an estimated centrality parameter or a severity ranking.
#'   `plot(fit, type = "rater_severity", diagnostics = diagnostics)` instead
#'   orders facet locations and returns SE sources and interval-use labels.
#'   Its whiskers are normal approximations, not likelihood-profile intervals
#'   or pairwise tests. See [plot_rater_severity_profile()] for interpretation
#'   and [mfrmr_plot_migration] for the public `plot()` entry point.
#'   [rater_halo_network_analysis()] instead compares observed-score profiles
#'   across shared contexts. Its graph centrality is a network property;
#'   neither its edge filters nor high same-rater correlations establish
#'   causal halo or residual independence. See the rater-effects section below
#'   for the distinction from centrality and the horn effect.
#'
#' - **What does a cut score imply?** Posterior ability probabilities concern a Person
#'   now; classification accuracy/consistency concerns a specified future rating plan and
#'   population. [score_mfrm_persons()], [predict_mfrm_units()] and
#'   [classification_table()]; the latter gives Rudner (2001, 2005) and Lee (2010), with
#'   DOI links.
#'
#' - **Would another assessment plan be more dependable?** G-studies estimate
#'   observed-score variance/covariance components; D-studies project specified future
#'   plans and score weights. [mfrm_g_study()] and
#'   [mfrm_d_study()]; Brennan (2001) and the cited design-specific methods.
#'   These are separate from MFRM separation reliability.
#'
#' - **How do I retain the analysis and its caveats?** Saved results, figures, notes and
#'   exports preserve the calculation and its interpretation. [mfrm_results()],
#'   [mfrm_report()], [plot_data()] and [export_mfrm_results()]. See
#'   [mfrm_calibration_workflow] for eligible calibration reuse.
#'   Saving an analysis supports later review; a frozen calibration supports
#'   scoring new Persons at known facet levels; saving a score result supports
#'   later tables and plots without rescoring. These are different objects.
#'   Freezing makes a checked calibration available for scoring. It does not
#'   establish that the model or its prior suits another cohort. The calibration help links
#'   Bock and Mislevy (1982) to conditional EAP scoring, separately from
#'   calibration-estimation and interval-performance claims.
#'
#' @section Use results as evidence about a construct:
#' A construct is the attribute a score is intended to represent, such as
#' speaking proficiency. Begin with a claim about that score and its use:
#' "Content and Organization describe related aspects of performance across
#' the specified tasks and raters." State the population, rubric, score units,
#' intended composite or separate scores, and an expected pattern before
#' inspecting results. A different claim, such as diagnosing an individual's
#' weakness in Organization, needs evidence about that specific score or
#' difference, not just the dependability of a total.
#'
#' The *Standards* (AERA, APA, & NCME, 2014, Standards 1.13--1.14) require
#' evidence for the proposed internal structure, subscore interpretation and
#' composite rationale. Kane (2013) connects score-use claims to their
#' assumptions and supporting evidence. The following is a guide to applying
#' those principles with mfrmr, not a validity test implemented by the package.
#'
#' | Proposed claim | Inspect with mfrmr | Alternative explanation and warranted conclusion |
#' | --- | --- | --- |
#' | Related criterion scores reflect shared differences between people. | G-study Person covariance, its uncertainty, and the other source matrices; [plot.mfrm_multivariate_gstudy()]. | Positive Person covariance is compatible with the claim under the model. It does not identify one latent dimension. A near-zero or negative estimate with a broad interval leaves the proposed relation unresolved, rather than establishing distinct constructs. |
#' | The association is not solely a common scoring condition. | Rater, Task, Person-by-Rater and Person-by-Task covariances alongside Person covariance. | Rater covariance describes judges shifting both scores together; task covariance describes tasks shifting both. Interactions describe person-specific shifts. Shared rubric interpretation, genuine performance demands and impressions can compete as explanations. No covariance sign identifies halo or horn. |
#' | The rubric and model represent the intended ordering. | Owner-specific category counts/curves, Wright-map targeting, and matching Infit/Outfit or residual diagnostics; [rating_scale_table()], [plot.mfrm_fit()]. | Sparse categories, poor targeting or unexplained residuals can challenge the proposed measurement. Good fit alone does not establish the meaning of the trait; an unexpected pattern calls for inspection of the rubric and performances before recoding or excluding ratings. |
#' | A proposed score is dependable enough for its stated decision. | G/Phi and SEM for the prespecified score, population and future plan; [summary.mfrm_multivariate_d_study()]. | A heterogeneous Person population can raise a coefficient without reducing error in score units. High G/Phi concerns the declared replication, not construct validity. High MFRM rater separation means distinguishable severities, not agreement or evidence against shared scoring bias. |
#'
#' **Work through conflicting results.** If G is high but Person covariance
#' between the criteria is uncertain, the composite may be dependable while
#' the proposed shared relation remains unresolved. If Person covariance is
#' positive but Person-by-Task variation is substantial, consider whether the
#' intended generalization spans those tasks. A task-specific strength may be
#' substantive; it is error only relative to the chosen across-task target.
#' If fit looks satisfactory but the claimed proficiency is poorly represented
#' by the tasks, fit does not repair content underrepresentation.
#'
#' Inspect all components rather than choosing favorable ones. Covariances
#' depend on units, and a raw component is not its contribution to a future
#' average: D-studies apply the relevant count divisors. Review negative or
#' indefinite component estimates and interval qualifications. Separate
#' pointwise intervals do not test equality between components. An interval
#' containing zero neither proves independence nor establishes equivalence.
#'
#' **Choose reporting and follow-up from the claim.** Prespecify composite
#' weights for substantive reasons; optimizing G/Phi after inspecting the data
#' can change the score being interpreted. Dependable subscores do not by
#' themselves establish dependable or useful differences between them. MFRM
#' and multivariate G-theory provide complementary evidence, not a joint
#' fitted model or a latent-dimensionality test. Agreement between outputs
#' based on the same ratings is not independent replication.
#'
#' Report the claim, design, observed pattern with uncertainty, competing
#' explanation, justified conclusion and unresolved evidence. Where the claim
#' requires it, examine task content, how raters reach judgments, or relations
#' with an appropriate external measure. These require substantive evidence
#' beyond the returned tables; the package does not certify validity or
#' prescribe another analysis for every assessment. See the
#' [G-theory worked example](../articles/mfrmr-gtheory.html#examine-a-construct-claim)
#' and the [figure reading guide](../articles/mfrmr-visual-diagnostics.html#connect-a-figure-to-a-claim).
#'
#' @section Choose an agreement target and keep its denominator:
#' Start with the same rated response, identified by all relevant keys such as
#' Person, Task and Criterion. Matching several responses from one Person does
#' not make them independent sampling units. The runnable
#' `vignette("mfrmr-rating-agreement", package = "mfrmr")` shows matching,
#' exclusions, interpretation and saving without fitting a measurement model.
#'
#' - **Identical scores and model expectations:** [interrater_agreement_table()]
#'   reports exact agreement on matched contexts. Its `ExpectedExact` averages
#'   products of fitted category probabilities, with calibration fixed. This
#'   is a conditional model expectation, not the chance-agreement term in kappa
#'   or AC2. MML uses a plug-in EAP, not posterior integration. Unavailable
#'   primary estimates withhold expectations with an `ExpectedReason`;
#'   observed agreement is retained. Overall exact agreement counts every matched pair-context
#'   opportunity; a context rated by more judges supplies more pairs.
#' - **Category agreement for two raters:** [pairwise_kappa_table()] uses the
#'   same common complete units for both rater margins. QWK penalizes category
#'   differences quadratically. Numeric categories use their actual spacing;
#'   character categories use their declared positions. Each pair may use a
#'   different sample, so these coefficients are not a panel average.
#' - **Pooled agreement:** [krippendorff_alpha_table()] excludes singleton
#'   units from both disagreements and weights pairable units by rating count.
#'   Its ordinal distances depend on pairable-category frequencies.
#'   [gwet_ac2_table()] weights pairable units equally for observed agreement;
#'   singletons enter its expected-agreement margins. Its declared unused
#'   categories also matter. Choose the definition before seeing which value
#'   is larger; neither numerical missing-data convention corrects selective
#'   assignment or missingness.
#' - **Reliability on an interval score scale:** [pairwise_icc_table()] requires
#'   choices about selected versus randomly sampled raters, absolute agreement
#'   versus consistency, and one rating versus the mean of two. Consistency
#'   excludes systematic rater offsets; agreement includes them. The default
#'   concerns a single rating from selected raters. An average-rating pairwise
#'   ICC concerns two ratings, not the whole panel; numeric category codes
#'   alone do not justify an interval scale.
#' - **A score source against a supplied reference:** [score_comparison_table()]
#'   separates QWK, correlation and standardized mean difference. Adding a
#'   second reference uses common complete records for all statistics and can
#'   change the sample. A chosen reference is not automatically correct, and
#'   the optional Williamson et al. (2012) screen is an application-specific
#'   choice rather than a universal quality standard.
#'
#' The raw-score helpers above return descriptive point estimates with their
#' counts, cases, settings and notes; they do not provide sampling intervals.
#' Save their **complete lists** with [base::saveRDS()] and export chosen tables
#' with [utils::write.csv()]. They are not fitted-model result bundles for
#' [export_mfrm_results()], and their custom tutorial figures are not new
#' `plot()` methods. The separate `mfrm_interrater` bundle supports `summary()`
#' and `plot()`; see [interrater_agreement_table()] for its different repeated-
#' cell handling. Low residual association does not establish raw agreement,
#' and high agreement does not establish accuracy against an expert standard.
#'
#' @section Distinguish rater effects before choosing a diagnostic:
#' Consider a judge scoring content, organization and language on a four-category
#' rubric. These questions concern different features of the ratings:
#'
#' - **Severity:** does the judge generally give lower scores after accounting
#'   for the modeled Person and criterion effects? Review fitted Rater locations.
#' - **Centrality:** does the judge compress ratings toward middle categories?
#'   [rater_centrality_table()] combines the supplied RSM residual/fit and
#'   rater-specific PCM threshold descriptions. A severity near zero does not
#'   answer this question; neither does graph-theoretic centrality.
#' - **Halo or horn:** does a favorable or unfavorable impression spill over
#'   into ratings on other criteria? For example, a poor impression of language
#'   could lower content and organization ratings too. The latter pattern is
#'   called a horn effect; it is not defined by a negative criterion correlation.
#'   [rater_halo_network_analysis()] describes observed-score associations and
#'   cannot determine their psychological cause or classify halo versus horn.
#'   Genuine relationships among performance criteria also produce correlations.
#' - **Local dependence:** do ratings share variation beyond the chosen model?
#'   A supported [fit_mfrm_testlet()] model can account for a declared shared
#'   effect, but its variance does not identify halo or horn.
#'
#' Infit/Outfit summarize squared residuals. Low values can have several causes
#' and do not identify which effect occurred or its positive/negative direction.
#' Inspect shared-rating counts, category use and the actual performances before
#' drawing conclusions about a rater. Agreement among several indices computed
#' from the same ratings is not independent confirmation.
#'
#' Direct models answer a further question: how would estimated measurements
#' change when a specified rater effect is included in the response model?
#' The FM-SC severity/centrality model and the MRFM-H halo mixture models are
#' **not implemented**. See the original references and model distinctions in
#' [rater_centrality_table()] and [rater_halo_network_analysis()]. These models
#' require their own identifying assumptions and uncertainty evaluation;
#' adding GPCM slopes or a testlet variance does not implement them.
#'
#' @section Choose the G-theory target and uncertainty:
#' Start with the rating design and the intended score, not the word
#' "reliability" alone. The multivariate route also accepts **one score**;
#' multiple score columns are fixed score components, not additional random
#' facets or latent dimensions. All methods below require explicit
#' `assumption = "normal"`, prespecified scores/weights and the documented
#' sampling/assignment assumptions. An incomplete source is not permission
#' to use a sparse future plan: current D-studies project complete balanced plans.
#'
#' First choose the quantity whose uncertainty you need:
#'
#' | Question | Interface | Meaning of the uncertainty |
#' | --- | --- | --- |
#' | How uncertain is a source variance or covariance? | [confint.mfrm_multivariate_gstudy()] | Component interval; variance units are squared score units, covariance units are products of score units. |
#' | How uncertain are the original component estimates jointly? | [vcov.mfrm_multivariate_gstudy()] | Sampling covariance and SEs of ANOVA/MINQUE components, including identifiable incomplete sources. Not covariance between individual ratings. |
#' | How uncertain is a future plan's G, Phi or SEM? | [confint.mfrm_multivariate_d_study()] | One selected population coefficient or measurement-error target at prespecified counts. Coefficient SE and measurement SEM are different quantities. |
#' | How uncertain is the change between two plans? | `confint(d, reference = 1, ...)` or [mfrm_multivariate_d_compare()] | Paired comparison using the shared source; never subtraction of separate marginal endpoints. |
#'
#' Then select a method allowed by the **observed source**. Here `g` is a
#' G-study and `d` its D-study. The table describes available calculations,
#' not a ranking by interval width or a guarantee of nominal coverage.
#'
#' | Target and source | Default or explicit choices | Qualification |
#' | --- | --- | --- |
#' | Component variance, balanced ANOVA | `confint(g)`: `auto`/`mls` | Exact residual chi-square; other component MLS intervals are approximate. |
#' | Signed covariance between two scores, balanced ANOVA | Explicit `gci` or `regression_region` | Generalized approximate coverage or conservative enclosure; no automatic covariance default. |
#' | One-facet G/Phi, balanced ANOVA | `confint(d)`: `exact_f` for G, `gci` for Phi | Normal-model exact F versus generalized approximate coverage. |
#' | Two-crossed-facet G/Phi, balanced ANOVA | Default `mls`; explicit `delta`, `approximate_f`, `mean_square_region`, `reml_profile`; G also permits `score` | Distinct approximations/conservative enclosures. Read formula-completion flags for MLS; these are not interchangeable calibrations. |
#' | Nested G/Phi, balanced ANOVA | With unchanged source child count per parent, the one-facet defaults apply. Otherwise explicitly select `delta`, `approximate_f`, `mean_square_region` or `reml_profile`. | Both facets remain random. The default reduction is unavailable when the child count changes. |
#' | Relative/absolute SEM, balanced crossed/nested ANOVA | Default `mls`; two-facet sources also permit `mean_square_region` | Error-variance intervals transformed by square root; exactness depends on the contributing mean squares. |
#' | G/Phi, identifiable MINQUE source | Explicit `delta` or `reml_profile`; one-facet G also permits `spectral_f` | Delta concerns the original coefficient; REML has a separate center. Spectral F has an exact normal-model pivot with its own residual-information requirements. |
#' | Component variance/covariance or SEM, identifiable MINQUE source | Explicit `reml_profile`; `residual_contrast` additionally supplies residual components and one-facet RelativeSEM | REML is experimental; residual contrasts have an exact variance law or conservative signed covariance enclosure. Residual-only inference does not supply other components or AbsoluteSEM. |
#'
#' `delta` is also an explicit G/Phi option for balanced one-facet sources;
#' `reml_profile` accepts every supported balanced source as well as identifiable
#' incomplete/unequal sources. `residual_contrast` and one-facet `spectral_f`
#' also accept retained complete ratings. These methods do not become an
#' automatic fallback when another calculation fails. MINQUE requests require
#' an explicit interval method even when the source happens to be balanced.
#' For derivations, required inputs, optional dependencies and references, use
#' [confint.mfrm_multivariate_gstudy()] and
#' [confint.mfrm_multivariate_d_study()]. Almehrizi's (2021) coefficient-SE
#' formula and the limits of its extension are in the latter's *Delta coefficient
#' SEs and intervals* section; SE accuracy does not establish Wald coverage.
#'
#' For **paired plans**, choose the relevant source law separately:
#'
#' - Balanced one-facet G/Phi map the whole common F/generalized set. The nested
#'   G/Phi reduction requires every future child count to equal the source count.
#'   One-facet SEM and nested SEM with a common future child count map one
#'   error-SD interval; that future child count may differ from the source.
#' - General balanced two-facet comparisons have explicit `mean_square_region`
#'   enclosures from one shared source region. These retain dependence but can
#'   be wider than the tight difference image. Two crossed random facets also
#'   permit [mfrm_multivariate_d_compare()] for a separately qualified paired
#'   delta approximation, including identifiable MINQUE sources.
#' - Explicit `reml_profile` supplies G/Phi/SEM comparisons for supported
#'   complete or incomplete sources. Proportional one-facet/common-child nested
#'   plans map a common profile; other plans profile the difference directly.
#'   Calibration remains approximate, especially at zero-contrast singularities.
#' - For retained one-facet ratings, `spectral_f` G and `residual_contrast`
#'   RelativeSEM differences map their shared confidence sets. Marginal
#'   `method = "delta"` in `confint(d)` does not supply paired differences.
#'
#' **Read and save the result.** [summary.mfrm_multivariate_d_study()] with
#' `intervals = TRUE` collects supported methods, original points, any separate
#' `IntervalCenter`, source/future counts and weights. Select `parm` explicitly
#' when the chosen method supports only some metrics. Ordinary `summary(d)`
#' returns the coefficient table. `confint(d, method = c(...))` compares compatible
#' methods for one target; `confint(g)` selects one component method at a time.
#'
#' Wide computable intervals remain results. Delta limits outside `[0,1]`,
#' negative original components, empty sets, failed numerical checks and absent
#' methods are different states. Read `Status`, availability and saved notes.
#' A REML or residual-projection center must not be relabeled as MINQUE.
#' General bootstrap intervals are not implemented. Nonnormal robustness,
#' selection-adjusted inference and sparse/unequal future plans are not supplied.
#' These limits do not make an implemented interval calculation nonexistent.
#'
#' There is no universal acceptable width or G/Phi cutoff. Define the intended
#' score use, a meaningful improvement or threshold, and task/rater costs before
#' choosing a plan. An interval spanning a decision threshold describes
#' uncertainty about that decision; selecting the narrowest method after seeing
#' its result does not resolve it. A method comparison on one dataset is
#' sensitivity analysis, not a coverage study. Neither high G/Phi nor a narrow
#' interval alone establishes construct validity or fairness.
#'
#' Use `plot(d)` for point curves, `intervals = TRUE` for an explicit interval
#' request, or supply saved intervals for reuse. `plot(ci)` compares D-study
#' intervals on one axis; `plot(g, intervals = TRUE, ...)` shows component
#' intervals and retains them in the returned plot object. Save complete
#' results with `saveRDS()`; `as_ggplot()` renders saved plot data without fitting.
#' Ribbons are pointwise, not simultaneous bands; overlap is not a paired test.
#' The [G-theory tutorial](../articles/mfrmr-gtheory.html) supplies a bundled-data
#' reading path, and the function help holds mathematical details and references.
#' The simpler [mfrm_generalizability()] / [mfrm_d_study()] route is a separate
#' main-effects approximation; its residual sensitivity is not a confidence interval.
#'
#' @section Read the method before interpreting the number:
#' **Model:** RSM/PCM use unit discrimination. One-family GPCM assigns slopes
#' to one selected facet; MML allows a different step owner, whereas JML
#' requires the same owner. The experimental two-family MML route has exactly
#' two non-Person facets and additional restrictions. Check
#' [gpcm_capability_matrix()] before transferring an ordinary-model analysis.
#' The original response-model papers do not establish the performance of
#' every mfrmr estimator, extension or interval.
#'
#' **Estimation and scoring:** MML uses numerical integration of the marginal
#' likelihood. JML includes Person parameters in the joint likelihood. EAP
#' scoring subsequently averages over a conditional ability distribution;
#' its posterior SD is not a sampling SE for a structural parameter. Ordinary
#' RSM/PCM MML fixes N(0,1) when `population_formula = NULL`; estimating a
#' population model is an explicit choice. A histogram of EAP estimates is
#' not itself an estimate of the latent population density. Numerical
#' integration, population assumptions and the source of uncertainty are
#' documented separately in [fit_mfrm()] and the scoring help pages.
#'
#' **Fit:** Infit weights squared standardized residuals by their model
#' variance; Outfit gives each retained residual equal weight when observation
#' weights are all one. Values above one indicate more residual variation than
#' the reference model expects. A flag is a reason to inspect ratings and
#' assumptions, not a verdict on a rater. Read the selected residual basis,
#' weighting and ZSTD convention in [diagnose_mfrm()]. Posterior diagnostics
#' for extended models do not inherit ordinary fit cutoffs.
#'
#' **Reliability and strata:** the ordinary summaries adjust the spread of
#' finite measures using their matching squared SEs. Invalid or missing SEs
#' can make the summary unavailable. High rater separation indicates
#' distinguishable severities, not agreement; low separation can also reflect
#' poor precision. Strata is a hypothetical separation index, not a discovered
#' number of groups or a set of pass/fail boundaries. Model and fit-adjusted
#' summaries are not confidence bounds. See the
#' [Facets explanation](https://www.winsteps.com/facetman64/reliability.htm)
#' alongside mfrmr's estimator-specific conventions in [diagnose_mfrm()].
#' These marginal-SE summaries omit covariance among calibrated measures.
#' Applying variance subtraction to EAP estimates is not posterior-variance
#' EAP reliability, and neither quantity is a G-theory coefficient.
#'
#' **Post-fit shrinkage:** [apply_empirical_bayes_shrinkage()] adds a
#' descriptive zero-centered normal-model adjustment while retaining the
#' fitted estimates. Its plug-in SD conditions on the chosen prior variance;
#' it is not a calibrated sampling SE or a general improvement guarantee.
#' Additional Person shrinkage requires ordinary JML sampling SEs. MML EAP
#' estimates and posterior SDs must not be fed through a second such update.
#'
#' **Categories and missingness:** use the intended rubric, not only the
#' observed score range. `category_policy = "preserve"` retains that ladder;
#' `"collapse"` changes a gapped ladder and records the mapping. Missing
#' assigned scores, unassigned ratings and unavailable rubric categories are
#' different situations. The current ordinary models use a common category
#' range; owner-specific ranges are not implied by PCM step ownership. See
#' [describe_mfrm_data()] and [fit_mfrm()] for support checks and retained rows.
#'
#' **Uncertainty:** distinguish the estimate, its numerical status and the
#' interpretation allowed by its assumptions. A finite experimental interval
#' is not evidence of nominal coverage. A missing endpoint is not infinity.
#' Conditional bias screens and fixed-calibration Person intervals omit
#' calibration uncertainty. Use [mfrmr_interval_guide()] to select the target
#' before choosing a method.
#'
#' @section Find and cite the supporting literature:
#' Open the function's **Details**, named method sections and **References**;
#' report its selected settings as well as the paper. For example,
#' `help("diagnose_mfrm")` connects the Infit/Outfit definitions to the ZSTD
#' conventions and reliability formulas; `help("classification_table")`
#' explains how the Lee recursion differs from the Rudner approximation.
#' The foundational sources below identify the model/estimation framework.
#' Specialist help pages provide the references for their own formulas and
#' adaptations. A citation to a general framework is not validation of a new
#' implementation or of your data's assumptions.
#'
#' Use `citation("mfrmr")` for the software citation and
#' `packageVersion("mfrmr")` for the installed version. Cite the relevant
#' methodological sources as well. A DOI identifies a paper even if access to
#' its full text requires a library subscription.
#'
#' @references
#' American Educational Research Association, American Psychological
#' Association, & National Council on Measurement in Education. (2014).
#' *Standards for educational and psychological testing*. American Educational
#' Research Association. Internal structure and score-combination rationale:
#' Standards 1.13--1.14, pp. 26--27; reliability/precision, Chapter 2.
#' \url{https://www.testingstandards.net/uploads/7/6/6/4/76643089/standards_2014edition.pdf}.
#'
#' Kane, M. T. (2013). Validating the interpretations and uses of test scores.
#' *Journal of Educational Measurement, 50*(1), 1--73.
#' \doi{10.1111/jedm.12000}. A framework for claims, assumptions and evidence;
#' not validation of a particular mfrmr calculation.
#'
#' Andrich, D. (1978). A rating formulation for ordered response categories.
#' *Psychometrika, 43*, 561-573. \doi{10.1007/BF02293814}. RSM foundation.
#'
#' Masters, G. N. (1982). A Rasch model for partial credit scoring.
#' *Psychometrika, 47*, 149-174. \doi{10.1007/BF02296272}. PCM foundation.
#'
#' Muraki, E. (1992). A generalized partial credit model: Application of an
#' EM algorithm. *Applied Psychological Measurement, 16*(2), 159-176.
#' \doi{10.1177/014662169201600206}. GPCM foundation.
#'
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation
#' of item parameters: Application of an EM algorithm. *Psychometrika, 46*,
#' 443-459. \doi{10.1007/BF02293801}. Marginal-likelihood estimation framework;
#' this citation does not mean every mfrmr engine uses EM.
#'
#' Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#' Many-facet formulation and interpretation.
#'
#' Wright, B. D., & Masters, G. N. (1982). *Rating scale analysis*.
#' MESA Press. Fit, separation and reliability conventions.
#'
#' Linacre, J. M. (2002). What do Infit and Outfit, Mean-square and
#' Standardized mean? *Rasch Measurement Transactions, 16*(2), 878.
#' \url{https://www.rasch.org/rmt/rmt162f.htm}. An introductory fit explanation.
#'
#' Brennan, R. L. (2003). *Coefficients and indices in generalizability theory*
#' (CASMA Research Report No. 1). University of Iowa. Pages 3--4 introduce
#' random-facet G/D error and coefficient targets.
#' \url{https://education.uiowa.edu/sites/education.uiowa.edu/files/2026-04/casma-research-report-1-archived.pdf}.
#'
#' @seealso [mfrmr_workflow_methods], [mfrmr_output_guide()],
#'   [mfrmr_interval_guide()], [gpcm_capability_matrix()]
#' @concept calculation methods
#' @concept methodological references
#' @name mfrmr_methods_reference
NULL
