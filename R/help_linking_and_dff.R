#' mfrmr Linking and DFF Guide
#'
#' @description
#' Package-native guide to checking connectedness, building anchor-based links,
#' monitoring drift, and screening differential facet functioning (DFF) in
#' `mfrmr`.
#'
#' @section What does an anchor do?:
#' Suppose a previous calibration placed Rater R01 at 0.35 logits. A direct
#' anchor tells the new fit to keep that value exactly; the new ratings do not
#' re-estimate it. This requires the same rater identity, model orientation and
#' a defensible stability assumption, not just the same text label.
#'
#' Keep four ideas separate:
#' - **Common ratings:** overlapping Persons, tasks or raters connect observed
#'   ratings. Their parameter values need not be known or fixed.
#' - **Scale identification:** a constraint such as a zero mean sets an origin
#'   for otherwise interchangeable coordinates. It does not establish invariance.
#' - **Direct anchors:** `Facet`, `Level`, `Anchor` fixes a location. For
#'   example, R01 = 0.35 means exactly 0.35 on the internal logit scale.
#' - **Group anchors:** `Facet`, `Level`, `Group`, `GroupValue` fixes a mean.
#'   For R01 and R02 in a group with target 0.20, their locations must sum to
#'   0.40; neither location is individually fixed unless also directly anchored.
#'   This is an unweighted mean across declared levels, including fixed members.
#'
#' @section What the current implementation supports:
#' Non-Person location constraints are available in supported ordinary MML/JML
#' fits. Person direct/group constraints require JML; MML integrates ability
#' and rejects them. This is an implementation boundary, not a claim that
#' fixed-Person MML is mathematically impossible. Two-family GPCM has separate
#' restrictions; see [gpcm_capability_matrix()].
#'
#' With `centering_scope = "facet"` (the ordinary-fit default), centered facets
#' sum to zero across all levels, including direct anchors and group totals.
#' A single ungrouped remainder therefore balances those totals. Use
#' `centering_scope = "remainder"` to reproduce earlier fits that centered only
#' the levels outside direct/group constraints. Ordinary additive RSM/PCM JML
#' also offers `centering_scope = "anchor_origin"`: direct or group anchors replace
#' their facet's centering, with an unanchored, active noncenter facet.
#' Group targets are repeated arithmetic means, not per-element contributions.
#' The default remains `"facet"`. Read [fit_mfrm()] before linking
#' a cohort with new raters or tasks; population and `noncenter_facet` settings
#' also determine the scale. Contradictory constraints stop before estimation.
#'
#' [anchor_to_baseline()] transfers selected facet locations, not a whole
#' calibration: thresholds/slopes are re-estimated, and most settings must be
#' supplied again. RSM/PCM MML defaults to a fixed standard-normal population;
#' anchoring does not automatically estimate a changed ability distribution.
#' Use [extract_mfrm_calibration()] and [score_mfrm_calibration()] to score against a complete
#' fixed calibration within that separate workflow's supported scope.
#'
#' @section Annual assessments with common recordings:
#' Suppose an English speech contest has new contestants each year, mostly
#' the same judges, and a common three-minute format. Judges also rate the
#' same recorded speeches each year. These recordings can connect annual
#' ratings if their performance meaning and scoring conditions remain
#' comparable. Returning judges or unchanged criterion names alone do not
#' establish that their severity or interpretation stayed constant.
#'
#' Give the same recording the same performance ID across years. Give a new
#' speech a new performance ID, even if its speaker competed before. Keep
#' judge identity and rating year/session separately: the question is how a
#' judge's severity changes, not whether that judge is the same person.
#' Contestant-group differences are distinct from individual growth.
#'
#' Two analysis strategies have different assumptions:
#' - **Joint calibration:** share an unknown performance measure for each
#'   common recording across years and estimate judge-year effects. The
#'   recording values need not be known or passed to `anchors`. Shared
#'   unknown measures still require an identified model and common scale.
#' - **Fixed-reference calibration:** fix previously calibrated recording
#'   measures when fitting a later year, while estimating its judges' effects.
#'   When recordings occupy the Person role, those direct anchors require
#'   JML in the current implementation. The receiving fit treats their values
#'   as exact; their source estimation error is excluded.
#'
#' Do not fix the returning judges' severities merely because they return.
#' Also, centering the panel to the same mean separately each year can remove
#' the panel-wide change being studied. Annual rater effects, criterion/step
#' invariance and population assumptions must be specified together. Separate
#' fixed-standard-normal MML fits do not estimate changing cohort means.
#'
#' Current helper output answers narrower questions. In an anchored result,
#' `$drift` compares retained fixed non-Person values with their supplied
#' values; zero differences do not show stability. In [detect_anchor_drift()],
#' `$drift_table$Drift` subtracts `LinkOffset`. If judges are the linking set,
#' a common shift in all their estimates is removed; the remainder cannot
#' establish that the panel was stable. Nor does `LinkOffset` alone identify
#' a real panel change, since it can also reflect different coordinates.
#'
#' Preserve individual ratings before discussion or disclosure of reference
#' scores, training phase, and video/live mode. Check coverage across performance
#' levels and criteria, memory of familiar scores, shared assignments and
#' dependence among repeated ratings. Training agreement alone does not
#' demonstrate stable operational scoring. A complete common-recording
#' workflow and formal year-difference intervals are not yet verified in mfrmr;
#' conditional anchor SEs do not supply cross-year covariance or source-anchor
#' uncertainty. The companion linking vignette develops this design example
#' and explains what the existing returned tables can establish.
#'
#' @section Two schools connected by common recordings:
#' Suppose School A has 30 students and five raters, School B has 40 students
#' and three raters, and both use three tasks. If each rater scores only their
#' own school's students, common tasks alone do not separate a school-specific
#' shift in student ability from an equal shift in its panel's severity.
#' Having all eight raters independently score the same three recordings can
#' supply that missing link. In a joint calibration, the common performances
#' share IDs and unknown measures; their scores need not be known anchors.
#' Unique local student IDs and distinct rater IDs prevent accidental overlap.
#' Record the task actually performed in each recording; do not replicate a
#' single recording across tasks that were never observed. Define whether
#' the measured unit is a student across tasks or a particular performance.
#'
#' This supports a model-based common scale under the linking assumptions,
#' not automatic interchangeability of school raw scores or a causal school
#' effect. Three recordings are three shared performances, not 24 independent
#' linking performances merely because eight raters assess each. Their ability
#' range, category coverage and stability across panels matter. Review whether
#' judgments were independent, reference scores were disclosed, or recording
#' versus live presentation changes severity. A connected graph is necessary
#' evidence to inspect, not proof of identification, precision or invariance.
#'
#' Specify each rater-by-school comparison separately. A rater's A-minus-B
#' effect relative to the entire panel is different from the difference
#' between two raters' A-minus-B effects. With incomplete school coverage,
#' the second comparison can be identified even when the first is not.
#' State the reference raters and their weights; dropping an unobserved rater
#' changes a whole-panel comparison. Under PCM, also state the step owner:
#' a location contrast alone need not describe every adjacent-category
#' response contrast when the compared raters own different step ladders.
#'
#' Preserve a recording's school of origin when another school rates it.
#' Shared neutral recordings can connect the scale without supplying the
#' missing A-versus-B observations for each rater. A common latent-population
#' assumption can supply information absent from an unrestricted-Person
#' comparison; allowing school-specific population means can restore that
#' ambiguity. These are different assumptions, not evidence that MML always
#' resolves an unidentified JML comparison. Use [describe_mfrm_data()] with
#' `comparisons` to review explicitly weighted RSM/PCM interaction contrasts
#' before fitting. Its conditional-map results and MML limitations remain
#' separate from marginal identification and from residual/refit DFF effects.
#' Its `IdentificationBasis` column distinguishes a target identified in the
#' unrestricted-Person response model from one whose conditional identification
#' requires the declared population restrictions. Neither makes MML
#' distribution-free or proves global marginal identification. CML instead
#' records identification from conditional support at the observed totals.
#' [analyze_dff()] can retain that review as `comparison_review` after checking
#' the retained responses and origin labels; it does not automatically certify
#' the identification of its own residual/refit contrasts.
#' For a matching joint interaction fit, pass the review to
#' [interaction_effect_table()] or `summary(fit, comparison_review = review)`.
#' These return descriptive weighted interaction estimates after full
#' conditional-rank, model/data alignment and retained-point numerical checks.
#' Unavailable rows retain reasons and separate optimizer traces. Positive
#' MML conditional ranks still do not certify global marginal identification;
#' no partial fit is run by extraction. For target-only identification, explicitly
#' use `fit_mfrm(fit_scope = "comparisons", comparison_review = review)` with the
#' review's model arguments. Its separate comparison result preserves the observed
#' likelihood and invariant targets, checks JML recession and numerical adequacy,
#' and retains unidentified/failed targets as missing. MML conclusions remain
#' conditional/local; JML/MML comparison fits supply no formal uncertainty.
#' Declaring `method = "CML"` in the review enables conditional-information
#' SEs and pointwise 95 percent Wald intervals for identified RSM/PCM contrasts.
#' `mfrm_facet_intervals(fit)` reuses these saved comparisons for interval plots
#' and APA tables, retaining unidentified rows without estimates. Its `level`
#' argument changes the normal critical value without refitting.
#' Identification uses conditional support at the actual Person totals, and
#' finite-MLE and information checks determine availability. Unidentified targets remain
#' missing. This route needs no ability distribution and estimates no population
#' SD; its asymptotic intervals do not establish finite-sample coverage or a
#' calibrated DFF classification. See [fit_mfrm()] and the linking vignette.
#'
#' For matching full RSM/PCM MML fits, [mfrm_facet_intervals()] with
#' `comparison_review` adds model-based interaction intervals after joint
#' information and finer-grid checks. `scale = "standardized"` divides by
#' the same fit's estimated intercept-only normal population SD, with the full
#' joint delta method. That common-population assumption differs from `~Group`;
#' the latter supports native intervals but its residual SD is not treated as
#' marginal population SD. Numerical eligibility does not certify global
#' marginal identification or finite-sample coverage.
#'
#' GPCM comparisons also require the slope owner and ability scale. A positive
#' location difference does not ensure that the expected-score difference has
#' the same sign across ability: unequal slopes can produce crossing response
#' curves. Specify the ability range and compare category probabilities or
#' expected scores with the actual step ladders. One-family slopes have
#' geometric mean one; under default MML the population standard deviation
#' supplies the common discrimination scale. Fixing both that standard deviation
#' and the slope geometric mean to one imposes a narrower model.
#'
#' With two slope families, distinguish each observed effective slope product
#' from its individual facet components. Disconnected facet crossings can
#' preserve observed products while reversing a component comparison; the
#' current fit route rejects that design. Predicting an unobserved crossing
#' depends on the multiplicative slope model. A full-rank product map alone
#' does not certify the full marginal response model, and a local numerical
#' rank check is not a global identification or precision guarantee. The
#' supported two-family route uses fixed-standard-normal MML with exactly two
#' measurement facets; it does not supply a rater-by-school interaction model.
#'
#' Prefer the joint route when no externally established reference values are
#' available. Fixing previously estimated recording values in the Person role
#' instead uses the JML-only anchor route and excludes their source uncertainty.
#' Do not force each school's rater panel to have the same mean to manufacture
#' a link, or add School as a free additive response facet alongside unrestricted
#' student effects without checking identification. A population school-mean
#' contrast needs its own population model and handling of common recordings
#' outside those cohorts; it is not an ordinary comparison of EAP means.
#' The existing common-recording workflow has not yet been verified end to end
#' for this school design or its school-difference uncertainty.
#'
#' @section Can I obtain intervals when I use anchors?:
#' For ordinary unit-weight RSM/PCM MML with fixed standard-normal ability and
#' fixed-grid integration, [mfrm_facet_intervals()] can return model-information
#' or cluster-sandwich intervals for free facet locations and contrasts, subject
#' to its checks. These intervals condition on the supplied anchors. A fixed
#' location has SE zero and no inferential bounds; its source uncertainty has
#' not been estimated by the receiving fit.
#'
#' Estimating the population mean and variance with anchors is a different
#' configuration. Its point estimates are supported, but the native population
#' and location interval helpers currently exclude anchors and constraint
#' groups. An `EstimabilityState` of `"not_evaluated"` denotes an incomplete
#' identification check, not failed optimization or proof of nonidentifiability.
#' Inspect `fit$readiness$fit` for numerical status and
#' `fit$data_review$estimability$nonlinear_local_estimability` for the separate
#' local-rank result. Passing that local check alone does not authorize intervals
#' or anchor export for reuse. Do not choose a fixed population solely to obtain
#' intervals when its assumptions are inappropriate.
#'
#' One-family GPCM MML has a separate slope/curve route: [confint.mfrm_fit()]
#' with `parm = "slopes"` and [mfrm_curve_intervals()] can condition on fixed anchors when their
#' own checks pass. Anchored native-location/profile-slope intervals and portable
#' GPCM extraction remain unavailable. An available interval does not certify
#' the separate integration checks required by [predict_mfrm_units()].
#'
#' @section Reusing an analysis saved before the centering change:
#' - **Read the original result:** `readRDS()` preserves its coefficients and
#'   constraints. `summary(fit, compute = "never")$settings_overview` and
#'   [specifications_report()] identify the centering convention. A native fit
#'   without convention metadata uses the legacy `"remainder"` interpretation.
#' - **Repeat the original analysis:** [build_mfrm_replay_script()] embeds its
#'   anchors and retains that convention. Running the script refits the model;
#'   it does not merely read a saved result.
#' - **Adopt whole-facet centering:** refit from the original responses with
#'   `centering_scope = "facet"` and the intended anchor, population, rubric,
#'   uncentered-facet and owner settings. Do not edit `fit$config` or recenter
#'   an output table to convert the fit. Keep the old and new analyses separately.
#' - **Refresh derived output:** recompute interval objects from the unchanged
#'   fit for revised fixed-contrast derivatives or anchor explanations. Rerun
#'   affected [analyze_dff()] subgroup analyses to retain target/Person anchors;
#'   this does require subgroup refitting. Loading old interval/DFF objects
#'   does not update them.
#' - **Keep a frozen calibration:** [load_mfrm_calibration()] preserves its
#'   original labels and values. Re-extract from the source fit, satisfying the
#'   current validation requirements, for explicit centering labels. If only
#'   the fit survives without the required RSM/PCM MML quadrature review,
#'   obtaining that review needs the original data and additional grid fits.
#'   If only
#'   the calibration file survives, do not infer a missing convention or relabel
#'   it; keep its original scoring interpretation. Changing its calibration
#'   requires the source analysis and a new validated, frozen artifact.
#'
#' The linking vignette provides a runnable migration example. Source-anchor
#' uncertainty remains excluded; a new label or file format adds no uncertainty
#' calculation. Older MML analyses containing Person anchors need substantive
#' review: those values were ignored, whereas current fitting rejects them.
#' Merely selecting the legacy centering mode does not repair such an analysis.
#'
#' @section Comparing anchors with TAM and ConQuest:
#' The same anchor number does not necessarily fix the same comparison.
#' Check the response model, parameter coding and remaining constraints before
#' transferring values between programs. A parameter index can change when
#' the model or retained facet levels change; match the actual design and labels,
#' not just the row number in an earlier output.
#'
#' TAM's `xsi.fixed` fixes indexed coefficients of its response design, which
#' need not be standalone facet locations. ConQuest's `anchor_parameters`
#' (`anchor_xsi`) serves that role. Population regression and covariance have
#' separate controls: TAM's `beta.fixed`/`variance.fixed` and ConQuest's
#' `anchor_beta`/`anchor_sigma`. Slopes also require their own coding; neither
#' an mfrmr location table nor an arbitrary list of external coefficients
#' represents a complete calibration.
#'
#' For a concrete RSM example, fix one of four raters at 0.70 and leave the
#' others outside group constraints. Whole-facet centering makes the other
#' three raters sum to -0.70: R01 is about 0.933 above their mean. This matches
#' the rater constraint in TAM's `tam.mml.mfr(formulaA = ~ item + rater + step,
#' constraint = "cases")` and ConQuest's `model task+rater+step` with
#' `lconstraints=cases`. Those settings also leave task locations uncentered;
#' use the corresponding `noncenter_facet` and match the population explicitly.
#' mfrmr's legacy `"remainder"` mode instead makes the other raters sum to zero,
#' putting R01 0.70 above their mean. Changing labels cannot reconcile these
#' restrictions. Group means likewise require matching linear constraints,
#' not copying them as individual fixed values.
#'
#' Match the population mean/variance, scoring, estimator and integration
#' accuracy as well. Equal numbers of integration nodes need not mean equal
#' rules or accuracy. ConQuest documents `anchor_theta` as ignored under MML;
#' mfrmr instead rejects Person anchors there. TAM's `tam.jml()` documents
#' `theta.fixed` for `version = 1`. Common Persons and fixed Person abilities
#' remain distinct concepts in all these comparisons.
#'
#' Agreement of point estimates does not establish agreement of SEs. Compare
#' uncertainty for the same free parameter or contrast, using the appropriate
#' covariance transformation. In ordinary fixed-anchor estimation, uncertainty
#' is conditional on the supplied values; a zero or missing SE for an anchor
#' does not account for its source calibration error. See the companion
#' vignette for a comparison table and the official manuals below.
#'
#' @section Start with the linking question:
#' - "Is the design connected enough to support a common scale?"
#'   Use [subset_connectivity_report()] and `plot(..., type = "design_matrix")`.
#' - "Which elements can I export as anchors from an existing fit?"
#'   Use [make_anchor_table()] and [review_mfrm_anchors()].
#' - "How do I anchor a new administration to a baseline?"
#'   Use [anchor_to_baseline()].
#' - "Have common elements drifted across separately fitted waves?"
#'   Use [detect_anchor_drift()] and [plot(x)][plot.mfrm_anchor_drift].
#' - "Can I synthesize anchor review, drift, and chain evidence into one review?"
#'   Use [build_linking_review()].
#' - "Do specific facet levels function differently across groups?"
#'   Use [analyze_dff()], [plot(x, type = "heatmap")][plot.mfrm_dif], and [plot(x, type = "summary")][plot.mfrm_dif].
#'
#' @section Recommended linking route:
#' Before fitting, use [describe_mfrm_data()] on the rating table and inspect
#' `design_connectivity`. Supply `expected_design` if the planned assignment
#' is available, so unassigned cells and missing assigned ratings are distinct.
#' A connected graph alone does not establish identification of all model
#' parameters, adequate precision, or absence of informative assignment.
#'
#' 1. Fit with [fit_mfrm()] and diagnose with [diagnose_mfrm()].
#' 2. Check connectedness with [subset_connectivity_report()].
#' 3. Build or review anchors with [make_anchor_table()] and
#'    [review_mfrm_anchors()].
#' 4. Use [anchor_to_baseline()] when you need to place raw new data onto a
#'    baseline scale.
#' 5. Use [build_equating_chain()] only as a screened linking aid across
#'    already fitted waves.
#' 6. Use [detect_anchor_drift()] for stability monitoring on separately fitted
#'    waves.
#' 7. Use [build_linking_review()] when you need one operational synthesis
#'    object rather than separate anchor/drift/chain tables.
#' 8. Run [analyze_dff()] only after checking connectivity and common-scale
#'    evidence.
#'
#' @section Keep anchor roles separate:
#' - A common-element or common-rating link is observed overlap in the design.
#' - `anchors` are direct equality constraints on selected parameter values.
#' - `group_anchors` constrain declared group means and require an externally
#'   defensible target or equal-mean assumption.
#' - [make_anchor_table()] exports candidate direct constraints mechanically;
#'   by default it requires a current inference-ready source fit, but it does
#'   not select invariant elements. `readiness_policy = "review"` is for
#'   review-only extraction, not anchor reuse.
#'
#' Constraint-based coordinate transfer does not create empirical overlap.
#' Before reuse, also verify consistent model/score/orientation/population
#' conventions, cross-run element identity, and the relevant assignment
#' connectedness. [review_mfrm_anchors()] checks syntax and receiving-data
#' support, not those substantive assumptions.
#'
#' @section Common persons, rotating assignments, and random assignments:
#' These are ways to collect ratings, not instructions to fix a parameter.
#' For example, asking all eight raters to score the same four performances
#' creates common-person links. Their abilities can still be estimated;
#' do not pass them to `anchors` unless known values really are intended.
#'
#' Existing [build_mfrm_sim_spec()] choices let you inspect these structures:
#' - Common persons: `assignment = "sparse_linked"` with
#'   `sparse_controls = list(link_persons = 4, link_raters_per_person = 8)`.
#'   Set `n_rater = 8` and specify the ordinary `raters_per_person` separately.
#'   The counts illustrate syntax, not a recommended number of linking persons.
#' - Overlapping rotation: `assignment = "rotating"` with
#'   `raters_per_person = 2` assigns successive subsets of raters.
#' - Random subsets: `assignment = "sparse_linked"` with
#'   `sparse_controls = list(link_persons = 0, assignment_mode = "random")`.
#'   This does not guarantee connectedness or equal workloads for every draw.
#' - Task-specific or externally planned assignments: `assignment = "skeleton"`
#'   accepts an explicit Person-by-facet template. The ordinary rotating
#'   generator assigns the chosen raters across all generated criteria; it
#'   is not a general implementation of task-nested spiral designs.
#'
#' Generate ratings with [simulate_mfrm_data()] and inspect them with
#' [describe_mfrm_data()]. Use the linking vignette for runnable examples.
#' A generator can represent a design without qualifying an estimator for it.
#' Compare both the amount and placement of overlap, and what happens when
#' a critical common person or rater is unavailable. Keep rating costs equal
#' when the question is which allocation uses a fixed budget better.
#'
#' @section Design the link, not only the anchor count:
#' There is no universal adequate anchor count or percentage. Linking quality
#' also depends on where overlap occurs, whether links are distributed or rely
#' on one critical bridge, the rating workload and coverage, and model--data
#' fit within the linking set. Inspect both Rater-centered and Person/task-
#' centered assignment graphs when the design is sparse. The package's
#' five-element `LinkSupportAdequate` screen is therefore a local warning, not
#' a design recommendation or proof of a common scale.
#'
#' Fixed anchor values also carry uncertainty from their source calibration.
#' Current anchor, drift, and chain summaries do not propagate all source-fit,
#' offset, or cross-fit covariance. Compare substantively defensible anchor
#' sets and report sensitivity when conclusions depend on the selected set.
#' The bounded fixed-SD RSM/PCM anchor-origin route now offers explicit
#' independent-source covariance propagation through [mfrm_facet_intervals()]
#' and [mfrm_population_intervals()]. This does not change drift or chain
#' summaries, correct noninvariant anchors, or handle dependent samples.
#'
#' @section Reading the theoretical background:
#' Uto (2021, pp. 1442--1444) explains fixed common-rater/task linking. Its
#' simulation fixes generated common parameters and removes a redundant rater
#' identification constraint; its Bayesian estimation and constraints differ
#' from mfrmr. It does not validate every centering convention above or provide
#' uncertainty propagation for anchors estimated in a previous mfrmr fit.
#'
#' Robitzsch (2024, pp. 1181--1182, 1187--1188) studies bias and linking error
#' from random uniform DIF in fixed-item 2PL calibration. It separates this
#' issue from sampling Persons and does not provide a variance estimator for
#' the fitted population parameters. It motivates checking anchor invariance,
#' not treating current mfrmr SEs as including every source of linking error.
#' For the distinction between an origin constraint and a fixed external value,
#' see the FACETS centering and element documentation cited below; its constraint
#' convention is not identical to mfrmr's.
#'
#' @section Drift and chain inference boundary:
#' [detect_anchor_drift()] and [build_equating_chain()] currently estimate one
#' pooled offset across all selected common elements and facets. Their
#' element-level SE ratios omit estimated-offset uncertainty and cross-fit
#' covariance, and chain output does not propagate uncertainty across adjacent
#' links. `Offset_SD` is residual spread, not an offset SE. Treat flags,
#' support counts, and cumulative offsets as review screens; when a common
#' shift across facet blocks is not defensible, analyze one coherent linking
#' facet at a time.
#' `ScreeningFallback` indicates that preliminary screening retained nothing
#' and the offset instead used all finite differences. `N_Contributing` counts
#' the rows actually entering the final average, which can differ from
#' `N_Retained`. Neither count certifies invariant anchors.
#' Missing offsets/spread are explained by `OffsetStatus`, `OffsetReason`
#' and, for chains, `OffsetSDStatus`. A single retained element has undefined
#' residual sample SD, not zero uncertainty. Drift `FlagBasis` distinguishes
#' the magnitude-only screen from a screen with an available SE ratio.
#'
#' @section Which helper answers which task:
#' \describe{
#'   \item{[subset_connectivity_report()]}{Summarizes connected subsets,
#'   bottleneck facets, and design-matrix coverage.}
#'   \item{[make_anchor_table()]}{Extracts candidate direct-anchor values from a
#'   fit without certifying them for reuse.}
#'   \item{[anchor_to_baseline()]}{Anchors new raw data to a baseline fit and
#'   returns anchored diagnostics plus a consistency check against the baseline
#'   scale.}
#'   \item{[detect_anchor_drift()]}{Compares fitted waves directly to flag
#'   unstable anchor elements.}
#'   \item{[build_equating_chain()]}{Accumulates screened pairwise links across
#'   a series of administrations or forms.}
#'   \item{[build_linking_review()]}{Synthesizes anchor review, drift, and
#'   screened-chain evidence into one operational review surface.}
#'   \item{[analyze_dff()]}{Describes group residual differences or linked
#'   subgroup facet estimates. Linked refit point contrasts remain screening-only
#'   because their uncertainty is conditional on baseline anchors.}
#'   \item{[mantel_haenszel_table()]}{Compares binary responses between two
#'   groups within declared matching strata, with one response per Person and
#'   target. Common odds-ratio inference does not certify matching validity
#'   or support pooled repeated ratings.}
#'   \item{[ordinal_logistic_table()]}{Compares ordered responses using nested
#'   cumulative logistic models for independent Persons. Separates uniform,
#'   nonuniform and joint group tests from a proportional-odds diagnostic;
#'   retains numerical failures and common-range effect contrasts.}
#' }
#'
#' @section Practical linking rules:
#' - Check connectedness before interpreting subgroup or wave differences.
#' - Use [mfrm_network_analysis()] for assignment/co-observation connectedness.
#'   Agreement, severity-direction, and halo networks describe score relations;
#'   they do not establish an empirical link or a common measurement scale.
#' - Residual differences do not isolate differential functioning and have no
#'   tests or binary classifications. Even with an
#'   adequate link, refit SEs condition on baseline anchors and omit their
#'   uncertainty and cross-refit covariance; formal inference remains unavailable.
#' - Always name the facet, facet level, and group pair involved in a DFF
#'   contrast. A generic "DIF exists" statement is not interpretable in a
#'   many-facet design.
#' - Linked refit contrasts remain descriptive; neither method assigns ETS
#'   A/B/C labels or supports a formal differential-functioning decision.
#' - Read a DFF figure with `plot_data(p, "pairs")`: summary plots retain
#'   all source comparisons with a `PlotStatus`, including those omitted by
#'   `top_n` or lacking a contrast. Both views retain method settings and,
#'   for residual results, the response-matching `cases` table. Positive
#'   residual and refit contrasts have different meanings; see [plot.mfrm_dif].
#' - For sparse binary cells, exact MH conditions on the stratum margins;
#'   intervals may be wide and p-values discrete. This does not correct invalid
#'   matching or dependence. Ordinal-logistic LR tests remain asymptotic:
#'   successful numerical convergence does not establish small-sample validity.
#' - MH or ordinal-logistic comparisons require explicit matching,
#'   and residual DFF's `$cases` records the actual requested-response matches.
#'   Inspect exclusions before comparing methods. Residual `N` counts rows;
#'   `NPersons` counts distinct Persons, without establishing independence.
#'   All methods require
#'   comparable rating plans and an appropriate unit of independence. A raw
#'   total across different plans or multiple rows from the same Person does
#'   not establish these conditions. See [analyze_dff()] and the companion
#'   vignette for differences in targets and response-model assumptions.
#'   Its "Compare an explicitly shared response roster" example aligns group
#'   direction, retains method-specific case counts and exclusions, and saves
#'   complete residual/MH/logistic results together. It uses a common binary
#'   item roster; it does not support pooling repeated Person-target ratings.
#' - Treat drift flags as prompts for review, not automatic evidence that an
#'   anchor must be removed.
#' - In drift/link-review results, `LinkSupportAdequate = FALSE` warns that
#'   at least one linking facet retained fewer than 5 common elements after
#'   screening. Refit DFF uses a different screen: by default at least five
#'   anchored levels in total across the non-target linking facets in each
#'   subgroup. Inspect `ScaleLinkStatus`, `LinkingThreshold`,
#'   `LinkingAnchoredLevels1`/`LinkingAnchoredLevels2`, and
#'   `LinkingDetail1`/`LinkingDetail2`. Neither screen proves invariance or
#'   supplies a universally adequate anchor count.
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
#'   fit each wave separately -> [detect_anchor_drift()] ->
#'   [build_linking_review()] -> [plot(x)][plot.mfrm_anchor_drift].
#' - Group comparison route:
#'   [subset_connectivity_report()] -> [analyze_dff()] ->
#'   [dif_report()] -> [plot(x, type = "heatmap")][plot.mfrm_dif] / [plot(x, type = "summary")][plot.mfrm_dif].
#'
#' @section Companion guides:
#' - For visual follow-up, see [mfrmr_visual_diagnostics].
#' - For report/table selection, see [mfrmr_reports_and_tables].
#' - For end-to-end routes, see [mfrmr_workflow_methods].
#' - For a longer walkthrough, see
#'   `vignette("mfrmr-linking-and-dff", package = "mfrmr")`.
#'
#' @references
#' TAM documentation: \href{https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html}{MML estimation and fixed parameters}
#' and \href{https://alexanderrobitzsch.github.io/TAM/reference/tam.jml.html}{JML and fixed Person parameters}.
#'
#' ACER ConQuest manual: \href{https://conquestmanual.acer.org/s4-00.html}{Command reference},
#' especially `import` (anchor parameter types) and `set` (location and scale
#' constraints). Defaults and supported options depend on the program version.
#'
#' FACETS documentation: \href{https://www.winsteps.com/facetman/centerfacet.htm}{Centering facets}
#' and \href{https://www.winsteps.com/facetman/elements.htm}{Elements and anchors}.
#'
#' Linacre, J. M. (2012). \href{https://www.winsteps.com/a/ftutorial4.pdf}{Facets Tutorial 4: Anchoring},
#' p. 9, discusses rater drift and the assumption involved in fixing a group
#' mean across sessions. It does not establish the common-recording workflow
#' described here or its uncertainty performance.
#'
#' Eckes, T. (2023). *Introduction to Many-Facet Rasch Measurement*,
#' Section 9.1, pp. 151--155, Table 9.1. Peter Lang.
#' ISBN 978-3-631-90305-6.
#'
#' Myford, C. M., & Wolfe, E. W. (2000). Strengthening the ties that bind:
#' Improving the linking network in sparsely connected rating designs.
#' *ETS Research Report Series*, 2000(1).
#' \doi{10.1002/j.2333-8504.2000.tb01832.x}
#'
#' Uto, M. (2021). Accuracy of performance-test linking based on a many-facet
#' Rasch model. *Behavior Research Methods*, 53(4), 1440--1454.
#' \doi{10.3758/s13428-020-01498-x}
#'
#' Wind, S. A., & Jones, E. (2018). The stabilizing influences of linking set
#' size and model--data fit in sparse rater-mediated assessment networks.
#' *Educational and Psychological Measurement*, 78(4), 679--707.
#' \doi{10.1177/0013164417703733}
#'
#' Robitzsch, A. (2024). Bias and linking error in fixed item parameter
#' calibration. *AppliedMath*, 4(3), 1181--1191.
#' \doi{10.3390/appliedmath4030063}
#'
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("mfrmr_linking_and_dff", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_core")
#' people <- unique(ratings$Person)
#' baseline_data <- ratings[ratings$Person %in% people[1:24], , drop = FALSE]
#' new_data <- ratings[ratings$Person %in% people[25:48], , drop = FALSE]
#' baseline <- fit_mfrm(
#'   baseline_data, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", rating_min = 1, rating_max = 4,
#'   category_policy = "preserve"
#' )
#' anchors <- make_anchor_table(baseline, facets = "Criterion")
#' review <- review_mfrm_anchors(
#'   new_data, "Person", c("Rater", "Criterion"), "Score", anchors = anchors,
#'   rating_min = 1, rating_max = 4, category_policy = "preserve"
#' )
#' review$issue_counts
#' linked <- anchor_to_baseline(
#'   new_data, baseline, "Person", c("Rater", "Criterion"), "Score",
#'   anchor_facets = "Criterion", anchor_policy = "error",
#'   rating_min = 1, rating_max = 4, category_policy = "preserve"
#' )
#' linked$fit$facets$others
#' # Fixed-anchor differences are zero by construction, not evidence of stability.
#' linked$drift[, c("Facet", "Level", "Drift")]
#'
#' }
#' @name mfrmr_linking_and_dff
NULL
