#' Confidence intervals for dependability and SEM in D-study plans
#'
#' Quantify uncertainty in relative G, absolute dependability Phi or measurement SEM for
#' specified future counts of common tasks, raters or other random facets.
#' Most methods use a complete balanced ANOVA source. An explicitly selected
#' experimental REML profile also accepts identifiable incomplete MINQUE sources
#' for G/Phi and SEM; future plans still describe complete balanced ratings.
#' For one-facet relative SEM, residual contrasts also accept incomplete sources
#' without likelihood fitting, under the normal residual model.
#' For one-facet G, explicit spectral F intervals accept incomplete sources
#' under independent normal Person effects and residuals.
#' Explicit delta-method coefficient SEs and normal-reference intervals also
#' use the original ANOVA/MINQUE estimator for supported crossed/nested sources.
#'
#' @param object A result from [mfrm_d_study()] based on a
#'   balanced ANOVA G-study with one or two crossed random measurement facets,
#'   or the supported balanced nested design. Nested G/Phi defaults require
#'   unchanged source child counts per parent; changed counts require an
#'   explicit experimental method. See "Nested plans" below. For G/Phi and SEM,
#'   `method = "reml_profile"` also accepts an incomplete or unequal MINQUE
#'   source with its original used ratings retained in the object.
#'   `method = "residual_contrast"` accepts retained one-facet ANOVA/MINQUE
#'   ratings for `parm = "RelativeSEM"` only, including plan differences.
#'   `method = "spectral_f"` accepts these sources for `parm = "G"`,
#'   including plan differences, with positive residual degrees of freedom.
#'   `method = "delta"` accepts ANOVA or identifiable MINQUE sources for
#'   individual-plan G/Phi, including changed nested child counts. MINQUE
#'   requires retained assignments and moment inputs for its sampling covariance.
#' @param parm One of `"G"` (default), `"Phi"`, `"RelativeSEM"` or
#'   `"AbsoluteSEM"`. With `reference`, return differences in that metric
#'   for the selected source and supported method.
#' @param level Confidence level strictly between zero and one.
#' @param score One original score name.
#' @param composite One named composite, or `"Composite"` for vector weights.
#'   Select either `score` or `composite`. With neither specified, select the
#'   sole composite, or the first original score if no composite is present.
#'   Multiple composites require a selection. Weights must be prespecified.
#' @param assumption Required: `"normal"` asserts independent normal random
#'   effects with constant component variances. The source persons and facet
#'   levels are sampled from their intended populations. This is not tested.
#' @param draws Number of independent pivotal draws for one-facet Phi GCI,
#'   default 100000.
#'   Must be an integer at least 1000, with at least 100 expected draws in
#'   each tail at the requested `level`. Increase it when the returned endpoint
#'   Monte Carlo bands are too wide for the intended decision. Unused by other methods.
#' @param seed Nonnegative integer seed for one-facet Phi GCI, default 1. The
#'   caller's random seed is preserved. The current R RNG kind is used and
#'   recorded. Unused by other methods.
#' @param reference Optional row number in `object$design_grid`. Return each
#'   other plan's selected metric minus this reference plan's value.
#'   The default supports one-facet sources and nested plans with a common
#'   future child count (unchanged from the source for G/Phi). Map the whole F set for G or generalized
#'   interval for Phi, preserving their different coverage qualifications.
#'   Explicit `method = "reml_profile"` also supports crossed/nested plans
#'   changing both facet counts, including incomplete/unequal sources.
#'   One-facet and common-child nested plans map a shared profile interval;
#'   other comparisons constrain the difference directly, using optional nloptr.
#'   Its calibration remains experimental, especially at zero differences.
#'   For general two-facet crossed/nested comparisons, explicitly select
#'   `method = "mean_square_region"` for a conservative joint-region enclosure
#'   for G, Phi or SEM. Approximate-F differences are not supplied.
#'   At least two plans are required. `NULL` returns individual-plan intervals.
#' @param method `"auto"` (default) selects `"exact_f"` for one-facet G and
#'   `"gci"` for one-facet Phi, including the nested reduction below.
#'   For individual G/Phi with two crossed random facets, it selects `"mls"`,
#'   an approximate signed-MLS null inversion with conservative completion
#'   of undefined widths; see "Two-crossed-facet MLS intervals" below.
#'   To compare other choices, or for nested plans with changed child counts, select
#'   `"approximate_f"`, `"mean_square_region"`, or a vector containing both
#'   to compare methods on the same data and target. Neither is an automatic
#'   default. `"mls"` can also be compared with these methods for crossed G/Phi.
#'   `"mean_square_region"` requires the optional lpSolve package.
#'   Balanced two-crossed-facet G also accepts explicit `"score"`, an
#'   approximate score inversion with extended nuisance parameters. It uses
#'   stored mean squares, preserves original points and requires no optional
#'   optimizer package. It can be compared with the other supported G methods.
#'   It is not a Phi, SEM, nested, incomplete-source or `reference` method.
#'   For SEM, `"auto"` selects `"mls"`: a single contributing mean square
#'   uses its chi-square interval; positive sums use Graybill-Wang MLS and
#'   signed sums use Ting MLS. This is a specified calculation rule, not
#'   selection based on which interval looks favorable. `Method` records
#'   the calculation used for each plan.
#'   Two-facet crossed/nested SEM also accepts `"mean_square_region"`, using
#'   linear error-variance projections followed by square roots. Select one
#'   method at a time unless including `"reml_profile"`, or comparing supported
#'   experimental coefficient methods above. `"reml_profile"` is an explicit G/Phi or SEM choice for
#'   balanced or incomplete sources. It requires optional lme4 and minqa,
#'   fits a separate normal REML model, and retains the original point estimate.
#'   It can be compared with other methods for the same supported target in a vector;
#'   `reference` uses a shared interval image or a direct difference profile
#'   as described below. The latter additionally requires nloptr. Its chi-square
#'   calibration is approximate, especially at component boundaries. It is
#'   never selected automatically. No fallback follows a failed calculation.
#'   Explicit `"residual_contrast"` supplies one-facet RelativeSEM intervals
#'   and differences from retained ANOVA/MINQUE ratings without likelihood
#'   fitting or simulation. It maps one residual chi-square interval, preserves
#'   original points and reports the separate projected `IntervalCenter`.
#'   It can be compared with MLS or REML for this supported target. It is not
#'   an AbsoluteSEM, G/Phi or multifacet method; see *Residual-contrast SEM*.
#'   Explicit `"spectral_f"` supplies one-facet G intervals and paired
#'   differences from retained ANOVA/MINQUE ratings. It removes facet effects
#'   and inverts an exact normal-model F pivot without likelihood fitting or
#'   simulation. It preserves original points and can be compared with
#'   `"exact_f"` (balanced sources) or `"reml_profile"`. It is not an
#'   automatically selected, Phi, SEM, nested or multifacet method.
#'   Explicit `"delta"` supplies G/Phi sampling SEs and normal-reference
#'   intervals from the original component estimator and its full sampling
#'   covariance. It accepts supported one-facet, crossed and nested sources,
#'   including incomplete MINQUE, and can be included in a method vector.
#'   It requires no refitting or simulation. It is not a SEM or `reference`
#'   method. Defaults are unchanged; see *Delta coefficient SEs and intervals*.
#' @param ... Reserved; additional arguments are rejected.
#'
#' @details
#' For a source-by-target choice table, see *Choose the G-theory target and
#' uncertainty* in [mfrmr_methods_reference]. Incomplete-source calculations
#' are available through explicit methods; this does not qualify a generally
#' calibrated automatic default. Read the selected method's assumptions below.
#'
#' Each interval concerns a population coefficient or SEM of a future plan: uncertainty comes
#' from the source G-study, not from collecting the proposed future ratings.
#' GENOVA-adjusted D-studies can use the same population-parameter intervals.
#' `Estimate` and `PointStatus` describe the adjusted point; `IntervalEstimate`
#' and `IntervalPointStatus` retain the original source-study projection.
#' `ComponentTreatment` and `IntervalBasis` distinguish their bases. The
#' original source mean squares, data, bounds, diagnostics and unavailable
#' states are preserved; adjusted components never replace the source inputs.
#' The bounds are not recentered and need not contain the adjusted point.
#' An existing method-specific `IntervalCenter` remains separate from both.
#' Where a method supplies a coefficient sampling SE, it is named `IntervalSE`
#' for adjusted results. `SamplingVariance`, error-variance columns and other
#' diagnostics also retain their original method basis. None describes the
#' sampling distribution of the zero-substituted or EMS-recalculated estimator.
#' In the method descriptions below, references to the original `Estimate`
#' mean `IntervalEstimate` when component treatment is adjusted.
#'
#' This distinction follows from the common population target: changing a
#' reported point does not change the confidence-set event based on the source
#' study. Exact procedures keep their stated normal-model qualifications;
#' approximate or experimental methods do not become exact through adjustment.
#' GENOVA's unchanged component SEs cannot be interpreted as an adjusted
#' estimator's SE. Its truncation is nonsmooth at zero; an adjusted-estimator
#' delta interval would require additional derivation and boundary validation.
#' Saved intervals must match the D-study's treatment as well as its source.
#' Use `confint(d, reference = ...)` for supported population plan differences;
#' the separate [mfrm_multivariate_d_compare()] delta-estimator API still
#' requires raw components.
#'
#' For a balanced two-crossed-facet source, start with
#' `confint(d, parm = "Phi", assumption = "normal")` for absolute dependability,
#' or `parm = "G"` for relative ordering. Read `Estimate`, `Lower`, `Upper`,
#' `Status` and `FormulaComplete`. The default is an approximation: a returned
#' interval is not evidence that the model assumptions hold. Use
#' `plot(d, metric = "Phi", intervals = TRUE, assumption = "normal")` to add
#' its ribbon, and `attr(ci, "notes")` to read cautions saved with a result `ci`.
#'
#' @section Two-crossed-facet MLS intervals:
#' Under the balanced independent normal-effects model, mean squares `M[j]`
#' are independent `lambda[j] * chi-square(df[j])/df[j]`. Let `theta=A*lambda`
#' be the variance components, `U=theta[Person]`, and `E=error' theta` the
#' relative or absolute error variance for the prespecified future counts.
#' A candidate `rho=U/(U+E)` is equivalent to the linear null
#' `h(rho)'lambda=0`, with `h(rho)'=((1-rho)e_P' - rho*error')*A`.
#' The denominator must be positive; positive source mean-square expectations
#' and residual variance are assumed. No estimated ratio or component is clipped.
#'
#' Compute the unrestricted Ting MLS interval for this signed linear combination
#' and retain `rho` when it contains zero. The asymmetric chi-square/F corrections
#' use source degrees of freedom. This is mfrmr's inversion of the published
#' linear-combination construction; Ting et al. and Hess/Iyer do not directly
#' establish coverage for this G/Phi procedure. The inversion inherits approximate,
#' rather than exact, calibration. The checked two-term cases at the usual
#' 95% level reproduce F inversion. This is not an identity for arbitrarily
#' low confidence levels: squaring a reversed chi-square correction can lose
#' its direction. The sampling evaluations below concern 95% intervals.
#'
#' Between coefficient sign changes, the squared widths and squared null point
#' are quadratic in `rho`. All polynomial crossings in `[0,1]` are evaluated,
#' including sign-change points, where the MLS formula can be discontinuous.
#' Diagnostics retain accepted pieces and individual nodes; `Lower`/`Upper`
#' enclose their union. This avoids a finite grid or repeated model fitting.
#' The polynomials use shifted and scaled coordinates to resolve narrow sets
#' near a boundary. If root separation cannot be represented on the coefficient
#' scale, bounds are missing with an explicit `Status` and
#' `NumericalComplete=FALSE`; the point estimate is retained. This is distinct
#' from a genuine empty set or a formula-defined endpoint at zero or one.
#'
#' A squared MLS width can be negative for some candidate values, especially
#' with few degrees of freedom and highly unequal components. Such candidates
#' are conservatively retained rather than assigned a zero width or silently
#' dropped. `FormulaComplete=FALSE` identifies this completion, and each row's
#' `attr(ci, "diagnostics")` entry retains the `undefined` ranges. This can only
#' enlarge the retained set; it cannot make the approximate calibration exact.
#' A numerical inversion failure remains unavailable and is distinct from this
#' defined completion rule. An empty set also remains explicit. An unavailable
#' original ANOVA point does not prevent calculation of an interval.
#'
#' Normal-model checks illustrate both potential and limits. In a matched
#' 400-source study with 60 persons, 4 raters, 12 tasks and component variances
#' `(1,.2,.3,.3,.4,.1,1)`, the source-plan G interval had mean width .0945
#' and coverage .950; Phi had mean width .3287 and coverage .970.
#' The conservative mean-square region gave widths .1974/.7249 and coverage
#' 1/1 in those same samples. The component order is Person, Rater, Task,
#' Person:Rater, Person:Task, Rater:Task, residual.
#' A further eight-condition normal study used 2,000 sources per condition,
#' Person counts 5--80, source facet counts from 2-by-2 to 2-by-20, zero/near-zero
#' Person variance, dominant facets/interactions and changed future counts.
#' All-attempt inclusion was .9505--.9885 for G and .948--1 for Phi after the
#' stated completion. At .95, Monte Carlo SE is about .0049 with 2,000 sources.
#' Some small-source Phi conditions returned `[0,1]` in over 90% of available
#' cases. These are limited normal-model evaluations, not uniform coverage,
#' robustness, or a promise of useful precision in every design. Increasing
#' future counts does not compensate for a weak source design.
#'
#' This choice supplies pointwise individual-plan intervals for balanced two-
#' crossed-facet sources and prespecified scores/composites. It does not supply
#' simultaneous curve bands, plan-difference intervals, or incomplete/nested-
#' source inference. The exact one-facet G and generalized one-facet Phi
#' defaults retain their different constructions. The name `"mls"` for SEM
#' instead refers to a direct error-variance interval followed by a square root.
#'
#' @section Interpreting the target and other defaults:
#' More future tasks can improve G but do not increase the source degrees of
#' freedom. A narrow interval is not guaranteed by a large proposed task count.
#' G concerns relative ordering; these intervals do not give pass/fail accuracy,
#' or Person-rank intervals. Phi includes the random measurement facet's main
#' effect in absolute error; use it when absolute score levels matter.
#'
#' The F pivot is exact under the balanced independent normal-effects model,
#' with positive residual variance. It is not a distribution-robust interval
#' for arbitrary numeric ordinal scores. It applies to a selected original
#' score or a fixed linear composite, including prespecified difference scores.
#' General two-facet plans and incomplete-source MINQUE(0) do not have this
#' pivot. The nested unchanged-child-count subcase below does reduce to it.
#' Balanced F/GCI, MLS, score and mean-square-region methods use stored ANOVA
#' summaries. The score method fits null nuisance expectations to those
#' summaries. One-facet Phi GCI draws independent chi-square quantities; it
#' is not an exact frequentist-coverage guarantee. Spectral F and residual
#' contrasts require retained ratings; delta uses ANOVA summaries or retained
#' MINQUE assignments and moment inputs. REML profiles require used ratings
#' and refit a separate model. These input requirements also apply to saved
#' G-study objects; reopening already computed intervals does not refit them.
#'
#' @section Score intervals for two-crossed-facet G:
#' Explicit `method = "score"` uses the independent normal-theory mean-square
#' sampling law of a balanced Person-by-two-random-facet source. For proposed
#' G = \eqn{\rho}, write the population target as
#' \eqn{\theta_P/(\theta_P+e^T\theta)} and test
#' \eqn{h(\rho)^T\lambda=0}, where \eqn{\theta=A\lambda},
#' \eqn{h(\rho)^T=((1-\rho)e_P^T-\rho e^T)A}, and \eqn{\lambda}
#' contains the positive mean-square expectations. Nuisance components may be
#' negative as long as this restricted marginal covariance remains positive
#' definite. These computational values do not replace the ANOVA components.
#' The nuisance score is checked independently of optimizer status labels.
#' If \eqn{S_\lambda=\kappa h}, the expected-information statistic is
#' \eqn{T=\kappa^2\sum_j 2h_j^2\lambda_j^2/df_j}. Invert its approximate
#' chi-square(1) test over \eqn{0\leq\rho\leq1}.
#'
#' This construction is motivated by Shedden and Ekvall's extended-nuisance
#' score approach; their preprint does not establish G/Phi ratio intervals.
#' For the active independent strata, the spectral/Frobenius bound in their
#' asymptotic condition is \eqn{1/\sqrt{\min_j df_j}}. G uses Person-containing
#' strata, whose degrees of freedom increase with source Persons even with
#' fixed facet counts. This does not give exact finite-sample coverage.
#' A retained normal-model evaluation used 360 independent sources per
#' condition, six conditions and three future plans: true-value inclusion was
#' 93.1--97.8%, with Monte Carlo SEs of 0.8--1.3 percentage points. These checks
#' test the statistic at the truth, not complete endpoint coverage. The lowest
#' cell's exact binomial 95% interval was 89.9--95.5%. The method remains
#' experimental and is never selected automatically.
#'
#' The numerical inversion compares two finite grids, includes contrast
#' sign-change locations, and verifies cutoff crossings. Diagnostics retain
#' all scanned fits and accepted pieces; bounds are their envelope. These
#' checks do not prove global nuisance optima or exclude narrower undetected
#' pieces. `NumericalComplete` describes the completed numerical check.
#' `Empty = TRUE` means the scanned confidence set is empty: it is not
#' an estimate of zero or zero uncertainty. Empty sets and numerical failures
#' both have missing `Lower`/`Upper`, but distinct status and diagnostics.
#' `CandidateLower`/`CandidateUpper` retain exploratory bounds after failure.
#' Cutoff-crossing errors must be smaller than
#' `min(1e-5, 1e-5 * critical)`. If this tolerance falls below
#' `64 * .Machine$double.eps * max(1, critical)`, the cutoff is numerically
#' unresolved: bounds are missing, `NumericalComplete = FALSE` and
#' `Empty = FALSE`. Original points remain available; this is not an empty
#' confidence set or a minimum sample-size restriction. Diagnostics retain
#' `critical`, `cutoff_tolerance`, `numeric_resolution`, `root_tolerance`
#' and `grid_tolerance`. Root and grid-comparison tolerances tighten for
#' small critical values. These arithmetic checks do not certify global
#' nuisance roots or complete confidence-set topology.
#' Positive observed mean squares are required. Few source Persons can yield
#' broad intervals or inaccurate calibration. The method is pointwise for a
#' prespecified score/composite and future plan, with no distribution-robust
#' or simultaneous-coverage claim.
#'
#' Phi also depends on few-df facet strata. Extending this same score rule to
#' Phi produced very broad sets and changes between nuisance roots near zero;
#' it is therefore not offered for Phi. Separate score intervals cannot be
#' subtracted to obtain a paired-plan difference interval.
#'
#' @section REML profiles for balanced or incomplete sources:
#' For a fixed score or composite, this method fits the normal random-effects
#' model of the declared grouping assignment: one common mean and independent
#' Person, facet, interaction and residual effects. Source omissions must be
#' noninformative conditional on that fixed assignment. It does not model
#' outcome-dependent missingness, unequal future exposure or selection of
#' score weights from the observed data. The residual variance is positive.
#' Raw used ratings are needed; ANOVA/MINQUE component summaries alone are
#' insufficient for the incomplete-source likelihood.
#'
#' Write the variance ratios as `v`, with residual ratio one. For each plan,
#' let `c` contain its relative or absolute error divisors, with `c[Person] = 0`.
#' The coefficient is `rho = v[Person] / (v[Person] + sum(c * v))`.
#' Fixing a proposed `rho` eliminates the Person ratio as
#' `rho / (1 - rho) * sum(c * v)`. The remaining nonnegative ratios are
#' optimized; the common mean and residual scale are profiled analytically
#' using the restricted Gaussian likelihood. The accepted numerical set has
#' deviance increase no larger than `qchisq(level, df = 1)`.
#'
#' For SEM, constrain the error variance `tau = sum(c * theta)` directly,
#' where `theta = s2 * v` includes the positive residual variance `s2`.
#' The profiled residual scale at fixed ratios is `s2_hat(v)`.
#' Restore it to the criterion using
#' `D(v,s2) = D_profile(v) + (n-1)*(log(s2/s2_hat) + s2_hat/s2 - 1)`.
#' At fixed positive `tau`, substitute `s2 = tau / sum(c*v)` and optimize
#' the remaining ratios. Take square roots of the error-variance endpoints
#' and restore the original score units. This is a constrained likelihood
#' for the error variance, not addition of component bounds or a delta method.
#' The residual scale cannot be profiled freely while fixing an absolute
#' variance target. SEM intervals describe sampling uncertainty in the
#' population measurement SEM, not intervals for individual Person scores.
#'
#' This is an implementation of target-specific restricted-likelihood
#' inversion, not a MINQUE likelihood or a published exact G-theory pivot.
#' Bates et al. describe the REML criterion and sparse computation, not a
#' validation of these G/Phi or SEM intervals. It also differs from lme4's ordinary
#' parameter-profile intervals, which use an ML criterion. The chi-square
#' cutoff is asymptotic; a generic boundary mixture is not substituted.
#' The special-case results of Crainiceanu and Ruppert illustrate why boundary
#' calibration cannot be assumed for arbitrary multi-component designs.
#'
#' A limited normal-model feasibility check used 40 samples per assignment
#' and population condition. True-value likelihood-ratio acceptance ranged
#' from 92.5% to 100% with positive components and from 87.5% to 100% in the
#' checked zero-component populations, counting numerical failures as misses.
#' These imprecise checks do not establish 95% coverage or general robustness,
#' and true-value acceptance is not a full endpoint-performance study.
#' Numerical endpoint calculations were separately checked on selected data.
#' Use this explicit experimental option only with its qualifications;
#' agreement with another interval in one dataset does not establish coverage.
#'
#' A separate SEM feasibility check reused 40 samples per assignment and
#' population. Positive-component true-value acceptance ranged from 92.5%
#' to 100%; checked boundary populations ranged from 82.5% to 100%.
#' The lowest result was absolute SEM in an incomplete crossed source with
#' 30 persons, 4 raters and 5 tasks and zero facet/interaction variances;
#' numerical availability there was 97.5% and the inclusion Monte Carlo SE
#' was 6.0 percentage points. These are coarse warning signals, not evidence
#' that nominal 95% inference is reliable. They concern acceptance at the
#' true target, not coverage of all returned interval endpoints.
#'
#' Three starting points, an independent BOBYQA optimization, and finite
#' differences with boundary checks verify numerical solutions. For G/Phi, a grid
#' including the REML estimate locates cutoff crossings; root finding refines
#' them. Disconnected accepted pieces are retained in `diagnostics`, while
#' `Lower/Upper` give their envelope. If the REML coefficient center rounds to
#' one, the current coefficient-scale search is unavailable: bounds remain
#' missing with an explicit numeric-resolution `Status`. The source
#' optimization may still have passed
#' its checks; the rounded center does not establish perfect reliability or
#' zero uncertainty, and the uncomputed interval need not have zero width.
#' This interval search also does not certify a global optimum or exclude
#' narrow pieces missed between grid points. An unresolved
#' profile or upper tail is flagged, with exploratory `CandidateLower/Upper`
#' retained separately; it is not replaced by a zero-width or `[0,1]` interval.
#' For a very small requested confidence level, a likelihood-ratio cutoff
#' indistinguishable from numerical error leaves the interval unavailable.
#' Endpoint tolerance is `min(1e-5, 1e-5 * cutoff)` and must exceed the
#' floating-point resolution check `64 * .Machine$double.eps *
#' max(1, abs(source_deviance))`. Row-aligned `diagnostics[[i]]$cutoff`
#' retains these values. This check does not certify optimizer accuracy.
#' A positive interval width or nonzero endpoint lost on restoring score units
#' also gives missing bounds with an explicit `Status`, keeping the original
#' point and representable REML center. These checks apply to component,
#' covariance, G/Phi, SEM and direct-difference profiles; an identical-plan
#' direct contrast remains an exact zero when its source checks pass.
#' Separate score and weight scaling preserves a composite under reciprocal
#' changes of units.
#'
#' For SEM, `log1p` variance-ratio coordinates retain exact zero components
#' while reducing scale differences. Search outward from the REML center
#' until the first checked cutoff crossing in each direction, then refine
#' the roots and check their likelihood-ratio residuals. This returns a
#' numerical interval around the center; it does not enumerate disconnected
#' remote accepted regions. A finite scan cannot prove global nuisance
#' optima or rule out narrower undetected crossings. Unresolved required
#' evaluations, roots or tails retain unavailable bounds and separate candidates;
#' reaching a search limit does not justify an infinite endpoint.
#' Proportional error-weight vectors reuse the same profile by scaling the
#' error-variance targets. This is an exact transformation, without new fits.
#'
#' `Estimate` and `PointStatus` always retain the source ANOVA/MINQUE result.
#' `IntervalCenter` is the REML estimate of the same population coefficient or SEM.
#' A source point can be unavailable or outside the profile interval; fitting
#' nonnegative REML variances does not repair the original component matrices.
#' Forests and D-study curves retain the original source points. Save the
#' interval object to avoid fitting again when plotting or reopening it.
#' Individual-plan intervals are pointwise. The shared-profile differences
#' below have a separate common-base justification. The normal assumption and a successful optimizer do not imply
#' accurate finite-sample calibration.
#'
#' @section Plan differences from one shared REML interval:
#' With `method = "reml_profile"` and `reference`, compare one-facet plans
#' or nested plans with the same future child count per parent. The child
#' count may differ from the source. The observed assignment, including
#' incomplete cells and unequal child counts, stays in the full REML model;
#' it is not averaged or replaced by balanced data.
#'
#' These future plans have error variance `e0 / m`, where `m` counts the
#' measurement-facet levels or nested parents. For G/Phi, let `t = U / e0`,
#' where `U` is Person variance. The coefficient at count `m` is
#' `m*t / (1+m*t)`. Calculate one base profile at count one, then map its
#' whole interval through the comparison-minus-reference function. With
#' counts `b` and `a`, its stationary point is `t = 1/sqrt(a*b)`;
#' this must be included when inside the base interval. Subtracting marginal
#' confidence endpoints or transforming only the two base endpoints is wrong.
#' For SEM, multiply the shared base interval for `sqrt(e0)` by
#' `1/sqrt(b) - 1/sqrt(a)`, reversing limits when this factor is negative.
#'
#' Whenever the base interval contains its population target, all returned
#' images contain their true differences for this fixed score/composite.
#' This event transfers the base method's coverage; it does not establish
#' exact nominal coverage or remove numerical-search limitations. Mapping
#' an envelope can be conservative, and SEM still uses the interval around
#' its REML center. No additional model is fitted for each comparison.
#' Original ANOVA/MINQUE differences remain in `Estimate`, with the REML
#' differences in `IntervalCenter`. An unavailable base remains unavailable,
#' even for identical plans; candidate images are exploratory only.
#'
#' `attr(result, "base_interval")` retains the complete count-one interval
#' and its profile diagnostics. Plot the comparison with `plot(result)`;
#' difference intervals cannot be used as individual-plan curve ribbons.
#'
#' @section Direct REML profiles for other plan differences:
#' For crossed plans or nested plans with different future child counts,
#' the same call uses a direct difference constraint on the shared source
#' likelihood. It does not subtract marginal confidence endpoints.
#' Let `a` and `b` be reference and comparison error weights, `v` the variance
#' ratios with residual ratio one, `ea = sum(a*v)`, `eb = sum(b*v)` and
#' `U = v[1]`. The G/Phi difference is
#' `U*(ea-eb)/((U+ea)*(U+eb))`. A nonzero target can have two positive
#' solutions for `U`, on either side of `sqrt(ea*eb)`; both seed the search.
#' At zero difference, fit both the zero-Person and equal-error branches.
#' The equal-error solution also seeds nearby targets to avoid following
#' only the zero-Person branch. Multiple starts and independent constraint
#' and stationarity checks still cannot certify global optima.
#'
#' The SEM difference is `sqrt(s2)*(sqrt(eb)-sqrt(ea))`, in standardized
#' score units. The constraint therefore restores the residual scale in the
#' full REML criterion. At zero SEM difference the errors must be equal;
#' zero Person variance alone does not satisfy this null. Final SEM bounds
#' are multiplied by the source score scale, retaining signed score units.
#'
#' A grid and checked cutoff crossings retain detected accepted pieces;
#' `Lower/Upper` enclose those pieces, not a guaranteed enumeration of the
#' global set. `attr(result, "diagnostics")` retains each comparison's
#' `sets`, `scan`, `profiles`, `roots`, numerical failures and candidate bounds.
#' SEM diagnostic targets and sets are standardized; multiply by
#' `attr(result, "reml_source")$score_scale` for the original score units.
#' Required failed evaluations or unresolved tails leave accepted bounds
#' unavailable. A search limit is not an infinite confidence bound.
#'
#' These comparisons are **pointwise and experimental**, not simultaneous.
#' Their chi-square(1) calibration and finite-sample coverage are not
#' established generally. In particular, the two branches of the zero
#' coefficient null intersect at a singularity. Drton (2009) explains why
#' singularities and boundaries can invalidate ordinary likelihood-ratio
#' calibration; that paper is not a derivation of the present G-theory method.
#' Neither numerical completion nor a warning supplies a coverage guarantee.
#' Prefer a method with justified coverage for the design and inferential
#' purpose; use this explicit choice to examine the fitted likelihood while
#' retaining its qualifications. It is never selected automatically.
#'
#' @section Sampling uncertainty in measurement SEM:
#' SEM describes measurement error in score units; its confidence interval
#' describes how uncertain that population SEM is given the source study.
#' It is not a confidence interval for a Person's score, a sampling SE of G,
#' or an interval obtained by adding component confidence bounds.
#'
#' Write the raw component estimates as `A %*% MS`. If `w` contains the
#' future-count divisors for relative or absolute error, the error-variance
#' estimate is `w' A MS`. Apply the normal-model interval to this one linear
#' combination of independent mean squares, intersect it with nonnegative
#' variances, then take square roots. This monotone transformation preserves
#' the original interval's coverage event. It does not improve an approximate
#' MLS interval's coverage. Component formulas and nested child-count units
#' are given in [mfrm_d_study()].
#' With `method = "mean_square_region"`, project the common simultaneous
#' mean-square region onto each error variance, then take square roots. This
#' retains conservative normal-model coverage, including joint coverage over
#' future plans for the fixed score. It can be much wider than MLS and has
#' no MLS sampling-SE columns; `diagnostics` retains variance-scale bounds.
#'
#' The single-mean-square interval is exact under the stated normal model
#' when its expectation is positive. For example, one-facet relative error
#' is `MS_residual / future_count`. A nested relative error also reduces to
#' one mean square when the future number of children per parent equals
#' that in the source. Otherwise positive combinations use Graybill-Wang
#' and signed combinations use Ting MLS, as set out by Hess and Iyer (2002).
#' Betebenner (1998) applies these constructions to G-theory error variances.
#' Both are approximations, not distribution-robust or simultaneous intervals.
#' Their squared chi-square corrections can lose their usual direction at
#' very low confidence levels; the sampling evaluation below concerns 95% limits.
#' If a positive variance interval or its square-root image loses its width
#' through numerical rounding, the SEM limits are missing and `Status` explains
#' the loss of resolution. The original SEM and any calculable error-variance SE
#' are retained. This also applies to residual-contrast SEM. An unavailable
#' common source interval does not become a zero-width plan-difference interval,
#' even for identical counts. Recompute `confint(d, ...)` from a saved D-study
#' to apply current numerical checks; a G-study refit is not required.
#'
#' A bounded check used nine normal-model source conditions with 2,000
#' replications each: one-facet, crossed and nested small/larger designs,
#' including zero components. Each had three future plans and both SEMs.
#' Across the 51 positive-error targets, nominal 95% coverage ranged from
#' 93.45% to 99.00% (Monte Carlo SE 0.22--0.55 percentage points), and
#' availability from 99.75% to 100%. The lowest coverage concerned absolute
#' SEM with 10 persons, 2 raters and 3 tasks, projected to one rater and one
#' task. The all-zero relative-error cases returned unavailable intervals.
#' This limited check identifies undercoverage and conservatism; it is not
#' a guarantee for other populations, designs or nonnormal scores.
#'
#' For a one-facet source, every SEM is `s / sqrt(m)` for the same
#' single-rating error SD `s` and future count `m`. With reference count `a`
#' and comparison `b`, transform the whole interval for `s` by
#' `1/sqrt(b) - 1/sqrt(a)`, reversing endpoints when negative. This retains
#' dependence between plans. The common base coverage event protects all
#' returned differences for that fixed score; relative SEM has exact
#' normal-model coverage for positive residual variance, while absolute SEM
#' retains approximate MLS coverage. Lower SEM differences favor the
#' comparison plan. Identical counts give `[0,0]` only when the base interval
#' is available. Nested comparisons with a common future child count use the
#' same mapping in the parent count; their base may require Ting MLS.
#' [mfrm_multivariate_d_compare()] retains its separate paired-delta scope.
#'
#' @section Spectral F intervals for one-facet G:
#' Use `method = "spectral_f"` when a one-facet source is incomplete or
#' unequal and the independent normal-effects model is appropriate. The
#' future plan must still be complete and balanced. Unlike residual-only SEM
#' inference, this method requires normal Person effects as well as normal
#' residuals. Assignment and omissions must preserve the assumed model.
#'
#' Write the scalar or fixed-composite response as
#' `y = Z_P p + Z_F f + e`, with Person variance `sigma_P^2` and positive
#' residual variance `sigma_E^2`. The facet incidence contains the intercept.
#' Removing its column space also removes every realized facet effect without
#' changing the population G target. Let `U diag(sqrt(lambda)) V'` be the
#' thin singular decomposition of the projected Person incidence, keeping
#' its `r` positive singular values. Then `z = U' y` has independent normal
#' coordinates with variances `sigma_E^2 + lambda_i * sigma_P^2`.
#' Projecting off both Person and facet incidence supplies independent
#' `S_E / sigma_E^2 ~ chi-square(nu)`, where
#' `nu = n - rank([Z_P, Z_F])`.
#'
#' For `t = sigma_P^2 / sigma_E^2`,
#' `F(t) = {sum(z_i^2 / (1 + lambda_i*t)) / r} / (S_E / nu)`
#' has an F distribution with `r, nu` degrees of freedom at the true ratio.
#' This is a specialization of the spectral pivot family in E, Hannig and
#' Iyer (2008), equations 3--7, using all positive Person eigenvalues in the
#' numerator and pure residual directions in the denominator. It is not their
#' fiducial procedure or Burch and Iyer's optimal choice among partitions.
#' The pivot decreases with `t`. Invert its central probability interval on
#' `t >= 0`; an empty intersection remains an empty set. Analytic brackets
#' using the smallest and largest eigenvalues locate endpoints, followed by
#' checked one-dimensional root finding on the log-ratio scale. There is no
#' likelihood optimizer, arbitrary upper parameter limit or pivotal simulation.
#'
#' Each future count `m` maps the same set by `G_m = m*t/(1+m*t)`.
#' With `reference`, map the entire set through the G difference, including
#' any interior extremum; see *Differences in G between plans*. The common
#' event protects all returned plans/differences for this prespecified score.
#' G intervals have exact normal-model coverage; difference images may be
#' conservative. This is not a joint guarantee across scores or metrics,
#' or protection against nonnormality, residual dependence or selected weights.
#'
#' In a balanced source all positive eigenvalues equal the source facet
#' count, recovering the usual `"exact_f"` limits. In an incomplete source,
#' removing facet directions can discard useful information. Exact coverage
#' does not establish optimal precision: sparse, small sources can give very
#' wide intervals. No remaining residual contrasts, unresolved projection
#' rank/variation and failed roots are reported explicitly. The residual
#' combines Person-by-facet interaction and observation error.
#' Distinct log-ratio limits that round to the same G value are also marked
#' unavailable on that numeric scale, rather than displayed as zero uncertainty.
#' This also applies to nonconstant plan-difference images whose widths or
#' nonzero endpoints are lost to rounding. Identical-count differences are
#' exactly zero when their base set is available; they cannot bypass a failed
#' base calculation.
#' A genuine singleton `{0}` from intersecting the F set with nonnegative
#' ratios retains `Status = "Upper endpoint at zero"`; it is distinct from
#' coincident finite bounds caused by lost numerical resolution.
#' Extremely small confidence levels can make the two F quantiles numerically
#' indistinguishable. On the log-F scale, endpoint tolerance is the smaller
#' of `1e-8` and `1e-5` times their separation. If it falls below
#' `64 * .Machine$double.eps * max(1, abs(log_quantiles))`, the interval is
#' unavailable with an explicit reason. This is an arithmetic-resolution
#' check, not a small-sample restriction or a test of the model assumptions.
#'
#' `Estimate` retains the original ANOVA/MINQUE G or difference, including
#' unavailable points; no replacement center or sampling SE is supplied.
#' `LogVarianceRatioLower/Upper` concern log(Person variance/residual variance),
#' where `-Inf` denotes zero. These differ from the mean-square-ratio fields
#' of the balanced F method. The `pivot` attribute saves eigenvalues,
#' normalized squared coordinates, residual mean square, df and root checks.
#' `endpoint_tolerance` and `quantile_resolution` record the numerical checks.
#' Source-matched summaries and ribbons can reuse saved results.
#'
#' @section Residual-contrast SEM:
#' For a one-facet source, relative measurement SEM at future count `m` is
#' `sigma_E / sqrt(m)`. Thus `method = "residual_contrast"` first obtains
#' one residual variance interval using the projection in
#' [confint.mfrm_multivariate_gstudy()], then takes square roots and divides
#' by `sqrt(m)`. This applies to incomplete or unequal observed assignments
#' that retain residual degrees of freedom. The future plan is still a
#' complete balanced assignment; its counts do not describe the sparse source.
#'
#' Under independent homoskedastic normal residuals and positive residual
#' variance, the residual chi-square law is exact. Every plan interval is a
#' monotone image of that same event. With `reference`, multiply the common
#' residual-SD interval by `1/sqrt(comparison count)-1/sqrt(reference count)`;
#' reverse endpoints for a negative multiplier. The common event protects all
#' returned plans/differences jointly for this one prespecified score/composite.
#' Identical-count differences are zero whenever the base interval is available.
#' This does not provide a joint guarantee across different scores, metrics
#' or adaptively chosen weights, nor robustness to selective omissions,
#' nonnormal residuals or residual dependence. The source, not the future
#' count, supplies the degrees of freedom and information.
#'
#' `Estimate` and `PointStatus` retain the original ANOVA/MINQUE SEM or
#' difference. `IntervalCenter` is the separately projected SEM or difference;
#' both and the endpoints use score units. A negative original residual
#' estimate can leave `Estimate` unavailable while the projected interval
#' remains available. The original point is not replaced by the projected
#' center or zero. No sampling SE of the original SEM estimator is supplied.
#' The `base_interval` attribute retains the residual **variance** interval
#' in squared score units, while `residual_projection` retains df and checks.
#' No remaining contrasts, numerically unresolved variation or a failed base
#' interval leaves limits unavailable, including for identical-count pairs.
#'
#' Absolute SEM and multifacet SEM involve additional components; this
#' transformation does not supply their intervals or G/Phi intervals. Use
#' `summary(d, intervals = TRUE, parm = "RelativeSEM", method =
#' "residual_contrast", assumption = "normal")` to collect supported results,
#' or the same method in `plot(d, type = "sem", metric = "RelativeSEM",
#' intervals = TRUE, assumption = "normal")` for curve ribbons. Save the
#' complete result to reuse it without repeating the projection.
#'
#' @section Nested plans:
#' Suppose raters are nested within tasks: different raters serve each task,
#' and every person is rated by all of them. If the source has three raters
#' per task and every future plan also has three, task counts may vary while
#' using `confint()` for G, Phi and their differences. This does not make the
#' raters or tasks fixed. New random tasks and their random raters remain the
#' generalization population.
#'
#' Let `R` be the balanced source child count per parent. Within-cell means
#' over the children have effective components `vP`, `vParent + vChild/R`,
#' and `vP:Parent + vResidual/R`. Their mean squares are
#' `(MS_P, MS_Parent, MS_P:Parent)/R`, with df
#' `(n-1, k-1, (n-1)*(k-1))`, where `k` is the source parent count. They have
#' the same independent normal-model sampling law as a one-facet G-study.
#' Reuse the F/GCI calculations on this triplet and vary the parent count;
#' no raw-data averaging or refitting is performed by `confint()`.
#' For G, the residual-MS term in Wong and McGraw's (1999, p.277, equation
#' 2.1.3) vanishes when the source/future child-count ratio is one. The GCI
#' route is a reduction to the existing one-facet procedure, not a nested
#' GCI attributed to that paper. G remains exact under the positive-error
#' normal model; Phi retains approximate coverage and Monte Carlo diagnostics.
#'
#' Changing the child count alters these effective components, so this
#' coefficient reduction no longer applies. Individual SEM intervals still
#' support changed child counts. For SEM differences, plans need only keep
#' the same future child count as each other, which may differ from `R`:
#' their common error SD scales by the inverse square root of the parent
#' count. The base interval may then use signed-sum MLS, retaining its
#' approximate coverage. Comparisons changing both nested counts are not
#' supported by this default. Explicit `method = "mean_square_region"`
#' comparisons permit both counts to change using the enclosure below.
#' `ReferenceCount` in the reduction counts parents; general region comparisons
#' instead retain all reference counts in `attr(result, "reference")$counts`.
#' `nesting_context` retains children per parent and the original source design.
#'
#' @section General nested coefficient calculations:
#' Both counts can change with explicit `"approximate_f"` or
#' `"mean_square_region"`. Use the five components in order Person, parent,
#' child(parent), Person:parent, residual. For future child/parent counts
#' `R,T`, relative error weights are `(0,0,0,1/T,1/(R*T))`;
#' absolute error adds `(0,1/T,1/(R*T),0,0)`.
#' The mean-square region has five independent mean squares rather than seven.
#'
#' For the F calculation let `n,r,t` be source persons, children per parent
#' and parents, and set `m = r/R`, `h = t/T`. The statistic is
#' `F(rho) = (1-rho)*MS_P/D(rho)`. For G,
#' `D = (1-rho+rho*h)*MS_P:parent + rho*h*(m-1)*MS_residual`.
#' For Phi, the coefficients of `MS_parent, MS_child(parent), MS_P:parent,
#' MS_residual` are respectively
#' `rho*h/n`, `rho*h*(m-1)/n`, `1-rho+rho*h*(1-1/n)` and
#' `rho*h*(m-1)*(1-1/n)`. These follow by substituting the expected mean squares
#' into the population coefficient, then isolating `MS_P`.
#'
#' This agrees with Wong and McGraw's (1999, p.277) relative statistic.
#' Their printed absolute formulas on p.278 omit the `-rho*h/n` term in the
#' Person:parent coefficient and `+rho*h/n` in the residual coefficient.
#' The implementation uses the expected-mean-square derivation above, checked
#' against the component target and the unchanged-child one-facet reduction,
#' rather than transcribing those printed formulas or their divided endpoints.
#' No author-confirmed correction is claimed. Freeze the Satterthwaite df at
#' the unrestricted ANOVA coefficient and invert the F inequalities, as for
#' the crossed case. Signed denominator terms, negative components and invalid
#' moments remain substantive limitations, not reasons to clip estimates.
#'
#' A bounded normal-model evaluation reused 2,000 mean-square samples for each
#' of three nested sources: 10 persons/2 children/3 parents, 60/4/8, and a
#' 10/2/3 boundary model with only residual variance. Each had three future
#' plans, with decreased, unchanged and increased child/parent counts.
#' In the first two sources, nominal 95% approximate-F coverage was
#' 93.55--96.45%, with 98.9--100% availability. In the boundary source projected
#' to 4 children and 6 parents, coverage and availability both fell to 62.85%
#' for G and 61.95% for Phi. Coverage retains unavailable replications in its
#' denominator. A separate 200-sample region check gave 99.5--100% availability;
#' some median coefficient widths exceeded .99. These limited results do not
#' qualify either method as a generally useful automatic default.
#'
#' @section Mean-square region calculation:
#' For `J` independent normal-model ANOVA strata, each
#' `df[j]*MS[j]/lambda[j]` has a chi-square distribution. Use marginal
#' confidence level `level^(1/J)` for every expected mean square `lambda[j]`.
#' Independence makes the joint rectangular region's coverage `level`.
#' Intersect it with `lambda = L * theta`, where `L` is the expected-mean-square
#' map and the scalar variance components `theta` are nonnegative.
#' This retains the true component vector whenever the original region does.
#'
#' Minimize and maximize `theta[Person]/(theta[Person] + error' theta)`
#' over that common region for G/Phi. Positive residual limits and positive
#' future counts keep this denominator positive. The change of variables
#' `z = theta/denominator`, `t = 1/denominator` converts each endpoint
#' problem to a linear program. SEM instead uses the extrema of the linear
#' error variance, followed by square roots. The component scaling used
#' internally changes coordinates, not the feasible set or target.
#'
#' These projections can only increase coverage relative to the joint event
#' for one prespecified score/composite. This does not guarantee narrow
#' intervals. Empty intersections count as noncoverage in repeated sampling;
#' coverage among only the returned intervals is a different quantity.
#' The theoretical statement presumes correct projections; solver failures
#' and finite-precision limitations remain separate availability conditions.
#'
#' Extremely small confidence levels can make F quantiles or chi-square
#' region limits numerically indistinguishable. The approximate-F calculation
#' requires separation of log-F quantiles greater than
#' `64 * .Machine$double.eps * max(1, abs(log_quantiles))`.
#' A region requires positive bounds with relative width `(upper-lower)/upper`
#' greater than `64 * .Machine$double.eps` for each expected mean square
#' before optimization. Missing limits and `Status` explain failures;
#' `diagnostics` retains F quantiles and resolution, or normalized
#' `mean_square_bounds`, `relative_width`, `region_resolution` and the
#' per-stratum `tail` probability.
#' A SEM transformation that loses a nonzero endpoint or positive width is
#' also unavailable. These checks retain the original point estimates and
#' do not impose a minimum observation count.
#'
#' @section Differences from a shared mean-square region:
#' With `reference` and `method = "mean_square_region"`, project the **same**
#' simultaneous source region onto the reference and comparison targets.
#' If their intervals are `[L_a,U_a]` and `[L_b,U_b]`, report the enclosure
#' `[L_b-U_a,U_b-L_a]`. The coverage justification is their common joint event,
#' not independence or marginal 95% coverage. This is generally wider than
#' the tight image of the difference over the region; it is not labeled exact.
#'
#' Intersect with any sign restriction implied by componentwise ordering of
#' the error weights. For example, if every error weight decreases, G/Phi
#' cannot decrease and SEM cannot increase.
#' For G/Phi, additionally intersect with the analytic population range.
#' For positive reference/comparison errors `ea, eb`, the signed extremum of
#' `U*(ea-eb)/((U+ea)*(U+eb))` occurs at `U=sqrt(ea*eb)` and equals
#' `tanh(log(ea/eb)/4)`. Across nonnegative components, `ea/eb` lies between
#' the smallest and largest active error-weight ratios `a[j]/b[j]`.
#' Transform those limits and include zero to obtain `FeasibleLower/Upper`.
#' These are parameter-space limits, not extra source information. Their
#' intersection preserves the common-region coverage argument, including
#' when an extremum is approached only as a limit.
#'
#' Equal error weights give `[0,0]`
#' only when the region projections are available. Preserve missing intervals
#' and diagnostics even for identical plans. The common region event gives
#' coverage of at least `level` jointly over the returned comparisons for
#' one fixed score/composite under positive mean-square expectations and
#' independent normal effects. It gives no distribution-robust guarantee,
#' simultaneous guarantee across selected scores, or tightness guarantee.
#' The returned `individual_intervals` and `reference` retain the projections.
#' Original ANOVA/MINQUE points remain unchanged, even if outside the
#' population range implied by nonnegative components.
#'
#' @section Point estimates and intervals at boundaries:
#' `PointStatus` describes the reported point estimate; `Status` describes
#' the interval calculation. A negative ANOVA estimate of Person variance
#' can make the point coefficient unavailable while a confidence interval
#' still exists. Inspect both statuses; do not remove such an interval merely
#' because `Estimate` is `NA`. An empty confidence set and a numerical failure
#' are different outcomes, and neither means zero uncertainty.
#'
#' A wide interval is retained: it may reflect limited information about the
#' sampled persons, raters or tasks, as well as a conservative method. More
#' numerical draws cannot supply information missing from the rating design.
#' Clipping negative variance estimates or rejecting simulation draws changes
#' a statistical procedure and does not by itself establish valid coverage.
#'
#' @section Choosing a two-facet method:
#' The approximate-F, mean-square-region, score and REML alternatives remain
#' explicit choices alongside the default MLS method. Their coverage theory,
#' precision and finite-sample evidence differ; see the method-specific sections.
#' `Available` describes whether endpoints were computed; `Validation` and
#' `NegativeComponent` describe separate limitations. A narrower interval is
#' not automatically better. Select methods before inspecting which is most
#' favorable, and retain unavailable rows when comparing performance.
#'
#' - Wong and McGraw (1999, pp. 275--276) give approximate F intervals for
#'   a fully crossed Person-by-two-facet random design. Both G and Phi use
#'   estimated Satterthwaite degrees of freedom. The authors advise against
#'   applying the method when variance-component estimates are negative.
#'   `method = "approximate_f"` implements the crossed-design calculation and
#'   the independently derived nested version above.
#'   Calculable intervals with negative components are retained as an
#'   experimental extension beyond that advice, with the flag set. Small
#'   normal-model evaluations found limitations in coverage, interval
#'   width and availability, particularly at zero components and when future
#'   facet counts change. See the G-theory vignette for conditions and results.
#' - `method = "mean_square_region"` projects a simultaneous chi-square
#'   confidence region for the seven crossed or five nested independent expected mean squares onto
#'   nonnegative variance components. It is a conservative comparison based
#'   on normal ANOVA distribution theory, not a method attributed to Wong and
#'   McGraw. It can return almost the whole `[0,1]` range. It requires lpSolve;
#'   zero mean squares, infeasibility and numerical failures remain explicit.
#' - `method = "reml_profile"` instead profiles the G/Phi target in a normal
#'   restricted likelihood, with approximate chi-square calibration. It also
#'   accepts identifiable incomplete/unequal sources. See the REML section
#'   above for the separate center, diagnostics and boundary limitations.
#'
#' Use `plot(result)` for a forest comparison and `plot_data()` for its notes.
#' - For balanced crossed G, `method = "score"` uses the extended-nuisance
#'   score construction above. Its treatment of nuisance component boundaries
#'   differs from the nonnegative REML profile; neither method is uniformly
#'   superior in the available finite-sample checks. Phi is unsupported.
#'
#' To show one method along a planning curve, use
#' `plot(d, metric = "Phi", intervals = result)` with matching score/composite
#' selection; see [plot.mfrm_multivariate_d_study()].
#' Agreement of methods on one dataset does not establish robustness or
#' frequentist coverage. The G-theory vignette includes small evaluation
#' results with their denominators and Monte Carlo uncertainty.
#'
#' @section Delta coefficient SEs and intervals:
#' Use `method = "delta"` to quantify uncertainty in the original G/Phi
#' coefficient without fitting a different estimator. Let \eqn{U} be estimated
#' universe variance and \eqn{E=a^T\hat\theta} the relative or absolute error
#' variance for the selected plan, where \eqn{\hat\theta} contains the original
#' component estimates. For \eqn{\hat\rho=U/(U+E)}, the gradient is
#' \eqn{g=\{(1-\hat\rho)e_U-\hat\rho a\}/(U+E)}, where \eqn{e_U}
#' selects the universe component. The estimated sampling variance is
#' \eqn{g^T\hat C g}, with \eqn{\hat C} from [vcov.mfrm_multivariate_gstudy()].
#' The full covariance matters: estimated universe and error variances are
#' not independent. Algebra uses normalized component scales to preserve
#' invariance to score units and nonzero multipliers of fixed composite weights.
#'
#' Balanced scalar ANOVA uses `A %*% diag(2*MS^2/(df+2)) %*% t(A)` for
#' \eqn{\hat C}. This agrees with the random crossed derivation of Almehrizi
#' (2021, equations 21--25). Nested and incomplete sources use the same
#' differentiable ratio with their own design-specific component covariance.
#' Those extensions follow the delta calculation; they were not evaluated in
#' that paper. Its simulations assess SE accuracy, not Wald interval coverage.
#' The covariance correction does not make a nonlinear coefficient SE unbiased.
#' A normal-reference justification also requires growing information for the
#' components that affect the target; more source Persons alone need not suffice.
#' For example, with a fixed number of random raters and positive rater variance,
#' increasing Persons does not consistently estimate the rater-population
#' variance entering Phi. The rater mean-square degrees of freedom stay fixed.
#' This limitation differs from the exact F result for one-facet G, which
#' eliminates that nuisance component under its stated normal model.
#'
#' `SE` is coefficient sampling uncertainty, not measurement SEM. Bounds are
#' `Estimate +/- qnorm((1 + level)/2) * SE`. They are returned without clipping
#' to `[0,1]`, even when wide, and summaries, forests and ribbons retain them.
#' `OutsideParameterSpace` and `AtBoundary` identify limitations of the normal
#' approximation; they do not mean the population coefficient can exceed its
#' parameter space. Small facet pools, weak information, boundaries and
#' nonnormal effects can impair coverage. Increasing future plan counts does
#' not increase the source study's information. Intervals are pointwise for
#' prespecified plans and weights, not simultaneous or selection-adjusted.
#'
#' The original point must be defined. A nonfinite or negative projected
#' sampling variance retains its status and does not supply an SE or interval.
#' Zero first-order variance retains `SE = 0`, but does not justify a degenerate
#' normal-reference interval. Positive SEs are also retained when the requested
#' interval width rounds to zero; the limits are then `NA` with a
#' numeric-resolution status. An indefinite estimated component covariance is
#' flagged; a finite positive projection still returns an explicit approximation.
#' Nothing is repaired to obtain a preferred interval. For paired two-crossed-
#' facet differences, [mfrm_multivariate_d_compare()] retains its separately
#' documented plug-in covariance, rather than the corrected covariance used here.
#'
#' A limited normal-model evaluation reuses 2,000 independent samples per
#' source below, with a fixed outcome-independent omission pattern retaining
#' about 8/11 of the cells. Person and residual variances are 1; all other
#' component variances are 0.2. Future plans are complete balanced ratings at
#' the source pool counts; nested child counts are per parent. Entries are
#' percentages. Coverage uses **all 2,000 attempts**, counting unavailable
#' intervals as not covering; return rates are shown separately.
#'
#' | Source persons x measurement counts | G returned | G covered | Phi returned | Phi covered |
#' | --- | ---: | ---: | ---: | ---: |
#' | One facet: 10 x 3 | 92.00 | 83.55 | 92.60 | 83.60 |
#' | One facet: 60 x 8 | 100.00 | 95.65 | 100.00 | 95.25 |
#' | Crossed: 10 x 2 x 3 | 78.95 | 78.30 | 86.45 | 82.90 |
#' | Crossed: 30 x 4 x 5 | 100.00 | 93.90 | 100.00 | 91.10 |
#' | Nested: 10 x 2 within 3 | 97.30 | 86.40 | 97.85 | 87.10 |
#' | Nested: 40 x 3 within 6 | 100.00 | 94.55 | 100.00 | 91.75 |
#'
#' Coverage Monte Carlo SEs range from 0.46 to 0.92 percentage points.
#' Conditional coverage can differ materially: crossed 10 x 2 x 3 G covers
#' in 99.18% of its 1,579 returned intervals, despite only 78.30% across all
#' attempts. These conditions demonstrate why the approximation is explicit,
#' not a generally calibrated default; they do not establish a sample-size
#' cutoff or robustness to nonnormality or selective omissions.
#'
#' @section Other methods and interpretation limits:
#' Data-level bootstrap intervals are not yet options in this function.
#' Incomplete-source G/Phi can use explicit `"delta"` or `"reml_profile"`,
#' with different estimators and approximations. One-facet G also supports
#' explicit `"spectral_f"` when residual degrees of freedom are positive.
#' These routes do not supply a general interval for every original MINQUE
#' component or validate all other incomplete-source targets.
#'
#' Parameter-distribution Monte Carlo and data-generating bootstrap methods
#' require a specified estimator, joint component covariance or refitting
#' procedure, and handling of boundaries and failed calculations. A result
#' for another estimator or sampling design cannot be transferred solely
#' because its coefficient is also named G or Phi.
#'
#' A fixed set of tasks scored by every sampled rater is a different target
#' from generalization to new tasks. Representing those tasks as fixed score
#' columns with prespecified weights can use the supported one-random-facet
#' route; see [mfrm_g_study()]. This corresponds to the crossed
#' fixed-facet G construction of Zhou et al. (2011). Merely holding a random
#' task count constant does not make that facet fixed.
#'
#' An unsupported method/design combination does not imply
#' that the coefficient has no uncertainty or that interval estimation is
#' mathematically impossible. For supported calculations, inspect `Status`
#' and the returned notes: a computable approximation can carry limitations,
#' whereas an undefined ratio, empty confidence set or failed numerical
#' calculation cannot be repaired by adding a warning. Preserve available
#' estimates and diagnostics; do not substitute measurement SEM, combine
#' separate component endpoints, or drop a random facet to obtain an interval.
#' The paired two-facet intervals from [mfrm_multivariate_d_compare()] concern
#' differences between plans and have their own approximation limits. They
#' do not supply an individual-plan coefficient interval.
#'
#' @section Two-facet calculations:
#' Write the seven ANOVA component estimates as `theta_hat = A * MS`, in
#' order Person, facet 1, facet 2, Person:1, Person:2, 1:2, residual. For future
#' counts `r,t`, relative error weights are `(0,0,0,1/r,1/t,0,1/(r*t))`;
#' absolute error adds `(0,1/r,1/t,0,0,1/(r*t),0)`. The population target is
#' `rho = theta_P / (theta_P + e' theta)`, with nonnegative components and
#' positive denominator. Fixed score weights are applied before this calculation.
#'
#' For approximate F, let `a` be the Person row of A and `d = e' A`.
#' Excluding the Person mean square, set `b0 = -a/a_P` and
#' `b1 = (a+d)/a_P`. The denominator is `D(rho) = sum((b0+rho*b1)*MS)`.
#' Estimate its Satterthwaite degrees of freedom at the unrestricted ANOVA
#' `rho_hat`: `nu = D(rho_hat)^2 / sum(((b0+rho_hat*b1)*MS)^2/df)`.
#' Freeze nu and invert `F(rho) = (1-rho)*MS_P/D(rho)` between the F quantiles
#' with Person df and nu. Intersect the original linear inequalities, including
#' positive D, with `[0,1]`. Dividing endpoint formulas and clipping can be
#' wrong when future counts change the divisor sign. `diagnostics` retains
#' the unrestricted roots and df; roots alone are not confidence limits.
#'
#' For the mean-square region, let `L = solve(A)` and `lambda = L*theta`.
#' Each `df*MS/lambda` has an independent chi-square distribution. Set each
#' central mean-square interval's coverage to `level^(1/K)`, where K is seven
#' for crossed or five for nested designs, making joint
#' coverage exactly `level` under positive expectations and normal effects.
#' Intersect this rectangular region for lambda with `theta >= 0`, then
#' minimize and maximize rho by a Charnes-Cooper linear-fractional
#' transformation. Feasibility and objective residuals are checked after
#' optimization. Its image has coverage at least `level`, jointly over future
#' plans for this fixed score; it is often conservative. This guarantee does
#' not extend to nonnormal effects or data-selected weights, or establish
#' numerical validity of a failed solve. Negative ANOVA estimates are retained
#' in the source and are not individually clipped before either method.
#'
#' @section Calculation of G:
#' Let `k` be the observed number of facet levels and `m` a future plan's count.
#' Person and residual mean squares `MS_P` and `MS_E` are independent under
#' the stated model, with degrees of freedom `n - 1` and `(n - 1)*(k - 1)`.
#' Their expectations are `k * sigma_P^2 + sigma_E^2` and `sigma_E^2`.
#' Thus, writing `R = 1 + k * sigma_P^2 / sigma_E^2`,
#' `(MS_P / MS_E) / R` follows that F distribution. Divide the observed ratio
#' by the upper and lower F quantiles to obtain the interval for `R`.
#' Intersect it with `R >= 1`, then transform its endpoints using
#' `G_m = m * (R - 1) / (m * (R - 1) + k)`.
#' These steps follow directly from the independent mean-square distributions
#' in Mondal et al. (2025, Section 2, Table 1); the transformation here expresses
#' the one-facet G-theory target at the requested future count. At `m = k`,
#' the usual estimate is `1 - MS_E / MS_P`.
#'
#' Intersect **before** transforming: outside `R >= 1`, the transformation
#' can have a pole when `m > k`. Clipping already transformed endpoints to
#' `[0, 1]` would therefore be incorrect. `RawRatioLower` and `RawRatioUpper`
#' retain unrestricted ratio endpoints; they are not G endpoints.
#'
#' An interval reaching `R = 1` includes G = 0. If the ratio interval is
#' entirely below 1, the valid intersection is empty and both G endpoints
#' are `NA`, with `Status = "Empty parameter-space interval"`. Negative
#' Person estimates are retained in the source; `Estimate` and `PointStatus`
#' follow the existing D-study rules. A point estimate may be unavailable
#' while a nonempty interval still exists. Zero residual mean square leaves
#' the interval unavailable, even if the point estimate is 1: the F pivot's
#' positive-variance condition cannot be checked from a zero denominator.
#' If a positive width rounds away on the G scale, endpoints are unavailable
#' with a numeric-resolution status. The retained log-ratio bounds can still
#' support representable plan-difference intervals. If the base log-ratio
#' endpoints themselves coincide through rounding (for example at an extremely
#' small `level`), both G and plan-difference intervals are unavailable. This
#' numerical failure differs from a positive-width base set whose valid
#' intersection is the single boundary value G = 0.
#'
#' All rows use the same ratio interval for the selected score. Under this
#' model their coverage events are identical; they are not independent
#' intervals, and their overlap is not a test of a plan difference. This
#' property does not provide simultaneous inference over different scores,
#' adaptively chosen weights, other designs or other metrics. Existing
#' D-study plots can display saved intervals through their `intervals` argument.
#'
#' @section Differences in G between plans:
#' Use `reference = 1` to ask how much G would change from the first plan.
#' The plans use the same estimated components: subtracting independent
#' marginal confidence limits or adding marginal variances is inappropriate.
#' This method instead maps the **entire** confidence set for the common
#' variance ratio to the G difference; it does not use a delta approximation.
#'
#' Put `t = sigma_P^2 / sigma_E^2 = (R - 1)/k`. For reference count `a` and
#' comparison count `b`, the difference is
#' `h(t) = (b-a)*t / ((1+a*t)*(1+b*t))`.
#' Its derivative has the sign of `(b-a)*(1-a*b*t^2)`. Therefore the range
#' over the nonnegative ratio interval is obtained from both endpoints and,
#' when included, the interior stationary point `t = 1/sqrt(a*b)`.
#' Mapping only the endpoints would miss the extremum. At zero and at the
#' limit `t -> Inf`, the difference is zero. Equal counts give zero throughout
#' a nonempty confidence set. An empty base set or a zero residual mean square
#' retains its unavailable status, including for repeated counts.
#' A nonconstant image whose width or nonzero endpoint is lost to numerical
#' rounding has missing limits with `Status` explaining the loss; its original
#' point remains unchanged. This numerical rule also applies to shared Phi
#' and REML coefficient-set images, with each method's own coverage scope.
#'
#' Whenever the base F confidence set contains the true variance ratio, all
#' its mapped difference intervals contain the corresponding true differences.
#' Thus coverage is **at least** the requested level under the stated model,
#' including jointly over the returned prespecified counts for this one fixed
#' score/composite. It can be conservative because `h` is not one-to-one.
#' This follows from confidence-set inclusion and the derivative above, not
#' from a new asymptotic distribution. It does not supply inference for other
#' scores, adaptive score weights or other facet designs. Phi uses the separate
#' generalized construction below, with no exact nominal-coverage guarantee.
#' The sign of the change from adding levels already follows from this model;
#' the interval quantifies its size, not whether more levels are cost-effective.
#' A point difference can be unavailable while its interval remains available.
#'
#' @section Calculation of Phi:
#' With `n` persons and `k` observed facet levels, write `Q_P`, `Q_F`, `Q_E`
#' for `df * MS / X`, where the three `X` values are independent chi-square
#' draws with degrees of freedom `n-1`, `k-1`, `(n-1)*(k-1)`, respectively.
#' The single-rating agreement-ICC generalized pivot is
#' `max(0, Q_P - Q_E) / (Q_P + (k/n)*Q_F + (k - 1 - k/n)*Q_E)`.
#' This is Equation 11 of Mondal et al. (2025), including their correction
#' to an earlier published denominator. Its central quantiles give the
#' single-rating GCI. For a future count `m`, transform both endpoints by
#' `m * rho / (1 + (m - 1) * rho)`, since
#' `Phi_m = sigma_P^2 / (sigma_P^2 + (sigma_F^2 + sigma_E^2)/m)`.
#' The observed `k`, not future `m`, determines the sampling distributions.
#'
#' Truncation of the pivot's numerator at zero is part of this GCI method;
#' it does not truncate the source component estimates or change the point
#' estimate. A negative Person estimate can therefore have an unavailable
#' point estimate and an available interval. A zero upper endpoint can also
#' occur; `Status` and the proportion of zero pivotal draws flag this boundary
#' result, which must not be read as proof of zero sampling uncertainty.
#' Any zero mean square makes this implementation return an unavailable
#' interval with `Status = "Degenerate mean squares"`.
#'
#' The default 100000 draws follow the numerical recommendation discussed in
#' Mondal et al. (Section 3.3). `attr(result, "monte_carlo")` contains draw/seed
#' settings, RNG kind, the zero proportion, and separate 95% Monte Carlo bands
#' for each endpoint, obtained from binomial order-statistic ranks. These bands
#' describe numerical quantile error, not another confidence interval for Phi.
#' They do not include uncertainty from violated model assumptions or approximate
#' frequentist coverage. All future counts share the same sampled interval;
#' rows are dependent and do not test plan differences. Fixed composites are
#' supported, but adaptive weights and simultaneous inference over scores are not.
#'
#' @section Differences in Phi between plans:
#' With `parm = "Phi", reference = 1`, ask how much absolute dependability
#' changes from the first plan. Write `t = sigma_P^2/(sigma_F^2+sigma_E^2)`.
#' Then `Phi_m = m*t/(1+m*t)`, so the count-difference function `h(t)` and
#' its interior extremum are the same as for G, with a different variance ratio.
#' Transform the whole single-rating generalized interval by
#' `t = rho/(1-rho)` and take the range of `h`. Zero and unit endpoints
#' correspond to `t = 0` and its infinite limit. This set-image extension
#' follows from the displayed algebra; it is not a new published difference
#' pivot or the quantiles of differences calculated draw by draw.
#'
#' Whenever the base generalized interval contains the true ratio, its image
#' contains every returned true plan difference for this fixed score. Mapping
#' can increase coverage through conservatism but does not make the base GCI's
#' approximate frequentist coverage exact. Normal effects and prespecified
#' counts/weights remain required. Do not infer the difference interval by
#' subtracting two marginal confidence limits. Use `plot(result)` for a forest
#' view on the difference scale; these results cannot be used as coefficient
#' ribbons in `plot(d)`.
#'
#' One pivotal sample supplies all comparisons. The `monte_carlo` attribute
#' retains `base_endpoint_bands` and at least 95% **joint numerical** bands
#' for both difference endpoints across the returned plans. Each base endpoint
#' has a Bonferroni-protected 97.5% order-statistic band. Map the outer base
#' interval and their common inner interval to bound the image endpoints;
#' when the inner interval is empty, use the outer image as a coarse bound
#' for each endpoint. These bands concern numerical quantile error, not the
#' probability that a coefficient interval covers its population target.
#' Degenerate/numerically unavailable base intervals remain unavailable even
#' for identical counts. A base `Upper endpoint at zero` status is retained
#' for unequal counts; it is not proof that the difference is precisely zero.
#'
#' @return A data frame of class `mfrm_dstudy_intervals`, with scenario/facet
#'   counts, `Kind`, `Score`, `Metric`,
#'   `Estimate`, `PointStatus`, `Lower`, `Upper`, unrestricted `RawRatioLower`
#'   and `RawRatioUpper` (balanced F/GCI coefficient methods), `Level`, `Method`, interval
#'   `Status`, `MethodCode`, `Available`, `NegativeComponent`, and `Validation`.
#'   Each selected method supplies a row per plan, including unavailable rows.
#'   Delta results additionally supply `SE`, `SamplingVariance`, `SEStatus`,
#'   `AtBoundary`, `OutsideParameterSpace` and `SamplingCovariancePSD`.
#'   Attribute `sampling_covariance` retains the component covariance with
#'   its own diagnostics and weights divided by their largest absolute value;
#'   `weights` retains the original score definition. Row-aligned `diagnostics` records normalized gradients,
#'   error coefficients and scales. MINQUE results retain `source_data` for
#'   source-identity checks. Method comparisons including delta retain the
#'   complete per-method objects in `method_results`.
#'   Spectral F additionally supplies `LogVarianceRatioLower/Upper` and
#'   projection/root diagnostics in `pivot`; it does not supply a replacement
#'   point estimate or SE.
#'   Coefficient MLS additionally retains `NumericalComplete`, `FormulaComplete`
#'   and `Empty`. Its row-aligned diagnostics contain the null problem, quadratic
#'   coefficients, accepted pieces/nodes, and any conservatively retained
#'   `undefined` ranges. Formula completeness and numerical availability are
#'   separate from the approximate statistical calibration.
#'   Score intervals additionally retain `NumericalComplete` and `Empty` to
#'   distinguish a resolved empty numerical set from a failed calculation;
#'   either has missing accepted endpoints. Scanned pieces and root checks
#'   remain in the row-aligned `diagnostics` attribute.
#'   `plot()` compares saved intervals without recomputation. `design_grid`
#'   retains future counts. For two facets, `diagnostics` is a row-aligned list
#'   of F moments/roots or solver status/residuals; `mean_squares` stores the
#'   normalized source mean squares, their scale, df and component map.
#'   Attributes `weights`, `assumption`, `source_design`, `source_rows`,
#'   `pivot` (G's log ratio/endpoints and degrees of freedom; Phi's single-rating
#'   bounds and degrees of freedom), and `notes` retain context. For Phi,
#'   `RawRatioLower/Upper` are `NA` because there is no G-type F ratio;
#'   `monte_carlo` retains numerical precision diagnostics. Phi endpoints and
#'   their numerical bands are withheld when a positive width rounds away.
#'   If the base pivotal interval has lost resolution, plan differences also
#'   remain unavailable, including identical-count comparisons. The pivot
#'   attribute retains the diagnostic base bounds. A sampled atom at zero
#'   retains its separate `"Upper endpoint at zero"` status.
#'   Save the whole
#'   object with `saveRDS()`. `source_summary` retains the design, covariance
#'   components, mean products and degrees of freedom for checking that a curve
#'   overlay uses the same source. REML profiles additionally retain
#'   `source_data` to check the used ratings, `reml_source` with the fit and
#'   numerical diagnostics, and row-aligned profile `diagnostics` with scans,
#'   accepted pieces, optimizer attempts and gradient checks. No live fitting
#'   environment is saved. `IntervalCenter` is the separate REML estimate;
#'   `CandidateLower/Upper` can remain finite when accepted bounds are unavailable.
#'   Shared-image REML `reference` results retain the common fit and scans in
#'   `base_interval`, with row-aligned transformation data in `diagnostics`.
#'   Direct difference profiles instead retain per-comparison scans and
#'   constraints in `diagnostics`, and all reference counts in
#'   `attr(result, "reference")$counts`.
#'   Residual-contrast SEM also retains `source_data` and a separate
#'   `IntervalCenter`; `base_interval` contains the residual variance interval,
#'   `residual_projection` its numerical checks, and row-aligned `diagnostics`
#'   the count multipliers. Method comparisons that include REML or residual
#'   contrasts retain each complete single-method object in `method_results`,
#'   so distinct bases and centers are not conflated by the combined table.
#'   For MLS SEM, individual-plan tables additionally retain `ErrorVariance`,
#'   `ErrorVarianceSE`, `RawVarianceLower` and `RawVarianceUpper` in squared
#'   score units. The squared `ErrorVarianceSE` uses the unbiased
#'   `df + 2` sampling-variance formula described in
#'   [confint.mfrm_multivariate_gstudy()]; no delta-method SEM SE is supplied.
#'   `diagnostics` retains component and mean-square weights indexed by plan,
#'   and `mean_squares` retains the scale and source df.
#'   With `reference`, rows exclude the reference plan, `Metric` is
#'   the selected metric followed by `" difference"`, and `Estimate` is comparison minus reference.
#'   `Reference`, `ReferenceCount`, `ReferenceValue`, and `Value` retain the
#'   two plans' context. Attribute `reference` also retains its point status.
#'   For MLS SEM differences, it retains the reference SEM interval, raw
#'   error-variance bounds, estimate and SE in place of per-plan variance columns.
#'   Nested results also retain `nesting_context` with facet identities, the
#'   source children per parent, and both counts for every future plan.
#'   A nested reduction/MLS `reference` includes `count_facet` and
#'   `children_per_parent`.
#'   General region comparisons instead retain all reference `counts`,
#'   its interval, status and diagnostics; `individual_intervals` stores
#'   the planwise projections before forming the difference enclosure.
#'   Their `diagnostics` list aligns with comparison rows. No scalar
#'   `ReferenceCount` is invented for a comparison changing both facets.
#'   No sampling SE or covariance matrix is inferred from the interval width.
#'   Stored G-study summaries are used and point projections are recomputed
#'   without refitting for the summary-based methods; REML profiles fit the
#'   stored used ratings. Every method leaves the input unchanged.
#' @references
#' E, L., Hannig, J., & Iyer, H. (2008). Fiducial intervals for variance
#' components in an unbalanced two-component normal mixed linear model.
#' *Journal of the American Statistical Association, 103*(482), 854--865.
#' \doi{10.1198/016214508000000229}. Sections 2 and 2.3, equations 3--7,
#' give the independent quadratic forms and existing F-pivot family used
#' by `"spectral_f"`; the fiducial method developed later is not implemented.
#'
#' Shedden, M., & Ekvall, K. O. (2026). Score-based confidence intervals for
#' variance-covariance parameters in linear mixed models. Preprint.
#' \doi{10.48550/arXiv.2610.04181}. The ratio-target adaptation above is not a
#' G-theory coverage result supplied by this paper.
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software, 67*(1),
#' 1--48. \doi{10.18637/jss.v067.i01}. Restricted likelihood computation.
#'
#' Drton, M. (2009). Likelihood ratio tests and singularities.
#' *The Annals of Statistics, 37*(2), 979--1012.
#' \doi{10.1214/07-AOS571}. General cautions about singular hypotheses;
#' not a calibration theorem for the present crossed REML contrasts.
#'
#' Crainiceanu, C. M., & Ruppert, D. (2004). Likelihood ratio tests in linear
#' mixed models with one variance component. *Journal of the Royal Statistical
#' Society: Series B, 66*(1), 165--185.
#' \doi{10.1111/j.1467-9868.2004.00438.x}. Special-case boundary distributions;
#' not a calibration result for general G/Phi targets.
#'
#' Betebenner, D. W. (1998). *Improved confidence interval estimation for
#' variance components and error variances in generalizability theory*.
#' Paper presented at AERA, San Diego. ERIC ED423304.
#' \url{https://files.eric.ed.gov/fulltext/ED423304.pdf}.
#'
#' Hess, A., & Iyer, H. (2002). A SAS/IML macro for computation of confidence
#' intervals for variance components of mixed models. *Proceedings of SUGI 27*,
#' Paper 252-27, pp. 2--3.
#' \url{https://support.sas.com/resources/papers/proceedings/proceedings/sugi27/p252-27.pdf}.
#'
#' Ting, N., Burdick, R. K., Graybill, F. A., Jeyaratnam, S., & Lu, T.-F. C.
#' (1990). Confidence intervals on linear combinations of variance components
#' that are unrestricted in sign. *Journal of Statistical Computation and
#' Simulation, 35*(3--4), 135--143. \doi{10.1080/00949659008811240}.
#' The G/Phi null inversion and conservative completion above are mfrmr's
#' adaptations, with their own limited finite-sample evaluation.
#'
#' Graybill, F. A., & Wang, C.-M. (1980). Confidence intervals on nonnegative
#' linear combinations of variances. *Journal of the American Statistical
#' Association, 75*(372), 869--873.
#'
#' Mondal, D., Candel, M. J. J. M., Cassese, A., & Vanbelle, S. (2025).
#' Confidence intervals and sample size for the ICC in two-way ANOVA models.
#' *Statistics in Medicine, 44*(10--12), e70106. \doi{10.1002/sim.70106}.
#' Section 2 and Table 1 give the independent normal-theory mean squares.
#' Its later comparison of agreement-ICC methods does not establish an exact
#' Phi interval or a multi-facet G interval.
#' Section 3.3, Equation 11 and Appendix C give the GCI used for Phi here.
#'
#' Tian, L., & Cappelleri, J. C. (2004). A new approach for interval estimation
#' and hypothesis testing of a certain intraclass correlation coefficient:
#' The generalized variable method. *Statistics in Medicine, 23*, 2125--2135.
#' \doi{10.1002/sim.1782}.
#'
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#' \doi{10.1007/978-1-4757-3456-0}. G-study and D-study targets.
#'
#' Wong, S. P., & McGraw, K. O. (1999). Confidence intervals and F tests for
#' intraclass correlations based on three-way random effects models.
#' *Educational and Psychological Measurement, 59*, 270--288.
#' \doi{10.1177/00131649921969848}. Crossed two-random-facet approximate F
#' method, with the negative-component extension flagged above. Page 277
#' gives the nested relative-G pivot used in the unchanged-child-count reduction.
#'
#' Zhou, H., Muellerleile, P., Ingram, D., & Wong, S. P. (2011). Confidence
#' intervals and F tests for intraclass correlation coefficients based on
#' three-way mixed effects models. *Journal of Educational and Behavioral
#' Statistics, 36*, 638--671. \doi{10.3102/1076998610381399}.
#' The fully crossed fixed-facet relative-G target agrees with the fixed
#' score-composite, one-random-facet representation described above.
#'
#' Almehrizi, R. S. (2021). Standard errors of variance components,
#' measurement errors and generalizability coefficients for crossed designs.
#' *Journal of Educational Measurement, 58*, 179--210.
#' \doi{10.1111/jedm.12277}. Coefficient SEs and their evaluation, not a
#' validation of two-facet coefficient interval coverage.
#' @seealso [mfrm_d_study()], [confint.mfrm_multivariate_gstudy()],
#'   [mfrm_multivariate_d_compare()]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("confint.mfrm_multivariate_d_study", run.dontrun = TRUE)
#' \dontrun{
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#' d <- mfrm_d_study(g, data.frame(Tasks = c(6, 12)),
#'   weights = c(V = .5, W = .5))
#' ci <- confint(d, assumption = "normal")
#' ci[, c("Tasks", "Estimate", "Lower", "Upper", "Status")]
#' attr(ci, "notes")
#' # Optional normal approximation: SE describes the coefficient, not measurement SEM.
#' delta_ci <- confint(d, method = "delta", assumption = "normal")
#' delta_ci[, c("Tasks", "Estimate", "SE", "Lower", "Upper", "SEStatus")]
#' phi_ci <- confint(d, parm = "Phi", assumption = "normal", seed = 1)
#' phi_ci[, c("Tasks", "Estimate", "Lower", "Upper", "Status")]
#' attr(phi_ci, "monte_carlo")$endpoint_bands
#' # How much would G improve by changing from six to twelve tasks?
#' change <- confint(d, parm = "G", reference = 1, assumption = "normal")
#' change[, c("Tasks", "ReferenceCount", "Estimate", "Lower", "Upper", "Status")]
#' # The same planning question for absolute dependability; coverage is approximate.
#' phi_change <- confint(d, parm = "Phi", reference = 1, assumption = "normal")
#' phi_change[, c("Tasks", "ReferenceCount", "Estimate", "Lower", "Upper", "Status")]
#' attr(phi_change, "monte_carlo")$endpoint_bands
#' # Sampling uncertainty in measurement SEM, rather than in G or Phi.
#' sem_ci <- confint(d, parm = "AbsoluteSEM", assumption = "normal")
#' sem_ci[, c("Tasks", "Estimate", "Lower", "Upper", "Method", "Status")]
#' sem_change <- confint(d, parm = "RelativeSEM", reference = 1,
#'                       assumption = "normal")
#' sem_change[, c("Tasks", "Estimate", "Lower", "Upper", "Status")]
#' sem_plot <- plot(d, type = "sem", metric = "AbsoluteSEM", intervals = TRUE,
#'                  assumption = "normal", draw = FALSE)
#' plot_data(sem_plot, "notes")
#' # Different raters serve each task; retain three raters per task.
#' data("mfrmr_example_gtheory_nested", package = "mfrmr")
#' gn <- mfrm_g_study(mfrmr_example_gtheory_nested,
#'   c("Content", "Organization"), nesting = c(Rater = "Task"))
#' dn <- mfrm_d_study(gn,
#'   data.frame(Raters = 3, Tasks = c(6, 12)))
#' nested_change <- confint(dn, parm = "G", reference = 1, assumption = "normal")
#' nested_change[, c("Raters", "Tasks", "ReferenceCount", "Estimate", "Lower", "Upper")]
#' attr(nested_change, "nesting_context")
#' # Incomplete one-facet ratings: relative SEM without likelihood fitting.
#' observed <- tasks[(tasks$Person + tasks$Task) %% 3 != 0, ]
#' sparse_g <- mfrm_g_study(observed, "V", rater = NULL, method = "minque0")
#' sparse_d <- mfrm_d_study(sparse_g, data.frame(Tasks = c(6, 12)))
#' residual_sem <- confint(sparse_d, parm = "RelativeSEM",
#'   method = "residual_contrast", assumption = "normal")
#' residual_sem[, c("Tasks", "Estimate", "IntervalCenter", "Lower", "Upper", "Status")]
#' residual_change <- confint(sparse_d, parm = "RelativeSEM", reference = 1,
#'   method = "residual_contrast", assumption = "normal")
#' residual_change[, c("Tasks", "Estimate", "IntervalCenter", "Lower", "Upper")]
#' # Relative G also uses Person variation; normal Person effects are required.
#' spectral_g <- confint(sparse_d, method = "spectral_f", assumption = "normal")
#' spectral_g[, c("Tasks", "Estimate", "Lower", "Upper", "Status")]
#' spectral_change <- confint(sparse_d, method = "spectral_f", reference = 1,
#'   assumption = "normal")
#' spectral_change[, c("Tasks", "Estimate", "Lower", "Upper")]
#' # A separate experimental REML model can address G/Phi or absolute SEM.
#' if (requireNamespace("lme4", quietly = TRUE) && requireNamespace("minqa", quietly = TRUE)) {
#'   profile_ci <- confint(sparse_d, method = "reml_profile", assumption = "normal")
#'   profile_ci[, c("Tasks", "Estimate", "IntervalCenter", "Lower", "Upper", "Validation")]
#'   attr(profile_ci, "notes")
#'   # Uncertainty in absolute measurement SEM, in score units.
#'   profile_sem <- confint(sparse_d, parm = "AbsoluteSEM",
#'                          method = "reml_profile", assumption = "normal")
#'   plot(sparse_d, type = "sem", metric = "AbsoluteSEM", intervals = profile_sem,
#'        draw = FALSE)
#'   # How much does doubling the task count change G?
#'   profile_change <- confint(sparse_d, reference = 1,
#'     method = "reml_profile", assumption = "normal")
#'   profile_change[, c("Tasks", "Estimate", "IntervalCenter", "Lower", "Upper")]
#'   attr(profile_change, "base_interval")
#'   plot(profile_change, draw = FALSE)
#' }
#' # Start with the default two-crossed-facet approximation on bundled data.
#' data("mfrmr_example_gtheory_crossed", package = "mfrmr")
#' crossed <- mfrm_g_study(mfrmr_example_gtheory_crossed,
#'   c("Content", "Organization"))
#' future <- mfrm_d_study(crossed,
#'   data.frame(Raters = c(2, 4), Tasks = c(6, 12)),
#'   weights = c(Content = .5, Organization = .5))
#' crossed_ci <- confint(future, parm = "Phi", assumption = "normal")
#' crossed_ci[, c("Raters", "Tasks", "Estimate", "Lower", "Upper",
#'                "FormulaComplete", "Status")]
#' attr(crossed_ci, "notes")
#' # The score method remains an explicit experimental comparison for G.
#' score_ci <- confint(future, method = "score", assumption = "normal")
#' score_ci[, c("Scenario", "Estimate", "Lower", "Upper", "Empty", "Status")]
#' # Vary task counts within each rater count to obtain planning-curve ribbons.
#' curve_plans <- mfrm_d_study(crossed,
#'   expand.grid(Raters = c(2, 4), Tasks = c(3, 6, 12)),
#'   weights = c(Content = .5, Organization = .5))
#' score_curve <- plot(curve_plans, metric = "G", intervals = TRUE,
#'   method = "score", assumption = "normal", draw = FALSE)
#' plot_data(score_curve, "notes")
#' if (requireNamespace("lpSolve", quietly = TRUE)) {
#'   comparison <- confint(future, parm = "Phi", assumption = "normal",
#'     method = c("mls", "approximate_f", "mean_square_region"))
#'   comparison[, c("Scenario", "Method", "Lower", "Upper", "Status")]
#'   figure <- plot(comparison, draw = FALSE)
#'   plot_data(figure, "notes")
#'   if (requireNamespace("ggplot2", quietly = TRUE)) as_ggplot(figure)
#'   # General nested comparison: both counts change, with eighteen ratings.
#'   general_nested <- mfrm_d_study(gn,
#'     data.frame(Raters = c(3, 2, 6), Tasks = c(6, 9, 3)))
#'   nested_region <- confint(general_nested, parm = "Phi", reference = 1,
#'     method = "mean_square_region", assumption = "normal")
#'   nested_region[, c("Raters", "Tasks", "Estimate", "Lower", "Upper", "Status")]
#'   attr(nested_region, "reference")$counts
#' }
#' }
#' @export
confint.mfrm_multivariate_d_study <- function(object, parm = "G", level = .95,
                                             score = NULL, composite = NULL,
                                             assumption, draws = 100000L,
                                             seed = 1L, reference = NULL,
                                             method = "auto", ...) {
  rlang::check_dots_empty()
  if (missing(assumption) || !identical(assumption, "normal"))
    stop("This interval requires `assumption = 'normal'`; the assumption is not tested.", call. = FALSE)
  if (!is.character(parm) || length(parm) != 1L || is.na(parm) ||
      !parm %in% c("G", "Phi", "RelativeSEM", "AbsoluteSEM"))
    stop("Select one of 'G', 'Phi', 'RelativeSEM' or 'AbsoluteSEM'.", call. = FALSE)
  sem <- parm %in% c("RelativeSEM", "AbsoluteSEM")
  if (!is.numeric(level) || is.object(level) || is.complex(level) ||
      !is.null(dim(level)) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be strictly between zero and one.", call. = FALSE)
  if (!inherits(object, "mfrm_multivariate_d_study"))
    stop("`object` must be a multivariate D-study result.", call. = FALSE)
  if (!identical(.mfrm_mvds_policy(object), "raw")) {
    d <- .mfrm_mvds_replay(object)
    raw <- mfrm_d_study(d$gstudy, d$design_grid, d$weights)
    ci <- confint.mfrm_multivariate_d_study(raw, parm = parm, level = level,
      score = score, composite = composite, assumption = assumption,
      draws = draws, seed = seed, reference = reference, method = method)
    return(.mfrm_mvds_adjusted_intervals(ci, d))
  }
  g <- object$gstudy
  d <- mfrm_d_study(g, object$design_grid, object$weights)
  anova <- identical(g$design$method, "Balanced multivariate ANOVA") &&
    isTRUE(g$design$balanced) && isTRUE(g$design$complete)
  nested <- !is.null(g$design$nesting)
  child <- if (nested) d$design_grid[[g$design$count_columns[[names(g$design$nesting)]]]]
  facets <- length(g$design$counts) - 1L
  one_target <- facets == 1L || (nested && all(child == g$design$child_counts[1L]))
  if (!is.character(method) || !length(method) || anyNA(method) ||
      anyDuplicated(method) || !is.null(dim(method)))
    stop("`method` must contain distinct method names.", call. = FALSE)
  if (identical(method, "auto")) {
    if (!anova)
      stop("Automatic intervals require a balanced ANOVA source. For a MINQUE source, explicitly select 'delta' for G/Phi, 'reml_profile', 'spectral_f' for one-facet G, or 'residual_contrast' for one-facet RelativeSEM, after reviewing the method's limitations.", call. = FALSE)
    if (!sem && !one_target && !is.null(reference))
      stop("Automatic coefficient plan-difference intervals require one random measurement facet or unchanged nested child count per parent. Otherwise explicitly select 'mean_square_region' or 'reml_profile' after reviewing their limitations.", call. = FALSE)
    if (!sem && !one_target && nested)
      stop("Automatic nested coefficient intervals require unchanged child count per parent. Otherwise explicitly select 'delta', 'approximate_f', 'mean_square_region' or 'reml_profile' after reviewing their limitations.", call. = FALSE)
    method <- if (sem || !one_target) "mls" else if (parm == "G") "exact_f" else "gci"
  }
  supported <- if (!anova) "reml_profile" else if (sem) c("mls", if (facets == 2L) "mean_square_region", "reml_profile") else
    c(if (one_target) if (parm == "G") "exact_f" else "gci",
      if (facets == 2L) c("approximate_f", "mean_square_region"),
      if (facets == 2L && !nested) c("mls", if (parm == "G") "score"), "reml_profile")
  if (facets == 1L && parm == "G") supported <- c(supported, "spectral_f")
  if (facets == 1L && parm == "RelativeSEM") supported <- c(supported, "residual_contrast")
  if (!sem) supported <- c(supported, "delta")
  if (any(!method %in% supported))
    stop("Unsupported method for this design and metric. Available: ",
      paste(supported, collapse = ", "), ".", call. = FALSE)
  if (length(method) > 1L && ((sem && !any(c("reml_profile", "residual_contrast") %in% method)) ||
      (!any(c("reml_profile", "spectral_f", "delta") %in% method) && any(method %in% c("exact_f", "gci")))))
    stop("Select one method, or compare 'approximate_f' and 'mean_square_region' for coefficients.", call. = FALSE)
  if (!is.null(reference)) {
    if ("delta" %in% method)
      stop("Delta `reference` intervals are not supplied here; mfrm_multivariate_d_compare() provides separately qualified paired delta intervals for two crossed facets.", call. = FALSE)
    if (!sem && "mls" %in% method)
      stop("Coefficient MLS `reference` intervals are not supported; marginal null sets do not give a plan-difference interval.", call. = FALSE)
    if ("score" %in% method)
      stop("Score `reference` intervals are not supported; marginal score sets do not give a plan-difference interval.", call. = FALSE)
    if ("approximate_f" %in% method)
      stop("Approximate-F `reference` intervals are not supported; select 'mean_square_region' for a conservative joint-region difference enclosure.", call. = FALSE)
    if (identical(method, "mls") && facets != 1L &&
        (!nested || length(unique(child)) != 1L))
      stop("MLS `reference` comparisons require one facet or nested plans with the same future child count per parent across plans. Otherwise select 'mean_square_region'.", call. = FALSE)
  }
  if (!is.null(reference)) {
    if (!is.numeric(reference) || is.object(reference) || is.complex(reference) ||
        !is.null(dim(reference)) || length(reference) != 1L || !is.finite(reference) ||
        reference != floor(reference) || nrow(d$design_grid) < 2L ||
        reference < 1L || reference > nrow(d$design_grid))
      stop("`reference` must name a valid row in a design grid with at least two plans.", call. = FALSE)
  }
  tab <- d$coefficients
  if (!is.null(score) && !is.null(composite))
    stop("Choose either `score` or `composite`.", call. = FALSE)
  if (is.null(score) && is.null(composite)) {
    choices <- unique(tab$Score[tab$Kind == "Composite"])
    if (length(choices) > 1L) stop("Select one `composite` or original `score`.", call. = FALSE)
    if (length(choices)) composite <- choices else score <- g$design$scores[1L]
  }
  kind <- if (is.null(composite)) "Score" else "Composite"
  selected <- if (is.null(composite)) score else composite
  if (!is.character(selected) || length(selected) != 1L || is.na(selected) ||
      !any(tab$Kind == kind & tab$Score == selected))
    stop("Select one score or composite present in the D-study.", call. = FALSE)
  w <- if (kind == "Score") setNames(as.numeric(g$design$scores == selected), g$design$scores) else
    if (is.matrix(d$weights)) d$weights[, selected] else d$weights
  tab <- tab[tab$Kind == kind & tab$Score == selected, , drop = FALSE]
  if (length(method) > 1L && any(c("reml_profile", "residual_contrast", "spectral_f", "delta") %in% method))
    return(.mfrm_mvgt_profile_bind(lapply(method, function(code)
      confint.mfrm_multivariate_d_study(d, parm = parm, level = level, score = score,
        composite = composite, assumption = assumption, draws = draws, seed = seed,
        reference = reference, method = code))))
  if (identical(method, "delta"))
    return(.mfrm_mvgt_delta_table(d, tab, w, parm, level))
  if (identical(method, "spectral_f"))
    return(.mfrm_mvgt_spectral_table(d, tab, w, level, reference))
  if (identical(method, "residual_contrast"))
    return(.mfrm_mvgt_residual_sem(d, tab, w, level, reference))
  if ("reml_profile" %in% method) {
    if (!is.null(reference)) {
      if (facets == 1L || (nested && length(unique(child)) == 1L))
        return(.mfrm_mvgt_profile_difference(d, tab, w, parm, level, reference))
      return(.mfrm_mvgt_contrast_table(d, tab, w, parm, level, reference))
    }
    return(.mfrm_mvgt_profile_table(d, tab, w, parm, level))
  }
  p <- .mfrm_mvgt_anova_projection(g, weights = w)
  if (!sem && "mls" %in% method) {
    if (length(method) > 1L) return(.mfrm_mvgt_profile_bind(lapply(method, function(code)
      confint.mfrm_multivariate_d_study(d, parm = parm, level = level, score = score,
        composite = composite, assumption = assumption, method = code))))
    return(.mfrm_mvgt_mls_coefficient_table(d, tab, p, parm, level))
  }
  if ("score" %in% method) {
    if (length(method) > 1L) return(.mfrm_mvgt_profile_bind(lapply(method, function(code)
      confint.mfrm_multivariate_d_study(d, parm = parm, level = level, score = score,
        composite = composite, assumption = assumption, draws = draws, seed = seed,
        method = code))))
    return(.mfrm_mvgt_score_table(d, tab, p, level))
  }
  if (identical(method, "mls")) return(.mfrm_mvgt_sem_table(d, tab, p, parm, level, reference))
  if (all(method %in% c("approximate_f", "mean_square_region")))
    return(.mfrm_mvgt_two_facet_intervals(d, tab, p, parm, level, method, reference))
  one <- .mfrm_mvgt_one_facet_target(d, p)
  if (parm == "Phi") return(.mfrm_mvgt_interval_result(
    .mfrm_mvgt_phi_table(d, tab, p, level, draws, seed, reference, one), d, p, method))
  interval <- .mfrm_mvgt_g_f_interval(one$ms[1L], one$ms[3L],
    one$df[c(1L, 3L)], one$observed, one$future, level)
  out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
  out$Metric <- "G"; out$Estimate <- tab$G; out$PointStatus <- tab$GStatus
  out$Lower <- interval$lower; out$Upper <- interval$upper
  out$RawRatioLower <- exp(interval$log_bounds[1L])
  out$RawRatioUpper <- exp(interval$log_bounds[2L])
  out$Level <- level; out$Method <- "Exact F under normal effects"
  out$Status <- interval$projection_status
  rownames(out) <- NULL
  attr(out, "weights") <- p$weights
  attr(out, "assumption") <- assumption
  attr(out, "source_design") <- g$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "pivot") <- list(log_ratio = interval$log_ratio,
    log_bounds = interval$log_bounds, df = one$df[c(1L, 3L)])
  attr(out, "notes") <- c(
    "Exact F coverage requires the balanced independent normal-effects model and positive residual variance. It is not distribution robust.",
    "Intervals describe population G for each future count, using uncertainty and degrees of freedom from the source study. They are not Person-rank, pass/fail, Phi or plan-difference intervals.",
    "All counts for this score use the same ratio interval. Do not treat rows as independent or infer a difference from interval overlap.",
    "Intersect the mean-square expectation-ratio interval with R >= 1 before transforming to G. An empty intersection is not zero uncertainty. RawRatioLower/Upper are ratio endpoints, not G endpoints.",
    "PointStatus and Status describe separate outputs. A positive-width ratio interval that rounds to identical G endpoints is unavailable on that scale; retained log-ratio bounds may still resolve plan differences. If the base log-ratio width itself rounds away, neither interval is available. SaveRDS() retains fixed weights, assumptions and source context.")
  if (!is.null(reference)) {
    counts <- one$future
    bounds <- .mfrm_mvgt_g_difference_range(interval, one$observed,
      counts, counts[reference])
    out$Reference <- as.integer(reference)
    out$ReferenceCount <- counts[reference]
    out$ReferenceValue <- tab$G[reference]; out$Value <- tab$G
    out$Metric <- "G difference"; out$Estimate <- out$Value - out$ReferenceValue
    out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
    out$Lower <- bounds$lower; out$Upper <- bounds$upper; out$Status <- bounds$status
    out$Method <- "F confidence-set image under normal effects"
    out <- out[-reference, , drop = FALSE]
    rownames(out) <- NULL
    attr(out,"reference") <- list(scenario=as.integer(reference),count=counts[reference],
      value=tab$G[reference],point_status=tab$GStatus[reference])
    attr(out,"notes") <- c(
      "G differences map the entire nonnegative F confidence set, including any interior extremum. They do not subtract marginal endpoints or use a delta approximation.",
      "Coverage is at least the requested level under balanced independent normal effects with positive residual variance; it can be conservative. The same event protects the returned comparisons jointly for this one fixed score/composite.",
      "Positive differences favor the comparison plan; the sign from adding levels follows from the model. The interval quantifies size, not cost-effectiveness, pass/fail accuracy or changes in Phi.",
      "RawRatioLower/Upper retain the base ratio interval, not difference endpoints. PointStatus and Status describe separate outputs. An empty base set or unresolved base log-ratio width remains unavailable, including for repeated counts.",
      "SaveRDS() retains the reference, fixed weights, assumptions and source context. Existing D-study plots do not display these difference intervals.")
  }
  .mfrm_mvgt_interval_result(out, d, p, method)
}

# For unchanged nested child counts, averaging within Person-by-parent cells
# gives the one-facet law. G/Phi are scale invariant: no division of the
# three mean squares by the common child count is needed for their pivots.
.mfrm_mvgt_one_facet_target <- function(d, p) {
  g <- d$gstudy
  if (is.null(g$design$nesting)) return(list(ms = p$ms, df = p$design$df,
    observed = unname(g$design$counts[2L]), future = d$design_grid[[1L]]))
  parent <- unname(g$design$nesting)
  list(ms = p$ms[c(1L, 2L, 4L)], df = p$design$df[c(1L, 2L, 4L)],
    observed = unname(g$design$counts[[parent]]),
    future = d$design_grid[[g$design$count_columns[[parent]]]])
}

.mfrm_mvgt_interval_result <- function(out, d, p, method) {
  out$MethodCode <- rep(method, length.out = nrow(out))
  out$Available <- is.finite(out$Lower) & is.finite(out$Upper)
  out$NegativeComponent <- any(p$estimate < 0)
  if (is.null(out$Validation)) out$Validation <- ifelse(out$MethodCode == "approximate_f",
    "Experimental approximate F; finite-sample undercoverage documented",
    ifelse(out$MethodCode == "mean_square_region",
      "Experimental numerical projection; conservative joint coverage under the normal model",
      ifelse(out$MethodCode == "gci",
        "Generalized interval; frequentist coverage not exact", "Normal-model F confidence set")))
  attr(out, "design_grid") <- d$design_grid
  # Retain the sufficient summaries so curve overlays cannot mix source studies.
  attr(out, "source_summary") <- d$gstudy[c("design", "components", "mean_products", "degrees_of_freedom")]
  nesting <- d$gstudy$design$nesting
  if (!is.null(nesting)) {
    g <- d$gstudy
    child <- d$design_grid[[g$design$count_columns[[names(nesting)]]]]
    parent <- d$design_grid[[g$design$count_columns[[unname(nesting)]]]]
    attr(out, "source_design") <- g$design[c("facets", "counts", "nesting", "child_counts",
      "method", "balanced", "complete", "score_convention")]
    attr(out, "nesting_context") <- list(child = names(nesting), parent = unname(nesting),
      source_children_per_parent = if (length(unique(g$design$child_counts)) == 1L)
        unname(g$design$child_counts[1L]) else g$design$child_counts,
      future = data.frame(Scenario = seq_along(child), ChildrenPerParent = child, Parents = parent))
    if (!is.null(attr(out, "reference"))) {
      ref <- attr(out, "reference")
      if (!is.null(ref$count)) ref$count_facet <- unname(nesting)
      ref$children_per_parent <- child[ref$scenario]
      attr(out, "reference") <- ref
    }
    note <- if (all(method == "reml_profile") && !is.null(attr(out, "base_interval")))
      "Nested REML comparisons keep a common future child count per parent and vary the parent count. The source may have different or unequal child counts; its full observed assignment remains in the likelihood." else if (all(method == "reml_profile"))
      "Nested counts mean children per parent and numbers of parents; both may change. REML uses the five grouping kernels of the observed nested assignment, including interactions, without requiring balanced mean-square pivots." else if (all(method == "mls"))
      "Nested counts mean children per parent and numbers of parents. For SEM comparisons, every plan must keep the same future child count; the common base SEM then scales with the inverse square root of the parent count. Exact versus MLS coverage follows the reported Method and Validation, including when the future child count differs from the source." else if (all(method %in% c("exact_f", "gci")))
      "The nested child count per parent equals the source count in every plan. Averaging within each Person-by-parent cell gives the one-facet sampling law with effective parent variance v_parent + v_child/R and error variance v_P:parent + v_E/R. This retains both random facets; it does not make a facet fixed. The original source and child counts are retained. ReferenceCount, when present, counts parents."
      else if (all(method == "delta"))
        "Nested counts mean children per parent and numbers of parents; both may change. The delta projection retains the five original ANOVA/MINQUE components and their sampling covariance; no unchanged-child reduction is used."
      else "Nested counts mean children per parent and numbers of parents; both may change. The five-component mean-square map retains the full random nested design. Explicit experimental methods do not use the unchanged-child one-facet reduction."
    attr(out, "notes") <- c(note, attr(out, "notes"))
  }
  class(out) <- c("mfrm_dstudy_intervals", "data.frame")
  out
}

.mfrm_mvgt_g_difference_range <- function(interval, observed, counts, reference) {
  x <- pmax(interval$log_bounds,0)
  # log((R-1)/k), retaining zero and infinite parameter-space endpoints.
  log_t <- x + log(-expm1(-x)) - log(observed)
  .mfrm_mvgt_count_difference_range(log_t, counts, reference, interval$status)
}

.mfrm_mvgt_count_difference_range <- function(log_t, counts, reference, base_status) {
  lower <- upper <- rep(NA_real_,length(counts))
  status <- rep(base_status,length(counts))
  if (base_status %in% c("Available","Includes zero","Upper endpoint at zero")) {
    softplus <- function(x) pmax(x,0) + log1p(exp(-abs(x)))
    for (i in seq_along(counts)) {
      b <- counts[i]; a <- reference
      values <- c(0,0)
      finite <- is.finite(log_t)
      if (a != b) values[finite] <- sign(b-a)*exp(log(abs(b-a))+log_t[finite] -
        softplus(log(a)+log_t[finite])-softplus(log(b)+log_t[finite]))
      lost <- a != b && any(values[finite] == 0)
      stationary <- -(log(a)+log(b))/2
      if (stationary >= log_t[1L] && stationary <= log_t[2L])
        values <- c(values,sign(b-a)*tanh(log1p(abs(b-a)/min(a,b))/4))
      lower[i] <- min(values); upper[i] <- max(values)
      status[i] <- if (a == b) "Identical counts" else
        if (base_status == "Upper endpoint at zero") base_status else
          if (lower[i] <= 0 && upper[i] >= 0) "Includes zero" else "Available"
      if (lost || (a != b && log_t[1L] < log_t[2L] && lower[i] == upper[i])) {
        lower[i] <- upper[i] <- NA_real_
        status[i] <- "Difference interval below numeric resolution"
      }
    }
  }
  list(lower=lower,upper=upper,status=status)
}

.mfrm_mvgt_phi_table <- function(d, tab, p, level, draws, seed, reference = NULL,
                               one = .mfrm_mvgt_one_facet_target(d, p)) {
  integer_scalar <- function(x) is.numeric(x) && !is.object(x) && !is.complex(x) &&
    is.null(dim(x)) && length(x) == 1L && is.finite(x) && x == floor(x)
  tail <- (1-level)/2
  if (!integer_scalar(draws) || draws < 1000 || draws > .Machine$integer.max || draws*tail < 100)
    stop("`draws` must be an integer >= 1000 with at least 100 expected draws in each tail.", call. = FALSE)
  if (!integer_scalar(seed) || seed < 0 || seed > .Machine$integer.max)
    stop("`seed` must be a nonnegative integer.", call. = FALSE)
  n <- unname(d$gstudy$design$counts[1L]); k <- one$observed
  m <- one$future
  transform <- function(rho) rho / (rho + (1-rho)/m)
  bounds <- rep(NA_real_, 2L); mc_bounds <- rep(NA_real_, 4L)
  zero <- NA_real_; status <- "Degenerate mean squares"
  if (all(one$ms > 0)) {
    z <- with_preserved_rng_seed(seed, .mfrm_mvgt_phi_pivots(one$ms, n, k, draws))
    status <- "Numerical interval failure"
    if (all(is.finite(z) & z >= 0 & z <= 1)) {
      # Reuse one ordered sample for both endpoint bands and all future counts.
      z <- sort(z)
      bounds <- unname(stats::quantile(z, c(tail,1-tail), type=7))
      # Differences depend on both base endpoints: protect them jointly.
      band_tail <- if (is.null(reference)) .025 else .0125
      rank <- as.integer(c(stats::qbinom(band_tail,draws,tail),
        stats::qbinom(1-band_tail,draws,tail)+1,
        stats::qbinom(band_tail,draws,1-tail),stats::qbinom(1-band_tail,draws,1-tail)+1))
      mc_bounds <- c(0,z,1)[rank+1L]
      zero <- mean(z == 0)
      status <- if (bounds[2L] == 0) "Upper endpoint at zero" else
        if (bounds[1L] == bounds[2L]) "Pivotal interval below numeric resolution" else
          if (bounds[1L] == 0) "Includes zero" else "Available"
    }
  }
  out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
  out$Metric <- "Phi"; out$Estimate <- tab$Phi; out$PointStatus <- tab$PhiStatus
  out$Lower <- transform(bounds[1L]); out$Upper <- transform(bounds[2L])
  out$RawRatioLower <- out$RawRatioUpper <- NA_real_
  out$Level <- level; out$Method <- "Generalized CI under normal effects"; out$Status <- status
  # A positive-width base interval can collapse after future-count scaling.
  collapsed <- is.finite(out$Lower) & is.finite(out$Upper) & out$Lower == out$Upper &
    (bounds[2L] > bounds[1L] | bounds[2L] > 0)
  out$Lower[collapsed] <- out$Upper[collapsed] <- NA_real_
  out$Status[collapsed] <- "Phi interval below numeric resolution"
  rownames(out) <- NULL
  attr(out,"weights") <- p$weights; attr(out,"assumption") <- "normal"
  attr(out,"source_design") <- d$gstudy$design[c("facets","counts","method","balanced","complete","score_convention")]
  attr(out,"source_rows") <- d$gstudy$data_usage$counts
  attr(out,"pivot") <- list(single_rating_bounds=bounds,df=one$df)
  attr(out,"monte_carlo") <- list(draws=draws,
    draws_used=if (all(one$ms > 0)) draws else 0L,seed=seed,rng_kind=RNGkind(),
    quantile_type=7L,zero_proportion=zero,band_level=.95,
    endpoint_bands=data.frame(Scenario=out$Scenario,
      LowerMCLower=transform(mc_bounds[1L]),LowerMCUpper=transform(mc_bounds[2L]),
      UpperMCLower=transform(mc_bounds[3L]),UpperMCUpper=transform(mc_bounds[4L])))
  attr(out, "monte_carlo")$endpoint_bands[collapsed, -1L] <- NA_real_
  attr(out,"notes") <- c(
    "Phi includes the random facet main effect in absolute error. Independent normal effects and a complete balanced one-facet source are required.",
    "This generalized confidence interval does not guarantee exact frequentist coverage. More pivotal draws reduce numerical error, not model uncertainty or coverage error.",
    "The GCI truncates its pivotal numerator at zero, not the source variance estimates. PointStatus and Status describe separate outputs; a zero endpoint is not proof of zero uncertainty.",
    "The monte_carlo attribute contains separate 95% numerical bands for the endpoints. Increase draws if these bands affect the intended decision; save the whole result with saveRDS().",
    "All future counts transform the same single-rating interval and use source degrees of freedom. Rows are dependent; overlap does not test a plan difference. These are not pass/fail accuracy intervals.",
    "Widths lost to floating-point resolution leave endpoints and their numerical bands unavailable. A base pivotal interval with lost resolution cannot supply difference bounds, including for identical plans. Diagnostic single-rating bounds remain in pivot; the sampled zero atom retains its separate status.")
  if (!is.null(reference)) {
    image <- function(rho, state = status) .mfrm_mvgt_count_difference_range(
      stats::qlogis(rho), m, m[reference], state)
    mapped <- image(bounds)
    outer <- image(mc_bounds[c(1L,4L)])
    # A common inner set bounds both image endpoints. If it is empty, the
    # outer image still bounds each endpoint without inventing a common set.
    inner <- if (all(is.finite(mc_bounds)) && mc_bounds[2L] <= mc_bounds[3L])
      image(mc_bounds[c(2L,3L)]) else list(lower=outer$upper,upper=outer$lower)
    out$Reference <- as.integer(reference); out$ReferenceCount <- m[reference]
    out$ReferenceValue <- tab$Phi[reference]; out$Value <- tab$Phi
    out$Metric <- "Phi difference"; out$Estimate <- out$Value-out$ReferenceValue
    out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
    out$Lower <- mapped$lower; out$Upper <- mapped$upper; out$Status <- mapped$status
    out$Method <- "Generalized confidence-set image under normal effects"
    out <- out[-reference, , drop=FALSE]; rownames(out) <- NULL
    attr(out,"reference") <- list(scenario=as.integer(reference),count=m[reference],
      value=tab$Phi[reference],point_status=tab$PhiStatus[reference])
    mc <- attr(out,"monte_carlo")
    mc$base_status <- status
    mc$base_endpoint_bands <- matrix(mc_bounds, nrow=2L, byrow=TRUE,
      dimnames=list(c("Lower","Upper"),c("MCLower","MCUpper")))
    mc$band_scope <- "At least 95% jointly over both numerical endpoints and returned comparisons; not coefficient coverage"
    mc$endpoint_bands <- data.frame(Scenario=out$Scenario,
      LowerMCLower=outer$lower[-reference],LowerMCUpper=inner$lower[-reference],
      UpperMCLower=inner$upper[-reference],UpperMCUpper=outer$upper[-reference])
    attr(out,"monte_carlo") <- mc
    attr(out,"notes") <- c(
      "Phi differences map the whole single-rating generalized interval, including any interior extremum. They do not subtract marginal endpoints or take quantiles of transformed differences.",
      "Whenever the base generalized set contains the true ratio, all mapped comparisons contain their true differences. Approximate frequentist coverage of the base method is retained; there is no exact nominal-coverage guarantee.",
      "Positive differences favor the comparison plan. Adding levels already determines the direction under this model; the interval describes size, not cost-effectiveness or pass/fail accuracy.",
      "The monte_carlo attribute gives at least 95% joint numerical bands from Bonferroni-protected base endpoints and their full set images. These describe quantile error, not sampling coverage. More draws do not add source information.",
      "PointStatus and Status remain separate. The base pivot bounds, zero proportion and reference context are retained; a sampled upper endpoint at zero is not proof of zero uncertainty. Base widths lost to numeric resolution leave all difference intervals unavailable. Save the full result with saveRDS().")
  }
  out
}

.mfrm_mvgt_phi_pivots <- function(ms, n, k, draws) {
  facet <- (k-1)*ms[2L]/stats::rchisq(draws,k-1)
  person <- (n-1)*ms[1L]/stats::rchisq(draws,n-1)
  residual <- (n-1)*(k-1)*ms[3L]/stats::rchisq(draws,(n-1)*(k-1))
  # Mondal et al. (2025), Eq. 11, simplified without changing truncation.
  pmax(0,person-residual)/(person+(k/n)*facet+(k-1-k/n)*residual)
}

.mfrm_mvgt_g_f_interval <- function(person_ms, residual_ms, df, observed, future, level) {
  lower <- upper <- rep(NA_real_, length(future))
  log_bounds <- rep(NA_real_, 2L); log_ratio <- NA_real_
  status <- "Zero residual mean square"
  if (residual_ms > 0) {
    log_ratio <- log(person_ms) - log(residual_ms)
    tail <- (1 - level) / 2
    log_bounds <- log_ratio - log(c(stats::qf(tail, df[1L], df[2L], lower.tail = FALSE),
      stats::qf(tail, df[1L], df[2L])))
    status <- "Numerical interval failure"
    if (!anyNA(log_bounds) && log_bounds[1L] <= log_bounds[2L]) {
      status <- "Empty parameter-space interval"
      if (is.finite(log_bounds[2L]) && log_bounds[2L] >= 0 &&
          log_bounds[1L] == log_bounds[2L]) {
        status <- "F ratio interval below numeric resolution"
      } else if (log_bounds[2L] >= 0) {
        bounded <- pmax(log_bounds, 0)
        transform <- function(x) {
          u <- -expm1(-x)
          u / (u + (observed / future) * exp(-x))
        }
        lower <- transform(bounded[1L]); upper <- transform(bounded[2L])
        status <- if (bounded[1L] == 0) "Includes zero" else "Available"
      }
    }
  }
  projection_status <- rep(status, length(future))
  collapsed <- is.finite(lower) & is.finite(upper) & lower == upper &
    log_bounds[2L] > max(0, log_bounds[1L])
  lower[collapsed] <- upper[collapsed] <- NA_real_
  projection_status[collapsed] <- "G interval below numeric resolution"
  # Keep the base ratio status: its log bounds can still identify differences
  # even when individual G endpoints both round to one.
  list(lower = lower, upper = upper, log_ratio = log_ratio,
    log_bounds = log_bounds, status = status, projection_status = projection_status)
}
