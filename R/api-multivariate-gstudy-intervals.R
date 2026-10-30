#' Sampling uncertainty in G-study variance and covariance components
#'
#' Add standard errors and confidence intervals to the variance components
#' of one observed score, one prespecified weighted composite, or the covariance
#' between two original scores. The default uses a complete balanced ANOVA
#' G-study without refitting. Explicit experimental REML profiles
#' also accept identifiable incomplete sources, fitting a separate normal model.
#' Explicit residual contrasts provide residual-only intervals from retained
#' used ratings, including incomplete sources, without likelihood fitting.
#'
#' @param object A result from [mfrm_g_study()] fitted with
#'   `method = "anova"`. One or two crossed random facets, or the supported
#'   balanced nested design, are allowed. With `method = "reml_profile"`,
#'   variance and two-score covariance intervals also accept incomplete or
#'   unequal MINQUE(0) sources with their raw used ratings retained.
#'   `method = "residual_contrast"` accepts these sources for Residual only.
#' @param parm Component names or integer positions to return; `NULL` selects
#'   all components. Names are those in `object$components`.
#' @param level Pointwise confidence level strictly between zero and one.
#' @param score One original score name for variances (the default is the first
#'   score), or exactly two distinct original names for their covariance.
#'   With two names, explicitly select a covariance `method` below.
#' @param weights Optional named numeric vector specifying one composite.
#'   Name every score once, with finite weights and at least one nonzero
#'   value. Weights are used as supplied, including signed weights for a
#'   difference score. Do not supply both `score` and `weights`.
#' @param assumption Required: `"normal"` asserts independent normal random
#'   effects with constant component variances for the selected score or
#'   composite. Persons and the selected measurement facets are random.
#'   The function does not test these assumptions.
#' @param method `"auto"` or `"mls"` for a single score/composite from balanced
#'   ANOVA. Explicit `"reml_profile"` computes approximate component intervals
#'   from raw used ratings, for balanced or identifiable incomplete sources;
#'   it requires optional lme4 and minqa. It retains original points and reports
#'   a separate REML center, but does not supply a sampling SE. Its chi-square
#'   calibration is experimental and it is never chosen automatically. For two-score
#'   covariance, explicitly choose `"gci"` (generalized interval with approximate
#'   coverage), `"regression_region"` (conservative normal-model enclosure), or
#'   `"reml_profile"` (experimental joint normal-model profile). The first two
#'   require balanced ANOVA. No covariance method is silently selected or
#'   substituted after failure.
#'   `"residual_contrast"` supplies only the Residual interval: exact chi-square
#'   for one score/composite, or a conservative enclosure for two-score
#'   covariance, under the normal residual model. It requires raw used ratings,
#'   permits incomplete sources, reports a separate `IntervalCenter`, and does
#'   not supply the original estimator's SE. Other components retain their
#'   estimates with unavailable intervals. See *Residual contrasts* below.
#' @param draws Number of pivotal draws for covariance `method = "gci"`;
#'   default 20000. At least 1000 and at least 100 expected draws in each tail
#'   are required. Inspect the numerical endpoint bands, not just this minimum.
#' @param seed Nonnegative integer seed for covariance GCI; the caller's random
#'   state is restored. `draws` and `seed` are rejected for other methods.
#' @param ... Reserved; additional arguments are rejected.
#'
#' @details
#' For a source-by-target choice table, see *Choose the G-theory target and
#' uncertainty* in [mfrmr_methods_reference]. Component intervals, G/Phi
#' intervals, measurement SEM and paired plan differences have distinct targets.
#'
#' A component estimates how much scores vary because of persons, conditions,
#' or their interactions. Its SE describes sampling uncertainty in that
#' estimate. Both the component and its SE have squared score units. This SE
#' is different from a D-study's measurement SEM, which is in score units.
#'
#' Default balanced intervals use the modified large-sample (MLS) method of Ting et al. (1990)
#' for signed linear combinations of independent ANOVA mean squares, as
#' detailed by Hess and Iyer (2002, pp. 2--3). The residual component, which
#' is a single mean square, has the usual chi-square interval. Other
#' component intervals are approximate, asymmetric and pointwise, not
#' simultaneous. Normal effects, a balanced design and prespecified score
#' weights matter: the method is not robust to arbitrary score distributions
#' or weights chosen after inspecting the data. A numeric ordinal score
#' does not by itself justify the normal model.
#' At very low requested confidence levels, squared chi-square corrections
#' can lose their usual direction; a small `level` does not guarantee narrow
#' or well-calibrated MLS limits.
#'
#' Negative component estimates are retained. `RawLower` and `RawUpper` are
#' the unrestricted MLS endpoints. `Lower` and `Upper` intersect that interval
#' with the nonnegative variance parameter space. A negative lower endpoint
#' therefore becomes zero, but the point estimate is not changed. If the
#' entire raw interval is below zero, the intersection is empty: both bounds
#' are `NA` with `Status = "Empty nonnegative interval"`, while the raw
#' endpoints remain visible. It is not replaced by an interval of zero width.
#' Inspect the design, assumptions and raw estimates in that case.
#'
#' For MLS, `Status = "Includes zero"` denotes a nonempty interval containing zero;
#' `"Available"` denotes a strictly positive interval. If all mean squares
#' contributing to a component are zero, its estimated SE is zero but bounds
#' are unavailable (`"Degenerate mean squares"`); this is not evidence of
#' precise population knowledge. `"Numerical interval failure"` retains the
#' estimate and SE when the MLS width cannot be evaluated reliably.
#' A positive interval whose endpoints coincide after rounding instead has
#' `Status = "Variance interval below numeric resolution"` (or the corresponding
#' component-scale status), with missing `Lower`/`Upper`. This can occur at an
#' extremely small requested confidence level. Positive endpoints that overflow
#' or underflow also remain unavailable. Raw endpoints are retained for diagnosis;
#' equal raw endpoints in such a row do not establish zero uncertainty.
#'
#' This method estimates variances of the selected score, including all
#' covariances needed for a fixed composite. The separate two-score covariance
#' calculation is described below. Neither supplies G/Phi or plan-difference
#' intervals. Explicit REML profiles and residual contrasts accept incomplete
#' sources, with the different scopes described below.
#' [mfrm_multivariate_d_compare()] provides a separate
#' approximation for certain paired plan differences. Use
#' `plot(g, intervals = TRUE, assumption = "normal")` to calculate and show
#' these component intervals directly, with `score` or `weights` selecting
#' the desired score. Its returned data retain the complete interval tables.
#'
#' @section Residual contrasts:
#' Use `method = "residual_contrast"` to ask how uncertain the combined
#' residual variance or covariance is. It supports the same crossed and nested
#' source designs, including identifiable incomplete MINQUE(0) sources, when
#' the observed assignment leaves residual degrees of freedom. It does not
#' supply intervals for Person, facet or lower-order interaction components
#' or G/Phi. For a one-facet D-study, [confint.mfrm_multivariate_d_study()]
#' can map this interval to relative SEM and its plan differences; the
#' component result itself remains in squared score units. With one rating per cell, the residual
#' combines the highest interaction and observation error; it cannot separate them.
#'
#' Let `Z` contain all nonresidual incidence columns, including interactions
#' and the mean. With `nu = n - rank(Z) > 0`, choose an orthonormal basis `H`
#' satisfying `H'Z = 0`. Under independent homoskedastic normal residuals,
#' `H'Y = H'E` and `Y'HH'Y` has a Wishart distribution with `nu` degrees of
#' freedom and residual covariance `Sigma_E`. Thus for prespecified weights
#' `w`, the residual sum of squares divided by `w'Sigma_E w` is chi-square
#' with `nu` degrees of freedom. Its quantiles give exact pointwise limits
#' for a positive residual variance. The calculation uses QR residuals,
#' checks rank against singular values and verifies orthogonality, without
#' forming an observation-by-observation projection matrix. This applies the
#' standard multivariate regression residual law to the G-study incidence model.
#'
#' For two scores, the same Wishart block feeds the conservative regression
#' enclosure described in *Between-score covariance elements* below. One-degree-of-freedom
#' or singular pairs instead use marginal variance bounds and Cauchy-Schwarz.
#' Signed covariance bounds are retained. Coverage is conservative under the
#' model, not exact equal-tail coverage, and the enclosure can be very wide.
#' A sparse crossed arrangement can leave few residual contrasts despite
#' many ratings. Inspect `attr(ci, "residual_projection")$df` and the limits.
#'
#' The assignment is treated as fixed; omissions must preserve the residual
#' model. This does not correct selective missingness, heteroskedasticity or
#' dependence between remaining residuals. Normal residuals are essential.
#' Claims are pointwise for prespecified scores/weights, not simultaneous or
#' adjusted for selecting a result after looking at the data.
#'
#' `Estimate` retains the original ANOVA/MINQUE component. `IntervalCenter`
#' is the separate residual-contrast estimate; the two agree in balanced
#' ANOVA but can differ for incomplete sources. The original point can lie
#' outside the interval. `SE = NA` does not mean zero
#' uncertainty: this method does not calculate the original estimator's SE.
#' The scalar chi-square calculation shares the numerical-width checks above;
#' an unresolved interval retains its original point and projected center.
#' Use [vcov.mfrm_multivariate_gstudy()] for its separately qualified sampling
#' covariance. Both the component and these limits have squared score units.
#' Zero residual degrees of freedom, numerically unresolved residual variation
#' or failed projection checks give unavailable limits with explicit `Status`.
#' Unsupported component rows remain visible. No unavailable row is a
#' zero-width interval. Raw ratings and stored components are checked for
#' agreement; recreating a modified G-study is preferable to editing its fields.
#'
#' @section REML component profiles:
#' For example, uncertainty in the Task variance asks how strongly randomly
#' sampled tasks differ in their average scores. `method = "reml_profile"`
#' fits a separate normal random-effects model to the observed assignment:
#' one common mean and independent Person, facet, interaction and residual
#' effects. Omissions must be noninformative conditional on that assignment.
#' It does not correct selective missingness or replace the original MINQUE
#' estimator. `Estimate` retains the original component, including negative
#' values; `IntervalCenter` is the REML estimate of the same population variance.
#' The original point can be outside the interval. Neither component matrices
#' nor original estimates are clipped or repaired.
#'
#' Write `theta = s2 * v`, with residual ratio one. For target component `j`,
#' constrain `tau = theta[j]` by setting `s2 = tau / v[j]`. Restore residual
#' scale to the profiled restricted criterion using
#' `D(v,s2) = D_profile(v) + (n-1)*(log(s2/s2_hat(v)) + s2_hat(v)/s2 - 1)`.
#' Optimize the remaining ratios. At a zero nonresidual target, set that
#' component exactly to zero and profile the others; the residual remains
#' strictly positive. The interval uses `D - D_min <= qchisq(level, 1)`.
#' See the REML section of [confint.mfrm_multivariate_d_study()] for the
#' underlying normal model, numerical checks and relationship to ordinary ML
#' profile intervals. Bates et al. supply the REML criterion, not a finite-sample
#' validation of this particular G-theory procedure.
#'
#' Search outward from the REML center and refine the first cutoff crossing
#' in each direction. This returns a numerical interval around that center;
#' it does not enumerate remote disconnected accepted regions. Multiple starts,
#' a second optimizer, and finite-difference checks including zero boundaries
#' do not certify global optima or rule out narrower undetected crossings.
#' Unresolved required evaluations, roots or tails leave accepted bounds
#' unavailable, retaining any resolved `CandidateLower/Upper` separately.
#' Reaching a search limit does not justify an infinite bound. A constant
#' selected score is degenerate, not evidence of zero-width uncertainty.
#' Weighted scores are scaled together with their weights before fitting;
#' reciprocal changes of score units and weights preserve the selected composite.
#' At very small confidence levels, the likelihood-ratio cutoff may be too
#' small to distinguish from numerical error. Such intervals have missing
#' bounds and an explicit `Status`, retaining the original estimate and
#' representable REML center. For each component, `diagnostics[[i]]$cutoff`
#' records the cutoff, endpoint tolerance and deviance resolution. Nonzero
#' bounds lost when restoring score units are also unavailable; normalized
#' profile diagnostics remain available. These checks also apply to joint
#' covariance profiles. They are numerical restrictions, not minimum sample
#' sizes or evidence of zero uncertainty.
#'
#' `Available` describes numerical completion, not statistical coverage.
#' A limited check reused 40 samples per condition: component true-value
#' acceptance was 92.5--100% with positive components and 90--100% for checked
#' zero components, counting numerical failures as misses. The Monte Carlo SE
#' can reach 4.7 percentage points. These are coarse feasibility results, not
#' coverage of full endpoint searches or support for a generally recommended
#' method. The chi-square cutoff is asymptotic; a universal boundary-mixture
#' correction is not assumed. Intervals are pointwise for prespecified scores
#' and components, neither simultaneous nor distribution robust.
#'
#' `SE = NA` and `SEStatus = "Not supplied by profile inversion"` mean that
#' this method did not calculate a sampling SE. The interval width must not
#' be converted into one by dividing by a normal critical value. `RawLower`
#' and `RawUpper` are inapplicable MLS fields and remain `NA`; they are not
#' failed or unrestricted profile endpoints. All reported estimates and bounds
#' use squared score units. Scalar component intervals cannot be added to
#' obtain SEM, G/Phi or covariance-element intervals.
#' Use [vcov.mfrm_multivariate_gstudy()] separately for sampling SEs of the
#' original ANOVA/MINQUE component estimates. Those SEs use normal product
#' moments and do not describe the separate REML estimator or calibrate its interval.
#'
#' Use `plot(g, intervals = TRUE, method = "reml_profile", assumption = "normal")`
#' to compute the intervals and draw the original points in one call. Its
#' returned data retain the full interval table, including REML centers and
#' diagnostics. Save that plot object with `saveRDS()` and use `as_ggplot()`
#' after reloading to render it without another fit.
#'
#' @section Calculation of the SE:
#' Let `MS` be the independent ANOVA mean squares, `df` their degrees of
#' freedom, and `A` the linear map from mean squares to raw components.
#' Under the normal model, `df[i] * MS[i] / E(MS[i])` is chi-square.
#' The sampling covariance is estimated by
#' `A %*% diag(2 * MS^2 / (df + 2)) %*% t(A)`; `SE` is the square root
#' of its diagonal. The squared SE is unbiased for the sampling variance
#' under this model; the SE itself need not be unbiased. This follows from
#' `E(MS[i]^2) = E(MS[i])^2 * (df[i] + 2) / df[i]` and matches the
#' urGENOVA SEVC convention. The existing paired D-study comparison instead
#' uses the plug-in denominator `df`. Its method is unchanged. MLS intervals
#' are calculated from the mean squares and chi-square/F quantiles, not by
#' adding and subtracting a normal critical value times this SE.
#'
#' @section Between-score covariance elements:
#' For example, `score = c("V", "W")` estimates how the V and W scores vary
#' together for each source. A negative covariance can be meaningful; its
#' interval is signed and is never intersected with nonnegative variances.
#' Covariance units are the product of the two score units. These are not
#' correlations, and large covariance alone does not establish a halo effect.
#' Select the pair before examining its results. Raw negative estimates and
#' non-positive-semidefinite component matrices are retained without repair.
#'
#' For each independent mean-product matrix `S` with df `nu`, `nu*S` has a
#' Wishart distribution under the balanced normal model. Write its entries
#' as `x = S[1,1]`, `c = S[1,2]`, `y = S[2,2]`. For `nu > 1`, an unbiased
#' estimator of the sampling variance of `c` is
#' `(nu*x*y + (nu-2)*c^2)/((nu-1)*(nu+2))`. Sum these estimates with squared
#' ANOVA-map coefficients to obtain the component's squared SE. Its square
#' root need not be unbiased. If a contributing block has df one, this SE
#' formula is unavailable; `SEStatus` says so independently of interval status.
#' Computation scales the scores separately, so a change in their relative
#' units does not erase the smaller marginal variance. The `mean_products`
#' attribute retains normalized entries `(x, c, y)`, the two `score_scale`
#' divisors and their product `scale`; multiply the normalized covariance
#' column by `scale` to restore covariance units (use the respective squared
#' divisors for the two variance columns).
#' A squared SE can be outside floating-point range even when its square root
#' is representable; the calculation retains that representable SE.
#' Unrepresentable endpoints or a positive interval width lost to rounding
#' produce missing limits with an explicit `Status`, retaining the original
#' estimate and any calculable SE. This also applies to residual-contrast
#' covariance intervals, which retain their separate `IntervalCenter`.
#'
#' **Generalized interval.** Krishnamoorthy and Xia (2007, Section 2) construct
#' a covariance-matrix pivot `C %*% solve(V)`, where `C` is the lower Cholesky
#' factor of the observed sum-of-products matrix and `V` has independent
#' Bartlett entries. Its covariance element simplifies to
#' `nu * (c - sqrt(x*y-c^2)*Z/sqrt(U2)) / U1`, with independent
#' `U1 ~ chi-square(nu)`, `U2 ~ chi-square(nu-1)` and `Z ~ N(0,1)`.
#' Our extension generates independent pivots for the ANOVA blocks, applies
#' the signed component map, and takes type-1 empirical quantiles. It does
#' not filter negative draws or impose component-matrix constraints. Each
#' contributing block must have df greater than one and a nonsingular pair.
#' The paper's favorable results for a simple correlation do not establish
#' exact frequentist coverage for these G-theory covariance components.
#'
#' The `monte_carlo` attribute gives 95% joint numerical bands for the two
#' endpoints of each component interval, using binomial order-statistic bounds.
#' These bands concern simulation error in the quantiles. They do not measure
#' source-study uncertainty or establish sampling coverage. More draws reduce
#' numerical error; they do not resolve insufficient facet levels.
#' At a very small requested confidence level, the two quantile probabilities
#' or empirical ranks may coincide. In that case the interval is unavailable
#' and `Status` distinguishes numerical probability resolution from an
#' unresolved interval at the current draw count. Increasing `draws` may
#' resolve ranks but cannot restore distinct probabilities lost to rounding.
#'
#' **Conservative enclosure.** For each component with `K` contributing
#' independent mean-product blocks, give each block coverage `gamma = level^(1/K)`.
#' When `nu > 1` and the pair is nonsingular, the regression-slope t pivot
#' (df `nu-1`) is independent of the predictor sum-of-squares chi-square pivot
#' (df `nu`). Give each pivot coverage `sqrt(gamma)`. The covariance is the
#' slope times the predictor variance; enclose their product over that rectangle.
#' In terms of `S`, take all products of
#' `c +/- t_quantile*sqrt((x*y-c^2)/(nu-1))` and the two ratios
#' `nu / chi_square_quantile`, with equal tails `(1-sqrt(gamma))/2`.
#' This expression is symmetric in the scores and respects changes of units.
#'
#' For a one-df or singular block, use two one-sided marginal variance upper
#' limits, each with tail `(1-gamma)/2`, and Cauchy-Schwarz to bound covariance
#' between plus/minus the square root of their product. Bonferroni protects
#' this block. Propagate block intervals through the signed ANOVA map. The
#' independent block events give at least `level` coverage for the component
#' under positive marginal variances and normal effects. This is our conservative
#' construction from standard pivots, not a published G-theory formula attributed
#' to Krishnamoorthy and Xia. It is pointwise, not simultaneous over components
#' or score pairs. It can be extremely wide: enclosing nuisance regions adds
#' conservatism beyond the source design's limited information.
#'
#' A limited normal-model check covered six small/larger one-facet, crossed
#' and nested sources. A 200-replication GCI screen was followed by four
#' exploratory checks using 2000 saved samples and 8000 pivotal draws each.
#' Their nominal 95% coverages were 94.65%, 94.80%, 99.65% and 96.50%
#' (Monte Carlo SE 0.13--0.50 percentage points). These concern respectively
#' one-facet residual (10 persons, 2 tasks), crossed residual and Task
#' (10 persons, 2 raters, 3 tasks), and crossed Rater (60,4,8) covariances.
#' The small Task interval remained very wide; covariance GCIs involving a
#' one-df block were unavailable. This check does not establish universal
#' coverage or qualify a generally preferred covariance method. Neither method
#' is distribution robust or applies to incomplete MINQUE(0) sources.
#'
#' `plot(g, type = "covariance", score = c("V", "W"), intervals = TRUE,
#' method = "gci", assumption = "normal")` displays one signed interval plot
#' across sources. Without intervals, `type = "covariance"` retains its matrix
#' heatmaps. Saved plot data retain the full interval table and numerical bands.
#'
#' @section Joint REML covariance profiles:
#' For incomplete sources, explicitly select `method = "reml_profile"` with
#' two original scores. This fits both scores jointly, with one mean per score,
#' independent bivariate normal random effects, a positive-semidefinite matrix
#' for each component and a positive-definite residual matrix. Omissions must
#' be noninformative conditional on the observed assignment. This is a different
#' inference procedure from MINQUE: the original covariance remains `Estimate`,
#' while `IntervalCenter` is the separate joint REML covariance. Both target
#' the same population quantity. The original point can lie outside the interval.
#'
#' With score-major observations, let `V = sum(Sigma[k] %x% K[k])`, where each
#' `K[k]` records shared levels, and `X = I[2] %x% rep(1,n)`.
#' The restricted deviance is `log|V| + log|X'V^-1 X| + y'Py` plus a constant,
#' with `P = V^-1 - V^-1 X (X'V^-1 X)^-1 X'V^-1`.
#' Fix one off-diagonal component element, optimize all other variances and
#' covariances, and invert `D - D_min <= qchisq(level, 1)`. No marginal interval
#' endpoints are combined. Internally, residual whitening permits sparse lme4
#' evaluation after restoring residual scale and the transformation determinant.
#' Independent direct-matrix calculations verify that criterion.
#'
#' Zero covariance is permitted exactly; it does not imply zero marginal
#' variances. Signed matrix factors preserve positive semidefiniteness during
#' optimization. Covariance-coordinate derivatives check boundary optimality
#' as well as ordinary parameter gradients. This numerical parameterization
#' does not remove the statistical boundary problem (Wu & Neale, 2013).
#' Multiple starts and these necessary local checks do not prove global optima.
#' If the initially selected candidate fails these checks, already computed
#' converged candidates within `1e-7` deviance of the smallest recorded value
#' are checked in objective order. The first candidate passing every unchanged
#' gradient, covariance-cone and target-constraint check is used. More distant
#' local solutions are not substituted. Selected indices and inspected checks
#' are retained in each fit's `selection`; no further optimization is triggered.
#'
#' The returned interval surrounds the REML center and ends at the first
#' checked cutoff crossing in each direction; it does not enumerate remote
#' accepted regions. Unresolved source/constraint fits, roots or tails leave
#' accepted bounds unavailable, preserving numerical candidates separately.
#' A search limit is not an infinite endpoint. Constant or numerically collinear
#' score pairs are degenerate, not evidence of zero uncertainty.
#'
#' This method is experimental and never selected automatically. Limited
#' balanced normal-model checks used 40 samples per condition and found
#' true-value LR inclusion of 82.5--100%, counting numerical failures as misses
#' (Monte Carlo SE up to 6.0 percentage points). The lowest result concerned
#' Rater covariance with 60 persons, four raters and eight tasks. More persons
#' alone do not increase the number of independent rater effects. These checks
#' do not establish full endpoint-search coverage; incomplete-source coverage
#' is not established. Do not interpret numerical availability as a coverage
#' guarantee or replace the existing balanced methods merely because this
#' method returns a narrower interval. No universal boundary-mixture correction
#' is assumed. The intervals are pointwise and not distribution robust.
#'
#' `SE` is `NA`: profile inversion does not supply sampling SEs for either the
#' original point or REML center. `RawLower/RawUpper` are inapplicable here;
#' `CandidateLower/Upper` retain exploratory endpoints. All reported covariance
#' estimates and bounds use products of the score units, retaining negative
#' signs. Saved diagnostics use standardized units; multiply by both entries
#' of `attr(result, "reml_source")$score_scale` to restore score-product units.
#'
#' @return A data frame with `Source`, selected `Score`, raw `Estimate`, `SE`,
#'   nonnegative-space `Lower` and `Upper`, unrestricted `RawLower` and
#'   `RawUpper`, `Level`, `Method`, and `Status`. Attributes `weights`, `kind`,
#'   `assumption`, `source_design`, `source_rows`, and `notes` retain context.
#'   Read `attr(result, "notes")` and save the full result with `saveRDS()`;
#'   a CSV alone does not retain weights and assumptions. The original
#'   G-study is unchanged. Stored mean products, degrees of freedom and design
#'   metadata are required for MLS; raw ratings are not needed for that calculation.
#'   REML profiles instead require raw used ratings and retain `PointStatus`,
#'   `IntervalCenter`, `SEStatus`, `CandidateLower/Upper`, `MethodCode`,
#'   `Available`, `Validation`, and portable `reml_source` and `diagnostics`
#'   attributes. Diagnostics are named in the selected component order, in
#'   standardized-score units; multiply by `reml_source$score_scale` twice for
#'   a variance or by the two separate scale entries for a covariance.
#'   No live fitting environment is saved. `PointStatus` is for scalar variances.
#'   Two-score results additionally include `Score1`, `Score2`, `SEStatus`,
#'   `score_pair` and method-specific `diagnostics`. Balanced GCI/region results
#'   retain normalized `mean_products` and optional `monte_carlo` attributes;
#'   their `RawLower`/`RawUpper` equal the signed `Lower`/`Upper`.
#'   REML results instead retain the separate center, candidates and profiles.
#'   Covariance intervals have no nonnegative-space intersection.
#'   Residual-contrast results retain `IntervalCenter`, `SEStatus`, `MethodCode`,
#'   `Available` and `Validation`, plus a `residual_projection` attribute with
#'   residual degrees of freedom, rank, scale and numerical checks. They require
#'   retained ratings but save no live fitting environment or projection matrix.
#' @references
#' Stanford University, STATS 305C. [Multivariate regression, Distribution results](https://web.stanford.edu/class/stats305c/lectures/Multivariate_regression.html#distribution-results).
#' Normal-model residual Wishart law; the G-study projection is described above.
#'
#' Wu, H., & Neale, M. C. (2013). On the likelihood ratio tests in bivariate
#' ACDE models. *Psychometrika, 78*, 441--463.
#' \doi{10.1007/s11336-012-9304-2}. Boundary-inference cautions, not a
#' calibration formula for these G-theory covariance-element intervals.
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software, 67*(1),
#' 1--48. \doi{10.18637/jss.v067.i01}. Restricted likelihood and sparse computation.
#'
#' Krishnamoorthy, K., & Xia, Y. (2007). Inferences on correlation coefficients:
#' One-sample, independent and correlated cases. *Journal of Statistical Planning
#' and Inference*, 137, 2362--2379. \doi{10.1016/j.jspi.2006.08.002}.
#' Section 2 supplies the Wishart matrix pivot; component propagation is the
#' extension documented above, with separate finite-sample evaluation.
#'
#' R Core Team. [Random Wishart distributed matrices](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/rWishart.html).
#' Normal sampling representation and elementwise sampling variance.
#'
#' Ting, N., Burdick, R. K., Graybill, F. A., Jeyaratnam, S., & Lu, T.-F. C.
#' (1990). Confidence intervals on linear combinations of variance components
#' that are unrestricted in sign. *Journal of Statistical Computation and
#' Simulation, 35*(3--4), 135--143. \doi{10.1080/00949659008811240}.
#'
#' Hess, A., & Iyer, H. (2002). A SAS/IML macro for computation of confidence
#' intervals for variance components of mixed models. *Proceedings of SUGI 27*,
#' Paper 252-27, pp. 2--3.
#' \url{https://support.sas.com/resources/papers/proceedings/proceedings/sugi27/p252-27.pdf}.
#'
#' Brennan, R. L. (2001). *Manual for urGENOVA, Version 2.1*. SEVC and
#' interval definitions, pp. 11--12; Appendix C numerical example, p. 26.
#' Available from [CASMA](https://education.uiowa.edu/casma/computer-programs).
#' @seealso [mfrm_g_study()], [plot.mfrm_multivariate_gstudy()],
#'   [mfrm_multivariate_d_compare()], [mfrmr_methods_reference]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("confint.mfrm_multivariate_gstudy", run.dontrun = TRUE)
#' \dontrun{
#' # Published synthetic data: 10 persons, 6 common tasks, scores V and W.
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#' ci <- confint(g, weights = c(V = .5, W = .5), assumption = "normal")
#' ci[, c("Source", "Estimate", "SE", "Lower", "Upper", "Status")]
#' attr(ci, "notes")
#' # One original score and one component; no new G-study is fitted.
#' confint(g, parm = "Person", score = "V", assumption = "normal", level = .90)
#' # Covariance between the two observed scores, retaining its sign.
#' covariance_ci <- confint(g, score = c("V", "W"), method = "gci",
#'   assumption = "normal")
#' covariance_ci[, c("Source", "Estimate", "SE", "Lower", "Upper", "Status")]
#' attr(covariance_ci, "monte_carlo")$endpoint_bands
#' # Incomplete ratings: a residual-only interval without likelihood fitting.
#' observed <- tasks[(tasks$Person + tasks$Task) %% 3 != 0, ]
#' sparse_g <- mfrm_g_study(observed, c("V", "W"), rater = NULL,
#'                                    method = "minque0")
#' residual_ci <- confint(sparse_g, parm = "Residual", score = "V",
#'   method = "residual_contrast", assumption = "normal")
#' residual_ci[, c("Estimate", "IntervalCenter", "Lower", "Upper", "Status")]
#' attr(residual_ci, "residual_projection")$df
#' # A separate experimental REML model can address other components.
#' if (requireNamespace("lme4", quietly = TRUE) && requireNamespace("minqa", quietly = TRUE)) {
#'   profile_ci <- confint(sparse_g, score = "V", method = "reml_profile",
#'                         assumption = "normal")
#'   profile_ci[, c("Source", "Estimate", "IntervalCenter", "Lower", "Upper", "SEStatus")]
#'   attr(profile_ci, "notes")
#'   # Two scores require a joint model; this is covariance, not correlation.
#'   joint_ci <- confint(sparse_g, parm = "Person", score = c("V", "W"),
#'                       method = "reml_profile", assumption = "normal")
#'   joint_ci[, c("Source", "Estimate", "IntervalCenter", "Lower", "Upper", "Validation")]
#' }
#' }
#' @export
confint.mfrm_multivariate_gstudy <- function(object, parm = NULL, level = .95,
                                            score = NULL, weights = NULL,
                                            assumption, method = "auto",
                                            draws = 20000L, seed = 20261016L, ...) {
  rlang::check_dots_empty()
  if (missing(assumption) || !identical(assumption, "normal"))
    stop("This method requires `assumption = 'normal'`; it does not test the random-effects assumption.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || is.object(level) ||
      !is.null(dim(level)) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be a number strictly between zero and one.", call. = FALSE)
  if (!is.character(method) || length(method) != 1L || is.na(method))
    stop("`method` must name one interval method.", call. = FALSE)
  if (method == "residual_contrast") {
    if (!missing(draws) || !missing(seed))
      stop("`draws` and `seed` apply only to `method = 'gci'`.", call. = FALSE)
    return(.mfrm_gstudy_interval_context(.mfrm_mvgt_residual_intervals(object, parm, score, weights, level), object))
  }
  if (is.character(score) && length(score) == 2L && is.null(weights)) {
    if (!method %in% c("gci", "regression_region", "reml_profile"))
      stop("For two-score covariance intervals, explicitly select `method = 'gci', 'regression_region', 'reml_profile' or residual-only 'residual_contrast'.", call. = FALSE)
    if (method != "gci" && (!missing(draws) || !missing(seed)))
      stop("`draws` and `seed` apply only to `method = 'gci'`.", call. = FALSE)
    if (method == "reml_profile")
      return(.mfrm_gstudy_interval_context(.mfrm_mvgt_covariance_profile(object, parm, score, level), object))
    out <- .mfrm_mvgt_covariance_intervals(object, score, level, method, draws, seed)
    return(.mfrm_gstudy_interval_context(out[.mfrm_mvgt_component_selection(parm, out$Source), , drop = FALSE], object))
  }
  if (!method %in% c("auto", "mls", "reml_profile") || !missing(draws) || !missing(seed))
    stop("Single-score or composite component intervals use `method = 'auto', 'mls', 'reml_profile' or residual-only 'residual_contrast', without draws or seed.", call. = FALSE)
  if (method == "reml_profile")
    return(.mfrm_gstudy_interval_context(.mfrm_mvgt_component_profile(object, parm, score, weights, level), object))
  projection <- .mfrm_mvgt_anova_projection(object, score, weights)
  design <- projection$design
  sources <- rownames(design$A)
  weights <- projection$weights; score <- projection$score; kind <- projection$kind
  ms <- projection$ms; estimate <- projection$estimate
  restore <- projection$restore
  selected <- .mfrm_mvgt_component_selection(parm, sources)
  intervals <- lapply(seq_along(sources), function(i)
    .mfrm_mvgt_linear_variance_interval(ms, design$df, design$A[i, ], level))
  se <- vapply(intervals, `[[`, numeric(1), "se") * restore
  raw <- t(vapply(intervals, `[[`, numeric(2), "raw")) * restore
  status <- vapply(intervals, `[[`, character(1), "status")
  for (i in which(status %in% c("Available", "Includes zero"))) {
    original <- intervals[[i]]$raw
    if (any(!is.finite(raw[i, ])) || any(raw[i, ] == 0 & original > 0)) {
      status[i] <- "Unrepresentable component endpoints"
    } else if (raw[i, 1L] == raw[i, 2L] && original[1L] < original[2L]) {
      status[i] <- "Component interval below numeric resolution"
    }
  }
  estimate <- estimate * restore
  if (any(!is.finite(c(estimate, se))))
    stop("Estimates or SEs exceed numeric range; rescale the scores or weights.", call. = FALSE)
  lower <- pmax(raw[, 1L], 0)
  upper <- raw[, 2L]
  unavailable <- !status %in% c("Available", "Includes zero")
  lower[unavailable] <- upper[unavailable] <- NA_real_
  out <- data.frame(Source = sources, Score = score, Estimate = unname(estimate),
    SE = unname(se), Lower = lower, Upper = upper, RawLower = raw[, 1L],
    RawUpper = raw[, 2L], Level = level,
    Method = vapply(intervals, `[[`, character(1), "method"), Status = status)
  out <- out[selected, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "weights") <- weights
  attr(out, "kind") <- kind
  attr(out, "assumption") <- assumption
  attr(out, "source_design") <- object$design[c("facets", "counts", "nesting",
    "child_counts", "balanced", "complete", "method", "score_convention")]
  attr(out, "source_rows") <- object$data_usage$counts
  attr(out, "notes") <- c(
    "Pointwise normal-random-effects intervals for a prespecified score or fixed composite; MLS intervals are approximate. Not simultaneous or distribution robust.",
    "Estimates and SEs have squared score units. SE is sampling uncertainty in a component, not measurement SEM or uncertainty in G/Phi.",
    "Raw estimates and endpoints are retained. Lower/Upper intersect the interval with nonnegative variances; consult Status for empty or unavailable intervals.",
    "Positive intervals whose width or endpoints cannot be represented remain unavailable with an explicit Status; a rounded zero width is not exact population knowledge.",
    "SE uses the df + 2 moment correction. Weights are used as supplied; saveRDS() preserves assumptions and source context.")
  .mfrm_gstudy_interval_context(out, object)
}

# Validate and project stored ANOVA summaries for both uncertainty methods.
.mfrm_mvgt_anova_projection <- function(object, score = NULL, weights = NULL) {
  design <- .mfrm_mvgt_anova_summary(object)
  sources <- rownames(design$A)
  scores <- object$design$scores
  selection <- .mfrm_mvgt_interval_score(scores, score, weights)
  weights <- selection$weights; score <- selection$score; kind <- selection$kind
  wscale <- max(abs(weights))
  w <- weights / wscale
  # Normalize matrices and weights before multiplication, preserving score units.
  active <- which(w != 0)
  mscale <- max(abs(unlist(lapply(object$mean_products,
    function(m) m[active, active, drop = FALSE]))))
  unit <- if (mscale == 0) 1 else mscale
  restore <- unit * wscale * wscale
  if (!is.finite(restore) || restore == 0)
    stop("Projected variances exceed numeric range; rescale the scores or weights.", call. = FALSE)
  project <- function(m) as.numeric(crossprod(w[active],
    (m[active, active, drop = FALSE] / unit) %*% w[active]))
  ms <- vapply(object$mean_products, project, numeric(1))
  estimate <- vapply(object$components, project, numeric(1))
  # Gram matrices may yield a tiny negative quadratic form through cancellation.
  tol <- 64 * .Machine$double.eps * length(scores)^2
  if (any(!is.finite(c(ms, estimate))) || any(ms < -tol))
    stop("Stored mean products give invalid mean squares for the selected score.", call. = FALSE)
  ms <- pmax(ms, 0)
  mapped <- drop(design$A %*% ms)
  if (any(abs(estimate - mapped) > tol * pmax(1, rowSums(abs(design$A)))))
    stop("Stored components disagree with their ANOVA mean squares; recreate the G-study.", call. = FALSE)
  list(design = design, weights = weights, score = score, kind = kind,
    ms = ms, estimate = estimate, mapped = mapped, restore = restore)
}

.mfrm_mvgt_interval_score <- function(scores, score, weights) {
  if (!is.null(score) && !is.null(weights))
    stop("Choose either `score` or `weights`, not both.", call. = FALSE)
  kind <- if (is.null(weights)) "Score" else "Composite"
  if (is.null(weights)) {
    if (is.null(score)) score <- scores[1L]
    if (!is.character(score) || length(score) != 1L || is.na(score) || !score %in% scores)
      stop("`score` must name one original score.", call. = FALSE)
    weights <- setNames(as.numeric(scores == score), scores)
  } else {
    if (!is.numeric(weights) || is.object(weights) || is.complex(weights) ||
        !is.null(dim(weights)) || length(weights) != length(scores) ||
        is.null(names(weights)) || anyNA(names(weights)) || anyDuplicated(names(weights)) ||
        !setequal(names(weights), scores) || any(!is.finite(weights)) || !any(weights != 0))
      stop("`weights` must name every score once, with finite values and at least one nonzero weight.", call. = FALSE)
    weights <- weights[scores]
    score <- "Composite"
  }
  list(weights = weights, score = score, kind = kind)
}

.mfrm_mvgt_component_selection <- function(parm, sources) {
  selected <- seq_along(sources)
  if (!is.null(parm)) {
    if (is.character(parm)) selected <- match(parm, sources) else {
      if (!is.numeric(parm) || is.complex(parm) || is.object(parm) ||
          !is.null(dim(parm)) || any(!is.finite(parm)) || any(parm != floor(parm)))
        stop("`parm` must contain component names or integer positions.", call. = FALSE)
      selected <- parm
    }
    if (!length(selected) || anyNA(selected) || anyDuplicated(selected) ||
        any(selected < 1 | selected > length(sources)))
      stop("`parm` must select distinct existing components.", call. = FALSE)
  }
  selected
}

.mfrm_mvgt_anova_summary <- function(object) {
  if (!inherits(object, "mfrm_multivariate_gstudy") ||
      !identical(object$design$method, "Balanced multivariate ANOVA") ||
      !isTRUE(object$design$balanced) || !isTRUE(object$design$complete))
    stop("MLS component intervals require a complete balanced ANOVA G-study; for MINQUE(0), explicitly select 'reml_profile' or residual-only 'residual_contrast' after reviewing its limitations.", call. = FALSE)
  design <- .mfrm_mvgt_anova_interval_design(object$design)
  sources <- rownames(design$A)
  scores <- object$design$scores
  valid_matrices <- function(m) is.list(m) && identical(names(m), sources) &&
    all(vapply(m, function(a) is.matrix(a) && is.numeric(a) && !is.complex(a) &&
      identical(dimnames(a), list(scores, scores)) && all(is.finite(a)) &&
      identical(a, t(a)), logical(1)))
  if (!is.character(scores) || !length(scores) || anyNA(scores) ||
      any(!nzchar(scores)) || anyDuplicated(scores) ||
      !valid_matrices(object$components) || !valid_matrices(object$mean_products) ||
      !identical(object$degrees_of_freedom$Source, sources) ||
      !isTRUE(all.equal(object$degrees_of_freedom$DF, unname(design$df), tolerance = 0)))
    stop("Stored ANOVA matrices or degrees of freedom are incomplete or inconsistent with the design; recreate the G-study.", call. = FALSE)
  design
}

# Mean-square maps for the balanced designs supported by the point estimator.
.mfrm_mvgt_anova_interval_design <- function(design) {
  n <- design$counts
  factors <- names(n)
  if (!is.numeric(n) || is.complex(n) || !is.null(dim(n)) ||
      !length(n) || any(!is.finite(n)) || any(n < 2 | n != floor(n)) ||
      !identical(factors[1L], "Person") || !.mfrm_mvgt_valid_facets(factors[-1L]))
    stop("The G-study factor counts or identities are incomplete or altered.", call. = FALSE)
  nesting <- .mfrm_mvgt_nesting(design$nesting, factors[-1L])
  sources <- .mfrm_mvgt_sources(factors, nesting)
  if (!is.null(nesting)) {
    P <- n[1L]; T <- n[unname(nesting)]
    children <- design$child_counts
    if (!is.numeric(children) || length(children) != T || any(!is.finite(children)) ||
        any(children < 2 | children != floor(children)) || length(unique(children)) != 1L ||
        sum(children) != n[names(nesting)])
      stop("The balanced nested child counts are incomplete or inconsistent.", call. = FALSE)
    R <- children[1L]
    A <- rbind(c(1, 0, 0, -1, 0) / (R * T), c(0, 1, -1, -1, 1) / (P * R),
      c(0, 0, 1, 0, -1) / P, c(0, 0, 0, 1, -1) / R, c(0, 0, 0, 0, 1))
    df <- c(P - 1, T - 1, T * (R - 1), (P - 1) * (T - 1), T * (P - 1) * (R - 1))
  } else {
    subsets <- if (length(n) == 2L) list(1L, 2L, 1:2) else
      list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
    L <- matrix(0, length(subsets), length(subsets))
    for (i in seq_along(subsets)) for (j in seq_along(subsets)) {
      if (all(subsets[[i]] %in% subsets[[j]]))
        L[i, j] <- prod(n[setdiff(seq_along(n), subsets[[j]])])
    }
    A <- solve(L)
    df <- vapply(subsets, function(s) prod(n[s] - 1), numeric(1))
  }
  dimnames(A) <- list(sources, sources)
  list(A = A, df = unname(df))
}

# Ting MLS lower width (Hess & Iyer, 2002, pp. 2--3).
# The upper width is the lower width for the negated linear combination.
.mfrm_mvgt_mls_width <- function(ms, df, a, tail) {
  p <- which(a > 0); n <- which(a < 0)
  z <- abs(a) * ms
  scale <- max(z)
  if (scale == 0) return(0)
  z <- z / scale
  G <- 1 - df / stats::qchisq(tail, df, lower.tail = FALSE)
  H <- df / stats::qchisq(tail, df) - 1
  terms <- c((G[p] * z[p])^2, (H[n] * z[n])^2)
  for (i in p) for (j in n) {
    f <- stats::qf(tail, df[i], df[j], lower.tail = FALSE)
    cross <- f - 2 + 1/f - G[i]^2 * f - H[j]^2/f
    terms <- c(terms, cross * z[i] * z[j])
  }
  if (length(p) > 1L && length(n)) {
    for (j in seq.int(2L, length(p))) for (i in seq_len(j - 1L)) {
      s <- p[i]; t <- p[j]; total <- df[s] + df[t]
      joint <- 1 - total / stats::qchisq(tail, total, lower.tail = FALSE)
      cross <- (joint^2 * total^2 / (df[s] * df[t]) -
        G[s]^2 * df[s]/df[t] - G[t]^2 * df[t]/df[s]) / (length(p) - 1L)
      terms <- c(terms, cross * z[s] * z[t])
    }
  }
  value <- sum(terms)
  if (!is.finite(value) || value < -64 * .Machine$double.eps * sum(abs(terms)))
    return(NA_real_)
  sqrt(max(0, value)) * scale
}
