#' Pointwise intervals for fixed facet estimates and contrasts
#'
#' Describe uncertainty about a fixed facet estimate, such as rater severity,
#' or a specified difference between two raters. Start from a supported fitted
#' model; the point estimates stay the same when you change the interval method.
#'
#' Compare ordinary observed-information intervals with a sandwich covariance
#' that treats a person's complete response vector, or an explicitly declared
#' larger cluster, as the independent sampling unit. Estimates are not refitted.
#'
#' @param fit An inference-ready RSM/PCM MML fit from [fit_mfrm()], using a
#'   fixed standard-normal person distribution, fixed quadrature and unit
#'   observation weights. A fit whose category review concerns singleton counts
#'   can additionally use `method = "model"` after the separate checks described
#'   under "Small samples and singleton categories". Experimental two-family GPCM MML also supports
#'   `method = "model"` with fixed-grid EM or direct adaptive integration;
#'   see "Two-family location targets" below. Direct MML native
#'   locations additionally support estimated normal RSM/PCM populations,
#'   including conditional-normal latent regression,
#'   and one-family GPCM with fixed N(0,1) or an estimated intercept-only normal
#'   population, using fixed or adaptive integration; see "Native locations"
#'   below for their scope and separate numerical checks. Native RSM/PCM CML
#'   results additionally support `method = "model", scale = "native"`, including
#'   the saved targets from `fit_scope = "comparisons"`.
#'   Matching full RSM/PCM MML interaction fits use the `comparison_review`
#'   extension below; their explicit model and source checks remain separate.
#' @param facet A non-person facet, for example `"Rater"`. May be omitted
#'   when `comparison_review` supplies a facet-by-origin interaction target,
#'   or for a CML comparison-only fit, which already retains its targets.
#' @param contrasts Optional numeric matrix with distinct target row names
#'   and columns named by every level of `facet`. Columns are aligned by name.
#'   A row with coefficients `c(1, -1, 0)` estimates the first level minus the
#'   second. By default, each level is reported.
#' @param method `"model"` (default) uses ordinary observed information;
#'   `"sandwich"` uses independent-cluster marginal-likelihood scores.
#'   The native-location and experimental two-family extensions support only
#'   `"model"`; sandwich support remains limited to the ordinary fixed-population
#'   RSM/PCM route.
#' @param clusters For sandwich inference, an optional data frame with exactly
#'   one row per fitted person and complete `Person` and `Cluster` identifiers.
#'   Without it, persons are the independent units. A larger cluster could be
#'   a school or clinic if different schools or clinics are independent.
#'   A person cannot be split across clusters. Supply original person labels,
#'   even when the input person column has another name.
#' @param adjust Logical; with `method = "sandwich"`, optionally multiply
#'   the covariance by `G/(G-1)`, where `G` is the cluster count. Default `FALSE`.
#'   This scaling does not guarantee small-sample coverage.
#' @param level Pointwise confidence level, strictly between zero and one.
#' @param scale `"native"` (default) retains the fitted ability units.
#'   `"standardized"` divides an explicit within-facet zero-sum contrast, or
#'   the interaction contrasts declared by `comparison_review`, by
#'   the population SD, using full joint delta-method uncertainty. Requires
#'   `method = "model"` and a fixed N(0,1) or estimated intercept-only normal
#'   population; individual locations and latent regression are not supported.
#' @param comparison_review Optional [describe_mfrm_data()] review or its
#'   `$comparison_review` for a matching full native RSM/PCM MML interaction fit.
#'   Uses its declared facet-by-origin contrasts instead of facet locations;
#'   omit `contrasts` and optionally `facet`. Requires `method = "model"`.
#'   See "MML interaction contrasts" for scope and assumptions.
#'   Omit this argument for a CML comparison-only fit; its saved review is reused.
#'
#' @param interaction Optional exact name of one fitted two-facet interaction,
#'   for example `"Rater:Criterion"`. Omit `facet` and `comparison_review`.
#'   This selects jointly fitted cell effects rather than fixed-baseline bias
#'   screens. Requires `method = "model"`, `scale = "native"`, full fixed-grid
#'   RSM/PCM MML (direct, EM or hybrid) with fixed N(0,1), one interaction,
#'   no anchors, groups, dummy facets or shrinkage. Every engine must pass the
#'   common numerical and full-joint information checks; EM stopping alone
#'   does not qualify intervals. See the joint-cell section below.
#' @param anchor_covariance Optional full sampling covariance of the supplied
#'   anchor estimates, including their correlations, in native ability units
#'   squared. Row and column names must match every directly anchored level
#'   on the one noncentered facet. Requires the fixed-grid direct RSM/PCM
#'   anchor-origin route with an intercept-only population and fixed SD.
#'   `NULL` retains the existing conditional-on-anchors calculation.
#' @param independent_anchors Logical, default `FALSE`. Set `TRUE` with
#'   `anchor_covariance` only when source calibration and receiving responses
#'   are independent. Shared respondents or dependent clusters require
#'   cross-sample covariance and are unsupported. This is an explicit
#'   assumption, not a check inferred from different Person identifiers.
#'
#' @return An `mfrm_facet_intervals` object with `table`, selected target
#'   `covariance`, `model_covariance`, free-parameter `parameter_covariance`,
#'   `model_parameter_covariance`, `person_scores`, `cluster_scores`, the
#'   `clusters` mapping, exact `contrasts`, `settings`, and source `fit`.
#'   The model-based SE and interval remain alongside the selected method.
#'   If ordinary MML joint information is unavailable, free-parameter
#'   covariances are `NULL` and target covariances contain `NA` (zero for
#'   constraint-fixed rows/columns); `unavailable_reason` and `checks` explain
#'   the failure. Sandwich scores and their rank are then uncomputed, not zero.
#'   Score rows are derivatives of a person's marginal log likelihood, not
#'   observed category scores or ability estimates. Scores are absent for
#'   `method = "model"`.
#'   When weak information passes numerical refinement, `cautions` and
#'   `information_review` retain the warning and checks. `InferenceCaution`
#'   also accompanies the interval table when applicable.
#'   Each row retains `Status`, `CIEligible` and a plain-language
#'   `InferenceReview`. `CIEligible` concerns the selected interval, not the
#'   validity of the point estimate or a guarantee of coverage. `fixed` means
#'   the target is fixed by constraints, with no inferential interval.
#'   `insufficient_cluster_rank` withholds sandwich bounds because independent
#'   cluster scores do not span all free parameters. The point and separate
#'   `ModelSE`/`ModelLower`/`ModelUpper` are retained; those model bounds are
#'   not a replacement robust interval. Review the independent sampling units
#'   and model dimension; do not split dependent clusters to obtain an interval.
#'   `nonpositive_variance` or `nonfinite_bounds` denotes an unusable
#'   target-specific calculation. Missing bounds are not zero uncertainty.
#'   Two-family results additionally retain `checks`, `numerical_checks`,
#'   `score_rank`, source identity, scale references and interval refusal reasons.
#'   Their table includes `CIEligible`; failed numerical checks leave bounds
#'   and SEs missing while preserving the fitted point estimates.
#'   Native-location results retain checks, numerical comparisons, full and
#'   target covariance, target Jacobian and source identity, without an
#'   empirical Person-score rank requirement.
#'   Standardized results also retain `NativeEstimate`, the population SD,
#'   whether it was estimated, and the full transformation Jacobian.
#'   `CIEligible` describes numerical availability; `CoverageQualified = FALSE`
#'   does not claim finite-sample coverage for the transformed interval.
#'
#' @section Jointly fitted interaction cells and contrasts:
#' Fit a prespecified two-facet interaction with `fit_mfrm(facet_interactions =
#' "Rater:Criterion")`, then use `mfrm_facet_intervals(fit,
#' interaction = "Rater:Criterion")`. The default returns all fitted cells,
#' including the dependent last row and column. Cell labels, counts, sparse
#' flags and identification constraints are retained in `interaction_cells`.
#' This bounded route uses a fixed standard-normal ability population; latent
#' regression, estimated population variance, adaptive integration, JML, CML,
#' GPCM and multiple fitted interaction terms are outside this extension.
#'
#' An interaction effect is an additive predictor shift after jointly fitting
#' the main effects. Positive values increase expected coded ratings; they do
#' not automatically imply favorable performance or unfairness. Effects have
#' zero row sums and zero column sums across the declared levels, without
#' weighting those constraints by cell counts. Their values depend on that
#' definition and the included main effects. This differs from [estimate_bias()],
#' which holds the fitted calibration and person measures fixed while estimating
#' an additional residual shift. Neither result proves substantive rater bias.
#'
#' Write free interaction coordinates as gamma and their expansion to all cells
#' as G gamma. Cell covariance is G V_gamma G', where V_gamma is the interaction
#' block of the inverse full marginal information, including main effects and
#' thresholds. It is not the inverse of only the interaction information block.
#' Dependent cells consequently have correlated estimates and a singular joint
#' cell covariance; this alone does not invalidate individual pointwise intervals.
#' For a prespecified cell contrast C, the covariance is C G V_gamma G' C'.
#' Supply every cell `Target` as a named column of `contrasts`, including zeros.
#' A contrast fixed at zero by the row/column constraints has `Status = "fixed"`
#' and no inferential interval; it is not evidence for a null interaction.
#'
#' Source consistency, identification/boundary records, full information,
#' stationarity and finer-quadrature checks must pass. Unavailable intervals
#' retain finite point estimates and reasons. Sparse-cell flags are descriptive,
#' not a sample-size guarantee. Full-design identification and numerical success
#' do not establish nominal coverage. These are pointwise Wald intervals, not
#' simultaneous bounds or multiplicity-adjusted cell decisions. Saved tables,
#' plots, [apa_table()] and [mfrm_results()] retain the same values and cautions.
#'
#' A bundled-data example uses the existing estimation and interval APIs:
#' ```r
#' ratings <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   model = "PCM", step_facet = "Criterion", noncenter_facet = "Criterion",
#'   facet_interactions = "Rater:Criterion", quad_points = 61)
#' cells <- mfrm_facet_intervals(fit, interaction = "Rater:Criterion")
#' summary(cells)[, c("Target", "Estimate", "SE", "Lower", "Upper", "Status")]
#' plot(cells)
#' # Compare the first two displayed cells; choose the scientific contrast first.
#' C <- matrix(0, 1, nrow(cells$table),
#'   dimnames = list("First minus second", cells$table$Target))
#' C[1, 1:2] <- c(1, -1)
#' difference <- mfrm_facet_intervals(fit, interaction = "Rater:Criterion",
#'   contrasts = C)
#' vcov(difference) # Uses covariance between both cells
#' ```
#'
#' @section Uncertainty from independently estimated anchors:
#' The default intervals condition on the exact anchor values. Supplying
#' `anchor_covariance` and `independent_anchors = TRUE` adds source-calibration
#' uncertainty without changing any point estimates. Use the full covariance,
#' not only a vector of SEs or a diagonal replacement for correlated anchors.
#' Its numerical validity is checked; its statistical validity, units and
#' correspondence to the supplied anchor values remain the analyst's responsibility.
#'
#' Write the free receiving parameters as u, supplied anchors as a, and the
#' negative log likelihood Hessian blocks as H_uu and H_ua. Implicit
#' differentiation gives S = -inverse(H_uu) H_ua. For a reported target t(u,a),
#' its source derivative is D = t_u S + t_a. The selected target covariance is
#' t_u V_conditional t_u' + D V_anchor D'. The direct derivative t_a matters
#' for targets containing anchored elements; adding variance only to the free
#' parameters would omit that contribution and their cross-covariance.
#' This is a first-order approximation for independent samples, not joint
#' recalibration, a bootstrap or a correction for anchor drift or biased anchors.
#'
#' Source derivatives use central differences of the analytic marginal score,
#' checked by halving the difference step and refining the quadrature. Failed
#' source-derivative checks preserve point estimates and the conditional
#' results, with unavailable selected uncertainty. `ConditionalSE`,
#' `ConditionalLower`, `ConditionalUpper` and `ConditionalStatus` retain that
#' separate calculation; facet `Model*` columns also retain it. `AnchorSourceSE`
#' is the square root of the added variance contribution, not a separate
#' confidence interval. The selected `SE`, bounds and `vcov()` include both
#' components. A conditionally fixed anchor may therefore have a selected
#' interval for its source-estimated value. A zero covariance contribution does
#' not create an interval for a target fixed under both specifications.
#'
#' `anchor_uncertainty` retains the supplied and numerically symmetrized/PSD
#' covariance, any floating-point adjustment, anchor values, sensitivity,
#' independent-sample assumption and derivative checks. Singular PSD matrices
#' are allowed. The free-parameter covariance includes source propagation;
#' it does not include anchored values as extra optimizer coordinates.
#' Saved population intervals can reuse this sensitivity and the conditional
#' information: `mfrm_population_intervals(rater_ci)` propagates it once, without
#' refitting or differentiating again. Reports and exported tables retain the
#' supplied covariance and source-derivative checks; RDS retains the full result.
#'
#' Only fixed-grid direct RSM/PCM MML with unit weights, no interactions,
#' shrinkage or groups, an intercept-only population, fixed SD and all anchors
#' on its noncentered facet is supported. Adaptive integration, free population
#' variance, GPCM, dependent calibration samples and uncertainty in the supplied
#' SD are outside this extension. No nominal coverage guarantee follows from
#' numerical agreement or positive covariance. Other facets retain their
#' centering; changing an item pool still requires a defensible scale link.
#'
#' A complete interface example starts from bundled ratings. This split is
#' illustrative: independent Persons and stable anchors must be justified for
#' the actual sampling design, and different IDs alone do not establish them.
#' ```r
#' ratings <- load_mfrmr_data("example_core")
#' ids <- unique(ratings$Person)
#' source <- ratings[ratings$Person %in% ids[1:20], ]
#' receiving <- ratings[!ratings$Person %in% ids[1:20], ]
#' baseline <- fit_mfrm(source, "Person", c("Rater", "Criterion"), "Score",
#'   model = "RSM", noncenter_facet = "Rater", quad_points = 61)
#' baseline_ci <- mfrm_facet_intervals(baseline, "Rater")
#' common <- c("R01", "R02")
#' anchors <- data.frame(Facet = "Rater", Level = common,
#'   Anchor = baseline_ci$table$Estimate[match(common, baseline_ci$table$Target)])
#' cohort <- fit_mfrm(receiving, "Person", c("Rater", "Criterion"), "Score",
#'   model = "RSM", anchors = anchors, noncenter_facet = "Rater",
#'   population_formula = ~1, person_data = unique(receiving["Person"]),
#'   population_sd = 1, quad_points = 61)
#' rater_ci <- mfrm_facet_intervals(cohort, "Rater",
#'   anchor_covariance = vcov(baseline_ci)[common, common], independent_anchors = TRUE)
#' summary(rater_ci) # Selected and conditional intervals side by side
#' population_ci <- mfrm_population_intervals(rater_ci)
#' ```
#'
#' @details For cluster score `s_g` and observed negative-log-likelihood
#'   Hessian `H`, the unadjusted sandwich is
#'   `solve(H) %*% sum_g(s_g %*% t(s_g)) %*% solve(H)`.
#'   Scores aggregate all responses of a person before forming the outer
#'   products; grouping individual rating rows would be a different and
#'   incorrect calculation for this marginal likelihood. The full covariance
#'   is transformed through the fitted constraints and requested contrasts.
#'   Both methods use normal critical values. No automatic method selection,
#'   multiplicity adjustment or significance flag is supplied.
#'
#'   The sandwich requires many independent sampling units and appropriate
#'   regularity. It permits dependence within the declared unit but does not
#'   establish independence between units. A small number of units can give
#'   poor intervals even with nonsingular covariance. If their scores do not
#'   span the free-parameter space, sandwich intervals are unavailable; point
#'   estimates and the reason remain. For ordinary fixed-population RSM/PCM,
#'   an ineligible source fit still causes an error. If the source is eligible
#'   but joint information is unavailable, nonpositive or regularized, finite
#'   points remain with `Status = "information_unavailable"`,
#'   `CIEligible = FALSE`, missing model and sandwich uncertainty, and the
#'   explanation in `InferenceReview`. This includes information-workspace
#'   limits. A sandwich also needs the inverse information; switching methods
#'   cannot repair its absence. Failed information can limit confidence in
#'   point accuracy along weak directions. No inverse is regularized for this
#'   output. Invalid arguments and cluster mappings remain errors.
#'   Targets fixed by constraints
#'   have no inferential interval. Known anchors exclude anchor uncertainty.
#'   An ill-conditioned but numerically verified unregularized inverse can be
#'   used with a warning. Successful inversion does not establish reliable
#'   normal intervals. Review interval width, boundary proximity and quadrature
#'   sensitivity; changing to sandwich covariance does not remove this concern.
#'
#' @section Small samples and singleton categories:
#' A small sample or a wide interval is not, by itself, a reason to suppress an
#' estimate or an otherwise available asymptotic interval. RSM/PCM model
#' intervals can retain singleton categories (one observed response in a
#' fitted or local category cell) when all declared categories are represented
#' in the fitted and local scopes. They require fresh category review, source
#' and convergence checks, positive unregularized joint information and a
#' finer-quadrature score/covariance comparison. Empty categories, unsupported
#' coordinates, known identification failures and failed numerical checks are
#' not excused by this route. Failed interval checks retain the point with
#' unavailable uncertainty and a reason.
#'
#' The fit keeps its category-review state and can still have
#' `InferenceReady = FALSE`. Inspect this interval object's `CIEligible`,
#' `checks`, `numerical_checks` and `InferenceCaution` for the specific calculation.
#' This does not qualify other analyses or guarantee small-sample coverage.
#' The same model covariance is used by `vcov(fit)` and default `confint(fit)`.
#' Singleton counts do not gain a sandwich or GPCM Wald extension here.
#' For supported one-Person RSM/PCM calibration, model-information intervals
#' may be available after the same numerical checks. Returned `cautions` and
#' `InferenceCaution` state that calibration relies on one response vector and
#' its assumed conditional independence; `CoverageQualified` is `FALSE`.
#' A finite point or information matrix does not establish small-sample
#' coverage or between-Person variation. One Person cannot provide sandwich
#' uncertainty; its ratings must not be relabeled as independent Persons.
#'
#' @section MML assumptions and numerical approximation:
#' MML integrates each Person's joint response probability under the declared
#' ability distribution. Model-information intervals require the response and
#' population models to hold for the actual rating assignments and covariates,
#' including local independence given ability. An overall normal ability
#' distribution does not by itself justify using the same distribution within
#' rosters selected according to ability. Response-dependent assignment,
#' stopping or missingness requires a separate justification of the likelihood.
#'
#' The regular asymptotic interpretation requires independent Persons, an
#' identifiable fixed-dimensional marginal model, a consistent interior
#' solution, smooth likelihood derivatives and information growing in every
#' free direction. With bounded ratings per Person and stable design
#' proportions, a positive-definite limiting information per Person is a
#' sufficient information-growth condition. Correctly specified MML can then
#' support inference with fixed exposure per Person; increasing Person count
#' alone cannot repair misspecification or an unidentified design. A positive
#' observed Hessian at one fitted vector does not prove global identification.
#'
#' The continuous-normal likelihood is evaluated numerically. A fixed finite
#' quadrature order can retain approximation error as Person count grows.
#' Optimizer and integration error must be negligible relative to sampling
#' uncertainty for the continuous-model Wald interpretation. Agreement at two
#' quadrature orders is a sensitivity diagnostic, not an error bound or a
#' proof of asymptotic equivalence. Standardization does not remove either
#' numerical approximation error or bias from a misspecified population model.
#' The ordinary sandwich route can estimate variation around a working-model
#' target under its own regularity conditions; that target need not be the
#' true structural parameter.
#'
#' @section Conditional likelihood:
#' Native `fit_mfrm(method = "CML")` results use the inverse full conditional
#' information, transformed through the fitted facet constraints and requested
#' contrasts. Other facets and steps remain estimated. Only `method = "model"`,
#' `clusters = NULL`, `adjust = FALSE` and `scale = "native"` are supported.
#' CML estimates no ability population SD: `scale = "standardized"` cannot be
#' computed from that fit. For MML with an estimated intercept-only normal
#' population, that option divides by the same fit's estimated SD and includes
#' its covariance with the contrast; see "Population-SD-standardized contrasts".
#' A full CML fit can retain finite numerical point estimates when its
#' identification, finite-MLE and stationarity checks pass but conditional
#' information fails. Its non-fixed targets have `Status = "information_unavailable"`,
#' missing SE/bounds, `CIEligible = FALSE` and the reason in `InferenceReview`.
#' The fit's `checks$information_status` and `checks$information_reason` record
#' the information decision. Constraint-fixed targets retain zero SE and no
#' inferential interval. An information warning also limits confidence in the
#' numerical accuracy of points along weak directions; retaining a point does
#' not qualify it for inference.
#' CML intervals are asymptotic, requiring the correct conditional response law,
#' independent Persons and information growing in every free direction of a
#' fixed-dimensional identified model. More Persons alone is insufficient;
#' extreme totals add no information. Response-dependent assignment, stopping
#' or missingness requires justification that the conditional law is preserved.
#' See the sufficient regular regime in [fit_mfrm()]. Neither numerical
#' acceptance nor a distribution-free ability treatment ensures small-sample
#' coverage. `plot(ci)`, `summary(ci)` and RDS reuse apply to CML interval objects;
#' attaching them to ordinary `mfrm_results()` is not supported.
#' Non-fixed CML targets require finite, strictly positive representable variance.
#' Rescale contrast coefficients if their units cause variance underflow or
#' overflow; numerical zero variance is not a constraint-fixed target.
#'
#' For `fit_scope = "comparisons"`, call `mfrm_facet_intervals(fit)` without
#' new `contrasts` or `comparison_review`. The saved estimates, joint target
#' covariance and qualifications are reused without optimization or information
#' recalculation. `level` changes only the pointwise normal critical value.
#' Unidentified, boundary or unqualified-stationarity comparisons keep missing
#' estimates and intervals. If only information failed, identified numerical
#' points remain, with missing uncertainty, their reasons and
#' `CIEligible = FALSE`; no unavailable value is zero-filled.
#' `parameter_covariance` is the covariance of the fit's identified subspace,
#' not of separately identified cell effects. Use `summary(ci)`, `plot(ci)`,
#' `apa_table(ci)` and `saveRDS(ci)`; the saved review retains group order,
#' contrast weights and observed/empty cells. JML/MML comparison-only fits
#' remain descriptive and are not covered by this interval route.
#'
#' @section Anchors and group means:
#'   The ordinary fixed-population RSM/PCM route supports location anchors and
#'   group means. Write the fitted facet locations as \eqn{b=b_0+J\eta} and a
#'   requested comparison as \eqn{Cb}. Its covariance is
#'   \eqn{CJ V_\eta J^{\mathsf T}C^{\mathsf T}}, where \eqn{V_\eta} is extracted
#'   after computing the full joint covariance. Other facets and thresholds
#'   remain estimated; inverting only the selected facet's information would
#'   omit their estimation uncertainty. The affine offset \eqn{b_0} affects
#'   estimates, but has no variance when its supplied values are treated as exact.
#'
#'   A comparison with \eqn{CJ=0} is fixed: its `Status` is `"fixed"`, `SE`
#'   is zero, and interval bounds are missing. This includes values implied by
#'   centering, not only rows in an anchor table. For example, with four raters,
#'   whole-facet centering and a fixed mean of 0.10 for the first three, the
#'   fourth rater is fixed at -0.30. The three group members can still differ
#'   unless other constraints fix them; their mean has no inferential interval.
#'   A difference between a free rater and an anchored rater can remain uncertain.
#'   These are conditional statements, not evidence that the source calibration
#'   or group-mean assumption is error-free. Native-location anchors use only
#'   the fixed-SD RSM/PCM extension below; groups remain excluded there.
#'
#' @section Native locations:
#'   For estimated normal RSM/PCM and intercept-only normal one-family GPCM MML,
#'   `method = "model"` supplies pointwise Wald intervals for
#'   locations or contrasts on the model's native ability scale. One-family
#'   GPCM also admits the explicit fixed-N(0,1) restriction
#'   (`gpcm_mml_identification = "fixed_standard_normal"`). Shared or separate
#'   slope/step owners are supported, with geometric-mean-one relative slopes.
#'   Use direct MML, fixed or adaptive integration, unit weights, additive
#'   centered facets and centered steps. RSM/PCM also accepts a fitted full-rank
#'   conditional-normal regression, retaining its coefficient and residual-variance
#'   cross-information when computing structural uncertainty. The saved population
#'   table identifies coefficients and residual variance; it does not label
#'   regression coefficients as marginal population means. GPCM latent
#'   regression and sandwich covariance remain unavailable here. The native
#'   location procedure performs its own qualification even when the source
#'   fit has `InferenceReady = FALSE`. Read `CIEligible`, `Status`, `checks`,
#'   `numerical_checks` and cautions on this result; interval availability does
#'   not change the source fit's overall readiness or qualify other procedures.
#'   Interactions use the separate `comparison_review` route below.
#'   Population-coefficient and variance intervals are separate targets and are
#'   not returned by this helper.
#'   RSM/PCM populations may fix residual SD with `population_sd`. For an
#'   intercept-only population (`population_formula = ~1`), a direct anchor on
#'   the selected non-Person `noncenter_facet` identifies the origin while free
#'   levels on that facet remain uncentered. These anchored intervals require
#'   no groups or dummy facets and condition on the supplied anchors and SD;
#'   their source-calibration uncertainty is not propagated. Other facets keep
#'   their declared centered constraints. Fixed SD has no variance coordinate
#'   in the joint information.
#'
#'   The covariance transforms the inverse full joint marginal information,
#'   including all estimated population, slope and step coordinates. A
#'   location-only inverse or independent printed SEs cannot replace it.
#'   The source, category support, local solution and unregularized information
#'   must pass checks; covariance and score changes are also checked at a finer
#'   integration order at the same fitted parameter vector. This model-information
#'   procedure has no empirical Person-score rank requirement.
#'   For RSM/PCM populations with estimated variance, a zero-residual-variance candidate at
#'   the retained nuisance estimates is compared with the fitted likelihood.
#'   A higher boundary likelihood beyond the numerical allowance prevents
#'   intervals. This check does not optimize or exclude the entire boundary.
#'   A numerical refusal preserves the native point and its reason with
#'   unavailable SE/bounds. Constraint-fixed targets have no inferential interval.
#'
#'   The asymptotic justification requires independent Persons, a correctly
#'   specified response and population model, an identified interior true
#'   parameter, a fixed number of structural parameters, nonsingular limiting
#'   information and sufficiently accurate integration and optimization.
#'   Information must grow in all target directions as Persons are added;
#'   connected facet labels alone do not establish this condition.
#'   These numerical checks do not establish finite-sample coverage, global
#'   identification or robustness to misspecification. Native locations are
#'   not divided by the estimated population SD unless `scale = "standardized"`
#'   is requested for a zero-sum contrast (see below). A location
#'   difference need not produce a uniform rating difference with unequal slopes
#'   or steps. The calculation does not refit or change the saved fit/readiness.
#'   Save and attach `ci` with `mfrm_results(fit, intervals = list(locations = ci),
#'   compute = "never")`; select `plot(res, type = "facet_locations")` for
#'   RSM/PCM or `"gpcm_locations"` for one-family GPCM. Reports and exports retain
#'   the native scale, inferential qualifications and unavailable outputs.
#'
#' @section Population-SD-standardized contrasts:
#'   Set `scale = "standardized"` and supply explicit contrast rows whose
#'   coefficients sum to zero. For native contrast \eqn{d}, population SD
#'   \eqn{\sigma} and free log-variance \eqn{\lambda = \log(\sigma^2)}, the
#'   target is \eqn{t = d/\sigma}. Its Jacobian has location entries
#'   \eqn{J_d/\sigma} and log-variance entry \eqn{-t/2}. The covariance is
#'   transformed from the inverse full joint observed information, including
#'   the covariance between the contrast and log-variance. Dividing native
#'   interval endpoints by the fitted SD omits this uncertainty. In scalar form,
#'   \deqn{\widehat{\mathrm{Var}}(\widehat t) =
#'     \frac{\widehat{\mathrm{Var}}(\widehat d)
#'       + \widehat d^2\widehat{\mathrm{Var}}(\widehat\lambda)/4
#'       - \widehat d\widehat{\mathrm{Cov}}(\widehat d,\widehat\lambda)}
#'       {\widehat\sigma^2}.}
#'   The covariance term can increase or decrease the SE relative to dividing
#'   the native SE by the fitted SD. A zero estimated contrast can still have
#'   positive variance. With fixed
#'   N(0,1), the result equals the native contrast because the known SD is one.
#'   With an explicit fixed `population_sd`, divide the native contrast and SE
#'   by that SD; no estimated-SD derivative or uncertainty is introduced.
#'   Conditional-normal regression is excluded: its residual SD is not a
#'   marginal population SD. Standardization does not equate expected ratings
#'   across models with different slopes or category steps.
#'
#'   The existing source and numerical checks still apply. These are pointwise
#'   normal delta-method intervals, not profile intervals or finite-sample
#'   coverage guarantees. The regular delta method assumes a fixed true
#'   population SD strictly above zero and an identifiable, interior joint
#'   solution under the stated marginal model. A positive fitted SD alone does
#'   not establish these conditions. At zero SD the target is undefined; near
#'   zero the ratio and normal approximation can be unstable. No universal SD
#'   cutoff certifies validity. Standardization does not resolve weak
#'   identification or boundary problems; see the MML assumptions above.
#'
#' @section MML interaction contrasts:
#' Supply `comparison_review` to obtain native or population-SD-standardized
#' intervals for its explicitly weighted facet-by-origin interaction contrasts.
#' The review must match the fitted responses, group identities, constraints,
#' model, PCM step owner, signs and population design. Its model has one
#' interaction, zero origin main facet and centered remaining facets and
#' interaction margins. This extension requires a full, unadjusted RSM/PCM
#' direct MML fit, unit weights and preserved categories, without shrinkage or
#' additional anchors. Comparison-only MML fits, JML and GPCM are not admitted.
#' CML comparison intervals remain available through [fit_mfrm()] but cannot
#' estimate a population SD.
#'
#' Full conditional-map identification is necessary here. Source tables,
#' categories, local solution, unregularized full information and finer-grid
#' stability must pass the native-location checks above. A failed numerical
#' check retains the point estimate with missing SE/bounds and a reason.
#' Positive conditional rank and local information do not establish global
#' marginal identification or finite-sample interval coverage; the output
#' retains `MarginalIdentificationState` and `CoverageQualified = FALSE`.
#' These are model-based asymptotic intervals, not a calibrated DFF classification.
#'
#' For `scale = "standardized"`, the denominator is the SD estimated in the
#' **same intercept-only normal MML fit** (`population_formula = ~1`). Such a
#' fit assumes one common ability distribution across the declared origin
#' groups. It is a different model from `~Group`, which allows group means.
#' For `~Group`, native interaction intervals are supported, but its residual
#' SD is not substituted for a marginal population SD. A fixed N(0,1) fit has
#' known SD one; no estimated-SD uncertainty is introduced in that case.
#'
#' The joint delta method uses interaction derivatives divided by the fitted
#' SD and derivative `-Estimate/2` for log population variance. All population,
#' step and structural nuisance coordinates remain in the full information
#' inverse. It includes the covariance between the numerator and denominator;
#' dividing the native endpoints or SE by the fitted SD is insufficient.
#'
#' Read `table` with `comparison_review`: group subtraction order, panel weights,
#' observed/empty cells and assumptions are retained. Use `plot(ci)`, `saveRDS(ci)`
#' to retain `IdentificationBasis` and its explanation with the intervals.
#' A free-Person response-map basis does not establish distribution-free MML
#' inference; population-dependent conditional identification does not establish
#' global marginal identification. Display the result directly
#' or attach it as `mfrm_results(fit, intervals = list(school = ci), compute = "never")`.
#' Plotting and reports reuse the saved uncertainty; they do not refit or
#' recalculate information. These results do not extend the separate
#' [mfrm_population_intervals()] consumer to models with fitted interactions.
#'
#' @section Two-family location targets:
#'   For native two-family GPCM MML, select either modeled non-Person facet.
#'   The first slope owner's locations sum to zero; the second owner's locations
#'   are uncentered on the fixed N(0,1) ability scale. These locations are measured
#'   before multiplication by both slope components. Increasing a location
#'   lowers expected ratings when ability, slopes and other effects stay fixed.
#'   A contrast between actual raters need not imply a uniform difference in
#'   expected ratings: their slopes and category-step offsets may also differ.
#'   Neither zero nor a location difference defines rater quality or accuracy.
#'
#'   Let `J` expand the free location coordinates under their fitted constraints
#'   and `C` be the requested contrast matrix. The target covariance is
#'   `C J V J' C'`, where `V` is the location block extracted *after inverting
#'   the full joint marginal information*. Slopes, steps and the other facet
#'   remain estimated nuisance parameters. Inverting a location-only Hessian,
#'   adding independent printed SEs, or using the EM conditional-objective
#'   Hessian would omit their joint estimation uncertainty.
#'
#'   The source fit must pass the same identity, category support, score
#'   derivative, stationarity, full-information and finer-quadrature checks used for
#'   experimental two-family slope intervals in [confint.mfrm_fit()]. The full
#'   covariance is checked; deficient empirical Person-score rank is a retained
#'   diagnostic and warning, not a veto on the inverse-observed-information
#'   method. Passing these numerical checks does not establish
#'   repeated-sample coverage for locations or contrasts. Evidence for component
#'   slope intervals does not qualify these different targets. These are
#'   reported separately as `CIEligible` (numerical interval availability) and
#'   `CoverageQualified = FALSE` in `table`. The latter means that nominal
#'   sampling coverage has not been established for the target and design;
#'   it is neither a coverage estimate of zero nor a numerical-failure flag.
#'   Saving, reporting or exporting the result retains both fields. These are
#'   experimental pointwise normal approximations, with no sandwich, profile,
#'   multiplicity or automatic significance decision. A failed check retains
#'   estimates and reasons with unavailable intervals, even if the coarse-grid
#'   information could be inverted. Constraint-fixed targets have no interval.
#'   The calculation does not refit or change quadrature in the saved fit.
#'
#'   For example, use `ci <- mfrm_facet_intervals(fit, "Rater")`, or supply a
#'   named zero-sum row of `contrasts` for a difference. Plot and save `ci`
#'   directly, or attach it with
#'   `mfrm_results(fit, intervals = list(locations = ci))` and select
#'   `plot(res, type = "gpcm_locations")`. The model must actually contain the
#'   specified facet. Location intervals do not add category-step/curve bands,
#'   Person uncertainty, Wright maps or tests of rater competence. Individual
#'   two-family rater sheets currently display component-slope intervals only;
#'   use the analyst report for these location and contrast results.
#'
#' @section What robustness means here:
#'   Under model misspecification, the sandwich describes sampling variation
#'   around the working model's limiting parameter (its pseudo-true target).
#'   That target need not equal the generating rater severity or criterion
#'   difficulty. Changing the SE does not correct a biased estimate, informative
#'   assignment, unmodeled population differences or MNAR nonresponse. Check
#'   the model and assignment before interpreting a severity contrast.
#'
#'   These intervals condition on the observed fixed facet levels. They do not
#'   generalize to replacement raters sampled from a rater population. They
#'   are not multiway crossed-cluster, G/D-study, variance-boundary, EAP or
#'   multiple-imputation intervals. Shared random rater/task effects spanning
#'   the declared clusters violate this one-way independence assumption.
#'   Numerical success is not a general coverage or rater-diagnosis guarantee.
#'   Review quadrature sensitivity separately; the helper reuses the fitted grid.
#'
#' @section Bounded evaluation:
#'   For the RSM/PCM route, a 1,600-dataset comparison used 80/320 independent persons, three
#'   fixed raters and two criteria. In its combined skewed-ability/sparse-design
#'   scenario, sandwich coverage of generating contrasts ranged from 87.0 to
#'   95.5 percent despite all intervals being available. Coverage of the
#'   independently calculated working-model targets ranged from 93.0 to 97.5
#'   percent. These are ranges across conditions/contrasts, not uncertainty
#'   bounds or universal operating characteristics. At 200 datasets per
#'   condition, MCSE near 95 percent is about 1.54 percentage points.
#'   See \code{vignette("mfrmr-facet-intervals", package = "mfrmr")} for the
#'   design, interpretation and a complete rater-feedback example.
#'
#' @references White, H. (1982). Maximum likelihood estimation of misspecified
#'   models. *Econometrica*, 50(1), 1--25.
#'   Zeileis, A. (2006). Object-oriented computation of sandwich
#'   estimators. *Journal of Statistical Software*, 16(9), 1--16.
#'   \doi{10.18637/jss.v016.i09}.
#'   Zeileis, A., Koell, S. and Graham, N. (2020). Various versatile variances:
#'   An object-oriented implementation of clustered covariances in R.
#'   *Journal of Statistical Software*, 95(1), 1--36.
#'   \doi{10.18637/jss.v095.i01}.
#' @section Save, display and report the selected intervals:
#'   Use [plot()], [as_ggplot()] and [plot_data()] to display or extract the
#'   saved result, and [apa_table()] for tables. Set `title = NULL`,
#'   `subtitle = NULL` or `caption = NULL` in the plot to omit that text.
#'   For supported RSM/PCM MML fits, attach one result or a named list to
#'   [mfrm_results()], for example
#'   `mfrm_results(fit, intervals = list(raters = intervals),
#'   include = c("fit", "plots"), compute = "never")`.
#'   The route `plot(res, type = "facet_raters")` shows the selected intervals.
#'   For two-family GPCM, the corresponding plot type is `"gpcm_raters"`.
#'   [mfrm_report()] and [export_mfrm_results()] retain the method, level,
#'   contrast coefficients, cluster mapping and unavailable outcomes.
#'   The source fit must match; replay reloads the saved results without
#'   refitting or recomputing covariance. These intervals do not replace
#'   uncertainty in ordinary Wright maps, fit diagnostics or Person scores.
#'   For CML, keep the interval object itself: use `summary(intervals)`,
#'   `plot(intervals)`, `apa_table(intervals)` and `saveRDS(intervals, file)`.
#'   Reopen it with `readRDS(file)`. Ordinary [mfrm_results()], [mfrm_report()]
#'   and [export_mfrm_results()] do not accept CML fits or their intervals.
#'
#' @seealso [analyze_facet_equivalence()], [pool_mfrm_imputed()],
#'   [plot.mfrm_facet_intervals()], [mfrm_population_intervals()]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("mfrm_facet_intervals", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   model = "RSM", method = "MML", rating_min = 1, rating_max = 4,
#'   category_policy = "preserve")
#' intervals <- mfrm_facet_intervals(fit, "Rater", method = "sandwich")
#' summary(intervals)
#' plot(intervals)
#'
#' # Compare two selected raters, rather than each rater against the origin.
#' # Here choose the first two example labels; prespecify your own comparison.
#' raters <- unique(as.character(ratings$Rater))
#' contrast <- matrix(0, 1, length(raters), dimnames = list(
#'   paste(raters[1], "minus", raters[2]), raters))
#' contrast[1, 1:2] <- c(1, -1)
#' difference <- mfrm_facet_intervals(fit, "Rater", contrasts = contrast)
#' summary(difference)[, c("Target", "Estimate", "SE", "Lower", "Upper", "Status")]
#' # With the default negative rater sign, a positive contrast means the first
#' # rater is stricter. Units are native logits. An interval containing zero
#' # leaves the direction uncertain; it does not establish equivalence.
#'
#' # Estimate the population SD in this same model before standardizing.
#' # Even an intercept-only population model needs one row per person.
#' persons <- unique(ratings["Person"])
#' population_fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"),
#'   "Score", model = "RSM", method = "MML", population_formula = ~ 1,
#'   person_data = persons, person_id = "Person",
#'   rating_min = 1, rating_max = 4, category_policy = "preserve", quad_points = 61)
#' standardized <- mfrm_facet_intervals(population_fit, "Rater",
#'   contrasts = contrast, scale = "standardized")
#' summary(standardized)
#' # Estimate is now in fitted ability-population SD units; its delta-method
#' # interval includes uncertainty in that SD. It is not a Wald z statistic.
#'
#' # CML estimates the native contrast without an ability distribution.
#' if (requireNamespace("lpSolve", quietly = TRUE)) {
#'   cml_fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'     model = "RSM", method = "CML", rating_min = 1, rating_max = 4,
#'     category_policy = "preserve")
#'   cml_ci <- mfrm_facet_intervals(cml_fit, "Rater", contrasts = contrast)
#'   summary(cml_ci)
#'   apa_table(cml_ci)
#'   file <- tempfile(fileext = ".rds") # Choose a lasting path for your analysis
#'   saveRDS(cml_ci, file)
#'   saved <- readRDS(file)
#'   summary(saved) # Same conditional-information result, without refitting
#'   # Use plot(saved) for a figure; ordinary mfrm_results() does not accept CML.
#' }
#' }
#' @export
#' @order 1
mfrm_facet_intervals <- function(fit, facet = NULL, contrasts = NULL,
                                 method = c("model", "sandwich"),
                                 clusters = NULL, adjust = FALSE, level = .95,
                                 scale = c("native", "standardized"),
                                 comparison_review = NULL,
                                 anchor_covariance = NULL, independent_anchors = FALSE,
                                 interaction = NULL) {
  joint <- !is.null(interaction)
  if (joint && (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_comparison_fit") ||
      !fit$config$model %in% c("RSM", "PCM") || !identical(fit$config$method, "MML")))
    stop("Joint cell intervals require a full native RSM/PCM MML fit.", call. = FALSE)
  if (joint && (!is.null(facet) || !is.null(comparison_review) ||
      match.arg(method) != "model" || match.arg(scale) != "native" ||
      !is.null(anchor_covariance) || !identical(independent_anchors, FALSE)))
    stop("Use interaction with method = 'model', scale = 'native', and omit facet, comparison_review and anchor-source arguments.", call. = FALSE)
  anchor_request <- mfrmr_anchor_uncertainty_request(fit, anchor_covariance, independent_anchors)
  if (!is.null(anchor_request) && (match.arg(method) != "model" || !is.null(comparison_review)))
    stop("Anchor-source intervals require method = 'model' without comparison_review.", call. = FALSE)
  if (inherits(fit, "mfrm_comparison_fit") && identical(fit$config$method, "CML")) return(mfrmr_cml_comparison_intervals(fit,
    facet, contrasts, match.arg(method), clusters, adjust, level, match.arg(scale), comparison_review))
  comparison <- !is.null(comparison_review)
  if (comparison && (!inherits(fit, "mfrm_fit") || !identical(fit$config$method, "MML") ||
      !fit$config$model %in% c("RSM", "PCM")))
    stop("Interaction intervals require a full native RSM/PCM MML fit; CML estimates no population SD and comparison-only MML fits have no interval route here.", call. = FALSE)
  if (inherits(fit, "mfrm_cml_fit")) return(mfrmr_cml_facet_intervals(fit, facet,
    contrasts, match.arg(method), clusters, adjust, level, match.arg(scale)))
  stop_if_jml_adjustment(fit, "mfrm_facet_intervals()")
  product <- mfrm_has_product_slopes(fit)
  native <- !product && (joint || comparison || isTRUE(fit$config$population_spec$active) ||
    identical(fit$config$model, "GPCM"))
  method <- match.arg(method)
  singleton <- !product && !native && identical(method, "model") &&
    inherits(fit, "mfrm_fit") && identical(fit$config$method, "MML") &&
    fit$config$model %in% c("RSM", "PCM") && !mfrm_inference_ready(fit) &&
    identical(as.character(mfrmr_get_readiness_record(fit)$fit$CategoryState[1]), "weak_information")
  standardized <- match.arg(scale) == "standardized"
  if (standardized && (method != "model" || (!comparison && is.null(contrasts))))
    stop("Standardized intervals require explicit zero-sum contrasts and method = 'model'.", call. = FALSE)
  population_active <- isTRUE(fit$config$population_spec$active)
  population_estimated <- population_active && is.null(fit$config$population_spec$fixed_sd)
  if (standardized && population_active &&
      !identical(fit$config$population_spec$design_columns, "(Intercept)"))
    stop("Standardized contrasts require an intercept-only normal population; regression residual SD is not marginal population SD.", call. = FALSE)
  if (!product && !native && (!inherits(fit, "mfrm_fit") || !identical(fit$config$method, "MML") ||
      !fit$config$model %in% c("RSM", "PCM") ||
      isTRUE(fit$config$population_spec$active) || mfrmr_adaptive_integration(fit$config) ||
      any(!is.finite(fit$prep$data$Weight)) || any(fit$prep$data$Weight != 1) ||
      (!mfrm_inference_ready(fit) && !singleton))) {
    stop("Use an inference-ready RSM/PCM MML fit with unit weights, fixed quadrature and a fixed standard-normal person distribution.", call. = FALSE)
  }
  if (!is.logical(adjust) || length(adjust) != 1L || is.na(adjust) ||
      !is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("Supply a logical `adjust` and 0 < level < 1.", call. = FALSE)
  }
  level <- unname(level)
  if (method == "model" && (!is.null(clusters) || adjust)) {
    stop("`clusters` and `adjust` apply only to method = 'sandwich'.", call. = FALSE)
  }
  if (product && (method != "model" || !inherits(fit, "mfrm_fit") ||
      inherits(fit, "mfrm_imported_fit") || !identical(fit$config$method, "MML"))) {
    stop("Two-family location intervals require native MML and method = 'model'; sandwich intervals are unavailable.", call. = FALSE)
  }
  if (native && method != "model") {
    stop("Native population/one-family GPCM location intervals require method = 'model'; sandwich intervals are unavailable.", call. = FALSE)
  }
  params <- if (native) mfrm_native_location_source(fit, comparison_review, interaction) else NULL
  sizes <- build_param_sizes(fit$config)
  if (!native) params <- expand_params(fit$opt$par, sizes, fit$config)
  if (joint) {
    facet <- interaction
    target <- mfrm_joint_interaction_targets(fit, interaction, contrasts)
    estimate <- drop(target$jacobian %*% fit$opt$par)
  } else if (comparison) {
    review <- mfrmr_match_comparison_model(fit$prep, fit$config, comparison_review)
    s <- review$specification
    if (!is.null(contrasts) || (!is.null(facet) && !identical(facet, s$facet)))
      stop("Use the comparison review's contrast weights and facet; omit contrasts and any different facet.", call. = FALSE)
    facet <- paste(s$facet, s$group, sep = ":")
    jacobian <- matrix(0, nrow(review$weights), length(fit$opt$par))
    jacobian[, review$parameter_map$OptimizerIndex] <- review$target_map
    target <- list(contrasts = review$weights, jacobian = jacobian,
      fixed = rowSums(abs(jacobian)) == 0, constant = rep(0, nrow(jacobian)))
    estimate <- drop(jacobian %*% fit$opt$par)
  } else {
    target <- mfrm_facet_contrasts(fit, facet, contrasts)
    estimate <- drop(target$contrasts %*% params$facets[[facet]])
    estimate[target$fixed] <- target$constant[target$fixed]
  }
  if (any(!is.finite(estimate)))
    stop("Contrast estimates must be finite; rescale excessively large coefficients.", call. = FALSE)
  slice <- if (comparison || joint) seq_along(fit$opt$par) else build_param_slices(sizes)[[facet]]
  if (standardized) {
    # Scale rows only for the zero-sum check, avoiding overflow without changing C.
    unit_rows <- target$contrasts / apply(abs(target$contrasts), 1L, max)
    if (any(abs(rowSums(unit_rows)) > 64 * .Machine$double.eps * rowSums(abs(unit_rows))))
      stop("Standardized contrasts must have zero-sum coefficient rows.", call. = FALSE)
    population_sd <- if (population_active) sqrt(params$population$sigma2) else 1
    if (!is.finite(population_sd) || population_sd <= 0)
      stop("Standardized contrasts require a finite positive population SD.", call. = FALSE)
    native_estimate <- estimate
    estimate <- estimate / population_sd
    target_jacobian <- matrix(0, nrow(target$contrasts), length(fit$opt$par),
      dimnames = list(rownames(target$contrasts), mfrm_checkpoint_parameter_names(sizes)))
    target_jacobian[, slice] <- target$jacobian / population_sd
    if (population_estimated)
      target_jacobian[, build_param_slices(sizes)[["log_sigma2"]]] <- -estimate / 2
    if (any(!is.finite(estimate)) || any(!is.finite(target_jacobian)))
      stop("Standardized estimates and derivatives must be finite; rescale excessively large coefficients.", call. = FALSE)
    target$fixed <- rowSums(abs(target_jacobian)) == 0
  }
  information_reason <- NULL
  if (method == "sandwich") clusters <- mfrm_sandwich_clusters(fit$prep$levels$Person, clusters)
  if (product) {
    expected <- build_other_facet_table(fit$config, fit$prep, params)
    keys <- c("Facet", "Level", "Estimate")
    if (!isTRUE(all.equal(fit$facets$others[keys], expected[keys],
        tolerance = 1e-10, check.attributes = FALSE)))
      stop("Saved locations must match the fitted two-family parameters.", call. = FALSE)
    inference <- mfrm_gpcm_product_inference(fit)
    information <- inference$information
    model_cov <- inference$covariance
    if (is.null(model_cov)) model_cov <- matrix(NA_real_, length(fit$opt$par), length(fit$opt$par))
  } else if (native || singleton) {
    inference <- mfrm_native_location_inference(fit)
    information <- inference$information
    model_cov <- inference$covariance
  } else {
    information <- compute_mml_parameter_covariance(fit)
    model_cov <- information$cov
    if (!identical(information$status, "ok") || isTRUE(information$regularized) ||
        !is.matrix(model_cov) || !is.numeric(model_cov) || is.complex(model_cov) ||
        !identical(dim(model_cov), rep(length(fit$opt$par), 2L)) ||
        any(!is.finite(model_cov)) ||
        is.null(tryCatch(chol(model_cov), error = function(e) NULL))) {
      information_reason <- paste("Unregularized positive-definite joint covariance is unavailable.",
        paste(information$detail, collapse = " "),
        "Point estimates are retained; model and sandwich uncertainty are withheld. Information failure can also limit confidence in point accuracy along weak directions.")
      # Do not retain a rejected inverse or allocate a dense placeholder after
      # an information-workspace refusal. Target tables still carry NA values.
      model_cov <- NULL
    }
  }
  cov <- model_cov
  person_scores <- cluster_scores <- NULL
  nclusters <- rank <- NA_integer_
  factor <- 1
  if (method == "sandwich" && is.null(information_reason)) {
    person_scores <- mfrm_person_likelihood_scores(fit)
    sandwich <- mfrm_cluster_sandwich(model_cov, person_scores, clusters, adjust)
    clusters <- sandwich$clusters; cluster_scores <- sandwich$cluster_scores
    nclusters <- nrow(cluster_scores); rank <- sandwich$rank
    factor <- sandwich$factor; cov <- sandwich$covariance
  } else if (method == "sandwich") {
    nclusters <- length(unique(clusters$Cluster))
    factor <- if (adjust) nclusters / (nclusters - 1) else 1
  }
  transform <- function(v) {
    z <- if (is.null(v)) matrix(NA_real_, nrow(target$contrasts), nrow(target$contrasts),
      dimnames = rep(list(rownames(target$contrasts)), 2L)) else if (standardized)
      symmetrize_matrix(target_jacobian %*% v %*% t(target_jacobian)) else
      symmetrize_matrix(target$jacobian %*% v[slice, slice, drop = FALSE] %*% t(target$jacobian))
    z[target$fixed, ] <- 0; z[, target$fixed] <- 0
    z
  }
  selected_cov <- transform(cov)
  model_target_cov <- transform(model_cov)
  se <- covariance_diag_se(selected_cov)
  model_se <- covariance_diag_se(model_target_cov)
  status <- rep("available", length(estimate))
  if ((product || native || singleton) && !isTRUE(inference$check$eligible)) status[] <- "numerical_review_failed"
  if (method == "sandwich" && !is.na(rank) && rank < ncol(cluster_scores)) status[] <- "insufficient_cluster_rank"
  status[status != "numerical_review_failed" & (!is.finite(se) | se <= 0)] <- "nonpositive_variance"
  if (!is.null(information_reason)) status[] <- "information_unavailable"
  status[target$fixed] <- "fixed"
  se[!status %in% c("available", "fixed")] <- NA_real_
  critical <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  margin <- critical * se
  margin[target$fixed] <- NA_real_
  model_margin <- critical * model_se
  model_margin[target$fixed] <- NA_real_
  tab <- data.frame(Target = rownames(target$contrasts), Estimate = estimate,
    SE = se, Lower = estimate - margin, Upper = estimate + margin,
    ModelSE = model_se, ModelLower = estimate - model_margin,
    ModelUpper = estimate + model_margin, Status = status, row.names = NULL)
  tab$CIEligible <- status == "available"
  tab$InferenceReview <- "Pointwise normal approximation under the selected covariance method; numerical availability does not guarantee coverage."
  tab$InferenceReview[status == "fixed"] <- "Constraint-fixed target; SE is zero by construction and no inferential interval is supplied."
  tab$InferenceReview[status == "insufficient_cluster_rank"] <- paste0(
    "Independent cluster scores do not span all free parameters (rank ", rank,
    " of ", length(fit$opt$par), "). Sandwich intervals are withheld; the point and separate model intervals remain. Review independent sampling units and model dimension; do not split dependent clusters.")
  tab$InferenceReview[status == "nonpositive_variance"] <- "The target variance is nonpositive or nonfinite; the point remains but the selected interval is unavailable. Inspect the covariance and contrast scaling."
  if (!is.null(information_reason)) tab$InferenceReview[status == "information_unavailable"] <- information_reason
  if (native || product || standardized || singleton) {
    invalid <- status == "available" & (!is.finite(tab$Lower) | !is.finite(tab$Upper))
    status[invalid] <- "nonfinite_bounds"
    tab$Status <- status
    tab[!status %in% c("available", "fixed"),
      c("SE", "Lower", "Upper", "ModelSE", "ModelLower", "ModelUpper")] <- NA_real_
  }
  out <- list(table = tab, covariance = selected_cov, model_covariance = model_target_cov,
    parameter_covariance = cov, model_parameter_covariance = model_cov,
    person_scores = person_scores, cluster_scores = cluster_scores, clusters = clusters,
    contrasts = target$contrasts,
    settings = list(facet = facet, method = method, level = level, adjust = adjust,
      adjustment_factor = factor, persons = fit$prep$n_person, clusters = nclusters,
      cluster_score_rank = rank, free_parameters = length(fit$opt$par),
      reference = "pointwise normal", quad_points = fit$config$estimation_control$quad_points),
    fit = fit)
  class(out) <- "mfrm_facet_intervals"
  out$cautions <- unique(c(mfrm_mml_information_caution(information),
    mfrm_single_person_caution(fit$prep)))
  if (fit$prep$n_person == 1L) out$table$CoverageQualified <- FALSE
  out$information_review <- information$solution_information$inverse_review
  if (!is.null(information_reason)) {
    out$unavailable_reason <- information_reason
    out$checks <- data.frame(Check = "Joint information", Passed = FALSE, Detail = information_reason)
    out$cautions <- unique(c(out$cautions, information_reason))
  }
  if (singleton) {
    out$checks <- inference$checks
    out$numerical_checks <- inference$numerical_checks
    out$comparison_information_review <- inference$comparison_information_review
    out$table$CoverageQualified <- FALSE
    out$table$InferenceReview[status == "available"] <-
      "Pointwise model-information normal approximation; singleton categories retained after joint-information and finer-quadrature checks. Small-sample coverage is not guaranteed."
    out$table$InferenceReview[status == "numerical_review_failed"] <- inference$check$review
    out$cautions <- unique(c(out$cautions, inference$check$caution))
  }
  if (product) {
    out$source <- mfrm_gpcm_inference_source(fit)
    out$settings$two_family <- TRUE
    out$settings$scale <- "fixed_standard_normal"
    out$settings$location_reference <- if (facet == fit$config$slope_facet[1L])
      "sum_to_zero_across_levels" else "uncentered_on_fixed_N01"
    out$settings$integration <- if (mfrmr_adaptive_integration(fit$config)) "adaptive" else "fixed"
    out$checks <- inference$checks
    out$numerical_checks <- inference$numerical_checks
    out$score_rank <- inference$score_rank
    out$table$CIEligible <- status == "available"
    out$table$CoverageQualified <- FALSE
    out$table$ScaleReference <- out$settings$location_reference
    out$table$InferenceReview <- inference$check$review
    out$table$InferenceReview[status == "fixed"] <- "Constraint-fixed native target; no inferential interval."
    out$table$InferenceReview[status %in% c("nonpositive_variance", "nonfinite_bounds")] <-
      "The target requires finite positive variance and finite normal bounds."
    out$cautions <- unique(c(out$cautions, inference$check$caution,
      "Experimental two-family location intervals: sampling coverage is not qualified. Location differences do not imply uniform differences in expected ratings when slopes or category steps differ."))
  }
  if (native) {
    out$settings$procedure <- "mml_native_location_model_v1"
    out$settings$scale <- "native_ability"
    out$settings$location_reference <- if (!is.null(fit$config$population_spec$fixed_sd) &&
        identical(facet, fit$config$noncenter_facet)) "direct_anchor_origin" else
      "sum_to_zero_across_levels"
    out$settings$integration <- if (mfrmr_adaptive_integration(fit$config)) "adaptive" else "fixed"
    out$settings$comparison_quad_points <- 2 * out$settings$quad_points - 1
    out$settings$population <- if (!isTRUE(fit$config$population_spec$active)) "fixed_standard_normal" else
      if (identical(fit$config$population_spec$design_columns, "(Intercept)"))
        "estimated_intercept_normal" else "conditional_normal_regression"
    out$settings$slope_owner <- fit$config$slope_facet %||% NA_character_
    out$settings$step_owner <- fit$config$step_facet %||% "shared"
    out$source <- mfrm_native_location_identity(fit)
    out$checks <- inference$checks
    out$numerical_checks <- inference$numerical_checks
    out$comparison_information_review <- inference$comparison_information_review
    out$target_jacobian <- matrix(0, nrow(target$contrasts), length(fit$opt$par),
      dimnames = list(rownames(target$contrasts), mfrm_checkpoint_parameter_names(sizes)))
    out$target_jacobian[, slice] <- target$jacobian
    out$population <- if (out$settings$population == "conditional_normal_regression")
      data.frame(Term = fit$config$population_spec$design_columns,
        Coefficient = unname(params$population$coefficients),
        ResidualVariance = params$population$sigma2) else
      data.frame(Mean = if (is.null(params$population)) 0 else
        unname(params$population$coefficients), Variance = params$population$sigma2 %||% 1)
    out$table$CIEligible <- status == "available"
    out$table$ScaleReference <- if (out$settings$location_reference == "direct_anchor_origin")
      "native ability units; supplied anchor origin" else
      "native ability units; centered facet locations"
    out$table$InferenceReview <- inference$check$review
    out$table$InferenceReview[status == "fixed"] <- "Constraint-fixed native target; no inferential interval."
    out$table$InferenceReview[status %in% c("nonpositive_variance", "nonfinite_bounds")] <-
      "The target requires finite positive variance and finite normal bounds."
    out$cautions <- unique(c(out$cautions, inference$check$caution,
      mfrm_gpcm_anchor_note(fit),
      if (!is.null(fit$config$population_spec$fixed_sd))
        "Intervals condition on the supplied population SD; its uncertainty is not propagated.",
      if (!standardized) "Native location intervals: pointwise Wald inference assumes a correctly specified, identified regular MML model; numerical eligibility does not guarantee finite-sample coverage. Native ability units are not standardized by the fitted population SD; unequal slopes or steps can make expected-rating differences vary with ability."))
  }
  if (standardized) {
    out$settings$scale <- "standardized_population_sd"
    out$settings$transformation <- "population_sd_delta_v1"
    out$settings$population_sd <- population_sd
    out$settings$population_sd_estimated <- population_estimated
    out$settings$location_reference <- "zero_sum_contrast_in_population_sd_units"
    out$target_jacobian <- target_jacobian
    out$table$NativeEstimate <- native_estimate
    out$table$CIEligible <- status == "available"
    out$table$CoverageQualified <- FALSE
    out$table$ScaleReference <- "within-facet zero-sum contrast / population SD"
    if (is.null(out$table$InferenceReview)) out$table$InferenceReview <- ""
    out$table$InferenceReview[status == "available"] <-
      "Population-SD-standardized pointwise delta-method interval; source numerical checks passed, sampling coverage is unqualified."
    out$table$InferenceReview[status == "fixed"] <- "Constraint-fixed standardized target; no inferential interval."
    out$table$InferenceReview[status %in% c("nonpositive_variance", "nonfinite_bounds")] <-
      "The standardized target requires finite positive variance and finite normal bounds."
    out$cautions <- unique(c(out$cautions,
      "Population-SD-standardized contrasts use full joint delta-method uncertainty; sampling coverage is not qualified. Regular inference requires a correctly specified, identified marginal model and a true population SD strictly above zero. A positive fitted SD alone does not establish these conditions; near-zero SD can make the ratio unstable. Quadrature agreement is a numerical diagnostic, not an error bound. Unequal slopes or steps can make expected-rating differences vary with ability."))
  }
  if (comparison) {
    out$settings$procedure <- "mml_interaction_comparison_model_v1"
    out$settings$target <- "facet_by_origin_interaction_contrast"
    out$settings$location_reference <- if (standardized) "interaction_contrast_in_population_sd_units" else
      "interaction_contrast_in_native_ability_units"
    out$comparison_review <- review
    out$table$Group1 <- s$groups[1L]; out$table$Group2 <- s$groups[2L]
    out$table$ConditionalState <- review$table$ConditionalState
    out$table$MarginalIdentificationState <- review$table$State
    out$table$IdentificationBasis <- review$table$IdentificationBasis
    out$table$IdentificationDetail <- review$table$IdentificationDetail
    out$table$CoverageQualified <- FALSE
    out$table$ScaleReference <- if (standardized) "facet-by-origin interaction contrast / population SD" else
      "facet-by-origin interaction contrast in native ability units"
    out$cautions <- unique(c(out$cautions,
      "Interaction contrasts use the declared group direction and panel weights. Positive conditional rank and local information do not prove global marginal identification; intervals are model-based asymptotic approximations, not calibrated DFF classifications."))
  }
  if (joint) {
    out$settings$procedure <- "mml_joint_interaction_model_v1"
    out$settings$target <- if (is.null(contrasts)) "joint_interaction_cells" else "joint_interaction_contrasts"
    out$settings$location_reference <- "two_way_sum_to_zero_margins"
    out$interaction_cells <- target$cells
    out$table$CoverageQualified <- FALSE
    out$table$ScaleReference <- "additive predictor shift; fixed N(0,1) ability units"
    if (is.null(contrasts)) out$table <- cbind(out$table,
      target$cells[setdiff(names(target$cells), c("Target", "Estimate"))])
    out$cautions <- unique(c(out$cautions,
      "Jointly fitted interaction effects use zero row and column sums across the declared facet levels. Positive values increase the score predictor relative to the jointly fitted main effects. These are not fixed-baseline estimate_bias() screens or proof of rater unfairness. Intervals are pointwise, not simultaneous; no multiplicity-adjusted decision is supplied. Sparse or missing cells require design and numerical review; finite estimates do not establish coverage."))
  }
  out <- mfrmr_anchor_uncertainty_intervals(out, anchor_request)
  if (length(out$cautions)) {
    out$table$InferenceCaution <- paste(out$cautions, collapse = " ")
    warning(paste(out$cautions, collapse = " "), call. = FALSE)
  }
  out
}

mfrm_facet_contrasts <- function(fit, facet, contrasts) {
  if (!is.character(facet) || length(facet) != 1L || is.na(facet) ||
      !facet %in% fit$config$facet_names) {
    stop("`facet` must name a non-person model facet.", call. = FALSE)
  }
  spec <- fit$config$facet_specs[[facet]]
  labels <- as.character(spec$levels)
  if (is.null(contrasts)) {
    contrasts <- diag(length(labels))
    dimnames(contrasts) <- list(labels, labels)
  }
  if (!is.matrix(contrasts) || !is.numeric(contrasts) || is.complex(contrasts) ||
      !nrow(contrasts) || ncol(contrasts) != length(labels) || any(!is.finite(contrasts)) ||
      is.null(colnames(contrasts)) || anyDuplicated(colnames(contrasts)) ||
      !setequal(colnames(contrasts), labels) || is.null(rownames(contrasts)) ||
      anyNA(rownames(contrasts)) || any(!nzchar(rownames(contrasts))) ||
      anyDuplicated(rownames(contrasts)) || any(rowSums(abs(contrasts)) == 0)) {
    stop("`contrasts` must have distinct target row names, all facet-level column names, and finite nonzero coefficient rows.", call. = FALSE)
  }
  contrasts <- contrasts[, labels, drop = FALSE]
  jac <- contrasts %*% constraint_jacobian(spec)
  list(labels = labels, contrasts = contrasts, jacobian = jac,
    fixed = if (ncol(jac)) rowSums(abs(jac)) == 0 else rep(TRUE, nrow(jac)),
    constant = drop(contrasts %*% expand_facet_with_constraints(numeric(spec$n_params), spec)))
}

# Reuse posterior moments and the same constraint projections as the likelihood
# gradient. The loop is over Persons after one joint quadrature calculation.
mfrm_person_likelihood_scores <- function(fit) {
  config <- fit$config
  sizes <- build_param_sizes(config)
  params <- expand_params(fit$opt$par, sizes, config)
  idx <- build_indices(fit$prep, step_facet = config$step_facet,
    interaction_specs = config$interaction_specs)
  quad <- gauss_hermite_normal(config$estimation_control$quad_points)
  bundle <- mfrm_mml_logprob_bundle(idx, config, quad, params,
    compute_base_eta(idx, params, config), include_probs = TRUE)
  posterior <- mfrm_mml_posterior_bundle(bundle)
  expected <- mfrm_mml_expected_category_bundle(bundle, posterior, include_p_geq = TRUE)
  residual <- idx$score_k - expected$expected_k
  step_residual <- expected$p_geq - outer(idx$score_k, seq_len(config$n_cat - 1L), ">=")
  rows <- split(seq_along(idx$person), factor(idx$person, levels = seq_len(config$n_person)))
  scores <- t(vapply(rows, function(at) {
    facet_scores <- unlist(lapply(config$facet_names, function(facet) {
      expanded <- numeric(length(params$facets[[facet]]))
      sums <- rowsum(matrix(residual[at], ncol = 1), idx$facets[[facet]][at], reorder = FALSE)
      expanded[as.integer(rownames(sums))] <- sums[, 1] * config$facet_signs[[facet]]
      constraint_grad_project(expanded, config$facet_specs[[facet]])
    }))
    sub_idx <- list(interactions = lapply(idx$interactions, `[`, at))
    interaction_scores <- compute_interaction_gradient_free(residual[at], sub_idx, config)
    if (config$model == "RSM") {
      step_scores <- colSums(step_residual[at, , drop = FALSE])
    } else {
      step_scores <- matrix(0, nrow(params$steps_mat), ncol(params$steps_mat))
      sums <- rowsum(step_residual[at, , drop = FALSE], idx$step_idx[at], reorder = FALSE)
      step_scores[as.integer(rownames(sums)), ] <- sums
    }
    c(facet_scores, interaction_scores, project_typed_step_gradient(step_scores, config))
  }, numeric(length(fit$opt$par))))
  dimnames(scores) <- list(fit$prep$levels$Person, mfrm_checkpoint_parameter_names(sizes))
  aggregate_gradient <- mfrm_grad_mml(fit$opt$par, idx, config, sizes, quad)
  if (any(!is.finite(scores)) ||
      max(abs(colSums(scores) + aggregate_gradient)) > 1e-8 * max(1, max(colSums(abs(scores))))) {
    stop("Person likelihood scores do not reproduce the fitted objective gradient.", call. = FALSE)
  }
  scores
}

#' @rdname mfrm_facet_intervals
#' @param object An object returned by `mfrm_facet_intervals()`.
#' @param ... For `print.mfrm_facet_intervals()`: Unused by print and summary methods.
#'
#' For `plot()`: Unused.
#' @export
print.mfrm_facet_intervals <- function(x, ...) {
  cat(x$settings$facet, "intervals:", x$settings$method, "covariance\n")
  if (x$settings$method == "sandwich") cat(x$settings$clusters,
    "declared independent clusters;", x$settings$persons, "persons\n")
  print(x$table[setdiff(names(x$table), "InferenceCaution")], row.names = FALSE)
  for (caution in x$cautions) print_wrapped_line(paste0("Caution: ", caution))
  cat("Pointwise normal intervals; changing covariance does not correct biased estimates.\n")
  invisible(x)
}

#' @rdname mfrm_facet_intervals
#' @name summary.mfrm_facet_intervals
#' @aliases summary.mfrm_facet_intervals
#' @usage \method{summary}{mfrm_facet_intervals}(object, ...)
#' @rawNamespace S3method(summary,mfrm_facet_intervals,summary.mfrm_pooled)
NULL

#' @rdname mfrm_facet_intervals
#' @param parm For `confint()`, optional distinct target names from
#'   `object$table$Target`. Extracts the stored limits without recalculation;
#'   `level` must match the saved level. Use `mfrm_facet_intervals()` to
#'   calculate a different level. Unavailable limits remain `NA`; inspect
#'   `summary(object)` for their eligibility and reasons.
#' @export
confint.mfrm_facet_intervals <- function(object, parm = NULL,
                                        level = object$settings$level, ...) {
  rlang::check_dots_empty()
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level != object$settings$level)
    stop("`level` must match the saved level; use `mfrm_facet_intervals()` to recalculate.", call. = FALSE)
  labels <- object$table$Target
  if (is.null(parm)) parm <- labels
  if (!is.character(parm) || anyNA(parm) || anyDuplicated(parm) || !all(parm %in% labels))
    stop("`parm` must contain distinct saved Target names.", call. = FALSE)
  out <- as.matrix(object$table[match(parm, labels), c("Lower", "Upper"), drop = FALSE])
  rownames(out) <- parm
  out
}

#' @rdname mfrm_facet_intervals
#' @importFrom stats vcov
#' @export
vcov.mfrm_facet_intervals <- function(object, ...) object$covariance

# Tables share the same saved values across APA output, reports and exports.
mfrm_facet_interval_tables <- function(x) {
  tab <- x$table
  tab$Facet <- x$settings$facet
  tab$Method <- x$settings$method
  tab$ConfidenceLevel <- paste0(format(100 * x$settings$level, trim = TRUE), "%")
  tab$Adjustment <- "Pointwise"
  settings <- as.data.frame(x$settings, stringsAsFactors = FALSE)
  settings$ConfidenceLevel <- paste0(format(100 * settings$level, trim = TRUE), "%")
  settings$level <- NULL
  tables <- list(intervals = tab,
    settings = settings,
    contrasts = data.frame(Comparison = rownames(x$contrasts),
      x$contrasts, row.names = NULL, check.names = FALSE))
  if (!is.null(x$interaction_cells)) tables$interaction_cells <- x$interaction_cells
  if (!is.null(x$comparison_review)) {
    tables$comparison_review <- x$comparison_review$table
    tables$comparison_cells <- x$comparison_review$cells
  }
  for (name in c("checks", "numerical_checks", "score_rank", "population", "comparison_information_review")) {
    if (!is.null(x[[name]])) tables[[name]] <- x[[name]]
  }
  tables <- c(tables, mfrmr_anchor_uncertainty_tables(x))
  if (!is.null(x$clusters)) tables$clusters <- x$clusters
  if (!is.null(x$information_review)) tables$information_review <- x$information_review
  if (!is.null(x$target_jacobian)) tables$target_jacobian <- data.frame(
    Target = rownames(x$target_jacobian), x$target_jacobian, row.names = NULL, check.names = FALSE)
  tables
}

mfrm_facet_results_inputs <- function(fit, intervals) {
  if (is.null(intervals)) return(NULL)
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !fit$config$model %in% c("RSM", "PCM")) {
    stop("Saved fixed-facet intervals require their native RSM/PCM fit.", call. = FALSE)
  }
  if (inherits(intervals, "mfrm_facet_intervals")) intervals <- list(inference = intervals)
  if (!is.list(intervals) || !length(intervals) || is.null(names(intervals)) ||
      anyNA(names(intervals)) || anyDuplicated(tolower(names(intervals))) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", names(intervals))) ||
      !all(vapply(intervals, inherits, logical(1), what = "mfrm_facet_intervals"))) {
    stop("`intervals` must be saved fixed-facet intervals or a named list of them; use distinct letter/number/underscore names.", call. = FALSE)
  }
  names(intervals) <- tolower(names(intervals))
  # The existing signature also describes ordinary native RSM/PCM fits.
  source <- mfrm_gpcm_inference_source(fit)
  for (x in intervals) {
    if (is.null(x$fit) || !identical(source, mfrm_gpcm_inference_source(x$fit))) {
      stop("Saved fixed-facet intervals must match the fitted parameters, data, constraints, population and integration settings.", call. = FALSE)
    }
    if ((x$settings$procedure %||% "") %in% c("mml_native_location_model_v1", "mml_interaction_comparison_model_v1", "mml_joint_interaction_model_v1") &&
        (!identical(mfrm_native_location_identity(fit), x$source) ||
          !identical(mfrm_native_location_identity(x$fit), x$source))) {
      stop("Saved native location intervals must match the source coordinates, tables and qualification metadata.", call. = FALSE)
    }
  }
  intervals
}

mfrm_facet_results_attach <- function(out, inputs) {
  if (is.null(inputs)) return(out)
  out$facet_intervals <- inputs
  for (name in names(inputs)) {
    key <- paste0("facet_", name)
    tables <- mfrm_facet_interval_tables(inputs[[name]])
    names(tables) <- paste(key, names(tables), sep = "_")
    out$tables <- c(out$tables, tables)
    out$components[[key]] <- inputs[[name]]
    out$status <- rbind(out$status, mfrm_results_status_row(key, "review",
      "Saved fixed-facet pointwise intervals; method, confidence level, contrasts and unavailable outcomes retained. No automatic rater-quality decision."))
    out$plot_map <- dplyr::bind_rows(out$plot_map, data.frame(Type = key,
      Available = "plots" %in% out$include, RequiredArtifact = FALSE,
      Route = paste0('plot(res, type = "', key, '")'),
      Detail = "Saved fixed-facet uncertainty; no refitting or covariance calculation.",
      InterpretationStatus = "approximate_inference", InterpretationReady = FALSE,
      ReadinessRoute = paste0("res$tables$", names(tables)[1])))
  }
  out$table_index <- mfrm_results_table_index(out$tables)
  out$notes <- unique(c(out$notes,
    unlist(lapply(inputs, `[[`, "cautions"), use.names = FALSE),
    "Fixed-facet pointwise intervals condition on the observed facet levels. Sandwich covariance does not correct biased estimates or imply rater quality."))
  out
}
