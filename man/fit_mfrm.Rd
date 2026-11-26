\name{fit_mfrm}
\alias{fit_mfrm}
\title{Estimate Many-Facet Ordered-Response Measurement Models (RSM, PCM, and GMFRM)}
\description{
Primary estimation entry point for Rating Scale (\code{"RSM"}), Partial Credit (\code{"PCM"}), and Generalized Many-Facet Rasch Models (GMFRM, selected with \code{model = "GPCM"}) across Marginal (\code{"MML"}), Joint (\code{"JML"}), Jeffreys-penalized Joint (\code{"JML_jeffreys"}), and Conditional (\code{"CML"}) Maximum Likelihood in C++. Supports element and group-mean anchors, fixed or shared \code{step_constraints}, declarative \code{category_support} for mixed binary/polytomous rubrics, explanatory \code{facet_formulas} with element covariates, joint \code{facet_interactions}, empirical Bayes or Laplace \code{facet_shrinkage}, Warm's (1989) \code{WLE} scoring, and automatic pre-estimation connectivity and identifiability checks.
}
\usage{
fit_mfrm(
  data, person = NULL, facets = NULL, score = NULL,
  model = c("RSM", "PCM", "GPCM"),
  method = c("MML", "JML", "JMLE", "JML_jeffreys", "CML"),
  step_facet = NULL, slope_facet = NULL, rating_min = NULL,
  rating_max = NULL, weight = NULL, anchors = NULL, group_anchors = NULL,
  population_formula = NULL, person_data = NULL, facet_interactions = NULL,
  step_constraints = NULL, category_support = NULL, facet_formulas = NULL,
  facet_data = NULL, person_score = c("auto", "EAP", "WLE"),
  warn_connectivity = TRUE, spec = NULL, control = list(), ...)
}
\arguments{
  \item{data}{A \code{data.frame} of observed ratings in long format (one row per observed rating; missing assignments are absent rows rather than zeros).}
  \item{person}{Character string naming the person (examinee/target) identifier column.}
  \item{facets}{Character vector naming one or more non-person facet columns (e.g., \code{c("Rater", "Criterion")}). By default, non-person facets carry negative sign (\code{Sign = -1}) so larger logit estimates represent higher rater severity or higher task/criterion difficulty.}
  \item{score}{Character string naming the integer or numeric rubric score column.
    Convert factor labels to numeric rubric values before using the unified engine.}
  \item{model}{\code{"GMFRM"} is accepted as an alias for \code{"GPCM"}; fitted objects retain the code \code{"GPCM"}. Measurement model family: \code{"RSM"} (shared step thresholds across elements), \code{"PCM"} (facet-specific step thresholds governed by \code{step_facet}), or \code{"GPCM"} (GMFRM with facet-specific discrimination and steps governed by \code{slope_facet} and \code{step_facet}).}
  \item{method}{Estimation method: \code{"MML"} (Gauss-Hermite or Laplace-adaptive marginal ML), \code{"JML"} / \code{"JMLE"} (joint ML), \code{"JML_jeffreys"} (Jeffreys/Firth-penalized JML with finite extreme-score estimates and Block-Schur standard errors), or \code{"CML"} (conditional ML conditioning on person total scores for \code{"RSM"} and \code{"PCM"}).}
  \item{step_facet}{Character string naming the facet in \code{facets} that owns category step thresholds when \code{model} is \code{"PCM"} or \code{"GPCM"}. For a GMFRM with common category steps on a common rubric, use \code{step_facet = "Common"}.}
  \item{slope_facet}{Character vector of length 1 or 2 naming the facet(s) in \code{facets} with element-specific discrimination slopes when \code{model = "GPCM"}.}
  \item{rating_min}{Optional integer lower bound of the rubric scale. Inferred from observed scores if \code{NULL}.}
  \item{rating_max}{Optional integer upper bound of the rubric scale. Inferred from observed scores if \code{NULL}.}
  \item{weight}{Optional character string naming a nonnegative observation case-weight column in \code{data}.}
  \item{anchors}{Optional \code{data.frame} with columns \code{Facet}, \code{Level}, and \code{Anchor} fixing individual facet elements to known logit values (see \code{\link{make_anchor_table}}).}
  \item{group_anchors}{Optional \code{data.frame} with columns \code{Facet}, \code{Level}, \code{Group}, and \code{GroupValue} constraining subset means of facet elements.}
  \item{population_formula}{Optional one-sided latent-regression formula (e.g., \code{~ Group}) for person ability means under \code{method = "MML"}.}
  \item{person_data}{Optional \code{data.frame} of person-level covariates with one row per person ID when \code{population_formula} is supplied.}
  \item{facet_interactions}{Optional list of length-2 character vectors (e.g., \code{list(c("Rater", "Criterion"))}) estimating sum-to-zero 2-way interaction terms jointly with main effects.}
  \item{step_constraints}{Optional \code{data.frame} with columns \code{Level}, \code{From}, \code{To}, and \code{Value} or \code{Group} fixing or tying step transitions across levels while preserving within-ladder sum-to-zero identification.}
  \item{category_support}{Optional \code{data.frame} with columns \code{Level}, \code{Min}, and \code{Max} declaring inclusive per-level category bounds on the original rubric scale for mixed binary/polytomous rubrics.}
  \item{facet_formulas}{Optional named list of one-sided covariate formulas (e.g., \code{list(Rater = ~ Experience)}) decomposing facet locations via explanatory linear predictors.}
  \item{facet_data}{Optional named list of element-level covariate \code{data.frame}s corresponding to \code{facet_formulas}, with one row per element ID in each table.}
  \item{person_score}{Primary person ability estimator stored in \code{fit$person$Estimate}: \code{"auto"} (\code{"EAP"} for \code{MML}, \code{"WLE"} for \code{CML}, \code{"JML"} for \code{JML}), \code{"EAP"}, or \code{"WLE"} (Warm, 1989).}
  \item{warn_connectivity}{Logical; if \code{TRUE} (default), audits multi-facet graph connectivity, design rank, and interior category support prior to optimization and emits a warning if disconnected subsets or rank deficiencies are detected.}
  \item{spec}{Optional named list of advanced model specification settings, including \code{slope_reference} (reference facet constrained to geometric mean 1 when two facets have discrimination parameters), \code{slope_type} (\code{"full_predictor"} or \code{"ability_only"}), \code{gpcm_mml_identification} (\code{"free_population"}, \code{"fixed_standard_normal"}, or \code{"unconstrained_slopes"}), \code{noncenter_facet}, \code{positive_facets}, \code{dummy_facets}, and \code{centering_scope} (\code{"remainder"} or \code{"facet"}).}
  \item{control}{Optional named list of numerical and regularization controls, including \code{quad_points} (default \code{31}), \code{mml_integration} (\code{"fixed"}, \code{"adaptive"}, or \code{"auto"}), \code{facet_shrinkage} (\code{"none"}, \code{"empirical_bayes"}, or \code{"laplace"}), \code{facet_prior_sd}, \code{slope_prior_sd}, \code{xtreme}, \code{whexact}, \code{maxit}, and \code{reltol}.}
  \item{...}{Additional supported specification or control arguments, supplied by name. In particular, \code{missing_codes} recodes specified values before retaining complete rating rows: an explicit vector applies to the score and identifier columns; \code{TRUE} or \code{"default"} applies the default codes (99, 999, -1, N, NA, n/a, period, blank) to scores only.}
}
\value{
An object of class \code{"mfrm_fit"} containing:
\itemize{
  \item \code{person}: \code{data.frame} of person ability estimates (\code{Estimate}, \code{ModelSE}, \code{RealSE}, \code{WLE}, \code{WLE_SE}, \code{FairM}, \code{FairZ}, \code{InfitMnSq}, \code{OutfitMnSq}, \code{ExtremeFlag}).
  \item \code{facets}: \code{data.frame} of non-person facet location estimates, \code{ModelSE}, \code{RealSE}, \code{SandwichSE}, \code{FairM}, and fit statistics.
  \item \code{steps}: \code{data.frame} of Andrich step threshold estimates and \code{ModelSE}.
  \item \code{slopes}: \code{data.frame} of GPCM discrimination estimates (\code{Estimate}, \code{LogEstimate}, \code{ModelSE}, \code{SandwichSE}) when \code{model = "GPCM"}.
  \item \code{population}, \code{connectivity}, \code{summary}: latent distribution summary, connectivity audit, and model fit indices (\code{loglik}, \code{aic}, \code{bic}).
}
Supports \code{summary()}, \code{plot()}, \code{predict()}, \code{coef()}, \code{vcov()}, \code{confint()}, \code{anova()}, \code{logLik()}, \code{AIC()}, \code{BIC()}, \code{fitted()}, and \code{residuals()}.
}
\details{
Use the named main arguments for the data and model, \code{spec} for
advanced model assumptions, and \code{control} for numerical settings.
Additional named arguments in \code{...} are retained for compatibility;
do not specify the same setting in both places.

Non-unit observation weights leave point estimates available but withhold
inferential reporting because the supported weight interpretation has not
been validated. A deficit in the constrained facet-location design after
projecting out centered person-location contrasts is also flagged: MML comparisons then
rely on population assumptions when effective anchors do not resolve the
deficit. This conditional location check is not a proof of global
identifiability for all nonlinear slope or threshold structures.

All estimators return the same \code{mfrm_fit} schema and work with
\code{mfrm_workflow(fit = ...)}. The implementation is selected internally;
there is no engine-selection argument. Refit saved objects from older versions
using their original data and model specification. Relabelling their class
does not convert the underlying parameter and covariance layout.

With \code{mml_integration = "auto"}, fixed quadrature fits at \eqn{q} and
\eqn{2q+1} points are compared. The finer fit is returned. Relative free-parameter
changes, absolute expected-score changes, and log-likelihood changes per rating
are checked against \code{1e-4}, \code{1e-4}, and \code{1e-6}, respectively.
The results and selected order are retained in \code{quadrature_check}.
A failed check withholds inferential reporting; increase \code{quad_points} and
inspect stability. Agreement at these two orders is a sensitivity check, not a
bound on integration error, global-optimum guarantee, or coverage validation.

The Generalized Many-Facet Rasch Model (GMFRM) is a unidimensional
adjacent-category model that estimates discrimination. The name follows Uto
and Ueno (2020), who use GMFRM for an extension of the many-facet Rasch model.
In the ordinary many-facet Rasch model, discrimination is held constant.
In this package, GMFRM covers a many-facet generalized partial credit model
(GPCM; Muraki, 1992) with facet-specific category steps and a rating-scale
form (many-facet RS-GPCM) with common category steps. These descriptions refer
to the fitted step structure. The response equation depends on the step
ownership and slope placement selected for the analysis.

Reports name the facets that carry discrimination. When two facets
carry discrimination, their slopes are multiplied. Report the step facet,
slope placement, scale identification, and estimation method when
describing the fitted model. \code{summary(fit)$model_label} gives the
display label for a fitted model; the model code is \code{"GPCM"}.

\code{Converged} reports the optimizer stopping condition. Separately,
\code{summary$InferenceReady} and \code{readiness$reason_codes} report failed
connectivity, convergence, finite-parameter, or positive-definite information checks.
For unified RSM and PCM MML, the population location terms must be
identifiable together with the constrained facet effects. A rank-deficient
joint design cannot identify these effects separately. Recheck the rating
assignments, fixed reference values and population formula.
For unpenalized RSM/PCM MML, \code{boundary_audit} uses the required dependency \pkg{lpSolve}
to seek a structural direction that improves every observed category probability
without changing person ability. A certified direction implies no finite structural
MLE and makes inference unavailable. Failure to find this sufficient condition is
not a general proof of existence; GPCM boundaries and all population-boundary
directions are not covered by this particular test. An unsuccessful numerical boundary check keeps inference unavailable. These checks do not establish finite-sample coverage.

Person and facet fit statistics, expected averages, FairM/FairZ, and local-bias
diagnostics use the primary person score. Selecting MML \code{person_score = "WLE"}
therefore changes these derived quantities as well as \code{person$Estimate};
the EAP and posterior SD remain available in separate columns.
}
\references{
Muraki, E. (1992). A generalized partial credit model: Application of an EM
algorithm. \emph{Applied Psychological Measurement, 16}(2), 159--176.
\doi{10.1177/014662169201600206}.

Uto, M., and Ueno, M. (2020). A generalized many-facet Rasch model and its
Bayesian estimation using Hamiltonian Monte Carlo. \emph{Behaviormetrika, 47},
469--496. \doi{10.1007/s41237-020-00115-7}.
}
\seealso{
\code{\link{mfrm_analysis}}, \code{\link{check_mfrm_connectivity}}, \code{\link{diagnose_mfrm}}, \code{\link{estimate_bias}}, \code{\link{analyze_dff}}, \code{\link{plot.mfrm_fit}}
}
\examples{
# Load bundled 4-rater x 4-criterion analytic rubric ratings (1 row = 1 rating)
ratings <- load_mfrmr_data("example_core")

# 1. Rating Scale Model (RSM) via MML with Warm's WLE person scoring
fit_rsm <- fit_mfrm(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  model = "RSM",
  method = "MML",
  person_score = "WLE"
)
fit_rsm$facets[, c("Facet", "Level", "Estimate", "RealSE", "FairM", "InfitMnSq")]

# 2. GMFRM with criterion-specific discrimination
fit_gpcm <- fit_mfrm(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  model = "GPCM",
  step_facet = "Criterion",
  slope_facet = "Criterion",
  control = list(quad_points = 21L)
)
fit_gpcm$slopes[, c("Facet", "Level", "Estimate", "ModelSE")]
predict(fit_gpcm, newdata = data.frame(theta = 0, Rater = "R01", Criterion = "Content"))
}
