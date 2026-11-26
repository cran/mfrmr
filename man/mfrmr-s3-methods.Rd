\name{mfrmr-s3-methods}
\alias{print.mfrm_fair_average}
\alias{mfrmr-s3-methods}
\alias{AIC.mfrm_fit}
\alias{anova.mfrm_fit}
\alias{as.data.frame.mfrm_drift}
\alias{as.data.frame.mfrm_fit}
\alias{BIC.mfrm_fit}
\alias{coef.mfrm_drift}
\alias{coef.mfrm_fit}
\alias{confint.mfrm_drift}
\alias{confint.mfrm_fair_average}
\alias{confint.mfrm_fit}
\alias{confint.mfrm_gtheory}
\alias{confint.mfrm_multivariate_d_study}
\alias{fitted.mfrm_fit}
\alias{logLik.mfrm_fit}
\alias{mfrm_d_study.default}
\alias{mfrm_d_study.mfrm_fit}
\alias{mfrm_d_study.mfrm_gtheory}
\alias{nobs.mfrm_fit}
\alias{predict.mfrm_drift}
\alias{predict.mfrm_fit}
\alias{print.mfrm_agreement}
\alias{print.mfrm_anchor_drift}
\alias{print.mfrm_anchor_review}
\alias{print.mfrm_anchored_fit}
\alias{print.mfrm_bias}
\alias{print.mfrm_bundle}
\alias{print.mfrm_comparison}
\alias{print.mfrm_data_description}
\alias{print.mfrm_dff}
\alias{print.mfrm_diagnostics}
\alias{print.mfrm_drift}
\alias{print.mfrm_fit}
\alias{print.mfrm_gtheory}
\alias{print.mfrm_multivariate_d_comparison}
\alias{print.mfrm_multivariate_d_study}
\alias{print.mfrm_person_fit_indices}
\alias{print.summary.mfrm_bias}
\alias{print.summary.mfrm_data_description}
\alias{print.summary.mfrm_dff}
\alias{print.summary.mfrm_diagnostics}
\alias{print.summary.mfrm_drift}
\alias{print.summary.mfrm_multivariate_d_study}
\alias{print.summary.mfrm_person_fit_indices}
\alias{residuals.mfrm_fit}
\alias{summary.mfrm_anchor_drift}
\alias{summary.mfrm_anchor_review}
\alias{summary.mfrm_anchored_fit}
\alias{summary.mfrm_bias}
\alias{summary.mfrm_bundle}
\alias{summary.mfrm_comparison}
\alias{summary.mfrm_data_description}
\alias{summary.mfrm_dff}
\alias{summary.mfrm_diagnostics}
\alias{summary.mfrm_drift}
\alias{summary.mfrm_fit}
\alias{summary.mfrm_gtheory}
\alias{summary.mfrm_multivariate_d_comparison}
\alias{summary.mfrm_multivariate_d_study}
\alias{summary.mfrm_person_fit_indices}
\alias{vcov.mfrm_drift}
\alias{vcov.mfrm_fit}
\alias{vcov.mfrm_gtheory}
\title{Standard S3 Extractor, Prediction, Inference, and Comparison Methods for mfrmr Objects}
\description{
Standard R S3 methods (\code{predict}, \code{coef}, \code{vcov}, \code{confint}, \code{anova}, \code{fitted}, \code{residuals}, \code{logLik}, \code{AIC}, \code{BIC}, \code{nobs}, \code{summary}, \code{print}, and \code{as.data.frame}) for objects returned by \pkg{mfrmr}. In particular, \code{predict()} evaluates expected scores, category probability matrices, linear predictors, \code{Fair(M)}, or new-examinee \code{WLE} scores under arbitrary user-supplied \code{newdata} conditions.
}
\usage{
\method{predict}{mfrm_fit}(object, newdata = NULL,
  type = c("response", "probs", "prob", "link", "fair", "wle", "eap"),
  person_data = NULL, facet_data = NULL, se.fit = FALSE, ...)
\method{coef}{mfrm_fit}(object,
  type = c("all", "facets", "locations", "steps", "slopes", "person", "persons",
    "explanatory", "interactions", "population", "free"), ...)
\method{vcov}{mfrm_fit}(object,
  component = c("design", "facets", "steps", "slopes", "explanatory", "population"),
  se_type = c("model", "sandwich"), facet = NULL, ...)
\method{confint}{mfrm_fit}(object, parm = NULL, level = 0.95,
  component = c("facets", "steps", "slopes", "person", "explanatory",
    "interactions", "population", "free"),
  se_type = c("model", "real", "sandwich", "profile"), ...)
\method{anova}{mfrm_fit}(object, ..., nested = TRUE)
\method{confint}{mfrm_gtheory}(object, parm = NULL, level = 0.95, ...)
}
\arguments{
\item{nested}{Logical; request nested-model likelihood-ratio inference after eligibility checks.}
  \item{object}{A fitted measurement model or G-Study object.}
  \item{newdata}{Rating conditions containing the fitted base facet columns.
    An explicit numeric \code{theta} column sets ability (\code{Theta} is an alias when it is not a model facet; conflicting columns are rejected); otherwise known person
    IDs use the primary fitted score. Unknown persons use the fitted population
    mean if available, or zero. \code{NULL} uses retained rating rows.}
  \item{type}{For prediction: \code{"response"} returns expected rubric scores,
    \code{"probs"} (alias \code{"prob"}) category probabilities,
    \code{"link"} the unscaled ability-plus-facet predictor,
    \code{"fair"} stored person FairM, and \code{"eap"}/\code{"wle"}
    person score tables. EAP scores are retrieved from the fit with
    \code{newdata = NULL}. New-data WLE requires a Person ID and a finite
    integer score within the fitted category support on each row.
    For \code{coef()}, selects a parameter table; \code{"free"} returns
    the identified free coordinates.}
  \item{person_data,facet_data}{Optional covariate tables for unseen persons or
    explanatory facet levels, using the fitted model's column definitions.}
  \item{se.fit}{For \code{"response"} or \code{"link"}, return a list with
    \code{fit} and conditional structural delta-method \code{se.fit}.}
  \item{component}{Parameter block. \code{vcov(..., component = "design")}
    returns the full free-coordinate covariance, including estimated population
    coordinates. For G-Study intervals, \code{"coefficients"} selects G/Phi.}
  \item{se_type}{\code{"model"} uses inverse information; \code{"sandwich"}
    uses person-clustered covariance for MML (model covariance otherwise).
    \code{"real"} uses fit-inflated SEs where available.
    \code{"profile"} profiles selected free structural coordinates.}
  \item{facet}{Optional facet name for covariance extraction.}
  \item{parm}{Optional parameter names or indices selecting interval rows.}
  \item{level}{Confidence level in \eqn{(0,1)}.}
  \item{...}{Additional fitted models for \code{anova()}; otherwise method-specific arguments.}
}
\details{
\strong{Prediction uncertainty.} Response/link SEs use the full covariance of
free structural coordinates, including cross-facet covariances, estimated steps,
interactions, explanatory coefficients, and GPCM log-slopes as applicable.
Ability is held fixed at the supplied or plug-in value: person-score uncertainty,
uncertainty in a predicted population mean, and future rating variance are excluded.
These are SEs of conditional model predictions, not prediction intervals for
individual ratings. An unusable calibration gives \code{NA} SEs.
\code{"link"} is the unscaled predictor; its gradient excludes steps and slopes.

\strong{Measurement-model intervals.} Model and sandwich intervals are local
normal approximations; slope intervals are formed on log-slopes then exponentiated.
For primary EAP scores the person SE is posterior SD conditional on calibration,
so the normal interval is an approximate conditional posterior summary, not a
calibration-inclusive frequentist confidence interval. Fit-inflated \code{"real"}
intervals are descriptive. Numerically ineligible calibrations are rejected.
For WLE, person intervals use conditional score uncertainty with calibration
parameters treated as estimated constants. They do not propagate calibration
uncertainty. MML sandwich SEs cluster scores by person; they still require
independent persons and enough informative clusters for the target parameter.
A large total number of ratings does not compensate for a facet level observed
on very few persons. Sandwich SEs do not repair biased point estimates or
guarantee small-sample coverage. Inspect interval availability separately from
the coverage properties established for the intended model and assignment.

\strong{Nested likelihood comparisons.} \code{anova()} requires the same retained
responses, person/facet assignments, weights, category support and numerical
integration settings, unpenalized MML or CML, converged eligible fits, and nested
identified parameter spaces. It rejects materially decreasing optimized
likelihoods instead of reporting a zero statistic. A failed check leaves
\code{lrt = NULL} and an unavailable \code{comparison_basis$lrt_status}.
The chi-square reference remains an asymptotic regular-interior approximation;
passing computational checks does not guarantee small-sample calibration.
G-Study covariance and interval assumptions are described in
\code{\link{mfrm_g_study}}.
}
\value{
Standard S3 return values: numeric vectors or matrices for \code{predict}, \code{coef}, \code{vcov}, \code{confint}, \code{fitted}, and \code{residuals}; \code{logLik} / scalar criteria for \code{logLik}, \code{AIC}, \code{BIC}, and \code{nobs}; \code{mfrm_comparison} tables for \code{anova}; and structured summaries for \code{summary} and \code{print}.
}
\examples{
ratings <- load_mfrmr_data("example_core")
fit_rsm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score", model = "RSM")
fit_pcm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
                    model = "PCM", step_facet = "Criterion")

coef(fit_rsm, type = "facets")
confint(fit_rsm, component = "facets", se_type = "real")
vcov(fit_rsm, facet = "Rater")
predict(fit_rsm, newdata = data.frame(theta = 0, Rater = "R01", Criterion = "Content"))
anova(fit_rsm, fit_pcm)
}
