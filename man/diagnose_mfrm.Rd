\name{diagnose_mfrm}
\alias{diagnose_mfrm}
\title{Compute Element Fit Diagnostics, Reliability, Separation, and Unexpected Cells}
\description{
Summarizes facet- and person-level measurement quality from a fitted Many-Facet model, including \code{Infit} and \code{Outfit} mean-square (\code{MnSq}) and standardized (\code{ZSTD}) fit statistics, \code{ModelSE}- vs \code{RealSE}-based separation, reliability, Strata, homogeneity chi-square tests, and standardized unexpected response cells (\code{|StdResidual| >= 2}).
}
\usage{
diagnose_mfrm(
  fit,
  residual_pca = c("none", "overall", "facet", "both"),
  pca_max_factors = 10L)
}
\arguments{
  \item{fit}{A fitted model object returned by \code{\link{fit_mfrm}}.}
  \item{residual_pca}{Optional residual principal-components mode (\code{"none"}, \code{"overall"}, \code{"facet"}, or \code{"both"}).}
  \item{pca_max_factors}{Maximum number of residual PCA factors extracted when \code{residual_pca != "none"}.}
}
\value{
An S3 diagnostic bundle of class \code{"mfrm_diagnostics"} containing \code{reliability} (per-facet \code{ModelReliability}, \code{RealReliability}, \code{ModelSeparation}, \code{RealSeparation}, \code{Strata}, and \code{HomogeneityChiSq}), \code{measures} (element-level \code{Estimate}, \code{ModelSE}, \code{RealSE}, \code{FairM}, \code{Infit}, \code{InfitZSTD}, \code{Outfit}, \code{OutfitZSTD}), and \code{unexpected} (flagged rating rows sorted by \code{|StdResidual|}).
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{compute_person_fit_indices}}, \code{\link{fair_average_table}}, \code{\link{plot.mfrm_fit}}
}
\details{
Unified diagnostics include \code{readiness}, independently of optimizer
convergence. When unready, homogeneity tests are withheld (\code{NA}) and plots
omit interval bars; descriptive estimates and fit statistics remain available.
Summary/print output and plot notes retain the reasons, including after saving
and reloading the diagnostic object. Numerical readiness is not a guarantee
of interval coverage or fit-screening accuracy.

Residual PCA averages standardized residuals for each Person and selected facet
level or joint context. Results are exploratory. Missing or non-positive
semidefinite pairwise correlations are reported as unavailable. For example,
use \code{plot(diagnose_mfrm(fit, residual_pca = "facet"), facet = "Criterion",
type = "residual_pca")} to see the first contrast.

For primary EAP person scores, model reliability is
\eqn{v/(v+s)}, where \eqn{v} is the sample variance of finite EAP estimates
(denominator \eqn{n-1}) and \eqn{s} is the mean squared posterior SD for the same
persons. Model separation is \eqn{\sqrt{v/s}}. This describes posterior score
separation conditional on the fitted calibration; it is not repeated-sample
classification accuracy. \code{RealReliability} substitutes the fit-inflated
\code{RealSE} and is a descriptive extension of that formula.
For WLE/JML person scores and facet estimates, true-score variance is estimated
as \eqn{\max(0,v-s)} and reliability as \eqn{\max(0,v-s)/v}.
\code{ReliabilityBasis} identifies which construction is used. Undefined
quantities remain \code{NA}.
}
\examples{
ratings <- load_mfrmr_data("example_operational")
fit <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score")
diag <- diagnose_mfrm(fit)

# Inspect reliability/separation per facet and top unexpected rating rows
diag$reliability
head(diag$unexpected, 5)
summary(diag)
}
