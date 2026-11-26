\name{compare_mfrm}
\alias{compare_mfrm}
\title{Compare Many-Facet Models via Information Criteria and Affine-Subspace LRT}
\description{
Compares two or more fitted Many-Facet models by log-likelihood, \code{AIC}, \code{BIC}, and (when \code{nested = TRUE} or called via \code{anova()}) formal likelihood-ratio tests after algebraically verifying affine parameter-subspace nesting across facet locations, step constraints, and GPCM log-slopes.
}
\usage{
compare_mfrm(
  ..., labels = NULL, warn_constraints = TRUE, nested = FALSE)
}
\arguments{
  \item{...}{Two or more fitted model objects from \code{\link{fit_mfrm}}.}
  \item{labels}{Optional character vector of model labels.}
  \item{warn_constraints}{Logical; whether to warn when compared models differ in estimation method, sample size, or non-nested parameter subspaces.}
  \item{nested}{Logical; if \code{TRUE}, verifies affine subspace containment via SVD and computes sequential chi-square likelihood-ratio tests for nested \code{MML} or \code{CML} models.}
}
\value{
An S3 object of class \code{"mfrm_comparison"} containing \code{table} (\code{Model}, \code{Family}, \code{Method}, \code{N}, \code{FreeParams}, \code{LogLik}, \code{AIC}, \code{BIC}, \code{DeltaAIC}, \code{DeltaBIC}), \code{lrt} (sequential LRT table when \code{nested = TRUE} and models are verified nested), and \code{comparison_basis}.
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{mfrmr-s3-methods}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
fit_rsm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score", model = "RSM")
fit_pcm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
                    model = "PCM", step_facet = "Criterion")

cmp <- compare_mfrm(fit_rsm, fit_pcm, labels = c("RSM", "PCM"), nested = TRUE)
cmp$table
cmp$lrt
}
