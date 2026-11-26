\name{apply_empirical_bayes_shrinkage}
\alias{apply_empirical_bayes_shrinkage}
\title{Apply Empirical Bayes Shrinkage to Facet and Person Estimates}
\description{
Applies post-estimation precision-weighted Empirical Bayes shrinkage toward each facet's mean location using estimated between-element variance \eqn{\hat{\tau}^2 = \max(0, s_{\hat{\beta}}^2 - \overline{\mathrm{SE}^2})} (or a user-supplied \code{facet_prior_sd}), stabilizing noisy estimates for low-exposure raters or tasks.
}
\usage{
apply_empirical_bayes_shrinkage(
  fit, facet_prior_sd = NULL, shrink_person = FALSE)
}
\arguments{
  \item{fit}{A fitted model object returned by \code{\link{fit_mfrm}}.}
  \item{facet_prior_sd}{Optional positive numeric prior standard deviation (in logits) overriding empirical \eqn{\hat{\tau}}.}
  \item{shrink_person}{Logical; if \code{TRUE}, also computes shrunk person estimates in \code{fit$person} (default \code{FALSE}).}
}
\value{
An updated fit object retaining original \code{Estimate} and adding \code{ShrunkEstimate} and \code{ShrunkSE} columns in \code{fit$facets} (and \code{fit$person} when \code{shrink_person = TRUE}).
}
\seealso{
\code{\link{fit_mfrm}}
}
\examples{
ratings <- load_mfrmr_data("example_operational")
fit <- fit_mfrm(
  ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score"
)

shrunk <- apply_empirical_bayes_shrinkage(fit)
cols <- c("Facet", "Level", "ObsCount", "Estimate", "ShrunkEstimate", "ShrunkSE")
shrunk$facets[, cols]
}
