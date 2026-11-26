\name{compute_information}
\alias{compute_information}
\title{Compute Category Probability and Fisher Test Information Curves}
\description{
Evaluates Fisher test information \eqn{I(\theta)}, conditional standard error \eqn{\mathrm{SE}(\theta) = 1 / \sqrt{I(\theta)}}, expected total score, and reference category probability curves across a grid of latent ability values for \code{"RSM"}, \code{"PCM"}, and \code{"GPCM"} calibrations.
}
\usage{
compute_information(
  fit, theta_range = c(-6, 6), theta_points = 201L, rating_plan = NULL,
  exposure_basis = c("observed_total", "per_person_average"))
}
\arguments{
  \item{fit}{A fitted model object returned by \code{\link{fit_mfrm}}.}
  \item{theta_range}{Numeric vector of length 2 giving the lower and upper logit bounds of the \code{theta} grid (default \code{c(-6, 6)}).}
  \item{theta_points}{Integer number of evaluation points along \code{theta_range} (default \code{201L}).}
  \item{rating_plan}{Optional \code{data.frame} with every non-Person facet and positive integer \code{Count} for one person's planned ratings. Supply calibrated levels only and combine repeated combinations in \code{Count}.}
  \item{exposure_basis}{With observed ratings, use \code{"observed_total"} for all retained exposure or \code{"per_person_average"} to divide by the number of retained persons. With the argument omitted, the reference curve counts each observed context once.}
}
\value{
An S3 object of class \code{"mfrm_information"} containing \code{test_information} (\code{Theta}, \code{Information}, \code{SE}, \code{ExpectedScore}) and \code{category_curves} (\code{ScaleLevel}, \code{Theta}, \code{Category}, \code{Probability}). Supports \code{\link{plot.mfrm_fit}}.
}
\details{
Information is additive for conditionally independent ratings with fixed
calibration. Conditional SE does not propagate calibration uncertainty.
Fitted calibrations report the exposure basis and counts and flag new
combinations of calibrated levels in
\code{exposures$ObservedContext}. An unsupported or numerically unready
calibration cannot be used for a precision curve.
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{classification_table}}, \code{\link{plot.mfrm_fit}}
}
\examples{
# Saved JML calibration of example_core; see fit_mfrm() for estimation.
saved <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
fit <- saved$fit

info <- compute_information(fit, theta_range = c(-4, 4), theta_points = 81L)
head(info$test_information)
plot(info, type = "information")
}
