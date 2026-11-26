\name{detect_anchor_drift}
\alias{detect_anchor_drift}
\title{Screen Element Parameter Displacement and Anchor Drift Across Waves}
\description{
Compares common facet element location estimates across two or more separately calibrated administration waves against a designated reference wave, flagging elements with large absolute displacement (\code{|Drift| > drift_threshold}) or standardized displacement (\code{|z| > flag_se_ratio}).
}
\usage{
detect_anchor_drift(
  fits, facets = NULL, drift_threshold = 0.5, flag_se_ratio = 2,
  reference = 1L, include_person = FALSE)
}
\arguments{
  \item{fits}{Named or unnamed list of two or more fitted \code{\link{fit_mfrm}} objects representing administration waves.}
  \item{facets}{Optional character vector restricting drift screening to selected facets.}
  \item{drift_threshold}{Absolute logit displacement threshold for flagging (default \code{0.5}).}
  \item{flag_se_ratio}{Standardized \code{|Drift| / SE} threshold for flagging (default \code{2}).}
  \item{reference}{Name or integer index of the reference wave in \code{fits} (default \code{1L}).}
  \item{include_person}{Logical; whether to include common persons (default \code{FALSE}).}
}
\value{
An S3 object of class \code{"mfrm_anchor_drift"} containing \code{table} (per-element \code{EstRef}, \code{EstWave}, \code{Drift}, \code{SE}, \code{z}, \code{pValue}, \code{Flagged}, and \code{Status}) and \code{reference}.
}
\details{
Estimates must already share a meaningful scale, model and origin. This screen
does not estimate an equating transformation. The standard error is the square
root of the sum of the two marginal variances and assumes independent wave
estimates; repeated persons or overlapping ratings can violate that assumption.
When either calibration is not ready for inference, displacement is retained but
standard errors, tests and flags are unavailable. With \code{include_person = TRUE},
common persons are included unless excluded by \code{facets}.
}
\seealso{
\code{\link{anchor_to_baseline}}, \code{\link{fit_mfrm_drift}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
f1 <- fit_mfrm(subset(ratings, Group == "A"), "Person", c("Rater", "Criterion"), "Score")
f2 <- fit_mfrm(subset(ratings, Group == "B"), "Person", c("Rater", "Criterion"), "Score")

drift_chk <- detect_anchor_drift(list(Wave1 = f1, Wave2 = f2), facets = "Rater")
drift_chk$table
}
