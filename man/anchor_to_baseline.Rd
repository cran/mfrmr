\name{anchor_to_baseline}
\alias{anchor_to_baseline}
\title{Link and Calibrate New Ratings Against Baseline Facet Anchors}
\description{
Extracts element anchors from a baseline calibration (\code{baseline_fit}), calibrates a new rating cohort (\code{new_data}) on the baseline logit scale, and optionally compares benchmark elements to estimates linked through declared stable reference anchors.
}
\usage{
anchor_to_baseline(
  new_data, baseline_fit, person, facets, score, anchor_facets = NULL,
  include_person = FALSE, weight = NULL, model = NULL, method = NULL,
  anchor_policy = "warn", mml_integration = NULL, centering_scope = NULL,
  ..., drift_reference = NULL)
}
\arguments{
  \item{new_data}{A \code{data.frame} of observed ratings from the new administration or wave.}
  \item{baseline_fit}{A fitted baseline model object from \code{\link{fit_mfrm}}.}
  \item{person}{Character string naming the person identifier column in \code{new_data}.}
  \item{facets}{Character vector naming non-person facet columns in \code{new_data}.}
  \item{score}{Character string naming the numeric score column in \code{new_data}.}
  \item{anchor_facets}{Optional subset of \code{facets} to anchor to \code{baseline_fit} (defaults to all common non-person facets).}
  \item{include_person}{Logical; whether to anchor common persons as well (default \code{FALSE}).}
  \item{weight}{Optional observation weight column name in \code{new_data}.}
  \item{model}{Optional model family override; defaults to \code{baseline_fit}'s model.}
  \item{method}{Optional estimation method override; defaults to \code{baseline_fit}'s method.}
  \item{anchor_policy}{Policy for anchors absent from \code{new_data} (\code{"warn"}, \code{"error"}, or \code{"silent"}).}
  \item{mml_integration}{Optional MML quadrature mode override.}
  \item{centering_scope}{Optional centering scope (\code{"remainder"} or \code{"facet"}).}
  \item{drift_reference}{Optional data frame of stable baseline anchor IDs with columns \code{Facet} and \code{Level}, for the unified engine. Reference elements set the origin for non-reference elements in the same facet.}
  \item{...}{Additional arguments forwarded to \code{\link{fit_mfrm}}.}
}
\value{
An S3 object of class \code{"mfrm_anchored_fit"} containing \code{fit} (new cohort calibrated on the baseline scale), \code{diagnostics}, \code{baseline_anchors}, and \code{drift} (element-by-element change on the linked scale).
For unified fits without an independent reference, \code{Drift} and \code{Flag}
are unavailable. In supported additive RSM/PCM fits, the comparison is conditional
on the supplied baseline reference and uses the new-fit standard error.
It is distinct from the local conditional displacement reported for a fixed
element or group in a single fitted calibration. Supports \code{summary()}.
}
\seealso{
\code{\link{make_anchor_table}}, \code{\link{detect_anchor_drift}}, \code{\link{fit_mfrm}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
wave1 <- subset(ratings, Group == "A")
wave2 <- subset(ratings, Group == "B")

base_fit <- fit_mfrm(wave1, person = "Person", facets = c("Rater", "Criterion"), score = "Score")
linked <- anchor_to_baseline(
  wave2,
  baseline_fit = base_fit,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  anchor_facets = "Criterion"
)
linked$drift
}
