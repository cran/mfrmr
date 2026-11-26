\name{fair_average_table}
\alias{fair_average_table}
\title{Compute Reference-Standardized Fair(M) and Zero-Reference Fair(Z) Expected Scores}
\description{
Translates latent logit parameter estimates for persons and non-person facets back onto the original rubric score scale by computing model-expected average ratings against a cohort-mean reference context (\code{FairM}) and a zero-logit reference context (\code{FairZ}), alongside raw observed averages (\code{ObsAvg}).
}
\usage{
fair_average_table(fit, diagnostics = NULL, facets = NULL, umean = 0,
  uscale = 1, fair_se = FALSE, ci_level = 0.95, reference_plan = NULL)
}
\arguments{
  \item{fit}{A fitted model object returned by \code{\link{fit_mfrm}}.}
  \item{diagnostics}{Optional precomputed \code{\link{diagnose_mfrm}} object from the same calibration.}
  \item{facets}{Optional character vector filtering which facets (including \code{"Person"}) to return.}
  \item{umean}{Additive origin shift applied to reported logit measures \code{umean + uscale * Estimate} (default \code{0}).}
  \item{uscale}{Multiplicative scale factor applied to reported logit measures and \code{SE} (default \code{1}).}
  \item{fair_se}{Logical; request uncertainty. Supply a common \code{reference_plan}.}
  \item{ci_level}{Confidence level in \code{(0, 1)} when intervals are requested (default \code{0.95}).}
  \item{reference_plan}{Optional common rating plan \code{data.frame} containing exactly all fitted base facet columns plus optional \code{Weight}. Levels must be known, cells unique and weights finite and strictly positive. Weights are normalized to sum to one.}
}
\value{
An S3 object of class \code{"mfrm_fair_average"} containing \code{stacked} (combined \code{data.frame} with \code{Facet}, \code{Level}, \code{ObsCount}, \code{ObsAvg}, \code{ExpAvg}, \code{FairM}, \code{FairZ}, \code{Measure}, and \code{SE}) and \code{raw_by_facet} (named list of per-facet \code{data.frame}s).
When a plan is supplied, also returns \code{reference_plan} with
\code{NormalizedWeight}, and \code{reference_scores} with one person per row.
}
\details{
Printing shows the common-plan \code{reference_scores}, including conditional
intervals when requested, before the observed-design Fair(M) tables.
Observed-design Fair(M) uses finite facet-level means, an observation-weighted
geometric reference discrimination for GMFRM, and each person's weighted
step-group mix. A geometric-mean-one identification constraint is not itself a
reference-panel definition. Use an explicit \code{reference_plan} to compare
people or models against the same specified rating plan.

\code{confint(x, level = 0.95)} returns these conditional person-score
intervals for a fair-average table with a reference plan. These intervals hold
calibration parameters fixed; they do not propagate estimated-anchor or
calibration uncertainty. Legacy refit-bootstrap and joint-delta arguments are
not supported by this method.

For the unified engine, \code{reference_scores$Estimate} is the weighted mean
model-expected rating under the same supplied plan for every person, conditional
on that person's primary ability estimate. It is on the original rubric scale
(\code{Units = "original_score"}). The plan preserves the fitted locations,
steps, interactions and slopes; its scores can therefore differ from FairM/FairZ.
An EAP plug-in score is not the posterior mean of the transformed score.
With \code{fair_se = TRUE}, the unified engine transforms the endpoints of
a normal interval for the primary person measure under the supplied plan. For
EAP these are conditional normal posterior approximations. Calibration and plan
weights are held fixed. The point estimate remains a plug-in expected score.
Unavailable bounds remain \code{NA} with an uncertainty status. The \code{SE}
column in the facet table stays on the latent measure scale.
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{diagnose_mfrm}}
}
\examples{
# Compare raw observed averages (ObsAvg) against rater-adjusted Fair(M) on incomplete panels
ratings <- load_mfrmr_data("example_operational")
fit <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score")

fair <- fair_average_table(fit, facets = c("Person", "Rater"))
head(fair$raw_by_facet$Person, 5)
fair$raw_by_facet$Rater

# Put every person on the same weighted panel (original score scale).
plan <- data.frame(Rater = c("R01", "R04"),
                   Criterion = c("Content", "Language"), Weight = c(1, 3))
head(fair_average_table(fit, reference_plan = plan)$reference_scores)
}
