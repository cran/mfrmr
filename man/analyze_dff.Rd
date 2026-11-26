\name{analyze_dff}
\alias{analyze_dff}
\title{Differential Facet Functioning (Parallel Split-Element LRT, RajuSA-Welch, MH, and Logistic)}
\description{
Screens Differential Facet Functioning (DFF) of raters, tasks, or criteria across reference and focal groups via warm-started parallel C++ split-element \code{MML} LRT (with person-clustered sandwich SEs and latent-regression group impact control), \code{CML} LRT, FACETS 4.5.3 \code{RajuSA-Welch} contrasts, stratified Mantel-Haenszel, or cumulative-logit LRT.
}
\usage{
analyze_dff(
  fit, facet, group, data = NULL, focal = NULL,
  method = c("mml_lrt", "cml_lrt", "raju_welch", "mh", "logistic"),
  df_method = c("facets", "classical"), min_obs = 10, p_adjust = "holm",
  threads = 1L)
}
\arguments{
  \item{fit}{A fitted model object returned by \code{\link{fit_mfrm}}.}
  \item{facet}{Character string naming the model facet to screen for differential functioning (e.g., \code{"Rater"} or \code{"Criterion"}).}
  \item{group}{Character string naming the grouping variable in \code{data} (e.g., examinee cohort or demographic group).}
  \item{data}{Optional \code{data.frame} containing \code{group} alongside the rating columns; defaults to the retained calibration data.}
  \item{focal}{Optional character string specifying the focal group level in \code{group}. Positive \code{Contrast} indicates higher severity/difficulty against the focal group.}
  \item{method}{DFF test family: \code{"mml_lrt"} (default C++ split-element MML likelihood-ratio test), \code{"cml_lrt"} (C++ split-element CML LRT), \code{"raju_welch"} (fixed-calibration pairwise contrast), \code{"mh"} (stratified mean-score Mantel-Haenszel), \code{"logistic"} (cumulative-logit uniform/non-uniform LRT).}
  \item{df_method}{Degrees-of-freedom convention when \code{method = "raju_welch"} (\code{"facets"} or \code{"classical"}).}
  \item{min_obs}{Minimum required ratings per group for an element to be tested (default \code{10}).}
  \item{p_adjust}{Multiple-testing adjustment method passed to \code{\link[stats]{p.adjust}} (default \code{"holm"}).}
  \item{threads}{Integer number of C++ worker threads for parallel split-element refits (default \code{1L}).}
}
\value{
An S3 object of class \code{"mfrm_dff"} containing \code{table} (\code{Facet}, \code{Level}, \code{ReferenceGroup}, \code{FocalGroup}, \code{ObsRef}, \code{ObsFocal}, \code{Contrast}, \code{SE}, \code{Statistic}, \code{df}, \code{pValue}, \code{pAdj}, \code{Flagged}), \code{method}, \code{reference}, and \code{focal}. Supports \code{summary()} and \code{plot()}.
}
\details{
When \env{_R_CHECK_LIMIT_CORES_} requests a core limit, native worker counts are capped at two. The default remains one worker.

Split-likelihood tests support unified, unpenalized RSM/PCM calibrations without
anchors, custom facet/step constraints, interactions or restricted category
support. Unsupported specifications are rejected rather than dropped.
Weights, facet signs and integration/optimizer controls are retained.
The selected test method determines the refitting estimator. With
\code{method = "mml_lrt"}, both models are refitted by MML even if the
supplied calibration was estimated by JML.
MML refits retain source population covariates and add person-level group
impact when linearly independent; a fixed-normal source is explicitly extended
to an estimated population. Both the baseline and split models use this same
population design. A group varying within person cannot be a person-level
impact covariate. CML tests condition on person totals and do not use a
population distribution.

Both likelihood refits must pass numerical readiness and likelihood checks.
Otherwise \code{SE}, \code{Statistic}, \code{pValue}, \code{pAdj} and
\code{Flagged} are \code{NA}, with \code{InferenceStatus}/\code{InferenceReason}.
\code{ElementsTested} in summaries counts available tests;
\code{ElementsUnavailable} counts unavailable rows. An all-unavailable
\code{Flagged} count is \code{NA}. These checks do not establish finite-sample
test size or power.

RajuSA-Welch analysis holds the source calibration fixed, including GPCM
slopes, and appends a zero-valued group context without recalibration. Other
fixed-calibration screens also withhold inferential output if the source is
not ready. Supplied \code{data} must contain exactly the retained calibration
ratings and weights; reordered rows are aligned before adding group labels.
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{estimate_bias}}, \code{\link{plot.mfrm_fit}}
}
\examples{
# Saved JML calibration of example_bias, with planted Group x Criterion DFF.
saved <- readRDS(system.file("examples", "bias-jml.rds", package = "mfrmr"))
ratings <- load_mfrmr_data(saved$dataset)
fit <- saved$fit
# Fixed-calibration screen: no split-model refitting.
dff <- analyze_dff(fit, facet = "Criterion", group = "Group", data = ratings,
                   method = "raju_welch")
summary(dff)
plot(dff)
}
