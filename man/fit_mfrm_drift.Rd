\name{fit_mfrm_drift}
\alias{fit_mfrm_drift}
\title{Estimate Dynamic Multi-Wave and Markov Rater Severity Drift Models}
\description{
Estimates non-Bayesian dynamic rater-severity drift across ordered scoring waves, sessions, or days without MCMC sampling (generalized from Uto, 2023). Integrates person abilities by normal quadrature, computes analytical Louis (1982) observed information for inner parameters, and estimates rater-specific random-walk drift standard deviations by Laplace-marginalized Empirical Bayes with Kass-Steffey (1989) hyperparameter variance propagation (\code{TotalSE} vs \code{ConditionalSE}).
}
\usage{
fit_mfrm_drift(
  data, person, rater, time, score, facets = character(),
  score_levels = NULL, model = c("PCM", "RSM"),
  drift_mode = c("markov_hetero", "markov_homo", "independent", "static"),
  eb_method = c("em_laplace", "full_laplace"), mu_sigma = -2,
  delta_sigma = 1, beta1_sd = 1, step_sd = 1, person_sd = 1,
  logistic_scale = 1, quad_points = 81L, maxit = 60L, tol = 1e-06,
  missing = c("fail", "omit"), drift_threshold = 0.3, flag_se_ratio = 2)
}
\arguments{
  \item{data}{A \code{data.frame} of observed ratings across ordered waves or occasions.}
  \item{person}{Character string naming the person identifier column.}
  \item{rater}{Character string naming the rater identifier column whose severity may drift over \code{time}.}
  \item{time}{Character string naming the ordered scoring wave, session, or day column.}
  \item{score}{Character string naming the ordered rubric score column.}
  \item{facets}{Optional character vector of additional time-invariant non-person facets (e.g., \code{"Criterion"}).}
  \item{score_levels}{Optional sorted numeric/integer vector declaring the full ordered rubric categories.}
  \item{model}{Category threshold structure: \code{"PCM"} (rater-specific step thresholds) or \code{"RSM"} (shared step thresholds).}
  \item{drift_mode}{Temporal severity specification: \code{"markov_hetero"} (rater-specific Markov drift SDs), \code{"markov_homo"} (common Markov drift SD), \code{"independent"} (unsmoothed per-wave severities), or \code{"static"} (time-invariant baseline).}
  \item{eb_method}{Hyperparameter optimization method: \code{"em_laplace"} (default) or \code{"full_laplace"}.}
  \item{mu_sigma}{Prior mean for log-drift standard deviations (default \code{-2}).}
  \item{delta_sigma}{Prior SD for log-drift standard deviations (default \code{1}).}
  \item{beta1_sd}{Prior SD for initial Wave-1 rater severities (default \code{1}).}
  \item{step_sd}{Weak regularizing prior SD on step parameters (default \code{1}).}
  \item{person_sd}{Standard deviation of the normal person ability distribution (default \code{1}).}
  \item{logistic_scale}{Multiplicative scale constant (default \code{1}; use \code{1.7} for normal-ogive scaling).}
  \item{quad_points}{Number of normal quadrature grid points for person integration (default \code{81L}).}
  \item{maxit}{Maximum outer optimization iterations (default \code{60L}).}
  \item{tol}{Relative convergence tolerance (default \code{1e-06}).}
  \item{missing}{Missing-value handling: \code{"fail"} or \code{"omit"}.}
  \item{drift_threshold}{Absolute logit change threshold used to flag substantive wave-to-wave drift (default \code{0.3}).}
  \item{flag_se_ratio}{Standardized \code{|Drift| / TotalSE} ratio threshold for flagging drift contrasts (default \code{2}).}
}
\value{
An S3 object of class \code{"mfrm_drift"} supporting \code{summary()}, \code{plot()}, \code{confint()}, \code{coef()}, \code{vcov()}, \code{predict()}, and \code{as.data.frame(x, component = "severity")} (returning rater-by-wave severity estimates, \code{ConditionalSE}, \code{TotalSE}, and drift flags).
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{detect_anchor_drift}}, \code{\link{plot.mfrm_fit}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
ratings$Wave <- ifelse(as.integer(factor(ratings$Person)) <= 24L, "Day1", "Day2")

drift_fit <- fit_mfrm_drift(
  ratings,
  person = "Person",
  rater = "Rater",
  time = "Wave",
  score = "Score",
  facets = "Criterion",
  model = "RSM"
)
as.data.frame(drift_fit, component = "severity")
plot(drift_fit, type = "trajectory")
}
