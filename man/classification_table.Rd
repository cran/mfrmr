\name{classification_table}
\alias{classification_table}
\title{Compute Cut-Score Classification Accuracy and Consistency Indices}
\description{
Computes conditional and population-weighted classification accuracy, consistency, overclassification, and underclassification probabilities under a specified multi-facet rating plan using both Lee's (2010) exact summed-score convolution method (\code{"lee_sum"}) and Rudner's (2005) asymptotic normal ability-estimate approximation (\code{"rudner_normal"}).
}
\usage{
classification_table(
  fit, rating_plan = NULL, cut_scores = 0, score_cuts = NULL,
  theta = seq(-3, 3, by = 0.5), theta_weights = NULL,
  operational_cuts = cut_scores, facets = NULL,
  readiness_policy = c("error", "review"))
}
\arguments{
  \item{fit}{A converged \code{"RSM"} or \code{"PCM"} calibration returned by \code{\link{fit_mfrm}}.}
  \item{rating_plan}{A \code{data.frame} listing the planned facet cells (one row per rating in one examinee's assessment plan). Defaults to the unique observed facet combinations on C++ fits.}
  \item{cut_scores}{Sorted numeric vector of true latent ability cutoffs in logits.}
  \item{score_cuts}{Nondecreasing integer vector of summed-score cutoffs (on \code{0..max_sum}), one per \code{cut_scores} entry, for Lee's summed-score rule. Computed proportionally if \code{NULL} on C++ fits.}
  \item{theta}{Numeric vector of latent ability evaluation points.}
  \item{theta_weights}{Optional nonnegative numeric weights matching \code{theta} to compute marginal population accuracy/consistency and the joint classification matrix.}
  \item{operational_cuts}{Sorted numeric vector of operational logit ability cutoffs used by Rudner's decision rule (defaults to \code{cut_scores}).}
  \item{facets}{Optional character vector mapping \code{rating_plan} column names to model facets.}
  \item{readiness_policy}{\code{"error"} (rejects calibrations with boundary exclusions or non-convergence) or \code{"review"} (retains conditional indices tagged \code{SourceStatus = "review_required"}).}
}
\value{
A list containing \code{table} (conditional accuracy, consistency, over/underclassification, expected sum, information, and SEM at each \code{theta}), \code{marginal} (population-weighted summary when \code{theta_weights} is supplied), \code{classification_matrix} (joint true-by-observed category probabilities), \code{cuts}, \code{probabilities}, \code{scores}, and \code{notes}.
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{compute_information}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
fit <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score")

# Evaluate a 2-rater x 4-criterion (8-rating, max category sum = 24) pass/fail plan at 0 logits
plan <- expand.grid(
  Rater = c("R01", "R02"),
  Criterion = c("Content", "Organization", "Language", "Accuracy"),
  stringsAsFactors = FALSE
)
th_grid <- seq(-3, 3, by = 0.5)
cls <- classification_table(
  fit,
  rating_plan = plan,
  cut_scores = 0,
  score_cuts = 12,
  theta = th_grid,
  theta_weights = stats::dnorm(th_grid)
)
cls$marginal
}
