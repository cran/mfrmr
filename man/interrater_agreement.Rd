\name{interrater_agreement}
\alias{interrater_agreement}
\title{Compute Exact, Weighted Kappa, Gwet AC2, Krippendorff Alpha, and ICC Agreement}
\description{
Consolidates exact and adjacent agreement, pairwise weighted kappa, Gwet's ordinal \eqn{AC_2}, Krippendorff's \eqn{\alpha}, and two-way random intraclass correlation (\code{ICC(2,1)} and \code{ICC(2,k)}) from raw rating tables or fitted \code{mfrm_fit} objects.
}
\usage{
interrater_agreement(
  data, person = "Person", rater = "Rater", score = "Score",
  context_facets = NULL, unit = NULL, categories = NULL,
  metric = c("all", "exact", "kappa", "gwet_ac2", "krippendorff_alpha", "icc"),
  weights = c("quadratic", "linear", "unweighted"), model_expected = FALSE
)
}
\arguments{
  \item{data}{A \code{data.frame} or fitted \code{mfrm_fit} object of observed ratings.}
  \item{person}{Character string naming the person identifier column.}
  \item{rater}{Character string naming the rater identifier column.}
  \item{score}{Character string naming the numeric rating column.}
  \item{context_facets}{Optional character vector of additional facets identifying rating units.}
  \item{unit}{Optional character vector specifying the full unit identifier columns (e.g., \code{c("Person", "Task")}), overriding \code{person} and \code{context_facets}.}
  \item{categories}{Optional numeric vector of declared rubric score categories (e.g., \code{0:3}).}
  \item{metric}{Agreement family: \code{"all"}, \code{"exact"}, \code{"kappa"}, \code{"gwet_ac2"}, \code{"krippendorff_alpha"}, or \code{"icc"}.}
  \item{model_expected}{Logical. For a unified-engine fit, also compute observed and expected rates and model-based exact kappa on the same usable pairs.}
  \item{weights}{Weighting scheme: \code{"quadratic"}, \code{"linear"}, or \code{"unweighted"}.}
}
\value{
A structured \pkg{mfrmr} S3 object or \code{data.frame} containing parameter estimates, standard errors, and diagnostic summaries.
}
\details{
Model expectation is the product of the two fitted category probability vectors,
conditional on each fitted person's primary score and independent ratings.
Complete unit-rater pairs with valid probabilities contribute equally. The model
table records pair counts and availability. The ordinary weighted kappa uses
marginal category frequencies.

Each unit-rater combination must have at most one observed rating. Duplicate
combinations are rejected to prevent row-order-dependent overwriting. For repeated
occasions, include the occasion column in \code{unit} or \code{context_facets}
so that the comparison refers to the intended rating occasion.
}
\examples{
ratings <- load_mfrmr_data("example_core")
interrater_agreement(ratings, person = "Person", rater = "Rater",
                     score = "Score", context_facets = "Criterion")
}
