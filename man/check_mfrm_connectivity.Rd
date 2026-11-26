\name{check_mfrm_connectivity}
\alias{check_mfrm_connectivity}
\title{Audit Multi-Facet Rating Graph Connectivity and Category Support in C++}
\description{
Performs a fast pre-estimation C++ union-find connectivity and category-support audit on a rating table before fitting a Many-Facet Rasch model. Identifies disconnected rating subsets (islands lacking shared persons, raters, tasks, anchors, or explanatory covariates) and checks whether interior rubric categories are observed.
}
\usage{
check_mfrm_connectivity(
  data, person, facets, score, model = c("RSM", "PCM", "GPCM"),
  step_facet = NULL, rating_min = NULL, rating_max = NULL, anchors = NULL,
  group_anchors = NULL, dummy_facets = NULL, category_support = NULL,
  facet_formulas = NULL, facet_data = NULL, rater_facet = NULL,
  network = FALSE, rater_panels = NULL)
}
\arguments{
  \item{data}{A \code{data.frame} of observed ratings (one row per rating).}
  \item{person}{Character string naming the person identifier column.}
  \item{facets}{Character vector naming non-person facet columns (e.g., \code{c("Rater", "Criterion")}).}
  \item{score}{Character string naming the numeric rubric score column.}
  \item{model}{Measurement model family (\code{"RSM"}, \code{"PCM"}, or \code{"GPCM"}).}
  \item{step_facet}{Facet owning category thresholds when \code{model} is \code{"PCM"} or \code{"GPCM"}.}
  \item{rating_min}{Optional integer minimum rubric score.}
  \item{rating_max}{Optional integer maximum rubric score.}
  \item{anchors}{Optional individual element anchor \code{data.frame} (\code{Facet}, \code{Level}, \code{Anchor}) linking subsets to a common logit origin.}
  \item{group_anchors}{Optional group-mean anchor \code{data.frame} (\code{Facet}, \code{Level}, \code{Group}, \code{GroupValue}).}
  \item{dummy_facets}{Optional character vector of facets fixed at 0 logits.}
  \item{category_support}{Optional \code{data.frame} with \code{Level}, \code{Min}, and \code{Max} declaring inclusive per-level bounds on the original rubric scale.}
  \item{facet_formulas}{Optional named list of explanatory facet covariate formulas.}
  \item{facet_data}{Optional named list of element covariate \code{data.frame}s.}
  \item{network}{Logical. Compute the Fiedler value of the unnormalized observed Person-by-rater Laplacian. Large graphs use the optional RSpectra package.}
  \item{rater_panels}{Optional named character vector assigning observed rater IDs to two or more panels. Compute resistance between equal-weight panel mean contrasts.}
  \item{rater_facet}{Facet for the observed Person-by-rater comparison graph. If \code{NULL}, infer a rater-like name when possible.}
}
\value{
A list containing the multi-facet graph status and counts, per-person score
and subset summaries, facet subsets, and the interior category support flag.
When a rater facet is available, \code{rater_linkage} gives the observed
Person-by-rater component table and the persons whose individual removal
separates the remaining observed panels. This graph uses shared persons,
independently of anchors and shared levels of other facets. Each observed
Person-by-rater pair contributes one edge regardless of the number of criteria
or repeated scores. The Fiedler value is zero for a disconnected graph.
A panel contrast that requires current to cross disconnected components has
infinite effective resistance. These statistics describe assignment structure.
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{describe_mfrm_data}}
}
\examples{
# Check incomplete two-rater operational assignment (48 persons, 6 raters, 3 criteria)
ratings <- load_mfrmr_data("example_operational")
conn <- check_mfrm_connectivity(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score"
)
conn$connected
conn$n_subsets
conn$rater_linkage
}
