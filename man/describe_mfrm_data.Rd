\name{describe_mfrm_data}
\alias{describe_mfrm_data}
\title{Summarize Rating Design Structure, Sparsity, Overlap, and Missing Codes}
\description{
Audits raw rating tables prior to calibration: summarizes category frequencies, element exposure counts, rater-pair shared-unit overlap, connectivity, duplicate cells, declared missing-value codes, and optional structural omissions relative to a planned assignment roster (\code{expected_design}).
}
\usage{
describe_mfrm_data(
  data, person = NULL, facets = NULL, score = NULL, weight = NULL,
  rating_min = NULL, rating_max = NULL, keep_original = FALSE,
  missing_codes = NULL, include_person_facet = FALSE,
  include_agreement = TRUE, rater_facet = NULL, context_facets = NULL,
  agreement_top_n = NULL, expected_design = NULL, min_linking_persons = 2L,
  category_policy = NULL, repeated_ratings = c("review", "independent"),
  comparisons = NULL, spec = NULL)
}
\arguments{
  \item{data}{A \code{data.frame} of observed ratings.}
  \item{person}{Character string naming the person identifier column.}
  \item{facets}{Character vector naming non-person facet columns.}
  \item{score}{Character string naming the numeric rating column.}
  \item{weight}{Optional character string naming an observation weight column.}
  \item{rating_min}{Optional integer minimum rubric score.}
  \item{rating_max}{Optional integer maximum rubric score.}
  \item{keep_original}{Logical; whether to retain original uncleaned columns in internal audit tables.}
  \item{missing_codes}{Optional vector of values in \code{score} to recode as missing (\code{NA}) prior to summarizing.}
  \item{include_person_facet}{Logical; whether to include individual person exposure rows in element tables.}
  \item{include_agreement}{Logical; whether to compute pairwise rater overlap and exact agreement summaries (default \code{TRUE}).}
  \item{rater_facet}{Optional character string identifying the rater facet for overlap summaries.}
  \item{context_facets}{Optional character vector of non-rater context facets defining rating units.}
  \item{agreement_top_n}{Optional integer limiting displayed pairwise overlap rows.}
  \item{expected_design}{Optional score-free roster \code{data.frame} of planned facet cells used to distinguish planned omissions from unassigned cells.}
  \item{min_linking_persons}{Minimum shared persons treated as adequate pairwise rater linkage (default \code{2L}).}
  \item{category_policy}{Optional category-gap handling policy (\code{"preserve"} or \code{"collapse"}).}
  \item{repeated_ratings}{Policy for duplicate person-by-facet cells (\code{"review"} or \code{"independent"}).}
  \item{comparisons}{Optional subgroup comparison specification.}
  \item{spec}{Optional named specification list shared with \code{\link{fit_mfrm}}.}
}
\value{
An S3 object of class \code{"mfrm_data_description"} supporting \code{print()} and \code{summary()}, containing design counts, category distributions, element exposure summaries, rater overlap tables, and structural missingness checks.
}
\seealso{
\code{\link{check_mfrm_connectivity}}, \code{\link{interrater_agreement}}, \code{\link{fit_mfrm}}
}
\examples{
ratings <- load_mfrmr_data("example_operational")
desc <- describe_mfrm_data(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score"
)
summary(desc)$overview
}
