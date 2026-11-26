\name{make_anchor_table}
\alias{make_anchor_table}
\title{Extract Direct Calibration Anchor Tables from a Fitted Model}
\description{
Extracts converged facet element location estimates (and optionally person estimates) from a baseline \code{\link{fit_mfrm}} object in the \code{data.frame(Facet, Level, Anchor)} format accepted by the \code{anchors} argument of subsequent \code{\link{fit_mfrm}} calibrations.
}
\usage{
make_anchor_table(
  fit, facets = NULL, include_person = FALSE, digits = NULL,
  readiness_policy = c("error", "review"))
}
\arguments{
  \item{fit}{A fitted baseline model object returned by \code{\link{fit_mfrm}}.}
  \item{facets}{Optional character vector of non-person facet names to include; defaults to all non-person facets.}
  \item{include_person}{Logical; if \code{TRUE}, also exports finite person ability estimates as \code{Facet = "Person"} rows.}
  \item{digits}{Optional integer rounding precision applied to \code{Anchor} values.}
  \item{readiness_policy}{\code{"error"} (rejects non-converged or unready legacy calibrations) or \code{"review"}.}
}
\value{
A \code{data.frame} (or \code{tibble}) with columns \code{Facet}, \code{Level}, and \code{Anchor}.
}
\seealso{
\code{\link{anchor_to_baseline}}, \code{\link{detect_anchor_drift}}, \code{\link{fit_mfrm}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
fit_base <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score")

anchors <- make_anchor_table(fit_base, facets = "Criterion")
anchors
}
