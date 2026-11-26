\name{mfrm_d_study}
\alias{mfrm_d_study}
\title{Project Decision-Study Generalizability (G), Dependability (Phi), and SEM}
\description{
Projects universe-score variance, relative error variance, absolute error variance, relative and absolute \code{SEM}, generalizability coefficient (\eqn{G}), and dependability index (\eqn{\Phi}) across candidate facet sample sizes (\code{design_grid}) from an estimated \code{\link{mfrm_g_study}} or \code{\link{fit_mfrm}} object.
}
\usage{
mfrm_d_study(
  x, design_grid = NULL, ...)
}
\arguments{
  \item{x}{An estimated G-Study object returned by \code{\link{mfrm_g_study}} (or a fitted \code{\link{fit_mfrm}} object).}
  \item{design_grid}{Optional \code{data.frame} specifying positive integer sample sizes, with one column per facet (e.g., \code{expand.grid(Rater = 1:4, Criterion = 2:5)}). For nested facets, counts are per parent level. For Henderson, \code{NULL} returns the default grid around observed harmonic facet counts. For balanced complete ANOVA/MINQUE designs it projects the observed counts; incomplete or unequal designs require an explicit future grid.}
  \item{...}{Additional arguments forwarded to methods (such as \code{weights} for \code{"anova"}/\code{"minque0"} multivariate G-Study objects).}
}
\value{
A \code{data.frame} (for \code{method = "urgenova"} G-Studies) or \code{"mfrm_multivariate_d_study"} object reporting candidate facet counts, \code{TotalRatingsPerPerson}, \code{UniverseVar}, \code{RelErrorVar}, \code{AbsErrorVar}, \code{RelSEM}, \code{AbsSEM}, \code{G}, and \code{Phi}.
}
\details{
For \code{urgenova}, \code{GStatus}/\code{PhiStatus} identify availability.
Zero total variance makes the coefficient undefined (\code{NA}), not zero;
its \code{MeetsG80}/\code{MeetsPhi80} flag is also \code{NA}.
Negative universe or projected error variance likewise makes the associated
coefficient unavailable. A zero universe variance with positive error variance
is a defined coefficient of zero. Negative error variance yields \code{NA} SEM.
Default and explicit plans use the same calculations, including for saved fits.
}
\seealso{
\code{\link{mfrm_g_study}}, \code{\link{plot.mfrm_fit}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
gt <- mfrm_g_study(
  ratings,
  scores = "Score",
  person = "Person",
  facets = c("Rater", "Criterion"),
  method = "urgenova"
)

# Compare 1-to-3 raters crossed with 2 or 4 criteria
mfrm_d_study(gt, design_grid = expand.grid(Rater = 1:3, Criterion = c(2, 4)))
}
