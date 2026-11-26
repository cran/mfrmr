\name{compute_person_fit_indices}
\alias{compute_person_fit_indices}
\title{Compute Parametric Person-Fit Statistics (lz and Outfit/Infit Screening)}
\description{
Computes the standardized log-likelihood person-fit index \eqn{l_z} (Drasgow, Levine, & Williams, 1985) alongside person-level \code{InfitMnSq} and \code{OutfitMnSq} from a fitted Many-Facet model or diagnostic bundle to screen aberrant examinee rating patterns.
}
\usage{
compute_person_fit_indices(
  diagnostics, fit = NULL)
}
\arguments{
  \item{diagnostics}{Either a fitted model object from \code{\link{fit_mfrm}} or a diagnostic bundle from \code{\link{diagnose_mfrm}}.}
  \item{fit}{Optional matching \code{\link{fit_mfrm}} object when \code{diagnostics} is a legacy diagnostic bundle.}
}
\value{
A data frame of class \code{mfrm_person_fit_indices} for either engine.
\code{TotalN}, \code{N} and \code{UnavailableN} count retained, usable and
unavailable response terms. \code{LogLik}, \code{lz}, \code{lz_star} and
\code{lz_star_status} describe the statistics and correction eligibility.
\code{ReportIndex}, \code{ReportValue} and \code{ReportFlag} use the corrected
index when available, otherwise uncorrected \eqn{l_z}. \code{ReportCaveat}
records the interpretation limits. Infit and Outfit are separate descriptive
indices. Compatibility aliases \code{ObsCount}, \code{L0} and \code{Flagged}
refer to \code{N}, \code{LogLik} and \code{ReportFlag}; \code{x$table} returns
the same table as a plain data frame. Use \code{summary(x)} for an overview.
}
\details{
For each retained response, category probabilities give the log-probability
mean and variance. If any response term for a person is unavailable, the
person statistic is withheld rather than computed on a selected subset.

For converged, unpenalized JML person estimates with unit weights, the
conditional ability correction projects centered log likelihood off the
person likelihood score. With total log-likelihood variance \eqn{V},
score information \eqn{I}, log-likelihood/score covariance \eqn{C},
centered log likelihood \eqn{L-E(L)} and observed score derivative \eqn{U},
the reported correction is
\deqn{l_z^* = \frac{L-E(L)-(C/I)U}{\sqrt{V-C^2/I}}.}
These moments sum over the person's observed rating design and condition
on the fitted non-person calibration. They do not include calibration
uncertainty. The correction is not supplied for EAP, WLE, penalized or
adjusted person scores, non-unit weights, nonconvergence, or insufficient
information. Their status remains explicit.

The two-sided normal-reference cutoffs \eqn{1.96} and \eqn{2.58} are
screening references, not calibrated individual or familywise error rates.
Short, sparse or dependent rating patterns need design-specific validation.
The plot uses the same reported index and flags; Outfit above 1.50 is
shown separately. Neither flag alone establishes aberrant behavior.
}
\references{
Snijders, T. A. B. (2001). Asymptotic null distribution of person fit
statistics with estimated person parameter. Psychometrika, 66, 331--342.
\doi{10.1007/BF02294437}.

Sinharay, S. (2016). Asymptotically correct standardization of person-fit
statistics beyond dichotomous items. Psychometrika, 81, 992--1013.
\doi{10.1007/s11336-015-9465-x}.
}
\seealso{
\code{\link{diagnose_mfrm}}, \code{\link{plot.mfrm_fit}}
}
\examples{
ratings <- load_mfrmr_data("example_core")
fit <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score")

pfit <- compute_person_fit_indices(fit)
summary(pfit)
plot(fit, type = "person_fit")
}
