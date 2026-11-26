\name{mfrm_g_study}
\alias{mfrm_g_study}
\title{Estimate Univariate and Multivariate Generalizability Study Variance Components}
\description{
Estimates variance and covariance components for balanced or unbalanced crossed and nested rating designs using Henderson Method 1 (\code{method = "urgenova"}), ANOVA, or \code{MINQUE(0)}.
}
\usage{
mfrm_g_study(
  data, scores = NULL, person = "Person", rater = "Rater", task = "Task",
  method = c("anova", "minque0", "urgenova"), missing = c("error", "omit"),
  facets = NULL, nesting = NULL, score = NULL, weights = NULL,
  d_study_grid = NULL, truncate_negative = TRUE)
}
\arguments{
  \item{data}{A \code{data.frame} of observed ratings.}
  \item{scores}{Character vector naming one or more score columns.}
  \item{person}{Character string naming the person identifier column.}
  \item{rater}{Character string naming the rater column when \code{facets = NULL}.}
  \item{task}{Character string naming the task column when \code{facets = NULL}.}
  \item{method}{Variance-component estimator (\code{"anova"}, \code{"minque0"}, or \code{"urgenova"}).}
  \item{missing}{Missing-value policy (\code{"error"} or \code{"omit"}).}
  \item{facets}{Character vector of non-person facet column names.}
  \item{nesting}{Optional named character vector \code{c(Child = "Parent")} declaring nested facets.}
  \item{score}{Alias for \code{scores} in univariate designs.}
  \item{weights}{Optional nonnegative composite weights matching \code{scores}. Weights are used as supplied: \code{c(1, 1)} forms a sum, while \code{c(.5, .5)} forms an average. Named weights are matched to score names; unnamed weights follow score order. The default is an equal-weight average. Requires explicit \code{method = "urgenova"}; weights never change the selected estimator. For ANOVA or MINQUE(0), supply composite weights to \code{mfrm_d_study()} after fitting the separate scores.}
  \item{d_study_grid}{Optional candidate D-Study design grid \code{data.frame}, with one column per facet and positive integer counts only. Nested-facet counts are per parent level.}
  \item{truncate_negative}{Logical; whether to clamp negative variance estimates at zero.}
}
\value{
An \code{mfrm_gtheory} object for every estimator, with \code{method_code},
\code{design}, raw score covariance matrices in \code{components}, composite
\code{variance_components}, \code{g_study_summary}, and a \code{d_study} table.
Method-specific numerical ingredients and unavailable summaries are recorded
explicitly. Incomplete ANOVA/MINQUE designs require a future \code{design_grid};
observed pool sizes are not silently interpreted as planned replication.
Use \code{summary(g)$component_diagnostics} to inspect covariance definiteness.
For Henderson, \code{vc_vcov} describes the supplied composite. For ANOVA/MINQUE,
\code{vcov(g, score = ..., assumption = "normal")} selects score or composite
sampling uncertainty. Balanced ANOVA supports \code{confint(g, method = "mls",
assumption = "normal")}; this interval construction is not applied to unbalanced
Henderson estimates. \code{tidy(g)} reports covariance entries for ANOVA/MINQUE;
Henderson's composite table additionally includes interval status.

}
\details{
Henderson Method 1 solves the expected mean-product equations for the declared
crossed or nested random-effects design. A rank-deficient coefficient matrix is
rejected: the data cannot distinguish the requested variance components.
For \code{urgenova}, unidentified components and too few complete observations
signal an error of class \code{"mfrm_gstudy_unavailable"}.
Negative raw estimates are retained in \code{RawEstimate}; \code{truncate_negative}
controls their use in projected point estimates.
For \code{urgenova}, zero total variance or negative universe/projected error
variance gives an unavailable G/Phi (\code{NA}), with the reason in
\code{GStatus}/\code{PhiStatus}; see \code{\link{mfrm_d_study}}.

The \code{urgenova} sampling covariance uses the Gaussian quadratic-form identity
\eqn{\operatorname{Cov}(y'Ay,y'By)=2\operatorname{tr}(AVBV)} with the actual
assignment pattern, including imbalance and replication. The working covariance
\eqn{V} substitutes nonnegative component estimates (negative raw estimates are
zeroed for this calculation). This is a plug-in covariance, not an exact
finite-sample distribution or a model-free uncertainty estimate.
\code{confint()} uses Satterthwaite chi-square approximations for positive
components and logit delta intervals for G/Phi; boundary coefficients may have
unavailable intervals (\code{NA}). The approximation is particularly limited
near zero components and in sparse designs. As a conservative reporting guard,
positive components with Satterthwaite degrees of freedom below 1 have both
endpoints set to \code{NA}. This cutoff is not a coverage guarantee for other
components. The named \code{status} attribute of the interval matrix records
\code{unstable_satterthwaite_df}, the approximation used, or numerical
unavailability. \code{tidy(gt)} and workflow reports carry these statuses and
95\% component intervals. MLS intervals for the balanced ANOVA representation
use different mean-product information and are not substituted for this
unbalanced Henderson approximation.
}
\examples{
ratings <- load_mfrmr_data("example_core")
gt <- mfrm_g_study(ratings, scores = "Score", person = "Person",
                   facets = c("Rater", "Criterion"), method = "urgenova")
mfrm_d_study(gt)
}
