\name{estimate_bias}
\alias{estimate_bias}
\title{Estimate 2-Way and All-Pairs Facet Bias with FACETS 4.5.3 RajuSA-Welch Contrasts}
\description{
Computes fast C++ fixed-measure 2-way local bias tables (FACETS Table 13) and pairwise contrasts (Table 14) for one facet pair or \code{interaction_facets = "all"}, defaulting to the FACETS 4.5.2/4.5.3 \code{RajuSA-Welch} integer-truncated Welch-Satterthwaite degrees of freedom (von Davier, 2025; Raju, 1990), or fits a joint C++ interaction model with omnibus likelihood-ratio test when \code{method = "joint"}.
}
\usage{
estimate_bias(
  fit, facet_a = NULL, facet_b = NULL,
  interaction_facets = NULL, method = c("local", "joint"),
  df_method = c("facets", "classical"), p_adjust = "holm")
}
\arguments{
  \item{fit}{A fitted model object returned by \code{\link{fit_mfrm}}.}
  \item{facet_a}{Optional character string naming the first interaction facet.}
  \item{facet_b}{Optional character string naming the second interaction facet.}
  \item{interaction_facets}{Character vector of two facet names, or \code{"all"} to screen all 2-way non-person facet pairs.}
  \item{method}{\code{"local"} (fixed-calibration post-hoc bias and pairwise contrasts) or \code{"joint"} (simultaneous C++ main-effect + sum-to-zero interaction refit with omnibus LRT).}
  \item{df_method}{Degrees-of-freedom convention for pairwise contrasts: \code{"facets"} (FACETS 4.5.3 \code{RajuSA-Welch} integer truncation with \code{df >= 1} floor) or \code{"classical"} (untruncated Welch-Satterthwaite).}
  \item{p_adjust}{Multiple-testing adjustment method passed to \code{\link[stats]{p.adjust}} (default \code{"holm"}).}
}
\value{
An S3 object of class \code{c("mfrm_bias", "mfrm_bundle")} containing \code{table} (cell-level bias estimates, \code{BiasSE}, \code{t}, \code{df}, \code{pValue}, \code{pAdj}), \code{pairwise} (pairwise context contrasts with \code{TestType}, \code{df}, \code{dfUnrounded}, \code{pAdj}), \code{summary}, and (when \code{method = "joint"}) \code{lrt} and \code{joint_fit}.
}
\details{
Local GPCM calculations retain the fitted observation discrimination products
and the chosen full-predictor or ability-only formula. Their cell curvature SEs
condition on the calibration and person estimates; pairwise SEs combine these
conditional cell errors. They do not propagate calibration uncertainty or
guarantee finite-sample coverage. RajuSA-Welch labels identify the contrast/df
convention, not independent validation of every GPCM extension.

Joint refits append only the declared zero-margin interaction to the source
encoding, retaining weights, anchors/constraints, population covariates,
discrimination specification and numerical controls. An already included
interaction or an empirically shrunk source is rejected. The omnibus LRT is
returned only when the common nested-likelihood checks used by \code{anova()}
are satisfied; otherwise \code{lrt} is \code{NULL} and
\code{comparison_basis} explains the restriction.

\code{InferenceStatus} and \code{InferenceReason} accompany the screening rows.
Unready calibrations retain descriptive estimates but return \code{NA} for
interval SEs and p-values. All-unavailable screen counts are \code{NA}, not zero.
Numerical readiness is not a guarantee of coverage or diagnostic accuracy.
Older saved calibrations require refitting with the current fit_mfrm().
}
\seealso{
\code{\link{fit_mfrm}}, \code{\link{analyze_dff}}, \code{\link{plot.mfrm_fit}}
}
\examples{
# Screen planted Rater x Criterion interaction (R04 x Accuracy) in example_bias
ratings <- load_mfrmr_data("example_bias")
fit <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score")

bias <- estimate_bias(fit, interaction_facets = c("Rater", "Criterion"), df_method = "facets")
head(bias$table[order(-abs(bias$table$t)),
                c("LevelA", "LevelB", "ObsExpAvg", "BiasSize", "BiasSE", "t", "pAdj")], 5)
plot(bias, type = "forest")
}
