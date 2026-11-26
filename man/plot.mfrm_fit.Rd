\name{plot.mfrm_fit}
\alias{plot.mfrm_bias}
\alias{plot.mfrm_diagnostics}
\alias{plot.mfrm_dff}
\alias{plot.mfrm_drift}
\alias{plot.mfrm_fit}
\alias{plot.mfrm_gtheory}
\alias{plot.mfrm_information}
\alias{plot.mfrm_multivariate_d_study}
\title{Color-Vision-Safe ggplot2 Visualizations for Many-Facet Models and Diagnostics}
\description{
Unified \pkg{ggplot2} \code{plot()} methods for fitted MFRM models, diagnostics, G/D-Studies, pairwise \code{RajuSA-Welch} bias tables, DFF screens, test information curves, and rater drift trajectories. Every multi-group aesthetic pairs color (\code{"okabe_ito"}, \code{"viridis"}, or strict \code{R == G == B} \code{"monochrome"}) with distinct point \code{shape} and \code{linetype} encodings. Canvas titles and footnotes are omitted by default and stored in \code{attr(p, "mfrm_notes")}.
}
\usage{
\method{plot}{mfrm_bias}(
  x, type = c("forest", "pairwise"), top_n = 20L,
  palette = c("okabe_ito", "viridis", "monochrome"), base_size = 11, ...)

\method{plot}{mfrm_diagnostics}(
  x, facet = NULL, palette = c("okabe_ito", "viridis", "monochrome"),
  base_size = 11, ..., type = c("pathway", "residual_pca"))

\method{plot}{mfrm_dff}(
  x, top_n = 25L, palette = c("okabe_ito", "viridis", "monochrome"),
  base_size = 11, ...)

\method{plot}{mfrm_drift}(
  x, type = c("trajectory", "contrasts"),
  palette = c("okabe_ito", "viridis", "monochrome"), base_size = 11, ...)

\method{plot}{mfrm_fit}(
  x,
  type = c("wright", "pathway", "steps", "category_curves", "persons",
    "raters", "person_fit", "residual_qq", "slopes", "qc"),
  palette = c("okabe_ito", "viridis", "monochrome"), facet = NULL,
  base_size = 11, ..., fit_statistic = c("infit", "outfit", "both"))

\method{plot}{mfrm_gtheory}(
  x, type = NULL,
  palette = c("okabe_ito", "viridis", "monochrome"), base_size = 11, ...)

\method{plot}{mfrm_information}(
  x, type = c("information", "ccc"),
  palette = c("okabe_ito", "viridis", "monochrome"), base_size = 11, ...)

\method{plot}{mfrm_multivariate_d_study}(
  x, palette = c("okabe_ito", "viridis", "monochrome"), base_size = 11, ...)



}
\arguments{
  \item{fit_statistic}{Horizontal mean-square statistic for a pathway plot. Use \code{"both"} for separate Infit and Outfit panels on the same scale.}
  \item{x}{An \pkg{mfrmr} model, diagnostic, bias, DFF, information, drift, G/D-Study object.}
  \item{type}{Character string selecting the visualization panel: for diagnostic objects, \code{"pathway"} or \code{"residual_pca"} (first residual contrast); for fitted models, \code{"wright"}, \code{"pathway"}, \code{"steps"}, \code{"category_curves"}, \code{"persons"}, \code{"raters"}, \code{"person_fit"}, \code{"residual_qq"}, \code{"slopes"} (GPCM), or \code{"qc"}; for G-Studies, \code{"dstudy"} or \code{"variance"} (\code{NULL} chooses D-Study for Henderson and component variance for ANOVA/MINQUE); for bias, \code{"forest"} or \code{"pairwise"}; for information, \code{"information"} or \code{"ccc"}; for drift, \code{"trajectory"} or \code{"contrasts"}.}
  \item{top_n}{Maximum number of ranked bias cells or DFF elements to display.}
  \item{palette}{Color-vision-safe palette: \code{"okabe_ito"} (default), \code{"viridis"}, or strict \code{"monochrome"} grayscale.}
  \item{base_size}{Base font size in points (default \code{11}).}
  \item{...}{Additional plotting options.}
  \item{facet}{Optional facet name filter for element pathway and rater diagnostic plots.}
}
\value{
A \code{ggplot} object with diagnostic and scale notes attached in \code{attr(p, "mfrm_notes")}.
}
\details{
The \code{"person_fit"} panel plots Outfit against the reported person-fit
index from \code{compute_person_fit_indices()}: conditional \eqn{l_z^*}
when available, otherwise \eqn{l_z}. Its two-sided flags agree with
\code{summary()} of that result. These are screening references, not
calibrated significance tests.

For unified-engine fits, slope bars use the same log-scale 95 percent model
intervals as \code{confint(x, component = "slopes")}. Person and facet diagnostic
bars use descriptive Real-SE intervals. When the calibration fails inference
readiness, these interval bars are withheld while finite descriptive points are
retained; \code{attr(p, "mfrm_notes")} records the reasons. Neither visible bars
nor a ready status establish nominal finite-sample coverage.
}
\examples{
# Saved JML calibration of example_core; see fit_mfrm() for estimation.
saved <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
fit <- saved$fit
p1 <- plot(fit, type = "wright", palette = "okabe_ito")
p2 <- plot(fit, type = "category_curves", palette = "viridis")
p3 <- plot(fit, type = "residual_qq", palette = "monochrome")
attr(p1, "mfrm_notes")
}
