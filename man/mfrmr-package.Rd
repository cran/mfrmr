\name{mfrmr-package}
\alias{mfrmr}
\alias{mfrmr-package}
\docType{package}
\title{Fast Estimation, Diagnostics, and Generalizability for Many-Facet Models}
\description{
Estimates Many-Facet Rasch Measurement models (\code{RSM}, \code{PCM}, and \code{GPCM}) including step constraints, declarative category support, explanatory facet covariates, and joint interactions via high-speed C++ Joint Maximum Likelihood, Jeffreys-penalized JML with Block-Schur standard errors, Gauss-Hermite and Laplace-adaptive Marginal Maximum Likelihood (with \code{EAP}, Warm's \code{WLE}, latent regression, and sandwich covariances), Conditional Maximum Likelihood, and Markov rater drift models. Also provides multi-threaded batch estimation, graph connectivity checks, \code{Fair(M)}, \code{Real SE}, \code{RajuSA-Welch} pairwise bias and parallel split-element DFF, color-vision-safe \pkg{ggplot2} graphics, and Generalizability (\code{G-Study}) and Decision (\code{D-Study}) planning.
}
\examples{
ratings <- load_mfrmr_data("example_core")
fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score")
summary(fit)$overview
}
