\name{tidymfrm-s3-methods}
\alias{tidymfrm-s3-methods}
\alias{anova.tidymfrm}
\alias{coef.tidymfrm}
\alias{confint.tidymfrm}
\alias{predict.tidymfrm}
\alias{print.tidymfrm}
\alias{print.tidymfrm_dashboard}
\alias{summary.tidymfrm}
\alias{vcov.tidymfrm}
\title{Standard S3 Methods for tidymfrm Workflow Objects}
\description{
Standard S3 methods (\code{anova}, \code{predict}, \code{coef}, \code{vcov}, \code{confint},
\code{summary}, and \code{print}) for \code{tidymfrm} workflow objects.
}
\usage{
\method{anova}{tidymfrm}(object, ...)
}
\arguments{
  \item{object}{A calibrated \code{tidymfrm} workflow.}
  \item{...}{For \code{anova()}, additional calibrated \code{tidymfrm} workflows.}
}
\details{
\code{anova()} compares separately fitted models using the underlying
\code{mfrmr} likelihood and nesting checks. Compare fits with the same retained
ratings, weights, likelihood method, population and integration settings.
Inspect \code{comparison_basis} when \code{lrt} is unavailable. Information
criteria are interpretable as differences only for comparable likelihoods.
\code{compare_methods} in \code{step_mfrm_measure()} addresses sensitivity to
the estimator within a specified model.
}
\value{
Standard S3 return values delegated from the underlying fitted \code{mfrmr}
model or workflow summary.
}
\examples{
utils::data("mfrmr_example_core", package = "mfrmr")
wf <- mfrm_workflow(mfrmr_example_core, "Person", c("Rater", "Criterion"), "Score") |>
  step_mfrm_measure()
coef(wf)
confint(wf)
head(predict(wf, type = "response"))
}
