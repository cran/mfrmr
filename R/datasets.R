#' Simulated MFRM datasets based on Eckes and Jin (2021)
#'
#' Synthetic many-facet rating datasets in long format.
#' All datasets include one row per observed rating.
#'
#' Available data objects:
#' - `mfrmr_example_core`
#' - `mfrmr_example_bias`
#' - `ej2021_study1`
#' - `ej2021_study2`
#' - `ej2021_combined`
#' - `ej2021_study1_itercal`
#' - `ej2021_study2_itercal`
#' - `ej2021_combined_itercal`
#'
#' @format A data.frame with 5 columns:
#' \describe{
#'   \item{Study}{Study label (`"Study1"` or `"Study2"`).}
#'   \item{Person}{Person/respondent identifier.}
#'   \item{Rater}{Rater identifier.}
#'   \item{Criterion}{Criterion facet label.}
#'   \item{Score}{Observed category score.}
#' }
#' @source Simulated for this package with design settings informed by Eckes and Jin (2021).
#' @details
#' Naming convention:
#' - `study1` / `study2`: separate simulation studies
#' - `combined`: row-bind of study1 and study2
#' - `_itercal`: iterative-calibration variant
#'
#' Use [load_mfrmr_data()] for programmatic selection by key.
#'
#' @section Data dimensions:
#' \tabular{lrrrr}{
#'   \strong{Dataset} \tab \strong{Rows} \tab \strong{Persons} \tab \strong{Raters} \tab \strong{Criteria} \cr
#'   study1 \tab 1842 \tab 307 \tab 18 \tab 3 \cr
#'   study2 \tab 3287 \tab 206 \tab 12 \tab 9 \cr
#'   combined \tab 5129 \tab 307 \tab 18 \tab 12 \cr
#'   study1_itercal \tab 1842 \tab 307 \tab 18 \tab 3 \cr
#'   study2_itercal \tab 3341 \tab 206 \tab 12 \tab 9 \cr
#'   combined_itercal \tab 5183 \tab 307 \tab 18 \tab 12
#' }
#' Score range: 1--4 (four-category rating scale).
#'
#' @section Simulation design:
#' Person ability is drawn from N(0, 1).  Rater severity effects span
#' approximately -0.5 to +0.5 logits.  Criterion difficulty effects span
#' approximately -0.3 to +0.3 logits.  Scores are generated from the
#' resulting linear predictor plus Gaussian noise, then discretized into
#' four categories.  The `_itercal` variants use a second iteration of
#' calibrated rater severity parameters.
#'
#' @section Interpreting output:
#' Each dataset is already in long format and can be passed directly to
#' [fit_mfrm()] after confirming column-role mapping.
#'
#' @section Typical workflow:
#' 1. Inspect available datasets with [list_mfrmr_data()].
#' 2. Load one dataset using [load_mfrmr_data()].
#' 3. Fit and diagnose with [fit_mfrm()] and [diagnose_mfrm()].
#'
#' @examples
#' data("ej2021_study1", package = "mfrmr")
#' head(ej2021_study1)
#' table(ej2021_study1$Study)
#' @name ej2021_data
#' @aliases ej2021_study1 ej2021_study2 ej2021_combined
#' @aliases ej2021_study1_itercal ej2021_study2_itercal ej2021_combined_itercal
NULL

#' Purpose-built example datasets for package help pages
#'
#' Compact synthetic many-facet datasets designed for documentation examples.
#' Both datasets are large enough to avoid tiny-sample toy behavior while
#' remaining fast in `R CMD check` examples.
#'
#' Available data objects:
#' - `mfrmr_example_core`
#' - `mfrmr_example_bias`
#'
#' @format A data.frame with 6 columns:
#' \describe{
#'   \item{Study}{Example dataset label (`"ExampleCore"` or `"ExampleBias"`).}
#'   \item{Person}{Person/respondent identifier.}
#'   \item{Rater}{Rater identifier.}
#'   \item{Criterion}{Criterion facet label.}
#'   \item{Score}{Observed category score on a four-category scale (`1`--`4`).}
#'   \item{Group}{Balanced grouping variable used in DFF/DIF examples (`"A"` / `"B"`).}
#' }
#' @source Synthetic documentation data generated from rating-scale Rasch facet
#'   designs with fixed seeds in `data-raw/make-example-data.R`.
#' @details
#' `mfrmr_example_core` is generated from a single latent trait plus rater and
#' criterion main effects, making it suitable for general fitting, plotting, and
#' reporting examples.
#'
#' `mfrmr_example_bias` starts from the same basic design but adds:
#' - a known `Group x Criterion` effect (`Group B` is advantaged on `Language`)
#' - a known `Rater x Criterion` interaction (`R04 x Accuracy`)
#'
#' This lets differential-functioning and bias-analysis help pages demonstrate non-null findings.
#'
#' @section Data dimensions:
#' \tabular{lrrrrr}{
#'   \strong{Dataset} \tab \strong{Rows} \tab \strong{Persons} \tab \strong{Raters} \tab \strong{Criteria} \tab \strong{Groups} \cr
#'   example_core \tab 768 \tab 48 \tab 4 \tab 4 \tab 2 \cr
#'   example_bias \tab 384 \tab 48 \tab 4 \tab 4 \tab 2
#' }
#'
#' @section Suggested usage:
#' - Use `mfrmr_example_core` for fitting, diagnostics, design-weighted precision curves,
#'   and generic plots/reports.
#' - Use `mfrmr_example_bias` for [analyze_dff()], [analyze_dif()], [dif_interaction_table()],
#'   [plot_dif_heatmap()], and [estimate_bias()].
#'
#' Both objects can be loaded either with [load_mfrmr_data()] or directly via
#' `data("mfrmr_example_core", package = "mfrmr")` /
#' `data("mfrmr_example_bias", package = "mfrmr")`.
#'
#' @examples
#' data("mfrmr_example_core", package = "mfrmr")
#' table(mfrmr_example_core$Score)
#' table(mfrmr_example_core$Group)
#' @name mfrmr_example_data
#' @aliases mfrmr_example_core mfrmr_example_bias
NULL

#' List packaged simulation datasets
#'
#' @return Character vector of dataset keys accepted by [load_mfrmr_data()].
#' @details
#' Use this helper when you want to select packaged data programmatically
#' (e.g., inside scripts, loops, or shiny/streamlit wrappers).
#'
#' Typical pattern:
#' 1. call `list_mfrmr_data()` to see available keys.
#' 2. pass one key to [load_mfrmr_data()].
#'
#' @section Interpreting output:
#' Returned values are canonical dataset keys accepted by [load_mfrmr_data()].
#'
#' @section Typical workflow:
#' 1. Capture keys in a script (`keys <- list_mfrmr_data()`).
#' 2. Select one key by index or name.
#' 3. Load data via [load_mfrmr_data()] and continue analysis.
#'
#' @seealso [load_mfrmr_data()], [ej2021_data]
#' @examples
#' keys <- list_mfrmr_data()
#' keys
#' d <- load_mfrmr_data(keys[1])
#' head(d)
#' @export
list_mfrmr_data <- function() {
  c(
    "example_core",
    "example_bias",
    "study1",
    "study2",
    "combined",
    "study1_itercal",
    "study2_itercal",
    "combined_itercal"
  )
}

#' Load a packaged simulation dataset
#'
#' @param name Dataset key. One of values from [list_mfrmr_data()].
#'
#' @return A data.frame in long format.
#' @details
#' This helper is useful in scripts/functions where you want to choose a dataset
#' by string key instead of calling `data()` manually.
#'
#' All returned datasets include the core long-format columns
#' `Study`, `Person`, `Rater`, `Criterion`, and `Score`.
#' Some datasets, such as the packaged documentation examples, also include
#' auxiliary variables like `Group` for DIF/bias demonstrations.
#'
#' @section Interpreting output:
#' The return value is a plain long-format `data.frame`, ready for direct use
#' in [fit_mfrm()] without additional reshaping.
#'
#' @section Typical workflow:
#' 1. list valid names with [list_mfrmr_data()].
#' 2. load one dataset key with `load_mfrmr_data(name)`.
#' 3. fit a model with [fit_mfrm()] and inspect with `summary()` / `plot()`.
#'
#' @seealso [list_mfrmr_data()], [ej2021_data]
#' @examples
#' data("mfrmr_example_core", package = "mfrmr")
#' head(mfrmr_example_core)
#'
#' d <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   data = d,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 25
#' )
#' summary(fit)
#' @export
load_mfrmr_data <- function(name = c(
                            "example_core",
                            "example_bias",
                            "study1",
                            "study2",
                            "combined",
                            "study1_itercal",
                            "study2_itercal",
                            "combined_itercal"
                          )) {
  key <- match.arg(tolower(name), choices = list_mfrmr_data())

  obj_name <- switch(
    key,
    example_core = "mfrmr_example_core",
    example_bias = "mfrmr_example_bias",
    study1 = "ej2021_study1",
    study2 = "ej2021_study2",
    combined = "ej2021_combined",
    study1_itercal = "ej2021_study1_itercal",
    study2_itercal = "ej2021_study2_itercal",
    combined_itercal = "ej2021_combined_itercal"
  )

  utils::data(list = obj_name, package = "mfrmr", envir = environment())
  get(obj_name, envir = environment(), inherits = FALSE)
}
