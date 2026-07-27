#' Legacy synthetic MFRM datasets inspired by Eckes and Jin (2021)
#'
#' Synthetic many-facet rating datasets in long format.
#' All datasets include one row per observed rating.
#'
#' Available data objects:
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
#' @source Simulated for this package with design settings informed by Eckes
#' and Jin (2021). The Eckes & Jin (2021) Method section reports the
#' following design parameters that motivated the synthetic versions
#' shipped here: Study 1 had 307 examinees (149 males, 158 females), 18
#' raters (4 males, 14 females), and 3 criteria (global impression, task
#' fulfillment, linguistic realization) on a 4-category rating scale (TDN
#' levels rescored 1-4); Study 2 had 206 examinees (66 males, 140 females),
#' 12 raters (1 male, 11 females), and 9 criteria on the same 4-category
#' scale. The packaged datasets reproduce these
#' (examinees, raters, criteria, categories) shapes but use simulated
#' responses, so they are not the real TestDaF data.
#'
#' @references
#' Eckes, T., & Jin, K.-Y. (2021). Measuring rater centrality effects in
#' writing assessment: A Bayesian facets modeling approach.
#' \emph{Psychological Test and Assessment Modeling, 63}(1), 65--94.
#' @details
#' Naming convention:
#' - `study1` / `study2`: separate simulation studies
#' - `combined`: row-bind of study1 and study2
#' - `_itercal`: legacy synthetic sensitivity variant. These objects can differ
#'   in observed rows as well as scores and should not be interpreted as a
#'   controlled one-parameter recalibration.
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
#' For the combined rows, `Persons` and `Raters` count unique raw labels.
#' Treating the two Study labels as distinct namespaces would instead give 513
#' person labels and 30 rater labels, but would leave two unlinked components.
#'
#' @section Provenance and limits:
#' These are legacy synthetic datasets whose stored responses reproduce the
#' dimensions described above. The exact response-generation code and random
#' seed are not available, so the objects must not be used as parameter-recovery
#' evidence or as evidence for a particular generating distribution. The
#' separately stored `_itercal` objects can differ in observed rows as well as
#' scores; they are legacy variants, not empirical calibration standards or a
#' controlled one-parameter recalibration.
#'
#' @section Combined-data caution:
#' The `combined` objects reuse `P001`--`P206` and `R01`--`R12` across the two
#' study labels. A combined analysis is meaningful only when those labels encode
#' an intended cross-study identity and an explicit anchor or linking design
#' establishes a common scale. Prefixing `Person` and `Rater` by `Study` removes
#' accidental label collisions, but it creates disconnected study components
#' and does not by itself establish a common scale. Analyze the studies
#' separately unless the linking design has been specified and reviewed. The
#' combined datasets are not beginner workflow examples or direct-fit examples.
#'
#' @section Interpreting output:
#' The study-specific datasets are in long format and can be passed to
#' [fit_mfrm()] after confirming column-role mapping. The `combined` objects
#' require a reviewed identity and anchor/linking design before a joint fit;
#' row-binding or prefixing identifiers alone does not create a common scale.
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

#' Synthetic many-facet rating examples
#'
#' Compact synthetic many-facet datasets sized to exercise the documented
#' workflow. Their dimensions are not evidence of sample-size adequacy for an
#' applied study.
#'
#' Available data objects:
#' - `mfrmr_example_operational`
#' - `mfrmr_example_operational_design` (documented separately below)
#' - `mfrmr_example_core`
#' - `mfrmr_example_bias`
#'
#' @format A data.frame with 6 columns:
#' \describe{
#'   \item{Study}{Example dataset label (`"OperationalExample"`,
#'     `"ExampleCore"`, or `"ExampleBias"`).}
#'   \item{Person}{Person/respondent identifier.}
#'   \item{Rater}{Rater identifier.}
#'   \item{Criterion}{Criterion facet label.}
#'   \item{Score}{Observed category score on a four-category scale (`1`--`4`).}
#'   \item{Group}{Balanced grouping label (`"A"` / `"B"`). It is neutral in
#'     the operational and core examples; the bias example has the planted
#'     group structure described below.}
#' }
#' @source Synthetic documentation data generated from rating-scale Rasch facet
#'   designs with fixed seeds. Generator scripts are maintained under
#'   `data-raw/` for this release in the public source repository.
#'   These are synthetic examples, not empirical records.
#' @details
#' `mfrmr_example_operational` is the primary applied teaching example. It has
#' a connected but incomplete two-rater assignment, moderately unequal rater
#' workloads, and six planned criterion-level omissions. Scores are sampled
#' directly from stated RSM category probabilities. The six unobserved ratings
#' are absent rows in the long data rather than `NA` scores. Each group contains
#' 24 persons, and three omissions in each group leave 141 observed rows for
#' Group A and 141 for Group B. The balanced `Group` variable has no effect in
#' the generating model; random observed differences may still occur.
#'
#' `mfrmr_example_core` is an idealized complete crossing generated from a
#' single latent trait plus rater and criterion main effects. It is useful as a
#' fast deterministic example, but it is not representative of routine
#' incomplete operational assignment.
#'
#' `mfrmr_example_bias` instead uses a balanced partial two-rater assignment.
#' Group A and B latent means are -0.1 and 0.1 logits, respectively, with a
#' common SD of 0.9. It also contains:
#' - a planted `Group x Criterion` effect (`Group B` is advantaged by 1.2
#'   logits on `Language`)
#' - a planted `Rater x Criterion` interaction (`R04 x Accuracy` lowers the
#'   linear predictor by 1.2 logits)
#'
#' This lets differential-functioning and bias-analysis help pages demonstrate non-null findings.
#'
#' @section Data dimensions:
#' \tabular{lrrrrr}{
#'   \strong{Dataset} \tab \strong{Rows} \tab \strong{Persons} \tab \strong{Raters} \tab \strong{Criteria} \tab \strong{Groups} \cr
#'   example_operational \tab 282 \tab 48 \tab 6 \tab 3 \tab 2 \cr
#'   example_core \tab 768 \tab 48 \tab 4 \tab 4 \tab 2 \cr
#'   example_bias \tab 384 \tab 48 \tab 4 \tab 4 \tab 2
#' }
#'
#' @section Suggested usage:
#' - Use `mfrmr_example_operational` for the beginner data-to-report workflow
#'   and for inspecting incomplete but connected assignment.
#' - Use `mfrmr_example_core` for fast, idealized checks and examples that
#'   specifically require complete crossing.
#' - Use `mfrmr_example_bias` for [analyze_dff()], [analyze_dif()], [dif_interaction_table()],
#'   [plot_dif_heatmap()], and [estimate_bias()].
#'
#' All three objects can be loaded either with [load_mfrmr_data()] or directly
#' with `data()`, for example
#' `data("mfrmr_example_operational", package = "mfrmr")`.
#'
#' @examples
#' data("mfrmr_example_operational", package = "mfrmr")
#' table(mfrmr_example_operational$Score)
#' table(mfrmr_example_operational$Group)
#' @name mfrmr_example_data
#' @aliases mfrmr_example_operational mfrmr_example_core mfrmr_example_bias
NULL

#' Planned assignment roster for the operational example
#'
#' A score-free roster declaring all Person x Rater x Criterion cells planned
#' for `mfrmr_example_operational`. Pass it to the `expected_design` argument of
#' [describe_mfrm_data()] to distinguish the six expected-but-unobserved cells
#' from combinations that were never assigned.
#'
#' @format A data.frame with 288 rows and 5 columns:
#' \describe{
#'   \item{Study}{Example dataset label (`"OperationalExample"`).}
#'   \item{Person}{Person/respondent identifier.}
#'   \item{Rater}{Planned rater identifier.}
#'   \item{Criterion}{Planned criterion label.}
#'   \item{Group}{Balanced grouping label (`"A"` / `"B"`).}
#' }
#' @source Synthetic assignment roster generated with the operational example
#'   by `data-raw/make-operational-example.R`. It contains no empirical records
#'   and no fabricated scores.
#' @details
#' The roster contains 288 planned cells. The observed
#' `mfrmr_example_operational` table contains 282 rows, so an explicit design
#' comparison identifies six planned omissions. Extra roster columns such as
#' `Study` and `Group` are ignored unless they are named as model facets.
#'
#' @examples
#' data("mfrmr_example_operational", package = "mfrmr")
#' data("mfrmr_example_operational_design", package = "mfrmr")
#' review <- describe_mfrm_data(
#'   mfrmr_example_operational,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   expected_design = mfrmr_example_operational_design
#' )
#' summary(review)$structural_missingness
#' @name mfrmr_example_operational_design
NULL

#' List packaged simulation datasets
#'
#' @param details If `FALSE` (default), return the backward-compatible
#'   character vector of dataset keys. If `TRUE`, return a catalog describing
#'   the intended teaching role and design of every dataset.
#'
#' @return A character vector of dataset keys accepted by [load_mfrmr_data()],
#'   or, when `details = TRUE`, a data.frame with `Key`, `Rows`, `Persons`,
#'   `Raters`, `Criteria`, `CountBasis`, `PrimaryUse`, `Design`, and
#'   `Empirical`.
#' @details
#' Use this helper when you want to select packaged data programmatically
#' (e.g., inside scripts, loops, or interactive-application wrappers).
#'
#' Typical pattern:
#' 1. call `list_mfrmr_data()` to see available keys.
#' 2. pass one key to [load_mfrmr_data()].
#'
#' @section Interpreting output:
#' With `details = FALSE`, returned values are canonical dataset keys accepted
#' by [load_mfrmr_data()]. With `details = TRUE`, use `PrimaryUse` and
#' `Design` to distinguish the applied teaching example from idealized,
#' planted-effect diagnostic, and larger sparse-design datasets. `CountBasis`
#' states whether person/rater counts use raw labels or Study-prefixed labels.
#' Every bundled dataset is synthetic rather than empirical.
#'
#' @section Typical workflow:
#' 1. Capture keys in a script (`keys <- list_mfrmr_data()`).
#' 2. Select one key by index or name.
#' 3. Load data via [load_mfrmr_data()] and continue analysis. Treat a
#'    `combined` key as a design-review object, not as a direct-fit example.
#'
#' @seealso [load_mfrmr_data()], [ej2021_data]
#' @examples
#' keys <- list_mfrmr_data()
#' keys
#' list_mfrmr_data(details = TRUE)[, c(
#'   "Key", "PrimaryUse", "Design", "CountBasis"
#' )]
#' d <- load_mfrmr_data("example_operational")
#' head(d)
#' @export
list_mfrmr_data <- function(details = FALSE) {
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("`details` must be TRUE or FALSE.", call. = FALSE)
  }
  keys <- c(
    "example_core",
    "example_bias",
    "example_operational",
    "study1",
    "study2",
    "combined",
    "study1_itercal",
    "study2_itercal",
    "combined_itercal"
  )
  if (!isTRUE(details)) return(keys)

  catalog <- data.frame(
    Key = keys,
    Rows = c(768L, 384L, 282L, 1842L, 3287L, 5129L, 1842L, 3341L, 5183L),
    Persons = c(48L, 48L, 48L, 307L, 206L, 307L, 307L, 206L, 307L),
    Raters = c(4L, 4L, 6L, 18L, 12L, 18L, 18L, 12L, 18L),
    Criteria = c(4L, 4L, 3L, 3L, 9L, 12L, 3L, 9L, 12L),
    CountBasis = c(
      rep("unique labels", 5L),
      "raw labels; 513 persons and 30 raters when Study-prefixed",
      rep("unique labels", 2L),
      "raw labels; 513 persons and 30 raters when Study-prefixed"
    ),
    PrimaryUse = c(
      "Idealized fast examples",
      "DFF and bias demonstrations with planted effects",
      "Beginner applied workflow",
      "Unequal-workload sparse-design review",
      "Larger sparse-design review",
      "Identity/linking design review; not direct fit",
      "Legacy synthetic variant review",
      "Legacy synthetic variant review",
      "Identity/linking sensitivity review; not direct fit"
    ),
    Design = c(
      "Complete crossing; no planned omissions",
      "Balanced two-rater assignment; planted non-null effects",
      "Connected two-rater assignment; six planned omissions",
      "Two raters per person; highly unequal rater workloads",
      "Two raters per person; incomplete criterion coverage",
      "Overlapping IDs; requires explicit anchors/linking for a common scale",
      "Legacy Study 1 variant; rows and scores can differ",
      "Legacy Study 2 variant; rows and scores can differ",
      "Overlapping IDs; requires explicit anchors/linking for a common scale"
    ),
    Empirical = FALSE,
    stringsAsFactors = FALSE
  )
  catalog
}

#' Load a packaged simulation dataset
#'
#' @param name Dataset key. One of the values from [list_mfrmr_data()]. If
#'   omitted, the backward-compatible default is `"example_core"`; new code
#'   should pass a key explicitly.
#'
#' @return A data.frame in long format.
#' @details
#' `load_mfrmr_data("<key>")` is the canonical loader for the packaged
#' datasets and the entry point used across the package help and
#' vignettes. The equivalent base-R alternative
#' `data("<object-name>", package = "mfrmr")` remains available for users
#' who prefer the full `data()` spelling; both paths return identical
#' long-format data frames.
#'
#' All returned datasets include the core long-format columns
#' `Study`, `Person`, `Rater`, `Criterion`, and `Score`.
#' Some datasets, such as the packaged documentation examples, also include
#' auxiliary variables like `Group` for DIF/bias demonstrations.
#'
#' @section Interpreting output:
#' The return value is a plain long-format `data.frame`. The example and
#' study-specific keys are ready for [fit_mfrm()] after checking role and score
#' mappings. The `combined` keys are design-review objects: overlapping IDs or
#' simple Study-based prefixes do not establish a common measurement scale, so
#' an explicit identity, anchor, or linking design is required before a joint
#' fit is interpretable.
#'
#' @section Typical workflow:
#' 1. list valid names with [list_mfrmr_data()].
#' 2. load one dataset key with `load_mfrmr_data(name)`.
#' 3. fit a model with [fit_mfrm()] and inspect with `summary()` / `plot()`.
#'
#' @seealso [list_mfrmr_data()], [ej2021_data]
#' @examples
#' data("mfrmr_example_operational", package = "mfrmr")
#' head(mfrmr_example_operational)
#'
#' d <- load_mfrmr_data("example_operational")
#' table(d$Rater)
#' table(d$Criterion, d$Score)
#' @export
load_mfrmr_data <- function(name = c(
                            "example_core",
                            "example_bias",
                            "example_operational",
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
    example_operational = "mfrmr_example_operational",
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
