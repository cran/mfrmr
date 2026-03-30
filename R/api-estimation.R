#' Fit a many-facet Rasch model with a flexible number of facets
#'
#' This is the package entry point. It wraps `mfrm_estimate()` and defaults to
#' `method = "MML"`. Any number of facet columns can be supplied via `facets`.
#'
#' @param data A data.frame with one row per observation.
#' @param person Column name for the person (character scalar).
#' @param facets Character vector of facet column names.
#' @param score Column name for observed category score.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param weight Optional weight column name.
#' @param keep_original Keep original category values.
#' @param model `"RSM"` or `"PCM"`.
#' @param method `"MML"` (default) or `"JML"` / `"JMLE"`.
#' @param step_facet Step facet for PCM.
#' @param anchors Optional anchor table.
#' @param group_anchors Optional group-anchor table.
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#' @param positive_facets Facets with positive orientation.
#' @param anchor_policy How to handle anchor-audit issues: `"warn"` (default),
#'   `"error"`, or `"silent"`.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   anchor-audit recommendations.
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in anchor-audit recommendations.
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in anchor-audit recommendations.
#' @param quad_points Quadrature points for MML.
#' @param maxit Maximum optimizer iterations.
#' @param reltol Optimization tolerance.
#'
#' @details
#' Data must be in **long format** (one row per observed rating event).
#'
#' @section Model:
#' `fit_mfrm()` estimates the many-facet Rasch model (Linacre, 1989).
#' For a two-facet design (rater \eqn{j}, criterion \eqn{i}) the model is:
#'
#' \deqn{\ln\frac{P(X_{nij} = k)}{P(X_{nij} = k-1)} =
#'   \theta_n - \delta_j - \beta_i - \tau_k}
#'
#' where \eqn{\theta_n} is person ability, \eqn{\delta_j} rater severity,
#' \eqn{\beta_i} criterion difficulty, and \eqn{\tau_k} the \eqn{k}-th
#' Rasch-Andrich threshold.  Any number of facets may be specified via the
#' `facets` argument; each enters as an additive term in the linear
#' predictor \eqn{\eta}.
#'
#' With `model = "RSM"`, thresholds \eqn{\tau_k} are shared across all
#' levels of all facets.
#' With `model = "PCM"`, each level of `step_facet` receives its own
#' threshold vector \eqn{\tau_{i,k}} on the package's shared observed
#' score scale.
#'
#' With `method = "MML"`, person parameters are integrated out using
#' Gauss-Hermite quadrature and EAP estimates are computed post-hoc.
#' With `method = "JML"`, all parameters are estimated jointly as fixed
#' effects.  See the "Estimation methods" section of [mfrmr-package] for
#' details.
#'
#' @section Input requirements:
#' Minimum required columns are:
#' - person identifier (`person`)
#' - one or more facet identifiers (`facets`)
#' - observed score (`score`)
#'
#' Scores are treated as ordered categories.
#' If your observed categories do not start at 0, set `rating_min`/`rating_max`
#' explicitly to avoid unintended recoding assumptions.
#'
#' Supported model/estimation combinations:
#' - `model = "RSM"` with `method = "MML"` or `"JML"/"JMLE"`
#' - `model = "PCM"` with a designated `step_facet` (defaults to first facet)
#'
#' Anchor inputs are optional:
#' - `anchors` should contain facet/level/fixed-value information.
#' - `group_anchors` should contain facet/level/group/group-value information.
#' Both are normalized internally, so column names can be flexible
#' (`facet`, `level`, `anchor`, `group`, `groupvalue`, etc.).
#'
#' Anchor audit behavior:
#' - `fit_mfrm()` runs an internal anchor audit.
#' - invalid rows are removed before estimation.
#' - duplicate rows keep the last occurrence for each key.
#' - `anchor_policy` controls whether detected issues are warned, treated as
#'   errors, or kept silent.
#'
#' Facet sign orientation:
#' - facets listed in `positive_facets` are treated as `+1`
#' - all other facets are treated as `-1`
#' This affects interpretation of reported facet measures.
#'
#' @section Performance tips:
#' For exploratory work, `method = "JML"` is usually faster than `method = "MML"`,
#' but it may require a larger `maxit` to converge on larger datasets.
#'
#' For MML runs, `quad_points` is the main accuracy/speed trade-off:
#' - `quad_points = 7` is a good lightweight default for quick iteration.
#' - `quad_points = 15` gives a more stable approximation for final reporting.
#'
#' Downstream diagnostics can also be staged:
#' - use `diagnose_mfrm(fit, residual_pca = "none")` for a quick first pass
#' - add residual PCA only when you need exploratory residual-structure evidence
#'
#' Downstream diagnostics report `ModelSE` / `RealSE` columns and related
#' reliability indices. For `MML`, non-person facet `ModelSE` values are based
#' on the observed information of the marginal log-likelihood and person rows
#' use posterior SDs from EAP scoring. For `JML`, these quantities remain
#' exploratory approximations and should not be treated as equally formal.
#'
#' @section Interpreting output:
#' A typical first-pass read is:
#' 1. `fit$summary` for convergence and global fit indicators.
#' 2. `summary(fit)` for human-readable overviews.
#' 3. `diagnose_mfrm(fit)` for element-level fit, approximate
#'    separation/reliability, and warning tables.
#'
#' @section Typical workflow:
#' 1. Fit the model with `fit_mfrm(...)`.
#' 2. Validate convergence and scale structure with `summary(fit)`.
#' 3. Run [diagnose_mfrm()] and proceed to reporting with [build_apa_outputs()].
#'
#' @return
#' An object of class `mfrm_fit` (named list) with:
#' - `summary`: one-row model summary (`LogLik`, `AIC`, `BIC`, convergence)
#' - `facets$person`: person estimates (`Estimate`; plus `SD` for MML)
#' - `facets$others`: facet-level estimates for each facet
#' - `steps`: estimated threshold/step parameters
#' - `config`: resolved model configuration used for estimation
#'   (includes `config$anchor_audit`)
#' - `prep`: preprocessed data/level metadata
#' - `opt`: raw optimizer result from [stats::optim()]
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()], [build_apa_outputs()],
#'   [mfrmr_workflow_methods], [mfrmr_reporting_and_apa]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#'
#' fit <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   model = "RSM",
#'   maxit = 25
#' )
#' fit$summary
#' s_fit <- summary(fit)
#' s_fit$overview[, c("Model", "Method", "Converged")]
#' p_fit <- plot(fit, draw = FALSE)
#' class(p_fit)
#'
#' # MML is the default:
#' fit_mml <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   model = "RSM",
#'   quad_points = 7,
#'   maxit = 25
#' )
#' summary(fit_mml)
#'
#' # Next steps after fitting:
#' diag_mml <- diagnose_mfrm(fit_mml, residual_pca = "none")
#' chk <- reporting_checklist(fit_mml, diagnostics = diag_mml)
#' head(chk$checklist[, c("Section", "Item", "DraftReady")])
#' }
#' @export
fit_mfrm <- function(data,
                     person,
                     facets,
                     score,
                     rating_min = NULL,
                     rating_max = NULL,
                     weight = NULL,
                     keep_original = FALSE,
                     model = c("RSM", "PCM"),
                     method = c("MML", "JML", "JMLE"),
                     step_facet = NULL,
                     anchors = NULL,
                     group_anchors = NULL,
                     noncenter_facet = "Person",
                     dummy_facets = NULL,
                     positive_facets = NULL,
                     anchor_policy = c("warn", "error", "silent"),
                     min_common_anchors = 5L,
                     min_obs_per_element = 30,
                     min_obs_per_category = 10,
                     quad_points = 15,
                     maxit = 400,
                     reltol = 1e-6) {
  # -- input validation --
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame. Got: ", class(data)[1], ". ",
         "Convert with as.data.frame() if needed.", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("`data` has zero rows. ",
         "Supply a data.frame with at least one observation.", call. = FALSE)
  }
  if (!is.character(person) || length(person) != 1 || !nzchar(person)) {
    stop("`person` must be a single non-empty character string ",
         "naming the person column.", call. = FALSE)
  }
  if (!is.character(facets) || length(facets) == 0) {
    stop("`facets` must be a character vector of one or more facet column names.",
         call. = FALSE)
  }
  if (!is.character(score) || length(score) != 1 || !nzchar(score)) {
    stop("`score` must be a single non-empty character string ",
         "naming the score column.", call. = FALSE)
  }
  if (!is.null(weight) && (!is.character(weight) || length(weight) != 1)) {
    stop("`weight` must be NULL or a single character string ",
         "naming the weight column.", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1) {
    stop("`maxit` must be a positive integer. Got: ", deparse(maxit), ".",
         call. = FALSE)
  }
  if (!is.numeric(reltol) || length(reltol) != 1 || reltol <= 0) {
    stop("`reltol` must be a positive number. Got: ", deparse(reltol), ".",
         call. = FALSE)
  }
  if (!is.numeric(quad_points) || length(quad_points) != 1 || quad_points < 1) {
    stop("`quad_points` must be a positive integer. Got: ", deparse(quad_points), ".",
         call. = FALSE)
  }

  model <- toupper(match.arg(model))
  method_input <- toupper(match.arg(method))
  method <- ifelse(method_input == "JML", "JMLE", method_input)
  anchor_policy <- tolower(match.arg(anchor_policy))

  anchor_audit <- audit_mfrm_anchors(
    data = data,
    person = person,
    facets = facets,
    score = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight = weight,
    keep_original = keep_original,
    anchors = anchors,
    group_anchors = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  anchors <- anchor_audit$anchors
  group_anchors <- anchor_audit$group_anchors

  issue_counts <- anchor_audit$issue_counts
  issue_total <- if (is.null(issue_counts) || nrow(issue_counts) == 0) 0L else sum(issue_counts$N, na.rm = TRUE)
  if (issue_total > 0) {
    msg <- format_anchor_audit_message(anchor_audit)
    if (anchor_policy == "error") {
      stop(msg, call. = FALSE)
    } else if (anchor_policy == "warn") {
      warning(msg, call. = FALSE)
    }
  }

  fit <- mfrm_estimate(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    model = model,
    method = method,
    step_facet = step_facet,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol
  )

  fit$config$anchor_audit <- anchor_audit
  fit$config$method_input <- method_input

  class(fit) <- c("mfrm_fit", class(fit))
  fit
}

format_anchor_audit_message <- function(anchor_audit) {
  if (is.null(anchor_audit$issue_counts) || nrow(anchor_audit$issue_counts) == 0) {
    return("Anchor audit detected no issues.")
  }
  nonzero <- anchor_audit$issue_counts |>
    dplyr::filter(.data$N > 0)
  if (nrow(nonzero) == 0) {
    return("Anchor audit detected no issues.")
  }
  labels <- paste0(nonzero$Issue, "=", nonzero$N)
  paste0(
    "Anchor audit detected ", sum(nonzero$N), " issue row(s): ",
    paste(labels, collapse = "; "),
    ". Invalid rows were removed; duplicate keys keep the last row."
  )
}

summarize_linkage_by_facet <- function(df, facet) {
  by_level <- df |>
    dplyr::group_by(.data[[facet]]) |>
    dplyr::summarize(
      PersonsPerLevel = dplyr::n_distinct(.data$Person),
      Observations = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    )

  by_person <- df |>
    dplyr::group_by(.data$Person) |>
    dplyr::summarize(
      LevelsPerPerson = dplyr::n_distinct(.data[[facet]]),
      .groups = "drop"
    )

  tibble::tibble(
    Facet = facet,
    Levels = nrow(by_level),
    MinPersonsPerLevel = min(by_level$PersonsPerLevel, na.rm = TRUE),
    MedianPersonsPerLevel = stats::median(by_level$PersonsPerLevel, na.rm = TRUE),
    MinLevelsPerPerson = min(by_person$LevelsPerPerson, na.rm = TRUE),
    MedianLevelsPerPerson = stats::median(by_person$LevelsPerPerson, na.rm = TRUE)
  )
}

normalize_compare_signature <- function(fit) {
  cfg <- fit$config %||% list()
  anchor_tables <- extract_anchor_tables(cfg)

  anchor_tbl <- as.data.frame(anchor_tables$anchors %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(anchor_tbl) > 0) {
    anchor_tbl <- anchor_tbl[, intersect(c("Facet", "Level", "Anchor"), names(anchor_tbl)), drop = FALSE]
    anchor_tbl <- anchor_tbl[do.call(order, c(anchor_tbl[intersect(c("Facet", "Level"), names(anchor_tbl))], na.last = TRUE)), , drop = FALSE]
  }

  group_tbl <- as.data.frame(anchor_tables$groups %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(group_tbl) > 0) {
    group_tbl <- group_tbl[, intersect(c("Facet", "Level", "Group", "GroupValue"), names(group_tbl)), drop = FALSE]
    group_tbl <- group_tbl[do.call(order, c(group_tbl[intersect(c("Facet", "Level", "Group"), names(group_tbl))], na.last = TRUE)), , drop = FALSE]
  }

  list(
    model = as.character(cfg$model %||% NA_character_),
    method = as.character(cfg$method %||% NA_character_),
    person = as.character(cfg$person_col %||% NA_character_),
    facets = sort(as.character(cfg$facet_cols %||% fit$prep$facet_names %||% character(0))),
    score = as.character(cfg$score_col %||% NA_character_),
    weight = as.character(cfg$weight_col %||% NA_character_),
    step_facet = as.character(cfg$step_facet %||% NA_character_),
    noncenter_facet = as.character(cfg$noncenter_facet %||% "Person"),
    dummy_facets = sort(as.character(cfg$dummy_facets %||% character(0))),
    positive_facets = sort(as.character(cfg$positive_facets %||% character(0))),
    anchors = anchor_tbl,
    group_anchors = group_tbl
  )
}

same_signature_component <- function(x, y) {
  if (is.data.frame(x) || is.data.frame(y)) {
    return(identical(
      as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE),
      as.data.frame(y %||% data.frame(), stringsAsFactors = FALSE)
    ))
  }
  identical(x, y)
}

audit_compare_mfrm_nesting <- function(fits, labels) {
  if (length(fits) != 2L) {
    return(list(
      eligible = FALSE,
      reason = "Likelihood-ratio tests are only audited for two-model comparisons.",
      simpler = NA_character_,
      complex = NA_character_,
      relation = "unsupported"
    ))
  }

  sigs <- lapply(fits, normalize_compare_signature)
  lbls <- as.character(labels)

  same_components <- c(
    person = same_signature_component(sigs[[1]]$person, sigs[[2]]$person),
    facets = same_signature_component(sigs[[1]]$facets, sigs[[2]]$facets),
    score = same_signature_component(sigs[[1]]$score, sigs[[2]]$score),
    weight = same_signature_component(sigs[[1]]$weight, sigs[[2]]$weight),
    noncenter_facet = same_signature_component(sigs[[1]]$noncenter_facet, sigs[[2]]$noncenter_facet),
    dummy_facets = same_signature_component(sigs[[1]]$dummy_facets, sigs[[2]]$dummy_facets),
    positive_facets = same_signature_component(sigs[[1]]$positive_facets, sigs[[2]]$positive_facets),
    anchors = same_signature_component(sigs[[1]]$anchors, sigs[[2]]$anchors),
    group_anchors = same_signature_component(sigs[[1]]$group_anchors, sigs[[2]]$group_anchors)
  )

  if (!all(same_components)) {
    mismatch <- names(same_components)[!same_components]
    return(list(
      eligible = FALSE,
      reason = paste0(
        "Models differ in structural comparison settings: ",
        paste(mismatch, collapse = ", "),
        "."
      ),
      simpler = NA_character_,
      complex = NA_character_,
      relation = "unsupported"
    ))
  }

  model_pair <- toupper(c(sigs[[1]]$model, sigs[[2]]$model))
  if (identical(model_pair[1], model_pair[2])) {
    return(list(
      eligible = FALSE,
      reason = "Both fits use the same model family, so there is no supported nested restriction to test.",
      simpler = lbls[1],
      complex = lbls[2],
      relation = "same_model"
    ))
  }

  if (setequal(model_pair, c("RSM", "PCM"))) {
    idx_simple <- which(model_pair == "RSM")[1]
    idx_complex <- which(model_pair == "PCM")[1]
    step_facet <- sigs[[idx_complex]]$step_facet
    if (!is.na(step_facet) && nzchar(step_facet)) {
      return(list(
        eligible = TRUE,
        reason = paste0(
          "Supported nesting audit passed: shared design/constraints with RSM nested inside PCM on step facet '",
          step_facet,
          "'."
        ),
        simpler = lbls[idx_simple],
        complex = lbls[idx_complex],
        relation = "RSM_in_PCM"
      ))
    }
  }

  list(
    eligible = FALSE,
    reason = "Automatic nesting audit currently supports only RSM nested inside PCM under shared design and constraints.",
    simpler = NA_character_,
    complex = NA_character_,
    relation = "unsupported"
  )
}

#' Summarize MFRM input data (TAM-style descriptive snapshot)
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param keep_original Keep original category values.
#' @param include_person_facet If `TRUE`, include person-level rows in
#'   `facet_level_summary`.
#' @param include_agreement If `TRUE`, include an observed-score inter-rater
#'   agreement bundle (summary/pairs/settings) in the output.
#' @param rater_facet Optional rater facet name used for agreement summaries.
#'   If `NULL`, inferred from facet names.
#' @param context_facets Optional facets used to define matched contexts for
#'   agreement. If `NULL`, all remaining facets (including `Person`) are used.
#' @param agreement_top_n Optional maximum number of agreement pair rows.
#'
#' @details
#' This function provides a compact descriptive bundle similar to the
#' pre-fit summaries commonly checked in TAM workflows:
#' sample size, score distribution, per-facet coverage, and linkage counts.
#' `psych::describe()` is used for numeric descriptives of score and weight.
#'
#' **Key data-quality checks to perform before fitting:**
#' - *Sparse categories*: any score category with fewer than 10 weighted
#'   observations may produce unstable threshold estimates
#'   (Linacre, 2002).  Consider collapsing adjacent categories.
#' - *Unlinked elements*: if a facet level has zero overlap with one or
#'   more levels of another facet, the design is disconnected and
#'   parameters cannot be placed on a common scale.  Check
#'   `linkage_summary` for low connectivity.
#' - *Extreme scores*: persons or facet levels with all-minimum or
#'   all-maximum scores yield infinite logit estimates under JML;
#'   they are handled via Bayesian shrinkage under MML.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `overview`: confirms sample size, facet count, and category span.
#'   The `MinWeightedN` column shows the smallest weighted observation
#'   count across all facet levels; values below 30 may lead to
#'   unstable parameter estimates.
#' - `missing_by_column`: identifies immediate data-quality risks.
#'   Any non-zero count warrants investigation before fitting.
#' - `score_distribution`: checks sparse/unused score categories.
#'   Balanced usage across categories is ideal; heavily skewed
#'   distributions may compress the measurement range.
#' - `facet_level_summary` and `linkage_summary`: checks per-level
#'   support and person-facet connectivity.  Low linkage ratios
#'   indicate sparse or disconnected design blocks.
#' - `agreement`: optional observed inter-rater consistency summary
#'   (exact agreement, correlation, mean differences per rater pair).
#'
#' @section Typical workflow:
#' 1. Run `describe_mfrm_data()` on long-format input.
#' 2. Review `summary(ds)` and `plot(ds, ...)`.
#' 3. Resolve missingness/sparsity issues before [fit_mfrm()].
#'
#' @return A list of class `mfrm_data_description` with:
#' - `overview`: one-row run-level summary
#' - `missing_by_column`: missing counts in selected input columns
#' - `score_descriptives`: output from [psych::describe()] for score
#' - `weight_descriptives`: output from [psych::describe()] for weight
#' - `score_distribution`: weighted and raw score frequencies
#' - `facet_level_summary`: per-level usage and score summaries
#' - `linkage_summary`: person-facet connectivity diagnostics
#' - `agreement`: observed-score inter-rater agreement bundle
#'
#' @seealso [fit_mfrm()], [audit_mfrm_anchors()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score"
#' )
#' s_ds <- summary(ds)
#' s_ds$overview
#' p_ds <- plot(ds, draw = FALSE)
#' class(p_ds)
#' @export
describe_mfrm_data <- function(data,
                               person,
                               facets,
                               score,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               include_person_facet = FALSE,
                               include_agreement = TRUE,
                               rater_facet = NULL,
                               context_facets = NULL,
                               agreement_top_n = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original
  )

  df <- prep$data |>
    dplyr::mutate(
      Person = as.character(.data$Person),
      dplyr::across(dplyr::all_of(prep$facet_names), as.character)
    )

  selected_cols <- unique(c(person, facets, score, if (!is.null(weight)) weight))
  missing_by_column <- tibble::tibble(
    Column = selected_cols,
    Missing = vapply(selected_cols, function(col) sum(is.na(data[[col]])), integer(1))
  )

  score_desc <- psych::describe(df$Score, fast = TRUE)
  weight_desc <- psych::describe(df$Weight, fast = TRUE)

  total_weight <- sum(df$Weight, na.rm = TRUE)
  score_distribution <- df |>
    dplyr::group_by(.data$Score) |>
    dplyr::summarize(
      RawN = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      Percent = ifelse(total_weight > 0, 100 * .data$WeightedN / total_weight, NA_real_),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Score)

  report_facets <- prep$facet_names
  if (isTRUE(include_person_facet)) {
    report_facets <- c("Person", report_facets)
  }

  facet_level_summary <- purrr::map_dfr(report_facets, function(facet) {
    df |>
      dplyr::group_by(.data[[facet]]) |>
      dplyr::summarize(
        RawN = dplyr::n(),
        WeightedN = sum(.data$Weight, na.rm = TRUE),
        MeanScore = weighted_mean(.data$Score, .data$Weight),
        SDScore = stats::sd(.data$Score, na.rm = TRUE),
        MinScore = min(.data$Score, na.rm = TRUE),
        MaxScore = max(.data$Score, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::rename(Level = dplyr::all_of(facet)) |>
      dplyr::mutate(
        Facet = facet,
        Level = as.character(.data$Level),
        .before = 1
      )
  })

  linkage_summary <- if (length(prep$facet_names) == 0) {
    tibble::tibble()
  } else {
    purrr::map_dfr(prep$facet_names, function(facet) summarize_linkage_by_facet(df, facet))
  }

  agreement_bundle <- list(
    summary = data.frame(),
    pairs = data.frame(),
    settings = list(
      included = FALSE,
      rater_facet = NA_character_,
      context_facets = character(0),
      expected_exact_from_model = FALSE,
      top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n))
    )
  )

  if (isTRUE(include_agreement) && length(prep$facet_names) > 0) {
    known_facets <- c("Person", prep$facet_names)
    if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
      rater_facet <- infer_default_rater_facet(prep$facet_names)
    } else {
      rater_facet <- as.character(rater_facet[1])
    }
    if (is.null(rater_facet) || !rater_facet %in% known_facets) {
      stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
    }
    if (identical(rater_facet, "Person")) {
      stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
    }

    if (is.null(context_facets)) {
      facet_cols <- known_facets
      resolved_context <- setdiff(facet_cols, rater_facet)
    } else {
      context_facets <- unique(as.character(context_facets))
      unknown <- setdiff(context_facets, known_facets)
      if (length(unknown) > 0) {
        stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "))
      }
      resolved_context <- setdiff(context_facets, rater_facet)
      if (length(resolved_context) == 0) {
        stop("`context_facets` must include at least one facet different from `rater_facet`.")
      }
      facet_cols <- c(rater_facet, resolved_context)
    }

    obs_agreement <- df |>
      dplyr::select(dplyr::all_of(unique(c("Person", prep$facet_names, "Score", "Weight")))) |>
      dplyr::rename(Observed = "Score")

    agreement <- calc_interrater_agreement(
      obs_df = obs_agreement,
      facet_cols = facet_cols,
      rater_facet = rater_facet,
      res = NULL
    )
    agreement_pairs <- as.data.frame(agreement$pairs, stringsAsFactors = FALSE)
    if (!is.null(agreement_top_n) && nrow(agreement_pairs) > 0) {
      agreement_pairs <- agreement_pairs |>
        dplyr::slice_head(n = max(1L, as.integer(agreement_top_n)))
    }

    agreement_bundle <- list(
      summary = as.data.frame(agreement$summary, stringsAsFactors = FALSE),
      pairs = agreement_pairs,
      settings = list(
        included = TRUE,
        rater_facet = rater_facet,
        context_facets = resolved_context,
        expected_exact_from_model = FALSE,
        top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n))
      )
    )
  }

  overview <- tibble::tibble(
    Observations = nrow(df),
    TotalWeight = total_weight,
    Persons = length(prep$levels$Person),
    Facets = length(prep$facet_names),
    Categories = prep$rating_max - prep$rating_min + 1,
    RatingMin = prep$rating_min,
    RatingMax = prep$rating_max
  )

  out <- list(
    overview = overview,
    missing_by_column = missing_by_column,
    score_descriptives = score_desc,
    weight_descriptives = weight_desc,
    score_distribution = score_distribution,
    facet_level_summary = facet_level_summary,
    linkage_summary = linkage_summary,
    agreement = agreement_bundle
  )
  class(out) <- c("mfrm_data_description", class(out))
  out
}

#' @export
print.mfrm_data_description <- function(x, ...) {
  cat("mfrm data description\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    print(x$overview, row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(x$score_distribution, row.names = FALSE)
  }
  if (!is.null(x$agreement$summary) && nrow(x$agreement$summary) > 0) {
    cat("\nInter-rater agreement (observed)\n")
    print(x$agreement$summary, row.names = FALSE)
  }
  invisible(x)
}

#' Summarize a data-description object
#'
#' @param object Output from [describe_mfrm_data()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in preview blocks.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is intended as a compact pre-fit quality snapshot for
#' manuscripts and analysis logs.
#'
#' @section Interpreting output:
#' Recommended read order:
#' - `overview`: sample size, persons/facets/categories.
#' - `missing`: missingness hotspots by selected input columns.
#' - `score_distribution`: category usage balance.
#' - `facet_overview`: coverage per facet (minimum/maximum weighted counts).
#' - `agreement`: observed-score inter-rater agreement (when available).
#'
#' Very low `MinWeightedN` in `facet_overview` is a practical warning for
#' unstable downstream facet estimates.
#'
#' @section Typical workflow:
#' 1. Run [describe_mfrm_data()] on raw long-format data.
#' 2. Inspect `summary(ds)` before model fitting.
#' 3. Resolve sparse/missing issues, then run [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_data_description`.
#' @seealso [describe_mfrm_data()], [summary.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(ds)
#' @export
summary.mfrm_data_description <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
  missing_tbl <- as.data.frame(object$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(missing_tbl) > 0 && all(c("Column", "Missing") %in% names(missing_tbl))) {
    missing_tbl <- missing_tbl |>
      dplyr::arrange(dplyr::desc(.data$Missing), .data$Column) |>
      dplyr::slice_head(n = top_n)
  }

  score_dist <- as.data.frame(object$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(score_dist) > 0) {
    score_dist <- utils::head(score_dist, n = top_n)
  }

  facet_tbl <- as.data.frame(object$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
  facet_overview <- data.frame()
  if (nrow(facet_tbl) > 0 && all(c("Facet", "Level", "WeightedN") %in% names(facet_tbl))) {
    facet_overview <- facet_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n_distinct(.data$Level),
        TotalWeightedN = sum(.data$WeightedN, na.rm = TRUE),
        MeanWeightedN = mean(.data$WeightedN, na.rm = TRUE),
        MinWeightedN = min(.data$WeightedN, na.rm = TRUE),
        MaxWeightedN = max(.data$WeightedN, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Facet) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  agreement_tbl <- as.data.frame(object$agreement$summary %||% data.frame(), stringsAsFactors = FALSE)
  notes <- if (nrow(missing_tbl) > 0 && any(suppressWarnings(as.numeric(missing_tbl$Missing)) > 0, na.rm = TRUE)) {
    "Missing values were detected in one or more input columns."
  } else {
    "No missing values were detected in selected input columns."
  }

  out <- list(
    overview = overview,
    missing = missing_tbl,
    score_distribution = score_dist,
    facet_overview = facet_overview,
    agreement = agreement_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_data_description"
  out
}

#' @export
print.summary.mfrm_data_description <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("mfrm Data Description Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$missing) && nrow(x$missing) > 0) {
    cat("\nMissing by column\n")
    print(round_numeric_df(as.data.frame(x$missing), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(round_numeric_df(as.data.frame(x$score_distribution), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_overview) && nrow(x$facet_overview) > 0) {
    cat("\nFacet coverage\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$agreement) && nrow(x$agreement) > 0) {
    cat("\nInter-rater agreement\n")
    print(round_numeric_df(as.data.frame(x$agreement), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$notes) && nzchar(x$notes)) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Plot a data-description object
#'
#' @param x Output from [describe_mfrm_data()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"score_distribution"`, `"facet_levels"`, or `"missing"`.
#' @param main Optional title override.
#' @param palette Optional named colors (`score`, `facet`, `missing`).
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method draws quick pre-fit quality views from [describe_mfrm_data()]:
#' - score distribution balance
#' - facet-level structure size
#' - missingness by selected columns
#'
#' @section Interpreting output:
#' - `"score_distribution"`: bar chart of weighted observation counts per
#'   score category.  Y-axis is `WeightedN` (sum of weights for each
#'   category).  Categories with very few observations (< 10) may produce
#'   unstable threshold estimates.  A roughly uniform or unimodal
#'   distribution is ideal; heavy floor/ceiling effects compress the
#'   measurement range.
#' - `"facet_levels"`: bar chart showing the number of distinct levels
#'   per facet.  Useful for verifying that the design structure matches
#'   expectations (e.g., expected number of raters or criteria).  Very
#'   large numbers of levels increase computation time and may require
#'   higher `maxit` in [fit_mfrm()].
#' - `"missing"`: bar chart of missing-value counts per input column.
#'   Columns with non-zero counts should be investigated before
#'   fitting---rows with missing scores, persons, or facet IDs are
#'   dropped during estimation.
#'
#' @section Typical workflow:
#' 1. Run [describe_mfrm_data()] before fitting.
#' 2. Inspect `summary(ds)` and `plot(ds, type = "missing")`.
#' 3. Check category/facet balance with other plot types.
#' 4. Fit model after resolving obvious data issues.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [describe_mfrm_data()], `plot()`
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(ds, draw = FALSE)
#' @export
plot.mfrm_data_description <- function(x,
                                       y = NULL,
                                       type = c("score_distribution", "facet_levels", "missing"),
                                       main = NULL,
                                       palette = NULL,
                                       label_angle = 45,
                                       draw = TRUE,
                                       ...) {
  type <- match.arg(tolower(as.character(type[1])), c("score_distribution", "facet_levels", "missing"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      score = "#2b8cbe",
      facet = "#31a354",
      missing = "#756bb1"
    )
  )

  if (type == "score_distribution") {
    tbl <- as.data.frame(x$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Score", "WeightedN") %in% names(tbl))) {
      stop("Score distribution is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$WeightedN)),
        labels = as.character(tbl$Score),
        col = pal["score"],
        main = if (is.null(main)) "Score distribution" else as.character(main[1]),
        ylab = "Weighted N",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_description",
      list(plot = "score_distribution", table = tbl)
    )))
  }

  if (type == "facet_levels") {
    tbl <- as.data.frame(x$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "Level") %in% names(tbl))) {
      stop("Facet level summary is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
    }
    agg <- tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(Levels = dplyr::n_distinct(.data$Level), .groups = "drop") |>
      dplyr::arrange(.data$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(agg$Levels)),
        labels = as.character(agg$Facet),
        col = pal["facet"],
        main = if (is.null(main)) "Facet levels" else as.character(main[1]),
        ylab = "Levels",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_description",
      list(plot = "facet_levels", table = as.data.frame(agg, stringsAsFactors = FALSE))
    )))
  }

  tbl <- as.data.frame(x$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Column", "Missing") %in% names(tbl))) {
    stop("Missing-by-column table is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
  }
  if (isTRUE(draw)) {
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl$Missing)),
      labels = as.character(tbl$Column),
      col = pal["missing"],
      main = if (is.null(main)) "Missing values by column" else as.character(main[1]),
      ylab = "Missing",
      label_angle = label_angle,
      mar_bottom = 8.0
    )
  }
  invisible(new_mfrm_plot_data(
    "data_description",
    list(plot = "missing", table = tbl)
  ))
}

#' Audit and normalize anchor/group-anchor tables
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param anchors Optional anchor table (Facet, Level, Anchor).
#' @param group_anchors Optional group-anchor table
#'   (Facet, Level, Group, GroupValue).
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param keep_original Keep original category values.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   recommendations (default `5`).
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in recommendations (default `30`).
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in recommendations (default `10`).
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#'
#' @details
#' **Anchoring** (also called "fixing" or scale linking) constrains selected
#' parameter estimates to pre-specified values, placing the current
#' analysis on a previously established scale.  This is essential when
#' comparing results across administrations, linking test forms, or
#' monitoring rater drift over time.
#'
#' This function applies the same preprocessing and key-resolution rules
#' as `fit_mfrm()`, but returns an audit object so constraints can be
#' checked *before* estimation.  Running the audit first helps avoid
#' estimation failures caused by misspecified or data-incompatible
#' anchors.
#'
#' **Anchor types:**
#' - *Direct anchors* fix individual element measures to specific logit
#'   values (e.g., Rater R1 anchored at 0.35 logits).
#' - *Group anchors* constrain the mean of a set of elements to a
#'   target value, allowing individual elements to vary freely around
#'   that mean.
#' - When both types overlap for the same element, the direct anchor
#'   takes precedence.
#'
#' **Design checks** verify that each anchored element has at least
#' `min_obs_per_element` weighted observations (default 30) and each
#' score category has at least `min_obs_per_category` (default 10).
#' These thresholds follow standard Rasch sample-size recommendations
#' (Linacre, 1994).
#'
#' @section Interpreting output:
#' - `issue_counts`/`issues`: concrete data or specification problems.
#' - `facet_summary`: constraint coverage by facet.
#' - `design_checks`: whether anchor targets have enough observations.
#' - `recommendations`: action items before estimation.
#'
#' @section Typical workflow:
#' 1. Build candidate anchors (e.g., with [make_anchor_table()]).
#' 2. Run `audit_mfrm_anchors(...)`.
#' 3. Resolve issues, then fit with [fit_mfrm()].
#'
#' @return A list of class `mfrm_anchor_audit` with:
#' - `anchors`: cleaned anchor table used by estimation
#' - `group_anchors`: cleaned group-anchor table used by estimation
#' - `facet_summary`: counts of levels, constrained levels, and free levels
#' - `design_checks`: observation-count checks by level/category
#' - `thresholds`: active threshold settings used for recommendations
#' - `issue_counts`: issue-type counts
#' - `issues`: list of issue tables
#' - `recommendations`: package-native anchor guidance strings
#'
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [make_anchor_table()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#'
#' anchors <- data.frame(
#'   Facet = c("Rater", "Rater"),
#'   Level = c("R1", "R1"),
#'   Anchor = c(0, 0.1),
#'   stringsAsFactors = FALSE
#' )
#' aud <- audit_mfrm_anchors(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   anchors = anchors
#' )
#' aud$issue_counts
#' summary(aud)
#' p_aud <- plot(aud, draw = FALSE)
#' class(p_aud)
#' @export
audit_mfrm_anchors <- function(data,
                               person,
                               facets,
                               score,
                               anchors = NULL,
                               group_anchors = NULL,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               min_common_anchors = 5L,
                               min_obs_per_element = 30,
                               min_obs_per_category = 10,
                               noncenter_facet = "Person",
                               dummy_facets = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original
  )

  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)

  audit <- audit_anchor_tables(
    prep = prep,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )
  class(audit) <- c("mfrm_anchor_audit", class(audit))
  audit
}

#' @export
print.mfrm_anchor_audit <- function(x, ...) {
  issue_total <- if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) sum(x$issue_counts$N) else 0
  cat("mfrm anchor audit\n")
  cat("  issue rows: ", issue_total, "\n", sep = "")

  if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
    nonzero <- x$issue_counts |>
      dplyr::filter(.data$N > 0)
    if (nrow(nonzero) > 0) {
      cat("\nIssue counts\n")
      print(nonzero, row.names = FALSE)
    }
  }

  if (!is.null(x$facet_summary) && nrow(x$facet_summary) > 0) {
    cat("\nFacet summary\n")
    print(x$facet_summary, row.names = FALSE)
  }

  if (!is.null(x$design_checks) &&
      !is.null(x$design_checks$level_observation_summary) &&
      nrow(x$design_checks$level_observation_summary) > 0) {
    cat("\nLevel observation summary\n")
    print(x$design_checks$level_observation_summary, row.names = FALSE)
  }

  invisible(x)
}

#' Summarize an anchor-audit object
#'
#' @param object Output from [audit_mfrm_anchors()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in issue previews.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary provides a compact pre-estimation audit of anchor and
#' group-anchor specifications.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `issue_counts`: primary triage table (non-zero issues first).
#' - `facet_summary`: anchored/grouped/free-level balance by facet.
#' - `level_observation_summary` and `category_counts`: sparse-cell diagnostics.
#' - `recommendations`: concrete remediation suggestions.
#'
#' If `issue_counts` is non-empty, treat anchor constraints as provisional and
#' resolve issues before final estimation.
#'
#' @section Typical workflow:
#' 1. Run [audit_mfrm_anchors()] with intended anchors/group anchors.
#' 2. Review `summary(aud)` and recommendations.
#' 3. Revise anchor tables, then call [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_anchor_audit`.
#' @seealso [audit_mfrm_anchors()], [fit_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' aud <- audit_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(aud)
#' @export
summary.mfrm_anchor_audit <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  issue_counts <- as.data.frame(object$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(issue_counts) > 0 && all(c("Issue", "N") %in% names(issue_counts))) {
    issue_counts <- issue_counts |>
      dplyr::filter(.data$N > 0) |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue) |>
      dplyr::slice_head(n = top_n)
  }

  facet_summary <- as.data.frame(object$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  level_summary <- as.data.frame(object$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  category_summary <- as.data.frame(object$design_checks$category_counts %||% data.frame(), stringsAsFactors = FALSE)

  recommendations <- as.character(object$recommendations %||% character(0))
  if (length(recommendations) > top_n) {
    recommendations <- recommendations[seq_len(top_n)]
  }

  notes <- if (nrow(issue_counts) > 0) {
    "Anchor-audit issues were detected. Review issue counts and recommendations."
  } else {
    "No anchor-table issue rows were detected."
  }

  out <- list(
    issue_counts = issue_counts,
    facet_summary = facet_summary,
    level_observation_summary = level_summary,
    category_counts = category_summary,
    recommendations = recommendations,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_anchor_audit"
  out
}

#' @export
print.summary.mfrm_anchor_audit <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Anchor Audit Summary\n")
  if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
    cat("\nIssue counts\n")
    print(round_numeric_df(as.data.frame(x$issue_counts), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_summary) && nrow(x$facet_summary) > 0) {
    cat("\nFacet summary\n")
    print(round_numeric_df(as.data.frame(x$facet_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$level_observation_summary) && nrow(x$level_observation_summary) > 0) {
    cat("\nLevel observation summary\n")
    print(round_numeric_df(as.data.frame(x$level_observation_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$category_counts) && nrow(x$category_counts) > 0) {
    cat("\nCategory counts\n")
    print(round_numeric_df(as.data.frame(x$category_counts), digits = digits), row.names = FALSE)
  }
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations\n")
    for (line in x$recommendations) cat(" - ", line, "\n", sep = "")
  }
  if (!is.null(x$notes) && nzchar(x$notes)) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Plot an anchor-audit object
#'
#' @param x Output from [audit_mfrm_anchors()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"issue_counts"`, `"facet_constraints"`,
#'   or `"level_observations"`.
#' @param main Optional title override.
#' @param palette Optional named colors.
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Base-R visualization helper for anchor audit outputs.
#'
#' @section Interpreting output:
#' - `"issue_counts"`: volume of each issue class.
#' - `"facet_constraints"`: anchored/grouped/free mix by facet.
#' - `"level_observations"`: observation support across levels.
#'
#' @section Typical workflow:
#' 1. Run [audit_mfrm_anchors()].
#' 2. Start with `plot(aud, type = "issue_counts")`.
#' 3. Inspect constraint and support plots before fitting.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [audit_mfrm_anchors()], [make_anchor_table()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' aud <- audit_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(aud, draw = FALSE)
#' @export
plot.mfrm_anchor_audit <- function(x,
                                   y = NULL,
                                   type = c("issue_counts", "facet_constraints", "level_observations"),
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45,
                                   draw = TRUE,
                                   ...) {
  type <- match.arg(tolower(as.character(type[1])), c("issue_counts", "facet_constraints", "level_observations"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      issues = "#cb181d",
      anchored = "#756bb1",
      grouped = "#9ecae1",
      levels = "#2b8cbe"
    )
  )

  if (type == "issue_counts") {
    tbl <- as.data.frame(x$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Issue", "N") %in% names(tbl))) {
      stop("Issue-count table is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    tbl <- tbl |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$N)),
        labels = as.character(tbl$Issue),
        col = pal["issues"],
        main = if (is.null(main)) "Anchor-audit issue counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_audit",
      list(plot = "issue_counts", table = tbl)
    )))
  }

  if (type == "facet_constraints") {
    tbl <- as.data.frame(x$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "AnchoredLevels", "GroupedLevels", "FreeLevels") %in% names(tbl))) {
      stop("Facet summary with constraint columns is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar), add = TRUE)
      mar <- old_mar
      mar[1] <- max(mar[1], 8.8)
      graphics::par(mar = mar)
      mat <- rbind(
        Anchored = suppressWarnings(as.numeric(tbl$AnchoredLevels)),
        Grouped = suppressWarnings(as.numeric(tbl$GroupedLevels)),
        Free = suppressWarnings(as.numeric(tbl$FreeLevels))
      )
      mids <- graphics::barplot(
        height = mat,
        beside = FALSE,
        names.arg = FALSE,
        col = c(pal["anchored"], pal["grouped"], "#d9d9d9"),
        border = "white",
        ylab = "Levels",
        main = if (is.null(main)) "Constraint profile by facet" else as.character(main[1])
      )
      draw_rotated_x_labels(
        at = mids,
        labels = as.character(tbl$Facet),
        srt = label_angle,
        cex = 0.82,
        line_offset = 0.085
      )
      graphics::legend(
        "topright",
        legend = c("Anchored", "Grouped", "Free"),
        fill = c(pal["anchored"], pal["grouped"], "#d9d9d9"),
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_audit",
      list(plot = "facet_constraints", table = tbl)
    )))
  }

  tbl <- as.data.frame(x$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Facet", "MinObsPerLevel") %in% names(tbl))) {
    stop("Level observation summary is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
  }
  if (isTRUE(draw)) {
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl$MinObsPerLevel)),
      labels = as.character(tbl$Facet),
      col = pal["levels"],
      main = if (is.null(main)) "Minimum observations per level" else as.character(main[1]),
      ylab = "Min observations",
      label_angle = label_angle,
      mar_bottom = 8.0
    )
    if ("RecommendedMinObs" %in% names(tbl)) {
      r <- suppressWarnings(as.numeric(tbl$RecommendedMinObs))
      r <- r[is.finite(r)]
      if (length(r) > 0) graphics::abline(h = unique(r)[1], lty = 2, col = "gray45")
    }
  }
  invisible(new_mfrm_plot_data(
    "anchor_audit",
    list(plot = "level_observations", table = tbl)
  ))
}

#' Build an anchor table from fitted estimates
#'
#' @param fit Output from [fit_mfrm()].
#' @param facets Optional subset of facets to include.
#' @param include_person Include person estimates as anchors.
#' @param digits Rounding digits for anchor values.
#'
#' @details
#' This function exports estimated facet parameters as an anchor table
#' for use in subsequent calibrations.  This is the standard approach
#' for **linking** across administrations: a reference
#' run establishes the measurement scale, and anchored re-analyses
#' place new data on that same scale.
#'
#' Anchor values should be exported from a well-fitting reference run
#' with adequate sample size.  If the reference model has convergence
#' issues or large misfit, the exported anchors may propagate
#' instability.  Re-run [audit_mfrm_anchors()] on the receiving data
#' to verify compatibility before estimation.
#'
#' The `digits` parameter controls rounding precision.  Use at least 4
#' digits for research applications; excessive rounding (e.g., 1 digit)
#' can introduce avoidable calibration error.
#'
#' @section Interpreting output:
#' - `Facet`: facet name to be anchored in later runs.
#' - `Level`: specific element/level name inside that facet.
#' - `Anchor`: fixed logit value (rounded by `digits`).
#'
#' @section Typical workflow:
#' 1. Fit a reference run with [fit_mfrm()].
#' 2. Export anchors with `make_anchor_table(fit)`.
#' 3. Pass selected rows back into `fit_mfrm(..., anchors = ...)`.
#'
#' @return A data.frame with `Facet`, `Level`, and `Anchor`.
#' @seealso [fit_mfrm()], [audit_mfrm_anchors()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' anchors_tbl <- make_anchor_table(fit)
#' head(anchors_tbl)
#' summary(anchors_tbl$Anchor)
#' @export
make_anchor_table <- function(fit,
                              facets = NULL,
                              include_person = FALSE,
                              digits = 6) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }

  digits <- max(0L, as.integer(digits))
  out <- tibble::tibble()

  if (isTRUE(include_person) && !is.null(fit$facets$person) && nrow(fit$facets$person) > 0) {
    per <- tibble::as_tibble(fit$facets$person)
    if ("Person" %in% names(per) && "Estimate" %in% names(per)) {
      out <- dplyr::bind_rows(
        out,
        per |>
          dplyr::transmute(
            Facet = "Person",
            Level = as.character(.data$Person),
            Anchor = round(as.numeric(.data$Estimate), digits = digits)
          )
      )
    }
  }

  others <- tibble::as_tibble(fit$facets$others)
  if (nrow(others) > 0 && all(c("Facet", "Level", "Estimate") %in% names(others))) {
    out <- dplyr::bind_rows(
      out,
      others |>
        dplyr::transmute(
          Facet = as.character(.data$Facet),
          Level = as.character(.data$Level),
          Anchor = round(as.numeric(.data$Estimate), digits = digits)
        )
    )
  }

  if (!is.null(facets)) {
    keep <- as.character(facets)
    out <- out |>
      dplyr::filter(.data$Facet %in% keep)
  }

  out |>
    dplyr::arrange(.data$Facet, .data$Level)
}

#' Compute diagnostics for an `mfrm_fit` object
#'
#' @param fit Output from [fit_mfrm()].
#' @param interaction_pairs Optional list of facet pairs.
#' @param top_n_interactions Number of top interactions.
#' @param whexact Use exact ZSTD transformation.
#' @param residual_pca Residual PCA mode: `"none"`, `"overall"`, `"facet"`, or `"both"`.
#' @param pca_max_factors Maximum number of PCA factors to retain per matrix.
#'
#' @details
#' This function computes a diagnostic bundle used by downstream reporting.
#' It calculates element-level fit statistics, approximate facet
#' separation/reliability summaries, residual-based QC diagnostics, and
#' optionally residual PCA for
#' exploratory residual-structure screening.
#'
#' **Key fit statistics computed for each element:**
#' - **Infit MnSq**: information-weighted mean-square residual; sensitive
#'   to on-target misfitting patterns.  Expected value = 1.0.
#' - **Outfit MnSq**: unweighted mean-square residual; sensitive to
#'   off-target outliers.  Expected value = 1.0.
#' - **ZSTD**: Wilson-Hilferty cube-root transformation of MnSq to an
#'   approximate standard normal deviate.
#' - **PTMEA**: point-measure correlation (item-rest correlation in MFRM
#'   context); positive values confirm alignment with the latent trait.
#'
#' **Misfit flagging guidelines (Bond & Fox, 2015):**
#' - MnSq < 0.5: overfit (too predictable; may inflate reliability)
#' - MnSq 0.5--1.5: productive for measurement
#' - MnSq > 1.5: underfit (noise degrades measurement)
#' - \eqn{|\mathrm{ZSTD}| > 2}: statistically significant misfit (5\%)
#'
#' When Infit and Outfit disagree, Infit is generally more informative
#' because it downweights extreme observations.  Large Outfit with
#' acceptable Infit typically indicates a few outlying responses rather
#' than systematic misfit.
#'
#' `interaction_pairs` controls which facet interactions are summarized.
#' Each element can be:
#' - a length-2 character vector such as `c("Rater", "Criterion")`, or
#' - omitted (`NULL`) to let the function select top interactions automatically.
#'
#' Residual PCA behavior:
#' - `"none"`: skip PCA (fastest; recommended for initial exploration)
#' - `"overall"`: compute overall residual PCA across all facets
#' - `"facet"`: compute facet-specific residual PCA for each facet
#' - `"both"`: compute both overall and facet-specific PCA
#'
#' Overall PCA examines the person \eqn{\times} combined-facet residual
#' matrix; facet-specific PCA examines person \eqn{\times} facet-level
#' matrices. These summaries are exploratory screens for residual
#' structure, not standalone proofs for or against unidimensionality.
#' Facet-specific PCA can help localise where a stronger residual signal
#' is concentrated.
#'
#' @section Reading key components:
#' Practical interpretation often starts with:
#' - `overall_fit`: global infit/outfit and degrees of freedom.
#' - `reliability`: facet-level model/real separation and reliability. `MML`
#'   uses model-based `ModelSE` values where available; `JML` keeps these
#'   quantities as exploratory approximations.
#' - `fit`: element-level misfit scan (`Infit`, `Outfit`, `ZSTD`).
#' - `unexpected`, `fair_average`, `displacement`: targeted QC bundles.
#' - `approximation_notes`: method notes for SE/CI/reliability summaries.
#'
#' @section Interpreting output:
#' Start with `overall_fit` and `reliability`, then move to element-level
#' diagnostics (`fit`) and targeted bundles (`unexpected`, `fair_average`,
#' `displacement`, `interrater`, `facets_chisq`).
#'
#' Consistent signals across multiple components are typically more robust than
#' a single isolated warning.  For example, an element flagged for both high
#' Outfit and high displacement is more concerning than one flagged on a
#' single criterion.
#'
#' `SE` is kept as a compatibility alias for `ModelSE`. `RealSE` is a
#' fit-adjusted companion defined as `ModelSE * sqrt(max(Infit, 1))`.
#' Reliability tables report model and fit-adjusted bounds from observed
#' variance, error variance, and true variance; `JML` entries should still be
#' treated as exploratory.
#'
#' @section Typical workflow:
#' 1. Run `diagnose_mfrm(fit, residual_pca = "none")` for baseline diagnostics.
#' 2. Inspect `summary(diag)` and targeted tables/plots.
#' 3. If needed, rerun with residual PCA (`"overall"` or `"both"`).
#'
#' @return
#' An object of class `mfrm_diagnostics` including:
#' - `obs`: observed/expected/residual-level table
#' - `measures`: facet/person fit table (`Infit`, `Outfit`, `ZSTD`, `PTMEA`)
#' - `overall_fit`: overall fit summary
#' - `fit`: element-level fit diagnostics
#' - `reliability`: facet-level model/real separation and reliability
#' - `precision_profile`: one-row summary of the active precision tier and its
#'   recommended use
#' - `precision_audit`: package-native checks for SE, CI, and reliability
#' - `facet_precision`: facet-level precision summary by distribution basis and
#'   SE mode
#' - `facets_chisq`: fixed/random facet variability summary
#' - `interactions`: top interaction diagnostics
#' - `interrater`: inter-rater agreement bundle (`summary`, `pairs`) including
#'   agreement and rater-severity spread indices
#' - `unexpected`: unexpected-response bundle
#' - `fair_average`: adjusted-score reference bundle
#' - `displacement`: displacement diagnostics bundle
#' - `approximation_notes`: method notes for SE/CI/reliability summaries
#' - `residual_pca_overall`: optional overall PCA object
#' - `residual_pca_by_facet`: optional facet PCA objects
#'
#' @seealso [fit_mfrm()], [analyze_residual_pca()], [build_visual_summaries()],
#'   [mfrmr_visual_diagnostics], [mfrmr_reporting_and_apa]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' s_diag <- summary(diag)
#' s_diag$overview[, c("Observations", "Facets", "Categories")]
#' p_qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE)
#' class(p_qc)
#'
#' # Optional: include residual PCA in the diagnostic bundle
#' diag_pca <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag_pca, mode = "overall")
#' head(pca$overall_table)
#'
#' # Reporting route:
#' prec <- precision_audit_report(fit, diagnostics = diag)
#' summary(prec)
#' }
#' @export
diagnose_mfrm <- function(fit,
                          interaction_pairs = NULL,
                          top_n_interactions = 20,
                          whexact = FALSE,
                          residual_pca = c("none", "overall", "facet", "both"),
                          pca_max_factors = 10L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))

  out <- mfrm_diagnostics(
    fit,
    interaction_pairs = interaction_pairs,
    top_n_interactions = top_n_interactions,
    whexact = whexact,
    residual_pca = residual_pca,
    pca_max_factors = pca_max_factors
  )
  class(out) <- c("mfrm_diagnostics", class(out))
  out
}

#' Compare two or more fitted MFRM models
#'
#' Produce a side-by-side comparison of multiple [fit_mfrm()] results using
#' information criteria, log-likelihood, and parameter counts. When exactly
#' two nested models are supplied, a likelihood-ratio test is included.
#'
#' @param ... Two or more `mfrm_fit` objects to compare.
#' @param labels Optional character vector of labels for each model.
#'   If `NULL`, labels are generated from model/method combinations.
#' @param warn_constraints Logical. If `TRUE` (the default), emit a warning
#'   when models use different centering constraints (`noncenter_facet` or
#'   `dummy_facets`), which can make information-criterion comparisons
#'   misleading.
#' @param nested Logical. Set to `TRUE` only when the supplied models are
#'   known to be nested and fitted with the same likelihood basis on the same
#'   observations. The default is `FALSE`, in which case no likelihood-ratio
#'   test is reported. When `TRUE`, the function still runs a conservative
#'   structural audit and computes the LRT only for supported nesting patterns.
#'
#' @details
#' Models should be fit to the **same data** (same rows, same person/facet
#' columns) for the comparison to be meaningful. The function checks that
#' observation counts match and warns otherwise.
#'
#' Information-criterion ranking is reported only when all candidate models
#' use the package's `MML` estimation path, analyze the same observations, and
#' converge successfully. Raw `AIC` and `BIC` values are still shown for each
#' model, but `Delta_*`, weights, and preferred-model summaries are suppressed
#' when the likelihood basis is not comparable enough for primary reporting.
#'
#' **Nesting**: Two models are *nested* when one is a special case of the
#' other obtained by imposing equality constraints.  The most common
#' nesting in MFRM is RSM (shared thresholds) inside PCM
#' (item-specific thresholds).  Models that differ only in estimation
#' method (MML vs JML) on the same specification are not nested in the
#' usual sense---use information criteria rather than LRT for that
#' comparison.
#'
#' The **likelihood-ratio test (LRT)** is reported only when exactly two
#' models are supplied, `nested = TRUE`, the structural audit passes, and the
#' difference in the number of parameters is positive:
#'
#' \deqn{\Lambda = -2 (\ell_{\mathrm{restricted}} - \ell_{\mathrm{full}})
#'   \sim \chi^2_{\Delta p}}
#'
#' The LRT is asymptotically valid when models are nested and the data
#' are independent.  With small samples or boundary conditions (e.g.,
#' variance components near zero), treat p-values as approximate.
#'
#' @section Information-criterion diagnostics:
#' In addition to raw AIC and BIC values, the function computes:
#' - **Delta_AIC / Delta_BIC**: difference from the best (minimum) value.
#'   A Delta < 2 is typically considered negligible; 4--7 suggests
#'   moderate evidence; > 10 indicates strong evidence against the
#'   higher-scoring model (Burnham & Anderson, 2002).
#' - **AkaikeWeight / BICWeight**: model probabilities derived from
#'   `exp(-0.5 * Delta)`, normalised across the candidate set.  An
#'   Akaike weight of 0.90 means the model has a 90\% probability of
#'   being the best in the candidate set.
#' - **Evidence ratios**: pairwise ratios of Akaike weights, quantifying
#'   the relative evidence for one model over another (e.g., an
#'   evidence ratio of 5 means the preferred model is 5 times more
#'   likely).
#'
#' AIC penalises complexity less than BIC; when they disagree, AIC
#' favours the more complex model and BIC the simpler one.
#'
#' @section What this comparison means:
#' `compare_mfrm()` is a same-basis model-comparison helper. Its strongest
#' claims apply only when the models were fit to the same response data,
#' under a compatible likelihood basis, and with compatible constraint
#' structure.
#'
#' @section What this comparison does not justify:
#' - Do not treat AIC/BIC differences as primary evidence when
#'   `table$ICComparable` is `FALSE`.
#' - Do not interpret the LRT unless `nested = TRUE` and the structural audit
#'   in `comparison_basis$nesting_audit` passes.
#' - Do not compare models fit to different datasets, different score codings,
#'   or materially different constraint systems as if they were commensurate.
#'
#' @section Interpreting output:
#' - Lower AIC/BIC values indicate better parsimony-accuracy trade-off only
#'   when `table$ICComparable` is `TRUE`.
#' - A significant LRT p-value suggests the more complex model provides a
#'   meaningfully better fit only when the nesting assumption truly holds.
#' - `preferred` indicates the model preferred by each criterion.
#' - `evidence_ratios` gives pairwise Akaike-weight ratios (returned only
#'   when Akaike weights can be computed for at least two models).
#' - When comparing more than two models, interpret evidence ratios
#'   cautiously---they do not adjust for multiple comparisons.
#'
#' @section How to read the main outputs:
#' - `table`: first-pass comparison table; start with `ICComparable`,
#'   `Model`, `Method`, `AIC`, and `BIC`.
#' - `comparison_basis`: records whether IC and LRT claims are defensible for
#'   the supplied models.
#' - `lrt`: nested-model test summary, present only when the requested and
#'   audited conditions are met.
#' - `preferred`: candidate preferred by each criterion when those summaries
#'   are available.
#'
#' @section Recommended next step:
#' Inspect `comparison_basis` before writing conclusions. If comparability is
#' weak, treat the result as descriptive and revise the model setup (for
#' example, explicit `step_facet`, common data, or common constraints) before
#' using IC or LRT results in reporting.
#'
#' @section Typical workflow:
#' 1. Fit two models with [fit_mfrm()] (e.g., RSM and PCM).
#' 2. Compare with `compare_mfrm(fit_rsm, fit_pcm)`.
#' 3. Inspect `summary(comparison)` for AIC/BIC diagnostics and, when
#'    appropriate, an LRT.
#'
#' @return
#' An object of class `mfrm_comparison` (named list) with:
#' - `table`: data.frame of model-level statistics (LogLik, AIC, BIC,
#'   Delta_AIC, AkaikeWeight, Delta_BIC, BICWeight, npar, nobs, Model,
#'   Method, Converged, ICComparable).
#' - `lrt`: data.frame with likelihood-ratio test result (only when two models
#'   are supplied and `nested = TRUE`). Contains `ChiSq`, `df`, `p_value`.
#' - `evidence_ratios`: data.frame of pairwise Akaike-weight ratios (Model1,
#'   Model2, EvidenceRatio). `NULL` when weights cannot be computed.
#' - `preferred`: named list with the preferred model label by each criterion.
#' - `comparison_basis`: list describing whether IC and LRT comparisons were
#'   considered comparable. Includes a conservative `nesting_audit`.
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#'
#' fit_rsm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "MML", model = "RSM", maxit = 25)
#' fit_pcm <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "MML", model = "PCM",
#'                      step_facet = "Criterion", maxit = 25)
#' comp <- compare_mfrm(fit_rsm, fit_pcm, labels = c("RSM", "PCM"))
#' comp$table
#' comp$evidence_ratios
#' }
#' @export
compare_mfrm <- function(..., labels = NULL, warn_constraints = TRUE, nested = FALSE) {
  fits <- list(...)
  if (length(fits) < 2) {
    stop("`compare_mfrm()` requires at least two `mfrm_fit` objects.",
         call. = FALSE)
  }
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "mfrm_fit")) {
      stop("Argument ", i, " is not an `mfrm_fit` object. Got: ",
           class(fits[[i]])[1], ".", call. = FALSE)
    }
  }

  # -- build labels --
  if (is.null(labels)) {
    labels <- vapply(fits, function(f) {
      m <- if (!is.null(f$config$model)) toupper(f$config$model[1]) else "?"
      e <- if (!is.null(f$config$method)) toupper(f$config$method[1]) else "?"
      paste0(m, "/", e)
    }, character(1))
    if (anyDuplicated(labels)) {
      labels <- paste0(labels, " (", seq_along(labels), ")")
    }
  }
  if (length(labels) != length(fits)) {
    stop("`labels` must have the same length as the number of models (",
         length(fits), ").", call. = FALSE)
  }

  normalize_compare_col <- function(col) {
    if (is.factor(col)) return(as.character(col))
    if (inherits(col, "POSIXt")) return(format(col, tz = "UTC", usetz = TRUE))
    if (is.numeric(col)) return(format(col, digits = 17, trim = TRUE, scientific = FALSE))
    if (is.logical(col)) return(ifelse(is.na(col), "NA", ifelse(col, "TRUE", "FALSE")))
    if (is.character(col)) return(col)
    as.character(col)
  }

  canonicalize_compare_data <- function(dat) {
    if (is.null(dat) || !is.data.frame(dat)) return(NULL)
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
    dat <- dat[, sort(names(dat)), drop = FALSE]
    dat[] <- lapply(dat, normalize_compare_col)
    if (nrow(dat) > 0 && ncol(dat) > 0) {
      ord <- do.call(order, c(dat, list(na.last = TRUE, method = "radix")))
      dat <- dat[ord, , drop = FALSE]
    }
    rownames(dat) <- NULL
    dat
  }

  compare_data_equal <- function(fit_a, fit_b) {
    dat_a <- canonicalize_compare_data(fit_a$prep$data %||% NULL)
    dat_b <- canonicalize_compare_data(fit_b$prep$data %||% NULL)
    if (is.null(dat_a) || is.null(dat_b)) return(FALSE)
    identical(dat_a, dat_b)
  }

  # -- constraint compatibility check --
  if (isTRUE(warn_constraints) && length(fits) >= 2) {
    noncenter <- vapply(fits, function(f) {
      nc <- f$config$noncenter_facet
      if (is.null(nc)) "Person" else as.character(nc[1])
    }, character(1))
    if (length(unique(noncenter)) > 1) {
      warning("Models use different centering constraints (",
              paste(unique(noncenter), collapse = ", "),
              "). IC comparisons may be misleading.", call. = FALSE)
    }
    dummy_sets <- lapply(fits, function(f) {
      d <- f$config$dummy_facets
      if (is.null(d)) character(0) else sort(as.character(d))
    })
    dummy_sigs <- vapply(dummy_sets, paste, character(1), collapse = ",")
    if (length(unique(dummy_sigs)) > 1) {
      warning("Models use different dummy-facet constraints (",
              paste(unique(dummy_sigs), collapse = " vs "),
              "). IC comparisons may be misleading.", call. = FALSE)
    }
  }

  # -- extract summary statistics --
  rows <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    s <- f$summary
    nobs <- if (!is.null(f$prep$data)) nrow(f$prep$data) else NA_integer_
    npar <- if (!is.null(f$opt$par)) length(f$opt$par) else NA_integer_
    method <- if (!is.null(f$config$method)) toupper(f$config$method[1]) else NA_character_
    method <- ifelse(identical(method, "JMLE"), "JML", method)
    tibble(
      Label     = labels[i],
      Model     = if (!is.null(f$config$model)) toupper(f$config$model[1]) else NA_character_,
      Method    = method,
      nobs      = nobs,
      npar      = npar,
      LogLik    = if ("LogLik" %in% names(s)) s$LogLik[1] else NA_real_,
      AIC       = if ("AIC" %in% names(s)) s$AIC[1] else NA_real_,
      BIC       = if ("BIC" %in% names(s)) s$BIC[1] else NA_real_,
      Converged = if ("Converged" %in% names(s)) s$Converged[1] else NA
    )
  })
  tbl <- bind_rows(rows)

  method_vals <- tbl$Method[!is.na(tbl$Method)]
  same_method <- length(unique(method_vals)) <= 1
  all_mml <- length(method_vals) > 0 && all(method_vals == "MML")
  obs_vals <- tbl$nobs[is.finite(tbl$nobs)]
  same_nobs <- length(unique(obs_vals)) <= 1
  same_data <- if (length(fits) >= 2) {
    all(vapply(fits[-1], function(f) compare_data_equal(fits[[1]], f), logical(1)))
  } else {
    TRUE
  }
  conv_vals <- tbl$Converged
  all_converged <- length(conv_vals) > 0 && all(!is.na(conv_vals) & as.logical(conv_vals))
  ic_comparable <- same_method && same_nobs && same_data && all_converged && all_mml

  if (!isTRUE(nested)) {
    if (!same_method) {
      warning(
        "Models use different estimation methods (",
        paste(unique(method_vals), collapse = ", "),
        "). Raw AIC/BIC values are shown, but cross-method deltas, weights, ",
        "and automatic preferences are suppressed.",
        call. = FALSE
      )
    } else if (!all_mml && all_converged) {
      warning(
        "Information-criterion ranking is limited to converged MML fits in this package. ",
        "Raw AIC/BIC values are shown for ",
        paste(unique(method_vals), collapse = ", "),
        " models, but deltas, weights, automatic preferences, and LRT were suppressed.",
        call. = FALSE
      )
    }

    # -- warn if observation counts differ --
    if (!same_nobs) {
      warning("Models were fit to different numbers of observations (",
              paste(obs_vals, collapse = ", "),
              "). Raw AIC/BIC values are shown, but cross-sample deltas, weights, ",
              "and automatic preferences are suppressed.", call. = FALSE)
    }
    if (!same_data) {
      warning(
        "Models were not fit to the same prepared response data. Raw AIC/BIC values are shown, ",
        "but deltas, weights, automatic preferences, and likelihood-ratio testing were suppressed.",
        call. = FALSE
      )
    }
    if (!all_converged) {
      warning(
        "At least one compared model did not converge. Raw AIC/BIC values are shown, ",
        "but IC ranking, weights, and likelihood-ratio testing were suppressed.",
        call. = FALSE
      )
    }
  }

  tbl$ICComparable <- ic_comparable
  preferred <- list()

  # -- Delta AIC and Akaike Weights --
  if (ic_comparable && any(is.finite(tbl$AIC))) {
    min_aic <- min(tbl$AIC, na.rm = TRUE)
    tbl$Delta_AIC <- tbl$AIC - min_aic
    raw_w <- exp(-0.5 * tbl$Delta_AIC)
    tbl$AkaikeWeight <- ifelse(is.finite(raw_w),
                               raw_w / sum(raw_w, na.rm = TRUE),
                               NA_real_)
    preferred$AIC <- tbl$Label[which.min(tbl$AIC)]
  } else {
    tbl$Delta_AIC <- NA_real_
    tbl$AkaikeWeight <- NA_real_
  }

  # -- Delta BIC and BIC Weights --
  if (ic_comparable && any(is.finite(tbl$BIC))) {
    min_bic <- min(tbl$BIC, na.rm = TRUE)
    tbl$Delta_BIC <- tbl$BIC - min_bic
    raw_bw <- exp(-0.5 * tbl$Delta_BIC)
    tbl$BICWeight <- ifelse(is.finite(raw_bw),
                            raw_bw / sum(raw_bw, na.rm = TRUE),
                            NA_real_)
    preferred$BIC <- tbl$Label[which.min(tbl$BIC)]
  } else {
    tbl$Delta_BIC <- NA_real_
    tbl$BICWeight <- NA_real_
  }

  # -- likelihood-ratio test (two models only) --
  lrt <- NULL
  nesting_audit <- audit_compare_mfrm_nesting(fits, labels = labels)
  if (isTRUE(nested) && length(fits) == 2) {
    if (!same_method || !same_nobs || !same_data || !all_converged || !all_mml) {
      warning(
        "`nested = TRUE` was requested, but the models do not share the same ",
        "formal MML likelihood basis, observation set, and convergence status. LRT was not computed.",
        call. = FALSE
      )
    } else if (!isTRUE(nesting_audit$eligible)) {
      warning(
        "`nested = TRUE` was requested, but the structural nesting audit did not pass. ",
        nesting_audit$reason,
        " LRT was not computed.",
        call. = FALSE
      )
    } else {
      ll <- tbl$LogLik
      np <- tbl$npar
      if (all(is.finite(ll)) && all(is.finite(np)) && np[1] != np[2]) {
        idx_simple <- which.min(np)
        idx_complex <- which.max(np)
        chi_sq <- 2 * (ll[idx_complex] - ll[idx_simple])
        df_diff <- np[idx_complex] - np[idx_simple]
        if (chi_sq >= 0 && df_diff > 0) {
          p_val <- stats::pchisq(chi_sq, df = df_diff, lower.tail = FALSE)
          lrt <- tibble(
            Simple   = tbl$Label[idx_simple],
            Complex  = tbl$Label[idx_complex],
            ChiSq    = chi_sq,
            df       = df_diff,
            p_value  = p_val
          )
          preferred$LRT <- if (p_val < 0.05) tbl$Label[idx_complex] else tbl$Label[idx_simple]
        }
      }
    }
  }

  # -- evidence ratios (pairwise Akaike-weight ratios) --
  evidence_ratios <- NULL
  if (ic_comparable && "AkaikeWeight" %in% names(tbl) && nrow(tbl) >= 2) {
    er_rows <- list()
    for (i in 1:(nrow(tbl) - 1)) {
      for (j in (i + 1):nrow(tbl)) {
        w_i <- tbl$AkaikeWeight[i]
        w_j <- tbl$AkaikeWeight[j]
        er <- if (is.finite(w_i) && is.finite(w_j) && w_j > 0) {
          w_i / w_j
        } else {
          NA_real_
        }
        er_rows[[length(er_rows) + 1]] <- tibble(
          Model1        = tbl$Label[i],
          Model2        = tbl$Label[j],
          EvidenceRatio = er
        )
      }
    }
    evidence_ratios <- bind_rows(er_rows)
  }

  comparison_basis <- list(
    same_method = same_method,
    all_mml = all_mml,
    same_nobs = same_nobs,
    same_data = same_data,
    all_converged = all_converged,
    ic_comparable = ic_comparable,
    nested_requested = isTRUE(nested),
    nesting_audit = nesting_audit
  )

  out <- list(
    table           = tbl,
    lrt             = lrt,
    evidence_ratios = evidence_ratios,
    preferred       = preferred,
    comparison_basis = comparison_basis
  )
  class(out) <- c("mfrm_comparison", class(out))
  out
}

#' @export
summary.mfrm_comparison <- function(object, ...) {
  out <- list(
    table           = object$table,
    lrt             = object$lrt,
    evidence_ratios = object$evidence_ratios,
    preferred       = object$preferred,
    comparison_basis = object$comparison_basis
  )
  class(out) <- "summary.mfrm_comparison"
  out
}

#' @export
print.summary.mfrm_comparison <- function(x, ...) {
  cat("--- MFRM Model Comparison ---\n\n")

  # -- main comparison table --
  tbl <- x$table
  # Format weight columns with 4 decimal places for readability
  fmt_tbl <- as.data.frame(tbl)
  weight_cols <- intersect(c("AkaikeWeight", "BICWeight"), names(fmt_tbl))
  for (wc in weight_cols) {
    fmt_tbl[[wc]] <- ifelse(is.na(fmt_tbl[[wc]]), NA_character_,
                            sprintf("%.4f", fmt_tbl[[wc]]))
  }
  delta_cols <- intersect(c("Delta_AIC", "Delta_BIC"), names(fmt_tbl))
  for (dc in delta_cols) {
    fmt_tbl[[dc]] <- ifelse(is.na(fmt_tbl[[dc]]), NA_character_,
                            sprintf("%.2f", fmt_tbl[[dc]]))
  }
  print(fmt_tbl, row.names = FALSE)

  if (!isTRUE(x$comparison_basis$ic_comparable)) {
    cat("\nInformation-criterion ranking was suppressed because the models do not share\n")
    cat("a comparable formal MML likelihood basis, observation set, and convergence status.\n")
  }

  if (!is.null(x$lrt)) {
    cat("\nLikelihood-ratio test:\n")
    cat(sprintf("  Chi-sq = %.3f, df = %d, p = %.4f\n",
                x$lrt$ChiSq[1], x$lrt$df[1], x$lrt$p_value[1]))
    cat(sprintf("  %s vs %s\n", x$lrt$Simple[1], x$lrt$Complex[1]))
  } else if (isTRUE(x$comparison_basis$nested_requested)) {
    cat("\nLikelihood-ratio test was not reported.\n")
    audit <- x$comparison_basis$nesting_audit %||% list()
    if (!is.null(audit$reason) && nzchar(audit$reason)) {
      cat("  Nesting audit:", audit$reason, "\n")
    }
  }

  # -- evidence ratios --
  if (!is.null(x$evidence_ratios) && nrow(x$evidence_ratios) > 0) {
    cat("\nEvidence ratios (Akaike weights):\n")
    er <- x$evidence_ratios
    for (k in seq_len(nrow(er))) {
      ratio_str <- if (is.finite(er$EvidenceRatio[k])) {
        sprintf("%.2f", er$EvidenceRatio[k])
      } else {
        "NA"
      }
      cat(sprintf("  %s / %s = %s\n",
                  er$Model1[k], er$Model2[k], ratio_str))
    }
  }

  if (length(x$preferred) > 0) {
    cat("\nPreferred model:\n")
    for (nm in names(x$preferred)) {
      cat(sprintf("  By %s: %s\n", nm, x$preferred[[nm]]))
    }
  }
  invisible(x)
}

#' @export
print.mfrm_comparison <- function(x, ...) {
  print(summary(x))
  invisible(x)
}
