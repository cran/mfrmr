#' Build an explicit simulation specification for MFRM design studies
#'
#' @param n_person Number of persons/respondents to generate.
#' @param n_rater Number of rater facet levels to generate.
#' @param n_criterion Number of criterion/item facet levels to generate.
#' @param raters_per_person Number of raters assigned to each person.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread used to generate equally spaced thresholds when
#'   `thresholds = NULL`.
#' @param thresholds Optional threshold specification. Use either a numeric
#'   vector of common thresholds or a data frame with columns `StepFacet`,
#'   `Step`/`StepIndex`, and `Estimate`.
#' @param model Measurement model recorded in the simulation specification.
#' @param step_facet Step facet used when `model = "PCM"` and threshold values
#'   vary across levels.
#' @param assignment Assignment design. `"crossed"` means every person sees
#'   every rater; `"rotating"` uses a balanced rotating subset; `"resampled"`
#'   reuses empirical person-level rater-assignment profiles; `"skeleton"`
#'   reuses an observed person-by-facet design skeleton.
#' @param latent_distribution Latent-value generator. `"normal"` samples from
#'   centered normal distributions using the supplied standard deviations.
#'   `"empirical"` resamples centered support values from
#'   `empirical_person`/`empirical_rater`/`empirical_criterion`.
#' @param empirical_person Optional numeric support values used when
#'   `latent_distribution = "empirical"`.
#' @param empirical_rater Optional numeric support values used when
#'   `latent_distribution = "empirical"`.
#' @param empirical_criterion Optional numeric support values used when
#'   `latent_distribution = "empirical"`.
#' @param assignment_profiles Optional data frame with columns
#'   `TemplatePerson` and `Rater` (optionally `Group`) describing empirical
#'   person-level rater-assignment profiles used when `assignment = "resampled"`.
#' @param design_skeleton Optional data frame with columns `TemplatePerson`,
#'   `Rater`, and `Criterion` (optionally `Group` and `Weight`) describing an
#'   observed response skeleton used when `assignment = "skeleton"`.
#' @param group_levels Optional character vector of group labels.
#' @param dif_effects Optional data frame of true group-linked DIF effects.
#' @param interaction_effects Optional data frame of true interaction effects.
#'
#' @details
#' `build_mfrm_sim_spec()` creates an explicit, portable simulation
#' specification that can be passed to [simulate_mfrm_data()]. The goal is to
#' make the data-generating mechanism inspectable and reusable rather than
#' relying only on ad hoc scalar arguments.
#'
#' The resulting object records:
#' - design counts (`n_person`, `n_rater`, `n_criterion`, `raters_per_person`)
#' - latent spread assumptions (`theta_sd`, `rater_sd`, `criterion_sd`)
#' - optional empirical latent support values for semi-parametric simulation
#' - threshold structure (`threshold_table`)
#' - assignment design (`assignment`)
#' - optional empirical assignment profiles (`assignment_profiles`) with
#'   optional person-level `Group` labels
#' - optional observed response skeleton (`design_skeleton`)
#'   with optional person-level `Group` labels and observation-level `Weight`
#'   values
#' - optional signal tables for DIF and interaction bias
#'
#' The current generator still targets the package's standard person x rater x
#' criterion workflow. When threshold values are provided by `StepFacet`, the
#' supported step facets are the generated `Rater` or `Criterion` levels.
#'
#' @section Interpreting output:
#' This object does not contain simulated data. It is a data-generating
#' specification that tells [simulate_mfrm_data()] how to generate them.
#'
#' @return An object of class `mfrm_sim_spec`.
#' @seealso [extract_mfrm_sim_spec()], [simulate_mfrm_data()]
#' @examples
#' spec <- build_mfrm_sim_spec(
#'   n_person = 30,
#'   n_rater = 4,
#'   n_criterion = 3,
#'   raters_per_person = 2,
#'   assignment = "rotating"
#' )
#' spec$model
#' spec$assignment
#' spec$threshold_table
#' @export
build_mfrm_sim_spec <- function(n_person = 50,
                                n_rater = 4,
                                n_criterion = 4,
                                raters_per_person = n_rater,
                                score_levels = 4,
                                theta_sd = 1,
                                rater_sd = 0.35,
                                criterion_sd = 0.25,
                                noise_sd = 0,
                                step_span = 1.4,
                                thresholds = NULL,
                                model = c("RSM", "PCM"),
                                step_facet = "Criterion",
                                assignment = c("crossed", "rotating", "resampled", "skeleton"),
                                latent_distribution = c("normal", "empirical"),
                                empirical_person = NULL,
                                empirical_rater = NULL,
                                empirical_criterion = NULL,
                                assignment_profiles = NULL,
                                design_skeleton = NULL,
                                group_levels = NULL,
                                dif_effects = NULL,
                                interaction_effects = NULL) {
  model <- match.arg(toupper(as.character(model[1])), c("RSM", "PCM"))
  assignment <- match.arg(tolower(as.character(assignment[1])), c("crossed", "rotating", "resampled", "skeleton"))
  latent_distribution <- match.arg(tolower(as.character(latent_distribution[1])), c("normal", "empirical"))

  n_person <- simulation_validate_count(n_person, "n_person", min_value = 2L)
  n_rater <- simulation_validate_count(n_rater, "n_rater", min_value = 2L)
  n_criterion <- simulation_validate_count(n_criterion, "n_criterion", min_value = 2L)
  raters_per_person <- simulation_validate_count(raters_per_person, "raters_per_person", min_value = 1L)
  score_levels <- simulation_validate_count(score_levels, "score_levels", min_value = 2L)

  if (raters_per_person > n_rater) {
    stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  }
  if (assignment == "crossed" && raters_per_person != n_rater) {
    stop("`assignment = \"crossed\"` requires `raters_per_person == n_rater`.", call. = FALSE)
  }

  step_facet <- as.character(step_facet[1] %||% "Criterion")
  if (!step_facet %in% c("Criterion", "Rater")) {
    stop("`step_facet` must currently be either \"Criterion\" or \"Rater\" for simulation.", call. = FALSE)
  }

  theta_sd <- simulation_validate_numeric(theta_sd, "theta_sd", lower = 0)
  rater_sd <- simulation_validate_numeric(rater_sd, "rater_sd", lower = 0)
  criterion_sd <- simulation_validate_numeric(criterion_sd, "criterion_sd", lower = 0)
  noise_sd <- simulation_validate_numeric(noise_sd, "noise_sd", lower = 0)
  step_span <- simulation_validate_numeric(step_span, "step_span", lower = 0)

  if (!is.null(group_levels)) {
    group_levels <- as.character(group_levels)
    group_levels <- unique(group_levels[!is.na(group_levels) & nzchar(group_levels)])
    if (length(group_levels) < 1L) {
      stop("`group_levels` must contain at least one non-empty label.", call. = FALSE)
    }
  }

  threshold_table <- simulation_build_threshold_table(
    thresholds = thresholds,
    score_levels = score_levels,
    step_span = step_span,
    model = model
  )

  empirical_support <- simulation_build_empirical_support(
    latent_distribution = latent_distribution,
    empirical_person = empirical_person,
    empirical_rater = empirical_rater,
    empirical_criterion = empirical_criterion
  )
  assignment_profiles <- simulation_normalize_assignment_profiles(
    assignment_profiles = assignment_profiles,
    assignment = assignment,
    n_rater = n_rater
  )
  design_skeleton <- simulation_normalize_design_skeleton(
    design_skeleton = design_skeleton,
    assignment = assignment,
    n_rater = n_rater,
    n_criterion = n_criterion
  )

  dif_effects <- simulation_normalize_effects(
    effects = dif_effects,
    arg_name = "dif_effects",
    allowed_cols = c("Group", "Person", "Rater", "Criterion")
  )
  interaction_effects <- simulation_normalize_effects(
    effects = interaction_effects,
    arg_name = "interaction_effects",
    allowed_cols = c("Group", "Person", "Rater", "Criterion")
  )

  structure(
    list(
      n_person = n_person,
      n_rater = n_rater,
      n_criterion = n_criterion,
      raters_per_person = raters_per_person,
      score_levels = score_levels,
      theta_sd = theta_sd,
      rater_sd = rater_sd,
      criterion_sd = criterion_sd,
      noise_sd = noise_sd,
      step_span = step_span,
      model = model,
      step_facet = step_facet,
      assignment = assignment,
      latent_distribution = latent_distribution,
      empirical_support = empirical_support,
      assignment_profiles = assignment_profiles,
      design_skeleton = design_skeleton,
      threshold_table = threshold_table,
      group_levels = group_levels,
      dif_effects = dif_effects,
      interaction_effects = interaction_effects,
      source = "manual"
    ),
    class = "mfrm_sim_spec"
  )
}

#' Derive a simulation specification from a fitted MFRM object
#'
#' @param fit Output from [fit_mfrm()].
#' @param assignment Assignment design to record in the returned specification.
#'   Use `"resampled"` to reuse empirical person-level rater-assignment
#'   profiles from the fitted data, or `"skeleton"` to reuse the observed
#'   person-by-facet design skeleton from the fitted data.
#' @param latent_distribution Latent-value generator to record in the returned
#'   specification. `"normal"` stores spread summaries for parametric draws;
#'   `"empirical"` additionally activates centered empirical resampling from the
#'   fitted person/rater/criterion estimates.
#' @param source_data Optional original source data used to recover additional
#'   non-calibration columns, currently person-level `group` labels, when
#'   building a fit-derived observed response skeleton.
#' @param person Optional person column name in `source_data`. Defaults to the
#'   person column recorded in `fit`.
#' @param group Optional group column name in `source_data` to merge into the
#'   returned `design_skeleton` as person-level metadata.
#'
#' @details
#' `extract_mfrm_sim_spec()` uses a fitted model as a practical starting point
#' for later simulation studies. It extracts:
#' - design counts from the fitted data
#' - empirical spread of person and facet estimates
#' - optional empirical support values for semi-parametric draws
#' - fitted threshold values
#' - either a simplified assignment summary (`"crossed"` / `"rotating"`),
#'   empirical resampled assignment profiles (`"resampled"`), or an observed
#'   response skeleton (`"skeleton"`, optionally carrying `Group`/`Weight`)
#'
#' This is intended as a **fit-derived parametric starting point**, not as a
#' claim that the fitted object perfectly recovers the true data-generating
#' mechanism. Users should review and, if necessary, edit the returned
#' specification before using it for design planning.
#'
#' If you want to carry person-level group labels into a fit-derived observed
#' response skeleton, provide the original `source_data` together with
#' `person` and `group`. Group labels are treated as person-level metadata and
#' are checked for one-label-per-person consistency before being merged.
#'
#' @section Interpreting output:
#' The returned object is a simulation specification, not a prediction about one
#' future sample. It captures one convenient approximation to the observed
#' design and estimated spread in the fitted run.
#'
#' @return An object of class `mfrm_sim_spec`.
#' @seealso [build_mfrm_sim_spec()], [simulate_mfrm_data()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 15)
#' spec <- extract_mfrm_sim_spec(fit, latent_distribution = "empirical")
#' spec$assignment
#' spec$model
#' head(spec$threshold_table)
#' @export
extract_mfrm_sim_spec <- function(fit,
                                  assignment = c("auto", "crossed", "rotating", "resampled", "skeleton"),
                                  latent_distribution = c("normal", "empirical"),
                                  source_data = NULL,
                                  person = NULL,
                                  group = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be output from fit_mfrm().", call. = FALSE)
  }

  assignment <- match.arg(tolower(as.character(assignment[1])), c("auto", "crossed", "rotating", "resampled", "skeleton"))
  latent_distribution <- match.arg(tolower(as.character(latent_distribution[1])), c("normal", "empirical"))
  prep <- fit$prep %||% NULL
  if (is.null(prep) || is.null(prep$data) || !is.data.frame(prep$data)) {
    stop("`fit` does not contain the prepared data needed to derive a simulation specification.", call. = FALSE)
  }

  required_facets <- c("Rater", "Criterion")
  if (!all(required_facets %in% names(prep$levels))) {
    stop("`extract_mfrm_sim_spec()` currently supports fitted models with `Rater` and `Criterion` facets.", call. = FALSE)
  }

  person_tbl <- fit$facets$person %||% tibble::tibble()
  other_tbl <- fit$facets$others %||% tibble::tibble()
  step_tbl <- fit$steps %||% tibble::tibble()

  rater_levels <- as.character(prep$levels$Rater)
  criterion_levels <- as.character(prep$levels$Criterion)
  person_levels <- as.character(prep$levels$Person)

  assignment_counts <- prep$data |>
    dplyr::distinct(.data$Person, .data$Rater) |>
    dplyr::count(.data$Person, name = "RatersPerPerson")
  raters_per_person <- as.integer(round(stats::median(assignment_counts$RatersPerPerson)))
  assignment_profiles <- prep$data |>
    dplyr::distinct(.data$Person, .data$Rater) |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$Person),
      Rater = as.character(.data$Rater)
    ) |>
    dplyr::arrange(.data$TemplatePerson, .data$Rater)
  keep_weight <- "Weight" %in% names(prep$data) && any(is.finite(prep$data$Weight) & prep$data$Weight != 1, na.rm = TRUE)
  design_keep <- c("Person", "Rater", "Criterion", intersect("Group", names(prep$data)))
  if (keep_weight) design_keep <- c(design_keep, "Weight")
  design_skeleton <- prep$data |>
    dplyr::select(dplyr::all_of(design_keep)) |>
    dplyr::distinct() |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$Person),
      Rater = as.character(.data$Rater),
      Criterion = as.character(.data$Criterion),
      Group = if ("Group" %in% names(prep$data)) as.character(.data$Group) else NA_character_,
      Weight = if (keep_weight) suppressWarnings(as.numeric(.data$Weight)) else NA_real_
    )
  if (!"Group" %in% names(prep$data)) {
    design_skeleton <- dplyr::select(design_skeleton, -dplyr::all_of("Group"))
  }
  if (!keep_weight) {
    design_skeleton <- dplyr::select(design_skeleton, -dplyr::all_of("Weight"))
  }
  person_col <- as.character(person[1] %||% prep$source_columns$person %||% "Person")
  group_col <- if (is.null(group)) NULL else as.character(group[1])
  if (!is.null(source_data) && !is.null(group_col)) {
    group_map <- simulation_extract_group_map(
      source_data = source_data,
      person_col = person_col,
      group_col = group_col,
      target_people = person_levels
    )
    assignment_profiles <- assignment_profiles |>
      dplyr::left_join(group_map, by = "TemplatePerson")
    if ("Group" %in% names(design_skeleton)) {
      design_skeleton <- design_skeleton |>
        dplyr::left_join(group_map, by = "TemplatePerson", suffix = c(".x", ".y")) |>
        dplyr::mutate(Group = dplyr::coalesce(.data$Group.y, .data$Group.x)) |>
        dplyr::select(-dplyr::any_of(c("Group.x", "Group.y")))
    } else {
      design_skeleton <- design_skeleton |>
        dplyr::left_join(group_map, by = "TemplatePerson")
    }
  }
  inferred_assignment <- if (all(assignment_counts$RatersPerPerson == length(rater_levels))) {
    "crossed"
  } else {
    "rotating"
  }
  if (assignment == "auto") assignment <- inferred_assignment

  score_values <- sort(unique(prep$data$Score))
  score_levels <- length(score_values)

  thresholds <- simulation_extract_thresholds_from_fit(
    step_tbl = step_tbl,
    model = as.character(fit$summary$Model[1] %||% fit$config$model %||% "RSM")
  )

  spec <- build_mfrm_sim_spec(
    n_person = length(person_levels),
    n_rater = length(rater_levels),
    n_criterion = length(criterion_levels),
    raters_per_person = raters_per_person,
    score_levels = score_levels,
    theta_sd = stats::sd(suppressWarnings(as.numeric(person_tbl$Estimate)), na.rm = TRUE),
    rater_sd = stats::sd(suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == "Rater"])), na.rm = TRUE),
    criterion_sd = stats::sd(suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == "Criterion"])), na.rm = TRUE),
    noise_sd = 0,
    step_span = if (is.null(thresholds) || nrow(thresholds) == 0) 0 else diff(range(thresholds$Estimate, na.rm = TRUE)) / 2,
    thresholds = thresholds,
    model = as.character(fit$summary$Model[1] %||% fit$config$model %||% "RSM"),
    step_facet = as.character(fit$config$step_facet %||% "Criterion"),
    assignment = assignment,
    latent_distribution = latent_distribution,
    empirical_person = suppressWarnings(as.numeric(person_tbl$Estimate)),
    empirical_rater = suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == "Rater"])),
    empirical_criterion = suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == "Criterion"])),
    assignment_profiles = assignment_profiles,
    design_skeleton = design_skeleton,
    group_levels = if ("Group" %in% names(design_skeleton)) sort(unique(as.character(design_skeleton$Group))) else NULL
  )

  spec$source <- "fit_mfrm"
  spec$source_summary <- list(
    observed_raters_per_person = assignment_counts,
    inferred_assignment = inferred_assignment,
    observed_score_values = score_values
  )
  spec
}

simulation_validate_count <- function(x, arg_name, min_value = 1L) {
  value <- as.integer(x[1])
  if (!is.finite(value) || value < min_value) {
    stop("`", arg_name, "` must be >= ", min_value, ".", call. = FALSE)
  }
  value
}

simulation_extract_group_map <- function(source_data,
                                         person_col,
                                         group_col,
                                         target_people) {
  if (!is.data.frame(source_data)) {
    stop("`source_data` must be a data.frame when `group` is supplied.", call. = FALSE)
  }
  required <- c(person_col, group_col)
  missing_cols <- setdiff(required, names(source_data))
  if (length(missing_cols) > 0) {
    stop("`source_data` is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }

  group_map <- tibble::as_tibble(source_data) |>
    dplyr::transmute(
      TemplatePerson = as.character(.data[[person_col]]),
      Group = as.character(.data[[group_col]])
    ) |>
    dplyr::filter(!is.na(.data$TemplatePerson), nzchar(.data$TemplatePerson),
                  !is.na(.data$Group), nzchar(.data$Group)) |>
    dplyr::distinct()

  if (nrow(group_map) == 0) {
    stop("`source_data` did not contain any valid person/group rows.", call. = FALSE)
  }

  ambiguity <- group_map |>
    dplyr::count(.data$TemplatePerson, name = "n_groups") |>
    dplyr::filter(.data$n_groups > 1L)
  if (nrow(ambiguity) > 0) {
    stop("`source_data` must assign at most one `group` label per person.", call. = FALSE)
  }

  missing_people <- setdiff(as.character(target_people), unique(group_map$TemplatePerson))
  if (length(missing_people) > 0) {
    stop(
      "`source_data` is missing `group` labels for fitted persons: ",
      paste(utils::head(missing_people, 5), collapse = ", "),
      if (length(missing_people) > 5) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  group_map
}

simulation_validate_numeric <- function(x, arg_name, lower = -Inf) {
  value <- as.numeric(x[1])
  if (!is.finite(value) || value < lower) {
    stop("`", arg_name, "` must be a finite numeric value >= ", lower, ".", call. = FALSE)
  }
  value
}

simulation_validate_empirical_vector <- function(x, arg_name) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    stop("`", arg_name, "` must contain at least two finite numeric values when `latent_distribution = \"empirical\"`.",
         call. = FALSE)
  }
  x
}

simulation_build_empirical_support <- function(latent_distribution,
                                               empirical_person,
                                               empirical_rater,
                                               empirical_criterion) {
  if (!identical(latent_distribution, "empirical")) {
    return(NULL)
  }

  list(
    person = simulation_validate_empirical_vector(empirical_person, "empirical_person"),
    rater = simulation_validate_empirical_vector(empirical_rater, "empirical_rater"),
    criterion = simulation_validate_empirical_vector(empirical_criterion, "empirical_criterion")
  )
}

simulation_normalize_assignment_profiles <- function(assignment_profiles,
                                                     assignment,
                                                     n_rater) {
  if (!identical(assignment, "resampled")) {
    return(NULL)
  }
  if (is.null(assignment_profiles)) {
    stop("`assignment = \"resampled\"` requires `assignment_profiles`.", call. = FALSE)
  }

  if (!is.data.frame(assignment_profiles)) {
    stop("`assignment_profiles` must be a data.frame with `TemplatePerson` and `Rater` columns.", call. = FALSE)
  }

  tbl <- tibble::as_tibble(assignment_profiles)
  if (!all(c("TemplatePerson", "Rater") %in% names(tbl))) {
    stop("`assignment_profiles` must include `TemplatePerson` and `Rater` columns.", call. = FALSE)
  }

  tbl <- tbl |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$TemplatePerson),
      Rater = as.character(.data$Rater),
      Group = if ("Group" %in% names(tbl)) as.character(.data$Group) else NA_character_
    ) |>
    dplyr::filter(!is.na(.data$TemplatePerson), nzchar(.data$TemplatePerson),
                  !is.na(.data$Rater), nzchar(.data$Rater)) |>
    dplyr::distinct()

  if (nrow(tbl) == 0) {
    stop("`assignment_profiles` did not contain any valid TemplatePerson/Rater rows.", call. = FALSE)
  }
  if ("Group" %in% names(tbl)) {
    group_check <- tbl |>
      dplyr::filter(!is.na(.data$Group), nzchar(.data$Group)) |>
      dplyr::distinct(.data$TemplatePerson, .data$Group) |>
      dplyr::count(.data$TemplatePerson, name = "n_groups")
    if (nrow(group_check) > 0 && any(group_check$n_groups > 1L)) {
      stop("`assignment_profiles` must assign at most one `Group` label per `TemplatePerson`.", call. = FALSE)
    }
    if (all(is.na(tbl$Group) | !nzchar(tbl$Group))) {
      tbl <- dplyr::select(tbl, -dplyr::all_of("Group"))
    }
  }
  observed_raters <- unique(tbl$Rater)
  if (length(observed_raters) != n_rater) {
    stop("`assignment_profiles` must reference exactly ", n_rater, " distinct rater levels.", call. = FALSE)
  }
  tbl
}

simulation_normalize_design_skeleton <- function(design_skeleton,
                                                 assignment,
                                                 n_rater,
                                                 n_criterion) {
  if (!identical(assignment, "skeleton")) {
    return(NULL)
  }
  if (is.null(design_skeleton)) {
    stop("`assignment = \"skeleton\"` requires `design_skeleton`.", call. = FALSE)
  }

  if (!is.data.frame(design_skeleton)) {
    stop("`design_skeleton` must be a data.frame.", call. = FALSE)
  }

  tbl <- tibble::as_tibble(design_skeleton)
  required <- c("TemplatePerson", "Rater", "Criterion")
  if (!all(required %in% names(tbl))) {
    stop("`design_skeleton` must include `TemplatePerson`, `Rater`, and `Criterion` columns.", call. = FALSE)
  }

  keep_cols <- c(required, intersect(c("Group", "Weight"), names(tbl)))
  tbl <- tbl |>
    dplyr::transmute(!!!rlang::syms(keep_cols)) |>
    dplyr::mutate(
      TemplatePerson = as.character(.data$TemplatePerson),
      Rater = as.character(.data$Rater),
      Criterion = as.character(.data$Criterion)
    ) |>
    dplyr::filter(!is.na(.data$TemplatePerson), nzchar(.data$TemplatePerson),
                  !is.na(.data$Rater), nzchar(.data$Rater),
                  !is.na(.data$Criterion), nzchar(.data$Criterion)) |>
    dplyr::distinct()

  if (nrow(tbl) == 0) {
    stop("`design_skeleton` did not contain any valid rows.", call. = FALSE)
  }
  if ("Group" %in% names(tbl)) {
    group_check <- tbl |>
      dplyr::distinct(.data$TemplatePerson, .data$Group) |>
      dplyr::count(.data$TemplatePerson, name = "n_groups")
    if (any(group_check$n_groups > 1L)) {
      stop("`design_skeleton` must assign at most one `Group` label per `TemplatePerson`.", call. = FALSE)
    }
  }
  if ("Weight" %in% names(tbl)) {
    tbl$Weight <- suppressWarnings(as.numeric(tbl$Weight))
    if (any(!is.finite(tbl$Weight) | tbl$Weight <= 0)) {
      stop("`design_skeleton$Weight` must contain positive finite values.", call. = FALSE)
    }
  }
  if (length(unique(tbl$Rater)) != n_rater) {
    stop("`design_skeleton` must reference exactly ", n_rater, " distinct rater levels.", call. = FALSE)
  }
  if (length(unique(tbl$Criterion)) != n_criterion) {
    stop("`design_skeleton` must reference exactly ", n_criterion, " distinct criterion levels.", call. = FALSE)
  }
  tbl
}

simulation_build_threshold_table <- function(thresholds, score_levels, step_span, model) {
  if (is.null(thresholds)) {
    est <- if (score_levels == 2L) 0 else seq(-abs(step_span), abs(step_span), length.out = score_levels - 1L)
    return(tibble::tibble(
      StepFacet = "Common",
      StepIndex = seq_along(est),
      Step = paste0("Step_", seq_along(est)),
      Estimate = as.numeric(est)
    ))
  }

  if (is.numeric(thresholds)) {
    est <- as.numeric(thresholds)
    if (length(est) != score_levels - 1L) {
      stop("Numeric `thresholds` must have length `score_levels - 1`.", call. = FALSE)
    }
    return(tibble::tibble(
      StepFacet = "Common",
      StepIndex = seq_along(est),
      Step = paste0("Step_", seq_along(est)),
      Estimate = est
    ))
  }

  if (!is.data.frame(thresholds)) {
    stop("`thresholds` must be NULL, a numeric vector, or a data.frame.", call. = FALSE)
  }

  tbl <- tibble::as_tibble(thresholds)
  if (!"StepFacet" %in% names(tbl) || !"Estimate" %in% names(tbl)) {
    stop("Threshold data frames must include `StepFacet` and `Estimate` columns.", call. = FALSE)
  }

  if (!"StepIndex" %in% names(tbl)) {
    if ("Step" %in% names(tbl)) {
      tbl$StepIndex <- suppressWarnings(as.integer(gsub("[^0-9]+", "", as.character(tbl$Step))))
    } else {
      stop("Threshold data frames must include either `StepIndex` or `Step`.", call. = FALSE)
    }
  }
  tbl$StepIndex <- suppressWarnings(as.integer(tbl$StepIndex))
  tbl$Estimate <- suppressWarnings(as.numeric(tbl$Estimate))
  tbl$StepFacet <- as.character(tbl$StepFacet)

  if (any(!is.finite(tbl$StepIndex)) || any(tbl$StepIndex < 1L)) {
    stop("`thresholds$StepIndex` must be positive integers.", call. = FALSE)
  }
  if (any(!is.finite(tbl$Estimate))) {
    stop("`thresholds$Estimate` must be finite numeric values.", call. = FALSE)
  }

  tbl <- tbl |>
    dplyr::mutate(Step = paste0("Step_", .data$StepIndex)) |>
    dplyr::arrange(.data$StepFacet, .data$StepIndex) |>
    dplyr::distinct(.data$StepFacet, .data$StepIndex, .keep_all = TRUE)

  counts <- tbl |>
    dplyr::count(.data$StepFacet, name = "n_steps")
  if (length(unique(counts$n_steps)) != 1L) {
    stop("Each `StepFacet` must supply the same number of step thresholds.", call. = FALSE)
  }
  expected_steps <- score_levels - 1L
  if (unique(counts$n_steps) != expected_steps) {
    stop("Threshold table implies ", unique(counts$n_steps) + 1L, " score levels, but `score_levels` = ", score_levels, ".", call. = FALSE)
  }

  if (identical(model, "RSM") && length(unique(tbl$StepFacet)) > 1L) {
    stop("`model = \"RSM\"` accepts only one common threshold set.", call. = FALSE)
  }

  tbl
}

simulation_extract_thresholds_from_fit <- function(step_tbl, model) {
  if (!is.data.frame(step_tbl) || nrow(step_tbl) == 0) {
    return(NULL)
  }
  tbl <- tibble::as_tibble(step_tbl)
  if (!"StepFacet" %in% names(tbl)) {
    tbl$StepFacet <- "Common"
  }
  if (!"StepIndex" %in% names(tbl)) {
    tbl$StepIndex <- suppressWarnings(as.integer(gsub("[^0-9]+", "", as.character(tbl$Step))))
  }
  tbl |>
    dplyr::transmute(
      StepFacet = as.character(.data$StepFacet),
      StepIndex = as.integer(.data$StepIndex),
      Step = paste0("Step_", .data$StepIndex),
      Estimate = as.numeric(.data$Estimate)
    ) |>
    dplyr::arrange(.data$StepFacet, .data$StepIndex)
}

simulation_spec_threshold_mode <- function(sim_spec) {
  tbl <- sim_spec$threshold_table %||% NULL
  if (!is.data.frame(tbl) || nrow(tbl) == 0) {
    return("implicit_common")
  }
  step_facets <- unique(as.character(tbl$StepFacet))
  if (length(step_facets) == 1L && identical(step_facets, "Common")) {
    "common"
  } else {
    "step_facet_specific"
  }
}

simulation_resolve_assignment <- function(base_assignment, n_rater, raters_per_person) {
  if (identical(base_assignment, "skeleton")) {
    return("skeleton")
  }
  if (identical(base_assignment, "resampled")) {
    return("resampled")
  }
  if (isTRUE(raters_per_person >= n_rater) && identical(base_assignment, "crossed")) {
    return("crossed")
  }
  "rotating"
}

simulation_override_spec_design <- function(sim_spec,
                                            n_person,
                                            n_rater,
                                            n_criterion,
                                            raters_per_person,
                                            group_levels = sim_spec$group_levels,
                                            dif_effects = sim_spec$dif_effects,
                                            interaction_effects = sim_spec$interaction_effects) {
  if (!inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must inherit from `mfrm_sim_spec`.", call. = FALSE)
  }

  n_person <- simulation_validate_count(n_person, "n_person", min_value = 2L)
  n_rater <- simulation_validate_count(n_rater, "n_rater", min_value = 2L)
  n_criterion <- simulation_validate_count(n_criterion, "n_criterion", min_value = 2L)
  raters_per_person <- simulation_validate_count(raters_per_person, "raters_per_person", min_value = 1L)
  if (raters_per_person > n_rater) {
    stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  }

  out <- sim_spec
  out$n_person <- n_person
  out$n_rater <- n_rater
  out$n_criterion <- n_criterion
  out$raters_per_person <- raters_per_person
  out$assignment <- simulation_resolve_assignment(sim_spec$assignment, n_rater = n_rater, raters_per_person = raters_per_person)
  out$group_levels <- group_levels
  out$dif_effects <- dif_effects
  out$interaction_effects <- interaction_effects

  if (identical(out$assignment, "resampled")) {
    if (!identical(n_rater, sim_spec$n_rater) || !identical(raters_per_person, sim_spec$raters_per_person)) {
      stop(
        "`assignment = \"resampled\"` reuses empirical person-level rater profiles. ",
        "It currently supports changing `n_person` only; keep `n_rater` and `raters_per_person` equal to the base specification.",
        call. = FALSE
      )
    }
  }
  if (identical(out$assignment, "skeleton")) {
    if (!identical(n_rater, sim_spec$n_rater) ||
        !identical(n_criterion, sim_spec$n_criterion) ||
        !identical(raters_per_person, sim_spec$raters_per_person)) {
      stop(
        "`assignment = \"skeleton\"` reuses the observed person-by-facet design skeleton. ",
        "It currently supports changing `n_person` only; keep `n_rater`, `n_criterion`, and `raters_per_person` equal to the base specification.",
        call. = FALSE
      )
    }
  }

  threshold_mode <- simulation_spec_threshold_mode(out)
  if (identical(threshold_mode, "step_facet_specific")) {
    expected_levels <- switch(
      out$step_facet,
      Criterion = sprintf("C%02d", seq_len(n_criterion)),
      Rater = sprintf("R%02d", seq_len(n_rater)),
      character(0)
    )
    observed_levels <- sort(unique(as.character(out$threshold_table$StepFacet)))
    if (!setequal(observed_levels, expected_levels)) {
      varying_arg <- if (identical(out$step_facet, "Criterion")) "n_criterion" else "n_rater"
      stop(
        "`sim_spec` contains step-facet-specific thresholds for `", out$step_facet,
        "` levels {", paste(observed_levels, collapse = ", "), "}. ",
        "Varying `", varying_arg, "` away from that threshold structure is not currently supported. ",
        "Use common thresholds or build a design-specific simulation specification.",
        call. = FALSE
      )
    }
  }

  out
}
