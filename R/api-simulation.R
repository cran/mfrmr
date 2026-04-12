#' Simulate long-format many-facet Rasch data for design studies
#'
#' @param n_person Number of persons/respondents.
#' @param n_rater Number of rater facet levels.
#' @param n_criterion Number of criterion/item facet levels.
#' @param raters_per_person Number of raters assigned to each person.
#' @param design Optional named design override supplied as a named list,
#'   named vector, or one-row data frame. When `sim_spec = NULL`, names may use
#'   canonical variables (`n_person`, `n_rater`, `n_criterion`,
#'   `raters_per_person`) or role keywords (`person`, `rater`, `criterion`,
#'   `assignment`). For the currently exposed facet keys, the schema-only
#'   future branch input `design$facets = c(person = ..., rater = ...,
#'   criterion = ...)` is also accepted. Do not specify the same variable
#'   through both `design` and the scalar count arguments.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread of step thresholds on the logit scale.
#' @param group_levels Optional character vector of group labels. When supplied,
#'   a balanced `Group` column is added to the simulated data.
#' @param dif_effects Optional data.frame describing true group-linked DIF
#'   effects. Must include `Group`, at least one design column such as
#'   `Criterion`, and numeric `Effect`.
#' @param interaction_effects Optional data.frame describing true non-group
#'   interaction effects. Must include at least one design column such as
#'   `Rater` or `Criterion`, plus numeric `Effect`.
#' @param seed Optional random seed.
#' @param model Measurement model recorded in the simulation setup. The current
#'   public generator supports `RSM`, `PCM`, and bounded `GPCM`.
#' @param step_facet Step facet used when `model = "PCM"` and threshold values
#'   vary across levels. Currently `"Criterion"` and `"Rater"` are supported.
#' @param slope_facet Slope facet used when `model = "GPCM"`. The current
#'   bounded `GPCM` branch requires `slope_facet == step_facet`.
#' @param thresholds Optional threshold specification. Use either a numeric
#'   vector of common thresholds or a data frame with columns `StepFacet`,
#'   `Step`/`StepIndex`, and `Estimate`.
#' @param slopes Optional slope specification used when `model = "GPCM"`.
#'   Use either a numeric vector aligned to the generated slope-facet levels or
#'   a data frame with columns `SlopeFacet` and `Estimate`. When omitted,
#'   slopes default to 1 for every slope-facet level, giving an exact `PCM`
#'   reduction.
#' @param assignment Assignment design. `"crossed"` means every person sees
#'   every rater; `"rotating"` uses a balanced rotating subset; `"resampled"`
#'   reuses person-level rater-assignment profiles stored in `sim_spec`;
#'   `"skeleton"` reuses an observed response skeleton stored in `sim_spec`,
#'   including optional `Group`/`Weight` columns when available. When omitted,
#'   the function chooses `"crossed"` if
#'   `raters_per_person == n_rater`, otherwise `"rotating"`.
#' @param sim_spec Optional output from [build_mfrm_sim_spec()] or
#'   [extract_mfrm_sim_spec()]. When supplied, it defines the generator setup;
#'   direct scalar arguments are treated as legacy inputs and should generally
#'   be left at their defaults except for `seed`. Any custom public two-facet
#'   names recorded in `sim_spec$facet_names` are also carried into the
#'   simulated output and downstream planning helpers. If `sim_spec` stores an
#'   active latent-regression population generator, the returned object also
#'   carries the generated one-row-per-person background-data table needed to
#'   refit that population model later.
#'
#' @details
#' This function generates synthetic MFRM data from the Rasch model.
#' The data-generating process is:
#'
#' 1. Draw person abilities: \eqn{\theta_n \sim N(0, \texttt{theta\_sd}^2)}
#' 2. Draw rater severities: \eqn{\delta_j \sim N(0, \texttt{rater\_sd}^2)}
#' 3. Draw criterion difficulties: \eqn{\beta_i \sim N(0, \texttt{criterion\_sd}^2)}
#' 4. Generate evenly-spaced step thresholds spanning \eqn{\pm}\code{step_span/2}
#' 5. For each observation, compute the linear predictor
#'    \eqn{\eta = \theta_n - \delta_j - \beta_i + \epsilon} where
#'    \eqn{\epsilon \sim N(0, \texttt{noise\_sd}^2)} (optional)
#' 6. Compute category probabilities under the recorded measurement model
#'    (`RSM`, `PCM`, or bounded `GPCM`) and sample the response
#'
#' Latent-value generation is explicit:
#' - `latent_distribution = "normal"` draws centered normal person/rater/
#'   criterion values using the supplied standard deviations
#' - `latent_distribution = "empirical"` resamples centered support values
#'   recorded in `sim_spec$empirical_support`
#' - if `sim_spec$population$active = TRUE`, person measures are generated from
#'   the stored latent-regression population model and template person
#'   covariates rather than from `theta_sd`
#'
#' When `dif_effects` is supplied, the specified logit shift is added to
#' \eqn{\eta} for the focal group on the target facet level, creating a
#' known DIF signal.  Similarly, `interaction_effects` injects a known
#' bias into specific facet-level combinations.
#'
#' The generator targets the common two-facet rating design (persons
#' \eqn{\times} raters \eqn{\times} criteria).  `raters_per_person`
#' controls the incomplete-block structure: when less than `n_rater`,
#' each person is assigned a rotating subset of raters to keep coverage
#' balanced and reproducible.
#'
#' Threshold handling is intentionally explicit:
#' - if `thresholds = NULL`, common equally spaced thresholds are generated
#'   from `step_span`
#' - if `thresholds` is a numeric vector, it is used as one common threshold set
#' - if `thresholds` is a data frame, threshold values may vary by `StepFacet`
#'   (currently `Criterion` or `Rater`)
#'
#' For bounded `GPCM`, the generator now requires an explicit slope
#' contract in parallel with the threshold table. The current public branch
#' keeps `slope_facet == step_facet` and uses the internal `category_prob_gpcm()`
#' helper for
#' response sampling. Broader design-planning helpers remain restricted until
#' that slope-aware contract is generalized beyond direct data generation.
#'
#' Assignment handling is also explicit:
#' - `"crossed"` uses the full person x rater x criterion design
#' - `"rotating"` assigns a deterministic rotating subset of raters per person
#' - `"resampled"` reuses empirical person-level rater profiles stored in
#'   `sim_spec$assignment_profiles`, optionally carrying over person-level
#'   `Group`
#' - `"skeleton"` reuses an observed person-by-rater-by-criterion response
#'   skeleton stored in `sim_spec$design_skeleton`, optionally carrying over
#'   `Group` and `Weight`
#'
#' For more controlled workflows, build a reusable simulation specification
#' first via [build_mfrm_sim_spec()] or derive one from an observed fit with
#' [extract_mfrm_sim_spec()], then pass it through `sim_spec`.
#'
#' Returned data include attributes:
#' - `mfrm_truth`: simulated true parameters (for parameter-recovery checks)
#' - `mfrm_truth$signals`: injected DIF and interaction signal tables
#' - `mfrm_truth$slope_table`: simulated discrimination table for bounded
#'   `GPCM`
#' - `mfrm_population_data`: generated one-row-per-person background data when
#'   the simulation specification stores an active latent-regression generator,
#'   including model-matrix xlevel and contrast provenance for categorical
#'   covariates
#' - `mfrm_simulation_spec`: generation settings (for reproducibility)
#'
#' @section Interpreting output:
#' - Higher `theta` values in `mfrm_truth$person` indicate higher person measures.
#' - Higher values in `mfrm_truth$facets$Rater` indicate more severe raters.
#' - Higher values in `mfrm_truth$facets$Criterion` indicate more difficult criteria.
#' - `mfrm_truth$signals$dif_effects` and `mfrm_truth$signals$interaction_effects`
#'   record any injected detection targets.
#'
#' @section Typical workflow:
#' 1. Generate one design with `simulate_mfrm_data()`.
#' 2. Fit with [fit_mfrm()] and diagnose with [diagnose_mfrm()].
#' 3. For repeated design studies, use [evaluate_mfrm_design()].
#'
#' @return A long-format `data.frame` with core columns `Study`, `Person`,
#'   two simulated non-person facet columns, and `Score`. By default those
#'   facet columns are `Rater` and `Criterion`; when `sim_spec` records custom
#'   public names, those names are used instead. If group labels are simulated
#'   or reused from an observed response skeleton, a `Group` column is
#'   included. If a weighted response skeleton is reused, a `Weight` column is
#'   also included.
#' @seealso [evaluate_mfrm_design()], [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' sim <- simulate_mfrm_data(
#'   n_person = 40,
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   seed = 123
#' )
#' head(sim)
#' names(attr(sim, "mfrm_truth"))
#' @export
simulate_mfrm_data <- function(n_person = 50,
                               n_rater = 4,
                               n_criterion = 4,
                               raters_per_person = n_rater,
                               design = NULL,
                               score_levels = 4,
                               theta_sd = 1,
                               rater_sd = 0.35,
                               criterion_sd = 0.25,
                               noise_sd = 0,
                               step_span = 1.4,
                               group_levels = NULL,
                               dif_effects = NULL,
                               interaction_effects = NULL,
                               seed = NULL,
                               model = c("RSM", "PCM", "GPCM"),
                               step_facet = "Criterion",
                               slope_facet = NULL,
                               thresholds = NULL,
                               slopes = NULL,
                               assignment = NULL,
                               sim_spec = NULL) {
  if (!is.null(sim_spec)) {
    if (!inherits(sim_spec, "mfrm_sim_spec")) {
      stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
    }
    n_person <- sim_spec$n_person
    n_rater <- sim_spec$n_rater
    n_criterion <- sim_spec$n_criterion
    raters_per_person <- sim_spec$raters_per_person
    score_levels <- sim_spec$score_levels
    theta_sd <- sim_spec$theta_sd
    rater_sd <- sim_spec$rater_sd
    criterion_sd <- sim_spec$criterion_sd
    noise_sd <- sim_spec$noise_sd
    step_span <- sim_spec$step_span
    group_levels <- sim_spec$group_levels
    dif_effects <- sim_spec$dif_effects
    interaction_effects <- sim_spec$interaction_effects
    model <- sim_spec$model
    step_facet <- sim_spec$step_facet
    slope_facet <- sim_spec$slope_facet %||% NULL
    thresholds <- sim_spec$threshold_table
    slopes <- sim_spec$slope_table %||% NULL
    assignment <- sim_spec$assignment
    latent_distribution <- as.character(sim_spec$latent_distribution %||% "normal")
    empirical_support <- sim_spec$empirical_support %||% NULL
    assignment_profiles <- sim_spec$assignment_profiles %||% NULL
    design_skeleton <- sim_spec$design_skeleton %||% NULL
    facet_names <- simulation_spec_output_facet_names(sim_spec)
    population_spec <- sim_spec$population %||% simulation_empty_population_spec()
  } else {
    supplied_counts <- intersect(
      names(as.list(match.call(expand.dots = FALSE))[-1]),
      c("n_person", "n_rater", "n_criterion", "raters_per_person")
    )
    score_levels <- as.integer(score_levels[1])
    model <- match.arg(toupper(as.character(model[1])), c("RSM", "PCM", "GPCM"))
    facet_names <- simulation_default_output_facet_names()
    design_counts <- simulation_resolve_design_counts(
      sim_spec = list(facet_names = facet_names),
      n_person = n_person,
      n_rater = n_rater,
      n_criterion = n_criterion,
      raters_per_person = raters_per_person,
      design = design,
      defaults = list(
        n_person = 50L,
        n_rater = 4L,
        n_criterion = 4L
      ),
      design_arg = "design",
      explicit_scalar_names = supplied_counts
    )
    n_person <- as.integer(design_counts$n_person[1])
    n_rater <- as.integer(design_counts$n_rater[1])
    n_criterion <- as.integer(design_counts$n_criterion[1])
    raters_per_person <- as.integer(design_counts$raters_per_person[1])
    assignment <- if (is.null(assignment)) {
      if (identical(raters_per_person, n_rater)) "crossed" else "rotating"
    } else {
      match.arg(tolower(as.character(assignment[1])), c("crossed", "rotating"))
    }
    if (identical(model, "GPCM")) {
      resolved_facets <- resolve_step_and_slope_facets(
        model = model,
        step_facet = step_facet[1] %||% "Criterion",
        slope_facet = slope_facet,
        facet_names = unname(facet_names)
      )
      step_facet <- resolved_facets$step_facet
      slope_facet <- resolved_facets$slope_facet
    } else {
      step_facet <- as.character(step_facet[1] %||% "Criterion")
      slope_facet <- NULL
    }
    thresholds <- simulation_build_threshold_table(
      thresholds = thresholds,
      score_levels = score_levels,
      step_span = as.numeric(step_span[1]),
      model = model
    )
    slopes <- simulation_build_slope_table(
      slopes = slopes,
      model = model,
      slope_facet = slope_facet,
      facet_names = stats::setNames(facet_names, c("rater", "criterion")),
      n_rater = n_rater,
      n_criterion = n_criterion
    )
    latent_distribution <- "normal"
    empirical_support <- NULL
    assignment_profiles <- NULL
    design_skeleton <- NULL
    population_spec <- simulation_empty_population_spec()
  }

  if (!is.finite(n_person) || n_person < 2L) stop("`n_person` must be >= 2.", call. = FALSE)
  if (!is.finite(n_rater) || n_rater < 2L) stop("`n_rater` must be >= 2.", call. = FALSE)
  if (!is.finite(n_criterion) || n_criterion < 2L) stop("`n_criterion` must be >= 2.", call. = FALSE)
  if (!is.finite(raters_per_person) || raters_per_person < 1L) stop("`raters_per_person` must be >= 1.", call. = FALSE)
  if (raters_per_person > n_rater) stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  if (!is.finite(score_levels) || score_levels < 2L) stop("`score_levels` must be >= 2.", call. = FALSE)
  if (!step_facet %in% facet_names) {
    stop(
      "`step_facet` must match one of the simulated facet columns: ",
      paste(facet_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!identical(assignment, "crossed") && !identical(assignment, "rotating") &&
      !identical(assignment, "resampled") && !identical(assignment, "skeleton")) {
    stop("`assignment` must be one of \"crossed\", \"rotating\", \"resampled\", or \"skeleton\".", call. = FALSE)
  }
  if (identical(assignment, "crossed") && raters_per_person != n_rater) {
    stop("`assignment = \"crossed\"` requires `raters_per_person == n_rater`.", call. = FALSE)
  }
  if (!identical(latent_distribution, "normal") && !identical(latent_distribution, "empirical")) {
    stop("`latent_distribution` must be either \"normal\" or \"empirical\".", call. = FALSE)
  }
  if (identical(latent_distribution, "empirical") && is.null(empirical_support)) {
    stop("`latent_distribution = \"empirical\"` requires empirical support values in `sim_spec`.", call. = FALSE)
  }
  if (identical(assignment, "resampled") && is.null(assignment_profiles)) {
    stop("`assignment = \"resampled\"` requires assignment profiles in `sim_spec`.", call. = FALSE)
  }
  if (identical(assignment, "skeleton") && is.null(design_skeleton)) {
    stop("`assignment = \"skeleton\"` requires a design skeleton in `sim_spec`.", call. = FALSE)
  }

  if (!is.null(group_levels)) {
    group_levels <- as.character(group_levels)
    group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
    group_levels <- unique(group_levels)
    if (length(group_levels) < 1L) stop("`group_levels` must contain at least one non-empty label.", call. = FALSE)
  }

  with_preserved_rng_seed(seed, {
    assignment_facet <- facet_names[1]
    criterion_facet <- facet_names[2]
    expand_design <- function(person_values, assignment_values, criterion_values) {
      expand.grid(
        setNames(
          list(person_values, assignment_values, criterion_values),
          c("Person", assignment_facet, criterion_facet)
        ),
        stringsAsFactors = FALSE
      )
    }
    person_ids <- sprintf("P%03d", seq_len(n_person))
    rater_ids <- if (identical(assignment, "skeleton") && !is.null(design_skeleton)) {
      sort(unique(as.character(design_skeleton$Rater)))
    } else if (identical(assignment, "resampled") && !is.null(assignment_profiles)) {
      sort(unique(as.character(assignment_profiles$Rater)))
    } else {
      simulation_spec_role_levels(sim_spec, "rater", count = n_rater)
    }
    criterion_ids <- if (identical(assignment, "skeleton") && !is.null(design_skeleton)) {
      sort(unique(as.character(design_skeleton$Criterion)))
    } else {
      simulation_spec_role_levels(sim_spec, "criterion", count = n_criterion)
    }

    population_bundle <- if (isTRUE(population_spec$active)) {
      simulation_generate_population_person_data(
        population_spec = population_spec,
        person_ids = person_ids
      )
    } else {
      NULL
    }
    theta <- if (isTRUE(population_spec$active)) {
      population_bundle$theta
    } else {
      stats::rnorm(n_person, mean = 0, sd = as.numeric(theta_sd[1]))
    }
    names(theta) <- person_ids

    rater_effects <- sort(stats::rnorm(n_rater, mean = 0, sd = as.numeric(rater_sd[1])))
    names(rater_effects) <- rater_ids

    criterion_effects <- sort(stats::rnorm(n_criterion, mean = 0, sd = as.numeric(criterion_sd[1])))
    names(criterion_effects) <- criterion_ids

    if (identical(latent_distribution, "empirical")) {
      if (!isTRUE(population_spec$active)) {
        theta <- simulation_sample_empirical_support(empirical_support$person, n_person)
        names(theta) <- person_ids
      }
      rater_effects <- simulation_sample_empirical_support(empirical_support$rater, n_rater)
      names(rater_effects) <- rater_ids
      criterion_effects <- simulation_sample_empirical_support(empirical_support$criterion, n_criterion)
      names(criterion_effects) <- criterion_ids
    }

    if (score_levels == 2L) {
      steps <- thresholds$Estimate[thresholds$StepFacet %in% c("Common", unique(thresholds$StepFacet)[1])][1]
    }

    if (assignment == "crossed") {
      dat <- expand_design(person_ids, rater_ids, criterion_ids)
    } else if (assignment == "skeleton") {
      dat <- simulation_generate_skeleton_assignment(
        person_ids = person_ids,
        design_skeleton = design_skeleton,
        facet_names = facet_names
      )
    } else if (assignment == "resampled") {
      dat <- simulation_generate_resampled_assignment(
        person_ids = person_ids,
        criterion_ids = criterion_ids,
        assignment_profiles = assignment_profiles,
        facet_names = facet_names
      )
    } else {
      rows <- vector("list", length(person_ids))
      for (i in seq_along(person_ids)) {
        assigned <- rater_ids[((i - 1L) + seq_len(raters_per_person) - 1L) %% n_rater + 1L]
        rows[[i]] <- expand_design(person_ids[i], assigned, criterion_ids)
      }
      dat <- dplyr::bind_rows(rows)
    }

  threshold_table <- tibble::as_tibble(thresholds) |>
      dplyr::mutate(
        StepFacet = as.character(.data$StepFacet),
        StepIndex = as.integer(.data$StepIndex),
        Step = as.character(.data$Step),
        Estimate = as.numeric(.data$Estimate)
      ) |>
      dplyr::arrange(.data$StepFacet, .data$StepIndex)
    threshold_lookup <- split(threshold_table$Estimate, threshold_table$StepFacet)
    if (!"Common" %in% names(threshold_lookup) && identical(model, "RSM")) {
      stop("RSM simulation requires one common threshold set.", call. = FALSE)
    }
    slope_table <- if (identical(model, "GPCM")) {
      tibble::as_tibble(slopes) |>
        dplyr::mutate(
          SlopeFacet = as.character(.data$SlopeFacet),
          Estimate = as.numeric(.data$Estimate)
        ) |>
        dplyr::arrange(.data$SlopeFacet)
    } else {
      NULL
    }
    slope_lookup <- if (identical(model, "GPCM")) {
      stats::setNames(as.numeric(slope_table$Estimate), as.character(slope_table$SlopeFacet))
    } else {
      NULL
    }

    if ("Group" %in% names(dat)) {
      group_assign_tbl <- dat |>
        dplyr::distinct(.data$Person, .data$Group)
      group_assign <- stats::setNames(as.character(group_assign_tbl$Group),
                                      as.character(group_assign_tbl$Person))
    } else if (!is.null(group_levels)) {
      group_assign <- rep(group_levels, length.out = n_person)
      group_assign <- sample(group_assign, size = n_person, replace = FALSE)
      names(group_assign) <- person_ids
      dat$Group <- unname(group_assign[dat$Person])
    } else {
      group_assign <- NULL
    }

    allowed_effect_cols <- intersect(c("Group", "Person", facet_names), names(dat))
    dif_effects <- simulation_normalize_effects(
      effects = dif_effects,
      arg_name = "dif_effects",
      allowed_cols = allowed_effect_cols
    )
    if (nrow(dif_effects) > 0) {
      if (!"Group" %in% names(dif_effects)) {
        stop("`dif_effects` must include a `Group` column.", call. = FALSE)
      }
      if (!any(c("Person", facet_names) %in% names(dif_effects))) {
        stop(
          "`dif_effects` must include at least one of `Person`, `",
          paste(facet_names, collapse = "`, or `"),
          "`.",
          call. = FALSE
        )
      }
    }

    interaction_effects <- simulation_normalize_effects(
      effects = interaction_effects,
      arg_name = "interaction_effects",
      allowed_cols = allowed_effect_cols
    )
    if (nrow(interaction_effects) > 0 && !any(c("Person", facet_names) %in% names(interaction_effects))) {
      stop(
        "`interaction_effects` must include at least one of `Person`, `",
        paste(facet_names, collapse = "`, or `"),
        "`.",
        call. = FALSE
      )
    }

    eta <- theta[dat$Person] - rater_effects[dat[[assignment_facet]]] - criterion_effects[dat[[criterion_facet]]]
    eta <- eta + simulation_apply_effects(dat, dif_effects)
    eta <- eta + simulation_apply_effects(dat, interaction_effects)
    if (isTRUE(is.finite(noise_sd)) && as.numeric(noise_sd[1]) > 0) {
      eta <- eta + stats::rnorm(length(eta), mean = 0, sd = as.numeric(noise_sd[1]))
    }

    threshold_key <- if ("Common" %in% names(threshold_lookup)) {
      rep("Common", nrow(dat))
    } else {
      as.character(dat[[step_facet]])
    }
    if (!all(threshold_key %in% names(threshold_lookup))) {
      missing_keys <- setdiff(unique(threshold_key), names(threshold_lookup))
      stop(
        "Threshold specification is missing step-facet levels required for simulation: ",
        paste(missing_keys, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if (identical(model, "PCM") || identical(model, "GPCM")) {
      step_levels <- if ("Common" %in% names(threshold_lookup)) {
        sort(unique(as.character(dat[[step_facet]])))
      } else {
        unique(as.character(threshold_table$StepFacet))
      }
      step_cum_mat <- t(vapply(step_levels, function(level) {
        step_vec <- if ("Common" %in% names(threshold_lookup)) {
          threshold_lookup[["Common"]]
        } else {
          threshold_lookup[[level]]
        }
        c(0, cumsum(step_vec))
      }, numeric(score_levels)))
      criterion_idx <- match(as.character(dat[[step_facet]]), step_levels)
      if (anyNA(criterion_idx)) {
        stop("PCM simulation could not align observations to step-facet thresholds.", call. = FALSE)
      }
      criterion_splits <- split(seq_along(criterion_idx), criterion_idx)
      if (identical(model, "GPCM")) {
        slope_levels <- sort(unique(as.character(dat[[slope_facet]])))
        if (!all(slope_levels %in% names(slope_lookup))) {
          missing_levels <- setdiff(slope_levels, names(slope_lookup))
          stop(
            "Slope specification is missing slope-facet levels required for simulation: ",
            paste(missing_levels, collapse = ", "),
            ".",
            call. = FALSE
          )
        }
        slope_vals <- unname(slope_lookup[slope_levels])
        slope_idx <- match(as.character(dat[[slope_facet]]), slope_levels)
        if (anyNA(slope_idx)) {
          stop("GPCM simulation could not align observations to slope-facet values.", call. = FALSE)
        }
        prob_mat <- category_prob_gpcm(
          eta = eta,
          step_cum_mat = step_cum_mat,
          criterion_idx = criterion_idx,
          slopes = slope_vals,
          slope_idx = slope_idx
        )
      } else {
        prob_mat <- category_prob_pcm(
          eta = eta,
          step_cum_mat = step_cum_mat,
          criterion_idx = criterion_idx,
          criterion_splits = criterion_splits
        )
      }
    } else {
      step_cum <- c(0, cumsum(threshold_lookup[["Common"]]))
      prob_mat <- category_prob_rsm(eta, step_cum)
    }
    dat$Score <- apply(prob_mat, 1, function(p) sample.int(score_levels, size = 1L, prob = p))
    dat$Study <- "SimulatedDesign"
    keep_cols <- c("Study", "Person", facet_names, "Score")
    if ("Group" %in% names(dat)) keep_cols <- c(keep_cols, "Group")
    if ("Weight" %in% names(dat)) keep_cols <- c(keep_cols, "Weight")
    dat <- dat[, keep_cols]

    truth_facets <- setNames(list(rater_effects, criterion_effects), facet_names)
    attr(dat, "mfrm_truth") <- list(
      person = theta,
      facets = truth_facets,
      steps = if ("Common" %in% names(threshold_lookup)) threshold_lookup$Common else threshold_table,
      step_table = threshold_table,
      slopes = if (identical(model, "GPCM")) unname(slope_lookup[as.character(slope_table$SlopeFacet)]) else NULL,
      slope_table = slope_table,
      population = if (isTRUE(population_spec$active)) {
        list(
          formula = population_spec$formula,
          coefficients = population_spec$coefficients,
          sigma2 = population_spec$sigma2,
          design_columns = population_spec$design_columns,
          xlevels = population_spec$xlevels,
          contrasts = population_spec$contrasts,
          linear_predictor = population_bundle$linear_predictor,
          person_data = population_bundle$person_data
        )
      } else {
        NULL
      },
      groups = group_assign,
      signals = list(
        dif_effects = dif_effects,
        interaction_effects = interaction_effects
      )
    )
    attr(dat, "mfrm_population_data") <- if (isTRUE(population_spec$active)) {
      list(
        active = TRUE,
        person_data = population_bundle$person_data,
        population_formula = population_spec$formula,
        person_id = "Person",
        population_policy = "error",
        design_columns = population_spec$design_columns,
        xlevels = population_spec$xlevels,
        contrasts = population_spec$contrasts,
        coefficients = population_spec$coefficients,
        sigma2 = population_spec$sigma2
      )
    } else {
      NULL
    }
    attr(dat, "mfrm_simulation_spec") <- list(
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
      group_levels = group_levels,
      model = model,
      step_facet = step_facet,
      slope_facet = slope_facet,
      facet_names = stats::setNames(facet_names, c("rater", "criterion")),
      facet_levels = list(rater = rater_ids, criterion = criterion_ids),
      assignment = assignment,
      latent_distribution = latent_distribution,
      threshold_table = threshold_table,
      slope_table = slope_table,
      design_skeleton = design_skeleton,
      population = population_spec
    )
    dat
  })
}

simulation_center_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(numeric(0))
  x - mean(x)
}

simulation_sample_empirical_support <- function(values, size) {
  centered <- simulation_center_numeric(values)
  if (length(centered) < 2L) {
    stop("Empirical latent support must contain at least two finite values.", call. = FALSE)
  }
  sample(centered, size = size, replace = TRUE)
}

simulation_generate_population_person_data <- function(population_spec, person_ids) {
  if (!isTRUE(population_spec$active)) {
    return(NULL)
  }

  template <- as.data.frame(population_spec$covariate_template %||% data.frame(), stringsAsFactors = FALSE)
  template_id <- as.character(population_spec$person_id[1] %||% "TemplatePerson")
  if (!template_id %in% names(template)) {
    stop("Stored population covariate template is missing its person-ID column.", call. = FALSE)
  }
  draw_idx <- sample.int(nrow(template), size = length(person_ids), replace = TRUE)
  sampled <- template[draw_idx, , drop = FALSE]
  sampled[[template_id]] <- person_ids
  person_data <- sampled
  names(person_data)[names(person_data) == template_id] <- "Person"

  model_frame_args <- list(
    formula = population_spec$formula,
    data = person_data,
    na.action = stats::na.pass
  )
  if (!is.null(population_spec$xlevels)) {
    model_frame_args$xlev <- population_spec$xlevels
  }
  mf <- do.call(stats::model.frame, model_frame_args)
  for (nm in intersect(names(mf), names(person_data))) {
    person_data[[nm]] <- mf[[nm]]
  }
  terms_obj <- attr(mf, "terms") %||% stats::terms(population_spec$formula)
  x_mat <- stats::model.matrix(
    terms_obj,
    data = mf,
    contrasts.arg = population_spec$contrasts %||% NULL
  )
  coeff <- as.numeric(population_spec$coefficients[population_spec$design_columns])
  mu <- drop(x_mat %*% coeff)
  theta <- mu + stats::rnorm(length(person_ids), mean = 0, sd = sqrt(population_spec$sigma2))
  names(theta) <- person_ids

  list(
    theta = theta,
    linear_predictor = stats::setNames(as.numeric(mu), person_ids),
    person_data = person_data
  )
}

simulation_generate_resampled_assignment <- function(person_ids,
                                                     criterion_ids,
                                                     assignment_profiles,
                                                     facet_names = simulation_default_output_facet_names()) {
  if (!is.data.frame(assignment_profiles) || nrow(assignment_profiles) == 0) {
    stop("Empirical assignment profiles are required for `assignment = \"resampled\"`.",
         call. = FALSE)
  }

  assignment_facet <- facet_names[1]
  criterion_facet <- facet_names[2]
  profiles <- tibble::as_tibble(assignment_profiles) |>
    dplyr::distinct(.data$TemplatePerson, .data$Rater) |>
    dplyr::arrange(.data$TemplatePerson, .data$Rater)
  group_map <- NULL
  if ("Group" %in% names(assignment_profiles)) {
    group_map <- tibble::as_tibble(assignment_profiles) |>
      dplyr::filter(!is.na(.data$Group), nzchar(.data$Group)) |>
      dplyr::distinct(.data$TemplatePerson, .data$Group)
  }
  template_people <- unique(profiles$TemplatePerson)
  profile_lookup <- split(
    profiles$Rater,
    factor(profiles$TemplatePerson, levels = template_people)
  )
  group_lookup <- if (is.null(group_map)) {
    NULL
  } else {
    split(
      group_map$Group,
      factor(group_map$TemplatePerson, levels = template_people)
    )
  }
  sampled_templates <- sample(template_people, size = length(person_ids), replace = TRUE)

  rows <- vector("list", length(person_ids))
  for (i in seq_along(person_ids)) {
    assigned <- profile_lookup[[sampled_templates[i]]]
    rows[[i]] <- expand.grid(
      setNames(
        list(person_ids[i], assigned, criterion_ids),
        c("Person", assignment_facet, criterion_facet)
      ),
      stringsAsFactors = FALSE
    )
    if (!is.null(group_lookup)) {
      matched_group <- group_lookup[[sampled_templates[i]]][1]
      if (is.character(matched_group) && nzchar(matched_group)) {
        rows[[i]]$Group <- matched_group
      }
    }
  }
  dplyr::bind_rows(rows)
}

simulation_generate_skeleton_assignment <- function(person_ids,
                                                   design_skeleton,
                                                   facet_names = simulation_default_output_facet_names()) {
  if (!is.data.frame(design_skeleton) || nrow(design_skeleton) == 0) {
    stop("Observed design skeleton is required for `assignment = \"skeleton\"`.",
         call. = FALSE)
  }

  assignment_facet <- facet_names[1]
  criterion_facet <- facet_names[2]
  skeleton <- tibble::as_tibble(design_skeleton) |>
    dplyr::distinct()
  template_people <- unique(skeleton$TemplatePerson)
  template_lookup <- split(
    skeleton,
    factor(skeleton$TemplatePerson, levels = template_people)
  )
  sampled_templates <- sample(template_people, size = length(person_ids), replace = TRUE)
  keep_group <- "Group" %in% names(skeleton)
  keep_weight <- "Weight" %in% names(skeleton)
  select_cols <- c("Rater", "Criterion")
  if (keep_group) select_cols <- c(select_cols, "Group")
  if (keep_weight) select_cols <- c(select_cols, "Weight")

  rows <- vector("list", length(person_ids))
  for (i in seq_along(person_ids)) {
    template_rows <- tibble::as_tibble(template_lookup[[sampled_templates[i]]])
    rows[[i]] <- dplyr::mutate(
      template_rows[, select_cols, drop = FALSE],
      Person = person_ids[i],
      .before = 1
    )
    names(rows[[i]])[names(rows[[i]]) == "Rater"] <- assignment_facet
    names(rows[[i]])[names(rows[[i]]) == "Criterion"] <- criterion_facet
  }
  dplyr::bind_rows(rows)
}

simulation_normalize_effects <- function(effects, arg_name, allowed_cols) {
  if (is.null(effects)) return(tibble::tibble())
  if (!is.data.frame(effects)) {
    stop("`", arg_name, "` must be a data.frame with an `Effect` column.", call. = FALSE)
  }
  eff <- tibble::as_tibble(effects)
  if (nrow(eff) == 0L && ncol(eff) == 0L) {
    return(tibble::tibble())
  }
  if (!"Effect" %in% names(eff)) {
    stop("`", arg_name, "` must include an `Effect` column.", call. = FALSE)
  }
  unknown <- setdiff(names(eff), c(allowed_cols, "Effect"))
  if (length(unknown) > 0) {
    stop("`", arg_name, "` includes unsupported columns: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  key_cols <- setdiff(names(eff), "Effect")
  if (length(key_cols) == 0) {
    stop("`", arg_name, "` must include at least one matching design column besides `Effect`.", call. = FALSE)
  }
  eff$Effect <- suppressWarnings(as.numeric(eff$Effect))
  if (any(!is.finite(eff$Effect))) {
    stop("`", arg_name, "` has non-finite values in `Effect`.", call. = FALSE)
  }
  eff
}

simulation_apply_effects <- function(dat, effects) {
  if (is.null(effects) || !is.data.frame(effects) || nrow(effects) == 0) {
    return(rep(0, nrow(dat)))
  }
  key_cols <- setdiff(names(effects), "Effect")
  adj <- numeric(nrow(dat))
  for (i in seq_len(nrow(effects))) {
    mask <- rep(TRUE, nrow(dat))
    for (col in key_cols) {
      value <- effects[[col]][i]
      if (is.na(value) || !nzchar(as.character(value))) next
      mask <- mask & as.character(dat[[col]]) == as.character(value)
    }
    adj[mask] <- adj[mask] + as.numeric(effects$Effect[i])
  }
  adj
}

design_eval_extract_truth <- function(truth, facet) {
  if (is.null(truth) || !is.list(truth)) return(NULL)
  if (identical(facet, "Person")) return(truth$person)
  if (is.null(truth$facets) || !facet %in% names(truth$facets)) return(NULL)
  truth$facets[[facet]]
}

design_eval_safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

design_eval_safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

simulation_mcse_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(n)
}

simulation_mcse_proportion <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n <= 1L) return(NA_real_)
  p <- mean(x, na.rm = TRUE)
  sqrt(p * (1 - p) / n)
}

simulation_threshold_summary <- function(sim_spec) {
  if (is.null(sim_spec)) {
    return(list(mode = "generated_common", step_facet = NA_character_))
  }
  list(
    mode = simulation_spec_threshold_mode(sim_spec),
    step_facet = as.character(sim_spec$step_facet %||% NA_character_)
  )
}

simulation_resolve_fit_step_facet <- function(model, step_facet, generator_step_facet) {
  if (!identical(model, "PCM")) {
    return(NA_character_)
  }
  explicit <- as.character(step_facet[1] %||% NA_character_)
  if (is.na(explicit) || !nzchar(explicit)) {
    inherited <- as.character(generator_step_facet[1] %||% NA_character_)
    if (!is.na(inherited) && nzchar(inherited)) {
      return(inherited)
    }
    return("Criterion")
  }
  explicit
}

simulation_recovery_contract <- function(generator_model,
                                         generator_step_facet,
                                         fitted_model,
                                         fitted_step_facet) {
  same_model <- identical(as.character(generator_model), as.character(fitted_model))
  if (!same_model) {
    return(list(
      comparable = FALSE,
      basis = "generator_fit_model_mismatch"
    ))
  }
  if (identical(as.character(fitted_model), "PCM")) {
    same_step_facet <- identical(
      as.character(generator_step_facet %||% NA_character_),
      as.character(fitted_step_facet %||% NA_character_)
    )
    if (!same_step_facet) {
      return(list(
        comparable = FALSE,
        basis = "generator_fit_step_facet_mismatch"
      ))
    }
  }
  list(
    comparable = TRUE,
    basis = "generator_fit_contract_aligned"
  )
}

simulation_design_variable_slug <- function(label) {
  value <- as.character(label[1] %||% "")
  value <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", value, perl = TRUE)
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  value <- gsub("_+", "_", value)
  value <- tolower(gsub("^_|_$", "", value))
  if (!nzchar(value)) {
    return("facet")
  }
  value
}

simulation_design_variable_aliases <- function(sim_spec = NULL) {
  facet_names <- simulation_spec_output_facet_names(sim_spec)
  rater_alias <- if (identical(facet_names[1], "Rater")) {
    "n_rater"
  } else {
    paste0("n_", simulation_design_variable_slug(facet_names[1]))
  }
  criterion_alias <- if (identical(facet_names[2], "Criterion")) {
    "n_criterion"
  } else {
    paste0("n_", simulation_design_variable_slug(facet_names[2]))
  }
  assignment_alias <- if (identical(facet_names[1], "Rater")) {
    "raters_per_person"
  } else {
    paste0(simulation_design_variable_slug(facet_names[1]), "_per_person")
  }
  stats::setNames(
    c("n_person", rater_alias, criterion_alias, assignment_alias),
    c("n_person", "n_rater", "n_criterion", "raters_per_person")
  )
}

simulation_design_descriptor <- function(sim_spec = NULL) {
  facet_names <- simulation_spec_output_facet_names(sim_spec)
  aliases <- simulation_design_variable_aliases(sim_spec)
  tibble::tibble(
    role = c("person", "rater", "criterion", "assignment"),
    canonical = c("n_person", "n_rater", "n_criterion", "raters_per_person"),
    alias = unname(aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")]),
    facet = c("Person", facet_names[1], facet_names[2], facet_names[1]),
    quantity = c("persons", "facet_levels", "facet_levels", "assignments_per_person"),
    description = c(
      "Number of persons in the design condition.",
      paste0("Number of ", facet_names[1], " levels in the design condition."),
      paste0("Number of ", facet_names[2], " levels in the design condition."),
      paste0("Number of ", facet_names[1], " assignments per person in the design condition.")
    )
  )
}

simulation_facet_manifest <- function(sim_spec = NULL, facet_names = NULL) {
  if (is.null(facet_names)) {
    facet_names <- simulation_spec_output_facet_names(sim_spec)
    n_person <- suppressWarnings(as.integer(sim_spec$n_person %||% NA_integer_))
    n_rater <- suppressWarnings(as.integer(sim_spec$n_rater %||% NA_integer_))
    n_criterion <- suppressWarnings(as.integer(sim_spec$n_criterion %||% NA_integer_))
    rater_count <- if (is.finite(n_rater)) n_rater else NULL
    criterion_count <- if (is.finite(n_criterion)) n_criterion else NULL
    rater_levels <- simulation_spec_role_levels(sim_spec, role = "rater", count = rater_count)
    criterion_levels <- simulation_spec_role_levels(sim_spec, role = "criterion", count = criterion_count)
    if (length(rater_levels) > 0L) {
      n_rater <- length(rater_levels)
    }
    if (length(criterion_levels) > 0L) {
      n_criterion <- length(criterion_levels)
    }
    count_values <- c(
      n_person = n_person,
      n_rater = n_rater,
      n_criterion = n_criterion
    )
    aliases <- simulation_design_variable_aliases(sim_spec)
  } else {
    facet_names <- simulation_validate_output_facet_names(facet_names)
    count_values <- c(n_person = NA_integer_, n_rater = NA_integer_, n_criterion = NA_integer_)
    aliases <- simulation_design_variable_aliases(list(facet_names = facet_names))
  }
  facet_names <- unname(as.character(facet_names))

  tibble::tibble(
    planner_role = c("person", "rater", "criterion"),
    facet = c("Person", facet_names[1], facet_names[2]),
    facet_kind = c("person", "non_person", "non_person"),
    level_count = unname(count_values[c("n_person", "n_rater", "n_criterion")]),
    planning_count_variable = c("n_person", "n_rater", "n_criterion"),
    planning_count_alias = unname(aliases[c("n_person", "n_rater", "n_criterion")]),
    current_planner_role_supported = TRUE,
    arbitrary_facet_branch_candidate = TRUE
  )
}

simulation_future_facet_table <- function(sim_spec = NULL, facet_names = NULL) {
  facet_manifest <- if (is.null(facet_names)) {
    simulation_facet_manifest(sim_spec)
  } else {
    simulation_facet_manifest(facet_names = facet_names)
  }
  facets <- unname(as.character(facet_manifest$facet))
  facet_kinds <- unname(as.character(facet_manifest$facet_kind))

  tibble::tibble(
    facet = facets,
    facet_kind = facet_kinds,
    level_count = unname(facet_manifest$level_count),
    future_facet_key = unname(vapply(facets, simulation_design_variable_slug, character(1))),
    future_axis_class = ifelse(
      facet_kinds == "person",
      "person_count",
      "facet_level_count"
    ),
    current_planning_count_variable = as.character(facet_manifest$planning_count_variable),
    current_planning_count_alias = as.character(facet_manifest$planning_count_alias),
    current_planner_role_supported = as.logical(facet_manifest$current_planner_role_supported),
    arbitrary_facet_branch_candidate = as.logical(facet_manifest$arbitrary_facet_branch_candidate),
    branch_stage = "schema_only"
  )
}

simulation_future_design_template <- function(sim_spec = NULL, facet_names = NULL) {
  future_facet_table <- if (is.null(facet_names)) {
    simulation_future_facet_table(sim_spec)
  } else {
    simulation_future_facet_table(facet_names = facet_names)
  }

  facet_counts <- stats::setNames(
    as.list(unname(future_facet_table$level_count %||% rep(NA_integer_, nrow(future_facet_table)))),
    as.character(future_facet_table$future_facet_key)
  )

  assignment_value <- NA_integer_
  if (!is.null(sim_spec)) {
    assignment_value <- as.integer(sim_spec$raters_per_person %||% NA_integer_)
  }

  list(
    facets = facet_counts,
    assignment = assignment_value,
    note = paste(
      "Schema-only stub mirroring the current design through the future-branch",
      "`design$facets(named counts)` contract while keeping the current",
      "assignment axis at top level."
    )
  )
}

simulation_future_branch_design_schema <- function(sim_spec = NULL, facet_names = NULL) {
  future_facet_table <- if (is.null(facet_names)) {
    simulation_future_facet_table(sim_spec)
  } else {
    simulation_future_facet_table(facet_names = facet_names)
  }
  future_design_template <- if (is.null(facet_names)) {
    simulation_future_design_template(sim_spec)
  } else {
    simulation_future_design_template(facet_names = facet_names)
  }

  alias_source <- if (!is.null(sim_spec)) {
    sim_spec
  } else if (!is.null(facet_names)) {
    list(facet_names = simulation_validate_output_facet_names(facet_names))
  } else {
    NULL
  }

  facet_axes <- future_facet_table |>
    dplyr::transmute(
      input_key = .data$future_facet_key,
      facet = .data$facet,
      facet_kind = .data$facet_kind,
      axis_class = .data$future_axis_class,
      level_count = .data$level_count,
      canonical_design_variable = .data$current_planning_count_variable,
      public_design_alias = .data$current_planning_count_alias,
      current_planner_role_supported = .data$current_planner_role_supported,
      arbitrary_facet_branch_candidate = .data$arbitrary_facet_branch_candidate
    )

  assignment_axis <- list(
    current_planning_count_variable = "raters_per_person",
    current_planning_count_alias = unname(simulation_design_variable_aliases(alias_source)[["raters_per_person"]]),
    future_input_key = "assignment",
    axis_class = "assignment_count",
    depends_on_input_key = as.character(facet_axes$input_key[match("facet_level_count", facet_axes$axis_class)] %||% "rater"),
    depends_on_canonical_design_variable = as.character(
      facet_axes$canonical_design_variable[match("facet_level_count", facet_axes$axis_class)] %||% "n_rater"
    ),
    depends_on_public_design_alias = as.character(
      facet_axes$public_design_alias[match("facet_level_count", facet_axes$axis_class)] %||% "n_rater"
    )
  )

  schema <- list(
    schema_contract = "arbitrary_facet_design_schema",
    schema_stage = "schema_only",
    facet_axes = facet_axes,
    assignment_axis = assignment_axis,
    input_keys = c(as.character(facet_axes$input_key), assignment_axis$future_input_key),
    canonical_design_variables = c(
      as.character(facet_axes$canonical_design_variable),
      assignment_axis$current_planning_count_variable
    ),
    default_design = future_design_template,
    note = paste(
      "Schema-only design-schema object for the future arbitrary-facet branch,",
      "bundling stable facet-count axes and the assignment axis in one",
      "machine-readable contract."
    )
  )
  schema$grid_semantics <- simulation_future_branch_grid_semantics(schema)
  schema
}

simulation_future_branch_grid_semantics <- function(design_schema) {
  design_schema <- simulation_coerce_future_branch_design_schema(design_schema)
  facet_axes <- design_schema$facet_axes
  assignment_axis <- design_schema$assignment_axis
  assignment_var <- as.character(
    assignment_axis$current_planning_count_variable %||% "raters_per_person"
  )
  assignment_key <- as.character(
    assignment_axis$future_input_key %||% "assignment"
  )
  canonical_columns <- c(
    "design_id",
    as.character(facet_axes$canonical_design_variable),
    assignment_var
  )
  public_columns <- c(
    "design_id",
    as.character(facet_axes$public_design_alias),
    as.character(assignment_axis$current_planning_count_alias)
  )
  branch_columns <- c(
    "design_id",
    as.character(facet_axes$input_key),
    assignment_key
  )

  list(
    semantics_contract = "arbitrary_facet_design_grid_semantics",
    semantics_stage = as.character(design_schema$schema_stage %||% "schema_only"),
    id_variable = "design_id",
    canonical_columns = canonical_columns,
    public_columns = public_columns,
    branch_columns = branch_columns,
    feasibility_rule = paste0(
      assignment_var,
      " <= ",
      as.character(assignment_axis$depends_on_canonical_design_variable %||% "n_rater")
    ),
    row_meaning = paste(
      "Each row is one schema-only future-branch design condition formed by",
      "crossing facet-count axes and the assignment axis after applying the",
      "current feasibility rule."
    ),
    default_id_prefix = "F"
  )
}

simulation_coerce_future_branch_design_schema <- function(future_branch_schema) {
  required_fields <- c(
    "schema_contract",
    "schema_stage",
    "facet_axes",
    "assignment_axis",
    "input_keys",
    "canonical_design_variables",
    "default_design",
    "note"
  )
  attach_grid_semantics <- function(x) {
    if (is.list(x$grid_semantics) &&
        identical(x$grid_semantics$semantics_contract, "arbitrary_facet_design_grid_semantics")) {
      return(c(x[required_fields], list(grid_semantics = x$grid_semantics)))
    }

    facet_axes <- x$facet_axes
    assignment_axis <- x$assignment_axis
    assignment_var <- as.character(
      assignment_axis$current_planning_count_variable %||% "raters_per_person"
    )
    assignment_key <- as.character(
      assignment_axis$future_input_key %||% "assignment"
    )

    c(
      x[required_fields],
      list(
        grid_semantics = list(
          semantics_contract = "arbitrary_facet_design_grid_semantics",
          semantics_stage = as.character(x$schema_stage %||% "schema_only"),
          id_variable = "design_id",
          canonical_columns = c(
            "design_id",
            as.character(facet_axes$canonical_design_variable),
            assignment_var
          ),
          public_columns = c(
            "design_id",
            as.character(facet_axes$public_design_alias),
            as.character(assignment_axis$current_planning_count_alias)
          ),
          branch_columns = c(
            "design_id",
            as.character(facet_axes$input_key),
            assignment_key
          ),
          feasibility_rule = paste0(
            assignment_var,
            " <= ",
            as.character(assignment_axis$depends_on_canonical_design_variable %||% "n_rater")
          ),
          row_meaning = paste(
            "Each row is one schema-only future-branch design condition formed by",
            "crossing facet-count axes and the assignment axis after applying the",
            "current feasibility rule."
          ),
          default_id_prefix = "F"
        )
      )
    )
  }
  if (is.list(future_branch_schema) &&
      all(required_fields %in% names(future_branch_schema)) &&
      is.data.frame(future_branch_schema$facet_axes) &&
      is.list(future_branch_schema$assignment_axis) &&
      is.list(future_branch_schema$default_design)) {
    return(attach_grid_semantics(future_branch_schema))
  }

  design_schema <- future_branch_schema$design_schema %||% NULL
  if (is.list(design_schema) &&
      all(required_fields %in% names(design_schema)) &&
      is.data.frame(design_schema$facet_axes) &&
      is.list(design_schema$assignment_axis) &&
      is.list(design_schema$default_design)) {
    return(attach_grid_semantics(design_schema))
  }

  if (!is.list(future_branch_schema) ||
      !is.data.frame(future_branch_schema$facet_table) ||
      !is.list(future_branch_schema$assignment_axis) ||
      !is.list(future_branch_schema$design_template)) {
    stop("`future_branch_schema` is not a valid schema-only arbitrary-facet branch contract.",
         call. = FALSE)
  }

  facet_axes <- future_branch_schema$facet_table |>
    dplyr::transmute(
      input_key = .data$future_facet_key,
      facet = .data$facet,
      facet_kind = .data$facet_kind,
      axis_class = .data$future_axis_class,
      level_count = .data$level_count,
      canonical_design_variable = .data$current_planning_count_variable,
      public_design_alias = .data$current_planning_count_alias,
      current_planner_role_supported = .data$current_planner_role_supported,
      arbitrary_facet_branch_candidate = .data$arbitrary_facet_branch_candidate
    )
  assignment_axis <- future_branch_schema$assignment_axis

  attach_grid_semantics(list(
    schema_contract = "arbitrary_facet_design_schema",
    schema_stage = as.character(future_branch_schema$planner_stage %||% "schema_only"),
    facet_axes = facet_axes,
    assignment_axis = assignment_axis,
    input_keys = c(as.character(facet_axes$input_key), assignment_axis$future_input_key),
    canonical_design_variables = c(
      as.character(facet_axes$canonical_design_variable),
      assignment_axis$current_planning_count_variable
    ),
    default_design = future_branch_schema$design_template,
    note = paste(
      "Compatibility coercion of a schema-only future arbitrary-facet branch",
      "contract into the bundled design-schema object."
    )
  ))
}

simulation_build_future_branch_design_grid <- function(design_schema,
                                                       design = NULL,
                                                       id_prefix = "F") {
  design_schema <- simulation_coerce_future_branch_design_schema(design_schema)
  parsed <- simulation_parse_future_branch_design(
    design = design,
    future_branch_schema = list(design_schema = design_schema),
    arg_name = "design"
  )

  defaults <- design_schema$default_design %||% list()
  default_facets <- defaults$facets %||% list()
  facet_axes <- design_schema$facet_axes
  assignment_axis <- design_schema$assignment_axis
  grid_semantics <- design_schema$grid_semantics

  values <- list()
  for (i in seq_len(nrow(facet_axes))) {
    canonical <- as.character(facet_axes$canonical_design_variable[i])
    input_key <- as.character(facet_axes$input_key[i])
    value <- parsed[[canonical]] %||% default_facets[[input_key]]
    values[[canonical]] <- simulation_validate_count_values(value, canonical, min_value = 2L)
  }

  assignment_var <- as.character(assignment_axis$current_planning_count_variable %||% "raters_per_person")
  assignment_key <- as.character(assignment_axis$future_input_key %||% "assignment")
  assignment_default <- defaults[[assignment_key]] %||%
    values[[as.character(assignment_axis$depends_on_canonical_design_variable %||% "n_rater")]]
  values[[assignment_var]] <- simulation_validate_count_values(
    parsed[[assignment_var]] %||% assignment_default,
    assignment_var,
    min_value = 1L
  )

  canonical_order <- setdiff(
    as.character(grid_semantics$canonical_columns %||% c("design_id", as.character(facet_axes$canonical_design_variable), assignment_var)),
    "design_id"
  )
  design_grid <- expand.grid(
    values[canonical_order],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  dependency_var <- as.character(
    assignment_axis$depends_on_canonical_design_variable %||% "n_rater"
  )
  design_grid <- design_grid[
    design_grid[[assignment_var]] <= design_grid[[dependency_var]],
    ,
    drop = FALSE
  ]
  if (nrow(design_grid) == 0L) {
    stop(
      "No valid future-branch design rows remain after enforcing `",
      assignment_var,
      " <= ",
      dependency_var,
      "`.",
      call. = FALSE
    )
  }

  prefix <- as.character(grid_semantics$default_id_prefix %||% "F")
  design_grid$design_id <- sprintf("%s%02d", as.character(id_prefix[1] %||% prefix), seq_len(nrow(design_grid)))
  design_grid <- tibble::as_tibble(design_grid[, c("design_id", canonical_order), drop = FALSE])

  branch_grid <- tibble::tibble(design_id = design_grid$design_id)
  for (i in seq_len(nrow(facet_axes))) {
    branch_grid[[as.character(facet_axes$input_key[i])]] <-
      design_grid[[as.character(facet_axes$canonical_design_variable[i])]]
  }
  branch_grid[[assignment_key]] <- design_grid[[assignment_var]]

  aliases <- stats::setNames(
    c(as.character(facet_axes$public_design_alias), as.character(assignment_axis$current_planning_count_alias)),
    canonical_order
  )

  list(
    canonical = design_grid,
    public = simulation_append_design_alias_columns(design_grid, aliases),
    branch = branch_grid,
    design_schema = design_schema,
    grid_semantics = grid_semantics
  )
}

simulation_future_branch_preview <- function(future_branch_schema,
                                             design = NULL,
                                             id_prefix = "F") {
  grid_contract <- simulation_future_branch_grid_contract(
    future_branch_schema = future_branch_schema,
    design = design,
    id_prefix = id_prefix
  )
  preview_fields <- c(
    "preview_available",
    "reason",
    "default_design",
    "canonical",
    "public",
    "branch",
    "design_schema",
    "grid_semantics"
  )

  names(grid_contract)[names(grid_contract) == "contract"] <- "preview_contract"
  names(grid_contract)[names(grid_contract) == "stage"] <- "preview_stage"
  grid_contract[c("preview_contract", "preview_stage", preview_fields)]
}

simulation_future_branch_grid_contract <- function(future_branch_schema,
                                                   design = NULL,
                                                   id_prefix = "F") {
  design_schema <- simulation_coerce_future_branch_design_schema(future_branch_schema)
  default_design <- design_schema$default_design %||% list()
  grid_semantics <- design_schema$grid_semantics
  facet_defaults <- default_design$facets %||% list()
  assignment_key <- as.character(
    design_schema$assignment_axis$future_input_key %||% "assignment"
  )
  assignment_default <- default_design[[assignment_key]] %||% NULL

  facet_complete <- length(facet_defaults) == nrow(design_schema$facet_axes) &&
    all(vapply(
      facet_defaults,
      function(x) length(x) == 1L && is.numeric(x) && is.finite(x),
      logical(1)
    ))
  assignment_complete <- length(assignment_default) == 1L &&
    is.numeric(assignment_default) && is.finite(assignment_default)

  preview_design <- design %||% default_design
  if (!facet_complete || !assignment_complete) {
    return(list(
      contract = "arbitrary_facet_design_grid_contract",
      stage = as.character(design_schema$schema_stage %||% "schema_only"),
      planner_contract = as.character(
        future_branch_schema$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      input_contract = as.character(
        future_branch_schema$input_contract %||% "design$facets(named counts)"
      ),
      preview_available = FALSE,
      reason = paste(
        "The schema-only future arbitrary-facet branch does not yet have",
        "complete default facet counts and assignment counts, so no",
        "schema-only branch grid can be materialized."
      ),
      default_design = default_design,
      design_schema = design_schema,
      grid_semantics = grid_semantics,
      note = paste(
        "Schema-only future-branch grid contract is unavailable until the",
        "default nested design carries finite facet and assignment counts."
      )
    ))
  }

  preview_grid <- simulation_build_future_branch_design_grid(
    design_schema = design_schema,
    design = preview_design,
    id_prefix = id_prefix
  )

  list(
    contract = "arbitrary_facet_design_grid_contract",
    stage = as.character(design_schema$schema_stage %||% "schema_only"),
    planner_contract = as.character(
      future_branch_schema$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    input_contract = as.character(
      future_branch_schema$input_contract %||% "design$facets(named counts)"
    ),
    preview_available = TRUE,
    reason = paste(
      "Schema-only branch grid built from the future arbitrary-facet",
      "branch design schema and its current default design."
    ),
    default_design = default_design,
    canonical = preview_grid$canonical,
    public = preview_grid$public,
    branch = preview_grid$branch,
    design_schema = design_schema,
    grid_semantics = preview_grid$grid_semantics,
    note = paste(
      "Schema-only future-branch grid contract bundling the default branch",
      "grid, design schema, and grid semantics without activating arbitrary-",
      "facet planner logic."
    )
  )
}

simulation_coerce_future_branch_grid_contract <- function(x) {
  required_fields <- c(
    "contract",
    "stage",
    "planner_contract",
    "input_contract",
    "preview_available",
    "reason",
    "default_design",
    "design_schema",
    "grid_semantics",
    "note"
  )
  optional_grid_fields <- c("canonical", "public", "branch")

  if (is.list(x) &&
      identical(x$contract %||% "", "arbitrary_facet_design_grid_contract") &&
      all(required_fields %in% names(x)) &&
      is.list(x$design_schema) &&
      is.list(x$grid_semantics)) {
    out <- x[unique(c(required_fields, optional_grid_fields[optional_grid_fields %in% names(x)]))]
    return(out)
  }

  nested <- x$grid_contract %||% x$future_branch_grid_contract %||% NULL
  if (is.list(nested)) {
    return(simulation_coerce_future_branch_grid_contract(nested))
  }

  branch <- x$future_branch_schema %||% x
  if (!is.list(branch)) {
    stop("`x` is not a valid schema-only future-branch grid contract.", call. = FALSE)
  }

  simulation_future_branch_grid_contract(branch)
}

simulation_build_future_branch_design_grid_from_contract <- function(grid_contract,
                                                                     design = NULL,
                                                                     id_prefix = NULL) {
  grid_contract <- simulation_coerce_future_branch_grid_contract(grid_contract)
  prefix <- as.character(id_prefix %||% grid_contract$grid_semantics$default_id_prefix %||% "F")

  if (is.null(design) &&
      isTRUE(grid_contract$preview_available) &&
      all(c("canonical", "public", "branch") %in% names(grid_contract))) {
    return(list(
      canonical = grid_contract$canonical,
      public = grid_contract$public,
      branch = grid_contract$branch,
      design_schema = grid_contract$design_schema,
      grid_semantics = grid_contract$grid_semantics,
      grid_contract = grid_contract
    ))
  }

  built <- simulation_build_future_branch_design_grid(
    design_schema = grid_contract$design_schema,
    design = design %||% grid_contract$default_design,
    id_prefix = prefix
  )
  built$grid_contract <- grid_contract
  built
}

simulation_materialize_future_branch_grid <- function(x,
                                                      design = NULL,
                                                      id_prefix = NULL) {
  grid_contract <- NULL

  if (is.list(x) &&
      identical(x$contract %||% "", "arbitrary_facet_design_grid_contract")) {
    grid_contract <- x
  }

  if (is.null(grid_contract) && is.list(x)) {
    nested <- x$future_branch_grid_contract %||%
      x$grid_contract %||%
      x$future_branch_schema %||%
      x$planning_schema %||%
      x$settings$planning_schema %||%
      x$sim_spec$planning_schema %||%
      x$ademp$data_generating_mechanism$planning_schema %||%
      NULL
    if (is.list(nested)) {
      grid_contract <- simulation_coerce_future_branch_grid_contract(nested)
    }
  }

  if (is.null(grid_contract)) {
    schema <- tryCatch(
      simulation_object_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$future_branch_grid_contract)) {
      grid_contract <- schema$future_branch_grid_contract
    }
  }

  if (is.null(grid_contract)) {
    schema <- tryCatch(
      simulation_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$future_branch_grid_contract)) {
      grid_contract <- schema$future_branch_grid_contract
    }
  }

  if (is.null(grid_contract)) {
    stop(
      "`x` is not a valid planning object, planning schema, sim spec, or ",
      "future-branch grid contract.",
      call. = FALSE
    )
  }

  simulation_build_future_branch_design_grid_from_contract(
    grid_contract = grid_contract,
    design = design,
    id_prefix = id_prefix
  )
}

simulation_future_branch_grid_bundle <- function(x,
                                                 design = NULL,
                                                 id_prefix = NULL) {
  grid_contract <- simulation_coerce_future_branch_grid_contract(x)

  if (is.null(design) && !isTRUE(grid_contract$preview_available)) {
    return(list(
      bundle_contract = "arbitrary_facet_design_grid_bundle",
      bundle_stage = as.character(grid_contract$stage %||% "schema_only"),
      planner_contract = as.character(
        grid_contract$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      input_contract = as.character(
        grid_contract$input_contract %||% "design$facets(named counts)"
      ),
      grid_available = FALSE,
      reason = as.character(
        grid_contract$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      default_design = grid_contract$default_design,
      design_schema = grid_contract$design_schema,
      grid_semantics = grid_contract$grid_semantics,
      grid_contract = grid_contract,
      note = paste(
        "Schema-only future-branch bundle carrying the branch grid contract,",
        "schema, and semantics without a materialized grid because default",
        "counts are not yet available."
      )
    ))
  }

  built <- simulation_build_future_branch_design_grid_from_contract(
    grid_contract = grid_contract,
    design = design,
    id_prefix = id_prefix
  )

  list(
    bundle_contract = "arbitrary_facet_design_grid_bundle",
    bundle_stage = as.character(grid_contract$stage %||% "schema_only"),
    planner_contract = as.character(
      grid_contract$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    input_contract = as.character(
      grid_contract$input_contract %||% "design$facets(named counts)"
    ),
    grid_available = TRUE,
    reason = paste(
      "Schema-only future-branch bundle materialized from the branch grid",
      "contract and matching design schema."
    ),
    canonical = built$canonical,
    public = built$public,
    branch = built$branch,
    default_design = grid_contract$default_design,
    design_schema = built$design_schema,
    grid_semantics = built$grid_semantics,
    grid_contract = grid_contract,
    note = paste(
      "Schema-only future-branch bundle that materializes canonical, public,",
      "and branch-facing design grids from one authoritative branch-grid",
      "contract."
    )
  )
}

simulation_coerce_future_branch_grid_bundle <- function(x) {
  required_fields <- c(
    "bundle_contract",
    "bundle_stage",
    "planner_contract",
    "input_contract",
    "grid_available",
    "reason",
    "default_design",
    "design_schema",
    "grid_semantics",
    "grid_contract",
    "note"
  )
  optional_grid_fields <- c("canonical", "public", "branch")

  if (is.list(x) &&
      identical(x$bundle_contract %||% "", "arbitrary_facet_design_grid_bundle") &&
      all(required_fields %in% names(x)) &&
      is.list(x$design_schema) &&
      is.list(x$grid_semantics) &&
      is.list(x$grid_contract)) {
    return(x[unique(c(required_fields, optional_grid_fields[optional_grid_fields %in% names(x)]))])
  }

  nested <- x$grid_bundle %||% x$future_branch_grid_bundle %||% NULL
  if (is.list(nested)) {
    return(simulation_coerce_future_branch_grid_bundle(nested))
  }

  simulation_future_branch_grid_bundle(x)
}

simulation_materialize_future_branch_grid_bundle <- function(x,
                                                             design = NULL,
                                                             id_prefix = NULL) {
  grid_bundle <- NULL

  if (is.list(x) &&
      identical(x$bundle_contract %||% "", "arbitrary_facet_design_grid_bundle")) {
    grid_bundle <- x
  }

  if (is.null(grid_bundle) && is.list(x)) {
    nested <- x$future_branch_grid_bundle %||%
      x$grid_bundle %||%
      x$future_branch_schema %||%
      x$planning_schema %||%
      x$settings$planning_schema %||%
      x$sim_spec$planning_schema %||%
      x$ademp$data_generating_mechanism$planning_schema %||%
      NULL
    if (is.list(nested)) {
      grid_bundle <- simulation_coerce_future_branch_grid_bundle(nested)
    }
  }

  if (is.null(grid_bundle)) {
    schema <- tryCatch(
      simulation_object_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$future_branch_grid_bundle)) {
      grid_bundle <- schema$future_branch_grid_bundle
    }
  }

  if (is.null(grid_bundle)) {
    schema <- tryCatch(
      simulation_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$future_branch_grid_bundle)) {
      grid_bundle <- schema$future_branch_grid_bundle
    }
  }

  if (is.null(grid_bundle)) {
    stop(
      "`x` is not a valid planning object, planning schema, sim spec, or ",
      "future-branch grid bundle.",
      call. = FALSE
    )
  }

  if (is.null(design) && isTRUE(grid_bundle$grid_available)) {
    return(grid_bundle)
  }

  simulation_future_branch_grid_bundle(
    x = grid_bundle$grid_contract %||% grid_bundle,
    design = design,
    id_prefix = id_prefix
  )
}

simulation_future_branch_grid_view <- function(x,
                                               view = c("canonical", "public", "branch"),
                                               design = NULL,
                                               id_prefix = NULL) {
  view <- match.arg(view)
  grid_bundle <- simulation_materialize_future_branch_grid_bundle(
    x = x,
    design = design,
    id_prefix = id_prefix
  )

  if (!isTRUE(grid_bundle$grid_available) || !is.data.frame(grid_bundle[[view]])) {
    stop(
      "The requested `", view, "` future-branch grid view is not currently ",
      "available from this schema-only bundle.",
      call. = FALSE
    )
  }

  grid_bundle[[view]]
}

simulation_future_branch_grid_context <- function(x,
                                                  design = NULL,
                                                  id_prefix = NULL) {
  grid_bundle <- simulation_materialize_future_branch_grid_bundle(
    x = x,
    design = design,
    id_prefix = id_prefix
  )

  design_schema <- grid_bundle$design_schema
  grid_semantics <- grid_bundle$grid_semantics
  facet_axes <- design_schema$facet_axes
  assignment_axis <- design_schema$assignment_axis

  axis_table <- dplyr::bind_rows(
    facet_axes |>
      dplyr::transmute(
        axis_source = "facet",
        input_key = .data$input_key,
        canonical_design_variable = .data$canonical_design_variable,
        public_design_alias = .data$public_design_alias,
        axis_class = .data$axis_class,
        facet = .data$facet
      ),
    tibble::tibble(
      axis_source = "assignment",
      input_key = as.character(assignment_axis$future_input_key %||% "assignment"),
      canonical_design_variable = as.character(
        assignment_axis$current_planning_count_variable %||% "raters_per_person"
      ),
      public_design_alias = as.character(
        assignment_axis$current_planning_count_alias %||% "raters_per_person"
      ),
      axis_class = as.character(assignment_axis$axis_class %||% "assignment_count"),
      facet = NA_character_
    )
  )

  if (!isTRUE(grid_bundle$grid_available) || !is.data.frame(grid_bundle$canonical)) {
    axis_table$values <- vector("list", nrow(axis_table))
    axis_table$n_values <- rep(NA_integer_, nrow(axis_table))
    axis_table$varying <- rep(NA, nrow(axis_table))
    axis_table$fixed_value <- rep(NA_integer_, nrow(axis_table))

    return(list(
      context_contract = "arbitrary_facet_design_grid_context",
      context_stage = as.character(grid_bundle$bundle_stage %||% "schema_only"),
      grid_available = FALSE,
      reason = as.character(grid_bundle$reason %||% "No future-branch grid is available."),
      id_variable = as.character(grid_semantics$id_variable %||% "design_id"),
      canonical_columns = as.character(grid_semantics$canonical_columns %||% character(0)),
      public_columns = as.character(grid_semantics$public_columns %||% character(0)),
      branch_columns = as.character(grid_semantics$branch_columns %||% character(0)),
      axis_table = axis_table,
      varying_canonical = character(0),
      fixed_canonical = character(0),
      varying_input_keys = character(0),
      fixed_input_keys = character(0),
      grid_bundle = grid_bundle,
      note = paste(
        "Schema-only future-branch context exposes axis metadata even when a",
        "materialized grid is not yet available."
      )
    ))
  }

  canonical_grid <- grid_bundle$canonical
  axis_values <- lapply(axis_table$canonical_design_variable, function(var) {
    sort(unique(as.integer(canonical_grid[[var]])))
  })
  n_values <- vapply(axis_values, length, integer(1))
  varying <- n_values > 1L
  fixed_value <- vapply(axis_values, function(vals) {
    if (length(vals) == 1L) vals[[1]] else NA_integer_
  }, integer(1))

  axis_table$values <- axis_values
  axis_table$n_values <- as.integer(n_values)
  axis_table$varying <- as.logical(varying)
  axis_table$fixed_value <- as.integer(fixed_value)

  list(
    context_contract = "arbitrary_facet_design_grid_context",
    context_stage = as.character(grid_bundle$bundle_stage %||% "schema_only"),
    grid_available = TRUE,
    reason = as.character(
      grid_bundle$reason %||% "Future-branch grid context built from a materialized bundle."
    ),
    id_variable = as.character(grid_semantics$id_variable %||% "design_id"),
    canonical_columns = as.character(grid_semantics$canonical_columns %||% names(canonical_grid)),
    public_columns = as.character(grid_semantics$public_columns %||% names(grid_bundle$public %||% canonical_grid)),
    branch_columns = as.character(grid_semantics$branch_columns %||% names(grid_bundle$branch %||% canonical_grid)),
    axis_table = axis_table,
    varying_canonical = as.character(axis_table$canonical_design_variable[axis_table$varying]),
    fixed_canonical = as.character(axis_table$canonical_design_variable[!axis_table$varying]),
    varying_input_keys = as.character(axis_table$input_key[axis_table$varying]),
    fixed_input_keys = as.character(axis_table$input_key[!axis_table$varying]),
    grid_bundle = grid_bundle,
    note = paste(
      "Schema-only future-branch context summarizing which branch axes vary",
      "within the currently materialized design grid."
    )
  )
}

simulation_future_branch_grid_summary <- function(x,
                                                  design = NULL,
                                                  id_prefix = NULL) {
  grid_context <- simulation_future_branch_grid_context(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  axis_table <- tibble::as_tibble(grid_context$axis_table)
  varying_rows <- !is.na(axis_table$varying) & as.logical(axis_table$varying)
  fixed_rows <- !is.na(axis_table$varying) & !as.logical(axis_table$varying)
  fixed_values <- stats::setNames(
    as.list(axis_table$fixed_value[fixed_rows]),
    as.character(axis_table$canonical_design_variable[fixed_rows])
  )

  if (!isTRUE(grid_context$grid_available)) {
    return(list(
      summary_contract = "arbitrary_facet_design_grid_summary",
      summary_stage = as.character(grid_context$context_stage %||% "schema_only"),
      planner_contract = as.character(
        grid_context$grid_bundle$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      input_contract = as.character(
        grid_context$grid_bundle$input_contract %||% "design$facets(named counts)"
      ),
      grid_available = FALSE,
      reason = as.character(
        grid_context$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      n_designs = 0L,
      n_varying_axes = sum(varying_rows),
      n_fixed_axes = sum(fixed_rows),
      varying_canonical = as.character(grid_context$varying_canonical),
      fixed_canonical = as.character(grid_context$fixed_canonical),
      varying_input_keys = as.character(grid_context$varying_input_keys),
      fixed_input_keys = as.character(grid_context$fixed_input_keys),
      fixed_values = fixed_values,
      axis_table = axis_table,
      grid_context = grid_context,
      grid_bundle = grid_context$grid_bundle,
      note = paste(
        "Schema-only future-branch summary is unavailable because no",
        "materialized branch grid exists yet, but axis metadata are still",
        "returned for branch-side planning helpers."
      )
    ))
  }

  canonical <- tibble::as_tibble(grid_context$grid_bundle$canonical)
  public <- tibble::as_tibble(grid_context$grid_bundle$public)
  branch <- tibble::as_tibble(grid_context$grid_bundle$branch)

  list(
    summary_contract = "arbitrary_facet_design_grid_summary",
    summary_stage = as.character(grid_context$context_stage %||% "schema_only"),
    planner_contract = as.character(
      grid_context$grid_bundle$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    input_contract = as.character(
      grid_context$grid_bundle$input_contract %||% "design$facets(named counts)"
    ),
    grid_available = TRUE,
    reason = as.character(
      grid_context$reason %||%
        "Schema-only future-branch grid is materialized."
    ),
    n_designs = nrow(canonical),
    n_varying_axes = sum(varying_rows),
    n_fixed_axes = sum(fixed_rows),
    varying_canonical = as.character(grid_context$varying_canonical),
    fixed_canonical = as.character(grid_context$fixed_canonical),
    varying_input_keys = as.character(grid_context$varying_input_keys),
    fixed_input_keys = as.character(grid_context$fixed_input_keys),
    fixed_values = fixed_values,
    axis_table = axis_table,
    canonical = canonical,
    public = public,
    branch = branch,
    grid_context = grid_context,
    grid_bundle = grid_context$grid_bundle,
    note = paste(
      "Schema-only future-branch summary bundling the materialized grid with",
      "axis-level varying/fixed metadata for later arbitrary-facet branch",
      "helpers."
    )
  )
}

simulation_future_branch_axis_lookup <- function(axis_table) {
  axis_table <- tibble::as_tibble(axis_table)
  lookup <- c(
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$canonical_design_variable)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$public_design_alias)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$input_key)
    )
  )
  lookup <- lookup[!is.na(names(lookup)) & nzchar(names(lookup))]
  lookup[!duplicated(names(lookup))]
}

simulation_resolve_future_branch_axis <- function(value,
                                                  axis_table,
                                                  arg_name = "axis",
                                                  allow_null = FALSE) {
  if (allow_null && is.null(value)) {
    return(NULL)
  }
  lookup <- simulation_future_branch_axis_lookup(axis_table)
  value <- as.character(value[1] %||% "")
  resolved <- unname(lookup[[value]])
  if (!is.null(resolved) && nzchar(resolved)) {
    return(resolved)
  }
  stop(
    "`", arg_name, "` must be one of: ",
    paste(sort(unique(names(lookup))), collapse = ", "),
    ".",
    call. = FALSE
  )
}

simulation_future_branch_axis_label <- function(canonical,
                                                axis_table,
                                                view = c("public", "canonical", "branch")) {
  view <- match.arg(view)
  axis_table <- tibble::as_tibble(axis_table)
  row <- axis_table[axis_table$canonical_design_variable == canonical, , drop = FALSE]
  if (nrow(row) == 0L) {
    return(as.character(canonical))
  }
  if (identical(view, "canonical")) {
    return(as.character(canonical))
  }
  if (identical(view, "branch")) {
    label <- as.character(row$input_key[[1]])
    if (nzchar(label)) return(label)
  }
  label <- as.character(row$public_design_alias[[1]])
  if (nzchar(label)) return(label)
  as.character(canonical)
}

simulation_future_branch_grid_recommendation <- function(x,
                                                         design = NULL,
                                                         prefer = NULL,
                                                         id_prefix = NULL) {
  grid_summary <- simulation_future_branch_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )

  axis_table <- tibble::as_tibble(grid_summary$axis_table)
  lookup <- c(
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$canonical_design_variable)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$public_design_alias)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$input_key)
    )
  )
  lookup <- lookup[!duplicated(names(lookup))]

  if (!isTRUE(grid_summary$grid_available)) {
    return(list(
      recommendation_contract = "arbitrary_facet_design_grid_recommendation",
      recommendation_stage = as.character(grid_summary$summary_stage %||% "schema_only"),
      planner_contract = as.character(
        grid_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      recommendation_available = FALSE,
      reason = as.character(
        grid_summary$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      selection_rule = paste(
        "Unavailable because no schema-only branch grid exists yet.",
        "This helper does not fabricate facet counts."
      ),
      prefer = character(0),
      rank_order = character(0),
      grid_summary = grid_summary,
      note = paste(
        "Schema-only future-branch recommendation is unavailable until the",
        "underlying branch grid is materialized."
      )
    ))
  }

  canonical_order <- setdiff(
    as.character(
      grid_summary$grid_bundle$grid_semantics$canonical_columns %||%
        names(grid_summary$canonical)
    ),
    "design_id"
  )
  default_prefer <- if (length(grid_summary$varying_canonical) > 0L) {
    intersect(canonical_order, grid_summary$varying_canonical)
  } else {
    canonical_order
  }

  if (is.null(prefer)) {
    resolved_prefer <- default_prefer
  } else {
    prefer <- unique(as.character(prefer))
    resolved_prefer <- unname(lookup[prefer])
    resolved_prefer <- unique(resolved_prefer[!is.na(resolved_prefer) & nzchar(resolved_prefer)])
    if (length(resolved_prefer) == 0L) {
      stop(
        "`prefer` must resolve to at least one future-branch design axis. Valid names: ",
        paste(sort(unique(names(lookup))), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  rank_order <- unique(c(resolved_prefer, setdiff(canonical_order, resolved_prefer)))

  ranked <- dplyr::arrange(grid_summary$canonical, !!!rlang::syms(rank_order))
  recommended_canonical <- dplyr::slice_head(ranked, n = 1)
  recommended_id <- as.character(recommended_canonical$design_id[[1]])
  recommended_public <- dplyr::filter(grid_summary$public, .data$design_id == recommended_id)
  recommended_branch <- dplyr::filter(grid_summary$branch, .data$design_id == recommended_id)

  list(
    recommendation_contract = "arbitrary_facet_design_grid_recommendation",
    recommendation_stage = as.character(grid_summary$summary_stage %||% "schema_only"),
    planner_contract = as.character(
      grid_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    recommendation_available = TRUE,
    reason = paste(
      "Schema-only future-branch baseline design selected by deterministic",
      "lexicographic ordering over the requested branch axes."
    ),
    selection_rule = paste(
      "Pick the lexicographically smallest feasible design row after sorting",
      "the schema-only canonical grid by `prefer`, then by the remaining",
      "canonical design axes. This is a deterministic baseline pick, not a",
      "performance-based recommendation."
    ),
    prefer = resolved_prefer,
    rank_order = rank_order,
    recommended_design_id = recommended_id,
    recommended_canonical = recommended_canonical,
    recommended_public = recommended_public,
    recommended_branch = recommended_branch,
    grid_summary = grid_summary,
    note = paste(
      "Schema-only future-branch baseline pick derived from the currently",
      "materialized branch grid without activating arbitrary-facet planner",
      "logic or performance criteria."
    )
  )
}

simulation_future_branch_grid_table <- function(x,
                                                design = NULL,
                                                prefer = NULL,
                                                view = c("public", "canonical", "branch"),
                                                id_prefix = NULL) {
  view <- match.arg(view)
  grid_summary <- simulation_future_branch_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  recommendation <- simulation_future_branch_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )

  if (!isTRUE(grid_summary$grid_available)) {
    return(list(
      table_contract = "arbitrary_facet_design_grid_table",
      table_stage = as.character(grid_summary$summary_stage %||% "schema_only"),
      planner_contract = as.character(
        grid_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      grid_available = FALSE,
      reason = as.character(
        grid_summary$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      view = view,
      view_label = view,
      table = tibble::tibble(),
      recommended_design_id = character(0),
      grid_summary = grid_summary,
      grid_recommendation = recommendation,
      note = paste(
        "Schema-only future-branch table is unavailable until the",
        "underlying branch grid is materialized."
      )
    ))
  }

  table <- switch(
    view,
    canonical = tibble::as_tibble(grid_summary$canonical),
    public = tibble::as_tibble(grid_summary$public),
    branch = tibble::as_tibble(grid_summary$branch)
  )
  table$recommended <- table$design_id == recommendation$recommended_design_id

  list(
    table_contract = "arbitrary_facet_design_grid_table",
    table_stage = as.character(grid_summary$summary_stage %||% "schema_only"),
    planner_contract = as.character(
      grid_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    grid_available = TRUE,
    reason = as.character(
      grid_summary$reason %||%
        "Schema-only future-branch grid is materialized."
    ),
    view = view,
    view_label = view,
    table = table,
    recommended_design_id = recommendation$recommended_design_id,
    grid_summary = grid_summary,
    grid_recommendation = recommendation,
    note = paste(
      "Schema-only future-branch table exposing one grid view together with",
      "the deterministic baseline design flag."
    )
  )
}

simulation_future_branch_grid_plot_payload <- function(x,
                                                       design = NULL,
                                                       x_var = NULL,
                                                       group_var = NULL,
                                                       prefer = NULL,
                                                       view = c("public", "canonical", "branch"),
                                                       id_prefix = NULL) {
  view <- match.arg(view)
  grid_summary <- simulation_future_branch_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  recommendation <- simulation_future_branch_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )
  axis_table <- tibble::as_tibble(grid_summary$axis_table)

  if (!isTRUE(grid_summary$grid_available)) {
    return(new_mfrm_plot_data(
      "future_branch_grid_schema",
      list(
        title = "Schema-Only Future-Branch Grid",
        subtitle = as.character(
          grid_summary$reason %||%
            "No schema-only branch grid is currently materialized."
        ),
        plot_available = FALSE,
        reason = as.character(
          grid_summary$reason %||%
            "No schema-only branch grid is currently materialized."
        ),
        view = view,
        x_var = NULL,
        x_label = NULL,
        group_var = NULL,
        group_label = NULL,
        recommended_design_id = character(0),
        data = tibble::tibble(),
        axis_table = axis_table,
        grid_summary = grid_summary,
        grid_recommendation = recommendation,
        note = paste(
          "Draw-free plotting payload for the schema-only future branch is",
          "unavailable until the branch grid can be materialized."
        )
      )
    ))
  }

  lookup <- simulation_future_branch_axis_lookup(axis_table)
  varying <- as.character(grid_summary$varying_canonical)
  canonical_default_x <- if (length(varying) > 0L) varying[[1]] else as.character(axis_table$canonical_design_variable[[1]])
  x_canonical <- if (is.null(x_var)) {
    canonical_default_x
  } else {
    simulation_resolve_future_branch_axis(x_var, axis_table, arg_name = "x_var")
  }

  group_default <- setdiff(varying, x_canonical)
  group_canonical <- if (is.null(group_var)) {
    if (length(group_default) > 0L) group_default[[1]] else NULL
  } else {
    simulation_resolve_future_branch_axis(group_var, axis_table, arg_name = "group_var", allow_null = TRUE)
  }
  if (!is.null(group_canonical) && identical(group_canonical, x_canonical)) {
    stop("`group_var` must differ from `x_var`.", call. = FALSE)
  }

  display_table <- switch(
    view,
    canonical = tibble::as_tibble(grid_summary$canonical),
    public = tibble::as_tibble(grid_summary$public),
    branch = tibble::as_tibble(grid_summary$branch)
  )
  canonical_axes <- tibble::as_tibble(grid_summary$canonical[, c("design_id", unique(c(x_canonical, group_canonical))), drop = FALSE])
  names(canonical_axes)[-1] <- paste0(names(canonical_axes)[-1], "__canonical")
  plot_data <- dplyr::left_join(display_table, canonical_axes, by = "design_id")
  plot_data$x_value <- plot_data[[paste0(x_canonical, "__canonical")]]
  if (!is.null(group_canonical)) {
    plot_data$group_value <- as.character(plot_data[[paste0(group_canonical, "__canonical")]])
  } else {
    plot_data$group_value <- "All designs"
  }
  plot_data$recommended <- plot_data$design_id == recommendation$recommended_design_id
  plot_data <- dplyr::arrange(plot_data, .data$x_value, .data$group_value, dplyr::desc(.data$recommended))

  x_label <- simulation_future_branch_axis_label(x_canonical, axis_table, view = view)
  group_label <- if (is.null(group_canonical)) "Design set" else simulation_future_branch_axis_label(group_canonical, axis_table, view = view)
  fixed_rows <- !is.na(axis_table$varying) & !as.logical(axis_table$varying)
  fixed_text <- if (any(fixed_rows)) {
    fixed_bits <- vapply(which(fixed_rows), function(i) {
      paste0(
        simulation_future_branch_axis_label(axis_table$canonical_design_variable[[i]], axis_table, view = view),
        "=",
        axis_table$fixed_value[[i]]
      )
    }, character(1))
    paste(fixed_bits, collapse = ", ")
  } else {
    "No fixed design axes"
  }

  legend <- if (!is.null(group_canonical)) {
    group_levels <- unique(as.character(plot_data$group_value))
    new_plot_legend(
      label = group_levels,
      role = rep("group", length(group_levels)),
      aesthetic = rep("line", length(group_levels)),
      value = group_levels
    )
  } else {
    new_plot_legend(
      label = "All designs",
      role = "group",
      aesthetic = "line",
      value = "All designs"
    )
  }

  new_mfrm_plot_data(
    "future_branch_grid_schema",
    list(
      title = "Schema-Only Future-Branch Grid",
      subtitle = fixed_text,
      plot_available = TRUE,
      reason = as.character(
        grid_summary$reason %||%
          "Schema-only future-branch grid is materialized."
      ),
      view = view,
      x_var = x_canonical,
      x_label = x_label,
      group_var = group_canonical,
      group_label = group_label,
      recommended_design_id = recommendation$recommended_design_id,
      selection_rule = recommendation$selection_rule,
      data = plot_data,
      axis_table = axis_table,
      grid_summary = grid_summary,
      grid_recommendation = recommendation,
      legend = legend,
      note = paste(
        "Draw-free plotting payload for the schema-only future branch. It",
        "summarizes the current branch grid and the deterministic baseline",
        "pick without implying an active arbitrary-facet planner."
      )
    )
  )
}

simulation_future_branch_cached_default <- function(x,
                                                    field,
                                                    design = NULL,
                                                    prefer = NULL,
                                                    view = NULL,
                                                    view_default = NULL,
                                                    mode = NULL,
                                                    mode_default = NULL,
                                                    surface = NULL,
                                                    surface_default = NULL,
                                                    table_component = NULL,
                                                    table_component_default = NULL,
                                                    component = NULL,
                                                    component_default = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  if (!is.list(x) || !is.list(x[[field]])) {
    return(NULL)
  }
  if (!is.null(design) || !is.null(prefer) || !is.null(x_var) ||
      !is.null(group_var) || !is.null(id_prefix)) {
    return(NULL)
  }
  if (!is.null(view_default) &&
      !identical(as.character(view[[1]] %||% view_default), as.character(view_default))) {
    return(NULL)
  }
  if (!is.null(mode_default) &&
      !identical(as.character(mode[[1]] %||% mode_default), as.character(mode_default))) {
    return(NULL)
  }
  if (!is.null(surface_default) &&
      !identical(as.character(surface[[1]] %||% surface_default), as.character(surface_default))) {
    return(NULL)
  }
  if (!is.null(table_component_default) &&
      !identical(
        as.character(table_component[[1]] %||% table_component_default),
        as.character(table_component_default)
      )) {
    return(NULL)
  }
  if (!is.null(component_default) &&
      !identical(as.character(component[[1]] %||% component_default), as.character(component_default))) {
    return(NULL)
  }
  x[[field]]
}

simulation_future_branch_report_bundle <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_bundle",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  grid_summary <- simulation_future_branch_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  recommendation <- simulation_future_branch_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )
  tables <- list(
    canonical = simulation_future_branch_grid_table(
      x = x,
      design = design,
      prefer = prefer,
      view = "canonical",
      id_prefix = id_prefix
    ),
    public = simulation_future_branch_grid_table(
      x = x,
      design = design,
      prefer = prefer,
      view = "public",
      id_prefix = id_prefix
    ),
    branch = simulation_future_branch_grid_table(
      x = x,
      design = design,
      prefer = prefer,
      view = "branch",
      id_prefix = id_prefix
    )
  )
  plots <- list(
    canonical = simulation_future_branch_grid_plot_payload(
      x = x,
      design = design,
      x_var = x_var,
      group_var = group_var,
      prefer = prefer,
      view = "canonical",
      id_prefix = id_prefix
    ),
    public = simulation_future_branch_grid_plot_payload(
      x = x,
      design = design,
      x_var = x_var,
      group_var = group_var,
      prefer = prefer,
      view = "public",
      id_prefix = id_prefix
    ),
    branch = simulation_future_branch_grid_plot_payload(
      x = x,
      design = design,
      x_var = x_var,
      group_var = group_var,
      prefer = prefer,
      view = "branch",
      id_prefix = id_prefix
    )
  )

  overview_table <- tibble::as_tibble(grid_summary$axis_table)

  if (!isTRUE(grid_summary$grid_available)) {
    return(list(
      report_contract = "arbitrary_facet_design_grid_report_bundle",
      report_stage = as.character(grid_summary$summary_stage %||% "schema_only"),
      planner_contract = as.character(
        grid_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      report_available = FALSE,
      reason = as.character(
        grid_summary$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      overview_table = overview_table,
      grid_summary = grid_summary,
      grid_recommendation = recommendation,
      tables = tables,
      plots = plots,
      note = paste(
        "Schema-only future-branch report bundle is unavailable until the",
        "underlying branch grid is materialized. The contract still bundles",
        "summary metadata, branch-side table contracts, and draw-free plot",
        "payloads in one internal object."
      )
    ))
  }

  list(
    report_contract = "arbitrary_facet_design_grid_report_bundle",
    report_stage = as.character(grid_summary$summary_stage %||% "schema_only"),
    planner_contract = as.character(
      grid_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    report_available = TRUE,
    reason = as.character(
      grid_summary$reason %||%
        "Schema-only future-branch grid is materialized."
    ),
    recommended_design_id = as.character(recommendation$recommended_design_id %||% character(0)),
    overview_table = overview_table,
    grid_summary = grid_summary,
    grid_recommendation = recommendation,
    tables = tables,
    plots = plots,
    note = paste(
      "Schema-only future-branch report bundle combining the internal",
      "summary, deterministic baseline recommendation, canonical/public/",
      "branch table views, and draw-free plot payloads without implying an",
      "active arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_summary <- function(x,
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_summary",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  report_bundle <- simulation_future_branch_report_bundle(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  table_views <- names(report_bundle$tables %||% list())
  plot_views <- names(report_bundle$plots %||% list())
  table_available <- vapply(
    report_bundle$tables %||% list(),
    function(obj) isTRUE(obj$grid_available),
    logical(1)
  )
  plot_available <- vapply(
    report_bundle$plots %||% list(),
    function(obj) isTRUE(obj$data$plot_available),
    logical(1)
  )
  first_or_na <- function(value) {
    if (length(value) == 0L) {
      return(NA_character_)
    }
    as.character(value[[1]] %||% NA_character_)
  }

  table_index <- tibble::tibble(
    component_type = "table",
    view = table_views,
    available = as.logical(table_available),
    n_rows = vapply(
      report_bundle$tables %||% list(),
      function(obj) if (is.data.frame(obj$table)) nrow(obj$table) else 0L,
      integer(1)
    ),
    x_var = NA_character_,
    group_var = NA_character_,
    recommended_design_id = vapply(
      report_bundle$tables %||% list(),
      function(obj) first_or_na(obj$recommended_design_id),
      character(1)
    )
  )
  plot_index <- tibble::tibble(
    component_type = "plot",
    view = plot_views,
    available = as.logical(plot_available),
    n_rows = vapply(
      report_bundle$plots %||% list(),
      function(obj) if (is.data.frame(obj$data$data)) nrow(obj$data$data) else 0L,
      integer(1)
    ),
    x_var = vapply(
      report_bundle$plots %||% list(),
      function(obj) first_or_na(obj$data$x_var),
      character(1)
    ),
    group_var = vapply(
      report_bundle$plots %||% list(),
      function(obj) first_or_na(obj$data$group_var),
      character(1)
    ),
    recommended_design_id = vapply(
      report_bundle$plots %||% list(),
      function(obj) first_or_na(obj$data$recommended_design_id),
      character(1)
    )
  )
  component_index <- dplyr::bind_rows(
    tibble::tibble(
      component_type = "summary",
      view = NA_character_,
      available = isTRUE(report_bundle$grid_summary$grid_available),
      n_rows = as.integer(report_bundle$grid_summary$n_designs %||% 0L),
      x_var = NA_character_,
      group_var = NA_character_,
      recommended_design_id = NA_character_
    ),
    tibble::tibble(
      component_type = "recommendation",
      view = NA_character_,
      available = isTRUE(report_bundle$grid_recommendation$recommendation_available),
      n_rows = NA_integer_,
      x_var = NA_character_,
      group_var = NA_character_,
      recommended_design_id = first_or_na(report_bundle$grid_recommendation$recommended_design_id)
    ),
    table_index,
    plot_index
  )

  if (!isTRUE(report_bundle$report_available)) {
    return(list(
      report_summary_contract = "arbitrary_facet_design_grid_report_summary",
      report_summary_stage = as.character(report_bundle$report_stage %||% "schema_only"),
      planner_contract = as.character(
        report_bundle$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      report_available = FALSE,
      reason = as.character(
        report_bundle$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      n_designs = 0L,
      available_table_views = as.character(table_views[table_available]),
      available_plot_views = as.character(plot_views[plot_available]),
      component_index = component_index,
      report_bundle = report_bundle,
      note = paste(
        "Schema-only future-branch report summary is unavailable until the",
        "underlying branch grid is materialized. The component index still",
        "describes which branch-side report surfaces would be populated."
      )
    ))
  }

  list(
    report_summary_contract = "arbitrary_facet_design_grid_report_summary",
    report_summary_stage = as.character(report_bundle$report_stage %||% "schema_only"),
    planner_contract = as.character(
      report_bundle$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    report_available = TRUE,
    reason = as.character(
      report_bundle$reason %||%
        "Schema-only future-branch grid is materialized."
    ),
    n_designs = as.integer(report_bundle$grid_summary$n_designs %||% 0L),
    recommended_design_id = first_or_na(report_bundle$recommended_design_id),
    available_table_views = as.character(table_views[table_available]),
    available_plot_views = as.character(plot_views[plot_available]),
    component_index = component_index,
    report_bundle = report_bundle,
    note = paste(
      "Schema-only future-branch report summary exposing a compact component",
      "index over the current report bundle without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_overview_table <- function(x,
                                                           design = NULL,
                                                           prefer = NULL,
                                                           x_var = NULL,
                                                           group_var = NULL,
                                                           id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_overview_table",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  report_summary <- simulation_future_branch_report_summary(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  report_bundle <- report_summary$report_bundle
  overview_source <- tibble::as_tibble(report_bundle$overview_table %||% tibble::tibble())
  has_fixed_value <- "fixed_value" %in% names(overview_source)

  metrics_table <- tibble::tibble(
    report_available = isTRUE(report_summary$report_available),
    n_designs = as.integer(report_summary$n_designs %||% 0L),
    recommended_design_id = as.character(report_summary$recommended_design_id %||% NA_character_),
    n_table_views = length(report_summary$available_table_views %||% character(0)),
    n_plot_views = length(report_summary$available_plot_views %||% character(0)),
    available_table_views = paste(report_summary$available_table_views %||% character(0), collapse = ", "),
    available_plot_views = paste(report_summary$available_plot_views %||% character(0), collapse = ", ")
  )

  axis_overview_table <- if (nrow(overview_source) == 0L) {
    overview_source
  } else {
    overview_source |>
      dplyr::mutate(
        axis_state = dplyr::case_when(
          is.na(.data$varying) ~ "unavailable",
          .data$varying ~ "varying",
          TRUE ~ "fixed"
        ),
        value_summary = vapply(.data$values, function(vals) {
          if (length(vals) == 0L || all(is.na(vals))) {
            return(NA_character_)
          }
          paste(as.character(vals), collapse = ", ")
        }, character(1)),
        fixed_value = if (has_fixed_value) as.integer(.data$fixed_value) else NA_integer_
      ) |>
      dplyr::select(
        dplyr::all_of(c(
          "axis_source",
          "input_key",
          "canonical_design_variable",
          "public_design_alias",
          "axis_class",
          "facet",
          "axis_state",
          "n_values",
          "value_summary",
          "fixed_value"
        ))
      )
  }

  if (!isTRUE(report_summary$report_available)) {
    return(list(
      overview_contract = "arbitrary_facet_design_report_overview_table",
      overview_stage = as.character(report_summary$report_summary_stage %||% "schema_only"),
      planner_contract = as.character(
        report_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      report_available = FALSE,
      reason = as.character(
        report_summary$reason %||%
          "No schema-only future-branch grid is currently materialized."
      ),
      metrics_table = metrics_table,
      axis_overview_table = axis_overview_table,
      component_index = tibble::as_tibble(report_summary$component_index %||% tibble::tibble()),
      report_summary = report_summary,
      report_bundle = report_bundle,
      note = paste(
        "Schema-only future-branch overview table is unavailable until the",
        "underlying branch grid is materialized. Metrics and axis metadata",
        "are still returned from the current report summary contract."
      )
    ))
  }

  list(
    overview_contract = "arbitrary_facet_design_report_overview_table",
    overview_stage = as.character(report_summary$report_summary_stage %||% "schema_only"),
    planner_contract = as.character(
      report_summary$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    report_available = TRUE,
    reason = as.character(
      report_summary$reason %||%
        "Schema-only future-branch grid is materialized."
    ),
    metrics_table = metrics_table,
    axis_overview_table = axis_overview_table,
    component_index = tibble::as_tibble(report_summary$component_index %||% tibble::tibble()),
    report_summary = report_summary,
    report_bundle = report_bundle,
    note = paste(
      "Schema-only future-branch overview table exposing compact report",
      "metrics and axis-level state summaries without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_overview_view <- function(x,
                                                          component = c("metrics", "axes", "components"),
                                                          design = NULL,
                                                          prefer = NULL,
                                                          x_var = NULL,
                                                          group_var = NULL,
                                                          id_prefix = NULL) {
  component <- match.arg(component)
  overview <- simulation_future_branch_report_overview_table(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  table <- switch(
    component,
    metrics = tibble::as_tibble(overview$metrics_table %||% tibble::tibble()),
    axes = tibble::as_tibble(overview$axis_overview_table %||% tibble::tibble()),
    components = tibble::as_tibble(overview$component_index %||% tibble::tibble())
  )
  component_label <- switch(
    component,
    metrics = "report metrics",
    axes = "axis overview",
    components = "component index"
  )

  list(
    overview_view_contract = "arbitrary_facet_design_report_overview_view",
    overview_stage = as.character(overview$overview_stage %||% "schema_only"),
    planner_contract = as.character(
      overview$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    component = component,
    component_label = component_label,
    report_available = isTRUE(overview$report_available),
    reason = as.character(
      overview$reason %||%
        "No schema-only future-branch report overview is currently available."
    ),
    table = table,
    report_overview = overview,
    note = paste(
      "Schema-only future-branch overview view exposing the selected",
      component_label,
      "table from the compact report-overview contract without implying an",
      "active arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_catalog <- function(x,
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_catalog",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  component_ids <- c("metrics", "axes", "components")
  views <- stats::setNames(vector("list", length(component_ids)), component_ids)

  for (component in component_ids) {
    views[[component]] <- simulation_future_branch_report_overview_view(
      x = x,
      component = component,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  }

  overview <- views[[1]]$report_overview %||% simulation_future_branch_report_overview_table(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  recommended_design_id <- as.character(
    overview$metrics_table$recommended_design_id[[1]] %||% NA_character_
  )
  surface_index <- tibble::tibble(
    component = component_ids,
    component_label = vapply(views, function(obj) {
      as.character(obj$component_label %||% NA_character_)
    }, character(1)),
    available = vapply(views, function(obj) {
      isTRUE(obj$report_available)
    }, logical(1)),
    n_rows = vapply(views, function(obj) {
      if (is.data.frame(obj$table)) nrow(obj$table) else 0L
    }, integer(1)),
    recommended_design_id = rep(recommended_design_id, length(component_ids))
  )

  list(
    catalog_contract = "arbitrary_facet_design_report_catalog",
    catalog_stage = as.character(overview$overview_stage %||% "schema_only"),
    planner_contract = as.character(
      overview$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    report_available = isTRUE(overview$report_available),
    reason = as.character(
      overview$reason %||%
        "No schema-only future-branch report catalog is currently available."
    ),
    recommended_design_id = recommended_design_id,
    surface_index = surface_index,
    views = views,
    report_overview = overview,
    note = paste(
      "Schema-only future-branch report catalog enumerating the compact",
      "overview surfaces that branch-side reporting code can request from",
      "the current overview contract without implying an active arbitrary-",
      "facet planner."
    )
  )
}

simulation_future_branch_report_digest <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_digest",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  report_catalog <- simulation_future_branch_report_catalog(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  report_overview <- report_catalog$report_overview %||% simulation_future_branch_report_overview_table(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  metrics_table <- tibble::as_tibble(report_overview$metrics_table %||% tibble::tibble())
  axis_table <- tibble::as_tibble(report_overview$axis_overview_table %||% tibble::tibble())
  surface_index <- tibble::as_tibble(report_catalog$surface_index %||% tibble::tibble())

  first_chr <- function(value, default = NA_character_) {
    if (length(value) == 0L) {
      return(as.character(default))
    }
    as.character(value[[1]] %||% default)
  }
  first_int <- function(value, default = 0L) {
    if (length(value) == 0L) {
      return(as.integer(default))
    }
    as.integer(value[[1]] %||% default)
  }

  available_surfaces <- if (nrow(surface_index) == 0L) {
    character(0)
  } else {
    as.character(surface_index$component[replace(as.logical(surface_index$available), is.na(surface_index$available), FALSE)])
  }
  varying_axes <- if (nrow(axis_table) == 0L) {
    character(0)
  } else {
    as.character(axis_table$canonical_design_variable[axis_table$axis_state %in% "varying"])
  }
  fixed_axes <- if (nrow(axis_table) == 0L) {
    character(0)
  } else {
    as.character(axis_table$canonical_design_variable[axis_table$axis_state %in% "fixed"])
  }

  digest_table <- tibble::tibble(
    report_available = isTRUE(report_catalog$report_available),
    n_designs = first_int(metrics_table$n_designs, 0L),
    recommended_design_id = first_chr(metrics_table$recommended_design_id),
    n_available_surfaces = length(available_surfaces),
    available_surfaces = paste(available_surfaces, collapse = ", "),
    varying_axes = paste(varying_axes, collapse = ", "),
    fixed_axes = paste(fixed_axes, collapse = ", ")
  )

  list(
    digest_contract = "arbitrary_facet_design_report_digest",
    digest_stage = as.character(report_catalog$catalog_stage %||% "schema_only"),
    planner_contract = as.character(
      report_catalog$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    report_available = isTRUE(report_catalog$report_available),
    reason = as.character(
      report_catalog$reason %||%
        "No schema-only future-branch report digest is currently available."
    ),
    recommended_design_id = first_chr(metrics_table$recommended_design_id),
    available_surfaces = available_surfaces,
    varying_axes = varying_axes,
    fixed_axes = fixed_axes,
    digest_table = digest_table,
    report_catalog = report_catalog,
    report_overview = report_overview,
    note = paste(
      "Schema-only future-branch report digest exposing headline report",
      "availability, baseline design metadata, and compact surface/axis",
      "summaries from one internal contract without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_surface <- function(x,
                                                    surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  surface <- match.arg(surface)

  source_object <- switch(
    surface,
    digest = simulation_future_branch_report_digest(
      x = x,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    catalog = simulation_future_branch_report_catalog(
      x = x,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    metrics = simulation_future_branch_report_overview_view(
      x = x,
      component = "metrics",
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    axes = simulation_future_branch_report_overview_view(
      x = x,
      component = "axes",
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    components = simulation_future_branch_report_overview_view(
      x = x,
      component = "components",
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  )

  surface_label <- switch(
    surface,
    digest = "report digest",
    catalog = "report catalog",
    metrics = "report metrics",
    axes = "axis overview",
    components = "component index"
  )
  source_contract <- switch(
    surface,
    digest = as.character(source_object$digest_contract %||% "arbitrary_facet_design_report_digest"),
    catalog = as.character(source_object$catalog_contract %||% "arbitrary_facet_design_report_catalog"),
    metrics = as.character(source_object$overview_view_contract %||% "arbitrary_facet_design_report_overview_view"),
    axes = as.character(source_object$overview_view_contract %||% "arbitrary_facet_design_report_overview_view"),
    components = as.character(source_object$overview_view_contract %||% "arbitrary_facet_design_report_overview_view")
  )
  recommended_design_id <- switch(
    surface,
    digest = as.character(source_object$recommended_design_id %||% NA_character_),
    catalog = as.character(source_object$recommended_design_id %||% NA_character_),
    metrics = as.character(source_object$table$recommended_design_id[[1]] %||% NA_character_),
    axes = as.character(source_object$report_overview$metrics_table$recommended_design_id[[1]] %||% NA_character_),
    components = as.character(source_object$report_overview$metrics_table$recommended_design_id[[1]] %||% NA_character_)
  )
  table <- switch(
    surface,
    digest = tibble::as_tibble(source_object$digest_table %||% tibble::tibble()),
    catalog = {
      tbl <- tibble::as_tibble(source_object$surface_index %||% tibble::tibble())
      if ("component" %in% names(tbl) && !("surface" %in% names(tbl))) {
        names(tbl)[names(tbl) == "component"] <- "surface"
      }
      tbl
    },
    metrics = tibble::as_tibble(source_object$table %||% tibble::tibble()),
    axes = tibble::as_tibble(source_object$table %||% tibble::tibble()),
    components = tibble::as_tibble(source_object$table %||% tibble::tibble())
  )
  report_available <- switch(
    surface,
    digest = isTRUE(source_object$report_available),
    catalog = isTRUE(source_object$report_available),
    metrics = isTRUE(source_object$report_available),
    axes = isTRUE(source_object$report_available),
    components = isTRUE(source_object$report_available)
  )
  reason <- switch(
    surface,
    digest = as.character(source_object$reason %||% NA_character_),
    catalog = as.character(source_object$reason %||% NA_character_),
    metrics = as.character(source_object$reason %||% NA_character_),
    axes = as.character(source_object$reason %||% NA_character_),
    components = as.character(source_object$reason %||% NA_character_)
  )
  surface_stage <- switch(
    surface,
    digest = as.character(source_object$digest_stage %||% "schema_only"),
    catalog = as.character(source_object$catalog_stage %||% "schema_only"),
    metrics = as.character(source_object$overview_stage %||% "schema_only"),
    axes = as.character(source_object$overview_stage %||% "schema_only"),
    components = as.character(source_object$overview_stage %||% "schema_only")
  )
  planner_contract <- switch(
    surface,
    digest = as.character(source_object$planner_contract %||% "arbitrary_facet_planning_scaffold"),
    catalog = as.character(source_object$planner_contract %||% "arbitrary_facet_planning_scaffold"),
    metrics = as.character(source_object$planner_contract %||% "arbitrary_facet_planning_scaffold"),
    axes = as.character(source_object$planner_contract %||% "arbitrary_facet_planning_scaffold"),
    components = as.character(source_object$planner_contract %||% "arbitrary_facet_planning_scaffold")
  )

  list(
    surface_contract = "arbitrary_facet_design_report_surface",
    surface_stage = surface_stage,
    planner_contract = planner_contract,
    surface = surface,
    surface_label = surface_label,
    source_contract = source_contract,
    report_available = report_available,
    reason = reason,
    recommended_design_id = recommended_design_id,
    table = table,
    source_object = source_object,
    note = paste(
      "Schema-only future-branch report surface exposing the selected",
      surface_label,
      "from one compact branch-side operation without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_surface_registry <- function(x,
                                                             design = NULL,
                                                             prefer = NULL,
                                                             x_var = NULL,
                                                             group_var = NULL,
                                                             id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_surface_registry",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface_ids <- c("digest", "catalog", "metrics", "axes", "components")
  surfaces <- stats::setNames(vector("list", length(surface_ids)), surface_ids)

  for (surface in surface_ids) {
    surfaces[[surface]] <- simulation_future_branch_report_surface(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  }

  recommended_design_id <- as.character(
    surfaces$digest$recommended_design_id %||%
      surfaces$catalog$recommended_design_id %||%
      NA_character_
  )
  surface_index <- tibble::tibble(
    surface = surface_ids,
    surface_label = vapply(surfaces, function(obj) {
      as.character(obj$surface_label %||% NA_character_)
    }, character(1)),
    source_contract = vapply(surfaces, function(obj) {
      as.character(obj$source_contract %||% NA_character_)
    }, character(1)),
    available = vapply(surfaces, function(obj) {
      isTRUE(obj$report_available)
    }, logical(1)),
    n_rows = vapply(surfaces, function(obj) {
      if (is.data.frame(obj$table)) nrow(obj$table) else 0L
    }, integer(1)),
    recommended_design_id = rep(recommended_design_id, length(surface_ids))
  )

  list(
    registry_contract = "arbitrary_facet_design_report_surface_registry",
    registry_stage = as.character(surfaces[[1]]$surface_stage %||% "schema_only"),
    planner_contract = as.character(
      surfaces[[1]]$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    report_available = isTRUE(surfaces[[1]]$report_available),
    reason = as.character(
      surfaces[[1]]$reason %||%
        "No schema-only future-branch report surfaces are currently available."
    ),
    recommended_design_id = recommended_design_id,
    surface_index = surface_index,
    surfaces = surfaces,
    note = paste(
      "Schema-only future-branch report surface registry exposing digest,",
      "catalog, and compact overview surfaces from one internal contract",
      "without implying an active arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_panel <- function(x,
                                                  surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                  design = NULL,
                                                  prefer = NULL,
                                                  x_var = NULL,
                                                  group_var = NULL,
                                                  id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_panel",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  registry <- simulation_future_branch_report_surface_registry(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  selected_surface <- registry$surfaces[[surface]] %||% simulation_future_branch_report_surface(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  digest_table <- tibble::as_tibble(
    registry$surfaces$digest$table %||% tibble::tibble()
  )

  list(
    panel_contract = "arbitrary_facet_design_report_panel",
    panel_stage = as.character(registry$registry_stage %||% "schema_only"),
    planner_contract = as.character(
      registry$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    surface = surface,
    surface_label = as.character(selected_surface$surface_label %||% surface),
    report_available = isTRUE(selected_surface$report_available),
    reason = as.character(
      selected_surface$reason %||%
        "No schema-only future-branch report panel is currently available."
    ),
    recommended_design_id = as.character(
      selected_surface$recommended_design_id %||%
        registry$recommended_design_id %||%
        NA_character_
    ),
    digest_table = digest_table,
    surface_index = tibble::as_tibble(registry$surface_index %||% tibble::tibble()),
    selected_table = tibble::as_tibble(selected_surface$table %||% tibble::tibble()),
    selected_surface = selected_surface,
    registry = registry,
    note = paste(
      "Schema-only future-branch report panel exposing one selected compact",
      "surface together with headline digest metadata and the surface index",
      "from one internal contract without implying an active arbitrary-",
      "facet planner."
    )
  )
}

simulation_future_branch_report_operation <- function(x,
                                                      surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                      design = NULL,
                                                      prefer = NULL,
                                                      x_var = NULL,
                                                      group_var = NULL,
                                                      id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_operation",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  panel <- simulation_future_branch_report_panel(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  registry <- panel$registry %||% simulation_future_branch_report_surface_registry(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  digest_surface <- registry$surfaces$digest %||% simulation_future_branch_report_surface(
    x = x,
    surface = "digest",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  metrics_surface <- registry$surfaces$metrics %||% simulation_future_branch_report_surface(
    x = x,
    surface = "metrics",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  axes_surface <- registry$surfaces$axes %||% simulation_future_branch_report_surface(
    x = x,
    surface = "axes",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  components_surface <- registry$surfaces$components %||% simulation_future_branch_report_surface(
    x = x,
    surface = "components",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  list(
    operation_contract = "arbitrary_facet_design_report_operation",
    operation_stage = as.character(
      panel$panel_stage %||% registry$registry_stage %||% "schema_only"
    ),
    planner_contract = as.character(
      panel$planner_contract %||%
        registry$planner_contract %||%
        "arbitrary_facet_planning_scaffold"
    ),
    surface = surface,
    surface_label = as.character(panel$surface_label %||% surface),
    report_available = isTRUE(panel$report_available),
    reason = as.character(
      panel$reason %||%
        "No schema-only future-branch report operation is currently available."
    ),
    recommended_design_id = as.character(
      panel$recommended_design_id %||%
        registry$recommended_design_id %||%
        NA_character_
    ),
    digest_table = tibble::as_tibble(digest_surface$table %||% tibble::tibble()),
    surface_index = tibble::as_tibble(registry$surface_index %||% tibble::tibble()),
    metrics_table = tibble::as_tibble(metrics_surface$table %||% tibble::tibble()),
    axis_overview_table = tibble::as_tibble(axes_surface$table %||% tibble::tibble()),
    component_index = tibble::as_tibble(components_surface$table %||% tibble::tibble()),
    selected_table = tibble::as_tibble(panel$selected_table %||% tibble::tibble()),
    selected_surface = panel$selected_surface,
    report_panel = panel,
    report_surface_registry = registry,
    note = paste(
      "Schema-only future-branch report operation exposing digest, compact",
      "overview tables, the current surface registry, and one selected",
      "surface from one internal contract without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_future_branch_report_snapshot <- function(x,
                                                     surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                     design = NULL,
                                                     prefer = NULL,
                                                     x_var = NULL,
                                                     group_var = NULL,
                                                     id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_snapshot",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  operation <- simulation_future_branch_report_operation(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  first_chr <- function(value, default = NA_character_) {
    if (length(value) == 0L) {
      return(as.character(default))
    }
    as.character(value[[1]] %||% default)
  }

  digest_tbl <- tibble::as_tibble(operation$digest_table %||% tibble::tibble())
  surface_index <- tibble::as_tibble(operation$surface_index %||% tibble::tibble())
  selected_tbl <- tibble::as_tibble(operation$selected_table %||% tibble::tibble())

  available_surfaces <- if (nrow(surface_index) == 0L) {
    character(0)
  } else {
    as.character(surface_index$surface[replace(as.logical(surface_index$available), is.na(surface_index$available), FALSE)])
  }
  varying_axes <- if (!("varying_axes" %in% names(digest_tbl))) {
    character(0)
  } else {
    strsplit(first_chr(digest_tbl$varying_axes, ""), ", ", fixed = TRUE)[[1]]
  }
  varying_axes <- varying_axes[nzchar(varying_axes)]
  fixed_axes <- if (!("fixed_axes" %in% names(digest_tbl))) {
    character(0)
  } else {
    strsplit(first_chr(digest_tbl$fixed_axes, ""), ", ", fixed = TRUE)[[1]]
  }
  fixed_axes <- fixed_axes[nzchar(fixed_axes)]

  list(
    snapshot_contract = "arbitrary_facet_design_report_snapshot",
    snapshot_stage = as.character(operation$operation_stage %||% "schema_only"),
    planner_contract = as.character(
      operation$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    surface = surface,
    surface_label = as.character(operation$surface_label %||% surface),
    report_available = isTRUE(operation$report_available),
    reason = as.character(
      operation$reason %||%
        "No schema-only future-branch report snapshot is currently available."
    ),
    recommended_design_id = as.character(
      operation$recommended_design_id %||% NA_character_
    ),
    n_designs = if ("n_designs" %in% names(digest_tbl)) as.integer(digest_tbl$n_designs[[1]] %||% 0L) else 0L,
    available_surfaces = available_surfaces,
    varying_axes = varying_axes,
    fixed_axes = fixed_axes,
    digest_table = digest_tbl,
    surface_index = surface_index,
    selected_table = selected_tbl,
    metrics_table = tibble::as_tibble(operation$metrics_table %||% tibble::tibble()),
    axis_overview_table = tibble::as_tibble(operation$axis_overview_table %||% tibble::tibble()),
    component_index = tibble::as_tibble(operation$component_index %||% tibble::tibble()),
    note = paste(
      "Schema-only future-branch report snapshot exposing compact headline",
      "metadata, compact overview tables, and one selected surface without",
      "carrying the deeper nested operation graph."
    )
  )
}

simulation_future_branch_report_brief <- function(x,
                                                  surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                  design = NULL,
                                                  prefer = NULL,
                                                  x_var = NULL,
                                                  group_var = NULL,
                                                  id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_brief",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  snapshot <- simulation_future_branch_report_snapshot(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  headline_table <- tibble::tibble(
    report_available = isTRUE(snapshot$report_available),
    n_designs = as.integer(snapshot$n_designs %||% 0L),
    recommended_design_id = as.character(
      snapshot$recommended_design_id %||% NA_character_
    ),
    surface = as.character(snapshot$surface %||% surface),
    surface_label = as.character(snapshot$surface_label %||% surface),
    available_surfaces = paste(snapshot$available_surfaces %||% character(0), collapse = ", "),
    varying_axes = paste(snapshot$varying_axes %||% character(0), collapse = ", "),
    fixed_axes = paste(snapshot$fixed_axes %||% character(0), collapse = ", ")
  )

  list(
    brief_contract = "arbitrary_facet_design_report_brief",
    brief_stage = as.character(snapshot$snapshot_stage %||% "schema_only"),
    planner_contract = as.character(
      snapshot$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    surface = surface,
    surface_label = as.character(snapshot$surface_label %||% surface),
    report_available = isTRUE(snapshot$report_available),
    reason = as.character(
      snapshot$reason %||%
        "No schema-only future-branch report brief is currently available."
    ),
    headline_table = headline_table,
    selected_table = tibble::as_tibble(snapshot$selected_table %||% tibble::tibble()),
    surface_index = tibble::as_tibble(snapshot$surface_index %||% tibble::tibble()),
    note = paste(
      "Schema-only future-branch report brief exposing one selected surface",
      "together with a headline table and surface index, without carrying",
      "the broader snapshot or nested operation graph."
    )
  )
}

simulation_future_branch_report_consume <- function(x,
                                                    mode = c("brief", "snapshot", "operation"),
                                                    surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "report_consumer",
    design = design,
    prefer = prefer,
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  mode <- match.arg(mode)
  surface <- match.arg(surface)

  payload <- switch(
    mode,
    brief = simulation_future_branch_report_brief(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    snapshot = simulation_future_branch_report_snapshot(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    operation = simulation_future_branch_report_operation(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  )

  payload_contract <- switch(
    mode,
    brief = as.character(payload$brief_contract %||% NA_character_),
    snapshot = as.character(payload$snapshot_contract %||% NA_character_),
    operation = as.character(payload$operation_contract %||% NA_character_)
  )
  stage <- switch(
    mode,
    brief = as.character(payload$brief_stage %||% "schema_only"),
    snapshot = as.character(payload$snapshot_stage %||% "schema_only"),
    operation = as.character(payload$operation_stage %||% "schema_only")
  )
  selected_table <- switch(
    mode,
    brief = tibble::as_tibble(payload$selected_table %||% tibble::tibble()),
    snapshot = tibble::as_tibble(payload$selected_table %||% tibble::tibble()),
    operation = tibble::as_tibble(payload$selected_table %||% tibble::tibble())
  )
  surface_index <- switch(
    mode,
    brief = tibble::as_tibble(payload$surface_index %||% tibble::tibble()),
    snapshot = tibble::as_tibble(payload$surface_index %||% tibble::tibble()),
    operation = tibble::as_tibble(payload$surface_index %||% tibble::tibble())
  )
  recommended_design_id <- switch(
    mode,
    brief = as.character(payload$headline_table$recommended_design_id[[1]] %||% NA_character_),
    snapshot = as.character(payload$recommended_design_id %||% NA_character_),
    operation = as.character(payload$recommended_design_id %||% NA_character_)
  )
  n_designs <- switch(
    mode,
    brief = as.integer(payload$headline_table$n_designs[[1]] %||% 0L),
    snapshot = as.integer(payload$n_designs %||% 0L),
    operation = if ("n_designs" %in% names(payload$digest_table)) {
      as.integer(payload$digest_table$n_designs[[1]] %||% 0L)
    } else {
      0L
    }
  )

  list(
    consumer_contract = "arbitrary_facet_design_report_consumer",
    consumer_stage = stage,
    planner_contract = as.character(
      payload$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    mode = mode,
    payload_contract = payload_contract,
    surface = surface,
    surface_label = as.character(payload$surface_label %||% surface),
    report_available = isTRUE(payload$report_available),
    reason = as.character(
      payload$reason %||%
        "No schema-only future-branch report payload is currently available."
    ),
    recommended_design_id = recommended_design_id,
    n_designs = n_designs,
    selected_table = selected_table,
    surface_index = surface_index,
    payload = payload,
    note = paste(
      "Schema-only future-branch report consumer dispatching to the selected",
      mode,
      "contract so branch-side code can switch report payload weight without",
      "re-implementing the dispatch logic."
    )
  )
}

simulation_future_branch_report_mode_registry <- function(x) {
  cached <- if (is.list(x)) x$report_mode_registry %||% NULL else NULL
  if (is.list(cached) || (is.data.frame(cached) && nrow(cached) >= 0L)) {
    return(cached)
  }

  tibble::tibble(
    mode = c("brief", "snapshot", "operation"),
    payload_contract = c(
      "arbitrary_facet_design_report_brief",
      "arbitrary_facet_design_report_snapshot",
      "arbitrary_facet_design_report_operation"
    ),
    default_surface = "digest",
    carries_nested_graph = c(FALSE, FALSE, TRUE),
    note = c(
      "Selected surface plus headline table and surface index.",
      "Compact headline metadata plus one selected surface.",
      "Full compact report operation with selected surface and nested registry."
    )
  )
}

simulation_future_branch_pilot <- function(x,
                                           design = NULL,
                                           prefer = NULL,
                                           view = c("public", "canonical", "branch"),
                                           mode = c("brief", "snapshot", "operation"),
                                           surface = c("digest", "catalog", "metrics", "axes", "components"),
                                           x_var = NULL,
                                           group_var = NULL,
                                           id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "pilot",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  view <- match.arg(view)
  mode <- match.arg(mode)
  surface <- match.arg(surface)

  grid_summary <- simulation_future_branch_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  grid_recommendation <- simulation_future_branch_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )
  grid_table <- simulation_future_branch_grid_table(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    id_prefix = id_prefix
  )
  plot_payload <- simulation_future_branch_grid_plot_payload(
    x = x,
    design = design,
    x_var = x_var,
    group_var = group_var,
    prefer = prefer,
    view = view,
    id_prefix = id_prefix
  )
  report_consumer <- simulation_future_branch_report_consume(
    x = x,
    mode = mode,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  pilot_available <- isTRUE(grid_summary$grid_available) && isTRUE(report_consumer$report_available)

  list(
    pilot_contract = "arbitrary_facet_planning_pilot",
    pilot_stage = if (pilot_available) "pilot_active" else "schema_only",
    planner_contract = as.character(
      grid_summary$planner_contract %||%
        report_consumer$planner_contract %||%
        "arbitrary_facet_planning_scaffold"
    ),
    pilot_available = pilot_available,
    reason = as.character(
      if (pilot_available) {
        "Schema-only future-branch pilot is materialized from the current grid and report consumer contracts."
      } else {
        grid_summary$reason %||%
          report_consumer$reason %||%
          "No schema-only future-branch pilot is currently available."
      }
    ),
    view = view,
    mode = mode,
    surface = surface,
    recommended_design_id = as.character(
      report_consumer$recommended_design_id %||%
        grid_recommendation$recommended_design_id %||%
        NA_character_
    ),
    n_designs = as.integer(
      report_consumer$n_designs %||%
        grid_summary$n_designs %||%
        0L
    ),
    grid_table = tibble::as_tibble(grid_table$table %||% tibble::tibble()),
    plot_payload = plot_payload,
    report_consumer = report_consumer,
    grid_summary = grid_summary,
    grid_recommendation = grid_recommendation,
    note = paste(
      "Internal schema-only future-branch pilot bundling one materialized grid",
      "view, one draw-free plot payload, and one mode-selected report consumer.",
      "This is the first active branch-side object, but it remains a",
      "deterministic scaffold rather than a performance-based arbitrary-facet planner."
    )
  )
}

simulation_future_branch_pilot_summary <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   view = c("public", "canonical", "branch"),
                                                   mode = c("brief", "snapshot", "operation"),
                                                   surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "pilot_summary",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  pilot <- simulation_future_branch_pilot(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  headline_table <- tibble::tibble(
    pilot_available = isTRUE(pilot$pilot_available),
    n_designs = as.integer(pilot$n_designs %||% 0L),
    recommended_design_id = as.character(pilot$recommended_design_id %||% NA_character_),
    view = as.character(pilot$view %||% NA_character_),
    mode = as.character(pilot$mode %||% NA_character_),
    surface = as.character(pilot$surface %||% NA_character_),
    grid_rows = if (is.data.frame(pilot$grid_table)) nrow(pilot$grid_table) else 0L,
    plot_available = isTRUE(pilot$plot_payload$data$plot_available %||% FALSE),
    report_payload_contract = as.character(
      pilot$report_consumer$payload_contract %||% NA_character_
    )
  )

  list(
    pilot_summary_contract = "arbitrary_facet_planning_pilot_summary",
    pilot_stage = as.character(pilot$pilot_stage %||% "schema_only"),
    planner_contract = as.character(
      pilot$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    pilot_available = isTRUE(pilot$pilot_available),
    reason = as.character(
      pilot$reason %||%
        "No schema-only future-branch pilot summary is currently available."
    ),
    headline_table = headline_table,
    pilot = pilot,
    note = paste(
      "Compact summary for the schema-only future-branch pilot, exposing one",
      "headline table over the current grid/report scaffold without implying",
      "an active performance-based arbitrary-facet planner."
    )
  )
}

simulation_future_branch_pilot_table <- function(x,
                                                 component = c("grid", "report", "surface_index"),
                                                 design = NULL,
                                                 prefer = NULL,
                                                 view = c("public", "canonical", "branch"),
                                                 mode = c("brief", "snapshot", "operation"),
                                                 surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                 x_var = NULL,
                                                 group_var = NULL,
                                                 id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "pilot_table",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    component = component,
    component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  component <- match.arg(component)
  pilot <- simulation_future_branch_pilot(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  table <- switch(
    component,
    grid = tibble::as_tibble(pilot$grid_table %||% tibble::tibble()),
    report = tibble::as_tibble(pilot$report_consumer$selected_table %||% tibble::tibble()),
    surface_index = tibble::as_tibble(pilot$report_consumer$surface_index %||% tibble::tibble())
  )

  list(
    pilot_table_contract = "arbitrary_facet_planning_pilot_table",
    pilot_stage = as.character(pilot$pilot_stage %||% "schema_only"),
    planner_contract = as.character(
      pilot$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    pilot_available = isTRUE(pilot$pilot_available),
    reason = as.character(
      pilot$reason %||%
        "No schema-only future-branch pilot table is currently available."
    ),
    component = component,
    table = table,
    pilot = pilot,
    note = paste(
      "Selected table view from the schema-only future-branch pilot.",
      "This exposes either the grid table, the selected report table, or",
      "the compact surface index from the current pilot scaffold."
    )
  )
}

simulation_future_branch_pilot_plot <- function(x,
                                                design = NULL,
                                                prefer = NULL,
                                                view = c("public", "canonical", "branch"),
                                                mode = c("brief", "snapshot", "operation"),
                                                surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                x_var = NULL,
                                                group_var = NULL,
                                                id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "pilot_plot",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  pilot <- simulation_future_branch_pilot(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  list(
    pilot_plot_contract = "arbitrary_facet_planning_pilot_plot",
    pilot_stage = as.character(pilot$pilot_stage %||% "schema_only"),
    planner_contract = as.character(
      pilot$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    pilot_available = isTRUE(pilot$pilot_available),
    reason = as.character(
      pilot$reason %||%
        "No schema-only future-branch pilot plot is currently available."
    ),
    plot = pilot$plot_payload,
    pilot = pilot,
    note = paste(
      "Draw-free plotting payload from the schema-only future-branch pilot.",
      "This preserves the existing plot contract while routing access through",
      "the pilot scaffold."
    )
  )
}

simulation_future_branch_active_branch <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   view = c("public", "canonical", "branch"),
                                                   mode = c("brief", "snapshot", "operation"),
                                                   surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                   table_component = c("grid", "report", "surface_index"),
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  if (inherits(x, "mfrm_future_branch_active_branch") ||
      identical(as.character(x$branch_contract %||% NA_character_), "arbitrary_facet_planning_active_branch")) {
    return(x)
  }

  view <- match.arg(view)
  mode <- match.arg(mode)
  surface <- match.arg(surface)
  table_component <- match.arg(table_component)

  pilot_summary <- simulation_future_branch_pilot_summary(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  pilot_table <- simulation_future_branch_pilot_table(
    x = x,
    component = table_component,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  pilot_plot <- simulation_future_branch_pilot_plot(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  branch_available <- isTRUE(pilot_summary$pilot_available) &&
    isTRUE(pilot_table$pilot_available) &&
    isTRUE(pilot_plot$pilot_available)
  canonical_grid <- tibble::as_tibble(
    pilot_summary$pilot$grid_summary$canonical %||% tibble::tibble()
  )

  structure(list(
    branch_contract = "arbitrary_facet_planning_active_branch",
    branch_stage = if (branch_available) "pilot_active" else "schema_only",
    planner_contract = as.character(
      pilot_summary$planner_contract %||%
        pilot_table$planner_contract %||%
        pilot_plot$planner_contract %||%
        "arbitrary_facet_planning_scaffold"
    ),
    branch_available = branch_available,
    reason = as.character(
      if (branch_available) {
        "Schema-only future-branch active branch is materialized from the pilot summary, table, and plot contracts."
      } else {
        pilot_summary$reason %||%
          pilot_table$reason %||%
          pilot_plot$reason %||%
          "No schema-only future-branch active branch is currently available."
      }
    ),
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    recommended_design_id = as.character(
      pilot_summary$headline_table$recommended_design_id[[1]] %||% NA_character_
    ),
    n_designs = as.integer(
      pilot_summary$headline_table$n_designs[[1]] %||% 0L
    ),
    summary = pilot_summary,
    table = pilot_table,
    plot = pilot_plot,
    canonical_grid = canonical_grid,
    note = paste(
      "Minimal active future-branch object bundling pilot-level summary, one",
      "selected table component, and one draw-free plot payload from the",
      "current arbitrary-facet scaffold."
    )
  ), class = c("mfrm_future_branch_active_branch", "list"))
}

simulation_future_branch_active_branch_canonical_grid <- function(branch) {
  candidates <- list(
    branch$canonical_grid,
    branch$summary$pilot$grid_summary$canonical,
    branch$table$pilot$grid_summary$canonical,
    branch$plot$pilot$grid_summary$canonical
  )

  for (candidate in candidates) {
    if (is.data.frame(candidate) && nrow(candidate) > 0L) {
      return(tibble::as_tibble(candidate))
    }
  }

  tibble::tibble()
}

simulation_future_branch_active_branch_metric_registry <- function() {
  tibble::tibble(
    metric = c(
      "total_observations",
      "observations_per_person",
      "observations_per_criterion",
      "expected_observations_per_rater",
      "assignment_fraction"
    ),
    basis_class = c(
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "balanced_expectation",
      "density_ratio"
    ),
    formula = c(
      "n_person * raters_per_person * n_criterion",
      "raters_per_person * n_criterion",
      "n_person * raters_per_person",
      "(n_person * raters_per_person * n_criterion) / n_rater",
      "raters_per_person / n_rater"
    ),
    interpretation = c(
      "Total number of scored observations implied by the current design row.",
      "Number of scored observations contributed by one person under the current design row.",
      "Number of scored observations contributed to one criterion across persons.",
      "Average expected rater load under balanced assignment, not a guaranteed realized count.",
      "Fraction of available raters assigned to each person."
    ),
    psychometric = FALSE
  )
}

simulation_future_branch_active_branch_load_balance_registry <- function() {
  tibble::tibble(
    metric = c(
      "expected_observations_per_rater",
      "expected_person_assignments_per_rater",
      "observation_load_floor",
      "observation_load_ceiling",
      "observation_load_remainder",
      "perfect_integer_observation_balance"
    ),
    basis_class = c(
      "balanced_expectation",
      "balanced_expectation",
      "integer_balance_bound",
      "integer_balance_bound",
      "exact_identity",
      "exact_identity"
    ),
    formula = c(
      "(n_person * raters_per_person * n_criterion) / n_rater",
      "(n_person * raters_per_person) / n_rater",
      "floor((n_person * raters_per_person * n_criterion) / n_rater)",
      "ceiling((n_person * raters_per_person * n_criterion) / n_rater)",
      "(n_person * raters_per_person * n_criterion) %% n_rater",
      "as.integer(((n_person * raters_per_person * n_criterion) %% n_rater) == 0)"
    ),
    interpretation = c(
      "Average observation load per rater under balanced assignment, not a guaranteed realized count.",
      "Average person-assignment load per rater under balanced assignment, not a guaranteed realized count.",
      "Lower integer observation count per rater under the most even observation-level split implied by the current design row.",
      "Upper integer observation count per rater under the most even observation-level split implied by the current design row.",
      "Number of observations left over after equal integer splitting across raters.",
      "Indicator that the total observation count is divisible by the number of raters."
    ),
    psychometric = FALSE
  )
}

simulation_future_branch_active_branch_coverage_registry <- function() {
  tibble::tibble(
    metric = c(
      "person_criterion_cells",
      "criterion_replications_per_person",
      "rater_pair_overlap_per_cell",
      "total_rater_pair_overlaps",
      "pair_coverage_fraction_per_cell",
      "redundant_scoring"
    ),
    basis_class = c(
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "exact_ratio",
      "exact_identity"
    ),
    formula = c(
      "n_person * n_criterion",
      "raters_per_person",
      "choose(raters_per_person, 2)",
      "n_person * n_criterion * choose(raters_per_person, 2)",
      "ifelse(choose(n_rater, 2) > 0, choose(raters_per_person, 2) / choose(n_rater, 2), 0)",
      "as.integer(raters_per_person > 1)"
    ),
    interpretation = c(
      "Number of person-by-criterion cells implied by the current design row.",
      "Number of ratings attached to each person-by-criterion cell.",
      "Number of distinct rater pairs overlapping within one scored cell.",
      "Total number of rater-pair overlaps implied across all person-by-criterion cells.",
      "Fraction of available rater pairs represented within one scored cell; set to 0 when fewer than two raters are available.",
      "Indicator that the design uses more than one rater per scored cell."
    ),
    psychometric = FALSE
  )
}

simulation_future_branch_active_branch_profile <- function(x,
                                                           design = NULL,
                                                           prefer = NULL,
                                                           view = c("public", "canonical", "branch"),
                                                           mode = c("brief", "snapshot", "operation"),
                                                           surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                           table_component = c("grid", "report", "surface_index"),
                                                           x_var = NULL,
                                                           group_var = NULL,
                                                           id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_profile",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  branch <- simulation_future_branch_active_branch(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  canonical <- simulation_future_branch_active_branch_canonical_grid(branch)

  if (!isTRUE(branch$branch_available) || nrow(canonical) == 0L) {
    return(list(
      profile_contract = "arbitrary_facet_planning_active_branch_profile",
      profile_stage = as.character(branch$branch_stage %||% "schema_only"),
      planner_contract = as.character(
        branch$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      branch_available = FALSE,
      reason = as.character(
        branch$reason %||%
          "No schema-only future-branch active profile is currently available."
      ),
      metric_registry = simulation_future_branch_active_branch_metric_registry(),
      profile_summary_table = tibble::tibble(),
      profile_table = tibble::tibble(),
      active_branch = branch,
      note = paste(
        "Deterministic design-profile metrics are unavailable until the",
        "active future-branch scaffold can be materialized."
      )
    ))
  }

  profile_table <- canonical |>
    dplyr::mutate(
      total_observations = as.numeric(.data$n_person) * as.numeric(.data$raters_per_person) * as.numeric(.data$n_criterion),
      observations_per_person = as.numeric(.data$raters_per_person) * as.numeric(.data$n_criterion),
      observations_per_criterion = as.numeric(.data$n_person) * as.numeric(.data$raters_per_person),
      expected_observations_per_rater = .data$total_observations / as.numeric(.data$n_rater),
      assignment_fraction = as.numeric(.data$raters_per_person) / as.numeric(.data$n_rater),
      recommended = .data$design_id == as.character(branch$recommended_design_id %||% NA_character_)
    )
  metric_registry <- simulation_future_branch_active_branch_metric_registry()
  recommended_row <- profile_table[profile_table$recommended %in% TRUE, , drop = FALSE]
  profile_summary_table <- metric_registry |>
    dplyr::mutate(
      min = vapply(.data$metric, function(metric) {
        min(profile_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$metric, function(metric) {
        max(profile_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$metric, function(metric) {
        mean(profile_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$metric, function(metric) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[metric]][[1]])
      }, numeric(1))
    )

  list(
    profile_contract = "arbitrary_facet_planning_active_branch_profile",
    profile_stage = as.character(branch$branch_stage %||% "schema_only"),
    planner_contract = as.character(
      branch$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = TRUE,
    reason = as.character(
      branch$reason %||%
        "Schema-only future-branch active profile is materialized."
    ),
    metric_registry = metric_registry,
    profile_summary_table = profile_summary_table,
    profile_table = profile_table,
    active_branch = branch,
    note = paste(
      "Deterministic design-bookkeeping profile for the active future-branch",
      "scaffold. These quantities summarize observation counts, expected rater",
      "load, and assignment density; they are not psychometric performance estimates."
    )
  )
}

simulation_future_branch_active_branch_load_balance <- function(x,
                                                                design = NULL,
                                                                prefer = NULL,
                                                                view = c("public", "canonical", "branch"),
                                                                mode = c("brief", "snapshot", "operation"),
                                                                surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                table_component = c("grid", "report", "surface_index"),
                                                                x_var = NULL,
                                                                group_var = NULL,
                                                                id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_load_balance",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  profile <- simulation_future_branch_active_branch_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  profile_table <- tibble::as_tibble(profile$profile_table %||% tibble::tibble())
  metric_registry <- simulation_future_branch_active_branch_load_balance_registry()

  if (!isTRUE(profile$branch_available) || nrow(profile_table) == 0L) {
    return(list(
      diagnostics_contract = "arbitrary_facet_planning_active_branch_load_balance",
      diagnostics_stage = as.character(profile$profile_stage %||% "schema_only"),
      planner_contract = as.character(
        profile$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      branch_available = FALSE,
      reason = as.character(
        profile$reason %||%
          "No schema-only future-branch load/balance diagnostics are currently available."
      ),
      metric_registry = metric_registry,
      diagnostic_summary_table = tibble::tibble(),
      diagnostic_table = tibble::tibble(),
      active_branch_profile = profile,
      note = paste(
        "Deterministic load/balance diagnostics are unavailable until the",
        "active future-branch scaffold can be materialized."
      )
    ))
  }

  diagnostic_table <- profile_table |>
    dplyr::mutate(
      expected_person_assignments_per_rater = (as.numeric(.data$n_person) * as.numeric(.data$raters_per_person)) / as.numeric(.data$n_rater),
      observation_load_floor = floor(as.numeric(.data$total_observations) / as.numeric(.data$n_rater)),
      observation_load_ceiling = ceiling(as.numeric(.data$total_observations) / as.numeric(.data$n_rater)),
      observation_load_remainder = as.numeric(.data$total_observations) %% as.numeric(.data$n_rater),
      perfect_integer_observation_balance = as.integer(.data$observation_load_remainder == 0)
    )
  recommended_row <- diagnostic_table[diagnostic_table$recommended %in% TRUE, , drop = FALSE]
  diagnostic_summary_table <- metric_registry |>
    dplyr::mutate(
      min = vapply(.data$metric, function(metric) {
        min(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$metric, function(metric) {
        max(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$metric, function(metric) {
        mean(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$metric, function(metric) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[metric]][[1]])
      }, numeric(1))
    )

  list(
    diagnostics_contract = "arbitrary_facet_planning_active_branch_load_balance",
    diagnostics_stage = as.character(profile$profile_stage %||% "schema_only"),
    planner_contract = as.character(
      profile$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = TRUE,
    reason = as.character(
      profile$reason %||%
        "Schema-only future-branch load/balance diagnostics are materialized."
    ),
    metric_registry = metric_registry,
    diagnostic_summary_table = diagnostic_summary_table,
    diagnostic_table = diagnostic_table,
    active_branch_profile = profile,
    note = paste(
      "Deterministic load/balance diagnostics for the active future-branch",
      "scaffold. Balanced-load quantities are labeled as expectations, while",
      "integer split diagnostics are exact combinatorial summaries of the",
      "current observation counts."
    )
  )
}

simulation_future_branch_active_branch_overview <- function(x,
                                                            design = NULL,
                                                            prefer = NULL,
                                                            view = c("public", "canonical", "branch"),
                                                            mode = c("brief", "snapshot", "operation"),
                                                            surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                            table_component = c("grid", "report", "surface_index"),
                                                            x_var = NULL,
                                                            group_var = NULL,
                                                            id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  branch <- simulation_future_branch_active_branch(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  profile <- simulation_future_branch_active_branch_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  headline_table <- tibble::tibble(
    branch_available = isTRUE(branch$branch_available),
    n_designs = as.integer(branch$n_designs %||% 0L),
    recommended_design_id = as.character(branch$recommended_design_id %||% NA_character_),
    view = as.character(branch$view %||% NA_character_),
    mode = as.character(branch$mode %||% NA_character_),
    surface = as.character(branch$surface %||% NA_character_),
    table_component = as.character(branch$table_component %||% NA_character_),
    selected_table_rows = if (is.data.frame(branch$table$table)) nrow(branch$table$table) else 0L,
    plot_available = isTRUE(branch$plot$plot$data$plot_available %||% FALSE),
    n_metrics = if (is.data.frame(profile$metric_registry)) nrow(profile$metric_registry) else 0L
  )

  list(
    overview_contract = "arbitrary_facet_planning_active_branch_overview",
    overview_stage = as.character(branch$branch_stage %||% "schema_only"),
    planner_contract = as.character(
      branch$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = isTRUE(branch$branch_available) && isTRUE(profile$branch_available),
    reason = as.character(
      branch$reason %||%
        profile$reason %||%
        "No future-branch active overview is currently available."
    ),
    headline_table = headline_table,
    metric_registry = tibble::as_tibble(profile$metric_registry %||% tibble::tibble()),
    metric_summary_table = tibble::as_tibble(profile$profile_summary_table %||% tibble::tibble()),
    active_branch = branch,
    active_branch_profile = profile,
    note = paste(
      "Compact overview for the active future-branch scaffold, combining the",
      "branch headline with metric-basis-aware deterministic design summaries."
    )
  )
}

simulation_future_branch_active_branch_load_balance_overview <- function(x,
                                                                         design = NULL,
                                                                         prefer = NULL,
                                                                         view = c("public", "canonical", "branch"),
                                                                         mode = c("brief", "snapshot", "operation"),
                                                                         surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                         table_component = c("grid", "report", "surface_index"),
                                                                         x_var = NULL,
                                                                         group_var = NULL,
                                                                         id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_load_balance_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  diagnostics <- simulation_future_branch_active_branch_load_balance(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  diagnostic_table <- tibble::as_tibble(diagnostics$diagnostic_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    branch_available = isTRUE(diagnostics$branch_available),
    n_designs = nrow(diagnostic_table),
    n_metrics = if (is.data.frame(diagnostics$metric_registry)) nrow(diagnostics$metric_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(diagnostic_table) &&
                                any(diagnostic_table$recommended %in% TRUE)) {
      as.character(diagnostic_table$design_id[which(diagnostic_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_perfect_integer_balance = if ("perfect_integer_observation_balance" %in% names(diagnostic_table)) {
      sum(diagnostic_table$perfect_integer_observation_balance %in% 1L, na.rm = TRUE)
    } else {
      0L
    },
    n_nondivisible_designs = if ("perfect_integer_observation_balance" %in% names(diagnostic_table)) {
      sum(diagnostic_table$perfect_integer_observation_balance %in% 0L, na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_active_branch_load_balance_overview",
    overview_stage = as.character(diagnostics$diagnostics_stage %||% "schema_only"),
    planner_contract = as.character(
      diagnostics$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = isTRUE(diagnostics$branch_available),
    reason = as.character(
      diagnostics$reason %||%
        "No future-branch load/balance overview is currently available."
    ),
    headline_table = headline_table,
    metric_registry = tibble::as_tibble(diagnostics$metric_registry %||% tibble::tibble()),
    diagnostic_summary_table = tibble::as_tibble(diagnostics$diagnostic_summary_table %||% tibble::tibble()),
    active_branch_load_balance = diagnostics,
    note = paste(
      "Compact load/balance overview for the active future-branch scaffold,",
      "combining a one-row divisibility headline with basis-aware deterministic",
      "observation-load diagnostics."
    )
  )
}

simulation_future_branch_active_branch_coverage <- function(x,
                                                            design = NULL,
                                                            prefer = NULL,
                                                            view = c("public", "canonical", "branch"),
                                                            mode = c("brief", "snapshot", "operation"),
                                                            surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                            table_component = c("grid", "report", "surface_index"),
                                                            x_var = NULL,
                                                            group_var = NULL,
                                                            id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_coverage",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  profile <- simulation_future_branch_active_branch_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  profile_table <- tibble::as_tibble(profile$profile_table %||% tibble::tibble())
  metric_registry <- simulation_future_branch_active_branch_coverage_registry()

  if (!isTRUE(profile$branch_available) || nrow(profile_table) == 0L) {
    return(list(
      diagnostics_contract = "arbitrary_facet_planning_active_branch_coverage",
      diagnostics_stage = as.character(profile$profile_stage %||% "schema_only"),
      planner_contract = as.character(
        profile$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      branch_available = FALSE,
      reason = as.character(
        profile$reason %||%
          "No schema-only future-branch coverage/connectivity diagnostics are currently available."
      ),
      metric_registry = metric_registry,
      diagnostic_summary_table = tibble::tibble(),
      diagnostic_table = tibble::tibble(),
      active_branch_profile = profile,
      note = paste(
        "Deterministic coverage/connectivity diagnostics are unavailable until",
        "the active future-branch scaffold can be materialized."
      )
    ))
  }

  diagnostic_table <- profile_table |>
    dplyr::mutate(
      person_criterion_cells = as.numeric(.data$n_person) * as.numeric(.data$n_criterion),
      criterion_replications_per_person = as.numeric(.data$raters_per_person),
      available_rater_pairs = choose(as.numeric(.data$n_rater), 2),
      rater_pair_overlap_per_cell = choose(as.numeric(.data$raters_per_person), 2),
      total_rater_pair_overlaps = as.numeric(.data$person_criterion_cells) * as.numeric(.data$rater_pair_overlap_per_cell),
      pair_coverage_fraction_per_cell = dplyr::if_else(
        .data$available_rater_pairs > 0,
        as.numeric(.data$rater_pair_overlap_per_cell) / as.numeric(.data$available_rater_pairs),
        0
      ),
      redundant_scoring = as.integer(as.numeric(.data$raters_per_person) > 1)
    ) |>
    dplyr::select(-"available_rater_pairs")

  recommended_row <- diagnostic_table[diagnostic_table$recommended %in% TRUE, , drop = FALSE]
  diagnostic_summary_table <- metric_registry |>
    dplyr::mutate(
      min = vapply(.data$metric, function(metric) {
        min(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$metric, function(metric) {
        max(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$metric, function(metric) {
        mean(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$metric, function(metric) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[metric]][[1]])
      }, numeric(1))
    )

  list(
    diagnostics_contract = "arbitrary_facet_planning_active_branch_coverage",
    diagnostics_stage = as.character(profile$profile_stage %||% "schema_only"),
    planner_contract = as.character(
      profile$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = TRUE,
    reason = as.character(
      profile$reason %||%
        "Schema-only future-branch coverage/connectivity diagnostics are materialized."
    ),
    metric_registry = metric_registry,
    diagnostic_summary_table = diagnostic_summary_table,
    diagnostic_table = diagnostic_table,
    active_branch_profile = profile,
    note = paste(
      "Deterministic coverage/connectivity diagnostics for the active future-branch",
      "scaffold. These quantities summarize scored-cell counts and rater-pair",
      "overlap identities without implying psychometric performance."
    )
  )
}

simulation_future_branch_active_branch_coverage_overview <- function(x,
                                                                     design = NULL,
                                                                     prefer = NULL,
                                                                     view = c("public", "canonical", "branch"),
                                                                     mode = c("brief", "snapshot", "operation"),
                                                                     surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                     table_component = c("grid", "report", "surface_index"),
                                                                     x_var = NULL,
                                                                     group_var = NULL,
                                                                     id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_coverage_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  diagnostics <- simulation_future_branch_active_branch_coverage(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  diagnostic_table <- tibble::as_tibble(diagnostics$diagnostic_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    branch_available = isTRUE(diagnostics$branch_available),
    n_designs = nrow(diagnostic_table),
    n_metrics = if (is.data.frame(diagnostics$metric_registry)) nrow(diagnostics$metric_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(diagnostic_table) &&
                                any(diagnostic_table$recommended %in% TRUE)) {
      as.character(diagnostic_table$design_id[which(diagnostic_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_redundant_designs = if ("redundant_scoring" %in% names(diagnostic_table)) {
      sum(diagnostic_table$redundant_scoring %in% 1L, na.rm = TRUE)
    } else {
      0L
    },
    n_single_rater_designs = if ("redundant_scoring" %in% names(diagnostic_table)) {
      sum(diagnostic_table$redundant_scoring %in% 0L, na.rm = TRUE)
    } else {
      0L
    },
    n_pair_connected_designs = if ("rater_pair_overlap_per_cell" %in% names(diagnostic_table)) {
      sum(diagnostic_table$rater_pair_overlap_per_cell > 0, na.rm = TRUE)
    } else {
      0L
    },
    n_zero_pair_overlap_designs = if ("rater_pair_overlap_per_cell" %in% names(diagnostic_table)) {
      sum(diagnostic_table$rater_pair_overlap_per_cell <= 0, na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_active_branch_coverage_overview",
    overview_stage = as.character(diagnostics$diagnostics_stage %||% "schema_only"),
    planner_contract = as.character(
      diagnostics$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = isTRUE(diagnostics$branch_available),
    reason = as.character(
      diagnostics$reason %||%
        "No future-branch coverage/connectivity overview is currently available."
    ),
    headline_table = headline_table,
    metric_registry = tibble::as_tibble(diagnostics$metric_registry %||% tibble::tibble()),
    diagnostic_summary_table = tibble::as_tibble(diagnostics$diagnostic_summary_table %||% tibble::tibble()),
    active_branch_coverage = diagnostics,
    note = paste(
      "Compact coverage/connectivity overview for the active future-branch",
      "scaffold, combining a one-row redundancy headline with basis-aware",
      "deterministic rater-overlap diagnostics."
    )
  )
}

simulation_future_branch_active_branch_guardrail_registry <- function() {
  tibble::tibble(
    guardrail = c(
      "rater_linking_regime",
      "pair_coverage_regime",
      "integer_balance_regime",
      "redundancy_regime"
    ),
    basis_class = c(
      "exact_classification",
      "exact_classification",
      "exact_classification",
      "exact_classification"
    ),
    rule = c(
      "ifelse(raters_per_person <= 1, 'single_rater', ifelse(raters_per_person >= n_rater, 'fully_crossed', 'partial_overlap'))",
      "ifelse(pair_coverage_fraction_per_cell <= 0, 'no_pair_coverage', ifelse(pair_coverage_fraction_per_cell >= 1, 'full_pair_coverage', 'partial_pair_coverage'))",
      "ifelse(perfect_integer_observation_balance == 1, 'integer_balanced', 'integer_unbalanced')",
      "ifelse(redundant_scoring == 1, 'redundant', 'single_rater_only')"
    ),
    interpretation = c(
      "Exact cell-level linking regime implied by the current assignments-per-person relative to the available rater count.",
      "Exact pair-coverage regime implied by the current within-cell rater overlap fraction.",
      "Exact integer-balance regime induced by the current total observation count and rater count.",
      "Exact indicator of whether each scored cell uses more than one rater."
    ),
    psychometric = FALSE
  )
}

simulation_future_branch_active_branch_guardrails <- function(x,
                                                              design = NULL,
                                                              prefer = NULL,
                                                              view = c("public", "canonical", "branch"),
                                                              mode = c("brief", "snapshot", "operation"),
                                                              surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                              table_component = c("grid", "report", "surface_index"),
                                                              x_var = NULL,
                                                              group_var = NULL,
                                                              id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_guardrails",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  coverage <- simulation_future_branch_active_branch_coverage(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  load_balance <- simulation_future_branch_active_branch_load_balance(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  coverage_table <- tibble::as_tibble(coverage$diagnostic_table %||% tibble::tibble())
  load_balance_table <- tibble::as_tibble(load_balance$diagnostic_table %||% tibble::tibble())
  guardrail_registry <- simulation_future_branch_active_branch_guardrail_registry()

  if (!isTRUE(coverage$branch_available) || nrow(coverage_table) == 0L ||
      !isTRUE(load_balance$branch_available) || nrow(load_balance_table) == 0L) {
    return(list(
      guardrail_contract = "arbitrary_facet_planning_active_branch_guardrails",
      guardrail_stage = as.character(
        coverage$diagnostics_stage %||% load_balance$diagnostics_stage %||% "schema_only"
      ),
      planner_contract = as.character(
        coverage$planner_contract %||%
          load_balance$planner_contract %||%
          "arbitrary_facet_planning_scaffold"
      ),
      branch_available = FALSE,
      reason = as.character(
        coverage$reason %||%
          load_balance$reason %||%
          "No schema-only future-branch guardrail classifications are currently available."
      ),
      guardrail_registry = guardrail_registry,
      guardrail_summary_table = tibble::tibble(),
      guardrail_table = tibble::tibble(),
      active_branch_coverage = coverage,
      active_branch_load_balance = load_balance,
      note = paste(
        "Deterministic guardrail classifications are unavailable until the",
        "active future-branch scaffold can be materialized."
      )
    ))
  }

  join_cols <- c("design_id", "recommended")
  joined <- dplyr::left_join(
    coverage_table,
    load_balance_table[, c(join_cols, "perfect_integer_observation_balance"), drop = FALSE],
    by = join_cols
  )

  guardrail_table <- joined |>
    dplyr::mutate(
      rater_linking_regime = dplyr::case_when(
        as.numeric(.data$criterion_replications_per_person) <= 1 ~ "single_rater",
        as.numeric(.data$criterion_replications_per_person) >= as.numeric(.data$n_rater) ~ "fully_crossed",
        TRUE ~ "partial_overlap"
      ),
      pair_coverage_regime = dplyr::case_when(
        as.numeric(.data$pair_coverage_fraction_per_cell) <= 0 ~ "no_pair_coverage",
        as.numeric(.data$pair_coverage_fraction_per_cell) >= 1 ~ "full_pair_coverage",
        TRUE ~ "partial_pair_coverage"
      ),
      integer_balance_regime = dplyr::if_else(
        as.integer(.data$perfect_integer_observation_balance) == 1L,
        "integer_balanced",
        "integer_unbalanced"
      ),
      redundancy_regime = dplyr::if_else(
        as.integer(.data$redundant_scoring) == 1L,
        "redundant",
        "single_rater_only"
      )
    )

  guardrail_summary_table <- guardrail_registry |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_levels = dplyr::n_distinct(guardrail_table[[.data$guardrail]]),
      recommended_level = {
        row <- guardrail_table[guardrail_table$recommended %in% TRUE, , drop = FALSE]
        if (nrow(row) == 0L) NA_character_ else as.character(row[[.data$guardrail]][[1]])
      },
      observed_levels = paste(sort(unique(as.character(guardrail_table[[.data$guardrail]]))), collapse = ", ")
    ) |>
    dplyr::ungroup()

  list(
    guardrail_contract = "arbitrary_facet_planning_active_branch_guardrails",
    guardrail_stage = as.character(
      coverage$diagnostics_stage %||% load_balance$diagnostics_stage %||% "schema_only"
    ),
    planner_contract = as.character(
      coverage$planner_contract %||%
        load_balance$planner_contract %||%
        "arbitrary_facet_planning_scaffold"
    ),
    branch_available = TRUE,
    reason = "Schema-only future-branch guardrail classifications are materialized.",
    guardrail_registry = guardrail_registry,
    guardrail_summary_table = guardrail_summary_table,
    guardrail_table = guardrail_table,
    active_branch_coverage = coverage,
    active_branch_load_balance = load_balance,
    note = paste(
      "Deterministic guardrail classifications for the active future-branch",
      "scaffold. These are exact design-regime labels derived from overlap and",
      "integer-balance structure; they do not estimate psychometric performance."
    )
  )
}

simulation_future_branch_active_branch_guardrail_overview <- function(x,
                                                                      design = NULL,
                                                                      prefer = NULL,
                                                                      view = c("public", "canonical", "branch"),
                                                                      mode = c("brief", "snapshot", "operation"),
                                                                      surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                      table_component = c("grid", "report", "surface_index"),
                                                                      x_var = NULL,
                                                                      group_var = NULL,
                                                                      id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_guardrail_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  guardrails <- simulation_future_branch_active_branch_guardrails(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  guardrail_table <- tibble::as_tibble(guardrails$guardrail_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    branch_available = isTRUE(guardrails$branch_available),
    n_designs = nrow(guardrail_table),
    n_guardrails = if (is.data.frame(guardrails$guardrail_registry)) nrow(guardrails$guardrail_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(guardrail_table) &&
                                any(guardrail_table$recommended %in% TRUE)) {
      as.character(guardrail_table$design_id[which(guardrail_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_single_rater_designs = if ("rater_linking_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$rater_linking_regime %in% "single_rater", na.rm = TRUE)
    } else {
      0L
    },
    n_partial_overlap_designs = if ("rater_linking_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$rater_linking_regime %in% "partial_overlap", na.rm = TRUE)
    } else {
      0L
    },
    n_fully_crossed_designs = if ("rater_linking_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$rater_linking_regime %in% "fully_crossed", na.rm = TRUE)
    } else {
      0L
    },
    n_integer_balanced_designs = if ("integer_balance_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$integer_balance_regime %in% "integer_balanced", na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_active_branch_guardrail_overview",
    overview_stage = as.character(guardrails$guardrail_stage %||% "schema_only"),
    planner_contract = as.character(
      guardrails$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = isTRUE(guardrails$branch_available),
    reason = as.character(
      guardrails$reason %||%
        "No future-branch guardrail overview is currently available."
    ),
    headline_table = headline_table,
    guardrail_registry = tibble::as_tibble(guardrails$guardrail_registry %||% tibble::tibble()),
    guardrail_summary_table = tibble::as_tibble(guardrails$guardrail_summary_table %||% tibble::tibble()),
    active_branch_guardrails = guardrails,
    note = paste(
      "Compact guardrail overview for the active future-branch scaffold,",
      "combining exact design-regime counts with a basis-aware guardrail index."
    )
  )
}

simulation_future_branch_active_branch_readiness_registry <- function() {
  tibble::tibble(
    indicator = c(
      "supports_multi_rater_cells",
      "supports_pair_overlap",
      "supports_integer_balanced_load",
      "supports_full_pair_coverage"
    ),
    basis_class = c(
      "exact_indicator",
      "exact_indicator",
      "exact_indicator",
      "exact_indicator"
    ),
    rule = c(
      "as.integer(redundancy_regime == 'redundant')",
      "as.integer(pair_coverage_regime != 'no_pair_coverage')",
      "as.integer(integer_balance_regime == 'integer_balanced')",
      "as.integer(pair_coverage_regime == 'full_pair_coverage')"
    ),
    interpretation = c(
      "Exact indicator that each scored cell uses more than one rater.",
      "Exact indicator that the current design yields at least one within-cell rater pair overlap.",
      "Exact indicator that the total observation count is evenly divisible across raters.",
      "Exact indicator that every available rater pair is represented within each scored cell."
    ),
    psychometric = FALSE
  )
}

simulation_future_branch_active_branch_readiness <- function(x,
                                                             design = NULL,
                                                             prefer = NULL,
                                                             view = c("public", "canonical", "branch"),
                                                             mode = c("brief", "snapshot", "operation"),
                                                             surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                             table_component = c("grid", "report", "surface_index"),
                                                             x_var = NULL,
                                                             group_var = NULL,
                                                             id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_readiness",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  guardrails <- simulation_future_branch_active_branch_guardrails(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  guardrail_table <- tibble::as_tibble(guardrails$guardrail_table %||% tibble::tibble())
  indicator_registry <- simulation_future_branch_active_branch_readiness_registry()

  if (!isTRUE(guardrails$branch_available) || nrow(guardrail_table) == 0L) {
    return(list(
      readiness_contract = "arbitrary_facet_planning_active_branch_readiness",
      readiness_stage = as.character(guardrails$guardrail_stage %||% "schema_only"),
      planner_contract = as.character(
        guardrails$planner_contract %||% "arbitrary_facet_planning_scaffold"
      ),
      branch_available = FALSE,
      reason = as.character(
        guardrails$reason %||%
          "No schema-only future-branch structural readiness summary is currently available."
      ),
      indicator_registry = indicator_registry,
      readiness_summary_table = tibble::tibble(),
      readiness_table = tibble::tibble(),
      active_branch_guardrails = guardrails,
      note = paste(
        "Deterministic structural readiness indicators are unavailable until",
        "the active future-branch scaffold can be materialized."
      )
    ))
  }

  readiness_table <- guardrail_table |>
    dplyr::mutate(
      supports_multi_rater_cells = as.integer(.data$redundancy_regime == "redundant"),
      supports_pair_overlap = as.integer(.data$pair_coverage_regime != "no_pair_coverage"),
      supports_integer_balanced_load = as.integer(.data$integer_balance_regime == "integer_balanced"),
      supports_full_pair_coverage = as.integer(.data$pair_coverage_regime == "full_pair_coverage"),
      structural_tier = dplyr::case_when(
        .data$supports_pair_overlap <= 0 ~ "single_rater_only",
        .data$supports_full_pair_coverage >= 1 & .data$supports_integer_balanced_load >= 1 ~ "full_overlap_balanced",
        .data$supports_full_pair_coverage >= 1 ~ "full_overlap_unbalanced",
        .data$supports_integer_balanced_load >= 1 ~ "partial_overlap_balanced",
        TRUE ~ "partial_overlap_unbalanced"
      )
    )

  recommended_row <- readiness_table[readiness_table$recommended %in% TRUE, , drop = FALSE]
  readiness_summary_table <- indicator_registry |>
    dplyr::mutate(
      min = vapply(.data$indicator, function(indicator) {
        min(readiness_table[[indicator]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$indicator, function(indicator) {
        max(readiness_table[[indicator]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$indicator, function(indicator) {
        mean(readiness_table[[indicator]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$indicator, function(indicator) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[indicator]][[1]])
      }, numeric(1))
    )

  list(
    readiness_contract = "arbitrary_facet_planning_active_branch_readiness",
    readiness_stage = as.character(guardrails$guardrail_stage %||% "schema_only"),
    planner_contract = as.character(
      guardrails$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = TRUE,
    reason = "Schema-only future-branch structural readiness summary is materialized.",
    indicator_registry = indicator_registry,
    readiness_summary_table = readiness_summary_table,
    readiness_table = readiness_table,
    active_branch_guardrails = guardrails,
    note = paste(
      "Deterministic structural readiness indicators for the active future-branch",
      "scaffold. These summarize exact overlap and balance preconditions only;",
      "they do not estimate psychometric performance."
    )
  )
}

simulation_future_branch_active_branch_readiness_overview <- function(x,
                                                                      design = NULL,
                                                                      prefer = NULL,
                                                                      view = c("public", "canonical", "branch"),
                                                                      mode = c("brief", "snapshot", "operation"),
                                                                      surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                      table_component = c("grid", "report", "surface_index"),
                                                                      x_var = NULL,
                                                                      group_var = NULL,
                                                                      id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_readiness_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  readiness <- simulation_future_branch_active_branch_readiness(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  readiness_table <- tibble::as_tibble(readiness$readiness_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    branch_available = isTRUE(readiness$branch_available),
    n_designs = nrow(readiness_table),
    n_indicators = if (is.data.frame(readiness$indicator_registry)) nrow(readiness$indicator_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(readiness_table) &&
                                any(readiness_table$recommended %in% TRUE)) {
      as.character(readiness_table$design_id[which(readiness_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_single_rater_only_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% "single_rater_only", na.rm = TRUE)
    } else {
      0L
    },
    n_partial_overlap_balanced_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% "partial_overlap_balanced", na.rm = TRUE)
    } else {
      0L
    },
    n_partial_overlap_unbalanced_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% "partial_overlap_unbalanced", na.rm = TRUE)
    } else {
      0L
    },
    n_full_overlap_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% c("full_overlap_balanced", "full_overlap_unbalanced"), na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_active_branch_readiness_overview",
    overview_stage = as.character(readiness$readiness_stage %||% "schema_only"),
    planner_contract = as.character(
      readiness$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    branch_available = isTRUE(readiness$branch_available),
    reason = as.character(
      readiness$reason %||%
        "No future-branch structural readiness overview is currently available."
    ),
    headline_table = headline_table,
    indicator_registry = tibble::as_tibble(readiness$indicator_registry %||% tibble::tibble()),
    readiness_summary_table = tibble::as_tibble(readiness$readiness_summary_table %||% tibble::tibble()),
    active_branch_readiness = readiness,
    note = paste(
      "Compact structural readiness overview for the active future-branch",
      "scaffold, combining exact overlap/balance tiers with a basis-aware",
      "indicator index."
    )
  )
}

simulation_future_branch_active_branch_recommendation <- function(x,
                                                                  design = NULL,
                                                                  prefer = NULL,
                                                                  view = c("public", "canonical", "branch"),
                                                                  mode = c("brief", "snapshot", "operation"),
                                                                  surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                  table_component = c("grid", "report", "surface_index"),
                                                                  x_var = NULL,
                                                                  group_var = NULL,
                                                                  id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_recommendation",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  readiness <- simulation_future_branch_active_branch_readiness(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  profile <- simulation_future_branch_active_branch_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  readiness_table <- tibble::as_tibble(readiness$readiness_table %||% tibble::tibble())
  profile_table <- tibble::as_tibble(profile$profile_table %||% tibble::tibble())

  if (!isTRUE(readiness$branch_available) || nrow(readiness_table) == 0L ||
      !isTRUE(profile$branch_available) || nrow(profile_table) == 0L) {
    return(list(
      recommendation_contract = "arbitrary_facet_planning_active_branch_recommendation",
      recommendation_stage = as.character(
        readiness$readiness_stage %||% profile$profile_stage %||% "schema_only"
      ),
      planner_contract = as.character(
        readiness$planner_contract %||%
          profile$planner_contract %||%
          "arbitrary_facet_planning_scaffold"
      ),
      recommendation_available = FALSE,
      reason = as.character(
        readiness$reason %||%
          profile$reason %||%
          "No schema-only future-branch structural recommendation is currently available."
      ),
      ranking_rule = paste(
        "Unavailable because the active future-branch readiness and profile",
        "contracts are not both materialized."
      ),
      recommended_design_id = character(0),
      recommended_tier = character(0),
      recommendation_table = tibble::tibble(),
      active_branch_readiness = readiness,
      active_branch_profile = profile,
      note = paste(
        "Structural recommendation is unavailable until the active future-branch",
        "scaffold can be materialized."
      )
    ))
  }

  tier_priority <- c(
    full_overlap_balanced = 1L,
    full_overlap_unbalanced = 2L,
    partial_overlap_balanced = 3L,
    partial_overlap_unbalanced = 4L,
    single_rater_only = 5L
  )

  recommendation_table <- readiness_table
  required_profile_cols <- c(
    "total_observations", "n_person", "n_rater", "n_criterion", "raters_per_person"
  )
  missing_profile_cols <- setdiff(required_profile_cols, names(recommendation_table))
  if (length(missing_profile_cols) > 0L) {
    recommendation_table <- dplyr::left_join(
      recommendation_table,
      profile_table[, c("design_id", missing_profile_cols), drop = FALSE],
      by = "design_id"
    )
  }

  recommendation_table <- recommendation_table |>
    dplyr::mutate(
      structural_priority = unname(tier_priority[as.character(.data$structural_tier)]),
      structural_priority = dplyr::coalesce(.data$structural_priority, length(tier_priority) + 1L)
    ) |>
    dplyr::arrange(
      .data$structural_priority,
      .data$total_observations,
      .data$n_person,
      .data$n_criterion,
      .data$raters_per_person,
      .data$n_rater,
      .data$design_id
    ) |>
    dplyr::mutate(recommended = dplyr::row_number() == 1L)

  recommended_row <- dplyr::slice_head(recommendation_table, n = 1)
  recommended_id <- as.character(recommended_row$design_id[[1]] %||% NA_character_)
  recommended_tier <- as.character(recommended_row$structural_tier[[1]] %||% NA_character_)

  list(
    recommendation_contract = "arbitrary_facet_planning_active_branch_recommendation",
    recommendation_stage = as.character(
      readiness$readiness_stage %||% profile$profile_stage %||% "schema_only"
    ),
    planner_contract = as.character(
      readiness$planner_contract %||%
        profile$planner_contract %||%
        "arbitrary_facet_planning_scaffold"
    ),
    recommendation_available = TRUE,
    reason = paste(
      "Structural recommendation is available from the active future-branch",
      "readiness and profile contracts."
    ),
    ranking_rule = paste(
      "Rank by structural tier in the order",
      "full_overlap_balanced > full_overlap_unbalanced > partial_overlap_balanced >",
      "partial_overlap_unbalanced > single_rater_only, then minimize total_observations",
      "and remaining canonical count variables."
    ),
    recommended_design_id = recommended_id,
    recommended_tier = recommended_tier,
    recommendation_table = recommendation_table,
    active_branch_readiness = readiness,
    active_branch_profile = profile,
    note = paste(
      "Conservative structural recommendation for the active future-branch",
      "scaffold. This ranking uses exact overlap/balance tiers and total",
      "observation counts only; it is not a psychometric optimization."
    )
  )
}

simulation_future_branch_active_branch_recommendation_overview <- function(x,
                                                                           design = NULL,
                                                                           prefer = NULL,
                                                                           view = c("public", "canonical", "branch"),
                                                                           mode = c("brief", "snapshot", "operation"),
                                                                           surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                           table_component = c("grid", "report", "surface_index"),
                                                                           x_var = NULL,
                                                                           group_var = NULL,
                                                                           id_prefix = NULL) {
  cached <- simulation_future_branch_cached_default(
    x = x,
    field = "active_branch_recommendation_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  recommendation <- simulation_future_branch_active_branch_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  recommendation_table <- tibble::as_tibble(recommendation$recommendation_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    recommendation_available = isTRUE(recommendation$recommendation_available),
    n_designs = nrow(recommendation_table),
    recommended_design_id = as.character(recommendation$recommended_design_id %||% NA_character_),
    recommended_tier = as.character(recommendation$recommended_tier %||% NA_character_),
    n_top_tier_candidates = if ("structural_priority" %in% names(recommendation_table)) {
      min_priority <- min(recommendation_table$structural_priority, na.rm = TRUE)
      sum(recommendation_table$structural_priority == min_priority, na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_active_branch_recommendation_overview",
    overview_stage = as.character(recommendation$recommendation_stage %||% "schema_only"),
    planner_contract = as.character(
      recommendation$planner_contract %||% "arbitrary_facet_planning_scaffold"
    ),
    recommendation_available = isTRUE(recommendation$recommendation_available),
    reason = as.character(
      recommendation$reason %||%
        "No future-branch structural recommendation overview is currently available."
    ),
    headline_table = headline_table,
    recommendation_table = recommendation_table,
    active_branch_recommendation = recommendation,
    note = paste(
      "Compact structural recommendation overview for the active future-branch",
      "scaffold, summarizing the conservative exact-tier ranking outcome."
    )
  )
}

future_branch_active_table_index <- function(active,
                                             overview,
                                             load_balance,
                                             coverage,
                                             guardrails,
                                             readiness,
                                             recommendation) {
  rows_or_zero <- function(tbl) {
    if (!is.data.frame(tbl)) return(0L)
    as.integer(nrow(tbl))
  }

  tibble::tibble(
    Table = c(
      "overview",
      "profile_summary",
      "load_balance_summary",
      "coverage_summary",
      "guardrail_summary",
      "readiness_summary",
      "recommendation_table"
    ),
    Rows = c(
      rows_or_zero(overview$headline_table),
      rows_or_zero(overview$metric_summary_table),
      rows_or_zero(load_balance$diagnostic_summary_table),
      rows_or_zero(coverage$diagnostic_summary_table),
      rows_or_zero(guardrails$guardrail_summary_table),
      rows_or_zero(readiness$readiness_summary_table),
      rows_or_zero(recommendation$recommendation_table)
    ),
    Role = c(
      "overview",
      "profile",
      "load_balance",
      "coverage",
      "guardrails",
      "readiness",
      "recommendation"
    ),
    Description = c(
      "Headline active-branch configuration and selected branch surface.",
      "Deterministic observation/load profile for the current design grid.",
      "Balanced-expectation and integer split summaries.",
      "Coverage/connectivity summaries implied by the current design grid.",
      "Exact overlap/balance regime summaries.",
      "Exact structural readiness indicators and tier summaries.",
      "Conservative deterministic recommendation ranking."
    ),
    stringsAsFactors = FALSE
  )
}

future_branch_appendix_selection_tables <- function(bundle,
                                                    label = "future_branch_active_branch",
                                                    presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  presets <- unique(as.character(presets))
  original_bundles <- stats::setNames(list(bundle), label)

  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  selection_catalog <- bind_tables(lapply(presets, function(preset) {
    selected <- export_select_summary_table_bundles_for_appendix(
      original_bundles,
      preset = preset
    )
    export_summary_table_selection_catalog(
      original_bundles = original_bundles,
      selected_bundles = selected,
      preset = preset
    )
  }))

  list(
    selection_catalog = selection_catalog,
    selection_summary = bind_tables(lapply(presets, function(preset) {
      export_summary_table_selection_summary(selection_catalog, preset = preset)
    })),
    selection_role_summary = bind_tables(lapply(presets, function(preset) {
      export_summary_table_selection_role_summary(selection_catalog, preset = preset)
    })),
    selection_section_summary = bind_tables(lapply(presets, function(preset) {
      export_summary_table_selection_section_summary(selection_catalog, preset = preset)
    }))
  )
}

future_branch_selection_table_summary <- function(selection_catalog,
                                                  preset = NULL) {
  tbl <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    return(data.frame())
  }

  if (!is.null(preset)) {
    tbl <- tbl[as.character(tbl$Preset %||% "") %in% unique(as.character(preset)), , drop = FALSE]
  }
  tbl <- tbl[tbl$Selected %in% TRUE, , drop = FALSE]
  if (nrow(tbl) == 0L || !"Table" %in% names(tbl)) {
    return(data.frame())
  }

  compact_unique <- function(x, max_n = 4L) {
    summary_table_bundle_compact_labels(unique(as.character(x %||% character(0))), max_n = max_n)
  }
  first_non_missing <- function(x, default = NA_character_) {
    x <- as.character(x %||% character(0))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
      return(default)
    }
    x[[1]]
  }
  first_numeric <- function(x, default = NA_real_) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(default)
    }
    x[[1]]
  }

  split_tbl <- split(tbl, as.character(tbl$Table))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(table_nm) {
      part <- split_tbl[[table_nm]]
      data.frame(
        Table = as.character(table_nm),
        PresetsSelected = length(unique(as.character(part$Preset %||% ""))),
        Presets = compact_unique(part$Preset, max_n = 7L),
        Rows = as.integer(first_numeric(part$Rows, default = 0)),
        Role = first_non_missing(part$Role),
        AppendixSection = first_non_missing(part$AppendixSection),
        PreferredAppendixOrder = as.integer(first_numeric(part$PreferredAppendixOrder, default = 9999)),
        PlotReady = any(part$PlotReady %in% TRUE, na.rm = TRUE),
        ExportReady = any(part$ExportReady %in% TRUE, na.rm = TRUE),
        ApaTableReady = any(part$ApaTableReady %in% TRUE, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$PreferredAppendixOrder, out$Table, na.last = TRUE), , drop = FALSE]
}

future_branch_selection_table_preset_summary <- function(selection_catalog,
                                                         presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_table_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

future_branch_selection_handoff_table_summary <- function(selection_catalog,
                                                          presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_table_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

future_branch_selection_handoff_summary <- function(selection_catalog,
                                                    presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  tbl <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    return(data.frame())
  }

  presets <- unique(as.character(presets))
  tbl <- tbl[tbl$Selected %in% TRUE & as.character(tbl$Preset %||% "") %in% presets, , drop = FALSE]
  if (nrow(tbl) == 0L || !"AppendixSection" %in% names(tbl)) {
    return(data.frame())
  }

  compact_unique <- function(x, max_n = 4L) {
    summary_table_bundle_compact_labels(unique(as.character(x %||% character(0))), max_n = max_n)
  }

  split_tbl <- split(tbl, paste(as.character(tbl$Preset), as.character(tbl$AppendixSection), sep = "\r"))
  out <- do.call(
    rbind,
    lapply(split_tbl, function(part) {
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = as.character(part$Preset[[1]] %||% ""),
        AppendixSection = as.character(part$AppendixSection[[1]] %||% ""),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        RolesCovered = compact_unique(part$Role, max_n = 4L),
        KeyTables = compact_unique(part$Table, max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$AppendixSection), , drop = FALSE]
}

future_branch_selection_handoff_bundle_summary <- function(selection_catalog,
                                                           presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_bundle_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

future_branch_selection_handoff_preset_summary <- function(selection_catalog,
                                                           presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_preset_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

future_branch_selection_handoff_role_summary <- function(selection_catalog,
                                                         presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_role_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

future_branch_selection_handoff_role_section_summary <- function(selection_catalog,
                                                                 presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_role_section_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

future_branch_active_plot_index_from_bundle <- function(plot_index) {
  idx <- as.data.frame(plot_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(idx) == 0L) {
    return(idx)
  }
  if (!all(c("Table", "DefaultPlotTypes") %in% names(idx))) {
    return(idx)
  }

  direct_routes <- c(
    future_branch_profile = "profile_metrics",
    future_branch_load_balance = "load_balance",
    future_branch_coverage = "coverage",
    future_branch_readiness = "readiness_tiers",
    future_branch_appendix_roles = "appendix_roles",
    future_branch_appendix_sections = "appendix_sections",
    future_branch_appendix_presets = "appendix_presets",
    future_branch_selection_table_presets = "selection_tables",
    future_branch_selection_handoff_presets = "selection_handoff_presets",
    future_branch_selection_handoff = "selection_handoff",
    future_branch_selection_handoff_bundles = "selection_handoff_bundles",
    future_branch_selection_handoff_roles = "selection_handoff_roles",
    future_branch_selection_handoff_role_sections = "selection_handoff_role_sections",
    future_branch_selection_tables = "selection_tables",
    future_branch_selection_summary = "selection_bundles",
    future_branch_selection_roles = "selection_roles",
    future_branch_selection_sections = "selection_sections"
  )

  idx$DefaultPlotTypes <- vapply(seq_len(nrow(idx)), function(i) {
    tbl <- as.character(idx$Table[i] %||% "")
    base_types <- strsplit(as.character(idx$DefaultPlotTypes[i] %||% ""), ",", fixed = TRUE)[[1]]
    base_types <- trimws(base_types)
    base_types <- base_types[nzchar(base_types)]
    extra <- unname(direct_routes[tbl])
    extra <- extra[!is.na(extra) & nzchar(extra)]
    paste(unique(c(extra, base_types)), collapse = ", ")
  }, character(1))

  idx
}

#' Summarize a future arbitrary-facet planning active branch
#'
#' @param object Output from the future-branch active planning scaffold stored
#'   in `planning_schema$future_branch_active_branch`.
#' @param digits Number of digits used in numeric summaries.
#' @param top_n Maximum number of recommendation rows to print in the preview.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is intentionally conservative. It aggregates only deterministic
#' branch-side quantities already validated in the schema-first arbitrary-facet
#' planning scaffold: observation bookkeeping, load/balance, coverage,
#' guardrails, structural readiness, and conservative recommendation ranking.
#' It also exposes the same manuscript-facing table/appendix metadata used by
#' [build_summary_table_bundle()] so the future branch can be reviewed directly
#' without first routing through planning summaries. In addition to bundle-level
#' appendix presets and section counts, it includes export-like appendix
#' selection summaries by preset, reporting role, manuscript section,
#' bundle-aware handoff summaries, preset-specific table surface, and a
#' table-level handoff crosswalk, plus direct `role_summary` / `table_profile`
#' surfaces for table-shape review.
#' It does not report psychometric recovery or Monte Carlo performance.
#'
#' @return An object of class `summary.mfrm_future_branch_active_branch`.
#' @seealso [summary.mfrm_design_evaluation()], [plot.mfrm_future_branch_active_branch()]
#' @export
summary.mfrm_future_branch_active_branch <- function(object, digits = 3, top_n = 8, ...) {
  active <- simulation_future_branch_active_branch(object)
  overview <- simulation_future_branch_active_branch_overview(active)
  load_balance <- simulation_future_branch_active_branch_load_balance_overview(active)
  coverage <- simulation_future_branch_active_branch_coverage_overview(active)
  guardrails <- simulation_future_branch_active_branch_guardrail_overview(active)
  readiness <- simulation_future_branch_active_branch_readiness_overview(active)
  recommendation <- simulation_future_branch_active_branch_recommendation_overview(active)

  digits <- max(0L, as.integer(digits[1]))
  top_n <- max(1L, as.integer(top_n[1]))

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  headline <- tibble::tibble(
    branch_available = isTRUE(active$branch_available),
    n_designs = as.integer(active$n_designs %||% 0L),
    recommended_design_id = as.character(active$recommended_design_id %||% NA_character_),
    view = as.character(active$view %||% NA_character_),
    mode = as.character(active$mode %||% NA_character_),
    surface = as.character(active$surface %||% NA_character_),
    table_component = as.character(active$table_component %||% NA_character_)
  )

  recommendation_table <- tibble::as_tibble(recommendation$recommendation_table %||% tibble::tibble())
  if (nrow(recommendation_table) > 0L) {
    recommendation_table <- utils::head(recommendation_table, n = top_n)
  }

  notes <- unique(stats::na.omit(c(
    as.character(active$note %||% character(0)),
    as.character(overview$note %||% character(0)),
    as.character(load_balance$note %||% character(0)),
    as.character(coverage$note %||% character(0)),
    as.character(guardrails$note %||% character(0)),
    as.character(readiness$note %||% character(0)),
    as.character(recommendation$note %||% character(0))
  )))

  out <- list(
    overview = round_df(headline),
    table_index = future_branch_active_table_index(
      active = active,
      overview = overview,
      load_balance = load_balance,
      coverage = coverage,
      guardrails = guardrails,
      readiness = readiness,
      recommendation = recommendation
    ),
    profile_summary = round_df(tibble::as_tibble(overview$metric_summary_table %||% tibble::tibble())),
    load_balance_summary = round_df(tibble::as_tibble(load_balance$diagnostic_summary_table %||% tibble::tibble())),
    coverage_summary = round_df(tibble::as_tibble(coverage$diagnostic_summary_table %||% tibble::tibble())),
    guardrail_summary = round_df(tibble::as_tibble(guardrails$guardrail_summary_table %||% tibble::tibble())),
    readiness_summary = round_df(tibble::as_tibble(readiness$readiness_summary_table %||% tibble::tibble())),
    recommendation_table = round_df(recommendation_table),
    notes = notes,
    digits = digits
  )
  temp_core <- out
  class(temp_core) <- "summary.mfrm_future_branch_active_branch"
  bundle_core <- build_summary_table_bundle(temp_core, include_empty = TRUE)
  selection_tables <- future_branch_appendix_selection_tables(bundle_core)
  out$selection_table_summary <- future_branch_selection_table_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_table_preset_summary <- future_branch_selection_table_preset_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_table_summary <- future_branch_selection_handoff_table_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_preset_summary <- future_branch_selection_handoff_preset_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_summary <- future_branch_selection_handoff_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_bundle_summary <- future_branch_selection_handoff_bundle_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_role_summary <- future_branch_selection_handoff_role_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_role_section_summary <- future_branch_selection_handoff_role_section_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_summary <- selection_tables$selection_summary
  out$selection_role_summary <- selection_tables$selection_role_summary
  out$selection_section_summary <- selection_tables$selection_section_summary
  out$selection_catalog <- selection_tables$selection_catalog

  temp_final <- out
  class(temp_final) <- "summary.mfrm_future_branch_active_branch"
  bundle_final <- build_summary_table_bundle(temp_final, include_empty = TRUE)
  out$table_index <- as.data.frame(bundle_final$table_index %||% data.frame(), stringsAsFactors = FALSE)
  out$plot_index <- future_branch_active_plot_index_from_bundle(
    as.data.frame(bundle_final$plot_index %||% data.frame(), stringsAsFactors = FALSE)
  )
  out$table_catalog <- summary_table_bundle_catalog(bundle_final)
  out$table_profile <- summary_table_bundle_profile(bundle_final)
  if (nrow(out$table_profile) > 0L) {
    ord <- order(out$table_profile$Rows, out$table_profile$Cols, decreasing = TRUE, na.last = TRUE)
    out$table_profile <- out$table_profile[ord, , drop = FALSE]
    out$table_profile <- utils::head(out$table_profile, n = top_n)
  }
  out$role_summary <- data.frame()
  if (nrow(out$table_index) > 0L && "Role" %in% names(out$table_index)) {
    roles <- split(out$table_index, out$table_index$Role %||% "")
    out$role_summary <- do.call(
      rbind,
      lapply(names(roles), function(role_nm) {
        part <- roles[[role_nm]]
        data.frame(
          Role = as.character(role_nm),
          Tables = nrow(part),
          TotalRows = sum(suppressWarnings(as.numeric(part$Rows)), na.rm = TRUE),
          TotalCols = sum(suppressWarnings(as.numeric(part$Cols)), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
    out$role_summary <- out$role_summary[
      order(out$role_summary$Tables, out$role_summary$Role, decreasing = TRUE),
      ,
      drop = FALSE
    ]
    rownames(out$role_summary) <- NULL
  }
  out$appendix_presets <- summary_table_bundle_appendix_presets(out$table_catalog)
  out$appendix_role_summary <- summary_table_bundle_appendix_role_summary(out$table_catalog)
  out$appendix_section_summary <- summary_table_bundle_appendix_section_summary(out$table_catalog)
  out$reporting_map <- summary_table_bundle_reporting_map(bundle_final, out$table_catalog)
  out$overview$RecommendedAppendixTables <- sum(out$table_catalog$RecommendedAppendix %in% TRUE, na.rm = TRUE)
  out$overview$CompactAppendixTables <- sum(out$table_catalog$CompactAppendix %in% TRUE, na.rm = TRUE)
  out$overview$NumericTables <- sum(out$table_profile$NumericColumns > 0, na.rm = TRUE)
  out$overview$AnyNumericTable <- nrow(out$table_profile) > 0L &&
    any(out$table_profile$NumericColumns > 0, na.rm = TRUE)
  class(out) <- "summary.mfrm_future_branch_active_branch"
  out
}

#' @export
print.summary.mfrm_future_branch_active_branch <- function(x, ...) {
  digits <- max(0L, as.integer(x$digits %||% 3L))

  cat("mfrmr Future Arbitrary-Facet Planning Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$table_index) && nrow(x$table_index) > 0L) {
    cat("\nTable index\n")
    print(as.data.frame(x$table_index), row.names = FALSE)
  }
  if (!is.null(x$plot_index) && nrow(x$plot_index) > 0L) {
    cat("\nPlot index\n")
    print(as.data.frame(x$plot_index), row.names = FALSE)
  }
  if (!is.null(x$table_catalog) && nrow(x$table_catalog) > 0L) {
    cat("\nTable catalog\n")
    print(round_numeric_df(as.data.frame(x$table_catalog), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$role_summary) && nrow(x$role_summary) > 0L) {
    cat("\nRole summary\n")
    print(round_numeric_df(as.data.frame(x$role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$table_profile) && nrow(x$table_profile) > 0L) {
    cat("\nTable profile\n")
    print(round_numeric_df(as.data.frame(x$table_profile), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$profile_summary) && nrow(x$profile_summary) > 0L) {
    cat("\nProfile summary\n")
    print(round_numeric_df(as.data.frame(x$profile_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$load_balance_summary) && nrow(x$load_balance_summary) > 0L) {
    cat("\nLoad/balance summary\n")
    print(round_numeric_df(as.data.frame(x$load_balance_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$coverage_summary) && nrow(x$coverage_summary) > 0L) {
    cat("\nCoverage summary\n")
    print(round_numeric_df(as.data.frame(x$coverage_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$guardrail_summary) && nrow(x$guardrail_summary) > 0L) {
    cat("\nGuardrail summary\n")
    print(round_numeric_df(as.data.frame(x$guardrail_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$readiness_summary) && nrow(x$readiness_summary) > 0L) {
    cat("\nReadiness summary\n")
    print(round_numeric_df(as.data.frame(x$readiness_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$recommendation_table) && nrow(x$recommendation_table) > 0L) {
    cat("\nRecommendation table\n")
    print(round_numeric_df(as.data.frame(x$recommendation_table), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_presets) && nrow(x$appendix_presets) > 0L) {
    cat("\nAppendix presets\n")
    print(round_numeric_df(as.data.frame(x$appendix_presets), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_role_summary) && nrow(x$appendix_role_summary) > 0L) {
    cat("\nAppendix role summary\n")
    print(round_numeric_df(as.data.frame(x$appendix_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_section_summary) && nrow(x$appendix_section_summary) > 0L) {
    cat("\nAppendix section summary\n")
    print(round_numeric_df(as.data.frame(x$appendix_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_summary) && nrow(x$selection_summary) > 0L) {
    cat("\nSelection summary\n")
    print(round_numeric_df(as.data.frame(x$selection_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_table_summary) && nrow(x$selection_table_summary) > 0L) {
    cat("\nSelection table summary\n")
    print(round_numeric_df(as.data.frame(x$selection_table_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_table_preset_summary) && nrow(x$selection_table_preset_summary) > 0L) {
    cat("\nSelection table preset summary\n")
    print(round_numeric_df(as.data.frame(x$selection_table_preset_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_table_summary) && nrow(x$selection_handoff_table_summary) > 0L) {
    cat("\nSelection handoff table summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_table_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_preset_summary) && nrow(x$selection_handoff_preset_summary) > 0L) {
    cat("\nSelection handoff preset summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_preset_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_summary) && nrow(x$selection_handoff_summary) > 0L) {
    cat("\nSelection handoff summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_bundle_summary) && nrow(x$selection_handoff_bundle_summary) > 0L) {
    cat("\nSelection handoff bundle summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_bundle_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_role_summary) && nrow(x$selection_handoff_role_summary) > 0L) {
    cat("\nSelection handoff role summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_role_section_summary) && nrow(x$selection_handoff_role_section_summary) > 0L) {
    cat("\nSelection handoff role-section summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_role_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_role_summary) && nrow(x$selection_role_summary) > 0L) {
    cat("\nSelection role summary\n")
    print(round_numeric_df(as.data.frame(x$selection_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_section_summary) && nrow(x$selection_section_summary) > 0L) {
    cat("\nSelection section summary\n")
    print(round_numeric_df(as.data.frame(x$selection_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$reporting_map) && nrow(x$reporting_map) > 0L) {
    cat("\nReporting map\n")
    print(as.data.frame(x$reporting_map), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

simulation_compact_future_branch_active_summary <- function(x,
                                                            digits = 3,
                                                            top_n = 6L) {
  planning_schema <- simulation_object_planning_schema(x)
  active_branch <- planning_schema$future_branch_active_branch %||% NULL
  if (!is.list(active_branch)) {
    return(NULL)
  }

  out <- tryCatch(
    summary.mfrm_future_branch_active_branch(
      active_branch,
      digits = digits,
      top_n = top_n
    ),
    error = function(e) NULL
  )
  if (!inherits(out, "summary.mfrm_future_branch_active_branch")) {
    return(NULL)
  }

  out
}

print_compact_future_branch_active_summary <- function(x,
                                                       digits = 3,
                                                       heading = "Future arbitrary-facet planning scaffold") {
  if (!inherits(x, "summary.mfrm_future_branch_active_branch")) {
    return(invisible(NULL))
  }

  cat("\n", heading, "\n", sep = "")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$readiness_summary) && nrow(x$readiness_summary) > 0L) {
    cat("\nReadiness summary\n")
    print(round_numeric_df(as.data.frame(x$readiness_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$recommendation_table) && nrow(x$recommendation_table) > 0L) {
    cat("\nRecommendation table\n")
    print(round_numeric_df(as.data.frame(x$recommendation_table), digits = digits), row.names = FALSE)
  }

  invisible(x)
}

#' Plot a future arbitrary-facet planning active branch
#'
#' @param x Output from the future-branch active planning scaffold stored in
#'   `planning_schema$future_branch_active_branch`.
#' @param y Unused placeholder for generic compatibility.
#' @param type Plot type: `"profile_metrics"` for recommended deterministic
#'   profile values by metric, `"load_balance"` for recommended load/balance
#'   values by metric, `"coverage"` for recommended coverage/connectivity values
#'   by metric, `"readiness_tiers"` for counts of structural tiers across the
#'   current active-branch design grid, `"table_rows"` / `"role_tables"` /
#'   `"appendix_roles"` for summary-table bundle QC,
#'   `"appendix_sections"` / `"appendix_presets"` for manuscript-facing
#'   appendix selection counts, `"selection_handoff_presets"` for preset-level
#'   appendix handoff counts, `"selection_tables"` for appendix-selected
#'   future-branch tables ranked by row count within a preset,
#'   `"selection_handoff"` for section-aware plot-ready appendix handoff counts,
#'   `"selection_handoff_bundles"` for section-and-bundle plot-ready appendix
#'   handoff counts,
#'   `"selection_handoff_roles"` for role-aware plot-ready appendix handoff
#'   counts, `"selection_handoff_role_sections"` for role-by-section plot-ready
#'   appendix handoff counts,
#'   or `"selection_bundles"` / `"selection_roles"` / `"selection_sections"`
#'   for preset-filtered appendix selection summaries.
#' @param appendix_preset Appendix preset used for `selection_*` plot types.
#' @param selection_value For `selection_*` plot types, whether to plot exact
#'   counts (`"count"`) or the matching exact fraction (`"fraction"`) when that
#'   surface exposes one. `selection_tables` remains count-only because it
#'   represents table row counts rather than a normalized selection surface.
#' @param draw If `TRUE`, draw with base graphics; otherwise return plotting data.
#' @param main Optional title override.
#' @param palette Optional named color overrides.
#' @param label_angle Axis-label rotation angle.
#' @param ... Reserved for generic compatibility.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [summary.mfrm_future_branch_active_branch()]
#' @export
plot.mfrm_future_branch_active_branch <- function(x,
                                                  y = NULL,
                                                  type = c("profile_metrics", "load_balance", "coverage", "readiness_tiers", "table_rows", "role_tables", "appendix_roles", "appendix_sections", "appendix_presets", "selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections"),
                                                  appendix_preset = c("recommended", "compact", "all", "methods", "results", "diagnostics", "reporting"),
                                                  selection_value = c("count", "fraction"),
                                                  draw = TRUE,
                                                  main = NULL,
                                                  palette = NULL,
                                                  label_angle = 45,
                                                  ...) {
  type <- match.arg(type)
  appendix_preset <- match.arg(
    tolower(as.character(appendix_preset[1])),
    choices = c("recommended", "compact", "all", "methods", "results", "diagnostics", "reporting")
  )
  selection_value <- match.arg(selection_value)
  active <- simulation_future_branch_active_branch(x)

  if (type %in% c("table_rows", "role_tables", "appendix_roles", "appendix_sections", "appendix_presets")) {
    bundle <- build_summary_table_bundle(active)
    return(plot.mfrm_summary_table_bundle(
      bundle,
      type = type,
      selection_value = selection_value,
      main = main,
      palette = palette,
      label_angle = label_angle,
      draw = draw,
      ...
    ))
  }

  if (type %in% c("selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections")) {
    sx <- summary.mfrm_future_branch_active_branch(active)
    if (type == "selection_handoff_presets") {
      tbl <- as.data.frame(sx$selection_handoff_preset_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Preset", "PlotReadyTables") %in% names(tbl))) {
        stop("No future-branch appendix handoff-preset summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Preset)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by preset for `", appendix_preset, "`")
      default_palette <- c(selection_handoff_presets = "#f4a261", grid = "#ececec")
      plot_name <- "selection_handoff_presets"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_tables") {
      tbl <- future_branch_selection_table_summary(
        selection_catalog = sx$selection_catalog,
        preset = appendix_preset
      )
      if (nrow(tbl) == 0L || !all(c("Table", "Rows") %in% names(tbl))) {
        stop("No future-branch appendix table selection is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Table)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Selected appendix tables for preset `", appendix_preset, "`")
      default_palette <- c(selection_tables = "#e76f51", grid = "#ececec")
      plot_name <- "selection_tables"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff") {
      tbl <- as.data.frame(sx$selection_handoff_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "PlotReadyTables") %in% names(tbl))) {
        stop("No future-branch appendix handoff summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$AppendixSection)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by section for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff = "#ff9f1c", grid = "#ececec")
      plot_name <- "selection_handoff"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff_bundles") {
      tbl <- as.data.frame(sx$selection_handoff_bundle_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Bundle", "PlotReadyTables") %in% names(tbl))) {
        stop("No future-branch appendix handoff-bundle summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Bundle))
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by section and bundle for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff_bundles = "#5c677d", grid = "#ececec")
      plot_name <- "selection_handoff_bundles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff_roles") {
      tbl <- as.data.frame(sx$selection_handoff_role_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Role", "PlotReadyTables") %in% names(tbl))) {
        stop("No future-branch appendix handoff-role summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Role)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by role for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff_roles = "#9c6644", grid = "#ececec")
      plot_name <- "selection_handoff_roles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff_role_sections") {
      tbl <- as.data.frame(sx$selection_handoff_role_section_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Role", "PlotReadyTables") %in% names(tbl))) {
        stop("No future-branch appendix handoff role-section summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Role))
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by section and role for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff_role_sections = "#7f5539", grid = "#ececec")
      plot_name <- "selection_handoff_role_sections"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_bundles") {
      tbl <- as.data.frame(sx$selection_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Bundle", "TablesSelected") %in% names(tbl))) {
        stop("No future-branch appendix selection summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Bundle)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Appendix tables by bundle for preset `", appendix_preset, "`")
      default_palette <- c(selection_bundles = "#54a24b", grid = "#ececec")
      plot_name <- "selection_bundles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_roles") {
      tbl <- as.data.frame(sx$selection_role_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Role", "Tables") %in% names(tbl))) {
        stop("No future-branch appendix role summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Role)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Selected appendix roles for preset `", appendix_preset, "`")
      default_palette <- c(selection_roles = "#b279a2", grid = "#ececec")
      plot_name <- "selection_roles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else {
      tbl <- as.data.frame(sx$selection_section_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Tables") %in% names(tbl))) {
        stop("No future-branch appendix section summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$AppendixSection)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Selected appendix sections for preset `", appendix_preset, "`")
      default_palette <- c(selection_sections = "#2a9d8f", grid = "#ececec")
      plot_name <- "selection_sections"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    }

    keep <- is.finite(values) & nzchar(labels)
    if (!any(keep)) {
      stop("Selected future-branch appendix plot has no finite values to display.", call. = FALSE)
    }
    labels <- labels[keep]
    values <- values[keep]
    tbl <- tbl[keep, , drop = FALSE]
    pal <- resolve_palette(palette = palette, defaults = default_palette)
    plot_title <- if (is.null(main)) paste("Future branch", gsub("_", " ", plot_name)) else as.character(main[1])

    if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = pal[plot_name],
        main = plot_title,
        ylab = ylab,
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }

    return(invisible(new_mfrm_plot_data(
      "future_branch_active_branch",
      list(
        plot = plot_name,
        selection_value = measure$selection_value,
        appendix_preset = appendix_preset,
        table = tbl,
        title = plot_title,
        subtitle = subtitle,
        legend = new_plot_legend(legend_label, "future_branch", "bar", pal[plot_name]),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (type == "profile_metrics") {
    tbl <- simulation_future_branch_active_branch_profile(active)$profile_summary_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !all(c("metric", "recommended_value") %in% names(tbl))) {
      stop("No future-branch profile metrics are available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$metric)
    values <- suppressWarnings(as.numeric(tbl$recommended_value))
    subtitle <- "Recommended deterministic profile values"
    default_palette <- c(profile_metrics = "#3a7ca5", grid = "#ececec")
    plot_name <- "profile_metrics"
    legend_label <- "Recommended value"
  } else if (type == "load_balance") {
    tbl <- simulation_future_branch_active_branch_load_balance(active)$diagnostic_summary_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !all(c("metric", "recommended_value") %in% names(tbl))) {
      stop("No future-branch load/balance metrics are available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$metric)
    values <- suppressWarnings(as.numeric(tbl$recommended_value))
    subtitle <- "Recommended deterministic load/balance values"
    default_palette <- c(load_balance = "#2a9d8f", grid = "#ececec")
    plot_name <- "load_balance"
    legend_label <- "Recommended value"
  } else if (type == "coverage") {
    tbl <- simulation_future_branch_active_branch_coverage(active)$diagnostic_summary_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !all(c("metric", "recommended_value") %in% names(tbl))) {
      stop("No future-branch coverage metrics are available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$metric)
    values <- suppressWarnings(as.numeric(tbl$recommended_value))
    subtitle <- "Recommended deterministic coverage/connectivity values"
    default_palette <- c(coverage = "#f4a261", grid = "#ececec")
    plot_name <- "coverage"
    legend_label <- "Recommended value"
  } else {
    tbl <- simulation_future_branch_active_branch_readiness(active)$readiness_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !"structural_tier" %in% names(tbl)) {
      stop("No future-branch readiness tiers are available for plotting.", call. = FALSE)
    }
    counts <- sort(table(as.character(tbl$structural_tier)), decreasing = TRUE)
    labels <- names(counts)
    values <- as.numeric(counts)
    subtitle <- "Structural tier counts across active-branch designs"
    default_palette <- c(readiness_tiers = "#b279a2", grid = "#ececec")
    plot_name <- "readiness_tiers"
    legend_label <- "Designs"
  }

  keep <- is.finite(values) & nzchar(labels)
  if (!any(keep)) {
    stop("Selected future-branch plot has no finite values to display.", call. = FALSE)
  }
  labels <- labels[keep]
  values <- values[keep]
  pal <- resolve_palette(palette = palette, defaults = default_palette)
  color_name <- names(default_palette)[1]
  plot_title <- if (is.null(main)) "Future arbitrary-facet planning profile" else as.character(main[1])

  if (isTRUE(draw)) {
    barplot_rot45(
      height = values,
      labels = labels,
      col = pal[color_name],
      main = plot_title,
      ylab = legend_label,
      label_angle = label_angle,
      mar_bottom = 9
    )
    graphics::abline(h = 0, col = pal["grid"], lty = 2)
  }

  invisible(new_mfrm_plot_data(
    "future_branch_active_branch",
    list(
      plot = plot_name,
      label = labels,
      value = values,
      title = plot_title,
      subtitle = subtitle,
      recommended_design_id = as.character(active$recommended_design_id %||% NA_character_),
      legend = new_plot_legend(legend_label, "future_branch", "bar", pal[color_name]),
      reference_lines = new_reference_lines("h", 0, "Zero reference", "dashed", "reference")
    )
  ))
}

simulation_compact_future_branch_report_cache <- function(x) {
  if (!is.list(x) || is.data.frame(x) || inherits(x, "mfrm_plot_data")) {
    return(x)
  }

  nested_report_fields <- c(
    "report_bundle",
    "report_summary",
    "report_overview",
    "report_catalog",
    "report_surface_registry",
    "report_panel",
    "source_object",
    "registry",
    "selected_surface"
  )
  x[intersect(names(x), nested_report_fields)] <- NULL

  for (nm in names(x)) {
    if (is.list(x[[nm]]) && !is.data.frame(x[[nm]]) && !inherits(x[[nm]], "mfrm_plot_data")) {
      x[[nm]] <- simulation_compact_future_branch_report_cache(x[[nm]])
    }
  }

  x
}

simulation_compact_future_branch_schema_report_cache <- function(branch_schema) {
  report_fields <- c(
    "report_summary",
    "report_overview_table",
    "report_catalog",
    "report_digest",
    "report_surface_registry",
    "report_panel",
    "report_operation",
    "report_snapshot",
    "report_brief",
    "report_consumer"
  )
  for (nm in intersect(report_fields, names(branch_schema))) {
    branch_schema[[nm]] <- simulation_compact_future_branch_report_cache(branch_schema[[nm]])
  }
  branch_schema
}

simulation_compact_future_branch_schema_active_cache <- function(branch_schema) {
  if (is.list(branch_schema$active_branch)) {
    if (is.list(branch_schema$active_branch$summary)) {
      branch_schema$active_branch$summary$pilot <- NULL
    }
    if (is.list(branch_schema$active_branch$table)) {
      branch_schema$active_branch$table$pilot <- NULL
    }
    if (is.list(branch_schema$active_branch$plot)) {
      branch_schema$active_branch$plot$pilot <- NULL
    }
  }

  drop_map <- list(
    pilot_summary = "pilot",
    pilot_table = "pilot",
    pilot_plot = "pilot",
    active_branch_profile = "active_branch",
    active_branch_load_balance = "active_branch_profile",
    active_branch_overview = c("active_branch", "active_branch_profile"),
    active_branch_load_balance_overview = "active_branch_load_balance",
    active_branch_coverage = "active_branch_profile",
    active_branch_coverage_overview = "active_branch_coverage",
    active_branch_guardrails = c("active_branch_coverage", "active_branch_load_balance"),
    active_branch_guardrail_overview = "active_branch_guardrails",
    active_branch_readiness = "active_branch_guardrails",
    active_branch_readiness_overview = "active_branch_readiness",
    active_branch_recommendation = c("active_branch_readiness", "active_branch_profile"),
    active_branch_recommendation_overview = "active_branch_recommendation"
  )

  for (nm in intersect(names(drop_map), names(branch_schema))) {
    branch_schema[[nm]][intersect(drop_map[[nm]], names(branch_schema[[nm]]))] <- NULL
  }

  branch_schema
}

simulation_future_branch_schema <- function(sim_spec = NULL, facet_names = NULL) {
  future_facet_table <- if (is.null(facet_names)) {
    simulation_future_facet_table(sim_spec)
  } else {
    simulation_future_facet_table(facet_names = facet_names)
  }
  future_design_template <- if (is.null(facet_names)) {
    simulation_future_design_template(sim_spec)
  } else {
    simulation_future_design_template(facet_names = facet_names)
  }
  design_schema <- if (is.null(facet_names)) {
    simulation_future_branch_design_schema(sim_spec)
  } else {
    simulation_future_branch_design_schema(facet_names = facet_names)
  }
  branch_schema <- list(
    planner_contract = "arbitrary_facet_planning_scaffold",
    planner_stage = "schema_only",
    input_contract = "design$facets(named counts)",
    facet_table = future_facet_table,
    design_template = future_design_template,
    assignment_axis = design_schema$assignment_axis,
    design_schema = design_schema,
    grid_semantics = design_schema$grid_semantics,
    note = paste(
      "Schema-only future-branch contract bundling the stable facet-count table",
      "and matching `design$facets(named counts)` template for a later",
      "arbitrary-facet planner, together with a preview-ready nested",
      "design schema."
    )
  )
  branch_schema$grid_contract <- simulation_future_branch_grid_contract(branch_schema)
  branch_schema$preview <- simulation_future_branch_preview(branch_schema)
  branch_schema$grid_bundle <- simulation_future_branch_grid_bundle(branch_schema)
  branch_schema$grid_context <- simulation_future_branch_grid_context(branch_schema)
  branch_schema$report_bundle <- simulation_future_branch_report_bundle(branch_schema)
  branch_schema$report_summary <- simulation_future_branch_report_summary(branch_schema)
  branch_schema$report_overview_table <- simulation_future_branch_report_overview_table(branch_schema)
  branch_schema$report_catalog <- simulation_future_branch_report_catalog(branch_schema)
  branch_schema$report_digest <- simulation_future_branch_report_digest(branch_schema)
  branch_schema$report_surface_registry <- simulation_future_branch_report_surface_registry(branch_schema)
  branch_schema$report_panel <- simulation_future_branch_report_panel(branch_schema)
  branch_schema$report_operation <- simulation_future_branch_report_operation(branch_schema)
  branch_schema$report_snapshot <- simulation_future_branch_report_snapshot(branch_schema)
  branch_schema$report_brief <- simulation_future_branch_report_brief(branch_schema)
  branch_schema$report_mode_registry <- simulation_future_branch_report_mode_registry(branch_schema)
  branch_schema$report_consumer <- simulation_future_branch_report_consume(branch_schema)
  branch_schema$pilot <- simulation_future_branch_pilot(branch_schema)
  branch_schema$pilot_summary <- simulation_future_branch_pilot_summary(branch_schema)
  branch_schema$pilot_table <- simulation_future_branch_pilot_table(branch_schema)
  branch_schema$pilot_plot <- simulation_future_branch_pilot_plot(branch_schema)
  branch_schema$active_branch <- simulation_future_branch_active_branch(branch_schema)
  branch_schema$active_branch_profile <- simulation_future_branch_active_branch_profile(branch_schema)
  branch_schema$active_branch_load_balance <- simulation_future_branch_active_branch_load_balance(branch_schema)
  branch_schema$active_branch_overview <- simulation_future_branch_active_branch_overview(branch_schema)
  branch_schema$active_branch_load_balance_overview <- simulation_future_branch_active_branch_load_balance_overview(branch_schema)
  branch_schema$active_branch_coverage <- simulation_future_branch_active_branch_coverage(branch_schema)
  branch_schema$active_branch_coverage_overview <- simulation_future_branch_active_branch_coverage_overview(branch_schema)
  branch_schema$active_branch_guardrails <- simulation_future_branch_active_branch_guardrails(branch_schema)
  branch_schema$active_branch_guardrail_overview <- simulation_future_branch_active_branch_guardrail_overview(branch_schema)
  branch_schema$active_branch_readiness <- simulation_future_branch_active_branch_readiness(branch_schema)
  branch_schema$active_branch_readiness_overview <- simulation_future_branch_active_branch_readiness_overview(branch_schema)
  branch_schema$active_branch_recommendation <- simulation_future_branch_active_branch_recommendation(branch_schema)
  branch_schema$active_branch_recommendation_overview <- simulation_future_branch_active_branch_recommendation_overview(branch_schema)
  branch_schema <- simulation_compact_future_branch_schema_report_cache(branch_schema)
  branch_schema <- simulation_compact_future_branch_schema_active_cache(branch_schema)
  branch_schema
}

simulation_planning_scope <- function(sim_spec = NULL, facet_names = NULL) {
  if (is.null(facet_names)) {
    facet_names <- simulation_spec_output_facet_names(sim_spec)
    aliases <- simulation_design_variable_aliases(sim_spec)
    facet_manifest <- simulation_facet_manifest(sim_spec)
  } else {
    facet_names <- simulation_validate_output_facet_names(facet_names)
    aliases <- simulation_design_variable_aliases(list(facet_names = facet_names))
    facet_manifest <- simulation_facet_manifest(facet_names = facet_names)
  }
  list(
    planner_contract = "role_based_two_non_person_facets",
    planner_stage = "first_release",
    supports_arbitrary_facet_planning = FALSE,
    supports_arbitrary_facet_estimation = TRUE,
    supported_non_person_roles = c("rater", "criterion"),
    supported_non_person_facet_count = 2L,
    role_labels = unname(facet_names),
    design_variables = c("n_person", "n_rater", "n_criterion", "raters_per_person"),
    design_variable_aliases = aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")],
    facet_manifest = facet_manifest,
    future_planner_contract = "arbitrary_facet_planning_scaffold",
    future_planner_stage = "schema_only",
    future_branch_input_contract = "design$facets(named counts)",
    note = paste0(
      "Current planning helpers vary one person count and exactly two non-person facet roles (",
      facet_names[1], " and ", facet_names[2], "). ",
      "The estimation core supports arbitrary facet counts, but planning/forecasting remain role-based until a fully arbitrary-facet planner is validated. ",
      "A facet manifest is now exposed so a future arbitrary-facet branch can reuse the same public facet labels without changing the current planner contract."
    )
  )
}

simulation_planning_scope_note <- function(scope) {
  if (is.list(scope) && is.character(scope$note) && length(scope$note) > 0L) {
    note <- as.character(scope$note[1])
    if (nzchar(note)) {
      return(note)
    }
  }
  character(0)
}

simulation_planning_constraints <- function(sim_spec = NULL) {
  all_vars <- c("n_person", "n_rater", "n_criterion", "raters_per_person")
  if (is.null(sim_spec) || !inherits(sim_spec, "mfrm_sim_spec")) {
    return(list(
      mutable_design_variables = all_vars,
      locked_design_variables = character(0),
      lock_reasons = stats::setNames(character(0), character(0)),
      feasibility_rule = "`raters_per_person <= n_rater`",
      note = "Current scalar-argument planning paths allow `n_person`, `n_rater`, `n_criterion`, and `raters_per_person` to vary subject to `raters_per_person <= n_rater`."
    ))
  }

  mutable <- all_vars
  reasons <- character(0)
  assignment <- as.character(sim_spec$assignment %||% "rotating")

  if (identical(assignment, "resampled")) {
    reasons["n_rater"] <- paste0(
      "`assignment = \"resampled\"` reuses empirical person-level ",
      simulation_spec_output_facet_names(sim_spec)[1],
      " profiles."
    )
    reasons["raters_per_person"] <- reasons[["n_rater"]]
  }
  if (identical(assignment, "skeleton")) {
    skeleton_reason <- "`assignment = \"skeleton\"` reuses the observed person-by-facet response skeleton."
    reasons["n_rater"] <- skeleton_reason
    reasons["n_criterion"] <- skeleton_reason
    reasons["raters_per_person"] <- skeleton_reason
  }

  threshold_mode <- simulation_spec_threshold_mode(sim_spec)
  if (identical(threshold_mode, "step_facet_specific")) {
    role <- simulation_step_facet_role(sim_spec, step_facet = sim_spec$step_facet)
    lock_var <- if (identical(role, "criterion")) "n_criterion" else if (identical(role, "rater")) "n_rater" else NULL
    if (!is.null(lock_var)) {
      reasons[lock_var] <- paste0(
        "`sim_spec` contains step-facet-specific thresholds for `",
        sim_spec$step_facet,
        "`."
      )
    }
  }

  slope_mode <- simulation_spec_slope_mode(sim_spec)
  if (identical(slope_mode, "slope_facet_specific")) {
    role <- simulation_step_facet_role(sim_spec, step_facet = sim_spec$slope_facet)
    lock_var <- if (identical(role, "criterion")) "n_criterion" else if (identical(role, "rater")) "n_rater" else NULL
    if (!is.null(lock_var)) {
      reasons[lock_var] <- paste0(
        "First-release `GPCM` stores slope values for `",
        sim_spec$slope_facet,
        "` levels."
      )
    }
  }

  reason_names <- names(reasons)
  if (is.null(reason_names)) {
    reason_names <- character(0)
  }
  locked <- intersect(all_vars, reason_names)
  mutable <- setdiff(all_vars, locked)
  note <- if (length(locked) == 0L) {
    "All current design variables remain mutable subject to `raters_per_person <= n_rater`."
  } else {
    paste0(
      "Current planning path allows changing ",
      paste(mutable, collapse = ", "),
      "; locked variables are ",
      paste(locked, collapse = ", "),
      "."
    )
  }

  list(
    mutable_design_variables = mutable,
    locked_design_variables = locked,
    lock_reasons = reasons[locked],
    feasibility_rule = "`raters_per_person <= n_rater`",
    note = note
  )
}

simulation_planning_constraints_note <- function(constraints) {
  if (is.list(constraints) && is.character(constraints$note) && length(constraints$note) > 0L) {
    note <- as.character(constraints$note[1])
    if (nzchar(note)) {
      return(note)
    }
  }
  character(0)
}

simulation_planning_schema <- function(sim_spec = NULL, facet_names = NULL) {
  descriptor <- if (is.null(facet_names)) {
    simulation_design_descriptor(sim_spec)
  } else {
    simulation_design_descriptor(list(facet_names = simulation_validate_output_facet_names(facet_names)))
  }
  scope <- if (is.null(facet_names)) {
    simulation_planning_scope(sim_spec)
  } else {
    simulation_planning_scope(facet_names = facet_names)
  }
  constraints <- simulation_planning_constraints(sim_spec)
  lock_lookup <- constraints$lock_reasons %||% stats::setNames(character(0), character(0))

  role_table <- descriptor |>
    dplyr::mutate(
      axis_class = c("person_count", "facet_level_count", "facet_level_count", "assignment_count"),
      depends_on_role = c(NA_character_, NA_character_, NA_character_, "rater"),
      mutable = .data$canonical %in% constraints$mutable_design_variables,
      locked = .data$canonical %in% constraints$locked_design_variables,
      lock_reason = unname(lock_lookup[.data$canonical])
    )
  role_table$lock_reason[is.na(role_table$lock_reason)] <- ""
  facet_manifest <- if (is.null(facet_names)) {
    simulation_facet_manifest(sim_spec)
  } else {
    simulation_facet_manifest(facet_names = facet_names)
  }
  future_facet_table <- if (is.null(facet_names)) {
    simulation_future_facet_table(sim_spec)
  } else {
    simulation_future_facet_table(facet_names = facet_names)
  }
  future_design_template <- if (is.null(facet_names)) {
    simulation_future_design_template(sim_spec)
  } else {
    simulation_future_design_template(facet_names = facet_names)
  }
  future_branch_schema <- if (is.null(facet_names)) {
    simulation_future_branch_schema(sim_spec)
  } else {
    simulation_future_branch_schema(facet_names = facet_names)
  }

  list(
    planner_contract = scope$planner_contract,
    planner_stage = scope$planner_stage,
    supports_arbitrary_facet_planning = scope$supports_arbitrary_facet_planning,
    supports_arbitrary_facet_estimation = scope$supports_arbitrary_facet_estimation,
    role_labels = scope$role_labels,
    design_variables = scope$design_variables,
    design_variable_aliases = scope$design_variable_aliases,
    facet_manifest = facet_manifest,
    future_planner_contract = scope$future_planner_contract,
    future_planner_stage = scope$future_planner_stage,
    future_branch_input_contract = scope$future_branch_input_contract,
    future_facet_table = future_facet_table,
    future_design_template = future_design_template,
    future_branch_schema = future_branch_schema,
    future_branch_preview = future_branch_schema$preview,
    future_branch_grid_semantics = future_branch_schema$grid_semantics,
    future_branch_grid_contract = future_branch_schema$grid_contract,
    future_branch_grid_bundle = future_branch_schema$grid_bundle,
    future_branch_grid_context = future_branch_schema$grid_context,
    future_branch_report_bundle = future_branch_schema$report_bundle,
    future_branch_report_summary = future_branch_schema$report_summary,
    future_branch_report_overview_table = future_branch_schema$report_overview_table,
    future_branch_report_catalog = future_branch_schema$report_catalog,
    future_branch_report_digest = future_branch_schema$report_digest,
    future_branch_report_surface_registry = future_branch_schema$report_surface_registry,
    future_branch_report_panel = future_branch_schema$report_panel,
    future_branch_report_operation = future_branch_schema$report_operation,
    future_branch_report_snapshot = future_branch_schema$report_snapshot,
    future_branch_report_brief = future_branch_schema$report_brief,
    future_branch_report_mode_registry = future_branch_schema$report_mode_registry,
    future_branch_report_consumer = future_branch_schema$report_consumer,
    future_branch_pilot = future_branch_schema$pilot,
    future_branch_pilot_summary = future_branch_schema$pilot_summary,
    future_branch_pilot_table = future_branch_schema$pilot_table,
    future_branch_pilot_plot = future_branch_schema$pilot_plot,
    future_branch_active_branch = future_branch_schema$active_branch,
    future_branch_active_branch_profile = future_branch_schema$active_branch_profile,
    future_branch_active_branch_load_balance = future_branch_schema$active_branch_load_balance,
    future_branch_active_branch_overview = future_branch_schema$active_branch_overview,
    future_branch_active_branch_load_balance_overview = future_branch_schema$active_branch_load_balance_overview,
    future_branch_active_branch_coverage = future_branch_schema$active_branch_coverage,
    future_branch_active_branch_coverage_overview = future_branch_schema$active_branch_coverage_overview,
    future_branch_active_branch_guardrails = future_branch_schema$active_branch_guardrails,
    future_branch_active_branch_guardrail_overview = future_branch_schema$active_branch_guardrail_overview,
    future_branch_active_branch_readiness = future_branch_schema$active_branch_readiness,
    future_branch_active_branch_readiness_overview = future_branch_schema$active_branch_readiness_overview,
    future_branch_active_branch_recommendation = future_branch_schema$active_branch_recommendation,
    future_branch_active_branch_recommendation_overview = future_branch_schema$active_branch_recommendation_overview,
    role_table = role_table,
    feasibility_rules = constraints$feasibility_rule,
    note = paste(
      "Current planning schema exposes one person-count axis, two non-person facet-count axes,",
      "and one assignments-per-person axis under the first-release role-based planner,",
      "while exposing a facet manifest and a nested schema-only future-branch",
      "contract for arbitrary-facet planning, including machine-readable",
      "future-branch preview, grid, bundle, context, report-bundle,",
      "report-summary, report-overview, report-catalog, and report-digest metadata",
      "plus a report-surface registry and compact report panel when default",
      "counts are available, alongside one combined report operation object",
      "plus one lightweight report snapshot, one selected-surface report brief,",
      "one mode-based report consumer, one internal active pilot object,",
      "compact pilot-level summary/table/plot consumers, and one bundled",
      "active-branch object plus deterministic active-branch profile,",
      "load/balance diagnostics, coverage/connectivity diagnostics,",
      "guardrail classifications, structural readiness summaries,",
      "and conservative recommendation/overview contracts."
    )
  )
}

simulation_planning_schema_note <- function(schema) {
  if (is.list(schema) && is.character(schema$note) && length(schema$note) > 0L) {
    note <- as.character(schema$note[1])
    if (nzchar(note)) {
      return(note)
    }
  }
  character(0)
}

simulation_object_design_variable_aliases <- function(x) {
  aliases <- x$design_variable_aliases %||% x$settings$design_variable_aliases %||% NULL
  if (is.character(aliases) &&
      identical(sort(names(aliases)), sort(c("n_person", "n_rater", "n_criterion", "raters_per_person")))) {
    return(aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")])
  }
  sim_spec <- x$settings$sim_spec %||% NULL
  simulation_design_variable_aliases(sim_spec)
}

simulation_object_design_descriptor <- function(x) {
  descriptor <- x$design_descriptor %||% x$settings$design_descriptor %||% x$ademp$data_generating_mechanism$design_descriptor %||% NULL
  required_cols <- c("role", "canonical", "alias", "facet", "quantity", "description")
  if (is.data.frame(descriptor) && all(required_cols %in% names(descriptor))) {
    descriptor <- tibble::as_tibble(descriptor)
    return(descriptor[, required_cols])
  }
  sim_spec <- x$settings$sim_spec %||% NULL
  simulation_design_descriptor(sim_spec)
}

simulation_object_planning_scope <- function(x) {
  scope <- x$planning_scope %||% x$settings$planning_scope %||%
    x$ademp$data_generating_mechanism$planning_scope %||%
    x$sim_spec$planning_scope %||% NULL
  required_fields <- c(
    "planner_contract",
    "planner_stage",
    "supports_arbitrary_facet_planning",
    "supports_arbitrary_facet_estimation",
    "supported_non_person_roles",
    "supported_non_person_facet_count",
    "role_labels",
    "design_variables",
    "design_variable_aliases",
    "facet_manifest",
    "future_planner_contract",
    "future_planner_stage",
    "future_branch_input_contract",
    "note"
  )
  if (is.list(scope) &&
      all(required_fields %in% names(scope)) &&
      is.data.frame(scope$facet_manifest)) {
    return(scope[required_fields])
  }
  sim_spec <- x$settings$sim_spec %||% x$sim_spec %||% NULL
  simulation_planning_scope(sim_spec)
}

simulation_object_planning_constraints <- function(x) {
  constraints <- x$planning_constraints %||% x$settings$planning_constraints %||%
    x$ademp$data_generating_mechanism$planning_constraints %||%
    x$sim_spec$planning_constraints %||% NULL
  required_fields <- c(
    "mutable_design_variables",
    "locked_design_variables",
    "lock_reasons",
    "feasibility_rule",
    "note"
  )
  if (is.list(constraints) && all(required_fields %in% names(constraints))) {
    return(constraints[required_fields])
  }
  sim_spec <- x$settings$sim_spec %||% x$sim_spec %||% NULL
  simulation_planning_constraints(sim_spec)
}

simulation_object_planning_schema <- function(x) {
  schema <- x$planning_schema %||% x$settings$planning_schema %||%
    x$ademp$data_generating_mechanism$planning_schema %||%
    x$sim_spec$planning_schema %||% NULL
  required_fields <- c(
    "planner_contract",
    "planner_stage",
    "supports_arbitrary_facet_planning",
    "supports_arbitrary_facet_estimation",
    "role_labels",
    "design_variables",
    "design_variable_aliases",
    "facet_manifest",
    "future_planner_contract",
    "future_planner_stage",
    "future_branch_input_contract",
    "future_facet_table",
    "future_design_template",
    "future_branch_schema",
    "future_branch_preview",
    "future_branch_grid_semantics",
    "future_branch_grid_contract",
    "future_branch_grid_bundle",
    "future_branch_grid_context",
    "future_branch_report_bundle",
    "future_branch_report_summary",
    "future_branch_report_overview_table",
    "future_branch_report_catalog",
    "future_branch_report_digest",
    "future_branch_report_surface_registry",
    "future_branch_report_panel",
    "future_branch_report_operation",
    "future_branch_report_snapshot",
    "future_branch_report_brief",
    "future_branch_report_mode_registry",
    "future_branch_report_consumer",
    "future_branch_pilot",
    "future_branch_pilot_summary",
    "future_branch_pilot_table",
    "future_branch_pilot_plot",
    "future_branch_active_branch",
    "future_branch_active_branch_profile",
    "future_branch_active_branch_load_balance",
    "future_branch_active_branch_overview",
    "future_branch_active_branch_load_balance_overview",
    "future_branch_active_branch_coverage",
    "future_branch_active_branch_coverage_overview",
    "future_branch_active_branch_guardrails",
    "future_branch_active_branch_guardrail_overview",
    "future_branch_active_branch_readiness",
    "future_branch_active_branch_readiness_overview",
    "future_branch_active_branch_recommendation",
    "future_branch_active_branch_recommendation_overview",
    "role_table",
    "feasibility_rules",
    "note"
  )
  if (is.list(schema) &&
      all(required_fields %in% names(schema)) &&
      is.data.frame(schema$role_table) &&
      is.data.frame(schema$facet_manifest) &&
      is.data.frame(schema$future_facet_table) &&
      is.list(schema$future_design_template) &&
      is.list(schema$future_branch_schema) &&
      is.list(schema$future_branch_preview) &&
      is.list(schema$future_branch_grid_semantics) &&
      is.list(schema$future_branch_grid_contract) &&
      is.list(schema$future_branch_grid_bundle) &&
      is.list(schema$future_branch_grid_context) &&
      is.list(schema$future_branch_report_bundle) &&
      is.list(schema$future_branch_report_summary) &&
      is.list(schema$future_branch_report_overview_table) &&
      is.list(schema$future_branch_report_catalog) &&
      is.list(schema$future_branch_report_digest) &&
      is.list(schema$future_branch_report_surface_registry) &&
      is.list(schema$future_branch_report_panel) &&
      is.list(schema$future_branch_report_operation) &&
      is.list(schema$future_branch_report_snapshot) &&
      is.list(schema$future_branch_report_brief) &&
      is.data.frame(schema$future_branch_report_mode_registry) &&
      is.list(schema$future_branch_report_consumer) &&
      is.list(schema$future_branch_pilot) &&
      is.list(schema$future_branch_pilot_summary) &&
      is.list(schema$future_branch_pilot_table) &&
      is.list(schema$future_branch_pilot_plot) &&
      is.list(schema$future_branch_active_branch) &&
      is.list(schema$future_branch_active_branch_profile) &&
      is.list(schema$future_branch_active_branch_load_balance) &&
      is.list(schema$future_branch_active_branch_overview) &&
      is.list(schema$future_branch_active_branch_load_balance_overview) &&
      is.list(schema$future_branch_active_branch_coverage) &&
      is.list(schema$future_branch_active_branch_coverage_overview) &&
      is.list(schema$future_branch_active_branch_guardrails) &&
      is.list(schema$future_branch_active_branch_guardrail_overview) &&
      is.list(schema$future_branch_active_branch_readiness) &&
      is.list(schema$future_branch_active_branch_readiness_overview) &&
      is.list(schema$future_branch_active_branch_recommendation) &&
      is.list(schema$future_branch_active_branch_recommendation_overview)) {
    return(schema[required_fields])
  }
  sim_spec <- x$settings$sim_spec %||% x$sim_spec %||% NULL
  simulation_planning_schema(sim_spec)
}

simulation_design_variable_choices <- function(aliases, descriptor = NULL) {
  role_names <- if (is.data.frame(descriptor) && "role" %in% names(descriptor)) as.character(descriptor$role) else character(0)
  unique(c(names(aliases), unname(aliases), role_names))
}

simulation_design_variable_lookup <- function(aliases, descriptor = NULL) {
  canonical <- names(aliases)
  lookup <- c(stats::setNames(canonical, canonical), stats::setNames(canonical, unname(aliases)))
  if (is.data.frame(descriptor) && all(c("role", "canonical") %in% names(descriptor))) {
    role_lookup <- stats::setNames(as.character(descriptor$canonical), as.character(descriptor$role))
    lookup <- c(lookup, role_lookup)
  }
  lookup[!duplicated(names(lookup))]
}

simulation_resolve_design_variable <- function(value, aliases, arg_name, descriptor = NULL) {
  value <- as.character(value[1])
  lookup <- simulation_design_variable_lookup(aliases, descriptor = descriptor)
  resolved <- unname(lookup[[value]])
  if (!is.null(resolved) && nzchar(resolved)) {
    return(resolved)
  }
  stop(
    "`", arg_name, "` must be one of: ",
    paste(simulation_design_variable_choices(aliases, descriptor = descriptor), collapse = ", "),
    ".",
    call. = FALSE
  )
}

simulation_resolve_design_variable_vector <- function(values, aliases, descriptor = NULL) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    return(character(0))
  }
  lookup <- simulation_design_variable_lookup(aliases, descriptor = descriptor)
  resolved <- unname(lookup[values])
  unique(resolved[!is.na(resolved) & nzchar(resolved)])
}

simulation_design_variable_label <- function(value, aliases) {
  label <- unname(aliases[[value]])
  if (is.null(label) || is.na(label) || !nzchar(label)) {
    return(value)
  }
  label
}

simulation_parse_future_facet_input <- function(facets,
                                                future_facet_table,
                                                arg_name = "design$facets") {
  if (is.null(facets)) {
    return(list())
  }
  if (is.data.frame(facets)) {
    if (nrow(facets) != 1L) {
      stop("`", arg_name, "` must be a one-row data frame when supplied as a data frame.",
           call. = FALSE)
    }
    facets <- as.list(facets[1, , drop = FALSE])
  } else if (is.atomic(facets) && !is.list(facets)) {
    facets <- as.list(facets)
  } else if (!is.list(facets)) {
    stop(
      "`", arg_name, "` must be a named list, named vector, or one-row data frame keyed by future facet names.",
      call. = FALSE
    )
  }

  raw_names <- names(facets)
  valid_names <- as.character(future_facet_table$future_facet_key)
  if (is.null(raw_names) || any(!nzchar(raw_names))) {
    stop(
      "`", arg_name, "` must use names such as ",
      paste(valid_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  lookup <- stats::setNames(
    as.character(future_facet_table$current_planning_count_variable),
    valid_names
  )
  resolved <- unname(lookup[raw_names])
  if (any(is.na(resolved) | !nzchar(resolved))) {
    bad <- unique(raw_names[is.na(resolved) | !nzchar(resolved)])
    stop(
      "`", arg_name, "` must use names such as ",
      paste(valid_names, collapse = ", "),
      ". Invalid names: ",
      paste(bad, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(resolved)) {
    dup <- unique(resolved[duplicated(resolved)])
    stop(
      "`", arg_name, "` supplies the same facet-count variable more than once: ",
      paste(dup, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  stats::setNames(unname(facets), resolved)
}

simulation_parse_future_branch_design <- function(design,
                                                  future_branch_schema,
                                                  arg_name = "design") {
  if (is.null(future_branch_schema) || !is.list(future_branch_schema)) {
    return(list())
  }
  design_schema <- simulation_coerce_future_branch_design_schema(future_branch_schema)

  parsed <- list()
  if ("facets" %in% names(design)) {
    parsed <- c(
      parsed,
      simulation_parse_future_facet_input(
        facets = design[["facets"]],
        future_facet_table = dplyr::transmute(
          design_schema$facet_axes,
          future_facet_key = .data$input_key,
          current_planning_count_variable = .data$canonical_design_variable
        ),
        arg_name = paste0(arg_name, "$facets")
      )
    )
  }

  assignment_axis <- design_schema$assignment_axis
  assignment_key <- as.character(assignment_axis$future_input_key %||% "assignment")
  assignment_var <- as.character(
    assignment_axis$current_planning_count_variable %||% "raters_per_person"
  )
  if (assignment_key %in% names(design)) {
    parsed[[assignment_var]] <- design[[assignment_key]]
  }

  parsed
}

simulation_parse_design_input <- function(design,
                                          aliases,
                                          descriptor = NULL,
                                          future_facet_table = NULL,
                                          future_branch_schema = NULL,
                                          arg_name = "design") {
  if (is.null(design)) {
    return(list())
  }

  if (is.data.frame(design)) {
    if (nrow(design) != 1L) {
      stop("`", arg_name, "` must be a one-row data frame when supplied as a data frame.",
           call. = FALSE)
    }
    design <- as.list(design[1, , drop = FALSE])
  } else if (is.atomic(design) && !is.list(design)) {
    design <- as.list(design)
  } else if (!is.list(design)) {
    stop(
      "`", arg_name, "` must be a named list, named vector, or one-row data frame. Valid names: ",
      paste(simulation_design_variable_choices(aliases, descriptor = descriptor), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  parsed_future_branch <- list()
  if ("facets" %in% names(design)) {
    if (!is.null(future_branch_schema) && is.list(future_branch_schema)) {
      parsed_future_branch <- simulation_parse_future_branch_design(
        design = design,
        future_branch_schema = future_branch_schema,
        arg_name = arg_name
      )
      assignment_key <- as.character(
        future_branch_schema$assignment_axis$future_input_key %||% "assignment"
      )
      design[["facets"]] <- NULL
      if (assignment_key %in% names(design)) {
        design[[assignment_key]] <- NULL
      }
    } else {
      if (is.null(future_facet_table) || !is.data.frame(future_facet_table)) {
        stop("`", arg_name, "$facets` is not available for this planning object.", call. = FALSE)
      }
      parsed_future_branch <- simulation_parse_future_facet_input(
        facets = design[["facets"]],
        future_facet_table = future_facet_table,
        arg_name = paste0(arg_name, "$facets")
      )
      design[["facets"]] <- NULL
    }
  } else if (!is.null(future_branch_schema) && is.list(future_branch_schema)) {
    assignment_key <- as.character(
      future_branch_schema$assignment_axis$future_input_key %||% "assignment"
    )
    if (assignment_key %in% names(design)) {
      parsed_future_branch <- simulation_parse_future_branch_design(
        design = design,
        future_branch_schema = future_branch_schema,
        arg_name = arg_name
      )
      design[[assignment_key]] <- NULL
    }
  }

  raw_names <- names(design)
  if (length(design) == 0L) {
    raw_names <- character(0)
  } else if (is.null(raw_names) || any(!nzchar(raw_names))) {
    stop(
      "`", arg_name, "` must use names such as ",
      paste(simulation_design_variable_choices(aliases, descriptor = descriptor), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  parsed_top_level <- list()
  if (length(raw_names) > 0L) {
    canonical_names <- vapply(
      raw_names,
      simulation_resolve_design_variable,
      aliases = aliases,
      arg_name = arg_name,
      descriptor = descriptor,
      FUN.VALUE = character(1)
    )
    if (anyDuplicated(canonical_names)) {
      dup <- unique(canonical_names[duplicated(canonical_names)])
      stop(
        "`", arg_name, "` supplies the same design variable more than once: ",
        paste(dup, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    parsed_top_level <- stats::setNames(unname(design), canonical_names)
  }

  overlap <- intersect(names(parsed_top_level), names(parsed_future_branch))
  if (length(overlap) > 0L) {
    stop(
      "`", arg_name, "` supplies the same design variable through both top-level names and `",
      arg_name,
      "$facets`: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  c(parsed_top_level, parsed_future_branch)
}

simulation_resolve_design_counts <- function(sim_spec = NULL,
                                             n_person = NULL,
                                             n_rater = NULL,
                                             n_criterion = NULL,
                                             raters_per_person = NULL,
                                             design = NULL,
                                             defaults = NULL,
                                             design_arg = "design",
                                             explicit_scalar_names = NULL) {
  aliases <- simulation_design_variable_aliases(sim_spec)
  descriptor <- simulation_design_descriptor(sim_spec)
  parsed_design <- simulation_parse_design_input(
    design = design,
    aliases = aliases,
    descriptor = descriptor,
    future_facet_table = simulation_future_facet_table(sim_spec),
    future_branch_schema = simulation_future_branch_schema(sim_spec),
    arg_name = design_arg
  )

  explicit_scalars <- list(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person
  )
  if (is.null(explicit_scalar_names)) {
    explicit_scalars <- explicit_scalars[!vapply(explicit_scalars, is.null, logical(1))]
  } else {
    explicit_scalar_names <- intersect(names(explicit_scalars), explicit_scalar_names)
    explicit_scalars <- explicit_scalars[explicit_scalar_names]
    explicit_scalars <- explicit_scalars[!vapply(explicit_scalars, is.null, logical(1))]
  }
  overlap <- intersect(names(parsed_design), names(explicit_scalars))
  if (length(overlap) > 0L) {
    stop(
      "Do not supply the same design variable through both scalar arguments and `",
      design_arg,
      "`: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  values <- as.list(defaults %||% list())
  values[names(parsed_design)] <- parsed_design
  values[names(explicit_scalars)] <- explicit_scalars

  required <- c("n_person", "n_rater", "n_criterion")
  missing_required <- required[!required %in% names(values)]
  if (length(missing_required) > 0L) {
    stop(
      "Missing required design counts: ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!"raters_per_person" %in% names(values) || is.null(values$raters_per_person)) {
    values$raters_per_person <- values$n_rater
  }

  out <- tibble::tibble(
    n_person = simulation_validate_count(values$n_person, "n_person", min_value = 2L),
    n_rater = simulation_validate_count(values$n_rater, "n_rater", min_value = 2L),
    n_criterion = simulation_validate_count(values$n_criterion, "n_criterion", min_value = 2L),
    raters_per_person = simulation_validate_count(values$raters_per_person, "raters_per_person", min_value = 1L)
  )
  if (out$raters_per_person > out$n_rater) {
    stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  }
  out
}

simulation_validate_count_values <- function(x, arg_name, min_value = 1L) {
  values <- as.integer(x)
  if (length(values) < 1L || any(!is.finite(values)) || any(values < min_value)) {
    stop("`", arg_name, "` must contain integer values >= ", min_value, ".", call. = FALSE)
  }
  values
}

simulation_resolve_design_grid_values <- function(sim_spec = NULL,
                                                  n_person = NULL,
                                                  n_rater = NULL,
                                                  n_criterion = NULL,
                                                  raters_per_person = NULL,
                                                  design = NULL,
                                                  defaults = NULL,
                                                  design_arg = "design",
                                                  explicit_scalar_names = NULL) {
  aliases <- simulation_design_variable_aliases(sim_spec)
  descriptor <- simulation_design_descriptor(sim_spec)
  parsed_design <- simulation_parse_design_input(
    design = design,
    aliases = aliases,
    descriptor = descriptor,
    future_facet_table = simulation_future_facet_table(sim_spec),
    future_branch_schema = simulation_future_branch_schema(sim_spec),
    arg_name = design_arg
  )

  if (is.null(explicit_scalar_names)) {
    overlap_candidates <- c("n_person", "n_rater", "n_criterion", "raters_per_person")
  } else {
    overlap_candidates <- intersect(
      c("n_person", "n_rater", "n_criterion", "raters_per_person"),
      explicit_scalar_names
    )
  }
  overlap <- intersect(names(parsed_design), overlap_candidates)
  if (length(overlap) > 0L) {
    stop(
      "Do not supply the same design variable through both scalar arguments and `",
      design_arg,
      "`: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  values <- as.list(defaults %||% list(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person
  ))
  values[names(parsed_design)] <- parsed_design

  required <- c("n_person", "n_rater", "n_criterion")
  missing_required <- required[!required %in% names(values)]
  if (length(missing_required) > 0L) {
    stop(
      "Missing required design values: ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!"raters_per_person" %in% names(values) || is.null(values$raters_per_person)) {
    values$raters_per_person <- values$n_rater
  }

  list(
    n_person = simulation_validate_count_values(values$n_person, "n_person", min_value = 2L),
    n_rater = simulation_validate_count_values(values$n_rater, "n_rater", min_value = 2L),
    n_criterion = simulation_validate_count_values(values$n_criterion, "n_criterion", min_value = 2L),
    raters_per_person = simulation_validate_count_values(values$raters_per_person, "raters_per_person", min_value = 1L)
  )
}

simulation_design_canonical_variables <- function(descriptor = NULL) {
  if (is.data.frame(descriptor) && "canonical" %in% names(descriptor)) {
    vals <- as.character(descriptor$canonical)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) > 0L) {
      return(vals)
    }
  }
  c("n_person", "n_rater", "n_criterion", "raters_per_person")
}

simulation_design_group_variables <- function(descriptor = NULL, include_design_id = TRUE) {
  vars <- simulation_design_canonical_variables(descriptor)
  if (isTRUE(include_design_id)) {
    c("design_id", vars)
  } else {
    vars
  }
}

simulation_append_design_alias_columns <- function(tbl, aliases) {
  tbl <- tibble::as_tibble(tbl)
  if (!is.character(aliases) || length(aliases) == 0L) {
    return(tbl)
  }
  for (canonical in names(aliases)) {
    alias <- as.character(aliases[[canonical]])
    if (!canonical %in% names(tbl)) next
    if (is.na(alias) || !nzchar(alias) || identical(alias, canonical) || alias %in% names(tbl)) next
    tbl[[alias]] <- tbl[[canonical]]
  }
  tbl
}

simulation_build_design_grid <- function(n_person,
                                         n_rater,
                                         n_criterion,
                                         raters_per_person,
                                         design = NULL,
                                         sim_spec = NULL,
                                         id_prefix = "D",
                                         explicit_scalar_names = NULL) {
  aliases <- simulation_design_variable_aliases(sim_spec)
  descriptor <- simulation_design_descriptor(sim_spec)
  canonical_order <- as.character(descriptor$canonical)
  design_values <- simulation_resolve_design_grid_values(
    sim_spec = sim_spec,
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    defaults = list(
      n_person = n_person,
      n_rater = n_rater,
      n_criterion = n_criterion,
      raters_per_person = raters_per_person
    ),
    design_arg = "design",
    explicit_scalar_names = explicit_scalar_names
  )

  design_grid <- expand.grid(
    n_person = design_values$n_person,
    n_rater = design_values$n_rater,
    n_criterion = design_values$n_criterion,
    raters_per_person = design_values$raters_per_person,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design_grid <- design_grid[design_grid$raters_per_person <= design_grid$n_rater, , drop = FALSE]
  if (nrow(design_grid) == 0) {
    msg <- "No valid design rows remain after enforcing `raters_per_person <= n_rater`."
    public_rule <- paste0("`", unname(aliases[["raters_per_person"]]), " <= ", unname(aliases[["n_rater"]]), "`")
    if (!identical(public_rule, "`raters_per_person <= n_rater`")) {
      msg <- paste0(msg, " Public aliases: ", public_rule, ".")
    }
    stop(msg, call. = FALSE)
  }
  design_grid$design_id <- sprintf("%s%02d", as.character(id_prefix[1]), seq_len(nrow(design_grid)))
  design_grid <- tibble::as_tibble(design_grid[, c("design_id", canonical_order), drop = FALSE])

  list(
    canonical = design_grid,
    public = simulation_append_design_alias_columns(design_grid, aliases),
    aliases = aliases,
    descriptor = descriptor
  )
}

simulation_build_ademp <- function(purpose,
                                   design_grid,
                                   generator_model,
                                   generator_step_facet,
                                   generator_assignment,
                                   sim_spec = NULL,
                                   estimands,
                                   analysis_methods,
                                   performance_measures) {
  threshold_info <- simulation_threshold_summary(sim_spec)
  design_variable_aliases <- simulation_design_variable_aliases(sim_spec)
  design_descriptor <- simulation_design_descriptor(sim_spec)
  planning_scope <- simulation_planning_scope(sim_spec)
  planning_constraints <- simulation_planning_constraints(sim_spec)
  planning_schema <- simulation_planning_schema(sim_spec)
  design_variables <- simulation_design_canonical_variables(design_descriptor)
  list(
    aims = purpose,
    data_generating_mechanism = list(
      source = if (is.null(sim_spec)) "scalar_arguments" else as.character(sim_spec$source %||% "mfrm_sim_spec"),
      model = generator_model,
      step_facet = generator_step_facet,
      assignment = generator_assignment,
      latent_distribution = if (is.null(sim_spec)) "normal" else as.character(sim_spec$latent_distribution %||% "normal"),
      threshold_mode = threshold_info$mode,
      threshold_step_facet = threshold_info$step_facet,
      design_variables = design_variables,
      design_variable_aliases = design_variable_aliases,
      design_descriptor = design_descriptor,
      planning_scope = planning_scope,
      planning_constraints = planning_constraints,
      planning_schema = planning_schema
    ),
    estimands = estimands,
    methods = analysis_methods,
    performance_measures = performance_measures
  )
}

design_eval_recovery_metrics <- function(est_levels, est_values, truth_vec) {
  idx <- match(as.character(est_levels), names(truth_vec))
  ok <- is.finite(idx) & is.finite(est_values)
  if (!any(ok)) {
    return(list(
      raw_rmse = NA_real_,
      raw_bias = NA_real_,
      aligned_rmse = NA_real_,
      aligned_bias = NA_real_
    ))
  }

  truth_vals <- as.numeric(truth_vec[idx[ok]])
  est_vals <- as.numeric(est_values[ok])
  diffs <- est_vals - truth_vals
  shift <- mean(diffs, na.rm = TRUE)
  aligned_diffs <- diffs - shift

  list(
    raw_rmse = sqrt(mean(diffs^2, na.rm = TRUE)),
    raw_bias = mean(diffs, na.rm = TRUE),
    aligned_rmse = sqrt(mean(aligned_diffs^2, na.rm = TRUE)),
    aligned_bias = mean(aligned_diffs, na.rm = TRUE)
  )
}

design_eval_match_metric <- function(metric) {
  switch(
    metric,
    separation = "MeanSeparation",
    reliability = "MeanReliability",
    infit = "MeanInfit",
    outfit = "MeanOutfit",
    misfitrate = "MeanMisfitRate",
    severityrmse = "MeanSeverityRMSE",
    severitybias = "MeanSeverityBias",
    convergencerate = "ConvergenceRate",
    elapsedsec = "MeanElapsedSec",
    mincategorycount = "MeanMinCategoryCount",
    stop("Unknown metric: ", metric, call. = FALSE)
  )
}

design_eval_build_notes <- function(summary_tbl) {
  notes <- character(0)
  if (!is.data.frame(summary_tbl) || nrow(summary_tbl) == 0) return(notes)
  if (any(summary_tbl$ConvergenceRate < 1, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions did not converge in every replication.")
  }
  if (any(summary_tbl$MeanMinCategoryCount < 10, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions produced sparse score categories (< 10 observations).")
  }
  if (any(summary_tbl$MeanSeparation < 2, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions yielded low facet separation (< 2.0).")
  }
  if (any(grepl("^Mcse", names(summary_tbl)))) {
    notes <- c(notes, "MCSE columns summarize finite-replication uncertainty around the reported means and rates.")
  }
  if ("RecoveryComparableRate" %in% names(summary_tbl) &&
      any(summary_tbl$RecoveryComparableRate < 1, na.rm = TRUE)) {
    notes <- c(
      notes,
      "Recovery metrics are reported only for design rows where generator and fitted model contracts align; rows with generator-fit mismatches set recovery fields to NA."
    )
  }
  notes
}

design_eval_summarize_results <- function(results, rep_overview, design_variable_aliases = NULL) {
  results_tbl <- tibble::as_tibble(results)
  rep_tbl <- tibble::as_tibble(rep_overview)
  design_descriptor <- simulation_object_design_descriptor(list(results = results_tbl, rep_overview = rep_tbl, design_variable_aliases = design_variable_aliases))
  design_vars <- simulation_design_canonical_variables(design_descriptor)
  grouping_vars <- c("design_id", "Facet", design_vars)
  arrange_vars <- c("Facet", design_vars)

  design_summary <- tibble::tibble()
  if (nrow(results_tbl) > 0) {
    design_summary <- results_tbl |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) |>
      dplyr::summarize(
        Reps = dplyr::n(),
        ConvergenceRate = mean(.data$Converged, na.rm = TRUE),
        McseConvergenceRate = simulation_mcse_proportion(.data$Converged),
        MeanSeparation = mean(.data$Separation, na.rm = TRUE),
        SdSeparation = design_eval_safe_sd(.data$Separation),
        McseSeparation = simulation_mcse_mean(.data$Separation),
        MeanReliability = mean(.data$Reliability, na.rm = TRUE),
        McseReliability = simulation_mcse_mean(.data$Reliability),
        MeanInfit = mean(.data$MeanInfit, na.rm = TRUE),
        McseInfit = simulation_mcse_mean(.data$MeanInfit),
        MeanOutfit = mean(.data$MeanOutfit, na.rm = TRUE),
        McseOutfit = simulation_mcse_mean(.data$MeanOutfit),
        MeanMisfitRate = mean(.data$MisfitRate, na.rm = TRUE),
        McseMisfitRate = simulation_mcse_mean(.data$MisfitRate),
        MeanSeverityRMSE = mean(.data$SeverityRMSE, na.rm = TRUE),
        McseSeverityRMSE = simulation_mcse_mean(.data$SeverityRMSE),
        MeanSeverityBias = mean(.data$SeverityBias, na.rm = TRUE),
        McseSeverityBias = simulation_mcse_mean(.data$SeverityBias),
        MeanSeverityRMSERaw = mean(.data$SeverityRMSERaw, na.rm = TRUE),
        McseSeverityRMSERaw = simulation_mcse_mean(.data$SeverityRMSERaw),
        MeanSeverityBiasRaw = mean(.data$SeverityBiasRaw, na.rm = TRUE),
        McseSeverityBiasRaw = simulation_mcse_mean(.data$SeverityBiasRaw),
        RecoveryComparableRate = mean(.data$RecoveryComparable, na.rm = TRUE),
        RecoveryBasis = paste(sort(unique(as.character(.data$RecoveryBasis))), collapse = "; "),
        MeanElapsedSec = mean(.data$ElapsedSec, na.rm = TRUE),
        McseElapsedSec = simulation_mcse_mean(.data$ElapsedSec),
        MeanMinCategoryCount = mean(.data$MinCategoryCount, na.rm = TRUE),
        McseMinCategoryCount = simulation_mcse_mean(.data$MinCategoryCount),
        .groups = "drop"
      ) |>
      dplyr::arrange(!!!rlang::syms(arrange_vars))
    design_summary <- simulation_append_design_alias_columns(design_summary, design_variable_aliases)
  }

  overview <- tibble::tibble(
    Designs = dplyr::n_distinct(rep_tbl$design_id),
    Replications = nrow(rep_tbl),
    SuccessfulRuns = sum(rep_tbl$RunOK, na.rm = TRUE),
    ConvergedRuns = sum(rep_tbl$Converged, na.rm = TRUE),
    MeanElapsedSec = design_eval_safe_mean(rep_tbl$ElapsedSec)
  )

  list(
    overview = overview,
    design_summary = design_summary,
    notes = design_eval_build_notes(design_summary)
  )
}

#' Evaluate MFRM design conditions by repeated simulation
#'
#' @param n_person Vector of person counts to evaluate.
#' @param n_rater Vector of rater counts to evaluate.
#' @param n_criterion Vector of criterion counts to evaluate.
#' @param raters_per_person Vector of rater assignments per person.
#' @param design Optional named design-grid override supplied as a named list,
#'   named vector, or one-row data frame. Names may use canonical variables
#'   (`n_person`, `n_rater`, `n_criterion`, `raters_per_person`), current
#'   public aliases implied by `sim_spec` (for example `n_judge`, `n_task`,
#'   `judge_per_person`), or role keywords (`person`, `rater`, `criterion`,
#'   `assignment`). Values may be vectors. The schema-only future branch input
#'   `design$facets = c(person = ..., judge = ..., task = ...)` is also
#'   accepted for the currently exposed facet keys. Do not specify the same
#'   variable through both `design` and the scalar design-grid arguments.
#' @param reps Number of replications per design condition.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread of step thresholds on the logit scale.
#' @param fit_method Estimation method passed to [fit_mfrm()].
#' @param model Measurement model passed to [fit_mfrm()]. The current design
#'   evaluator supports `RSM` and `PCM`; bounded `GPCM` is accepted only
#'   to produce an explicit unsupported-path error.
#' @param step_facet Step facet passed to [fit_mfrm()] when `model = "PCM"`.
#'   When left `NULL`, the function inherits the generator step facet from
#'   `sim_spec` when available and otherwise defaults to `"Criterion"`.
#' @param maxit Maximum iterations passed to [fit_mfrm()].
#' @param quad_points Quadrature points for `fit_method = "MML"`.
#' @param residual_pca Residual PCA mode passed to [diagnose_mfrm()].
#' @param sim_spec Optional output from [build_mfrm_sim_spec()] or
#'   [extract_mfrm_sim_spec()] used as the base data-generating mechanism.
#'   When supplied, the design grid still varies `n_person`, `n_rater`,
#'   `n_criterion`, and `raters_per_person`, but latent-spread assumptions,
#'   thresholds, and other generator settings come from `sim_spec`.
#'   If `sim_spec` contains step-facet-specific thresholds, the design grid may
#'   not vary the number of levels for that step facet away from the
#'   specification. If `sim_spec` stores an active latent-regression population
#'   generator, this helper currently requires `fit_method = "MML"` so each
#'   replication can refit the population model.
#' @param seed Optional seed for reproducible replications.
#'
#' @details
#' This helper runs a compact Monte Carlo design study for common rater-by-item
#' many-facet settings.
#'
#' For each design condition, the function:
#' 1. generates synthetic data with [simulate_mfrm_data()]
#' 2. fits the requested MFRM with [fit_mfrm()]
#' 3. computes diagnostics with [diagnose_mfrm()]
#' 4. stores recovery and precision summaries by facet
#'
#' The result is intended for planning questions such as:
#' - how many raters are needed for stable rater separation?
#' - how does `raters_per_person` affect severity recovery?
#' - when do category counts become too sparse for comfortable interpretation?
#'
#' This is a **parametric simulation study**. It does not take one observed
#' design (for example, 4 raters x 30 persons x 3 criteria) and analytically
#' extrapolate what would happen under a different design (for example,
#' 2 raters x 40 persons x 5 criteria). Instead, you specify a design grid and
#' data-generating assumptions (latent spread, facet spread, thresholds, noise,
#' and scoring structure), and the function repeatedly generates synthetic data
#' under those assumptions.
#'
#' When you want the simulated conditions to resemble an existing study, use
#' substantive knowledge or estimates from that study to choose
#' `theta_sd`, `rater_sd`, `criterion_sd`, `score_levels`, and related
#' settings before running the design evaluation.
#'
#' When `sim_spec` is supplied, the function uses it as the explicit
#' data-generating mechanism. This is the recommended route when you want a
#' design study to stay close to a previously fitted run while still varying the
#' candidate sample sizes or rater-assignment counts.
#'
#' If that specification also stores a latent-regression population generator,
#' each replication carries forward the simulated one-row-per-person background
#' data and refits the MML population-model branch. This remains a scenario
#' study under explicit assumptions; it is not a closed-form predictive
#' distribution for one future administration.
#'
#' First-release `GPCM` is not yet available in this design-evaluation helper.
#' The missing pieces are not just software wiring: the current package still
#' needs a validated slope-generating simulation contract and downstream
#' diagnostics compatible with the generalized ordered kernel. More broadly,
#' the current planning layer is still role-based for exactly two non-person
#' facets (`rater`-like and `criterion`-like), even though the estimation core
#' supports arbitrary facet counts.
#'
#' Recovery metrics are reported only when the generator and fitted model target
#' the same facet-parameter contract. In practice this means the same
#' `model`, and for `PCM`, the same `step_facet`. When these do not align,
#' recovery fields are set to `NA` and the output records the reason. Even when
#' these contract checks pass, the recovery summaries still assume compatible
#' orientation and anchoring conventions across the generator and fitted model.
#'
#' @section Reported metrics:
#' Facet-level simulation results include:
#' - `Separation` (\eqn{G = \mathrm{SD_{adj}} / \mathrm{RMSE}}):
#'   how many statistically distinct strata the facet resolves.
#' - `Reliability` (\eqn{G^2 / (1 + G^2)}): analogous to Cronbach's
#'   \eqn{\alpha} for the reproducibility of element ordering.
#' - `Strata` (\eqn{(4G + 1) / 3}): number of distinguishable groups.
#' - Mean `Infit` and `Outfit`: average fit mean-squares across elements.
#' - `MisfitRate`: share of elements with \eqn{|\mathrm{ZSTD}| > 2}.
#' - `SeverityRMSE`: root-mean-square error of recovered parameters vs
#'   the known truth **after facet-wise mean alignment**, so that the
#'   usual Rasch/MFRM location indeterminacy does not inflate recovery
#'   error. This quantity is reported only when the generator and fitted model
#'   target the same facet-parameter contract.
#' - `SeverityBias`: mean signed recovery error after the same alignment;
#'   values near zero are expected. This is likewise omitted when the
#'   generator/fitted-model contract does not align.
#'
#' @section Interpreting output:
#' Start with `summary(x)$design_summary`, then plot one focal metric at a time
#' (for example rater `Separation` or criterion `SeverityRMSE`).
#'
#' Higher separation/reliability is generally better, whereas lower
#' `SeverityRMSE`, `MeanMisfitRate`, and `MeanElapsedSec` are preferable.
#'
#' When choosing among designs, look for the point where increasing
#' `n_person` or `raters_per_person` yields diminishing returns in
#' separation and RMSE---this identifies the cost-effective design
#' frontier.  `ConvergedRuns / reps` should be near 1.0; low
#' convergence rates indicate the design is too small for the chosen
#' estimation method.
#'
#' @section References:
#' The simulation logic follows the general Monte Carlo / operating-characteristic
#' framework described by Morris, White, and Crowther (2019) and the
#' ADEMP-oriented planning/reporting guidance summarized for psychology by
#' Siepe et al. (2024). In `mfrmr`, `evaluate_mfrm_design()` is a practical
#' many-facet design-planning wrapper rather than a direct reproduction of one
#' published simulation study.
#'
#' - Morris, T. P., White, I. R., & Crowther, M. J. (2019).
#'   *Using simulation studies to evaluate statistical methods*.
#'   Statistics in Medicine, 38(11), 2074-2102.
#' - Siepe, B. S., Bartos, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W.,
#'   & Pawel, S. (2024). *Simulation studies for methodological research in
#'   psychology: A standardized template for planning, preregistration, and
#'   reporting*. Psychological Methods.
#'
#' @return An object of class `mfrm_design_evaluation` with components:
#' - `design_grid`: evaluated design conditions. When `sim_spec` carries custom
#'   public facet names, matching design-variable alias columns are included
#'   alongside the canonical internal columns.
#' - `results`: facet-level replicate results, with the same design-variable
#'   alias columns when applicable.
#' - `rep_overview`: run-level status and timing, with the same design-variable
#'   alias columns when applicable.
#' - `design_descriptor`: role-based design-variable metadata used by planning
#'   summaries and plots
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of which design variables remain
#'   mutable under the current simulation specification
#' - `planning_schema`: combined planner-schema contract bundling the role
#'   descriptor, scope boundary, and current mutability map
#' - `settings`: simulation settings
#' - `ademp`: simulation-study metadata (aims, DGM, estimands, methods, performance measures)
#' @seealso [simulate_mfrm_data()], [summary.mfrm_design_evaluation], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- suppressWarnings(evaluate_mfrm_design(
#'   design = list(person = c(8, 12), rater = 2, criterion = 2, assignment = 1),
#'   reps = 1,
#'   maxit = 8,
#'   seed = 123
#' ))
#' s_eval <- summary(sim_eval)
#' s_eval$design_summary[, c("Facet", "n_person", "MeanSeparation", "MeanSeverityRMSE")]
#' p_eval <- plot(sim_eval, facet = "Rater", metric = "separation", x_var = "n_person", draw = FALSE)
#' names(p_eval)
#' }
#' @export
evaluate_mfrm_design <- function(n_person = c(30, 50, 100),
                                 n_rater = c(3, 5),
                                 n_criterion = c(3, 5),
                                 raters_per_person = n_rater,
                                 design = NULL,
                                 reps = 10,
                                 score_levels = 4,
                                 theta_sd = 1,
                                 rater_sd = 0.35,
                                 criterion_sd = 0.25,
                                 noise_sd = 0,
                                 step_span = 1.4,
                                 fit_method = c("JML", "MML"),
                                 model = c("RSM", "PCM", "GPCM"),
                                 step_facet = NULL,
                                 maxit = 25,
                                 quad_points = 7,
                                 residual_pca = c("none", "overall", "facet", "both"),
                                 sim_spec = NULL,
                                 seed = NULL) {
  fit_method <- match.arg(fit_method)
  model <- match.arg(model)
  if (identical(model, "GPCM")) {
    stop(
      "`evaluate_mfrm_design()` does not yet support bounded `GPCM`. ",
      "Design evaluation still depends on simulation and diagnostics layers that remain validated only for `RSM` / `PCM`. ",
      gpcm_planning_scope_rationale(),
      call. = FALSE
    )
  }
  residual_pca <- match.arg(residual_pca)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!is.null(sim_spec) && identical(as.character(sim_spec$model %||% NA_character_), "GPCM")) {
    stop(
      "`evaluate_mfrm_design()` does not yet support bounded `GPCM` simulation specifications. ",
      "Direct data generation is available, but design evaluation still depends on diagnostics and recovery layers validated only for `RSM` / `PCM`. ",
      gpcm_planning_scope_rationale(),
      call. = FALSE
    )
  }
  if (!is.null(sim_spec) &&
      isTRUE((sim_spec$population %||% simulation_empty_population_spec())$active) &&
      !identical(fit_method, "MML")) {
    stop(
      "Simulation specifications with an active latent-regression population generator currently require `fit_method = \"MML\"`.",
      call. = FALSE
    )
  }
  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) stop("`reps` must be >= 1.", call. = FALSE)
  supplied_counts <- intersect(
    names(as.list(match.call(expand.dots = FALSE))[-1]),
    c("n_person", "n_rater", "n_criterion", "raters_per_person")
  )

  design_meta <- simulation_build_design_grid(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    sim_spec = sim_spec,
    id_prefix = "D",
    explicit_scalar_names = supplied_counts
  )
  design_grid <- design_meta$canonical
  design_grid_public <- design_meta$public

  generator_model <- if (is.null(sim_spec)) model else sim_spec$model
  base_facet_names <- if (is.null(sim_spec)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(sim_spec)
  generator_step_facet <- if (is.null(sim_spec)) if (identical(generator_model, "PCM")) base_facet_names[2] else NA_character_ else sim_spec$step_facet
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
  design_variable_aliases <- design_meta$aliases
  design_descriptor <- design_meta$descriptor
  planning_scope <- simulation_planning_scope(sim_spec)
  planning_constraints <- simulation_planning_constraints(sim_spec)
  planning_schema <- simulation_planning_schema(sim_spec)
  fit_step_facet <- simulation_resolve_fit_step_facet(model, step_facet, generator_step_facet)
  recovery_contract <- simulation_recovery_contract(
    generator_model = generator_model,
    generator_step_facet = generator_step_facet,
    fitted_model = model,
    fitted_step_facet = fit_step_facet
  )

  seeds <- with_preserved_rng_seed(
    seed,
    sample.int(.Machine$integer.max, size = nrow(design_grid) * reps, replace = FALSE)
  )
  seed_idx <- 0L

  result_rows <- vector("list", nrow(design_grid) * reps * 3L)
  rep_rows <- vector("list", nrow(design_grid) * reps)
  result_idx <- 0L
  rep_idx <- 0L

  for (i in seq_len(nrow(design_grid))) {
    design <- design_grid[i, , drop = FALSE]
    row_spec <- if (is.null(sim_spec)) {
      NULL
    } else {
      simulation_override_spec_design(
        sim_spec,
        n_person = design$n_person,
        n_rater = design$n_rater,
        n_criterion = design$n_criterion,
        raters_per_person = design$raters_per_person
      )
    }
    row_score_levels <- if (is.null(row_spec)) score_levels else row_spec$score_levels
    row_facet_names <- if (is.null(row_spec)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(row_spec)
    for (rep in seq_len(reps)) {
      seed_idx <- seed_idx + 1L
      sim <- if (is.null(row_spec)) {
        simulate_mfrm_data(
          n_person = design$n_person,
          n_rater = design$n_rater,
          n_criterion = design$n_criterion,
          raters_per_person = design$raters_per_person,
          score_levels = score_levels,
          theta_sd = theta_sd,
          rater_sd = rater_sd,
          criterion_sd = criterion_sd,
          noise_sd = noise_sd,
          step_span = step_span,
          seed = seeds[seed_idx]
        )
      } else {
        simulate_mfrm_data(sim_spec = row_spec, seed = seeds[seed_idx])
      }

      t0 <- proc.time()[["elapsed"]]
      fit_args <- list(
        data = sim,
        person = "Person",
        facets = row_facet_names,
        score = "Score",
        method = fit_method,
        model = model,
        maxit = maxit
      )
      if (identical(model, "PCM")) fit_args$step_facet <- fit_step_facet
      if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
      if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points
      sim_population <- attr(sim, "mfrm_population_data")
      if (is.list(sim_population) && isTRUE(sim_population$active)) {
        fit_args$population_formula <- sim_population$population_formula
        fit_args$person_data <- sim_population$person_data
        fit_args$person_id <- sim_population$person_id
        fit_args$population_policy <- sim_population$population_policy
      }

      fit <- tryCatch(do.call(fit_mfrm, fit_args), error = function(e) e)
      diag <- if (inherits(fit, "error")) fit else {
        tryCatch(
          diagnose_mfrm(fit, residual_pca = residual_pca),
          error = function(e) e
        )
      }
      elapsed <- proc.time()[["elapsed"]] - t0
      truth <- attr(sim, "mfrm_truth")
      min_category_count <- min(tabulate(sim$Score, nbins = row_score_levels))

      rep_idx <- rep_idx + 1L
      rep_row <- tibble::tibble(
        design_id = design$design_id,
        rep = rep,
        n_person = design$n_person,
        n_rater = design$n_rater,
        n_criterion = design$n_criterion,
        raters_per_person = design$raters_per_person,
        Observations = nrow(sim),
        MinCategoryCount = min_category_count,
        ElapsedSec = elapsed,
        RunOK = FALSE,
        Converged = FALSE,
        Error = NA_character_
      )

      if (inherits(fit, "error")) {
        rep_row$Error <- conditionMessage(fit)
        rep_rows[[rep_idx]] <- rep_row
        next
      }
      if (inherits(diag, "error")) {
        rep_row$Error <- conditionMessage(diag)
        rep_rows[[rep_idx]] <- rep_row
        next
      }

      converged <- isTRUE(as.logical(fit$summary$Converged[1]))
      rep_row$RunOK <- TRUE
      rep_row$Converged <- converged
      rep_rows[[rep_idx]] <- rep_row

      reliability_tbl <- tibble::as_tibble(diag$reliability)
      fit_tbl <- tibble::as_tibble(diag$fit)
      measure_tbl <- tibble::as_tibble(diag$measures)
      reliability_split <- split(reliability_tbl, reliability_tbl$Facet)
      fit_split <- split(fit_tbl, fit_tbl$Facet)
      measure_split <- split(measure_tbl, measure_tbl$Facet)

      for (facet in names(reliability_split)) {
        rel_row <- tibble::as_tibble(reliability_split[[facet]])
        facet_fit <- tibble::as_tibble(fit_split[[facet]] %||% tibble::tibble())
        facet_meas <- tibble::as_tibble(measure_split[[facet]] %||% tibble::tibble())
        truth_vec <- design_eval_extract_truth(truth, facet)

        severity_rmse <- NA_real_
        severity_bias <- NA_real_
        severity_rmse_raw <- NA_real_
        severity_bias_raw <- NA_real_
        if (isTRUE(recovery_contract$comparable) &&
            !is.null(truth_vec) && nrow(facet_meas) > 0 &&
            "Level" %in% names(facet_meas) && "Estimate" %in% names(facet_meas)) {
          recovery <- design_eval_recovery_metrics(
            est_levels = facet_meas$Level,
            est_values = suppressWarnings(as.numeric(facet_meas$Estimate)),
            truth_vec = truth_vec
          )
          severity_rmse <- recovery$aligned_rmse
          severity_bias <- recovery$aligned_bias
          severity_rmse_raw <- recovery$raw_rmse
          severity_bias_raw <- recovery$raw_bias
        }

        misfit_rate <- NA_real_
        if (nrow(facet_fit) > 0) {
          z_in <- suppressWarnings(as.numeric(facet_fit$InfitZSTD))
          z_out <- suppressWarnings(as.numeric(facet_fit$OutfitZSTD))
          misfit_rate <- mean(abs(z_in) > 2 | abs(z_out) > 2, na.rm = TRUE)
        }

        result_idx <- result_idx + 1L
        result_rows[[result_idx]] <- tibble::tibble(
          design_id = design$design_id,
          rep = rep,
          Facet = facet,
          n_person = design$n_person,
          n_rater = design$n_rater,
          n_criterion = design$n_criterion,
          raters_per_person = design$raters_per_person,
          Observations = nrow(sim),
          MinCategoryCount = min_category_count,
          ElapsedSec = elapsed,
          Converged = converged,
          GeneratorModel = generator_model,
          GeneratorStepFacet = generator_step_facet,
          FitModel = model,
          FitStepFacet = fit_step_facet,
          RecoveryComparable = recovery_contract$comparable,
          RecoveryBasis = recovery_contract$basis,
          Levels = suppressWarnings(as.integer(rel_row$Levels[1])),
          Separation = suppressWarnings(as.numeric(rel_row$Separation[1])),
          Strata = suppressWarnings(as.numeric(rel_row$Strata[1])),
          Reliability = suppressWarnings(as.numeric(rel_row$Reliability[1])),
          MeanInfit = suppressWarnings(as.numeric(rel_row$MeanInfit[1])),
          MeanOutfit = suppressWarnings(as.numeric(rel_row$MeanOutfit[1])),
          MisfitRate = misfit_rate,
          SeverityRMSE = severity_rmse,
          SeverityBias = severity_bias,
          SeverityRMSERaw = severity_rmse_raw,
          SeverityBiasRaw = severity_bias_raw
        )
      }
    }
  }

  results <- dplyr::bind_rows(result_rows[seq_len(result_idx)])
  rep_overview <- dplyr::bind_rows(rep_rows[seq_len(rep_idx)])
  results <- simulation_append_design_alias_columns(results, design_variable_aliases)
  rep_overview <- simulation_append_design_alias_columns(rep_overview, design_variable_aliases)
  ademp <- simulation_build_ademp(
    purpose = "Assess many-facet design conditions via repeated parametric simulation under explicit data-generating assumptions.",
    design_grid = design_grid,
    generator_model = generator_model,
    generator_step_facet = generator_step_facet,
    generator_assignment = generator_assignment,
    sim_spec = sim_spec,
    estimands = c(
      "Facet separation, reliability, and strata",
      "Mean infit/outfit and misfit rate",
      "Aligned facet-parameter recovery RMSE and bias",
      "Convergence rate and elapsed time"
    ),
    analysis_methods = list(
      fit_method = fit_method,
      fitted_model = model,
      maxit = maxit,
      quad_points = if (identical(fit_method, "MML")) quad_points else NA_integer_,
      residual_pca = residual_pca
    ),
    performance_measures = c(
      "Mean performance across replications",
      "MCSE for means and rates",
      "Convergence rate",
      "Sparse-category warning rate"
    )
  )

  structure(
    list(
      design_grid = design_grid_public,
      results = results,
      rep_overview = rep_overview,
      design_descriptor = design_descriptor,
      planning_scope = planning_scope,
      planning_constraints = planning_constraints,
      planning_schema = planning_schema,
      settings = list(
        reps = reps,
        score_levels = score_levels,
        theta_sd = theta_sd,
        rater_sd = rater_sd,
        criterion_sd = criterion_sd,
        noise_sd = noise_sd,
        step_span = step_span,
        fit_method = fit_method,
        model = model,
        step_facet = fit_step_facet,
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        sim_spec = sim_spec,
        facet_names = stats::setNames(base_facet_names, c("rater", "criterion")),
        design_variable_aliases = design_variable_aliases,
        design_descriptor = design_descriptor,
        planning_scope = planning_scope,
        planning_constraints = planning_constraints,
        planning_schema = planning_schema,
        generator_model = generator_model,
        generator_step_facet = generator_step_facet,
        generator_assignment = generator_assignment,
        recovery_comparable = recovery_contract$comparable,
        recovery_basis = recovery_contract$basis,
        seed = seed
      ),
      ademp = ademp
    ),
    class = "mfrm_design_evaluation"
  )
}

#' Summarize a design-simulation study
#'
#' @param object Output from [evaluate_mfrm_design()].
#' @param digits Number of digits used in the returned numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' The summary emphasizes condition-level averages that are useful for practical
#' design planning, especially:
#' - convergence rate
#' - separation and reliability by facet
#' - severity recovery RMSE
#' - mean misfit rate
#'
#' @return An object of class `summary.mfrm_design_evaluation` with components:
#' - `overview`: run-level overview
#' - `design_summary`: aggregated design-by-facet metrics, with design-variable
#'   alias columns when applicable
#' - `ademp`: simulation-study metadata carried forward from the original object
#' - `facet_names`: public facet labels carried from the simulation specification
#' - `design_variable_aliases`: accepted public aliases for design variables
#' - `design_descriptor`: role-based design-variable metadata
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract
#' - `future_branch_active_summary`: compact deterministic summary of the
#'   schema-only future arbitrary-facet planning branch embedded in the current
#'   planning schema
#' - `notes`: short interpretation notes
#' @seealso [evaluate_mfrm_design()], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- suppressWarnings(evaluate_mfrm_design(
#'   n_person = c(8, 12),
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 8,
#'   seed = 123
#' ))
#' s <- summary(sim_eval)
#' s$overview
#' head(s$design_summary)
#' }
#' @export
summary.mfrm_design_evaluation <- function(object, digits = 3, ...) {
  if (!is.list(object) || is.null(object$results) || is.null(object$rep_overview)) {
    stop("`object` must be output from evaluate_mfrm_design().")
  }
  digits <- max(0L, as.integer(digits[1]))
  out <- design_eval_summarize_results(
    object$results,
    object$rep_overview,
    design_variable_aliases = simulation_object_design_variable_aliases(object)
  )

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out$overview <- round_df(out$overview)
  out$design_summary <- round_df(out$design_summary)
  out$ademp <- object$ademp %||% NULL
  out$facet_names <- object$settings$facet_names %||% stats::setNames(simulation_default_output_facet_names(), c("rater", "criterion"))
  out$design_variable_aliases <- simulation_object_design_variable_aliases(object)
  out$design_descriptor <- simulation_object_design_descriptor(object)
  out$planning_scope <- simulation_object_planning_scope(object)
  out$planning_constraints <- simulation_object_planning_constraints(object)
  out$planning_schema <- simulation_object_planning_schema(object)
  out$future_branch_active_summary <- simulation_compact_future_branch_active_summary(
    object,
    digits = digits
  )
  out$digits <- digits
  scope_note <- simulation_planning_scope_note(out$planning_scope)
  if (length(scope_note) > 0L && !scope_note %in% out$notes) {
    out$notes <- c(out$notes, scope_note)
  }
  constraint_note <- simulation_planning_constraints_note(out$planning_constraints)
  if (length(constraint_note) > 0L && !constraint_note %in% out$notes) {
    out$notes <- c(out$notes, constraint_note)
  }
  schema_note <- simulation_planning_schema_note(out$planning_schema)
  if (length(schema_note) > 0L && !schema_note %in% out$notes) {
    out$notes <- c(out$notes, schema_note)
  }
  if (inherits(out$future_branch_active_summary, "summary.mfrm_future_branch_active_branch")) {
    out$notes <- c(
      out$notes,
      "A deterministic future arbitrary-facet planning scaffold is embedded in `future_branch_active_summary`; it reports structural bookkeeping and conservative recommendation logic, not Monte Carlo performance."
    )
  }
  out$notes <- unique(out$notes)
  class(out) <- "summary.mfrm_design_evaluation"
  out
}

#' @export
print.summary.mfrm_design_evaluation <- function(x, ...) {
  digits <- max(0L, as.integer(x$digits %||% 3L))
  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }
  preview_df <- function(df, n = 10L) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    utils::head(df, n = n)
  }

  cat("mfrmr Design Evaluation Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_df(as.data.frame(x$overview)), row.names = FALSE)
  }
  if (!is.null(x$design_summary) && nrow(x$design_summary) > 0) {
    cat("\nDesign summary (preview)\n")
    print(round_df(as.data.frame(preview_df(x$design_summary))), row.names = FALSE)
  }
  print_compact_future_branch_active_summary(
    x$future_branch_active_summary %||% NULL,
    digits = digits
  )
  if (is.list(x$ademp) && length(x$ademp) > 0L) {
    cat("\nADEMP metadata\n")
    cat(" - aims\n")
    cat(" - data_generating_mechanism\n")
    cat(" - estimands\n")
    cat(" - methods\n")
    cat(" - performance_measures\n")
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Plot a design-simulation study
#'
#' @param x Output from [evaluate_mfrm_design()].
#' @param facet Facet to visualize.
#' @param metric Metric to plot.
#' @param x_var Design variable used on the x-axis. When `x` was generated from
#'   a `sim_spec` with custom public facet names, the corresponding aliases
#'   (for example `n_judge`, `n_task`, `judge_per_person`) are also accepted.
#'   Role keywords (`person`, `rater`, `criterion`, `assignment`) are accepted
#'   as an abstraction over the current two-facet schema.
#' @param group_var Optional design variable used for separate lines. The same
#'   alias rules as `x_var` apply.
#' @param draw If `TRUE`, draw with base graphics; otherwise return plotting data.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method is designed for quick design-planning scans rather than polished
#' publication graphics.
#'
#' Useful first plots are:
#' - rater `metric = "separation"` against `x_var = "n_person"`
#' - criterion `metric = "severityrmse"` against `x_var = "n_person"`
#'   when you want aligned recovery error rather than raw location shifts
#' - rater `metric = "convergencerate"` against `x_var = "raters_per_person"`
#'
#' @return If `draw = TRUE`, invisibly returns a plotting-data list. If
#'   `draw = FALSE`, returns that list directly. The returned list includes
#'   resolved canonical variables (`x_var`, `group_var`) together with public
#'   labels (`x_label`, `group_label`), `design_variable_aliases`, and
#'   `design_descriptor`, plus `planning_scope`, `planning_constraints`, and
#'   `planning_schema`.
#' @seealso [evaluate_mfrm_design()], [summary.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- suppressWarnings(evaluate_mfrm_design(
#'   n_person = c(8, 12),
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 8,
#'   seed = 123
#' ))
#' p <- plot(sim_eval, facet = "Rater", metric = "separation", x_var = "n_person", draw = FALSE)
#' c(p$facet, p$x_var)
#' }
#' @export
plot.mfrm_design_evaluation <- function(x,
                                        facet = c("Rater", "Criterion", "Person"),
                                        metric = c("separation", "reliability", "infit", "outfit",
                                                   "misfitrate", "severityrmse", "severitybias",
                                                   "convergencerate", "elapsedsec", "mincategorycount"),
                                        x_var = c("n_person", "n_rater", "n_criterion", "raters_per_person"),
                                        group_var = NULL,
                                        draw = TRUE,
                                        ...) {
  if (!is.list(x) || is.null(x$results) || is.null(x$rep_overview)) {
    stop("`x` must be output from evaluate_mfrm_design().")
  }
  available_facets <- unique(as.character(tibble::as_tibble(x$results)$Facet))
  if (missing(facet)) {
    ordered_facets <- c(setdiff(available_facets, "Person"), intersect("Person", available_facets))
    facet <- ordered_facets[1]
  } else {
    facet <- as.character(facet[1])
    if (!facet %in% available_facets) {
      stop("`facet` must be one of: ", paste(available_facets, collapse = ", "), ".")
    }
  }
  metric <- match.arg(metric)
  design_variable_aliases <- simulation_object_design_variable_aliases(x)
  design_descriptor <- simulation_object_design_descriptor(x)
  x_var <- if (missing(x_var)) {
    "n_person"
  } else {
    simulation_resolve_design_variable(x_var, design_variable_aliases, "x_var", descriptor = design_descriptor)
  }
  x_label <- simulation_design_variable_label(x_var, design_variable_aliases)

  sum_obj <- design_eval_summarize_results(
    x$results,
    x$rep_overview,
    design_variable_aliases = design_variable_aliases
  )
  plot_tbl <- tibble::as_tibble(sum_obj$design_summary)
  if (nrow(plot_tbl) == 0) stop("No design-summary rows available for plotting.")
  plot_tbl <- plot_tbl[plot_tbl$Facet == facet, , drop = FALSE]
  if (nrow(plot_tbl) == 0) stop("No rows available for facet `", facet, "`.")

  metric_col <- design_eval_match_metric(metric)
  varying <- simulation_design_canonical_variables(design_descriptor)
  varying <- varying[varying != x_var]
  if (is.null(group_var)) {
    cand <- varying[vapply(plot_tbl[varying], function(col) length(unique(col)) > 1L, logical(1))]
    group_var <- if (length(cand) > 0) cand[1] else NULL
  } else {
    group_var <- simulation_resolve_design_variable(group_var, design_variable_aliases, "group_var", descriptor = design_descriptor)
    if (identical(group_var, x_var)) {
      stop("`group_var` must differ from `x_var`.")
    }
  }
  group_label <- if (is.null(group_var)) NULL else simulation_design_variable_label(group_var, design_variable_aliases)

  if (is.null(group_var)) {
    agg_tbl <- plot_tbl |>
      dplyr::group_by(.data[[x_var]]) |>
      dplyr::summarize(y = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data[[x_var]]) |>
      dplyr::mutate(group = "All designs")
  } else {
    agg_tbl <- plot_tbl |>
      dplyr::group_by(.data[[x_var]], .data[[group_var]]) |>
      dplyr::summarize(y = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data[[x_var]], .data[[group_var]]) |>
      dplyr::rename(group = dplyr::all_of(group_var))
  }

  out <- list(
    plot = "design_evaluation",
    facet = facet,
    metric = metric,
    metric_col = metric_col,
    x_var = x_var,
    x_label = x_label,
    group_var = group_var,
    group_label = group_label,
    design_variable_aliases = design_variable_aliases,
    design_descriptor = design_descriptor,
    planning_scope = simulation_object_planning_scope(x),
    planning_constraints = simulation_object_planning_constraints(x),
    planning_schema = simulation_object_planning_schema(x),
    data = agg_tbl
  )

  if (!isTRUE(draw)) return(out)

  groups <- unique(as.character(agg_tbl$group))
  cols <- grDevices::hcl.colors(max(1L, length(groups)), "Set 2")
  x_vals <- sort(unique(agg_tbl[[x_var]]))
  y_range <- range(agg_tbl$y, na.rm = TRUE)
  if (!all(is.finite(y_range))) stop("Selected metric has no finite values to plot.")

  graphics::plot(
    x = x_vals,
    y = rep(NA_real_, length(x_vals)),
    type = "n",
    xlab = x_label,
    ylab = metric_col,
    main = paste("Design simulation:", facet, metric_col),
    ylim = y_range
  )
  for (i in seq_along(groups)) {
    sub <- agg_tbl[as.character(agg_tbl$group) == groups[i], , drop = FALSE]
    sub <- sub[order(sub[[x_var]]), , drop = FALSE]
    graphics::lines(sub[[x_var]], sub$y, type = "b", lwd = 2, pch = 16 + (i - 1L) %% 5L, col = cols[i])
  }
  if (length(groups) > 1L) {
    graphics::legend("topleft", legend = groups, col = cols, lty = 1, lwd = 2, pch = 16 + (seq_along(groups) - 1L) %% 5L, bty = "n")
  }

  invisible(out)
}

#' Recommend a design condition from simulation results
#'
#' @param x Output from [evaluate_mfrm_design()] or [summary.mfrm_design_evaluation()].
#' @param facets Facets that must satisfy the planning thresholds.
#' @param min_separation Minimum acceptable mean separation.
#' @param min_reliability Minimum acceptable mean reliability.
#' @param max_severity_rmse Maximum acceptable severity recovery RMSE.
#' @param max_misfit_rate Maximum acceptable mean misfit rate.
#' @param min_convergence_rate Minimum acceptable convergence rate.
#' @param prefer Ranking priority among design variables. Earlier entries are
#'   optimized first when multiple designs pass. Custom public aliases from
#'   `sim_spec` are also accepted, as are the role keywords `person`, `rater`,
#'   `criterion`, and `assignment`.
#'
#' @details
#' This helper converts a design-study summary into a simple planning table.
#'
#' A design is marked as recommended when all requested facets satisfy all
#' selected thresholds simultaneously.
#' If multiple designs pass, the helper returns the smallest one according to
#' `prefer` (by default: fewer persons first, then fewer ratings per person,
#' then fewer raters, then fewer criteria).
#'
#' @section Typical workflow:
#' 1. Run [evaluate_mfrm_design()].
#' 2. Review [summary.mfrm_design_evaluation()] and [plot.mfrm_design_evaluation()].
#' 3. Use `recommend_mfrm_design(...)` to identify the smallest acceptable design.
#'
#' @return A list of class `mfrm_design_recommendation` with:
#' - `facet_table`: facet-level threshold checks, including design-variable
#'   alias columns when applicable
#' - `design_table`: design-level aggregated checks, including design-variable
#'   alias columns when applicable
#' - `recommended`: the first passing design after ranking
#' - `thresholds`: thresholds used in the recommendation
#' - `design_variable_aliases`: accepted public aliases for design variables
#' - `design_descriptor`: role-based design-variable metadata
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract
#' @seealso [evaluate_mfrm_design()], [summary.mfrm_design_evaluation], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- suppressWarnings(evaluate_mfrm_design(
#'   n_person = c(8, 12),
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 8,
#'   seed = 123
#' ))
#' rec <- recommend_mfrm_design(sim_eval)
#' rec$recommended
#' }
#' @export
recommend_mfrm_design <- function(x,
                                  facets = c("Rater", "Criterion"),
                                  min_separation = 2,
                                  min_reliability = 0.8,
                                  max_severity_rmse = 0.5,
                                  max_misfit_rate = 0.10,
                                  min_convergence_rate = 1,
                                  prefer = c("n_person", "raters_per_person", "n_rater", "n_criterion")) {
  if (inherits(x, "mfrm_design_evaluation")) {
    design_summary <- summary.mfrm_design_evaluation(x, digits = 6)$design_summary
  } else if (inherits(x, "summary.mfrm_design_evaluation")) {
    design_summary <- x$design_summary
  } else {
    stop("`x` must be output from evaluate_mfrm_design() or summary.mfrm_design_evaluation().")
  }

  design_summary <- tibble::as_tibble(design_summary)
  if (nrow(design_summary) == 0) stop("No design summary rows available.")

  if (missing(facets)) {
    facets <- setdiff(unique(as.character(design_summary$Facet)), "Person")
  } else {
    facets <- unique(as.character(facets))
  }
  missing_facets <- setdiff(facets, unique(design_summary$Facet))
  if (length(missing_facets) > 0) {
    stop("Requested facets not found in the design summary: ", paste(missing_facets, collapse = ", "))
  }

  design_variable_aliases <- simulation_object_design_variable_aliases(x)
  design_descriptor <- simulation_object_design_descriptor(x)
  prefer <- if (missing(prefer)) {
    c("n_person", "raters_per_person", "n_rater", "n_criterion")
  } else {
    simulation_resolve_design_variable_vector(prefer, design_variable_aliases, descriptor = design_descriptor)
  }
  if (length(prefer) == 0) {
    stop(
      "`prefer` must contain at least one valid design variable. Valid names: ",
      paste(simulation_design_variable_choices(design_variable_aliases, descriptor = design_descriptor), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  facet_table <- design_summary |>
    dplyr::filter(.data$Facet %in% facets) |>
    dplyr::mutate(
      SeparationPass = is.finite(.data$MeanSeparation) & .data$MeanSeparation >= min_separation,
      ReliabilityPass = is.finite(.data$MeanReliability) & .data$MeanReliability >= min_reliability,
      SeverityRMSEPass = is.finite(.data$MeanSeverityRMSE) & .data$MeanSeverityRMSE <= max_severity_rmse,
      MisfitRatePass = is.finite(.data$MeanMisfitRate) & .data$MeanMisfitRate <= max_misfit_rate,
      ConvergencePass = is.finite(.data$ConvergenceRate) & .data$ConvergenceRate >= min_convergence_rate,
      Pass = .data$SeparationPass & .data$ReliabilityPass & .data$SeverityRMSEPass &
        .data$MisfitRatePass & .data$ConvergencePass
    ) |>
    dplyr::arrange(.data$Facet, !!!rlang::syms(simulation_design_canonical_variables(design_descriptor)))

  design_table <- facet_table |>
    dplyr::group_by(dplyr::across(dplyr::all_of(simulation_design_group_variables(design_descriptor)))) |>
    dplyr::summarize(
      FacetsChecked = paste(.data$Facet, collapse = ", "),
      MinSeparation = min(.data$MeanSeparation, na.rm = TRUE),
      MinReliability = min(.data$MeanReliability, na.rm = TRUE),
      MaxSeverityRMSE = max(.data$MeanSeverityRMSE, na.rm = TRUE),
      MaxMisfitRate = max(.data$MeanMisfitRate, na.rm = TRUE),
      MinConvergenceRate = min(.data$ConvergenceRate, na.rm = TRUE),
      FacetsPassing = sum(.data$Pass, na.rm = TRUE),
      FacetsRequired = dplyr::n(),
      Pass = all(.data$Pass),
      .groups = "drop"
    )
  design_table <- simulation_append_design_alias_columns(design_table, design_variable_aliases)

  rank_vars <- c("Pass", prefer)
  design_table <- design_table |>
    dplyr::arrange(dplyr::desc(.data$Pass), !!!rlang::syms(prefer))

  recommended <- design_table |>
    dplyr::filter(.data$Pass) |>
    dplyr::slice_head(n = 1)

  structure(
    list(
      facet_table = facet_table,
      design_table = design_table,
      recommended = recommended,
      thresholds = list(
        facets = facets,
        min_separation = min_separation,
        min_reliability = min_reliability,
        max_severity_rmse = max_severity_rmse,
        max_misfit_rate = max_misfit_rate,
        min_convergence_rate = min_convergence_rate,
        prefer = prefer
      ),
      design_variable_aliases = design_variable_aliases,
      design_descriptor = design_descriptor,
      planning_scope = simulation_object_planning_scope(x),
      planning_constraints = simulation_object_planning_constraints(x),
      planning_schema = simulation_object_planning_schema(x)
    ),
    class = "mfrm_design_recommendation"
  )
}

diagnostic_screening_resolve_dependence_facet <- function(local_dependence_facet, facet_names) {
  facet_names <- as.character(facet_names %||% character(0))
  if (length(facet_names) < 2L) {
    stop("Local-dependence screening requires at least two non-person facets.", call. = FALSE)
  }

  if (is.null(local_dependence_facet)) {
    return(facet_names[2])
  }

  key <- as.character(local_dependence_facet[1] %||% "")
  role_lookup <- c(rater = facet_names[1], criterion = facet_names[2])

  if (key %in% facet_names) {
    return(key)
  }
  if (tolower(key) %in% names(role_lookup)) {
    return(unname(role_lookup[[tolower(key)]]))
  }

  stop(
    "`local_dependence_facet` must be one of `rater`, `criterion`, or the active facet names: ",
    paste(facet_names, collapse = ", "),
    ".",
    call. = FALSE
  )
}

diagnostic_screening_build_local_dependence_effects <- function(dat,
                                                                dependence_facet,
                                                                local_dependence_sd,
                                                                seed = NULL) {
  if (!is.data.frame(dat) || nrow(dat) == 0L) {
    return(tibble::tibble(Person = character(0), Effect = numeric(0)))
  }
  if (!"Person" %in% names(dat) || !dependence_facet %in% names(dat)) {
    stop("Local-dependence effect builder could not find `Person` and the requested dependence facet in the simulated data.", call. = FALSE)
  }
  if (!is.finite(local_dependence_sd) || local_dependence_sd < 0) {
    stop("`local_dependence_sd` must be a single non-negative numeric value.", call. = FALSE)
  }

  cells <- unique(dat[, c("Person", dependence_facet), drop = FALSE])
  if (identical(local_dependence_sd, 0)) {
    cells$Effect <- 0
    return(tibble::as_tibble(cells))
  }

  effects <- with_preserved_rng_seed(
    seed,
    stats::rnorm(nrow(cells), mean = 0, sd = local_dependence_sd)
  )
  effects <- effects - mean(effects, na.rm = TRUE)
  effects <- effects - stats::ave(
    effects,
    as.character(cells[[dependence_facet]]),
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  cells$Effect <- as.numeric(effects)
  tibble::as_tibble(cells)
}

diagnostic_screening_rescale_support <- function(values, target_sd) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (length(values) < 2L) {
    stop("Empirical latent support generation requires at least two finite values.", call. = FALSE)
  }

  values <- values - mean(values, na.rm = TRUE)
  target_sd <- as.numeric(target_sd[1] %||% 0)
  if (!is.finite(target_sd) || target_sd <= 0) {
    return(rep(0, length(values)))
  }

  current_sd <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(current_sd) || identical(current_sd, 0)) {
    return(rep(0, length(values)))
  }

  values / current_sd * target_sd
}

diagnostic_screening_build_latent_empirical_support <- function(theta_sd,
                                                                rater_sd,
                                                                criterion_sd,
                                                                support_n = 4000L,
                                                                person_shift = 1.5,
                                                                seed = NULL) {
  support_n <- as.integer(support_n[1] %||% 4000L)
  if (!is.finite(support_n) || support_n < 8L) {
    stop("`support_n` must be >= 8 for latent-misspecification screening.", call. = FALSE)
  }
  person_shift <- as.numeric(person_shift[1] %||% 1.5)
  if (!is.finite(person_shift) || person_shift <= 0) {
    stop("`person_shift` must be a positive numeric value.", call. = FALSE)
  }

  with_preserved_rng_seed(
    seed,
    {
      left_n <- support_n %/% 2L
      right_n <- support_n - left_n
      person_draws <- c(
        stats::rnorm(left_n, mean = -person_shift, sd = 1),
        stats::rnorm(right_n, mean = person_shift, sd = 1)
      )
      rater_draws <- stats::rnorm(support_n, mean = 0, sd = 1)
      criterion_draws <- stats::rnorm(support_n, mean = 0, sd = 1)

      list(
        person = diagnostic_screening_rescale_support(person_draws, theta_sd),
        rater = diagnostic_screening_rescale_support(rater_draws, rater_sd),
        criterion = diagnostic_screening_rescale_support(criterion_draws, criterion_sd)
      )
    }
  )
}

diagnostic_screening_build_latent_misspecification_spec <- function(row_spec_base,
                                                                    design_row,
                                                                    score_levels,
                                                                    theta_sd,
                                                                    rater_sd,
                                                                    criterion_sd,
                                                                    noise_sd,
                                                                    step_span,
                                                                    model,
                                                                    fit_step_facet,
                                                                    support_n = 4000L,
                                                                    person_shift = 1.5,
                                                                    seed = NULL) {
  empirical_support <- diagnostic_screening_build_latent_empirical_support(
    theta_sd = if (is.null(row_spec_base)) theta_sd else row_spec_base$theta_sd,
    rater_sd = if (is.null(row_spec_base)) rater_sd else row_spec_base$rater_sd,
    criterion_sd = if (is.null(row_spec_base)) criterion_sd else row_spec_base$criterion_sd,
    support_n = support_n,
    person_shift = person_shift,
    seed = seed
  )

  if (is.null(row_spec_base)) {
    spec_args <- list(
      n_person = design_row$n_person,
      n_rater = design_row$n_rater,
      n_criterion = design_row$n_criterion,
      raters_per_person = design_row$raters_per_person,
      score_levels = score_levels,
      theta_sd = theta_sd,
      rater_sd = rater_sd,
      criterion_sd = criterion_sd,
      noise_sd = noise_sd,
      step_span = step_span,
      model = model,
      assignment = if (identical(design_row$raters_per_person, design_row$n_rater)) "crossed" else "rotating",
      latent_distribution = "empirical",
      empirical_person = empirical_support$person,
      empirical_rater = empirical_support$rater,
      empirical_criterion = empirical_support$criterion
    )
    if (identical(model, "PCM")) {
      spec_args$step_facet <- fit_step_facet
    }
    return(do.call(build_mfrm_sim_spec, spec_args))
  }

  if (isTRUE((row_spec_base$population %||% list(active = FALSE))$active)) {
    stop(
      "Latent-misspecification screening does not yet support `sim_spec` objects with an active latent-regression population generator.",
      call. = FALSE
    )
  }

  out <- row_spec_base
  out$latent_distribution <- "empirical"
  out$empirical_support <- empirical_support
  out$planning_constraints <- simulation_planning_constraints(out)
  out$planning_schema <- simulation_planning_schema(out)
  out
}

diagnostic_screening_resolve_misspecified_step_facet <- function(fit_model,
                                                                 fit_step_facet,
                                                                 facet_names) {
  facet_names <- as.character(facet_names %||% character(0))
  if (length(facet_names) < 2L) {
    stop("Step-structure screening requires at least two non-person facets.", call. = FALSE)
  }
  if (identical(fit_model, "RSM")) {
    return(facet_names[2])
  }

  fit_step_facet <- as.character(fit_step_facet[1] %||% facet_names[2])
  alternate <- setdiff(facet_names, fit_step_facet)
  if (length(alternate) < 1L) {
    stop("Step-structure screening could not resolve an alternate step facet.", call. = FALSE)
  }
  alternate[1]
}

diagnostic_screening_build_step_threshold_table <- function(step_levels,
                                                            score_levels,
                                                            step_span) {
  step_levels <- as.character(step_levels)
  step_levels <- step_levels[!is.na(step_levels) & nzchar(step_levels)]
  if (length(step_levels) < 1L) {
    stop("Step-structure screening requires at least one generator step level.", call. = FALSE)
  }

  score_levels <- as.integer(score_levels[1] %||% 4L)
  if (!is.finite(score_levels) || score_levels < 2L) {
    stop("`score_levels` must be >= 2 for step-structure screening.", call. = FALSE)
  }

  base_steps <- if (score_levels == 2L) {
    0
  } else {
    seq(-abs(step_span), abs(step_span), length.out = score_levels - 1L)
  }
  offset_grid <- seq(-0.8, 0.8, length.out = length(step_levels))
  scale_grid <- seq(0.7, 1.35, length.out = length(step_levels))

  rows <- lapply(seq_along(step_levels), function(i) {
    tibble::tibble(
      StepFacet = step_levels[i],
      StepIndex = seq_len(score_levels - 1L),
      Step = paste0("Step_", seq_len(score_levels - 1L)),
      Estimate = base_steps * scale_grid[i] + offset_grid[i]
    )
  })
  dplyr::bind_rows(rows)
}

diagnostic_screening_build_step_misspecification_spec <- function(row_spec_base,
                                                                  design_row,
                                                                  score_levels,
                                                                  theta_sd,
                                                                  rater_sd,
                                                                  criterion_sd,
                                                                  noise_sd,
                                                                  step_span,
                                                                  fit_model,
                                                                  fit_step_facet) {
  if (!identical(fit_model, "RSM") && !identical(fit_model, "PCM")) {
    stop("Step-structure screening is currently scoped to `RSM` and `PCM` fits.", call. = FALSE)
  }
  if (!is.null(row_spec_base) && isTRUE((row_spec_base$population %||% list(active = FALSE))$active)) {
    stop(
      "Step-structure screening does not yet support `sim_spec` objects with an active latent-regression population generator.",
      call. = FALSE
    )
  }

  facet_names <- if (is.null(row_spec_base)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(row_spec_base)
  generator_step_facet <- diagnostic_screening_resolve_misspecified_step_facet(
    fit_model = fit_model,
    fit_step_facet = fit_step_facet,
    facet_names = facet_names
  )
  generator_step_role <- if (identical(generator_step_facet, facet_names[1])) "rater" else "criterion"
  step_levels <- if (is.null(row_spec_base)) {
    if (identical(generator_step_role, "rater")) {
      simulation_generated_role_levels(facet_names[1], "R", design_row$n_rater)
    } else {
      simulation_generated_role_levels(facet_names[2], "C", design_row$n_criterion)
    }
  } else {
    simulation_spec_role_levels(
      row_spec_base,
      role = generator_step_role,
      count = if (identical(generator_step_role, "rater")) design_row$n_rater else design_row$n_criterion
    )
  }
  threshold_table <- diagnostic_screening_build_step_threshold_table(
    step_levels = step_levels,
    score_levels = if (is.null(row_spec_base)) score_levels else row_spec_base$score_levels,
    step_span = if (is.null(row_spec_base)) step_span else row_spec_base$step_span
  )

  if (is.null(row_spec_base)) {
    return(build_mfrm_sim_spec(
      n_person = design_row$n_person,
      n_rater = design_row$n_rater,
      n_criterion = design_row$n_criterion,
      raters_per_person = design_row$raters_per_person,
      score_levels = score_levels,
      theta_sd = theta_sd,
      rater_sd = rater_sd,
      criterion_sd = criterion_sd,
      noise_sd = noise_sd,
      step_span = step_span,
      thresholds = threshold_table,
      model = "PCM",
      step_facet = generator_step_facet,
      assignment = if (identical(design_row$raters_per_person, design_row$n_rater)) "crossed" else "rotating",
      facet_names = facet_names
    ))
  }

  out <- row_spec_base
  out$model <- "PCM"
  out$step_facet <- generator_step_facet
  out$threshold_table <- tibble::as_tibble(threshold_table)
  out$planning_constraints <- simulation_planning_constraints(out)
  out$planning_schema <- simulation_planning_schema(out)
  out
}

diagnostic_screening_collect_metrics <- function(diag) {
  fit_tbl <- tibble::as_tibble(diag$fit %||% tibble::tibble())
  abs_z <- if (nrow(fit_tbl) > 0 && all(c("InfitZSTD", "OutfitZSTD") %in% names(fit_tbl))) {
    pmax(
      abs(suppressWarnings(as.numeric(fit_tbl$InfitZSTD))),
      abs(suppressWarnings(as.numeric(fit_tbl$OutfitZSTD))),
      na.rm = TRUE
    )
  } else {
    numeric(0)
  }
  abs_z <- abs_z[is.finite(abs_z)]

  marginal_summary <- tibble::as_tibble(diag$marginal_fit$summary %||% tibble::tibble())
  marginal_available <- isTRUE(diag$marginal_fit$available) && nrow(marginal_summary) > 0

  list(
    LegacyMeanAbsZ = if (length(abs_z) == 0L) NA_real_ else mean(abs_z, na.rm = TRUE),
    LegacyFlaggedLevels = if (length(abs_z) == 0L) NA_real_ else sum(abs_z >= 2, na.rm = TRUE),
    MarginalAvailable = marginal_available,
    MarginalOverallRMSD = if (marginal_available) suppressWarnings(as.numeric(marginal_summary$OverallRMSD[1])) else NA_real_,
    MarginalMaxAbsStdResidual = if (marginal_available) suppressWarnings(as.numeric(marginal_summary$OverallMaxAbsStdResidual[1])) else NA_real_,
    MarginalFlaggedGroups = if (marginal_available) {
      sum(
        suppressWarnings(as.numeric(marginal_summary$StepGroupsFlagged[1] %||% 0)),
        suppressWarnings(as.numeric(marginal_summary$FacetLevelsFlagged[1] %||% 0)),
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    PairwiseAvailable = isTRUE(diag$marginal_fit$pairwise$available),
    PairwiseFlaggedLevelPairs = if (marginal_available) suppressWarnings(as.numeric(marginal_summary$PairwiseFlaggedLevelPairs[1])) else NA_real_
  )
}

diagnostic_screening_safe_agreement <- function(a, b) {
  a <- as.logical(a)
  b <- as.logical(b)
  ok <- !is.na(a) & !is.na(b)
  if (!any(ok)) {
    return(NA_real_)
  }
  mean(a[ok] == b[ok], na.rm = TRUE)
}

diagnostic_screening_build_performance_summary <- function(results, design_variable_aliases = NULL) {
  results_tbl <- tibble::as_tibble(results)
  design_descriptor <- simulation_object_design_descriptor(
    list(results = results_tbl, design_variable_aliases = design_variable_aliases)
  )
  design_vars <- simulation_design_canonical_variables(design_descriptor)
  grouping_vars <- c("design_id", "Scenario", "ScenarioClass", "Model", "DependenceFacet", design_vars)

  flags_tbl <- results_tbl |>
    dplyr::mutate(
      LegacyAnyFlag = suppressWarnings(as.numeric(.data$LegacyFlaggedLevels)) > 0,
      MarginalAnyFlag = suppressWarnings(as.numeric(.data$MarginalFlaggedGroups)) > 0,
      PairwiseAnyFlag = suppressWarnings(as.numeric(.data$PairwiseFlaggedLevelPairs)) > 0,
      StrictAnyFlag = .data$MarginalAnyFlag | .data$PairwiseAnyFlag,
      ElapsedSecPer100Obs = dplyr::if_else(
        is.finite(.data$Observations) & .data$Observations > 0,
        .data$ElapsedSec / .data$Observations * 100,
        NA_real_
      )
    )

  out <- flags_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) |>
    dplyr::summarize(
      Reps = dplyr::n(),
      MeanElapsedSec = design_eval_safe_mean(.data$ElapsedSec),
      McseElapsedSec = simulation_mcse_mean(.data$ElapsedSec),
      MeanElapsedSecPer100Obs = design_eval_safe_mean(.data$ElapsedSecPer100Obs),
      LegacyAnyFlagRate = mean(.data$LegacyAnyFlag, na.rm = TRUE),
      McseLegacyAnyFlagRate = simulation_mcse_proportion(.data$LegacyAnyFlag),
      MarginalAnyFlagRate = mean(.data$MarginalAnyFlag, na.rm = TRUE),
      McseMarginalAnyFlagRate = simulation_mcse_proportion(.data$MarginalAnyFlag),
      PairwiseAnyFlagRate = mean(.data$PairwiseAnyFlag, na.rm = TRUE),
      McsePairwiseAnyFlagRate = simulation_mcse_proportion(.data$PairwiseAnyFlag),
      StrictAnyFlagRate = mean(.data$StrictAnyFlag, na.rm = TRUE),
      McseStrictAnyFlagRate = simulation_mcse_proportion(.data$StrictAnyFlag),
      LegacyVsMarginalAgreement = diagnostic_screening_safe_agreement(.data$LegacyAnyFlag, .data$MarginalAnyFlag),
      LegacyVsPairwiseAgreement = diagnostic_screening_safe_agreement(.data$LegacyAnyFlag, .data$PairwiseAnyFlag),
      LegacyVsStrictAgreement = diagnostic_screening_safe_agreement(.data$LegacyAnyFlag, .data$StrictAnyFlag),
      MarginalVsPairwiseAgreement = diagnostic_screening_safe_agreement(.data$MarginalAnyFlag, .data$PairwiseAnyFlag),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      EvaluationUse = dplyr::if_else(.data$Scenario == "well_specified", "type_I_proxy", "sensitivity_proxy"),
      LegacyTypeIProxy = dplyr::if_else(.data$Scenario == "well_specified", .data$LegacyAnyFlagRate, NA_real_),
      StrictTypeIProxy = dplyr::if_else(.data$Scenario == "well_specified", .data$StrictAnyFlagRate, NA_real_),
      LegacySensitivityProxy = dplyr::if_else(.data$Scenario != "well_specified", .data$LegacyAnyFlagRate, NA_real_),
      StrictSensitivityProxy = dplyr::if_else(.data$Scenario != "well_specified", .data$StrictAnyFlagRate, NA_real_),
      DeltaStrictMinusLegacyFlagRate = .data$StrictAnyFlagRate - .data$LegacyAnyFlagRate
    )
  simulation_append_design_alias_columns(out, design_variable_aliases)
}

diagnostic_screening_summarize_results <- function(results, design_variable_aliases = NULL) {
  results_tbl <- tibble::as_tibble(results)
  design_descriptor <- simulation_object_design_descriptor(
    list(results = results_tbl, design_variable_aliases = design_variable_aliases)
  )
  design_vars <- simulation_design_canonical_variables(design_descriptor)
  grouping_vars <- c("design_id", "Scenario", "ScenarioClass", "Model", "DependenceFacet", design_vars)

  scenario_summary <- results_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) |>
    dplyr::summarize(
      Reps = dplyr::n(),
      RunOKRate = mean(.data$RunOK, na.rm = TRUE),
      ConvergenceRate = mean(.data$Converged, na.rm = TRUE),
      MeanElapsedSec = design_eval_safe_mean(.data$ElapsedSec),
      MeanLegacyMeanAbsZ = design_eval_safe_mean(.data$LegacyMeanAbsZ),
      MeanLegacyFlaggedLevels = design_eval_safe_mean(.data$LegacyFlaggedLevels),
      LegacyAnyFlagRate = mean(.data$LegacyFlaggedLevels > 0, na.rm = TRUE),
      MeanMarginalOverallRMSD = design_eval_safe_mean(.data$MarginalOverallRMSD),
      MeanMarginalMaxAbsStdResidual = design_eval_safe_mean(.data$MarginalMaxAbsStdResidual),
      MeanMarginalFlaggedGroups = design_eval_safe_mean(.data$MarginalFlaggedGroups),
      MarginalAnyFlagRate = mean(.data$MarginalFlaggedGroups > 0, na.rm = TRUE),
      MeanPairwiseFlaggedLevelPairs = design_eval_safe_mean(.data$PairwiseFlaggedLevelPairs),
      PairwiseAnyFlagRate = mean(.data$PairwiseFlaggedLevelPairs > 0, na.rm = TRUE),
      PairwiseAvailabilityRate = mean(as.numeric(.data$PairwiseAvailable), na.rm = TRUE),
      .groups = "drop"
    )
  scenario_summary <- simulation_append_design_alias_columns(scenario_summary, design_variable_aliases)
  performance_summary <- diagnostic_screening_build_performance_summary(
    results_tbl,
    design_variable_aliases = design_variable_aliases
  )

  scenario_contrast <- tibble::tibble()
  if ("well_specified" %in% scenario_summary$Scenario &&
      any(scenario_summary$Scenario != "well_specified")) {
    by_vars <- c("design_id", "Model", "DependenceFacet", design_vars)
    contrast_cols <- c(
      "MeanLegacyMeanAbsZ",
      "MeanLegacyFlaggedLevels",
      "MeanMarginalOverallRMSD",
      "MeanMarginalMaxAbsStdResidual",
      "MeanMarginalFlaggedGroups",
      "MeanPairwiseFlaggedLevelPairs",
      "PairwiseAnyFlagRate"
    )
    well_tbl <- scenario_summary |>
      dplyr::filter(.data$Scenario == "well_specified") |>
      dplyr::select(dplyr::all_of(c(by_vars, contrast_cols)))
    alt_tbl <- scenario_summary |>
      dplyr::filter(.data$Scenario != "well_specified") |>
      dplyr::select(dplyr::all_of(c("Scenario", "ScenarioClass", by_vars, contrast_cols)))

    scenario_contrast <- dplyr::inner_join(
      alt_tbl,
      well_tbl,
      by = by_vars,
      suffix = c("_Scenario", "_WellSpecified")
    ) |>
      dplyr::mutate(
        DeltaLegacyMeanAbsZ = .data$MeanLegacyMeanAbsZ_Scenario - .data$MeanLegacyMeanAbsZ_WellSpecified,
        DeltaLegacyFlaggedLevels = .data$MeanLegacyFlaggedLevels_Scenario - .data$MeanLegacyFlaggedLevels_WellSpecified,
        DeltaMarginalOverallRMSD = .data$MeanMarginalOverallRMSD_Scenario - .data$MeanMarginalOverallRMSD_WellSpecified,
        DeltaMarginalMaxAbsStdResidual = .data$MeanMarginalMaxAbsStdResidual_Scenario - .data$MeanMarginalMaxAbsStdResidual_WellSpecified,
        DeltaMarginalFlaggedGroups = .data$MeanMarginalFlaggedGroups_Scenario - .data$MeanMarginalFlaggedGroups_WellSpecified,
        DeltaPairwiseFlaggedLevelPairs = .data$MeanPairwiseFlaggedLevelPairs_Scenario - .data$MeanPairwiseFlaggedLevelPairs_WellSpecified,
        DeltaPairwiseAnyFlagRate = .data$PairwiseAnyFlagRate_Scenario - .data$PairwiseAnyFlagRate_WellSpecified,
        StrictSignalImproved = (.data$DeltaMarginalFlaggedGroups > 0) | (.data$DeltaPairwiseFlaggedLevelPairs > 0),
        StrictSignalDominatesLegacy = (.data$DeltaMarginalFlaggedGroups > .data$DeltaLegacyFlaggedLevels) |
          (.data$DeltaPairwiseFlaggedLevelPairs > .data$DeltaLegacyFlaggedLevels)
      )
    scenario_contrast <- simulation_append_design_alias_columns(scenario_contrast, design_variable_aliases)
  }

  notes <- character(0)
  if (nrow(performance_summary) > 0 &&
      any(performance_summary$Scenario == "well_specified", na.rm = TRUE)) {
    notes <- c(
      notes,
      "Well-specified rows report screening-oriented Type I proxies from any-flag rates; they are not calibrated inferential alpha estimates."
    )
  }
  if (nrow(performance_summary) > 0 &&
      any(performance_summary$Scenario != "well_specified", na.rm = TRUE)) {
    notes <- c(
      notes,
      "Misspecification rows report screening-oriented sensitivity proxies from any-flag rates; they summarize detection behavior rather than formal power."
    )
  }
  if (nrow(scenario_contrast) > 0 && any(scenario_contrast$StrictSignalImproved, na.rm = TRUE)) {
    notes <- c(notes, "At least one misspecification scenario increased a strict screening signal relative to the well-specified baseline for an evaluated design row.")
  }
  if (nrow(scenario_contrast) > 0 && any(scenario_contrast$DeltaPairwiseFlaggedLevelPairs > 0, na.rm = TRUE)) {
    notes <- c(notes, "At least one misspecification scenario increased strict pairwise flagging relative to the well-specified baseline for an evaluated design row.")
  }
  if (nrow(scenario_contrast) > 0 && any(scenario_contrast$StrictSignalDominatesLegacy, na.rm = TRUE)) {
    notes <- c(notes, "At least one strict screening signal reacted more strongly than the legacy |ZSTD| screen for an evaluated design row.")
  }
  if (length(notes) == 0L) {
    notes <- "Compare `scenario_summary` and `scenario_contrast` to judge whether strict marginal diagnostics separate the misfit scenarios more clearly than legacy residual screens."
  }

  list(
    scenario_summary = scenario_summary,
    performance_summary = performance_summary,
    scenario_contrast = scenario_contrast,
    notes = unique(notes)
  )
}

#' Evaluate legacy and strict marginal diagnostic screening under controlled misfit scenarios
#'
#' @param n_person Vector of person counts to evaluate.
#' @param n_rater Vector of rater counts to evaluate.
#' @param n_criterion Vector of criterion counts to evaluate.
#' @param raters_per_person Vector of rater assignments per person.
#' @param design Optional named design-grid override supplied as a named list,
#'   named vector, or one-row data frame. Names may use canonical variables
#'   (`n_person`, `n_rater`, `n_criterion`, `raters_per_person`), current
#'   public aliases implied by `sim_spec`, or role keywords
#'   (`person`, `rater`, `criterion`, `assignment`). Values may be vectors.
#' @param reps Number of replications per design condition and scenario.
#' @param scenarios Screening scenarios to evaluate. The current first release
#'   supports `"well_specified"`, `"local_dependence"`, and
#'   `"latent_misspecification"`, plus
#'   `"step_structure_misspecification"`.
#' @param local_dependence_sd Standard deviation of the shared context effect
#'   injected in the `"local_dependence"` scenario.
#' @param local_dependence_facet Facet that receives the shared
#'   `Person x facet` dependence effect. Use `"criterion"`, `"rater"`, or an
#'   active public facet name. Defaults to the criterion-like facet.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread of step thresholds on the logit scale.
#' @param model Measurement model passed to [fit_mfrm()]. The current helper
#'   supports `RSM` and `PCM`; bounded `GPCM` is accepted only to
#'   produce an explicit unsupported-path error.
#' @param step_facet Step facet passed to [fit_mfrm()] when `model = "PCM"`.
#' @param maxit Maximum iterations passed to [fit_mfrm()].
#' @param quad_points Quadrature points for the internal `MML` fit.
#' @param residual_pca Residual PCA mode passed to [diagnose_mfrm()].
#' @param sim_spec Optional output from [build_mfrm_sim_spec()] or
#'   [extract_mfrm_sim_spec()] used as the base data-generating mechanism.
#' @param seed Optional seed for reproducible replications.
#'
#' @details
#' This helper performs a compact Monte Carlo validation study for the package's
#' current diagnostic architecture.
#'
#' For each design condition and scenario, the function:
#' 1. generates synthetic data with [simulate_mfrm_data()]
#' 2. fits the model with `method = "MML"`
#' 3. computes diagnostics with `diagnostic_mode = "both"`
#' 4. stores legacy residual-screen metrics and strict marginal-fit metrics
#' 5. aggregates the results into `scenario_summary` and `scenario_contrast`
#'
#' The `"well_specified"` scenario uses the ordinary generator with no injected
#' extra structure. The `"local_dependence"` scenario adds a shared
#' `Person x facet` random effect, centered within the selected facet levels, so
#' responses in the same context become correlated without changing the
#' facet-level mean effect contract. The `"latent_misspecification"` scenario
#' keeps the same marginal spread targets but replaces the normal person
#' distribution with a centered bimodal empirical support distribution, while
#' leaving the non-person facets on the original scale contract. The
#' `"step_structure_misspecification"` scenario uses a `PCM` generator with
#' facet-specific threshold tables that intentionally mismatch the fitted step
#' contract: `RSM` fits receive criterion-specific thresholds, and `PCM` fits
#' receive thresholds indexed by the opposite non-person facet.
#'
#' This function is intentionally screening-oriented. The strict marginal branch
#' remains exploratory in the current release, so the returned summaries should
#' be used to compare relative sensitivity across scenarios rather than to claim
#' calibrated inferential power.
#'
#' @return An object of class `mfrm_diagnostic_screening` with:
#' - `design_grid`: evaluated design conditions, including public alias columns
#'   when applicable
#' - `results`: replicate-level screening metrics for each design and scenario
#' - `scenario_summary`: aggregated scenario-by-design screening summaries
#' - `performance_summary`: scenario-by-design screening-performance summary
#'   including runtime, agreement, Type I proxy, and sensitivity proxy columns
#' - `scenario_contrast`: each misspecification scenario minus the
#'   well-specified baseline when the baseline scenario was evaluated
#' - `design_descriptor`: role-based design-variable metadata
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract
#' - `settings`: simulation and fitting settings
#' - `ademp`: simulation-study metadata
#' - `notes`: short interpretation notes
#'
#' @seealso [simulate_mfrm_data()], [evaluate_mfrm_design()], [diagnose_mfrm()]
#' @examples
#' \donttest{
#' diag_eval <- evaluate_mfrm_diagnostic_screening(
#'   design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
#'   reps = 1,
#'   maxit = 6,
#'   seed = 123
#' )
#' diag_eval$scenario_summary
#' diag_eval$scenario_contrast
#' }
#' @export
evaluate_mfrm_diagnostic_screening <- function(n_person = c(30, 50, 100),
                                               n_rater = c(4),
                                               n_criterion = c(4),
                                               raters_per_person = n_rater,
                                               design = NULL,
                                               reps = 10,
                                               scenarios = c("well_specified", "local_dependence"),
                                               local_dependence_sd = 0.8,
                                               local_dependence_facet = NULL,
                                               score_levels = 4,
                                               theta_sd = 1,
                                               rater_sd = 0.35,
                                               criterion_sd = 0.25,
                                               noise_sd = 0,
                                               step_span = 1.4,
                                               model = c("RSM", "PCM", "GPCM"),
                                               step_facet = NULL,
                                               maxit = 25,
                                               quad_points = 7,
                                               residual_pca = c("none", "overall", "facet", "both"),
                                               sim_spec = NULL,
                                               seed = NULL) {
  model <- match.arg(model)
  if (identical(model, "GPCM")) {
    stop(
      "`evaluate_mfrm_diagnostic_screening()` does not yet support bounded `GPCM`. ",
      "Strict marginal validation is currently scoped to `RSM` / `PCM`.",
      call. = FALSE
    )
  }
  residual_pca <- match.arg(residual_pca)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!is.null(sim_spec) && identical(as.character(sim_spec$model %||% NA_character_), "GPCM")) {
    stop(
      "`evaluate_mfrm_diagnostic_screening()` does not yet support bounded `GPCM` simulation specifications. ",
      "Strict marginal validation remains scoped to `RSM` / `PCM`.",
      call. = FALSE
    )
  }

  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) {
    stop("`reps` must be >= 1.", call. = FALSE)
  }

  scenarios <- unique(as.character(scenarios))
  valid_scenarios <- c(
    "well_specified",
    "local_dependence",
    "latent_misspecification",
    "step_structure_misspecification"
  )
  if (length(scenarios) == 0L || !all(scenarios %in% valid_scenarios)) {
    stop(
      "`scenarios` must be drawn from: ",
      paste(valid_scenarios, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  supplied_counts <- intersect(
    names(as.list(match.call(expand.dots = FALSE))[-1]),
    c("n_person", "n_rater", "n_criterion", "raters_per_person")
  )
  design_meta <- simulation_build_design_grid(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    sim_spec = sim_spec,
    id_prefix = "V",
    explicit_scalar_names = supplied_counts
  )
  design_grid <- design_meta$canonical
  design_grid_public <- design_meta$public
  design_variable_aliases <- design_meta$aliases
  design_descriptor <- design_meta$descriptor
  planning_scope <- simulation_planning_scope(sim_spec)
  planning_constraints <- simulation_planning_constraints(sim_spec)
  planning_schema <- simulation_planning_schema(sim_spec)

  generator_model <- if (is.null(sim_spec)) model else sim_spec$model
  base_facet_names <- if (is.null(sim_spec)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(sim_spec)
  generator_step_facet <- if (is.null(sim_spec)) {
    if (identical(generator_model, "PCM")) base_facet_names[2] else NA_character_
  } else {
    sim_spec$step_facet
  }
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
  fit_step_facet <- simulation_resolve_fit_step_facet(model, step_facet, generator_step_facet)

  seeds <- with_preserved_rng_seed(
    seed,
    sample.int(.Machine$integer.max, size = nrow(design_grid) * reps * length(scenarios), replace = FALSE)
  )
  seed_idx <- 0L
  out_idx <- 0L
  result_rows <- vector("list", nrow(design_grid) * reps * length(scenarios))

  for (i in seq_len(nrow(design_grid))) {
    design_row <- design_grid[i, , drop = FALSE]
    row_spec_base <- if (is.null(sim_spec)) {
      NULL
    } else {
      simulation_override_spec_design(
        sim_spec,
        n_person = design_row$n_person,
        n_rater = design_row$n_rater,
        n_criterion = design_row$n_criterion,
        raters_per_person = design_row$raters_per_person
      )
    }
    row_score_levels <- if (is.null(row_spec_base)) score_levels else row_spec_base$score_levels
    row_facet_names <- if (is.null(row_spec_base)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(row_spec_base)
    dependence_facet_name <- diagnostic_screening_resolve_dependence_facet(local_dependence_facet, row_facet_names)

    for (scenario in scenarios) {
      for (rep in seq_len(reps)) {
        seed_idx <- seed_idx + 1L
        sim_seed <- seeds[seed_idx]
        scenario_class <- dplyr::case_when(
          identical(scenario, "well_specified") ~ "null_reference",
          identical(scenario, "local_dependence") ~ "context_shared_person_by_facet_effect",
          identical(scenario, "latent_misspecification") ~ "bimodal_person_distribution_under_normal_fit",
          identical(scenario, "step_structure_misspecification") ~ "facet_specific_thresholds_under_misspecified_step_contract",
          TRUE ~ "other"
        )

        interaction_effects <- NULL
        row_spec <- row_spec_base
        if (identical(scenario, "local_dependence")) {
          skeleton <- if (is.null(row_spec_base)) {
            sim_args <- list(
              n_person = design_row$n_person,
              n_rater = design_row$n_rater,
              n_criterion = design_row$n_criterion,
              raters_per_person = design_row$raters_per_person,
              score_levels = score_levels,
              theta_sd = theta_sd,
              rater_sd = rater_sd,
              criterion_sd = criterion_sd,
              noise_sd = noise_sd,
              step_span = step_span,
              model = model,
              seed = sim_seed
            )
            if (identical(model, "PCM")) {
              sim_args$step_facet <- fit_step_facet
            }
            do.call(simulate_mfrm_data, sim_args)
          } else {
            simulate_mfrm_data(sim_spec = row_spec_base, seed = sim_seed)
          }
          interaction_effects <- diagnostic_screening_build_local_dependence_effects(
            dat = skeleton,
            dependence_facet = dependence_facet_name,
            local_dependence_sd = local_dependence_sd,
            seed = sim_seed + 1L
          )
          if (!is.null(row_spec_base)) {
            row_spec <- simulation_override_spec_design(
              row_spec_base,
              interaction_effects = interaction_effects
            )
          }
        } else if (identical(scenario, "latent_misspecification")) {
          row_spec <- diagnostic_screening_build_latent_misspecification_spec(
            row_spec_base = row_spec_base,
            design_row = design_row,
            score_levels = score_levels,
            theta_sd = theta_sd,
            rater_sd = rater_sd,
            criterion_sd = criterion_sd,
            noise_sd = noise_sd,
            step_span = step_span,
            model = model,
            fit_step_facet = fit_step_facet,
            seed = sim_seed + 1L
          )
        } else if (identical(scenario, "step_structure_misspecification")) {
          row_spec <- diagnostic_screening_build_step_misspecification_spec(
            row_spec_base = row_spec_base,
            design_row = design_row,
            score_levels = score_levels,
            theta_sd = theta_sd,
            rater_sd = rater_sd,
            criterion_sd = criterion_sd,
            noise_sd = noise_sd,
            step_span = step_span,
            fit_model = model,
            fit_step_facet = fit_step_facet
          )
        }

        sim <- if (is.null(row_spec)) {
          sim_args <- list(
            n_person = design_row$n_person,
            n_rater = design_row$n_rater,
            n_criterion = design_row$n_criterion,
            raters_per_person = design_row$raters_per_person,
            score_levels = score_levels,
            theta_sd = theta_sd,
            rater_sd = rater_sd,
            criterion_sd = criterion_sd,
            noise_sd = noise_sd,
            step_span = step_span,
            interaction_effects = interaction_effects,
            model = model,
            seed = sim_seed
          )
          if (identical(model, "PCM")) {
            sim_args$step_facet <- fit_step_facet
          }
          do.call(simulate_mfrm_data, sim_args)
        } else {
          simulate_mfrm_data(sim_spec = row_spec, seed = sim_seed)
        }

        t0 <- proc.time()[["elapsed"]]
        fit_args <- list(
          data = sim,
          person = "Person",
          facets = row_facet_names,
          score = "Score",
          method = "MML",
          model = model,
          maxit = maxit,
          quad_points = quad_points
        )
        if (identical(model, "PCM")) {
          fit_args$step_facet <- fit_step_facet
        }
        if ("Weight" %in% names(sim)) {
          fit_args$weight <- "Weight"
        }
        sim_population <- attr(sim, "mfrm_population_data")
        if (is.list(sim_population) && isTRUE(sim_population$active)) {
          fit_args$population_formula <- sim_population$population_formula
          fit_args$person_data <- sim_population$person_data
          fit_args$person_id <- sim_population$person_id
          fit_args$population_policy <- sim_population$population_policy
        }

        fit <- tryCatch(do.call(fit_mfrm, fit_args), error = function(e) e)
        diag <- if (inherits(fit, "error")) fit else {
          tryCatch(
            diagnose_mfrm(
              fit,
              diagnostic_mode = "both",
              residual_pca = residual_pca
            ),
            error = function(e) e
          )
        }
        elapsed <- proc.time()[["elapsed"]] - t0

        out_idx <- out_idx + 1L
        result_row <- tibble::tibble(
          design_id = design_row$design_id,
          Scenario = scenario,
          ScenarioClass = scenario_class,
          DependenceFacet = dependence_facet_name,
          Model = model,
          rep = rep,
          n_person = design_row$n_person,
          n_rater = design_row$n_rater,
          n_criterion = design_row$n_criterion,
          raters_per_person = design_row$raters_per_person,
          Observations = nrow(sim),
          MinCategoryCount = min(tabulate(sim$Score, nbins = row_score_levels)),
          ElapsedSec = elapsed,
          RunOK = FALSE,
          Converged = FALSE,
          Error = NA_character_,
          LegacyMeanAbsZ = NA_real_,
          LegacyFlaggedLevels = NA_real_,
          MarginalAvailable = FALSE,
          MarginalOverallRMSD = NA_real_,
          MarginalMaxAbsStdResidual = NA_real_,
          MarginalFlaggedGroups = NA_real_,
          PairwiseAvailable = FALSE,
          PairwiseFlaggedLevelPairs = NA_real_
        )

        if (inherits(fit, "error")) {
          result_row$Error <- conditionMessage(fit)
          result_rows[[out_idx]] <- result_row
          next
        }
        if (inherits(diag, "error")) {
          result_row$Error <- conditionMessage(diag)
          result_rows[[out_idx]] <- result_row
          next
        }

        metrics <- diagnostic_screening_collect_metrics(diag)
        result_row$RunOK <- TRUE
        result_row$Converged <- isTRUE(as.logical(fit$summary$Converged[1]))
        result_row$LegacyMeanAbsZ <- metrics$LegacyMeanAbsZ
        result_row$LegacyFlaggedLevels <- metrics$LegacyFlaggedLevels
        result_row$MarginalAvailable <- metrics$MarginalAvailable
        result_row$MarginalOverallRMSD <- metrics$MarginalOverallRMSD
        result_row$MarginalMaxAbsStdResidual <- metrics$MarginalMaxAbsStdResidual
        result_row$MarginalFlaggedGroups <- metrics$MarginalFlaggedGroups
        result_row$PairwiseAvailable <- metrics$PairwiseAvailable
        result_row$PairwiseFlaggedLevelPairs <- metrics$PairwiseFlaggedLevelPairs
        result_rows[[out_idx]] <- result_row
      }
    }
  }

  results <- dplyr::bind_rows(result_rows[seq_len(out_idx)])
  results <- simulation_append_design_alias_columns(results, design_variable_aliases)
  summary_bundle <- diagnostic_screening_summarize_results(
    results,
    design_variable_aliases = design_variable_aliases
  )
  ademp <- simulation_build_ademp(
    purpose = "Compare legacy residual screens and strict marginal-fit screens under well-specified and local-dependence simulation scenarios.",
    design_grid = design_grid,
    generator_model = generator_model,
    generator_step_facet = generator_step_facet,
    generator_assignment = generator_assignment,
    sim_spec = sim_spec,
    estimands = c(
      "Legacy residual-screen magnitude and flagged-level counts",
      "Strict marginal overall RMSD and flagged-group counts",
      "Strict pairwise local-dependence flagged-pair counts",
      "Scenario contrasts between local dependence and well-specified runs"
    ),
    analysis_methods = list(
      fit_method = "MML",
      fitted_model = model,
      maxit = maxit,
      quad_points = quad_points,
      diagnostic_mode = "both",
      residual_pca = residual_pca,
      scenarios = scenarios,
      local_dependence_sd = local_dependence_sd,
      local_dependence_facet = local_dependence_facet %||% "criterion"
    ),
    performance_measures = c(
      "Scenario-specific mean screening metrics",
      "Run and convergence rates",
      "Scenario contrasts in strict versus legacy screens"
    )
  )

  structure(
    list(
      design_grid = design_grid_public,
      results = results,
      scenario_summary = summary_bundle$scenario_summary,
      performance_summary = summary_bundle$performance_summary,
      scenario_contrast = summary_bundle$scenario_contrast,
      design_descriptor = design_descriptor,
      planning_scope = planning_scope,
      planning_constraints = planning_constraints,
      planning_schema = planning_schema,
      settings = list(
        reps = reps,
        scenarios = scenarios,
        local_dependence_sd = local_dependence_sd,
        local_dependence_facet = local_dependence_facet %||% "criterion",
        score_levels = score_levels,
        theta_sd = theta_sd,
        rater_sd = rater_sd,
        criterion_sd = criterion_sd,
        noise_sd = noise_sd,
        step_span = step_span,
        model = model,
        step_facet = fit_step_facet,
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        sim_spec = sim_spec,
        facet_names = stats::setNames(base_facet_names, c("rater", "criterion")),
        design_variable_aliases = design_variable_aliases,
        design_descriptor = design_descriptor,
        planning_scope = planning_scope,
        planning_constraints = planning_constraints,
        planning_schema = planning_schema,
        generator_model = generator_model,
        generator_step_facet = generator_step_facet,
        generator_assignment = generator_assignment,
        seed = seed
      ),
      ademp = ademp,
      notes = summary_bundle$notes
    ),
    class = "mfrm_diagnostic_screening"
  )
}

signal_eval_resolve_level <- function(value, prefix, n_levels, arg_name) {
  if (is.null(value)) {
    return(sprintf("%s%02d", prefix, n_levels))
  }
  if (is.numeric(value)) {
    idx <- as.integer(value[1])
    if (!is.finite(idx) || idx < 1L || idx > n_levels) {
      stop("`", arg_name, "` must be between 1 and ", n_levels, ".")
    }
    return(sprintf("%s%02d", prefix, idx))
  }
  value <- as.character(value[1])
  if (!nzchar(value) || is.na(value)) {
    stop("`", arg_name, "` must be a non-empty label or index.")
  }
  value
}

signal_eval_resolve_level_from_choices <- function(value, available_levels, arg_name) {
  available_levels <- as.character(available_levels)
  available_levels <- available_levels[!is.na(available_levels) & nzchar(available_levels)]
  if (length(available_levels) == 0L) {
    stop("No candidate levels are available for `", arg_name, "`.", call. = FALSE)
  }
  if (is.null(value)) {
    return(available_levels[length(available_levels)])
  }
  if (is.numeric(value)) {
    idx <- as.integer(value[1])
    if (!is.finite(idx) || idx < 1L || idx > length(available_levels)) {
      stop("`", arg_name, "` must be between 1 and ", length(available_levels), ".", call. = FALSE)
    }
    return(available_levels[idx])
  }
  value <- as.character(value[1])
  if (!nzchar(value) || is.na(value)) {
    stop("`", arg_name, "` must be a non-empty label or index.", call. = FALSE)
  }
  if (!value %in% available_levels) {
    stop(
      "`", arg_name, "` must match one of: ",
      paste(available_levels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  value
}

signal_eval_find_dif_row <- function(tbl, level, group1, group2) {
  tbl <- tibble::as_tibble(tbl)
  if (nrow(tbl) == 0) return(tbl)
  out <- tbl |>
    dplyr::filter(
      .data$Level == level,
      (.data$Group1 == group1 & .data$Group2 == group2) |
        (.data$Group1 == group2 & .data$Group2 == group1)
    ) |>
    dplyr::slice_head(n = 1)
  out
}

signal_eval_find_bias_row <- function(tbl, rater_level, criterion_level) {
  tbl <- tibble::as_tibble(tbl)
  if (nrow(tbl) == 0) return(tbl)
  out <- tbl |>
    dplyr::filter(
      .data$FacetA_Level == rater_level,
      .data$FacetB_Level == criterion_level
    ) |>
    dplyr::slice_head(n = 1)
  out
}

signal_eval_get_p <- function(tbl, adjusted_col = "p_adjusted", raw_col = "p_value") {
  tbl <- tibble::as_tibble(tbl)
  if (nrow(tbl) == 0) return(NA_real_)
  p_adj <- suppressWarnings(as.numeric(tbl[[adjusted_col]][1]))
  if (is.finite(p_adj)) return(p_adj)
  p_raw <- suppressWarnings(as.numeric(tbl[[raw_col]][1]))
  if (is.finite(p_raw)) return(p_raw)
  NA_real_
}

signal_eval_false_positive_rate <- function(flag_vec) {
  if (length(flag_vec) == 0) return(NA_real_)
  mean(flag_vec, na.rm = TRUE)
}

signal_eval_summary <- function(results, rep_overview, design_variable_aliases = NULL) {
  results_tbl <- tibble::as_tibble(results)
  rep_tbl <- tibble::as_tibble(rep_overview)
  design_descriptor <- simulation_object_design_descriptor(list(
    results = results_tbl,
    rep_overview = rep_tbl,
    design_variable_aliases = design_variable_aliases
  ))
  design_group_vars <- simulation_design_group_variables(design_descriptor)
  design_arrange_vars <- simulation_design_canonical_variables(design_descriptor)

  detection_summary <- tibble::tibble()
  if (nrow(results_tbl) > 0) {
    detection_summary <- results_tbl |>
      dplyr::group_by(
        dplyr::across(dplyr::all_of(design_group_vars)),
        .data$DIFTargetLevel,
        .data$BiasTargetRater,
        .data$BiasTargetCriterion
      ) |>
      dplyr::summarize(
        Reps = dplyr::n(),
        ConvergenceRate = mean(.data$Converged, na.rm = TRUE),
        McseConvergenceRate = simulation_mcse_proportion(.data$Converged),
        DIFPower = mean(.data$DIFDetected, na.rm = TRUE),
        McseDIFPower = simulation_mcse_proportion(.data$DIFDetected),
        DIFClassificationPower = mean(.data$DIFClassDetected, na.rm = TRUE),
        McseDIFClassificationPower = simulation_mcse_proportion(.data$DIFClassDetected),
        MeanTargetContrast = mean(.data$DIFContrast, na.rm = TRUE),
        McseTargetContrast = simulation_mcse_mean(.data$DIFContrast),
        MeanTargetContrastAbs = mean(abs(.data$DIFContrast), na.rm = TRUE),
        McseTargetContrastAbs = simulation_mcse_mean(abs(.data$DIFContrast)),
        DIFFalsePositiveRate = mean(.data$DIFFalsePositiveRate, na.rm = TRUE),
        McseDIFFalsePositiveRate = simulation_mcse_mean(.data$DIFFalsePositiveRate),
        BiasScreenRate = mean(.data$BiasDetected, na.rm = TRUE),
        McseBiasScreenRate = simulation_mcse_proportion(.data$BiasDetected),
        MeanTargetBias = mean(.data$BiasSize, na.rm = TRUE),
        McseTargetBias = simulation_mcse_mean(.data$BiasSize),
        MeanAbsTargetBias = mean(abs(.data$BiasSize), na.rm = TRUE),
        McseAbsTargetBias = simulation_mcse_mean(abs(.data$BiasSize)),
        MeanTargetBiasT = mean(.data$BiasT, na.rm = TRUE),
        McseTargetBiasT = simulation_mcse_mean(.data$BiasT),
        BiasScreenMetricAvailabilityRate = mean(.data$BiasScreenMetricAvailable, na.rm = TRUE),
        McseBiasScreenMetricAvailabilityRate = simulation_mcse_proportion(.data$BiasScreenMetricAvailable),
        BiasScreenFalsePositiveRate = mean(.data$BiasScreenFalsePositiveRate, na.rm = TRUE),
        McseBiasScreenFalsePositiveRate = simulation_mcse_mean(.data$BiasScreenFalsePositiveRate),
        MeanElapsedSec = mean(.data$ElapsedSec, na.rm = TRUE),
        McseElapsedSec = simulation_mcse_mean(.data$ElapsedSec),
        .groups = "drop"
      ) |>
      dplyr::arrange(!!!rlang::syms(design_arrange_vars))
    detection_summary <- simulation_append_design_alias_columns(detection_summary, design_variable_aliases)
  }

  overview <- tibble::tibble(
    Designs = dplyr::n_distinct(rep_tbl$design_id),
    Replications = nrow(rep_tbl),
    SuccessfulRuns = sum(rep_tbl$RunOK, na.rm = TRUE),
    ConvergedRuns = sum(rep_tbl$Converged, na.rm = TRUE),
    MeanElapsedSec = design_eval_safe_mean(rep_tbl$ElapsedSec)
  )

  notes <- character(0)
  if (nrow(detection_summary) > 0 && any(detection_summary$ConvergenceRate < 1, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions did not converge in every replication.")
  }
  if (nrow(detection_summary) > 0 && any(detection_summary$DIFPower < 0.8, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions showed DIF power below 0.80.")
  }
  if (nrow(detection_summary) > 0 && any(detection_summary$BiasScreenRate < 0.8, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions showed bias-screen hit rates below 0.80.")
  }
  if (nrow(detection_summary) > 0 && any(detection_summary$BiasScreenMetricAvailabilityRate < 1, na.rm = TRUE)) {
    notes <- c(notes, "Some design conditions did not yield usable bias-screening t/p metrics in every replication.")
  }
  notes <- c(
    notes,
    "Bias-side rates are screening summaries derived from `estimate_bias()` output and should not be interpreted as formal power or alpha-calibrated false-positive rates.",
    "MCSE columns summarize finite-replication uncertainty around the reported means and rates."
  )

  list(
    overview = overview,
    detection_summary = detection_summary,
    notes = notes
  )
}

signal_eval_metric_col <- function(signal, metric) {
  if (identical(signal, "dif")) {
    switch(
      metric,
      power = "DIFPower",
      false_positive = "DIFFalsePositiveRate",
      estimate = "MeanTargetContrast"
    )
  } else {
    switch(
      metric,
      power = "BiasScreenRate",
      false_positive = "BiasScreenFalsePositiveRate",
      estimate = "MeanTargetBias"
    )
  }
}

#' Evaluate DIF power and bias-screening behavior under known simulated signals
#'
#' @param n_person Vector of person counts to evaluate.
#' @param n_rater Vector of rater counts to evaluate.
#' @param n_criterion Vector of criterion counts to evaluate.
#' @param raters_per_person Vector of rater assignments per person.
#' @param design Optional named design-grid override supplied as a named list,
#'   named vector, or one-row data frame. Names may use canonical variables
#'   (`n_person`, `n_rater`, `n_criterion`, `raters_per_person`), current
#'   public aliases implied by `sim_spec` (for example `n_judge`, `n_task`,
#'   `judge_per_person`), or role keywords (`person`, `rater`, `criterion`,
#'   `assignment`). Values may be vectors. The schema-only future branch input
#'   `design$facets = c(person = ..., judge = ..., task = ...)` is also
#'   accepted for the currently exposed facet keys. Do not specify the same
#'   variable through both `design` and the scalar design-grid arguments.
#' @param reps Number of replications per design condition.
#' @param group_levels Group labels used for DIF simulation. The first two levels
#'   define the default reference and focal groups.
#' @param reference_group Optional reference group label used when extracting the
#'   target DIF contrast.
#' @param focal_group Optional focal group label used when extracting the target
#'   DIF contrast.
#' @param dif_level Target criterion level for the true DIF effect. Can be an
#'   integer index or a criterion label such as `"C04"`. Defaults to the last
#'   criterion level in each design.
#' @param dif_effect True DIF effect size added to the focal group on the target
#'   criterion.
#' @param bias_rater Target rater level for the true interaction-bias effect.
#'   Can be an integer index or a label such as `"R04"`. Defaults to the last
#'   rater level in each design.
#' @param bias_criterion Target criterion level for the true interaction-bias
#'   effect. Can be an integer index or a criterion label. Defaults to the last
#'   criterion level in each design.
#' @param bias_effect True interaction-bias effect added to the target
#'   `Rater x Criterion` cell.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread of step thresholds on the logit scale.
#' @param fit_method Estimation method passed to [fit_mfrm()].
#' @param model Measurement model passed to [fit_mfrm()]. The current
#'   signal-detection evaluator supports `RSM` and `PCM`; bounded `GPCM`
#'   is accepted only to produce an explicit unsupported-path error.
#' @param step_facet Step facet passed to [fit_mfrm()] when `model = "PCM"`.
#'   When left `NULL`, the function inherits the generator step facet from
#'   `sim_spec` when available and otherwise defaults to `"Criterion"`.
#' @param maxit Maximum iterations passed to [fit_mfrm()].
#' @param quad_points Quadrature points for `fit_method = "MML"`.
#' @param residual_pca Residual PCA mode passed to [diagnose_mfrm()].
#' @param sim_spec Optional output from [build_mfrm_sim_spec()] or
#'   [extract_mfrm_sim_spec()] used as the base data-generating mechanism.
#'   When supplied, the design grid still varies `n_person`, `n_rater`,
#'   `n_criterion`, and `raters_per_person`, but latent spread, thresholds,
#'   and other generator settings come from `sim_spec`. The target DIF and
#'   interaction-bias signals specified in this function override any signal
#'   tables stored in `sim_spec`. If `sim_spec` stores an active
#'   latent-regression population generator, this helper currently requires
#'   `fit_method = "MML"` so each replication can refit the population model.
#' @param dif_method Differential-functioning method passed to [analyze_dff()].
#' @param dif_min_obs Minimum observations per group cell for [analyze_dff()].
#' @param dif_p_adjust P-value adjustment method passed to [analyze_dff()].
#' @param dif_p_cut P-value cutoff for counting a target DIF detection.
#' @param dif_abs_cut Optional absolute contrast cutoff used when counting a
#'   target DIF detection. When omitted, the effective default is `0.43` for
#'   `dif_method = "refit"` and `0` (no additional magnitude cutoff) for
#'   `dif_method = "residual"`.
#' @param bias_max_iter Maximum iterations passed to [estimate_bias()].
#' @param bias_p_cut P-value cutoff for counting a target bias screen-positive result.
#' @param bias_abs_t Absolute t cutoff for counting a target bias screen-positive result.
#' @param seed Optional seed for reproducible replications.
#'
#' @details
#' This function performs Monte Carlo design screening for two related tasks:
#' DIF detection via [analyze_dff()] and interaction-bias screening via
#' [estimate_bias()].
#'
#' For each design condition (combination of `n_person`, `n_rater`,
#' `n_criterion`, `raters_per_person`), the function:
#' 1. Generates synthetic data with [simulate_mfrm_data()]
#' 2. Injects one known Group \eqn{\times} Criterion DIF effect
#'    (`dif_effect` logits added to the focal group on the target criterion)
#' 3. Injects one known Rater \eqn{\times} Criterion interaction-bias
#'    effect (`bias_effect` logits)
#' 4. Fits and diagnoses the MFRM
#' 5. Runs [analyze_dff()] and [estimate_bias()]
#' 6. Records whether the injected signals were detected or screen-positive
#'
#' **Detection criteria**:
#' A DIF signal is counted as "detected" when the target contrast has
#' \eqn{p <} `dif_p_cut` **and**, when an absolute contrast cutoff is in
#' force, \eqn{|\mathrm{Contrast}| \ge} `dif_abs_cut`. For
#' `dif_method = "refit"`, `dif_abs_cut` is interpreted on the logit scale.
#' For `dif_method = "residual"`, the residual-contrast screening result is
#' used and the default is to rely on the significance test alone.
#'
#' Bias results are different: [estimate_bias()] reports `t` and `Prob.` as
#' screening metrics rather than formal inferential quantities. Here, a bias
#' cell is counted as **screen-positive** only when those screening metrics are
#' available and satisfy
#'
#' First-release `GPCM` is not yet available in this helper because its signal-
#' detection path still depends on simulation and diagnostics layers validated
#' only for `RSM` / `PCM`. More broadly, the current planning layer is still
#' role-based for exactly two non-person facets (`rater`-like and
#' `criterion`-like), even though the estimation core supports arbitrary facet
#' counts.
#' \eqn{p <} `bias_p_cut` **and** \eqn{|t| \ge} `bias_abs_t`.
#'
#' **Power** is the proportion of replications in which the target signal
#' was correctly detected. For DIF this is a conventional power summary.
#' For bias, the primary summary is `BiasScreenRate`, a screening hit rate
#' rather than formal inferential power.
#'
#' **False-positive rate** is the proportion of non-target cells that were
#' incorrectly flagged. For DIF this is interpreted in the usual testing
#' sense. For bias, `BiasScreenFalsePositiveRate` is a screening rate and
#' should not be read as a calibrated inferential alpha level.
#'
#' **Default effect sizes**: `dif_effect = 0.6` logits corresponds to a
#' moderate criterion-linked differential-functioning effect; `bias_effect = -0.8`
#' logits represents a substantial rater-criterion interaction.  Adjust
#' these to match the smallest effect size of practical concern for your
#' application.
#'
#' This is again a **parametric simulation study**. The function does not
#' estimate a new design directly from one observed dataset. Instead, it
#' evaluates detection or screening behavior under user-specified design
#' conditions and known injected signals.
#'
#' If you want to approximate a real study, choose the design grid and
#' simulation settings so that they reflect the empirical context of interest.
#' For example, you may set `n_person`, `n_rater`, `n_criterion`,
#' `raters_per_person`, and the latent-spread arguments to values motivated by
#' an existing assessment program, then study how operating characteristics
#' change as those design settings vary.
#'
#' When `sim_spec` is supplied, the function uses it as the explicit
#' data-generating mechanism for the latent spreads, thresholds, and assignment
#' archetype, while still injecting the requested target DIF and bias effects
#' for each design condition.
#'
#' If that specification also stores a latent-regression population generator,
#' each replication carries simulated one-row-per-person background data into
#' the MML fit. This remains a screening-oriented Monte Carlo study; it is not
#' a person-level posterior prediction for one observed sample.
#'
#' @section References:
#' The simulation logic follows the general Monte Carlo / operating-characteristic
#' framework described by Morris, White, and Crowther (2019) and the
#' ADEMP-oriented planning/reporting guidance summarized for psychology by
#' Siepe et al. (2024). In `mfrmr`, `evaluate_mfrm_signal_detection()` is a
#' many-facet screening helper specialized to DIF and interaction-bias use
#' cases; it is not a direct implementation of one published many-facet Rasch
#' simulation design.
#'
#' - Morris, T. P., White, I. R., & Crowther, M. J. (2019).
#'   *Using simulation studies to evaluate statistical methods*.
#'   Statistics in Medicine, 38(11), 2074-2102.
#' - Siepe, B. S., Bartos, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W.,
#'   & Pawel, S. (2024). *Simulation studies for methodological research in
#'   psychology: A standardized template for planning, preregistration, and
#'   reporting*. Psychological Methods.
#'
#' @return An object of class `mfrm_signal_detection` with:
#' - `design_grid`: evaluated design conditions. When `sim_spec` carries custom
#'   public facet names, matching design-variable alias columns are included
#'   alongside the canonical internal columns.
#' - `results`: replicate-level detection results, with the same
#'   design-variable alias columns when applicable.
#' - `rep_overview`: run-level status and timing, with the same design-variable
#'   alias columns when applicable.
#' - `design_descriptor`: role-based design-variable metadata used by planning
#'   summaries and plots
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of which design variables remain
#'   mutable under the current simulation specification
#' - `planning_schema`: combined planner-schema contract bundling the role
#'   descriptor, scope boundary, and current mutability map
#' - `settings`: signal-analysis settings
#' - `ademp`: simulation-study metadata (aims, DGM, estimands, methods, performance measures)
#' @seealso [simulate_mfrm_data()], [evaluate_mfrm_design()], [analyze_dff()], [analyze_dif()], [estimate_bias()]
#' @examples
#' \donttest{
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   design = list(person = 8, rater = 2, criterion = 2, assignment = 1),
#'   reps = 1,
#'   maxit = 5,
#'   bias_max_iter = 1,
#'   seed = 123
#' ))
#' s_sig <- summary(sig_eval)
#' s_sig$overview
#' }
#' @export
evaluate_mfrm_signal_detection <- function(n_person = c(30, 50, 100),
                                           n_rater = c(4),
                                           n_criterion = c(4),
                                           raters_per_person = n_rater,
                                           design = NULL,
                                           reps = 10,
                                           group_levels = c("A", "B"),
                                           reference_group = NULL,
                                           focal_group = NULL,
                                           dif_level = NULL,
                                           dif_effect = 0.6,
                                           bias_rater = NULL,
                                           bias_criterion = NULL,
                                           bias_effect = -0.8,
                                           score_levels = 4,
                                           theta_sd = 1,
                                           rater_sd = 0.35,
                                           criterion_sd = 0.25,
                                           noise_sd = 0,
                                           step_span = 1.4,
                                           fit_method = c("JML", "MML"),
                                           model = c("RSM", "PCM", "GPCM"),
                                           step_facet = NULL,
                                           maxit = 25,
                                           quad_points = 7,
                                           residual_pca = c("none", "overall", "facet", "both"),
                                           sim_spec = NULL,
                                           dif_method = c("residual", "refit"),
                                           dif_min_obs = 10,
                                           dif_p_adjust = "holm",
                                           dif_p_cut = 0.05,
                                           dif_abs_cut = 0.43,
                                           bias_max_iter = 2,
                                           bias_p_cut = 0.05,
                                           bias_abs_t = 2,
                                           seed = NULL) {
  dif_abs_cut_missing <- missing(dif_abs_cut)
  fit_method <- match.arg(fit_method)
  model <- match.arg(model)
  if (identical(model, "GPCM")) {
    stop(
      "`evaluate_mfrm_signal_detection()` does not yet support bounded `GPCM`. ",
      "Signal-detection studies still rely on simulation and diagnostics layers that remain validated only for `RSM` / `PCM`. ",
      gpcm_planning_scope_rationale(),
      call. = FALSE
    )
  }
  residual_pca <- match.arg(residual_pca)
  dif_method <- match.arg(dif_method)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!is.null(sim_spec) && identical(as.character(sim_spec$model %||% NA_character_), "GPCM")) {
    stop(
      "`evaluate_mfrm_signal_detection()` does not yet support bounded `GPCM` simulation specifications. ",
      "Direct data generation is available, but signal-detection studies still depend on diagnostics layers validated only for `RSM` / `PCM`. ",
      gpcm_planning_scope_rationale(),
      call. = FALSE
    )
  }
  if (!is.null(sim_spec) &&
      isTRUE((sim_spec$population %||% simulation_empty_population_spec())$active) &&
      !identical(fit_method, "MML")) {
    stop(
      "Simulation specifications with an active latent-regression population generator currently require `fit_method = \"MML\"`.",
      call. = FALSE
    )
  }
  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) stop("`reps` must be >= 1.", call. = FALSE)
  supplied_counts <- intersect(
    names(as.list(match.call(expand.dots = FALSE))[-1]),
    c("n_person", "n_rater", "n_criterion", "raters_per_person")
  )
  dif_abs_cut <- as.numeric(dif_abs_cut[1])
  if (!is.finite(dif_abs_cut) || dif_abs_cut < 0) {
    stop("`dif_abs_cut` must be a single non-negative numeric value.", call. = FALSE)
  }
  dif_abs_cut_effective <- if (dif_abs_cut_missing && identical(dif_method, "residual")) {
    0
  } else {
    dif_abs_cut
  }

  group_levels <- unique(as.character(group_levels))
  group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
  if (length(group_levels) < 2L) {
    stop("`group_levels` must contain at least two non-empty labels.", call. = FALSE)
  }
  if (is.null(reference_group)) reference_group <- group_levels[1]
  if (is.null(focal_group)) focal_group <- group_levels[2]
  if (!reference_group %in% group_levels || !focal_group %in% group_levels) {
    stop("`reference_group` and `focal_group` must be members of `group_levels`.", call. = FALSE)
  }
  if (identical(reference_group, focal_group)) {
    stop("`reference_group` and `focal_group` must differ.", call. = FALSE)
  }

  design_meta <- simulation_build_design_grid(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    sim_spec = sim_spec,
    id_prefix = "S",
    explicit_scalar_names = supplied_counts
  )
  design_grid <- design_meta$canonical
  design_grid_public <- design_meta$public

  generator_model <- if (is.null(sim_spec)) model else sim_spec$model
  base_facet_names <- if (is.null(sim_spec)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(sim_spec)
  generator_step_facet <- if (is.null(sim_spec)) if (identical(generator_model, "PCM")) base_facet_names[2] else NA_character_ else sim_spec$step_facet
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
  design_variable_aliases <- design_meta$aliases
  design_descriptor <- design_meta$descriptor
  planning_scope <- simulation_planning_scope(sim_spec)
  planning_constraints <- simulation_planning_constraints(sim_spec)
  planning_schema <- simulation_planning_schema(sim_spec)
  fit_step_facet <- simulation_resolve_fit_step_facet(model, step_facet, generator_step_facet)

  seeds <- with_preserved_rng_seed(
    seed,
    sample.int(.Machine$integer.max, size = nrow(design_grid) * reps, replace = FALSE)
  )
  seed_idx <- 0L

  result_rows <- vector("list", nrow(design_grid) * reps)
  rep_rows <- vector("list", nrow(design_grid) * reps)
  out_idx <- 0L

  for (i in seq_len(nrow(design_grid))) {
    design <- design_grid[i, , drop = FALSE]
    row_facet_names <- if (is.null(sim_spec)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(sim_spec)
    dif_levels <- simulation_spec_role_levels(sim_spec, "criterion", count = design$n_criterion)
    bias_rater_levels <- simulation_spec_role_levels(sim_spec, "rater", count = design$n_rater)
    bias_criterion_levels <- simulation_spec_role_levels(sim_spec, "criterion", count = design$n_criterion)
    dif_target <- signal_eval_resolve_level_from_choices(dif_level, dif_levels, "dif_level")
    bias_rater_target <- signal_eval_resolve_level_from_choices(bias_rater, bias_rater_levels, "bias_rater")
    bias_criterion_target <- signal_eval_resolve_level_from_choices(bias_criterion, bias_criterion_levels, "bias_criterion")
    dif_tbl <- tibble::tibble(
      Group = focal_group,
      Effect = as.numeric(dif_effect[1])
    )
    dif_tbl[[row_facet_names[2]]] <- dif_target
    dif_tbl <- tibble::as_tibble(dif_tbl)
    dif_tbl <- dif_tbl[, c("Group", row_facet_names[2], "Effect")]
    bias_tbl <- tibble::tibble(Effect = as.numeric(bias_effect[1]))
    bias_tbl[[row_facet_names[1]]] <- bias_rater_target
    bias_tbl[[row_facet_names[2]]] <- bias_criterion_target
    bias_tbl <- tibble::as_tibble(bias_tbl)
    bias_tbl <- bias_tbl[, c(row_facet_names, "Effect")]
    row_spec <- if (is.null(sim_spec)) {
      NULL
    } else {
      simulation_override_spec_design(
        sim_spec,
        n_person = design$n_person,
        n_rater = design$n_rater,
        n_criterion = design$n_criterion,
        raters_per_person = design$raters_per_person,
        group_levels = group_levels,
        dif_effects = if (isTRUE(abs(as.numeric(dif_effect[1])) > 0)) dif_tbl else NULL,
        interaction_effects = if (isTRUE(abs(as.numeric(bias_effect[1])) > 0)) bias_tbl else NULL
      )
    }
    row_score_levels <- if (is.null(row_spec)) score_levels else row_spec$score_levels
    for (rep in seq_len(reps)) {
      seed_idx <- seed_idx + 1L
      sim <- if (is.null(row_spec)) {
        simulate_mfrm_data(
          n_person = design$n_person,
          n_rater = design$n_rater,
          n_criterion = design$n_criterion,
          raters_per_person = design$raters_per_person,
          score_levels = score_levels,
          theta_sd = theta_sd,
          rater_sd = rater_sd,
          criterion_sd = criterion_sd,
          noise_sd = noise_sd,
          step_span = step_span,
          group_levels = group_levels,
          dif_effects = if (isTRUE(abs(as.numeric(dif_effect[1])) > 0)) dif_tbl else NULL,
          interaction_effects = if (isTRUE(abs(as.numeric(bias_effect[1])) > 0)) bias_tbl else NULL,
          seed = seeds[seed_idx]
        )
      } else {
        simulate_mfrm_data(sim_spec = row_spec, seed = seeds[seed_idx])
      }

      t0 <- proc.time()[["elapsed"]]
      fit_args <- list(
        data = sim,
        person = "Person",
        facets = row_facet_names,
        score = "Score",
        method = fit_method,
        model = model,
        maxit = maxit
      )
      if (identical(model, "PCM")) fit_args$step_facet <- fit_step_facet
      if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
      if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points
      sim_population <- attr(sim, "mfrm_population_data")
      if (is.list(sim_population) && isTRUE(sim_population$active)) {
        fit_args$population_formula <- sim_population$population_formula
        fit_args$person_data <- sim_population$person_data
        fit_args$person_id <- sim_population$person_id
        fit_args$population_policy <- sim_population$population_policy
      }

      fit <- tryCatch(do.call(fit_mfrm, fit_args), error = function(e) e)
      diag <- if (inherits(fit, "error")) fit else {
        tryCatch(diagnose_mfrm(fit, residual_pca = residual_pca), error = function(e) e)
      }
      dif <- if (inherits(diag, "error")) diag else {
        tryCatch(
          analyze_dff(
            fit, diag,
            facet = row_facet_names[2],
            group = "Group",
            data = sim,
            method = dif_method,
            min_obs = dif_min_obs,
            p_adjust = dif_p_adjust
          ),
          error = function(e) e
        )
      }
      bias <- if (inherits(diag, "error")) diag else {
        tryCatch(
          estimate_bias(fit, diag, facet_a = row_facet_names[1], facet_b = row_facet_names[2], max_iter = bias_max_iter),
          error = function(e) e
        )
      }
      elapsed <- proc.time()[["elapsed"]] - t0
      converged <- !inherits(fit, "error") && isTRUE(as.logical(fit$summary$Converged[1]))

      err_msg <- character(0)
      if (inherits(fit, "error")) err_msg <- c(err_msg, conditionMessage(fit))
      if (inherits(diag, "error")) err_msg <- c(err_msg, conditionMessage(diag))
      if (inherits(dif, "error")) err_msg <- c(err_msg, conditionMessage(dif))
      if (inherits(bias, "error")) err_msg <- c(err_msg, conditionMessage(bias))
      run_ok <- length(err_msg) == 0L

      out_idx <- out_idx + 1L
      rep_rows[[out_idx]] <- tibble::tibble(
        design_id = design$design_id,
        rep = rep,
        n_person = design$n_person,
        n_rater = design$n_rater,
        n_criterion = design$n_criterion,
        raters_per_person = design$raters_per_person,
        Observations = nrow(sim),
        MinCategoryCount = min(tabulate(sim$Score, nbins = row_score_levels)),
        ElapsedSec = elapsed,
        RunOK = run_ok,
        Converged = converged,
        Error = if (length(err_msg) == 0L) NA_character_ else paste(unique(err_msg), collapse = " | ")
      )

      dif_target_row <- if (!inherits(dif, "error")) {
        signal_eval_find_dif_row(dif[["dif_table"]], dif_target, reference_group, focal_group)
      } else {
        tibble::tibble()
      }
      dif_target_p <- signal_eval_get_p(dif_target_row)
      dif_target_contrast <- if (nrow(dif_target_row) > 0) suppressWarnings(as.numeric(dif_target_row$Contrast[1])) else NA_real_
      dif_target_ets <- if (nrow(dif_target_row) > 0) as.character(dif_target_row$ETS[1]) else NA_character_
      dif_target_class <- if (nrow(dif_target_row) > 0) as.character(dif_target_row$Classification[1]) else NA_character_
      dif_target_class_system <- if (nrow(dif_target_row) > 0) as.character(dif_target_row$ClassificationSystem[1]) else NA_character_
      dif_detected <- is.finite(dif_target_p) && dif_target_p <= dif_p_cut &&
        is.finite(dif_target_contrast) && abs(dif_target_contrast) >= dif_abs_cut_effective
      dif_class_detected <- if (identical(dif_target_class_system, "ETS")) {
        !is.na(dif_target_ets) && dif_target_ets %in% c("B", "C")
      } else {
        identical(dif_target_class, "Screen positive")
      }

      dif_fp_rate <- NA_real_
      if (!inherits(dif, "error")) {
        dif_non_target <- tibble::as_tibble(dif[["dif_table"]]) |>
          dplyr::filter(
            .data$Level != dif_target,
            (.data$Group1 == reference_group & .data$Group2 == focal_group) |
              (.data$Group1 == focal_group & .data$Group2 == reference_group)
          ) |>
          dplyr::mutate(
            p_eval = dplyr::if_else(is.finite(.data$p_adjusted), .data$p_adjusted, .data$p_value),
            Flag = is.finite(.data$p_eval) & .data$p_eval <= dif_p_cut &
              is.finite(.data$Contrast) & abs(.data$Contrast) >= dif_abs_cut_effective
          )
        dif_fp_rate <- signal_eval_false_positive_rate(dif_non_target$Flag)
      }

      bias_target_row <- if (!inherits(bias, "error")) {
        signal_eval_find_bias_row(bias[["table"]], bias_rater_target, bias_criterion_target)
      } else {
        tibble::tibble()
      }
      bias_target_p <- if (nrow(bias_target_row) > 0) suppressWarnings(as.numeric(bias_target_row$`Prob.`[1])) else NA_real_
      bias_target_t <- if (nrow(bias_target_row) > 0) suppressWarnings(as.numeric(bias_target_row$t[1])) else NA_real_
      bias_target_size <- if (nrow(bias_target_row) > 0) suppressWarnings(as.numeric(bias_target_row$`Bias Size`[1])) else NA_real_
      bias_metric_available <- is.finite(bias_target_p) && is.finite(bias_target_t)
      bias_detected <- is.finite(bias_target_p) && bias_target_p <= bias_p_cut &&
        is.finite(bias_target_t) && abs(bias_target_t) >= bias_abs_t

      bias_fp_rate <- NA_real_
      if (!inherits(bias, "error")) {
        bias_non_target <- tibble::as_tibble(bias[["table"]]) |>
          dplyr::filter(!(.data$FacetA_Level == bias_rater_target & .data$FacetB_Level == bias_criterion_target)) |>
          dplyr::mutate(
            Flag = is.finite(.data$`Prob.`) & .data$`Prob.` <= bias_p_cut &
              is.finite(.data$t) & abs(.data$t) >= bias_abs_t
          )
        bias_fp_rate <- signal_eval_false_positive_rate(bias_non_target$Flag)
      }

      result_rows[[out_idx]] <- tibble::tibble(
        design_id = design$design_id,
        rep = rep,
        n_person = design$n_person,
        n_rater = design$n_rater,
        n_criterion = design$n_criterion,
        raters_per_person = design$raters_per_person,
        Observations = nrow(sim),
        ElapsedSec = elapsed,
        Converged = converged,
        DIFTargetLevel = dif_target,
        DIFContrast = dif_target_contrast,
        DIFP = dif_target_p,
        DIFClassificationSystem = dif_target_class_system,
        DIFClassification = dif_target_class,
        DIFETS = dif_target_ets,
        DIFDetected = dif_detected,
        DIFClassDetected = dif_class_detected,
        DIFFalsePositiveRate = dif_fp_rate,
        BiasTargetRater = bias_rater_target,
        BiasTargetCriterion = bias_criterion_target,
        BiasSize = bias_target_size,
        BiasP = bias_target_p,
        BiasT = bias_target_t,
        BiasScreenMetricAvailable = bias_metric_available,
        BiasDetected = bias_detected,
        BiasScreenFalsePositiveRate = bias_fp_rate
      )
    }
  }

  structure(
    list(
      design_grid = design_grid_public,
      results = simulation_append_design_alias_columns(dplyr::bind_rows(result_rows), design_variable_aliases),
      rep_overview = simulation_append_design_alias_columns(dplyr::bind_rows(rep_rows), design_variable_aliases),
      design_descriptor = design_descriptor,
      planning_scope = planning_scope,
      planning_constraints = planning_constraints,
      planning_schema = planning_schema,
      settings = list(
        group_levels = group_levels,
        reference_group = reference_group,
        focal_group = focal_group,
        dif_level = dif_level,
        dif_effect = dif_effect,
        bias_rater = bias_rater,
        bias_criterion = bias_criterion,
        bias_effect = bias_effect,
        reps = reps,
        score_levels = score_levels,
        theta_sd = theta_sd,
        rater_sd = rater_sd,
        criterion_sd = criterion_sd,
        noise_sd = noise_sd,
        step_span = step_span,
        fit_method = fit_method,
        model = model,
        step_facet = fit_step_facet,
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        dif_method = dif_method,
        dif_min_obs = dif_min_obs,
        dif_p_adjust = dif_p_adjust,
        dif_p_cut = dif_p_cut,
        dif_abs_cut = dif_abs_cut_effective,
        dif_abs_cut_input = dif_abs_cut,
        bias_max_iter = bias_max_iter,
        bias_p_cut = bias_p_cut,
        bias_abs_t = bias_abs_t,
        sim_spec = sim_spec,
        facet_names = stats::setNames(base_facet_names, c("rater", "criterion")),
        design_variable_aliases = design_variable_aliases,
        design_descriptor = design_descriptor,
        planning_scope = planning_scope,
        planning_constraints = planning_constraints,
        planning_schema = planning_schema,
        generator_model = generator_model,
        generator_step_facet = generator_step_facet,
        generator_assignment = generator_assignment,
        seed = seed
      ),
      ademp = simulation_build_ademp(
        purpose = "Assess DIF detection and interaction-bias screening behavior under repeated parametric many-facet simulations with known injected targets.",
        design_grid = design_grid,
        generator_model = generator_model,
        generator_step_facet = generator_step_facet,
        generator_assignment = generator_assignment,
        sim_spec = sim_spec,
        estimands = c(
          "DIF target-flag rate and non-target flag rate",
          "Bias screening hit rate and screening false-positive rate",
          "Target contrast and target bias summaries",
          "Convergence rate and elapsed time"
        ),
        analysis_methods = list(
          fit_method = fit_method,
          fitted_model = model,
          dif_method = dif_method,
          bias_method = "estimate_bias_screening",
          maxit = maxit,
          quad_points = if (identical(fit_method, "MML")) quad_points else NA_integer_,
          residual_pca = residual_pca
        ),
        performance_measures = c(
          "Mean detection/screening summaries across replications",
          "MCSE for means and rates",
          "Convergence rate",
          "Bias-screen metric availability rate"
        )
      )
    ),
    class = "mfrm_signal_detection"
  )
}

#' Summarize a DIF/bias screening simulation
#'
#' @param object Output from [evaluate_mfrm_signal_detection()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_signal_detection` with:
#' - `overview`: run-level overview
#' - `detection_summary`: aggregated detection rates by design, with
#'   design-variable alias columns when applicable
#' - `ademp`: simulation-study metadata carried forward from the original object
#' - `facet_names`: public facet labels carried from the simulation specification
#' - `design_variable_aliases`: accepted public aliases for design variables
#' - `design_descriptor`: role-based design-variable metadata
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract
#' - `future_branch_active_summary`: compact deterministic summary of the
#'   schema-only future arbitrary-facet planning branch embedded in the current
#'   planning schema
#' - `notes`: short interpretation notes, including the bias-side screening caveat
#' @seealso [evaluate_mfrm_signal_detection()], [plot.mfrm_signal_detection]
#' @examples
#' \dontrun{
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 8,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 5,
#'   bias_max_iter = 1,
#'   seed = 123
#' ))
#' summary(sig_eval)
#' }
#' @export
summary.mfrm_signal_detection <- function(object, digits = 3, ...) {
  if (!is.list(object) || is.null(object$results) || is.null(object$rep_overview)) {
    stop("`object` must be output from evaluate_mfrm_signal_detection().")
  }
  digits <- max(0L, as.integer(digits[1]))
  out <- signal_eval_summary(
    object$results,
    object$rep_overview,
    design_variable_aliases = simulation_object_design_variable_aliases(object)
  )

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out$overview <- round_df(out$overview)
  out$detection_summary <- round_df(out$detection_summary)
  out$ademp <- object$ademp %||% NULL
  out$facet_names <- object$settings$facet_names %||% stats::setNames(simulation_default_output_facet_names(), c("rater", "criterion"))
  out$design_variable_aliases <- simulation_object_design_variable_aliases(object)
  out$design_descriptor <- simulation_object_design_descriptor(object)
  out$planning_scope <- simulation_object_planning_scope(object)
  out$planning_constraints <- simulation_object_planning_constraints(object)
  out$planning_schema <- simulation_object_planning_schema(object)
  out$future_branch_active_summary <- simulation_compact_future_branch_active_summary(
    object,
    digits = digits
  )
  out$digits <- digits
  scope_note <- simulation_planning_scope_note(out$planning_scope)
  if (length(scope_note) > 0L && !scope_note %in% out$notes) {
    out$notes <- c(out$notes, scope_note)
  }
  constraint_note <- simulation_planning_constraints_note(out$planning_constraints)
  if (length(constraint_note) > 0L && !constraint_note %in% out$notes) {
    out$notes <- c(out$notes, constraint_note)
  }
  schema_note <- simulation_planning_schema_note(out$planning_schema)
  if (length(schema_note) > 0L && !schema_note %in% out$notes) {
    out$notes <- c(out$notes, schema_note)
  }
  if (inherits(out$future_branch_active_summary, "summary.mfrm_future_branch_active_branch")) {
    out$notes <- c(
      out$notes,
      "A deterministic future arbitrary-facet planning scaffold is embedded in `future_branch_active_summary`; it reports structural bookkeeping and conservative recommendation logic, not DIF/bias detection power."
    )
  }
  out$notes <- unique(out$notes)
  class(out) <- "summary.mfrm_signal_detection"
  out
}

#' @export
print.summary.mfrm_signal_detection <- function(x, ...) {
  digits <- max(0L, as.integer(x$digits %||% 3L))
  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }
  preview_df <- function(df, n = 10L) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    utils::head(df, n = n)
  }

  cat("mfrmr Signal Detection Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_df(as.data.frame(x$overview)), row.names = FALSE)
  }
  if (!is.null(x$detection_summary) && nrow(x$detection_summary) > 0) {
    cat("\nDetection summary (preview)\n")
    print(round_df(as.data.frame(preview_df(x$detection_summary))), row.names = FALSE)
  }
  print_compact_future_branch_active_summary(
    x$future_branch_active_summary %||% NULL,
    digits = digits
  )
  if (is.list(x$ademp) && length(x$ademp) > 0L) {
    cat("\nADEMP metadata\n")
    cat(" - aims\n")
    cat(" - data_generating_mechanism\n")
    cat(" - estimands\n")
    cat(" - methods\n")
    cat(" - performance_measures\n")
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

signal_detection_metric_label <- function(signal, metric_col) {
  if (identical(signal, "dif")) {
    return(switch(
      metric_col,
      DIFPower = "DIF target-flag rate",
      DIFFalsePositiveRate = "DIF non-target flag rate",
      MeanDIFEstimate = "Mean target contrast",
      metric_col
    ))
  }
  switch(
    metric_col,
    BiasScreenRate = "Bias screening hit rate",
    BiasScreenFalsePositiveRate = "Bias screening false-positive rate",
    MeanBiasEstimate = "Mean bias estimate",
    metric_col
  )
}

#' Plot DIF/bias screening simulation results
#'
#' @param x Output from [evaluate_mfrm_signal_detection()].
#' @param signal Whether to plot DIF or bias screening results.
#' @param metric Metric to plot. For `signal = "bias"`, prefer
#'   `metric = "screen_rate"` for the screening hit rate. The older
#'   `metric = "power"` spelling is retained as a backwards-compatible alias
#'   that maps to `BiasScreenRate`.
#' @param x_var Design variable used on the x-axis. When `x` was generated from
#'   a `sim_spec` with custom public facet names, the corresponding aliases
#'   (for example `n_judge`, `n_task`, `judge_per_person`) are also accepted.
#'   Role keywords (`person`, `rater`, `criterion`, `assignment`) are accepted
#'   as an abstraction over the current two-facet schema.
#' @param group_var Optional design variable used for separate lines. The same
#'   alias rules as `x_var` apply.
#' @param draw If `TRUE`, draw with base graphics; otherwise return plotting data.
#' @param ... Reserved for generic compatibility.
#'
#' @return If `draw = TRUE`, invisibly returns plotting data. If `draw = FALSE`,
#'   returns that plotting-data list directly. The returned list includes
#'   resolved canonical variables (`x_var`, `group_var`) together with public
#'   labels (`x_label`, `group_label`), `design_variable_aliases`,
#'   `design_descriptor`, `planning_scope`, `planning_constraints`,
#'   `planning_schema`,
#'   `display_metric`, and `interpretation_note` so
#'   callers can label bias-side plots as screening summaries rather than
#'   formal power/error-rate displays.
#' @seealso [evaluate_mfrm_signal_detection()], [summary.mfrm_signal_detection]
#' @examples
#' \dontrun{
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 8,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 5,
#'   bias_max_iter = 1,
#'   seed = 123
#' ))
#' plot(sig_eval, signal = "dif", metric = "power", x_var = "n_person", draw = FALSE)
#' }
#' @export
plot.mfrm_signal_detection <- function(x,
                                       signal = c("dif", "bias"),
                                       metric = c("power", "false_positive", "estimate",
                                                  "screen_rate", "screen_false_positive"),
                                       x_var = c("n_person", "n_rater", "n_criterion", "raters_per_person"),
                                       group_var = NULL,
                                       draw = TRUE,
                                       ...) {
  if (!is.list(x) || is.null(x$results) || is.null(x$rep_overview)) {
    stop("`x` must be output from evaluate_mfrm_signal_detection().")
  }
  signal <- match.arg(signal)
  metric <- match.arg(metric)
  if (identical(signal, "bias")) {
    metric <- switch(
      metric,
      screen_rate = "power",
      screen_false_positive = "false_positive",
      metric
    )
  }
  design_variable_aliases <- simulation_object_design_variable_aliases(x)
  design_descriptor <- simulation_object_design_descriptor(x)
  x_var <- if (missing(x_var)) {
    "n_person"
  } else {
    simulation_resolve_design_variable(x_var, design_variable_aliases, "x_var", descriptor = design_descriptor)
  }
  x_label <- simulation_design_variable_label(x_var, design_variable_aliases)

  sum_obj <- signal_eval_summary(
    x$results,
    x$rep_overview,
    design_variable_aliases = design_variable_aliases
  )
  plot_tbl <- tibble::as_tibble(sum_obj$detection_summary)
  if (nrow(plot_tbl) == 0) stop("No detection-summary rows available for plotting.")

  metric_col <- signal_eval_metric_col(signal, metric)
  varying <- simulation_design_canonical_variables(design_descriptor)
  varying <- varying[varying != x_var]
  if (is.null(group_var)) {
    cand <- varying[vapply(plot_tbl[varying], function(col) length(unique(col)) > 1L, logical(1))]
    group_var <- if (length(cand) > 0) cand[1] else NULL
  } else {
    group_var <- simulation_resolve_design_variable(group_var, design_variable_aliases, "group_var", descriptor = design_descriptor)
    if (identical(group_var, x_var)) {
      stop("`group_var` must differ from `x_var`.")
    }
  }
  group_label <- if (is.null(group_var)) NULL else simulation_design_variable_label(group_var, design_variable_aliases)

  if (is.null(group_var)) {
    agg_tbl <- plot_tbl |>
      dplyr::group_by(.data[[x_var]]) |>
      dplyr::summarize(y = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data[[x_var]]) |>
      dplyr::mutate(group = "All designs")
  } else {
    agg_tbl <- plot_tbl |>
      dplyr::group_by(.data[[x_var]], .data[[group_var]]) |>
      dplyr::summarize(y = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data[[x_var]], .data[[group_var]]) |>
      dplyr::rename(group = dplyr::all_of(group_var))
  }

  out <- list(
    plot = "signal_detection",
    signal = signal,
    metric = metric,
    metric_col = metric_col,
    display_metric = signal_detection_metric_label(signal, metric_col),
    interpretation_note = if (identical(signal, "bias")) {
      "Bias-side rates summarize screening behavior from estimate_bias(); they are not formal inferential power or alpha estimates."
    } else {
      "DIF-side rates summarize target/non-target flagging behavior under the selected DFF method and threshold settings."
    },
    x_var = x_var,
    x_label = x_label,
    group_var = group_var,
    group_label = group_label,
    design_variable_aliases = design_variable_aliases,
    design_descriptor = design_descriptor,
    planning_scope = simulation_object_planning_scope(x),
    planning_constraints = simulation_object_planning_constraints(x),
    planning_schema = simulation_object_planning_schema(x),
    data = agg_tbl
  )
  if (!isTRUE(draw)) return(out)

  groups <- unique(as.character(agg_tbl$group))
  cols <- grDevices::hcl.colors(max(1L, length(groups)), "Set 2")
  x_vals <- sort(unique(agg_tbl[[x_var]]))
  y_range <- range(agg_tbl$y, na.rm = TRUE)
  if (!all(is.finite(y_range))) stop("Selected metric has no finite values to plot.")

  graphics::plot(
    x = x_vals,
    y = rep(NA_real_, length(x_vals)),
    type = "n",
    xlab = x_label,
    ylab = out$display_metric,
    main = if (identical(signal, "bias")) {
      paste("Bias screening simulation:", out$display_metric)
    } else {
      paste("DIF screening simulation:", out$display_metric)
    },
    ylim = y_range
  )
  for (i in seq_along(groups)) {
    sub <- agg_tbl[as.character(agg_tbl$group) == groups[i], , drop = FALSE]
    sub <- sub[order(sub[[x_var]]), , drop = FALSE]
    graphics::lines(sub[[x_var]], sub$y, type = "b", lwd = 2, pch = 16 + (i - 1L) %% 5L, col = cols[i])
  }
  if (length(groups) > 1L) {
    graphics::legend("topleft", legend = groups, col = cols, lty = 1, lwd = 2, pch = 16 + (seq_along(groups) - 1L) %% 5L, bty = "n")
  }

  invisible(out)
}
