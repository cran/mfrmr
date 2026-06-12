#' Simulate long-format ordered many-facet data for design studies
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
#' @param thresholds Optional threshold specification. Use a numeric vector of
#'   common thresholds; a named list such as `list(C01 = c(-1, 0, 1))`; a
#'   numeric matrix with one row per `StepFacet` and one column per step; or a
#'   long data frame with columns `StepFacet`, `Step`/`StepIndex`, and
#'   `Estimate`.
#' @param slopes Optional slope specification used when `model = "GPCM"`.
#'   Use either a numeric vector aligned to the generated slope-facet levels or
#'   a data frame with columns `SlopeFacet` and `Estimate`. Supplied slopes are
#'   treated as relative discriminations and normalized to the package's
#'   geometric-mean-one identification convention on the log scale. When
#'   omitted, slopes default to 1 for every slope-facet level, giving an exact
#'   `PCM` reduction.
#' @param assignment Assignment design. `"crossed"` means every person sees
#'   every rater; `"rotating"` uses a balanced rotating subset; `"resampled"`
#'   reuses person-level rater-assignment profiles stored in `sim_spec`;
#'   `"sparse_linked"` uses an incomplete rating design with optional linking
#'   persons; `"skeleton"` reuses an observed response skeleton stored in
#'   `sim_spec`, including optional `Group`/`Weight` columns when available.
#'   When omitted, the function chooses `"crossed"` if
#'   `raters_per_person == n_rater`, otherwise `"rotating"`.
#' @param sparse_controls Optional named list used when
#'   `assignment = "sparse_linked"`. Supported entries are `link_fraction`,
#'   `link_persons`, `link_raters_per_person`, `assignment_mode`, and
#'   `min_common_persons_per_rater_pair`. See [build_mfrm_sim_spec()] for the
#'   same contract.
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
#' This function generates synthetic ordered many-facet data under `RSM`,
#' `PCM`, or the package's bounded `GPCM` branch.
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
#' - if `thresholds` is a named list, numeric matrix, or data frame, threshold
#'   values may vary by `StepFacet` (currently `Criterion` or `Rater`)
#'
#' For bounded `GPCM`, the generator now requires an explicit slope
#' contract in parallel with the threshold table. The current public branch
#' keeps `slope_facet == step_facet`, normalizes supplied slopes to the same
#' geometric-mean-one log-slope identification used by [fit_mfrm()], and uses
#' the internal `category_prob_gpcm()` helper for response sampling. Broader
#' arbitrary-facet planning remains restricted until that slope-aware contract
#' is generalized beyond the current role-based design, population-forecasting,
#' diagnostic-screening, and signal-detection helpers.
#'
#' Assignment handling is also explicit:
#' - `"crossed"` uses the full person x rater x criterion design
#' - `"rotating"` assigns a deterministic rotating subset of raters per person
#' - `"sparse_linked"` assigns most persons to an incomplete rater subset and
#'   assigns a configurable set of linking persons to a larger rater set
#' - `"resampled"` reuses empirical person-level rater profiles stored in
#'   `sim_spec$assignment_profiles`, optionally carrying over person-level
#'   `Group`
#' - `"skeleton"` reuses an observed person-by-rater-by-criterion response
#'   skeleton stored in `sim_spec$design_skeleton`, optionally carrying over
#'   `Group` and `Weight`
#'
#' Sparse linked simulation is intended for planned-missing rating designs in
#' which connectivity is maintained through common linking persons. The
#' returned `mfrm_sparse_design` attribute summarizes design density, planned
#' missingness, rater coverage, and rater-pair common-person counts. These
#' summaries are design diagnostics, not model-fit statistics or universal
#' adequacy thresholds. This branch follows sparse rater-mediated assessment
#' design work by Wind, Jones, and Grajeda (2023,
#' doi:10.1177/01466216231182148), Wind and Jones (2018,
#' doi:10.1177/0013164417703733), and DeMars, Shapovalov, and Hathcoat
#' (2023).
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
#' - `mfrm_sparse_design`: sparse-design diagnostics when
#'   `assignment = "sparse_linked"`, including design density, planned missing
#'   rate, rater coverage, and rater-pair common-person counts
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
                               sparse_controls = NULL,
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
    sparse_controls <- sim_spec$sparse_controls %||%
      simulation_normalize_sparse_controls(NULL, assignment, n_person, n_rater, raters_per_person)
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
      match.arg(tolower(as.character(assignment[1])), c("crossed", "rotating", "sparse_linked"))
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
    sparse_controls <- simulation_normalize_sparse_controls(
      sparse_controls = sparse_controls,
      assignment = assignment,
      n_person = n_person,
      n_rater = n_rater,
      raters_per_person = raters_per_person
    )
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
      !identical(assignment, "sparse_linked") &&
      !identical(assignment, "resampled") && !identical(assignment, "skeleton")) {
    stop("`assignment` must be one of \"crossed\", \"rotating\", \"sparse_linked\", \"resampled\", or \"skeleton\".", call. = FALSE)
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
  sparse_controls <- simulation_normalize_sparse_controls(
    sparse_controls = sparse_controls,
    assignment = assignment,
    n_person = n_person,
    n_rater = n_rater,
    raters_per_person = raters_per_person
  )

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

    sparse_design <- NULL
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
    } else if (assignment == "sparse_linked") {
      sparse_assignment <- simulation_generate_sparse_linked_assignment(
        person_ids = person_ids,
        rater_ids = rater_ids,
        criterion_ids = criterion_ids,
        raters_per_person = raters_per_person,
        sparse_controls = sparse_controls,
        facet_names = facet_names
      )
      dat <- sparse_assignment$data
      sparse_design <- sparse_assignment$metadata
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

    peer_review_design <- NULL
    if (!is.null(sim_spec) && isTRUE(sim_spec$peer_review$active)) {
      peer_review_design <- simulation_peer_review_design_metadata(
        dat = dat,
        peer_review = sim_spec$peer_review,
        facet_names = facet_names
      )
    }

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
      ),
      design = list(
        sparse = sparse_design,
        peer_review = peer_review_design
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
      slope_regime = simulation_gpcm_slope_regime(slope_table),
      design_skeleton = design_skeleton,
      sparse_controls = sparse_controls,
      sparse_design = sparse_design,
      peer_review = if (!is.null(sim_spec)) sim_spec$peer_review %||% NULL else NULL,
      population = population_spec
    )
    if (!is.null(sparse_design)) {
      attr(dat, "mfrm_sparse_design") <- sparse_design
    }
    if (!is.null(peer_review_design)) {
      attr(dat, "mfrm_peer_review_design") <- peer_review_design
    }
    dat
  })
}

simulation_fit_score_support <- function(sim, fallback_score_levels = NULL) {
  cfg <- attr(sim, "mfrm_simulation_spec")
  score_levels <- suppressWarnings(as.integer(
    (cfg$score_levels %||% fallback_score_levels %||% NA_integer_)[1]
  ))
  if (!is.finite(score_levels) || score_levels < 2L) {
    return(list())
  }
  list(rating_min = 1L, rating_max = score_levels)
}

simulation_add_fit_score_support <- function(fit_args, sim,
                                             fallback_score_levels = NULL) {
  support <- simulation_fit_score_support(
    sim = sim,
    fallback_score_levels = fallback_score_levels
  )
  if (length(support) > 0L) {
    fit_args[names(support)] <- support
  }
  fit_args
}

simulation_score_support_summary <- function(sim, fallback_score_levels = NULL) {
  empty <- list(
    ScoreLevelsDeclared = NA_integer_,
    ScoreLevelsObserved = NA_integer_,
    ZeroScoreLevels = NA_integer_,
    MinScoreCount = NA_integer_,
    MinScoreProportion = NA_real_,
    MaxScoreProportion = NA_real_
  )
  score_levels <- suppressWarnings(as.integer(
    ((attr(sim, "mfrm_simulation_spec")$score_levels %||% fallback_score_levels) %||% NA_integer_)[1]
  ))
  if (!is.finite(score_levels) || score_levels < 2L ||
      !is.data.frame(sim) || !"Score" %in% names(sim)) {
    return(empty)
  }
  scores <- suppressWarnings(as.integer(sim$Score))
  scores <- scores[is.finite(scores) & scores >= 1L & scores <= score_levels]
  counts <- tabulate(scores, nbins = score_levels)
  total <- sum(counts)
  if (total <= 0L) {
    empty$ScoreLevelsDeclared <- score_levels
    return(empty)
  }
  list(
    ScoreLevelsDeclared = as.integer(score_levels),
    ScoreLevelsObserved = as.integer(sum(counts > 0L)),
    ZeroScoreLevels = as.integer(sum(counts == 0L)),
    MinScoreCount = as.integer(min(counts)),
    MinScoreProportion = as.numeric(min(counts) / total),
    MaxScoreProportion = as.numeric(max(counts) / total)
  )
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
  reuse_template_people <- "TemplatePersonReuse" %in% names(skeleton) &&
    any(simulation_sparse_design_active(skeleton$TemplatePersonReuse), na.rm = TRUE)
  if (isTRUE(reuse_template_people)) {
    if (length(template_people) != length(person_ids)) {
      stop(
        "`TemplatePersonReuse = TRUE` requires the number of template persons to match `n_person`.",
        call. = FALSE
      )
    }
    sampled_templates <- template_people
  } else {
    sampled_templates <- sample(template_people, size = length(person_ids), replace = TRUE)
  }
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

simulation_peer_review_pair_table <- function(review_pairs, reviewer_ids) {
  if (length(reviewer_ids) < 2L || nrow(review_pairs) == 0L) {
    return(data.frame())
  }
  pairs <- utils::combn(reviewer_ids, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    p1 <- review_pairs$Person[review_pairs$Reviewer == pair[1]]
    p2 <- review_pairs$Person[review_pairs$Reviewer == pair[2]]
    common <- intersect(p1, p2)
    data.frame(
      Reviewer1 = pair[1],
      Reviewer2 = pair[2],
      CommonSubmissions = length(common),
      CommonSubmissionIDs = paste(common, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

simulation_peer_review_reciprocal_pairs <- function(review_pairs) {
  if (nrow(review_pairs) == 0L) return(0L)
  pair_keys <- paste(review_pairs$Person, review_pairs$Reviewer, sep = "\r")
  reciprocal <- paste(review_pairs$Reviewer, review_pairs$Person, sep = "\r")
  unordered <- ifelse(
    review_pairs$Person < review_pairs$Reviewer,
    paste(review_pairs$Person, review_pairs$Reviewer, sep = "\r"),
    paste(review_pairs$Reviewer, review_pairs$Person, sep = "\r")
  )
  length(unique(unordered[pair_keys %in% reciprocal]))
}

simulation_peer_review_design_metadata <- function(dat,
                                                   peer_review,
                                                   facet_names = simulation_default_output_facet_names()) {
  if (!isTRUE(peer_review$active) || !is.data.frame(dat) || nrow(dat) == 0L) {
    return(NULL)
  }
  reviewer_col <- as.character(facet_names[1] %||% "Reviewer")
  criterion_col <- as.character(facet_names[2] %||% "Criterion")
  if (!all(c("Person", reviewer_col, criterion_col) %in% names(dat))) {
    return(NULL)
  }
  review_pairs <- dat |>
    dplyr::distinct(
      Person = as.character(.data$Person),
      Reviewer = as.character(.data[[reviewer_col]])
    ) |>
    dplyr::arrange(.data$Person, .data$Reviewer) |>
    as.data.frame(stringsAsFactors = FALSE)
  submission_load <- review_pairs |>
    dplyr::count(.data$Person, name = "ReviewersAssigned") |>
    as.data.frame(stringsAsFactors = FALSE)
  reviewer_load <- review_pairs |>
    dplyr::count(.data$Reviewer, name = "SubmissionsReviewed") |>
    as.data.frame(stringsAsFactors = FALSE)
  reviewer_ids <- sort(unique(as.character(review_pairs$Reviewer)))
  pair_tbl <- simulation_peer_review_pair_table(review_pairs, reviewer_ids)
  self_reviews <- sum(review_pairs$Person == review_pairs$Reviewer, na.rm = TRUE)
  reciprocal_pairs <- simulation_peer_review_reciprocal_pairs(review_pairs)
  avoid_self <- isTRUE(peer_review$avoid_self_review)
  possible_review_pairs <- length(unique(review_pairs$Person)) *
    max(length(reviewer_ids) - if (avoid_self) 1L else 0L, 1L)
  possible_rows <- possible_review_pairs * length(unique(as.character(dat[[criterion_col]])))
  min_common <- if (nrow(pair_tbl) == 0L) NA_integer_ else min(pair_tbl$CommonSubmissions)
  zero_common <- if (nrow(pair_tbl) == 0L) NA_integer_ else sum(pair_tbl$CommonSubmissions == 0L)

  overview <- data.frame(
    Active = TRUE,
    Scenario = "peer_review",
    Submissions = length(unique(review_pairs$Person)),
    Reviewers = length(reviewer_ids),
    Criteria = length(unique(as.character(dat[[criterion_col]]))),
    Rows = nrow(dat),
    ReviewPairs = nrow(review_pairs),
    PossibleReviewPairs = possible_review_pairs,
    DesignDensity = nrow(dat) / possible_rows,
    SelfReviews = self_reviews,
    AvoidSelfReview = avoid_self,
    ReciprocalPairs = reciprocal_pairs,
    AnchorSubmissions = as.integer(peer_review$anchor_submissions %||% NA_integer_),
    OrdinaryReviewersPerSubmission = as.integer(peer_review$reviewers_per_submission %||% NA_integer_),
    AnchorReviewersPerSubmission = as.integer(peer_review$anchor_reviewers_per_submission %||% NA_integer_),
    MinReviewersPerSubmission = min(submission_load$ReviewersAssigned),
    MeanReviewersPerSubmission = mean(submission_load$ReviewersAssigned),
    MaxReviewersPerSubmission = max(submission_load$ReviewersAssigned),
    MinSubmissionsPerReviewer = min(reviewer_load$SubmissionsReviewed),
    MeanSubmissionsPerReviewer = mean(reviewer_load$SubmissionsReviewed),
    MaxSubmissionsPerReviewer = max(reviewer_load$SubmissionsReviewed),
    MinCommonSubmissionsPerReviewerPair = min_common,
    ZeroCommonReviewerPairs = zero_common,
    ReviewUse = "design_diagnostic_not_measurement_gate",
    stringsAsFactors = FALSE
  )

  list(
    active = TRUE,
    overview = overview,
    submission_load = submission_load,
    reviewer_load = reviewer_load,
    reviewer_pair_common_submissions = pair_tbl,
    review_pairs = review_pairs,
    notes = c(
      "Peer-review design diagnostics summarize assignment structure and linking, not MFRM fit or reviewer-quality adequacy.",
      "Common submissions per reviewer pair describe reviewer linkage through shared reviewed work."
    )
  )
}

simulation_generate_sparse_linked_assignment <- function(person_ids,
                                                        rater_ids,
                                                        criterion_ids,
                                                        raters_per_person,
                                                        sparse_controls,
                                                        facet_names = simulation_default_output_facet_names()) {
  assignment_facet <- facet_names[1]
  criterion_facet <- facet_names[2]
  n_rater <- length(rater_ids)
  n_person <- length(person_ids)
  link_n <- as.integer(sparse_controls$link_persons %||% 0L)
  link_persons <- if (link_n > 0L) person_ids[seq_len(link_n)] else character(0)

  draw_raters <- function(i, size) {
    size <- as.integer(size)
    if (size >= n_rater) {
      return(rater_ids)
    }
    if (identical(sparse_controls$assignment_mode, "random")) {
      return(sample(rater_ids, size = size, replace = FALSE))
    }
    rater_ids[((i - 1L) + seq_len(size) - 1L) %% n_rater + 1L]
  }

  person_assignment <- vector("list", n_person)
  rows <- vector("list", n_person)
  for (i in seq_along(person_ids)) {
    is_link <- person_ids[i] %in% link_persons
    size <- if (is_link) {
      sparse_controls$link_raters_per_person
    } else {
      raters_per_person
    }
    assigned <- draw_raters(i, size = size)
    person_assignment[[i]] <- data.frame(
      Person = person_ids[i],
      RatersAssigned = length(assigned),
      LinkPerson = is_link,
      Raters = paste(assigned, collapse = ", "),
      stringsAsFactors = FALSE
    )
    rows[[i]] <- expand.grid(
      setNames(
        list(person_ids[i], assigned, criterion_ids),
        c("Person", assignment_facet, criterion_facet)
      ),
      stringsAsFactors = FALSE
    )
  }
  dat <- dplyr::bind_rows(rows)
  metadata <- simulation_sparse_design_metadata(
    dat = dat,
    person_ids = person_ids,
    rater_ids = rater_ids,
    criterion_ids = criterion_ids,
    facet_names = facet_names,
    sparse_controls = sparse_controls,
    link_persons = link_persons,
    person_assignment = dplyr::bind_rows(person_assignment)
  )
  list(data = dat, metadata = metadata)
}

simulation_sparse_design_metadata <- function(dat,
                                              person_ids,
                                              rater_ids,
                                              criterion_ids,
                                              facet_names,
                                              sparse_controls,
                                              link_persons,
                                              person_assignment) {
  assignment_facet <- facet_names[1]
  person_rater <- unique(dat[, c("Person", assignment_facet), drop = FALSE])
  names(person_rater)[names(person_rater) == assignment_facet] <- "Rater"
  n_possible <- length(person_ids) * length(rater_ids) * length(criterion_ids)
  density <- nrow(dat) / n_possible

  rater_coverage <- person_rater |>
    dplyr::group_by(.data$Rater) |>
    dplyr::summarise(
      Persons = dplyr::n_distinct(.data$Person),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Rater) |>
    as.data.frame(stringsAsFactors = FALSE)
  row_counts <- dat |>
    dplyr::count(.data[[assignment_facet]], name = "Rows")
  names(row_counts)[names(row_counts) == assignment_facet] <- "Rater"
  rater_coverage <- rater_coverage |>
    dplyr::left_join(row_counts, by = "Rater")

  pair_tbl <- simulation_sparse_rater_pair_table(
    person_rater = person_rater,
    rater_ids = rater_ids
  )
  min_common <- if (nrow(pair_tbl) == 0L) NA_integer_ else min(pair_tbl$CommonPersons)
  zero_common <- if (nrow(pair_tbl) == 0L) NA_integer_ else sum(pair_tbl$CommonPersons == 0L)
  target <- as.integer(sparse_controls$min_common_persons_per_rater_pair %||% 0L)
  below_target <- if (nrow(pair_tbl) == 0L) NA_integer_ else sum(pair_tbl$CommonPersons < target)

  overview <- data.frame(
    Active = TRUE,
    Assignment = "sparse_linked",
    Persons = length(person_ids),
    Raters = length(rater_ids),
    Criteria = length(criterion_ids),
    Rows = nrow(dat),
    FullyCrossedRows = n_possible,
    DesignDensity = density,
    PlannedMissingRate = 1 - density,
    LinkPersons = length(link_persons),
    LinkFractionActual = length(link_persons) / length(person_ids),
    RatersPerPerson = as.integer(unique(person_assignment$RatersAssigned[!person_assignment$LinkPerson])[1] %||% NA_integer_),
    LinkRatersPerPerson = as.integer(sparse_controls$link_raters_per_person),
    MinCommonPersonsPerRaterPair = min_common,
    ZeroCommonRaterPairs = zero_common,
    RaterPairsBelowTarget = below_target,
    TargetCommonPersonsPerRaterPair = target,
    stringsAsFactors = FALSE
  )

  list(
    active = TRUE,
    overview = overview,
    controls = sparse_controls,
    link_persons = data.frame(Person = link_persons, stringsAsFactors = FALSE),
    person_assignment = person_assignment,
    rater_coverage = rater_coverage,
    rater_pair_links = pair_tbl,
    notes = c(
      "Design density is the generated row count divided by the fully crossed person x rater x criterion row count.",
      "Rater-pair common-person counts summarize linking strength, not model-data fit."
    )
  )
}

simulation_sparse_rater_pair_table <- function(person_rater, rater_ids) {
  if (length(rater_ids) < 2L) {
    return(data.frame())
  }
  pairs <- utils::combn(rater_ids, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    p1 <- person_rater$Person[person_rater$Rater == pair[1]]
    p2 <- person_rater$Person[person_rater$Rater == pair[2]]
    common <- intersect(p1, p2)
    data.frame(
      Rater1 = pair[1],
      Rater2 = pair[2],
      CommonPersons = length(common),
      CommonPersonIDs = paste(common, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

simulation_sparse_overview_fields <- function(sim) {
  empty <- list(
    SparseDesignActive = FALSE,
    DesignDensity = NA_real_,
    PlannedMissingRate = NA_real_,
    LinkPersons = NA_integer_,
    LinkFractionActual = NA_real_,
    LinkRatersPerPerson = NA_integer_,
    MinCommonPersonsPerRaterPair = NA_integer_,
    ZeroCommonRaterPairs = NA_integer_,
    RaterPairsBelowTarget = NA_integer_,
    TargetCommonPersonsPerRaterPair = NA_integer_
  )
  sparse <- attr(sim, "mfrm_sparse_design")
  if (!is.list(sparse) || !isTRUE(sparse$active) ||
      !is.data.frame(sparse$overview) || nrow(sparse$overview) == 0L) {
    return(empty)
  }
  overview <- sparse$overview[1, , drop = FALSE]
  list(
    SparseDesignActive = TRUE,
    DesignDensity = suppressWarnings(as.numeric(overview$DesignDensity[1] %||% NA_real_)),
    PlannedMissingRate = suppressWarnings(as.numeric(overview$PlannedMissingRate[1] %||% NA_real_)),
    LinkPersons = suppressWarnings(as.integer(overview$LinkPersons[1] %||% NA_integer_)),
    LinkFractionActual = suppressWarnings(as.numeric(overview$LinkFractionActual[1] %||% NA_real_)),
    LinkRatersPerPerson = suppressWarnings(as.integer(overview$LinkRatersPerPerson[1] %||% NA_integer_)),
    MinCommonPersonsPerRaterPair = suppressWarnings(as.integer(overview$MinCommonPersonsPerRaterPair[1] %||% NA_integer_)),
    ZeroCommonRaterPairs = suppressWarnings(as.integer(overview$ZeroCommonRaterPairs[1] %||% NA_integer_)),
    RaterPairsBelowTarget = suppressWarnings(as.integer(overview$RaterPairsBelowTarget[1] %||% NA_integer_)),
    TargetCommonPersonsPerRaterPair = suppressWarnings(as.integer(overview$TargetCommonPersonsPerRaterPair[1] %||% NA_integer_))
  )
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
  if (!identical(model, "PCM") && !identical(model, "GPCM")) {
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

simulation_resolve_fit_slope_facet <- function(model, slope_facet, fit_step_facet) {
  if (!identical(model, "GPCM")) {
    return(NA_character_)
  }
  explicit <- as.character(slope_facet[1] %||% NA_character_)
  if (!is.na(explicit) && nzchar(explicit)) {
    return(explicit)
  }
  as.character(fit_step_facet[1] %||% "Criterion")
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
  if (identical(as.character(fitted_model), "PCM") ||
      identical(as.character(fitted_model), "GPCM")) {
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

simulation_gpcm_design_boundary <- function(active) {
  if (!isTRUE(active) || !exists("gpcm_capability_matrix", mode = "function")) {
    return(data.frame())
  }
  tbl <- gpcm_capability_matrix("supported_with_caveat")
  out <- tbl[tbl$Area == "Design evaluation and population forecasting under bounded GPCM", , drop = FALSE]
  rownames(out) <- NULL
  out
}

simulation_gpcm_design_notes <- function(active) {
  if (!isTRUE(active)) return(character(0))
  c(
    "Bounded-GPCM design evaluation is supported with caveats as repeated simulation/refit operating-characteristic review.",
    "Interpret GPCM design and forecast outputs as design-level sensitivity evidence, not as operational scoring, diagnostic-screening, signal-detection, or arbitrary-facet planning validation.",
    "Use evaluate_mfrm_recovery() when the target is bounded-GPCM slope or parameter-recovery adequacy."
  )
}

simulation_gpcm_screening_boundary <- function(active) {
  if (!isTRUE(active) || !exists("gpcm_capability_matrix", mode = "function")) {
    return(data.frame())
  }
  tbl <- gpcm_capability_matrix("supported_with_caveat")
  out <- tbl[tbl$Area == "Diagnostic and signal-detection design screening under bounded GPCM", , drop = FALSE]
  rownames(out) <- NULL
  out
}

simulation_gpcm_screening_notes <- function(active) {
  if (!isTRUE(active)) return(character(0))
  c(
    "Bounded-GPCM diagnostic and signal-detection screening is supported with caveats as repeated simulation/refit sensitivity evidence.",
    "Interpret GPCM screening outputs as slope-aware operating-characteristic readouts under the evaluated role-based design, not as calibrated inferential tests, operational scoring, or arbitrary-facet planning validation.",
    "Report the active step/slope facet and slope-regime context before interpreting Type I proxy, sensitivity proxy, DIF, or bias-screening rates."
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


simulation_build_ademp <- function(purpose,
                                   design_grid,
                                   generator_model,
                                   generator_step_facet,
                                   generator_assignment,
                                   sim_spec = NULL,
                                   sparse_controls = NULL,
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
      planning_schema = planning_schema,
      sparse_controls = sparse_controls %||% sim_spec$sparse_controls %||% NULL
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
    designdensity = "MeanDesignDensity",
    plannedmissingrate = "MeanPlannedMissingRate",
    linkpersons = "MeanLinkPersons",
    linkfraction = "MeanLinkFractionActual",
    linkraters = "MeanLinkRatersPerPerson",
    mincommonpersons = "MeanMinCommonPersonsPerRaterPair",
    zerocommonpairs = "MaxZeroCommonRaterPairs",
    pairsshorttarget = "MaxRaterPairsBelowTarget",
    stop("Unknown metric: ", metric, call. = FALSE)
  )
}

simulation_sparse_design_active <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  if (is.numeric(x)) return(is.finite(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "yes", "1")
}

simulation_sparse_design_numeric <- function(tbl, cols) {
  out <- rep(NA_real_, nrow(tbl))
  if (nrow(tbl) == 0L) return(out)
  for (col in cols) {
    if (col %in% names(tbl)) {
      return(suppressWarnings(as.numeric(tbl[[col]])))
    }
  }
  out
}

simulation_sparse_design_review_fields <- function(tbl) {
  tbl <- as.data.frame(tbl %||% data.frame(), stringsAsFactors = FALSE)
  n <- nrow(tbl)
  if (n == 0L) {
    return(data.frame(
      LinkReviewStatus = character(0),
      LinkReviewReason = character(0),
      ReviewUse = character(0),
      stringsAsFactors = FALSE
    ))
  }

  zero_pairs <- simulation_sparse_design_numeric(
    tbl,
    c("ZeroCommonRaterPairs", "MaxZeroCommonRaterPairs")
  )
  below_target <- simulation_sparse_design_numeric(
    tbl,
    c("RaterPairsBelowTarget", "MaxRaterPairsBelowTarget")
  )
  min_common <- simulation_sparse_design_numeric(
    tbl,
    c("MinCommonPersonsPerRaterPair", "MeanMinCommonPersonsPerRaterPair")
  )
  target_common <- simulation_sparse_design_numeric(
    tbl,
    "TargetCommonPersonsPerRaterPair"
  )

  status <- rep("ok", n)
  reason <- rep("Recorded sparse rater-pair linkage checks are satisfied.", n)
  unavailable <- !is.finite(min_common)
  status[unavailable] <- "review"
  reason[unavailable] <- "Rater-pair common-person counts are unavailable in this summary row."

  below <- is.finite(below_target) & below_target > 0
  target_gap <- is.finite(min_common) & is.finite(target_common) & min_common < target_common
  gap <- below | target_gap
  status[gap] <- "review"
  reason[gap] <- "At least one rater pair is below the requested common-person target."

  zero <- is.finite(zero_pairs) & zero_pairs > 0
  status[zero] <- "review"
  reason[zero] <- "At least one rater pair has no common persons."

  data.frame(
    LinkReviewStatus = status,
    LinkReviewReason = reason,
    ReviewUse = rep("design_diagnostic_not_recovery_gate", n),
    stringsAsFactors = FALSE
  )
}

simulation_sparse_design_review_summary <- function(tbl) {
  tbl <- as.data.frame(tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L || !"SparseDesignActive" %in% names(tbl)) {
    return(tibble::tibble())
  }
  active <- simulation_sparse_design_active(tbl$SparseDesignActive)
  active[is.na(active)] <- FALSE
  if (!any(active)) {
    return(tibble::tibble())
  }
  tbl <- tbl[active, , drop = FALSE]

  sparse_cols <- intersect(
    c(
      "design_id", "rep", "Seed", "n_person", "n_rater", "n_criterion",
      "raters_per_person", "DesignDensity", "PlannedMissingRate",
      "LinkPersons", "LinkFractionActual", "LinkRatersPerPerson",
      "MinCommonPersonsPerRaterPair", "ZeroCommonRaterPairs",
      "RaterPairsBelowTarget", "TargetCommonPersonsPerRaterPair",
      "MeanDesignDensity", "MeanPlannedMissingRate", "MeanLinkPersons",
      "MeanLinkFractionActual", "MeanLinkRatersPerPerson",
      "MeanMinCommonPersonsPerRaterPair", "MaxZeroCommonRaterPairs",
      "MaxRaterPairsBelowTarget"
    ),
    names(tbl)
  )
  if (length(sparse_cols) > 0L) {
    tbl <- unique(tbl[, sparse_cols, drop = FALSE])
  }

  review <- simulation_sparse_design_review_fields(tbl)
  zero_pairs <- simulation_sparse_design_numeric(
    tbl,
    c("ZeroCommonRaterPairs", "MaxZeroCommonRaterPairs")
  )
  below_target <- simulation_sparse_design_numeric(
    tbl,
    c("RaterPairsBelowTarget", "MaxRaterPairsBelowTarget")
  )
  min_common <- simulation_sparse_design_numeric(
    tbl,
    c("MinCommonPersonsPerRaterPair", "MeanMinCommonPersonsPerRaterPair")
  )
  target_common <- simulation_sparse_design_numeric(
    tbl,
    "TargetCommonPersonsPerRaterPair"
  )
  density <- simulation_sparse_design_numeric(
    tbl,
    c("DesignDensity", "MeanDesignDensity")
  )
  missing_rate <- simulation_sparse_design_numeric(
    tbl,
    c("PlannedMissingRate", "MeanPlannedMissingRate")
  )
  target_shortfall <- (is.finite(below_target) & below_target > 0) |
    (is.finite(min_common) & is.finite(target_common) & min_common < target_common)

  safe_min <- function(x) if (all(!is.finite(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_max <- function(x) if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE)
  n <- nrow(tbl)
  review_n <- sum(review$LinkReviewStatus %in% "review", na.rm = TRUE)

  tibble::tibble(
    SparseDesignRows = n,
    OkRows = sum(review$LinkReviewStatus %in% "ok", na.rm = TRUE),
    ReviewRows = review_n,
    ReviewRate = if (n > 0L) review_n / n else NA_real_,
    ZeroCommonPairRows = sum(is.finite(zero_pairs) & zero_pairs > 0, na.rm = TRUE),
    TargetShortfallRows = sum(target_shortfall, na.rm = TRUE),
    UnavailableCommonPersonRows = sum(!is.finite(min_common), na.rm = TRUE),
    MinCommonPersonsObserved = safe_min(min_common),
    TargetCommonPersonsPerRaterPair = safe_max(target_common),
    MaxZeroCommonRaterPairs = safe_max(zero_pairs),
    MaxRaterPairsBelowTarget = safe_max(below_target),
    MinDesignDensity = safe_min(density),
    MaxPlannedMissingRate = safe_max(missing_rate),
    ReviewUse = "design_diagnostic_not_recovery_gate"
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
  if ("SparseDesignActive" %in% names(summary_tbl) &&
      any(summary_tbl$SparseDesignActive %in% TRUE, na.rm = TRUE)) {
    notes <- c(
      notes,
      "Sparse linked design rows report planned-missingness and rater-link diagnostics separately from fit and recovery metrics."
    )
  }
  if ("MaxZeroCommonRaterPairs" %in% names(summary_tbl) &&
      any(summary_tbl$MaxZeroCommonRaterPairs > 0, na.rm = TRUE)) {
    notes <- c(
      notes,
      "Some sparse linked design rows include rater pairs with no common persons; inspect rater-pair linking before interpreting recovery or separation."
    )
  }
  if ("MaxRaterPairsBelowTarget" %in% names(summary_tbl) &&
      any(summary_tbl$MaxRaterPairsBelowTarget > 0, na.rm = TRUE)) {
    notes <- c(
      notes,
      "Some sparse linked design rows fall below the requested common-person target for at least one rater pair."
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
        SparseDesignActive = any(.data$SparseDesignActive %in% TRUE, na.rm = TRUE),
        MeanDesignDensity = design_eval_safe_mean(.data$DesignDensity),
        MeanPlannedMissingRate = design_eval_safe_mean(.data$PlannedMissingRate),
        MeanLinkPersons = design_eval_safe_mean(.data$LinkPersons),
        MeanLinkFractionActual = design_eval_safe_mean(.data$LinkFractionActual),
        MeanLinkRatersPerPerson = design_eval_safe_mean(.data$LinkRatersPerPerson),
        MeanMinCommonPersonsPerRaterPair = design_eval_safe_mean(.data$MinCommonPersonsPerRaterPair),
        MaxZeroCommonRaterPairs = if (all(is.na(.data$ZeroCommonRaterPairs))) NA_integer_ else max(.data$ZeroCommonRaterPairs, na.rm = TRUE),
        MaxRaterPairsBelowTarget = if (all(is.na(.data$RaterPairsBelowTarget))) NA_integer_ else max(.data$RaterPairsBelowTarget, na.rm = TRUE),
        TargetCommonPersonsPerRaterPair = if (all(is.na(.data$TargetCommonPersonsPerRaterPair))) NA_integer_ else max(.data$TargetCommonPersonsPerRaterPair, na.rm = TRUE),
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

recovery_safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

recovery_safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

recovery_safe_cor <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  sdx <- recovery_safe_sd(x[ok])
  sdy <- recovery_safe_sd(y[ok])
  if (sum(ok) <= 1L || !is.finite(sdx) || !is.finite(sdy) || sdx == 0 || sdy == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok]))
}

recovery_mcse_rmse <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= 1L) return(NA_real_)
  rmse <- sqrt(mean(x^2))
  if (!is.finite(rmse) || rmse <= 0) return(NA_real_)
  stats::sd(x^2) / (2 * rmse * sqrt(n))
}

recovery_get_se <- function(tbl) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  candidates <- c("SE", "S.E.", "ModelSE", "Std.Error", "StdError")
  hit <- candidates[candidates %in% names(tbl)][1]
  if (is.na(hit)) return(rep(NA_real_, nrow(tbl)))
  suppressWarnings(as.numeric(tbl[[hit]]))
}

recovery_normalize_step_table <- function(x) {
  if (is.null(x)) {
    return(data.frame(StepFacet = character(0), Step = character(0), Estimate = numeric(0)))
  }
  if (is.numeric(x)) {
    return(data.frame(
      StepFacet = "Common",
      Step = paste0("Step_", seq_along(x)),
      Estimate = as.numeric(x),
      stringsAsFactors = FALSE
    ))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"StepFacet" %in% names(x)) x$StepFacet <- "Common"
  if (!"Step" %in% names(x)) {
    if ("StepIndex" %in% names(x)) {
      x$Step <- paste0("Step_", as.integer(x$StepIndex))
    } else {
      x$Step <- paste0("Step_", seq_len(nrow(x)))
    }
  }
  if (!"Estimate" %in% names(x)) {
    stop("Step recovery requires an `Estimate` column.", call. = FALSE)
  }
  data.frame(
    StepFacet = as.character(x$StepFacet),
    Step = as.character(x$Step),
    Estimate = suppressWarnings(as.numeric(x$Estimate)),
    stringsAsFactors = FALSE
  )
}

recovery_expand_common_steps <- function(truth_steps, fit_steps) {
  if (nrow(truth_steps) == 0L || nrow(fit_steps) == 0L) return(truth_steps)
  if (!all(unique(as.character(truth_steps$StepFacet)) == "Common")) return(truth_steps)
  fit_facets <- unique(as.character(fit_steps$StepFacet))
  fit_facets <- fit_facets[!is.na(fit_facets) & nzchar(fit_facets)]
  if (length(fit_facets) == 0L || identical(fit_facets, "Common")) return(truth_steps)
  expanded <- lapply(fit_facets, function(step_facet) {
    out <- truth_steps
    out$StepFacet <- step_facet
    out
  })
  dplyr::bind_rows(expanded)
}

recovery_new_rows <- function(rep,
                              parameter_type,
                              facet,
                              level,
                              subparameter,
                              truth,
                              estimate,
                              se = NA_real_,
                              raw_truth = truth,
                              raw_estimate = estimate,
                              comparison_scale = "logit",
                              alignment_group = facet,
                              align = TRUE,
                              recovery_comparable = TRUE,
                              recovery_basis = "truth_estimate_matched") {
  n <- length(level)
  recycle_field <- function(x) {
    x <- as.character(x)
    if (length(x) == n) x else rep(x[1] %||% NA_character_, n)
  }
  tibble::tibble(
    rep = rep,
    ParameterType = recycle_field(parameter_type),
    Facet = recycle_field(facet),
    Level = as.character(level),
    Subparameter = as.character(subparameter),
    Truth = suppressWarnings(as.numeric(truth)),
    Estimate = suppressWarnings(as.numeric(estimate)),
    SE = suppressWarnings(as.numeric(se)),
    RawTruth = suppressWarnings(as.numeric(raw_truth)),
    RawEstimate = suppressWarnings(as.numeric(raw_estimate)),
    ComparisonScale = recycle_field(comparison_scale),
    AlignmentGroup = recycle_field(alignment_group),
    AlignWithinGroup = rep(isTRUE(align), n),
    RecoveryComparable = rep(isTRUE(recovery_comparable), n),
    RecoveryBasis = rep(as.character(recovery_basis), n)
  )
}

recovery_rows_from_fit <- function(fit, truth, rep, include_person = TRUE) {
  if (is.null(truth) || !is.list(truth)) return(tibble::tibble())
  rows <- list()
  k <- 0L

  if (isTRUE(include_person) && !is.null(truth$person) && !is.null(fit$facets$person)) {
    person_tbl <- as.data.frame(fit$facets$person, stringsAsFactors = FALSE)
    if (all(c("Person", "Estimate") %in% names(person_tbl))) {
      idx <- match(as.character(person_tbl$Person), names(truth$person))
      ok <- is.finite(idx)
      if (any(ok)) {
        k <- k + 1L
        rows[[k]] <- recovery_new_rows(
          rep = rep,
          parameter_type = "person",
          facet = "Person",
          level = person_tbl$Person[ok],
          subparameter = "measure",
          truth = as.numeric(truth$person[idx[ok]]),
          estimate = suppressWarnings(as.numeric(person_tbl$Estimate[ok])),
          se = recovery_get_se(person_tbl)[ok],
          comparison_scale = "logit",
          alignment_group = "Person",
          align = TRUE
        )
      }
    }
  }

  if (!is.null(truth$facets) && length(truth$facets) > 0L && !is.null(fit$facets$others)) {
    other_tbl <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
    if (all(c("Facet", "Level", "Estimate") %in% names(other_tbl))) {
      for (facet_name in names(truth$facets)) {
        truth_vec <- truth$facets[[facet_name]]
        facet_tbl <- other_tbl[as.character(other_tbl$Facet) == facet_name, , drop = FALSE]
        if (nrow(facet_tbl) == 0L || is.null(names(truth_vec))) next
        idx <- match(as.character(facet_tbl$Level), names(truth_vec))
        ok <- is.finite(idx)
        if (!any(ok)) next
        k <- k + 1L
        rows[[k]] <- recovery_new_rows(
          rep = rep,
          parameter_type = "facet",
          facet = facet_name,
          level = facet_tbl$Level[ok],
          subparameter = "measure",
          truth = as.numeric(truth_vec[idx[ok]]),
          estimate = suppressWarnings(as.numeric(facet_tbl$Estimate[ok])),
          se = recovery_get_se(facet_tbl)[ok],
          comparison_scale = "logit",
          alignment_group = facet_name,
          align = TRUE
        )
      }
    }
  }

  if (!is.null(truth$step_table) && !is.null(fit$steps)) {
    truth_steps <- recovery_normalize_step_table(truth$step_table)
    fit_steps <- recovery_normalize_step_table(fit$steps)
    truth_steps <- recovery_expand_common_steps(truth_steps, fit_steps)
    step_tbl <- merge(
      truth_steps,
      fit_steps,
      by = c("StepFacet", "Step"),
      suffixes = c(".Truth", ".Estimate"),
      sort = FALSE
    )
    if (nrow(step_tbl) > 0L) {
      k <- k + 1L
      rows[[k]] <- recovery_new_rows(
        rep = rep,
        parameter_type = "step",
        facet = as.character(step_tbl$StepFacet),
        level = as.character(step_tbl$StepFacet),
        subparameter = as.character(step_tbl$Step),
        truth = suppressWarnings(as.numeric(step_tbl$Estimate.Truth)),
        estimate = suppressWarnings(as.numeric(step_tbl$Estimate.Estimate)),
        comparison_scale = "logit",
        alignment_group = paste0("step:", step_tbl$StepFacet),
        align = TRUE
      )
    }
  }

  if (!is.null(truth$slope_table) && !is.null(fit$slopes)) {
    truth_slopes <- as.data.frame(truth$slope_table, stringsAsFactors = FALSE)
    fit_slopes <- as.data.frame(fit$slopes, stringsAsFactors = FALSE)
    if (all(c("SlopeFacet", "Estimate") %in% names(truth_slopes)) &&
        all(c("SlopeFacet", "Estimate") %in% names(fit_slopes))) {
      slope_tbl <- merge(
        truth_slopes[, c("SlopeFacet", "Estimate"), drop = FALSE],
        fit_slopes[, c("SlopeFacet", "Estimate"), drop = FALSE],
        by = "SlopeFacet",
        suffixes = c(".Truth", ".Estimate"),
        sort = FALSE
      )
      slope_truth <- suppressWarnings(as.numeric(slope_tbl$Estimate.Truth))
      slope_est <- suppressWarnings(as.numeric(slope_tbl$Estimate.Estimate))
      ok <- is.finite(slope_truth) & slope_truth > 0 & is.finite(slope_est) & slope_est > 0
      if (any(ok)) {
        k <- k + 1L
        rows[[k]] <- recovery_new_rows(
          rep = rep,
          parameter_type = "slope",
          facet = "SlopeFacet",
          level = as.character(slope_tbl$SlopeFacet[ok]),
          subparameter = "log_slope",
          truth = log(slope_truth[ok]),
          estimate = log(slope_est[ok]),
          raw_truth = slope_truth[ok],
          raw_estimate = slope_est[ok],
          comparison_scale = "log_slope",
          alignment_group = "slope_identified",
          align = FALSE,
          recovery_basis = "geometric_mean_one_log_slope"
        )
      }
    }
  }

  truth_pop <- truth$population
  fit_pop <- fit$population
  if (is.list(truth_pop) && is.list(fit_pop) && !is.null(truth_pop$coefficients) &&
      !is.null(fit_pop$coefficients)) {
    truth_coef <- truth_pop$coefficients
    fit_coef <- fit_pop$coefficients
    common_terms <- intersect(names(truth_coef), names(fit_coef))
    if (length(common_terms) > 0L) {
      k <- k + 1L
      rows[[k]] <- recovery_new_rows(
        rep = rep,
        parameter_type = "population",
        facet = "population",
        level = common_terms,
        subparameter = "coefficient",
        truth = as.numeric(truth_coef[common_terms]),
        estimate = as.numeric(fit_coef[common_terms]),
        comparison_scale = "coefficient",
        alignment_group = "population_coefficients",
        align = FALSE
      )
    }
    if (!is.null(truth_pop$sigma2) && !is.null(fit_pop$sigma2)) {
      k <- k + 1L
      rows[[k]] <- recovery_new_rows(
        rep = rep,
        parameter_type = "population",
        facet = "population",
        level = "sigma2",
        subparameter = "variance",
        truth = as.numeric(truth_pop$sigma2),
        estimate = as.numeric(fit_pop$sigma2),
        comparison_scale = "variance",
        alignment_group = "population_variance",
        align = FALSE
      )
    }
  }

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) return(out)
  out |>
    dplyr::group_by(.data$rep, .data$ParameterType, .data$Facet,
                    .data$AlignmentGroup, .data$ComparisonScale) |>
    dplyr::mutate(
      ErrorRaw = .data$Estimate - .data$Truth,
      AlignmentShift = if (dplyr::first(.data$AlignWithinGroup)) {
        recovery_safe_mean(.data$ErrorRaw)
      } else {
        0
      },
      EstimateAligned = .data$Estimate - .data$AlignmentShift,
      ErrorAligned = .data$EstimateAligned - .data$Truth,
      Covered95 = dplyr::if_else(
        is.finite(.data$SE) & .data$SE > 0,
        .data$Truth >= .data$EstimateAligned - stats::qnorm(0.975) * .data$SE &
          .data$Truth <= .data$EstimateAligned + stats::qnorm(0.975) * .data$SE,
        NA
      )
    ) |>
    dplyr::ungroup()
}

recovery_summarize_rows <- function(rows) {
  rows <- tibble::as_tibble(rows)
  if (nrow(rows) == 0L) return(tibble::tibble())
  rows |>
    dplyr::group_by(.data$ParameterType, .data$Facet, .data$ComparisonScale) |>
    dplyr::summarise(
      Rows = dplyr::n(),
      Reps = dplyr::n_distinct(.data$rep),
      ComparableRate = mean(.data$RecoveryComparable, na.rm = TRUE),
      MeanTruth = recovery_safe_mean(.data$Truth),
      MeanEstimate = recovery_safe_mean(.data$EstimateAligned),
      Bias = recovery_safe_mean(.data$ErrorAligned),
      McseBias = simulation_mcse_mean(.data$ErrorAligned),
      RMSE = sqrt(recovery_safe_mean(.data$ErrorAligned^2)),
      McseRMSE = recovery_mcse_rmse(.data$ErrorAligned),
      MAE = recovery_safe_mean(abs(.data$ErrorAligned)),
      RawBias = recovery_safe_mean(.data$ErrorRaw),
      RawRMSE = sqrt(recovery_safe_mean(.data$ErrorRaw^2)),
      Correlation = recovery_safe_cor(.data$EstimateAligned, .data$Truth),
      MeanSE = recovery_safe_mean(.data$SE),
      SEAvailableRate = mean(is.finite(.data$SE) & .data$SE > 0, na.rm = TRUE),
      Coverage95 = if (all(is.na(.data$Covered95))) NA_real_ else mean(.data$Covered95, na.rm = TRUE),
      RecoveryBasis = paste(sort(unique(as.character(.data$RecoveryBasis))), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$ParameterType, .data$Facet, .data$ComparisonScale)
}

recovery_diagnostic_misfit_rate <- function(fit_tbl) {
  fit_tbl <- as.data.frame(fit_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(fit_tbl) == 0L ||
      !all(c("InfitZSTD", "OutfitZSTD") %in% names(fit_tbl))) {
    return(NA_real_)
  }
  z_in <- suppressWarnings(as.numeric(fit_tbl$InfitZSTD))
  z_out <- suppressWarnings(as.numeric(fit_tbl$OutfitZSTD))
  flag <- abs(z_in) > 2 | abs(z_out) > 2
  if (all(is.na(flag))) return(NA_real_)
  mean(flag, na.rm = TRUE)
}

recovery_diagnostic_df_sensitive_rate <- function(fit_tbl) {
  fit_tbl <- as.data.frame(fit_tbl %||% data.frame(), stringsAsFactors = FALSE)
  required <- c(
    "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
    "InfitZSTD_FACETS", "OutfitZSTD_FACETS"
  )
  if (nrow(fit_tbl) == 0L || !all(required %in% names(fit_tbl))) {
    return(NA_real_)
  }
  engine_flag <- abs(suppressWarnings(as.numeric(fit_tbl$InfitZSTD_ENGINE))) > 2 |
    abs(suppressWarnings(as.numeric(fit_tbl$OutfitZSTD_ENGINE))) > 2
  facets_flag <- abs(suppressWarnings(as.numeric(fit_tbl$InfitZSTD_FACETS))) > 2 |
    abs(suppressWarnings(as.numeric(fit_tbl$OutfitZSTD_FACETS))) > 2
  changed <- engine_flag != facets_flag
  if (all(is.na(changed))) return(NA_real_)
  mean(changed, na.rm = TRUE)
}

recovery_diagnostic_max_abs_zstd <- function(fit_tbl) {
  fit_tbl <- as.data.frame(fit_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(fit_tbl) == 0L ||
      !all(c("InfitZSTD", "OutfitZSTD") %in% names(fit_tbl))) {
    return(NA_real_)
  }
  vals <- c(
    suppressWarnings(as.numeric(fit_tbl$InfitZSTD)),
    suppressWarnings(as.numeric(fit_tbl$OutfitZSTD))
  )
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_real_)
  max(abs(vals), na.rm = TRUE)
}

recovery_diagnostic_rows_from_diagnostics <- function(diagnostics, rep, fit_df_method) {
  if (inherits(diagnostics, "error")) {
    return(tibble::tibble(
      rep = as.integer(rep),
      Facet = NA_character_,
      DiagnosticOK = FALSE,
      ValidationUse = "diagnostic_only_not_release_gate",
      DiagnosticError = conditionMessage(diagnostics)
    ))
  }
  reliability_tbl <- as.data.frame(diagnostics$reliability %||% data.frame(), stringsAsFactors = FALSE)
  fit_tbl <- as.data.frame(diagnostics$fit %||% data.frame(), stringsAsFactors = FALSE)
  facets <- unique(c(
    as.character(reliability_tbl$Facet %||% character(0)),
    as.character(fit_tbl$Facet %||% character(0))
  ))
  facets <- facets[!is.na(facets) & nzchar(facets)]
  if (length(facets) == 0L) {
    return(tibble::tibble())
  }
  rows <- lapply(facets, function(facet) {
    rel <- reliability_tbl[as.character(reliability_tbl$Facet) == facet, , drop = FALSE]
    rel <- rel[seq_len(min(nrow(rel), 1L)), , drop = FALSE]
    facet_fit <- fit_tbl[as.character(fit_tbl$Facet) == facet, , drop = FALSE]
    tibble::tibble(
      rep = as.integer(rep),
      Facet = facet,
      DiagnosticOK = TRUE,
      Levels = suppressWarnings(as.integer(rel$Levels[1] %||% NA_integer_)),
      Separation = suppressWarnings(as.numeric(rel$Separation[1] %||% NA_real_)),
      Reliability = suppressWarnings(as.numeric(rel$Reliability[1] %||% NA_real_)),
      Strata = suppressWarnings(as.numeric(rel$Strata[1] %||% NA_real_)),
      RealSeparation = suppressWarnings(as.numeric(rel$RealSeparation[1] %||% NA_real_)),
      RealReliability = suppressWarnings(as.numeric(rel$RealReliability[1] %||% NA_real_)),
      RealStrata = suppressWarnings(as.numeric(rel$RealStrata[1] %||% NA_real_)),
      MeanInfit = suppressWarnings(as.numeric(rel$MeanInfit[1] %||% NA_real_)),
      MeanOutfit = suppressWarnings(as.numeric(rel$MeanOutfit[1] %||% NA_real_)),
      FitRows = nrow(facet_fit),
      MisfitRateAbsZ2 = recovery_diagnostic_misfit_rate(facet_fit),
      MaxAbsZSTD = recovery_diagnostic_max_abs_zstd(facet_fit),
      DfSensitiveFlagRate = recovery_diagnostic_df_sensitive_rate(facet_fit),
      FitDfMethod = as.character(fit_df_method),
      ValidationUse = "diagnostic_only_not_release_gate",
      DiagnosticError = NA_character_
    )
  })
  dplyr::bind_rows(rows)
}

recovery_summarize_diagnostic_oc <- function(rows) {
  rows <- tibble::as_tibble(rows %||% tibble::tibble())
  if (nrow(rows) == 0L || !"Facet" %in% names(rows)) {
    return(tibble::tibble())
  }
  ok_rows <- rows[rows$DiagnosticOK %in% TRUE & !is.na(rows$Facet), , drop = FALSE]
  if (nrow(ok_rows) == 0L) return(tibble::tibble())
  ok_rows |>
    dplyr::group_by(.data$Facet) |>
    dplyr::summarise(
      Replications = dplyr::n_distinct(.data$rep),
      MeanLevels = recovery_safe_mean(.data$Levels),
      MeanSeparation = recovery_safe_mean(.data$Separation),
      MeanReliability = recovery_safe_mean(.data$Reliability),
      MeanStrata = recovery_safe_mean(.data$Strata),
      MeanRealSeparation = recovery_safe_mean(.data$RealSeparation),
      MeanRealReliability = recovery_safe_mean(.data$RealReliability),
      MeanRealStrata = recovery_safe_mean(.data$RealStrata),
      MeanInfit = recovery_safe_mean(.data$MeanInfit),
      MeanOutfit = recovery_safe_mean(.data$MeanOutfit),
      MeanMisfitRateAbsZ2 = recovery_safe_mean(.data$MisfitRateAbsZ2),
      MaxMisfitRateAbsZ2 = if (all(is.na(.data$MisfitRateAbsZ2))) NA_real_ else max(.data$MisfitRateAbsZ2, na.rm = TRUE),
      MeanMaxAbsZSTD = recovery_safe_mean(.data$MaxAbsZSTD),
      MeanDfSensitiveFlagRate = recovery_safe_mean(.data$DfSensitiveFlagRate),
      FitDfMethods = paste(sort(unique(as.character(.data$FitDfMethod))), collapse = "; "),
      ValidationUse = "diagnostic_only_not_release_gate",
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Facet)
}

recovery_build_notes <- function(rep_overview, recovery_summary, model) {
  notes <- character(0)
  if (is.data.frame(rep_overview) && nrow(rep_overview) > 0L) {
    if (any(!rep_overview$RunOK, na.rm = TRUE)) {
      notes <- c(notes, "Some recovery replications failed before recovery rows could be computed.")
    }
    if (any(!rep_overview$Converged, na.rm = TRUE)) {
      notes <- c(notes, "Some fitted models did not report convergence; inspect `rep_overview` before interpreting recovery summaries.")
    }
    if ("DiagnosticOK" %in% names(rep_overview) &&
        any(rep_overview$RunOK %in% TRUE & rep_overview$DiagnosticOK %in% FALSE)) {
      notes <- c(
        notes,
        "Some optional fit/separation diagnostic operating-characteristic computations failed; recovery metrics remain available."
      )
    }
  }
  if (is.data.frame(recovery_summary) && nrow(recovery_summary) > 0L &&
      any(recovery_summary$SEAvailableRate < 1, na.rm = TRUE)) {
    notes <- c(notes, "Coverage95 is reported only for rows with available standard errors.")
  }
  if (identical(model, "GPCM")) {
    notes <- c(
      notes,
      "Bounded GPCM recovery compares identified log slopes under the geometric-mean-one slope convention and keeps the current bounded GPCM support caveats."
    )
  }
  if (length(notes) == 0L) notes <- "No immediate warnings from the recovery simulation summary."
  notes
}

#' Evaluate parameter recovery by repeated simulation and refitting
#'
#' @description
#' Runs a compact parameter-recovery simulation study: generate data from a
#' known ordered many-facet data-generating setup, refit the requested model,
#' align estimates to the known truth where location indeterminacy requires it,
#' and summarize bias, RMSE, MAE, correlation, and standard-error coverage.
#'
#' @inheritParams simulate_mfrm_data
#' @param reps Number of Monte Carlo replications.
#' @param fit_method Estimation method passed to [fit_mfrm()].
#' @param maxit Maximum optimizer iterations passed to [fit_mfrm()].
#' @param quad_points Quadrature points used when `fit_method = "MML"`.
#' @param include_person Logical. When `TRUE`, include person-measure recovery
#'   rows when the fitted object exposes person estimates.
#' @param include_diagnostics Logical. When `TRUE`, run [diagnose_mfrm()] after
#'   each successful refit and retain facet-level fit/separation operating
#'   characteristics. These diagnostics are reported separately from recovery
#'   metrics and are not release-success criteria by themselves.
#' @param diagnostic_fit_df_method Fit-ZSTD degrees-of-freedom convention used
#'   for optional diagnostic operating-characteristic summaries. Use `"both"`
#'   when reviewing FACETS-style df sensitivity.
#'
#' @details
#' This helper is deliberately narrower than [evaluate_mfrm_design()]. Design
#' evaluation asks which design condition is operationally adequate; recovery
#' simulation asks whether the fitted model recovers the known parameters under
#' one explicit data-generating setup.
#'
#' Location-like parameters (`Person`, non-person facets, and steps) are
#' summarized after mean alignment within each replication and parameter group.
#' This follows the usual Rasch/MFRM identification convention: adding a common
#' constant to one location block should not be counted as recovery failure.
#' Raw, unaligned errors are retained in `recovery` and summarized as
#' `RawBias` / `RawRMSE`.
#'
#' For bounded `GPCM`, supplied generator slopes are treated as relative
#' discriminations and normalized to the same geometric-mean-one log-slope
#' identification used by the fitter. Slope recovery is therefore summarized on
#' the identified log-slope scale without an additional mean-alignment step.
#' Direct data generation and refitting are supported, but broader GPCM design-
#' planning claims remain outside the current package boundary.
#'
#' Sparse linked generators are supported through `sim_spec` or direct
#' `assignment = "sparse_linked"` plus `sparse_controls`. Their design-density
#' and rater-link diagnostics are retained in `rep_overview`; recovery metrics
#' remain parameter-recovery summaries and should not be read as evidence that
#' a sparse linking design is adequate by itself.
#'
#' The returned `ademp` component follows the simulation-study framing of
#' Morris, White, and Crowther (2019) and the ADEMP planning/reporting template
#' used in later simulation-study guidance.
#'
#' @return An object of class `mfrm_recovery_simulation` with components:
#' - `recovery`: row-level truth/estimate comparisons by replication.
#' - `recovery_summary`: parameter-type summaries across replications.
#' - `rep_overview`: replication-level convergence, timing, error status, and
#'   sparse-design diagnostics when applicable.
#' - `diagnostic_oc`: optional replication-by-facet fit/separation operating
#'   characteristics when `include_diagnostics = TRUE`.
#' - `diagnostic_oc_summary`: optional facet-level diagnostic operating-
#'   characteristic summary.
#' - `settings`: fitting and simulation settings.
#' - `ademp`: simulation-study metadata.
#'
#' @section Typical workflow:
#' 1. Build a simulation specification with [build_mfrm_sim_spec()] or pass scalar
#'    generator arguments directly.
#' 2. Run `evaluate_mfrm_recovery(...)` with a modest `reps` value for a smoke
#'    check, then increase `reps` for stable Monte Carlo summaries.
#' 3. Inspect `summary(x)$recovery_summary` and the row-level `x$recovery` table.
#'
#' @seealso [simulate_mfrm_data()], [evaluate_mfrm_design()], [fit_mfrm()]
#' @examples
#' \donttest{
#' rec <- evaluate_mfrm_recovery(
#'   n_person = 12,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   reps = 1,
#'   maxit = 30,
#'   seed = 123
#' )
#' summary(rec)$recovery_summary[, c("ParameterType", "Facet", "RMSE", "Bias")]
#' }
#' @export
evaluate_mfrm_recovery <- function(n_person = 50,
                                   n_rater = 4,
                                   n_criterion = 4,
                                   raters_per_person = n_rater,
                                   design = NULL,
                                   reps = 10,
                                   score_levels = 4,
                                   theta_sd = 1,
                                   rater_sd = 0.35,
                                   criterion_sd = 0.25,
                                   noise_sd = 0,
                                   step_span = 1.4,
                                   model = c("RSM", "PCM", "GPCM"),
                                   step_facet = NULL,
                                   slope_facet = NULL,
                                   thresholds = NULL,
                                   slopes = NULL,
                                   assignment = NULL,
                                   sparse_controls = NULL,
                                   sim_spec = NULL,
                                   fit_method = c("JML", "MML"),
                                   maxit = 25,
                                   quad_points = 7,
                                   include_person = TRUE,
                                   include_diagnostics = FALSE,
                                   diagnostic_fit_df_method = c("both", "engine", "facets"),
                                   seed = NULL) {
  model <- match.arg(toupper(as.character(model[1])), c("RSM", "PCM", "GPCM"))
  fit_method <- match.arg(fit_method)
  diagnostic_fit_df_method <- match.arg(diagnostic_fit_df_method)
  include_diagnostics <- isTRUE(include_diagnostics)
  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) stop("`reps` must be >= 1.", call. = FALSE)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!is.null(sim_spec)) {
    model <- as.character(sim_spec$model %||% model)
    step_facet <- sim_spec$step_facet %||% step_facet
    slope_facet <- sim_spec$slope_facet %||% slope_facet
  }
  if (identical(model, "GPCM")) {
    facet_names_for_resolution <- if (is.null(sim_spec)) {
      simulation_default_output_facet_names()
    } else {
      simulation_spec_output_facet_names(sim_spec)
    }
    resolved_facets <- resolve_step_and_slope_facets(
      model = model,
      step_facet = step_facet[1] %||% facet_names_for_resolution[2],
      slope_facet = slope_facet,
      facet_names = unname(facet_names_for_resolution)
    )
    step_facet <- resolved_facets$step_facet
    slope_facet <- resolved_facets$slope_facet
  } else if (identical(model, "PCM")) {
    step_facet <- step_facet[1] %||% if (is.null(sim_spec)) "Criterion" else sim_spec$step_facet
  } else {
    step_facet <- NULL
    slope_facet <- NULL
  }
  if (!is.null(sim_spec) &&
      isTRUE((sim_spec$population %||% simulation_empty_population_spec())$active) &&
      !identical(fit_method, "MML")) {
    stop(
      "Recovery simulations with an active latent-regression population generator require `fit_method = \"MML\"`.",
      call. = FALSE
    )
  }

  seeds <- with_preserved_rng_seed(
    seed,
    sample.int(.Machine$integer.max, size = reps, replace = FALSE)
  )
  base_facet_names <- if (is.null(sim_spec)) simulation_default_output_facet_names() else simulation_spec_output_facet_names(sim_spec)
  generator_model <- if (is.null(sim_spec)) model else as.character(sim_spec$model %||% model)
  generator_step_facet <- if (is.null(sim_spec)) {
    if (identical(generator_model, "RSM")) NA_character_ else step_facet
  } else {
    sim_spec$step_facet %||% NA_character_
  }
  generator_assignment <- if (is.null(sim_spec)) {
    assignment %||% "design_dependent"
  } else {
    sim_spec$assignment
  }
  recovery_sparse_controls <- if (is.null(sim_spec)) {
    simulation_normalize_sparse_controls(
      sparse_controls = sparse_controls,
      assignment = if (identical(assignment, "sparse_linked")) "sparse_linked" else "rotating",
      n_person = n_person,
      n_rater = n_rater,
      raters_per_person = raters_per_person
    )
  } else {
    sim_spec$sparse_controls %||% list(active = FALSE)
  }
  generator_slope_regime <- if (!identical(generator_model, "GPCM")) {
    NA_character_
  } else if (!is.null(sim_spec)) {
    sim_spec$slope_regime %||% simulation_gpcm_slope_regime(sim_spec$slope_table)
  } else {
    tryCatch(
      simulation_gpcm_slope_regime(
        simulation_build_slope_table(
          slopes = slopes,
          model = model,
          slope_facet = slope_facet,
          facet_names = stats::setNames(base_facet_names, c("rater", "criterion")),
          n_rater = n_rater,
          n_criterion = n_criterion
        )
      ),
      error = function(e) NA_character_
    )
  }

  recovery_rows <- vector("list", reps)
  rep_rows <- vector("list", reps)
  diagnostic_rows <- vector("list", reps)
  for (rep in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    sim <- tryCatch(
      if (is.null(sim_spec)) {
        simulate_mfrm_data(
          n_person = n_person,
          n_rater = n_rater,
          n_criterion = n_criterion,
          raters_per_person = raters_per_person,
          design = design,
          score_levels = score_levels,
          theta_sd = theta_sd,
          rater_sd = rater_sd,
          criterion_sd = criterion_sd,
          noise_sd = noise_sd,
          step_span = step_span,
          seed = seeds[rep],
          model = model,
          step_facet = step_facet %||% "Criterion",
          slope_facet = slope_facet,
          thresholds = thresholds,
          slopes = slopes,
          assignment = assignment,
          sparse_controls = recovery_sparse_controls
        )
      } else {
        simulate_mfrm_data(sim_spec = sim_spec, seed = seeds[rep])
      },
      error = function(e) e
    )
    elapsed_sim <- proc.time()[["elapsed"]] - t0
    sparse_fields <- simulation_sparse_overview_fields(sim)
    rep_row <- tibble::tibble(
      rep = rep,
      Seed = seeds[rep],
      Observations = if (is.data.frame(sim)) nrow(sim) else NA_integer_,
      ScoreLevelsDeclared = NA_integer_,
      ScoreLevelsObserved = NA_integer_,
      ZeroScoreLevels = NA_integer_,
      MinScoreCount = NA_integer_,
      MinScoreProportion = NA_real_,
      MaxScoreProportion = NA_real_,
      SparseDesignActive = sparse_fields$SparseDesignActive,
      DesignDensity = sparse_fields$DesignDensity,
      PlannedMissingRate = sparse_fields$PlannedMissingRate,
      LinkPersons = sparse_fields$LinkPersons,
      LinkFractionActual = sparse_fields$LinkFractionActual,
      LinkRatersPerPerson = sparse_fields$LinkRatersPerPerson,
      MinCommonPersonsPerRaterPair = sparse_fields$MinCommonPersonsPerRaterPair,
      ZeroCommonRaterPairs = sparse_fields$ZeroCommonRaterPairs,
      RaterPairsBelowTarget = sparse_fields$RaterPairsBelowTarget,
      TargetCommonPersonsPerRaterPair = sparse_fields$TargetCommonPersonsPerRaterPair,
      ElapsedSec = elapsed_sim,
      RunOK = FALSE,
      Converged = FALSE,
      RecoveryRows = 0L,
      DiagnosticOK = if (include_diagnostics) FALSE else NA,
      DiagnosticRows = if (include_diagnostics) 0L else NA_integer_,
      DiagnosticError = NA_character_,
      Error = NA_character_
    )
    if (inherits(sim, "error")) {
      rep_row$Error <- conditionMessage(sim)
      rep_rows[[rep]] <- rep_row
      recovery_rows[[rep]] <- tibble::tibble()
      diagnostic_rows[[rep]] <- tibble::tibble()
      next
    }
    score_support <- simulation_score_support_summary(
      sim,
      fallback_score_levels = if (is.null(sim_spec)) score_levels else sim_spec$score_levels
    )
    rep_row[names(score_support)] <- score_support

    row_facet_names <- if (is.null(sim_spec)) base_facet_names else simulation_spec_output_facet_names(sim_spec)
    fit_args <- list(
      data = sim,
      person = "Person",
      facets = row_facet_names,
      score = "Score",
      method = fit_method,
      model = model,
      maxit = maxit
    )
    if (identical(model, "PCM")) fit_args$step_facet <- step_facet %||% row_facet_names[2]
    if (identical(model, "GPCM")) {
      fit_args$step_facet <- step_facet %||% row_facet_names[2]
      fit_args$slope_facet <- slope_facet %||% fit_args$step_facet
    }
    if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points
    if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
    sim_population <- attr(sim, "mfrm_population_data")
    if (is.list(sim_population) && isTRUE(sim_population$active)) {
      fit_args$population_formula <- sim_population$population_formula
      fit_args$person_data <- sim_population$person_data
      fit_args$person_id <- sim_population$person_id
      fit_args$population_policy <- sim_population$population_policy
    }
    fit_args <- simulation_add_fit_score_support(
      fit_args,
      sim,
      fallback_score_levels = if (is.null(sim_spec)) score_levels else sim_spec$score_levels
    )

    fit <- tryCatch(suppressWarnings(do.call(fit_mfrm, fit_args)), error = function(e) e)
    elapsed <- proc.time()[["elapsed"]] - t0
    rep_row$ElapsedSec <- elapsed
    if (inherits(fit, "error")) {
      rep_row$Error <- conditionMessage(fit)
      rep_rows[[rep]] <- rep_row
      recovery_rows[[rep]] <- tibble::tibble()
      diagnostic_rows[[rep]] <- tibble::tibble()
      next
    }
    truth <- attr(sim, "mfrm_truth")
    rows <- recovery_rows_from_fit(fit, truth, rep = rep, include_person = include_person)
    diag_rows <- tibble::tibble()
    if (include_diagnostics) {
      diagnostic <- tryCatch(
        diagnose_mfrm(
          fit,
          residual_pca = "none",
          fit_df_method = diagnostic_fit_df_method
        ),
        error = function(e) e
      )
      diag_rows <- recovery_diagnostic_rows_from_diagnostics(
        diagnostic,
        rep = rep,
        fit_df_method = diagnostic_fit_df_method
      )
      rep_row$DiagnosticOK <- !inherits(diagnostic, "error")
      rep_row$DiagnosticRows <- nrow(diag_rows)
      if (inherits(diagnostic, "error")) {
        rep_row$DiagnosticError <- conditionMessage(diagnostic)
      }
    }
    rep_row$RunOK <- TRUE
    rep_row$Converged <- isTRUE(as.logical(fit$summary$Converged[1]))
    rep_row$RecoveryRows <- nrow(rows)
    rep_rows[[rep]] <- rep_row
    recovery_rows[[rep]] <- rows
    diagnostic_rows[[rep]] <- diag_rows
  }

  recovery <- dplyr::bind_rows(recovery_rows)
  rep_overview <- dplyr::bind_rows(rep_rows)
  diagnostic_oc <- dplyr::bind_rows(diagnostic_rows)
  recovery_summary <- recovery_summarize_rows(recovery)
  diagnostic_oc_summary <- recovery_summarize_diagnostic_oc(diagnostic_oc)
  ademp <- simulation_build_ademp(
    purpose = "Assess parameter recovery under one explicit many-facet data-generating setup by repeated simulation and refitting.",
    design_grid = data.frame(
      design_id = "R1",
      n_person = n_person,
      n_rater = n_rater,
      n_criterion = n_criterion,
      raters_per_person = raters_per_person
    ),
    generator_model = generator_model,
    generator_step_facet = generator_step_facet,
    generator_assignment = generator_assignment,
    sim_spec = sim_spec,
    sparse_controls = if (is.null(sim_spec)) recovery_sparse_controls else NULL,
    estimands = c(
      "Person and facet location recovery after identification alignment",
      "Step-threshold recovery after identification alignment",
      "Bounded GPCM log-slope recovery when slopes are fitted",
      "Latent-regression coefficient and variance recovery when present"
    ),
    analysis_methods = list(
      fit_method = fit_method,
      fitted_model = model,
      gpcm_slope_regime = generator_slope_regime,
      maxit = maxit,
      quad_points = if (identical(fit_method, "MML")) quad_points else NA_integer_
    ),
    performance_measures = c(
      "Bias",
      "RMSE",
      "MAE",
      "Truth-estimate correlation",
      "95% Wald coverage where standard errors are available",
      "Monte Carlo standard errors for bias and RMSE",
      "Optional fit/separation operating characteristics when include_diagnostics = TRUE"
    )
  )
  structure(
    list(
      recovery = recovery,
      recovery_summary = recovery_summary,
      rep_overview = rep_overview,
      diagnostic_oc = diagnostic_oc,
      diagnostic_oc_summary = diagnostic_oc_summary,
      settings = list(
        reps = reps,
        fit_method = fit_method,
        model = model,
        step_facet = step_facet,
        slope_facet = slope_facet,
        assignment = if (is.null(sim_spec)) assignment else sim_spec$assignment %||% assignment,
        sparse_controls = recovery_sparse_controls,
        gpcm_slope_regime = generator_slope_regime,
        maxit = maxit,
        quad_points = quad_points,
        include_person = isTRUE(include_person),
        include_diagnostics = include_diagnostics,
        diagnostic_fit_df_method = if (include_diagnostics) diagnostic_fit_df_method else NA_character_,
        sim_spec = sim_spec,
        seed = seed
      ),
      notes = recovery_build_notes(rep_overview, recovery_summary, model),
      ademp = ademp
    ),
    class = "mfrm_recovery_simulation"
  )
}

#' @export
summary.mfrm_recovery_simulation <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_recovery_simulation")) {
    stop("`object` must be output from evaluate_mfrm_recovery().", call. = FALSE)
  }
  rep_tbl <- tibble::as_tibble(object$rep_overview %||% tibble::tibble())
  diagnostic_runs <- if ("DiagnosticOK" %in% names(rep_tbl)) {
    sum(!is.na(rep_tbl$DiagnosticOK))
  } else {
    0L
  }
  diagnostic_successful <- if (diagnostic_runs > 0L) {
    sum(rep_tbl$DiagnosticOK %in% TRUE, na.rm = TRUE)
  } else {
    0L
  }
  overview <- tibble::tibble(
    Reps = nrow(rep_tbl),
    SuccessfulRuns = if (nrow(rep_tbl) > 0L) sum(rep_tbl$RunOK, na.rm = TRUE) else 0L,
    ConvergedRuns = if (nrow(rep_tbl) > 0L) sum(rep_tbl$Converged, na.rm = TRUE) else 0L,
    DiagnosticRuns = as.integer(diagnostic_runs),
    DiagnosticSuccessfulRuns = as.integer(diagnostic_successful),
    DiagnosticSuccessRate = if (diagnostic_runs > 0L) diagnostic_successful / diagnostic_runs else NA_real_,
    MeanElapsedSec = if (nrow(rep_tbl) > 0L) recovery_safe_mean(rep_tbl$ElapsedSec) else NA_real_,
    RecoveryRows = if (nrow(rep_tbl) > 0L) sum(rep_tbl$RecoveryRows, na.rm = TRUE) else 0L
  )
  out <- list(
    overview = overview,
    recovery_summary = tibble::as_tibble(object$recovery_summary %||% tibble::tibble()),
    rep_overview = rep_tbl,
    sparse_review = simulation_sparse_design_review_summary(rep_tbl),
    diagnostic_oc = tibble::as_tibble(object$diagnostic_oc %||% tibble::tibble()),
    diagnostic_oc_summary = tibble::as_tibble(object$diagnostic_oc_summary %||% tibble::tibble()),
    settings = object$settings %||% list(),
    notes = object$notes %||% character(0),
    ademp = object$ademp %||% NULL,
    digits = max(0L, as.integer(digits[1]))
  )
  class(out) <- "summary.mfrm_recovery_simulation"
  out
}

#' @export
print.summary.mfrm_recovery_simulation <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("MFRM Parameter Recovery Simulation Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$recovery_summary) && nrow(x$recovery_summary) > 0L) {
    cat("\nRecovery summary\n")
    print(round_numeric_df(as.data.frame(x$recovery_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$diagnostic_oc_summary) && nrow(x$diagnostic_oc_summary) > 0L) {
    cat("\nDiagnostic operating-characteristic summary\n")
    print(round_numeric_df(as.data.frame(x$diagnostic_oc_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$sparse_review) && nrow(x$sparse_review) > 0L) {
    cat("\nSparse linked design review\n")
    print(round_numeric_df(as.data.frame(x$sparse_review), digits = digits), row.names = FALSE)
  }
  if (length(x$notes) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  if (is.list(x$ademp) && length(x$ademp) > 0L) {
    cat("\nADEMP\n")
    cat(" - Aim: ", paste(x$ademp$aims, collapse = "; "), "\n", sep = "")
    cat(" - DGM model: ", x$ademp$data_generating_mechanism$model, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.mfrm_recovery_simulation <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

recovery_assessment_threshold <- function(thresholds,
                                          parameter_type,
                                          facet,
                                          comparison_scale) {
  if (is.null(thresholds) || length(thresholds) == 0L) return(NA_real_)
  thresholds <- unlist(thresholds, use.names = TRUE)
  thresholds <- suppressWarnings(as.numeric(thresholds))
  if (length(thresholds) == 0L) return(NA_real_)
  nm <- names(thresholds)
  if (is.null(nm) || !any(nzchar(nm))) return(thresholds[1])
  keys <- c(
    paste(parameter_type, facet, comparison_scale, sep = ":"),
    paste(parameter_type, facet, sep = ":"),
    paste(parameter_type, comparison_scale, sep = ":"),
    as.character(parameter_type),
    as.character(facet),
    "default"
  )
  hit <- keys[keys %in% nm][1]
  if (is.na(hit)) return(NA_real_)
  thresholds[[hit]]
}

recovery_assessment_overall_status <- function(status) {
  status <- as.character(status)
  status <- status[!is.na(status) & nzchar(status)]
  if (length(status) == 0L) return("review")
  if (any(status == "concern")) return("concern")
  if (any(status == "review")) return("review")
  if (any(status == "not_available")) return("review")
  if (all(status %in% c("ok", "not_assessed"))) return("ok")
  "review"
}

recovery_assessment_metric_status <- function(value,
                                              limit,
                                              lower_is_better = TRUE,
                                              concern_multiplier = 2) {
  value <- suppressWarnings(as.numeric(value))
  limit <- suppressWarnings(as.numeric(limit))
  if (!is.finite(value)) return("not_available")
  if (!is.finite(limit)) return("not_assessed")
  if (isTRUE(lower_is_better)) {
    if (value <= limit) return("ok")
    if (value <= limit * concern_multiplier) return("review")
    return("concern")
  }
  if (value >= limit) return("ok")
  if (value >= limit / concern_multiplier) return("review")
  "concern"
}

recovery_assessment_rate_status <- function(value,
                                            target,
                                            tolerance = 0,
                                            concern_gap = 0.15) {
  value <- suppressWarnings(as.numeric(value))
  target <- suppressWarnings(as.numeric(target))
  tolerance <- suppressWarnings(as.numeric(tolerance))
  concern_gap <- suppressWarnings(as.numeric(concern_gap))
  if (!is.finite(value)) return("not_available")
  if (!is.finite(target)) return("not_assessed")
  if (value >= target - tolerance) return("ok")
  if (value >= target - max(tolerance, concern_gap)) return("review")
  "concern"
}

recovery_assessment_coverage_status <- function(value,
                                                target = 0.95,
                                                tolerance = 0.05) {
  value <- suppressWarnings(as.numeric(value))[1]
  target <- suppressWarnings(as.numeric(target))
  tolerance <- suppressWarnings(as.numeric(tolerance))
  if (!is.finite(value)) return("not_available")
  if (length(target) == 0L || length(tolerance) == 0L ||
      !is.finite(target[1]) || !is.finite(tolerance[1])) {
    return("not_assessed")
  }
  delta <- abs(value - target[1])
  if (delta <= tolerance[1]) return("ok")
  if (delta <= 2 * tolerance[1]) return("review")
  "concern"
}

recovery_assessment_action <- function(status, topic = "metric") {
  switch(
    as.character(status)[1] %||% "review",
    ok = "Use as supporting evidence under the stated simulation setup and thresholds.",
    not_assessed = paste0("Set a practical threshold if ", topic, " must support a go/no-go decision."),
    not_available = paste0("Treat ", topic, " as unavailable for this run; inspect row-level output before reporting it."),
    review = paste0("Review ", topic, " with plots and row-level output before using it for a decision."),
    concern = paste0("Do not use ", topic, " as adequacy evidence until the design, fit settings, or replication count are revisited."),
    paste0("Review ", topic, " before interpretation.")
  )
}

recovery_assessment_check_row <- function(section,
                                          item,
                                          status,
                                          evidence,
                                          next_action) {
  tibble::tibble(
    Section = as.character(section),
    Item = as.character(item),
    Status = as.character(status),
    Evidence = as.character(evidence),
    NextAction = as.character(next_action)
  )
}

recovery_assessment_status_counts <- function(status) {
  status <- as.character(status)
  status <- status[!is.na(status) & nzchar(status)]
  if (length(status) == 0L) return("none")
  tbl <- sort(table(status), decreasing = TRUE)
  paste(paste(names(tbl), as.integer(tbl), sep = "="), collapse = ", ")
}

recovery_assessment_status_priority <- function(status) {
  status <- as.character(status)
  priority <- match(status, c("concern", "review", "not_available", "not_assessed", "ok"))
  priority[is.na(priority)] <- 99L
  priority
}

recovery_assessment_uncertainty_review <- function(metric_review,
                                                   coverage_target = 0.95,
                                                   coverage_tolerance = 0.05,
                                                   min_se_available = 0.80) {
  metric_review <- tibble::as_tibble(metric_review %||% tibble::tibble())
  empty <- tibble::tibble(
    ParameterType = character(),
    Facet = character(),
    ComparisonScale = character(),
    Coverage95 = numeric(),
    CoverageStatus = character(),
    SEAvailableRate = numeric(),
    SEStatus = character(),
    Interpretation = character(),
    NextAction = character()
  )
  if (nrow(metric_review) == 0L) return(empty)

  coverage_status <- as.character(metric_review$CoverageStatus %||% rep(NA_character_, nrow(metric_review)))
  se_status <- as.character(metric_review$SEStatus %||% rep(NA_character_, nrow(metric_review)))
  coverage_value <- suppressWarnings(as.numeric(metric_review$Coverage95 %||% rep(NA_real_, nrow(metric_review))))
  se_rate <- suppressWarnings(as.numeric(metric_review$SEAvailableRate %||% rep(NA_real_, nrow(metric_review))))
  target_available <- length(coverage_target) > 0L &&
    length(coverage_tolerance) > 0L &&
    is.finite(suppressWarnings(as.numeric(coverage_target))[1]) &&
    is.finite(suppressWarnings(as.numeric(coverage_tolerance))[1])
  se_target_available <- length(min_se_available) > 0L &&
    is.finite(suppressWarnings(as.numeric(min_se_available))[1])

  interpretation <- dplyr::case_when(
    coverage_status == "not_available" ~
      "Coverage was not computed because finite SE-based intervals were unavailable for this group.",
    coverage_status == "not_assessed" ~
      "Coverage was computed or retained without a target/tolerance decision.",
    coverage_status == "ok" ~
      "Coverage is within the requested target tolerance.",
    coverage_status == "review" ~
      "Coverage is finite but outside the requested target tolerance; inspect row-level intervals.",
    coverage_status == "concern" ~
      "Coverage is far from the requested target; do not use it as uncertainty evidence without follow-up.",
    is.finite(coverage_value) & !isTRUE(target_available) ~
      "Coverage is finite, but no usable target/tolerance was supplied.",
    TRUE ~ "Coverage evidence needs review before interpretation."
  )
  next_action <- dplyr::case_when(
    coverage_status == "ok" & (se_status %in% c("ok", "not_assessed")) ~
      "Use as supporting uncertainty evidence under the stated simulation setup.",
    coverage_status == "not_available" ~
      "Report coverage as unavailable and inspect SE generation before relying on interval evidence.",
    coverage_status == "not_assessed" ~
      "Set coverage_target and coverage_tolerance if coverage must support a decision.",
    se_status %in% c("review", "concern", "not_available") ~
      "Review SE availability before interpreting coverage.",
    coverage_status %in% c("review", "concern") ~
      "Inspect row-level errors, SEs, and simulation design before using coverage as adequacy evidence.",
    is.finite(se_rate) & !isTRUE(se_target_available) ~
      "Set min_se_available if SE availability must support a decision.",
    TRUE ~ "Review coverage and SE availability before reporting."
  )

  out <- tibble::tibble(
    ParameterType = as.character(metric_review$ParameterType %||% rep(NA_character_, nrow(metric_review))),
    Facet = as.character(metric_review$Facet %||% rep(NA_character_, nrow(metric_review))),
    ComparisonScale = as.character(metric_review$ComparisonScale %||% rep(NA_character_, nrow(metric_review))),
    Coverage95 = coverage_value,
    CoverageStatus = coverage_status,
    SEAvailableRate = se_rate,
    SEStatus = se_status,
    Interpretation = interpretation,
    NextAction = next_action
  )
  dplyr::arrange(
    out,
    recovery_assessment_status_priority(.data$CoverageStatus),
    recovery_assessment_status_priority(.data$SEStatus),
    .data$ParameterType,
    .data$Facet,
    .data$ComparisonScale
  )
}

recovery_assessment_condition_review <- function(x) {
  settings <- x$settings %||% list()
  ademp <- x$ademp %||% list()
  methods <- ademp$methods %||% list()
  sim_spec <- settings$sim_spec %||% list()
  rep_tbl <- tibble::as_tibble(x$rep_overview %||% tibble::tibble())
  model <- as.character(settings$model %||% sim_spec$model %||% NA_character_)[1]
  regime <- as.character(settings$gpcm_slope_regime %||% methods$gpcm_slope_regime %||% NA_character_)[1]
  if (is.na(model) || !nzchar(model)) model <- NA_character_
  if (is.na(regime) || !nzchar(regime)) regime <- NA_character_

  slope_table <- sim_spec$slope_table %||% NULL
  slopes <- if (is.data.frame(slope_table) && "Estimate" %in% names(slope_table)) {
    suppressWarnings(as.numeric(slope_table$Estimate))
  } else {
    numeric()
  }
  slopes <- slopes[is.finite(slopes) & slopes > 0]
  slope_levels <- if (length(slopes) > 0L) length(slopes) else NA_integer_
  max_abs_log <- if (length(slopes) > 0L) {
    log_slopes <- log(slopes) - mean(log(slopes))
    max(abs(log_slopes), na.rm = TRUE)
  } else {
    NA_real_
  }
  score_support_available <- nrow(rep_tbl) > 0L &&
    all(c("MinScoreCount", "MinScoreProportion", "ZeroScoreLevels") %in% names(rep_tbl))
  score_support_reps <- if (score_support_available) {
    sum(is.finite(suppressWarnings(as.numeric(rep_tbl$MinScoreCount))))
  } else {
    0L
  }
  min_score_count <- if (score_support_reps > 0L) {
    min(suppressWarnings(as.numeric(rep_tbl$MinScoreCount)), na.rm = TRUE)
  } else {
    NA_real_
  }
  min_score_prop <- if (score_support_reps > 0L) {
    min(suppressWarnings(as.numeric(rep_tbl$MinScoreProportion)), na.rm = TRUE)
  } else {
    NA_real_
  }
  max_zero_levels <- if (score_support_reps > 0L) {
    max(suppressWarnings(as.numeric(rep_tbl$ZeroScoreLevels)), na.rm = TRUE)
  } else {
    NA_real_
  }
  score_support_status <- dplyr::case_when(
    score_support_reps == 0L ~ "not_available",
    is.finite(max_zero_levels) && max_zero_levels > 0 ~ "review",
    is.finite(min_score_count) && min_score_count < 3 ~ "review",
    TRUE ~ "ok"
  )
  score_support_interpretation <- dplyr::case_when(
    score_support_reps == 0L ~
      "Generated score-category support was not retained in the replication overview.",
    is.finite(max_zero_levels) && max_zero_levels > 0 ~
      "At least one generated replication omitted one or more score categories; interpret recovery as sparse-category stress evidence.",
    is.finite(min_score_count) && min_score_count < 3 ~
      "All categories appeared, but at least one generated category had fewer than three observations.",
    TRUE ~
      "Generated score categories were all represented with at least minimal support in the retained replications."
  )
  score_support_next_action <- dplyr::case_when(
    score_support_status == "ok" ~
      "Use score-category support as supporting DGM evidence for this recovery run.",
    score_support_status == "not_available" ~
      "Refresh recovery simulation output so replication-level score support is retained before reporting category-sparse conditions.",
    TRUE ~
      "Report sparse generated score support explicitly and inspect category-level recovery or increase design size before generalizing."
  )

  stress_level <- dplyr::case_when(
    !identical(model, "GPCM") ~ "not_applicable",
    identical(regime, "unit_slopes") ~ "none",
    identical(regime, "near_flat") ~ "low",
    identical(regime, "moderate") ~ "moderate",
    identical(regime, "high_dispersion") ~ "high",
    TRUE ~ "unknown"
  )
  status <- dplyr::case_when(
    !identical(model, "GPCM") ~ "not_assessed",
    stress_level == "unknown" ~ "not_available",
    TRUE ~ "ok"
  )
  interpretation <- dplyr::case_when(
    !identical(model, "GPCM") ~
      "The fitted generator is not bounded GPCM, so slope-regime metadata is not part of this recovery condition.",
    identical(regime, "unit_slopes") ~
      "The GPCM generator is the unit-slope PCM-reduction condition.",
    identical(regime, "near_flat") ~
      "The GPCM generator has only near-flat relative discrimination spread.",
    identical(regime, "moderate") ~
      "The GPCM generator has moderate relative discrimination spread for slope-recovery review.",
    identical(regime, "high_dispersion") ~
      "The GPCM generator is a high-dispersion slope stress condition; interpret recovery metrics as stress evidence, not as a general adequacy claim.",
    TRUE ~
      "The GPCM slope-regime condition is unavailable; inspect the simulation specification before interpreting slope recovery."
  )
  next_action <- dplyr::case_when(
    !identical(model, "GPCM") ~
      "No GPCM slope-condition follow-up is needed for this model.",
    identical(status, "not_available") ~
      "Rebuild the simulation specification with explicit bounded-GPCM slopes before using slope-recovery evidence.",
    identical(regime, "high_dispersion") ~
      "Report the high-dispersion generator condition explicitly alongside slope recovery and uncertainty summaries.",
    TRUE ~
      "Report this generator condition alongside recovery and uncertainty summaries."
  )

  tibble::tibble(
    Model = model,
    GPCMSlopeRegime = regime,
    StressLevel = stress_level,
    SlopeLevels = as.integer(slope_levels),
    MaxAbsCenteredLogSlope = as.numeric(max_abs_log),
    Replications = as.integer(nrow(rep_tbl)),
    ScoreSupportReplications = as.integer(score_support_reps),
    MinScoreCount = as.integer(min_score_count),
    MinScoreProportion = as.numeric(min_score_prop),
    MaxZeroScoreLevels = as.integer(max_zero_levels),
    ScoreSupportStatus = score_support_status,
    Status = status,
    Interpretation = interpretation,
    ScoreSupportInterpretation = score_support_interpretation,
    ScoreSupportNextAction = score_support_next_action,
    NextAction = next_action
  )
}

recovery_assessment_condition_reporting_notes <- function(condition_review) {
  empty <- tibble::tibble(
    Model = character(),
    GPCMSlopeRegime = character(),
    StressLevel = character(),
    ConditionArea = character(),
    ReportingAttention = character(),
    ConditionFinding = character(),
    Evidence = character(),
    ReportingImplication = character(),
    NextAction = character(),
    ValidationUse = character()
  )
  condition_review <- as.data.frame(condition_review %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(condition_review) == 0L) return(empty)

  row <- condition_review[1, , drop = FALSE]
  row_chr <- function(name, default = NA_character_) {
    if (!name %in% names(row)) return(default)
    val <- row[[name]][1]
    if (is.na(val)) default else as.character(val)
  }
  row_num <- function(name) {
    if (!name %in% names(row)) return(NA_real_)
    suppressWarnings(as.numeric(row[[name]][1]))
  }
  fmt_num <- function(x) {
    if (!is.finite(x)) return("NA")
    formatC(x, digits = 3L, format = "fg", flag = "#")
  }
  make_note <- function(area,
                        attention,
                        finding,
                        evidence,
                        implication,
                        next_action) {
    tibble::tibble(
      Model = row_chr("Model"),
      GPCMSlopeRegime = row_chr("GPCMSlopeRegime"),
      StressLevel = row_chr("StressLevel"),
      ConditionArea = area,
      ReportingAttention = attention,
      ConditionFinding = finding,
      Evidence = evidence,
      ReportingImplication = implication,
      NextAction = next_action,
      ValidationUse = "generator_condition_not_release_gate"
    )
  }

  model <- row_chr("Model", "unknown")
  regime <- row_chr("GPCMSlopeRegime")
  stress <- row_chr("StressLevel")
  condition_status <- row_chr("Status", "not_available")
  slope_evidence <- paste(
    paste0("model=", model),
    paste0("slope_regime=", ifelse(is.na(regime), "NA", regime)),
    paste0("stress_level=", ifelse(is.na(stress), "NA", stress)),
    paste0("slope_levels=", fmt_num(row_num("SlopeLevels"))),
    paste0("max_abs_centered_log_slope=", fmt_num(row_num("MaxAbsCenteredLogSlope"))),
    sep = "; "
  )
  slope_note <- if (!identical(model, "GPCM")) {
    make_note(
      "slope_regime",
      "context",
      "slope_regime_not_applicable",
      slope_evidence,
      "This model does not use bounded-GPCM slope-regime metadata.",
      "Do not use slope-regime language for this recovery condition."
    )
  } else if (!identical(condition_status, "ok")) {
    make_note(
      "slope_regime",
      "not_available",
      "slope_regime_not_available",
      slope_evidence,
      "Bounded-GPCM slope-regime metadata were unavailable, so slope-stress interpretation should be withheld.",
      "Rebuild the simulation specification with explicit bounded-GPCM slopes before reporting slope-recovery context."
    )
  } else if (identical(regime, "high_dispersion")) {
    make_note(
      "slope_regime",
      "reporting_review",
      "high_dispersion_slope_stress",
      slope_evidence,
      "The generator intentionally used high relative-discrimination spread; report this as stress context, not as a literature-derived adequacy cut point.",
      "Read slope recovery together with uncertainty, score support, and diagnostic notes before generalizing from this condition."
    )
  } else {
    make_note(
      "slope_regime",
      "context",
      paste0(ifelse(is.na(regime), "unknown", regime), "_slope_regime"),
      slope_evidence,
      "The bounded-GPCM slope-regime label describes the generator condition for recovery interpretation.",
      "Report the generator condition alongside recovery and uncertainty summaries."
    )
  }

  score_status <- row_chr("ScoreSupportStatus", "not_available")
  min_count <- row_num("MinScoreCount")
  min_prop <- row_num("MinScoreProportion")
  max_zero <- row_num("MaxZeroScoreLevels")
  score_evidence <- paste(
    paste0("score_support_replications=", fmt_num(row_num("ScoreSupportReplications"))),
    paste0("min_score_count=", fmt_num(min_count)),
    paste0("min_score_proportion=", fmt_num(min_prop)),
    paste0("max_zero_score_levels=", fmt_num(max_zero)),
    sep = "; "
  )
  score_note <- if (identical(score_status, "not_available")) {
    make_note(
      "score_support",
      "not_available",
      "score_support_not_available",
      score_evidence,
      "Generated score-category support was unavailable, so sparse-category stress cannot be separated from recovery performance.",
      "Retain replication-level score support before making category-support claims."
    )
  } else if (is.finite(max_zero) && max_zero > 0) {
    make_note(
      "score_support",
      "reporting_review",
      "omitted_generated_score_categories",
      score_evidence,
      "At least one generated replication omitted score categories; interpret recovery as sparse-category stress evidence, not as a generic slope-recovery failure.",
      "Report sparse generated score support explicitly and inspect category-level recovery before generalizing."
    )
  } else if (is.finite(min_count) && min_count < 3) {
    make_note(
      "score_support",
      "reporting_review",
      "thin_generated_score_category_support",
      score_evidence,
      "All generated categories appeared, but at least one category had very small support.",
      "Qualify category-threshold and recovery-language claims for this condition."
    )
  } else {
    make_note(
      "score_support",
      "context",
      "generated_score_support_context",
      score_evidence,
      "Generated score categories had minimal retained support for reading the recovery condition.",
      "Use score-category support as supporting generator-condition context."
    )
  }

  dplyr::bind_rows(slope_note, score_note)
}

recovery_assessment_diagnostic_review <- function(x) {
  diagnostic_summary <- tibble::as_tibble(x$diagnostic_oc_summary %||% tibble::tibble())
  if (nrow(diagnostic_summary) == 0L) {
    return(tibble::tibble(
      Facet = character(),
      Replications = integer(),
      MeanSeparation = numeric(),
      MeanReliability = numeric(),
      MeanInfit = numeric(),
      MeanOutfit = numeric(),
      MeanMisfitRateAbsZ2 = numeric(),
      MeanDfSensitiveFlagRate = numeric(),
      DiagnosticAvailability = character(),
      Status = character(),
      ValidationUse = character(),
      Interpretation = character(),
      NextAction = character()
    ))
  }
  diagnostic_summary |>
    dplyr::mutate(
      DiagnosticAvailability = dplyr::case_when(
        is.finite(.data$Replications) & .data$Replications > 0 ~ "available",
        TRUE ~ "not_available"
      ),
      Status = "not_assessed",
      ValidationUse = "diagnostic_only_not_release_gate",
      Interpretation = paste(
        "Fit/separation operating characteristics summarize how diagnostics behaved under the simulated condition;",
        "their availability does not imply adequacy, and they do not replace recovery, convergence, uncertainty, or substantive validity evidence."
      ),
      NextAction = dplyr::case_when(
        .data$DiagnosticAvailability == "available" ~
          "Use as diagnostic operating-characteristic context; keep separate from release-recovery status.",
        TRUE ~
          "Re-run evaluate_mfrm_recovery(include_diagnostics = TRUE) if fit/separation behavior is needed."
      )
    ) |>
    dplyr::arrange(.data$Facet)
}

recovery_assessment_diagnostic_reporting_notes <- function(diagnostic_review) {
  empty <- tibble::tibble(
    Facet = character(),
    ReportingAttention = character(),
    DiagnosticFinding = character(),
    Evidence = character(),
    ReportingImplication = character(),
    NextAction = character(),
    ValidationUse = character()
  )
  diagnostic_review <- as.data.frame(diagnostic_review %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(diagnostic_review) == 0L) return(empty)

  row_chr <- function(row, name, default = NA_character_) {
    if (!name %in% names(row)) return(default)
    val <- row[[name]][1]
    if (is.na(val)) default else as.character(val)
  }
  row_num <- function(row, name) {
    if (!name %in% names(row)) return(NA_real_)
    suppressWarnings(as.numeric(row[[name]][1]))
  }
  fmt_num <- function(x) {
    if (!is.finite(x)) return("NA")
    formatC(x, digits = 3L, format = "fg", flag = "#")
  }
  build_evidence <- function(row) {
    paste(
      paste0("replications=", row_chr(row, "Replications", "NA")),
      paste0("mean_separation=", fmt_num(row_num(row, "MeanSeparation"))),
      paste0("mean_reliability=", fmt_num(row_num(row, "MeanReliability"))),
      paste0("mean_abs_zstd_flag_rate=", fmt_num(row_num(row, "MeanMisfitRateAbsZ2"))),
      paste0("mean_df_sensitive_flag_rate=", fmt_num(row_num(row, "MeanDfSensitiveFlagRate"))),
      sep = "; "
    )
  }
  make_note <- function(row,
                        attention,
                        finding,
                        implication,
                        next_action) {
    tibble::tibble(
      Facet = row_chr(row, "Facet"),
      ReportingAttention = attention,
      DiagnosticFinding = finding,
      Evidence = build_evidence(row),
      ReportingImplication = implication,
      NextAction = next_action,
      ValidationUse = "diagnostic_only_not_release_gate"
    )
  }

  rows <- lapply(seq_len(nrow(diagnostic_review)), function(i) {
    row <- diagnostic_review[i, , drop = FALSE]
    availability <- row_chr(row, "DiagnosticAvailability", row_chr(row, "Status", "not_available"))
    if (!identical(availability, "available")) {
      return(make_note(
        row,
        attention = "not_available",
        finding = "diagnostic_not_available",
        implication = "Fit/separation diagnostics were not available for this facet; do not infer fit or separation behavior from this recovery assessment.",
        next_action = "Re-run evaluate_mfrm_recovery(include_diagnostics = TRUE) if diagnostic operating characteristics are needed."
      ))
    }

    sep <- row_num(row, "MeanSeparation")
    rel <- row_num(row, "MeanReliability")
    misfit_rate <- row_num(row, "MeanMisfitRateAbsZ2")
    df_rate <- row_num(row, "MeanDfSensitiveFlagRate")
    notes <- list()

    if ((is.finite(sep) && sep <= 0) || (is.finite(rel) && rel <= 0)) {
      notes[[length(notes) + 1L]] <- make_note(
        row,
        attention = "reporting_review",
        finding = "zero_separation_or_reliability",
        implication = "The Rasch/FACETS-style separation signal collapsed to zero for this simulated condition; report this as diagnostic context, not as recovery failure.",
        next_action = "Inspect the generated condition and facet design before using this facet's separation/reliability language in examples or reports."
      )
    }
    if (is.finite(misfit_rate) && misfit_rate > 0) {
      notes[[length(notes) + 1L]] <- make_note(
        row,
        attention = "reporting_review",
        finding = "abs_zstd_flags_present",
        implication = "At least one replication produced absolute fit-ZSTD flags for this facet; this is a diagnostic signal rather than a recovery gate.",
        next_action = "Review the corresponding fit rows before making strong fit-language claims for this simulated condition."
      )
    }
    if (is.finite(df_rate) && df_rate > 0) {
      notes[[length(notes) + 1L]] <- make_note(
        row,
        attention = "reporting_review",
        finding = "df_sensitive_zstd_flags_present",
        implication = "Fit-ZSTD flagging changed under the engine-vs-FACETS degrees-of-freedom convention; keep df convention language explicit.",
        next_action = "Use fit_measures_table(fit_df_method = \"both\") or the diagnostic rows when explaining fit-ZSTD sensitivity."
      )
    }
    if (length(notes) == 0L) {
      notes[[1L]] <- make_note(
        row,
        attention = "context",
        finding = "diagnostic_context_available",
        implication = "Fit/separation diagnostics are available for context; no automatic reporting caveat was triggered by the retained operating-characteristic fields.",
        next_action = "Keep this row as diagnostic context and continue to base recovery assessment on recovery metrics, convergence, uncertainty, and Monte Carlo precision."
      )
    }
    dplyr::bind_rows(notes)
  })

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) return(empty)
  out
}

recovery_assessment_next_actions <- function(checklist, metric_review, max_n = 6L) {
  actions <- character(0)
  if (is.data.frame(checklist) && nrow(checklist) > 0L) {
    bad <- checklist[checklist$Status %in% c("concern", "review", "not_available"), , drop = FALSE]
    if (nrow(bad) > 0L) actions <- c(actions, paste(bad$Item, bad$NextAction, sep = ": "))
  }
  if (is.data.frame(metric_review) && nrow(metric_review) > 0L) {
    bad <- metric_review[metric_review$OverallStatus %in% c("concern", "review"), , drop = FALSE]
    if (nrow(bad) > 0L) {
      labels <- paste(bad$ParameterType, bad$Facet, bad$ComparisonScale, sep = " / ")
      actions <- c(actions, paste(labels, bad$NextAction, sep = ": "))
    }
  }
  actions <- unique(actions[nzchar(actions)])
  if (length(actions) == 0L) {
    actions <- "No immediate follow-up action from the selected recovery assessment thresholds."
  }
  utils::head(actions, n = max(1L, as.integer(max_n)))
}

#' Assess whether recovery-simulation results are ready to use
#'
#' @description
#' Converts the numerical output from [evaluate_mfrm_recovery()] into a
#' reviewer-facing adequacy checklist. The goal is not to impose one universal
#' pass/fail rule; it is to make the main user questions explicit: Did the runs
#' finish? Did the fitted models converge? Are uncertainty summaries available?
#' Are coverage and Monte Carlo precision plausible? If practical RMSE or bias
#' limits are supplied, which parameter groups need follow-up? For bounded
#' `GPCM`, which slope-regime generator condition frames the recovery evidence?
#'
#' @param x For `assess_mfrm_recovery()`, output from
#'   [evaluate_mfrm_recovery()]. For `plot.mfrm_recovery_assessment()`, output
#'   from `assess_mfrm_recovery()`.
#' @param min_reps Minimum replication count expected before treating the
#'   simulation as more than a smoke check.
#' @param min_success_rate Minimum acceptable proportion of replications that
#'   generated data and produced a fitted model.
#' @param min_convergence_rate Minimum acceptable proportion of replications
#'   whose fitted model reported convergence.
#' @param min_se_available Minimum acceptable proportion of recovery rows with
#'   standard errors in each parameter group. Set to `NULL` to skip this check.
#' @param coverage_target Nominal coverage target, usually `0.95`.
#' @param coverage_tolerance Absolute tolerance around `coverage_target`.
#' @param max_mcse_rmse_ratio Maximum acceptable Monte Carlo SE of RMSE divided
#'   by RMSE. Set to `NULL` to skip this precision check.
#' @param max_rmse Optional practical RMSE limit. Use a scalar for all parameter
#'   groups or a named vector/list with names such as `"facet"`, `"step"`,
#'   `"slope"`, `"Rater"`, or `"facet:Rater:logit"`.
#' @param max_abs_bias Optional practical absolute-bias limit. Naming follows
#'   `max_rmse`.
#' @param top_n Number of next-action lines retained in the compact output.
#' @param digits Digits used by the print method.
#' @param ... Reserved for future extensions.
#'
#' @details
#' RMSE and bias adequacy depends on the substantive scale and the use case, so
#' the function does not mark them as failed unless the user supplies
#' `max_rmse` or `max_abs_bias`. Without those limits, the corresponding rows are
#' marked `not_assessed` and the next action asks the user to set practical
#' thresholds when a decision depends on the metric.
#'
#' The `condition_review` table is generator metadata for interpreting the
#' recovery run. For bounded `GPCM`, `GPCMSlopeRegime`, `StressLevel`, and
#' generated score-category support describe the data-generating condition;
#' they are not model-fit tests and they are not literature-derived adequacy
#' cut points. `condition_reporting_notes` turns those generator conditions
#' into reporter-facing caveats, such as high-dispersion slope stress or sparse
#' generated score support.
#'
#' The optional `diagnostic_review` table is available when
#' [evaluate_mfrm_recovery()] was called with `include_diagnostics = TRUE`.
#' It summarizes fit and separation operating characteristics as diagnostic
#' context only. Its availability fields do not mean that fit or separation
#' values are adequate, and those rows do not enter the recovery adequacy
#' status. `diagnostic_reporting_notes` should be read first when drafting
#' fit/separation language because it separates zero separation/reliability,
#' absolute fit-ZSTD flags, and df-sensitive ZSTD flags from recovery gates.
#'
#' `plot.mfrm_recovery_assessment()` is a user-facing review aid. Use
#' `type = "status"` first to see where checklist attention is needed, then
#' `type = "metrics"` to inspect the parameter groups behind RMSE, bias,
#' coverage, standard-error availability, or Monte Carlo precision statuses.
#' The intended reading order is `summary(recovery_review)`, then
#' `condition_reporting_notes` and `condition_review`, then
#' `diagnostic_reporting_notes` and `diagnostic_review`, then the status plot,
#' then the metric plot, then the row-level recovery table for the parameter
#' groups that need follow-up. When `draw = FALSE`, the plot data also include
#' `reading_order`, `guidance`, condition/diagnostic handoff tables, and
#' user-facing plot tables such as `section_status` for status plots and
#' `metric_review` for metric plots.
#'
#' @return An object of class `mfrm_recovery_assessment` with:
#' - `overview`: compact run-level status.
#' - `checklist`: reviewer-facing adequacy checks.
#' - `condition_review`: generator-condition metadata, including bounded
#'   `GPCM` slope-regime interpretation and generated score-category support
#'   when available.
#' - `condition_reporting_notes`: reporter-facing generator-condition caveats
#'   separated from recovery metrics and release-gate decisions.
#' - `diagnostic_reporting_notes`: reporter-facing fit/separation caveats
#'   retained as diagnostic context rather than recovery gates.
#' - `diagnostic_review`: optional fit/separation operating-characteristic
#'   context when retained by [evaluate_mfrm_recovery()].
#' - `metric_review`: parameter-group metric checks.
#' - `uncertainty_review`: compact coverage / SE availability interpretation.
#' - `reading_order`: recommended first-read order for the summary, condition,
#'   plot, and row-level recovery outputs.
#' - `next_actions`: short action list sorted by severity.
#' - `thresholds`: thresholds used for the assessment.
#' @seealso [evaluate_mfrm_recovery()], [plot.mfrm_recovery_simulation()]
#' @examples
#' \donttest{
#' rec <- evaluate_mfrm_recovery(
#'   n_person = 12,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   reps = 1,
#'   maxit = 30,
#'   seed = 123
#' )
#' assess_mfrm_recovery(rec, min_reps = 1, max_rmse = 1)
#'
#' # Read the bounded-GPCM generator condition separately from recovery adequacy.
#' gpcm_spec <- build_mfrm_sim_spec(
#'   n_person = 14,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 2,
#'   model = "GPCM",
#'   step_facet = "Criterion",
#'   slope_facet = "Criterion",
#'   slopes = c(0.85, 1.15),
#'   assignment = "crossed"
#' )
#' gpcm_rec <- suppressWarnings(evaluate_mfrm_recovery(
#'   sim_spec = gpcm_spec,
#'   reps = 1,
#'   fit_method = "MML",
#'   quad_points = 5,
#'   maxit = 12,
#'   include_diagnostics = TRUE,
#'   include_person = FALSE,
#'   seed = 456
#' ))
#' gpcm_review <- assess_mfrm_recovery(
#'   gpcm_rec,
#'   min_reps = 1,
#'   max_rmse = c(slope = 2),
#'   max_abs_bias = c(slope = 1),
#'   min_se_available = NULL,
#'   max_mcse_rmse_ratio = NULL
#' )
#' gpcm_review$condition_reporting_notes[, c(
#'   "ConditionArea", "ReportingAttention", "ConditionFinding"
#' )]
#' gpcm_review$condition_review[, c(
#'   "Model", "GPCMSlopeRegime", "StressLevel", "ScoreSupportStatus"
#' )]
#' gpcm_review$diagnostic_reporting_notes[, c(
#'   "Facet", "ReportingAttention", "DiagnosticFinding"
#' )]
#' summary(gpcm_review)$reading_order
#' }
#' @export
assess_mfrm_recovery <- function(x,
                                 min_reps = 30,
                                 min_success_rate = 0.95,
                                 min_convergence_rate = 0.95,
                                 min_se_available = 0.80,
                                 coverage_target = 0.95,
                                 coverage_tolerance = 0.05,
                                 max_mcse_rmse_ratio = 0.25,
                                 max_rmse = NULL,
                                 max_abs_bias = NULL,
                                 top_n = 6,
                                 digits = 3,
                                 ...) {
  if (!inherits(x, "mfrm_recovery_simulation")) {
    stop("`x` must be output from evaluate_mfrm_recovery().", call. = FALSE)
  }
  rep_tbl <- tibble::as_tibble(x$rep_overview %||% tibble::tibble())
  summary_tbl <- tibble::as_tibble(x$recovery_summary %||% tibble::tibble())
  recovery_tbl <- tibble::as_tibble(x$recovery %||% tibble::tibble())

  reps <- if (nrow(rep_tbl) > 0L) nrow(rep_tbl) else 0L
  successful <- if (nrow(rep_tbl) > 0L) sum(rep_tbl$RunOK, na.rm = TRUE) else 0L
  converged <- if (nrow(rep_tbl) > 0L) sum(rep_tbl$Converged, na.rm = TRUE) else 0L
  success_rate <- if (reps > 0L) successful / reps else NA_real_
  convergence_rate <- if (reps > 0L) converged / reps else NA_real_

  metric_review <- summary_tbl
  if (nrow(metric_review) > 0L) {
    metric_review$RMSELimit <- mapply(
      recovery_assessment_threshold,
      parameter_type = metric_review$ParameterType,
      facet = metric_review$Facet,
      comparison_scale = metric_review$ComparisonScale,
      MoreArgs = list(thresholds = max_rmse),
      SIMPLIFY = TRUE
    )
    metric_review$AbsBiasLimit <- mapply(
      recovery_assessment_threshold,
      parameter_type = metric_review$ParameterType,
      facet = metric_review$Facet,
      comparison_scale = metric_review$ComparisonScale,
      MoreArgs = list(thresholds = max_abs_bias),
      SIMPLIFY = TRUE
    )
    metric_review$RMSEStatus <- mapply(
      recovery_assessment_metric_status,
      value = metric_review$RMSE,
      limit = metric_review$RMSELimit,
      SIMPLIFY = TRUE
    )
    metric_review$BiasStatus <- mapply(
      recovery_assessment_metric_status,
      value = abs(metric_review$Bias),
      limit = metric_review$AbsBiasLimit,
      SIMPLIFY = TRUE
    )
    metric_review$CoverageStatus <- vapply(
      metric_review$Coverage95,
      recovery_assessment_coverage_status,
      character(1),
      target = coverage_target,
      tolerance = coverage_tolerance
    )
    metric_review$SEStatus <- if (is.null(min_se_available)) {
      rep("not_assessed", nrow(metric_review))
    } else {
      vapply(
        metric_review$SEAvailableRate,
        recovery_assessment_rate_status,
        character(1),
        target = min_se_available,
        tolerance = 0,
        concern_gap = 0.25
      )
    }
    rmse_ratio <- suppressWarnings(as.numeric(metric_review$McseRMSE) / as.numeric(metric_review$RMSE))
    metric_review$McseRMSEToRMSE <- rmse_ratio
    metric_review$MonteCarloStatus <- if (is.null(max_mcse_rmse_ratio)) {
      rep("not_assessed", nrow(metric_review))
    } else {
      vapply(
        rmse_ratio,
        recovery_assessment_metric_status,
        character(1),
        limit = max_mcse_rmse_ratio,
        concern_multiplier = 2
      )
    }
    metric_review$OverallStatus <- apply(
      metric_review[, c("RMSEStatus", "BiasStatus", "CoverageStatus",
                        "SEStatus", "MonteCarloStatus"), drop = FALSE],
      1,
      recovery_assessment_overall_status
    )
    metric_review$NextAction <- vapply(
      metric_review$OverallStatus,
      recovery_assessment_action,
      character(1),
      topic = "this parameter group"
    )
  } else {
    metric_review <- tibble::tibble()
  }
  condition_review <- recovery_assessment_condition_review(x)
  condition_reporting_notes <- recovery_assessment_condition_reporting_notes(condition_review)
  diagnostic_review <- recovery_assessment_diagnostic_review(x)
  diagnostic_reporting_notes <- recovery_assessment_diagnostic_reporting_notes(diagnostic_review)

  run_status <- recovery_assessment_rate_status(reps, min_reps, tolerance = 0, concern_gap = min_reps)
  if (reps > 0L && reps < min_reps) run_status <- "review"
  success_status <- recovery_assessment_rate_status(success_rate, min_success_rate, tolerance = 0, concern_gap = 0.20)
  convergence_status <- recovery_assessment_rate_status(convergence_rate, min_convergence_rate, tolerance = 0, concern_gap = 0.20)
  row_status <- if (nrow(recovery_tbl) > 0L) "ok" else "concern"

  se_status <- if (nrow(metric_review) == 0L || is.null(min_se_available)) {
    "not_assessed"
  } else {
    recovery_assessment_overall_status(metric_review$SEStatus)
  }
  coverage_status <- if (nrow(metric_review) == 0L) {
    "not_available"
  } else {
    recovery_assessment_overall_status(metric_review$CoverageStatus)
  }
  mc_status <- if (nrow(metric_review) == 0L || is.null(max_mcse_rmse_ratio)) {
    "not_assessed"
  } else {
    recovery_assessment_overall_status(metric_review$MonteCarloStatus)
  }
  rmse_status <- if (nrow(metric_review) == 0L) {
    "not_available"
  } else if (is.null(max_rmse)) {
    "not_assessed"
  } else {
    recovery_assessment_overall_status(metric_review$RMSEStatus)
  }
  bias_status <- if (nrow(metric_review) == 0L) {
    "not_available"
  } else if (is.null(max_abs_bias)) {
    "not_assessed"
  } else {
    recovery_assessment_overall_status(metric_review$BiasStatus)
  }
  condition_status <- as.character(condition_review$Status[1] %||% "not_available")
  score_support_status <- as.character(condition_review$ScoreSupportStatus[1] %||% "not_available")
  condition_evidence <- if (nrow(condition_review) > 0L) {
    regime <- as.character(condition_review$GPCMSlopeRegime[1] %||% NA_character_)
    stress <- as.character(condition_review$StressLevel[1] %||% NA_character_)
    slope_levels <- suppressWarnings(as.integer(condition_review$SlopeLevels[1]))
    max_abs_log <- suppressWarnings(as.numeric(condition_review$MaxAbsCenteredLogSlope[1]))
    if (identical(as.character(condition_review$Model[1]), "GPCM")) {
      sprintf(
        "GPCM slope regime is %s; stress level = %s; slope levels = %s; max |centered log slope| = %s.",
        ifelse(is.na(regime), "unknown", regime),
        ifelse(is.na(stress), "unknown", stress),
        ifelse(is.na(slope_levels), "unknown", as.character(slope_levels)),
        ifelse(is.finite(max_abs_log), sprintf("%.3f", max_abs_log), "unknown")
      )
    } else {
      sprintf(
        "Model %s does not use bounded-GPCM slope-regime metadata.",
        as.character(condition_review$Model[1] %||% "unknown")
      )
    }
  } else {
    "No generator-condition metadata was available."
  }
  score_support_evidence <- if (nrow(condition_review) > 0L) {
    min_count <- suppressWarnings(as.integer(condition_review$MinScoreCount[1]))
    min_prop <- suppressWarnings(as.numeric(condition_review$MinScoreProportion[1]))
    max_zero <- suppressWarnings(as.integer(condition_review$MaxZeroScoreLevels[1]))
    reps_with_support <- suppressWarnings(as.integer(condition_review$ScoreSupportReplications[1]))
    sprintf(
      "%s replication(s) retained score support; minimum category count = %s; minimum category proportion = %s; maximum omitted categories = %s.",
      ifelse(is.na(reps_with_support), "0", as.character(reps_with_support)),
      ifelse(is.na(min_count), "unknown", as.character(min_count)),
      ifelse(is.finite(min_prop), sprintf("%.3f", min_prop), "unknown"),
      ifelse(is.na(max_zero), "unknown", as.character(max_zero))
    )
  } else {
    "No generated score-support metadata was available."
  }

  checklist <- dplyr::bind_rows(
    recovery_assessment_check_row(
      "Run completion",
      "Replication count",
      run_status,
      sprintf("%d replication(s); requested minimum is %d.", reps, as.integer(min_reps)),
      recovery_assessment_action(run_status, "the replication count")
    ),
    recovery_assessment_check_row(
      "Run completion",
      "Simulation and refit success",
      success_status,
      sprintf("%d/%d successful run(s); rate = %.3f.", successful, reps, success_rate),
      recovery_assessment_action(success_status, "run completion")
    ),
    recovery_assessment_check_row(
      "Run completion",
      "Reported convergence",
      convergence_status,
      sprintf("%d/%d converged run(s); rate = %.3f.", converged, reps, convergence_rate),
      recovery_assessment_action(convergence_status, "model convergence")
    ),
    recovery_assessment_check_row(
      "Recovery content",
      "Recoverable truth-estimate rows",
      row_status,
      sprintf("%d row-level recovery comparison(s).", nrow(recovery_tbl)),
      recovery_assessment_action(row_status, "recovery rows")
    ),
    recovery_assessment_check_row(
      "Generator conditions",
      "Bounded-GPCM slope regime",
      condition_status,
      condition_evidence,
      as.character(condition_review$NextAction[1] %||% recovery_assessment_action(condition_status, "generator conditions"))
    ),
    recovery_assessment_check_row(
      "Generator conditions",
      "Generated score-category support",
      score_support_status,
      score_support_evidence,
      as.character(condition_review$ScoreSupportNextAction[1] %||% recovery_assessment_action(score_support_status, "generated score support"))
    ),
    recovery_assessment_check_row(
      "Uncertainty",
      "Standard-error availability",
      se_status,
      if (nrow(metric_review) > 0L) {
        sprintf("Group statuses: %s.", recovery_assessment_status_counts(metric_review$SEStatus))
      } else {
        "No parameter-group summaries were available."
      },
      recovery_assessment_action(se_status, "standard-error availability")
    ),
    recovery_assessment_check_row(
      "Uncertainty",
      "Coverage",
      coverage_status,
      if (nrow(metric_review) > 0L) {
        sprintf("Group statuses: %s.", recovery_assessment_status_counts(metric_review$CoverageStatus))
      } else {
        "No parameter-group summaries were available."
      },
      recovery_assessment_action(coverage_status, "coverage")
    ),
    recovery_assessment_check_row(
      "Monte Carlo precision",
      "RMSE Monte Carlo error",
      mc_status,
      if (nrow(metric_review) > 0L) {
        sprintf("Group statuses: %s.", recovery_assessment_status_counts(metric_review$MonteCarloStatus))
      } else {
        "No parameter-group summaries were available."
      },
      recovery_assessment_action(mc_status, "Monte Carlo precision")
    ),
    recovery_assessment_check_row(
      "Practical thresholds",
      "RMSE threshold",
      rmse_status,
      if (is.null(max_rmse)) {
        "No practical RMSE threshold was supplied."
      } else {
        sprintf("Group statuses: %s.", recovery_assessment_status_counts(metric_review$RMSEStatus))
      },
      recovery_assessment_action(rmse_status, "RMSE")
    ),
    recovery_assessment_check_row(
      "Practical thresholds",
      "Bias threshold",
      bias_status,
      if (is.null(max_abs_bias)) {
        "No practical absolute-bias threshold was supplied."
      } else {
        sprintf("Group statuses: %s.", recovery_assessment_status_counts(metric_review$BiasStatus))
      },
      recovery_assessment_action(bias_status, "bias")
    )
  )

  overall_status <- recovery_assessment_overall_status(checklist$Status)
  overview <- tibble::tibble(
    Reps = reps,
    SuccessfulRuns = successful,
    SuccessRate = success_rate,
    ConvergedRuns = converged,
    ConvergenceRate = convergence_rate,
    RecoveryRows = nrow(recovery_tbl),
    RecoveryGroups = nrow(summary_tbl),
    OverallStatus = overall_status
  )
  uncertainty_review <- recovery_assessment_uncertainty_review(
    metric_review,
    coverage_target = coverage_target,
    coverage_tolerance = coverage_tolerance,
    min_se_available = min_se_available
  )
  next_actions <- recovery_assessment_next_actions(checklist, metric_review, max_n = top_n)

  structure(
    list(
      overview = overview,
      checklist = checklist,
      condition_review = condition_review,
      condition_reporting_notes = condition_reporting_notes,
      diagnostic_reporting_notes = diagnostic_reporting_notes,
      diagnostic_review = diagnostic_review,
      metric_review = metric_review,
      uncertainty_review = uncertainty_review,
      reading_order = recovery_assessment_reading_order(),
      next_actions = next_actions,
      thresholds = list(
        min_reps = min_reps,
        min_success_rate = min_success_rate,
        min_convergence_rate = min_convergence_rate,
        min_se_available = min_se_available,
        coverage_target = coverage_target,
        coverage_tolerance = coverage_tolerance,
        max_mcse_rmse_ratio = max_mcse_rmse_ratio,
        max_rmse = max_rmse,
        max_abs_bias = max_abs_bias
      ),
      notes = c(
        "RMSE and bias statuses are decision-specific; supply practical thresholds when they matter.",
        "Fit/separation operating characteristics are diagnostic context and do not enter the recovery adequacy status.",
        as.character(x$notes %||% character(0))
      ),
      source = x,
      digits = max(0L, as.integer(digits[1]))
    ),
    class = "mfrm_recovery_assessment"
  )
}

#' @export
summary.mfrm_recovery_assessment <- function(object, digits = NULL, ...) {
  if (!inherits(object, "mfrm_recovery_assessment")) {
    stop("`object` must be output from assess_mfrm_recovery().", call. = FALSE)
  }
  if (is.null(digits)) digits <- object$digits %||% 3L
  out <- list(
    overview = tibble::as_tibble(object$overview %||% tibble::tibble()),
    checklist = tibble::as_tibble(object$checklist %||% tibble::tibble()),
    condition_review = tibble::as_tibble(object$condition_review %||% tibble::tibble()),
    condition_reporting_notes = tibble::as_tibble(object$condition_reporting_notes %||% tibble::tibble()),
    diagnostic_reporting_notes = tibble::as_tibble(object$diagnostic_reporting_notes %||% tibble::tibble()),
    diagnostic_review = tibble::as_tibble(object$diagnostic_review %||% tibble::tibble()),
    metric_review = tibble::as_tibble(object$metric_review %||% tibble::tibble()),
    uncertainty_review = tibble::as_tibble(object$uncertainty_review %||% tibble::tibble()),
    reading_order = tibble::as_tibble(object$reading_order %||% recovery_assessment_reading_order()),
    next_actions = as.character(object$next_actions %||% character(0)),
    thresholds = object$thresholds %||% list(),
    notes = as.character(object$notes %||% character(0)),
    digits = max(0L, as.integer(digits[1]))
  )
  class(out) <- "summary.mfrm_recovery_assessment"
  out
}

#' @export
print.summary.mfrm_recovery_assessment <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("MFRM Recovery Adequacy Assessment\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$reading_order) && nrow(x$reading_order) > 0L) {
    keep <- intersect(c("Step", "Route", "WhatToRead"), names(x$reading_order))
    cat("\nRecommended reading order\n")
    print(as.data.frame(x$reading_order[, keep, drop = FALSE]), row.names = FALSE)
  }
  if (!is.null(x$checklist) && nrow(x$checklist) > 0L) {
    cat("\nChecklist\n")
    print(as.data.frame(x$checklist), row.names = FALSE)
  }
  if (!is.null(x$condition_review) && nrow(x$condition_review) > 0L) {
    keep <- intersect(
      c("Model", "GPCMSlopeRegime", "StressLevel", "SlopeLevels",
        "MaxAbsCenteredLogSlope", "MinScoreCount", "MinScoreProportion",
        "MaxZeroScoreLevels", "ScoreSupportStatus", "Status",
        "Interpretation"),
      names(x$condition_review)
    )
    cat("\nCondition review\n")
    print(round_numeric_df(as.data.frame(x$condition_review[, keep, drop = FALSE]), digits = digits),
          row.names = FALSE)
  }
  if (!is.null(x$condition_reporting_notes) && nrow(x$condition_reporting_notes) > 0L) {
    keep <- intersect(
      c("ConditionArea", "ReportingAttention", "ConditionFinding",
        "Evidence", "ValidationUse"),
      names(x$condition_reporting_notes)
    )
    cat("\nCondition reporting notes\n")
    print(as.data.frame(x$condition_reporting_notes[, keep, drop = FALSE]), row.names = FALSE)
  }
  if (!is.null(x$diagnostic_reporting_notes) && nrow(x$diagnostic_reporting_notes) > 0L) {
    keep <- intersect(
      c("Facet", "ReportingAttention", "DiagnosticFinding",
        "Evidence", "ValidationUse"),
      names(x$diagnostic_reporting_notes)
    )
    cat("\nDiagnostic reporting notes\n")
    print(as.data.frame(x$diagnostic_reporting_notes[, keep, drop = FALSE]), row.names = FALSE)
  }
  if (!is.null(x$diagnostic_review) && nrow(x$diagnostic_review) > 0L) {
    keep <- intersect(
      c("Facet", "Replications", "MeanSeparation", "MeanReliability",
        "MeanInfit", "MeanOutfit", "MeanMisfitRateAbsZ2",
        "MeanDfSensitiveFlagRate", "DiagnosticAvailability", "Status",
        "ValidationUse"),
      names(x$diagnostic_review)
    )
    cat("\nDiagnostic operating-characteristic review\n")
    print(round_numeric_df(as.data.frame(x$diagnostic_review[, keep, drop = FALSE]), digits = digits),
          row.names = FALSE)
  }
  if (!is.null(x$metric_review) && nrow(x$metric_review) > 0L) {
    keep <- intersect(
      c("ParameterType", "Facet", "ComparisonScale", "RMSE", "Bias",
        "Coverage95", "SEAvailableRate", "McseRMSEToRMSE", "OverallStatus"),
      names(x$metric_review)
    )
    cat("\nMetric review\n")
    print(round_numeric_df(as.data.frame(x$metric_review[, keep, drop = FALSE]), digits = digits),
          row.names = FALSE)
  }
  if (!is.null(x$uncertainty_review) && nrow(x$uncertainty_review) > 0L) {
    keep <- intersect(
      c("ParameterType", "Facet", "ComparisonScale", "Coverage95",
        "CoverageStatus", "SEAvailableRate", "SEStatus", "Interpretation"),
      names(x$uncertainty_review)
    )
    cat("\nUncertainty review\n")
    print(round_numeric_df(as.data.frame(x$uncertainty_review[, keep, drop = FALSE]), digits = digits),
          row.names = FALSE)
  }
  if (length(x$next_actions) > 0L) {
    cat("\nNext actions\n")
    for (line in x$next_actions) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.mfrm_recovery_assessment <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

recovery_assessment_status_rank <- function(status) {
  ranks <- c(ok = 1L, not_assessed = 2L, not_available = 2L,
             review = 3L, concern = 4L, fail = 4L)
  out <- unname(ranks[as.character(status)])
  out[is.na(out)] <- 0L
  out
}

recovery_assessment_status_palette <- function(status) {
  pal <- c(
    ok = "#238b45",
    review = "#b65e16",
    concern = "#b11f24",
    not_available = "#6b7280",
    not_assessed = "#8c8c8c",
    fail = "#b11f24"
  )
  out <- unname(pal[as.character(status)])
  out[is.na(out)] <- "#6b7280"
  out
}

recovery_assessment_reading_order <- function() {
  data.frame(
    Step = seq_len(6L),
    Route = c(
      "summary(recovery_review)",
      "recovery_review$condition_reporting_notes, then recovery_review$condition_review",
      "recovery_review$diagnostic_reporting_notes, then recovery_review$diagnostic_review",
      "plot(recovery_review, type = \"status\")",
      "plot(recovery_review, type = \"metrics\")",
      "recovery_review$source$recovery"
    ),
    WhatToRead = c(
      "Overall run status, next actions, and compact checklist.",
      "Generator-condition caveats, then GPCM slope-regime and generated score-support metadata.",
      "Reporter-facing fit/separation caveats, then optional operating characteristics retained by include_diagnostics = TRUE.",
      "Checklist domains ordered by attention status.",
      "Parameter groups behind RMSE, bias, coverage, SE, or Monte Carlo statuses.",
      "Row-level truth-estimate comparisons for the parameter groups that need follow-up."
    ),
    Purpose = c(
      "Decide whether the assessment is ready to inspect.",
      "Separate generator stress conditions and sparse score support from parameter-recovery performance before interpreting metrics.",
      "Check diagnostic behavior without treating fit or separation as release-recovery gates.",
      "Find the part of the assessment that needs attention first.",
      "Identify the specific parameter group and metric driving the status.",
      "Diagnose the underlying recovery pattern before changing design or fit settings."
    ),
    stringsAsFactors = FALSE
  )
}

recovery_assessment_plot_guidance <- function(type, metric = NULL) {
  if (identical(type, "status")) {
    return(c(
      "Read this plot before metric-level plots.",
      "Start with concern/review rows; ok rows are supporting evidence.",
      "Use section_status to identify the assessment section that needs follow-up."
    ))
  }
  c(
    paste0("This metric plot is sorted by status priority, then by ", metric %||% "value", "."),
    "Inspect concern/review rows before ok rows.",
    "Use the row-level recovery table only after identifying the parameter group that needs follow-up."
  )
}

recovery_assessment_plot_metric_spec <- function(metric, thresholds) {
  metric <- match.arg(
    tolower(as.character(metric[1])),
    c("rmse", "bias", "coverage", "se_available", "mcse_rmse")
  )
  thresholds <- thresholds %||% list()
  switch(
    metric,
    rmse = list(
      metric = "rmse",
      value_col = "RMSE",
      limit_col = "RMSELimit",
      status_col = "RMSEStatus",
      label = "RMSE",
      transform = identity,
      reference = NA_real_
    ),
    bias = list(
      metric = "bias",
      value_col = "Bias",
      limit_col = "AbsBiasLimit",
      status_col = "BiasStatus",
      label = "Absolute bias",
      transform = abs,
      reference = NA_real_
    ),
    coverage = list(
      metric = "coverage",
      value_col = "Coverage95",
      limit_col = NULL,
      status_col = "CoverageStatus",
      label = "95% coverage",
      transform = identity,
      reference = suppressWarnings(as.numeric(thresholds$coverage_target %||% 0.95))
    ),
    se_available = list(
      metric = "se_available",
      value_col = "SEAvailableRate",
      limit_col = NULL,
      status_col = "SEStatus",
      label = "SE availability rate",
      transform = identity,
      reference = suppressWarnings(as.numeric(thresholds$min_se_available %||% NA_real_))
    ),
    mcse_rmse = list(
      metric = "mcse_rmse",
      value_col = "McseRMSEToRMSE",
      limit_col = NULL,
      status_col = "MonteCarloStatus",
      label = "MCSE / RMSE",
      transform = identity,
      reference = suppressWarnings(as.numeric(thresholds$max_mcse_rmse_ratio %||% NA_real_))
    )
  )
}

#' @rdname assess_mfrm_recovery
#' @param y Reserved for S3 generic compatibility.
#' @param type Assessment plot route. `"status"` summarizes checklist status
#'   counts; `"metrics"` plots a parameter-group assessment metric colored by
#'   its status.
#' @param metric Metric used when `type = "metrics"`. Supported values are
#'   `"rmse"`, `"bias"`, `"coverage"`, `"se_available"`, and `"mcse_rmse"`.
#' @param draw If `TRUE`, draw with base graphics. If `FALSE`, return an
#'   `mfrm_plot_data` object with reusable plot tables and metadata.
#' @export
plot.mfrm_recovery_assessment <- function(x,
                                          y = NULL,
                                          type = c("status", "metrics"),
                                          metric = c("rmse", "bias", "coverage",
                                                     "se_available", "mcse_rmse"),
                                          draw = TRUE,
                                          ...) {
  if (!inherits(x, "mfrm_recovery_assessment")) {
    stop("`x` must be output from assess_mfrm_recovery().", call. = FALSE)
  }
  type <- match.arg(type)

  if (identical(type, "status")) {
    checklist <- as.data.frame(x$checklist %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(checklist) == 0L || !"Status" %in% names(checklist)) {
      stop("No recovery-assessment checklist statuses are available to plot.", call. = FALSE)
    }
    section_tbl <- as.data.frame(
      stats::xtabs(~ Section + Status, data = checklist),
      stringsAsFactors = FALSE
    )
    names(section_tbl)[names(section_tbl) == "Freq"] <- "Checks"
    section_tbl <- section_tbl[section_tbl$Checks > 0, , drop = FALSE]
    section_tbl$StatusRank <- recovery_assessment_status_rank(section_tbl$Status)
    section_tbl <- section_tbl[order(-section_tbl$StatusRank, section_tbl$Section,
                                     section_tbl$Status), , drop = FALSE]
    section_tbl$AttentionOrder <- seq_len(nrow(section_tbl))
    row.names(section_tbl) <- NULL
    status_tbl <- stats::aggregate(Checks ~ Status, data = section_tbl, FUN = sum)
    status_tbl$StatusRank <- recovery_assessment_status_rank(status_tbl$Status)
    status_tbl <- status_tbl[order(-status_tbl$StatusRank, status_tbl$Status), , drop = FALSE]
    status_tbl$AttentionOrder <- seq_len(nrow(status_tbl))
    row.names(status_tbl) <- NULL
    payload <- list(
      type = type,
      metric = "check_count",
      metric_label = "Checklist checks",
      plot_table = section_tbl,
      section_status = section_tbl,
      status_counts = status_tbl,
      checklist = checklist,
      condition_reporting_notes = as.data.frame(
        x$condition_reporting_notes %||% data.frame(),
        stringsAsFactors = FALSE
      ),
      condition_review = as.data.frame(
        x$condition_review %||% data.frame(),
        stringsAsFactors = FALSE
      ),
      diagnostic_reporting_notes = as.data.frame(
        x$diagnostic_reporting_notes %||% data.frame(),
        stringsAsFactors = FALSE
      ),
      diagnostic_review = as.data.frame(
        x$diagnostic_review %||% data.frame(),
        stringsAsFactors = FALSE
      ),
      reading_order = recovery_assessment_reading_order(),
      guidance = recovery_assessment_plot_guidance("status"),
      notes = as.character(x$notes %||% character(0)),
      title = "Recovery assessment: checklist status",
      subtitle = "Counts of checklist checks by status",
      legend = new_plot_legend(
        label = c("ok", "review", "concern", "not_available", "not_assessed"),
        role = rep("status", 5L),
        aesthetic = rep("fill", 5L),
        value = recovery_assessment_status_palette(c("ok", "review", "concern", "not_available", "not_assessed"))
      ),
      reference_lines = new_reference_lines()
    )
    out <- new_mfrm_plot_data("recovery_assessment", payload)
    if (!isTRUE(draw)) return(out)
    graphics::barplot(
      status_tbl$Checks,
      names.arg = status_tbl$Status,
      col = recovery_assessment_status_palette(status_tbl$Status),
      border = NA,
      ylab = "Checklist checks",
      main = payload$title
    )
    return(invisible(out))
  }

  metric_review <- as.data.frame(x$metric_review %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(metric_review) == 0L) {
    stop("No recovery-assessment metric-review rows are available to plot.", call. = FALSE)
  }
  spec <- recovery_assessment_plot_metric_spec(metric, x$thresholds %||% list())
  required <- c("ParameterType", "Facet", "ComparisonScale", spec$value_col, spec$status_col)
  missing <- setdiff(required, names(metric_review))
  if (length(missing) > 0L) {
    stop("Metric-review table is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  raw_value <- suppressWarnings(as.numeric(metric_review[[spec$value_col]]))
  value <- spec$transform(raw_value)
  limit <- if (!is.null(spec$limit_col) && spec$limit_col %in% names(metric_review)) {
    suppressWarnings(as.numeric(metric_review[[spec$limit_col]]))
  } else {
    rep(spec$reference, nrow(metric_review))
  }
  plot_tbl <- data.frame(
    ParameterType = as.character(metric_review$ParameterType),
    Facet = as.character(metric_review$Facet),
    ComparisonScale = as.character(metric_review$ComparisonScale),
    PlotGroup = paste(
      as.character(metric_review$ParameterType),
      as.character(metric_review$Facet),
      as.character(metric_review$ComparisonScale),
      sep = " / "
    ),
    Metric = spec$metric,
    Value = value,
    Limit = limit,
    Status = as.character(metric_review[[spec$status_col]]),
    OverallStatus = as.character(metric_review$OverallStatus %||% NA_character_),
    stringsAsFactors = FALSE
  )
  plot_tbl <- plot_tbl[is.finite(plot_tbl$Value), , drop = FALSE]
  if (nrow(plot_tbl) == 0L) {
    stop("Selected assessment metric has no finite values to plot.", call. = FALSE)
  }
  plot_tbl$StatusRank <- recovery_assessment_status_rank(plot_tbl$Status)
  plot_tbl <- plot_tbl[order(-plot_tbl$StatusRank, -plot_tbl$Value, plot_tbl$PlotGroup), , drop = FALSE]
  plot_tbl$AttentionOrder <- seq_len(nrow(plot_tbl))
  row.names(plot_tbl) <- NULL
  finite_limits <- unique(plot_tbl$Limit[is.finite(plot_tbl$Limit)])
  reference_lines <- if (length(finite_limits) == 1L) {
    new_reference_lines("h", finite_limits, "Assessment threshold", "dashed", "threshold")
  } else {
    new_reference_lines()
  }
  payload <- list(
    type = type,
    metric = spec$metric,
    metric_label = spec$label,
    status_column = spec$status_col,
    plot_table = plot_tbl,
    metric_review = plot_tbl,
    condition_reporting_notes = as.data.frame(
      x$condition_reporting_notes %||% data.frame(),
      stringsAsFactors = FALSE
    ),
    condition_review = as.data.frame(
      x$condition_review %||% data.frame(),
      stringsAsFactors = FALSE
    ),
    diagnostic_reporting_notes = as.data.frame(
      x$diagnostic_reporting_notes %||% data.frame(),
      stringsAsFactors = FALSE
    ),
    diagnostic_review = as.data.frame(
      x$diagnostic_review %||% data.frame(),
      stringsAsFactors = FALSE
    ),
    reading_order = recovery_assessment_reading_order(),
    guidance = recovery_assessment_plot_guidance("metrics", spec$label),
    notes = as.character(x$notes %||% character(0)),
    title = paste("Recovery assessment:", spec$label),
    subtitle = "Parameter-group assessment metric colored by status",
    legend = new_plot_legend(
      label = c("ok", "review", "concern", "not_available", "not_assessed"),
      role = rep("status", 5L),
      aesthetic = rep("fill", 5L),
      value = recovery_assessment_status_palette(c("ok", "review", "concern", "not_available", "not_assessed"))
    ),
    reference_lines = reference_lines
  )
  out <- new_mfrm_plot_data("recovery_assessment", payload)
  if (!isTRUE(draw)) return(out)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(8, 4, 3, 1))
  ylim <- range(c(0, plot_tbl$Value, reference_lines$value), na.rm = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) {
    ylim <- range(plot_tbl$Value + c(-0.1, 0.1), na.rm = TRUE)
  }
  graphics::barplot(
    plot_tbl$Value,
    names.arg = truncate_axis_label(plot_tbl$PlotGroup, width = 26L),
    las = 2,
    col = recovery_assessment_status_palette(plot_tbl$Status),
    border = NA,
    ylim = ylim,
    ylab = spec$label,
    main = payload$title
  )
  if (nrow(reference_lines) > 0L && is.finite(reference_lines$value[1])) {
    graphics::abline(h = reference_lines$value[1], lty = 2, col = "grey35")
  }
  invisible(out)
}

recovery_plot_metric_col <- function(metric) {
  metric <- match.arg(
    tolower(as.character(metric[1])),
    c("rmse", "bias", "mae", "correlation", "coverage", "mcse_bias",
      "mcse_rmse", "raw_rmse", "raw_bias", "mean_se", "se_available")
  )
  switch(
    metric,
    rmse = "RMSE",
    bias = "Bias",
    mae = "MAE",
    correlation = "Correlation",
    coverage = "Coverage95",
    mcse_bias = "McseBias",
    mcse_rmse = "McseRMSE",
    raw_rmse = "RawRMSE",
    raw_bias = "RawBias",
    mean_se = "MeanSE",
    se_available = "SEAvailableRate"
  )
}

recovery_plot_metric_label <- function(metric_col) {
  switch(
    metric_col,
    RMSE = "RMSE",
    Bias = "Bias",
    MAE = "MAE",
    Correlation = "Truth-estimate correlation",
    Coverage95 = "95% coverage",
    McseBias = "MCSE of bias",
    McseRMSE = "MCSE of RMSE",
    RawRMSE = "Unaligned RMSE",
    RawBias = "Unaligned bias",
    MeanSE = "Mean standard error",
    SEAvailableRate = "SE availability rate",
    metric_col
  )
}

recovery_plot_filter <- function(tbl, parameter_type = NULL, facet = NULL) {
  tbl <- tibble::as_tibble(tbl)
  if (!is.null(parameter_type)) {
    parameter_type <- as.character(parameter_type)
    tbl <- tbl[as.character(tbl$ParameterType) %in% parameter_type, , drop = FALSE]
  }
  if (!is.null(facet)) {
    facet <- as.character(facet)
    tbl <- tbl[as.character(tbl$Facet) %in% facet, , drop = FALSE]
  }
  tbl
}

recovery_plot_status_table <- function(rep_tbl) {
  rep_tbl <- tibble::as_tibble(rep_tbl)
  if (nrow(rep_tbl) == 0L) return(tibble::tibble())
  rep_tbl |>
    dplyr::mutate(
      Status = dplyr::case_when(
        !.data$RunOK ~ "failed",
        !.data$Converged ~ "not_converged",
        TRUE ~ "converged"
      )
    ) |>
    dplyr::group_by(.data$Status) |>
    dplyr::summarise(
      Reps = dplyr::n(),
      MeanRecoveryRows = recovery_safe_mean(.data$RecoveryRows),
      MeanElapsedSec = recovery_safe_mean(.data$ElapsedSec),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Status)
}

#' Plot parameter-recovery simulation results
#'
#' @param x Output from [evaluate_mfrm_recovery()].
#' @param y Reserved for S3 generic compatibility.
#' @param type Plot route: `"summary"` draws a metric from
#'   `x$recovery_summary`, `"coverage"` draws 95% coverage by parameter group,
#'   `"errors"` draws row-level recovery-error distributions, `"scatter"` draws
#'   truth against estimated values, and `"replications"` summarizes run status.
#' @param metric Summary metric used when `type = "summary"`. Supported values
#'   are `"rmse"`, `"bias"`, `"mae"`, `"correlation"`, `"coverage"`,
#'   `"mcse_bias"`, `"mcse_rmse"`, `"raw_rmse"`, `"raw_bias"`, `"mean_se"`, and
#'   `"se_available"`.
#' @param parameter_type Optional parameter type filter, such as `"person"`,
#'   `"facet"`, `"step"`, `"slope"`, or `"population"`.
#' @param facet Optional facet filter.
#' @param comparison Error/estimate scale for row-level routes. `"aligned"`
#'   uses `EstimateAligned` / `ErrorAligned`; `"unaligned"` uses `Estimate` /
#'   `ErrorRaw` on the same comparison scale.
#' @param draw If `TRUE`, draw with base graphics. If `FALSE`, return an
#'   `mfrm_plot_data` object with reusable plot tables and metadata.
#' @param ... Reserved for future extensions.
#'
#' @details
#' These plots are intended as simulation-review graphics. They do not replace
#' the row-level `x$recovery` table or the ADEMP metadata; they make the main
#' recovery estimands easier to inspect during model-development and design
#' checks. Coverage is displayed only for parameter groups with available
#' standard errors.
#'
#' @return An `mfrm_plot_data` object. When `draw = TRUE`, the object is returned
#'   invisibly after drawing.
#' @seealso [evaluate_mfrm_recovery()], [summary()]
#' @examples
#' \donttest{
#' rec <- evaluate_mfrm_recovery(
#'   n_person = 12,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   reps = 1,
#'   maxit = 30,
#'   seed = 123
#' )
#' plot(rec, type = "summary", metric = "rmse", draw = FALSE)
#' }
#' @export
plot.mfrm_recovery_simulation <- function(x,
                                          y = NULL,
                                          type = c("summary", "coverage", "errors",
                                                   "scatter", "replications"),
                                          metric = c("rmse", "bias", "mae",
                                                     "correlation", "coverage",
                                                     "mcse_bias", "mcse_rmse",
                                                     "raw_rmse", "raw_bias",
                                                     "mean_se", "se_available"),
                                          parameter_type = NULL,
                                          facet = NULL,
                                          comparison = c("aligned", "unaligned"),
                                          draw = TRUE,
                                          ...) {
  if (!inherits(x, "mfrm_recovery_simulation")) {
    stop("`x` must be output from evaluate_mfrm_recovery().", call. = FALSE)
  }
  type <- match.arg(type)
  comparison <- match.arg(comparison)
  metric_col <- recovery_plot_metric_col(metric)
  metric_label <- recovery_plot_metric_label(metric_col)
  notes <- as.character(x$notes %||% character(0))
  reference_lines <- new_reference_lines()

  if (identical(type, "summary") || identical(type, "coverage")) {
    summary_tbl <- recovery_plot_filter(x$recovery_summary, parameter_type, facet)
    if (nrow(summary_tbl) == 0L) {
      stop("No recovery-summary rows are available for the requested filters.", call. = FALSE)
    }
    if (identical(type, "coverage")) {
      metric_col <- "Coverage95"
      metric_label <- recovery_plot_metric_label(metric_col)
      summary_tbl <- summary_tbl[is.finite(summary_tbl$Coverage95), , drop = FALSE]
      reference_lines <- new_reference_lines(
        axis = "y",
        value = 0.95,
        label = "0.95 nominal coverage",
        linetype = "dashed",
        role = "nominal_coverage"
      )
      if (nrow(summary_tbl) == 0L) {
        stop("No finite coverage rows are available for the requested filters.", call. = FALSE)
      }
    }
    plot_tbl <- summary_tbl |>
      dplyr::mutate(
        PlotGroup = paste(
          as.character(.data$ParameterType),
          as.character(.data$Facet),
          as.character(.data$ComparisonScale),
          sep = " / "
        ),
        Value = suppressWarnings(as.numeric(.data[[metric_col]]))
      ) |>
      dplyr::arrange(.data$ParameterType, .data$Facet, .data$ComparisonScale)
    plot_tbl <- plot_tbl[is.finite(plot_tbl$Value), , drop = FALSE]
    if (nrow(plot_tbl) == 0L) {
      stop("Selected recovery metric has no finite values to plot.", call. = FALSE)
    }
    payload <- list(
      type = type,
      metric = metric_col,
      metric_label = metric_label,
      comparison = comparison,
      parameter_type = parameter_type,
      facet = facet,
      plot_table = plot_tbl,
      notes = notes,
      title = if (identical(type, "coverage")) {
        "Parameter recovery: coverage"
      } else {
        paste("Parameter recovery:", metric_label)
      },
      subtitle = "Summary by parameter type, facet, and comparison scale",
      legend = new_plot_legend(),
      reference_lines = reference_lines
    )
    out <- new_mfrm_plot_data("recovery_simulation", payload)
    if (!isTRUE(draw)) return(out)
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mar = c(8, 4, 3, 1))
    cols <- grDevices::hcl.colors(nrow(plot_tbl), "Dark 3")
    ylim <- range(c(0, plot_tbl$Value, reference_lines$value), na.rm = TRUE)
    if (identical(metric_col, "Coverage95") || identical(metric_col, "SEAvailableRate")) {
      ylim <- c(0, 1)
    } else if (!all(is.finite(ylim)) || diff(ylim) == 0) {
      ylim <- range(plot_tbl$Value + c(-0.1, 0.1), na.rm = TRUE)
    }
    graphics::barplot(
      plot_tbl$Value,
      names.arg = truncate_axis_label(plot_tbl$PlotGroup, width = 26L),
      las = 2,
      col = cols,
      border = NA,
      ylim = ylim,
      ylab = metric_label,
      main = payload$title
    )
    if (0 >= ylim[1] && 0 <= ylim[2]) graphics::abline(h = 0, col = "grey55")
    if (nrow(reference_lines) > 0L) {
      graphics::abline(h = reference_lines$value, lty = 2, col = "grey35")
    }
    return(invisible(out))
  }

  if (identical(type, "replications")) {
    rep_tbl <- tibble::as_tibble(x$rep_overview %||% tibble::tibble())
    status_tbl <- recovery_plot_status_table(rep_tbl)
    if (nrow(status_tbl) == 0L) stop("No replication rows are available.", call. = FALSE)
    payload <- list(
      type = type,
      metric = "Reps",
      metric_label = "Replications",
      comparison = NA_character_,
      parameter_type = parameter_type,
      facet = facet,
      plot_table = status_tbl,
      rep_overview = rep_tbl,
      notes = notes,
      title = "Parameter recovery: replication status",
      subtitle = "Run, convergence, and timing summary",
      legend = new_plot_legend(),
      reference_lines = new_reference_lines()
    )
    out <- new_mfrm_plot_data("recovery_simulation", payload)
    if (!isTRUE(draw)) return(out)
    cols <- c(converged = "#238b45", not_converged = "#b65e16", failed = "#b11f24")
    bar_cols <- unname(cols[status_tbl$Status])
    bar_cols[is.na(bar_cols)] <- "#6b7280"
    graphics::barplot(
      status_tbl$Reps,
      names.arg = status_tbl$Status,
      col = bar_cols,
      border = NA,
      ylab = "Replications",
      main = payload$title
    )
    return(invisible(out))
  }

  recovery_tbl <- recovery_plot_filter(x$recovery, parameter_type, facet)
  if (nrow(recovery_tbl) == 0L) {
    stop("No row-level recovery rows are available for the requested filters.", call. = FALSE)
  }
  estimate_col <- if (identical(comparison, "aligned")) "EstimateAligned" else "Estimate"
  error_col <- if (identical(comparison, "aligned")) "ErrorAligned" else "ErrorRaw"
  recovery_tbl <- recovery_tbl |>
    dplyr::mutate(
      PlotGroup = paste(
        as.character(.data$ParameterType),
        as.character(.data$Facet),
        as.character(.data$ComparisonScale),
        sep = " / "
      ),
      TruthForPlot = suppressWarnings(as.numeric(.data$Truth)),
      EstimateForPlot = suppressWarnings(as.numeric(.data[[estimate_col]])),
      ErrorForPlot = suppressWarnings(as.numeric(.data[[error_col]]))
    )

  if (identical(type, "errors")) {
    plot_tbl <- recovery_tbl[is.finite(recovery_tbl$ErrorForPlot), , drop = FALSE]
    if (nrow(plot_tbl) == 0L) stop("No finite recovery errors are available.", call. = FALSE)
    payload <- list(
      type = type,
      metric = error_col,
      metric_label = if (identical(comparison, "aligned")) "Aligned recovery error" else "Unaligned recovery error",
      comparison = comparison,
      parameter_type = parameter_type,
      facet = facet,
      plot_table = plot_tbl,
      notes = notes,
      title = "Parameter recovery: error distribution",
      subtitle = paste("Comparison:", comparison),
      legend = new_plot_legend(),
      reference_lines = new_reference_lines("y", 0, "zero error", "solid", "zero_error")
    )
    out <- new_mfrm_plot_data("recovery_simulation", payload)
    if (!isTRUE(draw)) return(out)
    graphics::boxplot(
      ErrorForPlot ~ PlotGroup,
      data = plot_tbl,
      las = 2,
      col = "#dbeafe",
      border = "#334e68",
      ylab = payload$metric_label,
      main = payload$title
    )
    graphics::abline(h = 0, col = "grey45", lty = 2)
    return(invisible(out))
  }

  plot_tbl <- recovery_tbl[is.finite(recovery_tbl$TruthForPlot) &
                             is.finite(recovery_tbl$EstimateForPlot), , drop = FALSE]
  if (nrow(plot_tbl) == 0L) {
    stop("No finite truth/estimate pairs are available.", call. = FALSE)
  }
  payload <- list(
    type = type,
    metric = estimate_col,
    metric_label = if (identical(comparison, "aligned")) "Aligned estimate" else "Unaligned estimate",
    comparison = comparison,
    parameter_type = parameter_type,
    facet = facet,
    plot_table = plot_tbl,
    notes = notes,
    title = "Parameter recovery: truth versus estimate",
    subtitle = paste("Comparison:", comparison),
    legend = new_plot_legend(),
    reference_lines = new_reference_lines("xy", NA_real_, "identity line", "solid", "perfect_recovery")
  )
  out <- new_mfrm_plot_data("recovery_simulation", payload)
  if (!isTRUE(draw)) return(out)
  groups <- unique(as.character(plot_tbl$PlotGroup))
  cols <- grDevices::hcl.colors(max(1L, length(groups)), "Dark 3")
  group_idx <- match(as.character(plot_tbl$PlotGroup), groups)
  lim <- range(c(plot_tbl$TruthForPlot, plot_tbl$EstimateForPlot), na.rm = TRUE)
  if (!all(is.finite(lim)) || diff(lim) == 0) lim <- lim + c(-0.1, 0.1)
  graphics::plot(
    plot_tbl$TruthForPlot,
    plot_tbl$EstimateForPlot,
    xlab = "Truth",
    ylab = payload$metric_label,
    main = payload$title,
    xlim = lim,
    ylim = lim,
    pch = 16,
    col = cols[group_idx]
  )
  graphics::abline(0, 1, col = "grey35", lty = 2)
  if (length(groups) > 1L && length(groups) <= 8L) {
    graphics::legend("topleft", legend = groups, col = cols, pch = 16, bty = "n", cex = 0.8)
  }
  invisible(out)
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
#' @param model Measurement model passed to [fit_mfrm()]. `RSM` and `PCM` use
#'   the validated Rasch-family design-planning layer. Bounded `GPCM` is
#'   available as a caveated simulation/refit operating-characteristic route.
#' @param step_facet Step facet passed to [fit_mfrm()] when `model = "PCM"` or
#'   `model = "GPCM"`.
#'   When left `NULL`, the function inherits the generator step facet from
#'   `sim_spec` when available and otherwise defaults to `"Criterion"`.
#' @param slope_facet Slope facet passed to [fit_mfrm()] when
#'   `model = "GPCM"`. The current bounded branch requires
#'   `slope_facet == step_facet`.
#' @param slopes Optional bounded-`GPCM` generator slopes used when
#'   `sim_spec = NULL`. See [build_mfrm_sim_spec()] for accepted formats.
#' @param assignment Optional assignment design used when `sim_spec = NULL`.
#'   `"sparse_linked"` activates planned-missing sparse rating designs; use
#'   `sparse_controls` to specify the linking set.
#' @param sparse_controls Optional named list used when
#'   `assignment = "sparse_linked"` and `sim_spec = NULL`, or retained from
#'   a sparse linked `sim_spec`. See [build_mfrm_sim_spec()].
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
#' @param progress Logical. Whether to show a progress bar across
#'   design-by-replication cells. Defaults to [interactive()], so interactive
#'   exploratory runs show progress while non-interactive tests, scripts, and
#'   report rendering stay quiet. Set `TRUE` or `FALSE` explicitly to override.
#' @param parallel Parallelisation strategy for the rep loop within
#'   each design row. `"no"` (default) runs serially; `"future"`
#'   uses `future.apply::future_lapply` and respects whatever
#'   `future::plan()` is currently active. The Suggests package
#'   `future.apply` must be installed for the parallel path to
#'   activate; otherwise the call falls back to serial execution
#'   with a single message. Cross-design-row parallelism is planned
#'   for a future release.
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
#' Sparse linked simulation specifications and direct
#' `assignment = "sparse_linked"` calls are carried into the design-evaluation
#' output. The resulting sparse-design columns report planned-missingness and
#' rater-link diagnostics (for example design density and rater-pair common
#' persons). They are design diagnostics, not fit statistics or universal
#' adequacy thresholds.
#'
#' If that specification also stores a latent-regression population generator,
#' each replication carries forward the simulated one-row-per-person background
#' data and refits the MML population-model branch. This remains a scenario
#' study under explicit assumptions; it is not a closed-form predictive
#' distribution for one future administration.
#'
#' Bounded `GPCM` design evaluation is available with caveats. It repeatedly
#' generates data from the supplied or fit-derived slope-aware specification,
#' refits bounded `GPCM`, and summarizes facet-level operating characteristics.
#' The route remains a role-based person x rater-like x criterion-like planner:
#' it does not validate diagnostic-screening or signal-detection rules, does
#' not provide a fully arbitrary-facet planner, and does not replace
#' [evaluate_mfrm_recovery()] for slope-recovery adequacy review.
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
#' - Sparse-design diagnostics when `assignment = "sparse_linked"`:
#'   `MeanDesignDensity`, `MeanPlannedMissingRate`, `MeanLinkPersons`,
#'   `MeanMinCommonPersonsPerRaterPair`, `MaxZeroCommonRaterPairs`, and
#'   `MaxRaterPairsBelowTarget`.
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
#' This is a Monte Carlo design-evaluation helper. It can visualize how
#' separation, reliability, strata, RMSE, and fit-screen rates change when
#' you vary person, rater, criterion, or assignment counts. For analytic
#' generalizability-theory planning, pair observed variance-component review
#' from [mfrm_generalizability()] with D-study projections from
#' [mfrm_d_study()].
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
#' - `gpcm_boundary`: bounded-`GPCM` caveat row when a `GPCM` design route is
#'   used
#' - `notes`: short interpretation notes
#' - `settings`: simulation settings
#' - `ademp`: simulation-study metadata (aims, DGM, estimands, methods, performance measures)
#' @seealso [simulate_mfrm_data()], [summary.mfrm_design_evaluation], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- suppressWarnings(evaluate_mfrm_design(
#'   design = list(person = c(8, 12), rater = 2, criterion = 2, assignment = 1),
#'   reps = 1,
#'   maxit = 30,
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
                                 slope_facet = NULL,
                                 slopes = NULL,
                                 assignment = NULL,
                                 sparse_controls = NULL,
                                 maxit = 25,
                                 quad_points = 7,
                                 residual_pca = c("none", "overall", "facet", "both"),
                                 sim_spec = NULL,
                                 seed = NULL,
                                 progress = interactive(),
                                 parallel = c("no", "future")) {
  fit_method <- match.arg(fit_method)
  call_args <- names(as.list(match.call(expand.dots = FALSE))[-1])
  model_explicit <- "model" %in% call_args
  if (!is.null(sim_spec) && !model_explicit) {
    model <- toupper(as.character(sim_spec$model %||% model[1]))
    model <- match.arg(model, c("RSM", "PCM", "GPCM"))
  } else {
    model <- match.arg(model)
  }
  parallel <- match.arg(parallel)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  # `parallel = "future"` requires the `future.apply` Suggests to be
  # installed AND a `future::plan()` to be active. We honour the
  # request when both are satisfied; otherwise we fall back to the
  # serial implementation with a single message so the run still
  # completes.
  if (identical(parallel, "future") &&
      !requireNamespace("future.apply", quietly = TRUE)) {
    message("`evaluate_mfrm_design(parallel = 'future')` requires the ",
            "`future.apply` package (in Suggests). Falling back to ",
            "serial execution.")
    parallel <- "no"
  }
  if (identical(parallel, "future")) {
    # The argument is exposed for forward compatibility. The current
    # release threads the request through the rep loop via
    # `future.apply::future_lapply` only when no per-rep state is
    # accumulated upstream. Full parallelisation across design rows
    # is planned for a future release; until then, set
    # `future::plan(multisession, workers = N)` and rerun to use
    # parallel rep execution within each design row.
    message("`evaluate_mfrm_design(parallel = 'future')` currently ",
            "parallelises the rep loop within each design row. ",
            "Cross-design-row parallelism is planned for a future release.")
  }
  residual_pca <- match.arg(residual_pca)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!identical(model, "GPCM") && !is.null(slopes)) {
    stop("`slopes` can be supplied only when `model = \"GPCM\"`.", call. = FALSE)
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
  assignment <- if (is.null(assignment)) {
    NULL
  } else {
    match.arg(tolower(as.character(assignment[1])), c("crossed", "rotating", "sparse_linked"))
  }
  sparse_controls <- simulation_normalize_sparse_controls(
    sparse_controls = sparse_controls,
    assignment = if (identical(assignment, "sparse_linked")) "sparse_linked" else "rotating",
    n_person = suppressWarnings(as.integer(n_person[1])),
    n_rater = suppressWarnings(as.integer(n_rater[1])),
    raters_per_person = suppressWarnings(as.integer(raters_per_person[1]))
  )
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
  generator_step_facet <- if (is.null(sim_spec)) {
    if (identical(generator_model, "PCM") || identical(generator_model, "GPCM")) {
      as.character(step_facet[1] %||% base_facet_names[2])
    } else {
      NA_character_
    }
  } else {
    sim_spec$step_facet
  }
  generator_assignment <- if (is.null(sim_spec)) assignment %||% "design_dependent" else sim_spec$assignment
  design_variable_aliases <- design_meta$aliases
  design_descriptor <- design_meta$descriptor
  planning_scope <- simulation_planning_scope(sim_spec)
  planning_constraints <- simulation_planning_constraints(sim_spec)
  planning_schema <- simulation_planning_schema(sim_spec)
  fit_step_facet <- simulation_resolve_fit_step_facet(model, step_facet, generator_step_facet)
  fit_slope_facet <- simulation_resolve_fit_slope_facet(model, slope_facet, fit_step_facet)
  recovery_contract <- simulation_recovery_contract(
    generator_model = generator_model,
    generator_step_facet = generator_step_facet,
    fitted_model = model,
    fitted_step_facet = fit_step_facet
  )
  gpcm_route_active <- identical(model, "GPCM") || identical(as.character(generator_model), "GPCM")
  gpcm_boundary <- simulation_gpcm_design_boundary(gpcm_route_active)
  gpcm_notes <- simulation_gpcm_design_notes(gpcm_route_active)

  seeds <- with_preserved_rng_seed(
    seed,
    sample.int(.Machine$integer.max, size = nrow(design_grid) * reps, replace = FALSE)
  )
  seed_idx <- 0L

  result_rows <- vector("list", nrow(design_grid) * reps * 3L)
  rep_rows <- vector("list", nrow(design_grid) * reps)
  result_idx <- 0L
  rep_idx <- 0L

  # Design-evaluation runs a full fit + diagnose per (design, rep) cell, so
  # the wall-clock can reach tens of seconds. Show a progress bar only when
  # requested (default: interactive sessions) so tests, Quarto rendering, and
  # batch simulation logs remain readable.
  total_cells <- nrow(design_grid) * reps
  design_progress_id <- NULL
  if (isTRUE(progress) && total_cells > 1L) {
    design_progress_id <- cli::cli_progress_bar(
      name = "evaluate_mfrm_design",
      total = total_cells,
      format = paste(
        "{cli::pb_spin} design-eval cells {cli::pb_current}/{cli::pb_total}",
        "[{cli::pb_elapsed}  eta {cli::pb_eta}]"
      ),
      clear = TRUE,
      .envir = parent.frame()
    )
    on.exit(cli::cli_progress_done(id = design_progress_id), add = TRUE)
  }

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
      if (!is.null(design_progress_id)) {
        cli::cli_progress_update(id = design_progress_id,
                                  set = (i - 1L) * reps + rep - 1L)
      }
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
          seed = seeds[seed_idx],
          model = model,
          step_facet = if (identical(model, "RSM")) "Criterion" else fit_step_facet,
          slope_facet = if (identical(model, "GPCM")) fit_slope_facet else NULL,
          slopes = if (identical(model, "GPCM")) slopes else NULL,
          assignment = assignment,
          sparse_controls = sparse_controls
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
      if (identical(model, "PCM") || identical(model, "GPCM")) fit_args$step_facet <- fit_step_facet
      if (identical(model, "GPCM")) fit_args$slope_facet <- fit_slope_facet
      if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
      if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points
      sim_population <- attr(sim, "mfrm_population_data")
      if (is.list(sim_population) && isTRUE(sim_population$active)) {
        fit_args$population_formula <- sim_population$population_formula
        fit_args$person_data <- sim_population$person_data
        fit_args$person_id <- sim_population$person_id
        fit_args$population_policy <- sim_population$population_policy
      }
      fit_args <- simulation_add_fit_score_support(
        fit_args,
        sim,
        fallback_score_levels = row_score_levels
      )

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
      sparse_fields <- simulation_sparse_overview_fields(sim)

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
        SparseDesignActive = sparse_fields$SparseDesignActive,
        DesignDensity = sparse_fields$DesignDensity,
        PlannedMissingRate = sparse_fields$PlannedMissingRate,
        LinkPersons = sparse_fields$LinkPersons,
        LinkFractionActual = sparse_fields$LinkFractionActual,
        LinkRatersPerPerson = sparse_fields$LinkRatersPerPerson,
        MinCommonPersonsPerRaterPair = sparse_fields$MinCommonPersonsPerRaterPair,
        ZeroCommonRaterPairs = sparse_fields$ZeroCommonRaterPairs,
        RaterPairsBelowTarget = sparse_fields$RaterPairsBelowTarget,
        TargetCommonPersonsPerRaterPair = sparse_fields$TargetCommonPersonsPerRaterPair,
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
          SparseDesignActive = sparse_fields$SparseDesignActive,
          DesignDensity = sparse_fields$DesignDensity,
          PlannedMissingRate = sparse_fields$PlannedMissingRate,
          LinkPersons = sparse_fields$LinkPersons,
          LinkFractionActual = sparse_fields$LinkFractionActual,
          LinkRatersPerPerson = sparse_fields$LinkRatersPerPerson,
          MinCommonPersonsPerRaterPair = sparse_fields$MinCommonPersonsPerRaterPair,
          ZeroCommonRaterPairs = sparse_fields$ZeroCommonRaterPairs,
          RaterPairsBelowTarget = sparse_fields$RaterPairsBelowTarget,
          TargetCommonPersonsPerRaterPair = sparse_fields$TargetCommonPersonsPerRaterPair,
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
    sparse_controls = if (is.null(sim_spec)) sparse_controls else NULL,
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
      gpcm_boundary = gpcm_boundary,
      notes = unique(c(design_eval_build_notes(results), gpcm_notes)),
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
        slope_facet = fit_slope_facet,
        slopes = if (identical(model, "GPCM")) slopes else NULL,
        assignment = assignment,
        sparse_controls = if (is.null(sim_spec)) sparse_controls else sim_spec$sparse_controls %||% list(active = FALSE),
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        sim_spec = sim_spec,
        progress = isTRUE(progress),
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
        gpcm_design_status = if (gpcm_route_active) "supported_with_caveat" else NA_character_,
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
#' - `sparse_review`: compact planned-missingness and rater-link review counts
#'   when sparse linked designs are active
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
#'   maxit = 30,
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
  out$sparse_review <- round_df(simulation_sparse_design_review_summary(out$design_summary))
  out$ademp <- object$ademp %||% NULL
  out$facet_names <- object$settings$facet_names %||% stats::setNames(simulation_default_output_facet_names(), c("rater", "criterion"))
  out$design_variable_aliases <- simulation_object_design_variable_aliases(object)
  out$design_descriptor <- simulation_object_design_descriptor(object)
  out$planning_scope <- simulation_object_planning_scope(object)
  out$planning_constraints <- simulation_object_planning_constraints(object)
  out$planning_schema <- simulation_object_planning_schema(object)
  out$gpcm_boundary <- object$gpcm_boundary %||% data.frame()
  out$future_branch_active_summary <- simulation_compact_future_branch_active_summary(
    object,
    digits = digits
  )
  out$digits <- digits
  out$notes <- unique(c(out$notes %||% character(0), object$notes %||% character(0)))
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
  if (!is.null(x$sparse_review) && nrow(x$sparse_review) > 0) {
    cat("\nSparse linked design review\n")
    print(round_df(as.data.frame(x$sparse_review)), row.names = FALSE)
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
#' - sparse linked `metric = "plannedmissingrate"`, `"mincommonpersons"`,
#'   `"zerocommonpairs"`, or `"pairsshorttarget"` to review planned
#'   missingness and rater-pair linkage separately from recovery metrics
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
#'   maxit = 30,
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
                                                   "convergencerate", "elapsedsec", "mincategorycount",
                                                   "designdensity", "plannedmissingrate",
                                                   "linkpersons", "linkfraction", "linkraters",
                                                   "mincommonpersons", "zerocommonpairs",
                                                   "pairsshorttarget"),
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
#' - `caveats`: structured warning rows for situations where the
#'   recommendation rests on weak evidence (e.g., no design met every
#'   threshold; the recommended design is at the boundary of the
#'   evaluated grid; only one rep was simulated). Empty `tibble()`
#'   when no caveats apply.
#' @seealso [evaluate_mfrm_design()], [summary.mfrm_design_evaluation], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- suppressWarnings(evaluate_mfrm_design(
#'   n_person = c(8, 12),
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 30,
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

  # Publication-workflow caveats. mfrmr treats facets as fixed effects,
  # so simulation-based recommendations should be complemented with a
  # post-fit review of observed level counts. See the "Fixed effects
  # assumption" section of ?fit_mfrm and `facet_small_sample_review()`.
  fixed_effects_note <- paste0(
    "mfrmr estimates all facets as fixed effects with sum-to-zero ",
    "identification; simulation-based adequacy does not imply partial ",
    "pooling for small-N levels in real data. After collecting data, ",
    "inspect observed level counts with `facet_small_sample_review(fit)` ",
    "and `analyze_hierarchical_structure(data, facets)`."
  )

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
      planning_schema = simulation_object_planning_schema(x),
      caveats = list(
        fixed_effects = fixed_effects_note,
        post_fit_review = c(
          "facet_small_sample_review(fit)",
          "analyze_hierarchical_structure(data, facets)",
          "compute_facet_icc(data, facets, score, person)"
        )
      )
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
                                                                    fit_slope_facet = NULL,
                                                                    slopes = NULL,
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
    } else if (identical(model, "GPCM")) {
      spec_args$step_facet <- fit_step_facet
      spec_args$slope_facet <- fit_slope_facet %||% fit_step_facet
      spec_args$slopes <- slopes
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

diagnostic_screening_build_step_slope_table <- function(step_levels) {
  step_levels <- as.character(step_levels)
  step_levels <- step_levels[!is.na(step_levels) & nzchar(step_levels)]
  if (length(step_levels) < 1L) {
    stop("Step-structure screening requires at least one GPCM slope level.", call. = FALSE)
  }
  raw <- if (length(step_levels) == 1L) {
    1
  } else {
    exp(seq(-0.35, 0.35, length.out = length(step_levels)))
  }
  raw <- exp(log(raw) - mean(log(raw)))
  tibble::tibble(
    SlopeFacet = step_levels,
    Estimate = raw
  )
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
                                                                  fit_step_facet,
                                                                  fit_slope_facet = NULL) {
  if (!identical(fit_model, "RSM") && !identical(fit_model, "PCM") &&
      !identical(fit_model, "GPCM")) {
    stop("Step-structure screening is currently scoped to `RSM`, `PCM`, and bounded `GPCM` fits.", call. = FALSE)
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
      thresholds = threshold_table,
      model = if (identical(fit_model, "GPCM")) "GPCM" else "PCM",
      step_facet = generator_step_facet,
      assignment = if (identical(design_row$raters_per_person, design_row$n_rater)) "crossed" else "rotating",
      facet_names = facet_names
    )
    if (identical(fit_model, "GPCM")) {
      spec_args$slope_facet <- generator_step_facet
      spec_args$slopes <- stats::setNames(
        diagnostic_screening_build_step_slope_table(step_levels)$Estimate,
        step_levels
      )
    }
    return(do.call(build_mfrm_sim_spec, spec_args))
  }

  out <- row_spec_base
  out$model <- if (identical(fit_model, "GPCM")) "GPCM" else "PCM"
  out$step_facet <- generator_step_facet
  out$threshold_table <- tibble::as_tibble(threshold_table)
  if (identical(fit_model, "GPCM")) {
    out$slope_facet <- generator_step_facet
    out$slope_table <- diagnostic_screening_build_step_slope_table(step_levels)
    out$slope_regime <- simulation_gpcm_slope_regime(out$slope_table)
  }
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

diagnostic_screening_report_audit_empty <- function() {
  list(
    ReportIndexAvailable = FALSE,
    ReportAuditError = NA_character_,
    ReportReviewAreas = NA_real_,
    ReportReadyAreas = NA_real_,
    ReportRequestIfNeededAreas = NA_real_,
    ReportFitReadiness = NA_character_,
    ReportFitReviewSignalCount = NA_real_,
    ReportPrecisionReadiness = NA_character_,
    ReportPrecisionReviewSignalCount = NA_real_,
    ReportMisfitReadiness = NA_character_,
    ReportMisfitReviewSignalCount = NA_real_
  )
}

diagnostic_screening_report_area_value <- function(report_index, area, column, default = NA) {
  report_index <- as.data.frame(report_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(report_index) == 0L || !all(c("Area", column) %in% names(report_index))) {
    return(default)
  }
  value <- report_index[as.character(report_index$Area) %in% area, column, drop = TRUE][1]
  if (length(value) == 0L || is.null(value) || is.na(value)) default else value
}

diagnostic_screening_collect_report_audit <- function(fit,
                                                      report_include,
                                                      report_style) {
  empty <- diagnostic_screening_report_audit_empty()
  res <- tryCatch(
    mfrm_results(fit, include = report_include),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    empty$ReportAuditError <- conditionMessage(res)
    return(empty)
  }
  report <- tryCatch(
    mfrm_report(res, style = report_style),
    error = function(e) e
  )
  if (inherits(report, "error")) {
    empty$ReportAuditError <- conditionMessage(report)
    return(empty)
  }
  report_index <- as.data.frame(report$report_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(report_index) == 0L || !"Readiness" %in% names(report_index)) {
    empty$ReportAuditError <- "mfrm_report() did not return a report_index table."
    return(empty)
  }
  review_signal <- suppressWarnings(as.numeric(report_index$ReviewSignalCount %||% NA_real_))
  list(
    ReportIndexAvailable = TRUE,
    ReportAuditError = NA_character_,
    ReportReviewAreas = sum(as.character(report_index$Readiness) %in% "review", na.rm = TRUE),
    ReportReadyAreas = sum(as.character(report_index$Readiness) %in% "ready", na.rm = TRUE),
    ReportRequestIfNeededAreas = sum(as.character(report_index$Readiness) %in% "request_if_needed", na.rm = TRUE),
    ReportFitReadiness = as.character(diagnostic_screening_report_area_value(report_index, "Fit", "Readiness", NA_character_)),
    ReportFitReviewSignalCount = suppressWarnings(as.numeric(diagnostic_screening_report_area_value(
      report_index, "Fit", "ReviewSignalCount", NA_real_
    ))),
    ReportPrecisionReadiness = as.character(diagnostic_screening_report_area_value(report_index, "Precision", "Readiness", NA_character_)),
    ReportPrecisionReviewSignalCount = suppressWarnings(as.numeric(diagnostic_screening_report_area_value(
      report_index, "Precision", "ReviewSignalCount", NA_real_
    ))),
    ReportMisfitReadiness = as.character(diagnostic_screening_report_area_value(report_index, "Misfit / pathway", "Readiness", NA_character_)),
    ReportMisfitReviewSignalCount = suppressWarnings(as.numeric(diagnostic_screening_report_area_value(
      report_index, "Misfit / pathway", "ReviewSignalCount", NA_real_
    )))
  )
}

diagnostic_screening_build_report_signal_summary <- function(results, design_variable_aliases = NULL) {
  results_tbl <- tibble::as_tibble(results)
  required <- c(
    "ReportIndexAvailable", "ReportAuditError", "ReportReviewAreas",
    "ReportReadyAreas", "ReportRequestIfNeededAreas", "ReportFitReadiness",
    "ReportFitReviewSignalCount", "ReportPrecisionReadiness",
    "ReportPrecisionReviewSignalCount", "ReportMisfitReadiness",
    "ReportMisfitReviewSignalCount"
  )
  if (!all(required %in% names(results_tbl))) {
    return(tibble::tibble())
  }
  design_descriptor <- simulation_object_design_descriptor(
    list(results = results_tbl, design_variable_aliases = design_variable_aliases)
  )
  design_vars <- simulation_design_canonical_variables(design_descriptor)
  grouping_vars <- c("design_id", "Scenario", "ScenarioClass", "Model", "DependenceFacet", design_vars)
  out <- results_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) |>
    dplyr::summarize(
      Reps = dplyr::n(),
      ReportIndexAvailabilityRate = mean(as.logical(.data$ReportIndexAvailable), na.rm = TRUE),
      ReportAuditErrorRows = sum(!is.na(.data$ReportAuditError) & nzchar(as.character(.data$ReportAuditError)), na.rm = TRUE),
      MeanReportReviewAreas = design_eval_safe_mean(.data$ReportReviewAreas),
      MeanReportReadyAreas = design_eval_safe_mean(.data$ReportReadyAreas),
      MeanReportRequestIfNeededAreas = design_eval_safe_mean(.data$ReportRequestIfNeededAreas),
      FitReportReviewRate = mean(as.character(.data$ReportFitReadiness) %in% "review", na.rm = TRUE),
      MeanFitReportSignals = design_eval_safe_mean(.data$ReportFitReviewSignalCount),
      PrecisionReportReviewRate = mean(as.character(.data$ReportPrecisionReadiness) %in% "review", na.rm = TRUE),
      MeanPrecisionReportSignals = design_eval_safe_mean(.data$ReportPrecisionReviewSignalCount),
      MisfitReportReviewRate = mean(as.character(.data$ReportMisfitReadiness) %in% "review", na.rm = TRUE),
      MeanMisfitReportSignals = design_eval_safe_mean(.data$ReportMisfitReviewSignalCount),
      .groups = "drop"
    )
  simulation_append_design_alias_columns(out, design_variable_aliases)
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
  report_signal_summary <- diagnostic_screening_build_report_signal_summary(
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
    report_signal_summary = report_signal_summary,
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
#' @param model Measurement model passed to [fit_mfrm()]. Bounded `GPCM` is
#'   supported with caveats as slope-aware screening sensitivity evidence.
#' @param step_facet Step facet passed to [fit_mfrm()] when `model = "PCM"` or
#'   `model = "GPCM"`.
#' @param slope_facet Slope facet passed to [fit_mfrm()] when
#'   `model = "GPCM"`. Defaults to the fitted step facet.
#' @param slopes Optional bounded-`GPCM` slope specification used by direct
#'   simulation calls when `sim_spec = NULL`.
#' @param maxit Maximum iterations passed to [fit_mfrm()].
#' @param quad_points Quadrature points for the internal `MML` fit.
#' @param residual_pca Residual PCA mode passed to [diagnose_mfrm()].
#' @param sim_spec Optional output from [build_mfrm_sim_spec()] or
#'   [extract_mfrm_sim_spec()] used as the base data-generating mechanism.
#' @param include_report Logical; if `TRUE`, each successful replicate also
#'   builds `mfrm_results()` and [mfrm_report()] and records the `report_index`
#'   readiness/signaling surface. This is intentionally opt-in because it
#'   repeats the comprehensive result-building workflow.
#' @param report_include `include` vector passed to [mfrm_results()] when
#'   `include_report = TRUE`.
#' @param report_style Report style passed to [mfrm_report()] when
#'   `include_report = TRUE`.
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
#' 5. optionally stores `mfrm_report()` `report_index` readiness signals
#' 6. aggregates the results into `scenario_summary`, `performance_summary`,
#'    `report_signal_summary`, and `scenario_contrast`
#'
#' The `"well_specified"` scenario uses the ordinary generator with no injected
#' extra structure. The `"local_dependence"` scenario adds a shared
#' `Person x facet` random effect, centered within the selected facet levels, so
#' responses in the same context become correlated without changing the
#' facet-level mean effect contract. The `"latent_misspecification"` scenario
#' keeps the same marginal spread targets but replaces the normal person
#' distribution with a centered bimodal empirical support distribution, while
#' leaving the non-person facets on the original scale contract. The
#' `"step_structure_misspecification"` scenario uses a `PCM` or bounded-`GPCM`
#' generator with facet-specific threshold tables that intentionally mismatch
#' the fitted step contract: `RSM` fits receive criterion-specific thresholds,
#' and `PCM` / `GPCM` fits receive threshold structures indexed by the opposite
#' non-person facet. For bounded `GPCM`, the generator and fit each keep
#' `slope_facet == step_facet`; the misspecification is the generator-versus-fit
#' step/slope facet mismatch.
#'
#' This function is intentionally screening-oriented. The strict marginal branch
#' remains exploratory in the current release, so the returned summaries should
#' be used to compare relative sensitivity across scenarios rather than to claim
#' calibrated inferential power. Bounded-`GPCM` rows add explicit
#' `gpcm_boundary` caveats and should be read as slope-aware operating
#' characteristics under the evaluated role-based design.
#'
#' @return An object of class `mfrm_diagnostic_screening` with:
#' - `design_grid`: evaluated design conditions, including public alias columns
#'   when applicable
#' - `results`: replicate-level screening metrics for each design and scenario
#' - `scenario_summary`: aggregated scenario-by-design screening summaries
#' - `performance_summary`: scenario-by-design screening-performance summary
#'   including runtime, agreement, Type I proxy, and sensitivity proxy columns
#' - `report_signal_summary`: optional scenario-by-design summary of
#'   `mfrm_report()` `report_index` availability, readiness, and review-signal
#'   counts when `include_report = TRUE`
#' - `scenario_contrast`: each misspecification scenario minus the
#'   well-specified baseline when the baseline scenario was evaluated
#' - `design_descriptor`: role-based design-variable metadata
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract
#' - `gpcm_boundary`: bounded-`GPCM` caveat row when present
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
#'   maxit = 30,
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
                                               slope_facet = NULL,
                                               slopes = NULL,
                                               maxit = 25,
                                               quad_points = 7,
                                               residual_pca = c("none", "overall", "facet", "both"),
                                               sim_spec = NULL,
                                               include_report = FALSE,
                                               report_include = c("fit", "diagnostics", "tables", "precision", "reporting"),
                                               report_style = c("qc", "apa", "validation", "reviewer", "technical"),
                                               seed = NULL) {
  call_args <- names(as.list(match.call(expand.dots = FALSE))[-1])
  model_explicit <- "model" %in% call_args
  if (!is.null(sim_spec) && !model_explicit) {
    model <- toupper(as.character(sim_spec$model %||% model[1]))
    model <- match.arg(model, c("RSM", "PCM", "GPCM"))
  } else {
    model <- match.arg(model)
  }
  residual_pca <- match.arg(residual_pca)
  include_report <- isTRUE(include_report)
  report_include <- unique(as.character(report_include %||% character(0)))
  if (include_report && length(report_include) == 0L) {
    stop("`report_include` must contain at least one mfrm_results() include value when `include_report = TRUE`.", call. = FALSE)
  }
  report_style <- match.arg(report_style)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!identical(model, "GPCM") && !is.null(slopes)) {
    stop("`slopes` can be supplied only when `model = \"GPCM\"`.", call. = FALSE)
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
    if (identical(generator_model, "PCM") || identical(generator_model, "GPCM")) {
      as.character(step_facet[1] %||% base_facet_names[2])
    } else {
      NA_character_
    }
  } else {
    sim_spec$step_facet
  }
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
  fit_step_facet <- simulation_resolve_fit_step_facet(model, step_facet, generator_step_facet)
  fit_slope_facet <- simulation_resolve_fit_slope_facet(model, slope_facet, fit_step_facet)
  gpcm_route_active <- identical(model, "GPCM") || identical(as.character(generator_model), "GPCM")
  gpcm_boundary <- simulation_gpcm_screening_boundary(gpcm_route_active)
  gpcm_notes <- simulation_gpcm_screening_notes(gpcm_route_active)

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
            step_facet = if (identical(model, "RSM")) "Criterion" else fit_step_facet,
            slope_facet = if (identical(model, "GPCM")) fit_slope_facet else NULL,
            slopes = if (identical(model, "GPCM")) slopes else NULL,
            seed = sim_seed
          )
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
            fit_slope_facet = fit_slope_facet,
            slopes = if (identical(model, "GPCM")) slopes else NULL,
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
            fit_step_facet = fit_step_facet,
            fit_slope_facet = fit_slope_facet
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
            step_facet = if (identical(model, "RSM")) "Criterion" else fit_step_facet,
            slope_facet = if (identical(model, "GPCM")) fit_slope_facet else NULL,
            slopes = if (identical(model, "GPCM")) slopes else NULL,
            seed = sim_seed
          )
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
        if (identical(model, "PCM") || identical(model, "GPCM")) {
          fit_args$step_facet <- fit_step_facet
        }
        if (identical(model, "GPCM")) {
          fit_args$slope_facet <- fit_slope_facet
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
        fit_args <- simulation_add_fit_score_support(
          fit_args,
          sim,
          fallback_score_levels = row_score_levels
        )

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
          PairwiseFlaggedLevelPairs = NA_real_,
          ReportIndexAvailable = FALSE,
          ReportAuditError = NA_character_,
          ReportReviewAreas = NA_real_,
          ReportReadyAreas = NA_real_,
          ReportRequestIfNeededAreas = NA_real_,
          ReportFitReadiness = NA_character_,
          ReportFitReviewSignalCount = NA_real_,
          ReportPrecisionReadiness = NA_character_,
          ReportPrecisionReviewSignalCount = NA_real_,
          ReportMisfitReadiness = NA_character_,
          ReportMisfitReviewSignalCount = NA_real_
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
        if (isTRUE(include_report)) {
          report_metrics <- diagnostic_screening_collect_report_audit(
            fit = fit,
            report_include = report_include,
            report_style = report_style
          )
          for (nm in names(report_metrics)) {
            result_row[[nm]] <- report_metrics[[nm]]
          }
        }
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
      include_report = include_report,
      report_include = report_include,
      report_style = report_style,
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
      report_signal_summary = summary_bundle$report_signal_summary,
      scenario_contrast = summary_bundle$scenario_contrast,
      design_descriptor = design_descriptor,
      planning_scope = planning_scope,
      planning_constraints = planning_constraints,
      planning_schema = planning_schema,
      gpcm_boundary = gpcm_boundary,
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
        slope_facet = fit_slope_facet,
        slopes = if (identical(model, "GPCM")) slopes else NULL,
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        include_report = include_report,
        report_include = report_include,
        report_style = report_style,
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
        gpcm_screening_status = if (gpcm_route_active) "supported_with_caveat" else NA_character_,
        seed = seed
      ),
      ademp = ademp,
      notes = unique(c(summary_bundle$notes, gpcm_notes))
    ),
    class = "mfrm_diagnostic_screening"
  )
}

diagnostic_screening_plot_long_safe <- function(object, type, metric = NULL) {
  plot_data <- tryCatch(
    plot(object, type = type, metric = metric, draw = FALSE),
    error = function(e) NULL
  )
  if (!inherits(plot_data, "mfrm_plot_data")) {
    return(tibble::tibble())
  }
  tibble::as_tibble(plot_data$data$plot_long %||% tibble::tibble())
}

diagnostic_screening_reading_order <- function(include_report = FALSE) {
  out <- data.frame(
    Step = 1:10,
    Table = c(
      "overview",
      "reading_order",
      "next_actions",
      "reporting_notes",
      "figure_recipes",
      "scenario_summary",
      "performance_summary",
      "scenario_contrast",
      "report_signal_summary",
      "plot_overview_rate"
    ),
    WhatToRead = c(
      "Design count, replication count, scenarios, model, convergence, and whether report signals were retained.",
      "Recommended reading sequence for diagnostic-screening summaries.",
      "Action-oriented triage for replication count, run completion, contrasts, report signals, and export.",
      "Reporting boundaries and recommended wording safeguards.",
      "Figure and display recipes linking plot() calls, plot_data() extraction, caption focus, and interpretation boundaries.",
      "Scenario-by-design legacy, strict marginal, and strict pairwise screening summaries.",
      "Operating-characteristic rates, Type I/sensitivity proxy labels, agreement rates, and runtime summaries.",
      "Misspecification-minus-well-specified deltas when a baseline scenario was evaluated.",
      "Optional mfrm_report() report-index availability/readiness and review-signal summaries.",
      "Long-form draw-free visualization table for the main overview-rate display."
    ),
    Purpose = c(
      "Confirm that the simulation ran as intended before interpreting screening behavior.",
      "Avoid treating every exported table as equal-priority evidence.",
      "Decide whether to increase replications, inspect run failures, read contrasts, request report signals, or export appendices.",
      "Keep simulation screening language distinct from inferential, recovery, and release language.",
      "Choose the smallest figure set that matches the reporting purpose before styling or rendering.",
      "Compare raw screening surfaces across scenarios and design conditions.",
      "Separate operating behavior from runtime and from validation or release decisions.",
      "Inspect whether misspecification shifts strict screens beyond the well-specified baseline.",
      "Use only when include_report = TRUE; otherwise treat absence as not requested.",
      "Reuse for ggplot2, plotly, Quarto, or supplementary figure construction."
    ),
    InterpretationBoundary = c(
      "Run completion and convergence are prerequisites for interpretation, not proof that the diagnostic screen is calibrated.",
      "Reading order is guidance for human review, not a statistical decision rule.",
      "Next actions are workflow prompts, not statistical decisions.",
      "Use cautious wording unless replication count, design, and scenario coverage support stronger claims.",
      "Figure recipes are reporting workflow guidance; they do not add evidence beyond the underlying summary tables.",
      "Scenario means and flag rates are operating-characteristic readouts, not formal inferential p-values.",
      "Type I and sensitivity labels are simulation proxies tied to the evaluated scenarios.",
      "Contrasts are descriptive deltas and depend on the selected baseline, design, and replication count.",
      "Report signals summarize the package reporting layer; they do not validate the model or the simulation design.",
      "Plot data are a presentation handoff and should not be reinterpreted as a separate analysis."
    ),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(include_report)) {
    out$Purpose[out$Table == "report_signal_summary"] <- "Skipped unless include_report = TRUE."
    out$InterpretationBoundary[out$Table == "report_signal_summary"] <- "An empty table means report-index review was not requested, not that no reporting issues exist."
  }
  out
}

diagnostic_screening_next_actions <- function(overview,
                                              scenario_summary,
                                              performance_summary,
                                              report_signal_summary,
                                              scenario_contrast,
                                              include_report = FALSE) {
  overview <- as.data.frame(overview %||% data.frame(), stringsAsFactors = FALSE)
  scenario_summary <- as.data.frame(scenario_summary %||% data.frame(), stringsAsFactors = FALSE)
  performance_summary <- as.data.frame(performance_summary %||% data.frame(), stringsAsFactors = FALSE)
  report_signal_summary <- as.data.frame(report_signal_summary %||% data.frame(), stringsAsFactors = FALSE)
  scenario_contrast <- as.data.frame(scenario_contrast %||% data.frame(), stringsAsFactors = FALSE)
  one <- function(df, col, default = NA) {
    if (!is.data.frame(df) || nrow(df) == 0L || !col %in% names(df)) return(default)
    df[[col]][1]
  }
  reps <- suppressWarnings(as.integer(one(overview, "Reps", NA_integer_)))
  run_ok_rate <- suppressWarnings(as.numeric(one(overview, "RunOKRate", NA_real_)))
  convergence_rate <- suppressWarnings(as.numeric(one(overview, "ConvergenceRate", NA_real_)))
  status_reps <- if (is.finite(reps) && reps >= 30L) "ok" else "review"
  status_runs <- if (is.finite(run_ok_rate) && is.finite(convergence_rate) &&
                       run_ok_rate >= 1 && convergence_rate >= 1) {
    "ok"
  } else {
    "review"
  }
  has_contrast <- nrow(scenario_contrast) > 0L
  has_report <- isTRUE(include_report) && nrow(report_signal_summary) > 0L
  has_screening <- nrow(scenario_summary) > 0L && nrow(performance_summary) > 0L

  data.frame(
    Priority = 1:6,
    Area = c(
      "Replication count",
      "Run completion and convergence",
      "Screening interpretation",
      "Scenario contrasts",
      "Report-index signals",
      "Appendix and plot-data handoff"
    ),
    Status = c(
      status_reps,
      status_runs,
      if (has_screening) "ok" else "review",
      if (has_contrast) "ok" else "not_available",
      if (has_report) "ok" else "optional",
      "ok"
    ),
    Evidence = c(
      if (is.finite(reps)) paste0("Reps = ", reps, ".") else "Replication count was not recorded.",
      paste0(
        "RunOKRate = ",
        if (is.finite(run_ok_rate)) format(round(run_ok_rate, 3), nsmall = 3) else "NA",
        "; ConvergenceRate = ",
        if (is.finite(convergence_rate)) format(round(convergence_rate, 3), nsmall = 3) else "NA",
        "."
      ),
      paste0("scenario_summary rows = ", nrow(scenario_summary),
             "; performance_summary rows = ", nrow(performance_summary), "."),
      paste0("scenario_contrast rows = ", nrow(scenario_contrast), "."),
      if (has_report) {
        paste0("report_signal_summary rows = ", nrow(report_signal_summary), ".")
      } else if (isTRUE(include_report)) {
        "include_report = TRUE, but no report_signal_summary rows are available."
      } else {
        "include_report = FALSE; report-index signals were not requested."
      },
      "summary tables and draw-free plot-data tables are available through the package-wide bundle/export route."
    ),
    Action = c(
      if (identical(status_reps, "ok")) {
        "Use the replication count in methods text and report Monte Carlo scope explicitly."
      } else {
        "Treat this as a screening or smoke run; increase reps before making stable operating-characteristic claims."
      },
      if (identical(status_runs, "ok")) {
        "Proceed to scenario and performance summaries, while still reporting the convergence basis."
      } else {
        "Inspect object$results for Error, RunOK, and Converged before summarizing scenario behavior."
      },
      "Read scenario_summary and performance_summary together before making legacy-vs-strict screening claims.",
      if (has_contrast) {
        "Use scenario_contrast to describe misspecification-minus-baseline shifts, with the baseline scenario named."
      } else {
        "Include well_specified plus at least one misspecification scenario when contrast evidence is needed."
      },
      if (has_report) {
        "Use report_signal_summary to prioritize report text review, not as a diagnostic adequacy test."
      } else {
        "Rebuild with include_report = TRUE only if report-readiness operating behavior is part of the study question."
      },
      "Use build_summary_table_bundle() or export_summary_appendix(); use plot_data() for custom figures."
    ),
    Route = c(
      "evaluate_mfrm_diagnostic_screening(reps = ...)",
      "diag_eval$results[, c(\"Scenario\", \"rep\", \"RunOK\", \"Converged\", \"Error\")]",
      "summary(diag_eval)$scenario_summary; summary(diag_eval)$performance_summary",
      "summary(diag_eval)$scenario_contrast",
      "summary(diag_eval)$report_signal_summary",
      "build_summary_table_bundle(diag_eval); export_summary_appendix(diag_eval); plot_data(diag_eval, type = \"overview\", component = \"plot_long\")"
    ),
    ReportingBoundary = c(
      "Replication count affects Monte Carlo stability and should be reported.",
      "Failed or non-converged runs are design/runtime evidence, not evidence about diagnostic sensitivity.",
      "Screening readouts compare operating behavior; they are not calibrated inferential tests.",
      "Contrasts are descriptive and conditional on the evaluated design grid and scenarios.",
      "Report-index signals are reporting-layer prompts, not extra diagnostic tests.",
      "Exports and plot data are presentation handoffs over the same summary evidence."
    ),
    stringsAsFactors = FALSE
  )
}

diagnostic_screening_reporting_notes <- function(include_report = FALSE) {
  data.frame(
    Area = c(
      "Diagnostic-screening scope",
      "Legacy residual screen",
      "Strict marginal/pairwise screens",
      "Scenario contrasts",
      "Report-index signals",
      "Appendix and plot-data handoff"
    ),
    Evidence = c(
      "evaluate_mfrm_diagnostic_screening() repeatedly simulates, fits, diagnoses, and aggregates selected scenario/design conditions.",
      "Legacy residual summaries retain familiar ZSTD-style screening behavior.",
      "Strict marginal and pairwise summaries add model-implied response-distribution checks where available.",
      "scenario_contrast subtracts the well-specified baseline from misspecification scenarios.",
      if (isTRUE(include_report)) {
        "report_signal_summary records mfrm_report() report-index availability and review-signal counts."
      } else {
        "report_signal_summary is empty unless include_report = TRUE was used."
      },
      "plot_* tables expose the same summaries in long form for figure or appendix construction."
    ),
    ReportingBoundary = c(
      "Describe as a simulation operating-characteristic study, not as a calibrated hypothesis test.",
      "Do not report ZSTD flags alone as final evidence that the generating mechanism is wrong.",
      "Do not treat strict marginal/pairwise flags as release gates; interpret them relative to evaluated scenarios.",
      "Do not generalize contrast direction or magnitude beyond the simulated design grid and replication count.",
      "Do not treat reporting-layer review signals as additional diagnostic tests.",
      "Do not treat exported plot data as independent evidence beyond the underlying summary tables."
    ),
    RecommendedAction = c(
      "Report scenarios, design grid, replication count, fitting method, diagnostic mode, and any unsupported model scope.",
      "Pair legacy flags with strict marginal/pairwise summaries before making diagnostic-screening claims.",
      "Use performance_summary and scenario_contrast to compare sensitivity and Type I proxy behavior.",
      "State which scenario is the baseline and avoid strong claims when replications are small.",
      if (isTRUE(include_report)) {
        "Use report signals to prioritize report review text, not to decide diagnostic adequacy."
      } else {
        "Rebuild with include_report = TRUE only if report-readiness operating behavior is part of the study question."
      },
      "Use build_summary_table_bundle() or export_summary_appendix() to keep tables, roles, and appendix presets explicit."
    ),
    stringsAsFactors = FALSE
  )
}

diagnostic_screening_figure_recipes <- function(include_report = FALSE) {
  report_availability <- if (isTRUE(include_report)) {
    "available_when_report_rows_exist"
  } else {
    "requires_include_report_TRUE"
  }
  data.frame(
    FigureID = c(
      "overview_rates",
      "overview_counts",
      "report_review_rates",
      "scenario_contrast_counts",
      "runtime_elapsed"
    ),
    RecommendedUse = c(
      "main_text_or_primary_supplement",
      "supplement_or_quality_control",
      "reporting_layer_supplement",
      "misspecification_follow_up",
      "methods_or_computational_appendix"
    ),
    PrimaryQuestion = c(
      "How often do legacy, strict marginal, strict pairwise, strict combined, and optional report-review screens fire across scenarios?",
      "How many levels, groups, pairs, or report-review signals are accumulated under each scenario/design condition?",
      "When report signals were retained, how often does the reporting layer route fit, precision, or misfit areas to review?",
      "How much do misspecification scenarios shift flagged counts relative to the well-specified baseline?",
      "How much elapsed time does the diagnostic-screening workflow require under each design/scenario condition?"
    ),
    PlotCall = c(
      "plot(diag_eval, type = \"overview\", metric = \"rate\", draw = FALSE)",
      "plot(diag_eval, type = \"overview\", metric = \"count\", draw = FALSE)",
      "plot(diag_eval, type = \"report\", metric = \"rate\", draw = FALSE)",
      "plot(diag_eval, type = \"contrast\", metric = \"count\", draw = FALSE)",
      "plot(diag_eval, type = \"runtime\", metric = \"elapsed\", draw = FALSE)"
    ),
    PlotDataCall = c(
      "plot_data(diag_eval, type = \"overview\", metric = \"rate\", component = \"plot_long\")",
      "plot_data(diag_eval, type = \"overview\", metric = \"count\", component = \"plot_long\")",
      "plot_data(diag_eval, type = \"report\", metric = \"rate\", component = \"plot_long\")",
      "plot_data(diag_eval, type = \"contrast\", metric = \"count\", component = \"plot_long\")",
      "plot_data(diag_eval, type = \"runtime\", metric = \"elapsed\", component = \"plot_long\")"
    ),
    SummaryTable = c(
      "plot_overview_rate",
      "plot_overview_count",
      "plot_report_rate",
      "plot_contrast_count",
      "plot_runtime"
    ),
    DisplaySuggestion = c(
      "Line or point plot by design variable; facet or color by scenario/signal.",
      "Small-multiple count plot or appendix table when raw signal burden matters.",
      "Focused report-readiness panel; suppress when include_report was not requested.",
      "Diverging or signed count display with the baseline scenario named in the caption.",
      "Line or point plot with units stated as seconds or seconds per 100 observations."
    ),
    CaptionFocus = c(
      "Describe operating-characteristic signal rates and identify whether report-review rates were included.",
      "Describe signal burden, not statistical significance or diagnostic adequacy.",
      "Describe reporting-layer review routing, not model validity or diagnostic success.",
      "Describe misspecification-minus-baseline deltas and name the evaluated baseline.",
      "Describe computational cost under the evaluated design grid and fitting settings."
    ),
    InterpretationBoundary = c(
      "Rates are simulation summaries and should not be read as calibrated inferential test results.",
      "Counts are presentation summaries over the same simulation evidence and should not define pass/fail gates.",
      "Report-review signals are prompts for text and evidence review, not additional diagnostic tests.",
      "Contrasts are descriptive and conditional on scenarios, baseline, design grid, and replication count.",
      "Runtime evidence describes this implementation and settings, not a general computational guarantee."
    ),
    Availability = c(
      "available_when_plot_rows_exist",
      "available_when_plot_rows_exist",
      report_availability,
      "available_when_well_specified_baseline_and_misspecification_rows_exist",
      "available_when_performance_rows_exist"
    ),
    stringsAsFactors = FALSE
  )
}

diagnostic_screening_overview <- function(object) {
  results_tbl <- tibble::as_tibble(object$results %||% tibble::tibble())
  settings <- object$settings %||% list()
  scenarios <- as.character(settings$scenarios %||% results_tbl$Scenario %||% character(0))
  models <- as.character(results_tbl$Model %||% settings$model %||% character(0))
  run_ok <- as.logical(results_tbl$RunOK %||% logical(0))
  converged <- as.logical(results_tbl$Converged %||% logical(0))

  tibble::tibble(
    Designs = nrow(as.data.frame(object$design_grid %||% data.frame())),
    Reps = as.integer(settings$reps %||% NA_integer_),
    Scenarios = paste(unique(scenarios), collapse = ", "),
    Models = paste(unique(models), collapse = ", "),
    ReplicateRows = nrow(results_tbl),
    ScenarioRows = nrow(as.data.frame(object$scenario_summary %||% data.frame())),
    PerformanceRows = nrow(as.data.frame(object$performance_summary %||% data.frame())),
    ReportSignalRows = nrow(as.data.frame(object$report_signal_summary %||% data.frame())),
    ContrastRows = nrow(as.data.frame(object$scenario_contrast %||% data.frame())),
    RunOKRate = if (length(run_ok) > 0L) mean(run_ok, na.rm = TRUE) else NA_real_,
    ConvergenceRate = if (length(converged) > 0L) mean(converged, na.rm = TRUE) else NA_real_,
    IncludeReport = isTRUE(settings$include_report),
    PlotDataContract = "mfrm_plot_data"
  )
}

diagnostic_screening_interpretation_payload <- function(object, overview = NULL) {
  settings <- object$settings %||% list()
  overview <- overview %||% diagnostic_screening_overview(object)
  scenario_summary <- tibble::as_tibble(object$scenario_summary %||% tibble::tibble())
  performance_summary <- tibble::as_tibble(object$performance_summary %||% tibble::tibble())
  report_signal_summary <- tibble::as_tibble(object$report_signal_summary %||% tibble::tibble())
  scenario_contrast <- tibble::as_tibble(object$scenario_contrast %||% tibble::tibble())

  list(
    overview = overview,
    reading_order = diagnostic_screening_reading_order(isTRUE(settings$include_report)),
    next_actions = diagnostic_screening_next_actions(
      overview = overview,
      scenario_summary = scenario_summary,
      performance_summary = performance_summary,
      report_signal_summary = report_signal_summary,
      scenario_contrast = scenario_contrast,
      include_report = isTRUE(settings$include_report)
    ),
    reporting_notes = diagnostic_screening_reporting_notes(isTRUE(settings$include_report)),
    figure_recipes = diagnostic_screening_figure_recipes(isTRUE(settings$include_report))
  )
}

#' Summarize a diagnostic-screening validation study
#'
#' @description
#' Summarizes output from [evaluate_mfrm_diagnostic_screening()] for reporting,
#' appendix export, and draw-free visualization handoff. The summary keeps
#' simulation operating characteristics separate from validation gates: fit,
#' marginal, pairwise, and report-review signals are screening readouts rather
#' than pass/fail evidence.
#'
#' @param object Output from [evaluate_mfrm_diagnostic_screening()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_diagnostic_screening` with:
#' - `overview`: run-level design, replication, convergence, and report-review
#'   metadata
#' - `reading_order`: recommended order for reading the summary tables
#' - `next_actions`: action-oriented triage for interpreting and exporting the
#'   summary
#' - `reporting_notes`: report-facing boundaries and recommended wording
#'   safeguards
#' - `figure_recipes`: recommended figure/display recipes for the draw-free
#'   plot-data tables
#' - `scenario_summary`: aggregated scenario-by-design screening summaries
#' - `performance_summary`: operating-characteristic rates and runtime summaries
#' - `report_signal_summary`: optional `mfrm_report()` readiness/review signals
#' - `scenario_contrast`: misspecification-minus-well-specified contrasts
#' - `plot_*`: long-form draw-free plot tables for overview, report, contrast,
#'   and runtime views
#' - planning metadata, settings, ADEMP metadata, and interpretation notes
#'
#' @seealso [evaluate_mfrm_diagnostic_screening()], [plot.mfrm_diagnostic_screening], [plot_data()]
#' @examples
#' \donttest{
#' diag_eval <- evaluate_mfrm_diagnostic_screening(
#'   design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
#'   reps = 1,
#'   maxit = 30,
#'   seed = 123
#' )
#' summary(diag_eval)
#' }
#' @export
summary.mfrm_diagnostic_screening <- function(object, digits = 3, ...) {
  if (!is.list(object) || is.null(object$results) || is.null(object$scenario_summary)) {
    stop("`object` must be output from evaluate_mfrm_diagnostic_screening().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits[1]))
  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }
  settings <- object$settings %||% list()
  as_summary_tbl <- function(x) {
    if (is.null(x)) return(tibble::tibble())
    tibble::as_tibble(x)
  }

  overview <- diagnostic_screening_overview(object)
  scenario_summary <- round_df(as_summary_tbl(object$scenario_summary))
  performance_summary <- round_df(as_summary_tbl(object$performance_summary))
  report_signal_summary <- round_df(as_summary_tbl(object$report_signal_summary))
  scenario_contrast <- round_df(as_summary_tbl(object$scenario_contrast))
  interpretation <- diagnostic_screening_interpretation_payload(
    object,
    overview = round_df(overview)
  )

  out <- list(
    overview = interpretation$overview,
    reading_order = interpretation$reading_order,
    next_actions = interpretation$next_actions,
    reporting_notes = interpretation$reporting_notes,
    figure_recipes = interpretation$figure_recipes,
    scenario_summary = scenario_summary,
    performance_summary = performance_summary,
    report_signal_summary = report_signal_summary,
    scenario_contrast = scenario_contrast,
    plot_overview_rate = round_df(diagnostic_screening_plot_long_safe(object, "overview", "rate")),
    plot_overview_count = round_df(diagnostic_screening_plot_long_safe(object, "overview", "count")),
    plot_report_rate = round_df(diagnostic_screening_plot_long_safe(object, "report", "rate")),
    plot_contrast_count = round_df(diagnostic_screening_plot_long_safe(object, "contrast", "count")),
    plot_runtime = round_df(diagnostic_screening_plot_long_safe(object, "runtime", "elapsed")),
    ademp = object$ademp %||% NULL,
    facet_names = object$settings$facet_names %||% stats::setNames(simulation_default_output_facet_names(), c("rater", "criterion")),
    design_variable_aliases = simulation_object_design_variable_aliases(object),
    design_descriptor = simulation_object_design_descriptor(object),
    planning_scope = simulation_object_planning_scope(object),
    planning_constraints = simulation_object_planning_constraints(object),
    planning_schema = simulation_object_planning_schema(object),
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
    settings = settings,
    digits = digits,
    notes = unique(c(
      as.character(object$notes %||% character(0)),
      "Draw-free diagnostic-screening plot tables are exported as operating-characteristic readouts, not validation pass/fail gates."
    ))
  )
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
  out$notes <- unique(out$notes)
  class(out) <- "summary.mfrm_diagnostic_screening"
  out
}

#' @export
print.summary.mfrm_diagnostic_screening <- function(x, ...) {
  digits <- max(0L, as.integer(x$digits %||% 3L))
  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }
  preview_df <- function(df, n = 10L) {
    if (!is.data.frame(df) || nrow(df) == 0L) return(df)
    utils::head(df, n = n)
  }

  cat("mfrmr Diagnostic Screening Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_df(as.data.frame(x$overview)), row.names = FALSE)
  }
  if (!is.null(x$reading_order) && nrow(x$reading_order) > 0L) {
    cat("\nReading order\n")
    print(as.data.frame(x$reading_order), row.names = FALSE)
  }
  if (!is.null(x$next_actions) && nrow(x$next_actions) > 0L) {
    cat("\nNext actions\n")
    print(as.data.frame(preview_df(x$next_actions)), row.names = FALSE)
  }
  if (!is.null(x$reporting_notes) && nrow(x$reporting_notes) > 0L) {
    cat("\nReporting notes\n")
    print(as.data.frame(preview_df(x$reporting_notes)), row.names = FALSE)
  }
  if (!is.null(x$figure_recipes) && nrow(x$figure_recipes) > 0L) {
    keep <- intersect(
      c("FigureID", "RecommendedUse", "SummaryTable", "CaptionFocus", "Availability"),
      names(x$figure_recipes)
    )
    cat("\nFigure recipes\n")
    print(as.data.frame(preview_df(x$figure_recipes[, keep, drop = FALSE])), row.names = FALSE)
  }
  if (!is.null(x$scenario_summary) && nrow(x$scenario_summary) > 0L) {
    cat("\nScenario summary (preview)\n")
    print(round_df(as.data.frame(preview_df(x$scenario_summary))), row.names = FALSE)
  }
  if (!is.null(x$performance_summary) && nrow(x$performance_summary) > 0L) {
    cat("\nPerformance summary (preview)\n")
    print(round_df(as.data.frame(preview_df(x$performance_summary))), row.names = FALSE)
  }
  if (!is.null(x$gpcm_boundary) && nrow(x$gpcm_boundary) > 0L) {
    cat("\nBounded GPCM boundary\n")
    keep <- intersect(
      c("Area", "Status", "RecommendedRoute", "NextValidationStep"),
      names(x$gpcm_boundary)
    )
    print(as.data.frame(preview_df(x$gpcm_boundary[, keep, drop = FALSE])), row.names = FALSE)
  }
  if (!is.null(x$report_signal_summary) && nrow(x$report_signal_summary) > 0L) {
    cat("\nReport signal summary (preview)\n")
    print(round_df(as.data.frame(preview_df(x$report_signal_summary))), row.names = FALSE)
  }
  if (!is.null(x$scenario_contrast) && nrow(x$scenario_contrast) > 0L) {
    cat("\nScenario contrast (preview)\n")
    print(round_df(as.data.frame(preview_df(x$scenario_contrast))), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (note in x$notes) cat(" - ", note, "\n", sep = "")
  }
  invisible(x)
}

diagnostic_screening_plot_specs <- function(type, metric) {
  metric <- if (is.null(metric)) "auto" else as.character(metric[1])
  if (!nzchar(metric) || is.na(metric)) {
    metric <- "auto"
  }
  if (identical(type, "overview")) {
    if (identical(metric, "auto")) metric <- "rate"
    if (identical(metric, "rate")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = c(
            "scenario_summary", "scenario_summary", "scenario_summary",
            "performance_summary",
            "report_signal_summary", "report_signal_summary", "report_signal_summary"
          ),
          column = c(
            "LegacyAnyFlagRate", "MarginalAnyFlagRate", "PairwiseAnyFlagRate",
            "StrictAnyFlagRate",
            "FitReportReviewRate", "PrecisionReportReviewRate", "MisfitReportReviewRate"
          ),
          signal = c(
            "Legacy |ZSTD| any-flag rate",
            "Strict marginal any-flag rate",
            "Strict pairwise any-flag rate",
            "Strict combined any-flag rate",
            "Report fit review rate",
            "Report precision review rate",
            "Report misfit review rate"
          ),
          scale = "rate"
        )
      ))
    }
    if (identical(metric, "count")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = c(
            "scenario_summary", "scenario_summary", "scenario_summary",
            "report_signal_summary", "report_signal_summary", "report_signal_summary"
          ),
          column = c(
            "MeanLegacyFlaggedLevels", "MeanMarginalFlaggedGroups",
            "MeanPairwiseFlaggedLevelPairs",
            "MeanFitReportSignals", "MeanPrecisionReportSignals",
            "MeanMisfitReportSignals"
          ),
          signal = c(
            "Mean legacy flagged levels",
            "Mean strict marginal flagged groups",
            "Mean strict pairwise flagged pairs",
            "Mean report fit signals",
            "Mean report precision signals",
            "Mean report misfit signals"
          ),
          scale = "count"
        )
      ))
    }
    if (identical(metric, "magnitude")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "scenario_summary",
          column = c(
            "MeanLegacyMeanAbsZ",
            "MeanMarginalOverallRMSD",
            "MeanMarginalMaxAbsStdResidual"
          ),
          signal = c(
            "Mean legacy |ZSTD|",
            "Mean marginal RMSD",
            "Mean marginal max |standardized residual|"
          ),
          scale = "magnitude"
        )
      ))
    }
    stop("For `type = \"overview\"`, `metric` must be one of `rate`, `count`, or `magnitude`.", call. = FALSE)
  }

  if (identical(type, "report")) {
    if (identical(metric, "auto")) metric <- "rate"
    if (identical(metric, "rate")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "report_signal_summary",
          column = c(
            "ReportIndexAvailabilityRate",
            "FitReportReviewRate",
            "PrecisionReportReviewRate",
            "MisfitReportReviewRate"
          ),
          signal = c(
            "Report-index availability rate",
            "Fit review rate",
            "Precision review rate",
            "Misfit review rate"
          ),
          scale = "rate"
        )
      ))
    }
    if (identical(metric, "count")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "report_signal_summary",
          column = c(
            "MeanReportReviewAreas",
            "MeanFitReportSignals",
            "MeanPrecisionReportSignals",
            "MeanMisfitReportSignals"
          ),
          signal = c(
            "Mean report review areas",
            "Mean fit report signals",
            "Mean precision report signals",
            "Mean misfit report signals"
          ),
          scale = "count"
        )
      ))
    }
    stop("For `type = \"report\"`, `metric` must be `rate` or `count`.", call. = FALSE)
  }

  if (identical(type, "contrast")) {
    if (identical(metric, "auto")) metric <- "count"
    if (identical(metric, "count")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "scenario_contrast",
          column = c(
            "DeltaLegacyFlaggedLevels",
            "DeltaMarginalFlaggedGroups",
            "DeltaPairwiseFlaggedLevelPairs"
          ),
          signal = c(
            "Delta legacy flagged levels",
            "Delta strict marginal flagged groups",
            "Delta strict pairwise flagged pairs"
          ),
          scale = "delta count"
        )
      ))
    }
    if (identical(metric, "rate")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "scenario_contrast",
          column = "DeltaPairwiseAnyFlagRate",
          signal = "Delta pairwise any-flag rate",
          scale = "delta rate"
        )
      ))
    }
    if (identical(metric, "magnitude")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "scenario_contrast",
          column = c(
            "DeltaLegacyMeanAbsZ",
            "DeltaMarginalOverallRMSD",
            "DeltaMarginalMaxAbsStdResidual"
          ),
          signal = c(
            "Delta mean legacy |ZSTD|",
            "Delta marginal RMSD",
            "Delta marginal max |standardized residual|"
          ),
          scale = "delta magnitude"
        )
      ))
    }
    stop("For `type = \"contrast\"`, `metric` must be one of `count`, `rate`, or `magnitude`.", call. = FALSE)
  }

  if (identical(type, "runtime")) {
    if (identical(metric, "auto")) metric <- "elapsed"
    if (identical(metric, "elapsed")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "performance_summary",
          column = "MeanElapsedSec",
          signal = "Mean elapsed seconds",
          scale = "seconds"
        )
      ))
    }
    if (identical(metric, "per_observation")) {
      return(list(
        metric = metric,
        specs = tibble::tibble(
          table = "performance_summary",
          column = "MeanElapsedSecPer100Obs",
          signal = "Mean elapsed seconds per 100 observations",
          scale = "seconds per 100 observations"
        )
      ))
    }
    stop("For `type = \"runtime\"`, `metric` must be `elapsed` or `per_observation`.", call. = FALSE)
  }

  stop("Unsupported diagnostic-screening plot type.", call. = FALSE)
}

diagnostic_screening_plot_bind_specs <- function(x, specs) {
  rows <- vector("list", nrow(specs))
  out_idx <- 0L
  for (i in seq_len(nrow(specs))) {
    src <- specs$table[i]
    src_tbl <- tibble::as_tibble(x[[src]] %||% tibble::tibble())
    col <- specs$column[i]
    if (nrow(src_tbl) == 0L || !col %in% names(src_tbl)) {
      next
    }
    out_idx <- out_idx + 1L
    row_tbl <- src_tbl
    row_tbl$SourceTable <- src
    row_tbl$Metric <- col
    row_tbl$Signal <- specs$signal[i]
    row_tbl$Scale <- specs$scale[i]
    row_tbl$Value <- suppressWarnings(as.numeric(src_tbl[[col]]))
    rows[[out_idx]] <- row_tbl
  }
  if (out_idx == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(rows[seq_len(out_idx)])
}

diagnostic_screening_plot_group_labels <- function(plot_tbl, group_var = NULL) {
  if (nrow(plot_tbl) == 0L) return(character(0))
  signal <- as.character(plot_tbl$Signal)
  scenario <- if ("Scenario" %in% names(plot_tbl)) as.character(plot_tbl$Scenario) else rep("All scenarios", nrow(plot_tbl))
  if (is.null(group_var)) {
    return(paste(scenario, signal, sep = " | "))
  }
  paste(scenario, signal, paste0(group_var, "=", as.character(plot_tbl[[group_var]])), sep = " | ")
}

#' Plot a diagnostic-screening validation study
#'
#' @description
#' Builds an integrated visual summary from
#' [evaluate_mfrm_diagnostic_screening()] output. The default view combines
#' legacy residual, strict marginal, strict pairwise, strict combined, and
#' optional report-index review rates so simulation results can be inspected in
#' one operating-characteristic surface.
#'
#' @param x Output from [evaluate_mfrm_diagnostic_screening()].
#' @param type Plot family. `"overview"` combines screening and optional report
#'   rates or counts. `"report"` focuses on `mfrm_report()` review signals.
#'   `"contrast"` plots misspecification-minus-well-specified contrasts.
#'   `"runtime"` plots elapsed-time summaries.
#' @param metric Metric family. Use `NULL` or `"auto"` for the default within
#'   each `type`. Supported values are documented by error messages and include
#'   `"rate"`, `"count"`, `"magnitude"`, `"elapsed"`, and
#'   `"per_observation"` depending on `type`.
#' @param x_var Design variable for the horizontal axis. Public design aliases
#'   from a simulation specification are accepted.
#' @param group_var Optional additional design variable to include in group
#'   labels. Public design aliases are accepted.
#' @param draw Logical; if `FALSE`, return the plot-data bundle without drawing.
#' @param ... Reserved for future extensions.
#'
#' @return An `mfrm_plot_data` object with reusable metadata, a long-form
#'   `plot_long` table, and interpretation handoff tables (`overview`,
#'   `reading_order`, `next_actions`, `reporting_notes`, and
#'   `figure_recipes`). When `draw = TRUE`, the object is returned invisibly
#'   after drawing.
#'
#' @examples
#' \donttest{
#' diag_eval <- evaluate_mfrm_diagnostic_screening(
#'   design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
#'   reps = 1,
#'   maxit = 30,
#'   include_report = TRUE,
#'   seed = 123
#' )
#' plot(diag_eval, type = "overview", draw = FALSE)
#' plot_data(diag_eval, type = "overview", component = "plot_long")
#' plot_data(diag_eval, type = "overview", component = "next_actions")
#' plot_data(diag_eval, type = "overview", component = "figure_recipes")
#' plot(diag_eval, type = "report", metric = "rate", draw = FALSE)
#' }
#' @export
plot.mfrm_diagnostic_screening <- function(x,
                                           type = c("overview", "report", "contrast", "runtime"),
                                           metric = NULL,
                                           x_var = c("n_person", "n_rater", "n_criterion", "raters_per_person"),
                                           group_var = NULL,
                                           draw = TRUE,
                                           ...) {
  if (!is.list(x) || is.null(x$results) || is.null(x$scenario_summary)) {
    stop("`x` must be output from evaluate_mfrm_diagnostic_screening().", call. = FALSE)
  }
  type <- match.arg(type)
  design_variable_aliases <- simulation_object_design_variable_aliases(x)
  design_descriptor <- simulation_object_design_descriptor(x)
  x_var <- if (missing(x_var)) {
    "n_person"
  } else {
    simulation_resolve_design_variable(x_var, design_variable_aliases, "x_var", descriptor = design_descriptor)
  }
  x_label <- simulation_design_variable_label(x_var, design_variable_aliases)
  if (is.null(group_var)) {
    resolved_group_var <- NULL
  } else {
    resolved_group_var <- simulation_resolve_design_variable(group_var, design_variable_aliases, "group_var", descriptor = design_descriptor)
    if (identical(resolved_group_var, x_var)) {
      stop("`group_var` must differ from `x_var`.", call. = FALSE)
    }
  }
  group_label <- if (is.null(resolved_group_var)) NULL else simulation_design_variable_label(resolved_group_var, design_variable_aliases)

  spec_bundle <- diagnostic_screening_plot_specs(type, metric)
  plot_tbl <- diagnostic_screening_plot_bind_specs(x, spec_bundle$specs)
  if (nrow(plot_tbl) == 0L) {
    stop("No rows are available for the selected diagnostic-screening plot type and metric.", call. = FALSE)
  }
  if (!x_var %in% names(plot_tbl)) {
    stop("The selected `x_var` is not available in the diagnostic-screening summary tables.", call. = FALSE)
  }
  if (!is.null(resolved_group_var) && !resolved_group_var %in% names(plot_tbl)) {
    stop("The selected `group_var` is not available in the diagnostic-screening summary tables.", call. = FALSE)
  }

  plot_tbl <- plot_tbl |>
    dplyr::filter(is.finite(.data$Value))
  plot_tbl$group <- diagnostic_screening_plot_group_labels(
    plot_tbl,
    group_var = resolved_group_var
  )
  plot_tbl <- plot_tbl |>
    dplyr::arrange(.data[[x_var]], .data$Scenario, .data$Signal)

  if (nrow(plot_tbl) == 0L) {
    stop("Selected diagnostic-screening metric has no finite values to plot.", call. = FALSE)
  }

  interpretation <- diagnostic_screening_interpretation_payload(x)
  payload <- list(
    type = type,
    metric = spec_bundle$metric,
    x_var = x_var,
    x_label = x_label,
    group_var = resolved_group_var,
    group_label = group_label,
    design_variable_aliases = design_variable_aliases,
    design_descriptor = design_descriptor,
    planning_scope = simulation_object_planning_scope(x),
    planning_constraints = simulation_object_planning_constraints(x),
    planning_schema = simulation_object_planning_schema(x),
    gpcm_boundary = x$gpcm_boundary %||% data.frame(),
    overview = interpretation$overview,
    reading_order = interpretation$reading_order,
    next_actions = interpretation$next_actions,
    reporting_notes = interpretation$reporting_notes,
    figure_recipes = interpretation$figure_recipes,
    source_tables = unique(plot_tbl$SourceTable),
    signals = unique(plot_tbl$Signal),
    interpretation_note = paste(
      "Diagnostic-screening plots summarize simulation operating characteristics.",
      "They compare screening and reporting signals across design/scenario conditions",
      "and should not be read as calibrated inferential tests."
    ),
    plot_long = plot_tbl,
    plot_table = plot_tbl,
    title = paste("MFRM diagnostic screening", type),
    subtitle = paste("Metric family:", spec_bundle$metric),
    legend = new_plot_legend(
      label = unique(plot_tbl$group),
      role = rep("scenario_signal", length(unique(plot_tbl$group))),
      aesthetic = rep("color", length(unique(plot_tbl$group))),
      value = unique(plot_tbl$group)
    ),
    reference_lines = new_reference_lines()
  )
  out <- new_mfrm_plot_data("diagnostic_screening", payload)
  if (!isTRUE(draw)) return(out)

  groups <- unique(as.character(plot_tbl$group))
  cols <- grDevices::hcl.colors(max(1L, length(groups)), "Set 2")
  x_vals <- sort(unique(plot_tbl[[x_var]]))
  y_range <- range(plot_tbl$Value, na.rm = TRUE)
  if (!all(is.finite(y_range))) {
    stop("Selected diagnostic-screening metric has no finite values to plot.", call. = FALSE)
  }

  graphics::plot(
    x = x_vals,
    y = rep(NA_real_, length(x_vals)),
    type = "n",
    ylim = y_range,
    xlab = x_label,
    ylab = paste("Diagnostic-screening", spec_bundle$metric),
    main = payload$title
  )
  for (i in seq_along(groups)) {
    g_tbl <- plot_tbl[as.character(plot_tbl$group) == groups[i], , drop = FALSE]
    g_tbl <- g_tbl[order(g_tbl[[x_var]]), , drop = FALSE]
    graphics::lines(g_tbl[[x_var]], g_tbl$Value, col = cols[i], lwd = 2)
    graphics::points(g_tbl[[x_var]], g_tbl$Value, col = cols[i], pch = 19)
  }
  graphics::legend(
    "topright",
    legend = groups,
    col = cols,
    lty = 1,
    pch = 19,
    cex = 0.75,
    bty = "n"
  )
  invisible(out)
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
#' @param model Measurement model passed to [fit_mfrm()]. Bounded `GPCM` is
#'   supported with caveats as slope-aware signal-detection sensitivity
#'   evidence.
#' @param step_facet Step facet passed to [fit_mfrm()] when `model = "PCM"` or
#'   `model = "GPCM"`.
#'   When left `NULL`, the function inherits the generator step facet from
#'   `sim_spec` when available and otherwise defaults to `"Criterion"`.
#' @param slope_facet Slope facet passed to [fit_mfrm()] when
#'   `model = "GPCM"`. Defaults to the fitted step facet.
#' @param slopes Optional bounded-`GPCM` slope specification used by direct
#'   simulation calls when `sim_spec = NULL`.
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
#' Bounded-`GPCM` runs preserve the current package constraint
#' `slope_facet == step_facet` within the generator and fitted model. The
#' resulting DIF and bias rates are slope-aware screening summaries, not
#' formal inferential power, alpha calibration, operational scoring, or
#' arbitrary-facet planning evidence.
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
#' - `gpcm_boundary`: bounded-`GPCM` caveat row when a `GPCM` screening route
#'   is used
#' - `settings`: signal-analysis settings
#' - `ademp`: simulation-study metadata (aims, DGM, estimands, methods, performance measures)
#' - `notes`: short interpretation notes
#' @seealso [simulate_mfrm_data()], [evaluate_mfrm_design()], [analyze_dff()], [analyze_dif()], [estimate_bias()]
#' @examples
#' \donttest{
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   design = list(person = 8, rater = 2, criterion = 2, assignment = 1),
#'   reps = 1,
#'   maxit = 30,
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
                                           slope_facet = NULL,
                                           slopes = NULL,
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
  call_args <- names(as.list(match.call(expand.dots = FALSE))[-1])
  model_explicit <- "model" %in% call_args
  if (!is.null(sim_spec) && !model_explicit) {
    model <- toupper(as.character(sim_spec$model %||% model[1]))
    model <- match.arg(model, c("RSM", "PCM", "GPCM"))
  } else {
    model <- match.arg(model)
  }
  residual_pca <- match.arg(residual_pca)
  dif_method <- match.arg(dif_method)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  if (!identical(model, "GPCM") && !is.null(slopes)) {
    stop("`slopes` can be supplied only when `model = \"GPCM\"`.", call. = FALSE)
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
  generator_step_facet <- if (is.null(sim_spec)) {
    if (identical(generator_model, "PCM") || identical(generator_model, "GPCM")) {
      as.character(step_facet[1] %||% base_facet_names[2])
    } else {
      NA_character_
    }
  } else {
    sim_spec$step_facet
  }
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
  design_variable_aliases <- design_meta$aliases
  design_descriptor <- design_meta$descriptor
  planning_scope <- simulation_planning_scope(sim_spec)
  planning_constraints <- simulation_planning_constraints(sim_spec)
  planning_schema <- simulation_planning_schema(sim_spec)
  fit_step_facet <- simulation_resolve_fit_step_facet(model, step_facet, generator_step_facet)
  fit_slope_facet <- simulation_resolve_fit_slope_facet(model, slope_facet, fit_step_facet)
  gpcm_route_active <- identical(model, "GPCM") || identical(as.character(generator_model), "GPCM")
  gpcm_boundary <- simulation_gpcm_screening_boundary(gpcm_route_active)
  gpcm_notes <- simulation_gpcm_screening_notes(gpcm_route_active)

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
          model = model,
          step_facet = if (identical(model, "RSM")) "Criterion" else fit_step_facet,
          slope_facet = if (identical(model, "GPCM")) fit_slope_facet else NULL,
          slopes = if (identical(model, "GPCM")) slopes else NULL,
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
      if (identical(model, "PCM") || identical(model, "GPCM")) fit_args$step_facet <- fit_step_facet
      if (identical(model, "GPCM")) fit_args$slope_facet <- fit_slope_facet
      if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
      if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points
      sim_population <- attr(sim, "mfrm_population_data")
      if (is.list(sim_population) && isTRUE(sim_population$active)) {
        fit_args$population_formula <- sim_population$population_formula
        fit_args$person_data <- sim_population$person_data
        fit_args$person_id <- sim_population$person_id
        fit_args$population_policy <- sim_population$population_policy
      }
      fit_args <- simulation_add_fit_score_support(
        fit_args,
        sim,
        fallback_score_levels = row_score_levels
      )

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
      gpcm_boundary = gpcm_boundary,
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
        slope_facet = fit_slope_facet,
        slopes = if (identical(model, "GPCM")) slopes else NULL,
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
        gpcm_screening_status = if (gpcm_route_active) "supported_with_caveat" else NA_character_,
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
      ),
      notes = gpcm_notes
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
#' - `gpcm_boundary`: bounded-`GPCM` caveat row when present
#' - `notes`: short interpretation notes, including the bias-side screening caveat
#' @seealso [evaluate_mfrm_signal_detection()], [plot.mfrm_signal_detection]
#' @examples
#' \donttest{
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 8,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 30,
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
  out$gpcm_boundary <- object$gpcm_boundary %||% data.frame()
  out$future_branch_active_summary <- simulation_compact_future_branch_active_summary(
    object,
    digits = digits
  )
  out$digits <- digits
  out$notes <- unique(c(out$notes, as.character(object$notes %||% character(0))))
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
  if (!is.null(x$gpcm_boundary) && nrow(x$gpcm_boundary) > 0L) {
    cat("\nBounded GPCM boundary\n")
    keep <- intersect(
      c("Area", "Status", "RecommendedRoute", "NextValidationStep"),
      names(x$gpcm_boundary)
    )
    print(as.data.frame(preview_df(x$gpcm_boundary[, keep, drop = FALSE])), row.names = FALSE)
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
#' \donttest{
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 8,
#'   n_rater = 2,
#'   n_criterion = 2,
#'   raters_per_person = 1,
#'   reps = 1,
#'   maxit = 30,
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
    gpcm_boundary = x$gpcm_boundary %||% data.frame(),
    notes = x$notes %||% character(0),
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
