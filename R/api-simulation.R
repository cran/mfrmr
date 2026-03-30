#' Simulate long-format many-facet Rasch data for design studies
#'
#' @param n_person Number of persons/respondents.
#' @param n_rater Number of rater facet levels.
#' @param n_criterion Number of criterion/item facet levels.
#' @param raters_per_person Number of raters assigned to each person.
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
#' @param model Measurement model recorded in the simulation setup.
#' @param step_facet Step facet used when `model = "PCM"` and threshold values
#'   vary across levels. Currently `"Criterion"` and `"Rater"` are supported.
#' @param thresholds Optional threshold specification. Use either a numeric
#'   vector of common thresholds or a data frame with columns `StepFacet`,
#'   `Step`/`StepIndex`, and `Estimate`.
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
#'   be left at their defaults except for `seed`.
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
#'    (`RSM` or `PCM`) and sample the response
#'
#' Latent-value generation is explicit:
#' - `latent_distribution = "normal"` draws centered normal person/rater/
#'   criterion values using the supplied standard deviations
#' - `latent_distribution = "empirical"` resamples centered support values
#'   recorded in `sim_spec$empirical_support`
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
#'   `Rater`, `Criterion`, and `Score`. If group labels are simulated or
#'   reused from an observed response skeleton, a `Group` column is included.
#'   If a weighted response skeleton is reused, a `Weight` column is also
#'   included.
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
                               model = c("RSM", "PCM"),
                               step_facet = "Criterion",
                               thresholds = NULL,
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
    thresholds <- sim_spec$threshold_table
    assignment <- sim_spec$assignment
    latent_distribution <- as.character(sim_spec$latent_distribution %||% "normal")
    empirical_support <- sim_spec$empirical_support %||% NULL
    assignment_profiles <- sim_spec$assignment_profiles %||% NULL
    design_skeleton <- sim_spec$design_skeleton %||% NULL
  } else {
    n_person <- as.integer(n_person[1])
    n_rater <- as.integer(n_rater[1])
    n_criterion <- as.integer(n_criterion[1])
    raters_per_person <- as.integer(raters_per_person[1])
    score_levels <- as.integer(score_levels[1])
    model <- match.arg(toupper(as.character(model[1])), c("RSM", "PCM"))
    step_facet <- as.character(step_facet[1] %||% "Criterion")
    assignment <- if (is.null(assignment)) {
      if (identical(raters_per_person, n_rater)) "crossed" else "rotating"
    } else {
      match.arg(tolower(as.character(assignment[1])), c("crossed", "rotating"))
    }
    thresholds <- simulation_build_threshold_table(
      thresholds = thresholds,
      score_levels = score_levels,
      step_span = as.numeric(step_span[1]),
      model = model
    )
    latent_distribution <- "normal"
    empirical_support <- NULL
    assignment_profiles <- NULL
    design_skeleton <- NULL
  }

  if (!is.finite(n_person) || n_person < 2L) stop("`n_person` must be >= 2.", call. = FALSE)
  if (!is.finite(n_rater) || n_rater < 2L) stop("`n_rater` must be >= 2.", call. = FALSE)
  if (!is.finite(n_criterion) || n_criterion < 2L) stop("`n_criterion` must be >= 2.", call. = FALSE)
  if (!is.finite(raters_per_person) || raters_per_person < 1L) stop("`raters_per_person` must be >= 1.", call. = FALSE)
  if (raters_per_person > n_rater) stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  if (!is.finite(score_levels) || score_levels < 2L) stop("`score_levels` must be >= 2.", call. = FALSE)
  if (!step_facet %in% c("Criterion", "Rater")) {
    stop("`step_facet` must currently be either \"Criterion\" or \"Rater\" for simulation.", call. = FALSE)
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
    person_ids <- sprintf("P%03d", seq_len(n_person))
    rater_ids <- if (identical(assignment, "skeleton") && !is.null(design_skeleton)) {
      sort(unique(as.character(design_skeleton$Rater)))
    } else if (identical(assignment, "resampled") && !is.null(assignment_profiles)) {
      sort(unique(as.character(assignment_profiles$Rater)))
    } else {
      sprintf("R%02d", seq_len(n_rater))
    }
    criterion_ids <- if (identical(assignment, "skeleton") && !is.null(design_skeleton)) {
      sort(unique(as.character(design_skeleton$Criterion)))
    } else {
      sprintf("C%02d", seq_len(n_criterion))
    }

    theta <- stats::rnorm(n_person, mean = 0, sd = as.numeric(theta_sd[1]))
    names(theta) <- person_ids

    rater_effects <- sort(stats::rnorm(n_rater, mean = 0, sd = as.numeric(rater_sd[1])))
    names(rater_effects) <- rater_ids

    criterion_effects <- sort(stats::rnorm(n_criterion, mean = 0, sd = as.numeric(criterion_sd[1])))
    names(criterion_effects) <- criterion_ids

    if (identical(latent_distribution, "empirical")) {
      theta <- simulation_sample_empirical_support(empirical_support$person, n_person)
      names(theta) <- person_ids
      rater_effects <- simulation_sample_empirical_support(empirical_support$rater, n_rater)
      names(rater_effects) <- rater_ids
      criterion_effects <- simulation_sample_empirical_support(empirical_support$criterion, n_criterion)
      names(criterion_effects) <- criterion_ids
    }

    if (score_levels == 2L) {
      steps <- thresholds$Estimate[thresholds$StepFacet %in% c("Common", unique(thresholds$StepFacet)[1])][1]
    }

    if (assignment == "crossed") {
      dat <- expand.grid(
        Person = person_ids,
        Rater = rater_ids,
        Criterion = criterion_ids,
        stringsAsFactors = FALSE
      )
    } else if (assignment == "skeleton") {
      dat <- simulation_generate_skeleton_assignment(
        person_ids = person_ids,
        design_skeleton = design_skeleton
      )
    } else if (assignment == "resampled") {
      dat <- simulation_generate_resampled_assignment(
        person_ids = person_ids,
        criterion_ids = criterion_ids,
        assignment_profiles = assignment_profiles
      )
    } else {
      rows <- vector("list", length(person_ids))
      for (i in seq_along(person_ids)) {
        assigned <- rater_ids[((i - 1L) + seq_len(raters_per_person) - 1L) %% n_rater + 1L]
        rows[[i]] <- expand.grid(
          Person = person_ids[i],
          Rater = assigned,
          Criterion = criterion_ids,
          stringsAsFactors = FALSE
        )
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

    allowed_effect_cols <- intersect(c("Group", "Person", "Rater", "Criterion"), names(dat))
    dif_effects <- simulation_normalize_effects(
      effects = dif_effects,
      arg_name = "dif_effects",
      allowed_cols = allowed_effect_cols
    )
    if (nrow(dif_effects) > 0) {
      if (!"Group" %in% names(dif_effects)) {
        stop("`dif_effects` must include a `Group` column.", call. = FALSE)
      }
      if (!any(c("Person", "Rater", "Criterion") %in% names(dif_effects))) {
        stop("`dif_effects` must include at least one of `Person`, `Rater`, or `Criterion`.", call. = FALSE)
      }
    }

    interaction_effects <- simulation_normalize_effects(
      effects = interaction_effects,
      arg_name = "interaction_effects",
      allowed_cols = allowed_effect_cols
    )
    if (nrow(interaction_effects) > 0 && !any(c("Person", "Rater", "Criterion") %in% names(interaction_effects))) {
      stop("`interaction_effects` must include at least one of `Person`, `Rater`, or `Criterion`.", call. = FALSE)
    }

    eta <- theta[dat$Person] - rater_effects[dat$Rater] - criterion_effects[dat$Criterion]
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
    if (identical(model, "PCM")) {
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
      prob_mat <- category_prob_pcm(
        eta = eta,
        step_cum_mat = step_cum_mat,
        criterion_idx = criterion_idx,
        criterion_splits = criterion_splits
      )
    } else {
      step_cum <- c(0, cumsum(threshold_lookup[["Common"]]))
      prob_mat <- category_prob_rsm(eta, step_cum)
    }
    dat$Score <- apply(prob_mat, 1, function(p) sample.int(score_levels, size = 1L, prob = p))
    dat$Study <- "SimulatedDesign"
    keep_cols <- c("Study", "Person", "Rater", "Criterion", "Score")
    if ("Group" %in% names(dat)) keep_cols <- c(keep_cols, "Group")
    if ("Weight" %in% names(dat)) keep_cols <- c(keep_cols, "Weight")
    dat <- dat[, keep_cols]

    attr(dat, "mfrm_truth") <- list(
      person = theta,
      facets = list(
        Rater = rater_effects,
        Criterion = criterion_effects
      ),
      steps = if ("Common" %in% names(threshold_lookup)) threshold_lookup$Common else threshold_table,
      step_table = threshold_table,
      groups = group_assign,
      signals = list(
        dif_effects = dif_effects,
        interaction_effects = interaction_effects
      )
    )
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
      assignment = assignment,
      latent_distribution = latent_distribution,
      threshold_table = threshold_table,
      design_skeleton = design_skeleton
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

simulation_generate_resampled_assignment <- function(person_ids,
                                                     criterion_ids,
                                                     assignment_profiles) {
  if (!is.data.frame(assignment_profiles) || nrow(assignment_profiles) == 0) {
    stop("Empirical assignment profiles are required for `assignment = \"resampled\"`.",
         call. = FALSE)
  }

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
  sampled_templates <- sample(template_people, size = length(person_ids), replace = TRUE)

  rows <- vector("list", length(person_ids))
  for (i in seq_along(person_ids)) {
    assigned <- profiles |>
      dplyr::filter(.data$TemplatePerson == sampled_templates[i]) |>
      dplyr::pull(.data$Rater)
    rows[[i]] <- expand.grid(
      Person = person_ids[i],
      Rater = assigned,
      Criterion = criterion_ids,
      stringsAsFactors = FALSE
    )
    if (!is.null(group_map)) {
      matched_group <- group_map$Group[group_map$TemplatePerson == sampled_templates[i]][1]
      if (is.character(matched_group) && nzchar(matched_group)) {
        rows[[i]]$Group <- matched_group
      }
    }
  }
  dplyr::bind_rows(rows)
}

simulation_generate_skeleton_assignment <- function(person_ids, design_skeleton) {
  if (!is.data.frame(design_skeleton) || nrow(design_skeleton) == 0) {
    stop("Observed design skeleton is required for `assignment = \"skeleton\"`.",
         call. = FALSE)
  }

  skeleton <- tibble::as_tibble(design_skeleton) |>
    dplyr::distinct()
  template_people <- unique(skeleton$TemplatePerson)
  sampled_templates <- sample(template_people, size = length(person_ids), replace = TRUE)
  keep_group <- "Group" %in% names(skeleton)
  keep_weight <- "Weight" %in% names(skeleton)

  rows <- vector("list", length(person_ids))
  for (i in seq_along(person_ids)) {
    template_rows <- skeleton |>
      dplyr::filter(.data$TemplatePerson == sampled_templates[i]) |>
      dplyr::mutate(Person = person_ids[i], .before = 1)
    select_cols <- c("Person", "Rater", "Criterion")
    if (keep_group) select_cols <- c(select_cols, "Group")
    if (keep_weight) select_cols <- c(select_cols, "Weight")
    rows[[i]] <- dplyr::select(template_rows, dplyr::all_of(select_cols))
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
      design_variables = names(design_grid)[names(design_grid) %in% c("n_person", "n_rater", "n_criterion", "raters_per_person")]
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

design_eval_summarize_results <- function(results, rep_overview) {
  results_tbl <- tibble::as_tibble(results)
  rep_tbl <- tibble::as_tibble(rep_overview)

  design_summary <- tibble::tibble()
  if (nrow(results_tbl) > 0) {
    design_summary <- results_tbl |>
      dplyr::group_by(
        .data$design_id,
        .data$Facet,
        .data$n_person,
        .data$n_rater,
        .data$n_criterion,
        .data$raters_per_person
      ) |>
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
      dplyr::arrange(.data$Facet, .data$n_person, .data$n_rater, .data$n_criterion, .data$raters_per_person)
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
#' @param reps Number of replications per design condition.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread of step thresholds on the logit scale.
#' @param fit_method Estimation method passed to [fit_mfrm()].
#' @param model Measurement model passed to [fit_mfrm()].
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
#'   specification.
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
#' Recovery metrics are reported only when the generator and fitted model target
#' the same facet-parameter contract. In practice this means the same
#' `model`, and for `PCM`, the same `step_facet`. When these do not align,
#' recovery fields are set to `NA` and the output records the reason.
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
#' - `design_grid`: evaluated design conditions
#' - `results`: facet-level replicate results
#' - `rep_overview`: run-level status and timing
#' - `settings`: simulation settings
#' - `ademp`: simulation-study metadata (aims, DGM, estimands, methods, performance measures)
#' @seealso [simulate_mfrm_data()], [summary.mfrm_design_evaluation], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- evaluate_mfrm_design(
#'   n_person = c(30, 50),
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 15,
#'   seed = 123
#' )
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
                                 reps = 10,
                                 score_levels = 4,
                                 theta_sd = 1,
                                 rater_sd = 0.35,
                                 criterion_sd = 0.25,
                                 noise_sd = 0,
                                 step_span = 1.4,
                                 fit_method = c("JML", "MML"),
                                 model = c("RSM", "PCM"),
                                 step_facet = NULL,
                                 maxit = 25,
                                 quad_points = 7,
                                 residual_pca = c("none", "overall", "facet", "both"),
                                 sim_spec = NULL,
                                 seed = NULL) {
  fit_method <- match.arg(fit_method)
  model <- match.arg(model)
  residual_pca <- match.arg(residual_pca)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) stop("`reps` must be >= 1.", call. = FALSE)

  design_grid <- expand.grid(
    n_person = as.integer(n_person),
    n_rater = as.integer(n_rater),
    n_criterion = as.integer(n_criterion),
    raters_per_person = as.integer(raters_per_person),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design_grid <- design_grid[design_grid$raters_per_person <= design_grid$n_rater, , drop = FALSE]
  if (nrow(design_grid) == 0) {
    stop("No valid design rows remain after enforcing `raters_per_person <= n_rater`.", call. = FALSE)
  }
  design_grid$design_id <- sprintf("D%02d", seq_len(nrow(design_grid)))
  design_grid <- design_grid[, c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person")]

  generator_model <- if (is.null(sim_spec)) model else sim_spec$model
  generator_step_facet <- if (is.null(sim_spec)) if (identical(generator_model, "PCM")) "Criterion" else NA_character_ else sim_spec$step_facet
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
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
    for (rep in seq_len(reps)) {
      seed_idx <- seed_idx + 1L
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
        facets = c("Rater", "Criterion"),
        score = "Score",
        method = fit_method,
        model = model,
        maxit = maxit
      )
      if (identical(model, "PCM")) fit_args$step_facet <- fit_step_facet
      if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
      if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points

      fit <- tryCatch(do.call(fit_mfrm, fit_args), error = function(e) e)
      diag <- if (inherits(fit, "error")) fit else {
        tryCatch(
          diagnose_mfrm(fit, residual_pca = residual_pca),
          error = function(e) e
        )
      }
      elapsed <- proc.time()[["elapsed"]] - t0
      truth <- attr(sim, "mfrm_truth")

      rep_idx <- rep_idx + 1L
      rep_row <- tibble::tibble(
        design_id = design$design_id,
        rep = rep,
        n_person = design$n_person,
        n_rater = design$n_rater,
        n_criterion = design$n_criterion,
        raters_per_person = design$raters_per_person,
        Observations = nrow(sim),
        MinCategoryCount = min(tabulate(sim$Score, nbins = row_score_levels)),
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

      for (facet in reliability_tbl$Facet) {
        rel_row <- reliability_tbl[reliability_tbl$Facet == facet, , drop = FALSE]
        facet_fit <- fit_tbl[fit_tbl$Facet == facet, , drop = FALSE]
        facet_meas <- measure_tbl[measure_tbl$Facet == facet, , drop = FALSE]
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
          MinCategoryCount = min(tabulate(sim$Score, nbins = row_score_levels)),
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
      design_grid = design_grid,
      results = results,
      rep_overview = rep_overview,
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
#' - `design_summary`: aggregated design-by-facet metrics
#' - `ademp`: simulation-study metadata carried forward from the original object
#' - `notes`: short interpretation notes
#' @seealso [evaluate_mfrm_design()], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- evaluate_mfrm_design(
#'   n_person = c(30, 50),
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 15,
#'   seed = 123
#' )
#' summary(sim_eval)
#' }
#' @export
summary.mfrm_design_evaluation <- function(object, digits = 3, ...) {
  if (!is.list(object) || is.null(object$results) || is.null(object$rep_overview)) {
    stop("`object` must be output from evaluate_mfrm_design().")
  }
  digits <- max(0L, as.integer(digits[1]))
  out <- design_eval_summarize_results(object$results, object$rep_overview)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out$overview <- round_df(out$overview)
  out$design_summary <- round_df(out$design_summary)
  out$ademp <- object$ademp %||% NULL
  class(out) <- "summary.mfrm_design_evaluation"
  out
}

#' Plot a design-simulation study
#'
#' @param x Output from [evaluate_mfrm_design()].
#' @param facet Facet to visualize.
#' @param metric Metric to plot.
#' @param x_var Design variable used on the x-axis.
#' @param group_var Optional design variable used for separate lines.
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
#'   `draw = FALSE`, returns that list directly.
#' @seealso [evaluate_mfrm_design()], [summary.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- evaluate_mfrm_design(
#'   n_person = c(30, 50),
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 15,
#'   seed = 123
#' )
#' plot(sim_eval, facet = "Rater", metric = "separation", x_var = "n_person", draw = FALSE)
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
  facet <- match.arg(facet)
  metric <- match.arg(metric)
  x_var <- match.arg(x_var)

  sum_obj <- design_eval_summarize_results(x$results, x$rep_overview)
  plot_tbl <- tibble::as_tibble(sum_obj$design_summary)
  if (nrow(plot_tbl) == 0) stop("No design-summary rows available for plotting.")
  plot_tbl <- plot_tbl[plot_tbl$Facet == facet, , drop = FALSE]
  if (nrow(plot_tbl) == 0) stop("No rows available for facet `", facet, "`.")

  metric_col <- design_eval_match_metric(metric)
  varying <- c("n_person", "n_rater", "n_criterion", "raters_per_person")
  varying <- varying[varying != x_var]
  if (is.null(group_var)) {
    cand <- varying[vapply(plot_tbl[varying], function(col) length(unique(col)) > 1L, logical(1))]
    group_var <- if (length(cand) > 0) cand[1] else NULL
  } else {
    if (!group_var %in% c("n_person", "n_rater", "n_criterion", "raters_per_person")) {
      stop("`group_var` must be one of the evaluated design variables.")
    }
    if (identical(group_var, x_var)) {
      stop("`group_var` must differ from `x_var`.")
    }
  }

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
    group_var = group_var,
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
    xlab = x_var,
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
#'   optimized first when multiple designs pass.
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
#' - `facet_table`: facet-level threshold checks
#' - `design_table`: design-level aggregated checks
#' - `recommended`: the first passing design after ranking
#' - `thresholds`: thresholds used in the recommendation
#' @seealso [evaluate_mfrm_design()], [summary.mfrm_design_evaluation], [plot.mfrm_design_evaluation]
#' @examples
#' \donttest{
#' sim_eval <- evaluate_mfrm_design(
#'   n_person = c(30, 50),
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 15,
#'   seed = 123
#' )
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

  facets <- unique(as.character(facets))
  missing_facets <- setdiff(facets, unique(design_summary$Facet))
  if (length(missing_facets) > 0) {
    stop("Requested facets not found in the design summary: ", paste(missing_facets, collapse = ", "))
  }

  prefer <- intersect(prefer, c("n_person", "raters_per_person", "n_rater", "n_criterion"))
  if (length(prefer) == 0) {
    stop("`prefer` must contain at least one valid design variable.")
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
    dplyr::arrange(.data$Facet, .data$n_person, .data$n_rater, .data$n_criterion, .data$raters_per_person)

  design_table <- facet_table |>
    dplyr::group_by(
      .data$design_id,
      .data$n_person,
      .data$n_rater,
      .data$n_criterion,
      .data$raters_per_person
    ) |>
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
      )
    ),
    class = "mfrm_design_recommendation"
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

signal_eval_summary <- function(results, rep_overview) {
  results_tbl <- tibble::as_tibble(results)
  rep_tbl <- tibble::as_tibble(rep_overview)

  detection_summary <- tibble::tibble()
  if (nrow(results_tbl) > 0) {
    detection_summary <- results_tbl |>
      dplyr::group_by(
        .data$design_id,
        .data$n_person,
        .data$n_rater,
        .data$n_criterion,
        .data$raters_per_person,
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
      dplyr::arrange(.data$n_person, .data$n_rater, .data$n_criterion, .data$raters_per_person)
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
#' @param model Measurement model passed to [fit_mfrm()].
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
#'   tables stored in `sim_spec`.
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
#' - `design_grid`: evaluated design conditions
#' - `results`: replicate-level detection results
#' - `rep_overview`: run-level status and timing
#' - `settings`: signal-analysis settings
#' - `ademp`: simulation-study metadata (aims, DGM, estimands, methods, performance measures)
#' @seealso [simulate_mfrm_data()], [evaluate_mfrm_design()], [analyze_dff()], [analyze_dif()], [estimate_bias()]
#' @examples
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 20,
#'   n_rater = 3,
#'   n_criterion = 3,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 10,
#'   bias_max_iter = 1,
#'   seed = 123
#' ))
#' s_sig <- summary(sig_eval)
#' s_sig$detection_summary[, c("n_person", "DIFPower", "BiasScreenRate")]
#' @export
evaluate_mfrm_signal_detection <- function(n_person = c(30, 50, 100),
                                           n_rater = c(4),
                                           n_criterion = c(4),
                                           raters_per_person = n_rater,
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
                                           model = c("RSM", "PCM"),
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
  residual_pca <- match.arg(residual_pca)
  dif_method <- match.arg(dif_method)
  if (!is.null(sim_spec) && !inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
  }
  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) stop("`reps` must be >= 1.", call. = FALSE)
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

  design_grid <- expand.grid(
    n_person = as.integer(n_person),
    n_rater = as.integer(n_rater),
    n_criterion = as.integer(n_criterion),
    raters_per_person = as.integer(raters_per_person),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design_grid <- design_grid[design_grid$raters_per_person <= design_grid$n_rater, , drop = FALSE]
  if (nrow(design_grid) == 0) {
    stop("No valid design rows remain after enforcing `raters_per_person <= n_rater`.", call. = FALSE)
  }
  design_grid$design_id <- sprintf("S%02d", seq_len(nrow(design_grid)))
  design_grid <- design_grid[, c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person")]

  generator_model <- if (is.null(sim_spec)) model else sim_spec$model
  generator_step_facet <- if (is.null(sim_spec)) if (identical(generator_model, "PCM")) "Criterion" else NA_character_ else sim_spec$step_facet
  generator_assignment <- if (is.null(sim_spec)) "design_dependent" else sim_spec$assignment
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
    dif_target <- signal_eval_resolve_level(dif_level, "C", design$n_criterion, "dif_level")
    bias_rater_target <- signal_eval_resolve_level(bias_rater, "R", design$n_rater, "bias_rater")
    bias_criterion_target <- signal_eval_resolve_level(bias_criterion, "C", design$n_criterion, "bias_criterion")

    for (rep in seq_len(reps)) {
      seed_idx <- seed_idx + 1L
      dif_tbl <- tibble::tibble(
        Group = focal_group,
        Criterion = dif_target,
        Effect = as.numeric(dif_effect[1])
      )
      bias_tbl <- tibble::tibble(
        Rater = bias_rater_target,
        Criterion = bias_criterion_target,
        Effect = as.numeric(bias_effect[1])
      )

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
        facets = c("Rater", "Criterion"),
        score = "Score",
        method = fit_method,
        model = model,
        maxit = maxit
      )
      if (identical(model, "PCM")) fit_args$step_facet <- fit_step_facet
      if ("Weight" %in% names(sim)) fit_args$weight <- "Weight"
      if (identical(fit_method, "MML")) fit_args$quad_points <- quad_points

      fit <- tryCatch(do.call(fit_mfrm, fit_args), error = function(e) e)
      diag <- if (inherits(fit, "error")) fit else {
        tryCatch(diagnose_mfrm(fit, residual_pca = residual_pca), error = function(e) e)
      }
      dif <- if (inherits(diag, "error")) diag else {
        tryCatch(
          analyze_dff(
            fit, diag,
            facet = "Criterion",
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
          estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = bias_max_iter),
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
      design_grid = design_grid,
      results = dplyr::bind_rows(result_rows),
      rep_overview = dplyr::bind_rows(rep_rows),
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
#' - `detection_summary`: aggregated detection rates by design
#' - `ademp`: simulation-study metadata carried forward from the original object
#' - `notes`: short interpretation notes, including the bias-side screening caveat
#' @seealso [evaluate_mfrm_signal_detection()], [plot.mfrm_signal_detection]
#' @examples
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 20,
#'   n_rater = 3,
#'   n_criterion = 3,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 10,
#'   bias_max_iter = 1,
#'   seed = 123
#' ))
#' summary(sig_eval)
#' @export
summary.mfrm_signal_detection <- function(object, digits = 3, ...) {
  if (!is.list(object) || is.null(object$results) || is.null(object$rep_overview)) {
    stop("`object` must be output from evaluate_mfrm_signal_detection().")
  }
  digits <- max(0L, as.integer(digits[1]))
  out <- signal_eval_summary(object$results, object$rep_overview)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out$overview <- round_df(out$overview)
  out$detection_summary <- round_df(out$detection_summary)
  out$ademp <- object$ademp %||% NULL
  class(out) <- "summary.mfrm_signal_detection"
  out
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
#' @param metric Metric to plot. For `signal = "bias"`, `metric = "power"`
#'   maps to the screening hit rate (`BiasScreenRate`).
#' @param x_var Design variable used on the x-axis.
#' @param group_var Optional design variable used for separate lines.
#' @param draw If `TRUE`, draw with base graphics; otherwise return plotting data.
#' @param ... Reserved for generic compatibility.
#'
#' @return If `draw = TRUE`, invisibly returns plotting data. If `draw = FALSE`,
#'   returns that plotting-data list directly. The returned list includes
#'   `display_metric` and `interpretation_note` so callers can label bias-side
#'   plots as screening summaries rather than formal power/error-rate displays.
#' @seealso [evaluate_mfrm_signal_detection()], [summary.mfrm_signal_detection]
#' @examples
#' sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
#'   n_person = 20,
#'   n_rater = 3,
#'   n_criterion = 3,
#'   raters_per_person = 2,
#'   reps = 1,
#'   maxit = 10,
#'   bias_max_iter = 1,
#'   seed = 123
#' ))
#' plot(sig_eval, signal = "dif", metric = "power", x_var = "n_person", draw = FALSE)
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
  x_var <- match.arg(x_var)

  sum_obj <- signal_eval_summary(x$results, x$rep_overview)
  plot_tbl <- tibble::as_tibble(sum_obj$detection_summary)
  if (nrow(plot_tbl) == 0) stop("No detection-summary rows available for plotting.")

  metric_col <- signal_eval_metric_col(signal, metric)
  varying <- c("n_person", "n_rater", "n_criterion", "raters_per_person")
  varying <- varying[varying != x_var]
  if (is.null(group_var)) {
    cand <- varying[vapply(plot_tbl[varying], function(col) length(unique(col)) > 1L, logical(1))]
    group_var <- if (length(cand) > 0) cand[1] else NULL
  } else {
    if (!group_var %in% c("n_person", "n_rater", "n_criterion", "raters_per_person")) {
      stop("`group_var` must be one of the evaluated design variables.")
    }
    if (identical(group_var, x_var)) {
      stop("`group_var` must differ from `x_var`.")
    }
  }

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
    group_var = group_var,
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
    xlab = x_var,
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
