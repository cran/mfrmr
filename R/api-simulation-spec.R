#' Build an explicit simulation specification for MFRM design studies
#'
#' @param n_person Number of persons/respondents to generate.
#' @param n_rater Number of rater facet levels to generate.
#' @param n_criterion Number of criterion/item facet levels to generate.
#' @param raters_per_person Number of raters assigned to each person.
#' @param design Optional named design override supplied as a named list,
#'   named vector, or one-row data frame. Names may use canonical variables
#'   (`n_person`, `n_rater`, `n_criterion`, `raters_per_person`), current
#'   public aliases implied by `facet_names` (for example `n_judge`,
#'   `n_task`, `judge_per_person`), or role keywords (`person`, `rater`,
#'   `criterion`, `assignment`). The schema-only future branch input
#'   `design$facets = c(person = ..., judge = ..., task = ...)` is also
#'   accepted for the currently exposed facet keys. Do not specify the same
#'   variable through both `design` and the scalar count arguments.
#' @param score_levels Number of ordered score categories.
#' @param theta_sd Standard deviation of simulated person measures.
#' @param rater_sd Standard deviation of simulated rater severities.
#' @param criterion_sd Standard deviation of simulated criterion difficulties.
#' @param noise_sd Optional observation-level noise added to the linear predictor.
#' @param step_span Spread used to generate equally spaced thresholds when
#'   `thresholds = NULL`.
#' @param thresholds Optional threshold specification. Use a numeric vector of
#'   common thresholds; a named list such as `list(C01 = c(-1, 0, 1))`; a
#'   numeric matrix with one row per `StepFacet` and one column per step; or a
#'   long data frame with columns `StepFacet`, `Step`/`StepIndex`, and
#'   `Estimate`.
#' @param model Measurement model recorded in the simulation specification.
#' @param step_facet Step facet used when `model = "PCM"` and threshold values
#'   vary across levels.
#' @param slope_facet Slope facet used when `model = "GPCM"`. The current
#'   bounded `GPCM` branch requires `slope_facet == step_facet`.
#' @param slopes Optional slope specification for `model = "GPCM"`. Use either
#'   a numeric vector aligned to the generated slope-facet levels or a data
#'   frame with columns `SlopeFacet` and `Estimate`. Supplied slopes are treated
#'   as relative discriminations and normalized to the package's geometric-mean-
#'   one identification convention on the log scale. When omitted, slopes
#'   default to 1 for every slope-facet level, giving an exact `PCM`
#'   reduction. The `GPCM` model form follows Muraki's generalized partial
#'   credit model; the package's slope-regime labels are validation stress
#'   labels, not published psychometric cut points.
#' @param facet_names Optional public names for the two simulated non-person
#'   facet columns. Supply either an unnamed character vector of length 2
#'   in rater-like / criterion-like order, or a named vector with names
#'   `c("rater", "criterion")`.
#' @param assignment Assignment design. `"crossed"` means every person sees
#'   every rater; `"rotating"` uses a balanced rotating subset;
#'   `"sparse_linked"` uses an incomplete rating design with optional linking
#'   persons; `"resampled"` reuses empirical person-level rater-assignment
#'   profiles; `"skeleton"` reuses an observed person-by-facet design skeleton.
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
#'   `TemplatePerson` and the public rater-like facet column (optionally
#'   `Group`) describing empirical person-level rater-assignment profiles used
#'   when `assignment = "resampled"`. The canonical name `Rater` is also
#'   accepted.
#' @param design_skeleton Optional data frame with columns `TemplatePerson`,
#'   the public rater-like facet column, and the public criterion-like facet
#'   column (optionally `Group`, `Weight`, and `TemplatePersonReuse`)
#'   describing an observed response skeleton used when
#'   `assignment = "skeleton"`. The canonical names `Rater` and `Criterion`
#'   are also accepted. `TemplatePersonReuse = TRUE` asks
#'   [simulate_mfrm_data()] to keep the template-person order instead of
#'   resampling templates, which is used by [build_peer_review_sim_spec()].
#' @param sparse_controls Optional named list used when
#'   `assignment = "sparse_linked"`. Supported entries are `link_fraction`
#'   (default `0.1`), `link_persons` (overrides `link_fraction`),
#'   `link_raters_per_person` (default `n_rater`), `assignment_mode`
#'   (`"balanced"` or `"random"`), and
#'   `min_common_persons_per_rater_pair` (diagnostic target, default `1`).
#' @param group_levels Optional character vector of group labels.
#' @param dif_effects Optional data frame of true group-linked DIF effects.
#' @param interaction_effects Optional data frame of true interaction effects.
#' @param population_formula Optional one-sided formula describing a
#'   person-level latent-regression population model used when generating
#'   person measures, for example `~ X + G`. When supplied, person measures are
#'   generated from `X %*% beta + e` rather than from `N(0, theta_sd^2)`.
#' @param population_coefficients Optional numeric vector of latent-regression
#'   coefficients corresponding to the design matrix implied by
#'   `population_formula`.
#' @param population_sigma2 Optional residual variance for the latent-regression
#'   person distribution.
#' @param population_covariates Optional template data frame containing one row
#'   per template person and the background variables referenced by
#'   `population_formula`. Numeric/logical and categorical factor/character
#'   variables are expanded through the same `stats::model.matrix()` contract
#'   used by latent-regression fitting. During simulation, template rows are
#'   resampled to the requested `n_person`.
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
#' - optional discrimination structure for bounded `GPCM` (`slope_table`) and
#'   its identified log-slope spread label (`slope_regime`)
#' - assignment design (`assignment`)
#' - optional sparse linked-design controls (`sparse_controls`) when
#'   `assignment = "sparse_linked"`
#' - optional empirical assignment profiles (`assignment_profiles`) with
#'   optional person-level `Group` labels
#' - optional observed response skeleton (`design_skeleton`)
#'   with optional person-level `Group` labels and observation-level `Weight`
#'   values
#' - optional person-level latent-regression population metadata including
#'   `population_formula`, `population_coefficients`, `population_sigma2`, and
#'   a reusable template of person-level covariates, including model-matrix
#'   xlevel/contrast provenance for categorical covariates
#' - `planning_scope`, an explicit record that the current planning/forecasting
#'   helpers target the role-based person x rater-like x criterion-like
#'   design contract rather than a fully arbitrary-facet planner
#' - `planning_constraints`, an explicit record of which design variables can
#'   currently be changed from that specification without rebuilding it
#' - `planning_schema`, a combined schema contract bundling the role descriptor,
#'   scope boundary, current mutability map, a `facet_manifest`, a
#'   schema-only `future_facet_table`, and a matching
#'   `future_design_template`, plus a nested `future_branch_schema` scaffold
#'   for a future arbitrary-facet planning branch
#' - the current `design$facets(...)` parser now normalizes nested facet-count
#'   input through that bundled `future_branch_schema`, whose nested
#'   `design_schema` is now the authoritative schema-only branch object
#' - optional signal tables for DIF and interaction bias
#'
#' The current generator targets the package's standard person x rater x
#' criterion workflow, but the public output names for those two facet roles
#' can now be customized with `facet_names`. This naming layer improves public
#' ergonomics; it does not yet turn the generator into a fully arbitrary-facet
#' simulator. Internally, helper objects keep canonical role mappings so that
#' planning functions can treat the first non-person facet as rater-like and
#' the second as criterion-like. When threshold values are provided by
#' `StepFacet`, the supported step facets are the generated levels of the
#' chosen public rater-like or criterion-like column. For convenience,
#' step-facet-specific thresholds can be supplied as a named list or as a
#' numeric matrix whose row names are `StepFacet` labels.
#' When `model = "GPCM"`, the same public facet naming rules apply to the
#' slope table; the current bounded branch keeps `slope_facet` equal to
#' `step_facet`.
#' The `slope_regime` field summarizes the centered log-slope spread so
#' recovery simulations can be read against the intended generator stress
#' level. Its labels (`unit_slopes`, `near_flat`, `moderate`, and
#' `high_dispersion`) are package validation labels; they are not model-fit
#' decisions and should not be interpreted as literature-derived adequacy
#' thresholds. The GPCM data-generating form follows Muraki (1992,
#' doi:10.1177/014662169201600206), while information-function interpretation
#' follows Muraki (1993, doi:10.1177/014662169301700403). The explicit
#' simulation-specification metadata is intended to support ADEMP-style
#' simulation reporting as described by Morris, White, and Crowther (2019,
#' doi:10.1002/sim.8086).
#'
#' The `assignment = "sparse_linked"` branch follows sparse
#' rater-mediated assessment work in treating incomplete rater assignment as
#' planned missingness rather than incidental nonresponse. Its `link_persons`
#' and `link_raters_per_person` controls emulate a common linking set so users
#' can inspect rater-pair common-person counts before using the generated data
#' in recovery or design studies. This is a design generator and diagnostic
#' metadata layer; it does not impose a universal minimum-linking cutoff.
#' Sparse-design motivation follows Wind, Jones, and Grajeda (2023,
#' doi:10.1177/01466216231182148), Wind and Jones (2018,
#' doi:10.1177/0013164417703733), and DeMars, Shapovalov, and Hathcoat
#' (2023).
#'
#' If `population_formula` is supplied, the simulation specification carries a
#' first-version person-level latent-regression generator. This affects only the
#' person distribution. The current implementation keeps the non-person facets
#' in the existing many-facet Rasch generator and resamples rows from
#' `population_covariates` to the requested design size before computing
#' \eqn{\theta_n = x_n^\top \beta + \varepsilon_n} with
#' \eqn{\varepsilon_n \sim N(0, \sigma^2)}.
#'
#' @section Interpreting output:
#' This object does not contain simulated data. It is a data-generating
#' specification that tells [simulate_mfrm_data()] how to generate them.
#'
#' @return An object of class `mfrm_sim_spec`.
#' @seealso [extract_mfrm_sim_spec()], [simulate_mfrm_data()]
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   design = list(person = 8, rater = 2, criterion = 2, assignment = 1),
#'   assignment = "rotating"
#' )
#' spec$model
#' spec$assignment
#' nrow(spec$threshold_table)
#' }
#' @export
build_mfrm_sim_spec <- function(n_person = 50,
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
                                thresholds = NULL,
                                model = c("RSM", "PCM", "GPCM"),
                                step_facet = NULL,
                                slope_facet = NULL,
                                slopes = NULL,
                                facet_names = NULL,
                                assignment = c("crossed", "rotating", "sparse_linked", "resampled", "skeleton"),
                                latent_distribution = c("normal", "empirical"),
                                empirical_person = NULL,
                                empirical_rater = NULL,
                                empirical_criterion = NULL,
                                assignment_profiles = NULL,
                                design_skeleton = NULL,
                                sparse_controls = NULL,
                                group_levels = NULL,
                                dif_effects = NULL,
                                interaction_effects = NULL,
                                population_formula = NULL,
                                population_coefficients = NULL,
                                population_sigma2 = NULL,
                                population_covariates = NULL) {
  model <- match.arg(toupper(as.character(model[1])), c("RSM", "PCM", "GPCM"))
  assignment <- match.arg(tolower(as.character(assignment[1])), c("crossed", "rotating", "sparse_linked", "resampled", "skeleton"))
  latent_distribution <- match.arg(tolower(as.character(latent_distribution[1])), c("normal", "empirical"))

  facet_names <- simulation_validate_output_facet_names(facet_names)
  supplied_counts <- intersect(
    names(as.list(match.call(expand.dots = FALSE))[-1]),
    c("n_person", "n_rater", "n_criterion", "raters_per_person")
  )
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
  n_person <- simulation_validate_count(design_counts$n_person, "n_person", min_value = 2L)
  n_rater <- simulation_validate_count(design_counts$n_rater, "n_rater", min_value = 2L)
  n_criterion <- simulation_validate_count(design_counts$n_criterion, "n_criterion", min_value = 2L)
  raters_per_person <- simulation_validate_count(design_counts$raters_per_person, "raters_per_person", min_value = 1L)
  score_levels <- simulation_validate_count(score_levels, "score_levels", min_value = 2L)

  if (raters_per_person > n_rater) {
    stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  }
  if (assignment == "crossed" && raters_per_person != n_rater) {
    stop("`assignment = \"crossed\"` requires `raters_per_person == n_rater`.", call. = FALSE)
  }

  if (identical(model, "GPCM")) {
    resolved_facets <- resolve_step_and_slope_facets(
      model = model,
      step_facet = step_facet[1] %||% facet_names[["criterion"]],
      slope_facet = slope_facet,
      facet_names = unname(facet_names)
    )
    step_facet <- resolved_facets$step_facet
    slope_facet <- resolved_facets$slope_facet
  } else {
    step_facet <- simulation_validate_step_facet_name(step_facet[1] %||% facet_names[["criterion"]])
    slope_facet <- NULL
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
  slope_table <- simulation_build_slope_table(
    slopes = slopes,
    model = model,
    slope_facet = slope_facet,
    facet_names = facet_names,
    n_rater = n_rater,
    n_criterion = n_criterion
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
    n_rater = n_rater,
    facet_names = facet_names
  )
  design_skeleton <- simulation_normalize_design_skeleton(
    design_skeleton = design_skeleton,
    assignment = assignment,
    n_rater = n_rater,
    n_criterion = n_criterion,
    facet_names = facet_names
  )
  sparse_controls <- simulation_normalize_sparse_controls(
    sparse_controls = sparse_controls,
    assignment = assignment,
    n_person = n_person,
    n_rater = n_rater,
    raters_per_person = raters_per_person
  )

  dif_effects <- simulation_normalize_effects(
    effects = dif_effects,
    arg_name = "dif_effects",
    allowed_cols = c("Group", "Person", unname(facet_names))
  )
  interaction_effects <- simulation_normalize_effects(
    effects = interaction_effects,
    arg_name = "interaction_effects",
    allowed_cols = c("Group", "Person", unname(facet_names))
  )
  population <- simulation_build_population_spec(
    population_formula = population_formula,
    population_coefficients = population_coefficients,
    population_sigma2 = population_sigma2,
    population_covariates = population_covariates
  )

  spec <- list(
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
      slope_facet = slope_facet,
      assignment = assignment,
      latent_distribution = latent_distribution,
      empirical_support = empirical_support,
      assignment_profiles = assignment_profiles,
      design_skeleton = design_skeleton,
      sparse_controls = sparse_controls,
      threshold_table = threshold_table,
      slope_table = slope_table,
      slope_regime = simulation_gpcm_slope_regime(slope_table),
      facet_names = facet_names,
      facet_levels = list(
        rater = simulation_generated_role_levels(facet_names[["rater"]], "R", n_rater),
        criterion = simulation_generated_role_levels(facet_names[["criterion"]], "C", n_criterion)
      ),
      group_levels = group_levels,
      dif_effects = dif_effects,
      interaction_effects = interaction_effects,
      population = population,
      planning_scope = NULL,
      planning_constraints = NULL,
      planning_schema = NULL,
      source = "manual"
    )
  spec$planning_scope <- simulation_planning_scope(structure(spec, class = "mfrm_sim_spec"))
  spec$planning_constraints <- simulation_planning_constraints(structure(spec, class = "mfrm_sim_spec"))
  spec$planning_schema <- simulation_planning_schema(structure(spec, class = "mfrm_sim_spec"))
  structure(
    spec,
    class = "mfrm_sim_spec"
  )
}

peer_review_validate_scalar_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg_name, "` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }
  x
}

peer_review_resolve_anchor_submissions <- function(n_submission,
                                                   anchor_fraction,
                                                   anchor_submissions) {
  if (is.null(anchor_submissions)) {
    anchor_fraction <- suppressWarnings(as.numeric(anchor_fraction[1]))
    if (!is.finite(anchor_fraction) || anchor_fraction < 0 || anchor_fraction > 1) {
      stop("`anchor_fraction` must be a finite number in [0, 1].", call. = FALSE)
    }
    return(as.integer(ceiling(n_submission * anchor_fraction)))
  }
  value <- suppressWarnings(as.integer(anchor_submissions[1]))
  if (!is.finite(value) || value < 0L || value > n_submission) {
    stop("`anchor_submissions` must be an integer between 0 and `n_submission`.", call. = FALSE)
  }
  value
}

peer_review_rotate_candidates <- function(candidates, offset) {
  if (length(candidates) <= 1L) return(candidates)
  idx <- ((seq_along(candidates) + offset - 2L) %% length(candidates)) + 1L
  candidates[idx]
}

peer_review_choose_reviewers <- function(candidates,
                                         k,
                                         reviewer_load,
                                         offset,
                                         assignment_mode) {
  if (identical(assignment_mode, "random")) {
    return(sample(candidates, size = k, replace = FALSE))
  }
  rotated <- peer_review_rotate_candidates(candidates, offset)
  load <- reviewer_load[rotated]
  rotated[order(load, seq_along(rotated), na.last = TRUE)][seq_len(k)]
}

peer_review_build_design_skeleton <- function(n_submission,
                                              n_criterion,
                                              reviewers_per_submission,
                                              anchor_submissions,
                                              anchor_reviewers_per_submission,
                                              avoid_self_review,
                                              assignment_mode,
                                              seed = NULL) {
  submission_ids <- sprintf("P%03d", seq_len(n_submission))
  reviewer_ids <- submission_ids
  criterion_ids <- simulation_generated_levels("C", n_criterion)
  reviewer_load <- stats::setNames(rep(0L, length(reviewer_ids)), reviewer_ids)
  anchor_ids <- if (anchor_submissions > 0L) {
    submission_ids[seq_len(anchor_submissions)]
  } else {
    character(0)
  }

  assignment_rows <- with_preserved_rng_seed(seed, {
    rows <- vector("list", length(submission_ids))
    for (i in seq_along(submission_ids)) {
      submission <- submission_ids[i]
      candidates <- reviewer_ids
      if (isTRUE(avoid_self_review)) {
        candidates <- setdiff(candidates, submission)
      }
      k <- if (submission %in% anchor_ids) {
        anchor_reviewers_per_submission
      } else {
        reviewers_per_submission
      }
      chosen <- peer_review_choose_reviewers(
        candidates = candidates,
        k = k,
        reviewer_load = reviewer_load,
        offset = i,
        assignment_mode = assignment_mode
      )
      reviewer_load[chosen] <- reviewer_load[chosen] + 1L
      rows[[i]] <- expand.grid(
        TemplatePerson = submission,
        Reviewer = chosen,
        Criterion = criterion_ids,
        stringsAsFactors = FALSE
      )
    }
    dplyr::bind_rows(rows)
  })
  assignment_rows$TemplatePersonReuse <- TRUE
  assignment_rows
}

#' Build a peer-review simulation specification
#'
#' @param n_submission Number of submissions/authors to generate.
#' @param n_criterion Number of rubric criteria.
#' @param reviewers_per_submission Number of peer reviewers assigned to each
#'   ordinary submission.
#' @param anchor_fraction Fraction of submissions treated as common-link
#'   anchor submissions when `anchor_submissions` is not supplied.
#' @param anchor_submissions Optional number of common-link anchor
#'   submissions. Anchor submissions receive
#'   `anchor_reviewers_per_submission` reviewers.
#' @param anchor_reviewers_per_submission Number of reviewers assigned to each
#'   anchor submission. Defaults to all eligible peers when anchors are used
#'   and self-review is disallowed; recorded as 0 when no anchor submissions
#'   are requested.
#' @param avoid_self_review Logical; if `TRUE`, a reviewer is never assigned to
#'   review their own submission.
#' @param assignment_mode Assignment algorithm. `"balanced"` assigns reviewers
#'   with the lowest current load using deterministic rotating tie-breaks.
#'   `"random"` samples eligible reviewers without replacement.
#' @param seed Optional seed used only for random peer-review assignment when
#'   `assignment_mode = "random"`.
#' @param score_levels,theta_sd,reviewer_sd,criterion_sd,noise_sd,step_span
#'   Generator settings passed to [build_mfrm_sim_spec()]. `reviewer_sd` maps
#'   to the standard MFRM rater-severity spread.
#' @param model,step_facet,thresholds Measurement-model settings passed to
#'   [build_mfrm_sim_spec()]. The first public facet is `Reviewer`; the second
#'   is `Criterion`.
#' @param group_levels,dif_effects,interaction_effects Optional signal settings
#'   passed to [build_mfrm_sim_spec()].
#'
#' @details
#' `build_peer_review_sim_spec()` creates a fixed person-by-reviewer-by-rubric
#' skeleton for peer-assessment or peer-review studies. Submissions and peer
#' reviewers share the same ID universe (`P001`, `P002`, ...), so
#' self-review can be structurally excluded and checked in the generated data.
#' The specification uses the existing `assignment = "skeleton"` generator and
#' records peer-review metadata; it does not introduce a new measurement
#' model. MFRM still estimates person/submission measures, reviewer severity,
#' and criterion difficulty, while design-network review can inspect whether
#' the peer-review graph is sufficiently linked.
#'
#' The common-link anchor controls follow the same logic used in sparse
#' rater-mediated designs: when most submissions receive only a few peer
#' reviews, assigning all or many reviewers to a small anchor set can strengthen
#' links among reviewers. The helper labels these rows as design diagnostics,
#' not universal adequacy thresholds for fit, separation, or recovery.
#'
#' @section References:
#' - Farrokhi, F., Esfandiari, R., & Schaefer, E. (2012). A many-facet Rasch
#'   measurement of differential rater severity/leniency in three types of
#'   assessment. *JALT Journal*, 34(1), 79-102.
#'   doi:10.37546/JALTJJ34.1-3.
#' - Uto, M., & Ueno, M. (2020). A generalized many-facet Rasch model and its
#'   Bayesian estimation using Hamiltonian Monte Carlo. *Behaviormetrika*,
#'   47, 469-496. doi:10.1007/s41237-020-00115-7.
#' - DeMars, C. E., Shapovalov, Y. A., & Hathcoat, J. D. (2023).
#'   *Many-Facet Rasch Designs: How Should Raters be Assigned to Examinees?*
#'   NCME presentation.
#'
#' @return An object of class `mfrm_sim_spec` with `peer_review` metadata and
#'   a fixed peer-review design skeleton.
#' @seealso [simulate_mfrm_data()], [build_mfrm_network_review()],
#'   [build_mfrm_sim_spec()]
#' @examples
#' peer_spec <- build_peer_review_sim_spec(
#'   n_submission = 12,
#'   n_criterion = 3,
#'   reviewers_per_submission = 2,
#'   anchor_submissions = 2
#' )
#' peer_spec$peer_review$overview
#' @export
build_peer_review_sim_spec <- function(n_submission = 50,
                                       n_criterion = 4,
                                       reviewers_per_submission = 3,
                                       anchor_fraction = 0.1,
                                       anchor_submissions = NULL,
                                       anchor_reviewers_per_submission = NULL,
                                       avoid_self_review = TRUE,
                                       assignment_mode = c("balanced", "random"),
                                       seed = NULL,
                                       score_levels = 4,
                                       theta_sd = 1,
                                       reviewer_sd = 0.45,
                                       criterion_sd = 0.25,
                                       noise_sd = 0,
                                       step_span = 1.4,
                                       model = c("RSM", "PCM", "GPCM"),
                                       step_facet = "Criterion",
                                       thresholds = NULL,
                                       group_levels = NULL,
                                       dif_effects = NULL,
                                       interaction_effects = NULL) {
  n_submission <- simulation_validate_count(n_submission, "n_submission", min_value = 3L)
  n_criterion <- simulation_validate_count(n_criterion, "n_criterion", min_value = 2L)
  reviewers_per_submission <- simulation_validate_count(
    reviewers_per_submission,
    "reviewers_per_submission",
    min_value = 1L
  )
  avoid_self_review <- peer_review_validate_scalar_logical(
    avoid_self_review,
    "avoid_self_review"
  )
  assignment_mode <- match.arg(tolower(as.character(assignment_mode[1])), c("balanced", "random"))
  max_eligible <- if (isTRUE(avoid_self_review)) n_submission - 1L else n_submission
  if (reviewers_per_submission > max_eligible) {
    stop(
      "`reviewers_per_submission` cannot exceed the number of eligible reviewers per submission.",
      call. = FALSE
    )
  }
  anchor_submissions <- peer_review_resolve_anchor_submissions(
    n_submission = n_submission,
    anchor_fraction = anchor_fraction,
    anchor_submissions = anchor_submissions
  )
  if (is.null(anchor_reviewers_per_submission)) {
    anchor_reviewers_per_submission <- if (anchor_submissions > 0L) max_eligible else 0L
  }
  anchor_reviewers_per_submission <- simulation_validate_count(
    anchor_reviewers_per_submission,
    "anchor_reviewers_per_submission",
    min_value = if (anchor_submissions > 0L) 1L else 0L
  )
  if (anchor_reviewers_per_submission > max_eligible) {
    stop(
      "`anchor_reviewers_per_submission` cannot exceed the number of eligible reviewers per submission.",
      call. = FALSE
    )
  }
  if (anchor_submissions > 0L && anchor_reviewers_per_submission < reviewers_per_submission) {
    stop(
      "`anchor_reviewers_per_submission` must be >= `reviewers_per_submission` when anchor submissions are used.",
      call. = FALSE
    )
  }

  skeleton <- peer_review_build_design_skeleton(
    n_submission = n_submission,
    n_criterion = n_criterion,
    reviewers_per_submission = reviewers_per_submission,
    anchor_submissions = anchor_submissions,
    anchor_reviewers_per_submission = anchor_reviewers_per_submission,
    avoid_self_review = avoid_self_review,
    assignment_mode = assignment_mode,
    seed = seed
  )
  reviewer_levels <- sprintf("P%03d", seq_len(n_submission))
  criterion_levels <- simulation_generated_levels("C", n_criterion)

  spec <- build_mfrm_sim_spec(
    n_person = n_submission,
    n_rater = n_submission,
    n_criterion = n_criterion,
    raters_per_person = reviewers_per_submission,
    score_levels = score_levels,
    theta_sd = theta_sd,
    rater_sd = reviewer_sd,
    criterion_sd = criterion_sd,
    noise_sd = noise_sd,
    step_span = step_span,
    thresholds = thresholds,
    model = model,
    step_facet = step_facet,
    facet_names = c(reviewer = "Reviewer", criterion = "Criterion"),
    assignment = "skeleton",
    design_skeleton = skeleton,
    group_levels = group_levels,
    dif_effects = dif_effects,
    interaction_effects = interaction_effects
  )
  spec$facet_levels <- list(
    rater = reviewer_levels,
    criterion = criterion_levels
  )
  ordinary_rows <- skeleton[!skeleton$TemplatePerson %in%
                              sprintf("P%03d", seq_len(anchor_submissions)), , drop = FALSE]
  anchor_rows <- skeleton[skeleton$TemplatePerson %in%
                            sprintf("P%03d", seq_len(anchor_submissions)), , drop = FALSE]
  spec$peer_review <- list(
    active = TRUE,
    scenario = "peer_review",
    avoid_self_review = avoid_self_review,
    assignment_mode = assignment_mode,
    reviewers_per_submission = reviewers_per_submission,
    anchor_submissions = anchor_submissions,
    anchor_reviewers_per_submission = anchor_reviewers_per_submission,
    reviewer_ids = reviewer_levels,
    criterion_ids = criterion_levels,
    overview = data.frame(
      Scenario = "peer_review",
      Submissions = n_submission,
      Reviewers = n_submission,
      Criteria = n_criterion,
      OrdinaryReviewersPerSubmission = reviewers_per_submission,
      AnchorSubmissions = anchor_submissions,
      AnchorReviewersPerSubmission = anchor_reviewers_per_submission,
      AvoidSelfReview = avoid_self_review,
      AssignmentMode = assignment_mode,
      SkeletonRows = nrow(skeleton),
      OrdinaryRows = nrow(ordinary_rows),
      AnchorRows = nrow(anchor_rows),
      ReviewUse = "design_diagnostic_not_measurement_gate",
      stringsAsFactors = FALSE
    ),
    notes = c(
      "Peer-review simulation uses a fixed skeleton so submissions and reviewers can share the same ID universe.",
      "Common-link anchor submissions are a design device for reviewer linkage, not a fit or recovery adequacy threshold."
    )
  )
  spec$source <- "peer_review"
  spec$source_summary <- c(
    spec$source_summary %||% list(),
    list(peer_review = spec$peer_review$overview)
  )
  spec$planning_scope <- simulation_planning_scope(structure(spec, class = "mfrm_sim_spec"))
  spec$planning_constraints <- simulation_planning_constraints(structure(spec, class = "mfrm_sim_spec"))
  spec$planning_schema <- simulation_planning_schema(structure(spec, class = "mfrm_sim_spec"))
  structure(spec, class = "mfrm_sim_spec")
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
#' - when the fit used the latent-regression branch, the fitted
#'   `population_formula`, coefficient vector, residual variance, and the
#'   stored person-level covariate table, including model-matrix xlevel and
#'   contrast provenance for categorical covariates
#'
#' This is intended as a **fit-derived parametric starting point**, not as a
#' claim that the fitted object perfectly recovers the true data-generating
#' mechanism. Users should review and, if necessary, edit the returned
#' specification before using it for design planning.
#'
#' First-release `GPCM` fits are now supported here for direct data generation
#' and parameter-recovery checks, provided that the returned simulation
#' specification stores both a threshold table and a parallel slope table.
#' The same fit-derived specification can feed caveated role-based design
#' evaluation, population forecasting, and fit-based report/export bundles.
#' Diagnostic/signal-detection design screening, full FACETS score-side
#' contract review, posterior predictive checks, and heavy backend extensions
#' remain outside the bounded-`GPCM` boundary until those downstream contracts
#' are widened explicitly.
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
#' \donttest{
#' toy <- simulate_mfrm_data(
#'   n_person = 8,
#'   n_rater = 3,
#'   n_criterion = 2,
#'   seed = 123
#' )
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' spec <- extract_mfrm_sim_spec(fit, latent_distribution = "empirical")
#' spec$assignment
#' spec$model
#' head(spec$threshold_table)
#' }
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
  fit_model <- as.character(fit$summary$Model[1] %||% fit$config$model %||% "RSM")

  facet_names <- as.character(fit$config$facet_names %||% setdiff(names(prep$levels), "Person"))
  facet_names <- facet_names[!is.na(facet_names) & nzchar(facet_names)]
  if (length(facet_names) != 2L) {
    stop(
      "`extract_mfrm_sim_spec()` currently supports fitted models with exactly two non-person facets. ",
      "This fit has ", length(facet_names), ".",
      call. = FALSE
    )
  }
  assignment_facet <- facet_names[1]
  criterion_facet <- facet_names[2]

  person_tbl <- fit$facets$person %||% tibble::tibble()
  other_tbl <- fit$facets$others %||% tibble::tibble()
  step_tbl <- fit$steps %||% tibble::tibble()

  rater_levels <- as.character(prep$levels[[assignment_facet]])
  criterion_levels <- as.character(prep$levels[[criterion_facet]])
  person_levels <- as.character(prep$levels$Person)

  assignment_counts <- prep$data |>
    dplyr::distinct(.data$Person, .data[[assignment_facet]]) |>
    dplyr::count(.data$Person, name = "RatersPerPerson")
  raters_per_person <- as.integer(round(stats::median(assignment_counts$RatersPerPerson)))
  assignment_profiles <- prep$data |>
    dplyr::distinct(.data$Person, .data[[assignment_facet]]) |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$Person),
      Rater = as.character(.data[[assignment_facet]])
    ) |>
    dplyr::arrange(.data$TemplatePerson, .data$Rater)
  keep_weight <- "Weight" %in% names(prep$data) && any(is.finite(prep$data$Weight) & prep$data$Weight != 1, na.rm = TRUE)
  design_keep <- c("Person", assignment_facet, criterion_facet, intersect("Group", names(prep$data)))
  if (keep_weight) design_keep <- c(design_keep, "Weight")
  design_skeleton <- prep$data |>
    dplyr::select(dplyr::all_of(design_keep)) |>
    dplyr::distinct() |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$Person),
      Rater = as.character(.data[[assignment_facet]]),
      Criterion = as.character(.data[[criterion_facet]]),
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
    model = fit_model
  )
  slopes <- simulation_extract_slopes_from_fit(
    slope_tbl = fit$slopes %||% NULL,
    model = fit_model
  )
  slope_input <- if (identical(fit_model, "GPCM")) {
    slope_match <- match(criterion_levels, slopes$SlopeFacet)
    if (anyNA(slope_match)) {
      stop(
        "`extract_mfrm_sim_spec()` could not align fitted `GPCM` slopes to the observed `",
        criterion_facet,
        "` levels.",
        call. = FALSE
      )
    }
    slopes <- slopes[slope_match, , drop = FALSE]
    as.numeric(slopes$Estimate)
  } else {
    slopes
  }

  spec <- build_mfrm_sim_spec(
    n_person = length(person_levels),
    n_rater = length(rater_levels),
    n_criterion = length(criterion_levels),
    raters_per_person = raters_per_person,
    score_levels = score_levels,
    theta_sd = stats::sd(suppressWarnings(as.numeric(person_tbl$Estimate)), na.rm = TRUE),
    rater_sd = stats::sd(suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == assignment_facet])), na.rm = TRUE),
    criterion_sd = stats::sd(suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == criterion_facet])), na.rm = TRUE),
    noise_sd = 0,
    step_span = if (is.null(thresholds) || nrow(thresholds) == 0) 0 else diff(range(thresholds$Estimate, na.rm = TRUE)) / 2,
    thresholds = thresholds,
    model = fit_model,
    step_facet = simulation_validate_step_facet_name(fit$config$step_facet %||% assignment_facet),
    slope_facet = fit$config$slope_facet %||% NULL,
    slopes = slope_input,
    assignment = assignment,
    latent_distribution = latent_distribution,
    empirical_person = suppressWarnings(as.numeric(person_tbl$Estimate)),
    empirical_rater = suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == assignment_facet])),
    empirical_criterion = suppressWarnings(as.numeric(other_tbl$Estimate[other_tbl$Facet == criterion_facet])),
    assignment_profiles = assignment_profiles,
    design_skeleton = design_skeleton,
    group_levels = if ("Group" %in% names(design_skeleton)) sort(unique(as.character(design_skeleton$Group))) else NULL
  )

  spec$facet_names <- stats::setNames(facet_names, c("rater", "criterion"))
  spec$facet_levels <- list(
    rater = rater_levels,
    criterion = criterion_levels
  )
  if (identical(fit_model, "GPCM")) {
    spec$slope_table <- tibble::as_tibble(slopes)
    spec$slope_regime <- simulation_gpcm_slope_regime(spec$slope_table)
  }
  fit_population <- fit$population %||% list()
  if (isTRUE(fit_population$active)) {
    template <- as.data.frame(fit_population$person_table %||% data.frame(), stringsAsFactors = FALSE)
    template_id <- as.character(fit_population$person_id[1] %||% "Person")
    if (nrow(template) > 0 && template_id %in% names(template)) {
      names(template)[names(template) == template_id] <- "TemplatePerson"
    }
    spec$population <- simulation_build_population_spec(
      population_formula = fit_population$formula,
      population_coefficients = fit_population$coefficients,
      population_sigma2 = fit_population$sigma2,
      population_covariates = template,
      population_xlevels = fit_population$xlevels,
      population_contrasts = fit_population$contrasts
    )
  } else {
    spec$population <- simulation_empty_population_spec()
  }
  spec$source <- "fit_mfrm"
  spec$source_summary <- list(
    observed_raters_per_person = assignment_counts,
    inferred_assignment = inferred_assignment,
    observed_score_values = score_values,
    facet_names = facet_names
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

simulation_default_output_facet_names <- function() {
  c("Rater", "Criterion")
}

simulation_reserved_output_names <- function() {
  c("Study", "Person", "Score", "Group", "Weight", "TemplatePerson")
}

simulation_validate_output_facet_names <- function(facet_names = NULL) {
  if (is.null(facet_names)) {
    return(stats::setNames(simulation_default_output_facet_names(), c("rater", "criterion")))
  }

  value <- facet_names
  if (!is.null(names(value)) && all(c("rater", "criterion") %in% names(value))) {
    value <- value[c("rater", "criterion")]
  }
  value <- unname(as.character(value))

  if (length(value) != 2L || any(is.na(value) | !nzchar(value)) || anyDuplicated(value)) {
    stop(
      "`facet_names` must be a character vector of length 2 naming the rater-like and criterion-like simulated facet columns.",
      call. = FALSE
    )
  }

  reserved <- intersect(value, simulation_reserved_output_names())
  if (length(reserved) > 0L) {
    stop(
      "`facet_names` cannot reuse reserved output names: ",
      paste(reserved, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  stats::setNames(value, c("rater", "criterion"))
}

simulation_validate_step_facet_name <- function(step_facet) {
  value <- as.character(step_facet[1] %||% NA_character_)
  if (is.na(value) || !nzchar(value)) {
    stop("`step_facet` must be a single non-empty character string.", call. = FALSE)
  }
  value
}

simulation_generated_levels <- function(prefix, n_levels) {
  sprintf("%s%02d", prefix, seq_len(n_levels))
}

simulation_generated_role_levels <- function(label, fallback_prefix, n_levels) {
  clean <- gsub("[^A-Za-z0-9]+", "", as.character(label[1] %||% ""))
  prefix <- toupper(substr(clean, 1L, 1L))
  if (!nzchar(prefix)) {
    prefix <- fallback_prefix
  }
  simulation_generated_levels(prefix, n_levels)
}

simulation_spec_output_facet_names <- function(sim_spec = NULL) {
  fallback <- simulation_default_output_facet_names()
  facet_names <- sim_spec$facet_names %||% NULL
  if (is.null(facet_names)) {
    return(fallback)
  }
  facet_names <- unname(as.character(facet_names))
  facet_names <- facet_names[!is.na(facet_names) & nzchar(facet_names)]
  if (length(facet_names) != 2L || anyDuplicated(facet_names)) {
    return(fallback)
  }
  facet_names
}

simulation_step_facet_role <- function(sim_spec, step_facet = sim_spec$step_facet) {
  facet_names <- simulation_spec_output_facet_names(sim_spec)
  match_idx <- match(as.character(step_facet[1] %||% NA_character_), facet_names)
  if (is.na(match_idx)) {
    return(NA_character_)
  }
  c("rater", "criterion")[match_idx]
}

simulation_spec_role_levels <- function(sim_spec, role = c("rater", "criterion"), count = NULL) {
  role <- match.arg(role)
  stored <- sim_spec$facet_levels[[role]] %||% NULL
  stored <- as.character(stored)
  stored <- stored[!is.na(stored) & nzchar(stored)]
  if (!is.null(count) && length(stored) == count && length(stored) > 0L) {
    return(stored)
  }
  if (is.null(count) && length(stored) > 0L) {
    return(stored)
  }
  if (is.null(count)) {
    return(character(0))
  }
  facet_names <- simulation_spec_output_facet_names(sim_spec)
  label <- if (identical(role, "rater")) facet_names[1] else facet_names[2]
  fallback_prefix <- if (identical(role, "rater")) "R" else "C"
  simulation_generated_role_levels(label, fallback_prefix, count)
}

simulation_expected_step_levels <- function(sim_spec, step_facet, n_rater, n_criterion) {
  role <- simulation_step_facet_role(sim_spec, step_facet = step_facet)
  if (identical(role, "rater")) {
    return(simulation_spec_role_levels(sim_spec, "rater", count = n_rater))
  }
  if (identical(role, "criterion")) {
    return(simulation_spec_role_levels(sim_spec, "criterion", count = n_criterion))
  }
  character(0)
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

simulation_empty_population_spec <- function() {
  list(
    active = FALSE,
    formula = NULL,
    coefficients = NULL,
    sigma2 = NULL,
    design_columns = NULL,
    xlevels = NULL,
    contrasts = NULL,
    person_id = NULL,
    covariate_template = NULL,
    notes = "No latent-regression population generator is stored in this simulation specification."
  )
}

simulation_build_population_spec <- function(population_formula = NULL,
                                             population_coefficients = NULL,
                                             population_sigma2 = NULL,
                                             population_covariates = NULL,
                                             population_xlevels = NULL,
                                             population_contrasts = NULL) {
  if (is.null(population_formula)) {
    if (!is.null(population_coefficients) || !is.null(population_sigma2) || !is.null(population_covariates)) {
      stop(
        "`population_coefficients`, `population_sigma2`, and `population_covariates` can be supplied only when `population_formula` is also supplied.",
        call. = FALSE
      )
    }
    return(simulation_empty_population_spec())
  }

  if (!is.data.frame(population_covariates)) {
    stop(
      "`population_covariates` must be a data.frame with one row per template person when `population_formula` is supplied.",
      call. = FALSE
    )
  }

  template <- as.data.frame(population_covariates, stringsAsFactors = FALSE)
  if (!"TemplatePerson" %in% names(template)) {
    template$TemplatePerson <- sprintf("TP%03d", seq_len(nrow(template)))
  }
  template$TemplatePerson <- as.character(template$TemplatePerson)
  if (anyNA(template$TemplatePerson) || any(!nzchar(template$TemplatePerson))) {
    stop("`population_covariates$TemplatePerson` must contain non-missing non-empty IDs.", call. = FALSE)
  }
  if (anyDuplicated(template$TemplatePerson)) {
    stop("`population_covariates` must contain one row per `TemplatePerson`.", call. = FALSE)
  }
  if (nrow(template) < 2L) {
    stop("`population_covariates` must contain at least two template persons.", call. = FALSE)
  }

  scaffold_tbl <- template
  names(scaffold_tbl)[names(scaffold_tbl) == "TemplatePerson"] <- "Person"
  scaffold <- prepare_mfrm_population_scaffold(
    data = data.frame(Person = scaffold_tbl$Person, stringsAsFactors = FALSE),
    person = "Person",
    population_formula = population_formula,
    person_data = scaffold_tbl,
    person_id = "Person",
    population_policy = "error",
    require_full_rank = TRUE,
    population_xlevels = population_xlevels,
    population_contrasts = population_contrasts
  )

  coeff <- suppressWarnings(as.numeric(population_coefficients))
  if (length(coeff) != ncol(scaffold$design_matrix) || any(!is.finite(coeff))) {
    stop(
      "`population_coefficients` must be a finite numeric vector with length ",
      ncol(scaffold$design_matrix),
      " to match the design matrix implied by `population_formula`.",
      call. = FALSE
    )
  }
  coeff_names <- names(population_coefficients %||% NULL)
  if (!is.null(coeff_names) && length(coeff_names) == length(coeff)) {
    if (!setequal(coeff_names, scaffold$design_columns)) {
      stop(
        "Named `population_coefficients` must match the design columns implied by `population_formula`: ",
        paste(scaffold$design_columns, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    coeff <- coeff[match(scaffold$design_columns, coeff_names)]
  }
  names(coeff) <- scaffold$design_columns

  sigma2 <- as.numeric(population_sigma2[1])
  if (!is.finite(sigma2) || sigma2 < 0) {
    stop("`population_sigma2` must be a finite numeric value >= 0.", call. = FALSE)
  }

  list(
    active = TRUE,
    formula = scaffold$formula,
    coefficients = coeff,
    sigma2 = sigma2,
    design_columns = scaffold$design_columns,
    xlevels = scaffold$xlevels,
    contrasts = scaffold$contrasts,
    person_id = "TemplatePerson",
    covariate_template = template,
    notes = c(
      "This simulation specification stores a first-version latent-regression person generator.",
      "Template person rows are resampled to the requested design size before computing theta = X beta + e."
    )
  )
}

simulation_normalize_assignment_profiles <- function(assignment_profiles,
                                                     assignment,
                                                     n_rater,
                                                     facet_names = simulation_default_output_facet_names()) {
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
  assignment_facet <- as.character(facet_names[1] %||% "Rater")
  rater_col <- c(assignment_facet, "Rater")
  rater_col <- unique(rater_col[rater_col %in% names(tbl)])[1] %||% NA_character_
  if (!"TemplatePerson" %in% names(tbl) || is.na(rater_col) || !nzchar(rater_col)) {
    stop(
      "`assignment_profiles` must include `TemplatePerson` and `", assignment_facet,
      "` columns. The canonical `Rater` name is also accepted.",
      call. = FALSE
    )
  }

  tbl <- tbl |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$TemplatePerson),
      Rater = as.character(.data[[rater_col]]),
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
                                                 n_criterion,
                                                 facet_names = simulation_default_output_facet_names()) {
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
  assignment_facet <- as.character(facet_names[1] %||% "Rater")
  criterion_facet <- as.character(facet_names[2] %||% "Criterion")
  rater_col <- c(assignment_facet, "Rater")
  rater_col <- unique(rater_col[rater_col %in% names(tbl)])[1] %||% NA_character_
  criterion_col <- c(criterion_facet, "Criterion")
  criterion_col <- unique(criterion_col[criterion_col %in% names(tbl)])[1] %||% NA_character_
  if (!"TemplatePerson" %in% names(tbl) ||
      is.na(rater_col) || !nzchar(rater_col) ||
      is.na(criterion_col) || !nzchar(criterion_col)) {
    stop(
      "`design_skeleton` must include `TemplatePerson`, `", assignment_facet,
      "`, and `", criterion_facet, "` columns. Canonical `Rater` / `Criterion` names are also accepted.",
      call. = FALSE
    )
  }

  keep_group <- "Group" %in% names(tbl)
  keep_weight <- "Weight" %in% names(tbl)
  keep_reuse <- "TemplatePersonReuse" %in% names(tbl)
  tbl <- tbl |>
    dplyr::transmute(
      TemplatePerson = as.character(.data$TemplatePerson),
      Rater = as.character(.data[[rater_col]]),
      Criterion = as.character(.data[[criterion_col]]),
      Group = if (keep_group) as.character(.data$Group) else NA_character_,
      Weight = if (keep_weight) .data$Weight else NA_real_,
      TemplatePersonReuse = if (keep_reuse) .data$TemplatePersonReuse else FALSE
    ) |>
    dplyr::mutate(
      Weight = suppressWarnings(as.numeric(.data$Weight)),
      TemplatePersonReuse = simulation_sparse_design_active(.data$TemplatePersonReuse)
    ) |>
    dplyr::filter(!is.na(.data$TemplatePerson), nzchar(.data$TemplatePerson),
                  !is.na(.data$Rater), nzchar(.data$Rater),
                  !is.na(.data$Criterion), nzchar(.data$Criterion)) |>
    dplyr::distinct()
  if (!keep_group) {
    tbl <- dplyr::select(tbl, -dplyr::all_of("Group"))
  }
  if (!keep_weight) {
    tbl <- dplyr::select(tbl, -dplyr::all_of("Weight"))
  }
  if (!keep_reuse) {
    tbl <- dplyr::select(tbl, -dplyr::all_of("TemplatePersonReuse"))
  }

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

simulation_normalize_sparse_controls <- function(sparse_controls,
                                                 assignment,
                                                 n_person,
                                                 n_rater,
                                                 raters_per_person) {
  if (!identical(assignment, "sparse_linked")) {
    return(list(active = FALSE))
  }
  if (is.null(sparse_controls)) {
    sparse_controls <- list()
  }
  if (!is.list(sparse_controls)) {
    stop("`sparse_controls` must be a named list when `assignment = \"sparse_linked\"`.",
         call. = FALSE)
  }

  link_fraction <- suppressWarnings(as.numeric((sparse_controls$link_fraction %||% 0.1)[1]))
  if (!is.finite(link_fraction) || link_fraction < 0 || link_fraction > 1) {
    stop("`sparse_controls$link_fraction` must be a finite number in [0, 1].",
         call. = FALSE)
  }
  link_persons_source <- as.character(sparse_controls$link_persons_source %||% NA_character_)[1]
  use_link_fraction <- is.null(sparse_controls$link_persons) ||
    identical(link_persons_source, "link_fraction")
  link_persons <- if (use_link_fraction) {
    as.integer(ceiling(n_person * link_fraction))
  } else {
    suppressWarnings(as.integer(sparse_controls$link_persons[1]))
  }
  if (!is.finite(link_persons) || link_persons < 0L || link_persons > n_person) {
    stop("`sparse_controls$link_persons` must be an integer between 0 and `n_person`.",
         call. = FALSE)
  }

  link_raters_source <- as.character(sparse_controls$link_raters_per_person_source %||% NA_character_)[1]
  use_all_raters_link <- is.null(sparse_controls$link_raters_per_person) ||
    identical(link_raters_source, "n_rater")
  link_raters_per_person <- if (use_all_raters_link) {
    as.integer(n_rater)
  } else {
    suppressWarnings(as.integer(sparse_controls$link_raters_per_person[1]))
  }
  if (!is.finite(link_raters_per_person) || link_raters_per_person < 1L ||
      link_raters_per_person > n_rater) {
    stop("`sparse_controls$link_raters_per_person` must be an integer between 1 and `n_rater`.",
         call. = FALSE)
  }
  if (link_persons > 0L && link_raters_per_person < raters_per_person) {
    stop("`sparse_controls$link_raters_per_person` must be >= `raters_per_person` when linking persons are used.",
         call. = FALSE)
  }

  assignment_mode <- tolower(as.character(sparse_controls$assignment_mode %||% "balanced")[1])
  assignment_mode <- match.arg(assignment_mode, c("balanced", "random"))
  min_common <- suppressWarnings(as.integer(
    (sparse_controls$min_common_persons_per_rater_pair %||% 1L)[1]
  ))
  if (!is.finite(min_common) || min_common < 0L) {
    stop("`sparse_controls$min_common_persons_per_rater_pair` must be an integer >= 0.",
         call. = FALSE)
  }

  list(
    active = TRUE,
    link_fraction = as.numeric(link_fraction),
    link_persons = as.integer(link_persons),
    link_persons_source = if (use_link_fraction) "link_fraction" else "link_persons",
    link_raters_per_person = as.integer(link_raters_per_person),
    link_raters_per_person_source = if (use_all_raters_link) "n_rater" else "link_raters_per_person",
    assignment_mode = assignment_mode,
    min_common_persons_per_rater_pair = as.integer(min_common),
    notes = c(
      "Sparse linked simulation creates planned missingness through incomplete rater assignment.",
      "Linking persons are a design device for preserving common rater links, not a separate population."
    )
  )
}

simulation_threshold_list_to_table <- function(thresholds) {
  if (length(thresholds) == 0L) {
    stop("`thresholds` list must contain at least one threshold set.", call. = FALSE)
  }

  step_facets <- names(thresholds)
  if (is.null(step_facets)) {
    if (length(thresholds) == 1L) {
      step_facets <- "Common"
    } else {
      stop("`thresholds` lists with multiple entries must name each `StepFacet`.", call. = FALSE)
    }
  } else {
    step_facets <- as.character(step_facets)
    missing_names <- !nzchar(step_facets)
    if (length(thresholds) == 1L && all(missing_names)) {
      step_facets <- "Common"
    } else if (any(missing_names)) {
      stop("Every entry in a `thresholds` list must have a non-empty `StepFacet` name.", call. = FALSE)
    }
  }
  if (anyDuplicated(step_facets)) {
    stop("`thresholds` list names must be unique `StepFacet` labels.", call. = FALSE)
  }

  pieces <- vector("list", length(thresholds))
  for (i in seq_along(thresholds)) {
    est <- thresholds[[i]]
    if (!is.numeric(est) || is.matrix(est) || length(dim(est)) > 0L) {
      stop("Each entry in a `thresholds` list must be a numeric vector.", call. = FALSE)
    }
    pieces[[i]] <- tibble::tibble(
      StepFacet = step_facets[[i]],
      StepIndex = seq_along(est),
      Step = paste0("Step_", seq_along(est)),
      Estimate = as.numeric(est)
    )
  }

  dplyr::bind_rows(pieces)
}

simulation_threshold_matrix_to_table <- function(thresholds) {
  if (!is.numeric(thresholds)) {
    stop("`thresholds` matrices must be numeric.", call. = FALSE)
  }
  if (nrow(thresholds) == 0L || ncol(thresholds) == 0L) {
    stop("`thresholds` matrices must have at least one row and one column.", call. = FALSE)
  }

  step_facets <- rownames(thresholds)
  if (is.null(step_facets)) {
    if (nrow(thresholds) == 1L) {
      step_facets <- "Common"
    } else {
      stop("`thresholds` matrices with multiple rows must have row names for `StepFacet` labels.", call. = FALSE)
    }
  } else {
    step_facets <- as.character(step_facets)
    if (any(!nzchar(step_facets))) {
      stop("Every row in a `thresholds` matrix must have a non-empty `StepFacet` name.", call. = FALSE)
    }
  }
  if (anyDuplicated(step_facets)) {
    stop("`thresholds` matrix row names must be unique `StepFacet` labels.", call. = FALSE)
  }

  step_index <- rep(seq_len(ncol(thresholds)), times = nrow(thresholds))
  tibble::tibble(
    StepFacet = rep(step_facets, each = ncol(thresholds)),
    StepIndex = step_index,
    Step = paste0("Step_", step_index),
    Estimate = as.numeric(t(thresholds))
  )
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

  if (is.matrix(thresholds)) {
    thresholds <- simulation_threshold_matrix_to_table(thresholds)
  } else if (is.list(thresholds) && !is.data.frame(thresholds)) {
    thresholds <- simulation_threshold_list_to_table(thresholds)
  } else if (is.numeric(thresholds)) {
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
    stop("`thresholds` must be NULL, a numeric vector, a named list, a numeric matrix, or a data.frame.", call. = FALSE)
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

simulation_expected_slope_levels <- function(facet_names, slope_facet, n_rater, n_criterion) {
  if (is.null(slope_facet) || !nzchar(as.character(slope_facet[1]))) {
    return(character(0))
  }
  slope_facet <- as.character(slope_facet[1])
  if (identical(slope_facet, facet_names[["rater"]])) {
    return(simulation_generated_role_levels(facet_names[["rater"]], "R", n_rater))
  }
  if (identical(slope_facet, facet_names[["criterion"]])) {
    return(simulation_generated_role_levels(facet_names[["criterion"]], "C", n_criterion))
  }
  character(0)
}

simulation_build_slope_table <- function(slopes,
                                         model,
                                         slope_facet,
                                         facet_names,
                                         n_rater,
                                         n_criterion) {
  if (!identical(model, "GPCM")) {
    if (!is.null(slopes)) {
      stop("`slopes` can be supplied only when `model = \"GPCM\"`.", call. = FALSE)
    }
    return(NULL)
  }

  expected_levels <- simulation_expected_slope_levels(
    facet_names = facet_names,
    slope_facet = slope_facet,
    n_rater = n_rater,
    n_criterion = n_criterion
  )
  if (length(expected_levels) == 0L) {
    stop("Could not resolve expected slope-facet levels for bounded `GPCM`.", call. = FALSE)
  }

  normalize_slopes <- function(x) {
    x <- as.numeric(x)
    if (any(!is.finite(x)) || any(x <= 0)) {
      stop("`GPCM` slopes must contain strictly positive finite values.", call. = FALSE)
    }
    exp(log(x) - mean(log(x)))
  }

  if (is.null(slopes)) {
    return(tibble::tibble(
      SlopeFacet = expected_levels,
      Estimate = rep(1, length(expected_levels))
    ))
  }

  if (is.numeric(slopes)) {
    slope_vals <- as.numeric(slopes)
    if (!is.null(names(slopes)) && setequal(names(slopes), expected_levels)) {
      slope_vals <- slope_vals[match(expected_levels, names(slopes))]
    }
    if (length(slope_vals) != length(expected_levels)) {
      stop(
        "Numeric `slopes` must have length equal to the number of `", slope_facet,
        "` levels (", length(expected_levels), ").",
        call. = FALSE
      )
    }
    slope_vals <- normalize_slopes(slope_vals)
    return(tibble::tibble(
      SlopeFacet = expected_levels,
      Estimate = slope_vals
    ))
  }

  if (!is.data.frame(slopes)) {
    stop("`slopes` must be NULL, a numeric vector, or a data.frame.", call. = FALSE)
  }

  tbl <- tibble::as_tibble(slopes)
  if (!all(c("SlopeFacet", "Estimate") %in% names(tbl))) {
    stop("Slope data frames must include `SlopeFacet` and `Estimate` columns.", call. = FALSE)
  }
  tbl <- tbl |>
    dplyr::transmute(
      SlopeFacet = as.character(.data$SlopeFacet),
      Estimate = suppressWarnings(as.numeric(.data$Estimate))
    ) |>
    dplyr::filter(!is.na(.data$SlopeFacet), nzchar(.data$SlopeFacet)) |>
    dplyr::distinct(.data$SlopeFacet, .keep_all = TRUE) |>
    dplyr::arrange(.data$SlopeFacet)
  if (nrow(tbl) == 0L) {
    stop("`slopes` did not contain any valid slope rows.", call. = FALSE)
  }
  tbl$Estimate <- normalize_slopes(tbl$Estimate)
  if (!setequal(tbl$SlopeFacet, expected_levels)) {
    stop(
      "`slopes` must cover exactly the generated `", slope_facet, "` levels: ",
      paste(expected_levels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  tbl[match(expected_levels, tbl$SlopeFacet), , drop = FALSE]
}

simulation_gpcm_slope_regime <- function(slope_table,
                                         near_flat_log = log(1.05),
                                         high_dispersion_log = log(1.50)) {
  # Operational simulation stress bins, not literature-derived model-fit cutoffs.
  if (!is.data.frame(slope_table) || nrow(slope_table) == 0L ||
      !"Estimate" %in% names(slope_table)) {
    return(NA_character_)
  }
  slopes <- suppressWarnings(as.numeric(slope_table$Estimate))
  slopes <- slopes[is.finite(slopes) & slopes > 0]
  if (length(slopes) == 0L) return(NA_character_)
  log_slopes <- log(slopes) - mean(log(slopes))
  max_abs_log <- max(abs(log_slopes), na.rm = TRUE)
  if (!is.finite(max_abs_log)) return(NA_character_)
  if (max_abs_log <= sqrt(.Machine$double.eps)) return("unit_slopes")
  if (max_abs_log <= near_flat_log) return("near_flat")
  if (max_abs_log <= high_dispersion_log) return("moderate")
  "high_dispersion"
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

simulation_extract_slopes_from_fit <- function(slope_tbl, model) {
  if (!identical(model, "GPCM")) {
    return(NULL)
  }
  if (!is.data.frame(slope_tbl) || nrow(slope_tbl) == 0) {
    stop(
      "`extract_mfrm_sim_spec()` requires a non-empty `fit$slopes` table for bounded `GPCM` fits.",
      call. = FALSE
    )
  }
  tbl <- tibble::as_tibble(slope_tbl)
  if (!all(c("SlopeFacet", "Estimate") %in% names(tbl))) {
    stop(
      "`extract_mfrm_sim_spec()` requires `fit$slopes` to contain `SlopeFacet` and `Estimate` columns.",
      call. = FALSE
    )
  }
  tbl |>
    dplyr::transmute(
      SlopeFacet = as.character(.data$SlopeFacet),
      Estimate = as.numeric(.data$Estimate)
    ) |>
    dplyr::arrange(.data$SlopeFacet)
}

simulation_spec_slope_mode <- function(sim_spec) {
  tbl <- sim_spec$slope_table %||% NULL
  if (!is.data.frame(tbl) || nrow(tbl) == 0) {
    return("none")
  }
  "slope_facet_specific"
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
  if (identical(base_assignment, "sparse_linked")) {
    return("sparse_linked")
  }
  if (isTRUE(raters_per_person >= n_rater) && identical(base_assignment, "crossed")) {
    return("crossed")
  }
  "rotating"
}

simulation_override_spec_design <- function(sim_spec,
                                            n_person = NULL,
                                            n_rater = NULL,
                                            n_criterion = NULL,
                                            raters_per_person = NULL,
                                            design = NULL,
                                            group_levels = sim_spec$group_levels,
                                            dif_effects = sim_spec$dif_effects,
                                            interaction_effects = sim_spec$interaction_effects) {
  if (!inherits(sim_spec, "mfrm_sim_spec")) {
    stop("`sim_spec` must inherit from `mfrm_sim_spec`.", call. = FALSE)
  }

  design_counts <- simulation_resolve_design_counts(
    sim_spec = sim_spec,
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    defaults = list(
      n_person = sim_spec$n_person,
      n_rater = sim_spec$n_rater,
      n_criterion = sim_spec$n_criterion,
      raters_per_person = sim_spec$raters_per_person
    ),
    design_arg = "design"
  )
  n_person <- design_counts$n_person
  n_rater <- design_counts$n_rater
  n_criterion <- design_counts$n_criterion
  raters_per_person <- design_counts$raters_per_person

  out <- sim_spec
  out$n_person <- n_person
  out$n_rater <- n_rater
  out$n_criterion <- n_criterion
  out$raters_per_person <- raters_per_person
  out$assignment <- simulation_resolve_assignment(sim_spec$assignment, n_rater = n_rater, raters_per_person = raters_per_person)
  out$group_levels <- group_levels
  out$dif_effects <- dif_effects
  out$interaction_effects <- interaction_effects
  out$sparse_controls <- simulation_normalize_sparse_controls(
    sparse_controls = sim_spec$sparse_controls %||% NULL,
    assignment = out$assignment,
    n_person = n_person,
    n_rater = n_rater,
    raters_per_person = raters_per_person
  )
  out$planning_constraints <- simulation_planning_constraints(out)
  out$planning_schema <- simulation_planning_schema(out)

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
    expected_levels <- simulation_expected_step_levels(
      out,
      step_facet = out$step_facet,
      n_rater = n_rater,
      n_criterion = n_criterion
    )
    observed_levels <- sort(unique(as.character(out$threshold_table$StepFacet)))
    if (!setequal(observed_levels, expected_levels)) {
      role <- simulation_step_facet_role(out, step_facet = out$step_facet)
      varying_arg <- if (identical(role, "criterion")) "n_criterion" else if (identical(role, "rater")) "n_rater" else "design dimensions"
      stop(
        "`sim_spec` contains step-facet-specific thresholds for `", out$step_facet,
        "` levels {", paste(observed_levels, collapse = ", "), "}. ",
        "Varying `", varying_arg, "` away from that threshold structure is not currently supported. ",
        "Use common thresholds or build a design-specific simulation specification.",
        call. = FALSE
      )
    }
  }

  slope_mode <- simulation_spec_slope_mode(out)
  if (identical(slope_mode, "slope_facet_specific")) {
    expected_levels <- simulation_expected_step_levels(
      out,
      step_facet = out$slope_facet,
      n_rater = n_rater,
      n_criterion = n_criterion
    )
    observed_levels <- sort(unique(as.character(out$slope_table$SlopeFacet)))
    if (!setequal(observed_levels, expected_levels)) {
      role <- simulation_step_facet_role(out, step_facet = out$slope_facet)
      varying_arg <- if (identical(role, "criterion")) "n_criterion" else if (identical(role, "rater")) "n_rater" else "design dimensions"
      stop(
        "`sim_spec` contains `GPCM` slopes for `", out$slope_facet,
        "` levels {", paste(observed_levels, collapse = ", "), "}. ",
        "Varying `", varying_arg, "` away from that slope structure is not currently supported. ",
        "Build a design-specific simulation specification instead.",
        call. = FALSE
      )
    }
  }

  out$planning_constraints <- simulation_planning_constraints(out)
  out$planning_schema <- simulation_planning_schema(out)

  out
}
