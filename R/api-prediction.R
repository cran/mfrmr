#' Forecast population-level MFRM operating characteristics for one future design
#'
#' @param fit Optional output from [fit_mfrm()] used to derive a fit-based
#'   simulation specification.
#' @param sim_spec Optional output from [build_mfrm_sim_spec()] or
#'   [extract_mfrm_sim_spec()]. Supply exactly one of `fit` or `sim_spec`.
#' @param n_person Number of persons/respondents in the future design. Defaults
#'   to the value stored in the base simulation specification.
#' @param n_rater Number of rater facet levels in the future design. Defaults to
#'   the value stored in the base simulation specification.
#' @param n_criterion Number of criterion/item facet levels in the future
#'   design. Defaults to the value stored in the base simulation specification.
#' @param raters_per_person Number of raters assigned to each person in the
#'   future design. Defaults to the value stored in the base simulation
#'   specification.
#' @param design Optional named design override supplied as a named list,
#'   named vector, or one-row data frame. Names may use canonical variables
#'   (`n_person`, `n_rater`, `n_criterion`, `raters_per_person`), current
#'   public aliases (for example `n_judge`, `n_task`, `judge_per_person`), or
#'   role keywords (`person`, `rater`, `criterion`, `assignment`). The
#'   schema-only future branch input `design$facets = c(person = ..., judge =
#'   ..., task = ...)` is also accepted for the currently exposed facet keys.
#'   Do not specify the same variable through both `design` and the scalar
#'   count arguments.
#' @param reps Number of replications used in the forecast simulation.
#' @param fit_method Estimation method used inside the forecast simulation. When
#'   `fit` is supplied, defaults to that fit's estimation method; otherwise
#'   defaults to `"MML"`.
#' @param model Measurement model used when refitting the forecasted design.
#'   Defaults to the model recorded in the base simulation specification.
#' @param maxit Maximum iterations passed to [fit_mfrm()] in each replication.
#' @param quad_points Quadrature points for `fit_method = "MML"`.
#' @param residual_pca Residual PCA mode passed to [diagnose_mfrm()].
#' @param seed Optional seed for reproducible replications.
#'
#' @details
#' `predict_mfrm_population()` is a **scenario-level forecasting helper** built
#' on top of [evaluate_mfrm_design()]. It is intended for questions such as:
#' - what separation/reliability would we expect if the next administration had
#'   60 persons, 4 raters, and 2 ratings per person?
#' - how much Monte Carlo uncertainty remains around those expected summaries?
#'
#' The function deliberately returns **aggregate operating characteristics**
#' (for example mean separation, reliability, recovery RMSE, convergence rate)
#' rather than future individual true values for one respondent or one rater.
#'
#' If `fit` is supplied, the function first constructs a fit-derived parametric
#' starting point with [extract_mfrm_sim_spec()] and then evaluates the
#' requested future design under that explicit data-generating mechanism. This
#' should be interpreted as a fit-based forecast under modeling assumptions, not
#' as a guaranteed out-of-sample prediction.
#'
#' When that fit-derived or manually built simulation specification stores an
#' active latent-regression population generator, the helper still operates at
#' the **design / operating-characteristic** level. It repeatedly simulates
#' person-level covariates and responses, refits the MML population-model
#' branch, and summarizes the resulting facet-level behavior. This is distinct
#' from the fitted-model posterior scoring provided by [predict_mfrm_units()].
#'
#' Bounded `GPCM` forecasts are available with caveats through the same
#' repeated simulation/refit design route used by [evaluate_mfrm_design()].
#' They summarize design-level operating characteristics under the supplied or
#' fit-derived slope-aware specification; they do not validate operational
#' scoring, diagnostic-screening or signal-detection rules, slope-recovery
#' adequacy, or arbitrary-facet planning.
#'
#' @section Interpreting output:
#' - `forecast` contains facet-level expected summaries for the requested
#'   future design.
#' - `Mcse*` columns quantify Monte Carlo uncertainty from using a finite number
#'   of replications.
#' - `design_variable_aliases` and `design_descriptor` carry the same public
#'   naming metadata used by the underlying planning object. They rename the
#'   standard two non-person facet roles for presentation, but they do not turn
#'   the current planner into a fully arbitrary-facet simulator.
#' - If `sim_spec$population$active = TRUE`, the forecast summarizes repeated
#'   latent-regression MML refits under that stored person-level generator; it
#'   is still a scenario forecast rather than direct posterior scoring for one
#'   observed sample.
#' - `simulation` stores the full design-evaluation object in case you want to
#'   inspect replicate-level behavior.
#'
#' @section What this does not justify:
#' This helper does not produce definitive future person measures or rater
#' severities for one concrete sample. It forecasts design-level behavior under
#' the supplied or derived parametric assumptions.
#'
#' @section References:
#' The forecast is implemented as a one-scenario Monte Carlo / operating-
#' characteristic study following the general guidance of Morris, White, and
#' Crowther (2019) and the ADEMP-oriented reporting framework discussed by
#' Siepe et al. (2024). In `mfrmr`, this function is a practical wrapper for
#' future-design planning rather than a direct implementation of a published
#' many-facet forecasting procedure.
#'
#' - Morris, T. P., White, I. R., & Crowther, M. J. (2019).
#'   *Using simulation studies to evaluate statistical methods*.
#'   Statistics in Medicine, 38(11), 2074-2102.
#' - Siepe, B. S., Bartos, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W.,
#'   & Pawel, S. (2024). *Simulation studies for methodological research in
#'   psychology: A standardized template for planning, preregistration, and
#'   reporting*. Psychological Methods.
#'
#' @return An object of class `mfrm_population_prediction` with components:
#' - `design`: requested future design
#' - `forecast`: facet-level forecast table
#' - `overview`: run-level overview
#' - `simulation`: underlying [evaluate_mfrm_design()] result
#' - `sim_spec`: simulation specification used for the forecast
#' - `facet_names`: public non-person facet names carried by the simulation
#'   specification
#' - `design_variable_aliases`: public aliases for
#'   `n_person`/`n_rater`/`n_criterion`/`raters_per_person`
#' - `design_descriptor`: role-based description of design variables carried
#'   from the underlying planning object
#' - `planning_scope`: explicit record of the current planning contract,
#'   including a `facet_manifest` and future-planner scaffold marker
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract carrying the role
#'   table, current boundary, mutability map, facet manifest, and a
#'   schema-only future facet-count table
#' - `gpcm_boundary`: bounded-`GPCM` caveat row when a `GPCM` forecast route is
#'   used
#' - `settings`: forecasting settings
#' - `ademp`: simulation-study metadata
#' - `notes`: interpretation notes
#' @seealso [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()],
#'   [evaluate_mfrm_design()], [summary.mfrm_population_prediction]
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   n_person = 16,
#'   n_rater = 3,
#'   n_criterion = 2,
#'   raters_per_person = 2,
#'   assignment = "rotating"
#' )
#' pred <- predict_mfrm_population(
#'   sim_spec = spec,
#'   design = list(person = 18),
#'   reps = 1,
#'   maxit = 30,
#'   seed = 123
#' )
#' s_pred <- summary(pred)
#' s_pred$forecast[, c("Facet", "MeanSeparation", "McseSeparation")]
#' }
#' @export
predict_mfrm_population <- function(fit = NULL,
                                    sim_spec = NULL,
                                    n_person = NULL,
                                    n_rater = NULL,
                                    n_criterion = NULL,
                                    raters_per_person = NULL,
                                    design = NULL,
                                    reps = 50,
                                    fit_method = NULL,
                                    model = NULL,
                                    maxit = 25,
                                    quad_points = 7,
                                    residual_pca = c("none", "overall", "facet", "both"),
                                    seed = NULL) {
  residual_pca <- match.arg(residual_pca)
  has_fit <- !is.null(fit)
  has_spec <- !is.null(sim_spec)
  if (identical(has_fit, has_spec)) {
    stop("Supply exactly one of `fit` or `sim_spec`.", call. = FALSE)
  }

  if (has_fit) {
    if (!inherits(fit, "mfrm_fit")) {
      stop("`fit` must be output from fit_mfrm().", call. = FALSE)
    }
    default_fit_method <- as.character(fit$summary$Method[1])
    if (!is.character(default_fit_method) || !nzchar(default_fit_method)) {
      default_fit_method <- "MML"
    }
    if (identical(default_fit_method, "JMLE")) default_fit_method <- "JML"
    default_model <- as.character(fit$summary$Model[1])
    if (!is.character(default_model) || !nzchar(default_model)) {
      default_model <- as.character(fit$config$model)
    }
    if (identical(default_model, "GPCM")) {
      # GPCM fit-derived forecasts are caveated design-level simulations, not
      # operational score predictions.
    }
    base_spec <- extract_mfrm_sim_spec(fit)
  } else {
    if (!inherits(sim_spec, "mfrm_sim_spec")) {
      stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
    }
    base_spec <- sim_spec
    default_fit_method <- "MML"
    default_model <- as.character(base_spec$model)
    if (identical(default_model, "GPCM")) {
      # GPCM sim-spec forecasts are caveated design-level simulations.
    }
  }

  fit_method <- toupper(as.character(fit_method[1] %||% default_fit_method))
  fit_method <- match.arg(fit_method, c("JML", "MML"))
  model <- toupper(as.character(model[1] %||% default_model))
  model <- match.arg(model, c("RSM", "PCM", "GPCM"))

  design <- simulation_resolve_design_counts(
    sim_spec = base_spec,
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    defaults = list(
      n_person = base_spec$n_person,
      n_rater = base_spec$n_rater,
      n_criterion = base_spec$n_criterion,
      raters_per_person = base_spec$raters_per_person
    ),
    design_arg = "design"
  )

  forecast_spec <- simulation_override_spec_design(
    base_spec,
    design = as.list(design[1, , drop = FALSE])
  )

  sim_eval <- evaluate_mfrm_design(
    n_person = design$n_person,
    n_rater = design$n_rater,
    n_criterion = design$n_criterion,
    raters_per_person = design$raters_per_person,
    reps = reps,
    fit_method = fit_method,
    model = model,
    maxit = maxit,
    quad_points = quad_points,
    residual_pca = residual_pca,
    sim_spec = forecast_spec,
    seed = seed
  )

  sim_summary <- summary(sim_eval, digits = 6)
  design_variable_aliases <- simulation_object_design_variable_aliases(sim_eval)
  design_descriptor <- simulation_object_design_descriptor(sim_eval)
  planning_scope <- simulation_object_planning_scope(sim_eval)
  planning_constraints <- simulation_object_planning_constraints(sim_eval)
  planning_schema <- simulation_object_planning_schema(sim_eval)
  gpcm_boundary <- sim_eval$gpcm_boundary %||% data.frame()
  facet_names <- simulation_spec_output_facet_names(forecast_spec)
  design_public <- simulation_append_design_alias_columns(design, design_variable_aliases)
  notes <- c(
    "This forecast summarizes expected design-level behavior under the supplied or fit-derived simulation specification.",
    "MCSE columns quantify Monte Carlo uncertainty from using a finite number of replications.",
    "Do not interpret this output as deterministic future person/rater true values."
  )
  notes <- unique(c(notes, sim_eval$notes %||% character(0)))
  scope_note <- simulation_planning_scope_note(planning_scope)
  if (length(scope_note) > 0L && !scope_note %in% notes) {
    notes <- c(notes, scope_note)
  }
  constraint_note <- simulation_planning_constraints_note(planning_constraints)
  if (length(constraint_note) > 0L && !constraint_note %in% notes) {
    notes <- c(notes, constraint_note)
  }
  schema_note <- simulation_planning_schema_note(planning_schema)
  if (length(schema_note) > 0L && !schema_note %in% notes) {
    notes <- c(notes, schema_note)
  }

  structure(
    list(
      design = design_public,
      forecast = tibble::as_tibble(sim_summary$design_summary),
      overview = tibble::as_tibble(sim_summary$overview),
      simulation = sim_eval,
      sim_spec = forecast_spec,
      facet_names = facet_names,
      design_variable_aliases = design_variable_aliases,
      design_descriptor = design_descriptor,
      planning_scope = planning_scope,
      planning_constraints = planning_constraints,
      planning_schema = planning_schema,
      gpcm_boundary = gpcm_boundary,
      settings = list(
        reps = as.integer(reps[1]),
        fit_method = fit_method,
        model = model,
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        source = if (has_fit) "fit_mfrm" else "mfrm_sim_spec",
        seed = seed,
        facet_names = facet_names,
        design_variable_aliases = design_variable_aliases,
        design_descriptor = design_descriptor,
        planning_scope = planning_scope,
        planning_constraints = planning_constraints,
        planning_schema = planning_schema,
        gpcm_design_status = if (identical(model, "GPCM") ||
          identical(as.character(forecast_spec$model %||% NA_character_), "GPCM")) {
          "supported_with_caveat"
        } else {
          NA_character_
        }
      ),
      ademp = sim_eval$ademp,
      notes = notes
    ),
    class = "mfrm_population_prediction"
  )
}

#' Summarize a population-level design forecast
#'
#' @param object Output from [predict_mfrm_population()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_population_prediction` with:
#' - `design`: requested future design
#' - `overview`: run-level overview
#' - `forecast`: facet-level forecast table
#' - `facet_names`: public non-person facet names used in the forecast
#' - `design_variable_aliases`: public aliases for design variables
#' - `design_descriptor`: role-based description of design variables
#' - `planning_scope`: explicit record of the current planning contract
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: combined planner-schema contract
#' - `gpcm_boundary`: bounded-`GPCM` caveat row when present
#' - `future_branch_active_summary`: compact deterministic summary of the
#'   schema-only future arbitrary-facet planning branch embedded in the current
#'   planning schema
#' - `ademp`: simulation-study metadata
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_population()]
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   n_person = 16,
#'   n_rater = 3,
#'   n_criterion = 2,
#'   raters_per_person = 2,
#'   assignment = "rotating"
#' )
#' pred <- predict_mfrm_population(
#'   sim_spec = spec,
#'   design = list(person = 18),
#'   reps = 1,
#'   maxit = 30,
#'   seed = 123
#' )
#' s <- summary(pred)
#' s$overview
#' s$forecast[, c("Facet", "MeanSeparation", "McseSeparation")]
#' }
#' @method summary mfrm_population_prediction
#' @export
summary.mfrm_population_prediction <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_population_prediction")) {
    stop("`object` must be output from predict_mfrm_population().", call. = FALSE)
  }
  digits <- prediction_validate_integer(digits, "digits", min_value = 0L, positive = FALSE)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out <- list(
    design = round_df(object$design),
    overview = round_df(object$overview),
    forecast = round_df(object$forecast),
    facet_names = object$facet_names %||% object$settings$facet_names %||% simulation_spec_output_facet_names(object$sim_spec),
    design_variable_aliases = object$design_variable_aliases %||% object$settings$design_variable_aliases %||% simulation_design_variable_aliases(object$sim_spec),
    design_descriptor = object$design_descriptor %||% object$settings$design_descriptor %||% simulation_design_descriptor(object$sim_spec),
    planning_scope = object$planning_scope %||% object$settings$planning_scope %||% simulation_planning_scope(object$sim_spec),
    planning_constraints = object$planning_constraints %||% object$settings$planning_constraints %||% simulation_planning_constraints(object$sim_spec),
    planning_schema = object$planning_schema %||% object$settings$planning_schema %||% simulation_planning_schema(object$sim_spec),
    gpcm_boundary = object$gpcm_boundary %||% data.frame(),
    future_branch_active_summary = simulation_compact_future_branch_active_summary(
      object,
      digits = digits
    ),
    ademp = object$ademp %||% NULL,
    notes = object$notes %||% character(0),
    digits = digits
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
  if (inherits(out$future_branch_active_summary, "summary.mfrm_future_branch_active_branch")) {
    out$notes <- c(
      out$notes,
      "A deterministic future arbitrary-facet planning scaffold is embedded in `future_branch_active_summary`; it reports structural bookkeeping and conservative recommendation logic, not forecast uncertainty."
    )
  }
  out$notes <- unique(out$notes)
  class(out) <- "summary.mfrm_population_prediction"
  out
}

#' @export
print.summary.mfrm_population_prediction <- function(x, ...) {
  digits <- prediction_validate_integer(x$digits %||% 3L, "digits", min_value = 0L, positive = FALSE)
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

  cat("mfrmr Population Prediction Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_df(as.data.frame(x$overview)), row.names = FALSE)
  }
  if (!is.null(x$design) && nrow(x$design) > 0) {
    cat("\nDesign grid (preview)\n")
    print(round_df(as.data.frame(preview_df(x$design))), row.names = FALSE)
  }
  if (!is.null(x$forecast) && nrow(x$forecast) > 0) {
    cat("\nForecast (preview)\n")
    print(round_df(as.data.frame(preview_df(x$forecast))), row.names = FALSE)
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

resolve_prediction_facets <- function(fit, facets = NULL) {
  fit_facets <- as.character(fit$config$facet_names %||% character(0))
  if (length(fit_facets) == 0) {
    stop("`fit` does not contain any calibrated facet names.", call. = FALSE)
  }

  facets <- facets %||% (fit$config$source_columns$facets %||% fit_facets)
  if (!is.character(facets) || length(facets) != length(fit_facets)) {
    stop("`facets` must be a character vector naming one column per calibrated facet: ",
         paste(fit_facets, collapse = ", "), ".", call. = FALSE)
  }

  if (!is.null(names(facets)) && any(nzchar(names(facets)))) {
    if (!setequal(names(facets), fit_facets)) {
      stop("Named `facets` must use the calibrated facet names: ",
           paste(fit_facets, collapse = ", "), ".", call. = FALSE)
    }
    facets <- facets[fit_facets]
  } else {
    names(facets) <- fit_facets
  }

  if (any(!nzchar(unname(facets)))) {
    stop("`facets` contains an empty column name.", call. = FALSE)
  }

  facets
}

prediction_validate_integer <- function(x,
                                        arg,
                                        min_value = 0L,
                                        positive = FALSE) {
  if (length(x) != 1L) {
    stop("`", arg, "` must be a single integer value.", call. = FALSE)
  }

  x_num <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x_num) || is.na(x_num) || x_num > .Machine$integer.max) {
    if (positive) {
      stop("`", arg, "` must be a positive integer.", call. = FALSE)
    }
    stop("`", arg, "` must be a non-negative integer.", call. = FALSE)
  }

  rounded <- round(x_num)
  if (abs(x_num - rounded) > sqrt(.Machine$double.eps)) {
    if (positive) {
      stop("`", arg, "` must be a positive integer.", call. = FALSE)
    }
    stop("`", arg, "` must be a non-negative integer.", call. = FALSE)
  }

  x_int <- as.integer(rounded)
  if (positive && x_int <= 0L) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }
  if (!positive && x_int < min_value) {
    stop("`", arg, "` must be a non-negative integer.", call. = FALSE)
  }

  x_int
}

prepare_mfrm_prediction_data <- function(fit,
                                         new_data,
                                         person = NULL,
                                         facets = NULL,
                                         score = NULL,
                                         weight = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be output from fit_mfrm().", call. = FALSE)
  }
  if (!is.data.frame(new_data)) {
    stop("`new_data` must be a data.frame.", call. = FALSE)
  }
  if (nrow(new_data) == 0) {
    stop("`new_data` has zero rows.", call. = FALSE)
  }

  source_cols <- fit$config$source_columns %||% fit$prep$source_columns %||% list()
  facet_map <- resolve_prediction_facets(fit, facets = facets)
  person_col <- as.character(person[1] %||% source_cols$person %||% "Person")
  score_col <- as.character(score[1] %||% source_cols$score %||% "Score")
  weight_col <- if (is.null(weight) && !is.null(source_cols$weight)) {
    as.character(source_cols$weight[1])
  } else if (is.null(weight)) {
    NULL
  } else {
    as.character(weight[1])
  }

  required <- c(person_col, unname(facet_map), score_col, weight_col)
  required <- required[!is.na(required) & nzchar(required)]
  if (length(unique(required)) != length(required)) {
    stop("Prediction columns must be distinct. Check `person`, `facets`, `score`, and `weight`.",
         call. = FALSE)
  }
  missing_cols <- setdiff(required, names(new_data))
  if (length(missing_cols) > 0) {
    stop("Prediction columns not found in `new_data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }

  df <- as.data.frame(new_data[, required, drop = FALSE], stringsAsFactors = FALSE)
  names(df) <- c("Person", names(facet_map), "Score", if (!is.null(weight_col)) "Weight")

  blank_to_na <- function(x) {
    x <- as.character(x)
    x[trimws(x) == ""] <- NA_character_
    x
  }

  raw_score <- as.character(df$Score)
  raw_weight <- if ("Weight" %in% names(df)) as.character(df$Weight) else NULL
  score_num <- suppressWarnings(as.numeric(raw_score))
  weight_num <- if (is.null(raw_weight)) {
    rep(1, nrow(df))
  } else {
    suppressWarnings(as.numeric(raw_weight))
  }

  df$Person <- blank_to_na(df$Person)
  for (facet in names(facet_map)) {
    df[[facet]] <- blank_to_na(df[[facet]])
  }

  bad_score <- is.na(score_num) & !is.na(raw_score) & nzchar(trimws(raw_score))
  bad_weight <- if (is.null(raw_weight)) {
    rep(FALSE, nrow(df))
  } else {
    is.na(weight_num) & !is.na(raw_weight) & nzchar(trimws(raw_weight))
  }
  missing_required <- is.na(df$Person) | is.na(score_num)
  for (facet in names(facet_map)) {
    missing_required <- missing_required | is.na(df[[facet]])
  }
  nonpositive_weight <- is.na(weight_num) | weight_num <= 0
  drop_rows <- missing_required | bad_score | nonpositive_weight

  row_review <- tibble::tibble(
    InputRows = nrow(df),
    KeptRows = sum(!drop_rows),
    DroppedRows = sum(drop_rows),
    DroppedMissing = sum(missing_required),
    DroppedBadScore = sum(bad_score),
    DroppedBadWeight = sum(bad_weight),
    DroppedNonpositiveWeight = sum(nonpositive_weight & !is.na(weight_num))
  )

  if (any(drop_rows)) {
    warning(
      "Dropped ", sum(drop_rows), " row(s) from `new_data` before posterior scoring due to missing, non-numeric, or non-positive values.",
      call. = FALSE
    )
  }

  df <- df[!drop_rows, , drop = FALSE]
  score_num <- score_num[!drop_rows]
  weight_num <- weight_num[!drop_rows]
  if (nrow(df) == 0) {
    stop("No valid rows remain in `new_data` after removing missing/invalid observations.",
      call. = FALSE)
  }

  input_data <- df
  input_data$Score <- as.numeric(score_num)
  input_data$Weight <- as.numeric(weight_num)
  input_data <- as.data.frame(input_data, stringsAsFactors = FALSE)

  score_map <- fit$prep$score_map %||% tibble::tibble(
    OriginalScore = seq(fit$prep$rating_min, fit$prep$rating_max),
    InternalScore = seq(fit$prep$rating_min, fit$prep$rating_max)
  )
  internal_score <- score_map$InternalScore[match(score_num, score_map$OriginalScore)]
  unknown_scores <- sort(unique(score_num[is.na(internal_score)]))
  if (length(unknown_scores) > 0) {
    stop(
      "Prediction scores are outside the calibration score support: ",
      paste(unknown_scores, collapse = ", "),
      ". Use the same observed score coding used during model fitting.",
      call. = FALSE
    )
  }

  calibration_levels <- fit$prep$levels[names(facet_map)]
  for (facet in names(facet_map)) {
    unknown_levels <- sort(setdiff(unique(df[[facet]]), calibration_levels[[facet]]))
    if (length(unknown_levels) > 0) {
      stop(
        "Prediction data contain unseen levels for facet `", facet, "`: ",
        paste(unknown_levels, collapse = ", "),
        ". Score future units only against previously calibrated non-person facet levels.",
        call. = FALSE
      )
    }
  }

  pred_df <- df
  pred_df$Person <- factor(pred_df$Person)
  for (facet in names(facet_map)) {
    pred_df[[facet]] <- factor(pred_df[[facet]], levels = calibration_levels[[facet]])
  }
  pred_df$Score <- as.integer(internal_score)
  pred_df$Weight <- as.numeric(weight_num)
  pred_df$score_k <- pred_df$Score - fit$prep$rating_min

  prep <- list(
    data = pred_df,
    n_obs = nrow(pred_df),
    weighted_n = sum(pred_df$Weight, na.rm = TRUE),
    n_person = length(levels(pred_df$Person)),
    rating_min = fit$prep$rating_min,
    rating_max = fit$prep$rating_max,
    score_map = score_map,
    facet_names = names(facet_map),
    levels = c(list(Person = levels(pred_df$Person)), calibration_levels),
    weight_col = if (!is.null(weight_col)) weight_col else NULL,
    keep_original = isTRUE(fit$prep$keep_original),
    source_columns = list(
      person = person_col,
      facets = unname(facet_map),
      score = score_col,
      weight = weight_col
    )
  )

  list(prep = prep, row_review = row_review, input_data = input_data)
}

filter_mfrm_prediction_persons <- function(prepared, keep_persons) {
  keep_persons <- as.character(keep_persons %||% character(0))
  keep_mask <- as.character(prepared$prep$data$Person) %in% keep_persons

  prep_data <- prepared$prep$data[keep_mask, , drop = FALSE]
  prep_data$Person <- droplevels(prep_data$Person)
  input_data <- prepared$input_data[
    as.character(prepared$input_data$Person) %in% keep_persons,
    ,
    drop = FALSE
  ]

  if (nrow(prep_data) == 0) {
    stop("No valid scored persons remain after applying the population-model person-data policy.",
         call. = FALSE)
  }

  updated_prep <- prepared$prep
  updated_prep$data <- prep_data
  updated_prep$n_obs <- nrow(prep_data)
  updated_prep$weighted_n <- sum(prep_data$Weight, na.rm = TRUE)
  updated_prep$n_person <- length(unique(as.character(prep_data$Person)))
  updated_prep$levels$Person <- levels(prep_data$Person)

  prepared$prep <- updated_prep
  prepared$input_data <- input_data
  prepared
}

prepare_mfrm_prediction_population <- function(fit,
                                               prepared,
                                               person_data = NULL,
                                               person_id = NULL,
                                               population_policy = c("error", "omit")) {
  population_policy <- match.arg(population_policy)
  fit_population <- fit$population %||% list()
  target_columns <- as.character(
    fit_population$design_columns %||%
      names(fit_population$coefficients %||% numeric(0))
  )
  covariate_columns <- setdiff(target_columns, "(Intercept)")
  posterior_basis <- as.character(
    fit$config$posterior_basis %||%
      fit_population$posterior_basis %||%
      "legacy_mml"
  )
  active <- isTRUE(fit_population$active) || identical(posterior_basis, "population_model")
  auto_person_data <- FALSE

  if (!active) {
    return(list(
      active = FALSE,
      prepared = prepared,
      scaffold = NULL,
      spec = NULL,
      population_review = NULL,
      input_data = NULL,
      auto_person_data = FALSE
    ))
  }

  if (!is.data.frame(person_data)) {
    if (length(covariate_columns) == 0L) {
      person_data <- data.frame(
        Person = unique(as.character(prepared$input_data$Person)),
        stringsAsFactors = FALSE
      )
      person_id <- "Person"
      auto_person_data <- TRUE
    } else {
      stop(
        "`person_data` must be supplied when scoring a latent-regression fit with background covariates. ",
        "Provide one row per scored person with the background variables used in the fitted `population_formula`.",
        call. = FALSE
      )
    }
  }

  scaffold <- prepare_mfrm_population_scaffold(
    data = prepared$input_data,
    person = "Person",
    population_formula = fit_population$formula %||% fit$config$population_formula,
    person_data = person_data,
    person_id = person_id,
    population_policy = population_policy,
    population_xlevels = fit_population$xlevels %||% fit$config$population_spec$xlevels,
    population_contrasts = fit_population$contrasts %||% fit$config$population_spec$contrasts,
    require_full_rank = FALSE
  )

  current_columns <- colnames(scaffold$design_matrix %||% matrix(numeric(0), nrow = 0L))
  if (length(target_columns) > 0L) {
    if (setequal(current_columns, target_columns)) {
      scaffold$design_matrix <- scaffold$design_matrix[, target_columns, drop = FALSE]
    } else {
      stop(
        "The scored `person_data` do not reproduce the latent-regression design matrix used in the fitted model. ",
        "Check that `person_data`, `person_id`, the fitted `population_formula`, and any categorical predictor levels use the same model-matrix coding.",
        call. = FALSE
      )
    }
  }

  scaffold$coefficients <- fit_population$coefficients
  scaffold$sigma2 <- fit_population$sigma2
  scaffold$converged <- isTRUE(fit_population$converged)
  scaffold$posterior_basis <- "population_model"

  if (length(scaffold$omitted_persons) > 0L) {
    prepared <- filter_mfrm_prediction_persons(prepared, scaffold$included_persons)
  }

  population_review <- tibble::tibble(
    InputPersons = length(unique(as.character(prepared$input_data$Person))) + length(scaffold$omitted_persons),
    RetainedPersons = length(scaffold$included_persons),
    OmittedPersons = length(scaffold$omitted_persons),
    RetainedRows = as.integer(scaffold$response_rows_retained),
    OmittedRows = as.integer(scaffold$response_rows_omitted),
    Policy = as.character(scaffold$policy %||% population_policy),
    PosteriorBasis = "population_model"
  )

  spec <- compact_population_spec(
    scaffold,
    person_levels = prepared$prep$levels$Person
  )

  list(
    active = TRUE,
    prepared = prepared,
    scaffold = scaffold,
    spec = spec,
    population_review = population_review,
    input_data = as.data.frame(person_data, stringsAsFactors = FALSE),
    auto_person_data = auto_person_data
  )
}

posterior_quantile_on_grid <- function(nodes, probs, p) {
  ord <- order(nodes)
  nodes_ord <- nodes[ord]
  probs_ord <- probs[ord]
  cum <- cumsum(probs_ord)
  hit <- which(cum >= p)[1]
  if (is.na(hit)) {
    return(utils::tail(nodes_ord, 1))
  }
  nodes_ord[hit]
}

compute_person_posterior_summary <- function(idx,
                                             config,
                                             params,
                                             quad,
                                             person_labels,
                                             population_spec = NULL,
                                             interval_level = 0.95,
                                             n_draws = 0,
                                             seed = NULL) {
  n <- length(idx$score_k)
  if (n == 0) {
    empty_tbl <- tibble::tibble(
      Person = character(0),
      Estimate = numeric(0),
      SD = numeric(0),
      Lower = numeric(0),
      Upper = numeric(0),
      Observations = integer(0),
      WeightedN = numeric(0)
    )
    return(list(estimates = empty_tbl, draws = tibble::tibble()))
  }

  base_eta <- compute_base_eta(idx, params, config)
  person_int <- idx$person
  n_nodes <- length(quad$nodes)
  score_k <- idx$score_k
  quad_basis <- resolve_person_quadrature_basis(
    quad = quad,
    population_spec = population_spec,
    person_count = length(person_labels)
  )
  person_nodes <- quad_basis$nodes

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    k_cat <- length(step_cum)
    step_cum_row <- matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
    obs_idx <- cbind(seq_len(n), score_k + 1L)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + person_nodes[person_int, q]
      eta_mat <- outer(eta_q, 0:(k_cat - 1))
      log_num <- eta_mat - step_cum_row
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  } else if (identical(config$model, "GPCM")) {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    k_cat <- ncol(step_cum_mat)
    obs_idx <- cbind(seq_len(n), score_k + 1L)
    step_cum_obs <- step_cum_mat[idx$step_idx, , drop = FALSE]
    slope_obs <- matrix(params$slopes[idx$slope_idx], nrow = n, ncol = k_cat)
    k_vals <- 0:(k_cat - 1)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + person_nodes[person_int, q]
      linear_part <- outer(eta_q, k_vals) - step_cum_obs
      log_num <- linear_part * slope_obs
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    k_cat <- ncol(step_cum_mat)
    obs_idx <- cbind(seq_len(n), score_k + 1L)
    step_cum_obs <- step_cum_mat[idx$step_idx, , drop = FALSE]
    k_vals <- 0:(k_cat - 1)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + person_nodes[person_int, q]
      log_num <- outer(eta_q, k_vals) - step_cum_obs
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  }

  ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)
  person_ids <- as.integer(rownames(ll_by_person))
  n_persons <- nrow(ll_by_person)
  aligned_person_labels <- as.character(person_labels[person_ids])
  log_post <- quad_basis$log_weights[person_ids, , drop = FALSE] + ll_by_person
  row_max <- log_post[cbind(seq_len(n_persons), max.col(log_post))]
  log_norm <- row_max + log(rowSums(exp(log_post - row_max)))
  post_w <- exp(log_post - log_norm)

  nodes_mat <- quad_basis$nodes[person_ids, , drop = FALSE]
  eap <- rowSums(nodes_mat * post_w)
  sd_eap <- sqrt(rowSums((nodes_mat - eap)^2 * post_w))
  alpha <- max(min((1 - interval_level) / 2, 0.5), 0)
  lower <- vapply(
    seq_len(n_persons),
    function(i) posterior_quantile_on_grid(nodes_mat[i, ], post_w[i, ], alpha),
    numeric(1)
  )
  upper <- vapply(
    seq_len(n_persons),
    function(i) posterior_quantile_on_grid(nodes_mat[i, ], post_w[i, ], 1 - alpha),
    numeric(1)
  )

  obs_n <- as.integer(rowsum(rep(1L, n), person_int, reorder = FALSE)[, 1])
  weight_n <- as.numeric(rowsum(if (is.null(idx$weight)) rep(1, n) else idx$weight,
                                person_int, reorder = FALSE)[, 1])

  estimates <- tibble::tibble(
    Person = aligned_person_labels,
    Estimate = eap,
    SD = sd_eap,
    Lower = lower,
    Upper = upper,
    Observations = obs_n,
    WeightedN = weight_n
  )

  draws_tbl <- tibble::tibble()
  if (n_draws > 0) {
    draws_tbl <- with_preserved_rng_seed(seed, {
      draw_list <- lapply(seq_len(n_persons), function(i) {
        tibble::tibble(
          Person = aligned_person_labels[i],
          Draw = seq_len(n_draws),
          Value = sample(nodes_mat[i, ], size = n_draws, replace = TRUE, prob = post_w[i, ])
        )
      })
      dplyr::bind_rows(draw_list)
    })
  }

  list(estimates = estimates, draws = draws_tbl)
}

prediction_resolve_fit_method <- function(fit) {
  method <- as.character(
    fit$config$method_input %||%
      fit$config$method %||%
      fit$summary$Method[1] %||%
      NA_character_
  )
  method <- toupper(method[1] %||% NA_character_)
  if (identical(method, "JMLE")) {
    method <- "JML"
  }
  method
}

#' Score future or partially observed units under the fitted scoring basis
#'
#' @param fit Output from [fit_mfrm()] estimated with `method = "MML"` or
#'   `method = "JML"`. When `fit` uses the latent-regression MML branch
#'   (`posterior_basis = "population_model"`), score the target persons with
#'   the same background-variable contract via `person_data`.
#' @param new_data Long-format data for the future or partially observed units
#'   to be scored.
#' @param person Optional person column in `new_data`. Defaults to the person
#'   column recorded in `fit`.
#' @param facets Optional facet-column mapping for `new_data`. Supply either an
#'   unnamed character vector in the calibrated facet order or a named vector
#'   whose names are the calibrated facet names and whose values are the column
#'   names in `new_data`.
#' @param score Optional score column in `new_data`. Defaults to the score
#'   column recorded in `fit`.
#' @param weight Optional weight column in `new_data`. Defaults to the weight
#'   column recorded in `fit`, if any.
#' @param person_data Optional one-row-per-person data.frame with the
#'   background variables required by a latent-regression fit. Ignored for
#'   ordinary fixed-calibration scoring. For intercept-only latent-regression
#'   fits (`population_formula = ~ 1`), `mfrmr` reconstructs the minimal
#'   one-row-per-person table internally from the scored person IDs. This is the
#'   scoring-time table for `new_data`, not the fit object's replay/export
#'   provenance table. For categorical background variables, supply values on
#'   the same coding scale used at fit time; the fitted factor levels and
#'   contrasts are reused when building the scoring design matrix.
#' @param person_id Optional person-ID column in `person_data`. Defaults to
#'   `person` when that column exists, otherwise `"Person"` for the canonical
#'   scoring layout.
#' @param population_policy How missing background data are handled when
#'   `fit` uses the latent-regression branch. `"error"` (default) requires
#'   complete person-level covariates for all scored persons; `"omit"` drops
#'   scored persons lacking complete covariates and records that omission in
#'   `population_review`.
#' @param interval_level Posterior interval level returned in `Lower`/`Upper`.
#' @param n_draws Optional number of quadrature-grid posterior draws to return
#'   per scored person. Use 0 to skip draws.
#' @param seed Optional seed for reproducible posterior draws.
#'
#' @details
#' `predict_mfrm_units()` is the **individual-unit companion** to
#' [predict_mfrm_population()]. It uses the fitted calibration and, when
#' available, the fitted one-dimensional population model to score new or
#' partially observed persons via Expected A Posteriori (EAP) summaries on a
#' quadrature grid.
#'
#' When the original fit uses ordinary `method = "MML"`, the posterior
#' summaries are taken under that fitted MML calibration. When the original fit
#' uses the latent-regression MML branch, the scoring prior is the fitted
#' conditional normal population model \eqn{\theta \mid x \sim
#' N(x^\top\hat\beta, \hat\sigma^2)}, so the returned summaries are
#' population-model-aware posterior EAP estimates. When the original fit uses
#' `method = "JML"`, `mfrmr` applies the fitted facet/step parameters with a
#' standard normal reference prior on the quadrature grid, so the returned
#' person scores remain fixed-calibration EAP summaries rather than direct JML
#' estimates from the fitting step.
#'
#' When the fitted population model is intercept-only (`population_formula = ~
#' 1`), `predict_mfrm_units()` still uses the fitted population-model basis,
#' but it can reconstruct the minimal scored-person table internally because no
#' background covariates are needed beyond the person IDs in `new_data`.
#'
#' The current bounded `GPCM` branch is included in this scoring layer,
#' so fitted `GPCM` objects can be used for the same fixed-calibration
#' posterior summaries. This does not imply that every downstream diagnostic or
#' reporting helper has already been generalized to `GPCM`.
#'
#' This is appropriate for questions such as:
#' - what posterior location/uncertainty do these partially observed new
#'   respondents have under the existing calibration?
#' - how uncertain are those scores, given the observed response pattern?
#'
#' All non-person facet levels in `new_data` must already exist in the fitted
#' calibration. The function does **not** recalibrate the model, update facet
#' estimates, or treat overlapping person IDs as the same latent units from the
#' training data. Person IDs in `new_data` are treated as labels for the rows
#' being scored.
#'
#' When `n_draws > 0`, the returned `draws` component contains discrete
#' quadrature-grid posterior draws that can be used as approximate plausible
#' values under the fitted scoring basis. They should be interpreted as
#' posterior uncertainty summaries, not as deterministic future truth values.
#'
#' For `JML` fits, this scoring stage is intentionally post hoc: `mfrmr` uses
#' the fitted facet and step parameters from the joint-likelihood fit, then
#' adds a standard normal reference prior only for the scoring layer so that
#' new or partially observed units can be summarized on a quadrature grid.
#' This is a practical fixed-calibration EAP procedure, not a claim that the
#' original `JML` fit itself estimated a population model.
#'
#' @section Interpreting output:
#' - `estimates` contains posterior EAP summaries for each person in
#'   `new_data`.
#' - `Lower` and `Upper` are quadrature-grid posterior interval bounds at the
#'   requested `interval_level`.
#' - `SD` is posterior uncertainty under the fitted scoring basis used for
#'   scoring.
#' - `draws`, when requested, contains approximate plausible values on the
#'   fitted quadrature grid.
#' - `population_review`, when present, records whether scored persons were
#'   omitted because their background data were incomplete for a
#'   latent-regression fit.
#'
#' @section What this does not justify:
#' This helper does not update the original calibration, estimate new non-person
#' facet levels, or produce deterministic future person true values. It scores
#' new response patterns under the fitted calibration and, when applicable, the
#' fitted one-dimensional population model.
#'
#' @section References:
#' The posterior summaries follow the usual quadrature-based EAP scoring
#' framework used in item response modeling under calibrated parameters
#' (for example Bock & Aitkin, 1981). When `fit` uses the latent-regression
#' branch, `mfrmr` scores under the fitted conditional normal population model
#' in the general plausible-values spirit discussed by Mislevy (1991). Optional
#' posterior draws are exposed as quadrature-grid plausible-value-style
#' summaries for practical many-facet scoring rather than as a claim of full
#' ConQuest numerical equivalence. When the source fit is `JML`, the same
#' literature supports
#' the quadrature-based scoring layer, but the standard normal prior is a
#' package-level reference prior introduced for post hoc scoring rather than an
#' estimated population distribution.
#'
#' - Bock, R. D., & Aitkin, M. (1981). *Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm*. Psychometrika, 46(4),
#'   443-459.
#' - Mislevy, R. J. (1991). *Randomization-based inference about latent
#'   variables from complex samples*. Psychometrika, 56(2), 177-196.
#' - Muraki, E. (1992). *A generalized partial credit model: Application of an
#'   EM algorithm*. Applied Psychological Measurement, 16(2), 159-176.
#'
#' @return An object of class `mfrm_unit_prediction` with components:
#' - `estimates`: posterior summaries by person
#' - `draws`: optional quadrature-grid posterior draws
#' - `row_review`: row-level preparation review for `new_data`
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `input_data`: cleaned canonical scoring rows retained from `new_data`
#' - `person_data`: cleaned or supplied person-level background data used for
#'   latent-regression scoring; `NULL` otherwise
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_population()], [fit_mfrm()],
#'   [summary.mfrm_unit_prediction]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- suppressWarnings(
#'   fit_mfrm(
#'     toy[toy$Person %in% keep_people, , drop = FALSE],
#'     "Person", c("Rater", "Criterion"), "Score",
#'     method = "MML",
#'     quad_points = 5,
#'     maxit = 30
#'   )
#' )
#' raters <- unique(toy$Rater)[1:2]
#' criteria <- unique(toy$Criterion)[1:2]
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01", "NEW02", "NEW02"),
#'   Rater = c(raters[1], raters[2], raters[1], raters[2]),
#'   Criterion = c(criteria[1], criteria[2], criteria[1], criteria[2]),
#'   Score = c(2, 3, 2, 4)
#' )
#' pred_units <- predict_mfrm_units(toy_fit, new_units, n_draws = 0)
#' summary(pred_units)$estimates[, c("Person", "Estimate", "Lower", "Upper")]
#' @export
predict_mfrm_units <- function(fit,
                               new_data,
                               person = NULL,
                               facets = NULL,
                               score = NULL,
                               weight = NULL,
                               person_data = NULL,
                               person_id = NULL,
                               population_policy = c("error", "omit"),
                               interval_level = 0.95,
                               n_draws = 0,
                               seed = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be output from fit_mfrm().", call. = FALSE)
  }
  fit_method <- prediction_resolve_fit_method(fit)
  if (!fit_method %in% c("MML", "JML")) {
    stop("`predict_mfrm_units()` currently supports only fits estimated with method = 'MML' or 'JML'.",
         call. = FALSE)
  }
  interval_level <- as.numeric(interval_level[1])
  if (!is.finite(interval_level) || interval_level <= 0 || interval_level >= 1) {
    stop("`interval_level` must be a single number in (0, 1).", call. = FALSE)
  }
  n_draws <- prediction_validate_integer(n_draws[1] %||% 0L, "n_draws", min_value = 0L, positive = FALSE)

  prepared <- prepare_mfrm_prediction_data(
    fit = fit,
    new_data = new_data,
    person = person,
    facets = facets,
    score = score,
    weight = weight
  )
  population_ready <- prepare_mfrm_prediction_population(
    fit = fit,
    prepared = prepared,
    person_data = person_data,
    person_id = person_id,
    population_policy = population_policy
  )
  prepared <- population_ready$prepared

  idx <- build_indices(
    prepared$prep,
    step_facet = fit$config$step_facet,
    slope_facet = fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  sizes <- build_param_sizes(fit$config)
  params <- expand_params(fit$opt$par, sizes, fit$config)
  quad_points <- fit$config$estimation_control$quad_points %||% 15L
  quad <- gauss_hermite_normal(as.integer(quad_points))

  scored <- compute_person_posterior_summary(
    idx = idx,
    config = fit$config,
    params = params,
    quad = quad,
    person_labels = prepared$prep$levels$Person,
    population_spec = population_ready$spec,
    interval_level = interval_level,
    n_draws = n_draws,
    seed = seed
  )

  calibration_note <- if (isTRUE(population_ready$active)) {
    "Posterior summaries are computed under the fitted MML calibration together with the fitted conditional normal population model for the scored persons."
  } else if (identical(fit_method, "MML")) {
    "Posterior summaries are computed under the fixed fitted MML calibration."
  } else {
    "Posterior summaries are computed under the fixed fitted JML calibration using a standard normal reference prior on the quadrature grid."
  }

  notes <- c(
    calibration_note,
    "Non-person facets in `new_data` must already exist in the fitted calibration.",
    "Overlapping person IDs are treated as labels in `new_data`; the original fitted person estimates are not updated."
  )
  if (identical(fit_method, "JML")) {
    notes <- c(
      notes,
      "For JML fits, the returned person scores are post hoc fixed-calibration EAP summaries rather than direct JML estimates from the original fit."
    )
  }
  if (isTRUE(population_ready$active)) {
    if (isTRUE(population_ready$auto_person_data)) {
      notes <- c(
        notes,
        "This latent-regression fit uses an intercept-only population model, so `mfrmr` reconstructed the minimal scored-person table from the person IDs in `new_data`."
      )
    } else {
      notes <- c(
        notes,
        "For latent-regression fits with background covariates, supply one-row-per-person background data for the scored units so the fitted population model can be used during posterior scoring."
      )
    }
    if (nrow(population_ready$population_review %||% data.frame()) > 0 &&
        isTRUE(population_ready$population_review$OmittedPersons[1] > 0)) {
      notes <- c(
        notes,
        paste0(
          "Population-model scoring omitted ",
          population_ready$population_review$OmittedPersons[1],
          " person(s) and ",
          population_ready$population_review$OmittedRows[1],
          " response row(s) under `population_policy = '",
          population_ready$population_review$Policy[1],
          "'`."
        )
      )
    }
  }
  if (n_draws > 0) {
    notes <- c(
      notes,
      "The `draws` component contains quadrature-grid posterior draws that can be used as approximate plausible-value summaries."
    )
  }

  structure(
    list(
      estimates = scored$estimates,
      draws = scored$draws,
      row_review = prepared$row_review,
      population_review = population_ready$population_review,
      input_data = prepared$input_data,
      person_data = population_ready$input_data,
      settings = list(
        interval_level = interval_level,
        n_draws = n_draws,
        quad_points = as.integer(quad_points),
        seed = seed,
        method = fit_method,
        source_columns = prepared$prep$source_columns,
        posterior_basis = if (isTRUE(population_ready$active)) {
          "population_model"
        } else {
          as.character(fit$config$posterior_basis %||% "legacy_mml")
        },
        person_id = if (isTRUE(population_ready$active)) {
          as.character(population_ready$scaffold$person_id %||% person_id %||% "Person")
        } else {
          NULL
        },
        population_policy = if (isTRUE(population_ready$active)) {
          as.character(population_ready$scaffold$policy %||% population_policy)
        } else {
          NULL
        },
        population_formula = if (isTRUE(population_ready$active)) {
          paste(deparse(population_ready$scaffold$formula), collapse = " ")
        } else {
          NULL
        }
      ),
      notes = notes
    ),
    class = "mfrm_unit_prediction"
  )
}

#' Summarize posterior unit scoring output
#'
#' @param object Output from [predict_mfrm_units()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_unit_prediction` with:
#' - `estimates`: posterior summaries by person
#' - `row_review`: row-preparation review
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_units()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- suppressWarnings(
#'   fit_mfrm(
#'     toy[toy$Person %in% keep_people, , drop = FALSE],
#'     "Person", c("Rater", "Criterion"), "Score",
#'     method = "MML",
#'     quad_points = 5,
#'     maxit = 30
#'   )
#' )
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01"),
#'   Rater = unique(toy$Rater)[1],
#'   Criterion = unique(toy$Criterion)[1:2],
#'   Score = c(2, 3)
#' )
#' pred_units <- predict_mfrm_units(toy_fit, new_units)
#' summary(pred_units)
#' @method summary mfrm_unit_prediction
#' @export
summary.mfrm_unit_prediction <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_unit_prediction")) {
    stop("`object` must be output from predict_mfrm_units().", call. = FALSE)
  }
  digits <- prediction_validate_integer(digits, "digits", min_value = 0L, positive = FALSE)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out <- list(
    estimates = round_df(object$estimates),
    row_review = round_df(object$row_review),
    population_review = round_df(object$population_review),
    settings = object$settings,
    notes = object$notes %||% character(0),
    digits = digits
  )
  class(out) <- "summary.mfrm_unit_prediction"
  out
}

#' @export
print.summary.mfrm_unit_prediction <- function(x, ...) {
  digits <- prediction_validate_integer(x$digits %||% 3L, "digits", min_value = 0L, positive = FALSE)
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

  cat("mfrmr Unit Prediction Summary\n")
  if (!is.null(x$estimates) && nrow(x$estimates) > 0) {
    cat("\nPosterior estimates\n")
    print(round_df(as.data.frame(preview_df(x$estimates))), row.names = FALSE)
  }
  if (!is.null(x$row_review) && nrow(x$row_review) > 0) {
    cat("\nRow preparation review\n")
    print(round_df(as.data.frame(x$row_review)), row.names = FALSE)
  }
  if (!is.null(x$population_review) && nrow(x$population_review) > 0) {
    cat("\nPopulation-model omission review\n")
    print(round_df(as.data.frame(x$population_review)), row.names = FALSE)
  }
  if (!is.null(x$settings) && length(x$settings) > 0L) {
    cat("\nSettings\n")
    print(bundle_settings_table(x$settings), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Sample approximate plausible values under fitted posterior scoring
#'
#' @param fit Output from [fit_mfrm()] estimated with `method = "MML"` or
#'   `method = "JML"`.
#' @param new_data Long-format data for the future or partially observed units
#'   to be scored.
#' @param person Optional person column in `new_data`. Defaults to the person
#'   column recorded in `fit`.
#' @param facets Optional facet-column mapping for `new_data`. Supply either an
#'   unnamed character vector in the calibrated facet order or a named vector
#'   whose names are the calibrated facet names and whose values are the column
#'   names in `new_data`.
#' @param score Optional score column in `new_data`. Defaults to the score
#'   column recorded in `fit`.
#' @param weight Optional weight column in `new_data`. Defaults to the weight
#'   column recorded in `fit`, if any.
#' @param person_data Optional one-row-per-person data.frame with the
#'   background variables required by a latent-regression fit. Ignored for
#'   ordinary fixed-calibration scoring. Intercept-only latent-regression fits
#'   can reconstruct the minimal scored-person table internally. This is the
#'   scoring-time table for `new_data`, not the fit object's replay/export
#'   provenance table. For categorical background variables, supply values on
#'   the same coding scale used at fit time; the fitted factor levels and
#'   contrasts are reused when building the scoring design matrix.
#' @param person_id Optional person-ID column in `person_data`.
#' @param population_policy How missing background data are handled when
#'   `fit` uses the latent-regression branch. `"error"` (default) requires
#'   complete person-level covariates; `"omit"` drops scored persons lacking
#'   complete covariates and records that omission in `population_review`.
#' @param n_draws Number of posterior draws per person. Must be a positive
#'   integer.
#' @param interval_level Posterior interval level passed to
#'   [predict_mfrm_units()] for the accompanying EAP summary table.
#' @param seed Optional seed for reproducible posterior draws.
#'
#' @details
#' `sample_mfrm_plausible_values()` is a thin public wrapper around
#' [predict_mfrm_units()] that exposes the fixed-calibration posterior draws as
#' a standalone object. It is useful when downstream workflows want repeated
#' latent-value imputations rather than just one posterior EAP summary.
#'
#' In the current `mfrmr` implementation these are **approximate plausible
#' values** drawn from the fitted quadrature-grid posterior under the scoring
#' basis implied by `fit`. For ordinary `MML` fits this is the fitted marginal
#' calibration; for latent-regression `MML` fits it is the fitted conditional
#' normal population model for the scored persons; for `JML` fits it is the
#' fixed facet/step calibration together with a standard normal reference prior
#' on the quadrature grid. They should be interpreted as posterior uncertainty
#' summaries for the scored persons, not as deterministic future truth values
#' and not as a claim of full many-facet plausible-values equivalence with
#' population-model software.
#'
#' In other words, the `JML` path here is a practical scoring approximation
#' layered on top of the fitted joint-likelihood calibration, whereas the
#' latent-regression `MML` path uses the fitted one-dimensional conditional
#' normal population model. Neither path should be described as a full
#' many-facet plausible-values system with all ConQuest-style extensions.
#'
#' @section Interpreting output:
#' - `values` contains one row per person per draw.
#' - `estimates` contains the companion posterior EAP summaries from
#'   [predict_mfrm_units()].
#' - `summary()` reports draw counts and empirical draw summaries by person.
#'
#' @section What this does not justify:
#' This helper does not update the calibration, estimate new non-person facet
#' levels, or provide exact future true values. It samples from the fixed-grid
#' posterior implied by the existing fixed calibration.
#'
#' @section References:
#' The underlying posterior scoring follows the usual quadrature-based EAP
#' framework of Bock and Aitkin (1981). The interpretation of multiple
#' posterior draws as plausible-value-style summaries follows the general logic
#' discussed by Mislevy (1991), while the current implementation remains a
#' practical fixed-calibration approximation rather than a full published
#' many-facet plausible-values method. For `JML` source fits, the quadrature
#' posterior uses a package-level standard normal reference prior for this
#' post hoc scoring layer.
#'
#' - Bock, R. D., & Aitkin, M. (1981). *Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm*. Psychometrika, 46(4),
#'   443-459.
#' - Mislevy, R. J. (1991). *Randomization-based inference about latent
#'   variables from complex samples*. Psychometrika, 56(2), 177-196.
#'
#' @return An object of class `mfrm_plausible_values` with components:
#' - `values`: one row per person per draw
#' - `estimates`: companion posterior EAP summaries
#' - `row_review`: row-preparation review
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `input_data`: cleaned canonical scoring rows retained from `new_data`
#' - `person_data`: cleaned or supplied person-level background data used for
#'   latent-regression scoring; `NULL` otherwise
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_units()], [summary.mfrm_plausible_values]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- suppressWarnings(
#'   fit_mfrm(
#'     toy[toy$Person %in% keep_people, , drop = FALSE],
#'     "Person", c("Rater", "Criterion"), "Score",
#'     method = "MML",
#'     quad_points = 5,
#'     maxit = 30
#'   )
#' )
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01"),
#'   Rater = unique(toy$Rater)[1],
#'   Criterion = unique(toy$Criterion)[1:2],
#'   Score = c(2, 3)
#' )
#' pv <- sample_mfrm_plausible_values(toy_fit, new_units, n_draws = 3, seed = 1)
#' summary(pv)$draw_summary
#' @export
sample_mfrm_plausible_values <- function(fit,
                                         new_data,
                                         person = NULL,
                                         facets = NULL,
                                         score = NULL,
                                         weight = NULL,
                                         person_data = NULL,
                                         person_id = NULL,
                                         population_policy = c("error", "omit"),
                                         n_draws = 5,
                                         interval_level = 0.95,
                                         seed = NULL) {
  n_draws <- prediction_validate_integer(n_draws[1] %||% 0L, "n_draws", positive = TRUE)

  pred <- predict_mfrm_units(
    fit = fit,
    new_data = new_data,
    person = person,
    facets = facets,
    score = score,
    weight = weight,
    person_data = person_data,
    person_id = person_id,
    population_policy = population_policy,
    interval_level = interval_level,
    n_draws = n_draws,
    seed = seed
  )

  fit_method <- prediction_resolve_fit_method(fit)
  draw_note <- if (identical(pred$settings$posterior_basis %||% "", "population_model")) {
    "These draws are sampled from the fitted population-model posterior implied by the latent-regression MML branch."
  } else if (identical(fit_method, "MML")) {
    "These draws are sampled from the fixed quadrature-grid posterior under the existing MML calibration."
  } else {
    "These draws are sampled from a fixed-calibration quadrature-grid posterior built from the fitted JML parameters and a standard normal reference prior."
  }

  notes <- c(
    draw_note,
    "Use them as approximate plausible-value summaries for posterior uncertainty, not as deterministic future truth values."
  )

  structure(
    list(
      values = pred$draws,
      estimates = pred$estimates,
      row_review = pred$row_review,
      population_review = pred$population_review,
      input_data = pred$input_data,
      person_data = pred$person_data,
      settings = pred$settings,
      notes = notes
    ),
    class = "mfrm_plausible_values"
  )
}

#' Summarize approximate plausible values from posterior scoring
#'
#' @param object Output from [sample_mfrm_plausible_values()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_plausible_values` with:
#' - `draw_summary`: empirical summaries of the sampled values by person
#' - `estimates`: companion posterior EAP summaries
#' - `row_review`: row-preparation review
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [sample_mfrm_plausible_values()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- suppressWarnings(
#'   fit_mfrm(
#'     toy[toy$Person %in% keep_people, , drop = FALSE],
#'     "Person", c("Rater", "Criterion"), "Score",
#'     method = "MML",
#'     quad_points = 5,
#'     maxit = 30
#'   )
#' )
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01"),
#'   Rater = unique(toy$Rater)[1],
#'   Criterion = unique(toy$Criterion)[1:2],
#'   Score = c(2, 3)
#' )
#' pv <- sample_mfrm_plausible_values(toy_fit, new_units, n_draws = 3, seed = 1)
#' summary(pv)
#' @method summary mfrm_plausible_values
#' @export
summary.mfrm_plausible_values <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_plausible_values")) {
    stop("`object` must be output from sample_mfrm_plausible_values().", call. = FALSE)
  }
  digits <- prediction_validate_integer(digits, "digits", min_value = 0L, positive = FALSE)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  draw_summary <- object$values |>
    dplyr::group_by(.data$Person) |>
    dplyr::summarise(
      Draws = dplyr::n(),
      MeanValue = mean(.data$Value),
      SDValue = stats::sd(.data$Value),
      LowerValue = stats::quantile(.data$Value, probs = 0.025, names = FALSE, type = 1),
      UpperValue = stats::quantile(.data$Value, probs = 0.975, names = FALSE, type = 1),
      .groups = "drop"
    )

  out <- list(
    draw_summary = round_df(draw_summary),
    estimates = round_df(object$estimates),
    row_review = round_df(object$row_review),
    population_review = round_df(object$population_review),
    settings = object$settings,
    notes = object$notes %||% character(0),
    digits = digits
  )
  class(out) <- "summary.mfrm_plausible_values"
  out
}

#' @export
print.summary.mfrm_plausible_values <- function(x, ...) {
  digits <- prediction_validate_integer(x$digits %||% 3L, "digits", min_value = 0L, positive = FALSE)
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

  cat("mfrmr Plausible Values Summary\n")
  if (!is.null(x$draw_summary) && nrow(x$draw_summary) > 0) {
    cat("\nDraw summary\n")
    print(round_df(as.data.frame(preview_df(x$draw_summary))), row.names = FALSE)
  }
  if (!is.null(x$estimates) && nrow(x$estimates) > 0) {
    cat("\nCompanion estimates\n")
    print(round_df(as.data.frame(preview_df(x$estimates))), row.names = FALSE)
  }
  if (!is.null(x$row_review) && nrow(x$row_review) > 0) {
    cat("\nRow preparation review\n")
    print(round_df(as.data.frame(x$row_review)), row.names = FALSE)
  }
  if (!is.null(x$population_review) && nrow(x$population_review) > 0) {
    cat("\nPopulation-model omission review\n")
    print(round_df(as.data.frame(x$population_review)), row.names = FALSE)
  }
  if (!is.null(x$settings) && length(x$settings) > 0L) {
    cat("\nSettings\n")
    print(bundle_settings_table(x$settings), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}
