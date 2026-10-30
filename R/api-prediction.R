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
#'   role keywords (`person`, `rater`, `criterion`, `assignment`). The nested
#'   named-facet form `design$facets = c(person = ..., judge = ..., task = ...)`
#'   is also accepted for the supported person/rater/criterion design.
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
#' Generating parameters are held fixed at the values in that specification.
#' Monte Carlo SEs describe simulation precision conditional on these values;
#' they do not include source-fit uncertainty or uncertainty about whether the
#' generating assumptions describe a future administration.
#'
#' When that fit-derived or manually built simulation specification stores an
#' active latent-regression population generator, the helper still operates at
#' the **design / operating-characteristic** level. It repeatedly simulates
#' person-level covariates and responses, refits the MML population-model
#' branch, and summarizes the resulting facet-level behavior. This is distinct
#' from the fitted-model posterior scoring provided by [predict_mfrm_units()].
#'
#' `GPCM` forecasts are available with caveats through the same
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
#' - `fit_identity`: source-fit metadata for fit-derived forecasts; `NULL` for
#'   an explicitly supplied simulation specification
#' - `sim_spec`: simulation specification used for the forecast
#' - `facet_names`: public non-person facet names carried by the simulation
#'   specification
#' - `design_variable_aliases`: public aliases for
#'   `n_person`/`n_rater`/`n_criterion`/`raters_per_person`
#' - `design_descriptor`: role-based description of design variables carried
#'   from the underlying planning object
#' - `planning_scope`: explicit record of the supported planning contract,
#'   including a `facet_manifest`
#' - `planning_constraints`: explicit record of mutable/locked design variables
#' - `planning_schema`: structured planning metadata carrying the role
#'   table, supported boundary, mutability map, and facet manifest
#' - `gpcm_boundary`: `GPCM` caveat row when a `GPCM` forecast route is
#'   used
#' - `settings`: forecasting settings
#' - `ademp`: simulation-study metadata
#' - `notes`: interpretation notes
#' @seealso [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()],
#'   [evaluate_mfrm_design()], [summary.mfrm_population_prediction]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("predict_mfrm_population", run.dontrun = TRUE)
#' \dontrun{
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
#' @section Saved-result reuse:
#' Fit-derived forecasts must match the source fit when supplied to archive or
#' replay helpers. Explicit `sim_spec` forecasts remain separate scenarios;
#' attaching them does not establish that their assumptions match the fit.
#' Fit-derived forecasts with a fixed `population_sd` are unavailable: the
#' design-simulation refits do not preserve that variance constraint or an
#' anchor-defined origin. An explicit `sim_spec` can describe a separate
#' forecasting scenario, but does not reproduce those source-fit constraints.
#' See [mfrmr_compatibility_layer] for migration and replay behavior.
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
    if (!is.null(fit$population$fixed_sd %||% fit$config$population_spec$fixed_sd)) {
      stop("Fit-derived population forecasts do not yet preserve fixed `population_sd` or an anchor-defined origin during simulation refits. Supply an explicit `sim_spec` only for a separately declared forecasting scenario; it does not reproduce those source-fit constraints.", call. = FALSE)
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
    if (has_fit) "This forecast uses a simulation specification derived from the recorded source fit." else
      "This forecast uses an explicitly supplied simulation specification; it is not asserted to derive from a fit stored alongside it.",
    "MCSE columns quantify Monte Carlo uncertainty from using a finite number of replications.",
    "Generating parameters are held fixed; MCSEs do not include source-fit uncertainty or uncertainty about the generating assumptions.",
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
      fit_identity = if (has_fit) prediction_fit_identity(fit) else NULL,
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
#' - `planning_schema`: structured planning metadata
#' - `gpcm_boundary`: `GPCM` caveat row when present
#' - `structural_design_review`: deterministic structural review of the
#'   named-facet design grid; it is not
#'   a forecast-uncertainty result
#' - `ademp`: simulation-study metadata
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_population()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("summary.mfrm_population_prediction", run.dontrun = TRUE)
#' \dontrun{
#' # One saved replication illustrates output structure, not precise rates.
#' # Regeneration recipe: examples/simulation-results.R in the package.
#' saved <- readRDS(system.file("examples", "simulation-population.rds",
#'   package = "mfrmr", mustWork = TRUE))
#' # Inspect saved$generation_warnings and the result's run/status tables.
#' pred <- saved$result
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
    structural_design_review = simulation_compact_structural_design_review_summary(
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
  if (inherits(out$structural_design_review, "summary.mfrm_structural_design_review")) {
    out$notes <- c(
      out$notes,
      "The structural design review reports deterministic bookkeeping and conservative design guidance, not forecast uncertainty."
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
  print_compact_structural_design_review_summary(
    x$structural_design_review %||% NULL,
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

  # Two-family parameter order follows slope owners, which may differ from
  # the facets argument in the original call. Their native names are identities.
  if (is.null(facets) && mfrm_has_product_slopes(fit)) facets <- setNames(fit_facets, fit_facets)
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

  if (!is.null(raw_weight)) {
    invalid_weight <- is.na(weight_num) | !is.finite(weight_num) |
      weight_num <= 0
    if (any(invalid_weight)) {
      stop(
        "Scoring weights must be finite and strictly positive. Invalid ",
        "value(s) were found at row(s): ",
        paste(utils::head(which(invalid_weight), 10L), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
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

# Scoring also depends on future-population coding and source admission states.
prediction_fit_identity <- function(fit) {
  pop <- fit$population %||% list()
  spec <- fit$config$population_spec %||% list()
  active <- isTRUE(pop$active) || identical(fit$config$posterior_basis, "population_model")
  terms <- if (active) prediction_population_terms(fit) else NULL
  # The fitted transformation basis belongs to the identity, its workspace does not.
  if (!is.null(terms)) attr(terms, ".Environment") <- NULL
  list(calibration = mfrm_result_fit_identity(fit), scoring = mfrm_checkpoint_fingerprint(list(
    posterior_basis = fit$config$posterior_basis,
    terms = terms, xlevels = pop$xlevels %||% spec$xlevels,
    contrasts = pop$contrasts %||% spec$contrasts,
    design_columns = pop$design_columns %||% names(pop$coefficients),
    source_columns = fit$config$source_columns,
    readiness = mfrmr_get_readiness_record(fit),
    numerical = fit$opt$optimizer_diagnostics,
    repeated_ratings = fit$config$repeated_ratings,
    preparation_notes = fit$prep$preparation_notes %||% fit$data_review$preparation_notes,
    estimability = fit$config$estimability_audit,
    data_estimability = fit$data_review$estimability,
    boundary = fit$config$boundary_audit,
    adjustment = fit$jml_adjustment[c("estimator", "point", "policy")]
  ), canonical = TRUE))
}

prediction_population_terms <- function(fit) {
  pop <- fit$population %||% list()
  spec <- fit$config$population_spec %||% list()
  stored <- pop$terms %||% spec$terms
  if (!is.null(stored)) return(stored)

  formula <- pop$formula %||% fit$config$population_formula
  terms <- stats::terms(formula)
  # Plain variables need only the retained factor levels/contrasts. Calls may
  # have learned a basis (scale, poly, splines); never relearn it on new Persons.
  variables <- as.list(attr(terms, "variables"))[-1L]
  if (!any(vapply(variables, is.call, logical(1)))) return(terms)

  training <- pop$person_table_replay %||% pop$person_table
  person_id <- pop$person_id %||% spec$person_id
  if (!is.data.frame(training) || length(person_id) != 1L ||
      !person_id %in% names(training)) {
    stop("Transformed population predictors require fitted terms or retained training person data. Refit the calibration before scoring.", call. = FALSE)
  }
  rebuilt <- prepare_mfrm_population_scaffold(
    data = training, person = person_id, population_formula = formula,
    person_data = training, person_id = person_id,
    population_policy = pop$policy %||% "error",
    population_xlevels = pop$xlevels %||% spec$xlevels,
    population_contrasts = pop$contrasts %||% spec$contrasts,
    require_full_rank = FALSE
  )
  design <- pop$design_matrix %||% spec$design_matrix
  if (is.null(design) || !identical(colnames(rebuilt$design_matrix), colnames(design)) ||
      !isTRUE(all.equal(unname(rebuilt$design_matrix), unname(design), tolerance = 1e-12))) {
    stop("Retained person data do not reproduce the fitted population training design. Refit the calibration before scoring.", call. = FALSE)
  }
  rebuilt$terms
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
    require_full_rank = FALSE,
    population_terms = prediction_population_terms(fit)
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
  fixed_sd <- fit_population$fixed_sd %||% fit$config$population_spec$fixed_sd
  if (!is.null(fixed_sd)) scaffold$fixed_sd <- fixed_sd
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

compute_person_posterior_summary <- function(idx,
                                             config,
                                             params,
                                             quad,
                                             person_labels,
                                             population_spec = NULL,
                                             interval_level = 0.95,
                                             n_draws = 0,
                                             seed = NULL, cut_scores = NULL) {
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
    out <- list(estimates = empty_tbl, draws = tibble::tibble())
    if (!is.null(cut_scores)) out$decision_probabilities <- matrix(numeric(),0L,2L*length(cut_scores)+1L)
    return(out)
  }

  base_eta <- compute_base_eta(idx, params, config)
  person_int <- idx$person
  quad_basis <- if (mfrmr_adaptive_integration(config)) {
    mfrmr_adaptive_quadrature_basis(idx, config, params, quad, base_eta,
      population_spec = population_spec, person_count = length(person_labels))
  } else resolve_person_quadrature_basis(
    quad = quad, population_spec = population_spec, person_count = length(person_labels)
  )
  # Use the fitted-model kernels with this scoring prior/person grid.
  quad$person_basis <- quad_basis
  logprob <- mfrm_mml_logprob_bundle_r(idx, config, quad, params, base_eta)
  posterior <- mfrm_mml_person_bundle(logprob$log_prob_mat, person_int,
    quad_basis, include_posterior = TRUE)
  person_ids <- posterior$person_ids
  n_persons <- length(person_ids)
  aligned_person_labels <- as.character(person_labels[person_ids])
  post_w <- posterior$posterior

  nodes_mat <- quad_basis$nodes[person_ids, , drop = FALSE]
  eap <- rowSums(nodes_mat * post_w)
  sd_eap <- sqrt(rowSums((nodes_mat - eap)^2 * post_w))
  intervals <- mfrmr_person_posterior_intervals(idx, config, params, base_eta,
    person_labels, population_spec, interval_level, cut_scores)

  obs_n <- as.integer(rowsum(rep(1L, n), person_int, reorder = FALSE)[, 1])
  weight_n <- as.numeric(rowsum(if (is.null(idx$weight)) rep(1, n) else idx$weight,
                                person_int, reorder = FALSE)[, 1])

  estimates <- tibble::tibble(
    Person = aligned_person_labels,
    Estimate = eap,
    SD = sd_eap,
    Lower = intervals[, "Lower"],
    Upper = intervals[, "Upper"],
    PriorMean = if (isTRUE(quad_basis$transformed)) quad_basis$mu[person_ids] else 0,
    PriorSD = if (isTRUE(quad_basis$transformed)) quad_basis$sigma else 1,
    WeightedLikelihood = as.numeric(rowsum(as.integer((idx$weight %||% rep(1, n)) != 1), person_int,
                                           reorder = FALSE)[, 1]) > 0,
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

  out <- list(estimates = estimates, draws = draws_tbl)
  if (!is.null(cut_scores)) out$decision_probabilities <- intervals[, -(1:2), drop=FALSE]
  out
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

# Compare the actual reported EAP/SD with two adaptive reference orders.
prediction_score_integration_review <- function(estimates, review) {
  rows <- lapply(seq_len(nrow(estimates)), function(i) {
    d <- review[review$Person == estimates$Person[i], , drop = FALSE]
    d <- d[order(d$AdaptiveNodes), , drop = FALSE]
    out <- data.frame(Person = estimates$Person[i], EAPChange = NA_real_,
      SDChange = NA_real_, ReferenceChange = NA_real_, Passed = FALSE)
    if (nrow(d) < 2L || any(d$Status != "computed")) return(out)
    last <- tail(d, 2L)
    out$EAPChange <- abs(estimates$Estimate[i] - last$AdaptiveEAP[2])
    out$SDChange <- abs(estimates$SD[i] - last$AdaptivePosteriorSD[2])
    out$ReferenceChange <- max(abs(diff(last$AdaptiveEAP)),
      abs(diff(last$AdaptivePosteriorSD)), abs(diff(last$AdaptiveLogMarginal)))
    out$Passed <- all(is.finite(unlist(out[-1]))) &&
      max(out$EAPChange, out$SDChange, out$ReferenceChange) <= 1e-5
    out
  })
  do.call(rbind, rows)
}

# Output-specific local calibration checks; never changes global readiness.
prediction_local_calibration_review <- function(fit, covariance = NULL) {
  out <- list(eligible = FALSE, basis = "local_mml_and_integration_v1",
    review = "Local calibration checks were not completed.")
  tryCatch({
    if (isTRUE(fit$config$population_spec$active) &&
        !identical(fit$config$population_spec$design_columns, "(Intercept)"))
      mfrm_native_location_source(fit)
    # The promise is evaluated only after the existing local-source guards.
    check <- mfrm_ic_fit_check(fit, mfrm_extract_fit_ic_contract(fit),
      information = if (is.null(covariance)) NULL else tryCatch({
        info <- covariance(fit)
        info$solution_information %||% list(status = "unavailable", detail = info$detail)
      }, error = function(e) list(status = "unavailable", detail = conditionMessage(e))))
    out$review <- check$review
    out$caution <- check$caution %||% character(0)
    if (!isTRUE(check$eligible)) return(out)
    cfg <- fit$config
    order <- as.integer(cfg$estimation_control$quad_points)
    higher <- max(order + 10L, 2L * order - 1L)
    if (inherits(tryCatch(gauss_hermite_normal(higher), error = identity), "error")) {
      higher <- order + 10L
    }
    sizes <- build_param_sizes(cfg)
    idx <- build_indices(fit$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
    evaluate <- function(q) {
      cfg$estimation_control$quad_points <- q
      evaluator <- make_mfrm_direct_evaluator("MML",
        make_param_cache(sizes, cfg, idx, is_mml = TRUE), idx, cfg, sizes,
        gauss_hermite_normal(q))
      list(value = evaluator$value(fit$opt$par), gradient = evaluator$gradient(fit$opt$par))
    }
    original <- evaluate(order); refined <- evaluate(higher)
    out$integration <- data.frame(OriginalOrder = order, ComparisonOrder = higher,
      NLLChange = abs(original$value - refined$value),
      GradientChange = max(abs(original$gradient - refined$gradient)),
      ComparisonGradient = max(abs(refined$gradient)))
    out$eligible <- all(is.finite(unlist(out$integration))) &&
      out$integration$NLLChange <= 1e-5 && out$integration$GradientChange <= 1e-4 &&
      out$integration$ComparisonGradient <= 1e-4
    if (!out$eligible) out$review <- paste(
      "Calibration integration changes exceed the scoring review tolerance.",
      "Refit with more quadrature points and inspect mml_quadrature_sensitivity().")
    out
  }, error = function(e) {
    out$review <- paste("Local calibration checks are unavailable:", conditionMessage(e))
    out
  })
}

# Local JML checks retain the original global identification/boundary states.
prediction_jml_calibration_review <- function(fit) {
  out <- list(eligible = FALSE, basis = "local_jml_and_curvature_v1",
    review = "Local joint-likelihood checks were not completed.")
  tryCatch({
    cfg <- fit$config
    if (!identical(cfg$model, "GPCM") || !identical(cfg$method, "JML") ||
        !identical(cfg$slope_facet, cfg$step_facet) ||
        length(cfg$interaction_specs) || nrow(mfrmr_calibration_extract_anchors(fit)) ||
        any(!is.na(cfg$theta_spec$anchors)) || length(cfg$theta_spec$group_values) ||
        isTRUE(cfg$population_spec$active) || anyNA(fit$prep$data$Weight) ||
        any(fit$prep$data$Weight != 1))
      stop("Conditional GPCM JML scoring requires shared owners, unit weights, no anchors/interactions and no fitted population model.")
    audit <- cfg$boundary_audit
    person_status <- audit$parameter_status
    if (!is.data.frame(person_status) || !nrow(person_status) ||
        !"ParameterStatus" %in% names(person_status))
      stop("The source must retain its Person boundary audit.")
    slope_status <- mfrmr_readiness_gpcm_slope_parameters(cfg)$ParameterStatus
    states <- vapply(audit[c("structural_additive", "joint_additive",
      "gpcm_slope_boundary", "gpcm_joint_boundary")],
      function(x) as.character(x$state %||% "missing"), character(1))
    if (any(states == "missing")) stop("The source must retain its additive and slope boundary audits.")
    if (any(grepl("^unbounded", c(person_status$ParameterStatus, slope_status))) ||
        any(startsWith(states, "certified_")))
      stop("A retained Person, additive or slope boundary certificate prevents automatic conditional scoring; finite optimizer traces do not replace it.")
    sizes <- build_param_sizes(cfg)
    idx <- build_indices(fit$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
    evaluator <- make_mfrm_direct_evaluator("JML",
      make_param_cache(sizes, cfg, idx, is_mml = FALSE), idx, cfg, sizes)
    par <- fit$opt$par
    objective_difference <- abs(evaluator$value(par) - fit$opt$value)
    gradient <- evaluator$gradient(par)
    if (length(objective_difference) != 1L || !is.finite(objective_difference) ||
        objective_difference > 1e-6 || !length(gradient) || any(!is.finite(gradient)) ||
        max(abs(gradient)) > 1e-4)
      stop("The joint objective must agree within 1e-6 and its maximum absolute gradient must be <= 1e-4.")
    # Includes all free Person and structural parameters; no regularization or inversion.
    information <- optimHess(par, evaluator$value, evaluator$gradient)
    if (!identical(dim(information), rep(length(par), 2L)) || any(!is.finite(information)))
      stop("Joint JML curvature is unavailable.")
    chol(information)
    condition <- rcond(information)
    if (!is.finite(condition) || condition <= 1e-10)
      stop("Joint JML curvature is numerically singular.")
    out$checks <- list(objective_difference = objective_difference,
      max_abs_gradient = max(abs(gradient)), dimension = as.integer(length(par)),
      positive_definite = TRUE, reciprocal_condition = condition,
      known_boundary_absent = TRUE)
    out$eligible <- TRUE
    out$review <- "Fresh joint-likelihood and unregularized full-curvature checks passed for conditional scoring."
    out$caution <- "Local JML checks do not establish a unique global maximum, complete the global boundary audit or qualify formal slope intervals."
    out
  }, error = function(e) {
    out$review <- paste("Local JML calibration checks are unavailable:", conditionMessage(e))
    out
  })
}

# Fixed-calibration scoring needs the saved point solution, not a covariance
# certificate. Reuse its identity check, then compare the fitted integration
# with a finer rule without fitting or changing the source readiness record.
prediction_gmfrm_source_identity <- function(fit, legacy = FALSE) {
  # An exact saved-source match, not a new numerical check or authenticity claim.
  paste0(if (legacy) "gmfrm_scoring_source_v1:" else "gmfrm_scoring_source_v2:",
    mfrm_checkpoint_fingerprint(list(
    source = mfrm_gpcm_inference_source(fit), config = fit$config,
    specification = fit$gmfrm, score_map = fit$prep$score_map,
    source_columns = fit$prep$source_columns,
    locations = fit$facets$others, steps = fit$steps, slopes = fit$slopes,
    converged = fit$summary$Converged, numerical = fit$opt$optimizer_diagnostics,
    readiness = mfrmr_get_readiness_record(fit)), canonical = !legacy))
}

prediction_gmfrm_calibration <- function(fit, params) {
  list(slope_facets = fit$config$slope_facet, step_facet = fit$config$step_facet,
    scale_reference = "fixed_standard_normal", levels = fit$prep$levels[fit$config$slope_facet],
    score_map = fit$prep$score_map, locations = fit$facets$others,
    steps = fit$steps, slope_components = params$slope_components,
    observed_contexts = unique(fit$prep$data[fit$config$slope_facet]))
}

prediction_gmfrm_scoring_readiness <- function(fit) {
  input <- mfrm_gmfrm_response_input(fit, operation = "new-Person scoring")
  cfg <- input$config
  record <- mfrmr_get_readiness_record(fit)$fit
  field <- function(name) if (is.data.frame(record) && nrow(record) == 1L &&
      name %in% names(record)) as.character(record[[name]]) else "missing"
  eligible <- identical(field("ReadinessContractVersion"), mfrmr_readiness_contract_version()) &&
    field("InputState") == "pass" && field("CategoryState") == "adequate" &&
    field("NumericalState") == "ready" &&
    field("EstimabilityState") %in% c("identified", "not_evaluated") &&
    field("BoundaryState") %in% c("finite", "not_applicable", "not_evaluated") &&
    !isTRUE(cfg$estimability_audit$nonlinear_local_estimability$local_first_order_rank_deficient)
  local <- list(eligible = FALSE, basis = "two_family_point_calibration_v1",
    review = "The source readiness record requires review.",
    caution = "Two-family scoring is experimental. Numerical checks do not establish a global maximum, complete identification or boundary audits, or validate future-population scores.")
  if (isTRUE(eligible)) local <- tryCatch({
    q <- as.integer(cfg$estimation_control$quad_points)
    higher <- max(q + 10L, 2L * q - 1L)
    if (inherits(tryCatch(gauss_hermite_normal(higher), error = identity), "error")) higher <- q + 10L
    n <- length(input$prep$levels$Person)
    if (mfrmr_adaptive_integration(cfg)) {
      refined <- mfrmr_make_adaptive_mml_evaluator(input$indices, cfg,
        build_param_sizes(cfg), higher)(fit$opt$par)
      value <- refined$value; gradient <- refined$gradient / n
    } else {
      specification <- fit$gmfrm$specification
      specification$quadrature <- gauss_hermite_normal(higher)
      refined <- do.call(mfrm_gmfrm_problem, specification)$marginal(fit$opt$par)
      value <- -refined$logLik; gradient <- refined$gradient
    }
    local$integration <- data.frame(OriginalOrder = q, ComparisonOrder = higher,
      NLLChangePerPerson = abs(value - fit$opt$value) / n,
      ComparisonMeanGradient = max(abs(gradient)))
    local$eligible <- all(is.finite(unlist(local$integration))) &&
      local$integration$NLLChangePerPerson <= 1e-6 &&
      local$integration$ComparisonMeanGradient <= 1e-6
    local$review <- if (local$eligible) {
      "Saved two-family calibration identity and convergence passed; finer-rule NLL change and maximum absolute gradient per Person are at most 1e-6."
    } else "The two-family calibration needs finer integration or further optimization; numerical source scoring checks did not pass."
    local
  }, error = function(e) {
    local$review <- paste("Two-family source integration checks are unavailable:", conditionMessage(e))
    local
  })
  local$source_identity <- prediction_gmfrm_source_identity(fit)
  state <- mfrm_convergence_state(fit)
  list(ready = isTRUE(local$eligible),
    status = if (isTRUE(local$eligible)) "conditional" else "review_only",
    policy_basis = "two_family_conditional_eap_v1", local_calibration_review = local,
    parameter_ready = TRUE,
    reason_codes = if (isTRUE(local$eligible)) character(0) else "local_calibration_check_failed",
    source_audit_states = c(identification = field("EstimabilityState"),
      boundary = field("BoundaryState"), numerical = field("NumericalState")),
    fit_readiness = state$fit_readiness, inference_ready = isTRUE(state$inference_ready),
    readiness_contract_version = field("ReadinessContractVersion"))
}

prediction_source_scoring_readiness <- function(fit, covariance = NULL) {
  if (mfrm_has_product_slopes(fit)) return(prediction_gmfrm_scoring_readiness(fit))
  state <- mfrm_convergence_state(fit)
  record <- mfrmr_get_readiness_record(fit)
  readiness <- as.data.frame(record$fit %||% data.frame(),
                             stringsAsFactors = FALSE)
  value <- function(field, default = "") {
    if (nrow(readiness) == 1L && field %in% names(readiness)) {
      as.character(readiness[[field]][1])
    } else {
      as.character(default)
    }
  }

  contract_ready <- nrow(readiness) == 1L && identical(
    value("ReadinessContractVersion"),
    mfrmr_readiness_contract_version()
  )

  input_state <- value("InputState", "legacy_unknown")
  preparation_notes <- as.data.frame(
    fit$prep$preparation_notes %||%
      fit$data_review$preparation_notes %||%
      data.frame(),
    stringsAsFactors = FALSE
  )
  review_conditions <- if (
    nrow(preparation_notes) > 0L &&
      all(c("Condition", "Severity") %in% names(preparation_notes))
  ) {
    review_rows <- tolower(as.character(preparation_notes$Severity)) %in%
      c("review", "warning", "warn", "error")
    unique(as.character(preparation_notes$Condition[review_rows]))
  } else {
    character(0)
  }
  review_conditions <- review_conditions[
    !is.na(review_conditions) & nzchar(trimws(review_conditions))
  ]
  # A recorded, deterministic score recoding is part of the scoring schema and
  # does not by itself invalidate fixed-parameter scoring. Other input review
  # states remain fail-closed.
  input_ready <- identical(input_state, "pass") ||
    (identical(input_state, "review") &&
       length(review_conditions) > 0L &&
       all(review_conditions %in% "score_categories_recoded"))

  population_active <- isTRUE(fit$config$population_spec$active) ||
    isTRUE(fit$config$population_active) || isTRUE(fit$population$active) ||
    identical(fit$config$posterior_basis, "population_model")
  estimability_state <- value("EstimabilityState", "legacy_unknown")
  estimability_ready <- estimability_state %in%
    c("identified", "population_assumption_linked")
  category_state <- value("CategoryState", "legacy_unknown")
  category_ready <- category_state %in% c("adequate", "not_applicable")
  boundary_state <- value("BoundaryState", "legacy_unknown")
  boundary_ready <- boundary_state %in% c("finite", "not_applicable")
  numerical_state <- value("NumericalState", "legacy_unknown")
  numerical_ready <- identical(numerical_state, "ready")

  sizes <- tryCatch(build_param_sizes(fit$config), error = function(e) NULL)
  par <- fit$opt$par %||% NULL
  parameter_ready <- is.numeric(par) && length(par) > 0L &&
    all(is.finite(par)) && !is.null(sizes) &&
    length(par) == sum(as.integer(unlist(sizes, use.names = FALSE)))

  # Finite optimizer coordinates alone do not ensure finite expanded slopes,
  # anchors or population variance. Invalid calibrations cannot be reviewed into use.
  expanded <- if (parameter_ready) tryCatch(expand_params(par, sizes, fit$config),
    error = function(e) NULL) else NULL
  parameter_ready <- parameter_ready && !is.null(expanded) &&
    all(is.finite(unlist(expanded, use.names = FALSE)))
  pop <- materialize_population_spec(fit$config, expanded)
  if (population_active) parameter_ready <- parameter_ready &&
    isTRUE(pop$active) && length(pop$sigma2) == 1L && is.finite(pop$sigma2) &&
    pop$sigma2 > 0 && length(pop$coefficients) > 0L && all(is.finite(pop$coefficients)) &&
    identical(fit$population$fixed_sd, pop$fixed_sd) &&
    isTRUE(all.equal(fit$population$sigma2, pop$sigma2, tolerance = 1e-10)) &&
    isTRUE(all.equal(as.numeric(fit$population$coefficients), as.numeric(pop$coefficients), tolerance = 1e-10)) &&
    identical(as.character(fit$population$design_columns), names(pop$coefficients))
  rank_deficient <- isTRUE(fit$config$estimability_audit$nonlinear_local_estimability$local_first_order_rank_deficient) ||
    isTRUE(fit$data_review$estimability$nonlinear_local_estimability$local_first_order_rank_deficient)
  if (rank_deficient) estimability_ready <- FALSE
  ready <- !population_active && contract_ready && input_ready && estimability_ready &&
    category_ready && boundary_ready && numerical_ready && parameter_ready
  local <- NULL
  normal_intercept <- !population_active || (is.matrix(pop$design_matrix) &&
    ncol(pop$design_matrix) == 1L && all(is.finite(pop$design_matrix)) &&
    all(pop$design_matrix == 1))
  normal_scoring <- normal_intercept || (population_active &&
    fit$config$model %in% c("RSM", "PCM"))
  if (!isTRUE(ready) && identical(fit$config$method, "MML") &&
      (population_active || identical(fit$config$model, "GPCM")) &&
      contract_ready && input_ready && category_ready && numerical_ready && parameter_ready &&
      !rank_deficient && normal_scoring &&
      estimability_state %in% c("identified", "population_assumption_linked", "not_evaluated") &&
      boundary_state %in% c("finite", "not_applicable", "not_evaluated")) {
    local <- prediction_local_calibration_review(fit, covariance = covariance)
    ready <- isTRUE(local$eligible)
  }
  if (!isTRUE(ready) && identical(fit$config$method, "JML") &&
      identical(fit$config$model, "GPCM") && !population_active &&
      contract_ready && input_ready && category_ready && numerical_ready && parameter_ready &&
      !rank_deficient &&
      estimability_state %in% c("identified", "not_evaluated") &&
      boundary_state %in% c("finite", "not_applicable", "not_evaluated")) {
    local <- prediction_jml_calibration_review(fit)
    ready <- isTRUE(local$eligible)
  }
  reason_codes <- character(0)
  if (population_active) {
    if (!isTRUE(local$eligible)) reason_codes <- c(reason_codes,
      if (normal_scoring) "population_calibration_requires_review" else "population_scoring_scope_requires_review")
  }
  if (!contract_ready) {
    reason_codes <- c(reason_codes, "readiness_contract_not_current")
  }
  if (!input_ready) {
    reason_codes <- c(
      reason_codes, paste0("input_not_scoring_ready:", input_state)
    )
  }
  if (!estimability_ready) {
    reason_codes <- c(
      reason_codes,
      paste0("estimability_not_scoring_ready:", estimability_state)
    )
  }
  if (!category_ready) {
    reason_codes <- c(
      reason_codes, paste0("category_not_scoring_ready:", category_state)
    )
  }
  if (!boundary_ready) {
    reason_codes <- c(
      reason_codes, paste0("boundary_not_scoring_ready:", boundary_state)
    )
  }
  if (!numerical_ready) {
    reason_codes <- c(
      reason_codes, paste0("numerical_not_scoring_ready:", numerical_state)
    )
  }
  if (!parameter_ready) {
    reason_codes <- c(reason_codes, "calibration_parameter_layout_invalid")
  }
  if (rank_deficient) reason_codes <- c(reason_codes, "calibration_rank_deficient")
  if (!is.null(local) && !isTRUE(local$eligible)) reason_codes <- c(reason_codes, "local_calibration_check_failed")
  list(
    ready = isTRUE(ready),
    status = if (!parameter_ready) "unavailable" else if (isTRUE(local$eligible))
      "conditional" else if (isTRUE(ready)) "ready" else "review_only",
    policy_basis = "conditional_scoring_checks_v3",
    local_calibration_review = local,
    parameter_ready = isTRUE(parameter_ready),
    reason_codes = if (isTRUE(ready)) character(0) else unique(reason_codes),
    source_audit_states = c(identification = estimability_state, boundary = boundary_state,
      numerical = numerical_state),
    fit_readiness = as.character(state$fit_readiness %||% "unknown"),
    inference_ready = isTRUE(state$inference_ready),
    readiness_contract_version = as.character(
      state$readiness_contract_version %||% ""
    )
  )
}

prediction_scoring_reasons <- function(codes) {
  labels <- c(
    population_scoring_validity_not_evaluated = "Scoring with an estimated population distribution still requires review.",
    population_calibration_requires_review = "The estimated population distribution requires calibration review.",
    population_scoring_scope_requires_review = "Scoring with population covariates requires explicit review.",
    calibration_rank_deficient = "The source calibration has a detected local identification deficiency.",
    local_calibration_check_failed = "The source calibration did not pass the local solution and integration checks.",
    readiness_contract_not_current = "The source fit does not retain current scoring checks.",
    input_not_scoring_ready = "The source response data require review.",
    estimability_not_scoring_ready = "The source calibration's identification requires review.",
    category_not_scoring_ready = "The source score categories require review.",
    boundary_not_scoring_ready = "A source estimate may lie at a boundary or be unbounded.",
    numerical_not_scoring_ready = "Source numerical convergence requires review.",
    calibration_parameter_layout_invalid = "The retained calibration parameters are incomplete or incompatible."
  )
  codes <- as.character(codes)
  result <- unname(labels[sub(":.*$", "", codes)])
  result[codes == "boundary_not_scoring_ready:not_evaluated"] <-
    "The source calibration's boundary behavior has not been evaluated."
  result[codes == "estimability_not_scoring_ready:not_evaluated"] <-
    "The source calibration's identification has not been evaluated."
  result[is.na(result)] <- "Inspect the source calibration diagnostics before interpreting these scores."
  unique(result)
}

prediction_validate_scoring_prior <- function(prior) {
  if (is.null(prior)) return(NULL)
  scalar <- function(z) is.numeric(z) && is.null(dim(z)) && length(z) == 1L && is.finite(z)
  if (!is.list(prior) || is.data.frame(prior) || length(prior) != 2L ||
      !setequal(names(prior), c("mean", "sd")) ||
      !scalar(prior$mean) || !scalar(prior$sd) || prior$sd <= 0 ||
      !is.finite(prior$sd^2) || prior$sd^2 <= 0) {
    stop("`scoring_prior` must be NULL or list(mean = <finite number>, sd = <positive finite number>) with a representable positive variance.", call. = FALSE)
  }
  list(mean = unname(prior$mean), sd = unname(prior$sd))
}

prediction_validate_interval_record <- function(x) {
  level <- x$settings$interval_level
  algorithm <- x$settings$scoring_algorithm
  valid <- is.numeric(level) && length(level) == 1L && is.finite(level) &&
    level > 0 && level < 1 && is.character(algorithm) && length(algorithm) == 1L &&
    algorithm %in% c("quadrature_eap_v2", "adaptive_quadrature_eap_v2") &&
    is.data.frame(x$estimates) &&
    all(c("IntervalLevel", "ScoringAlgorithm") %in% names(x$estimates)) &&
    identical(unname(x$estimates$IntervalLevel), rep(unname(level), nrow(x$estimates))) &&
    identical(unname(x$estimates$ScoringAlgorithm), rep(unname(algorithm), nrow(x$estimates)))
  if (!isTRUE(valid))
    stop("Scoring output has inconsistent interval records. Re-score to change the interval level or algorithm; saved bounds cannot be relabelled.", call. = FALSE)
  invisible(x)
}

prediction_estimate_table <- function(x) {
  table <- x$estimates
  if (!is.data.frame(table) || !nrow(table)) return(table)
  settings <- x$settings
  supplied <- identical(settings$posterior_basis, "user_supplied_normal")
  population <- identical(settings$retained_posterior_basis %||% settings$posterior_basis, "population_model")
  table$IntervalLevel <- settings$interval_level %||% NA_real_
  table$ScoringAlgorithm <- settings$scoring_algorithm %||% NA_character_
  corrected <- identical(settings$source_scoring_policy_basis,"corrected_jml_reference_eap_v1")
  table$CalibrationMethod <- if (corrected) "Corrected JML" else settings$method %||% NA_character_
  if (corrected) table$CorrectionOrder <- settings$local_calibration_review$correction_order
  table$EstimateBasis <- "Posterior EAP conditional on point calibration and scoring prior"
  table$Prior <- if (supplied) "User-supplied normal" else if (population) "Fitted conditional normal" else "Standard normal N(0,1)"
  if (!is.null(settings$prior_comparison)) {
    table$PriorSource <- if (supplied) "user_supplied" else "retained"
    table$RetainedPrior <- if (population) "Fitted conditional normal" else "Standard normal reference N(0,1)"
  }
  if (!"PriorMean" %in% names(table)) table$PriorMean <- if (population) NA_real_ else 0
  if (!"PriorSD" %in% names(table)) table$PriorSD <- if (population) NA_real_ else 1
  if (!"WeightedLikelihood" %in% names(table)) table$WeightedLikelihood <- NA
  table$UncertaintyBasis <- "Conditional posterior uncertainty; calibration and population estimation uncertainty excluded"
  algorithm <- settings$scoring_algorithm %||% ""
  table$IntervalBasis <- if (endsWith(algorithm, "_v2")) {
    "Continuous equal-tail posterior quantiles"
  } else if (endsWith(algorithm, "_v1")) {
    "Quadrature-grid endpoints; posterior mass may differ from requested level"
  } else "Interval calculation not recorded; rescore to establish its interpretation"
  table
}

prediction_draw_table <- function(draws, estimates) {
  if (!is.data.frame(draws) || !nrow(draws)) return(draws)
  columns <- intersect(c("CalibrationMethod", "CorrectionOrder", "Prior", "PriorMean", "PriorSD",
    "PriorSource", "RetainedPrior", "RetainedPriorMean", "RetainedPriorSD",
    "WeightedLikelihood", "UncertaintyBasis", "ScoringAlgorithm",
    "SourceScoringReady", "ScoreIntegrationReady", "EstimateUse"), names(estimates))
  for (name in columns) draws[[name]] <- estimates[[name]][match(draws$Person, estimates$Person)]
  draws$DrawBasis <- "Discrete quadrature posterior; conditional on point calibration and scoring prior"
  draws
}

prediction_output_notes <- function(x) {
  notes <- x$notes %||% character(0)
  notes <- notes[!grepl("Reason codes:", notes, fixed = TRUE) &
                   !startsWith(notes, "Intervals invert the continuous posterior CDF")]
  notes <- c(notes, if (identical(x$settings$posterior_basis, "user_supplied_normal")) {
    paste0("Scoring prior: user-supplied normal, mean = ", x$settings$scoring_prior$mean,
      ", SD = ", x$settings$scoring_prior$sd,
      ". The calibration and its original population assumptions are unchanged; this prior is an analyst assumption, not a newly estimated population.")
  } else if (identical(x$settings$posterior_basis, "population_model")) {
    "Scoring prior: fitted conditional normal distribution, held at its estimated parameters."
  } else "Scoring prior: standard normal N(0,1).")
  notes <- c(notes, mfrmr_calibration_score_interval_note(x$settings),
    "Posterior SDs and intervals condition on the point calibration and prior; uncertainty from estimating either is excluded.")
  if (identical(x$settings$source_scoring_status, "conditional")) {
    notes <- c(notes,
      if (identical(x$settings$source_scoring_policy_basis, "two_family_conditional_eap_v1")) {
        "The saved two-family calibration passed identity, convergence and finer-integration checks. Scores condition on these estimates; this is not a global-maximum or future-population validation."
      } else "The calibration passed local solution and integration checks. Scores condition on these estimates; this is not a global-maximum or future-population validation.",
      x$settings$local_calibration_review$caution %||% character(0))
  }
  if (!is.null(x$settings$score_integration_review) &&
      !all(x$settings$score_integration_review$Passed)) {
    notes <- c(notes, "Flagged scores require integration review; increase scoring_quad_points before interpreting them.")
  }
  if (identical(x$settings$source_scoring_ready, FALSE)) {
    notes <- c(notes, "The source fit requires review; these scores are review-only and should be treated as exploratory.",
      prediction_scoring_reasons(x$settings$source_scoring_reason_codes),
      x$settings$local_calibration_review$review %||% character(0))
  }
  unique(notes)
}

prediction_display_estimates <- function(table) {
  columns <- intersect(c("Person", "Estimate", "SD", "Lower", "Upper", "IntervalLevel",
                          "PriorMean", "PriorSD", "Observations"), names(table))
  out <- as.data.frame(table[columns])
  if ("SourceScoringReady" %in% names(table)) {
    out$Review <- ifelse(is.na(table$SourceScoringReady), "Not recorded",
                        ifelse(table$SourceScoringReady, "Source scoring checks passed", "Required"))
  }
  if ("ScoreIntegrationReady" %in% names(table)) {
    out$Review[!table$ScoreIntegrationReady] <- "Integration review required"
  }
  out
}

prediction_print_summary <- function(x, plausible = FALSE) {
  prediction_validate_population_output(x)
  cat(if (plausible) "mfrmr Plausible Values Summary\n" else "mfrmr Unit Prediction Summary\n")
  estimates <- prediction_estimate_table(x)
  print_wrapped_line(paste0("Calibration estimated by ", unique(estimates$CalibrationMethod)[1] %||% "an unrecorded method",
    "; scoring uses posterior EAP. Prior: ", unique(estimates$Prior)[1], "."))
  print_wrapped_line(mfrmr_calibration_score_interval_note(x$settings))
  print_wrapped_line("Posterior SDs and intervals condition on point estimates of the calibration and prior; their estimation uncertainty is excluded.")
  if (identical(x$settings$posterior_basis, "user_supplied_normal")) {
    print_wrapped_line(paste0("User-supplied scoring prior: mean = ", x$settings$scoring_prior$mean,
      ", SD = ", x$settings$scoring_prior$sd,
      ". The retained prior is recorded separately; the new prior was not estimated from these responses."))
  }
  if (!is.null(x$quadrature_overview)) {
    cat("\nFixed-parameter integration review (adaptive minus fixed)\n")
    print(x$quadrature_overview, row.names = FALSE)
  }
  if (plausible) {
    cat("\nEmpirical draw summaries (first 10)\n")
    print(utils::head(as.data.frame(x$draw_summary[c("Person", "Draws", "MeanValue",
      "SDValue", "LowerValue", "UpperValue")]), 10L), row.names = FALSE)
    print_wrapped_line("Draw limits are empirical quantiles at the requested level. With few draws they are coarse; use the companion posterior interval and its stated calculation method for interval reporting.")
  }
  cat("\nPosterior estimates (first 10)\n")
  shown <- prediction_display_estimates(estimates)
  shown <- shown[setdiff(names(shown), c("IntervalLevel", "PriorMean", "PriorSD"))]
  print(utils::head(shown, 10L), row.names = FALSE)
  if (!is.null(x$decisions)) {
    mfrm_validate_decisions(x)
    cat("\nConditional decision probabilities (first 10)\n")
    print(utils::head(x$decisions,10L), row.names=FALSE)
    print_wrapped_line(mfrm_decision_note())
  }
  if (nrow(x$row_review %||% data.frame())) {
    cat("\nResponse rows\n"); print(as.data.frame(x$row_review), row.names = FALSE)
  }
  if (nrow(x$population_review %||% data.frame())) {
    cat("\nBackground-data omissions\n")
    columns <- intersect(c("InputPersons", "RetainedPersons", "OmittedPersons", "RetainedRows", "OmittedRows"), names(x$population_review))
    print(as.data.frame(x$population_review[columns]), row.names = FALSE)
  }
  notes <- prediction_output_notes(x)
  notes <- notes[!grepl("^Posterior summaries|^Scoring prior:|^[0-9.]+% intervals:|^Intervals:|^Posterior SDs and intervals", notes)]
  for (note in notes) print_wrapped_line(note)
  invisible(x)
}

prediction_validate_population_output <- function(x) {
  if (inherits(x, "summary.mfrm_plausible_values") &&
      !all(c("IntervalLevel", "DrawSummaryBasis") %in% names(x$draw_summary))) {
    stop("Recreate this summary with summary(original_plausible_values) to label the draw quantiles and interval level; no refitting or resampling is needed.", call. = FALSE)
  }
  if (!inherits(x, c("mfrm_unit_prediction", "mfrm_plausible_values",
                      "summary.mfrm_unit_prediction", "summary.mfrm_plausible_values"))) return(invisible(x))
  product <- identical(x$settings$source_scoring_policy_basis, "two_family_conditional_eap_v1") ||
    !is.null(x$settings$two_family_calibration)
  if (product) {
    prediction_validate_interval_record(x)
    calibration <- x$settings$two_family_calibration
    owners <- calibration$slope_facets
    valid <- identical(x$settings$source_scoring_policy_basis, "two_family_conditional_eap_v1") &&
      x$settings$source_scoring_status %in% c("conditional", "review_only") &&
      identical(isTRUE(x$settings$source_scoring_ready),
        identical(x$settings$source_scoring_status, "conditional")) &&
      identical(isTRUE(x$settings$source_scoring_ready),
        isTRUE(x$settings$local_calibration_review$eligible)) &&
      identical(x$settings$method, "MML") &&
      identical(x$settings$posterior_basis, "fixed_standard_normal") &&
      is.null(x$settings$scoring_prior) &&
      identical(calibration$scale_reference, "fixed_standard_normal") &&
      is.character(owners) && length(owners) == 2L && !anyNA(owners) && !anyDuplicated(owners) &&
      identical(calibration$step_facet, owners[2L]) &&
      identical(names(calibration$levels), owners) &&
      identical(names(calibration$slope_components), owners) &&
      all(vapply(calibration$slope_components, function(tab)
        is.data.frame(tab) && is.numeric(tab$Slope) && length(tab$Slope) > 1L &&
          all(is.finite(tab$Slope) & tab$Slope > 0), TRUE)) &&
      is.data.frame(calibration$score_map) &&
      all(c("OriginalScore", "InternalScore") %in% names(calibration$score_map)) &&
      identical(calibration$score_map$OriginalScore, calibration$score_map$InternalScore)
    if (isTRUE(x$settings$source_scoring_ready)) {
      check <- x$settings$local_calibration_review$integration
      valid <- valid && is.data.frame(check) && nrow(check) == 1L &&
        all(c("NLLChangePerPerson", "ComparisonMeanGradient") %in% names(check)) &&
        all(is.finite(unlist(check))) && check$NLLChangePerPerson >= 0 &&
        check$NLLChangePerPerson <= 1e-6 && check$ComparisonMeanGradient >= 0 &&
        check$ComparisonMeanGradient <= 1e-6
    }
    if (!isTRUE(valid)) stop("Two-family scoring output has missing or inconsistent calibration records. Re-score from the fitted model before summary or export.", call. = FALSE)
  }
  supplied <- identical(x$settings$posterior_basis, "user_supplied_normal")
  if (supplied || !is.null(x$settings$scoring_prior) || !is.null(x$settings$prior_comparison) ||
      any(x$estimates[["PriorSource"]] %in% "user_supplied")) {
    prior <- prediction_validate_scoring_prior(x$settings$scoring_prior)
    p <- x$settings$prior_comparison
    columns <- c("RetainedPriorMean", "RetainedPriorSD", "PriorMean", "PriorSD")
    valid <- identical(supplied, !is.null(prior)) &&
      is.character(x$settings$retained_posterior_basis) &&
      length(x$settings$retained_posterior_basis) == 1L && !is.na(x$settings$retained_posterior_basis) &&
      x$settings$retained_posterior_basis != "user_supplied_normal" &&
      is.data.frame(p) && all(c("Person", columns) %in% names(p)) &&
      all(c(columns, "PriorSource") %in% names(x$estimates)) &&
      identical(as.character(p$Person), as.character(x$estimates$Person)) &&
      all(vapply(p[columns], is.numeric, TRUE)) &&
      all(is.finite(as.matrix(p[columns]))) && all(p$RetainedPriorSD > 0 & p$PriorSD > 0)
    if (valid) valid <- all(vapply(columns, function(nm)
      identical(as.numeric(p[[nm]]), as.numeric(x$estimates[[nm]])), TRUE)) &&
      all(x$estimates$PriorSource == if (supplied) "user_supplied" else "retained") &&
      if (supplied) {
        all(p$PriorMean == prior$mean) &&
          isTRUE(all.equal(p$PriorSD, rep(prior$sd, nrow(p)), tolerance = 1e-14))
      } else all(p$PriorMean == p$RetainedPriorMean & p$PriorSD == p$RetainedPriorSD)
    if (!isTRUE(valid)) stop("Scoring output has missing or inconsistent retained/scoring prior records. Re-score before summary or export.", call. = FALSE)
  }
  conditional <- identical(x$settings$source_scoring_status, "conditional")
  population <- identical(x$settings$posterior_basis, "population_model") ||
    identical(x$settings$retained_posterior_basis, "population_model")
  if (population && !conditional &&
      (!identical(x$settings$source_scoring_ready, FALSE) ||
       !identical(x$settings$readiness_policy, "review") ||
       !identical(x$settings$source_scoring_status, "review_only") ||
       !all(c("SourceScoringReady", "EstimateUse") %in% names(x$estimates)) ||
       anyNA(x$estimates$SourceScoringReady) || any(x$estimates$SourceScoringReady) ||
       anyNA(x$estimates$EstimateUse) ||
       any(x$estimates$EstimateUse != "review_only_nonready_source"))) {
    stop("Population-model scoring output requires explicit review-only eligibility. ",
         "Regenerate it with `readiness_policy = \"review\"` before summary or export.",
         call. = FALSE)
  }
  if (conditional || supplied || product || isTRUE(x$settings$score_integration_required) ||
      !is.null(x$settings$score_integration_review) || "ScoreIntegrationReady" %in% names(x$estimates)) {
    review <- x$settings$score_integration_review
    valid_review <- is.data.frame(review) && nrow(review) == nrow(x$estimates) &&
      all(c("Person", "Passed", "EAPChange", "SDChange", "ReferenceChange") %in% names(review)) &&
      identical(as.character(review$Person), as.character(x$estimates$Person)) &&
      !anyNA(review$Passed)
    if (valid_review) {
      actual <- apply(as.matrix(review[c("EAPChange", "SDChange", "ReferenceChange")]), 1L,
        function(z) all(is.finite(z)) && all(z >= 0 & z <= 1e-5))
      valid_review <- identical(unname(actual), review$Passed)
    }
    ready <- isTRUE(x$settings$source_scoring_ready)
    valid <- valid_review &&
      (!conditional || (x$settings$source_scoring_policy_basis %in% c("conditional_scoring_checks_v3","corrected_jml_reference_eap_v1", "two_family_conditional_eap_v1") &&
        ready && isTRUE(x$settings$local_calibration_review$eligible))) &&
      all(c("SourceScoringReady", "ScoreIntegrationReady", "EstimateUse") %in% names(x$estimates)) &&
      all(!is.na(x$estimates$SourceScoringReady) & x$estimates$SourceScoringReady == ready) &&
      identical(x$estimates$ScoreIntegrationReady, review$Passed) &&
      identical(x$estimates$EstimateUse, if (!ready) rep("review_only_nonready_source", nrow(review)) else
        ifelse(review$Passed, "fitted_object_scoring", "review_only_scoring_integration")) &&
      ((ready && all(review$Passed)) || identical(x$settings$readiness_policy, "review"))
    if (conditional && identical(x$settings$method, "JML")) {
      evidence <- list(policy_basis = x$settings$source_scoring_policy_basis,
        status = x$settings$source_scoring_status,
        local_calibration_review = x$settings$local_calibration_review,
        source_audit_states = x$settings$source_audit_states,
        inference_ready = x$settings$source_inference_ready)
      corrected <- identical(evidence$policy_basis,"corrected_jml_reference_eap_v1")
      valid <- valid && (if (corrected) mfrm_jml_scoring_evidence_valid else mfrmr_calibration_gpcm_jml_evidence_valid)(evidence)
      if (corrected) valid <- valid && all(c("CalibrationMethod","CorrectionOrder") %in% names(x$estimates)) &&
        isTRUE(all(x$estimates$CalibrationMethod=="Corrected JML" & x$estimates$CorrectionOrder==evidence$local_calibration_review$correction_order))
    }
    if (!valid) stop("Scoring output has missing or inconsistent numerical review records. Re-score with the fitted model before summary or export.", call. = FALSE)
    return(invisible(x))
  }
  invisible(x)
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
#'   ordinary fitted-object posterior scoring. For intercept-only latent-regression
#'   fits (`population_formula = ~ 1`), `mfrmr` reconstructs the minimal
#'   one-row-per-person table internally from the scored person IDs. This is the
#'   scoring-time table for `new_data`, not the fit object's replay/export
#'   provenance table. For categorical background variables, supply values on
#'   the same coding scale used at fit time; the fitted factor levels and
#'   contrasts are reused when building the scoring design matrix.
#'   Fitted standardization and polynomial/spline bases are also reused;
#'   they are not re-estimated from the scoring cohort.
#' @param person_id Optional person-ID column in `person_data`. Defaults to
#'   `person` when that column exists, otherwise `"Person"` for the canonical
#'   scoring layout.
#' @param population_policy How missing background data are handled when
#'   `fit` uses the latent-regression branch. `"error"` (default) requires
#'   complete person-level covariates for all scored persons; `"omit"` drops
#'   scored persons lacking complete covariates and records that omission in
#'   `population_review`.
#' @param interval_level Posterior interval level returned in `Lower`/`Upper`.
#' @param scoring_quad_points Number of Gauss-Hermite nodes used only for this
#'   scoring call. It is independent of the quadrature order used while fitting
#'   `fit`; the default is 31 and values below 2 are refused.
#'   The fixed or adaptive integration mode is inherited from `fit`.
#' @param readiness_policy How a source fit that is not scoring-ready is
#'   handled. `"error"` (default) refuses scoring. `"review"` permits an
#'   explicitly review-only fitted-object calculation and labels the returned
#'   estimates and settings accordingly; it does not make the source fit ready.
#'   GPCM MML and RSM/PCM normal population models may pass the
#'   conditional-scoring checks described below even when their global inference
#'   audit remains incomplete. GPCM population models with background covariates
#'   still require explicit review. Invalid calibration parameters or inconsistent
#'   stored prior parameters are refused under either policy. Older population
#'   results without the new check records retain their review-only restriction.
#' @param n_draws Optional number of quadrature-grid posterior draws to return
#'   per scored person. Use 0 to skip draws.
#' @param seed Optional seed for reproducible posterior draws.
#' @param adaptive_quad_points Optional vector of at least two distinct integer
#'   orders >= 3, for example `c(31, 61)`. Adds `quadrature_review`, comparing
#'   a fixed-prior grid with grids centered and scaled to each Person's
#'   posterior. Calibration parameters and the prior are held fixed. Inspect
#'   both fixed/adaptive differences and movement between adaptive orders;
#'   this diagnostic does not replace estimates, intervals or draws. Every
#'   scoring call automatically adds reference orders and uses them to check
#'   numerical score accuracy, including ordinary RSM/PCM with the retained
#'   prior; see below.
#' @param scoring_prior Optional `list(mean = ..., sd = ...)` specifying one
#'   normal scoring prior shared by all scored persons, on the unchanged
#'   calibration ability scale. `sd` is a standard deviation, not a variance.
#'   `NULL` (default) retains the MML scoring prior; JML uses a post-hoc
#'   standard-normal reference prior. Both values must be
#'   finite; SD and its representable variance must be positive. Overrides of
#'   population models with background covariates are not supported. A supplied
#'   prior does not refit the calibration, estimate a new population, or bypass
#'   source checks. Numerical scoring checks run under the supplied prior.
#' @param cut_scores Optional finite, strictly increasing ability cuts with
#'   optional distinct names. Requires ordinary RSM/PCM and unit weights.
#'   A `decisions` table returns continuous probabilities above cuts and within
#'   bands (`Lower < ability <= Upper`), conditional on the point calibration
#'   and actual scoring prior. Calibration/population and cut-score estimation
#'   uncertainty are excluded. This is not a Person-difference test or a
#'   probability of passing a future administration. Cut labels alone do not
#'   establish CEFR alignment. Tiny tails have absolute, not relative, numerical
#'   accuracy. Use `apa_table(x, which = "decisions")` or export the saved table.
#'
#' @details
#' `predict_mfrm_units()` is the **individual-unit companion** to
#' [predict_mfrm_population()]. It uses the fitted calibration and, when
#' available, the fitted one-dimensional population model to score new or
#' partially observed persons via Expected A Posteriori (EAP) summaries on a
#' quadrature grid.
#'
#' With the default `scoring_prior = NULL`, when the original fit uses ordinary
#' `method = "MML"`, the posterior
#' summaries use that fitted MML calibration and a standard normal scoring
#' prior, unless a population model was fitted. When the original fit
#' uses the latent-regression MML branch, the scoring prior is the fitted
#' conditional normal population model \eqn{\theta \mid x \sim
#' N(x^\top\hat\beta, \hat\sigma^2)}, so the returned summaries are
#' population-model-aware posterior EAP estimates. When the original fit uses
#' `method = "JML"`, `mfrmr` applies the fitted facet/step parameters with a
#' standard normal reference prior on the quadrature grid, so the returned
#' person scores remain fitted-object EAP summaries rather than direct JML
#' estimates from the fitting step.
#'
#' When the fitted population model is intercept-only (`population_formula = ~
#' 1`), `predict_mfrm_units()` still uses the fitted population-model basis,
#' but it can reconstruct the minimal scored-person table internally because no
#' background covariates are needed beyond the person IDs in `new_data`.
#'
#' The current `GPCM` branch is included in this scoring layer,
#' so fitted `GPCM` objects can be used for the same fitted-object
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
#' adds a standard normal reference prior by default (or the explicit
#' `scoring_prior`) only for the scoring layer so that
#' new or partially observed units can be summarized on a quadrature grid.
#' This is a practical fitted-object EAP procedure, not a claim that the
#' original `JML` fit itself estimated a population model. It does not reproduce
#' the training Person joint maximum-likelihood estimates or provide ML/WLE
#' scoring for new Persons. Posterior uncertainty conditions on the fitted
#' calibration; it does not include calibration-parameter uncertainty.
#'
#' @section Experimental two-family scoring:
#' Two-family GPCM MML supports new Persons with known non-Person facet levels,
#' unit weights and the retained fixed N(0,1) prior. Both slope components,
#' facet locations and the second owner's category steps remain fixed. New
#' combinations of known levels are allowed under the fitted product model;
#' this does not validate those combinations empirically. Repeated Person-facet
#' combinations are rejected. Missing responses are omitted with a warning and
#' row review; Persons without valid responses receive no prior-only score.
#' Unknown levels or scores outside the original category coding are rejected.
#' Prior overrides and person-level population inputs are not supported.
#'
#' The saved specification, point estimates, objective and convergence must
#' agree before either readiness policy permits scoring. A finer rule (normally
#' `max(q + 10, 2 * q - 1)`) must change the source negative log likelihood per
#' Person by at most 1e-6 and have maximum absolute gradient per Person at most
#' 1e-6 for conditional source scoring. Source-readiness or refinement failures
#' require explicit `readiness_policy = "review"`; invalid calibration records
#' or failed convergence cannot be enabled by review. No covariance is required.
#' Each new batch also receives the EAP/SD integration checks described below.
#' These tolerances are numerical checks, not statistical validity thresholds.
#'
#' EAP, posterior SD and continuous equal-tail intervals condition on the point
#' calibration and fixed prior, excluding calibration uncertainty. Numerical
#' checks do not complete the identification/boundary audits or establish
#' sampling coverage, population transport or decision accuracy. Settings retain
#' both owners, component slopes, locations, steps, category coding, observed
#' calibration combinations and separate source/batch checks through saveRDS()
#' and summary(). [extract_mfrm_calibration()] provides the corresponding
#' portable format 6 route after passing source checks. Use
#' `mfrm_results(fit, scores = scores)` to attach matching saved native or
#' format-6 scores for reports and exports without rescoring. The saved source
#' identity, calibration and scoring records must agree with the supplied fit;
#' source/batch checks, omitted rows and review-only labels remain visible.
#'
#' @section Conditional scoring and source checks:
#' A finite calibration can define a person's posterior score without proving
#' that the calibration is the unique global optimum. For one-family GPCM MML and models
#' with an estimated RSM/PCM normal population, this function checks the
#' current likelihood, gradient and unregularized joint information using the
#' same local-solution machinery as [compare_mfrm()]. Detected identification
#' deficiencies, inadequate categories and failed convergence do not pass.
#' RSM/PCM regressions additionally require the consistent full-rank source,
#' unit weights, direct MML and centered unanchored facets/steps used by
#' [mfrm_population_intervals()]. Scoring reuses the fitted predictor coding
#' and conditional residual variance; no population uncertainty is propagated.
#' Boundary and identification audits labelled "not evaluated" remain so;
#' they are not relabelled as completed by these local checks.
#'
#' Shared-owner GPCM JML can also use conditional scoring with unit weights and
#' no anchors/interactions. Fresh joint objective/gradient and unregularized
#' full-curvature checks must pass; known Person/additive/slope boundary
#' certificates are not overridden by finite optimizer coordinates. Incomplete
#' global identification/boundary audits remain incomplete. Source curvature
#' includes free Person parameters and uses a dense matrix; saved portable
#' calibration avoids repeating it. No MML integration check or formal slope
#' interval is inferred from this JML calculation.
#'
#' One-family MML calibration checks also compare the retained solution at its fitting order
#' and a higher order (normally `2 * q - 1`, at least `q + 10`; `q + 10` is
#' used if the larger rule is not representable). They require NLL movement
#' <= 1e-5, gradient movement <= 1e-4 and a higher-order gradient <= 1e-4,
#' in addition to the existing local-solution criteria. These are numerical
#' screening tolerances, not statistical coverage or global-existence guarantees.
#'
#' For every native scoring route, after a calibration passes, the actual reported EAP and SD are compared
#' with adaptive reference orders under the same calibration and scoring prior.
#' The reference includes the scoring order (at least 3 for this comparison)
#' and a higher order chosen by the same rule. EAP/SD discrepancies and changes
#' in adaptive EAP, SD and log marginal likelihood must each be <= 1e-5.
#' Increase `scoring_quad_points` if this check fails. `readiness_policy =
#' "review"` retains and labels those scores for inspection instead of stopping.
#' Nonfinite posterior calculations or invalid parameters cannot be enabled by review.
#'
#' `settings$local_calibration_review` records the source checks and
#' `settings$score_integration_review` records the per-Person numerical checks.
#' New outputs set `settings$score_integration_required = TRUE`; their saved
#' checks must remain with the scores. Older outputs without these checks can
#' still be summarized, but are not retrospectively numerically qualified.
#' `SourceScoringReady` concerns the calibration; `ScoreIntegrationReady`
#' concerns the reported scores. `EstimateUse` identifies any scores requiring
#' integration review. These records accompany summaries, draws and exports.
#' The conditional route has `source_scoring_status = "conditional"` and does
#' not change the fit's global `InferenceReady` or boundary/identification status.
#' Information-matrix work respects `mfrmr.max_information_bytes`; unavailable
#' checks require explicit review, rather than accepting a stale covariance.
#'
#' For example, a rubric calibration that passes these checks can score new
#' respondents while holding its estimated normal prior fixed. Passing does
#' not establish that the next cohort has the same ability distribution.
#' Calibration and prior estimation uncertainty remain excluded from posterior
#' intervals, and no portable GPCM calibration artifact is created here.
#'
#' @section Scoring a later cohort:
#' By default the scoring prior is retained from the fitted model; this function does not
#' estimate a new cohort's population mean or spread from `new_data`. For
#' example, a calibration from an advanced class may not supply a suitable
#' prior for a beginning class, even when both classes use the same rubric.
#' With only a few ratings, scores can depend substantially on that prior.
#' The number of ratings alone is not an information measure: their calibrated
#' difficulties, discriminations and observed categories also matter.
#' Inspect extreme response patterns separately; increasing the rating count
#' does not guarantee less prior sensitivity for every person.
#'
#' A narrower prior can yield narrower posterior intervals without any new
#' rating evidence. This is greater certainty under a stronger assumption,
#' not evidence that the assessment became more informative. Changing the
#' prior while fixing calibration is a sensitivity analysis; estimating the
#' calibration again changes a different part of the analysis. Neither action
#' alone establishes coverage for a new cohort.
#'
#' Use `scoring_prior = list(mean = 0.5, sd = 1.5)`, for example, to inspect
#' another normal-prior assumption while retaining the calibration. Choose
#' plausible values for your setting; these example numbers are not a default
#' recommendation. Compare the same response rows under both priors. The
#' supplied prior is labelled as an analyst assumption, even if its numbers
#' equal the retained prior. It does not repair a poorly estimated calibration.
#'
#' `PriorMean`/`PriorSD` describe the prior actually used for scoring;
#' `RetainedPriorMean`/`RetainedPriorSD` describe the original scoring prior.
#' `PriorSource` distinguishes `"retained"` from `"user_supplied"`.
#' `RetainedPrior` distinguishes an estimated normal population from the
#' standard normal reference (including post hoc JML scoring). These unrounded
#' columns accompany estimates, posterior draws and their summaries.
#' `settings$prior_comparison` records both priors by person, and
#' `settings$retained_posterior_basis` records the original basis. Do not edit
#' `fit$population`: its consistency with the retained calibration is checked
#' independently of `scoring_prior`. Saving a prediction with `saveRDS()` stores
#' a result; it does not create a portable GPCM calibration for future scoring.
#'
#' @section Interpreting output:
#' - `estimates` contains posterior EAP summaries for each person in
#'   `new_data`.
#' - `Lower` and `Upper` are continuous equal-tail posterior interval bounds at
#'   the requested `interval_level`, computed by numerical CDF inversion.
#'   EAP, SD, and optional draws still use the selected quadrature rule.
#'   These intervals condition on the fitted calibration and scoring prior,
#'   including any estimated population coefficients and variance. They exclude
#'   uncertainty from estimating calibration or population parameters; they are
#'   not confidence intervals at each fixed true Person ability.
#'   Their interpretation also depends on the scoring prior: a new population
#'   with a different ability mean, spread or shape can have different coverage.
#'   More quadrature points check numerical approximation under the same prior;
#'   they do not establish that this prior matches the new population.
#' - `SD` is posterior uncertainty under the fitted scoring basis used for
#'   scoring.
#' - Estimate tables retain `IntervalLevel` without rounding, the prior form,
#'   per-Person `PriorMean` and `PriorSD`, calibration method, interval method
#'   and uncertainty interpretation. `WeightedLikelihood` indicates whether
#'   any response contribution for that Person was raised to a non-unit weight.
#'   Such weights do not by themselves establish equivalent independent ratings
#'   or frequentist coverage. Draw tables retain their discrete-grid basis and
#'   the same prior and calibration interpretation.
#' - `draws`, when requested, contains approximate plausible values on the
#'   fitted quadrature grid.
#' - `population_review`, when present, records whether scored persons were
#'   omitted because their background data were incomplete for a
#'   latent-regression fit.
#' - `quadrature_review`, when requested, retains unrounded per-Person
#'   fixed/adaptive log-marginal, EAP and posterior-SD differences, changes
#'   between adaptive orders, and computation status/reasons. `summary()` also
#'   supplies a compact `quadrature_overview`. A `computed` status does not
#'   certify accuracy; inspect the differences and any unavailable rows.
#'
#' Re-summarize older results to recover recorded interval settings and readable
#' notes without changing numerical scores. Prior means/SDs for older
#' estimated-population results may be unavailable because those parameters
#' were not retained in the result. Re-score with the existing fitted model
#' to retain them; no calibration refit is needed. An unrecorded interval
#' algorithm is labelled unavailable, not assumed to use continuous quantiles.
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
#' (Bock & Mislevy, 1982, pp. 432-433; Bock & Aitkin, 1981).
#' When `fit` uses the latent-regression
#' branch, `mfrmr` scores under the fitted conditional normal population model
#' in the general plausible-values spirit discussed by Mislevy (1991). Optional
#' posterior draws are exposed as quadrature-grid plausible-value-style
#' summaries for practical many-facet scoring rather than as a claim of full
#' ConQuest numerical equivalence. When the source fit is `JML`, the same
#' literature supports
#' the quadrature-based scoring layer, but the standard normal prior is a
#' package-level reference prior introduced for post hoc scoring rather than an
#' estimated population distribution.
#' The local JML objective, gradient and curvature cutoffs are package
#' numerical criteria, not literature-derived guarantees of calibration
#' accuracy or interval coverage. The scoring literature does not establish
#' the completeness of the source boundary audit.
#'
#' - Bock, R. D., & Mislevy, R. J. (1982). *Adaptive EAP estimation of ability
#'   in a microcomputer environment*. Applied Psychological Measurement,
#'   6(4), 431-444. \doi{10.1177/014662168200600405}.
#' - Bock, R. D., & Aitkin, M. (1981). *Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm*. Psychometrika, 46(4),
#'   443-459.
#' - Mislevy, R. J. (1991). *Randomization-based inference about latent
#'   variables from complex samples*. Psychometrika, 56(2), 177-196.
#' - Muraki, E. (1992). *A generalized partial credit model: Application of an
#'   EM algorithm*. Applied Psychological Measurement, 16(2), 159-176.
#'
#' @section Corrected JML:
#' An experimental shared-owner corrected GPCM calibration can score new
#' Persons after fresh checks of its adjusted equation, full Jacobian and
#' saved parameter identities. These checks never apply an unadjusted
#' likelihood-stationarity requirement to the corrected root. Unit weights
#' and known facet levels are required; unresolved roots cannot be enabled by
#' `readiness_policy = "review"`. No model is refitted. The returned EAPs
#' use the corrected point calibration and a separate normal reference prior,
#' not the original profiled Person estimates. Continuous posterior intervals
#' exclude calibration uncertainty and residual calibration bias. The correction
#' order and estimator remain recorded in the score table and saved settings.
#' [extract_mfrm_calibration()] provides the corresponding portable route.
#'
#' @section Conditional ability decisions:
#' Supply `cut_scores` on the calibration's ability scale to obtain
#' `decisions`, with one row per Person and requested event. `Type = "above_cut"`
#' gives \eqn{P(\theta > c \mid y, \widehat\psi)}; `Type = "ability_band"`
#' gives the probability between successive cuts, including the two outer
#' bands. `Lower < ability <= Upper` defines every row. `Label` uses the
#' supplied cut names (or `Cut1`, ...) and `Band1`, ... in increasing order.
#' These events overlap: only the ability-band probabilities sum to one.
#'
#' `PriorMean` and `PriorSD` record the actual Person-specific normal prior;
#' `Status` and `Reason` retain scoring availability. `available_conditional`
#' means the applicable numerical checks passed, not that the model or the
#' cut's substantive meaning has been validated. `review_only` requires
#' inspection before interpretation. Source-roster scoring may return
#' `prior_only` for missing-only Persons; portable scoring instead retains
#' them as `unavailable` with `NA` probabilities. Fitted-object scoring omits
#' Persons with no retained rows; inspect its row and population reviews.
#' A Person absent from the input is never invented.
#'
#' Probabilities integrate the continuous posterior, rather than counting
#' [sample_mfrm_plausible_values()] draws or applying a normal approximation
#' to EAP/SD. They condition on point calibration, population/prior parameters
#' and supplied cuts, as recorded in `Uncertainty`. They do not propagate
#' estimation error, compare two Persons or predict passing a future test.
#' CEFR labels require externally justified cuts on this same scale.
#' Very small tails have absolute, not relative, numerical accuracy; a
#' returned zero or one is not a claim of impossibility or certainty.
#'
#' Review the scoring population as well as the calibration. An estimated
#' normal prior is treated as known during scoring; its sampling uncertainty
#' can matter especially with a small calibration cohort. A larger calibration
#' cohort can improve parameter precision without correcting a prior that
#' mismatches a new population's mean, spread or shape. Passing numerical checks
#' does not establish that probabilities are accurate for that new population.
#' A fixed reference prior is not a universal remedy: justify it for the target
#' population on the calibration scale. If a prior override is supported,
#' comparing externally justified alternatives is a sensitivity analysis, not
#' evidence that the alternative giving the preferred classification is valid.
#' Population-averaged accuracy also does not guarantee accuracy near each cut,
#' within every subgroup or for each individual Person.
#'
#' Use [score_mfrm_persons()] for eligible ordinary RSM MML source Persons,
#' [predict_mfrm_units()] for ordinary RSM/PCM MML or JML fitted-object scoring, and
#' [score_mfrm_calibration()] for eligible continuous-posterior RSM/PCM
#' artifacts. Each route retains its existing source and population
#' restrictions. Unit weights are required; GPCM, shared-rater and testlet
#' cut probabilities are not supported. JML uses a separate scoring prior,
#' not a population estimated by JML. `NULL` leaves the existing output unchanged.
#'
#' `summary(x)$decisions` and [apa_table()] with `which = "decisions"` reuse
#' saved probabilities without rounding them. Fitted-object and portable
#' score batches also support [build_summary_table_bundle()] and
#' [export_summary_appendix()]. Save the complete result with [saveRDS()]
#' to retain cuts, priors and review states for later reporting. Existing
#' score plots display EAPs and intervals, not these event probabilities.
#'
#' @return An object of class `mfrm_unit_prediction` with components:
#' - `estimates`: posterior summaries by person
#' - `decisions`: optional cut-score and ability-band probabilities, with
#'   boundaries, labels, priors, availability and uncertainty basis
#' - `draws`: optional quadrature-grid posterior draws
#' - `fit_identity`: source calibration, fitted population coding and scoring
#'   status used to check reuse with fit-level exports
#' - `row_review`: row-level preparation review for `new_data`
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `quadrature_review`, `quadrature_overview`: optional unrounded numerical
#'   integration comparison and its compact overview
#' - `input_data`: cleaned canonical scoring rows retained from `new_data`
#' - `person_data`: cleaned or supplied person-level background data used for
#'   latent-regression scoring; `NULL` otherwise
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [score_mfrm_persons()], [score_mfrm_calibration()], [classification_table()], [apa_table()],
#'   [predict_mfrm_population()], [fit_mfrm()], [summary.mfrm_unit_prediction]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- fit_mfrm(
#'   toy[toy$Person %in% keep_people, , drop = FALSE],
#'   "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML",
#'   quad_points = 61
#' )
#' raters <- unique(toy$Rater)[1:2]
#' criteria <- unique(toy$Criterion)[1:2]
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01", "NEW02", "NEW02"),
#'   Rater = c(raters[1], raters[2], raters[1], raters[2]),
#'   Criterion = c(criteria[1], criteria[2], criteria[1], criteria[2]),
#'   Score = c(2, 3, 2, 4)
#' )
#' # Illustrative boundaries only: these labels do not establish CEFR levels.
#' cuts <- c(Basic = -0.5, Advanced = 0.8)
#' pred_units <- predict_mfrm_units(toy_fit, new_units,
#'   scoring_quad_points = 121, cut_scores = cuts)
#' summary(pred_units)$estimates[, c("Person", "Estimate", "Lower", "Upper")]
#' pred_units$decisions
#'
#' # Optional prior-sensitivity comparison performs another scoring calculation.
#' # example("predict_mfrm_units", run.dontrun = TRUE)
#' \dontrun{
#' apa_table(pred_units, which = "decisions")
#' # Compare assumptions on the same response rows, without refitting.
#' alternative <- predict_mfrm_units(
#'   toy_fit, new_units, scoring_prior = list(mean = 0.5, sd = 1.5),
#'   scoring_quad_points = 121, cut_scores = cuts
#' )
#' comparison <- dplyr::bind_rows(
#'   "Retained prior" = pred_units$estimates,
#'   "Alternative prior" = alternative$estimates, .id = "Scenario"
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(comparison, ggplot2::aes(Scenario, Estimate)) +
#'   ggplot2::geom_pointrange(ggplot2::aes(ymin = Lower, ymax = Upper)) +
#'   ggplot2::facet_wrap(~ Person) + ggplot2::theme_minimal() +
#'   ggplot2::labs(x = NULL, y = "Posterior EAP",
#'     caption = "95% posterior intervals conditional on calibration and each prior")
#' }
#' # Retain the complete result; summaries and tables do not score again.
#' saved <- tempfile(fileext = ".rds")
#' saveRDS(pred_units, saved)
#' apa_table(readRDS(saved), which = "decisions")
#' unlink(saved)
#' }
#' @section Saved-result reuse:
#' Fit-level export and replay helpers check that supplied scores describe the
#' same calibration, fitted population coding and scoring status. Recreate
#' older scores from the saved fit and original scoring inputs/settings; the
#' calibration does not need refitting. See [mfrmr_compatibility_layer].
#' @export
#' @order 1
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
                               scoring_quad_points = 31L,
                               readiness_policy = c("error", "review"),
                               n_draws = 0,
                               seed = NULL,
                               adaptive_quad_points = NULL,
                               scoring_prior = NULL,
                               cut_scores = NULL) {
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit")) {
    stop("`fit` must be a native fit_mfrm() result; use the source package to score an imported model.", call. = FALSE)
  }
  # Retain the caller's calibration before temporary corrected-JML expansion.
  fit_identity <- prediction_fit_identity(fit)
  cut_scores <- mfrm_validate_cut_scores(cut_scores)
  adjusted <- if (mfrm_has_jml_adjustment(fit)) mfrm_jml_scoring_components(fit) else NULL
  if (!is.null(adjusted)) fit$config <- adjusted$config
  product <- mfrm_has_product_slopes(fit)
  if (product && (!is.null(scoring_prior) || !is.null(person_data) || !is.null(person_id))) {
    stop("Two-family scoring retains the fitted N(0,1) prior; scoring_prior and person-level population inputs are not supported.", call. = FALSE)
  }
  scoring_prior <- prediction_validate_scoring_prior(scoring_prior)
  adaptive_quad_points <- mfrmr_validate_adaptive_quad_points(adaptive_quad_points)
  if (!is.null(cut_scores) && (!fit$config$model %in% c("RSM","PCM") || mfrm_has_jml_adjustment(fit)))
    stop("Cut-score probabilities currently require ordinary RSM/PCM scoring.", call. = FALSE)
  if (!is.null(scoring_prior) && length(setdiff(
      fit$population$design_columns %||% names(fit$population$coefficients), "(Intercept)"))) {
    stop("`scoring_prior` cannot override a population model with background covariates; use the fitted person_data contract without a prior override.", call. = FALSE)
  }
  fit_method <- prediction_resolve_fit_method(fit)
  if (!fit_method %in% c("MML", "JML")) {
    stop("`predict_mfrm_units()` currently supports only fits estimated with method = 'MML' or 'JML'.",
         call. = FALSE)
  }
  readiness_policy <- match.arg(readiness_policy)
  source_readiness <- if (is.null(adjusted)) prediction_source_scoring_readiness(fit) else adjusted$source
  if (identical(source_readiness$status, "unavailable")) {
    stop("The retained calibration or scoring prior is invalid; review cannot enable scoring.", call. = FALSE)
  }
  if (!isTRUE(source_readiness$ready) &&
      identical(readiness_policy, "error")) {
    stop(
      "`fit` is not ready for fitted-object scoring. ",
      paste(prediction_scoring_reasons(source_readiness$reason_codes), collapse = " "),
      if (!is.null(source_readiness$local_calibration_review))
        paste0(" ", source_readiness$local_calibration_review$review),
      " Resolve the source fit or use `readiness_policy = \"review\"` ",
      "for an explicitly review-only calculation.",
      call. = FALSE
    )
  }
  if (!is.numeric(interval_level) || length(interval_level) != 1L ||
      !is.finite(interval_level) || interval_level <= 0 || interval_level >= 1) {
    stop("`interval_level` must be a single number in (0, 1).", call. = FALSE)
  }
  n_draws <- prediction_validate_integer(n_draws[1] %||% 0L, "n_draws", min_value = 0L, positive = FALSE)
  scoring_quad_points <- prediction_validate_integer(
    scoring_quad_points[1] %||% 31L,
    "scoring_quad_points",
    positive = TRUE
  )
  if (scoring_quad_points < 2L) {
    stop("`scoring_quad_points` must be an integer greater than or equal to 2.",
         call. = FALSE)
  }

  higher <- max(scoring_quad_points + 10L, 2L * scoring_quad_points - 1L)
  if (inherits(tryCatch(gauss_hermite_normal(higher), error = identity), "error")) {
    higher <- scoring_quad_points + 10L
  }
  adaptive_quad_points <- mfrmr_validate_adaptive_quad_points(
    unique(c(adaptive_quad_points, max(3L, scoring_quad_points), higher)))

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
  scoring_population <- population_ready$spec
  labels <- prepared$prep$levels$Person
  retained_basis <- if (isTRUE(population_ready$active)) "population_model" else
    as.character(fit$config$posterior_basis %||% "legacy_mml")
  retained <- resolve_person_quadrature_basis(gauss_hermite_normal(2L),
    scoring_population, person_count = length(labels))
  prior_comparison <- data.frame(Person = as.character(labels),
    RetainedPriorMean = if (isTRUE(retained$transformed)) retained$mu else 0,
    RetainedPriorSD = if (isTRUE(retained$transformed)) retained$sigma else 1)
  if (!is.null(scoring_prior)) {
    if (isTRUE(population_ready$active) &&
        (ncol(scoring_population$design_matrix) != 1L ||
          !all(scoring_population$design_matrix == 1))) {
      stop("`scoring_prior` cannot override a population model with background covariates; use the fitted person_data contract without a prior override.", call. = FALSE)
    }
    scoring_population <- list(active = TRUE,
      design_matrix = matrix(1, length(labels), 1L),
      person_lookup = seq_along(labels),
      coefficients = c(`(Intercept)` = scoring_prior$mean),
      sigma2 = scoring_prior$sd^2, posterior_basis = "user_supplied_normal")
    if (!is.null(population_ready$population_review)) {
      population_ready$population_review$PosteriorBasis <- "user_supplied_normal"
    }
  }

  idx <- build_indices(
    prepared$prep,
    step_facet = fit$config$step_facet,
    slope_facet = if (product) NULL else fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  params <- if (is.null(adjusted)) {
    sizes <- build_param_sizes(fit$config)
    expand_params(fit$opt$par, sizes, fit$config)
  } else adjusted$params
  if (product) {
    if (any(prepared$prep$data$Weight != 1))
      stop("Two-family scoring requires unit observation weights.", call. = FALSE)
    if (anyDuplicated(prepared$prep$data[c("Person", fit$config$slope_facet)]))
      stop("Repeated person-facet combinations require a separate observation model.", call. = FALSE)
    # Assemble from the fitted components, including previously unobserved
    # combinations of known levels. Never index the fitted crossings by batch order.
    params$slopes <- Reduce(`*`, lapply(fit$config$slope_facet, function(owner)
      params$slope_components[[owner]]$Slope[idx$facets[[owner]]]))
    if (any(!is.finite(params$slopes) | params$slopes <= 0))
      stop("Predicted slope products must be finite and positive.", call. = FALSE)
    idx$slope_idx <- seq_along(idx$score_k)
  }
  if (!is.null(adjusted) && any(prepared$prep$data$Weight != 1))
    stop("Corrected-JML scoring requires unit observation weights.",call.=FALSE)
  quad_points <- as.integer(scoring_quad_points)
  quad <- gauss_hermite_normal(quad_points)

  if (!is.null(cut_scores) && any((idx$weight %||% 1) != 1))
    stop("Cut-score probabilities require unit observation weights.", call. = FALSE)
  scored <- compute_person_posterior_summary(
    idx = idx,
    config = fit$config,
    params = params,
    quad = quad,
    person_labels = prepared$prep$levels$Person,
    population_spec = scoring_population,
    interval_level = interval_level,
    n_draws = n_draws,
    seed = seed,
    cut_scores = cut_scores
  )
  values <- as.matrix(scored$estimates[c("Estimate", "SD", "Lower", "Upper", "PriorMean", "PriorSD")])
  if (any(!is.finite(values)) || any(scored$estimates$SD < 0) ||
      any(scored$estimates$PriorSD <= 0)) {
    stop("Scoring could not produce finite posterior summaries under the supplied calibration and prior.", call. = FALSE)
  }
  prior_comparison <- prior_comparison[match(scored$estimates$Person, prior_comparison$Person), , drop = FALSE]
  prior_comparison$PriorMean <- scored$estimates$PriorMean
  prior_comparison$PriorSD <- scored$estimates$PriorSD
  scored$estimates$RetainedPriorMean <- prior_comparison$RetainedPriorMean
  scored$estimates$RetainedPriorSD <- prior_comparison$RetainedPriorSD
  scored$estimates$SourceScoringReady <- isTRUE(source_readiness$ready)
  scored$estimates$EstimateUse <- if (isTRUE(source_readiness$ready)) {
    "fitted_object_scoring"
  } else {
    "review_only_nonready_source"
  }

  calibration_note <- if (!is.null(scoring_prior)) {
    "Posterior summaries use the unchanged fitted calibration and the explicitly supplied normal scoring prior. The original population assumptions remain recorded separately."
  } else if (isTRUE(population_ready$active)) {
    "Posterior summaries are computed under the fitted MML calibration together with the fitted conditional normal population model for the scored persons."
  } else if (identical(fit_method, "MML")) {
    "Posterior summaries are computed under the fixed fitted MML calibration with a standard normal N(0,1) scoring prior."
  } else {
    "Posterior summaries are computed under the fixed fitted JML calibration using a standard normal reference prior on the quadrature grid."
  }

  notes <- c(
    calibration_note,
    "Intervals invert the continuous posterior CDF and condition on the point calibration and scoring prior; uncertainty from estimating calibration or population parameters is excluded.",
    "Non-person facets in `new_data` must already exist in the fitted calibration.",
    "Overlapping person IDs are treated as labels in `new_data`; the original fitted person estimates are not updated."
  )
  if (identical(fit_method, "JML")) {
    notes <- c(
      notes,
      "For JML fits, the returned person scores are post hoc fitted-object EAP summaries rather than direct JML estimates from the original fit."
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
  if (any(scored$estimates$WeightedLikelihood)) {
    notes <- c(notes, "Non-unit response weights exponentiate likelihood contributions. Posterior uncertainty describes that weighting rule; weights do not establish equivalent independent ratings or frequentist coverage.")
  }

  out <- structure(
    list(
      estimates = scored$estimates,
      fit_identity = fit_identity,
      draws = scored$draws,
      row_review = prepared$row_review,
      population_review = population_ready$population_review,
      input_data = prepared$input_data,
      person_data = population_ready$input_data,
      settings = list(
        scoring_prior = scoring_prior,
        retained_posterior_basis = retained_basis,
        prior_comparison = prior_comparison,
        interval_level = interval_level,
        n_draws = n_draws,
        quad_points = as.integer(quad_points),
        scoring_quad_points = as.integer(quad_points),
        score_integration_required = TRUE,
        scoring_algorithm = if (mfrmr_adaptive_integration(fit$config)) {
          "adaptive_quadrature_eap_v2"
        } else "quadrature_eap_v2",
        readiness_policy = readiness_policy,
        source_scoring_ready = isTRUE(source_readiness$ready),
        source_repeated_ratings = fit$config$repeated_ratings %||% "review",
        source_scoring_status = source_readiness$status,
        source_scoring_policy_basis = source_readiness$policy_basis,
        local_calibration_review = source_readiness$local_calibration_review,
        source_audit_states = source_readiness$source_audit_states,
        source_scoring_reason_codes = source_readiness$reason_codes,
        source_fit_readiness = source_readiness$fit_readiness,
        source_inference_ready = source_readiness$inference_ready,
        source_readiness_contract_version =
          source_readiness$readiness_contract_version,
        seed = seed,
        method = fit_method,
        source_columns = prepared$prep$source_columns,
        posterior_basis = if (!is.null(scoring_prior)) {
          "user_supplied_normal"
        } else if (isTRUE(population_ready$active)) {
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
  out$quadrature_review <- mfrmr_adaptive_quadrature_review(
    idx, fit$config, params, quad, prepared$prep$levels$Person,
    adaptive_quad_points, population_spec = scoring_population
  )
  out$settings$adaptive_quad_points <- adaptive_quad_points
  if (is.null(scoring_prior) && !is.null(population_ready$spec$fixed_sd))
    out$settings$population_fixed_sd <- population_ready$spec$fixed_sd
  out$notes <- c(out$notes,
    "Scoring integration compares the reported EAP and SD with adaptive reference orders under the same calibration and prior. Passing does not validate the scoring prior for another population.")
  out$settings$score_integration_review <- prediction_score_integration_review(
    out$estimates, out$quadrature_review)
  passed <- out$settings$score_integration_review$Passed
  out$estimates$ScoreIntegrationReady <- passed
  if (!all(passed)) {
    if (identical(readiness_policy, "error")) {
      stop("Posterior scoring integration did not pass the numerical accuracy check. Increase `scoring_quad_points` or use `readiness_policy = \"review\"` to inspect the flagged scores.", call. = FALSE)
    }
    if (isTRUE(source_readiness$ready)) {
      out$estimates$EstimateUse[!passed] <- "review_only_scoring_integration"
    }
  }
  if (!is.null(adjusted)) out$notes <- c(out$notes,mfrm_jml_scoring_note(adjusted$source))
  if (product) {
    out$settings$two_family_calibration <- prediction_gmfrm_calibration(fit, params)
    out$notes <- c(out$notes,
      "Experimental two-family EAP uses both fitted slope components and the original step coding on the fixed N(0,1) scale. Known levels may form new combinations; their adequacy is a model assumption, not an observed calibration result.",
      "These conditional scores do not establish global identification, boundary absence, sampling coverage or validity for another population. extract_mfrm_calibration() provides the separately validated portable route.")
  }
  out$estimates <- prediction_estimate_table(out)
  out$draws <- prediction_draw_table(out$draws, out$estimates)
  out$notes <- prediction_output_notes(out)
  if (!is.null(cut_scores)) {
    out$settings$cut_scores <- cut_scores
    out$decisions <- mfrm_decision_table(scored$decision_probabilities,
      out$estimates$Person, cut_scores, out$estimates$PriorMean, out$estimates$PriorSD,
      ifelse(out$estimates$SourceScoringReady & out$estimates$ScoreIntegrationReady,
        "available_conditional", "review_only"),
      ifelse(out$estimates$SourceScoringReady & out$estimates$ScoreIntegrationReady,
        "", out$estimates$EstimateUse))
    out$notes <- c(out$notes, mfrm_decision_note())
  }
  out
}

#' @rdname predict_mfrm_units
#' @details `summary()` method: Summarize posterior unit scoring output
#' @param object Output from [predict_mfrm_units()].
#' @param digits Number of digits used in numeric summaries. Decision
#'   probabilities, cut boundaries and scoring priors are not rounded.
#' @param ... Reserved for generic compatibility.
#' @return For `summary()`: An object of class `summary.mfrm_unit_prediction` with:
#' - `estimates`: posterior summaries by person
#' - `decisions`: unrounded conditional probabilities, when requested
#' - `row_review`: row-preparation review
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_units()] for an executable scoring and reporting
#'   example; [apa_table()], [export_summary_appendix()]
#' @method summary mfrm_unit_prediction
#' @export
summary.mfrm_unit_prediction <- function(object, digits = 3, ...) {
  mfrm_validate_decisions(object)
  if (!inherits(object, "mfrm_unit_prediction")) {
    stop("`object` must be output from predict_mfrm_units().", call. = FALSE)
  }
  prediction_validate_population_output(object)
  digits <- prediction_validate_integer(digits, "digits", min_value = 0L, positive = FALSE)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1)) & !names(df) %in% c("IntervalLevel", "PriorMean", "PriorSD", "RetainedPriorMean", "RetainedPriorSD")
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  out <- list(
    estimates = round_df(prediction_estimate_table(object)),
    row_review = round_df(object$row_review),
    population_review = round_df(object$population_review),
    settings = object$settings,
    notes = prediction_output_notes(object),
    digits = digits
  )
  if (!is.null(object$quadrature_review)) {
    out$quadrature_review <- object$quadrature_review
    out$quadrature_overview <- mfrmr_adaptive_quadrature_overview(object$quadrature_review)
  }
  class(out) <- "summary.mfrm_unit_prediction"
  if (!is.null(object$decisions)) out$decisions <- object$decisions
  out
}

#' @export
print.summary.mfrm_unit_prediction <- function(x, ...) {
  prediction_print_summary(x)
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
#'   ordinary fitted-object posterior scoring. Intercept-only latent-regression fits
#'   can reconstruct the minimal scored-person table internally. This is the
#'   scoring-time table for `new_data`, not the fit object's replay/export
#'   provenance table. For categorical background variables, supply values on
#'   the same coding scale used at fit time; the fitted factor levels and
#'   contrasts are reused when building the scoring design matrix.
#'   Fitted standardization and polynomial/spline bases are also reused;
#'   they are not re-estimated from the scoring cohort.
#' @param person_id Optional person-ID column in `person_data`.
#' @param population_policy How missing background data are handled when
#'   `fit` uses the latent-regression branch. `"error"` (default) requires
#'   complete person-level covariates; `"omit"` drops scored persons lacking
#'   complete covariates and records that omission in `population_review`.
#' @param n_draws Number of posterior draws per person. Must be a positive
#'   integer.
#' @param interval_level Posterior interval level passed to
#'   [predict_mfrm_units()] for the accompanying EAP summary table. The same
#'   level selects the empirical draw quantiles reported by `summary()`.
#' @param scoring_quad_points Number of Gauss-Hermite nodes used only for this
#'   scoring call. Passed to [predict_mfrm_units()] and independent of the
#'   fit-time quadrature order. Fixed or adaptive integration is inherited
#'   from `fit`.
#' @param readiness_policy Source-fit readiness policy passed to
#'   [predict_mfrm_units()]. Conditional calibration and scoring-integration
#'   checks are shared with that function and retained in the returned settings.
#' @param seed Optional seed for reproducible posterior draws.
#' @param scoring_prior Optional normal scoring prior, `list(mean = ..., sd = ...)`,
#'   passed to [predict_mfrm_units()]. `NULL` retains the fitted scoring prior.
#'   Estimates and draws retain both the original and the supplied prior;
#'   checks of the fitted calibration and numerical accuracy are unchanged.
#'
#' @details
#' `sample_mfrm_plausible_values()` uses [predict_mfrm_units()] to return
#' draws from each Person's posterior distribution as
#' a standalone object. It is useful when downstream workflows want repeated
#' latent-value imputations rather than just one posterior EAP summary.
#'
#' In the current `mfrmr` implementation these are **approximate plausible
#' values** drawn from the fitted quadrature-grid posterior under the scoring
#' basis implied by `fit`, unless `scoring_prior` explicitly supplies another
#' normal prior. With the default prior, for ordinary `MML` fits this is the fitted marginal
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
#'   `LowerValue` and `UpperValue` use the requested `interval_level` and are
#'   labelled with `IntervalLevel` and `DrawSummaryBasis`. They are empirical
#'   quantiles of a finite sample of discrete draws, not the continuous posterior
#'   limits in the companion estimates table. With few draws they are coarse;
#'   one draw has an unavailable empirical SD.
#'   Recreate older summaries from the original plausible-values object to
#'   update these limits and labels; no refitting or resampling is needed.
#'
#' @section What this does not justify:
#' This helper does not update the calibration, estimate new non-person facet
#' levels, or provide exact future true values. It samples from the quadrature-
#' grid posterior implied by the existing fitted-model scoring basis, using
#' fixed or Person-specific adaptive nodes according to the fit's setting.
#' Calibration and population parameters remain fixed. These draws alone do
#' not validate downstream group comparisons or regressions; those analyses
#' require a compatible conditioning model and sampling design.
#'
#' @section References:
#' The underlying posterior scoring follows the usual quadrature-based EAP
#' framework of Bock and Aitkin (1981). The interpretation of multiple
#' posterior draws as plausible-value-style summaries follows the general logic
#' discussed by Mislevy (1991), while the current implementation remains a
#' practical fitted-object posterior approximation rather than a full published
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
#' - `fit_identity`: source calibration, fitted population coding and scoring
#'   status retained from [predict_mfrm_units()]
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
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("sample_mfrm_plausible_values", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- fit_mfrm(
#'   toy[toy$Person %in% keep_people, , drop = FALSE],
#'   "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML",
#'   quad_points = 5,
#'   maxit = 30
#' )
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01"),
#'   Rater = unique(toy$Rater)[1],
#'   Criterion = unique(toy$Criterion)[1:2],
#'   Score = c(2, 3)
#' )
#' pv <- sample_mfrm_plausible_values(toy_fit, new_units, n_draws = 3, seed = 1)
#' summary(pv)$draw_summary
#' }
#' @section Saved-result reuse:
#' Fit-level export and replay helpers check the source calibration and scoring
#' status of supplied draws. Recreate older results from the saved fit and
#' original scoring inputs/settings, including the seed. This repeats drawing
#' values without refitting the calibration. See [mfrmr_compatibility_layer].
#' @export
#' @order 1
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
                                         scoring_quad_points = 31L,
                                         readiness_policy = c("error", "review"),
                                         seed = NULL,
                                         scoring_prior = NULL) {
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
    scoring_quad_points = scoring_quad_points,
    readiness_policy = readiness_policy,
    n_draws = n_draws,
    seed = seed,
    scoring_prior = scoring_prior
  )

  fit_method <- prediction_resolve_fit_method(fit)
  draw_note <- if (identical(pred$settings$posterior_basis, "user_supplied_normal")) {
    "These draws use the user-supplied normal scoring prior and unchanged fitted calibration; the prior is an analyst assumption."
  } else if (identical(pred$settings$posterior_basis %||% "", "population_model")) {
    "These draws are sampled from the fitted population-model posterior implied by the latent-regression MML branch."
  } else if (identical(fit_method, "MML")) {
    "These draws are sampled from the quadrature-grid posterior under the existing MML calibration and its fixed or adaptive integration setting."
  } else {
    "These draws are sampled from a fitted-object quadrature-grid posterior built from the fitted JML parameters and a standard normal reference prior."
  }

  notes <- c(
    draw_note,
    "Use them as approximate plausible-value summaries for posterior uncertainty, not as deterministic future truth values.",
    "Draws alone do not validate downstream group comparisons or regressions; check the conditioning model and sampling design for the intended analysis.",
    pred$notes
  )

  structure(
    list(
      values = pred$draws,
      fit_identity = pred$fit_identity,
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

#' @rdname sample_mfrm_plausible_values
#' @details `summary()` method: Summarize approximate plausible values from posterior scoring
#' @param object Output from [sample_mfrm_plausible_values()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#' @return For `summary()`: An object of class `summary.mfrm_plausible_values` with:
#' - `draw_summary`: empirical summaries of the sampled values by person
#' - `estimates`: companion posterior EAP summaries
#' - `row_review`: row-preparation review
#' - `population_review`: optional person-level omission review for
#'   latent-regression scoring
#' - `settings`: scoring settings
#' - `notes`: interpretation notes
#' @seealso [sample_mfrm_plausible_values()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("summary.mfrm_plausible_values", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy$Person)[1:18]
#' toy_fit <- fit_mfrm(
#'   toy[toy$Person %in% keep_people, , drop = FALSE],
#'   "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML",
#'   quad_points = 5,
#'   maxit = 30
#' )
#' new_units <- data.frame(
#'   Person = c("NEW01", "NEW01"),
#'   Rater = unique(toy$Rater)[1],
#'   Criterion = unique(toy$Criterion)[1:2],
#'   Score = c(2, 3)
#' )
#' pv <- sample_mfrm_plausible_values(toy_fit, new_units, n_draws = 3, seed = 1)
#' summary(pv)
#' }
#' @method summary mfrm_plausible_values
#' @export
summary.mfrm_plausible_values <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_plausible_values")) {
    stop("`object` must be output from sample_mfrm_plausible_values().", call. = FALSE)
  }
  prediction_validate_population_output(object)
  digits <- prediction_validate_integer(digits, "digits", min_value = 0L, positive = FALSE)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    num_cols <- vapply(df, is.numeric, logical(1)) & !names(df) %in% c("IntervalLevel", "PriorMean", "PriorSD", "RetainedPriorMean", "RetainedPriorSD")
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  level <- object$settings$interval_level %||% NA_real_
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("The requested interval level is unavailable; recreate the scoring result with an explicit interval_level before summarizing draw quantiles.", call. = FALSE)
  }
  alpha <- (1 - level) / 2
  draw_summary <- object$values |>
    dplyr::group_by(.data$Person) |>
    dplyr::summarise(
      Draws = dplyr::n(),
      MeanValue = mean(.data$Value),
      SDValue = stats::sd(.data$Value),
      LowerValue = stats::quantile(.data$Value, probs = alpha, names = FALSE, type = 1),
      UpperValue = stats::quantile(.data$Value, probs = 1 - alpha, names = FALSE, type = 1),
      .groups = "drop"
    )

  draw_summary$IntervalLevel <- level
  draw_summary$DrawSummaryBasis <- "Empirical quantiles of discrete posterior draws; not continuous posterior interval bounds"
  draw_summary <- prediction_draw_table(draw_summary, prediction_estimate_table(object))

  out <- list(
    draw_summary = round_df(draw_summary),
    estimates = round_df(prediction_estimate_table(object)),
    row_review = round_df(object$row_review),
    population_review = round_df(object$population_review),
    settings = object$settings,
    notes = prediction_output_notes(object),
    digits = digits
  )
  class(out) <- "summary.mfrm_plausible_values"
  out
}

#' @export
print.summary.mfrm_plausible_values <- function(x, ...) {
  prediction_print_summary(x, plausible = TRUE)
}
