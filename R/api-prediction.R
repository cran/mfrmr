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
#' @section Interpreting output:
#' - `forecast` contains facet-level expected summaries for the requested
#'   future design.
#' - `Mcse*` columns quantify Monte Carlo uncertainty from using a finite number
#'   of replications.
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
#' - `settings`: forecasting settings
#' - `ademp`: simulation-study metadata
#' - `notes`: interpretation notes
#' @seealso [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()],
#'   [evaluate_mfrm_design()], [summary.mfrm_population_prediction]
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   n_person = 40,
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   assignment = "rotating"
#' )
#' pred <- predict_mfrm_population(
#'   sim_spec = spec,
#'   n_person = 60,
#'   reps = 2,
#'   maxit = 10,
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
    base_spec <- extract_mfrm_sim_spec(fit)
    default_fit_method <- as.character(fit$summary$Method[1])
    if (!is.character(default_fit_method) || !nzchar(default_fit_method)) {
      default_fit_method <- "MML"
    }
    if (identical(default_fit_method, "JMLE")) default_fit_method <- "JML"
    default_model <- as.character(fit$summary$Model[1])
    if (!is.character(default_model) || !nzchar(default_model)) {
      default_model <- as.character(fit$config$model)
    }
  } else {
    if (!inherits(sim_spec, "mfrm_sim_spec")) {
      stop("`sim_spec` must be output from build_mfrm_sim_spec() or extract_mfrm_sim_spec().", call. = FALSE)
    }
    base_spec <- sim_spec
    default_fit_method <- "MML"
    default_model <- as.character(base_spec$model)
  }

  fit_method <- toupper(as.character(fit_method[1] %||% default_fit_method))
  fit_method <- match.arg(fit_method, c("JML", "MML"))
  model <- toupper(as.character(model[1] %||% default_model))
  model <- match.arg(model, c("RSM", "PCM"))

  design <- tibble::tibble(
    n_person = simulation_validate_count(if (is.null(n_person)) base_spec$n_person else n_person, "n_person", min_value = 2L),
    n_rater = simulation_validate_count(if (is.null(n_rater)) base_spec$n_rater else n_rater, "n_rater", min_value = 2L),
    n_criterion = simulation_validate_count(if (is.null(n_criterion)) base_spec$n_criterion else n_criterion, "n_criterion", min_value = 2L),
    raters_per_person = simulation_validate_count(if (is.null(raters_per_person)) base_spec$raters_per_person else raters_per_person, "raters_per_person", min_value = 1L)
  )
  if (design$raters_per_person > design$n_rater) {
    stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  }

  forecast_spec <- simulation_override_spec_design(
    base_spec,
    n_person = design$n_person,
    n_rater = design$n_rater,
    n_criterion = design$n_criterion,
    raters_per_person = design$raters_per_person
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
  notes <- c(
    "This forecast summarizes expected design-level behavior under the supplied or fit-derived simulation specification.",
    "MCSE columns quantify Monte Carlo uncertainty from using a finite number of replications.",
    "Do not interpret this output as deterministic future person/rater true values."
  )

  structure(
    list(
      design = design,
      forecast = tibble::as_tibble(sim_summary$design_summary),
      overview = tibble::as_tibble(sim_summary$overview),
      simulation = sim_eval,
      sim_spec = forecast_spec,
      settings = list(
        reps = as.integer(reps[1]),
        fit_method = fit_method,
        model = model,
        maxit = maxit,
        quad_points = quad_points,
        residual_pca = residual_pca,
        source = if (has_fit) "fit_mfrm" else "mfrm_sim_spec",
        seed = seed
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
#' - `ademp`: simulation-study metadata
#' - `notes`: interpretation notes
#' @seealso [predict_mfrm_population()]
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   n_person = 40,
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   assignment = "rotating"
#' )
#' pred <- predict_mfrm_population(
#'   sim_spec = spec,
#'   n_person = 60,
#'   reps = 2,
#'   maxit = 10,
#'   seed = 123
#' )
#' summary(pred)
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
    ademp = object$ademp %||% NULL,
    notes = object$notes %||% character(0)
  )
  class(out) <- "summary.mfrm_population_prediction"
  out
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

  audit <- tibble::tibble(
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

  list(prep = prep, audit = audit, input_data = input_data)
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
  log_w <- log(quad$weights)
  n_nodes <- length(quad$nodes)
  score_k <- idx$score_k

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    k_cat <- length(step_cum)
    step_cum_row <- matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
    obs_idx <- cbind(seq_len(n), score_k + 1L)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      eta_mat <- outer(eta_q, 0:(k_cat - 1))
      log_num <- eta_mat - step_cum_row
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
      eta_q <- base_eta + quad$nodes[q]
      log_num <- outer(eta_q, k_vals) - step_cum_obs
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  }

  ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)
  n_persons <- nrow(ll_by_person)
  log_post <- matrix(log_w, nrow = n_persons, ncol = n_nodes, byrow = TRUE) + ll_by_person
  row_max <- log_post[cbind(seq_len(n_persons), max.col(log_post))]
  log_norm <- row_max + log(rowSums(exp(log_post - row_max)))
  post_w <- exp(log_post - log_norm)

  nodes_mat <- matrix(quad$nodes, nrow = n_persons, ncol = n_nodes, byrow = TRUE)
  eap <- rowSums(nodes_mat * post_w)
  sd_eap <- sqrt(rowSums((nodes_mat - eap)^2 * post_w))
  alpha <- max(min((1 - interval_level) / 2, 0.5), 0)
  lower <- apply(post_w, 1, posterior_quantile_on_grid, nodes = quad$nodes, p = alpha)
  upper <- apply(post_w, 1, posterior_quantile_on_grid, nodes = quad$nodes, p = 1 - alpha)

  obs_n <- as.integer(rowsum(rep(1L, n), person_int, reorder = FALSE)[, 1])
  weight_n <- as.numeric(rowsum(if (is.null(idx$weight)) rep(1, n) else idx$weight,
                                person_int, reorder = FALSE)[, 1])

  estimates <- tibble::tibble(
    Person = as.character(person_labels),
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
          Person = as.character(person_labels[i]),
          Draw = seq_len(n_draws),
          Value = sample(quad$nodes, size = n_draws, replace = TRUE, prob = post_w[i, ])
        )
      })
      dplyr::bind_rows(draw_list)
    })
  }

  list(estimates = estimates, draws = draws_tbl)
}

#' Score future or partially observed units under a fixed MML calibration
#'
#' @param fit Output from [fit_mfrm()] with `method = "MML"`.
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
#' @param interval_level Posterior interval level returned in `Lower`/`Upper`.
#' @param n_draws Optional number of quadrature-grid posterior draws to return
#'   per scored person. Use 0 to skip draws.
#' @param seed Optional seed for reproducible posterior draws.
#'
#' @details
#' `predict_mfrm_units()` is the **individual-unit companion** to
#' [predict_mfrm_population()]. It uses a fixed MML calibration and scores new
#' or partially observed persons via Expected A Posteriori (EAP) summaries on
#' the fitted quadrature grid.
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
#' values under the fixed calibration. They should be interpreted as posterior
#' uncertainty summaries, not as deterministic future truth values.
#'
#' @section Interpreting output:
#' - `estimates` contains posterior EAP summaries for each person in
#'   `new_data`.
#' - `Lower` and `Upper` are quadrature-grid posterior interval bounds at the
#'   requested `interval_level`.
#' - `SD` is posterior uncertainty under the fixed MML calibration.
#' - `draws`, when requested, contains approximate plausible values on the
#'   fitted quadrature grid.
#'
#' @section What this does not justify:
#' This helper does not update the original calibration, estimate new non-person
#' facet levels, or produce deterministic future person true values. It scores
#' new response patterns under a fixed calibration model.
#'
#' @section References:
#' The posterior summaries follow the usual MML/EAP scoring framework used in
#' item response modeling under fixed calibrated parameters (for example Bock &
#' Aitkin, 1981). Optional posterior draws are exposed as quadrature-grid
#' plausible-value-style summaries in the spirit of Mislevy (1991), but here
#' they are offered as practical uncertainty summaries for fixed-calibration
#' many-facet scoring rather than as a direct implementation of a published
#' many-facet plausible-values procedure.
#'
#' - Bock, R. D., & Aitkin, M. (1981). *Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm*. Psychometrika, 46(4),
#'   443-459.
#' - Mislevy, R. J. (1991). *Randomization-based inference about latent
#'   variables from complex samples*. Psychometrika, 56(2), 177-196.
#'
#' @return An object of class `mfrm_unit_prediction` with components:
#' - `estimates`: posterior summaries by person
#' - `draws`: optional quadrature-grid posterior draws
#' - `audit`: row-level preparation audit for `new_data`
#' - `input_data`: cleaned canonical scoring rows retained from `new_data`
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
#'     maxit = 15
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
                               interval_level = 0.95,
                               n_draws = 0,
                               seed = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be output from fit_mfrm().", call. = FALSE)
  }
  if (!identical(as.character(fit$config$method), "MML")) {
    stop("`predict_mfrm_units()` currently supports only fits estimated with method = 'MML'.",
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

  idx <- build_indices(prepared$prep, step_facet = fit$config$step_facet)
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
    interval_level = interval_level,
    n_draws = n_draws,
    seed = seed
  )

  notes <- c(
    "Posterior summaries are computed under the fixed fitted MML calibration.",
    "Non-person facets in `new_data` must already exist in the fitted calibration.",
    "Overlapping person IDs are treated as labels in `new_data`; the original fitted person estimates are not updated."
  )
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
      audit = prepared$audit,
      input_data = prepared$input_data,
      settings = list(
        interval_level = interval_level,
        n_draws = n_draws,
        quad_points = as.integer(quad_points),
        seed = seed,
        method = "MML",
        source_columns = prepared$prep$source_columns
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
#' - `audit`: row-preparation audit
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
#'     maxit = 15
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
    audit = round_df(object$audit),
    settings = object$settings,
    notes = object$notes %||% character(0)
  )
  class(out) <- "summary.mfrm_unit_prediction"
  out
}

#' Sample approximate plausible values under a fixed MML calibration
#'
#' @param fit Output from [fit_mfrm()] with `method = "MML"`.
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
#' values** drawn from the fitted quadrature-grid posterior under a fixed MML
#' calibration. They should be interpreted as posterior uncertainty summaries
#' for the scored persons, not as deterministic future truth values and not as
#' a full many-facet plausible-values procedure with additional population
#' modeling.
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
#' posterior implied by the existing MML fit.
#'
#' @section References:
#' The underlying posterior scoring follows the usual MML/EAP framework of
#' Bock and Aitkin (1981). The interpretation of multiple posterior draws as
#' plausible-value-style summaries follows the general logic discussed by
#' Mislevy (1991), while the current implementation remains a practical fixed-
#' calibration approximation rather than a full published many-facet plausible-
#' values method.
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
#' - `audit`: row-preparation audit
#' - `input_data`: cleaned canonical scoring rows retained from `new_data`
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
#'     maxit = 15
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
    interval_level = interval_level,
    n_draws = n_draws,
    seed = seed
  )

  notes <- c(
    "These draws are sampled from the fixed quadrature-grid posterior under the existing MML calibration.",
    "Use them as approximate plausible-value summaries for posterior uncertainty, not as deterministic future truth values."
  )

  structure(
    list(
      values = pred$draws,
      estimates = pred$estimates,
      audit = pred$audit,
      input_data = pred$input_data,
      settings = pred$settings,
      notes = notes
    ),
    class = "mfrm_plausible_values"
  )
}

#' Summarize approximate plausible values from fixed-calibration scoring
#'
#' @param object Output from [sample_mfrm_plausible_values()].
#' @param digits Number of digits used in numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @return An object of class `summary.mfrm_plausible_values` with:
#' - `draw_summary`: empirical summaries of the sampled values by person
#' - `estimates`: companion posterior EAP summaries
#' - `audit`: row-preparation audit
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
#'     maxit = 15
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
    audit = round_df(object$audit),
    settings = object$settings,
    notes = object$notes %||% character(0)
  )
  class(out) <- "summary.mfrm_plausible_values"
  out
}
