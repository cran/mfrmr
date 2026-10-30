# Shared adaptive integration and explicit fixed-parameter integration review.

mfrmr_adaptive_integration <- function(config) {
  mode <- config$estimation_control$mml_integration %||% "fixed"
  if (!is.character(mode) || length(mode) != 1L || !mode %in% c("fixed", "adaptive")) {
    stop("Unrecognized MML integration mode in the retained configuration.", call. = FALSE)
  }
  identical(mode, "adaptive")
}

new_adaptive_quadrature_numeric_error <- function(message) {
  structure(list(message = message, call = NULL),
    class = c("mfrmr_adaptive_quadrature_numeric_error", "error", "condition"))
}

mfrmr_adaptive_posterior_location <- function(evaluate) {
  mode <- stats::uniroot(function(z) evaluate(z)$score, c(-1, 1),
                         extendInt = "downX", check.conv = TRUE, tol = 1e-12)$root
  at_mode <- evaluate(mode)
  scale <- 1 / sqrt(at_mode$information)
  if (!is.finite(scale) || scale <= 0) {
    stop(new_adaptive_quadrature_numeric_error("Invalid local posterior scale."))
  }
  list(mode = mode, scale = scale, score = at_mode$score)
}

mfrmr_adaptive_quadrature_basis <- function(idx, config, params, quad, base_eta,
                                            population_spec = NULL,
                                            person_count = config$n_person) {
  basis <- resolve_person_quadrature_basis(
    quad, population_spec %||% materialize_population_spec(config, params),
    person_count = person_count
  )
  mu <- if (isTRUE(basis$transformed)) basis$mu else rep(0, person_count)
  sigma <- if (isTRUE(basis$transformed)) basis$sigma else 1
  weights <- idx$weight %||% rep(1, length(idx$score_k))
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Adaptive quadrature requires finite non-negative observation weights.", call. = FALSE)
  }
  persons <- sort(unique(idx$person))
  groups <- split(seq_along(idx$person), factor(idx$person, levels = persons))
  for (i in seq_along(groups)) {
    person <- persons[i]
    rows <- groups[[i]]
    local_idx <- list(person = rep(1L, length(rows)), score_k = idx$score_k[rows],
      weight = weights[rows], step_idx = idx$step_idx[rows], slope_idx = idx$slope_idx[rows])
    evaluate <- mfrmr_adaptive_person_kernel(
      local_idx, config, params, base_eta[rows], mu[person], sigma
    )
    location <- mfrmr_adaptive_posterior_location(evaluate)
    z <- location$mode + location$scale * quad$nodes
    basis$nodes[person, ] <- mu[person] + sigma * z
    basis$log_weights[person, ] <- log(quad$weights) + stats::dnorm(z, log = TRUE) -
      stats::dnorm(quad$nodes, log = TRUE) + log(location$scale)
  }
  basis
}

# Shared one-Person kernel, with coordinates standardized to the original prior.
mfrmr_adaptive_person_kernel <- function(idx, config, params, base_eta, mu, sigma) {
  config$n_person <- 1L
  config$population_spec <- list(active = FALSE)
  categories <- seq_len(config$n_cat) - 1L
  slope <- if (identical(config$model, "GPCM")) params$slopes[idx$slope_idx] else
    rep(1, length(idx$score_k))
  if (any(!is.finite(slope)) || any(slope <= 0)) {
    stop("The retained slopes are not positive finite values.", call. = FALSE)
  }
  # Thresholds are fixed for this parameter vector throughout mode finding,
  # quadrature and CDF integration; reuse the same cumulative sums.
  step_cum <- if (config$model == "RSM") {
    c(0, cumsum(params$steps))
  } else {
    t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
  }
  n <- length(idx$score_k)
  rows <- seq_len(n)
  step_rows <- if (config$model == "RSM") {
    matrix(step_cum, n, length(categories), byrow = TRUE)
  } else step_cum[idx$step_idx, , drop = FALSE]
  function(z, include_linear_part = FALSE, include_moments = TRUE) {
    # CDF integration needs only likelihoods; retain the R moment/gradient path.
    if (!include_moments && mfrm_use_cpp11_backend(config, include_linear_part)) {
      bundle <- mfrm_mml_logprob_bundle_cpp11(
        idx, config, list(nodes = mu + sigma * z, weights = rep(1, length(z))),
        params, base_eta, step_cum = step_cum, include_probs = FALSE
      )
      return(list(log_likelihood = colSums(bundle$log_prob_mat)))
    }
    # Stack nodes within this Person so fixed row/category work is not repeated.
    node_rows <- rep(rows, length(z))
    eta <- base_eta[node_rows] + rep(mu + sigma * z, each = n)
    linear <- outer(eta, categories) - step_rows[node_rows, , drop = FALSE]
    log_num <- if (identical(config$model, "GPCM")) linear * slope[node_rows] else linear
    flat_rows <- seq_along(node_rows)
    row_max <- log_num[cbind(flat_rows, max.col(log_num, ties.method = "first"))]
    log_num <- log_num - row_max
    log_denom <- log(rowSums(exp(log_num)))
    log_prob <- log_num[cbind(flat_rows, idx$score_k[node_rows] + 1L)] - log_denom
    if (!is.null(idx$weight)) log_prob <- log_prob * idx$weight[node_rows]
    log_likelihood <- colSums(matrix(log_prob, n, length(z)))
    if (!include_moments) return(list(log_likelihood = log_likelihood))
    probability <- exp(log_num - log_denom)
    node_indices <- lapply(seq_along(z), function(q) rows + (q - 1L) * n)
    prob_list <- lapply(node_indices, function(i) probability[i, , drop = FALSE])
    linear_part_list <- if (include_linear_part) {
      if (identical(config$model, "GPCM")) {
        lapply(node_indices, function(i) linear[i, , drop = FALSE])
      } else vector("list", length(z))
    } else NULL
    # Rows are already stacked by node. Compute all category moments together
    # instead of dispatching the same small matrix operations once per node.
    expected <- as.vector(probability %*% categories)
    variance <- rowSums(probability *
      (rep(categories, each = nrow(probability)) - expected)^2)
    moments <- rbind(
      score = sigma * colSums(matrix(
        (idx$weight * slope)[node_rows] * (idx$score_k[node_rows] - expected), n)),
      information = sigma^2 * colSums(matrix(
        (idx$weight * slope^2)[node_rows] * variance, n)))
    if (any(!is.finite(moments))) {
      stop(new_adaptive_quadrature_numeric_error("Non-finite adaptive posterior score or information."))
    }
    list(log_likelihood = log_likelihood,
         score = moments[1L, ] - z, information = moments[2L, ] + 1,
         prob_list = prob_list, linear_part_list = linear_part_list)
  }
}

mfrmr_validate_adaptive_quad_points <- function(points) {
  if (is.null(points)) return(NULL)
  if (!is.numeric(points) || length(points) < 2L || any(!is.finite(points)) ||
      any(points < 3 | points > .Machine$integer.max | points != floor(points)) ||
      length(unique(points)) < 2L) {
    stop("`adaptive_quad_points` must be NULL or at least two distinct integer orders >= 3.",
         call. = FALSE)
  }
  points <- sort(unique(as.integer(points)))
  # Validate representability once, before fitting or scoring.
  invisible(lapply(points, gauss_hermite_normal))
  points
}

mfrmr_adaptive_quadrature_review <- function(
    idx, config, params, quad, person_labels, adaptive_quad_points,
    population_spec = NULL, base_eta = NULL) {
  rules <- lapply(adaptive_quad_points, gauss_hermite_normal)
  basis <- resolve_person_quadrature_basis(
    list(nodes = c(0, 1), weights = c(1, 1)),
    population_spec = population_spec %||% materialize_population_spec(config, params),
    person_count = length(person_labels)
  )
  mu <- if (isTRUE(basis$transformed)) basis$mu else rep(0, length(person_labels))
  sigma <- if (isTRUE(basis$transformed)) basis$sigma else 1
  base_eta <- base_eta %||% compute_base_eta(idx, params, config)
  weights <- idx$weight %||% rep(1, length(idx$score_k))
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Adaptive quadrature review requires finite non-negative observation weights.",
         call. = FALSE)
  }
  person_ids <- sort(unique(idx$person))
  groups <- split(seq_along(idx$person), factor(idx$person, levels = person_ids))

  rows <- lapply(seq_along(groups), function(index) {
    person <- person_ids[index]
    observations <- groups[[index]]
    local_idx <- list(
      person = rep(1L, length(observations)),
      score_k = idx$score_k[observations], weight = weights[observations],
      step_idx = idx$step_idx[observations], slope_idx = idx$slope_idx[observations]
    )
    out <- data.frame(
      Person = rep(as.character(person_labels[person]), length(rules)),
      FixedNodes = length(quad$nodes), AdaptiveNodes = adaptive_quad_points,
      FixedLogMarginal = NA_real_, AdaptiveLogMarginal = NA_real_,
      LogMarginalChange = NA_real_,
      FixedEAP = NA_real_, AdaptiveEAP = NA_real_, EAPChange = NA_real_,
      FixedPosteriorSD = NA_real_, AdaptivePosteriorSD = NA_real_,
      PosteriorSDChange = NA_real_,
      AdaptiveLogMarginalChangeFromPrevious = NA_real_,
      AdaptiveEAPChangeFromPrevious = NA_real_,
      AdaptiveSDChangeFromPrevious = NA_real_,
      PosteriorMode = NA_real_, LocalPosteriorSD = NA_real_,
      ModeScore = NA_real_, Status = "unavailable", Detail = "",
      stringsAsFactors = FALSE
    )
    tryCatch({
      evaluate <- mfrmr_adaptive_person_kernel(
        local_idx, config, params, base_eta[observations], mu[person], sigma
      )
      summarize <- function(z, likelihood, log_weights) {
        offset <- max(likelihood)
        joint <- (likelihood - offset) + log_weights
        if (any(!is.finite(joint))) {
          stop(new_adaptive_quadrature_numeric_error("Non-finite integration terms."))
        }
        normalizer <- logsumexp(joint)
        probability <- exp(joint - max(joint))
        probability <- probability / sum(probability)
        center <- sum(probability * z)
        value <- c(log_marginal = offset + normalizer, eap = mu[person] + sigma * center,
                   sd = sigma * sqrt(sum(probability * (z - center)^2)))
        if (any(!is.finite(value))) stop("Non-finite posterior moments.", call. = FALSE)
        value
      }
      fixed <- summarize(quad$nodes,
        evaluate(quad$nodes, include_moments = FALSE)$log_likelihood, log(quad$weights))
      # For non-negative weights and positive slopes the normal-prior log
      # posterior is strictly concave: information = 1 + sum(w a^2 Var(K)).
      location <- mfrmr_adaptive_posterior_location(evaluate)
      mode <- location$mode
      scale <- location$scale
      adaptive <- vapply(rules, function(rule) {
        z <- mode + scale * rule$nodes
        # The density ratio and Jacobian retain the original prior/integral.
        log_weights <- log(rule$weights) +
          stats::dnorm(z, log = TRUE) - stats::dnorm(rule$nodes, log = TRUE) + log(scale)
        summarize(z, evaluate(z, include_moments = FALSE)$log_likelihood, log_weights)
      }, c(log_marginal = 0, eap = 0, sd = 0))
      out$FixedLogMarginal <- fixed["log_marginal"]
      out$AdaptiveLogMarginal <- adaptive["log_marginal", ]
      out$LogMarginalChange <- out$AdaptiveLogMarginal - out$FixedLogMarginal
      out$FixedEAP <- fixed["eap"]
      out$AdaptiveEAP <- adaptive["eap", ]
      out$EAPChange <- out$AdaptiveEAP - out$FixedEAP
      out$FixedPosteriorSD <- fixed["sd"]
      out$AdaptivePosteriorSD <- adaptive["sd", ]
      out$PosteriorSDChange <- out$AdaptivePosteriorSD - out$FixedPosteriorSD
      out$AdaptiveLogMarginalChangeFromPrevious <- c(NA_real_, diff(out$AdaptiveLogMarginal))
      out$AdaptiveEAPChangeFromPrevious <- c(NA_real_, diff(out$AdaptiveEAP))
      out$AdaptiveSDChangeFromPrevious <- c(NA_real_, diff(out$AdaptivePosteriorSD))
      out$PosteriorMode <- mu[person] + sigma * mode
      out$LocalPosteriorSD <- sigma * scale
      out$ModeScore <- location$score
      out$Status <- "computed"
      out
    }, error = function(condition) {
      out$Detail <- conditionMessage(condition)
      out
    })
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mfrmr_adaptive_quadrature_note <- function() {
  paste(
    "The quadrature_review table holds calibration parameters and the prior fixed",
    "and compares a fixed-prior grid with mode/curvature-adapted grids.",
    "Inspect both fixed/adaptive differences and changes between adaptive orders.",
    "Computed rows are numerical diagnostics, not certified accuracy or revised scores;",
    "unavailable rows retain the reason. Fit and score readiness are unchanged."
  )
}

mfrmr_adaptive_quadrature_overview <- function(review) {
  if (is.null(review) || nrow(review) == 0L) return(data.frame())
  groups <- split(review, interaction(review$FixedNodes, review$AdaptiveNodes, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(rows) {
    maximum <- function(column) {
      values <- abs(rows[[column]][rows$Status == "computed"])
      if (length(values) == 0L) NA_real_ else max(values, na.rm = TRUE)
    }
    data.frame(FixedNodes = rows$FixedNodes[1L], AdaptiveNodes = rows$AdaptiveNodes[1L],
      Persons = nrow(rows), Unavailable = sum(rows$Status != "computed"),
      MaxAbsLogMarginalChange = maximum("LogMarginalChange"),
      MaxAbsEAPChange = maximum("EAPChange"),
      MaxAbsSDChange = maximum("PosteriorSDChange"))
  }))
  rownames(out) <- NULL
  out
}

# The gradient differentiates the finite adaptive sum, including mode, scale
# and Jacobian, in the optimizer's actual free coordinates.
mfrmr_make_adaptive_mml_evaluator <- function(idx, config, sizes, quad_points,
                                               reuse_identical_patterns = TRUE,
                                               include_person_summaries = FALSE) {
  if (!identical(config$method, "MML")) {
    stop("Adaptive MML evaluation requires an MML configuration.", call. = FALSE)
  }
  rule <- gauss_hermite_normal(quad_points)
  n <- length(idx$score_k)
  p <- sum(unlist(sizes, use.names = FALSE))
  categories <- seq_len(config$n_cat) - 1L
  weights <- idx$weight %||% rep(1, n)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Adaptive MML requires finite non-negative observation weights.", call. = FALSE)
  }
  slices <- build_param_slices(sizes)
  adjacent <- mfrmr_estimability_adjacent_design(
    NULL, idx, config, sizes, include_person = FALSE, include_population_beta = FALSE
  )
  mapping <- Matrix::sparseMatrix(
    i = seq_len(ncol(adjacent$design)), j = adjacent$map$OptimizerIndex, x = 1,
    dims = c(ncol(adjacent$design), p)
  )
  adjacent <- adjacent$design %*% mapping
  cumulative <- lapply(seq_len(config$n_cat - 1L), function(k) {
    adjacent[seq_len(n) + (k - 1L) * n, , drop = FALSE]
  })
  cumulative <- do.call(rbind, Reduce(`+`, cumulative, accumulate = TRUE))
  log_slope_design <- mfrmr_empty_sparse_matrix(n, p)
  if (length(slices[["log_slopes"]])) {
    log_slope_design[, slices[["log_slopes"]]] <-
      gpcm_log_slope_design(config$gpcm_spec)[idx$slope_idx, , drop = FALSE]
  }
  patterns <- mfrmr_mml_person_design_groups(idx, config,
    reuse_identical_designs = reuse_identical_patterns, include_responses = TRUE)
  person_ids <- patterns$person_ids[patterns$representative_positions]
  groups <- patterns$canonical_rows[patterns$representative_positions]
  multiplicity <- patterns$group_sizes
  # Keep category designs sparse. Only the current Person's probabilities are dense.
  designs <- lapply(groups, function(rows) {
    list(cumulative = cumulative[as.vector(outer(rows, (categories[-1L] - 1L) * n, `+`)),
                                 , drop = FALSE],
         log_slope = log_slope_design[rows, , drop = FALSE])
  })
  function(par) {
    if (length(par) != p || any(!is.finite(par))) {
      stop("Adaptive MML requires a finite free-parameter vector of the expected length.", call. = FALSE)
    }
    params <- expand_params(par, sizes, config)
    base <- compute_base_eta(idx, params, config)
    population <- materialize_population_spec(config, params)
    basis <- resolve_person_quadrature_basis(
      list(nodes = 0, weights = 1), population, person_count = config$n_person
    )
    mu <- if (isTRUE(basis$transformed)) basis$mu else rep(0, config$n_person)
    sigma <- if (isTRUE(basis$transformed)) basis$sigma else 1
    value <- 0
    gradient <- numeric(p)
    if (include_person_summaries) {
      person_summaries <- matrix(NA_real_, length(patterns$person_ids), 3L,
        dimnames = list(NULL, c("LogMarginal", "EAP", "PosteriorSD")))
    }
    for (i in seq_along(groups)) {
      person <- person_ids[i]
      rows <- groups[[i]]
      local_idx <- list(person = rep(1L, length(rows)), score_k = idx$score_k[rows],
        weight = weights[rows], step_idx = idx$step_idx[rows], slope_idx = idx$slope_idx[rows])
      slope <- if (identical(config$model, "GPCM")) params$slopes[local_idx$slope_idx] else
        rep(1, length(rows))
      evaluate <- mfrmr_adaptive_person_kernel(
        local_idx, config, params, base[rows], mu[person], sigma
      )
      location <- mfrmr_adaptive_posterior_location(evaluate)
      mode <- location$mode
      at_mode <- evaluate(mode, include_linear_part = TRUE)
      scale <- location$scale
      mean_derivative <- variance_derivative <- numeric(p)
      if (length(slices[["beta"]])) {
        mean_derivative[slices[["beta"]]] <-
          population$design_matrix[population$person_lookup[person], ]
      }
      if (length(slices[["log_sigma2"]])) variance_derivative[slices[["log_sigma2"]]] <- 1
      # Projects category-logit derivatives into the optimizer's actual free
      # coordinates, including anchored steps and signed facet interactions.
      project <- function(coefficient, linear) {
        answer <- as.numeric(Matrix::crossprod(designs[[i]]$cumulative,
          as.vector(coefficient[, -1L, drop = FALSE] * slope)))
        if (length(slices[["log_slopes"]])) {
          answer <- answer + as.numeric(Matrix::crossprod(designs[[i]]$log_slope,
            rowSums(coefficient * linear) * slope))
        }
        answer
      }
      probability <- at_mode$prob_list[[1L]]
      expected <- as.vector(probability %*% categories)
      centered <- rep(categories, each = length(rows)) - expected
      variance <- rowSums(probability * centered^2)
      w <- local_idx$weight
      linear <- at_mode$linear_part_list[[1L]]
      # In physical theta coordinates: dm = h_theta,psi/H and
      # d log(s) = -(H_psi + H_theta dm)/(2H), H = -h_theta,theta.
      mixed <- -project(probability * centered * (w * slope), linear) +
        as.numeric(Matrix::crossprod(designs[[i]]$log_slope,
          w * slope * (local_idx$score_k - expected))) +
        mean_derivative / sigma^2 + mode / sigma * variance_derivative
      curvature_derivative <- project(
        probability * (centered^2 - variance) * (w * slope^2), linear
      ) + as.numeric(Matrix::crossprod(designs[[i]]$log_slope, 2 * w * slope^2 * variance)) -
        variance_derivative / sigma^2
      curvature_theta <- sum(w * slope^3 * rowSums(probability * centered^3))
      curvature <- at_mode$information / sigma^2
      mode_derivative <- mixed / curvature
      log_scale_derivative <- -(curvature_derivative + curvature_theta * mode_derivative) /
        (2 * curvature)
      z <- mode + scale * rule$nodes
      nodes <- evaluate(z, include_linear_part = TRUE)
      offset <- max(nodes$log_likelihood)
      joint <- log(rule$weights) + (nodes$log_likelihood - offset) + stats::dnorm(z, log = TRUE) -
        stats::dnorm(rule$nodes, log = TRUE) + log(scale)
      if (any(!is.finite(joint))) {
        stop(new_adaptive_quadrature_numeric_error("Non-finite integration terms."))
      }
      normalizer <- logsumexp(joint)
      posterior <- exp(joint - max(joint))
      posterior <- posterior / sum(posterior)
      # Independent Persons multiply marginal likelihoods; their count belongs
      # outside the integral, never in the observation weights.
      value <- value - multiplicity[i] * (offset + normalizer)
      if (include_person_summaries) {
        center <- sum(posterior * z)
        moments <- c(offset + normalizer, mu[person] + sigma * center,
          sigma * sqrt(sum(posterior * (z - center)^2)))
        # Expand exact pattern reuse back to every independent Person.
        person_summaries[patterns$group_positions[[i]], ] <-
          matrix(moments, multiplicity[i], 3L, byrow = TRUE)
      }
      # Projection is linear: average category derivatives before multiplying
      # by the same sparse designs. Keep the slope's residual-times-logit
      # average separate; a product of averages would change the gradient.
      residual_mean <- matrix(0, length(rows), config$n_cat)
      slope_mean <- numeric(length(rows))
      observed <- cbind(seq_along(rows), local_idx$score_k + 1L)
      for (q in seq_along(rule$nodes)) {
        residual <- -nodes$prob_list[[q]]
        residual[observed] <- residual[observed] + 1
        residual_mean <- residual_mean + posterior[q] * residual
        if (length(slices[["log_slopes"]])) {
          slope_mean <- slope_mean + posterior[q] *
            rowSums(residual * nodes$linear_part_list[[q]])
        }
      }
      fixed_derivative <- as.numeric(Matrix::crossprod(designs[[i]]$cumulative,
        as.vector(residual_mean[, -1L, drop = FALSE] * (w * slope))))
      if (length(slices[["log_slopes"]])) {
        fixed_derivative <- fixed_derivative + as.numeric(Matrix::crossprod(
          designs[[i]]$log_slope, slope_mean * (w * slope)))
      }
      fixed_derivative <- fixed_derivative +
        sum(posterior * z) / sigma * mean_derivative +
        sum(posterior * (z^2 - 1)) / 2 * variance_derivative
      gradient <- gradient - multiplicity[i] * (fixed_derivative +
        sum(posterior * nodes$score) / sigma * mode_derivative +
        sum(posterior * (scale * nodes$score * rule$nodes + 1)) * log_scale_derivative)
    }
    if (!is.finite(value) || any(!is.finite(gradient))) {
      stop(new_adaptive_quadrature_numeric_error("Non-finite adaptive MML objective or gradient."))
    }
    out <- list(value = value, gradient = gradient)
    if (include_person_summaries) {
      if (any(!is.finite(person_summaries))) {
        stop(new_adaptive_quadrature_numeric_error("Non-finite posterior moments."))
      }
      out$per_person <- data.frame(PersonIndex = patterns$person_ids, person_summaries)
    }
    out
  }
}
