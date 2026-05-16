# Core MFRM computation functions extracted from app.R (non-UI section)
# Source: ../../app.R
# NOTE: This file intentionally keeps function names aligned with the app implementation.
# ---- math helpers ----

# Numerically stable log-sum-exp: log(sum(exp(x))).
# Subtracts max(x) before exponentiation to avoid overflow/underflow.
#   logsumexp(x) = max(x) + log(sum(exp(x - max(x))))
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Weighted mean that silently drops non-finite values and zero weights.
weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# row_max_fast() and the three category_prob_* helpers now live in
# R/core-category-probabilities.R. They are loaded with the rest of
# the package namespace before any function is called, so downstream
# code in this file continues to work without change.

mfrm_cpp11_backend_available <- function() {
  exists("mfrm_cpp_rsm_logprob_bundle", mode = "function") &&
    exists("mfrm_cpp_pcm_logprob_bundle", mode = "function") &&
    exists("mfrm_cpp_expected_category_bundle", mode = "function") &&
    exists("mfrm_cpp_compute_p_geq", mode = "function")
}

mfrm_use_cpp11_backend <- function(config, include_linear_part = FALSE) {
  isTRUE(getOption("mfrmr.use_cpp11_backend", FALSE)) &&
    mfrm_cpp11_backend_available() &&
    !isTRUE(include_linear_part) &&
    identical(config$model %in% c("RSM", "PCM"), TRUE)
}

get_weights <- function(df) {
  if (!is.null(df) && "Weight" %in% names(df)) {
    w <- suppressWarnings(as.numeric(df$Weight))
    w <- ifelse(is.finite(w) & w > 0, w, 0)
    return(w)
  }
  rep(1, nrow(df))
}

stop_if_gpcm_out_of_scope <- function(fit,
                                      helper,
                                      supported = "fitting, core summary output, fixed-calibration posterior scoring, compute_information(), direct simulation and recovery checks, residual-based diagnostics, the curve/report helpers, and the slope-aware element-conditional fair_average_table() and estimate_bias() (with the SE caveats documented in their help pages)") {
  if (!inherits(fit, "mfrm_fit")) return(invisible(NULL))
  model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  if (identical(model, "GPCM")) {
    stop(
      "`", helper, "` does not support `GPCM` fits in this release. ",
      "Supported `GPCM` workflows are: ",
      supported,
      ". See `gpcm_capability_matrix()` for the current scope.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

gpcm_fair_average_rationale <- function() {
  paste0(
    "The diagnostics/QC dashboard keeps the embedded fair-average panel ",
    "disabled for bounded `GPCM` fits so it does not present score-side ",
    "adjustments as a fully validated FACETS-equivalent reporting surface. ",
    "Use `fair_average_table()` directly for the supported slope-aware ",
    "element-conditional GPCM fair averages, and keep its documented SE ",
    "caveat in force."
  )
}

gpcm_planning_scope_rationale <- function() {
  paste0(
    "The current planning layer targets the role-based person x ",
    "rater-like x criterion-like design contract rather than a fully ",
    "arbitrary-facet planner."
  )
}

# Gauss-Hermite nodes/weights for standard normal integration
# Based on Golub-Welsch for Hermite polynomials
# Returns nodes and weights for phi(x) dx
#   integral f(x) phi(x) dx \approx sum w_i f(x_i)
gauss_hermite_normal <- function(n) {
  if (n < 1) stop("Gauss-Hermite quadrature requires n >= 1 quadrature points. ",
                   "Check the 'quad_points' argument.", call. = FALSE)
  if (n == 1) {
    return(list(nodes = 0, weights = 1))
  }
  i <- seq_len(n - 1)
  a <- rep(0, n)
  b <- sqrt(i / 2)
  jmat <- matrix(0, nrow = n, ncol = n)
  diag(jmat) <- a
  jmat[cbind(seq_len(n - 1), seq_len(n - 1) + 1)] <- b
  jmat[cbind(seq_len(n - 1) + 1, seq_len(n - 1))] <- b
  eig <- eigen(jmat, symmetric = TRUE)
  nodes <- eig$values
  weights <- sqrt(pi) * (eig$vectors[1, ]^2)
  # Convert from exp(-x^2) to standard normal
  nodes <- sqrt(2) * nodes
  weights <- weights / sqrt(pi)
  list(nodes = nodes, weights = weights)
}

compact_population_spec <- function(population = NULL, person_levels = character(0)) {
  pop <- population %||% list()
  person_levels <- as.character(person_levels %||% character(0))
  person_table <- pop$person_table
  person_id <- as.character(pop$person_id %||% "Person")
  person_lookup <- integer(0)

  if (!is.null(person_table) &&
      is.data.frame(person_table) &&
      length(person_levels) > 0 &&
      person_id %in% names(person_table)) {
    person_lookup <- match(person_levels, as.character(person_table[[person_id]]))
  }

  list(
    active = isTRUE(pop$active),
    posterior_basis = as.character(
      pop$posterior_basis %||%
        if (isTRUE(pop$active)) "population_model" else "legacy_mml"
    ),
    formula = pop$formula %||% NULL,
    person_id = person_id,
    design_matrix = pop$design_matrix %||% NULL,
    design_columns = pop$design_columns %||%
      if (!is.null(pop$design_matrix)) colnames(pop$design_matrix) else NULL,
    xlevels = pop$xlevels %||% NULL,
    contrasts = pop$contrasts %||% NULL,
    coefficients = pop$coefficients %||% NULL,
    sigma2 = pop$sigma2 %||% NULL,
    converged = isTRUE(pop$converged),
    policy = pop$policy %||% NULL,
    person_lookup = person_lookup,
    person_rows = if (!is.null(person_table) && is.data.frame(person_table)) {
      nrow(person_table)
    } else {
      0L
    },
    included_persons = as.character(pop$included_persons %||% character(0)),
    omitted_persons = as.character(pop$omitted_persons %||% character(0)),
    response_rows_retained = as.integer(pop$response_rows_retained %||% NA_integer_),
    response_rows_omitted = as.integer(pop$response_rows_omitted %||% NA_integer_),
    notes = pop$notes %||% character(0)
  )
}

materialize_population_spec <- function(config, params = NULL) {
  pop <- config$population_spec %||% list()
  if (!is.null(params) && !is.null(params$population)) {
    if (!is.null(params$population$coefficients)) {
      pop$coefficients <- params$population$coefficients
    }
    if (!is.null(params$population$sigma2)) {
      pop$sigma2 <- params$population$sigma2
    }
  }
  pop
}

resolve_person_quadrature_basis <- function(quad,
                                            population_spec = NULL,
                                            person_ids = NULL,
                                            person_count = NULL) {
  if (is.null(person_ids)) {
    n_persons <- as.integer(person_count %||% 0L)
    person_ids <- seq_len(max(0L, n_persons))
  } else {
    person_ids <- as.integer(person_ids)
    n_persons <- length(person_ids)
  }
  n_nodes <- length(quad$nodes)
  log_weights <- matrix(log(quad$weights), nrow = n_persons, ncol = n_nodes, byrow = TRUE)
  default_nodes <- matrix(quad$nodes, nrow = n_persons, ncol = n_nodes, byrow = TRUE)
  pop <- population_spec %||% list()

  if (!isTRUE(pop$active)) {
    return(list(
      nodes = default_nodes,
      log_weights = log_weights,
      posterior_basis = as.character(pop$posterior_basis %||% "legacy_mml"),
      transformed = FALSE
    ))
  }

  if (is.null(pop$design_matrix) || is.null(pop$coefficients) || is.null(pop$sigma2)) {
    stop("Population-model quadrature requested before regression coefficients and variance were estimated.",
         call. = FALSE)
  }
  if (length(pop$person_lookup) == 0L || anyNA(pop$person_lookup[person_ids])) {
    stop("Population-model quadrature could not align all fitted persons to the latent-regression person table.",
         call. = FALSE)
  }

  design_matrix <- pop$design_matrix[pop$person_lookup[person_ids], , drop = FALSE]
  beta <- as.numeric(pop$coefficients)
  if (length(beta) != ncol(design_matrix)) {
    stop("Population-model coefficient length does not match the latent-regression design matrix.",
         call. = FALSE)
  }
  sigma2 <- as.numeric(pop$sigma2[1])
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("Population-model residual variance `sigma2` must be a single positive finite value.",
         call. = FALSE)
  }

  mu <- as.numeric(design_matrix %*% beta)
  sigma <- sqrt(sigma2)
  nodes <- outer(mu, quad$nodes, FUN = function(m, z) m + sigma * z)

  list(
    nodes = nodes,
    log_weights = log_weights,
    posterior_basis = as.character(pop$posterior_basis %||% "population_model"),
    transformed = TRUE,
    mu = mu,
    sigma = sigma
  )
}

center_sum_zero <- function(x) {
  if (length(x) == 0) return(x)
  x - mean(x)
}

expand_facet <- function(free, n_levels) {
  if (n_levels <= 1) return(rep(0, n_levels))
  c(free, -sum(free))
}

sum_zero_param_count <- function(n_levels) {
  max(as.integer(n_levels %||% 0L) - 1L, 0L)
}

expand_sum_zero_vector <- function(free, n_levels) {
  expand_facet(as.numeric(free %||% numeric(0)), as.integer(n_levels %||% 0L))
}

initial_sum_zero_free <- function(x) {
  x <- center_sum_zero(as.numeric(x %||% numeric(0)))
  if (length(x) <= 1L) numeric(0) else x[seq_len(length(x) - 1L)]
}

expand_step_matrix <- function(free, n_levels, n_steps) {
  n_levels <- as.integer(n_levels %||% 0L)
  n_steps <- as.integer(n_steps %||% 0L)
  if (n_levels <= 0L || n_steps <= 0L) {
    return(matrix(0, nrow = max(n_levels, 0L), ncol = max(n_steps, 0L)))
  }
  row_free_n <- sum_zero_param_count(n_steps)
  out <- matrix(0, nrow = n_levels, ncol = n_steps)
  if (row_free_n == 0L) {
    return(out)
  }
  free <- as.numeric(free %||% numeric(0))
  expected <- n_levels * row_free_n
  if (length(free) != expected) {
    stop("Step parameter vector has length ", length(free),
         " but the identified step parameterization requires ", expected, ".",
         call. = FALSE)
  }
  pos <- 1L
  for (i in seq_len(n_levels)) {
    out[i, ] <- expand_sum_zero_vector(free[pos:(pos + row_free_n - 1L)], n_steps)
    pos <- pos + row_free_n
  }
  out
}

project_step_matrix_gradient <- function(grad_step_mat) {
  if (!is.matrix(grad_step_mat)) {
    grad_step_mat <- as.matrix(grad_step_mat)
  }
  if (nrow(grad_step_mat) == 0L || ncol(grad_step_mat) <= 1L) {
    return(numeric(0))
  }
  unlist(
    lapply(seq_len(nrow(grad_step_mat)), function(i) {
      project_sum_zero_gradient(grad_step_mat[i, ])
    }),
    use.names = FALSE
  )
}

build_gpcm_slope_spec <- function(levels,
                                  slope_facet,
                                  step_facet) {
  slope_levels <- as.character(levels %||% character(0))
  list(
    active = TRUE,
    slope_facet = as.character(slope_facet %||% NA_character_),
    step_facet = as.character(step_facet %||% NA_character_),
    levels = slope_levels,
    n_params = max(length(slope_levels) - 1L, 0L),
    identification = "sum_to_zero_log_slopes",
    scale_reference = "geometric_mean_one",
    reduction_reference = "PCM when all slopes equal 1"
  )
}

expand_gpcm_log_slopes <- function(free,
                                   spec) {
  n_levels <- length(spec$levels %||% character(0))
  if (n_levels == 0L) {
    return(list(
      log_slopes = numeric(0),
      slopes = numeric(0)
    ))
  }

  log_slopes <- if ((spec$n_params %||% 0L) == 0L) {
    rep(0, n_levels)
  } else {
    expand_facet(as.numeric(free), n_levels)
  }

  list(
    log_slopes = log_slopes,
    slopes = exp(log_slopes)
  )
}

project_sum_zero_gradient <- function(grad_expanded) {
  grad_expanded <- as.numeric(grad_expanded %||% numeric(0))
  n <- length(grad_expanded)
  if (n <= 1L) {
    return(numeric(0))
  }
  grad_expanded[seq_len(n - 1L)] - grad_expanded[n]
}

facet_interactions_active <- function(x) {
  specs <- if (is.list(x) && !is.null(x$interaction_specs)) {
    x$interaction_specs
  } else {
    x
  }
  specs <- specs %||% list()
  if (!is.list(specs) || length(specs) == 0L) return(FALSE)
  any(vapply(specs, function(spec) {
    if (!is.list(spec) || is.null(spec$n_params)) return(FALSE)
    as.integer(spec$n_params %||% 0L) > 0L
  }, logical(1)))
}

normalize_facet_interactions <- function(facet_interactions, facet_names) {
  if (is.null(facet_interactions) || length(facet_interactions) == 0L) {
    return(list())
  }
  facet_names <- as.character(facet_names %||% character(0))

  terms <- if (is.character(facet_interactions)) {
    lapply(facet_interactions, function(term) {
      term <- trimws(as.character(term))
      if (is.na(term) || !nzchar(term)) {
        stop("`facet_interactions` cannot contain empty strings.", call. = FALSE)
      }
      parts <- trimws(strsplit(term, ":", fixed = TRUE)[[1]])
      if (length(parts) != 2L || any(!nzchar(parts))) {
        stop(
          "`facet_interactions` must use explicit two-way terms such as ",
          "`\"Rater:Criterion\"`. Problem term: ", shQuote(term), ".",
          call. = FALSE
        )
      }
      parts
    })
  } else if (is.list(facet_interactions) && !is.data.frame(facet_interactions)) {
    lapply(facet_interactions, function(term) {
      parts <- trimws(as.character(term))
      parts <- parts[!is.na(parts) & nzchar(parts)]
      if (length(parts) != 2L) {
        stop(
          "List entries in `facet_interactions` must each contain exactly ",
          "two facet names, for example `list(c(\"Rater\", \"Criterion\"))`.",
          call. = FALSE
        )
      }
      parts
    })
  } else {
    stop(
      "`facet_interactions` must be NULL, a character vector of `A:B` terms, ",
      "or a list of length-two character vectors.",
      call. = FALSE
    )
  }

  canonical <- character(length(terms))
  names_out <- character(length(terms))
  for (i in seq_along(terms)) {
    parts <- terms[[i]]
    if (any(parts == "Person")) {
      stop(
        "`facet_interactions` currently supports only non-person facet pairs. ",
        "Remove `Person` from ", shQuote(paste(parts, collapse = ":")), ".",
        call. = FALSE
      )
    }
    if (!all(parts %in% facet_names)) {
      missing <- setdiff(parts, facet_names)
      stop(
        "`facet_interactions` references facet(s) not supplied in `facets`: ",
        paste(shQuote(missing), collapse = ", "), ". Available non-person ",
        "facets are: ", paste(shQuote(facet_names), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (length(unique(parts)) != 2L) {
      stop(
        "`facet_interactions` terms must name two distinct facets. Problem term: ",
        shQuote(paste(parts, collapse = ":")), ".",
        call. = FALSE
      )
    }
    canonical[i] <- paste(sort(parts), collapse = "\r")
    names_out[i] <- paste(parts, collapse = ":")
  }
  if (any(duplicated(canonical))) {
    dup <- names_out[duplicated(canonical)]
    stop(
      "`facet_interactions` contains duplicate facet pairs: ",
      paste(shQuote(dup), collapse = ", "), ". Each two-way pair may be ",
      "specified only once.",
      call. = FALSE
    )
  }

  names(terms) <- names_out
  terms
}

build_facet_interaction_specs <- function(prep,
                                          facet_interactions = NULL,
                                          min_obs_per_interaction = 10,
                                          interaction_policy = c("warn", "error", "silent")) {
  interaction_policy <- match.arg(interaction_policy)
  if (!is.numeric(min_obs_per_interaction) ||
      length(min_obs_per_interaction) != 1L ||
      !is.finite(min_obs_per_interaction) ||
      min_obs_per_interaction < 0) {
    stop(
      "`min_obs_per_interaction` must be a single non-negative finite number.",
      call. = FALSE
    )
  }
  min_obs_per_interaction <- as.numeric(min_obs_per_interaction)

  terms <- normalize_facet_interactions(facet_interactions, prep$facet_names)
  if (length(terms) == 0L) return(list())

  specs <- lapply(names(terms), function(name) {
    facets <- terms[[name]]
    facet_a <- facets[1]
    facet_b <- facets[2]
    levels_a <- prep$levels[[facet_a]]
    levels_b <- prep$levels[[facet_b]]
    n_a <- length(levels_a)
    n_b <- length(levels_b)
    if (n_a < 2L || n_b < 2L) {
      stop(
        "Facet interaction ", shQuote(name), " is not estimable because ",
        "both facets must have at least two observed levels.",
        call. = FALSE
      )
    }

    a_idx <- as.integer(prep$data[[facet_a]])
    b_idx <- as.integer(prep$data[[facet_b]])
    cell_id <- a_idx + (b_idx - 1L) * n_a
    n_cells <- n_a * n_b
    cell_n <- matrix(tabulate(cell_id, nbins = n_cells), nrow = n_a, ncol = n_b)

    w <- suppressWarnings(as.numeric(prep$data$Weight))
    cell_weighted <- numeric(n_cells)
    weight_sum <- rowsum(matrix(w, ncol = 1), cell_id, reorder = FALSE)
    weight_ids <- as.integer(rownames(weight_sum))
    cell_weighted[weight_ids] <- as.vector(weight_sum)
    cell_weighted <- matrix(cell_weighted, nrow = n_a, ncol = n_b)

    sparse_cells <- cell_n < min_obs_per_interaction
    sparse_count <- sum(sparse_cells)
    zero_count <- sum(cell_n == 0L)
    if (sparse_count > 0L && !identical(interaction_policy, "silent")) {
      msg <- paste0(
        "Facet interaction ", shQuote(name), " has ", sparse_count,
        " of ", n_cells, " cells below `min_obs_per_interaction = ",
        min_obs_per_interaction, "` (", zero_count, " zero-count cells). ",
        "Interaction estimates in sparse cells are weakly identified; ",
        "use the term only when the pair is substantively justified."
      )
      if (identical(interaction_policy, "error")) {
        stop(msg, call. = FALSE)
      }
      warning(msg, call. = FALSE)
    }

    list(
      name = name,
      facets = facets,
      facet_a = facet_a,
      facet_b = facet_b,
      levels_a = levels_a,
      levels_b = levels_b,
      n_a = n_a,
      n_b = n_b,
      n_cells = n_cells,
      n_params = (n_a - 1L) * (n_b - 1L),
      cell_n = cell_n,
      cell_weighted_n = cell_weighted,
      sparse_count = sparse_count,
      zero_count = zero_count,
      min_obs_per_interaction = min_obs_per_interaction,
      identification = "two_way_sum_to_zero_margins",
      interpretation = paste(
        "Positive estimates indicate scores higher than expected under the",
        "additive main-effects model for that facet-level combination."
      )
    )
  })
  names(specs) <- names(terms)
  specs
}

expand_two_way_interaction <- function(free, spec) {
  n_a <- as.integer(spec$n_a %||% 0L)
  n_b <- as.integer(spec$n_b %||% 0L)
  mat <- matrix(0, nrow = n_a, ncol = n_b)
  if (n_a <= 1L || n_b <= 1L) return(mat)
  free <- as.numeric(free %||% numeric(0))
  expected <- (n_a - 1L) * (n_b - 1L)
  if (length(free) != expected) {
    stop(
      "Internal interaction parameter length mismatch for ",
      shQuote(spec$name %||% "interaction"), ": expected ", expected,
      " free parameter(s), got ", length(free), ".",
      call. = FALSE
    )
  }
  core <- matrix(free, nrow = n_a - 1L, ncol = n_b - 1L)
  mat[seq_len(n_a - 1L), seq_len(n_b - 1L)] <- core
  mat[seq_len(n_a - 1L), n_b] <- -rowSums(core)
  mat[n_a, seq_len(n_b - 1L)] <- -colSums(core)
  mat[n_a, n_b] <- sum(core)
  dimnames(mat) <- list(spec$levels_a, spec$levels_b)
  mat
}

project_two_way_interaction_gradient <- function(grad_expanded, spec) {
  n_a <- as.integer(spec$n_a %||% nrow(grad_expanded))
  n_b <- as.integer(spec$n_b %||% ncol(grad_expanded))
  if (n_a <= 1L || n_b <= 1L) return(numeric(0))
  grad_expanded <- matrix(as.numeric(grad_expanded), nrow = n_a, ncol = n_b)
  core <- grad_expanded[seq_len(n_a - 1L), seq_len(n_b - 1L), drop = FALSE]
  last_col <- grad_expanded[seq_len(n_a - 1L), n_b]
  last_row <- grad_expanded[n_a, seq_len(n_b - 1L)]
  corner <- grad_expanded[n_a, n_b]
  projected <- core -
    matrix(last_col, nrow = n_a - 1L, ncol = n_b - 1L) -
    matrix(last_row, nrow = n_a - 1L, ncol = n_b - 1L, byrow = TRUE) +
    corner
  as.vector(projected)
}

expand_interaction_params <- function(free, specs) {
  specs <- specs %||% list()
  if (length(specs) == 0L) return(list())
  out <- vector("list", length(specs))
  names(out) <- names(specs)
  pos <- 1L
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    k <- as.integer(spec$n_params %||% 0L)
    seg <- if (k > 0L) {
      free[pos:(pos + k - 1L)]
    } else {
      numeric(0)
    }
    out[[nm]] <- expand_two_way_interaction(seg, spec)
    pos <- pos + k
  }
  out
}

compute_interaction_eta <- function(idx, params, config) {
  specs <- config$interaction_specs %||% list()
  if (length(specs) == 0L) return(rep(0, length(idx$score_k)))
  interactions <- params$interactions %||% list()
  eta <- rep(0, length(idx$score_k))
  for (nm in names(specs)) {
    mat <- interactions[[nm]]
    cell_id <- idx$interactions[[nm]]
    if (is.null(mat) || is.null(cell_id)) next
    eta <- eta + as.vector(mat)[cell_id]
  }
  eta
}

compute_interaction_gradient_free <- function(residual, idx, config) {
  specs <- config$interaction_specs %||% list()
  if (length(specs) == 0L) return(numeric(0))
  residual <- as.numeric(residual)
  out <- lapply(names(specs), function(nm) {
    spec <- specs[[nm]]
    k <- as.integer(spec$n_params %||% 0L)
    if (k == 0L) return(numeric(0))
    cell_id <- idx$interactions[[nm]]
    if (is.null(cell_id)) {
      stop("Internal interaction index not found for ", shQuote(nm), ".", call. = FALSE)
    }
    grad_vec <- numeric(spec$n_cells)
    rs <- rowsum(matrix(residual, ncol = 1), cell_id, reorder = FALSE)
    ids <- as.integer(rownames(rs))
    grad_vec[ids] <- as.vector(rs)
    grad_mat <- matrix(grad_vec, nrow = spec$n_a, ncol = spec$n_b)
    project_two_way_interaction_gradient(grad_mat, spec)
  })
  unlist(out, use.names = FALSE)
}

build_interaction_effect_table <- function(config, prep, params) {
  specs <- config$interaction_specs %||% list()
  if (length(specs) == 0L) return(tibble())
  interactions <- params$interactions %||% list()
  tbls <- lapply(names(specs), function(nm) {
    spec <- specs[[nm]]
    mat <- interactions[[nm]]
    if (is.null(mat)) {
      mat <- matrix(0, nrow = spec$n_a, ncol = spec$n_b)
    }
    grid <- expand.grid(
      FacetA_Level = spec$levels_a,
      FacetB_Level = spec$levels_b,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    tibble(
      Interaction = spec$name,
      FacetA = spec$facet_a,
      FacetA_Level = as.character(grid$FacetA_Level),
      FacetB = spec$facet_b,
      FacetB_Level = as.character(grid$FacetB_Level),
      Estimate = as.vector(mat),
      N = as.integer(as.vector(spec$cell_n)),
      WeightedN = as.numeric(as.vector(spec$cell_weighted_n)),
      Sparse = as.vector(spec$cell_n) < spec$min_obs_per_interaction,
      Identification = spec$identification
    )
  })
  bind_rows(tbls)
}

build_facet_constraint <- function(levels,
                                   anchors = NULL,
                                   groups = NULL,
                                   group_values = NULL,
                                   centered = TRUE) {
  lvl <- as.character(levels)
  anchors_vec <- rep(NA_real_, length(lvl))
  names(anchors_vec) <- lvl
  if (!is.null(anchors)) {
    anchors <- anchors[!is.na(names(anchors))]
    anchors_vec[names(anchors)] <- as.numeric(anchors)
  }

  groups_vec <- rep(NA_character_, length(lvl))
  names(groups_vec) <- lvl
  if (!is.null(groups)) {
    groups <- groups[!is.na(names(groups))]
    groups_vec[names(groups)] <- as.character(groups)
  }

  group_values_map <- numeric(0)
  if (!is.null(group_values)) {
    group_values <- group_values[!is.na(names(group_values))]
    group_values_map <- as.numeric(group_values)
    names(group_values_map) <- names(group_values)
  }

  spec <- list(
    levels = lvl,
    anchors = anchors_vec,
    groups = groups_vec,
    group_values = group_values_map,
    centered = centered
  )
  spec$n_params <- count_facet_params(spec)
  spec
}

count_facet_params <- function(spec) {
  anchors <- spec$anchors
  groups <- spec$groups
  free_idx <- which(is.na(anchors))
  if (length(free_idx) == 0) return(0)

  n_params <- 0
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(anchors[group_levels])]
      k <- length(free_in_group)
      if (k > 1) n_params <- n_params + (k - 1)
    }
  }

  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m == 0) return(n_params)
  if (spec$centered) {
    n_params <- n_params + max(m - 1, 0)
  } else {
    n_params <- n_params + m
  }
  n_params
}

expand_facet_with_constraints <- function(free, spec) {
  out <- spec$anchors
  groups <- spec$groups
  group_values <- spec$group_values
  centered <- isTRUE(spec$centered)
  free_idx <- which(is.na(out))
  if (length(free_idx) == 0) return(out)

  used <- 0
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(out[group_levels])]
      if (length(free_in_group) == 0) next
      group_value <- if (gid %in% names(group_values)) group_values[[gid]] else 0
      anchor_sum <- sum(out[group_levels], na.rm = TRUE)
      target_sum <- group_value * length(group_levels)
      k <- length(free_in_group)
      if (k == 1) {
        out[free_in_group] <- target_sum - anchor_sum
      } else {
        seg <- free[(used + 1):(used + k - 1)]
        used <- used + (k - 1)
        last_val <- target_sum - anchor_sum - sum(seg)
        out[free_in_group] <- c(seg, last_val)
      }
    }
  }

  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m == 0) return(out)
  if (centered) {
    if (m == 1) {
      out[ungrouped_idx] <- 0
    } else {
      seg <- free[(used + 1):(used + m - 1)]
      used <- used + (m - 1)
      out[ungrouped_idx] <- c(seg, -sum(seg))
    }
  } else {
    seg <- free[(used + 1):(used + m)]
    used <- used + m
    out[ungrouped_idx] <- seg
  }
  out
}

# ---- Constraint gradient projection ----
# Reverse of expand_facet_with_constraints: projects gradient from expanded
# (full-length) parameter space back to the free (optimizable) parameter space.
# This computes J^T * grad_expanded, where J is the Jacobian of
# expand_facet_with_constraints w.r.t. the free parameters.
constraint_grad_project <- function(grad_expanded, spec) {
  anchors <- spec$anchors
  groups <- spec$groups
  centered <- isTRUE(spec$centered)
  free_idx <- which(is.na(anchors))
  if (length(free_idx) == 0 || spec$n_params == 0) return(numeric(0))

  grad_free <- numeric(spec$n_params)
  pos <- 1L

  # Grouped levels: last free in group = target - anchor_sum - sum(others)
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(anchors[group_levels])]
      k <- length(free_in_group)
      if (k <= 1) next
      last_g <- grad_expanded[free_in_group[k]]
      for (j in seq_len(k - 1)) {
        grad_free[pos] <- grad_expanded[free_in_group[j]] - last_g
        pos <- pos + 1L
      }
    }
  }

  # Ungrouped levels
  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m > 0) {
    if (centered) {
      if (m >= 2) {
        last_g <- grad_expanded[ungrouped_idx[m]]
        for (j in seq_len(m - 1)) {
          grad_free[pos] <- grad_expanded[ungrouped_idx[j]] - last_g
          pos <- pos + 1L
        }
      }
      # m == 1: constrained to 0, no free params
    } else {
      for (j in seq_len(m)) {
        grad_free[pos] <- grad_expanded[ungrouped_idx[j]]
        pos <- pos + 1L
      }
    }
  }

  grad_free
}

build_param_sizes <- function(config) {
  n_steps <- max(config$n_cat - 1, 0)
  sizes <- list(
    theta = if (config$method == "JMLE") config$theta_spec$n_params else 0
  )
  for (facet in config$facet_names) {
    sizes[[facet]] <- config$facet_specs[[facet]]$n_params
  }
  interaction_n <- sum(vapply(config$interaction_specs %||% list(), function(spec) {
    as.integer(spec$n_params %||% 0L)
  }, integer(1)))
  if (interaction_n > 0L) {
    sizes$interactions <- interaction_n
  }
  if (config$model == "RSM") {
    sizes$steps <- sum_zero_param_count(n_steps)
  } else {
    if (is.null(config$step_facet) || !config$step_facet %in% config$facet_names) {
      stop("PCM model requires 'step_facet' to name one of the facet columns: ",
           paste(config$facet_names, collapse = ", "), ". ",
           "Supply step_facet = '<name>'.", call. = FALSE)
    }
    sizes$steps <- length(config$facet_levels[[config$step_facet]]) *
      sum_zero_param_count(n_steps)
  }
  if (identical(config$model, "GPCM")) {
    gpcm_spec <- config$gpcm_spec %||% list()
    sizes$log_slopes <- as.integer(gpcm_spec$n_params %||% 0L)
  }
  if (identical(config$method, "MML") && isTRUE(config$population_spec$active)) {
    design_matrix <- config$population_spec$design_matrix
    if (is.null(design_matrix) || !is.matrix(design_matrix)) {
      stop("Active latent-regression MML requires a numeric person-level design matrix.",
           call. = FALSE)
    }
    sizes$beta <- ncol(design_matrix)
    sizes$log_sigma2 <- 1L
  }
  sizes
}

split_params <- function(par, sizes) {
  out <- list()
  idx <- 1
  for (nm in names(sizes)) {
    k <- sizes[[nm]]
    if (k == 0) {
      out[[nm]] <- numeric(0)
    } else {
      out[[nm]] <- par[idx:(idx + k - 1)]
      idx <- idx + k
    }
  }
  out
}

center_step_matrix_rows <- function(step_mat) {
  if (!is.matrix(step_mat)) {
    step_mat <- as.matrix(step_mat)
  }
  if (nrow(step_mat) == 0 || ncol(step_mat) == 0) {
    return(step_mat)
  }
  centered_rows <- lapply(
    seq_len(nrow(step_mat)),
    function(i) center_sum_zero(step_mat[i, ])
  )
  matrix(
    unlist(centered_rows, use.names = FALSE),
    nrow = nrow(step_mat),
    ncol = ncol(step_mat),
    byrow = TRUE
  )
}

expand_params <- function(par, sizes, config) {
  parts <- split_params(par, sizes)
  theta <- if (config$method == "JMLE") {
    expand_facet_with_constraints(parts$theta, config$theta_spec)
  } else {
    numeric(0)
  }

  facets <- lapply(config$facet_names, function(facet) {
    expand_facet_with_constraints(parts[[facet]], config$facet_specs[[facet]])
  })
  names(facets) <- config$facet_names
  interactions <- expand_interaction_params(
    free = parts$interactions %||% numeric(0),
    specs = config$interaction_specs %||% list()
  )

  if (config$model == "RSM") {
    steps <- expand_sum_zero_vector(parts$steps, config$n_cat - 1L)
    steps_mat <- NULL
  } else {
    n_levels <- length(config$facet_levels[[config$step_facet]])
    if (n_levels == 0 || config$n_cat <= 1) {
      steps_mat <- matrix(0, nrow = n_levels, ncol = max(config$n_cat - 1, 0))
    } else {
      steps_mat <- expand_step_matrix(
        free = parts$steps,
        n_levels = n_levels,
        n_steps = config$n_cat - 1L
      )
    }
    steps <- NULL
  }
  log_slopes <- NULL
  slopes <- NULL
  if (identical(config$model, "GPCM")) {
    gpcm_expanded <- expand_gpcm_log_slopes(parts$log_slopes, config$gpcm_spec %||% list())
    log_slopes <- gpcm_expanded$log_slopes
    slopes <- gpcm_expanded$slopes
    if (!all(is.finite(slopes)) || any(slopes <= 0)) {
      stop("Expanded GPCM slopes must be finite and strictly positive after log-scale identification.",
           call. = FALSE)
    }
  }
  population <- NULL
  if (identical(config$method, "MML") && isTRUE(config$population_spec$active)) {
    beta <- as.numeric(parts$beta)
    beta_names <- config$population_spec$design_columns %||%
      colnames(config$population_spec$design_matrix) %||%
      paste0("beta_", seq_along(beta))
    names(beta) <- beta_names
    log_sigma2 <- as.numeric(parts$log_sigma2[1] %||% 0)
    population <- list(
      coefficients = beta,
      log_sigma2 = log_sigma2,
      sigma2 = exp(log_sigma2)
    )
  }

  list(
    theta = theta,
    facets = facets,
    interactions = interactions,
    steps = steps,
    steps_mat = steps_mat,
    log_slopes = log_slopes,
    slopes = slopes,
    population = population
  )
}



compute_base_eta <- function(idx, params, config) {
  eta <- rep(0, length(idx$score_k))
  facet_signs <- config$facet_signs
  for (facet in config$facet_names) {
    sign <- if (!is.null(facet_signs) && !is.null(facet_signs[[facet]])) {
      facet_signs[[facet]]
    } else {
      -1
    }
    eta <- eta + sign * params$facets[[facet]][idx$facets[[facet]]]
  }
  eta <- eta + compute_interaction_eta(idx, params, config)
  eta
}

compute_eta <- function(idx, params, config, theta_override = NULL) {
  theta <- if (!is.null(theta_override)) theta_override else params$theta
  eta <- if (length(theta) == 0) {
    rep(0, length(idx$score_k))
  } else {
    theta[idx$person]
  }
  eta + compute_base_eta(idx, params, config)
}

# Joint Maximum Likelihood Estimation (JMLE) negative log-likelihood.
# Estimates person abilities and facet parameters simultaneously.
# Returns -LL for minimization by optim().
mfrm_loglik_jmle <- function(par, idx, config, sizes) {
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    ll <- loglik_rsm(eta, idx$score_k, step_cum, weight = idx$weight)
  } else if (identical(config$model, "GPCM")) {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    ll <- loglik_gpcm(
      eta = eta,
      score_k = idx$score_k,
      step_cum_mat = step_cum_mat,
      criterion_idx = idx$step_idx,
      slopes = params$slopes,
      slope_idx = idx$slope_idx,
      weight = idx$weight
    )
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    ll <- loglik_pcm(eta, idx$score_k, step_cum_mat, idx$step_idx, weight = idx$weight,
                     criterion_splits = idx$criterion_splits)
  }
  -ll
}

# Cached version of JMLE log-likelihood (uses pre-computed params/eta/step_cum)
mfrm_loglik_jmle_cached <- function(cache, idx, config) {
  eta <- cache$eta()
  step_cum <- cache$step_cum()
  if (config$model == "RSM") {
    ll <- loglik_rsm(eta, idx$score_k, step_cum, weight = idx$weight)
  } else if (identical(config$model, "GPCM")) {
    ll <- loglik_gpcm(
      eta = eta,
      score_k = idx$score_k,
      step_cum_mat = step_cum,
      criterion_idx = idx$step_idx,
      slopes = cache$params()$slopes,
      slope_idx = idx$slope_idx,
      weight = idx$weight
    )
  } else {
    ll <- loglik_pcm(eta, idx$score_k, step_cum, idx$step_idx, weight = idx$weight,
                     criterion_splits = idx$criterion_splits)
  }
  -ll
}

mfrm_mml_logprob_bundle_r <- function(idx,
                                      config,
                                      quad,
                                      params,
                                      base_eta,
                                      step_cum = NULL,
                                      include_probs = FALSE,
                                      include_linear_part = FALSE) {
  n <- length(idx$score_k)
  score_k <- idx$score_k
  person_int <- idx$person
  n_nodes <- length(quad$nodes)
  obs_idx <- cbind(seq_len(n), score_k + 1L)
  prob_list <- if (isTRUE(include_probs)) vector("list", n_nodes) else NULL
  linear_part_list <- if (isTRUE(include_linear_part)) vector("list", n_nodes) else NULL

  if (is.null(step_cum)) {
    if (config$model == "RSM") {
      step_cum <- c(0, cumsum(params$steps))
    } else {
      step_cum <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    }
  }

  quad_basis <- resolve_person_quadrature_basis(
    quad = quad,
    population_spec = materialize_population_spec(config, params),
    person_count = config$n_person
  )
  person_nodes <- quad_basis$nodes

  if (config$model == "RSM") {
    k_cat <- length(step_cum)
    step_cum_row <- matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
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
      if (isTRUE(include_probs)) {
        log_denom_mat <- matrix(log_denom, nrow = n, ncol = k_cat)
        prob_list[[q]] <- exp(log_num - log_denom_mat)
      }
    }
  } else if (identical(config$model, "GPCM")) {
    k_cat <- ncol(step_cum)
    crit <- idx$step_idx
    slope_idx <- idx$slope_idx
    step_cum_obs <- step_cum[crit, , drop = FALSE]
    slope_obs <- matrix(params$slopes[slope_idx], nrow = n, ncol = k_cat)
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
      if (isTRUE(include_probs)) {
        log_denom_mat <- matrix(log_denom, nrow = n, ncol = k_cat)
        prob_list[[q]] <- exp(log_num - log_denom_mat)
      }
      if (isTRUE(include_linear_part)) {
        linear_part_list[[q]] <- linear_part
      }
    }
  } else {
    k_cat <- ncol(step_cum)
    step_cum_obs <- step_cum[idx$step_idx, , drop = FALSE]
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
      if (isTRUE(include_probs)) {
        log_denom_mat <- matrix(log_denom, nrow = n, ncol = k_cat)
        prob_list[[q]] <- exp(log_num - log_denom_mat)
      }
    }
  }

  out <- list(
    log_prob_mat = log_prob_mat,
    quad_basis = quad_basis,
    person_int = person_int
  )
  if (isTRUE(include_probs)) {
    out$prob_list <- prob_list
  }
  if (isTRUE(include_linear_part)) {
    out$linear_part_list <- linear_part_list
  }
  out
}

mfrm_mml_logprob_bundle_cpp11 <- function(idx,
                                          config,
                                          quad,
                                          params,
                                          base_eta,
                                          step_cum = NULL,
                                          include_probs = FALSE) {
  n <- length(idx$score_k)
  person_int <- idx$person

  if (is.null(step_cum)) {
    if (config$model == "RSM") {
      step_cum <- c(0, cumsum(params$steps))
    } else {
      step_cum <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    }
  }

  quad_basis <- resolve_person_quadrature_basis(
    quad = quad,
    population_spec = materialize_population_spec(config, params),
    person_count = config$n_person
  )
  person_nodes <- quad_basis$nodes
  weight <- idx$weight %||% NULL

  compiled <- if (config$model == "RSM") {
    mfrm_cpp_rsm_logprob_bundle(
      score_k = as.integer(idx$score_k),
      person_int = as.integer(person_int),
      base_eta = as.numeric(base_eta),
      person_nodes = as.matrix(person_nodes),
      step_cum = as.numeric(step_cum),
      weight = weight,
      include_probs = isTRUE(include_probs)
    )
  } else {
    step_cum_obs <- step_cum[idx$step_idx, , drop = FALSE]
    mfrm_cpp_pcm_logprob_bundle(
      score_k = as.integer(idx$score_k),
      person_int = as.integer(person_int),
      base_eta = as.numeric(base_eta),
      person_nodes = as.matrix(person_nodes),
      step_cum_obs = as.matrix(step_cum_obs),
      weight = weight,
      include_probs = isTRUE(include_probs)
    )
  }

  out <- list(
    log_prob_mat = compiled$log_prob_mat,
    quad_basis = quad_basis,
    person_int = person_int
  )
  if (isTRUE(include_probs)) {
    out$prob_list <- compiled$prob_list
  }
  out
}

mfrm_mml_logprob_bundle <- function(idx,
                                    config,
                                    quad,
                                    params,
                                    base_eta,
                                    step_cum = NULL,
                                    include_probs = FALSE,
                                    include_linear_part = FALSE) {
  if (mfrm_use_cpp11_backend(config, include_linear_part = include_linear_part)) {
    return(mfrm_mml_logprob_bundle_cpp11(
      idx = idx,
      config = config,
      quad = quad,
      params = params,
      base_eta = base_eta,
      step_cum = step_cum,
      include_probs = include_probs
    ))
  }

  mfrm_mml_logprob_bundle_r(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta,
    step_cum = step_cum,
    include_probs = include_probs,
    include_linear_part = include_linear_part
  )
}

mfrm_mml_person_bundle <- function(log_prob_mat,
                                   person_int,
                                   quad_basis,
                                   include_posterior = FALSE) {
  ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)
  person_ids <- as.integer(rownames(ll_by_person))
  log_joint <- quad_basis$log_weights[person_ids, , drop = FALSE] + ll_by_person
  row_max <- row_max_fast(log_joint)
  log_marginal <- row_max + log(rowSums(exp(log_joint - row_max)))

  out <- list(
    ll_by_person = ll_by_person,
    person_ids = person_ids,
    log_joint = log_joint,
    log_marginal = log_marginal
  )

  if (isTRUE(include_posterior)) {
    log_posterior <- log_joint - log_marginal
    out$log_posterior <- log_posterior
    out$posterior <- exp(log_posterior)
  }

  out
}

mfrm_mml_posterior_bundle <- function(logprob_bundle) {
  person_bundle <- mfrm_mml_person_bundle(
    log_prob_mat = logprob_bundle$log_prob_mat,
    person_int = logprob_bundle$person_int,
    quad_basis = logprob_bundle$quad_basis,
    include_posterior = TRUE
  )
  person_to_row <- integer(nrow(logprob_bundle$quad_basis$nodes))
  person_to_row[person_bundle$person_ids] <- seq_along(person_bundle$person_ids)
  obs_person_row <- person_to_row[logprob_bundle$person_int]

  list(
    person_bundle = person_bundle,
    obs_person_row = obs_person_row,
    obs_posterior = person_bundle$posterior[obs_person_row, , drop = FALSE]
  )
}

mfrm_mml_expected_category_bundle_r <- function(logprob_bundle,
                                                posterior_bundle,
                                                include_p_geq = FALSE) {
  prob_list <- logprob_bundle$prob_list
  if (is.null(prob_list) || length(prob_list) == 0L) {
    stop("Expected-category quantities require `mfrm_mml_logprob_bundle(..., include_probs = TRUE)`.",
         call. = FALSE)
  }

  obs_posterior <- posterior_bundle$obs_posterior
  n <- nrow(obs_posterior)
  k_cat <- ncol(prob_list[[1]])
  posterior_prob <- matrix(0, nrow = n, ncol = k_cat)

  for (q in seq_along(prob_list)) {
    posterior_prob <- posterior_prob + prob_list[[q]] * obs_posterior[, q]
  }

  k_vals <- 0:(k_cat - 1)
  expected_k <- as.vector(posterior_prob %*% k_vals)
  var_k <- as.vector(posterior_prob %*% (k_vals^2)) - expected_k^2

  out <- list(
    posterior_prob = posterior_prob,
    expected_k = expected_k,
    var_k = var_k
  )
  if (isTRUE(include_p_geq) && k_cat > 1L) {
    out$p_geq <- compute_P_geq(posterior_prob)
  }
  out
}

mfrm_mml_expected_category_bundle <- function(logprob_bundle,
                                              posterior_bundle,
                                              include_p_geq = FALSE) {
  prob_list <- logprob_bundle$prob_list
  if (is.null(prob_list) || length(prob_list) == 0L) {
    stop("Expected-category quantities require `mfrm_mml_logprob_bundle(..., include_probs = TRUE)`.",
         call. = FALSE)
  }

  if (mfrm_cpp11_backend_available()) {
    return(mfrm_cpp_expected_category_bundle(
      prob_list = prob_list,
      obs_posterior = posterior_bundle$obs_posterior,
      include_p_geq = isTRUE(include_p_geq)
    ))
  }

  mfrm_mml_expected_category_bundle_r(
    logprob_bundle = logprob_bundle,
    posterior_bundle = posterior_bundle,
    include_p_geq = include_p_geq
  )
}

# Cached version of MML log-likelihood (uses pre-computed params/base_eta/step_cum)
mfrm_loglik_mml_cached <- function(cache, idx, config, quad) {
  params <- cache$params()
  base_eta <- cache$base_eta()
  step_cum <- cache$step_cum()
  n <- length(idx$score_k)
  if (n == 0) return(0)

  logprob_bundle <- mfrm_mml_logprob_bundle(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta,
    step_cum = step_cum
  )
  person_bundle <- mfrm_mml_person_bundle(
    log_prob_mat = logprob_bundle$log_prob_mat,
    person_int = logprob_bundle$person_int,
    quad_basis = logprob_bundle$quad_basis
  )
  -sum(person_bundle$log_marginal)
}

# Cached version of JMLE gradient (uses pre-computed params/eta/step_cum)
mfrm_grad_jmle_cached <- function(cache, idx, config, sizes) {
  mfrm_grad_jmle_core(cache$params(), cache$eta(), idx, config, sizes)
}

# Cached version of MML gradient (uses pre-computed params/base_eta/step_cum)
mfrm_grad_mml_cached <- function(cache, idx, config, sizes, quad) {
  mfrm_grad_mml_core(
    cache$params(),
    cache$base_eta(),
    idx,
    config,
    sizes,
    quad,
    step_cum = cache$step_cum()
  )
}

# ---- Analytical gradient for JMLE ----
# Returns gradient of -LL (negative log-likelihood) w.r.t. free parameters.
# Key derivation:
#   d(LL)/d(eta_i)      = k_i - E[k|eta_i]  (observed minus expected score)
#   d(LL)/d(delta_s)    = P(X >= s) - I(X >= s)  (for step parameters)
# Chain rule maps these through constraint Jacobians to free parameters.
mfrm_grad_jmle <- function(par, idx, config, sizes) {
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)
  mfrm_grad_jmle_core(params, eta, idx, config, sizes)
}

mfrm_grad_jmle_core <- function(params, eta, idx, config, sizes) {
  n <- length(eta)
  score_k <- idx$score_k
  weight <- idx$weight

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
    k_cat <- ncol(probs)
    n_steps <- k_cat - 1

    # Residuals: observed - expected
    expected <- as.vector(probs %*% (0:(k_cat - 1)))
    residual <- score_k - expected
    if (!is.null(weight)) residual <- residual * weight

    # Gradient w.r.t. expanded theta
    n_theta <- length(params$theta)
    grad_theta_exp <- numeric(n_theta)
    if (n_theta > 0) {
      rs <- rowsum(matrix(residual, ncol = 1), idx$person)
      grad_theta_exp[as.integer(rownames(rs))] <- as.vector(rs)
    }

    # Gradient w.r.t. expanded facets
    grad_facets_exp <- list()
    for (facet in config$facet_names) {
      sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
      n_lev <- length(params$facets[[facet]])
      g <- numeric(n_lev)
      rs <- rowsum(matrix(sign_f * residual, ncol = 1), idx$facets[[facet]])
      g[as.integer(rownames(rs))] <- as.vector(rs)
      grad_facets_exp[[facet]] <- g
    }
    grad_interaction_free <- compute_interaction_gradient_free(residual, idx, config)

    # Gradient w.r.t. centered step parameters
    P_geq <- compute_P_geq(probs)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_resid <- P_geq - I_geq
    if (!is.null(weight)) step_resid <- step_resid * weight
    grad_step_centered <- colSums(step_resid)

    # Map to the identified step parameters: last centered step is
    # constrained to minus the sum of the preceding free steps.
    grad_step_free <- project_sum_zero_gradient(grad_step_centered)
    grad_log_slope_free <- numeric(0)

  } else if (identical(config$model, "GPCM")) {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_gpcm(
      eta = eta,
      step_cum_mat = step_cum_mat,
      criterion_idx = idx$step_idx,
      slopes = params$slopes,
      slope_idx = idx$slope_idx
    )
    k_cat <- ncol(probs)
    n_steps <- k_cat - 1
    n_criteria <- nrow(params$steps_mat)
    k_vals <- 0:(k_cat - 1)
    slope_obs <- params$slopes[idx$slope_idx]
    linear_part <- outer(eta, k_vals) - step_cum_mat[idx$step_idx, , drop = FALSE]

    expected <- as.vector(probs %*% k_vals)
    residual <- slope_obs * (score_k - expected)
    if (!is.null(weight)) residual <- residual * weight

    n_theta <- length(params$theta)
    grad_theta_exp <- numeric(n_theta)
    if (n_theta > 0) {
      rs <- rowsum(matrix(residual, ncol = 1), idx$person)
      grad_theta_exp[as.integer(rownames(rs))] <- as.vector(rs)
    }

    grad_facets_exp <- list()
    for (facet in config$facet_names) {
      sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
      n_lev <- length(params$facets[[facet]])
      g <- numeric(n_lev)
      rs <- rowsum(matrix(sign_f * residual, ncol = 1), idx$facets[[facet]])
      g[as.integer(rownames(rs))] <- as.vector(rs)
      grad_facets_exp[[facet]] <- g
    }
    grad_interaction_free <- compute_interaction_gradient_free(residual, idx, config)

    P_geq <- compute_P_geq(probs)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_resid <- (P_geq - I_geq) * slope_obs
    if (!is.null(weight)) step_resid <- step_resid * weight
    grad_step_mat <- matrix(0, n_criteria, n_steps)
    rs_step <- rowsum(step_resid, idx$step_idx)
    rs_ids <- as.integer(rownames(rs_step))
    grad_step_mat[rs_ids, ] <- rs_step
    grad_step_free <- project_step_matrix_gradient(grad_step_mat)

    obs_linear <- linear_part[cbind(seq_len(n), score_k + 1L)]
    expected_linear <- rowSums(probs * linear_part)
    slope_score <- slope_obs * (obs_linear - expected_linear)
    if (!is.null(weight)) slope_score <- slope_score * weight
    grad_log_slope_exp <- numeric(length(params$slopes))
    rs_slope <- rowsum(matrix(slope_score, ncol = 1), idx$slope_idx)
    slope_ids <- as.integer(rownames(rs_slope))
    grad_log_slope_exp[slope_ids] <- as.vector(rs_slope)
    grad_log_slope_free <- project_sum_zero_gradient(grad_log_slope_exp)
  } else {
    # PCM
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                               criterion_splits = idx$criterion_splits)
    k_cat <- ncol(probs)
    n_steps <- k_cat - 1
    n_criteria <- nrow(params$steps_mat)

    expected <- as.vector(probs %*% (0:(k_cat - 1)))
    residual <- score_k - expected
    if (!is.null(weight)) residual <- residual * weight

    n_theta <- length(params$theta)
    grad_theta_exp <- numeric(n_theta)
    if (n_theta > 0) {
      rs <- rowsum(matrix(residual, ncol = 1), idx$person)
      grad_theta_exp[as.integer(rownames(rs))] <- as.vector(rs)
    }

    grad_facets_exp <- list()
    for (facet in config$facet_names) {
      sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
      n_lev <- length(params$facets[[facet]])
      g <- numeric(n_lev)
      rs <- rowsum(matrix(sign_f * residual, ncol = 1), idx$facets[[facet]])
      g[as.integer(rownames(rs))] <- as.vector(rs)
      grad_facets_exp[[facet]] <- g
    }
    grad_interaction_free <- compute_interaction_gradient_free(residual, idx, config)

    # Step gradients per criterion
    P_geq <- compute_P_geq(probs)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_resid <- P_geq - I_geq
    if (!is.null(weight)) step_resid <- step_resid * weight

    # Vectorized step gradient aggregation via rowsum
    grad_step_mat <- matrix(0, n_criteria, n_steps)
    rs_step <- rowsum(step_resid, idx$step_idx)
    rs_ids <- as.integer(rownames(rs_step))
    grad_step_mat[rs_ids, ] <- rs_step

    grad_step_free <- project_step_matrix_gradient(grad_step_mat)
    grad_log_slope_free <- numeric(0)
  }

  # Project through constraints
  grad_theta_free <- constraint_grad_project(grad_theta_exp, config$theta_spec)
  grad_facet_free <- unlist(lapply(config$facet_names, function(f) {
    constraint_grad_project(grad_facets_exp[[f]], config$facet_specs[[f]])
  }))

  # Return negative gradient (minimizing -LL)
  -c(grad_theta_free, grad_facet_free, grad_interaction_free,
     grad_step_free, grad_log_slope_free)
}

# Marginal Maximum Likelihood (MML) negative log-likelihood.
# Integrates over the person ability distribution using Gauss-Hermite quadrature:
#   L = prod_p integral P(X_p | theta, facets) phi(theta) d(theta)
#     ~ prod_p sum_q w_q P(X_p | theta_q, facets)
# Person parameters are not estimated directly; instead, EAP post-hoc.
mfrm_loglik_mml <- function(par, idx, config, sizes, quad) {
  params <- expand_params(par, sizes, config)
  n <- length(idx$score_k)
  if (n == 0) return(0)

  base_eta <- compute_base_eta(idx, params, config)
  logprob_bundle <- mfrm_mml_logprob_bundle(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta
  )
  person_bundle <- mfrm_mml_person_bundle(
    log_prob_mat = logprob_bundle$log_prob_mat,
    person_int = logprob_bundle$person_int,
    quad_basis = logprob_bundle$quad_basis
  )

  -sum(person_bundle$log_marginal)
}

# ---- Analytical gradient for MML ----
# Returns gradient of -LL w.r.t. free parameters under marginal ML.
# Key formula:
#   d(LL)/d(param) = sum_p sum_q posterior_pq * d(ll_pq)/d(param)
# where posterior_pq = w_q * L_p(theta_q) / sum_q' w_q' * L_p(theta_q')
# and d(ll_pq)/d(param) uses the same residual/step derivatives as JMLE
# but evaluated at each quadrature node theta_q.
mfrm_grad_mml <- function(par, idx, config, sizes, quad) {
  params <- expand_params(par, sizes, config)
  base_eta <- compute_base_eta(idx, params, config)
  mfrm_grad_mml_core(params, base_eta, idx, config, sizes, quad)
}

mfrm_grad_mml_core <- function(params, base_eta, idx, config, sizes, quad, step_cum = NULL) {
  n <- length(idx$score_k)
  if (n == 0) return(rep(0, sum(unlist(sizes))))
  score_k <- idx$score_k
  weight <- idx$weight
  population_spec <- materialize_population_spec(config, params)
  pop_active <- isTRUE(population_spec$active)
  person_ids_template <- seq_len(config$n_person)
  if (pop_active) {
    sigma_pop <- sqrt(as.numeric(population_spec$sigma2[1]))
    design_rows <- population_spec$person_lookup[person_ids_template]
    design_matrix <- population_spec$design_matrix[design_rows, , drop = FALSE]
    grad_beta <- numeric(ncol(design_matrix))
    grad_log_sigma2 <- 0
  } else {
    sigma_pop <- NA_real_
    design_matrix <- NULL
  }
  logprob_bundle <- mfrm_mml_logprob_bundle(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta,
    step_cum = step_cum,
    include_probs = TRUE,
    include_linear_part = identical(config$model, "GPCM")
  )
  posterior_bundle <- mfrm_mml_posterior_bundle(logprob_bundle)
  person_bundle <- posterior_bundle$person_bundle
  posterior <- person_bundle$posterior
  obs_posterior <- posterior_bundle$obs_posterior
  person_ids <- person_bundle$person_ids
  n_nodes <- ncol(posterior)
  design_matrix_person <- if (pop_active) {
    design_matrix[person_ids, , drop = FALSE]
  } else {
    NULL
  }
  grad_facets_exp <- lapply(config$facet_names, function(f) numeric(length(params$facets[[f]])))
  names(grad_facets_exp) <- config$facet_names
  grad_interaction_free <- rep(
    0,
    sum(vapply(config$interaction_specs %||% list(), function(spec) {
      as.integer(spec$n_params %||% 0L)
    }, integer(1)))
  )

  if (config$model == "RSM") {
    k_cat <- ncol(logprob_bundle$prob_list[[1]])
    n_steps <- k_cat - 1
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    k_vals <- 0:(k_cat - 1)
    grad_step_centered <- numeric(n_steps)

    for (q in seq_len(n_nodes)) {
      probs_q <- logprob_bundle$prob_list[[q]]
      expected_q <- as.vector(probs_q %*% k_vals)
      residual_q <- score_k - expected_q
      if (!is.null(weight)) residual_q <- residual_q * weight

      if (pop_active) {
        rs_person <- rowsum(matrix(residual_q, ncol = 1), logprob_bundle$person_int, reorder = FALSE)
        person_score_q <- as.numeric(rs_person[, 1])
        grad_beta <- grad_beta +
          as.numeric(crossprod(design_matrix_person, posterior[, q] * person_score_q))
        grad_log_sigma2 <- grad_log_sigma2 +
          sum(posterior[, q] * person_score_q * quad$nodes[q]) * (0.5 * sigma_pop)
      }

      obs_post_q <- obs_posterior[, q]
      w_residual <- residual_q * obs_post_q

      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]], reorder = FALSE)
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }
      grad_interaction_free <- grad_interaction_free +
        compute_interaction_gradient_free(w_residual, idx, config)

      P_geq <- compute_P_geq(probs_q)
      step_resid <- (P_geq - I_geq) * obs_post_q
      if (!is.null(weight)) step_resid <- step_resid * weight
      grad_step_centered <- grad_step_centered + colSums(step_resid)
    }

    grad_step_free <- project_sum_zero_gradient(grad_step_centered)
    grad_log_slope_free <- numeric(0)

  } else if (identical(config$model, "GPCM")) {
    k_cat <- ncol(logprob_bundle$prob_list[[1]])
    n_steps <- k_cat - 1
    n_criteria <- nrow(params$steps_mat)
    k_vals <- 0:(k_cat - 1)
    slope_idx_int <- idx$slope_idx
    slope_obs <- params$slopes[slope_idx_int]
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    grad_step_mat <- matrix(0, n_criteria, n_steps)
    grad_log_slope_exp <- numeric(length(params$slopes))
    obs_idx <- cbind(seq_len(n), score_k + 1L)

    for (q in seq_len(n_nodes)) {
      probs_q <- logprob_bundle$prob_list[[q]]
      linear_part_q <- logprob_bundle$linear_part_list[[q]]
      expected_q <- as.vector(probs_q %*% k_vals)
      residual_q <- slope_obs * (score_k - expected_q)
      if (!is.null(weight)) residual_q <- residual_q * weight

      if (pop_active) {
        rs_person <- rowsum(matrix(residual_q, ncol = 1), logprob_bundle$person_int, reorder = FALSE)
        person_score_q <- as.numeric(rs_person[, 1])
        grad_beta <- grad_beta +
          as.numeric(crossprod(design_matrix_person, posterior[, q] * person_score_q))
        grad_log_sigma2 <- grad_log_sigma2 +
          sum(posterior[, q] * person_score_q * quad$nodes[q]) * (0.5 * sigma_pop)
      }

      obs_post_q <- obs_posterior[, q]
      w_residual <- residual_q * obs_post_q

      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]], reorder = FALSE)
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }
      grad_interaction_free <- grad_interaction_free +
        compute_interaction_gradient_free(w_residual, idx, config)

      step_resid <- (compute_P_geq(probs_q) - I_geq) * slope_obs * obs_post_q
      if (!is.null(weight)) step_resid <- step_resid * weight
      rs_step <- rowsum(step_resid, idx$step_idx, reorder = FALSE)
      rs_ids <- as.integer(rownames(rs_step))
      grad_step_mat[rs_ids, ] <- grad_step_mat[rs_ids, ] + rs_step

      obs_linear <- linear_part_q[obs_idx]
      expected_linear <- rowSums(probs_q * linear_part_q)
      slope_score <- slope_obs * (obs_linear - expected_linear) * obs_post_q
      if (!is.null(weight)) slope_score <- slope_score * weight
      rs_slope <- rowsum(matrix(slope_score, ncol = 1), slope_idx_int, reorder = FALSE)
      slope_ids <- as.integer(rownames(rs_slope))
      grad_log_slope_exp[slope_ids] <- grad_log_slope_exp[slope_ids] + as.vector(rs_slope)
    }

    grad_step_free <- project_step_matrix_gradient(grad_step_mat)
    grad_log_slope_free <- project_sum_zero_gradient(grad_log_slope_exp)
  } else {
    k_cat <- ncol(logprob_bundle$prob_list[[1]])
    n_steps <- k_cat - 1
    n_criteria <- nrow(params$steps_mat)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_idx_int <- idx$step_idx
    k_vals <- 0:(k_cat - 1)
    grad_step_mat <- matrix(0, n_criteria, n_steps)

    for (q in seq_len(n_nodes)) {
      probs_q <- logprob_bundle$prob_list[[q]]
      expected_q <- as.vector(probs_q %*% k_vals)
      residual_q <- score_k - expected_q
      if (!is.null(weight)) residual_q <- residual_q * weight

      if (pop_active) {
        rs_person <- rowsum(matrix(residual_q, ncol = 1), logprob_bundle$person_int, reorder = FALSE)
        person_score_q <- as.numeric(rs_person[, 1])
        grad_beta <- grad_beta +
          as.numeric(crossprod(design_matrix_person, posterior[, q] * person_score_q))
        grad_log_sigma2 <- grad_log_sigma2 +
          sum(posterior[, q] * person_score_q * quad$nodes[q]) * (0.5 * sigma_pop)
      }

      obs_post_q <- obs_posterior[, q]
      w_residual <- residual_q * obs_post_q

      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]], reorder = FALSE)
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }
      grad_interaction_free <- grad_interaction_free +
        compute_interaction_gradient_free(w_residual, idx, config)

      P_geq <- compute_P_geq(probs_q)
      step_resid <- (P_geq - I_geq) * obs_post_q
      if (!is.null(weight)) step_resid <- step_resid * weight
      rs_step <- rowsum(step_resid, step_idx_int, reorder = FALSE)
      rs_ids <- as.integer(rownames(rs_step))
      grad_step_mat[rs_ids, ] <- grad_step_mat[rs_ids, ] + rs_step
    }

    grad_step_free <- project_step_matrix_gradient(grad_step_mat)
    grad_log_slope_free <- numeric(0)
  }

  # Project facet gradients through constraints
  grad_facet_free <- unlist(lapply(config$facet_names, function(f) {
    constraint_grad_project(grad_facets_exp[[f]], config$facet_specs[[f]])
  }))

  # Return negative gradient (minimizing -LL); no theta in MML
  if (pop_active) {
    return(-c(grad_facet_free, grad_interaction_free, grad_step_free,
              grad_log_slope_free, grad_beta, grad_log_sigma2))
  }
  -c(grad_facet_free, grad_interaction_free, grad_step_free, grad_log_slope_free)
}

# Expected A Posteriori (EAP) person ability estimates under MML.
# For each person p:
#   theta_EAP = sum_q theta_q * w_q * L(X_p|theta_q) / sum_q w_q * L(X_p|theta_q)
#   SD_EAP   = sqrt(E[theta^2] - (E[theta])^2)  (posterior standard deviation)
compute_person_eap <- function(idx, config, params, quad) {
  n <- length(idx$score_k)
  if (n == 0) {
    return(tibble(Person = character(0), Estimate = numeric(0), SD = numeric(0)))
  }
  base_eta <- compute_base_eta(idx, params, config)
  logprob_bundle <- mfrm_mml_logprob_bundle(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta
  )
  person_bundle <- mfrm_mml_person_bundle(
    log_prob_mat = logprob_bundle$log_prob_mat,
    person_int = logprob_bundle$person_int,
    quad_basis = logprob_bundle$quad_basis,
    include_posterior = TRUE
  )

  # EAP and SD
  nodes_mat <- logprob_bundle$quad_basis$nodes[person_bundle$person_ids, , drop = FALSE]
  eap <- rowSums(nodes_mat * person_bundle$posterior)
  sd_eap <- sqrt(rowSums((nodes_mat - eap)^2 * person_bundle$posterior))

  tibble(Estimate = eap, SD = sd_eap)
}

prepare_constraint_specs <- function(prep,
                                     anchor_df = NULL,
                                     group_anchor_df = NULL,
                                     noncenter_facet = "Person",
                                     dummy_facets = character(0)) {
  facet_names <- prep$facet_names
  all_facets <- c("Person", facet_names)
  anchor_review <- audit_anchor_tables(
    prep = prep,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  anchor_df <- anchor_review$anchors
  group_anchor_df <- anchor_review$group_anchors

  anchor_map <- setNames(vector("list", length(all_facets)), all_facets)
  group_map <- setNames(vector("list", length(all_facets)), all_facets)
  group_values <- setNames(vector("list", length(all_facets)), all_facets)

  for (facet in all_facets) {
    levels <- prep$levels[[facet]]
    if (!is.null(levels)) {
      if (facet %in% dummy_facets) {
        anchor_map[[facet]] <- setNames(rep(0, length(levels)), levels)
        group_map[facet] <- list(NULL)
        group_values[[facet]] <- numeric(0)
        next
      }
      if (nrow(anchor_df) > 0) {
        df <- filter(anchor_df, Facet == facet)
        if (nrow(df) > 0) {
          anchors <- setNames(df$Anchor, df$Level)
          anchors <- anchors[names(anchors) %in% levels]
          if (length(anchors) > 0) anchor_map[[facet]] <- anchors
        }
      }

      if (nrow(group_anchor_df) > 0) {
        df <- filter(group_anchor_df, Facet == facet)
        if (nrow(df) > 0) {
          groups <- setNames(df$Group, df$Level)
          groups <- groups[names(groups) %in% levels]
          if (length(groups) > 0) group_map[[facet]] <- groups

          group_vals <- df |>
            select(Group, GroupValue) |>
            distinct() |>
            filter(!is.na(Group)) |>
            group_by(Group) |>
            summarize(
              GroupValue = {
                vals <- GroupValue[!is.na(GroupValue)]
                if (length(vals) == 0) NA_real_ else vals[1]
              },
              .groups = "drop"
            ) |>
            mutate(GroupValue = ifelse(is.na(GroupValue), 0, GroupValue))
          if (nrow(group_vals) > 0) {
            group_values[[facet]] <- setNames(group_vals$GroupValue, group_vals$Group)
          }
        }
      }
    }
  }

  theta_spec <- build_facet_constraint(
    levels = prep$levels$Person,
    anchors = anchor_map$Person,
    groups = group_map$Person,
    group_values = group_values$Person,
    centered = !(identical(noncenter_facet, "Person"))
  )

  facet_specs <- lapply(facet_names, function(facet) {
    build_facet_constraint(
      levels = prep$levels[[facet]],
      anchors = anchor_map[[facet]],
      groups = group_map[[facet]],
      group_values = group_values[[facet]],
      centered = !identical(noncenter_facet, facet)
    )
  })
  names(facet_specs) <- facet_names

  anchor_summary <- tibble(
    Facet = all_facets,
    AnchoredLevels = vapply(anchor_map, function(x) length(x), integer(1)),
    GroupAnchors = vapply(group_map, function(x) length(unique(x)), integer(1)),
    DummyFacet = all_facets %in% dummy_facets
  )

  list(
    theta_spec = theta_spec,
    facet_specs = facet_specs,
    anchor_summary = anchor_summary,
    anchor_review = anchor_review
  )
}


resolve_step_and_slope_facets <- function(model,
                                          step_facet,
                                          slope_facet,
                                          facet_names) {
  if (identical(model, "RSM")) {
    return(list(
      step_facet = NULL,
      slope_facet = NULL
    ))
  }

  if (identical(model, "GPCM") && is.null(step_facet)) {
    stop("The current bounded `GPCM` branch requires an explicit `step_facet`.",
         call. = FALSE)
  }

  resolved_step <- if (is.null(step_facet)) facet_names[1] else as.character(step_facet[1])
  if (!resolved_step %in% facet_names) {
    stop("step_facet = '", resolved_step, "' is not among the declared facets: ",
         paste(facet_names, collapse = ", "), ". ",
         "Supply a valid facet name.", call. = FALSE)
  }

  if (identical(model, "PCM")) {
    return(list(
      step_facet = resolved_step,
      slope_facet = NULL
    ))
  }

  resolved_slope <- if (is.null(slope_facet)) resolved_step else as.character(slope_facet[1])
  if (!resolved_slope %in% facet_names) {
    stop("slope_facet = '", resolved_slope, "' is not among the declared facets: ",
         paste(facet_names, collapse = ", "), ". ",
         "Supply a valid facet name.", call. = FALSE)
  }
  if (!identical(resolved_step, resolved_slope)) {
    stop("The current bounded `GPCM` branch requires `slope_facet == step_facet`.",
         call. = FALSE)
  }

  list(
    step_facet = resolved_step,
    slope_facet = resolved_slope
  )
}

resolve_pcm_step_facet <- function(model, step_facet, facet_names) {
  resolve_step_and_slope_facets(
    model = model,
    step_facet = step_facet,
    slope_facet = NULL,
    facet_names = facet_names
  )$step_facet
}

stop_if_first_release_gpcm_downstream <- function(fit,
                                                  helper,
                                                  supported = c(
                                                    "fitting",
                                                    "core summary output",
                                                    "fixed-calibration posterior scoring"
                                                  )) {
  model <- as.character(fit$config$model %||% NA_character_)
  if (!identical(model, "GPCM")) {
    return(invisible(NULL))
  }

  helper <- as.character(helper[1] %||% "This helper")
  supported <- as.character(supported %||% character(0))
  supported <- supported[!is.na(supported) & nzchar(supported)]
  if (length(supported) == 0L) {
    supported_text <- "fitting"
  } else if (length(supported) == 1L) {
    supported_text <- supported
  } else {
    supported_text <- paste0(
      paste(utils::head(supported, -1L), collapse = ", "),
      ", and ",
      utils::tail(supported, 1L)
    )
  }

  stop(
    "`", helper, "` is not yet validated for bounded `GPCM` fits. ",
    "Current source-backed `GPCM` support is limited to ",
    supported_text,
    ". Generalized diagnostics, information, simulation, and reporting remain outside the current validated boundary.",
    call. = FALSE
  )
}

sanitize_noncenter_facet <- function(noncenter_facet, facet_names) {
  if (is.null(noncenter_facet)) return("Person")
  valid_set <- c("Person", facet_names)
  if (length(noncenter_facet) == 1L && noncenter_facet %in% valid_set) {
    return(noncenter_facet)
  }
  warning(
    "`noncenter_facet` value ", shQuote(noncenter_facet),
    " is not one of the available facets (",
    paste(shQuote(valid_set), collapse = ", "),
    "); falling back to 'Person'. Check for typos.",
    call. = FALSE
  )
  "Person"
}

sanitize_dummy_facets <- function(dummy_facets, facet_names) {
  if (is.null(dummy_facets)) return(character(0))
  valid_set <- c("Person", facet_names)
  kept <- intersect(dummy_facets, valid_set)
  dropped <- setdiff(dummy_facets, valid_set)
  if (length(dropped) > 0L) {
    warning(
      "`dummy_facets` contains value(s) not present in the model: ",
      paste(shQuote(dropped), collapse = ", "),
      " (available: ", paste(shQuote(valid_set), collapse = ", "),
      "). Unknown entries are ignored; check for typos.",
      call. = FALSE
    )
  }
  kept
}

build_facet_signs <- function(facet_names, positive_facets = character(0)) {
  positives <- intersect(positive_facets, facet_names)
  dropped <- setdiff(positive_facets, facet_names)
  if (length(dropped) > 0L) {
    warning(
      "`positive_facets` contains value(s) not present in the model: ",
      paste(shQuote(dropped), collapse = ", "),
      " (available: ", paste(shQuote(facet_names), collapse = ", "),
      "). Unknown entries are ignored; check for typos (case-sensitive).",
      call. = FALSE
    )
  }
  signs <- setNames(ifelse(facet_names %in% positives, 1, -1), facet_names)
  list(signs = signs, positive_facets = positives)
}

audit_interaction_orientation <- function(config, interaction_facets) {
  interaction_facets <- unique(as.character(interaction_facets %||% character(0)))
  interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
  if (length(interaction_facets) == 0) {
    return(list(
      table = tibble(),
      mixed_sign = FALSE,
      direction_note = "No interaction facets were supplied for orientation review.",
      recommended_action = "None."
    ))
  }

  sign_map <- config$facet_signs %||% list()
  sign_vals <- vapply(interaction_facets, function(facet) {
    sign_map[[facet]] %||% ifelse(identical(facet, "Person"), 1, -1)
  }, numeric(1))
  positive_facets <- as.character(config$positive_facets %||% character(0))
  orientation_tbl <- tibble(
    Facet = interaction_facets,
    Sign = sign_vals,
    Orientation = ifelse(sign_vals >= 0, "positive", "negative"),
    InPositiveFacets = interaction_facets %in% positive_facets
  )

  mixed_sign <- length(unique(sign_vals[is.finite(sign_vals)])) > 1L
  if (mixed_sign) {
    direction_note <- paste(
      "Selected interaction facets mix positive and negative score orientations.",
      "Report cells as higher-than-expected or lower-than-expected rather than using directional labels tied to a single facet."
    )
    recommended_action <- paste(
      "For clearer interpretation, either align the selected facets through `positive_facets`",
      "or keep the mixed orientation and use neutral higher/lower-than-expected wording in reports."
    )
  } else {
    direction_note <- "Selected interaction facets share the same score orientation."
    recommended_action <- "Directional wording is stable across the selected interaction facets."
  }

  list(
    table = orientation_tbl,
    mixed_sign = mixed_sign,
    direction_note = direction_note,
    recommended_action = recommended_action
  )
}

build_estimation_config <- function(prep,
                                    model,
                                    method,
                                    step_facet,
                                    slope_facet,
                                    interaction_specs = NULL,
                                    weight_col,
                                    facet_signs,
                                    positive_facets,
                                    noncenter_facet,
                                    dummy_facets,
                                    anchor_df,
                                    group_anchor_df,
                                    population = NULL) {
  config <- list(
    model = model,
    method = method,
    n_person = length(prep$levels$Person),
    person_levels = prep$levels$Person,
    n_cat = prep$rating_max - prep$rating_min + 1,
    facet_names = prep$facet_names,
    facet_levels = prep$levels[prep$facet_names],
    step_facet = step_facet,
    slope_facet = slope_facet,
    interaction_specs = interaction_specs %||% list(),
    facet_interactions = names(interaction_specs %||% list()),
    keep_original = isTRUE(prep$keep_original),
    rating_min = prep$rating_min,
    rating_max = prep$rating_max,
    rating_range_source = prep$rating_range_source %||% NA_character_,
    rating_min_source = prep$rating_min_source %||% NA_character_,
    rating_max_source = prep$rating_max_source %||% NA_character_,
    score_map = prep$score_map
  )
  config$weight_col <- if (!is.null(weight_col)) weight_col else NULL
  config$positive_facets <- positive_facets
  config$facet_signs <- facet_signs

  constraint_specs <- prepare_constraint_specs(
    prep = prep,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  config$theta_spec <- constraint_specs$theta_spec
  config$facet_specs <- constraint_specs$facet_specs
  config$noncenter_facet <- noncenter_facet
  config$dummy_facets <- dummy_facets
  config$anchor_summary <- constraint_specs$anchor_summary
  config$anchor_review <- constraint_specs$anchor_review
  config$source_columns <- prep$source_columns
  config$population_spec <- compact_population_spec(population, prep$levels$Person)
  config$gpcm_spec <- if (identical(model, "GPCM")) {
    build_gpcm_slope_spec(
      levels = prep$levels[[slope_facet]],
      slope_facet = slope_facet,
      step_facet = step_facet
    )
  } else {
    NULL
  }

  list(
    config = config,
    sizes = build_param_sizes(config)
  )
}

build_initial_param_vector <- function(config, sizes) {
  n_cat <- config$n_cat
  step_init <- if (n_cat > 1) {
    seq(-1, 1, length.out = n_cat - 1)
  } else {
    numeric(0)
  }

  facet_starts <- unlist(lapply(config$facet_names, function(f) rep(0, sizes[[f]])))
  step_init_free <- initial_sum_zero_free(step_init)
  base <- c(
    rep(0, sizes$theta),
    facet_starts,
    if (facet_interactions_active(config)) {
      rep(0, sizes$interactions %||% 0L)
    } else {
      numeric(0)
    },
    if (config$model == "RSM") {
      step_init_free
    } else {
      rep(step_init_free, length(config$facet_levels[[config$step_facet]]))
    },
    if (identical(config$model, "GPCM")) {
      rep(0, sizes$log_slopes %||% 0L)
    } else {
      numeric(0)
    }
  )
  if (identical(config$method, "MML") && isTRUE(config$population_spec$active)) {
    base <- c(
      base,
      rep(0, sizes$beta),
      0
    )
  }
  base
}

# Parameter cache shared between fn and gr to avoid redundant expand_params /
# compute_eta calls within the same optim iteration (BFGS calls fn then gr
# with identical par).
make_param_cache <- function(sizes, config, idx, is_mml = FALSE) {
  cached_par <- NULL
  cached_params <- NULL
  cached_eta <- NULL      # full eta for JMLE
  cached_base_eta <- NULL  # base_eta for MML
  cached_step_cum <- NULL  # step_cum (RSM) or step_cum_mat (PCM)

  list(
    ensure = function(par) {
      if (identical(par, cached_par)) return(invisible(NULL))
      cached_par <<- par
      cached_params <<- expand_params(par, sizes, config)
      if (is_mml) {
        cached_base_eta <<- compute_base_eta(idx, cached_params, config)
      } else {
        cached_eta <<- compute_eta(idx, cached_params, config)
      }
      if (config$model == "RSM") {
        cached_step_cum <<- c(0, cumsum(cached_params$steps))
      } else {
        cached_step_cum <<- t(apply(cached_params$steps_mat, 1,
                                    function(x) c(0, cumsum(x))))
      }
      invisible(NULL)
    },
    params = function() cached_params,
    eta = function() cached_eta,
    base_eta = function() cached_base_eta,
    step_cum = function() cached_step_cum
  )
}

compute_gradient_metrics <- function(gradient) {
  grad <- as.numeric(gradient %||% numeric(0))
  grad <- grad[is.finite(grad)]
  if (length(grad) == 0L) {
    return(list(
      TerminalGradientSupNorm = NA_real_,
      TerminalGradientRMS = NA_real_
    ))
  }

  list(
    TerminalGradientSupNorm = max(abs(grad)),
    TerminalGradientRMS = sqrt(mean(grad^2))
  )
}

build_optimizer_diagnostics <- function(opt,
                                        gradient = NULL,
                                        reltol = 1e-6,
                                        maxit = NA_integer_,
                                        optimizer_method = "BFGS",
                                        convergence_basis = c("optimizer_gradient",
                                                              "relative_loglik")) {
  convergence_basis <- match.arg(convergence_basis)
  code <- as.integer(opt$convergence %||% NA_integer_)
  counts <- opt$counts %||% c()
  fn_evals <- as.integer(unname(counts[["function"]] %||% NA_integer_))
  gr_evals <- as.integer(unname(counts[["gradient"]] %||% NA_integer_))
  message_raw <- opt$message
  message <- if (length(message_raw) == 0L || is.null(message_raw)) {
    NA_character_
  } else {
    as.character(message_raw[1])
  }

  grad_metrics <- compute_gradient_metrics(gradient)
  grad_tol <- if (is.finite(reltol)) max(1e-4, 10 * reltol) else 1e-4
  if (identical(convergence_basis, "relative_loglik")) {
    reviewable_warning <- FALSE
    if (is.na(code)) {
      status <- "unknown"
      reason <- "missing_code"
      severity <- "review"
      detail <- "EM did not return a convergence code."
    } else if (identical(code, 0L)) {
      status <- "converged"
      reason <- "relative_loglik_tolerance_met"
      severity <- "pass"
      detail <- "EM relative log-likelihood change met the stopping tolerance."
    } else if (identical(code, 1L)) {
      status <- "iteration_limit"
      reason <- "relative_loglik_iteration_limit"
      severity <- "fail"
      detail <- "EM reached the iteration limit before the relative log-likelihood change met the stopping tolerance."
    } else {
      status <- "optimizer_warning"
      reason <- "nonzero_code"
      severity <- "fail"
      detail <- "EM returned a nonzero convergence code that requires follow-up."
    }
  } else {
    reviewable_warning <- !is.na(code) &&
      code != 0L &&
      is.finite(grad_metrics$TerminalGradientSupNorm) &&
      grad_metrics$TerminalGradientSupNorm <= grad_tol
    precision_warning <- is.character(message) &&
      !is.na(message) &&
      grepl("precision loss", message, ignore.case = TRUE)

    if (is.na(code)) {
      status <- "unknown"
      reason <- "missing_code"
      severity <- "review"
      detail <- "Optimizer did not return a convergence code."
    } else if (identical(code, 0L)) {
      status <- "converged"
      reason <- "tolerance_met"
      severity <- "pass"
      detail <- "Optimizer returned convergence code 0."
    } else if (reviewable_warning && precision_warning) {
      status <- "reviewable_warning"
      reason <- "precision_warning_small_gradient"
      severity <- "review"
      detail <- "Optimizer reported precision loss, but the terminal gradient was already within the review tolerance."
    } else if (reviewable_warning && identical(code, 1L)) {
      status <- "reviewable_warning"
      reason <- "iteration_limit_small_gradient"
      severity <- "review"
      detail <- "Optimizer reached the iteration limit, but the terminal gradient was already within the review tolerance."
    } else if (reviewable_warning) {
      status <- "reviewable_warning"
      reason <- "nonzero_code_small_gradient"
      severity <- "review"
      detail <- "Optimizer returned a nonzero code, but the terminal gradient was already within the review tolerance."
    } else if (identical(code, 1L)) {
      status <- "iteration_limit"
      reason <- "iteration_limit_large_gradient"
      severity <- "fail"
      detail <- "Optimizer reached the iteration limit before the terminal gradient became small enough for review-only acceptance."
    } else {
      status <- "optimizer_warning"
      reason <- "nonzero_code"
      severity <- "fail"
      detail <- "Optimizer returned a nonzero convergence code that requires follow-up."
    }
  }

  list(
    OptimizerMethod = as.character(optimizer_method),
    ConvergenceCode = code,
    ConvergenceBasis = as.character(convergence_basis),
    ConvergenceStatus = status,
    ConvergenceReason = reason,
    ConvergenceSeverity = severity,
    ConvergenceMessage = message,
    ConvergenceDetail = detail,
    ReviewableWarning = reviewable_warning,
    GradientReviewTolerance = grad_tol,
    FunctionEvaluations = fn_evals,
    GradientEvaluations = gr_evals,
    TerminalGradientSupNorm = grad_metrics$TerminalGradientSupNorm,
    TerminalGradientRMS = grad_metrics$TerminalGradientRMS
  )
}

normalize_mml_engine <- function(mml_engine = "direct") {
  match.arg(
    tolower(as.character(mml_engine %||% "direct")[1]),
    c("direct", "em", "hybrid")
  )
}

resolve_mml_engine_plan <- function(method,
                                    model,
                                    requested = "direct",
                                    population_active = FALSE,
                                    interaction_active = FALSE) {
  requested <- normalize_mml_engine(requested)

  if (!identical(method, "MML")) {
    return(list(
      Requested = requested,
      Used = "not_applicable",
      Fallback = FALSE,
      Detail = "MML engine selection applies only to method = 'MML'."
    ))
  }

  if (identical(requested, "direct")) {
    return(list(
      Requested = requested,
      Used = requested,
      Fallback = FALSE,
      Detail = "Direct marginal likelihood optimization."
    ))
  }

  if (identical(model, "GPCM")) {
    return(list(
      Requested = requested,
      Used = "direct",
      Fallback = TRUE,
      Detail = "EM and hybrid MML are currently implemented only for RSM/PCM; falling back to direct optimization for GPCM."
    ))
  }

  if (isTRUE(interaction_active)) {
    return(list(
      Requested = requested,
      Used = "direct",
      Fallback = TRUE,
      Detail = "EM and hybrid MML are currently implemented only for additive RSM/PCM; falling back to direct optimization for model-estimated facet interactions."
    ))
  }

  if (isTRUE(population_active)) {
    return(list(
      Requested = requested,
      Used = "direct",
      Fallback = TRUE,
      Detail = "EM and hybrid MML currently require `population = NULL`; falling back to direct optimization for the latent-regression branch."
    ))
  }

  list(
    Requested = requested,
    Used = requested,
    Fallback = FALSE,
    Detail = if (identical(requested, "em")) {
      "Expectation-maximization over the marginal likelihood."
    } else {
      "EM warm start followed by direct marginal likelihood optimization."
    }
  )
}

compute_hybrid_em_maxit <- function(maxit) {
  maxit <- as.integer(maxit %||% 0L)
  max(2L, min(3L, maxit %/% 10L))
}

compute_hybrid_em_mstep_maxit <- function(em_maxit) {
  em_maxit <- as.integer(em_maxit %||% 0L)
  max(2L, min(4L, em_maxit + 1L))
}

compute_hybrid_em_reltol <- function(reltol) {
  reltol <- as.numeric(reltol %||% 1e-6)
  max(1e-3, sqrt(reltol))
}


build_person_table <- function(method, idx, config, params, prep, quad_points) {
  # Flag extreme persons (all observed at rating_min or rating_max).
  # Under JMLE the corresponding theta diverges to +-Inf, so surfacing
  # the flag early lets downstream diagnostics and reporting avoid
  # presenting these runaway estimates as substantively meaningful.
  # See Wright (1998) on the JML extreme-score problem.
  extreme_flag <- .flag_extreme_persons(prep, rating_min = prep$rating_min,
                                        rating_max = prep$rating_max)

  if (method == "MML") {
    quad <- gauss_hermite_normal(quad_points)
    tbl <- compute_person_eap(idx, config, params, quad) |>
      mutate(Person = prep$levels$Person) |>
      select(Person, Estimate, SD)
    # The `SD` column is actually the posterior standard deviation of
    # theta under the Gauss-Hermite prior. Expose it under the more
    # precise name `PosteriorSD` while keeping `SD` as an alias for
    # backward compatibility. `SE` mirrors PosteriorSD so callers that
    # reach for the Rasch-convention SE label find it too.
    tbl$PosteriorSD <- tbl$SD
    tbl$SE <- tbl$SD
    tbl$Extreme <- extreme_flag[match(tbl$Person, names(extreme_flag))]
    return(tbl)
  }

  tbl <- tibble(
    Person = prep$levels$Person,
    Estimate = params$theta
  )
  # JML does not produce a posterior SD; per-person SE comes from the
  # observation information inverted in diagnose_mfrm(). Populate a NA
  # SE column so the schema matches MML; callers that need JML SEs
  # should consult `diagnose_mfrm(fit)$measures`.
  tbl$SE <- NA_real_
  tbl$Extreme <- extreme_flag[match(tbl$Person, names(extreme_flag))]
  tbl
}

#' @keywords internal
#' @noRd
.flag_extreme_persons <- function(prep, rating_min, rating_max) {
  if (is.null(prep$data) || nrow(prep$data) == 0L) {
    return(setNames(character(0), character(0)))
  }
  by_person <- split(prep$data$Score, prep$data$Person)
  flag <- vapply(by_person, function(s) {
    s <- s[is.finite(s)]
    if (length(s) == 0L) return(NA_character_)
    if (all(s == rating_max, na.rm = TRUE)) return("high")
    if (all(s == rating_min, na.rm = TRUE)) return("low")
    "none"
  }, character(1))
  flag
}

build_other_facet_table <- function(config, prep, params) {
  facet_tbls <- lapply(config$facet_names, function(facet) {
    tibble(Level = prep$levels[[facet]], Estimate = params$facets[[facet]]) |>
      mutate(Facet = facet, .before = 1)
  })
  bind_rows(facet_tbls)
}

build_step_table <- function(config, prep, params) {
  if (config$model == "RSM") {
    return(
      tibble(
        Step = paste0("Step_", seq_len(config$n_cat - 1)),
        Estimate = params$steps
      )
    )
  }

  expand_grid(
    StepFacet = prep$levels[[config$step_facet]],
    Step = paste0("Step_", seq_len(config$n_cat - 1))
  ) |>
    mutate(Estimate = as.vector(t(params$steps_mat)))
}

build_slope_table <- function(config, prep, params) {
  if (!identical(config$model, "GPCM")) {
    return(tibble())
  }

  tibble(
    SlopeFacet = prep$levels[[config$slope_facet]],
    LogEstimate = params$log_slopes,
    Estimate = params$slopes
  )
}

build_estimation_summary <- function(model, method, prep, config, sizes, opt) {
  k_params <- sum(unlist(sizes))
  loglik <- -opt$value
  n_obs <- if (!is.null(config$weight_col) && "Weight" %in% names(prep$data)) {
    sum(prep$data$Weight, na.rm = TRUE)
  } else {
    nrow(prep$data)
  }
  aic <- 2 * k_params - 2 * loglik
  bic <- log(n_obs) * k_params - 2 * loglik
  interaction_specs <- config$interaction_specs %||% list()
  interaction_count <- length(interaction_specs)
  interaction_params <- sum(vapply(interaction_specs, function(spec) {
    as.integer(spec$n_params %||% 0L)
  }, integer(1)))
  interaction_cells <- sum(vapply(interaction_specs, function(spec) {
    as.integer(spec$n_cells %||% 0L)
  }, integer(1)))
  interaction_sparse_cells <- sum(vapply(interaction_specs, function(spec) {
    as.integer(spec$sparse_count %||% 0L)
  }, integer(1)))

  # Sample-size flag: worst Linacre band across all non-person facet levels.
  # Bands are adapted from Linacre (1994, RMT 7:4): the 30-level band
  # preserves Linacre's approximately +-1.0 logit at 95% CI line, while
  # the `< 10 sparse` floor and the `< 50 standard` watermark are
  # mfrmr-specific screening choices below Linacre's 30-examinee minimum.
  # See facet_small_sample_review() for per-level documentation.
  facet_sample_flag <- NA_character_
  facet_sample_min_n <- NA_integer_
  facet_sample_sparse_count <- 0L
  if (!is.null(config$facet_names) && length(config$facet_names) > 0L &&
      !is.null(prep$data) && nrow(prep$data) > 0L) {
    facet_min_ns <- vapply(config$facet_names, function(f) {
      if (!f %in% names(prep$data)) return(NA_integer_)
      counts <- tabulate(as.integer(factor(prep$data[[f]])))
      if (length(counts) == 0L) NA_integer_ else as.integer(min(counts))
    }, integer(1))
    facet_sample_min_n <- suppressWarnings(min(facet_min_ns, na.rm = TRUE))
    if (!is.finite(facet_sample_min_n)) facet_sample_min_n <- NA_integer_
    if (is.finite(facet_sample_min_n)) {
      facet_sample_flag <- if (facet_sample_min_n < 10L) "sparse"
      else if (facet_sample_min_n < 30L) "marginal"
      else if (facet_sample_min_n < 50L) "standard"
      else "strong"
      facet_sample_sparse_count <- sum(facet_min_ns < 10L, na.rm = TRUE)
    }
  }
  optimizer_diag <- opt$optimizer_diagnostics %||%
    build_optimizer_diagnostics(opt = opt)
  mml_engine <- opt$mml_engine %||% list()
  mml_used <- if (identical(method, "MML")) {
    as.character(mml_engine$Used %||% config$estimation_control$mml_engine_used %||% "direct")
  } else {
    NA_character_
  }
  iterations <- if (identical(method, "MML") && identical(mml_used, "em")) {
    as.integer(mml_engine$EMIterations %||% optimizer_diag$FunctionEvaluations)
  } else {
    optimizer_diag$FunctionEvaluations
  }
  iterations_basis <- if (identical(method, "MML") && identical(mml_used, "em")) {
    "em_iterations"
  } else {
    "function_evaluations"
  }

  tibble(
    Model = model,
    Method = public_mfrm_method_label(method),
    MethodUsed = method,
    N = n_obs,
    Persons = config$n_person,
    Facets = length(config$facet_names),
    FacetInteractions = interaction_count,
    InteractionParameters = interaction_params,
    InteractionCells = interaction_cells,
    InteractionSparseCells = interaction_sparse_cells,
    Categories = config$n_cat,
    LogLik = loglik,
    AIC = aic,
    BIC = bic,
    Converged = opt$convergence == 0,
    Iterations = iterations,
    IterationsBasis = iterations_basis,
    MMLEngineRequested = if (identical(method, "MML")) {
      as.character(mml_engine$Requested %||% config$estimation_control$mml_engine_requested %||% "direct")
    } else {
      NA_character_
    },
    MMLEngineUsed = mml_used,
    MMLEngineDetail = if (identical(method, "MML")) {
      as.character(mml_engine$Detail %||% NA_character_)
    } else {
      NA_character_
    },
    EMIterations = if (identical(method, "MML")) {
      as.integer(mml_engine$EMIterations %||% NA_integer_)
    } else {
      NA_integer_
    },
    EMConverged = if (identical(method, "MML")) {
      as.logical(mml_engine$EMConverged %||% NA)
    } else {
      NA
    },
    EMRelativeChange = if (identical(method, "MML")) {
      as.numeric(mml_engine$EMRelativeChange %||% NA_real_)
    } else {
      NA_real_
    },
    OptimizerMethod = optimizer_diag$OptimizerMethod,
    ConvergenceCode = optimizer_diag$ConvergenceCode,
    ConvergenceBasis = optimizer_diag$ConvergenceBasis,
    ConvergenceStatus = optimizer_diag$ConvergenceStatus,
    ConvergenceReason = optimizer_diag$ConvergenceReason,
    ConvergenceSeverity = optimizer_diag$ConvergenceSeverity,
    ConvergenceMessage = optimizer_diag$ConvergenceMessage,
    ConvergenceDetail = optimizer_diag$ConvergenceDetail,
    ReviewableWarning = optimizer_diag$ReviewableWarning,
    GradientReviewTolerance = optimizer_diag$GradientReviewTolerance,
    FunctionEvaluations = optimizer_diag$FunctionEvaluations,
    GradientEvaluations = optimizer_diag$GradientEvaluations,
    TerminalGradientSupNorm = optimizer_diag$TerminalGradientSupNorm,
    TerminalGradientRMS = optimizer_diag$TerminalGradientRMS,
    # Summary of facet sample-size adequacy. See facet_small_sample_review()
    # for a full per-level report.
    FacetSampleSizeFlag = facet_sample_flag,
    FacetMinLevelN = facet_sample_min_n,
    FacetSparseCount = facet_sample_sparse_count,
    # Extreme-person counts (0.1.6 polish). Under JMLE, extreme persons
    # produce +-Inf theta; MML returns a finite EAP but the information
    # is small. See `fit$facets$person$Extreme` for the per-person flag.
    ExtremeHighN = .count_extreme_persons(prep, which = "high"),
    ExtremeLowN = .count_extreme_persons(prep, which = "low")
  )
}

#' @keywords internal
#' @noRd
.count_extreme_persons <- function(prep, which = c("high", "low")) {
  which <- match.arg(which)
  flag <- .flag_extreme_persons(prep,
                                rating_min = prep$rating_min,
                                rating_max = prep$rating_max)
  as.integer(sum(flag == which, na.rm = TRUE))
}

# ---- estimation wrapper ----
mfrm_estimate <- function(data, person_col, facet_cols, score_col,
                          rating_min = NULL, rating_max = NULL,
                          weight_col = NULL, keep_original = FALSE,
                          missing_codes = NULL,
                          model = c("RSM", "PCM", "GPCM"), method = c("JMLE", "MML"),
                          step_facet = NULL,
                          slope_facet = NULL,
                          facet_interactions = NULL,
                          min_obs_per_interaction = 10,
                          interaction_policy = c("warn", "error", "silent"),
                          anchor_df = NULL,
                          group_anchor_df = NULL,
                          noncenter_facet = "Person",
                          dummy_facets = character(0),
                          positive_facets = character(0),
                          population = NULL,
                          quad_points = 15, maxit = 400, reltol = 1e-6,
                          mml_engine = "direct",
                          checkpoint = NULL) {
  # Stage 1: Normalize model options and input data.
  model <- match.arg(model)
  method <- match.arg(method)
  interaction_policy <- match.arg(interaction_policy)

  prep <- prepare_mfrm_data(
    data,
    person_col = person_col,
    facet_cols = facet_cols,
    score_col = score_col,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight_col,
    keep_original = keep_original,
    missing_codes = missing_codes
  )

  # Stage 2: Resolve facet-level modeling choices.
  resolved_families <- resolve_step_and_slope_facets(
    model = model,
    step_facet = step_facet,
    slope_facet = slope_facet,
    facet_names = prep$facet_names
  )
  step_facet <- resolved_families$step_facet
  slope_facet <- resolved_families$slope_facet
  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)
  sign_info <- build_facet_signs(prep$facet_names, positive_facets)
  if (!is.null(facet_interactions) && length(facet_interactions) > 0L &&
      identical(model, "GPCM")) {
    stop(
      "`facet_interactions` currently supports only `model = \"RSM\"` or ",
      "`model = \"PCM\"`. Re-fit without interactions or use a Rasch-family ",
      "model before extending the bounded `GPCM` branch.",
      call. = FALSE
    )
  }
  interaction_specs <- build_facet_interaction_specs(
    prep = prep,
    facet_interactions = facet_interactions,
    min_obs_per_interaction = min_obs_per_interaction,
    interaction_policy = interaction_policy
  )

  # Stage 3: Build reusable structures for optimization.
  idx <- build_indices(
    prep,
    step_facet = step_facet,
    slope_facet = slope_facet,
    interaction_specs = interaction_specs
  )
  cfg <- build_estimation_config(
    prep = prep,
    model = model,
    method = method,
    step_facet = step_facet,
    slope_facet = slope_facet,
    interaction_specs = interaction_specs,
    weight_col = weight_col,
    facet_signs = sign_info$signs,
    positive_facets = sign_info$positive_facets,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df,
    population = population
  )
  config <- cfg$config
  config$estimation_control <- list(
    maxit = as.integer(maxit),
    reltol = as.numeric(reltol),
    quad_points = as.integer(quad_points),
    mml_engine_requested = normalize_mml_engine(mml_engine),
    checkpoint = checkpoint
  )
  sizes <- cfg$sizes
  start <- build_initial_param_vector(config, sizes)

  # Stage 4: Optimize model parameters.
  opt <- run_mfrm_optimization(
    start = start,
    method = method,
    idx = idx,
    config = config,
    sizes = sizes,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol,
    checkpoint = config$estimation_control$checkpoint
  )
  config$estimation_control$mml_engine_used <- as.character(opt$mml_engine$Used %||% NA_character_)
  config$estimation_control$mml_engine_detail <- as.character(opt$mml_engine$Detail %||% NA_character_)

  # Stage 5: Build human-readable output tables.
  params <- expand_params(opt$par, sizes, config)
  person_tbl <- build_person_table(method, idx, config, params, prep, quad_points)
  facet_tbl <- build_other_facet_table(config, prep, params)
  interaction_tbl <- build_interaction_effect_table(config, prep, params)
  step_tbl <- build_step_table(config, prep, params)
  slope_tbl <- build_slope_table(config, prep, params)
  summary_tbl <- build_estimation_summary(model, method, prep, config, sizes, opt)

  list(
    summary = summary_tbl,
    facets = list(
      person = person_tbl,
      others = facet_tbl
    ),
    interactions = list(
      effects = interaction_tbl,
      specs = config$interaction_specs
    ),
    steps = step_tbl,
    slopes = slope_tbl,
    config = config,
    prep = prep,
    opt = opt
  )
}

expected_score_table <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)

  prob_bundle <- compute_response_probability_bundle(config, idx, params, eta)
  tibble(
    Observed = prep$data$Score,
    Expected = prep$rating_min + prob_bundle$expected_k
  )
}

compute_obs_table <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  person_levels <- prep$levels$Person
  person_measure <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate[match(person_levels, res$facets$person$Person)]
  }
  person_measure_by_row <- person_measure[idx$person]
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  prob_bundle <- compute_response_probability_bundle(config, idx, params, eta)
  resid_k <- idx$score_k - prob_bundle$expected_k
  std_sq <- resid_k^2 / prob_bundle$var_k

  # Per-observation categorical-distribution moments needed by
  # `compute_person_fit_indices()` to compute Drasgow et al. (1985) lz
  # in its proper polytomous form. `pr_obs` is the model probability of
  # the observed category; `item_entropy` and `item_var_logp` are the
  # per-item expectation and variance of log P(X|theta) under the model
  # categorical, which are the centering and per-item-variance terms in
  # Drasgow's lz.
  probs_mat <- prob_bundle$probs
  n_rows <- nrow(probs_mat)
  if (n_rows > 0L) {
    score_k_int <- as.integer(idx$score_k)
    obs_col <- pmin(pmax(score_k_int + 1L, 1L), ncol(probs_mat))
    pr_obs <- probs_mat[cbind(seq_len(n_rows), obs_col)]
    eps <- .Machine$double.eps
    log_probs_mat <- log(pmax(probs_mat, eps))
    item_entropy <- as.numeric(rowSums(probs_mat * log_probs_mat))
    item_e_logp_sq <- as.numeric(rowSums(probs_mat * log_probs_mat^2))
    item_var_logp <- pmax(item_e_logp_sq - item_entropy^2, 0)
    k_mat <- matrix(
      seq_len(ncol(probs_mat)) - 1L,
      nrow = n_rows,
      ncol = ncol(probs_mat),
      byrow = TRUE
    )
    score_deriv_mat <- sweep(k_mat, 1L, prob_bundle$expected_k, FUN = "-") *
      matrix(prob_bundle$slope_obs, nrow = n_rows, ncol = ncol(probs_mat))
    item_logp_score_cov <- as.numeric(rowSums(probs_mat * log_probs_mat * score_deriv_mat))
    observed_score_deriv <- prob_bundle$slope_obs * (idx$score_k - prob_bundle$expected_k)
  } else {
    pr_obs <- numeric(0)
    item_entropy <- numeric(0)
    item_var_logp <- numeric(0)
    item_logp_score_cov <- numeric(0)
    observed_score_deriv <- numeric(0)
  }

  prep$data |>
    mutate(
      PersonMeasure = person_measure_by_row,
      Observed = prep$rating_min + idx$score_k,
      Expected = prep$rating_min + prob_bundle$expected_k,
      Var = prob_bundle$var_k,
      FourthCentralMoment = prob_bundle$fourth_central_moment,
      ScoreInformation = prob_bundle$score_information,
      ScoreSlope = prob_bundle$slope_obs,
      ObservedScoreDerivative = observed_score_deriv,
      Residual = Observed - Expected,
      StdResidual = Residual / sqrt(Var),
      StdSq = std_sq,
      PrObserved = pr_obs,
      ItemEntropy = item_entropy,
      ItemVarLogP = item_var_logp,
      ItemLogPScoreCov = item_logp_score_cov
    )
}

extract_bias_facet_spec <- function(bias_results, data_cols = NULL) {
  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(NULL)
  }
  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)

  facets <- character(0)
  if (!is.null(bias_results$interaction_facets)) {
    facets <- c(facets, as.character(bias_results$interaction_facets))
  }
  if (!is.null(bias_results$facet_a) && length(bias_results$facet_a) > 0) {
    facets <- c(facets, as.character(bias_results$facet_a[1]))
  }
  if (!is.null(bias_results$facet_b) && length(bias_results$facet_b) > 0) {
    facets <- c(facets, as.character(bias_results$facet_b[1]))
  }
  facets <- facets[!is.na(facets) & nzchar(facets)]
  facets <- unique(facets)

  level_cols <- grep("^Facet[0-9]+_Level$", names(tbl), value = TRUE)
  if (length(level_cols) > 0) {
    ord <- suppressWarnings(as.integer(sub("^Facet([0-9]+)_Level$", "\\1", level_cols)))
    ok <- is.finite(ord)
    level_cols <- level_cols[ok]
    ord <- ord[ok]
    if (length(level_cols) == 0) return(NULL)
    o <- order(ord)
    ord <- ord[o]
    level_cols <- level_cols[o]

    facet_name_cols <- paste0("Facet", ord)
    if (all(facet_name_cols %in% names(tbl)) && nrow(tbl) > 0) {
      facets_from_tbl <- vapply(facet_name_cols, function(col) as.character(tbl[[col]][1]), character(1))
      if (all(!is.na(facets_from_tbl) & nzchar(facets_from_tbl))) {
        facets <- facets_from_tbl
      }
    }
    if (length(facets) < length(level_cols)) {
      return(NULL)
    }
    facets <- facets[seq_along(level_cols)]
    index_cols <- paste0("Facet", ord, "_Index")
    measure_cols <- paste0("Facet", ord, "_Measure")
    se_cols <- paste0("Facet", ord, "_SE")
  } else if (all(c("FacetA_Level", "FacetB_Level") %in% names(tbl))) {
    if (length(facets) < 2 && all(c("FacetA", "FacetB") %in% names(tbl)) && nrow(tbl) > 0) {
      facets <- c(as.character(tbl$FacetA[1]), as.character(tbl$FacetB[1]))
    }
    facets <- facets[!is.na(facets) & nzchar(facets)]
    facets <- unique(facets)
    if (length(facets) < 2) return(NULL)
    facets <- facets[1:2]
    level_cols <- c("FacetA_Level", "FacetB_Level")
    index_cols <- c("FacetA_Index", "FacetB_Index")
    measure_cols <- c("FacetA_Measure", "FacetB_Measure")
    se_cols <- c("FacetA_SE", "FacetB_SE")
  } else {
    return(NULL)
  }

  if (!is.null(data_cols) && !all(facets %in% data_cols)) {
    return(NULL)
  }

  list(
    facets = facets,
    level_cols = level_cols,
    index_cols = index_cols,
    measure_cols = measure_cols,
    se_cols = se_cols,
    interaction_order = length(facets),
    interaction_mode = ifelse(length(facets) > 2, "higher_order", "pairwise")
  )
}

compute_bias_adjustment_vector <- function(res, bias_results = NULL) {
  n <- nrow(res$prep$data)
  if (n == 0) return(numeric(0))
  adj <- rep(0, n)

  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(adj)
  }

  data <- res$prep$data
  spec <- extract_bias_facet_spec(bias_results, data_cols = names(data))
  if (is.null(spec) || length(spec$facets) < 2) {
    return(adj)
  }

  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  req_cols <- c(spec$level_cols, "Bias Size")
  if (!all(req_cols %in% names(tbl))) {
    return(adj)
  }

  level_cols <- spec$level_cols
  tbl <- tbl |>
    mutate(
      across(all_of(level_cols), as.character),
      `Bias Size` = suppressWarnings(as.numeric(.data$`Bias Size`))
    ) |>
    filter(is.finite(.data$`Bias Size`))
  if (nrow(tbl) == 0) return(adj)

  # Duplicate rows can appear after manual edits; average as a safe default.
  lookup <- tbl |>
    group_by(across(all_of(level_cols))) |>
    summarize(Bias = mean(.data$`Bias Size`, na.rm = TRUE), .groups = "drop")
  lookup_key_parts <- lapply(level_cols, function(col) as.character(lookup[[col]]))
  lookup$Key <- do.call(paste, c(lookup_key_parts, sep = "||"))

  data_key_parts <- lapply(spec$facets, function(f) as.character(data[[f]]))
  names(data_key_parts) <- spec$facets
  keys <- do.call(paste, c(data_key_parts, sep = "||"))

  idx_hit <- match(keys, lookup$Key)
  hit <- !is.na(idx_hit)
  if (any(hit)) {
    adj[hit] <- lookup$Bias[idx_hit[hit]]
  }
  adj
}

compute_prob_matrix_with_bias <- function(res, bias_results = NULL) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") params$theta else res$facets$person$Estimate
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  bias_adj <- compute_bias_adjustment_vector(res, bias_results = bias_results)
  if (length(bias_adj) == length(eta)) {
    eta <- eta + bias_adj
  }
  compute_response_probability_bundle(config, idx, params, eta)$probs
}

compute_obs_table_with_bias <- function(res, bias_results = NULL) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  person_levels <- prep$levels$Person
  person_measure <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate[match(person_levels, res$facets$person$Person)]
  }
  person_measure_by_row <- person_measure[idx$person]
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  bias_adj <- compute_bias_adjustment_vector(res, bias_results = bias_results)
  if (length(bias_adj) == length(eta)) {
    eta <- eta + bias_adj
  } else {
    bias_adj <- rep(0, length(eta))
  }
  prob_bundle <- compute_response_probability_bundle(config, idx, params, eta)
  resid_k <- idx$score_k - prob_bundle$expected_k
  std_sq <- resid_k^2 / prob_bundle$var_k

  prep$data |>
    mutate(
      PersonMeasure = person_measure_by_row,
      Observed = prep$rating_min + idx$score_k,
      Expected = prep$rating_min + prob_bundle$expected_k,
      Var = prob_bundle$var_k,
      FourthCentralMoment = prob_bundle$fourth_central_moment,
      ScoreInformation = prob_bundle$score_information,
      ScoreSlope = prob_bundle$slope_obs,
      Residual = Observed - Expected,
      StdResidual = Residual / sqrt(Var),
      StdSq = std_sq,
      BiasAdjustment = bias_adj
    )
}

calc_unexpected_response_table <- function(obs_df,
                                           probs,
                                           facet_names,
                                           rating_min,
                                           abs_z_min = 2,
                                           prob_max = 0.30,
                                           top_n = 100,
                                           rule = c("either", "both")) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  if (is.null(obs_df) || nrow(obs_df) == 0 || is.null(probs) || nrow(probs) != nrow(obs_df)) {
    return(tibble())
  }
  if (!"score_k" %in% names(obs_df) || !"StdResidual" %in% names(obs_df)) {
    return(tibble())
  }

  score_k <- suppressWarnings(as.integer(obs_df$score_k))
  n_cat <- ncol(probs)
  valid <- is.finite(score_k) & score_k >= 0 & score_k < n_cat
  obs_prob <- rep(NA_real_, nrow(obs_df))
  row_idx <- which(valid)
  if (length(row_idx) > 0) {
    obs_prob[row_idx] <- probs[cbind(row_idx, score_k[row_idx] + 1L)]
  }

  max_idx <- max.col(probs, ties.method = "first")
  most_likely_k <- max_idx - 1L
  most_likely <- rating_min + most_likely_k
  most_likely_prob <- probs[cbind(seq_len(nrow(probs)), max_idx)]

  abs_std <- abs(obs_df$StdResidual)
  surprise <- -log10(pmax(obs_prob, .Machine$double.xmin))
  cat_gap <- abs(obs_df$Observed - most_likely)
  low_prob <- is.finite(obs_prob) & obs_prob <= prob_max
  high_resid <- is.finite(abs_std) & abs_std >= abs_z_min
  flagged <- if (rule == "both") {
    low_prob & high_resid
  } else {
    low_prob | high_resid
  }
  if (!any(flagged, na.rm = TRUE)) return(tibble())

  out <- obs_df |>
    mutate(
      Row = dplyr::row_number(),
      ObsProb = obs_prob,
      MostLikely = most_likely,
      MostLikelyProb = most_likely_prob,
      Surprise = surprise,
      AbsStdResidual = abs_std,
      CategoryGap = cat_gap,
      FlagLowProbability = low_prob,
      FlagLargeResidual = high_resid,
      Unexpected = flagged,
      Direction = dplyr::case_when(
        .data$Residual > 0 ~ "Higher than expected",
        .data$Residual < 0 ~ "Lower than expected",
        TRUE ~ "As expected"
      ),
      Severity = .data$AbsStdResidual + .data$Surprise + 0.5 * .data$CategoryGap
    ) |>
    filter(.data$Unexpected)

  id_cols <- c("Person", facet_names)
  if ("Weight" %in% names(obs_df)) id_cols <- c(id_cols, "Weight")
  id_cols <- unique(id_cols)
  keep_cols <- c(
    "Row",
    id_cols,
    "Score",
    "Observed",
    "Expected",
    "Residual",
    "StdResidual",
    "ObsProb",
    "MostLikely",
    "MostLikelyProb",
    "CategoryGap",
    "Surprise",
    "Direction",
    "FlagLowProbability",
    "FlagLargeResidual",
    "Severity"
  )
  keep_cols <- keep_cols[keep_cols %in% names(out)]

  top_n <- max(1L, as.integer(top_n))
  out |>
    arrange(desc(.data$Severity), desc(abs(.data$StdResidual)), .data$ObsProb) |>
    select(dplyr::all_of(keep_cols)) |>
    slice_head(n = top_n)
}

summarize_unexpected_response_table <- function(unexpected_tbl,
                                                total_observations,
                                                abs_z_min = 2,
                                                prob_max = 0.30,
                                                rule = "either") {
  if (is.null(unexpected_tbl) || nrow(unexpected_tbl) == 0) {
    return(tibble(
      TotalObservations = total_observations,
      UnexpectedN = 0L,
      UnexpectedPercent = 0,
      LowProbabilityN = 0L,
      LargeResidualN = 0L,
      Rule = rule,
      AbsZThreshold = abs_z_min,
      ProbThreshold = prob_max
    ))
  }
  low_n <- if ("FlagLowProbability" %in% names(unexpected_tbl)) {
    sum(unexpected_tbl$FlagLowProbability, na.rm = TRUE)
  } else {
    NA_integer_
  }
  resid_n <- if ("FlagLargeResidual" %in% names(unexpected_tbl)) {
    sum(unexpected_tbl$FlagLargeResidual, na.rm = TRUE)
  } else {
    NA_integer_
  }
  tibble(
    TotalObservations = total_observations,
    UnexpectedN = nrow(unexpected_tbl),
    UnexpectedPercent = ifelse(total_observations > 0, 100 * nrow(unexpected_tbl) / total_observations, NA_real_),
    LowProbabilityN = low_n,
    LargeResidualN = resid_n,
    Rule = rule,
    AbsZThreshold = abs_z_min,
    ProbThreshold = prob_max
  )
}

compute_prob_matrix <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  compute_response_probability_bundle(config, idx, params, eta)$probs
}

# Displacement threshold defaults follow the convention documented in
# Linacre (2026), `A User's Guide to Winsteps 5.11.0`, where
# `Displacement > 0.5 logit`
# together with `|t| > 2` is used to flag anchor instability. The combined
# rule keeps the flag rate conservative for small-sample designs where a
# 0.5-logit shift alone may still be within sampling error.
calc_displacement_table <- function(obs_df,
                                    res,
                                    measures = NULL,
                                    abs_displacement_warn = 0.5,
                                    abs_t_warn = 2) {
  if (is.null(obs_df) || nrow(obs_df) == 0 || is.null(res$config)) {
    return(tibble())
  }

  facet_cols <- c("Person", res$config$facet_names)
  obs_df <- obs_df |>
    mutate(.Weight = get_weights(obs_df))
  info_col <- if ("ScoreInformation" %in% names(obs_df)) "ScoreInformation" else "Var"

  displacement <- purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        WeightedN = sum(.data$.Weight, na.rm = TRUE),
        ResidualSum = sum(.data$Residual * .data$.Weight, na.rm = TRUE),
        Information = sum(.data[[info_col]] * .data$.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Facet = facet,
        Level = as.character(.data[[facet]])
      ) |>
      select(
        "Facet",
        "Level",
        "WeightedN",
        "ResidualSum",
        "Information"
      )
  })

  if (nrow(displacement) == 0) return(tibble())

  displacement <- displacement |>
    mutate(
      Displacement = ifelse(.data$Information > 0, .data$ResidualSum / .data$Information, NA_real_),
      DisplacementSE = ifelse(.data$Information > 0, 1 / sqrt(.data$Information), NA_real_),
      DisplacementT = ifelse(
        is.finite(.data$DisplacementSE) & .data$DisplacementSE > 0,
        .data$Displacement / .data$DisplacementSE,
        NA_real_
      )
    )

  if (!is.null(measures) && nrow(measures) > 0) {
    displacement <- displacement |>
      left_join(
        measures |>
          select("Facet", "Level", "Estimate", "SE", "N"),
        by = c("Facet", "Level")
      )
  } else {
    displacement <- displacement |>
      mutate(Estimate = NA_real_, SE = NA_real_, N = .data$WeightedN)
  }

  anchor_tbl <- extract_anchor_tables(res$config)$anchors
  if (is.null(anchor_tbl) || nrow(anchor_tbl) == 0) {
    anchor_tbl <- tibble(Facet = character(0), Level = character(0), AnchorValue = numeric(0))
  } else {
    anchor_tbl <- anchor_tbl |>
      transmute(
        Facet = as.character(.data$Facet),
        Level = as.character(.data$Level),
        AnchorValue = as.numeric(.data$Anchor)
      )
  }

  status_tbl <- purrr::map_dfr(unique(displacement$Facet), function(facet) {
    lv <- displacement |>
      filter(.data$Facet == facet) |>
      pull(.data$Level)
    tibble(
      Facet = facet,
      Level = lv,
      AnchorStatus = facet_anchor_status(facet, lv, res$config)
    )
  })

  displacement |>
    left_join(anchor_tbl, by = c("Facet", "Level")) |>
    left_join(status_tbl, by = c("Facet", "Level")) |>
    mutate(
      AnchorType = case_when(
        .data$AnchorStatus == "A" ~ "Anchor",
        .data$AnchorStatus == "G" ~ "Group",
        TRUE ~ "Free"
      ),
      ReleasedEstimate = ifelse(
        is.finite(.data$Estimate) & is.finite(.data$Displacement),
        .data$Estimate + .data$Displacement,
        NA_real_
      ),
      AnchorGap = ifelse(
        is.finite(.data$AnchorValue) & is.finite(.data$ReleasedEstimate),
        .data$ReleasedEstimate - .data$AnchorValue,
        NA_real_
      ),
      FlagDisplacement = is.finite(.data$Displacement) & abs(.data$Displacement) >= abs_displacement_warn,
      FlagT = is.finite(.data$DisplacementT) & abs(.data$DisplacementT) >= abs_t_warn,
      Flag = .data$FlagDisplacement | .data$FlagT
    ) |>
    arrange(desc(abs(.data$Displacement)), desc(abs(.data$DisplacementT)))
}

summarize_displacement_table <- function(displacement_tbl,
                                         abs_displacement_warn = 0.5,
                                         abs_t_warn = 2) {
  if (is.null(displacement_tbl) || nrow(displacement_tbl) == 0) {
    return(tibble(
      Levels = 0L,
      AnchoredLevels = 0L,
      FlaggedLevels = 0L,
      FlaggedAnchoredLevels = 0L,
      MaxAbsDisplacement = NA_real_,
      MaxAbsDisplacementT = NA_real_,
      AbsDisplacementThreshold = abs_displacement_warn,
      AbsTThreshold = abs_t_warn
    ))
  }
  is_anchored <- displacement_tbl$AnchorType %in% c("Anchor", "Group")
  flagged <- if ("Flag" %in% names(displacement_tbl)) {
    displacement_tbl$Flag
  } else {
    rep(FALSE, nrow(displacement_tbl))
  }
  abs_disp <- abs(displacement_tbl$Displacement)
  abs_t <- abs(displacement_tbl$DisplacementT)
  tibble(
    Levels = nrow(displacement_tbl),
    AnchoredLevels = sum(is_anchored, na.rm = TRUE),
    FlaggedLevels = sum(flagged, na.rm = TRUE),
    FlaggedAnchoredLevels = sum(flagged & is_anchored, na.rm = TRUE),
    MaxAbsDisplacement = if (any(is.finite(abs_disp))) max(abs_disp, na.rm = TRUE) else NA_real_,
    MaxAbsDisplacementT = if (any(is.finite(abs_t))) max(abs_t, na.rm = TRUE) else NA_real_,
    AbsDisplacementThreshold = abs_displacement_warn,
    AbsTThreshold = abs_t_warn
  )
}

compute_scorefile <- function(res) {
  obs <- compute_obs_table(res)
  if (nrow(obs) == 0) return(tibble())
  probs <- compute_prob_matrix(res)
  if (is.null(probs) || nrow(probs) != nrow(obs)) return(obs)
  cat_vals <- seq(res$prep$rating_min, res$prep$rating_max)
  prob_df <- as_tibble(probs, .name_repair = "minimal")
  names(prob_df) <- paste0("P_", cat_vals)
  max_idx <- apply(probs, 1, which.max)
  max_prob <- probs[cbind(seq_len(nrow(probs)), max_idx)]
  most_likely <- cat_vals[max_idx]
  base_cols <- c("Person", res$config$facet_names)
  if ("Weight" %in% names(obs)) {
    base_cols <- c(base_cols, "Weight")
  }
  base_cols <- c(base_cols, "Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "PersonMeasure")
  bind_cols(
    obs |>
      select(all_of(base_cols)),
    prob_df
  ) |>
    mutate(
      MostLikely = most_likely,
      MaxProb = max_prob
    )
}

compute_residual_file <- function(res) {
  obs <- compute_obs_table(res)
  if (nrow(obs) == 0) return(tibble())
  base_cols <- c("Person", res$config$facet_names)
  if ("Weight" %in% names(obs)) {
    base_cols <- c(base_cols, "Weight")
  }
  base_cols <- c(base_cols, "Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "StdSq", "PersonMeasure")
  obs |>
    select(all_of(base_cols))
}

fit_df_method_choices <- c("engine", "facets", "both")

match_fit_df_method <- function(fit_df_method = "engine") {
  match.arg(as.character(fit_df_method %||% "engine"), fit_df_method_choices)
}

fit_zstd_transform_label <- function(whexact = FALSE) {
  if (isTRUE(whexact)) "linear normal approximation" else "Wilson-Hilferty"
}

compute_facets_fit_df <- function(sum_var_w, sum_w, denom_infit, denom_outfit) {
  denom_infit <- suppressWarnings(as.numeric(denom_infit))
  denom_outfit <- suppressWarnings(as.numeric(denom_outfit))
  sum_var_w <- suppressWarnings(as.numeric(sum_var_w))
  sum_w <- suppressWarnings(as.numeric(sum_w))
  df_infit <- if (is.finite(denom_infit) && denom_infit > 1e-12 &&
                   is.finite(sum_var_w) && sum_var_w > 0) {
    2 * sum_var_w^2 / denom_infit
  } else {
    NA_real_
  }
  df_outfit <- if (is.finite(denom_outfit) && denom_outfit > 1e-12 &&
                    is.finite(sum_w) && sum_w > 0) {
    2 * sum_w^2 / denom_outfit
  } else {
    NA_real_
  }
  c(DF_Infit_FACETS = df_infit, DF_Outfit_FACETS = df_outfit)
}

fit_df_terms <- function(var, fourth, weight) {
  var <- suppressWarnings(as.numeric(var))
  fourth <- suppressWarnings(as.numeric(fourth))
  weight <- suppressWarnings(as.numeric(weight))
  if (length(weight) == 0L) weight <- rep(1, length(var))
  safe_var <- ifelse(is.finite(var) & var > 1e-12, var, NA_real_)
  valid <- is.finite(safe_var) & is.finite(fourth) & is.finite(weight) & weight > 0
  denom_infit <- sum(weight[valid] * (fourth[valid] - safe_var[valid]^2), na.rm = TRUE)
  denom_outfit <- sum(
    weight[valid] * (fourth[valid] / (safe_var[valid]^2) - 1),
    na.rm = TRUE
  )
  c(FACETSDenom_Infit = denom_infit, FACETSDenom_Outfit = denom_outfit)
}

apply_fit_df_method <- function(tbl,
                                fit_df_method = "engine",
                                whexact = FALSE,
                                facets_zstd_cap = 9) {
  fit_df_method <- match_fit_df_method(fit_df_method)
  if (!all(c("Infit", "Outfit", "DF_Infit", "DF_Outfit",
             "DF_Infit_FACETS", "DF_Outfit_FACETS") %in% names(tbl))) {
    return(tbl)
  }

  engine_cols <- c(
    DF_Infit_ENGINE = "DF_Infit",
    DF_Outfit_ENGINE = "DF_Outfit",
    InfitZSTD_ENGINE = "InfitZSTD",
    OutfitZSTD_ENGINE = "OutfitZSTD"
  )
  for (nm in names(engine_cols)) {
    src <- unname(engine_cols[nm])
    tbl[[nm]] <- tbl[[src]]
  }

  tbl$InfitZSTD_FACETS <- zstd_from_mnsq_facets(
    tbl$Infit,
    tbl$DF_Infit_FACETS,
    whexact = whexact,
    cap = facets_zstd_cap
  )
  tbl$OutfitZSTD_FACETS <- zstd_from_mnsq_facets(
    tbl$Outfit,
    tbl$DF_Outfit_FACETS,
    whexact = whexact,
    cap = facets_zstd_cap
  )

  if (identical(fit_df_method, "facets")) {
    tbl$DF_Infit <- tbl$DF_Infit_FACETS
    tbl$DF_Outfit <- tbl$DF_Outfit_FACETS
    tbl$InfitZSTD <- tbl$InfitZSTD_FACETS
    tbl$OutfitZSTD <- tbl$OutfitZSTD_FACETS
  }

  tbl$FitDfMethod <- dplyr::case_when(
    identical(fit_df_method, "facets") ~ "facets_wright_masters",
    identical(fit_df_method, "both") ~ "engine_primary_facets_available",
    TRUE ~ "engine"
  )
  tbl$FitZSTDTransform <- fit_zstd_transform_label(whexact)
  tbl$FitZSTDCap <- if (identical(fit_df_method, "engine")) {
    NA_real_
  } else {
    suppressWarnings(as.numeric(facets_zstd_cap))
  }

  if (identical(fit_df_method, "engine")) {
    drop_cols <- c(
      names(engine_cols), "DF_Infit_FACETS", "DF_Outfit_FACETS",
      "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
      "FitDfMethod", "FitZSTDTransform", "FitZSTDCap"
    )
    tbl <- tbl[, setdiff(names(tbl), drop_cols), drop = FALSE]
  }

  tbl
}

# Overall model fit: weighted infit and outfit mean-square statistics.
# Infit (information-weighted): sum(StdSq * Var * w) / sum(Var * w)
#   Sensitive to unexpected responses near the person's ability level.
# Outfit (unweighted): sum(StdSq * w) / sum(w)
#   Sensitive to outlying unexpected responses far from the ability level.
# The primary ZSTD uses the package's engine df convention unless the caller
# requests the FACETS/Wright-Masters df convention explicitly.
calc_overall_fit <- function(obs_df,
                             whexact = FALSE,
                             fit_df_method = "engine",
                             facets_zstd_cap = 9) {
  fit_df_method <- match_fit_df_method(fit_df_method)
  w <- get_weights(obs_df)
  infit <- sum(obs_df$StdSq * obs_df$Var * w, na.rm = TRUE) / sum(obs_df$Var * w, na.rm = TRUE)
  outfit <- sum(obs_df$StdSq * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  df_infit <- sum(obs_df$Var * w, na.rm = TRUE)
  df_outfit <- sum(w, na.rm = TRUE)
  fourth <- if ("FourthCentralMoment" %in% names(obs_df)) {
    obs_df$FourthCentralMoment
  } else {
    rep(NA_real_, nrow(obs_df))
  }
  df_terms <- fit_df_terms(obs_df$Var, fourth, w)
  facets_df <- compute_facets_fit_df(
    sum_var_w = df_infit,
    sum_w = df_outfit,
    denom_infit = df_terms["FACETSDenom_Infit"],
    denom_outfit = df_terms["FACETSDenom_Outfit"]
  )
  out <- tibble(
    Infit = infit,
    Outfit = outfit,
    InfitZSTD = zstd_from_mnsq(infit, df_infit, whexact = whexact),
    OutfitZSTD = zstd_from_mnsq(outfit, df_outfit, whexact = whexact),
    DF_Infit = df_infit,
    DF_Outfit = df_outfit,
    DF_Infit_FACETS = unname(facets_df["DF_Infit_FACETS"]),
    DF_Outfit_FACETS = unname(facets_df["DF_Outfit_FACETS"])
  )
  apply_fit_df_method(out, fit_df_method = fit_df_method, whexact = whexact,
                      facets_zstd_cap = facets_zstd_cap)
}

calc_facet_fit <- function(obs_df,
                           facet_cols,
                           whexact = FALSE,
                           fit_df_method = "engine",
                           facets_zstd_cap = 9) {
  fit_df_method <- match_fit_df_method(fit_df_method)
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  if (!"FourthCentralMoment" %in% names(obs_df)) {
    obs_df$FourthCentralMoment <- NA_real_
  }
  purrr::map_dfr(facet_cols, function(facet) {
    df <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        Infit = sum(StdSq * Var * .Weight, na.rm = TRUE) / sum(Var * .Weight, na.rm = TRUE),
        Outfit = sum(StdSq * .Weight, na.rm = TRUE) / sum(.Weight, na.rm = TRUE),
        DF_Infit = sum(Var * .Weight, na.rm = TRUE),
        DF_Outfit = sum(.Weight, na.rm = TRUE),
        FACETSDenom_Infit = sum(
          .Weight * (.data$FourthCentralMoment - .data$Var^2),
          na.rm = TRUE
        ),
        FACETSDenom_Outfit = sum(
          .Weight * (.data$FourthCentralMoment / pmax(.data$Var, 1e-12)^2 - 1),
          na.rm = TRUE
        ),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      )
    df |>
      mutate(
        DF_Infit_FACETS = ifelse(
          is.finite(.data$FACETSDenom_Infit) & .data$FACETSDenom_Infit > 1e-12 &
            is.finite(.data$DF_Infit) & .data$DF_Infit > 0,
          2 * .data$DF_Infit^2 / .data$FACETSDenom_Infit,
          NA_real_
        ),
        DF_Outfit_FACETS = ifelse(
          is.finite(.data$FACETSDenom_Outfit) & .data$FACETSDenom_Outfit > 1e-12 &
            is.finite(.data$DF_Outfit) & .data$DF_Outfit > 0,
          2 * .data$DF_Outfit^2 / .data$FACETSDenom_Outfit,
          NA_real_
        ),
        InfitZSTD = zstd_from_mnsq(Infit, DF_Infit, whexact = whexact),
        OutfitZSTD = zstd_from_mnsq(Outfit, DF_Outfit, whexact = whexact)
      ) |>
      apply_fit_df_method(fit_df_method = fit_df_method, whexact = whexact,
                          facets_zstd_cap = facets_zstd_cap) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      select(-dplyr::all_of(c(facet, "FACETSDenom_Infit", "FACETSDenom_Outfit"))) |>
      dplyr::relocate(Facet, Level, N, Infit, Outfit, InfitZSTD, OutfitZSTD, DF_Infit, DF_Outfit)
  })
}

# Approximate facet-level SE from summed observation information.
calc_facet_se <- function(obs_df, facet_cols) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  info_col <- if ("ScoreInformation" %in% names(obs_df)) "ScoreInformation" else "Var"
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        Info = sum(.data[[info_col]] * .Weight, na.rm = TRUE),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Facet = facet,
        Level = .data[[facet]],
        SE = ifelse(Info > 0, 1 / sqrt(Info), NA_real_)
      ) |>
      select(Facet, Level, N, SE)
  })
}

calc_bias_facet <- function(obs_df, facet_cols) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        ObservedAverage = weighted_mean(Observed, .Weight),
        ExpectedAverage = weighted_mean(Expected, .Weight),
        MeanResidual = weighted_mean(Residual, .Weight),
        MeanStdResidual = weighted_mean(StdResidual, .Weight),
        MeanAbsStdResidual = weighted_mean(abs(StdResidual), .Weight),
        ChiSq = sum((StdResidual^2) * .Weight, na.rm = TRUE),
        SE_Residual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, sqrt(sum(Var * .Weight, na.rm = TRUE)) / n_w, NA_real_)
        },
        SE_StdResidual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, 1 / sqrt(n_w), NA_real_)
        },
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      mutate(
        Bias = MeanResidual,
        DF = ifelse(N > 1, N - 1, NA_real_),
        t_Residual = ifelse(is.finite(SE_Residual) & SE_Residual > 0, MeanResidual / SE_Residual, NA_real_),
        t_StdResidual = ifelse(is.finite(SE_StdResidual) & SE_StdResidual > 0, MeanStdResidual / SE_StdResidual, NA_real_),
        p_Residual = ifelse(is.finite(DF) & is.finite(t_Residual), 2 * stats::pt(-abs(t_Residual), df = DF), NA_real_),
        p_StdResidual = ifelse(is.finite(DF) & is.finite(t_StdResidual), 2 * stats::pt(-abs(t_StdResidual), df = DF), NA_real_),
        ChiDf = DF,
        ChiP = ifelse(is.finite(ChiSq) & is.finite(ChiDf) & ChiDf > 0, 1 - stats::pchisq(ChiSq, df = ChiDf), NA_real_)
      ) |>
      select(Facet, Level, N, ObservedAverage, ExpectedAverage, Bias, MeanResidual, MeanStdResidual,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP, SE_Residual, t_Residual, p_Residual,
             SE_StdResidual, t_StdResidual, p_StdResidual, DF)
  })
}

calc_bias_interactions <- function(obs_df, facet_cols, pairs = NULL, top_n = 20) {
  if (length(facet_cols) < 2) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  if (is.null(pairs)) {
    combos <- combn(facet_cols, 2, simplify = FALSE)
  } else if (length(pairs) == 0) {
    return(tibble())
  } else {
    combos <- pairs
  }
  out <- purrr::map_dfr(combos, function(pair) {
    pair1 <- pair[1]
    pair2 <- pair[2]
    obs_df |>
      group_by(.data[[pair1]], .data[[pair2]]) |>
      summarize(
        ObservedAverage = weighted_mean(Observed, .Weight),
        ExpectedAverage = weighted_mean(Expected, .Weight),
        MeanResidual = weighted_mean(Residual, .Weight),
        MeanStdResidual = weighted_mean(StdResidual, .Weight),
        MeanAbsStdResidual = weighted_mean(abs(StdResidual), .Weight),
        ChiSq = sum((StdResidual^2) * .Weight, na.rm = TRUE),
        SE_Residual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, sqrt(sum(Var * .Weight, na.rm = TRUE)) / n_w, NA_real_)
        },
        SE_StdResidual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, 1 / sqrt(n_w), NA_real_)
        },
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Pair = paste(pair1, "x", pair2),
        Level = paste(.data[[pair1]], .data[[pair2]], sep = " | ")
      ) |>
      mutate(
        Bias = MeanResidual,
        DF = ifelse(N > 1, N - 1, NA_real_),
        t_Residual = ifelse(is.finite(SE_Residual) & SE_Residual > 0, MeanResidual / SE_Residual, NA_real_),
        t_StdResidual = ifelse(is.finite(SE_StdResidual) & SE_StdResidual > 0, MeanStdResidual / SE_StdResidual, NA_real_),
        p_Residual = ifelse(is.finite(DF) & is.finite(t_Residual), 2 * stats::pt(-abs(t_Residual), df = DF), NA_real_),
        p_StdResidual = ifelse(is.finite(DF) & is.finite(t_StdResidual), 2 * stats::pt(-abs(t_StdResidual), df = DF), NA_real_),
        ChiDf = DF,
        ChiP = ifelse(is.finite(ChiSq) & is.finite(ChiDf) & ChiDf > 0, 1 - stats::pchisq(ChiSq, df = ChiDf), NA_real_)
      ) |>
      select(Pair, Level, N, ObservedAverage, ExpectedAverage, Bias, MeanResidual, MeanStdResidual,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP, SE_Residual, t_Residual, p_Residual,
             SE_StdResidual, t_StdResidual, p_StdResidual, DF)
  })
  out |>
    mutate(AbsStd = abs(MeanStdResidual)) |>
    arrange(desc(AbsStd)) |>
    select(-AbsStd) |>
    slice_head(n = top_n)
}

safe_cor <- function(x, y, w = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (is.null(w)) {
    if (!any(ok)) return(NA_real_)
    x <- x[ok]
    y <- y[ok]
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(NA_real_)
    }
    return(suppressWarnings(stats::cor(x, y, use = "complete.obs")))
  }
  ok <- ok & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  w_sum <- sum(w)
  if (w_sum <= 0) return(NA_real_)
  mx <- sum(w * x) / w_sum
  my <- sum(w * y) / w_sum
  vx <- sum(w * (x - mx)^2) / w_sum
  vy <- sum(w * (y - my)^2) / w_sum
  if (vx <= 0 || vy <= 0) return(NA_real_)
  cov <- sum(w * (x - mx) * (y - my)) / w_sum
  cov / sqrt(vx * vy)
}

weighted_mean_safe <- function(x, w) {
  weighted_mean(x, w)
}

infer_default_rater_facet <- function(facet_names) {
  facet_names <- as.character(facet_names)
  if (length(facet_names) == 0) return(NULL)

  lowered <- tolower(facet_names)
  patterns <- c("rater", "judge", "grader", "reader", "scorer", "assessor", "evaluator")
  hits <- vapply(
    lowered,
    function(x) any(vapply(patterns, function(p) grepl(p, x, fixed = TRUE), logical(1))),
    logical(1)
  )
  if (any(hits)) {
    return(facet_names[which(hits)[1]])
  }
  facet_names[1]
}

calc_interrater_agreement <- function(obs_df, facet_cols, rater_facet, res = NULL) {
  if (is.null(obs_df) || nrow(obs_df) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }
  if (is.null(rater_facet) || !rater_facet %in% facet_cols) {
    return(list(summary = tibble(), pairs = tibble()))
  }
  context_cols <- setdiff(facet_cols, rater_facet)
  if (length(context_cols) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  df <- obs_df |>
    mutate(across(all_of(context_cols), as.character)) |>
    tidyr::unite(".context", all_of(context_cols), sep = "|", remove = FALSE) |>
    select(.context, !!rlang::sym(rater_facet), Observed, any_of("Weight"))
  df$.Weight <- get_weights(df)

  df <- df |>
    group_by(.context, !!rlang::sym(rater_facet)) |>
    summarize(Score = weighted_mean(Observed, .Weight), .groups = "drop")

  if (nrow(df) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  wide <- tryCatch(
    tidyr::pivot_wider(
      df,
      id_cols = .context,
      names_from = !!rlang::sym(rater_facet),
      values_from = Score
    ),
    error = function(e) NULL
  )
  if (is.null(wide)) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  rater_cols <- setdiff(names(wide), ".context")
  if (length(rater_cols) < 2) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  # `prob_map` was historically a `list()` keyed by "ctx||rater" which
  # gave O(N) average-case lookup against tens of thousands of keys.
  # Swap to an `environment` (hash-backed for character keys) so the
  # per-pair loop below is O(N_pairs) rather than O(N_pairs * N_keys).
  prob_map <- new.env(parent = emptyenv(), hash = TRUE)
  prob_map_size <- 0L
  if (!is.null(res)) {
    probs <- compute_prob_matrix(res)
    if (!is.null(probs) && nrow(probs) == nrow(obs_df)) {
      prob_cols <- paste0(".p", seq_len(ncol(probs)))
      prob_df <- obs_df |>
        mutate(across(all_of(context_cols), as.character)) |>
        tidyr::unite(".context", all_of(context_cols), sep = "|", remove = FALSE) |>
        select(.context, !!rlang::sym(rater_facet), any_of("Weight"))
      prob_df[prob_cols] <- probs
      prob_df$.Weight <- get_weights(prob_df)

      prob_avg <- prob_df |>
        group_by(.context, !!rlang::sym(rater_facet)) |>
        summarize(
          across(all_of(prob_cols), ~ weighted_mean(.x, .Weight)),
          .groups = "drop"
        )

      if (nrow(prob_avg) > 0) {
        prob_avg_mat <- as.matrix(prob_avg[, prob_cols, drop = FALSE])
        ctx_vec <- as.character(prob_avg$.context)
        rater_vec <- as.character(prob_avg[[rater_facet]])
        for (i in seq_len(nrow(prob_avg))) {
          key <- paste(ctx_vec[i], rater_vec[i], sep = "||")
          assign(key, prob_avg_mat[i, ], envir = prob_map)
        }
        prob_map_size <- nrow(prob_avg)
      }
    }
  }

  pairs <- combn(rater_cols, 2, simplify = FALSE)
  pair_tbl <- purrr::map_dfr(pairs, function(pair) {
    sub <- wide |>
      select(.context, !!rlang::sym(pair[1]), !!rlang::sym(pair[2])) |>
      filter(is.finite(.data[[pair[1]]]), is.finite(.data[[pair[2]]]))
    n_ok <- nrow(sub)
    if (n_ok == 0) {
      return(tibble(
        Rater1 = pair[1],
        Rater2 = pair[2],
        N = 0,
        OpportunityCount = 0,
        ExactCount = 0,
        ExpectedExactCount = NA_real_,
        AdjacentCount = 0,
        Exact = NA_real_,
        ExpectedExact = NA_real_,
        Adjacent = NA_real_,
        MeanDiff = NA_real_,
        MAD = NA_real_,
        Corr = NA_real_
      ))
    }
    v1 <- sub[[pair[1]]]
    v2 <- sub[[pair[2]]]
    diff <- v1 - v2
    exact_count <- sum(diff == 0, na.rm = TRUE)

    # Preallocate `exp_vals` and fill by index. `c(exp_vals, ...)`
    # inside the loop allocated O(N^2) memory and dominated runtime
    # for diagnose_mfrm() on N >= 50k observations (Round 4 audit
    # measured a 643x slowdown vs the no-Person path).
    contexts <- as.character(sub$.context)
    n_ctx <- length(contexts)
    exp_vals <- if (prob_map_size > 0L && n_ctx > 0L) {
      tmp <- numeric(n_ctx)
      kept <- 0L
      key1_vec <- paste(contexts, pair[1], sep = "||")
      key2_vec <- paste(contexts, pair[2], sep = "||")
      for (k in seq_len(n_ctx)) {
        p1 <- get0(key1_vec[k], envir = prob_map, inherits = FALSE)
        p2 <- get0(key2_vec[k], envir = prob_map, inherits = FALSE)
        if (is.null(p1) || is.null(p2)) next
        if (any(!is.finite(p1)) || any(!is.finite(p2))) next
        kept <- kept + 1L
        tmp[kept] <- sum(p1 * p2)
      }
      if (kept > 0L) tmp[seq_len(kept)] else numeric(0)
    } else {
      numeric(0)
    }
    exp_mean <- if (length(exp_vals) > 0) mean(exp_vals) else NA_real_

    tibble(
      Rater1 = pair[1],
      Rater2 = pair[2],
      N = n_ok,
      OpportunityCount = n_ok,
      ExactCount = exact_count,
      ExpectedExactCount = if (length(exp_vals) > 0) sum(exp_vals) else NA_real_,
      AdjacentCount = sum(abs(diff) <= 1, na.rm = TRUE),
      Exact = exact_count / n_ok,
      ExpectedExact = exp_mean,
      Adjacent = mean(abs(diff) <= 1, na.rm = TRUE),
      MeanDiff = mean(diff, na.rm = TRUE),
      MAD = mean(abs(diff), na.rm = TRUE),
      Corr = safe_cor(v1, v2)
    )
  })

  if (nrow(pair_tbl) == 0L) {
    pair_tbl <- data.frame(
      Rater1 = character(0),
      Rater2 = character(0),
      N = integer(0),
      OpportunityCount = integer(0),
      ExactCount = integer(0),
      ExpectedExactCount = numeric(0),
      AdjacentCount = integer(0),
      Exact = numeric(0),
      ExpectedExact = numeric(0),
      Adjacent = numeric(0),
      MeanDiff = numeric(0),
      MAD = numeric(0),
      Corr = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  contexts_with_pairs <- sum(rowSums(!is.na(wide[rater_cols])) >= 2)
  total_pairs <- sum(pair_tbl$N, na.rm = TRUE)
  total_exact <- sum(pair_tbl$Exact * pair_tbl$N, na.rm = TRUE)
  expected_available <- any(is.finite(pair_tbl$ExpectedExact))
  total_expected <- if (expected_available) {
    sum(pair_tbl$ExpectedExact * pair_tbl$N, na.rm = TRUE)
  } else {
    NA_real_
  }
  summary_tbl <- tibble(
    RaterFacet = rater_facet,
    Raters = length(rater_cols),
    Pairs = nrow(pair_tbl),
    Contexts = contexts_with_pairs,
    TotalPairs = total_pairs,
    OpportunityCount = total_pairs,
    ExactAgreements = total_exact,
    ExpectedAgreements = ifelse(expected_available, total_expected, NA_real_),
    ExactAgreement = ifelse(total_pairs > 0, total_exact / total_pairs, NA_real_),
    ExpectedExactAgreement = ifelse(expected_available && total_pairs > 0, total_expected / total_pairs, NA_real_),
    AgreementMinusExpected = ifelse(
      expected_available && total_pairs > 0,
      (total_exact - total_expected) / total_pairs,
      NA_real_
    ),
    AdjacentAgreements = sum(pair_tbl$AdjacentCount, na.rm = TRUE),
    AdjacentAgreement = weighted_mean_safe(pair_tbl$Adjacent, pair_tbl$N),
    MeanAbsDiff = weighted_mean_safe(pair_tbl$MAD, pair_tbl$N),
    MeanCorr = weighted_mean_safe(pair_tbl$Corr, pair_tbl$N)
  )

  list(summary = summary_tbl, pairs = pair_tbl)
}

augment_interrater_with_precision <- function(agreement, reliability_tbl, rater_facet = NULL) {
  if (is.null(agreement) || !is.list(agreement)) return(agreement)
  summary_tbl <- as.data.frame(agreement$summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_tbl) == 0) return(agreement)

  if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
    if ("RaterFacet" %in% names(summary_tbl)) {
      rater_facet <- as.character(summary_tbl$RaterFacet[1])
    }
  } else {
    rater_facet <- as.character(rater_facet[1])
  }

  rel_tbl <- as.data.frame(reliability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (!"Facet" %in% names(rel_tbl) || !nzchar(rater_facet)) {
    agreement$summary <- summary_tbl
    return(agreement)
  }

  rel_row <- rel_tbl[as.character(rel_tbl$Facet) == rater_facet, , drop = FALSE]
  if (nrow(rel_row) == 0) {
    agreement$summary <- summary_tbl
    return(agreement)
  }

  summary_tbl$RaterSeparation <- suppressWarnings(as.numeric(rel_row$Separation[1]))
  summary_tbl$RaterStrata <- suppressWarnings(as.numeric(rel_row$Strata[1]))
  summary_tbl$RaterReliability <- suppressWarnings(as.numeric(rel_row$Reliability[1]))
  summary_tbl$RaterRealSeparation <- suppressWarnings(as.numeric(rel_row$RealSeparation[1]))
  summary_tbl$RaterRealReliability <- suppressWarnings(as.numeric(rel_row$RealReliability[1]))
  agreement$summary <- summary_tbl
  agreement
}

calc_ptmea <- function(obs_df, facet_cols) {
  facet_cols <- setdiff(facet_cols, "Person")
  if (length(facet_cols) == 0) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        PTMEA = safe_cor(Observed, PersonMeasure, w = .Weight),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      select(Facet, Level, PTMEA, N)
  })
}

expected_score_from_eta <- function(eta, step_cum, rating_min) {
  if (!is.finite(eta) || length(step_cum) == 0) return(NA_real_)
  probs <- category_prob_rsm(eta, step_cum)
  k_vals <- 0:(length(step_cum) - 1)
  rating_min + sum(probs * k_vals)
}

# Slope-aware element-conditional expected score for GPCM fair-averages.
# Computes
#     E[X | eta, a, step_cum] = sum_k k * exp(a * (k * eta - step_cum_k))
#                                / sum_r exp(a * (r * eta - step_cum_r))
# in the log-space-shifted form for numerical stability. Reduces exactly
# to `expected_score_from_eta()` when `slope == 1` (the PCM/RSM special
# case). Degenerate or non-finite slopes return `NA` rather than falling back
# to the PCM/RSM special case, because a silent slope reset can hide a broken
# bounded-GPCM fit object.
expected_score_from_eta_gpcm <- function(eta, step_cum, slope, rating_min) {
  if (!is.finite(eta) || length(step_cum) == 0) return(NA_real_)
  if (!is.finite(slope) || slope <= 0) return(NA_real_)
  k_vals <- 0:(length(step_cum) - 1)
  log_num <- slope * (k_vals * eta - step_cum)
  m <- max(log_num)
  if (!is.finite(m)) return(NA_real_)
  probs <- exp(log_num - m)
  probs <- probs / sum(probs)
  rating_min + sum(probs * k_vals)
}

estimate_eta_from_target <- function(target, step_cum, rating_min, rating_max,
                                     slope = 1) {
  if (!is.finite(target) || length(step_cum) == 0) return(NA_real_)
  if (target <= rating_min) return(-Inf)
  if (target >= rating_max) return(Inf)
  # Uses the slope-aware helper. At slope = 1 (the default and the
  # RSM/PCM case) this is identical to `expected_score_from_eta()` to
  # machine precision (verified in tests/testthat/test-gpcm-fair-
  # average.R), so the existing PCM/RSM xtreme-eta path is unchanged.
  f <- function(eta) expected_score_from_eta_gpcm(eta, step_cum, slope, rating_min) - target
  lower <- -10
  upper <- 10
  f_low <- f(lower)
  f_up <- f(upper)
  if (!is.finite(f_low) || !is.finite(f_up)) return(NA_real_)
  if (f_low * f_up > 0) {
    lower <- -20
    upper <- 20
    f_low <- f(lower)
    f_up <- f(upper)
    if (!is.finite(f_low) || !is.finite(f_up) || f_low * f_up > 0) return(NA_real_)
  }
  uniroot(f, lower = lower, upper = upper)$root
}

facet_anchor_status <- function(facet, levels, config) {
  spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
  if (is.null(spec)) return(rep("", length(levels)))
  anchors <- spec$anchors
  groups <- spec$groups
  status <- rep("", length(levels))
  if (!is.null(anchors)) {
    anchor_vals <- anchors[match(levels, names(anchors))]
    status[is.finite(anchor_vals)] <- "A"
  }
  if (!is.null(groups)) {
    group_vals <- groups[match(levels, names(groups))]
    status[status == "" & !is.na(group_vals) & group_vals != ""] <- "G"
  }
  status
}

calc_facets_report_tbls <- function(res,
                                    diagnostics,
                                    totalscore = TRUE,
                                    umean = 0,
                                    uscale = 1,
                                    udecimals = 2,
                                    omit_unobserved = FALSE,
                                    xtreme = 0) {
  # Stage 1: Validate required inputs and extract core model objects.
  if (is.null(res) || is.null(diagnostics)) return(list())
  obs_df <- diagnostics$obs
  measures <- diagnostics$measures
  if (nrow(obs_df) == 0 || nrow(measures) == 0) return(list())

  prep <- res$prep
  config <- res$config
  rating_min <- prep$rating_min
  rating_max <- prep$rating_max
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  theta_mean <- if (length(theta_hat) > 0) mean(theta_hat, na.rm = TRUE) else 0
  facet_means <- purrr::map_dbl(config$facet_names, function(f) {
    mean(params$facets[[f]], na.rm = TRUE)
  })
  names(facet_means) <- config$facet_names
  facet_signs <- config$facet_signs
  if (is.null(facet_signs) || length(facet_signs) == 0) {
    facet_signs <- setNames(rep(-1, length(config$facet_names)), config$facet_names)
  }

  # Stage 2: Build step structures used for fair-average calculations.
  if (config$model == "RSM") {
    step_cum_common <- c(0, cumsum(params$steps))
    step_cum_mean <- step_cum_common
  } else {
    step_mat <- params$steps_mat
    if (is.null(step_mat) || length(step_mat) == 0) {
      step_cum_common <- numeric(0)
      step_cum_mean <- numeric(0)
    } else {
      step_mean <- colMeans(step_mat, na.rm = TRUE)
      step_cum_common <- t(apply(step_mat, 1, function(x) c(0, cumsum(x))))
      step_cum_mean <- c(0, cumsum(step_mean))
    }
  }

  facet_names <- c("Person", config$facet_names)
  facet_levels_all <- lapply(facet_names, function(facet) {
    if (facet == "Person") {
      prep$levels$Person
    } else {
      prep$levels[[facet]]
    }
  })
  names(facet_levels_all) <- facet_names

  # Stage 3: Detect extreme-only levels and cache row-level flags.
  extreme_levels <- purrr::map(facet_names, function(facet) {
    if (!facet %in% names(obs_df)) return(character(0))
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(
        (MinScore == prep$rating_min & MaxScore == prep$rating_min) |
          (MinScore == prep$rating_max & MaxScore == prep$rating_max)
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      pull(Level)
  })
  names(extreme_levels) <- facet_names

  extreme_flags <- list()
  for (facet in facet_names) {
    if (facet %in% names(obs_df)) {
      extreme_flags[[facet]] <- obs_df[[facet]] %in% extreme_levels[[facet]]
    }
  }
  if (length(extreme_flags) > 0) {
    extreme_flag_df <- as_tibble(extreme_flags)
    extreme_count <- rowSums(extreme_flag_df, na.rm = TRUE)
  } else {
    extreme_count <- rep(0, nrow(obs_df))
  }

  # Stage 4: Build one compatibility table per facet.
  out <- purrr::map(facet_names, function(facet) {
    if (!facet %in% names(obs_df)) return(tibble())
    status_tbl <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        TotalCountAll = n(),
        .groups = "drop"
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      select(Level, MinScore, MaxScore, TotalCountAll)

    if (isTRUE(totalscore)) {
      score_source <- obs_df
    } else {
      # Legacy printer-style rule: remove rows with multiple extreme flags,
      # but keep rows where the focal facet is the only extreme flag.
      flag <- extreme_flags[[facet]]
      if (is.null(flag)) {
        flag <- rep(FALSE, nrow(obs_df))
      }
      active_mask <- (extreme_count == 0) | ((extreme_count == 1) & flag)
      score_source <- obs_df[active_mask, , drop = FALSE]
    }
    score_tbl <- score_source |>
      group_by(.data[[facet]]) |>
      summarize(
        TotalScore = sum(Observed, na.rm = TRUE),
        TotalCount = n(),
        WeightdScore = sum(Observed * Weight, na.rm = TRUE),
        WeightdCount = sum(Weight, na.rm = TRUE),
        ObservedAverage = ifelse(sum(Weight, na.rm = TRUE) > 0,
                                 sum(Observed * Weight, na.rm = TRUE) / sum(Weight, na.rm = TRUE),
                                 NA_real_),
        .groups = "drop"
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      select(Level, TotalScore, TotalCount, WeightdScore, WeightdCount, ObservedAverage)

    level_tbl <- tibble(Level = as.character(facet_levels_all[[facet]]))
    score_tbl <- level_tbl |>
      left_join(score_tbl, by = "Level") |>
      mutate(
        TotalScore = replace_na(TotalScore, 0),
        TotalCount = replace_na(TotalCount, 0),
        WeightdScore = replace_na(WeightdScore, 0),
        WeightdCount = replace_na(WeightdCount, 0),
        ObservedAverage = ifelse(WeightdCount > 0, ObservedAverage, NA_real_)
      )

    meas_tbl <- measures |>
      filter(Facet == facet) |>
      mutate(Level = as.character(Level)) |>
      select(Level, Estimate, SE, ModelSE, RealSE, Infit, Outfit, InfitZSTD, OutfitZSTD, PTMEA)

    tbl <- level_tbl |>
      left_join(score_tbl, by = "Level") |>
      left_join(status_tbl, by = "Level") |>
      left_join(meas_tbl, by = "Level")

    anchor_status <- facet_anchor_status(facet, tbl$Level, config)
    status <- dplyr::case_when(
      tbl$TotalCountAll == 0 ~ "No data",
      tbl$MinScore == tbl$MaxScore & tbl$MinScore == rating_min ~ "Minimum",
      tbl$MinScore == tbl$MaxScore & tbl$MaxScore == rating_max ~ "Maximum",
      tbl$TotalCountAll == 1 ~ "One datum",
      TRUE ~ ""
    )

    sign_vec <- facet_signs[names(facet_means)]
    if (facet == "Person") {
      other_sum <- sum(sign_vec * facet_means, na.rm = TRUE)
      eta_m <- tbl$Estimate + other_sum
      eta_z <- tbl$Estimate
    } else {
      sign <- if (!is.null(facet_signs[[facet]])) facet_signs[[facet]] else -1
      other_sum <- sum(sign_vec[names(facet_means) != facet] *
                         facet_means[names(facet_means) != facet], na.rm = TRUE)
      eta_m <- theta_mean + other_sum + sign * tbl$Estimate
      eta_z <- sign * tbl$Estimate
    }

    if (config$model %in% c("PCM", "GPCM") && !is.null(config$step_facet)) {
      step_levels <- prep$levels[[config$step_facet]]
      if (facet == config$step_facet && length(step_levels) > 0 && length(step_cum_common) > 0) {
        step_cum_list <- purrr::map(tbl$Level, function(lvl) {
          idx <- match(lvl, step_levels)
          if (is.na(idx) || idx < 1 || idx > nrow(step_cum_common)) {
            step_cum_mean
          } else {
            step_cum_common[idx, ]
          }
        })
      } else {
        step_cum_list <- rep(list(step_cum_mean), nrow(tbl))
      }
    } else {
      step_cum_list <- rep(list(step_cum_common), nrow(tbl))
    }

    # Per-row slope vector for the slope-aware GPCM fair-average
    # construction. For the slope-facet's own rows, use that level's
    # own slope from `params$slopes`; for all other rows (Person,
    # Rater, ..., and any non-GPCM fit) the slope is 1. Setting the
    # non-slope-facet slope to 1 is the geometric-mean-one
    # identification convention: it represents the "average slope"
    # across slope-facet elements and makes
    # `expected_score_from_eta_gpcm()` reduce exactly to the PCM/RSM
    # Linacre fair-average. Net effect: only the slope-facet element
    # rows carry visible slope heterogeneity in the fair-average
    # column; non-slope rows remain continuous with the standard PCM
    # Linacre construction.
    slope_facet <- config$slope_facet %||% config$step_facet
    if (config$model == "GPCM" && !is.null(slope_facet) &&
        !is.null(params$slopes) && facet == slope_facet) {
      slope_levels <- prep$levels[[slope_facet]] %||% character(0)
      slopes_vec <- as.numeric(params$slopes)
      slope_list <- purrr::map_dbl(tbl$Level, function(lvl) {
        idx <- match(lvl, slope_levels)
        if (is.na(idx) || idx < 1L || idx > length(slopes_vec)) {
          1
        } else {
          slopes_vec[idx]
        }
      })
    } else {
      slope_list <- rep(1, nrow(tbl))
    }

    if (config$model == "GPCM") {
      fair_m <- purrr::pmap_dbl(
        list(eta_m, step_cum_list, slope_list),
        function(e, s, a) expected_score_from_eta_gpcm(e, s, a, rating_min)
      )
      fair_z <- purrr::pmap_dbl(
        list(eta_z, step_cum_list, slope_list),
        function(e, s, a) expected_score_from_eta_gpcm(e, s, a, rating_min)
      )
    } else {
      fair_m <- purrr::map2_dbl(eta_m, step_cum_list, ~ expected_score_from_eta(.x, .y, rating_min))
      fair_z <- purrr::map2_dbl(eta_z, step_cum_list, ~ expected_score_from_eta(.x, .y, rating_min))
    }

    xtreme_target <- ifelse(
      status == "Minimum", rating_min + xtreme,
      ifelse(status == "Maximum", rating_max - xtreme, NA_real_)
    )
    xtreme_eta <- purrr::pmap_dbl(
      list(xtreme_target, step_cum_list, slope_list),
      function(t, s, a) {
        if (!is.finite(t) || xtreme <= 0) return(NA_real_)
        estimate_eta_from_target(t, s, rating_min, rating_max, slope = a)
      }
    )

    # Convert any xtreme-adjusted eta back to facet-specific measure units.
    measure_logit <- tbl$Estimate
    if (any(is.finite(xtreme_eta))) {
      if (facet == "Person") {
        measure_logit <- ifelse(
          is.finite(xtreme_eta),
          xtreme_eta - other_sum,
          measure_logit
        )
      } else {
        sign <- if (!is.null(facet_signs[[facet]])) facet_signs[[facet]] else -1
        measure_logit <- ifelse(
          is.finite(xtreme_eta),
          (xtreme_eta - theta_mean - other_sum) / sign,
          measure_logit
        )
      }
    }

    scale_factor <- ifelse(is.finite(uscale), uscale, 1)
    scale_origin <- ifelse(is.finite(umean), umean, 0)

    tbl <- tbl |>
      mutate(
        Anchor = anchor_status,
        Status = status,
        FairM = fair_m,
        FairZ = fair_z,
        Measure = ifelse(is.finite(measure_logit), measure_logit * scale_factor + scale_origin, NA_real_),
        ModelSE = ifelse(
          is.finite(dplyr::coalesce(ModelSE, SE)),
          abs(scale_factor) * dplyr::coalesce(ModelSE, SE),
          NA_real_
        ),
        RealSE = ifelse(
          is.finite(dplyr::coalesce(RealSE, ModelSE, SE)),
          abs(scale_factor) * dplyr::coalesce(
            RealSE,
            ModelSE * sqrt(pmax(dplyr::coalesce(Infit, 1), 1)),
            SE * sqrt(pmax(dplyr::coalesce(Infit, 1), 1))
          ),
          NA_real_
        ),
        Infit = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, Infit),
        Outfit = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, Outfit),
        InfitZSTD = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, InfitZSTD),
        OutfitZSTD = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, OutfitZSTD)
      ) |>
      transmute(
        TotalScore,
        TotalCount,
        WeightdScore,
        WeightdCount,
        ObservedAverage,
        FairM,
        FairZ,
        Measure,
        ModelSE,
        RealSE,
        InfitMnSq = Infit,
        InfitZStd = InfitZSTD,
        OutfitMnSq = Outfit,
        OutfitZStd = OutfitZSTD,
        PtMeaCorr = ifelse(facet == "Person", NA_real_, PTMEA),
        Anchor,
        Status,
        Level,
        TotalCountAll
      ) |>
      arrange(desc(Measure), desc(TotalCount))

    if (isTRUE(omit_unobserved)) {
      tbl <- tbl |> filter(TotalCountAll > 0)
    }

    tbl |>
      select(-TotalCountAll)
  })

  names(out) <- facet_names
  out
}

format_facets_report_gt <- function(tbl,
                                    facet,
                                    decimals = 2,
                                    totalscore = TRUE,
                                    reference = c("both", "mean", "zero"),
                                    label_style = c("both", "native", "legacy")) {
  reference <- match.arg(reference)
  label_style <- match.arg(label_style)
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(data.frame(Message = "No facet report available.", stringsAsFactors = FALSE))
  }
  out <- as.data.frame(tbl, stringsAsFactors = FALSE)
  label_map <- c(
    TotalScore = if (isTRUE(totalscore)) "Total Score" else "Obsvd Score",
    TotalCount = if (isTRUE(totalscore)) "Total Count" else "Obsvd Count",
    WeightdScore = "Weightd Score",
    WeightdCount = "Weightd Count",
    ObservedAverage = "Obsvd Average",
    FairM = "Fair(M) Average",
    FairZ = "Fair(Z) Average",
    FairMSE = "Fair(M) S.E.",
    FairM_CI_Lower = "Fair(M) CI Lower",
    FairM_CI_Upper = "Fair(M) CI Upper",
    FairM_CI_Level = "Fair(M) CI Level",
    FairM_SE_Method = "Fair(M) S.E. Method",
    FairM_SE_Status = "Fair(M) S.E. Status",
    FairM_SE_Detail = "Fair(M) S.E. Detail",
    FairZSE = "Fair(Z) S.E.",
    FairZ_CI_Lower = "Fair(Z) CI Lower",
    FairZ_CI_Upper = "Fair(Z) CI Upper",
    FairZ_CI_Level = "Fair(Z) CI Level",
    FairZ_SE_Method = "Fair(Z) S.E. Method",
    FairZ_SE_Status = "Fair(Z) S.E. Status",
    FairZ_SE_Detail = "Fair(Z) S.E. Detail",
    Measure = "Measure",
    ModelSE = "Model S.E.",
    RealSE = "Real S.E.",
    InfitMnSq = "Infit MnSq",
    InfitZStd = "Infit ZStd",
    OutfitMnSq = "Outfit MnSq",
    OutfitZStd = "Outfit ZStd",
    PtMeaCorr = "PtMea Corr",
    Anchor = "Anch",
    Status = "Status",
    Level = "Element"
  )
  keep <- intersect(names(out), names(label_map))
  names(out)[match(keep, names(out))] <- unname(label_map[keep])
  if ("Element" %in% names(out)) {
    out <- out[, c(setdiff(names(out), "Element"), "Element"), drop = FALSE]
  }

  count_cols <- intersect(c("Total Score", "Total Count", "Weightd Score", "Weightd Count"), names(out))
  for (col in count_cols) out[[col]] <- round(as.numeric(out[[col]]), digits = 0)
  value_cols <- intersect(
    c("Obsvd Average", "Fair(M) Average", "Fair(Z) Average", "Measure", "Model S.E.", "Real S.E.",
      "Fair(M) S.E.", "Fair(M) CI Lower", "Fair(M) CI Upper", "Fair(M) CI Level",
      "Fair(Z) S.E.", "Fair(Z) CI Lower", "Fair(Z) CI Upper", "Fair(Z) CI Level",
      "Infit MnSq", "Infit ZStd", "Outfit MnSq", "Outfit ZStd", "PtMea Corr"),
    names(out)
  )
  for (col in value_cols) out[[col]] <- round(as.numeric(out[[col]]), digits = decimals)

  # Preserve the legacy display labels while exposing package-native aliases.
  if ("Obsvd Average" %in% names(out) && !"ObservedAverage" %in% names(out)) {
    out$ObservedAverage <- out[["Obsvd Average"]]
  }
  if ("Fair(M) Average" %in% names(out) && !"AdjustedAverage" %in% names(out)) {
    out$AdjustedAverage <- out[["Fair(M) Average"]]
  }
  if ("Fair(Z) Average" %in% names(out) && !"StandardizedAdjustedAverage" %in% names(out)) {
    out$StandardizedAdjustedAverage <- out[["Fair(Z) Average"]]
  }
  if ("Fair(M) S.E." %in% names(out) && !"AdjustedAverageSE" %in% names(out)) {
    out$AdjustedAverageSE <- out[["Fair(M) S.E."]]
  }
  if ("Fair(M) CI Lower" %in% names(out) && !"AdjustedAverageCI_Lower" %in% names(out)) {
    out$AdjustedAverageCI_Lower <- out[["Fair(M) CI Lower"]]
  }
  if ("Fair(M) CI Upper" %in% names(out) && !"AdjustedAverageCI_Upper" %in% names(out)) {
    out$AdjustedAverageCI_Upper <- out[["Fair(M) CI Upper"]]
  }
  if ("Fair(M) CI Level" %in% names(out) && !"AdjustedAverageCI_Level" %in% names(out)) {
    out$AdjustedAverageCI_Level <- out[["Fair(M) CI Level"]]
  }
  if ("Fair(M) S.E. Method" %in% names(out) && !"AdjustedAverageSEMethod" %in% names(out)) {
    out$AdjustedAverageSEMethod <- out[["Fair(M) S.E. Method"]]
  }
  if ("Fair(M) S.E. Status" %in% names(out) && !"AdjustedAverageSEStatus" %in% names(out)) {
    out$AdjustedAverageSEStatus <- out[["Fair(M) S.E. Status"]]
  }
  if ("Fair(M) S.E. Detail" %in% names(out) && !"AdjustedAverageSEDetail" %in% names(out)) {
    out$AdjustedAverageSEDetail <- out[["Fair(M) S.E. Detail"]]
  }
  if ("Fair(Z) S.E." %in% names(out) && !"StandardizedAdjustedAverageSE" %in% names(out)) {
    out$StandardizedAdjustedAverageSE <- out[["Fair(Z) S.E."]]
  }
  if ("Fair(Z) CI Lower" %in% names(out) && !"StandardizedAdjustedAverageCI_Lower" %in% names(out)) {
    out$StandardizedAdjustedAverageCI_Lower <- out[["Fair(Z) CI Lower"]]
  }
  if ("Fair(Z) CI Upper" %in% names(out) && !"StandardizedAdjustedAverageCI_Upper" %in% names(out)) {
    out$StandardizedAdjustedAverageCI_Upper <- out[["Fair(Z) CI Upper"]]
  }
  if ("Fair(Z) CI Level" %in% names(out) && !"StandardizedAdjustedAverageCI_Level" %in% names(out)) {
    out$StandardizedAdjustedAverageCI_Level <- out[["Fair(Z) CI Level"]]
  }
  if ("Fair(Z) S.E. Method" %in% names(out) && !"StandardizedAdjustedAverageSEMethod" %in% names(out)) {
    out$StandardizedAdjustedAverageSEMethod <- out[["Fair(Z) S.E. Method"]]
  }
  if ("Fair(Z) S.E. Status" %in% names(out) && !"StandardizedAdjustedAverageSEStatus" %in% names(out)) {
    out$StandardizedAdjustedAverageSEStatus <- out[["Fair(Z) S.E. Status"]]
  }
  if ("Fair(Z) S.E. Detail" %in% names(out) && !"StandardizedAdjustedAverageSEDetail" %in% names(out)) {
    out$StandardizedAdjustedAverageSEDetail <- out[["Fair(Z) S.E. Detail"]]
  }
  if ("Model S.E." %in% names(out) && !"ModelBasedSE" %in% names(out)) {
    out$ModelBasedSE <- out[["Model S.E."]]
  }
  if ("Real S.E." %in% names(out) && !"FitAdjustedSE" %in% names(out)) {
    out$FitAdjustedSE <- out[["Real S.E."]]
  }

  if (identical(reference, "mean")) {
    drop_cols <- intersect(c(
      "Fair(Z) Average", "StandardizedAdjustedAverage",
      "Fair(Z) S.E.", "Fair(Z) CI Lower", "Fair(Z) CI Upper", "Fair(Z) CI Level",
      "Fair(Z) S.E. Method", "Fair(Z) S.E. Status", "Fair(Z) S.E. Detail",
      "StandardizedAdjustedAverageSE", "StandardizedAdjustedAverageCI_Lower",
      "StandardizedAdjustedAverageCI_Upper", "StandardizedAdjustedAverageCI_Level",
      "StandardizedAdjustedAverageSEMethod", "StandardizedAdjustedAverageSEStatus",
      "StandardizedAdjustedAverageSEDetail"
    ), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  } else if (identical(reference, "zero")) {
    drop_cols <- intersect(c(
      "Fair(M) Average", "AdjustedAverage",
      "Fair(M) S.E.", "Fair(M) CI Lower", "Fair(M) CI Upper", "Fair(M) CI Level",
      "Fair(M) S.E. Method", "Fair(M) S.E. Status", "Fair(M) S.E. Detail",
      "AdjustedAverageSE", "AdjustedAverageCI_Lower", "AdjustedAverageCI_Upper",
      "AdjustedAverageCI_Level", "AdjustedAverageSEMethod",
      "AdjustedAverageSEStatus", "AdjustedAverageSEDetail"
    ), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  }

  if (identical(label_style, "native")) {
    drop_cols <- intersect(c(
      "Obsvd Average", "Fair(M) Average", "Fair(Z) Average",
      "Fair(M) S.E.", "Fair(M) CI Lower", "Fair(M) CI Upper", "Fair(M) CI Level",
      "Fair(M) S.E. Method", "Fair(M) S.E. Status", "Fair(M) S.E. Detail",
      "Fair(Z) S.E.", "Fair(Z) CI Lower", "Fair(Z) CI Upper", "Fair(Z) CI Level",
      "Fair(Z) S.E. Method", "Fair(Z) S.E. Status", "Fair(Z) S.E. Detail",
      "Model S.E.", "Real S.E."
    ), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  } else if (identical(label_style, "legacy")) {
    drop_cols <- intersect(c(
      "ObservedAverage", "AdjustedAverage", "StandardizedAdjustedAverage",
      "AdjustedAverageSE", "AdjustedAverageCI_Lower", "AdjustedAverageCI_Upper",
      "AdjustedAverageCI_Level", "AdjustedAverageSEMethod",
      "AdjustedAverageSEStatus", "AdjustedAverageSEDetail",
      "StandardizedAdjustedAverageSE", "StandardizedAdjustedAverageCI_Lower",
      "StandardizedAdjustedAverageCI_Upper", "StandardizedAdjustedAverageCI_Level",
      "StandardizedAdjustedAverageSEMethod", "StandardizedAdjustedAverageSEStatus",
      "StandardizedAdjustedAverageSEDetail",
      "ModelBasedSE", "FitAdjustedSE"
    ), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  }

  attr(out, "facet") <- facet
  attr(out, "totalscore") <- isTRUE(totalscore)
  attr(out, "reference") <- reference
  attr(out, "label_style") <- label_style
  out
}

compute_fair_average_value_from_par <- function(res,
                                                par,
                                                facet,
                                                level,
                                                metric = c("FairM", "FairZ")) {
  metric <- match.arg(metric)
  config <- res$config
  prep <- res$prep
  if (!identical(config$model, "GPCM")) return(NA_real_)
  facet <- as.character(facet[1])
  level <- as.character(level[1])
  if (identical(facet, "Person") || !facet %in% config$facet_names) {
    return(NA_real_)
  }

  sizes <- build_param_sizes(config)
  params <- expand_params(par, sizes, config)
  facet_levels <- as.character(prep$levels[[facet]] %||% character(0))
  facet_idx <- match(level, facet_levels)
  if (is.na(facet_idx) || facet_idx < 1L) return(NA_real_)

  theta_hat <- res$facets$person$Estimate
  theta_mean <- if (length(theta_hat) > 0L) mean(theta_hat, na.rm = TRUE) else 0
  if (!is.finite(theta_mean)) theta_mean <- 0

  facet_means <- purrr::map_dbl(config$facet_names, function(f) {
    vals <- params$facets[[f]]
    if (length(vals) == 0L) return(0)
    m <- mean(vals, na.rm = TRUE)
    if (is.finite(m)) m else 0
  })
  names(facet_means) <- config$facet_names

  facet_signs <- config$facet_signs
  if (is.null(facet_signs) || length(facet_signs) == 0L) {
    facet_signs <- stats::setNames(rep(-1, length(config$facet_names)), config$facet_names)
  }
  sign_vec <- suppressWarnings(as.numeric(facet_signs[names(facet_means)]))
  sign_vec[!is.finite(sign_vec)] <- -1
  names(sign_vec) <- names(facet_means)
  sign <- suppressWarnings(as.numeric(facet_signs[[facet]] %||% -1))
  if (!is.finite(sign)) sign <- -1

  estimate <- params$facets[[facet]][facet_idx]
  if (!is.finite(estimate)) return(NA_real_)
  other_sum <- sum(sign_vec[names(facet_means) != facet] *
                     facet_means[names(facet_means) != facet], na.rm = TRUE)
  eta <- if (identical(metric, "FairM")) {
    theta_mean + other_sum + sign * estimate
  } else {
    sign * estimate
  }
  if (!is.finite(eta)) return(NA_real_)

  step_mat <- params$steps_mat
  if (is.null(step_mat) || length(step_mat) == 0L) return(NA_real_)
  step_profile <- colMeans(step_mat, na.rm = TRUE)
  step_facet <- config$step_facet %||% NA_character_
  if (!is.na(step_facet) && identical(facet, step_facet)) {
    step_levels <- as.character(prep$levels[[step_facet]] %||% character(0))
    step_idx <- match(level, step_levels)
    if (!is.na(step_idx) && step_idx >= 1L && step_idx <= nrow(step_mat)) {
      step_profile <- step_mat[step_idx, ]
    }
  }
  if (!all(is.finite(step_profile))) return(NA_real_)
  step_cum <- c(0, cumsum(step_profile))

  slope <- 1
  slope_facet <- config$slope_facet %||% config$step_facet
  if (!is.null(slope_facet) && identical(facet, slope_facet) &&
      !is.null(params$slopes)) {
    slope_levels <- as.character(prep$levels[[slope_facet]] %||% character(0))
    slope_idx <- match(level, slope_levels)
    if (!is.na(slope_idx) && slope_idx >= 1L && slope_idx <= length(params$slopes)) {
      slope <- params$slopes[slope_idx]
    }
  }

  expected_score_from_eta_gpcm(
    eta = eta,
    step_cum = step_cum,
    slope = slope,
    rating_min = prep$rating_min
  )
}

finite_difference_gradient <- function(fn, par, rel_step = 1e-5) {
  par <- as.numeric(par %||% numeric(0))
  if (length(par) == 0L) return(numeric(0))
  grad <- rep(NA_real_, length(par))
  base <- tryCatch(fn(par), error = function(e) NA_real_)
  base <- suppressWarnings(as.numeric(base[1]))
  for (j in seq_along(par)) {
    h <- rel_step * max(1, abs(par[j]))
    par_p <- par
    par_m <- par
    par_p[j] <- par_p[j] + h
    par_m[j] <- par_m[j] - h
    f_p <- tryCatch(fn(par_p), error = function(e) NA_real_)
    f_m <- tryCatch(fn(par_m), error = function(e) NA_real_)
    f_p <- suppressWarnings(as.numeric(f_p[1]))
    f_m <- suppressWarnings(as.numeric(f_m[1]))
    if (is.finite(f_p) && is.finite(f_m)) {
      grad[j] <- (f_p - f_m) / (2 * h)
    } else if (is.finite(f_p) && is.finite(base)) {
      grad[j] <- (f_p - base) / h
    } else if (is.finite(f_m) && is.finite(base)) {
      grad[j] <- (base - f_m) / h
    }
  }
  grad
}

add_gpcm_fair_average_delta_se <- function(raw_tbls,
                                           res,
                                           covariance = NULL,
                                           ci_level = 0.95) {
  add_unavailable_columns <- function(tbl, status, detail) {
    if (is.null(tbl) || nrow(tbl) == 0L) return(tbl)
    for (prefix in c("FairM", "FairZ")) {
      tbl[[paste0(prefix, "SE")]] <- NA_real_
      tbl[[paste0(prefix, "_CI_Lower")]] <- NA_real_
      tbl[[paste0(prefix, "_CI_Upper")]] <- NA_real_
      tbl[[paste0(prefix, "_CI_Level")]] <- ci_level
      tbl[[paste0(prefix, "_SE_Method")]] <- "not available"
      tbl[[paste0(prefix, "_SE_Status")]] <- status
      tbl[[paste0(prefix, "_SE_Detail")]] <- detail
    }
    tbl
  }

  if (!identical(res$config$model, "GPCM")) {
    return(raw_tbls)
  }
  covariance <- covariance %||% compute_mml_parameter_covariance(res)
  if (!covariance$status %in% c("ok", "regularized") || is.null(covariance$cov)) {
    return(lapply(raw_tbls, add_unavailable_columns,
                  status = covariance$status,
                  detail = covariance$detail))
  }

  par <- as.numeric(res$opt$par %||% numeric(0))
  cov_mat <- symmetrize_matrix(covariance$cov)
  z <- stats::qnorm((1 + ci_level) / 2)
  rating_min <- suppressWarnings(as.numeric(res$prep$rating_min %||% NA_real_))
  rating_max <- suppressWarnings(as.numeric(res$prep$rating_max %||% NA_real_))
  se_method <- if (identical(covariance$status, "regularized")) {
    "Structural delta method (MML observed information; regularized Hessian)"
  } else {
    "Structural delta method (MML observed information)"
  }
  se_detail <- paste(
    "Propagates the structural covariance of facet, step, and slope",
    "parameters; MML person EAP estimates are conditioned on rather than",
    "included in the Hessian."
  )

  compute_one <- function(facet, level, estimate, metric) {
    if (!is.finite(estimate)) {
      return(list(se = NA_real_, lower = NA_real_, upper = NA_real_,
                  status = "not available",
                  detail = "Fair-average estimate is not finite."))
    }
    if (identical(facet, "Person")) {
      return(list(se = NA_real_, lower = NA_real_, upper = NA_real_,
                  status = "not available",
                  detail = paste(
                    "Person rows are not assigned structural fair-average SEs",
                    "because MML person EAPs are not part of the structural",
                    "observed-information Hessian."
                  )))
    }
    fn <- function(p) {
      compute_fair_average_value_from_par(
        res = res,
        par = p,
        facet = facet,
        level = level,
        metric = metric
      )
    }
    grad <- finite_difference_gradient(fn, par)
    if (length(grad) != ncol(cov_mat) || !any(is.finite(grad))) {
      return(list(se = NA_real_, lower = NA_real_, upper = NA_real_,
                  status = "not available",
                  detail = "Finite-difference gradient was not available."))
    }
    grad[!is.finite(grad)] <- 0
    var <- as.numeric(t(grad) %*% cov_mat %*% grad)
    if (is.finite(var) && var < 0 && abs(var) < 1e-10) var <- 0
    if (!is.finite(var) || var < 0) {
      return(list(se = NA_real_, lower = NA_real_, upper = NA_real_,
                  status = "not available",
                  detail = "Delta-method variance was not positive."))
    }
    se <- sqrt(var)
    lower <- estimate - z * se
    upper <- estimate + z * se
    if (is.finite(rating_min)) lower <- max(lower, rating_min)
    if (is.finite(rating_max)) upper <- min(upper, rating_max)
    list(se = se, lower = lower, upper = upper,
         status = covariance$status, detail = se_detail)
  }

  out <- lapply(names(raw_tbls), function(facet) {
    tbl <- raw_tbls[[facet]]
    if (is.null(tbl) || nrow(tbl) == 0L) return(tbl)
    facet_chr <- as.character(facet)
    for (prefix in c("FairM", "FairZ")) {
      estimate <- suppressWarnings(as.numeric(tbl[[prefix]]))
      vals <- lapply(seq_len(nrow(tbl)), function(i) {
        compute_one(
          facet = facet_chr,
          level = as.character(tbl$Level[i]),
          estimate = estimate[i],
          metric = prefix
        )
      })
      tbl[[paste0(prefix, "SE")]] <- vapply(vals, `[[`, numeric(1), "se")
      tbl[[paste0(prefix, "_CI_Lower")]] <- vapply(vals, `[[`, numeric(1), "lower")
      tbl[[paste0(prefix, "_CI_Upper")]] <- vapply(vals, `[[`, numeric(1), "upper")
      tbl[[paste0(prefix, "_CI_Level")]] <- ci_level
      tbl[[paste0(prefix, "_SE_Method")]] <- se_method
      tbl[[paste0(prefix, "_SE_Status")]] <- vapply(vals, `[[`, character(1), "status")
      tbl[[paste0(prefix, "_SE_Detail")]] <- vapply(vals, `[[`, character(1), "detail")
    }
    tbl
  })
  stats::setNames(out, names(raw_tbls))
}

calc_fair_average_bundle <- function(res,
                                     diagnostics,
                                     facets = NULL,
                                     totalscore = TRUE,
                                     umean = 0,
                                     uscale = 1,
                                     udecimals = 2,
                                     reference = c("both", "mean", "zero"),
                                     label_style = c("both", "native", "legacy"),
                                     omit_unobserved = FALSE,
                                     xtreme = 0,
                                     fair_se = FALSE,
                                     ci_level = 0.95) {
  reference <- match.arg(reference)
  label_style <- match.arg(label_style)
  raw_tbls <- calc_facets_report_tbls(
    res = res,
    diagnostics = diagnostics,
    totalscore = totalscore,
    umean = umean,
    uscale = uscale,
    udecimals = udecimals,
    omit_unobserved = omit_unobserved,
    xtreme = xtreme
  )
  if (is.null(raw_tbls) || length(raw_tbls) == 0) {
    return(list(raw_by_facet = list(), by_facet = list(), stacked = tibble()))
  }

  keep_names <- names(raw_tbls)
  if (!is.null(facets)) {
    keep_names <- intersect(keep_names, as.character(facets))
  }
  if (length(keep_names) == 0) {
    return(list(raw_by_facet = list(), by_facet = list(), stacked = tibble()))
  }
  raw_tbls <- raw_tbls[keep_names]
  if (isTRUE(fair_se)) {
    raw_tbls <- add_gpcm_fair_average_delta_se(
      raw_tbls = raw_tbls,
      res = res,
      covariance = NULL,
      ci_level = ci_level
    )
  }

  by_facet <- lapply(names(raw_tbls), function(facet) {
    tbl <- raw_tbls[[facet]]
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(data.frame())
    }
    format_facets_report_gt(
      tbl = tbl,
      facet = facet,
      decimals = udecimals,
      totalscore = totalscore,
      reference = reference,
      label_style = label_style
    )
  })
  names(by_facet) <- names(raw_tbls)
  by_facet <- by_facet[vapply(by_facet, nrow, integer(1)) > 0]

  stacked <- if (length(by_facet) == 0) {
    tibble()
  } else {
    dplyr::bind_rows(
      lapply(names(by_facet), function(facet) {
        df <- by_facet[[facet]]
        df$Facet <- facet
        df
      })
    ) |>
      select("Facet", everything())
  }

  list(
    raw_by_facet = raw_tbls,
    by_facet = by_facet,
    stacked = stacked
  )
}

calc_expected_category_counts <- function(res) {
  if (is.null(res)) return(tibble())
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  probs <- compute_response_probability_bundle(config, idx, params, eta)$probs
  if (length(probs) == 0) return(tibble())
  w <- idx$weight
  exp_counts <- if (is.null(w)) {
    colSums(probs, na.rm = TRUE)
  } else {
    colSums(probs * w, na.rm = TRUE)
  }
  total_exp <- sum(exp_counts)
  cat_vals <- seq(prep$rating_min, prep$rating_max)
  tibble(
    Category = cat_vals,
    ExpectedCount = exp_counts,
    ExpectedPercent = if (total_exp > 0) 100 * exp_counts / total_exp else rep(NA_real_, length(exp_counts))
  )
}

compute_mml_expected_category_diagnostics <- function(res,
                                                      include_p_geq = FALSE) {
  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  model <- as.character(res$config$model %||% res$summary$Model[1] %||% NA_character_)

  if (!identical(method, "MML")) {
    return(list(
      available = FALSE,
      reason = "Strict marginal diagnostics currently require an MML fit."
    ))
  }
  if (!model %in% c("RSM", "PCM", "GPCM")) {
    return(list(
      available = FALSE,
      reason = "Strict marginal diagnostics are currently implemented only for RSM, PCM, and bounded GPCM."
    ))
  }

  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  quad_points <- max(1L, as.integer(config$estimation_control$quad_points %||% 15L))
  quad <- gauss_hermite_normal(quad_points)
  base_eta <- compute_base_eta(idx, params, config)
  logprob_bundle <- mfrm_mml_logprob_bundle(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta,
    include_probs = TRUE
  )
  posterior_bundle <- mfrm_mml_posterior_bundle(logprob_bundle)
  expected_bundle <- mfrm_mml_expected_category_bundle(
    logprob_bundle = logprob_bundle,
    posterior_bundle = posterior_bundle,
    include_p_geq = include_p_geq
  )
  weights <- idx$weight
  if (is.null(weights) || length(weights) == 0L) {
    weights <- rep(1, nrow(prep$data))
  }
  weights <- ifelse(is.finite(weights) & weights > 0, weights, 0)

  list(
    available = TRUE,
    prep = prep,
    config = config,
    idx = idx,
    params = params,
    quad = quad,
    logprob_bundle = logprob_bundle,
    posterior_bundle = posterior_bundle,
    expected_bundle = expected_bundle,
    weights = weights,
    observed_cat = idx$score_k + prep$rating_min,
    categories = seq(prep$rating_min, prep$rating_max)
  )
}

summarize_marginal_fit_grid <- function(group_df,
                                        observed_cat,
                                        posterior_prob,
                                        weights,
                                        categories,
                                        abs_z_warn = 2,
                                        rmsd_warn = 0.05) {
  n <- length(observed_cat)
  if (is.null(group_df) || ncol(as.data.frame(group_df)) == 0L) {
    group_df <- data.frame(Scope = rep("All", n), stringsAsFactors = FALSE)
  } else {
    group_df <- as.data.frame(group_df, stringsAsFactors = FALSE)
  }
  if (nrow(group_df) != n) {
    stop("`group_df` must have one row per observation.", call. = FALSE)
  }

  group_cols <- names(group_df)
  group_key <- do.call(
    paste,
    c(lapply(group_df[group_cols], function(x) ifelse(is.na(x), "<NA>", as.character(x))), sep = "\r")
  )
  split_idx <- split(seq_len(n), group_key, drop = TRUE)
  exemplar_idx <- vapply(split_idx, `[`, integer(1), 1)
  group_index_tbl <- group_df[exemplar_idx, group_cols, drop = FALSE]

  cell_rows <- vector("list", length(split_idx))
  for (i in seq_along(split_idx)) {
    idx_g <- split_idx[[i]]
    probs_g <- posterior_prob[idx_g, , drop = FALSE]
    w_g <- weights[idx_g]
    total_w <- sum(w_g, na.rm = TRUE)

    cell_rows[[i]] <- bind_rows(lapply(seq_along(categories), function(j) {
      obs_count <- sum(w_g[observed_cat[idx_g] == categories[j]], na.rm = TRUE)
      exp_count <- sum(w_g * probs_g[, j], na.rm = TRUE)
      var_count <- sum((w_g^2) * probs_g[, j] * pmax(1 - probs_g[, j], 0), na.rm = TRUE)
      obs_prop <- if (is.finite(total_w) && total_w > 0) obs_count / total_w else NA_real_
      exp_prop <- if (is.finite(total_w) && total_w > 0) exp_count / total_w else NA_real_
      prop_diff <- obs_prop - exp_prop
      std_resid <- if (is.finite(var_count) && var_count > 0) {
        (obs_count - exp_count) / sqrt(var_count)
      } else {
        NA_real_
      }

      bind_cols(
        group_index_tbl[i, , drop = FALSE],
        tibble(
          Category = categories[j],
          GroupCount = total_w,
          ObservedCount = obs_count,
          ExpectedCount = exp_count,
          VarianceCount = var_count,
          ResidualCount = obs_count - exp_count,
          ObservedProp = obs_prop,
          ExpectedProp = exp_prop,
          PropDiff = prop_diff,
          SquaredPropDiff = prop_diff^2,
          StdResidual = std_resid,
          FlaggedAbsZ = is.finite(std_resid) && abs(std_resid) >= abs_z_warn
        )
      )
    }))
  }

  cell_stats <- bind_rows(cell_rows)
  summary_stats <- cell_stats |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      GroupCount = dplyr::first(.data$GroupCount),
      ObservedCountTotal = sum(.data$ObservedCount, na.rm = TRUE),
      ExpectedCountTotal = sum(.data$ExpectedCount, na.rm = TRUE),
      RMSD = sqrt(mean(.data$SquaredPropDiff, na.rm = TRUE)),
      MaxAbsStdResidual = if (all(!is.finite(.data$StdResidual))) {
        NA_real_
      } else {
        max(abs(.data$StdResidual), na.rm = TRUE)
      },
      FlaggedCellCount = sum(.data$FlaggedAbsZ, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      FlaggedRMSD = is.finite(.data$RMSD) & .data$RMSD >= rmsd_warn,
      Flagged = .data$FlaggedCellCount > 0 | .data$FlaggedRMSD
    )

  list(
    cell_stats = cell_stats,
    summary_stats = summary_stats
  )
}

build_marginal_fit_guidance <- function(pairwise_available = FALSE) {
  rows <- list(
    tibble(
      Component = "first_order_category_counts",
      StatisticLabel = "strict marginal category residual",
      LiteratureSeries = "limited_information_inspired_marginal_screen",
      PrimaryRole = "category / facet screening",
      InferenceTier = "exploratory",
      SupportsFormalInference = FALSE,
      FormalInferenceEligible = FALSE,
      PrimaryReportingEligible = FALSE,
      ClassificationSystem = "screening",
      ReportingUse = "screening_only",
      InterpretationNote = "Use as a latent-integrated limited-information screen; do not treat isolated flags as definitive inferential evidence."
    )
  )

  rows[[length(rows) + 1L]] <- tibble(
    Component = "pairwise_local_dependence",
    StatisticLabel = "strict pairwise agreement residual",
    LiteratureSeries = if (isTRUE(pairwise_available)) {
      "agreement_screen_informed_by_generalized_residual_logic"
    } else {
      "agreement_screen_informed_by_generalized_residual_logic_not_available"
    },
    PrimaryRole = "local dependence follow-up",
    InferenceTier = "exploratory",
    SupportsFormalInference = FALSE,
    FormalInferenceEligible = FALSE,
    PrimaryReportingEligible = FALSE,
    ClassificationSystem = "screening",
    ReportingUse = "screening_only",
    InterpretationNote = if (isTRUE(pairwise_available)) {
      "Use as an exploratory local-dependence follow-up after first-order marginal screening; this is informed by generalized-residual logic but is not a formal Haberman-Sinharay generalized residual test."
    } else {
      "Reserved for exploratory local-dependence follow-up when pairwise comparisons are available."
    }
  )

  rows[[length(rows) + 1L]] <- tibble(
    Component = "posterior_predictive_follow_up",
    StatisticLabel = "posterior predictive follow-up",
    LiteratureSeries = "posterior_predictive_model_check_follow_up",
    PrimaryRole = "misfit corroboration / practical-significance follow-up",
    InferenceTier = "exploratory",
    SupportsFormalInference = FALSE,
    FormalInferenceEligible = FALSE,
    PrimaryReportingEligible = FALSE,
    ClassificationSystem = "screening",
    ReportingUse = "screening_only",
    InterpretationNote = "Reserve for corroborating strict marginal flags and assessing practical significance; the current release labels this follow-up path but does not yet compute posterior predictive checks."
  )

  bind_rows(rows)
}

build_diagnostic_basis_guide <- function(diagnostic_mode = c("legacy", "marginal_fit", "both"),
                                         marginal_fit = NULL) {
  diagnostic_mode <- match.arg(diagnostic_mode, c("legacy", "marginal_fit", "both"))
  legacy_requested <- diagnostic_mode %in% c("legacy", "both")
  marginal_requested <- diagnostic_mode %in% c("marginal_fit", "both")
  marginal_available <- is.list(marginal_fit) && isTRUE(marginal_fit$available)
  pairwise_available <- marginal_available &&
    is.list(marginal_fit$pairwise) &&
    isTRUE(marginal_fit$pairwise$available)

  tibble(
    DiagnosticPath = c(
      "legacy_residual_fit",
      "strict_marginal_fit",
      "strict_pairwise_local_dependence",
      "posterior_predictive_follow_up"
    ),
    Component = c(
      "element_residual_fit",
      "first_order_category_counts",
      "pairwise_local_dependence",
      "posterior_predictive_follow_up"
    ),
    Status = c(
      if (legacy_requested) "computed" else "available_but_not_requested",
      if (marginal_requested && marginal_available) {
        "computed"
      } else if (marginal_requested) {
        "requested_not_available"
      } else {
        "not_requested"
      },
      if (marginal_requested && pairwise_available) {
        "computed"
      } else if (marginal_requested && marginal_available) {
        "not_available_for_run"
      } else if (marginal_requested) {
        "requested_not_available"
      } else {
        "not_requested"
      },
      "planned_not_implemented"
    ),
    Basis = c(
      "plugin_residuals_and_eap_tables",
      "latent_integrated_first_order_counts",
      "latent_integrated_second_order_agreement",
      "posterior_predictive_replication"
    ),
    PrimaryStatistics = c(
      "Infit / Outfit / ZSTD / PTMEA / residual QC bundles",
      "category residuals / standardized residuals / RMSD",
      "exact and adjacent agreement residuals",
      "replicated-data discrepancy checks"
    ),
    ReportingUse = c(
      "legacy_compatibility_screen",
      "screening_only",
      "screening_only",
      "screening_only"
    ),
    InterpretationNote = c(
      "Legacy fit statistics use plug-in residual machinery and should not be interpreted as latent-integrated marginal evidence.",
      "Strict marginal fit integrates over the latent distribution and is the preferred strict screen for category-level misfit in the current release.",
      "Strict pairwise local-dependence checks are exploratory follow-ups to first-order marginal flags, not standalone inferential tests.",
      "Posterior predictive checking is reserved for corroborating strict marginal flags and practical-significance review in a later release."
    )
  )
}

calc_marginal_pairwise_bundle <- function(expected_core,
                                          abs_z_warn = 2,
                                          exact_gap_warn = 0.10,
                                          adjacent_gap_warn = 0.10) {
  prep <- expected_core$prep
  config <- expected_core$config
  prob_list <- expected_core$logprob_bundle$prob_list
  obs_posterior <- expected_core$posterior_bundle$obs_posterior
  observed_cat <- expected_core$observed_cat
  weights <- expected_core$weights
  n <- length(observed_cat)
  k_cat <- length(expected_core$categories)

  if (length(config$facet_names) == 0L || length(prob_list) == 0L || n == 0L) {
    return(list(
      available = FALSE,
      reason = "No facet-wise pairwise contexts were available for strict second-order diagnostics."
    ))
  }

  adjacent_mask <- abs(row(matrix(0, k_cat, k_cat)) - col(matrix(0, k_cat, k_cat))) <= 1L
  pair_rows <- vector("list", 0L)
  pair_idx <- 0L

  for (facet in config$facet_names) {
    context_cols <- c("Person", setdiff(config$facet_names, facet))
    if (!all(context_cols %in% names(prep$data)) || !facet %in% names(prep$data)) {
      next
    }

    context_df <- prep$data[, context_cols, drop = FALSE]
    context_df[] <- lapply(context_df, function(x) ifelse(is.na(x), "<NA>", as.character(x)))
    group_key <- do.call(paste, c(unname(context_df), sep = "\r"))
    split_idx <- split(seq_len(n), group_key, drop = TRUE)

    for (idx_g in split_idx) {
      if (length(idx_g) < 2L) next
      level_vals <- as.character(prep$data[[facet]][idx_g])
      if (length(unique(level_vals)) < 2L) next

      pair_mat <- utils::combn(seq_along(idx_g), 2L)
      if (length(pair_mat) == 0L) next

      for (col_idx in seq_len(ncol(pair_mat))) {
        row_i <- idx_g[pair_mat[1L, col_idx]]
        row_j <- idx_g[pair_mat[2L, col_idx]]
        level_i <- as.character(prep$data[[facet]][row_i])
        level_j <- as.character(prep$data[[facet]][row_j])
        if (!nzchar(level_i) || !nzchar(level_j) || identical(level_i, level_j)) next

        level_pair <- sort(c(level_i, level_j), method = "radix")
        pair_weight <- weights[row_i] * weights[row_j]
        if (!is.finite(pair_weight) || pair_weight <= 0) next

        posterior <- obs_posterior[row_i, ]
        if (!all(is.finite(posterior))) next

        exact_q <- numeric(length(prob_list))
        adjacent_q <- numeric(length(prob_list))
        for (q in seq_along(prob_list)) {
          prob_i <- prob_list[[q]][row_i, ]
          prob_j <- prob_list[[q]][row_j, ]
          exact_q[q] <- sum(prob_i * prob_j)
          adjacent_q[q] <- sum(outer(prob_i, prob_j)[adjacent_mask])
        }

        exp_exact <- sum(posterior * exact_q)
        exp_adjacent <- sum(posterior * adjacent_q)
        obs_exact <- as.numeric(observed_cat[row_i] == observed_cat[row_j])
        obs_adjacent <- as.numeric(abs(observed_cat[row_i] - observed_cat[row_j]) <= 1)

        pair_idx <- pair_idx + 1L
        pair_rows[[pair_idx]] <- tibble(
          Facet = facet,
          Level1 = level_pair[1],
          Level2 = level_pair[2],
          PairWeight = pair_weight,
          ObservedExactCount = pair_weight * obs_exact,
          ExpectedExactCount = pair_weight * exp_exact,
          ExactVarianceCount = (pair_weight^2) * exp_exact * pmax(1 - exp_exact, 0),
          ObservedAdjacentCount = pair_weight * obs_adjacent,
          ExpectedAdjacentCount = pair_weight * exp_adjacent,
          AdjacentVarianceCount = (pair_weight^2) * exp_adjacent * pmax(1 - exp_adjacent, 0),
          SquaredExactDiff = (obs_exact - exp_exact)^2,
          SquaredAdjacentDiff = (obs_adjacent - exp_adjacent)^2
        )
      }
    }
  }

  if (length(pair_rows) == 0L) {
    return(list(
      available = FALSE,
      reason = "No within-context pairwise comparisons were available for strict second-order diagnostics."
    ))
  }

  pair_rows_df <- bind_rows(pair_rows)
  pair_stats <- pair_rows_df |>
    group_by(.data$Facet, .data$Level1, .data$Level2) |>
    summarize(
      LevelPairCount = dplyr::n(),
      OpportunityWeight = sum(.data$PairWeight, na.rm = TRUE),
      ObservedExactCount = sum(.data$ObservedExactCount, na.rm = TRUE),
      ExpectedExactCount = sum(.data$ExpectedExactCount, na.rm = TRUE),
      ExactVarianceCount = sum(.data$ExactVarianceCount, na.rm = TRUE),
      ObservedAdjacentCount = sum(.data$ObservedAdjacentCount, na.rm = TRUE),
      ExpectedAdjacentCount = sum(.data$ExpectedAdjacentCount, na.rm = TRUE),
      AdjacentVarianceCount = sum(.data$AdjacentVarianceCount, na.rm = TRUE),
      ExactRMSD = sqrt(stats::weighted.mean(.data$SquaredExactDiff, .data$PairWeight, na.rm = TRUE)),
      AdjacentRMSD = sqrt(stats::weighted.mean(.data$SquaredAdjacentDiff, .data$PairWeight, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(
      ExactAgreement = ifelse(.data$OpportunityWeight > 0, .data$ObservedExactCount / .data$OpportunityWeight, NA_real_),
      ExpectedExactAgreement = ifelse(.data$OpportunityWeight > 0, .data$ExpectedExactCount / .data$OpportunityWeight, NA_real_),
      ExactGap = .data$ExactAgreement - .data$ExpectedExactAgreement,
      ExactStdResidual = ifelse(
        is.finite(.data$ExactVarianceCount) & .data$ExactVarianceCount > 0,
        (.data$ObservedExactCount - .data$ExpectedExactCount) / sqrt(.data$ExactVarianceCount),
        NA_real_
      ),
      AdjacentAgreement = ifelse(.data$OpportunityWeight > 0, .data$ObservedAdjacentCount / .data$OpportunityWeight, NA_real_),
      ExpectedAdjacentAgreement = ifelse(.data$OpportunityWeight > 0, .data$ExpectedAdjacentCount / .data$OpportunityWeight, NA_real_),
      AdjacentGap = .data$AdjacentAgreement - .data$ExpectedAdjacentAgreement,
      AdjacentStdResidual = ifelse(
        is.finite(.data$AdjacentVarianceCount) & .data$AdjacentVarianceCount > 0,
        (.data$ObservedAdjacentCount - .data$ExpectedAdjacentCount) / sqrt(.data$AdjacentVarianceCount),
        NA_real_
      ),
      FlaggedExact = (is.finite(.data$ExactStdResidual) & abs(.data$ExactStdResidual) >= abs_z_warn) |
        (is.finite(.data$ExactGap) & abs(.data$ExactGap) >= exact_gap_warn),
      FlaggedAdjacent = (is.finite(.data$AdjacentStdResidual) & abs(.data$AdjacentStdResidual) >= abs_z_warn) |
        (is.finite(.data$AdjacentGap) & abs(.data$AdjacentGap) >= adjacent_gap_warn),
      Flagged = .data$FlaggedExact | .data$FlaggedAdjacent
    ) |>
    mutate(
      InferenceTier = "exploratory",
      SupportsFormalInference = FALSE,
      FormalInferenceEligible = FALSE,
      PrimaryReportingEligible = FALSE,
      ClassificationSystem = "screening",
      ReportingUse = "screening_only",
      StatisticLabel = "strict pairwise agreement residual",
      LiteratureSeries = "agreement_screen_informed_by_generalized_residual_logic"
    )

  facet_summary <- pair_stats |>
    group_by(.data$Facet) |>
    summarize(
      LevelPairs = dplyr::n(),
      ContextPairs = sum(.data$LevelPairCount, na.rm = TRUE),
      MeanAbsExactGap = stats::weighted.mean(abs(.data$ExactGap), .data$OpportunityWeight, na.rm = TRUE),
      MeanAbsAdjacentGap = stats::weighted.mean(abs(.data$AdjacentGap), .data$OpportunityWeight, na.rm = TRUE),
      MaxAbsExactStdResidual = if (all(!is.finite(.data$ExactStdResidual))) {
        NA_real_
      } else {
        max(abs(.data$ExactStdResidual), na.rm = TRUE)
      },
      MaxAbsAdjacentStdResidual = if (all(!is.finite(.data$AdjacentStdResidual))) {
        NA_real_
      } else {
        max(abs(.data$AdjacentStdResidual), na.rm = TRUE)
      },
      OpportunityWeight = sum(.data$OpportunityWeight, na.rm = TRUE),
      FlaggedLevelPairs = sum(.data$Flagged, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      InferenceTier = "exploratory",
      SupportsFormalInference = FALSE,
      FormalInferenceEligible = FALSE,
      PrimaryReportingEligible = FALSE,
      ClassificationSystem = "screening",
      ReportingUse = "screening_only",
      StatisticLabel = "strict pairwise agreement residual",
      LiteratureSeries = "agreement_screen_informed_by_generalized_residual_logic"
    )

  top_pairs <- pair_stats |>
    mutate(
      AbsExactStdResidual = abs(.data$ExactStdResidual),
      AbsAdjacentStdResidual = abs(.data$AdjacentStdResidual)
    ) |>
    arrange(
      desc(.data$AbsExactStdResidual),
      desc(.data$AbsAdjacentStdResidual),
      desc(abs(.data$ExactGap)),
      desc(abs(.data$AdjacentGap))
    ) |>
    slice_head(n = 20)

  nonunit_weights <- any(abs(weights - 1) > sqrt(.Machine$double.eps), na.rm = TRUE)

  list(
    available = TRUE,
    pair_stats = pair_stats,
    facet_summary = facet_summary,
    top_pairs = top_pairs,
    thresholds = list(
      abs_z_warn = abs_z_warn,
      exact_gap_warn = exact_gap_warn,
      adjacent_gap_warn = adjacent_gap_warn
    ),
    notes = c(
      "Strict second-order diagnostics compare observed and posterior-integrated expected exact/adjacent agreement within Person x remaining-facets contexts.",
      "These agreement screens are informed by generalized-residual logic but are not formal Haberman-Sinharay generalized residual tests.",
      if (isTRUE(nonunit_weights)) {
        "Pairwise opportunity counts use the product of row weights."
      } else {
        "Pairwise opportunity counts are unweighted row-pair counts."
      }
    )
  )
}

calc_marginal_fit_bundle <- function(res,
                                     abs_z_warn = 2,
                                     rmsd_warn = 0.05) {
  expected_core <- compute_mml_expected_category_diagnostics(res)
  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  model <- as.character(res$config$model %||% res$summary$Model[1] %||% NA_character_)

  if (!isTRUE(expected_core$available)) {
    return(list(
      available = FALSE,
      summary = tibble(
        Available = FALSE,
        Method = method,
        Model = model,
        Basis = "latent_integrated_first_order_counts",
        Reason = as.character(expected_core$reason %||% "Marginal fit diagnostics were not available.")
      ),
      thresholds = list(abs_z_warn = abs_z_warn, rmsd_warn = rmsd_warn),
      notes = as.character(expected_core$reason %||% "Marginal fit diagnostics were not available.")
    ))
  }

  n <- length(expected_core$observed_cat)
  weights <- expected_core$weights
  posterior_prob <- expected_core$expected_bundle$posterior_prob
  categories <- expected_core$categories
  prep <- expected_core$prep
  config <- expected_core$config

  overall <- summarize_marginal_fit_grid(
    group_df = data.frame(CellType = rep("overall", n), Scope = rep("All", n), stringsAsFactors = FALSE),
    observed_cat = expected_core$observed_cat,
    posterior_prob = posterior_prob,
    weights = weights,
    categories = categories,
    abs_z_warn = abs_z_warn,
    rmsd_warn = rmsd_warn
  )
  step_or_scale <- summarize_marginal_fit_grid(
    group_df = data.frame(
      CellType = rep("step_facet", n),
      StepFacet = if (identical(config$model, "PCM")) {
        as.character(prep$data[[config$step_facet]])
      } else {
        rep("Common", n)
      },
      stringsAsFactors = FALSE
    ),
    observed_cat = expected_core$observed_cat,
    posterior_prob = posterior_prob,
    weights = weights,
    categories = categories,
    abs_z_warn = abs_z_warn,
    rmsd_warn = rmsd_warn
  )
  facet_parts <- lapply(config$facet_names, function(facet) {
    summarize_marginal_fit_grid(
      group_df = data.frame(
        CellType = rep("facet_level", n),
        Facet = rep(facet, n),
        Level = as.character(prep$data[[facet]]),
        stringsAsFactors = FALSE
      ),
      observed_cat = expected_core$observed_cat,
      posterior_prob = posterior_prob,
      weights = weights,
      categories = categories,
      abs_z_warn = abs_z_warn,
      rmsd_warn = rmsd_warn
    )
  })
  facet_level <- list(
    cell_stats = bind_rows(lapply(facet_parts, `[[`, "cell_stats")),
    summary_stats = bind_rows(lapply(facet_parts, `[[`, "summary_stats"))
  )
  pairwise <- calc_marginal_pairwise_bundle(
    expected_core = expected_core,
    abs_z_warn = abs_z_warn
  )
  guidance_tbl <- build_marginal_fit_guidance(pairwise_available = isTRUE(pairwise$available))
  top_cells <- bind_rows(
    step_or_scale$cell_stats,
    facet_level$cell_stats
  ) |>
    mutate(AbsStdResidual = abs(.data$StdResidual)) |>
    arrange(desc(.data$AbsStdResidual), desc(abs(.data$PropDiff))) |>
    slice_head(n = 20)

  overall_summary <- overall$summary_stats |>
    slice_head(n = 1)
  summary_tbl <- tibble(
    Available = TRUE,
    Method = method,
    Model = model,
    Basis = "latent_integrated_first_order_counts",
    ObservationCount = n,
    CategoryCount = length(categories),
    OverallRMSD = overall_summary$RMSD[1] %||% NA_real_,
    OverallMaxAbsStdResidual = overall_summary$MaxAbsStdResidual[1] %||% NA_real_,
    StepGroupsFlagged = sum(step_or_scale$summary_stats$Flagged, na.rm = TRUE),
    FacetLevelsFlagged = sum(facet_level$summary_stats$Flagged, na.rm = TRUE),
    PairwiseAvailable = isTRUE(pairwise$available),
    PairwiseFlaggedLevelPairs = if (isTRUE(pairwise$available)) {
      sum(pairwise$pair_stats$Flagged, na.rm = TRUE)
    } else {
      NA_integer_
    },
    PosteriorPredictiveFollowUp = "planned_not_implemented",
    InferenceTier = "exploratory",
    SupportsFormalInference = FALSE,
    FormalInferenceEligible = FALSE,
    PrimaryReportingEligible = FALSE,
    ClassificationSystem = "screening",
    ReportingUse = "screening_only",
    StatisticLabel = "strict marginal fit screen",
    LiteratureSeries = if (isTRUE(pairwise$available)) {
      "limited_information_inspired_marginal_screen + agreement_screen_informed_by_generalized_residual_logic"
    } else {
      "limited_information_inspired_marginal_screen"
    },
    InterpretationBasis = "latent_integrated",
    FollowUpRecommendation = "Use as screening evidence and corroborate with companion diagnostics, model comparisons, and substantive design review."
  )

  list(
    available = TRUE,
    summary = summary_tbl,
    overall = overall,
    step_or_scale = step_or_scale,
    facet_level = facet_level,
    pairwise = pairwise,
    guidance = guidance_tbl,
    top_cells = top_cells,
    thresholds = list(abs_z_warn = abs_z_warn, rmsd_warn = rmsd_warn),
    notes = c(
      "Strict marginal diagnostics use latent-integrated posterior-expected first-order category counts.",
      "The current first release summarizes overall, step-facet/scale, facet-level category residuals, and pairwise local-dependence checks.",
      "Posterior predictive checking is currently a planned corroborating follow-up path rather than a default computed statistic."
    )
  )
}

calc_category_stats <- function(obs_df, res = NULL, whexact = FALSE) {
  if (nrow(obs_df) == 0) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  total_n <- sum(obs_df$.Weight, na.rm = TRUE)
  obs_summary <- obs_df |>
    group_by(Observed) |>
    summarize(
      Count = sum(.Weight, na.rm = TRUE),
      AvgPersonMeasure = weighted_mean(PersonMeasure, .Weight),
      ExpectedAverage = weighted_mean(Expected, .Weight),
      Infit = sum(StdSq * Var * .Weight, na.rm = TRUE) / sum(Var * .Weight, na.rm = TRUE),
      Outfit = sum(StdSq * .Weight, na.rm = TRUE) / sum(.Weight, na.rm = TRUE),
      MeanResidual = weighted_mean(Residual, .Weight),
      DF_Infit = sum(Var * .Weight, na.rm = TRUE),
      DF_Outfit = sum(.Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(Category = Observed)

  all_categories <- if (!is.null(res)) {
    seq(res$prep$rating_min, res$prep$rating_max)
  } else {
    sort(unique(obs_df$Observed))
  }

  cat_tbl <- tibble(Category = all_categories) |>
    left_join(obs_summary, by = "Category") |>
    mutate(
      Count = replace_na(Count, 0),
      Percent = if (total_n > 0) 100 * Count / total_n else NA_real_,
      InfitZSTD = zstd_from_mnsq(Infit, DF_Infit, whexact = whexact),
      OutfitZSTD = zstd_from_mnsq(Outfit, DF_Outfit, whexact = whexact)
    )

  exp_tbl <- calc_expected_category_counts(res)
  if (nrow(exp_tbl) > 0) {
    cat_tbl <- cat_tbl |>
      left_join(exp_tbl, by = "Category") |>
      mutate(
        DiffCount = ifelse(is.finite(ExpectedCount), Count - ExpectedCount, NA_real_),
        DiffPercent = ifelse(is.finite(ExpectedPercent), Percent - ExpectedPercent, NA_real_)
      )
  }

  cat_tbl |>
    mutate(
      LowCount = Count < 10,
      InfitFlag = !is.na(Infit) & (Infit < 0.5 | Infit > 1.5),
      OutfitFlag = !is.na(Outfit) & (Outfit < 0.5 | Outfit > 1.5),
      ZSTDFlag = (is.finite(InfitZSTD) & abs(InfitZSTD) >= 2) |
        (is.finite(OutfitZSTD) & abs(OutfitZSTD) >= 2)
    ) |>
    arrange(Category)
}

make_union_find <- function(nodes) {
  parent <- setNames(nodes, nodes)
  # Iterative `find_root` with path compression. The earlier recursive
  # implementation hit `options(expressions)` (default 5,000) on
  # degenerate union chains in large designs (Round 4 K-audit).
  find_root <- function(x) {
    if (is.null(parent[[x]])) return(NA_character_)
    # Walk up to the root.
    root <- x
    repeat {
      pr <- parent[[root]]
      if (is.null(pr) || pr == root) break
      root <- pr
    }
    # Path compression: every node on the walked path now points at root.
    cur <- x
    while (!identical(parent[[cur]], root)) {
      nxt <- parent[[cur]]
      parent[[cur]] <<- root
      cur <- nxt
      if (is.null(cur) || cur == root) break
    }
    root
  }
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (is.na(ra) || is.na(rb) || ra == rb) return(NULL)
    parent[[rb]] <<- ra
  }
  list(find = find_root, union = union_nodes)
}

calc_subsets <- function(obs_df, facet_cols) {
  if (is.null(obs_df) || nrow(obs_df) == 0 || length(facet_cols) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }
  df <- obs_df |>
    select(all_of(facet_cols)) |>
    filter(if_all(everything(), ~ !is.na(.)))
  if (nrow(df) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }

  nodes <- unique(unlist(lapply(facet_cols, function(facet) {
    paste0(facet, ":", as.character(df[[facet]]))
  })))
  if (length(nodes) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }
  uf <- make_union_find(nodes)

  for (i in seq_len(nrow(df))) {
    row_vals <- unlist(df[i, facet_cols, drop = FALSE], use.names = FALSE)
    row_nodes <- paste0(facet_cols, ":", as.character(row_vals))
    row_nodes <- row_nodes[!is.na(row_nodes)]
    if (length(row_nodes) < 2) next
    base <- row_nodes[1]
    for (node in row_nodes[-1]) {
      uf$union(base, node)
    }
  }

  comp_ids <- vapply(nodes, uf$find, character(1))
  comp_levels <- unique(comp_ids)
  comp_index <- setNames(seq_along(comp_levels), comp_levels)

  node_tbl <- tibble(
    Node = nodes,
    Component = comp_ids,
    Subset = comp_index[comp_ids],
    Facet = stringr::str_extract(nodes, "^[^:]+"),
    Level = stringr::str_replace(nodes, "^[^:]+:", "")
  )

  facet_counts <- node_tbl |>
    group_by(Subset, Facet) |>
    summarize(Levels = n_distinct(Level), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Facet, values_from = Levels, values_fill = 0)

  row_subset <- vapply(seq_len(nrow(df)), function(i) {
    node <- paste0(facet_cols[1], ":", as.character(df[[facet_cols[1]]][i]))
    comp_index[uf$find(node)]
  }, integer(1))
  obs_counts <- tibble(Subset = row_subset) |>
    count(Subset, name = "Observations")

  summary_tbl <- facet_counts |>
    left_join(obs_counts, by = "Subset") |>
    mutate(Observations = replace_na(Observations, 0)) |>
    arrange(desc(Observations))

  list(summary = summary_tbl, nodes = node_tbl)
}

calc_step_order <- function(step_tbl) {
  if (is.null(step_tbl) || nrow(step_tbl) == 0) return(tibble())
  step_tbl <- step_tbl |>
    mutate(StepIndex = suppressWarnings(as.integer(stringr::str_extract(Step, "\\d+"))))
  if (!"StepFacet" %in% names(step_tbl)) {
    step_tbl <- mutate(step_tbl, StepFacet = "Common")
  }
  step_tbl |>
    arrange(StepFacet, StepIndex) |>
    group_by(StepFacet) |>
    mutate(
      Spacing = Estimate - lag(Estimate),
      Ordered = ifelse(is.na(Spacing), NA, Spacing > 0)
    ) |>
    ungroup()
}

category_warnings_text <- function(cat_tbl, step_tbl = NULL) {
  if (is.null(cat_tbl) || nrow(cat_tbl) == 0) return("No category diagnostics available.")
  msgs <- character()
  unused <- cat_tbl |>
    filter(Count == 0)
  if (nrow(unused) > 0) {
    msgs <- c(msgs, paste0("Unused categories: ", paste(unused$Category, collapse = ", ")))
  }
  low_counts <- cat_tbl |>
    filter(Count < 10)
  if (nrow(low_counts) > 0) {
    msgs <- c(msgs, paste0("Low category counts (<10): ", paste(low_counts$Category, collapse = ", ")))
  }
  if ("DiffPercent" %in% names(cat_tbl)) {
    diff_bad <- cat_tbl |>
      filter(is.finite(DiffPercent), abs(DiffPercent) >= 5)
    if (nrow(diff_bad) > 0) {
      msgs <- c(msgs, paste0("Observed vs expected % differs by >= 5: ", paste(diff_bad$Category, collapse = ", ")))
    }
  }
  if ("InfitZSTD" %in% names(cat_tbl) && "OutfitZSTD" %in% names(cat_tbl)) {
    zstd_bad <- cat_tbl |>
      filter((is.finite(InfitZSTD) & abs(InfitZSTD) >= 2) | (is.finite(OutfitZSTD) & abs(OutfitZSTD) >= 2))
    if (nrow(zstd_bad) > 0) {
      msgs <- c(msgs, paste0("Category |ZSTD| >= 2: ", paste(zstd_bad$Category, collapse = ", ")))
    }
  }
  avg_tbl <- cat_tbl |>
    filter(is.finite(AvgPersonMeasure)) |>
    arrange(Category)
  if (nrow(avg_tbl) >= 3 && is.unsorted(avg_tbl$AvgPersonMeasure, strictly = FALSE)) {
    msgs <- c(msgs, "Category averages are not monotonic (Avg Measure by category).")
  }
  if (!is.null(step_tbl) && nrow(step_tbl) > 0) {
    disordered <- step_tbl |>
      filter(!is.na(Ordered), Ordered == FALSE)
    if (nrow(disordered) > 0) {
      bad <- disordered |>
        mutate(Label = paste0(StepFacet, ":", Step)) |>
        pull(Label)
      msgs <- c(msgs, paste0("Disordered thresholds detected: ", paste(bad, collapse = ", ")))
    }
  }
  if (length(msgs) == 0) {
    "No major category warnings detected."
  } else {
    paste(msgs, collapse = "\n")
  }
}

get_extreme_levels <- function(obs_df, facet_names, rating_min, rating_max) {
  out <- list()
  for (facet in facet_names) {
    if (!facet %in% names(obs_df)) {
      out[[facet]] <- character(0)
      next
    }
    stat <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        .groups = "drop"
      )
    extreme <- stat |>
      filter(
        (MinScore == rating_min & MaxScore == rating_min) |
          (MinScore == rating_max & MaxScore == rating_max)
      )
    out[[facet]] <- as.character(extreme[[facet]])
  }
  out
}

estimate_bias_interaction <- function(res,
                                      diagnostics,
                                      facet_a = NULL,
                                      facet_b = NULL,
                                      interaction_facets = NULL,
                                      max_abs = 10,
                                      omit_extreme = TRUE,
                                      max_iter = 4,
                                      tol = 1e-3) {
  if (is.null(res) || is.null(diagnostics)) return(list())
  obs_df <- diagnostics$obs
  if (is.null(obs_df) || nrow(obs_df) == 0) return(list())

  facet_names <- c("Person", res$config$facet_names)
  selected_facets <- if (!is.null(interaction_facets)) {
    as.character(interaction_facets)
  } else if (!is.null(facet_a) && !is.null(facet_b)) {
    c(as.character(facet_a[1]), as.character(facet_b[1]))
  } else {
    character(0)
  }
  selected_facets <- selected_facets[!is.na(selected_facets) & nzchar(selected_facets)]
  selected_facets <- unique(selected_facets)
  if (length(selected_facets) < 2) return(list())
  if (!all(selected_facets %in% facet_names)) return(list())
  if (!all(selected_facets %in% names(obs_df))) return(list())

  prep <- res$prep
  config <- res$config
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  idx <- build_indices(prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  theta_hat <- if (config$method == "JMLE") params$theta else res$facets$person$Estimate
  eta_base <- compute_eta(idx, params, config, theta_override = theta_hat)
  score_k <- idx$score_k
  weight <- idx$weight
  step_idx <- idx$step_idx
  # Per-observation slope index for GPCM. For non-GPCM fits these stay
  # NULL and the dispatch below falls through to the RSM/PCM branch.
  slope_idx <- if (identical(config$model, "GPCM")) {
    idx$slope_idx %||% idx$step_idx
  } else {
    NULL
  }
  slopes_full <- if (identical(config$model, "GPCM")) {
    as.numeric(params$slopes)
  } else {
    NULL
  }

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    step_cum_mat <- NULL
  } else {
    step_cum <- NULL
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
  }

  if (omit_extreme) {
    extreme_levels <- get_extreme_levels(
      obs_df,
      selected_facets,
      prep$rating_min,
      prep$rating_max
    )
  } else {
    extreme_levels <- list()
  }

  measures <- diagnostics$measures
  meas_map <- list()
  se_map <- list()
  if (!is.null(measures) && nrow(measures) > 0) {
    for (i in seq_len(nrow(measures))) {
      key <- paste(measures$Facet[i], as.character(measures$Level[i]), sep = "||")
      meas_map[[key]] <- measures$Estimate[i]
      se_map[[key]] <- measures$SE[i]
    }
  }

  level_map <- lapply(facet_names, function(facet) {
    if (facet == "Person") {
      as.character(prep$levels$Person)
    } else {
      as.character(prep$levels[[facet]])
    }
  })
  names(level_map) <- facet_names

  interaction_parts <- lapply(selected_facets, function(f) as.character(obs_df[[f]]))
  names(interaction_parts) <- selected_facets
  group_keys <- do.call(interaction, c(interaction_parts, list(drop = TRUE, sep = "||", lex.order = TRUE)))
  group_indices <- split(seq_len(nrow(obs_df)), group_keys)

  groups <- list()
  for (key in names(group_indices)) {
    lvl_vals <- strsplit(as.character(key), "\\|\\|")[[1]]
    if (length(lvl_vals) < length(selected_facets)) next
    lvl_vals <- as.character(lvl_vals[seq_along(selected_facets)])
    names(lvl_vals) <- selected_facets
    if (isTRUE(omit_extreme)) {
      is_extreme <- vapply(selected_facets, function(f) {
        lvl_vals[[f]] %in% (extreme_levels[[f]] %||% character(0))
      }, logical(1))
      if (any(is_extreme)) next
    }
    idx_rows <- group_indices[[key]]
    if (length(idx_rows) == 0) next
    groups[[as.character(key)]] <- list(levels = lvl_vals, idx = idx_rows)
  }
  if (length(groups) == 0) return(list())

  orientation_review <- audit_interaction_orientation(config, selected_facets)
  has_bias_error <- function(x) {
    is.character(x) && length(x) > 0L && !is.na(x[1]) && isTRUE(nzchar(x[1]))
  }

  make_bias_nll <- function(eta_sub,
                            score_k_sub,
                            weight_sub = NULL,
                            step_idx_sub = NULL,
                            slope_idx_sub = NULL) {
    if (config$model == "RSM") {
      return(function(b) -loglik_rsm(eta_sub + b, score_k_sub, step_cum,
                                     weight = weight_sub))
    }
    if (identical(config$model, "GPCM")) {
      return(function(b) -loglik_gpcm(eta_sub + b, score_k_sub, step_cum_mat,
                                      criterion_idx = step_idx_sub,
                                      slopes = slopes_full,
                                      slope_idx = slope_idx_sub,
                                      weight = weight_sub))
    }
    function(b) -loglik_pcm(eta_sub + b, score_k_sub, step_cum_mat,
                            step_idx_sub, weight = weight_sub)
  }

  safe_nll_value <- function(nll, b) {
    val <- tryCatch(nll(b), error = function(e) NA_real_)
    val <- suppressWarnings(as.numeric(val[1]))
    if (is.finite(val)) val else NA_real_
  }

  profile_bias_ci <- function(nll,
                              estimate,
                              nll_min,
                              max_abs,
                              level = 0.95) {
    empty <- list(
      lower = NA_real_,
      upper = NA_real_,
      level = level,
      status = "not available"
    )
    if (!is.finite(estimate) || !is.finite(nll_min) ||
        !is.finite(max_abs) || max_abs <= 0 ||
        !is.finite(level) || level <= 0 || level >= 1) {
      return(empty)
    }
    cutoff <- stats::qchisq(level, df = 1)
    target <- function(b) {
      val <- safe_nll_value(nll, b)
      if (!is.finite(val)) return(NA_real_)
      2 * (val - nll_min) - cutoff
    }
    root_side <- function(bound) {
      if (abs(bound - estimate) <= sqrt(.Machine$double.eps)) {
        return(list(value = bound, status = "limited by search range"))
      }
      target_bound <- target(bound)
      if (!is.finite(target_bound)) {
        return(list(value = NA_real_, status = "not available"))
      }
      if (target_bound <= 0) {
        return(list(value = bound, status = "limited by search range"))
      }
      lo <- min(bound, estimate)
      hi <- max(bound, estimate)
      root <- tryCatch(
        stats::uniroot(function(z) target(z), lower = lo, upper = hi,
                       tol = 1e-8)$root,
        error = function(e) NA_real_
      )
      if (is.finite(root)) {
        list(value = root, status = "ok")
      } else {
        list(value = NA_real_, status = "not available")
      }
    }

    lower <- root_side(-abs(max_abs))
    upper <- root_side(abs(max_abs))
    statuses <- unique(c(lower$status, upper$status))
    status <- if (all(statuses == "ok")) {
      "ok"
    } else {
      paste(statuses, collapse = "; ")
    }
    list(
      lower = lower$value,
      upper = upper$value,
      level = level,
      status = status
    )
  }

  estimate_bias_for_group <- function(idx_rows) {
    eta_sub <- eta_base[idx_rows]
    score_k_sub <- score_k[idx_rows]
    weight_sub <- if (!is.null(weight)) weight[idx_rows] else NULL
    step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL
    slope_idx_sub <- if (!is.null(slope_idx)) slope_idx[idx_rows] else NULL
    nll <- make_bias_nll(eta_sub, score_k_sub, weight_sub,
                         step_idx_sub, slope_idx_sub)
    opt <- tryCatch(
      stats::optimize(nll, interval = c(-max_abs, max_abs)),
      error = function(e) e
    )
    if (inherits(opt, "error")) {
      return(list(value = NA_real_, error = conditionMessage(opt)))
    }
    list(value = opt$minimum, error = NULL)
  }

  iteration_metrics <- function(bias_map) {
    max_resid <- 0
    max_resid_pct <- NA_real_
    for (key in names(groups)) {
      g <- groups[[key]]
      idx_rows <- g$idx
      bias <- bias_map[[key]]
      if (!is.finite(bias) || has_bias_error(bias_error_map[[key]])) next
      eta_sub <- eta_base[idx_rows] + ifelse(is.finite(bias), bias, 0)
      score_k_sub <- score_k[idx_rows]
      step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL
      slope_idx_sub <- if (!is.null(slope_idx)) slope_idx[idx_rows] else NULL
      probs <- if (config$model == "RSM") {
        category_prob_rsm(eta_sub, step_cum)
      } else if (identical(config$model, "GPCM")) {
        category_prob_gpcm(eta_sub, step_cum_mat,
                            criterion_idx = step_idx_sub,
                            slopes = slopes_full,
                            slope_idx = slope_idx_sub)
      } else {
        category_prob_pcm(eta_sub, step_cum_mat, step_idx_sub)
      }
      k_vals <- 0:(ncol(probs) - 1)
      expected_k <- as.vector(probs %*% k_vals)
      expected_score <- prep$rating_min + expected_k
      obs_score <- obs_df$Observed[idx_rows]
      if ("Weight" %in% names(obs_df)) {
        w <- obs_df$Weight[idx_rows]
        obs_score <- obs_score * w
        expected_score <- expected_score * w
      }
      obs_sum <- sum(obs_score, na.rm = TRUE)
      exp_sum <- sum(expected_score, na.rm = TRUE)
      resid_sum <- obs_sum - exp_sum
      if (abs(resid_sum) >= abs(max_resid)) {
        max_resid <- resid_sum
        max_resid_pct <- ifelse(exp_sum != 0, resid_sum / exp_sum * 100, NA_real_)
      }
    }
    list(
      max_resid = max_resid,
      max_resid_pct = max_resid_pct,
      max_resid_categories = NA_real_
    )
  }

  bias_map <- stats::setNames(as.list(rep(0, length(groups))), names(groups))
  bias_error_map <- stats::setNames(as.list(rep(NA_character_, length(groups))), names(groups))

  iter_rows <- list()
  for (it in seq_len(max_iter)) {
    max_change_abs <- 0
    max_change_signed <- 0
    changes <- numeric(0)
    for (key in names(groups)) {
      g <- groups[[key]]
      bias_result <- estimate_bias_for_group(g$idx)
      bias_hat <- bias_result$value
      prev <- bias_map[[key]]
      if (is.finite(bias_hat) && is.finite(prev)) {
        delta <- bias_hat - prev
        if (abs(delta) >= max_change_abs) {
          max_change_abs <- abs(delta)
          max_change_signed <- delta
        }
        changes <- c(changes, abs(delta))
      } else {
        changes <- c(changes, NA_real_)
      }
      bias_map[[key]] <- bias_hat
      bias_error_map[[key]] <- if (has_bias_error(bias_result$error)) bias_result$error else NA_character_
    }
    resid_info <- iteration_metrics(bias_map)
    iter_rows[[it]] <- tibble(
      Iteration = it,
      MaxScoreResidual = resid_info$max_resid,
      MaxScoreResidualPct = resid_info$max_resid_pct,
      MaxScoreResidualCategories = resid_info$max_resid_categories,
      MaxLogitChange = max_change_signed,
      BiasCells = sum(changes > tol, na.rm = TRUE)
    )
    if (max_change_abs < tol) break
  }

  interaction_label <- paste(selected_facets, collapse = " x ")
  interaction_order <- length(selected_facets)
  interaction_mode <- ifelse(interaction_order > 2, "higher_order", "pairwise")

  rows <- list()
  seq_id <- 1
  for (key in names(groups)) {
    g <- groups[[key]]
    idx_rows <- g$idx
    levels <- g$levels
    bias_hat <- bias_map[[key]]
    bias_error <- bias_error_map[[key]]
    bias_ok <- is.finite(bias_hat) && !has_bias_error(bias_error)
    eta_sub <- eta_base[idx_rows]
    score_k_sub <- score_k[idx_rows]
    weight_sub <- if (!is.null(weight)) weight[idx_rows] else rep(1, length(idx_rows))
    step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL
    slope_idx_sub <- if (!is.null(slope_idx)) slope_idx[idx_rows] else NULL

    if (bias_ok) {
      probs <- if (config$model == "RSM") {
        category_prob_rsm(eta_sub + bias_hat, step_cum)
      } else if (identical(config$model, "GPCM")) {
        category_prob_gpcm(eta_sub + bias_hat, step_cum_mat,
                            criterion_idx = step_idx_sub,
                            slopes = slopes_full,
                            slope_idx = slope_idx_sub)
      } else {
        category_prob_pcm(eta_sub + bias_hat, step_cum_mat, step_idx_sub)
      }
      k_vals <- 0:(ncol(probs) - 1)
      expected_k <- as.vector(probs %*% k_vals)
      var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
      var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
      resid_k <- score_k_sub - expected_k
      std_sq <- resid_k^2 / var_k

      # The bias parameter is an additive shift on eta.  Under GPCM,
      # d log L / d b = a * (X - E[X]), so the conditional information
      # for b is a^2 * Var(X), not Var(X).  RSM/PCM are the a = 1 case.
      info_var <- var_k
      if (identical(config$model, "GPCM")) {
        slope_obs_sub <- if (!is.null(slope_idx_sub) && !is.null(slopes_full)) {
          as.numeric(slopes_full[slope_idx_sub])
        } else {
          rep(NA_real_, length(var_k))
        }
        info_var <- ifelse(
          is.finite(slope_obs_sub) & slope_obs_sub > 0,
          (slope_obs_sub^2) * var_k,
          NA_real_
        )
      }
      info <- sum(info_var * weight_sub, na.rm = TRUE)
      se <- ifelse(is.finite(info) && info > 0, 1 / sqrt(info), NA_real_)
      infit <- ifelse(sum(var_k * weight_sub, na.rm = TRUE) > 0,
                      sum(std_sq * var_k * weight_sub, na.rm = TRUE) / sum(var_k * weight_sub, na.rm = TRUE),
                      NA_real_)
      outfit <- ifelse(sum(weight_sub, na.rm = TRUE) > 0,
                       sum(std_sq * weight_sub, na.rm = TRUE) / sum(weight_sub, na.rm = TRUE),
                       NA_real_)
    } else {
      expected_k <- rep(NA_real_, length(score_k_sub))
      se <- NA_real_
      infit <- NA_real_
      outfit <- NA_real_
    }

    lr_chisq <- NA_real_
    lr_df <- NA_real_
    lr_prob <- NA_real_
    profile_ci <- list(
      lower = NA_real_,
      upper = NA_real_,
      level = 0.95,
      status = "not available"
    )
    likelihood_basis <- NA_character_
    if (bias_ok && identical(config$model, "GPCM")) {
      nll_sub <- make_bias_nll(eta_sub, score_k_sub, weight_sub,
                               step_idx_sub, slope_idx_sub)
      nll_hat <- safe_nll_value(nll_sub, bias_hat)
      nll_null <- safe_nll_value(nll_sub, 0)
      if (is.finite(nll_hat) && is.finite(nll_null)) {
        lr_chisq <- max(0, 2 * (nll_null - nll_hat))
        lr_df <- 1
        lr_prob <- stats::pchisq(lr_chisq, df = lr_df, lower.tail = FALSE)
        profile_ci <- profile_bias_ci(
          nll = nll_sub,
          estimate = bias_hat,
          nll_min = nll_hat,
          max_abs = max_abs,
          level = 0.95
        )
        likelihood_basis <- paste(
          "conditional profile likelihood for one additive GPCM bias shift;",
          "theta, steps, slopes, and other facet estimates held fixed"
        )
      }
    }

    obs_slice <- obs_df[idx_rows, , drop = FALSE]
    w_obs <- if ("Weight" %in% names(obs_slice)) obs_slice$Weight else rep(1, nrow(obs_slice))
    obs_score <- sum(obs_slice$Observed * w_obs, na.rm = TRUE)
    exp_score <- if (bias_ok) sum((prep$rating_min + expected_k) * w_obs, na.rm = TRUE) else NA_real_
    obs_count <- sum(w_obs, na.rm = TRUE)
    obs_exp_avg <- ifelse(bias_ok && obs_count > 0 && is.finite(exp_score), (obs_score - exp_score) / obs_count, NA_real_)

    n_obs <- nrow(obs_slice)
    df_t <- ifelse(is.finite(obs_count), max(obs_count - 1, 0), NA_real_)
    t_val <- ifelse(is.finite(bias_hat) && is.finite(se) && se > 0, bias_hat / se, NA_real_)
    p_val <- ifelse(is.finite(t_val) && df_t > 0, 2 * stats::pt(-abs(t_val), df = df_t), NA_real_)

    row <- list(
      Sq = seq_id,
      `Observd Score` = obs_score,
      `Expctd Score` = exp_score,
      `Observd Count` = obs_count,
      `Obs-Exp Average` = obs_exp_avg,
      `Bias Size` = bias_hat,
      `S.E.` = se,
      t = t_val,
      `d.f.` = df_t,
      `Prob.` = p_val,
      InferenceTier = "screening",
      SupportsFormalInference = FALSE,
      FormalInferenceEligible = FALSE,
      PrimaryReportingEligible = FALSE,
      ReportingUse = "screening_only",
      SEBasis = "conditional plug-in information",
      DFBasis = "observed-count minus 1 approximation",
      StatisticLabel = "screening t",
      ProbabilityMetric = "screening tail area",
      Infit = infit,
      Outfit = outfit,
      ObsN = n_obs,
      OptimizationStatus = if (bias_ok) "ok" else "failed",
      OptimizationDetail = if (bias_ok) NA_character_ else bias_error,
      InteractionFacets = interaction_label,
      InteractionOrder = interaction_order,
      InteractionMode = interaction_mode
    )
    if (identical(config$model, "GPCM")) {
      row$`LR ChiSq` <- lr_chisq
      row$`LR d.f.` <- lr_df
      row$`LR Prob.` <- lr_prob
      row$`Profile CI Lower` <- profile_ci$lower
      row$`Profile CI Upper` <- profile_ci$upper
      row$`Profile CI Level` <- profile_ci$level
      row$`Profile CI Status` <- profile_ci$status
      row$`Likelihood Basis` <- likelihood_basis
    }

    for (j in seq_along(selected_facets)) {
      facet_j <- selected_facets[j]
      level_j <- as.character(levels[[facet_j]])
      index_j <- match(level_j, level_map[[facet_j]])
      row[[paste0("Facet", j)]] <- facet_j
      row[[paste0("Facet", j, "_Level")]] <- level_j
      row[[paste0("Facet", j, "_Index")]] <- ifelse(is.na(index_j), NA_real_, index_j)
      row[[paste0("Facet", j, "_Measure")]] <- meas_map[[paste(facet_j, level_j, sep = "||")]]
      row[[paste0("Facet", j, "_SE")]] <- se_map[[paste(facet_j, level_j, sep = "||")]]
    }

    # Backward compatibility for existing 2-way workflows and helper APIs.
    if (interaction_order >= 2) {
      row$FacetA <- row$Facet1
      row$FacetA_Level <- row$Facet1_Level
      row$FacetA_Index <- row$Facet1_Index
      row$FacetA_Measure <- row$Facet1_Measure
      row$FacetA_SE <- row$Facet1_SE

      row$FacetB <- row$Facet2
      row$FacetB_Level <- row$Facet2_Level
      row$FacetB_Index <- row$Facet2_Index
      row$FacetB_Measure <- row$Facet2_Measure
      row$FacetB_SE <- row$Facet2_SE
    }

    rows[[seq_id]] <- tibble::as_tibble(row)
    seq_id <- seq_id + 1
  }

  bias_tbl <- bind_rows(rows)
  if (nrow(bias_tbl) == 0) return(list())

  numeric_cols <- c("Observd Score", "Expctd Score", "Observd Count", "Obs-Exp Average", "Bias Size", "S.E.")
  pop_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    m <- mean(x)
    sqrt(mean((x - m)^2))
  }
  mean_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)))
  sd_pop_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ pop_sd(.x)))
  sd_sample_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ stats::sd(.x, na.rm = TRUE)))
  summary_tbl <- bind_rows(
    mean_row,
    sd_pop_row,
    sd_sample_row
  ) |>
    mutate(
      Statistic = c(paste0("Mean (Count: ", nrow(bias_tbl), ")"), "S.D. (Population)", "S.D. (Sample)"),
      InteractionFacets = interaction_label,
      InteractionOrder = interaction_order,
      InteractionMode = interaction_mode,
      MixedSign = orientation_review$mixed_sign,
      .before = 1
    )

  se_bias <- bias_tbl$`S.E.`
  bias_vals <- bias_tbl$`Bias Size`
  w_chi <- ifelse(is.finite(se_bias) & se_bias > 0, 1 / (se_bias^2), NA_real_)
  ok <- is.finite(w_chi) & is.finite(bias_vals)
  fixed_chi <- if (sum(ok) >= 2) {
    sum(w_chi[ok] * bias_vals[ok]^2) - (sum(w_chi[ok] * bias_vals[ok])^2) / sum(w_chi[ok])
  } else {
    NA_real_
  }
  fixed_df <- max(nrow(bias_tbl) - 1, 0)
  fixed_prob <- ifelse(is.finite(fixed_chi) && fixed_df > 0, 1 - stats::pchisq(fixed_chi, df = fixed_df), NA_real_)
  optimization_failures <- purrr::imap_dfr(groups, function(g, key) {
    err <- bias_error_map[[key]]
    if (!has_bias_error(err)) {
      return(tibble::tibble())
    }
    out <- tibble::tibble(
      InteractionKey = key,
      Error = err
    )
    for (facet_j in selected_facets) {
      out[[facet_j]] <- as.character(g$levels[[facet_j]])
    }
    out
  })

  chi_tbl <- tibble(
    FixedChiSq = fixed_chi,
    FixedDF = fixed_df,
    FixedProb = fixed_prob,
    InferenceTier = "screening",
    SupportsFormalInference = FALSE,
    FormalInferenceEligible = FALSE,
    PrimaryReportingEligible = FALSE,
    ReportingUse = "screening_only",
    TestBasis = "conditional plug-in heterogeneity screen",
    InteractionFacets = interaction_label,
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode
  )

  list(
    facet_a = selected_facets[1],
    facet_b = selected_facets[2],
    interaction_facets = selected_facets,
    interaction_order = interaction_order,
    interaction_mode = interaction_mode,
    table = bias_tbl,
    summary = summary_tbl,
    chi_sq = chi_tbl,
    iteration = bind_rows(iter_rows),
    orientation_review = orientation_review$table,
    mixed_sign = orientation_review$mixed_sign,
    direction_note = orientation_review$direction_note,
    recommended_action = orientation_review$recommended_action,
    inference_tier = "screening",
    optimization_failures = optimization_failures
  )
}

calc_bias_pairwise <- function(bias_tbl, target_facet, context_facet) {
  if (is.null(bias_tbl) || nrow(bias_tbl) == 0) return(tibble())
  if ("InteractionOrder" %in% names(bias_tbl)) {
    ord <- suppressWarnings(as.numeric(bias_tbl$InteractionOrder[1]))
    if (is.finite(ord) && ord != 2) return(tibble())
  }

  use_a <- bias_tbl$FacetA[1] == target_facet
  use_b <- bias_tbl$FacetB[1] == target_facet
  if (!(use_a || use_b)) return(tibble())

  target_prefix <- if (use_a) "FacetA" else "FacetB"
  context_prefix <- if (use_a) "FacetB" else "FacetA"

  sub <- bias_tbl |>
    filter(.data[[context_prefix]] %in% context_facet)
  if (nrow(sub) == 0) sub <- bias_tbl

  rows <- list()
  for (tgt_level in unique(sub[[paste0(target_prefix, "_Level")]])) {
    group_df <- sub |> filter(.data[[paste0(target_prefix, "_Level")]] == tgt_level)
    contexts <- unique(group_df[[paste0(context_prefix, "_Level")]])
    if (length(contexts) < 2) next
    ctx_pairs <- combn(contexts, 2, simplify = FALSE)
    for (pair in ctx_pairs) {
      c1 <- pair[1]
      c2 <- pair[2]
      r1 <- group_df |> filter(.data[[paste0(context_prefix, "_Level")]] == c1) |> slice(1)
      r2 <- group_df |> filter(.data[[paste0(context_prefix, "_Level")]] == c2) |> slice(1)

      tgt_measure <- r1[[paste0(target_prefix, "_Measure")]]
      tgt_se <- r1[[paste0(target_prefix, "_SE")]]
      tgt_index <- r1[[paste0(target_prefix, "_Index")]]

      bias1 <- r1$`Bias Size`
      bias2 <- r2$`Bias Size`
      bias_se1 <- r1$`S.E.`
      bias_se2 <- r2$`S.E.`

      local1 <- ifelse(is.finite(tgt_measure) & is.finite(bias1), tgt_measure + bias1, NA_real_)
      local2 <- ifelse(is.finite(tgt_measure) & is.finite(bias2), tgt_measure + bias2, NA_real_)
      se1 <- ifelse(is.finite(tgt_se) & is.finite(bias_se1), sqrt(tgt_se^2 + bias_se1^2), NA_real_)
      se2 <- ifelse(is.finite(tgt_se) & is.finite(bias_se2), sqrt(tgt_se^2 + bias_se2^2), NA_real_)

      contrast <- ifelse(is.finite(local1) & is.finite(local2), local1 - local2, NA_real_)
      se_contrast <- ifelse(
        is.finite(bias_se1) & is.finite(bias_se2),
        sqrt(bias_se1^2 + bias_se2^2),
        NA_real_
      )
      t_val <- ifelse(is.finite(contrast) & is.finite(se_contrast) & se_contrast > 0, contrast / se_contrast, NA_real_)

      n1 <- ifelse(is.finite(r1$ObsN), r1$ObsN, 0)
      n2 <- ifelse(is.finite(r2$ObsN), r2$ObsN, 0)
      df_t <- welch_satterthwaite_df(c(bias_se1^2, bias_se2^2), c(n1 - 1, n2 - 1))
      p_val <- ifelse(is.finite(t_val) & is.finite(df_t) & df_t > 0,
                      2 * stats::pt(-abs(t_val), df = df_t),
                      NA_real_)

      rows[[length(rows) + 1]] <- tibble(
        Target = tgt_level,
        `Target N` = tgt_index,
        `Target Measure` = tgt_measure,
        `Target S.E.` = tgt_se,
        Context1 = c1,
        `Context1 N` = r1[[paste0(context_prefix, "_Index")]],
        `Local Measure1` = local1,
        SE1 = se1,
        `Obs-Exp Avg1` = r1$`Obs-Exp Average`,
        Count1 = r1$`Observd Count`,
        ObsN1 = n1,
        Context2 = c2,
        `Context2 N` = r2[[paste0(context_prefix, "_Index")]],
        `Local Measure2` = local2,
        SE2 = se2,
        `Obs-Exp Avg2` = r2$`Obs-Exp Average`,
        Count2 = r2$`Observd Count`,
        ObsN2 = n2,
        Contrast = contrast,
        SE = se_contrast,
        t = t_val,
        `d.f.` = df_t,
        `Prob.` = p_val,
        InferenceTier = "screening",
        SupportsFormalInference = FALSE,
        FormalInferenceEligible = FALSE,
        PrimaryReportingEligible = FALSE,
        ReportingUse = "screening_only",
        ContrastBasis = "difference between local target measures across contexts (target term cancels to a bias contrast)",
        SEBasis = "combined context-specific bias standard errors",
        StatisticLabel = "Bias-contrast Welch screening t",
        ProbabilityMetric = "screening tail area",
        DFBasis = "Welch-Satterthwaite approximation"
      )
    }
  }
  bind_rows(rows)
}

extract_anchor_tables <- function(config) {
  facet_names <- c("Person", config$facet_names)
  anchor_tbl <- purrr::map_dfr(facet_names, function(facet) {
    spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
    anchors <- spec$anchors
    if (is.null(anchors)) return(tibble())
    df <- tibble(Facet = facet, Level = names(anchors), Anchor = as.numeric(anchors)) |>
      filter(is.finite(Anchor))
    if (nrow(df) == 0) return(tibble())
    df |>
      mutate(Source = ifelse(facet %in% config$dummy_facets, "Dummy facet", "Anchor"))
  })

  group_tbl <- purrr::map_dfr(facet_names, function(facet) {
    spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
    groups <- spec$groups
    if (is.null(groups)) return(tibble())
    df <- tibble(Facet = facet, Level = names(groups), Group = as.character(groups)) |>
      filter(!is.na(Group), Group != "")
    if (nrow(df) == 0) return(tibble())
    group_values <- spec$group_values
    df |>
      mutate(GroupValue = ifelse(Group %in% names(group_values), group_values[Group], NA_real_))
  })

  list(anchors = anchor_tbl, groups = group_tbl)
}

calc_facets_chisq <- function(measure_df) {
  if (is.null(measure_df) || nrow(measure_df) == 0) return(tibble())
  measure_df |>
    group_by(Facet) |>
    summarize(
      Levels = n(),
      MeanMeasure = mean(Estimate, na.rm = TRUE),
      SD = sd(Estimate, na.rm = TRUE),
      EstimateVec = list(Estimate),
      SEVec = list(SE),
      FixedChiSq = {
        w <- ifelse(is.finite(SE) & SE > 0, 1 / (SE^2), NA_real_)
        d <- Estimate
        ok <- is.finite(w) & is.finite(d)
        if (sum(ok) < 2) NA_real_ else sum(w[ok] * d[ok]^2) - (sum(w[ok] * d[ok])^2) / sum(w[ok])
      },
      FixedDF = pmax(Levels - 1, 0),
      RandomVar = {
        d <- Estimate
        se2 <- SE^2
        ok <- is.finite(d) & is.finite(se2)
        if (sum(ok) < 2) NA_real_ else (sum((d[ok] - mean(d[ok]))^2) / (sum(ok) - 1)) - (sum(se2[ok]) / sum(ok))
      },
      .groups = "drop"
    ) |>
    mutate(
      FixedProb = ifelse(is.finite(FixedChiSq) & FixedDF > 0,
                         1 - stats::pchisq(FixedChiSq, df = FixedDF),
                         NA_real_),
      RandomVar = ifelse(is.finite(RandomVar) & RandomVar > 0, RandomVar, NA_real_)
    ) |>
    rowwise() |>
    mutate(
      RandomChiSq = {
        if (!is.finite(RandomVar) || RandomVar <= 0) {
          NA_real_
        } else {
          d <- unlist(EstimateVec)
          se <- unlist(SEVec)
          w <- ifelse(is.finite(se) & se > 0, 1 / (RandomVar + se^2), NA_real_)
          ok <- is.finite(w) & is.finite(d)
          if (sum(ok) < 2) NA_real_ else sum(w[ok] * d[ok]^2) - (sum(w[ok] * d[ok])^2) / sum(w[ok])
        }
      },
      RandomDF = pmax(Levels - 2, 0),
      RandomProb = ifelse(is.finite(RandomChiSq) & RandomDF > 0,
                          1 - stats::pchisq(RandomChiSq, df = RandomDF),
                          NA_real_)
    ) |>
    ungroup() |>
    select(-EstimateVec, -SEVec)
}

build_param_slices <- function(sizes) {
  out <- list()
  idx <- 1L
  for (nm in names(sizes)) {
    k <- as.integer(sizes[[nm]] %||% 0L)
    if (k > 0L) {
      out[[nm]] <- idx:(idx + k - 1L)
      idx <- idx + k
    } else {
      out[[nm]] <- integer(0)
    }
  }
  out
}

constraint_jacobian <- function(spec) {
  n_levels <- length(spec$levels %||% character(0))
  n_params <- as.integer(spec$n_params %||% 0L)
  out <- matrix(0, nrow = n_levels, ncol = n_params)
  rownames(out) <- as.character(spec$levels %||% character(0))
  if (n_levels == 0L || n_params == 0L) {
    return(out)
  }

  baseline <- expand_facet_with_constraints(numeric(n_params), spec)
  for (j in seq_len(n_params)) {
    free_vec <- numeric(n_params)
    free_vec[j] <- 1
    out[, j] <- expand_facet_with_constraints(free_vec, spec) - baseline
  }
  out
}

sum_zero_jacobian <- function(n_levels) {
  n_levels <- as.integer(n_levels %||% 0L)
  n_params <- sum_zero_param_count(n_levels)
  out <- matrix(0, nrow = max(n_levels, 0L), ncol = n_params)
  if (n_levels <= 1L || n_params == 0L) {
    return(out)
  }
  out[seq_len(n_params), seq_len(n_params)] <- diag(n_params)
  out[n_levels, ] <- -1
  out
}

block_diag_sum_zero_jacobian <- function(n_blocks, n_levels) {
  n_blocks <- as.integer(n_blocks %||% 0L)
  n_levels <- as.integer(n_levels %||% 0L)
  one <- sum_zero_jacobian(n_levels)
  n_params <- ncol(one)
  out <- matrix(0, nrow = max(n_blocks, 0L) * max(n_levels, 0L),
                ncol = max(n_blocks, 0L) * n_params)
  if (n_blocks <= 0L || n_levels <= 0L || n_params == 0L) {
    return(out)
  }
  for (i in seq_len(n_blocks)) {
    row_idx <- ((i - 1L) * n_levels + 1L):(i * n_levels)
    col_idx <- ((i - 1L) * n_params + 1L):(i * n_params)
    out[row_idx, col_idx] <- one
  }
  out
}

symmetrize_matrix <- function(mat) {
  if (is.null(mat)) return(NULL)
  (mat + t(mat)) / 2
}

invert_information_matrix <- function(info_mat) {
  if (is.null(info_mat)) {
    return(list(cov = NULL, regularized = FALSE, rank = 0L))
  }

  info_mat <- suppressWarnings(as.matrix(info_mat))
  if (!is.numeric(info_mat) || length(info_mat) == 0L || any(!is.finite(info_mat))) {
    return(list(cov = NULL, regularized = FALSE, rank = 0L))
  }

  info_mat <- symmetrize_matrix(info_mat)
  eig <- tryCatch(
    eigen(info_mat, symmetric = TRUE),
    error = function(e) NULL
  )
  if (is.null(eig)) {
    return(list(cov = NULL, regularized = FALSE, rank = 0L))
  }

  values <- eig$values
  max_val <- suppressWarnings(max(abs(values), na.rm = TRUE))
  tol <- if (is.finite(max_val) && max_val > 0) {
    max_val * sqrt(.Machine$double.eps)
  } else {
    sqrt(.Machine$double.eps)
  }
  rank <- sum(values > tol)
  regularized <- any(values <= tol)
  safe_values <- pmax(values, tol)
  cov_mat <- eig$vectors %*% diag(1 / safe_values, nrow = length(safe_values)) %*% t(eig$vectors)
  cov_mat <- symmetrize_matrix(cov_mat)

  list(
    cov = cov_mat,
    regularized = regularized,
    rank = rank
  )
}

compute_mml_parameter_covariance <- function(res) {
  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  if (!identical(method, "MML")) {
    return(list(
      cov = NULL,
      hessian = NULL,
      sizes = build_param_sizes(res$config),
      param_slices = build_param_slices(build_param_sizes(res$config)),
      status = "not_applicable",
      detail = "MML observed-information covariance is only computed for MML fits.",
      regularized = FALSE,
      rank = 0L
    ))
  }

  if (is.null(res$opt$par) || length(res$opt$par) == 0L) {
    sizes <- build_param_sizes(res$config)
    return(list(
      cov = NULL,
      hessian = NULL,
      sizes = sizes,
      param_slices = build_param_slices(sizes),
      status = "fallback",
      detail = "Optimized parameter vector was not available; observed-information covariance was not computed.",
      regularized = FALSE,
      rank = 0L
    ))
  }

  config <- res$config
  sizes <- build_param_sizes(config)
  idx <- build_indices(res$prep, step_facet = config$step_facet,
                       slope_facet = config$slope_facet,
                       interaction_specs = config$interaction_specs)
  quad_points <- max(1L, as.integer(config$estimation_control$quad_points %||% 15L))
  quad <- gauss_hermite_normal(quad_points)
  cache <- make_param_cache(sizes, config, idx, is_mml = TRUE)

  fn <- function(par, idx, config, sizes, quad) {
    cache$ensure(par)
    mfrm_loglik_mml_cached(cache, idx, config, quad)
  }
  gr <- function(par, idx, config, sizes, quad) {
    cache$ensure(par)
    mfrm_grad_mml_cached(cache, idx, config, sizes, quad)
  }

  hess <- tryCatch(
    stats::optimHess(
      par = res$opt$par,
      fn = fn,
      gr = gr,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad
    ),
    error = function(e) e
  )

  if (inherits(hess, "error") || is.null(hess)) {
    msg <- if (inherits(hess, "error")) conditionMessage(hess) else "Unknown Hessian error."
    return(list(
      cov = NULL,
      hessian = NULL,
      sizes = sizes,
      param_slices = build_param_slices(sizes),
      status = "fallback",
      detail = paste0("Observed-information Hessian was unavailable; covariance was not computed. ", msg),
      regularized = FALSE,
      rank = 0L
    ))
  }

  inv_info <- invert_information_matrix(hess)
  cov_free <- inv_info$cov
  if (is.null(cov_free)) {
    return(list(
      cov = NULL,
      hessian = hess,
      sizes = sizes,
      param_slices = build_param_slices(sizes),
      status = "fallback",
      detail = "Observed-information matrix could not be inverted; covariance was not computed.",
      regularized = FALSE,
      rank = inv_info$rank
    ))
  }

  detail <- if (isTRUE(inv_info$regularized)) {
    "MML parameter covariance uses the observed information of the marginal log-likelihood; a near-singular Hessian was regularized during inversion."
  } else {
    "MML parameter covariance uses the observed information of the marginal log-likelihood."
  }

  list(
    cov = cov_free,
    hessian = hess,
    sizes = sizes,
    param_slices = build_param_slices(sizes),
    status = if (isTRUE(inv_info$regularized)) "regularized" else "ok",
    detail = detail,
    regularized = isTRUE(inv_info$regularized),
    rank = inv_info$rank
  )
}

covariance_diag_se <- function(cov_mat) {
  if (is.null(cov_mat)) return(numeric(0))
  diag_var <- diag(symmetrize_matrix(cov_mat))
  diag_var[diag_var < 0 & abs(diag_var) < 1e-10] <- 0
  ifelse(diag_var >= 0, sqrt(diag_var), NA_real_)
}

compute_mml_facet_model_se <- function(res, covariance = NULL) {
  covariance <- covariance %||% compute_mml_parameter_covariance(res)
  if (!covariance$status %in% c("ok", "regularized") || is.null(covariance$cov)) {
    return(list(
      table = tibble(),
      status = covariance$status,
      detail = covariance$detail
    ))
  }

  config <- res$config
  cov_free <- covariance$cov
  param_slices <- covariance$param_slices
  facet_tbl <- purrr::map_dfr(config$facet_names, function(facet) {
    spec <- config$facet_specs[[facet]]
    levels <- as.character(spec$levels %||% res$prep$levels[[facet]] %||% character(0))
    slice <- param_slices[[facet]] %||% integer(0)

    if (length(levels) == 0L) {
      return(tibble(Facet = character(0), Level = character(0), ModelSE = numeric(0)))
    }

    if (length(slice) == 0L || is.null(spec) || (spec$n_params %||% 0L) == 0L) {
      return(tibble(
        Facet = facet,
        Level = levels,
        ModelSE = NA_real_
      ))
    }

    jac <- constraint_jacobian(spec)
    cov_block <- cov_free[slice, slice, drop = FALSE]
    cov_expanded <- symmetrize_matrix(jac %*% cov_block %*% t(jac))

    tibble(
      Facet = facet,
      Level = levels,
      ModelSE = covariance_diag_se(cov_expanded)
    )
  })

  detail <- if (identical(covariance$status, "regularized")) {
    "MML facet ModelSE values use the observed information of the marginal log-likelihood; a near-singular Hessian was regularized during inversion."
  } else {
    "MML facet ModelSE values use the observed information of the marginal log-likelihood."
  }

  list(
    table = facet_tbl,
    status = covariance$status,
    detail = detail
  )
}

compute_mml_structural_parameter_se <- function(res,
                                                covariance = NULL,
                                                ci_level = 0.95) {
  step_tbl <- as.data.frame(res$steps %||% tibble(), stringsAsFactors = FALSE)
  slope_tbl <- as.data.frame(res$slopes %||% tibble(), stringsAsFactors = FALSE)
  covariance <- covariance %||% compute_mml_parameter_covariance(res)
  z <- stats::qnorm((1 + ci_level) / 2)

  add_unavailable <- function(tbl, kind) {
    if (!is.data.frame(tbl) || nrow(tbl) == 0L) return(tbl)
    tbl$SE <- NA_real_
    tbl$CI_Lower <- NA_real_
    tbl$CI_Upper <- NA_real_
    tbl$CI_Level <- ci_level
    tbl$SE_Method <- if (identical(covariance$status, "not_applicable")) {
      "not applicable"
    } else {
      "not available"
    }
    tbl$SE_Status <- covariance$status
    tbl$SE_Detail <- covariance$detail
    if (identical(kind, "slope")) {
      tbl$LogSE <- NA_real_
      tbl$LogCI_Lower <- NA_real_
      tbl$LogCI_Upper <- NA_real_
    }
    tbl
  }

  if (!covariance$status %in% c("ok", "regularized") || is.null(covariance$cov)) {
    return(list(
      steps = add_unavailable(step_tbl, "step"),
      slopes = add_unavailable(slope_tbl, "slope"),
      summary = tibble(
        Component = c("steps", "slopes"),
        Available = FALSE,
        Status = covariance$status,
        Detail = covariance$detail
      ),
      status = covariance$status,
      detail = covariance$detail
    ))
  }

  config <- res$config
  cov_free <- covariance$cov
  param_slices <- covariance$param_slices
  status <- covariance$status
  se_method <- if (identical(status, "regularized")) {
    "Observed information (MML; regularized Hessian)"
  } else {
    "Observed information (MML)"
  }

  if (nrow(step_tbl) > 0L) {
    step_slice <- param_slices$steps %||% integer(0)
    n_steps <- max(config$n_cat - 1L, 0L)
    if (length(step_slice) > 0L && n_steps > 0L) {
      if (identical(config$model, "RSM")) {
        jac <- sum_zero_jacobian(n_steps)
      } else {
        n_step_levels <- length(config$facet_levels[[config$step_facet]] %||% character(0))
        jac <- block_diag_sum_zero_jacobian(n_step_levels, n_steps)
      }
      cov_step <- symmetrize_matrix(jac %*% cov_free[step_slice, step_slice, drop = FALSE] %*% t(jac))
      step_se <- covariance_diag_se(cov_step)
      step_tbl$SE <- if (length(step_se) == nrow(step_tbl)) step_se else NA_real_
    } else {
      step_tbl$SE <- NA_real_
    }
    step_tbl$CI_Lower <- ifelse(is.finite(step_tbl$SE), step_tbl$Estimate - z * step_tbl$SE, NA_real_)
    step_tbl$CI_Upper <- ifelse(is.finite(step_tbl$SE), step_tbl$Estimate + z * step_tbl$SE, NA_real_)
    step_tbl$CI_Level <- ci_level
    step_tbl$SE_Method <- se_method
    step_tbl$SE_Status <- status
    step_tbl$SE_Detail <- covariance$detail
  }

  if (nrow(slope_tbl) > 0L && identical(config$model, "GPCM")) {
    slope_slice <- param_slices$log_slopes %||% integer(0)
    n_slopes <- nrow(slope_tbl)
    if (length(slope_slice) > 0L && n_slopes > 0L) {
      jac <- sum_zero_jacobian(n_slopes)
      cov_log <- symmetrize_matrix(jac %*% cov_free[slope_slice, slope_slice, drop = FALSE] %*% t(jac))
      log_se <- covariance_diag_se(cov_log)
      slope_vals <- as.numeric(slope_tbl$Estimate)
      cov_slope <- symmetrize_matrix(diag(slope_vals, nrow = length(slope_vals)) %*%
                                       cov_log %*%
                                       diag(slope_vals, nrow = length(slope_vals)))
      slope_se <- covariance_diag_se(cov_slope)
      slope_tbl$LogSE <- if (length(log_se) == nrow(slope_tbl)) log_se else NA_real_
      slope_tbl$SE <- if (length(slope_se) == nrow(slope_tbl)) slope_se else NA_real_
    } else {
      slope_tbl$LogSE <- NA_real_
      slope_tbl$SE <- NA_real_
    }
    slope_tbl$LogCI_Lower <- ifelse(is.finite(slope_tbl$LogSE), slope_tbl$LogEstimate - z * slope_tbl$LogSE, NA_real_)
    slope_tbl$LogCI_Upper <- ifelse(is.finite(slope_tbl$LogSE), slope_tbl$LogEstimate + z * slope_tbl$LogSE, NA_real_)
    slope_tbl$CI_Lower <- ifelse(is.finite(slope_tbl$LogCI_Lower), exp(slope_tbl$LogCI_Lower), NA_real_)
    slope_tbl$CI_Upper <- ifelse(is.finite(slope_tbl$LogCI_Upper), exp(slope_tbl$LogCI_Upper), NA_real_)
    slope_tbl$CI_Level <- ci_level
    slope_tbl$SE_Method <- se_method
    slope_tbl$SE_Status <- status
    slope_tbl$SE_Detail <- covariance$detail
  }

  list(
    steps = tibble::as_tibble(step_tbl),
    slopes = tibble::as_tibble(slope_tbl),
    summary = tibble(
      Component = c("steps", "slopes"),
      Available = c(
        nrow(step_tbl) > 0L && "SE" %in% names(step_tbl) && any(is.finite(step_tbl$SE)),
        nrow(slope_tbl) > 0L && "SE" %in% names(slope_tbl) && any(is.finite(slope_tbl$SE))
      ),
      Status = status,
      Detail = covariance$detail
    ),
    status = status,
    detail = covariance$detail,
    covariance = list(
      parameters = length(res$opt$par %||% numeric(0)),
      rank = covariance$rank,
      regularized = covariance$regularized
    )
  )
}

build_measure_se_table <- function(res, obs_df, facet_cols, fit_tbl, covariance = NULL) {
  approx_tbl <- calc_facet_se(obs_df, facet_cols) |>
    rename(ApproxSE = SE)

  base_tbl <- bind_rows(
    res$facets$person |>
      transmute(Facet = "Person", Level = as.character(Person)),
    res$facets$others |>
      transmute(Facet = as.character(Facet), Level = as.character(Level))
  ) |>
    distinct(Facet, Level)

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  converged <- isTRUE(as.logical(res$summary$Converged[1] %||% FALSE))
  fallback_label <- if (identical(method, "MML")) {
    "Fallback observation-table information"
  } else {
    "Observation-table information"
  }

  person_model_tbl <- if (identical(method, "MML") && "SD" %in% names(res$facets$person)) {
    res$facets$person |>
      transmute(
        Facet = "Person",
        Level = as.character(Person),
        ModelSE = as.numeric(SD),
        SE_Method = "Posterior SD (EAP)"
      )
  } else {
    tibble(
      Facet = character(0),
      Level = character(0),
      ModelSE = numeric(0),
      SE_Method = character(0)
    )
  }

  facet_model_bundle <- compute_mml_facet_model_se(res, covariance = covariance)
  facet_model_tbl <- facet_model_bundle$table |>
    mutate(SE_Method = "Observed information (MML)")

  base_tbl |>
    left_join(approx_tbl, by = c("Facet", "Level")) |>
    left_join(
      bind_rows(person_model_tbl, facet_model_tbl) |>
        rename(ModelSE_Raw = ModelSE, SE_Method_Raw = SE_Method),
      by = c("Facet", "Level")
    ) |>
    left_join(
      fit_tbl |>
        select("Facet", "Level", "Infit"),
      by = c("Facet", "Level")
    ) |>
    mutate(
      PrecisionTier = dplyr::case_when(
        identical(method, "MML") & !is.na(.data$SE_Method_Raw) ~ "model_based",
        identical(method, "MML") ~ "hybrid",
        TRUE ~ "exploratory"
      ),
      Converged = converged,
      SupportsFormalInference = .data$PrecisionTier == "model_based" & .data$Converged,
      ModelSE = dplyr::coalesce(.data$ModelSE_Raw, .data$ApproxSE),
      RealSE = ifelse(
        is.finite(.data$ModelSE),
        .data$ModelSE * sqrt(pmax(dplyr::coalesce(.data$Infit, 1), 1)),
        NA_real_
      ),
      SE = .data$ModelSE,
      SE_Method = dplyr::coalesce(.data$SE_Method_Raw, fallback_label),
      SEUse = dplyr::case_when(
        .data$PrecisionTier == "model_based" & .data$Converged ~ "primary_reporting",
        .data$PrecisionTier == "model_based" ~ "review_before_reporting",
        .data$PrecisionTier == "hybrid" ~ "review_before_reporting",
        TRUE ~ "screening_only"
      ),
      CIBasis = dplyr::case_when(
        .data$PrecisionTier == "model_based" & .data$Converged ~ "Normal interval from model-based SE",
        .data$PrecisionTier == "model_based" ~ "Normal interval from model-based SE; optimizer convergence review required",
        .data$PrecisionTier == "hybrid" ~ "Normal interval from fallback observation-table SE",
        TRUE ~ "Normal interval from exploratory observation-table SE"
      ),
      CIUse = dplyr::case_when(
        .data$PrecisionTier == "model_based" & .data$Converged ~ "primary_reporting",
        .data$PrecisionTier == "model_based" ~ "review_before_reporting",
        .data$PrecisionTier == "hybrid" ~ "review_before_reporting",
        TRUE ~ "screening_only"
      )
    ) |>
    select("Facet", "Level", "N", "SE", "ModelSE", "RealSE", "SE_Method",
           "Converged",
           "PrecisionTier", "SupportsFormalInference", "SEUse", "CIBasis", "CIUse") |>
    arrange(.data$Facet, .data$Level) |>
    structure(mml_se_status = facet_model_bundle$status, mml_se_detail = facet_model_bundle$detail)
}

# Separation, reliability, and strata indices used in package diagnostics.
# Conventions follow Wright & Masters (1982) "Rating Scale Analysis"
# (Chapter 5) and Linacre's Facets/Winsteps documentation:
#
#   Separation  G = TrueSD / RMSE
#   Reliability R = TrueVariance / ObservedVariance = G^2 / (1 + G^2)
#   Strata      H = (4 * G + 1) / 3
#
# Model statistics use ModelSE; fit-adjusted statistics use RealSE
# when available.
summarize_precision_basis <- function(df, se_col, distribution_basis = c("sample", "population")) {
  distribution_basis <- match.arg(distribution_basis)

  est <- suppressWarnings(as.numeric(df$Estimate))
  se_vals <- if (se_col %in% names(df)) suppressWarnings(as.numeric(df[[se_col]])) else rep(NA_real_, nrow(df))
  ok_est <- is.finite(est)
  ok_se <- is.finite(se_vals) & se_vals >= 0
  n_est <- sum(ok_est)

  observed_mean <- if (n_est > 0L) mean(est[ok_est]) else NA_real_
  observed_var <- if (n_est == 0L) {
    NA_real_
  } else if (identical(distribution_basis, "sample")) {
    if (n_est >= 2L) stats::var(est[ok_est]) else NA_real_
  } else {
    mean((est[ok_est] - observed_mean)^2)
  }
  observed_sd <- if (is.finite(observed_var)) sqrt(observed_var) else NA_real_
  rmse <- if (any(ok_se)) sqrt(mean(se_vals[ok_se]^2, na.rm = TRUE)) else NA_real_
  error_var <- if (is.finite(rmse)) rmse^2 else NA_real_
  true_var <- if (is.finite(observed_var) && is.finite(error_var)) {
    pmax(observed_var - error_var, 0)
  } else {
    NA_real_
  }
  true_sd <- if (is.finite(true_var)) sqrt(true_var) else NA_real_
  separation <- if (is.finite(rmse) && rmse > 0 && is.finite(true_sd)) true_sd / rmse else NA_real_
  reliability <- if (is.finite(observed_var) && observed_var > 0 && is.finite(true_var)) true_var / observed_var else NA_real_
  strata <- if (is.finite(separation)) (4 * separation + 1) / 3 else NA_real_

  list(
    ObservedMean = observed_mean,
    ObservedSD = observed_sd,
    RMSE = rmse,
    TrueSD = true_sd,
    ObservedVariance = observed_var,
    ErrorVariance = error_var,
    TrueVariance = true_var,
    Separation = separation,
    Strata = strata,
    Reliability = reliability,
    SEAvailable = sum(ok_se),
    MeanSE = if (any(ok_se)) mean(se_vals[ok_se], na.rm = TRUE) else NA_real_,
    MedianSE = if (any(ok_se)) stats::median(se_vals[ok_se], na.rm = TRUE) else NA_real_
  )
}

build_facet_precision_summary <- function(measure_df, variability_tbl = NULL) {
  if (is.null(measure_df) || nrow(measure_df) == 0) return(tibble())

  variability_tbl <- as.data.frame(variability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  variability_keep <- intersect(
    c("Facet", "FixedChiSq", "FixedDF", "FixedProb", "RandomVar", "RandomChiSq", "RandomDF", "RandomProb"),
    names(variability_tbl)
  )
  variability_tbl <- variability_tbl[, variability_keep, drop = FALSE]

  mode_map <- c(model = "ModelSE", fit_adjusted = "RealSE")
  rows <- list()
  row_id <- 1L

  facets <- unique(as.character(measure_df$Facet))
  facets <- facets[!is.na(facets) & nzchar(facets)]
  for (facet in facets) {
    sub <- measure_df[as.character(measure_df$Facet) == facet, , drop = FALSE]
    for (mode_name in names(mode_map)) {
      se_col <- mode_map[[mode_name]]
      if (!se_col %in% names(sub)) next
      for (basis_name in c("sample", "population")) {
        stats <- summarize_precision_basis(sub, se_col = se_col, distribution_basis = basis_name)
        rows[[row_id]] <- tibble(
          Facet = facet,
          Levels = nrow(sub),
          DistributionBasis = basis_name,
          SEMode = mode_name,
          SEColumn = se_col,
          ObservedMean = stats$ObservedMean,
          ObservedSD = stats$ObservedSD,
          RMSE = stats$RMSE,
          TrueSD = stats$TrueSD,
          ObservedVariance = stats$ObservedVariance,
          ErrorVariance = stats$ErrorVariance,
          TrueVariance = stats$TrueVariance,
          Separation = stats$Separation,
          Strata = stats$Strata,
          Reliability = stats$Reliability,
          SEAvailable = stats$SEAvailable,
          MeanSE = stats$MeanSE,
          MedianSE = stats$MedianSE,
          MeanInfit = if ("Infit" %in% names(sub)) mean(sub$Infit, na.rm = TRUE) else NA_real_,
          MeanOutfit = if ("Outfit" %in% names(sub)) mean(sub$Outfit, na.rm = TRUE) else NA_real_
        )
        row_id <- row_id + 1L
      }
    }
  }

  out <- bind_rows(rows)
  if (nrow(out) == 0) {
    return(tibble())
  }
  if (nrow(variability_tbl) > 0 && "Facet" %in% names(variability_tbl)) {
    out <- left_join(out, variability_tbl, by = "Facet")
  }
  out |>
    arrange(.data$Facet, .data$DistributionBasis, .data$SEMode)
}

build_precision_profile <- function(res, measure_df, reliability_tbl, facet_precision_tbl) {
  measure_df <- as.data.frame(measure_df %||% data.frame(), stringsAsFactors = FALSE)
  reliability_tbl <- as.data.frame(reliability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  facet_precision_tbl <- as.data.frame(facet_precision_tbl %||% data.frame(), stringsAsFactors = FALSE)

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  method <- ifelse(identical(method, "JMLE"), "JML", method)
  converged <- isTRUE(as.logical(res$summary$Converged[1] %||% FALSE))

  person_labels <- character(0)
  nonperson_labels <- character(0)
  all_labels <- character(0)
  if (nrow(measure_df) > 0 && "SE_Method" %in% names(measure_df) && "Facet" %in% names(measure_df)) {
    all_labels <- unique(as.character(stats::na.omit(measure_df$SE_Method)))
    person_labels <- unique(as.character(stats::na.omit(measure_df$SE_Method[measure_df$Facet == "Person"])))
    nonperson_labels <- unique(as.character(stats::na.omit(measure_df$SE_Method[measure_df$Facet != "Person"])))
  }

  has_fallback <- identical(method, "MML") &&
    any(grepl("Fallback observation-table information", all_labels))
  precision_tier <- if (!identical(method, "MML")) {
    "exploratory"
  } else if (has_fallback) {
    "hybrid"
  } else {
    "model_based"
  }
  supports_formal <- identical(precision_tier, "model_based") && converged

  has_fit_adjusted <- nrow(measure_df) > 0 &&
    "RealSE" %in% names(measure_df) &&
    any(is.finite(suppressWarnings(as.numeric(measure_df$RealSE))))

  coverage_complete <- FALSE
  if (nrow(facet_precision_tbl) > 0 &&
      all(c("Facet", "DistributionBasis", "SEMode") %in% names(facet_precision_tbl))) {
    precision_facets <- unique(as.character(facet_precision_tbl$Facet))
    precision_facets <- precision_facets[!is.na(precision_facets) & nzchar(precision_facets)]
    if (length(precision_facets) > 0) {
      expected <- expand.grid(
        DistributionBasis = c("sample", "population"),
        SEMode = c("model", "fit_adjusted"),
        stringsAsFactors = FALSE
      )
      coverage_complete <- all(vapply(precision_facets, function(facet_name) {
        sub <- facet_precision_tbl[facet_precision_tbl$Facet == facet_name, c("DistributionBasis", "SEMode"), drop = FALSE]
        all(apply(expected, 1, function(row) {
          any(sub$DistributionBasis == row[["DistributionBasis"]] & sub$SEMode == row[["SEMode"]])
        }))
      }, logical(1)))
    }
  }

  tibble(
    Method = method,
    Converged = converged,
    PrecisionTier = precision_tier,
    SupportsFormalInference = supports_formal,
    HasFallbackSE = has_fallback,
    PersonSEBasis = if (length(person_labels) > 0) paste(sort(person_labels), collapse = "; ") else NA_character_,
    NonPersonSEBasis = if (length(nonperson_labels) > 0) paste(sort(nonperson_labels), collapse = "; ") else NA_character_,
    CIBasis = if (identical(precision_tier, "model_based")) {
      if (isTRUE(converged)) {
        "Normal interval from model-based SE"
      } else {
        "Normal interval from model-based SE; optimizer convergence review required"
      }
    } else if (identical(precision_tier, "hybrid")) {
      "Normal interval from mixed model-based and fallback SE"
    } else {
      "Normal interval from exploratory SE"
    },
    ReliabilityBasis = if (identical(precision_tier, "model_based")) {
      if (isTRUE(converged)) {
        "Observed variance with model-based and fit-adjusted error bounds"
      } else {
        "Observed variance with model-based and fit-adjusted error bounds; optimizer convergence review required"
      }
    } else if (identical(precision_tier, "hybrid")) {
      "Observed variance with mixed model-based and fallback error bounds"
    } else {
      "Exploratory variance summary with model-based and fit-adjusted error bounds"
    },
    HasFitAdjustedSE = has_fit_adjusted,
    HasSamplePopulationCoverage = coverage_complete,
    RecommendedUse = if (identical(precision_tier, "model_based") && isTRUE(converged)) {
      "Use for primary reporting of SE, CI, and reliability in this package."
    } else if (identical(precision_tier, "model_based")) {
      "This run reached the model-based precision path, but optimizer convergence should be reviewed before primary reporting."
    } else if (identical(precision_tier, "hybrid")) {
      "Use model-based rows for primary reporting, but review levels that fell back to observation-table information before treating the whole run as formal inference."
    } else {
      "Use for screening and calibration triage; confirm formal SE, CI, and reliability with an MML fit."
    }
  )
}

audit_precision_outputs <- function(res, measure_df, reliability_tbl, facet_precision_tbl, precision_profile_tbl = NULL) {
  measure_df <- as.data.frame(measure_df %||% data.frame(), stringsAsFactors = FALSE)
  reliability_tbl <- as.data.frame(reliability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  facet_precision_tbl <- as.data.frame(facet_precision_tbl %||% data.frame(), stringsAsFactors = FALSE)
  precision_profile_tbl <- as.data.frame(precision_profile_tbl %||% data.frame(), stringsAsFactors = FALSE)

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  method <- ifelse(identical(method, "JMLE"), "JML", method)
  converged <- isTRUE(as.logical(res$summary$Converged[1] %||% FALSE))

  finite_model <- if ("ModelSE" %in% names(measure_df)) is.finite(suppressWarnings(as.numeric(measure_df$ModelSE))) else logical(0)
  model_share <- if (length(finite_model) > 0) mean(finite_model) else NA_real_

  real_ok <- TRUE
  if (all(c("ModelSE", "RealSE") %in% names(measure_df))) {
    model_vals <- suppressWarnings(as.numeric(measure_df$ModelSE))
    real_vals <- suppressWarnings(as.numeric(measure_df$RealSE))
    ok <- is.finite(model_vals) & is.finite(real_vals)
    if (any(ok)) {
      real_ok <- all(real_vals[ok] >= model_vals[ok])
    }
  }

  rel_ok <- TRUE
  if (all(c("Reliability", "RealReliability") %in% names(reliability_tbl))) {
    rel_vals <- suppressWarnings(as.numeric(reliability_tbl$Reliability))
    real_rel_vals <- suppressWarnings(as.numeric(reliability_tbl$RealReliability))
    ok <- is.finite(rel_vals) & is.finite(real_rel_vals)
    if (any(ok)) {
      rel_ok <- all(real_rel_vals[ok] <= rel_vals[ok])
    }
  }

  coverage_ok <- FALSE
  if (nrow(precision_profile_tbl) > 0 && "HasSamplePopulationCoverage" %in% names(precision_profile_tbl)) {
    coverage_ok <- isTRUE(precision_profile_tbl$HasSamplePopulationCoverage[1])
  }

  label_status <- "review"
  label_detail <- "SE source labels were not available for audit."
  if (nrow(measure_df) > 0 && all(c("Facet", "SE_Method") %in% names(measure_df))) {
    labels <- as.character(measure_df$SE_Method)
    person_labels <- labels[measure_df$Facet == "Person"]
    nonperson_labels <- labels[measure_df$Facet != "Person"]
    if (identical(method, "MML")) {
      has_person_eap <- length(person_labels) == 0 || any(grepl("Posterior SD \\(EAP\\)", person_labels))
      has_nonperson_info <- length(nonperson_labels) == 0 || any(grepl("Observed information \\(MML\\)", nonperson_labels))
      has_fallback <- any(grepl("Fallback observation-table information", labels))
      if (has_person_eap && has_nonperson_info && !has_fallback) {
        label_status <- "pass"
        label_detail <- "Person and non-person SE labels match the MML precision path."
      } else if (has_person_eap && (has_nonperson_info || has_fallback)) {
        label_status <- "review"
        label_detail <- "MML labels were found, but at least one level fell back to observation-table information."
      } else {
        label_status <- "warn"
        label_detail <- "MML SE labels did not consistently identify posterior-SD/person and observed-information/facet sources."
      }
    } else {
      uses_obs_info <- all(grepl("Observation-table information", labels))
      if (uses_obs_info) {
        label_status <- "pass"
        label_detail <- "JML SE labels consistently identify observation-table information."
      } else {
        label_status <- "warn"
        label_detail <- "JML SE labels were expected to indicate observation-table information throughout."
      }
    }
  }

  tibble(
    Check = c(
      "Precision tier",
      "Optimizer convergence",
      "ModelSE availability",
      "Fit-adjusted SE ordering",
      "Reliability ordering",
      "Facet precision coverage",
      "SE source labels"
    ),
    Status = c(
      if (nrow(precision_profile_tbl) > 0) {
        if (identical(as.character(precision_profile_tbl$PrecisionTier[1]), "model_based")) "pass" else "review"
      } else if (identical(method, "MML")) {
        "review"
      } else {
        "review"
      },
      if (isTRUE(converged)) "pass" else "review",
      if (!is.finite(model_share)) "warn" else if (model_share >= 0.99) "pass" else if (model_share >= 0.90) "review" else "warn",
      if (isTRUE(real_ok)) "pass" else "warn",
      if (isTRUE(rel_ok)) "pass" else "warn",
      if (isTRUE(coverage_ok)) "pass" else "review",
      label_status
    ),
    Detail = c(
      if (nrow(precision_profile_tbl) > 0 && identical(as.character(precision_profile_tbl$PrecisionTier[1]), "model_based")) {
        if (isTRUE(converged)) {
          "This run uses the package's model-based precision path."
        } else {
          "This run uses the package's model-based precision path, but optimizer convergence should be reviewed before formal reporting."
        }
      } else if (nrow(precision_profile_tbl) > 0 && identical(as.character(precision_profile_tbl$PrecisionTier[1]), "hybrid")) {
        "This run mixes model-based SE with fallback observation-table SE for at least one level; review before treating the full run as formal inference."
      } else {
        "This run uses the package's exploratory precision path; prefer MML for formal SE, CI, and reliability reporting."
      },
      if (isTRUE(converged)) {
        "The optimizer reported convergence."
      } else {
        "The optimizer did not report convergence; keep SE, CI, and reliability in review mode."
      },
      if (!is.finite(model_share)) {
        "ModelSE coverage could not be computed."
      } else {
        paste0("Finite ModelSE values were available for ", sprintf("%.1f", 100 * model_share), "% of rows.")
      },
      if (isTRUE(real_ok)) {
        "Fit-adjusted SE values were not smaller than their paired ModelSE values."
      } else {
        "At least one RealSE value was smaller than ModelSE."
      },
      if (isTRUE(rel_ok)) {
        "Conservative reliability values were not larger than the model-based values."
      } else {
        "At least one RealReliability value exceeded Reliability."
      },
      if (isTRUE(coverage_ok)) {
        "Each facet had sample/population summaries for both model and fit-adjusted SE modes."
      } else {
        "At least one facet is missing a sample/population or model/fit-adjusted precision combination."
      },
      label_detail
    )
  )
}

calc_reliability <- function(measure_df) {
  measure_df |>
    dplyr::group_by(Facet) |>
    dplyr::group_modify(function(.x, .y) {
      model_stats <- summarize_precision_basis(.x, if ("ModelSE" %in% names(.x)) "ModelSE" else "SE", distribution_basis = "sample")
      real_stats <- summarize_precision_basis(.x, if ("RealSE" %in% names(.x)) "RealSE" else if ("ModelSE" %in% names(.x)) "ModelSE" else "SE", distribution_basis = "sample")
      tier_vals <- unique(as.character(stats::na.omit(.x$PrecisionTier %||% character(0))))
      converged_vals <- unique(as.logical(stats::na.omit(.x$Converged %||% logical(0))))
      facet_converged <- if (length(converged_vals) == 0) {
        FALSE
      } else {
        all(converged_vals)
      }
      facet_tier <- if (length(tier_vals) == 0) {
        NA_character_
      } else if (all(tier_vals == "model_based")) {
        "model_based"
      } else if (any(tier_vals == "exploratory")) {
        "exploratory"
      } else {
        "hybrid"
      }
      supports_formal <- identical(facet_tier, "model_based") && isTRUE(facet_converged)

      tibble(
        Levels = nrow(.x),
        Converged = facet_converged,
        PrecisionTier = facet_tier,
        SupportsFormalInference = supports_formal,
        ReliabilityUse = dplyr::case_when(
          identical(facet_tier, "model_based") && isTRUE(facet_converged) ~ "primary_reporting",
          identical(facet_tier, "model_based") ~ "review_before_reporting",
          identical(facet_tier, "hybrid") ~ "review_before_reporting",
          identical(facet_tier, "exploratory") ~ "screening_only",
          TRUE ~ NA_character_
        ),
        SD = model_stats$ObservedSD,
        RMSE = model_stats$RMSE,
        TrueSD = model_stats$TrueSD,
        ObservedVariance = model_stats$ObservedVariance,
        ModelErrorVariance = model_stats$ErrorVariance,
        ModelTrueVariance = model_stats$TrueVariance,
        Separation = model_stats$Separation,
        Strata = model_stats$Strata,
        Reliability = model_stats$Reliability,
        ModelRMSE = model_stats$RMSE,
        ModelTrueSD = model_stats$TrueSD,
        ModelSeparation = model_stats$Separation,
        ModelStrata = model_stats$Strata,
        ModelReliability = model_stats$Reliability,
        RealRMSE = real_stats$RMSE,
        RealTrueSD = real_stats$TrueSD,
        RealErrorVariance = real_stats$ErrorVariance,
        RealTrueVariance = real_stats$TrueVariance,
        RealSeparation = real_stats$Separation,
        RealStrata = real_stats$Strata,
        RealReliability = real_stats$Reliability,
        MeanInfit = if ("Infit" %in% names(.x)) mean(.x$Infit, na.rm = TRUE) else NA_real_,
        MeanOutfit = if ("Outfit" %in% names(.x)) mean(.x$Outfit, na.rm = TRUE) else NA_real_
      )
    }) |>
    dplyr::ungroup()
}

ensure_positive_definite <- function(mat) {
  eig <- tryCatch(eigen(mat, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
  if (is.null(eig)) return(mat)
  if (any(eig < .Machine$double.eps)) {
    smoothed <- tryCatch(suppressWarnings(psych::cor.smooth(mat)), error = function(e) NULL)
    if (!is.null(smoothed)) {
      if (is.list(smoothed) && !is.null(smoothed$R)) {
        mat <- smoothed$R
      } else {
        mat <- smoothed
      }
    }
  }
  mat
}

compute_principal_components <- function(cor_matrix, n_factors) {
  warnings <- character(0)
  pca <- withCallingHandlers(
    tryCatch(
      psych::principal(cor_matrix, nfactors = n_factors, rotate = "none"),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(
    pca = pca,
    warning = if (length(warnings) == 0) NULL else paste(unique(warnings), collapse = "; ")
  )
}

compute_pca_overall <- function(obs_df, facet_names, max_factors = 10L) {
  pca_failure <- function(message, residual_matrix = NULL, cor_matrix = NULL, warning = NULL) {
    list(pca = NULL, residual_matrix = residual_matrix, cor_matrix = cor_matrix, error = message, warning = warning)
  }
  if (length(facet_names) == 0) return(NULL)
  df_aug <- obs_df |>
    mutate(
      Person = as.character(Person),
      item_combination = paste(!!!rlang::syms(facet_names), sep = "_")
    ) |>
    select(Person, item_combination, StdResidual)

  residual_matrix_prep <- df_aug |>
    group_by(Person, item_combination) |>
    summarize(std_residual = mean(StdResidual, na.rm = TRUE), .groups = "drop")

  residual_matrix_wide <- tryCatch({
    residual_matrix_prep |>
      tidyr::pivot_wider(
        id_cols = Person,
        names_from = item_combination,
        values_from = std_residual,
        values_fill = list(std_residual = NA)
      ) |>
      tibble::column_to_rownames("Person")
  }, error = function(e) e)

  if (inherits(residual_matrix_wide, "error")) {
    return(pca_failure(paste0("Could not reshape the residual matrix: ", conditionMessage(residual_matrix_wide))))
  }
  if (is.null(residual_matrix_wide) || nrow(residual_matrix_wide) < 2 || ncol(residual_matrix_wide) < 2) {
    return(pca_failure("Residual PCA requires at least two persons and two combined-facet columns."))
  }

  residual_matrix_clean <- residual_matrix_wide[, colSums(is.na(residual_matrix_wide)) < nrow(residual_matrix_wide), drop = FALSE]
  if (ncol(residual_matrix_clean) < 2) {
    return(pca_failure("Residual PCA requires at least two combined-facet columns with observed residuals.", residual_matrix = residual_matrix_wide))
  }

  cor_matrix <- tryCatch(suppressWarnings(stats::cor(residual_matrix_clean, use = "pairwise.complete.obs")), error = function(e) e)
  if (inherits(cor_matrix, "error")) {
    return(pca_failure(
      paste0("Could not compute the residual correlation matrix: ", conditionMessage(cor_matrix)),
      residual_matrix = residual_matrix_wide
    ))
  }
  if (is.null(cor_matrix)) {
    return(pca_failure("Residual correlation matrix was not available.", residual_matrix = residual_matrix_wide))
  }
  cor_matrix[is.na(cor_matrix)] <- 0
  diag(cor_matrix) <- 1
  cor_matrix <- ensure_positive_definite(cor_matrix)

  max_cap <- if (is.na(max_factors) || !is.finite(max_factors)) {
    min(10L, ncol(cor_matrix) - 1L, nrow(cor_matrix) - 1L)
  } else {
    as.integer(max_factors)
  }
  n_factors <- max(1, min(max_cap, ncol(cor_matrix) - 1, nrow(cor_matrix) - 1))
  pca_run <- compute_principal_components(cor_matrix, n_factors)
  pca_result <- pca_run$pca
  if (inherits(pca_result, "error")) {
    return(pca_failure(
      paste0("Principal components could not be computed: ", conditionMessage(pca_result)),
      residual_matrix = residual_matrix_wide,
      cor_matrix = cor_matrix,
      warning = pca_run$warning
    ))
  }
  list(pca = pca_result, residual_matrix = residual_matrix_wide, cor_matrix = cor_matrix, error = NULL, warning = pca_run$warning)
}

compute_pca_by_facet <- function(obs_df, facet_names, max_factors = 10L) {
  out <- list()
  for (facet in facet_names) {
    pca_failure <- function(message, residual_matrix = NULL, cor_matrix = NULL, warning = NULL) {
      list(pca = NULL, residual_matrix = residual_matrix, cor_matrix = cor_matrix, error = message, warning = warning)
    }
    facet_sym <- rlang::sym(facet)
    prep <- obs_df |>
      mutate(Person = as.character(Person), .Level = as.character(!!facet_sym)) |>
      select(Person, .Level, StdResidual) |>
      group_by(Person, .Level) |>
      summarize(std_residual = mean(StdResidual, na.rm = TRUE), .groups = "drop")

    wide <- tryCatch({
      tidyr::pivot_wider(
        prep,
        id_cols = Person,
        names_from = .Level,
        values_from = std_residual,
        values_fill = list(std_residual = NA)
      ) |>
        tibble::column_to_rownames("Person")
    }, error = function(e) e)

    if (inherits(wide, "error")) {
      out[[facet]] <- pca_failure(paste0("Could not reshape the residual matrix: ", conditionMessage(wide)))
      next
    }
    if (is.null(wide) || nrow(wide) < 2 || ncol(wide) < 2) {
      out[[facet]] <- pca_failure("Residual PCA requires at least two persons and two facet-level columns.")
      next
    }

    keep <- colSums(is.na(wide)) < nrow(wide)
    wide <- wide[, keep, drop = FALSE]
    if (ncol(wide) < 2) {
      out[[facet]] <- pca_failure("Residual PCA requires at least two facet-level columns with observed residuals.", residual_matrix = wide)
      next
    }

    cor_mat <- tryCatch(suppressWarnings(stats::cor(wide, use = "pairwise.complete.obs")), error = function(e) e)
    if (inherits(cor_mat, "error")) {
      out[[facet]] <- pca_failure(
        paste0("Could not compute the residual correlation matrix: ", conditionMessage(cor_mat)),
        residual_matrix = wide
      )
      next
    }
    if (is.null(cor_mat)) {
      out[[facet]] <- pca_failure("Residual correlation matrix was not available.", residual_matrix = wide)
      next
    }
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1
    cor_mat <- ensure_positive_definite(cor_mat)

    max_cap <- if (is.na(max_factors) || !is.finite(max_factors)) {
      min(10L, ncol(cor_mat) - 1L, nrow(cor_mat) - 1L)
    } else {
      as.integer(max_factors)
    }
    n_factors <- max(1, min(max_cap, ncol(cor_mat) - 1, nrow(cor_mat) - 1))
    pca_run <- compute_principal_components(cor_mat, n_factors)
    pca_obj <- pca_run$pca
    if (inherits(pca_obj, "error")) {
      out[[facet]] <- pca_failure(
        paste0("Principal components could not be computed: ", conditionMessage(pca_obj)),
        residual_matrix = wide,
        cor_matrix = cor_mat,
        warning = pca_run$warning
      )
      next
    }
    out[[facet]] <- list(pca = pca_obj, cor_matrix = cor_mat, residual_matrix = wide, error = NULL, warning = pca_run$warning)
  }
  out
}

mfrm_diagnostics <- function(res,
                             interaction_pairs = NULL,
                             top_n_interactions = 20,
                             whexact = FALSE,
                             fit_df_method = c("engine", "facets", "both"),
                             diagnostic_mode = c("legacy", "marginal_fit", "both"),
                             residual_pca = c("none", "overall", "facet", "both"),
                             pca_max_factors = 10L) {
  fit_df_method <- match_fit_df_method(fit_df_method)
  diagnostic_mode <- match.arg(diagnostic_mode, c("legacy", "marginal_fit", "both"))
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))
  obs_df <- compute_obs_table(res)
  facet_cols <- c("Person", res$config$facet_names)
  overall_fit <- calc_overall_fit(
    obs_df,
    whexact = whexact,
    fit_df_method = fit_df_method
  )
  fit_tbl <- calc_facet_fit(
    obs_df,
    facet_cols,
    whexact = whexact,
    fit_df_method = fit_df_method
  )
  parameter_covariance <- compute_mml_parameter_covariance(res)
  se_tbl <- build_measure_se_table(res, obs_df, facet_cols, fit_tbl,
                                   covariance = parameter_covariance)
  structural_uncertainty <- compute_mml_structural_parameter_se(
    res,
    covariance = parameter_covariance
  )
  bias_tbl <- calc_bias_facet(obs_df, facet_cols)
  interaction_tbl <- calc_bias_interactions(obs_df, facet_cols, pairs = interaction_pairs,
                                             top_n = top_n_interactions)
  ptmea_tbl <- calc_ptmea(obs_df, facet_cols)
  subset_tbls <- calc_subsets(obs_df, facet_cols)

  person_tbl <- res$facets$person |>
    mutate(
      Facet = "Person",
      Level = Person
    )
  facet_tbl <- res$facets$others |>
    mutate(Level = as.character(Level))

  measure_ci_level <- 0.95
  measure_ci_z <- stats::qnorm(1 - (1 - measure_ci_level) / 2)

  measures <- bind_rows(
    person_tbl |>
      select(Facet, Level, Estimate),
    facet_tbl |>
      select(Facet, Level, Estimate)
  ) |>
    left_join(se_tbl, by = c("Facet", "Level")) |>
    left_join(fit_tbl, by = c("Facet", "Level")) |>
    left_join(bias_tbl, by = c("Facet", "Level")) |>
    left_join(ptmea_tbl, by = c("Facet", "Level")) |>
    mutate(
      CI_Lower = ifelse(is.finite(SE), Estimate - measure_ci_z * SE, NA_real_),
      CI_Upper = ifelse(is.finite(SE), Estimate + measure_ci_z * SE, NA_real_),
      CI_Level = measure_ci_level,
      CI_Method = "Normal approximation",
      CIEligible = dplyr::coalesce(.data$SupportsFormalInference, FALSE),
      CILabel = dplyr::case_when(
        .data$PrecisionTier == "model_based" ~ "Model-based normal interval",
        .data$PrecisionTier == "hybrid" ~ "Approximate interval; review fallback SE",
        .data$PrecisionTier == "exploratory" ~ "Approximate interval; screening only",
        TRUE ~ "Approximate interval"
      )
    )

  if (!"N" %in% names(measures)) {
    if ("N.x" %in% names(measures) || "N.y" %in% names(measures)) {
      measures <- measures |>
        mutate(N = dplyr::coalesce(.data$N.x, .data$N.y))
    } else {
      measures <- measures |>
        mutate(N = NA_real_)
    }
  }

  reliability_tbl <- calc_reliability(measures)
  facets_chisq_tbl <- calc_facets_chisq(measures)
  facet_precision_tbl <- build_facet_precision_summary(measures, facets_chisq_tbl)
  default_rater_facet <- infer_default_rater_facet(res$config$facet_names)
  interrater_tbl <- augment_interrater_with_precision(
    calc_interrater_agreement(
      obs_df = obs_df,
      facet_cols = facet_cols,
      rater_facet = default_rater_facet,
      res = res
    ),
    reliability_tbl = reliability_tbl,
    rater_facet = default_rater_facet
  )

  unexpected_abs_z_min <- 2
  unexpected_prob_max <- 0.30
  unexpected_rule <- "either"
  unexpected_tbl <- calc_unexpected_response_table(
    obs_df = obs_df,
    probs = compute_prob_matrix(res),
    facet_names = res$config$facet_names,
    rating_min = res$prep$rating_min,
    abs_z_min = unexpected_abs_z_min,
    prob_max = unexpected_prob_max,
    top_n = 100,
    rule = unexpected_rule
  )
  unexpected_summary <- summarize_unexpected_response_table(
    unexpected_tbl = unexpected_tbl,
    total_observations = nrow(obs_df),
    abs_z_min = unexpected_abs_z_min,
    prob_max = unexpected_prob_max,
    rule = unexpected_rule
  )

  fair_average <- if (identical(res$config$model, "GPCM")) {
    list(
      raw_by_facet = list(),
      by_facet = list(),
      stacked = tibble(),
      available = FALSE,
      reason = gpcm_fair_average_rationale()
    )
  } else {
    calc_fair_average_bundle(
      res = res,
      diagnostics = list(obs = obs_df, measures = measures),
      totalscore = TRUE,
      umean = 0,
      uscale = 1,
      udecimals = 2,
      omit_unobserved = FALSE,
      xtreme = 0
    )
  }

  displacement_abs_warn <- 0.5
  displacement_abs_t_warn <- 2
  displacement_tbl <- calc_displacement_table(
    obs_df = obs_df,
    res = res,
    measures = measures,
    abs_displacement_warn = displacement_abs_warn,
    abs_t_warn = displacement_abs_t_warn
  )
  displacement_summary <- summarize_displacement_table(
    displacement_tbl = displacement_tbl,
    abs_displacement_warn = displacement_abs_warn,
    abs_t_warn = displacement_abs_t_warn
  )

  pca_overall <- NULL
  pca_by_facet <- NULL
  if (residual_pca %in% c("overall", "both")) {
    pca_overall <- compute_pca_overall(
      obs_df = obs_df,
      facet_names = res$config$facet_names,
      max_factors = pca_max_factors
    )
  }
  if (residual_pca %in% c("facet", "both")) {
    pca_by_facet <- compute_pca_by_facet(
      obs_df = obs_df,
      facet_names = res$config$facet_names,
      max_factors = pca_max_factors
    )
  }

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  se_detail <- attr(se_tbl, "mml_se_detail", exact = TRUE) %||%
    "SE values default to observation-table information approximations."
  precision_profile_tbl <- build_precision_profile(res, measures, reliability_tbl, facet_precision_tbl)
  precision_review_tbl <- audit_precision_outputs(
    res,
    measures,
    reliability_tbl,
    facet_precision_tbl,
    precision_profile_tbl = precision_profile_tbl
  )
  approximation_notes <- tibble::tibble(
    Component = c("SE", "RealSE", "CI", "Reliability"),
    Method = c(
      if (identical(method, "MML")) "Model-based SE (MML)" else "Model-based SE (exploratory JML)",
      "Fit-adjusted SE",
      "Normal interval",
      "Variance decomposition"
    ),
    Detail = c(
      paste0(
        "The `SE` column is kept as a compatibility alias for `ModelSE`. ",
        if (identical(method, "MML")) {
          "Persons use posterior SD from EAP estimation; non-person facets use observed-information SEs from the marginal log-likelihood. "
        } else {
          "JML uses observation-table information as an exploratory approximation. "
        },
        se_detail,
        " Row-level `PrecisionTier`, `Converged`, `SupportsFormalInference`, and `SEUse` indicate whether a given estimate is suitable for primary reporting."
      ),
      "RealSE is a fit-adjusted companion to ModelSE: ModelSE * sqrt(max(Infit, 1)), so it is never smaller than ModelSE.",
      "CI_Lower and CI_Upper are symmetric 95% normal bands computed as Estimate +/- qnorm(0.975) * SE. `CI_Level` records 0.95 and `CI_Method` records the interval family. Use `CIEligible`, `Converged`, `CIBasis`, and `CIUse` to distinguish primary-reporting intervals from review or screening approximations.",
      "Reliability tables report model and real bounds using observed variance, error variance, and true variance (Observed variance - mean SE^2). `Reliability`/`Separation` remain compatibility aliases for the model-based values, while `PrecisionTier`, `Converged`, `SupportsFormalInference`, and `ReliabilityUse` indicate how strongly each facet summary supports formal reporting."
    )
  )
  marginal_fit <- if (diagnostic_mode %in% c("marginal_fit", "both")) {
    calc_marginal_fit_bundle(res)
  } else {
    list(
      available = FALSE,
      summary = tibble(
        Available = FALSE,
        Method = method,
        Model = as.character(res$config$model %||% res$summary$Model[1] %||% NA_character_),
        Basis = "latent_integrated_first_order_counts",
        Reason = "Not computed; use `diagnostic_mode = \"marginal_fit\"` or `\"both\"`."
      ),
      thresholds = list(abs_z_warn = 2, rmsd_warn = 0.05),
      notes = "Strict marginal diagnostics were not requested for this run."
    )
  }
  diagnostic_basis_tbl <- build_diagnostic_basis_guide(
    diagnostic_mode = diagnostic_mode,
    marginal_fit = marginal_fit
  )
  fit_standardization <- tibble::tibble(
    PrimaryFitDfMethod = fit_df_method,
    PrimaryZSTDColumns = if (identical(fit_df_method, "facets")) {
      "InfitZSTD / OutfitZSTD use FACETS-style df"
    } else {
      "InfitZSTD / OutfitZSTD use engine df"
    },
    EngineDf = "DF_Infit = sum(Var * Weight); DF_Outfit = sum(Weight)",
    FacetsDf = "DF = 2 / q^2 using the Wright-Masters/FACETS fourth-moment approximation",
    FacetsColumnsAvailable = identical(fit_df_method, "both") || identical(fit_df_method, "facets"),
    ZSTDTransform = fit_zstd_transform_label(whexact),
    FacetsZSTDCap = if (identical(fit_df_method, "engine")) NA_real_ else 9
  )

  list(
    obs = obs_df,
    facet_names = res$config$facet_names,
    diagnostic_mode = diagnostic_mode,
    diagnostic_basis = diagnostic_basis_tbl,
    fit_standardization = fit_standardization,
    overall_fit = overall_fit,
    measures = measures,
    fit = fit_tbl,
    reliability = reliability_tbl,
    precision_profile = precision_profile_tbl,
    precision_review = precision_review_tbl,
    parameter_uncertainty = structural_uncertainty,
    facet_precision = facet_precision_tbl,
    facets_chisq = facets_chisq_tbl,
    bias = bias_tbl,
    interactions = interaction_tbl,
    interrater = interrater_tbl,
    unexpected = list(
      table = unexpected_tbl,
      summary = unexpected_summary,
      thresholds = list(
        abs_z_min = unexpected_abs_z_min,
        prob_max = unexpected_prob_max,
        rule = unexpected_rule
      )
    ),
    fair_average = fair_average,
    displacement = list(
      table = displacement_tbl,
      summary = displacement_summary,
      thresholds = list(
        abs_displacement_warn = displacement_abs_warn,
        abs_t_warn = displacement_abs_t_warn
      )
    ),
    approximation_notes = approximation_notes,
    marginal_fit = marginal_fit,
    subsets = subset_tbls,
    residual_pca_mode = residual_pca,
    residual_pca_overall = pca_overall,
    residual_pca_by_facet = pca_by_facet
  )
}
