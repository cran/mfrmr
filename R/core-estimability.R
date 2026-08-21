# ==============================================================================
# Estimator-specific constrained estimability audit
# ==============================================================================
#
# The audit works in the exact free-coordinate parameterization used by the
# optimizer. For RSM/PCM it builds the sparse adjacent-category-logit design,
# including Person coordinates for JML and excluding them for MML. Anchors,
# group constraints, facet signs, two-way interactions, and step constraints
# enter through their actual expansion Jacobians.
#
# GPCM log slopes and latent-regression residual variance are nonlinear
# coordinates. The additive block is audited before fitting. Retained-point
# conditional or marginal score maps can supply a sufficient local first-order
# full-rank certificate, but neither that certificate nor a bounded fitted-
# information instrument establishes global identification or weak-information
# readiness. This distinction is retained in `Complete` and `NonlinearBlocks`.

mfrmr_estimability_contract_version <- function() {
  mfrmr_readiness_contract_version()
}

mfrmr_empty_sparse_matrix <- function(nrow = 0L, ncol = 0L) {
  Matrix::sparseMatrix(
    i = integer(0), j = integer(0), x = numeric(0),
    dims = c(as.integer(nrow), as.integer(ncol))
  )
}

mfrmr_estimability_map <- function(Block = character(0),
                                    Coordinate = character(0),
                                    Facet = character(0),
                                    Level = character(0),
                                    ReferenceLevel = character(0),
                                    Constraint = character(0),
                                    OptimizerIndex = integer(0)) {
  data.frame(
    Block = as.character(Block),
    Coordinate = as.character(Coordinate),
    Facet = as.character(Facet),
    Level = as.character(Level),
    ReferenceLevel = as.character(ReferenceLevel),
    Constraint = as.character(Constraint),
    OptimizerIndex = as.integer(OptimizerIndex),
    stringsAsFactors = FALSE
  )
}

mfrmr_constraint_jacobian_sparse <- function(spec, block) {
  levels <- as.character(spec$levels %||% character(0))
  n_levels <- length(levels)
  anchors <- as.numeric(spec$anchors %||% rep(NA_real_, n_levels))
  groups <- as.character(spec$groups %||% rep(NA_character_, n_levels))
  if (length(anchors) != n_levels || length(groups) != n_levels) {
    stop("Internal estimability audit found a malformed facet constraint.",
         call. = FALSE)
  }

  free_idx <- which(is.na(anchors))
  row_i <- integer(0)
  col_j <- integer(0)
  values <- numeric(0)
  map_rows <- list()
  col_cursor <- 0L

  add_column <- function(level_idx, reference_idx, constraint, group = "") {
    col_cursor <<- col_cursor + 1L
    row_i <<- c(row_i, level_idx, reference_idx)
    col_j <<- c(col_j, col_cursor, col_cursor)
    values <<- c(values, 1, -1)
    suffix <- if (nzchar(group)) paste0("|group=", group) else ""
    map_rows[[col_cursor]] <<- mfrmr_estimability_map(
      Block = block,
      Coordinate = paste0(block, ":", levels[level_idx], suffix),
      Facet = block,
      Level = levels[level_idx],
      ReferenceLevel = levels[reference_idx],
      Constraint = constraint,
      OptimizerIndex = NA_integer_
    )
  }

  group_ids <- unique(stats::na.omit(groups[free_idx]))
  group_ids <- group_ids[nzchar(group_ids)]
  grouped_idx <- integer(0)
  if (length(group_ids) > 0L) {
    for (group_id in group_ids) {
      group_levels <- which(groups == group_id)
      free_in_group <- group_levels[is.na(anchors[group_levels])]
      grouped_idx <- c(grouped_idx, free_in_group)
      if (length(free_in_group) <= 1L) next
      reference <- free_in_group[length(free_in_group)]
      for (level_idx in free_in_group[-length(free_in_group)]) {
        add_column(
          level_idx = level_idx,
          reference_idx = reference,
          constraint = "group_sum",
          group = group_id
        )
      }
    }
  }

  ungrouped_idx <- setdiff(free_idx, grouped_idx)
  if (length(ungrouped_idx) > 0L) {
    if (isTRUE(spec$centered)) {
      if (length(ungrouped_idx) > 1L) {
        reference <- ungrouped_idx[length(ungrouped_idx)]
        for (level_idx in ungrouped_idx[-length(ungrouped_idx)]) {
          add_column(
            level_idx = level_idx,
            reference_idx = reference,
            constraint = "sum_zero"
          )
        }
      }
    } else {
      for (level_idx in ungrouped_idx) {
        col_cursor <- col_cursor + 1L
        row_i <- c(row_i, level_idx)
        col_j <- c(col_j, col_cursor)
        values <- c(values, 1)
        map_rows[[col_cursor]] <- mfrmr_estimability_map(
          Block = block,
          Coordinate = paste0(block, ":", levels[level_idx]),
          Facet = block,
          Level = levels[level_idx],
          ReferenceLevel = "",
          Constraint = "free_location",
          OptimizerIndex = NA_integer_
        )
      }
    }
  }

  expected <- as.integer(spec$n_params %||% 0L)
  if (col_cursor != expected) {
    stop(
      "Internal estimability audit/free-parameter mismatch for ",
      shQuote(block), ": audit built ", col_cursor,
      " coordinate(s), optimizer expects ", expected, ".",
      call. = FALSE
    )
  }

  jacobian <- if (col_cursor == 0L) {
    mfrmr_empty_sparse_matrix(n_levels, 0L)
  } else {
    Matrix::sparseMatrix(
      i = row_i, j = col_j, x = values,
      dims = c(n_levels, col_cursor)
    )
  }
  map <- if (length(map_rows) == 0L) {
    mfrmr_estimability_map()
  } else {
    do.call(rbind, map_rows)
  }
  rownames(map) <- NULL
  dimnames(jacobian) <- list(
    levels,
    if (nrow(map) > 0L) map$Coordinate else character(0)
  )
  list(jacobian = jacobian, map = map)
}

mfrmr_interaction_jacobian_sparse <- function(spec) {
  n_a <- as.integer(spec$n_a %||% 0L)
  n_b <- as.integer(spec$n_b %||% 0L)
  n_params <- as.integer(spec$n_params %||% 0L)
  n_cells <- n_a * n_b
  if (n_params == 0L) {
    return(list(
      jacobian = mfrmr_empty_sparse_matrix(n_cells, 0L),
      map = mfrmr_estimability_map()
    ))
  }

  row_i <- integer(0)
  col_j <- integer(0)
  values <- numeric(0)
  map_rows <- vector("list", n_params)
  col_cursor <- 0L
  cell <- function(a, b) a + (b - 1L) * n_a
  for (b in seq_len(n_b - 1L)) {
    for (a in seq_len(n_a - 1L)) {
      col_cursor <- col_cursor + 1L
      row_i <- c(row_i, cell(a, b), cell(a, n_b),
                 cell(n_a, b), cell(n_a, n_b))
      col_j <- c(col_j, rep(col_cursor, 4L))
      values <- c(values, 1, -1, -1, 1)
      map_rows[[col_cursor]] <- mfrmr_estimability_map(
        Block = "interactions",
        Coordinate = paste0(
          spec$name, ":", spec$levels_a[a], "*", spec$levels_b[b]
        ),
        Facet = spec$name,
        Level = paste(spec$levels_a[a], spec$levels_b[b], sep = "*"),
        ReferenceLevel = paste(spec$levels_a[n_a], spec$levels_b[n_b],
                               sep = "*"),
        Constraint = "two_way_sum_zero_margins",
        OptimizerIndex = NA_integer_
      )
    }
  }
  if (col_cursor != n_params) {
    stop("Internal interaction estimability dimension mismatch.", call. = FALSE)
  }
  list(
    jacobian = Matrix::sparseMatrix(
      i = row_i, j = col_j, x = values,
      dims = c(n_cells, n_params)
    ),
    map = do.call(rbind, map_rows)
  )
}

mfrmr_step_jacobian_sparse <- function(config, sizes) {
  n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
  n_free_per_scope <- sum_zero_param_count(n_steps)
  n_scope <- if (identical(config$model, "RSM")) {
    1L
  } else {
    length(config$facet_levels[[config$step_facet]] %||% character(0))
  }
  n_rows <- n_scope * n_steps
  n_cols <- n_scope * n_free_per_scope
  if (n_cols == 0L) {
    return(list(
      jacobian = mfrmr_empty_sparse_matrix(n_rows, 0L),
      map = mfrmr_estimability_map()
    ))
  }

  row_i <- integer(0)
  col_j <- integer(0)
  values <- numeric(0)
  map_rows <- vector("list", n_cols)
  col_cursor <- 0L
  scope_levels <- if (identical(config$model, "RSM")) {
    "shared"
  } else {
    as.character(config$facet_levels[[config$step_facet]])
  }
  for (scope in seq_len(n_scope)) {
    row_offset <- (scope - 1L) * n_steps
    for (transition in seq_len(n_free_per_scope)) {
      col_cursor <- col_cursor + 1L
      row_i <- c(row_i, row_offset + transition, row_offset + n_steps)
      col_j <- c(col_j, col_cursor, col_cursor)
      values <- c(values, 1, -1)
      map_rows[[col_cursor]] <- mfrmr_estimability_map(
        Block = "steps",
        Coordinate = paste0(
          "steps:", scope_levels[scope], ":transition", transition
        ),
        Facet = if (identical(config$model, "RSM")) {
          "shared_rating_scale"
        } else {
          config$step_facet
        },
        Level = scope_levels[scope],
        ReferenceLevel = paste0("transition", n_steps),
        Constraint = "within_ladder_sum_zero",
        OptimizerIndex = NA_integer_
      )
    }
  }
  list(
    jacobian = Matrix::sparseMatrix(
      i = row_i, j = col_j, x = values,
      dims = c(n_rows, n_cols)
    ),
    map = do.call(rbind, map_rows)
  )
}

mfrmr_estimability_set_optimizer_index <- function(map, block, slices) {
  if (nrow(map) == 0L) return(map)
  idx <- slices[[block]] %||% integer(0)
  if (length(idx) == nrow(map)) {
    map$OptimizerIndex <- as.integer(idx)
  }
  map
}

mfrmr_estimability_eta_design <- function(prep, idx, config, sizes,
                                           include_person,
                                           include_population_beta) {
  n_obs <- length(idx$score_k)
  slices <- build_param_slices(sizes)
  blocks <- list()
  maps <- list()

  append_block <- function(name, matrix, map) {
    if (ncol(matrix) == 0L) return(invisible(NULL))
    blocks[[length(blocks) + 1L]] <<- matrix
    maps[[length(maps) + 1L]] <<- map
    invisible(NULL)
  }

  if (isTRUE(include_person)) {
    theta <- mfrmr_constraint_jacobian_sparse(config$theta_spec, "Person")
    theta$map <- mfrmr_estimability_set_optimizer_index(
      theta$map, "theta", slices
    )
    append_block("theta", theta$jacobian[idx$person, , drop = FALSE],
                 theta$map)
  }

  for (facet in config$facet_names) {
    facet_block <- mfrmr_constraint_jacobian_sparse(
      config$facet_specs[[facet]], facet
    )
    sign <- as.numeric(config$facet_signs[[facet]] %||% -1)
    facet_block$map <- mfrmr_estimability_set_optimizer_index(
      facet_block$map, facet, slices
    )
    append_block(
      facet,
      sign * facet_block$jacobian[idx$facets[[facet]], , drop = FALSE],
      facet_block$map
    )
  }

  interaction_specs <- config$interaction_specs %||% list()
  if (length(interaction_specs) > 0L) {
    interaction_matrices <- list()
    interaction_maps <- list()
    for (name in names(interaction_specs)) {
      interaction <- mfrmr_interaction_jacobian_sparse(
        interaction_specs[[name]]
      )
      interaction_matrices[[length(interaction_matrices) + 1L]] <-
        interaction$jacobian[idx$interactions[[name]], , drop = FALSE]
      interaction_maps[[length(interaction_maps) + 1L]] <- interaction$map
    }
    interaction_matrix <- do.call(cbind, interaction_matrices)
    interaction_map <- do.call(rbind, interaction_maps)
    rownames(interaction_map) <- NULL
    interaction_map <- mfrmr_estimability_set_optimizer_index(
      interaction_map, "interactions", slices
    )
    append_block("interactions", interaction_matrix, interaction_map)
  }

  if (isTRUE(include_population_beta) &&
      isTRUE(config$population_spec$active) &&
      as.integer(sizes$beta %||% 0L) > 0L) {
    pop <- config$population_spec
    lookup <- as.integer(pop$person_lookup %||% integer(0))
    if (length(lookup) != config$n_person || anyNA(lookup[idx$person])) {
      stop("Internal estimability audit could not align population-model persons.",
           call. = FALSE)
    }
    beta_matrix <- as.matrix(pop$design_matrix)[lookup[idx$person], , drop = FALSE]
    beta_names <- colnames(beta_matrix)
    if (is.null(beta_names)) beta_names <- paste0("beta_", seq_len(ncol(beta_matrix)))
    beta_map <- mfrmr_estimability_map(
      Block = rep("beta", ncol(beta_matrix)),
      Coordinate = paste0("beta:", beta_names),
      Facet = rep("PersonPopulation", ncol(beta_matrix)),
      Level = beta_names,
      ReferenceLevel = rep("", ncol(beta_matrix)),
      Constraint = rep("population_location", ncol(beta_matrix)),
      OptimizerIndex = as.integer(slices$beta %||% rep(NA_integer_, ncol(beta_matrix)))
    )
    append_block(
      "beta", Matrix::Matrix(beta_matrix, sparse = TRUE), beta_map
    )
  }

  design <- if (length(blocks) == 0L) {
    mfrmr_empty_sparse_matrix(n_obs, 0L)
  } else {
    do.call(cbind, blocks)
  }
  map <- if (length(maps) == 0L) {
    mfrmr_estimability_map()
  } else {
    do.call(rbind, maps)
  }
  rownames(map) <- NULL
  list(design = design, map = map)
}

mfrmr_estimability_adjacent_design <- function(prep, idx, config, sizes,
                                                include_person = NULL,
                                                include_population_beta = NULL) {
  if (is.null(include_person)) {
    include_person <- identical(config$method, "JML")
  }
  if (is.null(include_population_beta)) {
    include_population_beta <- identical(config$method, "MML")
  }
  eta <- mfrmr_estimability_eta_design(
    prep = prep,
    idx = idx,
    config = config,
    sizes = sizes,
    include_person = include_person,
    include_population_beta = include_population_beta
  )
  n_obs <- nrow(eta$design)
  n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
  if (n_steps <= 0L) {
    stop("Estimability audit requires at least two response categories.",
         call. = FALSE)
  }
  transition <- rep(seq_len(n_steps), each = n_obs)
  obs_index <- rep(seq_len(n_obs), times = n_steps)
  eta_repeated <- eta$design[obs_index, , drop = FALSE]

  step <- mfrmr_step_jacobian_sparse(config, sizes)
  step$map <- mfrmr_estimability_set_optimizer_index(
    step$map, "steps", build_param_slices(sizes)
  )
  expanded_step_index <- if (identical(config$model, "RSM")) {
    transition
  } else {
    scope <- rep(as.integer(idx$step_idx), times = n_steps)
    (scope - 1L) * n_steps + transition
  }
  step_selected <- step$jacobian[expanded_step_index, , drop = FALSE]
  design <- cbind(eta_repeated, -step_selected)
  map <- rbind(eta$map, step$map)
  rownames(map) <- NULL
  if (ncol(design) != nrow(map)) {
    stop("Internal estimability design/map dimension mismatch.", call. = FALSE)
  }
  list(
    design = methods::as(design, "dgCMatrix"),
    map = map,
    transition_rows = nrow(design),
    observation_rows = n_obs,
    transitions = n_steps
  )
}

mfrmr_estimability_null_directions <- function(scaled_design,
                                                parameter_map,
                                                rank,
                                                dense_element_limit = 2e6,
                                                max_directions = 5L,
                                                max_coordinates = 8L) {
  n_elements <- as.double(nrow(scaled_design)) * as.double(ncol(scaled_design))
  nullity <- ncol(scaled_design) - rank
  if (nullity <= 0L || n_elements > dense_element_limit ||
      ncol(scaled_design) == 0L) {
    return(data.frame())
  }
  dense <- as.matrix(scaled_design)
  sv <- tryCatch(
    base::svd(dense, nu = 0L, nv = ncol(dense)),
    error = function(e) NULL
  )
  if (is.null(sv) || is.null(sv$v) || ncol(sv$v) < ncol(dense)) {
    return(data.frame())
  }
  direction_idx <- seq.int(rank + 1L, ncol(dense))
  direction_idx <- utils::head(direction_idx, max_directions)
  rows <- lapply(seq_along(direction_idx), function(direction_number) {
    loadings <- sv$v[, direction_idx[direction_number]]
    keep <- order(abs(loadings), decreasing = TRUE)
    keep <- utils::head(keep, max_coordinates)
    data.frame(
      Direction = direction_number,
      Coordinate = parameter_map$Coordinate[keep],
      Block = parameter_map$Block[keep],
      Loading = as.numeric(loadings[keep]),
      AbsLoading = abs(as.numeric(loadings[keep])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mfrmr_estimability_rank_audit <- function(design, parameter_map,
                                           tolerances = c(1e-12, 1e-10, 1e-8),
                                           structural_tolerance = 1e-12) {
  design <- methods::as(design, "dgCMatrix")
  n_params <- ncol(design)
  if (n_params == 0L) {
    return(list(
      Rank = 0L,
      Nullity = 0L,
      State = "identified",
      ToleranceSensitive = FALSE,
      ToleranceRanks = data.frame(
        Tolerance = tolerances,
        Rank = 0L,
        Nullity = 0L
      ),
      ColumnNormMin = NA_real_,
      ColumnNormMax = NA_real_,
      SmallestSingularValue = NA_real_,
      ConditionIndex = NA_real_,
      ZeroCoordinates = character(0),
      NullDirections = data.frame()
    ))
  }

  norm_sq <- as.numeric(Matrix::colSums(design^2))
  norms <- sqrt(pmax(norm_sq, 0))
  zero <- !is.finite(norms) | norms == 0
  inverse_norm <- ifelse(zero, 1, 1 / norms)
  scaled <- design %*% Matrix::Diagonal(x = inverse_norm)
  scaled <- Matrix::drop0(methods::as(scaled, "dgCMatrix"))

  ranks <- vapply(tolerances, function(tolerance) {
    as.integer(Matrix::rankMatrix(
      scaled, method = "qr", tol = as.numeric(tolerance)
    ))
  }, integer(1))
  structural_index <- which.min(abs(tolerances - structural_tolerance))
  structural_rank <- ranks[structural_index]
  nullity <- n_params - structural_rank
  tolerance_sensitive <- length(unique(ranks)) > 1L
  state <- if (nullity > 0L) "structurally_unidentified" else "identified"

  smallest <- NA_real_
  condition_index <- NA_real_
  n_elements <- as.double(nrow(scaled)) * as.double(ncol(scaled))
  if (n_elements <= 2e6 && ncol(scaled) > 0L) {
    singular_values <- tryCatch(
      base::svd(as.matrix(scaled), nu = 0L, nv = 0L)$d,
      error = function(e) numeric(0)
    )
    if (length(singular_values) > 0L) {
      smallest <- min(singular_values)
      positive <- singular_values[singular_values > .Machine$double.eps]
      if (length(positive) > 0L) {
        condition_index <- max(singular_values) / min(positive)
      }
    }
  }

  list(
    Rank = as.integer(structural_rank),
    Nullity = as.integer(nullity),
    State = state,
    ToleranceSensitive = isTRUE(tolerance_sensitive),
    ToleranceRanks = data.frame(
      Tolerance = as.numeric(tolerances),
      Rank = as.integer(ranks),
      Nullity = as.integer(n_params - ranks),
      stringsAsFactors = FALSE
    ),
    ColumnNormMin = if (all(zero)) 0 else min(norms[!zero]),
    ColumnNormMax = if (all(zero)) 0 else max(norms[!zero]),
    SmallestSingularValue = smallest,
    ConditionIndex = condition_index,
    ZeroCoordinates = parameter_map$Coordinate[zero],
    NullDirections = mfrmr_estimability_null_directions(
      scaled_design = scaled,
      parameter_map = parameter_map,
      rank = structural_rank
    )
  )
}

mfrmr_information_rank_ladder <- function(hessian,
                                           tolerances = c(1e-10, 1e-8, 1e-6)) {
  hessian <- suppressWarnings(as.matrix(hessian))
  if (!is.numeric(hessian) || nrow(hessian) != ncol(hessian) ||
      length(hessian) == 0L || any(!is.finite(hessian))) {
    return(list(
      valid = FALSE,
      rank_ladder = data.frame(),
      eigenvalues = numeric(0),
      scale = NA_real_,
      smallest_eigenvalue = NA_real_,
      largest_eigenvalue = NA_real_
    ))
  }

  symmetric <- (hessian + t(hessian)) / 2
  eigenvalues <- tryCatch(
    as.numeric(eigen(symmetric, symmetric = TRUE, only.values = TRUE)$values),
    error = function(e) numeric(0)
  )
  if (length(eigenvalues) != nrow(symmetric) || any(!is.finite(eigenvalues))) {
    return(list(
      valid = FALSE,
      rank_ladder = data.frame(),
      eigenvalues = eigenvalues,
      scale = NA_real_,
      smallest_eigenvalue = NA_real_,
      largest_eigenvalue = NA_real_
    ))
  }

  scale <- max(abs(eigenvalues))
  thresholds <- if (is.finite(scale) && scale > 0) {
    scale * as.numeric(tolerances)
  } else {
    as.numeric(tolerances)
  }
  rank_ladder <- data.frame(
    RelativeTolerance = as.numeric(tolerances),
    AbsoluteThreshold = thresholds,
    PositiveRank = vapply(
      thresholds, function(threshold) sum(eigenvalues > threshold), integer(1)
    ),
    NegativeCount = vapply(
      thresholds, function(threshold) sum(eigenvalues < -threshold), integer(1)
    ),
    NearZeroCount = vapply(
      thresholds,
      function(threshold) sum(abs(eigenvalues) <= threshold),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  list(
    valid = TRUE,
    rank_ladder = rank_ladder,
    eigenvalues = eigenvalues,
    scale = scale,
    smallest_eigenvalue = min(eigenvalues),
    largest_eigenvalue = max(eigenvalues)
  )
}

mfrmr_fitted_information_block_summary <- function(hessian,
                                                    sizes,
                                                    blocks) {
  slices <- build_param_slices(sizes)
  rows <- lapply(as.character(blocks), function(block) {
    slice <- as.integer(slices[[block]] %||% integer(0))
    diagonal <- if (length(slice) > 0L) {
      diag(hessian)[slice]
    } else {
      numeric(0)
    }
    data.frame(
      Block = block,
      FreeCoordinates = length(slice),
      OptimizerIndexStart = if (length(slice) > 0L) min(slice) else NA_integer_,
      OptimizerIndexEnd = if (length(slice) > 0L) max(slice) else NA_integer_,
      FiniteDiagonal = length(diagonal) > 0L && all(is.finite(diagonal)),
      MinimumDiagonal = if (length(diagonal) > 0L) min(diagonal) else NA_real_,
      MaximumDiagonal = if (length(diagonal) > 0L) max(diagonal) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(data.frame(
      Block = character(0),
      FreeCoordinates = integer(0),
      OptimizerIndexStart = integer(0),
      OptimizerIndexEnd = integer(0),
      FiniteDiagonal = logical(0),
      MinimumDiagonal = numeric(0),
      MaximumDiagonal = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

mfrmr_numeric_transformation_jacobian <- function(transform,
                                                   free,
                                                   relative_step = 1e-6) {
  free <- as.numeric(free)
  relative_step <- as.numeric(relative_step[1])
  baseline <- tryCatch(
    as.numeric(transform(free)),
    error = function(e) numeric(0)
  )
  if (length(free) == 0L || length(baseline) == 0L ||
      !is.finite(relative_step) || relative_step <= 0 ||
      any(!is.finite(free)) || any(!is.finite(baseline))) {
    return(list(
      valid = FALSE,
      jacobian = matrix(numeric(0), nrow = length(baseline), ncol = length(free)),
      relative_step = relative_step,
      absolute_steps = numeric(0)
    ))
  }

  absolute_steps <- relative_step * pmax(1, abs(free))
  jacobian <- matrix(
    NA_real_, nrow = length(baseline), ncol = length(free)
  )
  for (column in seq_along(free)) {
    high <- low <- free
    high[column] <- high[column] + absolute_steps[column]
    low[column] <- low[column] - absolute_steps[column]
    high_value <- tryCatch(as.numeric(transform(high)),
                           error = function(e) numeric(0))
    low_value <- tryCatch(as.numeric(transform(low)),
                          error = function(e) numeric(0))
    if (length(high_value) != length(baseline) ||
        length(low_value) != length(baseline) ||
        any(!is.finite(high_value)) || any(!is.finite(low_value))) {
      next
    }
    jacobian[, column] <-
      (high_value - low_value) / (2 * absolute_steps[column])
  }

  list(
    valid = all(is.finite(jacobian)),
    jacobian = jacobian,
    relative_step = relative_step,
    absolute_steps = absolute_steps
  )
}

mfrmr_transformation_rank_ladder <- function(jacobian,
                                             tolerances = c(1e-10, 1e-8, 1e-6)) {
  jacobian <- suppressWarnings(as.matrix(jacobian))
  if (!is.numeric(jacobian) || ncol(jacobian) == 0L ||
      nrow(jacobian) == 0L || any(!is.finite(jacobian))) {
    return(list(
      valid = FALSE,
      rank_ladder = data.frame(),
      singular_values = numeric(0),
      condition_index = NA_real_
    ))
  }
  singular_values <- tryCatch(
    as.numeric(base::svd(jacobian, nu = 0L, nv = 0L)$d),
    error = function(e) numeric(0)
  )
  if (length(singular_values) != min(dim(jacobian)) ||
      any(!is.finite(singular_values))) {
    return(list(
      valid = FALSE,
      rank_ladder = data.frame(),
      singular_values = singular_values,
      condition_index = NA_real_
    ))
  }
  scale <- max(singular_values)
  thresholds <- if (is.finite(scale) && scale > 0) {
    scale * as.numeric(tolerances)
  } else {
    as.numeric(tolerances)
  }
  ranks <- vapply(
    thresholds,
    function(threshold) sum(singular_values > threshold),
    integer(1)
  )
  positive <- singular_values[singular_values > 0]
  list(
    valid = TRUE,
    rank_ladder = data.frame(
      RelativeTolerance = as.numeric(tolerances),
      AbsoluteThreshold = thresholds,
      Rank = ranks,
      Nullity = as.integer(ncol(jacobian) - ranks),
      stringsAsFactors = FALSE
    ),
    singular_values = singular_values,
    condition_index = if (length(positive) > 0L) {
      max(singular_values) / min(positive)
    } else {
      Inf
    }
  )
}

mfrmr_empty_transformation_summary <- function() {
  data.frame(
    Block = character(0),
    CoordinateSystem = character(0),
    FreeCoordinates = integer(0),
    ExpandedCoordinates = integer(0),
    ExpectedRank = integer(0),
    LogCoordinateRank = integer(0),
    NaturalCoordinateRank = integer(0),
    LogFullColumnRank = logical(0),
    NaturalFullColumnRank = logical(0),
    NaturalMinimum = numeric(0),
    NaturalMaximum = numeric(0),
    InvariantResidual = numeric(0),
    MaxAbsLogJacobianDifference = numeric(0),
    MaxAbsNaturalJacobianDifference = numeric(0),
    MaxScaledNaturalJacobianDifference = numeric(0),
    LogConditionIndex = numeric(0),
    NaturalConditionIndex = numeric(0),
    Finite = logical(0),
    Status = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_invalid_transformation_summary <- function(block,
                                                  coordinate_system,
                                                  free_coordinates,
                                                  expanded_coordinates) {
  out <- mfrmr_empty_transformation_summary()[NA_integer_, , drop = FALSE]
  out$Block <- as.character(block)
  out$CoordinateSystem <- as.character(coordinate_system)
  out$FreeCoordinates <- as.integer(free_coordinates)
  out$ExpandedCoordinates <- as.integer(expanded_coordinates)
  out$ExpectedRank <- as.integer(free_coordinates)
  out$Finite <- FALSE
  out$Status <- "invalid_configuration"
  out
}

mfrmr_transformation_block_audit <- function(
    block,
    coordinate_system,
    rank_labels,
    free,
    log_transform,
    natural_transform,
    analytic_log_jacobian,
    analytic_natural_jacobian,
    natural_values,
    invariant_residual,
    relative_step,
    tolerances) {
  numeric_log <- mfrmr_numeric_transformation_jacobian(
    log_transform, free, relative_step
  )
  numeric_natural <- mfrmr_numeric_transformation_jacobian(
    natural_transform, free, relative_step
  )
  log_rank <- mfrmr_transformation_rank_ladder(
    analytic_log_jacobian, tolerances
  )
  natural_rank <- mfrmr_transformation_rank_ladder(
    analytic_natural_jacobian, tolerances
  )
  primary_rank <- function(x) {
    if (isTRUE(x$valid)) as.integer(x$rank_ladder$Rank[1]) else NA_integer_
  }
  difference <- function(analytic, numeric) {
    if (!isTRUE(numeric$valid)) return(NA_real_)
    max(abs(analytic - numeric$jacobian))
  }
  scaled_difference <- function(analytic, numeric) {
    if (!isTRUE(numeric$valid)) return(NA_real_)
    max(
      abs(analytic - numeric$jacobian) /
        pmax(1, abs(analytic), abs(numeric$jacobian))
    )
  }
  log_rank_primary <- primary_rank(log_rank)
  natural_rank_primary <- primary_rank(natural_rank)
  free_coordinates <- length(free)
  finite <- all(is.finite(natural_values)) && all(natural_values > 0) &&
    isTRUE(numeric_log$valid) && isTRUE(numeric_natural$valid) &&
    isTRUE(log_rank$valid) && isTRUE(natural_rank$valid)
  summary <- data.frame(
    Block = block,
    CoordinateSystem = coordinate_system,
    FreeCoordinates = free_coordinates,
    ExpandedCoordinates = length(natural_values),
    ExpectedRank = free_coordinates,
    LogCoordinateRank = log_rank_primary,
    NaturalCoordinateRank = natural_rank_primary,
    LogFullColumnRank = isTRUE(log_rank_primary == free_coordinates),
    NaturalFullColumnRank = isTRUE(
      natural_rank_primary == free_coordinates
    ),
    NaturalMinimum = if (length(natural_values) > 0L) {
      min(natural_values)
    } else {
      NA_real_
    },
    NaturalMaximum = if (length(natural_values) > 0L) {
      max(natural_values)
    } else {
      NA_real_
    },
    InvariantResidual = as.numeric(invariant_residual),
    MaxAbsLogJacobianDifference = difference(
      analytic_log_jacobian, numeric_log
    ),
    MaxAbsNaturalJacobianDifference = difference(
      analytic_natural_jacobian, numeric_natural
    ),
    MaxScaledNaturalJacobianDifference = scaled_difference(
      analytic_natural_jacobian, numeric_natural
    ),
    LogConditionIndex = log_rank$condition_index,
    NaturalConditionIndex = natural_rank$condition_index,
    Finite = finite,
    Status = if (finite) "evaluated_diagnostic_only" else "unavailable",
    stringsAsFactors = FALSE
  )
  ladder_rows <- lapply(seq_along(rank_labels), function(index) {
    rank <- list(log_rank, natural_rank)[[index]]
    if (!isTRUE(rank$valid)) return(NULL)
    value <- rank$rank_ladder
    value$Block <- block
    value$CoordinateSystem <- rank_labels[index]
    value[, c(
      "Block", "CoordinateSystem", "RelativeTolerance",
      "AbsoluteThreshold", "Rank", "Nullity"
    )]
  })
  ladder_rows <- Filter(Negate(is.null), ladder_rows)
  list(
    summary = summary,
    rank_ladder = if (length(ladder_rows) > 0L) {
      do.call(rbind, ladder_rows)
    } else {
      data.frame()
    }
  )
}

mfrmr_nonlinear_transformation_audit <- function(
    par,
    sizes,
    config,
    nonlinear_blocks,
    relative_step = 1e-6,
    tolerances = c(1e-10, 1e-8, 1e-6)) {
  nonlinear_blocks <- as.character(nonlinear_blocks %||% character(0))
  base <- list(
    role = "free_to_expanded_nonlinear_parameter_transformation",
    evaluation_point = "retained_optimizer_free_coordinate_vector",
    attempted = FALSE,
    status = "not_required",
    nonlinear_blocks = nonlinear_blocks,
    numerical_differentiation = list(
      method = "coordinate_scaled_central_difference",
      relative_step = as.numeric(relative_step[1])
    ),
    block_summary = mfrmr_empty_transformation_summary(),
    rank_ladder = data.frame(),
    parameterization_only = TRUE,
    likelihood_jacobian_evaluated = FALSE,
    structural_identification_classified = FALSE,
    readiness_effect = "none_parameterization_audit_only",
    detail = paste(
      "No nonlinear free-coordinate transformation was required for this fit."
    )
  )
  if (length(nonlinear_blocks) == 0L) return(base)

  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  if (is.null(par) || length(par) != free_dimension || any(!is.finite(par))) {
    base$status <- "not_evaluated_parameter_vector"
    base$detail <- paste(
      "The nonlinear transformation audit was not evaluated because the",
      "retained free parameter vector was unavailable, non-finite, or",
      "dimensionally inconsistent."
    )
    return(base)
  }

  slices <- build_param_slices(sizes)
  results <- list()

  if ("log_slopes" %in% nonlinear_blocks) {
    slice <- as.integer(slices$log_slopes %||% integer(0))
    free <- as.numeric(par[slice])
    levels <- as.character(config$gpcm_spec$levels %||% character(0))
    n_levels <- length(levels)
    valid_config <- identical(config$model, "GPCM") &&
      length(slice) > 0L && n_levels == length(free) + 1L
    if (valid_config) {
      transformed <- expand_gpcm_log_slopes(free, config$gpcm_spec)
      log_jacobian <- sum_zero_jacobian(n_levels)
      natural_jacobian <- diag(
        transformed$slopes, nrow = n_levels
      ) %*% log_jacobian
      results[[length(results) + 1L]] <- mfrmr_transformation_block_audit(
        block = "log_slopes",
        coordinate_system = paste(
          "free_log_slopes -> sum_zero_log_slopes ->",
          "positive_geometric_mean_one_slopes"
        ),
        rank_labels = c(
          "expanded_sum_zero_log_slopes",
          "positive_geometric_mean_one_slopes"
        ),
        free = free,
        log_transform = function(value) c(value, -sum(value)),
        natural_transform = function(value) exp(c(value, -sum(value))),
        analytic_log_jacobian = log_jacobian,
        analytic_natural_jacobian = natural_jacobian,
        natural_values = transformed$slopes,
        invariant_residual = abs(mean(log(transformed$slopes))),
        relative_step = relative_step,
        tolerances = tolerances
      )
    } else {
      results[[length(results) + 1L]] <- list(
        summary = mfrmr_invalid_transformation_summary(
          block = "log_slopes",
          coordinate_system = "invalid_GPCM_slope_parameterization",
          free_coordinates = length(free),
          expanded_coordinates = n_levels
        ),
        rank_ladder = data.frame()
      )
    }
  }

  if ("log_sigma2" %in% nonlinear_blocks) {
    slice <- as.integer(slices$log_sigma2 %||% integer(0))
    free <- as.numeric(par[slice])
    valid_config <- identical(config$method, "MML") &&
      isTRUE(config$population_spec$active) && length(free) == 1L
    if (valid_config) {
      sigma2 <- exp(free)
      log_jacobian <- matrix(1, nrow = 1L, ncol = 1L)
      natural_jacobian <- matrix(sigma2, nrow = 1L, ncol = 1L)
      results[[length(results) + 1L]] <- mfrmr_transformation_block_audit(
        block = "log_sigma2",
        coordinate_system =
          "free_log_sigma2 -> positive_residual_variance",
        rank_labels = c(
          "log_residual_variance", "positive_residual_variance"
        ),
        free = free,
        log_transform = identity,
        natural_transform = exp,
        analytic_log_jacobian = log_jacobian,
        analytic_natural_jacobian = natural_jacobian,
        natural_values = sigma2,
        invariant_residual = 0,
        relative_step = relative_step,
        tolerances = tolerances
      )
    } else {
      results[[length(results) + 1L]] <- list(
        summary = mfrmr_invalid_transformation_summary(
          block = "log_sigma2",
          coordinate_system = "invalid_latent_variance_parameterization",
          free_coordinates = length(free),
          expanded_coordinates = 1L
        ),
        rank_ladder = data.frame()
      )
    }
  }

  base$attempted <- TRUE
  base$block_summary <- if (length(results) > 0L) {
    value <- do.call(rbind, lapply(results, `[[`, "summary"))
    rownames(value) <- NULL
    value
  } else {
    mfrmr_empty_transformation_summary()
  }
  ladder_rows <- Filter(
    function(x) is.data.frame(x) && nrow(x) > 0L,
    lapply(results, `[[`, "rank_ladder")
  )
  base$rank_ladder <- if (length(ladder_rows) > 0L) {
    value <- do.call(rbind, ladder_rows)
    rownames(value) <- NULL
    value
  } else {
    data.frame()
  }
  base$status <- if (nrow(base$block_summary) > 0L &&
                     all(base$block_summary$Status ==
                           "evaluated_diagnostic_only")) {
    "evaluated_diagnostic_only"
  } else {
    "partially_unavailable"
  }
  base$detail <- paste(
    "This record verifies the free-to-expanded nonlinear parameter",
    "transformations and their numerical derivatives only. It is not a",
    "response-likelihood Jacobian, structural-identification classification,",
    "weak-information rule, or readiness decision."
  )
  base
}

mfrmr_gpcm_log_slope_map <- function(config, sizes) {
  levels <- as.character(config$gpcm_spec$levels %||% character(0))
  n_free <- as.integer(sizes$log_slopes %||% 0L)
  if (n_free == 0L) return(mfrmr_estimability_map())
  if (length(levels) != n_free + 1L) {
    stop(
      "Internal GPCM response-kernel audit found an inconsistent slope map.",
      call. = FALSE
    )
  }
  slices <- build_param_slices(sizes)
  mfrmr_estimability_map(
    Block = rep("log_slopes", n_free),
    Coordinate = paste0("log_slopes:", levels[seq_len(n_free)]),
    Facet = rep(as.character(config$slope_facet), n_free),
    Level = levels[seq_len(n_free)],
    ReferenceLevel = rep(levels[length(levels)], n_free),
    Constraint = rep("sum_zero_log_slopes", n_free),
    OptimizerIndex = as.integer(slices$log_slopes)
  )
}

mfrmr_gpcm_selected_log_slope_jacobian_sparse <- function(slope_index,
                                                            n_levels) {
  slope_index <- as.integer(slope_index)
  n_levels <- as.integer(n_levels)
  n_free <- max(n_levels - 1L, 0L)
  if (n_free == 0L) {
    return(mfrmr_empty_sparse_matrix(length(slope_index), 0L))
  }
  if (anyNA(slope_index) || any(slope_index < 1L) ||
      any(slope_index > n_levels)) {
    stop("The GPCM response-kernel audit found an invalid slope index.",
         call. = FALSE)
  }
  direct_rows <- which(slope_index <= n_free)
  reference_rows <- which(slope_index == n_levels)
  Matrix::sparseMatrix(
    i = c(
      direct_rows,
      rep(reference_rows, each = n_free)
    ),
    j = c(
      slope_index[direct_rows],
      rep(seq_len(n_free), times = length(reference_rows))
    ),
    x = c(
      rep(1, length(direct_rows)),
      rep(-1, length(reference_rows) * n_free)
    ),
    dims = c(length(slope_index), n_free)
  )
}

mfrmr_gpcm_adjacent_logits <- function(par, idx, config, sizes) {
  if (!identical(config$model, "GPCM") ||
      !identical(config$method, "JML")) {
    stop(
      "The retained adjacent-logit response kernel is defined here only for JML GPCM.",
      call. = FALSE
    )
  }
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)
  n_obs <- length(eta)
  n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
  if (n_obs == 0L || n_steps <= 0L) {
    stop("The GPCM response-kernel audit requires retained observations and transitions.",
         call. = FALSE)
  }
  step_idx <- as.integer(idx$step_idx)
  slope_idx <- as.integer(idx$slope_idx)
  if (length(step_idx) != n_obs || length(slope_idx) != n_obs ||
      anyNA(step_idx) || anyNA(slope_idx)) {
    stop("The GPCM response-kernel audit could not align step or slope indices.",
         call. = FALSE)
  }
  step_obs <- params$steps_mat[step_idx, , drop = FALSE]
  eta_minus_step <- matrix(eta, nrow = n_obs, ncol = n_steps) - step_obs
  slope_obs <- as.numeric(params$slopes[slope_idx])
  as.vector(eta_minus_step * matrix(
    slope_obs, nrow = n_obs, ncol = n_steps
  ))
}

mfrmr_gpcm_response_kernel_design <- function(prep, idx, config, sizes, par) {
  if (!identical(config$model, "GPCM") ||
      !identical(config$method, "JML")) {
    stop(
      "The retained response-kernel Jacobian is available only for JML GPCM.",
      call. = FALSE
    )
  }
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  if (is.null(par) || length(par) != free_dimension || any(!is.finite(par))) {
    stop("The retained GPCM free parameter vector is unavailable or malformed.",
         call. = FALSE)
  }

  additive <- mfrmr_estimability_adjacent_design(
    prep = prep,
    idx = idx,
    config = config,
    sizes = sizes,
    include_person = TRUE,
    include_population_beta = FALSE
  )
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)
  n_obs <- length(eta)
  n_steps <- additive$transitions
  obs_index <- rep(seq_len(n_obs), times = n_steps)
  transition <- rep(seq_len(n_steps), each = n_obs)
  step_idx <- as.integer(idx$step_idx)
  slope_idx <- as.integer(idx$slope_idx)
  step_value <- params$steps_mat[cbind(step_idx[obs_index], transition)]
  eta_minus_step <- eta[obs_index] - step_value
  slope_repeated <- as.numeric(params$slopes[slope_idx[obs_index]])

  additive_jacobian <- Matrix::Diagonal(x = slope_repeated) %*%
    additive$design
  n_slope_levels <- length(params$slopes)
  slope_map <- mfrmr_gpcm_log_slope_map(config, sizes)
  slope_jacobian <- if (n_slope_levels <= 1L) {
    mfrmr_empty_sparse_matrix(nrow(additive$design), 0L)
  } else {
    selected <- mfrmr_gpcm_selected_log_slope_jacobian_sparse(
      slope_index = slope_idx[obs_index],
      n_levels = n_slope_levels
    )
    Matrix::Diagonal(x = slope_repeated * eta_minus_step) %*% selected
  }
  jacobian <- methods::as(
    cbind(additive_jacobian, slope_jacobian), "dgCMatrix"
  )
  map <- rbind(additive$map, slope_map)
  rownames(map) <- NULL
  if (ncol(jacobian) != free_dimension || nrow(map) != free_dimension ||
      !identical(as.integer(map$OptimizerIndex), seq_len(free_dimension))) {
    stop(
      "Internal GPCM response-kernel Jacobian does not match optimizer coordinates.",
      call. = FALSE
    )
  }
  list(
    jacobian = jacobian,
    parameter_map = map,
    adjacent_logits = mfrmr_gpcm_adjacent_logits(par, idx, config, sizes),
    observation_rows = n_obs,
    transition_rows = nrow(jacobian),
    transitions = n_steps,
    additive_free_dimension = ncol(additive$design),
    log_slope_free_dimension = ncol(slope_jacobian),
    eta_minus_step = as.numeric(eta_minus_step),
    slope_repeated = as.numeric(slope_repeated)
  )
}

audit_mfrm_gpcm_response_kernel <- function(
    prep,
    idx,
    config,
    sizes,
    par,
    max_numeric_free_dimension = 80L,
    max_numeric_elements = 2e6,
    relative_step = 1e-6,
    tolerances = c(1e-12, 1e-10, 1e-8)) {
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  base <- list(
    role = "retained_gpcm_adjacent_category_logit_response_kernel_jacobian",
    response_scale = "log(P[k]/P[k-1])",
    evaluation_point = "retained_optimizer_free_coordinate_vector",
    method = as.character(config$method),
    model = as.character(config$model),
    attempted = FALSE,
    status = if (identical(config$model, "GPCM")) {
      "pending_after_optimization"
    } else {
      "not_applicable_model"
    },
    kernel_scope = if (identical(config$method, "JML")) {
      "conditional_response_given_retained_person_coordinates"
    } else {
      "person_integrated_response_pattern_required"
    },
    free_dimension = as.integer(free_dimension),
    observation_rows = 0L,
    transition_rows = 0L,
    transitions = max(as.integer(config$n_cat %||% 0L) - 1L, 0L),
    nonzero_entries = 0L,
    local_rank = NA_integer_,
    local_nullity = NA_integer_,
    local_rank_state = "not_evaluated",
    tolerance_sensitive = FALSE,
    rank_ladder = data.frame(),
    parameter_blocks = data.frame(
      Block = character(0), FreeCoordinates = integer(0)
    ),
    parameter_map = mfrmr_estimability_map(),
    zero_coordinates = character(0),
    null_directions = data.frame(),
    numerical_differentiation = list(
      status = "not_evaluated",
      method = "coordinate_scaled_central_difference",
      relative_step = as.numeric(relative_step[1]),
      free_dimension_limit = as.integer(max_numeric_free_dimension),
      element_limit = as.double(max_numeric_elements),
      max_abs_difference = NA_real_,
      max_scaled_difference = NA_real_
    ),
    conditional_response_kernel_jacobian_evaluated = FALSE,
    marginal_person_pattern_jacobian_evaluated = FALSE,
    structural_identification_classified = FALSE,
    readiness_effect = "none_pending_property_and_marginal_model_audits",
    detail = "This response-kernel audit is not applicable to the fitted model."
  )
  if (!identical(config$model, "GPCM")) return(base)

  if (identical(config$method, "MML")) {
    base$status <- "not_evaluated_marginal_person_pattern_required"
    base$detail <- paste(
      "MML integrates Person coordinates, so a conditional observation-level",
      "adjacent-logit Jacobian cannot certify the marginal person-pattern",
      "model. The fitted-information Hessian is a separate local likelihood",
      "diagnostic. The separate observed-pattern and exhaustive retained-",
      "design expected-information audits remain local and do not supply a",
      "global marginal structural-identification result."
    )
    return(base)
  }
  if (!identical(config$method, "JML")) {
    base$status <- "not_evaluated_estimator_scope"
    base$detail <- "The estimator is outside the implemented response-kernel scope."
    return(base)
  }

  built <- tryCatch(
    mfrmr_gpcm_response_kernel_design(prep, idx, config, sizes, par),
    error = function(e) e
  )
  base$attempted <- TRUE
  if (inherits(built, "error")) {
    base$status <- "unavailable"
    base$detail <- paste0(
      "The retained JML GPCM response-kernel Jacobian was unavailable: ",
      conditionMessage(built), "."
    )
    return(base)
  }

  rank <- mfrmr_estimability_rank_audit(
    built$jacobian, built$parameter_map,
    tolerances = tolerances,
    structural_tolerance = tolerances[1]
  )
  block_summary <- if (nrow(built$parameter_map) == 0L) {
    base$parameter_blocks
  } else {
    stats::aggregate(
      rep(1L, nrow(built$parameter_map)),
      by = list(Block = built$parameter_map$Block),
      FUN = sum
    ) |>
      stats::setNames(c("Block", "FreeCoordinates"))
  }
  base$status <- "evaluated_local_diagnostic_only"
  base$observation_rows <- as.integer(built$observation_rows)
  base$transition_rows <- as.integer(built$transition_rows)
  base$transitions <- as.integer(built$transitions)
  base$nonzero_entries <- length(built$jacobian@x)
  base$local_rank <- rank$Rank
  base$local_nullity <- rank$Nullity
  base$local_rank_state <- if (rank$Nullity == 0L) {
    "locally_full_column_rank"
  } else {
    "locally_rank_deficient"
  }
  base$tolerance_sensitive <- isTRUE(rank$ToleranceSensitive)
  base$rank_ladder <- rank$ToleranceRanks
  base$parameter_blocks <- block_summary
  base$parameter_map <- built$parameter_map
  base$zero_coordinates <- rank$ZeroCoordinates
  base$null_directions <- rank$NullDirections
  base$conditional_response_kernel_jacobian_evaluated <- TRUE

  numeric_elements <- as.double(built$transition_rows) *
    as.double(free_dimension)
  if (free_dimension <= as.integer(max_numeric_free_dimension) &&
      numeric_elements <= as.double(max_numeric_elements)) {
    numeric <- mfrmr_numeric_transformation_jacobian(
      transform = function(value) {
        mfrmr_gpcm_adjacent_logits(value, idx, config, sizes)
      },
      free = par,
      relative_step = relative_step
    )
    if (isTRUE(numeric$valid) &&
        identical(dim(numeric$jacobian), dim(built$jacobian))) {
      analytic <- as.matrix(built$jacobian)
      absolute_difference <- abs(analytic - numeric$jacobian)
      scaled_difference <- absolute_difference /
        pmax(1, abs(analytic), abs(numeric$jacobian))
      base$numerical_differentiation$status <- "evaluated"
      base$numerical_differentiation$max_abs_difference <-
        max(absolute_difference)
      base$numerical_differentiation$max_scaled_difference <-
        max(scaled_difference)
    } else {
      base$numerical_differentiation$status <- "unavailable"
    }
  } else {
    base$numerical_differentiation$status <-
      "not_evaluated_execution_limit"
  }
  base$detail <- paste(
    "The analytic Jacobian combines constrained Person/facet/interaction/step",
    "coordinates with sum-zero log-slope coordinates for every retained JML",
    "GPCM adjacent-category logit. Its local rank and independent derivative",
    "check are diagnostic only; they do not yet classify structural or weak",
    "identification or alter readiness."
  )
  base
}

mfrmr_subset_observation_indices <- function(idx, rows) {
  rows <- as.integer(rows)
  n_obs <- length(idx$score_k)
  if (anyNA(rows) || any(rows < 1L) || any(rows > n_obs)) {
    stop("Internal MML pattern-score audit received invalid observation rows.",
         call. = FALSE)
  }
  subset_vector <- function(value) {
    if (is.null(value)) return(NULL)
    if (length(value) == n_obs) value[rows] else value
  }
  out <- idx
  out$person <- subset_vector(idx$person)
  out$score_k <- subset_vector(idx$score_k)
  out$weight <- subset_vector(idx$weight)
  out$step_idx <- subset_vector(idx$step_idx)
  out$slope_idx <- subset_vector(idx$slope_idx)
  out$facets <- lapply(idx$facets %||% list(), subset_vector)
  out$interactions <- lapply(idx$interactions %||% list(), subset_vector)
  out$criterion_splits <- if (!is.null(out$step_idx)) {
    split(seq_along(out$step_idx), out$step_idx)
  } else {
    NULL
  }
  out$slope_splits <- if (!is.null(out$slope_idx)) {
    split(seq_along(out$slope_idx), out$slope_idx)
  } else {
    NULL
  }
  out
}

mfrmr_mml_person_log_marginals <- function(par, idx, config, sizes, quad) {
  if (!identical(config$method, "MML")) {
    stop("Person-pattern marginal contributions are defined here only for MML.",
         call. = FALSE)
  }
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  if (length(par) != free_dimension || any(!is.finite(par))) {
    stop("The MML free parameter vector is unavailable or malformed.",
         call. = FALSE)
  }
  params <- expand_params(par, sizes, config)
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
  list(
    person_ids = as.integer(person_bundle$person_ids),
    log_marginal = as.numeric(person_bundle$log_marginal)
  )
}

mfrmr_mml_observed_person_score_matrix <- function(par, idx, config, sizes,
                                                    quad) {
  contribution <- mfrmr_mml_person_log_marginals(
    par = par, idx = idx, config = config, sizes = sizes, quad = quad
  )
  params <- expand_params(par, sizes, config)
  score_rows <- lapply(contribution$person_ids, function(person_id) {
    rows <- which(as.integer(idx$person) == person_id)
    person_idx <- mfrmr_subset_observation_indices(idx, rows)
    person_base_eta <- compute_base_eta(person_idx, params, config)
    # `mfrm_grad_mml_core()` returns the gradient of the negative marginal
    # log-likelihood. Negating it gives the conventional log-likelihood score
    # for this one observed Person response pattern.
    -mfrm_grad_mml_core(
      params = params,
      base_eta = person_base_eta,
      idx = person_idx,
      config = config,
      sizes = sizes,
      quad = quad
    )
  })
  score <- do.call(rbind, score_rows)
  if (is.null(dim(score))) {
    score <- matrix(score, nrow = length(contribution$person_ids), byrow = TRUE)
  }
  score <- unname(as.matrix(score))
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  if (nrow(score) != length(contribution$person_ids) ||
      ncol(score) != free_dimension || any(!is.finite(score))) {
    stop("Internal MML Person-pattern score matrix is malformed.",
         call. = FALSE)
  }
  list(
    score = score,
    person_ids = contribution$person_ids,
    log_marginal = contribution$log_marginal
  )
}

mfrmr_enumerate_response_patterns <- function(n_observations, n_categories) {
  n_observations <- suppressWarnings(as.integer(n_observations)[1])
  n_categories <- suppressWarnings(as.integer(n_categories)[1])
  if (!is.finite(n_observations) || n_observations < 0L) {
    stop("`n_observations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(n_categories) || n_categories < 2L) {
    stop("`n_categories` must be one integer of at least two.", call. = FALSE)
  }
  if (n_observations == 0L) {
    return(matrix(integer(0), nrow = 1L, ncol = 0L))
  }
  pattern_count <- n_categories^n_observations
  if (!is.finite(pattern_count) || pattern_count > .Machine$integer.max) {
    stop("The requested response-pattern grid is too large to enumerate.",
         call. = FALSE)
  }
  code <- 0:(as.integer(pattern_count) - 1L)
  radix <- n_categories^(seq_len(n_observations) - 1L)
  patterns <- outer(
    code, radix,
    FUN = function(value, place) floor(value / place) %% n_categories
  )
  storage.mode(patterns) <- "integer"
  unname(patterns)
}

mfrmr_mml_evaluate_person_patterns <- function(
    par,
    person_idx,
    config,
    sizes,
    quad,
    patterns,
    include_scores = TRUE) {
  if (!identical(config$method, "MML")) {
    stop("All-pattern marginal evaluation is defined here only for MML.",
         call. = FALSE)
  }
  patterns <- as.matrix(patterns)
  n_obs <- length(person_idx$score_k)
  if (n_obs <= 0L || ncol(patterns) != n_obs || nrow(patterns) <= 0L) {
    stop("The Person response-pattern matrix is empty or dimensionally invalid.",
         call. = FALSE)
  }
  person_ids <- unique(as.integer(person_idx$person))
  if (length(person_ids) != 1L || anyNA(person_ids)) {
    stop("All-pattern marginal evaluation requires exactly one observed Person design.",
         call. = FALSE)
  }
  n_cat <- as.integer(config$n_cat)
  if (any(!is.finite(patterns)) || any(patterns < 0L) ||
      any(patterns >= n_cat) || any(patterns != floor(patterns))) {
    stop("The response-pattern matrix contains an invalid category index.",
         call. = FALSE)
  }
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  if (length(par) != free_dimension || any(!is.finite(par))) {
    stop("The MML free parameter vector is unavailable or malformed.",
         call. = FALSE)
  }

  params <- expand_params(par, sizes, config)
  base_eta <- compute_base_eta(person_idx, params, config)
  log_marginal <- numeric(nrow(patterns))
  score <- if (isTRUE(include_scores)) {
    matrix(NA_real_, nrow = nrow(patterns), ncol = free_dimension)
  } else {
    matrix(numeric(0), nrow = 0L, ncol = free_dimension)
  }

  for (pattern_row in seq_len(nrow(patterns))) {
    pattern_idx <- person_idx
    pattern_idx$score_k <- as.integer(patterns[pattern_row, ])
    logprob_bundle <- mfrm_mml_logprob_bundle(
      idx = pattern_idx,
      config = config,
      quad = quad,
      params = params,
      base_eta = base_eta,
      include_probs = isTRUE(include_scores),
      include_linear_part = isTRUE(include_scores) &&
        identical(config$model, "GPCM")
    )
    posterior_bundle <- mfrm_mml_posterior_bundle(logprob_bundle)
    marginal <- posterior_bundle$person_bundle$log_marginal
    if (length(marginal) != 1L || !is.finite(marginal)) {
      stop("A Person response-pattern marginal probability was non-finite.",
           call. = FALSE)
    }
    log_marginal[pattern_row] <- marginal
    if (isTRUE(include_scores)) {
      score[pattern_row, ] <- -mfrm_grad_mml_core(
        params = params,
        base_eta = base_eta,
        idx = pattern_idx,
        config = config,
        sizes = sizes,
        quad = quad,
        logprob_bundle = logprob_bundle,
        posterior_bundle = posterior_bundle
      )
    }
  }

  if (isTRUE(include_scores) && any(!is.finite(score))) {
    stop("An all-pattern marginal score was non-finite.", call. = FALSE)
  }
  list(
    person_id = person_ids,
    patterns = patterns,
    log_marginal = log_marginal,
    score = score
  )
}

mfrmr_mml_person_design_groups <- function(
    idx,
    config,
    reuse_identical_designs = TRUE) {
  person_ids <- sort(unique(as.integer(idx$person)))
  person_ids <- person_ids[is.finite(person_ids)]
  canonical_rows <- vector("list", length(person_ids))
  signatures <- character(length(person_ids))
  population_spec <- config$population_spec %||% list(active = FALSE)

  for (person_position in seq_along(person_ids)) {
    person_id <- person_ids[person_position]
    rows <- which(as.integer(idx$person) == person_id)
    order_keys <- list()
    if (!is.null(idx$step_idx)) {
      order_keys[[length(order_keys) + 1L]] <- as.integer(idx$step_idx[rows])
    }
    if (!is.null(idx$slope_idx)) {
      order_keys[[length(order_keys) + 1L]] <- as.integer(idx$slope_idx[rows])
    }
    for (facet in names(idx$facets %||% list())) {
      order_keys[[length(order_keys) + 1L]] <-
        as.integer(idx$facets[[facet]][rows])
    }
    for (interaction in names(idx$interactions %||% list())) {
      order_keys[[length(order_keys) + 1L]] <-
        as.integer(idx$interactions[[interaction]][rows])
    }
    within_person_order <- if (length(order_keys) > 0L) {
      do.call(
        order,
        c(order_keys, list(na.last = TRUE, method = "radix"))
      )
    } else {
      seq_along(rows)
    }
    rows <- rows[within_person_order]
    canonical_rows[[person_position]] <- as.integer(rows)
    person_idx <- mfrmr_subset_observation_indices(idx, rows)

    population_row <- NULL
    if (isTRUE(population_spec$active)) {
      lookup <- as.integer(population_spec$person_lookup %||% integer(0))
      design <- population_spec$design_matrix
      if (length(lookup) < person_id || is.na(lookup[person_id]) ||
          is.null(design) || !is.matrix(design) ||
          lookup[person_id] < 1L || lookup[person_id] > nrow(design)) {
        stop("MML all-pattern design reuse could not align a population-design row.",
             call. = FALSE)
      }
      population_row <- as.numeric(design[lookup[person_id], , drop = TRUE])
    }
    components <- list(
      observation_rows = length(rows),
      facets = lapply(person_idx$facets %||% list(), as.integer),
      interactions = lapply(
        person_idx$interactions %||% list(), as.integer
      ),
      step_idx = as.integer(person_idx$step_idx %||% integer(0)),
      slope_idx = as.integer(person_idx$slope_idx %||% integer(0)),
      population_design_row = population_row
    )
    serialized <- serialize(components, connection = NULL, version = 2)
    signatures[person_position] <- paste(as.integer(serialized), collapse = ".")
    if (!isTRUE(reuse_identical_designs)) {
      signatures[person_position] <- paste0(
        signatures[person_position], "|person_position=", person_position
      )
    }
  }

  unique_signatures <- unique(signatures)
  group_index <- match(signatures, unique_signatures)
  group_positions <- split(
    seq_along(person_ids),
    factor(group_index, levels = seq_along(unique_signatures))
  )
  list(
    person_ids = as.integer(person_ids),
    canonical_rows = canonical_rows,
    group_positions = unname(group_positions),
    representative_positions = as.integer(vapply(
      group_positions, `[[`, integer(1), 1L
    )),
    group_sizes = as.integer(lengths(group_positions)),
    unique_person_designs = as.integer(length(group_positions)),
    reuse_identical_designs = isTRUE(reuse_identical_designs)
  )
}

mfrmr_mml_all_pattern_workload <- function(idx, n_categories) {
  person_ids <- sort(unique(as.integer(idx$person)))
  person_ids <- person_ids[is.finite(person_ids)]
  observation_counts <- vapply(person_ids, function(person_id) {
    sum(as.integer(idx$person) == person_id)
  }, integer(1))
  log_pattern_counts <- observation_counts * log(as.numeric(n_categories))
  pattern_counts <- ifelse(
    log_pattern_counts <= log(.Machine$double.xmax),
    exp(log_pattern_counts),
    Inf
  )
  pattern_counts <- round(pattern_counts)
  list(
    person_ids = as.integer(person_ids),
    observation_counts = as.integer(observation_counts),
    pattern_counts = as.double(pattern_counts),
    total_patterns = sum(pattern_counts)
  )
}

mfrmr_mml_all_pattern_expected_information <- function(
    par,
    idx,
    config,
    sizes,
    quad,
    reuse_identical_designs = TRUE) {
  if (!identical(config$method, "MML")) {
    stop("All-pattern expected information is defined here only for MML.",
         call. = FALSE)
  }
  weights <- idx$weight
  if (!is.null(weights) &&
      (any(!is.finite(weights)) || any(abs(weights - 1) > 1e-12))) {
    stop("All-pattern expected information currently requires unit row weights.",
         call. = FALSE)
  }
  workload <- mfrmr_mml_all_pattern_workload(idx, config$n_cat)
  design_groups <- mfrmr_mml_person_design_groups(
    idx = idx,
    config = config,
    reuse_identical_designs = reuse_identical_designs
  )
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  by_person <- vector("list", length(workload$person_ids))
  weighted_score <- vector("list", design_groups$unique_person_designs)
  expected_score <- matrix(
    0,
    nrow = length(workload$person_ids),
    ncol = free_dimension
  )
  probability_mass <- numeric(length(workload$person_ids))
  evaluated_pattern_counts <- numeric(design_groups$unique_person_designs)

  for (group_position in seq_len(design_groups$unique_person_designs)) {
    member_positions <- design_groups$group_positions[[group_position]]
    representative_position <- member_positions[1L]
    person_id <- workload$person_ids[representative_position]
    rows <- design_groups$canonical_rows[[representative_position]]
    person_idx <- mfrmr_subset_observation_indices(idx, rows)
    patterns <- mfrmr_enumerate_response_patterns(
      n_observations = length(rows),
      n_categories = config$n_cat
    )
    evaluated <- mfrmr_mml_evaluate_person_patterns(
      par = par,
      person_idx = person_idx,
      config = config,
      sizes = sizes,
      quad = quad,
      patterns = patterns,
      include_scores = TRUE
    )
    probability <- exp(evaluated$log_marginal)
    if (any(!is.finite(probability)) || any(probability < 0)) {
      stop("An all-pattern marginal probability was invalid.", call. = FALSE)
    }
    mass <- sum(probability)
    one_person_expected_score <- colSums(evaluated$score * probability)
    probability_mass[member_positions] <- mass
    expected_score[member_positions, ] <- matrix(
      one_person_expected_score,
      nrow = length(member_positions),
      ncol = free_dimension,
      byrow = TRUE
    )
    weighted_score[[group_position]] <- evaluated$score *
      sqrt(probability * length(member_positions))
    evaluated_pattern_counts[group_position] <- nrow(evaluated$patterns)
    evaluated$probability <- probability
    for (person_position in member_positions) {
      person_evaluated <- evaluated
      person_evaluated$person_id <- workload$person_ids[person_position]
      by_person[[person_position]] <- person_evaluated
    }
  }

  weighted_score <- do.call(rbind, weighted_score)
  expected_information <- crossprod(weighted_score)
  evaluated_patterns <- sum(evaluated_pattern_counts)
  if (nrow(weighted_score) != evaluated_patterns ||
      ncol(weighted_score) != free_dimension ||
      any(!is.finite(expected_information))) {
    stop("The all-pattern expected-information assembly is malformed.",
         call. = FALSE)
  }
  list(
    person_ids = workload$person_ids,
    observation_counts = workload$observation_counts,
    pattern_counts = workload$pattern_counts,
    total_patterns = workload$total_patterns,
    evaluated_pattern_counts = as.double(evaluated_pattern_counts),
    evaluated_patterns = as.double(evaluated_patterns),
    unique_person_designs = design_groups$unique_person_designs,
    design_group_sizes = design_groups$group_sizes,
    canonical_rows = design_groups$canonical_rows,
    reuse_identical_designs = design_groups$reuse_identical_designs,
    probability_mass = probability_mass,
    expected_score = expected_score,
    weighted_score = unname(weighted_score),
    expected_information = unname(expected_information),
    by_person = by_person
  )
}

mfrmr_mml_optimizer_parameter_map <- function(prep, idx, config, sizes) {
  if (!identical(config$method, "MML")) {
    stop("The MML optimizer map requires method = 'MML'.", call. = FALSE)
  }
  eta <- mfrmr_estimability_eta_design(
    prep = prep,
    idx = idx,
    config = config,
    sizes = sizes,
    include_person = FALSE,
    include_population_beta = TRUE
  )
  step <- mfrmr_step_jacobian_sparse(config, sizes)
  step$map <- mfrmr_estimability_set_optimizer_index(
    step$map, "steps", build_param_slices(sizes)
  )
  maps <- list(eta$map, step$map)
  if (identical(config$model, "GPCM")) {
    maps[[length(maps) + 1L]] <- mfrmr_gpcm_log_slope_map(config, sizes)
  }
  if (as.integer(sizes$log_sigma2 %||% 0L) > 0L) {
    slices <- build_param_slices(sizes)
    maps[[length(maps) + 1L]] <- mfrmr_estimability_map(
      Block = "log_sigma2",
      Coordinate = "log_sigma2:residual_variance",
      Facet = "PersonPopulation",
      Level = "residual_variance",
      ReferenceLevel = "",
      Constraint = "positive_via_log_scale",
      OptimizerIndex = as.integer(slices$log_sigma2)
    )
  }
  map <- do.call(rbind, maps)
  rownames(map) <- NULL
  map <- map[order(map$OptimizerIndex), , drop = FALSE]
  rownames(map) <- NULL
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  if (nrow(map) != free_dimension ||
      !identical(as.integer(map$OptimizerIndex), seq_len(free_dimension))) {
    stop("Internal MML Person-pattern parameter map is incomplete.",
         call. = FALSE)
  }
  map
}

audit_mfrm_mml_observed_pattern_score <- function(
    opt,
    prep,
    idx,
    config,
    sizes,
    quad_points,
    nonlinear_blocks,
    max_persons = 200L,
    max_numeric_persons = 100L,
    max_free_dimension = 80L,
    max_score_elements = 8000,
    relative_step = 1e-6,
    tolerances = c(1e-12, 1e-10, 1e-8)) {
  nonlinear_blocks <- as.character(nonlinear_blocks %||% character(0))
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  person_rows <- length(unique(as.integer(idx$person)))
  weights <- as.numeric(idx$weight %||% rep(1, length(idx$score_k)))
  unit_weights <- length(weights) == length(idx$score_k) &&
    all(is.finite(weights)) && all(abs(weights - 1) <= 1e-12)
  severity <- as.character(
    opt$optimizer_diagnostics$ConvergenceSeverity %||% NA_character_
  )
  base <- list(
    role = "observed_person_log_marginal_pattern_score_jacobian",
    evaluation_point = "retained_optimizer_free_coordinate_vector",
    method = as.character(config$method),
    model = as.character(config$model),
    nonlinear_blocks = nonlinear_blocks,
    optimizer_convergence_severity = severity,
    attempted = FALSE,
    status = if (identical(config$method, "MML") &&
                 length(nonlinear_blocks) > 0L) {
      "pending_after_optimization"
    } else if (identical(config$method, "MML")) {
      "not_required_linear_preflight_scope"
    } else {
      "not_applicable_estimator"
    },
    person_rows = as.integer(person_rows),
    free_dimension = as.integer(free_dimension),
    score_elements = as.double(person_rows) * as.double(free_dimension),
    execution_limits = list(
      persons = as.integer(max_persons),
      numeric_persons = as.integer(max_numeric_persons),
      free_dimension = as.integer(max_free_dimension),
      score_elements = as.double(max_score_elements)
    ),
    quadrature_points = as.integer(quad_points),
    unit_row_weights_observed = isTRUE(unit_weights),
    finite_parameter_vector_observed = FALSE,
    observed_patterns_only = TRUE,
    all_possible_response_patterns_evaluated = FALSE,
    conditional_jml_kernel_reused = FALSE,
    score_rank = NA_integer_,
    score_nullity = NA_integer_,
    score_rank_state = "not_evaluated",
    tolerance_sensitive = FALSE,
    rank_ladder = data.frame(),
    parameter_blocks = data.frame(
      Block = character(0), FreeCoordinates = integer(0)
    ),
    parameter_map = mfrmr_estimability_map(),
    zero_coordinates = character(0),
    null_directions = data.frame(),
    score_row_norm_summary = data.frame(
      Minimum = NA_real_, Median = NA_real_, Maximum = NA_real_
    ),
    reconstruction = data.frame(
      NegativeLogLikelihoodDifference = NA_real_,
      ScoreSumVsNegativeGradientMaxAbs = NA_real_,
      stringsAsFactors = FALSE
    ),
    numerical_differentiation = list(
      status = "not_evaluated",
      method = "coordinate_scaled_central_difference_of_person_log_marginals",
      relative_step = as.numeric(relative_step[1]),
      max_abs_difference = NA_real_,
      max_scaled_difference = NA_real_
    ),
    structural_identification_classified = FALSE,
    weak_information_classified = FALSE,
    readiness_effect = "none_observed_patterns_diagnostic_only",
    detail = paste(
      "No nonlinear MML block required the observed Person-pattern score",
      "audit in this implementation slice."
    )
  )
  if (!identical(config$method, "MML") ||
      length(nonlinear_blocks) == 0L) return(base)

  if (person_rows <= 0L || free_dimension <= 0L ||
      person_rows > as.integer(max_persons) ||
      free_dimension > as.integer(max_free_dimension) ||
      as.double(person_rows) * as.double(free_dimension) >
        as.double(max_score_elements)) {
    base$status <- "not_evaluated_execution_limit"
    base$detail <- paste(
      "The observed Person-pattern score audit exceeded a recorded bounded",
      "execution limit. This is a computational state, not an estimability",
      "or readiness classification."
    )
    return(base)
  }
  if (is.null(opt$par) || length(opt$par) != free_dimension ||
      any(!is.finite(opt$par))) {
    base$status <- "not_evaluated_parameter_vector"
    base$detail <- paste(
      "The retained MML free parameter vector was unavailable, non-finite,",
      "or dimensionally inconsistent."
    )
    return(base)
  }

  quad <- gauss_hermite_normal(max(1L, as.integer(quad_points)))
  built <- tryCatch(
    mfrmr_mml_observed_person_score_matrix(
      par = opt$par,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad
    ),
    error = function(e) e
  )
  base$attempted <- TRUE
  if (inherits(built, "error")) {
    base$status <- "unavailable"
    base$detail <- paste0(
      "The observed Person-pattern score matrix was unavailable: ",
      conditionMessage(built), "."
    )
    return(base)
  }

  parameter_map <- tryCatch(
    mfrmr_mml_optimizer_parameter_map(prep, idx, config, sizes),
    error = function(e) e
  )
  if (inherits(parameter_map, "error")) {
    base$status <- "unavailable"
    base$detail <- paste0(
      "The observed Person-pattern parameter map was unavailable: ",
      conditionMessage(parameter_map), "."
    )
    return(base)
  }
  rank <- suppressWarnings(
    mfrmr_estimability_rank_audit(
      design = Matrix::Matrix(built$score, sparse = TRUE),
      parameter_map = parameter_map,
      tolerances = tolerances,
      structural_tolerance = tolerances[1]
    )
  )
  full_objective <- mfrm_loglik_mml(opt$par, idx, config, sizes, quad)
  full_negative_gradient <- mfrm_grad_mml(
    opt$par, idx, config, sizes, quad
  )
  objective_difference <- -sum(built$log_marginal) - full_objective
  gradient_difference <- colSums(built$score) + full_negative_gradient
  row_norm <- sqrt(rowSums(built$score^2))

  base$status <- "evaluated_observed_pattern_diagnostic_only"
  base$finite_parameter_vector_observed <- TRUE
  base$score_rank <- rank$Rank
  base$score_nullity <- rank$Nullity
  base$score_rank_state <- if (rank$Nullity == 0L) {
    "observed_pattern_full_column_rank"
  } else {
    "observed_pattern_rank_deficient"
  }
  base$tolerance_sensitive <- isTRUE(rank$ToleranceSensitive)
  base$rank_ladder <- rank$ToleranceRanks
  base$parameter_map <- parameter_map
  base$zero_coordinates <- rank$ZeroCoordinates
  base$null_directions <- rank$NullDirections
  base$parameter_blocks <- stats::aggregate(
    rep(1L, nrow(parameter_map)),
    by = list(Block = parameter_map$Block),
    FUN = sum
  ) |>
    stats::setNames(c("Block", "FreeCoordinates"))
  base$score_row_norm_summary <- data.frame(
    Minimum = min(row_norm),
    Median = stats::median(row_norm),
    Maximum = max(row_norm),
    stringsAsFactors = FALSE
  )
  base$reconstruction <- data.frame(
    NegativeLogLikelihoodDifference = objective_difference,
    ScoreSumVsNegativeGradientMaxAbs = max(abs(gradient_difference)),
    stringsAsFactors = FALSE
  )

  if (person_rows <= as.integer(max_numeric_persons)) {
    numeric <- mfrmr_numeric_transformation_jacobian(
      transform = function(value) {
        mfrmr_mml_person_log_marginals(
          par = value,
          idx = idx,
          config = config,
          sizes = sizes,
          quad = quad
        )$log_marginal
      },
      free = opt$par,
      relative_step = relative_step
    )
    if (isTRUE(numeric$valid) &&
        identical(dim(numeric$jacobian), dim(built$score))) {
      difference <- abs(built$score - numeric$jacobian)
      scaled <- difference /
        pmax(1, abs(built$score), abs(numeric$jacobian))
      base$numerical_differentiation$status <- "evaluated"
      base$numerical_differentiation$max_abs_difference <- max(difference)
      base$numerical_differentiation$max_scaled_difference <- max(scaled)
    } else {
      base$numerical_differentiation$status <- "unavailable"
    }
  } else {
    base$numerical_differentiation$status <-
      "not_evaluated_execution_limit"
  }
  base$detail <- paste(
    "Rows are conventional log-likelihood scores for the observed marginal",
    "Person response patterns. Their sum reconstructs the full MML score and",
    "eligible bounded cases are checked against a central difference of that",
    "Person's log marginal contribution. The separate local-classification",
    "layer may use full column rank of unit-weight observed-pattern scores as",
    "a sufficient subset-of-support certificate for positive fixed-quadrature",
    "Fisher information. Rank deficiency is inconclusive without the all-",
    "pattern audit, and no result establishes global identification, weak",
    "information, or readiness."
  )
  base
}

audit_mfrm_mml_all_pattern_information <- function(
    opt,
    prep,
    idx,
    config,
    sizes,
    quad_points,
    nonlinear_blocks,
    max_persons = 100L,
    max_patterns_per_person = 4096,
    max_total_patterns = 5000,
    max_free_dimension = 80L,
    max_score_elements = 4e5,
    relative_step = 1e-6,
    tolerances = c(1e-12, 1e-10, 1e-8)) {
  nonlinear_blocks <- as.character(nonlinear_blocks %||% character(0))
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  n_cat <- as.integer(config$n_cat %||% 0L)
  workload <- mfrmr_mml_all_pattern_workload(idx, max(n_cat, 2L))
  design_groups <- tryCatch(
    mfrmr_mml_person_design_groups(
      idx = idx,
      config = config,
      reuse_identical_designs = TRUE
    ),
    error = function(e) e
  )
  person_designs <- length(workload$person_ids)
  total_patterns <- workload$total_patterns
  unique_person_designs <- if (inherits(design_groups, "error")) {
    NA_integer_
  } else {
    design_groups$unique_person_designs
  }
  evaluated_patterns <- if (inherits(design_groups, "error")) {
    Inf
  } else {
    sum(workload$pattern_counts[design_groups$representative_positions])
  }
  score_elements <- evaluated_patterns * as.double(free_dimension)
  weights <- idx$weight
  unit_weights <- is.null(weights) ||
    (all(is.finite(weights)) && all(abs(weights - 1) <= 1e-12))
  severity <- as.character(
    opt$optimizer_diagnostics$ConvergenceSeverity %||% NA_character_
  )
  empty_map <- mfrmr_estimability_map()
  base <- list(
    role = "all_response_patterns_expected_marginal_score_information",
    evaluation_point = "retained_optimizer_free_coordinate_vector",
    method = as.character(config$method),
    model = as.character(config$model),
    nonlinear_blocks = nonlinear_blocks,
    optimizer_convergence_severity = severity,
    attempted = FALSE,
    status = if (identical(config$method, "MML") &&
                 length(nonlinear_blocks) > 0L) {
      "pending_after_optimization"
    } else if (identical(config$method, "MML")) {
      "not_required_linear_preflight_scope"
    } else {
      "not_applicable_estimator"
    },
    person_designs = as.integer(person_designs),
    retained_observation_rows = length(idx$score_k),
    retained_observations_per_person = data.frame(
      Minimum = if (length(workload$observation_counts) > 0L) {
        min(workload$observation_counts)
      } else {
        NA_integer_
      },
      Median = if (length(workload$observation_counts) > 0L) {
        stats::median(workload$observation_counts)
      } else {
        NA_real_
      },
      Maximum = if (length(workload$observation_counts) > 0L) {
        max(workload$observation_counts)
      } else {
        NA_integer_
      },
      stringsAsFactors = FALSE
    ),
    categories = n_cat,
    total_response_patterns = as.double(total_patterns),
    evaluated_response_patterns = as.double(evaluated_patterns),
    unique_person_designs = as.integer(unique_person_designs),
    design_reuse = data.frame(
      PersonDesigns = as.integer(person_designs),
      UniquePersonDesigns = as.integer(unique_person_designs),
      ReusedPersonDesigns = if (is.finite(unique_person_designs)) {
        as.integer(person_designs - unique_person_designs)
      } else {
        NA_integer_
      },
      LargestDesignGroup = if (!inherits(design_groups, "error") &&
                               length(design_groups$group_sizes) > 0L) {
        max(design_groups$group_sizes)
      } else {
        NA_integer_
      },
      ConceptualToEvaluatedPatternRatio = if (is.finite(evaluated_patterns) &&
                                               evaluated_patterns > 0) {
        total_patterns / evaluated_patterns
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    ),
    free_dimension = as.integer(free_dimension),
    score_elements = as.double(score_elements),
    execution_limits = list(
      persons = as.integer(max_persons),
      patterns_per_person = as.double(max_patterns_per_person),
      evaluated_patterns = as.double(max_total_patterns),
      free_dimension = as.integer(max_free_dimension),
      score_elements = as.double(max_score_elements)
    ),
    quadrature_points = as.integer(quad_points),
    unit_row_weights_required = TRUE,
    unit_row_weights_observed = isTRUE(unit_weights),
    retained_observation_designs = TRUE,
    missing_rows_imputed = FALSE,
    observed_patterns_only = FALSE,
    all_possible_response_patterns_evaluated = FALSE,
    expected_information_evaluated = FALSE,
    conditional_jml_kernel_reused = FALSE,
    local_rank = NA_integer_,
    local_nullity = NA_integer_,
    local_rank_state = "not_evaluated",
    tolerance_sensitive = FALSE,
    rank_ladder = data.frame(),
    parameter_blocks = data.frame(
      Block = character(0), FreeCoordinates = integer(0)
    ),
    parameter_map = empty_map,
    zero_coordinates = character(0),
    null_directions = data.frame(),
    probability_normalization = data.frame(
      MinimumMass = NA_real_,
      MedianMass = NA_real_,
      MaximumMass = NA_real_,
      MaximumAbsoluteMassError = NA_real_,
      stringsAsFactors = FALSE
    ),
    expected_score_identity = data.frame(
      MaximumPersonScoreAbs = NA_real_,
      MaximumAggregateScoreAbs = NA_real_,
      stringsAsFactors = FALSE
    ),
    expected_information_summary = data.frame(
      MinimumEigenvalue = NA_real_,
      MaximumEigenvalue = NA_real_,
      MaximumAsymmetry = NA_real_,
      Trace = NA_real_,
      stringsAsFactors = FALSE
    ),
    numerical_differentiation = list(
      status = "not_evaluated",
      method = "coordinate_scaled_central_difference_of_selected_all_pattern_log_marginals",
      relative_step = as.numeric(relative_step[1]),
      selected_person_designs = 0L,
      selected_patterns = 0L,
      max_abs_difference = NA_real_,
      max_scaled_difference = NA_real_
    ),
    structural_identification_classified = FALSE,
    weak_information_classified = FALSE,
    readiness_effect = "none_all_patterns_local_diagnostic_only",
    detail = paste(
      "No nonlinear MML block required the all-response-pattern expected",
      "information audit in this implementation slice."
    )
  )
  if (!identical(config$method, "MML") ||
      length(nonlinear_blocks) == 0L) return(base)

  if (!isTRUE(unit_weights)) {
    base$status <- "not_evaluated_nonunit_row_weights"
    base$detail <- paste(
      "The exhaustive response-pattern probability identity currently",
      "requires unit row weights. Nonunit likelihood weights do not define",
      "the same normalized finite response-pattern distribution and are",
      "therefore not silently reinterpreted."
    )
    return(base)
  }
  if (inherits(design_groups, "error")) {
    base$status <- "not_evaluated_design_grouping_unavailable"
    base$detail <- paste0(
      "The all-response-pattern Person-design grouping was unavailable: ",
      conditionMessage(design_groups), "."
    )
    return(base)
  }
  maximum_person_patterns <- if (length(workload$pattern_counts) > 0L) {
    max(workload$pattern_counts)
  } else {
    Inf
  }
  if (person_designs <= 0L || free_dimension <= 0L || n_cat < 2L ||
      person_designs > as.integer(max_persons) ||
      maximum_person_patterns > as.double(max_patterns_per_person) ||
      evaluated_patterns > as.double(max_total_patterns) ||
      free_dimension > as.integer(max_free_dimension) ||
      score_elements > as.double(max_score_elements)) {
    base$status <- "not_evaluated_execution_limit"
    base$detail <- paste(
      "The all-response-pattern expected-information audit exceeded a",
      "recorded evaluated-pattern or dimensional execution limit after exact",
      "reuse of identical Person designs. This is a",
      "computational state, not an estimability or readiness classification."
    )
    return(base)
  }
  if (is.null(opt$par) || length(opt$par) != free_dimension ||
      any(!is.finite(opt$par))) {
    base$status <- "not_evaluated_parameter_vector"
    base$detail <- paste(
      "The retained MML free parameter vector was unavailable, non-finite,",
      "or dimensionally inconsistent."
    )
    return(base)
  }

  quad <- gauss_hermite_normal(max(1L, as.integer(quad_points)))
  built <- tryCatch(
    mfrmr_mml_all_pattern_expected_information(
      par = opt$par,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad,
      reuse_identical_designs = TRUE
    ),
    error = function(e) e
  )
  base$attempted <- TRUE
  if (inherits(built, "error")) {
    base$status <- "unavailable"
    base$detail <- paste0(
      "The all-response-pattern expected information was unavailable: ",
      conditionMessage(built), "."
    )
    return(base)
  }
  parameter_map <- tryCatch(
    mfrmr_mml_optimizer_parameter_map(prep, idx, config, sizes),
    error = function(e) e
  )
  if (inherits(parameter_map, "error")) {
    base$status <- "unavailable"
    base$detail <- paste0(
      "The all-response-pattern parameter map was unavailable: ",
      conditionMessage(parameter_map), "."
    )
    return(base)
  }

  rank <- suppressWarnings(
    mfrmr_estimability_rank_audit(
      design = Matrix::Matrix(built$weighted_score, sparse = TRUE),
      parameter_map = parameter_map,
      tolerances = tolerances,
      structural_tolerance = tolerances[1]
    )
  )
  information <- built$expected_information
  information_symmetric <- (information + t(information)) / 2
  information_eigenvalues <- eigen(
    information_symmetric, symmetric = TRUE, only.values = TRUE
  )$values
  expected_aggregate <- colSums(built$expected_score)
  probability_error <- abs(built$probability_mass - 1)

  base$status <- "evaluated_all_patterns_local_diagnostic_only"
  base$all_possible_response_patterns_evaluated <- TRUE
  base$expected_information_evaluated <- TRUE
  base$evaluated_response_patterns <- built$evaluated_patterns
  base$unique_person_designs <- built$unique_person_designs
  base$local_rank <- rank$Rank
  base$local_nullity <- rank$Nullity
  base$local_rank_state <- if (rank$Nullity == 0L) {
    "all_patterns_locally_full_column_rank"
  } else {
    "all_patterns_locally_rank_deficient"
  }
  base$tolerance_sensitive <- isTRUE(rank$ToleranceSensitive)
  base$rank_ladder <- rank$ToleranceRanks
  base$parameter_map <- parameter_map
  base$zero_coordinates <- rank$ZeroCoordinates
  base$null_directions <- rank$NullDirections
  base$parameter_blocks <- stats::aggregate(
    rep(1L, nrow(parameter_map)),
    by = list(Block = parameter_map$Block),
    FUN = sum
  ) |>
    stats::setNames(c("Block", "FreeCoordinates"))
  base$probability_normalization <- data.frame(
    MinimumMass = min(built$probability_mass),
    MedianMass = stats::median(built$probability_mass),
    MaximumMass = max(built$probability_mass),
    MaximumAbsoluteMassError = max(probability_error),
    stringsAsFactors = FALSE
  )
  base$expected_score_identity <- data.frame(
    MaximumPersonScoreAbs = max(abs(built$expected_score)),
    MaximumAggregateScoreAbs = max(abs(expected_aggregate)),
    stringsAsFactors = FALSE
  )
  base$expected_information_summary <- data.frame(
    MinimumEigenvalue = min(information_eigenvalues),
    MaximumEigenvalue = max(information_eigenvalues),
    MaximumAsymmetry = max(abs(information - t(information))),
    Trace = sum(diag(information)),
    stringsAsFactors = FALSE
  )

  selected_positions <- unique(c(1L, length(built$person_ids)))
  selected <- lapply(selected_positions, function(person_position) {
    evaluated <- built$by_person[[person_position]]
    pattern_rows <- unique(c(
      1L,
      as.integer(ceiling(nrow(evaluated$patterns) / 2)),
      nrow(evaluated$patterns)
    ))
    person_idx <- mfrmr_subset_observation_indices(
      idx, built$canonical_rows[[person_position]]
    )
    list(
      person_idx = person_idx,
      patterns = evaluated$patterns[pattern_rows, , drop = FALSE],
      score = evaluated$score[pattern_rows, , drop = FALSE]
    )
  })
  analytic_selected <- do.call(rbind, lapply(selected, `[[`, "score"))
  numeric <- mfrmr_numeric_transformation_jacobian(
    transform = function(value) {
      unlist(lapply(selected, function(selection) {
        mfrmr_mml_evaluate_person_patterns(
          par = value,
          person_idx = selection$person_idx,
          config = config,
          sizes = sizes,
          quad = quad,
          patterns = selection$patterns,
          include_scores = FALSE
        )$log_marginal
      }), use.names = FALSE)
    },
    free = opt$par,
    relative_step = relative_step
  )
  base$numerical_differentiation$selected_person_designs <-
    length(selected)
  base$numerical_differentiation$selected_patterns <-
    nrow(analytic_selected)
  if (isTRUE(numeric$valid) &&
      identical(dim(numeric$jacobian), dim(analytic_selected))) {
    difference <- abs(analytic_selected - numeric$jacobian)
    scaled <- difference /
      pmax(1, abs(analytic_selected), abs(numeric$jacobian))
    base$numerical_differentiation$status <- "evaluated"
    base$numerical_differentiation$max_abs_difference <- max(difference)
    base$numerical_differentiation$max_scaled_difference <- max(scaled)
  } else {
    base$numerical_differentiation$status <- "unavailable"
  }
  base$detail <- paste(
    "For each observed Person design after row omission, every finite",
    "category-response pattern is enumerated under unit weights. Marginal",
    "probability mass, the zero expected-score identity, the score-outer-",
    "product expected information, and selected independent derivatives are",
    "checked in optimizer coordinates. Mathematically identical Person",
    "designs, including any latent-regression design row, are evaluated once",
    "and reconstructed by exact group multiplicity. This retained-point, observed-design",
    "geometry is stronger than an observed-pattern rank but remains a local",
    "diagnostic; it does not yet classify structural identification, weak",
    "information, readiness, or external-comparison eligibility."
  )
  base
}

mfrmr_nonlinear_local_estimability_contract_version <- function() {
  "mfrmr-nonlinear-local-estimability-0.2.3-v1"
}

mfrmr_nonlinear_local_map_complete <- function(source) {
  free_dimension <- as.integer(source$free_dimension %||% NA_integer_)
  parameter_map <- as.data.frame(
    source$parameter_map %||% data.frame(), stringsAsFactors = FALSE
  )
  is.finite(free_dimension) && free_dimension > 0L &&
    nrow(parameter_map) == free_dimension &&
    "OptimizerIndex" %in% names(parameter_map) &&
    identical(as.integer(parameter_map$OptimizerIndex),
              seq_len(free_dimension))
}

mfrmr_classify_nonlinear_local_estimability <- function(audit, config) {
  method <- as.character(config$method %||% "unknown")
  model <- as.character(config$model %||% "unknown")
  nonlinear_blocks <- as.character(
    audit$nonlinear_blocks %||% character(0)
  )
  base <- list(
    contract_version =
      mfrmr_nonlinear_local_estimability_contract_version(),
    method = method,
    model = model,
    nonlinear_blocks = nonlinear_blocks,
    state = if (length(nonlinear_blocks) == 0L) {
      "not_required"
    } else {
      "not_evaluated"
    },
    evidence_basis = if (length(nonlinear_blocks) == 0L) {
      "linear_preflight_complete"
    } else {
      "none"
    },
    probability_model_scope = if (identical(method, "JML")) {
      "conditional_response_given_retained_person_coordinates"
    } else if (identical(method, "MML")) {
      "implemented_fixed_quadrature_marginal_model"
    } else {
      "unsupported_estimator"
    },
    local_rank = NA_integer_,
    local_nullity = NA_integer_,
    free_dimension = NA_integer_,
    tolerance_sensitive = NA,
    parameter_map_complete = FALSE,
    local_first_order_classified = FALSE,
    local_full_rank_sufficient = FALSE,
    local_first_order_rank_deficient = FALSE,
    local_nonidentifiability_established = FALSE,
    global_identification_classified = FALSE,
    continuous_integral_identification_classified = FALSE,
    weak_information_classified = FALSE,
    boundary_classified = FALSE,
    numerical_derivative_status = "not_evaluated",
    quadrature_points = NA_integer_,
    reason_codes = if (length(nonlinear_blocks) == 0L) {
      "nonlinear_coordinates_not_present"
    } else {
      "local_nonlinear_map_not_evaluated"
    },
    readiness_effect = "none_local_property_only",
    detail = paste(
      "No nonlinear local first-order classification was required for this",
      "fit."
    )
  )
  if (length(nonlinear_blocks) == 0L) return(base)

  classify_source <- function(source, basis, rank_field, nullity_field,
                              derivative_status, full_rank_reason,
                              deficient_reason) {
    rank <- as.integer(source[[rank_field]] %||% NA_integer_)
    nullity <- as.integer(source[[nullity_field]] %||% NA_integer_)
    free_dimension <- as.integer(source$free_dimension %||% NA_integer_)
    map_complete <- mfrmr_nonlinear_local_map_complete(source)
    tolerance_sensitive <- isTRUE(source$tolerance_sensitive)
    valid <- is.finite(rank) && is.finite(nullity) &&
      is.finite(free_dimension) && free_dimension > 0L &&
      rank >= 0L && nullity >= 0L &&
      rank + nullity == free_dimension && map_complete
    if (!valid) return(NULL)
    state <- if (tolerance_sensitive) {
      "indeterminate_tolerance_sensitive"
    } else if (nullity == 0L && rank == free_dimension) {
      "locally_full_rank_sufficient"
    } else {
      "locally_first_order_rank_deficient"
    }
    list(
      state = state,
      evidence_basis = basis,
      local_rank = rank,
      local_nullity = nullity,
      free_dimension = free_dimension,
      tolerance_sensitive = tolerance_sensitive,
      parameter_map_complete = map_complete,
      local_first_order_classified = !tolerance_sensitive,
      local_full_rank_sufficient = identical(
        state, "locally_full_rank_sufficient"
      ),
      local_first_order_rank_deficient = identical(
        state, "locally_first_order_rank_deficient"
      ),
      quadrature_points = as.integer(
        source$quadrature_points %||% NA_integer_
      ),
      numerical_derivative_status = derivative_status,
      reason_codes = if (tolerance_sensitive) {
        "local_rank_tolerance_sensitive;global_identification_not_classified"
      } else if (nullity == 0L && rank == free_dimension) {
        paste0(full_rank_reason,
               ";global_identification_not_classified")
      } else {
        paste0(deficient_reason,
               ";local_nonidentifiability_not_established;",
               "global_identification_not_classified")
      }
    )
  }

  classified <- NULL
  if (identical(method, "JML") && identical(model, "GPCM")) {
    source <- audit$gpcm_response_kernel %||% list()
    if (identical(as.character(source$status %||% ""),
                  "evaluated_local_diagnostic_only") &&
        isTRUE(source$conditional_response_kernel_jacobian_evaluated)) {
      classified <- classify_source(
        source = source,
        basis = "conditional_adjacent_logit_jacobian_full_free_coordinates",
        rank_field = "local_rank",
        nullity_field = "local_nullity",
        derivative_status = as.character(
          source$numerical_differentiation$status %||% "not_evaluated"
        ),
        full_rank_reason =
          "conditional_adjacent_logit_map_locally_full_column_rank",
        deficient_reason =
          "conditional_adjacent_logit_map_first_order_rank_deficient"
      )
    }
  } else if (identical(method, "MML")) {
    observed <- audit$mml_observed_pattern_score %||% list()
    observed_eligible <- identical(
      as.character(observed$status %||% ""),
      "evaluated_observed_pattern_diagnostic_only"
    ) && isTRUE(observed$unit_row_weights_observed) &&
      isTRUE(observed$finite_parameter_vector_observed)
    if (observed_eligible) {
      observed_classification <- classify_source(
        source = observed,
        basis = "observed_pattern_subset_score_span_fixed_quadrature",
        rank_field = "score_rank",
        nullity_field = "score_nullity",
        derivative_status = as.character(
          observed$numerical_differentiation$status %||% "not_evaluated"
        ),
        full_rank_reason = paste0(
          "unit_weight_observed_pattern_scores_span_all_free_coordinates_",
          "with_positive_pattern_probability"
        ),
        deficient_reason =
          "observed_pattern_subset_score_rank_deficient"
      )
      if (!is.null(observed_classification) &&
          identical(observed_classification$state,
                    "locally_full_rank_sufficient")) {
        classified <- observed_classification
      }
    }

    # A deficient or tolerance-sensitive observed subset is inconclusive.
    # Exhaustive all-pattern information, when computationally available,
    # provides the complete fixed-quadrature first-order classification.
    if (is.null(classified)) {
      all_pattern <- audit$mml_all_pattern_information %||% list()
      if (identical(as.character(all_pattern$status %||% ""),
                    "evaluated_all_patterns_local_diagnostic_only") &&
          isTRUE(all_pattern$all_possible_response_patterns_evaluated) &&
          isTRUE(all_pattern$expected_information_evaluated) &&
          isTRUE(all_pattern$unit_row_weights_observed)) {
        classified <- classify_source(
          source = all_pattern,
          basis = "all_pattern_expected_score_information_fixed_quadrature",
          rank_field = "local_rank",
          nullity_field = "local_nullity",
          derivative_status = as.character(
            all_pattern$numerical_differentiation$status %||%
              "not_evaluated"
          ),
          full_rank_reason =
            "all_pattern_expected_score_information_locally_full_rank",
          deficient_reason =
            "all_pattern_expected_score_information_first_order_rank_deficient"
        )
      }
    }
    if (is.null(classified) && observed_eligible) {
      classified <- observed_classification
      if (!is.null(classified) &&
          identical(classified$state,
                    "locally_first_order_rank_deficient")) {
        classified$state <- "not_evaluated"
        classified$local_first_order_classified <- FALSE
        classified$local_first_order_rank_deficient <- FALSE
        classified$reason_codes <- paste0(
          "observed_pattern_subset_rank_deficient_is_inconclusive;",
          "all_pattern_information_not_evaluated;",
          "global_identification_not_classified"
        )
      }
    }
  }

  if (is.null(classified)) {
    base$detail <- paste(
      "The estimator-specific nonlinear probability map did not supply a",
      "complete retained-point first-order classification. This may reflect",
      "an execution limit, nonunit weights, an unavailable parameter map, or",
      "an estimator outside the implemented scope."
    )
    return(base)
  }
  for (field in names(classified)) base[[field]] <- classified[[field]]
  base$detail <- if (identical(method, "JML")) {
    paste(
      "Full column rank of the smooth conditional adjacent-logit map in all",
      "optimizer free coordinates is a sufficient retained-point local",
      "first-order certificate. Rank deficiency is only a first-order",
      "diagnosis. Neither result establishes global identification, boundary",
      "interiority, weak information, or inference readiness."
    )
  } else {
    paste(
      "For the implemented fixed-quadrature marginal model, unit-weight",
      "observed-pattern score vectors are a subset of the strictly positive",
      "finite response-pattern support. If they span all optimizer free",
      "coordinates, the corresponding expected score information is positive",
      "definite. A deficient observed subset is inconclusive and is classified",
      "only when exhaustive all-pattern information is available. This does",
      "not establish continuous-integral or global identification, boundary",
      "interiority, weak information, or inference readiness."
    )
  }
  base
}

audit_mfrm_fitted_information <- function(opt,
                                          idx,
                                          config,
                                          sizes,
                                          quad_points,
                                          nonlinear_blocks,
                                          max_free_dimension = 80L,
                                          tolerances = c(1e-10, 1e-8, 1e-6)) {
  nonlinear_blocks <- as.character(nonlinear_blocks %||% character(0))
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  base <- list(
    role = "observed_negative_loglikelihood_hessian_at_retained_solution",
    evaluation_point = "retained_optimizer_free_coordinate_vector",
    method = as.character(config$method),
    model = as.character(config$model),
    nonlinear_blocks = nonlinear_blocks,
    attempted = FALSE,
    status = "not_required",
    free_dimension = as.integer(free_dimension),
    dimension_limit = as.integer(max_free_dimension),
    hessian_elements = as.double(free_dimension)^2,
    numerical_differentiation = list(
      method = "central_difference_of_analytical_gradient_via_stats_optimHess",
      parameter_scale = 1,
      free_coordinate_step = 1e-3
    ),
    evaluation_summary = data.frame(
      RetainedObjective = NA_real_,
      ReevaluatedObjective = NA_real_,
      ObjectiveDifference = NA_real_,
      GradientMaxAbs = NA_real_,
      MaximumHessianAsymmetry = NA_real_,
      stringsAsFactors = FALSE
    ),
    rank_ladder = data.frame(),
    block_summary = mfrmr_fitted_information_block_summary(
      matrix(numeric(0), nrow = 0L, ncol = 0L), sizes, character(0)
    ),
    eigenvalue_summary = data.frame(
      Smallest = NA_real_,
      Largest = NA_real_,
      AbsoluteScale = NA_real_,
      stringsAsFactors = FALSE
    ),
    weak_information_classified = FALSE,
    readiness_effect = "none_pending_pilot_calibrated_rule",
    detail = "No nonlinear free-coordinate block required this fitted-information slice."
  )

  if (length(nonlinear_blocks) == 0L) return(base)

  severity <- as.character(
    opt$optimizer_diagnostics$ConvergenceSeverity %||% NA_character_
  )
  if (!identical(severity, "pass")) {
    base$status <- "not_evaluated_nonstationary"
    base$detail <- paste(
      "Fitted information was not evaluated because the retained optimizer",
      "state did not pass the terminal numerical stationarity gate."
    )
    return(base)
  }
  if (free_dimension <= 0L || is.null(opt$par) ||
      length(opt$par) != free_dimension || any(!is.finite(opt$par))) {
    base$status <- "not_evaluated_parameter_vector"
    base$detail <- paste(
      "Fitted information was not evaluated because the retained free",
      "parameter vector was unavailable, non-finite, or dimensionally inconsistent."
    )
    return(base)
  }
  if (free_dimension > as.integer(max_free_dimension)) {
    base$status <- "not_evaluated_dimension_limit"
    base$detail <- paste0(
      "Fitted information requires a dense ", free_dimension, " x ",
      free_dimension, " Hessian, exceeding the instrumentation limit of ",
      as.integer(max_free_dimension), ". This is an execution limit, not an ",
      "estimability classification."
    )
    return(base)
  }

  quad <- gauss_hermite_normal(max(1L, as.integer(quad_points)))
  cache <- make_param_cache(
    sizes = sizes,
    config = config,
    idx = idx,
    is_mml = identical(config$method, "MML")
  )
  evaluator <- make_mfrm_direct_evaluator(
    method = config$method,
    cache = cache,
    idx = idx,
    config = config,
    sizes = sizes,
    quad = quad
  )
  evaluation <- tryCatch(
    list(
      value = as.numeric(evaluator$value(opt$par)),
      gradient = as.numeric(evaluator$gradient(opt$par))
    ),
    error = function(e) list(value = NA_real_, gradient = numeric(0))
  )
  parameter_scale <- rep(1, free_dimension)
  difference_step <- rep(1e-3, free_dimension)
  hessian <- tryCatch(
    stats::optimHess(
      par = opt$par,
      fn = evaluator$value,
      gr = evaluator$gradient,
      control = list(
        fnscale = 1,
        parscale = parameter_scale,
        ndeps = difference_step
      )
    ),
    error = function(e) e
  )
  base$attempted <- TRUE
  if (inherits(hessian, "error") || is.null(hessian)) {
    base$status <- "unavailable"
    base$detail <- paste0(
      "Fitted-information Hessian was unavailable: ",
      if (inherits(hessian, "error")) conditionMessage(hessian) else "unknown error",
      "."
    )
    return(base)
  }

  hessian <- suppressWarnings(as.matrix(hessian))
  retained_objective <- as.numeric(opt$value %||% NA_real_)
  reevaluated_objective <- as.numeric(evaluation$value[1] %||% NA_real_)
  objective_difference <- if (is.finite(retained_objective) &&
                              is.finite(reevaluated_objective)) {
    reevaluated_objective - retained_objective
  } else {
    NA_real_
  }
  gradient_max_abs <- if (length(evaluation$gradient) > 0L &&
                          all(is.finite(evaluation$gradient))) {
    max(abs(evaluation$gradient))
  } else {
    NA_real_
  }
  maximum_asymmetry <- if (nrow(hessian) == ncol(hessian) &&
                           all(is.finite(hessian))) {
    max(abs(hessian - t(hessian)))
  } else {
    NA_real_
  }
  ladder <- mfrmr_information_rank_ladder(
    hessian = hessian,
    tolerances = tolerances
  )
  if (!isTRUE(ladder$valid)) {
    base$status <- "unavailable"
    base$detail <- paste(
      "Fitted-information Hessian or its symmetric eigenvalue decomposition",
      "was non-finite or unavailable."
    )
    return(base)
  }

  base$status <- "evaluated_diagnostic_only"
  base$evaluation_summary <- data.frame(
    RetainedObjective = retained_objective,
    ReevaluatedObjective = reevaluated_objective,
    ObjectiveDifference = objective_difference,
    GradientMaxAbs = gradient_max_abs,
    MaximumHessianAsymmetry = maximum_asymmetry,
    stringsAsFactors = FALSE
  )
  base$rank_ladder <- ladder$rank_ladder
  base$block_summary <- mfrmr_fitted_information_block_summary(
    hessian = hessian,
    sizes = sizes,
    blocks = nonlinear_blocks
  )
  base$eigenvalue_summary <- data.frame(
    Smallest = ladder$smallest_eigenvalue,
    Largest = ladder$largest_eigenvalue,
    AbsoluteScale = ladder$scale,
    stringsAsFactors = FALSE
  )
  base$detail <- paste(
    "The retained solution was evaluated through a dense numerical Hessian",
    "of the same negative log-likelihood and analytical gradient used by the",
    "direct optimizer, using a recorded free-coordinate difference step.",
    "The tolerance ladder is diagnostic only; it does not",
    "classify weak information or change readiness before pilot calibration."
  )
  base
}

mfrmr_estimability_nonlinear_blocks <- function(sizes) {
  candidates <- intersect(
    names(sizes), c("log_slopes", "log_sigma2")
  )
  candidates[vapply(
    sizes[candidates],
    function(x) as.integer(x %||% 0L) > 0L,
    logical(1)
  )]
}

audit_mfrm_estimability <- function(prep, idx, config, sizes) {
  primary <- mfrmr_estimability_adjacent_design(
    prep = prep, idx = idx, config = config, sizes = sizes
  )
  rank_audit <- mfrmr_estimability_rank_audit(
    primary$design, primary$map
  )

  nonlinear_blocks <- mfrmr_estimability_nonlinear_blocks(sizes)
  complete <- length(nonlinear_blocks) == 0L

  counterfactual <- NULL
  population_assumption_linked <- FALSE
  if (identical(config$method, "MML")) {
    counter_design <- mfrmr_estimability_adjacent_design(
      prep = prep,
      idx = idx,
      config = config,
      sizes = sizes,
      include_person = TRUE,
      include_population_beta = FALSE
    )
    counter_rank <- mfrmr_estimability_rank_audit(
      counter_design$design, counter_design$map
    )
    population_assumption_linked <-
      identical(rank_audit$State, "identified") &&
      identical(counter_rank$State, "structurally_unidentified")
    counterfactual <- list(
      Role = "same fixed-effect structure with free JML Person coordinates",
      Rank = counter_rank$Rank,
      Nullity = counter_rank$Nullity,
      State = counter_rank$State,
      FreeDimension = ncol(counter_design$design),
      ToleranceRanks = counter_rank$ToleranceRanks
    )
  }

  subsets <- calc_subsets(prep$data, c("Person", prep$facet_names))
  components <- nrow(as.data.frame(subsets$summary %||% data.frame()))
  disconnected <- is.finite(components) && components > 1L

  state <- rank_audit$State
  reason_codes <- character(0)
  if (identical(state, "structurally_unidentified")) {
    reason_codes <- c(reason_codes, "design_rank_deficient")
    if (disconnected) {
      reason_codes <- c(reason_codes, "disconnected_without_link")
    }
  } else if (isTRUE(population_assumption_linked)) {
    state <- "population_assumption_linked"
    reason_codes <- c(reason_codes, "population_assumption_linked")
  }
  if (!isTRUE(complete)) {
    reason_codes <- c(reason_codes, "design_rank_not_evaluated")
  }
  reason_codes <- unique(reason_codes)

  optimizer_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  block_summary <- if (nrow(primary$map) == 0L) {
    data.frame(Block = character(0), FreeCoordinates = integer(0))
  } else {
    stats::aggregate(
      rep(1L, nrow(primary$map)),
      by = list(Block = primary$map$Block),
      FUN = sum
    ) |>
      stats::setNames(c("Block", "FreeCoordinates"))
  }
  null_blocks <- if (nrow(rank_audit$NullDirections) == 0L) {
    data.frame(Block = character(0), Appearances = integer(0))
  } else {
    stats::aggregate(
      rep(1L, nrow(rank_audit$NullDirections)),
      by = list(Block = rank_audit$NullDirections$Block),
      FUN = sum
    ) |>
      stats::setNames(c("Block", "Appearances"))
  }

  readiness <- data.frame(
    ReadinessContractVersion = mfrmr_estimability_contract_version(),
    ReadinessScope = "fit",
    EstimabilityState = state,
    ReasonCodes = paste(reason_codes, collapse = ";"),
    Complete = isTRUE(complete),
    AuditedFreeDimension = ncol(primary$design),
    OptimizerFreeDimension = optimizer_dimension,
    Rank = rank_audit$Rank,
    Nullity = rank_audit$Nullity,
    PopulationAssumptionLinked = isTRUE(population_assumption_linked),
    stringsAsFactors = FALSE
  )

  list(
    contract_version = mfrmr_estimability_contract_version(),
    method = config$method,
    model = config$model,
    readiness = readiness,
    complete = isTRUE(complete),
    nonlinear_blocks = nonlinear_blocks,
    design = list(
      role = "adjacent_category_logit_constrained_free_coordinate_design",
      observation_rows = primary$observation_rows,
      transition_rows = primary$transition_rows,
      transitions = primary$transitions,
      free_dimension = ncol(primary$design),
      nonzero_entries = length(primary$design@x),
      rank = rank_audit$Rank,
      nullity = rank_audit$Nullity,
      state = rank_audit$State,
      tolerance_sensitive = isTRUE(rank_audit$ToleranceSensitive),
      tolerance_ranks = rank_audit$ToleranceRanks,
      column_norm_min = rank_audit$ColumnNormMin,
      column_norm_max = rank_audit$ColumnNormMax,
      smallest_singular_value = rank_audit$SmallestSingularValue,
      condition_index = rank_audit$ConditionIndex
    ),
    parameter_blocks = block_summary,
    parameter_map = primary$map,
    zero_coordinates = rank_audit$ZeroCoordinates,
    null_directions = rank_audit$NullDirections,
    null_blocks = null_blocks,
    observed_components = as.integer(components),
    counterfactual_jml = counterfactual,
    population_assumption_linked = isTRUE(population_assumption_linked),
    nonlinear_transformation = list(
      status = if (length(nonlinear_blocks) > 0L) {
        "pending_after_optimization"
      } else {
        "not_required"
      },
      attempted = FALSE,
      nonlinear_blocks = nonlinear_blocks,
      parameterization_only = TRUE,
      likelihood_jacobian_evaluated = FALSE,
      structural_identification_classified = FALSE,
      readiness_effect = "none_parameterization_audit_only"
    ),
    gpcm_response_kernel = list(
      role = "retained_gpcm_adjacent_category_logit_response_kernel_jacobian",
      status = if (identical(config$model, "GPCM")) {
        "pending_after_optimization"
      } else {
        "not_applicable_model"
      },
      attempted = FALSE,
      conditional_response_kernel_jacobian_evaluated = FALSE,
      marginal_person_pattern_jacobian_evaluated = FALSE,
      structural_identification_classified = FALSE,
      readiness_effect = "none_pending_property_and_marginal_model_audits"
    ),
    mml_observed_pattern_score = list(
      role = "observed_person_log_marginal_pattern_score_jacobian",
      status = if (identical(config$method, "MML") &&
                   length(nonlinear_blocks) > 0L) {
        "pending_after_optimization"
      } else if (identical(config$method, "MML")) {
        "not_required_linear_preflight_scope"
      } else {
        "not_applicable_estimator"
      },
      attempted = FALSE,
      observed_patterns_only = TRUE,
      all_possible_response_patterns_evaluated = FALSE,
      structural_identification_classified = FALSE,
      weak_information_classified = FALSE,
      readiness_effect = "none_observed_patterns_diagnostic_only"
    ),
    mml_all_pattern_information = list(
      role = "all_response_patterns_expected_marginal_score_information",
      status = if (identical(config$method, "MML") &&
                   length(nonlinear_blocks) > 0L) {
        "pending_after_optimization"
      } else if (identical(config$method, "MML")) {
        "not_required_linear_preflight_scope"
      } else {
        "not_applicable_estimator"
      },
      attempted = FALSE,
      observed_patterns_only = FALSE,
      all_possible_response_patterns_evaluated = FALSE,
      expected_information_evaluated = FALSE,
      retained_observation_designs = TRUE,
      missing_rows_imputed = FALSE,
      structural_identification_classified = FALSE,
      weak_information_classified = FALSE,
      readiness_effect = "none_all_patterns_local_diagnostic_only"
    ),
    nonlinear_local_estimability = list(
      contract_version =
        mfrmr_nonlinear_local_estimability_contract_version(),
      method = config$method,
      model = config$model,
      nonlinear_blocks = nonlinear_blocks,
      state = if (length(nonlinear_blocks) > 0L) {
        "pending_after_optimization"
      } else {
        "not_required"
      },
      evidence_basis = if (length(nonlinear_blocks) > 0L) {
        "none"
      } else {
        "linear_preflight_complete"
      },
      local_full_rank_sufficient = FALSE,
      global_identification_classified = FALSE,
      continuous_integral_identification_classified = FALSE,
      weak_information_classified = FALSE,
      boundary_classified = FALSE,
      readiness_effect = "none_local_property_only"
    ),
    fitted_information = list(
      status = if (length(nonlinear_blocks) > 0L) {
        "pending_after_optimization"
      } else {
        "pending_weak_information_calibration"
      },
      attempted = FALSE,
      nonlinear_blocks = nonlinear_blocks,
      weak_information_classified = FALSE,
      readiness_effect = "none_pending_pilot_calibrated_rule"
    )
  )
}

mfrmr_estimability_condition <- function(message, audit,
                                          type = c("error", "warning")) {
  type <- match.arg(type)
  subclass <- if (identical(type, "error")) {
    "mfrmr_estimability_error"
  } else {
    "mfrmr_estimability_warning"
  }
  structure(
    list(
      message = as.character(message),
      call = NULL,
      readiness = audit$readiness,
      estimability = audit
    ),
    class = c(
      subclass,
      "mfrmr_readiness_condition",
      type,
      "condition"
    )
  )
}

mfrmr_stop_unidentified <- function(audit) {
  affected <- as.character(audit$null_blocks$Block %||% character(0))
  affected <- affected[!is.na(affected) & nzchar(affected)]
  affected_text <- if (length(affected) > 0L) {
    paste0(" Affected free-coordinate blocks include ",
           paste(affected, collapse = ", "), ".")
  } else {
    ""
  }
  message <- paste0(
    "The estimator-specific constrained design is structurally unidentified ",
    "(rank ", audit$design$rank, " of ", audit$design$free_dimension,
    "; nullity ", audit$design$nullity, "). Optimization was not run.",
    affected_text,
    " Change the observation design or supply defensible identifying anchors."
  )
  stop(mfrmr_estimability_condition(message, audit, type = "error"))
}

mfrmr_warn_estimability <- function(audit) {
  state <- as.character(audit$readiness$EstimabilityState[1])
  message <- if (identical(state, "population_assumption_linked")) {
    paste(
      "The MML fixed-effect design is full rank, but the corresponding",
      "free-Person JML design is rank deficient. Cross-panel identification",
      "therefore relies on the common latent-population assumption; review",
      "assignment and population invariance before reporting facet contrasts."
    )
  } else {
    paste0(
      "The constrained design is full rank only at the strictest audited ",
      "tolerance. Treat this fit as weak-information review rather than as ",
      "an exact nonidentification failure."
    )
  }
  warning(mfrmr_estimability_condition(message, audit, type = "warning"))
  invisible(audit)
}
