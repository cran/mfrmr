# Constrained JML additive recession audits
# ==============================================================================
#
# The Person-fixed audit asks a narrower question than the sufficient-score
# rule: is there a direction in the retained additive structural parameter
# space along which every observed category is at least as well supported as
# every alternative category and at least one comparison improves?  Its joint
# companion lets retained Person and structural coordinates move together and
# screens the global cone before target enumeration.  Both reuse the exact
# optimizer expansion, so facet signs, direct/group anchors, sum-zero
# constraints, interactions, and steps enter without being reconstructed from
# response totals.
#
# The result is deliberately an internal certificate instrument.  It does not
# yet overwrite the finite optimizer iterate in public facet, interaction, or
# step tables.  That promotion belongs to the cross-surface propagation work
# after the certificate and target mapping have been stress-tested.

mfrmr_jml_recession_empty_target_table <- function() {
  data.frame(
    ParameterId = character(0),
    ParameterClass = character(0),
    Facet = character(0),
    Level = character(0),
    OptimizerEstimate = numeric(0),
    PositiveRecession = logical(0),
    NegativeRecession = logical(0),
    CandidateStatus = character(0),
    EvaluationState = character(0),
    ReasonCodes = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_recession_empty_certificate_table <- function() {
  data.frame(
    ParameterId = character(0),
    RequestedDirection = character(0),
    SolverStatus = integer(0),
    TargetCapacity = numeric(0),
    TargetChange = numeric(0),
    MinimumContrastMargin = numeric(0),
    PositiveContrastMargin = numeric(0),
    StrictContrastRows = integer(0),
    DirectionL1 = numeric(0),
    Certified = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_recession_empty_loading_table <- function() {
  data.frame(
    ParameterId = character(0),
    RequestedDirection = character(0),
    OptimizerIndex = integer(0),
    Coordinate = character(0),
    Loading = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_observed_contrast_design <- function(adjacent_design,
                                                score_k,
                                                n_obs,
                                                n_steps,
                                                implementation = c(
                                                  "preallocated", "reference"
                                                )) {
  adjacent_design <- methods::as(adjacent_design, "dgCMatrix")
  score_k <- suppressWarnings(as.integer(score_k))
  n_obs <- as.integer(n_obs)
  n_steps <- as.integer(n_steps)
  if (n_obs < 1L || n_steps < 1L || length(score_k) != n_obs ||
      nrow(adjacent_design) != n_obs * n_steps ||
      anyNA(score_k) || any(score_k < 0L | score_k > n_steps)) {
    stop("Internal JML recession audit received a malformed adjacent design.",
         call. = FALSE)
  }
  implementation <- match.arg(implementation)

  rows_per_observation <- n_steps
  n_contrasts <- n_obs * rows_per_observation
  if (identical(implementation, "reference")) {
    row_i <- integer(0)
    col_j <- integer(0)
    values <- numeric(0)
    contrast_cursor <- 0L

    for (observation in seq_len(n_obs)) {
      observed <- score_k[observation]
      alternatives <- setdiff(0:n_steps, observed)
      for (alternative in alternatives) {
        contrast_cursor <- contrast_cursor + 1L
        if (alternative < observed) {
          transitions <- seq.int(alternative + 1L, observed)
          signs <- rep(1, length(transitions))
        } else {
          transitions <- seq.int(observed + 1L, alternative)
          signs <- rep(-1, length(transitions))
        }
        adjacent_rows <- (transitions - 1L) * n_obs + observation
        row_i <- c(row_i, rep.int(contrast_cursor, length(adjacent_rows)))
        col_j <- c(col_j, adjacent_rows)
        values <- c(values, signs)
      }
    }
  } else {
    templates <- lapply(0:n_steps, function(observed) {
      alternatives <- setdiff(0:n_steps, observed)
      widths <- abs(alternatives - observed)
      transitions <- unlist(lapply(alternatives, function(alternative) {
        if (alternative < observed) {
          seq.int(alternative + 1L, observed)
        } else {
          seq.int(observed + 1L, alternative)
        }
      }), use.names = FALSE)
      list(
        contrast_row = rep.int(seq_along(alternatives), widths),
        transition = as.integer(transitions),
        value = rep.int(ifelse(alternatives < observed, 1, -1), widths)
      )
    })
    template_sizes <- vapply(templates, function(template) {
      length(template$transition)
    }, integer(1))
    stored_entries <- sum(as.double(template_sizes[score_k + 1L]))
    if (!is.finite(stored_entries) || stored_entries > .Machine$integer.max) {
      stop("Internal JML recession contrast design exceeded sparse index limits.",
           call. = FALSE)
    }

    row_i <- integer(stored_entries)
    col_j <- integer(stored_entries)
    values <- numeric(stored_entries)
    entry_cursor <- 0L
    for (observation in seq_len(n_obs)) {
      template <- templates[[score_k[observation] + 1L]]
      width <- length(template$transition)
      slots <- entry_cursor + seq_len(width)
      row_i[slots] <-
        (observation - 1L) * rows_per_observation + template$contrast_row
      col_j[slots] <-
        (template$transition - 1L) * n_obs + observation
      values[slots] <- template$value
      entry_cursor <- entry_cursor + width
    }
    if (entry_cursor != stored_entries) {
      stop("Internal JML recession contrast allocation was inconsistent.",
           call. = FALSE)
    }
  }

  contrast_operator <- Matrix::sparseMatrix(
    i = row_i,
    j = col_j,
    x = values,
    dims = c(n_contrasts, nrow(adjacent_design))
  )
  methods::as(contrast_operator %*% adjacent_design, "dgCMatrix")
}

mfrmr_jml_structural_target_system <- function(config, sizes, params,
                                                include_person = FALSE) {
  slices <- build_param_slices(sizes)
  total_parameters <- sum(as.integer(unlist(sizes, use.names = FALSE)))
  target_rows <- list()
  triplet_i <- integer(0)
  triplet_j <- integer(0)
  triplet_x <- numeric(0)
  row_cursor <- 0L

  append_targets <- function(metadata, jacobian, optimizer_slice) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    jacobian <- methods::as(jacobian, "dgCMatrix")
    optimizer_slice <- as.integer(optimizer_slice)
    if (nrow(metadata) != nrow(jacobian) ||
        ncol(jacobian) != length(optimizer_slice)) {
      stop("Internal JML recession target-map dimension mismatch.",
           call. = FALSE)
    }
    if (nrow(metadata) == 0L) return(invisible(NULL))
    target_rows[[length(target_rows) + 1L]] <<- metadata
    nonzero <- Matrix::summary(jacobian)
    if (nrow(nonzero) > 0L) {
      triplet_i <<- c(triplet_i, row_cursor + nonzero$i)
      triplet_j <<- c(triplet_j, optimizer_slice[nonzero$j])
      triplet_x <<- c(triplet_x, nonzero$x)
    }
    row_cursor <<- row_cursor + nrow(metadata)
    invisible(NULL)
  }

  if (isTRUE(include_person)) {
    levels <- as.character(config$theta_spec$levels %||% character(0))
    jacobian <- mfrmr_constraint_jacobian_sparse(
      config$theta_spec, "Person"
    )$jacobian
    append_targets(
      data.frame(
        ParameterId = paste0("Person:", levels),
        ParameterClass = "person",
        Facet = "Person",
        Level = levels,
        OptimizerEstimate = as.numeric(params$theta),
        stringsAsFactors = FALSE
      ),
      jacobian,
      slices$theta %||% integer(0)
    )
  }

  for (facet in config$facet_names) {
    levels <- as.character(config$facet_specs[[facet]]$levels %||%
                             config$facet_levels[[facet]] %||% character(0))
    jacobian <- mfrmr_constraint_jacobian_sparse(
      config$facet_specs[[facet]], facet
    )$jacobian
    append_targets(
      data.frame(
        ParameterId = paste0("Facet:", facet, ":", levels),
        ParameterClass = "facet",
        Facet = facet,
        Level = levels,
        OptimizerEstimate = as.numeric(params$facets[[facet]]),
        stringsAsFactors = FALSE
      ),
      jacobian,
      slices[[facet]] %||% integer(0)
    )
  }

  interaction_specs <- config$interaction_specs %||% list()
  interaction_slice <- as.integer(slices$interactions %||% integer(0))
  interaction_cursor <- 0L
  if (length(interaction_specs) > 0L) {
    for (name in names(interaction_specs)) {
      spec <- interaction_specs[[name]]
      built <- mfrmr_interaction_jacobian_sparse(spec)
      local_n <- ncol(built$jacobian)
      local_slice <- if (local_n > 0L) {
        interaction_slice[interaction_cursor + seq_len(local_n)]
      } else {
        integer(0)
      }
      interaction_cursor <- interaction_cursor + local_n
      grid <- expand.grid(
        FacetA_Level = as.character(spec$levels_a),
        FacetB_Level = as.character(spec$levels_b),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      cell_level <- paste(grid$FacetA_Level, grid$FacetB_Level, sep = "*")
      append_targets(
        data.frame(
          ParameterId = paste0("Interaction:", name, ":", cell_level),
          ParameterClass = "interaction",
          Facet = name,
          Level = cell_level,
          OptimizerEstimate = as.numeric(params$interactions[[name]]),
          stringsAsFactors = FALSE
        ),
        built$jacobian,
        local_slice
      )
    }
  }

  step_built <- mfrmr_step_jacobian_sparse(config, sizes)
  n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
  step_scopes <- if (identical(config$model, "RSM")) {
    "shared"
  } else {
    as.character(config$facet_levels[[config$step_facet]] %||% character(0))
  }
  step_grid <- if (length(step_scopes) > 0L && n_steps > 0L) {
    expand.grid(
      Step = paste0("Step_", seq_len(n_steps)),
      Scope = step_scopes,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )[, c("Scope", "Step"), drop = FALSE]
  } else {
    data.frame(Scope = character(0), Step = character(0))
  }
  step_estimates <- if (identical(config$model, "RSM")) {
    as.numeric(params$steps)
  } else {
    as.numeric(t(params$steps_mat))
  }
  append_targets(
    data.frame(
      ParameterId = paste0("Step:", step_grid$Scope, ":", step_grid$Step),
      ParameterClass = "step",
      Facet = if (identical(config$model, "RSM")) {
        "shared_rating_scale"
      } else {
        config$step_facet
      },
      Level = paste(step_grid$Scope, step_grid$Step, sep = "*"),
      OptimizerEstimate = step_estimates,
      stringsAsFactors = FALSE
    ),
    step_built$jacobian,
    slices$steps %||% integer(0)
  )

  metadata <- if (length(target_rows) == 0L) {
    mfrmr_jml_recession_empty_target_table()[, c(
      "ParameterId", "ParameterClass", "Facet", "Level",
      "OptimizerEstimate"
    )]
  } else {
    do.call(rbind, target_rows)
  }
  rownames(metadata) <- NULL
  expansion <- Matrix::sparseMatrix(
    i = triplet_i,
    j = triplet_j,
    x = triplet_x,
    dims = c(nrow(metadata), total_parameters)
  )
  list(metadata = metadata, expansion = methods::as(expansion, "dgCMatrix"))
}

mfrmr_jml_recession_lp_base <- function(
    contrast_design,
    representation = c("sparse_triplet", "dense_reference")) {
  representation <- match.arg(representation)
  contrast_design <- methods::as(contrast_design, "dgCMatrix")
  n_parameters <- ncol(contrast_design)
  n_contrasts <- nrow(contrast_design)
  constraint_matrix <- NULL
  constraint_triplet <- matrix(
    numeric(0), nrow = 0L, ncol = 3L,
    dimnames = list(NULL, c("constraint", "variable", "value"))
  )
  if (identical(representation, "dense_reference")) {
    dense_contrast <- as.matrix(contrast_design)
    constraint_matrix <- rbind(
      cbind(dense_contrast, -dense_contrast),
      cbind(diag(n_parameters), diag(n_parameters))
    )
  } else {
    nonzero <- Matrix::summary(contrast_design)
    contrast_triplet <- if (nrow(nonzero) > 0L) {
      rbind(
        cbind(nonzero$i, nonzero$j, nonzero$x),
        cbind(nonzero$i, n_parameters + nonzero$j, -nonzero$x)
      )
    } else {
      matrix(numeric(0), nrow = 0L, ncol = 3L)
    }
    box_rows <- n_contrasts + seq_len(n_parameters)
    box_triplet <- rbind(
      cbind(box_rows, seq_len(n_parameters), 1),
      cbind(box_rows, n_parameters + seq_len(n_parameters), 1)
    )
    constraint_triplet <- rbind(contrast_triplet, box_triplet)
    storage.mode(constraint_triplet) <- "double"
    colnames(constraint_triplet) <- c("constraint", "variable", "value")
  }
  list(
    contrast_design = contrast_design,
    representation = representation,
    constraint_matrix = constraint_matrix,
    constraint_triplet = constraint_triplet,
    constraint_direction = c(
      rep(">=", n_contrasts),
      rep("<=", n_parameters)
    ),
    constraint_rhs = c(
      rep(0, n_contrasts),
      rep(1, n_parameters)
    ),
    n_parameters = as.integer(n_parameters),
    n_constraints = as.integer(n_contrasts + n_parameters),
    stored_constraint_nonzeros = if (identical(
      representation, "sparse_triplet"
    )) {
      as.integer(nrow(constraint_triplet))
    } else {
      as.integer(sum(constraint_matrix != 0))
    }
  )
}

mfrmr_jml_recession_shared_geometry <- function(prep, idx, config, sizes,
                                                  params) {
  finish <- function(state, detail, adjacent = NULL,
                     target_system = NULL, contrast = NULL) {
    list(
      contract_version = "mfrmr-jml-recession-shared-geometry-v1",
      state = state,
      ready = identical(state, "ready"),
      adjacent = adjacent,
      target_system = target_system,
      contrast = contrast,
      detail = detail
    )
  }

  adjacent <- tryCatch(
    mfrmr_estimability_adjacent_design(
      prep = prep,
      idx = idx,
      config = config,
      sizes = sizes,
      include_person = TRUE,
      include_population_beta = FALSE
    ),
    error = function(e) e
  )
  if (inherits(adjacent, "error")) {
    return(finish(
      "not_evaluated_design",
      paste0("Adjacent design construction failed: ",
             conditionMessage(adjacent))
    ))
  }

  target_system <- tryCatch(
    mfrmr_jml_structural_target_system(
      config, sizes, params, include_person = TRUE
    ),
    error = function(e) e
  )
  if (inherits(target_system, "error")) {
    return(finish(
      "not_evaluated_mapping",
      paste0("Expanded target mapping failed: ",
             conditionMessage(target_system)),
      adjacent = adjacent
    ))
  }

  contrast <- tryCatch(
    mfrmr_jml_observed_contrast_design(
      adjacent_design = adjacent$design,
      score_k = idx$score_k,
      n_obs = adjacent$observation_rows,
      n_steps = adjacent$transitions
    ),
    error = function(e) e
  )
  if (inherits(contrast, "error")) {
    return(finish(
      "not_evaluated_contrast",
      paste0("Observed-category contrast construction failed: ",
             conditionMessage(contrast)),
      adjacent = adjacent,
      target_system = target_system
    ))
  }

  optimizer_index <- as.integer(adjacent$map$OptimizerIndex)
  mapping_valid <-
    is.data.frame(adjacent$map) &&
    ncol(adjacent$design) == nrow(adjacent$map) &&
    ncol(contrast) == ncol(adjacent$design) &&
    is.data.frame(target_system$metadata) &&
    nrow(target_system$metadata) == nrow(target_system$expansion) &&
    !anyNA(optimizer_index) && !anyDuplicated(optimizer_index) &&
    ncol(target_system$expansion) >= max(c(0L, optimizer_index))
  if (!isTRUE(mapping_valid)) {
    return(finish(
      "not_evaluated_mapping",
      "Shared recession geometry returned inconsistent coordinate mappings.",
      adjacent = adjacent,
      target_system = target_system,
      contrast = contrast
    ))
  }

  finish(
    "ready",
    paste(
      "One full joint adjacent design, target system, and observed contrast",
      "is available for exact structural column projection."
    ),
    adjacent = adjacent,
    target_system = target_system,
    contrast = methods::as(contrast, "dgCMatrix")
  )
}

mfrmr_jml_recession_fit_policy <- function() {
  list(
    contract_version = "mfrmr-jml-recession-fit-policy-v1",
    scope = "additive_structural_and_joint_recession",
    native_timeout_seconds = 10L,
    attempts_per_stage = 1L,
    retry_after_unaccepted = FALSE,
    maximum_native_seconds_per_positive_target = 20L
  )
}

mfrmr_jml_recession_run_lp <- function(lp_base,
                                        objective,
                                        extra_constraint = NULL,
                                        extra_direction = NULL,
                                        extra_rhs = NULL,
                                        timeout = (
                                          mfrmr_jml_recession_fit_policy()$
                                            native_timeout_seconds
                                        )) {
  objective <- as.numeric(objective)
  constraint_direction <- lp_base$constraint_direction
  constraint_rhs <- lp_base$constraint_rhs
  constraint_matrix <- lp_base$constraint_matrix
  constraint_triplet <- lp_base$constraint_triplet

  if (!is.null(extra_constraint)) {
    extra_constraint <- as.numeric(extra_constraint)
    if (length(extra_constraint) != length(objective) ||
        length(extra_direction) != 1L || length(extra_rhs) != 1L) {
      stop("Internal JML recession LP augmentation is malformed.",
           call. = FALSE)
    }
    constraint_direction <- c(
      constraint_direction, as.character(extra_direction)
    )
    constraint_rhs <- c(constraint_rhs, as.numeric(extra_rhs))
    if (identical(lp_base$representation, "dense_reference")) {
      constraint_matrix <- rbind(constraint_matrix, extra_constraint)
    } else {
      keep <- which(extra_constraint != 0)
      if (length(keep) > 0L) {
        extra_triplet <- cbind(
          lp_base$n_constraints + 1L,
          keep,
          extra_constraint[keep]
        )
        storage.mode(extra_triplet) <- "double"
        constraint_triplet <- rbind(constraint_triplet, extra_triplet)
      }
    }
  }

  arguments <- list(
    direction = "max",
    objective.in = objective,
    const.dir = constraint_direction,
    const.rhs = constraint_rhs,
    timeout = as.integer(timeout)
  )
  if (identical(lp_base$representation, "dense_reference")) {
    arguments$const.mat <- constraint_matrix
  } else {
    arguments$dense.const <- constraint_triplet
  }
  do.call(lpSolve::lp, arguments)
}

mfrmr_jml_recession_target_lp <- function(lp_base,
                                           target,
                                           objective_tolerance = 1e-7,
                                           certificate_tolerance = 1e-7,
                                           timeout = (
                                             mfrmr_jml_recession_fit_policy()$
                                               native_timeout_seconds
                                           )) {
  target <- as.numeric(target)
  n_parameters <- length(target)
  signed_target <- c(target, -target)
  capacity_fit <- tryCatch(
    mfrmr_jml_recession_run_lp(
      lp_base = lp_base,
      objective = signed_target,
      timeout = as.integer(timeout)
    ),
    error = function(e) e
  )
  if (inherits(capacity_fit, "error") ||
      !identical(as.integer(capacity_fit$status %||% -1L), 0L)) {
    return(list(
      evaluated = FALSE,
      certified = FALSE,
      solver_status = if (inherits(capacity_fit, "error")) {
        NA_integer_
      } else {
        as.integer(capacity_fit$status)
      },
      target_capacity = NA_real_,
      target_change = NA_real_,
      minimum_margin = NA_real_,
      positive_margin = NA_real_,
      strict_rows = NA_integer_,
      direction = rep(NA_real_, n_parameters),
      lp_calls = 1L,
      reason = "linear_program_capacity_failed"
    ))
  }

  capacity <- as.numeric(capacity_fit$objval)
  if (!is.finite(capacity) || capacity <= 10 * objective_tolerance) {
    return(list(
      evaluated = TRUE,
      certified = FALSE,
      solver_status = as.integer(capacity_fit$status),
      target_capacity = capacity,
      target_change = 0,
      minimum_margin = 0,
      positive_margin = 0,
      strict_rows = 0L,
      direction = rep(0, n_parameters),
      lp_calls = 1L,
      reason = "no_target_recession_direction"
    ))
  }

  target_floor <- max(objective_tolerance * 2, capacity * 1e-5)
  strict_objective <- as.numeric(Matrix::colSums(lp_base$contrast_design))
  strict_fit <- tryCatch(
    mfrmr_jml_recession_run_lp(
      lp_base = lp_base,
      objective = c(strict_objective, -strict_objective),
      extra_constraint = signed_target,
      extra_direction = ">=",
      extra_rhs = target_floor,
      timeout = as.integer(timeout)
    ),
    error = function(e) e
  )
  if (inherits(strict_fit, "error") ||
      !identical(as.integer(strict_fit$status %||% -1L), 0L)) {
    return(list(
      evaluated = FALSE,
      certified = FALSE,
      solver_status = if (inherits(strict_fit, "error")) {
        NA_integer_
      } else {
        as.integer(strict_fit$status)
      },
      target_capacity = capacity,
      target_change = NA_real_,
      minimum_margin = NA_real_,
      positive_margin = NA_real_,
      strict_rows = NA_integer_,
      direction = rep(NA_real_, n_parameters),
      lp_calls = 2L,
      reason = "linear_program_strictness_failed"
    ))
  }

  solution <- as.numeric(strict_fit$solution)
  direction <- solution[seq_len(n_parameters)] -
    solution[n_parameters + seq_len(n_parameters)]
  margins <- as.numeric(lp_base$contrast_design %*% direction)
  target_change <- sum(target * direction)
  minimum_margin <- min(margins)
  positive_margin <- sum(pmax(margins, 0))
  strict_rows <- sum(margins > certificate_tolerance)
  certified <- is.finite(target_change) &&
    target_change > objective_tolerance &&
    is.finite(minimum_margin) &&
    minimum_margin >= -certificate_tolerance &&
    is.finite(positive_margin) &&
    positive_margin > objective_tolerance &&
    strict_rows > 0L

  list(
    evaluated = TRUE,
    certified = isTRUE(certified),
    solver_status = as.integer(strict_fit$status),
    target_capacity = capacity,
    target_change = target_change,
    minimum_margin = minimum_margin,
    positive_margin = positive_margin,
    strict_rows = as.integer(strict_rows),
    direction = direction,
    lp_calls = 2L,
    reason = if (isTRUE(certified)) {
      "certified_additive_recession_direction"
    } else {
      "candidate_failed_postsolve_certificate"
    }
  )
}

mfrmr_jml_recession_target_nullspace_screen <- function(
    contrast_design,
    target_expansion,
    tolerances = c(1e-12, 1e-10, 1e-8),
    max_nonzeros = 2e6,
    max_rows = 100000L,
    max_coordinates = 750L,
    max_target_rows = 200L) {
  contrast_design <- methods::as(contrast_design, "dgCMatrix")
  target_expansion <- methods::as(target_expansion, "dgCMatrix")
  tolerances <- as.numeric(tolerances)
  finish <- function(state, evaluated, safe_to_skip,
                     rank_ladder = data.frame(), detail = "") {
    list(
      contract_version =
        "mfrmr-jml-selected-target-nullspace-screen-v1",
      state = state,
      evaluated = isTRUE(evaluated),
      safe_to_skip = isTRUE(safe_to_skip),
      target_null_movement = if (identical(
        state, "target_changes_quotient_nullspace"
      )) {
        TRUE
      } else if (isTRUE(safe_to_skip)) {
        FALSE
      } else {
        NA
      },
      tolerance_sensitive = if (nrow(rank_ladder) > 0L) {
        length(unique(rank_ladder$ContrastRank)) > 1L ||
          length(unique(rank_ladder$AugmentedRank)) > 1L ||
          length(unique(rank_ladder$RankIncrement)) > 1L
      } else {
        NA
      },
      contrast_rows = as.integer(nrow(contrast_design)),
      target_rows = as.integer(nrow(target_expansion)),
      coordinates = as.integer(ncol(contrast_design)),
      stored_nonzeros = as.double(
        length(contrast_design@x) + length(target_expansion@x)
      ),
      tolerances = tolerances,
      rank_ladder = rank_ladder,
      detail = detail
    )
  }

  if (ncol(contrast_design) != ncol(target_expansion) ||
      length(tolerances) < 1L || any(!is.finite(tolerances)) ||
      any(tolerances <= 0) || any(!is.finite(contrast_design@x)) ||
      any(!is.finite(target_expansion@x))) {
    return(finish(
      "not_evaluated_mapping", FALSE, FALSE,
      detail = "Contrast and target coordinates or tolerances were malformed."
    ))
  }
  free_target <- as.numeric(Matrix::rowSums(abs(target_expansion))) > 0
  target_expansion <- target_expansion[free_target, , drop = FALSE]
  if (nrow(target_expansion) == 0L) {
    return(finish(
      "no_free_selected_targets", TRUE, TRUE,
      detail = "No selected target changes in the retained coordinates."
    ))
  }
  stored_nonzeros <- as.double(
    length(contrast_design@x) + length(target_expansion@x)
  )
  if (stored_nonzeros > as.double(max_nonzeros) ||
      nrow(contrast_design) > as.integer(max_rows) ||
      ncol(contrast_design) > as.integer(max_coordinates) ||
      nrow(target_expansion) > as.integer(max_target_rows)) {
    return(finish(
      "not_evaluated_size_limit", FALSE, FALSE,
      detail = paste0(
        "The common-scale rank comparison exceeded its execution limit (",
        format(stored_nonzeros, scientific = FALSE), " stored nonzeros; ",
        nrow(contrast_design), " contrast rows; ",
        nrow(target_expansion), " target rows; ",
        ncol(contrast_design), " coordinates)."
      )
    ))
  }
  if (ncol(contrast_design) == 0L) {
    return(finish(
      "no_retained_coordinates", TRUE, TRUE,
      detail = "No retained coordinate can change a selected target."
    ))
  }

  combined <- rbind(contrast_design, target_expansion)
  norm_sq <- as.numeric(Matrix::colSums(combined^2))
  norms <- sqrt(pmax(norm_sq, 0))
  if (any(!is.finite(norms))) {
    return(finish(
      "not_evaluated_rank", FALSE, FALSE,
      detail = "Common column scaling was non-finite."
    ))
  }
  inverse_norm <- ifelse(norms == 0, 1, 1 / norms)
  scaling <- Matrix::Diagonal(x = inverse_norm)
  contrast_scaled <- Matrix::drop0(methods::as(
    contrast_design %*% scaling, "dgCMatrix"
  ))
  target_scaled <- Matrix::drop0(methods::as(
    target_expansion %*% scaling, "dgCMatrix"
  ))
  augmented_scaled <- rbind(contrast_scaled, target_scaled)
  rank_one <- function(x, tolerance) {
    if (nrow(x) == 0L || ncol(x) == 0L) return(0L)
    tryCatch(
      suppressWarnings(as.integer(Matrix::rankMatrix(
        x, method = "qr", tol = as.numeric(tolerance)
      ))),
      error = function(e) NA_integer_
    )
  }
  contrast_rank <- vapply(
    tolerances, function(tolerance) rank_one(contrast_scaled, tolerance),
    integer(1)
  )
  augmented_rank <- vapply(
    tolerances, function(tolerance) rank_one(augmented_scaled, tolerance),
    integer(1)
  )
  rank_ladder <- data.frame(
    Tolerance = tolerances,
    ContrastRank = contrast_rank,
    AugmentedRank = augmented_rank,
    RankIncrement = as.integer(augmented_rank - contrast_rank),
    stringsAsFactors = FALSE
  )
  if (anyNA(contrast_rank) || anyNA(augmented_rank) ||
      any(rank_ladder$RankIncrement < 0L)) {
    return(finish(
      "not_evaluated_rank", FALSE, FALSE, rank_ladder,
      "Sparse QR did not return a coherent common-scale rank ladder."
    ))
  }
  tolerance_sensitive <-
    length(unique(contrast_rank)) > 1L ||
    length(unique(augmented_rank)) > 1L ||
    length(unique(rank_ladder$RankIncrement)) > 1L
  if (any(rank_ladder$RankIncrement > 0L)) {
    return(finish(
      "target_changes_quotient_nullspace", TRUE, FALSE, rank_ladder,
      paste0(
        "At least one selected target row increases quotient rank by ",
        max(rank_ladder$RankIncrement), "."
      )
    ))
  }
  if (tolerance_sensitive) {
    return(finish(
      "tolerance_sensitive", TRUE, FALSE, rank_ladder,
      "Rank was target-invariant but changed across guarded tolerances."
    ))
  }
  finish(
    "no_target_change_in_quotient_nullspace", TRUE, TRUE, rank_ladder,
    "Every selected target row is in the quotient contrast row space."
  )
}

audit_mfrm_jml_additive_recession <- function(prep,
                                               idx,
                                               config,
                                               sizes,
                                               params,
                                               coordinate_scope = c(
                                                 "structural_fixed_person",
                                                 "joint_person_structural"
                                               ),
                                               target_person_ids = character(0),
                                               profiled_person_ids = character(0),
                                               shared_geometry = NULL,
                                               screen_global_cone = FALSE,
                                               max_lp_elements = NULL,
                                               max_lp_nonzeros = 2e6,
                                               max_lp_constraints = 100000L,
                                               max_structural_coordinates = 500L,
                                               max_person_coordinates = 500L,
                                               max_additive_coordinates = 750L,
                                               max_target_directions = 200L,
                                               lp_timeout = (
                                                 mfrmr_jml_recession_fit_policy()$
                                                   native_timeout_seconds
                                               ),
                                               cone_objective_tolerance = 1e-10,
                                               cone_certificate_tolerance = 1e-7,
                                               nullspace_tolerances = c(
                                                 1e-12, 1e-10, 1e-8
                                               ),
                                               lp_representation = c(
                                                 "sparse_triplet",
                                                 "dense_reference"
                                               )) {
  coordinate_scope <- match.arg(coordinate_scope)
  lp_representation <- match.arg(lp_representation)
  joint_scope <- identical(coordinate_scope, "joint_person_structural")
  target_person_ids <- unique(as.character(target_person_ids))
  target_person_ids <- target_person_ids[
    !is.na(target_person_ids) & nzchar(target_person_ids)
  ]
  profiled_person_ids <- unique(as.character(profiled_person_ids))
  profiled_person_ids <- profiled_person_ids[
    !is.na(profiled_person_ids) & nzchar(profiled_person_ids)
  ]
  method <- as.character(config$method %||% NA_character_)
  model <- as.character(config$model %||% NA_character_)
  empty_targets <- mfrmr_jml_recession_empty_target_table()
  empty_certificates <- mfrmr_jml_recession_empty_certificate_table()
  empty_loadings <- mfrmr_jml_recession_empty_loading_table()
  prescreen <- list(
    contract_version = "mfrmr-jml-global-cone-prescreen-v1",
    requested = isTRUE(screen_global_cone),
    scope = if (joint_scope) {
      "joint_person_structural"
    } else {
      "structural_fixed_person"
    },
    state = if (isTRUE(screen_global_cone)) "not_reached" else "not_requested",
    evaluated = FALSE,
    cone_certified = NA,
    target_enumeration_skipped = FALSE,
    cone_lp_calls = 0L,
    relevance_lp_calls = 0L,
    target_directions_evaluated = 0L,
    target_lp_calls = 0L,
    total_lp_calls = 0L,
    objective_tolerance = as.numeric(cone_objective_tolerance),
    certificate_tolerance = as.numeric(cone_certificate_tolerance)
  )
  relevance <- list(
    contract_version = "mfrmr-jml-known-person-quotient-prescreen-v1",
    requested = joint_scope && length(profiled_person_ids) > 0L,
    state = if (joint_scope && length(profiled_person_ids) > 0L) {
      "not_reached"
    } else {
      "not_applicable"
    },
    evaluated = FALSE,
    mapping_valid = NA,
    profiled_persons = as.integer(length(profiled_person_ids)),
    removed_observations = 0L,
    removed_contrast_rows = 0L,
    removed_coordinates = 0L,
    boundary_rays_valid = NA,
    quotient_cone_evaluated = FALSE,
    quotient_cone_certified = NA,
    quotient_cone_lp_calls = 0L,
    quotient_cone_capacity = NA_real_,
    nullspace_screen = list(),
    selected_target_exclusion_certified = FALSE,
    detail = ""
  )

  finish <- function(state, complete, detail,
                     target_status = empty_targets,
                     certificates = empty_certificates,
                     loadings = empty_loadings,
                     dimensions = data.frame(),
                     cone_certificate = empty_certificates,
                     cone_loadings = empty_loadings) {
    prescreen_out <- prescreen
    prescreen_out$total_lp_calls <- as.integer(
      prescreen_out$cone_lp_calls + prescreen_out$relevance_lp_calls +
        prescreen_out$target_lp_calls
    )
    list(
      contract_version = mfrmr_boundary_contract_version(),
      method = method,
      model = model,
      scope = if (joint_scope) {
        "JML joint Person-structural additive coordinates"
      } else {
        "JML structural additive coordinates with Person fixed"
      },
      state = state,
      complete = isTRUE(complete),
      target_status = target_status,
      certificates = certificates,
      cone_certificate = cone_certificate,
      cone_direction_loadings = cone_loadings,
      direction_loadings = loadings,
      dimensions = dimensions,
      prescreen = prescreen_out,
      relevance_screen = relevance,
      detail = detail,
      limitations = if (joint_scope) {
        paste(
          "This bounded audit certifies joint additive recession directions",
          "for structural targets and unresolved constraint-coupled extreme",
          "Person targets. Ordinary free extreme Persons remain governed by",
          "the sufficient-score audit. GPCM log-slope recession and public-",
          "table propagation are not part of this implementation slice."
        )
      } else {
        paste(
          "This bounded audit certifies additive structural recession directions",
          "with Person coordinates fixed. GPCM log-slope recession and public-",
          "table propagation are not part of this implementation slice; joint",
          "Person movement is evaluated by the companion additive audit."
        )
      }
    )
  }

  if (!identical(method, "JML")) {
    prescreen$state <- "not_applicable_mml"
    relevance$state <- "not_applicable_mml"
    return(finish(
      "not_applicable_mml", TRUE,
      "Additive fixed-effect recession certification is scoped to JML."
    ))
  }
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    prescreen$state <- "not_evaluated_dependency"
    return(finish(
      "not_evaluated_dependency", FALSE,
      "Package 'lpSolve' is required for the bounded internal certificate."
    ))
  }

  use_shared_geometry <- isTRUE(tryCatch(
    is.list(shared_geometry) &&
      identical(
        shared_geometry$contract_version,
        "mfrmr-jml-recession-shared-geometry-v1"
      ) &&
      isTRUE(shared_geometry$ready) &&
      identical(shared_geometry$state, "ready") &&
      is.list(shared_geometry$adjacent) &&
      inherits(shared_geometry$adjacent$design, "Matrix") &&
      is.data.frame(shared_geometry$adjacent$map) &&
      inherits(shared_geometry$contrast, "Matrix") &&
      ncol(shared_geometry$contrast) ==
        ncol(shared_geometry$adjacent$design) &&
      is.list(shared_geometry$target_system) &&
      is.data.frame(shared_geometry$target_system$metadata) &&
      inherits(shared_geometry$target_system$expansion, "Matrix") &&
      nrow(shared_geometry$target_system$metadata) ==
        nrow(shared_geometry$target_system$expansion),
    error = function(e) FALSE
  ))
  adjacent <- if (use_shared_geometry) {
    shared_geometry$adjacent
  } else {
    tryCatch(
      mfrmr_estimability_adjacent_design(
        prep = prep,
        idx = idx,
        config = config,
        sizes = sizes,
        include_person = TRUE,
        include_population_beta = FALSE
      ),
      error = function(e) e
    )
  }
  if (inherits(adjacent, "error")) {
    prescreen$state <- "not_evaluated_design"
    return(finish(
      "not_evaluated_design", FALSE,
      paste0("Adjacent design construction failed: ", conditionMessage(adjacent))
    ))
  }

  structural_columns <- if (joint_scope) {
    seq_len(ncol(adjacent$design))
  } else {
    which(adjacent$map$Block != "Person")
  }
  structural_map <- adjacent$map[structural_columns, , drop = FALSE]
  structural_optimizer_index <- as.integer(structural_map$OptimizerIndex)
  target_system <- if (use_shared_geometry) {
    shared_geometry$target_system
  } else {
    tryCatch(
      mfrmr_jml_structural_target_system(
        config, sizes, params, include_person = joint_scope
      ),
      error = function(e) e
    )
  }
  if (inherits(target_system, "error")) {
    prescreen$state <- "not_evaluated_mapping"
    return(finish(
      "not_evaluated_mapping", FALSE,
      paste0("Expanded target mapping failed: ", conditionMessage(target_system))
    ))
  }
  targets <- target_system$metadata
  if (joint_scope) {
    target_keep <- targets$ParameterClass != "person" |
      targets$ParameterId %in% target_person_ids
    targets <- targets[target_keep, , drop = FALSE]
    target_system$expansion <- target_system$expansion[
      target_keep, , drop = FALSE
    ]
    rownames(targets) <- NULL
  } else if (use_shared_geometry) {
    target_keep <- targets$ParameterClass != "person"
    targets <- targets[target_keep, , drop = FALSE]
    target_system$expansion <- target_system$expansion[
      target_keep, , drop = FALSE
    ]
    rownames(targets) <- NULL
  }
  if (length(structural_columns) == 0L) {
    prescreen$state <- "not_needed_no_free_coordinates"
    targets$PositiveRecession <- FALSE
    targets$NegativeRecession <- FALSE
    targets$CandidateStatus <- "fixed"
    targets$EvaluationState <- "evaluated"
    targets$ReasonCodes <- "no_free_coordinate"
    dimensions <- data.frame(
      Observations = as.integer(adjacent$observation_rows),
      Categories = as.integer(adjacent$transitions + 1L),
      ContrastRows = as.integer(adjacent$transition_rows),
      StructuralFreeCoordinates = 0L,
      PersonFreeCoordinates = 0L,
      AdditiveFreeCoordinates = 0L,
      ExpandedTargets = as.integer(nrow(targets)),
      FreeTargets = 0L,
      TargetDirections = 0L,
      LPVariables = 0L,
      LPConstraints = 0L,
      SparseLPNonzeros = 0,
      DenseReferenceElements = 0,
      LPRepresentation = lp_representation,
      stringsAsFactors = FALSE
    )
    return(finish(
      if (joint_scope) {
        "no_free_additive_coordinates"
      } else {
        "no_free_structural_coordinates"
      },
      TRUE,
      if (joint_scope) {
        "No additive Person or structural free coordinates were present."
      } else {
        "No additive structural free coordinates were present."
      },
      target_status = targets,
      dimensions = dimensions
    ))
  }
  if (anyNA(structural_optimizer_index) ||
      anyDuplicated(structural_optimizer_index)) {
    prescreen$state <- "not_evaluated_mapping"
    return(finish(
      "not_evaluated_mapping", FALSE,
      "Audited optimizer-coordinate mapping was incomplete or duplicated."
    ))
  }

  contrast <- if (use_shared_geometry) {
    shared_geometry$contrast[, structural_columns, drop = FALSE]
  } else {
    tryCatch(
      mfrmr_jml_observed_contrast_design(
        adjacent_design = adjacent$design[, structural_columns, drop = FALSE],
        score_k = idx$score_k,
        n_obs = adjacent$observation_rows,
        n_steps = adjacent$transitions
      ),
      error = function(e) e
    )
  }
  if (inherits(contrast, "error")) {
    prescreen$state <- "not_evaluated_contrast"
    return(finish(
      "not_evaluated_contrast", FALSE,
      paste0("Observed-category contrast construction failed: ",
             conditionMessage(contrast))
    ))
  }

  target_expansion <- target_system$expansion[
    , structural_optimizer_index, drop = FALSE
  ]
  if (nrow(targets) != nrow(target_expansion)) {
    prescreen$state <- "not_evaluated_mapping"
    return(finish(
      "not_evaluated_mapping", FALSE,
      "Expanded target metadata and Jacobian rows did not align."
    ))
  }

  free_target <- as.numeric(Matrix::rowSums(abs(target_expansion))) > 0
  target_directions <- 2L * sum(free_target)
  additive_coordinates <- ncol(contrast)
  person_coordinates <- sum(structural_map$Block == "Person")
  structural_coordinate_count <- sum(structural_map$Block != "Person")
  lp_constraints <- nrow(contrast) + additive_coordinates + 1L
  lp_variables <- 2L * additive_coordinates
  maximum_target_nonzeros <- if (any(free_target)) {
    max(as.numeric(Matrix::rowSums(
      target_expansion[free_target, , drop = FALSE] != 0
    )))
  } else {
    0
  }
  cone_target_nonzeros <- if (isTRUE(screen_global_cone)) {
    sum(as.numeric(Matrix::colSums(contrast)) != 0)
  } else {
    0
  }
  maximum_extra_nonzeros <- max(maximum_target_nonzeros, cone_target_nonzeros)
  sparse_lp_nonzeros <- 2 * length(contrast@x) +
    2 * additive_coordinates + 2 * maximum_extra_nonzeros
  dense_reference_elements <- as.double(lp_constraints) *
    as.double(lp_variables)
  dimensions <- data.frame(
    Observations = as.integer(adjacent$observation_rows),
    Categories = as.integer(adjacent$transitions + 1L),
    ContrastRows = as.integer(nrow(contrast)),
    StructuralFreeCoordinates = as.integer(structural_coordinate_count),
    PersonFreeCoordinates = as.integer(person_coordinates),
    AdditiveFreeCoordinates = as.integer(additive_coordinates),
    ExpandedTargets = as.integer(nrow(targets)),
    FreeTargets = as.integer(sum(free_target)),
    TargetDirections = as.integer(target_directions),
    LPVariables = as.integer(lp_variables),
    LPConstraints = as.integer(lp_constraints),
    SparseLPNonzeros = as.double(sparse_lp_nonzeros),
    DenseReferenceElements = as.double(dense_reference_elements),
    LPRepresentation = lp_representation,
    stringsAsFactors = FALSE
  )
  legacy_dense_limit_hit <- !is.null(max_lp_elements) &&
    length(max_lp_elements) == 1L && is.finite(max_lp_elements) &&
    dense_reference_elements > as.double(max_lp_elements)
  if (legacy_dense_limit_hit ||
      sparse_lp_nonzeros > as.double(max_lp_nonzeros) ||
      lp_constraints > as.integer(max_lp_constraints) ||
      structural_coordinate_count > as.integer(max_structural_coordinates) ||
      person_coordinates > as.integer(max_person_coordinates) ||
      additive_coordinates > as.integer(max_additive_coordinates) ||
      (!joint_scope &&
       target_directions > as.integer(max_target_directions))) {
    prescreen$state <- "not_evaluated_size_limit"
    targets$PositiveRecession <- NA
    targets$NegativeRecession <- NA
    targets$CandidateStatus <- ifelse(
      free_target, "not_evaluated_size_limit", "fixed"
    )
    targets$EvaluationState <- ifelse(
      free_target, "not_evaluated", "evaluated"
    )
    targets$ReasonCodes <- ifelse(
      free_target, "bounded_lp_size_limit", "no_free_coordinate"
    )
    return(finish(
      "not_evaluated_size_limit", FALSE,
      paste0(
        "The bounded LP audit exceeded its current execution limit (",
        format(sparse_lp_nonzeros, scientific = FALSE),
        " sparse nonzeros; ", lp_constraints, " constraints; ",
        person_coordinates, " Person coordinates; ",
        structural_coordinate_count, " structural coordinates; ",
        additive_coordinates, " additive coordinates; ",
        target_directions, " target directions; dense-reference equivalent ",
        format(dense_reference_elements, scientific = FALSE), " elements)."
      ),
      target_status = targets,
      dimensions = dimensions
    ))
  }

  lp_base <- mfrmr_jml_recession_lp_base(
    contrast,
    representation = lp_representation
  )
  cone_certificate <- empty_certificates
  cone_loadings <- empty_loadings
  cone_certified <- FALSE
  if (isTRUE(screen_global_cone)) {
    cone_solved <- mfrmr_jml_recession_target_lp(
      lp_base = lp_base,
      target = as.numeric(Matrix::colSums(contrast)),
      objective_tolerance = as.numeric(cone_objective_tolerance),
      certificate_tolerance = as.numeric(cone_certificate_tolerance),
      timeout = lp_timeout
    )
    prescreen$evaluated <- isTRUE(cone_solved$evaluated)
    prescreen$cone_certified <- isTRUE(cone_solved$certified)
    prescreen$cone_lp_calls <- as.integer(cone_solved$lp_calls %||% 0L)
    cone_certificate <- data.frame(
      ParameterId = if (joint_scope) {
        "__joint_additive_cone__"
      } else {
        "__structural_additive_cone__"
      },
      RequestedDirection = "likelihood_improving",
      SolverStatus = as.integer(cone_solved$solver_status),
      TargetCapacity = as.numeric(cone_solved$target_capacity),
      TargetChange = as.numeric(cone_solved$target_change),
      MinimumContrastMargin = as.numeric(cone_solved$minimum_margin),
      PositiveContrastMargin = as.numeric(cone_solved$positive_margin),
      StrictContrastRows = as.integer(cone_solved$strict_rows),
      DirectionL1 = sum(abs(as.numeric(cone_solved$direction)), na.rm = TRUE),
      Certified = isTRUE(cone_solved$certified),
      stringsAsFactors = FALSE
    )
    cone_certified <- isTRUE(cone_solved$certified)
    if (cone_certified) {
      keep <- which(abs(cone_solved$direction) > 1e-9)
      if (length(keep) > 0L) {
        cone_loadings <- data.frame(
          ParameterId = if (joint_scope) {
            "__joint_additive_cone__"
          } else {
            "__structural_additive_cone__"
          },
          RequestedDirection = "likelihood_improving",
          OptimizerIndex = structural_optimizer_index[keep],
          Coordinate = structural_map$Coordinate[keep],
          Loading = as.numeric(cone_solved$direction[keep]),
          stringsAsFactors = FALSE
        )
      }
    }
    if (!isTRUE(cone_solved$evaluated)) {
      prescreen$state <- "not_evaluated_solver"
      targets$PositiveRecession <- ifelse(free_target, NA, FALSE)
      targets$NegativeRecession <- ifelse(free_target, NA, FALSE)
      targets$CandidateStatus <- ifelse(
        free_target, "not_evaluated_solver", "fixed"
      )
      targets$EvaluationState <- ifelse(
        free_target, "not_evaluated", "evaluated"
      )
      targets$ReasonCodes <- ifelse(
        free_target, "linear_program_failed", "no_free_coordinate"
      )
      return(finish(
        "not_evaluated_solver", FALSE,
        "The global joint additive recession cone could not be evaluated.",
        target_status = targets,
        dimensions = dimensions,
        cone_certificate = cone_certificate,
        cone_loadings = cone_loadings
      ))
    }
    if (!cone_certified) {
      prescreen$state <- "negative_cone_enumeration_skipped"
      prescreen$target_enumeration_skipped <- TRUE
      targets$PositiveRecession <- FALSE
      targets$NegativeRecession <- FALSE
      targets$CandidateStatus <- ifelse(
        free_target, "finite_in_audited_subspace", "fixed"
      )
      targets$EvaluationState <- "evaluated"
      targets$ReasonCodes <- ifelse(
        free_target,
        "no_additive_recession_direction_certified",
        "no_free_coordinate"
      )
      return(finish(
        "none_certified", !identical(model, "GPCM"),
        if (joint_scope) {
          paste0(
            "No joint additive recession direction was certified for ",
            sum(free_target), " selected free expanded targets."
          )
        } else {
          paste0(
            "No additive structural recession direction was certified for ",
            sum(free_target),
            " free expanded targets with Person coordinates fixed."
          )
        },
        target_status = targets,
        dimensions = dimensions,
        cone_certificate = cone_certificate,
        cone_loadings = cone_loadings
      ))
    }
    prescreen$state <- "positive_cone_target_enumeration"
    if (isTRUE(relevance$requested)) {
      profiled_column <- match(
        profiled_person_ids, as.character(structural_map$Coordinate)
      )
      observation_person_id <- paste0(
        "Person:", prep$levels$Person[idx$person]
      )
      profiled_observation <- observation_person_id %in% profiled_person_ids
      target_uses_profiled_coordinate <- if (
        length(profiled_column) > 0L && !anyNA(profiled_column)
      ) {
        length(target_expansion[, profiled_column, drop = FALSE]@x) > 0L
      } else {
        TRUE
      }
      relevance$mapping_valid <-
        length(profiled_column) == length(profiled_person_ids) &&
        !anyNA(profiled_column) && !anyDuplicated(profiled_column) &&
        all(structural_map$Block[profiled_column] == "Person") &&
        setequal(
          as.character(structural_map$Coordinate[profiled_column]),
          profiled_person_ids
        ) &&
        all(profiled_person_ids %in% observation_person_id) &&
        !any(profiled_person_ids %in% targets$ParameterId) &&
        !target_uses_profiled_coordinate
      if (!isTRUE(relevance$mapping_valid)) {
        relevance$state <- "not_evaluated_mapping_fallback"
        relevance$detail <- paste(
          "Known free extreme Person rows, coordinates, and selected-target",
          "mapping did not support quotient profiling; target enumeration",
          "was retained."
        )
      } else {
        boundary_ray_valid <- vapply(
          seq_along(profiled_person_ids),
          function(person_index) {
            person_id <- profiled_person_ids[person_index]
            observation_keep <- observation_person_id == person_id
            scores <- as.integer(idx$score_k[observation_keep])
            direction_sign <- if (
              length(scores) > 0L && all(scores == 0L)
            ) {
              -1
            } else if (
              length(scores) > 0L &&
                all(scores == as.integer(adjacent$transitions))
            ) {
              1
            } else {
              NA_real_
            }
            if (!is.finite(direction_sign)) return(FALSE)
            contrast_keep <- rep(
              observation_keep, each = as.integer(adjacent$transitions)
            )
            ray_margin <- direction_sign * as.numeric(
              contrast[, profiled_column[person_index], drop = FALSE]
            )
            own_margin <- ray_margin[contrast_keep]
            other_margin <- ray_margin[!contrast_keep]
            length(own_margin) > 0L &&
              all(is.finite(own_margin)) &&
              min(own_margin) > as.numeric(cone_objective_tolerance) &&
              (
                length(other_margin) == 0L ||
                  max(abs(other_margin)) <=
                    as.numeric(cone_certificate_tolerance)
              )
          },
          logical(1)
        )
        relevance$boundary_rays_valid <- all(boundary_ray_valid)
        if (!isTRUE(relevance$boundary_rays_valid)) {
          relevance$state <- "not_evaluated_boundary_ray_fallback"
          relevance$detail <- paste(
            "At least one profiled Person coordinate was not a verified",
            "strict one-sided ray confined to that Person's observations;",
            "target enumeration was retained."
          )
        } else {
          retained_observation <- !profiled_observation
          retained_contrast_row <- rep(
            retained_observation, each = as.integer(adjacent$transitions)
          )
          retained_coordinate <- rep(TRUE, ncol(contrast))
          retained_coordinate[profiled_column] <- FALSE
          quotient_contrast <- contrast[
            retained_contrast_row, retained_coordinate, drop = FALSE
          ]
          quotient_target_expansion <- target_expansion[
            , retained_coordinate, drop = FALSE
          ]
          relevance$removed_observations <- as.integer(
            sum(profiled_observation)
          )
          relevance$removed_contrast_rows <- as.integer(
            sum(!retained_contrast_row)
          )
          relevance$removed_coordinates <- as.integer(
            sum(!retained_coordinate)
          )
          quotient_cone <- if (ncol(quotient_contrast) == 0L) {
            list(
              evaluated = TRUE, certified = FALSE, solver_status = 0L,
              target_capacity = 0, target_change = 0,
              minimum_margin = 0, positive_margin = 0, strict_rows = 0L,
              direction = numeric(0), lp_calls = 0L,
              reason = "no_retained_coordinate"
            )
          } else {
            mfrmr_jml_recession_target_lp(
              lp_base = mfrmr_jml_recession_lp_base(
                quotient_contrast, representation = lp_representation
              ),
              target = as.numeric(Matrix::colSums(quotient_contrast)),
              objective_tolerance = as.numeric(cone_objective_tolerance),
              certificate_tolerance =
                as.numeric(cone_certificate_tolerance),
              timeout = lp_timeout
            )
          }
          relevance$quotient_cone_evaluated <-
            isTRUE(quotient_cone$evaluated)
          relevance$quotient_cone_certified <-
            isTRUE(quotient_cone$certified)
          relevance$quotient_cone_lp_calls <-
            as.integer(quotient_cone$lp_calls %||% 0L)
          relevance$quotient_cone_capacity <-
            as.numeric(quotient_cone$target_capacity %||% NA_real_)
          prescreen$relevance_lp_calls <- as.integer(
            prescreen$relevance_lp_calls +
              relevance$quotient_cone_lp_calls
          )
          if (!isTRUE(quotient_cone$evaluated)) {
            relevance$state <- "not_evaluated_solver_fallback"
            relevance$detail <- paste(
              "The quotient cone solver did not complete; target",
              "enumeration was retained."
            )
          } else if (isTRUE(quotient_cone$certified)) {
            relevance$state <-
              "positive_quotient_cone_target_enumeration"
            relevance$evaluated <- TRUE
            relevance$detail <- paste(
              "A strictly improving quotient cone remained after known",
              "Person boundaries were profiled out; target enumeration",
              "was retained."
            )
          } else {
            target_row_limit <- if (
              length(max_target_directions) == 1L &&
                is.finite(max_target_directions)
            ) {
              max(0L, as.integer(floor(max_target_directions / 2)))
            } else {
              nrow(quotient_target_expansion)
            }
            nullspace <- mfrmr_jml_recession_target_nullspace_screen(
              contrast_design = quotient_contrast,
              target_expansion = quotient_target_expansion,
              tolerances = nullspace_tolerances,
              max_nonzeros = max_lp_nonzeros,
              max_rows = max_lp_constraints,
              max_coordinates = max_additive_coordinates,
              max_target_rows = target_row_limit
            )
            relevance$nullspace_screen <- nullspace
            relevance$evaluated <- isTRUE(nullspace$evaluated)
            if (isTRUE(nullspace$safe_to_skip)) {
              relevance$state <-
                "negative_quotient_no_target_null_movement"
              relevance$selected_target_exclusion_certified <- TRUE
              relevance$detail <- paste(
                "The quotient strict cone was negative and every selected",
                "target annihilated its nullspace."
              )
              prescreen$state <-
                "positive_known_person_cone_target_exclusion"
              prescreen$target_enumeration_skipped <- TRUE
              targets$PositiveRecession <- FALSE
              targets$NegativeRecession <- FALSE
              targets$CandidateStatus <- ifelse(
                free_target, "finite_in_audited_subspace", "fixed"
              )
              targets$EvaluationState <- "evaluated"
              targets$ReasonCodes <- ifelse(
                free_target,
                "no_additive_recession_direction_certified",
                "no_free_coordinate"
              )
              return(finish(
                "certified_recession", !identical(model, "GPCM"),
                paste0(
                  "Certified a global joint additive recession cone fully ",
                  "explained by ", length(profiled_person_ids),
                  " ordinary free extreme Person boundary direction(s); ",
                  "the quotient cone and selected-target nullspace screen ",
                  "excluded recession for ", sum(free_target),
                  " selected free expanded targets."
                ),
                target_status = targets,
                dimensions = dimensions,
                cone_certificate = cone_certificate,
                cone_loadings = cone_loadings
              ))
            }
            relevance$state <- paste0(
              "negative_quotient_", nullspace$state, "_fallback"
            )
            relevance$detail <- paste(
              "The quotient strict cone was negative, but selected-target",
              "nullspace exclusion was not stable; target enumeration was",
              "retained."
            )
          }
        }
      }
    }
    if (target_directions > as.integer(max_target_directions)) {
      prescreen$state <- "positive_cone_target_limit"
      targets$PositiveRecession <- ifelse(free_target, NA, FALSE)
      targets$NegativeRecession <- ifelse(free_target, NA, FALSE)
      targets$CandidateStatus <- ifelse(
        free_target, "not_evaluated_size_limit", "fixed"
      )
      targets$EvaluationState <- ifelse(
        free_target, "not_evaluated", "evaluated"
      )
      targets$ReasonCodes <- ifelse(
        free_target, "bounded_target_direction_limit", "no_free_coordinate"
      )
      return(finish(
        "certified_recession", FALSE,
        paste0(
          "Certified a global joint additive recession cone, but ",
          target_directions, " selected target directions exceeded the ",
          "bounded enumeration limit."
        ),
        target_status = targets,
        dimensions = dimensions,
        cone_certificate = cone_certificate,
        cone_loadings = cone_loadings
      ))
    }
  }
  certificate_rows <- list()
  loading_rows <- list()
  positive <- negative <- rep(FALSE, nrow(targets))
  evaluated_positive <- evaluated_negative <- !free_target
  certificate_cursor <- 0L
  loading_cursor <- 0L

  for (target_index in which(free_target)) {
    target_vector <- as.numeric(target_expansion[target_index, ])
    for (requested in c("positive", "negative")) {
      signed_vector <- if (identical(requested, "positive")) {
        target_vector
      } else {
        -target_vector
      }
      solved <- mfrmr_jml_recession_target_lp(
        lp_base = lp_base,
        target = signed_vector,
        timeout = lp_timeout
      )
      prescreen$target_directions_evaluated <- as.integer(
        prescreen$target_directions_evaluated + 1L
      )
      prescreen$target_lp_calls <- as.integer(
        prescreen$target_lp_calls + (solved$lp_calls %||% 0L)
      )
      if (identical(requested, "positive")) {
        evaluated_positive[target_index] <- isTRUE(solved$evaluated)
        positive[target_index] <- isTRUE(solved$certified)
      } else {
        evaluated_negative[target_index] <- isTRUE(solved$evaluated)
        negative[target_index] <- isTRUE(solved$certified)
      }
      certificate_cursor <- certificate_cursor + 1L
      certificate_rows[[certificate_cursor]] <- data.frame(
        ParameterId = targets$ParameterId[target_index],
        RequestedDirection = requested,
        SolverStatus = as.integer(solved$solver_status),
        TargetCapacity = as.numeric(solved$target_capacity),
        TargetChange = as.numeric(solved$target_change),
        MinimumContrastMargin = as.numeric(solved$minimum_margin),
        PositiveContrastMargin = as.numeric(solved$positive_margin),
        StrictContrastRows = as.integer(solved$strict_rows),
        DirectionL1 = sum(abs(as.numeric(solved$direction)), na.rm = TRUE),
        Certified = isTRUE(solved$certified),
        stringsAsFactors = FALSE
      )
      if (isTRUE(solved$certified)) {
        keep <- which(abs(solved$direction) > 1e-9)
        if (length(keep) > 0L) {
          loading_cursor <- loading_cursor + 1L
          loading_rows[[loading_cursor]] <- data.frame(
            ParameterId = targets$ParameterId[target_index],
            RequestedDirection = requested,
            OptimizerIndex = structural_optimizer_index[keep],
            Coordinate = structural_map$Coordinate[keep],
            Loading = as.numeric(solved$direction[keep]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  evaluated <- evaluated_positive & evaluated_negative
  targets$PositiveRecession <- ifelse(evaluated_positive, positive, NA)
  targets$NegativeRecession <- ifelse(evaluated_negative, negative, NA)
  targets$CandidateStatus <- ifelse(
    !free_target,
    "fixed",
    ifelse(
      !evaluated,
      "not_evaluated_solver",
      ifelse(
        positive & negative,
        "unbounded_direction_ambiguous",
        ifelse(
          positive,
          "unbounded_high",
          ifelse(negative, "unbounded_low", "finite_in_audited_subspace")
        )
      )
    )
  )
  targets$EvaluationState <- ifelse(evaluated, "evaluated", "not_evaluated")
  targets$ReasonCodes <- ifelse(
    !free_target,
    "no_free_coordinate",
    ifelse(
      !evaluated,
      "linear_program_failed",
      ifelse(
        positive | negative,
        "certified_additive_recession_direction",
        "no_additive_recession_direction_certified"
      )
    )
  )
  certificates <- if (length(certificate_rows) == 0L) {
    empty_certificates
  } else {
    do.call(rbind, certificate_rows)
  }
  loadings <- if (length(loading_rows) == 0L) {
    empty_loadings
  } else {
    do.call(rbind, loading_rows)
  }
  any_target_certified <- any(positive | negative)
  any_certified <- any_target_certified || cone_certified
  all_evaluated <- all(evaluated)
  state <- if (!all_evaluated) {
    "not_evaluated_solver"
  } else if (any_certified) {
    "certified_recession"
  } else {
    "none_certified"
  }
  detail <- if (any_certified && joint_scope) {
    paste0(
      "Certified a joint additive recession cone; ",
      sum(positive | negative), " of ", sum(free_target),
      " selected free expanded targets have a certified direction."
    )
  } else if (any_certified) {
    paste0(
      "Certified additive structural recession for ",
      sum(positive | negative), " of ", sum(free_target),
      " free expanded targets with Person coordinates fixed."
    )
  } else if (all_evaluated) {
    paste0(
      "No additive structural recession direction was certified for ",
      sum(free_target), " free expanded targets with Person coordinates fixed."
    )
  } else {
    "At least one target-direction linear program was not evaluated successfully."
  }
  finish(
    state = state,
    complete = all_evaluated && !identical(model, "GPCM"),
    detail = detail,
    target_status = targets,
    certificates = certificates,
    loadings = loadings,
    dimensions = dimensions,
    cone_certificate = cone_certificate,
    cone_loadings = cone_loadings
  )
}

audit_mfrm_jml_structural_recession <- function(prep, idx, config, sizes,
                                                 params,
                                                 shared_geometry = NULL,
                                                 screen_global_cone = TRUE,
                                                 ...) {
  audit_mfrm_jml_additive_recession(
    prep = prep,
    idx = idx,
    config = config,
    sizes = sizes,
    params = params,
    coordinate_scope = "structural_fixed_person",
    shared_geometry = shared_geometry,
    screen_global_cone = screen_global_cone,
    ...
  )
}

audit_mfrm_jml_joint_recession <- function(prep, idx, config, sizes, params,
                                            person_boundary_audit,
                                            shared_geometry = NULL, ...) {
  person_status <- as.data.frame(
    person_boundary_audit$parameter_status %||% data.frame(),
    stringsAsFactors = FALSE
  )
  target_person_ids <- character(0)
  profiled_person_ids <- character(0)
  if (nrow(person_status) > 0L &&
      all(c("ParameterId", "ParameterStatus", "ResponseExtreme") %in%
          names(person_status))) {
    unresolved <- person_status$ParameterStatus == "weak_information" &
      person_status$ResponseExtreme %in% c("low", "high")
    target_person_ids <- as.character(person_status$ParameterId[unresolved])
    profiled <- person_status$ParameterStatus %in%
      c("unbounded_low", "unbounded_high") &
      person_status$ResponseExtreme %in% c("low", "high")
    profiled_person_ids <- as.character(person_status$ParameterId[profiled])
  }
  audit_mfrm_jml_additive_recession(
    prep = prep,
    idx = idx,
    config = config,
    sizes = sizes,
    params = params,
    coordinate_scope = "joint_person_structural",
    target_person_ids = target_person_ids,
    profiled_person_ids = profiled_person_ids,
    shared_geometry = shared_geometry,
    screen_global_cone = TRUE,
    ...
  )
}
