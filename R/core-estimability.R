mfrmr_empty_sparse_matrix <- function(nrow = 0L, ncol = 0L) {
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(as.integer(nrow), as.integer(ncol)))
}

mfrmr_estimability_map <- function(Block = character(0), Coordinate = character(0), Facet = character(0),
    Level = character(0), ReferenceLevel = character(0), Constraint = character(0), OptimizerIndex = integer(0)) {
    data.frame(Block = as.character(Block), Coordinate = as.character(Coordinate), Facet = as.character(Facet),
        Level = as.character(Level), ReferenceLevel = as.character(ReferenceLevel), Constraint = as.character(Constraint),
        OptimizerIndex = as.integer(OptimizerIndex), stringsAsFactors = FALSE)
}

mfrmr_constraint_jacobian_sparse <- function(spec, block) {
    levels <- as.character(spec$levels %||% character(0))
    n_levels <- length(levels)
    anchors <- as.numeric(spec$anchors %||% rep(NA_real_, n_levels))
    groups <- as.character(spec$groups %||% rep(NA_character_, n_levels))
    if (length(anchors) != n_levels || length(groups) != n_levels) {
        stop("Internal estimability audit found a malformed facet constraint.", call. = FALSE)
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
        suffix <- if (nzchar(group))
            paste0("|group=", group)
        else ""
        map_rows[[col_cursor]] <<- mfrmr_estimability_map(Block = block, Coordinate = paste0(block, ":",
            levels[level_idx], suffix), Facet = block, Level = levels[level_idx], ReferenceLevel = levels[reference_idx],
            Constraint = constraint, OptimizerIndex = NA_integer_)
    }
    group_ids <- unique(stats::na.omit(groups[free_idx]))
    group_ids <- group_ids[nzchar(group_ids)]
    grouped_idx <- integer(0)
    if (length(group_ids) > 0L) {
        for (group_id in group_ids) {
            group_levels <- which(groups == group_id)
            free_in_group <- group_levels[is.na(anchors[group_levels])]
            grouped_idx <- c(grouped_idx, free_in_group)
            if (length(free_in_group) <= 1L)
                next
            reference <- free_in_group[length(free_in_group)]
            for (level_idx in free_in_group[-length(free_in_group)]) {
                add_column(level_idx = level_idx, reference_idx = reference, constraint = "group_sum",
                  group = group_id)
            }
        }
    }
    ungrouped_idx <- setdiff(free_idx, grouped_idx)
    if (length(ungrouped_idx) > 0L) {
        if (isTRUE(spec$centered)) {
            if (length(ungrouped_idx) > 1L) {
                reference <- ungrouped_idx[length(ungrouped_idx)]
                for (level_idx in ungrouped_idx[-length(ungrouped_idx)]) {
                  add_column(level_idx = level_idx, reference_idx = reference, constraint = if (identical(spec$centering_scope,
                    "facet")) {
                    "whole_facet_sum_zero"
                  }
                  else {
                    "sum_zero"
                  })
                }
            }
        }
        else {
            for (level_idx in ungrouped_idx) {
                col_cursor <- col_cursor + 1L
                row_i <- c(row_i, level_idx)
                col_j <- c(col_j, col_cursor)
                values <- c(values, 1)
                map_rows[[col_cursor]] <- mfrmr_estimability_map(Block = block, Coordinate = paste0(block,
                  ":", levels[level_idx]), Facet = block, Level = levels[level_idx], ReferenceLevel = "",
                  Constraint = "free_location", OptimizerIndex = NA_integer_)
            }
        }
    }
    expected <- as.integer(spec$n_params %||% 0L)
    if (col_cursor != expected) {
        stop("Internal estimability audit/free-parameter mismatch for ", shQuote(block), ": audit built ",
            col_cursor, " coordinate(s), optimizer expects ", expected, ".", call. = FALSE)
    }
    jacobian <- if (col_cursor == 0L) {
        mfrmr_empty_sparse_matrix(n_levels, 0L)
    }
    else {
        Matrix::sparseMatrix(i = row_i, j = col_j, x = values, dims = c(n_levels, col_cursor))
    }
    map <- if (length(map_rows) == 0L) {
        mfrmr_estimability_map()
    }
    else {
        do.call(rbind, map_rows)
    }
    rownames(map) <- NULL
    dimnames(jacobian) <- list(levels, if (nrow(map) > 0L) map$Coordinate else character(0))
    list(jacobian = jacobian, map = map)
}

mfrmr_interaction_jacobian_sparse <- function(spec) {
    n_a <- as.integer(spec$n_a %||% 0L)
    n_b <- as.integer(spec$n_b %||% 0L)
    n_params <- as.integer(spec$n_params %||% 0L)
    n_cells <- n_a * n_b
    if (n_params == 0L) {
        return(list(jacobian = mfrmr_empty_sparse_matrix(n_cells, 0L), map = mfrmr_estimability_map()))
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
            row_i <- c(row_i, cell(a, b), cell(a, n_b), cell(n_a, b), cell(n_a, n_b))
            col_j <- c(col_j, rep(col_cursor, 4L))
            values <- c(values, 1, -1, -1, 1)
            map_rows[[col_cursor]] <- mfrmr_estimability_map(Block = "interactions", Coordinate = paste0(spec$name,
                ":", spec$levels_a[a], "*", spec$levels_b[b]), Facet = spec$name, Level = paste(spec$levels_a[a],
                spec$levels_b[b], sep = "*"), ReferenceLevel = paste(spec$levels_a[n_a], spec$levels_b[n_b],
                sep = "*"), Constraint = "two_way_sum_zero_margins", OptimizerIndex = NA_integer_)
        }
    }
    if (col_cursor != n_params) {
        stop("Internal interaction estimability dimension mismatch.", call. = FALSE)
    }
    list(jacobian = Matrix::sparseMatrix(i = row_i, j = col_j, x = values, dims = c(n_cells, n_params)),
        map = do.call(rbind, map_rows))
}

mfrmr_step_jacobian_sparse <- function(config, sizes) {
    n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
    specs <- mfrmr_step_specs(config)
    n_scope <- length(specs)
    n_rows <- n_scope * n_steps
    step_map <- config$step_constraint_map
    if (!is.null(step_map)) {
        rows <- step_map$free
        map <- if (!length(rows))
            mfrmr_estimability_map()
        else {
            scope <- names(specs)[(rows - 1L)%/%n_steps + 1L]
            transition <- (rows - 1L)%%n_steps + 1L
            mfrmr_estimability_map(Block = "steps", Coordinate = paste0("steps:", scope, ":transition",
                transition), Facet = if (config$model == "RSM")
                "shared_rating_scale"
            else config$step_facet, Level = scope, ReferenceLevel = NA_character_, Constraint = "fixed_or_shared_centered_steps",
                OptimizerIndex = NA_integer_)
        }
        return(list(jacobian = step_map$basis, map = map))
    }
    n_cols <- sum(vapply(specs, function(spec) {
        as.integer(spec$n_params)
    }, integer(1)))
    if (n_cols == 0L) {
        return(list(jacobian = mfrmr_empty_sparse_matrix(n_rows, 0L), map = mfrmr_estimability_map()))
    }
    row_i <- integer(0)
    col_j <- integer(0)
    values <- numeric(0)
    map_rows <- vector("list", n_cols)
    col_cursor <- 0L
    scope_levels <- if (identical(config$model, "RSM")) {
        "shared"
    }
    else {
        as.character(config$facet_levels[[config$step_facet]])
    }
    for (scope in seq_len(n_scope)) {
        row_offset <- (scope - 1L) * n_steps
        free_index <- as.integer(specs[[scope]]$free_index)
        if (length(free_index) <= 1L)
            next
        reference_transition <- free_index[length(free_index)]
        estimated_transitions <- free_index[-length(free_index)]
        for (transition in estimated_transitions) {
            col_cursor <- col_cursor + 1L
            row_i <- c(row_i, row_offset + transition, row_offset + reference_transition)
            col_j <- c(col_j, col_cursor, col_cursor)
            values <- c(values, 1, -1)
            map_rows[[col_cursor]] <- mfrmr_estimability_map(Block = "steps", Coordinate = paste0("steps:",
                scope_levels[scope], ":transition", transition), Facet = if (identical(config$model,
                "RSM")) {
                "shared_rating_scale"
            }
            else {
                config$step_facet
            }, Level = scope_levels[scope], ReferenceLevel = paste0("transition", reference_transition),
                Constraint = if (any(is.finite(specs[[scope]]$anchors))) {
                  "within_ladder_sum_zero_with_fixed_anchors"
                }
                else {
                  "within_ladder_sum_zero"
                }, OptimizerIndex = NA_integer_)
        }
    }
    list(jacobian = Matrix::sparseMatrix(i = row_i, j = col_j, x = values, dims = c(n_rows, n_cols)),
        map = do.call(rbind, map_rows))
}

mfrmr_estimability_set_optimizer_index <- function(map, block, slices) {
    if (nrow(map) == 0L)
        return(map)
    idx <- slices[[block]] %||% integer(0)
    if (length(idx) == nrow(map)) {
        map$OptimizerIndex <- as.integer(idx)
    }
    map
}

mfrmr_estimability_eta_design <- function(prep, idx, config, sizes, include_person, include_population_beta) {
    n_obs <- length(idx$score_k)
    slices <- build_param_slices(sizes)
    blocks <- list()
    maps <- list()
    append_block <- function(name, matrix, map) {
        if (ncol(matrix) == 0L)
            return(invisible(NULL))
        blocks[[length(blocks) + 1L]] <<- matrix
        maps[[length(maps) + 1L]] <<- map
        invisible(NULL)
    }
    if (isTRUE(include_person)) {
        theta <- mfrmr_constraint_jacobian_sparse(config$theta_spec, "Person")
        theta$map <- mfrmr_estimability_set_optimizer_index(theta$map, "theta", slices)
        append_block("theta", theta$jacobian[idx$person, , drop = FALSE], theta$map)
    }
    for (facet in config$facet_names) {
        facet_block <- mfrmr_constraint_jacobian_sparse(config$facet_specs[[facet]], facet)
        sign <- as.numeric(config$facet_signs[[facet]] %||% -1)
        facet_block$map <- mfrmr_estimability_set_optimizer_index(facet_block$map, facet, slices)
        append_block(facet, sign * facet_block$jacobian[idx$facets[[facet]], , drop = FALSE], facet_block$map)
    }
    interaction_specs <- config$interaction_specs %||% list()
    if (length(interaction_specs) > 0L) {
        interaction_matrices <- list()
        interaction_maps <- list()
        for (name in names(interaction_specs)) {
            interaction <- mfrmr_interaction_jacobian_sparse(interaction_specs[[name]])
            interaction_matrices[[length(interaction_matrices) + 1L]] <- interaction$jacobian[idx$interactions[[name]],
                , drop = FALSE]
            interaction_maps[[length(interaction_maps) + 1L]] <- interaction$map
        }
        interaction_matrix <- do.call(cbind, interaction_matrices)
        interaction_map <- do.call(rbind, interaction_maps)
        rownames(interaction_map) <- NULL
        interaction_map <- mfrmr_estimability_set_optimizer_index(interaction_map, "interactions", slices)
        append_block("interactions", interaction_matrix, interaction_map)
    }
    if (isTRUE(include_population_beta) && isTRUE(config$population_spec$active) && as.integer(sizes[["beta"]] %||%
        0L) > 0L) {
        pop <- config$population_spec
        lookup <- as.integer(pop$person_lookup %||% integer(0))
        if (length(lookup) != config$n_person || anyNA(lookup[idx$person])) {
            stop("Internal estimability audit could not align population-model persons.", call. = FALSE)
        }
        beta_matrix <- as.matrix(pop$design_matrix)[lookup[idx$person], , drop = FALSE]
        beta_names <- colnames(beta_matrix)
        if (is.null(beta_names))
            beta_names <- paste0("beta_", seq_len(ncol(beta_matrix)))
        beta_map <- mfrmr_estimability_map(Block = rep("beta", ncol(beta_matrix)), Coordinate = paste0("beta:",
            beta_names), Facet = rep("PersonPopulation", ncol(beta_matrix)), Level = beta_names, ReferenceLevel = rep("",
            ncol(beta_matrix)), Constraint = rep("population_location", ncol(beta_matrix)), OptimizerIndex = as.integer(slices[["beta"]] %||%
            rep(NA_integer_, ncol(beta_matrix))))
        append_block("beta", Matrix::Matrix(beta_matrix, sparse = TRUE), beta_map)
    }
    design <- if (length(blocks) == 0L) {
        mfrmr_empty_sparse_matrix(n_obs, 0L)
    }
    else {
        do.call(cbind, blocks)
    }
    map <- if (length(maps) == 0L) {
        mfrmr_estimability_map()
    }
    else {
        do.call(rbind, maps)
    }
    rownames(map) <- NULL
    list(design = design, map = map)
}

mfrmr_estimability_adjacent_design <- function(prep, idx, config, sizes, include_person = NULL, include_population_beta = NULL) {
    if (is.null(include_person)) {
        include_person <- identical(config$method, "JML")
    }
    if (is.null(include_population_beta)) {
        include_population_beta <- identical(config$method, "MML")
    }
    eta <- mfrmr_estimability_eta_design(prep = prep, idx = idx, config = config, sizes = sizes, include_person = include_person,
        include_population_beta = include_population_beta)
    n_obs <- nrow(eta$design)
    n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
    if (n_steps <= 0L) {
        stop("Estimability audit requires at least two response categories.", call. = FALSE)
    }
    transition <- rep(seq_len(n_steps), each = n_obs)
    obs_index <- rep(seq_len(n_obs), times = n_steps)
    eta_repeated <- eta$design[obs_index, , drop = FALSE]
    step <- mfrmr_step_jacobian_sparse(config, sizes)
    step$map <- mfrmr_estimability_set_optimizer_index(step$map, "steps", build_param_slices(sizes))
    expanded_step_index <- if (identical(config$model, "RSM")) {
        transition
    }
    else {
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
    list(design = methods::as(design, "dgCMatrix"), map = map, transition_rows = nrow(design), observation_rows = n_obs,
        transitions = n_steps)
}
