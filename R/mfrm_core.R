
weighted_mean <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok))
        return(NA_real_)
    sum(x[ok] * w[ok])/sum(w[ok])
}

get_weights <- function(df) {
    if (!is.null(df) && "Weight" %in% names(df)) {
        w <- suppressWarnings(as.numeric(df$Weight))
        w <- ifelse(is.finite(w) & w > 0, w, 0)
        return(w)
    }
    rep(1, nrow(df))
}

compact_population_spec <- function(population = NULL, person_levels = character(0)) {
    pop <- population %||% list()
    person_levels <- as.character(person_levels %||% character(0))
    person_table <- pop$person_table
    person_id <- as.character(pop$person_id %||% "Person")
    person_lookup <- integer(0)
    if (!is.null(person_table) && is.data.frame(person_table) && length(person_levels) > 0 && person_id %in%
        names(person_table)) {
        person_lookup <- match(person_levels, as.character(person_table[[person_id]]))
    }
    out <- list(active = isTRUE(pop$active), source = as.character(pop$source %||% "none"), identification_role = as.character(pop$identification_role %||%
        "not_applicable"), posterior_basis = as.character(pop$posterior_basis %||% if (isTRUE(pop$active)) "population_model" else "legacy_mml"),
        formula = pop$formula %||% NULL, terms = pop$terms %||% NULL, person_id = person_id, design_matrix = pop$design_matrix %||%
            NULL, design_columns = pop$design_columns %||% if (!is.null(pop$design_matrix)) colnames(pop$design_matrix) else NULL,
        xlevels = pop$xlevels %||% NULL, contrasts = pop$contrasts %||% NULL, coefficients = pop$coefficients %||%
            NULL, sigma2 = pop$sigma2 %||% NULL, converged = isTRUE(pop$converged), estimation_converged = isTRUE(pop$estimation_converged),
        inference_ready = isTRUE(pop$inference_ready), converged_basis = as.character(pop$converged_basis %||%
            "legacy_unspecified"), policy = pop$policy %||% NULL, person_lookup = person_lookup, person_rows = if (!is.null(person_table) &&
            is.data.frame(person_table)) {
            nrow(person_table)
        } else {
            0L
        }, included_persons = as.character(pop$included_persons %||% character(0)), omitted_persons = as.character(pop$omitted_persons %||%
            character(0)), response_rows_retained = as.integer(pop$response_rows_retained %||% NA_integer_),
        response_rows_omitted = as.integer(pop$response_rows_omitted %||% NA_integer_), notes = pop$notes %||%
            character(0))
    if (!is.null(pop$fixed_sd))
        out$fixed_sd <- pop$fixed_sd
    out
}

build_step_constraint <- function(n_steps, anchors = NULL, target = 0, scope = "shared") {
    n_steps <- as.integer(n_steps %||% 0L)
    if (length(n_steps) != 1L || is.na(n_steps) || n_steps < 0L) {
        stop("`n_steps` must be one non-negative integer.", call. = FALSE)
    }
    target <- as.numeric(target %||% NA_real_)
    if (length(target) != 1L || !is.finite(target)) {
        stop("Step-constraint target must be one finite number.", call. = FALSE)
    }
    anchor_values <- rep(NA_real_, n_steps)
    names(anchor_values) <- as.character(seq_len(n_steps))
    if (!is.null(anchors) && length(anchors) > 0L) {
        if (is.null(names(anchors)) || anyNA(names(anchors)) || any(!nzchar(names(anchors))) || anyDuplicated(names(anchors))) {
            stop("Step anchors must have unique transition-index names.", call. = FALSE)
        }
        anchor_index <- suppressWarnings(as.integer(names(anchors)))
        if (anyNA(anchor_index) || any(anchor_index < 1L | anchor_index > n_steps) || any(as.character(anchor_index) !=
            names(anchors))) {
            stop("Step-anchor names must be canonical transition indices in range.", call. = FALSE)
        }
        anchor_numeric <- as.numeric(anchors)
        if (any(!is.finite(anchor_numeric))) {
            stop("Step-anchor values must be finite.", call. = FALSE)
        }
        anchor_values[anchor_index] <- anchor_numeric
    }
    free_index <- which(is.na(anchor_values))
    if (length(free_index) == 0L && !isTRUE(all.equal(sum(anchor_values), target, tolerance = 1e-12))) {
        stop("Fully anchored step coordinates are incompatible with the declared sum constraint.", call. = FALSE)
    }
    reference_index <- if (length(free_index) > 0L) {
        free_index[length(free_index)]
    }
    else {
        NA_integer_
    }
    list(scope = as.character(scope %||% "shared")[1], n_steps = n_steps, anchors = anchor_values, target = target,
        free_index = as.integer(free_index), reference_index = as.integer(reference_index), n_params = as.integer(max(length(free_index) -
            1L, 0L)), identification = "within_ladder_sum_to_target_with_fixed_anchors")
}

mfrmr_default_step_specs <- function(config) {
    n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
    scopes <- if (identical(config$model, "RSM")) {
        "shared"
    }
    else {
        step_facet <- as.character(config$step_facet %||% NA_character_)
        if (length(step_facet) != 1L || is.na(step_facet) || !step_facet %in% as.character(config$facet_names %||%
            character(0))) {
            stop("PCM model requires `step_facet` to name a declared facet.", call. = FALSE)
        }
        as.character(config$facet_levels[[step_facet]] %||% character(0))
    }
    specs <- lapply(scopes, function(scope) {
        build_step_constraint(n_steps = n_steps, scope = scope)
    })
    names(specs) <- scopes
    specs
}

mfrmr_step_specs <- function(config) {
    specs <- config$step_specs %||% mfrmr_default_step_specs(config)
    if (!is.list(specs)) {
        stop("Internal step-constraint specification must be a list.", call. = FALSE)
    }
    expected_scopes <- if (identical(config$model, "RSM")) {
        "shared"
    }
    else {
        as.character(config$facet_levels[[config$step_facet]] %||% character(0))
    }
    if (!identical(names(specs), expected_scopes)) {
        stop("Step-constraint scopes do not match the declared model ownership.", call. = FALSE)
    }
    n_steps <- max(as.integer(config$n_cat %||% 0L) - 1L, 0L)
    valid <- vapply(seq_along(specs), function(i) {
        spec <- specs[[i]]
        is.list(spec) && identical(as.integer(spec$n_steps), n_steps) && identical(as.character(spec$scope),
            expected_scopes[i]) && length(spec$anchors) == n_steps && identical(as.integer(spec$free_index),
            unname(which(is.na(spec$anchors)))) && identical(as.integer(spec$n_params), as.integer(max(length(spec$free_index) -
            1L, 0L)))
    }, logical(1))
    if (!all(valid)) {
        stop("Step-constraint specification is inconsistent with the model.", call. = FALSE)
    }
    specs
}

build_gpcm_slope_spec <- function(levels, slope_facet, step_facet) {
    slope_levels <- as.character(levels %||% character(0))
    list(active = TRUE, slope_facet = as.character(slope_facet %||% NA_character_), step_facet = as.character(step_facet %||%
        NA_character_), levels = slope_levels, n_params = max(length(slope_levels) - 1L, 0L), identification = "sum_to_zero_log_slopes",
        scale_reference = "geometric_mean_one", reduction_reference = "PCM when all slopes equal 1")
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
                stop("`facet_interactions` must use explicit two-way terms such as ", "`\"Rater:Criterion\"`. Problem term: ",
                  shQuote(term), ".", call. = FALSE)
            }
            parts
        })
    }
    else if (is.list(facet_interactions) && !is.data.frame(facet_interactions)) {
        lapply(facet_interactions, function(term) {
            parts <- trimws(as.character(term))
            parts <- parts[!is.na(parts) & nzchar(parts)]
            if (length(parts) != 2L) {
                stop("List entries in `facet_interactions` must each contain exactly ", "two facet names, for example `list(c(\"Rater\", \"Criterion\"))`.",
                  call. = FALSE)
            }
            parts
        })
    }
    else {
        stop("`facet_interactions` must be NULL, a character vector of `A:B` terms, ", "or a list of length-two character vectors.",
            call. = FALSE)
    }
    canonical <- character(length(terms))
    names_out <- character(length(terms))
    for (i in seq_along(terms)) {
        parts <- terms[[i]]
        if (any(parts == "Person")) {
            stop("`facet_interactions` currently supports only non-person facet pairs. ", "Remove `Person` from ",
                shQuote(paste(parts, collapse = ":")), ".", call. = FALSE)
        }
        if (!all(parts %in% facet_names)) {
            missing <- setdiff(parts, facet_names)
            stop("`facet_interactions` references facet(s) not supplied in `facets`: ", paste(shQuote(missing),
                collapse = ", "), ". Available non-person ", "facets are: ", paste(shQuote(facet_names),
                collapse = ", "), ".", call. = FALSE)
        }
        if (length(unique(parts)) != 2L) {
            stop("`facet_interactions` terms must name two distinct facets. Problem term: ", shQuote(paste(parts,
                collapse = ":")), ".", call. = FALSE)
        }
        canonical[i] <- paste(sort(parts), collapse = "\r")
        names_out[i] <- paste(parts, collapse = ":")
    }
    if (any(duplicated(canonical))) {
        dup <- names_out[duplicated(canonical)]
        stop("`facet_interactions` contains duplicate facet pairs: ", paste(shQuote(dup), collapse = ", "),
            ". Each two-way pair may be ", "specified only once.", call. = FALSE)
    }
    names(terms) <- names_out
    terms
}

build_facet_interaction_specs <- function(prep, facet_interactions = NULL, min_obs_per_interaction = 10,
    interaction_policy = c("warn", "error", "silent")) {
    interaction_policy <- match.arg(interaction_policy)
    if (!is.numeric(min_obs_per_interaction) || length(min_obs_per_interaction) != 1L || !is.finite(min_obs_per_interaction) ||
        min_obs_per_interaction < 0) {
        stop("`min_obs_per_interaction` must be a single non-negative finite number.", call. = FALSE)
    }
    min_obs_per_interaction <- as.numeric(min_obs_per_interaction)
    terms <- normalize_facet_interactions(facet_interactions, prep$facet_names)
    if (length(terms) == 0L)
        return(list())
    specs <- lapply(names(terms), function(name) {
        facets <- terms[[name]]
        facet_a <- facets[1]
        facet_b <- facets[2]
        levels_a <- prep$levels[[facet_a]]
        levels_b <- prep$levels[[facet_b]]
        n_a <- length(levels_a)
        n_b <- length(levels_b)
        if (n_a < 2L || n_b < 2L) {
            stop("Facet interaction ", shQuote(name), " is not estimable because ", "both facets must have at least two observed levels.",
                call. = FALSE)
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
            msg <- paste0("Facet interaction ", shQuote(name), " has ", sparse_count, " of ", n_cells,
                " cells below `min_obs_per_interaction = ", min_obs_per_interaction, "` (", zero_count,
                " zero-count cells). ", "Interaction estimates in sparse cells are weakly identified; ",
                "use the term only when the pair is substantively justified.")
            if (identical(interaction_policy, "error")) {
                stop(msg, call. = FALSE)
            }
            warning(msg, call. = FALSE)
        }
        list(name = name, facets = facets, facet_a = facet_a, facet_b = facet_b, levels_a = levels_a,
            levels_b = levels_b, n_a = n_a, n_b = n_b, n_cells = n_cells, n_params = (n_a - 1L) * (n_b -
                1L), cell_n = cell_n, cell_weighted_n = cell_weighted, sparse_count = sparse_count, zero_count = zero_count,
            min_obs_per_interaction = min_obs_per_interaction, identification = "two_way_sum_to_zero_margins",
            interpretation = paste("Positive estimates indicate scores higher than expected under the",
                "additive main-effects model for that facet-level combination."))
    })
    names(specs) <- names(terms)
    specs
}

build_facet_constraint <- function(levels, anchors = NULL, groups = NULL, group_values = NULL, centered = TRUE,
    centering_scope = c("remainder", "facet")) {
    centering_scope <- match.arg(centering_scope)
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
    spec <- list(levels = lvl, anchors = anchors_vec, groups = groups_vec, group_values = group_values_map,
        centered = centered)
    if (identical(centering_scope, "facet")) {
        spec$centering_scope <- "facet"
        spec$groups[!is.na(spec$groups) & !nzchar(spec$groups)] <- NA_character_
        grouped <- !is.na(spec$groups)
        fixed_total <- sum(spec$anchors[!grouped], na.rm = TRUE)
        fixed_magnitude <- sum(abs(spec$anchors[!grouped]), na.rm = TRUE)
        for (gid in unique(spec$groups[grouped])) {
            members <- which(spec$groups == gid)
            value <- if (gid %in% names(spec$group_values))
                spec$group_values[[gid]]
            else 0
            if (length(value) != 1L || !is.finite(value)) {
                stop("A group-mean constraint requires a finite target.", call. = FALSE)
            }
            spec$group_values[gid] <- value
            target <- value * length(members)
            if (all(is.finite(spec$anchors[members]))) {
                residual <- sum(spec$anchors[members]) - target
                scale <- max(1, abs(target), sum(abs(spec$anchors[members])))
                if (abs(residual) > 1e-12 * scale) {
                  stop("Fixed anchors conflict with the mean of group '", gid, "'.", call. = FALSE)
                }
            }
            fixed_total <- fixed_total + target
            fixed_magnitude <- fixed_magnitude + abs(target)
        }
        spec$remainder_sum <- if (isTRUE(centered))
            -fixed_total
        else 0
        remaining <- which(is.na(spec$anchors) & !grouped)
        if (isTRUE(centered) && !length(remaining) && abs(fixed_total) > 1e-12 * max(1, fixed_magnitude)) {
            stop("Fixed anchors and group means conflict with whole-facet centering. ", "Revise the declared constraints or the facet's centering choice.",
                call. = FALSE)
        }
    }
    spec$n_params <- count_facet_params(spec)
    spec
}

count_facet_params <- function(spec) {
    anchors <- spec$anchors
    groups <- spec$groups
    free_idx <- which(is.na(anchors))
    if (length(free_idx) == 0)
        return(0)
    n_params <- 0
    group_ids <- unique(na.omit(groups[free_idx]))
    if (length(group_ids) > 0) {
        for (gid in group_ids) {
            group_levels <- which(groups == gid)
            free_in_group <- group_levels[is.na(anchors[group_levels])]
            k <- length(free_in_group)
            if (k > 1)
                n_params <- n_params + (k - 1)
        }
    }
    ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
    m <- length(ungrouped_idx)
    if (m == 0)
        return(n_params)
    if (spec$centered) {
        n_params <- n_params + max(m - 1, 0)
    }
    else {
        n_params <- n_params + m
    }
    n_params
}


constraint_grad_project <- function(grad_expanded, spec) {
    anchors <- spec$anchors
    groups <- spec$groups
    centered <- isTRUE(spec$centered)
    free_idx <- which(is.na(anchors))
    if (length(free_idx) == 0 || spec$n_params == 0)
        return(numeric(0))
    grad_free <- numeric(spec$n_params)
    pos <- 1L
    group_ids <- unique(na.omit(groups[free_idx]))
    if (length(group_ids) > 0) {
        for (gid in group_ids) {
            group_levels <- which(groups == gid)
            free_in_group <- group_levels[is.na(anchors[group_levels])]
            k <- length(free_in_group)
            if (k <= 1)
                next
            last_g <- grad_expanded[free_in_group[k]]
            for (j in seq_len(k - 1)) {
                grad_free[pos] <- grad_expanded[free_in_group[j]] - last_g
                pos <- pos + 1L
            }
        }
    }
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
        }
        else {
            for (j in seq_len(m)) {
                grad_free[pos] <- grad_expanded[ungrouped_idx[j]]
                pos <- pos + 1L
            }
        }
    }
    grad_free
}

build_param_sizes <- function(config) {
    mfrm_validate_facet_names(config$facet_names)
    n_steps <- max(config$n_cat - 1, 0)
    sizes <- list(theta = if (config$method == "JML") config$theta_spec$n_params else 0)
    for (facet in config$facet_names) {
        sizes[[facet]] <- config$facet_specs[[facet]]$n_params
    }
    interaction_n <- sum(vapply(config$interaction_specs %||% list(), function(spec) {
        as.integer(spec$n_params %||% 0L)
    }, integer(1)))
    if (interaction_n > 0L) {
        sizes$interactions <- interaction_n
    }
    if (config$model != "RSM") {
        if (is.null(config$step_facet) || !config$step_facet %in% config$facet_names) {
            stop("PCM model requires 'step_facet' to name one of the facet columns: ", paste(config$facet_names,
                collapse = ", "), ". ", "Supply step_facet = '<name>'.", call. = FALSE)
        }
    }
    step_specs <- mfrmr_step_specs(config)
    if (!is.null(config$step_constraint_map) && !isNamespaceLoaded("Matrix"))
        loadNamespace("Matrix")
    sizes$steps <- ncol(config$step_constraint_map$basis) %||% sum(vapply(step_specs, function(spec) {
        as.integer(spec$n_params)
    }, integer(1)))
    if (identical(config$model, "GPCM")) {
        gpcm_spec <- config$gpcm_spec %||% list()
        sizes$log_slopes <- as.integer(gpcm_spec$n_params %||% 0L)
    }
    if (identical(config$method, "MML") && isTRUE(config$population_spec$active)) {
        design_matrix <- config$population_spec$design_matrix
        if (is.null(design_matrix) || !is.matrix(design_matrix)) {
            stop("Active latent-regression MML requires a numeric person-level design matrix.", call. = FALSE)
        }
        sizes$beta <- ncol(design_matrix)
        sizes$log_sigma2 <- if (is.null(config$population_spec$fixed_sd))
            1L
        else 0L
    }
    sizes
}

prepare_constraint_specs <- function(prep, anchor_df = NULL, group_anchor_df = NULL, noncenter_facet = "Person",
    dummy_facets = character(0), centering_scope = "remainder", anchor_review = NULL) {
    centering_scope <- match.arg(centering_scope, c("facet", "remainder", "anchor_origin"))
    facet_names <- prep$facet_names
    all_facets <- c("Person", facet_names)
    anchor_review <- anchor_review %||% audit_anchor_tables(prep = prep, anchor_df = anchor_df, group_anchor_df = group_anchor_df,
        noncenter_facet = noncenter_facet, dummy_facets = dummy_facets, centering_scope = centering_scope)
    anchor_df <- anchor_review$anchors
    group_anchor_df <- anchor_review$group_anchors
    anchor_origin <- identical(centering_scope, "anchor_origin")
    anchored_facets <- unique(c(as.character(anchor_df$Facet), as.character(group_anchor_df$Facet)))
    if (anchor_origin) {
        if (noncenter_facet %in% c(anchored_facets, dummy_facets)) {
            stop("`centering_scope = 'anchor_origin'` requires an unanchored, active `noncenter_facet`. Choose a facet without direct anchors, group anchors or dummy constraints.",
                call. = FALSE)
        }
    }
    map_scope <- if (anchor_origin)
        "facet"
    else centering_scope
    is_centered <- function(facet) {
        !identical(noncenter_facet, facet) && !(anchor_origin && facet %in% anchored_facets)
    }
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
                df <- dplyr::filter(anchor_df, Facet == facet)
                if (nrow(df) > 0) {
                  anchors <- setNames(df$Anchor, df$Level)
                  anchors <- anchors[names(anchors) %in% levels]
                  if (length(anchors) > 0)
                    anchor_map[[facet]] <- anchors
                }
            }
            if (nrow(group_anchor_df) > 0) {
                df <- dplyr::filter(group_anchor_df, Facet == facet)
                if (nrow(df) > 0) {
                  groups <- setNames(df$Group, df$Level)
                  groups <- groups[names(groups) %in% levels]
                  if (length(groups) > 0)
                    group_map[[facet]] <- groups
                  group_vals <- dplyr::mutate(dplyr::summarize(dplyr::group_by(dplyr::filter(dplyr::distinct(dplyr::select(df, Group, GroupValue)),
                    !is.na(Group)), Group), GroupValue = {
                    vals <- GroupValue[!is.na(GroupValue)]
                    if (length(vals) == 0)
                      NA_real_
                    else vals[1]
                  }, .groups = "drop"), GroupValue = ifelse(is.na(GroupValue), 0, GroupValue))
                  if (nrow(group_vals) > 0) {
                    group_values[[facet]] <- setNames(group_vals$GroupValue, group_vals$Group)
                  }
                }
            }
        }
    }
    theta_spec <- build_facet_constraint(levels = prep$levels$Person, anchors = anchor_map$Person, groups = group_map$Person,
        group_values = group_values$Person, centered = is_centered("Person"), centering_scope = map_scope)
    facet_specs <- lapply(facet_names, function(facet) {
        build_facet_constraint(levels = prep$levels[[facet]], anchors = anchor_map[[facet]], groups = group_map[[facet]],
            group_values = group_values[[facet]], centered = is_centered(facet), centering_scope = map_scope)
    })
    names(facet_specs) <- facet_names
    anchor_summary <- tibble::tibble(Facet = all_facets, AnchoredLevels = vapply(anchor_map, function(x) length(x),
        integer(1)), GroupAnchors = vapply(group_map, function(x) length(unique(x)), integer(1)), DummyFacet = all_facets %in%
        dummy_facets)
    list(theta_spec = theta_spec, facet_specs = facet_specs, anchor_summary = anchor_summary, anchor_review = anchor_review)
}

build_facet_signs <- function(facet_names, positive_facets = character(0)) {
    positives <- intersect(positive_facets, facet_names)
    dropped <- setdiff(positive_facets, facet_names)
    if (length(dropped) > 0L) {
        warning("`positive_facets` contains value(s) not present in the model: ", paste(shQuote(dropped),
            collapse = ", "), " (available: ", paste(shQuote(facet_names), collapse = ", "), "). Unknown entries are ignored; check for typos (case-sensitive).",
            call. = FALSE)
    }
    signs <- setNames(ifelse(facet_names %in% positives, 1, -1), facet_names)
    list(signs = signs, positive_facets = positives)
}


build_estimation_config <- function(prep, model, method, step_facet, slope_facet, interaction_specs = NULL,
    weight_col, facet_signs, positive_facets, noncenter_facet, dummy_facets, anchor_df, group_anchor_df,
    population = NULL, centering_scope = "remainder", anchor_review = NULL, constraint_specs = NULL,
    step_constraints = NULL) {
    if (identical(centering_scope, "anchor_origin") && (!identical(method, "JML") || !model %in% c("RSM",
        "PCM") || length(interaction_specs) > 0L)) {
        stop("`centering_scope = 'anchor_origin'` requires ordinary additive RSM/PCM JML.", call. = FALSE)
    }
    config <- list(model = model, method = method, n_person = length(prep$levels$Person), person_levels = prep$levels$Person,
        n_cat = prep$rating_max - prep$rating_min + 1, facet_names = prep$facet_names, facet_levels = prep$levels[prep$facet_names],
        step_facet = step_facet, slope_facet = slope_facet, interaction_specs = interaction_specs %||%
            list(), facet_interactions = names(interaction_specs %||% list()), keep_original = isTRUE(prep$keep_original),
        rating_min = prep$rating_min, rating_max = prep$rating_max, rating_range_source = prep$rating_range_source %||%
            NA_character_, rating_min_source = prep$rating_min_source %||% NA_character_, rating_max_source = prep$rating_max_source %||%
            NA_character_, score_map = prep$score_map)
    config$weight_col <- if (!is.null(weight_col))
        weight_col
    else NULL
    config$positive_facets <- positive_facets
    config$facet_signs <- facet_signs
    constraint_specs <- constraint_specs %||% prepare_constraint_specs(prep = prep, anchor_df = anchor_df,
        group_anchor_df = group_anchor_df, noncenter_facet = noncenter_facet, dummy_facets = dummy_facets,
        centering_scope = centering_scope, anchor_review = anchor_review)
    theta_spec <- constraint_specs$theta_spec
    if (identical(method, "MML") && (any(is.finite(theta_spec$anchors)) || any(!is.na(theta_spec$groups)))) {
        stop("Person anchors and Person group anchors are not supported with MML. ", "MML integrates Person ability rather than fixing individual values or ",
            "their group mean. Common Persons in the rating design do not require ", "anchors. Fixed Person constraints are available with JML, whose ",
            "estimation assumptions differ.", call. = FALSE)
    }
    config$theta_spec <- constraint_specs$theta_spec
    config$facet_specs <- constraint_specs$facet_specs
    config$noncenter_facet <- noncenter_facet
    config$centering_scope <- centering_scope
    config$dummy_facets <- dummy_facets
    config$anchor_summary <- constraint_specs$anchor_summary
    config$anchor_review <- constraint_specs$anchor_review
    config$step_specs <- mfrmr_default_step_specs(config)
    if (!is.null(step_constraints)) {
        config$step_constraint_map <- build_mfrm_step_map(step_constraints, config, prep)
    }
    config$source_columns <- prep$source_columns
    config$repeated_ratings <- prep$repeated_ratings %||% "review"
    config$population_spec <- compact_population_spec(population, prep$levels$Person)
    config$gpcm_spec <- if (identical(model, "GPCM")) {
        build_gpcm_slope_spec(levels = prep$levels[[slope_facet]], slope_facet = slope_facet, step_facet = step_facet)
    }
    else {
        NULL
    }
    list(config = config, sizes = build_param_sizes(config))
}



safe_cor <- function(x, y, w = NULL) {
    ok <- is.finite(x) & is.finite(y)
    if (is.null(w)) {
        if (!any(ok))
            return(NA_real_)
        x <- x[ok]
        y <- y[ok]
        if (length(unique(x)) < 2 || length(unique(y)) < 2) {
            return(NA_real_)
        }
        return(suppressWarnings(stats::cor(x, y, use = "complete.obs")))
    }
    ok <- ok & is.finite(w) & w > 0
    if (!any(ok))
        return(NA_real_)
    x <- x[ok]
    y <- y[ok]
    w <- w[ok]
    w_sum <- sum(w)
    if (w_sum <= 0)
        return(NA_real_)
    mx <- sum(w * x)/w_sum
    my <- sum(w * y)/w_sum
    vx <- sum(w * (x - mx)^2)/w_sum
    vy <- sum(w * (y - my)^2)/w_sum
    if (vx <= 0 || vy <= 0)
        return(NA_real_)
    cov <- sum(w * (x - mx) * (y - my))/w_sum
    cov/sqrt(vx * vy)
}

weighted_mean_safe <- function(x, w) {
    weighted_mean(x, w)
}

infer_default_rater_facet <- function(facet_names, fallback_first = TRUE) {
    facet_names <- as.character(facet_names)
    if (length(facet_names) == 0)
        return(NULL)
    lowered <- tolower(facet_names)
    patterns <- c("rater", "judge", "grader", "reader", "scorer", "assessor", "evaluator")
    hits <- vapply(lowered, function(x) any(vapply(patterns, function(p) grepl(p, x, fixed = TRUE), logical(1))),
        logical(1))
    if (any(hits)) {
        return(facet_names[which(hits)[1]])
    }
    if (isTRUE(fallback_first))
        facet_names[1]
    else NULL
}

calc_interrater_agreement <- function(obs_df, facet_cols, rater_facet, res = NULL) {
    if (is.null(obs_df) || nrow(obs_df) == 0) {
        return(list(summary = tibble::tibble(), pairs = tibble::tibble()))
    }
    if (is.null(rater_facet) || !rater_facet %in% facet_cols) {
        return(list(summary = tibble::tibble(), pairs = tibble::tibble()))
    }
    context_cols <- setdiff(facet_cols, rater_facet)
    if (length(context_cols) == 0) {
        return(list(summary = tibble::tibble(), pairs = tibble::tibble()))
    }
    obs_df$.context <- dplyr::group_indices(dplyr::group_by(obs_df, dplyr::across(dplyr::all_of(context_cols))))
    df <- dplyr::select(obs_df, .context, !!rlang::sym(rater_facet), Observed, dplyr::any_of("Weight"))
    df$.Weight <- get_weights(df)
    df <- dplyr::summarize(dplyr::group_by(df, .context, !!rlang::sym(rater_facet)), Score = if (all(is.finite(Observed)))
        weighted_mean(Observed, .Weight)
    else NA_real_, Ratings = dplyr::n(), .groups = "drop")
    if (nrow(df) == 0) {
        return(list(summary = tibble::tibble(), pairs = tibble::tibble()))
    }
    rater_chr <- ifelse(is.na(df[[rater_facet]]), "NA", as.character(df[[rater_facet]]))
    rater_cols <- unique(rater_chr)
    n_raters <- length(rater_cols)
    if (n_raters < 2L) {
        return(list(summary = tibble::tibble(), pairs = tibble::tibble()))
    }
    prob_map <- new.env(parent = emptyenv(), hash = TRUE)
    prob_map_size <- 0L
    if (!is.null(res))
        stop("Use interrater_agreement(fit, ...) for model-expected agreement.", call. = FALSE)
    r_id <- match(rater_chr, rater_cols)
    fin_mask <- is.finite(df$Score)
    r_fin <- r_id[fin_mask]
    ctx_fin <- as.integer(df$.context[fin_mask])
    sc_fin <- df$Score[fin_mask]
    ctx_tab <- if (length(ctx_fin) > 0L)
        tabulate(ctx_fin)
    else integer(0)
    contexts_with_pairs <- sum(ctx_tab >= 2L)
    pair_mat_all <- utils::combn(n_raters, 2L)
    n_all_pairs <- ncol(pair_mat_all)
    r1_names <- rater_cols[pair_mat_all[1L, ]]
    r2_names <- rater_cols[pair_mat_all[2L, ]]
    n_vec <- integer(n_all_pairs)
    exact_cnt_vec <- integer(n_all_pairs)
    adj_cnt_vec <- integer(n_all_pairs)
    exp_exact_cnt_vec <- rep(NA_real_, n_all_pairs)
    exp_avail_vec <- integer(n_all_pairs)
    exp_unavail_vec <- integer(n_all_pairs)
    exp_reason_vec <- rep("no_matched_scores", n_all_pairs)
    exact_vec <- rep(NA_real_, n_all_pairs)
    exp_exact_vec <- rep(NA_real_, n_all_pairs)
    adj_vec <- rep(NA_real_, n_all_pairs)
    mean_diff_vec <- rep(NA_real_, n_all_pairs)
    mad_vec <- rep(NA_real_, n_all_pairs)
    corr_vec <- rep(NA_real_, n_all_pairs)
    if (contexts_with_pairs > 0L) {
        keep_multi <- ctx_tab[ctx_fin] >= 2L
        ctx_m <- ctx_fin[keep_multi]
        r_m <- r_fin[keep_multi]
        sc_m <- sc_fin[keep_multi]
        if (is.unsorted(ctx_m)) {
            ord <- order(ctx_m)
            ctx_m <- ctx_m[ord]
            r_m <- r_m[ord]
            sc_m <- sc_m[ord]
        }
        n_m <- length(ctx_m)
        run_ends <- c(which(ctx_m[-1L] != ctx_m[-n_m]), n_m)
        run_starts <- c(1L, run_ends[-length(run_ends)] + 1L)
        run_lens <- run_ends - run_starts + 1L
        total_co <- as.integer(sum(as.double(run_lens) * as.double(run_lens - 1L)/2))
        co_pid <- integer(total_co)
        co_ctx <- integer(total_co)
        co_v1 <- numeric(total_co)
        co_v2 <- numeric(total_co)
        pos <- 0L
        for (g in seq_along(run_starts)) {
            s_idx <- run_starts[g]
            m_len <- run_lens[g]
            c_val <- ctx_m[s_idx]
            for (a in seq_len(m_len - 1L)) {
                ia <- s_idx + a - 1L
                ra <- r_m[ia]
                va <- sc_m[ia]
                for (b in (a + 1L):m_len) {
                  ib <- s_idx + b - 1L
                  rb <- r_m[ib]
                  vb <- sc_m[ib]
                  pos <- pos + 1L
                  co_ctx[pos] <- c_val
                  if (ra < rb) {
                    co_pid[pos] <- as.integer((ra - 1) * n_raters - ((ra - 1) * ra)/2 + (rb - ra))
                    co_v1[pos] <- va
                    co_v2[pos] <- vb
                  }
                  else {
                    co_pid[pos] <- as.integer((rb - 1) * n_raters - ((rb - 1) * rb)/2 + (ra - rb))
                    co_v1[pos] <- vb
                    co_v2[pos] <- va
                  }
                }
            }
        }
        by_pid <- split(seq_len(pos), co_pid)
        for (nm in names(by_pid)) {
            pid <- as.integer(nm)
            rows_p <- by_pid[[nm]]
            n_ok <- length(rows_p)
            v1 <- co_v1[rows_p]
            v2 <- co_v2[rows_p]
            diff <- v1 - v2
            exact_count <- sum(diff == 0, na.rm = TRUE)
            contexts <- as.character(co_ctx[rows_p])
            expected_reasons <- rep("fitted_probabilities_unavailable", n_ok)
            exp_vals <- if (prob_map_size > 0L && n_ok > 0L) {
                tmp <- numeric(n_ok)
                kept <- 0L
                key1_vec <- paste(contexts, r1_names[pid], sep = "||")
                key2_vec <- paste(contexts, r2_names[pid], sep = "||")
                for (k in seq_len(n_ok)) {
                  p1 <- get0(key1_vec[k], envir = prob_map, inherits = FALSE)
                  p2 <- get0(key2_vec[k], envir = prob_map, inherits = FALSE)
                  if (is.null(p1) || is.null(p2))
                    next
                  if (!p1$primary_available || !p2$primary_available) {
                    expected_reasons[k] <- "primary_estimate_unavailable"
                    next
                  }
                  if (p1$repeated || p2$repeated) {
                    expected_reasons[k] <- "repeated_ratings"
                    next
                  }
                  p1 <- p1$probability
                  p2 <- p2$probability
                  if (any(!is.finite(p1)) || any(!is.finite(p2)) || any(p1 < 0 | p1 > 1) || any(p2 < 0 | p2 >
                    1) || abs(sum(p1) - 1) > 1e-08 || abs(sum(p2) - 1) > 1e-08)
                    next
                  kept <- kept + 1L
                  tmp[kept] <- sum(p1 * p2)
                  expected_reasons[k] <- "available"
                }
                if (kept > 0L)
                  tmp[seq_len(kept)]
                else numeric(0)
            }
            else {
                numeric(0)
            }
            n_exp <- length(exp_vals)
            n_vec[pid] <- n_ok
            exact_cnt_vec[pid] <- exact_count
            adj_cnt_vec[pid] <- sum(abs(diff) <= 1, na.rm = TRUE)
            exp_exact_cnt_vec[pid] <- if (n_exp == n_ok)
                sum(exp_vals)
            else NA_real_
            exp_avail_vec[pid] <- n_exp
            exp_unavail_vec[pid] <- n_ok - n_exp
            exp_reason_vec[pid] <- if (n_exp == n_ok)
                "available"
            else paste(sort(unique(expected_reasons[expected_reasons != "available"])), collapse = ";")
            exact_vec[pid] <- exact_count/n_ok
            exp_exact_vec[pid] <- if (n_exp == n_ok)
                mean(exp_vals)
            else NA_real_
            adj_vec[pid] <- mean(abs(diff) <= 1, na.rm = TRUE)
            mean_diff_vec[pid] <- mean(diff, na.rm = TRUE)
            mad_vec[pid] <- mean(abs(diff), na.rm = TRUE)
            corr_vec[pid] <- safe_cor(v1, v2)
        }
    }
    count_cast <- if (any(n_vec == 0L))
        as.numeric
    else as.integer
    pair_tbl <- tibble::tibble(Rater1 = r1_names, Rater2 = r2_names, N = count_cast(n_vec), OpportunityCount = count_cast(n_vec),
        ExactCount = count_cast(exact_cnt_vec), ExpectedExactCount = exp_exact_cnt_vec, ExpectedAvailableContexts = exp_avail_vec,
        ExpectedUnavailableContexts = exp_unavail_vec, ExpectedReason = exp_reason_vec, AdjacentCount = count_cast(adj_cnt_vec),
        Exact = exact_vec, ExpectedExact = exp_exact_vec, Adjacent = adj_vec, MeanDiff = mean_diff_vec, MAD = mad_vec,
        Corr = corr_vec)
    total_pairs <- sum(pair_tbl$N, na.rm = TRUE)
    total_exact <- sum(pair_tbl$Exact * pair_tbl$N, na.rm = TRUE)
    expected_available <- total_pairs > 0 && all(is.finite(pair_tbl$ExpectedExact[pair_tbl$N > 0]))
    total_expected <- if (expected_available) {
        sum(pair_tbl$ExpectedExact * pair_tbl$N, na.rm = TRUE)
    }
    else {
        NA_real_
    }
    summary_tbl <- tibble::tibble(RaterFacet = rater_facet, Raters = length(rater_cols), Pairs = nrow(pair_tbl), AvailablePairs = sum(pair_tbl$N >
        0), UnavailablePairs = sum(pair_tbl$N == 0), CorrelationPairs = sum(is.finite(pair_tbl$Corr)), UnavailableCorrelationPairs = sum(!is.finite(pair_tbl$Corr)),
        RepeatedCells = sum(df$Ratings > 1L), ExpectedAvailableContexts = sum(pair_tbl$ExpectedAvailableContexts),
        ExpectedUnavailableContexts = sum(pair_tbl$ExpectedUnavailableContexts), ExpectedReason = if (expected_available)
            "available"
        else if (total_pairs == 0)
            "no_matched_scores"
        else paste(sort(unique(unlist(strsplit(pair_tbl$ExpectedReason[pair_tbl$N > 0 & !is.finite(pair_tbl$ExpectedExact)],
            ";", fixed = TRUE)))), collapse = ";"), Contexts = contexts_with_pairs, TotalPairs = total_pairs,
        OpportunityCount = total_pairs, ExactAgreements = total_exact, ExpectedAgreements = ifelse(expected_available,
            total_expected, NA_real_), ExactAgreement = ifelse(total_pairs > 0, total_exact/total_pairs, NA_real_),
        ExpectedExactAgreement = ifelse(expected_available && total_pairs > 0, total_expected/total_pairs,
            NA_real_), AgreementMinusExpected = ifelse(expected_available && total_pairs > 0, (total_exact -
            total_expected)/total_pairs, NA_real_), AdjacentAgreements = sum(pair_tbl$AdjacentCount, na.rm = TRUE),
        AdjacentAgreement = weighted_mean_safe(pair_tbl$Adjacent, pair_tbl$N), MeanAbsDiff = weighted_mean_safe(pair_tbl$MAD,
            pair_tbl$N), MeanCorr = weighted_mean_safe(pair_tbl$Corr, pair_tbl$N))
    list(summary = summary_tbl, pairs = pair_tbl)
}

build_param_slices <- function(sizes) {
    out <- list()
    idx <- 1L
    for (nm in names(sizes)) {
        k <- as.integer(sizes[[nm]] %||% 0L)
        if (k > 0L) {
            out[[nm]] <- idx:(idx + k - 1L)
            idx <- idx + k
        }
        else {
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
    for (i in seq_len(n_levels)) {
        unit <- numeric(n_levels)
        unit[i] <- 1
        out[i, ] <- constraint_grad_project(unit, spec)
    }
    out
}

residual_pca_correlation <- function(residual_matrix) {
    mat <- as.matrix(residual_matrix)
    mat[!is.finite(mat)] <- NA_real_
    correlation <- tryCatch(suppressWarnings(stats::cor(mat, use = "pairwise.complete.obs")), error = function(e) NULL)
    error <- NULL
    if (is.null(correlation) || any(!is.finite(correlation))) {
        error <- paste("Residual correlations are unavailable for one or more columns or pairs.", "Check shared-person counts and residual variation; missing correlations are not replaced with zero.")
    }
    else {
        values <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
        if (min(values) < -sqrt(.Machine$double.eps) * max(1, max(values))) {
            error <- paste("Pairwise residual correlations do not form a positive semidefinite matrix.",
                "Review missingness and shared-person coverage; the matrix is not automatically smoothed.")
        }
    }
    list(cor_matrix = correlation, error = error)
}

compute_principal_components <- function(cor_matrix, n_factors) {
    result <- tryCatch({
        decomposition <- eigen(cor_matrix, symmetric = TRUE)
        loadings <- sweep(decomposition$vectors[, seq_len(n_factors), drop = FALSE],
            2L, sqrt(pmax(0, decomposition$values[seq_len(n_factors)])), "*")
        orientation <- ifelse(colSums(loadings) < 0, -1, 1)
        loadings <- sweep(loadings, 2L, orientation, "*")
        dimnames(loadings) <- list(colnames(cor_matrix), paste0("PC", seq_len(n_factors)))
        list(values = decomposition$values, loadings = structure(loadings, class = "loadings"))
    }, error = identity)
    list(pca = result, warning = NULL)
}

compute_pca_overall <- function(obs_df, facet_names, max_factors = 10L) {
    pca_failure <- function(message, residual_matrix = NULL, cor_matrix = NULL, warning = NULL) {
        list(pca = NULL, residual_matrix = residual_matrix, cor_matrix = cor_matrix, error = message,
            warning = warning, calculation_version = 2L, max_factors = as.integer(max_factors))
    }
    if (length(facet_names) == 0)
        return(NULL)
    combinations <- unique(obs_df[, facet_names, drop = FALSE])
    if (anyDuplicated(do.call(paste, c(combinations, sep = "_")))) {
        return(pca_failure("Combined-facet column labels are ambiguous. Relabel facet levels to distinguish the combinations."))
    }
    df_aug <- dplyr::select(dplyr::mutate(obs_df, Person = as.character(Person), item_combination = paste(!!!rlang::syms(facet_names),
        sep = "_")), Person, item_combination, StdResidual)
    residual_matrix_prep <- dplyr::summarize(dplyr::group_by(df_aug, Person, item_combination), std_residual = mean(StdResidual,
        na.rm = TRUE), .groups = "drop")
    residual_matrix_wide <- tryCatch({
        tibble::column_to_rownames(tidyr::pivot_wider(residual_matrix_prep, id_cols = Person, names_from = item_combination,
            values_from = std_residual, values_fill = list(std_residual = NA)), "Person")
    }, error = function(e) e)
    if (inherits(residual_matrix_wide, "error")) {
        return(pca_failure(paste0("Could not reshape the residual matrix: ", conditionMessage(residual_matrix_wide))))
    }
    if (is.null(residual_matrix_wide) || nrow(residual_matrix_wide) < 2 || ncol(residual_matrix_wide) <
        2) {
        return(pca_failure("Residual PCA requires at least two persons and two combined-facet columns."))
    }
    correlation <- residual_pca_correlation(residual_matrix_wide)
    cor_matrix <- correlation$cor_matrix
    if (!is.null(correlation$error)) {
        return(pca_failure(correlation$error, residual_matrix_wide, cor_matrix))
    }
    max_cap <- if (is.na(max_factors) || !is.finite(max_factors)) {
        min(10L, ncol(cor_matrix) - 1L, nrow(cor_matrix) - 1L)
    }
    else {
        as.integer(max_factors)
    }
    n_factors <- max(1, min(max_cap, ncol(cor_matrix) - 1, nrow(cor_matrix) - 1))
    pca_run <- compute_principal_components(cor_matrix, n_factors)
    pca_result <- pca_run$pca
    if (inherits(pca_result, "error")) {
        return(pca_failure(paste0("Principal components could not be computed: ", conditionMessage(pca_result)),
            residual_matrix = residual_matrix_wide, cor_matrix = cor_matrix, warning = pca_run$warning))
    }
    list(pca = pca_result, residual_matrix = residual_matrix_wide, cor_matrix = cor_matrix, error = NULL,
        warning = pca_run$warning, calculation_version = 2L, max_factors = as.integer(max_factors))
}

compute_pca_by_facet <- function(obs_df, facet_names, max_factors = 10L) {
    out <- list()
    for (facet in facet_names) {
        pca_failure <- function(message, residual_matrix = NULL, cor_matrix = NULL, warning = NULL) {
            list(pca = NULL, residual_matrix = residual_matrix, cor_matrix = cor_matrix, error = message,
                warning = warning, calculation_version = 2L, max_factors = as.integer(max_factors))
        }
        facet_sym <- rlang::sym(facet)
        prep <- dplyr::summarize(dplyr::group_by(dplyr::select(dplyr::mutate(obs_df, Person = as.character(Person), .Level = as.character(!!facet_sym)),
            Person, .Level, StdResidual), Person, .Level), std_residual = mean(StdResidual, na.rm = TRUE),
            .groups = "drop")
        wide <- tryCatch({
            tibble::column_to_rownames(tidyr::pivot_wider(prep, id_cols = Person, names_from = .Level,
                values_from = std_residual, values_fill = list(std_residual = NA)), "Person")
        }, error = function(e) e)
        if (inherits(wide, "error")) {
            out[[facet]] <- pca_failure(paste0("Could not reshape the residual matrix: ", conditionMessage(wide)))
            next
        }
        if (is.null(wide) || nrow(wide) < 2 || ncol(wide) < 2) {
            out[[facet]] <- pca_failure("Residual PCA requires at least two persons and two facet-level columns.")
            next
        }
        correlation <- residual_pca_correlation(wide)
        cor_mat <- correlation$cor_matrix
        if (!is.null(correlation$error)) {
            out[[facet]] <- pca_failure(correlation$error, wide, cor_mat)
            next
        }
        max_cap <- if (is.na(max_factors) || !is.finite(max_factors)) {
            min(10L, ncol(cor_mat) - 1L, nrow(cor_mat) - 1L)
        }
        else {
            as.integer(max_factors)
        }
        n_factors <- max(1, min(max_cap, ncol(cor_mat) - 1, nrow(cor_mat) - 1))
        pca_run <- compute_principal_components(cor_mat, n_factors)
        pca_obj <- pca_run$pca
        if (inherits(pca_obj, "error")) {
            out[[facet]] <- pca_failure(paste0("Principal components could not be computed: ", conditionMessage(pca_obj)),
                residual_matrix = wide, cor_matrix = cor_mat, warning = pca_run$warning)
            next
        }
        out[[facet]] <- list(pca = pca_obj, cor_matrix = cor_mat, residual_matrix = wide, error = NULL,
            warning = pca_run$warning, calculation_version = 2L, max_factors = as.integer(max_factors))
    }
    out
}
