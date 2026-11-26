
.mfrm_fit_spec_keys <- function() c(c("person", "facets", "score", "model", "method", "step_facet", "slope_facet",
    "slope_reference", "slope_type", "gpcm_mml_identification", "rating_min", "rating_max", "weight", "anchors",
    "group_anchors", "dummy_facets", "noncenter_facet", "positive_facets", "centering_scope", "population_formula",
    "person_data", "person_id", "population_sd", "step_constraints", "category_support", "facet_formulas",
    "facet_data", "facet_interactions"), "missing_codes")

.mfrm_fit_control_keys <- function() c("quad_points", "mml_integration", "maxit", "reltol", "facet_shrinkage",
"facet_prior_sd", "slope_prior_sd", "jml_penalty", "person_score",
"xtreme", "whexact", "bias_facets", "warn_connectivity")

.mfrm_fit_defaults <- function() list(data = NULL, person = NULL, facets = NULL, score = NULL, model = "RSM",
    method = "MML", step_facet = NULL, slope_facet = NULL, slope_reference = NULL, slope_type = "full_predictor",
    gpcm_mml_identification = "free_population", person_score = "auto", rating_min = NULL, rating_max = NULL,
    weight = NULL, anchors = NULL, group_anchors = NULL, dummy_facets = NULL, noncenter_facet = "Person",
    positive_facets = NULL, centering_scope = "facet", population_formula = NULL, person_data = NULL, person_id = NULL,
    population_sd = 1, step_constraints = NULL, category_support = NULL, facet_formulas = NULL, facet_data = NULL,
    facet_interactions = NULL, facet_shrinkage = "none", facet_prior_sd = NULL, slope_prior_sd = NULL, quad_points = 31L,
    mml_integration = "fixed", xtreme = 0, jml_penalty = "none", whexact = FALSE, bias_facets = NULL, warn_connectivity = TRUE,
    maxit = 400L, reltol = 1e-09, missing_codes = NULL)

.mfrm_fit_names <- function() {
    names(.mfrm_fit_defaults())
}

.mfrm_argument_hint <- function(bad_names) {
    canonical <- setdiff(.mfrm_fit_names(), "data")
    spec_keys <- .mfrm_fit_spec_keys()
    hints <- character(0)
    for (nm in bad_names) {
        if (!nzchar(nm))
            next
        d <- utils::adist(nm, canonical, ignore.case = TRUE)[1L, ]
        best_i <- which.min(d)
        if (length(best_i) == 1L && is.finite(d[best_i]) && d[best_i] <= 3) {
            cand <- canonical[best_i]
            bucket <- if (cand %in% spec_keys)
                "`spec = list(...)`"
            else "`control = list(...)`"
            hints <- c(hints, sprintf("Did you mean `%s` (in %s)?", cand, bucket))
        }
    }
    if (length(hints))
        paste0(" ", paste(unique(hints), collapse = " "))
    else ""
}

.mfrm_validate_spec_list <- function(spec, caller = "fit_mfrm()") {
    if (is.null(spec))
        return(list())
    if (!is.list(spec) || is.data.frame(spec)) {
        stop("`spec` in ", caller, " must be NULL or a named list of measurement-model specification settings.",
            call. = FALSE)
    }
    if (!length(spec))
        return(list())
    nm <- names(spec)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
        stop("`spec` in ", caller, " must have unique non-empty names.", call. = FALSE)
    }
    extra <- setdiff(nm, .mfrm_fit_spec_keys())
    if (length(extra) > 0L) {
        in_ctrl <- intersect(extra, .mfrm_fit_control_keys())
        if (length(in_ctrl) > 0L) {
            stop("`spec` in ", caller, " contains estimation/numerical control entries (", paste(in_ctrl,
                collapse = ", "), "); pass them in `control = list(...)` instead of `spec`.", call. = FALSE)
        }
        stop("Unsupported or unknown entries in `spec`: ", paste(extra, collapse = ", "), ".", .mfrm_argument_hint(extra),
            call. = FALSE)
    }
    if (!identical(caller, "fit_mfrm()") && !is.null(spec$step_constraints)) {
        stop(caller, " does not support `step_constraints`. Use fit_mfrm().", call. = FALSE)
    }
    spec
}

.mfrm_validate_control_list <- function(control, caller = "fit_mfrm()") {
    if (is.null(control))
        return(list())
    if (!is.list(control) || is.data.frame(control)) {
        stop("`control` in ", caller, " must be a named list of estimation/numerical control settings.",
            call. = FALSE)
    }
    if (!length(control))
        return(list())
    nm <- names(control)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
        stop("`control` in ", caller, " must have unique non-empty names.", call. = FALSE)
    }
    extra <- setdiff(nm, .mfrm_fit_control_keys())
    if (length(extra) > 0L) {
        in_spec <- intersect(extra, .mfrm_fit_spec_keys())
        if (length(in_spec) > 0L) {
            stop("`control` in ", caller, " contains model/data specification entries (", paste(in_spec,
                collapse = ", "), "); pass them as top-level arguments or in `spec = list(...)` instead of `control`.",
                call. = FALSE)
        }
        stop("Unsupported or unknown entries in `control`: ", paste(extra, collapse = ", "), ".", .mfrm_argument_hint(extra),
            call. = FALSE)
    }
    control
}

fit_mfrm <- function(data, person = NULL, facets = NULL, score = NULL, model = c("RSM", "PCM", "GPCM"), method = c("MML",
    "JML", "JMLE", "JML_jeffreys", "CML"), step_facet = NULL, slope_facet = NULL, rating_min = NULL, rating_max = NULL,
    weight = NULL, anchors = NULL, group_anchors = NULL, population_formula = NULL, person_data = NULL, facet_interactions = NULL,
    step_constraints = NULL, category_support = NULL, facet_formulas = NULL, facet_data = NULL, person_score = c("auto",
        "EAP", "WLE"), warn_connectivity = TRUE, spec = NULL, control = list(), ...) {
    mc <- match.call(expand.dots = FALSE)
    formal_supplied <- setdiff(names(as.list(mc)[-1L]), c("...", "spec", "control"))
    spec_list <- .mfrm_validate_spec_list(spec, "fit_mfrm()")
    ctrl_list <- .mfrm_validate_control_list(control, "fit_mfrm()")
    core_formal_names <- c("data", "person", "facets", "score", "model", "method", "step_facet", "slope_facet",
        "rating_min", "rating_max", "weight", "anchors", "group_anchors", "population_formula", "person_data",
        "facet_interactions", "step_constraints", "category_support", "facet_formulas", "facet_data", "person_score",
        character(0), "warn_connectivity")
    advanced_names <- setdiff(.mfrm_fit_names(), core_formal_names)
    dot_len <- ...length()
    dots <- if (dot_len > 0L) {
        dot_nms <- ...names()
        if (is.null(dot_nms) || length(dot_nms) != dot_len || anyNA(dot_nms) || any(!nzchar(dot_nms))) {
            stop("Unnamed extra arguments in `...` are not supported in `fit_mfrm()`. ", "Supply core arguments by name, or pass advanced settings via `spec = list(...)` or `control = list(...)`.",
                call. = FALSE)
        }
        if (anyDuplicated(dot_nms)) {
            stop("Duplicate argument(s) in `...`: ", paste(unique(dot_nms[duplicated(dot_nms)]), collapse = ", "),
                ".", call. = FALSE)
        }
        bad_dots <- setdiff(dot_nms, advanced_names)
        if (length(bad_dots) > 0L) {
            stop("Unknown argument(s) in `fit_mfrm()`: ", paste(bad_dots, collapse = ", "), ".", .mfrm_argument_hint(bad_dots),
                call. = FALSE)
        }
        list(...)
    }
    else {
        list()
    }
    dup_spec_dots <- intersect(names(spec_list), names(dots))
    if (length(dup_spec_dots) > 0L) {
        stop("Duplicate specification in both `spec` and `...`: ", paste(dup_spec_dots, collapse = ", "),
            ". Specify each setting in one place.", call. = FALSE)
    }
    dup_ctrl_dots <- intersect(names(ctrl_list), names(dots))
    if (length(dup_ctrl_dots) > 0L) {
        stop("Duplicate specification in both `control` and `...`: ", paste(dup_ctrl_dots, collapse = ", "),
            ". Specify each control setting in one place.", call. = FALSE)
    }
    resolved <- .mfrm_fit_defaults()
    if (length(spec_list) > 0L)
        resolved[names(spec_list)] <- spec_list
    if (length(ctrl_list) > 0L)
        resolved[names(ctrl_list)] <- ctrl_list
    if (length(dots) > 0L)
        resolved[names(dots)] <- dots
    for (nm in setdiff(formal_supplied, "data")) {
        resolved[nm] <- list(get(nm, envir = environment(), inherits = FALSE))
    }
    resolved$data <- data
    resolved$model <- .mfrm_model_code(resolved$model)
    missing_codes <- resolved$missing_codes
    if (!is.null(missing_codes) && !isFALSE(missing_codes)) {
        defaults <- isTRUE(missing_codes) || identical(tolower(as.character(missing_codes)), "default")
        resolved$data <- recode_missing_codes(data, columns = if (defaults)
            resolved$score
        else c(resolved$person, resolved$facets, resolved$score), codes = if (defaults)
            c("99", "999", "-1", "N", "NA", "n/a", ".", "")
        else missing_codes)
    }
    resolved$missing_codes <- NULL
    do.call(.mfrm_estimate, resolved)
}

prepare_mfrm_population_scaffold <- function(data, person, population_formula = NULL, person_data = NULL,
    person_id = NULL, population_policy = c("error", "omit"), population_xlevels = NULL, population_contrasts = NULL,
    require_full_rank = TRUE, population_terms = NULL) {
    if (is.null(population_formula)) {
        return(list(active = FALSE, formula = NULL, person_id = if (is.null(person_id)) person else person_id,
            person_table = NULL, person_table_replay = NULL, person_table_replay_scope = NULL, design_matrix = NULL,
            design_columns = NULL, xlevels = NULL, contrasts = NULL, coefficients = NULL, sigma2 = NULL,
            converged = FALSE, posterior_basis = "legacy_mml", policy = NULL, included_persons = character(0),
            omitted_persons = character(0), response_rows_retained = nrow(data), response_rows_omitted = 0L,
            notes = "No population model was requested; this fit uses an unconditional normal person distribution."))
    }
    population_policy <- match.arg(population_policy)
    if (inherits(population_formula, "terms")) {
        population_contrasts <- population_contrasts %||% attr(population_formula, "mfrm_population_contrasts",
            exact = TRUE)
        population_xlevels <- population_xlevels %||% attr(population_formula, "mfrm_population_xlevels",
            exact = TRUE)
    }
    population_formula <- tryCatch(stats::as.formula(population_formula), error = function(e) {
        stop("`population_formula` must be coercible to a one-sided formula such as `~ grade + ses`. ",
            "Original error: ", conditionMessage(e), call. = FALSE)
    })
    if (length(population_formula) != 2L) {
        stop("`population_formula` must be a one-sided formula such as `~ grade + ses`.", call. = FALSE)
    }
    if (!is.data.frame(person_data)) {
        stop("`person_data` must be a data.frame with one row per person when `population_formula` is supplied.",
            call. = FALSE)
    }
    if (is.null(person_id)) {
        if (!person %in% names(person_data)) {
            stop("`person_data` must contain the person column `", person, "` or you must supply `person_id` explicitly when `population_formula` is used.",
                call. = FALSE)
        }
        person_id <- person
    }
    if (!person_id %in% names(person_data)) {
        stop("`person_id = '", person_id, "'` is not a column in `person_data`.", call. = FALSE)
    }
    observed_persons <- unique(as.character(data[[person]]))
    person_tbl <- as.data.frame(person_data, stringsAsFactors = FALSE)
    person_tbl[[person_id]] <- as.character(person_tbl[[person_id]])
    if (anyNA(person_tbl[[person_id]]) || any(!nzchar(person_tbl[[person_id]]))) {
        stop("`person_data` contains missing or empty person IDs in column `", person_id, "`.", call. = FALSE)
    }
    dup_ids <- unique(person_tbl[[person_id]][duplicated(person_tbl[[person_id]])])
    if (length(dup_ids) > 0) {
        preview <- paste(utils::head(dup_ids, 5), collapse = ", ")
        stop("`person_data` must contain one unique row per person. Duplicate IDs detected in `", person_id,
            "`: ", preview, if (length(dup_ids) > 5)
                ", ..."
            else ".", call. = FALSE)
    }
    missing_persons <- setdiff(observed_persons, person_tbl[[person_id]])
    if (length(missing_persons) > 0) {
        preview <- paste(utils::head(missing_persons, 5), collapse = ", ")
        stop("`person_data` is missing ", length(missing_persons), " person(s) observed in `data`: ",
            preview, if (length(missing_persons) > 5)
                ", ..."
            else ".", call. = FALSE)
    }
    keep_idx <- match(observed_persons, person_tbl[[person_id]])
    person_tbl <- person_tbl[keep_idx, , drop = FALSE]
    rownames(person_tbl) <- NULL
    person_tbl_replay <- person_tbl
    model_frame_args <- list(formula = population_terms %||% population_formula, data = person_tbl, na.action = stats::na.pass)
    if (!is.null(population_xlevels)) {
        model_frame_args$xlev <- population_xlevels
    }
    mf <- tryCatch(do.call(stats::model.frame, model_frame_args), error = function(e) {
        stop("Could not build the latent-regression model frame from `person_data`. ", "Original error: ",
            conditionMessage(e), call. = FALSE)
    })
    mf_names <- names(mf)
    if (length(mf_names) > 0) {
        unsupported <- mf_names[vapply(mf[mf_names], function(x) is.list(x) || is.data.frame(x), logical(1))]
        if (length(unsupported) > 0) {
            stop("Variables referenced in `population_formula` must be atomic columns usable by stats::model.matrix(). ",
                "Numeric, logical, factor, ordered factor, and character predictors are supported; unsupported columns: ",
                paste(unsupported, collapse = ", "), ".", call. = FALSE)
        }
    }
    complete_mask <- if (ncol(mf) == 0)
        rep(TRUE, nrow(person_tbl))
    else stats::complete.cases(mf)
    omitted_persons <- person_tbl[[person_id]][!complete_mask]
    if (length(omitted_persons) > 0 && identical(population_policy, "error")) {
        preview <- paste(utils::head(omitted_persons, 5), collapse = ", ")
        stop("`person_data` contains missing covariate values for persons referenced by `population_formula`. ",
            "Use `population_policy = 'omit'` to fit the complete-case analysis set instead. Affected IDs: ",
            preview, if (length(omitted_persons) > 5)
                ", ..."
            else ".", call. = FALSE)
    }
    if (length(omitted_persons) > 0) {
        person_tbl <- person_tbl[complete_mask, , drop = FALSE]
        mf <- mf[complete_mask, , drop = FALSE]
        rownames(person_tbl) <- NULL
        rownames(mf) <- NULL
    }
    if (nrow(person_tbl) == 0) {
        stop("No persons remain in `person_data` after applying the latent-regression covariate policy.",
            call. = FALSE)
    }
    included_persons <- as.character(person_tbl[[person_id]])
    response_keep <- as.character(data[[person]]) %in% included_persons
    response_rows_retained <- sum(response_keep)
    response_rows_omitted <- sum(!response_keep)
    terms_obj <- attr(mf, "terms") %||% stats::terms(population_formula)
    mm <- tryCatch(stats::model.matrix(terms_obj, data = mf, contrasts.arg = population_contrasts), error = function(e) {
        stop("Could not build the latent-regression design matrix from `person_data`. ", "Original error: ",
            conditionMessage(e), call. = FALSE)
    })
    xlevels <- tryCatch(stats::.getXlevels(terms_obj, mf), error = function(e) NULL)
    mm_contrasts <- attr(mm, "contrasts", exact = TRUE)
    qr_mm <- qr(mm)
    if (isTRUE(require_full_rank) && qr_mm$rank < ncol(mm)) {
        stop("The latent-regression design matrix is rank-deficient. ", "Adjust `population_formula` or the coding of `person_data` before fitting.",
            call. = FALSE)
    }
    mm_rcond <- tryCatch(suppressWarnings(rcond(mm)), error = function(e) NA_real_)
    if (isTRUE(require_full_rank) && is.finite(mm_rcond) && mm_rcond < 1e-08) {
        warning(sprintf(paste0("Latent-regression design matrix is near-singular (rcond = %.1e). ", "Estimated coefficients may be unstable; consider dropping ",
            "redundant covariates or rescaling `person_data`."), mm_rcond), call. = FALSE)
    }
    notes <- c("Latent-regression covariates are active.", "Person-level covariates and their design matrix were prepared for MML estimation.")
    if (length(omitted_persons) > 0) {
        notes <- c(notes, paste0("Complete-case preparation retained ", nrow(person_tbl), " person(s) and omitted ",
            length(omitted_persons), " person(s) under `population_policy = 'omit'`."))
    }
    if (response_rows_omitted > 0) {
        notes <- c(notes, paste0("For latent-regression estimation, ", response_rows_omitted, " response row(s) would be excluded because their persons lack complete background data."))
    }
    list(active = TRUE, formula = population_formula, terms = terms_obj, person_id = person_id, person_table = person_tbl,
        person_table_replay = person_tbl_replay, person_table_replay_scope = "observed_person_subset_pre_omit",
        design_matrix = mm, design_columns = colnames(mm), xlevels = xlevels, contrasts = mm_contrasts,
        coefficients = NULL, sigma2 = NULL, converged = FALSE, posterior_basis = "population_model",
        policy = population_policy, included_persons = included_persons, omitted_persons = as.character(omitted_persons),
        response_rows_retained = as.integer(response_rows_retained), response_rows_omitted = as.integer(response_rows_omitted),
        notes = notes)
}


normalize_expected_mfrm_design <- function(expected_design, person, facets) {
    if (is.null(expected_design))
        return(NULL)
    if (!is.data.frame(expected_design)) {
        stop("`expected_design` must be a data.frame.", call. = FALSE)
    }
    required <- c(person, facets)
    missing_cols <- setdiff(required, names(expected_design))
    if (length(missing_cols) > 0L) {
        stop("`expected_design` is missing required column(s): ", paste(missing_cols, collapse = ", "),
            ".", call. = FALSE)
    }
    if (nrow(expected_design) == 0L) {
        stop("`expected_design` must contain at least one planned rating cell.", call. = FALSE)
    }
    out <- as.data.frame(expected_design[, required, drop = FALSE], stringsAsFactors = FALSE)
    names(out)[1L] <- "Person"
    out[] <- lapply(out, function(x) trimws(as.character(x)))
    invalid <- vapply(out, function(x) any(is.na(x) | !nzchar(x)), logical(1))
    if (any(invalid)) {
        stop("`expected_design` contains missing or blank IDs in: ", paste(names(out)[invalid], collapse = ", "),
            ".", call. = FALSE)
    }
    duplicate_rows <- duplicated(out) | duplicated(out, fromLast = TRUE)
    if (any(duplicate_rows)) {
        stop("`expected_design` must contain one row per planned Person x ", paste(facets, collapse = " x "),
            " cell; duplicate planned cells were found.", call. = FALSE)
    }
    out
}

audit_mfrm_expected_design <- function(df, expected_design, facets) {
    key_cols <- c("Person", facets)
    observed <- unique(as.data.frame(df[, key_cols, drop = FALSE], stringsAsFactors = FALSE))
    empty_cells <- observed[0, , drop = FALSE]
    if (is.null(expected_design)) {
        level_coverage <- dplyr::bind_rows(lapply(facets, function(facet) {
            tibble::tibble(Facet = facet, ExpectedLevels = NA_integer_, ObservedLevels = dplyr::n_distinct(observed[[facet]]),
                ExpectedOnlyLevels = NA_integer_, UnexpectedLevels = NA_integer_)
        }))
        return(list(summary = data.frame(Status = "not_declared", ExpectedCells = NA_integer_, ObservedCells = nrow(observed),
            MatchedCells = NA_integer_, MissingExpectedCells = NA_integer_, UnexpectedObservedCells = NA_integer_,
            CoverageRate = NA_real_, ExpectedOnlyPersons = NA_integer_, UnexpectedPersons = NA_integer_, stringsAsFactors = FALSE),
            missing_expected_cells = empty_cells, unexpected_observed_cells = empty_cells, level_coverage = as.data.frame(level_coverage,
                stringsAsFactors = FALSE), settings = list(declared = FALSE, key = key_cols, note = paste0("Structural missingness was not assessed because `expected_design` was not supplied. ",
                "Absent rows cannot be distinguished from cells that were never assigned."))))
    }
    matched <- dplyr::semi_join(expected_design, observed, by = key_cols)
    missing_expected <- dplyr::anti_join(expected_design, observed, by = key_cols)
    unexpected_observed <- dplyr::anti_join(observed, expected_design, by = key_cols)
    level_coverage <- dplyr::bind_rows(lapply(facets, function(facet) {
        expected_levels <- unique(expected_design[[facet]])
        observed_levels <- unique(observed[[facet]])
        tibble::tibble(Facet = facet, ExpectedLevels = length(expected_levels), ObservedLevels = length(observed_levels),
            ExpectedOnlyLevels = length(setdiff(expected_levels, observed_levels)), UnexpectedLevels = length(setdiff(observed_levels,
                expected_levels)))
    }))
    list(summary = data.frame(Status = "declared", ExpectedCells = nrow(expected_design), ObservedCells = nrow(observed),
        MatchedCells = nrow(matched), MissingExpectedCells = nrow(missing_expected), UnexpectedObservedCells = nrow(unexpected_observed),
        CoverageRate = nrow(matched)/nrow(expected_design), ExpectedOnlyPersons = length(setdiff(expected_design$Person,
            observed$Person)), UnexpectedPersons = length(setdiff(observed$Person, expected_design$Person)),
        stringsAsFactors = FALSE), missing_expected_cells = as.data.frame(missing_expected, stringsAsFactors = FALSE),
        unexpected_observed_cells = as.data.frame(unexpected_observed, stringsAsFactors = FALSE), level_coverage = as.data.frame(level_coverage,
            stringsAsFactors = FALSE), settings = list(declared = TRUE, key = key_cols, note = "Observed unique cells were compared with the explicitly declared assignment roster."))
}

summarize_duplicate_mfrm_cells <- function(df, facets) {
    key_cols <- c("Person", facets)
    detail <- as.data.frame(dplyr::arrange(dplyr::filter(dplyr::summarise(dplyr::group_by(df, dplyr::across(dplyr::all_of(key_cols))),
        Observations = dplyr::n(), WeightedN = sum(.data$Weight, na.rm = TRUE), .groups = "drop"), .data$Observations >
        1L), dplyr::desc(.data$Observations)), stringsAsFactors = FALSE)
    list(summary = data.frame(DuplicateCells = nrow(detail), RowsInDuplicateCells = if (nrow(detail) >
        0L) sum(detail$Observations) else 0L, ExtraRows = if (nrow(detail) > 0L) sum(detail$Observations -
        1L) else 0L, HasDuplicates = nrow(detail) > 0L, stringsAsFactors = FALSE), detail = detail)
}

mfrm_bipartite_components <- function(cells, facet, basis, include_persons = FALSE, include_membership = FALSE) {
    edges <- unique(as.data.frame(cells[, c("Person", facet), drop = FALSE], stringsAsFactors = FALSE))
    edges$Person <- as.character(edges$Person)
    edges[[facet]] <- as.character(edges[[facet]])
    persons_to_levels <- split(edges[[facet]], edges$Person)
    levels_to_persons <- split(edges$Person, edges[[facet]])
    unseen <- unique(edges$Person)
    components <- list()
    membership <- if (include_membership) stats::setNames(integer(length(levels_to_persons)), names(levels_to_persons))
    while (length(unseen) > 0L) {
        component_persons <- unseen[1L]
        component_levels <- character(0)
        repeat {
            next_levels <- unique(unname(unlist(persons_to_levels[component_persons], use.names = FALSE)))
            next_persons <- unique(unname(unlist(levels_to_persons[next_levels], use.names = FALSE)))
            expanded_persons <- union(component_persons, next_persons)
            expanded_levels <- union(component_levels, next_levels)
            if (setequal(expanded_persons, component_persons) && setequal(expanded_levels, component_levels))
                break
            component_persons <- expanded_persons
            component_levels <- expanded_levels
        }
        component_edges <- sum(edges$Person %in% component_persons & edges[[facet]] %in% component_levels)
        components[[length(components) + 1L]] <- data.frame(Basis = basis, Facet = facet, Component = length(components) +
            1L, Persons = length(component_persons), FacetLevels = length(component_levels), Edges = component_edges,
            LevelLabels = paste(sort(component_levels), collapse = ", "), PersonLabels = if (isTRUE(include_persons)) {
                paste(sort(component_persons), collapse = ", ")
            }
            else {
                ""
            }, stringsAsFactors = FALSE)
        if (include_membership) membership[component_levels] <- length(components)
        unseen <- setdiff(unseen, component_persons)
    }
    component_table <- do.call(rbind, components)
    largest_index <- which.max(component_table$Persons)
    summary <- data.frame(Basis = basis, Facet = facet, PersonNodes = length(unique(edges$Person)), FacetLevelNodes = length(unique(edges[[facet]])),
        Edges = nrow(edges), Components = nrow(component_table), LargestComponentPersons = component_table$Persons[largest_index],
        LargestComponentLevels = component_table$FacetLevels[largest_index], LargestComponentPercent = 100 *
            component_table$Persons[largest_index]/length(unique(edges$Person)), Connected = nrow(component_table) ==
            1L, stringsAsFactors = FALSE)
    list(summary = summary, components = component_table, membership = membership)
}

audit_mfrm_connectivity <- function(df, expected_design, facets, include_persons = FALSE) {
    bases <- list(observed = df)
    if (!is.null(expected_design))
        bases$declared_expected <- expected_design
    audits <- unlist(lapply(names(bases), function(basis) {
        lapply(facets, function(facet) {
            mfrm_bipartite_components(bases[[basis]], facet, basis, include_persons = include_persons)
        })
    }), recursive = FALSE)
    list(summary = do.call(rbind, lapply(audits, `[[`, "summary")), components = do.call(rbind, lapply(audits,
        `[[`, "components")))
}

summarize_linkage_by_facet <- function(df, facet, min_linking_persons = 2L) {
    by_level <- dplyr::summarize(dplyr::group_by(df, .data[[facet]]), PersonsPerLevel = dplyr::n_distinct(.data$Person),
        Observations = dplyr::n(), WeightedN = sum(.data$Weight, na.rm = TRUE), .groups = "drop")
    by_person <- dplyr::summarize(dplyr::group_by(df, .data$Person), LevelsPerPerson = dplyr::n_distinct(.data[[facet]]),
        .groups = "drop")
    tibble::tibble(Facet = facet, Levels = nrow(by_level), Persons = nrow(by_person), MinPersonsPerLevel = min(by_level$PersonsPerLevel,
        na.rm = TRUE), MedianPersonsPerLevel = stats::median(by_level$PersonsPerLevel, na.rm = TRUE),
        MinLevelsPerPerson = min(by_person$LevelsPerPerson, na.rm = TRUE), MedianLevelsPerPerson = stats::median(by_person$LevelsPerPerson,
            na.rm = TRUE), LinkingPersons = sum(by_person$LevelsPerPerson > 1L), LinkingPersonRate = mean(by_person$LevelsPerPerson >
            1L), SingleLevelPersons = sum(by_person$LevelsPerPerson == 1L), SparseLevels = sum(by_level$PersonsPerLevel <
            min_linking_persons), SparseLevelThreshold = min_linking_persons)
}

resolve_mfrm_category_policy <- function(keep_original, category_policy, keep_original_supplied, rating_min,
    rating_max) {
    if (!is.logical(keep_original) || length(keep_original) != 1L || is.na(keep_original)) {
        stop("`keep_original` must be a single logical value (TRUE or FALSE).", call. = FALSE)
    }
    if (is.null(category_policy)) {
        if (!keep_original_supplied && (!is.null(rating_min) || !is.null(rating_max)))
            return(TRUE)
        return(keep_original)
    }
    if (!is.character(category_policy) || length(category_policy) != 1L || is.na(category_policy)) {
        stop("`category_policy` must be NULL, 'collapse' or 'preserve'.", call. = FALSE)
    }
    category_policy <- match.arg(category_policy, c("collapse", "preserve"))
    preserve <- identical(category_policy, "preserve")
    if (keep_original_supplied && keep_original != preserve) {
        stop("`category_policy` and `keep_original` conflict; supply one choice or matching values.",
            call. = FALSE)
    }
    preserve
}

describe_mfrm_data <- function(data, person = NULL, facets = NULL, score = NULL, weight = NULL, rating_min = NULL,
    rating_max = NULL, keep_original = FALSE, missing_codes = NULL, include_person_facet = FALSE, include_agreement = TRUE,
    rater_facet = NULL, context_facets = NULL, agreement_top_n = NULL, expected_design = NULL, min_linking_persons = 2L,
    category_policy = NULL, repeated_ratings = c("review", "independent"), comparisons = NULL, spec = NULL) {
    mc <- match.call(expand.dots = FALSE)
    formal_supplied <- setdiff(names(as.list(mc)[-1L]), "spec")
    spec_list <- .mfrm_validate_spec_list(spec, "describe_mfrm_data()")
    desc_spec_keys <- c("person", "facets", "score", "weight", "rating_min", "rating_max", "keep_original",
        "missing_codes", "category_policy", "repeated_ratings")
    for (nm in intersect(names(spec_list), desc_spec_keys)) {
        if (!nm %in% formal_supplied) {
            assign(nm, spec_list[[nm]], envir = environment())
        }
    }
    if ("keep_original" %in% formal_supplied && !"category_policy" %in% formal_supplied) {
        category_policy <- NULL
    }
    keep_original_supplied <- ("keep_original" %in% formal_supplied) || ("keep_original" %in% names(spec_list) &&
        !"category_policy" %in% formal_supplied)
    if (!is.character(person) || length(person) != 1L || is.na(person) || !nzchar(person)) {
        stop("`person` must be a single non-empty character string naming the person column.", call. = FALSE)
    }
    if (!is.character(facets) || length(facets) == 0L || anyNA(facets) || any(!nzchar(facets))) {
        stop("`facets` must be a character vector of one or more facet column names.", call. = FALSE)
    }
    if (!is.character(score) || length(score) != 1L || is.na(score) || !nzchar(score)) {
        stop("`score` must be a single non-empty character string naming the score column.", call. = FALSE)
    }
    repeated_ratings <- match.arg(repeated_ratings, c("review", "independent"))
    keep_original <- resolve_mfrm_category_policy(keep_original, category_policy, keep_original_supplied,
        rating_min, rating_max)
    if (!is.numeric(min_linking_persons) || length(min_linking_persons) != 1L || !is.finite(min_linking_persons) ||
        min_linking_persons < 1L || abs(min_linking_persons - round(min_linking_persons)) > 1e-08) {
        stop("`min_linking_persons` must be a positive integer.", call. = FALSE)
    }
    min_linking_persons <- as.integer(round(min_linking_persons))
    expected_design_normalized <- normalize_expected_mfrm_design(expected_design = expected_design, person = person,
        facets = facets)
    prep <- prepare_mfrm_data(data = data, person_col = person, facet_cols = facets, score_col = score, rating_min = rating_min,
        rating_max = rating_max, weight_col = weight, keep_original = keep_original, missing_codes = missing_codes,
        repeated_ratings = repeated_ratings, allow_single_person = TRUE)
    df <- dplyr::mutate(prep$data, Person = as.character(.data$Person), dplyr::across(dplyr::all_of(prep$facet_names),
        as.character))
    selected_cols <- unique(c(person, facets, score, if (!is.null(weight)) weight))
    missing_by_column <- tibble::tibble(Column = selected_cols, Missing = vapply(selected_cols, function(col) sum(is.na(data[[col]])),
        integer(1)))
    score_desc <- .mfrm_describe_numeric(df$Score)
    weight_desc <- .mfrm_describe_numeric(df$Weight)
    total_weight <- sum(df$Weight, na.rm = TRUE)
    observed_score_distribution <- dplyr::arrange(dplyr::summarize(dplyr::group_by(df, .data$Score), RawN = dplyr::n(),
        WeightedN = sum(.data$Weight, na.rm = TRUE), Percent = ifelse(total_weight > 0, 100 * .data$WeightedN/total_weight,
            NA_real_), .groups = "drop"), .data$Score)
    score_distribution <- dplyr::mutate(dplyr::left_join(tibble::tibble(Score = seq(prep$rating_min, prep$rating_max)),
        observed_score_distribution, by = "Score"), RawN = tidyr::replace_na(.data$RawN, 0L), WeightedN = tidyr::replace_na(.data$WeightedN,
        0), Percent = tidyr::replace_na(.data$Percent, 0))
    score_distribution <- mfrm_original_score_labels(score_distribution, prep$score_map)
    report_facets <- prep$facet_names
    if (isTRUE(include_person_facet)) {
        report_facets <- c("Person", report_facets)
    }
    facet_level_summary <- dplyr::bind_rows(lapply(report_facets, function(facet) {
        dplyr::mutate(dplyr::rename(dplyr::summarize(dplyr::group_by(df, .data[[facet]]), RawN = dplyr::n(),
            WeightedN = sum(.data$Weight, na.rm = TRUE), MeanScore = weighted_mean(.data$Score, .data$Weight),
            SDScore = stats::sd(.data$Score, na.rm = TRUE), MinScore = min(.data$Score, na.rm = TRUE), MaxScore = max(.data$Score,
                na.rm = TRUE), .groups = "drop"), Level = dplyr::all_of(facet)), Facet = facet, Level = as.character(.data$Level),
            .before = 1)
    }))
    linkage_summary <- if (length(prep$facet_names) == 0) {
        tibble::tibble()
    }
    else {
        dplyr::bind_rows(lapply(prep$facet_names, function(facet) {
            summarize_linkage_by_facet(df, facet, min_linking_persons = min_linking_persons)
        }))
    }
    structural_missingness <- audit_mfrm_expected_design(df = df, expected_design = expected_design_normalized,
        facets = prep$facet_names)
    connectivity <- audit_mfrm_connectivity(df = df, expected_design = expected_design_normalized, facets = prep$facet_names,
        include_persons = include_person_facet)
    duplicate_cells <- summarize_duplicate_mfrm_cells(df, prep$facet_names)
    agreement_bundle <- list(summary = data.frame(), pairs = data.frame(), settings = list(included = FALSE,
        status = if (isTRUE(include_agreement)) "not_computed" else "not_requested", rater_facet = NA_character_,
        context_facets = character(0), expected_exact_from_model = FALSE, top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L,
            as.integer(agreement_top_n)), note = if (isTRUE(include_agreement)) {
            "No rater-like facet was detected. Set `rater_facet` explicitly to request facet-agreement summaries."
        } else {
            "Agreement summaries were not requested."
        }))
    if (isTRUE(include_agreement) && length(prep$facet_names) > 0) {
        known_facets <- c("Person", prep$facet_names)
        if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
            rater_facet <- infer_default_rater_facet(prep$facet_names, fallback_first = FALSE)
        }
        else {
            rater_facet <- as.character(rater_facet[1])
        }
        if (is.null(rater_facet)) {
            rater_facet <- NULL
        }
        else if (!rater_facet %in% known_facets) {
            stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
        }
        if (!is.null(rater_facet) && identical(rater_facet, "Person")) {
            stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
        }
        if (!is.null(rater_facet) && is.null(context_facets)) {
            facet_cols <- known_facets
            resolved_context <- setdiff(facet_cols, rater_facet)
        }
        else if (!is.null(rater_facet)) {
            context_facets <- unique(as.character(context_facets))
            unknown <- setdiff(context_facets, known_facets)
            if (length(unknown) > 0) {
                stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "))
            }
            resolved_context <- setdiff(context_facets, rater_facet)
            if (length(resolved_context) == 0) {
                stop("`context_facets` must include at least one facet different from `rater_facet`.")
            }
            facet_cols <- c(rater_facet, resolved_context)
        }
        if (!is.null(rater_facet)) {
            obs_agreement <- dplyr::rename(dplyr::select(df, dplyr::all_of(unique(c("Person", prep$facet_names,
                "Score", "Weight")))), Observed = "Score")
            agreement <- calc_interrater_agreement(obs_df = obs_agreement, facet_cols = facet_cols, rater_facet = rater_facet,
                res = NULL)
            agreement_pairs <- as.data.frame(agreement$pairs, stringsAsFactors = FALSE)
            if (!is.null(agreement_top_n) && nrow(agreement_pairs) > 0) {
                agreement_pairs <- dplyr::slice_head(agreement_pairs, n = max(1L, as.integer(agreement_top_n)))
            }
            agreement_bundle <- list(summary = as.data.frame(agreement$summary, stringsAsFactors = FALSE),
                pairs = agreement_pairs, settings = list(included = TRUE, status = "computed", rater_facet = rater_facet,
                  context_facets = resolved_context, expected_exact_from_model = FALSE, top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L,
                    as.integer(agreement_top_n)), note = "Observed-score agreement was computed for the selected scorer facet."))
        }
    }
    overview <- tibble::tibble(Observations = nrow(df), TotalWeight = total_weight, Persons = length(prep$levels$Person),
        Facets = length(prep$facet_names), Categories = prep$rating_max - prep$rating_min + 1, RatingMin = prep$rating_min,
        RatingMax = prep$rating_max, RatingRangeSource = prep$rating_range_source %||% "unknown", RatingMinSource = prep$rating_min_source %||%
            "unknown", RatingMaxSource = prep$rating_max_source %||% "unknown", OriginalRatingMin = prep$original_rating_min,
        OriginalRatingMax = prep$original_rating_max, !!!mfrm_category_settings(prep))
    facet_crosstabs <- list()
    if (length(prep$facet_names) >= 2L) {
        for (i in seq_len(length(prep$facet_names) - 1L)) {
            for (j in seq(i + 1L, length(prep$facet_names))) {
                a <- prep$facet_names[i]
                b <- prep$facet_names[j]
                ctab <- dplyr::mutate(dplyr::count(df, LevelA = as.character(.data[[a]]), LevelB = as.character(.data[[b]]),
                  name = "N"), FacetA = a, FacetB = b)
                facet_crosstabs[[paste(a, b, sep = "__")]] <- as.data.frame(ctab, stringsAsFactors = FALSE)
            }
        }
    }
    n_total <- nrow(data)
    missing_rate_summary <- data.frame(Column = missing_by_column$Column, Missing = as.integer(missing_by_column$Missing),
        NonMissing = n_total - as.integer(missing_by_column$Missing), MissingRate = ifelse(rep(n_total > 0,
            nrow(missing_by_column)), missing_by_column$Missing/n_total, NA_real_), stringsAsFactors = FALSE)
    out <- list(overview = overview, comparison_review = mfrmr_build_comparison_review(prep, comparisons),
        missing_by_column = missing_by_column, missing_rate_summary = missing_rate_summary, score_descriptives = score_desc,
        weight_descriptives = weight_desc, score_distribution = score_distribution, facet_level_summary = facet_level_summary,
        facet_crosstabs = facet_crosstabs, linkage_summary = linkage_summary, structural_missingness = structural_missingness,
        design_connectivity = as.data.frame(connectivity$summary, stringsAsFactors = FALSE), design_components = as.data.frame(connectivity$components,
            stringsAsFactors = FALSE), repeated_ratings = prep$repeated_ratings, duplicate_cell_summary = duplicate_cells$summary,
        duplicate_cell_detail = duplicate_cells$detail, agreement = agreement_bundle, row_retention = as.data.frame(prep$row_retention %||%
            data.frame(), stringsAsFactors = FALSE), preparation_notes = as.data.frame(prep$preparation_notes %||%
            data.frame(), stringsAsFactors = FALSE), missing_recoding = if (is.data.frame(prep$missing_recoding)) {
            as.data.frame(prep$missing_recoding, stringsAsFactors = FALSE)
        } else {
            prep$missing_recoding
        }, score_support = list(keep_original = prep$keep_original, data = data.frame(Score = sort(unique(df$Score))),
            rating_min = prep$rating_min, rating_max = prep$rating_max, original_rating_min = prep$original_rating_min,
            original_rating_max = prep$original_rating_max, rating_range_source = prep$rating_range_source %||%
                "unknown", rating_min_source = prep$rating_min_source %||% "unknown", rating_max_source = prep$rating_max_source %||%
                "unknown", unused_score_categories = prep$unused_score_categories, score_map = prep$score_map))
    class(out) <- c("mfrm_data_description", class(out))
    out
}

recode_missing_codes <- function(data, columns = NULL, codes = c("99", "999", "-1", "N", "NA", "n/a",
    ".", ""), numeric_codes = TRUE, verbose = FALSE) {
    if (!is.data.frame(data))
        stop("`data` must be a data.frame.", call. = FALSE)
    if (is.null(columns))
        columns <- names(data)
    bad_cols <- setdiff(columns, names(data))
    if (length(bad_cols) > 0L) {
        stop("Unknown column(s): ", paste(bad_cols, collapse = ", "), call. = FALSE)
    }
    codes <- as.character(codes)
    num_codes <- if (isTRUE(numeric_codes)) {
        suppressWarnings(as.numeric(codes))
    }
    else {
        numeric(0)
    }
    num_codes <- num_codes[is.finite(num_codes)]
    counts <- integer(length(columns))
    names(counts) <- columns
    for (col in columns) {
        x <- data[[col]]
        if (is.character(x)) {
            hit <- trimws(x) %in% codes
        }
        else if (is.numeric(x)) {
            hit <- x %in% num_codes
        }
        else {
            hit <- trimws(as.character(x)) %in% codes
        }
        n_hit <- sum(hit, na.rm = TRUE)
        counts[col] <- as.integer(n_hit)
        if (n_hit > 0L) {
            data[[col]][hit] <- NA
        }
    }
    attr(data, "mfrm_missing_recoding") <- data.frame(Column = columns, Replaced = as.integer(counts),
        stringsAsFactors = FALSE)
    if (isTRUE(verbose)) {
        total <- sum(counts, na.rm = TRUE)
        message(sprintf("recode_missing_codes(): replaced %d cell(s) across %d column(s).", total, sum(counts >
            0L, na.rm = TRUE)))
    }
    data
}

print.mfrm_data_description <- function(x, ...) {
    cat("mfrm data description\n")
    mfrmr_print_comparison_review(x$comparison_review)
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
        print(x$overview, row.names = FALSE)
    }
    structural <- as.data.frame(x$structural_missingness$summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(structural) > 0L) {
        cat("\nPlanned assignment coverage\n")
        print(structural, row.names = FALSE)
    }
    print_summary_sections(x, .mfrmr_summary_print_sections[["print.mfrm_data_description_01"]])
    if (!is.null(x$agreement$summary) && nrow(x$agreement$summary) > 0) {
        agreement_facet <- as.character(x$agreement$settings$rater_facet %||% "selected facet")
        cat("\nObserved agreement by ", agreement_facet, "\n", sep = "")
        print(x$agreement$summary, row.names = FALSE)
    }
    invisible(x)
}

collect_mfrm_design_caveats <- function(object) {
    out <- empty_mfrm_caveats
    add_caveat <- function(area, severity, condition, facets, message, action, details = "") {
        out <<- rbind(out, data.frame(Area = area, Severity = severity, Condition = condition, Categories = paste(as.character(facets),
            collapse = ", "), CategoryType = "design", Message = message, RecommendedAction = action,
            Details = details, stringsAsFactors = FALSE))
        invisible(NULL)
    }
    structural <- as.data.frame(object$structural_missingness$summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(structural) > 0L && identical(as.character(structural$Status[1L]), "declared")) {
        missing_n <- suppressWarnings(as.integer(structural$MissingExpectedCells[1L]))
        unexpected_n <- suppressWarnings(as.integer(structural$UnexpectedObservedCells[1L]))
        if (is.finite(missing_n) && missing_n > 0L) {
            add_caveat(area = "structural_missingness", severity = "review", condition = "planned_cells_not_observed",
                facets = object$structural_missingness$settings$key %||% character(0), message = paste0(missing_n,
                  " planned rating cell(s) in `expected_design` were not observed."), action = "Inspect the returned object's `structural_missingness$missing_expected_cells` table and document whether each absence was planned, unavailable, or lost.")
        }
        if (is.finite(unexpected_n) && unexpected_n > 0L) {
            add_caveat(area = "structural_missingness", severity = "warning", condition = "observed_cells_not_in_design",
                facets = object$structural_missingness$settings$key %||% character(0), message = paste0(unexpected_n,
                  " observed rating cell(s) were not present in `expected_design`."), action = "Reconcile the assignment roster and observed data before fitting.")
        }
    }
    connectivity <- as.data.frame(object$design_connectivity %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(connectivity) > 0L && all(c("Connected", "Basis", "Facet", "Components") %in% names(connectivity))) {
        disconnected <- connectivity[!as.logical(connectivity$Connected), , drop = FALSE]
        for (i in seq_len(nrow(disconnected))) {
            basis <- as.character(disconnected$Basis[i])
            facet <- as.character(disconnected$Facet[i])
            components <- suppressWarnings(as.integer(disconnected$Components[i]))
            add_caveat(area = "design_connectivity", severity = "warning", condition = paste0("disconnected_",
                basis, "_person_facet_graph"), facets = facet, message = paste0("The ", gsub("_", " ",
                basis), " Person-", facet, " graph has ", components, " disconnected components."), action = paste0("Inspect the returned object's `design_components` table for ",
                facet, " and add defensible linking observations or analyze components separately."))
        }
    }
    duplicates <- as.data.frame(object$duplicate_cell_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(duplicates) > 0L && isTRUE(duplicates$HasDuplicates[1L])) {
        add_caveat(area = "duplicate_cells", severity = if (identical(object$repeated_ratings, "independent"))
            "info"
        else "warning", condition = if (identical(object$repeated_ratings, "independent"))
            "independent_repeated_ratings_declared"
        else "duplicate_person_facet_cells", facets = object$structural_missingness$settings$key %||%
            character(0), message = paste0(duplicates$DuplicateCells[1L], " Person x facet cell(s) contain multiple observed rows."),
            action = if (identical(object$repeated_ratings, "independent"))
                "Conditional independence was declared by the analyst; the data check does not verify it."
            else "Review repeated events: correct accidental copies, model dependence, or explicitly declare conditional independence.")
    }
    linkage <- as.data.frame(object$linkage_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(linkage) > 0L && all(c("Facet", "SparseLevels", "SparseLevelThreshold") %in% names(linkage))) {
        sparse <- linkage[suppressWarnings(as.integer(linkage$SparseLevels)) > 0L, , drop = FALSE]
        for (i in seq_len(nrow(sparse))) {
            add_caveat(area = "design_linkage", severity = "review", condition = "sparse_facet_levels",
                facets = sparse$Facet[i], message = paste0(sparse$SparseLevels[i], " level(s) of ", sparse$Facet[i],
                  " were observed for fewer than ", sparse$SparseLevelThreshold[i], " distinct persons."),
                action = "Review the affected level support and intended inference; the threshold is descriptive, not an automatic exclusion rule.")
        }
    }
    out
}

summary.mfrm_data_description <- function(object, digits = 3, top_n = 10, ...) {
    digits <- max(0L, as.integer(digits))
    top_n <- max(1L, as.integer(top_n))
    overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(overview) > 0L) {
        category_settings <- mfrm_category_settings(object$score_support %||% object$prep)
        for (field in setdiff(names(category_settings), names(overview))) {
            overview[[field]] <- category_settings[[field]]
        }
    }
    missing_tbl <- as.data.frame(object$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(missing_tbl) > 0 && all(c("Column", "Missing") %in% names(missing_tbl))) {
        missing_tbl <- dplyr::slice_head(dplyr::arrange(missing_tbl, dplyr::desc(.data$Missing), .data$Column),
            n = top_n)
    }
    score_dist_full <- as.data.frame(object$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
    score_dist_full <- mfrm_original_score_labels(score_dist_full, object$score_support$score_map)
    score_dist <- score_dist_full
    if (nrow(score_dist) > 0) {
        score_dist <- utils::head(score_dist, n = top_n)
    }
    facet_tbl <- as.data.frame(object$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
    facet_overview <- data.frame()
    if (nrow(facet_tbl) > 0 && all(c("Facet", "Level", "WeightedN") %in% names(facet_tbl))) {
        facet_overview <- as.data.frame(dplyr::arrange(dplyr::summarise(dplyr::group_by(facet_tbl, .data$Facet),
            Levels = dplyr::n_distinct(.data$Level), TotalWeightedN = sum(.data$WeightedN, na.rm = TRUE),
            MeanWeightedN = mean(.data$WeightedN, na.rm = TRUE), MinWeightedN = min(.data$WeightedN,
                na.rm = TRUE), MaxWeightedN = max(.data$WeightedN, na.rm = TRUE), .groups = "drop"),
            .data$Facet), stringsAsFactors = FALSE)
    }
    agreement_tbl <- as.data.frame(object$agreement$summary %||% data.frame(), stringsAsFactors = FALSE)
    row_retention <- as.data.frame(object$row_retention %||% data.frame(), stringsAsFactors = FALSE)
    preparation_notes <- as.data.frame(object$preparation_notes %||% data.frame(), stringsAsFactors = FALSE)
    structural_missingness <- as.data.frame(object$structural_missingness$summary %||% data.frame(),
        stringsAsFactors = FALSE)
    structural_level_coverage <- as.data.frame(object$structural_missingness$level_coverage %||% data.frame(),
        stringsAsFactors = FALSE)
    design_connectivity <- as.data.frame(object$design_connectivity %||% data.frame(), stringsAsFactors = FALSE)
    design_components <- as.data.frame(object$design_components %||% data.frame(), stringsAsFactors = FALSE)
    linkage_summary <- as.data.frame(object$linkage_summary %||% data.frame(), stringsAsFactors = FALSE)
    duplicate_cell_summary <- as.data.frame(object$duplicate_cell_summary %||% data.frame(), stringsAsFactors = FALSE)
    reporting_map <- data.frame(Area = c("Sample / design counts", "Missingness review", "Planned assignment coverage",
        "Person-facet connectivity", "Score usage / category distribution", "Facet coverage", "Rater-facet agreement",
        "Fit / reliability / residual PCA"), CoveredHere = c("yes", "yes", if (nrow(structural_missingness) >
        0L && identical(as.character(structural_missingness$Status[1L]), "declared")) "yes" else "not declared",
        "yes", "yes", "yes", if (nrow(agreement_tbl) > 0) "yes" else "no", "no"), CompanionOutput = c("summary(describe_mfrm_data(...))",
        "summary(describe_mfrm_data(...))", "data_review$structural_missingness$missing_expected_cells",
        "data_review$design_components", "summary(describe_mfrm_data(...))", "summary(describe_mfrm_data(...))",
        "summary(describe_mfrm_data(...)) / plot(x, type = \"agreement\")", "summary(diagnose_mfrm(fit))"),
        stringsAsFactors = FALSE)
    notes <- character(0)
    if (nrow(missing_tbl) > 0 && any(suppressWarnings(as.numeric(missing_tbl$Missing)) > 0, na.rm = TRUE)) {
        notes <- c(notes, "Missing values were detected in one or more input columns.")
    }
    else {
        notes <- c(notes, "No missing values were detected in selected input columns.")
    }
    caveat_prep <- object$score_support %||% object$prep %||% NULL
    caveats <- collect_mfrm_caveats(prep = caveat_prep, score_distribution = score_dist_full, include_recode = TRUE,
        context = "data")
    if (nrow(caveats) > 0 && "Message" %in% names(caveats)) {
        notes <- c(notes, as.character(caveats$Message))
    }
    design_caveats <- collect_mfrm_design_caveats(object)
    caveats <- rbind(caveats, design_caveats)
    if (nrow(design_caveats) > 0L && "Message" %in% names(design_caveats)) {
        notes <- c(notes, as.character(design_caveats$Message))
    }
    structural_note <- as.character(object$structural_missingness$settings$note %||% "")
    if (nzchar(structural_note))
        notes <- c(notes, structural_note)
    prep_review_messages <- if (nrow(preparation_notes) > 0L && all(c("Severity", "Message") %in% names(preparation_notes))) {
        preparation_notes$Message[tolower(as.character(preparation_notes$Severity)) %in% c("review",
            "warning", "error")]
    }
    else {
        character(0)
    }
    notes <- c(notes, prep_review_messages)
    out <- list(overview = overview, comparison_review = object$comparison_review, missing = missing_tbl,
        score_distribution = score_dist, facet_overview = facet_overview, structural_missingness = structural_missingness,
        structural_level_coverage = structural_level_coverage, design_connectivity = design_connectivity,
        design_components = design_components, linkage_summary = linkage_summary, duplicate_cell_summary = duplicate_cell_summary,
        agreement = agreement_tbl, agreement_settings = object$agreement$settings %||% list(), row_retention = row_retention,
        preparation_notes = preparation_notes, reporting_map = reporting_map, caveats = caveats, notes = notes,
        digits = digits)
    class(out) <- "summary.mfrm_data_description"
    out
}

print.summary.mfrm_data_description <- function(x, ...) {
    digits <- as.integer(x$digits %||% 3L)
    if (!is.finite(digits))
        digits <- 3L
    cat("mfrm Data Description Summary\n")
    mfrmr_print_comparison_review(x$comparison_review)
    print_summary_sections(x, .mfrmr_summary_print_sections[["print.summary.mfrm_data_description_01"]],
        digits = digits)
    if (!is.null(x$row_retention) && nrow(x$row_retention) > 0 && "DroppedRows" %in% names(x$row_retention) &&
        any(suppressWarnings(as.numeric(x$row_retention$DroppedRows)) > 0, na.rm = TRUE)) {
        cat("\nRow retention\n")
        print(round_numeric_df(as.data.frame(x$row_retention), digits = digits), row.names = FALSE)
    }
    print_note_section(notes = x$preparation_notes, title = "Data preparation notes", preparation = TRUE)
    if (!is.null(x$agreement) && nrow(x$agreement) > 0) {
        agreement_facet <- as.character(x$agreement_settings$rater_facet %||% "selected facet")
        cat("\nObserved agreement by ", agreement_facet, "\n", sep = "")
        print(round_numeric_df(as.data.frame(x$agreement), digits = digits), row.names = FALSE)
    }
    print_note_section(notes = x$caveats, preparation = FALSE)
    print_summary_sections(x, .mfrmr_summary_print_sections[["print.summary.mfrm_data_description_02"]])
    note_lines <- as.character(x$notes %||% character(0))
    note_lines <- note_lines[nzchar(note_lines)]
    if (length(note_lines) > 0L) {
        cat("\nNotes\n")
        cat(paste0(" - ", note_lines, collapse = "\n"), "\n", sep = "")
    }
    invisible(x)
}

print.mfrm_anchor_review <- function(x, ...) {
    issue_total <- if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0)
        sum(x$issue_counts$N)
    else 0
    cat("mfrm anchor review\n")
    cat("  issue rows: ", issue_total, "\n", sep = "")
    if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
        nonzero <- dplyr::filter(x$issue_counts, .data$N > 0)
        if (nrow(nonzero) > 0) {
            cat("\nIssue counts\n")
            print(nonzero, row.names = FALSE)
        }
    }
    print_summary_sections(x, .mfrmr_summary_print_sections[["print.mfrm_anchor_review_01"]])
    if (!is.null(x$design_checks) && !is.null(x$design_checks$level_observation_summary) && nrow(x$design_checks$level_observation_summary) >
        0) {
        cat("\nLevel observation summary\n")
        print(x$design_checks$level_observation_summary, row.names = FALSE)
    }
    invisible(x)
}

summary.mfrm_anchor_review <- function(object, digits = 3, top_n = 10, ...) {
    digits <- max(0L, as.integer(digits))
    top_n <- max(1L, as.integer(top_n))
    issue_counts <- as.data.frame(object$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(issue_counts) > 0 && all(c("Issue", "N") %in% names(issue_counts))) {
        issue_counts <- dplyr::slice_head(dplyr::arrange(dplyr::filter(issue_counts, .data$N > 0), dplyr::desc(.data$N),
            .data$Issue), n = top_n)
    }
    facet_summary <- as.data.frame(object$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
    level_summary <- as.data.frame(object$design_checks$level_observation_summary %||% data.frame(),
        stringsAsFactors = FALSE)
    category_summary <- as.data.frame(object$design_checks$category_counts %||% data.frame(), stringsAsFactors = FALSE)
    recommendations <- as.character(object$recommendations %||% character(0))
    if (length(recommendations) > top_n) {
        recommendations <- recommendations[seq_len(top_n)]
    }
    notes <- if (nrow(issue_counts) > 0) {
        "Anchor-review issues were detected. Review issue counts and recommendations."
    }
    else {
        "No anchor-table issue rows were detected."
    }
    out <- list(issue_counts = issue_counts, facet_summary = facet_summary, level_observation_summary = level_summary,
        category_counts = category_summary, recommendations = recommendations, notes = notes, digits = digits)
    class(out) <- "summary.mfrm_anchor_review"
    out
}

make_anchor_table <- function(fit, facets = NULL, include_person = FALSE, digits = NULL, readiness_policy = c("error",
    "review")) {
    .mfrm_validate_fit(fit)
    readiness_policy <- match.arg(readiness_policy)
    if (!isTRUE(fit$converged) && identical(readiness_policy, "error")) {
        stop("`fit` did not converge; resolve convergence or set `readiness_policy = 'review'`.", call. = FALSE)
    }
    out <- fit$facets[, c("Facet", "Level", "Estimate"), drop = FALSE]
    names(out)[3L] <- "Anchor"
    if (isTRUE(include_person)) {
        p_anc <- data.frame(Facet = "Person", Level = fit$person$Person, Anchor = fit$person$Estimate, stringsAsFactors = FALSE)
        out <- rbind(p_anc, out)
    }
    if (!is.null(facets))
        out <- out[out$Facet %in% as.character(facets), , drop = FALSE]
    if (!is.null(digits))
        out$Anchor <- round(out$Anchor, digits = as.integer(digits))
    rownames(out) <- NULL
    return(tibble::as_tibble(out))
}

diagnose_mfrm <- function(fit, residual_pca = c("none", "overall", "facet", "both"), pca_max_factors = 10L) {
    .mfrm_validate_fit(fit)
    return(.mfrm_diagnose_fit(fit, residual_pca = match.arg(residual_pca), pca_max_factors = pca_max_factors))
}

compare_mfrm <- function(..., labels = NULL, warn_constraints = TRUE, nested = FALSE) {
    if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
        stop("`nested` must be a single logical value (TRUE or FALSE).", call. = FALSE)
    }
    fits <- list(...)
    plan_input <- vapply(fits, inherits, logical(1), what = "mfrm_multivariate_d_study")
    if (any(plan_input)) {
        if (!plan_input[1L] || sum(plan_input) != 1L)
            stop("Supply one D-study containing all plans to compare.", call. = FALSE)
        if (!missing(labels) || !missing(warn_constraints) || !missing(nested))
            stop("Model-comparison arguments do not apply to a D-study plan comparison.", call. = FALSE)
        x <- fits[[1L]]
        defaults <- list(reference = 1L, score = NULL, composite = NULL, assumption = NULL, level = 0.95)
        options <- fits[-1L]
        if (length(options) && (is.null(names(options)) || anyNA(names(options)) || any(!nzchar(names(options))) ||
            anyDuplicated(names(options)) || any(!names(options) %in% names(defaults))))
            stop("Name D-study options with `reference`, `score`, `composite`, `assumption` or `level`.",
                call. = FALSE)
        defaults[names(options)] <- options
        reference <- defaults$reference
        score <- defaults$score
        composite <- defaults$composite
        assumption <- defaults$assumption
        level <- defaults$level
        if (!identical(assumption, "normal")) {
            stop("This approximation requires `assumption = 'normal'` for the random effects; it does not test that assumption.",
                call. = FALSE)
        }
        .mfrm_mvds_require_raw(x)
        d <- mfrm_d_study(x$gstudy, x$design_grid, x$weights)
        g <- d$gstudy
        if (length(g$design$counts) != 3L || !is.null(g$design$nesting)) {
            stop("This interval method currently supports two common crossed random facets; nested-design intervals are not available.",
                call. = FALSE)
        }
        nplans <- nrow(d$design_grid)
        if (nplans < 2L || !is.numeric(reference) || is.complex(reference) || is.object(reference) || !is.null(dim(reference)) ||
            length(reference) != 1L || !is.finite(reference) || reference != floor(reference) || reference <
            1L || reference > nplans) {
            stop("Supply at least two plans and a valid reference row number.", call. = FALSE)
        }
        if (!is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) || length(level) !=
            1L || !is.finite(level) || level <= 0 || level >= 1) {
            stop("`level` must be a number strictly between zero and one.", call. = FALSE)
        }
        target <- .mfrm_mvds_interval_target(d, score, composite, comparison = TRUE)
        kind <- target$kind
        selected <- target$selected
        w <- target$weights
        tab <- target$table
        weight_scale <- max(abs(w))
        unit_w <- w/weight_scale
        theta <- vapply(g$components, function(a) drop(crossprod(unit_w, a %*% unit_w)), numeric(1))
        component_scale <- max(abs(theta))
        if (component_scale == 0)
            component_scale <- 1
        theta <- theta/component_scale
        sem_scale <- sqrt(component_scale) * weight_scale
        C <- .mfrm_mvgt_component_covariance(g, theta)
        eigenvalues <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
        covariance_psd <- min(eigenvalues) >= -sqrt(.Machine$double.eps) * max(abs(eigenvalues))
        metrics <- c("G", "Phi", "RelativeSEM", "AbsoluteSEM")
        gradients <- lapply(seq_len(nplans), function(i) {
            counts <- as.numeric(d$design_grid[i, ])
            r <- unname(.mfrm_mvgt_profile_error(d, i, "G", counts))
            a <- unname(.mfrm_mvgt_profile_error(d, i, "Phi", counts))
            J <- matrix(NA_real_, 4, 7)
            for (j in 1:2) {
                c <- if (j == 1L)
                  r
                else a
                e <- sum(c * theta)
                u <- theta[1L]
                if (u > 0 && e > 0) {
                  J[j, ] <- -u * c/(u + e)^2
                  J[j, 1L] <- e/(u + e)^2
                }
                if (e > 0)
                  J[j + 2L, ] <- sem_scale * c/(2 * sqrt(e))
            }
            J
        })
        others <- setdiff(seq_len(nplans), reference)
        J <- do.call(rbind, lapply(others, function(i) gradients[[i]] - gradients[[reference]]))
        table <- do.call(rbind, lapply(others, function(i) data.frame(Reference = reference, Scenario = i,
            Kind = kind, Score = selected, Metric = metrics, ReferenceValue = as.numeric(tab[reference, metrics]),
            Value = as.numeric(tab[i, metrics]), Difference = as.numeric(tab[i, metrics]) - as.numeric(tab[reference,
                metrics]), Status = ifelse(is.finite(as.numeric(tab[i, metrics])) & is.finite(as.numeric(tab[reference,
                metrics])), "Available", "Point projection unavailable"), row.names = NULL)))
        identical_plan <- vapply(others, function(i) all(d$design_grid[i, ] == d$design_grid[reference, ]),
            logical(1))
        exact <- rep(identical_plan, each = 4L) & table$Status == "Available"
        J[exact, ] <- 0
        covariance <- J %*% C %*% t(J)
        rownames(covariance) <- colnames(covariance) <- paste(table$Scenario, table$Metric, sep = "/")
        variance <- diag(covariance)
        boundary <- rowSums(!is.finite(J)) > 0
        table$Status[table$Status == "Available" & boundary] <- "Interval unavailable at a boundary"
        table$Status[table$Status == "Available" & (!is.finite(variance) | variance <= 0) & !exact] <- "Nonpositive or nonfinite difference variance"
        ok <- table$Status == "Available"
        table$SE <- sqrt(ifelse(ok, variance, NA_real_))
        critical <- stats::qnorm((1 - level)/2, lower.tail = FALSE)
        table$Lower <- table$Difference - critical * table$SE
        table$Upper <- table$Difference + critical * table$SE
        return(structure(list(comparisons = table, covariance = covariance, design_grid = d$design_grid, reference = as.integer(reference),
            score = selected, kind = kind, weights = w, level = level, method = "Normal-theory paired delta approximation",
            source_design = g$design[c("method", "counts", "rows", "complete", "balanced", "nesting", "score_convention")],
            source_rows = g$data_usage$counts, component_diagnostics = d$component_diagnostics, sampling_covariance_psd = covariance_psd),
            class = "mfrm_multivariate_d_comparison"))
    }
    .mfrm_compare_fits(fits, labels = labels, nested = nested, warn_constraints = warn_constraints)
}

summary.mfrm_comparison <- function(object, ...) {
    out <- list(table = object$table, lrt = object$lrt, evidence_ratios = object$evidence_ratios, preferred = object$preferred,
        comparison_basis = object$comparison_basis)
    class(out) <- "summary.mfrm_comparison"
    out
}


print.mfrm_comparison <- function(x, ...) {
    print(summary(x))
    invisible(x)
}
