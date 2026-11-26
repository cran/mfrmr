build_mfrm_step_map <- function(declarations, config, prep) {
    required <- c("Level", "From", "To")
    optional <- c("Value", "Group")
    if (!is.data.frame(declarations) || !nrow(declarations) || anyDuplicated(names(declarations)) ||
        !all(required %in% names(declarations)) || !any(optional %in% names(declarations)) || !all(names(declarations) %in%
        c(required, optional))) {
        stop("`step_constraints` needs a data frame with Level, From, To, and Value or Group.", call. = FALSE)
    }
    spec <- declarations
    owners <- names(config$step_specs)
    categories <- prep$score_map
    if (is.null(categories) || !identical(as.numeric(categories$OriginalScore), as.numeric(categories$InternalScore))) {
        stop("`step_constraints` requires consecutive, unrelabelled score categories.", call. = FALSE)
    }
    from <- seq.int(prep$rating_min, prep$rating_max - 1L)
    n_steps <- length(from)
    owner <- match(spec$Level, owners)
    transition <- match(spec$From, from)
    if (!is.character(spec$Level) || !is.null(dim(spec$Level)) || anyNA(owner) || !is.null(dim(spec$From)) ||
        !is.null(dim(spec$To)) || !is.numeric(spec$From) || is.complex(spec$From) || !is.numeric(spec$To) ||
        is.complex(spec$To) || anyNA(transition) || anyNA(spec$To) || any(spec$To != from[transition] +
        1L)) {
        stop("`step_constraints` needs a known Level and consecutive From/To scores in the rating range.",
            call. = FALSE)
    }
    row <- (owner - 1L) * n_steps + transition
    if (anyDuplicated(row))
        stop("Declare each Level and transition only once in `step_constraints`.", call. = FALSE)
    value <- spec$Value %||% rep(NA_real_, nrow(spec))
    group <- spec$Group %||% rep(NA_character_, nrow(spec))
    if (is.logical(value) && all(is.na(value)))
        value <- as.numeric(value)
    if (is.logical(group) && all(is.na(group)))
        group <- as.character(group)
    if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) || any(is.nan(value) | is.infinite(value)) ||
        !is.character(group) || !is.null(dim(group)) || any(!is.na(group) & !nzchar(group)) || any(is.na(value) &
        is.na(group))) {
        stop("Each step declaration needs a finite Value or a nonempty Group (use NA for omissions).",
            call. = FALSE)
    }
    n_rows <- length(owners) * n_steps
    representatives <- seq_len(n_rows)
    groups <- rep(NA_character_, n_rows)
    groups[row] <- group
    for (member in split(seq_len(nrow(spec)), group)) {
        if (length(member) < 2L || length(unique(transition[member])) != 1L) {
            stop("Each step Group must share the same transition across at least two levels.", call. = FALSE)
        }
        if (length(unique(value[member][!is.na(value[member])])) > 1L) {
            stop("Conflicting fixed Values in a shared step Group.", call. = FALSE)
        }
        representatives[row[member]] <- min(row[member])
    }
    cells <- sort(unique(representatives))
    cell <- match(representatives, cells)
    fixed <- rep(NA_real_, length(cells))
    fixed[cell[row[!is.na(value)]]] <- value[!is.na(value)]
    offset <- replace(fixed, is.na(fixed), 0)
    available <- which(is.na(fixed))
    centering <- matrix(0, length(owners), length(cells))
    centering[cbind(rep(seq_along(owners), each = n_steps), cell)] <- 1
    system <- centering[, rev(available), drop = FALSE]
    q <- qr(system)
    dependent <- rev(available)[q$pivot[seq_len(q$rank)]]
    free <- setdiff(available, dependent)
    basis <- matrix(0, length(cells), length(free))
    basis[cbind(free, seq_along(free))] <- 1
    if (length(dependent)) {
        calculated <- qr.coef(qr(centering[, dependent, drop = FALSE]), -cbind(centering %*% offset,
            centering[, free, drop = FALSE]))
        offset[dependent] <- calculated[, 1L]
        basis[dependent, ] <- calculated[, -1L, drop = FALSE]
    }
    basis[abs(basis) < 1e-12] <- 0
    lhs <- as.numeric(centering %*% offset)
    scale <- as.numeric(centering %*% abs(offset))
    if (any(!is.finite(lhs)) || any(!is.finite(scale)) || any(abs(lhs) > 1e-12 * pmax(1, scale)) || any(abs(centering %*%
        basis) > 1e-12)) {
        stop("Fixed or shared steps conflict with within-ladder sum-to-zero constraints.", call. = FALSE)
    }
    list(offset = offset[cell], basis = Matrix::Matrix(basis[cell, , drop = FALSE], sparse = TRUE), free = cells[free],
        groups = groups, fixed = !is.na(fixed[cell]))
}


mfrm_affine_subspace_contained <- function(sub, sup, tol = 1e-10) {
    if (length(sub$offset) != length(sup$offset))
        return(FALSE)
    delta_c <- sub$offset - sup$offset
    res_c <- if (sup$rank == 0L) {
        delta_c
    }
    else {
        delta_c - drop(sup$Q %*% crossprod(sup$Q, delta_c))
    }
    scale_c <- max(1, max(abs(sub$offset)), max(abs(sup$offset)))
    if (!all(is.finite(res_c)) || max(abs(res_c)) > tol * scale_c)
        return(FALSE)
    if (sub$rank == 0L)
        return(TRUE)
    if (sup$rank < sub$rank)
        return(FALSE)
    res_B <- sub$Q - sup$Q %*% crossprod(sup$Q, sub$Q)
    all(is.finite(res_B)) && max(abs(res_B)) <= tol
}
