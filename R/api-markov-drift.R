.drift_quadrature_normal <- function(n_quad = 81L, mean = 0, sd = 1, bound = 6) {
    z <- seq(-bound, bound, length.out = as.integer(n_quad))
    weights <- stats::dnorm(z)
    weights <- weights/sum(weights)
    list(nodes = mean + sd * z, weights = weights, bound = bound)
}

.drift_rowsum_by_id <- function(x, id, n_id) {
    rs <- rowsum(x, id, reorder = FALSE)
    idx <- as.integer(rownames(rs))
    if (is.matrix(x)) {
        out <- matrix(0, nrow = n_id, ncol = ncol(x))
        out[idx, ] <- rs
    }
    else {
        out <- numeric(n_id)
        out[idx] <- as.numeric(rs)
    }
    out
}

.drift_order_time_levels <- function(x) {
    if (is.factor(x)) {
        lv <- levels(droplevels(x))
        return(lv[nzchar(trimws(lv))])
    }
    if (is.numeric(x) || inherits(x, c("Date", "POSIXt"))) {
        u <- sort(unique(x))
        return(as.character(u))
    }
    u <- unique(trimws(as.character(x)))
    u <- u[!is.na(u) & nzchar(u)]
    u_num <- suppressWarnings(as.numeric(u))
    if (length(u) > 0L && !anyNA(u_num) && !anyDuplicated(u_num)) {
        return(u[order(u_num)])
    }
    sort(u)
}

.drift_order_id_levels <- function(x) {
    if (is.factor(x)) {
        lv <- levels(droplevels(x))
        return(lv[nzchar(trimws(lv))])
    }
    u <- unique(trimws(as.character(x)))
    u[!is.na(u) & nzchar(u)]
}

.drift_prepare_input <- function(data, person, rater, time, score, facets = character(), score_levels = NULL,
    missing = c("fail", "omit")) {
    missing <- match.arg(missing)
    if (!is.data.frame(data) || nrow(data) == 0L) {
        stop("`data` must be a non-empty data.frame.", call. = FALSE)
    }
    for (nm in list(person = person, rater = rater, time = time, score = score)) {
        if (!is.character(nm) || length(nm) != 1L || is.na(nm) || !nzchar(nm)) {
            stop("`person`, `rater`, `time`, and `score` must each be a single column name.", call. = FALSE)
        }
    }
    if (is.null(facets))
        facets <- character()
    if (!is.character(facets) || anyNA(facets) || any(!nzchar(facets))) {
        stop("`facets` must be a character vector of column names.", call. = FALSE)
    }
    facets <- unique(facets)
    role_cols <- c(person, rater, time, score, facets)
    if (anyDuplicated(c(person, rater, time, score)) || any(facets %in% c(person, rater, time, score))) {
        stop("`person`, `rater`, `time`, `score`, and `facets` must name distinct columns.", call. = FALSE)
    }
    missing_cols <- setdiff(role_cols, names(data))
    if (length(missing_cols) > 0L) {
        stop("Missing required column(s) in `data`: ", paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }
    df <- as.data.frame(data[, role_cols, drop = FALSE], stringsAsFactors = FALSE)
    n_raw <- nrow(df)
    has_na <- Reduce(`|`, lapply(df, function(col) {
        if (is.character(col) || is.factor(col)) {
            is.na(col) | !nzchar(trimws(as.character(col)))
        }
        else {
            is.na(col)
        }
    }))
    n_omitted <- sum(has_na)
    if (n_omitted > 0L) {
        if (identical(missing, "fail")) {
            stop(sprintf("`data` contains %d incomplete row(s); set `missing = 'omit'` to drop them explicitly.",
                n_omitted), call. = FALSE)
        }
        df <- df[!has_na, , drop = FALSE]
    }
    if (nrow(df) < 4L) {
        stop("At least 4 complete rating rows are required to fit a drift MFRM.", call. = FALSE)
    }
    raw_score <- df[[score]]
    if (!is.numeric(raw_score) || is.complex(raw_score) || any(!is.finite(raw_score)) || any(raw_score !=
        floor(raw_score))) {
        stop("`score` must contain finite integer category values.", call. = FALSE)
    }
    raw_score <- as.integer(raw_score)
    if (is.null(score_levels)) {
        s_min <- min(raw_score)
        s_max <- max(raw_score)
        if (s_max <= s_min) {
            stop("`score` must contain at least two distinct rating categories.", call. = FALSE)
        }
        score_levels <- seq.int(s_min, s_max)
    }
    else {
        if (!is.numeric(score_levels) || is.complex(score_levels) || length(score_levels) < 2L || any(!is.finite(score_levels)) ||
            any(score_levels != floor(score_levels))) {
            stop("`score_levels` must be a vector of at least 2 finite integers.", call. = FALSE)
        }
        score_levels <- sort(unique(as.integer(score_levels)))
        if (length(score_levels) < 2L || !all(diff(score_levels) == 1L)) {
            stop("`score_levels` must be consecutive integers (e.g., `0:4` or `1:5`).", call. = FALSE)
        }
        if (any(!raw_score %in% score_levels)) {
            stop("Observed `score` values fall outside `score_levels`.", call. = FALSE)
        }
    }
    K <- length(score_levels)
    x_idx <- match(raw_score, score_levels)
    person_levels <- .drift_order_id_levels(df[[person]])
    rater_levels <- .drift_order_id_levels(df[[rater]])
    time_levels <- .drift_order_time_levels(df[[time]])
    if (length(person_levels) < 2L)
        stop("At least 2 distinct persons are required.", call. = FALSE)
    if (length(rater_levels) < 2L)
        stop("At least 2 distinct raters are required.", call. = FALSE)
    if (length(time_levels) < 2L)
        stop("At least 2 distinct time points/waves are required.", call. = FALSE)
    j_idx <- match(trimws(as.character(df[[person]])), person_levels)
    r_idx <- match(trimws(as.character(df[[rater]])), rater_levels)
    t_idx <- match(if (is.factor(df[[time]]) || is.numeric(df[[time]]) || inherits(df[[time]], c("Date",
        "POSIXt"))) {
        as.character(df[[time]])
    }
    else {
        trimws(as.character(df[[time]]))
    }, time_levels)
    facet_levels <- list()
    facet_idx <- list()
    facet_basis <- list()
    active_facets <- character()
    fixed_separated <- FALSE
    for (fc in facets) {
        lv <- .drift_order_id_levels(df[[fc]])
        f_id <- match(trimws(as.character(df[[fc]])), lv)
        facet_levels[[fc]] <- lv
        facet_idx[[fc]] <- f_id
        if (length(lv) >= 2L) {
            active_facets <- c(active_facets, fc)
            B <- stats::contr.sum(length(lv))
            rownames(B) <- lv
            colnames(B) <- lv[-length(lv)]
            facet_basis[[fc]] <- B
            by_lv <- split(x_idx, factor(f_id, levels = seq_along(lv)))
            if (any(vapply(by_lv, function(v) length(v) == 0L || all(v == 1L) || all(v == K), logical(1)))) {
                fixed_separated <- TRUE
            }
        }
    }
    pw_mat <- table(j_idx, t_idx) > 0L
    person_waves <- rowSums(pw_mat)
    bridge_persons <- sum(person_waves >= 2L)
    Tn <- length(time_levels)
    adj_w <- crossprod(pw_mat) > 0L
    visited <- rep(FALSE, Tn)
    n_wave_comp <- 0L
    for (w0 in seq_len(Tn)) {
        if (!visited[w0]) {
            n_wave_comp <- n_wave_comp + 1L
            queue <- w0
            visited[w0] <- TRUE
            while (length(queue) > 0L) {
                cur <- queue[1L]
                queue <- queue[-1L]
                nbrs <- which(adj_w[cur, ] & !visited)
                if (length(nbrs) > 0L) {
                  visited[nbrs] <- TRUE
                  queue <- c(queue, nbrs)
                }
            }
        }
    }
    rc_counts <- matrix(0L, nrow = length(rater_levels), ncol = K, dimnames = list(rater_levels, as.character(score_levels)))
    rc_tab <- table(factor(r_idx, levels = seq_along(rater_levels)), factor(x_idx, levels = seq_len(K)))
    rc_counts[, ] <- as.integer(rc_tab)
    rt_counts <- matrix(0L, nrow = length(rater_levels), ncol = Tn, dimnames = list(rater_levels, time_levels))
    rt_tab <- table(factor(r_idx, levels = seq_along(rater_levels)), factor(t_idx, levels = seq_len(Tn)))
    rt_counts[, ] <- as.integer(rt_tab)
    list(data = df, j = as.integer(j_idx), r = as.integer(r_idx), t = as.integer(t_idx), x = as.integer(x_idx),
        raw_score = raw_score, N = nrow(df), N_raw = n_raw, N_omitted = n_omitted, J = length(person_levels),
        R = length(rater_levels), T = Tn, K = K, score_levels = score_levels, columns = list(person = person,
            rater = rater, time = time, score = score, facets = facets), levels = c(list(person = person_levels,
            rater = rater_levels, time = time_levels), facet_levels), active_facets = active_facets,
        facet_idx = facet_idx, facet_basis = facet_basis, fixed_separated = fixed_separated, bridge_persons = as.integer(bridge_persons),
        person_waves = as.integer(person_waves), wave_components = as.integer(n_wave_comp), rc_counts = rc_counts,
        rt_counts = rt_counts)
}

.drift_layout <- function(dat, time_mode = c("dynamic", "static"), step_mode = c("rater", "shared")) {
    time_mode <- match.arg(time_mode)
    step_mode <- match.arg(step_mode)
    R <- dat$R
    Tn <- dat$T
    K <- dat$K
    n_u <- K - 2L
    n_beta <- if (time_mode == "dynamic")
        R * Tn
    else R
    n_step <- if (step_mode == "rater")
        R * n_u
    else n_u
    f_names <- dat$active_facets
    f_sizes_tilde <- integer(length(f_names))
    f_sizes_free <- integer(length(f_names))
    for (i in seq_along(f_names)) {
        B <- dat$facet_basis[[f_names[i]]]
        f_sizes_tilde[i] <- nrow(B)
        f_sizes_free[i] <- ncol(B)
    }
    n_facet_tilde <- sum(f_sizes_tilde)
    n_facet_free <- sum(f_sizes_free)
    n_tilde <- n_beta + n_facet_tilde + n_step
    n_free <- n_beta + n_facet_free + n_step
    has_facets <- length(f_names) > 0L
    C <- NULL
    f_tilde_offsets <- integer(length(f_names))
    f_free_offsets <- integer(length(f_names))
    if (has_facets) {
        C <- matrix(0, nrow = n_tilde, ncol = n_free)
        C[cbind(seq_len(n_beta), seq_len(n_beta))] <- 1
        off_t <- n_beta
        off_f <- n_beta
        for (i in seq_along(f_names)) {
            B <- dat$facet_basis[[f_names[i]]]
            f_tilde_offsets[i] <- off_t
            f_free_offsets[i] <- off_f
            C[off_t + seq_len(nrow(B)), off_f + seq_len(ncol(B))] <- B
            off_t <- off_t + nrow(B)
            off_f <- off_f + ncol(B)
        }
        if (n_step > 0L) {
            C[cbind(off_t + seq_len(n_step), off_f + seq_len(n_step))] <- 1
        }
    }
    obs_ord <- order(dat$j)
    person_ptr <- as.integer(c(0L, cumsum(tabulate(dat$j, nbins = dat$J))))
    r_obs0 <- as.integer(dat$r[obs_ord] - 1L)
    x_obs0 <- as.integer(dat$x[obs_ord] - 1L)
    n_loc_blocks <- 1L + length(f_names)
    loc_cols <- vector("list", n_loc_blocks)
    beta_idx_obs <- if (time_mode == "dynamic")
        (dat$r - 1L) * Tn + dat$t
    else dat$r
    loc_cols[[1L]] <- as.integer(beta_idx_obs[obs_ord] - 1L)
    if (has_facets) {
        for (i in seq_along(f_names)) {
            fc <- f_names[i]
            loc_cols[[1L + i]] <- as.integer(f_tilde_offsets[i] + dat$facet_idx[[fc]][obs_ord] - 1L)
        }
    }
    loc_idx0 <- as.integer(t(do.call(cbind, loc_cols)))
    u_tilde_offset <- as.integer(n_beta + n_facet_tilde)
    list(time_mode = time_mode, step_mode = step_mode, n_u = n_u, n_beta = n_beta, n_facet_tilde = n_facet_tilde,
        n_facet_free = n_facet_free, n_step = n_step, n_tilde = n_tilde, n_free = n_free, has_facets = has_facets,
        f_names = f_names, f_sizes_tilde = f_sizes_tilde, f_sizes_free = f_sizes_free, f_tilde_offsets = f_tilde_offsets,
        f_free_offsets = f_free_offsets, C = C, person_ptr = person_ptr, r_obs0 = r_obs0, x_obs0 = x_obs0,
        loc_idx0 = loc_idx0, n_loc_blocks = n_loc_blocks, u_tilde_offset = u_tilde_offset)
}

.drift_mml_eval <- function(xi, dat, quad, layout, D = 1, need_hessian = TRUE) {
    tilde_psi <- if (layout$has_facets)
        as.numeric(layout$C %*% xi)
    else as.numeric(xi)
    res <- mfrm_cpp_drift_mml_eval(tilde_psi = tilde_psi, person_ptr = layout$person_ptr, r_obs0 = layout$r_obs0,
        x_obs0 = layout$x_obs0, loc_idx0 = layout$loc_idx0, n_loc_blocks = layout$n_loc_blocks, R = dat$R,
        K = dat$K, u_tilde_offset = layout$u_tilde_offset, step_rater = identical(layout$step_mode, "rater"),
        quad_nodes = as.numeric(quad$nodes), quad_weights = as.numeric(quad$weights), D = as.numeric(D),
        need_hessian = isTRUE(need_hessian))
    grad <- if (layout$has_facets)
        as.numeric(crossprod(layout$C, res$grad_tilde))
    else res$grad_tilde
    if (!isTRUE(need_hessian)) {
        return(list(ll = res$ll, grad = grad, pi_jm = res$pi_jm))
    }
    H_mml <- if (layout$has_facets) {
        crossprod(layout$C, res$H_mml_tilde %*% layout$C)
    }
    else {
        res$H_mml_tilde
    }
    H_mml <- 0.5 * (H_mml + t(H_mml))
    list(ll = res$ll, grad = grad, H_mml = H_mml, pi_jm = res$pi_jm, person_eap = res$person_eap, person_psd = res$person_psd)
}

.drift_prior_penalty <- function(xi, dat, layout, prior_mode = c("markov", "iid"), log_sigma_r = NULL,
    beta1_sd = 1, step_sd = 1) {
    prior_mode <- match.arg(prior_mode)
    R <- dat$R
    Tn <- dat$T
    n_u <- layout$n_u
    n_beta <- layout$n_beta
    n_facet_free <- layout$n_facet_free
    n_step <- layout$n_step
    n_free <- layout$n_free
    P <- matrix(0, nrow = n_free, ncol = n_free)
    log_norm <- 0
    if (layout$time_mode == "static" || prior_mode == "iid") {
        diag(P)[seq_len(n_beta)] <- 1/(beta1_sd^2)
        log_norm <- log_norm - n_beta * log(beta1_sd)
    }
    else {
        D1 <- matrix(0, nrow = Tn - 1L, ncol = Tn)
        for (t in seq_len(Tn - 1L)) {
            D1[t, t] <- -1
            D1[t, t + 1L] <- 1
        }
        DtD <- crossprod(D1)
        e1e1 <- matrix(0, Tn, Tn)
        e1e1[1L, 1L] <- 1/(beta1_sd^2)
        for (r in seq_len(R)) {
            idx <- (r - 1L) * Tn + seq_len(Tn)
            prec_rw <- exp(-2 * log_sigma_r[r])
            P[idx, idx] <- e1e1 + prec_rw * DtD
            log_norm <- log_norm - log(beta1_sd) - (Tn - 1L) * log_sigma_r[r]
        }
    }
    if (n_u > 0L) {
        A_u <- (diag(n_u) + matrix(1, n_u, n_u))/(step_sd^2)
        u_off <- n_beta + n_facet_free
        if (layout$step_mode == "rater") {
            for (r in seq_len(R)) {
                idx <- u_off + (r - 1L) * n_u + seq_len(n_u)
                P[idx, idx] <- A_u
            }
            log_norm <- log_norm + 0.5 * R * as.numeric(determinant(A_u, logarithm = TRUE)$modulus)
        }
        else {
            idx <- u_off + seq_len(n_u)
            P[idx, idx] <- A_u
            log_norm <- log_norm + 0.5 * as.numeric(determinant(A_u, logarithm = TRUE)$modulus)
        }
    }
    lp <- log_norm - 0.5 * sum(xi * as.numeric(P %*% xi))
    grad_p <- -as.numeric(P %*% xi)
    list(lp = lp, grad = grad_p, P = P)
}

.drift_solve_inner <- function(xi_init, dat, quad, layout, D = 1, prior_mode = "markov", log_sigma_r = NULL,
    beta1_sd = 1, step_sd = 1, tol = 1e-09, maxit = 45L) {
    xi <- xi_init
    pen <- .drift_prior_penalty(xi, dat, layout, prior_mode, log_sigma_r, beta1_sd, step_sd)
    ev <- .drift_mml_eval(xi, dat, quad, layout, D = D, need_hessian = TRUE)
    obj <- ev$ll + pen$lp
    g <- ev$grad + pen$grad
    inner_it <- 0L
    for (it in seq_len(maxit)) {
        if (max(abs(g)) < tol)
            break
        inner_it <- it
        H_tot <- ev$H_mml + pen$P
        eig_vals <- eigen(H_tot, symmetric = TRUE, only.values = TRUE)$values
        eig_min <- min(eig_vals)
        ridge <- if (!is.finite(eig_min) || eig_min < 1e-06)
            max(1e-04, 1e-06 - eig_min + 1e-04)
        else 0
        if (ridge > 0)
            diag(H_tot) <- diag(H_tot) + ridge
        step <- tryCatch(solve(H_tot, g), error = function(e) qr.solve(H_tot, g))
        alpha <- 1
        accepted <- FALSE
        while (alpha >= 1e-06) {
            xi_try <- xi + alpha * step
            ev_try <- .drift_mml_eval(xi_try, dat, quad, layout, D = D, need_hessian = FALSE)
            pen_try <- .drift_prior_penalty(xi_try, dat, layout, prior_mode, log_sigma_r, beta1_sd, step_sd)
            obj_try <- ev_try$ll + pen_try$lp
            if (is.finite(obj_try) && obj_try >= obj - 1e-12) {
                xi <- xi_try
                accepted <- TRUE
                break
            }
            alpha <- alpha * 0.5
        }
        if (!accepted)
            break
        pen <- .drift_prior_penalty(xi, dat, layout, prior_mode, log_sigma_r, beta1_sd, step_sd)
        ev <- .drift_mml_eval(xi, dat, quad, layout, D = D, need_hessian = TRUE)
        obj <- ev$ll + pen$lp
        g <- ev$grad + pen$grad
        if (max(abs(alpha * step)) < tol * 0.1 && max(abs(g)) < tol * 10)
            break
    }
    H_tot <- 0.5 * ((ev$H_mml + pen$P) + t(ev$H_mml + pen$P))
    eig_final <- eigen(H_tot, symmetric = TRUE, only.values = TRUE)$values
    info_positive <- all(is.finite(eig_final)) && min(eig_final) > max(eig_final) * 1e-10
    chol_H <- tryCatch(chol(H_tot), error = function(e) NULL)
    if (!is.null(chol_H)) {
        log_det_H <- 2 * sum(log(diag(chol_H)))
        V_xi <- chol2inv(chol_H)
    }
    else {
        info_positive <- FALSE
        eig_sym <- eigen(H_tot, symmetric = TRUE)
        pos_v <- pmax(eig_sym$values, 1e-06)
        log_det_H <- sum(log(pos_v))
        V_xi <- eig_sym$vectors %*% (t(eig_sym$vectors)/pos_v)
    }
    list(xi = xi, obj = obj, ll = ev$ll, grad = g, grad_mml = ev$grad, grad_max = max(abs(g)), iterations = inner_it,
        info_positive = info_positive, H_mml = ev$H_mml, P = pen$P, H_tot = H_tot, V_xi = V_xi, log_det_H = log_det_H,
        person_eap = ev$person_eap, person_psd = ev$person_psd, pi_jm = ev$pi_jm)
}

fit_mfrm_drift <- function(data, person, rater, time, score, facets = character(), score_levels = NULL,
    model = c("PCM", "RSM"), drift_mode = c("markov_hetero", "markov_homo", "independent", "static"),
    eb_method = c("em_laplace", "full_laplace"), mu_sigma = -2, delta_sigma = 1, beta1_sd = 1, step_sd = 1,
    person_sd = 1, logistic_scale = 1, quad_points = 81L, maxit = 60L, tol = 1e-06, missing = c("fail",
        "omit"), drift_threshold = 0.3, flag_se_ratio = 2) {
    t_start <- proc.time()
    model <- match.arg(model)
    drift_mode <- match.arg(drift_mode)
    eb_method <- match.arg(eb_method)
    missing <- match.arg(missing)
    for (nm_val in list(list("delta_sigma", delta_sigma, 0), list("beta1_sd", beta1_sd, 0), list("step_sd",
        step_sd, 0), list("person_sd", person_sd, 0), list("logistic_scale", logistic_scale, 0), list("tol",
        tol, 0))) {
        v <- nm_val[[2L]]
        if (!is.numeric(v) || is.complex(v) || length(v) != 1L || !is.finite(v) || v <= nm_val[[3L]]) {
            stop(sprintf("`%s` must be a single finite positive number.", nm_val[[1L]]), call. = FALSE)
        }
    }
    if (!is.numeric(mu_sigma) || is.complex(mu_sigma) || length(mu_sigma) != 1L || !is.finite(mu_sigma)) {
        stop("`mu_sigma` must be a single finite numeric value.", call. = FALSE)
    }
    if (!is.numeric(quad_points) || is.complex(quad_points) || length(quad_points) != 1L || !is.finite(quad_points) ||
        quad_points != floor(quad_points) || quad_points < 21L || quad_points > 201L) {
        stop("`quad_points` must be an integer between 21 and 201.", call. = FALSE)
    }
    if (!is.numeric(maxit) || is.complex(maxit) || length(maxit) != 1L || !is.finite(maxit) || maxit !=
        floor(maxit) || maxit < 5L || maxit > 300L) {
        stop("`maxit` must be an integer between 5 and 300.", call. = FALSE)
    }
    if (!is.numeric(drift_threshold) || length(drift_threshold) != 1L || !is.finite(drift_threshold) ||
        drift_threshold < 0) {
        stop("`drift_threshold` must be a single non-negative number.", call. = FALSE)
    }
    if (!is.numeric(flag_se_ratio) || length(flag_se_ratio) != 1L || !is.finite(flag_se_ratio) || flag_se_ratio <
        0) {
        stop("`flag_se_ratio` must be a single non-negative number.", call. = FALSE)
    }
    quad_points <- as.integer(quad_points)
    maxit <- as.integer(maxit)
    step_mode <- if (identical(model, "PCM"))
        "rater"
    else "shared"
    time_mode <- if (identical(drift_mode, "static"))
        "static"
    else "dynamic"
    dat <- .drift_prepare_input(data = data, person = person, rater = rater, time = time, score = score,
        facets = facets, score_levels = score_levels, missing = missing)
    layout <- .drift_layout(dat, time_mode = time_mode, step_mode = step_mode)
    quad <- .drift_quadrature_normal(n_quad = quad_points, mean = 0, sd = person_sd)
    est <- .drift_estimate_outer(dat = dat, quad = quad, layout = layout, drift_mode = drift_mode, eb_method = eb_method,
        mu_sigma = mu_sigma, delta_sigma = delta_sigma, beta1_sd = beta1_sd, step_sd = step_sd, logistic_scale = logistic_scale,
        maxit = maxit, tol = tol)
    inner <- est$inner
    check_quad_points <- quad_points + 40L
    quad_fine <- .drift_quadrature_normal(n_quad = check_quad_points, mean = 0, sd = person_sd)
    ev_fine <- .drift_mml_eval(inner$xi, dat, quad_fine, layout, D = logistic_scale, need_hessian = FALSE)
    ll_diff <- abs(ev_fine$ll - inner$ll)
    grad_diff <- max(abs(ev_fine$grad - inner$grad_mml))
    quadrature_stable <- isTRUE(ll_diff < 1e-05) && isTRUE(grad_diff < 1e-04)
    if (!quadrature_stable) {
        warning(sprintf("Person quadrature is not stable at `quad_points = %d` (|dLogLik| = %.2e, max |dGrad| = %.2e). Refit with larger `quad_points`.",
            quad_points, ll_diff, grad_diff), call. = FALSE)
    }
    info_ready <- isTRUE(inner$info_positive) && !dat$fixed_separated
    numerical_ready <- isTRUE(est$outer_converged) && isTRUE(inner$grad_max < 1e-04) && info_ready &&
        quadrature_stable
    if (!est$outer_converged) {
        warning("Outer Empirical Bayes iteration reached `maxit` before `tol` convergence; inspect `fit$checks`.",
            call. = FALSE)
    }
    if (dat$fixed_separated) {
        warning("One or more fixed-facet levels have only extreme or zero observations; fixed-facet MLE is separated.",
            call. = FALSE)
    }
    tabs <- .drift_build_tables(dat = dat, layout = layout, quad = quad, inner = inner, est = est, model = model,
        drift_mode = drift_mode, step_mode = step_mode, logistic_scale = logistic_scale, drift_threshold = drift_threshold,
        flag_se_ratio = flag_se_ratio)
    unobs_rc <- if (identical(model, "PCM"))
        sum(dat$rc_counts == 0L)
    else sum(colSums(dat$rc_counts) == 0L)
    unobs_rt <- sum(dat$rt_counts == 0L)
    checks <- data.frame(NumericalReady = numerical_ready, InformationPositive = info_ready, HyperCurvaturePositive = est$hyper_positive,
        OuterConverged = est$outer_converged, OuterIterations = as.integer(est$outer_iterations), MaxGradient = inner$grad_max,
        QuadraturePoints = quad_points, CheckPoints = check_quad_points, LogLikDifference = ll_diff,
        GradientDifference = grad_diff, QuadratureStable = quadrature_stable, FixedFacetSeparation = dat$fixed_separated,
        BridgePersons = dat$bridge_persons, WaveComponents = dat$wave_components, StationaryCohortRequired = dat$bridge_persons ==
            0L || dat$wave_components > 1L, UnobservedRaterWaveCells = as.integer(unobs_rt), UnobservedStepCells = as.integer(unobs_rc),
        stringsAsFactors = FALSE)
    elapsed_sec <- as.numeric((proc.time() - t_start)[3L])
    model_summary <- data.frame(Model = model, DriftMode = drift_mode, EBMethod = if (drift_mode %in%
        c("markov_hetero", "markov_homo"))
        eb_method
    else "fixed_prior", NObs = dat$N, NOmitted = dat$N_omitted, NPersons = dat$J, NRaters = dat$R, NWaves = dat$T,
        NCategories = dat$K, NFixedFacets = length(dat$active_facets), LogisticScale = logistic_scale,
        PersonSD = person_sd, MML_LogLik = inner$ll, LaplaceLogPost = est$laplace_obj, NominalParameters = as.integer(layout$n_free +
            tabs$n_hyper), EffectiveDF = tabs$eff_df, EffectiveAIC = -2 * inner$ll + 2 * tabs$eff_df,
        EAPReliability = tabs$eap_rel, FlaggedRaters = sum(tabs$rater_drift_table$Flagged), ElapsedSec = elapsed_sec,
        stringsAsFactors = FALSE)
    xi_hat <- inner$xi
    V_total <- est$V_total
    V_s <- est$V_s
    coord_names <- .drift_coordinate_names(dat, layout)
    names(xi_hat) <- coord_names
    dimnames(inner$V_xi) <- list(coord_names, coord_names)
    dimnames(V_total) <- list(coord_names, coord_names)
    if (!all(is.na(V_s)))
        dimnames(V_s) <- list(dat$levels$rater, dat$levels$rater)
    structure(list(model_summary = model_summary, severity_table = tabs$severity_table, rater_drift_table = tabs$rater_drift_table,
        step_table = tabs$step_table, facet_table = tabs$facet_table, person_table = tabs$person_table,
        checks = checks, matrices = tabs$matrices, coefficients = xi_hat, vcov_conditional = inner$V_xi,
        vcov_total = V_total, vcov_hyper = V_s, posterior = list(quad_nodes = quad$nodes, quad_weights = quad$weights,
            person_weights = inner$pi_jm), input = list(columns = dat$columns, levels = dat$levels, active_facets = dat$active_facets,
            facet_basis = dat$facet_basis, score_levels = dat$score_levels, n_raw = dat$N_raw, n_used = dat$N,
            n_omitted = dat$N_omitted), settings = list(model = model, drift_mode = drift_mode, eb_method = eb_method,
            mu_sigma = mu_sigma, delta_sigma = delta_sigma, beta1_sd = beta1_sd, step_sd = step_sd, person_sd = person_sd,
            logistic_scale = logistic_scale, quad_points = quad_points, maxit = maxit, tol = tol, missing = missing,
            drift_threshold = drift_threshold, flag_se_ratio = flag_se_ratio)), class = "mfrm_drift")
}

.drift_estimate_outer <- function(dat, quad, layout, drift_mode, eb_method, mu_sigma, delta_sigma, beta1_sd,
    step_sd, logistic_scale, maxit, tol) {
    R <- dat$R
    Tn <- dat$T
    n_facet_free <- layout$n_facet_free
    xi_init <- numeric(layout$n_free)
    outer_converged <- TRUE
    outer_iterations <- 1L
    hyper_positive <- TRUE
    if (drift_mode %in% c("static", "independent")) {
        inner <- .drift_solve_inner(xi_init = xi_init, dat = dat, quad = quad, layout = layout, D = logistic_scale,
            prior_mode = "iid", log_sigma_r = NULL, beta1_sd = beta1_sd, step_sd = step_sd, tol = 1e-10)
        return(list(inner = inner, outer_converged = TRUE, outer_iterations = 1L, hyper_positive = TRUE,
            log_sigma_hat = rep(NA_real_, R), sigma_map = rep(NA_real_, R), sigma_eap = rep(NA_real_,
                R), V_s = matrix(NA_real_, R, R), V_total = inner$V_xi, df_hyper = 0, laplace_obj = inner$obj -
                0.5 * inner$log_det_H + 0.5 * n_facet_free * log(2 * pi)))
    }
    D1 <- matrix(0, nrow = Tn - 1L, ncol = Tn)
    for (t in seq_len(Tn - 1L)) {
        D1[t, t] <- -1
        D1[t, t + 1L] <- 1
    }
    DtD <- crossprod(D1)
    s_vec <- rep(mu_sigma, R)
    xi_cur <- xi_init
    outer_converged <- FALSE
    for (out_it in seq_len(maxit)) {
        outer_iterations <- out_it
        inner <- .drift_solve_inner(xi_init = xi_cur, dat = dat, quad = quad, layout = layout, D = logistic_scale,
            prior_mode = "markov", log_sigma_r = s_vec, beta1_sd = beta1_sd, step_sd = step_sd, tol = 1e-09)
        xi_cur <- inner$xi
        S_r <- numeric(R)
        for (r in seq_len(R)) {
            idx <- (r - 1L) * Tn + seq_len(Tn)
            b_r <- xi_cur[idx]
            V_rr <- inner$V_xi[idx, idx, drop = FALSE]
            S_r[r] <- sum(diff(b_r)^2) + sum(DtD * V_rr)
        }
        s_new <- s_vec
        if (identical(drift_mode, "markov_hetero")) {
            grad_s <- -(Tn - 1L) + S_r * exp(-2 * s_vec) - (s_vec - mu_sigma)/(delta_sigma^2)
            if (max(abs(grad_s)) < tol) {
                outer_converged <- TRUE
                break
            }
            if (out_it <= 2L) {
                for (r in seq_len(R)) {
                  sr <- s_vec[r]
                  for (nit in 1:12) {
                    f1 <- -(Tn - 1L) + S_r[r] * exp(-2 * sr) - (sr - mu_sigma)/(delta_sigma^2)
                    f2 <- -2 * S_r[r] * exp(-2 * sr) - 1/(delta_sigma^2)
                    sr <- sr - f1/f2
                  }
                  s_new[r] <- sr
                }
            }
            else {
                h_diag <- numeric(R)
                for (r in seq_len(R)) {
                  idx <- (r - 1L) * Tn + seq_len(Tn)
                  V_rr <- inner$V_xi[idx, idx, drop = FALSE]
                  b_r <- xi_cur[idx]
                  A_r <- exp(-2 * s_vec[r]) * DtD
                  dphi_r <- 2 * as.numeric(A_r %*% b_r)
                  mode_term <- sum(dphi_r * as.numeric(V_rr %*% dphi_r))
                  tr_term <- 2 * sum((A_r %*% V_rr) * (V_rr %*% A_r))
                  h_r <- 2 * S_r[r] * exp(-2 * s_vec[r]) + 1/(delta_sigma^2) - mode_term - tr_term
                  em_curv <- 2 * S_r[r] * exp(-2 * s_vec[r]) + 1/(delta_sigma^2)
                  h_diag[r] <- max(h_r, 0.2 * em_curv)
                }
                step_s <- pmax(pmin(grad_s/h_diag, 1.5), -1.5)
                s_new <- s_vec + step_s
            }
        }
        else {
            sr <- s_vec[1L]
            S_tot <- sum(S_r)
            f1_cur <- -R * (Tn - 1L) + S_tot * exp(-2 * sr) - (sr - mu_sigma)/(delta_sigma^2)
            if (abs(f1_cur) < tol) {
                outer_converged <- TRUE
                break
            }
            for (nit in 1:15) {
                f1 <- -R * (Tn - 1L) + S_tot * exp(-2 * sr) - (sr - mu_sigma)/(delta_sigma^2)
                f2 <- -2 * S_tot * exp(-2 * sr) - 1/(delta_sigma^2)
                sr <- sr - f1/f2
            }
            s_new <- rep(sr, R)
        }
        max_ds <- max(abs(s_new - s_vec))
        s_vec <- s_new
        if (max_ds < tol) {
            outer_converged <- TRUE
            break
        }
    }
    if (identical(eb_method, "full_laplace")) {
        env_cache <- new.env(parent = emptyenv())
        env_cache$xi <- xi_cur
        eval_neg_laplace <- function(par_s) {
            s_full <- if (identical(drift_mode, "markov_hetero"))
                par_s
            else rep(par_s[1L], R)
            inn <- .drift_solve_inner(xi_init = env_cache$xi, dat = dat, quad = quad, layout = layout,
                D = logistic_scale, prior_mode = "markov", log_sigma_r = s_full, beta1_sd = beta1_sd,
                step_sd = step_sd, tol = 1e-10)
            env_cache$xi <- inn$xi
            hyper_lp <- if (identical(drift_mode, "markov_hetero")) {
                -0.5 * sum((par_s - mu_sigma)^2)/(delta_sigma^2) - R * log(delta_sigma)
            }
            else {
                -0.5 * ((par_s[1L] - mu_sigma)^2)/(delta_sigma^2) - log(delta_sigma)
            }
            -(inn$obj - 0.5 * inn$log_det_H + hyper_lp)
        }
        par_init <- if (identical(drift_mode, "markov_hetero"))
            s_vec
        else s_vec[1L]
        opt_lap <- stats::optim(par_init, eval_neg_laplace, method = "BFGS", control = list(reltol = 1e-08,
            maxit = 40L))
        s_vec <- if (identical(drift_mode, "markov_hetero"))
            opt_lap$par
        else rep(opt_lap$par[1L], R)
        outer_converged <- outer_converged || identical(opt_lap$convergence, 0L)
    }
    inner <- .drift_solve_inner(xi_init = xi_cur, dat = dat, quad = quad, layout = layout, D = logistic_scale,
        prior_mode = "markov", log_sigma_r = s_vec, beta1_sd = beta1_sd, step_sd = step_sd, tol = 1e-10)
    d2Phi_dxi_ds <- matrix(0, nrow = layout$n_free, ncol = R)
    S_r <- numeric(R)
    for (r in seq_len(R)) {
        idx <- (r - 1L) * Tn + seq_len(Tn)
        b_r <- inner$xi[idx]
        V_rr <- inner$V_xi[idx, idx, drop = FALSE]
        S_r[r] <- sum(diff(b_r)^2) + sum(DtD * V_rr)
        d2Phi_dxi_ds[idx, r] <- 2 * exp(-2 * s_vec[r]) * as.numeric(DtD %*% b_r)
    }
    J_xi_s <- inner$V_xi %*% d2Phi_dxi_ds
    H_s_exact <- matrix(0, R, R)
    for (r1 in seq_len(R)) {
        idx1 <- (r1 - 1L) * Tn + seq_len(Tn)
        A1 <- exp(-2 * s_vec[r1]) * DtD
        for (r2 in r1:R) {
            idx2 <- (r2 - 1L) * Tn + seq_len(Tn)
            A2 <- exp(-2 * s_vec[r2]) * DtD
            V12 <- inner$V_xi[idx1, idx2, drop = FALSE]
            tr_term <- 2 * sum((A1 %*% V12) * (V12 %*% A2))
            val <- -as.numeric(crossprod(d2Phi_dxi_ds[, r1], J_xi_s[, r2])) - tr_term
            if (r1 == r2) {
                val <- val + 2 * S_r[r1] * exp(-2 * s_vec[r1]) + 1/(delta_sigma^2)
            }
            H_s_exact[r1, r2] <- val
            H_s_exact[r2, r1] <- val
        }
    }
    if (identical(drift_mode, "markov_hetero")) {
        eig_hs <- eigen(H_s_exact, symmetric = TRUE)
        hyper_positive <- all(is.finite(eig_hs$values)) && min(eig_hs$values) > 0
        vals_pos <- pmax(eig_hs$values, 0.25/(delta_sigma^2))
        V_s <- eig_hs$vectors %*% (t(eig_hs$vectors)/vals_pos)
        log_det_Vs <- -sum(log(vals_pos))
        df_hyper <- sum(pmax(0, 1 - diag(V_s)/(delta_sigma^2)))
        hyper_quad <- -0.5 * sum((s_vec - mu_sigma)^2)/(delta_sigma^2)
        occam_hyper <- 0.5 * log_det_Vs - R * log(delta_sigma)
    }
    else {
        h_1_raw <- sum(H_s_exact)
        hyper_positive <- is.finite(h_1_raw) && h_1_raw > 0
        h_1 <- max(h_1_raw, 0.25/(delta_sigma^2))
        V_s <- matrix(1/h_1, R, R)
        log_det_Vs <- -log(h_1)
        df_hyper <- max(0, 1 - (1/h_1)/(delta_sigma^2))
        hyper_quad <- -0.5 * ((s_vec[1L] - mu_sigma)^2)/(delta_sigma^2)
        occam_hyper <- 0.5 * log_det_Vs - log(delta_sigma)
    }
    V_total <- inner$V_xi + J_xi_s %*% V_s %*% t(J_xi_s)
    V_total <- 0.5 * (V_total + t(V_total))
    laplace_obj <- inner$obj - 0.5 * inner$log_det_H + 0.5 * n_facet_free * log(2 * pi) + hyper_quad +
        occam_hyper
    list(inner = inner, outer_converged = outer_converged, outer_iterations = outer_iterations, hyper_positive = hyper_positive,
        log_sigma_hat = s_vec, sigma_map = exp(s_vec), sigma_eap = exp(s_vec + 0.5 * pmax(diag(V_s),
            0)), V_s = V_s, V_total = V_total, df_hyper = df_hyper, laplace_obj = laplace_obj)
}

.drift_build_tables <- function(dat, layout, quad, inner, est, model, drift_mode, step_mode, logistic_scale,
    drift_threshold, flag_se_ratio) {
    R <- dat$R
    Tn <- dat$T
    K <- dat$K
    n_u <- layout$n_u
    n_beta <- layout$n_beta
    n_facet_free <- layout$n_facet_free
    V_total <- est$V_total
    V_s <- est$V_s
    xi_hat <- inner$xi
    se_cond_all <- sqrt(pmax(diag(inner$V_xi), 0))
    se_tot_all <- sqrt(pmax(diag(V_total), 0))
    rater_levels <- dat$levels$rater
    time_levels <- dat$levels$time
    score_lvls <- dat$score_levels
    if (layout$time_mode == "dynamic") {
        beta_rt_mat <- matrix(xi_hat[seq_len(n_beta)], nrow = R, ncol = Tn, byrow = TRUE, dimnames = list(rater_levels,
            time_levels))
        se_beta_cond <- matrix(se_cond_all[seq_len(n_beta)], nrow = R, ncol = Tn, byrow = TRUE, dimnames = list(rater_levels,
            time_levels))
        se_beta_tot <- matrix(se_tot_all[seq_len(n_beta)], nrow = R, ncol = Tn, byrow = TRUE, dimnames = list(rater_levels,
            time_levels))
    }
    else {
        beta_rt_mat <- matrix(xi_hat[seq_len(R)], nrow = R, ncol = Tn, dimnames = list(rater_levels,
            time_levels))
        se_beta_cond <- matrix(se_cond_all[seq_len(R)], nrow = R, ncol = Tn, dimnames = list(rater_levels,
            time_levels))
        se_beta_tot <- matrix(se_tot_all[seq_len(R)], nrow = R, ncol = Tn, dimnames = list(rater_levels,
            time_levels))
    }
    sev_rows <- vector("list", R)
    rater_rows <- vector("list", R)
    log_sigma_se <- if (all(is.na(V_s)))
        rep(NA_real_, R)
    else sqrt(pmax(diag(V_s), 0))
    for (r in seq_len(R)) {
        b_r <- beta_rt_mat[r, ]
        chg <- b_r - b_r[1L]
        chg_se <- numeric(Tn)
        chg_se_cond <- numeric(Tn)
        if (layout$time_mode == "dynamic") {
            idx_r <- (r - 1L) * Tn + seq_len(Tn)
            Vr_tot <- V_total[idx_r, idx_r, drop = FALSE]
            Vr_cnd <- inner$V_xi[idx_r, idx_r, drop = FALSE]
            for (t in seq_len(Tn)) {
                if (t == 1L) {
                  chg_se[t] <- 0
                  chg_se_cond[t] <- 0
                }
                else {
                  v_diff <- Vr_tot[t, t] + Vr_tot[1L, 1L] - 2 * Vr_tot[t, 1L]
                  vc_diff <- Vr_cnd[t, t] + Vr_cnd[1L, 1L] - 2 * Vr_cnd[t, 1L]
                  chg_se[t] <- sqrt(max(v_diff, 0))
                  chg_se_cond[t] <- sqrt(max(vc_diff, 0))
                }
            }
        }
        z_chg <- ifelse(chg_se > 0, abs(chg)/chg_se, 0)
        flagged_wave <- (seq_len(Tn) > 1L) & (abs(chg) >= drift_threshold) & (if (flag_se_ratio > 0)
            z_chg >= flag_se_ratio
        else TRUE)
        n_rt <- as.integer(dat$rt_counts[r, ])
        sev_rows[[r]] <- data.frame(Rater = rep(rater_levels[r], Tn), Time = time_levels, WaveIndex = seq_len(Tn),
            NObs = n_rt, Observed = n_rt > 0L, Estimate = as.numeric(b_r), ConditionalSE = as.numeric(se_beta_cond[r,
                ]), TotalSE = as.numeric(se_beta_tot[r, ]), ChangeFromBaseline = as.numeric(chg), ChangeConditionalSE = as.numeric(chg_se_cond),
            ChangeSE = as.numeric(chg_se), AbsZ = as.numeric(z_chg), FlaggedDrift = as.logical(flagged_wave),
            stringsAsFactors = FALSE)
        rater_rows[[r]] <- data.frame(Rater = rater_levels[r], NObs = sum(n_rt), ActiveWaves = sum(n_rt >
            0L), LogSigma = est$log_sigma_hat[r], LogSigmaSE = log_sigma_se[r], DriftSD_MAP = est$sigma_map[r],
            DriftSD_EAP = est$sigma_eap[r], BaselineSeverity = b_r[1L], FinalSeverity = b_r[Tn], NetDrift = chg[Tn],
            NetDriftSE = chg_se[Tn], MaxAbsChange = max(abs(chg)), MaxAbsZ = max(z_chg), TrajectoryRange = diff(range(b_r)),
            Flagged = any(flagged_wave), stringsAsFactors = FALSE)
    }
    severity_table <- do.call(rbind, sev_rows)
    rownames(severity_table) <- NULL
    rater_drift_table <- do.call(rbind, rater_rows)
    rownames(rater_drift_table) <- NULL
    facet_rows <- list()
    if (layout$has_facets) {
        for (i in seq_along(layout$f_names)) {
            fc <- layout$f_names[i]
            B <- dat$facet_basis[[fc]]
            idx_free <- layout$f_free_offsets[i] + seq_len(layout$f_sizes_free[i])
            v_fc <- xi_hat[idx_free]
            est_fc <- as.numeric(B %*% v_fc)
            V_fc_cnd <- B %*% inner$V_xi[idx_free, idx_free, drop = FALSE] %*% t(B)
            V_fc_tot <- B %*% V_total[idx_free, idx_free, drop = FALSE] %*% t(B)
            lv_fc <- rownames(B)
            n_fc <- as.integer(table(factor(dat$facet_idx[[fc]], levels = seq_along(lv_fc))))
            facet_rows[[fc]] <- data.frame(Facet = rep(fc, length(lv_fc)), Level = lv_fc, NObs = n_fc,
                Estimate = est_fc, ConditionalSE = sqrt(pmax(diag(V_fc_cnd), 0)), TotalSE = sqrt(pmax(diag(V_fc_tot),
                  0)), Constraint = rep("sum_to_zero", length(lv_fc)), stringsAsFactors = FALSE)
        }
    }
    facet_table <- if (length(facet_rows) > 0L) {
        out_ft <- do.call(rbind, facet_rows)
        rownames(out_ft) <- NULL
        out_ft
    }
    else {
        data.frame(Facet = character(), Level = character(), NObs = integer(), Estimate = numeric(),
            ConditionalSE = numeric(), TotalSE = numeric(), Constraint = character(), stringsAsFactors = FALSE)
    }
    d_mat <- matrix(0, nrow = R, ncol = K - 1L, dimnames = list(rater_levels, paste0("Step_", seq.int(2L,
        K))))
    se_d_cond_mat <- matrix(0, nrow = R, ncol = K - 1L, dimnames = dimnames(d_mat))
    se_d_tot_mat <- matrix(0, nrow = R, ncol = K - 1L, dimnames = dimnames(d_mat))
    u_off <- n_beta + n_facet_free
    if (n_u > 0L) {
        B_step <- rbind(diag(n_u), rep(-1, n_u))
        if (identical(step_mode, "rater")) {
            for (r in seq_len(R)) {
                idx <- u_off + (r - 1L) * n_u + seq_len(n_u)
                u_r <- xi_hat[idx]
                d_mat[r, ] <- as.numeric(B_step %*% u_r)
                Vc_u <- B_step %*% inner$V_xi[idx, idx, drop = FALSE] %*% t(B_step)
                Vt_u <- B_step %*% V_total[idx, idx, drop = FALSE] %*% t(B_step)
                se_d_cond_mat[r, ] <- sqrt(pmax(diag(Vc_u), 0))
                se_d_tot_mat[r, ] <- sqrt(pmax(diag(Vt_u), 0))
            }
        }
        else {
            idx <- u_off + seq_len(n_u)
            u_vec <- xi_hat[idx]
            d_row <- as.numeric(B_step %*% u_vec)
            Vc_u <- B_step %*% inner$V_xi[idx, idx, drop = FALSE] %*% t(B_step)
            Vt_u <- B_step %*% V_total[idx, idx, drop = FALSE] %*% t(B_step)
            d_mat <- matrix(d_row, nrow = R, ncol = K - 1L, byrow = TRUE, dimnames = dimnames(d_mat))
            se_d_cond_mat <- matrix(sqrt(pmax(diag(Vc_u), 0)), nrow = R, ncol = K - 1L, byrow = TRUE,
                dimnames = dimnames(d_mat))
            se_d_tot_mat <- matrix(sqrt(pmax(diag(Vt_u), 0)), nrow = R, ncol = K - 1L, byrow = TRUE,
                dimnames = dimnames(d_mat))
        }
    }
    if (identical(step_mode, "rater")) {
        step_rows <- lapply(seq_len(R), function(r) {
            data.frame(Rater = rep(rater_levels[r], K - 1L), Step = paste0("d_", seq.int(2L, K)), FromCategory = score_lvls[seq_len(K -
                1L)], ToCategory = score_lvls[seq.int(2L, K)], CategoryCount = as.integer(dat$rc_counts[r,
                seq.int(2L, K)]), Estimate = as.numeric(d_mat[r, ]), ConditionalSE = as.numeric(se_d_cond_mat[r,
                ]), TotalSE = as.numeric(se_d_tot_mat[r, ]), stringsAsFactors = FALSE)
        })
        step_table <- do.call(rbind, step_rows)
    }
    else {
        tot_cat_counts <- colSums(dat$rc_counts)
        step_table <- data.frame(Rater = rep("(Shared)", K - 1L), Step = paste0("d_", seq.int(2L, K)),
            FromCategory = score_lvls[seq_len(K - 1L)], ToCategory = score_lvls[seq.int(2L, K)], CategoryCount = as.integer(tot_cat_counts[seq.int(2L,
                K)]), Estimate = as.numeric(d_mat[1L, ]), ConditionalSE = as.numeric(se_d_cond_mat[1L,
                ]), TotalSE = as.numeric(se_d_tot_mat[1L, ]), stringsAsFactors = FALSE)
    }
    rownames(step_table) <- NULL
    eff_df <- sum(inner$V_xi * inner$H_mml) + est$df_hyper
    n_hyper <- switch(drift_mode, markov_hetero = R, markov_homo = 1L, 0L)
    var_eap <- if (dat$J > 1L)
        stats::var(inner$person_eap)
    else 0
    mean_psd2 <- mean(inner$person_psd^2)
    eap_rel <- if ((var_eap + mean_psd2) > 0)
        var_eap/(var_eap + mean_psd2)
    else NA_real_
    fit_core <- list(beta_rt = beta_rt_mat, d_rm = d_mat, facet_table = facet_table, person_eap = inner$person_eap,
        person_psd = inner$person_psd, pi_jm = inner$pi_jm, quad = quad, score_levels = score_lvls, logistic_scale = logistic_scale,
        levels = dat$levels)
    fair_w1 <- .drift_compute_fair_scores(fit_core, ref_raters = seq_len(R), ref_waves = 1L)
    fair_all <- .drift_compute_fair_scores(fit_core, ref_raters = seq_len(R), ref_waves = seq_len(Tn))
    p_nobs <- as.integer(tabulate(dat$j, nbins = dat$J))
    p_mean_raw <- as.numeric(.drift_rowsum_by_id(as.numeric(dat$raw_score), dat$j, dat$J))/pmax(p_nobs,
        1L)
    person_table <- data.frame(Person = dat$levels$person, NObs = p_nobs, NWaves = dat$person_waves,
        MeanScore = p_mean_raw, ThetaEAP = inner$person_eap, ThetaPSD = inner$person_psd, FairM_Wave1 = fair_w1$post_expected,
        FairM_PanelMean = fair_all$post_expected, stringsAsFactors = FALSE)
    list(severity_table = severity_table, rater_drift_table = rater_drift_table, step_table = step_table,
        facet_table = facet_table, person_table = person_table, eff_df = eff_df, n_hyper = n_hyper, eap_rel = eap_rel,
        matrices = list(beta_rt = beta_rt_mat, se_beta_cond = se_beta_cond, se_beta_total = se_beta_tot,
            d_rm = d_mat, se_d_cond = se_d_cond_mat, se_d_total = se_d_tot_mat, rc_counts = dat$rc_counts,
            rt_counts = dat$rt_counts))
}

.drift_coordinate_names <- function(dat, layout) {
    r_lv <- dat$levels$rater
    t_lv <- dat$levels$time
    b_names <- if (layout$time_mode == "dynamic") {
        as.vector(t(outer(r_lv, t_lv, function(r, t) paste0("severity:", r, "@", t))))
    }
    else {
        paste0("severity:", r_lv)
    }
    f_names <- character()
    if (layout$has_facets) {
        for (fc in layout$f_names) {
            B <- dat$facet_basis[[fc]]
            f_names <- c(f_names, paste0("facet:", fc, ":", colnames(B)))
        }
    }
    s_names <- character()
    if (layout$n_u > 0L) {
        st_idx <- seq.int(2L, dat$K - 1L)
        if (identical(layout$step_mode, "rater")) {
            s_names <- as.vector(t(outer(r_lv, st_idx, function(r, m) paste0("step:", r, ":d", m))))
        }
        else {
            s_names <- paste0("step:d", st_idx)
        }
    }
    c(b_names, f_names, s_names)
}

.drift_compute_fair_scores <- function(fit_core, ref_raters, ref_waves) {
    K <- length(fit_core$score_levels)
    D <- fit_core$logistic_scale
    c_vec <- as.numeric(0:(K - 1L))
    s_vals <- as.numeric(fit_core$score_levels)
    eval_at_theta <- function(theta_vec) {
        n_t <- length(theta_vec)
        acc <- numeric(n_t)
        n_cells <- length(ref_raters) * length(ref_waves)
        for (r in ref_raters) {
            g_r <- c(0, cumsum(fit_core$d_rm[r, ]))
            for (w in ref_waves) {
                b_rw <- fit_core$beta_rt[r, w]
                eta <- matrix(0, nrow = n_t, ncol = K)
                for (k in 2:K) {
                  eta[, k] <- D * (c_vec[k] * (theta_vec - b_rw) - g_r[k])
                }
                mx <- apply(eta, 1L, max)
                pr <- exp(eta - mx)
                pr <- pr/rowSums(pr)
                acc <- acc + as.numeric(pr %*% s_vals)
            }
        }
        acc/n_cells
    }
    plug <- eval_at_theta(fit_core$person_eap)
    node_exp <- eval_at_theta(fit_core$quad$nodes)
    post_exp <- as.numeric(fit_core$pi_jm %*% node_exp)
    list(plugin = plug, post_expected = post_exp)
}

print.mfrm_drift <- function(x, digits = 3, ...) {
    ms <- x$model_summary
    cat("Markov Rater-Severity Drift MFRM (`mfrm_drift`)\n")
    cat(sprintf("  Model: %s | Drift mode: %s | Scale D = %s | Quad points: %d\n", ms$Model, ms$DriftMode,
        format(ms$LogisticScale, digits = 3), x$checks$QuadraturePoints))
    cat(sprintf("  Data: %d ratings | %d persons | %d raters | %d waves | %d categories (%d..%d)\n",
        ms$NObs, ms$NPersons, ms$NRaters, ms$NWaves, ms$NCategories, min(x$input$score_levels), max(x$input$score_levels)))
    cat(sprintf("  Fit : MML logLik = %.*f | Laplace logPost = %.*f | Eff.df = %.2f | EAP rel = %.*f\n",
        digits, ms$MML_LogLik, digits, ms$LaplaceLogPost, ms$EffectiveDF, digits, ms$EAPReliability))
    cat(sprintf("  Checks: NumericalReady = %s | QuadratureStable = %s | BridgePersons = %d | Flagged raters = %d/%d\n",
        x$checks$NumericalReady, x$checks$QuadratureStable, x$checks$BridgePersons, ms$FlaggedRaters,
        ms$NRaters))
    invisible(x)
}

summary.mfrm_drift <- function(object, level = 0.95, se_type = c("total", "conditional"), ...) {
    rlang::check_dots_empty()
    se_type <- match.arg(se_type)
    if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
        stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
    }
    ci_sev <- confint.mfrm_drift(object, parm = "severity", level = level, se_type = se_type)
    ci_sd <- if (object$settings$drift_mode %in% c("markov_hetero", "markov_homo")) {
        confint.mfrm_drift(object, parm = "drift_sd", level = level, se_type = se_type)
    }
    else {
        NULL
    }
    notes <- character()
    if (isTRUE(object$checks$StationaryCohortRequired)) {
        notes <- c(notes, "Cross-wave identification: waves are not fully linked by common bridge persons; common wave severity shifts rely on a stationary N(0, person_sd^2) cohort distribution.")
    }
    if (object$checks$UnobservedRaterWaveCells > 0L) {
        notes <- c(notes, sprintf("%d rater-by-wave cell(s) have 0 ratings (`Observed == FALSE` in `severity_table`); Markov estimates interpolate across adjacent waves.",
            object$checks$UnobservedRaterWaveCells))
    }
    if (identical(object$settings$model, "PCM") && object$checks$UnobservedStepCells > 0L) {
        notes <- c(notes, sprintf("%d rater-by-category cell(s) have 0 observed ratings; consider `model = 'RSM'` if category support is sparse.",
            object$checks$UnobservedStepCells))
    }
    structure(list(model_summary = object$model_summary, rater_drift = object$rater_drift_table, drift_sd_intervals = ci_sd,
        severity_intervals = ci_sev, facets = object$facet_table, checks = object$checks, level = level,
        se_type = se_type, notes = notes), class = "summary.mfrm_drift")
}

print.summary.mfrm_drift <- function(x, digits = 3, ...) {
    ms <- x$model_summary
    cat("Summary of Markov Rater-Severity Drift MFRM\n")
    cat(sprintf("Model: %s (%s) | N = %d | Persons = %d | Raters = %d | Waves = %d | Categories = %d\n",
        ms$Model, ms$DriftMode, ms$NObs, ms$NPersons, ms$NRaters, ms$NWaves, ms$NCategories))
    cat(sprintf("MML logLik = %.*f | Laplace logPost = %.*f | Effective df = %.2f | Effective AIC = %.*f | EAP rel = %.*f\n\n",
        digits, ms$MML_LogLik, digits, ms$LaplaceLogPost, ms$EffectiveDF, digits, ms$EffectiveAIC, digits,
        ms$EAPReliability))
    cat("Rater Drift Summary:\n")
    rd <- x$rater_drift
    num_cols <- vapply(rd, is.numeric, logical(1))
    rd_print <- rd
    for (nm in names(rd_print)[num_cols]) {
        if (!nm %in% c("NObs", "ActiveWaves")) {
            rd_print[[nm]] <- round(rd_print[[nm]], digits)
        }
    }
    print(rd_print, row.names = FALSE)
    if (nrow(x$facets) > 0L) {
        cat("\nFixed Facet Estimates (sum-to-zero):\n")
        ft_print <- x$facets
        for (nm in c("Estimate", "ConditionalSE", "TotalSE")) {
            ft_print[[nm]] <- round(ft_print[[nm]], digits)
        }
        print(ft_print, row.names = FALSE)
    }
    if (length(x$notes) > 0L) {
        cat("\nIdentification & Design Notes:\n")
        for (nt in x$notes) {
            cat(" - ", nt, "\n", sep = "")
        }
    }
    invisible(x)
}

coef.mfrm_drift <- function(object, parm = c("severity", "drift_sd", "contrasts", "steps", "facets",
    "all"), ...) {
    rlang::check_dots_empty()
    parm <- match.arg(parm)
    switch(parm, severity = {
        st <- object$severity_table
        stats::setNames(st$Estimate, paste0(st$Rater, "@", st$Time))
    }, drift_sd = {
        rd <- object$rater_drift_table
        stats::setNames(rd$DriftSD_EAP, rd$Rater)
    }, contrasts = {
        st <- object$severity_table[object$severity_table$WaveIndex > 1L, , drop = FALSE]
        stats::setNames(st$ChangeFromBaseline, paste0(st$Rater, ":", st$Time, "-", object$input$levels$time[1L]))
    }, steps = {
        sp <- object$step_table
        nm <- if (identical(object$settings$model, "PCM")) paste0(sp$Rater, ":", sp$Step) else sp$Step
        stats::setNames(sp$Estimate, nm)
    }, facets = {
        ft <- object$facet_table
        if (nrow(ft) == 0L) return(stats::setNames(numeric(), character()))
        stats::setNames(ft$Estimate, paste0(ft$Facet, ":", ft$Level))
    }, all = object$coefficients)
}

vcov.mfrm_drift <- function(object, se_type = c("total", "conditional", "hyper"), ...) {
    rlang::check_dots_empty()
    se_type <- match.arg(se_type)
    switch(se_type, total = object$vcov_total, conditional = object$vcov_conditional, hyper = object$vcov_hyper)
}

confint.mfrm_drift <- function(object, parm = c("severity", "contrasts", "drift_sd", "steps", "facets"),
    level = 0.95, se_type = c("total", "conditional"), ...) {
    rlang::check_dots_empty()
    parm <- match.arg(parm)
    se_type <- match.arg(se_type)
    if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
        stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
    }
    z <- stats::qnorm(0.5 + 0.5 * level)
    ready <- isTRUE(object$checks$NumericalReady) && isTRUE(object$checks$InformationPositive)
    if (identical(parm, "severity")) {
        st <- object$severity_table
        se <- if (identical(se_type, "total"))
            st$TotalSE
        else st$ConditionalSE
        if (!ready)
            se <- rep(NA_real_, length(se))
        return(data.frame(Parameter = "severity", Rater = st$Rater, Time = st$Time, WaveIndex = st$WaveIndex,
            Observed = st$Observed, Estimate = st$Estimate, SE = se, Lower = st$Estimate - z * se, Upper = st$Estimate +
                z * se, Level = level, SEType = se_type, stringsAsFactors = FALSE))
    }
    if (identical(parm, "contrasts")) {
        st <- object$severity_table[object$severity_table$WaveIndex > 1L, , drop = FALSE]
        se <- if (identical(se_type, "total"))
            st$ChangeSE
        else st$ChangeConditionalSE
        if (!ready)
            se <- rep(NA_real_, length(se))
        base_lbl <- object$input$levels$time[1L]
        return(data.frame(Parameter = "contrast", Rater = st$Rater, Time = st$Time, Baseline = rep(base_lbl,
            nrow(st)), Estimate = st$ChangeFromBaseline, SE = se, Lower = st$ChangeFromBaseline - z *
            se, Upper = st$ChangeFromBaseline + z * se, Level = level, SEType = se_type, stringsAsFactors = FALSE))
    }
    if (identical(parm, "drift_sd")) {
        rd <- object$rater_drift_table
        if (!object$settings$drift_mode %in% c("markov_hetero", "markov_homo")) {
            return(data.frame(Parameter = rep("drift_sd", nrow(rd)), Rater = rd$Rater, LogSigma = rep(NA_real_,
                nrow(rd)), LogSigmaSE = rep(NA_real_, nrow(rd)), DriftSD_MAP = rep(NA_real_, nrow(rd)),
                DriftSD_EAP = rep(NA_real_, nrow(rd)), Lower = rep(NA_real_, nrow(rd)), Upper = rep(NA_real_,
                  nrow(rd)), Level = level, stringsAsFactors = FALSE))
        }
        se_s <- if (ready && isTRUE(object$checks$HyperCurvaturePositive))
            rd$LogSigmaSE
        else rep(NA_real_, nrow(rd))
        return(data.frame(Parameter = rep("drift_sd", nrow(rd)), Rater = rd$Rater, LogSigma = rd$LogSigma,
            LogSigmaSE = se_s, DriftSD_MAP = rd$DriftSD_MAP, DriftSD_EAP = rd$DriftSD_EAP, Lower = exp(rd$LogSigma -
                z * se_s), Upper = exp(rd$LogSigma + z * se_s), Level = level, stringsAsFactors = FALSE))
    }
    if (identical(parm, "steps")) {
        sp <- object$step_table
        se <- if (identical(se_type, "total"))
            sp$TotalSE
        else sp$ConditionalSE
        if (!ready)
            se <- rep(NA_real_, length(se))
        return(data.frame(Parameter = "step", Rater = sp$Rater, Step = sp$Step, Estimate = sp$Estimate,
            SE = se, Lower = sp$Estimate - z * se, Upper = sp$Estimate + z * se, Level = level, SEType = se_type,
            stringsAsFactors = FALSE))
    }
    ft <- object$facet_table
    if (nrow(ft) == 0L) {
        return(data.frame(Parameter = character(), Facet = character(), LevelName = character(), Estimate = numeric(),
            SE = numeric(), Lower = numeric(), Upper = numeric(), Level = numeric(), SEType = character(),
            stringsAsFactors = FALSE))
    }
    se <- if (identical(se_type, "total"))
        ft$TotalSE
    else ft$ConditionalSE
    if (!ready)
        se <- rep(NA_real_, length(se))
    data.frame(Parameter = "facet", Facet = ft$Facet, LevelName = ft$Level, Estimate = ft$Estimate, SE = se,
        Lower = ft$Estimate - z * se, Upper = ft$Estimate + z * se, Level = level, SEType = se_type,
        stringsAsFactors = FALSE)
}

predict.mfrm_drift <- function(object, newdata = NULL, type = c("fair_score", "expected", "probabilities"),
    ref_time = NULL, ref_raters = NULL, ability = NULL, integration = c("posterior", "plugin"), ...) {
    rlang::check_dots_empty()
    if (!isTRUE(object$checks$NumericalReady) || !isTRUE(object$checks$InformationPositive)) {
        stop("Resolve numerical and information checks in `object$checks` before prediction.", call. = FALSE)
    }
    integration <- match.arg(integration)
    if (is.null(newdata)) {
        type <- if (missing(type))
            "fair_score"
        else match.arg(type)
        if (!identical(type, "fair_score")) {
            stop("Supply `newdata` when `type` is `'expected'` or `'probabilities'`.", call. = FALSE)
        }
        t_lvls <- object$input$levels$time
        r_lvls <- object$input$levels$rater
        if (is.null(ref_time)) {
            w_idx <- 1L
        }
        else if (is.numeric(ref_time)) {
            w_idx <- as.integer(ref_time)
            if (anyNA(w_idx) || any(w_idx < 1L | w_idx > length(t_lvls))) {
                stop("Numeric `ref_time` indices must lie within 1..", length(t_lvls), ".", call. = FALSE)
            }
        }
        else {
            w_idx <- match(as.character(ref_time), t_lvls)
            if (anyNA(w_idx)) {
                stop("Unknown `ref_time` wave label(s): ", paste(setdiff(as.character(ref_time), t_lvls),
                  collapse = ", "), ".", call. = FALSE)
            }
        }
        if (is.null(ref_raters)) {
            r_idx <- seq_along(r_lvls)
        }
        else {
            r_idx <- match(as.character(ref_raters), r_lvls)
            if (anyNA(r_idx)) {
                stop("Unknown `ref_raters` label(s): ", paste(setdiff(as.character(ref_raters), r_lvls),
                  collapse = ", "), ".", call. = FALSE)
            }
        }
        fit_core <- list(beta_rt = object$matrices$beta_rt, d_rm = object$matrices$d_rm, person_eap = object$person_table$ThetaEAP,
            pi_jm = object$posterior$person_weights, quad = list(nodes = object$posterior$quad_nodes),
            score_levels = object$input$score_levels, logistic_scale = object$settings$logistic_scale)
        fs <- .drift_compute_fair_scores(fit_core, ref_raters = unique(r_idx), ref_waves = unique(w_idx))
        fair_sel <- if (identical(integration, "posterior"))
            fs$post_expected
        else fs$plugin
        out <- data.frame(Person = object$person_table$Person, NObs = object$person_table$NObs, MeanScore = object$person_table$MeanScore,
            ThetaEAP = object$person_table$ThetaEAP, ThetaPSD = object$person_table$ThetaPSD, FairM = fair_sel,
            FairM_Posterior = fs$post_expected, FairM_PlugIn = fs$plugin, ReferenceTime = paste(t_lvls[unique(w_idx)],
                collapse = "+"), ReferenceRaters = if (length(unique(r_idx)) == length(r_lvls))
                "ALL"
            else paste(r_lvls[unique(r_idx)], collapse = "+"), stringsAsFactors = FALSE)
        return(out)
    }
    type <- if (missing(type))
        "expected"
    else match.arg(type)
    if (identical(type, "fair_score")) {
        stop("Use `newdata = NULL` for person-level `'fair_score'`, or choose `type = 'expected'` / `'probabilities'` with `newdata`.",
            call. = FALSE)
    }
    if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
        stop("`newdata` must be a non-empty data.frame.", call. = FALSE)
    }
    cols <- object$input$columns
    req <- c(cols$rater, cols$time, object$input$active_facets)
    miss_req <- setdiff(req, names(newdata))
    if (length(miss_req) > 0L) {
        stop("`newdata` is missing required column(s): ", paste(miss_req, collapse = ", "), ".", call. = FALSE)
    }
    n_nd <- nrow(newdata)
    r_id <- match(trimws(as.character(newdata[[cols$rater]])), object$input$levels$rater)
    if (anyNA(r_id))
        stop("`newdata` contains unknown rater level(s).", call. = FALSE)
    t_raw <- if (is.factor(newdata[[cols$time]]) || is.numeric(newdata[[cols$time]]) || inherits(newdata[[cols$time]],
        c("Date", "POSIXt"))) {
        as.character(newdata[[cols$time]])
    }
    else {
        trimws(as.character(newdata[[cols$time]]))
    }
    t_id <- match(t_raw, object$input$levels$time)
    if (anyNA(t_id))
        stop("`newdata` contains unknown time/wave level(s).", call. = FALSE)
    eta_offset <- object$matrices$beta_rt[cbind(r_id, t_id)]
    if (length(object$input$active_facets) > 0L) {
        ft <- object$facet_table
        for (fc in object$input$active_facets) {
            sub_ft <- ft[ft$Facet == fc, , drop = FALSE]
            f_id <- match(trimws(as.character(newdata[[fc]])), sub_ft$Level)
            if (anyNA(f_id))
                stop("`newdata` contains unknown level(s) for fixed facet `", fc, "`.", call. = FALSE)
            eta_offset <- eta_offset + sub_ft$Estimate[f_id]
        }
    }
    use_post <- FALSE
    if (!is.null(ability)) {
        if (is.character(ability) && length(ability) == 1L && ability %in% names(newdata)) {
            th_vec <- newdata[[ability]]
        }
        else {
            th_vec <- ability
        }
        if (!is.numeric(th_vec) || any(!is.finite(th_vec)) || !length(th_vec) %in% c(1L, n_nd)) {
            stop("`ability` must be finite numeric values of length 1 or `nrow(newdata)`.", call. = FALSE)
        }
        th_vec <- rep(as.numeric(th_vec), length.out = n_nd)
    }
    else if (cols$person %in% names(newdata)) {
        p_id <- match(trimws(as.character(newdata[[cols$person]])), object$input$levels$person)
        if (anyNA(p_id)) {
            stop("`newdata` contains unknown `Person` ID(s); supply explicit `ability` for new examinees.",
                call. = FALSE)
        }
        th_vec <- object$person_table$ThetaEAP[p_id]
        use_post <- identical(integration, "posterior")
    }
    else {
        stop("Supply either a known `", cols$person, "` column in `newdata` or explicit `ability`.",
            call. = FALSE)
    }
    K <- length(object$input$score_levels)
    s_vals <- as.numeric(object$input$score_levels)
    c_vec <- as.numeric(0:(K - 1L))
    D <- object$settings$logistic_scale
    g_mat <- cbind(0, t(apply(object$matrices$d_rm, 1L, cumsum)))
    if (K == 2L)
        g_mat <- cbind(0, object$matrices$d_rm[, 1L])
    g_nd <- g_mat[r_id, , drop = FALSE]
    prob_at <- function(theta_row) {
        eta <- matrix(0, nrow = n_nd, ncol = K)
        for (k in 2:K) {
            eta[, k] <- D * (c_vec[k] * (theta_row - eta_offset) - g_nd[, k])
        }
        mx <- apply(eta, 1L, max)
        pr <- exp(eta - mx)
        pr/rowSums(pr)
    }
    if (use_post) {
        nodes <- object$posterior$quad_nodes
        W_p <- object$posterior$person_weights[p_id, , drop = FALSE]
        probs <- matrix(0, nrow = n_nd, ncol = K)
        for (q in seq_along(nodes)) {
            probs <- probs + W_p[, q] * prob_at(rep(nodes[q], n_nd))
        }
    }
    else {
        probs <- prob_at(th_vec)
    }
    colnames(probs) <- as.character(object$input$score_levels)
    if (identical(type, "probabilities")) {
        return(probs)
    }
    out_exp <- newdata
    out_exp$ExpectedScore <- as.numeric(probs %*% s_vals)
    out_exp
}

as.data.frame.mfrm_drift <- function(x, row.names = NULL, optional = FALSE, component = c("severity",
    "rater_drift", "steps", "facets", "persons", "summary"), ...) {
    rlang::check_dots_empty()
    component <- match.arg(component)
    switch(component, severity = x$severity_table, rater_drift = x$rater_drift_table, steps = x$step_table,
        facets = x$facet_table, persons = x$person_table, summary = x$model_summary)
}
