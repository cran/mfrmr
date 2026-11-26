# ==============================================================================
# Ultra-Fast Exact Henderson Method 1 (urGENOVA) G-Study & D-Study Engine
# Supports Crossed & Nested Multi-Facet Designs + Multivariate Composite Weights
# Reference: Brennan, R. L. (2001). Generalizability Theory (Ch. 2, 5, 6, 7, 10)
# ==============================================================================

































.mfrm_gstudy_henderson <- function(data, person, facets, score = NULL, scores = NULL, nesting = NULL, weights = NULL,
    d_study_grid = NULL, truncate_negative = TRUE) {
    t0 <- proc.time()[[3L]]
    if (!is.data.frame(data) || nrow(data) == 0L) {
        stop("`data` must be a non-empty data.frame.", call. = FALSE)
    }
    score_cols <- as.character(scores %||% score)
    if (length(score_cols) < 1L || anyNA(score_cols) || any(!nzchar(score_cols))) {
        stop("Supply at least one score column name in `score` or `scores`.", call. = FALSE)
    }
    facets <- unname(as.character(facets))
    if (!is.character(facets) || length(facets) < 1L || anyNA(facets) || any(!nzchar(trimws(facets))) || anyDuplicated(facets) ||
        person %in% facets || any(score_cols %in% c(person, facets))) {
        stop("`facets` must name distinct non-person columns disjoint from `person` and `scores`.", call. = FALSE)
    }
    req <- unique(c(person, facets, score_cols))
    miss <- setdiff(req, names(data))
    if (length(miss) > 0L) {
        stop("Missing column(s) in `data`: ", paste(miss, collapse = ", "), ".", call. = FALSE)
    }
    if (length(facets) > 4L) {
        stop(errorCondition("The G-Study supports up to four non-person facets. Supply a subset with `g_facets` in mfrm_analysis(), or with `facets` in mfrm_g_study().",
            class = "mfrm_gstudy_unavailable"))
    }
    nest_child_idx <- integer(0)
    nest_parent_idx <- integer(0)
    all_factors <- c(person, facets)
    if (!is.null(nesting) && length(nesting) > 0L) {
        if (!is.character(nesting) || is.null(names(nesting)) || any(!nzchar(names(nesting)))) {
            stop("`nesting` must be a named character vector `c(Child = 'Parent')`.", call. = FALSE)
        }
        ch_nms <- names(nesting)
        pa_nms <- unname(as.character(nesting))
        if (any(!ch_nms %in% facets) || any(!pa_nms %in% all_factors) || any(ch_nms == pa_nms)) {
            stop("`nesting` children must be in `facets` and parents in `c(person, facets)`.", call. = FALSE)
        }
        nest_child_idx <- match(ch_nms, all_factors)
        nest_parent_idx <- match(pa_nms, all_factors)
    }
    P_sc <- length(score_cols)
    if (is.null(weights)) {
        w_vec <- rep(1/P_sc, P_sc)
    }
    else {
        if (!is.numeric(weights) || is.complex(weights))
            stop("`weights` must be numeric.", call. = FALSE)
        if (!is.null(names(weights))) {
            if (anyDuplicated(names(weights)) || !setequal(names(weights), score_cols))
                stop("Named `weights` must match `scores` exactly.", call. = FALSE)
            weights <- weights[score_cols]
        }
        w_vec <- as.numeric(weights)
        if (length(w_vec) != P_sc || any(!is.finite(w_vec)) || any(w_vec < 0) || sum(w_vec) <= 0) {
            stop("`weights` must be a nonnegative numeric vector matching `scores` with positive sum.", call. = FALSE)
        }
    }
    names(w_vec) <- score_cols
    Y_raw <- as.matrix(data[, score_cols, drop = FALSE])
    storage.mode(Y_raw) <- "double"
    keep <- rowSums(!is.finite(Y_raw)) == 0L & !is.na(data[[person]])
    for (fc in facets) keep <- keep & !is.na(data[[fc]])
    if (!all(keep)) {
        data <- data[keep, , drop = FALSE]
        Y_raw <- Y_raw[keep, , drop = FALSE]
    }
    N <- nrow(Y_raw)
    if (N < 2L) {
        stop(errorCondition("At least 2 complete observations are required for G-Study.", class = "mfrm_gstudy_unavailable"))
    }
    d <- length(all_factors)
    fac_int <- lapply(all_factors, function(col) as.integer(factor(data[[col]])))
    names(fac_int) <- all_factors
    fac_n_levels <- vapply(fac_int, max, integer(1))
    all_subsets <- list()
    for (sz in seq_len(d)) {
        cmbs <- utils::combn(seq_len(d), sz, simplify = FALSE)
        all_subsets <- c(all_subsets, cmbs)
    }
    is_admissible <- function(v) {
        if (length(nest_child_idx) == 0L)
            return(TRUE)
        for (kk in seq_along(nest_child_idx)) {
            if ((nest_child_idx[kk] %in% v) && !(nest_parent_idx[kk] %in% v))
                return(FALSE)
        }
        TRUE
    }
    subsets <- Filter(is_admissible, all_subsets)
    M <- length(subsets)
    format_subset_label <- function(v) {
        paste(all_factors[v], collapse = ":")
    }
    make_cell_id <- function(idx_vec) {
        if (length(idx_vec) == 1L)
            return(fac_int[[idx_vec]])
        id <- as.double(fac_int[[idx_vec[1L]]] - 1L)
        for (k in idx_vec[-1L]) {
            id <- id * as.double(fac_n_levels[k]) + as.double(fac_int[[k]] - 1L)
        }
        match(id, unique(id))
    }
    cell_ids <- lapply(subsets, make_cell_id)
    cell_counts <- lapply(cell_ids, tabulate)
    has_replicates <- (length(cell_counts[[M]]) < N)
    Y_cen <- sweep(Y_raw, 2L, colMeans(Y_raw), "-")
    y_comp <- as.numeric(Y_cen %*% w_vec)
    k_empty <- vapply(seq_len(M), function(u) {
        cnts <- as.numeric(cell_counts[[u]])
        sum(cnts * cnts)/N
    }, numeric(1))
    if (!has_replicates) {
        K_mat <- matrix(0, nrow = M, ncol = M)
        for (s in seq_len(M)) {
            cid_s <- cell_ids[[s]]
            n_s <- as.numeric(rowsum(rep(1, N), cid_s, reorder = FALSE)[, 1L])
            for (u in seq_len(M)) {
                if (all(subsets[[u]] %in% subsets[[s]])) {
                  k_su <- as.numeric(N)
                }
                else if (s == M) {
                  k_su <- as.numeric(N)
                }
                else {
                  su_idx <- sort(unique(c(subsets[[s]], subsets[[u]])))
                  su_pos <- Position(function(v) identical(v, su_idx), subsets)
                  cid_su <- cell_ids[[su_pos]]
                  first_obs <- match(seq_len(max(cid_su)), cid_su)
                  parent_s <- cid_s[first_obs]
                  cnt_su <- as.numeric(cell_counts[[su_pos]])
                  num_by_s <- rowsum(cnt_su * cnt_su, parent_s, reorder = FALSE)[, 1L]
                  k_su <- sum(num_by_s/n_s)
                }
                K_mat[s, u] <- k_su - k_empty[u]
            }
        }
        comp_names <- vapply(subsets, format_subset_label, character(1))
        comp_names[M] <- paste0(comp_names[M], ",e")
        contains_person <- vapply(subsets, function(v) 1L %in% v, logical(1))
        non_person_subsets <- lapply(subsets, function(v) setdiff(v, 1L) - 1L)
        eff_subsets <- subsets
    }
    else {
        M_ext <- M + 1L
        K_mat <- matrix(0, nrow = M_ext, ncol = M_ext)
        k_empty_ext <- c(k_empty, 1)
        for (s in seq_len(M)) {
            cid_s <- cell_ids[[s]]
            n_s <- as.numeric(rowsum(rep(1, N), cid_s, reorder = FALSE)[, 1L])
            for (u in seq_len(M)) {
                if (all(subsets[[u]] %in% subsets[[s]])) {
                  k_su <- as.numeric(N)
                }
                else {
                  su_idx <- sort(unique(c(subsets[[s]], subsets[[u]])))
                  su_pos <- Position(function(v) identical(v, su_idx), subsets)
                  cid_su <- cell_ids[[su_pos]]
                  first_obs <- match(seq_len(max(cid_su)), cid_su)
                  parent_s <- cid_s[first_obs]
                  cnt_su <- as.numeric(cell_counts[[su_pos]])
                  num_by_s <- rowsum(cnt_su * cnt_su, parent_s, reorder = FALSE)[, 1L]
                  k_su <- sum(num_by_s/n_s)
                }
                K_mat[s, u] <- k_su - k_empty_ext[u]
            }
            K_mat[s, M_ext] <- length(n_s) - 1
        }
        for (u in seq_len(M_ext)) K_mat[M_ext, u] <- as.numeric(N) - k_empty_ext[u]
        comp_names <- c(vapply(subsets, format_subset_label, character(1)), "Residual")
        contains_person <- c(vapply(subsets, function(v) 1L %in% v, logical(1)), TRUE)
        non_person_subsets <- c(lapply(subsets, function(v) setdiff(v, 1L) - 1L), list(seq_along(facets)))
        eff_subsets <- c(subsets, list(seq_len(d + 1L)))
    }
    sv_k <- svd(K_mat)
    keep_k <- sv_k$d > 1e-10 * max(1, sv_k$d[1L])
    if (!all(keep_k)) {
        stop(errorCondition("The G-Study variance components are not identifiable from these assignments. Add shared ratings or revise the facet/nesting design.",
            class = "mfrm_gstudy_unavailable"))
    }
    K_pinv <- solve(K_mat)
    M_tot <- nrow(K_mat)
    RHS_mat <- matrix(0, nrow = M_tot, ncol = P_sc * P_sc)
    tot_sum <- colSums(Y_cen)
    T_empty_mat <- outer(tot_sum, tot_sum)/N
    for (s in seq_len(M)) {
        S_s <- rowsum(Y_cen, cell_ids[[s]], reorder = FALSE)
        cnt_s <- as.numeric(rowsum(rep(1, N), cell_ids[[s]], reorder = FALSE)[, 1L])
        T_s_mat <- crossprod(S_s/sqrt(cnt_s))
        RHS_mat[s, ] <- as.numeric(T_s_mat - T_empty_mat)
    }
    if (has_replicates) {
        T_obs_mat <- crossprod(Y_cen)
        RHS_mat[M_tot, ] <- as.numeric(T_obs_mat - T_empty_mat)
    }
    VC_flat <- K_pinv %*% RHS_mat
    vcov_components <- stats::setNames(vector("list", M_tot), comp_names)
    raw_vc <- numeric(M_tot)
    for (u in seq_len(M_tot)) {
        Sig_u <- matrix(VC_flat[u, ], nrow = P_sc, ncol = P_sc, dimnames = list(score_cols, score_cols))
        Sig_u <- 0.5 * (Sig_u + t(Sig_u))
        vcov_components[[u]] <- Sig_u
        raw_vc[u] <- as.numeric(t(w_vec) %*% Sig_u %*% w_vec)
    }
    est_vc <- if (isTRUE(truncate_negative))
        pmax(0, raw_vc)
    else raw_vc
    tot_vc <- sum(pmax(0, est_vc))
    groups <- c(cell_ids, if (has_replicates) list(seq_len(N)) else list())
    vc_vcov <- .mfrm_gstudy_covariance(groups, pmax(0, raw_vc), K_pinv)
    dimnames(vc_vcov) <- list(comp_names, comp_names)
    vc_se <- sqrt(pmax(0, diag(vc_vcov)))
    vc_df_satt <- ifelse(est_vc > 1e-12 & diag(vc_vcov) > 1e-20, 2 * (est_vc^2)/diag(vc_vcov), NA_real_)
    vc_df <- data.frame(Component = comp_names, Estimate = est_vc, RawEstimate = raw_vc, SE = vc_se, df = vc_df_satt,
        PercentVariance = if (tot_vc > 0)
            100 * pmax(0, est_vc)/tot_vc
        else rep(NA_real_, length(est_vc)), ContainsPerson = contains_person, Truncated = raw_vc < 0 & isTRUE(truncate_negative),
        stringsAsFactors = FALSE)
    obs_n_per_facet <- vapply(seq_along(facets), function(f) {
        f_factor_idx <- f + 1L
        par_idx <- nest_parent_idx[nest_child_idx == f_factor_idx]
        base_cols <- unique(c(1L, par_idx))
        parent_cid <- make_cell_id(base_cols)
        joint_cid <- make_cell_id(unique(c(base_cols, f_factor_idx)))
        ub <- unique(cbind(parent_cid, joint_cid))
        cnt_p <- tabulate(ub[, 1L], nbins = max(parent_cid))
        cnt_p <- cnt_p[cnt_p > 0L]
        length(cnt_p)/sum(1/cnt_p)
    }, numeric(1))
    names(obs_n_per_facet) <- facets
    eval_design_row <- function(n_vec) {
        sigma2_p <- est_vc[1L]
        sigma2_rel <- 0
        sigma2_abs <- 0
        for (j in seq_along(est_vc)[-1L]) {
            np_idx <- non_person_subsets[[j]]
            denom <- if (length(np_idx) == 0L)
                1
            else prod(n_vec[np_idx])
            contrib <- est_vc[j]/denom
            sigma2_abs <- sigma2_abs + contrib
            if (contains_person[j]) {
                sigma2_rel <- sigma2_rel + contrib
            }
        }
        .mfrm_gstudy_metrics(sigma2_p, sigma2_rel, sigma2_abs)
    }
    g_obs_metrics <- eval_design_row(obs_n_per_facet)
    g_summary <- as.data.frame(c(as.list(stats::setNames(obs_n_per_facet, paste0("Harmonic_", facets))), as.list(g_obs_metrics)),
        stringsAsFactors = FALSE)
    if (is.null(d_study_grid)) {
        grids <- lapply(seq_along(facets), function(f) {
            n0 <- max(1L, as.integer(round(obs_n_per_facet[f])))
            sort(unique(c(1L, 2L, 3L, 4L, 6L, 8L, n0, 2L * n0)))
        })
        names(grids) <- facets
        d_study_grid <- do.call(expand.grid, c(grids, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
    }
    else {
        miss_g <- setdiff(facets, names(d_study_grid))
        if (length(miss_g) > 0L) {
            stop("`d_study_grid` must contain column(s): ", paste(miss_g, collapse = ", "), ".", call. = FALSE)
        }
    }
    d_study_grid <- .mfrm_design_grid(d_study_grid, facets)
    d_mat <- do.call(rbind, lapply(seq_len(nrow(d_study_grid)), function(i) eval_design_row(unlist(d_study_grid[i,
        facets, drop = FALSE], use.names = FALSE))))
    d_df <- cbind(d_study_grid[, facets, drop = FALSE], TotalRatingsPerPerson = apply(as.matrix(d_study_grid[,
        facets, drop = FALSE]), 1L, prod), as.data.frame(d_mat, stringsAsFactors = FALSE))
    d_df$MeetsG80 <- d_df$G >= 0.8
    d_df$MeetsPhi80 <- d_df$Phi >= 0.8
    rownames(d_df) <- NULL
    univar_df <- NULL
    if (P_sc > 1L) {
        u_rows <- lapply(seq_len(P_sc), function(p_i) {
            sc_nm <- score_cols[p_i]
            v_p <- vapply(vcov_components, function(S) S[p_i, p_i], numeric(1))
            if (isTRUE(truncate_negative))
                v_p <- pmax(0, v_p)
            sp <- v_p[1L]
            s_rel <- 0
            s_abs <- 0
            for (j in seq_along(v_p)[-1L]) {
                np_idx <- non_person_subsets[[j]]
                den <- if (length(np_idx) == 0L)
                  1
                else prod(obs_n_per_facet[np_idx])
                s_abs <- s_abs + v_p[j]/den
                if (contains_person[j])
                  s_rel <- s_rel + v_p[j]/den
            }
            cbind(data.frame(Score = sc_nm, Weight = unname(w_vec[p_i])), .mfrm_gstudy_metrics(sp, s_rel,
                s_abs))
        })
        univar_df <- do.call(rbind, u_rows)
        rownames(univar_df) <- NULL
    }
    .mfrm_new_gstudy(structure(list(person = person, facets = facets, score = if (P_sc == 1L) score_cols[1L] else score_cols,
        scores = score_cols, weights = w_vec, nesting = nesting, n_obs = N, n_person = fac_n_levels[1L], facet_levels = fac_n_levels[-1L],
        obs_n_per_facet = obs_n_per_facet, non_person_subsets = non_person_subsets, method = "Henderson Method 1 (urGENOVA)",
        truncate_negative = isTRUE(truncate_negative), variance_components = vc_df, vc_vcov = vc_vcov, covariance_basis = "Gaussian plug-in quadratic-form covariance; negative components set to zero for the working covariance",
        vcov_components = if (P_sc > 1L) vcov_components else NULL, univariate = univar_df, g_study_summary = g_summary,
        d_study = d_df, elapsed_sec = proc.time()[[3L]] - t0), class = "mfrm_gtheory"), "urgenova", data)
}


.mfrm_gstudy_vcov_henderson <- function(object, ...) {
  object$vc_vcov
}


.mfrm_gstudy_confint_henderson <- function(object,
                                      parm = NULL,
                                      level = 0.95,
                                      component = c("all", "variance_components", "coefficients"),
                                      ...) {
  if (!is.null(parm) && is.character(parm) && length(parm) == 1L &&
      parm %in% c("all", "variance_components", "coefficients")) {
    component <- parm
    parm <- NULL
  }
  rlang::check_dots_empty()
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one finite number strictly between zero and one.", call. = FALSE)
  }
  component <- match.arg(component)
  alpha <- 1.0 - level
  pct_lbls <- paste0(format(100 * c(alpha / 2, 1 - alpha / 2), trim = TRUE), "%")
  vc <- object$variance_components
  V <- object$vc_vcov

  # 1. Satterthwaite chi-square CIs for positive variance components (0 lower clamp if truncated)
  vc_ci <- matrix(NA_real_, nrow = nrow(vc), ncol = 2L, dimnames = list(vc$Component, pct_lbls))
  vc_status <- rep("unavailable", nrow(vc))
  for (i in seq_len(nrow(vc))) {
    est <- vc$Estimate[i]
    df_i <- vc$df[i]
    if (is.finite(est) && est > 1e-12 && is.finite(df_i) && df_i > 0) {
      # A conservative reporting guard, not a guarantee of coverage above df = 1.
      if (df_i < 1) {
        vc_status[i] <- "unstable_satterthwaite_df"
        next
      }
      q_hi <- stats::qchisq(1 - alpha / 2, df = df_i)
      q_lo <- stats::qchisq(alpha / 2, df = df_i)
      vc_ci[i, 1L] <- (df_i * est) / q_hi
      vc_ci[i, 2L] <- if (q_lo > 0) (df_i * est) / q_lo else Inf
      vc_status[i] <- "approximate_satterthwaite"
    } else if (is.finite(vc$RawEstimate[i]) && is.finite(vc$SE[i])) {
      z <- stats::qnorm(1 - alpha / 2)
      vc_ci[i, ] <- pmax(0.0, vc$RawEstimate[i] + c(-1, 1) * z * vc$SE[i])
      vc_status[i] <- "approximate_truncated_wald"
    }
    if (any(!is.finite(vc_ci[i, ]))) {
      vc_ci[i, ] <- NA_real_
      vc_status[i] <- "unavailable_numerical_range"
    }
  }

  # 2. Logit-transformed delta-method CIs on observed G and Phi
  n_obs_f <- object$obs_n_per_facet
  np_subs <- object$non_person_subsets
  w_rel <- numeric(nrow(vc))
  w_abs <- numeric(nrow(vc))
  for (j in seq_len(nrow(vc))[-1L]) {
    den <- if (length(np_subs[[j]]) == 0L) 1.0 else prod(n_obs_f[np_subs[[j]]])
    w_abs[j] <- if (vc$Estimate[j] > 0) 1.0 / den else 0.0
    if (isTRUE(vc$ContainsPerson[j])) w_rel[j] <- w_abs[j]
  }
  z <- stats::qnorm(1 - alpha / 2)
  logit_ci <- function(sp, s_err, w_err) {
    if (!is.finite(sp) || sp <= 1e-12 || !is.finite(s_err) || s_err <= 1e-12) return(c(NA_real_, NA_real_))
    # g = log(sp) - log(s_err); grad_1 = 1/sp, grad_j = -w_err_j / s_err
    grad <- -w_err / s_err
    grad[1L] <- 1.0 / sp
    v_logit <- as.numeric(t(grad) %*% V %*% grad)
    if (!is.finite(v_logit) || v_logit <= 0) return(c(NA_real_, NA_real_))
    eta <- log(sp / s_err)
    stats::plogis(eta + c(-1, 1) * z * sqrt(v_logit))
  }
  coef_ci <- rbind(
    G = logit_ci(object$g_study_summary$UniverseVar, object$g_study_summary$RelErrorVar, w_rel),
    Phi = logit_ci(object$g_study_summary$UniverseVar, object$g_study_summary$AbsErrorVar, w_abs)
  )
  colnames(coef_ci) <- pct_lbls

  out <- switch(
    component,
    variance_components = vc_ci,
    coefficients = coef_ci,
    all = rbind(vc_ci, coef_ci)
  )
  status <- stats::setNames(c(vc_status,
    ifelse(rowSums(is.finite(coef_ci)) == 2L, "approximate_logit_delta", "unavailable")),
    c(vc$Component, rownames(coef_ci)))
  if (!is.null(parm)) {
    sel <- if (is.character(parm)) match(parm, rownames(out)) else as.integer(parm)
    out <- out[sel[!is.na(sel)], , drop = FALSE]
  }
  attr(out, "status") <- status[rownames(out)]
  out
}


.mfrm_gstudy_print_henderson <- function(x, ...) {
  x <- .mfrm_refresh_gstudy(x)
  cat(sprintf(
    "Generalizability (G-Study) / Decision (D-Study) [%s | N = %d ratings, %d persons | %.4f s]\n",
    x$method, x$n_obs, x$n_person, x$elapsed_sec
  ))
  if (!is.null(x$nesting) && length(x$nesting) > 0L) {
    cat("  Nesting:", paste(paste0(names(x$nesting), " : ", unname(x$nesting)), collapse = ", "), "\n")
  }
  if (length(x$scores) > 1L) {
    cat("  Composite Scores:", paste(sprintf("%s (w=%.2f)", x$scores, x$weights), collapse = ", "), "\n")
  }
  cat("\nVariance Components:\n")
  vc_disp <- x$variance_components[, c("Component", "Estimate", "SE", "RawEstimate", "PercentVariance")]
  vc_disp$Estimate <- round(vc_disp$Estimate, 4)
  vc_disp$SE <- round(vc_disp$SE, 4)
  vc_disp$RawEstimate <- round(vc_disp$RawEstimate, 4)
  vc_disp$PercentVariance <- round(vc_disp$PercentVariance, 2)
  print(vc_disp, row.names = FALSE)
  cat(sprintf(
    "\nObserved Harmonic Design: G = %.4f (RelSEM = %.4f), Phi = %.4f (AbsSEM = %.4f)\n",
    x$g_study_summary$G, x$g_study_summary$RelSEM,
    x$g_study_summary$Phi, x$g_study_summary$AbsSEM
  ))
  cat("  G:", x$g_study_summary$GStatus, "| Phi:", x$g_study_summary$PhiStatus, "\n")
  invisible(x)
}


