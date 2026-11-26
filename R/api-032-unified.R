# ==============================================================================
# Unified 0.3.2 Fast Bias (RajuSA-Welch), Parallel Split DFF, Affine LRT,
# Unified Interrater Agreement, & Adapter Methods
# ==============================================================================

# Extract C++ encoding + parameter state from the standardized calibration
.mfrm_fit_state <- function(fit, data_override = NULL) {
    .mfrm_validate_fit(fit)
    enc <- fit$encoding
    if (!is.null(data_override))
        enc$retained_data <- .mfrm_diagnostic_data(enc, data_override)
    th <- fit$person$Estimate
    excl <- (fit$person$ObsCount <= 0) | (identical(fit$method, "JML") & fit$person$ExtremeFlag != 0L & !is.finite(th))
    beta_list <- lapply(seq_along(enc$facets), function(f) {
        fnm <- enc$facets[f]
        if (fnm %in% enc$base_facets) {
            sub <- fit$facets[fit$facets$Facet == fnm, ]
            sub$Estimate[match(enc$facet_levels[[f]], sub$Level)]
        }
        else {
            sub <- fit$interactions[fit$interactions$Interaction == fnm, ]
            sub$Estimate
        }
    })
    return(list(enc = enc, theta = as.double(th), theta_exclude = as.logical(excl), beta_list = beta_list,
        step_cum = fit$step_cumulative, alpha = fit$alpha_obs))
}

# Fast C++ 2-way & all-pairs bias + FACETS 4.5.3 RajuSA-Welch pairwise contrasts
.mfrm_estimate_bias_cpp <- function(fit,
                                    facet_a = NULL,
                                    facet_b = NULL,
                                    interaction_facets = NULL,
                                    method = c("local", "joint"),
                                    df_method = c("facets", "classical"),
                                    p_adjust = "holm") {
  method <- match.arg(method)
  df_method <- match.arg(df_method)
  st <- .mfrm_fit_state(fit)
  enc <- st$enc
  base_f <- enc$base_facets

  if (identical(interaction_facets, "all") || (is.null(interaction_facets) && is.null(facet_a) && is.null(facet_b))) {
    if (length(base_f) < 2L) {
      stop("At least two non-person facets are required for `interaction_facets = 'all'`.", call. = FALSE)
    }
    cmb <- utils::combn(seq_along(base_f), 2L)
    pair_fa <- as.integer(match(base_f[cmb[1L, ]], enc$facets) - 1L)
    pair_fb <- as.integer(match(base_f[cmb[2L, ]], enc$facets) - 1L)
  } else {
    if (is.null(interaction_facets)) interaction_facets <- c(facet_a, facet_b)
    interaction_facets <- as.character(interaction_facets)
    if (length(interaction_facets) != 2L || any(!interaction_facets %in% base_f)) {
      stop("`interaction_facets` must name two facets in `facets` (or `'all'`).", call. = FALSE)
    }
    pair_fa <- as.integer(match(interaction_facets[1L], enc$facets) - 1L)
    pair_fb <- as.integer(match(interaction_facets[2L], enc$facets) - 1L)
  }

  if (identical(method, "joint")) {
    if (length(pair_fa) != 1L) {
      stop("`method = 'joint'` requires a single 2-way facet pair.", call. = FALSE)
    }
    pair_nms <- enc$facets[c(pair_fa + 1L, pair_fb + 1L)]
    if (!is.null(fit$facets$ShrunkEstimate))
      stop("Joint bias refits require the original, unshrunk calibration.", call. = FALSE)
    m_base <- fit
    enc_int <- .mfrm_add_interaction(enc, pair_nms)
    ctrl <- fit$control
    ctrl$compute_diagnostics <- TRUE
    raw <- mfrm_cpp_fast_fit_spec(enc_int, ctrl)
    m_int <- .mfrm_format_result(enc_int, raw, fit$method, fit$score_method,
                                     ctrl$xtreme, ctrl, NA_real_)
    comparison <- .mfrm_compare_fits(list(m_base, m_int), c("Additive", "Joint"), nested = TRUE)
    readiness <- .mfrm_readiness(m_int)
    if (!.mfrm_readiness(m_base)$ready) {
      readiness$ready <- FALSE
      readiness$status <- "review_required"
      readiness$reason_codes <- unique(c(readiness$reason_codes, "source_calibration_not_ready"))
    }
    tb <- .mfrm_screen_table(m_int$interactions, readiness,
                                  c("ModelSE", "SandwichSE", "z", "pValue"))
    tb$pAdj <- stats::p.adjust(tb$pValue, method = p_adjust)
    lrt_df <- NULL
    if (!is.null(comparison$lrt)) {
      lr <- comparison$lrt[1L, ]
      lrt_df <- data.frame(Interaction = paste(pair_nms, collapse = ":"),
        LogLikAdditive = m_base$loglik, LogLikJoint = m_int$loglik,
        ChiSq = lr$ChiSq, df = lr$df, pValue = lr$pValue)
    }
    return(structure(
      list(
        method = "joint",
        df_method = "wald_normal",
        table = tb,
        pairwise = data.frame(),
        lrt = lrt_df,
        joint_fit = m_int,
        readiness = readiness,
        comparison_basis = comparison$comparison_basis
      ),
      class = c("mfrm_bias", "mfrm_bundle")
    ))
  }

  use_facets_df <- identical(df_method, "facets")
  cpp_out <- mfrm_cpp_fast_bias_pairwise(
    enc, st$theta, st$theta_exclude, st$beta_list, st$step_cum,
    pair_fa, pair_fb, use_facets_df, st$alpha
  )

  cell_dfs <- vector("list", length(pair_fa))
  pair_dfs <- vector("list", length(pair_fa))
  for (pr in seq_along(pair_fa)) {
    fa <- pair_fa[pr] + 1L
    fb <- pair_fb[pr] + 1L
    fa_nm <- enc$facets[fa]
    fb_nm <- enc$facets[fb]
    lvs_a <- enc$facet_levels[[fa]]
    lvs_b <- enc$facet_levels[[fb]]

    cl <- cpp_out$cells[[pr]]
    if (length(cl$level_a) > 0L) {
      p_cell <- 2.0 * stats::pt(-abs(cl$t_stat), df = pmax(1.0, cl$df))
      cell_dfs[[pr]] <- data.frame(
        Pair = paste0(fa_nm, ":", fb_nm),
        FacetA = fa_nm,
        LevelA = lvs_a[cl$level_a + 1L],
        FacetB = fb_nm,
        LevelB = lvs_b[cl$level_b + 1L],
        ObsCount = cl$obs_count,
        ObsScore = cl$obs_score,
        ExpScore = cl$exp_score,
        ObsExpAvg = (cl$obs_score - cl$exp_score) / pmax(1e-12, cl$obs_count),
        BaseMeasure = cl$base_measure,
        TargetMeasure = cl$target_measure,
        BiasSize = cl$bias_est,
        BiasSE = cl$bias_se,
        t = cl$t_stat,
        df = cl$df,
        pValue = p_cell,
        pAdj = stats::p.adjust(p_cell, method = p_adjust),
        stringsAsFactors = FALSE
      )
    }

    ct <- cpp_out$contrasts[[pr]]
    if (length(ct$level_a) > 0L) {
      p_con <- ifelse(
        is.finite(ct$t_stat) & is.finite(ct$df) & ct$df > 0,
        2.0 * stats::pt(-abs(ct$t_stat), df = ct$df),
        NA_real_
      )
      pair_dfs[[pr]] <- data.frame(
        Pair = paste0(fa_nm, ":", fb_nm),
        TargetFacet = fa_nm,
        TargetLevel = lvs_a[ct$level_a + 1L],
        ContextFacet = fb_nm,
        Context1 = lvs_b[ct$level_b1 + 1L],
        Context2 = lvs_b[ct$level_b2 + 1L],
        Measure1 = ct$target_1,
        SE1 = ct$se_1,
        Obs1 = ct$n_1,
        Measure2 = ct$target_2,
        SE2 = ct$se_2,
        Obs2 = ct$n_2,
        Contrast = ct$contrast,
        ContrastSE = ct$contrast_se,
        t = ct$t_stat,
        df = ct$df,
        dfUnrounded = ct$df_unrounded,
        TestType = if (use_facets_df) "RajuSA-Welch (FACETS 4.5.3)" else "Welch-Satterthwaite",
        pValue = p_con,
        pAdj = stats::p.adjust(p_con, method = p_adjust),
        stringsAsFactors = FALSE
      )
    }
  }

  cells_all <- do.call(rbind, cell_dfs)
  if (is.null(cells_all)) cells_all <- data.frame()
  rownames(cells_all) <- NULL

  pairs_all <- do.call(rbind, pair_dfs)
  if (is.null(pairs_all)) pairs_all <- data.frame()
  rownames(pairs_all) <- NULL

  readiness <- .mfrm_readiness(fit)
  cells_all <- .mfrm_screen_table(cells_all, readiness, c("BiasSE", "t", "pValue", "pAdj"))
  pairs_all <- .mfrm_screen_table(pairs_all, readiness, c("SE1", "SE2", "ContrastSE", "t", "pValue", "pAdj"))

  sum_df <- if (nrow(cells_all) > 0L) {
    do.call(rbind, lapply(split(cells_all, cells_all$Pair), function(d) {
      data.frame(
        Pair = d$Pair[1L],
        Cells = nrow(d),
        MaxAbsBias = max(abs(d$BiasSize), na.rm = TRUE),
        RMSBias = sqrt(mean(d$BiasSize^2, na.rm = TRUE)),
        FlaggedCells = if (all(is.na(d$pAdj))) NA_integer_ else sum(d$pAdj < 0.05, na.rm = TRUE),
        UnavailableCells = sum(is.na(d$pAdj)),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame()
  }
  rownames(sum_df) <- NULL

  structure(
    list(
      method = "local",
      df_method = if (use_facets_df) "facets" else "classical",
      test_label = if (use_facets_df) "RajuSA-Welch (FACETS 4.5.3 / von Davier 2025)" else "Classical Welch-Satterthwaite",
      table = cells_all,
      pairwise = pairs_all,
      summary = sum_df,
      readiness = readiness
    ),
    class = c("mfrm_bias", "mfrm_bundle")
  )
}

# Fast C++ DFF analyzer: parallel warm-started MML-LRT, CML-LRT, RajuSA-Welch, MH, & Logistic
.mfrm_analyze_dff <- function(fit,
                                   facet,
                                   group,
                                   data = NULL,
                                   focal = NULL,
                                   method = c("mml_lrt", "cml_lrt", "raju_welch", "mh", "logistic"),
                                   df_method = c("facets", "classical"),
                                   min_obs = 5L,
                                   p_adjust = "holm",
                                   threads = 1L) {
  method <- match.arg(method)
  df_method <- match.arg(df_method)
  st <- .mfrm_fit_state(fit, data_override = data)
  enc0 <- st$enc
  orig_df <- enc0$retained_data
  orig_df$.mfrm_dff_row <- seq_len(nrow(orig_df))
  if (!group %in% names(orig_df)) {
    stop("`group` column '", group, "' not found in `data`.", call. = FALSE)
  }
  if (!facet %in% enc0$base_facets) {
    stop("`facet` '", facet, "' is not one of the model facets: ",
         paste(enc0$base_facets, collapse = ", "), ".", call. = FALSE)
  }

  source_readiness <- .mfrm_readiness(fit)
  finish <- function(result) {
    if (!method %in% c("mml_lrt", "cml_lrt")) {
      result$table <- .mfrm_screen_table(result$table, source_readiness,
        c("SE", "SandwichSE", "WaldZ", "Statistic", "NonUniformStat", "NonUniformP", "pValue", "pAdj"))
    }
    result$table$Flagged <- ifelse(is.finite(result$table$pAdj), result$table$pAdj < 0.05, NA)
    result$source_readiness <- source_readiness
    result
  }
  if (method %in% c("mml_lrt", "cml_lrt")) {
    if (!inherits(fit, "mfrm_fit") || !fit$model %in% c("RSM", "PCM") ||
        any(enc0$facet_dummy != 0L) || any(enc0$facet_noncenter != 0L) ||
        any(enc0$has_custom_facet_map != 0L) || isTRUE(enc0$has_custom_step_map) ||
        any(vapply(enc0$anchor_fixed, function(x) any(x != 0L), logical(1))) ||
        !is.null(enc0$category_support) || !is.null(fit$facets$ShrunkEstimate) ||
        any(c(fit$control$facet_l2_pen, fit$control$slope_l2_pen) > 0)) {
      stop("Split DFF likelihood tests require an unpenalized unified RSM/PCM fit without anchors, custom facet/step constraints, interactions, or restricted category support. These settings cannot be dropped during refitting.", call. = FALSE)
    }
  }

  g_chr <- trimws(as.character(orig_df[[group]]))
  keep <- !is.na(g_chr) & nzchar(g_chr) & is.finite(as.numeric(orig_df[[enc0$score_col]]))
  d_work <- orig_df[keep, , drop = FALSE]
  g_chr <- g_chr[keep]
  g_lvls <- sort(unique(g_chr))
  if (length(g_lvls) < 2L) {
    stop("Grouping column '", group, "' must have at least 2 distinct levels.", call. = FALSE)
  }
  foc_lvl <- if (is.null(focal)) g_lvls[2L] else as.character(focal)[1L]
  if (!foc_lvl %in% g_lvls) stop("`focal` ('", foc_lvl, "') not found in `group`.", call. = FALSE)
  ref_lvl <- setdiff(g_lvls, foc_lvl)[1L]

  # Filter to Reference vs Focal comparison if >2 groups
  if (length(g_lvls) > 2L) {
    sub_2 <- g_chr %in% c(ref_lvl, foc_lvl)
    d_work <- d_work[sub_2, , drop = FALSE]
    g_chr <- g_chr[sub_2]
  }
  d_work$.dff_group <- factor(g_chr, levels = c(ref_lvl, foc_lvl))

  f_vals <- as.character(d_work[[facet]])
  all_elems <- .mfrm_sort_levels(f_vals)
  tab_eg <- table(factor(f_vals, levels = all_elems), d_work$.dff_group)
  cand_elems <- all_elems[tab_eg[, ref_lvl] >= min_obs & tab_eg[, foc_lvl] >= min_obs]
  if (length(cand_elems) == 0L) {
    stop("No elements of facet '", facet, "' have at least `min_obs` = ", min_obs,
         " ratings in both reference ('", ref_lvl, "') and focal ('", foc_lvl, "') groups.", call. = FALSE)
  }

  # 1. RajuSA-Welch (FACETS 4.5.3 pairwise fixed-calibration DFF)
  if (identical(method, "raju_welch")) {
    if (".dff_grp_dummy" %in% enc0$facets) stop("Reserved DFF context name conflicts with a model facet.", call. = FALSE)
    d_aug <- d_work
    d_aug$.dff_grp_dummy <- as.character(d_work$.dff_group)
    f_facets <- c(enc0$base_facets, ".dff_grp_dummy")
    # Add a zero-valued context to the fixed calibration, without refitting it.
    bfit <- fit
    eb <- enc0
    ii <- d_work$.mfrm_dff_row
    for (nm in c("person0", "score0", "weight", "step_group0")) eb[[nm]] <- eb[[nm]][ii]
    eb$facet_mat0 <- eb$facet_mat0[ii, , drop = FALSE]
    eb$retained_data <- d_aug
    eb$n_obs <- nrow(d_aug)
    eb$facets <- c(eb$facets, ".dff_grp_dummy")
    eb$base_facets <- c(eb$base_facets, ".dff_grp_dummy")
    eb$facet_levels[[".dff_grp_dummy"]] <- c(ref_lvl, foc_lvl)
    eb$facet_n_levels <- c(eb$facet_n_levels, 2L)
    eb$facet_signs <- c(eb$facet_signs, -1)
    eb$facet_dummy <- c(eb$facet_dummy, 1L)
    eb$facet_noncenter <- c(eb$facet_noncenter, 0L)
    eb$has_custom_facet_map <- c(eb$has_custom_facet_map, 0L)
    eb$custom_facet_q <- c(eb$custom_facet_q, 0L)
    for (nm in c("anchor_fixed", "custom_facet_anchored")) eb[[nm]][[length(eb$facets)]] <- integer(2L)
    for (nm in c("anchor_val", "custom_facet_d0")) eb[[nm]][[length(eb$facets)]] <- numeric(2L)
    eb$custom_facet_B[[length(eb$facets)]] <- matrix(0, 2L, 0L)
    eb$facet_mat0 <- cbind(eb$facet_mat0, as.integer(d_work$.dff_group) - 1L)
    dummy <- fit$facets[rep(1L, 2L), , drop = FALSE]
    dummy$Facet <- ".dff_grp_dummy"; dummy$Level <- c(ref_lvl, foc_lvl)
    dummy$Estimate <- 0; dummy$Sign <- -1
    bfit$facets <- rbind(fit$facets, dummy)
    bfit$encoding <- eb
    bfit$alpha_obs <- fit$alpha_obs[ii]
    b_res <- .mfrm_estimate_bias_cpp(
      bfit, facet_a = facet, facet_b = ".dff_grp_dummy",
      method = "local", df_method = df_method, p_adjust = p_adjust
    )
    pw <- b_res$pairwise
    pw <- pw[pw$TargetLevel %in% cand_elems, , drop = FALSE]
    # Orient contrast as Focal - Reference
    flip <- pw$Context1 == ref_lvl & pw$Context2 == foc_lvl
    est_ref <- ifelse(flip, pw$Measure1, pw$Measure2)
    est_foc <- ifelse(flip, pw$Measure2, pw$Measure1)
    se_ref <- ifelse(flip, pw$SE1, pw$SE2)
    se_foc <- ifelse(flip, pw$SE2, pw$SE1)
    n_ref <- ifelse(flip, pw$Obs1, pw$Obs2)
    n_foc <- ifelse(flip, pw$Obs2, pw$Obs1)
    delta <- est_foc - est_ref
    t_val <- delta / pw$ContrastSE
    p_val <- ifelse(is.finite(t_val) & is.finite(pw$df) & pw$df > 0,
                    2.0 * stats::pt(-abs(t_val), df = pw$df), NA_real_)
    p_adj <- stats::p.adjust(p_val, method = p_adjust)
    tbl <- data.frame(
      Facet = facet,
      Level = pw$TargetLevel,
      ReferenceGroup = ref_lvl,
      FocalGroup = foc_lvl,
      ObsRef = n_ref,
      ObsFocal = n_foc,
      EstRef = est_ref,
      EstFocal = est_foc,
      Contrast = delta,
      SE = pw$ContrastSE,
      Statistic = t_val,
      df = pw$df,
      dfUnrounded = pw$dfUnrounded,
      pValue = p_val,
      pAdj = p_adj,
      Flagged = !is.na(p_adj) & (p_adj < 0.05),
      Method = b_res$test_label,
      stringsAsFactors = FALSE
    )
    return(finish(structure(
      list(method = "raju_welch", facet = facet, group = group,
           reference = ref_lvl, focal = foc_lvl, table = tbl),
      class = "mfrm_dff"
    )))
  }

  # 2. Parallel Warm-Started C++ Split-Element MML-LRT & CML-LRT
  if (method %in% c("mml_lrt", "cml_lrt")) {
    est_method <- if (identical(method, "mml_lrt")) "MML" else "CML"
    m_code <- if (identical(method, "mml_lrt")) 2L else 3L

    # Check if group is nested within person (person-level impact covariate)
    pg_tab <- unique(d_work[, c(enc0$person_col, ".dff_group"), drop = FALSE])
    person_nested_group <- !anyDuplicated(pg_tab[[enc0$person_col]])
    pop_fml <- if (identical(est_method, "MML") && person_nested_group) ~ .dff_group else if (identical(est_method, "MML")) ~ 1 else NULL
    pdata_g <- if (person_nested_group && !is.null(pop_fml)) pg_tab else NULL

    base_enc <- .mfrm_encode_dataset(
      data = d_work, person = enc0$person_col, facets = enc0$base_facets,
      score = enc0$score_col, model = enc0$model, step_facet = enc0$step_facet,
      rating_min = enc0$rating_min, rating_max = enc0$rating_max,
      weight = enc0$weight_col, positive_facets = enc0$base_facets[enc0$facet_signs > 0],
      population_formula = pop_fml, person_data = pdata_g, person_id = enc0$person_col
    )
    pop_design <- function(e) {
      if (!identical(est_method, "MML")) return(e)
      old <- if (enc0$pop_mode == 1L) enc0$pop_X_mat[match(e$person_levels, enc0$person_levels), , drop = FALSE] else
        matrix(1, e$n_person, 1L, dimnames = list(NULL, "(Intercept)"))
      if (qr(old)$rank < ncol(old)) stop("The source population design is rank deficient in the selected DFF groups.", call. = FALSE)
      if (person_nested_group) {
        gv <- as.numeric(pdata_g$.dff_group[match(e$person_levels, pdata_g[[enc0$person_col]])] == foc_lvl)
        if (qr(cbind(old, gv))$rank > ncol(old)) {
          old <- cbind(old, gv)
          colnames(old)[ncol(old)] <- paste0("DFFGroup:", foc_lvl)
        }
      }
      e$pop_X_mat <- old; e$pop_col_names <- colnames(old)
      e$pop_mode <- 1L
      e
    }
    base_enc <- pop_design(base_enc)
    base_ctrl <- fit$control
    base_ctrl$method_code <- m_code
    base_ctrl$compute_diagnostics <- TRUE
    base_ctrl$compute_wle <- TRUE
    base_ctrl$wle_primary <- FALSE
    base_ctrl$bias_facet_a <- base_ctrl$bias_facet_b <- -1L
    refit_readiness <- function(e, rw) {
      f <- .mfrm_format_result(e, rw, est_method, "auto", base_ctrl$xtreme, base_ctrl, NA_real_)
      r <- .mfrm_readiness(f)
      if (identical(est_method, "MML")) {
        issue <- tryCatch({.mfrm_check_mml_population_boundary(e); NULL}, error = conditionMessage)
        if (!is.null(issue)) {
          r$ready <- FALSE; r$status <- "review_required"
          r$reason_codes <- c(r$reason_codes, paste0("population_boundary_review: ", issue))
        }
      }
      r
    }
    base_raw <- mfrm_cpp_fast_fit_spec(base_enc, base_ctrl)
    base_readiness <- refit_readiness(base_enc, base_raw)
    ll_base <- as.numeric(base_raw$loglik)

    # Prepare split-element specs for all candidate elements in parallel!
    f_pos <- match(facet, enc0$base_facets)
    spec_list <- vector("list", length(cand_elems))
    ctrl_list <- vector("list", length(cand_elems))

    for (k in seq_along(cand_elems)) {
      el <- cand_elems[k]
      d_sp <- d_work
      is_el <- as.character(d_sp[[facet]]) == el
      d_sp[[facet]] <- as.character(d_sp[[facet]])
      d_sp[[facet]][is_el & d_sp$.dff_group == ref_lvl] <- paste0(el, "__REF__")
      d_sp[[facet]][is_el & d_sp$.dff_group == foc_lvl] <- paste0(el, "__FOC__")
      # If split facet is also PCM step_facet, keep shared step thresholds across the split pair
      sc_decl <- NULL
      if (identical(enc0$model, "PCM") && identical(enc0$step_facet, facet)) {
        K_st <- enc0$n_cat - 1L
        sc_decl <- data.frame(
          Level = rep(c(paste0(el, "__REF__"), paste0(el, "__FOC__")), each = K_st),
          From = rep(enc0$rating_min + seq(0L, K_st - 1L), times = 2L),
          To = rep(enc0$rating_min + seq_len(K_st), times = 2L),
          Group = rep(paste0("sh_", seq_len(K_st)), times = 2L),
          stringsAsFactors = FALSE
        )
      }
      enc_sp <- .mfrm_encode_dataset(
        data = d_sp, person = enc0$person_col, facets = enc0$base_facets,
        score = enc0$score_col, model = enc0$model, step_facet = enc0$step_facet,
        rating_min = enc0$rating_min, rating_max = enc0$rating_max,
        weight = enc0$weight_col, positive_facets = enc0$base_facets[enc0$facet_signs > 0],
        population_formula = pop_fml, person_data = pdata_g, person_id = enc0$person_col,
        step_constraints = sc_decl
      )
      spec_list[[k]] <- pop_design(enc_sp)
      ctrl_list[[k]] <- base_ctrl
    }

    sp_raws <- mfrm_cpp_fast_batch_warm_list(spec_list, ctrl_list, .mfrm_threads(threads))

    rows_out <- lapply(seq_along(cand_elems), function(k) {
      el <- cand_elems[k]
      enc_sp <- spec_list[[k]]
      rw <- sp_raws[[k]]
      lvs_sp <- enc_sp$facet_levels[[f_pos]]
      i_ref <- match(paste0(el, "__REF__"), lvs_sp)
      i_foc <- match(paste0(el, "__FOC__"), lvs_sp)
      est_sp <- rw$facet_est[[f_pos]]
      cm_sp <- rw$facet_cov_model[[f_pos]]
      cs_sp <- rw$facet_cov_sandwich[[f_pos]]

      delta <- est_sp[i_foc] - est_sp[i_ref]
      v_mod <- cm_sp[i_foc, i_foc] + cm_sp[i_ref, i_ref] - 2.0 * cm_sp[i_foc, i_ref]
      se_mod <- if (is.finite(v_mod) && v_mod > 0) sqrt(v_mod) else NA_real_
      se_sw <- NA_real_
      if (identical(est_method, "MML") && nrow(cs_sp) >= max(i_ref, i_foc)) {
        v_sw <- cs_sp[i_foc, i_foc] + cs_sp[i_ref, i_ref] - 2.0 * cs_sp[i_foc, i_ref]
        if (is.finite(v_sw) && v_sw > 0) se_sw <- sqrt(v_sw)
      }
      readiness <- refit_readiness(enc_sp, rw)
      readiness$reason_codes <- unique(c(if (length(base_readiness$reason_codes)) paste0("base:", base_readiness$reason_codes) else character(), readiness$reason_codes))
      lr_raw <- 2 * (as.numeric(rw$loglik) - ll_base)
      df_k <- as.integer(rw$n_free_par) - as.integer(base_raw$n_free_par)
      if (!is.finite(lr_raw) || lr_raw < -1e-7 * max(1, abs(ll_base)) || df_k < 1L)
        readiness$reason_codes <- c(readiness$reason_codes, "invalid_nested_likelihood")
      readiness$ready <- isTRUE(base_readiness$ready) && isTRUE(readiness$ready) && !length(readiness$reason_codes)
      readiness$status <- if (readiness$ready) "ready" else "review_required"
      lr_stat <- if (readiness$ready) max(0, lr_raw) else NA_real_
      pval <- if (readiness$ready) stats::pchisq(lr_stat, df = df_k, lower.tail = FALSE) else NA_real_

      .mfrm_screen_table(data.frame(
        Facet = facet,
        Level = el,
        ReferenceGroup = ref_lvl,
        FocalGroup = foc_lvl,
        ObsRef = as.numeric(tab_eg[el, ref_lvl]),
        ObsFocal = as.numeric(tab_eg[el, foc_lvl]),
        EstRef = est_sp[i_ref],
        EstFocal = est_sp[i_foc],
        Contrast = delta,
        SE = se_mod,
        SandwichSE = se_sw,
        WaldZ = ifelse(is.finite(se_mod) & se_mod > 0, delta / se_mod, NA_real_),
        Statistic = lr_stat,
        df = df_k,
        pValue = pval,
        stringsAsFactors = FALSE
      ), readiness, c("SE", "SandwichSE", "WaldZ", "Statistic", "pValue"))
    })
    tbl <- do.call(rbind, rows_out)
    tbl$pAdj <- stats::p.adjust(tbl$pValue, method = p_adjust)
    tbl$Flagged <- !is.na(tbl$pAdj) & (tbl$pAdj < 0.05)
    tbl$Method <- toupper(method)
    rownames(tbl) <- NULL

    return(finish(structure(
      list(method = method, facet = facet, group = group,
           reference = ref_lvl, focal = foc_lvl, loglik_base = ll_base, base_readiness = base_readiness, table = tbl),
      class = "mfrm_dff"
    )))
  }

  # 3. Stratified Generalized Mantel-Haenszel (MH / Liu-Agresti mean-score statistic)
  p_map <- stats::setNames(st$theta, enc0$person_levels)
  th_obs <- p_map[as.character(d_work[[enc0$person_col]])]
  y_obs <- as.numeric(d_work[[enc0$score_col]])

  if (identical(method, "mh")) {
    # Stratify on person total raw score (or deciles of theta if unbalanced)
    p_tot <- stats::ave(y_obs, d_work[[enc0$person_col]], FUN = sum)
    strata <- if (length(unique(p_tot)) <= 25L) as.factor(p_tot) else cut(th_obs, breaks = unique(stats::quantile(th_obs, probs = seq(0, 1, 0.2), na.rm = TRUE)), include.lowest = TRUE)

    rows_mh <- lapply(cand_elems, function(el) {
      idx_e <- which(as.character(d_work[[facet]]) == el)
      y_e <- y_obs[idx_e]
      g_e <- as.integer(d_work$.dff_group[idx_e] == foc_lvl)
      s_e <- droplevels(as.factor(strata[idx_e]))
      num <- 0.0
      den <- 0.0
      w_diff <- 0.0
      w_sum <- 0.0
      for (st_lv in levels(s_e)) {
        ii <- which(s_e == st_lv)
        n_k <- length(ii)
        n1 <- sum(g_e[ii] == 1L)
        n0 <- n_k - n1
        if (n1 == 0L || n0 == 0L || n_k < 2L) next
        yk <- y_e[ii]
        mu_k <- mean(yk)
        v_k <- sum((yk - mu_k)^2) / (n_k - 1.0)
        if (v_k <= 1e-12) next
        s_foc <- sum(yk[g_e[ii] == 1L])
        num <- num + (s_foc - n1 * mu_k)
        den <- den + (as.double(n1) * as.double(n0) / n_k) * v_k
        wt <- (as.double(n1) * as.double(n0)) / n_k
        w_diff <- w_diff + wt * (mean(yk[g_e[ii] == 1L]) - mean(yk[g_e[ii] == 0L]))
        w_sum <- w_sum + wt
      }
      stat <- if (den > 1e-12) (num * num) / den else NA_real_
      pval <- if (is.finite(stat)) stats::pchisq(stat, df = 1L, lower.tail = FALSE) else NA_real_
      # Positive Contrast = harder/stricter against Focal (lower score in Focal)
      smd <- if (w_sum > 0) -(w_diff / w_sum) else NA_real_
      data.frame(
        Facet = facet, Level = el, ReferenceGroup = ref_lvl, FocalGroup = foc_lvl,
        ObsRef = as.numeric(tab_eg[el, ref_lvl]), ObsFocal = as.numeric(tab_eg[el, foc_lvl]),
        Contrast = smd, SE = if (den > 1e-12 && w_sum > 0) sqrt(den) / w_sum else NA_real_,
        Statistic = stat, df = 1, pValue = pval, stringsAsFactors = FALSE
      )
    })
    tbl <- do.call(rbind, rows_mh)
    tbl$pAdj <- stats::p.adjust(tbl$pValue, method = p_adjust)
    tbl$Flagged <- !is.na(tbl$pAdj) & (tbl$pAdj < 0.05)
    tbl$Method <- "Mantel-Haenszel (Mean-Score)"
    rownames(tbl) <- NULL
    return(finish(structure(
      list(method = "mh", facet = facet, group = group,
           reference = ref_lvl, focal = foc_lvl, table = tbl),
      class = "mfrm_dff"
    )))
  }

  # 4. Logistic / Cumulative-Logit Uniform & Non-Uniform DFF
  rows_log <- lapply(cand_elems, function(el) {
    idx_e <- which(as.character(d_work[[facet]]) == el & is.finite(th_obs))
    ye <- y_obs[idx_e]
    ge <- as.numeric(d_work$.dff_group[idx_e] == foc_lvl)
    te <- th_obs[idx_e]
    # Fit cumulative logit via IRLS / optim on ordered thresholds
    u_cats <- sort(unique(ye))
    if (length(u_cats) < 2L) {
      return(data.frame(
        Facet = facet, Level = el, ReferenceGroup = ref_lvl, FocalGroup = foc_lvl,
        ObsRef = as.numeric(tab_eg[el, ref_lvl]), ObsFocal = as.numeric(tab_eg[el, foc_lvl]),
        Contrast = NA_real_, SE = NA_real_, Statistic = NA_real_, df = 1,
        NonUniformStat = NA_real_, NonUniformP = NA_real_, pValue = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    y_int <- match(ye, u_cats)
    Kc <- length(u_cats) - 1L
    fit_clm <- function(X_mat) {
      q <- ncol(X_mat)
      p0 <- c(seq(-1, 1, length.out = Kc), rep(0, q))
      if (Kc > 1L) p0[2:Kc] <- log(diff(p0[1:Kc]))
      nll <- function(par) {
        cuts <- if (Kc == 1L) par[1L] else cumsum(c(par[1L], exp(par[2:Kc])))
        xb <- as.numeric(X_mat %*% par[(Kc + 1L):(Kc + q)])
        c_hi <- ifelse(y_int == Kc + 1L, Inf, cuts[pmin(Kc, y_int)] - xb)
        c_lo <- ifelse(y_int == 1L, -Inf, cuts[pmax(1L, y_int - 1L)] - xb)
        pr <- stats::plogis(c_hi) - stats::plogis(c_lo)
        if (any(!is.finite(pr) | pr <= 1e-15)) return(1e10)
        -sum(log(pr))
      }
      stats::optim(p0, nll, method = "BFGS", hessian = TRUE)
    }
    m0 <- fit_clm(cbind(te))
    m1 <- fit_clm(cbind(te, ge))
    m2 <- fit_clm(cbind(te, ge, te * ge))
    lr_u <- max(0.0, 2.0 * (m0$value - m1$value))
    lr_nu <- max(0.0, 2.0 * (m1$value - m2$value))
    beta_g <- m1$par[Kc + 2L]
    H_inv <- tryCatch(solve(m1$hessian), error = function(e) matrix(NA_real_, nrow(m1$hessian), ncol(m1$hessian)))
    se_g <- if (is.finite(H_inv[Kc + 2L, Kc + 2L]) && H_inv[Kc + 2L, Kc + 2L] > 0) sqrt(H_inv[Kc + 2L, Kc + 2L]) else NA_real_
    data.frame(
      Facet = facet, Level = el, ReferenceGroup = ref_lvl, FocalGroup = foc_lvl,
      ObsRef = as.numeric(tab_eg[el, ref_lvl]), ObsFocal = as.numeric(tab_eg[el, foc_lvl]),
      Contrast = -beta_g, SE = se_g,
      Statistic = lr_u, df = 1,
      NonUniformStat = lr_nu,
      NonUniformP = stats::pchisq(lr_nu, df = 1L, lower.tail = FALSE),
      pValue = stats::pchisq(lr_u, df = 1L, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, rows_log)
  tbl$pAdj <- stats::p.adjust(tbl$pValue, method = p_adjust)
  tbl$Flagged <- !is.na(tbl$pAdj) & (tbl$pAdj < 0.05)
  tbl$Method <- "Cumulative Logit LRT"
  rownames(tbl) <- NULL
  finish(structure(
    list(method = "logistic", facet = facet, group = group,
         reference = ref_lvl, focal = foc_lvl, table = tbl),
    class = "mfrm_dff"
  ))
}

# ==============================================================================
# Affine-Subspace Model Comparison & `anova()` for `mfrm_fit` & `mfrm_fit`
# ==============================================================================

.mfrm_full_affine_space <- function(fit, ref_step_owners = NULL) {
  enc <- fit$encoding
  F_base <- length(enc$base_facets)
  K <- enc$n_cat - 1L

  # Facet linear space [sum(L_f) x q_facets]
  f_d0 <- numeric(0)
  f_B_blocks <- list()
  for (f in seq_len(F_base)) {
    Lf <- enc$facet_n_levels[f]
    if (enc$facet_dummy[f] == 1L) {
      d0 <- numeric(Lf)
      B <- matrix(0.0, Lf, 0L)
    } else if (enc$has_custom_facet_map[f] == 1L) {
      d0 <- enc$custom_facet_d0[[f]]
      B <- enc$custom_facet_B[[f]]
    } else {
      d0 <- numeric(Lf)
      anc <- which(enc$anchor_fixed[[f]] == 1L)
      unanc <- which(enc$anchor_fixed[[f]] == 0L)
      if (length(anc) > 0L) d0[anc] <- enc$anchor_val[[f]][anc]
      targ <- if (enc$centering_scope == 1L) -sum(d0[anc]) else 0.0
      if (length(unanc) == 0L) {
        B <- matrix(0.0, Lf, 0L)
      } else if (enc$facet_noncenter[f] == 1L) {
        B <- matrix(0.0, Lf, length(unanc))
        B[cbind(unanc, seq_along(unanc))] <- 1.0
      } else {
        qf <- length(unanc) - 1L
        d0[unanc[length(unanc)]] <- targ
        B <- matrix(0.0, Lf, qf)
        if (qf > 0L) {
          B[cbind(unanc[seq_len(qf)], seq_len(qf))] <- 1.0
          B[unanc[length(unanc)], ] <- -1.0
        }
      }
    }
    f_d0 <- c(f_d0, d0)
    f_B_blocks[[f]] <- B
  }

  # Step linear space lifted to ref_step_owners if comparing RSM vs PCM
  G <- enc$n_step_groups
  if (isTRUE(enc$has_custom_step_map)) {
    s_d0 <- replace(enc$custom_step_d0, !is.finite(enc$custom_step_d0), 0.0)
    s_B <- enc$custom_step_B
  } else {
    qs <- fit$xi_layout$step_len
    s_d0 <- numeric(G * K)
    s_B <- matrix(0.0, G * K, qs)
    cc <- 0L
    for (g in seq_len(G)) {
      ag <- enc$step_cat_min[g]
      bg <- enc$step_cat_max[g]
      Kg <- bg - ag
      if (Kg >= 2L) {
        for (t in seq_len(Kg - 1L)) {
          s_B[(g - 1L) * K + ag + t, cc + t] <- 1.0
          s_B[(g - 1L) * K + bg, cc + t] <- -1.0
        }
        cc <- cc + (Kg - 1L)
      }
    }
  }
  if (identical(enc$model, "RSM") && !is.null(ref_step_owners) && length(ref_step_owners) > 1L) {
    lift <- kronecker(rep(1.0, length(ref_step_owners)), diag(1.0, K))
    s_d0 <- as.numeric(lift %*% s_d0)
    s_B <- lift %*% s_B
  }

  # Log-slope linear space [sum(L_f) x q_slopes] + log_sigma [1 x q_sig]
  sl_d0 <- numeric(0)
  sl_B_blocks <- list()
  sl_facets <- enc$slope_facets %||% character(0)
  for (f in seq_len(F_base)) {
    Lf <- enc$facet_n_levels[f]
    fnm <- enc$base_facets[f]
    if (fnm %in% sl_facets) {
      si <- match(fnm, sl_facets)
      sl_d0 <- c(sl_d0, enc$slope_d0[[si]])
      sl_B_blocks[[f]] <- enc$slope_B[[si]]
    } else {
      sl_d0 <- c(sl_d0, numeric(Lf))
      sl_B_blocks[[f]] <- matrix(0.0, Lf, 0L)
    }
  }
  sig_d0 <- if (enc$pop_mode == 1L) 0.0 else log(enc$fixed_sigma)
  sig_B <- if (enc$pop_mode == 1L) matrix(1.0, 1L, 1L) else matrix(0.0, 1L, 0L)

  # Combine block-diagonal basis
  all_blocks <- c(f_B_blocks, list(s_B), sl_B_blocks, list(sig_B))
  nr_tot <- sum(vapply(all_blocks, nrow, integer(1)))
  nc_tot <- sum(vapply(all_blocks, ncol, integer(1)))
  B_tot <- matrix(0.0, nr_tot, nc_tot)
  r_cur <- 0L
  c_cur <- 0L
  for (blk in all_blocks) {
    nr <- nrow(blk)
    nc <- ncol(blk)
    if (nr > 0L && nc > 0L) {
      B_tot[r_cur + seq_len(nr), c_cur + seq_len(nc)] <- blk
    }
    r_cur <- r_cur + nr
    c_cur <- c_cur + nc
  }
  offset_tot <- c(f_d0, s_d0, sl_d0, sig_d0)
  if (ncol(B_tot) == 0L || all(abs(B_tot) == 0)) {
    Q <- matrix(0.0, length(offset_tot), 0L)
  } else {
    sv <- svd(B_tot)
    keep <- sv$d > 1e-10 * max(1.0, sv$d[1L])
    Q <- sv$u[, keep, drop = FALSE]
  }
  list(
    offset = offset_tot,
    Q = Q,
    rank = ncol(Q),
    interactions = names(enc$interaction_meta)
  )
}

.mfrm_compare_fits <- function(fits, labels = NULL, nested = FALSE, warn_constraints = TRUE) {
  if (length(fits) < 2L) stop("Supply at least two fitted models to compare.", call. = FALSE)
  invisible(lapply(fits, .mfrm_validate_fit))
  if (!is.logical(warn_constraints) || length(warn_constraints) != 1L || is.na(warn_constraints))
    stop("warn_constraints must be TRUE or FALSE.", call. = FALSE)
  if (is.null(labels)) labels <- paste0("Model_", seq_along(fits))
  if (!is.character(labels) || length(labels) != length(fits) || anyNA(labels))
    stop("Supply one label per compared model.", call. = FALSE)

  tbl <- data.frame(
    Model = labels,
    Family = vapply(fits, `[[`, character(1), "model"),
    Method = vapply(fits, `[[`, character(1), "method"),
    N = vapply(fits, `[[`, integer(1), "n_obs"),
    FreeParams = vapply(fits, `[[`, integer(1), "n_free_par"),
    LogLik = vapply(fits, `[[`, numeric(1), "loglik"),
    AIC = vapply(fits, `[[`, numeric(1), "aic"),
    BIC = vapply(fits, `[[`, numeric(1), "bic"),
    Converged = vapply(fits, `[[`, logical(1), "converged"),
    stringsAsFactors = FALSE
  )
  tbl$DeltaAIC <- tbl$AIC - min(tbl$AIC, na.rm = TRUE)
  tbl$DeltaBIC <- tbl$BIC - min(tbl$BIC, na.rm = TRUE)

  lrt_df <- NULL
  lrt_status <- "not_requested"
  nesting_rel <- "unreviewed"

  if (isTRUE(nested)) {
    methods <- unique(tbl$Method)
    ns <- unique(tbl$N)
    same_data <- all(vapply(fits[-1L], function(f) .mfrm_same_observations(fits[[1L]], f), logical(1)))
    ready <- all(vapply(fits, function(f) .mfrm_readiness(f)$ready, logical(1)))
    penalized <- any(vapply(fits, function(f) any(c(f$control$facet_l2_pen, f$control$slope_l2_pen) > 0) ||
                             !is.null(f$facets$ShrunkEstimate), logical(1)))
    if (length(methods) != 1L || !methods %in% c("MML", "CML") || length(ns) != 1L || !same_data || !ready || penalized) {
      if (warn_constraints) warning("Nested LRT requires the same unpenalized MML or CML likelihood, identical observations and weights, and numerically eligible converged fits.", call. = FALSE)
      lrt_status <- "ineligible"
    } else {
      # Order models by FreeParams ascending
      ord <- order(tbl$FreeParams)
      ref_owners <- NULL
      for (f in fits) {
        if (f$model %in% c("PCM", "GPCM")) {
          ref_owners <- f$encoding$step_owners
          break
        }
      }
      spaces <- lapply(fits[ord], .mfrm_full_affine_space, ref_step_owners = ref_owners)
      lbl_ord <- labels[ord]
      fp_ord <- tbl$FreeParams[ord]
      ll_ord <- tbl$LogLik[ord]

      lrt_rows <- list()
      all_ok <- TRUE
      for (k in seq_len(length(ord) - 1L)) {
        s1 <- spaces[[k]]
        s2 <- spaces[[k + 1L]]
        in_12 <- mfrm_affine_subspace_contained(s1, s2) && all(s1$interactions %in% s2$interactions)
        in_21 <- mfrm_affine_subspace_contained(s2, s1) && all(s2$interactions %in% s1$interactions)
        df_diff <- fp_ord[k + 1L] - fp_ord[k]
        if (in_12 && in_21 && df_diff == 0L) {
          nesting_rel <- "equivalent"
          all_ok <- FALSE
          if (warn_constraints) warning("Models span equivalent parameter subspaces (0 restrictions).", call. = FALSE)
          break
        }
        if (!in_12 || df_diff < 1L) {
          nesting_rel <- "not_nested"
          all_ok <- FALSE
          if (warn_constraints) warning("Affine parameter subspaces of compared models are not nested.", call. = FALSE)
          break
        }
        lr <- 2 * (ll_ord[k + 1L] - ll_ord[k])
        if (!is.finite(lr) || lr < -1e-7 * max(1, abs(ll_ord[k]))) {
          all_ok <- FALSE
          if (warn_constraints) warning("The larger model has a lower maximized likelihood; LRT unavailable. Review optimization and integration.", call. = FALSE)
          break
        }
        chisq <- max(0, lr)
        lrt_rows[[k]] <- data.frame(
          Simpler = lbl_ord[k],
          Complex = lbl_ord[k + 1L],
          FreeParamsSimpler = fp_ord[k],
          FreeParamsComplex = fp_ord[k + 1L],
          df = df_diff,
          ChiSq = chisq,
          pValue = stats::pchisq(chisq, df = df_diff, lower.tail = FALSE),
          stringsAsFactors = FALSE
        )
      }
      if (all_ok && length(lrt_rows) > 0L) {
        lrt_df <- do.call(rbind, lrt_rows)
        lrt_status <- "computed"
        nesting_rel <- "nested"
      }
    }
  }

  structure(
    list(
      table = tbl,
      lrt = lrt_df,
      comparison_basis = list(
        lrt_status = lrt_status,
        nesting_review = list(relation = nesting_rel, eligible = identical(lrt_status, "computed"))
      )
    ),
    class = "mfrm_comparison"
  )
}


anova.mfrm_fit <- function(object, ..., nested = TRUE) {
 .mfrm_validate_fit(object)
 fits <- c(list(object),list(...))
 expressions <- c(list(substitute(object)),as.list(substitute(list(...)))[-1L])
 .mfrm_compare_fits(fits, labels=vapply(expressions,deparse1,character(1)),nested=nested)
}



# ==============================================================================
# Unified `interrater_agreement()` (Exact, Kappa, Gwet's AC2, Krippendorff's alpha, ICC)
# ==============================================================================

























interrater_agreement <- function(data, person = "Person", rater = "Rater", score = "Score", context_facets = NULL,
    unit = NULL, categories = NULL, metric = c("all", "exact", "kappa", "gwet_ac2", "krippendorff_alpha",
        "icc"), weights = c("quadratic", "linear", "unweighted"), model_expected = FALSE) {
    if (isTRUE(model_expected) && !inherits(data, "mfrm_fit"))
        stop("model_expected requires a fitted unified-engine calibration.", call. = FALSE)
    metric <- match.arg(metric)
    weights <- match.arg(weights)
    if (!is.null(unit) && length(unit) >= 1L) {
        person <- unit[1L]
        context_facets <- if (length(unit) >= 2L)
            unit[-1L]
        else NULL
    }
    if (inherits(data, "mfrm_fit")) {
        enc <- data$encoding
        person <- enc$person_col
        score <- enc$score_col
        if (!rater %in% enc$base_facets)
            rater <- enc$base_facets[1L]
        if (is.null(context_facets))
            context_facets <- setdiff(enc$base_facets, rater)
        df <- enc$retained_data
    }
    else if (is.data.frame(data)) {
        df <- data
    }
    else {
        stop("`data` must be a rating data.frame or a fitted `mfrm_fit` object.", call. = FALSE)
    }
    req <- c(person, rater, score, context_facets)
    miss <- setdiff(req, names(df))
    if (length(miss) > 0L) {
        stop("Missing column(s) in `data`: ", paste(miss, collapse = ", "), ".", call. = FALSE)
    }
    y <- as.numeric(df[[score]])
    keep <- is.finite(y) & !is.na(df[[person]]) & !is.na(df[[rater]])
    for (cf in context_facets) keep <- keep & !is.na(df[[cf]])
    df <- df[keep, , drop = FALSE]
    y <- y[keep]
    unit_cols <- c(person, context_facets)
    unit_id <- if (length(unit_cols) == 1L) {
        as.character(df[[unit_cols]])
    }
    else {
        as.character(.mfrm_mvgt_group(lapply(df[unit_cols], factor)))
    }
    r_chr <- as.character(df[[rater]])
    u_lvls <- unique(unit_id)
    r_lvls <- .mfrm_sort_levels(r_chr)
    cats <- if (!is.null(categories)) {
        sort(unique(c(as.numeric(categories), y)))
    }
    else {
        sort(unique(y))
    }
    cats <- cats[is.finite(cats)]
    Q <- length(cats)
    if (anyDuplicated(df[c(unit_cols, rater)])) {
        stop("Repeated ratings for the same unit and rater are ambiguous. Include occasion in `unit`/`context_facets` or resolve duplicates explicitly.",
            call. = FALSE)
    }
    mat <- matrix(NA_real_, nrow = length(u_lvls), ncol = length(r_lvls), dimnames = list(u_lvls, r_lvls))
    u_idx <- match(unit_id, u_lvls)
    r_idx <- match(r_chr, r_lvls)
    for (i in seq_along(y)) {
        mat[u_idx[i], r_idx[i]] <- y[i]
    }
    n_per_u <- rowSums(!is.na(mat))
    mat_multi <- mat[n_per_u >= 2L, , drop = FALSE]
    if (nrow(mat_multi) == 0L) {
        stop("No rating units were evaluated by 2 or more raters.", call. = FALSE)
    }
    W <- matrix(0, Q, Q)
    span <- max(1, max(cats) - min(cats))
    for (a in seq_len(Q)) {
        for (b in seq_len(Q)) {
            if (identical(weights, "unweighted")) {
                W[a, b] <- if (a == b)
                  1
                else 0
            }
            else if (identical(weights, "linear")) {
                W[a, b] <- 1 - abs(cats[a] - cats[b])/span
            }
            else {
                W[a, b] <- 1 - ((cats[a] - cats[b])/span)^2
            }
        }
    }
    pw_rows <- list()
    if (ncol(mat_multi) >= 2L) {
        cmbs <- utils::combn(ncol(mat_multi), 2L)
        for (cc in seq_len(ncol(cmbs))) {
            r1 <- cmbs[1L, cc]
            r2 <- cmbs[2L, cc]
            ok <- !is.na(mat_multi[, r1]) & !is.na(mat_multi[, r2])
            if (sum(ok) < 2L)
                next
            x1 <- mat_multi[ok, r1]
            x2 <- mat_multi[ok, r2]
            i1 <- match(x1, cats)
            i2 <- match(x2, cats)
            ex <- mean(x1 == x2)
            adj <- mean(abs(x1 - x2) <= 1)
            po <- mean(W[cbind(i1, i2)])
            p1 <- tabulate(i1, nbins = Q)/length(i1)
            p2 <- tabulate(i2, nbins = Q)/length(i2)
            pe <- sum(outer(p1, p2) * W)
            kap <- if (1 - pe > 1e-12)
                (po - pe)/(1 - pe)
            else NA_real_
            pw_rows[[length(pw_rows) + 1L]] <- data.frame(Rater1 = r_lvls[r1], Rater2 = r_lvls[r2], N = sum(ok),
                ExactAgreement = ex, AdjacentAgreement = adj, Kappa = kap, stringsAsFactors = FALSE)
        }
    }
    pw_df <- if (length(pw_rows) > 0L)
        do.call(rbind, pw_rows)
    else data.frame()
    tot_pairs <- 0
    ex_sum <- 0
    adj_sum <- 0
    po_w_sum <- 0
    cat_freq <- numeric(Q)
    for (u in seq_len(nrow(mat_multi))) {
        vals <- mat_multi[u, !is.na(mat_multi[u, ])]
        m <- length(vals)
        iv <- match(vals, cats)
        for (k in iv) cat_freq[k] <- cat_freq[k] + 1/m
        for (a in seq_len(m - 1L)) {
            for (b in seq.int(a + 1L, m)) {
                tot_pairs <- tot_pairs + 1
                if (vals[a] == vals[b])
                  ex_sum <- ex_sum + 1
                if (abs(vals[a] - vals[b]) <= 1)
                  adj_sum <- adj_sum + 1
                po_w_sum <- po_w_sum + W[iv[a], iv[b]]
            }
        }
    }
    exact_ov <- ex_sum/tot_pairs
    adj_ov <- adj_sum/tot_pairs
    po_ov <- po_w_sum/tot_pairs
    pi_q <- cat_freq/nrow(mat_multi)
    tw_sum <- sum(W)
    pe_gwet <- if (Q > 1L)
        (tw_sum/(Q * (Q - 1))) * sum(pi_q * (1 - pi_q))
    else 0
    gwet_ac2 <- if (1 - pe_gwet > 1e-12)
        (po_ov - pe_gwet)/(1 - pe_gwet)
    else NA_real_
    O_mat <- matrix(0, Q, Q)
    for (u in seq_len(nrow(mat_multi))) {
        vals <- mat_multi[u, !is.na(mat_multi[u, ])]
        m <- length(vals)
        iv <- match(vals, cats)
        for (a in seq_len(m)) {
            for (b in seq_len(m)) {
                if (a != b)
                  O_mat[iv[a], iv[b]] <- O_mat[iv[a], iv[b]] + 1/(m - 1)
            }
        }
    }
    n_tot_c <- sum(O_mat)
    n_c <- rowSums(O_mat)
    D_mat <- 1 - W
    Do <- sum(O_mat * D_mat)/n_tot_c
    De <- sum(outer(n_c, n_c) * D_mat)/(n_tot_c * (n_tot_c - 1))
    krip_alpha <- if (De > 1e-12)
        1 - (Do/De)
    else NA_real_
    d_ur <- data.frame(Unit = unit_id, Rater = r_chr, Score = y, stringsAsFactors = FALSE)
    gt_ur <- tryCatch(.mfrm_gstudy_henderson(d_ur, person = "Unit", facets = "Rater", score = "Score", truncate_negative = TRUE),
        error = function(e) NULL)
    icc_2_1 <- NA_real_
    icc_2_k <- NA_real_
    if (!is.null(gt_ur)) {
        vc <- gt_ur$variance_components$Estimate
        vp <- vc[1L]
        vr <- vc[2L]
        ve <- vc[3L]
        if (vp + vr + ve > 0)
            icc_2_1 <- vp/(vp + vr + ve)
        icc_2_k <- gt_ur$g_study_summary$Phi
    }
    mean_kap <- if (nrow(pw_df) > 0L)
        mean(pw_df$Kappa, na.rm = TRUE)
    else NA_real_
    sum_tbl <- data.frame(Metric = c("ExactAgreement", "AdjacentAgreement", "MeanWeightedKappa", "GwetAC2",
        "KrippendorffAlpha", "ICC_2_1", "ICC_2_k"), Estimate = c(exact_ov, adj_ov, mean_kap, gwet_ac2, krip_alpha,
        icc_2_1, icc_2_k), Weights = c("unweighted", "+/-1 category", weights, weights, weights, "variance",
        "variance"), Units = nrow(mat_multi), Raters = ncol(mat_multi), stringsAsFactors = FALSE)
    if (!identical(metric, "all")) {
        keep_m <- switch(metric, exact = c("ExactAgreement", "AdjacentAgreement"), kappa = "MeanWeightedKappa",
            gwet_ac2 = "GwetAC2", krippendorff_alpha = "KrippendorffAlpha", icc = c("ICC_2_1", "ICC_2_k"))
        sum_tbl <- sum_tbl[sum_tbl$Metric %in% keep_m, , drop = FALSE]
    }
    model_table <- NULL
    if (isTRUE(model_expected)) {
        ready <- .mfrm_readiness(data)
        eligible <- if (ready$ready) {
            mfrm_cpp_expected_agreement(split(seq_along(y), u_idx), y, predict(data, type = "probs")[keep,
                , drop = FALSE])
        }
        else rep(0, 6L)
        used <- eligible[2L]
        values <- if (used > 0)
            eligible[3L:6L]/used
        else rep(NA_real_, 4L)
        kappa <- if (used > 0 && values[3L] < 1 - 1e-12)
            (values[1L] - values[3L])/(1 - values[3L])
        else NA_real_
        model_table <- data.frame(ObservedExact = values[1L], ObservedWithinOne = values[2L], ExpectedExact = values[3L],
            ExpectedWithinOne = values[4L], ModelBasedKappa = kappa, Units = eligible[1L], Pairs = used, ObservedPairs = tot_pairs,
            Status = if (!ready$ready)
                paste(ready$reason_codes, collapse = "; ")
            else if (used == 0)
                "no_usable_model_pairs"
            else if (used < tot_pairs)
                "partially_available"
            else "available")
    }
    structure(list(summary = sum_tbl, model_expected = model_table, pairwise = pw_df, table = pw_df, settings = list(person = person,
        rater = rater, score = score, context_facets = context_facets, categories = cats, metric = metric,
        weights = weights)), class = "mfrm_agreement")
}


print.mfrm_agreement <- function(x, ...) {
  cat("Inter-Rater Agreement Summary:\n")
  disp <- x$summary
  disp$Estimate <- round(disp$Estimate, 4)
  print(disp, row.names = FALSE)
  invisible(x)
}

# ==============================================================================
# Adapters for `diagnose_mfrm`, `fair_average_table`, `compute_information`,
# `make_anchor_table`, `anchor_to_baseline`, `detect_anchor_drift`, & `mfrm_d_study`
# on `mfrm_fit` / `mfrm_gtheory`
# ==============================================================================

.mfrm_diagnose_fit <- function(fit, residual_threshold = 2.0,
                                   residual_pca = "none", pca_max_factors = 10L) {
  p_df <- fit$person
  f_df <- fit$facets

  calc_rel <- function(est, se, eap = FALSE) {
    ok <- is.finite(est) & is.finite(se) & (se > 0)
    if (sum(ok) < 2L) {
      return(list(obs_var = NA_real_, mse = NA_real_, true_var = NA_real_,
                  rel = NA_real_, sep = NA_real_, strata = NA_real_, chisq = NA_real_, df = 0L, p = NA_real_))
    }
    e <- est[ok]
    s <- se[ok]
    v_obs <- stats::var(e)
    mse <- mean(s * s)
    v_true <- if (eap) v_obs else max(0.0, v_obs - mse)
    total <- if (eap) v_obs + mse else v_obs
    rel <- if (total > 0) v_true / total else NA_real_
    sep <- if (mse > 0) sqrt(v_true / mse) else 0.0
    strata <- (4.0 * sep + 1.0) / 3.0
    w <- 1.0 / (s * s)
    mu_w <- sum(w * e) / sum(w)
    chisq <- sum(w * (e - mu_w)^2)
    df <- length(e) - 1L
    list(obs_var = v_obs, mse = mse, true_var = v_true, rel = rel, sep = sep,
         strata = strata, chisq = chisq, df = df, p = stats::pchisq(chisq, df = df, lower.tail = FALSE))
  }

  rel_rows <- list()
  eap <- identical(fit$score_method, "EAP")
  pr_m <- calc_rel(p_df$Estimate, p_df$ModelSE, eap)
  pr_r <- calc_rel(p_df$Estimate, p_df$RealSE, eap)
  rel_rows[[1L]] <- data.frame(
    Facet = "Person", N = sum(is.finite(p_df$Estimate)),
    ReliabilityBasis = if (eap) "EAP posterior variance decomposition" else "Error-corrected measure variance",
    ModelReliability = pr_m$rel, RealReliability = pr_r$rel,
    ModelSeparation = pr_m$sep, RealSeparation = pr_r$sep,
    Strata = pr_r$strata, HomogeneityChiSq = pr_m$chisq, df = pr_m$df, pValue = pr_m$p,
    stringsAsFactors = FALSE
  )
  for (fnm in unique(f_df$Facet)) {
    sub <- f_df[f_df$Facet == fnm, ]
    fr_m <- calc_rel(sub$Estimate, sub$ModelSE)
    fr_r <- calc_rel(sub$Estimate, sub$RealSE)
    rel_rows[[length(rel_rows) + 1L]] <- data.frame(
      Facet = fnm, N = nrow(sub),
      ReliabilityBasis = "Error-corrected measure variance",
      ModelReliability = fr_m$rel, RealReliability = fr_r$rel,
      ModelSeparation = fr_m$sep, RealSeparation = fr_r$sep,
      Strata = fr_r$strata, HomogeneityChiSq = fr_m$chisq, df = fr_m$df, pValue = fr_m$p,
      stringsAsFactors = FALSE
    )
  }
  readiness <- .mfrm_readiness(fit)
  rel_tbl <- .mfrm_screen_table(do.call(rbind, rel_rows), readiness, c("HomogeneityChiSq", "pValue"))

  p_meas <- data.frame(
    Facet = "Person", Level = p_df$Person, Estimate = p_df$Estimate,
    ModelSE = p_df$ModelSE, RealSE = p_df$RealSE, FairM = p_df$FairM,
    ObsCount = p_df$ObsCount, ObsAvg = p_df$ObsAvg, ExpAvg = p_df$ExpAvg,
    Infit = p_df$InfitMnSq, InfitZSTD = p_df$InfitZSTD,
    Outfit = p_df$OutfitMnSq, OutfitZSTD = p_df$OutfitZSTD,
    stringsAsFactors = FALSE
  )
  f_meas <- data.frame(
    Facet = f_df$Facet, Level = f_df$Level, Estimate = f_df$Estimate,
    ModelSE = f_df$ModelSE, RealSE = f_df$RealSE, FairM = f_df$FairM,
    ObsCount = f_df$ObsCount, ObsAvg = f_df$ObsAvg, ExpAvg = f_df$ExpAvg,
    Infit = f_df$InfitMnSq, InfitZSTD = f_df$InfitZSTD,
    Outfit = f_df$OutfitMnSq, OutfitZSTD = f_df$OutfitZSTD,
    stringsAsFactors = FALSE
  )
  measures_tbl <- rbind(p_meas, f_meas)
  rownames(measures_tbl) <- NULL

  std_res <- stats::residuals(fit, type = "pearson")
  exp_sc <- stats::fitted(fit, type = "response")
  enc <- fit$encoding
  unexp_idx <- which(is.finite(std_res) & abs(std_res) >= residual_threshold)
  unexp_df <- enc$retained_data[unexp_idx, c(enc$person_col, enc$base_facets, enc$score_col), drop = FALSE]
  if (nrow(unexp_df) > 0L) {
    unexp_df$Expected <- exp_sc[unexp_idx]
    unexp_df$StdResidual <- std_res[unexp_idx]
    unexp_df <- unexp_df[order(-abs(unexp_df$StdResidual)), , drop = FALSE]
  }

  pca_overall <- pca_by_facet <- NULL
  if (residual_pca != "none") {
    cols <- enc$base_facets
    if (any(c("Person", "StdResidual") %in% cols)) {
      stop("Person and StdResidual are reserved facet names for residual PCA.", call. = FALSE)
    }
    obs <- data.frame(Person = enc$person_levels[enc$person0 + 1L],
                      enc$retained_data[cols], StdResidual = std_res, check.names = FALSE)
    if (residual_pca %in% c("overall", "both")) pca_overall <- compute_pca_overall(obs, cols, pca_max_factors)
    if (residual_pca %in% c("facet", "both")) pca_by_facet <- compute_pca_by_facet(obs, cols, pca_max_factors)
  }
  structure(
    list(
      readiness = readiness,
      measures = measures_tbl,
      reliability = rel_tbl,
      unexpected = unexp_df,
      steps = fit$steps,
      slopes = fit$slopes,
      connectivity = fit$connectivity,
      residual_pca_mode = residual_pca,
      residual_pca_overall = pca_overall,
      residual_pca_by_facet = pca_by_facet,
      fit = fit
    ),
    class = "mfrm_diagnostics"
  )
}


print.mfrm_diagnostics <- function(x, ...) {
  x$readiness <- .mfrm_readiness(x$fit)
  x$reliability <- .mfrm_screen_table(x$reliability, x$readiness, c("HomogeneityChiSq", "pValue"))
  cat("Many-Facet Diagnostic Summary:\n")
  cat("  Numerical inference:", x$readiness$status, paste(x$readiness$reason_codes, collapse = "; "), "\n")
  r_disp <- x$reliability
  num_c <- vapply(r_disp, is.numeric, logical(1))
  r_disp[num_c] <- lapply(r_disp[num_c], round, 4)
  print(r_disp, row.names = FALSE)
  cat(sprintf("\nFlagged Unexpected Responses (|z| >= 2): %d\n", nrow(x$unexpected)))
  invisible(x)
}


summary.mfrm_diagnostics <- function(object, digits = 3L, top_n = 10L, include_person = FALSE, ...) {
  readiness <- .mfrm_readiness(object$fit)
  object$reliability <- .mfrm_screen_table(object$reliability, readiness, c("HomogeneityChiSq", "pValue"))
  meas <- object$measures
  if (!isTRUE(include_person)) {
    meas_sub <- meas[meas$Facet != "Person", , drop = FALSE]
  } else {
    meas_sub <- meas
  }
  flagged_fit <- meas_sub[
    (!is.na(meas_sub$Infit) & (meas_sub$Infit < 0.5 | meas_sub$Infit > 1.5)) |
    (!is.na(meas_sub$Outfit) & (meas_sub$Outfit < 0.5 | meas_sub$Outfit > 1.5)),
    , drop = FALSE
  ]
  top_unexp <- if (nrow(object$unexpected) > 0L) utils::head(object$unexpected, max(1L, as.integer(top_n))) else object$unexpected
  structure(
    list(
      readiness = readiness,
      reliability = object$reliability,
      flagged_elements = flagged_fit,
      unexpected = top_unexp,
      n_unexpected = nrow(object$unexpected),
      steps = object$steps,
      slopes = object$slopes,
      connectivity = object$connectivity,
      residual_pca_overall = object$residual_pca_overall,
      residual_pca_by_facet = object$residual_pca_by_facet,
      digits = as.integer(digits)
    ),
    class = "summary.mfrm_diagnostics"
  )
}


print.summary.mfrm_diagnostics <- function(x, ...) {
  d <- x$digits %||% 3L
  cat("Many-Facet Diagnostics Summary\n")
  cat("  Numerical inference:", x$readiness$status, paste(x$readiness$reason_codes, collapse = "; "), "\n")
  cat(sprintf("  Connected: %s (%d subset%s) | Unexpected responses: %d | Misfit-flagged elements: %d\n\n",
              if (isTRUE(x$connectivity$connected)) "YES" else "NO",
              x$connectivity$n_subsets %||% 1L,
              if (identical(x$connectivity$n_subsets, 1L)) "" else "s",
              x$n_unexpected, nrow(x$flagged_elements)))
  r_disp <- x$reliability
  num_c <- vapply(r_disp, is.numeric, logical(1))
  r_disp[num_c] <- lapply(r_disp[num_c], round, d)
  print(r_disp, row.names = FALSE)
  if (nrow(x$flagged_elements) > 0L) {
    cat("\nFlagged Elements (MnSq outside [0.5, 1.5]):\n")
    fl <- x$flagged_elements
    nc <- vapply(fl, is.numeric, logical(1))
    fl[nc] <- lapply(fl[nc], round, d)
    print(fl, row.names = FALSE)
  }
  invisible(x)
}

.mfrm_information_fit <- function(fit, theta_range, theta_points, rating_plan, exposure_basis) {
  check <- summary(fit)$readiness
  if (!check$ready) stop("Calibration is not ready for information: ",
                         paste(check$reason_codes, collapse = ", "), ".", call. = FALSE)
  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      any(!is.finite(theta_range)) || theta_range[1L] >= theta_range[2L] ||
      !is.numeric(theta_points) || length(theta_points) != 1L || !is.finite(theta_points) ||
      theta_points < 2L || theta_points != floor(theta_points) || theta_points > .Machine$integer.max) {
    stop("Supply increasing finite theta bounds and an integer theta_points of at least two.", call. = FALSE)
  }
  th_grid <- seq(theta_range[1L], theta_range[2L], length.out = as.integer(theta_points))
  enc <- fit$encoding
  C <- enc$n_cat
  sform <- enc$slope_formula %||% 0L
  facets <- enc$base_facets
  used <- enc$weight > 0
  observed <- enc$retained_data[used, facets, drop = FALSE]
  if (!is.null(rating_plan)) {
    if (identical(exposure_basis, "per_person_average")) {
      stop("A rating_plan already specifies one person's rating counts.", call. = FALSE)
    }
    if (any(c("Count", "ObservedContext") %in% facets) ||
        !is.data.frame(rating_plan) || !nrow(rating_plan) || anyDuplicated(names(rating_plan)) ||
        !setequal(names(rating_plan), c(facets, "Count"))) {
      stop("rating_plan must contain every non-Person facet and Count.", call. = FALSE)
    }
    contexts <- rating_plan[facets]
    counts <- rating_plan$Count
    if (!is.numeric(counts) || !is.null(dim(counts)) ||
        any(!is.finite(counts) | counts < 1 | counts != floor(counts)) || !is.finite(sum(counts))) {
      stop("Rating-plan Count must be positive finite integers.", call. = FALSE)
    }
    for (f in facets) {
      values <- contexts[[f]]
      levels <- enc$facet_levels[[match(f, facets)]]
      if (!is.atomic(values) || !is.null(dim(values)) || anyNA(values) ||
          !all(as.character(values) %in% levels)) {
        stop("Unknown or missing rating-plan level in facet ", f, ".", call. = FALSE)
      }
    }
    if (anyDuplicated(contexts)) stop("Combine repeated rating-plan contexts in Count.", call. = FALSE)
    exposure_basis <- "rating_plan"
  } else if (identical(exposure_basis, "unique_contexts")) {
    contexts <- unique(observed)
    counts <- rep(1, nrow(contexts))
  } else {
    grouped <- stats::aggregate(list(Count = enc$weight[used]), by = observed, FUN = sum)
    contexts <- grouped[seq_along(facets)]
    counts <- grouped[[length(facets) + 1L]]
    if (identical(exposure_basis, "per_person_average")) {
      counts <- counts / length(unique(enc$person0[used]))
    }
  }
  observed_context <- rep(TRUE, nrow(contexts))
  if (!is.null(rating_plan)) {
    observed_context <- vapply(seq_len(nrow(contexts)), function(i) {
      any(Reduce(`&`, lapply(facets, function(f) as.character(observed[[f]]) == as.character(contexts[[f]][i]))))
    }, logical(1))
  }
  indexes <- lapply(facets, function(f) match(as.character(contexts[[f]]), enc$facet_levels[[f]]))
  names(indexes) <- facets
  comb_eta <- numeric(nrow(contexts))
  for (f in seq_along(facets)) {
    name <- facets[f]
    tb <- fit$facets[fit$facets$Facet == name, , drop = FALSE]
    values <- tb$Estimate[match(enc$facet_levels[[name]], tb$Level)]
    comb_eta <- comb_eta + enc$facet_signs[f] * values[indexes[[name]]]
  }
  for (im in enc$interaction_meta) {
    i <- (indexes[[im$facet_a]] - 1L) * length(enc$facet_levels[[im$facet_b]]) + indexes[[im$facet_b]]
    tb <- fit$interactions[fit$interactions$Interaction == im$name, , drop = FALSE]
    comb_eta <- comb_eta + enc$facet_signs[im$facet_idx] * tb$Estimate[i]
  }
  if (any(!is.finite(comb_eta))) stop("A rating context has unavailable facet locations.", call. = FALSE)
  comb_sg <- if (enc$model_code == 1L) indexes[[enc$step_facet]] - 1L else integer(nrow(contexts))
  comb_alpha <- rep(1, nrow(contexts))
  for (f in enc$slope_facets) {
    sl <- fit$slopes[fit$slopes$Facet == f, , drop = FALSE]
    comb_alpha <- comb_alpha * sl$Estimate[match(as.character(contexts[[f]]), sl$Level)]
  }
  if (any(!is.finite(comb_alpha) | comb_alpha <= 0)) {
    stop("Information requires finite positive discrimination in every context.", call. = FALSE)
  }

  test_info <- numeric(length(th_grid))
  exp_total <- numeric(length(th_grid))
  for (t_i in seq_along(th_grid)) {
    theta <- th_grid[t_i]
    ev <- .mfrm_eval_row_predictions(theta + comb_eta, fit$step_cumulative, comb_sg,
      enc$step_cat_min, enc$step_cat_max, enc$rating_min,
      alpha_vec = if (length(enc$slope_facets)) comb_alpha else NULL,
      theta_vec = rep(theta, nrow(contexts)), base_eta_vec = comb_eta, slope_formula = sform)
    if (any(!is.finite(ev$variance) | !is.finite(ev$expected))) {
      stop("Category probabilities are unavailable in an information context.", call. = FALSE)
    }
    test_info[t_i] <- sum(counts * comb_alpha^2 * ev$variance)
    exp_total[t_i] <- sum(counts * ev$expected)
  }
  if (any(!is.finite(test_info) | !is.finite(exp_total))) {
    stop("The combined rating exposures exceed the numerical range for information.", call. = FALSE)
  }

  # Reference category probability curve at base_eta = 0 for each step owner
  cpc_rows <- list()
  for (g in seq_len(enc$n_step_groups)) {
    owner_nm <- enc$step_owners[g]
    a_g <- 1.0
    if (!is.null(fit$slopes) && nrow(fit$slopes) > 0L && !is.null(enc$step_facet) &&
        enc$step_facet %in% enc$slope_facets) {
      sl_g <- fit$slopes[fit$slopes$Facet == enc$step_facet & fit$slopes$Level == owner_nm, ]
      if (nrow(sl_g) == 1L) a_g <- sl_g$Estimate[1L]
    }
    ev_g <- .mfrm_eval_row_predictions(
      eta_vec = th_grid,
      step_cum = fit$step_cumulative,
      step_group0 = rep(g - 1L, length(th_grid)),
      step_cat_min = enc$step_cat_min,
      step_cat_max = enc$step_cat_max,
      rating_min = enc$rating_min,
      alpha_vec = if (!is.null(fit$slopes)) rep(a_g, length(th_grid)) else NULL,
      theta_vec = th_grid,
      base_eta_vec = numeric(length(th_grid)),
      slope_formula = sform
    )
    for (k in seq_len(C)) {
      cpc_rows[[length(cpc_rows) + 1L]] <- data.frame(
        ScaleLevel = owner_nm,
        Theta = th_grid,
        Category = colnames(ev_g$probs)[k],
        Probability = ev_g$probs[, k],
        stringsAsFactors = FALSE
      )
    }
  }

  exposures <- data.frame(contexts, Count = counts, ObservedContext = observed_context, check.names = FALSE)
  names(exposures) <- make.unique(names(exposures))
  structure(
    list(
      test_information = data.frame(
        Theta = th_grid,
        Information = test_info,
        SE = ifelse(test_info > 0, 1.0 / sqrt(test_info), NA_real_),
        ExpectedScore = exp_total,
        stringsAsFactors = FALSE
      ),
      category_curves = do.call(rbind, cpc_rows),
      exposure_basis = exposure_basis,
      exposures = exposures,
      notes = c("Conditional on fixed calibration and independent ratings.",
                if (exposure_basis %in% c("observed_total", "per_person_average") && any(enc$weight != 1))
                  "Observed exposures use calibration weights; conditional SE assumes frequency weights.",
                if (any(!observed_context)) "Some planned combinations were not observed in calibration.")
    ),
    class = "mfrm_information"
  )
}


.mfrm_d_study_henderson <- function(x, design_grid = NULL, ...) {
  if (is.null(design_grid)) design_grid <- x$d_study[x$facets]
  facets <- x$facets
  design_grid <- .mfrm_design_grid(design_grid, facets)
  vc <- x$variance_components
  sigma2_p <- vc$Estimate[1L]
  # Parse non-person factors in each component name
  comp_parts <- strsplit(sub(",e$", "", vc$Component), ":", fixed = TRUE)
  eval_row <- function(n_vec) {
    s_rel <- 0.0
    s_abs <- 0.0
    for (j in seq_len(nrow(vc))[-1L]) {
      f_in <- intersect(comp_parts[[j]], facets)
      if (identical(vc$Component[j], "Residual")) f_in <- facets
      den <- if (length(f_in) == 0L) 1.0 else prod(n_vec[match(f_in, facets)])
      c_val <- vc$Estimate[j] / den
      s_abs <- s_abs + c_val
      if (isTRUE(vc$ContainsPerson[j])) s_rel <- s_rel + c_val
    }
    .mfrm_gstudy_metrics(sigma2_p, s_rel, s_abs)
  }
  mat <- do.call(rbind, lapply(seq_len(nrow(design_grid)), function(i)
    eval_row(unlist(design_grid[i, facets, drop = FALSE], use.names = FALSE))))
  out <- cbind(
    design_grid[, facets, drop = FALSE],
    TotalRatingsPerPerson = apply(as.matrix(design_grid[, facets, drop = FALSE]), 1L, prod),
    as.data.frame(mat, stringsAsFactors = FALSE)
  )
  out$MeetsG80 <- out$G >= 0.80
  out$MeetsPhi80 <- out$Phi >= 0.80
  rownames(out) <- NULL
  out
}


mfrm_d_study.mfrm_fit <- function(x, design_grid = NULL, ...) {
    .mfrm_validate_fit(x)
    enc <- x$encoding
    gt <- .mfrm_gstudy_henderson(data = enc$retained_data, person = enc$person_col, facets = enc$base_facets[seq_len(min(3L,
        length(enc$base_facets)))], score = enc$score_col, d_study_grid = design_grid)
    gt$d_study
}


print.mfrm_dff <- function(x, ...) {
  cat(sprintf(
    "Differential Facet Functioning [%s | Facet: %s | Group: %s (%s vs %s)]\n",
    x$method, x$facet %||% "", x$group %||% "", x$reference %||% "Ref", x$focal %||% "Focal"
  ))
  if (is.data.frame(x$table) && nrow(x$table) > 0L) {
    disp <- x$table
    num_c <- vapply(disp, is.numeric, logical(1))
    disp[num_c] <- lapply(disp[num_c], round, 4)
    print(disp, row.names = FALSE)
  }
  invisible(x)
}
