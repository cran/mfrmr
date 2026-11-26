# ==============================================================================
# Ultra-Fast Native C++ MFRM Engine, Affine Design Compiler, & S3 Model API
# ==============================================================================

.mfrm_sort_levels <- function(v) {
  if (is.factor(v)) return(levels(droplevels(v)))
  u <- unique(as.character(v))
  u <- u[!is.na(u) & nzchar(u)]
  num <- suppressWarnings(as.numeric(u))
  if (length(u) > 0L && !anyNA(num) && !anyDuplicated(num)) {
    u[order(num)]
  } else {
    sort(u)
  }
}

.mfrm_build_cpp_step_map <- function(step_constraints,
                                     owners,
                                     rating_min,
                                     K,
                                     step_cat_min,
                                     step_cat_max) {
  req_cols <- c("Level", "From", "To")
  opt_cols <- c("Value", "Group")
  if (!is.data.frame(step_constraints) || nrow(step_constraints) == 0L ||
      anyDuplicated(names(step_constraints)) ||
      !all(req_cols %in% names(step_constraints)) ||
      !any(opt_cols %in% names(step_constraints)) ||
      !all(names(step_constraints) %in% c(req_cols, opt_cols))) {
    stop(
      "`step_constraints` must be a data.frame with columns `Level`, `From`, `To`, and `Value` or `Group`.",
      call. = FALSE
    )
  }

  G <- length(owners)
  lvl_in <- as.character(step_constraints$Level)
  if (G == 1L && identical(owners, "Common")) {
    lvl_in[lvl_in %in% c("shared", "Shared", "COMMON", "common")] <- "Common"
  }
  owner_idx <- match(lvl_in, owners)
  from_seq <- rating_min + seq(0L, K - 1L)
  trans_idx <- match(as.numeric(step_constraints$From), from_seq)
  to_num <- as.numeric(step_constraints$To)

  if (anyNA(owner_idx) || anyNA(trans_idx) || anyNA(to_num) ||
      any(to_num != from_seq[trans_idx] + 1)) {
    stop(
      "`step_constraints` needs a known `Level` and consecutive `From`/`To` scores within the declared rating range.",
      call. = FALSE
    )
  }

  # Verify declared transitions lie within each owner's supported category range
  for (r in seq_len(nrow(step_constraints))) {
    g <- owner_idx[r]
    t0 <- trans_idx[r] - 1L
    if (t0 < step_cat_min[g] || t0 >= step_cat_max[g]) {
      stop(
        sprintf(
          "`step_constraints` transition %g->%g for Level '%s' lies outside its declared category_support [%g, %g].",
          step_constraints$From[r], step_constraints$To[r], owners[g],
          rating_min + step_cat_min[g], rating_min + step_cat_max[g]
        ),
        call. = FALSE
      )
    }
  }

  row_id <- (owner_idx - 1L) * K + trans_idx
  if (anyDuplicated(row_id)) {
    stop("Declare each `Level` and transition only once in `step_constraints`.", call. = FALSE)
  }

  val_vec <- if ("Value" %in% names(step_constraints)) step_constraints$Value else rep(NA_real_, nrow(step_constraints))
  grp_vec <- if ("Group" %in% names(step_constraints)) step_constraints$Group else rep(NA_character_, nrow(step_constraints))
  if (is.logical(val_vec) && all(is.na(val_vec))) val_vec <- as.numeric(val_vec)
  if (is.logical(grp_vec) && all(is.na(grp_vec))) grp_vec <- as.character(grp_vec)
  if (!is.numeric(val_vec) || any(is.nan(val_vec) | is.infinite(val_vec)) ||
      !is.character(grp_vec) || any(!is.na(grp_vec) & !nzchar(grp_vec)) ||
      any(is.na(val_vec) & is.na(grp_vec))) {
    stop("Each row of `step_constraints` needs a finite `Value` or a non-empty `Group`.", call. = FALSE)
  }

  n_rows <- G * K
  active_row <- logical(n_rows)
  for (g in seq_len(G)) {
    ag <- step_cat_min[g]
    bg <- step_cat_max[g]
    if (bg > ag) {
      active_row[(g - 1L) * K + seq.int(ag + 1L, bg)] <- TRUE
    }
  }

  act_indices <- which(active_row)
  n_act <- length(act_indices)
  row_to_act <- integer(n_rows)
  row_to_act[act_indices] <- seq_len(n_act)

  reps <- seq_len(n_act)
  non_na_grp <- which(!is.na(grp_vec))
  if (length(non_na_grp) > 0L) {
    for (mem in split(non_na_grp, grp_vec[non_na_grp])) {
      if (length(mem) < 2L || length(unique(trans_idx[mem])) != 1L) {
        stop("Each shared step `Group` must tie the same transition across at least two levels.", call. = FALSE)
      }
      vals_m <- unique(val_vec[mem][!is.na(val_vec[mem])])
      if (length(vals_m) > 1L) {
        stop("Conflicting fixed `Value`s in a shared step `Group`.", call. = FALSE)
      }
      act_m <- row_to_act[row_id[mem]]
      reps[act_m] <- min(act_m)
    }
  }

  cells <- sort(unique(reps))
  cell_of_act <- match(reps, cells)
  n_cells <- length(cells)
  fixed_cell <- rep(NA_real_, n_cells)

  # Single-transition (binary) groups must have their sole step fixed at 0
  for (g in seq_len(G)) {
    if (step_cat_max[g] - step_cat_min[g] == 1L) {
      r_bin <- (g - 1L) * K + (step_cat_min[g] + 1L)
      c_bin <- cell_of_act[row_to_act[r_bin]]
      fixed_cell[c_bin] <- 0.0
    }
  }

  non_na_val <- which(!is.na(val_vec))
  if (length(non_na_val) > 0L) {
    for (r in non_na_val) {
      c_idx <- cell_of_act[row_to_act[row_id[r]]]
      if (!is.na(fixed_cell[c_idx]) && abs(fixed_cell[c_idx] - val_vec[r]) > 1e-12) {
        stop("Fixed or shared steps conflict with within-ladder sum-to-zero constraints.", call. = FALSE)
      }
      fixed_cell[c_idx] <- val_vec[r]
    }
  }

  offset_cell <- replace(fixed_cell, is.na(fixed_cell), 0.0)
  avail <- which(is.na(fixed_cell))

  # Within-ladder sum-to-zero matrix over active owners with >= 1 active transition
  active_owners <- which((step_cat_max - step_cat_min) >= 1L)
  centering <- matrix(0.0, nrow = length(active_owners), ncol = n_cells)
  for (oi in seq_along(active_owners)) {
    g <- active_owners[oi]
    ag <- step_cat_min[g]
    bg <- step_cat_max[g]
    r_seq <- (g - 1L) * K + seq.int(ag + 1L, bg)
    for (rw in r_seq) {
      c_idx <- cell_of_act[row_to_act[rw]]
      centering[oi, c_idx] <- centering[oi, c_idx] + 1.0
    }
  }

  if (length(avail) > 0L && nrow(centering) > 0L) {
    sys <- centering[, rev(avail), drop = FALSE]
    qr_sys <- qr(sys)
    dep <- if (qr_sys$rank > 0L) rev(avail)[qr_sys$pivot[seq_len(qr_sys$rank)]] else integer(0)
    free <- setdiff(avail, dep)
  } else {
    dep <- integer(0)
    free <- avail
  }

  basis_cell <- matrix(0.0, nrow = n_cells, ncol = length(free))
  if (length(free) > 0L) {
    basis_cell[cbind(free, seq_along(free))] <- 1.0
  }
  if (length(dep) > 0L) {
    calc <- qr.coef(
      qr(centering[, dep, drop = FALSE]),
      -cbind(centering %*% offset_cell, centering[, free, drop = FALSE])
    )
    offset_cell[dep] <- calc[, 1L]
    if (length(free) > 0L) {
      basis_cell[dep, ] <- calc[, -1L, drop = FALSE]
    }
  }
  basis_cell[abs(basis_cell) < 1e-12] <- 0.0

  if (nrow(centering) > 0L) {
    lhs <- as.numeric(centering %*% offset_cell)
    sc <- as.numeric(centering %*% abs(offset_cell))
    if (any(!is.finite(lhs)) || any(abs(lhs) > 1e-10 * pmax(1.0, sc)) ||
        (ncol(basis_cell) > 0L && any(abs(centering %*% basis_cell) > 1e-10))) {
      stop("Fixed or shared steps conflict with within-ladder sum-to-zero constraints.", call. = FALSE)
    }
  }

  d0_full <- rep(NA_real_, n_rows)
  B_full <- matrix(0.0, nrow = n_rows, ncol = ncol(basis_cell))
  fixed_full <- rep(1L, n_rows)

  d0_full[act_indices] <- offset_cell[cell_of_act]
  if (ncol(basis_cell) > 0L) {
    B_full[act_indices, ] <- basis_cell[cell_of_act, , drop = FALSE]
  }
  # A step row has sampling variance iff at least one basis column is nonzero
  if (ncol(B_full) > 0L) {
    has_var <- rowSums(abs(B_full) > 1e-14) > 0L
    fixed_full[act_indices[has_var]] <- 0L
  }

  list(
    q = as.integer(ncol(B_full)),
    d0 = as.double(d0_full),
    B = B_full,
    fixed = as.integer(fixed_full),
    declarations = step_constraints
  )
}

# Compile individual anchors + group_anchors + centering into exact affine map beta_f = d0 + B * xi
.mfrm_build_cpp_group_anchor_map <- function(levels_f,
                                             anchor_fixed_f,
                                             anchor_val_f,
                                             group_rows,
                                             is_noncenter,
                                             centering_scope) {
  Lf <- length(levels_f)
  fixed_val <- rep(NA_real_, Lf)
  anc_idx <- which(anchor_fixed_f == 1L)
  if (length(anc_idx) > 0L) fixed_val[anc_idx] <- anchor_val_f[anc_idx]

  l_col <- intersect(c("Level", "level", "Element", "element"), names(group_rows))[1L]
  g_col <- intersect(c("Group", "group"), names(group_rows))[1L]
  v_col <- intersect(c("GroupValue", "Value", "Anchor", "value", "anchor"), names(group_rows))[1L]
  if (is.na(l_col) || is.na(g_col) || is.na(v_col)) {
    stop("`group_anchors` must contain columns `Facet`, `Level`, `Group`, and `GroupValue` (or `Value`).", call. = FALSE)
  }

  C_rows <- list()
  c_rhs <- numeric(0)
  anchored_mask <- as.integer(anchor_fixed_f == 1L)

  for (ii in anc_idx) {
    r <- numeric(Lf); r[ii] <- 1.0
    C_rows[[length(C_rows) + 1L]] <- r
    c_rhs <- c(c_rhs, anchor_val_f[ii])
  }

  grp_split <- split(seq_len(nrow(group_rows)), as.character(group_rows[[g_col]]))
  for (gnm in names(grp_split)) {
    r_idx <- grp_split[[gnm]]
    lv_idx <- match(as.character(group_rows[[l_col]][r_idx]), levels_f)
    lv_idx <- unique(lv_idx[!is.na(lv_idx)])
    if (length(lv_idx) == 0L) next
    gvals <- unique(as.numeric(group_rows[[v_col]][r_idx]))
    gvals <- gvals[is.finite(gvals)]
    if (length(gvals) != 1L) {
      stop("Each `group_anchors` group must declare a single finite `GroupValue`.", call. = FALSE)
    }
    r <- numeric(Lf); r[lv_idx] <- 1.0
    C_rows[[length(C_rows) + 1L]] <- r
    c_rhs <- c(c_rhs, length(lv_idx) * gvals[1L])
    anchored_mask[lv_idx] <- 1L
  }

  if (!isTRUE(is_noncenter)) {
    if (identical(centering_scope, "facet")) {
      C_rows[[length(C_rows) + 1L]] <- rep(1.0, Lf)
      c_rhs <- c(c_rhs, 0.0)
    } else {
      rem <- which(anchored_mask == 0L)
      if (length(rem) > 0L) {
        r <- numeric(Lf); r[rem] <- 1.0
        C_rows[[length(C_rows) + 1L]] <- r
        c_rhs <- c(c_rhs, 0.0)
      }
    }
  }

  if (length(C_rows) == 0L) {
    return(list(q = Lf, d0 = numeric(Lf), B = diag(1.0, Lf), anchored_mask = anchored_mask))
  }
  C_mat <- do.call(rbind, C_rows)
  # Solve C_mat %*% beta = c_rhs via pivoted QR onrev(1:Lf) so trailing element is dependent
  ord <- rev(seq_len(Lf))
  C_ord <- C_mat[, ord, drop = FALSE]
  qr_c <- qr(C_ord)
  rnk <- qr_c$rank
  dep <- if (rnk > 0L) ord[qr_c$pivot[seq_len(rnk)]] else integer(0)
  free <- setdiff(seq_len(Lf), dep)

  d0 <- numeric(Lf)
  B <- matrix(0.0, nrow = Lf, ncol = length(free))
  if (length(free) > 0L) {
    B[cbind(free, seq_along(free))] <- 1.0
  }
  if (length(dep) > 0L) {
    # Keep independent constraint rows via QR of transpose
    qr_t <- qr(t(C_mat))
    ind_rows <- qr_t$pivot[seq_len(qr_t$rank)]
    C_ind <- C_mat[ind_rows, , drop = FALSE]
    rhs_ind <- c_rhs[ind_rows]
    sol <- qr.coef(
      qr(C_ind[, dep, drop = FALSE]),
      cbind(rhs_ind, -C_ind[, free, drop = FALSE])
    )
    if (is.null(dim(sol))) sol <- matrix(sol, nrow = 1L)
    d0[dep] <- sol[, 1L]
    if (length(free) > 0L) B[dep, ] <- sol[, -1L, drop = FALSE]
  }
  resid_c <- as.numeric(C_mat %*% d0 - c_rhs)
  if (any(!is.finite(resid_c)) || max(abs(resid_c)) > 1e-9) {
    stop("Conflicting `anchors`, `group_anchors`, and centering constraints.", call. = FALSE)
  }
  B[abs(B) < 1e-12] <- 0.0
  storage.mode(B) <- "double"
  list(q = as.integer(ncol(B)), d0 = as.double(d0), B = B, anchored_mask = as.integer(anchored_mask))
}

.mfrm_unique_covariate_ids <- function(ids, source) {
  ids <- as.character(ids)
  duplicate <- anyDuplicated(ids)
  if (duplicate) {
    stop(source, " must have one row per ID. Duplicate ID: '", ids[duplicate], "'.", call. = FALSE)
  }
  ids
}

.mfrm_encode_dataset <- function(data,
                                      person,
                                      facets,
                                      score,
                                      model = "RSM",
                                      step_facet = NULL,
                                      slope_facet = NULL,
                                      slope_reference = NULL,
                                      slope_type = "full_predictor",
                                      gpcm_mml_identification = "free_population",
                                      method_in = "MML",
                                      rating_min = NULL,
                                      rating_max = NULL,
                                      weight = NULL,
                                      anchors = NULL,
                                      group_anchors = NULL,
                                      dummy_facets = NULL,
                                      noncenter_facet = "Person",
                                      positive_facets = NULL,
                                      centering_scope = "remainder",
                                      population_formula = NULL,
                                      person_data = NULL,
                                      person_id = NULL,
                                      population_sd = 1.0,
                                      step_constraints = NULL,
                                      category_support = NULL,
                                      facet_formulas = NULL,
                                      facet_data = NULL,
                                      facet_interactions = NULL) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }
  if (!is.character(person) || length(person) != 1L || !nzchar(person)) {
    stop("`person` must be a single column name.", call. = FALSE)
  }
  if (!is.character(facets) || length(facets) < 1L || any(!nzchar(facets))) {
    stop("`facets` must be a non-empty character vector of facet column names.", call. = FALSE)
  }
  req_cols <- c(person, facets, score, weight)
  miss_cols <- setdiff(req_cols, names(data))
  if (length(miss_cols) > 0L) {
    stop("Missing column(s) in `data`: ", paste(miss_cols, collapse = ", "), ".", call. = FALSE)
  }

  if (is.factor(data[[score]])) {
    stop("The score column must contain numeric rating values. Convert factor labels to numeric scores before calibration.", call. = FALSE)
  }
  raw_y <- as.numeric(data[[score]])
  w_vec <- if (is.null(weight)) rep(1.0, nrow(data)) else as.numeric(data[[weight]])
  keep <- is.finite(raw_y) & !is.na(data[[person]]) & is.finite(w_vec) & (w_vec > 0)
  for (fc in facets) keep <- keep & !is.na(data[[fc]])
  if (!all(keep)) {
    data <- data[keep, , drop = FALSE]
    raw_y <- raw_y[keep]
    w_vec <- w_vec[keep]
  }
  if (nrow(data) == 0L) {
    stop("No complete finite rating rows remain after omitting missing values.", call. = FALSE)
  }

  rmin <- if (is.null(rating_min)) min(raw_y) else as.numeric(rating_min)[1L]
  rmax <- if (is.null(rating_max)) max(raw_y) else as.numeric(rating_max)[1L]
  if (!is.null(category_support) && is.data.frame(category_support)) {
    if ("Min" %in% names(category_support) && is.null(rating_min)) {
      rmin <- min(rmin, min(as.numeric(category_support$Min), na.rm = TRUE))
    }
    if ("Max" %in% names(category_support) && is.null(rating_max)) {
      rmax <- max(rmax, max(as.numeric(category_support$Max), na.rm = TRUE))
    }
  }
  if (!is.finite(rmin) || !is.finite(rmax) || rmax <= rmin ||
      abs((rmax - rmin) - round(rmax - rmin)) > 1e-8) {
    stop("`rating_min` and `rating_max` must define an integer span with `rating_max > rating_min`.", call. = FALSE)
  }
  if (any(abs(raw_y - round(raw_y)) > 1e-8)) {
    stop("Observed scores must be integer-valued categories.", call. = FALSE)
  }
  score0 <- as.integer(round(raw_y - rmin))
  n_cat <- as.integer(round(rmax - rmin) + 1L)
  K <- n_cat - 1L
  if (any(score0 < 0L | score0 >= n_cat)) {
    stop("Observed scores fall outside [`rating_min`, `rating_max`].", call. = FALSE)
  }

  p_chr <- as.character(data[[person]])
  p_levels <- .mfrm_sort_levels(data[[person]])
  person0 <- match(p_chr, p_levels) - 1L
  n_person <- length(p_levels)

  base_facets <- as.character(facets)
  n_base_facets <- length(base_facets)
  facet_levels <- vector("list", n_base_facets)
  names(facet_levels) <- base_facets
  facet_n_levels <- integer(n_base_facets)
  facet_mat0 <- matrix(0L, nrow = nrow(data), ncol = n_base_facets)
  for (f in seq_len(n_base_facets)) {
    fc <- base_facets[f]
    lv <- .mfrm_sort_levels(data[[fc]])
    facet_levels[[f]] <- lv
    facet_n_levels[f] <- length(lv)
    facet_mat0[, f] <- match(as.character(data[[fc]]), lv) - 1L
  }

  facet_signs <- ifelse(base_facets %in% as.character(positive_facets), 1.0, -1.0)
  facet_dummy <- as.integer(base_facets %in% as.character(dummy_facets))
  noncenter_chr <- if (is.null(noncenter_facet)) "Person" else as.character(noncenter_facet)
  facet_noncenter <- as.integer(base_facets %in% setdiff(noncenter_chr, "Person"))
  cscope_int <- if (identical(centering_scope, "facet")) 1L else 0L

  anchor_fixed <- vector("list", n_base_facets)
  anchor_val <- vector("list", n_base_facets)
  for (f in seq_len(n_base_facets)) {
    anchor_fixed[[f]] <- integer(facet_n_levels[f])
    anchor_val[[f]] <- numeric(facet_n_levels[f])
  }
  if (!is.null(anchors) && is.data.frame(anchors) && nrow(anchors) > 0L) {
    f_col <- intersect(c("Facet", "facet"), names(anchors))[1L]
    l_col <- intersect(c("Level", "level", "Element", "element"), names(anchors))[1L]
    v_col <- intersect(c("Anchor", "anchor", "Value", "value", "Estimate"), names(anchors))[1L]
    if (is.na(f_col) || is.na(l_col) || is.na(v_col)) {
      stop("`anchors` must contain columns `Facet`, `Level`, and `Anchor`.", call. = FALSE)
    }
    for (r in seq_len(nrow(anchors))) {
      fnm <- as.character(anchors[[f_col]][r])
      fidx <- match(fnm, base_facets)
      if (is.na(fidx)) next
      lidx <- match(as.character(anchors[[l_col]][r]), facet_levels[[fidx]])
      if (is.na(lidx)) next
      aval <- as.numeric(anchors[[v_col]][r])
      if (is.finite(aval)) {
        anchor_fixed[[fidx]][lidx] <- 1L
        anchor_val[[fidx]][lidx] <- aval
      }
    }
  }

  # Custom affine facet maps: Explanatory RSM/PCM/GPCM & group_anchors
  has_custom_facet_map <- integer(n_base_facets)
  custom_facet_q <- integer(n_base_facets)
  custom_facet_d0 <- lapply(facet_n_levels, numeric)
  custom_facet_B <- lapply(facet_n_levels, function(L) matrix(0.0, nrow = L, ncol = 0L))
  custom_facet_anchored <- vector("list", n_base_facets)
  group_anchor_meta <- list()
  explanatory_meta <- list()

  model_up <- toupper(model[1L])
  # For 2+ family GPCM under fixed N(0,1) without anchors or population_formula,
  # absorb grand location mean into the free-slope facet when noncenter_facet == "Person"
  if (identical(model_up, "GPCM") && identical(toupper(method_in[1L]), "MML") &&
      length(slope_facet) >= 2L && identical(noncenter_chr, "Person") &&
      is.null(population_formula) &&
      all(vapply(anchor_fixed, function(x) all(x == 0L), logical(1))) &&
      is.null(group_anchors)) {
    s_ref_tmp <- if (is.null(slope_reference)) slope_facet[1L] else as.character(slope_reference)[1L]
    s_free_tmp <- setdiff(as.character(slope_facet), s_ref_tmp)[1L]
    f_nc_auto <- match(s_free_tmp, base_facets)
    if (!is.na(f_nc_auto)) facet_noncenter[f_nc_auto] <- 1L
  }

  if (!is.null(group_anchors) && is.data.frame(group_anchors) && nrow(group_anchors) > 0L) {
    f_col_ga <- intersect(c("Facet", "facet"), names(group_anchors))[1L]
    if (is.na(f_col_ga)) {
      stop("`group_anchors` must contain column `Facet`.", call. = FALSE)
    }
    for (fnm in unique(as.character(group_anchors[[f_col_ga]]))) {
      fidx <- match(fnm, base_facets)
      if (is.na(fidx) || facet_dummy[fidx] == 1L) next
      if (!is.null(facet_formulas) && fnm %in% names(facet_formulas)) {
        stop("Cannot combine `group_anchors` and `facet_formulas` on the same facet ('", fnm, "').", call. = FALSE)
      }
      g_sub <- group_anchors[as.character(group_anchors[[f_col_ga]]) == fnm, , drop = FALSE]
      ga_map <- .mfrm_build_cpp_group_anchor_map(
        levels_f = facet_levels[[fidx]],
        anchor_fixed_f = anchor_fixed[[fidx]],
        anchor_val_f = anchor_val[[fidx]],
        group_rows = g_sub,
        is_noncenter = (facet_noncenter[fidx] == 1L),
        centering_scope = centering_scope
      )
      has_custom_facet_map[fidx] <- 1L
      custom_facet_q[fidx] <- ga_map$q
      custom_facet_d0[[fidx]] <- ga_map$d0
      custom_facet_B[[fidx]] <- ga_map$B
      custom_facet_anchored[[fidx]] <- ga_map$anchored_mask
      group_anchor_meta[[fnm]] <- list(facet = fnm, facet_idx = fidx, q = ga_map$q)
    }
  }

  if (!is.null(facet_formulas)) {
    if (!is.list(facet_formulas) || is.null(names(facet_formulas)) ||
        any(!names(facet_formulas) %in% base_facets)) {
      stop("`facet_formulas` must be a named list of one-sided formulas keyed by `facets`.", call. = FALSE)
    }
    for (fnm in names(facet_formulas)) {
      fidx <- match(fnm, base_facets)
      if (facet_dummy[fidx] == 1L) {
        stop("Cannot specify `facet_formulas` on dummy facet '", fnm, "'.", call. = FALSE)
      }
      if (any(anchor_fixed[[fidx]] == 1L)) {
        stop("Cannot combine fixed `anchors` and `facet_formulas` on the same facet ('", fnm, "').", call. = FALSE)
      }
      fml <- facet_formulas[[fnm]]
      if (!inherits(fml, "formula")) {
        stop("`facet_formulas[['", fnm, "']]` must be a one-sided formula.", call. = FALSE)
      }
      lvs <- facet_levels[[fidx]]
      Lf <- length(lvs)
      fdf <- if (!is.null(facet_data) && is.list(facet_data) && !is.data.frame(facet_data)) {
        facet_data[[fnm]]
      } else if (is.data.frame(facet_data)) {
        facet_data
      } else {
        NULL
      }
      if (is.null(fdf) || !is.data.frame(fdf)) {
        stop("`facet_data` must supply a data.frame of element covariates for facet '", fnm, "'.", call. = FALSE)
      }
      id_col <- intersect(c(fnm, "Level", "Element"), names(fdf))[1L]
      if (!is.na(id_col)) {
        m_idx <- match(lvs, .mfrm_unique_covariate_ids(fdf[[id_col]],
          paste0("`facet_data` for facet '", fnm, "'")))
      } else if (!is.null(rownames(fdf))) {
        m_idx <- match(lvs, rownames(fdf))
      } else {
        stop("`facet_data[['", fnm, "']]` must include column '", fnm, "' or `'Level'`.", call. = FALSE)
      }
      if (anyNA(m_idx)) {
        stop("`facet_data[['", fnm, "']]` is missing observed level(s): ",
             paste(lvs[is.na(m_idx)], collapse = ", "), ".", call. = FALSE)
      }
      fdf_ord <- fdf[m_idx, , drop = FALSE]
      tt <- stats::delete.response(stats::terms(fml))
      mm <- stats::model.matrix(tt, data = fdf_ord)
      if (any(!is.finite(mm))) {
        stop("Covariate design matrix for facet '", fnm, "' contains missing or non-finite values.", call. = FALSE)
      }
      is_centered <- (facet_noncenter[fidx] == 0L)
      col_means <- numeric(0)
      if (is_centered) {
        keep_c <- colnames(mm) != "(Intercept)"
        mm <- mm[, keep_c, drop = FALSE]
        if (ncol(mm) == 0L) {
          stop("Centered explanatory formula for facet '", fnm, "' has no non-intercept predictors.", call. = FALSE)
        }
        col_means <- colMeans(mm)
        mm <- sweep(mm, 2L, col_means, "-")
        max_q <- Lf - 1L
      } else {
        max_q <- Lf
      }
      qr_mm <- qr(mm)
      if (ncol(mm) > max_q || qr_mm$rank < ncol(mm)) {
        stop(
          sprintf(
            "Explanatory design matrix for facet '%s' is rank-deficient or overparameterized (ncol = %d, rank = %d, max identifiable = %d).",
            fnm, ncol(mm), qr_mm$rank, max_q
          ),
          call. = FALSE
        )
      }
      storage.mode(mm) <- "double"
      has_custom_facet_map[fidx] <- 1L
      custom_facet_q[fidx] <- ncol(mm)
      custom_facet_d0[[fidx]] <- numeric(Lf)
      custom_facet_B[[fidx]] <- mm
      explanatory_meta[[fnm]] <- list(
        facet = fnm,
        facet_idx = fidx,
        formula = fml,
        terms = tt,
        centered = is_centered,
        col_means = col_means,
        colnames = colnames(mm),
        B = mm
      )
    }
  }

  # Optional joint facet_interactions: append synthetic zero-margin Kronecker facets
  interaction_meta <- list()
  all_facets <- base_facets
  if (!is.null(facet_interactions)) {
    int_pairs <- if (is.character(facet_interactions) && length(facet_interactions) == 2L) {
      list(facet_interactions)
    } else if (is.list(facet_interactions)) {
      facet_interactions
    } else {
      stop("`facet_interactions` must be a list of length-2 facet name vectors.", call. = FALSE)
    }
    for (pr in int_pairs) {
      pr <- as.character(pr)
      if (length(pr) != 2L || pr[1L] == pr[2L] || any(!pr %in% base_facets)) {
        stop("Each entry in `facet_interactions` must name two distinct facets in `facets`.", call. = FALSE)
      }
      fa <- match(pr[1L], base_facets)
      fb <- match(pr[2L], base_facets)
      La <- facet_n_levels[fa]
      Lb <- facet_n_levels[fb]
      if (La < 2L || Lb < 2L) next
      sum_zero_basis <- function(L) {
        B <- rbind(diag(1.0, L - 1L), rep(-1.0, L - 1L))
        B
      }
      Ba <- sum_zero_basis(La)
      Bb <- sum_zero_basis(Lb)
      # Kronecker basis over cell index l = la * Lb + lb (0-based, row-major)
      B_ab <- kronecker(Ba, Bb)
      storage.mode(B_ab) <- "double"
      cell_levels <- as.vector(t(outer(facet_levels[[fa]], facet_levels[[fb]], paste, sep = ":")))
      cell_idx0 <- facet_mat0[, fa] * Lb + facet_mat0[, fb]
      int_name <- paste0(pr[1L], ":", pr[2L])

      all_facets <- c(all_facets, int_name)
      facet_levels[[int_name]] <- cell_levels
      facet_n_levels <- c(facet_n_levels, La * Lb)
      facet_mat0 <- cbind(facet_mat0, as.integer(cell_idx0))
      # Sign matches FacetA orientation so positive interaction = higher severity/difficulty when sign=-1
      facet_signs <- c(facet_signs, facet_signs[fa])
      facet_dummy <- c(facet_dummy, 0L)
      facet_noncenter <- c(facet_noncenter, 0L)
      anchor_fixed[[length(all_facets)]] <- integer(La * Lb)
      anchor_val[[length(all_facets)]] <- numeric(La * Lb)
      has_custom_facet_map <- c(has_custom_facet_map, 1L)
      custom_facet_q <- c(custom_facet_q, ncol(B_ab))
      custom_facet_d0[[length(all_facets)]] <- numeric(La * Lb)
      custom_facet_B[[length(all_facets)]] <- B_ab
      custom_facet_anchored[[length(all_facets)]] <- integer(La * Lb)

      interaction_meta[[int_name]] <- list(
        name = int_name,
        facet_a = pr[1L],
        facet_b = pr[2L],
        facet_idx = length(all_facets),
        levels_a = rep(facet_levels[[fa]], each = Lb),
        levels_b = rep(facet_levels[[fb]], times = La),
        cell_levels = cell_levels,
        q = ncol(B_ab),
        B = B_ab
      )
    }
  }

  # Step facet & GPCM / GMFRM slope design compilation
  is_gpcm <- identical(model_up, "GPCM")
  if (is_gpcm && toupper(method_in[1L]) %in% c("CML", "JML_JEFFREYS")) {
    stop("`model = 'GPCM'` does not have a valid parameter-free sufficient statistic for CML or Jeffreys JML; use `method = 'MML'` (recommended) or `'JML'`.", call. = FALSE)
  }

  n_slope_facets <- 0L
  slope_formula_code <- 0L
  slope_facets0 <- integer(0)
  slope_q <- integer(0)
  slope_d0 <- list()
  slope_B <- list()
  slope_fixed <- list()
  slope_facet_names <- character(0)
  slope_type <- match.arg(slope_type[1L], c("full_predictor", "ability_only"))
  gpcm_id <- match.arg(gpcm_mml_identification[1L],
                       c("free_population", "fixed_standard_normal", "unconstrained_slopes"))

  if (is_gpcm) {
    sf_in <- if (is.null(slope_facet)) {
      if (!is.null(step_facet) && step_facet[1L] %in% base_facets) as.character(step_facet)[1L] else base_facets[length(base_facets)]
    } else {
      unique(as.character(slope_facet))
    }
    if (length(sf_in) < 1L || any(!sf_in %in% base_facets)) {
      stop("`slope_facet` must name one or more facets declared in `facets`.", call. = FALSE)
    }
    if (any(facet_dummy[match(sf_in, base_facets)] == 1L)) {
      stop("Cannot estimate discrimination slopes on `dummy_facets`.", call. = FALSE)
    }
    slope_facet_names <- sf_in
    n_slope_facets <- length(sf_in)
    slope_facets0 <- as.integer(match(sf_in, base_facets) - 1L)
    slope_formula_code <- if (identical(slope_type, "ability_only")) 1L else 0L

    s_ref <- if (is.null(slope_reference)) sf_in[1L] else as.character(slope_reference)[1L]
    if (!s_ref %in% sf_in) {
      stop("`slope_reference` ('", s_ref, "') must be one of `slope_facet`.", call. = FALSE)
    }

    for (s in seq_len(n_slope_facets)) {
      fnm_s <- sf_in[s]
      f_idx_s <- slope_facets0[s] + 1L
      Ls <- facet_n_levels[f_idx_s]
      # Determine whether family s is constrained to GM(alpha) = 1 (sum(log_alpha) = 0) or free
      is_constrained_gm1 <- if (!identical(toupper(method_in[1L]), "MML")) {
        TRUE
      } else if (n_slope_facets == 1L) {
        !identical(gpcm_id, "unconstrained_slopes")
      } else {
        # For 2+ slope families under fixed N(0,1): reference family (and any 3rd+) has GM=1;
        # the primary non-reference family is free on the N(0,1) latent scale!
        non_ref <- setdiff(sf_in, s_ref)
        !(identical(fnm_s, non_ref[1L]) && is.null(population_formula))
      }
      slope_d0[[s]] <- numeric(Ls)
      if (is_constrained_gm1) {
        if (Ls < 2L) {
          slope_q[s] <- 0L
          slope_B[[s]] <- matrix(0.0, nrow = Ls, ncol = 0L)
          slope_fixed[[s]] <- rep(1L, Ls)
        } else {
          qs <- Ls - 1L
          slope_q[s] <- qs
          Bs <- rbind(diag(1.0, qs), rep(-1.0, qs))
          storage.mode(Bs) <- "double"
          slope_B[[s]] <- Bs
          slope_fixed[[s]] <- integer(Ls)
        }
      } else {
        slope_q[s] <- Ls
        Bs <- diag(1.0, Ls)
        storage.mode(Bs) <- "double"
        slope_B[[s]] <- Bs
        slope_fixed[[s]] <- integer(Ls)
      }
    }

    # For 1-family GPCM MML with default `free_population` and no explicit population_formula,
    # activate ~ 1 latent population (free mean + free SD, or free SD if a facet is noncentered/anchored)
    if (identical(toupper(method_in[1L]), "MML") && n_slope_facets == 1L &&
        identical(gpcm_id, "free_population") && is.null(population_formula)) {
      any_anc <- any(vapply(anchor_fixed, function(x) any(x == 1L), logical(1))) || length(group_anchor_meta) > 0L
      if (any(facet_noncenter == 1L) || any_anc) {
        population_formula <- ~ 0
      } else {
        population_formula <- ~ 1
      }
    }
  }

  use_pcm_steps <- identical(model_up, "PCM") ||
    (is_gpcm && !(identical(step_facet, "Common") || identical(step_facet, "RSM"))) ||
    (!is.null(category_support) && !is.null(step_facet))
  model_code <- if (use_pcm_steps) 1L else 0L
  step_f_idx <- -1L
  if (model_code == 1L) {
    sf <- if (is.null(step_facet)) {
      if (is_gpcm && length(slope_facet_names) > 0L) slope_facet_names[length(slope_facet_names)] else base_facets[length(base_facets)]
    } else {
      as.character(step_facet)[1L]
    }
    step_f_idx <- match(sf, base_facets) - 1L
    if (is.na(step_f_idx)) {
      stop("`step_facet` ('", sf, "') must be one of `facets`.", call. = FALSE)
    }
  }

  step_owners <- if (model_code == 1L) facet_levels[[step_f_idx + 1L]] else "Common"
  n_step_groups <- length(step_owners)
  step_group0 <- if (model_code == 1L) as.integer(facet_mat0[, step_f_idx + 1L]) else integer(nrow(data))

  # Declarative category_support: per-owner [Min, Max] integer bounds
  step_cat_min <- integer(n_step_groups)
  step_cat_max <- rep(as.integer(K), n_step_groups)
  if (!is.null(category_support)) {
    if (!is.data.frame(category_support) || !all(c("Min", "Max") %in% names(category_support))) {
      stop("`category_support` must be a data.frame containing columns `Level` (for PCM), `Min`, and `Max`.", call. = FALSE)
    }
    lvl_cs <- if ("Level" %in% names(category_support)) as.character(category_support$Level) else rep("Common", nrow(category_support))
    if (n_step_groups == 1L && identical(step_owners, "Common")) {
      lvl_cs[lvl_cs %in% c("shared", "Shared", "COMMON", "common")] <- "Common"
    }
    if (anyDuplicated(lvl_cs)) {
      stop("Duplicate `Level` entries in `category_support`.", call. = FALSE)
    }
    bad_lv <- setdiff(lvl_cs, step_owners)
    if (length(bad_lv) > 0L) {
      stop("Unknown `Level`(s) in `category_support`: ", paste(bad_lv, collapse = ", "), ".", call. = FALSE)
    }
    for (r in seq_len(nrow(category_support))) {
      g <- match(lvl_cs[r], step_owners)
      mn <- as.numeric(category_support$Min[r])
      mx <- as.numeric(category_support$Max[r])
      if (!is.finite(mn) || !is.finite(mx) || mx <= mn ||
          abs(mn - round(mn)) > 1e-8 || abs(mx - round(mx)) > 1e-8 ||
          mn < rmin || mx > rmax) {
        stop(
          sprintf("Invalid `category_support` bounds [%g, %g] for Level '%s' (must be integers inside [%g, %g] with Max > Min).",
                  mn, mx, step_owners[g], rmin, rmax),
          call. = FALSE
        )
      }
      step_cat_min[g] <- as.integer(round(mn - rmin))
      step_cat_max[g] <- as.integer(round(mx - rmin))
    }
    # Check observed ratings obey declared category support
    obs_min_req <- step_cat_min[step_group0 + 1L]
    obs_max_req <- step_cat_max[step_group0 + 1L]
    viol <- which(score0 < obs_min_req | score0 > obs_max_req)
    if (length(viol) > 0L) {
      i0 <- viol[1L]
      g0 <- step_group0[i0] + 1L
      stop(
        sprintf(
          "Observed score %g at row %d violates declared `category_support` [%g, %g] for Level '%s'.",
          raw_y[i0], i0, rmin + step_cat_min[g0], rmin + step_cat_max[g0], step_owners[g0]
        ),
        call. = FALSE
      )
    }
  }

  # Optional custom affine step map (`step_constraints`)
  has_custom_step_map <- FALSE
  custom_step_q <- 0L
  custom_step_d0 <- numeric(0)
  custom_step_B <- matrix(0.0, nrow = 0L, ncol = 0L)
  custom_step_fixed <- integer(0)
  step_map_info <- NULL

  if (!is.null(step_constraints)) {
    step_map_info <- .mfrm_build_cpp_step_map(
      step_constraints = step_constraints,
      owners = step_owners,
      rating_min = rmin,
      K = K,
      step_cat_min = step_cat_min,
      step_cat_max = step_cat_max
    )
    has_custom_step_map <- TRUE
    custom_step_q <- step_map_info$q
    custom_step_d0 <- step_map_info$d0
    custom_step_B <- step_map_info$B
    custom_step_fixed <- step_map_info$fixed
  }

  pop_mode <- 0L
  pop_col_names <- character(0)
  pop_X_mat <- matrix(0.0, nrow = 0L, ncol = 0L)
  pop_terms <- NULL
  if (!is.null(population_formula)) {
    pop_mode <- 1L
    if (is.null(person_data)) {
      pdata <- data.frame(.person_id = p_levels, stringsAsFactors = FALSE)
      names(pdata)[1L] <- person
    } else {
      pid_col <- if (is.null(person_id)) person else person_id
      if (!pid_col %in% names(person_data)) {
        stop("`person_data` must contain person ID column '", pid_col, "'.", call. = FALSE)
      }
      idx_m <- match(p_levels, .mfrm_unique_covariate_ids(person_data[[pid_col]], "`person_data`"))
      if (anyNA(idx_m)) {
        stop("`person_data` is missing one or more observed Person IDs.", call. = FALSE)
      }
      pdata <- person_data[idx_m, , drop = FALSE]
    }
    pop_terms <- stats::delete.response(stats::terms(population_formula))
    mm <- stats::model.matrix(pop_terms, data = pdata)
    if (any(facet_noncenter == 1L) && "(Intercept)" %in% colnames(mm)) {
      mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    }
    storage.mode(mm) <- "double"
    pop_X_mat <- mm
    pop_col_names <- colnames(mm)
  }

  enc <- list(
    person_col = person,
    score_col = score,
    weight_col = weight,
    n_person = as.integer(n_person),
    person_levels = p_levels,
    base_facets = base_facets,
    facets = all_facets,
    facet_levels = facet_levels,
    facet_n_levels = as.integer(facet_n_levels),
    facet_signs = as.double(facet_signs),
    facet_dummy = as.integer(facet_dummy),
    facet_noncenter = as.integer(facet_noncenter),
    centering_scope = as.integer(cscope_int),
    anchor_fixed = anchor_fixed,
    anchor_val = anchor_val,
    has_custom_facet_map = as.integer(has_custom_facet_map),
    custom_facet_q = as.integer(custom_facet_q),
    custom_facet_d0 = custom_facet_d0,
    custom_facet_B = custom_facet_B,
    custom_facet_anchored = custom_facet_anchored,
    group_anchor_meta = group_anchor_meta,
    group_anchors = group_anchors,
    explanatory_meta = explanatory_meta,
    interaction_meta = interaction_meta,
    n_slope_facets = as.integer(n_slope_facets),
    slope_formula = as.integer(slope_formula_code),
    slope_type = slope_type,
    slope_facets = slope_facet_names,
    slope_reference = if (length(slope_facet_names) > 1L) slope_reference,
    slope_facets0 = as.integer(slope_facets0),
    slope_q = as.integer(slope_q),
    slope_d0 = slope_d0,
    slope_B = slope_B,
    slope_fixed = slope_fixed,
    gpcm_mml_identification = gpcm_id,
    n_cat = as.integer(n_cat),
    rating_min = as.double(rmin),
    rating_max = as.double(rmax),
    model = model_up,
    model_code = as.integer(model_code),
    step_facet = if (model_code == 1L) base_facets[step_f_idx + 1L] else NULL,
    step_facet_idx = as.integer(step_f_idx),
    step_owners = step_owners,
    n_step_groups = as.integer(n_step_groups),
    step_group0 = as.integer(step_group0),
    step_cat_min = as.integer(step_cat_min),
    step_cat_max = as.integer(step_cat_max),
    has_custom_step_map = isTRUE(has_custom_step_map),
    custom_step_q = as.integer(custom_step_q),
    custom_step_d0 = as.double(custom_step_d0),
    custom_step_B = custom_step_B,
    custom_step_fixed = as.integer(custom_step_fixed),
    step_map_info = step_map_info,
    category_support = category_support,
    person0 = as.integer(person0),
    score0 = as.integer(score0),
    weight = as.double(w_vec),
    facet_mat0 = facet_mat0,
    pop_mode = as.integer(pop_mode),
    fixed_sigma = as.double(population_sd),
    pop_X_mat = pop_X_mat,
    pop_col_names = pop_col_names,
    pop_formula = population_formula,
    pop_terms = pop_terms,
    person_id = if (is.null(person_id)) person else person_id,
    retained_data = data,
    n_obs = nrow(data)
  )
  if (identical(toupper(method_in[1L]), "MML")) {
    rank_check <- .mfrm_population_location_rank(enc)
    if (!is.null(rank_check) && rank_check$rank < rank_check$columns) {
      stop(sprintf(paste(
        "The MML population and facet location design is rank-deficient (rank %d of %d).",
        "Population location and facet effects cannot be separated under the current constraints.",
        "Review the population formula, rating assignments, and identifying anchors."
      ), rank_check$rank, rank_check$columns), call. = FALSE)
    }
  }
  enc
}

.mfrm_emit_connectivity_warning <- function(conn_info, warn_connectivity = TRUE) {
  if (!isTRUE(warn_connectivity)) return(invisible(NULL))
  if (!isTRUE(conn_info$connected) && isTRUE(conn_info$n_subsets > 1L)) {
    warning(
      sprintf(
        "Disconnected rating design detected (%d disjoint subsets across Person and active facets). Relative locations across subsets are not jointly identified without common links or anchors; inspect check_mfrm_connectivity().",
        as.integer(conn_info$n_subsets)
      ),
      call. = FALSE
    )
  }
  if (!isTRUE(conn_info$interior_categories_supported)) {
    warning(
      "One or more interior rating categories have zero observed counts within their declared support; step thresholds may be weakly or non-finitely identified without step_constraints or category_support.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.mfrm_compute_xi_layout <- function(enc) {
  F <- length(enc$facets)
  K <- enc$n_cat - 1L
  f_start <- integer(F)
  f_len <- integer(F)
  f_names <- vector("list", F)
  cur <- 0L
  for (f in seq_len(F)) {
    f_start[f] <- cur
    fnm <- enc$facets[f]
    if (enc$facet_dummy[f] == 1L) {
      f_len[f] <- 0L
      f_names[[f]] <- character(0)
      next
    }
    if (enc$has_custom_facet_map[f] == 1L) {
      qf <- enc$custom_facet_q[f]
      f_len[f] <- qf
      if (!is.null(enc$explanatory_meta[[fnm]])) {
        f_names[[f]] <- paste0(fnm, ":", enc$explanatory_meta[[fnm]]$colnames)
      } else if (!is.null(enc$interaction_meta[[fnm]])) {
        f_names[[f]] <- paste0(fnm, "[c", seq_len(qf), "]")
      } else {
        f_names[[f]] <- paste0(fnm, "[", seq_len(qf), "]")
      }
      cur <- cur + qf
      next
    }
    unanc <- which(enc$anchor_fixed[[f]] == 0L)
    if (length(unanc) == 0L) {
      qf <- 0L
      nms <- character(0)
    } else if (enc$facet_noncenter[f] == 1L) {
      qf <- length(unanc)
      nms <- paste0(fnm, ":", enc$facet_levels[[f]][unanc])
    } else {
      qf <- max(0L, length(unanc) - 1L)
      nms <- if (qf > 0L) paste0(fnm, ":", enc$facet_levels[[f]][unanc[seq_len(qf)]]) else character(0)
    }
    f_len[f] <- qf
    f_names[[f]] <- nms
    cur <- cur + qf
  }
  step_start <- cur
  if (isTRUE(enc$has_custom_step_map)) {
    qs <- enc$custom_step_q
    s_names <- if (qs > 0L) paste0("StepBasis[", seq_len(qs), "]") else character(0)
  } else {
    qs <- 0L
    s_names <- character(0)
    for (g in seq_len(enc$n_step_groups)) {
      ag <- enc$step_cat_min[g]
      bg <- enc$step_cat_max[g]
      Kg <- bg - ag
      if (Kg >= 2L) {
        qs <- qs + (Kg - 1L)
        s_names <- c(s_names, paste0(enc$step_owners[g], ":Step_", seq_len(Kg - 1L)))
      }
    }
  }
  cur <- cur + qs

  slope_start <- cur
  q_sl_tot <- 0L
  sl_names <- character(0)
  sl_offset <- integer(enc$n_slope_facets %||% 0L)
  if (isTRUE((enc$n_slope_facets %||% 0L) > 0L)) {
    for (s in seq_len(enc$n_slope_facets)) {
      sl_offset[s] <- cur
      q_s <- enc$slope_q[s]
      fnm_s <- enc$slope_facets[s]
      f_idx_s <- enc$slope_facets0[s] + 1L
      lvs_s <- enc$facet_levels[[f_idx_s]]
      if (q_s > 0L) {
        sl_names <- c(sl_names, paste0("LogSlope:", fnm_s, ":", lvs_s[seq_len(q_s)]))
      }
      cur <- cur + q_s
      q_sl_tot <- q_sl_tot + q_s
    }
  }
  list(
    facet_start = f_start,
    facet_len = f_len,
    facet_xi_names = f_names,
    step_start = step_start,
    step_len = qs,
    step_xi_names = s_names,
    slope_start = slope_start,
    slope_len = q_sl_tot,
    slope_offset = sl_offset,
    slope_xi_names = sl_names,
    n_xi = cur,
    xi_names = c(unlist(f_names, use.names = FALSE), s_names, sl_names)
  )
}

.mfrm_format_result <- function(enc, raw, method_label, person_score, xtreme, ctrl_list, elapsed_sec = NA_real_) {
    use_wle_primary <- (identical(method_label, "MML") && identical(toupper(person_score), "WLE"))
    p_est <- if (use_wle_primary)
        raw$person_wle
    else raw$person_est
    p_se <- if (use_wle_primary)
        raw$person_wle_se
    else raw$person_se
    p_infit <- raw$person_diag$infit_mnsq
    p_rse <- if (use_wle_primary) {
        ifelse(is.finite(p_se) & is.finite(p_infit), p_se * sqrt(pmax(1, p_infit)), p_se)
    }
    else {
        raw$person_real_se
    }
    n_base <- length(enc$base_facets)
    p_fair_m <- raw$person_diag$fair_m
    score_lbl <- switch(method_label, JML = if (xtreme > 0) "JML_xtreme" else "JML", JML_jeffreys = "JML_jeffreys",
        MML = if (use_wle_primary) "WLE" else "EAP", CML = "WLE", method_label)
    person_df <- data.frame(Person = enc$person_levels, Estimate = p_est, ModelSE = p_se, RealSE = p_rse,
        ConditionalSE = raw$person_cond_se, EAP = if (length(raw$person_eap) == enc$n_person)
            raw$person_eap
        else rep(NA_real_, enc$n_person), PosteriorSD = if (length(raw$person_psd) == enc$n_person)
            raw$person_psd
        else rep(NA_real_, enc$n_person), WLE = if (length(raw$person_wle) == enc$n_person)
            raw$person_wle
        else rep(NA_real_, enc$n_person), WLE_SE = if (length(raw$person_wle_se) == enc$n_person)
            raw$person_wle_se
        else rep(NA_real_, enc$n_person), FairM = p_fair_m, FairZ = raw$person_diag$fair_z, ObsCount = raw$person_diag$obs_n,
        ObsAvg = ifelse(raw$person_diag$obs_n > 0, raw$person_diag$obs_score/raw$person_diag$obs_n, NA_real_),
        ExpAvg = ifelse(raw$person_diag$obs_n > 0, raw$person_diag$exp_score/raw$person_diag$obs_n, NA_real_),
        InfitMnSq = p_infit, InfitZSTD = raw$person_diag$infit_zstd, OutfitMnSq = raw$person_diag$outfit_mnsq,
        OutfitZSTD = raw$person_diag$outfit_zstd, ExtremeFlag = raw$person_extreme_flag, Subset = raw$person_subset,
        ScoreMethod = score_lbl, stringsAsFactors = FALSE)
    facet_rows <- vector("list", n_base)
    cov_model <- vector("list", length(enc$facets))
    cov_sandwich <- vector("list", length(enc$facets))
    names(cov_model) <- names(cov_sandwich) <- enc$facets
    for (f in seq_along(enc$facets)) {
        fnm <- enc$facets[f]
        lvs <- enc$facet_levels[[f]]
        cm <- raw$facet_cov_model[[f]]
        cs <- raw$facet_cov_sandwich[[f]]
        if (nrow(cm) == length(lvs) && ncol(cm) == length(lvs))
            dimnames(cm) <- list(lvs, lvs)
        if (nrow(cs) == length(lvs) && ncol(cs) == length(lvs))
            dimnames(cs) <- list(lvs, lvs)
        cov_model[[f]] <- cm
        cov_sandwich[[f]] <- cs
        if (f <= n_base) {
            fd <- raw$facet_diag[[f]]
            facet_rows[[f]] <- data.frame(Facet = fnm, Level = lvs, Sign = enc$facet_signs[f], Estimate = raw$facet_est[[f]],
                ModelSE = raw$facet_se[[f]], RealSE = raw$facet_real_se[[f]], SandwichSE = raw$facet_joint_se[[f]],
                FairM = fd$fair_m, FairZ = fd$fair_z, ObsCount = fd$obs_n, ObsAvg = ifelse(fd$obs_n > 0, fd$obs_score/fd$obs_n,
                  NA_real_), ExpAvg = ifelse(fd$obs_n > 0, fd$exp_score/fd$obs_n, NA_real_), InfitMnSq = fd$infit_mnsq,
                InfitZSTD = fd$infit_zstd, OutfitMnSq = fd$outfit_mnsq, OutfitZSTD = fd$outfit_zstd, Subset = raw$facet_subset[[f]],
                stringsAsFactors = FALSE)
        }
    }
    facets_df <- do.call(rbind, facet_rows)
    rownames(facets_df) <- NULL
    interactions_df <- NULL
    if (length(enc$interaction_meta) > 0L) {
        int_rows <- lapply(enc$interaction_meta, function(im) {
            f <- im$facet_idx
            est <- raw$facet_est[[f]]
            mse <- raw$facet_se[[f]]
            sse <- raw$facet_joint_se[[f]]
            fd <- raw$facet_diag[[f]]
            zval <- ifelse(is.finite(mse) & mse > 0, est/mse, NA_real_)
            data.frame(Interaction = im$name, FacetA = im$facet_a, LevelA = im$levels_a, FacetB = im$facet_b,
                LevelB = im$levels_b, Estimate = est, ModelSE = mse, SandwichSE = sse, ObsCount = fd$obs_n,
                z = zval, pValue = 2 * stats::pnorm(-abs(zval)), stringsAsFactors = FALSE)
        })
        interactions_df <- do.call(rbind, int_rows)
        rownames(interactions_df) <- NULL
    }
    layout <- .mfrm_compute_xi_layout(enc)
    xi_est <- as.numeric(raw$xi_est)
    xi_cov_mod <- raw$xi_cov_model
    xi_cov_sw <- raw$xi_cov_sandwich
    full_xi_names <- layout$xi_names
    if (identical(method_label, "MML") && enc$pop_mode == 1L) {
        full_xi_names <- c(full_xi_names, paste0("Pop:", enc$pop_col_names), "Pop:log_sd")
    }
    if (length(xi_est) == length(full_xi_names)) {
        names(xi_est) <- full_xi_names
        if (nrow(xi_cov_mod) == length(full_xi_names)) {
            dimnames(xi_cov_mod) <- list(full_xi_names, full_xi_names)
        }
        if (nrow(xi_cov_sw) == length(full_xi_names)) {
            dimnames(xi_cov_sw) <- list(full_xi_names, full_xi_names)
        }
    }
    explanatory_df <- NULL
    if (length(enc$explanatory_meta) > 0L) {
        exp_rows <- lapply(enc$explanatory_meta, function(em) {
            f <- em$facet_idx
            idx <- layout$facet_start[f] + seq_len(layout$facet_len[f])
            est <- unname(xi_est[idx])
            mse <- if (nrow(xi_cov_mod) >= max(idx))
                sqrt(pmax(0, diag(xi_cov_mod)[idx]))
            else rep(NA_real_, length(idx))
            sse <- if (nrow(xi_cov_sw) >= max(idx))
                sqrt(pmax(0, diag(xi_cov_sw)[idx]))
            else rep(NA_real_, length(idx))
            zval <- ifelse(is.finite(mse) & mse > 0, est/mse, NA_real_)
            data.frame(Facet = em$facet, Term = em$colnames, Estimate = est, ModelSE = mse, SandwichSE = sse,
                z = zval, pValue = 2 * stats::pnorm(-abs(zval)), stringsAsFactors = FALSE)
        })
        explanatory_df <- do.call(rbind, exp_rows)
        rownames(explanatory_df) <- NULL
    }
    K <- enc$n_cat - 1L
    n_sg <- nrow(raw$step_est)
    sg_labels <- enc$step_owners
    step_sw_se <- if (nrow(raw$step_cov_sandwich) == n_sg * K) {
        matrix(sqrt(pmax(0, diag(raw$step_cov_sandwich))), nrow = n_sg, ncol = K, byrow = TRUE)
    }
    else {
        matrix(NA_real_, nrow = n_sg, ncol = K)
    }
    step_list <- vector("list", n_sg)
    for (g in seq_len(n_sg)) {
        step_list[[g]] <- data.frame(StepFacet = if (enc$model_code == 1L)
            enc$step_facet
        else "Common", ScaleLevel = sg_labels[g], Step = paste0("Step_", seq_len(K)), StepIndex = seq_len(K),
            FromCategory = enc$rating_min + seq(0L, K - 1L), ToCategory = enc$rating_min + seq_len(K), Estimate = as.numeric(raw$step_est[g,
                ]), ModelSE = as.numeric(raw$step_se[g, ]), SandwichSE = as.numeric(step_sw_se[g, ]), Supported = (seq(0L,
                K - 1L) >= enc$step_cat_min[g]) & (seq(0L, K - 1L) < enc$step_cat_max[g]), stringsAsFactors = FALSE)
    }
    steps_df <- do.call(rbind, step_list)
    rownames(steps_df) <- NULL
    step_row_nms <- paste0(steps_df$ScaleLevel, ":", steps_df$Step)
    step_cov_mod <- raw$step_cov_model
    step_cov_sw <- raw$step_cov_sandwich
    if (nrow(step_cov_mod) == length(step_row_nms))
        dimnames(step_cov_mod) <- list(step_row_nms, step_row_nms)
    if (nrow(step_cov_sw) == length(step_row_nms))
        dimnames(step_cov_sw) <- list(step_row_nms, step_row_nms)
    slopes_df <- NULL
    slope_cov_mod <- list()
    slope_cov_sw <- list()
    if (isTRUE((enc$n_slope_facets %||% 0L) > 0L) && length(raw$slope_est) == enc$n_slope_facets) {
        sl_rows <- vector("list", enc$n_slope_facets)
        for (s in seq_len(enc$n_slope_facets)) {
            fnm_s <- enc$slope_facets[s]
            f_idx_s <- enc$slope_facets0[s] + 1L
            lvs_s <- enc$facet_levels[[f_idx_s]]
            fd_s <- raw$facet_diag[[f_idx_s]]
            cm_s <- raw$slope_cov_model[[s]]
            cs_s <- raw$slope_cov_sandwich[[s]]
            if (nrow(cm_s) == length(lvs_s))
                dimnames(cm_s) <- list(lvs_s, lvs_s)
            if (nrow(cs_s) == length(lvs_s))
                dimnames(cs_s) <- list(lvs_s, lvs_s)
            slope_cov_mod[[fnm_s]] <- cm_s
            slope_cov_sw[[fnm_s]] <- cs_s
            sl_rows[[s]] <- data.frame(Facet = fnm_s, SlopeFacet = lvs_s, Level = lvs_s, LogEstimate = as.numeric(raw$log_slope_est[[s]]),
                Estimate = as.numeric(raw$slope_est[[s]]), ModelSE = as.numeric(raw$slope_se[[s]]), SandwichSE = as.numeric(raw$slope_joint_se[[s]]),
                ObsCount = as.numeric(fd_s$obs_n), Formula = enc$slope_type, stringsAsFactors = FALSE)
        }
        slopes_df <- do.call(rbind, sl_rows)
        rownames(slopes_df) <- NULL
    }
    bias_df <- NULL
    if (length(raw$local_bias$level_a) > 0L) {
        bias_df <- data.frame(LevelA_Index = raw$local_bias$level_a + 1L, LevelB_Index = raw$local_bias$level_b +
            1L, ObsCount = raw$local_bias$obs_count, ObsScore = raw$local_bias$obs_score, ExpScore = raw$local_bias$exp_score,
            BiasEstimate = raw$local_bias$bias_est, BiasSE = raw$local_bias$bias_se, tStat = raw$local_bias$t_stat,
            df = raw$local_bias$df, pValue = 2 * stats::pt(-abs(raw$local_bias$t_stat), df = pmax(1, raw$local_bias$df)),
            stringsAsFactors = FALSE)
    }
    pop_se_mod <- as.numeric(raw$pop_se_model)
    pop_se_sw <- as.numeric(raw$pop_se_sandwich)
    pop_tbl <- NULL
    if (enc$pop_mode == 1L) {
        dX <- length(enc$pop_col_names)
        p_terms <- c(enc$pop_col_names, "SD")
        p_vals <- c(as.numeric(raw$pop_beta), raw$pop_sigma)
        se_m <- if (length(pop_se_mod) == dX + 1L)
            c(pop_se_mod[seq_len(dX)], raw$pop_sigma * pop_se_mod[dX + 1L])
        else rep(NA_real_, dX + 1L)
        se_s <- if (length(pop_se_sw) == dX + 1L)
            c(pop_se_sw[seq_len(dX)], raw$pop_sigma * pop_se_sw[dX + 1L])
        else rep(NA_real_, dX + 1L)
        pop_tbl <- data.frame(Term = p_terms, Estimate = p_vals, ModelSE = se_m, SandwichSE = se_s, stringsAsFactors = FALSE)
    }
    pop_list <- list(mode = if (enc$pop_mode == 1L) "estimated" else "fixed", sigma = raw$pop_sigma, variance = raw$pop_sigma^2,
        coefficients = stats::setNames(as.numeric(raw$pop_beta), enc$pop_col_names), table = pop_tbl, eap_reliability = {
            keep <- is.finite(raw$person_eap) & is.finite(raw$person_psd)
            v <- stats::var(raw$person_eap[keep])
            s <- mean(raw$person_psd[keep]^2)
            if (is.finite(v + s) && v + s > 0) v/(v + s) else NA_real_
        })
    n_fp <- as.integer(raw$n_free_par)
    ll_val <- as.numeric(raw$loglik)
    aic_val <- if (is.finite(ll_val))
        -2 * ll_val + 2 * n_fp
    else NA_real_
    bic_n <- if (identical(method_label, "MML"))
        enc$n_person
    else enc$n_obs
    bic_val <- if (is.finite(ll_val) && bic_n > 0L)
        -2 * ll_val + log(bic_n) * n_fp
    else NA_real_
    summary_tbl <- data.frame(Model = enc$model, Method = method_label, N = as.integer(enc$n_obs), Persons = as.integer(enc$n_person),
        FreeParams = n_fp, LogLik = ll_val, AIC = aic_val, BIC = bic_val, Converged = isTRUE(raw$converged),
        Iterations = as.integer(raw$iterations), stringsAsFactors = FALSE)
    out <- structure(list(method = method_label, model = enc$model, step_facet = enc$step_facet, score_method = score_lbl,
        converged = isTRUE(raw$converged), iterations = as.integer(raw$iterations), loglik = ll_val, aic = aic_val,
        bic = bic_val, objective = as.numeric(raw$objective), max_abs_grad = as.numeric(raw$max_abs_grad),
        n_obs = as.integer(enc$n_obs), n_person = as.integer(enc$n_person), n_free_par = n_fp, rating_min = enc$rating_min,
        rating_max = enc$rating_max, connectivity = list(connected = isTRUE(raw$connected), n_subsets = as.integer(raw$n_subsets),
            interior_categories_supported = isTRUE(raw$interior_categories_supported)), summary = summary_tbl,
        person = person_df, facets = facets_df, steps = steps_df, slopes = slopes_df, explanatory = explanatory_df,
        interactions = interactions_df, step_cumulative = raw$step_cum, alpha_obs = if (length(raw$alpha_obs) ==
            enc$n_obs) as.numeric(raw$alpha_obs) else rep(1, enc$n_obs), population = pop_list, covariance_model = cov_model,
        covariance_sandwich = cov_sandwich, step_covariance_model = step_cov_mod, step_covariance_sandwich = step_cov_sw,
        slope_covariance_model = slope_cov_mod, slope_covariance_sandwich = slope_cov_sw, xi = xi_est, xi_covariance_model = xi_cov_mod,
        xi_covariance_sandwich = xi_cov_sw, xi_layout = layout, local_bias = bias_df, encoding = enc, control = ctrl_list,
        elapsed_sec = elapsed_sec), class = "mfrm_fit")
    out$schema_version <- 1L
    out$boundary_audit <- .mfrm_structural_boundary(out)
    if (identical(method_label, "MML")) {
        out["population_location_rank"] <- list(.mfrm_population_location_rank(enc))
    }
    out$conditional_location_rank <- .mfrm_conditional_location_rank(out)
    out$readiness <- .mfrm_readiness(out)
    out$summary$InferenceReady <- out$readiness$ready
    out$summary$InferenceStatus <- out$readiness$status
    if (isTRUE(out$boundary_audit$detected))
        warning("No finite structural MML estimate: a likelihood-increasing boundary direction was detected. Inference and classification require review.",
            call. = FALSE)
    out
}























.mfrm_rater_linkage <- function(enc, facet, network = FALSE, rater_panels = NULL) {
  j <- match(facet, enc$base_facets)
  if (!length(j) || is.na(j)) return(NULL)
  edges <- unique(data.frame(Person = enc$person_levels[enc$person0 + 1L],
                             Rater = enc$facet_levels[[j]][enc$facet_mat0[, j] + 1L]))
  audit <- mfrm_bipartite_components(edges, "Rater", "observed",
                                     include_membership = !is.null(rater_panels))
  n_components <- audit$summary$Components
  bridge_persons <- character()
  if (n_components == 1L) {
    assigned <- split(match(edges$Rater, enc$facet_levels[[j]]), edges$Person)
    key <- vapply(assigned, function(x) paste(sort(x), collapse = ","), character(1))
    # Another person with the same assignment supplies every link in that assignment.
    candidates <- names(assigned)[lengths(assigned) > 1L & !duplicated(key) & !duplicated(key, fromLast = TRUE)]
    if (length(candidates)) {
      stable <- edges[!edges$Person %in% candidates, , drop = FALSE]
      stable_connected <- length(unique(stable$Rater)) == length(unique(edges$Rater)) &&
        mfrm_bipartite_components(stable, "Rater", "observed")$summary$Connected
      if (!stable_connected) {
        bridge_persons <- candidates[vapply(candidates, function(p) {
          without <- edges[edges$Person != p, , drop = FALSE]
          nrow(without) > 0L && mfrm_bipartite_components(without, "Rater", "observed")$summary$Components > 1L
        }, logical(1))]
      }
    }
  }
  audit$components$Facet <- facet
  robustness <- if (isTRUE(network) || !is.null(rater_panels))
    .mfrm_network_robustness(edges, audit, rater_panels)
  list(facet = facet, n_components = n_components, components = audit$components,
       bridge_persons = bridge_persons, robustness = robustness)
}

check_mfrm_connectivity <- function(data,
                                    person,
                                    facets,
                                    score,
                                    model = c("RSM", "PCM", "GPCM"),
                                    step_facet = NULL,
                                    rating_min = NULL,
                                    rating_max = NULL,
                                    anchors = NULL,
                                    group_anchors = NULL,
                                    dummy_facets = NULL,
                                    category_support = NULL,
                                    facet_formulas = NULL,
                                    facet_data = NULL,
                                    rater_facet = NULL, network = FALSE, rater_panels = NULL) {
  model <- .mfrm_model_code(model)
  enc <- .mfrm_encode_dataset(
    data = data, person = person, facets = facets, score = score,
    model = model, step_facet = step_facet,
    rating_min = rating_min, rating_max = rating_max,
    anchors = anchors, group_anchors = group_anchors, dummy_facets = dummy_facets,
    category_support = category_support,
    facet_formulas = facet_formulas, facet_data = facet_data
  )
  if (!is.null(rater_facet) && (length(rater_facet) != 1L || !rater_facet %in% enc$base_facets)) {
    stop("rater_facet must name one of the declared facets.", call. = FALSE)
  }
  rater_facet <- rater_facet %||% infer_default_rater_facet(enc$base_facets, fallback_first = FALSE)
  if (is.null(rater_facet) && !is.null(rater_panels))
    stop("Specify rater_facet when supplying rater_panels.", call. = FALSE)
  c(.mfrm_connectivity_summary(enc), list(
    rater_linkage = .mfrm_rater_linkage(enc, rater_facet, network, rater_panels)
  ))
}

.mfrm_connectivity_summary <- function(enc) {
  raw <- mfrm_cpp_fast_connectivity_spec(enc)
  f_sub <- lapply(seq_along(enc$base_facets), function(f) {
    data.frame(
      Facet = enc$base_facets[f],
      Level = enc$facet_levels[[f]],
      Subset = raw$facet_subset[[f]],
      stringsAsFactors = FALSE
    )
  })
  list(
    connected = isTRUE(raw$connected),
    n_subsets = as.integer(raw$n_subsets),
    interior_categories_supported = isTRUE(raw$interior_categories_supported),
    n_obs = enc$n_obs,
    n_person = enc$n_person,
    n_extreme_min = sum(raw$person_extreme_flag == -1L),
    n_extreme_max = sum(raw$person_extreme_flag == 1L),
    person_summary = data.frame(
      Person = enc$person_levels,
      ObsCount = raw$person_obs_count,
      RawScore = raw$person_raw_score + raw$person_obs_count * enc$rating_min,
      MinScore = raw$person_min_score + raw$person_obs_count * enc$rating_min,
      MaxScore = raw$person_max_score + raw$person_obs_count * enc$rating_min,
      ExtremeFlag = raw$person_extreme_flag,
      Subset = raw$person_subset,
      stringsAsFactors = FALSE
    ),
    facet_subsets = do.call(rbind, f_sub)
  )
}

# Self-contained directed reachability & LP recession-cone boundary check for MML + population_formula
.mfrm_has_strict_edge <- function(from, to, n) {
  edges <- unique(data.frame(from = as.integer(from), to = as.integer(to)))
  if (nrow(edges) == 0L) return(FALSE)
  parent <- seq_len(n)
  find_r <- function(x) {
    while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
    x
  }
  for (r in seq_len(nrow(edges))) {
    ra <- find_r(edges$from[r]); rb <- find_r(edges$to[r])
    if (ra != rb) parent[rb] <- ra
  }
  active <- unique(c(edges$from, edges$to))
  roots <- unique(vapply(active, find_r, integer(1)))
  reaches_all <- function(src, dst) {
    adj <- split(dst, factor(src, levels = seq_len(n)))
    reached <- logical(n)
    q <- integer(n)
    for (rt in roots) {
      hd <- 1L; tl <- 1L; q[1L] <- rt; reached[rt] <- TRUE
      while (hd <= tl) {
        nxt <- adj[[q[hd]]]
        nxt <- nxt[!reached[nxt]]
        if (length(nxt) > 0L) {
          reached[nxt] <- TRUE
          q[tl + seq_along(nxt)] <- nxt
          tl <- tl + length(nxt)
        }
        hd <- hd + 1L
      }
    }
    all(reached[active])
  }
  !reaches_all(edges$from, edges$to) || !reaches_all(edges$to, edges$from)
}

.mfrm_facet_basis <- function(enc, f) {
  Lf <- enc$facet_n_levels[f]
  if (enc$facet_dummy[f] == 1L) return(matrix(0.0, Lf, 0L))
  if (enc$has_custom_facet_map[f] == 1L) return(enc$custom_facet_B[[f]])
  unanc <- which(enc$anchor_fixed[[f]] == 0L)
  if (length(unanc) == 0L) return(matrix(0.0, Lf, 0L))
  if (enc$facet_noncenter[f] == 1L) {
    B <- matrix(0.0, Lf, length(unanc))
    B[cbind(unanc, seq_along(unanc))] <- 1.0
    return(B)
  }
  qf <- length(unanc) - 1L
  B <- matrix(0.0, Lf, qf)
  if (qf > 0L) {
    B[cbind(unanc[seq_len(qf)], seq_len(qf))] <- 1.0
    B[unanc[length(unanc)], ] <- -1.0
  }
  B
}

.mfrm_population_location_rank <- function(enc) {
  # In RSM/PCM, a change in the population mean adds to the same predictor as
  # a facet location change while steps and population variance stay fixed.
  if (!identical(enc$pop_mode, 1L) || !enc$model %in% c("RSM", "PCM") ||
      !ncol(enc$pop_X_mat)) return(NULL)

  population <- enc$pop_X_mat[enc$person0 + 1L, , drop = FALSE]
  rows <- !duplicated(cbind(enc$facet_mat0, population))
  levels <- enc$facet_mat0[rows, , drop = FALSE]
  facets <- lapply(seq_along(enc$facets), function(f) {
    basis <- .mfrm_facet_basis(enc, f)
    enc$facet_signs[f] * basis[levels[, f] + 1L, , drop = FALSE]
  })
  design <- do.call(cbind, c(list(population[rows, , drop = FALSE]), facets))
  scale <- apply(abs(design), 2L, max)
  design <- sweep(design, 2L, ifelse(scale > 0, scale, 1), "/")
  list(rank = qr(design, tol = 1e-10)$rank, columns = ncol(design))
}

.mfrm_check_mml_population_boundary <- function(enc) {
  if (!enc$model %in% c("RSM", "PCM") || enc$pop_mode != 1L ||
      any(enc$facet_noncenter == 1L) || length(enc$interaction_meta) > 0L) {
    return(invisible(NULL))
  }
  design <- enc$pop_X_mat
  if (!is.matrix(design) || ncol(design) == 0L ||
      (ncol(design) == 1L && length(unique(design[, 1L])) == 1L)) {
    return(invisible(NULL))
  }
  obs_p <- sort(unique(enc$person0 + 1L))
  p_des <- design[obs_p, , drop = FALSE]
  if (any(!is.finite(p_des))) return(invisible(NULL))

  prof_key <- apply(p_des, 1L, paste, collapse = "\r")
  u_keys <- unique(prof_key)
  n_groups <- length(u_keys)
  prof_of_obs <- match(prof_key, u_keys)
  group_row <- prof_of_obs[match(enc$person0 + 1L, obs_p)]
  profile_design <- p_des[match(u_keys, prof_key), , drop = FALSE]

  minimum <- enc$score0 == 0L
  maximum <- enc$score0 == (enc$n_cat - 1L)
  if (!any(minimum | maximum)) return(invisible(NULL))

  stop_boundary <- function(facets) {
    stop(
      "No finite maximum-likelihood estimate for this MML fit. ",
      "The population and facet location design (", paste(facets, collapse = ", "),
      ") allows shifts that keep increasing the likelihood. ",
      "Review the shared ratings and category support before reporting adjusted scores.",
      call. = FALSE
    )
  }

  cell_key <- if (length(enc$base_facets) == 1L) {
    as.character(enc$facet_mat0[, 1L] + 1L)
  } else {
    apply(enc$facet_mat0[, seq_along(enc$base_facets), drop = FALSE], 1L, paste, collapse = "\r")
  }
  u_cells <- unique(cell_key)
  cell_id <- match(cell_key, u_cells)
  strict_edge <- function(cid) {
    lvl <- n_groups + cid
    .mfrm_has_strict_edge(c(group_row[!maximum], lvl[!minimum]),
                               c(lvl[!maximum], group_row[!minimum]),
                               n_groups + max(cid))
  }
  if (!strict_edge(cell_id)) return(invisible(NULL))

  if (n_groups == ncol(design) && qr(profile_design)$rank == n_groups) {
    for (f in seq_along(enc$base_facets)) {
      fnm <- enc$base_facets[f]
      if (enc$facet_n_levels[f] >= 2L && enc$facet_noncenter[f] == 0L &&
          all(enc$anchor_fixed[[f]] == 0L) && is.null(enc$group_anchor_meta[[fnm]]) &&
          is.null(enc$explanatory_meta[[fnm]]) &&
          strict_edge(enc$facet_mat0[, f] + 1L)) {
        stop_boundary(fnm)
      }
    }
  }

  gc_key <- paste(group_row, cell_id, sep = "\r")
  u_gc <- unique(gc_key)
  idx_gc <- match(gc_key, u_gc)
  rep_rows <- match(seq_along(u_gc), idx_gc)

  loc <- profile_design[group_row[rep_rows], , drop = FALSE]
  sc <- apply(abs(loc), 2L, max)
  loc <- sweep(loc, 2L, ifelse(sc > 0, sc, 1.0), "/")
  for (f in seq_along(enc$base_facets)) {
    Bf <- .mfrm_facet_basis(enc, f)
    if (ncol(Bf) > 0L) {
      lv <- enc$facet_mat0[rep_rows, f] + 1L
      loc <- cbind(loc, enc$facet_signs[f] * Bf[lv, , drop = FALSE])
    }
  }
  contrast <- rbind(
    -loc[unique(idx_gc[!maximum]), , drop = FALSE],
    loc[unique(idx_gc[!minimum]), , drop = FALSE]
  )
  m_col <- ncol(contrast)
  if (m_col == 0L || nrow(contrast) == 0L) return(invisible(NULL))
  target <- colSums(contrast)
  const_mat <- rbind(
    cbind(contrast, -contrast),
    rep(1.0, 2L * m_col)
  )
  const_dir <- c(rep(">=", nrow(contrast)), "<=")
  const_rhs <- c(rep(0.0, nrow(contrast)), 1.0)
  lp <- tryCatch(
    lpSolve::lp("max", objective.in = c(target, -target),
                const.mat = const_mat, const.dir = const_dir, const.rhs = const_rhs),
    error = function(e) NULL
  )
  if (is.null(lp) || lp$status != 0L || !is.finite(lp$objval)) {
    stop("Could not verify the population and facet location boundary numerically. Review the population design and anchors before reporting adjusted scores.", call. = FALSE)
  }
  if (lp$objval <= 1e-11) return(invisible(NULL))
  dir_vec <- lp$solution[seq_len(m_col)] - lp$solution[m_col + seq_len(m_col)]
  margins <- as.numeric(contrast %*% dir_vec)
  roundoff <- 32 * .Machine$double.eps * as.numeric(abs(contrast) %*% abs(dir_vec))
  if (all(is.finite(margins)) && all(margins >= -roundoff) && any(margins > roundoff)) {
    stop_boundary(enc$base_facets)
  }
  invisible(NULL)
}

.mfrm_check_design_rank <- function(enc, warn_connectivity = TRUE) {
  if (!isTRUE(warn_connectivity) || length(enc$base_facets) < 2L) return(invisible(NULL))
  f_mat <- enc$facet_mat0[, seq_along(enc$base_facets), drop = FALSE]
  u_rows <- !duplicated(f_mat)
  f_u <- f_mat[u_rows, , drop = FALSE]
  blocks <- list()
  for (f in seq_along(enc$base_facets)) {
    Bf <- .mfrm_facet_basis(enc, f)
    if (ncol(Bf) > 0L) {
      blocks[[length(blocks) + 1L]] <- enc$facet_signs[f] * Bf[f_u[, f] + 1L, , drop = FALSE]
    }
  }
  if (length(blocks) == 0L) return(invisible(NULL))
  X_f <- do.call(cbind, blocks)
  if (ncol(X_f) > 0L) {
    rnk <- qr(X_f)$rank
    if (rnk < ncol(X_f)) {
      warning(
        sprintf(
          "Additive facet design matrix is rank-deficient across (%s): rank %d < %d free facet parameters (nullity = %d). Check for confounded facets (e.g., nested/panel-mode-cohort collinearity) or supply identifying anchors.",
          paste(enc$base_facets, collapse = ", "), rnk, ncol(X_f), ncol(X_f) - rnk
        ),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

# Estimator used by fit_mfrm(); result layout is independent of estimation method.
.mfrm_estimate <- function(data,
                          person,
                          facets,
                          score,
                          model = c("RSM", "PCM", "GPCM"),
                          method = c("MML", "JML", "JMLE", "JML_jeffreys", "CML"),
                          step_facet = NULL,
                          slope_facet = NULL,
                          slope_reference = NULL,
                          slope_type = c("full_predictor", "ability_only"),
                          gpcm_mml_identification = c("free_population", "fixed_standard_normal", "unconstrained_slopes"),
                          person_score = c("auto", "EAP", "WLE"),
                          rating_min = NULL,
                          rating_max = NULL,
                          weight = NULL,
                          anchors = NULL,
                          group_anchors = NULL,
                          dummy_facets = NULL,
                          noncenter_facet = "Person",
                          positive_facets = NULL,
                          centering_scope = c("remainder", "facet"),
                          population_formula = NULL,
                          person_data = NULL,
                          person_id = NULL,
                          population_sd = 1.0,
                          step_constraints = NULL,
                          category_support = NULL,
                          facet_formulas = NULL,
                          facet_data = NULL,
                          facet_interactions = NULL,
                          facet_shrinkage = c("none", "empirical_bayes", "laplace"),
                          facet_prior_sd = NULL,
                          slope_prior_sd = NULL,
                          quad_points = 31L,
                          mml_integration = c("fixed", "adaptive", "auto"),
                          xtreme = 0.0,
                          jml_penalty = c("none", "jeffreys"),
                          whexact = FALSE,
                          bias_facets = NULL,
                          warn_connectivity = TRUE,
                          maxit = 400L,
                          reltol = 1e-10) {
  t0 <- proc.time()[[3L]]
  model <- .mfrm_model_code(model)
  method_in <- match.arg(method[1L], c("MML", "JML", "JMLE", "JML_jeffreys", "CML"))
  jml_penalty <- match.arg(jml_penalty)
  if (method_in %in% c("JML", "JMLE") && identical(jml_penalty, "jeffreys")) {
    method_in <- "JML_jeffreys"
  }
  if (identical(method_in, "JMLE")) method_in <- "JML"
  person_score <- match.arg(person_score)
  centering_scope <- match.arg(centering_scope)
  mml_integration <- match.arg(mml_integration)
  if (mml_integration == "auto" && method_in != "MML") stop("Automatic quadrature applies only to MML.", call. = FALSE)
  slope_type <- match.arg(slope_type)
  gpcm_mml_identification <- match.arg(gpcm_mml_identification)
  facet_shrinkage <- match.arg(facet_shrinkage)

  method_code <- switch(
    method_in,
    JML = 0L,
    JML_jeffreys = 1L,
    MML = 2L,
    CML = 3L
  )

  enc <- .mfrm_encode_dataset(
    data = data, person = person, facets = facets, score = score,
    model = model, step_facet = step_facet,
    slope_facet = slope_facet, slope_reference = slope_reference,
    slope_type = slope_type, gpcm_mml_identification = gpcm_mml_identification,
    method_in = method_in,
    rating_min = rating_min, rating_max = rating_max, weight = weight,
    anchors = anchors, group_anchors = group_anchors, dummy_facets = dummy_facets,
    noncenter_facet = noncenter_facet, positive_facets = positive_facets,
    centering_scope = centering_scope,
    population_formula = population_formula, person_data = person_data,
    person_id = person_id, population_sd = population_sd,
    step_constraints = step_constraints,
    category_support = category_support,
    facet_formulas = facet_formulas,
    facet_data = facet_data,
    facet_interactions = facet_interactions
  )

  # L3 Safety Gate 1: MML population + facet separation non-existence audit
  if (identical(method_in, "MML") && !is.null(population_formula)) {
    .mfrm_check_mml_population_boundary(enc)
  }

  # L3 Safety Gate 2: Connectivity, multi-facet collinearity rank, & unused endpoint warning
  conn_pre <- mfrm_cpp_fast_connectivity_spec(enc)
  if (isTRUE(warn_connectivity)) {
    .mfrm_emit_connectivity_warning(conn_pre, warn_connectivity = TRUE)
    if (isTRUE(conn_pre$connected)) {
      .mfrm_check_design_rank(enc, warn_connectivity = TRUE)
      linkage <- .mfrm_rater_linkage(enc, infer_default_rater_facet(enc$base_facets, fallback_first = FALSE))
      if (!is.null(linkage) && linkage$n_components > 1L) {
        warning("The Person-by-", linkage$facet, " graph has ", linkage$n_components,
                " components with no shared persons between them. Cross-component comparisons rely on calibration assumptions or anchors.",
                call. = FALSE)
      } else if (length(linkage$bridge_persons)) {
        warning("The Person-by-", linkage$facet, " graph loses connectivity if an individual linking person is removed: ",
                paste(linkage$bridge_persons, collapse = ", "), ".", call. = FALSE)
      }
    }
    if (identical(model, "GPCM") && is.null(category_support)) {
      obs_cats <- unique(enc$score0)
      if (!(0L %in% obs_cats) || !((enc$n_cat - 1L) %in% obs_cats)) {
        warning(
          "One or more extreme rubric categories have zero observed counts in GPCM; unpenalized step/location parameters at unused endpoints may diverge without `category_support` or `step_constraints`.",
          call. = FALSE
        )
      }
    }
  }

  bf_a <- -1L
  bf_b <- -1L
  if (!is.null(bias_facets)) {
    if (length(bias_facets) != 2L || any(!bias_facets %in% enc$base_facets)) {
      stop("`bias_facets` must name two facets included in `facets`.", call. = FALSE)
    }
    bf_a <- match(bias_facets[1L], enc$facets) - 1L
    bf_b <- match(bias_facets[2L], enc$facets) - 1L
  }

  f_l2 <- 0.0
  if (identical(facet_shrinkage, "laplace")) {
    psd_f <- if (is.null(facet_prior_sd)) 2.0 else as.numeric(facet_prior_sd)[1L]
    if (is.finite(psd_f) && psd_f > 0) f_l2 <- 1.0 / (psd_f * psd_f)
  }
  s_l2 <- 0.0
  if (!is.null(slope_prior_sd)) {
    psd_s <- as.numeric(slope_prior_sd)[1L]
    if (is.finite(psd_s) && psd_s > 0) s_l2 <- 1.0 / (psd_s * psd_s)
  }

  ctrl_list <- list(
    method_code = as.integer(method_code),
    maxit = as.integer(maxit),
    reltol = as.double(reltol),
    quad_points = as.integer(quad_points),
    adaptive_quad = identical(mml_integration, "adaptive"),
    xtreme = as.double(xtreme),
    whexact = isTRUE(whexact),
    compute_diagnostics = TRUE,
    compute_wle = TRUE,
    wle_primary = identical(person_score, "WLE"),
    bias_facet_a = as.integer(bf_a),
    bias_facet_b = as.integer(bf_b),
    facet_l2_pen = as.double(f_l2),
    slope_l2_pen = as.double(s_l2)
  )

  raw <- mfrm_cpp_fast_fit_spec(enc, ctrl_list)
  out <- .mfrm_format_result(enc, raw, method_in, person_score, xtreme, ctrl_list, proc.time()[[3L]] - t0)
  if (identical(mml_integration, "auto")) out <- .mfrm_refine_quadrature(out)
  if (!is.null(bias_facets) && !is.null(out$local_bias) && nrow(out$local_bias) > 0L) {
    out$local_bias$FacetA <- bias_facets[1L]
    out$local_bias$LevelA <- enc$facet_levels[[bf_a + 1L]][out$local_bias$LevelA_Index]
    out$local_bias$FacetB <- bias_facets[2L]
    out$local_bias$LevelB <- enc$facet_levels[[bf_b + 1L]][out$local_bias$LevelB_Index]
  }
  if (identical(facet_shrinkage, "empirical_bayes")) {
    out <- apply_empirical_bayes_shrinkage(out, facet_prior_sd = facet_prior_sd)
    out$facets$MLEEstimate <- out$facets$Estimate
    out$facets$Estimate <- out$facets$ShrunkEstimate
  }
  out
}

# Internal helper retained for WLE scoring on any fit object











































mfrm_batch_fit <- function(data,
                           person,
                           facets,
                           score = "score",
                           score_matrix = NULL,
                           model = c("RSM", "PCM", "GPCM"),
                           method = c("MML", "JML", "JML_jeffreys", "CML"),
                           step_facet = NULL,
                           slope_facet = NULL,
                           slope_reference = NULL,
                           slope_type = c("full_predictor", "ability_only"),
                           gpcm_mml_identification = c("free_population", "fixed_standard_normal", "unconstrained_slopes"),
                           person_score = c("auto", "EAP", "WLE"),
                           rating_min = NULL,
                           rating_max = NULL,
                           anchors = NULL,
                           group_anchors = NULL,
                           dummy_facets = NULL,
                           population_formula = NULL,
                           person_data = NULL,
                           step_constraints = NULL,
                           category_support = NULL,
                           facet_formulas = NULL,
                           facet_data = NULL,
                           slope_prior_sd = NULL,
                           quad_points = 31L,
                           mml_integration = c("fixed", "adaptive"),
                           xtreme = 0.0,
                           whexact = FALSE,
                           bias_facets = NULL,
                           maxit = 400L,
                           reltol = 1e-10,
                           threads = 1L,
                           return_fits = TRUE) {
  t0 <- proc.time()[[3L]]
  model <- .mfrm_model_code(model)
  method <- match.arg(method[1L], c("MML", "JML", "JML_jeffreys", "CML"))
  person_score <- match.arg(person_score)
  mml_integration <- match.arg(mml_integration)
  slope_type <- match.arg(slope_type)
  gpcm_mml_identification <- match.arg(gpcm_mml_identification)
  threads <- .mfrm_threads(threads)

  method_code <- switch(method, JML = 0L, JML_jeffreys = 1L, MML = 2L, CML = 3L)
  bf_a <- -1L
  bf_b <- -1L
  if (!is.null(bias_facets)) {
    bf_a <- match(bias_facets[1L], facets) - 1L
    bf_b <- match(bias_facets[2L], facets) - 1L
  }
  s_l2 <- 0.0
  if (!is.null(slope_prior_sd)) {
    psd_s <- as.numeric(slope_prior_sd)[1L]
    if (is.finite(psd_s) && psd_s > 0) s_l2 <- 1.0 / (psd_s * psd_s)
  }
  ctrl_list <- list(
    method_code = as.integer(method_code),
    maxit = as.integer(maxit),
    reltol = as.double(reltol),
    quad_points = as.integer(quad_points),
    adaptive_quad = identical(mml_integration, "adaptive"),
    xtreme = as.double(xtreme),
    whexact = isTRUE(whexact),
    compute_diagnostics = TRUE,
    compute_wle = TRUE,
    wle_primary = identical(person_score, "WLE"),
    bias_facet_a = as.integer(bf_a),
    bias_facet_b = as.integer(bf_b),
    slope_l2_pen = as.double(s_l2)
  )

  if (!is.null(score_matrix)) {
    if (!is.data.frame(data)) {
      stop("When `score_matrix` is supplied, `data` must be a single template data.frame.", call. = FALSE)
    }
    score_matrix <- as.matrix(score_matrix)
    if (nrow(score_matrix) != nrow(data)) {
      stop("`nrow(score_matrix)` must equal `nrow(data)`.", call. = FALSE)
    }
    if (ncol(score_matrix) == 0L || !is.numeric(score_matrix) ||
        any(!is.finite(score_matrix)) || any(abs(score_matrix - round(score_matrix)) > 1e-8)) {
      stop("`score_matrix` must contain at least one column of finite integer-valued ratings.", call. = FALSE)
    }
    tmp_df <- data
    rmin <- if (is.null(rating_min)) min(score_matrix) else as.numeric(rating_min)[1L]
    rmax <- if (is.null(rating_max)) max(score_matrix) else as.numeric(rating_max)[1L]
    tmp_df[[score]] <- score_matrix[, 1L]
    enc <- .mfrm_encode_dataset(
      data = tmp_df, person = person, facets = facets, score = score,
      model = model, step_facet = step_facet,
      slope_facet = slope_facet, slope_reference = slope_reference,
      slope_type = slope_type, gpcm_mml_identification = gpcm_mml_identification,
      method_in = method,
      rating_min = rmin, rating_max = rmax,
      anchors = anchors, group_anchors = group_anchors,
      dummy_facets = dummy_facets,
      population_formula = population_formula, person_data = person_data,
      step_constraints = step_constraints, category_support = category_support,
      facet_formulas = facet_formulas, facet_data = facet_data
    )
    if (enc$n_obs != nrow(data)) {
      stop("The shared design must have complete Person and facet identifiers. Use a list of rating tables for designs with omitted rows.", call. = FALSE)
    }
    score_offset <- round(score_matrix - enc$rating_min)
    if (any(score_offset < enc$step_cat_min[enc$step_group0 + 1L] |
            score_offset > enc$step_cat_max[enc$step_group0 + 1L])) {
      stop("Scores in `score_matrix` fall outside the rubric bounds or declared `category_support`.", call. = FALSE)
    }
    smat0 <- matrix(as.integer(score_offset), nrow = nrow(score_matrix), ncol = ncol(score_matrix))
    raw_list <- mfrm_cpp_fast_batch_shared_design(enc, smat0, ctrl_list, threads)
    elapsed <- proc.time()[[3L]] - t0
    if (!isTRUE(return_fits)) {
      return(data.frame(
        Replicate = seq_along(raw_list),
        Converged = vapply(raw_list, function(x) isTRUE(x$converged), logical(1)),
        Iterations = vapply(raw_list, function(x) as.integer(x$iterations), integer(1)),
        LogLik = vapply(raw_list, function(x) as.numeric(x$loglik), numeric(1)),
        MaxAbsGrad = vapply(raw_list, function(x) as.numeric(x$max_abs_grad), numeric(1))
      ))
    }
    fits <- lapply(seq_along(raw_list), function(i) {
      replicate_enc <- enc
      replicate_enc$score0 <- smat0[, i]
      replicate_enc$retained_data[[score]] <- score_matrix[, i]
      .mfrm_format_result(replicate_enc, raw_list[[i]], method, person_score,
                               xtreme, ctrl_list, elapsed / length(raw_list))
    })
    attr(fits, "elapsed_sec") <- elapsed
    return(fits)
  }

  if (!is.list(data) || is.data.frame(data)) {
    stop("`data` must be a list of data.frames (or a template data.frame with `score_matrix`).", call. = FALSE)
  }

  enc_list <- lapply(seq_along(data), function(i) {
    pdata_i <- if (is.list(person_data) && !is.data.frame(person_data)) person_data[[i]] else person_data
    .mfrm_encode_dataset(
      data = data[[i]], person = person, facets = facets, score = score,
      model = model, step_facet = step_facet,
      slope_facet = slope_facet, slope_reference = slope_reference,
      slope_type = slope_type, gpcm_mml_identification = gpcm_mml_identification,
      method_in = method,
      rating_min = rating_min, rating_max = rating_max,
      anchors = anchors, group_anchors = group_anchors,
      dummy_facets = dummy_facets,
      population_formula = population_formula, person_data = pdata_i,
      step_constraints = step_constraints, category_support = category_support,
      facet_formulas = facet_formulas, facet_data = facet_data
    )
  })

  raw_list <- mfrm_cpp_fast_batch_list(enc_list, ctrl_list, threads)
  elapsed <- proc.time()[[3L]] - t0

  if (!isTRUE(return_fits)) {
    return(data.frame(
      Replicate = seq_along(raw_list),
      Converged = vapply(raw_list, function(x) isTRUE(x$converged), logical(1)),
      Iterations = vapply(raw_list, function(x) as.integer(x$iterations), integer(1)),
      LogLik = vapply(raw_list, function(x) as.numeric(x$loglik), numeric(1)),
      MaxAbsGrad = vapply(raw_list, function(x) as.numeric(x$max_abs_grad), numeric(1))
    ))
  }

  fits <- lapply(seq_along(raw_list), function(i) {
    .mfrm_format_result(enc_list[[i]], raw_list[[i]], method, person_score, xtreme, ctrl_list, elapsed / length(raw_list))
  })
  attr(fits, "elapsed_sec") <- elapsed
  fits
}

# ==============================================================================
# Standard S3 Model Methods for `mfrm_fit` (and `newdata` Prediction)
# ==============================================================================

.mfrm_eval_row_predictions <- function(eta_vec, step_cum, step_group0, step_cat_min, step_cat_max, rating_min,
                                       alpha_vec = NULL, theta_vec = NULL, base_eta_vec = NULL, slope_formula = 0L) {
  out <- mfrm_cpp_row_predictions(
    as.double(eta_vec), step_cum, as.integer(step_group0), as.integer(step_cat_min),
    as.integer(step_cat_max), rating_min, as.double(alpha_vec),
    as.double(theta_vec), as.double(base_eta_vec), slope_formula == 1L
  )
  colnames(out$probs) <- as.character(rating_min + seq_len(ncol(step_cum)) - 1L)
  out$eta <- eta_vec
  out
}

.mfrm_wle_newdata <- function(newdata, enc, step_group0) {
  score_col <- intersect(c(enc$score_col, "Score", "score"), names(newdata))[1L]
  person_col <- intersect(c(enc$person_col, "Person", "person"), names(newdata))[1L]
  if (is.na(score_col) || is.na(person_col)) {
    stop("`newdata` must contain person and score columns when `type = 'wle'`.", call. = FALSE)
  }
  persons <- as.character(newdata[[person_col]])
  if (anyNA(persons) || any(!nzchar(trimws(persons)))) {
    stop("New-data WLE requires a non-missing Person ID on every row.", call. = FALSE)
  }
  scores <- newdata[[score_col]]
  if (is.factor(scores)) scores <- as.character(scores)
  if (!is.null(dim(scores)) || !(is.character(scores) || is.numeric(scores)) || is.complex(scores)) {
    stop("New-data WLE scores must be finite, integer-valued categories.", call. = FALSE)
  }
  scores <- suppressWarnings(as.numeric(scores))
  if (any(!is.finite(scores)) || any(abs(scores - round(scores)) > 1e-8)) {
    stop("New-data WLE scores must be finite, integer-valued categories.", call. = FALSE)
  }
  offset <- round(scores - enc$rating_min)
  group <- step_group0 + 1L
  invalid <- which(offset < enc$step_cat_min[group] | offset > enc$step_cat_max[group])
  if (length(invalid)) {
    stop(sprintf("New-data WLE score at row %d falls outside the fitted rubric or category support.",
                 invalid[1L]), call. = FALSE)
  }
  levels <- .mfrm_sort_levels(newdata[[person_col]])
  list(levels = levels, person0 = match(persons, levels) - 1L, score0 = as.integer(offset))
}


predict.mfrm_fit <- function(object, newdata = NULL, type = c("response", "probs", "prob", "link", "fair",
    "wle", "eap"), person_data = NULL, facet_data = NULL, se.fit = FALSE, ...) {
    .mfrm_validate_fit(object)
    type <- match.arg(type)
    if (!is.null(newdata) && "Theta" %in% names(newdata) && !"Theta" %in% object$encoding$base_facets) {
        if ("theta" %in% names(newdata) && !isTRUE(all.equal(newdata$theta, newdata$Theta)))
            stop("theta and Theta specify different abilities; supply one ability column.", call. = FALSE)
        if (!"theta" %in% names(newdata)) newdata$theta <- newdata$Theta
    }
    if (identical(type, "prob"))
        type <- "probs"
    enc <- object$encoding
    if (identical(type, "fair")) {
        return(stats::setNames(object$person$FairM, object$person$Person))
    }
    if (identical(type, "eap")) {
        if (!is.null(newdata))
            stop("EAP scoring of new ratings is unavailable. Omit `newdata` for fitted EAP scores.", call. = FALSE)
        return(data.frame(Person = object$person$Person, Estimate = object$person$EAP, SE = object$person$PosteriorSD,
            stringsAsFactors = FALSE))
    }
    if (identical(type, "wle") && is.null(newdata)) {
        return(data.frame(Person = object$person$Person, Estimate = object$person$WLE, SE = object$person$WLE_SE,
            stringsAsFactors = FALSE))
    }
    if (is.null(newdata)) {
        th_row <- object$person$Estimate[enc$person0 + 1L]
        base_eta <- numeric(enc$n_obs)
        for (f in seq_along(enc$facets)) {
            fnm <- enc$facets[f]
            if (f <= length(enc$base_facets)) {
                f_sub <- object$facets[object$facets$Facet == fnm, ]
                b_vec <- f_sub$Estimate[match(enc$facet_levels[[f]], f_sub$Level)]
            }
            else {
                i_sub <- object$interactions[object$interactions$Interaction == fnm, ]
                b_vec <- i_sub$Estimate
            }
            base_eta <- base_eta + enc$facet_signs[f] * b_vec[enc$facet_mat0[, f] + 1L]
        }
        eta_vec <- th_row + base_eta
        sg0 <- enc$step_group0
    }
    else {
        if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
            stop("`newdata` must be a non-empty data.frame.", call. = FALSE)
        }
        miss_f <- setdiff(enc$base_facets, names(newdata))
        if (length(miss_f) > 0L) {
            stop("`newdata` is missing facet column(s): ", paste(miss_f, collapse = ", "), ".", call. = FALSE)
        }
        N_new <- nrow(newdata)
        base_eta <- numeric(N_new)
        f_lidx_list <- vector("list", length(enc$base_facets))
        for (f in seq_along(enc$base_facets)) {
            fnm <- enc$base_facets[f]
            vals_chr <- as.character(newdata[[fnm]])
            f_sub <- object$facets[object$facets$Facet == fnm, ]
            m_old <- match(vals_chr, f_sub$Level)
            b_row <- f_sub$Estimate[m_old]
            if (anyNA(m_old)) {
                em <- enc$explanatory_meta[[fnm]]
                if (!is.null(em) && !is.null(facet_data)) {
                  fdf <- if (is.list(facet_data) && !is.data.frame(facet_data))
                    facet_data[[fnm]]
                  else facet_data
                  id_c <- intersect(c(fnm, "Level", "Element"), names(fdf))[1L]
                  unseen <- which(is.na(m_old))
                  m_fd <- match(vals_chr[unseen], .mfrm_unique_covariate_ids(fdf[[id_c]], paste0("`facet_data` for facet '",
                    fnm, "'")))
                  if (anyNA(m_fd)) {
                    stop("Unseen level(s) in `newdata[['", fnm, "']]` not found in `facet_data`: ", paste(unique(vals_chr[unseen][is.na(m_fd)]),
                      collapse = ", "), ".", call. = FALSE)
                  }
                  mm_new <- stats::model.matrix(em$terms, data = fdf[m_fd, , drop = FALSE])
                  if (em$centered) {
                    mm_new <- mm_new[, colnames(mm_new) != "(Intercept)", drop = FALSE]
                    mm_new <- sweep(mm_new, 2L, em$col_means, "-")
                  }
                  exp_sub <- object$explanatory[object$explanatory$Facet == fnm, ]
                  gam <- exp_sub$Estimate[match(colnames(mm_new), exp_sub$Term)]
                  b_row[unseen] <- as.numeric(mm_new %*% gam)
                }
                else {
                  stop("Unknown level(s) in `newdata[['", fnm, "']]`: ", paste(unique(vals_chr[is.na(m_old)]),
                    collapse = ", "), ".", call. = FALSE)
                }
            }
            base_eta <- base_eta + enc$facet_signs[f] * b_row
            f_lidx_list[[f]] <- m_old - 1L
        }
        if (length(enc$interaction_meta) > 0L) {
            for (im in enc$interaction_meta) {
                fa <- match(im$facet_a, enc$base_facets)
                fb <- match(im$facet_b, enc$base_facets)
                Lb <- enc$facet_n_levels[fb]
                cidx <- f_lidx_list[[fa]] * Lb + f_lidx_list[[fb]] + 1L
                i_sub <- object$interactions[object$interactions$Interaction == im$name, ]
                base_eta <- base_eta + enc$facet_signs[im$facet_idx] * i_sub$Estimate[cidx]
            }
        }
        sg0 <- if (enc$model_code == 1L) {
            sf_chr <- as.character(newdata[[enc$step_facet]])
            m_sg <- match(sf_chr, enc$step_owners)
            if (anyNA(m_sg))
                stop("Unknown step_facet level in `newdata`.", call. = FALSE)
            as.integer(m_sg - 1L)
        }
        else {
            integer(N_new)
        }
        sform <- enc$slope_formula %||% 0L
        alpha_new <- NULL
        if (!is.null(object$slopes) && nrow(object$slopes) > 0L) {
            alpha_new <- rep(1, N_new)
            for (fnm_s in enc$slope_facets) {
                sl_sub <- object$slopes[object$slopes$Facet == fnm_s, ]
                m_sl <- match(as.character(newdata[[fnm_s]]), sl_sub$Level)
                alpha_new <- alpha_new * ifelse(is.na(m_sl), 1, sl_sub$Estimate[m_sl])
            }
        }
        if (identical(type, "wle")) {
            scoring <- .mfrm_wle_newdata(newdata, enc, sg0)
            p_new_lv <- scoring$levels
            p0_new <- scoring$person0
            sc0_new <- scoring$score0
            if (is.null(alpha_new)) {
                ans <- mfrm_cpp_fast_wle(as.integer(length(p_new_lv)), as.integer(enc$n_cat), as.integer(p0_new),
                  sc0_new, rep(1, N_new), as.double(base_eta), object$step_cumulative, sg0, enc$step_cat_min,
                  enc$step_cat_max)
                return(data.frame(Person = p_new_lv, Estimate = ans$theta, SE = ans$se, stringsAsFactors = FALSE))
            }
            th_out <- numeric(length(p_new_lv))
            se_out <- rep(NA_real_, length(p_new_lv))
            for (pp in seq_along(p_new_lv)) {
                ii <- which(p0_new == (pp - 1L))
                th <- 0
                for (it_w in seq_len(30L)) {
                  S_th <- 0
                  I_th <- 0
                  J_th <- 0
                  dJ_th <- 0
                  for (r in ii) {
                    g <- sg0[r] + 1L
                    kmin <- enc$step_cat_min[g]
                    kmax <- enc$step_cat_max[g]
                    ks <- kmin:kmax
                    al <- alpha_new[r]
                    S_vec <- object$step_cumulative[g, ks + 1L]
                    u <- if (sform == 1L) {
                      (ks - kmin) * (al * th + base_eta[r]) - S_vec
                    }
                    else {
                      al * ((ks - kmin) * (th + base_eta[r]) - S_vec)
                    }
                    u <- u - max(u)
                    pr <- exp(u)
                    pr <- pr/sum(pr)
                    m1 <- sum(ks * pr)
                    d_k <- ks - m1
                    v_r <- sum((d_k^2) * pr)
                    mu3 <- sum((d_k^3) * pr)
                    mu4 <- sum((d_k^4) * pr)
                    S_th <- S_th + al * (sc0_new[r] - m1)
                    I_th <- I_th + (al^2) * v_r
                    J_th <- J_th + (al^3) * mu3
                    dJ_th <- dJ_th + (al^4) * (mu4 - 3 * v_r^2)
                  }
                  if (I_th < 1e-12)
                    break
                  w_sc <- S_th + J_th/(2 * I_th)
                  h_eff <- max(1e-06, I_th - (dJ_th * I_th - J_th^2)/(2 * I_th^2))
                  dth <- max(-2, min(2, w_sc/h_eff))
                  th <- max(-15, min(15, th + dth))
                  if (abs(dth) < 1e-09)
                    break
                }
                th_out[pp] <- th
                I_fin <- 0
                for (r in ii) {
                  g <- sg0[r] + 1L
                  kmin <- enc$step_cat_min[g]
                  kmax <- enc$step_cat_max[g]
                  ks <- kmin:kmax
                  al <- alpha_new[r]
                  S_vec <- object$step_cumulative[g, ks + 1L]
                  u <- if (sform == 1L)
                    (ks - kmin) * (al * th + base_eta[r]) - S_vec
                  else al * ((ks - kmin) * (th + base_eta[r]) - S_vec)
                  u <- u - max(u)
                  pr <- exp(u)
                  pr <- pr/sum(pr)
                  m1 <- sum(ks * pr)
                  I_fin <- I_fin + (al^2) * sum(((ks - m1)^2) * pr)
                }
                if (I_fin > 1e-12)
                  se_out[pp] <- 1/sqrt(I_fin)
            }
            return(data.frame(Person = p_new_lv, Estimate = th_out, SE = se_out, stringsAsFactors = FALSE))
        }
        if ("theta" %in% names(newdata)) {
            th_row <- as.numeric(newdata$theta)
        }
        else if (enc$person_col %in% names(newdata)) {
            p_chr <- as.character(newdata[[enc$person_col]])
            m_p <- match(p_chr, object$person$Person)
            th_row <- object$person$Estimate[m_p]
            if (anyNA(m_p)) {
                unseen_p <- which(is.na(m_p))
                if (enc$pop_mode == 1L && !is.null(person_data)) {
                  m_pd <- match(p_chr[unseen_p], .mfrm_unique_covariate_ids(person_data[[enc$person_id]],
                    "`person_data`"))
                  if (anyNA(m_pd)) {
                    th_row[unseen_p] <- 0
                  }
                  else {
                    mm_p <- stats::model.matrix(enc$pop_terms, data = person_data[m_pd, , drop = FALSE])
                    mm_p <- mm_p[, enc$pop_col_names, drop = FALSE]
                    th_row[unseen_p] <- as.numeric(mm_p %*% object$population$coefficients)
                  }
                }
                else if (enc$pop_mode == 1L && "(Intercept)" %in% names(object$population$coefficients)) {
                  th_row[unseen_p] <- unname(object$population$coefficients["(Intercept)"])
                }
                else {
                  th_row[unseen_p] <- 0
                }
            }
        }
        else {
            th_row <- rep(0, N_new)
        }
        eta_vec <- th_row + base_eta
    }
    sform <- enc$slope_formula %||% 0L
    alpha_row <- NULL
    if (!is.null(object$slopes) && nrow(object$slopes) > 0L) {
        if (is.null(newdata)) {
            alpha_row <- object$alpha_obs
        }
        else {
            alpha_row <- rep(1, nrow(newdata))
            for (fnm_s in enc$slope_facets) {
                sl_sub <- object$slopes[object$slopes$Facet == fnm_s, ]
                m_sl <- match(as.character(newdata[[fnm_s]]), sl_sub$Level)
                alpha_row <- alpha_row * ifelse(is.na(m_sl), 1, sl_sub$Estimate[m_sl])
            }
        }
    }
    ev <- .mfrm_eval_row_predictions(eta_vec = eta_vec, step_cum = object$step_cumulative, step_group0 = sg0,
        step_cat_min = enc$step_cat_min, step_cat_max = enc$step_cat_max, rating_min = enc$rating_min, alpha_vec = alpha_row,
        theta_vec = th_row, base_eta_vec = base_eta, slope_formula = sform)
    if (identical(type, "probs"))
        return(ev$probs)
    val <- if (identical(type, "link"))
        ev$eta
    else ev$expected
    if (!isTRUE(se.fit))
        return(val)
    G <- .mfrm_prediction_gradient(object, newdata, ev, th_row, base_eta, sg0, alpha_row, type, facet_data)
    V <- object$xi_covariance_model
    se_out <- if (!.mfrm_readiness(object)$ready || any(!is.finite(V)))
        rep(NA_real_, length(val))
    else sqrt(pmax(0, rowSums((G %*% V) * G)))
    list(fit = val, se.fit = se_out)
}


fitted.mfrm_fit <- function(object, type = c("response", "probs", "prob"), ...) {
    .mfrm_validate_fit(object)
    type <- match.arg(type)
    predict.mfrm_fit(object, newdata = NULL, type = type, ...)
}


residuals.mfrm_fit <- function(object, type = c("response", "pearson", "standardized"), ...) {
    .mfrm_validate_fit(object)
    type <- match.arg(type)
    enc <- object$encoding
    th_row <- object$person$Estimate[enc$person0 + 1L]
    base_eta <- numeric(enc$n_obs)
    for (f in seq_along(enc$facets)) {
        fnm <- enc$facets[f]
        if (f <= length(enc$base_facets)) {
            f_sub <- object$facets[object$facets$Facet == fnm, ]
            b_vec <- f_sub$Estimate[match(enc$facet_levels[[f]], f_sub$Level)]
        }
        else {
            i_sub <- object$interactions[object$interactions$Interaction == fnm, ]
            b_vec <- i_sub$Estimate
        }
        base_eta <- base_eta + enc$facet_signs[f] * b_vec[enc$facet_mat0[, f] + 1L]
    }
    ev <- .mfrm_eval_row_predictions(eta_vec = th_row + base_eta, step_cum = object$step_cumulative, step_group0 = enc$step_group0,
        step_cat_min = enc$step_cat_min, step_cat_max = enc$step_cat_max, rating_min = enc$rating_min, alpha_vec = if (!is.null(object$slopes))
            object$alpha_obs
        else NULL, theta_vec = th_row, base_eta_vec = base_eta, slope_formula = enc$slope_formula %||% 0L)
    y_obs <- enc$rating_min + enc$score0
    raw_res <- y_obs - ev$expected
    if (identical(type, "response"))
        return(raw_res)
    ifelse(is.finite(ev$variance) & ev$variance > 1e-12, raw_res/sqrt(ev$variance), NA_real_)
}


logLik.mfrm_fit <- function(object, ...) {
    .mfrm_validate_fit(object)
    structure(object$loglik, df = object$n_free_par, nobs = object$n_obs, class = "logLik")
}


AIC.mfrm_fit <- function(object, ..., k = 2) {
    .mfrm_validate_fit(object)
    dots <- list(...)
    if (length(dots) == 0L) {
        return(-2 * object$loglik + k * object$n_free_par)
    }
    objs <- c(list(object), dots)
    data.frame(df = vapply(objs, function(x) as.numeric(attr(stats::logLik(x), "df")), numeric(1)), AIC = vapply(objs,
        function(x) -2 * as.numeric(stats::logLik(x)) + k * as.numeric(attr(stats::logLik(x), "df")), numeric(1)))
}


BIC.mfrm_fit <- function(object, ...) {
    .mfrm_validate_fit(object)
    n_eff <- if (identical(object$method, "MML"))
        object$n_person
    else object$n_obs
    -2 * object$loglik + log(n_eff) * object$n_free_par
}


nobs.mfrm_fit <- function(object, ...) {
    .mfrm_validate_fit(object)
    object$n_obs
}


coef.mfrm_fit <- function(object, type = c("all", "facets", "locations", "steps", "slopes", "person", "persons",
    "explanatory", "interactions", "population", "free"), ...) {
    .mfrm_validate_fit(object)
    type <- match.arg(type)
    if (identical(type, "locations"))
        type <- "facets"
    if (identical(type, "persons"))
        type <- "person"
    if (identical(type, "free"))
        return(object$xi)
    if (identical(type, "facets")) {
        return(stats::setNames(object$facets$Estimate, paste0(object$facets$Facet, ":", object$facets$Level)))
    }
    if (identical(type, "steps")) {
        return(stats::setNames(object$steps$Estimate, paste0(object$steps$ScaleLevel, ":", object$steps$Step)))
    }
    if (identical(type, "slopes")) {
        if (is.null(object$slopes) || nrow(object$slopes) == 0L)
            return(numeric(0))
        return(stats::setNames(object$slopes$Estimate, paste0(object$slopes$Facet, ":", object$slopes$Level)))
    }
    if (identical(type, "person")) {
        return(stats::setNames(object$person$Estimate, paste0("Person:", object$person$Person)))
    }
    if (identical(type, "explanatory")) {
        if (is.null(object$explanatory))
            return(numeric(0))
        return(stats::setNames(object$explanatory$Estimate, paste0(object$explanatory$Facet, ":", object$explanatory$Term)))
    }
    if (identical(type, "interactions")) {
        if (is.null(object$interactions))
            return(numeric(0))
        return(stats::setNames(object$interactions$Estimate, paste0(object$interactions$Interaction, ":",
            object$interactions$LevelA, ":", object$interactions$LevelB)))
    }
    if (identical(type, "population")) {
        if (is.null(object$population$table))
            return(c(SD = object$population$sigma))
        return(stats::setNames(object$population$table$Estimate, paste0("Pop:", object$population$table$Term)))
    }
    c(stats::setNames(object$facets$Estimate, paste0(object$facets$Facet, ":", object$facets$Level)), stats::setNames(object$steps$Estimate,
        paste0(object$steps$ScaleLevel, ":", object$steps$Step)), if (!is.null(object$slopes) && nrow(object$slopes) >
        0L) {
        stats::setNames(object$slopes$Estimate, paste0("Slope:", object$slopes$Facet, ":", object$slopes$Level))
    }, if (!is.null(object$explanatory)) {
        stats::setNames(object$explanatory$Estimate, paste0("Explanatory:", object$explanatory$Facet, ":",
            object$explanatory$Term))
    }, if (!is.null(object$population$table)) {
        stats::setNames(object$population$table$Estimate, paste0("Pop:", object$population$table$Term))
    })
}


vcov.mfrm_fit <- function(object, component = c("design", "facets", "steps", "slopes", "explanatory", "population"),
    se_type = c("model", "sandwich"), facet = NULL, ...) {
    .mfrm_validate_fit(object)
    if (!is.null(facet) && missing(component))
        component <- "facets"
    component <- match.arg(component)
    se_type <- match.arg(se_type)
    use_sw <- identical(se_type, "sandwich") && identical(object$method, "MML")
    if (identical(component, "design")) {
        return(if (use_sw) object$xi_covariance_sandwich else object$xi_covariance_model)
    }
    if (identical(component, "steps")) {
        return(if (use_sw) object$step_covariance_sandwich else object$step_covariance_model)
    }
    if (identical(component, "slopes")) {
        sl_cov <- if (use_sw)
            object$slope_covariance_sandwich
        else object$slope_covariance_model
        if (is.null(sl_cov) || length(sl_cov) == 0L)
            return(matrix(0, 0L, 0L))
        if (!is.null(facet)) {
            if (!facet %in% names(sl_cov)) {
                stop("`facet` must be one of: ", paste(names(sl_cov), collapse = ", "), ".", call. = FALSE)
            }
            return(sl_cov[[facet]])
        }
        if (length(sl_cov) == 1L)
            return(sl_cov[[1L]])
        return(sl_cov)
    }
    if (identical(component, "explanatory")) {
        V <- if (use_sw)
            object$xi_covariance_sandwich
        else object$xi_covariance_model
        idx <- unlist(lapply(object$encoding$explanatory_meta, function(em) {
            f <- em$facet_idx
            object$xi_layout$facet_start[f] + seq_len(object$xi_layout$facet_len[f])
        }))
        if (length(idx) == 0L || nrow(V) < max(idx))
            return(matrix(0, 0L, 0L))
        return(V[idx, idx, drop = FALSE])
    }
    if (identical(component, "population")) {
        V <- if (use_sw)
            object$xi_covariance_sandwich
        else object$xi_covariance_model
        n_xi <- object$xi_layout$n_xi
        if (nrow(V) <= n_xi)
            return(matrix(0, 0L, 0L))
        idx <- seq.int(n_xi + 1L, nrow(V))
        return(V[idx, idx, drop = FALSE])
    }
    cov_list <- if (use_sw)
        object$covariance_sandwich
    else object$covariance_model
    if (is.null(facet)) {
        if (length(object$encoding$base_facets) == 1L)
            return(cov_list[[1L]])
        return(cov_list[object$encoding$base_facets])
    }
    if (!facet %in% names(cov_list)) {
        stop("`facet` must be one of: ", paste(names(cov_list), collapse = ", "), ".", call. = FALSE)
    }
    cov_list[[facet]]
}

# 1D C++ profile-likelihood CI helper for any free coordinate j in 1..n_xi
.mfrm_profile_xi_ci <- function(object, xi_indices, level = 0.95) {
  enc <- object$encoding
  lay <- object$xi_layout
  ll_max <- object$loglik
  crit <- 0.5 * stats::qchisq(level, df = 1L)
  target_ll <- ll_max - crit
  xi_hat <- object$xi
  V_mod <- object$xi_covariance_model

  # Build explicit affine maps for ALL facets and steps so any single xi column can be fixed cleanly
  F <- length(enc$facets)
  K <- enc$n_cat - 1L
  G <- enc$n_step_groups

  base_f_d0 <- vector("list", F)
  base_f_B <- vector("list", F)
  for (f in seq_len(F)) {
    Lf <- enc$facet_n_levels[f]
    if (enc$facet_dummy[f] == 1L) {
      base_f_d0[[f]] <- numeric(Lf)
      base_f_B[[f]] <- matrix(0.0, Lf, 0L)
    } else if (enc$has_custom_facet_map[f] == 1L) {
      base_f_d0[[f]] <- enc$custom_facet_d0[[f]]
      base_f_B[[f]] <- enc$custom_facet_B[[f]]
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
      base_f_d0[[f]] <- d0
      base_f_B[[f]] <- B
    }
  }

  if (isTRUE(enc$has_custom_step_map)) {
    base_s_d0 <- enc$custom_step_d0
    base_s_B <- enc$custom_step_B
    base_s_fix <- enc$custom_step_fixed
  } else {
    qs <- lay$step_len
    base_s_d0 <- rep(NA_real_, G * K)
    base_s_B <- matrix(0.0, G * K, qs)
    base_s_fix <- rep(1L, G * K)
    cc <- 0L
    for (g in seq_len(G)) {
      ag <- enc$step_cat_min[g]
      bg <- enc$step_cat_max[g]
      Kg <- bg - ag
      if (Kg >= 1L) base_s_d0[(g - 1L) * K + seq.int(ag + 1L, bg)] <- 0.0
      if (Kg >= 2L) {
        for (t in seq_len(Kg - 1L)) {
          rt <- (g - 1L) * K + ag + t
          rl <- (g - 1L) * K + bg
          base_s_B[rt, cc + t] <- 1.0
          base_s_B[rl, cc + t] <- -1.0
          base_s_fix[rt] <- 0L
        }
        base_s_fix[(g - 1L) * K + bg] <- 0L
        cc <- cc + (Kg - 1L)
      }
    }
  }

  eval_profile_ll <- function(j_global, val) {
    ds_p <- enc
    ds_p$has_custom_facet_map <- rep(1L, F)
    ds_p$custom_facet_d0 <- base_f_d0
    ds_p$custom_facet_B <- base_f_B
    ds_p$custom_facet_q <- vapply(base_f_B, ncol, integer(1))
    ds_p$has_custom_step_map <- TRUE
    ds_p$custom_step_d0 <- base_s_d0
    ds_p$custom_step_B <- base_s_B
    ds_p$custom_step_q <- ncol(base_s_B)
    ds_p$custom_step_fixed <- base_s_fix

    if (j_global <= lay$step_start) {
      f_hit <- which(lay$facet_start < j_global & (lay$facet_start + lay$facet_len) >= j_global)[1L]
      c_loc <- j_global - lay$facet_start[f_hit]
      Bf <- ds_p$custom_facet_B[[f_hit]]
      ds_p$custom_facet_d0[[f_hit]] <- ds_p$custom_facet_d0[[f_hit]] + Bf[, c_loc] * val
      ds_p$custom_facet_B[[f_hit]] <- Bf[, -c_loc, drop = FALSE]
      ds_p$custom_facet_q[f_hit] <- ncol(ds_p$custom_facet_B[[f_hit]])
    } else if (j_global <= (lay$slope_start %||% (lay$step_start + lay$step_len))) {
      c_loc <- j_global - lay$step_start
      Bs <- ds_p$custom_step_B
      act <- is.finite(ds_p$custom_step_d0)
      ds_p$custom_step_d0[act] <- ds_p$custom_step_d0[act] + Bs[act, c_loc] * val
      ds_p$custom_step_B <- Bs[, -c_loc, drop = FALSE]
      ds_p$custom_step_q <- ncol(ds_p$custom_step_B)
    } else {
      sl_lens <- lay$slope_len
      sl_starts <- lay$slope_start + c(0L, cumsum(sl_lens)[-length(sl_lens)])
      s_hit <- which(sl_starts < j_global & (sl_starts + sl_lens) >= j_global)[1L]
      c_loc <- j_global - sl_starts[s_hit]
      Bsl <- ds_p$slope_B[[s_hit]]
      ds_p$slope_d0[[s_hit]] <- ds_p$slope_d0[[s_hit]] + Bsl[, c_loc] * val
      ds_p$slope_B[[s_hit]] <- Bsl[, -c_loc, drop = FALSE]
      ds_p$slope_q[s_hit] <- ncol(ds_p$slope_B[[s_hit]])
    }

    ctrl_p <- object$control
    ctrl_p$compute_diagnostics <- FALSE
    ctrl_p$compute_wle <- FALSE
    ctrl_p$bias_facet_a <- -1L
    ctrl_p$bias_facet_b <- -1L
    ctrl_p$warm_init <- as.double(xi_hat[-j_global])
    rw <- mfrm_cpp_fast_fit_spec(ds_p, ctrl_p)
    as.numeric(rw$loglik)
  }

  out_ci <- matrix(NA_real_, nrow = length(xi_indices), ncol = 2L)
  z_mult <- stats::qnorm(0.5 + 0.5 * level)
  for (ii in seq_along(xi_indices)) {
    j <- xi_indices[ii]
    est <- unname(xi_hat[j])
    se0 <- if (nrow(V_mod) >= j && is.finite(V_mod[j, j]) && V_mod[j, j] > 0) sqrt(V_mod[j, j]) else 0.25
    for (dir in c(-1.0, 1.0)) {
      step <- dir * max(1e-3, z_mult * se0)
      x_lo <- est
      x_hi <- est + step
      f_hi <- eval_profile_ll(j, x_hi) - target_ll
      for (k_exp in seq_len(8L)) {
        if (!is.finite(f_hi) || f_hi <= 0.0) break
        x_lo <- x_hi
        step <- step * 1.8
        x_hi <- est + step
        f_hi <- eval_profile_ll(j, x_hi) - target_ll
      }
      if (is.finite(f_hi) && f_hi <= 0.0) {
        iv <- sort(c(x_lo, x_hi))
        rt <- tryCatch(
          stats::uniroot(function(v) eval_profile_ll(j, v) - target_ll, interval = iv, tol = 1e-4)$root,
          error = function(e) NA_real_
        )
        out_ci[ii, if (dir < 0) 1L else 2L] <- rt
      }
    }
  }
  out_ci
}


confint.mfrm_fit <- function(object, parm = NULL, level = 0.95, component = c("facets", "steps", "slopes",
    "person", "explanatory", "interactions", "population", "free"), se_type = c("model", "real", "sandwich",
    "profile"), ...) {
    .mfrm_validate_fit(object)
    readiness <- .mfrm_readiness(object)
    if (!readiness$ready)
        stop("Calibration is not ready for intervals: ", paste(readiness$reason_codes, collapse = ", "), ".",
            call. = FALSE)
    if (!is.null(parm) && is.character(parm) && length(parm) == 1L && parm %in% c("facets", "locations", "steps",
        "slopes", "person", "persons", "explanatory", "interactions", "population", "free")) {
        component <- switch(parm, locations = "facets", persons = "person", parm)
        parm <- NULL
    }
    component <- match.arg(component)
    se_type <- match.arg(se_type)
    z <- stats::qnorm(0.5 + 0.5 * level)
    pct_lbls <- paste0(format(100 * c((1 - level)/2, (1 + level)/2), trim = TRUE), "%")
    if (identical(se_type, "profile")) {
        if (identical(component, "slopes")) {
            sl_start <- object$xi_layout$slope_start %||% object$xi_layout$n_xi
            n_xi <- object$xi_layout$n_xi
            if (sl_start >= n_xi)
                return(matrix(NA_real_, 0L, 2L, dimnames = list(character(0), pct_lbls)))
            idx <- seq.int(sl_start + 1L, n_xi)
            nms <- sub("^LogSlope:", "", names(object$xi)[idx])
            if (!is.null(parm)) {
                sel <- if (is.character(parm))
                  match(parm, nms)
                else as.integer(parm)
                sel <- sel[!is.na(sel) & sel >= 1L & sel <= length(idx)]
                idx <- idx[sel]
                nms <- nms[sel]
            }
            ci_mat <- exp(.mfrm_profile_xi_ci(object, idx, level = level))
            dimnames(ci_mat) <- list(nms, pct_lbls)
            return(ci_mat)
        }
        n_xi <- object$xi_layout$n_xi
        idx <- seq_len(n_xi)
        nms <- names(object$xi)[idx]
        if (!is.null(parm)) {
            idx <- if (is.character(parm))
                match(parm, nms)
            else as.integer(parm)
            idx <- idx[!is.na(idx) & idx >= 1L & idx <= n_xi]
        }
        ci_mat <- .mfrm_profile_xi_ci(object, idx, level = level)
        dimnames(ci_mat) <- list(nms[idx], pct_lbls)
        return(ci_mat)
    }
    pick_se <- function(mod_se, real_se = NULL, sw_se = NULL) {
        if (identical(se_type, "real") && !is.null(real_se))
            return(real_se)
        if (identical(se_type, "sandwich") && !is.null(sw_se) && identical(object$method, "MML"))
            return(sw_se)
        mod_se
    }
    if (identical(component, "facets")) {
        tb <- object$facets
        est <- tb$Estimate
        se <- pick_se(tb$ModelSE, tb$RealSE, tb$SandwichSE)
        nms <- paste0(tb$Facet, ":", tb$Level)
    }
    else if (identical(component, "steps")) {
        tb <- object$steps
        est <- tb$Estimate
        se <- pick_se(tb$ModelSE, NULL, tb$SandwichSE)
        nms <- paste0(tb$ScaleLevel, ":", tb$Step)
    }
    else if (identical(component, "slopes")) {
        tb <- object$slopes
        if (is.null(tb) || nrow(tb) == 0L)
            return(matrix(NA_real_, 0L, 2L, dimnames = list(character(0), pct_lbls)))
        est <- log(tb$Estimate)
        se <- pick_se(tb$ModelSE, NULL, tb$SandwichSE)/tb$Estimate
        nms <- paste0(tb$Facet, ":", tb$Level)
    }
    else if (identical(component, "person")) {
        tb <- object$person
        est <- tb$Estimate
        se <- pick_se(tb$ModelSE, tb$RealSE, NULL)
        nms <- paste0("Person:", tb$Person)
    }
    else if (identical(component, "explanatory")) {
        tb <- object$explanatory
        if (is.null(tb))
            return(matrix(NA_real_, 0L, 2L, dimnames = list(character(0), pct_lbls)))
        est <- tb$Estimate
        se <- pick_se(tb$ModelSE, NULL, tb$SandwichSE)
        nms <- paste0(tb$Facet, ":", tb$Term)
    }
    else if (identical(component, "interactions")) {
        tb <- object$interactions
        if (is.null(tb))
            return(matrix(NA_real_, 0L, 2L, dimnames = list(character(0), pct_lbls)))
        est <- tb$Estimate
        se <- pick_se(tb$ModelSE, NULL, tb$SandwichSE)
        nms <- paste0(tb$Interaction, ":", tb$LevelA, ":", tb$LevelB)
    }
    else if (identical(component, "population")) {
        tb <- object$population$table
        if (is.null(tb))
            return(matrix(NA_real_, 0L, 2L, dimnames = list(character(0), pct_lbls)))
        est <- tb$Estimate
        se <- pick_se(tb$ModelSE, NULL, tb$SandwichSE)
        nms <- paste0("Pop:", tb$Term)
    }
    else {
        est <- unname(object$xi)
        V <- if (identical(se_type, "sandwich") && identical(object$method, "MML"))
            object$xi_covariance_sandwich
        else object$xi_covariance_model
        se <- if (nrow(V) == length(est))
            sqrt(pmax(0, diag(V)))
        else rep(NA_real_, length(est))
        nms <- names(object$xi)
    }
    ok <- is.finite(est) & is.finite(se) & (se > 0)
    lo <- ifelse(ok, est - z * se, NA_real_)
    hi <- ifelse(ok, est + z * se, NA_real_)
    out <- cbind(lo, hi)
    if (identical(component, "slopes"))
        out <- exp(out)
    dimnames(out) <- list(nms, pct_lbls)
    if (!is.null(parm)) {
        sel <- if (is.character(parm))
            match(parm, nms)
        else as.integer(parm)
        if (anyNA(sel))
            stop("Unknown parameter(s) in `parm`.", call. = FALSE)
        out <- out[sel, , drop = FALSE]
    }
    out
}


as.data.frame.mfrm_fit <- function(x, row.names = NULL, optional = FALSE, component = c("facets", "person",
    "steps", "slopes"), ...) {
    .mfrm_validate_fit(x)
    component <- match.arg(component)
    x[[component]]
}

.mfrm_model_label <- function(model, slope_facets = character()) {
  if (!identical(model, "GPCM")) return(model)
  label <- "Generalized Many-Facet Rasch Model (GMFRM)"
  if (!length(slope_facets)) return(label)
  paste0(label, ", discrimination by ", paste(slope_facets, collapse = ", "))
}


print.mfrm_fit <- function(x, ...) {
    .mfrm_validate_fit(x)
    cat(sprintf("Many-Facet Measurement Fit [%s / %s | Score: %s]\n", x$method, .mfrm_model_label(x$model,
        unique(x$slopes$Facet)), x$score_method))
    cat(sprintf("  Ratings: %d | Persons: %d | Free par: %d | Connected: %s (%d subset%s)\n", x$n_obs, x$n_person,
        x$n_free_par, if (x$connectivity$connected)
            "YES"
        else "NO", x$connectivity$n_subsets, if (x$connectivity$n_subsets == 1L)
            ""
        else "s"))
    cat(sprintf("  Converged: %s (%d iter, max|grad| = %.2e, %.4f s) | LogLik = %.4f | AIC = %.2f\n", if (x$converged)
        "TRUE"
    else "FALSE", x$iterations, x$max_abs_grad, x$elapsed_sec, x$loglik, x$aic))
    readiness <- .mfrm_readiness(x)
    cat("  Inference:", readiness$status, if (length(readiness$reason_codes))
        paste0(" (", paste(readiness$reason_codes, collapse = ", "), ")")
    else "", "\n")
    if (identical(x$method, "MML")) {
        cat(sprintf("  Population SD = %.4f | EAP Reliability = %.4f\n", x$population$sigma, x$population$eap_reliability))
    }
    if (!is.null(x$slopes) && nrow(x$slopes) > 0L) {
        cat(sprintf("  GMFRM discrimination (%s): %d element(s), range [%.3f, %.3f]\n", paste(unique(x$slopes$Facet),
            collapse = ", "), nrow(x$slopes), min(x$slopes$Estimate, na.rm = TRUE), max(x$slopes$Estimate,
            na.rm = TRUE)))
    }
    if (!is.null(x$explanatory)) {
        cat("  Explanatory Covariate Terms:", nrow(x$explanatory), "\n")
    }
    invisible(x)
}


summary.mfrm_fit <- function(object, ...) {
    .mfrm_validate_fit(object)
    readiness <- .mfrm_readiness(object)
    overview <- object$summary
    overview$InferenceReady <- readiness$ready
    overview$InferenceStatus <- readiness$status
    list(overview = overview, model_label = .mfrm_model_label(object$model, unique(object$slopes$Facet)),
        readiness = readiness, connectivity = object$connectivity, population = object$population, facets = object$facets,
        steps = object$steps, slopes = object$slopes, explanatory = object$explanatory, interactions = object$interactions,
        person_summary = data.frame(MeanAbility = mean(object$person$Estimate[is.finite(object$person$Estimate)]),
            SDAbility = stats::sd(object$person$Estimate[is.finite(object$person$Estimate)]), MeanModelSE = mean(object$person$ModelSE,
                na.rm = TRUE), MeanRealSE = mean(object$person$RealSE, na.rm = TRUE), MeanInfit = mean(object$person$InfitMnSq,
                na.rm = TRUE), MeanOutfit = mean(object$person$OutfitMnSq, na.rm = TRUE), ExtremeMin = sum(object$person$ExtremeFlag ==
                -1L), ExtremeMax = sum(object$person$ExtremeFlag == 1L)))
}
