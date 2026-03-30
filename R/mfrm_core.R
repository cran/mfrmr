# Core MFRM computation functions extracted from app.R (non-UI section)
# Source: ../../app.R
# NOTE: This file intentionally keeps function names aligned with the app implementation.
# ---- math helpers ----

# Numerically stable log-sum-exp: log(sum(exp(x))).
# Subtracts max(x) before exponentiation to avoid overflow/underflow.
#   logsumexp(x) = max(x) + log(sum(exp(x - max(x))))
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Weighted mean that silently drops non-finite values and zero weights.
weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

row_max_fast <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(numeric(0))
  mat[cbind(seq_len(nrow(mat)), max.col(mat, ties.method = "first"))]
}

get_weights <- function(df) {
  if (!is.null(df) && "Weight" %in% names(df)) {
    w <- suppressWarnings(as.numeric(df$Weight))
    w <- ifelse(is.finite(w) & w > 0, w, 0)
    return(w)
  }
  rep(1, nrow(df))
}

# Gauss-Hermite nodes/weights for standard normal integration
# Based on Golub-Welsch for Hermite polynomials
# Returns nodes and weights for phi(x) dx
#   integral f(x) phi(x) dx \approx sum w_i f(x_i)
gauss_hermite_normal <- function(n) {
  if (n < 1) stop("Gauss-Hermite quadrature requires n >= 1 quadrature points. ",
                   "Check the 'quad_points' argument.", call. = FALSE)
  if (n == 1) {
    return(list(nodes = 0, weights = 1))
  }
  i <- seq_len(n - 1)
  a <- rep(0, n)
  b <- sqrt(i / 2)
  jmat <- matrix(0, nrow = n, ncol = n)
  diag(jmat) <- a
  jmat[cbind(seq_len(n - 1), seq_len(n - 1) + 1)] <- b
  jmat[cbind(seq_len(n - 1) + 1, seq_len(n - 1))] <- b
  eig <- eigen(jmat, symmetric = TRUE)
  nodes <- eig$values
  weights <- sqrt(pi) * (eig$vectors[1, ]^2)
  # Convert from exp(-x^2) to standard normal
  nodes <- sqrt(2) * nodes
  weights <- weights / sqrt(pi)
  list(nodes = nodes, weights = weights)
}

center_sum_zero <- function(x) {
  if (length(x) == 0) return(x)
  x - mean(x)
}

expand_facet <- function(free, n_levels) {
  if (n_levels <= 1) return(rep(0, n_levels))
  c(free, -sum(free))
}

build_facet_constraint <- function(levels,
                                   anchors = NULL,
                                   groups = NULL,
                                   group_values = NULL,
                                   centered = TRUE) {
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

  spec <- list(
    levels = lvl,
    anchors = anchors_vec,
    groups = groups_vec,
    group_values = group_values_map,
    centered = centered
  )
  spec$n_params <- count_facet_params(spec)
  spec
}

count_facet_params <- function(spec) {
  anchors <- spec$anchors
  groups <- spec$groups
  free_idx <- which(is.na(anchors))
  if (length(free_idx) == 0) return(0)

  n_params <- 0
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(anchors[group_levels])]
      k <- length(free_in_group)
      if (k > 1) n_params <- n_params + (k - 1)
    }
  }

  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m == 0) return(n_params)
  if (spec$centered) {
    n_params <- n_params + max(m - 1, 0)
  } else {
    n_params <- n_params + m
  }
  n_params
}

expand_facet_with_constraints <- function(free, spec) {
  out <- spec$anchors
  groups <- spec$groups
  group_values <- spec$group_values
  centered <- isTRUE(spec$centered)
  free_idx <- which(is.na(out))
  if (length(free_idx) == 0) return(out)

  used <- 0
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(out[group_levels])]
      if (length(free_in_group) == 0) next
      group_value <- if (gid %in% names(group_values)) group_values[[gid]] else 0
      anchor_sum <- sum(out[group_levels], na.rm = TRUE)
      target_sum <- group_value * length(group_levels)
      k <- length(free_in_group)
      if (k == 1) {
        out[free_in_group] <- target_sum - anchor_sum
      } else {
        seg <- free[(used + 1):(used + k - 1)]
        used <- used + (k - 1)
        last_val <- target_sum - anchor_sum - sum(seg)
        out[free_in_group] <- c(seg, last_val)
      }
    }
  }

  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m == 0) return(out)
  if (centered) {
    if (m == 1) {
      out[ungrouped_idx] <- 0
    } else {
      seg <- free[(used + 1):(used + m - 1)]
      used <- used + (m - 1)
      out[ungrouped_idx] <- c(seg, -sum(seg))
    }
  } else {
    seg <- free[(used + 1):(used + m)]
    used <- used + m
    out[ungrouped_idx] <- seg
  }
  out
}

# ---- Constraint gradient projection ----
# Reverse of expand_facet_with_constraints: projects gradient from expanded
# (full-length) parameter space back to the free (optimizable) parameter space.
# This computes J^T * grad_expanded, where J is the Jacobian of
# expand_facet_with_constraints w.r.t. the free parameters.
constraint_grad_project <- function(grad_expanded, spec) {
  anchors <- spec$anchors
  groups <- spec$groups
  centered <- isTRUE(spec$centered)
  free_idx <- which(is.na(anchors))
  if (length(free_idx) == 0 || spec$n_params == 0) return(numeric(0))

  grad_free <- numeric(spec$n_params)
  pos <- 1L

  # Grouped levels: last free in group = target - anchor_sum - sum(others)
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(anchors[group_levels])]
      k <- length(free_in_group)
      if (k <= 1) next
      last_g <- grad_expanded[free_in_group[k]]
      for (j in seq_len(k - 1)) {
        grad_free[pos] <- grad_expanded[free_in_group[j]] - last_g
        pos <- pos + 1L
      }
    }
  }

  # Ungrouped levels
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
      # m == 1: constrained to 0, no free params
    } else {
      for (j in seq_len(m)) {
        grad_free[pos] <- grad_expanded[ungrouped_idx[j]]
        pos <- pos + 1L
      }
    }
  }

  grad_free
}

build_param_sizes <- function(config) {
  n_steps <- max(config$n_cat - 1, 0)
  sizes <- list(
    theta = if (config$method == "JMLE") config$theta_spec$n_params else 0
  )
  for (facet in config$facet_names) {
    sizes[[facet]] <- config$facet_specs[[facet]]$n_params
  }
  if (config$model == "RSM") {
    sizes$steps <- n_steps
  } else {
    if (is.null(config$step_facet) || !config$step_facet %in% config$facet_names) {
      stop("PCM model requires 'step_facet' to name one of the facet columns: ",
           paste(config$facet_names, collapse = ", "), ". ",
           "Supply step_facet = '<name>'.", call. = FALSE)
    }
    sizes$steps <- length(config$facet_levels[[config$step_facet]]) * n_steps
  }
  sizes
}

split_params <- function(par, sizes) {
  out <- list()
  idx <- 1
  for (nm in names(sizes)) {
    k <- sizes[[nm]]
    if (k == 0) {
      out[[nm]] <- numeric(0)
    } else {
      out[[nm]] <- par[idx:(idx + k - 1)]
      idx <- idx + k
    }
  }
  out
}

expand_params <- function(par, sizes, config) {
  parts <- split_params(par, sizes)
  theta <- if (config$method == "JMLE") {
    expand_facet_with_constraints(parts$theta, config$theta_spec)
  } else {
    numeric(0)
  }

  facets <- lapply(config$facet_names, function(facet) {
    expand_facet_with_constraints(parts[[facet]], config$facet_specs[[facet]])
  })
  names(facets) <- config$facet_names

  if (config$model == "RSM") {
    steps <- center_sum_zero(parts$steps)
    steps_mat <- NULL
  } else {
    n_levels <- length(config$facet_levels[[config$step_facet]])
    if (n_levels == 0 || config$n_cat <= 1) {
      steps_mat <- matrix(0, nrow = n_levels, ncol = max(config$n_cat - 1, 0))
    } else {
      steps_mat <- matrix(parts$steps, nrow = n_levels, byrow = TRUE)
      steps_mat <- t(apply(steps_mat, 1, center_sum_zero))
    }
    steps <- NULL
  }

  list(theta = theta, facets = facets, steps = steps, steps_mat = steps_mat)
}

# ---- data preparation ----
prepare_mfrm_data <- function(data, person_col, facet_cols, score_col,
                              rating_min = NULL, rating_max = NULL,
                              weight_col = NULL, keep_original = FALSE) {
  required <- c(person_col, facet_cols, score_col)
  if (!is.null(weight_col)) {
    required <- c(required, weight_col)
  }
  if (length(unique(required)) != length(required)) {
    dup_names <- required[duplicated(required)]
    stop("The 'person', 'score', and 'facets' arguments must name distinct columns, ",
         "but duplicates were found: ", paste(dup_names, collapse = ", "), ". ",
         "Remove or rename the duplicated references.", call. = FALSE)
  }
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "), ". ",
         "Available columns: ", paste(names(data), collapse = ", "), ". ",
         "Check spelling of person/facets/score arguments.", call. = FALSE)
  }
  if (any(duplicated(names(data)))) {
    dupes <- unique(names(data)[duplicated(names(data))])
    if (any(required %in% dupes)) {
      stop("Selected columns include duplicate names in the data: ",
           paste(intersect(required, dupes), collapse = ", "), ". ",
           "Rename columns so each name is unique.", call. = FALSE)
    }
  }
  if (length(facet_cols) == 0) {
    stop("No facet columns were specified. ",
         "Supply at least one column name via 'facets' ",
         "(e.g., facets = c('Rater', 'Task')).", call. = FALSE)
  }

  cols <- c(person_col, facet_cols, score_col)
  if (!is.null(weight_col)) {
    cols <- c(cols, weight_col)
  }
  df <- data |>
    select(all_of(cols)) |>
    rename(
      Person = all_of(person_col),
      Score = all_of(score_col)
    )
  if (!is.null(weight_col)) {
    df <- df |> rename(Weight = all_of(weight_col))
  }

  raw_score <- as.character(df$Score)
  raw_weight <- if ("Weight" %in% names(df)) as.character(df$Weight) else NULL

  score_num <- suppressWarnings(as.numeric(raw_score))
  bad_score <- is.na(score_num) & !is.na(raw_score) & nzchar(trimws(raw_score))
  if (any(bad_score)) {
    warning(
      "`Score` contained ", sum(bad_score), " non-numeric value(s); affected row(s) will be removed before estimation.",
      call. = FALSE
    )
  }
  df <- df |>
    mutate(
      Person = as.character(Person),
      across(all_of(facet_cols), ~ as.character(.x)),
      Score = score_num
    )
  if (!"Weight" %in% names(df)) {
    df <- df |> mutate(Weight = 1)
  } else {
    weight_num <- suppressWarnings(as.numeric(raw_weight))
    bad_weight <- is.na(weight_num) & !is.na(raw_weight) & nzchar(trimws(raw_weight))
    if (any(bad_weight)) {
      warning(
        "`Weight` contained ", sum(bad_weight), " non-numeric value(s); affected row(s) will be removed before estimation.",
        call. = FALSE
      )
    }
    df <- df |> mutate(Weight = weight_num)
  }

  df <- df |>
    tidyr::drop_na() |>
    filter(Weight > 0)

  if (nrow(df) == 0) {
    stop("No valid observations remain after removing missing values and ",
         "zero-weight rows. Check that person, facet, score, and weight columns ",
         "contain valid (non-NA, non-empty) data.", call. = FALSE)
  }

  df <- df |>
    mutate(Score = as.integer(Score))

  observed_score_values <- sort(unique(df$Score))

  if (length(unique(df$Score)) < 2) {
    stop("Only one score category found in the data (Score = ",
         unique(df$Score), "). ",
         "MFRM requires at least two distinct response categories.", call. = FALSE)
  }

  if (is.null(rating_min)) rating_min <- min(df$Score, na.rm = TRUE)
  if (is.null(rating_max)) rating_max <- max(df$Score, na.rm = TRUE)

  if (!isTRUE(keep_original)) {
    score_vals <- sort(unique(df$Score))
    expected_vals <- seq(rating_min, rating_max)
    if (!identical(score_vals, expected_vals)) {
      df <- df |>
        mutate(Score = match(Score, score_vals) + rating_min - 1)
      rating_max <- rating_min + length(score_vals) - 1
    }
  }

  if (isTRUE(keep_original)) {
    score_map <- tibble(
      OriginalScore = seq(rating_min, rating_max),
      InternalScore = seq(rating_min, rating_max)
    )
  } else {
    score_map <- tibble(
      OriginalScore = observed_score_values,
      InternalScore = seq(rating_min, rating_max)
    )
  }

  df <- df |>
    mutate(score_k = Score - rating_min)

  df <- df |>
    mutate(
      Person = factor(Person),
      across(all_of(facet_cols), ~ factor(.x))
    )

  facet_names <- facet_cols
  facet_levels <- lapply(facet_names, function(f) levels(df[[f]]))
  names(facet_levels) <- facet_names

  list(
    data = df,
    n_obs = nrow(df),
    weighted_n = sum(df$Weight, na.rm = TRUE),
    n_person = length(levels(df$Person)),
    rating_min = rating_min,
    rating_max = rating_max,
    score_map = score_map,
    facet_names = facet_names,
    levels = c(list(Person = levels(df$Person)), facet_levels),
    weight_col = if (!is.null(weight_col)) weight_col else NULL,
    keep_original = isTRUE(keep_original),
    source_columns = list(
      person = person_col,
      facets = facet_cols,
      score = score_col,
      weight = if (!is.null(weight_col)) weight_col else NULL
    )
  )
}

build_indices <- function(prep, step_facet = NULL) {
  df <- prep$data
  facets_idx <- lapply(prep$facet_names, function(f) as.integer(df[[f]]))
  names(facets_idx) <- prep$facet_names
  step_idx <- if (!is.null(step_facet)) {
    as.integer(df[[step_facet]])
  } else {
    NULL
  }
  # Pre-split observation indices by criterion for PCM (avoids repeated which())
  criterion_splits <- if (!is.null(step_idx)) {
    split(seq_len(nrow(df)), step_idx)
  } else {
    NULL
  }
  list(
    person = as.integer(df$Person),
    facets = facets_idx,
    step_idx = step_idx,
    criterion_splits = criterion_splits,
    score_k = as.integer(df$score_k),
    weight = suppressWarnings(as.numeric(df$Weight))
  )
}

sample_mfrm_data <- function(seed = 20240131) {
  with_preserved_rng_seed(seed, {
    persons <- paste0("P", sprintf("%02d", 1:36))
    raters <- paste0("R", 1:3)
    tasks <- paste0("T", 1:4)
    criteria <- paste0("C", 1:3)
    df <- expand_grid(
      Person = persons,
      Rater = raters,
      Task = tasks,
      Criterion = criteria
    )
    ability <- rnorm(length(persons), 0, 1)
    rater_eff <- c(-0.4, 0, 0.4)
    task_eff <- seq(-0.5, 0.5, length.out = length(tasks))
    crit_eff <- c(-0.3, 0, 0.3)
    eta <- ability[match(df$Person, persons)] -
      rater_eff[match(df$Rater, raters)] -
      task_eff[match(df$Task, tasks)] -
      crit_eff[match(df$Criterion, criteria)]
    raw <- eta + rnorm(nrow(df), 0, 0.6)
    score <- as.integer(cut(
      raw,
      breaks = c(-Inf, -1.0, -0.3, 0.3, 1.0, Inf),
      labels = 1:5
    ))
    df$Score <- score
    df
  })
}

with_preserved_rng_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed[1]))
  force(expr)
}

format_tab_template <- function(df) {
  char_df <- df |> mutate(across(everything(), ~ replace_na(as.character(.x), "")))
  widths <- vapply(seq_along(char_df), function(i) {
    max(nchar(c(names(char_df)[i], char_df[[i]])), na.rm = TRUE)
  }, integer(1))
  format_row <- function(row_vec) {
    padded <- mapply(function(value, width) {
      value <- ifelse(is.na(value), "", value)
      stringr::str_pad(value, width = width, side = "right")
    }, row_vec, widths, SIMPLIFY = TRUE)
    paste(padded, collapse = "\t")
  }
  header <- format_row(names(char_df))
  rows <- apply(char_df, 1, format_row)
  paste(c(header, rows), collapse = "\n")
}

template_tab_source_demo <- sample_mfrm_data(seed = 20240131) |>
  slice_head(n = 24)
template_tab_source_toy <- sample_mfrm_data(seed = 20240131) |>
  slice_head(n = 8)
template_tab_text <- format_tab_template(template_tab_source_demo)
template_tab_text_toy <- format_tab_template(template_tab_source_toy)
template_header_text <- format_tab_template(template_tab_source_demo[0, ])
download_sample_data <- sample_mfrm_data(seed = 20240131)

guess_col <- function(cols, patterns, fallback = 1) {
  if (length(cols) == 0) return(character(0))
  hit <- which(stringr::str_detect(tolower(cols), paste(patterns, collapse = "|")))
  if (length(hit) > 0) return(cols[hit[1]])
  cols[min(fallback, length(cols))]
}

truncate_label <- function(x, width = 28) {
  stringr::str_trunc(as.character(x), width = width)
}

facet_report_id <- function(facet) {
  paste0("facet_report_", stringr::str_replace_all(as.character(facet), "[^A-Za-z0-9]", "_"))
}

# ---- likelihoods ----
# RSM log-likelihood: sum_i w_i log P(X_i = k_i | eta_i).
# Under the Rating Scale Model (Andrich, 1978):
#   P(X = k | eta) = exp(k*eta - tau_k) / sum_j exp(j*eta - tau_j)
# where tau_k = cumulative step parameters and eta = theta - sum(facets).
# Computation uses log-domain subtraction with logsumexp for stability.
loglik_rsm <- function(eta, score_k, step_cum, weight = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- length(step_cum)
  eta_mat <- outer(eta, 0:(k_cat - 1))
  log_num <- eta_mat - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  row_max <- row_max_fast(log_num)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  log_num_obs <- log_num[cbind(seq_len(n), score_k + 1)]
  diff <- log_num_obs - log_denom
  if (is.null(weight)) {
    sum(diff)
  } else {
    sum(diff * weight)
  }
}

# PCM log-likelihood: same structure as RSM but with criterion-specific steps.
# Under the Partial Credit Model (Masters, 1982):
#   P(X = k | eta, criterion c) = exp(k*eta - tau_{c,k}) / sum_j exp(j*eta - tau_{c,j})
# step_cum_mat has one row per criterion level, columns = cumulative thresholds.
loglik_pcm <- function(eta, score_k, step_cum_mat, criterion_idx, weight = NULL,
                       criterion_splits = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- ncol(step_cum_mat)
  total <- 0
  splits <- criterion_splits %||% split(seq_len(n), criterion_idx)
  for (ci in seq_along(splits)) {
    rows <- splits[[ci]]
    if (length(rows) == 0) next
    c_idx <- as.integer(names(splits)[ci])
    eta_c <- eta[rows]
    step_cum <- step_cum_mat[c_idx, ]
    nr <- length(rows)
    eta_mat <- outer(eta_c, 0:(k_cat - 1))
    log_num <- eta_mat - matrix(step_cum, nrow = nr, ncol = k_cat, byrow = TRUE)
    row_max <- log_num[cbind(seq_len(nr), max.col(log_num))]
    log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
    log_num_obs <- log_num[cbind(seq_len(nr), score_k[rows] + 1)]
    diff <- log_num_obs - log_denom
    if (is.null(weight)) {
      total <- total + sum(diff)
    } else {
      total <- total + sum(diff * weight[rows])
    }
  }
  total
}

# Category response probabilities under RSM.
# Returns an n x K matrix where K = number of categories.
# Each row sums to 1; probabilities are computed in log-domain for stability.
category_prob_rsm <- function(eta, step_cum) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = length(step_cum)))
  k_cat <- length(step_cum)
  eta_mat <- outer(eta, 0:(k_cat - 1))
  log_num <- eta_mat - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  row_max <- row_max_fast(log_num)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  exp(log_num - matrix(log_denom, nrow = n, ncol = k_cat))
}

# Category response probabilities under PCM (criterion-specific thresholds).
# Returns an n x K matrix; each row sums to 1.
category_prob_pcm <- function(eta, step_cum_mat, criterion_idx,
                              criterion_splits = NULL) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = ncol(step_cum_mat)))
  k_cat <- ncol(step_cum_mat)
  probs <- matrix(0, nrow = n, ncol = k_cat)
  splits <- criterion_splits %||% split(seq_len(n), criterion_idx)
  for (ci in seq_along(splits)) {
    rows <- splits[[ci]]
    if (length(rows) == 0) next
    c_idx <- as.integer(names(splits)[ci])
    step_cum <- step_cum_mat[c_idx, ]
    eta_c <- eta[rows]
    eta_mat <- outer(eta_c, 0:(k_cat - 1))
    log_num <- eta_mat - matrix(step_cum, nrow = length(rows), ncol = k_cat, byrow = TRUE)
    row_max <- row_max_fast(log_num)
    log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
    probs[rows, ] <- exp(log_num - matrix(log_denom, nrow = length(rows), ncol = k_cat))
  }
  probs
}

# Compute P(X >= s) matrix for s = 1,...,K-1 from category probabilities.
# Input: probs (n x K matrix of category probabilities, columns for k=0,...,K-1)
# Output: P_geq (n x (K-1) matrix), P_geq[i,s] = P(X_i >= s)
compute_P_geq <- function(probs) {
  k_cat <- ncol(probs)
  n_steps <- k_cat - 1
  if (n_steps == 0) return(matrix(0, nrow(probs), 0))
  n <- nrow(probs)
  P_geq <- matrix(0, n, n_steps)
  P_geq[, n_steps] <- probs[, k_cat]
  if (n_steps >= 2) {
    for (s in (n_steps - 1):1) {
      P_geq[, s] <- P_geq[, s + 1] + probs[, s + 1]
    }
  }
  P_geq
}

# Convert mean-square fit statistic to a standardized z-score (ZSTD).
# Default uses the Wilson-Hilferty (1931) cube-root approximation:
#   ZSTD = (MnSq^(1/3) - (1 - 2/(9*df))) / sqrt(2/(9*df))
# When whexact = TRUE, uses the simpler linear approximation:
#   ZSTD = (MnSq - 1) * sqrt(df / 2)
# Values near 0 indicate expected fit; |ZSTD| > 2 flags potential misfit.
zstd_from_mnsq <- function(mnsq, df, whexact = FALSE) {
  mnsq <- as.numeric(mnsq)
  df <- as.numeric(df)

  if (length(df) == 1L && length(mnsq) > 1L) {
    df <- rep(df, length(mnsq))
  } else if (length(mnsq) == 1L && length(df) > 1L) {
    mnsq <- rep(mnsq, length(df))
  }

  n <- min(length(mnsq), length(df))
  if (n == 0L) return(numeric(0))

  out <- rep(NA_real_, n)
  m <- mnsq[seq_len(n)]
  d <- df[seq_len(n)]
  ok <- is.finite(m) & is.finite(d) & (d > 0)

  if (isTRUE(whexact)) {
    out[ok] <- (m[ok] - 1) * sqrt(d[ok] / 2)
  } else {
    out[ok] <- (m[ok]^(1 / 3) - (1 - 2 / (9 * d[ok]))) / sqrt(2 / (9 * d[ok]))
  }
  out
}

compute_base_eta <- function(idx, params, config) {
  eta <- rep(0, length(idx$score_k))
  facet_signs <- config$facet_signs
  for (facet in config$facet_names) {
    sign <- if (!is.null(facet_signs) && !is.null(facet_signs[[facet]])) {
      facet_signs[[facet]]
    } else {
      -1
    }
    eta <- eta + sign * params$facets[[facet]][idx$facets[[facet]]]
  }
  eta
}

compute_eta <- function(idx, params, config, theta_override = NULL) {
  theta <- if (!is.null(theta_override)) theta_override else params$theta
  eta <- if (length(theta) == 0) {
    rep(0, length(idx$score_k))
  } else {
    theta[idx$person]
  }
  eta + compute_base_eta(idx, params, config)
}

# Joint Maximum Likelihood Estimation (JMLE) negative log-likelihood.
# Estimates person abilities and facet parameters simultaneously.
# Returns -LL for minimization by optim().
mfrm_loglik_jmle <- function(par, idx, config, sizes) {
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    ll <- loglik_rsm(eta, idx$score_k, step_cum, weight = idx$weight)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    ll <- loglik_pcm(eta, idx$score_k, step_cum_mat, idx$step_idx, weight = idx$weight,
                     criterion_splits = idx$criterion_splits)
  }
  -ll
}

# Cached version of JMLE log-likelihood (uses pre-computed params/eta/step_cum)
mfrm_loglik_jmle_cached <- function(cache, idx, config) {
  eta <- cache$eta()
  step_cum <- cache$step_cum()
  if (config$model == "RSM") {
    ll <- loglik_rsm(eta, idx$score_k, step_cum, weight = idx$weight)
  } else {
    ll <- loglik_pcm(eta, idx$score_k, step_cum, idx$step_idx, weight = idx$weight,
                     criterion_splits = idx$criterion_splits)
  }
  -ll
}

# Cached version of MML log-likelihood (uses pre-computed params/base_eta/step_cum)
mfrm_loglik_mml_cached <- function(cache, idx, config, quad) {
  params <- cache$params()
  base_eta <- cache$base_eta()
  step_cum <- cache$step_cum()
  n <- length(idx$score_k)
  if (n == 0) return(0)

  person_int <- idx$person
  log_w <- log(quad$weights)
  n_nodes <- length(quad$nodes)
  score_k <- idx$score_k

  if (config$model == "RSM") {
    k_cat <- length(step_cum)
    step_cum_row <- matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
    obs_idx <- cbind(seq_len(n), score_k + 1L)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      eta_mat <- outer(eta_q, 0:(k_cat - 1))
      log_num <- eta_mat - step_cum_row
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  } else {
    k_cat <- ncol(step_cum)
    crit <- idx$step_idx
    obs_idx <- cbind(seq_len(n), score_k + 1L)
    step_cum_obs <- step_cum[crit, , drop = FALSE]
    k_vals <- 0:(k_cat - 1)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      log_num <- outer(eta_q, k_vals) - step_cum_obs
      rm <- log_num[cbind(seq_len(n), max.col(log_num))]
      ld <- rm + log(rowSums(exp(log_num - rm)))
      lp <- log_num[obs_idx] - ld
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  }

  ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)
  log_w_mat <- matrix(log_w, nrow = nrow(ll_by_person), ncol = n_nodes, byrow = TRUE)
  combined <- log_w_mat + ll_by_person
  row_max <- combined[cbind(seq_len(nrow(combined)), max.col(combined))]
  ll_person <- row_max + log(rowSums(exp(combined - row_max)))
  -sum(ll_person)
}

# Cached version of JMLE gradient (uses pre-computed params/eta/step_cum)
mfrm_grad_jmle_cached <- function(cache, idx, config, sizes) {
  mfrm_grad_jmle_core(cache$params(), cache$eta(), idx, config, sizes)
}

# Cached version of MML gradient (uses pre-computed params/base_eta/step_cum)
mfrm_grad_mml_cached <- function(cache, idx, config, sizes, quad) {
  mfrm_grad_mml_core(cache$params(), cache$base_eta(), idx, config, sizes, quad)
}

# ---- Analytical gradient for JMLE ----
# Returns gradient of -LL (negative log-likelihood) w.r.t. free parameters.
# Key derivation:
#   d(LL)/d(eta_i)      = k_i - E[k|eta_i]  (observed minus expected score)
#   d(LL)/d(delta_s)    = P(X >= s) - I(X >= s)  (for step parameters)
# Chain rule maps these through constraint Jacobians to free parameters.
mfrm_grad_jmle <- function(par, idx, config, sizes) {
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)
  mfrm_grad_jmle_core(params, eta, idx, config, sizes)
}

mfrm_grad_jmle_core <- function(params, eta, idx, config, sizes) {
  n <- length(eta)
  score_k <- idx$score_k
  weight <- idx$weight

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
    k_cat <- ncol(probs)
    n_steps <- k_cat - 1

    # Residuals: observed - expected
    expected <- as.vector(probs %*% (0:(k_cat - 1)))
    residual <- score_k - expected
    if (!is.null(weight)) residual <- residual * weight

    # Gradient w.r.t. expanded theta
    n_theta <- length(params$theta)
    grad_theta_exp <- numeric(n_theta)
    if (n_theta > 0) {
      rs <- rowsum(matrix(residual, ncol = 1), idx$person)
      grad_theta_exp[as.integer(rownames(rs))] <- as.vector(rs)
    }

    # Gradient w.r.t. expanded facets
    grad_facets_exp <- list()
    for (facet in config$facet_names) {
      sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
      n_lev <- length(params$facets[[facet]])
      g <- numeric(n_lev)
      rs <- rowsum(matrix(sign_f * residual, ncol = 1), idx$facets[[facet]])
      g[as.integer(rownames(rs))] <- as.vector(rs)
      grad_facets_exp[[facet]] <- g
    }

    # Gradient w.r.t. centered step parameters
    P_geq <- compute_P_geq(probs)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_resid <- P_geq - I_geq
    if (!is.null(weight)) step_resid <- step_resid * weight
    grad_step_centered <- colSums(step_resid)

    # Map to free step params (centering: centered = raw - mean(raw))
    grad_step_free <- grad_step_centered - mean(grad_step_centered)

  } else {
    # PCM
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                               criterion_splits = idx$criterion_splits)
    k_cat <- ncol(probs)
    n_steps <- k_cat - 1
    n_criteria <- nrow(params$steps_mat)

    expected <- as.vector(probs %*% (0:(k_cat - 1)))
    residual <- score_k - expected
    if (!is.null(weight)) residual <- residual * weight

    n_theta <- length(params$theta)
    grad_theta_exp <- numeric(n_theta)
    if (n_theta > 0) {
      rs <- rowsum(matrix(residual, ncol = 1), idx$person)
      grad_theta_exp[as.integer(rownames(rs))] <- as.vector(rs)
    }

    grad_facets_exp <- list()
    for (facet in config$facet_names) {
      sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
      n_lev <- length(params$facets[[facet]])
      g <- numeric(n_lev)
      rs <- rowsum(matrix(sign_f * residual, ncol = 1), idx$facets[[facet]])
      g[as.integer(rownames(rs))] <- as.vector(rs)
      grad_facets_exp[[facet]] <- g
    }

    # Step gradients per criterion
    P_geq <- compute_P_geq(probs)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_resid <- P_geq - I_geq
    if (!is.null(weight)) step_resid <- step_resid * weight

    # Vectorized step gradient aggregation via rowsum
    grad_step_mat <- matrix(0, n_criteria, n_steps)
    rs_step <- rowsum(step_resid, idx$step_idx)
    rs_ids <- as.integer(rownames(rs_step))
    grad_step_mat[rs_ids, ] <- rs_step

    # Center each criterion row, then flatten to row-major vector
    grad_step_mat_free <- grad_step_mat - rowMeans(grad_step_mat)
    grad_step_free <- as.vector(t(grad_step_mat_free))
  }

  # Project through constraints
  grad_theta_free <- constraint_grad_project(grad_theta_exp, config$theta_spec)
  grad_facet_free <- unlist(lapply(config$facet_names, function(f) {
    constraint_grad_project(grad_facets_exp[[f]], config$facet_specs[[f]])
  }))

  # Return negative gradient (minimizing -LL)
  -c(grad_theta_free, grad_facet_free, grad_step_free)
}

# Marginal Maximum Likelihood (MML) negative log-likelihood.
# Integrates over the person ability distribution using Gauss-Hermite quadrature:
#   L = prod_p integral P(X_p | theta, facets) phi(theta) d(theta)
#     ~ prod_p sum_q w_q P(X_p | theta_q, facets)
# Person parameters are not estimated directly; instead, EAP post-hoc.
mfrm_loglik_mml <- function(par, idx, config, sizes, quad) {
  params <- expand_params(par, sizes, config)
  n <- length(idx$score_k)
  if (n == 0) return(0)

  base_eta <- compute_base_eta(idx, params, config)
  person_int <- idx$person
  log_w <- log(quad$weights)
  n_nodes <- length(quad$nodes)
  score_k <- idx$score_k

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    k_cat <- length(step_cum)
    step_cum_row <- matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
    obs_col <- score_k + 1L
    obs_idx <- cbind(seq_len(n), obs_col)

    # Vectorized: Q passes over all observations, then per-person aggregation
    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      eta_mat <- outer(eta_q, 0:(k_cat - 1))
      log_num <- eta_mat - step_cum_row
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }

    # Per-person sum via rowsum (highly optimized)
    ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)

    # Vectorized logsumexp across quadrature nodes
    log_w_mat <- matrix(log_w, nrow = nrow(ll_by_person), ncol = n_nodes,
                        byrow = TRUE)
    combined <- log_w_mat + ll_by_person
    row_max <- combined[cbind(seq_len(nrow(combined)),
                              max.col(combined))]
    ll_person <- row_max + log(rowSums(exp(combined - row_max)))

  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    k_cat <- ncol(step_cum_mat)
    crit <- idx$step_idx
    obs_col <- score_k + 1L
    obs_idx <- cbind(seq_len(n), obs_col)

    # Vectorized: hoist per-observation step parameters outside quadrature loop
    step_cum_obs <- step_cum_mat[crit, , drop = FALSE]
    k_vals <- 0:(k_cat - 1)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      log_num <- outer(eta_q, k_vals) - step_cum_obs
      rm <- log_num[cbind(seq_len(n), max.col(log_num))]
      ld <- rm + log(rowSums(exp(log_num - rm)))
      lp <- log_num[obs_idx] - ld
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }

    ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)
    log_w_mat <- matrix(log_w, nrow = nrow(ll_by_person), ncol = n_nodes,
                        byrow = TRUE)
    combined <- log_w_mat + ll_by_person
    row_max <- combined[cbind(seq_len(nrow(combined)),
                              max.col(combined))]
    ll_person <- row_max + log(rowSums(exp(combined - row_max)))
  }

  -sum(ll_person)
}

# ---- Analytical gradient for MML ----
# Returns gradient of -LL w.r.t. free parameters under marginal ML.
# Key formula:
#   d(LL)/d(param) = sum_p sum_q posterior_pq * d(ll_pq)/d(param)
# where posterior_pq = w_q * L_p(theta_q) / sum_q' w_q' * L_p(theta_q')
# and d(ll_pq)/d(param) uses the same residual/step derivatives as JMLE
# but evaluated at each quadrature node theta_q.
mfrm_grad_mml <- function(par, idx, config, sizes, quad) {
  params <- expand_params(par, sizes, config)
  base_eta <- compute_base_eta(idx, params, config)
  mfrm_grad_mml_core(params, base_eta, idx, config, sizes, quad)
}

mfrm_grad_mml_core <- function(params, base_eta, idx, config, sizes, quad) {
  n <- length(idx$score_k)
  if (n == 0) return(rep(0, sum(unlist(sizes))))
  person_int <- idx$person
  log_w <- log(quad$weights)
  n_nodes <- length(quad$nodes)
  score_k <- idx$score_k
  weight <- idx$weight
  n_persons <- max(person_int)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    k_cat <- length(step_cum)
    n_steps <- k_cat - 1

    # Phase 1: Compute log-prob matrix and category probs at each node
    log_prob_mat <- matrix(0, n, n_nodes)
    prob_list <- vector("list", n_nodes)

    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      probs_q <- category_prob_rsm(eta_q, step_cum)
      log_p <- log(pmax(probs_q[cbind(seq_len(n), score_k + 1L)], 1e-300))
      if (!is.null(weight)) log_p <- log_p * weight
      log_prob_mat[, q] <- log_p
      prob_list[[q]] <- probs_q
    }

    # Phase 2: Compute posterior weights
    ll_by_person <- rowsum(log_prob_mat, person_int)
    person_ids <- as.integer(rownames(ll_by_person))
    log_w_mat <- matrix(log_w, nrow = nrow(ll_by_person), ncol = n_nodes, byrow = TRUE)
    combined <- log_w_mat + ll_by_person
    row_max <- row_max_fast(combined)
    log_posterior <- combined - (row_max + log(rowSums(exp(combined - row_max))))
    posterior <- exp(log_posterior)

    # Map person_int to posterior row indices
    person_to_row <- integer(n_persons)
    person_to_row[person_ids] <- seq_along(person_ids)
    obs_person_row <- person_to_row[person_int]

    # Pre-compute I(X >= s) (doesn't depend on quadrature node)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0

    # Phase 3: Accumulate gradients
    grad_facets_exp <- lapply(config$facet_names, function(f) numeric(length(params$facets[[f]])))
    names(grad_facets_exp) <- config$facet_names
    grad_step_centered <- numeric(n_steps)

    for (q in seq_len(n_nodes)) {
      probs_q <- prob_list[[q]]
      expected_q <- as.vector(probs_q %*% (0:(k_cat - 1)))
      residual_q <- score_k - expected_q
      if (!is.null(weight)) residual_q <- residual_q * weight

      obs_post_q <- posterior[obs_person_row, q]
      w_residual <- residual_q * obs_post_q

      # Facet gradients
      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]])
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }

      # Step gradients
      P_geq <- compute_P_geq(probs_q)
      step_resid <- P_geq - I_geq
      if (!is.null(weight)) step_resid <- step_resid * weight
      step_resid_w <- step_resid * obs_post_q
      grad_step_centered <- grad_step_centered + colSums(step_resid_w)
    }

    grad_step_free <- grad_step_centered - mean(grad_step_centered)

  } else {
    # PCM case
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    k_cat <- ncol(step_cum_mat)
    n_steps <- k_cat - 1
    n_criteria <- nrow(params$steps_mat)

    log_prob_mat <- matrix(0, n, n_nodes)
    prob_list <- vector("list", n_nodes)

    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      probs_q <- category_prob_pcm(eta_q, step_cum_mat, idx$step_idx,
                                   criterion_splits = idx$criterion_splits)
      log_p <- log(pmax(probs_q[cbind(seq_len(n), score_k + 1L)], 1e-300))
      if (!is.null(weight)) log_p <- log_p * weight
      log_prob_mat[, q] <- log_p
      prob_list[[q]] <- probs_q
    }

    ll_by_person <- rowsum(log_prob_mat, person_int)
    person_ids <- as.integer(rownames(ll_by_person))
    log_w_mat <- matrix(log_w, nrow = nrow(ll_by_person), ncol = n_nodes, byrow = TRUE)
    combined <- log_w_mat + ll_by_person
    row_max <- row_max_fast(combined)
    log_posterior <- combined - (row_max + log(rowSums(exp(combined - row_max))))
    posterior <- exp(log_posterior)

    person_to_row <- integer(n_persons)
    person_to_row[person_ids] <- seq_along(person_ids)
    obs_person_row <- person_to_row[person_int]

    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    step_idx_int <- idx$step_idx

    grad_facets_exp <- lapply(config$facet_names, function(f) numeric(length(params$facets[[f]])))
    names(grad_facets_exp) <- config$facet_names
    grad_step_mat <- matrix(0, n_criteria, n_steps)

    for (q in seq_len(n_nodes)) {
      probs_q <- prob_list[[q]]
      expected_q <- as.vector(probs_q %*% (0:(k_cat - 1)))
      residual_q <- score_k - expected_q
      if (!is.null(weight)) residual_q <- residual_q * weight

      obs_post_q <- posterior[obs_person_row, q]
      w_residual <- residual_q * obs_post_q

      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]])
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }

      P_geq <- compute_P_geq(probs_q)
      step_resid <- P_geq - I_geq
      if (!is.null(weight)) step_resid <- step_resid * weight
      step_resid_w <- step_resid * obs_post_q

      # Vectorized step gradient aggregation via rowsum
      rs_step <- rowsum(step_resid_w, step_idx_int)
      rs_ids <- as.integer(rownames(rs_step))
      grad_step_mat[rs_ids, ] <- grad_step_mat[rs_ids, ] + rs_step
    }

    grad_step_mat_free <- grad_step_mat - rowMeans(grad_step_mat)
    grad_step_free <- as.vector(t(grad_step_mat_free))
  }

  # Project facet gradients through constraints
  grad_facet_free <- unlist(lapply(config$facet_names, function(f) {
    constraint_grad_project(grad_facets_exp[[f]], config$facet_specs[[f]])
  }))

  # Return negative gradient (minimizing -LL); no theta in MML
  -c(grad_facet_free, grad_step_free)
}

# Expected A Posteriori (EAP) person ability estimates under MML.
# For each person p:
#   theta_EAP = sum_q theta_q * w_q * L(X_p|theta_q) / sum_q w_q * L(X_p|theta_q)
#   SD_EAP   = sqrt(E[theta^2] - (E[theta])^2)  (posterior standard deviation)
compute_person_eap <- function(idx, config, params, quad) {
  n <- length(idx$score_k)
  if (n == 0) {
    return(tibble(Person = character(0), Estimate = numeric(0), SD = numeric(0)))
  }
  base_eta <- compute_base_eta(idx, params, config)
  person_int <- idx$person
  log_w <- log(quad$weights)
  n_nodes <- length(quad$nodes)
  score_k <- idx$score_k

  # Build per-obs log-prob matrix (n x n_nodes), same as mfrm_loglik_mml
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    k_cat <- length(step_cum)
    step_cum_row <- matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
    obs_idx <- cbind(seq_len(n), score_k + 1L)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      eta_mat <- outer(eta_q, 0:(k_cat - 1))
      log_num <- eta_mat - step_cum_row
      row_max <- log_num[cbind(seq_len(n), max.col(log_num))]
      log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
      lp <- log_num[obs_idx] - log_denom
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    k_cat <- ncol(step_cum_mat)
    crit <- idx$step_idx
    obs_col <- score_k + 1L
    obs_idx <- cbind(seq_len(n), obs_col)

    # Vectorized: hoist per-observation step parameters outside quadrature loop
    step_cum_obs <- step_cum_mat[crit, , drop = FALSE]
    k_vals <- 0:(k_cat - 1)

    log_prob_mat <- matrix(0, n, n_nodes)
    for (q in seq_len(n_nodes)) {
      eta_q <- base_eta + quad$nodes[q]
      log_num <- outer(eta_q, k_vals) - step_cum_obs
      rm <- log_num[cbind(seq_len(n), max.col(log_num))]
      ld <- rm + log(rowSums(exp(log_num - rm)))
      lp <- log_num[obs_idx] - ld
      if (!is.null(idx$weight)) lp <- lp * idx$weight
      log_prob_mat[, q] <- lp
    }
  }

  # Per-person sum of log-probs via rowsum
  ll_by_person <- rowsum(log_prob_mat, person_int, reorder = FALSE)
  n_persons <- nrow(ll_by_person)

  # Posterior weights: log_w + ll_nodes, normalized per person
  log_w_mat <- matrix(log_w, nrow = n_persons, ncol = n_nodes, byrow = TRUE)
  log_post <- log_w_mat + ll_by_person
  # Normalize each row
  row_max <- log_post[cbind(seq_len(n_persons), max.col(log_post))]
  log_norm <- row_max + log(rowSums(exp(log_post - row_max)))
  log_post <- log_post - log_norm
  post_w <- exp(log_post)

  # EAP and SD
  nodes_mat <- matrix(quad$nodes, nrow = n_persons, ncol = n_nodes, byrow = TRUE)
  eap <- rowSums(nodes_mat * post_w)
  sd_eap <- sqrt(rowSums((nodes_mat - eap)^2 * post_w))

  tibble(Estimate = eap, SD = sd_eap)
}

prepare_constraint_specs <- function(prep,
                                     anchor_df = NULL,
                                     group_anchor_df = NULL,
                                     noncenter_facet = "Person",
                                     dummy_facets = character(0)) {
  facet_names <- prep$facet_names
  all_facets <- c("Person", facet_names)
  anchor_audit <- audit_anchor_tables(
    prep = prep,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  anchor_df <- anchor_audit$anchors
  group_anchor_df <- anchor_audit$group_anchors

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
        df <- filter(anchor_df, Facet == facet)
        if (nrow(df) > 0) {
          anchors <- setNames(df$Anchor, df$Level)
          anchors <- anchors[names(anchors) %in% levels]
          if (length(anchors) > 0) anchor_map[[facet]] <- anchors
        }
      }

      if (nrow(group_anchor_df) > 0) {
        df <- filter(group_anchor_df, Facet == facet)
        if (nrow(df) > 0) {
          groups <- setNames(df$Group, df$Level)
          groups <- groups[names(groups) %in% levels]
          if (length(groups) > 0) group_map[[facet]] <- groups

          group_vals <- df |>
            select(Group, GroupValue) |>
            distinct() |>
            filter(!is.na(Group)) |>
            group_by(Group) |>
            summarize(
              GroupValue = {
                vals <- GroupValue[!is.na(GroupValue)]
                if (length(vals) == 0) NA_real_ else vals[1]
              },
              .groups = "drop"
            ) |>
            mutate(GroupValue = ifelse(is.na(GroupValue), 0, GroupValue))
          if (nrow(group_vals) > 0) {
            group_values[[facet]] <- setNames(group_vals$GroupValue, group_vals$Group)
          }
        }
      }
    }
  }

  theta_spec <- build_facet_constraint(
    levels = prep$levels$Person,
    anchors = anchor_map$Person,
    groups = group_map$Person,
    group_values = group_values$Person,
    centered = !(identical(noncenter_facet, "Person"))
  )

  facet_specs <- lapply(facet_names, function(facet) {
    build_facet_constraint(
      levels = prep$levels[[facet]],
      anchors = anchor_map[[facet]],
      groups = group_map[[facet]],
      group_values = group_values[[facet]],
      centered = !identical(noncenter_facet, facet)
    )
  })
  names(facet_specs) <- facet_names

  anchor_summary <- tibble(
    Facet = all_facets,
    AnchoredLevels = vapply(anchor_map, function(x) length(x), integer(1)),
    GroupAnchors = vapply(group_map, function(x) length(unique(x)), integer(1)),
    DummyFacet = all_facets %in% dummy_facets
  )

  list(
    theta_spec = theta_spec,
    facet_specs = facet_specs,
    anchor_summary = anchor_summary,
    anchor_audit = anchor_audit
  )
}

read_flexible_table <- function(text_value, file_input) {
  if (!is.null(file_input) && !is.null(file_input$datapath)) {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", file_input$name, ignore.case = TRUE), "\t", ",")
    return(read.csv(file_input$datapath,
                    sep = sep,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE))
  }
  if (is.null(text_value)) return(tibble())
  text_value <- trimws(text_value)
  if (!nzchar(text_value)) return(tibble())
  sep <- if (grepl("\t", text_value)) {
    "\t"
  } else if (grepl(";", text_value)) {
    ";"
  } else {
    ","
  }
  read.csv(text = text_value,
           sep = sep,
           header = TRUE,
           stringsAsFactors = FALSE,
           check.names = FALSE)
}

normalize_anchor_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Anchor = numeric(0)))
  }
  nm <- tolower(names(df))
  facet_col <- which(nm %in% c("facet", "facets"))
  level_col <- which(nm %in% c("level", "element", "label"))
  anchor_col <- which(nm %in% c("anchor", "value", "measure"))
  if (length(facet_col) == 0 || length(level_col) == 0 || length(anchor_col) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Anchor = numeric(0)))
  }
  tibble(
    Facet = as.character(df[[facet_col[1]]]),
    Level = as.character(df[[level_col[1]]]),
    Anchor = suppressWarnings(as.numeric(df[[anchor_col[1]]]))
  ) |>
    filter(!is.na(Facet), !is.na(Level))
}

normalize_group_anchor_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Group = character(0), GroupValue = numeric(0)))
  }
  nm <- tolower(names(df))
  facet_col <- which(nm %in% c("facet", "facets"))
  level_col <- which(nm %in% c("level", "element", "label"))
  group_col <- which(nm %in% c("group", "subset"))
  value_col <- which(nm %in% c("groupvalue", "value", "anchor"))
  if (length(facet_col) == 0 || length(level_col) == 0 || length(group_col) == 0 || length(value_col) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Group = character(0), GroupValue = numeric(0)))
  }
  tibble(
    Facet = as.character(df[[facet_col[1]]]),
    Level = as.character(df[[level_col[1]]]),
    Group = as.character(df[[group_col[1]]]),
    GroupValue = suppressWarnings(as.numeric(df[[value_col[1]]]))
  ) |>
    filter(!is.na(Facet), !is.na(Level), !is.na(Group))
}

collect_anchor_levels <- function(prep) {
  facets <- c("Person", prep$facet_names)
  rows <- lapply(facets, function(facet) {
    lv <- prep$levels[[facet]]
    if (is.null(lv) || length(lv) == 0) return(tibble())
    tibble(
      Facet = as.character(facet),
      Level = as.character(lv)
    )
  })
  bind_rows(rows)
}

safe_join_key <- function(facet, level) {
  paste0(as.character(facet), "\r", as.character(level))
}

build_anchor_issue_counts <- function(issue_tables) {
  issue_names <- names(issue_tables)
  tibble(
    Issue = issue_names,
    N = vapply(issue_tables, nrow, integer(1))
  )
}

build_anchor_recommendations <- function(facet_summary,
                                         issue_counts,
                                         design_checks = NULL,
                                         min_common_anchors = 5L,
                                         min_obs_per_element = 30,
                                         min_obs_per_category = 10,
                                         noncenter_facet = "Person",
                                         dummy_facets = character(0)) {
  rec <- character(0)

  if (!is.null(issue_counts) && nrow(issue_counts) > 0) {
    n_overlap <- issue_counts$N[issue_counts$Issue == "overlap_anchor_group"]
    n_missing <- issue_counts$N[issue_counts$Issue == "missing_group_values"]
    n_dup_anchor <- issue_counts$N[issue_counts$Issue == "duplicate_anchors"]
    n_dup_group <- issue_counts$N[issue_counts$Issue == "duplicate_group_assignments"]
    n_group_conf <- issue_counts$N[issue_counts$Issue == "group_value_conflicts"]

    if (length(n_overlap) > 0 && n_overlap > 0) {
      rec <- c(rec, "Levels listed in both anchor and group-anchor tables are directly anchored (fixed anchors take precedence).")
    }
    if (length(n_missing) > 0 && n_missing > 0) {
      rec <- c(rec, "Some group anchors had missing GroupValue; default 0 was applied using the legacy-compatible group-centering rule.")
    }
    if (length(n_dup_anchor) > 0 && n_dup_anchor > 0) {
      rec <- c(rec, "Duplicate anchors were detected; the last row per Facet-Level was retained.")
    }
    if (length(n_dup_group) > 0 && n_dup_group > 0) {
      rec <- c(rec, "A Facet-Level appeared in multiple groups; the last row per Facet-Level was retained.")
    }
    if (length(n_group_conf) > 0 && n_group_conf > 0) {
      rec <- c(rec, "Conflicting GroupValue settings were detected within the same Facet-Group; the most recent finite value was retained.")
    }
  }

  if (!is.null(facet_summary) && nrow(facet_summary) > 0) {
    link_tbl <- facet_summary |>
      filter(Facet != "Person", AnchoredLevels > 0, AnchoredLevels < min_common_anchors)
    if (nrow(link_tbl) > 0) {
      rec <- c(
        rec,
        paste0(
          "FACETS linking guideline: consider >= ", min_common_anchors,
          " common anchor levels per linking facet. Low-anchor facets: ",
          paste(link_tbl$Facet, collapse = ", "), "."
        )
      )
    }

    fixed_tbl <- facet_summary |>
      filter(FreeLevels <= 0, !Facet %in% dummy_facets)
    if (nrow(fixed_tbl) > 0) {
      rec <- c(
        rec,
        paste0(
          "Some facets are fully constrained (no free levels): ",
          paste(fixed_tbl$Facet, collapse = ", "),
          ". Verify this is intentional."
        )
      )
    }
  }

  if (!is.null(design_checks) && is.list(design_checks)) {
    if (!is.null(design_checks$low_observation_levels) && nrow(design_checks$low_observation_levels) > 0) {
      low_facets <- unique(design_checks$low_observation_levels$Facet)
      rec <- c(
        rec,
        paste0(
          "Linacre guideline: about ", fmt_count(min_obs_per_element),
          " observations per element are desirable. Low-observation facets: ",
          paste(low_facets, collapse = ", "), "."
        )
      )
    }
    if (!is.null(design_checks$low_categories) && nrow(design_checks$low_categories) > 0) {
      cats <- paste(design_checks$low_categories$Category, collapse = ", ")
      rec <- c(
        rec,
        paste0(
          "Linacre guideline: about ", fmt_count(min_obs_per_category),
          " observations per rating category are desirable. Low categories: ", cats, "."
        )
      )
    }
  }

  rec <- c(
    rec,
    "For linked analyses, keep Umean/Uscale from the source calibration so reporting origin and scaling stay consistent.",
    paste0("Current noncenter facet is '", noncenter_facet, "'. Other facets are centered unless constrained by anchors/group anchors.")
  )

  unique(rec)
}

audit_anchor_tables <- function(prep,
                                anchor_df = NULL,
                                group_anchor_df = NULL,
                                min_common_anchors = 5L,
                                min_obs_per_element = 30,
                                min_obs_per_category = 10,
                                noncenter_facet = "Person",
                                dummy_facets = character(0)) {
  all_facets <- c("Person", prep$facet_names)
  level_df <- collect_anchor_levels(prep)
  valid_keys <- safe_join_key(level_df$Facet, level_df$Level)
  anchor_schema_issue <- if (!is.null(anchor_df) && nrow(anchor_df) > 0) {
    nm <- tolower(names(anchor_df))
    has_facet <- any(nm %in% c("facet", "facets"))
    has_level <- any(nm %in% c("level", "element", "label"))
    has_anchor <- any(nm %in% c("anchor", "value", "measure"))
    if (!(has_facet && has_level && has_anchor)) {
      tibble::tibble(
        Columns = paste(names(anchor_df), collapse = ", "),
        Required = "facet + level + anchor/value/measure"
      )
    } else {
      tibble::tibble()
    }
  } else {
    tibble::tibble()
  }
  group_schema_issue <- if (!is.null(group_anchor_df) && nrow(group_anchor_df) > 0) {
    nm <- tolower(names(group_anchor_df))
    has_facet <- any(nm %in% c("facet", "facets"))
    has_level <- any(nm %in% c("level", "element", "label"))
    has_group <- any(nm %in% c("group", "subset"))
    has_value <- any(nm %in% c("groupvalue", "value", "anchor"))
    if (!(has_facet && has_level && has_group && has_value)) {
      tibble::tibble(
        Columns = paste(names(group_anchor_df), collapse = ", "),
        Required = "facet + level + group/subset + groupvalue/value/anchor"
      )
    } else {
      tibble::tibble()
    }
  } else {
    tibble::tibble()
  }

  anchor_in <- normalize_anchor_df(anchor_df) |>
    mutate(
      .Row = row_number(),
      Facet = trimws(as.character(Facet)),
      Level = trimws(as.character(Level)),
      .Key = safe_join_key(Facet, Level),
      .ValidFacet = Facet %in% all_facets,
      .ValidLevel = .Key %in% valid_keys,
      .ValidValue = is.finite(Anchor)
    )

  group_in <- normalize_group_anchor_df(group_anchor_df) |>
    mutate(
      .Row = row_number(),
      Facet = trimws(as.character(Facet)),
      Level = trimws(as.character(Level)),
      Group = trimws(as.character(Group)),
      .Key = safe_join_key(Facet, Level),
      .ValidFacet = Facet %in% all_facets,
      .ValidLevel = .Key %in% valid_keys,
      .ValidGroup = nzchar(Group),
      .FiniteGroupValue = is.finite(GroupValue)
    )

  issues <- list(
    anchor_schema_mismatch = anchor_schema_issue,
    group_anchor_schema_mismatch = group_schema_issue,
    unknown_anchor_facets = anchor_in |>
      filter(!.ValidFacet) |>
      select(Facet, Level, Anchor),
    unknown_anchor_levels = anchor_in |>
      filter(.ValidFacet, !.ValidLevel) |>
      select(Facet, Level, Anchor),
    invalid_anchor_values = anchor_in |>
      filter(.ValidFacet, .ValidLevel, !.ValidValue) |>
      select(Facet, Level, Anchor),
    duplicate_anchors = anchor_in |>
      filter(.ValidFacet, .ValidLevel, .ValidValue) |>
      group_by(Facet, Level) |>
      summarize(
        Rows = n(),
        DistinctValues = n_distinct(Anchor),
        Values = paste(unique(round(Anchor, 6)), collapse = ", "),
        .groups = "drop"
      ) |>
      filter(Rows > 1 | DistinctValues > 1),
    unknown_group_facets = group_in |>
      filter(!.ValidFacet) |>
      select(Facet, Level, Group, GroupValue),
    unknown_group_levels = group_in |>
      filter(.ValidFacet, !.ValidLevel) |>
      select(Facet, Level, Group, GroupValue),
    invalid_group_labels = group_in |>
      filter(.ValidFacet, .ValidLevel, !.ValidGroup) |>
      select(Facet, Level, Group, GroupValue),
    duplicate_group_assignments = group_in |>
      filter(.ValidFacet, .ValidLevel, .ValidGroup) |>
      group_by(Facet, Level) |>
      summarize(
        Rows = n(),
        DistinctGroups = n_distinct(Group),
        Groups = paste(unique(Group), collapse = ", "),
        .groups = "drop"
      ) |>
      filter(Rows > 1 | DistinctGroups > 1)
  )

  anchors_clean <- anchor_in |>
    filter(.ValidFacet, .ValidLevel, .ValidValue) |>
    arrange(.Row) |>
    group_by(Facet, Level) |>
    slice_tail(n = 1) |>
    ungroup() |>
    select(Facet, Level, Anchor)

  groups_clean <- group_in |>
    filter(.ValidFacet, .ValidLevel, .ValidGroup) |>
    arrange(.Row) |>
    group_by(Facet, Level) |>
    slice_tail(n = 1) |>
    ungroup() |>
    select(Facet, Level, Group, GroupValue, .FiniteGroupValue)

  group_value_tbl <- groups_clean |>
    arrange(Facet, Group) |>
    group_by(Facet, Group) |>
    summarize(
      .NFinite = sum(.FiniteGroupValue, na.rm = TRUE),
      ChosenGroupValue = if (any(.FiniteGroupValue)) dplyr::last(GroupValue[.FiniteGroupValue]) else 0,
      DistinctFiniteValues = n_distinct(GroupValue[.FiniteGroupValue]),
      FiniteValues = paste(unique(round(GroupValue[.FiniteGroupValue], 6)), collapse = ", "),
      .groups = "drop"
    )

  issues$missing_group_values <- group_value_tbl |>
    filter(.NFinite == 0) |>
    select(Facet, Group)

  issues$group_value_conflicts <- group_value_tbl |>
    filter(DistinctFiniteValues > 1) |>
    select(Facet, Group, DistinctFiniteValues, FiniteValues)

  groups_clean <- groups_clean |>
    select(Facet, Level, Group) |>
    left_join(group_value_tbl |> select(Facet, Group, ChosenGroupValue), by = c("Facet", "Group")) |>
    rename(GroupValue = ChosenGroupValue) |>
    mutate(GroupValue = ifelse(is.finite(GroupValue), GroupValue, 0))

  overlap_tbl <- inner_join(
    anchors_clean |> select(Facet, Level),
    groups_clean |> select(Facet, Level),
    by = c("Facet", "Level")
  )
  issues$overlap_anchor_group <- overlap_tbl

  constrained_counts <- bind_rows(
    anchors_clean |> select(Facet, Level),
    groups_clean |> select(Facet, Level)
  ) |>
    distinct(Facet, Level) |>
    group_by(Facet) |>
    summarize(ConstrainedLevels = n_distinct(Level), .groups = "drop")

  facet_counts <- level_df |>
    group_by(Facet) |>
    summarize(Levels = n_distinct(Level), .groups = "drop")

  anchor_counts <- anchors_clean |>
    group_by(Facet) |>
    summarize(AnchoredLevels = n_distinct(Level), .groups = "drop")

  group_counts <- groups_clean |>
    group_by(Facet) |>
    summarize(GroupedLevels = n_distinct(Level), GroupCount = n_distinct(Group), .groups = "drop")

  overlap_counts <- overlap_tbl |>
    group_by(Facet) |>
    summarize(OverlapLevels = n_distinct(Level), .groups = "drop")

  facet_summary <- facet_counts |>
    left_join(anchor_counts, by = "Facet") |>
    left_join(group_counts, by = "Facet") |>
    left_join(constrained_counts, by = "Facet") |>
    left_join(overlap_counts, by = "Facet") |>
    mutate(
      AnchoredLevels = tidyr::replace_na(AnchoredLevels, 0L),
      GroupedLevels = tidyr::replace_na(GroupedLevels, 0L),
      GroupCount = tidyr::replace_na(GroupCount, 0L),
      ConstrainedLevels = tidyr::replace_na(ConstrainedLevels, 0L),
      OverlapLevels = tidyr::replace_na(OverlapLevels, 0L),
      FreeLevels = pmax(Levels - ConstrainedLevels, 0L),
      Noncenter = Facet == noncenter_facet,
      DummyFacet = Facet %in% dummy_facets
    ) |>
    arrange(match(Facet, all_facets))

  design_df <- prep$data |>
    mutate(
      Person = as.character(Person),
      Score = as.numeric(Score),
      Weight = as.numeric(Weight),
      Weight = ifelse(is.finite(Weight) & Weight > 0, Weight, 0)
    )

  level_obs <- bind_rows(lapply(all_facets, function(facet) {
    if (!facet %in% names(design_df)) return(tibble())
    design_df |>
      mutate(Level = as.character(.data[[facet]])) |>
      group_by(Facet = facet, Level) |>
      summarize(
        RawN = n(),
        WeightedN = sum(Weight, na.rm = TRUE),
        .groups = "drop"
      )
  }))

  level_obs_summary <- if (nrow(level_obs) == 0) {
    tibble()
  } else {
    level_obs |>
      group_by(Facet) |>
      summarize(
        Levels = n(),
        MinObsPerLevel = min(WeightedN, na.rm = TRUE),
        MedianObsPerLevel = stats::median(WeightedN, na.rm = TRUE),
        RecommendedMinObs = as.numeric(min_obs_per_element),
        PassMinObs = all(WeightedN >= min_obs_per_element),
        .groups = "drop"
      )
  }

  low_observation_levels <- if (nrow(level_obs) == 0) {
    tibble()
  } else {
    level_obs |>
      filter(WeightedN < min_obs_per_element) |>
      arrange(Facet, WeightedN, Level)
  }

  category_counts <- design_df |>
    group_by(Category = Score) |>
    summarize(
      RawN = n(),
      WeightedN = sum(Weight, na.rm = TRUE),
      RecommendedMinObs = as.numeric(min_obs_per_category),
      PassMinObs = WeightedN >= min_obs_per_category,
      .groups = "drop"
    ) |>
    arrange(Category)

  low_categories <- category_counts |>
    filter(!PassMinObs)

  design_checks <- list(
    level_observation_summary = level_obs_summary,
    low_observation_levels = low_observation_levels,
    category_counts = category_counts,
    low_categories = low_categories
  )

  issue_counts <- build_anchor_issue_counts(issues)
  rec <- build_anchor_recommendations(
    facet_summary = facet_summary,
    issue_counts = issue_counts,
    design_checks = design_checks,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  thresholds <- list(
    min_common_anchors = as.integer(min_common_anchors),
    min_obs_per_element = as.numeric(min_obs_per_element),
    min_obs_per_category = as.numeric(min_obs_per_category)
  )

  list(
    anchors = anchors_clean,
    group_anchors = groups_clean,
    facet_summary = facet_summary,
    design_checks = design_checks,
    thresholds = thresholds,
    issues = issues,
    issue_counts = issue_counts,
    recommendations = rec
  )
}

resolve_pcm_step_facet <- function(model, step_facet, facet_names) {
  if (model != "PCM") return(NULL)
  resolved <- if (is.null(step_facet)) facet_names[1] else step_facet
  if (!resolved %in% facet_names) {
    stop("step_facet = '", resolved, "' is not among the declared facets: ",
         paste(facet_names, collapse = ", "), ". ",
         "Supply a valid facet name.", call. = FALSE)
  }
  resolved
}

sanitize_noncenter_facet <- function(noncenter_facet, facet_names) {
  if (!is.null(noncenter_facet) && noncenter_facet %in% c("Person", facet_names)) {
    return(noncenter_facet)
  }
  "Person"
}

sanitize_dummy_facets <- function(dummy_facets, facet_names) {
  if (is.null(dummy_facets)) return(character(0))
  intersect(dummy_facets, c("Person", facet_names))
}

build_facet_signs <- function(facet_names, positive_facets = character(0)) {
  positives <- intersect(positive_facets, facet_names)
  signs <- setNames(ifelse(facet_names %in% positives, 1, -1), facet_names)
  list(signs = signs, positive_facets = positives)
}

audit_interaction_orientation <- function(config, interaction_facets) {
  interaction_facets <- unique(as.character(interaction_facets %||% character(0)))
  interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
  if (length(interaction_facets) == 0) {
    return(list(
      table = tibble(),
      mixed_sign = FALSE,
      direction_note = "No interaction facets were supplied for orientation audit.",
      recommended_action = "None."
    ))
  }

  sign_map <- config$facet_signs %||% list()
  sign_vals <- vapply(interaction_facets, function(facet) {
    sign_map[[facet]] %||% ifelse(identical(facet, "Person"), 1, -1)
  }, numeric(1))
  positive_facets <- as.character(config$positive_facets %||% character(0))
  orientation_tbl <- tibble(
    Facet = interaction_facets,
    Sign = sign_vals,
    Orientation = ifelse(sign_vals >= 0, "positive", "negative"),
    InPositiveFacets = interaction_facets %in% positive_facets
  )

  mixed_sign <- length(unique(sign_vals[is.finite(sign_vals)])) > 1L
  if (mixed_sign) {
    direction_note <- paste(
      "Selected interaction facets mix positive and negative score orientations.",
      "Report cells as higher-than-expected or lower-than-expected rather than using directional labels tied to a single facet."
    )
    recommended_action <- paste(
      "For clearer interpretation, either align the selected facets through `positive_facets`",
      "or keep the mixed orientation and use neutral higher/lower-than-expected wording in reports."
    )
  } else {
    direction_note <- "Selected interaction facets share the same score orientation."
    recommended_action <- "Directional wording is stable across the selected interaction facets."
  }

  list(
    table = orientation_tbl,
    mixed_sign = mixed_sign,
    direction_note = direction_note,
    recommended_action = recommended_action
  )
}

build_estimation_config <- function(prep,
                                    model,
                                    method,
                                    step_facet,
                                    weight_col,
                                    facet_signs,
                                    positive_facets,
                                    noncenter_facet,
                                    dummy_facets,
                                    anchor_df,
                                    group_anchor_df) {
  config <- list(
    model = model,
    method = method,
    n_person = length(prep$levels$Person),
    n_cat = prep$rating_max - prep$rating_min + 1,
    facet_names = prep$facet_names,
    facet_levels = prep$levels[prep$facet_names],
    step_facet = step_facet,
    keep_original = isTRUE(prep$keep_original),
    rating_min = prep$rating_min,
    rating_max = prep$rating_max,
    score_map = prep$score_map
  )
  config$weight_col <- if (!is.null(weight_col)) weight_col else NULL
  config$positive_facets <- positive_facets
  config$facet_signs <- facet_signs

  constraint_specs <- prepare_constraint_specs(
    prep = prep,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  config$theta_spec <- constraint_specs$theta_spec
  config$facet_specs <- constraint_specs$facet_specs
  config$noncenter_facet <- noncenter_facet
  config$dummy_facets <- dummy_facets
  config$anchor_summary <- constraint_specs$anchor_summary
  config$anchor_audit <- constraint_specs$anchor_audit
  config$source_columns <- prep$source_columns

  list(
    config = config,
    sizes = build_param_sizes(config)
  )
}

build_initial_param_vector <- function(config, sizes) {
  n_cat <- config$n_cat
  step_init <- if (n_cat > 1) {
    seq(-1, 1, length.out = n_cat - 1)
  } else {
    numeric(0)
  }

  facet_starts <- unlist(lapply(config$facet_names, function(f) rep(0, sizes[[f]])))
  c(
    rep(0, sizes$theta),
    facet_starts,
    if (config$model == "RSM") {
      step_init
    } else {
      rep(step_init, length(config$facet_levels[[config$step_facet]]))
    }
  )
}

# Parameter cache shared between fn and gr to avoid redundant expand_params /
# compute_eta calls within the same optim iteration (BFGS calls fn then gr
# with identical par).
make_param_cache <- function(sizes, config, idx, is_mml = FALSE) {
  cached_par <- NULL
  cached_params <- NULL
  cached_eta <- NULL      # full eta for JMLE
  cached_base_eta <- NULL  # base_eta for MML
  cached_step_cum <- NULL  # step_cum (RSM) or step_cum_mat (PCM)

  list(
    ensure = function(par) {
      if (identical(par, cached_par)) return(invisible(NULL))
      cached_par <<- par
      cached_params <<- expand_params(par, sizes, config)
      if (is_mml) {
        cached_base_eta <<- compute_base_eta(idx, cached_params, config)
      } else {
        cached_eta <<- compute_eta(idx, cached_params, config)
      }
      if (config$model == "RSM") {
        cached_step_cum <<- c(0, cumsum(cached_params$steps))
      } else {
        cached_step_cum <<- t(apply(cached_params$steps_mat, 1,
                                    function(x) c(0, cumsum(x))))
      }
      invisible(NULL)
    },
    params = function() cached_params,
    eta = function() cached_eta,
    base_eta = function() cached_base_eta,
    step_cum = function() cached_step_cum
  )
}

run_mfrm_optimization <- function(start,
                                  method,
                                  idx,
                                  config,
                                  sizes,
                                  quad_points,
                                  maxit,
                                  reltol,
                                  suppress_convergence_warning = FALSE) {
  control <- list(maxit = maxit, reltol = reltol)
  quad <- NULL

  if (method == "JMLE") {
    cache <- make_param_cache(sizes, config, idx, is_mml = FALSE)
  } else {
    quad <- gauss_hermite_normal(quad_points)
    cache <- make_param_cache(sizes, config, idx, is_mml = TRUE)
  }

  fn <- function(par, idx, config, sizes, quad = NULL) {
    cache$ensure(par)
    if (method == "JMLE") {
      mfrm_loglik_jmle_cached(cache, idx, config)
    } else {
      mfrm_loglik_mml_cached(cache, idx, config, quad)
    }
  }

  gr <- function(par, idx, config, sizes, quad = NULL) {
    cache$ensure(par)
    if (method == "JMLE") {
      mfrm_grad_jmle_cached(cache, idx, config, sizes)
    } else {
      mfrm_grad_mml_cached(cache, idx, config, sizes, quad)
    }
  }

  opt <- tryCatch(
    optim(par = start, fn = fn, gr = gr, method = "BFGS",
          control = control, idx = idx, config = config, sizes = sizes,
          quad = quad),
    error = function(e) {
      stop("Model optimization failed: ", conditionMessage(e), ". ",
           "Possible causes: (1) insufficient data for the number of parameters, ",
           "(2) extreme score distributions, (3) near-constant responses. ",
           "Try reducing facets, increasing maxit, or checking data quality.",
           call. = FALSE)
    }
  )

  if (opt$convergence != 0 && !isTRUE(suppress_convergence_warning)) {
    warning("Optimizer did not fully converge (code = ", opt$convergence, "). ",
            "Consider increasing maxit (current: ", maxit, ") ",
            "or relaxing reltol (current: ", reltol, ").",
            call. = FALSE)
  }

  opt
}

build_person_table <- function(method, idx, config, params, prep, quad_points) {
  if (method == "MML") {
    quad <- gauss_hermite_normal(quad_points)
    return(
      compute_person_eap(idx, config, params, quad) |>
        mutate(Person = prep$levels$Person) |>
        select(Person, Estimate, SD)
    )
  }

  tibble(
    Person = prep$levels$Person,
    Estimate = params$theta
  )
}

build_other_facet_table <- function(config, prep, params) {
  facet_tbls <- lapply(config$facet_names, function(facet) {
    tibble(Level = prep$levels[[facet]], Estimate = params$facets[[facet]]) |>
      mutate(Facet = facet, .before = 1)
  })
  bind_rows(facet_tbls)
}

build_step_table <- function(config, prep, params) {
  if (config$model == "RSM") {
    return(
      tibble(
        Step = paste0("Step_", seq_len(config$n_cat - 1)),
        Estimate = params$steps
      )
    )
  }

  expand_grid(
    StepFacet = prep$levels[[config$step_facet]],
    Step = paste0("Step_", seq_len(config$n_cat - 1))
  ) |>
    mutate(Estimate = as.vector(t(params$steps_mat)))
}

build_estimation_summary <- function(model, method, prep, config, sizes, opt) {
  k_params <- sum(unlist(sizes))
  loglik <- -opt$value
  n_obs <- if (!is.null(config$weight_col) && "Weight" %in% names(prep$data)) {
    sum(prep$data$Weight, na.rm = TRUE)
  } else {
    nrow(prep$data)
  }
  aic <- 2 * k_params - 2 * loglik
  bic <- log(n_obs) * k_params - 2 * loglik

  tibble(
    Model = model,
    Method = method,
    N = n_obs,
    Persons = config$n_person,
    Facets = length(config$facet_names),
    Categories = config$n_cat,
    LogLik = loglik,
    AIC = aic,
    BIC = bic,
    Converged = opt$convergence == 0,
    Iterations = opt$counts[["function"]]
  )
}

# ---- estimation wrapper ----
mfrm_estimate <- function(data, person_col, facet_cols, score_col,
                          rating_min = NULL, rating_max = NULL,
                          weight_col = NULL, keep_original = FALSE,
                          model = c("RSM", "PCM"), method = c("JMLE", "MML"),
                          step_facet = NULL,
                          anchor_df = NULL,
                          group_anchor_df = NULL,
                          noncenter_facet = "Person",
                          dummy_facets = character(0),
                          positive_facets = character(0),
                          quad_points = 15, maxit = 400, reltol = 1e-6) {
  # Stage 1: Normalize model options and input data.
  model <- match.arg(model)
  method <- match.arg(method)

  prep <- prepare_mfrm_data(
    data,
    person_col = person_col,
    facet_cols = facet_cols,
    score_col = score_col,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight_col,
    keep_original = keep_original
  )

  # Stage 2: Resolve facet-level modeling choices.
  step_facet <- resolve_pcm_step_facet(model, step_facet, prep$facet_names)
  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)
  sign_info <- build_facet_signs(prep$facet_names, positive_facets)

  # Stage 3: Build reusable structures for optimization.
  idx <- build_indices(prep, step_facet = step_facet)
  cfg <- build_estimation_config(
    prep = prep,
    model = model,
    method = method,
    step_facet = step_facet,
    weight_col = weight_col,
    facet_signs = sign_info$signs,
    positive_facets = sign_info$positive_facets,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df
  )
  config <- cfg$config
  config$estimation_control <- list(
    maxit = as.integer(maxit),
    reltol = as.numeric(reltol),
    quad_points = as.integer(quad_points)
  )
  sizes <- cfg$sizes
  start <- build_initial_param_vector(config, sizes)

  # Stage 4: Optimize model parameters.
  opt <- run_mfrm_optimization(
    start = start,
    method = method,
    idx = idx,
    config = config,
    sizes = sizes,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol
  )

  # Stage 5: Build human-readable output tables.
  params <- expand_params(opt$par, sizes, config)
  person_tbl <- build_person_table(method, idx, config, params, prep, quad_points)
  facet_tbl <- build_other_facet_table(config, prep, params)
  step_tbl <- build_step_table(config, prep, params)
  summary_tbl <- build_estimation_summary(model, method, prep, config, sizes, opt)

  list(
    summary = summary_tbl,
    facets = list(
      person = person_tbl,
      others = facet_tbl
    ),
    steps = step_tbl,
    config = config,
    prep = prep,
    opt = opt
  )
}

expected_score_table <- function(res) {
  prep <- res$prep
  idx <- build_indices(prep, step_facet = res$config$step_facet)
  config <- res$config
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                                criterion_splits = idx$criterion_splits)
  }
  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  tibble(
    Observed = prep$data$Score,
    Expected = prep$rating_min + expected_k
  )
}

compute_obs_table <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  person_levels <- prep$levels$Person
  person_measure <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate[match(person_levels, res$facets$person$Person)]
  }
  person_measure_by_row <- person_measure[idx$person]
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                                criterion_splits = idx$criterion_splits)
  }

  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
  var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
  resid_k <- idx$score_k - expected_k
  std_sq <- resid_k^2 / var_k

  prep$data |>
    mutate(
      PersonMeasure = person_measure_by_row,
      Observed = prep$rating_min + idx$score_k,
      Expected = prep$rating_min + expected_k,
      Var = var_k,
      Residual = Observed - Expected,
      StdResidual = Residual / sqrt(Var),
      StdSq = std_sq
    )
}

extract_bias_facet_spec <- function(bias_results, data_cols = NULL) {
  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(NULL)
  }
  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)

  facets <- character(0)
  if (!is.null(bias_results$interaction_facets)) {
    facets <- c(facets, as.character(bias_results$interaction_facets))
  }
  if (!is.null(bias_results$facet_a) && length(bias_results$facet_a) > 0) {
    facets <- c(facets, as.character(bias_results$facet_a[1]))
  }
  if (!is.null(bias_results$facet_b) && length(bias_results$facet_b) > 0) {
    facets <- c(facets, as.character(bias_results$facet_b[1]))
  }
  facets <- facets[!is.na(facets) & nzchar(facets)]
  facets <- unique(facets)

  level_cols <- grep("^Facet[0-9]+_Level$", names(tbl), value = TRUE)
  if (length(level_cols) > 0) {
    ord <- suppressWarnings(as.integer(sub("^Facet([0-9]+)_Level$", "\\1", level_cols)))
    ok <- is.finite(ord)
    level_cols <- level_cols[ok]
    ord <- ord[ok]
    if (length(level_cols) == 0) return(NULL)
    o <- order(ord)
    ord <- ord[o]
    level_cols <- level_cols[o]

    facet_name_cols <- paste0("Facet", ord)
    if (all(facet_name_cols %in% names(tbl)) && nrow(tbl) > 0) {
      facets_from_tbl <- vapply(facet_name_cols, function(col) as.character(tbl[[col]][1]), character(1))
      if (all(!is.na(facets_from_tbl) & nzchar(facets_from_tbl))) {
        facets <- facets_from_tbl
      }
    }
    if (length(facets) < length(level_cols)) {
      return(NULL)
    }
    facets <- facets[seq_along(level_cols)]
    index_cols <- paste0("Facet", ord, "_Index")
    measure_cols <- paste0("Facet", ord, "_Measure")
    se_cols <- paste0("Facet", ord, "_SE")
  } else if (all(c("FacetA_Level", "FacetB_Level") %in% names(tbl))) {
    if (length(facets) < 2 && all(c("FacetA", "FacetB") %in% names(tbl)) && nrow(tbl) > 0) {
      facets <- c(as.character(tbl$FacetA[1]), as.character(tbl$FacetB[1]))
    }
    facets <- facets[!is.na(facets) & nzchar(facets)]
    facets <- unique(facets)
    if (length(facets) < 2) return(NULL)
    facets <- facets[1:2]
    level_cols <- c("FacetA_Level", "FacetB_Level")
    index_cols <- c("FacetA_Index", "FacetB_Index")
    measure_cols <- c("FacetA_Measure", "FacetB_Measure")
    se_cols <- c("FacetA_SE", "FacetB_SE")
  } else {
    return(NULL)
  }

  if (!is.null(data_cols) && !all(facets %in% data_cols)) {
    return(NULL)
  }

  list(
    facets = facets,
    level_cols = level_cols,
    index_cols = index_cols,
    measure_cols = measure_cols,
    se_cols = se_cols,
    interaction_order = length(facets),
    interaction_mode = ifelse(length(facets) > 2, "higher_order", "pairwise")
  )
}

compute_bias_adjustment_vector <- function(res, bias_results = NULL) {
  n <- nrow(res$prep$data)
  if (n == 0) return(numeric(0))
  adj <- rep(0, n)

  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(adj)
  }

  data <- res$prep$data
  spec <- extract_bias_facet_spec(bias_results, data_cols = names(data))
  if (is.null(spec) || length(spec$facets) < 2) {
    return(adj)
  }

  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  req_cols <- c(spec$level_cols, "Bias Size")
  if (!all(req_cols %in% names(tbl))) {
    return(adj)
  }

  level_cols <- spec$level_cols
  tbl <- tbl |>
    mutate(
      across(all_of(level_cols), as.character),
      `Bias Size` = suppressWarnings(as.numeric(.data$`Bias Size`))
    ) |>
    filter(is.finite(.data$`Bias Size`))
  if (nrow(tbl) == 0) return(adj)

  # Duplicate rows can appear after manual edits; average as a safe default.
  lookup <- tbl |>
    group_by(across(all_of(level_cols))) |>
    summarize(Bias = mean(.data$`Bias Size`, na.rm = TRUE), .groups = "drop")
  lookup_key_parts <- lapply(level_cols, function(col) as.character(lookup[[col]]))
  lookup$Key <- do.call(paste, c(lookup_key_parts, sep = "||"))

  data_key_parts <- lapply(spec$facets, function(f) as.character(data[[f]]))
  names(data_key_parts) <- spec$facets
  keys <- do.call(paste, c(data_key_parts, sep = "||"))

  idx_hit <- match(keys, lookup$Key)
  hit <- !is.na(idx_hit)
  if (any(hit)) {
    adj[hit] <- lookup$Bias[idx_hit[hit]]
  }
  adj
}

compute_prob_matrix_with_bias <- function(res, bias_results = NULL) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") params$theta else res$facets$person$Estimate
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  bias_adj <- compute_bias_adjustment_vector(res, bias_results = bias_results)
  if (length(bias_adj) == length(eta)) {
    eta <- eta + bias_adj
  }

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                                criterion_splits = idx$criterion_splits)
  }
}

compute_obs_table_with_bias <- function(res, bias_results = NULL) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  person_levels <- prep$levels$Person
  person_measure <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate[match(person_levels, res$facets$person$Person)]
  }
  person_measure_by_row <- person_measure[idx$person]
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  bias_adj <- compute_bias_adjustment_vector(res, bias_results = bias_results)
  if (length(bias_adj) == length(eta)) {
    eta <- eta + bias_adj
  } else {
    bias_adj <- rep(0, length(eta))
  }

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                                criterion_splits = idx$criterion_splits)
  }

  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
  var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
  resid_k <- idx$score_k - expected_k
  std_sq <- resid_k^2 / var_k

  prep$data |>
    mutate(
      PersonMeasure = person_measure_by_row,
      Observed = prep$rating_min + idx$score_k,
      Expected = prep$rating_min + expected_k,
      Var = var_k,
      Residual = Observed - Expected,
      StdResidual = Residual / sqrt(Var),
      StdSq = std_sq,
      BiasAdjustment = bias_adj
    )
}

calc_unexpected_response_table <- function(obs_df,
                                           probs,
                                           facet_names,
                                           rating_min,
                                           abs_z_min = 2,
                                           prob_max = 0.30,
                                           top_n = 100,
                                           rule = c("either", "both")) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  if (is.null(obs_df) || nrow(obs_df) == 0 || is.null(probs) || nrow(probs) != nrow(obs_df)) {
    return(tibble())
  }
  if (!"score_k" %in% names(obs_df) || !"StdResidual" %in% names(obs_df)) {
    return(tibble())
  }

  score_k <- suppressWarnings(as.integer(obs_df$score_k))
  n_cat <- ncol(probs)
  valid <- is.finite(score_k) & score_k >= 0 & score_k < n_cat
  obs_prob <- rep(NA_real_, nrow(obs_df))
  row_idx <- which(valid)
  if (length(row_idx) > 0) {
    obs_prob[row_idx] <- probs[cbind(row_idx, score_k[row_idx] + 1L)]
  }

  max_idx <- max.col(probs, ties.method = "first")
  most_likely_k <- max_idx - 1L
  most_likely <- rating_min + most_likely_k
  most_likely_prob <- probs[cbind(seq_len(nrow(probs)), max_idx)]

  abs_std <- abs(obs_df$StdResidual)
  surprise <- -log10(pmax(obs_prob, .Machine$double.xmin))
  cat_gap <- abs(obs_df$Observed - most_likely)
  low_prob <- is.finite(obs_prob) & obs_prob <= prob_max
  high_resid <- is.finite(abs_std) & abs_std >= abs_z_min
  flagged <- if (rule == "both") {
    low_prob & high_resid
  } else {
    low_prob | high_resid
  }
  if (!any(flagged, na.rm = TRUE)) return(tibble())

  out <- obs_df |>
    mutate(
      Row = dplyr::row_number(),
      ObsProb = obs_prob,
      MostLikely = most_likely,
      MostLikelyProb = most_likely_prob,
      Surprise = surprise,
      AbsStdResidual = abs_std,
      CategoryGap = cat_gap,
      FlagLowProbability = low_prob,
      FlagLargeResidual = high_resid,
      Unexpected = flagged,
      Direction = dplyr::case_when(
        .data$Residual > 0 ~ "Higher than expected",
        .data$Residual < 0 ~ "Lower than expected",
        TRUE ~ "As expected"
      ),
      Severity = .data$AbsStdResidual + .data$Surprise + 0.5 * .data$CategoryGap
    ) |>
    filter(.data$Unexpected)

  id_cols <- c("Person", facet_names)
  if ("Weight" %in% names(obs_df)) id_cols <- c(id_cols, "Weight")
  id_cols <- unique(id_cols)
  keep_cols <- c(
    "Row",
    id_cols,
    "Score",
    "Observed",
    "Expected",
    "Residual",
    "StdResidual",
    "ObsProb",
    "MostLikely",
    "MostLikelyProb",
    "CategoryGap",
    "Surprise",
    "Direction",
    "FlagLowProbability",
    "FlagLargeResidual",
    "Severity"
  )
  keep_cols <- keep_cols[keep_cols %in% names(out)]

  top_n <- max(1L, as.integer(top_n))
  out |>
    arrange(desc(.data$Severity), desc(abs(.data$StdResidual)), .data$ObsProb) |>
    select(dplyr::all_of(keep_cols)) |>
    slice_head(n = top_n)
}

summarize_unexpected_response_table <- function(unexpected_tbl,
                                                total_observations,
                                                abs_z_min = 2,
                                                prob_max = 0.30,
                                                rule = "either") {
  if (is.null(unexpected_tbl) || nrow(unexpected_tbl) == 0) {
    return(tibble(
      TotalObservations = total_observations,
      UnexpectedN = 0L,
      UnexpectedPercent = 0,
      LowProbabilityN = 0L,
      LargeResidualN = 0L,
      Rule = rule,
      AbsZThreshold = abs_z_min,
      ProbThreshold = prob_max
    ))
  }
  low_n <- if ("FlagLowProbability" %in% names(unexpected_tbl)) {
    sum(unexpected_tbl$FlagLowProbability, na.rm = TRUE)
  } else {
    NA_integer_
  }
  resid_n <- if ("FlagLargeResidual" %in% names(unexpected_tbl)) {
    sum(unexpected_tbl$FlagLargeResidual, na.rm = TRUE)
  } else {
    NA_integer_
  }
  tibble(
    TotalObservations = total_observations,
    UnexpectedN = nrow(unexpected_tbl),
    UnexpectedPercent = ifelse(total_observations > 0, 100 * nrow(unexpected_tbl) / total_observations, NA_real_),
    LowProbabilityN = low_n,
    LargeResidualN = resid_n,
    Rule = rule,
    AbsZThreshold = abs_z_min,
    ProbThreshold = prob_max
  )
}

compute_prob_matrix <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                                criterion_splits = idx$criterion_splits)
  }
  probs
}

calc_displacement_table <- function(obs_df,
                                    res,
                                    measures = NULL,
                                    abs_displacement_warn = 0.5,
                                    abs_t_warn = 2) {
  if (is.null(obs_df) || nrow(obs_df) == 0 || is.null(res$config)) {
    return(tibble())
  }

  facet_cols <- c("Person", res$config$facet_names)
  obs_df <- obs_df |>
    mutate(.Weight = get_weights(obs_df))

  displacement <- purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        WeightedN = sum(.data$.Weight, na.rm = TRUE),
        ResidualSum = sum(.data$Residual * .data$.Weight, na.rm = TRUE),
        Information = sum(.data$Var * .data$.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Facet = facet,
        Level = as.character(.data[[facet]])
      ) |>
      select(
        "Facet",
        "Level",
        "WeightedN",
        "ResidualSum",
        "Information"
      )
  })

  if (nrow(displacement) == 0) return(tibble())

  displacement <- displacement |>
    mutate(
      Displacement = ifelse(.data$Information > 0, .data$ResidualSum / .data$Information, NA_real_),
      DisplacementSE = ifelse(.data$Information > 0, 1 / sqrt(.data$Information), NA_real_),
      DisplacementT = ifelse(
        is.finite(.data$DisplacementSE) & .data$DisplacementSE > 0,
        .data$Displacement / .data$DisplacementSE,
        NA_real_
      )
    )

  if (!is.null(measures) && nrow(measures) > 0) {
    displacement <- displacement |>
      left_join(
        measures |>
          select("Facet", "Level", "Estimate", "SE", "N"),
        by = c("Facet", "Level")
      )
  } else {
    displacement <- displacement |>
      mutate(Estimate = NA_real_, SE = NA_real_, N = .data$WeightedN)
  }

  anchor_tbl <- extract_anchor_tables(res$config)$anchors
  if (is.null(anchor_tbl) || nrow(anchor_tbl) == 0) {
    anchor_tbl <- tibble(Facet = character(0), Level = character(0), AnchorValue = numeric(0))
  } else {
    anchor_tbl <- anchor_tbl |>
      transmute(
        Facet = as.character(.data$Facet),
        Level = as.character(.data$Level),
        AnchorValue = as.numeric(.data$Anchor)
      )
  }

  status_tbl <- purrr::map_dfr(unique(displacement$Facet), function(facet) {
    lv <- displacement |>
      filter(.data$Facet == facet) |>
      pull(.data$Level)
    tibble(
      Facet = facet,
      Level = lv,
      AnchorStatus = facet_anchor_status(facet, lv, res$config)
    )
  })

  displacement |>
    left_join(anchor_tbl, by = c("Facet", "Level")) |>
    left_join(status_tbl, by = c("Facet", "Level")) |>
    mutate(
      AnchorType = case_when(
        .data$AnchorStatus == "A" ~ "Anchor",
        .data$AnchorStatus == "G" ~ "Group",
        TRUE ~ "Free"
      ),
      ReleasedEstimate = ifelse(
        is.finite(.data$Estimate) & is.finite(.data$Displacement),
        .data$Estimate + .data$Displacement,
        NA_real_
      ),
      AnchorGap = ifelse(
        is.finite(.data$AnchorValue) & is.finite(.data$ReleasedEstimate),
        .data$ReleasedEstimate - .data$AnchorValue,
        NA_real_
      ),
      FlagDisplacement = is.finite(.data$Displacement) & abs(.data$Displacement) >= abs_displacement_warn,
      FlagT = is.finite(.data$DisplacementT) & abs(.data$DisplacementT) >= abs_t_warn,
      Flag = .data$FlagDisplacement | .data$FlagT
    ) |>
    arrange(desc(abs(.data$Displacement)), desc(abs(.data$DisplacementT)))
}

summarize_displacement_table <- function(displacement_tbl,
                                         abs_displacement_warn = 0.5,
                                         abs_t_warn = 2) {
  if (is.null(displacement_tbl) || nrow(displacement_tbl) == 0) {
    return(tibble(
      Levels = 0L,
      AnchoredLevels = 0L,
      FlaggedLevels = 0L,
      FlaggedAnchoredLevels = 0L,
      MaxAbsDisplacement = NA_real_,
      MaxAbsDisplacementT = NA_real_,
      AbsDisplacementThreshold = abs_displacement_warn,
      AbsTThreshold = abs_t_warn
    ))
  }
  is_anchored <- displacement_tbl$AnchorType %in% c("Anchor", "Group")
  flagged <- if ("Flag" %in% names(displacement_tbl)) {
    displacement_tbl$Flag
  } else {
    rep(FALSE, nrow(displacement_tbl))
  }
  tibble(
    Levels = nrow(displacement_tbl),
    AnchoredLevels = sum(is_anchored, na.rm = TRUE),
    FlaggedLevels = sum(flagged, na.rm = TRUE),
    FlaggedAnchoredLevels = sum(flagged & is_anchored, na.rm = TRUE),
    MaxAbsDisplacement = max(abs(displacement_tbl$Displacement), na.rm = TRUE),
    MaxAbsDisplacementT = max(abs(displacement_tbl$DisplacementT), na.rm = TRUE),
    AbsDisplacementThreshold = abs_displacement_warn,
    AbsTThreshold = abs_t_warn
  )
}

compute_scorefile <- function(res) {
  obs <- compute_obs_table(res)
  if (nrow(obs) == 0) return(tibble())
  probs <- compute_prob_matrix(res)
  if (is.null(probs) || nrow(probs) != nrow(obs)) return(obs)
  cat_vals <- seq(res$prep$rating_min, res$prep$rating_max)
  prob_df <- as_tibble(probs, .name_repair = "minimal")
  names(prob_df) <- paste0("P_", cat_vals)
  max_idx <- apply(probs, 1, which.max)
  max_prob <- probs[cbind(seq_len(nrow(probs)), max_idx)]
  most_likely <- cat_vals[max_idx]
  base_cols <- c("Person", res$config$facet_names)
  if ("Weight" %in% names(obs)) {
    base_cols <- c(base_cols, "Weight")
  }
  base_cols <- c(base_cols, "Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "PersonMeasure")
  bind_cols(
    obs |>
      select(all_of(base_cols)),
    prob_df
  ) |>
    mutate(
      MostLikely = most_likely,
      MaxProb = max_prob
    )
}

compute_residual_file <- function(res) {
  obs <- compute_obs_table(res)
  if (nrow(obs) == 0) return(tibble())
  base_cols <- c("Person", res$config$facet_names)
  if ("Weight" %in% names(obs)) {
    base_cols <- c(base_cols, "Weight")
  }
  base_cols <- c(base_cols, "Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "StdSq", "PersonMeasure")
  obs |>
    select(all_of(base_cols))
}

# Overall model fit: weighted infit and outfit mean-square statistics.
# Infit (information-weighted): sum(StdSq * Var * w) / sum(Var * w)
#   Sensitive to unexpected responses near the person's ability level.
# Outfit (unweighted): sum(StdSq * w) / sum(w)
#   Sensitive to outlying unexpected responses far from the ability level.
# Both transformed to ZSTD via Wilson-Hilferty for significance testing.
calc_overall_fit <- function(obs_df, whexact = FALSE) {
  w <- get_weights(obs_df)
  infit <- sum(obs_df$StdSq * obs_df$Var * w, na.rm = TRUE) / sum(obs_df$Var * w, na.rm = TRUE)
  outfit <- sum(obs_df$StdSq * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  df_infit <- sum(obs_df$Var * w, na.rm = TRUE)
  df_outfit <- sum(w, na.rm = TRUE)
  tibble(
    Infit = infit,
    Outfit = outfit,
    InfitZSTD = zstd_from_mnsq(infit, df_infit, whexact = whexact),
    OutfitZSTD = zstd_from_mnsq(outfit, df_outfit, whexact = whexact),
    DF_Infit = df_infit,
    DF_Outfit = df_outfit
  )
}

calc_facet_fit <- function(obs_df, facet_cols, whexact = FALSE) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    df <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        Infit = sum(StdSq * Var * .Weight, na.rm = TRUE) / sum(Var * .Weight, na.rm = TRUE),
        Outfit = sum(StdSq * .Weight, na.rm = TRUE) / sum(.Weight, na.rm = TRUE),
        DF_Infit = sum(Var * .Weight, na.rm = TRUE),
        DF_Outfit = sum(.Weight, na.rm = TRUE),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      )
    df |>
      mutate(
        InfitZSTD = zstd_from_mnsq(Infit, DF_Infit, whexact = whexact),
        OutfitZSTD = zstd_from_mnsq(Outfit, DF_Outfit, whexact = whexact)
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      select(Facet, Level, N, Infit, Outfit, InfitZSTD, OutfitZSTD, DF_Infit, DF_Outfit)
  })
}

# Approximate facet-level SE from summed observation information.
calc_facet_se <- function(obs_df, facet_cols) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        Info = sum(Var * .Weight, na.rm = TRUE),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Facet = facet,
        Level = .data[[facet]],
        SE = ifelse(Info > 0, 1 / sqrt(Info), NA_real_)
      ) |>
      select(Facet, Level, N, SE)
  })
}

calc_bias_facet <- function(obs_df, facet_cols) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        ObservedAverage = weighted_mean(Observed, .Weight),
        ExpectedAverage = weighted_mean(Expected, .Weight),
        MeanResidual = weighted_mean(Residual, .Weight),
        MeanStdResidual = weighted_mean(StdResidual, .Weight),
        MeanAbsStdResidual = weighted_mean(abs(StdResidual), .Weight),
        ChiSq = sum((StdResidual^2) * .Weight, na.rm = TRUE),
        SE_Residual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, sqrt(sum(Var * .Weight, na.rm = TRUE)) / n_w, NA_real_)
        },
        SE_StdResidual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, 1 / sqrt(n_w), NA_real_)
        },
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      mutate(
        Bias = MeanResidual,
        DF = ifelse(N > 1, N - 1, NA_real_),
        t_Residual = ifelse(is.finite(SE_Residual) & SE_Residual > 0, MeanResidual / SE_Residual, NA_real_),
        t_StdResidual = ifelse(is.finite(SE_StdResidual) & SE_StdResidual > 0, MeanStdResidual / SE_StdResidual, NA_real_),
        p_Residual = ifelse(is.finite(DF) & is.finite(t_Residual), 2 * stats::pt(-abs(t_Residual), df = DF), NA_real_),
        p_StdResidual = ifelse(is.finite(DF) & is.finite(t_StdResidual), 2 * stats::pt(-abs(t_StdResidual), df = DF), NA_real_),
        ChiDf = DF,
        ChiP = ifelse(is.finite(ChiSq) & is.finite(ChiDf) & ChiDf > 0, 1 - stats::pchisq(ChiSq, df = ChiDf), NA_real_)
      ) |>
      select(Facet, Level, N, ObservedAverage, ExpectedAverage, Bias, MeanResidual, MeanStdResidual,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP, SE_Residual, t_Residual, p_Residual,
             SE_StdResidual, t_StdResidual, p_StdResidual, DF)
  })
}

calc_bias_interactions <- function(obs_df, facet_cols, pairs = NULL, top_n = 20) {
  if (length(facet_cols) < 2) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  if (is.null(pairs)) {
    combos <- combn(facet_cols, 2, simplify = FALSE)
  } else if (length(pairs) == 0) {
    return(tibble())
  } else {
    combos <- pairs
  }
  out <- purrr::map_dfr(combos, function(pair) {
    pair1 <- pair[1]
    pair2 <- pair[2]
    obs_df |>
      group_by(.data[[pair1]], .data[[pair2]]) |>
      summarize(
        ObservedAverage = weighted_mean(Observed, .Weight),
        ExpectedAverage = weighted_mean(Expected, .Weight),
        MeanResidual = weighted_mean(Residual, .Weight),
        MeanStdResidual = weighted_mean(StdResidual, .Weight),
        MeanAbsStdResidual = weighted_mean(abs(StdResidual), .Weight),
        ChiSq = sum((StdResidual^2) * .Weight, na.rm = TRUE),
        SE_Residual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, sqrt(sum(Var * .Weight, na.rm = TRUE)) / n_w, NA_real_)
        },
        SE_StdResidual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, 1 / sqrt(n_w), NA_real_)
        },
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Pair = paste(pair1, "x", pair2),
        Level = paste(.data[[pair1]], .data[[pair2]], sep = " | ")
      ) |>
      mutate(
        Bias = MeanResidual,
        DF = ifelse(N > 1, N - 1, NA_real_),
        t_Residual = ifelse(is.finite(SE_Residual) & SE_Residual > 0, MeanResidual / SE_Residual, NA_real_),
        t_StdResidual = ifelse(is.finite(SE_StdResidual) & SE_StdResidual > 0, MeanStdResidual / SE_StdResidual, NA_real_),
        p_Residual = ifelse(is.finite(DF) & is.finite(t_Residual), 2 * stats::pt(-abs(t_Residual), df = DF), NA_real_),
        p_StdResidual = ifelse(is.finite(DF) & is.finite(t_StdResidual), 2 * stats::pt(-abs(t_StdResidual), df = DF), NA_real_),
        ChiDf = DF,
        ChiP = ifelse(is.finite(ChiSq) & is.finite(ChiDf) & ChiDf > 0, 1 - stats::pchisq(ChiSq, df = ChiDf), NA_real_)
      ) |>
      select(Pair, Level, N, ObservedAverage, ExpectedAverage, Bias, MeanResidual, MeanStdResidual,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP, SE_Residual, t_Residual, p_Residual,
             SE_StdResidual, t_StdResidual, p_StdResidual, DF)
  })
  out |>
    mutate(AbsStd = abs(MeanStdResidual)) |>
    arrange(desc(AbsStd)) |>
    select(-AbsStd) |>
    slice_head(n = top_n)
}

safe_cor <- function(x, y, w = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (is.null(w)) {
    if (!any(ok)) return(NA_real_)
    x <- x[ok]
    y <- y[ok]
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(NA_real_)
    }
    return(suppressWarnings(stats::cor(x, y, use = "complete.obs")))
  }
  ok <- ok & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  w_sum <- sum(w)
  if (w_sum <= 0) return(NA_real_)
  mx <- sum(w * x) / w_sum
  my <- sum(w * y) / w_sum
  vx <- sum(w * (x - mx)^2) / w_sum
  vy <- sum(w * (y - my)^2) / w_sum
  if (vx <= 0 || vy <= 0) return(NA_real_)
  cov <- sum(w * (x - mx) * (y - my)) / w_sum
  cov / sqrt(vx * vy)
}

weighted_mean_safe <- function(x, w) {
  weighted_mean(x, w)
}

infer_default_rater_facet <- function(facet_names) {
  facet_names <- as.character(facet_names)
  if (length(facet_names) == 0) return(NULL)

  lowered <- tolower(facet_names)
  patterns <- c("rater", "judge", "grader", "reader", "scorer", "assessor", "evaluator")
  hits <- vapply(
    lowered,
    function(x) any(vapply(patterns, function(p) grepl(p, x, fixed = TRUE), logical(1))),
    logical(1)
  )
  if (any(hits)) {
    return(facet_names[which(hits)[1]])
  }
  facet_names[1]
}

calc_interrater_agreement <- function(obs_df, facet_cols, rater_facet, res = NULL) {
  if (is.null(obs_df) || nrow(obs_df) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }
  if (is.null(rater_facet) || !rater_facet %in% facet_cols) {
    return(list(summary = tibble(), pairs = tibble()))
  }
  context_cols <- setdiff(facet_cols, rater_facet)
  if (length(context_cols) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  df <- obs_df |>
    mutate(across(all_of(context_cols), as.character)) |>
    tidyr::unite(".context", all_of(context_cols), sep = "|", remove = FALSE) |>
    select(.context, !!rlang::sym(rater_facet), Observed, any_of("Weight"))
  df$.Weight <- get_weights(df)

  df <- df |>
    group_by(.context, !!rlang::sym(rater_facet)) |>
    summarize(Score = weighted_mean(Observed, .Weight), .groups = "drop")

  if (nrow(df) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  wide <- tryCatch(
    tidyr::pivot_wider(
      df,
      id_cols = .context,
      names_from = !!rlang::sym(rater_facet),
      values_from = Score
    ),
    error = function(e) NULL
  )
  if (is.null(wide)) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  rater_cols <- setdiff(names(wide), ".context")
  if (length(rater_cols) < 2) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  prob_map <- list()
  if (!is.null(res)) {
    probs <- compute_prob_matrix(res)
    if (!is.null(probs) && nrow(probs) == nrow(obs_df)) {
      prob_cols <- paste0(".p", seq_len(ncol(probs)))
      prob_df <- obs_df |>
        mutate(across(all_of(context_cols), as.character)) |>
        tidyr::unite(".context", all_of(context_cols), sep = "|", remove = FALSE) |>
        select(.context, !!rlang::sym(rater_facet), any_of("Weight"))
      prob_df[prob_cols] <- probs
      prob_df$.Weight <- get_weights(prob_df)

      prob_avg <- prob_df |>
        group_by(.context, !!rlang::sym(rater_facet)) |>
        summarize(
          across(all_of(prob_cols), ~ weighted_mean(.x, .Weight)),
          .groups = "drop"
        )

      if (nrow(prob_avg) > 0) {
        for (i in seq_len(nrow(prob_avg))) {
          ctx <- prob_avg$.context[[i]]
          rater_val <- prob_avg[[rater_facet]][[i]]
          key <- paste(ctx, rater_val, sep = "||")
          prob_map[[key]] <- as.numeric(prob_avg[i, prob_cols, drop = TRUE])
        }
      }
    }
  }

  pairs <- combn(rater_cols, 2, simplify = FALSE)
  pair_tbl <- purrr::map_dfr(pairs, function(pair) {
    sub <- wide |>
      select(.context, !!rlang::sym(pair[1]), !!rlang::sym(pair[2])) |>
      filter(is.finite(.data[[pair[1]]]), is.finite(.data[[pair[2]]]))
    n_ok <- nrow(sub)
    if (n_ok == 0) {
      return(tibble(
        Rater1 = pair[1],
        Rater2 = pair[2],
        N = 0,
        Exact = NA_real_,
        ExpectedExact = NA_real_,
        Adjacent = NA_real_,
        MeanDiff = NA_real_,
        MAD = NA_real_,
        Corr = NA_real_
      ))
    }
    v1 <- sub[[pair[1]]]
    v2 <- sub[[pair[2]]]
    diff <- v1 - v2
    exact_count <- sum(diff == 0, na.rm = TRUE)

    exp_vals <- numeric(0)
    if (length(prob_map) > 0) {
      for (ctx in sub$.context) {
        key1 <- paste(ctx, pair[1], sep = "||")
        key2 <- paste(ctx, pair[2], sep = "||")
        p1 <- prob_map[[key1]]
        p2 <- prob_map[[key2]]
        if (is.null(p1) || is.null(p2)) next
        if (any(!is.finite(p1)) || any(!is.finite(p2))) next
        exp_vals <- c(exp_vals, sum(p1 * p2))
      }
    }
    exp_mean <- if (length(exp_vals) > 0) mean(exp_vals) else NA_real_

    tibble(
      Rater1 = pair[1],
      Rater2 = pair[2],
      N = n_ok,
      OpportunityCount = n_ok,
      ExactCount = exact_count,
      ExpectedExactCount = if (length(exp_vals) > 0) sum(exp_vals) else NA_real_,
      AdjacentCount = sum(abs(diff) <= 1, na.rm = TRUE),
      Exact = exact_count / n_ok,
      ExpectedExact = exp_mean,
      Adjacent = mean(abs(diff) <= 1, na.rm = TRUE),
      MeanDiff = mean(diff, na.rm = TRUE),
      MAD = mean(abs(diff), na.rm = TRUE),
      Corr = safe_cor(v1, v2)
    )
  })

  contexts_with_pairs <- sum(rowSums(!is.na(wide[rater_cols])) >= 2)
  total_pairs <- sum(pair_tbl$N, na.rm = TRUE)
  total_exact <- sum(pair_tbl$Exact * pair_tbl$N, na.rm = TRUE)
  expected_available <- any(is.finite(pair_tbl$ExpectedExact))
  total_expected <- if (expected_available) {
    sum(pair_tbl$ExpectedExact * pair_tbl$N, na.rm = TRUE)
  } else {
    NA_real_
  }
  summary_tbl <- tibble(
    RaterFacet = rater_facet,
    Raters = length(rater_cols),
    Pairs = nrow(pair_tbl),
    Contexts = contexts_with_pairs,
    TotalPairs = total_pairs,
    OpportunityCount = total_pairs,
    ExactAgreements = total_exact,
    ExpectedAgreements = ifelse(expected_available, total_expected, NA_real_),
    ExactAgreement = ifelse(total_pairs > 0, total_exact / total_pairs, NA_real_),
    ExpectedExactAgreement = ifelse(expected_available && total_pairs > 0, total_expected / total_pairs, NA_real_),
    AgreementMinusExpected = ifelse(
      expected_available && total_pairs > 0,
      (total_exact - total_expected) / total_pairs,
      NA_real_
    ),
    AdjacentAgreements = sum(pair_tbl$AdjacentCount, na.rm = TRUE),
    AdjacentAgreement = weighted_mean_safe(pair_tbl$Adjacent, pair_tbl$N),
    MeanAbsDiff = weighted_mean_safe(pair_tbl$MAD, pair_tbl$N),
    MeanCorr = weighted_mean_safe(pair_tbl$Corr, pair_tbl$N)
  )

  list(summary = summary_tbl, pairs = pair_tbl)
}

augment_interrater_with_precision <- function(agreement, reliability_tbl, rater_facet = NULL) {
  if (is.null(agreement) || !is.list(agreement)) return(agreement)
  summary_tbl <- as.data.frame(agreement$summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_tbl) == 0) return(agreement)

  if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
    if ("RaterFacet" %in% names(summary_tbl)) {
      rater_facet <- as.character(summary_tbl$RaterFacet[1])
    }
  } else {
    rater_facet <- as.character(rater_facet[1])
  }

  rel_tbl <- as.data.frame(reliability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (!"Facet" %in% names(rel_tbl) || !nzchar(rater_facet)) {
    agreement$summary <- summary_tbl
    return(agreement)
  }

  rel_row <- rel_tbl[as.character(rel_tbl$Facet) == rater_facet, , drop = FALSE]
  if (nrow(rel_row) == 0) {
    agreement$summary <- summary_tbl
    return(agreement)
  }

  summary_tbl$RaterSeparation <- suppressWarnings(as.numeric(rel_row$Separation[1]))
  summary_tbl$RaterStrata <- suppressWarnings(as.numeric(rel_row$Strata[1]))
  summary_tbl$RaterReliability <- suppressWarnings(as.numeric(rel_row$Reliability[1]))
  summary_tbl$RaterRealSeparation <- suppressWarnings(as.numeric(rel_row$RealSeparation[1]))
  summary_tbl$RaterRealReliability <- suppressWarnings(as.numeric(rel_row$RealReliability[1]))
  agreement$summary <- summary_tbl
  agreement
}

calc_ptmea <- function(obs_df, facet_cols) {
  facet_cols <- setdiff(facet_cols, "Person")
  if (length(facet_cols) == 0) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        PTMEA = safe_cor(Observed, PersonMeasure, w = .Weight),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      select(Facet, Level, PTMEA, N)
  })
}

expected_score_from_eta <- function(eta, step_cum, rating_min) {
  if (!is.finite(eta) || length(step_cum) == 0) return(NA_real_)
  probs <- category_prob_rsm(eta, step_cum)
  k_vals <- 0:(length(step_cum) - 1)
  rating_min + sum(probs * k_vals)
}

estimate_eta_from_target <- function(target, step_cum, rating_min, rating_max) {
  if (!is.finite(target) || length(step_cum) == 0) return(NA_real_)
  if (target <= rating_min) return(-Inf)
  if (target >= rating_max) return(Inf)
  f <- function(eta) expected_score_from_eta(eta, step_cum, rating_min) - target
  lower <- -10
  upper <- 10
  f_low <- f(lower)
  f_up <- f(upper)
  if (!is.finite(f_low) || !is.finite(f_up)) return(NA_real_)
  if (f_low * f_up > 0) {
    lower <- -20
    upper <- 20
    f_low <- f(lower)
    f_up <- f(upper)
    if (!is.finite(f_low) || !is.finite(f_up) || f_low * f_up > 0) return(NA_real_)
  }
  uniroot(f, lower = lower, upper = upper)$root
}

facet_anchor_status <- function(facet, levels, config) {
  spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
  if (is.null(spec)) return(rep("", length(levels)))
  anchors <- spec$anchors
  groups <- spec$groups
  status <- rep("", length(levels))
  if (!is.null(anchors)) {
    anchor_vals <- anchors[match(levels, names(anchors))]
    status[is.finite(anchor_vals)] <- "A"
  }
  if (!is.null(groups)) {
    group_vals <- groups[match(levels, names(groups))]
    status[status == "" & !is.na(group_vals) & group_vals != ""] <- "G"
  }
  status
}

calc_facets_report_tbls <- function(res,
                                    diagnostics,
                                    totalscore = TRUE,
                                    umean = 0,
                                    uscale = 1,
                                    udecimals = 2,
                                    omit_unobserved = FALSE,
                                    xtreme = 0) {
  # Stage 1: Validate required inputs and extract core model objects.
  if (is.null(res) || is.null(diagnostics)) return(list())
  obs_df <- diagnostics$obs
  measures <- diagnostics$measures
  if (nrow(obs_df) == 0 || nrow(measures) == 0) return(list())

  prep <- res$prep
  config <- res$config
  rating_min <- prep$rating_min
  rating_max <- prep$rating_max
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  theta_mean <- if (length(theta_hat) > 0) mean(theta_hat, na.rm = TRUE) else 0
  facet_means <- purrr::map_dbl(config$facet_names, function(f) {
    mean(params$facets[[f]], na.rm = TRUE)
  })
  names(facet_means) <- config$facet_names
  facet_signs <- config$facet_signs
  if (is.null(facet_signs) || length(facet_signs) == 0) {
    facet_signs <- setNames(rep(-1, length(config$facet_names)), config$facet_names)
  }

  # Stage 2: Build step structures used for fair-average calculations.
  if (config$model == "RSM") {
    step_cum_common <- c(0, cumsum(params$steps))
    step_cum_mean <- step_cum_common
  } else {
    step_mat <- params$steps_mat
    if (is.null(step_mat) || length(step_mat) == 0) {
      step_cum_common <- numeric(0)
      step_cum_mean <- numeric(0)
    } else {
      step_mean <- colMeans(step_mat, na.rm = TRUE)
      step_cum_common <- t(apply(step_mat, 1, function(x) c(0, cumsum(x))))
      step_cum_mean <- c(0, cumsum(step_mean))
    }
  }

  facet_names <- c("Person", config$facet_names)
  facet_levels_all <- lapply(facet_names, function(facet) {
    if (facet == "Person") {
      prep$levels$Person
    } else {
      prep$levels[[facet]]
    }
  })
  names(facet_levels_all) <- facet_names

  # Stage 3: Detect extreme-only levels and cache row-level flags.
  extreme_levels <- purrr::map(facet_names, function(facet) {
    if (!facet %in% names(obs_df)) return(character(0))
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(
        (MinScore == prep$rating_min & MaxScore == prep$rating_min) |
          (MinScore == prep$rating_max & MaxScore == prep$rating_max)
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      pull(Level)
  })
  names(extreme_levels) <- facet_names

  extreme_flags <- list()
  for (facet in facet_names) {
    if (facet %in% names(obs_df)) {
      extreme_flags[[facet]] <- obs_df[[facet]] %in% extreme_levels[[facet]]
    }
  }
  if (length(extreme_flags) > 0) {
    extreme_flag_df <- as_tibble(extreme_flags)
    extreme_count <- rowSums(extreme_flag_df, na.rm = TRUE)
  } else {
    extreme_count <- rep(0, nrow(obs_df))
  }

  # Stage 4: Build one compatibility table per facet.
  out <- purrr::map(facet_names, function(facet) {
    if (!facet %in% names(obs_df)) return(tibble())
    status_tbl <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        TotalCountAll = n(),
        .groups = "drop"
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      select(Level, MinScore, MaxScore, TotalCountAll)

    if (isTRUE(totalscore)) {
      score_source <- obs_df
    } else {
      # Legacy printer-style rule: remove rows with multiple extreme flags,
      # but keep rows where the focal facet is the only extreme flag.
      flag <- extreme_flags[[facet]]
      if (is.null(flag)) {
        flag <- rep(FALSE, nrow(obs_df))
      }
      active_mask <- (extreme_count == 0) | ((extreme_count == 1) & flag)
      score_source <- obs_df[active_mask, , drop = FALSE]
    }
    score_tbl <- score_source |>
      group_by(.data[[facet]]) |>
      summarize(
        TotalScore = sum(Observed, na.rm = TRUE),
        TotalCount = n(),
        WeightdScore = sum(Observed * Weight, na.rm = TRUE),
        WeightdCount = sum(Weight, na.rm = TRUE),
        ObservedAverage = ifelse(sum(Weight, na.rm = TRUE) > 0,
                                 sum(Observed * Weight, na.rm = TRUE) / sum(Weight, na.rm = TRUE),
                                 NA_real_),
        .groups = "drop"
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      select(Level, TotalScore, TotalCount, WeightdScore, WeightdCount, ObservedAverage)

    level_tbl <- tibble(Level = as.character(facet_levels_all[[facet]]))
    score_tbl <- level_tbl |>
      left_join(score_tbl, by = "Level") |>
      mutate(
        TotalScore = replace_na(TotalScore, 0),
        TotalCount = replace_na(TotalCount, 0),
        WeightdScore = replace_na(WeightdScore, 0),
        WeightdCount = replace_na(WeightdCount, 0),
        ObservedAverage = ifelse(WeightdCount > 0, ObservedAverage, NA_real_)
      )

    meas_tbl <- measures |>
      filter(Facet == facet) |>
      mutate(Level = as.character(Level)) |>
      select(Level, Estimate, SE, ModelSE, RealSE, Infit, Outfit, InfitZSTD, OutfitZSTD, PTMEA)

    tbl <- level_tbl |>
      left_join(score_tbl, by = "Level") |>
      left_join(status_tbl, by = "Level") |>
      left_join(meas_tbl, by = "Level")

    anchor_status <- facet_anchor_status(facet, tbl$Level, config)
    status <- dplyr::case_when(
      tbl$TotalCountAll == 0 ~ "No data",
      tbl$MinScore == tbl$MaxScore & tbl$MinScore == rating_min ~ "Minimum",
      tbl$MinScore == tbl$MaxScore & tbl$MaxScore == rating_max ~ "Maximum",
      tbl$TotalCountAll == 1 ~ "One datum",
      TRUE ~ ""
    )

    sign_vec <- facet_signs[names(facet_means)]
    if (facet == "Person") {
      other_sum <- sum(sign_vec * facet_means, na.rm = TRUE)
      eta_m <- tbl$Estimate + other_sum
      eta_z <- tbl$Estimate
    } else {
      sign <- if (!is.null(facet_signs[[facet]])) facet_signs[[facet]] else -1
      other_sum <- sum(sign_vec[names(facet_means) != facet] *
                         facet_means[names(facet_means) != facet], na.rm = TRUE)
      eta_m <- theta_mean + other_sum + sign * tbl$Estimate
      eta_z <- sign * tbl$Estimate
    }

    if (config$model == "PCM" && !is.null(config$step_facet)) {
      step_levels <- prep$levels[[config$step_facet]]
      if (facet == config$step_facet && length(step_levels) > 0 && length(step_cum_common) > 0) {
        step_cum_list <- purrr::map(tbl$Level, function(lvl) {
          idx <- match(lvl, step_levels)
          if (is.na(idx) || idx < 1 || idx > nrow(step_cum_common)) {
            step_cum_mean
          } else {
            step_cum_common[idx, ]
          }
        })
      } else {
        step_cum_list <- rep(list(step_cum_mean), nrow(tbl))
      }
    } else {
      step_cum_list <- rep(list(step_cum_common), nrow(tbl))
    }

    fair_m <- purrr::map2_dbl(eta_m, step_cum_list, ~ expected_score_from_eta(.x, .y, rating_min))
    fair_z <- purrr::map2_dbl(eta_z, step_cum_list, ~ expected_score_from_eta(.x, .y, rating_min))

    xtreme_target <- ifelse(
      status == "Minimum", rating_min + xtreme,
      ifelse(status == "Maximum", rating_max - xtreme, NA_real_)
    )
    xtreme_eta <- purrr::map2_dbl(xtreme_target, step_cum_list, ~ {
      if (!is.finite(.x) || xtreme <= 0) return(NA_real_)
      estimate_eta_from_target(.x, .y, rating_min, rating_max)
    })

    # Convert any xtreme-adjusted eta back to facet-specific measure units.
    measure_logit <- tbl$Estimate
    if (any(is.finite(xtreme_eta))) {
      if (facet == "Person") {
        measure_logit <- ifelse(
          is.finite(xtreme_eta),
          xtreme_eta - other_sum,
          measure_logit
        )
      } else {
        sign <- if (!is.null(facet_signs[[facet]])) facet_signs[[facet]] else -1
        measure_logit <- ifelse(
          is.finite(xtreme_eta),
          (xtreme_eta - theta_mean - other_sum) / sign,
          measure_logit
        )
      }
    }

    scale_factor <- ifelse(is.finite(uscale), uscale, 1)
    scale_origin <- ifelse(is.finite(umean), umean, 0)

    tbl <- tbl |>
      mutate(
        Anchor = anchor_status,
        Status = status,
        FairM = fair_m,
        FairZ = fair_z,
        Measure = ifelse(is.finite(measure_logit), measure_logit * scale_factor + scale_origin, NA_real_),
        ModelSE = ifelse(
          is.finite(dplyr::coalesce(ModelSE, SE)),
          abs(scale_factor) * dplyr::coalesce(ModelSE, SE),
          NA_real_
        ),
        RealSE = ifelse(
          is.finite(dplyr::coalesce(RealSE, ModelSE, SE)),
          abs(scale_factor) * dplyr::coalesce(
            RealSE,
            ModelSE * sqrt(pmax(dplyr::coalesce(Infit, 1), 1)),
            SE * sqrt(pmax(dplyr::coalesce(Infit, 1), 1))
          ),
          NA_real_
        ),
        Infit = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, Infit),
        Outfit = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, Outfit),
        InfitZSTD = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, InfitZSTD),
        OutfitZSTD = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, OutfitZSTD)
      ) |>
      transmute(
        TotalScore,
        TotalCount,
        WeightdScore,
        WeightdCount,
        ObservedAverage,
        FairM,
        FairZ,
        Measure,
        ModelSE,
        RealSE,
        InfitMnSq = Infit,
        InfitZStd = InfitZSTD,
        OutfitMnSq = Outfit,
        OutfitZStd = OutfitZSTD,
        PtMeaCorr = ifelse(facet == "Person", NA_real_, PTMEA),
        Anchor,
        Status,
        Level,
        TotalCountAll
      ) |>
      arrange(desc(Measure), desc(TotalCount))

    if (isTRUE(omit_unobserved)) {
      tbl <- tbl |> filter(TotalCountAll > 0)
    }

    tbl |>
      select(-TotalCountAll)
  })

  names(out) <- facet_names
  out
}

format_facets_report_gt <- function(tbl,
                                    facet,
                                    decimals = 2,
                                    totalscore = TRUE,
                                    reference = c("both", "mean", "zero"),
                                    label_style = c("both", "native", "legacy")) {
  reference <- match.arg(reference)
  label_style <- match.arg(label_style)
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(data.frame(Message = "No facet report available.", stringsAsFactors = FALSE))
  }
  out <- as.data.frame(tbl, stringsAsFactors = FALSE)
  label_map <- c(
    TotalScore = if (isTRUE(totalscore)) "Total Score" else "Obsvd Score",
    TotalCount = if (isTRUE(totalscore)) "Total Count" else "Obsvd Count",
    WeightdScore = "Weightd Score",
    WeightdCount = "Weightd Count",
    ObservedAverage = "Obsvd Average",
    FairM = "Fair(M) Average",
    FairZ = "Fair(Z) Average",
    Measure = "Measure",
    ModelSE = "Model S.E.",
    RealSE = "Real S.E.",
    InfitMnSq = "Infit MnSq",
    InfitZStd = "Infit ZStd",
    OutfitMnSq = "Outfit MnSq",
    OutfitZStd = "Outfit ZStd",
    PtMeaCorr = "PtMea Corr",
    Anchor = "Anch",
    Status = "Status",
    Level = "Element"
  )
  keep <- intersect(names(out), names(label_map))
  names(out)[match(keep, names(out))] <- unname(label_map[keep])
  if ("Element" %in% names(out)) {
    out <- out[, c(setdiff(names(out), "Element"), "Element"), drop = FALSE]
  }

  count_cols <- intersect(c("Total Score", "Total Count", "Weightd Score", "Weightd Count"), names(out))
  for (col in count_cols) out[[col]] <- round(as.numeric(out[[col]]), digits = 0)
  value_cols <- intersect(
    c("Obsvd Average", "Fair(M) Average", "Fair(Z) Average", "Measure", "Model S.E.", "Real S.E.",
      "Infit MnSq", "Infit ZStd", "Outfit MnSq", "Outfit ZStd", "PtMea Corr"),
    names(out)
  )
  for (col in value_cols) out[[col]] <- round(as.numeric(out[[col]]), digits = decimals)

  # Preserve the legacy display labels while exposing package-native aliases.
  if ("Obsvd Average" %in% names(out) && !"ObservedAverage" %in% names(out)) {
    out$ObservedAverage <- out[["Obsvd Average"]]
  }
  if ("Fair(M) Average" %in% names(out) && !"AdjustedAverage" %in% names(out)) {
    out$AdjustedAverage <- out[["Fair(M) Average"]]
  }
  if ("Fair(Z) Average" %in% names(out) && !"StandardizedAdjustedAverage" %in% names(out)) {
    out$StandardizedAdjustedAverage <- out[["Fair(Z) Average"]]
  }
  if ("Model S.E." %in% names(out) && !"ModelBasedSE" %in% names(out)) {
    out$ModelBasedSE <- out[["Model S.E."]]
  }
  if ("Real S.E." %in% names(out) && !"FitAdjustedSE" %in% names(out)) {
    out$FitAdjustedSE <- out[["Real S.E."]]
  }

  if (identical(reference, "mean")) {
    drop_cols <- intersect(c("Fair(Z) Average", "StandardizedAdjustedAverage"), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  } else if (identical(reference, "zero")) {
    drop_cols <- intersect(c("Fair(M) Average", "AdjustedAverage"), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  }

  if (identical(label_style, "native")) {
    drop_cols <- intersect(c("Obsvd Average", "Fair(M) Average", "Fair(Z) Average", "Model S.E.", "Real S.E."), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  } else if (identical(label_style, "legacy")) {
    drop_cols <- intersect(c("ObservedAverage", "AdjustedAverage", "StandardizedAdjustedAverage", "ModelBasedSE", "FitAdjustedSE"), names(out))
    out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
  }

  attr(out, "facet") <- facet
  attr(out, "totalscore") <- isTRUE(totalscore)
  attr(out, "reference") <- reference
  attr(out, "label_style") <- label_style
  out
}

calc_fair_average_bundle <- function(res,
                                     diagnostics,
                                     facets = NULL,
                                     totalscore = TRUE,
                                     umean = 0,
                                     uscale = 1,
                                     udecimals = 2,
                                     reference = c("both", "mean", "zero"),
                                     label_style = c("both", "native", "legacy"),
                                     omit_unobserved = FALSE,
                                     xtreme = 0) {
  reference <- match.arg(reference)
  label_style <- match.arg(label_style)
  raw_tbls <- calc_facets_report_tbls(
    res = res,
    diagnostics = diagnostics,
    totalscore = totalscore,
    umean = umean,
    uscale = uscale,
    udecimals = udecimals,
    omit_unobserved = omit_unobserved,
    xtreme = xtreme
  )
  if (is.null(raw_tbls) || length(raw_tbls) == 0) {
    return(list(raw_by_facet = list(), by_facet = list(), stacked = tibble()))
  }

  keep_names <- names(raw_tbls)
  if (!is.null(facets)) {
    keep_names <- intersect(keep_names, as.character(facets))
  }
  if (length(keep_names) == 0) {
    return(list(raw_by_facet = list(), by_facet = list(), stacked = tibble()))
  }
  raw_tbls <- raw_tbls[keep_names]

  by_facet <- lapply(names(raw_tbls), function(facet) {
    tbl <- raw_tbls[[facet]]
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(data.frame())
    }
    format_facets_report_gt(
      tbl = tbl,
      facet = facet,
      decimals = udecimals,
      totalscore = totalscore,
      reference = reference,
      label_style = label_style
    )
  })
  names(by_facet) <- names(raw_tbls)
  by_facet <- by_facet[vapply(by_facet, nrow, integer(1)) > 0]

  stacked <- if (length(by_facet) == 0) {
    tibble()
  } else {
    dplyr::bind_rows(
      lapply(names(by_facet), function(facet) {
        df <- by_facet[[facet]]
        df$Facet <- facet
        df
      })
    ) |>
      select("Facet", everything())
  }

  list(
    raw_by_facet = raw_tbls,
    by_facet = by_facet,
    stacked = stacked
  )
}

calc_expected_category_counts <- function(res) {
  if (is.null(res)) return(tibble())
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx,
                                criterion_splits = idx$criterion_splits)
  }
  if (length(probs) == 0) return(tibble())
  w <- idx$weight
  exp_counts <- if (is.null(w)) {
    colSums(probs, na.rm = TRUE)
  } else {
    colSums(probs * w, na.rm = TRUE)
  }
  total_exp <- sum(exp_counts)
  cat_vals <- seq(prep$rating_min, prep$rating_max)
  tibble(
    Category = cat_vals,
    ExpectedCount = exp_counts,
    ExpectedPercent = ifelse(total_exp > 0, 100 * exp_counts / total_exp, NA_real_)
  )
}

calc_category_stats <- function(obs_df, res = NULL, whexact = FALSE) {
  if (nrow(obs_df) == 0) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  total_n <- sum(obs_df$.Weight, na.rm = TRUE)
  obs_summary <- obs_df |>
    group_by(Observed) |>
    summarize(
      Count = sum(.Weight, na.rm = TRUE),
      AvgPersonMeasure = weighted_mean(PersonMeasure, .Weight),
      ExpectedAverage = weighted_mean(Expected, .Weight),
      Infit = sum(StdSq * Var * .Weight, na.rm = TRUE) / sum(Var * .Weight, na.rm = TRUE),
      Outfit = sum(StdSq * .Weight, na.rm = TRUE) / sum(.Weight, na.rm = TRUE),
      MeanResidual = weighted_mean(Residual, .Weight),
      DF_Infit = sum(Var * .Weight, na.rm = TRUE),
      DF_Outfit = sum(.Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(Category = Observed)

  all_categories <- if (!is.null(res)) {
    seq(res$prep$rating_min, res$prep$rating_max)
  } else {
    sort(unique(obs_df$Observed))
  }

  cat_tbl <- tibble(Category = all_categories) |>
    left_join(obs_summary, by = "Category") |>
    mutate(
      Count = replace_na(Count, 0),
      Percent = ifelse(total_n > 0, 100 * Count / total_n, NA_real_),
      InfitZSTD = zstd_from_mnsq(Infit, DF_Infit, whexact = whexact),
      OutfitZSTD = zstd_from_mnsq(Outfit, DF_Outfit, whexact = whexact)
    )

  exp_tbl <- calc_expected_category_counts(res)
  if (nrow(exp_tbl) > 0) {
    cat_tbl <- cat_tbl |>
      left_join(exp_tbl, by = "Category") |>
      mutate(
        DiffCount = ifelse(is.finite(ExpectedCount), Count - ExpectedCount, NA_real_),
        DiffPercent = ifelse(is.finite(ExpectedPercent), Percent - ExpectedPercent, NA_real_)
      )
  }

  cat_tbl |>
    mutate(
      LowCount = Count < 10,
      InfitFlag = !is.na(Infit) & (Infit < 0.5 | Infit > 1.5),
      OutfitFlag = !is.na(Outfit) & (Outfit < 0.5 | Outfit > 1.5),
      ZSTDFlag = (is.finite(InfitZSTD) & abs(InfitZSTD) >= 2) |
        (is.finite(OutfitZSTD) & abs(OutfitZSTD) >= 2)
    ) |>
    arrange(Category)
}

make_union_find <- function(nodes) {
  parent <- setNames(nodes, nodes)
  find_root <- function(x) {
    px <- parent[[x]]
    if (is.null(px)) return(NA_character_)
    if (px != x) {
      parent[[x]] <<- find_root(px)
    }
    parent[[x]]
  }
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (is.na(ra) || is.na(rb) || ra == rb) return(NULL)
    parent[[rb]] <<- ra
  }
  list(find = find_root, union = union_nodes)
}

calc_subsets <- function(obs_df, facet_cols) {
  if (is.null(obs_df) || nrow(obs_df) == 0 || length(facet_cols) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }
  df <- obs_df |>
    select(all_of(facet_cols)) |>
    filter(if_all(everything(), ~ !is.na(.)))
  if (nrow(df) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }

  nodes <- unique(unlist(lapply(facet_cols, function(facet) {
    paste0(facet, ":", as.character(df[[facet]]))
  })))
  if (length(nodes) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }
  uf <- make_union_find(nodes)

  for (i in seq_len(nrow(df))) {
    row_vals <- unlist(df[i, facet_cols, drop = FALSE], use.names = FALSE)
    row_nodes <- paste0(facet_cols, ":", as.character(row_vals))
    row_nodes <- row_nodes[!is.na(row_nodes)]
    if (length(row_nodes) < 2) next
    base <- row_nodes[1]
    for (node in row_nodes[-1]) {
      uf$union(base, node)
    }
  }

  comp_ids <- vapply(nodes, uf$find, character(1))
  comp_levels <- unique(comp_ids)
  comp_index <- setNames(seq_along(comp_levels), comp_levels)

  node_tbl <- tibble(
    Node = nodes,
    Component = comp_ids,
    Subset = comp_index[comp_ids],
    Facet = stringr::str_extract(nodes, "^[^:]+"),
    Level = stringr::str_replace(nodes, "^[^:]+:", "")
  )

  facet_counts <- node_tbl |>
    group_by(Subset, Facet) |>
    summarize(Levels = n_distinct(Level), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Facet, values_from = Levels, values_fill = 0)

  row_subset <- vapply(seq_len(nrow(df)), function(i) {
    node <- paste0(facet_cols[1], ":", as.character(df[[facet_cols[1]]][i]))
    comp_index[uf$find(node)]
  }, integer(1))
  obs_counts <- tibble(Subset = row_subset) |>
    count(Subset, name = "Observations")

  summary_tbl <- facet_counts |>
    left_join(obs_counts, by = "Subset") |>
    mutate(Observations = replace_na(Observations, 0)) |>
    arrange(desc(Observations))

  list(summary = summary_tbl, nodes = node_tbl)
}

calc_step_order <- function(step_tbl) {
  if (is.null(step_tbl) || nrow(step_tbl) == 0) return(tibble())
  step_tbl <- step_tbl |>
    mutate(StepIndex = suppressWarnings(as.integer(stringr::str_extract(Step, "\\d+"))))
  if (!"StepFacet" %in% names(step_tbl)) {
    step_tbl <- mutate(step_tbl, StepFacet = "Common")
  }
  step_tbl |>
    arrange(StepFacet, StepIndex) |>
    group_by(StepFacet) |>
    mutate(
      Spacing = Estimate - lag(Estimate),
      Ordered = ifelse(is.na(Spacing), NA, Spacing > 0)
    ) |>
    ungroup()
}

category_warnings_text <- function(cat_tbl, step_tbl = NULL) {
  if (is.null(cat_tbl) || nrow(cat_tbl) == 0) return("No category diagnostics available.")
  msgs <- character()
  unused <- cat_tbl |>
    filter(Count == 0)
  if (nrow(unused) > 0) {
    msgs <- c(msgs, paste0("Unused categories: ", paste(unused$Category, collapse = ", ")))
  }
  low_counts <- cat_tbl |>
    filter(Count < 10)
  if (nrow(low_counts) > 0) {
    msgs <- c(msgs, paste0("Low category counts (<10): ", paste(low_counts$Category, collapse = ", ")))
  }
  if ("DiffPercent" %in% names(cat_tbl)) {
    diff_bad <- cat_tbl |>
      filter(is.finite(DiffPercent), abs(DiffPercent) >= 5)
    if (nrow(diff_bad) > 0) {
      msgs <- c(msgs, paste0("Observed vs expected % differs by >= 5: ", paste(diff_bad$Category, collapse = ", ")))
    }
  }
  if ("InfitZSTD" %in% names(cat_tbl) && "OutfitZSTD" %in% names(cat_tbl)) {
    zstd_bad <- cat_tbl |>
      filter((is.finite(InfitZSTD) & abs(InfitZSTD) >= 2) | (is.finite(OutfitZSTD) & abs(OutfitZSTD) >= 2))
    if (nrow(zstd_bad) > 0) {
      msgs <- c(msgs, paste0("Category |ZSTD| >= 2: ", paste(zstd_bad$Category, collapse = ", ")))
    }
  }
  avg_tbl <- cat_tbl |>
    filter(is.finite(AvgPersonMeasure)) |>
    arrange(Category)
  if (nrow(avg_tbl) >= 3 && is.unsorted(avg_tbl$AvgPersonMeasure, strictly = FALSE)) {
    msgs <- c(msgs, "Category averages are not monotonic (Avg Measure by category).")
  }
  if (!is.null(step_tbl) && nrow(step_tbl) > 0) {
    disordered <- step_tbl |>
      filter(!is.na(Ordered), Ordered == FALSE)
    if (nrow(disordered) > 0) {
      bad <- disordered |>
        mutate(Label = paste0(StepFacet, ":", Step)) |>
        pull(Label)
      msgs <- c(msgs, paste0("Disordered thresholds detected: ", paste(bad, collapse = ", ")))
    }
  }
  if (length(msgs) == 0) {
    "No major category warnings detected."
  } else {
    paste(msgs, collapse = "\n")
  }
}

get_extreme_levels <- function(obs_df, facet_names, rating_min, rating_max) {
  out <- list()
  for (facet in facet_names) {
    if (!facet %in% names(obs_df)) {
      out[[facet]] <- character(0)
      next
    }
    stat <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        .groups = "drop"
      )
    extreme <- stat |>
      filter(
        (MinScore == rating_min & MaxScore == rating_min) |
          (MinScore == rating_max & MaxScore == rating_max)
      )
    out[[facet]] <- as.character(extreme[[facet]])
  }
  out
}

estimate_bias_interaction <- function(res,
                                      diagnostics,
                                      facet_a = NULL,
                                      facet_b = NULL,
                                      interaction_facets = NULL,
                                      max_abs = 10,
                                      omit_extreme = TRUE,
                                      max_iter = 4,
                                      tol = 1e-3) {
  if (is.null(res) || is.null(diagnostics)) return(list())
  obs_df <- diagnostics$obs
  if (is.null(obs_df) || nrow(obs_df) == 0) return(list())

  facet_names <- c("Person", res$config$facet_names)
  selected_facets <- if (!is.null(interaction_facets)) {
    as.character(interaction_facets)
  } else if (!is.null(facet_a) && !is.null(facet_b)) {
    c(as.character(facet_a[1]), as.character(facet_b[1]))
  } else {
    character(0)
  }
  selected_facets <- selected_facets[!is.na(selected_facets) & nzchar(selected_facets)]
  selected_facets <- unique(selected_facets)
  if (length(selected_facets) < 2) return(list())
  if (!all(selected_facets %in% facet_names)) return(list())
  if (!all(selected_facets %in% names(obs_df))) return(list())

  prep <- res$prep
  config <- res$config
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  idx <- build_indices(prep, step_facet = config$step_facet)
  theta_hat <- if (config$method == "JMLE") params$theta else res$facets$person$Estimate
  eta_base <- compute_eta(idx, params, config, theta_override = theta_hat)
  score_k <- idx$score_k
  weight <- idx$weight
  step_idx <- idx$step_idx

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    step_cum_mat <- NULL
  } else {
    step_cum <- NULL
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
  }

  if (omit_extreme) {
    extreme_levels <- get_extreme_levels(
      obs_df,
      selected_facets,
      prep$rating_min,
      prep$rating_max
    )
  } else {
    extreme_levels <- list()
  }

  measures <- diagnostics$measures
  meas_map <- list()
  se_map <- list()
  if (!is.null(measures) && nrow(measures) > 0) {
    for (i in seq_len(nrow(measures))) {
      key <- paste(measures$Facet[i], as.character(measures$Level[i]), sep = "||")
      meas_map[[key]] <- measures$Estimate[i]
      se_map[[key]] <- measures$SE[i]
    }
  }

  level_map <- lapply(facet_names, function(facet) {
    if (facet == "Person") {
      as.character(prep$levels$Person)
    } else {
      as.character(prep$levels[[facet]])
    }
  })
  names(level_map) <- facet_names

  interaction_parts <- lapply(selected_facets, function(f) as.character(obs_df[[f]]))
  names(interaction_parts) <- selected_facets
  group_keys <- do.call(interaction, c(interaction_parts, list(drop = TRUE, sep = "||", lex.order = TRUE)))
  group_indices <- split(seq_len(nrow(obs_df)), group_keys)

  groups <- list()
  for (key in names(group_indices)) {
    lvl_vals <- strsplit(as.character(key), "\\|\\|")[[1]]
    if (length(lvl_vals) < length(selected_facets)) next
    lvl_vals <- as.character(lvl_vals[seq_along(selected_facets)])
    names(lvl_vals) <- selected_facets
    if (isTRUE(omit_extreme)) {
      is_extreme <- vapply(selected_facets, function(f) {
        lvl_vals[[f]] %in% (extreme_levels[[f]] %||% character(0))
      }, logical(1))
      if (any(is_extreme)) next
    }
    idx_rows <- group_indices[[key]]
    if (length(idx_rows) == 0) next
    groups[[as.character(key)]] <- list(levels = lvl_vals, idx = idx_rows)
  }
  if (length(groups) == 0) return(list())

  orientation_audit <- audit_interaction_orientation(config, selected_facets)
  has_bias_error <- function(x) {
    is.character(x) && length(x) > 0L && isTRUE(nzchar(x[1]))
  }

  estimate_bias_for_group <- function(idx_rows) {
    eta_sub <- eta_base[idx_rows]
    score_k_sub <- score_k[idx_rows]
    weight_sub <- if (!is.null(weight)) weight[idx_rows] else NULL
    step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL

    if (config$model == "RSM") {
      nll <- function(b) -loglik_rsm(eta_sub + b, score_k_sub, step_cum, weight = weight_sub)
    } else {
      nll <- function(b) -loglik_pcm(eta_sub + b, score_k_sub, step_cum_mat, step_idx_sub, weight = weight_sub)
    }
    opt <- tryCatch(
      stats::optimize(nll, interval = c(-max_abs, max_abs)),
      error = function(e) e
    )
    if (inherits(opt, "error")) {
      return(list(value = NA_real_, error = conditionMessage(opt)))
    }
    list(value = opt$minimum, error = NULL)
  }

  iteration_metrics <- function(bias_map) {
    max_resid <- 0
    max_resid_pct <- NA_real_
    for (key in names(groups)) {
      g <- groups[[key]]
      idx_rows <- g$idx
      bias <- bias_map[[key]]
      if (!is.finite(bias) || has_bias_error(bias_error_map[[key]])) next
      eta_sub <- eta_base[idx_rows] + ifelse(is.finite(bias), bias, 0)
      score_k_sub <- score_k[idx_rows]
      step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL
      probs <- if (config$model == "RSM") {
        category_prob_rsm(eta_sub, step_cum)
      } else {
        category_prob_pcm(eta_sub, step_cum_mat, step_idx_sub)
      }
      k_vals <- 0:(ncol(probs) - 1)
      expected_k <- as.vector(probs %*% k_vals)
      expected_score <- prep$rating_min + expected_k
      obs_score <- obs_df$Observed[idx_rows]
      if ("Weight" %in% names(obs_df)) {
        w <- obs_df$Weight[idx_rows]
        obs_score <- obs_score * w
        expected_score <- expected_score * w
      }
      obs_sum <- sum(obs_score, na.rm = TRUE)
      exp_sum <- sum(expected_score, na.rm = TRUE)
      resid_sum <- obs_sum - exp_sum
      if (abs(resid_sum) >= abs(max_resid)) {
        max_resid <- resid_sum
        max_resid_pct <- ifelse(exp_sum != 0, resid_sum / exp_sum * 100, NA_real_)
      }
    }
    list(
      max_resid = max_resid,
      max_resid_pct = max_resid_pct,
      max_resid_categories = NA_real_
    )
  }

  bias_map <- stats::setNames(as.list(rep(0, length(groups))), names(groups))
  bias_error_map <- stats::setNames(as.list(rep(NA_character_, length(groups))), names(groups))

  iter_rows <- list()
  for (it in seq_len(max_iter)) {
    max_change_abs <- 0
    max_change_signed <- 0
    changes <- numeric(0)
    for (key in names(groups)) {
      g <- groups[[key]]
      bias_result <- estimate_bias_for_group(g$idx)
      bias_hat <- bias_result$value
      prev <- bias_map[[key]]
      if (is.finite(bias_hat) && is.finite(prev)) {
        delta <- bias_hat - prev
        if (abs(delta) >= max_change_abs) {
          max_change_abs <- abs(delta)
          max_change_signed <- delta
        }
        changes <- c(changes, abs(delta))
      } else {
        changes <- c(changes, NA_real_)
      }
      bias_map[[key]] <- bias_hat
      bias_error_map[[key]] <- if (has_bias_error(bias_result$error)) bias_result$error else NA_character_
    }
    resid_info <- iteration_metrics(bias_map)
    iter_rows[[it]] <- tibble(
      Iteration = it,
      MaxScoreResidual = resid_info$max_resid,
      MaxScoreResidualPct = resid_info$max_resid_pct,
      MaxScoreResidualCategories = resid_info$max_resid_categories,
      MaxLogitChange = max_change_signed,
      BiasCells = sum(changes > tol, na.rm = TRUE)
    )
    if (max_change_abs < tol) break
  }

  interaction_label <- paste(selected_facets, collapse = " x ")
  interaction_order <- length(selected_facets)
  interaction_mode <- ifelse(interaction_order > 2, "higher_order", "pairwise")

  rows <- list()
  seq_id <- 1
  for (key in names(groups)) {
    g <- groups[[key]]
    idx_rows <- g$idx
    levels <- g$levels
    bias_hat <- bias_map[[key]]
    bias_error <- bias_error_map[[key]]
    bias_ok <- is.finite(bias_hat) && !has_bias_error(bias_error)
    eta_sub <- eta_base[idx_rows]
    score_k_sub <- score_k[idx_rows]
    weight_sub <- if (!is.null(weight)) weight[idx_rows] else rep(1, length(idx_rows))
    step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL

    if (bias_ok) {
      probs <- if (config$model == "RSM") {
        category_prob_rsm(eta_sub + bias_hat, step_cum)
      } else {
        category_prob_pcm(eta_sub + bias_hat, step_cum_mat, step_idx_sub)
      }
      k_vals <- 0:(ncol(probs) - 1)
      expected_k <- as.vector(probs %*% k_vals)
      var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
      var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
      resid_k <- score_k_sub - expected_k
      std_sq <- resid_k^2 / var_k

      info <- sum(var_k * weight_sub, na.rm = TRUE)
      se <- ifelse(is.finite(info) && info > 0, 1 / sqrt(info), NA_real_)
      infit <- ifelse(sum(var_k * weight_sub, na.rm = TRUE) > 0,
                      sum(std_sq * var_k * weight_sub, na.rm = TRUE) / sum(var_k * weight_sub, na.rm = TRUE),
                      NA_real_)
      outfit <- ifelse(sum(weight_sub, na.rm = TRUE) > 0,
                       sum(std_sq * weight_sub, na.rm = TRUE) / sum(weight_sub, na.rm = TRUE),
                       NA_real_)
    } else {
      expected_k <- rep(NA_real_, length(score_k_sub))
      se <- NA_real_
      infit <- NA_real_
      outfit <- NA_real_
    }

    obs_slice <- obs_df[idx_rows, , drop = FALSE]
    w_obs <- if ("Weight" %in% names(obs_slice)) obs_slice$Weight else rep(1, nrow(obs_slice))
    obs_score <- sum(obs_slice$Observed * w_obs, na.rm = TRUE)
    exp_score <- if (bias_ok) sum((prep$rating_min + expected_k) * w_obs, na.rm = TRUE) else NA_real_
    obs_count <- sum(w_obs, na.rm = TRUE)
    obs_exp_avg <- ifelse(bias_ok && obs_count > 0 && is.finite(exp_score), (obs_score - exp_score) / obs_count, NA_real_)

    n_obs <- nrow(obs_slice)
    df_t <- ifelse(is.finite(obs_count), max(obs_count - 1, 0), NA_real_)
    t_val <- ifelse(is.finite(bias_hat) && is.finite(se) && se > 0, bias_hat / se, NA_real_)
    p_val <- ifelse(is.finite(t_val) && df_t > 0, 2 * stats::pt(-abs(t_val), df = df_t), NA_real_)

    row <- list(
      Sq = seq_id,
      `Observd Score` = obs_score,
      `Expctd Score` = exp_score,
      `Observd Count` = obs_count,
      `Obs-Exp Average` = obs_exp_avg,
      `Bias Size` = bias_hat,
      `S.E.` = se,
      t = t_val,
      `d.f.` = df_t,
      `Prob.` = p_val,
      InferenceTier = "screening",
      SupportsFormalInference = FALSE,
      FormalInferenceEligible = FALSE,
      PrimaryReportingEligible = FALSE,
      ReportingUse = "screening_only",
      SEBasis = "conditional plug-in information",
      DFBasis = "observed-count minus 1 approximation",
      StatisticLabel = "screening t",
      ProbabilityMetric = "screening tail area",
      Infit = infit,
      Outfit = outfit,
      ObsN = n_obs,
      OptimizationStatus = if (bias_ok) "ok" else "failed",
      OptimizationDetail = if (bias_ok) NA_character_ else bias_error,
      InteractionFacets = interaction_label,
      InteractionOrder = interaction_order,
      InteractionMode = interaction_mode
    )

    for (j in seq_along(selected_facets)) {
      facet_j <- selected_facets[j]
      level_j <- as.character(levels[[facet_j]])
      index_j <- match(level_j, level_map[[facet_j]])
      row[[paste0("Facet", j)]] <- facet_j
      row[[paste0("Facet", j, "_Level")]] <- level_j
      row[[paste0("Facet", j, "_Index")]] <- ifelse(is.na(index_j), NA_real_, index_j)
      row[[paste0("Facet", j, "_Measure")]] <- meas_map[[paste(facet_j, level_j, sep = "||")]]
      row[[paste0("Facet", j, "_SE")]] <- se_map[[paste(facet_j, level_j, sep = "||")]]
    }

    # Backward compatibility for existing 2-way workflows and helper APIs.
    if (interaction_order >= 2) {
      row$FacetA <- row$Facet1
      row$FacetA_Level <- row$Facet1_Level
      row$FacetA_Index <- row$Facet1_Index
      row$FacetA_Measure <- row$Facet1_Measure
      row$FacetA_SE <- row$Facet1_SE

      row$FacetB <- row$Facet2
      row$FacetB_Level <- row$Facet2_Level
      row$FacetB_Index <- row$Facet2_Index
      row$FacetB_Measure <- row$Facet2_Measure
      row$FacetB_SE <- row$Facet2_SE
    }

    rows[[seq_id]] <- tibble::as_tibble(row)
    seq_id <- seq_id + 1
  }

  bias_tbl <- bind_rows(rows)
  if (nrow(bias_tbl) == 0) return(list())

  numeric_cols <- c("Observd Score", "Expctd Score", "Observd Count", "Obs-Exp Average", "Bias Size", "S.E.")
  pop_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    m <- mean(x)
    sqrt(mean((x - m)^2))
  }
  mean_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)))
  sd_pop_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ pop_sd(.x)))
  sd_sample_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ stats::sd(.x, na.rm = TRUE)))
  summary_tbl <- bind_rows(
    mean_row,
    sd_pop_row,
    sd_sample_row
  ) |>
    mutate(
      Statistic = c(paste0("Mean (Count: ", nrow(bias_tbl), ")"), "S.D. (Population)", "S.D. (Sample)"),
      InteractionFacets = interaction_label,
      InteractionOrder = interaction_order,
      InteractionMode = interaction_mode,
      MixedSign = orientation_audit$mixed_sign,
      .before = 1
    )

  se_bias <- bias_tbl$`S.E.`
  bias_vals <- bias_tbl$`Bias Size`
  w_chi <- ifelse(is.finite(se_bias) & se_bias > 0, 1 / (se_bias^2), NA_real_)
  ok <- is.finite(w_chi) & is.finite(bias_vals)
  fixed_chi <- if (sum(ok) >= 2) {
    sum(w_chi[ok] * bias_vals[ok]^2) - (sum(w_chi[ok] * bias_vals[ok])^2) / sum(w_chi[ok])
  } else {
    NA_real_
  }
  fixed_df <- max(nrow(bias_tbl) - 1, 0)
  fixed_prob <- ifelse(is.finite(fixed_chi) && fixed_df > 0, 1 - stats::pchisq(fixed_chi, df = fixed_df), NA_real_)
  optimization_failures <- purrr::imap_dfr(groups, function(g, key) {
    err <- bias_error_map[[key]]
    if (!has_bias_error(err)) {
      return(tibble::tibble())
    }
    out <- tibble::tibble(
      InteractionKey = key,
      Error = err
    )
    for (facet_j in selected_facets) {
      out[[facet_j]] <- as.character(g$levels[[facet_j]])
    }
    out
  })

  chi_tbl <- tibble(
    FixedChiSq = fixed_chi,
    FixedDF = fixed_df,
    FixedProb = fixed_prob,
    InferenceTier = "screening",
    SupportsFormalInference = FALSE,
    FormalInferenceEligible = FALSE,
    PrimaryReportingEligible = FALSE,
    ReportingUse = "screening_only",
    TestBasis = "conditional plug-in heterogeneity screen",
    InteractionFacets = interaction_label,
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode
  )

  list(
    facet_a = selected_facets[1],
    facet_b = selected_facets[2],
    interaction_facets = selected_facets,
    interaction_order = interaction_order,
    interaction_mode = interaction_mode,
    table = bias_tbl,
    summary = summary_tbl,
    chi_sq = chi_tbl,
    iteration = bind_rows(iter_rows),
    orientation_audit = orientation_audit$table,
    mixed_sign = orientation_audit$mixed_sign,
    direction_note = orientation_audit$direction_note,
    recommended_action = orientation_audit$recommended_action,
    inference_tier = "screening",
    optimization_failures = optimization_failures
  )
}

calc_bias_pairwise <- function(bias_tbl, target_facet, context_facet) {
  if (is.null(bias_tbl) || nrow(bias_tbl) == 0) return(tibble())
  if ("InteractionOrder" %in% names(bias_tbl)) {
    ord <- suppressWarnings(as.numeric(bias_tbl$InteractionOrder[1]))
    if (is.finite(ord) && ord != 2) return(tibble())
  }

  use_a <- bias_tbl$FacetA[1] == target_facet
  use_b <- bias_tbl$FacetB[1] == target_facet
  if (!(use_a || use_b)) return(tibble())

  target_prefix <- if (use_a) "FacetA" else "FacetB"
  context_prefix <- if (use_a) "FacetB" else "FacetA"

  sub <- bias_tbl |>
    filter(.data[[context_prefix]] %in% context_facet)
  if (nrow(sub) == 0) sub <- bias_tbl

  rows <- list()
  for (tgt_level in unique(sub[[paste0(target_prefix, "_Level")]])) {
    group_df <- sub |> filter(.data[[paste0(target_prefix, "_Level")]] == tgt_level)
    contexts <- unique(group_df[[paste0(context_prefix, "_Level")]])
    if (length(contexts) < 2) next
    ctx_pairs <- combn(contexts, 2, simplify = FALSE)
    for (pair in ctx_pairs) {
      c1 <- pair[1]
      c2 <- pair[2]
      r1 <- group_df |> filter(.data[[paste0(context_prefix, "_Level")]] == c1) |> slice(1)
      r2 <- group_df |> filter(.data[[paste0(context_prefix, "_Level")]] == c2) |> slice(1)

      tgt_measure <- r1[[paste0(target_prefix, "_Measure")]]
      tgt_se <- r1[[paste0(target_prefix, "_SE")]]
      tgt_index <- r1[[paste0(target_prefix, "_Index")]]

      bias1 <- r1$`Bias Size`
      bias2 <- r2$`Bias Size`
      bias_se1 <- r1$`S.E.`
      bias_se2 <- r2$`S.E.`

      local1 <- ifelse(is.finite(tgt_measure) & is.finite(bias1), tgt_measure + bias1, NA_real_)
      local2 <- ifelse(is.finite(tgt_measure) & is.finite(bias2), tgt_measure + bias2, NA_real_)
      se1 <- ifelse(is.finite(tgt_se) & is.finite(bias_se1), sqrt(tgt_se^2 + bias_se1^2), NA_real_)
      se2 <- ifelse(is.finite(tgt_se) & is.finite(bias_se2), sqrt(tgt_se^2 + bias_se2^2), NA_real_)

      contrast <- ifelse(is.finite(local1) & is.finite(local2), local1 - local2, NA_real_)
      se_contrast <- ifelse(
        is.finite(bias_se1) & is.finite(bias_se2),
        sqrt(bias_se1^2 + bias_se2^2),
        NA_real_
      )
      t_val <- ifelse(is.finite(contrast) & is.finite(se_contrast) & se_contrast > 0, contrast / se_contrast, NA_real_)

      n1 <- ifelse(is.finite(r1$ObsN), r1$ObsN, 0)
      n2 <- ifelse(is.finite(r2$ObsN), r2$ObsN, 0)
      df_t <- welch_satterthwaite_df(c(bias_se1^2, bias_se2^2), c(n1 - 1, n2 - 1))
      p_val <- ifelse(is.finite(t_val) & is.finite(df_t) & df_t > 0,
                      2 * stats::pt(-abs(t_val), df = df_t),
                      NA_real_)

      rows[[length(rows) + 1]] <- tibble(
        Target = tgt_level,
        `Target N` = tgt_index,
        `Target Measure` = tgt_measure,
        `Target S.E.` = tgt_se,
        Context1 = c1,
        `Context1 N` = r1[[paste0(context_prefix, "_Index")]],
        `Local Measure1` = local1,
        SE1 = se1,
        `Obs-Exp Avg1` = r1$`Obs-Exp Average`,
        Count1 = r1$`Observd Count`,
        ObsN1 = n1,
        Context2 = c2,
        `Context2 N` = r2[[paste0(context_prefix, "_Index")]],
        `Local Measure2` = local2,
        SE2 = se2,
        `Obs-Exp Avg2` = r2$`Obs-Exp Average`,
        Count2 = r2$`Observd Count`,
        ObsN2 = n2,
        Contrast = contrast,
        SE = se_contrast,
        t = t_val,
        `d.f.` = df_t,
        `Prob.` = p_val,
        InferenceTier = "screening",
        SupportsFormalInference = FALSE,
        FormalInferenceEligible = FALSE,
        PrimaryReportingEligible = FALSE,
        ReportingUse = "screening_only",
        ContrastBasis = "difference between local target measures across contexts (target term cancels to a bias contrast)",
        SEBasis = "combined context-specific bias standard errors",
        StatisticLabel = "Bias-contrast Welch screening t",
        ProbabilityMetric = "screening tail area",
        DFBasis = "Welch-Satterthwaite approximation"
      )
    }
  }
  bind_rows(rows)
}

extract_anchor_tables <- function(config) {
  facet_names <- c("Person", config$facet_names)
  anchor_tbl <- purrr::map_dfr(facet_names, function(facet) {
    spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
    anchors <- spec$anchors
    if (is.null(anchors)) return(tibble())
    df <- tibble(Facet = facet, Level = names(anchors), Anchor = as.numeric(anchors)) |>
      filter(is.finite(Anchor))
    if (nrow(df) == 0) return(tibble())
    df |>
      mutate(Source = ifelse(facet %in% config$dummy_facets, "Dummy facet", "Anchor"))
  })

  group_tbl <- purrr::map_dfr(facet_names, function(facet) {
    spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
    groups <- spec$groups
    if (is.null(groups)) return(tibble())
    df <- tibble(Facet = facet, Level = names(groups), Group = as.character(groups)) |>
      filter(!is.na(Group), Group != "")
    if (nrow(df) == 0) return(tibble())
    group_values <- spec$group_values
    df |>
      mutate(GroupValue = ifelse(Group %in% names(group_values), group_values[Group], NA_real_))
  })

  list(anchors = anchor_tbl, groups = group_tbl)
}

calc_facets_chisq <- function(measure_df) {
  if (is.null(measure_df) || nrow(measure_df) == 0) return(tibble())
  measure_df |>
    group_by(Facet) |>
    summarize(
      Levels = n(),
      MeanMeasure = mean(Estimate, na.rm = TRUE),
      SD = sd(Estimate, na.rm = TRUE),
      EstimateVec = list(Estimate),
      SEVec = list(SE),
      FixedChiSq = {
        w <- ifelse(is.finite(SE) & SE > 0, 1 / (SE^2), NA_real_)
        d <- Estimate
        ok <- is.finite(w) & is.finite(d)
        if (sum(ok) < 2) NA_real_ else sum(w[ok] * d[ok]^2) - (sum(w[ok] * d[ok])^2) / sum(w[ok])
      },
      FixedDF = pmax(Levels - 1, 0),
      RandomVar = {
        d <- Estimate
        se2 <- SE^2
        ok <- is.finite(d) & is.finite(se2)
        if (sum(ok) < 2) NA_real_ else (sum((d[ok] - mean(d[ok]))^2) / (sum(ok) - 1)) - (sum(se2[ok]) / sum(ok))
      },
      .groups = "drop"
    ) |>
    mutate(
      FixedProb = ifelse(is.finite(FixedChiSq) & FixedDF > 0,
                         1 - stats::pchisq(FixedChiSq, df = FixedDF),
                         NA_real_),
      RandomVar = ifelse(is.finite(RandomVar) & RandomVar > 0, RandomVar, NA_real_)
    ) |>
    rowwise() |>
    mutate(
      RandomChiSq = {
        if (!is.finite(RandomVar) || RandomVar <= 0) {
          NA_real_
        } else {
          d <- unlist(EstimateVec)
          se <- unlist(SEVec)
          w <- ifelse(is.finite(se) & se > 0, 1 / (RandomVar + se^2), NA_real_)
          ok <- is.finite(w) & is.finite(d)
          if (sum(ok) < 2) NA_real_ else sum(w[ok] * d[ok]^2) - (sum(w[ok] * d[ok])^2) / sum(w[ok])
        }
      },
      RandomDF = pmax(Levels - 2, 0),
      RandomProb = ifelse(is.finite(RandomChiSq) & RandomDF > 0,
                          1 - stats::pchisq(RandomChiSq, df = RandomDF),
                          NA_real_)
    ) |>
    ungroup() |>
    select(-EstimateVec, -SEVec)
}

build_param_slices <- function(sizes) {
  out <- list()
  idx <- 1L
  for (nm in names(sizes)) {
    k <- as.integer(sizes[[nm]] %||% 0L)
    if (k > 0L) {
      out[[nm]] <- idx:(idx + k - 1L)
      idx <- idx + k
    } else {
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

  baseline <- expand_facet_with_constraints(numeric(n_params), spec)
  for (j in seq_len(n_params)) {
    free_vec <- numeric(n_params)
    free_vec[j] <- 1
    out[, j] <- expand_facet_with_constraints(free_vec, spec) - baseline
  }
  out
}

symmetrize_matrix <- function(mat) {
  if (is.null(mat)) return(NULL)
  (mat + t(mat)) / 2
}

invert_information_matrix <- function(info_mat) {
  if (is.null(info_mat)) {
    return(list(cov = NULL, regularized = FALSE, rank = 0L))
  }

  info_mat <- suppressWarnings(as.matrix(info_mat))
  if (!is.numeric(info_mat) || length(info_mat) == 0L || any(!is.finite(info_mat))) {
    return(list(cov = NULL, regularized = FALSE, rank = 0L))
  }

  info_mat <- symmetrize_matrix(info_mat)
  eig <- tryCatch(
    eigen(info_mat, symmetric = TRUE),
    error = function(e) NULL
  )
  if (is.null(eig)) {
    return(list(cov = NULL, regularized = FALSE, rank = 0L))
  }

  values <- eig$values
  max_val <- suppressWarnings(max(abs(values), na.rm = TRUE))
  tol <- if (is.finite(max_val) && max_val > 0) {
    max_val * sqrt(.Machine$double.eps)
  } else {
    sqrt(.Machine$double.eps)
  }
  rank <- sum(values > tol)
  regularized <- any(values <= tol)
  safe_values <- pmax(values, tol)
  cov_mat <- eig$vectors %*% diag(1 / safe_values, nrow = length(safe_values)) %*% t(eig$vectors)
  cov_mat <- symmetrize_matrix(cov_mat)

  list(
    cov = cov_mat,
    regularized = regularized,
    rank = rank
  )
}

compute_mml_facet_model_se <- function(res) {
  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  if (!identical(method, "MML")) {
    return(list(
      table = tibble(),
      status = "not_applicable",
      detail = "MML observed-information SEs are only computed for MML fits."
    ))
  }

  if (is.null(res$opt$par) || length(res$opt$par) == 0L) {
    return(list(
      table = tibble(),
      status = "fallback",
      detail = "Optimized parameter vector was not available; fell back to observation-table SEs."
    ))
  }

  config <- res$config
  sizes <- build_param_sizes(config)
  idx <- build_indices(res$prep, step_facet = config$step_facet)
  quad_points <- max(1L, as.integer(config$estimation_control$quad_points %||% 15L))
  quad <- gauss_hermite_normal(quad_points)
  cache <- make_param_cache(sizes, config, idx, is_mml = TRUE)

  fn <- function(par, idx, config, sizes, quad) {
    cache$ensure(par)
    mfrm_loglik_mml_cached(cache, idx, config, quad)
  }
  gr <- function(par, idx, config, sizes, quad) {
    cache$ensure(par)
    mfrm_grad_mml_cached(cache, idx, config, sizes, quad)
  }

  hess <- tryCatch(
    stats::optimHess(
      par = res$opt$par,
      fn = fn,
      gr = gr,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad
    ),
    error = function(e) e
  )

  if (inherits(hess, "error") || is.null(hess)) {
    msg <- if (inherits(hess, "error")) conditionMessage(hess) else "Unknown Hessian error."
    return(list(
      table = tibble(),
      status = "fallback",
      detail = paste0("Observed-information Hessian was unavailable; fell back to observation-table SEs. ", msg)
    ))
  }

  inv_info <- invert_information_matrix(hess)
  cov_free <- inv_info$cov
  if (is.null(cov_free)) {
    return(list(
      table = tibble(),
      status = "fallback",
      detail = "Observed-information matrix could not be inverted; fell back to observation-table SEs."
    ))
  }

  param_slices <- build_param_slices(sizes)
  facet_tbl <- purrr::map_dfr(config$facet_names, function(facet) {
    spec <- config$facet_specs[[facet]]
    levels <- as.character(spec$levels %||% res$prep$levels[[facet]] %||% character(0))
    slice <- param_slices[[facet]] %||% integer(0)

    if (length(levels) == 0L) {
      return(tibble(Facet = character(0), Level = character(0), ModelSE = numeric(0)))
    }

    if (length(slice) == 0L || is.null(spec) || (spec$n_params %||% 0L) == 0L) {
      return(tibble(
        Facet = facet,
        Level = levels,
        ModelSE = NA_real_
      ))
    }

    jac <- constraint_jacobian(spec)
    cov_block <- cov_free[slice, slice, drop = FALSE]
    cov_expanded <- symmetrize_matrix(jac %*% cov_block %*% t(jac))
    diag_var <- diag(cov_expanded)
    diag_var[diag_var < 0 & abs(diag_var) < 1e-10] <- 0

    tibble(
      Facet = facet,
      Level = levels,
      ModelSE = ifelse(diag_var >= 0, sqrt(diag_var), NA_real_)
    )
  })

  detail <- if (isTRUE(inv_info$regularized)) {
    "MML facet ModelSE values use the observed information of the marginal log-likelihood; a near-singular Hessian was regularized during inversion."
  } else {
    "MML facet ModelSE values use the observed information of the marginal log-likelihood."
  }

  list(
    table = facet_tbl,
    status = if (isTRUE(inv_info$regularized)) "regularized" else "ok",
    detail = detail
  )
}

build_measure_se_table <- function(res, obs_df, facet_cols, fit_tbl) {
  approx_tbl <- calc_facet_se(obs_df, facet_cols) |>
    rename(ApproxSE = SE)

  base_tbl <- bind_rows(
    res$facets$person |>
      transmute(Facet = "Person", Level = as.character(Person)),
    res$facets$others |>
      transmute(Facet = as.character(Facet), Level = as.character(Level))
  ) |>
    distinct(Facet, Level)

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  converged <- isTRUE(as.logical(res$summary$Converged[1] %||% FALSE))
  fallback_label <- if (identical(method, "MML")) {
    "Fallback observation-table information"
  } else {
    "Observation-table information"
  }

  person_model_tbl <- if (identical(method, "MML") && "SD" %in% names(res$facets$person)) {
    res$facets$person |>
      transmute(
        Facet = "Person",
        Level = as.character(Person),
        ModelSE = as.numeric(SD),
        SE_Method = "Posterior SD (EAP)"
      )
  } else {
    tibble(
      Facet = character(0),
      Level = character(0),
      ModelSE = numeric(0),
      SE_Method = character(0)
    )
  }

  facet_model_bundle <- compute_mml_facet_model_se(res)
  facet_model_tbl <- facet_model_bundle$table |>
    mutate(SE_Method = "Observed information (MML)")

  base_tbl |>
    left_join(approx_tbl, by = c("Facet", "Level")) |>
    left_join(
      bind_rows(person_model_tbl, facet_model_tbl) |>
        rename(ModelSE_Raw = ModelSE, SE_Method_Raw = SE_Method),
      by = c("Facet", "Level")
    ) |>
    left_join(
      fit_tbl |>
        select("Facet", "Level", "Infit"),
      by = c("Facet", "Level")
    ) |>
    mutate(
      PrecisionTier = dplyr::case_when(
        identical(method, "MML") & !is.na(.data$SE_Method_Raw) ~ "model_based",
        identical(method, "MML") ~ "hybrid",
        TRUE ~ "exploratory"
      ),
      Converged = converged,
      SupportsFormalInference = .data$PrecisionTier == "model_based" & .data$Converged,
      ModelSE = dplyr::coalesce(.data$ModelSE_Raw, .data$ApproxSE),
      RealSE = ifelse(
        is.finite(.data$ModelSE),
        .data$ModelSE * sqrt(pmax(dplyr::coalesce(.data$Infit, 1), 1)),
        NA_real_
      ),
      SE = .data$ModelSE,
      SE_Method = dplyr::coalesce(.data$SE_Method_Raw, fallback_label),
      SEUse = dplyr::case_when(
        .data$PrecisionTier == "model_based" & .data$Converged ~ "primary_reporting",
        .data$PrecisionTier == "model_based" ~ "review_before_reporting",
        .data$PrecisionTier == "hybrid" ~ "review_before_reporting",
        TRUE ~ "screening_only"
      ),
      CIBasis = dplyr::case_when(
        .data$PrecisionTier == "model_based" & .data$Converged ~ "Normal interval from model-based SE",
        .data$PrecisionTier == "model_based" ~ "Normal interval from model-based SE; optimizer convergence review required",
        .data$PrecisionTier == "hybrid" ~ "Normal interval from fallback observation-table SE",
        TRUE ~ "Normal interval from exploratory observation-table SE"
      ),
      CIUse = dplyr::case_when(
        .data$PrecisionTier == "model_based" & .data$Converged ~ "primary_reporting",
        .data$PrecisionTier == "model_based" ~ "review_before_reporting",
        .data$PrecisionTier == "hybrid" ~ "review_before_reporting",
        TRUE ~ "screening_only"
      )
    ) |>
    select("Facet", "Level", "N", "SE", "ModelSE", "RealSE", "SE_Method",
           "Converged",
           "PrecisionTier", "SupportsFormalInference", "SEUse", "CIBasis", "CIUse") |>
    arrange(.data$Facet, .data$Level) |>
    structure(mml_se_status = facet_model_bundle$status, mml_se_detail = facet_model_bundle$detail)
}

# Separation, reliability, and strata indices used in package diagnostics.
# Model statistics use ModelSE; fit-adjusted statistics use RealSE when available.
# Separation G = TrueSD / RMSE
# Reliability R = TrueVariance / ObservedVariance = G^2 / (1 + G^2)
# Strata H = (4*G + 1) / 3
summarize_precision_basis <- function(df, se_col, distribution_basis = c("sample", "population")) {
  distribution_basis <- match.arg(distribution_basis)

  est <- suppressWarnings(as.numeric(df$Estimate))
  se_vals <- if (se_col %in% names(df)) suppressWarnings(as.numeric(df[[se_col]])) else rep(NA_real_, nrow(df))
  ok_est <- is.finite(est)
  ok_se <- is.finite(se_vals) & se_vals >= 0
  n_est <- sum(ok_est)

  observed_mean <- if (n_est > 0L) mean(est[ok_est]) else NA_real_
  observed_var <- if (n_est == 0L) {
    NA_real_
  } else if (identical(distribution_basis, "sample")) {
    if (n_est >= 2L) stats::var(est[ok_est]) else NA_real_
  } else {
    mean((est[ok_est] - observed_mean)^2)
  }
  observed_sd <- if (is.finite(observed_var)) sqrt(observed_var) else NA_real_
  rmse <- if (any(ok_se)) sqrt(mean(se_vals[ok_se]^2, na.rm = TRUE)) else NA_real_
  error_var <- if (is.finite(rmse)) rmse^2 else NA_real_
  true_var <- if (is.finite(observed_var) && is.finite(error_var)) {
    pmax(observed_var - error_var, 0)
  } else {
    NA_real_
  }
  true_sd <- if (is.finite(true_var)) sqrt(true_var) else NA_real_
  separation <- if (is.finite(rmse) && rmse > 0 && is.finite(true_sd)) true_sd / rmse else NA_real_
  reliability <- if (is.finite(observed_var) && observed_var > 0 && is.finite(true_var)) true_var / observed_var else NA_real_
  strata <- if (is.finite(separation)) (4 * separation + 1) / 3 else NA_real_

  list(
    ObservedMean = observed_mean,
    ObservedSD = observed_sd,
    RMSE = rmse,
    TrueSD = true_sd,
    ObservedVariance = observed_var,
    ErrorVariance = error_var,
    TrueVariance = true_var,
    Separation = separation,
    Strata = strata,
    Reliability = reliability,
    SEAvailable = sum(ok_se),
    MeanSE = if (any(ok_se)) mean(se_vals[ok_se], na.rm = TRUE) else NA_real_,
    MedianSE = if (any(ok_se)) stats::median(se_vals[ok_se], na.rm = TRUE) else NA_real_
  )
}

build_facet_precision_summary <- function(measure_df, variability_tbl = NULL) {
  if (is.null(measure_df) || nrow(measure_df) == 0) return(tibble())

  variability_tbl <- as.data.frame(variability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  variability_keep <- intersect(
    c("Facet", "FixedChiSq", "FixedDF", "FixedProb", "RandomVar", "RandomChiSq", "RandomDF", "RandomProb"),
    names(variability_tbl)
  )
  variability_tbl <- variability_tbl[, variability_keep, drop = FALSE]

  mode_map <- c(model = "ModelSE", fit_adjusted = "RealSE")
  rows <- list()
  row_id <- 1L

  facets <- unique(as.character(measure_df$Facet))
  facets <- facets[!is.na(facets) & nzchar(facets)]
  for (facet in facets) {
    sub <- measure_df[as.character(measure_df$Facet) == facet, , drop = FALSE]
    for (mode_name in names(mode_map)) {
      se_col <- mode_map[[mode_name]]
      if (!se_col %in% names(sub)) next
      for (basis_name in c("sample", "population")) {
        stats <- summarize_precision_basis(sub, se_col = se_col, distribution_basis = basis_name)
        rows[[row_id]] <- tibble(
          Facet = facet,
          Levels = nrow(sub),
          DistributionBasis = basis_name,
          SEMode = mode_name,
          SEColumn = se_col,
          ObservedMean = stats$ObservedMean,
          ObservedSD = stats$ObservedSD,
          RMSE = stats$RMSE,
          TrueSD = stats$TrueSD,
          ObservedVariance = stats$ObservedVariance,
          ErrorVariance = stats$ErrorVariance,
          TrueVariance = stats$TrueVariance,
          Separation = stats$Separation,
          Strata = stats$Strata,
          Reliability = stats$Reliability,
          SEAvailable = stats$SEAvailable,
          MeanSE = stats$MeanSE,
          MedianSE = stats$MedianSE,
          MeanInfit = if ("Infit" %in% names(sub)) mean(sub$Infit, na.rm = TRUE) else NA_real_,
          MeanOutfit = if ("Outfit" %in% names(sub)) mean(sub$Outfit, na.rm = TRUE) else NA_real_
        )
        row_id <- row_id + 1L
      }
    }
  }

  out <- bind_rows(rows)
  if (nrow(out) == 0) {
    return(tibble())
  }
  if (nrow(variability_tbl) > 0 && "Facet" %in% names(variability_tbl)) {
    out <- left_join(out, variability_tbl, by = "Facet")
  }
  out |>
    arrange(.data$Facet, .data$DistributionBasis, .data$SEMode)
}

build_precision_profile <- function(res, measure_df, reliability_tbl, facet_precision_tbl) {
  measure_df <- as.data.frame(measure_df %||% data.frame(), stringsAsFactors = FALSE)
  reliability_tbl <- as.data.frame(reliability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  facet_precision_tbl <- as.data.frame(facet_precision_tbl %||% data.frame(), stringsAsFactors = FALSE)

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  method <- ifelse(identical(method, "JMLE"), "JML", method)
  converged <- isTRUE(as.logical(res$summary$Converged[1] %||% FALSE))

  person_labels <- character(0)
  nonperson_labels <- character(0)
  all_labels <- character(0)
  if (nrow(measure_df) > 0 && "SE_Method" %in% names(measure_df) && "Facet" %in% names(measure_df)) {
    all_labels <- unique(as.character(stats::na.omit(measure_df$SE_Method)))
    person_labels <- unique(as.character(stats::na.omit(measure_df$SE_Method[measure_df$Facet == "Person"])))
    nonperson_labels <- unique(as.character(stats::na.omit(measure_df$SE_Method[measure_df$Facet != "Person"])))
  }

  has_fallback <- identical(method, "MML") &&
    any(grepl("Fallback observation-table information", all_labels))
  precision_tier <- if (!identical(method, "MML")) {
    "exploratory"
  } else if (has_fallback) {
    "hybrid"
  } else {
    "model_based"
  }
  supports_formal <- identical(precision_tier, "model_based") && converged

  has_fit_adjusted <- nrow(measure_df) > 0 &&
    "RealSE" %in% names(measure_df) &&
    any(is.finite(suppressWarnings(as.numeric(measure_df$RealSE))))

  coverage_complete <- FALSE
  if (nrow(facet_precision_tbl) > 0 &&
      all(c("Facet", "DistributionBasis", "SEMode") %in% names(facet_precision_tbl))) {
    precision_facets <- unique(as.character(facet_precision_tbl$Facet))
    precision_facets <- precision_facets[!is.na(precision_facets) & nzchar(precision_facets)]
    if (length(precision_facets) > 0) {
      expected <- expand.grid(
        DistributionBasis = c("sample", "population"),
        SEMode = c("model", "fit_adjusted"),
        stringsAsFactors = FALSE
      )
      coverage_complete <- all(vapply(precision_facets, function(facet_name) {
        sub <- facet_precision_tbl[facet_precision_tbl$Facet == facet_name, c("DistributionBasis", "SEMode"), drop = FALSE]
        all(apply(expected, 1, function(row) {
          any(sub$DistributionBasis == row[["DistributionBasis"]] & sub$SEMode == row[["SEMode"]])
        }))
      }, logical(1)))
    }
  }

  tibble(
    Method = method,
    Converged = converged,
    PrecisionTier = precision_tier,
    SupportsFormalInference = supports_formal,
    HasFallbackSE = has_fallback,
    PersonSEBasis = if (length(person_labels) > 0) paste(sort(person_labels), collapse = "; ") else NA_character_,
    NonPersonSEBasis = if (length(nonperson_labels) > 0) paste(sort(nonperson_labels), collapse = "; ") else NA_character_,
    CIBasis = if (identical(precision_tier, "model_based")) {
      if (isTRUE(converged)) {
        "Normal interval from model-based SE"
      } else {
        "Normal interval from model-based SE; optimizer convergence review required"
      }
    } else if (identical(precision_tier, "hybrid")) {
      "Normal interval from mixed model-based and fallback SE"
    } else {
      "Normal interval from exploratory SE"
    },
    ReliabilityBasis = if (identical(precision_tier, "model_based")) {
      if (isTRUE(converged)) {
        "Observed variance with model-based and fit-adjusted error bounds"
      } else {
        "Observed variance with model-based and fit-adjusted error bounds; optimizer convergence review required"
      }
    } else if (identical(precision_tier, "hybrid")) {
      "Observed variance with mixed model-based and fallback error bounds"
    } else {
      "Exploratory variance summary with model-based and fit-adjusted error bounds"
    },
    HasFitAdjustedSE = has_fit_adjusted,
    HasSamplePopulationCoverage = coverage_complete,
    RecommendedUse = if (identical(precision_tier, "model_based") && isTRUE(converged)) {
      "Use for primary reporting of SE, CI, and reliability in this package."
    } else if (identical(precision_tier, "model_based")) {
      "This run reached the model-based precision path, but optimizer convergence should be reviewed before primary reporting."
    } else if (identical(precision_tier, "hybrid")) {
      "Use model-based rows for primary reporting, but review levels that fell back to observation-table information before treating the whole run as formal inference."
    } else {
      "Use for screening and calibration triage; confirm formal SE, CI, and reliability with an MML fit."
    }
  )
}

audit_precision_outputs <- function(res, measure_df, reliability_tbl, facet_precision_tbl, precision_profile_tbl = NULL) {
  measure_df <- as.data.frame(measure_df %||% data.frame(), stringsAsFactors = FALSE)
  reliability_tbl <- as.data.frame(reliability_tbl %||% data.frame(), stringsAsFactors = FALSE)
  facet_precision_tbl <- as.data.frame(facet_precision_tbl %||% data.frame(), stringsAsFactors = FALSE)
  precision_profile_tbl <- as.data.frame(precision_profile_tbl %||% data.frame(), stringsAsFactors = FALSE)

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  method <- ifelse(identical(method, "JMLE"), "JML", method)
  converged <- isTRUE(as.logical(res$summary$Converged[1] %||% FALSE))

  finite_model <- if ("ModelSE" %in% names(measure_df)) is.finite(suppressWarnings(as.numeric(measure_df$ModelSE))) else logical(0)
  model_share <- if (length(finite_model) > 0) mean(finite_model) else NA_real_

  real_ok <- TRUE
  if (all(c("ModelSE", "RealSE") %in% names(measure_df))) {
    model_vals <- suppressWarnings(as.numeric(measure_df$ModelSE))
    real_vals <- suppressWarnings(as.numeric(measure_df$RealSE))
    ok <- is.finite(model_vals) & is.finite(real_vals)
    if (any(ok)) {
      real_ok <- all(real_vals[ok] >= model_vals[ok])
    }
  }

  rel_ok <- TRUE
  if (all(c("Reliability", "RealReliability") %in% names(reliability_tbl))) {
    rel_vals <- suppressWarnings(as.numeric(reliability_tbl$Reliability))
    real_rel_vals <- suppressWarnings(as.numeric(reliability_tbl$RealReliability))
    ok <- is.finite(rel_vals) & is.finite(real_rel_vals)
    if (any(ok)) {
      rel_ok <- all(real_rel_vals[ok] <= rel_vals[ok])
    }
  }

  coverage_ok <- FALSE
  if (nrow(precision_profile_tbl) > 0 && "HasSamplePopulationCoverage" %in% names(precision_profile_tbl)) {
    coverage_ok <- isTRUE(precision_profile_tbl$HasSamplePopulationCoverage[1])
  }

  label_status <- "review"
  label_detail <- "SE source labels were not available for audit."
  if (nrow(measure_df) > 0 && all(c("Facet", "SE_Method") %in% names(measure_df))) {
    labels <- as.character(measure_df$SE_Method)
    person_labels <- labels[measure_df$Facet == "Person"]
    nonperson_labels <- labels[measure_df$Facet != "Person"]
    if (identical(method, "MML")) {
      has_person_eap <- length(person_labels) == 0 || any(grepl("Posterior SD \\(EAP\\)", person_labels))
      has_nonperson_info <- length(nonperson_labels) == 0 || any(grepl("Observed information \\(MML\\)", nonperson_labels))
      has_fallback <- any(grepl("Fallback observation-table information", labels))
      if (has_person_eap && has_nonperson_info && !has_fallback) {
        label_status <- "pass"
        label_detail <- "Person and non-person SE labels match the MML precision path."
      } else if (has_person_eap && (has_nonperson_info || has_fallback)) {
        label_status <- "review"
        label_detail <- "MML labels were found, but at least one level fell back to observation-table information."
      } else {
        label_status <- "warn"
        label_detail <- "MML SE labels did not consistently identify posterior-SD/person and observed-information/facet sources."
      }
    } else {
      uses_obs_info <- all(grepl("Observation-table information", labels))
      if (uses_obs_info) {
        label_status <- "pass"
        label_detail <- "JML SE labels consistently identify observation-table information."
      } else {
        label_status <- "warn"
        label_detail <- "JML SE labels were expected to indicate observation-table information throughout."
      }
    }
  }

  tibble(
    Check = c(
      "Precision tier",
      "Optimizer convergence",
      "ModelSE availability",
      "Fit-adjusted SE ordering",
      "Reliability ordering",
      "Facet precision coverage",
      "SE source labels"
    ),
    Status = c(
      if (nrow(precision_profile_tbl) > 0) {
        if (identical(as.character(precision_profile_tbl$PrecisionTier[1]), "model_based")) "pass" else "review"
      } else if (identical(method, "MML")) {
        "review"
      } else {
        "review"
      },
      if (isTRUE(converged)) "pass" else "review",
      if (!is.finite(model_share)) "warn" else if (model_share >= 0.99) "pass" else if (model_share >= 0.90) "review" else "warn",
      if (isTRUE(real_ok)) "pass" else "warn",
      if (isTRUE(rel_ok)) "pass" else "warn",
      if (isTRUE(coverage_ok)) "pass" else "review",
      label_status
    ),
    Detail = c(
      if (nrow(precision_profile_tbl) > 0 && identical(as.character(precision_profile_tbl$PrecisionTier[1]), "model_based")) {
        if (isTRUE(converged)) {
          "This run uses the package's model-based precision path."
        } else {
          "This run uses the package's model-based precision path, but optimizer convergence should be reviewed before formal reporting."
        }
      } else if (nrow(precision_profile_tbl) > 0 && identical(as.character(precision_profile_tbl$PrecisionTier[1]), "hybrid")) {
        "This run mixes model-based SE with fallback observation-table SE for at least one level; review before treating the full run as formal inference."
      } else {
        "This run uses the package's exploratory precision path; prefer MML for formal SE, CI, and reliability reporting."
      },
      if (isTRUE(converged)) {
        "The optimizer reported convergence."
      } else {
        "The optimizer did not report convergence; keep SE, CI, and reliability in review mode."
      },
      if (!is.finite(model_share)) {
        "ModelSE coverage could not be computed."
      } else {
        paste0("Finite ModelSE values were available for ", sprintf("%.1f", 100 * model_share), "% of rows.")
      },
      if (isTRUE(real_ok)) {
        "Fit-adjusted SE values were not smaller than their paired ModelSE values."
      } else {
        "At least one RealSE value was smaller than ModelSE."
      },
      if (isTRUE(rel_ok)) {
        "Conservative reliability values were not larger than the model-based values."
      } else {
        "At least one RealReliability value exceeded Reliability."
      },
      if (isTRUE(coverage_ok)) {
        "Each facet had sample/population summaries for both model and fit-adjusted SE modes."
      } else {
        "At least one facet is missing a sample/population or model/fit-adjusted precision combination."
      },
      label_detail
    )
  )
}

calc_reliability <- function(measure_df) {
  measure_df |>
    dplyr::group_by(Facet) |>
    dplyr::group_modify(function(.x, .y) {
      model_stats <- summarize_precision_basis(.x, if ("ModelSE" %in% names(.x)) "ModelSE" else "SE", distribution_basis = "sample")
      real_stats <- summarize_precision_basis(.x, if ("RealSE" %in% names(.x)) "RealSE" else if ("ModelSE" %in% names(.x)) "ModelSE" else "SE", distribution_basis = "sample")
      tier_vals <- unique(as.character(stats::na.omit(.x$PrecisionTier %||% character(0))))
      converged_vals <- unique(as.logical(stats::na.omit(.x$Converged %||% logical(0))))
      facet_converged <- if (length(converged_vals) == 0) {
        FALSE
      } else {
        all(converged_vals)
      }
      facet_tier <- if (length(tier_vals) == 0) {
        NA_character_
      } else if (all(tier_vals == "model_based")) {
        "model_based"
      } else if (any(tier_vals == "exploratory")) {
        "exploratory"
      } else {
        "hybrid"
      }
      supports_formal <- identical(facet_tier, "model_based") && isTRUE(facet_converged)

      tibble(
        Levels = nrow(.x),
        Converged = facet_converged,
        PrecisionTier = facet_tier,
        SupportsFormalInference = supports_formal,
        ReliabilityUse = dplyr::case_when(
          identical(facet_tier, "model_based") && isTRUE(facet_converged) ~ "primary_reporting",
          identical(facet_tier, "model_based") ~ "review_before_reporting",
          identical(facet_tier, "hybrid") ~ "review_before_reporting",
          identical(facet_tier, "exploratory") ~ "screening_only",
          TRUE ~ NA_character_
        ),
        SD = model_stats$ObservedSD,
        RMSE = model_stats$RMSE,
        TrueSD = model_stats$TrueSD,
        ObservedVariance = model_stats$ObservedVariance,
        ModelErrorVariance = model_stats$ErrorVariance,
        ModelTrueVariance = model_stats$TrueVariance,
        Separation = model_stats$Separation,
        Strata = model_stats$Strata,
        Reliability = model_stats$Reliability,
        ModelRMSE = model_stats$RMSE,
        ModelTrueSD = model_stats$TrueSD,
        ModelSeparation = model_stats$Separation,
        ModelStrata = model_stats$Strata,
        ModelReliability = model_stats$Reliability,
        RealRMSE = real_stats$RMSE,
        RealTrueSD = real_stats$TrueSD,
        RealErrorVariance = real_stats$ErrorVariance,
        RealTrueVariance = real_stats$TrueVariance,
        RealSeparation = real_stats$Separation,
        RealStrata = real_stats$Strata,
        RealReliability = real_stats$Reliability,
        MeanInfit = if ("Infit" %in% names(.x)) mean(.x$Infit, na.rm = TRUE) else NA_real_,
        MeanOutfit = if ("Outfit" %in% names(.x)) mean(.x$Outfit, na.rm = TRUE) else NA_real_
      )
    }) |>
    dplyr::ungroup()
}

ensure_positive_definite <- function(mat) {
  eig <- tryCatch(eigen(mat, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
  if (is.null(eig)) return(mat)
  if (any(eig < .Machine$double.eps)) {
    smoothed <- tryCatch(suppressWarnings(psych::cor.smooth(mat)), error = function(e) NULL)
    if (!is.null(smoothed)) {
      if (is.list(smoothed) && !is.null(smoothed$R)) {
        mat <- smoothed$R
      } else {
        mat <- smoothed
      }
    }
  }
  mat
}

compute_pca_overall <- function(obs_df, facet_names, max_factors = 10L) {
  pca_failure <- function(message, residual_matrix = NULL, cor_matrix = NULL) {
    list(pca = NULL, residual_matrix = residual_matrix, cor_matrix = cor_matrix, error = message)
  }
  if (length(facet_names) == 0) return(NULL)
  df_aug <- obs_df |>
    mutate(
      Person = as.character(Person),
      item_combination = paste(!!!rlang::syms(facet_names), sep = "_")
    ) |>
    select(Person, item_combination, StdResidual)

  residual_matrix_prep <- df_aug |>
    group_by(Person, item_combination) |>
    summarize(std_residual = mean(StdResidual, na.rm = TRUE), .groups = "drop")

  residual_matrix_wide <- tryCatch({
    residual_matrix_prep |>
      tidyr::pivot_wider(
        id_cols = Person,
        names_from = item_combination,
        values_from = std_residual,
        values_fill = list(std_residual = NA)
      ) |>
      tibble::column_to_rownames("Person")
  }, error = function(e) e)

  if (inherits(residual_matrix_wide, "error")) {
    return(pca_failure(paste0("Could not reshape the residual matrix: ", conditionMessage(residual_matrix_wide))))
  }
  if (is.null(residual_matrix_wide) || nrow(residual_matrix_wide) < 2 || ncol(residual_matrix_wide) < 2) {
    return(pca_failure("Residual PCA requires at least two persons and two combined-facet columns."))
  }

  residual_matrix_clean <- residual_matrix_wide[, colSums(is.na(residual_matrix_wide)) < nrow(residual_matrix_wide), drop = FALSE]
  if (ncol(residual_matrix_clean) < 2) {
    return(pca_failure("Residual PCA requires at least two combined-facet columns with observed residuals.", residual_matrix = residual_matrix_wide))
  }

  cor_matrix <- tryCatch(suppressWarnings(stats::cor(residual_matrix_clean, use = "pairwise.complete.obs")), error = function(e) e)
  if (inherits(cor_matrix, "error")) {
    return(pca_failure(
      paste0("Could not compute the residual correlation matrix: ", conditionMessage(cor_matrix)),
      residual_matrix = residual_matrix_wide
    ))
  }
  if (is.null(cor_matrix)) {
    return(pca_failure("Residual correlation matrix was not available.", residual_matrix = residual_matrix_wide))
  }
  cor_matrix[is.na(cor_matrix)] <- 0
  diag(cor_matrix) <- 1
  cor_matrix <- ensure_positive_definite(cor_matrix)

  n_factors <- max(1, min(as.integer(max_factors), ncol(cor_matrix) - 1, nrow(cor_matrix) - 1))
  pca_result <- tryCatch(psych::principal(cor_matrix, nfactors = n_factors, rotate = "none"), error = function(e) e)
  if (inherits(pca_result, "error")) {
    return(pca_failure(
      paste0("Principal components could not be computed: ", conditionMessage(pca_result)),
      residual_matrix = residual_matrix_wide,
      cor_matrix = cor_matrix
    ))
  }
  list(pca = pca_result, residual_matrix = residual_matrix_wide, cor_matrix = cor_matrix, error = NULL)
}

compute_pca_by_facet <- function(obs_df, facet_names, max_factors = 10L) {
  out <- list()
  for (facet in facet_names) {
    pca_failure <- function(message, residual_matrix = NULL, cor_matrix = NULL) {
      list(pca = NULL, residual_matrix = residual_matrix, cor_matrix = cor_matrix, error = message)
    }
    facet_sym <- rlang::sym(facet)
    prep <- obs_df |>
      mutate(Person = as.character(Person), .Level = as.character(!!facet_sym)) |>
      select(Person, .Level, StdResidual) |>
      group_by(Person, .Level) |>
      summarize(std_residual = mean(StdResidual, na.rm = TRUE), .groups = "drop")

    wide <- tryCatch({
      tidyr::pivot_wider(
        prep,
        id_cols = Person,
        names_from = .Level,
        values_from = std_residual,
        values_fill = list(std_residual = NA)
      ) |>
        tibble::column_to_rownames("Person")
    }, error = function(e) e)

    if (inherits(wide, "error")) {
      out[[facet]] <- pca_failure(paste0("Could not reshape the residual matrix: ", conditionMessage(wide)))
      next
    }
    if (is.null(wide) || nrow(wide) < 2 || ncol(wide) < 2) {
      out[[facet]] <- pca_failure("Residual PCA requires at least two persons and two facet-level columns.")
      next
    }

    keep <- colSums(is.na(wide)) < nrow(wide)
    wide <- wide[, keep, drop = FALSE]
    if (ncol(wide) < 2) {
      out[[facet]] <- pca_failure("Residual PCA requires at least two facet-level columns with observed residuals.", residual_matrix = wide)
      next
    }

    cor_mat <- tryCatch(suppressWarnings(stats::cor(wide, use = "pairwise.complete.obs")), error = function(e) e)
    if (inherits(cor_mat, "error")) {
      out[[facet]] <- pca_failure(
        paste0("Could not compute the residual correlation matrix: ", conditionMessage(cor_mat)),
        residual_matrix = wide
      )
      next
    }
    if (is.null(cor_mat)) {
      out[[facet]] <- pca_failure("Residual correlation matrix was not available.", residual_matrix = wide)
      next
    }
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1
    cor_mat <- ensure_positive_definite(cor_mat)

    n_factors <- max(1, min(as.integer(max_factors), ncol(cor_mat) - 1, nrow(cor_mat) - 1))
    pca_obj <- tryCatch(psych::principal(cor_mat, nfactors = n_factors, rotate = "none"), error = function(e) e)
    if (inherits(pca_obj, "error")) {
      out[[facet]] <- pca_failure(
        paste0("Principal components could not be computed: ", conditionMessage(pca_obj)),
        residual_matrix = wide,
        cor_matrix = cor_mat
      )
      next
    }
    out[[facet]] <- list(pca = pca_obj, cor_matrix = cor_mat, residual_matrix = wide, error = NULL)
  }
  out
}

mfrm_diagnostics <- function(res,
                             interaction_pairs = NULL,
                             top_n_interactions = 20,
                             whexact = FALSE,
                             residual_pca = c("none", "overall", "facet", "both"),
                             pca_max_factors = 10L) {
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))
  obs_df <- compute_obs_table(res)
  facet_cols <- c("Person", res$config$facet_names)
  overall_fit <- calc_overall_fit(obs_df, whexact = whexact)
  fit_tbl <- calc_facet_fit(obs_df, facet_cols, whexact = whexact)
  se_tbl <- build_measure_se_table(res, obs_df, facet_cols, fit_tbl)
  bias_tbl <- calc_bias_facet(obs_df, facet_cols)
  interaction_tbl <- calc_bias_interactions(obs_df, facet_cols, pairs = interaction_pairs,
                                             top_n = top_n_interactions)
  ptmea_tbl <- calc_ptmea(obs_df, facet_cols)
  subset_tbls <- calc_subsets(obs_df, facet_cols)

  person_tbl <- res$facets$person |>
    mutate(
      Facet = "Person",
      Level = Person
    )
  facet_tbl <- res$facets$others |>
    mutate(Level = as.character(Level))

  measures <- bind_rows(
    person_tbl |>
      select(Facet, Level, Estimate),
    facet_tbl |>
      select(Facet, Level, Estimate)
  ) |>
    left_join(se_tbl, by = c("Facet", "Level")) |>
    left_join(fit_tbl, by = c("Facet", "Level")) |>
    left_join(bias_tbl, by = c("Facet", "Level")) |>
    left_join(ptmea_tbl, by = c("Facet", "Level")) |>
    mutate(
      CI_Lower = ifelse(is.finite(SE), Estimate - 1.96 * SE, NA_real_),
      CI_Upper = ifelse(is.finite(SE), Estimate + 1.96 * SE, NA_real_),
      CIEligible = dplyr::coalesce(.data$SupportsFormalInference, FALSE),
      CILabel = dplyr::case_when(
        .data$PrecisionTier == "model_based" ~ "Model-based normal interval",
        .data$PrecisionTier == "hybrid" ~ "Approximate interval; review fallback SE",
        .data$PrecisionTier == "exploratory" ~ "Approximate interval; screening only",
        TRUE ~ "Approximate interval"
      )
    )

  if (!"N" %in% names(measures)) {
    if ("N.x" %in% names(measures) || "N.y" %in% names(measures)) {
      measures <- measures |>
        mutate(N = dplyr::coalesce(.data$N.x, .data$N.y))
    } else {
      measures <- measures |>
        mutate(N = NA_real_)
    }
  }

  reliability_tbl <- calc_reliability(measures)
  facets_chisq_tbl <- calc_facets_chisq(measures)
  facet_precision_tbl <- build_facet_precision_summary(measures, facets_chisq_tbl)
  default_rater_facet <- infer_default_rater_facet(res$config$facet_names)
  interrater_tbl <- augment_interrater_with_precision(
    calc_interrater_agreement(
      obs_df = obs_df,
      facet_cols = facet_cols,
      rater_facet = default_rater_facet,
      res = res
    ),
    reliability_tbl = reliability_tbl,
    rater_facet = default_rater_facet
  )

  unexpected_abs_z_min <- 2
  unexpected_prob_max <- 0.30
  unexpected_rule <- "either"
  unexpected_tbl <- calc_unexpected_response_table(
    obs_df = obs_df,
    probs = compute_prob_matrix(res),
    facet_names = res$config$facet_names,
    rating_min = res$prep$rating_min,
    abs_z_min = unexpected_abs_z_min,
    prob_max = unexpected_prob_max,
    top_n = 100,
    rule = unexpected_rule
  )
  unexpected_summary <- summarize_unexpected_response_table(
    unexpected_tbl = unexpected_tbl,
    total_observations = nrow(obs_df),
    abs_z_min = unexpected_abs_z_min,
    prob_max = unexpected_prob_max,
    rule = unexpected_rule
  )

  fair_average <- calc_fair_average_bundle(
    res = res,
    diagnostics = list(obs = obs_df, measures = measures),
    totalscore = TRUE,
    umean = 0,
    uscale = 1,
    udecimals = 2,
    omit_unobserved = FALSE,
    xtreme = 0
  )

  displacement_abs_warn <- 0.5
  displacement_abs_t_warn <- 2
  displacement_tbl <- calc_displacement_table(
    obs_df = obs_df,
    res = res,
    measures = measures,
    abs_displacement_warn = displacement_abs_warn,
    abs_t_warn = displacement_abs_t_warn
  )
  displacement_summary <- summarize_displacement_table(
    displacement_tbl = displacement_tbl,
    abs_displacement_warn = displacement_abs_warn,
    abs_t_warn = displacement_abs_t_warn
  )

  pca_overall <- NULL
  pca_by_facet <- NULL
  if (residual_pca %in% c("overall", "both")) {
    pca_overall <- compute_pca_overall(
      obs_df = obs_df,
      facet_names = res$config$facet_names,
      max_factors = pca_max_factors
    )
  }
  if (residual_pca %in% c("facet", "both")) {
    pca_by_facet <- compute_pca_by_facet(
      obs_df = obs_df,
      facet_names = res$config$facet_names,
      max_factors = pca_max_factors
    )
  }

  method <- as.character(res$summary$Method[1] %||% res$config$method %||% NA_character_)
  se_detail <- attr(se_tbl, "mml_se_detail", exact = TRUE) %||%
    "SE values default to observation-table information approximations."
  precision_profile_tbl <- build_precision_profile(res, measures, reliability_tbl, facet_precision_tbl)
  precision_audit_tbl <- audit_precision_outputs(
    res,
    measures,
    reliability_tbl,
    facet_precision_tbl,
    precision_profile_tbl = precision_profile_tbl
  )
  approximation_notes <- tibble::tibble(
    Component = c("SE", "RealSE", "CI", "Reliability"),
    Method = c(
      if (identical(method, "MML")) "Model-based SE (MML)" else "Model-based SE (exploratory JML)",
      "Fit-adjusted SE",
      "Normal interval",
      "Variance decomposition"
    ),
    Detail = c(
      paste0(
        "The `SE` column is kept as a compatibility alias for `ModelSE`. ",
        if (identical(method, "MML")) {
          "Persons use posterior SD from EAP estimation; non-person facets use observed-information SEs from the marginal log-likelihood. "
        } else {
          "JML uses observation-table information as an exploratory approximation. "
        },
        se_detail,
        " Row-level `PrecisionTier`, `Converged`, `SupportsFormalInference`, and `SEUse` indicate whether a given estimate is suitable for primary reporting."
      ),
      "RealSE is a fit-adjusted companion to ModelSE: ModelSE * sqrt(max(Infit, 1)), so it is never smaller than ModelSE.",
      "CI_Lower and CI_Upper are symmetric bands computed as Estimate +/- 1.96 * SE. Use `CIEligible`, `Converged`, `CIBasis`, and `CIUse` to distinguish primary-reporting intervals from review or screening approximations.",
      "Reliability tables report model and real bounds using observed variance, error variance, and true variance (Observed variance - mean SE^2). `Reliability`/`Separation` remain compatibility aliases for the model-based values, while `PrecisionTier`, `Converged`, `SupportsFormalInference`, and `ReliabilityUse` indicate how strongly each facet summary supports formal reporting."
    )
  )

  list(
    obs = obs_df,
    facet_names = res$config$facet_names,
    overall_fit = overall_fit,
    measures = measures,
    fit = fit_tbl,
    reliability = reliability_tbl,
    precision_profile = precision_profile_tbl,
    precision_audit = precision_audit_tbl,
    facet_precision = facet_precision_tbl,
    facets_chisq = facets_chisq_tbl,
    bias = bias_tbl,
    interactions = interaction_tbl,
    interrater = interrater_tbl,
    unexpected = list(
      table = unexpected_tbl,
      summary = unexpected_summary,
      thresholds = list(
        abs_z_min = unexpected_abs_z_min,
        prob_max = unexpected_prob_max,
        rule = unexpected_rule
      )
    ),
    fair_average = fair_average,
    displacement = list(
      table = displacement_tbl,
      summary = displacement_summary,
      thresholds = list(
        abs_displacement_warn = displacement_abs_warn,
        abs_t_warn = displacement_abs_t_warn
      )
    ),
    approximation_notes = approximation_notes,
    subsets = subset_tbls,
    residual_pca_mode = residual_pca,
    residual_pca_overall = pca_overall,
    residual_pca_by_facet = pca_by_facet
  )
}
