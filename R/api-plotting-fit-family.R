# Plot methods centered on fitted many-facet Rasch objects.

step_index_from_label <- function(step_labels) {
  step_labels <- as.character(step_labels)
  idx <- suppressWarnings(as.integer(gsub("[^0-9]+", "", step_labels)))
  if (all(is.na(idx))) return(seq_along(step_labels))
  if (any(is.na(idx))) {
    fill_start <- max(idx, na.rm = TRUE)
    idx[is.na(idx)] <- fill_start + seq_len(sum(is.na(idx)))
  }
  idx
}

build_step_curve_spec <- function(x) {
  step_tbl <- x$steps
  if (is.null(step_tbl) || nrow(step_tbl) == 0 || !"Estimate" %in% names(step_tbl)) {
    stop("Step estimates are required for pathway/CCC plots.")
  }
  step_tbl <- tibble::as_tibble(step_tbl)

  model <- toupper(as.character(x$config$model[1]))
  rating_min <- suppressWarnings(min(as.numeric(x$prep$data$Score), na.rm = TRUE))
  if (!is.finite(rating_min)) rating_min <- 0

  n_cat <- suppressWarnings(as.integer(x$config$n_cat[1]))
  if (!is.finite(n_cat) || n_cat < 2) {
    n_cat <- max(2L, nrow(step_tbl) + 1L)
  }
  categories <- rating_min + 0:(n_cat - 1L)

  groups <- list()
  step_points <- tibble::tibble()

  if (model == "RSM") {
    if (!"Step" %in% names(step_tbl)) {
      step_tbl$Step <- paste0("Step_", seq_len(nrow(step_tbl)))
    }
    ord <- order(step_index_from_label(step_tbl$Step))
    tau <- as.numeric(step_tbl$Estimate[ord])
    groups[["Common"]] <- list(name = "Common", step_cum = c(0, cumsum(tau)), tau = tau)
    step_points <- tibble::tibble(
      CurveGroup = "Common",
      Step = as.character(step_tbl$Step[ord]),
      StepIndex = seq_along(tau),
      Threshold = tau
    )
  } else {
    if (!all(c("StepFacet", "Step", "Estimate") %in% names(step_tbl))) {
      stop("PCM step table must include StepFacet, Step, and Estimate columns.")
    }
    step_facet <- as.character(x$config$step_facet[1])
    ordered_levels <- x$prep$levels[[step_facet]]
    if (is.null(ordered_levels)) {
      ordered_levels <- unique(as.character(step_tbl$StepFacet))
    }
    for (lvl in ordered_levels) {
      sub <- step_tbl[as.character(step_tbl$StepFacet) == as.character(lvl), , drop = FALSE]
      if (nrow(sub) == 0) next
      ord <- order(step_index_from_label(sub$Step))
      tau <- as.numeric(sub$Estimate[ord])
      groups[[as.character(lvl)]] <- list(
        name = as.character(lvl),
        step_cum = c(0, cumsum(tau)),
        tau = tau
      )
      step_points <- dplyr::bind_rows(
        step_points,
        tibble::tibble(
          CurveGroup = as.character(lvl),
          Step = as.character(sub$Step[ord]),
          StepIndex = seq_along(tau),
          Threshold = tau
        )
      )
    }
  }

  if (length(groups) == 0) stop("No step groups available for pathway/CCC plots.")

  list(
    model = model,
    categories = categories,
    rating_min = rating_min,
    n_cat = n_cat,
    groups = groups,
    step_points = step_points
  )
}

build_curve_tables <- function(curve_spec, theta_grid) {
  prob_tables <- list()
  exp_tables <- list()
  idx_prob <- 1L
  idx_exp <- 1L
  for (g in names(curve_spec$groups)) {
    grp <- curve_spec$groups[[g]]
    probs <- category_prob_rsm(theta_grid, grp$step_cum)
    k_vals <- as.numeric(curve_spec$categories)
    expected <- as.numeric(probs %*% matrix(k_vals, ncol = 1))
    exp_tables[[idx_exp]] <- tibble::tibble(
      Theta = theta_grid,
      ExpectedScore = expected,
      CurveGroup = grp$name
    )
    idx_exp <- idx_exp + 1L
    for (k in seq_len(ncol(probs))) {
      prob_tables[[idx_prob]] <- tibble::tibble(
        Theta = theta_grid,
        Probability = probs[, k],
        Category = as.character(curve_spec$categories[k]),
        CurveGroup = grp$name
      )
      idx_prob <- idx_prob + 1L
    }
  }
  list(
    probabilities = dplyr::bind_rows(prob_tables),
    expected = dplyr::bind_rows(exp_tables)
  )
}

compute_se_for_plot <- function(x, ci_level = 0.95) {
  tryCatch({
    obs_df <- compute_obs_table(x)
    facet_cols <- x$config$source_columns$facets
    if (is.null(facet_cols) || length(facet_cols) == 0) {
      facet_cols <- x$config$facet_names
    }
    se_tbl <- calc_facet_se(obs_df, facet_cols)
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    se_tbl$CI_Lower <- -z * se_tbl$SE
    se_tbl$CI_Upper <- z * se_tbl$SE
    se_tbl$CI_Level <- ci_level
    as.data.frame(se_tbl, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

build_wright_map_data <- function(x, top_n = 30L, se_tbl = NULL, include_steps = TRUE) {
  person_tbl <- tibble::as_tibble(x$facets$person)
  person_tbl <- person_tbl[is.finite(person_tbl$Estimate), , drop = FALSE]
  if (nrow(person_tbl) == 0) stop("Person estimates are not available for Wright map.")

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (!all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available for Wright map.")
  }
  facet_tbl <- facet_tbl[is.finite(facet_tbl$Estimate), , drop = FALSE]

  step_points <- if (isTRUE(include_steps)) {
    tryCatch(build_step_curve_spec(x)$step_points, error = function(e) tibble::tibble())
  } else {
    tibble::tibble()
  }
  step_tbl <- if (nrow(step_points) > 0) {
    tibble::tibble(
      PlotType = "Step threshold",
      Group = if ("CurveGroup" %in% names(step_points)) paste0("Step:", step_points$CurveGroup) else "Step",
      Label = step_points$Step,
      Estimate = step_points$Threshold
    )
  } else {
    tibble::tibble()
  }

  facet_pts <- facet_tbl |>
    dplyr::transmute(
      PlotType = "Facet level",
      Group = .data$Facet,
      Label = .data$Level,
      Estimate = .data$Estimate
    )
  if (!is.null(se_tbl) && is.data.frame(se_tbl) && nrow(se_tbl) > 0 &&
      all(c("Facet", "Level", "SE") %in% names(se_tbl))) {
    se_join <- se_tbl[, intersect(c("Facet", "Level", "SE", "CI_Lower", "CI_Upper"), names(se_tbl)), drop = FALSE]
    se_join$Level <- as.character(se_join$Level)
    facet_pts$Group <- as.character(facet_pts$Group)
    facet_pts$Label <- as.character(facet_pts$Label)
    facet_pts <- merge(
      facet_pts,
      se_join,
      by.x = c("Group", "Label"),
      by.y = c("Facet", "Level"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  point_tbl <- dplyr::bind_rows(facet_pts, step_tbl)
  if (nrow(point_tbl) == 0) stop("No facet/step locations available for Wright map.")

  if (nrow(point_tbl) > top_n) {
    point_tbl <- point_tbl |>
      dplyr::mutate(.Abs = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$.Abs)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(-.data$.Abs)
  }

  group_levels <- unique(point_tbl$Group)
  point_tbl <- point_tbl |>
    dplyr::mutate(Group = factor(.data$Group, levels = group_levels)) |>
    dplyr::group_by(.data$Group) |>
    dplyr::arrange(.data$Estimate, .by_group = TRUE) |>
    dplyr::mutate(
      XBase = as.numeric(.data$Group),
      X = if (dplyr::n() == 1) .data$XBase else .data$XBase + seq(-0.18, 0.18, length.out = dplyr::n())
    ) |>
    dplyr::ungroup()

  hist_data <- graphics::hist(person_tbl$Estimate, breaks = "FD", plot = FALSE)
  y_range <- range(c(point_tbl$Estimate, person_tbl$Estimate), finite = TRUE)
  if (!all(is.finite(y_range)) || diff(y_range) <= sqrt(.Machine$double.eps)) {
    center <- mean(c(point_tbl$Estimate, person_tbl$Estimate), na.rm = TRUE)
    y_range <- center + c(-0.5, 0.5)
  }

  label_tbl <- point_tbl |>
    dplyr::group_by(.data$Group) |>
    dplyr::slice(c(1L, dplyr::n())) |>
    dplyr::ungroup() |>
    dplyr::distinct(.data$Group, .data$Label, .keep_all = TRUE)

  group_summary <- point_tbl |>
    dplyr::group_by(.data$Group, .data$PlotType) |>
    dplyr::summarise(
      Min = min(.data$Estimate, na.rm = TRUE),
      Q1 = stats::quantile(.data$Estimate, 0.25, na.rm = TRUE, names = FALSE),
      Median = stats::median(.data$Estimate, na.rm = TRUE),
      Q3 = stats::quantile(.data$Estimate, 0.75, na.rm = TRUE, names = FALSE),
      Max = max(.data$Estimate, na.rm = TRUE),
      N = dplyr::n(),
      XBase = mean(.data$XBase, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(TargetGap = .data$Median - mean(person_tbl$Estimate, na.rm = TRUE))

  person_stats <- tibble::tibble(
    N = nrow(person_tbl),
    Mean = mean(person_tbl$Estimate, na.rm = TRUE),
    Median = stats::median(person_tbl$Estimate, na.rm = TRUE),
    SD = stats::sd(person_tbl$Estimate, na.rm = TRUE)
  )

  list(
    title = "Wright Map",
    person = person_tbl,
    person_hist = hist_data,
    person_stats = person_stats,
    locations = point_tbl,
    label_points = label_tbl,
    group_summary = group_summary,
    group_levels = group_levels,
    y_range = y_range
  )
}

draw_wright_map <- function(plot_data,
                            title = NULL,
                            palette = NULL,
                            label_angle = 45,
                            show_ci = FALSE,
                            ci_level = 0.95) {
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet_level = "#1b9e77",
      step_threshold = "#d95f02",
      person_hist = "gray80",
      grid = "#ececec",
      range = "#94a3b8",
      iqr = "#334155"
    )
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.15, 2.55))
  graphics::par(mgp = c(2.6, 0.85, 0))

  loc <- plot_data$locations
  yr <- plot_data$y_range
  pretty_y <- pretty(yr, n = 6)
  hist_data <- plot_data$person_hist
  hist_max <- max(hist_data$counts, na.rm = TRUE)
  if (!is.finite(hist_max) || hist_max <= 0) hist_max <- 1

  graphics::par(mar = c(8.8, 4.6, 3.2, 0.8))
  graphics::plot(
    x = c(0, hist_max * 1.12),
    y = yr,
    type = "n",
    xlab = "Persons",
    ylab = "Logit scale",
    main = "Person distribution",
    yaxt = "n"
  )
  graphics::axis(2, at = pretty_y, las = 1)
  graphics::abline(h = pretty_y, col = pal["grid"], lty = 1)
  for (i in seq_along(hist_data$counts)) {
    graphics::rect(
      xleft = 0,
      ybottom = hist_data$breaks[i],
      xright = hist_data$counts[i],
      ytop = hist_data$breaks[i + 1],
      col = grDevices::adjustcolor(pal["person_hist"], alpha.f = 0.85),
      border = "white"
    )
  }
  if (nrow(plot_data$person_stats) > 0) {
    graphics::abline(h = plot_data$person_stats$Median[1], col = grDevices::adjustcolor("gray25", alpha.f = 0.9), lty = 2, lwd = 1.2)
    graphics::abline(h = plot_data$person_stats$Mean[1], col = grDevices::adjustcolor("gray45", alpha.f = 0.9), lty = 3, lwd = 1.2)
    graphics::mtext(
      sprintf(
        "N=%d  Mean=%.2f  Median=%.2f",
        plot_data$person_stats$N[1],
        plot_data$person_stats$Mean[1],
        plot_data$person_stats$Median[1]
      ),
      side = 3,
      line = 0.2,
      cex = 0.78,
      adj = 0
    )
  }

  xr <- c(0.5, length(plot_data$group_levels) + 0.75)
  graphics::par(mar = c(8.8, 2.4, 3.2, 2.8))
  graphics::plot(
    x = loc$X,
    y = loc$Estimate,
    type = "n",
    xlim = xr,
    ylim = yr,
    xaxt = "n",
    xlab = "",
    yaxt = "n",
    ylab = "",
    main = if (is.null(title)) "Facet and step locations" else as.character(title[1])
  )
  for (i in seq_along(plot_data$group_levels)) {
    if (i %% 2 == 1) {
      graphics::rect(
        xleft = i - 0.5,
        ybottom = yr[1] - 10,
        xright = i + 0.5,
        ytop = yr[2] + 10,
        col = grDevices::adjustcolor(pal["grid"], alpha.f = 0.18),
        border = NA
      )
    }
  }
  graphics::abline(h = pretty_y, col = pal["grid"], lty = 1)
  graphics::abline(h = 0, col = grDevices::adjustcolor("gray35", alpha.f = 0.9), lty = 2)
  group_summary <- plot_data$group_summary
  if (!is.null(group_summary) && nrow(group_summary) > 0) {
    range_cols <- ifelse(
      group_summary$PlotType == "Step threshold",
      grDevices::adjustcolor(pal["step_threshold"], alpha.f = 0.42),
      grDevices::adjustcolor(pal["range"], alpha.f = 0.72)
    )
    iqr_cols <- ifelse(
      group_summary$PlotType == "Step threshold",
      grDevices::adjustcolor(pal["step_threshold"], alpha.f = 0.78),
      grDevices::adjustcolor(pal["iqr"], alpha.f = 0.92)
    )
    graphics::segments(
      x0 = group_summary$XBase,
      y0 = group_summary$Min,
      x1 = group_summary$XBase,
      y1 = group_summary$Max,
      col = range_cols,
      lwd = 1.6
    )
    graphics::segments(
      x0 = group_summary$XBase,
      y0 = group_summary$Q1,
      x1 = group_summary$XBase,
      y1 = group_summary$Q3,
      col = iqr_cols,
      lwd = 5.2
    )
    graphics::segments(
      x0 = group_summary$XBase - 0.13,
      y0 = group_summary$Median,
      x1 = group_summary$XBase + 0.13,
      y1 = group_summary$Median,
      col = grDevices::adjustcolor("gray15", alpha.f = 0.95),
      lwd = 1.4
    )
  }
  graphics::axis(4, at = pretty_y, labels = pretty_y, las = 1)
  x_at <- seq_along(plot_data$group_levels)
  x_lab <- truncate_axis_label(plot_data$group_levels, width = 16L)
  draw_rotated_x_labels(
    at = x_at,
    labels = x_lab,
    srt = label_angle,
    cex = 0.82,
    line_offset = 0.085
  )
  cols <- ifelse(loc$PlotType == "Step threshold", pal["step_threshold"], pal["facet_level"])
  pch <- ifelse(loc$PlotType == "Step threshold", 17, 16)
  if (isTRUE(show_ci) && "SE" %in% names(loc)) {
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    ci_ok <- is.finite(loc$SE) & loc$SE > 0
    if (any(ci_ok)) {
      ci_lo <- loc$Estimate[ci_ok] - z * loc$SE[ci_ok]
      ci_hi <- loc$Estimate[ci_ok] + z * loc$SE[ci_ok]
      graphics::arrows(
        x0 = loc$X[ci_ok],
        y0 = ci_lo,
        x1 = loc$X[ci_ok],
        y1 = ci_hi,
        angle = 90,
        code = 3,
        length = 0.04,
        col = grDevices::adjustcolor(cols[ci_ok], alpha.f = 0.6),
        lwd = 1.2
      )
    }
  }
  graphics::points(
    loc$X,
    loc$Estimate,
    pch = pch,
    col = cols,
    bg = grDevices::adjustcolor(cols, alpha.f = 0.22),
    cex = ifelse(loc$PlotType == "Step threshold", 1.08, 1)
  )
  label_tbl <- plot_data$label_points
  if (!is.null(label_tbl) && nrow(label_tbl) > 0) {
    graphics::text(
      x = pmin(label_tbl$X + 0.06, xr[2] - 0.02),
      y = label_tbl$Estimate,
      labels = truncate_axis_label(label_tbl$Label, width = 14L),
      pos = 4,
      cex = 0.72,
      xpd = NA,
      col = grDevices::adjustcolor("gray20", alpha.f = 0.95)
    )
  }
  graphics::legend(
    "topleft",
    legend = c("Facet level", "Step threshold", "IQR / range"),
    pch = c(16, 17, NA),
    lty = c(0, 0, 1),
    lwd = c(NA, NA, 4),
    col = c(pal["facet_level"], pal["step_threshold"], pal["iqr"]),
    bty = "n",
    cex = 0.9
  )
}

build_pathway_map_data <- function(x, theta_range = c(-6, 6), theta_points = 241L) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  cats <- as.character(curve_spec$categories)
  step_df <- curve_spec$step_points |>
    dplyr::mutate(
      PathY = curve_spec$rating_min + .data$StepIndex - 0.5,
      ThresholdLabel = ifelse(
        .data$StepIndex < length(cats),
        paste0(cats[pmax(1L, .data$StepIndex)], "/", cats[pmin(length(cats), .data$StepIndex + 1L)]),
        as.character(.data$Step)
      )
    )
  endpoint_labels <- curve_tbl$expected |>
    dplyr::group_by(.data$CurveGroup) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
  dominance_regions <- curve_tbl$probabilities |>
    dplyr::group_by(.data$CurveGroup, .data$Theta) |>
    dplyr::slice_max(.data$Probability, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$CurveGroup, .data$Theta) |>
    dplyr::group_by(.data$CurveGroup) |>
    dplyr::mutate(
      Region = cumsum(dplyr::coalesce(.data$Category != dplyr::lag(.data$Category), TRUE))
    ) |>
    dplyr::group_by(.data$CurveGroup, .data$Category, .data$Region) |>
    dplyr::summarise(
      ThetaStart = min(.data$Theta, na.rm = TRUE),
      ThetaEnd = max(.data$Theta, na.rm = TRUE),
      ThetaMid = mean(range(.data$Theta, na.rm = TRUE)),
      .groups = "drop"
    )
  list(
    title = "Pathway Map (Expected Score by Theta)",
    expected = curve_tbl$expected,
    steps = step_df,
    endpoint_labels = endpoint_labels,
    dominance_regions = dominance_regions,
    score_range = range(curve_tbl$expected$ExpectedScore, finite = TRUE)
  )
}

draw_pathway_map <- function(plot_data, title = NULL, palette = NULL) {
  exp_df <- plot_data$expected
  groups <- unique(exp_df$CurveGroup)
  defaults <- stats::setNames(
    grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)],
    groups
  )
  cols <- resolve_palette(palette = palette, defaults = defaults)
  y_rng <- range(exp_df$ExpectedScore, finite = TRUE)
  y_breaks <- pretty(y_rng, n = 6)
  dominance_regions <- plot_data$dominance_regions
  use_strips <- !is.null(dominance_regions) && nrow(dominance_regions) > 0 && length(groups) <= 6
  band_pad <- if (isTRUE(use_strips)) diff(y_rng) * (0.08 + 0.05 * length(groups)) else 0
  graphics::plot(
    x = range(exp_df$Theta, finite = TRUE),
    y = c(y_rng[1] - band_pad, y_rng[2]),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Expected score",
    main = if (is.null(title)) plot_data$title else as.character(title[1])
  )
  graphics::abline(h = y_breaks, col = grDevices::adjustcolor("#d9dde3", alpha.f = 0.8), lty = 1)
  if (isTRUE(use_strips)) {
    band_h <- max(diff(y_rng) * 0.035, 0.08)
    band_gap <- band_h * 0.35
    band_top <- y_rng[1] - band_h * 0.45
    group_index <- stats::setNames(seq_along(groups), groups)
    for (i in seq_len(nrow(dominance_regions))) {
      grp <- as.character(dominance_regions$CurveGroup[i])
      row <- group_index[[grp]]
      y_top <- band_top - (row - 1) * (band_h + band_gap)
      y_bottom <- y_top - band_h
      graphics::rect(
        xleft = dominance_regions$ThetaStart[i],
        ybottom = y_bottom,
        xright = dominance_regions$ThetaEnd[i],
        ytop = y_top,
        col = grDevices::adjustcolor(cols[grp], alpha.f = 0.14),
        border = "white"
      )
      graphics::text(
        x = dominance_regions$ThetaMid[i],
        y = (y_bottom + y_top) / 2,
        labels = as.character(dominance_regions$Category[i]),
        cex = 0.72,
        col = cols[grp]
      )
    }
    graphics::text(
      x = rep(min(exp_df$Theta, na.rm = TRUE), length(groups)),
      y = band_top - (seq_along(groups) - 1) * (band_h + band_gap) - band_h / 2,
      labels = truncate_axis_label(groups, width = 14L),
      pos = 2,
      xpd = NA,
      cex = 0.75,
      col = cols[groups]
    )
    graphics::mtext("Dominant category by theta", side = 1, line = 4.6, adj = 0)
  }
  if (!is.null(plot_data$steps) && nrow(plot_data$steps) > 0) {
    step_groups <- as.character(plot_data$steps$CurveGroup)
    step_cols <- cols[step_groups]
    for (i in seq_len(nrow(plot_data$steps))) {
      graphics::abline(
        v = plot_data$steps$Threshold[i],
        col = grDevices::adjustcolor(step_cols[i], alpha.f = 0.16),
        lty = 3
      )
    }
  }
  for (i in seq_along(groups)) {
    sub <- exp_df[exp_df$CurveGroup == groups[i], , drop = FALSE]
    graphics::lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2)
  }
  if (nrow(plot_data$steps) > 0) {
    step_groups <- as.character(plot_data$steps$CurveGroup)
    step_cols <- cols[step_groups]
    graphics::points(plot_data$steps$Threshold, plot_data$steps$PathY, pch = 18, col = step_cols)
    graphics::text(
      x = plot_data$steps$Threshold,
      y = plot_data$steps$PathY,
      labels = truncate_axis_label(plot_data$steps$ThresholdLabel, width = 6L),
      pos = 3,
      cex = 0.62,
      col = step_cols,
      offset = 0.45
    )
  }
  if (!is.null(plot_data$endpoint_labels) && nrow(plot_data$endpoint_labels) > 0 && length(groups) <= 5) {
    graphics::text(
      x = plot_data$endpoint_labels$Theta,
      y = plot_data$endpoint_labels$ExpectedScore,
      labels = truncate_axis_label(plot_data$endpoint_labels$CurveGroup, width = 14L),
      pos = 4,
      cex = 0.78,
      xpd = NA,
      col = cols[as.character(plot_data$endpoint_labels$CurveGroup)]
    )
  } else {
    graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n")
  }
}

build_ccc_data <- function(x, theta_range = c(-6, 6), theta_points = 241L) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  prob_df <- build_curve_tables(curve_spec, theta_grid)$probabilities
  list(
    title = "Category Characteristic Curves",
    probabilities = prob_df
  )
}

draw_ccc <- function(plot_data, title = NULL, palette = NULL) {
  prob_df <- plot_data$probabilities
  traces <- unique(paste(prob_df$CurveGroup, prob_df$Category, sep = " | Cat "))
  defaults <- stats::setNames(
    grDevices::hcl.colors(max(3L, length(traces)), "Dark 3")[seq_along(traces)],
    traces
  )
  cols <- resolve_palette(palette = palette, defaults = defaults)
  graphics::plot(
    x = range(prob_df$Theta, finite = TRUE),
    y = c(0, 1),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Probability",
    main = if (is.null(title)) plot_data$title else as.character(title[1])
  )
  for (i in seq_along(traces)) {
    parts <- strsplit(traces[i], " \\| Cat ", fixed = FALSE)[[1]]
    sub <- prob_df[prob_df$CurveGroup == parts[1] & prob_df$Category == parts[2], , drop = FALSE]
    graphics::lines(sub$Theta, sub$Probability, col = cols[traces[i]], lwd = 1.5)
  }
}

draw_person_plot <- function(person_tbl, bins, title = "Person measure distribution", palette = NULL) {
  pal <- resolve_palette(palette = palette, defaults = c(person_hist = "#2C7FB8"))
  graphics::hist(
    x = person_tbl$Estimate,
    breaks = bins,
    main = title,
    xlab = "Person estimate",
    ylab = "Count",
    col = pal["person_hist"],
    border = "white"
  )
}

draw_step_plot <- function(step_tbl,
                           title = "Step parameter estimates",
                           palette = NULL,
                           label_angle = 45) {
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      step_line = "#1B9E77",
      grid = "#ececec"
    )
  )
  x_idx <- step_index_from_label(step_tbl$Step)
  ord <- order(x_idx)
  step_tbl <- step_tbl[ord, , drop = FALSE]
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar <- old_mar
  mar[1] <- max(mar[1], 8.6)
  graphics::par(mar = mar)
  graphics::plot(
    x = seq_len(nrow(step_tbl)),
    y = step_tbl$Estimate,
    type = "b",
    pch = 16,
    xaxt = "n",
    xlab = "",
    ylab = "Estimate",
    main = title,
    col = pal["step_line"]
  )
  graphics::abline(h = pretty(step_tbl$Estimate, n = 5), col = pal["grid"], lty = 1)
  draw_rotated_x_labels(
    at = seq_len(nrow(step_tbl)),
    labels = truncate_axis_label(step_tbl$Step, width = 16L),
    srt = label_angle,
    cex = 0.84,
    line_offset = 0.085
  )
}

draw_facet_plot <- function(facet_tbl,
                            title,
                            palette = NULL,
                            show_ci = FALSE,
                            ci_level = 0.95) {
  pal <- resolve_palette(palette = palette, defaults = c(facet_bar = "gray70"))
  labels <- paste0(facet_tbl$Facet, ": ", facet_tbl$Level)
  mids <- graphics::barplot(
    height = facet_tbl$Estimate,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    main = title,
    xlab = "Estimate",
    col = pal["facet_bar"],
    border = "white"
  )
  graphics::abline(v = 0, lty = 2, col = "gray40")
  if (isTRUE(show_ci) && "SE" %in% names(facet_tbl)) {
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    ci_ok <- is.finite(facet_tbl$SE) & facet_tbl$SE > 0
    if (any(ci_ok)) {
      ci_lo <- facet_tbl$Estimate[ci_ok] - z * facet_tbl$SE[ci_ok]
      ci_hi <- facet_tbl$Estimate[ci_ok] + z * facet_tbl$SE[ci_ok]
      graphics::arrows(
        x0 = ci_lo,
        y0 = mids[ci_ok],
        x1 = ci_hi,
        y1 = mids[ci_ok],
        angle = 90,
        code = 3,
        length = 0.04,
        col = "gray30",
        lwd = 1.2
      )
    }
  }
}

#' Plot fitted MFRM results with base R
#'
#' @param x An `mfrm_fit` object from [fit_mfrm()].
#' @param type Plot type. Use `NULL`, `"bundle"`, or `"all"` for the
#'   three-part fit bundle; otherwise choose one of `"facet"`, `"person"`,
#'   `"step"`, `"wright"`, `"pathway"`, or `"ccc"`.
#' @param facet Optional facet name for `type = "facet"`.
#' @param top_n Maximum number of facet/step locations retained for
#'   compact displays.
#' @param theta_range Numeric length-2 range for pathway and CCC plots.
#' @param theta_points Number of theta grid points used for pathway and
#'   CCC plots.
#' @param title Optional custom title.
#' @param palette Optional color overrides.
#' @param label_angle Rotation angle for x-axis labels where applicable.
#' @param show_ci If `TRUE`, add approximate confidence intervals when
#'   available.
#' @param ci_level Confidence level used when `show_ci = TRUE`.
#' @param draw If `TRUE`, draw the plot with base graphics.
#' @param preset Visual preset (`"standard"`, `"publication"`, or
#'   `"compact"`).
#' @param ... Additional arguments ignored for S3 compatibility.
#'
#' @details
#' This S3 plotting method provides the core fit-family visuals for
#' `mfrmr`. When `type` is omitted, it returns a bundle containing a
#' Wright map, pathway map, and category characteristic curves. The
#' returned object still carries machine-readable metadata through the
#' `mfrm_plot_data` contract, even when the plot is drawn immediately.
#'
#' `type = "wright"` shows persons, facet levels, and step thresholds on
#' a shared logit scale. `type = "pathway"` shows expected score traces
#' and dominant-category regions across theta. `type = "ccc"` shows
#' category response probabilities. The remaining types provide compact
#' person, step, or facet-specific displays.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Use `plot(fit)` to inspect the three core fit-family visuals.
#' 3. Switch to `type = "wright"` or `type = "pathway"` when you need a
#'    single figure for reporting or manuscript preparation.
#'
#' @section Further guidance:
#' For a plot-selection guide and extended examples, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @return Invisibly, an `mfrm_plot_data` object or an `mfrm_plot_bundle`
#'   when `type` is omitted.
#' @seealso [fit_mfrm()], [plot_wright_unified()], [plot_bubble()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy,
#'   "Person",
#'   c("Rater", "Criterion"),
#'   "Score",
#'   method = "JML",
#'   model = "RSM",
#'   maxit = 25
#' )
#' bundle <- plot(fit, draw = FALSE)
#' names(bundle)
#' if (interactive()) {
#'   plot(
#'     fit,
#'     type = "wright",
#'     preset = "publication",
#'     title = "Customized Wright Map",
#'     show_ci = TRUE,
#'     label_angle = 45
#'   )
#'   plot(
#'     fit,
#'     type = "pathway",
#'     title = "Customized Pathway Map",
#'     palette = c("#1f78b4")
#'   )
#'   plot(
#'     fit,
#'     type = "ccc",
#'     title = "Customized Category Characteristic Curves",
#'     palette = c("#1b9e77", "#d95f02", "#7570b3")
#'   )
#' }
#' @export
plot.mfrm_fit <- function(x,
                          type = NULL,
                          facet = NULL,
                          top_n = 30,
                          theta_range = c(-6, 6),
                          theta_points = 241,
                          title = NULL,
                          palette = NULL,
                          label_angle = 45,
                          show_ci = FALSE,
                          ci_level = 0.95,
                          draw = TRUE,
                          preset = c("standard", "publication", "compact"),
                          ...) {
  if (!inherits(x, "mfrm_fit")) {
    stop("`x` must be an mfrm_fit object from fit_mfrm().")
  }
  top_n <- max(1L, as.integer(top_n))
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  style <- resolve_plot_preset(preset)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }

  as_plot_data <- function(name, data) new_mfrm_plot_data(name, data)

  se_tbl_ci <- if (isTRUE(show_ci)) compute_se_for_plot(x, ci_level = ci_level) else NULL

  default_bundle <- missing(type) || is.null(type)
  if (default_bundle || tolower(as.character(type[1])) %in% c("bundle", "all", "default")) {
    out <- list(
      wright_map = as_plot_data("wright_map", c(
        build_wright_map_data(x, top_n = top_n, se_tbl = se_tbl_ci),
        list(
          title = title %||% "Wright map",
          subtitle = "Shared logit scale for persons, facets, and thresholds",
          preset = style$name,
          legend = new_plot_legend(
            label = c("Person density", "Facet levels", "Step thresholds"),
            role = c("distribution", "location", "location"),
            aesthetic = c("histogram", "points", "points"),
            value = c(style$fill_muted, style$accent_tertiary, style$accent_secondary)
          ),
          reference_lines = new_reference_lines("h", 0, "Centered logit reference", "dashed", "reference")
        )
      )),
      pathway_map = as_plot_data("pathway_map", c(
        build_pathway_map_data(x, theta_range = theta_range, theta_points = theta_points),
        list(
          title = title %||% "Pathway map",
          subtitle = "Dominant score categories across the latent continuum",
          preset = style$name,
          legend = new_plot_legend(
            label = c("Step thresholds", "Dominant-category regions"),
            role = c("threshold", "region"),
            aesthetic = c("line", "band"),
            value = c(style$accent_secondary, style$fill_soft)
          ),
          reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
        )
      )),
      category_characteristic_curves = as_plot_data("category_characteristic_curves", c(
        build_ccc_data(x, theta_range = theta_range, theta_points = theta_points),
        list(
          title = title %||% "Category characteristic curves",
          subtitle = "Category response probabilities across theta",
          preset = style$name,
          legend = new_plot_legend(
            label = "Category curves",
            role = "probability",
            aesthetic = "line",
            value = "category_palette"
          ),
          reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
        )
      ))
    )
    class(out) <- c("mfrm_plot_bundle", class(out))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_wright_map(
        out$wright_map$data,
        title = title,
        palette = resolve_palette(
          palette = palette,
          defaults = c(
            facet_level = style$accent_tertiary,
            step_threshold = style$accent_secondary,
            person_hist = style$fill_muted,
            grid = style$grid
          )
        ),
        label_angle = label_angle,
        show_ci = show_ci,
        ci_level = ci_level
      )
      draw_pathway_map(
        out$pathway_map$data,
        title = title,
        palette = palette %||% c(style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn)
      )
      draw_ccc(
        out$category_characteristic_curves$data,
        title = title,
        palette = palette %||% c(style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn)
      )
    }
    return(invisible(out))
  }

  type <- match.arg(tolower(as.character(type[1])), c("facet", "person", "step", "wright", "pathway", "ccc"))

  if (type == "wright") {
    out <- as_plot_data("wright_map", c(
      build_wright_map_data(x, top_n = top_n, se_tbl = se_tbl_ci),
      list(
        title = title %||% "Wright map",
        subtitle = "Shared logit scale for persons, facets, and thresholds",
        preset = style$name,
        legend = new_plot_legend(
          label = c("Person density", "Facet levels", "Step thresholds"),
          role = c("distribution", "location", "location"),
          aesthetic = c("histogram", "points", "points"),
          value = c(style$fill_muted, style$accent_tertiary, style$accent_secondary)
        ),
        reference_lines = new_reference_lines("h", 0, "Centered logit reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_wright_map(
        out$data,
        title = title,
        palette = resolve_palette(
          palette = palette,
          defaults = c(
            facet_level = style$accent_tertiary,
            step_threshold = style$accent_secondary,
            person_hist = style$fill_muted,
            grid = style$grid
          )
        ),
        label_angle = label_angle,
        show_ci = show_ci,
        ci_level = ci_level
      )
    }
    return(invisible(out))
  }
  if (type == "pathway") {
    out <- as_plot_data("pathway_map", c(
      build_pathway_map_data(x, theta_range = theta_range, theta_points = theta_points),
      list(
        title = title %||% "Pathway map",
        subtitle = "Dominant score categories across the latent continuum",
        preset = style$name,
        legend = new_plot_legend(
          label = c("Step thresholds", "Dominant-category regions"),
          role = c("threshold", "region"),
          aesthetic = c("line", "band"),
          value = c(style$accent_secondary, style$fill_soft)
        ),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_pathway_map(out$data, title = title, palette = palette %||% c(
        style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn
      ))
    }
    return(invisible(out))
  }
  if (type == "ccc") {
    out <- as_plot_data("category_characteristic_curves", c(
      build_ccc_data(x, theta_range = theta_range, theta_points = theta_points),
      list(
        title = title %||% "Category characteristic curves",
        subtitle = "Category response probabilities across theta",
        preset = style$name,
        legend = new_plot_legend(
          label = "Category curves",
          role = "probability",
          aesthetic = "line",
          value = "category_palette"
        ),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_ccc(out$data, title = title, palette = palette %||% c(
        style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn
      ))
    }
    return(invisible(out))
  }
  if (type == "person") {
    person_tbl <- tibble::as_tibble(x$facets$person)
    person_tbl <- person_tbl[is.finite(person_tbl$Estimate), , drop = FALSE]
    if (nrow(person_tbl) == 0) stop("No finite person estimates available for plotting.")
    bins <- max(10L, min(35L, as.integer(round(sqrt(nrow(person_tbl))))))
    this_title <- if (is.null(title)) "Person measure distribution" else as.character(title[1])
    out <- as_plot_data("person", list(
      person = person_tbl,
      bins = bins,
      title = this_title,
      subtitle = "Distribution of person measures on the logit scale",
      legend = new_plot_legend("Person distribution", "distribution", "histogram", style$accent_primary),
      reference_lines = new_reference_lines(),
      preset = style$name
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_person_plot(
        person_tbl,
        bins = bins,
        title = this_title,
        palette = resolve_palette(palette = palette, defaults = c(person_hist = style$accent_primary))
      )
    }
    return(invisible(out))
  }
  if (type == "step") {
    step_tbl <- tibble::as_tibble(x$steps)
    if (nrow(step_tbl) == 0 || !all(c("Step", "Estimate") %in% names(step_tbl))) {
      stop("Step estimates are not available in this fit object.")
    }
    this_title <- if (is.null(title)) "Step parameter estimates" else as.character(title[1])
    out <- as_plot_data("step", list(
      steps = step_tbl,
      title = this_title,
      subtitle = "Estimated step thresholds on the shared logit scale",
      legend = new_plot_legend("Step thresholds", "threshold", "line-point", style$accent_tertiary),
      reference_lines = new_reference_lines("h", 0, "Centered step reference", "dashed", "reference"),
      preset = style$name
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_step_plot(
        step_tbl,
        title = this_title,
        palette = resolve_palette(
          palette = palette,
          defaults = c(step_line = style$accent_tertiary, grid = style$grid)
        ),
        label_angle = label_angle
      )
    }
    return(invisible(out))
  }

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (nrow(facet_tbl) == 0 || !all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available in this fit object.")
  }
  if (!is.null(facet)) {
    facet_tbl <- dplyr::filter(facet_tbl, .data$Facet == as.character(facet[1]))
    if (nrow(facet_tbl) == 0) stop("Requested `facet` not found in facet-level estimates.")
  }
  if (isTRUE(show_ci) && !is.null(se_tbl_ci) && is.data.frame(se_tbl_ci) &&
      nrow(se_tbl_ci) > 0 && all(c("Facet", "Level", "SE") %in% names(se_tbl_ci))) {
    se_join <- se_tbl_ci[, intersect(c("Facet", "Level", "SE"), names(se_tbl_ci)), drop = FALSE]
    se_join$Level <- as.character(se_join$Level)
    facet_tbl$Level <- as.character(facet_tbl$Level)
    facet_tbl <- merge(facet_tbl, se_join, by = c("Facet", "Level"), all.x = TRUE, sort = FALSE)
    facet_tbl <- tibble::as_tibble(facet_tbl)
  }
  facet_tbl <- dplyr::arrange(facet_tbl, .data$Estimate)
  if (nrow(facet_tbl) > top_n) {
    facet_tbl <- facet_tbl |>
      dplyr::mutate(AbsEstimate = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$AbsEstimate)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::arrange(.data$Estimate) |>
      dplyr::select(-"AbsEstimate")
  }
  facet_title <- if (is.null(facet)) "Facet-level estimates" else paste0("Facet-level estimates: ", as.character(facet[1]))
  if (!is.null(title)) facet_title <- as.character(title[1])
  out <- as_plot_data("facet", list(
    facets = facet_tbl,
    title = facet_title,
    subtitle = if (is.null(facet)) "Ranked facet-level measures" else paste0("Ranked measures for facet: ", as.character(facet[1])),
    legend = new_plot_legend("Facet estimates", "location", "bar", style$fill_soft),
    reference_lines = new_reference_lines("h", 0, "Centered facet reference", "dashed", "reference"),
    preset = style$name
  ))
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    draw_facet_plot(
      facet_tbl,
      title = facet_title,
      palette = resolve_palette(palette = palette, defaults = c(facet_bar = style$fill_soft)),
      show_ci = show_ci,
      ci_level = ci_level
    )
  }
  invisible(out)
}
