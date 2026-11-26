# ==============================================================================
# Unified Color-Vision-Safe ggplot2 `plot()` Methods for `mfrmr`
# Supports Okabe-Ito, Viridis, and Strict Monochrome + Redundant Shape/Linetype
# ==============================================================================

.mfrm_gg_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
      plot.background = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#EBEBEB", linewidth = 0.35),
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", colour = "#CCCCCC"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#222222"),
      legend.background = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
      legend.key = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

.mfrm_cvd_palette <- function(n, palette = c("okabe_ito", "viridis", "monochrome")) {
  palette <- match.arg(palette)
  if (n <= 0L) return(character(0))
  if (identical(palette, "monochrome")) {
    if (n == 1L) return("#222222")
    return(grDevices::gray(seq(0.12, 0.65, length.out = n)))
  }
  if (identical(palette, "viridis")) {
    if (n == 1L) return("#21918C")
    return(grDevices::hcl.colors(n, palette = "viridis", rev = FALSE))
  }
  oi <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#000000", "#F0E442"
  )
  if (n <= length(oi)) oi[seq_len(n)] else grDevices::colorRampPalette(oi)(n)
}

.mfrm_cvd_shapes <- function(n) {
  sh <- c(16L, 17L, 15L, 18L, 3L, 4L, 8L, 1L, 2L, 0L)
  rep(sh, length.out = max(1L, n))[seq_len(n)]
}

.mfrm_cvd_linetypes <- function(n) {
  lt <- c("solid", "dashed", "dotdash", "twodash", "longdash", "dotted")
  rep(lt, length.out = max(1L, n))[seq_len(n)]
}

.mfrm_gg_wright <- function(person_df, facets_df, steps_df, model_lbl, palette = "okabe_ito", base_size = 11) {
  pt <- person_df[is.finite(person_df$Estimate), , drop = FALSE]
  fc <- facets_df[is.finite(facets_df$Estimate), , drop = FALSE]
  st <- steps_df[is.finite(steps_df$Estimate), , drop = FALSE]

  p_df <- data.frame(Panel = "Person", Element = pt$Person, Logit = pt$Estimate, stringsAsFactors = FALSE)
  f_df <- data.frame(Panel = fc$Facet, Element = fc$Level, Logit = fc$Estimate, stringsAsFactors = FALSE)
  s_df <- if (nrow(st) > 0L) {
    lbl <- if (identical(model_lbl, "PCM")) paste0(st$ScaleLevel, ":", st$Step) else as.character(st$Step)
    data.frame(Panel = "Thresholds", Element = lbl, Logit = st$Estimate, stringsAsFactors = FALSE)
  } else {
    data.frame(Panel = character(0), Element = character(0), Logit = numeric(0))
  }

  panel_levels <- unique(c("Person", unique(fc$Facet), if (nrow(s_df) > 0L) "Thresholds"))
  all_df <- rbind(p_df, f_df, s_df)
  all_df$Panel <- factor(all_df$Panel, levels = panel_levels)

  cols <- stats::setNames(.mfrm_cvd_palette(length(panel_levels), palette), panel_levels)
  shapes <- stats::setNames(.mfrm_cvd_shapes(length(panel_levels)), panel_levels)

  non_p <- all_df[all_df$Panel != "Person", , drop = FALSE]
  p_only <- all_df[all_df$Panel == "Person", , drop = FALSE]

  hist_df <- data.frame(Panel = factor(character(0), levels = panel_levels),
                        ymin = numeric(0), ymax = numeric(0), xmin = numeric(0), xmax = numeric(0))
  if (nrow(p_only) >= 2L) {
    rng <- range(all_df$Logit, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) < 1e-6) rng <- rng + c(-1, 1)
    brks <- seq(rng[1L] - 0.05, rng[2L] + 0.05, length.out = 19L)
    h <- graphics::hist(p_only$Logit, breaks = brks, plot = FALSE)
    max_c <- max(h$counts, 1L)
    hist_df <- data.frame(
      Panel = factor("Person", levels = panel_levels),
      ymin = h$breaks[-length(h$breaks)],
      ymax = h$breaks[-1L],
      xmin = 1.0 - 0.38 * (h$counts / max_c),
      xmax = 1.0 + 0.38 * (h$counts / max_c),
      Count = h$counts,
      stringsAsFactors = FALSE
    )
    hist_df <- hist_df[hist_df$Count > 0L, , drop = FALSE]
  }
  non_p$XPos <- as.numeric(non_p$Panel)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888", linewidth = 0.4)
  if (nrow(hist_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = hist_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = cols[["Person"]], alpha = 0.30, colour = cols[["Person"]], linewidth = 0.3
    )
  }
  if (nrow(non_p) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data = non_p,
        ggplot2::aes(x = .data$XPos, y = .data$Logit, colour = .data$Panel, shape = .data$Panel),
        size = 2.7
      ) +
      ggplot2::geom_text(
        data = non_p,
        ggplot2::aes(x = .data$XPos + 0.08, y = .data$Logit, label = .data$Element),
        hjust = 0, size = 3.0, colour = "#222222", check_overlap = TRUE
      )
  }
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_along(panel_levels), labels = panel_levels,
      limits = c(0.45, length(panel_levels) + 0.65)
    ) +
    ggplot2::scale_colour_manual(values = cols, name = "Facet") +
    ggplot2::scale_shape_manual(values = shapes, name = "Facet") +
    ggplot2::labs(x = NULL, y = "Logit Scale") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- c(
    sprintf("Threshold structure: %s.", model_lbl),
    "Facet elements are encoded simultaneously by horizontal column, colour, and point shape."
  )
  p
}

.mfrm_gg_pathway <- function(facets_df, mnsq_range = c(0.70, 1.30), palette = "okabe_ito", base_size = 11,
                             fit_statistic = c("infit", "outfit", "both")) {
  fit_statistic <- match.arg(fit_statistic)
  statistics <- if (fit_statistic == "both") c("Infit", "Outfit") else tools::toTitleCase(fit_statistic)
  df <- do.call(rbind, lapply(statistics, function(statistic) {
    vals <- facets_df[[paste0(statistic, "MnSq")]]
    keep <- is.finite(facets_df$Estimate) & is.finite(vals)
    data.frame(facets_df[keep, , drop = FALSE],
               FitMnSq = vals[keep], Statistic = rep(statistic, sum(keep)))
  }))
  se_use <- if ("RealSE" %in% names(df)) df$RealSE else df$ModelSE
  df$Low <- df$Estimate - 1.96 * se_use
  df$High <- df$Estimate + 1.96 * se_use
  df$FitBand <- ifelse(
    df$FitMnSq > mnsq_range[2L], "Underfit (MnSq > Upper)",
    ifelse(df$FitMnSq < mnsq_range[1L], "Overfit / Muted (MnSq < Lower)", "Acceptable")
  )
  bands <- c("Acceptable", "Overfit / Muted (MnSq < Lower)", "Underfit (MnSq > Upper)")
  df$FitBand <- factor(df$FitBand, levels = bands)
  cols <- stats::setNames(.mfrm_cvd_palette(3L, palette), bands)
  shapes <- stats::setNames(c(16L, 15L, 17L), bands)
  ltypes <- stats::setNames(c("solid", "dashed", "dotdash"), bands)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$FitMnSq, y = .data$Estimate)) +
    ggplot2::geom_vline(xintercept = 1.0, linetype = "solid", colour = "#999999", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = mnsq_range, linetype = "dashed", colour = "#666666", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.0, linetype = "dotted", colour = "#888888", linewidth = 0.35) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$Low, ymax = .data$High, colour = .data$FitBand, linetype = .data$FitBand),
      width = 0.0, linewidth = 0.55, alpha = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$FitBand, shape = .data$FitBand),
      size = 2.9
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Level),
      nudge_y = 0.04 * if (nrow(df)) max(diff(range(df$Estimate)), 0.5) else 0.5,
      size = 2.9, colour = "#222222", check_overlap = TRUE
    ) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Fit Screening") +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE, name = "Fit Screening") +
    ggplot2::scale_linetype_manual(values = ltypes, drop = FALSE, name = "Fit Screening") +
    ggplot2::labs(x = if (length(statistics) == 1L) paste(statistics, "Mean-Square (MnSq)") else "Mean-Square (MnSq)",
                  y = "Measure (Logits)") +
    .mfrm_gg_theme(base_size = base_size)
  if (fit_statistic == "both" && nrow(df)) {
    p <- p + ggplot2::facet_grid(Facet ~ Statistic)
  } else if (length(unique(df$Facet)) > 1L) {
    p <- p + ggplot2::facet_wrap(~ Facet)
  }
  attr(p, "mfrm_notes") <- c(
    sprintf("Dashed vertical reference lines mark the displayed MnSq bounds [%.2f, %.2f].", mnsq_range[1L], mnsq_range[2L]),
    "Error bars show normal 95% intervals when available. The selected fit statistic determines the horizontal position."
  )
  p
}

.mfrm_gg_ccc <- function(step_cum, steps_df, rating_min, step_cat_min, step_cat_max,
                         slopes_df = NULL, step_facet = NULL, slope_formula = 0L,
                         palette = "okabe_ito", base_size = 11) {
  C <- ncol(step_cum)
  n_sg <- nrow(step_cum)
  sg_names <- unique(steps_df$ScaleLevel)
  if (length(sg_names) != n_sg) sg_names <- paste0("Group_", seq_len(n_sg))
  cat_labels <- as.character(rating_min + seq(0L, C - 1L))
  theta_grid <- seq(-4.5, 4.5, length.out = 121L)

  rows <- vector("list", n_sg)
  for (g in seq_len(n_sg)) {
    a_g <- 1.0
    if (!is.null(slopes_df) && nrow(slopes_df) > 0L && !is.null(step_facet)) {
      sl_hit <- slopes_df[slopes_df$Facet == step_facet & slopes_df$Level == sg_names[g], , drop = FALSE]
      if (nrow(sl_hit) == 1L && is.finite(sl_hit$Estimate[1L])) a_g <- sl_hit$Estimate[1L]
    }
    ev <- .mfrm_eval_row_predictions(
      eta_vec = theta_grid,
      step_cum = step_cum,
      step_group0 = rep(g - 1L, length(theta_grid)),
      step_cat_min = step_cat_min,
      step_cat_max = step_cat_max,
      rating_min = rating_min,
      alpha_vec = if (!is.null(slopes_df) && nrow(slopes_df) > 0L) rep(a_g, length(theta_grid)) else NULL,
      theta_vec = theta_grid,
      base_eta_vec = numeric(length(theta_grid)),
      slope_formula = slope_formula
    )
    df_g <- expand.grid(ThetaIdx = seq_along(theta_grid), CatIdx = seq_len(C), KEEP.OUT.ATTRS = FALSE)
    df_g$ScaleGroup <- sg_names[g]
    df_g$LinearMeasure <- theta_grid[df_g$ThetaIdx]
    df_g$Category <- factor(cat_labels[df_g$CatIdx], levels = cat_labels)
    df_g$Probability <- as.numeric(ev$probs)
    rows[[g]] <- df_g
  }
  ccc_df <- do.call(rbind, rows)
  # Subsample points along theta so distinct point shapes reinforce linetypes
  pts_df <- ccc_df[ccc_df$ThetaIdx %in% seq(1L, length(theta_grid), by = 15L), , drop = FALSE]

  cols <- stats::setNames(.mfrm_cvd_palette(C, palette), cat_labels)
  ltypes <- stats::setNames(.mfrm_cvd_linetypes(C), cat_labels)
  shapes <- stats::setNames(.mfrm_cvd_shapes(C), cat_labels)

  st_lines <- steps_df[is.finite(steps_df$Estimate), , drop = FALSE]
  st_lines$ScaleGroup <- st_lines$ScaleLevel

  p <- ggplot2::ggplot(
    ccc_df,
    ggplot2::aes(x = .data$LinearMeasure, y = .data$Probability,
                 colour = .data$Category, linetype = .data$Category)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(
      data = pts_df,
      ggplot2::aes(shape = .data$Category),
      size = 2.0
    ) +
    ggplot2::geom_vline(
      data = st_lines,
      ggplot2::aes(xintercept = .data$Estimate),
      linetype = "dotted", colour = "#666666", linewidth = 0.4
    ) +
    ggplot2::scale_colour_manual(values = cols, name = "Category") +
    ggplot2::scale_linetype_manual(values = ltypes, name = "Category") +
    ggplot2::scale_shape_manual(values = shapes, name = "Category") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    ggplot2::labs(x = "Latent Measure Relative to Element Location (Logits)",
                  y = "Category Response Probability") +
    .mfrm_gg_theme(base_size = base_size)
  if (n_sg > 1L) p <- p + ggplot2::facet_wrap(~ ScaleGroup)
  attr(p, "mfrm_notes") <- "Categories are redundantly encoded by colour, linetype, and marker shape."
  p
}

.mfrm_gg_person_fit <- function(fit_obj, palette = "okabe_ito", base_size = 11, pfi = NULL) {
  if (is.null(pfi)) pfi <- compute_person_fit_indices(fit_obj)
  df <- as.data.frame(pfi)
  df <- df[is.finite(df$OutfitMnSq) & is.finite(df$ReportValue), , drop = FALSE]
  lvls <- c("Not flagged", "Outfit > 1.50", "Person-fit flag", "Both flagged")
  st_chr <- ifelse(df$OutfitMnSq > 1.50 & df$ReportFlag, "Both flagged",
    ifelse(df$ReportFlag, "Person-fit flag", ifelse(df$OutfitMnSq > 1.50, "Outfit > 1.50", "Not flagged")))
  df$Status <- factor(st_chr, levels = lvls)
  cols <- stats::setNames(.mfrm_cvd_palette(4L, palette), lvls)
  shapes <- stats::setNames(c(16L, 17L, 15L, 8L), lvls)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$OutfitMnSq, y = .data$ReportValue)) +
    ggplot2::geom_vline(xintercept = 1.50, linetype = "dashed", colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = c(-1, 1) * stats::qnorm(.975), linetype = "dotdash", colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Status, shape = .data$Status), size = 2.6, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Person Fit") +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE, name = "Person Fit") +
    ggplot2::labs(x = "Person Outfit Mean-Square (MnSq)", y = "Reported person-fit index (lz* or lz)") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- c(
    "Uses the same ReportIndex, ReportValue and two-sided normal-reference flags as summary(compute_person_fit_indices(fit)).",
    "These screening cutoffs do not provide calibrated false-positive rates. Outfit > 1.50 is shown separately.",
    unique(pfi$ReportCaveat))
  p
}

.mfrm_gg_residual_qq <- function(fit_obj, palette = "okabe_ito", base_size = 11) {
  z_vec <- stats::residuals(fit_obj, type = "pearson")
  z_vec <- sort(z_vec[is.finite(z_vec)])
  n <- length(z_vec)
  if (n < 2L) stop("Insufficient finite standardized residuals for Q-Q plot.", call. = FALSE)
  pp <- (seq_len(n) - 0.5) / n
  theo <- stats::qnorm(pp)
  se_env <- sqrt(pp * (1.0 - pp) / n) / pmax(1e-8, stats::dnorm(theo))
  lvls <- c("|z| < 2.0", "|z| >= 2.0")
  df <- data.frame(
    Theoretical = theo,
    Observed = z_vec,
    EnvLow = theo - 1.96 * se_env,
    EnvHigh = theo + 1.96 * se_env,
    Band = factor(ifelse(abs(z_vec) >= 2.0, "|z| >= 2.0", "|z| < 2.0"), levels = lvls)
  )
  cols <- stats::setNames(.mfrm_cvd_palette(2L, palette), lvls)
  shapes <- stats::setNames(c(16L, 17L), lvls)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Theoretical, y = .data$Observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$EnvLow, ymax = .data$EnvHigh),
      fill = "#CCCCCC", alpha = 0.28, colour = "#888888", linetype = "dotted", linewidth = 0.35
    ) +
    ggplot2::geom_abline(intercept = 0.0, slope = 1.0, linetype = "dashed", colour = "#555555", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Band, shape = .data$Band), size = 1.9, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Residual Tail") +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE, name = "Residual Tail") +
    ggplot2::labs(x = "Theoretical Standard Normal Quantiles", y = "Standardized Response Residuals (z)") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- "Dotted band shows pointwise 95% normal order-statistic reference envelope; colour and shape distinguish |z| >= 2.0."
  p
}

.mfrm_gg_slopes <- function(slopes_df, palette = "okabe_ito", base_size = 11, intervals = NULL) {
  if (is.null(slopes_df) || nrow(slopes_df) == 0L) {
    stop("`type = 'slopes'` requires a fitted GMFRM (`model = 'GPCM'`) with estimated discrimination parameters.", call. = FALSE)
  }
  keep <- is.finite(slopes_df$Estimate)
  df <- slopes_df[keep, , drop = FALSE]
  df$Low <- if (is.null(intervals)) pmax(0.0, df$Estimate - 1.96 * df$ModelSE) else intervals[keep, 1L]
  df$High <- if (is.null(intervals)) df$Estimate + 1.96 * df$ModelSE else intervals[keep, 2L]
  lvls <- c("Low (< 0.80)", "Moderate (0.80 - 1.25)", "High (> 1.25)")
  df$Band <- factor(
    ifelse(df$Estimate < 0.80, "Low (< 0.80)", ifelse(df$Estimate > 1.25, "High (> 1.25)", "Moderate (0.80 - 1.25)")),
    levels = lvls
  )
  df$Elem <- factor(df$Level, levels = rev(unique(df$Level)))
  cols <- stats::setNames(.mfrm_cvd_palette(3L, palette), lvls)
  shapes <- stats::setNames(c(15L, 16L, 17L), lvls)
  ltypes <- stats::setNames(c("dashed", "solid", "dotdash"), lvls)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Estimate, y = .data$Elem)) +
    ggplot2::geom_vline(xintercept = 1.0, linetype = "dotted", colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$Low, xmax = .data$High, colour = .data$Band, linetype = .data$Band),
      orientation = "y", width = 0.0, linewidth = 0.6
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Band, shape = .data$Band), size = 2.8) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Discrimination") +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE, name = "Discrimination") +
    ggplot2::scale_linetype_manual(values = ltypes, drop = FALSE, name = "Discrimination") +
    ggplot2::labs(x = "Discrimination Parameter Alpha (95% Model-SE CI)", y = NULL) +
    .mfrm_gg_theme(base_size = base_size)
  if (length(unique(df$Facet)) > 1L) p <- p + ggplot2::facet_wrap(~ Facet, scales = "free_y")
  attr(p, "mfrm_notes") <- "Reference line at alpha = 1.0 denotes equal-slope Rasch discrimination; colour, shape, and linetype encode discrimination band."
  p
}

.mfrm_gg_qc <- function(fit_obj, palette = "okabe_ito", base_size = 11) {
  fc <- fit_obj$facets[is.finite(fit_obj$facets$InfitMnSq), , drop = FALSE]
  pt <- fit_obj$person[is.finite(fit_obj$person$InfitMnSq), , drop = FALSE]
  st <- fit_obj$steps[is.finite(fit_obj$steps$Estimate), , drop = FALSE]
  zr <- stats::residuals(fit_obj, type = "pearson")
  zr <- zr[is.finite(zr)]

  lvls <- c("Within Reference", "Flagged")
  d1 <- data.frame(
    Panel = "1. Facet Infit MnSq",
    X = seq_len(nrow(fc)),
    Y = fc$InfitMnSq,
    Status = ifelse(fc$InfitMnSq < 0.70 | fc$InfitMnSq > 1.30, "Flagged", "Within Reference"),
    stringsAsFactors = FALSE
  )
  d2 <- data.frame(
    Panel = "2. Person Outfit MnSq",
    X = pt$Estimate,
    Y = pt$OutfitMnSq,
    Status = ifelse(pt$OutfitMnSq > 1.50, "Flagged", "Within Reference"),
    stringsAsFactors = FALSE
  )
  d3 <- if (nrow(st) > 0L) {
    dis_vec <- ave(st$Estimate, st$ScaleLevel, FUN = function(v) {
      if (length(v) < 2L) FALSE else c(FALSE, diff(v) < 0)
    })
    data.frame(
      Panel = "3. Step Thresholds (Logits)",
      X = seq_len(nrow(st)),
      Y = st$Estimate,
      Status = ifelse(dis_vec > 0, "Flagged", "Within Reference"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Panel = character(0), X = numeric(0), Y = numeric(0), Status = character(0))
  }
  pp <- (seq_along(zr) - 0.5) / max(1L, length(zr))
  d4 <- data.frame(
    Panel = "4. Residual Q-Q",
    X = stats::qnorm(pp),
    Y = sort(zr),
    Status = ifelse(abs(sort(zr)) >= 2.0, "Flagged", "Within Reference"),
    stringsAsFactors = FALSE
  )
  all_qc <- rbind(d1, d2, d3, d4)
  all_qc$Status <- factor(all_qc$Status, levels = lvls)
  all_qc$Panel <- factor(all_qc$Panel, levels = unique(all_qc$Panel))
  cols <- stats::setNames(.mfrm_cvd_palette(2L, palette), lvls)
  shapes <- stats::setNames(c(16L, 17L), lvls)

  p <- ggplot2::ggplot(all_qc, ggplot2::aes(x = .data$X, y = .data$Y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Status, shape = .data$Status), size = 2.0, alpha = 0.85) +
    ggplot2::facet_wrap(~ Panel, scales = "free") +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "QC Screening") +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE, name = "QC Screening") +
    ggplot2::labs(x = "Index / Latent Scale / Quantile", y = "Diagnostic Statistic") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- "Four-panel quality-control summary; screening status is encoded by colour and point shape."
  p
}




















plot.mfrm_fit <- function(x, type = c("wright", "pathway", "steps", "category_curves", "persons", "raters",
    "person_fit", "residual_qq", "slopes", "qc"), palette = c("okabe_ito", "viridis", "monochrome"), facet = NULL,
    base_size = 11, ..., fit_statistic = c("infit", "outfit", "both")) {
    .mfrm_validate_fit(x)
    type <- match.arg(type)
    palette <- match.arg(palette)
    interval_ready <- TRUE
    interval_notes <- character()
    if (type %in% c("pathway", "raters", "persons", "slopes")) {
        readiness <- .mfrm_readiness(x)
        interval_ready <- readiness$ready
        if (!interval_ready) {
            interval_notes <- paste0("Intervals unavailable: ", paste(readiness$reason_codes, collapse = ", "),
                ". Points are descriptive estimates requiring review.")
            for (component in c("person", "facets", "slopes")) {
                if (is.data.frame(x[[component]]))
                  for (se in intersect(c("ModelSE", "RealSE"), names(x[[component]]))) x[[component]][[se]] <- NA_real_
            }
        }
    }
    with_interval_notes <- function(p) {
        attr(p, "mfrm_notes") <- c(attr(p, "mfrm_notes"), interval_notes)
        p
    }
    if (identical(type, "wright")) {
        label <- .mfrm_model_label(x$model, unique(x$slopes$Facet))
        return(.mfrm_gg_wright(x$person, x$facets, x$steps, label, palette = palette, base_size = base_size))
    }
    if (type %in% c("pathway", "raters")) {
        fc <- x$facets
        if (identical(type, "raters") && is.null(facet)) {
            hit <- grep("rater|judge|examiner|marker|grader", unique(fc$Facet), ignore.case = TRUE, value = TRUE)
            facet <- if (length(hit) > 0L)
                hit[1L]
            else unique(fc$Facet)[1L]
        }
        if (!is.null(facet))
            fc <- fc[fc$Facet %in% facet, , drop = FALSE]
        return(with_interval_notes(.mfrm_gg_pathway(fc, palette = palette, base_size = base_size, fit_statistic = match.arg(fit_statistic))))
    }
    if (type %in% c("steps", "category_curves")) {
        return(.mfrm_gg_ccc(x$step_cumulative, x$steps, x$rating_min, x$encoding$step_cat_min, x$encoding$step_cat_max,
            slopes_df = x$slopes, step_facet = x$encoding$step_facet, slope_formula = x$encoding$slope_formula %||%
                0L, palette = palette, base_size = base_size))
    }
    if (identical(type, "person_fit")) {
        return(.mfrm_gg_person_fit(x, palette = palette, base_size = base_size))
    }
    if (identical(type, "residual_qq")) {
        return(.mfrm_gg_residual_qq(x, palette = palette, base_size = base_size))
    }
    if (identical(type, "slopes")) {
        ci <- if (interval_ready)
            stats::confint(x, component = "slopes")
        else matrix(NA_real_, NROW(x$slopes), 2L)
        return(with_interval_notes(.mfrm_gg_slopes(x$slopes, palette = palette, base_size = base_size, intervals = ci)))
    }
    if (identical(type, "qc")) {
        return(.mfrm_gg_qc(x, palette = palette, base_size = base_size))
    }
    pt <- x$person[is.finite(x$person$Estimate) & is.finite(x$person$FairM), , drop = FALSE]
    pt$Low <- pt$Estimate - 1.96 * pt$RealSE
    pt$High <- pt$Estimate + 1.96 * pt$RealSE
    pt$Band <- factor(ifelse(is.finite(pt$InfitMnSq) & pt$InfitMnSq > 1.3, "Infit MnSq > 1.30", "Infit MnSq <= 1.30"),
        levels = c("Infit MnSq <= 1.30", "Infit MnSq > 1.30"))
    cols <- stats::setNames(.mfrm_cvd_palette(2L, palette), levels(pt$Band))
    p <- ggplot2::ggplot(pt, ggplot2::aes(x = .data$Estimate, y = .data$FairM)) + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$Low,
        xmax = .data$High, colour = .data$Band, linetype = .data$Band), orientation = "y", width = 0, alpha = 0.65,
        linewidth = 0.45) + ggplot2::geom_point(ggplot2::aes(colour = .data$Band, shape = .data$Band), size = 2.3) +
        ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Person Fit") + ggplot2::scale_shape_manual(values = c(`Infit MnSq <= 1.30` = 16L,
        `Infit MnSq > 1.30` = 17L), drop = FALSE, name = "Person Fit") + ggplot2::scale_linetype_manual(values = c(`Infit MnSq <= 1.30` = "solid",
        `Infit MnSq > 1.30` = "dashed"), drop = FALSE, name = "Person Fit") + ggplot2::labs(x = "Person Ability (Logits; 95% Real-SE CI)",
        y = "Fair Average Score Fair(M)") + .mfrm_gg_theme(base_size = base_size)
    attr(p, "mfrm_notes") <- "Person fit status is encoded via colour, point shape, and interval linetype."
    with_interval_notes(p)
}






plot.mfrm_diagnostics <- function(x,
                                  facet = NULL,
                                  palette = c("okabe_ito", "viridis", "monochrome"),
                                  base_size = 11,
                                  ..., type = c("pathway", "residual_pca")) {
  palette <- match.arg(palette)
  type <- match.arg(type)
  if (type == "residual_pca") {
    if (length(facet) > 1L) stop("Choose one facet for the residual PCA plot.", call. = FALSE)
    result <- if (is.null(facet)) x$residual_pca_overall else x$residual_pca_by_facet[[facet]]
    if (is.null(result)) stop("Request residual PCA in diagnose_mfrm() for the plotted grouping.", call. = FALSE)
    if (!is.null(result$error)) stop(result$error, call. = FALSE)
    loadings <- as.matrix(result$pca$loadings)
    d <- data.frame(Level = rownames(loadings), Loading = loadings[, 1L])
    d$Level <- stats::reorder(d$Level, d$Loading)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Loading, y = .data$Level)) +
      ggplot2::geom_vline(xintercept = 0, color = "grey65") +
      ggplot2::geom_point(shape = 16, color = .mfrm_cvd_palette(1L, palette), size = 2.8) +
      ggplot2::labs(x = "Loading on first residual contrast", y = facet) +
      .mfrm_gg_theme(base_size = base_size)
    overlap <- crossprod(is.finite(as.matrix(result$residual_matrix)) * 1L)
    attr(p, "mfrm_notes") <- c(
      sprintf("First contrast eigenvalue: %.2f. Pairwise shared persons: %d\u2013%d.",
              result$pca$values[1L], min(overlap[upper.tri(overlap)]), max(overlap[upper.tri(overlap)])),
      "Exploratory PCA of averaged standardized residuals. Contrast signs are arbitrary.",
      result$warning
    )
    return(p)
  }
  readiness <- .mfrm_readiness(x$fit)
  ms <- as.data.frame(x$measures)
  if ("Infit" %in% names(ms) && !"InfitMnSq" %in% names(ms)) ms$InfitMnSq <- ms$Infit
  ms <- ms[ms$Facet != "Person", , drop = FALSE]
  if (!is.null(facet)) ms <- ms[ms$Facet %in% facet, , drop = FALSE]
  if (!readiness$ready) ms$ModelSE <- ms$RealSE <- NA_real_
  p <- .mfrm_gg_pathway(ms, palette = palette, base_size = base_size)
  if (!readiness$ready) attr(p, "mfrm_notes") <- c(attr(p, "mfrm_notes"),
    paste("Intervals unavailable:", paste(readiness$reason_codes, collapse = "; ")))
  p
}



.mfrm_gstudy_plot_henderson <- function(x,
                                   type = c("dstudy", "variance"),
                                   palette = c("okabe_ito", "viridis", "monochrome"),
                                   base_size = 11,
                                   ...) {
  x <- .mfrm_refresh_gstudy(x)
  type <- match.arg(type)
  palette <- match.arg(palette)
  if (identical(type, "variance")) {
    vc <- x$variance_components
    vc$CompOrder <- factor(vc$Component, levels = rev(vc$Component))
    vc$Role <- factor(ifelse(vc$ContainsPerson, "Universe / Relative", "Facet Main / Absolute"),
                      levels = c("Universe / Relative", "Facet Main / Absolute"))
    cols <- stats::setNames(.mfrm_cvd_palette(2L, palette), levels(vc$Role))
    p <- ggplot2::ggplot(vc, ggplot2::aes(x = .data$PercentVariance, y = .data$CompOrder)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data$PercentVariance, y = .data$CompOrder, yend = .data$CompOrder,
                     colour = .data$Role, linetype = .data$Role),
        linewidth = 0.8
      ) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$Role, shape = .data$Role), size = 3.0) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Variance Role") +
      ggplot2::scale_shape_manual(values = c(16L, 17L), drop = FALSE, name = "Variance Role") +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE, name = "Variance Role") +
      ggplot2::labs(x = "Percentage of Total Variance (%)", y = NULL) +
      .mfrm_gg_theme(base_size = base_size)
    attr(p, "mfrm_notes") <- "G-Study variance components estimated via Henderson Method 1 (urGENOVA)."
    return(p)
  }

  ds <- as.data.frame(x$d_study)
  f_cols <- x$facets
  p_f <- f_cols[1L]
  s_f <- if (length(f_cols) >= 2L) f_cols[2L] else NULL
  if (length(f_cols) >= 3L) {
    for (fc3 in f_cols[-c(1L, 2L)]) {
      u3 <- sort(unique(ds[[fc3]]))
      ds <- ds[ds[[fc3]] == u3[ceiling(length(u3) / 2)], , drop = FALSE]
    }
  }
  grp <- if (!is.null(s_f)) paste0(s_f, " = ", ds[[s_f]]) else "Design"
  long_df <- rbind(
    data.frame(PrimaryN = ds[[p_f]], Group = grp, Metric = "G (Relative)", Coef = ds$G, stringsAsFactors = FALSE),
    data.frame(PrimaryN = ds[[p_f]], Group = grp, Metric = "Phi (Absolute)", Coef = ds$Phi, stringsAsFactors = FALSE)
  )
  grp_lv <- unique(long_df$Group)
  long_df$Group <- factor(long_df$Group, levels = grp_lv)
  cols <- stats::setNames(.mfrm_cvd_palette(length(grp_lv), palette), grp_lv)
  shapes <- stats::setNames(.mfrm_cvd_shapes(length(grp_lv)), grp_lv)

  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = .data$PrimaryN, y = .data$Coef,
                 colour = .data$Group, shape = .data$Group, linetype = .data$Metric)
  ) +
    ggplot2::geom_hline(yintercept = 0.80, linetype = "dotted", colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = cols, name = if (!is.null(s_f)) s_f else "Design") +
    ggplot2::scale_shape_manual(values = shapes, name = if (!is.null(s_f)) s_f else "Design") +
    ggplot2::scale_linetype_manual(values = c("G (Relative)" = "solid", "Phi (Absolute)" = "dashed"), name = "Index") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(x = sprintf("Planned Sample Size (%s)", p_f), y = "Generalizability / Dependability") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- c("Secondary facet count is encoded by colour + point shape; G vs Phi by linetype.",
    unique(c(ds$GStatus[is.na(ds$G)], ds$PhiStatus[is.na(ds$Phi)])))
  p
}



plot.mfrm_bias <- function(x,
                           type = c("forest", "pairwise"),
                           top_n = 20L,
                           palette = c("okabe_ito", "viridis", "monochrome"),
                           base_size = 11,
                           ...) {
  type <- match.arg(type)
  palette <- match.arg(palette)
  if (identical(type, "pairwise") && is.data.frame(x$pairwise) && nrow(x$pairwise) > 0L) {
    df <- x$pairwise
    df <- df[order(-abs(df$t)), , drop = FALSE]
    df <- utils::head(df, as.integer(top_n))
    df$Label <- paste0(df$TargetLevel, ": ", df$Context1, " - ", df$Context2)
    df$Label <- factor(df$Label, levels = rev(unique(df$Label)))
    crit <- stats::qt(0.975, df = pmax(1, df$df))
    df$Low <- df$Contrast - crit * df$ContrastSE
    df$High <- df$Contrast + crit * df$ContrastSE
    df$Status <- factor(ifelse(is.na(df$pAdj), "Unavailable", ifelse(df$pAdj < 0.05, "Adjusted p < .05", "Not Flagged")),
                        levels = c("Not Flagged", "Adjusted p < .05", "Unavailable"))
    cols <- stats::setNames(.mfrm_cvd_palette(3L, palette), levels(df$Status))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Contrast, y = .data$Label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "#666666") +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$Low, xmax = .data$High, colour = .data$Status, linetype = .data$Status),
        orientation = "y", width = 0.0, linewidth = 0.6, na.rm = TRUE
      ) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$Status, shape = .data$Status), size = 2.8) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Pairwise Bias") +
      ggplot2::scale_shape_manual(values = c(16L, 17L, 4L), drop = FALSE, name = "Pairwise Bias") +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotted"), drop = FALSE, name = "Pairwise Bias") +
      ggplot2::labs(x = "Pairwise Bias Contrast (Logits; 95% Welch CI)", y = NULL) +
      .mfrm_gg_theme(base_size = base_size)
    attr(p, "mfrm_notes") <- c(sprintf("Degrees of freedom convention: %s.", x$test_label %||% x$df_method), unique(df$InferenceReason[nzchar(df$InferenceReason)]))
    return(p)
  }

  tb <- as.data.frame(x$table)
  est_col <- intersect(c("BiasSize", "Estimate", "Bias"), names(tb))[1L]
  se_col <- intersect(c("BiasSE", "ModelSE", "SE"), names(tb))[1L]
  tb$Est <- tb[[est_col]]
  tb$SE <- tb[[se_col]]
  tb <- tb[is.finite(tb$Est), , drop = FALSE]
  tb <- tb[order(-abs(tb$Est / pmax(1e-8, tb$SE))), , drop = FALSE]
  tb <- utils::head(tb, as.integer(top_n))
  tb$Cell <- paste0(tb$LevelA, " x ", tb$LevelB)
  tb$Cell <- factor(tb$Cell, levels = rev(unique(tb$Cell)))
  tb$Low <- tb$Est - 1.96 * tb$SE
  tb$High <- tb$Est + 1.96 * tb$SE
  p_col <- intersect(c("pAdj", "pValue"), names(tb))[1L]
  tb$Status <- factor(ifelse(is.na(tb[[p_col]]), "Unavailable", ifelse(tb[[p_col]] < 0.05, "Flagged (p < .05)", "Not Flagged")),
                      levels = c("Not Flagged", "Flagged (p < .05)", "Unavailable"))
  cols <- stats::setNames(.mfrm_cvd_palette(3L, palette), levels(tb$Status))

  p <- ggplot2::ggplot(tb, ggplot2::aes(x = .data$Est, y = .data$Cell)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "#666666") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$Low, xmax = .data$High, colour = .data$Status, linetype = .data$Status),
      orientation = "y", width = 0.0, linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Status, shape = .data$Status), size = 2.8) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "Interaction Cell") +
    ggplot2::scale_shape_manual(values = c(16L, 17L, 4L), drop = FALSE, name = "Interaction Cell") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotted"), drop = FALSE, name = "Interaction Cell") +
    ggplot2::labs(x = "Local Bias / Interaction Estimate (Logits; 95% CI)", y = NULL) +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- c("Cells ordered by absolute standardized bias; colour, shape, and linetype encode screening status.", unique(tb$InferenceReason[nzchar(tb$InferenceReason)]))
  p
}



plot.mfrm_dff <- function(x,
                          top_n = 25L,
                          palette = c("okabe_ito", "viridis", "monochrome"),
                          base_size = 11,
                          ...) {
  palette <- match.arg(palette)
  tb <- as.data.frame(x$table)
  c_col <- intersect(c("Contrast", "Delta", "Estimate"), names(tb))[1L]
  s_col <- intersect(c("SandwichSE", "SE", "ModelSE"), names(tb))
  se_vec <- tb[[s_col[1L]]]
  if (all(is.na(se_vec)) && length(s_col) >= 2L) se_vec <- tb[[s_col[2L]]]
  tb$Est <- tb[[c_col]]
  tb$SE_use <- se_vec
  tb <- utils::head(tb[order(-abs(tb$Est)), , drop = FALSE], as.integer(top_n))
  lv_col <- intersect(c("Level", "Element"), names(tb))[1L]
  tb$Elem <- factor(tb[[lv_col]], levels = rev(unique(tb[[lv_col]])))
  tb$Low <- tb$Est - 1.96 * tb$SE_use
  tb$High <- tb$Est + 1.96 * tb$SE_use
  flag_v <- if ("Flagged" %in% names(tb)) tb$Flagged else ifelse(is.na(tb$pAdj), NA, tb$pAdj < 0.05)
  tb$Status <- factor(ifelse(is.na(flag_v), "Unavailable", ifelse(flag_v, "Flagged DFF (p_adj < .05)", "No Flag")),
                      levels = c("No Flag", "Flagged DFF (p_adj < .05)", "Unavailable"))
  cols <- stats::setNames(.mfrm_cvd_palette(3L, palette), levels(tb$Status))

  p <- ggplot2::ggplot(tb, ggplot2::aes(x = .data$Est, y = .data$Elem)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "#666666") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$Low, xmax = .data$High, colour = .data$Status, linetype = .data$Status),
      orientation = "y", width = 0.0, linewidth = 0.6
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Status, shape = .data$Status), size = 2.8) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = "DFF Status") +
    ggplot2::scale_shape_manual(values = c(16L, 17L, 4L), drop = FALSE, name = "DFF Status") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotted"), drop = FALSE, name = "DFF Status") +
    ggplot2::labs(x = sprintf("DFF Contrast (%s - %s)", x$focal %||% "Focal", x$reference %||% "Reference"),
                  y = x$facet %||% "Element") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- c(sprintf("DFF method: %s.", x$method), unique(tb$InferenceReason[nzchar(tb$InferenceReason)]))
  p
}



plot.mfrm_information <- function(x,
                                  type = c("information", "ccc"),
                                  palette = c("okabe_ito", "viridis", "monochrome"),
                                  base_size = 11,
                                  ...) {
  type <- match.arg(type)
  palette <- match.arg(palette)
  if (identical(type, "ccc") && !is.null(x$category_curves)) {
    cc <- x$category_curves
    cats <- unique(cc$Category)
    cc$Category <- factor(cc$Category, levels = cats)
    cols <- stats::setNames(.mfrm_cvd_palette(length(cats), palette), cats)
    ltypes <- stats::setNames(.mfrm_cvd_linetypes(length(cats)), cats)
    shapes <- stats::setNames(.mfrm_cvd_shapes(length(cats)), cats)
    pts <- cc[seq(1L, nrow(cc), by = 15L), , drop = FALSE]
    p <- ggplot2::ggplot(cc, ggplot2::aes(x = .data$Theta, y = .data$Probability,
                                          colour = .data$Category, linetype = .data$Category)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(data = pts, ggplot2::aes(shape = .data$Category), size = 2.0) +
      ggplot2::scale_colour_manual(values = cols, name = "Category") +
      ggplot2::scale_linetype_manual(values = ltypes, name = "Category") +
      ggplot2::scale_shape_manual(values = shapes, name = "Category") +
      ggplot2::labs(x = "Latent Ability Theta (Logits)", y = "Category Probability") +
      .mfrm_gg_theme(base_size = base_size)
    if (length(unique(cc$ScaleLevel)) > 1L) p <- p + ggplot2::facet_wrap(~ ScaleLevel)
    return(p)
  }

  ti <- x$test_information %||% x$total
  th_c <- intersect(c("Theta", "theta"), names(ti))[1L]
  in_c <- intersect(c("Information", "information", "TotalInformation"), names(ti))[1L]
  se_c <- intersect(c("SE", "se", "ConditionalSE"), names(ti))[1L]
  df_l <- rbind(
    data.frame(Theta = ti[[th_c]], Value = ti[[in_c]], Curve = "Test Information I(theta)", stringsAsFactors = FALSE),
    data.frame(Theta = ti[[th_c]], Value = ti[[se_c]], Curve = "Conditional SE(theta)", stringsAsFactors = FALSE)
  )
  df_l <- df_l[is.finite(df_l$Value), , drop = FALSE]
  df_l$Curve <- factor(df_l$Curve, levels = c("Test Information I(theta)", "Conditional SE(theta)"))
  cols <- stats::setNames(.mfrm_cvd_palette(2L, palette), levels(df_l$Curve))
  pts <- df_l[seq(1L, nrow(df_l), by = 18L), , drop = FALSE]
  p <- ggplot2::ggplot(df_l, ggplot2::aes(x = .data$Theta, y = .data$Value,
                                          colour = .data$Curve, linetype = .data$Curve)) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(data = pts, ggplot2::aes(shape = .data$Curve), size = 2.2) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(16L, 17L), name = NULL) +
    ggplot2::facet_wrap(~ Curve, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Person ability (logits)", y = NULL) +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- c("Test information and conditional SE use separate vertical scales.",
                             if (!is.null(x$exposure_basis)) paste("Exposure basis:", x$exposure_basis),
                             x$notes)
  p
}



plot.mfrm_drift <- function(x,
                            type = c("trajectory", "contrasts"),
                            palette = c("okabe_ito", "viridis", "monochrome"),
                            base_size = 11,
                            ...) {
  type <- match.arg(type)
  palette <- match.arg(palette)
  sev <- as.data.frame(x, component = "severity")
  r_col <- intersect(c("Rater", "rater"), names(sev))[1L]
  t_col <- intersect(c("Time", "Wave", "time"), names(sev))[1L]
  e_col <- intersect(c("Estimate", "Severity", "estimate"), names(sev))[1L]
  s_col <- intersect(c("SE", "TotalSE", "se"), names(sev))[1L]
  raters <- unique(as.character(sev[[r_col]]))
  times <- unique(as.character(sev[[t_col]]))
  sev$RaterF <- factor(sev[[r_col]], levels = raters)
  sev$TimeF <- factor(sev[[t_col]], levels = times)
  sev$Est <- sev[[e_col]]
  sev$SE_u <- if (!is.na(s_col)) sev[[s_col]] else 0.0

  cols <- stats::setNames(.mfrm_cvd_palette(length(raters), palette), raters)
  shapes <- stats::setNames(.mfrm_cvd_shapes(length(raters)), raters)
  ltypes <- stats::setNames(.mfrm_cvd_linetypes(length(raters)), raters)

  p <- ggplot2::ggplot(
    sev,
    ggplot2::aes(x = .data$TimeF, y = .data$Est, group = .data$RaterF,
                 colour = .data$RaterF, linetype = .data$RaterF, shape = .data$RaterF)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "#888888") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = cols, name = "Rater") +
    ggplot2::scale_linetype_manual(values = ltypes, name = "Rater") +
    ggplot2::scale_shape_manual(values = shapes, name = "Rater") +
    ggplot2::labs(x = "Occasion / Wave", y = "Dynamic Rater Severity (Logits)") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- "Each rater trajectory is encoded simultaneously by colour, linetype, and marker shape."
  p
}






.mfrm_gstudy_plot_moments <- function(x,
                                          palette = c("okabe_ito", "viridis", "monochrome"),
                                          base_size = 11,
                                          ...) {
  palette <- match.arg(palette)
  srcs <- names(x$components)
  rows <- lapply(srcs, function(nm) {
    mat <- as.matrix(x$components[[nm]])
    sc <- colnames(mat) %||% paste0("Score", seq_len(ncol(mat)))
    data.frame(Source = nm, Score = sc, Variance = pmax(0.0, diag(mat)), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  sc_lvls <- unique(df$Score)
  df$Source <- factor(df$Source, levels = rev(srcs))
  df$Score <- factor(df$Score, levels = sc_lvls)
  cols <- stats::setNames(.mfrm_cvd_palette(length(sc_lvls), palette), sc_lvls)
  shapes <- stats::setNames(.mfrm_cvd_shapes(length(sc_lvls)), sc_lvls)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Variance, y = .data$Source,
                                        colour = .data$Score, shape = .data$Score)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "#888888") +
    ggplot2::geom_point(size = 3.0, position = ggplot2::position_dodge(width = 0.45)) +
    ggplot2::scale_colour_manual(values = cols, name = "Score") +
    ggplot2::scale_shape_manual(values = shapes, name = "Score") +
    ggplot2::labs(x = "Estimated Variance Component (Nonnegative-Clipped)", y = "Source of Variation") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- "Diagonal variance components from multivariate G-Study."
  p
}



plot.mfrm_multivariate_d_study <- function(x,
                                           palette = c("okabe_ito", "viridis", "monochrome"),
                                           base_size = 11,
                                           ...) {
  palette <- match.arg(palette)
  cf <- x$coefficients
  cf <- cf[is.finite(cf$G) | is.finite(cf$Phi), , drop = FALSE]
  scen <- paste0("S", cf$Scenario)
  df_l <- rbind(
    data.frame(Scenario = scen, Score = cf$Score, Metric = "G (Relative)", Value = cf$G, stringsAsFactors = FALSE),
    data.frame(Scenario = scen, Score = cf$Score, Metric = "Phi (Absolute)", Value = cf$Phi, stringsAsFactors = FALSE)
  )
  df_l <- df_l[is.finite(df_l$Value), , drop = FALSE]
  df_l$Scenario <- factor(df_l$Scenario, levels = unique(scen))
  df_l$Metric <- factor(df_l$Metric, levels = c("G (Relative)", "Phi (Absolute)"))
  cols <- stats::setNames(.mfrm_cvd_palette(2L, palette), levels(df_l$Metric))

  p <- ggplot2::ggplot(df_l, ggplot2::aes(x = .data$Scenario, y = .data$Value,
                                          group = .data$Metric, colour = .data$Metric,
                                          linetype = .data$Metric, shape = .data$Metric)) +
    ggplot2::geom_hline(yintercept = 0.80, linetype = "dotted", colour = "#666666") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.7) +
    ggplot2::scale_colour_manual(values = cols, name = "Coefficient") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = "Coefficient") +
    ggplot2::scale_shape_manual(values = c(16L, 17L), name = "Coefficient") +
    ggplot2::facet_wrap(~ Score) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "D-Study Scenario", y = "Generalizability / Dependability") +
    .mfrm_gg_theme(base_size = base_size)
  attr(p, "mfrm_notes") <- "Projected G and Phi across candidate D-Study scenarios."
  p
}
