# ==============================================================================
# Publication-Ready Color-Vision-Safe ggplot2 Visualizations for `tidymfrm`
# Supports Okabe-Ito, Viridis, and Strict Monochrome + Redundant Shape/Linetype
# ==============================================================================

.tidymfrm_theme <- function(base_size = 11) {
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

.tidymfrm_palette <- function(n, palette = c("okabe_ito", "viridis", "monochrome")) {
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

.tidymfrm_shapes <- function(n) {
  sh <- c(16L, 17L, 15L, 18L, 3L, 4L, 8L, 1L, 2L, 0L)
  rep(sh, length.out = max(1L, n))[seq_len(n)]
}

.tidymfrm_linetypes <- function(n) {
  lt <- c("solid", "dashed", "dotdash", "twodash", "longdash", "dotted")
  rep(lt, length.out = max(1L, n))[seq_len(n)]
}

.tidymfrm_flag_styles <- function(palette = c("okabe_ito", "viridis", "monochrome")) {
  palette <- match.arg(palette)
  lvls <- c("Acceptable", "Markedly Severe", "Markedly Lenient", "Muted / Centrality", "Misfit (Unpredictable)", "Review required")
  cols <- if (identical(palette, "monochrome")) {
    c("#333333", "#111111", "#666666", "#777777", "#000000", "#555555")
  } else if (identical(palette, "viridis")) {
    grDevices::hcl.colors(6L, palette = "viridis")
  } else {
    c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#666666")
  }
  list(
    colours = stats::setNames(cols, lvls),
    shapes = stats::setNames(c(16L, 17L, 15L, 18L, 8L, 4L), lvls),
    linetypes = stats::setNames(c("solid", "dashed", "dotdash", "twodash", "longdash", "dotted"), lvls)
  )
}

.tidymfrm_plot_raters <- function(raters_df, mnsq_range = c(0.70, 1.30), palette = "okabe_ito", base_size = 11) {
  df <- raters_df[is.finite(raters_df$Estimate) & is.finite(raters_df$InfitMnSq), , drop = FALSE]
  if (!"DiagnosticFlag" %in% names(df)) {
    df$DiagnosticFlag <- ifelse(
      df$InfitMnSq > mnsq_range[2L], "Misfit (Unpredictable)",
      ifelse(df$InfitMnSq < mnsq_range[1L], "Muted / Centrality", "Acceptable")
    )
  }
  df$RealLow <- df$Estimate - 1.96 * df$RealSE
  df$RealHigh <- df$Estimate + 1.96 * df$RealSE

  unavailable <- df$DiagnosticFlag == "Review required"
  df$RealLow[unavailable] <- df$RealHigh[unavailable] <- NA_real_

  st <- .tidymfrm_flag_styles(palette = palette)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$InfitMnSq, y = .data$Estimate)) +
    ggplot2::geom_vline(xintercept = 1.0, linetype = "solid", colour = "#999999", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = mnsq_range, linetype = "dashed", colour = "#666666", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.0, linetype = "dotted", colour = "#888888", linewidth = 0.35) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$RealLow, ymax = .data$RealHigh,
                   colour = .data$DiagnosticFlag, linetype = .data$DiagnosticFlag),
      width = 0.0, linewidth = 0.55, alpha = 0.8, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$DiagnosticFlag, shape = .data$DiagnosticFlag),
      size = 3.0
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Level),
      nudge_y = 0.04 * max(diff(range(df$Estimate, na.rm = TRUE)), 0.5),
      size = 3.0, colour = "#222222", check_overlap = TRUE
    ) +
    ggplot2::scale_colour_manual(values = st$colours, name = "Diagnostic Status") +
    ggplot2::scale_shape_manual(values = st$shapes, name = "Diagnostic Status") +
    ggplot2::scale_linetype_manual(values = st$linetypes, name = "Diagnostic Status") +
    ggplot2::labs(
      x = "Infit Mean-Square (MnSq)",
      y = "Rater Severity (Logits; 95% Real-SE CI)"
    ) +
    .tidymfrm_theme(base_size = base_size)

  attr(p, "mfrm_notes") <- c(
    sprintf("Dashed vertical lines indicate Infit MnSq screening bounds [%.2f, %.2f].", mnsq_range[1L], mnsq_range[2L]),
    "Status is redundantly encoded via colour, point shape, and error-bar linetype.",
    unique(df$InferenceReason[nzchar(df$InferenceReason)])
  )
  p
}

.tidymfrm_plot_dstudy <- function(gt, palette = "okabe_ito", base_size = 11) {
  ds <- as.data.frame(gt$d_study)
  facet_cols <- gt$facets
  fixed <- character()
  primary_f <- facet_cols[1L]
  secondary_f <- if (length(facet_cols) >= 2L) facet_cols[2L] else NULL
  if (length(facet_cols) >= 3L) {
    for (fc3 in facet_cols[-c(1L, 2L)]) {
      u3 <- sort(unique(ds[[fc3]]))
      count <- u3[ceiling(length(u3) / 2)]
      ds <- ds[ds[[fc3]] == count, , drop = FALSE]
      fixed <- c(fixed, paste(fc3, "=", count))
    }
  }

  grp_chr <- if (!is.null(secondary_f)) paste0(secondary_f, " = ", ds[[secondary_f]]) else "Single-facet design"
  long_df <- rbind(
    data.frame(PrimaryN = ds[[primary_f]], SecondaryGroup = grp_chr,
               Metric = "G (Relative)", Coefficient = ds$G, stringsAsFactors = FALSE),
    data.frame(PrimaryN = ds[[primary_f]], SecondaryGroup = grp_chr,
               Metric = "Phi (Absolute)", Coefficient = ds$Phi, stringsAsFactors = FALSE)
  )
  grp_lvls <- unique(long_df$SecondaryGroup)
  long_df$SecondaryGroup <- factor(long_df$SecondaryGroup, levels = grp_lvls)
  cols <- stats::setNames(.tidymfrm_palette(length(grp_lvls), palette = palette), grp_lvls)
  shapes <- stats::setNames(.tidymfrm_shapes(length(grp_lvls)), grp_lvls)
  targets <- unname(gt$targets[c("G", "Phi")])

  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x = .data$PrimaryN,
      y = .data$Coefficient,
      colour = .data$SecondaryGroup,
      shape = .data$SecondaryGroup,
      linetype = .data$Metric
    )
  ) +
    ggplot2::geom_hline(yintercept = targets[1L], linetype = "dotted", colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.4, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = cols, name = if (!is.null(secondary_f)) secondary_f else "Design") +
    ggplot2::scale_shape_manual(values = shapes, name = if (!is.null(secondary_f)) secondary_f else "Design") +
    ggplot2::scale_linetype_manual(values = c("G (Relative)" = "solid", "Phi (Absolute)" = "dashed"), name = "Coefficient") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(
      x = sprintf("Sample Size for %s", primary_f),
      y = "Generalizability / Dependability"
    ) +
    .tidymfrm_theme(base_size = base_size)

  if (targets[1L] != targets[2L]) p <- p + ggplot2::geom_hline(
    yintercept = targets[2L], linetype = "dotdash", colour = "#666666", linewidth = 0.45)
  attr(p, "mfrm_notes") <- c(
    "G-Study and D-Study coefficients operate on the observed score scale (Brennan, 2001).",
    if (targets[1L] == targets[2L]) sprintf("Dotted line marks the shared G and Phi target %.2f.", targets[1L]) else
      sprintf("Dotted line marks G target %.2f. Dot-dash line marks Phi target %.2f.", targets[1L], targets[2L]),
    if (length(fixed)) paste("Displayed slice fixes", paste(fixed, collapse = ", ")),
    if (anyNA(ds$G) || anyNA(ds$Phi)) paste("Unavailable in this slice:",
      paste(unique(c(ds$GStatus[is.na(ds$G)], ds$PhiStatus[is.na(ds$Phi)])), collapse = ", "))
  )
  p
}

.tidymfrm_plot_context_curves <- function(x, newdata, palette, base_size) {
  facets <- x$spec$facets
  if (is.null(newdata)) newdata <- x$fit$encoding$retained_data[1L, facets, drop = FALSE]
  if (!is.data.frame(newdata) || !nrow(newdata) || !all(facets %in% names(newdata))) {
    stop("Provide one row per rating context in newdata with every model facet column.", call. = FALSE)
  }
  contexts <- unique(newdata[facets])
  labels <- make.unique(vapply(seq_len(nrow(contexts)), function(i) {
    paste(paste(facets, unlist(contexts[i, ], use.names = FALSE)), collapse = " / ")
  }, character(1)))
  theta <- seq(-4, 4, length.out = 121L)
  index <- rep(seq_len(nrow(contexts)), each = length(theta))
  nd <- contexts[index, , drop = FALSE]
  nd$theta <- rep(theta, times = nrow(contexts))
  probs <- stats::predict(x$fit, newdata = nd, type = "probs")
  df <- data.frame(theta = nd$theta, ContextID = index,
                   RatingContext = factor(labels[index], levels = labels),
                   Category = rep(colnames(probs), each = nrow(probs)),
                   Probability = as.vector(probs))
  df$Category <- factor(df$Category, levels = colnames(probs))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$theta, y = .data$Probability, colour = .data$Category,
                                       linetype = .data$Category)) +
    ggplot2::geom_line(linewidth = .8) +
    ggplot2::facet_wrap(~ RatingContext, labeller = ggplot2::label_wrap_gen(width = 28)) +
    ggplot2::scale_colour_manual(values = .tidymfrm_palette(ncol(probs), palette)) +
    ggplot2::scale_linetype_manual(values = .tidymfrm_linetypes(ncol(probs))) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Person ability (logits)", y = "Category probability") +
    .tidymfrm_theme(base_size)
  check <- summary(x$fit)$readiness
  attr(p, "mfrm_notes") <- c(
    "Curves use fitted steps, facet locations and all discrimination families for the stated rating contexts.",
    "These are conditional point predictions. Contexts can differ in severity, thresholds and discrimination.",
    if (!check$ready) paste("Review calibration:", paste(check$reason_codes, collapse = ", "))
  )
  p
}

















































plot.tidymfrm <- function(x,
                          type = c("wright", "raters", "steps", "category_curves", "dstudy", "persons",
                                   "person_fit", "residual_qq", "slopes", "drift", "qc", "dashboard",
                                   "dstudy_frontier", "gstudy", "pathway"),
                          palette = c("okabe_ito", "viridis", "monochrome"),
                          base_size = 11,
                          draw = TRUE,
                          ..., newdata = NULL) {
  x <- .tidymfrm_refresh(x)
  type <- match.arg(type)
  palette <- match.arg(palette)

  if (type == "dashboard") {
    if (is.null(x$fit)) stop("Run `step_mfrm_measure()` before plotting a dashboard.", call. = FALSE)
    p1 <- plot.tidymfrm(x, type = "wright", palette = palette, base_size = base_size * 0.9, draw = FALSE)
    p2 <- plot.tidymfrm(x, type = "raters", palette = palette, base_size = base_size * 0.9, draw = FALSE)
    p3 <- plot.tidymfrm(x, type = "steps", palette = palette, base_size = base_size * 0.9, draw = FALSE)
    p4 <- if (!is.null(x$gtheory)) {
      plot.tidymfrm(x, type = "dstudy", palette = palette, base_size = base_size * 0.9, draw = FALSE)
    } else {
      plot.tidymfrm(x, type = "persons", palette = palette, base_size = base_size * 0.9, draw = FALSE)
    }
    dash <- structure(
      list(wright = p1, raters = p2, steps = p3, fourth = p4),
      class = "tidymfrm_dashboard"
    )
    if (isTRUE(draw)) print(dash)
    return(invisible(dash))
  }

  p <- switch(
    type,
    pathway = {
      if (is.null(x$fit)) stop("Run step_mfrm_measure() before plotting a pathway.", call. = FALSE)
      graphics::plot(x$fit, type = "pathway", palette = palette, base_size = base_size, ...)
    },
    wright = {
      if (is.null(x$fit)) stop("Run `step_mfrm_measure()` before `type = 'wright'`.", call. = FALSE)
      graphics::plot(x$fit, type = "wright", palette = palette, base_size = base_size, ...)
    },
    raters = {
      if (!is.null(x$diagnostics$raters)) {
        .tidymfrm_plot_raters(x$diagnostics$raters, mnsq_range = x$diagnostics$mnsq_range,
                              palette = palette, base_size = base_size)
      } else if (!is.null(x$fit)) {
        graphics::plot(x$fit, type = "raters", facet = x$spec$rater_facet,
                       palette = palette, base_size = base_size, ...)
      } else {
        stop("Run `step_mfrm_measure()` before `type = 'raters'`.", call. = FALSE)
      }
    },
    steps = ,
    category_curves = {
      if (is.null(x$fit)) stop("Run `step_mfrm_measure()` before `type = 'steps'`.", call. = FALSE)
      if (identical(x$fit$model, "GPCM") || !is.null(newdata)) {
        .tidymfrm_plot_context_curves(x, newdata = newdata, palette = palette, base_size = base_size)
      } else graphics::plot(x$fit, type = "category_curves", palette = palette, base_size = base_size, ...)
    },
    person_fit = ,
    residual_qq = ,
    slopes = ,
    qc = {
      if (is.null(x$fit)) stop("Run `step_mfrm_measure()` before plotting.", call. = FALSE)
      graphics::plot(x$fit, type = type, palette = palette, base_size = base_size, ...)
    },
    drift = {
      if (is.null(x$drift)) stop("Supply `occasion` in `mfrm_analysis()` / `step_mfrm_measure()` before `type = 'drift'`.", call. = FALSE)
      graphics::plot(x$drift, palette = palette, base_size = base_size, ...)
    },
    dstudy = {
      if (is.null(x$gtheory)) stop("Run `step_mfrm_gtheory()` before `type = 'dstudy'`.", call. = FALSE)
      .tidymfrm_plot_dstudy(x$gtheory, palette = palette, base_size = base_size)
    },
    dstudy_frontier = {
      if (is.null(x$gtheory)) stop("Run step_mfrm_gtheory() before plotting D-Study plans.", call. = FALSE)
      .tidymfrm_plot_frontier(x$gtheory, palette, base_size)
    },
    gstudy = {
      if (is.null(x$gtheory)) stop("Run step_mfrm_gtheory() before plotting G-Study coefficients.", call. = FALSE)
      .tidymfrm_plot_gstudy(x$gtheory, palette, base_size)
    },
    persons = {
      if (is.null(x$fit)) stop("Run `step_mfrm_measure()` before `type = 'persons'`.", call. = FALSE)
      p_per <- graphics::plot(x$fit, type = "persons", palette = palette, base_size = base_size, ...)
      cuts <- if (!is.null(x$diagnostics$decisions)) x$diagnostics$decisions$CutFairM else NULL
      if (!is.null(cuts) && length(cuts) >= 1L) {
        p_per <- p_per + ggplot2::geom_hline(
          yintercept = as.numeric(cuts),
          linetype = "dashed", colour = "#555555", linewidth = 0.45
        )
      }
      p_per
    }
  )
  if (isTRUE(draw)) print(p)
  if (isTRUE(draw)) invisible(p) else p
}


print.tidymfrm_dashboard <- function(x, ...) {
  grid::grid.newpage()
  vp_layout <- grid::grid.layout(nrow = 2L, ncol = 2L)
  grid::pushViewport(grid::viewport(layout = vp_layout))
  positions <- list(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
  for (i in seq_along(x)) {
    pos <- positions[[i]]
    print(x[[i]], vp = grid::viewport(layout.pos.row = pos[1L], layout.pos.col = pos[2L]))
  }
  grid::popViewport()
  invisible(x)
}
