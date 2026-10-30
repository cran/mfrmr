# Shared presentation of saved extended-model estimates. No fitting or scoring.
mfrm_extended_estimate_plot <- function(table, labels, name, default_title, xlab,
    note, draw, settings, style = c("interval", "precision", "distribution"),
    sort = c("input", "estimate", "uncertainty"), palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, reference = 0, text_scale = 1, point_size = 2.5,
    extra = list()) {
  style <- match.arg(style); sort <- match.arg(sort); palette <- match.arg(palette)
  for (key in c("draw", "show_title", "show_notes", "show_labels")) {
    value <- get(key)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("`", key, "` must be TRUE or FALSE.", call. = FALSE)
  }
  for (key in c("title", "caption")) {
    value <- get(key)
    if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value))) stop("`", key, "` must be NULL or one string.", call. = FALSE)
  }
  for (key in c("text_scale", "point_size")) {
    value <- get(key)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) stop("`", key, "` must be positive and finite.", call. = FALSE)
  }
  if (!is.null(reference) && (!is.numeric(reference) || length(reference) != 1L || !is.finite(reference))) stop("`reference` must be NULL or one finite number.", call. = FALSE)
  prior <- if ("Status" %in% names(table)) table$Status %in% "prior_only" else rep(FALSE, nrow(table))
  unavailable <- if ("Status" %in% names(table)) table$Status %in% "unavailable" else rep(FALSE, nrow(table))
  width <- table$Upper - table$Lower
  rows <- switch(sort, input = seq_len(nrow(table)),
    estimate = order(table$Estimate, na.last = TRUE), uncertainty = order(width, na.last = TRUE))
  valid <- is.finite(table$Estimate) & !unavailable
  reason <- ifelse(valid, "", "Estimate unavailable")
  if (style == "precision") {
    reason[valid & !is.finite(width)] <- "Finite interval unavailable"
    valid <- valid & is.finite(width)
  } else if (style == "distribution") {
    reason[prior] <- "Prior only"
    valid <- valid & !prior
  }
  colours <- if (palette == "mono") c("#222222", "#222222") else c("#0072B2", "#9C4700")
  display_data <- data.frame(Label = as.character(labels), Estimate = table$Estimate,
    Width = width, Shape = ifelse(prior, 1, 16), Colour = colours[ifelse(prior, 2, 1)],
    Included = valid, Reason = reason, stringsAsFactors = FALSE)
  used <- if (style == "distribution") table$Estimate[valid] else numeric()
  values <- base::sort(unique(used))
  distribution <- data.frame(Estimate = values,
    Cumulative = if (length(values)) findInterval(values, base::sort(used)) / length(used) else numeric())
  limits <- range(c((reference %||% 0) + c(-.1, .1), table$Estimate,
    if (style == "interval") c(table$Lower, table$Upper,
      if (isTRUE(settings$comparison)) c(table$OrdinaryLower, table$OrdinaryUpper))), finite = TRUE)
  padding <- max(diff(limits) * .08, .1); limits <- limits + c(-padding, padding)
  distribution_path <- if (length(values)) data.frame(
    Estimate = c(limits[1], rep(values, each = 2), limits[2]),
    Cumulative = c(0, as.vector(rbind(c(0, head(distribution$Cumulative, -1)), distribution$Cumulative)), 1)
  ) else data.frame(Estimate = numeric(), Cumulative = numeric())
  view_note <- switch(style, interval = "",
    precision = "Width = upper minus lower bound; not a fit statistic or reliability coefficient.",
    distribution = paste0("Empirical distribution of point estimates; not a latent population distribution.",
      if ("Status" %in% names(table)) " Prior-only scores excluded." else ""))
  count_note <- paste0(sprintf("%d of %d estimates displayed; ",sum(valid),nrow(table)),
    if ("Status" %in% names(table)) sprintf("%d prior only; ",sum(prior)) else "",
    sprintf("%d omitted from this view.",sum(!valid)))
  interpretation <- paste(c(note, view_note[nzchar(view_note)], count_note[style != "interval" || any(!valid)]), collapse = "\n")
  caption_note <- if (style == "distribution" && grepl("Approximate|conditional intervals|normal approximation", note))
    "Point estimates only; no intervals displayed." else note
  annotation <- paste(c(caption_note, view_note[nzchar(view_note)], count_note[style != "interval" || any(!valid)]), collapse = "\n")
  payload <- c(list(table = table, settings = settings, labels = labels,
    title = title %||% default_title, xlab = xlab,
    ylab = switch(style, interval = "", precision = "Interval width (logits)", distribution = "Cumulative proportion"),
    caption = caption %||% annotation,
    notes = data.frame(Type = "Interpretation", Text = interpretation),
    alt_text = paste(default_title, paste0("(", style, " view)."), count_note, note, view_note),
    display = list(style = style, sort = sort, rows = rows, palette = palette,
      show_title = show_title, show_notes = show_notes, show_labels = show_labels,
      reference = reference, text_scale = text_scale, point_size = point_size,
      limits = limits, padding = padding), display_data = display_data,
    distribution = distribution, distribution_path = distribution_path), extra)
  out <- new_mfrm_plot_data(name, payload)
  if (draw) mfrm_draw_extended_estimates(payload)
  invisible(out)
}

mfrm_draw_extended_estimates <- function(payload) {
  opt <- payload$display; style <- opt$style
  tab <- payload$table[opt$rows, , drop = FALSE]
  dd <- payload$display_data[opt$rows, , drop = FALSE]
  colour <- if (opt$palette == "mono") "#222222" else "#0072B2"
  caption <- if (opt$show_notes) payload$caption else ""
  caption_lines <- unlist(lapply(strsplit(caption, "\n", fixed = TRUE)[[1]], strwrap, width = 85))
  left <- if (style == "interval" && opt$show_labels) max(5, min(18, max(nchar(dd$Label), 0) * .55 + 1)) else 5
  title <- if (opt$show_title) payload$title else ""
  title_lines <- unlist(lapply(strsplit(title, "\n", fixed = TRUE)[[1]], strwrap,
    width = max(20, floor(grDevices::dev.size("in")[1] * 7 / opt$text_scale))))
  old <- graphics::par(mar = c(5 + length(caption_lines) * .85, left,
    1.5 + length(title_lines) * 1.5, 2), cex = opt$text_scale)
  on.exit(graphics::par(old), add = TRUE)
  title <- paste(title_lines, collapse = "\n")
  y <- rev(seq_len(nrow(tab))); comparison <- isTRUE(payload$settings$comparison)
  offset <- if (comparison) .12 else 0
  ylim <- switch(style, interval = c(.5, max(1, nrow(tab)) + .5 + as.numeric(comparison)),
    precision = range(c(0, dd$Width[dd$Included]), finite = TRUE), distribution = c(0, 1))
  if (diff(ylim) == 0) ylim <- ylim + c(0, 1)
  if (style == "precision") ylim[2] <- ylim[2] * 1.15
  plot(NA_real_, NA_real_, type = "n", xlim = opt$limits, ylim = ylim,
    xlab = payload$xlab, ylab = payload$ylab, main = title,
    yaxt = if (style == "interval") "n" else "s", bty = "n")
  abline(h = if (style == "interval") y else graphics::axTicks(2), col = "grey92")
  if (!is.null(opt$reference)) abline(v = opt$reference, lty = 2, col = "grey45")
  if (style == "interval") {
    if (comparison) segments(tab$OrdinaryLower, y - offset, tab$OrdinaryUpper, y - offset,
      col = "grey35", lty = 2, lwd = 1.5)
    segments(pmax(opt$limits[1], tab$Lower), y + offset,
      pmin(opt$limits[2], tab$Upper), y + offset, col = dd$Colour, lwd = 1.5)
    for (i in seq_len(nrow(tab))) {
      if (is.infinite(tab$Lower[i])) arrows(opt$limits[1] + opt$padding, y[i] + offset,
        opt$limits[1], y[i] + offset, length = .08, col = dd$Colour[i])
      if (is.infinite(tab$Upper[i])) arrows(opt$limits[2] - opt$padding, y[i] + offset,
        opt$limits[2], y[i] + offset, length = .08, col = dd$Colour[i])
    }
    points(tab$Estimate[dd$Included], (y + offset)[dd$Included],
      pch = dd$Shape[dd$Included], col = dd$Colour[dd$Included], cex = opt$point_size / 2.5)
    if (opt$show_labels) axis(2, at = y, labels = dd$Label, las = 1, tick = FALSE)
    if (comparison) legend("top", horiz = TRUE, bty = "n", cex = .8,
      legend = c(paste("Bootstrap", payload$settings$method), "Ordinary normal"),
      col = c(colour, "grey35"), lty = c(1, 2))
  } else if (style == "precision") {
    keep <- dd$Included
    points(dd$Estimate[keep], dd$Width[keep], pch = dd$Shape[keep],
      col = dd$Colour[keep], cex = opt$point_size / 2.5)
    if (opt$show_labels && any(keep)) text(dd$Estimate[keep], dd$Width[keep], dd$Label[keep], pos = 3, cex = .8)
  } else if (nrow(payload$distribution)) {
    path <- payload$distribution_path; ecdf <- payload$distribution
    graphics::polygon(c(path$Estimate, rev(range(path$Estimate))),
      c(path$Cumulative, 0, 0), col = grDevices::adjustcolor(colour, alpha.f = .12), border = NA)
    lines(path, col = colour, lwd = 2)
    points(ecdf, col = colour, pch = 16, cex = opt$point_size / 2.5)
    graphics::rug(dd$Estimate[dd$Included], col = colour)
  }
  if (style != "interval" && !any(dd$Included)) text(mean(opt$limits), mean(ylim), "No eligible estimates")
  for (i in seq_along(caption_lines)) mtext(caption_lines[i], side = 1,
    line = 4 + (i - 1) * .85, cex = .7, adj = 0)
  invisible(payload)
}

.mfrmr_gg_extended_estimates <- function(payload) {
  opt <- payload$display; style <- opt$style
  df <- payload$table[opt$rows, , drop = FALSE]
  dd <- payload$display_data[opt$rows, , drop = FALSE]
  df$.Row <- rev(seq_len(nrow(df))); df$.Shape <- dd$Shape; df$.Colour <- dd$Colour
  df$.Label <- dd$Label; df$.Width <- dd$Width
  colour <- if (opt$palette == "mono") "#222222" else "#0072B2"
  comparison <- isTRUE(payload$settings$comparison); offset <- if (comparison) .12 else 0
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Estimate, y = .data$.Row))
  if (!is.null(opt$reference)) p <- p + ggplot2::geom_vline(xintercept = opt$reference, linetype = "dashed", colour = "grey45")
  if (style == "interval") {
    df$.Row <- df$.Row + offset
    ci <- df[!is.na(df$Lower) & !is.na(df$Upper), , drop = FALSE]
    ci$.Lower <- pmax(opt$limits[1], ci$Lower); ci$.Upper <- pmin(opt$limits[2], ci$Upper)
    ci$.Method <- rep(paste("Bootstrap", payload$settings$method), nrow(ci))
    df$.Method <- rep("Ordinary normal", nrow(df))
    interval_aes <- ggplot2::aes(x = .data$.Lower, xend = .data$.Upper, yend = .data$.Row, colour = .data$.Colour)
    if (comparison) interval_aes$linetype <- ggplot2::aes(linetype = .data$.Method)$linetype
    p <- p + ggplot2::geom_segment(data = ci, interval_aes, linewidth = .6) +
      ggplot2::geom_point(data = df[dd$Included, , drop = FALSE],
        ggplot2::aes(shape = .data$.Shape, colour = .data$.Colour), size = opt$point_size)
    if (comparison) {
      # The legend uses line types so comparison remains readable in monochrome.
      p <- p + ggplot2::geom_segment(data = df[is.finite(df$OrdinaryLower) & is.finite(df$OrdinaryUpper), , drop = FALSE],
        ggplot2::aes(x = .data$OrdinaryLower, xend = .data$OrdinaryUpper,
          y = .data$.Row - 2 * offset, yend = .data$.Row - 2 * offset, linetype = .data$.Method), colour = "grey35", linewidth = .6) +
        ggplot2::scale_linetype_manual(name = NULL,
          values = stats::setNames(c("solid", "dashed"), c(paste("Bootstrap", payload$settings$method), "Ordinary normal")))
    }
    for (side in c("Lower", "Upper")) {
      arrow <- ci[is.infinite(ci[[side]]), , drop = FALSE]
      if (!nrow(arrow)) next
      direction <- if (side == "Lower") 1 else -1
      arrow$.End <- opt$limits[if (side == "Lower") 1 else 2]
      arrow$.Start <- arrow$.End + direction * opt$padding
      p <- p + ggplot2::geom_segment(data = arrow,
        ggplot2::aes(x = .data$.Start, xend = .data$.End, yend = .data$.Row, colour = .data$.Colour),
        arrow = ggplot2::arrow(length = ggplot2::unit(.08, "inches")), linewidth = .6)
    }
    p <- p + ggplot2::scale_y_continuous(breaks = rev(seq_len(nrow(df))),
      labels = if (opt$show_labels) dd$Label else rep("", nrow(df)), limits = c(.5, max(1, nrow(df)) + .5))
  } else if (style == "precision") {
    points <- df[dd$Included, , drop = FALSE]
    p <- p + ggplot2::geom_point(data = points,
      ggplot2::aes(y = .data$.Width, shape = .data$.Shape, colour = .data$.Colour), size = opt$point_size) +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(.03, .15)))
    if (opt$show_labels) p <- p + ggplot2::geom_text(data = points,
      ggplot2::aes(y = .data$.Width, label = .data$.Label), vjust = -1, size = 3 * opt$text_scale)
  } else {
    path <- payload$distribution_path
    if (nrow(path)) {
      polygon <- rbind(path, data.frame(Estimate = rev(range(path$Estimate)), Cumulative = c(0, 0)))
      p <- p + ggplot2::geom_polygon(data = polygon,
        ggplot2::aes(x = .data$Estimate, y = .data$Cumulative), fill = colour, alpha = .12) +
        ggplot2::geom_path(data = path, ggplot2::aes(y = .data$Cumulative), colour = colour, linewidth = .8) +
        ggplot2::geom_point(data = payload$distribution, ggplot2::aes(y = .data$Cumulative), colour = colour, size = opt$point_size) +
        ggplot2::geom_rug(data = df[dd$Included, , drop = FALSE], ggplot2::aes(x = .data$Estimate),
          inherit.aes = FALSE, sides = "b", colour = colour)
    }
    p <- p + ggplot2::scale_y_continuous(limits = c(0, 1))
  }
  if (style != "interval" && !any(dd$Included)) p <- p + ggplot2::annotate("text",
    x = mean(opt$limits), y = .5, label = "No eligible estimates")
  p <- p + ggplot2::scale_shape_identity() + ggplot2::scale_colour_identity() +
    ggplot2::coord_cartesian(xlim = opt$limits) + .mfrmr_gg_theme() +
    ggplot2::theme(text = ggplot2::element_text(size = 11 * opt$text_scale),
      legend.position = if (comparison) "bottom" else "none")
  p <- .mfrmr_gg_labs(p, payload, x = payload$xlab, y = payload$ylab)
  p <- p + ggplot2::labs(alt = payload$alt_text)
  attr(p, "mfrmr_alt_text") <- payload$alt_text
  attr(p, "mfrmr_display_data") <- payload$display_data
  p
}
