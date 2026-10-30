# Error-variance contributions, projected from the same components as the D-study.
.mfrm_mvds_error_plot <- function(x, score, kind, weights, metric, ncol, style, draw) {
  metrics <- c("RelativeErrorVariance", "AbsoluteErrorVariance")
  if (!is.null(metric)) {
    if (!is.character(metric) || length(metric) != 1L || is.na(metric) || !metric %in% metrics)
      stop("For `error_components`, `metric` must be RelativeErrorVariance or AbsoluteErrorVariance.", call. = FALSE)
    metrics <- metric
  }
  if (!is.null(ncol) && (!is.numeric(ncol) || is.complex(ncol) || length(ncol) != 1L ||
      !is.finite(ncol) || ncol < 1 || ncol != floor(ncol)))
    stop("`ncol` must be a positive integer.", call. = FALSE)
  # Reuse the point-projection validation (including old count maps); no fit or
  # interval calculation. The displayed components must explain the saved result.
  d <- .mfrm_mvds_replay(x)
  totals <- x$coefficients[x$coefficients$Kind == kind & x$coefficients$Score %in% score, , drop = FALSE]
  g <- d$gstudy
  scores <- g$design$scores
  vectors <- if (kind == "Score") {
    v <- diag(length(scores))[, match(score, scores), drop = FALSE]
    dimnames(v) <- list(scores, score)
    v
  } else if (is.matrix(weights)) weights[scores, score, drop = FALSE] else
    matrix(weights[scores], ncol = 1L, dimnames = list(scores, score))
  sources <- names(g$components)
  source_labels <- gsub(":", " x ", sources, fixed = TRUE)
  nesting <- g$design$nesting
  highest <- if (is.null(nesting)) paste(names(g$design$counts), collapse = " x ") else
    paste0("Person x ", names(nesting), "(", nesting, ")")
  source_labels[sources == "Residual"] <- paste(highest, "+ residual")
  count_labels <- names(d$design_grid)
  if (!is.null(nesting)) {
    child <- match(g$design$count_columns[[names(nesting)]], count_labels)
    count_labels[child] <- paste(count_labels[child], "per", unname(nesting))
  }
  plans <- vapply(seq_len(nrow(d$design_grid)), function(i)
    paste0("Plan ", i, ": ", paste(count_labels, unlist(d$design_grid[i, ]), sep = " = ", collapse = ", ")),
    character(1))
  panels <- expand.grid(Score = score, Metric = metrics, stringsAsFactors = FALSE)
  panels$Panel <- seq_len(nrow(panels))
  panels$Label <- paste(panels$Score, ifelse(panels$Metric == "RelativeErrorVariance",
    "Relative error", "Absolute error"), sep = ": ")
  rows <- vector("list", nrow(panels)*length(plans))
  for (j in seq_len(nrow(panels))) {
    w <- vectors[, panels$Score[j]]
    scale <- max(abs(w))
    w <- w/scale
    # Separate diagonal and cross-score terms without subtracting two nearly
    # equal quadratic forms. Scaling matches the D-study projection.
    diagonal <- vapply(d$effective_components, function(a) sum(w^2*diag(a)), numeric(1))*scale^2
    cross <- vapply(d$effective_components, function(a)
      2*sum((outer(w, w)*a)[upper.tri(a)]), numeric(1))*scale^2
    projected <- vapply(d$effective_components, function(a)
      as.numeric(crossprod(w, a %*% w)), numeric(1))*scale^2
    for (i in seq_along(plans)) {
      multiplier <- .mfrm_mvgt_profile_error(d, i,
        if (panels$Metric[j] == "RelativeErrorVariance") "RelativeSEM" else "AbsoluteSEM")
      use <- multiplier != 0
      z <- data.frame(Scenario = i, d$design_grid[i, , drop = FALSE], Kind = kind,
        Score = panels$Score[j], Metric = panels$Metric[j], Panel = j,
        Plan = plans[i], Source = sources[use], SourceLabel = source_labels[use],
        Multiplier = unname(multiplier[use]), ComponentVariance = unname(projected[use]),
        DiagonalContribution = unname(diagonal[use]*multiplier[use]),
        CovarianceContribution = unname(cross[use]*multiplier[use]),
        Contribution = unname(projected[use]*multiplier[use]), row.names = NULL)
      total <- totals[totals$Scenario == i & totals$Score == panels$Score[j], , drop = FALSE]
      target <- total[[panels$Metric[j]]]
      error_scale <- max(abs(target), sum(abs(z$Contribution)), .Machine$double.xmin)
      if (any(!is.finite(as.matrix(z[c("ComponentVariance", "DiagonalContribution",
          "CovarianceContribution", "Contribution")]))))
        stop("Component contributions exceed numerical range; rescale scores or weights.", call. = FALSE)
      if (abs(sum(z$Contribution)-target) > 256*.Machine$double.eps*error_scale)
        stop("Component contributions do not reproduce the D-study error variance within numerical precision.", call. = FALSE)
      z$TotalErrorVariance <- target
      z$SEMStatus <- total[[if (panels$Metric[j] == "RelativeErrorVariance")
        "RelativeSEMStatus" else "AbsoluteSEMStatus"]]
      rows[[(j-1L)*length(plans)+i]] <- z
    }
  }
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  panel_columns <- min(nrow(panels), if (is.null(ncol)) 1L else ncol)
  cols <- if (style$name == "monochrome") rep(style$foreground, length(plans)) else
    grDevices::hcl.colors(length(plans), "Dark 3")
  shapes <- rep(c(16, 17, 15, 18, 1, 2, 0, 5), length.out = length(plans))
  out <- new_mfrm_plot_data("multivariate_d_study", list(
    plot = "error_components", table = table, series = table, totals = totals,
    design_grid = d$design_grid, weights = weights, nesting = nesting,
    negative_components = d$negative_components,
    component_diagnostics = d$component_diagnostics,
    score = score, kind = kind, metric = metrics, panels = panels, ncol = panel_columns,
    source_labels = source_labels[sources %in% table$Source],
    point_shapes = shapes, preset = style$name,
    legend = new_plot_legend(plans, rep("plan", length(plans)), rep("color", length(plans)), cols),
    cautions = c("Signed point contributions to error variance (squared score units), not SEM or confidence intervals.",
      "Contributions sum to each saved error variance; covariance terms are included for composites. Weights are not normalized.",
      .mfrm_mvds_component_note(d),
      "Plans assume future complete common ratings; equal rating counts need not mean equal cost. This plot does not select optimal counts or weights."),
    title = "Error sources under planned measurement conditions",
    x_label = "Contribution to error variance (squared score units)"))
  if (draw) .mfrm_mvds_draw_errors(out$data, style)
  invisible(out)
}

.mfrm_mvds_draw_errors <- function(p, style) {
  old <- graphics::par()[c("mfrow", "mar", "oma", "cex", "mex", "fig", "mfg")]
  on.exit(graphics::par(old), add = TRUE)
  apply_plot_preset(style)
  labels <- unique(p$table$SourceLabel)
  graphics::par(mfrow = c(ceiling(nrow(p$panels)/p$ncol), p$ncol),
    mar = c(4.3, 4, 2.5, 1), oma = c(2+length(p$legend$label), 0, 0, 0))
  margin <- graphics::par("mar")
  margin[2] <- max(4, max(graphics::strwidth(labels, units = "inches", cex = .8))/
    (graphics::par("csi")*graphics::par("mex"))+2)
  graphics::par(mar = margin)
  limits <- range(c(0, p$table$Contribution))
  if (diff(limits) == 0) limits <- limits + c(-1, 1)*.04
  offsets <- if (nrow(p$legend) == 1L) 0 else seq(-.24, .24, length.out = nrow(p$legend))
  for (j in p$panels$Panel) {
    z <- p$table[p$table$Panel == j, ]
    sources <- unique(z$SourceLabel)
    y <- match(z$SourceLabel, rev(sources))
    group <- match(z$Plan, p$legend$label)
    y <- y + offsets[group]
    plot(limits, c(.5, length(sources)+.5), type = "n",
      xlab = p$x_label, ylab = "", yaxt = "n", main = p$panels$Label[j])
    axis(2, at = seq_along(sources), labels = rev(sources), las = 1, cex.axis = .8)
    abline(v = 0, col = style$grid)
    segments(0, y, z$Contribution, y, col = p$legend$value[group])
    points(z$Contribution, y, col = p$legend$value[group], pch = p$point_shapes[group])
  }
  graphics::par(oma = rep(0, 4), fig = c(0, 1, 0, 1), mar = rep(0, 4), new = TRUE)
  plot.new()
  legend("bottom", legend = p$legend$label, col = p$legend$value,
    pch = p$point_shapes, bty = "n", cex = .8, xpd = NA)
}
