#' Compare fit-screening thresholds against known simulation truth
#'
#' See how detection and false-flag rates change when you change fit-screening
#' thresholds in a simulation with known truth. Apply specified mean-square
#' bands to saved Infit/Outfit values without refitting models. Keep underfit,
#' overfit and their union separate. Here sensitivity analysis means comparing
#' thresholds; the reported Sensitivity rate means detecting an affected target.
#' @inheritParams mfrm_screening_performance
#' @param measures Data frame with `Condition`, `Replicate`, `Target`, `Infit`
#'   and `Outfit`. Include `InfitZSTD` and `OutfitZSTD` when requesting the
#'   combined ZSTD rule. The selected statistic's columns are required;
#'   nonfinite values and omitted rows remain unavailable. Additional source
#'   columns, such as numerical readiness, are retained without filtering.
#' @param thresholds Data frame with unique nonempty `Profile` labels and
#'   numeric `Lower`, `Upper` bounds satisfying `0 < Lower < Upper`.
#'   Profiles are compared in supplied order, without selecting an optimum.
#' @param statistic `"either"` (default) flags Infit OR Outfit; `"infit"` or
#'   `"outfit"` uses only the specified index.
#' @param zstd_cut `NULL` (default) uses only mean squares. A positive finite
#'   number adds directional ZSTD flags with inclusive boundaries. This is an
#'   explicit alternative rule, not a significance calibration. Record the
#'   chosen df convention and residual definition in `rule`.
#' @details Mean-square boundaries are strict: equality is inside the band.
#'   Missing values use three-valued logic: TRUE OR NA is TRUE, but FALSE OR NA
#'   is NA. The family event counts any affected or any unaffected target once
#'   per independent replication. Correlated raters are not independent trials.
#'   The source roster fixes truth and planned denominators for every profile.
#'   `Affected` denotes the departure of interest; when comparing directional
#'   screens, an affected-target flag can have the wrong direction for that
#'   departure. Inspect directions together, not just the union's sensitivity.
#'
#'   Thresholds are review heuristics, not universal error-controlled tests.
#'   A low mean square describes low residual variability, not poor rater quality
#'   or proof of misconduct. ZSTD is sensitive to sample size and its df
#'   convention. Selecting a threshold on these results and reporting its same-
#'   sample performance is optimistic; use a separately designed validation.
#'   Pointwise Monte Carlo intervals are not simultaneous across profiles or
#'   conditions. Their overlap is not a paired test of rules applied to the
#'   same replications. No automatic exclusion or refitting is performed.
#'
#'   This helper does not qualify the supplied statistics. In particular,
#'   ordinary Rasch bands must not be transferred to extended-model posterior-
#'   predictive residual summaries that lack the same reference distribution.
#' @return An `mfrm_screening_sensitivity` object with `by_target`, `by_family`,
#'   aligned `outcomes`, replication-level `family_outcomes`, source `roster`
#'   and `measures`, `thresholds` and `settings`. Each result table identifies
#'   `Profile`, `Lower`, `Upper` and `Direction`; the remaining columns have
#'   the meanings documented by [mfrm_screening_performance()].
#' @seealso [fit_measures_table()], [mfrm_screening_performance()]
#' @references Linacre, J. M. (2003). Size vs. significance: standardized
#'   chi-square fit statistic. *Rasch Measurement Transactions*, 17(1), 918.
#'   \url{https://www.rasch.org/rmt/rmt171n.htm}.
#' @examples
#' roster <- expand.grid(Condition = "Null", Replicate = 1:4,
#'                       Target = c("R1", "R2"), stringsAsFactors = FALSE)
#' roster$Affected <- FALSE
#' measures <- roster[c("Condition", "Replicate", "Target")]
#' measures$Infit <- c(.45, .8, 1.1, NA, .9, 1.4, 1.2, NA)
#' measures$Outfit <- measures$Infit
#' bands <- data.frame(Profile = c("Broad", "Narrow"),
#'                     Lower = c(.5, .7), Upper = c(1.5, 1.3))
#' sensitivity <- mfrm_screening_sensitivity(roster, measures, bands,
#'   rule = "Illustrative fixed-rater residual screen; no exclusion/refit")
#' plot(sensitivity, direction = "overfit")
#' @export
#' @param draw For `plot()`: `FALSE` returns plot data without opening a device.
#' @order 1
mfrm_screening_sensitivity <- function(roster, measures, thresholds, rule,
    statistic = c("either", "infit", "outfit"), zstd_cut = NULL, level = .95) {
  statistic <- match.arg(statistic)
  columns <- switch(statistic, either = c("Infit", "Outfit"), infit = "Infit", outfit = "Outfit")
  if (!is.null(zstd_cut) && (!is.numeric(zstd_cut) || is.complex(zstd_cut) ||
      length(zstd_cut) != 1L || !is.finite(zstd_cut) || zstd_cut <= 0)) stop("Use NULL or a positive finite zstd_cut.", call. = FALSE)
  required <- c("Condition", "Replicate", "Target", columns,
    if (!is.null(zstd_cut)) paste0(columns, "ZSTD"))
  if (!is.data.frame(measures) || anyDuplicated(names(measures)) || !all(required %in% names(measures))) stop(
    "`measures` needs identifiers and the selected fit-statistic columns.", call. = FALSE)
  for (name in setdiff(required, c("Condition", "Replicate", "Target"))) {
    v <- measures[[name]]
    if (!is.numeric(v) || is.complex(v) || !is.null(dim(v)) ||
        (name %in% columns && any(v[is.finite(v)] < 0))) stop("Fit columns must be numeric; finite mean squares must be nonnegative.", call. = FALSE)
  }
  if (!is.data.frame(thresholds) || anyDuplicated(names(thresholds)) ||
      !all(c("Profile", "Lower", "Upper") %in% names(thresholds)) || !nrow(thresholds)) stop("Supply Profile, Lower and Upper in thresholds.", call. = FALSE)
  if (!is.character(thresholds$Profile) || anyNA(thresholds$Profile) ||
      any(!nzchar(trimws(thresholds$Profile))) || anyDuplicated(thresholds$Profile)) stop("Profile labels must be distinct nonempty strings.", call. = FALSE)
  for (key in c("Lower", "Upper")) if (!is.numeric(thresholds[[key]]) || is.complex(thresholds[[key]]) ||
      !is.null(dim(thresholds[[key]])) || any(!is.finite(thresholds[[key]]))) stop("Threshold bounds must be finite numbers.", call. = FALSE)
  if (any(thresholds$Lower <= 0 | thresholds$Upper <= thresholds$Lower)) stop("Use 0 < Lower < Upper for each profile.", call. = FALSE)
  parts <- list()
  for (i in seq_len(nrow(thresholds))) {
    band <- thresholds[i, ]
    flags <- mfrm_fit_flags(measures, band$Lower, band$Upper,
      statistic = if (statistic == "either") "either" else columns, zstd_cut = zstd_cut)
    for (direction in c("underfit", "overfit", "either")) {
      results <- measures[c("Condition", "Replicate", "Target")]; results$Flag <- flags[[direction]]
      part <- mfrm_screening_performance(roster, results, rule, level)
      for (name in c("by_target", "by_family", "outcomes", "family_outcomes")) {
        tab <- part[[name]]
        tab$Profile <- rep(band$Profile, nrow(tab)); tab$Lower <- rep(band$Lower, nrow(tab)); tab$Upper <- rep(band$Upper, nrow(tab))
        tab$Direction <- rep(direction, nrow(tab)); part[[name]] <- tab
      }
      parts[[length(parts) + 1L]] <- part
    }
  }
  out <- lapply(c("by_target", "by_family", "outcomes", "family_outcomes"),
    function(name) do.call(rbind, lapply(parts, `[[`, name)))
  names(out) <- c("by_target", "by_family", "outcomes", "family_outcomes")
  out$roster <- roster; out$measures <- measures; out$thresholds <- thresholds
  out$settings <- list(rule = rule, statistic = statistic, zstd_cut = zstd_cut,
    flag_basis = if (is.null(zstd_cut)) "mnsq" else "mnsq_or_zstd", level = level,
    uncertainty_unit = "Independent simulation replication within condition")
  structure(out, class = "mfrm_screening_sensitivity")
}

#' @rdname mfrm_screening_sensitivity
#' @param object An `mfrm_screening_sensitivity` object.
#' @param ... For `summary.mfrm_screening_sensitivity()`: Unused by print and summary.
#'
#' For `plot()`: Unused.
#' @export
summary.mfrm_screening_sensitivity <- function(object, ...) {
  list(by_family = object$by_family, by_target = object$by_target, settings = object$settings)
}

#' @rdname mfrm_screening_sensitivity
#' @export
print.mfrm_screening_sensitivity <- function(x, ...) {
  cat("Fit-screening threshold sensitivity\n", x$settings$rule, "\n")
  print(x$by_family[x$by_family$Planned > 0, ], row.names = FALSE)
  invisible(x)
}

#' @rdname mfrm_screening_sensitivity
#' @details `plot()` method: Plot directional screening rates across threshold profiles
#' @param x Output from [mfrm_screening_sensitivity()].
#'
#' Also for `summary()`: An `mfrm_screening_sensitivity` object.
#' @param style `"tiles"` shows a labelled rate matrix. `"curves"` shows
#'   profile comparisons with pointwise Monte Carlo intervals by condition.
#' @param direction One or more of `"either"`, `"underfit"`, `"overfit"`.
#' @param metric `"false_flags"` or `"detection"` chooses the family event.
#' @param condition Optional condition labels to display.
#' @param quantity `"rate"` (default) or `"unavailable"`, the proportion of
#'   planned family events whose outcome is unresolved.
#' @inheritParams plot.mfrm_extended_comparison
#' @param palette `"accessible"` uses a sequential blue scale for tiles and
#'   blue/orange/purple lines with different symbols. `"mono"` uses greys.
#' @param show_labels Show percentages and known/planned counts in tiles.
#' @return For `plot()`: Invisibly, `mfrm_plot_data` with exact selected tables, settings,
#'   notes and alternative text. [as_ggplot()] supports either style.
#' @details For `plot()`: Tile shading is supplemented by numeric labels; missing estimates
#'   are labelled NA. Profile order is supplied order, not an optimized ranking.
#'   Curves connect discrete profile choices, not a continuous threshold scale.
#'   Tile plots do not display uncertainty intervals: consult the retained table
#'   or curve view. Exact intervals are pointwise and conditional on known
#'   outcomes; all-trial bounds are retained, not confidence intervals.
#' @export
plot.mfrm_screening_sensitivity <- function(x, style = c("tiles", "curves"),
    direction = "either", metric = c("false_flags", "detection"), condition = NULL,
    quantity = c("rate", "unavailable"), draw = TRUE, palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, text_scale = 1, point_size = 2.5, ...) {
  rlang::check_dots_empty(); style <- match.arg(style); metric <- match.arg(metric)
  quantity <- match.arg(quantity); palette <- match.arg(palette)
  if (!is.character(direction) || !length(direction) || anyNA(direction) || anyDuplicated(direction) ||
      !all(direction %in% c("underfit", "overfit", "either"))) stop("Choose distinct underfit, overfit or either directions.", call. = FALSE)
  for (key in c("draw", "show_title", "show_notes", "show_labels")) if (!is.logical(get(key)) || length(get(key)) != 1L || is.na(get(key))) stop("Display switches must be TRUE or FALSE.", call. = FALSE)
  for (key in c("title", "caption")) if (!is.null(get(key)) && (!is.character(get(key)) || length(get(key)) != 1L || is.na(get(key)))) stop("Titles/captions must be NULL or one string.", call. = FALSE)
  for (key in c("text_scale", "point_size")) if (!is.numeric(get(key)) || length(get(key)) != 1L || !is.finite(get(key)) || get(key) <= 0) stop("Text and point sizes must be positive finite numbers.", call. = FALSE)
  conditions <- unique(x$by_family$Condition); condition <- condition %||% conditions
  if (!is.character(condition) || !length(condition) || anyNA(condition) || anyDuplicated(condition) || !all(condition %in% conditions)) stop("Choose recorded condition labels.", call. = FALSE)
  event <- if (metric == "false_flags") "Any unaffected target flagged" else "Any affected target flagged"
  tab <- x$by_family[x$by_family$Metric == event & x$by_family$Direction %in% direction &
    x$by_family$Condition %in% condition & x$by_family$Planned > 0, ]
  if (!nrow(tab)) stop("No planned targets for this event.", call. = FALSE)
  condition <- condition[condition %in% tab$Condition]
  tab$Profile <- factor(tab$Profile, levels = x$thresholds$Profile)
  tab$Condition <- factor(tab$Condition, levels = condition)
  tab$Direction <- factor(tab$Direction, levels = direction)
  tab$Value <- if (quantity == "rate") tab$Rate else tab$Unavailable / tab$Planned
  tab$Label <- vapply(tab$Value, function(v) if (is.finite(v))
    paste0(format(signif(100 * v, 3), trim = TRUE, scientific = FALSE), "%") else "NA", character(1))
  tab$Label <- paste0(tab$Label, "\n", tab$Available, "/", tab$Planned)
  note <- if (quantity == "rate") "Conditional rates; known/planned counts. Intervals are pointwise; unresolved-outcome bounds remain in the table." else
    "Unresolved family events / planned replications; a known positive remains known when another target is unavailable."
  payload <- list(table = tab, thresholds = x$thresholds, settings = x$settings,
    title = title %||% paste(event, "across thresholds"), caption = caption %||% note,
    alt_text = paste(event, "by declared threshold profile and condition.", paste(direction, collapse = ", "), note),
    notes = data.frame(Note = c(note, "No best threshold or automatic rater exclusion is selected.")),
    display = list(style = style, quantity = quantity, palette = palette, show_title = show_title,
      show_notes = show_notes, show_labels = show_labels, text_scale = text_scale, point_size = point_size))
  out <- new_mfrm_plot_data("screening_sensitivity", payload)
  if (draw) mfrm_draw_screening_sensitivity(payload)
  invisible(out)
}

mfrm_draw_screening_sensitivity <- function(payload) {
  opt <- payload$display; tab <- payload$table; profiles <- levels(tab$Profile)
  directions <- levels(tab$Direction); conditions <- levels(tab$Condition)
  tiles <- opt$style == "tiles"; panels <- if (tiles) directions else conditions
  old <- graphics::par(mfrow = c(ceiling(length(panels)/min(3,length(panels))), min(3,length(panels))),
    mar = if (tiles) c(5, 11, 3, 1) else c(5, 4, 3, 1),
    oma = c(if (opt$show_notes) 5 else 0, 0, if (opt$show_title) 3 else 0, 0), cex = opt$text_scale)
  on.exit(graphics::par(old), add = TRUE)
  colour <- if (opt$palette == "mono") c("#111111", "#666666", "#999999") else c("#0072B2", "#D55E00", "#7B3294")
  for (panel in panels) {
    d <- if (tiles) tab[tab$Direction == panel, ] else tab[tab$Condition == panel, ]
    if (tiles) {
      graphics::plot(NA_real_, NA_real_, xlim = c(.5,length(profiles)+.5), ylim = c(.5,length(conditions)+.5), axes = FALSE, xlab = "Threshold profile", ylab = "", main = panel)
      xx <- as.integer(d$Profile); yy <- length(conditions)+1-as.integer(d$Condition)
      colours <- mfrm_screening_tile_colours(d$Value, opt$palette)
      graphics::rect(xx-.48, yy-.48, xx+.48, yy+.48, col = colours$fill, border = "white")
      if (opt$show_labels) graphics::text(xx, yy, d$Label, cex = .7,
        col = colours$text)
      graphics::axis(1, seq_along(profiles), profiles, cex.axis = .7)
      graphics::axis(2, seq_along(conditions), rev(conditions), las = 1, tick = FALSE, cex.axis = .75)
    } else {
      graphics::plot(NA_real_, NA_real_, xlim = c(.8,length(profiles)+.2), ylim = c(0,1), xaxt = "n", xlab = "Threshold profile", ylab = "Proportion", main = panel, bty = "n")
      graphics::axis(1, seq_along(profiles), profiles, cex.axis = .7)
      for (j in seq_along(directions)) {
        a <- d[d$Direction == directions[j], ]; xx <- as.integer(a$Profile)
        if (opt$quantity == "rate") graphics::segments(xx,a$MCLower,xx,a$MCUpper,col=colour[j])
        graphics::lines(xx,a$Value,type="b",col=colour[j],pch=c(16,17,15)[j],lty=j,cex=opt$point_size/2.5)
      }
      graphics::legend("topleft", directions, col=colour[seq_along(directions)], pch=c(16,17,15)[seq_along(directions)],lty=seq_along(directions),bty="n",cex=.65)
    }
  }
  if (opt$show_title) graphics::mtext(payload$title,side=3,outer=TRUE,line=1)
  if (opt$show_notes) graphics::mtext(paste(strwrap(payload$caption,100),collapse="\n"),side=1,outer=TRUE,line=1,cex=.7)
}

mfrm_screening_tile_colours <- function(value, palette) {
  ramp <- grDevices::colorRamp(c("#F4F8FC", if (palette == "mono") "#222222" else "#084B83"), space = "Lab")
  fill <- grDevices::rgb(ramp(ifelse(is.finite(value), value, 0)), maxColorValue = 255)
  fill[!is.finite(value)] <- "grey85"
  rgb <- grDevices::col2rgb(fill)/255
  linear <- ifelse(rgb <= .04045, rgb/12.92, ((rgb + .055)/1.055)^2.4)
  luminance <- drop(c(.2126, .7152, .0722) %*% linear)
  list(fill = fill, text = ifelse(luminance <= sqrt(.0525)-.05, "white", "black"),
    contrast = pmax(1.05/(luminance+.05), (luminance+.05)/.05))
}

.mfrmr_gg_screening_sensitivity <- function(payload) {
  tab <- payload$table; opt <- payload$display
  p <- ggplot2::ggplot(tab)
  if (opt$style == "tiles") {
    p <- p + ggplot2::aes(x=.data$Profile,y=.data$Condition,fill=.data$Value) +
      ggplot2::geom_tile(colour="white") + ggplot2::facet_wrap(~Direction) +
      ggplot2::scale_y_discrete(limits=rev(levels(tab$Condition))) +
      ggplot2::scale_fill_gradient(low="#F4F8FC",high=if(opt$palette=="mono") "#222222" else "#084B83",limits=c(0,1),na.value="grey85",name="Proportion")
    if(opt$show_labels) p <- p + ggplot2::geom_text(ggplot2::aes(label=.data$Label),
      colour=mfrm_screening_tile_colours(tab$Value,opt$palette)$text,size=3*opt$text_scale)
  } else {
    p <- p + ggplot2::aes(x=.data$Profile,y=.data$Value,group=.data$Direction,colour=.data$Direction,shape=.data$Direction,linetype=.data$Direction)
    if(opt$quantity=="rate") p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$MCLower,ymax=.data$MCUpper),width=.08,na.rm=TRUE)
    p <- p + ggplot2::geom_line(na.rm=TRUE) + ggplot2::geom_point(size=opt$point_size,na.rm=TRUE) +
      ggplot2::facet_wrap(~Condition,ncol=min(3,nlevels(tab$Condition))) + ggplot2::scale_y_continuous(limits=c(0,1)) +
      ggplot2::scale_colour_manual(values=if(opt$palette=="mono") c("#111111","#666666","#999999") else c("#0072B2","#D55E00","#7B3294"))
  }
  p <- p + .mfrmr_gg_theme() + ggplot2::theme(text=ggplot2::element_text(size=11*opt$text_scale))
  p <- .mfrmr_gg_labs(p,payload,x="Threshold profile",y=if(opt$style=="tiles") NULL else "Proportion") + ggplot2::labs(alt=payload$alt_text)
  attr(p,"mfrmr_alt_text") <- payload$alt_text
  p
}
