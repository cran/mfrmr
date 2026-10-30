# ============================================================================
# Anchor Drift & Equating Chain Plots
# ============================================================================

# --- Internal plot helpers (not exported) ------------------------------------

.plot_drift_dot <- function(dt, config, draw = TRUE, style = resolve_plot_preset("standard"), ...) {
  out <- new_mfrm_plot_data(
    "anchor_drift",
    list(
      plot = "drift",
      table = dt,
      title = "Anchor drift",
      subtitle = paste0("Review threshold: |drift| > ", format(config$drift_threshold)),
      legend = new_plot_legend(
        label = c("Within review band", "Flagged drift"),
        role = c("status", "status"),
        aesthetic = c("point", "point"),
        value = c(style$accent_primary, style$warn)
      ),
      reference_lines = new_reference_lines(
        axis = c("v", "v", "v"),
        value = c(-config$drift_threshold, 0, config$drift_threshold),
        label = c("Drift review threshold", "Centered drift reference", "Drift review threshold"),
        linetype = c("dotted", "dashed", "dotted"),
        role = c("threshold", "reference", "threshold")
      ),
      preset = style$name
    )
  )
  if (!draw) return(invisible(out))

  opar <- graphics::par()["mar"]
  on.exit(graphics::par(opar), add = TRUE)

  dt <- dt |> dplyr::arrange(dplyr::desc(abs(.data$Drift)))
  if (!any(is.finite(dt$Drift))) {
    plot.new()
    text(.5, .5, "No finite drift estimates")
    return(invisible(out))
  }
  labels <- paste0(dt$Facet, ":", dt$Level, " @ ", dt$Wave)
  n <- nrow(dt)

  max_abs <- max(abs(dt$Drift), na.rm = TRUE) * 1.2
  if (max_abs == 0) max_abs <- 1

  label_width <- max(graphics::strwidth(labels, units = "inches", cex = .7))
  if (label_width + 1.5 >= graphics::par("fin")[1]) {
    stop("Anchor labels need more horizontal space; use a wider graphics device or inspect the table with `draw = FALSE`.",
      call. = FALSE)
  }
  margin <- label_width / (graphics::par("csi") * graphics::par("mex")) + 1.2
  graphics::par(mar = c(4, max(8, margin), 3, 1))
  plot(dt$Drift, seq_len(n), xlim = c(-max_abs, max_abs),
                 yaxt = "n", xlab = "Drift (logits)", ylab = "",
                 main = "", pch = 19,
                 col = ifelse(dt$Flag %in% TRUE, style$warn, style$accent_primary), ...)
  axis(2, at = seq_len(n), labels = labels, las = 1, cex.axis = 0.7)
  abline(v = pretty(c(-max_abs, max_abs), n = 5), col = grDevices::adjustcolor(style$grid, alpha.f = 0.85), lty = 1)
  abline(v = 0, lty = 2, col = grDevices::adjustcolor(style$foreground, alpha.f = 0.7))
  abline(v = c(-config$drift_threshold, config$drift_threshold),
                   lty = 3, col = grDevices::adjustcolor(style$warn, alpha.f = 0.9))

  invisible(out)
}

.plot_drift_heatmap <- function(dt, config, draw = TRUE, style = resolve_plot_preset("standard"), ...) {
  dt_wide <- dt |>
    dplyr::mutate(Element = paste0(.data$Facet, ":", .data$Level)) |>
    dplyr::select("Element", "Wave", "Drift")

  mat <- tryCatch({
    tidyr::pivot_wider(dt_wide, names_from = "Wave",
                        values_from = "Drift") |>
      tibble::column_to_rownames("Element") |>
      as.matrix()
  }, error = function(e) NULL)

  if (is.null(mat) || nrow(mat) == 0) {
    if (draw) message("Insufficient data for heatmap.")
    return(invisible(NULL))
  }

  out <- new_mfrm_plot_data(
    "anchor_drift",
    list(
      plot = "heatmap",
      table = dt,
      matrix = mat,
      title = "Anchor drift heatmap",
      subtitle = paste0("Wave-by-element drift; review threshold = ", format(config$drift_threshold)),
      legend = new_plot_legend(
        label = c("Negative drift", "Positive drift"),
        role = c("drift", "drift"),
        aesthetic = c("heatmap", "heatmap"),
        value = c(style$accent_primary, style$warn)
      ),
      reference_lines = new_reference_lines(),
      preset = style$name
    )
  )

  if (!draw) return(invisible(out))

  if (!any(is.finite(mat))) {
    plot.new()
    text(.5, .5, "No finite drift estimates")
    return(invisible(out))
  }
  max_abs <- max(abs(mat), na.rm = TRUE)
  if (max_abs == 0) max_abs <- 1
  n_colors <- 21
  breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)
  blues <- grDevices::colorRampPalette(c(style$accent_primary, "white", style$warn))(n_colors)

  opar <- graphics::par()["mar"]
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(5, 8, 3, 2))

  image(t(mat[nrow(mat):1, , drop = FALSE]),
                  axes = FALSE, col = blues, breaks = breaks,
                  main = "", ...)
  axis(1, at = seq(0, 1, length.out = ncol(mat)),
                 labels = colnames(mat), las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, 1, length.out = nrow(mat)),
                 labels = rev(rownames(mat)), las = 1, cex.axis = 0.7)

  invisible(out)
}

.equating_graph_data <- function(x) {
  links <- as.data.frame(x$links %||% data.frame(), stringsAsFactors = FALSE)
  detail <- as.data.frame(x$element_detail %||% data.frame(), stringsAsFactors = FALSE)
  if (!all(c("From", "To") %in% names(links))) {
    if (all(c("Wave1", "Wave2") %in% names(links))) {
      links$From <- links$Wave1
      links$To <- links$Wave2
    } else if (nrow(links)) {
      stop("Graph view requires From/To columns in links.", call. = FALSE)
    }
  }
  waves <- as.character(x$config$waves %||% x$cumulative$Wave %||%
                          unique(c(links$From, links$To)))
  if (!length(waves) || anyNA(waves) || any(!nzchar(waves))) {
    stop("Graph view requires non-missing, non-empty wave labels.", call. = FALSE)
  }
  links$LinkID <- as.character(links$Link %||% seq_len(nrow(links)))
  if (anyNA(links$LinkID) || anyDuplicated(links$LinkID)) {
    stop("Graph view requires unique, non-missing link IDs.", call. = FALSE)
  }
  if (!all(c("FromID", "ToID") %in% names(links))) {
    # Older chains have ordered adjacent links. Never split display labels.
    adjacent <- nrow(links) == length(waves) - 1L &&
      identical(as.character(links$From), head(waves, -1L)) &&
      identical(as.character(links$To), tail(waves, -1L))
    if (adjacent) {
      links$FromID <- seq_len(nrow(links))
      links$ToID <- seq_len(nrow(links)) + 1L
    } else {
      if (nrow(links) && anyDuplicated(waves)) {
        stop("Ambiguous wave labels; rebuild the chain with explicit wave IDs.", call. = FALSE)
      }
      links$FromID <- match(links$From, waves)
      links$ToID <- match(links$To, waves)
    }
  }
  for (endpoint in c("From", "To")) {
    id <- links[[paste0(endpoint, "ID")]]
    if (!is.numeric(id) || anyNA(id) || any(id != as.integer(id)) ||
        any(id < 1L | id > length(waves)) ||
        !identical(as.character(links[[endpoint]]), waves[id])) {
      if (nrow(links)) stop("Invalid wave IDs or endpoint labels in links.", call. = FALSE)
    }
  }
  if (nrow(detail)) {
    if (!all(c("Facet", "Level") %in% names(detail))) {
      stop("Graph view requires Facet/Level columns in element_detail.", call. = FALSE)
    }
    if ("LinkID" %in% names(detail)) {
      index <- match(as.character(detail$LinkID), links$LinkID)
    } else {
      keys <- paste(links$From, links$To, sep = " -> ")
      if (!"Link" %in% names(detail) || anyDuplicated(keys)) {
        stop("Ambiguous or missing element link IDs; rebuild the chain.", call. = FALSE)
      }
      index <- match(as.character(detail$Link), keys)
    }
    if (anyNA(index)) stop("Element links do not match the links table.", call. = FALSE)
  } else {
    detail$Facet <- detail$Level <- character()
    index <- integer()
  }
  detail$LinkID <- links$LinkID[index]
  detail$Facet <- as.character(detail$Facet)
  detail$Level <- as.character(detail$Level)
  if (anyNA(detail[, c("Facet", "Level")]) ||
      anyDuplicated(detail[, c("LinkID", "Facet", "Level")])) {
    stop("Graph elements require non-missing identities and unique rows per link.", call. = FALSE)
  }
  for (flag in c("Retained", "Flag")) {
    if (is.null(detail[[flag]])) detail[[flag]] <- rep(NA, nrow(detail))
    if (!is.logical(detail[[flag]])) {
      stop(sprintf("`element_detail$%s` must be logical.", flag), call. = FALSE)
    }
  }
  detail$Status <- ifelse(is.na(detail$Retained), "candidate",
    ifelse(!detail$Retained, "excluded",
      ifelse(is.na(detail$Flag) | detail$Flag, "retained_review", "retained")))
  detail$Reason <- ifelse(detail$Status == "candidate", "Retention not recorded",
    ifelse(detail$Status == "excluded", "Excluded from the retained set",
      ifelse(detail$Status == "retained_review", "Retained; drift flagged or flag unavailable",
        "Retained; no residual drift flag")))
  anchors <- unique(detail[, c("Facet", "Level"), drop = FALSE])
  anchors$AnchorId <- sprintf("anchor_%d", seq_len(nrow(anchors)))
  detail <- dplyr::left_join(detail, anchors, by = c("Facet", "Level"))
  # A node represents one element in one reviewed link. Merging the same
  # element across links would create paths that bypass pair-specific screening.
  detail$NodeId <- sprintf("element_%d", seq_len(nrow(detail)))
  wave_nodes <- data.frame(
    NodeId = paste0("wave_", seq_along(waves)), Label = waves,
    Kind = "wave", WaveID = seq_along(waves), stringsAsFactors = FALSE
  )
  element_nodes <- data.frame(
    NodeId = detail$NodeId,
    Label = sprintf("%s: %s [L%s]", detail$Facet, detail$Level, detail$LinkID),
    Kind = rep("element", nrow(detail)),
    detail[, c("AnchorId", "LinkID", "Facet", "Level", "Status"), drop = FALSE],
    stringsAsFactors = FALSE
  )
  nodes <- as.data.frame(dplyr::bind_rows(wave_nodes, element_nodes))
  edges <- dplyr::bind_rows(lapply(c("From", "To"), function(endpoint) {
    data.frame(
      from = sprintf("wave_%d", links[[paste0(endpoint, "ID")]][index]),
      to = detail$NodeId, Wave = waves[links[[paste0(endpoint, "ID")]][index]],
      detail[, c("AnchorId", "LinkID", "Facet", "Level", "Retained", "Flag", "Status", "Reason"), drop = FALSE],
      stringsAsFactors = FALSE
    )
  }))
  retained <- !is.na(edges$Retained) & edges$Retained
  nodes$RetainedDegree <- tabulate(match(c(edges$from[retained], edges$to[retained]), nodes$NodeId), nrow(nodes))
  links$from <- sprintf("wave_%d", links$FromID)
  links$to <- sprintf("wave_%d", links$ToID)
  notes <- .mfrmr_table_.equating_graph_data_notes
  isolated <- wave_nodes$Label[nodes$RetainedDegree[seq_along(waves)] == 0L]
  if (length(isolated)) notes <- rbind(notes, data.frame(Type = "warning", Text = paste0(
    "No retained common-element connection for: ", paste(isolated, collapse = ", "), ".")))
  weak <- links$LinkSupportAdequate %in% FALSE
  if (any(weak)) notes <- rbind(notes, data.frame(Type = "warning", Text = paste0(
    "Thin or absent linking support recorded for link(s): ", paste(links$LinkID[weak], collapse = ", "),
    ". Inspect the link and common_by_facet tables; retention alone does not imply adequate support.")))
  list(data = list(nodes = nodes, edges = edges, links = links, elements = detail,
    n_waves = length(waves), n_anchors = nrow(anchors), n_element_nodes = nrow(detail)), notes = notes)
}

.equating_link_topology <- function(graph_data, removal = FALSE) {
  data <- graph_data$data
  nodes <- data$nodes[data$nodes$Kind == "wave", c("NodeId", "Label", "WaveID"), drop = FALSE]
  links <- data$links
  elements <- data$elements
  li <- match(elements$LinkID, links$LinkID)
  retained <- elements$Retained %in% TRUE
  links$RetainedElements <- tabulate(li[retained], nrow(links))
  links$UnknownRetention <- tabulate(li[is.na(elements$Retained)], nrow(links))
  links$ReviewElements <- tabulate(li[retained & (is.na(elements$Flag) | elements$Flag)], nrow(links))
  links$Connected <- links$RetainedElements > 0L
  support <- links$LinkSupportAdequate %||% rep(NA, nrow(links))
  links$Status <- ifelse(!links$Connected, "no_retained",
    ifelse(links$ReviewElements > 0L | is.na(support) | !support,
      "retained_review", "retained"))
  membership <- function(connected) {
    uf <- make_union_find(nodes$NodeId)
    for (i in which(connected)) uf$union(links$from[i], links$to[i])
    root <- vapply(nodes$NodeId, uf$find, character(1))
    match(root, unique(root))
  }
  baseline <- membership(links$Connected)
  nodes$Component <- baseline
  connected_pairs <- function(component) {
    size <- as.numeric(table(component))
    sum(size * (size - 1) / 2)
  }
  baseline_pairs <- connected_pairs(baseline)
  result <- list(nodes = nodes, links = links, elements = elements,
    summary = data.frame(Waves = nrow(nodes), RecordedLinks = nrow(links),
      RetainedLinks = sum(links$Connected), Components = length(unique(baseline)),
      ConnectedPairs = baseline_pairs, UnknownRetention = sum(links$UnknownRetention)))
  notes <- rbind(graph_data$notes, data.frame(Type = "definition", Text =
    "A wave-to-wave edge exists when at least one common element is recorded as Retained = TRUE. Flagged retained elements are included; unknown retention is not counted. Component labels describe only this recorded graph."))
  if (any(is.na(elements$Retained))) notes <- rbind(notes, data.frame(Type = "warning", Text =
    "Retention is missing for some elements. Connectivity uses known retained rows only; absence of a recorded connection is not evidence of a disconnected rating design."))
  if (!removal) return(list(data = result, notes = notes))

  audit <- unique(elements[, c("AnchorId", "Facet", "Level"), drop = FALSE])
  audit$RetainedOccurrences <- audit$UnknownOccurrences <- audit$LostLinks <- audit$NewlyDisconnectedPairs <- integer(nrow(audit))
  audit$ComponentsAfter <- rep(length(unique(baseline)), nrow(audit))
  component_rows <- lost_link_rows <- lost_pair_rows <- list()
  # ponytail: one deletion per element; large wave sets may need lazy pair export.
  for (i in seq_len(nrow(audit))) {
    removed <- retained & elements$AnchorId == audit$AnchorId[i]
    remaining <- links$RetainedElements - tabulate(li[removed], nrow(links))
    lost <- links$Connected & remaining == 0L
    after <- if (any(lost)) membership(remaining > 0L) else baseline
    audit$RetainedOccurrences[i] <- sum(removed)
    audit$UnknownOccurrences[i] <- sum(is.na(elements$Retained) & elements$AnchorId == audit$AnchorId[i])
    audit$LostLinks[i] <- sum(lost)
    audit$ComponentsAfter[i] <- length(unique(after))
    audit$NewlyDisconnectedPairs[i] <- baseline_pairs - connected_pairs(after)
    component_rows[[i]] <- data.frame(AnchorId = audit$AnchorId[i],
      WaveID = nodes$WaveID, Wave = nodes$Label, Before = baseline, After = after)
    lost_link_rows[[i]] <- data.frame(AnchorId = rep(audit$AnchorId[i], sum(lost)),
      links[lost, c("LinkID", "FromID", "ToID", "From", "To"), drop = FALSE])
    # Enumerate only newly separated pairs, keeping existing disconnections out.
    pairs <- lapply(seq_len(nrow(nodes)), function(a) {
      b <- which(seq_len(nrow(nodes)) > a & baseline == baseline[a] & after != after[a])
      data.frame(AnchorId = rep(audit$AnchorId[i], length(b)),
        FromID = rep(nodes$WaveID[a], length(b)), ToID = nodes$WaveID[b],
        From = rep(nodes$Label[a], length(b)), To = nodes$Label[b])
    })
    lost_pair_rows[[i]] <- dplyr::bind_rows(pairs)
  }
  result$removal <- audit
  result$components_after <- if (length(component_rows)) dplyr::bind_rows(component_rows) else
    data.frame(AnchorId = character(), WaveID = integer(), Wave = character(), Before = integer(), After = integer())
  result$lost_links <- if (length(lost_link_rows)) dplyr::bind_rows(lost_link_rows) else
    data.frame(AnchorId = character(), LinkID = character(), FromID = integer(), ToID = integer(), From = character(), To = character())
  result$lost_pairs <- if (length(lost_pair_rows)) dplyr::bind_rows(lost_pair_rows) else
    data.frame(AnchorId = character(), FromID = integer(), ToID = integer(), From = character(), To = character())
  notes <- rbind(notes, data.frame(Type = "interpretation", Text =
    "Each scenario removes one facet/level identity from all recorded comparisons, starting from the unchanged baseline. Screening and other retained flags are held fixed. Newly disconnected pairs were connected at baseline, including indirect paths. Zero means no new graph disconnection, not negligible statistical influence. No offsets, estimates or SEs are recomputed."))
  list(data = result, notes = notes)
}

.equating_offset_sensitivity <- function(x) {
  source <- .equating_graph_data(x)$data
  links <- source$links
  elements <- source$elements
  waves <- source$nodes[source$nodes$Kind == "wave", c("WaveID", "Label"), drop = FALSE]
  threshold <- x$config$drift_threshold
  guideline <- x$config$min_common_per_facet
  if (!identical(x$config$method, "screened_common_element_alignment") ||
      !is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold < 0 ||
      !is.numeric(guideline) || length(guideline) != 1L || !is.finite(guideline) || guideline < 1 || guideline != floor(guideline)) {
    stop("Offset sensitivity requires the recorded screened-linking method, drift threshold and support guideline; rebuild the chain.", call. = FALSE)
  }
  if (nrow(waves) < 2L || nrow(links) != nrow(waves) - 1L ||
      !identical(as.integer(links$FromID), seq_len(nrow(links))) ||
      !identical(as.integer(links$ToID), seq_len(nrow(links)) + 1L)) {
    stop("Offset sensitivity requires an ordered adjacent chain.", call. = FALSE)
  }
  for (field in c("Est_From", "Est_To", "SE_From", "SE_To")) {
    if (!nrow(elements)) elements[[field]] <- numeric()
    if (!is.numeric(elements[[field]])) stop(
      "Offset sensitivity requires numeric Est_From/Est_To and SE_From/SE_To columns; rebuild the chain.", call. = FALSE)
  }
  if (!is.numeric(links$Offset)) stop("Offset sensitivity requires recorded numeric offsets.", call. = FALSE)
  by_link <- split(seq_len(nrow(elements)), factor(elements$LinkID, levels = links$LinkID))
  recalculate <- function(removed = NA_character_) {
    rows <- support <- detail <- vector("list", nrow(links))
    for (j in seq_len(nrow(links))) {
      original <- elements[by_link[[j]], , drop = FALSE]
      e <- original[is.na(removed) | original$AnchorId != removed, , drop = FALSE]
      diffs <- e$Est_To - e$Est_From
      info <- compute_equating_offset(diffs, e$SE_From, e$SE_To, threshold)
      s <- .summarise_link_support(e, info$retained, guideline)
      # Keep a baseline facet visible when deletion removes its last element.
      s <- dplyr::left_join(data.frame(Facet = unique(original$Facet)), s, by = "Facet")
      for (field in c("N_Common", "N_Retained")) s[[field]][is.na(s[[field]])] <- 0L
      s$GuidelineMinCommon <- rep(guideline, nrow(s))
      s$LinkSupportAdequate <- s$N_Retained >= guideline
      s$LinkID <- rep(links$LinkID[j], nrow(s))
      support[[j]] <- s
      status <- info$status
      rows[[j]] <- data.frame(LinkID = links$LinkID[j], FromID = links$FromID[j],
        ToID = links$ToID[j], From = links$From[j], To = links$To[j],
        N_Common = nrow(e), N_Retained = info$n_retained,
        N_Contributing = sum(info$contributing), Offset = if (status == "computed") info$offset else NA_real_,
        Offset_Prelim = info$offset_prelim, Offset_Method = info$weighting,
        ScreeningFallback = info$screening_fallback,
        LinkSupportAdequate = nrow(s) > 0L && all(s$LinkSupportAdequate),
        Status = status, Reason = info$reason)
      detail[[j]] <- data.frame(LinkID = e$LinkID, AnchorId = e$AnchorId,
        Facet = e$Facet, Level = e$Level, Est_From = e$Est_From, Est_To = e$Est_To,
        SE_From = e$SE_From, SE_To = e$SE_To, Retained = info$retained,
        Contributing = info$contributing, Flag = abs(info$residual) > threshold)
    }
    l <- dplyr::bind_rows(rows)
    list(links = l, support = dplyr::bind_rows(support), elements = dplyr::bind_rows(detail),
      cumulative = data.frame(WaveID = waves$WaveID, Wave = waves$Label,
        Offset = .linking_cumulative_offsets(l$Offset)))
  }
  baseline <- recalculate()
  same <- (is.na(links$Offset) & is.na(baseline$links$Offset)) |
    (is.finite(links$Offset) & is.finite(baseline$links$Offset) &
       abs(links$Offset - baseline$links$Offset) <= 1e-8 * pmax(1, abs(links$Offset)))
  if (anyNA(same) || !all(same)) stop(
    "Recomputed baseline does not match recorded offsets; rebuild the chain before sensitivity analysis.", call. = FALSE)
  audit <- unique(elements[, c("AnchorId", "Facet", "Level"), drop = FALSE])
  audit$MaxAbsLinkChange <- audit$MaxAbsCumulativeChange <- rep(NA_real_, nrow(audit))
  audit$UnavailableLinksAfter <- audit$NewlyUnavailableLinks <- audit$RetentionChanges <- audit$FallbackLinks <- integer(nrow(audit))
  audit$CompleteComparison <- rep(FALSE, nrow(audit))
  link_rows <- cumulative_rows <- changed_rows <- support_rows <- vector("list", nrow(audit))
  max_finite <- function(z) if (any(is.finite(z))) max(abs(z[is.finite(z)])) else NA_real_
  # ponytail: independent deletion reruns; retain full tables before adding batch optimizations.
  for (i in seq_len(nrow(audit))) {
    id <- audit$AnchorId[i]
    after <- recalculate(id)
    l <- after$links
    l$RemovedAnchorId <- rep(id, nrow(l))
    l$Offset_Before <- baseline$links$Offset
    l$Change <- l$Offset - l$Offset_Before
    link_rows[[i]] <- l
    c <- after$cumulative
    c$RemovedAnchorId <- rep(id, nrow(c))
    c$Offset_Before <- baseline$cumulative$Offset
    c$Change <- c$Offset - c$Offset_Before
    cumulative_rows[[i]] <- c
    after$support$RemovedAnchorId <- rep(id, nrow(after$support))
    support_rows[[i]] <- after$support
    fields <- c("LinkID", "AnchorId", "Facet", "Level", "Retained", "Contributing", "Flag")
    e <- dplyr::inner_join(baseline$elements[, fields], after$elements[, fields],
      by = c("LinkID", "AnchorId", "Facet", "Level"), suffix = c("_Before", "_After"))
    changed <- function(a, b) (is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
    retention_changed <- changed(e$Retained_Before, e$Retained_After)
    e <- e[retention_changed | changed(e$Contributing_Before, e$Contributing_After) |
      changed(e$Flag_Before, e$Flag_After), , drop = FALSE]
    e$RemovedAnchorId <- rep(id, nrow(e))
    changed_rows[[i]] <- e
    audit$MaxAbsLinkChange[i] <- max_finite(l$Change)
    audit$MaxAbsCumulativeChange[i] <- max_finite(c$Change[c$WaveID > 1L])
    audit$UnavailableLinksAfter[i] <- sum(!is.finite(l$Offset))
    audit$NewlyUnavailableLinks[i] <- sum(is.finite(l$Offset_Before) & !is.finite(l$Offset))
    audit$RetentionChanges[i] <- sum(retention_changed)
    audit$FallbackLinks[i] <- sum(l$ScreeningFallback)
    audit$CompleteComparison[i] <- all(is.finite(c$Change[c$WaveID > 1L]))
  }
  list(data = list(settings = data.frame(Method = x$config$method, DriftThreshold = threshold,
      MinCommonPerFacet = guideline, SourceEstimatesFixed = TRUE, SourceSEsFixed = TRUE, BaselineTolerance = 1e-8),
    removal = audit, baseline_links = baseline$links,
    baseline_cumulative = baseline$cumulative, baseline_elements = baseline$elements,
    links = dplyr::bind_rows(link_rows), cumulative = dplyr::bind_rows(cumulative_rows),
    retention_changes = dplyr::bind_rows(changed_rows), common_by_facet = dplyr::bind_rows(support_rows)),
    notes = .mfrmr_table_.equating_offset_sensitivity_data)
}

.draw_equating_offset_sensitivity <- function(out, style) {
  rows <- out$data$data$removal
  apply_plot_preset(style)
  old <- graphics::par(c("mar", "xpd"))
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = c(if (out$data$display$show_notes) 8 else 6.8, 8,
    if (out$data$display$show_title) 3 else 1, 1))
  title <- .mfrm_plot_display_title(out$data)
  if (!nrow(rows)) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "No common elements to remove.", cex = 0.85)
    return(invisible(out))
  }
  # Keep input order for ties at the precision printed next to the points.
  rows <- rows[order(signif(rows$MaxAbsCumulativeChange, 3), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  value <- rows$MaxAbsCumulativeChange
  available <- is.finite(value)
  position <- ifelse(available, value, 0)
  y <- rev(seq_len(nrow(rows)))
  key <- out$data$styles
  k <- ifelse(!available, 3L, ifelse(rows$CompleteComparison, 1L, 2L))
  plot(position, y, type = "n", xlim = c(0, max(0.1, position) * 1.4),
    ylim = c(0.5, nrow(rows) + 0.5), yaxt = "n", ylab = "", xaxs = "i",
    xaxt = if (any(available)) "s" else "n",
    xlab = if (any(available)) "Max. absolute cumulative offset change (logits)" else
      "No finite cumulative offset comparison", main = title, cex.lab = 0.8)
  labels <- sprintf("%s: %s", rows$Facet, rows$Level)
  duplicate <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate] <- paste0(labels[duplicate], " [", rows$AnchorId[duplicate], "]")
  labels <- truncate_axis_label(labels, width = 22L)
  wide <- graphics::strwidth(labels, cex = 0.7, units = "inches") > graphics::par("mai")[2] - 0.25
  if (any(wide)) {
    labels[wide] <- rows$AnchorId[wide]
    warning("Some sensitivity labels exceed the left margin; anchor IDs are shown. Full identities remain in the removal table.", call. = FALSE)
  }
  axis(2, at = y, labels = labels, las = 1, cex.axis = 0.7, tick = FALSE)
  if (any(available)) segments(0, y[available], value[available], y[available], col = style$neutral)
  points(position, y, pch = key$Shape[k], col = key$Color[k], cex = 1.1, xpd = NA)
  text(position, y, ifelse(available, format(signif(value, 3), trim = TRUE), "unavailable"),
    pos = 4, offset = 0.4, cex = 0.7)
  if (graphics::strheight("M", cex = 0.7) * 1.4 > 1) warning(
    "Sensitivity labels are crowded; increase the device height or use the returned removal table.", call. = FALSE)
  graphics::par(xpd = NA)
  usr <- graphics::par("usr")
  lh <- graphics::par("csi") * graphics::par("mex") * diff(usr[3:4]) / graphics::par("pin")[2]
  legend(mean(usr[1:2]), usr[3] - 3.9 * lh, xjust = 0.5, yjust = 1,
    legend = key$Label, pch = key$Shape, col = key$Color, bty = "n", cex = 0.7)
  if (out$data$display$show_notes) mtext(
    "Offsets recomputed; source fits and SEs fixed. No new uncertainty estimates.", side = 1, line = 6.7, cex = 0.6)
  invisible(out)
}

.draw_equating_topology <- function(out, style, removal = FALSE) {
  data <- out$data$data
  apply_plot_preset(style)
  old <- graphics::par(c("mar", "xpd"))
  on.exit(graphics::par(old), add = TRUE)
  title <- .mfrm_plot_display_title(out$data)
  bottom <- if (removal) {
    if (out$data$display$show_notes) 8 else 6.8
  } else if (out$data$display$show_notes) 6 else 4.5
  top <- if (out$data$display$show_title) 3 else 1
  key <- out$data$styles
  if (!removal) {
    graphics::par(mar = c(bottom, 1, top, 1))
    nodes <- data$nodes
    angle <- seq(0, 2 * pi, length.out = nrow(nodes) + 1L)[-1L] + pi / 2
    xy <- cbind(cos(angle), sin(angle))
    plot(xy, type = "n", xlim = c(-1.65, 1.65), ylim = c(-1.5, 1.5),
      axes = FALSE, xlab = "", ylab = "", main = title)
    edges <- data$links
    a <- match(edges$from, nodes$NodeId)
    b <- match(edges$to, nodes$NodeId)
    k <- match(edges$Status, key$Status)
    segments(xy[a, 1], xy[a, 2], xy[b, 1], xy[b, 2],
      lty = key$Linetype[k], col = key$Color[k], lwd = 1.7)
    # Counts are annotations of recorded links, not edge weights for distance.
    if (nrow(edges)) text((xy[a, 1] + xy[b, 1]) / 2, (xy[a, 2] + xy[b, 2]) / 2,
      labels = edges$RetainedElements, pos = 3, offset = 0.3, cex = 0.8)
    points(xy, pch = 22, bg = style$background, cex = 1.8)
    labels <- nodes$Label
    duplicate <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    labels[duplicate] <- paste0(labels[duplicate], " [W", nodes$WaveID[duplicate], "]")
    labels <- truncate_axis_label(labels, width = 18L)
    usr <- graphics::par("usr")
    too_wide <- graphics::strwidth(labels, cex = 0.75) >
      1.9 * pmin(xy[, 1] - usr[1], usr[2] - xy[, 1])
    if (any(too_wide)) {
      labels[too_wide] <- paste0("Wave ", nodes$WaveID[too_wide])
      warning("Some wave labels exceed the panel width; wave IDs are shown. Full labels remain in the node table.", call. = FALSE)
    }
    if (nrow(nodes) > 8L) warning("Wave-link labels may be crowded; enlarge the device or use the returned link table.", call. = FALSE)
    text(xy[, 1], xy[, 2], labels, pos = ifelse(xy[, 2] >= 0, 3, 1),
      offset = 0.8, cex = 0.75)
    graphics::par(xpd = NA)
    legend("bottom", inset = c(0, -0.32), legend = key$Label,
      lty = key$Linetype, col = key$Color, lwd = 1.7, bty = "n", cex = 0.7)
    footer <- "Numbers: retained elements. Connections do not establish invariance."
  } else {
    graphics::par(mar = c(bottom, 8, top, 1))
    rows <- data$removal
    rows <- rows[order(rows$NewlyDisconnectedPairs, decreasing = TRUE), , drop = FALSE]
    if (!nrow(rows)) {
      plot.new()
      title(main = title)
      text(0.5, 0.5, "No common elements to remove.", cex = 0.85)
    } else {
      value <- rows$NewlyDisconnectedPairs
      y <- rev(seq_len(nrow(rows)))
      k <- ifelse(value > 0, 1L, ifelse(rows$RetainedOccurrences > 0, 2L, 3L))
      plot(value, y, type = "n", xlim = c(0, max(1, value) * 1.15),
        ylim = c(0.5, nrow(rows) + 0.5), yaxt = "n", xaxt = "n", ylab = "",
        xlab = "Newly disconnected wave pairs", main = title, xaxs = "i")
      ticks <- pretty(c(0, max(1, value)), n = 4)
      ticks <- ticks[ticks >= 0 & ticks == floor(ticks) & ticks <= max(1, value)]
      axis(1, at = ticks, labels = ticks)
      labels <- sprintf("%s: %s", rows$Facet, rows$Level)
      duplicate <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
      labels[duplicate] <- paste0(labels[duplicate], " [", rows$AnchorId[duplicate], "]")
      labels <- truncate_axis_label(labels, width = 22L)
      too_wide <- graphics::strwidth(labels, cex = 0.7, units = "inches") > graphics::par("mai")[2] - 0.25
      if (any(too_wide)) {
        labels[too_wide] <- rows$AnchorId[too_wide]
        warning("Some removal labels exceed the left margin; anchor IDs are shown. Full labels remain in the removal table.", call. = FALSE)
      }
      axis(2, at = y, labels = labels, las = 1, cex.axis = 0.7, tick = FALSE)
      segments(0, y, value, y, col = style$neutral)
      points(value, y, pch = key$Shape[k], col = key$Color[k], cex = 1.1, xpd = NA)
      text(value, y, value, pos = 4, offset = 0.4, cex = 0.7)
      if (graphics::strheight("M", cex = 0.7) * 1.4 > 1) warning(
        "Removal labels are crowded; increase the device height or use the returned removal table.", call. = FALSE)
      graphics::par(xpd = NA)
      usr <- graphics::par("usr")
      line_height <- graphics::par("csi") * graphics::par("mex") * diff(usr[3:4]) / graphics::par("pin")[2]
      legend(mean(usr[1:2]), usr[3] - 3.9 * line_height,
        xjust = 0.5, yjust = 1, legend = key$Label,
        pch = key$Shape, col = key$Color, bty = "n", cex = 0.7)
    }
    footer <- "Topology only; screening fixed. Estimates and SEs are not recomputed."
  }
  if (out$data$display$show_notes) mtext(footer, side = 1,
    line = if (removal) 6.7 else 4.8, cex = 0.6)
  invisible(out)
}

.draw_equating_graph <- function(out, style) {
  data <- out$data$data
  nodes <- data$nodes
  edges <- data$edges
  nw <- data$n_waves
  ne <- data$n_element_nodes
  # ponytail: fixed two-column layout; add link panels if dense chains need them.
  y <- function(n) if (n == 1L) 0.5 else seq(0.95, 0.05, length.out = n)
  nodes$X <- c(rep(0.30, nw), rep(0.62, ne))
  nodes$Y <- c(y(nw), y(ne))
  key <- out$data$styles
  edge_style <- match(edges$Status, key$Status)
  apply_plot_preset(style)
  old <- graphics::par(c("mar", "xpd"))
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = c(if (out$data$display$show_notes) 6 else 4.5, 0.5,
                        if (out$data$display$show_title) 3 else 1, 0.5))
  plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  cex <- 0.72
  labels <- nodes$Label
  duplicate <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate] <- paste0(labels[duplicate], " [", nodes$NodeId[duplicate], "]")
  labels <- truncate_axis_label(labels, width = 30L)
  available <- ifelse(nodes$Kind == "wave", 0.27, 0.35)
  too_wide <- graphics::strwidth(labels, cex = cex) > available
  if (any(too_wide)) {
    # Keep stable IDs on shortened labels; full text remains in the payload.
    labels[too_wide] <- ifelse(nodes$Kind[too_wide] == "wave",
      paste0("Wave ", nodes$WaveID[too_wide]), nodes$NodeId[too_wide])
    warning("Some graph labels exceed the panel width; node IDs are shown. Widen the device or use the returned node table for full labels.", call. = FALSE)
  }
  if (max(nw, ne) > 1L && 0.9 / (max(nw, ne) - 1L) < 1.4 * graphics::strheight("M", cex = cex)) {
    warning("Graph labels are crowded; increase the device height or use the returned node/edge tables.", call. = FALSE)
  }
  a <- match(edges$from, nodes$NodeId)
  b <- match(edges$to, nodes$NodeId)
  segments(nodes$X[a], nodes$Y[a], nodes$X[b], nodes$Y[b],
    col = key$Color[edge_style], lty = key$Linetype[edge_style], lwd = 1.5)
  points(nodes$X, nodes$Y, pch = ifelse(nodes$Kind == "wave", 22, 21),
    bg = style$background, col = style$foreground, cex = 1.4)
  text(nodes$X, nodes$Y, labels, pos = ifelse(nodes$Kind == "wave", 2, 4),
    cex = cex, col = style$foreground, offset = 0.7)
  title(main = .mfrm_plot_display_title(out$data))
  graphics::par(xpd = NA)
  legend(0.5, -0.035, xjust = 0.5, yjust = 1, bty = "n", ncol = 2,
    legend = c("Wave / form", "Element in link", key$Label),
    pch = c(22, 21, rep(NA, nrow(key))), pt.bg = style$background,
    lty = c(NA, NA, key$Linetype), col = c(rep(style$foreground, 2), key$Color),
    lwd = 1.5, cex = 0.7)
  if (out$data$display$show_notes) mtext(
    "Recorded comparisons only. Retention does not establish invariance.",
    side = 1, line = 4.5, cex = 0.63)
  invisible(out)
}

.plot_equating_chain <- function(x, draw = TRUE, style = resolve_plot_preset("standard"),
                                 show_title = FALSE, show_notes = FALSE, ...) {
  cum <- x$cumulative
  out <- new_mfrm_plot_data(
    "anchor_drift",
    list(
      plot = "chain",
      table = cum,
      links = x$links,
      display = list(show_title = show_title, show_notes = show_notes),
      title = "Equating chain",
      subtitle = "Cumulative offsets across linked calibration waves",
      legend = new_plot_legend(
        label = c("Cumulative offset", "Centered chain reference"),
        role = c("offset", "reference"),
        aesthetic = c("line-point", "line"),
        value = c(style$accent_primary, style$foreground)
      ),
      reference_lines = new_reference_lines("h", 0, "Centered chain reference", "dashed", "reference"),
      preset = style$name
    )
  )
  out$data$notes <- .mfrm_plot_notes(out$data)
  if (!draw) return(invisible(out))

  opar <- graphics::par()["mar"]
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(5, 4, 3, 1))

  n <- nrow(cum)
  plot(seq_len(n), cum$Cumulative_Offset, type = "b",
                 pch = 19, col = style$accent_primary, lwd = 2,
                 xaxt = "n", xlab = "", ylab = "Cumulative Offset (logits)",
                 main = .mfrm_plot_display_title(out$data), ...)
  axis(1, at = seq_len(n), labels = cum$Wave, las = 2, cex.axis = 0.8)
  abline(h = pretty(cum$Cumulative_Offset, n = 5), col = grDevices::adjustcolor(style$grid, alpha.f = 0.85), lty = 1)
  abline(h = 0, lty = 2, col = grDevices::adjustcolor(style$foreground, alpha.f = 0.7))

  links <- x$links
  for (i in seq_len(nrow(links))) {
    mid_x <- i + 0.5
    mid_y <- (cum$Cumulative_Offset[i] + cum$Cumulative_Offset[i + 1]) / 2
    text(mid_x, mid_y, sprintf("n=%d", links$N_Common[i]),
                   cex = 0.7, col = grDevices::adjustcolor(style$foreground, alpha.f = 0.82))
  }

  invisible(out)
}

# --- Exported plot function --------------------------------------------------

#' Plot anchor drift or a screened linking chain
#'
#' Creates base-R plots for inspecting anchor drift across calibration waves
#' or visualising the cumulative offset in a screened linking chain.
#'
#' @param x An `mfrm_anchor_drift` or `mfrm_equating_chain` object.
#' @param type Plot type: `"drift"` (dot plot of element drift),
#'   `"chain"` (cumulative offset line plot), `"heatmap"`
#'   (wave-by-element drift heatmap), or `"forest"` (per-(Facet, Level,
#'   Wave) anchor estimate with `+/- z * SE` whiskers; requires
#'   `mfrm_anchor_drift`).
#' @param facet Optional character vector to filter drift plots to specific
#'   facets.
#' @param ci_level Confidence level used by `type = "forest"` for the
#'   anchor-estimate whiskers (default `0.95`). Ignored for other
#'   plot types.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
#' @param draw If `FALSE`, return the plot data invisibly without drawing.
#' @param ... Additional graphical parameters passed to base plotting
#'   functions.
#'
#' @details
#' Four plot types are supported:
#'
#' - **`"drift"`** (for `mfrm_anchor_drift` objects): A dot plot of each
#'   element's drift value. Vertical reference lines mark the drift threshold.
#'   The preset's warning colour indicates flagged elements.
#' - **`"heatmap"`** (for `mfrm_anchor_drift` objects): A wave-by-element
#'   heat matrix showing drift magnitude.  Darker cells represent larger
#'   absolute drift.  Useful for spotting systematic patterns (e.g., all
#'   criteria shifting in the same direction).
#' - **`"forest"`** (for `mfrm_anchor_drift` objects): Per-wave estimates
#'   with approximate normal whiskers from their available SEs. These are
#'   not intervals for between-wave differences; linking-offset uncertainty
#'   and cross-wave covariance are not included.
#' - **`"chain"`** (for `mfrm_equating_chain` objects): A line plot of
#'   cumulative offsets across the screened linking chain. A flatter line
#'   indicates smaller between-wave shifts; steep segments suggest larger
#'   link offsets that deserve review.
#'
#' @section Which plot should I use?:
#' - Use `type = "drift"` with an `mfrm_anchor_drift` object to review flagged
#'   elements directly.
#' - Use `type = "heatmap"` with an `mfrm_anchor_drift` object to spot
#'   wave-by-element patterns.
#' - Use `type = "chain"` with an `mfrm_equating_chain` object after
#'   [build_equating_chain()] to inspect cumulative offsets across waves.
#'
#' @section Interpreting plots:
#' **Drift** is the change in an element's estimated measure between
#' calibration waves, after accounting for the screened common-element link
#' offset. An
#' element is flagged when its absolute drift exceeds a threshold
#' (typically 0.5 logits) **and** the drift-to-SE ratio exceeds a
#' secondary criterion (typically 2.0). These are heuristic screening rules,
#' not a calibrated test of change. Subtracting a common offset can remove
#' a panel-wide shift; unflagged elements do not establish absence of annual drift.
#'
#' - Drift dots distinguish flagged elements. Heatmap shading shows magnitude,
#'   not whether both screening thresholds were exceeded. Possible explanations
#'   require substantive review, including changes in scoring or exposure.
#' - In chain plots, uneven changes between waves indicate differential
#'   shifts in the screened linking offsets. The \eqn{y}-axis shows cumulative
#'   logit-scale offsets; flatter segments indicate more stable adjacent links.
#'   Steep segments should be checked alongside `LinkSupportAdequate` and the
#'   retained common-element counts before making longitudinal claims.
#' - For drift objects, it is usually best to read `summary(x)` first
#'   and then use the plot to see where the flagged values sit.
#'
#' @section Typical workflow:
#' 1. Build a drift or screened-linking object with [detect_anchor_drift()] or
#'    [build_equating_chain()].
#' 2. Start with `draw = FALSE` if you want the plotting data for custom
#'    reporting.
#' 3. Use the base-R plot for quick screening and then inspect the underlying
#'    tables for exact values.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @return A plotting-data object of class `mfrm_plot_data`. With
#'   `draw = FALSE`, `result$data$table` contains the filtered drift or chain
#'   table, `result$data$matrix` contains the heatmap matrix when requested,
#'   and the returned plot data includes package-native `title`, `subtitle`,
#'   `legend`, and `reference_lines`.
#'
#' @seealso [detect_anchor_drift()], [build_equating_chain()],
#'   [plot(x, type = "heatmap")][plot.mfrm_dif], [plot(x, type = "bubble")][plot.mfrm_fit], [mfrmr_visual_diagnostics]
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept linking
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_anchor_drift <- function(x, type = c("drift", "chain", "heatmap", "forest"),
                              facet = NULL,
                              ci_level = 0.95,
                              preset = c("standard", "publication", "compact", "monochrome"),
                              draw = TRUE, ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  type <- match.arg(type)
  style <- resolve_plot_preset(preset)

  if (inherits(x, "mfrm_equating_chain")) {
    if (type == "chain") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_equating_chain(x, draw = draw, style = style, ...))
    }
  }

  if (inherits(x, "mfrm_anchor_drift")) {
    dt <- x$drift_table
    if (!is.null(facet)) dt <- dt |> dplyr::filter(.data$Facet %in% facet)

    if (nrow(dt) == 0) {
      if (draw) message("No drift data to plot.")
      return(invisible(NULL))
    }

    if (type == "drift") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_drift_dot(dt, x$config, draw = draw, style = style, ...))
    } else if (type == "heatmap") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_drift_heatmap(dt, x$config, draw = draw, style = style, ...))
    } else if (type == "forest") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_drift_forest(dt, x$config, ci_level = ci_level,
                                 draw = draw, style = style, ...))
    }
  }

  stop("Unsupported object class or plot type combination.", call. = FALSE)
}

# Forest-style anchor-drift visualisation. Each (Facet, Level, Wave)
# row is rendered as a horizontal CI whisker around the per-wave
# anchor estimate; rows are grouped by Facet level so reviewers see
# one ladder of CI bands per anchored element across waves.
.plot_drift_forest <- function(drift_tbl, config, ci_level = 0.95,
                               draw = TRUE, style = NULL, ...) {
  if (is.null(style)) style <- resolve_plot_preset("standard")
  dt <- as.data.frame(drift_tbl, stringsAsFactors = FALSE)
  est_col <- if ("Estimate" %in% names(dt)) "Estimate" else if ("Wave_Est" %in% names(dt)) "Wave_Est" else NA_character_
  if (is.na(est_col) || !all(c("Facet", "Level", "Wave") %in% names(dt))) {
    stop("Forest plot requires Facet, Level, Wave, and Estimate (or Wave_Est) columns.",
         call. = FALSE)
  }
  dt$Estimate <- dt[[est_col]]
  se_col <- if ("SE_Wave" %in% names(dt)) "SE_Wave" else if ("SE" %in% names(dt)) "SE" else if ("ModelSE" %in% names(dt)) "ModelSE" else NA_character_
  z_ci <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
  if (!is.na(se_col)) {
    dt$CI_Lower <- dt$Estimate - z_ci * dt[[se_col]]
    dt$CI_Upper <- dt$Estimate + z_ci * dt[[se_col]]
  } else {
    dt$CI_Lower <- NA_real_
    dt$CI_Upper <- NA_real_
  }
  dt$Label <- paste0(dt$Facet, ":", dt$Level, " @ ", dt$Wave)
  dt <- dt[order(dt$Facet, dt$Level, dt$Wave), , drop = FALSE]
  out <- new_mfrm_plot_data(
    "anchor_drift_forest",
    list(
      data = dt,
      ci_level = ci_level,
      title = "Anchor drift forest",
      subtitle = sprintf("Per-wave anchor estimates with %g%% CI",
                         round(100 * ci_level)),
      preset = style$name,
      legend = new_plot_legend(
        label = c("Anchor estimate", "CI whisker"),
        role = c("location", "uncertainty"),
        aesthetic = c("point", "segment"),
        value = c(style$accent_primary, style$accent_primary)
      ),
      reference_lines = new_reference_lines("v", 0,
                                             "Centred logit reference",
                                             "dashed", "reference")
    )
  )
  if (isTRUE(draw)) {
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    label_width <- max(graphics::strwidth(dt$Label, units = "inches", cex = .7))
    if (label_width + graphics::par("mai")[4] + 1.25 >= graphics::par("fin")[1]) {
      stop("Anchor labels need more horizontal space; use a wider graphics device or inspect the table with `draw = FALSE`.",
        call. = FALSE)
    }
    margin <- label_width / (graphics::par("csi") * graphics::par("mex")) + 1.2
    graphics::par(mar = c(old_mar[1], max(old_mar[2], margin), old_mar[3:4]))
    y <- seq_len(nrow(dt))
    xrange <- range(c(dt$Estimate, dt$CI_Lower, dt$CI_Upper, 0),
                    finite = TRUE, na.rm = TRUE)
    plot(
      x = dt$Estimate, y = y, type = "n",
      xlim = xrange, yaxt = "n",
      xlab = "Estimate (logit)", ylab = "",
      main = ""
    )
    abline(v = 0, lty = 2, col = style$neutral)
    valid <- is.finite(dt$CI_Lower) & is.finite(dt$CI_Upper)
    if (any(valid)) {
      segments(dt$CI_Lower[valid], y[valid],
                         dt$CI_Upper[valid], y[valid],
                         col = style$accent_primary, lwd = 2)
    }
    points(dt$Estimate, y, pch = 19, col = style$accent_primary)
    axis(2, at = y, labels = dt$Label, las = 1, cex.axis = 0.7)
  }
  invisible(out)
}

# ============================================================================
# QC Pipeline Plot
# ============================================================================

#' Plot QC pipeline results
#'
#' Visualizes the output from [run_qc_pipeline()] as either a traffic-light
#' bar chart or a panel containing each check's explanation.
#'
#' @param x Output from [run_qc_pipeline()].
#' @param type Plot type: `"traffic_light"` (default) or `"detail"`.
#' @param draw If `FALSE`, return plot data invisibly without drawing.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @details
#' Two plot types are provided for visual triage of QC results:
#'
#' - **`"traffic_light"`** (default): A horizontal bar chart with one row
#'   per QC check.  Bars are coloured green (Pass), amber (Warn), or red
#'   (Fail).  Provides an at-a-glance summary of the current QC review state.
#' - **`"detail"`**: A panel showing the saved `Detail` text beside each
#'   verdict. Read this explanation and the returned table to understand
#'   why a check was flagged or unavailable.
#'
#' @section Interpretation:
#' The rows are the actual checks returned by [run_qc_pipeline()]: convergence,
#' global fit, selected-facet reliability and separation, element misfit,
#' unexpected responses, category structure, connectivity, inter-rater
#' agreement, and the functioning/bias screen. Element misfit uses the
#' configured mean-square band; it is not a ZSTD test. PCA and displacement
#' are separate diagnostics and are not checks in this pipeline.
#'
#' Green (Pass), amber (Warn), and red (Fail) describe the selected screening
#' rules. Grey (Skip) marks unrequested or unavailable checks; `AffectsOverall`
#' in the returned data distinguishes their contribution to the overall result.
#' A missing result is not evidence of acceptable fit. High rater separation
#' is not high agreement, and no colour establishes statistical validity or
#' justifies deleting raters or collapsing categories without further review.
#'
#' @return Invisible verdicts tibble from the QC pipeline.
#'
#' @seealso [run_qc_pipeline()], [plot(x, type = "qc")][plot.mfrm_fit],
#'   [build_visual_summaries()], [mfrmr_visual_diagnostics]
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_qc_pipeline <- function(x, type = c("traffic_light", "detail"),
                             draw = TRUE, ...) {
  type <- match.arg(type)
  stopifnot(inherits(x, "mfrm_qc_pipeline"))
  validate_qc_pipeline_output(x)

  vt <- x$verdicts
  if (!draw) return(invisible(vt))

  n <- nrow(vt)
  cols <- ifelse(vt$Verdict == "Pass", "#2ca02c",
                 ifelse(vt$Verdict == "Warn", "#ff7f0e",
                        ifelse(vt$Verdict == "Fail", "#d62728", "#999999")))

  opar <- graphics::par()["mar"]
  on.exit(graphics::par(opar), add = TRUE)

  if (type == "traffic_light") {
    graphics::par(mar = c(3, 14, 3, 4))

    plot(NULL, xlim = c(0, 1), ylim = c(0.5, n + 0.5),
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = "", ...)

    for (i in seq_len(n)) {
      rect(0, i - 0.4, 1, i + 0.4, col = cols[i], border = NA)
      text(0.5, i, vt$Verdict[i], col = "white", font = 2, cex = 0.9)
    }

    axis(2, at = seq_len(n), labels = vt$Check,
                   las = 1, cex.axis = 0.8, tick = FALSE)

    for (i in seq_len(n)) {
      mtext(vt$Value[i], side = 4, at = i,
                      las = 1, cex = 0.6, line = 0.5)
    }
  } else {
    graphics::par(mar = c(3, 14, 3, 8))

    plot(NULL, xlim = c(0, 1), ylim = c(0.5, n + 0.5),
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = "", ...)

    for (i in seq_len(n)) {
      rect(0, i - 0.4, 0.15, i + 0.4, col = cols[i], border = NA)
      text(0.075, i, substr(vt$Verdict[i], 1, 1),
                     col = "white", font = 2, cex = 0.8)
      text(0.2, i, vt$Detail[i], adj = 0, cex = 0.7)
    }

    axis(2, at = seq_len(n), labels = vt$Check,
                   las = 1, cex.axis = 0.8, tick = FALSE)
  }

  invisible(vt)
}
