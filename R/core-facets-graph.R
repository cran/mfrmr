# Readers and descriptive consumers of exported FACETS Graphfile coordinates.
# They do not construct a native model or infer item/scale correspondences.
.facets_graph_number <- function(token) {
  token <- trimws(token)
  valid <- !is.na(token) & grepl("^[+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?$", token)
  value <- suppressWarnings(as.numeric(token))
  reason <- rep("available", length(token))
  reason[!valid] <- "invalid_source_token"
  reason[is.na(token) | token %in% c("", ".", "*", "NA", "N/A", "NaN", "Inf", "-Inf")] <-
    "missing_or_nonfinite"
  reason[!is.na(value) & !is.finite(value)] <- "missing_or_nonfinite"
  reason[!is.na(token) & grepl("[<>]", token)] <- "censored_source_value"
  value[!valid | !is.finite(value)] <- NA_real_
  mantissa <- sub("[eE].*$", "", token)
  decimals <- fit_review_reported_decimal_places(mantissa)
  exponent <- rep(0, length(token))
  scientific <- valid & grepl("[eE]", token)
  exponent[scientific] <- suppressWarnings(as.numeric(sub(".*[eE]", "", token[scientific])))
  resolution <- 10^(exponent - decimals)
  resolution[reason != "available"] <- NA_real_
  list(value = value, raw = token, resolution = resolution, reason = reason)
}

.facets_graph_read <- function(path, encoding, source_id) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path) || dir.exists(path))
    stop("`graph_file` must name one existing text Graphfile.", call. = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  hash <- unname(tools::md5sum(path))
  lines <- fit_review_read_lines(path, encoding)
  source_rows <- which(nzchar(trimws(lines)))
  lines <- lines[source_rows]
  if (length(lines) < 2L) stop("Graphfile needs headings and data rows.", call. = FALSE)
  sep <- if (grepl("\t", lines[1], fixed = TRUE)) "\t" else
    if (grepl(",", lines[1], fixed = TRUE)) "," else NULL
  if (!is.null(sep)) {
    con <- textConnection(lines)
    counts <- tryCatch(utils::count.fields(con, sep = sep, quote = '"', comment.char = "",
      blank.lines.skip = FALSE), finally = close(con))
    if (anyNA(counts) || any(counts != counts[1]))
      stop("Graphfile rows must have the same fields as their single-line headings.", call. = FALSE)
    raw <- utils::read.table(text = paste(lines, collapse = "\n"), header = FALSE,
      sep = sep, quote = '"', comment.char = "", colClasses = "character",
      na.strings = NULL, check.names = FALSE, fill = FALSE, encoding = "UTF-8")
    header <- trimws(as.character(raw[1, ]))
    raw <- raw[-1, , drop = FALSE]
    layout <- if (sep == ",") "csv" else "tab"
  } else {
    if (any(grepl("[^\x01-\x7F]", lines)))
      stop("Fixed-width Graphfiles must be ASCII; use named CSV/TSV for other text.", call. = FALSE)
    width <- 10L
    fields <- ceiling(nchar(lines[1], type = "chars") / width)
    if (any(nchar(lines, type = "chars") > fields * width))
      stop("Graphfile exceeds its fixed-width headings.", call. = FALSE)
    split_fields <- function(line) trimws(substring(line,
      (seq_len(fields) - 1L) * width + 1L, seq_len(fields) * width))
    header <- split_fields(lines[1])
    raw <- as.data.frame(do.call(rbind, lapply(lines[-1], split_fields)),
      stringsAsFactors = FALSE)
    layout <- "fixed_width"
  }
  if (length(header) < 5L || anyDuplicated(header) ||
      !identical(header[1:4], c("Scale", "Measure", "Expected", "ExpCat")) ||
      any(!grepl("^Prob:.+", header[-(1:4)])))
    stop("Graphfile requires unique Scale, Measure, Expected, ExpCat and Prob:<category> headings.", call. = FALSE)
  names(raw) <- header
  row.names(raw) <- NULL
  scale_num <- .facets_graph_number(raw$Scale)
  if (any(scale_num$reason != "available" | scale_num$value < 1 |
          scale_num$value != floor(scale_num$value)))
    stop("Graphfile scale numbers must be positive integers on every row.", call. = FALSE)
  scale <- as.character(scale_num$value)
  measure <- .facets_graph_number(raw$Measure)
  expected <- .facets_graph_number(raw$Expected)
  expcat <- .facets_graph_number(raw$ExpCat)
  finite <- is.finite(measure$value)
  if (anyDuplicated(data.frame(Scale = scale[finite], Measure = measure$value[finite])))
    stop("Duplicate Graphfile scale/measure coordinates; export distinct coordinates at sufficient precision.", call. = FALSE)
  for (s in unique(scale)) {
    delta <- diff(measure$value[scale == s & finite])
    if (length(delta) && !all(delta > 0) && !all(delta < 0))
      stop("Graphfile measures must be monotone within each scale; rows are not reordered.", call. = FALSE)
  }
  prob_names <- header[-(1:4)]
  categories <- sub("^Prob:", "", prob_names)
  probs <- lapply(raw[prob_names], .facets_graph_number)
  for (i in seq_along(probs)) {
    p <- probs[[i]]
    p$reason[is.finite(p$value) & (p$value < 0 | p$value > 1)] <- "probability_out_of_range"
    probs[[i]] <- p
  }
  available <- vapply(probs, function(p) p$reason == "available", logical(nrow(raw)))
  values <- vapply(probs, `[[`, numeric(nrow(raw)), "value")
  resolution <- vapply(probs, `[[`, numeric(nrow(raw)), "resolution")
  # Preserve the matrix dimensions for one-row or one-category exports.
  dim(available) <- dim(values) <- dim(resolution) <- c(nrow(raw), length(probs))
  complete <- rowSums(available) == length(probs)
  sums <- rowSums(values)
  allowance <- rowSums(resolution / 2)
  outside <- complete & abs(sums - 1) > allowance + 16 * .Machine$double.eps * length(probs)
  status <- ifelse(complete, ifelse(outside, "outside_rounding_allowance", "within_rounding_allowance"),
    "incomplete_or_invalid")
  diagnostics <- data.frame(SourceID = source_id, SourceRow = source_rows[-1], Scale = scale,
    Measure = measure$value, ProbabilitySum = sums, RoundingAllowance = allowance,
    ProbabilityStatus = status)
  series <- c(list(expected), probs)
  curves <- lapply(seq_along(series), function(i) {
    p <- series[[i]]
    reason <- p$reason
    if (i > 1L) reason[outside] <- "probability_sum_outside_rounding"
    reason[measure$reason != "available"] <- paste0("measure:", measure$reason[measure$reason != "available"])
    data.frame(SourceID = source_id, SourceRow = source_rows[-1], Scale = scale,
      ScaleRaw = raw$Scale, Measure = measure$value, MeasureRaw = measure$raw,
      MeasureResolution = measure$resolution,
      Quantity = if (i == 1L) "expected" else "probability",
      Category = if (i == 1L) NA_character_ else categories[i - 1L],
      Value = p$value, ValueRaw = p$raw, ValueResolution = p$resolution,
      ExpCat = expcat$value, ExpCatRaw = expcat$raw,
      Available = reason == "available", Reason = reason, ProbabilityStatus = status)
  })
  if (!identical(unname(tools::md5sum(path)), hash))
    stop("A Graphfile changed during import; retry after the export finishes.", call. = FALSE)
  list(table = raw, curves = do.call(rbind, curves), diagnostics = diagnostics,
    categories = categories, layout = layout,
    source = data.frame(SourceID = source_id, Path = path, MD5 = hash,
      Encoding = encoding, Input = "graph_file"))
}

.facets_results_add_graph <- function(x, path, encoding) {
  graph <- .facets_graph_read(path, encoding, nrow(x$sources) + 1L)
  if (graph$source$Path %in% x$sources$Path)
    stop("The Graphfile must be separate from the element files.", call. = FALSE)
  x$sources <- rbind(x$sources, graph$source)
  graph$source <- NULL
  graph$metadata <- list(software = "FACETS", version = x$metadata$version,
    coordinate = "Reported Graphfile Measure", measure_units = NA_character_,
    scale_identity = "Table 8 subtable", element_scale_alignment = "not_inferred")
  x$graph <- graph
  x$notes <- rbind(x$notes, data.frame(Type = c("graph_source", "graph_scale", "graph_availability"),
    Text = c("Reported Graphfile coordinates; no refitting, confidence intervals or new model predictions. Connecting segments are visual guides.",
      "Graphfile scales identify Table 8 subtables. Source category coding and Measure origins are preserved; no item mapping or alignment to Scorefile measures is inferred.",
      "Missing/invalid values break lines. Probability sums are checked against printed rounding precision without renormalization. Available values from incomplete vectors remain explicitly flagged; inspect graph diagnostics.")))
  x
}

.facets_graph_curves <- function(x) {
  if (is.null(x$graph)) stop("No Graphfile was imported; supply `graph_file` when reading results.", call. = FALSE)
  x$graph$curves
}

.facets_graph_summary <- function(x) {
  d <- .facets_graph_curves(x)
  do.call(rbind, lapply(unique(d$Scale), function(s) {
    z <- d[d$Scale == s, , drop = FALSE]
    diagnostics <- x$graph$diagnostics[x$graph$diagnostics$Scale == s, , drop = FALSE]
    data.frame(Scale = s, GridRows = nrow(diagnostics), Categories = length(x$graph$categories),
      ExpectedAvailable = sum(z$Quantity == "expected" & z$Available),
      ProbabilityAvailable = sum(z$Quantity == "probability" & z$Available),
      ProbabilityCheckIssues = sum(diagnostics$ProbabilityStatus != "within_rounding_allowance"))
  }))
}

.facets_graph_plot <- function(x, type, scales, draw, show_title, show_notes, show_legend) {
  d <- .facets_graph_curves(x)
  if (is.null(scales)) scales <- unique(d$Scale)
  if (!is.character(scales) && !is.numeric(scales))
    stop("`scales` must select Graphfile scale numbers.", call. = FALSE)
  scales <- as.character(scales)
  if (!length(scales) || anyNA(scales) || anyDuplicated(scales) || any(!scales %in% d$Scale))
    stop("`scales` must select distinct scale numbers in the Graphfile.", call. = FALSE)
  quantity <- if (type == "ccc") "probability" else "expected"
  d <- d[d$Scale %in% scales & d$Quantity == quantity, , drop = FALSE]
  d$Panel <- paste("Scale", d$Scale)
  d$Series <- if (type == "ccc") d$Category else "Expected score"
  labels <- unique(d$Series)
  colors <- if (type == "ccc") grDevices::hcl.colors(length(labels), "Dark 3") else "black"
  d$Color <- colors[match(d$Series, labels)]
  d$LineType <- (match(d$Series, labels) - 1L) %% 6L + 1L
  d$X <- ifelse(d$Available, d$Measure, NA_real_)
  d$Y <- ifelse(d$Available, d$Value, NA_real_)
  d$Segment <- NA_integer_
  for (s in scales) for (label in labels) {
    idx <- which(d$Scale == s & d$Series == label)
    segment <- cumsum(!d$Available[idx])
    d$Segment[idx] <- ifelse(d$Available[idx], segment, NA_integer_)
  }
  xlab <- "Graphfile measure"
  ylab <- if (type == "ccc") "Category probability" else "Reported expected score"
  out <- new_mfrm_plot_data(paste0("facets_results_", type), list(
    data = x$graph$table, plot_data = d, diagnostics = x$graph$diagnostics,
    sources = x$sources, metadata = x$graph$metadata, notes = x$notes,
    title = "Imported FACETS curves", xlab = xlab, ylab = ylab,
    legend = new_plot_legend(label = labels, role = if (type == "ccc") "source_category" else "expected_score",
      aesthetic = "color", value = colors),
    display = list(show_title = show_title, show_notes = show_notes, show_legend = show_legend)))
  if (draw) {
    nc <- min(2L, length(scales)); nr <- ceiling(length(scales) / nc)
    key <- show_legend && type == "ccc"
    key_space <- if (key) 1 + ceiling(length(labels) / 4) else 0
    old <- graphics::par(mfrow = c(nr, nc), mar = c(4.5, 4.5, if (show_title) 3 else 1, 1),
      oma = c(if (show_notes) 3 else 0, 0, 0, 0), xpd = FALSE)
    on.exit(graphics::par(old), add = TRUE)
    if (key) {
      # Reserve a legend region instead of freezing device-to-user coordinates.
      panel_cex <- graphics::par("cex")
      panel_layout <- matrix(c(seq_along(scales), rep(0L, nr * nc - length(scales))),
        nrow = nr, byrow = TRUE)
      graphics::layout(rbind(panel_layout, rep(length(scales) + 1L, nc)),
        heights = c(rep(1, nr), graphics::lcm(2.54 * graphics::par("csi") * key_space)))
      graphics::par(cex = panel_cex)
    }
    limits <- function(z) {
      z <- z[is.finite(z)]
      if (!length(z)) return(c(0, 1))
      r <- range(z)
      if (diff(r) == 0) r <- r + c(-.5, .5)
      r
    }
    for (s in scales) {
      z <- d[d$Scale == s, , drop = FALSE]
      graphics::plot(NA_real_, NA_real_, xlim = limits(z$X),
        ylim = if (type == "ccc") c(0, 1) else limits(z$Y),
        xlab = paste0(xlab, " [Scale ", s, "]"), ylab = ylab,
        main = if (show_title) out$data$title else "")
      for (label in labels) {
        v <- z[z$Series == label, , drop = FALSE]
        if (!nrow(v)) next
        graphics::lines(v$X, v$Y, col = v$Color[1], lty = v$LineType[1], lwd = 1.5)
        # Show isolated available coordinates even when no line can be drawn.
        isolated <- !is.na(v$Segment) & !duplicated(v$Segment) & !duplicated(v$Segment, fromLast = TRUE)
        graphics::points(v$X[isolated], v$Y[isolated], col = v$Color[1], pch = 16)
      }
      if (!any(z$Available)) graphics::text(mean(graphics::par("usr")[1:2]),
        mean(graphics::par("usr")[3:4]), "No available coordinates")
    }
    if (key) {
      graphics::par(mar = rep(0, 4))
      graphics::plot.new()
      graphics::legend("center", legend = labels, col = colors,
        lty = (seq_along(labels) - 1L) %% 6L + 1L,
        ncol = min(4L, length(labels)), cex = .8, bty = "n", xpd = NA, title = "Source category")
    }
    if (show_notes) graphics::mtext("Reported FACETS coordinates; inspect source scales and saved diagnostics.",
      side = 1, outer = TRUE, line = .5, cex = .7)
  }
  out
}
