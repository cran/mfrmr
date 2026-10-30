#' Read FACETS results for standalone visualization
#'
#' Import reported element measures, fit statistics and optional Graphfile
#' curves without refitting a model.
#' Use [plot.mfrm_facets_results()] to inspect the saved results.
#'
#' @inheritParams read_facets_fit_table
#' @param file Scorefile input as described below, or `NULL` when only
#'   `graph_file` is supplied.
#' @param graph_file Optional path to one FACETS Graphfile with headings, in
#'   CSV, tab-separated or validated ASCII fixed-width format. Supply an export
#'   from the same analysis when also importing element results. File identity
#'   is recorded, but the reader cannot certify that the exports belong together.
#' @param source_version Optional FACETS version recorded with these exports.
#'   Unknown versions remain `NA`; the reader does not infer a version.
#' @param measure_units Optional description of the reported measure units,
#'   after checking FACETS reporting transformations. Unknown units remain `NA`.
#' @param common_scale Logical; explicitly declare that the imported facets
#'   have compatible units and origins for a shared-axis Wright map. Defaults
#'   to false. This declaration does not transform, center or reverse measures.
#' @details
#' Each file must describe one analysis. Supply a vector of facet files or a
#' directory containing `score.N.txt` files, or a named delimited table. Comma,
#' tab and semicolon Scorefiles retain their native facet heading lines. Numeric
#' facet numbers default to `Facet1`, `Facet2`, etc.; supply `facet_map` to name
#' them. A data frame is also accepted. The supported layouts follow
#' [read_facets_fit_table()]. Use named delimited files for non-ASCII labels; fixed-width
#' imports require ASCII because byte and character widths can differ.
#' For text output from FACETS, request `Scorefile=score.txt`,
#' `Heading lines=Yes`, `QM=Double` and `CSV=Commas`, `Tabs` or `Semicolons`.
#' Check that the files were actually generated; a zero process exit code alone
#' does not establish a successful export. Specify the encoding of the existing
#' file; the reader does not configure FACETS input encoding or repair labels
#' already lost during export.
#' Element IDs are preserved separately from display
#' labels. Where no ID was exported, the label is the key; ambiguous duplicate
#' keys are rejected, not averaged. Unknown source status remains unknown.
#'
#' FACETS status codes distinguish anchor-flagged, estimated, extreme, no-data,
#' unmeasurable and single-observation elements. `StatusRaw = "1"` and the
#' retained `SourceStatus = "anchored"` mean an anchor flag, not proof that the
#' individual measure was fixed. FACETS 4.5.0 uses this flag for estimated
#' group-anchor members as well as direct anchors and singleton groups.
#' `GroupRaw` can also describe a reporting group; it does not establish the
#' anchor subtype. Neither a group's imported size nor a reported SE resolves
#' that ambiguity. Check the original `Labels=` specification for the anchor
#' mode and Table 7 for row annotations. A Table 7 G annotation alone need not
#' distinguish the G and X input modes. This reader does not parse those sources
#' or infer a subtype.
#' Extreme group members can instead carry Status -2/-3, so `AnchorFlagged`
#' does not count every member of an anchored group. Finite extreme measures
#' are retained as reported locations, not ordinary finite maximum-likelihood
#' estimates. The status alone does not identify their adjustment rule:
#' native X group anchoring can report a different location from G anchoring.
#' Reported fit values for extremes can be placeholders, so their fit plots are
#' unavailable. No-data and unmeasurable elements remain in the returned table.
#' Missing or censored values are never replaced with zero. Marginal source SEs
#' are retained for inspection, but do not supply joint covariance or validated
#' confidence intervals. The object contains no likelihood, response data,
#' inferred thresholds, new-person scoring or calibration-inclusive inference.
#'
#' Source measures, signs and rounding are retained. A positive or negative sign
#' describes the source model contribution; the plot does not multiply measures
#' by that sign. Check the original specification before interpreting severity,
#' ability or comparing facet origins. Reading results does not validate their
#' convergence, model fit or substantive interpretation.
#'
#' Graphfiles preserve the reported `Scale`, `Measure`, `Expected`, `ExpCat`
#' and `Prob:<category>` fields. Scale numbers identify Table 8 subtables, not
#' automatically items or facets. Category labels can differ from response
#' codes (including binary exports labelled 1/2); no recoding is performed.
#' `Expected` is used directly, not recomputed from probabilities or `ExpCat`.
#' The Graphfile Measure axis is not automatically a Person ability axis or
#' the coordinate system of Scorefile measures; `measure_units` and
#' `common_scale` apply only to element plots.
#'
#' Graph imports require unambiguous headings and strictly monotone distinct
#' finite measures within each scale. Missing coordinates/values are retained
#' with reasons and break curve lines. Probability sums are checked against
#' the sum of half the printed units of precision, plus floating-point tolerance;
#' sums outside that rounding allowance are withheld from probability plots.
#' Incomplete probability vectors remain flagged, with available individual
#' probabilities retained. Nothing is renormalized, smoothed or extrapolated.
#' These are source-consistency checks, not validation of the fitted model.
#' @return A `mfrm_facets_results` list with `table` (all imported rows, source
#'   tokens, IDs, status and per-quantity availability), `sources` (file paths,
#'   MD5 checksums and encoding; unavailable for data frames), `metadata` (source
#'   declarations) and `notes`. When supplied, `graph` contains the original
#'   character `table`, long `curves` with raw tokens, precision and availability,
#'   row-level `diagnostics`, `categories` and detected `layout`.
#'   Graph `metadata` identifies its separate coordinate system; its units
#'   remain unknown rather than inheriting element `measure_units`.
#'   `summary()` counts available rows by facet. `AnchorFlagged` counts
#'   `SourceStatus = "anchored"` rows, not individually fixed parameters; it
#'   replaces the earlier candidate's ambiguous `Anchored` summary column.
#'   `as.data.frame()` exports the
#'   element table. Use `component = "curves"` with either method to inspect
#'   Graphfile curves. Use `saveRDS()` to retain the full object.
#' @references FACETS Scorefile specification:
#'   <https://www.winsteps.com/facetman64/scoreoutput.htm>.
#'   FACETS Graphfile specification:
#'   <https://www.winsteps.com/facetman64/graphfile.htm>.
#'   FACETS element and Table 7 anchor annotations:
#'   <https://www.winsteps.com/facetman64/elements.htm> and
#'   <https://www.winsteps.com/facetman64/table7.htm>.
#' @examples
#' # Small synthetic binary-rating analysis exported by FACETS.
#' # Each row in result$table is one Person, Rater or Task element.
#' path <- system.file("extdata", "facets-example", package = "mfrmr")
#' result <- read_facets_results(path,
#'   facet_map = c("Person", "Rater", "Task"), source_version = "4.5.0")
#' summary(result) # Availability, not a model-adequacy assessment.
#' result$table[, c("Facet", "ElementID", "Label", "Estimate", "SourceStatus")]
#' plot(result, type = "measure", facets = "Rater", draw = FALSE)
#' @export
#' @order 1
read_facets_results <- function(file = NULL, facet = NULL, facet_map = NULL,
                                format = c("auto", "delimited", "scorefile"),
                                facet_col = NULL, level_col = NULL,
                                delimiter = NULL, encoding = "UTF-8",
                                source_version = NULL, measure_units = NULL,
                                common_scale = FALSE, graph_file = NULL) {
  format <- match.arg(format)
  .facets_results_flag(common_scale, "common_scale")
  text_option <- function(value, name) {
    if (is.null(value)) return(NA_character_)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value)))
      stop("`", name, "` must be one nonempty string or NULL.", call. = FALSE)
    value
  }
  source_version <- text_option(source_version, "source_version")
  measure_units <- text_option(measure_units, "measure_units")
  encoding <- text_option(encoding, "encoding")
  if (is.na(encoding)) stop("Supply an encoding.", call. = FALSE)
  if (is.null(file)) {
    if (is.null(graph_file)) stop("Supply `file` or `graph_file`.", call. = FALSE)
    x <- structure(list(table = data.frame(),
      sources = data.frame(SourceID = integer(), Path = character(), MD5 = character(),
        Encoding = character(), Input = character()),
      metadata = list(software = "FACETS", version = source_version,
        measure_units = measure_units, common_scale = common_scale),
      notes = data.frame(Type = character(), Text = character())),
      class = c("mfrm_facets_results", "list"))
    return(.facets_results_add_graph(x, graph_file, encoding))
  }
  if (is.data.frame(file)) {
    d <- fit_review_standardize_frame(file, facet = facet, facet_col = facet_col,
      level_col = level_col, source = "data_frame", preserve_rows = TRUE)
    if (!nrow(d)) stop("No FACETS result rows were supplied.", call. = FALSE)
    d$SourceID <- 1L
    sources <- data.frame(SourceID = 1L, Path = NA_character_, MD5 = NA_character_,
      Encoding = NA_character_, Input = "data_frame")
  } else {
    if (!is.character(file) || !length(file) || anyNA(file) || any(!nzchar(file)))
      stop("`file` must contain file paths, a score directory or a data frame.", call. = FALSE)
    if (!is.null(facet) && (!is.character(facet) || anyNA(facet) ||
        !length(facet) || any(!nzchar(facet)) || !length(facet) %in% c(1L, length(file))))
      stop("`facet` must name one facet or one facet per input path.", call. = FALSE)
    paths <- character(); facet_names <- character()
    for (i in seq_along(file)) {
      if (!file.exists(file[i])) stop("FACETS input not found: ", file[i], call. = FALSE)
      if (dir.exists(file[i])) {
        if (!is.null(facet)) stop("Use `facet_map`, not `facet`, with a directory.", call. = FALSE)
        p <- list.files(file[i], pattern = "^score[._-][0-9]+\\.txt$",
          full.names = TRUE, ignore.case = TRUE)
        if (!length(p)) stop("Directory contains no score.N.txt files.", call. = FALSE)
        f <- rep(NA_character_, length(p))
      } else {
        p <- file[i]
        f <- if (is.null(facet)) NA_character_ else facet[if (length(facet) == 1L) 1L else i]
      }
      paths <- c(paths, p); facet_names <- c(facet_names, f)
    }
    paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    if (anyDuplicated(paths)) stop("An input file was supplied more than once.", call. = FALSE)
    sources <- data.frame(SourceID = seq_along(paths), Path = paths,
      MD5 = unname(tools::md5sum(paths)), Encoding = encoding, Input = "file")
    pieces <- lapply(seq_along(paths), function(i) {
      part <- fit_review_read_one_path(paths[i],
        facet = if (is.na(facet_names[i])) NULL else facet_names[i],
        facet_map = facet_map %||% character(), format = format,
        facet_col = facet_col, level_col = level_col, delimiter = delimiter,
        encoding = encoding, preserve_rows = TRUE)
      part$SourceID <- i
      part
    })
    if (!identical(unname(tools::md5sum(paths)), sources$MD5))
      stop("A FACETS input changed during import; retry after the export finishes.", call. = FALSE)
    d <- dplyr::bind_rows(pieces)
  }
  if (!nrow(d)) stop("No FACETS result rows were supplied.", call. = FALSE)
  present <- function(z) !is.na(z) & nzchar(trimws(z))
  d$Label <- ifelse(present(d$ElementLabel), d$ElementLabel, d$Level)
  d$Label <- ifelse(present(d$Label), d$Label, d$ElementID)
  d$KeyType <- ifelse(present(d$ElementID), "element_id", "label")
  d$ElementKey <- ifelse(present(d$ElementID), d$ElementID, d$Label)
  if (any(!present(d$Facet) | !present(d$ElementKey)))
    stop("Every row needs a facet and an element ID or label; no rows were discarded.", call. = FALSE)
  if (anyDuplicated(d[c("Facet", "ElementKey")]))
    stop("Duplicate facet/element keys. Export distinct element IDs for duplicate labels, and import one analysis at a time.", call. = FALSE)
  status <- c("1" = "anchored", "0" = "no_data", "-1" = "estimated",
    "-2" = "extreme_low", "-3" = "extreme_high", "-4" = "unmeasurable",
    "-5" = "single_observation")
  d$SourceStatus <- unname(status[d$StatusRaw])
  d$SourceStatus[is.na(d$SourceStatus)] <- "unknown"
  d$Sign <- fit_review_clean_numeric(d$SignRaw)
  d$Sign[!d$Sign %in% c(-1, 0, 1)] <- NA_real_
  d$SourceVersion <- source_version
  d$MeasureUnits <- measure_units
  for (quantity in c("Estimate", "Infit", "Outfit")) {
    reason <- rep("available", nrow(d))
    reason[!is.finite(d[[quantity]])] <- "missing_or_nonfinite"
    token <- d[[paste0(quantity, "Raw")]]
    censored <- !is.na(token) & grepl("[<>]", token)
    reason[censored] <- "censored_source_value"
    if (quantity != "Estimate") {
      reason[is.finite(d[[quantity]]) & d[[quantity]] < 0] <- "invalid_mean_square"
      reason[d$SourceStatus %in% c("extreme_low", "extreme_high", "single_observation")] <-
        "source_fit_not_interpretable"
    }
    reason[d$SourceStatus %in% c("no_data", "unmeasurable")] <- "source_not_estimable"
    d[[paste0(quantity, "Available")]] <- reason == "available"
    d[[paste0(quantity, "Reason")]] <- reason
  }
  notes <- data.frame(Type = c("source", "scale", "uncertainty", "availability"),
    Text = c("Reported FACETS results; no model was fitted or validated by this import.",
      "Source units, origins and signs are preserved; larger values do not automatically mean better performance.",
      "Points only. Source marginal SEs do not supply joint covariance or calibration-inclusive intervals.",
      "Unavailable and extreme rows remain in the table. Extreme fit values can be placeholders; unknown status is not evidence of adequacy."))
  notes <- .facets_results_anchor_notes(notes)
  x <- structure(list(table = d, sources = sources,
    metadata = list(software = "FACETS", version = source_version,
      measure_units = measure_units, common_scale = common_scale), notes = notes),
    class = c("mfrm_facets_results", "list"))
  if (is.null(graph_file)) x else .facets_results_add_graph(x, graph_file, encoding)
}

.facets_results_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
}

.facets_results_anchor_notes <- function(notes) {
  # Rebuild the interpretation for saved imports too, without mutating their
  # source data or inferring A/G/X from a reporting group or a marginal SE.
  dplyr::bind_rows(notes[!notes$Type %in% "anchor_status", , drop = FALSE],
    data.frame(Type = "anchor_status", Text = paste(
      "Status 1 (SourceStatus = anchored) flags direct or group anchoring;",
      "the individual measure need not be fixed. GroupRaw, imported group size",
      "and reported SE do not identify the anchor subtype.",
      "Extreme group members can instead carry Status -2/-3; AnchorFlagged",
      "does not count every member of an anchored group. Their reported",
      "locations do not identify an extreme-score adjustment rule.",
      "Check the original Labels= specification for the anchor mode and",
      "Table 7 for row annotations; G alone need not distinguish G/X input.")))
}

#' @rdname read_facets_results
#' @param object A saved `mfrm_facets_results` object.
#' @param ... Must be empty.
#' @param component `"elements"` (default) or `"curves"` for Graphfile results.
#' @export
summary.mfrm_facets_results <- function(object, ..., component = c("elements", "curves")) {
  rlang::check_dots_empty()
  component <- match.arg(component)
  if (component == "curves") return(.facets_graph_summary(object))
  d <- object$table
  if (!nrow(d)) return(data.frame(Facet = character(), Rows = integer(),
    Measures = integer(), Infit = integer(), Outfit = integer(), AnchorFlagged = integer(),
    Extreme = integer(), UnknownStatus = integer()))
  do.call(rbind, lapply(unique(d$Facet), function(f) {
    z <- d[d$Facet == f, , drop = FALSE]
    data.frame(Facet = f, Rows = nrow(z), Measures = sum(z$EstimateAvailable),
      Infit = sum(z$InfitAvailable), Outfit = sum(z$OutfitAvailable),
      AnchorFlagged = sum(z$SourceStatus == "anchored"),
      Extreme = sum(z$SourceStatus %in% c("extreme_low", "extreme_high")),
      UnknownStatus = sum(z$SourceStatus == "unknown"))
  }))
}

#' @rdname read_facets_results
#' @export
print.mfrm_facets_results <- function(x, ...) {
  rlang::check_dots_empty()
  cat("Imported FACETS results (descriptive)\n")
  counts <- summary(x)
  print(counts, row.names = FALSE)
  if (any(counts$AnchorFlagged > 0L))
    cat("AnchorFlagged includes direct or group anchoring; the individual measure need not be fixed.\n")
  if (!is.null(x$graph)) {
    cat("Graphfile curves (reported coordinates)\n")
    print(summary(x, component = "curves"), row.names = FALSE)
  }
  invisible(x)
}

#' @rdname read_facets_results
#' @param row.names Optional row names for the exported data frame.
#' @param optional Passed to [as.data.frame()].
#' @export
as.data.frame.mfrm_facets_results <- function(x, row.names = NULL, optional = FALSE, ...,
                                            component = c("elements", "curves")) {
  rlang::check_dots_empty()
  component <- match.arg(component)
  d <- if (component == "curves") .facets_graph_curves(x) else x$table
  as.data.frame(d, row.names = row.names, optional = optional)
}

.facets_results_display_labels <- function(d, label_map) {
  labels <- d$Label
  if (is.null(label_map)) return(labels)
  fields <- c("Facet", "ElementID", "DisplayLabel")
  if (!is.data.frame(label_map) || anyDuplicated(names(label_map)) ||
      !all(fields %in% names(label_map)))
    stop("`label_map` must be a data frame with Facet, ElementID and DisplayLabel columns.", call. = FALSE)
  for (field in fields) {
    value <- label_map[[field]]
    if (!is.character(value) || anyNA(value) || any(!nzchar(trimws(value))))
      stop("`label_map$", field, "` must contain nonempty character values without NA; preserve IDs as text.", call. = FALSE)
  }
  if (anyDuplicated(label_map[c("Facet", "ElementID")]))
    stop("`label_map` contains duplicate Facet/ElementID keys.", call. = FALSE)
  # Match within facets, without concatenating keys or using source labels.
  for (facet in unique(label_map$Facet)) {
    target <- which(d$Facet == facet & !is.na(d$ElementID) & nzchar(d$ElementID))
    supplied <- which(label_map$Facet == facet)
    at <- match(label_map$ElementID[supplied], d$ElementID[target])
    if (anyNA(at))
      stop("Every `label_map` key must match an exported Facet/ElementID in the result.", call. = FALSE)
    labels[target[at]] <- label_map$DisplayLabel[supplied]
  }
  labels
}

#' @rdname read_facets_results
#' @details `plot()` method: Plot imported FACETS measures, fit statistics and curves
#' @param x Saved output of [read_facets_results()].
#'
#' Also for `summary()`: A saved `mfrm_facets_results` object.
#' @param type `"measure"` for labelled element locations by facet; `"wright"`
#'   for points by facet on one measure axis; `"fit"` for measure versus mean
#'   square; or `"infit_outfit"` for Infit versus Outfit. Wright maps require
#'   the explicit `common_scale = TRUE` declaration at import.
#'   `"ccc"` plots Graphfile category probabilities and `"expected"` plots
#'   Graphfile expected scores, each in separate scale panels.
#' @param facets Optional character vector of facet names to display.
#' @param scales Optional vector of Graphfile scale numbers for curve plots.
#'   Facet selectors do not apply to curves, and scale selectors do not apply
#'   to element plots.
#' @param statistic `"infit"` or `"outfit"`, used for `type = "fit"`.
#' @param draw Logical; draw with base graphics when true.
#' @param show_title,show_notes Logical; show title/interpretation notes in the
#'   plot, off by default. Notes always remain in the returned data.
#' @param show_legend Logical; display the source-status symbol key or category
#'   curve key below the panels. Defaults to true.
#' @param label_map Optional data frame with character columns `Facet`,
#'   `ElementID` and `DisplayLabel` for element plots. Keys must match exported
#'   IDs in `x`, including any leading zeros, and must be unique within a facet.
#'   A partial map is allowed; unmapped rows retain their source labels.
#'   Mapping is checked against the full result before `facets` selection.
#'   This changes display labels only, not source labels or element identity.
#'   Curve plots do not accept this map. Only `type = "measure"` draws element
#'   labels; the other element plots retain them in returned plot data.
#' @return For `plot()`: Invisibly an `mfrm_plot_data` object. `data$plot_data` retains every
#'   selected row with `X`, `Y`, `Available` and `Reason`; element plots also
#'   have `Symbol`, and curves have `Series`, `Color`, `LineType` and `Segment`.
#'   Unavailable
#'   numerical coordinates are `NA` (measure plots retain label positions).
#'   Original reported values/status remain in the same
#'   rows. `data$notes`, `data$sources` and `data$metadata` describe interpretation
#'   and origin. `data$label_map` retains a supplied display map. Filtering and
#'   display labels do not change the saved results.
#' @details For `plot()`: Each measure panel has its own axis; the Wright display shares one measure
#' axis without centering, reversing signs or inventing category thresholds.
#' Filled circles denote estimated elements, squares direct/group anchor flags, downward/upward
#' triangles extreme low/high source status, diamonds single observations and
#' open circles unknown source status. The symbol key is in `data$legend`.
#' A square marks source Status 1; it does not identify an individually fixed
#' measure. Group-anchor members can vary around a fixed group mean, and
#' `GroupRaw` may instead be a reporting group. Check the original specification
#' or Table 7 to distinguish direct and group anchors. Returned notes retain
#' this limitation, also when plotting imports saved by an earlier version;
#' the source object is not changed. Recreate saved plot-data objects with
#' `plot()` to obtain the revised legend and notes.
#' Duplicate display labels are disambiguated by element key in measure plots.
#' Use `label_map` to display local names while retaining the original FACETS
#' labels and IDs. It cannot recover names lost in an export; supply the known
#' correspondence yourself. Character preservation and font rendering are
#' separate: choose an installed font with the required glyphs on the graphics
#' device. With Cairo support, [grDevices::png()] with `type = "cairo"` and
#' [grDevices::cairo_pdf()] accept a `family` such as `"Meiryo"` on Windows
#' when installed. Font availability differs across systems. The plot does not
#' install fonts or change global font settings.
#' Infit/Outfit reference lines at one describe the model reference, not a
#' fairness or acceptability cutoff. No confidence bars are calculated.
#' Extreme/single-observation fit placeholders and unmeasurable/no-data elements
#' are retained as unavailable. A plot can have no usable points; inspect its
#' row-level reasons. Unknown status is displayed explicitly, not certified.
#' Graphfile curves join consecutive exported coordinates with straight line
#' segments; this does not calculate new model predictions. Missing or invalid
#' rows break the lines. The legend identifies source category labels. Expected
#' scores and category labels retain their source coding. Each panel has its
#' own reported Measure axis, independent of element-plot scale declarations.
#' For curves, `data$data` holds the raw Graphfile table, `data$plot_data` the
#' selected long rows, and `data$diagnostics` the probability-sum checks.
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_facets_results", run.dontrun = TRUE)
#' \dontrun{
#' path <- system.file("extdata", "facets-example", package = "mfrmr")
#' result <- read_facets_results(path,
#'   facet_map = c("Person", "Rater", "Task"), source_version = "4.5.0")
#' plot(result, type = "measure", facets = "Person")
#' p <- plot(result, type = "infit_outfit", facets = "Rater", draw = FALSE)
#' p$data$plot_data[, c("Label", "X", "Y", "Available", "Reason")]
#' p$data$notes
#' # Display names are joined by facet and exported ID, never by row order.
#' labels <- data.frame(Facet = "Person", ElementID = c("1", "2"),
#'   DisplayLabel = c("Participant A", "Participant B"))
#' labelled <- plot(result, facets = "Person", label_map = labels, draw = FALSE)
#' labelled$data$plot_data[, c("ElementID", "Label", "DisplayLabel")]
#' # A separate synthetic PCM analysis, with four exported scale curves.
#' graph_path <- system.file("extdata", "facets-graph-example", "graph.txt",
#'   package = "mfrmr")
#' curves <- read_facets_results(graph_file = graph_path, source_version = "4.5.0")
#' summary(curves, component = "curves")
#' plot(curves, type = "ccc", scales = 1)
#' expected <- plot_data(curves, "plot_data", type = "expected", scales = 1)
#' head(expected[, c("Scale", "Measure", "Value", "Available")])
#' }
#' @export
plot.mfrm_facets_results <- function(x, type = c("measure", "wright", "fit", "infit_outfit", "ccc", "expected"),
                                    facets = NULL, statistic = c("infit", "outfit"),
                                    draw = TRUE, show_title = FALSE, show_notes = FALSE,
                                    show_legend = TRUE, scales = NULL, label_map = NULL, ...) {
  rlang::check_dots_empty()
  type <- match.arg(type); statistic <- match.arg(statistic)
  for (name in c("draw", "show_title", "show_notes", "show_legend")) .facets_results_flag(get(name), name)
  if (type %in% c("ccc", "expected")) {
    if (!is.null(facets)) stop("Use `scales`, not `facets`, for Graphfile curves.", call. = FALSE)
    if (!is.null(label_map)) stop("`label_map` applies only to element plots.", call. = FALSE)
    return(invisible(.facets_graph_plot(x, type, scales, draw, show_title, show_notes, show_legend)))
  }
  if (!is.null(scales)) stop("`scales` applies only to Graphfile curves.", call. = FALSE)
  if (!nrow(x$table)) stop("No element results were imported; select a Graphfile curve type.", call. = FALSE)
  if (type == "wright" && !isTRUE(x$metadata$common_scale))
    stop("A Wright map requires compatible source units and origins. Verify them and import with `common_scale = TRUE`.", call. = FALSE)
  d <- as.data.frame(x$table)
  display_labels <- .facets_results_display_labels(d, label_map)
  if (is.null(facets)) facets <- unique(d$Facet)
  if (!is.character(facets) || !length(facets) || anyNA(facets) ||
      anyDuplicated(facets) || any(!facets %in% d$Facet))
    stop("`facets` must select distinct facet names in the result.", call. = FALSE)
  selected <- d$Facet %in% facets
  display_labels <- display_labels[selected]
  d <- d[selected, , drop = FALSE]
  stat <- if (statistic == "infit") "Infit" else "Outfit"
  quantities <- switch(type, measure = "Estimate", wright = "Estimate",
    fit = c("Estimate", stat), infit_outfit = c("Infit", "Outfit"))
  reasons <- lapply(quantities, function(q) ifelse(d[[paste0(q, "Available")]], "",
    paste0(q, ":", d[[paste0(q, "Reason")]])))
  d$Reason <- vapply(seq_len(nrow(d)), function(i) {
    z <- vapply(reasons, `[`, character(1), i)
    if (any(nzchar(z))) paste(z[nzchar(z)], collapse = "; ") else "available"
  }, character(1))
  d$Available <- d$Reason == "available"
  d$Panel <- if (type == "wright") "Wright map" else d$Facet
  d$X <- switch(type, measure = d$Estimate, wright = as.numeric(match(d$Facet, facets)),
    fit = d$Estimate, infit_outfit = d$Infit)
  d$Y <- switch(type, measure = as.numeric(ave(seq_len(nrow(d)), d$Facet,
    FUN = function(z) rev(seq_along(z)))), wright = d$Estimate,
    fit = d[[stat]], infit_outfit = d$Outfit)
  d$X[!d$Available] <- NA_real_
  if (type != "measure") d$Y[!d$Available] <- NA_real_
  d$DisplayLabel <- display_labels
  repeat {
    dup <- duplicated(d[c("Facet", "DisplayLabel")]) |
      duplicated(d[c("Facet", "DisplayLabel")], fromLast = TRUE)
    if (!any(dup)) break
    d$DisplayLabel[dup] <- paste0(d$DisplayLabel[dup], " [", d$ElementKey[dup], "]")
  }
  symbols <- c(estimated = 19L, anchored = 15L, extreme_low = 25L,
    extreme_high = 24L, single_observation = 5L, unknown = 1L,
    no_data = 4L, unmeasurable = 4L)
  d$Symbol <- unname(symbols[d$SourceStatus])
  symbol_labels <- c("Estimated", "Direct/group anchor", "Extreme low", "Extreme high",
    "Single observation", "Unknown status", "No data", "Unmeasurable")
  used_symbols <- which(names(symbols) %in% d$SourceStatus[d$Available])
  units <- x$metadata$measure_units
  measure_label <- if (is.na(units)) "Reported measure" else paste0("Measure (", units, ")")
  xlab <- switch(type, measure = measure_label, wright = "Facet", fit = measure_label,
    infit_outfit = "Infit mean square")
  ylab <- switch(type, measure = "", wright = measure_label,
    fit = paste0(stat, " mean square"), infit_outfit = "Outfit mean square")
  out <- new_mfrm_plot_data(paste0("facets_results_", type), list(
    data = as.data.frame(x$table), plot_data = d, sources = x$sources,
    metadata = x$metadata, notes = .facets_results_anchor_notes(x$notes), label_map = label_map,
    title = "Imported FACETS results",
    xlab = xlab, ylab = ylab,
    legend = new_plot_legend(label = symbol_labels, role = rep("source_status", length(symbols)),
      aesthetic = rep("pch", length(symbols)), value = unname(symbols)),
    reference_lines = if (type %in% c("fit", "infit_outfit")) new_reference_lines(
      axis = if (type == "fit") "y" else c("x", "y"), value = 1,
      label = "Model reference", linetype = "dashed", role = "mean_square_reference") else NULL,
    display = list(show_title = show_title, show_notes = show_notes, show_legend = show_legend)))
  if (draw) {
    panels <- if (type == "wright") "Wright map" else facets
    nc <- min(2L, length(panels)); nr <- ceiling(length(panels) / nc)
    legend_space <- if (show_legend && length(used_symbols)) 3 else 0
    old <- graphics::par(mfrow = c(nr, nc), mar = c(4.5, if (type == "measure") 6 else 4.5,
      if (show_title) 3 else 1, 1), oma = c(if (show_notes) 3 else 0, 0, 0, 0), xpd = FALSE)
    on.exit(graphics::par(old), add = TRUE)
    if (legend_space > 0) {
      # A separate region survives recorded-plot replay at a new device size.
      # Converting device coordinates to user units freezes the old margins.
      panel_cex <- graphics::par("cex")
      panel_layout <- matrix(c(seq_along(panels), rep(0L, nr * nc - length(panels))),
        nrow = nr, byrow = TRUE)
      graphics::layout(rbind(panel_layout, rep(length(panels) + 1L, nc)),
        heights = c(rep(1, nr), graphics::lcm(2.54 * graphics::par("csi") * legend_space)))
      graphics::par(cex = panel_cex)
    }
    limits <- function(z) {
      z <- z[is.finite(z)]
      if (!length(z)) return(c(0, 1))
      r <- range(z)
      if (diff(r) == 0) r <- r + c(-1, 1) * max(.1, abs(r[1]) * .05)
      r
    }
    for (panel in panels) {
      z <- d[d$Panel == panel, , drop = FALSE]
      graphics::plot(NA_real_, NA_real_, xlim = if (type == "wright") c(.5, length(facets) + .5) else limits(z$X),
        ylim = if (type == "measure") c(.5, nrow(z) + .5) else limits(z$Y),
        xlab = if (type == "wright") xlab else paste0(xlab, " [", panel, "]"),
        ylab = ylab, xaxt = if (type == "wright") "n" else "s",
        yaxt = if (type == "measure") "n" else "s",
        main = if (show_title) out$data$title else "")
      if (type == "measure") graphics::axis(2, at = z$Y, labels = z$DisplayLabel, las = 2, cex.axis = .8)
      if (type == "wright") graphics::axis(1, at = seq_along(facets), labels = facets)
      if (type %in% c("fit", "infit_outfit")) graphics::abline(h = 1, col = "grey60", lty = 2)
      if (type == "infit_outfit") graphics::abline(v = 1, col = "grey60", lty = 2)
      graphics::points(z$X, z$Y, pch = z$Symbol, bg = "white")
      if (!any(z$Available)) graphics::text(mean(graphics::par("usr")[1:2]),
        mean(graphics::par("usr")[3:4]), "No available points")
    }
    if (legend_space > 0) {
      graphics::par(mar = rep(0, 4))
      graphics::plot.new()
      graphics::legend("center", legend = symbol_labels[used_symbols], pch = unname(symbols[used_symbols]),
        pt.bg = "white", ncol = min(3L, length(used_symbols)), cex = .7, bty = "n", xpd = NA)
    }
    if (show_notes) graphics::mtext("Reported FACETS points; no confidence intervals. Inspect saved notes and source scales.",
      side = 1, outer = TRUE, line = .5, cex = .7)
  }
  invisible(out)
}
