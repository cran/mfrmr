# Shared numerical basis for the Q3-style view and residual exploration.
.mfrm_residual_associations <- function(fit, diagnostics, facet, min_pairs,
                                        missing = "pairwise") {
  if (!inherits(fit, "mfrm_fit")) stop("`fit` must be an mfrm_fit.", call. = FALSE)
  if (!is.character(facet) || length(facet) != 1L || is.na(facet) ||
      !facet %in% fit$config$facet_names)
    stop("`facet` must be one of: ", paste(fit$config$facet_names, collapse = ", "),
      "; select a single fitted facet.", call. = FALSE)
  if (!is.numeric(min_pairs) || is.complex(min_pairs) || length(min_pairs) != 1L ||
      !is.finite(min_pairs) || min_pairs < 3 || min_pairs != floor(min_pairs))
    stop("`min_pairs` must be a single integer of at least 3.", call. = FALSE)
  obs <- if (is.null(diagnostics)) {
    compute_obs_table(.mfrm_validate_diagnostic_fit(fit))
  } else {
    mfrm_results_validate_diagnostics_identity(fit, diagnostics,
      helper = "Residual associations")
    diagnostics$obs
  }
  obs <- as.data.frame(obs)
  if (!all(c("Person", facet, "StdResidual") %in% names(obs)))
    stop("Diagnostics require Person, the selected facet and StdResidual.", call. = FALSE)
  valid <- !is.na(obs$Person) & !is.na(obs[[facet]])
  persons <- sort(unique(as.character(obs$Person[valid])))
  levels <- sort(unique(as.character(obs[[facet]][valid])))
  if (length(levels) < 2L) stop("At least two distinct facet levels are required.", call. = FALSE)
  # Integer cell indices avoid collisions between arbitrary user labels.
  i <- match(as.character(obs$Person[valid]), persons)
  j <- match(as.character(obs[[facet]][valid]), levels)
  values <- obs$StdResidual[valid]
  values[!is.finite(values)] <- NA_real_
  cell <- i + (j - 1L) * length(persons)
  means <- tapply(values, cell, mean, na.rm = TRUE)
  residuals <- matrix(NA_real_, length(persons), length(levels),
    dimnames = list(persons, levels))
  residuals[as.integer(names(means))] <- means
  residuals[!is.finite(residuals)] <- NA_real_
  counts <- matrix(0L, length(persons), length(levels), dimnames = dimnames(residuals))
  finite_counts <- tapply(is.finite(values), cell, sum)
  counts[as.integer(names(finite_counts))] <- finite_counts
  used <- if (missing == "complete") stats::complete.cases(residuals) else
    rowSums(is.finite(residuals)) > 0L
  r <- residuals[used, , drop = FALSE]
  mat <- matrix(NA_real_, length(levels), length(levels), dimnames = list(levels, levels))
  diag(mat) <- vapply(seq_along(levels), function(k) {
    a <- r[, k]; a <- a[is.finite(a)]
    if (length(a) >= min_pairs && stats::sd(a) > 0) 1 else NA_real_
  }, numeric(1))
  pair <- utils::combn(seq_along(levels), 2L)
  pairs <- do.call(rbind, lapply(seq_len(ncol(pair)), function(k) {
    a <- r[, pair[1L, k]]; b <- r[, pair[2L, k]]
    ok <- is.finite(a) & is.finite(b); n <- sum(ok)
    value <- if (n >= min_pairs) suppressWarnings(stats::cor(a[ok], b[ok])) else NA_real_
    if (is.finite(value)) value <- max(-1, min(1, value))
    mat[pair[1L, k], pair[2L, k]] <<- mat[pair[2L, k], pair[1L, k]] <<- value
    data.frame(Level1 = levels[pair[1L, k]], Level2 = levels[pair[2L, k]],
      ResidualCor = value, N = n, Reason = if (n < min_pairs)
        "Insufficient shared persons" else if (!is.finite(value))
        "No residual variation in the shared persons" else "", stringsAsFactors = FALSE)
  }))
  list(matrix = mat, pairs = pairs, residuals = residuals, cell_counts = counts,
    persons_used = used, excluded_identifier_rows = sum(!valid),
    nonfinite_residual_rows = sum(!is.finite(obs$StdResidual)), input_rows = nrow(obs))
}

#' Explore residual associations and groups of facet levels
#'
#' After reviewing a fitted measurement model, locate facet levels whose
#' standardized residuals move together across persons. For example, explore
#' criteria that share unexplained variation, or raters with similar residual
#' profiles. This descriptive, signed correlation network complements
#' [q3_statistic()], residual PCA and case review.
#'
#' @param fit A supported [fit_mfrm()] result, with ordinary observation-level
#'   standardized residuals available from [diagnose_mfrm()].
#' @param diagnostics Optional matching diagnostics; computed once if omitted.
#' @param facet A single fitted facet. Its levels, not persons, are the nodes
#'   and the objects clustered. Defaults to `"Rater"`.
#' @param min_pairs Minimum shared persons for a correlation, at least three.
#' @param missing `"pairwise"` uses the available persons for each pair.
#'   `"complete"` uses only persons with finite aggregated residuals at every
#'   selected facet level. No residual is imputed. Clustering requires
#'   `"complete"` so all distances describe the same sample.
#' @param min_abs Absolute correlation threshold for displayed edges, in
#'   `[0, 1]`; default zero retains every nonzero available association.
#'   This is a display choice, not a significance threshold. It never changes
#'   the correlations or the clustering distance.
#' @param clusters Optional integer number of groups from 2 to the number of
#'   levels. `NULL` skips clustering. The number is supplied, not estimated.
#' @param linkage `"average"` (default) or `"complete"` agglomerative linkage.
#'
#' @section Statistical definition:
#' Finite `StdResidual` values are averaged without weighting within each
#' Person-by-facet-level cell. Pearson correlations across persons define the
#' network. Thus repeated ratings are not treated as independent persons.
#' `cell_counts` records differing amounts of evidence behind those means.
#' Missing overlap and zero variance produce unavailable (`NA`) associations,
#' with reasons in `pairs`; they never become zero edges. All nodes, including
#' isolates, remain in the output. Positive edges indicate co-movement and
#' negative edges counter-movement after the fitted model. These are marginal
#' residual correlations, not conditional/partial correlations or raw Yen Q3.
#'
#' Unequal assignments, aggregation, estimated Person effects and finite
#' samples can induce residual structure, including negative correlations.
#' Inspect overlap and the rating design before explaining an edge as a
#' substantive relationship. Threshold sensitivity is not a calibrated test.
#' No edge p-values, uncertainty intervals or stability estimates are supplied.
#'
#' @section Machine learning and interpretation:
#' Hierarchical agglomerative clustering is unsupervised learning: it merges
#' similar nodes without known outcome labels. With complete-person Pearson
#' correlations, `sqrt(2 * (1 - r))` is Euclidean distance between centered,
#' unit-length residual profiles. Average linkage uses the mean between-group
#' distance; complete linkage uses the maximum. Strong negative correlation
#' means dissimilar profiles, not similarity based on absolute correlation.
#' A supplied `clusters` cuts the resulting [stats::hclust()] tree with
#' [stats::cutree()]. Changing `min_abs` cannot change these groups.
#'
#' Clustering remains unavailable if any level lacks sufficient common-person
#' variation; the graph and failure reason are still returned. Pairwise
#' deletion can yield an indefinite correlation matrix, so it is not used
#' to manufacture clustering distances. Common-person deletion can select a
#' different population and greatly reduce the sample; inspect `person_usage`.
#'
#' Groups are exploratory residual-profile similarities, not latent classes,
#' demonstrated rubric dimensions, rater competence or automatic grounds for
#' excluding a rater. This procedure is not exploratory graph analysis (EGA):
#' Golino and Epskamp (2017) used a graphical-LASSO partial-correlation network
#' and community detection. It is also not a jointly estimated residual-network
#' measurement model. Compare substantive explanations, inspect cases and
#' check sensitivity to sample, linkage and group count. For validation, repeat
#' preprocessing, measurement fitting and clustering within independent or
#' appropriately resampled data; resampling saved residuals alone omits
#' calibration uncertainty. The appropriate resampling unit follows the rating
#' design, including dependencies across persons, tasks and raters.
#'
#' @section Where random forests belong:
#' Wesolowski (2019) first used k-means on differential-rater-functioning
#' indices to construct rater types, then random forests to predict those
#' labels. That study motivates a supervised follow-up, not random forests as
#' the clustering algorithm implemented here. Its features and targets differ
#' from these residual correlations. In this supervised use, a forest combines
#' predictions from bootstrap samples and random candidate feature subsets;
#' it requires a declared prediction target and out-of-sample assessment.
#'
#' A classifier predicting exploratory groups can reproduce the grouping rule
#' without establishing external validity. Fit feature transformations,
#' measurement models and cluster discovery within training data, keep a
#' separate test set and define how new units receive reference labels. Split
#' at the intended generalization unit (for example new raters, persons or
#' tasks); row-wise splitting of repeated ratings can leak information.
#' Out-of-bag error alone does not validate labels constructed on all data.
#' Report class balance, confusion matrices, discrimination/calibration where
#' appropriate and independent substantive outcomes. Feature importance is
#' predictive, not causal. This function does not fit a random forest; use
#' established ML packages for a separately designed prediction study.
#'
#' @return An `mfrm_residual_network` with `nodes`, all `pairs`, displayed
#'   `edges`, the correlation `matrix`, aggregated `residuals`, `cell_counts`,
#'   `person_usage`, `settings`, `clustering` (status, distance, tree and
#'   membership), source model metadata and interpretation `notes`.
#'   [plot()] returns saved figure data; [mfrm_report()] and
#'   [build_summary_table_bundle()] connect the result to reporting.
#'   Complete RDS archives include person identifiers: review sharing permissions.
#' @references
#' Wesolowski, B. C. (2019). Predicting operational rater-type classifications
#' using Rasch measurement theory and random forests: A music performance
#' assessment perspective. *Journal of Educational Measurement, 56*, 610--625.
#' \doi{10.1111/jedm.12227}.
#'
#' Golino, H. F., & Epskamp, S. (2017). Exploratory graph analysis: A new
#' approach for estimating the number of dimensions in psychological research.
#' *PLOS ONE, 12*, e0174035. \doi{10.1371/journal.pone.0174035}.
#' @seealso [mfrm_cluster_kmeans()], [mfrm_cluster_hierarchical()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("mfrm_residual_network", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' net <- mfrm_residual_network(fit, facet = "Criterion", missing = "complete",
#'                              clusters = 2)
#' net$pairs                 # association, shared-person count, unavailability
#' net$clustering$status     # check before interpreting groups
#' plot(net, type = "matrix", draw = FALSE)
#' report <- mfrm_report(net)
#' build_summary_table_bundle(report)$table_index
#' }
#' @export
#' @order 1
mfrm_residual_network <- function(fit, diagnostics = NULL, facet = "Rater",
    min_pairs = 5L, missing = c("pairwise", "complete"), min_abs = 0,
    clusters = NULL, linkage = c("average", "complete")) {
  missing <- match.arg(missing); linkage <- match.arg(linkage)
  if (!is.numeric(min_abs) || is.complex(min_abs) || length(min_abs) != 1L || !is.finite(min_abs) ||
      min_abs < 0 || min_abs > 1) stop("`min_abs` must lie in [0, 1].", call. = FALSE)
  if (!is.null(clusters) && (!is.numeric(clusters) || is.complex(clusters) || length(clusters) != 1L ||
      !is.finite(clusters) || clusters < 2 || clusters != floor(clusters)))
    stop("`clusters` must be an integer of at least 2 or NULL.", call. = FALSE)
  if (!is.null(clusters) && missing != "complete")
    stop("Clustering requires explicit `missing = 'complete'` for common-person distances.", call. = FALSE)
  a <- .mfrm_residual_associations(fit, diagnostics, facet, min_pairs, missing)
  if (!is.null(clusters) && clusters > ncol(a$matrix))
    stop("`clusters` cannot exceed the number of facet levels.", call. = FALSE)
  a$pairs$Displayed <- is.finite(a$pairs$ResidualCor) &
    abs(a$pairs$ResidualCor) > 0 & abs(a$pairs$ResidualCor) >= min_abs
  edges <- a$pairs[a$pairs$Displayed, , drop = FALSE]
  groups <- list(status = "not_requested", reason = "No cluster count supplied.",
    distance = NULL, tree = NULL, membership = NULL)
  if (!is.null(clusters)) {
    groups$status <- "unavailable"
    groups$reason <- "Insufficient common persons or no common-person residual variation at one or more levels."
    if (all(is.finite(a$matrix))) {
      groups$distance <- stats::as.dist(sqrt(pmax(2 * (1 - a$matrix), 0)))
      groups$tree <- stats::hclust(groups$distance, method = linkage)
      groups$membership <- stats::cutree(groups$tree, k = clusters)
      groups$status <- "available"; groups$reason <- "Exploratory groups; stability is not estimated."
    }
  }
  nodes <- data.frame(Level = colnames(a$matrix),
    Persons = colSums(is.finite(a$residuals[a$persons_used, , drop = FALSE])),
    DisplayedDegree = as.integer(table(factor(c(edges$Level1, edges$Level2),
      levels = colnames(a$matrix)))),
    Cluster = if (is.null(groups$membership)) NA_integer_ else unname(groups$membership),
    row.names = NULL)
  notes <- c(
    "Standardized residuals are averaged without weighting within Person/facet-level cells, then correlated across persons. Inspect cell_counts and pair overlap.",
    "Signed marginal residual correlations are descriptive; unavailable pairs are NA, not zero. Edge thresholds do not establish significance or local independence.",
    "Exploratory groups describe facet-level residual profiles, not confirmed dimensions, latent classes or rater quality. No stability or calibration uncertainty is estimated.",
    "Shared estimated Person effects, aggregation and assignments can induce residual associations. Review fit and individual cases before substantive interpretation.",
    "Complete-person deletion changes the analyzed sample. Cluster distances use sqrt(2*(1-r)); edge thresholds never determine membership.",
    "Higher coded scores need not mean better performance. Record rubric direction; a residual sign is relative to model expectation, not a quality label.")
  structure(list(nodes = nodes, pairs = a$pairs, edges = edges, matrix = a$matrix,
    residuals = a$residuals, cell_counts = a$cell_counts,
    person_usage = data.frame(Person = rownames(a$residuals), Used = a$persons_used,
      FiniteLevels = rowSums(is.finite(a$residuals)), row.names = NULL),
    settings = data.frame(Facet = facet, Missing = missing, MinPairs = min_pairs,
      MinAbs = min_abs, Clusters = clusters %||% NA_integer_, Linkage = linkage,
      Persons = nrow(a$residuals), UsedPersons = sum(a$persons_used),
      InputRows = a$input_rows, ExcludedIdentifierRows = a$excluded_identifier_rows,
      NonfiniteResidualRows = a$nonfinite_residual_rows),
    source = list(model = fit$config$model, method = fit$config$method,
      facets = fit$config$facet_names, facet_signs = fit$config$facet_signs),
    clustering = groups, notes = notes), class = "mfrm_residual_network")
}

#' @rdname mfrm_residual_network
#' @param object An `mfrm_residual_network`.
#' @param ... Reserved; unused arguments are rejected.
#' @export
summary.mfrm_residual_network <- function(object, ...) {
  rlang::check_dots_empty()
  list(settings = object$settings, nodes = object$nodes,
    pair_status = data.frame(Pairs = nrow(object$pairs),
      Available = sum(is.finite(object$pairs$ResidualCor)), Displayed = nrow(object$edges)),
    clustering = data.frame(Status = object$clustering$status, Reason = object$clustering$reason),
    notes = object$notes)
}

#' @rdname mfrm_residual_network
#' @param x For `print.mfrm_residual_network()`: An `mfrm_residual_network`.
#'
#' For `plot()`: A result from [mfrm_residual_network()].
#' @export
print.mfrm_residual_network <- function(x, ...) {
  rlang::check_dots_empty()
  cat("Residual association network (descriptive)\n")
  print(x$settings, row.names = FALSE)
  cat(nrow(x$nodes), "nodes;", nrow(x$edges), "displayed edges;",
    sum(!is.finite(x$pairs$ResidualCor)), "unavailable pairs\n")
  cat("Clustering:", x$clustering$status, "-", x$clustering$reason, "\n")
  invisible(x)
}

#' @rdname mfrm_residual_network
#' @details `plot()` method: Plot saved residual associations and exploratory groups
#'
#' Compare the full matrix with its thresholded network, then inspect the
#' clustering tree when requested. Network positions are a deterministic
#' circle: distance on the page is not a statistical distance. Positive edges
#' are solid red and negative edges dashed blue; thickness represents |r|.
#' Gray matrix cells denote unavailable correlations. Node colors do not
#' encode quality. Dendrogram labels append the chosen group number in square
#' brackets; these numbers are labels, not ranks. Titles and cautions are retained in [plot_data()], not
#' printed into the figure.
#' @param type `"network"`, `"matrix"` or `"dendrogram"`.
#' @param preset Visual preset; see [mfrmr_visual_diagnostics].
#' @param draw Whether to draw with base graphics.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` containing the selected view, saved
#'   matrix, nodes, pairs, edges, clustering, settings and interpretation notes.
#'   Plotting does not recompute residuals or clustering. Save the analysis
#'   object to redraw every view. The matrix supports [as_ggplot()].
#' @export
plot.mfrm_residual_network <- function(x, type = c("network", "matrix", "dendrogram"),
    preset = c("standard", "publication", "compact", "monochrome"), draw = TRUE, ...) {
  rlang::check_dots_empty(); type <- match.arg(type)
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  style <- resolve_plot_preset(preset)
  if (type == "dendrogram" && is.null(x$clustering$tree))
    stop("No clustering tree is available: ", x$clustering$reason, call. = FALSE)
  n <- nrow(x$nodes); angle <- 2 * pi * (seq_len(n) - 1L) / n
  layout <- data.frame(Level = x$nodes$Level, X = cos(angle), Y = sin(angle))
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old <- graphics::par(c("mar", "xpd")); on.exit(graphics::par(old), add = TRUE)
    if (type == "dendrogram") {
      graphics::par(mar = c(4, 2, 1, 10))
      tree <- stats::dendrapply(stats::as.dendrogram(x$clustering$tree), function(node) {
        if (stats::is.leaf(node)) {
          label <- attr(node, "label")
          attr(node, "label") <- paste0(label, " [", x$clustering$membership[[label]], "]")
        }
        node
      })
      graphics::plot(tree, horiz = TRUE, main = "", axes = FALSE,
        nodePar = list(lab.cex = .85, pch = NA),
        xlab = paste(x$settings$Linkage, "linkage distance"))
      graphics::axis(1)
    } else if (type == "matrix") {
      graphics::par(mar = c(6, 6, 3, 2))
      graphics::plot(NA, xlim = c(.5, n + .5), ylim = c(.5, n + .5),
        xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i", asp = 1)
      palette <- grDevices::hcl.colors(101, "RdBu", rev = TRUE)
      for (i in seq_len(n)) for (j in seq_len(n)) {
        r <- x$matrix[i, j]
        col <- if (is.finite(r)) palette[1L + round((r + 1) * 50)] else "grey80"
        graphics::rect(i - .5, j - .5, i + .5, j + .5, col = col, border = "white")
        if (n <= 12L && is.finite(r)) graphics::text(i, j, sprintf("%.2f", r),
          col = if (abs(r) > .65) "white" else "black")
      }
      graphics::axis(1, seq_len(n), x$nodes$Level, las = 2)
      graphics::axis(2, seq_len(n), x$nodes$Level, las = 1)
      graphics::legend("top", inset = c(0, -.12), legend = c("-1", "0", "+1", "Unavailable"),
        fill = c(palette[c(1, 51, 101)], "grey80"), horiz = TRUE,
        cex = .75, bty = "n", xpd = NA)
    } else {
      graphics::par(mar = c(1, 1, 2, 1))
      graphics::plot(layout$X, layout$Y, type = "n", axes = FALSE, xlab = "", ylab = "",
        xlim = c(-1.5, 1.5), ylim = c(-1.4, 1.4), asp = 1)
      for (i in seq_len(nrow(x$edges))) {
        a <- match(x$edges$Level1[i], layout$Level); b <- match(x$edges$Level2[i], layout$Level)
        r <- x$edges$ResidualCor[i]
        color <- if (style$name == "monochrome") "grey35" else if (r > 0) "#B2182B" else "#2166AC"
        graphics::segments(layout$X[a], layout$Y[a], layout$X[b], layout$Y[b],
          col = color, lty = if (r > 0) 1 else 2, lwd = 1 + 4 * abs(r))
      }
      graphics::points(layout$X, layout$Y, pch = 21, bg = "white", cex = 2)
      for (i in seq_len(n)) graphics::text(layout$X[i] * 1.12, layout$Y[i] * 1.12,
        layout$Level[i], cex = .85, adj = c(
          if (layout$X[i] > .7) 0 else if (layout$X[i] < -.7) 1 else .5,
          if (layout$Y[i] > .7) 0 else if (layout$Y[i] < -.7) 1 else .5))
      graphics::legend("top", legend = c("Positive", "Negative"), lty = c(1, 2),
        col = if (style$name == "monochrome") rep("grey35", 2) else c("#B2182B", "#2166AC"),
        bty = "n", horiz = TRUE)
    }
  }
  invisible(new_mfrm_plot_data(if (type == "matrix") "local_dependence_heatmap" else "residual_network",
    list(type = type, matrix = x$matrix, nodes = x$nodes, pairs = x$pairs,
      edges = x$edges, layout = layout, clustering = x$clustering,
      settings = x$settings, facet = x$settings$Facet, notes = x$notes,
      title = "Residual association exploration", preset = style$name)))
}
