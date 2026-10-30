# Saved planning results: describe the decision without rerunning simulations.
.mfrm_design_target <- function(facets) {
  data.frame(Facet = facets,
    Target = ifelse(facets == "Person", "Person differentiation", "Facet calibration"),
    Interpretation = ifelse(facets == "Person",
      "Distinguishing Persons in the simulated population; not classification accuracy at a cut score.",
      "Distinguishing facet-level locations; high separation reliability does not mean high inter-rater agreement."),
    stringsAsFactors = FALSE)
}

#' @rdname recommend_mfrm_design
#' @details `summary()` method: Review a simulation-based design recommendation
#' @param object Output from [recommend_mfrm_design()].
#' @details For `summary()`: Recommendations screen the saved simulation aggregates under the requested
#' criteria. The default targets non-Person facet calibration. For Person
#' differentiation, request `facets = "Person"` when recommending the design.
#' Neither target establishes classification accuracy, inter-rater agreement,
#' a universal minimum sample size or precision outside the simulated population.
#'
#' `checks` has one row per design, target facet and criterion, including
#' unavailable metrics and inactive constraints. `DesignRow` links these rows
#' to the complete, ranked `design_table`. `Selected` marks the first passing
#' plan under the stored `prefer` order; it does not establish an optimum over
#' unevaluated plans. `Pass` and selected plans are reused without recalculation.
#' Threshold sources are recorded for new results; older saved results show
#' `not_recorded`. No confidence intervals or selection-adjusted uncertainty
#' are calculated. Retain the simulation's generating assumptions and replication
#' evidence with the recommendation.
#' @return For `summary()`: `summary()` returns a `summary.mfrm_design_recommendation` list
#' containing `target`, `thresholds`, `checks`, `design_table`, `facet_table`,
#' `recommended`, and `notes`. Print methods return their input invisibly.
#' @seealso [recommend_mfrm_design()], [plot.mfrm_design_recommendation()],
#'   [mfrm_report()], [export_summary_appendix()]
#' @examples # Saved synthetic simulation; one replication illustrates structure only.
#' saved <- readRDS(system.file("examples", "simulation-design.rds",
#'   package = "mfrmr", mustWork = TRUE))
#' rec <- recommend_mfrm_design(saved$result, facets = "Person")
#' summary(rec)
#' report <- mfrm_report(rec)
#' report$tables$checks
#' @export
summary.mfrm_design_recommendation <- function(object, ...) {
  rlang::check_dots_empty()
  x <- object
  if (!is.data.frame(x$design_table) || !is.data.frame(x$facet_table) ||
      !is.list(x$thresholds) || !length(x$thresholds$facets))
    stop("The recommendation must retain its design/facet tables and thresholds.", call. = FALSE)
  designs <- as.data.frame(x$design_table)
  designs$DesignRow <- seq_len(nrow(designs))
  designs$Selected <- FALSE
  selected <- which(designs$Pass %in% TRUE)
  if (length(selected)) designs$Selected[selected[1L]] <- TRUE
  keys <- simulation_design_group_variables(x$design_descriptor)
  facets <- dplyr::left_join(x$facet_table, designs[c(keys, "DesignRow")], by = keys)
  criteria <- c(Separation = "MeanSeparation", Reliability = "MeanReliability",
    SeverityRMSE = "MeanSeverityRMSE", MisfitRate = "MeanMisfitRate",
    Convergence = "ConvergenceRate")
  limits <- c("min_separation", "min_reliability", "max_severity_rmse",
    "max_misfit_rate", "min_convergence_rate")
  checks <- lapply(seq_along(criteria), function(i) {
    value <- facets[[criteria[i]]]
    data.frame(DesignRow = facets$DesignRow, Facet = facets$Facet,
      Criterion = names(criteria)[i], Value = value,
      Relation = if (startsWith(limits[i], "min")) ">=" else "<=",
      Threshold = x$thresholds[[limits[i]]],
      Status = ifelse(!is.finite(value), "unavailable",
        ifelse(facets[[paste0(names(criteria)[i], "Pass")]] %in% TRUE, "pass", "fail")))
  })
  add <- function(name, value, limit, relation, pass, active = TRUE) {
    data.frame(DesignRow = designs$DesignRow, Facet = "All requested facets",
      Criterion = name, Value = value, Relation = relation,
      Threshold = if (is.null(limit)) NA_real_ else limit,
      Status = if (!active) rep("not_requested", nrow(designs)) else
        ifelse(!is.finite(value), "unavailable", ifelse(pass %in% TRUE, "pass", "fail")))
  }
  checks <- dplyr::bind_rows(checks,
    add("Facet coverage", vapply(designs$DesignRow, function(i)
      sum(x$thresholds$facets %in% facets$Facet[facets$DesignRow == i]), integer(1)),
      designs$FacetsRequired, "=", designs$FacetsMissing == ""),
    add("Total ratings", designs$MaxRatings, x$thresholds$max_ratings, "<=",
      designs$RatingsPass, !is.null(x$thresholds$max_ratings)),
    add("Ratings per rater", designs$MaxRatingsPerRater, x$thresholds$max_ratings_per_rater, "<=",
      designs$RaterWorkloadPass, !is.null(x$thresholds$max_ratings_per_rater)),
    add("Connected assignments", ifelse(designs$ConnectivityStatus == "not_assessed", NA_real_,
      as.numeric(designs$ConnectivityStatus == "connected")), 1, "=", designs$ConnectivityPass,
      isTRUE(x$thresholds$require_connected)))
  thresholds <- data.frame(Setting = names(x$thresholds),
    Value = vapply(x$thresholds, function(v) if (is.null(v)) "not_requested" else paste(v, collapse = ", "), character(1)),
    Source = unname(x$threshold_sources[names(x$thresholds)] %||%
      rep("not_recorded", length(x$thresholds))), stringsAsFactors = FALSE)
  thresholds$Source[is.na(thresholds$Source)] <- "not_recorded"
  notes <- c(unlist(x$caveats, use.names = FALSE),
    "Separation reliability describes reproducible differences within a facet, not inter-rater agreement. Low separation can reflect small differences or insufficient precision.",
    "Thresholds are planning choices, not universal standards. Passing describes the evaluated conditions only; no minimum sample size outside those conditions is established.",
    "Averaged simulation metrics do not guarantee every replication meets a threshold. Inspect the source evaluation's replication counts, failures and Monte Carlo uncertainty.",
    "No plan-selection uncertainty or confidence intervals are supplied. Sparse-overlap review is separate from Pass.")
  structure(list(target = .mfrm_design_target(x$thresholds$facets),
    thresholds = thresholds, checks = as.data.frame(checks), design_table = designs,
    facet_table = as.data.frame(x$facet_table), recommended = x$recommended, notes = notes),
    class = "summary.mfrm_design_recommendation")
}

#' @rdname recommend_mfrm_design
#' @export
print.summary.mfrm_design_recommendation <- function(x, ...) {
  rlang::check_dots_empty()
  cat("Simulation-based design recommendation\n")
  print(x$target, row.names = FALSE)
  cat("\n", sum(x$design_table$Pass %in% TRUE), "of", nrow(x$design_table), "evaluated designs pass.\n")
  if (nrow(x$recommended)) print(x$recommended) else cat("No evaluated design meets all requested criteria.\n")
  cat("Use $checks for all criterion outcomes and $thresholds for their source.\n")
  cat("This screen does not establish agreement, classification accuracy or a universal sample size.\n")
  invisible(x)
}

#' @rdname recommend_mfrm_design
#' @name print.mfrm_design_recommendation
#' @aliases print.mfrm_design_recommendation
#' @usage \method{print}{mfrm_design_recommendation}(x, ...)
#' @rawNamespace S3method(print,mfrm_design_recommendation,print.mfrm_bundle)
NULL

#' @rdname recommend_mfrm_design
#' @details `plot()` method: Plot criteria behind a design recommendation
#' @param draw Logical; draw the criterion-status matrix. Use `FALSE` to
#'   retrieve every criterion and source table without drawing.
#' @param ... Unused; extra arguments are rejected.
#' @return For `plot()`: An invisible `mfrm_plot_data` with `table` (long criterion outcomes),
#'   `designs`, `targets`, `thresholds` and interpretation notes.
#' @details For `plot()`: Columns are ranked candidates, labelled by `DesignRow`; an asterisk
#'   marks the selected candidate. Rows combine the facet and criterion.
#'   Grey means unavailable and white means not requested, neither a passing
#'   result. Selection follows the stored ordering, not a new utility function.
#'   For custom ggplot/raincloud graphics, use `plot_data(x, "table")` or obtain
#'   full design metrics with `plot_data(x, "designs", format = "long",
#'   value_cols = c("MinReliability", "MaxRatings"))`. These aggregate metrics
#'   are not replication-level distributions for a raincloud plot.
#' @seealso [summary.mfrm_design_recommendation()], [as_ggplot()], [plot_data()]
#' @examples saved <- readRDS(system.file("examples", "simulation-design.rds",
#'   package = "mfrmr", mustWork = TRUE))
#' rec <- recommend_mfrm_design(saved$result, facets = "Person")
#' plot(rec)
#' head(plot_data(rec, "table"))
#' @export
plot.mfrm_design_recommendation <- function(x, draw = TRUE, ...) {
  rlang::check_dots_empty()
  s <- summary(x)
  tab <- s$checks
  tab$Check <- paste(ifelse(tab$Facet == "All requested facets", "Design", tab$Facet),
    tab$Criterion, sep = ": ")
  out <- new_mfrm_plot_data("design_recommendation", list(table = tab,
    designs = s$design_table, targets = s$target, thresholds = s$thresholds,
    notes = data.frame(Type = "interpretation", Text = s$notes)))
  if (!isTRUE(draw)) return(out)
  levels <- c("pass", "fail", "unavailable", "not_requested")
  colors <- c("#2A9D8F", "#D55E00", "#999999", "#FFFFFF")
  criteria <- rev(unique(tab$Check))
  n <- nrow(s$design_table)
  if (!n) stop("No evaluated designs to plot.", call. = FALSE)
  left <- max(graphics::strwidth(criteria, units = "inches", cex = .7)) + .3
  old <- graphics::par(mai = c(.8, left, .15, .15))
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot(NA_real_, NA_real_, xlim = c(.5, n + .5),
    ylim = c(.5, length(criteria) + 2.5), xlab = "Ranked design row (* selected)", ylab = "", axes = FALSE)
  graphics::axis(1, seq_len(n), labels = paste0(seq_len(n), ifelse(s$design_table$Selected, " *", "")))
  graphics::axis(2, seq_along(criteria), labels = criteria, las = 1, cex.axis = .7)
  graphics::points(tab$DesignRow, match(tab$Check, criteria), pch = 22, cex = 1.5,
    bg = colors[match(tab$Status, levels)], col = "#444444")
  graphics::legend("topright", legend = levels, pch = 22, pt.bg = colors, bty = "n", cex = .7, ncol = 2)
  invisible(out)
}
