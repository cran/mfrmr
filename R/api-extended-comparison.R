# Descriptive, matched-event comparisons. Estimation is never performed here.
mfrm_compare_events <- function(data, columns) {
  if (!is.data.frame(data) || !all(columns %in% names(data))) stop(
    "Required rating-event columns are absent; supply complete saved fits and assigned-score data.", call. = FALSE)
  out <- as.data.frame(data[columns], stringsAsFactors = FALSE)
  out[] <- lapply(out, as.character)
  if (nrow(out)) out <- out[do.call(order, c(out, list(na.last = TRUE, method = "radix"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

mfrm_compare_extended <- function(fits, labels, nested, response_diagnostics = NULL, person_scores = NULL) {
  extended <- vapply(fits, mfrm_extended_fit, logical(1))
  if (length(fits) != 2L || sum(extended) != 1L || !inherits(fits[[which(!extended)[1]]], "mfrm_fit")) stop(
    "Compare exactly one ordinary mfrm_fit and one testlet or shared-rater fit.", call. = FALSE)
  if (!is.logical(nested) || length(nested) != 1L || is.na(nested) || nested) stop(
    "Extended-model comparisons are descriptive; nested = TRUE and ordinary chi-squared LRTs are not supported.", call. = FALSE)
  ei <- which(extended); oi <- which(!extended); ext <- fits[[ei]]; ordinary <- fits[[oi]]
  cfg <- ordinary$config; prep <- ordinary$prep; cols <- ext$input$columns
  testlet <- inherits(ext, "mfrm_testlet")
  facets <- if (testlet) cols$facets else c(cols$rater, cols$facets)
  source <- prep$source_columns %||% cfg$source_columns
  if (!identical(toupper(cfg$model), "RSM") || !identical(toupper(cfg$method), "MML") ||
      !setequal(source$facets, facets) || !identical(source$person, cols$person) ||
      !identical(source$score, cols$score)) stop(
    "Use an ordinary RSM MML fit with the same person/score columns and fixed facets (including the shared rater as a fixed facet).", call. = FALSE)
  anchors <- extract_anchor_tables(cfg)
  if (length(cfg$facet_interactions) || nrow(as.data.frame(anchors$anchors)) ||
      nrow(as.data.frame(anchors$groups)) || length(cfg$dummy_facets) ||
      !identical(cfg$noncenter_facet, "Person") ||
      !identical(cfg$facet_shrinkage %||% "none", "none") ||
      length(cfg$facet_signs[facets]) != length(facets) || anyNA(cfg$facet_signs[facets]) ||
      !all(cfg$facet_signs[facets] == -1) ||
      !is.data.frame(prep$data) || !all(is.finite(prep$data$Weight) & prep$data$Weight == 1)) stop(
    "Use unit weights, additive unanchored severity facets, Person as the noncentered facet, and no facet shrinkage or interactions.", call. = FALSE)
  categories <- ext$input$score_levels
  map <- prep$score_map
  if (!is.data.frame(map) || !all(c("OriginalScore", "InternalScore") %in% names(map))) stop(
    "The ordinary fit lacks a verified category mapping.", call. = FALSE)
  map <- map[order(map$InternalScore), , drop = FALSE]
  if (!identical(as.numeric(map$OriginalScore), as.numeric(categories)) ||
      any(diff(map$InternalScore) != 1) ||
      !all(prep$data$score_k == match(prep$data$Score, map$InternalScore) - 1L)) stop(
    "Both fits must retain the same ordered consecutive categories and adjacent-category coding.", call. = FALSE)
  pop <- ordinary$population %||% cfg$population_spec
  estimated <- isTRUE(pop$active)
  ext_estimated <- !is.null(ext$calibration$person_variance) && is.null(ext$settings$fixed_person_sd)
  if (estimated && (!identical(pop$design_columns, "(Intercept)") ||
      !is.matrix(cfg$population_spec$design_matrix) ||
      ncol(cfg$population_spec$design_matrix) != 1L ||
      any(!is.finite(cfg$population_spec$design_matrix)) ||
      any(cfg$population_spec$design_matrix != 1))) stop(
    "The ordinary reference must use an intercept-only normal population (population_formula = ~1).", call. = FALSE)
  if (estimated != ext_estimated || (!estimated && !isTRUE(all.equal(ext$calibration$person_sd %||% 1, 1)))) stop(
    "Match the population assumption: population_formula = ~1 for estimated ability SD, or the ordinary default and person_sd = 1 for known N(0,1).", call. = FALSE)
  if (isTRUE((pop$response_rows_omitted %||% 0) > 0)) stop(
    "Population-data omissions are outside this matched assigned-score comparison.", call. = FALSE)
  columns <- c(cols$person, sort(facets), cols$score)
  ord_data <- prep$data[c("Person", sort(facets), "Score")]
  names(ord_data) <- columns
  ord_data[[cols$score]] <- map$OriginalScore[match(prep$data$Score, map$InternalScore)]
  events <- mfrm_compare_events(ord_data, columns)
  if (!identical(events, mfrm_compare_events(ext$input$data, columns))) stop(
    "Observed rating events differ, including IDs, scores or repeated-event multiplicities; fit both models to the same retained rows.", call. = FALSE)
  retention <- prep$row_retention
  input_n <- retention$Rows[match("input_selected_columns", retention$Stage)]
  if (length(input_n) != 1L || !is.finite(input_n) || input_n < nrow(events) ||
      !identical(as.numeric(input_n), as.numeric(ext$input$input_rows)) ||
      input_n - nrow(events) != length(ext$input$omitted_rows)) stop(
    "Input and omitted-row accounting differ or are unavailable; use the same assigned-score roster.", call. = FALSE)
  omitted_count <- input_n - nrow(events)
  all_columns <- unique(c(columns, if (testlet) cols$testlet))
  omission_rows <- ext$input$data[FALSE, all_columns, drop = FALSE]
  if (omitted_count > 0) {
    omitted <- prep$omitted_data
    if (!is.data.frame(omitted) || nrow(omitted) != omitted_count || is.null(ext$input$assigned_data)) stop(
      "Omitted-event identities are unavailable in an older fit. Refit from the complete assigned-score roster to retain omission provenance.", call. = FALSE)
    if (anyNA(omitted[c("Person", facets)]) || any(!is.na(omitted$Score)) ||
        anyNA(omitted$Weight) || any(omitted$Weight != 1)) stop(
      "This comparison supports missing assigned scores only, not missing IDs or excluded weights.", call. = FALSE)
    omitted <- omitted[c("Person", sort(facets), "Score")]; names(omitted) <- columns
    omission_rows <- ext$input$assigned_data[is.na(ext$input$assigned_data[[cols$score]]), all_columns, drop = FALSE]
    if (!identical(mfrm_compare_events(omitted, columns), mfrm_compare_events(omission_rows, columns))) stop(
      "Omitted rating-event identities differ, despite the observed events or counts matching.", call. = FALSE)
  }
  defaults <- c("Ordinary fixed-facet RSM", if (testlet) "Testlet RSM" else "Shared-rater RSM")
  if (ei == 1L) defaults <- rev(defaults)
  labels <- labels %||% defaults
  if (!is.character(labels) || length(labels) != 2L || anyNA(labels) ||
      any(!nzchar(trimws(labels))) || anyDuplicated(labels)) stop("Use two distinct nonempty labels.", call. = FALSE)
  ord_ready <- identical(as.character(ordinary$summary$NumericalState), "ready") &&
    isTRUE(mfrm_convergence_state(ordinary)$code_converged)
  ext_ready <- isTRUE(ext$checks$NumericalReady) && isTRUE(ext$checks$InformationPositive)
  ready <- if (ei == 1L) c(ext_ready, ord_ready) else c(ord_ready, ext_ready)
  usable <- ready
  usable[oi] <- usable[oi] && !identical(ordinary$readiness$fit$FitReadiness, "blocked")
  ord_effects <- as.data.frame(ordinary$facets$others[intersect(c("Facet", "Level", "Estimate", "ParameterStatus"), names(ordinary$facets$others))])
  ext_effects <- ext$calibration_table[ext$calibration_table$Parameter == "Fixed facet", c("Facet", "Level", "Estimate"), drop = FALSE]
  if (!testlet) ext_effects <- rbind(ext_effects, data.frame(Facet = cols$rater,
    Level = ext$raters$Rater, Estimate = ext$raters$Estimate))
  tables <- if (ei == 1L) list(ext_effects, ord_effects) else list(ord_effects, ext_effects)
  effects <- do.call(rbind, lapply(sort(facets), function(facet) {
    levels <- sort(unique(as.character(ord_data[[facet]])))
    parts <- lapply(seq_along(tables), function(i) {
      table <- tables[[i]][tables[[i]]$Facet == facet, , drop = FALSE]
      if (anyDuplicated(table$Level)) stop("Facet estimates have duplicate level identities.", call. = FALSE)
      raw <- table$Estimate[match(levels, as.character(table$Level))]
      excluded_status <- "ParameterStatus" %in% names(table) &&
        any(!table$ParameterStatus[match(levels, as.character(table$Level))] %in% c("estimable", "fixed", "weak_information"))
      reason <- if (!usable[i]) "Source checks require review" else if (excluded_status)
        "Source parameter status excludes comparison" else if (any(!is.finite(raw)))
        "All matched levels need finite estimates for centering" else ""
      center <- if (!nzchar(reason)) mean(raw) else NA_real_
      list(raw = raw, center = center, value = raw - center, reason = reason)
    })
    data.frame(Facet = facet, Level = levels,
      SourceReference = parts[[1]]$raw, SourceComparison = parts[[2]]$raw,
      CenterReference = parts[[1]]$center, CenterComparison = parts[[2]]$center,
      Reference = parts[[1]]$value, Comparison = parts[[2]]$value,
      Mean = (parts[[1]]$value + parts[[2]]$value)/2,
      Difference = parts[[2]]$value - parts[[1]]$value,
      ReferenceKind = if (ei == 1L && !testlet && facet == cols$rater) "Conditional rater mode" else "Fixed facet estimate",
      ComparisonKind = if (ei == 2L && !testlet && facet == cols$rater) "Conditional rater mode" else "Fixed facet estimate",
      Status = if (all(vapply(parts, function(z) !nzchar(z$reason), logical(1)))) "available_descriptive" else "unavailable",
      Reason = paste(Filter(nzchar, unique(vapply(parts, `[[`, character(1), "reason"))), collapse = "; "),
      stringsAsFactors = FALSE)
  }))
  models <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]; is_ext <- extended[i]
    data.frame(Label = labels[i], Model = defaults[i],
      Input = input_n, Observed = nrow(events), Omitted = input_n - nrow(events), Persons = length(unique(events[[cols$person]])),
      AbilityMeanSource = if (is_ext) 0 else if (estimated) as.numeric(pop$coefficients[1]) else 0,
      AbilitySD = if (is_ext) f$calibration$person_sd %||% 1 else if (estimated) sqrt(pop$sigma2) else 1,
      AbilitySDModel = if (estimated) "Estimated common normal SD" else "Known N(0,1)",
      LogLik = if (is_ext) f$loglik else as.numeric(f$summary$LogLik),
      NumericalReady = ready[i], DescriptiveSourceAvailable = usable[i],
      InferenceReadiness = if (is_ext) "Target-specific; general qualification incomplete" else as.character(f$readiness$fit$FitReadiness %||% "unavailable"),
      Location = if (is_ext) "Mean-zero population; free step location" else if (estimated) "Estimated population intercept; centered steps" else "Known mean-zero population",
      Integration = if (is_ext) f$settings$method else paste("MML", f$summary$MMLIntegration),
      stringsAsFactors = FALSE)
  })
  checks <- data.frame(Check = c("Observed rating events", "Categories", "Assigned-score omissions", "Facet specification", "Ability population", "Unit and orientation", "Numerical checks"),
    Status = c(rep("matched", 6), if (all(ready)) "passed" else "review"),
    Detail = c("Identical multiset; repeated events retain their multiplicity.",
      "Identical ordered consecutive categories and adjacent logits.",
      if (input_n == nrow(events)) "Neither fit omitted assigned rows." else "Saved omitted-event identities agree; only assigned scores are missing.",
      if (testlet) "Same fixed facets; Person-specific local dependence is the changed structure." else "Same other fixed facets; observed raters change from fixed effects to shared random effects.",
      if (estimated) "Both estimate one common normal SD; raw population locations use different identification conventions." else "Both impose known N(0,1).",
      "Unit Rasch slope; positive facet effects mean greater severity/difficulty.",
      "Descriptive readiness only; does not establish adequacy or interval coverage."))
  checks <- rbind(checks, data.frame(Check = "Descriptive source checks",
    Status = if (all(usable)) "passed" else "review",
    Detail = "Blocked source readiness withholds differences; parameter exclusions remain in the effects table."))
  notes <- c("Effects are centered at the unweighted mean over the same complete set of levels within each facet. Pairwise level contrasts are preserved; this is not a scale calibration or a statistical test.",
    "Differences are comparison minus reference. Fixed rater coefficients and conditional random-rater modes have different estimation targets and shrinkage; a smaller spread does not establish improvement.",
    "No SE/CI for model differences, AIC weights, Person-count BIC, chi-squared LRT or preferred model is supplied. Log likelihoods retain their own integration basis and are not ranked.",
    "Person scores, step locations, predictive fit statistics and replacement-rater predictions are not compared by this facet-effect output.")
  responses <- if (!is.null(response_diagnostics)) mfrm_compare_responses(fits, response_diagnostics, columns) else NULL
  persons <- if (!is.null(person_scores)) mfrm_compare_persons(fits,person_scores) else NULL
  if (!is.null(persons)) notes <- c(notes,persons$note)
  if (!is.null(responses)) {
    checks <- rbind(checks, data.frame(Check = "Selected posterior predictive events", Status = "matched",
      Detail = "Identical selected-event multiset, probability target and group accounting; all observed source events condition each model."))
    notes <- c(notes, responses$note)
  }
  structure(list(models = do.call(rbind, models), checks = checks, effects = effects, responses = responses, persons = persons,
    notes = data.frame(Note = notes), source = mfrm_extended_prediction_source(ext),
    source_events = mfrm_compare_events(ext$input$data, unique(c(columns, if (testlet) cols$testlet))),
    omitted_events = mfrm_compare_events(omission_rows, all_columns), event_columns = all_columns, source_settings = list(ordinary = cfg$replay_inputs %||% cfg, extension = ext$settings),
    numerical_checks = list(ordinary = ordinary$readiness, extension = ext$checks),
    labels = labels, extension_index = ei), class = "mfrm_extended_comparison")
}

#' @rdname compare_mfrm
#' @export
summary.mfrm_extended_comparison <- function(object, ...) {
  list(models = object$models, checks = object$checks, effects = object$effects,
    responses = object$responses, persons = object$persons, notes = object$notes)
}

#' @rdname compare_mfrm
#' @export
print.mfrm_extended_comparison <- function(x, ...) {
  cat("Matched-event facet comparison (descriptive)\n")
  print(x$models, row.names = FALSE)
  print(x$effects[c("Facet", "Level", "Reference", "Comparison", "Difference", "Status")], row.names = FALSE)
  cat("Centered within each facet; no difference intervals or automatic model preference.\n")
  if (!is.null(x$responses)) {
    cat("Saved posterior predictive group comparisons (same data, descriptive only):\n")
    print(x$responses$measures, row.names = FALSE)
  }
  if (!is.null(x$persons)) {
    cat("Population-centered conditional Person scores (descriptive only):\n")
    print(x$persons$table, row.names=FALSE)
  }
  invisible(x)
}

#' @rdname compare_mfrm
#' @details `plot()` method: Plot matched facet effects or posterior predictions across MFRMs
#' @param x A descriptive extended-model result from [compare_mfrm()].
#'
#' Also for `compare_mfrm()`: An extended-model comparison returned by `compare_mfrm()`.
#' @param style `"paired"` plots comparison against reference with an equality
#'   line. `"difference"` plots comparison minus reference against their mean.
#' @param facet Panels to display; `NULL` includes all. For effects and group
#'   indices these are facet names; for row moments the panel is
#'   `"Selected ratings"`; for probabilities it is `"Score k"` for category k.
#' @param metric `"effects"` compares centered facet locations (default).
#'   `"expected_score"`, `"variance"`, `"infit"`, `"outfit"` and
#'   `"probability"` require saved `response_diagnostics` in [compare_mfrm()].
#'   These compare the same posterior predictive target in original score
#'   units, squared score units, descriptive mean squares or probabilities.
#'   `"person"` requires saved `person_scores` and compares EAPs centered at
#'   each fitted population mean, in logits. No difference interval is drawn.
#' @param category Numeric score categories to display with
#'   `metric = "probability"`; `NULL` includes all categories.
#' @inheritParams plot.mfrm_testlet_scores
#' @param palette `"accessible"` uses blue points and a dashed grey reference;
#'   `"mono"` uses dark grey points. Meaning does not depend on colour.
#' @param show_labels Show matched IDs beside the points. Defaults to `TRUE`
#'   for effects and group indices, `FALSE` for rating-level moments and
#'   probabilities. Consult the retained table when labels are hidden.
#' @details For `plot()`: Points are centered at the unweighted mean within each complete
#'   matched facet. These displays compare fitted summaries, including
#'   shrinkage, rather than testing differences or choosing a better model.
#'   The equality/zero line is a descriptive reference. No difference intervals
#'   or limits of agreement are computed. Unavailable rows remain in plot data.
#'   Select `facet` or use `show_labels = FALSE` for crowded displays.
#'
#'   Predictive quantities are not centered. Each row label maps through
#'   `data$response_events` and the comparison's `responses$rows`, which retain
#'   both original input row numbers. Infit/Outfit use selected rows in each
#'   group; posterior conditioning retains all observed source ratings.
#'   Smaller residual indices do not establish better prediction or model
#'   adequacy. There are no ordinary fit cutoffs, tests or automatic ranking.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object with source effects, checks,
#'   interpretation notes and a text alternative. [as_ggplot()] supports further
#'   styling without estimation. Both renderers respect the display controls.
#' @export
plot.mfrm_extended_comparison <- function(x, style = c("paired", "difference"),
    facet = NULL, draw = TRUE, palette = c("accessible", "mono"), title = NULL,
    caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = metric %in% c("effects", "infit", "outfit"),
    text_scale = 1, point_size = 2.5,
    metric = c("effects", "expected_score", "variance", "infit", "outfit", "probability", "person"), category = NULL, ...) {
  rlang::check_dots_empty(); style <- match.arg(style); palette <- match.arg(palette)
  metric <- match.arg(metric)
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
  response <- if (metric == "person") {
    if (is.null(x$persons)) stop("Supply matching saved person_scores to compare_mfrm() first.",call.=FALSE)
    if (!is.null(category)) stop("category applies only to probabilities.",call.=FALSE)
    d <- x$persons$table; d$Facet <- "Person"; d$Level <- d$Person
    d$Reference[d$Status!="available_descriptive"] <- d$Comparison[d$Status!="available_descriptive"] <- NA_real_
    d$Mean <- (d$Reference+d$Comparison)/2
    list(table=d,unit="population-centered ability (logits)")
  } else if (metric != "effects") mfrm_response_comparison_plot_table(x, metric, category) else NULL
  if (metric == "effects" && !is.null(category)) stop("`category` applies only to probability comparisons.", call. = FALSE)
  source_table <- if (is.null(response)) x$effects else response$table
  facets <- unique(source_table$Facet); facet <- facet %||% facets
  if (!is.character(facet) || !length(facet) || anyNA(facet) || anyDuplicated(facet) || !all(facet %in% facets)) stop("Choose available facet names or panel names for the selected metric.", call. = FALSE)
  tab <- source_table[source_table$Facet %in% facet, , drop = FALSE]
  tab$Facet <- factor(tab$Facet, levels = facet)
  tab$X <- if (style == "paired") tab$Reference else tab$Mean
  tab$Y <- if (style == "paired") tab$Comparison else tab$Difference
  note <- paste(if (metric == "effects") "Centered within each facet; descriptive, no difference intervals." else if(metric=="person")
    "Population-centered conditional EAPs; calibration fixed; no difference intervals or ranking." else
    "Same-data posterior predictions; calibration fixed; descriptive, no model ranking or fit cutoffs.",
    sprintf("%d of %d matched summaries displayed.", sum(is.finite(tab$X) & is.finite(tab$Y)), nrow(tab)))
  unit <- if (metric == "effects") "centered logits" else response$unit
  expand_range <- function(values, zero = FALSE) {
    values <- values[is.finite(values)]
    limits <- range(c(if (zero) 0, if (!length(values)) c(-.1, .1), values))
    limits + c(-1, 1) * max(.01, diff(limits) * .18)
  }
  xlim <- if (style == "paired") {
    if (metric == "probability") c(0, 1) else expand_range(c(tab$X, tab$Y))
  } else expand_range(tab$X)
  ylim <- if (style == "paired") xlim else expand_range(tab$Y, zero = TRUE)
  payload <- list(effects = x$effects, table = tab, models = x$models, checks = x$checks,
    notes = data.frame(Type = "Interpretation", Text = x$notes$Note),
    response_events = x$responses$events,
    xlim = xlim, ylim = ylim,
    title = title %||% if (metric == "effects") "Facet effects across models" else paste("Across models:", unit), caption = caption %||% note,
    xlab = if (style == "paired") paste0(x$labels[1], " (", unit, ")") else paste0("Mean across models (", unit, ")"),
    ylab = if (style == "paired") paste0(x$labels[2], " (", unit, ")") else paste0("Comparison - reference (", unit, ")"),
    alt_text = paste("Matched", if (metric == "effects") "facet effects:" else if(metric=="person") "conditional Person summaries:" else "posterior predictive summaries:", paste(facet, collapse = ", "), paste0("(", style, " view)."), note,
      if (metric == "effects") "Fixed facet estimates and conditional random-rater modes are different fitted summaries." else if(metric=="person") x$persons$note else x$responses$note),
    display = list(style = style, palette = palette, facet = facet, metric = metric,
      show_title = show_title, show_notes = show_notes, show_labels = show_labels,
      text_scale = text_scale, point_size = point_size))
  out <- new_mfrm_plot_data("extended_model_comparison", payload)
  if (draw) {
    cols <- min(3, ceiling(sqrt(length(facet)))); rows <- ceiling(length(facet)/cols)
    old <- graphics::par(mfrow = c(rows, cols), mar = c(5, 5, 3, 1),
      oma = c(if (show_notes) 4 else 0, 0, if (show_title) 3 else 0, 0), cex = text_scale)
    on.exit(graphics::par(old), add = TRUE)
    colour <- if (palette == "mono") "#222222" else "#0072B2"
    for (f in facet) {
      d <- tab[tab$Facet == f & is.finite(tab$X) & is.finite(tab$Y), , drop = FALSE]
      graphics::plot(NA_real_, NA_real_, xlim = xlim, ylim = ylim,
        asp = if (style == "paired") 1 else NA_real_, xlab = paste(strwrap(payload$xlab, 32), collapse = "\n"),
        ylab = paste(strwrap(payload$ylab, 32), collapse = "\n"), main = f, bty = "n")
      if (style == "paired") graphics::abline(a = 0, b = 1, col = "grey45", lty = 2) else graphics::abline(h = 0, col = "grey45", lty = 2)
      if (nrow(d)) {
        graphics::points(d$X, d$Y, pch = 16, col = colour, cex = point_size/2.5)
        if (show_labels) graphics::text(d$X, d$Y, d$Level, pos = 3, cex = .75)
      } else graphics::text(mean(xlim), mean(ylim), "No available summaries")
    }
    if (show_title) graphics::mtext(payload$title, side = 3, outer = TRUE, line = 1)
    if (show_notes) graphics::mtext(paste(strwrap(payload$caption, 85), collapse = "\n"), side = 1, outer = TRUE, line = 1, cex = .7)
  }
  invisible(out)
}

.mfrmr_gg_extended_comparison <- function(payload) {
  opt <- payload$display; tab <- payload$table
  points <- tab[is.finite(tab$X) & is.finite(tab$Y), , drop = FALSE]
  colour <- if (opt$palette == "mono") "#222222" else "#0072B2"
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$X, .data$Y))
  if (opt$style == "paired") {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
      ggplot2::coord_equal(xlim = payload$xlim, ylim = payload$ylim, expand = FALSE)
  } else p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .18)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = .18))
  p <- p + ggplot2::geom_point(data = points, colour = colour, size = opt$point_size)
  if (opt$show_labels) p <- p + ggplot2::geom_text(data = points, ggplot2::aes(label = .data$Level),
    vjust = -1, size = 3 * opt$text_scale)
  absent <- setdiff(levels(tab$Facet), as.character(points$Facet))
  if (length(absent)) p <- p + ggplot2::geom_text(data = data.frame(Facet = factor(absent, levels = levels(tab$Facet)),
    X = mean(payload$xlim %||% c(-.1, .1)), Y = mean(payload$ylim %||% c(-.1, .1))),
    label = "No available summaries")
  p <- p + ggplot2::facet_wrap(~ Facet, ncol = min(3, ceiling(sqrt(nlevels(tab$Facet))))) + .mfrmr_gg_theme() +
    ggplot2::theme(text = ggplot2::element_text(size = 11 * opt$text_scale),
      panel.spacing = ggplot2::unit(1.5, "lines"))
  p <- .mfrmr_gg_labs(p, payload, x = payload$xlab, y = payload$ylab) + ggplot2::labs(alt = payload$alt_text)
  attr(p, "mfrmr_alt_text") <- payload$alt_text
  p
}
