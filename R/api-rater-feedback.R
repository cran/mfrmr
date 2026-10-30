# A distributed sheet is a whitelist projection, never a wrapper around a fit.
mfrm_rater_feedback <- function(x, facet, rater, audience, label, max_cases, interval) {
  scalar_text <- function(value, name) {
    if (length(value) != 1L || is.na(value) || !is.character(value) || !nzchar(value)) {
      stop(paste0("`", name, "` must be one nonempty string."), call. = FALSE)
    }
    unname(as.character(value))
  }
  facet <- scalar_text(facet, "facet")
  rater <- scalar_text(rater, "rater")
  label <- if (is.null(label)) "Selected rater" else scalar_text(label, "label")
  audience <- match.arg(audience, c("rater", "researcher"))
  if (length(max_cases) != 1L || !is.numeric(max_cases) || is.na(max_cases) ||
      !is.finite(max_cases) || max_cases < 0 || max_cases != floor(max_cases)) {
    stop("`max_cases` must be a nonnegative whole number; use 0 to omit cases.", call. = FALSE)
  }
  fit <- x[["fit"]]
  if (mfrm_has_product_slopes(fit)) {
    if (!is.null(interval)) interval <- tolower(scalar_text(interval, "interval"))
    return(mfrm_gmfrm_rater_feedback(x, facet, rater, audience, label, max_cases, interval))
  }
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !isTRUE(fit$config$model %in% c("RSM", "PCM")) ||
      !isTRUE(fit$config$method %in% c("MML", "JML")) ||
      length(fit$config$interaction_specs) || length(fit$interactions$specs)) {
    stop("Rater sheets require a native additive RSM/PCM fit. GPCM, interaction, imported, testlet and random-rater models need their model-specific reports.", call. = FALSE)
  }
  if (!facet %in% fit$prep$facet_names) stop("`facet` must name a fitted non-Person facet.", call. = FALSE)
  coef <- fit$facets$others
  selected <- which(as.character(coef$Facet) == facet & as.character(coef$Level) == rater)
  if (length(selected) != 1L) stop("`rater` must select exactly one level of `facet`.", call. = FALSE)
  sign <- as.numeric(fit$config$facet_signs[facet])
  if (length(sign) != 1L || is.na(sign) || !sign %in% c(-1, 1)) {
    stop("The selected facet lacks a supported effect direction.", call. = FALSE)
  }
  estimate <- as.numeric(coef$Estimate[selected])
  ready <- isTRUE(.mfrm_fit_plot_readiness(fit)$ready)
  d <- fit$prep$data
  take <- which(as.character(d[[facet]]) == rater)
  if (!length(take)) stop("The selected rater has no retained rating rows.", call. = FALSE)
  weight <- as.numeric(d$Weight[take])
  if (any(!is.finite(weight) | weight < 0)) stop("Invalid stored rating weights.", call. = FALSE)
  tables <- list(
    severity = data.frame(Severity = -sign * estimate, Reference = 0),
    exposure = data.frame(RatingRows = length(take), Persons = length(unique(d$Person[take])),
                          WeightSum = sum(weight)),
    uncertainty = data.frame(), fit = data.frame(), categories = data.frame(), cases = data.frame()
  )
  notes <- c(
    severity = "Positive severity means lower scores are expected, with other modeled effects held constant. Zero is the fitted model's reference, which need not be the average of the other raters. This model scale is different from the original score scale. Severity is not rater quality.",
    exposure = "Counts cover retained rating rows only. They do not measure completion of planned assignments. Weight sum describes estimation weights; the other counts are unweighted.",
    uncertainty = "No saved individual-rater interval was supplied. An interval is not calculated when creating this sheet.",
    fit = "Saved ordinary-model diagnostics were not supplied. Fit statistics are not calculated when creating this sheet.",
    categories = "Category use describes these ratings and the people assigned to this rater. Unequal use alone does not show a problem with the rubric or the rater.",
    cases = "Saved ordinary-model diagnostics were not supplied. Unexpected ratings are not calculated when creating this sheet."
  )
  # Select an individual coefficient by its contrast vector, never its row label.
  intervals <- mfrm_facet_results_inputs(fit, x[["facet_intervals"]])
  interval_row <- function(ci) {
    if (!identical(ci$settings$facet, facet)) return(integer(0))
    j <- ci$contrasts
    levels <- as.character(coef$Level[as.character(coef$Facet) == facet])
    if (!is.matrix(j) || is.null(colnames(j)) || anyDuplicated(colnames(j)) ||
        !setequal(colnames(j), levels) || anyNA(j) || nrow(j) != nrow(ci$table)) return(integer(0))
    wanted <- as.numeric(colnames(j) == rater)
    which(apply(j, 1L, function(row) all(row == wanted)))
  }
  candidates <- lapply(intervals, interval_row)
  eligible <- names(candidates)[lengths(candidates) > 0L]
  if (!is.null(interval)) {
    interval <- tolower(scalar_text(interval, "interval"))
    if (!interval %in% eligible) stop("`interval` must name a saved interval containing this individual rater coefficient.", call. = FALSE)
  } else if (length(eligible) > 1L) {
    stop("Several saved intervals contain this rater. Choose one with `interval`.", call. = FALSE)
  } else if (length(eligible) == 1L) interval <- eligible
  if (!is.null(interval)) {
    ci <- intervals[[interval]]
    row <- candidates[[interval]]
    if (length(row) != 1L) stop("The saved interval repeats this individual rater coefficient.", call. = FALSE)
    tab <- ci$table[row, , drop = FALSE]
    if (!isTRUE(all.equal(as.numeric(tab$Estimate), estimate, tolerance = 1e-8))) {
      stop("The saved individual-rater estimate does not match the fit.", call. = FALSE)
    }
    if (!ci$settings$method %in% c("model", "sandwich")) stop("Unsupported saved interval method.", call. = FALSE)
    ends <- as.numeric(c(tab$Lower, tab$Upper))
    ends <- if (sign == 1) -rev(ends) else ends
    valid_bounds <- all(is.finite(ends)) && ends[1] <= ends[2]
    available <- identical(as.character(tab$Status), "available") && valid_bounds &&
      (!"CIEligible" %in% names(tab) || isTRUE(tab$CIEligible))
    reason <- mfrm_rater_interval_reason(available,
      if (identical(as.character(tab$Status), "available") && !valid_bounds)
        "invalid_bounds" else as.character(tab$Status))
    if (!available) ends[] <- NA_real_
    tables$uncertainty <- data.frame(Level = as.numeric(ci$settings$level),
      Lower = ends[1], Upper = ends[2], Method = as.character(ci$settings$method),
      Available = available, Reason = reason)
    notes["uncertainty"] <- paste(
      if (available) "Saved approximate pointwise interval for this fixed rater's severity." else
        paste("The saved individual-rater interval is unavailable.", reason),
      "It conditions on the fitted facet levels and does not classify rater quality or guarantee coverage.",
      if (ci$settings$method == "sandwich") "The sandwich calculation assumes independent declared clusters and does not correct biased estimates." else
        "The model-based calculation relies on the fitted model's assumptions."
    )
  }
  score <- as.numeric(d$Score[take])
  map <- fit$prep$score_map
  mapped <- is.data.frame(map) && all(c("InternalScore", "OriginalScore") %in% names(map))
  same_scale <- mapped && isTRUE(all.equal(as.numeric(map$InternalScore), as.numeric(map$OriginalScore)))
  categories <- if (mapped) as.numeric(map$InternalScore) else sort(unique(score))
  counts <- vapply(categories, function(k) sum(score == k), numeric(1))
  sums <- vapply(categories, function(k) sum(weight[score == k]), numeric(1))
  tables$categories <- data.frame(ModelScore = categories, Ratings = counts,
    Percent = 100 * counts / length(take), WeightSum = sums)
  if (mapped) tables$categories$OriginalScore <- as.numeric(map$OriginalScore)
  if (!same_scale) notes["categories"] <- paste(notes["categories"],
    "Model scores use the fitted category coding. Expected scores and residuals below stay on that coding; they are not converted to the original scale.")
  diagnostics <- x[["diagnostics"]]
  if (!is.null(diagnostics)) {
    mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "mfrm_report(style = 'rater')")
    measures <- diagnostics$measures
    at <- which(as.character(measures$Facet) == facet & as.character(measures$Level) == rater)
    if (length(at) == 1L && all(c("Infit", "Outfit") %in% names(measures))) {
      tables$fit <- data.frame(Statistic = c("Infit", "Outfit"),
        Value = as.numeric(c(measures$Infit[at], measures$Outfit[at])))
      notes["fit"] <- "These saved statistics compare ratings with ordinary-model expectations. Values above 1 describe more variation than expected; values below 1 describe less. Infit emphasizes informative ratings, while Outfit is more sensitive to unexpected ratings. Neither statistic is a stand-alone verdict on the rater; no automatic warning cutoff is applied. Missing values remain unavailable."
    } else notes["fit"] <- "The saved diagnostics contain no unique Infit/Outfit row for this rater."
    obs <- diagnostics$obs
    fields <- c("Observed", "Expected", "Residual", "StdResidual")
    if (all(c(facet, fields) %in% names(obs))) {
      rows <- which(as.character(obs[[facet]]) == rater & is.finite(obs$StdResidual))
      rows <- utils::head(rows[order(abs(obs$StdResidual[rows]), decreasing = TRUE)], max_cases)
      tables$cases <- data.frame(Case = seq_along(rows),
        Observed = as.numeric(obs$Observed[rows]), Expected = as.numeric(obs$Expected[rows]),
        Residual = as.numeric(obs$Residual[rows]), StandardizedResidual = as.numeric(obs$StdResidual[rows]))
      notes["cases"] <- if (length(rows)) paste(
        "Selected by largest absolute standardized residual, not by a warning threshold.",
        "Case numbers are local to this sheet and are not source row numbers.",
        "Positive residuals mean higher scores than expected. Review context before drawing conclusions."
      ) else "No finite saved standardized residuals are available for the requested case display."
    } else notes["cases"] <- "The saved diagnostics lack the observation fields needed for this section."
  }
  if (max_cases == 0) {
    tables$cases <- data.frame()
    notes["cases"] <- "Individual rating cases were omitted by request."
  }
  technical <- if (audience == "researcher") paste(
    fit$config$model, "model;", fit$config$method,
    "estimation. Severity is minus the facet sign times its fitted coefficient (logits).",
    "Fit and case summaries use saved plug-in expectations, not marginal or posterior predictive extended-model diagnostics.",
    "Reported intervals, when supplied, are fixed-facet pointwise approximations; fit statistics do not estimate diagnostic accuracy."
  ) else "Use this sheet to discuss scoring patterns with the assessment team. Differences may reflect the ratings assigned as well as the model's assumptions."
  if (audience == "rater") {
    notes["severity"] <- paste(
      "Positive values indicate a tendency toward lower scores; negative values indicate higher scores, after accounting for the other modeled effects.",
      "Zero is the model's reference, not necessarily the other raters' average. This model scale differs from the rubric's scores. Severity is not rater quality.")
    if (nrow(tables$uncertainty) && isTRUE(tables$uncertainty$Available)) {
      notes["uncertainty"] <- paste(
        "This approximate interval describes uncertainty in this rater's scoring tendency under the analysis assumptions.",
        "It is not a guarantee of accuracy or a judgement of scoring quality.",
        if (tables$uncertainty$Method == "sandwich") "The sandwich method assumes independent groups and cannot correct biased estimates." else
          "The model-based calculation depends on the model's assumptions.")
    }
    if (nrow(tables$fit)) notes["fit"] <- paste(
      "Values above 1 describe more variation than the model expects; values below 1 describe less.",
      "Infit gives more weight to ratings that help distinguish performance levels. Outfit is more influenced by unusual ratings.",
      "Neither is a verdict on the rater. No automatic warning cutoff is applied; missing values remain unavailable.")
    if (nrow(tables$cases)) notes["cases"] <- paste(
      "These are the largest departures from expected scores, taking expected variation into account.",
      "They are discussion examples, not automatically scoring errors. Positive residuals mean higher scores than expected.",
      "Case numbers refer only to this list, not the original data rows.")
  }
  # No source object, names, row names, attributes or free-form source notes cross this boundary.
  out <- list(title = "Rater feedback", label = label, audience = audience,
    review = if (ready) "Review the assessment context before sharing or acting on this sheet." else
      "The source fit needs review. Treat these values as provisional and discuss them with the analyst before acting.",
    guidance = technical, notes = notes, tables = tables)
  class(out) <- "mfrm_rater_feedback"
  out$markdown <- mfrm_rater_feedback_render(out, html = FALSE)
  out
}

# Two-family sheets use the same recipient boundary and renderer, with their
# own targets. No ordinary fit reference or location interval is transferred.
mfrm_gmfrm_rater_feedback <- function(x, facet, rater, audience, label, max_cases, interval) {
  fit <- x$fit; cfg <- fit$config
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !identical(cfg$model, "GPCM") || !identical(cfg$method, "MML") ||
      !identical(cfg$posterior_basis, "fixed_standard_normal") ||
      !identical(cfg$step_facet, cfg$slope_facet[2L]) ||
      !isTRUE(all(cfg$facet_signs[cfg$slope_facet] == -1)) ||
      length(cfg$interaction_specs) || isTRUE(fit$population$active))
    stop("Two-family rater sheets require the native fixed-N(0,1) GPCM MML fit.", call. = FALSE)
  roles <- mfrm_gpcm_slope_roles(cfg)
  if (!facet %in% roles$SlopeOwner) stop("`facet` must name a fitted non-Person slope facet.", call. = FALSE)
  if (!rater %in% cfg$facet_levels[[facet]]) stop("`rater` must select exactly one level of `facet`.", call. = FALSE)
  params <- expand_params(fit$opt$par, build_param_sizes(cfg), cfg)
  expected <- list(locations = build_other_facet_table(cfg, fit$prep, params),
    steps = build_step_table(cfg, fit$prep, params), slopes = build_slope_table(cfg, fit$prep, params))
  actual <- list(locations = fit$facets$others, steps = fit$steps, slopes = fit$slopes)
  actual$slopes$Estimate <- actual$slopes$OptimizerEstimate
  keys <- list(locations = c("Facet", "Level", "Estimate"),
    steps = c("StepFacet", "Step", "Estimate"), slopes = c("SlopeOwner", "SlopeFacet", "Estimate"))
  for (key in names(keys)) {
    if (!is.data.frame(actual[[key]]) || !all(keys[[key]] %in% names(actual[[key]])) ||
        !isTRUE(all.equal(actual[[key]][keys[[key]]], expected[[key]][keys[[key]]],
          tolerance = 1e-10, check.attributes = FALSE)))
      stop("Saved location, step and slope values must match the two-family fitted parameters.", call. = FALSE)
  }
  location <- as.numeric(params$facets[[facet]][match(rater, cfg$facet_levels[[facet]])])
  slope <- as.numeric(params$slope_components[[facet]]$Slope[match(rater, cfg$facet_levels[[facet]])])
  if (any(!is.finite(c(location, slope))) || slope <= 0)
    stop("The selected two-family point estimates are unavailable.", call. = FALSE)
  first <- identical(facet, cfg$slope_facet[1L])
  d <- fit$prep$data; take <- which(as.character(d[[facet]]) == rater)
  other <- setdiff(cfg$slope_facet, facet)
  if (!length(take) || any(!is.finite(d$Weight)) || any(d$Weight != 1))
    stop("Two-family rater sheets require retained unit-weight rating rows.", call. = FALSE)
  score <- as.numeric(d$Score[take]); map <- fit$prep$score_map
  if (!is.data.frame(map) || !identical(as.numeric(map$OriginalScore), as.numeric(map$InternalScore)) ||
      !identical(as.numeric(map$InternalScore), as.numeric(0:(cfg$n_cat - 1L))) ||
      anyNA(score) || !all(score %in% map$InternalScore))
    stop("Retain the original zero-based two-family category coding.", call. = FALSE)
  counts <- vapply(as.numeric(map$InternalScore), function(k) sum(score == k), integer(1))
  tables <- list(
    severity = data.frame(Location = location),
    slope = data.frame(ComponentSlope = slope),
    exposure = data.frame(RatingRows = length(take), Persons = length(unique(d$Person[take])),
      ObservedCounterpartLevels = length(unique(d[[other]][take])),
      FittedCounterpartLevels = length(cfg$facet_levels[[other]])),
    uncertainty = data.frame(),
    steps = if (!first) data.frame(Step = seq_len(cfg$n_cat - 1L),
      Offset = as.numeric(params$steps_mat[match(rater, cfg$facet_levels[[facet]]), ])) else data.frame(),
    categories = data.frame(ModelScore = as.numeric(map$InternalScore), Ratings = counts,
      Percent = 100 * counts / length(take)),
    residuals = data.frame(), cases = data.frame())
  notes <- c(
    severity = paste("Increasing this location lowers expected scores when ability, both slopes and all other modeled effects are held fixed.",
      if (first) "Locations in this family sum to zero across fitted levels." else
        "Locations in this family are not centered; the ability population has mean zero and SD one.",
      "Zero is not a quality standard or necessarily the other raters' average. Location intervals are not displayed in this individual sheet."),
    slope = paste("This component multiplies the other family's slope for each rating. Together they control how category probabilities change with ability.",
      if (first) "Components in this family have geometric mean one." else
        "Components in this family are free on the fixed standard-normal ability scale; their geometric mean need not equal one.",
      "One is not a competence benchmark. A larger slope does not establish scoring accuracy or training effectiveness."),
    exposure = "Counts cover retained ratings, people and observed levels of the other facet. Unobserved combinations are not evidence of the rater's behavior there. These counts do not measure planned-assignment completion or make different raters' assignments comparable.",
    uncertainty = "No saved interval for this component slope was supplied. Location and category-step intervals are not displayed in this sheet; no interval is calculated here.",
    steps = if (first) "Category steps belong to the other facet; there is no single set of steps specific to this selected level." else
      "These category-step offsets sum to zero within this level. An adjacent-category threshold also includes both facet locations; an offset alone is not a threshold. Step intervals are unavailable.",
    categories = "Percentages describe this rater's retained assignments on the fitted category scale. Unequal use or an unused category alone does not establish a rubric or rater problem.",
    residuals = "No saved posterior response diagnostics were supplied. Ordinary Infit/Outfit reference values, cutoffs and tests are unavailable for this model.",
    cases = "No saved posterior response diagnostics were supplied; rating cases are not calculated for this sheet.")
  intervals <- mfrm_gpcm_results_inputs(fit, x$gpcm_inference)
  candidates <- lapply(intervals, function(ci) {
    if (!inherits(ci, "mfrm_slope_intervals")) return(integer(0))
    tab <- attr(ci, "diagnostics")
    if (!all(c("SlopeOwner", "SlopeLevel") %in% names(tab))) return(integer(0))
    which(tab$SlopeOwner == facet & tab$SlopeLevel == rater)
  })
  eligible <- names(candidates)[lengths(candidates) > 0L]
  if (!is.null(interval) && !interval %in% eligible)
    stop("`interval` must name a saved interval for this individual component slope.", call. = FALSE)
  if (is.null(interval) && length(eligible) > 1L)
    stop("Several saved intervals contain this component slope. Choose one with `interval`.", call. = FALSE)
  if (is.null(interval) && length(eligible) == 1L) interval <- eligible
  if (!is.null(interval)) {
    ci <- intervals[[interval]]; at <- candidates[[interval]]
    tab <- attr(ci, "diagnostics")[at, , drop = FALSE]; settings <- attr(ci, "settings")
    level <- attr(ci, "level"); adjustment <- attr(ci, "simultaneous")
    if (length(at) != 1L ||
        !all(c("Estimate", "CI_Level", "CI_Lower", "CI_Upper", "CIEligible", "ScaleReference") %in% names(tab)) ||
        !isTRUE(all.equal(as.numeric(tab$Estimate), slope, tolerance = 1e-8)) ||
        !identical(settings$scale, "standardized") || !is.null(settings$contrasts) ||
        !isTRUE(settings$method %in% c("model", "profile")) ||
        !isTRUE(adjustment %in% c("none", "bonferroni")) ||
        !is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1 ||
        !identical(as.numeric(tab$CI_Level), as.numeric(level)) ||
        !identical(as.character(tab$ScaleReference), if (first) "geometric_mean_one" else "fixed_standard_normal"))
      stop("The saved interval must retain this component's estimate, scale and method.", call. = FALSE)
    available <- isTRUE(tab$CIEligible) && all(is.finite(c(tab$CI_Lower, tab$CI_Upper))) &&
      tab$CI_Lower > 0 && tab$CI_Lower <= tab$CI_Upper
    reason <- mfrm_rater_interval_reason(available,
      if (isTRUE(tab$CIEligible) && !available) "invalid_bounds" else NULL,
      checks = attr(ci, "checks"), endpoints = attr(ci, "profile")$endpoints)
    tables$uncertainty <- data.frame(Level = as.numeric(level),
      Lower = if (available) as.numeric(tab$CI_Lower) else NA_real_,
      Upper = if (available) as.numeric(tab$CI_Upper) else NA_real_,
      Method = if (settings$method == "model") "Log-Wald" else "Profile likelihood",
      Adjustment = if (adjustment == "none") "Pointwise" else "Bonferroni",
      Available = available, Reason = reason)
    notes["uncertainty"] <- paste(
      if (available) "Saved experimental interval for this component slope only." else
        paste("The saved component-slope interval is unavailable.", reason),
      "It does not describe uncertainty in location or category steps. Coverage is not established; this is not a test of rater quality or a comparison with another rater.",
      if (adjustment == "bonferroni") "The saved Bonferroni adjustment retains the original family of targets; selecting one recipient does not change it.")
  }
  diagnostics <- x$response_diagnostics
  if (!is.null(diagnostics)) {
    mfrm_validate_response_diagnostics(fit, diagnostics)
    tab <- diagnostics$rows
    if (!all(c("InputRow", "Score", "ExpectedScore", "StandardizedResidual", "Status") %in% names(tab)) ||
        anyNA(tab$InputRow) || anyDuplicated(tab$InputRow) ||
        any(tab$InputRow != floor(tab$InputRow) | tab$InputRow < 1 | tab$InputRow > nrow(diagnostics$source_data)))
      stop("Saved response diagnostics must retain distinct original row positions.", call. = FALSE)
    tab <- tab[as.character(diagnostics$source_data[[facet]][tab$InputRow]) == rater, , drop = FALSE]
    available <- tab$Status == "available_conditional" & is.finite(tab$ExpectedScore) &
      is.finite(tab$StandardizedResidual) & is.finite(tab$Score)
    available[is.na(available)] <- FALSE
    complete <- nrow(tab) > 0 && all(available)
    tables$residuals <- data.frame(SavedRows = nrow(tab), Available = sum(available),
      Unavailable = sum(!available), NotIncluded = length(take) - nrow(tab),
      MeanResidual = if (complete) mean(tab$Score - tab$ExpectedScore) else NA_real_,
      MeanSquaredStandardizedResidual = if (complete) mean(tab$StandardizedResidual^2) else NA_real_)
    notes["residuals"] <- paste("These summaries cover only the saved rows for this recipient; rows not included in the diagnostics are counted separately.",
      "Expectations use each person's full observed record and fixed calibration. They describe the same data used to fit the model, not independent predictions.",
      "A mean squared standardized residual has no expectation-one reference, warning cutoff or p-value here. Averages are unavailable if any saved row is unresolved; calibration uncertainty is excluded.")
    rows <- which(available)
    rows <- utils::head(rows[order(-abs(tab$StandardizedResidual[rows]), tab$InputRow[rows])], max_cases)
    tables$cases <- data.frame(Case = seq_along(rows), Observed = as.numeric(tab$Score[rows]),
      Expected = as.numeric(tab$ExpectedScore[rows]), Residual = as.numeric(tab$Score[rows] - tab$ExpectedScore[rows]),
      StandardizedResidual = as.numeric(tab$StandardizedResidual[rows]))
    notes["cases"] <- paste("The largest absolute standardized residuals among available saved rows are shown; unresolved and unselected rows remain counted above.",
      "Positive residuals mean higher scores than expected. These are discussion examples, not automatic scoring errors. Case numbers are local to this sheet, not original row numbers.")
  }
  if (max_cases == 0) {
    tables$cases <- data.frame(); notes["cases"] <- "Individual rating cases were omitted by request."
  }
  numerical <- isTRUE(fit$summary$Converged) && identical(fit$opt$convergence, 0L) &&
    identical(fit$opt$optimizer_diagnostics$ConvergenceSeverity, "pass")
  out <- structure(list(title = "Two-family GPCM rater feedback", label = label, audience = audience,
    review = paste(if (numerical) "Experimental descriptive feedback; numerical convergence does not establish model adequacy or reliable inference." else
      "The source fit needs review: numerical convergence is unresolved. Treat all displayed estimates as provisional.",
      if (isTRUE(tables$residuals$Unavailable > 0))
        "Saved response diagnostics contain unresolved rows; residual averages are unavailable.",
      if (isTRUE(tables$residuals$NotIncluded > 0)) "Residuals describe a selected subset of retained ratings.",
      "Review assignment differences and assessment context before sharing or acting."),
    guidance = if (audience == "researcher") paste(
      "Two-family GPCM MML with fixed N(0,1) ability. The product of component slopes multiplies the complete adjacent-category predictor.",
      "Locations and centered step offsets are on the ability scale; posterior replicate residuals condition on all observed ratings and fixed calibration.") else
      "Use scoring tendency, response sensitivity, category use and rating examples together to guide a discussion with the assessment team. This sheet does not rank raters or establish competence or training effects.",
    notes = notes, tables = tables), class = "mfrm_rater_feedback")
  out$markdown <- mfrm_rater_feedback_render(out, html = FALSE)
  out
}

# Only fixed descriptions cross the recipient boundary, never source free text.
mfrm_rater_interval_reason <- function(available, status = NULL, checks = NULL, endpoints = NULL) {
  if (available) return("Saved interval is numerically available under its stated assumptions; coverage is not guaranteed.")
  reasons <- c(
    fixed = "This value is fixed by a model constraint; no inferential interval is supplied.",
    insufficient_cluster_rank = "There is insufficient independent-group information for the requested sandwich interval; model-based bounds are not substituted.",
    nonpositive_variance = "The estimated variance is nonpositive or nonfinite; uncertainty cannot be calculated.",
    nonfinite_bounds = "The interval bounds cannot be represented as finite numbers.",
    invalid_bounds = "The saved bounds are missing, nonfinite or invalid for this interval.",
    numerical_review_failed = "The source estimation or information checks did not pass.")
  if (length(status) == 1L && !is.na(status) && status %in% names(reasons))
    return(unname(reasons[status]))
  profile_reasons <- c(
    source_profile_failed = "Profile checks at the fitted value failed.",
    numerical_or_start_failure = "Numerical checks or agreement between starting values failed.",
    nonmonotone_profile = "The searched likelihood did not change consistently away from the fitted value.",
    endpoint_search_failed = "An interval endpoint search failed.",
    endpoint_check_failed = "An interval endpoint did not meet the numerical accuracy checks.",
    crossing_not_found = "No likelihood cutoff crossing was found within the searched range; this does not establish an unbounded interval.",
    source_likelihood_improved = "The search found a better source likelihood; the original fit needs review.")
  if (is.data.frame(endpoints) && "Status" %in% names(endpoints)) {
    known <- unique(as.character(endpoints$Status))
    known <- known[known %in% names(profile_reasons)]
    if (length(known)) return(paste(unname(profile_reasons[known]), collapse = " "))
  }
  if (is.data.frame(checks) && all(c("Check", "Passed") %in% names(checks))) {
    failed <- as.character(checks$Check[!checks$Passed %in% TRUE])
    descriptions <- c(
      "Scope" = "The saved analysis does not meet the method's supported conditions.",
      "Model identity" = "The saved model and source data did not agree.",
      "Slope identity" = "The saved slope values and model constraints did not agree.",
      "EM convergence" = "Estimation did not meet the convergence checks.",
      "Direct convergence" = "Estimation did not meet the convergence checks.",
      "Category support" = "The observed categories do not support this interval approximation.",
      "Joint information" = "The joint information matrix could not support an unregularized uncertainty calculation.",
      "Likelihood identity" = "The saved likelihood did not match its numerical reevaluation.",
      "Score workspace" = "The numerical verification exceeded its configured memory allowance.",
      "Score derivatives" = "The numerical derivatives could not be verified.",
      "Score derivative stability" = "The numerical derivatives were not sufficiently stable.",
      "Stationary local maximum" = "The fitted solution did not meet the local optimum checks.",
      "Refined quadrature information" = "The more accurate integration did not provide usable information.",
      "Quadrature sensitivity" = "The numerical integration was not sufficiently stable.")
    known <- failed[failed %in% names(descriptions)]
    if (length(known)) return(paste(unique(unname(descriptions[known])), collapse = " "))
    if (length(failed)) return("The saved source checks did not pass; ask the analyst to review the full interval result.")
    if (nrow(checks)) return("The variance or interval bounds could not be calculated despite passing the source checks.")
  }
  "The saved interval is marked unavailable or lacks usable diagnostics; ask the analyst to review the full interval result."
}

mfrm_rater_feedback_render <- function(x, html = TRUE) {
  headings <- c(severity = "Scoring tendency", exposure = "Ratings included",
    uncertainty = "Uncertainty in scoring tendency", fit = "Consistency with the model",
    categories = "Category use", cases = "Ratings to discuss", slope = "Response sensitivity",
    steps = "Category-step offsets", residuals = "Saved response residuals",
    study_context = "Study and score meaning")
  product <- "slope" %in% names(x$tables)
  if (product) headings["uncertainty"] <- "Uncertainty in the component slope"
  headings <- headings[setdiff(names(x$tables), "interpretation_notes")]
  display <- lapply(x$tables, function(tab) {
    tab[] <- lapply(tab, function(col) {
      if (is.numeric(col)) ifelse(is.na(col), "Not available", format(round(col, 3), trim = TRUE)) else col
    })
    tab
  })
  names(display$severity) <- if (product) "Location" else c("Severity", "Model reference")
  names(display$exposure) <- if (product) c("Rating rows", "People rated", "Observed counterpart levels", "Fitted counterpart levels") else
    c("Rating rows", "People rated", "Weight total")
  if (nrow(display$uncertainty)) {
    # The same reason is shown beside the table and retained in its CSV fields.
    display$uncertainty$Reason <- NULL
    display$uncertainty$Level <- paste0(100 * x$tables$uncertainty$Level, "%")
    if (!product) display$uncertainty$Method <- ifelse(display$uncertainty$Method == "model", "Model-based", "Sandwich")
    display$uncertainty$Available <- ifelse(display$uncertainty$Available, "Available", "Unavailable")
    names(display$uncertainty) <- c("Confidence", "Lower", "Upper", "Calculation", if (product) "Adjustment", "Status")
  }
  names(display$categories) <- c("Model score", "Ratings", "Percent", if (!product) "Weight total",
    if ("OriginalScore" %in% names(display$categories)) "Original score")
  if (product) {
    names(display$slope) <- "Component slope"
    if (nrow(display$steps)) names(display$steps) <- c("Step", "Offset")
    if (nrow(display$residuals)) names(display$residuals) <- c("Saved rows", "Available", "Unavailable", "Not included",
      "Mean residual", "Mean squared standardized residual")
    # These are one-recipient summaries; vertical fields also fit printed cards.
    for (key in c("uncertainty", "residuals")) if (nrow(display[[key]]))
      display[[key]] <- data.frame(Measure = names(display[[key]]),
        Value = as.character(unlist(display[[key]][1, ], use.names = FALSE)))
  }
  if (nrow(display$cases)) names(display$cases) <- c("Case", "Observed", "Expected", "Residual", "Standardized residual")
  if (!html) return(paste(c(paste0("# ", x$title), html_escape(x$label), x$review, x$guidance,
    unlist(lapply(names(headings), function(key) c(paste0("## ", headings[key]),
      x$notes[key], if (nrow(display[[key]])) mfrm_report_markdown_table(display[[key]]))))), collapse = "\n\n"))
  sections <- vapply(names(headings), function(key) {
    table <- if (nrow(display[[key]])) sub("<table>", paste0('<table aria-label="', headings[key], '">'),
      gsub("<th>", '<th scope="col">', dataframe_to_html_table(display[[key]]), fixed = TRUE), fixed = TRUE) else ""
    bars <- if (key == "categories") paste0('<div class="bars" aria-hidden="true">',
      paste0('<div class="bar-row"><span>', html_escape(display$categories[["Model score"]]),
        '</span><span class="track"><span style="width:',
        format(x$tables$categories$Percent, scientific = FALSE, trim = TRUE),
        '%"></span></span></div>', collapse = ""), '</div>') else ""
    paste0('<section><h2>', headings[key], '</h2><p>', html_escape(x$notes[key]),
      '</p>', bars, '<div class="table-scroll" tabindex="0" role="region" aria-label="', headings[key], '">', table, '</div></section>')
  }, character(1))
  paste0('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1"><title>Rater feedback</title>',
    '<style>body{font:16px/1.5 system-ui,sans-serif;color:#192d3c;background:#f2f5f7;margin:0}',
    'main{max-width:920px;margin:auto;padding:28px}h1{margin:0;font-size:2rem}h2{font-size:1.1rem;margin:0 0 8px}',
    '.label{font-size:1.2rem;overflow-wrap:anywhere}.review{border-left:4px solid #245d75;padding:10px 14px;background:#e6eef3}',
    '.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}section{min-width:0;background:white;border:1px solid #bdcbd3;border-radius:8px;padding:16px}',
    'p{margin:0 0 12px}.table-scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:.86rem}',
    'th,td{padding:5px;text-align:left;border-bottom:1px solid #ccd5db}th{background:#edf2f5}',
    '.table-scroll:focus-visible{outline:3px solid #245d75;outline-offset:2px}',
    '.bars{margin:10px 0}.bar-row{display:flex;gap:10px;align-items:center}.bar-row>span:first-child{width:2em}',
    '.track{flex:1;background:#e5ebef;height:12px}.track>span{display:block;background:#245d75;height:12px}',
    '@media(max-width:700px){main{padding:14px}.grid{grid-template-columns:1fr}}',
    '@media print{@page{size:A4;margin:12mm}body{background:white;font-size:9pt;line-height:1.3}main{padding:0;max-width:none}',
    'h1{font-size:18pt}h2{font-size:11pt}.grid{display:grid;grid-template-columns:1fr 1fr;gap:8px}section{break-inside:avoid;margin:0;padding:8px}',
    'p{margin-bottom:6px}table{font-size:8pt}.bars{display:none}.table-scroll{overflow:visible}.table-scroll:focus-visible{outline:none}}',
    '</style></head><body><main><h1>', html_escape(x$title), '</h1><p class="label">', html_escape(x$label),
    '</p><p class="review">', html_escape(x$review), '</p><p>', html_escape(x$guidance),
    '</p><div class="grid">', paste(sections, collapse = ""), '</div></main></body></html>')
}

#' @export
print.mfrm_rater_feedback <- function(x, ...) {
  cat(x$title, "\n", x$label, "\n", x$review, "\n", sep = "")
  print(x$tables$severity, row.names = FALSE)
  cat("\n", x$notes[["uncertainty"]], "\n", sep = "")
  if (nrow(x$tables$uncertainty)) print(x$tables$uncertainty, row.names = FALSE)
  invisible(x)
}

mfrm_rater_feedback_html <- function(report) {
  html <- mfrm_rater_feedback_render(report)
  path <- tempfile("mfrmr_rater_", fileext = ".html")
  writeLines(enc2utf8(html), path, useBytes = TRUE)
  structure(list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
    report = report, html = html), class = "mfrm_report_html")
}
