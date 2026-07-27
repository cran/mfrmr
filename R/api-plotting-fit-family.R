# Plot methods centered on fitted many-facet Rasch objects.

.mfrm_fit_plot_readiness <- function(fit) {
  convergence <- mfrm_convergence_state(fit)
  data_review <- fit$data_review %||% list()
  review_status <- as.data.frame(
    data_review$status %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(review_status) == 0L && !is.null(fit$prep)) {
    data_review <- tryCatch(
      build_mfrm_data_review(fit$prep),
      error = function(e) list()
    )
    review_status <- as.data.frame(
      data_review$status %||% data.frame(),
      stringsAsFactors = FALSE
    )
  }
  domain_status <- function(domain, default = "not_assessed") {
    if (nrow(review_status) == 0L ||
        !all(c("Domain", "Status") %in% names(review_status))) return(default)
    value <- review_status$Status[match(domain, review_status$Domain)]
    if (length(value) == 0L || is.na(value[1L]) || !nzchar(value[1L])) {
      default
    } else {
      as.character(value[1L])
    }
  }
  convergence_severity <- tolower(trimws(as.character(
    convergence$severity %||% "review"
  )[1L]))
  numerical_status <- if (isTRUE(convergence$inference_ready)) {
    "pass"
  } else if (convergence_severity %in% c("fail", "error")) {
    "fail"
  } else {
    "review"
  }
  data_status <- domain_status("Data")
  design_status <- domain_status("Design")
  stability_status <- domain_status("Stability")
  interpretation_ready <- identical(numerical_status, "pass") &&
    identical(data_status, "pass") &&
    identical(design_status, "pass_linked") &&
    identical(stability_status, "pass")
  interpretation_status <- if (interpretation_ready) {
    "ready_for_diagnostic_interpretation"
  } else {
    "review_only"
  }
  statuses <- data.frame(
    Domain = c("Numerical", "Data", "Design", "Stability", "Plot"),
    Status = c(
      numerical_status,
      data_status,
      design_status,
      stability_status,
      interpretation_status
    ),
    stringsAsFactors = FALSE
  )
  detail <- if (interpretation_ready) {
    paste0(
      "Fit-level numerical, data-support, connectivity, and stability gates ",
      "passed. Treat this display as diagnostic evidence, not automatic ",
      "publication approval."
    )
  } else {
    paste0(
      "Review-only display: Numerical=", numerical_status,
      ", Data=", data_status,
      ", Design=", design_status,
      ", Stability=", stability_status,
      ". Inspect `summary(fit)$readiness` and `fit$data_review` before ",
      "substantive or cross-subset interpretation."
    )
  }
  list(
    status = interpretation_status,
    ready = interpretation_ready,
    detail = detail,
    table = statuses
  )
}

.mfrm_attach_plot_readiness <- function(plot_object, readiness) {
  plot_object$data$fit_readiness <- readiness$table
  plot_object$data$interpretation_status <- readiness$status
  plot_object$data$interpretation_note <- readiness$detail
  if (!isTRUE(readiness$ready)) {
    subtitle <- as.character(plot_object$data$subtitle %||% "")
    review_label <- "REVIEW ONLY"
    plot_object$data$subtitle <- if (length(subtitle) > 0L && nzchar(subtitle[1L])) {
      paste0(review_label, " - ", subtitle[1L])
    } else {
      review_label
    }
  }
  plot_object
}

.mfrm_plot_display_title <- function(plot_data, title = NULL,
                                     fallback = "mfrmr plot") {
  display_title <- title
  if (is.null(display_title) || length(display_title) == 0L ||
      is.na(display_title[1L]) || !nzchar(as.character(display_title[1L]))) {
    display_title <- plot_data$title %||% fallback
  }
  display_title <- as.character(display_title[1L])
  if (!nzchar(display_title)) display_title <- as.character(fallback[1L])
  review_only <- identical(
    as.character(plot_data$interpretation_status %||% "")[1L],
    "review_only"
  )
  if (isTRUE(review_only) && !startsWith(display_title, "REVIEW ONLY")) {
    display_title <- paste0("REVIEW ONLY - ", display_title)
  }
  display_title
}

step_index_from_label <- function(step_labels) {
  step_labels <- as.character(step_labels)
  idx <- suppressWarnings(as.integer(gsub("[^0-9]+", "", step_labels)))
  if (all(is.na(idx))) return(seq_along(step_labels))
  if (any(is.na(idx))) {
    fill_start <- max(idx, na.rm = TRUE)
    idx[is.na(idx)] <- fill_start + seq_len(sum(is.na(idx)))
  }
  idx
}

.build_wright_group_density <- function(x, group, group_data = NULL) {
  # Returns a long-format data.frame with columns Group, Theta, Density
  # for each subgroup of persons. The group label can come from:
  #  * a column in `group_data` (data.frame with a Person column;
  #    typical use for DIF screening when fit_mfrm was called without
  #    the Group column);
  #  * a column in `x$prep$data` (when the user passed it through the
  #    fit), or
  #  * a vector aligned with `x$facets$person$Person` (length match).
  # NA / extreme rows are dropped before density estimation.
  person_tbl <- as.data.frame(x$facets$person, stringsAsFactors = FALSE)
  if (nrow(person_tbl) == 0L) {
    stop("No person estimates available for group overlay.", call. = FALSE)
  }
  joined <- NULL
  if (is.character(group) && length(group) == 1L) {
    if (!is.null(group_data) && is.data.frame(group_data) &&
        "Person" %in% names(group_data) &&
        group %in% names(group_data)) {
      lookup <- unique(group_data[, c("Person", group), drop = FALSE])
      names(lookup)[2] <- "Group"
      joined <- merge(person_tbl[, c("Person", "Estimate")], lookup,
                      by = "Person", all.x = TRUE)
    } else {
      obs_df <- as.data.frame(x$prep$data %||% NULL, stringsAsFactors = FALSE)
      if (group %in% names(obs_df) && "Person" %in% names(obs_df)) {
        lookup <- unique(obs_df[, c("Person", group), drop = FALSE])
        names(lookup)[2] <- "Group"
        joined <- merge(person_tbl[, c("Person", "Estimate")], lookup,
                        by = "Person", all.x = TRUE)
      } else {
        stop(sprintf(
          "Group column '%s' not found. Pass `group_data = <data.frame>` ",
          group
        ), sprintf(
          "with a Person column and a '%s' column, or include the column in fit_mfrm()'s data argument.",
          group
        ), call. = FALSE)
      }
    }
  } else if (length(group) == nrow(person_tbl)) {
    joined <- data.frame(
      Person = person_tbl$Person,
      Estimate = person_tbl$Estimate,
      Group = as.character(group),
      stringsAsFactors = FALSE
    )
  } else {
    stop("`group` must be a column name (string) or a vector aligned ",
         "with fit$facets$person rows.", call. = FALSE)
  }
  joined <- joined[is.finite(joined$Estimate) & !is.na(joined$Group), ,
                    drop = FALSE]
  groups <- unique(as.character(joined$Group))
  if (length(groups) < 2L) {
    stop("Group overlay requires at least two distinct group levels.",
         call. = FALSE)
  }
  out <- do.call(rbind, lapply(groups, function(g) {
    sub <- joined$Estimate[as.character(joined$Group) == g]
    if (length(sub) < 5L) {
      return(NULL)
    }
    d <- stats::density(sub, n = 256L)
    data.frame(Group = g, Theta = d$x, Density = d$y,
               stringsAsFactors = FALSE)
  }))
  if (is.null(out) || nrow(out) == 0L) {
    stop("No subgroup had at least 5 finite person estimates.",
         call. = FALSE)
  }
  out
}

build_step_curve_spec <- function(x) {
  step_tbl <- x$steps
  if (is.null(step_tbl) || nrow(step_tbl) == 0 || !"Estimate" %in% names(step_tbl)) {
    stop("Step estimates are required for pathway/CCC plots.")
  }
  step_tbl <- tibble::as_tibble(step_tbl)

  model <- toupper(as.character(x$config$model[1]))
  rating_min <- suppressWarnings(as.numeric(x$prep$rating_min %||% NA_real_))
  if (!is.finite(rating_min)) {
    rating_min <- suppressWarnings(min(as.numeric(x$prep$data$Score), na.rm = TRUE))
  }
  if (!is.finite(rating_min)) rating_min <- 0

  n_cat <- suppressWarnings(as.integer(x$config$n_cat[1]))
  if (!is.finite(n_cat) || n_cat < 2) {
    n_cat <- max(2L, nrow(step_tbl) + 1L)
  }
  categories <- rating_min + 0:(n_cat - 1L)

  groups <- list()
  step_points <- tibble::tibble()

  if (model == "RSM") {
    if (!"Step" %in% names(step_tbl)) {
      step_tbl$Step <- paste0("Step_", seq_len(nrow(step_tbl)))
    }
    ord <- order(step_index_from_label(step_tbl$Step))
    tau <- as.numeric(step_tbl$Estimate[ord])
    groups[["Common"]] <- list(name = "Common", step_cum = c(0, cumsum(tau)), tau = tau)
    step_points <- tibble::tibble(
      CurveGroup = "Common",
      Step = as.character(step_tbl$Step[ord]),
      StepIndex = seq_along(tau),
      Threshold = tau
    )
  } else {
    if (!all(c("StepFacet", "Step", "Estimate") %in% names(step_tbl))) {
      stop("PCM step table must include StepFacet, Step, and Estimate columns.")
    }
    step_facet <- as.character(x$config$step_facet[1])
    ordered_levels <- x$prep$levels[[step_facet]]
    if (is.null(ordered_levels)) {
      ordered_levels <- unique(as.character(step_tbl$StepFacet))
    }
    slope_lookup <- NULL
    if (model == "GPCM") {
      slope_tbl <- tibble::as_tibble(x$slopes %||% tibble::tibble())
      if (!all(c("SlopeFacet", "Estimate") %in% names(slope_tbl))) {
        stop("GPCM pathway/CCC plots require a `fit$slopes` table with `SlopeFacet` and `Estimate` columns.")
      }
      slope_lookup <- stats::setNames(
        suppressWarnings(as.numeric(slope_tbl$Estimate)),
        as.character(slope_tbl$SlopeFacet)
      )
    }
    for (lvl in ordered_levels) {
      sub <- step_tbl[as.character(step_tbl$StepFacet) == as.character(lvl), , drop = FALSE]
      if (nrow(sub) == 0) next
      ord <- order(step_index_from_label(sub$Step))
      tau <- as.numeric(sub$Estimate[ord])
      slope_val <- if (model == "GPCM") {
        sval <- unname(slope_lookup[[as.character(lvl)]])
        if (!is.finite(sval) || sval <= 0) {
          stop("GPCM pathway/CCC plots require finite positive slopes for every `StepFacet` level.")
        }
        sval
      } else {
        1
      }
      groups[[as.character(lvl)]] <- list(
        name = as.character(lvl),
        step_cum = c(0, cumsum(tau)),
        tau = tau,
        slope = slope_val
      )
      step_points <- dplyr::bind_rows(
        step_points,
        tibble::tibble(
          CurveGroup = as.character(lvl),
          Step = as.character(sub$Step[ord]),
          StepIndex = seq_along(tau),
          Threshold = tau,
          Slope = slope_val
        )
      )
    }
  }

  if (length(groups) == 0) stop("No step groups available for pathway/CCC plots.")

  list(
    model = model,
    step_facet = if (model == "RSM") NA_character_ else as.character(x$config$step_facet[1]),
    categories = categories,
    rating_min = rating_min,
    n_cat = n_cat,
    groups = groups,
    step_points = step_points
  )
}

empty_pathway_fit_payload <- function(status = "not_available", message = "") {
  fit_measures <- data.frame(
    Facet = character(),
    Level = character(),
    Measure = numeric(),
    SE = numeric(),
    Infit = numeric(),
    Outfit = numeric(),
    InfitZSTD = numeric(),
    OutfitZSTD = numeric(),
    FitStatus = character(),
    Underfit = logical(),
    Overfit = logical(),
    ReviewReason = character(),
    stringsAsFactors = FALSE
  )
  curve_fit_status <- data.frame(
    CurveGroup = character(),
    Facet = character(),
    Level = character(),
    Measure = numeric(),
    SE = numeric(),
    Infit = numeric(),
    Outfit = numeric(),
    InfitZSTD = numeric(),
    OutfitZSTD = numeric(),
    FitStatus = character(),
    Underfit = logical(),
    Overfit = logical(),
    ReviewReason = character(),
    MatchedFitRow = logical(),
    stringsAsFactors = FALSE
  )
  curve_fit_annotations <- data.frame(
    AnnotationType = character(),
    CurveGroup = character(),
    Facet = character(),
    Level = character(),
    X = numeric(),
    Y = numeric(),
    Label = character(),
    Measure = numeric(),
    SE = numeric(),
    FitStatus = character(),
    Underfit = logical(),
    Overfit = logical(),
    ReviewReason = character(),
    stringsAsFactors = FALSE
  )
  list(
    fit_measures = fit_measures,
    fit_status = data.frame(Facet = character(), FitStatus = character(), Rows = integer(), stringsAsFactors = FALSE),
    curve_fit_status = curve_fit_status,
    curve_fit_annotations = curve_fit_annotations,
    fit_measure_status = data.frame(
      Available = identical(status, "available"),
      Status = status,
      Message = as.character(message %||% ""),
      stringsAsFactors = FALSE
    )
  )
}

pathway_fit_payload <- function(x,
                                diagnostics = NULL,
                                include_fit_measures = TRUE,
                                curve_groups = character(),
                                step_facet = NA_character_,
                                endpoint_labels = data.frame()) {
  curve_groups <- unique(as.character(curve_groups %||% character(0)))
  curve_groups <- curve_groups[nzchar(curve_groups)]
  if (!isTRUE(include_fit_measures)) {
    out <- empty_pathway_fit_payload("not_requested", "`include_fit_measures = FALSE`.")
    out$curve_fit_status <- data.frame(
      CurveGroup = curve_groups,
      Facet = NA_character_,
      Level = NA_character_,
      Measure = NA_real_,
      SE = NA_real_,
      Infit = NA_real_,
      Outfit = NA_real_,
      InfitZSTD = NA_real_,
      OutfitZSTD = NA_real_,
      FitStatus = "not_requested",
      Underfit = FALSE,
      Overfit = FALSE,
      ReviewReason = "Fit measures not requested",
      MatchedFitRow = FALSE,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  fm <- tryCatch(
    fit_measures_table(
      x,
      diagnostics = diagnostics,
      include_person = FALSE,
      sort_by = "facet",
      top_n = Inf
    ),
    error = function(e) e
  )
  if (inherits(fm, "error")) {
    out <- empty_pathway_fit_payload("error", conditionMessage(fm))
    out$curve_fit_status <- data.frame(
      CurveGroup = curve_groups,
      Facet = NA_character_,
      Level = NA_character_,
      Measure = NA_real_,
      SE = NA_real_,
      Infit = NA_real_,
      Outfit = NA_real_,
      InfitZSTD = NA_real_,
      OutfitZSTD = NA_real_,
      FitStatus = "not_available",
      Underfit = FALSE,
      Overfit = FALSE,
      ReviewReason = conditionMessage(fm),
      MatchedFitRow = FALSE,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  fit_tbl <- as.data.frame(fm$table %||% data.frame(), stringsAsFactors = FALSE)
  status_tbl <- as.data.frame(fm$status_summary %||% data.frame(), stringsAsFactors = FALSE)
  keep_cols <- c(
    "Facet", "Level", "Measure", "SE", "Infit", "Outfit", "InfitZSTD",
    "OutfitZSTD", "FitStatus", "Underfit", "Overfit", "ReviewReason"
  )
  for (nm in keep_cols) {
    if (!nm %in% names(fit_tbl)) {
      fit_tbl[[nm]] <- if (nm %in% c("Facet", "Level", "FitStatus", "ReviewReason")) {
        NA_character_
      } else if (nm %in% c("Underfit", "Overfit")) {
        NA
      } else {
        NA_real_
      }
    }
  }
  fit_tbl <- fit_tbl[, keep_cols, drop = FALSE]

  candidate_tbl <- fit_tbl
  step_facet <- as.character(step_facet %||% NA_character_)
  if (!is.na(step_facet) && nzchar(step_facet) && "Facet" %in% names(candidate_tbl)) {
    narrowed <- candidate_tbl[as.character(candidate_tbl$Facet) == step_facet, , drop = FALSE]
    if (nrow(narrowed) > 0L) candidate_tbl <- narrowed
  }
  match_idx <- match(curve_groups, as.character(candidate_tbl$Level))
  curve_fit_status <- data.frame(
    CurveGroup = curve_groups,
    Facet = as.character(candidate_tbl$Facet[match_idx]),
    Level = as.character(candidate_tbl$Level[match_idx]),
    Measure = suppressWarnings(as.numeric(candidate_tbl$Measure[match_idx])),
    SE = suppressWarnings(as.numeric(candidate_tbl$SE[match_idx])),
    Infit = suppressWarnings(as.numeric(candidate_tbl$Infit[match_idx])),
    Outfit = suppressWarnings(as.numeric(candidate_tbl$Outfit[match_idx])),
    InfitZSTD = suppressWarnings(as.numeric(candidate_tbl$InfitZSTD[match_idx])),
    OutfitZSTD = suppressWarnings(as.numeric(candidate_tbl$OutfitZSTD[match_idx])),
    FitStatus = as.character(candidate_tbl$FitStatus[match_idx]),
    Underfit = as.logical(candidate_tbl$Underfit[match_idx]),
    Overfit = as.logical(candidate_tbl$Overfit[match_idx]),
    ReviewReason = as.character(candidate_tbl$ReviewReason[match_idx]),
    MatchedFitRow = !is.na(match_idx),
    stringsAsFactors = FALSE
  )
  curve_fit_status$FitStatus[!curve_fit_status$MatchedFitRow] <- "not_matched"
  curve_fit_status$Underfit[is.na(curve_fit_status$Underfit)] <- FALSE
  curve_fit_status$Overfit[is.na(curve_fit_status$Overfit)] <- FALSE
  curve_fit_status$ReviewReason[!curve_fit_status$MatchedFitRow] <- "No fit-measure row matched this pathway curve group."

  endpoints <- as.data.frame(endpoint_labels %||% data.frame(), stringsAsFactors = FALSE)
  fit_annotations <- curve_fit_status[
    curve_fit_status$FitStatus %in% c("underfit", "overfit", "mixed"),
    ,
    drop = FALSE
  ]
  curve_fit_annotations <- empty_pathway_fit_payload()$curve_fit_annotations
  if (nrow(fit_annotations) > 0L && nrow(endpoints) > 0L &&
      all(c("CurveGroup", "Theta", "ExpectedScore") %in% names(endpoints))) {
    endpoint_idx <- match(fit_annotations$CurveGroup, as.character(endpoints$CurveGroup))
    ok <- !is.na(endpoint_idx)
    if (any(ok)) {
      fit_annotations <- fit_annotations[ok, , drop = FALSE]
      endpoint_idx <- endpoint_idx[ok]
      curve_fit_annotations <- data.frame(
        AnnotationType = "fit_status",
        CurveGroup = fit_annotations$CurveGroup,
        Facet = fit_annotations$Facet,
        Level = fit_annotations$Level,
        X = suppressWarnings(as.numeric(endpoints$Theta[endpoint_idx])),
        Y = suppressWarnings(as.numeric(endpoints$ExpectedScore[endpoint_idx])),
        Label = paste0(fit_annotations$CurveGroup, ": ", fit_annotations$FitStatus),
        Measure = fit_annotations$Measure,
        SE = fit_annotations$SE,
        FitStatus = fit_annotations$FitStatus,
        Underfit = fit_annotations$Underfit,
        Overfit = fit_annotations$Overfit,
        ReviewReason = fit_annotations$ReviewReason,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    fit_measures = fit_tbl,
    fit_status = status_tbl,
    curve_fit_status = curve_fit_status,
    curve_fit_annotations = curve_fit_annotations,
    fit_measure_status = data.frame(
      Available = TRUE,
      Status = "available",
      Message = "Fit measures are available for pathway annotation and custom plotting.",
      stringsAsFactors = FALSE
    )
  )
}

build_curve_tables <- function(curve_spec, theta_grid) {
  prob_tables <- list()
  exp_tables <- list()
  idx_prob <- 1L
  idx_exp <- 1L
  for (g in names(curve_spec$groups)) {
    grp <- curve_spec$groups[[g]]
    slope_val <- if (identical(curve_spec$model, "GPCM")) {
      as.numeric(grp$slope %||% NA_real_)
    } else {
      1
    }
    if (!is.finite(slope_val) || slope_val <= 0) {
      stop("Category-curve information requires finite positive slopes.")
    }
    probs <- if (identical(curve_spec$model, "GPCM")) {
      category_prob_gpcm(
        eta = theta_grid,
        step_cum_mat = matrix(grp$step_cum, nrow = 1L),
        criterion_idx = rep(1L, length(theta_grid)),
        slopes = slope_val,
        slope_idx = rep(1L, length(theta_grid))
      )
    } else {
      category_prob_rsm(theta_grid, grp$step_cum)
    }
    k_vals <- as.numeric(curve_spec$categories)
    expected <- as.numeric(probs %*% matrix(k_vals, ncol = 1))
    second <- as.numeric(probs %*% matrix(k_vals^2, ncol = 1))
    score_variance <- pmax(second - expected^2, 0)
    information <- (slope_val^2) * score_variance
    exp_tables[[idx_exp]] <- tibble::tibble(
      Theta = theta_grid,
      ExpectedScore = expected,
      ScoreVariance = score_variance,
      Information = information,
      Slope = slope_val,
      Model = curve_spec$model,
      CurveGroup = grp$name
    )
    idx_exp <- idx_exp + 1L
    for (k in seq_len(ncol(probs))) {
      category_value <- as.numeric(curve_spec$categories[k])
      category_information <- (slope_val^2) * probs[, k] * (category_value - expected)^2
      prob_tables[[idx_prob]] <- tibble::tibble(
        Theta = theta_grid,
        Probability = probs[, k],
        ExpectedScore = expected,
        ScoreVariance = score_variance,
        Information = information,
        CategoryInformation = category_information,
        CategoryInformationShare = ifelse(information > 0, category_information / information, NA_real_),
        Slope = slope_val,
        Model = curve_spec$model,
        Category = as.character(curve_spec$categories[k]),
        CurveGroup = grp$name
      )
      idx_prob <- idx_prob + 1L
    }
  }
  list(
    probabilities = dplyr::bind_rows(prob_tables),
    expected = dplyr::bind_rows(exp_tables)
  )
}

compute_se_for_plot <- function(x, ci_level = 0.95, diagnostics = NULL) {
  ci_level <- suppressWarnings(as.numeric(ci_level[1]))
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(diagnostics)) {
    measures <- as.data.frame(diagnostics$measures %||% data.frame(), stringsAsFactors = FALSE)
    if (!all(c("Facet", "Level", "Estimate") %in% names(measures))) {
      stop("`diagnostics$measures` must contain Facet, Level, and Estimate columns.", call. = FALSE)
    }
    if (!"SE" %in% names(measures)) measures$SE <- NA_real_
    keep <- intersect(
      c(
        "Facet", "Level", "Estimate", "SE", "SE_Method", "PrecisionTier",
        "SupportsFormalInference", "SEUse", "CIBasis", "CIUse"
      ),
      names(measures)
    )
    se_tbl <- measures[, keep, drop = FALSE]
    se_tbl$Facet <- as.character(se_tbl$Facet)
    se_tbl$Level <- as.character(se_tbl$Level)
    se_tbl$Estimate <- suppressWarnings(as.numeric(se_tbl$Estimate))
    se_tbl$SE <- suppressWarnings(as.numeric(se_tbl$SE))
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    se_tbl$CI_Lower <- ifelse(
      is.finite(se_tbl$Estimate) & is.finite(se_tbl$SE),
      se_tbl$Estimate - z * se_tbl$SE,
      NA_real_
    )
    se_tbl$CI_Upper <- ifelse(
      is.finite(se_tbl$Estimate) & is.finite(se_tbl$SE),
      se_tbl$Estimate + z * se_tbl$SE,
      NA_real_
    )
    se_tbl$CI_Level <- ci_level
    se_tbl$Measure_Source <- "diagnostics$measures"
    return(se_tbl)
  }
  tryCatch({
    obs_df <- compute_obs_table(x)
    facet_cols <- x$config$source_columns$facets
    if (is.null(facet_cols) || length(facet_cols) == 0) {
      facet_cols <- x$config$facet_names
    }
    se_tbl <- calc_facet_se(obs_df, facet_cols)
    estimates <- tibble::as_tibble(x$facets$others) |>
      dplyr::transmute(
        Facet = as.character(.data$Facet),
        Level = as.character(.data$Level),
        Estimate = as.numeric(.data$Estimate)
      )
    se_tbl <- se_tbl |>
      dplyr::mutate(
        Facet = as.character(.data$Facet),
        Level = as.character(.data$Level)
      ) |>
      dplyr::left_join(estimates, by = c("Facet", "Level"))
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    se_tbl$CI_Lower <- ifelse(
      is.finite(se_tbl$Estimate) & is.finite(se_tbl$SE),
      se_tbl$Estimate - z * se_tbl$SE,
      NA_real_
    )
    se_tbl$CI_Upper <- ifelse(
      is.finite(se_tbl$Estimate) & is.finite(se_tbl$SE),
      se_tbl$Estimate + z * se_tbl$SE,
      NA_real_
    )
    se_tbl$CI_Level <- ci_level
    se_tbl$SE_Method <- "Observation-table information"
    se_tbl$Measure_Source <- "fit + observation-table information"
    as.data.frame(se_tbl, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

.wright_native_legend_spec <- function(plot_core,
                                       show_ci = FALSE,
                                       ci_level = 0.95) {
  loc <- as.data.frame(plot_core$locations %||% data.frame(),
                       stringsAsFactors = FALSE)
  has_type <- function(value) {
    nrow(loc) > 0L && "PlotType" %in% names(loc) &&
      any(as.character(loc$PlotType) == value)
  }
  logical_column <- function(name) {
    if (name %in% names(loc)) {
      value <- as.logical(loc[[name]])
      value[is.na(value)] <- FALSE
      value
    } else {
      rep(FALSE, nrow(loc))
    }
  }
  has_interval <- if (all(c("OriginalCI_Lower", "OriginalCI_Upper") %in%
                          names(loc))) {
    is.finite(loc$OriginalCI_Lower) & is.finite(loc$OriginalCI_Upper)
  } else if ("SE" %in% names(loc)) {
    is.finite(loc$SE) & loc$SE > 0
  } else {
    rep(FALSE, nrow(loc))
  }
  suppressed <- logical_column("CISuppressed")
  clipped <- logical_column("CIClipped")
  boundary_end <- if ("BoundaryEnd" %in% names(loc)) {
    as.character(loc$BoundaryEnd)
  } else {
    rep("none", nrow(loc))
  }
  ci_available <- isTRUE(show_ci) && any(has_interval & !suppressed)
  clipped_available <- isTRUE(show_ci) && any(
    has_interval & !suppressed & clipped
  )
  boundary_available <- any(boundary_end %in% c("lower", "upper"))
  boundary_ci <- isTRUE(show_ci) && any(has_interval & suppressed)

  key <- c(
    "person_density",
    if (has_type("Facet level")) "facet_level",
    if (has_type("Step threshold")) "step_threshold",
    if (nrow(plot_core$group_summary %||% data.frame()) > 0L) "iqr_range",
    if (ci_available) "ci_whisker",
    if (clipped_available) "ci_continues",
    if (boundary_available) "boundary_endpoint"
  )
  label <- vapply(key, function(value) switch(
    value,
    person_density = "Person density",
    facet_level = "Facet level",
    step_threshold = "Step threshold",
    iqr_range = "IQR / range",
    ci_whisker = sprintf("%g%% mfrmr CI", round(100 * ci_level)),
    ci_continues = "CI continues beyond ruler",
    boundary_endpoint = if (boundary_ci) {
      "Boundary level / CI beyond ruler"
    } else {
      "Boundary level at ruler end"
    }
  ), character(1))
  role <- vapply(key, function(value) switch(
    value,
    person_density = "distribution",
    facet_level = "location",
    step_threshold = "location",
    iqr_range = "summary",
    ci_whisker = "uncertainty",
    ci_continues = "uncertainty",
    boundary_endpoint = "stability"
  ), character(1))
  aesthetic <- vapply(key, function(value) switch(
    value,
    person_density = "histogram",
    facet_level = "point",
    step_threshold = "point",
    iqr_range = "range",
    ci_whisker = "whisker",
    ci_continues = "endpoint",
    boundary_endpoint = "endpoint"
  ), character(1))
  data.frame(
    key = key,
    label = label,
    role = role,
    aesthetic = aesthetic,
    stringsAsFactors = FALSE
  )
}

.wright_map_legend <- function(plot_core,
                               style,
                               show_ci = FALSE,
                               ci_level = 0.95) {
  wright_style <- as.character(plot_core$wright_style %||% "native")[1L]
  if (!identical(wright_style, "native")) {
    ci_label <- if (isTRUE(show_ci)) {
      sprintf("%g%% mfrmr CI", round(100 * ci_level))
    } else {
      character(0)
    }
    return(new_plot_legend(
      label = c("Person density", "Facet levels", "Step thresholds", ci_label),
      role = c(
        "distribution", "location", "location",
        if (length(ci_label) > 0L) "uncertainty"
      ),
      aesthetic = c(
        "histogram", "points", "points",
        if (length(ci_label) > 0L) "whisker"
      ),
      value = c(
        style$fill_muted, style$accent_tertiary, style$accent_secondary,
        if (length(ci_label) > 0L) style$accent_tertiary
      )
    ))
  }

  spec <- .wright_native_legend_spec(
    plot_core,
    show_ci = show_ci,
    ci_level = ci_level
  )
  value <- vapply(spec$key, function(key) switch(
    key,
    person_density = style$fill_muted,
    facet_level = style$accent_tertiary,
    step_threshold = style$accent_secondary,
    iqr_range = style$foreground,
    ci_whisker = style$accent_tertiary,
    ci_continues = style$accent_tertiary,
    boundary_endpoint = style$warn
  ), character(1))
  new_plot_legend(
    label = spec$label,
    role = spec$role,
    aesthetic = spec$aesthetic,
    value = value
  )
}

build_wright_map_data <- function(x,
                                  top_n = 30L,
                                  se_tbl = NULL,
                                  include_steps = TRUE,
                                  wright_style = c("native", "facets_style"),
                                  renderer = NULL,
                                  category_labels = NULL,
                                  rows_per_logit = 2L,
                                  wright_range = NULL,
                                  extreme_placement = c("ends", "estimate"),
                                  persons_per_star = NULL) {
  top_n <- suppressWarnings(as.numeric(top_n[1]))
  if (length(top_n) != 1L || is.na(top_n) || top_n <= 0) {
    stop("`top_n` must be a positive number or `Inf`.", call. = FALSE)
  }
  if (is.finite(top_n)) {
    top_n <- as.integer(floor(top_n))
    if (top_n < 1L) stop("`top_n` must be at least 1 or `Inf`.", call. = FALSE)
  }
  wright_style <- match_wright_style(wright_style, renderer = renderer)
  person_tbl <- tibble::as_tibble(x$facets$person)
  person_tbl <- person_tbl[is.finite(person_tbl$Estimate), , drop = FALSE]
  if (nrow(person_tbl) == 0) stop("Person estimates are not available for Wright map.")

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (!all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available for Wright map.")
  }
  facet_tbl <- facet_tbl[is.finite(facet_tbl$Estimate), , drop = FALSE]

  step_points <- if (isTRUE(include_steps)) {
    tryCatch(build_step_curve_spec(x)$step_points, error = function(e) tibble::tibble())
  } else {
    tibble::tibble()
  }
  step_tbl <- if (nrow(step_points) > 0) {
    categories <- seq(x$prep$rating_min, x$prep$rating_max)
    score_labels <- wright_score_labels(
      x,
      categories = categories,
      category_labels = category_labels
    )
    step_index <- suppressWarnings(as.integer(step_points$StepIndex))
    if (length(step_index) != nrow(step_points) || any(!is.finite(step_index))) {
      step_index <- step_index_from_label(step_points$Step)
    }
    lower_label <- score_labels$DisplayLabel[step_index]
    upper_label <- score_labels$DisplayLabel[step_index + 1L]
    transition_label <- ifelse(
      !is.na(lower_label) & !is.na(upper_label),
      paste0(lower_label, " -> ", upper_label),
      as.character(step_points$Step)
    )
    tibble::tibble(
      PlotType = "Step threshold",
      Group = if ("CurveGroup" %in% names(step_points)) paste0("Step:", step_points$CurveGroup) else "Step",
      Label = transition_label,
      Step = as.character(step_points$Step),
      StepIndex = step_index,
      Estimate = step_points$Threshold
    )
  } else {
    tibble::tibble()
  }

  facet_pts <- facet_tbl |>
    dplyr::transmute(
      PlotType = "Facet level",
      Group = .data$Facet,
      Label = .data$Level,
      Estimate = .data$Estimate
    )
  if (!is.null(se_tbl) && is.data.frame(se_tbl) && nrow(se_tbl) > 0 &&
      all(c("Facet", "Level") %in% names(se_tbl)) &&
      "SE" %in% names(se_tbl)) {
    se_join <- se_tbl[, intersect(
      c(
        "Facet", "Level", "SE",
        "CI_Level", "SE_Method", "PrecisionTier", "SupportsFormalInference",
        "SEUse", "CIBasis", "CIUse", "Measure_Source"
      ),
      names(se_tbl)
    ), drop = FALSE]
    se_join$Level <- as.character(se_join$Level)
    facet_pts$Group <- as.character(facet_pts$Group)
    facet_pts$Label <- as.character(facet_pts$Label)
    facet_pts <- merge(
      facet_pts,
      se_join,
      by.x = c("Group", "Label"),
      by.y = c("Facet", "Level"),
      all.x = TRUE,
      sort = FALSE
    )
    if ("SE" %in% names(facet_pts)) {
      ci_level_value <- if ("CI_Level" %in% names(facet_pts)) {
        suppressWarnings(as.numeric(facet_pts$CI_Level))
      } else {
        rep(0.95, nrow(facet_pts))
      }
      ci_level_value[!is.finite(ci_level_value) | ci_level_value <= 0 | ci_level_value >= 1] <- 0.95
      z <- stats::qnorm(1 - (1 - ci_level_value) / 2)
      facet_pts$CI_Lower <- ifelse(
        is.finite(facet_pts$Estimate) & is.finite(facet_pts$SE),
        facet_pts$Estimate - z * facet_pts$SE,
        NA_real_
      )
      facet_pts$CI_Upper <- ifelse(
        is.finite(facet_pts$Estimate) & is.finite(facet_pts$SE),
        facet_pts$Estimate + z * facet_pts$SE,
        NA_real_
      )
    }
  }
  point_tbl <- dplyr::bind_rows(facet_pts, step_tbl)
  if (nrow(point_tbl) == 0) stop("No facet/step locations available for Wright map.")
  boundary_levels <- as.data.frame(
    x$data_review$boundary_levels %||% data.frame(),
    stringsAsFactors = FALSE
  )
  boundary_key <- character()
  if (nrow(boundary_levels) > 0L &&
      all(c("Facet", "Level") %in% names(boundary_levels))) {
    boundary_key <- unique(paste(
      as.character(boundary_levels$Facet),
      as.character(boundary_levels$Level),
      sep = "\r"
    ))
  }
  point_key <- paste(
    as.character(point_tbl$Group),
    as.character(point_tbl$Label),
    sep = "\r"
  )
  point_tbl$BoundarySeparated <- point_tbl$PlotType == "Facet level" &
    point_key %in% boundary_key

  native_display_scale <- NULL
  if (identical(wright_style, "native") && is.null(wright_range) &&
      any(point_tbl$BoundarySeparated)) {
    # The native map shares the FACETS-style automatic range policy only for
    # diagnosed boundary separation. Ordinary native maps retain their prior
    # estimate-based range exactly.
    native_display_scale <- build_wright_display_scale(
      fit = x,
      person_values = person_tbl$Estimate,
      location_tbl = point_tbl,
      category_labels = NULL,
      rows_per_logit = 2L,
      wright_range = NULL,
      include_steps = include_steps
    )
  }
  location_totals <- c(
    `Facet level` = sum(point_tbl$PlotType %in% "Facet level"),
    `Step threshold` = sum(point_tbl$PlotType %in% "Step threshold")
  )

  if (identical(wright_style, "native") && nrow(point_tbl) > top_n) {
    # Category transitions are part of the scale definition, so keep every
    # step label even when facet-level locations are compacted by `top_n`.
    step_keep <- point_tbl$PlotType %in% "Step threshold"
    step_rows <- point_tbl[step_keep, , drop = FALSE]
    facet_rows <- point_tbl[!step_keep, , drop = FALSE]
    facet_rows <- facet_rows[
      order(abs(facet_rows$Estimate), decreasing = TRUE),
      ,
      drop = FALSE
    ]
    facet_limit <- max(0L, top_n - nrow(step_rows))
    facet_rows <- utils::head(facet_rows, facet_limit)
    point_tbl <- dplyr::bind_rows(facet_rows, step_rows)
  }
  location_shown <- c(
    `Facet level` = sum(point_tbl$PlotType %in% "Facet level"),
    `Step threshold` = sum(point_tbl$PlotType %in% "Step threshold")
  )
  retention <- tibble::tibble(
    Component = c(names(location_totals), "All locations"),
    Shown = as.integer(c(location_shown, sum(location_shown))),
    Total = as.integer(c(location_totals, sum(location_totals))),
    Omitted = as.integer(c(location_totals - location_shown, sum(location_totals - location_shown))),
    RequestedTopN = rep(top_n, 3L)
  )
  retention$Complete <- retention$Omitted == 0L
  omitted_facets <- retention$Omitted[retention$Component == "Facet level"]
  retention_note <- if (length(omitted_facets) == 1L && omitted_facets > 0L) {
    paste0(
      "Compact native view: showing ",
      retention$Shown[retention$Component == "Facet level"],
      " of ", retention$Total[retention$Component == "Facet level"],
      " facet locations; ", omitted_facets,
      " omitted. Use `top_n = Inf` for the complete final map."
    )
  } else {
    ""
  }

  group_levels <- unique(point_tbl$Group)
  point_tbl <- point_tbl |>
    dplyr::mutate(Group = factor(.data$Group, levels = group_levels)) |>
    dplyr::group_by(.data$Group) |>
    dplyr::arrange(.data$Estimate, .by_group = TRUE) |>
    dplyr::mutate(
      XBase = as.numeric(.data$Group),
      X = if (dplyr::n() == 1) .data$XBase else .data$XBase + seq(-0.18, 0.18, length.out = dplyr::n())
    ) |>
    dplyr::ungroup()

  hist_data <- graphics::hist(person_tbl$Estimate, breaks = "FD", plot = FALSE)
  y_range <- range(c(point_tbl$Estimate, person_tbl$Estimate), finite = TRUE)
  if (!all(is.finite(y_range)) || diff(y_range) <= sqrt(.Machine$double.eps)) {
    center <- mean(c(point_tbl$Estimate, person_tbl$Estimate), na.rm = TRUE)
    y_range <- center + c(-0.5, 0.5)
  }
  auto_range_policy <- "all_fitted_locations"
  if (!is.null(native_display_scale) && identical(
    native_display_scale$auto_range_policy,
    "boundary_levels_at_ends"
  )) {
    y_range <- c(native_display_scale$lower, native_display_scale$upper)
    auto_range_policy <- native_display_scale$auto_range_policy
  }

  point_tbl$OriginalEstimate <- suppressWarnings(as.numeric(point_tbl$Estimate))
  point_tbl$BelowRange <- point_tbl$OriginalEstimate < y_range[1]
  point_tbl$AboveRange <- point_tbl$OriginalEstimate > y_range[2]
  point_tbl$DisplayEstimate <- pmin(
    y_range[2],
    pmax(y_range[1], point_tbl$OriginalEstimate)
  )
  point_tbl$DisplayLabel <- ifelse(
    point_tbl$BelowRange | point_tbl$AboveRange,
    paste0("(", as.character(point_tbl$Label), ")"),
    as.character(point_tbl$Label)
  )
  point_tbl <- tibble::as_tibble(add_wright_ci_display_metadata(
    point_tbl,
    lower = y_range[1],
    upper = y_range[2]
  ))

  step_label_tbl <- point_tbl |>
    dplyr::filter(.data$PlotType == "Step threshold")
  facet_label_tbl <- point_tbl |>
    dplyr::filter(.data$PlotType != "Step threshold") |>
    dplyr::group_by(.data$Group) |>
    dplyr::slice(c(1L, dplyr::n())) |>
    dplyr::ungroup() |>
    dplyr::distinct(.data$Group, .data$Label, .keep_all = TRUE)
  facet_label_limit <- max(0L, top_n - nrow(step_label_tbl))
  if (nrow(facet_label_tbl) > facet_label_limit) {
    facet_label_tbl <- utils::head(facet_label_tbl, facet_label_limit)
  }
  label_tbl <- dplyr::bind_rows(facet_label_tbl, step_label_tbl) |>
    dplyr::distinct(.data$Group, .data$Label, .keep_all = TRUE)

  group_summary <- point_tbl |>
    dplyr::group_by(.data$Group, .data$PlotType) |>
    dplyr::summarise(
      Min = min(.data$Estimate, na.rm = TRUE),
      Q1 = stats::quantile(.data$Estimate, 0.25, na.rm = TRUE, names = FALSE),
      Median = stats::median(.data$Estimate, na.rm = TRUE),
      Q3 = stats::quantile(.data$Estimate, 0.75, na.rm = TRUE, names = FALSE),
      Max = max(.data$Estimate, na.rm = TRUE),
      DisplayMin = min(.data$DisplayEstimate, na.rm = TRUE),
      DisplayQ1 = stats::quantile(.data$DisplayEstimate, 0.25, na.rm = TRUE, names = FALSE),
      DisplayMedian = stats::median(.data$DisplayEstimate, na.rm = TRUE),
      DisplayQ3 = stats::quantile(.data$DisplayEstimate, 0.75, na.rm = TRUE, names = FALSE),
      DisplayMax = max(.data$DisplayEstimate, na.rm = TRUE),
      N = dplyr::n(),
      XBase = mean(.data$XBase, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      TargetGap = .data$Median - mean(person_tbl$Estimate, na.rm = TRUE),
      DisplayTargetGap = .data$DisplayMedian - mean(person_tbl$Estimate, na.rm = TRUE)
    )

  person_stats <- tibble::tibble(
    N = nrow(person_tbl),
    Mean = mean(person_tbl$Estimate, na.rm = TRUE),
    Median = stats::median(person_tbl$Estimate, na.rm = TRUE),
    SD = stats::sd(person_tbl$Estimate, na.rm = TRUE)
  )

  display_settings <- tibble::tibble(
    Renderer = if (identical(wright_style, "native")) "native" else "facets",
    LowerLogit = y_range[1],
    UpperLogit = y_range[2],
    AutoRangePolicy = auto_range_policy,
    BoundaryLevelsAtEnds = sum(point_tbl$BoundaryEnd != "none"),
    CIClippedCount = sum(point_tbl$CIClipped),
    BoundaryCIEndpointCount = sum(point_tbl$CISuppressed),
    CIDisplayPolicy = paste0(
      "Exact intervals remain in CI_Lower/CI_Upper; clipped intervals use ",
      "endpoint markers and boundary-separated intervals may be omitted from ",
      "the ruler."
    )
  )

  out <- list(
    title = "Wright Map",
    wright_style = wright_style,
    renderer = if (identical(wright_style, "facets_style")) "facets" else "native",
    visual_contract = if (identical(wright_style, "facets_style")) {
      "FACETS Table 6-style visual layout; not FACETS numerical equivalence"
    } else {
      "mfrmr native Wright-map layout"
    },
    person = person_tbl,
    person_hist = hist_data,
    person_stats = person_stats,
    locations = point_tbl,
    label_points = label_tbl,
    group_summary = group_summary,
    group_levels = group_levels,
    y_range = y_range,
    display_settings = display_settings,
    label_limit = top_n,
    retention = retention,
    retention_note = retention_note
  )
  if (identical(wright_style, "facets_style")) {
    out$facets_style <- build_facets_style_wright_data(
      fit = x,
      plot_data = out,
      category_labels = category_labels,
      rows_per_logit = rows_per_logit,
      wright_range = wright_range,
      extreme_placement = extreme_placement,
      persons_per_star = persons_per_star,
      include_steps = include_steps
    )
    out$y_range <- c(
      out$facets_style$settings$LowerLogit[1],
      out$facets_style$settings$UpperLogit[1]
    )
    out$display_settings <- out$facets_style$settings
  }
  out
}

draw_wright_map <- function(plot_data,
                            title = NULL,
                            palette = NULL,
                            label_angle = 45,
                            show_ci = FALSE,
                            ci_level = 0.95) {
  facets_style <- identical(
    as.character(plot_data$wright_style %||% "native"),
    "facets_style"
  )
  default_title <- if (isTRUE(facets_style)) {
    "FACETS-style Wright map"
  } else {
    "Facet and step locations"
  }
  display_title <- .mfrm_plot_display_title(
    plot_data,
    title = if (is.null(title)) default_title else title,
    fallback = default_title
  )
  if (isTRUE(facets_style)) {
    return(draw_wright_facets_style(
      plot_data,
      title = display_title,
      palette = palette,
      show_ci = show_ci,
      ci_level = ci_level
    ))
  }
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet_level = "#1b9e77",
      step_threshold = "#d95f02",
      person_hist = "gray80",
      grid = "#ececec",
      range = "#94a3b8",
      iqr = "#334155"
    )
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.15, 2.55))
  graphics::par(mgp = c(2.6, 0.85, 0))

  loc <- plot_data$locations
  yr <- plot_data$y_range
  loc_y <- if ("DisplayEstimate" %in% names(loc)) {
    suppressWarnings(as.numeric(loc$DisplayEstimate))
  } else {
    suppressWarnings(as.numeric(loc$Estimate))
  }
  boundary_end <- if ("BoundaryEnd" %in% names(loc)) {
    as.character(loc$BoundaryEnd)
  } else {
    rep("none", nrow(loc))
  }
  ci_clipped <- if ("CIClipped" %in% names(loc)) {
    value <- as.logical(loc$CIClipped)
    value[is.na(value)] <- FALSE
    value
  } else {
    rep(FALSE, nrow(loc))
  }
  ci_suppressed <- if ("CISuppressed" %in% names(loc)) {
    value <- as.logical(loc$CISuppressed)
    value[is.na(value)] <- FALSE
    value
  } else {
    rep(FALSE, nrow(loc))
  }
  retention_note <- as.character(plot_data$retention_note %||% "")
  footer_notes <- character()
  if (length(retention_note) > 0L && nzchar(retention_note[1L])) {
    footer_notes <- c(footer_notes, retention_note[1L])
  }
  boundary_at_ends <- any(boundary_end %in% c("lower", "upper"))
  if (boundary_at_ends) {
    footer_notes <- c(footer_notes, if (isTRUE(show_ci) && any(ci_suppressed)) {
      paste0(
        "Boundary-separated levels use end triangles and their intervals are ",
        "omitted from the ruler; inspect OriginalEstimate and CI_Lower/CI_Upper."
      )
    } else {
      "Boundary-separated levels use end triangles; inspect OriginalEstimate for exact values."
    })
  } else if (isTRUE(show_ci) && any(ci_clipped)) {
    footer_notes <- c(
      footer_notes,
      "CI endpoint triangles identify clipping; exact bounds remain in CI_Lower/CI_Upper."
    )
  }
  footer_lines <- unlist(lapply(
    footer_notes,
    function(note) strwrap(note, width = 118L)
  ), use.names = FALSE)
  bottom_margin <- if (length(footer_lines) > 1L) {
    9.3 + 0.65 * (length(footer_lines) - 1L)
  } else {
    8.8
  }
  pretty_y <- pretty(yr, n = 6)
  hist_data <- plot_data$person_hist
  hist_max <- max(hist_data$counts, na.rm = TRUE)
  if (!is.finite(hist_max) || hist_max <= 0) hist_max <- 1

  graphics::par(mar = c(bottom_margin, 4.6, 3.2, 0.8))
  graphics::plot(
    x = c(0, hist_max * 1.12),
    y = yr,
    type = "n",
    xlab = "Persons",
    ylab = "Logit scale",
    main = "Person distribution",
    yaxt = "n"
  )
  graphics::axis(2, at = pretty_y, las = 1)
  graphics::abline(h = pretty_y, col = pal["grid"], lty = 1)
  for (i in seq_along(hist_data$counts)) {
    graphics::rect(
      xleft = 0,
      ybottom = hist_data$breaks[i],
      xright = hist_data$counts[i],
      ytop = hist_data$breaks[i + 1],
      col = grDevices::adjustcolor(pal["person_hist"], alpha.f = 0.85),
      border = "white"
    )
  }
  if (nrow(plot_data$person_stats) > 0) {
    graphics::abline(h = plot_data$person_stats$Median[1], col = grDevices::adjustcolor("gray25", alpha.f = 0.9), lty = 2, lwd = 1.2)
    graphics::abline(h = plot_data$person_stats$Mean[1], col = grDevices::adjustcolor("gray45", alpha.f = 0.9), lty = 3, lwd = 1.2)
    graphics::mtext(
      sprintf(
        "N=%d  Mean=%.2f  Median=%.2f",
        plot_data$person_stats$N[1],
        plot_data$person_stats$Mean[1],
        plot_data$person_stats$Median[1]
      ),
      side = 3,
      line = 0.2,
      cex = 0.78,
      adj = 0
    )
  }

  xr <- c(0.5, length(plot_data$group_levels) + 0.75)
  graphics::par(mar = c(bottom_margin, 2.4, 3.2, 2.8))
  graphics::plot(
    x = loc$X,
    y = loc_y,
    type = "n",
    xlim = xr,
    ylim = yr,
    xaxt = "n",
    xlab = "",
    yaxt = "n",
    ylab = "",
    main = display_title
  )
  for (i in seq_along(plot_data$group_levels)) {
    if (i %% 2 == 1) {
      graphics::rect(
        xleft = i - 0.5,
        ybottom = yr[1] - 10,
        xright = i + 0.5,
        ytop = yr[2] + 10,
        col = grDevices::adjustcolor(pal["grid"], alpha.f = 0.18),
        border = NA
      )
    }
  }
  graphics::abline(h = pretty_y, col = pal["grid"], lty = 1)
  graphics::abline(h = 0, col = grDevices::adjustcolor("gray35", alpha.f = 0.9), lty = 2)
  group_summary <- plot_data$group_summary
  if (!is.null(group_summary) && nrow(group_summary) > 0) {
    summary_value <- function(display_name, original_name) {
      if (display_name %in% names(group_summary)) {
        group_summary[[display_name]]
      } else {
        group_summary[[original_name]]
      }
    }
    range_cols <- ifelse(
      group_summary$PlotType == "Step threshold",
      grDevices::adjustcolor(pal["step_threshold"], alpha.f = 0.42),
      grDevices::adjustcolor(pal["range"], alpha.f = 0.72)
    )
    iqr_cols <- ifelse(
      group_summary$PlotType == "Step threshold",
      grDevices::adjustcolor(pal["step_threshold"], alpha.f = 0.78),
      grDevices::adjustcolor(pal["iqr"], alpha.f = 0.92)
    )
    graphics::segments(
      x0 = group_summary$XBase,
      y0 = summary_value("DisplayMin", "Min"),
      x1 = group_summary$XBase,
      y1 = summary_value("DisplayMax", "Max"),
      col = range_cols,
      lwd = 1.6
    )
    graphics::segments(
      x0 = group_summary$XBase,
      y0 = summary_value("DisplayQ1", "Q1"),
      x1 = group_summary$XBase,
      y1 = summary_value("DisplayQ3", "Q3"),
      col = iqr_cols,
      lwd = 5.2
    )
    graphics::segments(
      x0 = group_summary$XBase - 0.13,
      y0 = summary_value("DisplayMedian", "Median"),
      x1 = group_summary$XBase + 0.13,
      y1 = summary_value("DisplayMedian", "Median"),
      col = grDevices::adjustcolor("gray15", alpha.f = 0.95),
      lwd = 1.4
    )
  }
  graphics::axis(4, at = pretty_y, labels = pretty_y, las = 1)
  x_at <- seq_along(plot_data$group_levels)
  x_lab <- truncate_axis_label(plot_data$group_levels, width = 16L)
  draw_rotated_x_labels(
    at = x_at,
    labels = x_lab,
    srt = label_angle,
    cex = 0.82,
    line_offset = 0.085
  )
  cols <- ifelse(loc$PlotType == "Step threshold", pal["step_threshold"], pal["facet_level"])
  pch <- ifelse(loc$PlotType == "Step threshold", 17, 16)
  boundary_color <- "#9A3412"
  boundary_lower <- boundary_end == "lower"
  boundary_upper <- boundary_end == "upper"
  pch[boundary_lower] <- 25
  pch[boundary_upper] <- 24
  cols[boundary_lower | boundary_upper] <- boundary_color
  if (isTRUE(show_ci)) {
    if (all(c("OriginalCI_Lower", "OriginalCI_Upper",
              "DisplayCI_Lower", "DisplayCI_Upper") %in% names(loc))) {
      ci_ok <- is.finite(loc$OriginalCI_Lower) & is.finite(loc$OriginalCI_Upper)
      ci_lo <- loc$DisplayCI_Lower
      ci_hi <- loc$DisplayCI_Upper
    } else if ("SE" %in% names(loc)) {
      z <- stats::qnorm(1 - (1 - ci_level) / 2)
      ci_ok <- is.finite(loc$SE) & loc$SE > 0
      ci_lo <- pmin(yr[2], pmax(yr[1], loc$Estimate - z * loc$SE))
      ci_hi <- pmin(yr[2], pmax(yr[1], loc$Estimate + z * loc$SE))
    } else {
      ci_ok <- rep(FALSE, nrow(loc))
      ci_lo <- ci_hi <- rep(NA_real_, nrow(loc))
    }
    regular_ci <- ci_ok & !ci_suppressed
    if (any(regular_ci)) {
      graphics::arrows(
        x0 = loc$X[regular_ci],
        y0 = ci_lo[regular_ci],
        x1 = loc$X[regular_ci],
        y1 = ci_hi[regular_ci],
        angle = 90,
        code = 3,
        length = 0.04,
        col = grDevices::adjustcolor(cols[regular_ci], alpha.f = 0.6),
        lwd = 1.2
      )
      clipped_lower <- regular_ci &
        if ("CIClippedLower" %in% names(loc)) loc$CIClippedLower else FALSE
      clipped_upper <- regular_ci &
        if ("CIClippedUpper" %in% names(loc)) loc$CIClippedUpper else FALSE
      if (any(clipped_lower)) {
        graphics::points(
          loc$X[clipped_lower], rep(yr[1], sum(clipped_lower)),
          pch = 25, bg = "white", col = cols[clipped_lower], cex = 0.78
        )
      }
      if (any(clipped_upper)) {
        graphics::points(
          loc$X[clipped_upper], rep(yr[2], sum(clipped_upper)),
          pch = 24, bg = "white", col = cols[clipped_upper], cex = 0.78
        )
      }
    }
  }
  graphics::points(
    loc$X,
    loc_y,
    pch = pch,
    col = cols,
    bg = ifelse(boundary_lower | boundary_upper, "white",
                grDevices::adjustcolor(cols, alpha.f = 0.22)),
    cex = ifelse(loc$PlotType == "Step threshold", 1.08, 1)
  )
  label_tbl <- plot_data$label_points
  if (!is.null(label_tbl) && nrow(label_tbl) > 0) {
    place_left <- label_tbl$X > xr[1] + 0.72 * diff(xr)
    graphics::text(
      x = ifelse(place_left, label_tbl$X - 0.06, label_tbl$X + 0.06),
      y = if ("DisplayEstimate" %in% names(label_tbl)) {
        label_tbl$DisplayEstimate
      } else {
        label_tbl$Estimate
      },
      labels = truncate_axis_label(
        if ("DisplayLabel" %in% names(label_tbl)) {
          label_tbl$DisplayLabel
        } else {
          label_tbl$Label
        },
        width = 14L
      ),
      pos = ifelse(place_left, 2, 4),
      cex = 0.72,
      xpd = NA,
      col = grDevices::adjustcolor("gray20", alpha.f = 0.95)
    )
  }
  legend_spec <- .wright_native_legend_spec(
    plot_data,
    show_ci = show_ci,
    ci_level = ci_level
  )
  legend_pch <- vapply(legend_spec$key, function(key) switch(
    key,
    person_density = 15L,
    facet_level = 16L,
    step_threshold = 17L,
    iqr_range = NA_integer_,
    ci_whisker = NA_integer_,
    ci_continues = 24L,
    boundary_endpoint = 24L
  ), integer(1))
  legend_lty <- ifelse(
    legend_spec$key %in% c("iqr_range", "ci_whisker"),
    1,
    0
  )
  legend_lwd <- ifelse(
    legend_spec$key == "iqr_range",
    4,
    ifelse(legend_spec$key == "ci_whisker", 1.2, NA_real_)
  )
  legend_col <- vapply(legend_spec$key, function(key) switch(
    key,
    person_density = unname(pal["person_hist"]),
    facet_level = unname(pal["facet_level"]),
    step_threshold = unname(pal["step_threshold"]),
    iqr_range = unname(pal["iqr"]),
    ci_whisker = unname(grDevices::adjustcolor(pal["facet_level"], alpha.f = 0.6)),
    ci_continues = unname(pal["facet_level"]),
    boundary_endpoint = boundary_color
  ), character(1))
  graphics::legend(
    "topright",
    legend = legend_spec$label,
    pch = legend_pch,
    pt.bg = ifelse(
      legend_spec$key %in% c("ci_continues", "boundary_endpoint"),
      "white",
      legend_col
    ),
    lty = legend_lty,
    lwd = legend_lwd,
    col = legend_col,
    bty = "n",
    cex = 0.9
  )
  if (length(footer_lines) > 0L) {
    footer_col <- rep("#9A3412", length(footer_lines))
    for (i in seq_along(footer_lines)) {
      graphics::mtext(
        footer_lines[i],
        side = 1,
        line = 7.4 + 0.65 * (i - 1L),
        adj = 0,
        cex = 0.68,
        col = footer_col[i]
      )
    }
  }
}

build_pathway_map_data <- function(x,
                                   theta_range = c(-6, 6),
                                   theta_points = 241L,
                                   diagnostics = NULL,
                                   include_fit_measures = TRUE) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  cats <- as.character(curve_spec$categories)
  step_df <- curve_spec$step_points |>
    dplyr::mutate(
      PathY = curve_spec$rating_min + .data$StepIndex - 0.5,
      ThresholdLabel = ifelse(
        .data$StepIndex < length(cats),
        paste0(cats[pmax(1L, .data$StepIndex)], "/", cats[pmin(length(cats), .data$StepIndex + 1L)]),
        as.character(.data$Step)
      )
    )
  endpoint_labels <- curve_tbl$expected |>
    dplyr::group_by(.data$CurveGroup) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
  dominance_regions <- curve_tbl$probabilities |>
    dplyr::group_by(.data$CurveGroup, .data$Theta) |>
    dplyr::slice_max(.data$Probability, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$CurveGroup, .data$Theta) |>
    dplyr::group_by(.data$CurveGroup) |>
    dplyr::mutate(
      Region = cumsum(dplyr::coalesce(.data$Category != dplyr::lag(.data$Category), TRUE))
    ) |>
    dplyr::group_by(.data$CurveGroup, .data$Category, .data$Region) |>
    dplyr::summarise(
      ThetaStart = min(.data$Theta, na.rm = TRUE),
      ThetaEnd = max(.data$Theta, na.rm = TRUE),
      ThetaMid = mean(range(.data$Theta, na.rm = TRUE)),
      .groups = "drop"
    )
  pathway_long_expected <- data.frame(
    Layer = "expected_score",
    CurveGroup = as.character(curve_tbl$expected$CurveGroup),
    Theta = suppressWarnings(as.numeric(curve_tbl$expected$Theta)),
    Value = suppressWarnings(as.numeric(curve_tbl$expected$ExpectedScore)),
    ValueName = "ExpectedScore",
    Category = NA_character_,
    Step = NA_character_,
    StepIndex = NA_integer_,
    Label = as.character(curve_tbl$expected$CurveGroup),
    DisplayedByDefault = TRUE,
    Model = as.character(curve_tbl$expected$Model),
    Slope = suppressWarnings(as.numeric(curve_tbl$expected$Slope)),
    stringsAsFactors = FALSE
  )
  pathway_long_steps <- data.frame(
    Layer = "step_threshold",
    CurveGroup = as.character(step_df$CurveGroup),
    Theta = suppressWarnings(as.numeric(step_df$Threshold)),
    Value = suppressWarnings(as.numeric(step_df$PathY)),
    ValueName = "ThresholdPathY",
    Category = as.character(step_df$ThresholdLabel),
    Step = as.character(step_df$Step),
    StepIndex = suppressWarnings(as.integer(step_df$StepIndex)),
    Label = as.character(step_df$ThresholdLabel),
    DisplayedByDefault = TRUE,
    Model = curve_spec$model,
    Slope = if ("Slope" %in% names(step_df)) suppressWarnings(as.numeric(step_df$Slope)) else 1,
    stringsAsFactors = FALSE
  )
  pathway_long_endpoints <- data.frame(
    Layer = "endpoint_label",
    CurveGroup = as.character(endpoint_labels$CurveGroup),
    Theta = suppressWarnings(as.numeric(endpoint_labels$Theta)),
    Value = suppressWarnings(as.numeric(endpoint_labels$ExpectedScore)),
    ValueName = "ExpectedScore",
    Category = NA_character_,
    Step = NA_character_,
    StepIndex = NA_integer_,
    Label = as.character(endpoint_labels$CurveGroup),
    DisplayedByDefault = TRUE,
    Model = as.character(endpoint_labels$Model),
    Slope = suppressWarnings(as.numeric(endpoint_labels$Slope)),
    stringsAsFactors = FALSE
  )
  pathway_long <- rbind(pathway_long_expected, pathway_long_steps, pathway_long_endpoints)
  rownames(pathway_long) <- NULL

  step_annotations <- data.frame(
    AnnotationType = "step_threshold",
    CurveGroup = as.character(step_df$CurveGroup),
    Facet = if (!is.na(curve_spec$step_facet) && nzchar(curve_spec$step_facet)) curve_spec$step_facet else NA_character_,
    Level = ifelse(as.character(step_df$CurveGroup) == "Common", NA_character_, as.character(step_df$CurveGroup)),
    X = suppressWarnings(as.numeric(step_df$Threshold)),
    Y = suppressWarnings(as.numeric(step_df$PathY)),
    Label = as.character(step_df$ThresholdLabel),
    Measure = suppressWarnings(as.numeric(step_df$Threshold)),
    SE = NA_real_,
    FitStatus = NA_character_,
    Underfit = FALSE,
    Overfit = FALSE,
    ReviewReason = "Step threshold location",
    stringsAsFactors = FALSE
  )
  endpoint_annotations <- data.frame(
    AnnotationType = "endpoint_label",
    CurveGroup = as.character(endpoint_labels$CurveGroup),
    Facet = NA_character_,
    Level = as.character(endpoint_labels$CurveGroup),
    X = suppressWarnings(as.numeric(endpoint_labels$Theta)),
    Y = suppressWarnings(as.numeric(endpoint_labels$ExpectedScore)),
    Label = as.character(endpoint_labels$CurveGroup),
    Measure = NA_real_,
    SE = NA_real_,
    FitStatus = NA_character_,
    Underfit = FALSE,
    Overfit = FALSE,
    ReviewReason = "Curve endpoint label",
    stringsAsFactors = FALSE
  )
  fit_payload <- pathway_fit_payload(
    x,
    diagnostics = diagnostics,
    include_fit_measures = include_fit_measures,
    curve_groups = unique(as.character(curve_tbl$expected$CurveGroup)),
    step_facet = curve_spec$step_facet,
    endpoint_labels = endpoint_labels
  )
  pathway_annotations <- rbind(
    step_annotations,
    endpoint_annotations,
    fit_payload$curve_fit_annotations
  )
  rownames(pathway_annotations) <- NULL
  list(
    title = "Expected-score pathway (theta to score)",
    expected = curve_tbl$expected,
    steps = step_df,
    endpoint_labels = endpoint_labels,
    dominance_regions = dominance_regions,
    pathway_long = pathway_long,
    pathway_annotations = pathway_annotations,
    fit_measures = fit_payload$fit_measures,
    fit_status = fit_payload$fit_status,
    curve_fit_status = fit_payload$curve_fit_status,
    fit_measure_status = fit_payload$fit_measure_status,
    score_range = range(curve_tbl$expected$ExpectedScore, finite = TRUE)
  )
}

draw_pathway_map <- function(plot_data, title = NULL, palette = NULL) {
  exp_df <- plot_data$expected
  groups <- unique(exp_df$CurveGroup)
  defaults <- stats::setNames(
    grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)],
    groups
  )
  cols <- resolve_palette(palette = palette, defaults = defaults)
  y_rng <- range(exp_df$ExpectedScore, finite = TRUE)
  y_breaks <- pretty(y_rng, n = 6)
  dominance_regions <- plot_data$dominance_regions
  use_strips <- !is.null(dominance_regions) && nrow(dominance_regions) > 0 && length(groups) <= 6
  band_pad <- if (isTRUE(use_strips)) diff(y_rng) * (0.08 + 0.05 * length(groups)) else 0
  graphics::plot(
    x = range(exp_df$Theta, finite = TRUE),
    y = c(y_rng[1] - band_pad, y_rng[2]),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Expected score",
    main = .mfrm_plot_display_title(
      plot_data,
      title = title,
      fallback = "Expected-score pathway"
    )
  )
  graphics::abline(h = y_breaks, col = grDevices::adjustcolor("#d9dde3", alpha.f = 0.8), lty = 1)
  if (isTRUE(use_strips)) {
    band_h <- max(diff(y_rng) * 0.035, 0.08)
    band_gap <- band_h * 0.35
    band_top <- y_rng[1] - band_h * 0.45
    group_index <- stats::setNames(seq_along(groups), groups)
    for (i in seq_len(nrow(dominance_regions))) {
      grp <- as.character(dominance_regions$CurveGroup[i])
      row <- group_index[[grp]]
      y_top <- band_top - (row - 1) * (band_h + band_gap)
      y_bottom <- y_top - band_h
      graphics::rect(
        xleft = dominance_regions$ThetaStart[i],
        ybottom = y_bottom,
        xright = dominance_regions$ThetaEnd[i],
        ytop = y_top,
        col = grDevices::adjustcolor(cols[grp], alpha.f = 0.14),
        border = "white"
      )
      graphics::text(
        x = dominance_regions$ThetaMid[i],
        y = (y_bottom + y_top) / 2,
        labels = as.character(dominance_regions$Category[i]),
        cex = 0.72,
        col = cols[grp]
      )
    }
    graphics::text(
      x = rep(min(exp_df$Theta, na.rm = TRUE), length(groups)),
      y = band_top - (seq_along(groups) - 1) * (band_h + band_gap) - band_h / 2,
      labels = truncate_axis_label(groups, width = 14L),
      pos = 2,
      xpd = NA,
      cex = 0.75,
      col = cols[groups]
    )
    graphics::mtext("Dominant category by theta", side = 1, line = 4.6, adj = 0)
  }
  if (!is.null(plot_data$steps) && nrow(plot_data$steps) > 0) {
    step_groups <- as.character(plot_data$steps$CurveGroup)
    step_cols <- cols[step_groups]
    for (i in seq_len(nrow(plot_data$steps))) {
      graphics::abline(
        v = plot_data$steps$Threshold[i],
        col = grDevices::adjustcolor(step_cols[i], alpha.f = 0.16),
        lty = 3
      )
    }
  }
  for (i in seq_along(groups)) {
    sub <- exp_df[exp_df$CurveGroup == groups[i], , drop = FALSE]
    graphics::lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2)
  }
  if (nrow(plot_data$steps) > 0) {
    step_groups <- as.character(plot_data$steps$CurveGroup)
    step_cols <- cols[step_groups]
    graphics::points(plot_data$steps$Threshold, plot_data$steps$PathY, pch = 18, col = step_cols)
    graphics::text(
      x = plot_data$steps$Threshold,
      y = plot_data$steps$PathY,
      labels = truncate_axis_label(plot_data$steps$ThresholdLabel, width = 6L),
      pos = 3,
      cex = 0.62,
      col = step_cols,
      offset = 0.45
    )
  }
  if (!is.null(plot_data$endpoint_labels) && nrow(plot_data$endpoint_labels) > 0 && length(groups) <= 5) {
    graphics::text(
      x = plot_data$endpoint_labels$Theta,
      y = plot_data$endpoint_labels$ExpectedScore,
      labels = truncate_axis_label(plot_data$endpoint_labels$CurveGroup, width = 14L),
      pos = 4,
      cex = 0.78,
      xpd = NA,
      col = cols[as.character(plot_data$endpoint_labels$CurveGroup)]
    )
  } else {
    graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n")
  }
}

build_fit_pathway_data <- function(x,
                                   diagnostics = NULL,
                                   fit_stat = c("Infit", "Outfit"),
                                   fit_scale = c("mnsq", "zstd"),
                                   zstd_method = c("engine", "facets"),
                                   include_person = FALSE,
                                   person_subset = NULL,
                                   top_n_person = 30L,
                                   person_labels = c("flagged", "all", "none"),
                                   facet_labels = c("all", "flagged", "none"),
                                   panel = c("combined", "facet"),
                                   fit_range = c(0.5, 1.5),
                                   zstd_cut = 2,
                                   ci_level = 0.95) {
  fit_stat <- match.arg(fit_stat)
  fit_scale <- match.arg(tolower(as.character(fit_scale[1])), c("mnsq", "zstd"))
  zstd_method <- match.arg(tolower(as.character(zstd_method[1])), c("engine", "facets"))
  person_labels <- match.arg(tolower(as.character(person_labels[1])), c("flagged", "all", "none"))
  facet_labels <- match.arg(tolower(as.character(facet_labels[1])), c("all", "flagged", "none"))
  panel <- match.arg(tolower(as.character(panel[1])), c("combined", "facet"))
  if (!is.logical(include_person) || length(include_person) != 1L || is.na(include_person)) {
    stop("`include_person` must be TRUE or FALSE.", call. = FALSE)
  }
  fit_range <- suppressWarnings(as.numeric(fit_range))
  if (length(fit_range) != 2L || !all(is.finite(fit_range)) || fit_range[1] >= fit_range[2]) {
    stop("`fit_range` must be a finite increasing numeric length-2 vector.", call. = FALSE)
  }
  zstd_cut <- suppressWarnings(as.numeric(zstd_cut[1]))
  if (!is.finite(zstd_cut) || zstd_cut <= 0) {
    stop("`zstd_cut` must be a positive finite number.", call. = FALSE)
  }
  top_n_person <- suppressWarnings(as.numeric(top_n_person[1]))
  if (is.na(top_n_person) || top_n_person <= 0) {
    stop("`top_n_person` must be a positive number or `Inf`.", call. = FALSE)
  }

  diagnostics_source <- if (is.null(diagnostics)) "computed" else "supplied"
  diag_use <- diagnostics
  if (is.null(diag_use)) {
    diag_use <- diagnose_mfrm(
      x,
      residual_pca = "none",
      fit_df_method = if (identical(fit_scale, "zstd")) "both" else "engine"
    )
  }
  if (!is.list(diag_use) || is.null(diag_use$measures)) {
    stop("`diagnostics` must contain a `measures` table.", call. = FALSE)
  }
  if (identical(fit_scale, "zstd")) {
    selected_z_col <- paste0(fit_stat, "ZSTD_", toupper(zstd_method))
    if (!selected_z_col %in% names(diag_use$measures)) {
      diag_use <- diagnose_mfrm(x, residual_pca = "none", fit_df_method = "both")
      diagnostics_source <- paste0(diagnostics_source, "; recomputed for engine/FACETS ZSTD companions")
    }
  }

  fm <- fit_measures_table(
    x,
    diagnostics = diag_use,
    include_person = isTRUE(include_person),
    fit_df_method = if (identical(fit_scale, "zstd")) "both" else "engine",
    lower = fit_range[1],
    upper = fit_range[2],
    zstd_cut = zstd_cut,
    ci_level = ci_level,
    sort_by = "facet",
    top_n = Inf
  )
  tbl <- as.data.frame(fm$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    stop("No fit-measure rows are available for the fit pathway.", call. = FALSE)
  }

  raw_measures <- as.data.frame(diag_use$measures, stringsAsFactors = FALSE)
  metadata_cols <- intersect(
    c("SE_Method", "PrecisionTier", "SupportsFormalInference", "SEUse",
      "CIBasis", "CIUse", "Converged", "CI_Method", "CILabel"),
    names(raw_measures)
  )
  raw_key <- paste(as.character(raw_measures$Facet), as.character(raw_measures$Level), sep = "\r")
  tbl_key <- paste(as.character(tbl$Facet), as.character(tbl$Level), sep = "\r")
  meta_idx <- match(tbl_key, raw_key)
  for (nm in metadata_cols) tbl[[nm]] <- raw_measures[[nm]][meta_idx]
  for (nm in setdiff(
    c("SE_Method", "PrecisionTier", "SupportsFormalInference", "SEUse",
      "CIBasis", "CIUse", "Converged", "CI_Method", "CILabel"),
    names(tbl)
  )) {
    tbl[[nm]] <- NA
  }

  fit_col <- if (identical(fit_scale, "mnsq")) {
    fit_stat
  } else {
    candidate <- paste0(fit_stat, "ZSTD_", toupper(zstd_method))
    if (candidate %in% names(tbl)) candidate else paste0(fit_stat, "ZSTD")
  }
  if (!fit_col %in% names(tbl)) {
    stop("The selected fit statistic is not available: `", fit_col, "`.", call. = FALSE)
  }
  tbl$FitValue <- suppressWarnings(as.numeric(tbl[[fit_col]]))
  tbl$Measure <- suppressWarnings(as.numeric(tbl$Measure))
  tbl$SE <- suppressWarnings(as.numeric(tbl$SE))
  tbl$ElementType <- ifelse(as.character(tbl$Facet) == "Person", "Person", "Facet level")
  tbl$FitScale <- fit_scale
  tbl$FitStatistic <- fit_stat
  tbl$FitColumn <- fit_col
  tbl$FitDistance <- if (identical(fit_scale, "mnsq")) {
    abs(tbl$FitValue - 1)
  } else {
    abs(tbl$FitValue)
  }
  tbl$Flagged <- if (identical(fit_scale, "mnsq")) {
    tbl$FitValue < fit_range[1] | tbl$FitValue > fit_range[2]
  } else {
    abs(tbl$FitValue) >= zstd_cut
  }
  tbl$FitDirection <- ifelse(
    !is.finite(tbl$FitValue), "not_available",
    ifelse(
      identical(fit_scale, "mnsq") & tbl$FitValue > fit_range[2] |
        identical(fit_scale, "zstd") & tbl$FitValue > zstd_cut,
      "underfit",
      ifelse(
        identical(fit_scale, "mnsq") & tbl$FitValue < fit_range[1] |
          identical(fit_scale, "zstd") & tbl$FitValue < -zstd_cut,
        "overfit",
        "within_band"
      )
    )
  )
  tbl <- tbl[is.finite(tbl$FitValue) & is.finite(tbl$Measure), , drop = FALSE]
  if (nrow(tbl) == 0L) {
    stop("No finite fit/measure pairs are available for the fit pathway.", call. = FALSE)
  }

  person_rows <- tbl$ElementType == "Person"
  if (isTRUE(include_person) && !is.null(person_subset)) {
    keep_ids <- unique(as.character(person_subset))
    tbl <- tbl[!person_rows | as.character(tbl$Level) %in% keep_ids, , drop = FALSE]
    person_rows <- tbl$ElementType == "Person"
  }
  if (isTRUE(include_person) && is.finite(top_n_person) && sum(person_rows) > top_n_person) {
    person_idx <- which(person_rows)
    ord <- order(tbl$FitDistance[person_idx], decreasing = TRUE, na.last = NA)
    keep_person_idx <- person_idx[utils::head(ord, as.integer(top_n_person))]
    tbl <- tbl[!person_rows | seq_len(nrow(tbl)) %in% keep_person_idx, , drop = FALSE]
  }
  if (nrow(tbl) == 0L) {
    stop("No fit-pathway rows remain after person filtering.", call. = FALSE)
  }

  person_rows <- tbl$ElementType == "Person"
  tbl$Panel <- if (identical(panel, "facet")) as.character(tbl$Facet) else "All elements"
  tbl$Shape <- ifelse(person_rows, 15L, 21L)
  tbl$LabelText <- as.character(tbl$Level)
  if (identical(person_labels, "none")) {
    tbl$LabelText[person_rows] <- ""
  } else if (identical(person_labels, "flagged")) {
    tbl$LabelText[person_rows & !tbl$Flagged] <- ""
  }
  if (identical(facet_labels, "none")) {
    tbl$LabelText[!person_rows] <- ""
  } else if (identical(facet_labels, "flagged")) {
    tbl$LabelText[!person_rows & !tbl$Flagged] <- ""
  }

  uncertainty_cols <- c(
    "ElementType", "SE_Method", "PrecisionTier", "SupportsFormalInference",
    "SEUse", "CIBasis", "CIUse", "Converged"
  )
  uncertainty_basis <- unique(tbl[, uncertainty_cols, drop = FALSE])
  uncertainty_basis$CI_Level <- ci_level
  uncertainty_basis$WhiskerDefinition <- sprintf(
    "Measure +/- %.6g x SE (normal %.1f%% interval)",
    stats::qnorm(1 - (1 - ci_level) / 2), 100 * ci_level
  )
  rownames(uncertainty_basis) <- NULL

  reference_values <- if (identical(fit_scale, "mnsq")) {
    c(fit_range[1], 1, fit_range[2])
  } else {
    c(-zstd_cut, 0, zstd_cut)
  }
  reference_labels <- if (identical(fit_scale, "mnsq")) {
    c("Lower MnSq review line", "Expected MnSq", "Upper MnSq review line")
  } else {
    c("Lower ZSTD review line", "Expected ZSTD", "Upper ZSTD review line")
  }

  list(
    title = paste0("Fit pathway: ", fit_stat, " ", if (fit_scale == "mnsq") "MnSq" else "ZSTD", " by measure"),
    subtitle = if (identical(fit_scale, "zstd")) {
      paste0("Vertical logits with ", zstd_method, "-df ZSTD on the horizontal axis")
    } else {
      "Vertical logits with mean-square fit on the horizontal axis"
    },
    table = tbl,
    full_table = as.data.frame(fm$table, stringsAsFactors = FALSE),
    uncertainty_basis = uncertainty_basis,
    fit_stat = fit_stat,
    fit_scale = fit_scale,
    view = "fit_measure",
    fit_column = fit_col,
    zstd_method = zstd_method,
    fit_range = fit_range,
    zstd_cut = zstd_cut,
    include_person = isTRUE(include_person),
    person_subset = person_subset,
    top_n_person = top_n_person,
    person_labels = person_labels,
    facet_labels = facet_labels,
    panel = panel,
    ci_level = ci_level,
    diagnostics_source = diagnostics_source,
    legend = new_plot_legend(
      label = c("Facet level", "Person", "Measure uncertainty"),
      role = c("element", "element", "interval"),
      aesthetic = c("point", "point", "vertical whisker"),
      value = c("circle", "square", paste0(round(100 * ci_level), "% CI"))
    ),
    reference_lines = new_reference_lines(
      axis = rep("v", 3L),
      value = reference_values,
      label = reference_labels,
      linetype = c("dotted", "dashed", "dotted"),
      role = c("threshold", "reference", "threshold")
    )
  )
}

draw_fit_pathway <- function(plot_data,
                             title = NULL,
                             palette = NULL,
                             show_ci = TRUE) {
  tbl <- as.data.frame(plot_data$table, stringsAsFactors = FALSE)
  facets <- unique(as.character(tbl$Facet))
  default_cols <- stats::setNames(
    grDevices::hcl.colors(max(3L, length(facets)), "Dark 3")[seq_along(facets)],
    facets
  )
  cols <- resolve_palette(palette = palette, defaults = default_cols)
  refs <- as.data.frame(plot_data$reference_lines, stringsAsFactors = FALSE)
  x_ref <- suppressWarnings(as.numeric(refs$value))

  draw_panel <- function(panel_tbl, panel_title) {
    panel_title <- .mfrm_plot_display_title(
      plot_data,
      title = panel_title,
      fallback = "Fit pathway"
    )
    x_rng <- range(c(panel_tbl$FitValue, x_ref), finite = TRUE)
    y_values <- panel_tbl$Measure
    if (isTRUE(show_ci)) {
      y_values <- c(y_values, panel_tbl$CI_Lower, panel_tbl$CI_Upper)
    }
    y_rng <- range(y_values, finite = TRUE)
    if (diff(x_rng) <= sqrt(.Machine$double.eps)) x_rng <- x_rng[1] + c(-0.5, 0.5)
    if (diff(y_rng) <= sqrt(.Machine$double.eps)) y_rng <- y_rng[1] + c(-0.5, 0.5)
    x_rng <- x_rng + c(-0.08, 0.08) * diff(x_rng)
    y_rng <- y_rng + c(-0.08, 0.08) * diff(y_rng)
    graphics::plot(
      panel_tbl$FitValue,
      panel_tbl$Measure,
      type = "n",
      xlim = x_rng,
      ylim = y_rng,
      xlab = if (identical(plot_data$fit_scale, "mnsq")) {
        paste0(plot_data$fit_stat, " MnSq")
      } else {
        paste0(plot_data$fit_stat, " ZSTD (", plot_data$zstd_method, " df)")
      },
      ylab = "Measure (logits)",
      main = panel_title
    )
    if (identical(plot_data$fit_scale, "mnsq")) {
      graphics::rect(
        plot_data$fit_range[1], y_rng[1], plot_data$fit_range[2], y_rng[2],
        col = grDevices::adjustcolor("#E8F0FE", alpha.f = 0.55), border = NA
      )
    } else {
      graphics::rect(
        -plot_data$zstd_cut, y_rng[1], plot_data$zstd_cut, y_rng[2],
        col = grDevices::adjustcolor("#E8F0FE", alpha.f = 0.55), border = NA
      )
    }
    for (i in seq_len(nrow(refs))) {
      graphics::abline(
        v = refs$value[i],
        lty = if (refs$role[i] == "reference") 2 else 3,
        col = if (refs$role[i] == "reference") "gray35" else "gray65"
      )
    }
    graphics::abline(h = pretty(y_rng, n = 6), col = "gray92", lty = 1)
    point_cols <- unname(cols[as.character(panel_tbl$Facet)])
    if (isTRUE(show_ci)) {
      ci_ok <- is.finite(panel_tbl$CI_Lower) & is.finite(panel_tbl$CI_Upper)
      if (any(ci_ok)) {
        graphics::arrows(
          x0 = panel_tbl$FitValue[ci_ok], y0 = panel_tbl$CI_Lower[ci_ok],
          x1 = panel_tbl$FitValue[ci_ok], y1 = panel_tbl$CI_Upper[ci_ok],
          angle = 90, code = 3, length = 0.035,
          col = grDevices::adjustcolor(point_cols[ci_ok], alpha.f = 0.60)
        )
      }
    }
    graphics::points(
      panel_tbl$FitValue, panel_tbl$Measure,
      pch = panel_tbl$Shape,
      col = point_cols,
      bg = grDevices::adjustcolor(point_cols, alpha.f = 0.35),
      cex = ifelse(panel_tbl$ElementType == "Person", 0.82, 1)
    )
    label_rows <- nzchar(as.character(panel_tbl$LabelText))
    if (any(label_rows)) {
      graphics::text(
        panel_tbl$FitValue[label_rows], panel_tbl$Measure[label_rows],
        labels = truncate_axis_label(panel_tbl$LabelText[label_rows], width = 16L),
        pos = 4, cex = 0.65, offset = 0.35, xpd = NA,
        col = point_cols[label_rows]
      )
    }
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  if (identical(plot_data$panel, "facet") && length(unique(tbl$Panel)) > 1L) {
    panel_names <- unique(as.character(tbl$Panel))
    n_cols <- ceiling(sqrt(length(panel_names)))
    graphics::par(mfrow = c(ceiling(length(panel_names) / n_cols), n_cols))
    for (nm in panel_names) {
      draw_panel(tbl[as.character(tbl$Panel) == nm, , drop = FALSE], nm)
    }
  } else {
    draw_panel(tbl, title %||% plot_data$title)
    graphics::legend(
      "topleft", legend = facets, col = cols[facets], pch = 16,
      bty = "n", cex = 0.76
    )
    if (any(tbl$ElementType == "Person")) {
      graphics::legend(
        "bottomright", legend = c("Facet level", "Person"),
        pch = c(21, 15), col = "gray30", pt.bg = "gray80",
        bty = "n", cex = 0.76
      )
    }
  }
}

build_ccc_data <- function(x, theta_range = c(-6, 6), theta_points = 241L) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  prob_df <- build_curve_tables(curve_spec, theta_grid)$probabilities
  list(
    title = "Category Characteristic Curves",
    probabilities = prob_df
  )
}

build_ccc_surface_data <- function(x, theta_range = c(-6, 6), theta_points = 121L) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  prob_df <- as.data.frame(build_curve_tables(curve_spec, theta_grid)$probabilities, stringsAsFactors = FALSE)

  category_levels <- as.character(curve_spec$categories)
  prob_df$Category <- factor(as.character(prob_df$Category), levels = category_levels)
  prob_df$CategoryIndex <- as.integer(prob_df$Category)
  prob_df$CategoryScore <- suppressWarnings(as.numeric(as.character(prob_df$Category)))
  prob_df$SurfaceX <- prob_df$Theta
  prob_df$SurfaceY <- prob_df$CategoryIndex
  prob_df$SurfaceZ <- prob_df$Probability
  prob_df$Category <- as.character(prob_df$Category)
  prob_df <- prob_df[order(prob_df$CurveGroup, prob_df$CategoryIndex, prob_df$Theta), , drop = FALSE]

  data.frame_category <- data.frame(
    Category = category_levels,
    CategoryIndex = seq_along(category_levels),
    CategoryScore = suppressWarnings(as.numeric(category_levels)),
    stringsAsFactors = FALSE
  )
  data.frame_group <- data.frame(
    CurveGroup = names(curve_spec$groups),
    SurfacePanel = seq_along(curve_spec$groups),
    Model = curve_spec$model,
    stringsAsFactors = FALSE
  )
  observed_scores <- if (!is.null(x$prep$data) && "Score" %in% names(x$prep$data)) {
    as.character(x$prep$data$Score)
  } else {
    character()
  }
  observed_counts <- table(observed_scores, useNA = "no")
  category_counts <- unname(observed_counts[category_levels])
  category_counts[is.na(category_counts)] <- 0L
  category_support <- data.frame(
    Category = category_levels,
    CategoryIndex = seq_along(category_levels),
    CategoryScore = suppressWarnings(as.numeric(category_levels)),
    ObservedCount = as.integer(category_counts),
    ZeroObserved = as.integer(category_counts) == 0L,
    SupportRole = ifelse(
      as.integer(category_counts) == 0L,
      "retained score-support category with zero observed responses",
      "observed response category"
    ),
    stringsAsFactors = FALSE
  )
  axis_contract <- data.frame(
    Axis = c("x", "y", "z"),
    Column = c("SurfaceX", "SurfaceY", "SurfaceZ"),
    Label = c("Theta / Logit", "Observed category index", "Category probability"),
    stringsAsFactors = FALSE
  )
  renderer_contract <- data.frame(
    Renderer = c("base", "external 3D renderer"),
    Status = c("not rendered as 3D in base graphics", "plot data only; no plotly/rgl dependency"),
    RecommendedUse = c(
      "Use the 2D CCC/pathway views for default reports.",
      "Use the surface plot data for exploratory teaching, audit, or downstream interactive rendering."
    ),
    stringsAsFactors = FALSE
  )
  zero_count_message <- if (any(category_support$ZeroObserved)) {
    paste0(
      "Retained categories with zero observed responses: ",
      paste(category_support$Category[category_support$ZeroObserved], collapse = ", "),
      ". Keep these labels for scale support, but avoid overinterpreting adjacent threshold ridges."
    )
  } else {
    "All retained categories have at least one observed response in the fitted data."
  }
  interpretation_guide <- data.frame(
    Topic = c(
      "Axes",
      "Category dominance",
      "Zero-frequency categories",
      "Reporting use",
      "Renderer boundary"
    ),
    Guidance = c(
      "Read SurfaceX as theta/logit, SurfaceY as retained category index, and SurfaceZ as predicted category probability.",
      "A high SurfaceZ ridge marks the theta region where that category is most probable; compare with the 2D CCC/pathway view before making a reporting claim.",
      zero_count_message,
      "Use the 2D pathway or CCC figure as the default manuscript/report figure; use this surface as exploratory or teaching material.",
      "mfrmr returns the data contract only and does not execute plotly/rgl or any other interactive 3D renderer."
    ),
    stringsAsFactors = FALSE
  )
  reporting_policy <- data.frame(
    UseCase = c(
      "Manuscript core figure",
      "Appendix or teaching figure",
      "Interactive downstream rendering",
      "FACETS/ConQuest/SPSS comparison"
    ),
    Recommendation = c(
      "Prefer plot(fit, type = \"pathway\") or plot(fit, type = \"ccc\").",
      "Acceptable when explicitly described as exploratory category-probability support.",
      "Use surface, categories, groups, axis_contract, renderer_contract, and category_support as the renderer input contract.",
      "Do not present the surface as a direct FACETS, ConQuest, or SPSS-equivalent display."
    ),
    stringsAsFactors = FALSE
  )

  list(
    surface = prob_df,
    categories = data.frame_category,
    category_support = category_support,
    groups = data.frame_group,
    axis_contract = axis_contract,
    renderer_contract = renderer_contract,
    interpretation_guide = interpretation_guide,
    reporting_policy = reporting_policy,
    surface_role = "theta_by_category_by_probability",
    recommended_use = "Exploratory surface data only; keep 2D pathway/CCC views as the default reporting figures."
  )
}

draw_ccc <- function(plot_data, title = NULL, palette = NULL) {
  prob_df <- plot_data$probabilities
  traces <- unique(paste(prob_df$CurveGroup, prob_df$Category, sep = " | Cat "))
  defaults <- stats::setNames(
    grDevices::hcl.colors(max(3L, length(traces)), "Dark 3")[seq_along(traces)],
    traces
  )
  cols <- resolve_palette(palette = palette, defaults = defaults)
  graphics::plot(
    x = range(prob_df$Theta, finite = TRUE),
    y = c(0, 1),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Probability",
    main = .mfrm_plot_display_title(
      plot_data,
      title = title,
      fallback = "Category characteristic curves"
    )
  )
  for (i in seq_along(traces)) {
    parts <- strsplit(traces[i], " \\| Cat ", fixed = FALSE)[[1]]
    sub <- prob_df[prob_df$CurveGroup == parts[1] & prob_df$Category == parts[2], , drop = FALSE]
    graphics::lines(sub$Theta, sub$Probability, col = cols[traces[i]], lwd = 1.5)
  }
}

draw_person_plot <- function(person_tbl, bins, title = "Person measure distribution", palette = NULL) {
  pal <- resolve_palette(palette = palette, defaults = c(person_hist = "#2C7FB8"))
  graphics::hist(
    x = person_tbl$Estimate,
    breaks = bins,
    main = title,
    xlab = "Person estimate",
    ylab = "Count",
    col = pal["person_hist"],
    border = "white"
  )
}

draw_step_plot <- function(step_tbl,
                           title = "Step parameter estimates",
                           palette = NULL,
                           label_angle = 45) {
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      step_line = "#1B9E77",
      grid = "#ececec"
    )
  )
  x_idx <- step_index_from_label(step_tbl$Step)
  ord <- order(x_idx)
  step_tbl <- step_tbl[ord, , drop = FALSE]
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar <- old_mar
  mar[1] <- max(mar[1], 8.6)
  graphics::par(mar = mar)
  graphics::plot(
    x = seq_len(nrow(step_tbl)),
    y = step_tbl$Estimate,
    type = "b",
    pch = 16,
    xaxt = "n",
    xlab = "",
    ylab = "Estimate",
    main = title,
    col = pal["step_line"]
  )
  graphics::abline(h = pretty(step_tbl$Estimate, n = 5), col = pal["grid"], lty = 1)
  draw_rotated_x_labels(
    at = seq_len(nrow(step_tbl)),
    labels = truncate_axis_label(step_tbl$Step, width = 16L),
    srt = label_angle,
    cex = 0.84,
    line_offset = 0.085
  )
}

draw_facet_plot <- function(facet_tbl,
                            title,
                            palette = NULL,
                            show_ci = FALSE,
                            ci_level = 0.95) {
  pal <- resolve_palette(palette = palette, defaults = c(facet_bar = "gray70"))
  labels <- paste0(facet_tbl$Facet, ": ", facet_tbl$Level)
  mids <- graphics::barplot(
    height = facet_tbl$Estimate,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    main = title,
    xlab = "Estimate",
    col = pal["facet_bar"],
    border = "white"
  )
  graphics::abline(v = 0, lty = 2, col = "gray40")
  if (isTRUE(show_ci) && "SE" %in% names(facet_tbl)) {
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    ci_ok <- is.finite(facet_tbl$SE) & facet_tbl$SE > 0
    if (any(ci_ok)) {
      ci_lo <- facet_tbl$Estimate[ci_ok] - z * facet_tbl$SE[ci_ok]
      ci_hi <- facet_tbl$Estimate[ci_ok] + z * facet_tbl$SE[ci_ok]
      graphics::arrows(
        x0 = ci_lo,
        y0 = mids[ci_ok],
        x1 = ci_hi,
        y1 = mids[ci_ok],
        angle = 90,
        code = 3,
        length = 0.04,
        col = "gray30",
        lwd = 1.2
      )
    }
  }
}

#' Plot fitted MFRM results with base R
#'
#' @param x An `mfrm_fit` object from [fit_mfrm()].
#' @param type Plot type. Omit `type` (or use `NULL` / `"wright"`) for the
#'   primary native Wright map. Use `"bundle"`, `"all"`, or `"default"` for
#'   the three-part fit bundle; otherwise choose one of `"facet"`, `"person"`,
#'   `"step"`, `"wright"`, `"pathway"`, `"fit_pathway"`, `"ccc"`,
#'   `"ccc_surface"`, or `"category_surface"`. `"pathway"` is the
#'   theta-to-expected-score display; `"fit_pathway"` is the
#'   fit-statistic-to-measure display.
#' @param facet Optional facet name for `type = "facet"`.
#' @param top_n Maximum number of facet/step locations retained by the native
#'   Wright map for compact displays. Step transitions are always retained;
#'   any omitted facet locations are counted in the returned `retention` table
#'   and disclosed in the plot subtitle/note. Use `Inf` for a complete
#'   all-level final map. The FACETS-style payload always retains every fitted
#'   facet and step location.
#' @param theta_range Numeric length-2 range for pathway, CCC, and
#'   category-surface plot data.
#' @param theta_points Number of theta grid points used for pathway, CCC, and
#'   category-surface plot data.
#' @param title Optional custom title.
#' @param palette Optional color overrides.
#' @param label_angle Rotation angle for x-axis labels where applicable.
#' @param show_ci Whether to add approximate confidence intervals when
#'   available. `NULL` selects the display-specific default: `TRUE` for the
#'   primary native Wright map and `type = "fit_pathway"`, and `FALSE` for the
#'   FACETS-style ruler and other plot types. With `renderer = "facets"`,
#'   explicitly setting `show_ci = TRUE` creates a hybrid display: FACETS-style
#'   ruler grammar with mfrmr uncertainty intervals.
#' @param ci_level Confidence level used when `show_ci = TRUE`.
#' @param group Optional grouping for `type = "wright"` to overlay
#'   per-group person-density curves (DIF / DFF screening view).
#'   Either a column name (looked up first in `group_data` when
#'   supplied through `...`, then in `fit$prep$data`) or a vector
#'   aligned with `fit$facets$person`. Ignored for other `type`
#'   values and for `wright_style = "facets_style"`, whose person column is a
#'   single FACETS-style star frequency. To pass the source data alongside, use
#'   `plot(fit, type = "wright", group = "MyCol", group_data = <df>)`.
#' @param diagnostics Optional output from [diagnose_mfrm()]. When supplied,
#'   Wright-map standard errors and precision metadata reuse matching rows from
#'   `diagnostics$measures` without replacing fitted coordinates,
#'   while pathway plot data reuse `fit_measures`, `fit_status`, and
#'   `curve_fit_status` instead of recomputing diagnostics.
#' @param include_fit_measures If `TRUE` (default), pathway plot data include
#'   tidy fit-measure and fit-status tables for custom R graphics. Set to
#'   `FALSE` when only the curve coordinates are needed.
#' @param fit_stat For `type = "fit_pathway"`, the horizontal fit statistic:
#'   `"Infit"` (default) or `"Outfit"`.
#' @param fit_scale For `type = "fit_pathway"`, use mean squares (`"mnsq"`)
#'   or standardized fit (`"zstd"`) on the horizontal axis.
#' @param zstd_method For a ZSTD fit pathway, use the package engine degrees
#'   of freedom (`"engine"`, default) or the FACETS/Wright-Masters companion
#'   calculation (`"facets"`). The latter is a comparison aid, not a claim of
#'   complete FACETS output equivalence.
#' @param include_person If `TRUE`, include person rows in a fit pathway.
#' @param person_subset Optional character vector of person IDs retained in a
#'   fit pathway before ranking.
#' @param top_n_person Maximum person rows retained in a fit pathway, ranked
#'   by distance from expected fit. Use `Inf` for all selected persons.
#' @param person_labels Person-label policy for a fit pathway: `"flagged"`
#'   (default), `"all"`, or `"none"`.
#' @param facet_labels Facet-level label policy for a fit pathway: `"all"`
#'   (default), `"flagged"`, or `"none"`. This changes drawn labels only;
#'   all retained facet rows remain in the draw-free table.
#' @param panel Fit-pathway layout: `"combined"` (distinguish persons by
#'   shape) or `"facet"` (one base-graphics panel per facet).
#' @param fit_range Mean-square review lines for a fit pathway.
#' @param zstd_cut Absolute ZSTD review line for a fit pathway.
#' @param draw If `TRUE`, draw the plot with base graphics.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`,
#'   or `"monochrome"`).
#' @param renderer Wright-map renderer. Use `"native"` (default) or
#'   `"facets"` for the FACETS Table 6-style visual layout. This is the
#'   canonical selector for new code. The native map is the primary display
#'   and includes available facet uncertainty by default. The closest
#'   FACETS-style visual uses `show_ci = FALSE`.
#' @param wright_style Wright-map renderer. `"native"` preserves the package's
#'   histogram, point, range, and facet-SE display. `"facets_style"` adds a
#'   FACETS Table 6-style text ruler with person-frequency stars, signed facet
#'   headers, and labeled score-transition lines. It is a visual-layout option,
#'   not a claim of numerical equivalence with FACETS. This explicit style name
#'   is equivalent to `renderer = "facets"`. Adding `show_ci = TRUE` to that
#'   style is an mfrmr/FACETS hybrid rather than the closest FACETS-style view.
#' @param category_labels Optional score rubric labels used by
#'   `wright_style = "facets_style"`. Supply a named character vector keyed by
#'   original score, an unnamed vector with one label per retained category, or
#'   a data frame with `Score` and `Label` columns.
#' @param rows_per_logit Number of FACETS-style ruler rows per logit.
#' @param wright_range Optional finite increasing length-2 logit range for the
#'   FACETS-style ruler. `NULL` derives a range that contains the fitted data.
#' @param extreme_placement In the FACETS-style renderer, place extreme-score
#'   persons at the ruler `"ends"` or retain their fitted `"estimate"`.
#' @param persons_per_star Number of persons represented by one `*` in the
#'   FACETS-style frequency column. `NULL` chooses a compact value automatically.
#' @param ... Additional arguments ignored for S3 compatibility.
#'
#' @details
#' This S3 plotting method provides the core fit-family visuals for
#' `mfrmr`. When `type` is omitted, it returns the Wright map alone as
#' an `mfrm_plot_data` object (the most useful single figure for a
#' first inspection). Pass `type = "bundle"` (or `"all"` / `"default"`)
#' to obtain the legacy three-plot `mfrm_plot_bundle` containing a
#' Wright map, pathway map, and category characteristic curves. The
#' compact native default records any omitted facet locations in
#' `data$retention` and annotates the subtitle and drawn figure; use
#' `top_n = Inf` to retain every fitted location in the final Wright map.
#' Native text labels remain collision-aware, so a retained point may be
#' unlabeled; its exact level and estimate remain available in the returned
#' plot data. When the fit records boundary-separated facet levels and no
#' display range was supplied, the native and FACETS-style maps use the same
#' robust automatic range and place those levels at ruler ends. Exact fitted
#' values and intervals remain in `OriginalEstimate`, `CI_Lower`, and
#' `CI_Upper`; `DisplayEstimate`, `DisplayCI_Lower`, `DisplayCI_Upper`, and the
#' clipping-status columns describe only the rendered coordinates. Endpoint
#' triangles and the plot footer disclose any omitted or clipped interval.
#' The returned object always carries machine-readable metadata through
#' the `mfrm_plot_data` contract, even when the plot is drawn
#' immediately.
#'
#' Every fit-derived payload also carries `data$fit_readiness`,
#' `data$interpretation_status`, and `data$interpretation_note`. Availability
#' and interpretability are separate: when a numerical, data, design, or
#' stability gate requires review, the coordinates remain available for
#' diagnosis, but the call warns and marks the returned subtitle and drawn
#' title `REVIEW ONLY`.
#'
#' `type = "wright"` shows persons, facet levels, and step thresholds on
#' a shared logit scale. Estimates are plotted as fitted, so the sign
#' convention follows the fit: higher person values indicate higher
#' ability, and higher non-person facet values indicate greater
#' severity/difficulty under the default negative facet orientation.
#' Facets listed in `fit_mfrm(positive_facets = ...)` are reversed
#' (higher values raise expected scores); state the active orientation in
#' figure captions when reporting.
#'
#' Set `renderer = "facets"` (equivalently,
#' `wright_style = "facets_style"`) for a line-printer-inspired common
#' ruler. That mode retains every facet location in its data and shows actual
#' (original, when scores were internally recoded) adjacent score transitions.
#' Its `facets_style` payload contains tidy ruler, person-frequency, facet,
#' step, mean-half-score, header, category-label, and settings tables for custom
#' graphics. `RulerValue` records the nearest discrete ruler row for
#' line-printer reconstruction, whereas step and midpoint lines are drawn at
#' their exact fitted `DrawValue` and step labels print that logit value. The
#' styling emulates the semantics of FACETS Table 6; estimates
#' and standard errors remain those produced by `mfrmr`. Use
#' `show_ci = FALSE` for the closest FACETS-style rendering. If
#' `show_ci = TRUE` is requested, interpret the result as a hybrid that adds
#' mfrmr uncertainty intervals to FACETS-style ruler grammar.
#'
#' `type = "pathway"` shows expected score
#' traces and dominant-category regions across theta. This expected-score
#' display is distinct from the Bond-and-Fox-style fit-versus-measure
#' pathway: use `type = "fit_pathway"` for Infit/Outfit on the horizontal
#' axis and measure logits on the vertical axis. Its draw-free plot data
#' include person-selection settings, CI columns, and explicit SE-source
#' metadata. The expected-score pathway draw-free data also
#' includes `pathway_long`, `pathway_annotations`, `fit_measures`,
#' `fit_status`, and `curve_fit_status`, so R users can rebuild the pathway
#' map in ggplot2, plotly, or a report pipeline while keeping the same
#' underfit/overfit labels used by [fit_measures_table()]. `type = "ccc"` shows
#' category response probabilities. `type = "ccc_surface"` or
#' `type = "category_surface"` returns 3D-ready category-probability surface
#' data for external rendering; it deliberately does not add a
#' plotly/rgl dependency or replace the 2D CCC/pathway reporting figures. The
#' returned object includes `category_support`, `interpretation_guide`, and
#' `reporting_policy` tables so retained zero-frequency categories and
#' manuscript-use boundaries remain visible to beginners. The remaining
#' types (`"facet"`, `"person"`, `"step"`, `"shrinkage"`) provide
#' compact location-specific displays.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Use `plot(fit)` to inspect the Wright map at a glance.
#' 3. Switch to `type = "pathway"`, `"fit_pathway"`, `"ccc"`, or
#'    `"shrinkage"` for the
#'    relevant follow-up figure, or `type = "bundle"` for the
#'    three-plot overview when preparing a FACETS-style summary.
#'
#' @section Further guidance:
#' For a plot-selection guide and extended examples, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @return Invisibly, an `mfrm_plot_data` object (default and for any
#'   single `type`), or an `mfrm_plot_bundle` when
#'   `type = "bundle"` / `"all"` / `"default"`. Each returned fit plot
#'   includes domain-specific readiness and an interpretation status in its
#'   data payload.
#' @seealso [fit_mfrm()], [plot_wright_unified()], [plot_bubble()],
#'   [mfrmr_visual_diagnostics]
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept shrinkage
#' @examples
#' toy <- load_mfrmr_data("example_operational")
#' # Seven quadrature points keep this executable example short. For a final
#' # analysis, restore the default or a prespecified grid and review sensitivity.
#' fit <- fit_mfrm(
#'   toy,
#'   "Person",
#'   c("Rater", "Criterion"),
#'   "Score",
#'   method = "MML",
#'   model = "RSM",
#'   quad_points = 7,
#'   maxit = 30
#' )
#' wright <- plot(fit, draw = FALSE)
#' head(wright$data$locations)
#' # Look for: persons clustered against the facet / step rows on the
#' #   shared logit axis. Large gaps between the person density and
#' #   the step / facet rails indicate weak targeting; ceiling /
#' #   floor stripes mean the test is too easy / hard.
#' bundle <- plot(fit, type = "bundle", draw = FALSE)
#' bundle$wright_map$data$group_summary
#' # Look for: pathway curves rising in the expected order with
#' #   visible dominant-category bands; CCC curves peaking sequentially
#' #   without one category being completely overlapped by neighbours.
#' surface <- plot(fit, type = "ccc_surface", draw = FALSE)
#' head(surface$data$surface)
#' surface$data$category_support
#' # Look for: every retained category having `Observed > 0`; categories
#' #   with zero observations are returned as a zero-observation slice and
#' #   should not be interpreted as a real score region.
#' surface$data$interpretation_guide
#' if (interactive()) {
#'   plot(
#'     fit,
#'     type = "wright",
#'     preset = "publication",
#'     title = "Customized Wright Map",
#'     show_ci = TRUE,
#'     label_angle = 45
#'   )
#'   plot(
#'     fit,
#'     type = "pathway",
#'     title = "Customized Pathway Map",
#'     palette = c("#1f78b4")
#'   )
#'   plot(
#'     fit,
#'     type = "ccc",
#'     title = "Customized Category Characteristic Curves",
#'     palette = c("#1b9e77", "#d95f02", "#7570b3")
#'   )
#' }
#' @export
plot.mfrm_fit <- function(x,
                          type = NULL,
                          facet = NULL,
                          top_n = 30,
                          theta_range = c(-6, 6),
                          theta_points = 241,
                          title = NULL,
                          palette = NULL,
                          label_angle = 45,
                          show_ci = NULL,
                          ci_level = 0.95,
                          group = NULL,
                          diagnostics = NULL,
                          include_fit_measures = TRUE,
                          fit_stat = c("Infit", "Outfit"),
                          fit_scale = c("mnsq", "zstd"),
                          zstd_method = c("engine", "facets"),
                          include_person = FALSE,
                          person_subset = NULL,
                          top_n_person = 30,
                          person_labels = c("flagged", "all", "none"),
                          facet_labels = c("all", "flagged", "none"),
                          panel = c("combined", "facet"),
                          fit_range = c(0.5, 1.5),
                          zstd_cut = 2,
                          draw = TRUE,
                          preset = c("standard", "publication", "compact", "monochrome"),
                          renderer = NULL,
                          wright_style = c("native", "facets_style"),
                          category_labels = NULL,
                          rows_per_logit = 2L,
                          wright_range = NULL,
                          extreme_placement = c("ends", "estimate"),
                          persons_per_star = NULL,
                          ...) {
  if (!inherits(x, "mfrm_fit")) {
    stop("`x` must be an mfrm_fit object from fit_mfrm().")
  }
  show_ci_auto <- missing(show_ci) || is.null(show_ci)
  top_n <- suppressWarnings(as.numeric(top_n[1]))
  if (length(top_n) != 1L || is.na(top_n) || top_n <= 0) {
    stop("`top_n` must be a positive number or `Inf`.", call. = FALSE)
  }
  if (is.finite(top_n)) {
    top_n <- as.integer(floor(top_n))
    if (top_n < 1L) stop("`top_n` must be at least 1 or `Inf`.", call. = FALSE)
  }
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  style <- resolve_plot_preset(preset)
  wright_style <- match_wright_style(wright_style, renderer = renderer)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }

  # Resolve the primary route before deciding the uncertainty default. Native
  # Wright maps retain mfrmr's uncertainty advantage; the FACETS-style ruler
  # stays closest to its source grammar unless a hybrid is requested explicitly.
  if (missing(type) || is.null(type)) {
    type <- "wright"
  }
  requested_type <- tolower(as.character(type[1]))
  bundle_requested <- requested_type %in% c("bundle", "all", "default")
  if (isTRUE(show_ci_auto)) {
    show_ci <- isTRUE(
      (requested_type %in% "wright" || bundle_requested) &&
        identical(wright_style, "native")
    ) || requested_type %in% "fit_pathway"
  } else if (!is.logical(show_ci) || length(show_ci) != 1L || is.na(show_ci)) {
    stop("`show_ci` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }

  fit_plot_readiness <- .mfrm_fit_plot_readiness(x)
  readiness_warning_emitted <- FALSE

  as_plot_data <- function(name, data) {
    data_names <- names(data) %||% rep("", length(data))
    duplicated_named <- nzchar(data_names) & duplicated(data_names, fromLast = TRUE)
    out <- new_mfrm_plot_data(name, data[!duplicated_named])
    out <- .mfrm_attach_plot_readiness(out, fit_plot_readiness)
    if (!isTRUE(fit_plot_readiness$ready) && !isTRUE(readiness_warning_emitted)) {
      warning(fit_plot_readiness$detail, call. = FALSE)
      readiness_warning_emitted <<- TRUE
    }
    out
  }
  append_retention_subtitle <- function(plot_object) {
    note <- as.character(plot_object$data$retention_note %||% "")
    if (length(note) > 0L && nzchar(note[1])) {
      plot_object$data$subtitle <- paste(plot_object$data$subtitle, note[1], sep = " - ")
    }
    plot_object
  }
  wright_legend <- function(plot_core) {
    .wright_map_legend(
      plot_core,
      style = style,
      show_ci = show_ci,
      ci_level = ci_level
    )
  }

  se_tbl_ci <- if (isTRUE(show_ci) || !is.null(diagnostics)) {
    compute_se_for_plot(x, ci_level = ci_level, diagnostics = diagnostics)
  } else {
    NULL
  }

  # A bare plot(fit) call renders the primary native Wright map. The previous
  # three-plot bundle (Wright + pathway + CCC) remains explicit.
  if (isTRUE(bundle_requested)) {
    wright_core <- build_wright_map_data(
      x,
      top_n = top_n,
      se_tbl = se_tbl_ci,
      wright_style = wright_style,
      category_labels = category_labels,
      rows_per_logit = rows_per_logit,
      wright_range = wright_range,
      extreme_placement = extreme_placement,
      persons_per_star = persons_per_star
    )
    out <- list(
      wright_map = as_plot_data("wright_map", c(
        wright_core,
        list(
          title = title %||% "Wright map",
          subtitle = if (identical(wright_style, "facets_style")) {
            if (isTRUE(show_ci)) {
              "Hybrid: FACETS Table 6-style ruler with mfrmr uncertainty intervals"
            } else {
              "FACETS Table 6-style visual layout (mfrmr estimates)"
            }
          } else {
            "Shared logit scale for persons, facets, and thresholds"
          },
          show_ci = isTRUE(show_ci),
          uncertainty_display = if (identical(wright_style, "facets_style")) {
            if (isTRUE(show_ci)) "hybrid_mfrmr_ci" else "facets_style_no_ci"
          } else if (isTRUE(show_ci)) {
            "native_mfrmr_ci"
          } else {
            "native_no_ci"
          },
          preset = style$name,
          legend = wright_legend(wright_core),
          reference_lines = new_reference_lines("h", 0, "Centered logit reference", "dashed", "reference")
        )
      )),
      pathway_map = as_plot_data("pathway_map", c(
        build_pathway_map_data(
          x,
          theta_range = theta_range,
          theta_points = theta_points,
          diagnostics = diagnostics,
          include_fit_measures = include_fit_measures
        ),
        list(
          title = title %||% "Expected-score pathway (theta to score)",
          subtitle = "Expected score and dominant categories across the latent continuum",
          preset = style$name,
          legend = new_plot_legend(
            label = c("Step thresholds", "Dominant-category regions"),
            role = c("threshold", "region"),
            aesthetic = c("line", "band"),
            value = c(style$accent_secondary, style$fill_soft)
          ),
          reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
        )
      )),
      category_characteristic_curves = as_plot_data("category_characteristic_curves", c(
        build_ccc_data(x, theta_range = theta_range, theta_points = theta_points),
        list(
          title = title %||% "Category characteristic curves",
          subtitle = "Category response probabilities across theta",
          preset = style$name,
          legend = new_plot_legend(
            label = "Category curves",
            role = "probability",
            aesthetic = "line",
            value = "category_palette"
          ),
          reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
        )
      ))
    )
    out$wright_map <- append_retention_subtitle(out$wright_map)
    class(out) <- c("mfrm_plot_bundle", class(out))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_wright_map(
        out$wright_map$data,
        title = title,
        palette = resolve_palette(
          palette = palette,
          defaults = c(
            facet_level = style$accent_tertiary,
            step_threshold = style$accent_secondary,
            person_hist = style$fill_muted,
            grid = style$grid
          )
        ),
        label_angle = label_angle,
        show_ci = show_ci,
        ci_level = ci_level
      )
      draw_pathway_map(
        out$pathway_map$data,
        title = title,
        palette = palette %||% c(style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn)
      )
      draw_ccc(
        out$category_characteristic_curves$data,
        title = title,
        palette = palette %||% c(style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn)
      )
    }
    return(invisible(out))
  }

  # Locale-independent validation (match.arg() would produce a
  # translated error string on non-English locales).
  type_choices <- c("facet", "person", "step", "wright", "pathway", "fit_pathway",
                    "ccc", "ccc_overlay", "ccc_surface", "category_surface",
                    "shrinkage")
  type_in <- tolower(as.character(type[1]))
  if (!(type_in %in% type_choices)) {
    stop(sprintf(
      "`type = \"%s\"` is not a recognized plot type. Choose one of: %s.",
      type_in, paste(shQuote(type_choices), collapse = ", ")
    ), call. = FALSE)
  }
  type <- type_in
  if (identical(type, "category_surface")) type <- "ccc_surface"

  if (type == "shrinkage") {
    data_list <- .build_shrinkage_plot_data(x)
    if (isTRUE(show_ci)) {
      ci_level <- .validate_shrinkage_ci_level(ci_level)
      data_list$table <- .add_shrinkage_ci_columns(
        data_list$table,
        ci_level = ci_level
      )
    }
    legend_label <- c("Original estimate", "Shrunk estimate",
                      "Shrinkage direction", "Sum-to-zero (reference)")
    legend_role <- c("location", "location", "arrow", "reference")
    legend_aesthetic <- c("point", "point", "arrow", "line")
    legend_value <- c(style$accent_primary, style$accent_tertiary,
                      style$neutral, style$neutral)
    if (isTRUE(show_ci)) {
      legend_label <- c(legend_label, sprintf("%g%% CI whisker",
                                              round(100 * ci_level)))
      legend_role <- c(legend_role, "interval")
      legend_aesthetic <- c(legend_aesthetic, "line")
      legend_value <- c(legend_value, style$accent_primary)
    }
    out <- as_plot_data("shrinkage", list(
      data = data_list$table,
      shrinkage_report = data_list$report,
      mode = data_list$mode,
      show_ci = isTRUE(show_ci),
      ci_level = if (isTRUE(show_ci)) ci_level else NA_real_,
      title = title %||% "Empirical-Bayes shrinkage",
      subtitle = sprintf(
        "Original (filled) vs shrunk (open) estimates; mode = %s",
        data_list$mode
      ),
      preset = style$name,
      legend = new_plot_legend(
        label = legend_label,
        role = legend_role,
        aesthetic = legend_aesthetic,
        value = legend_value
      ),
      reference_lines = new_reference_lines(
        "v", 0, "Sum-to-zero reference", "dashed", "reference"
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      .draw_shrinkage_plot(data_list, style = style,
                           title = .mfrm_plot_display_title(
                             out$data,
                             title = title,
                             fallback = "Empirical-Bayes shrinkage"
                           ),
                           show_ci = isTRUE(show_ci),
                           ci_level = ci_level)
    }
    return(invisible(out))
  }

  if (type == "wright") {
    # Build optional subgroup density data when `group` names a
    # column on the prepared data. Each subgroup gets its own per-bin
    # density profile so DIF screening reads the same axis as the
    # main Wright map.
    group_payload <- NULL
    if (!is.null(group)) {
      if (identical(wright_style, "facets_style")) {
        warning(
          "`group` overlays are available only for `wright_style = \"native\"`; ",
          "the FACETS-style star column remains ungrouped.",
          call. = FALSE
        )
      } else {
        group_payload <- tryCatch(
          .build_wright_group_density(x, group = group,
                                       group_data = list(...)$group_data),
          error = function(e) {
            message(sprintf(
              "Wright group overlay disabled: %s",
              conditionMessage(e)
            ))
            NULL
          }
        )
      }
    }
    wright_core <- build_wright_map_data(
      x,
      top_n = top_n,
      se_tbl = se_tbl_ci,
      wright_style = wright_style,
      category_labels = category_labels,
      rows_per_logit = rows_per_logit,
      wright_range = wright_range,
      extreme_placement = extreme_placement,
      persons_per_star = persons_per_star
    )
    out <- as_plot_data("wright_map", c(
      wright_core,
      list(
        title = title %||% "Wright map",
        subtitle = if (identical(wright_style, "facets_style")) {
          if (isTRUE(show_ci)) {
            "Hybrid: FACETS Table 6-style ruler with mfrmr uncertainty intervals"
          } else {
            "FACETS Table 6-style visual layout (mfrmr estimates; not numerical equivalence)"
          }
        } else if (is.null(group_payload)) {
          "Shared logit scale for persons, facets, and thresholds"
        } else {
          sprintf(
            "Shared logit scale; person density split by `%s` (%d level(s))",
            group, length(unique(group_payload$Group))
          )
        },
        show_ci = isTRUE(show_ci),
        uncertainty_display = if (identical(wright_style, "facets_style")) {
          if (isTRUE(show_ci)) "hybrid_mfrmr_ci" else "facets_style_no_ci"
        } else if (isTRUE(show_ci)) {
          "native_mfrmr_ci"
        } else {
          "native_no_ci"
        },
        group = group_payload,
        preset = style$name,
        legend = wright_legend(wright_core),
        reference_lines = new_reference_lines("h", 0, "Centered logit reference", "dashed", "reference")
      )
    ))
    out <- append_retention_subtitle(out)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_wright_map(
        out$data,
        title = title,
        palette = resolve_palette(
          palette = palette,
          defaults = c(
            facet_level = style$accent_tertiary,
            step_threshold = style$accent_secondary,
            person_hist = style$fill_muted,
            grid = style$grid
          )
        ),
        label_angle = label_angle,
        show_ci = show_ci,
        ci_level = ci_level
      )
      if (identical(wright_style, "native") && !is.null(group_payload) && nrow(group_payload) > 0L) {
        # Overlay per-group density curves on the Wright map's left
        # density column, using the existing usr coordinate system.
        usr <- graphics::par("usr")
        groups <- unique(group_payload$Group)
        cols <- grDevices::hcl.colors(length(groups), "Dark 3")
        for (i in seq_along(groups)) {
          sub <- group_payload[group_payload$Group == groups[i], ,
                                drop = FALSE]
          if (nrow(sub) < 2L) next
          # Scale density to the left margin band (10% of the x-axis).
          x_band <- usr[1] + diff(usr[1:2]) * 0.10 * sub$Density /
            max(group_payload$Density, na.rm = TRUE)
          graphics::lines(x_band, sub$Theta, col = cols[i], lwd = 2)
        }
        graphics::legend(
          "topright", legend = groups, col = cols, lwd = 2,
          bty = "n", cex = 0.8, inset = 0.02
        )
      }
    }
    return(invisible(out))
  }
  if (type == "fit_pathway") {
    fit_pathway_show_ci <- isTRUE(show_ci)
    out <- as_plot_data("fit_pathway", c(
      build_fit_pathway_data(
        x,
        diagnostics = diagnostics,
        fit_stat = fit_stat,
        fit_scale = fit_scale,
        zstd_method = zstd_method,
        include_person = include_person,
        person_subset = person_subset,
        top_n_person = top_n_person,
        person_labels = person_labels,
        facet_labels = facet_labels,
        panel = panel,
        fit_range = fit_range,
        zstd_cut = zstd_cut,
        ci_level = ci_level
      ),
      list(
        show_ci = fit_pathway_show_ci,
        preset = style$name
      )
    ))
    if (!is.null(title)) out$data$title <- as.character(title[1])
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_fit_pathway(
        out$data,
        title = title,
        palette = palette,
        show_ci = fit_pathway_show_ci
      )
    }
    return(invisible(out))
  }
  if (type == "pathway") {
    out <- as_plot_data("pathway_map", c(
      build_pathway_map_data(
        x,
        theta_range = theta_range,
        theta_points = theta_points,
        diagnostics = diagnostics,
        include_fit_measures = include_fit_measures
      ),
      list(
        title = title %||% "Expected-score pathway (theta to score)",
        subtitle = "Expected score and dominant categories across the latent continuum",
        preset = style$name,
        legend = new_plot_legend(
          label = c("Step thresholds", "Dominant-category regions"),
          role = c("threshold", "region"),
          aesthetic = c("line", "band"),
          value = c(style$accent_secondary, style$fill_soft)
        ),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_pathway_map(out$data, title = title, palette = palette %||% c(
        style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn
      ))
    }
    return(invisible(out))
  }
  if (type == "ccc") {
    out <- as_plot_data("category_characteristic_curves", c(
      build_ccc_data(x, theta_range = theta_range, theta_points = theta_points),
      list(
        title = title %||% "Category characteristic curves",
        subtitle = "Category response probabilities across theta",
        preset = style$name,
        legend = new_plot_legend(
          label = "Category curves",
          role = "probability",
          aesthetic = "line",
          value = "category_palette"
        ),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_ccc(out$data, title = title, palette = palette %||% c(
        style$accent_primary, style$accent_secondary, style$accent_tertiary, style$warn
      ))
    }
    return(invisible(out))
  }
  if (type == "ccc_overlay") {
    # Empirical category proportions binned by person measure, overlaid
    # on the model CCC curves. Provides a visual model-data fit check.
    base <- as_plot_data("category_characteristic_curves", c(
      build_ccc_data(x, theta_range = theta_range,
                     theta_points = theta_points),
      list(preset = style$name)
    ))
    person_tbl <- as.data.frame(x$facets$person, stringsAsFactors = FALSE)
    obs_df <- as.data.frame(x$prep$data %||% NULL, stringsAsFactors = FALSE)
    overlay <- tryCatch({
      if (nrow(person_tbl) == 0L || nrow(obs_df) == 0L) {
        data.frame(Bin = integer(0), Theta = numeric(0),
                   Category = integer(0), Proportion = numeric(0),
                   N = integer(0), stringsAsFactors = FALSE)
      } else {
        merged <- merge(
          obs_df[, c("Person", "Score")],
          person_tbl[, c("Person", "Estimate")],
          by = "Person", all.x = TRUE
        )
        merged <- merged[is.finite(merged$Estimate) &
                           is.finite(merged$Score), , drop = FALSE]
        if (nrow(merged) == 0L) {
          data.frame(Bin = integer(0), Theta = numeric(0),
                     Category = integer(0), Proportion = numeric(0),
                     N = integer(0), stringsAsFactors = FALSE)
        } else {
          n_bins <- max(4L, min(12L, floor(nrow(merged) / 30)))
          breaks <- stats::quantile(
            merged$Estimate,
            probs = seq(0, 1, length.out = n_bins + 1L),
            na.rm = TRUE, names = FALSE
          )
          breaks <- unique(breaks)
          if (length(breaks) < 2L) {
            breaks <- range(merged$Estimate) + c(-0.05, 0.05)
          }
          merged$Bin <- as.integer(cut(merged$Estimate, breaks = breaks,
                                       include.lowest = TRUE))
          cats <- sort(unique(merged$Score))
          rows <- list()
          for (b in sort(unique(merged$Bin))) {
            sub <- merged[merged$Bin == b, , drop = FALSE]
            if (nrow(sub) == 0L) next
            theta_b <- mean(sub$Estimate, na.rm = TRUE)
            n_b <- nrow(sub)
            for (k in cats) {
              p_k <- mean(sub$Score == k, na.rm = TRUE)
              rows[[length(rows) + 1L]] <- data.frame(
                Bin = b, Theta = theta_b, Category = k,
                Proportion = p_k, N = n_b,
                stringsAsFactors = FALSE
              )
            }
          }
          do.call(rbind, rows) %||% data.frame()
        }
      }
    }, error = function(e) {
      data.frame(Bin = integer(0), Theta = numeric(0),
                 Category = integer(0), Proportion = numeric(0),
                 N = integer(0), stringsAsFactors = FALSE)
    })
    ccc_payload <- build_ccc_data(x, theta_range = theta_range,
                                   theta_points = theta_points)
    out <- as_plot_data("category_characteristic_curves_overlay", c(
      ccc_payload,
      list(
        overlay = overlay,
        title = title %||% "Category curves with empirical overlay",
        subtitle = sprintf(
          "Model curves vs observed proportions in %d theta bin(s)",
          length(unique(overlay$Bin))
        ),
        preset = style$name,
        legend = new_plot_legend(
          label = c("Model category curves", "Observed bin proportion"),
          role = c("probability", "empirical"),
          aesthetic = c("line", "point"),
          value = c(style$accent_primary, style$fail)
        ),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference",
                                               "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_ccc(out$data, title = title,
               palette = palette %||% c(
                 style$accent_primary, style$accent_secondary,
                 style$accent_tertiary, style$warn
               ))
      if (nrow(overlay) > 0L) {
        cats <- sort(unique(overlay$Category))
        cat_palette <- grDevices::hcl.colors(length(cats), "Dark 3")
        max_n <- max(overlay$N, na.rm = TRUE)
        for (i in seq_along(cats)) {
          sub <- overlay[overlay$Category == cats[i], , drop = FALSE]
          graphics::points(sub$Theta, sub$Proportion, pch = 21,
                           bg = cat_palette[i], col = "white",
                           cex = 0.8 + 0.6 * sqrt(sub$N / max(1, max_n)))
        }
      }
    }
    return(invisible(out))
  }
  if (type == "ccc_surface") {
    out <- as_plot_data("category_probability_surface", c(
      build_ccc_surface_data(x, theta_range = theta_range, theta_points = theta_points),
      list(
        title = title %||% "Category probability surface",
        subtitle = "3D-ready theta x category x probability data; no package-native 3D renderer",
        preset = style$name,
        legend = new_plot_legend(
          label = "Category probability surface",
          role = "probability_surface",
          aesthetic = "surface",
          value = "SurfaceZ"
        ),
        reference_lines = new_reference_lines("x", 0, "Centered theta reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      warning(
        "`type = \"ccc_surface\"` currently returns 3D-ready `mfrm_plot_data`; ",
        "use `draw = FALSE` and pass `x$data$surface` to an external renderer if needed.",
        call. = FALSE
      )
    }
    return(invisible(out))
  }
  if (type == "person") {
    person_tbl <- tibble::as_tibble(x$facets$person)
    person_tbl <- person_tbl[is.finite(person_tbl$Estimate), , drop = FALSE]
    if (nrow(person_tbl) == 0) stop("No finite person estimates available for plotting.")
    bins <- max(10L, min(35L, as.integer(round(sqrt(nrow(person_tbl))))))
    this_title <- if (is.null(title)) "Person measure distribution" else as.character(title[1])
    out <- as_plot_data("person", list(
      person = person_tbl,
      bins = bins,
      title = this_title,
      subtitle = "Distribution of person measures on the logit scale",
      legend = new_plot_legend("Person distribution", "distribution", "histogram", style$accent_primary),
      reference_lines = new_reference_lines(),
      preset = style$name
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_person_plot(
        person_tbl,
        bins = bins,
        title = .mfrm_plot_display_title(
          out$data,
          title = this_title,
          fallback = "Person measure distribution"
        ),
        palette = resolve_palette(palette = palette, defaults = c(person_hist = style$accent_primary))
      )
    }
    return(invisible(out))
  }
  if (type == "step") {
    step_tbl <- tibble::as_tibble(x$steps)
    if (nrow(step_tbl) == 0 || !all(c("Step", "Estimate") %in% names(step_tbl))) {
      stop("Step estimates are not available in this fit object.")
    }
    this_title <- if (is.null(title)) "Step parameter estimates" else as.character(title[1])
    out <- as_plot_data("step", list(
      steps = step_tbl,
      title = this_title,
      subtitle = "Estimated step thresholds on the shared logit scale",
      legend = new_plot_legend("Step thresholds", "threshold", "line-point", style$accent_tertiary),
      reference_lines = new_reference_lines("h", 0, "Centered step reference", "dashed", "reference"),
      preset = style$name
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_step_plot(
        step_tbl,
        title = .mfrm_plot_display_title(
          out$data,
          title = this_title,
          fallback = "Step parameter estimates"
        ),
        palette = resolve_palette(
          palette = palette,
          defaults = c(step_line = style$accent_tertiary, grid = style$grid)
        ),
        label_angle = label_angle
      )
    }
    return(invisible(out))
  }

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (nrow(facet_tbl) == 0 || !all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available in this fit object.")
  }
  if (!is.null(facet)) {
    facet_tbl <- dplyr::filter(facet_tbl, .data$Facet == as.character(facet[1]))
    if (nrow(facet_tbl) == 0) stop("Requested `facet` not found in facet-level estimates.")
  }
  if (isTRUE(show_ci) && !is.null(se_tbl_ci) && is.data.frame(se_tbl_ci) &&
      nrow(se_tbl_ci) > 0 && all(c("Facet", "Level", "SE") %in% names(se_tbl_ci))) {
    se_join <- se_tbl_ci[, intersect(c("Facet", "Level", "SE"), names(se_tbl_ci)), drop = FALSE]
    se_join$Level <- as.character(se_join$Level)
    facet_tbl$Level <- as.character(facet_tbl$Level)
    facet_tbl <- merge(facet_tbl, se_join, by = c("Facet", "Level"), all.x = TRUE, sort = FALSE)
    facet_tbl <- tibble::as_tibble(facet_tbl)
  }
  facet_tbl <- dplyr::arrange(facet_tbl, .data$Estimate)
  if (nrow(facet_tbl) > top_n) {
    facet_tbl <- facet_tbl |>
      dplyr::mutate(AbsEstimate = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$AbsEstimate)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::arrange(.data$Estimate) |>
      dplyr::select(-"AbsEstimate")
  }
  facet_title <- if (is.null(facet)) "Facet-level estimates" else paste0("Facet-level estimates: ", as.character(facet[1]))
  if (!is.null(title)) facet_title <- as.character(title[1])
  out <- as_plot_data("facet", list(
    facets = facet_tbl,
    title = facet_title,
    subtitle = if (is.null(facet)) "Ranked facet-level measures" else paste0("Ranked measures for facet: ", as.character(facet[1])),
    legend = new_plot_legend("Facet estimates", "location", "bar", style$fill_soft),
    reference_lines = new_reference_lines("h", 0, "Centered facet reference", "dashed", "reference"),
    preset = style$name
  ))
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    draw_facet_plot(
      facet_tbl,
      title = .mfrm_plot_display_title(
        out$data,
        title = facet_title,
        fallback = "Facet-level estimates"
      ),
      palette = resolve_palette(palette = palette, defaults = c(facet_bar = style$fill_soft)),
      show_ci = show_ci,
      ci_level = ci_level
    )
  }
  invisible(out)
}
