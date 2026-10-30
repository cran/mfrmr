# Plot methods centered on fitted many-facet Rasch objects.

.mfrm_fit_plot_readiness <- function(fit) {
  if (inherits(fit, "mfrm_imported_fit")) {
    .validate_imported_metric(fit)
    return(list(status = "review_only", ready = FALSE,
      detail = paste(fit$source$metric,
        "Imported point estimates retain the source scale; import does not establish the model's suitability for inference."),
      table = data.frame(Domain = "Source", Status = "imported_descriptive")))
  }
  convergence <- mfrm_convergence_state(fit)
  stored <- as.data.frame(
    mfrmr_get_readiness_record(fit)$fit, stringsAsFactors = FALSE
  )
  has_stored <- nrow(stored) == 1L &&
    "ReadinessContractVersion" %in% names(stored) &&
    identical(
      as.character(stored$ReadinessContractVersion[1]),
      mfrmr_readiness_contract_version()
    )
  stored_value <- function(field, default = NA_character_) {
    if (has_stored && field %in% names(stored)) stored[[field]][1] else default
  }
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
  numerical_status <- if (has_stored) {
    switch(
      as.character(stored_value("NumericalState", "legacy_unknown")),
      ready = "pass", review = "review", failed = "fail", not_run = "fail",
      legacy_unknown = "legacy_unknown", "not_assessed"
    )
  } else if (isTRUE(convergence$inference_ready)) {
    "pass"
  } else if (convergence_severity %in% c("fail", "error")) {
    "fail"
  } else {
    "review"
  }
  fit_status <- as.character(stored_value("FitReadiness", "legacy_unknown"))
  data_status <- if (has_stored) {
    switch(
      as.character(stored_value("InputState", "legacy_unknown")),
      pass = "pass", review = "review", blocked = "blocked",
      legacy_unknown = "legacy_unknown", "not_assessed"
    )
  } else {
    domain_status("Data")
  }
  design_status <- domain_status("Design")
  stability_status <- domain_status("Stability")
  interpretation_ready <- (!has_stored || identical(fit_status, "ready")) &&
    identical(numerical_status, "pass") &&
    identical(data_status, "pass") &&
    identical(design_status, "pass_linked") &&
    identical(stability_status, "pass")
  interpretation_status <- if (interpretation_ready) {
    "ready_for_diagnostic_interpretation"
  } else {
    "review_only"
  }
  statuses <- data.frame(
    Domain = c("Fit", "Numerical", "Data", "Design", "Stability", "Plot"),
    Status = c(
      fit_status,
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
      "Stored fit readiness plus numerical, data-support, connectivity, and stability checks ",
      "passed. Treat this display as diagnostic evidence, not automatic ",
      "publication approval."
    )
  } else {
    paste0(
      "Review-only display: Fit=", fit_status,
      ", Numerical=", numerical_status,
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
  if (!isTRUE(plot_data$display$show_title)) return("")
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
  if (isTRUE(review_only) && !identical(plot_data$display$show_notes, FALSE) &&
      !startsWith(display_title, "REVIEW ONLY")) {
    display_title <- paste0("REVIEW ONLY - ", display_title)
  }
  display_title
}

.mfrm_plot_notes <- function(data) {
  text <- c(data$interpretation_note, data$subtitle,
            data$curve_basis$Description, data$retention_note)
  type <- rep(c("interpretation", "description", "reference_profile", "retention"),
              lengths(list(data$interpretation_note, data$subtitle,
                           data$curve_basis$Description, data$retention_note)))
  if (identical(data$interpretation_status, "review_only")) type[1L] <- "warning"
  if (!is.null(data$locations)) {
    extra <- .wright_footer_notes(data, show_ci = isTRUE(data$show_ci))
    if (nrow(data$person_stats %||% data.frame()) > 0L) {
      extra <- c(extra, sprintf("Persons: N=%d; mean=%.2f; median=%.2f.",
        data$person_stats$N[1], data$person_stats$Mean[1], data$person_stats$Median[1]))
    }
    text <- c(text, extra)
    type <- c(type, rep("display", length(extra)))
  }
  out <- data.frame(Type = type, Text = as.character(text), stringsAsFactors = FALSE)
  out <- out[!is.na(out$Text) & nzchar(out$Text) & !duplicated(out$Text), , drop = FALSE]
  rownames(out) <- NULL
  out
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
  if (inherits(x, "mfrm_imported_fit")) {
    stop("Imported absolute thresholds do not supply a native response-model contract. Use the source package for probability curves.", call. = FALSE)
  }
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

  profiles <- data.frame(CurveGroup = "Common", StepLevel = "Common", SlopeLevel = NA_character_)
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
    separate_owners <- model == "GPCM" &&
      !identical(x$config$slope_facet, x$config$step_facet)
    profiles <- if (separate_owners) {
      expand.grid(StepLevel = ordered_levels,
        SlopeLevel = x$prep$levels[[x$config$slope_facet]], stringsAsFactors = FALSE)
    } else {
      data.frame(StepLevel = ordered_levels, SlopeLevel = ordered_levels,
        stringsAsFactors = FALSE)
    }
    profiles$CurveGroup <- if (separate_owners) make.unique(paste0(
      step_facet, " = ", profiles$StepLevel, "; ", x$config$slope_facet,
      " = ", profiles$SlopeLevel)) else profiles$StepLevel
    for (i in seq_len(nrow(profiles))) {
      lvl <- profiles$StepLevel[i]
      group <- profiles$CurveGroup[i]
      sub <- step_tbl[as.character(step_tbl$StepFacet) == as.character(lvl), , drop = FALSE]
      if (nrow(sub) == 0) next
      ord <- order(step_index_from_label(sub$Step))
      tau <- as.numeric(sub$Estimate[ord])
      slope_val <- if (model == "GPCM") {
        sval <- unname(slope_lookup[profiles$SlopeLevel[i]])
        if (!is.finite(sval) || sval <= 0) {
          stop("GPCM pathway/CCC plots require finite positive slopes for every slope-facet level.")
        }
        sval
      } else {
        1
      }
      groups[[group]] <- list(
        name = group,
        step_cum = c(0, cumsum(tau)),
        tau = tau,
        slope = slope_val
      )
      step_points <- dplyr::bind_rows(
        step_points,
        tibble::tibble(
          CurveGroup = group,
          StepLevel = as.character(lvl),
          Step = as.character(sub$Step[ord]),
          StepIndex = seq_along(tau),
          Threshold = tau,
          Slope = slope_val
        )
      )
    }
  }

  if (length(groups) == 0) stop("No step groups available for pathway/CCC plots.")
  profiles$StepOwner <- if (model == "RSM") "Common" else x$config$step_facet
  profiles$SlopeOwner <- if (model == "GPCM") x$config$slope_facet else NA_character_
  if (model != "GPCM") profiles$SlopeLevel <- NA_character_

  list(
    model = model,
    step_facet = if (model == "RSM") NA_character_ else as.character(x$config$step_facet[1]),
    categories = categories,
    rating_min = rating_min,
    n_cat = n_cat,
    profiles = profiles,
    curve_basis = "zero_additive_facet_profile",
    predictor_offset = 0,
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
                                endpoint_labels = data.frame(),
                                joint_profiles = FALSE) {
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
  match_idx <- if (joint_profiles) rep(NA_integer_, length(curve_groups)) else
    match(curve_groups, as.character(candidate_tbl$Level))
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
  curve_fit_status$ReviewReason[!curve_fit_status$MatchedFitRow] <- if (joint_profiles)
    "A joint slope/step profile has no single-facet fit statistic; inspect fit_measures separately." else
    "No fit-measure row matched this pathway curve group."

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
  curve_basis <- as.character(
    curve_spec$curve_basis %||% "zero_additive_facet_profile"
  )[1]
  predictor_offset <- suppressWarnings(as.numeric(
    curve_spec$predictor_offset %||% 0
  )[1])
  if (!is.finite(predictor_offset)) {
    stop("Category curves require a finite predictor offset.")
  }
  for (g in names(curve_spec$groups)) {
    grp <- curve_spec$groups[[g]]
    group_offset <- grp$predictor_offset %||% predictor_offset
    if (length(group_offset) != 1L || !is.finite(group_offset)) {
      stop("Category curves require a finite predictor offset.")
    }
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
        eta = theta_grid + group_offset,
        step_cum_mat = matrix(grp$step_cum, nrow = 1L),
        criterion_idx = rep(1L, length(theta_grid)),
        slopes = slope_val,
        slope_idx = rep(1L, length(theta_grid))
      )
    } else {
      category_prob_rsm(theta_grid + group_offset, grp$step_cum)
    }
    k_vals <- as.numeric(curve_spec$categories)
    expected <- as.numeric(probs %*% matrix(k_vals, ncol = 1))
    score_variance <- mfrm_category_variance(probs, k_vals, expected)
    information <- (slope_val^2) * score_variance
    exp_tables[[idx_exp]] <- tibble::tibble(
      Theta = theta_grid,
      ExpectedScore = expected,
      ScoreVariance = score_variance,
      Information = information,
      Slope = slope_val,
      Model = curve_spec$model,
      CurveGroup = grp$name,
      CurveBasis = curve_basis,
      PredictorOffset = group_offset
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
        CurveGroup = grp$name,
        CurveBasis = curve_basis,
        PredictorOffset = group_offset
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
    mfrm_results_validate_diagnostics_identity(x, diagnostics,
                                               helper = "fit plot")
    measures <- as.data.frame(diagnostics$measures %||% data.frame(), stringsAsFactors = FALSE)
    if (!all(c("Facet", "Level", "Estimate") %in% names(measures))) {
      stop("`diagnostics$measures` must contain Facet, Level, and Estimate columns.", call. = FALSE)
    }
    if (!"SE" %in% names(measures)) measures$SE <- NA_real_
    keep <- intersect(
      c(
        "Facet", "Level", "Estimate", "SE", "SE_Method", "PrecisionTier",
        "SupportsFormalInference", "SEUse", "CIBasis", "CIUse",
        "CIEligible", "CILabel", "CI_Method"
      ),
      names(measures)
    )
    se_tbl <- measures[, keep, drop = FALSE]
    se_tbl$Facet <- as.character(se_tbl$Facet)
    se_tbl$Level <- as.character(se_tbl$Level)
    se_tbl$Estimate <- suppressWarnings(as.numeric(se_tbl$Estimate))
    se_tbl$SE <- suppressWarnings(as.numeric(se_tbl$SE))
    z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
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
    return(apply_fixed_measure_precision(se_tbl, x$config))
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
    z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
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
    se_tbl$PrecisionTier <- "exploratory"
    se_tbl$SupportsFormalInference <- FALSE
    se_tbl$SEUse <- "screening_only"
    se_tbl$CIBasis <-
      "Normal interval from exploratory observation-table SE"
    se_tbl$CIUse <- "screening_only"
    se_tbl$CIEligible <- FALSE
    se_tbl$CILabel <- "Approximate interval; screening only"
    se_tbl$Measure_Source <- "fit + observation-table information"
    as.data.frame(apply_fixed_measure_precision(se_tbl, x$config), stringsAsFactors = FALSE)
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
    if (any(as.character(
      (plot_core$group_summary %||% data.frame())$PlotType %||% character(0)
    ) == "Facet level")) "iqr_range",
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
    step_threshold = "line-point",
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

.spread_wright_label_positions <- function(values, lower, upper, min_gap = 0) {
  values <- suppressWarnings(as.numeric(values))
  n <- length(values)
  if (n <= 1L || !is.finite(lower) || !is.finite(upper) || upper <= lower) {
    return(pmin(upper, pmax(lower, values)))
  }
  span <- upper - lower
  gap <- min(max(0.045 * span, 0.08, min_gap), 0.92 * span / (n - 1L))
  ord <- order(values, na.last = TRUE)
  placed <- pmin(upper, pmax(lower, values[ord]))
  for (i in 2:n) {
    placed[i] <- max(placed[i], placed[i - 1L] + gap)
  }
  if (placed[n] > upper) {
    placed <- placed - (placed[n] - upper)
  }
  for (i in seq.int(n - 1L, 1L)) {
    placed[i] <- min(placed[i], placed[i + 1L] - gap)
  }
  if (placed[1L] < lower) {
    placed <- placed + (lower - placed[1L])
  }
  out <- numeric(n)
  out[ord] <- pmin(upper, pmax(lower, placed))
  out
}

.place_plot_labels <- function(x, y, labels, cex, point_x, point_y, adj = 0.5,
                                xlim = graphics::par("usr")[1:2],
                                ylim = graphics::par("usr")[3:4]) {
  width <- graphics::strwidth(labels, cex = cex)
  height <- graphics::strheight(labels, cex = cex)
  pad_x <- graphics::strwidth("M", cex = cex) * 0.3
  pad_y <- graphics::strheight("M", cex = cex) * 0.3
  center_x <- x + (0.5 - adj) * width
  half_width <- width / 2 + pad_x
  half_height <- height / 2 + pad_y
  center_x <- pmax(xlim[1] + half_width, pmin(xlim[2] - half_width, center_x))
  out <- data.frame(X = center_x, Y = y,
                    Left = center_x - half_width, Right = center_x + half_width,
                    Bottom = y - half_height, Top = y + half_height)
  points <- data.frame(Left = point_x - pad_x, Right = point_x + pad_x,
                       Bottom = point_y - pad_y, Top = point_y + pad_y)
  unresolved <- logical(length(y))
  # ponytail: greedy vertical placement; an overcrowded canvas still needs more space.
  for (i in seq_along(y)) {
    previous <- rbind(points, out[seq_len(i - 1L), names(points), drop = FALSE])
    previous$Point <- seq_len(nrow(previous)) <= nrow(points)
    previous <- previous[previous$Left < out$Right[i] & previous$Right > out$Left[i], , drop = FALSE]
    candidates <- c(y[i], previous$Top + half_height[i],
                     previous$Bottom - half_height[i])
    candidates <- unique(pmax(ylim[1] + half_height[i],
                               pmin(ylim[2] - half_height[i], candidates)))
    candidates <- candidates[order(abs(candidates - y[i]))]
    overlap <- vapply(candidates, function(candidate) {
      amount <- pmax(0, pmin(candidate + half_height[i], previous$Top) -
                         pmax(candidate - half_height[i], previous$Bottom))
      c(sum(amount[previous$Point]), sum(amount[!previous$Point]))
    }, numeric(2))
    overlap[overlap < diff(ylim) * 1e-10] <- 0
    best <- order(overlap[1, ], overlap[2, ], seq_along(candidates))[1]
    unresolved[i] <- any(overlap[, best] > 0) ||
      2 * half_width[i] > diff(xlim) || 2 * half_height[i] > diff(ylim)
    out$Y[i] <- candidates[best]
    out$Bottom[i] <- out$Y[i] - half_height[i]
    out$Top[i] <- out$Y[i] + half_height[i]
  }
  if (any(unresolved)) {
    warning(sprintf(
      "Plot labels need more space (%d of %d could not be separated). Enlarge the plotting window or export a larger figure; plotted values are retained.",
      sum(unresolved), length(y)), call. = FALSE)
  }
  out
}

.build_wright_label_layout <- function(point_tbl, y_range) {
  out <- as.data.frame(point_tbl, stringsAsFactors = FALSE)
  if (nrow(out) == 0L) return(tibble::as_tibble(out))
  out$LabelY <- suppressWarnings(as.numeric(out$DisplayEstimate))
  out$LabelSide <- "right"
  out$LabelX <- suppressWarnings(as.numeric(out$X)) + 0.075
  out$LabelHjust <- 0
  out$LabelText <- as.character(out$DisplayLabel)

  groups <- unique(as.character(out$Group))
  for (group in groups) {
    idx <- which(as.character(out$Group) == group)
    ord <- idx[order(out$DisplayEstimate[idx], as.character(out$Label[idx]))]
    out$LabelY[ord] <- .spread_wright_label_positions(
      out$DisplayEstimate[ord], y_range[1L], y_range[2L]
    )
    is_step <- as.character(out$PlotType[ord]) == "Step threshold"
    side <- ifelse(
      is_step,
      rep(c("left", "right"), length.out = length(ord)),
      ifelse(out$X[ord] <= out$XBase[ord], "left", "right")
    )
    if (length(ord) == 1L && !is_step[1L]) side <- "right"
    out$LabelSide[ord] <- side
    out$LabelX[ord] <- out$X[ord] + ifelse(side == "right", 0.075, -0.075)
    out$LabelHjust[ord] <- ifelse(side == "right", 0, 1)
    out$LabelText[ord] <- ifelse(
      is_step,
      sprintf(
        "%s (%+.2f)",
        as.character(out$DisplayLabel[ord]),
        suppressWarnings(as.numeric(out$OriginalEstimate[ord]))
      ),
      as.character(out$DisplayLabel[ord])
    )
  }
  out$LabelDisplaced <- abs(out$LabelY - out$DisplayEstimate) >
    sqrt(.Machine$double.eps)
  tibble::as_tibble(out)
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
                                  persons_per_star = NULL,
                                  highlight = NULL) {
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
  source_fit_state <- if (
    "SourceFitReadiness" %in% names(person_tbl) && nrow(person_tbl) > 0L
  ) {
    as.character(person_tbl$SourceFitReadiness[1L])
  } else {
    fit_readiness <- as.data.frame(
      mfrmr_get_readiness_record(x)$fit %||% data.frame(),
      stringsAsFactors = FALSE
    )
    if (nrow(fit_readiness) == 1L &&
        "FitReadiness" %in% names(fit_readiness)) {
      as.character(fit_readiness$FitReadiness[1L])
    } else {
      "legacy_unknown"
    }
  }
  prior_regularized_extreme <- if ("ReasonCodes" %in% names(person_tbl)) {
    grepl(
      "(^|;)mml_extreme_response_prior_regularized($|;)",
      as.character(person_tbl$ReasonCodes)
    )
  } else {
    rep(FALSE, nrow(person_tbl))
  }
  blocked_extreme <- identical(source_fit_state, "blocked") &
    prior_regularized_extreme
  person_exclusions <- person_tbl[blocked_extreme, , drop = FALSE]
  if (nrow(person_exclusions) > 0L) {
    person_exclusions$PlotExclusionReason <-
      "prior_regularized_extreme_eap_from_blocked_source_fit"
  }
  typed_extreme <- if ("ParameterStatus" %in% names(person_tbl)) {
    as.character(person_tbl$ParameterStatus) %in%
      c("unbounded_low", "unbounded_high")
  } else {
    rep(FALSE, nrow(person_tbl))
  }
  keep_person <- !blocked_extreme & (is.finite(person_tbl$Estimate) |
    (identical(wright_style, "facets_style") & typed_extreme)
  )
  person_tbl <- person_tbl[keep_person, , drop = FALSE]
  if (nrow(person_tbl) == 0L) {
    if (nrow(person_exclusions) > 0L) {
      stop(
        paste(
          "A Wright map is unavailable because every Person estimate is a",
          "prior-regularized extreme EAP from a blocked source fit.",
          "Inspect `summary(fit)$decision` and refit before plotting."
        ),
        call. = FALSE
      )
    }
    stop("Person estimates are not available for Wright map.")
  }

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (!all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available for Wright map.")
  }
  highlight_keys <- character()
  if (!is.null(highlight)) {
    if (!is.list(highlight) || is.data.frame(highlight) ||
        length(highlight) == 0L || is.null(names(highlight)) ||
        anyNA(names(highlight)) || any(!nzchar(names(highlight))) ||
        anyDuplicated(names(highlight))) {
      stop("`highlight` must be NULL or a named list, e.g. list(Rater = c('R01', 'R02')).", call. = FALSE)
    }
    for (facet_name in names(highlight)) {
      selected <- highlight[[facet_name]]
      if (!is.character(selected) || length(selected) == 0L ||
          anyNA(selected) || any(!nzchar(selected))) {
        stop("Each `highlight` entry must contain non-empty character level IDs.", call. = FALSE)
      }
      available <- facet_tbl[facet_tbl$Facet == facet_name, , drop = FALSE]
      if (nrow(available) == 0L) {
        stop("Unknown highlight facet: ", facet_name, call. = FALSE)
      }
      if (any(!selected %in% available$Level)) {
        stop("Unknown highlight level(s) in ", facet_name, ": ",
          paste(setdiff(selected, available$Level), collapse = ", "), call. = FALSE)
      }
      if (any(!is.finite(available$Estimate[match(selected, available$Level)]))) {
        stop("Highlighted levels must have finite fitted locations.", call. = FALSE)
      }
      highlight_keys <- c(highlight_keys, paste(facet_name, unique(selected), sep = "\r"))
    }
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
        "Facet", "Level", "SE", "Fixed",
        "CI_Level", "SE_Method", "PrecisionTier", "SupportsFormalInference",
        "SEUse", "CIBasis", "CIUse", "CIEligible", "CILabel",
        "Measure_Source"
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
      z <- stats::qnorm((1 - ci_level_value) / 2, lower.tail = FALSE)
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
  point_tbl$Highlighted <- point_tbl$PlotType == "Facet level" &
    point_key %in% highlight_keys
  all_location_estimates <- point_tbl$Estimate

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
    facet_limit <- top_n
    # Keep the ordinary compact view and add requested levels if necessary.
    # Selection must not silently substitute one rater for another.
    facet_rows <- facet_rows[
      seq_len(nrow(facet_rows)) <= facet_limit | facet_rows$Highlighted,
      , drop = FALSE
    ]
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
  if (length(highlight_keys) > 0L) {
    retention_note <- paste(c(retention_note[nzchar(retention_note)],
      paste0("Labels prefixed with '> ' identify user-selected facet levels; ",
        "selection is not a statistical flag. Selected levels are retained in ",
        "addition to the usual top_n selection, so the count may exceed top_n. ",
        "Fitted locations and uncertainty are unchanged.")), collapse = " ")
  }
  if (nrow(person_exclusions) > 0L) {
    person_note <- paste0(
      nrow(person_exclusions),
      " prior-regularized extreme Person EAP(s) from the blocked source fit ",
      "were omitted from the plotted scale; exact values remain in ",
      "`fit$facets$person` and `person_exclusions`."
    )
    retention_note <- paste(
      c(retention_note[nzchar(retention_note)], person_note),
      collapse = " "
    )
  }

  group_levels <- unique(point_tbl$Group)
  point_tbl <- point_tbl |>
    dplyr::mutate(Group = factor(.data$Group, levels = group_levels)) |>
    dplyr::group_by(.data$Group) |>
    dplyr::arrange(.data$Estimate, .by_group = TRUE) |>
    dplyr::mutate(
      XBase = as.numeric(.data$Group),
      X = if (dplyr::n() == 1) .data$XBase else .data$XBase + seq(-0.18, 0.18, length.out = dplyr::n()),
      X = ifelse(.data$PlotType == "Step threshold", .data$XBase, .data$X)
    ) |>
    dplyr::ungroup()

  finite_person_estimates <- person_tbl$Estimate[is.finite(person_tbl$Estimate)]
  finite_person_center <- if (length(finite_person_estimates) > 0L) {
    mean(finite_person_estimates)
  } else {
    NA_real_
  }
  hist_data <- if (length(finite_person_estimates) > 0L) {
    # Automatic breaks for one estimate near numerical zero can be subpixel.
    hist(finite_person_estimates,
      breaks = if (length(finite_person_estimates) == 1L)
        finite_person_center + c(-0.5, 0.5) else "FD", plot = FALSE)
  } else {
    structure(
      list(
        breaks = c(-0.5, 0.5), counts = 0L, density = 0,
        intensities = 0, mids = 0, xname = "finite Person estimates",
        equidist = TRUE
      ),
      class = "histogram"
    )
  }
  if (length(finite_person_estimates) == 1L && wright_style == "native") retention_note <- paste(
    c(retention_note[nzchar(retention_note)],
      "The single-Person histogram bin has width one scale unit for display; it is not an uncertainty interval."),
    collapse = " ")
  finite_range_basis <- c(all_location_estimates, finite_person_estimates)
  finite_range_basis <- finite_range_basis[is.finite(finite_range_basis)]
  y_range <- if (length(finite_range_basis) > 0L) {
    range(finite_range_basis)
  } else {
    c(-0.5, 0.5)
  }
  if (!all(is.finite(y_range)) || diff(y_range) <= sqrt(.Machine$double.eps)) {
    center <- if (length(finite_range_basis) > 0L) {
      mean(finite_range_basis)
    } else {
      0
    }
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
  point_tbl$DisplayLabel[point_tbl$Highlighted] <- paste0(
    "> ", point_tbl$DisplayLabel[point_tbl$Highlighted]
  )
  point_tbl <- tibble::as_tibble(add_wright_ci_display_metadata(
    point_tbl,
    lower = y_range[1],
    upper = y_range[2]
  ))

  # Every retained location remains identifiable. Collision avoidance changes
  # label coordinates rather than silently dropping interior facet levels.
  label_tbl <- .build_wright_label_layout(
    point_tbl |>
      dplyr::distinct(.data$Group, .data$Label, .keep_all = TRUE),
    y_range = y_range
  )

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
      TargetGap = .data$Median - finite_person_center,
      DisplayTargetGap = .data$DisplayMedian - finite_person_center
    )

  person_stats <- tibble::tibble(
    N = nrow(person_tbl),
    ReviewExcludedN = nrow(person_exclusions),
    FiniteN = length(finite_person_estimates),
    BoundaryExcludedN = nrow(person_tbl) - length(finite_person_estimates),
    Mean = finite_person_center,
    Median = if (length(finite_person_estimates) > 0L) {
      stats::median(finite_person_estimates)
    } else {
      NA_real_
    },
    SD = if (length(finite_person_estimates) > 1L) {
      stats::sd(finite_person_estimates)
    } else {
      NA_real_
    }
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
    person_exclusions = person_exclusions,
    person_hist = hist_data,
    person_stats = person_stats,
    locations = point_tbl,
    label_points = label_tbl,
    group_summary = group_summary,
    group_levels = group_levels,
    y_range = y_range,
    display_settings = display_settings,
    label_limit = top_n,
    highlight = highlight,
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

.wright_footer_notes <- function(plot_data, show_ci = FALSE) {
  loc <- plot_data$locations
  boundary_end <- as.character(loc$BoundaryEnd %||% "none")
  ci_suppressed <- as.logical(loc$CISuppressed %||% FALSE)
  ci_clipped <- as.logical(loc$CIClipped %||% FALSE)
  ci_suppressed[is.na(ci_suppressed)] <- FALSE
  ci_clipped[is.na(ci_clipped)] <- FALSE
  retention_note <- as.character(plot_data$retention_note %||% "")
  footer_notes <- character()
  if (identical(plot_data$wright_style, "facets_style")) {
    settings <- plot_data$facets_style$settings
    footer_notes <- c(
      sprintf("* = %s person(s); rows/%s = %d", settings$PersonsPerStar[1],
              if (is.null(plot_data$axis_label)) "logit" else "scale unit", settings$RowsPerLogit[1]),
      if (is.null(plot_data$axis_label)) {
        "FACETS Table 6-style visual layout; estimates remain mfrmr estimates (not numerical equivalence)."
      } else {
        "FACETS Table 6-style visual layout; estimates retain the source scale."
      })
    if (isTRUE(show_ci)) footer_notes <- c(footer_notes,
      sprintf("Whiskers show %g%% mfrmr confidence intervals; triangles mark bounds beyond the displayed ruler.",
              round(100 * (plot_data$ci_level %||% 0.95))))
  }
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
  footer_notes
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
  default_title <- if (!is.null(plot_data$axis_label)) {
    "Source person and item/design-cell estimates"
  } else if (isTRUE(facets_style)) {
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
      facet_level = "#0072B2",
      step_threshold = "#D55E00",
      person_hist = "gray80",
      grid = "#ececec",
      range = "#94a3b8",
      iqr = "#334155"
    )
  )
  old_par <- graphics::par()[c("mfrow", "cex", "mex", "mar", "oma", "mgp")]
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.15, 2.55))
  device_scale <- min(1, graphics::par("din")[1] / 7, graphics::par("din")[2] / 5)
  graphics::par(mgp = c(2.6, 0.85, 0), cex = old_par$cex * device_scale,
                mex = device_scale)

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
  ci_suppressed <- if ("CISuppressed" %in% names(loc)) {
    value <- as.logical(loc$CISuppressed)
    value[is.na(value)] <- FALSE
    value
  } else {
    rep(FALSE, nrow(loc))
  }
  footer_notes <- if (!identical(plot_data$display$show_notes, FALSE)) {
    .wright_footer_notes(plot_data, show_ci = show_ci)
  } else character()
  footer_lines <- unlist(lapply(
    footer_notes,
    function(note) strwrap(note, width = max(40L, floor(graphics::par("din")[1] * 17)))
  ), use.names = FALSE)
  x_lab <- truncate_axis_label(plot_data$group_levels, width = 16L)
  angle <- label_angle * pi / 180
  label_height <- max(abs(sin(angle)) * graphics::strwidth(x_lab, units = "inches", cex = 0.82) +
    abs(cos(angle)) * graphics::strheight(x_lab, units = "inches", cex = 0.82))
  # Ambiguous abbreviations retain full labels, so their margin must grow too.
  bottom_margin <- max(6.6, (label_height + 0.1 * graphics::par("din")[2]) /
    (graphics::par("csi") * graphics::par("mex")))
  graphics::par(oma = c(0.6 + 0.7 * length(footer_lines), 0,
                        if (nzchar(display_title)) 2.8 else 0.2, 0))
  pretty_y <- pretty(yr, n = 6)
  hist_data <- plot_data$person_hist
  hist_max <- max(hist_data$counts, na.rm = TRUE)
  if (!is.finite(hist_max) || hist_max <= 0) hist_max <- 1

  graphics::par(mar = c(bottom_margin, 4.6, 3.0, 0.8))
  plot(
    x = c(0, hist_max * 1.12),
    y = yr,
    type = "n",
    xlab = "Persons",
    ylab = plot_data$axis_label %||% "Logit scale",
    main = "Persons",
    yaxt = "n"
  )
  axis(2, at = pretty_y, las = 1)
  abline(h = pretty_y, col = pal["grid"], lty = 1)
  for (i in seq_along(hist_data$counts)) {
    rect(
      xleft = 0,
      ybottom = hist_data$breaks[i],
      xright = hist_data$counts[i],
      ytop = hist_data$breaks[i + 1],
      col = grDevices::adjustcolor(pal["person_hist"], alpha.f = 0.85),
      border = "white"
    )
  }
  if (nrow(plot_data$person_stats) > 0) {
    abline(h = plot_data$person_stats$Median[1], col = grDevices::adjustcolor("gray25", alpha.f = 0.9), lty = 2, lwd = 1.2)
    abline(h = plot_data$person_stats$Mean[1], col = grDevices::adjustcolor("gray45", alpha.f = 0.9), lty = 3, lwd = 1.2)
    if (!identical(plot_data$display$show_notes, FALSE)) mtext(
      sprintf(
        "N=%d  Mean=%.2f  Median=%.2f",
        plot_data$person_stats$N[1],
        plot_data$person_stats$Mean[1],
        plot_data$person_stats$Median[1]
      ),
      side = 3,
      line = 0.2,
      cex = 0.65 * device_scale,
      adj = 0
    )
  }

  group_density <- plot_data$group
  if (!is.null(group_density) && nrow(group_density) > 0L) {
    groups <- unique(group_density$Group)
    cols <- .plot_series_colors(groups, plot_data$preset %||% "standard")
    line_types <- .plot_series_linetypes(groups)
    density_max <- max(group_density$Density, na.rm = TRUE)
    for (i in seq_along(groups)) {
      sub <- group_density[group_density$Group == groups[i], , drop = FALSE]
      if (nrow(sub) < 2L) next
      lines(hist_max * sub$Density / density_max, sub$Theta,
                      col = cols[i], lwd = 2, lty = line_types[i])
    }
    legend("topright", legend = groups, col = cols, lwd = 2, lty = line_types,
                     bty = "n", cex = 0.65, title = "Scaled density", title.cex = 0.6)
  }

  xr <- c(0.5, length(plot_data$group_levels) + 0.75)
  # Reserve a dedicated right-hand gutter for the legend. Keeping the legend
  # outside the measurement ruler prevents it from hiding high step thresholds
  # (especially when a GPCM supplies one step ladder per curve group).
  legend_spec <- .wright_native_legend_spec(plot_data, show_ci, ci_level)
  legend_margin <- max(graphics::strwidth(legend_spec$label, units = "inches", cex = 0.72)) /
    (graphics::par("csi") * graphics::par("mex")) + 4 / device_scale
  graphics::par(mar = c(bottom_margin, 2.4, 3.0, legend_margin))
  plot(
    x = loc$X,
    y = loc_y,
    type = "n",
    xlim = xr,
    ylim = yr,
    xaxt = "n",
    xlab = "",
    yaxt = "n",
    ylab = "",
    main = "Locations"
  )
  for (i in seq_along(plot_data$group_levels)) {
    if (i %% 2 == 1) {
      rect(
        xleft = i - 0.5,
        ybottom = yr[1] - 10,
        xright = i + 0.5,
        ytop = yr[2] + 10,
        col = grDevices::adjustcolor(pal["grid"], alpha.f = 0.18),
        border = NA
      )
    }
  }
  abline(h = pretty_y, col = pal["grid"], lty = 1)
  abline(h = 0, col = grDevices::adjustcolor("gray35", alpha.f = 0.9), lty = 2)
  group_summary <- plot_data$group_summary
  if (!is.null(group_summary) && nrow(group_summary) > 0) {
    group_summary <- group_summary[
      as.character(group_summary$PlotType) == "Facet level",
      ,
      drop = FALSE
    ]
  }
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
    segments(
      x0 = group_summary$XBase,
      y0 = summary_value("DisplayMin", "Min"),
      x1 = group_summary$XBase,
      y1 = summary_value("DisplayMax", "Max"),
      col = range_cols,
      lwd = 1.6
    )
    segments(
      x0 = group_summary$XBase,
      y0 = summary_value("DisplayQ1", "Q1"),
      x1 = group_summary$XBase,
      y1 = summary_value("DisplayQ3", "Q3"),
      col = iqr_cols,
      lwd = 5.2
    )
    segments(
      x0 = group_summary$XBase - 0.13,
      y0 = summary_value("DisplayMedian", "Median"),
      x1 = group_summary$XBase + 0.13,
      y1 = summary_value("DisplayMedian", "Median"),
      col = grDevices::adjustcolor("gray15", alpha.f = 0.95),
      lwd = 1.4
    )
  }
  axis(4, at = pretty_y, labels = pretty_y, las = 1)
  x_at <- seq_along(plot_data$group_levels)
  draw_rotated_x_labels(
    at = x_at,
    labels = x_lab,
    srt = label_angle,
    cex = 0.82,
    line_offset = 0.085
  )
  mtext(display_title, side = 3, outer = TRUE, line = 0.9,
                  font = 2, cex = 0.95 * device_scale)
  cols <- ifelse(loc$PlotType == "Step threshold", pal["step_threshold"], pal["facet_level"])
  pch <- ifelse(loc$PlotType == "Step threshold", 17, 16)
  step_loc <- loc[as.character(loc$PlotType) == "Step threshold", , drop = FALSE]
  if (nrow(step_loc) > 1L) {
    for (group in unique(as.character(step_loc$Group))) {
      sub <- step_loc[as.character(step_loc$Group) == group, , drop = FALSE]
      if (nrow(sub) > 1L) {
        segments(
          x0 = sub$XBase[1L], y0 = min(sub$DisplayEstimate, na.rm = TRUE),
          x1 = sub$XBase[1L], y1 = max(sub$DisplayEstimate, na.rm = TRUE),
          col = grDevices::adjustcolor(pal["step_threshold"], alpha.f = 0.48),
          lwd = 1.2
        )
      }
    }
  }
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
      z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
      ci_ok <- is.finite(loc$SE) & loc$SE > 0
      ci_lo <- pmin(yr[2], pmax(yr[1], loc$Estimate - z * loc$SE))
      ci_hi <- pmin(yr[2], pmax(yr[1], loc$Estimate + z * loc$SE))
    } else {
      ci_ok <- rep(FALSE, nrow(loc))
      ci_lo <- ci_hi <- rep(NA_real_, nrow(loc))
    }
    regular_ci <- ci_ok & !ci_suppressed
    if (any(regular_ci)) {
      arrows(
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
        points(
          loc$X[clipped_lower], rep(yr[1], sum(clipped_lower)),
          pch = 25, bg = "white", col = cols[clipped_lower], cex = 0.78
        )
      }
      if (any(clipped_upper)) {
        points(
          loc$X[clipped_upper], rep(yr[2], sum(clipped_upper)),
          pch = 24, bg = "white", col = cols[clipped_upper], cex = 0.78
        )
      }
    }
  }
  points(
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
    label_x <- suppressWarnings(as.numeric(label_tbl$LabelX %||% label_tbl$X))
    label_y <- suppressWarnings(as.numeric(
      label_tbl$LabelY %||% label_tbl$DisplayEstimate %||% label_tbl$Estimate
    ))
    label_text <- as.character(
      label_tbl$LabelText %||% label_tbl$DisplayLabel %||% label_tbl$Label
    )
    label_side <- as.character(
      label_tbl$LabelSide %||% rep("right", nrow(label_tbl))
    )
    point_y <- suppressWarnings(as.numeric(
      label_tbl$DisplayEstimate %||% label_tbl$Estimate
    ))
    label_col <- ifelse(
      as.character(label_tbl$PlotType) == "Step threshold",
      pal["step_threshold"],
      "gray20"
    )
    label_text <- truncate_axis_label(label_text, width = 22L)
    placed <- .place_plot_labels(label_x, label_y, label_text, cex = 0.70,
      point_x = loc$X, point_y = loc_y,
      adj = ifelse(label_side == "left", 1, 0), ylim = yr)
    segments(
      x0 = label_tbl$X, y0 = point_y,
      x1 = pmax(placed$Left, pmin(placed$Right, label_tbl$X)),
      y1 = pmax(placed$Bottom, pmin(placed$Top, point_y)),
      col = grDevices::adjustcolor(label_col, alpha.f = 0.62),
      lwd = 0.7,
      xpd = NA
    )
    text(placed$X, placed$Y, labels = label_text, cex = 0.70,
                    col = "gray15")
  }
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
    legend_spec$key %in% c("step_threshold", "iqr_range", "ci_whisker"),
    1,
    0
  )
  legend_lwd <- ifelse(
    legend_spec$key == "iqr_range",
    4,
    ifelse(legend_spec$key %in% c("step_threshold", "ci_whisker"), 1.2, NA_real_)
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
  ruler_usr <- graphics::par("usr")
  legend(
    x = graphics::grconvertX(
      graphics::grconvertX(ruler_usr[2], from = "user", to = "inches") + 0.26,
      from = "inches", to = "user"),
    y = ruler_usr[4],
    xjust = 0,
    yjust = 1,
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
    cex = 0.72,
    xpd = NA
  )
  if (length(footer_lines) > 0L) {
    footer_col <- rep("#9A3412", length(footer_lines))
    for (i in seq_along(footer_lines)) {
      mtext(
        footer_lines[i],
        side = 1,
        line = 0.1 + 0.7 * (i - 1L),
        outer = TRUE,
        adj = 0,
        cex = 0.68 * device_scale,
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
    Level = ifelse(as.character(step_df$CurveGroup) == "Common", NA_character_,
      as.character(step_df$StepLevel %||% step_df$CurveGroup)),
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
    endpoint_labels = endpoint_labels,
    joint_profiles = identical(x$config$model, "GPCM") &&
      !identical(x$config$slope_facet, x$config$step_facet)
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
    curve_basis = .category_curve_basis(curve_spec),
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
  defaults <- .plot_series_colors(groups, plot_data$preset %||% "standard")
  line_types <- .plot_series_linetypes(groups)
  cols <- resolve_palette(palette = palette, defaults = defaults)
  y_rng <- range(exp_df$ExpectedScore, finite = TRUE)
  y_breaks <- pretty(y_rng, n = 6)
  dominance_regions <- plot_data$dominance_regions
  use_strips <- !is.null(dominance_regions) && nrow(dominance_regions) > 0 && length(groups) <= 6
  band_pad <- if (isTRUE(use_strips)) diff(y_rng) * (0.08 + 0.05 * length(groups)) else 0
  old_par <- graphics::par()[c("mar", "cex", "mex")]
  on.exit(graphics::par(old_par), add = TRUE)
  device_scale <- min(1, graphics::par("fin")[1] / 7, graphics::par("fin")[2] / 5)
  graphics::par(cex = old_par$cex * device_scale, mex = old_par$mex * device_scale)
  label_width <- max(graphics::strwidth(truncate_axis_label(groups, width = 14L),
                                       units = "inches", cex = 0.78))
  label_margin <- label_width / (graphics::par("csi") * graphics::par("mex")) + 1.5
  basis_note <- as.character(plot_data$curve_basis$Description %||% "")
  if (identical(plot_data$display$show_notes, FALSE)) basis_note <- ""
  basis_lines <- strwrap(basis_note,
    width = max(30L, floor(graphics::par("fin")[1] * 13 / device_scale)))
  note_start <- if (use_strips) 5.8 else 4.1
  graphics::par(mar = c(if (identical(plot_data$display$show_notes, FALSE)) 5.1 else max(if (use_strips) 6.2 else 5.1,
                            note_start + 0.7 * length(basis_lines) + 1.4),
                        max(4.6, label_margin),
                        if (identical(plot_data$display$show_title, FALSE)) 1 else 3.2,
                        if (length(groups) <= 5) max(2.1, label_margin) else 2.1))
  plot(
    x = range(exp_df$Theta, finite = TRUE),
    y = c(y_rng[1] - band_pad, y_rng[2]),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Expected score",
    yaxt = "n",
    main = .mfrm_plot_display_title(
      plot_data,
      title = title,
      fallback = "Expected-score pathway"
    )
  )
  axis(2, at = y_breaks[y_breaks >= y_rng[1] & y_breaks <= y_rng[2]], las = 1)
  abline(h = y_breaks, col = grDevices::adjustcolor("#d9dde3", alpha.f = 0.8), lty = 1)
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
      rect(
        xleft = dominance_regions$ThetaStart[i],
        ybottom = y_bottom,
        xright = dominance_regions$ThetaEnd[i],
        ytop = y_top,
        col = grDevices::adjustcolor(cols[grp], alpha.f = 0.14),
        border = "white"
      )
      text(
        x = dominance_regions$ThetaMid[i],
        y = (y_bottom + y_top) / 2,
        labels = as.character(dominance_regions$Category[i]),
        cex = 0.72,
        col = "gray15"
      )
    }
    text(
      x = rep(min(exp_df$Theta, na.rm = TRUE), length(groups)),
      y = band_top - (seq_along(groups) - 1) * (band_h + band_gap) - band_h / 2,
      labels = truncate_axis_label(groups, width = 14L),
      pos = 2,
      xpd = NA,
      cex = 0.75,
      col = "gray15"
    )
    if (!identical(plot_data$display$show_notes, FALSE)) {
      mtext("Dominant category by theta", side = 1, line = 4.6, adj = 0)
    }
  }
  if (!is.null(plot_data$steps) && nrow(plot_data$steps) > 0) {
    step_groups <- as.character(plot_data$steps$CurveGroup)
    step_cols <- cols[step_groups]
    for (i in seq_len(nrow(plot_data$steps))) {
      abline(
        v = plot_data$steps$Threshold[i],
        col = grDevices::adjustcolor(step_cols[i], alpha.f = 0.16),
        lty = 3
      )
    }
  }
  for (i in seq_along(groups)) {
    sub <- exp_df[exp_df$CurveGroup == groups[i], , drop = FALSE]
    lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2, lty = line_types[groups[i]])
  }
  visible_steps <- plot_data$steps[
    is.finite(plot_data$steps$Threshold) &
      is.finite(plot_data$steps$PathY) &
      plot_data$steps$Threshold >= min(exp_df$Theta) &
      plot_data$steps$Threshold <= max(exp_df$Theta) &
      plot_data$steps$PathY >= y_rng[1] & plot_data$steps$PathY <= y_rng[2], , drop = FALSE
  ]
  if (nrow(visible_steps) > 0) {
    step_groups <- as.character(visible_steps$CurveGroup)
    step_cols <- cols[step_groups]
    points(visible_steps$Threshold, visible_steps$PathY, pch = 18, col = step_cols)
    step_labels <- truncate_axis_label(visible_steps$ThresholdLabel, width = 6L)
    placed <- .place_plot_labels(visible_steps$Threshold,
      visible_steps$PathY + graphics::strheight("M", cex = 0.62) * 1.2,
      step_labels, cex = 0.62, point_x = visible_steps$Threshold,
      point_y = visible_steps$PathY, ylim = y_rng)
    segments(visible_steps$Threshold, visible_steps$PathY,
      pmax(placed$Left, pmin(placed$Right, visible_steps$Threshold)),
      pmax(placed$Bottom, pmin(placed$Top, visible_steps$PathY)),
      col = step_cols, lwd = 0.7)
    text(
      x = placed$X,
      y = placed$Y,
      labels = step_labels,
      cex = 0.62,
      col = "gray15"
    )
  }
  if (!is.null(plot_data$endpoint_labels) && nrow(plot_data$endpoint_labels) > 0 && length(groups) <= 5) {
    label_y <- .spread_wright_label_positions(
      plot_data$endpoint_labels$ExpectedScore, y_rng[1], y_rng[2],
      min_gap = 1.4 * graphics::strheight("M", cex = 0.78)
    )
    label_x <- max(exp_df$Theta) + diff(range(exp_df$Theta)) * 0.03
    segments(plot_data$endpoint_labels$Theta,
                       plot_data$endpoint_labels$ExpectedScore,
                       label_x, label_y, xpd = NA,
                       col = cols[as.character(plot_data$endpoint_labels$CurveGroup)])
    text(
      x = label_x,
      y = label_y,
      labels = truncate_axis_label(plot_data$endpoint_labels$CurveGroup, width = 14L),
      pos = 4,
      cex = 0.78,
      xpd = NA,
      col = "gray15"
    )
  } else {
    legend("topleft", legend = truncate_axis_label(groups),
                      col = cols[groups], lty = line_types[groups], lwd = 2, bty = "n")
  }
  if (length(basis_lines)) mtext(basis_lines, side = 1,
    line = note_start + 0.7 * seq_along(basis_lines), cex = 0.68)
}

build_fit_pathway_data <- function(x,
                                   diagnostics = NULL,
                                   fit_stat = c("Infit", "Outfit"),
                                   fit_scale = c("mnsq", "zstd"),
                                   zstd_method = c("facets", "engine"),
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
      fit_df_method = if (identical(fit_scale, "zstd")) "both" else "facets"
    )
  }
  if (!is.list(diag_use) || is.null(diag_use$measures)) {
    stop("`diagnostics` must contain a `measures` table.", call. = FALSE)
  }
  if (identical(fit_scale, "zstd")) {
    selected_z_col <- paste0(fit_stat, "ZSTD_", toupper(zstd_method))
    if (!selected_z_col %in% names(diag_use$measures)) {
      replay <- diag_use$replay_inputs %||% list()
      replay$residual_pca <- "none"
      replay$fit_df_method <- "both"
      diag_use <- do.call(diagnose_mfrm, c(list(fit = x), replay))
      diagnostics_source <- paste0(diagnostics_source, "; recomputed for engine/FACETS ZSTD companions")
    }
  }

  fm <- fit_measures_table(
    x,
    diagnostics = diag_use,
    include_person = isTRUE(include_person),
    fit_df_method = if (identical(fit_scale, "zstd")) "both" else "facets",
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
    stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE), 100 * ci_level
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
  default_cols <- .plot_series_colors(facets, plot_data$preset %||% "standard")
  cols <- resolve_palette(palette = palette, defaults = default_cols)
  refs <- as.data.frame(plot_data$reference_lines, stringsAsFactors = FALSE)
  x_ref <- suppressWarnings(as.numeric(refs$value))

  draw_panel <- function(panel_tbl, panel_title) {
    if (!(identical(plot_data$display$show_title, FALSE) &&
          identical(plot_data$panel, "facet") && length(unique(tbl$Panel)) > 1L)) {
      panel_title <- .mfrm_plot_display_title(
        plot_data, title = panel_title, fallback = "Fit pathway")
    }
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
    plot(
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
      rect(
        plot_data$fit_range[1], y_rng[1], plot_data$fit_range[2], y_rng[2],
        col = grDevices::adjustcolor("#E8F0FE", alpha.f = 0.55), border = NA
      )
    } else {
      rect(
        -plot_data$zstd_cut, y_rng[1], plot_data$zstd_cut, y_rng[2],
        col = grDevices::adjustcolor("#E8F0FE", alpha.f = 0.55), border = NA
      )
    }
    for (i in seq_len(nrow(refs))) {
      abline(
        v = refs$value[i],
        lty = if (refs$role[i] == "reference") 2 else 3,
        col = if (refs$role[i] == "reference") "gray35" else "gray65"
      )
    }
    abline(h = pretty(y_rng, n = 6), col = "gray92", lty = 1)
    point_cols <- unname(cols[as.character(panel_tbl$Facet)])
    if (isTRUE(show_ci)) {
      ci_ok <- is.finite(panel_tbl$CI_Lower) & is.finite(panel_tbl$CI_Upper)
      if (any(ci_ok)) {
        arrows(
          x0 = panel_tbl$FitValue[ci_ok], y0 = panel_tbl$CI_Lower[ci_ok],
          x1 = panel_tbl$FitValue[ci_ok], y1 = panel_tbl$CI_Upper[ci_ok],
          angle = 90, code = 3, length = 0.035,
          col = grDevices::adjustcolor(point_cols[ci_ok], alpha.f = 0.60)
        )
      }
    }
    points(
      panel_tbl$FitValue, panel_tbl$Measure,
      pch = panel_tbl$Shape,
      col = point_cols,
      bg = grDevices::adjustcolor(point_cols, alpha.f = 0.35),
      cex = ifelse(panel_tbl$ElementType == "Person", 0.82, 1)
    )
    label_rows <- nzchar(as.character(panel_tbl$LabelText))
    if (any(label_rows)) {
      text(
        panel_tbl$FitValue[label_rows], panel_tbl$Measure[label_rows],
        labels = truncate_axis_label(panel_tbl$LabelText[label_rows], width = 16L),
        pos = 4, cex = 0.65, offset = 0.35, xpd = NA,
        col = "gray15"
      )
    }
  }

  if (identical(plot_data$panel, "facet") && length(unique(tbl$Panel)) > 1L) {
    old_par <- graphics::par()[c("mfrow", "cex", "mex")]
    on.exit(graphics::par(old_par), add = TRUE)
    panel_names <- unique(as.character(tbl$Panel))
    n_cols <- ceiling(sqrt(length(panel_names)))
    graphics::par(mfrow = c(ceiling(length(panel_names) / n_cols), n_cols))
    for (nm in panel_names) {
      draw_panel(tbl[as.character(tbl$Panel) == nm, , drop = FALSE], nm)
    }
  } else {
    draw_panel(tbl, title %||% plot_data$title)
    legend(
      "topleft", legend = facets, col = cols[facets], pch = 16,
      bty = "n", cex = 0.76
    )
    if (any(tbl$ElementType == "Person")) {
      legend(
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
    probabilities = prob_df,
    curve_basis = .category_curve_basis(curve_spec)
  )
}

.category_curve_basis <- function(curve_spec) {
  data.frame(
      CurveBasis = curve_spec$curve_basis,
      PredictorOffset = curve_spec$predictor_offset,
      Description = "Fitted steps (GPCM: also slopes); additive effects and interactions fixed at zero.",
      stringsAsFactors = FALSE
  )
}

.ccc_reference_subtitle <- function(prefix = "Category response probabilities across theta") {
  paste0(
    prefix,
    "; reference profile fixes additive facet effects and fitted interactions at zero"
  )
}

.plot_series_colors <- function(categories, preset = "standard") {
  categories <- as.character(categories)
  colors <- if (identical(as.character(preset)[1L], "monochrome")) {
    grDevices::gray.colors(max(3L, length(categories)), start = 0.12, end = 0.55)
  } else if (length(categories) > 8L) {
    # ponytail: many series use the dark half of viridis; labels/panels remain
    # necessary when colour and the six available line types are insufficient.
    grDevices::hcl.colors(2L * length(categories), "viridis")[seq_along(categories)]
  } else {
    # Okabe-Ito foundation; omit bright yellow and darken pink/amber/cyan
    # entries for lines on white. Small text uses neutral ink, not these colours.
    c("#0072B2", "#D55E00", "#009E73", "#B56794",
      "#222222", "#666666", "#A66F00", "#007F9E")
  }
  stats::setNames(colors[seq_along(categories)], categories)
}

.plot_series_linetypes <- function(series) {
  stats::setNames(rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                      length.out = length(series)), as.character(series))
}

.ccc_plot_legend <- function(plot_core, style, palette = NULL, overlay = FALSE) {
  prob_df <- as.data.frame(plot_core$probabilities %||% data.frame(), stringsAsFactors = FALSE)
  categories <- unique(as.character(prob_df$Category %||% character(0)))
  defaults <- .plot_series_colors(categories, preset = style$name)
  colors <- resolve_palette(
    palette = palette,
    defaults = defaults
  )
  legend <- new_plot_legend(
    label = paste("Category", categories),
    role = rep("probability", length(categories)),
    aesthetic = rep("line", length(categories)),
    value = unname(colors[categories])
  )
  if (isTRUE(overlay)) {
    legend <- rbind(
      legend,
      new_plot_legend(
        label = "Observed bin proportion",
        role = "empirical",
        aesthetic = "point",
        value = "matching_category_color"
      )
    )
  }
  legend
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
  reporting_policy <- .mfrmr_table_build_ccc_surface_data_reporting_policy

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
  groups <- unique(as.character(prob_df$CurveGroup))
  categories <- unique(as.character(prob_df$Category))
  defaults <- .plot_series_colors(
    categories,
    preset = as.character(plot_data$preset %||% "standard")
  )
  cols <- resolve_palette(palette = palette, defaults = defaults)
  line_types <- .plot_series_linetypes(categories)
  point_shapes <- stats::setNames(rep(c(21, 22, 24, 23, 25),
    length.out = length(categories)), categories)
  display_title <- .mfrm_plot_display_title(
    plot_data,
    title = title,
    fallback = "Category characteristic curves"
  )
  curve_basis <- as.data.frame(plot_data$curve_basis %||% data.frame(), stringsAsFactors = FALSE)
  basis_note <- if (nrow(curve_basis) > 0L && "Description" %in% names(curve_basis)) {
    as.character(curve_basis$Description[1L])
  } else {
    "Reference profile: additive facet effects and fitted interactions are fixed at zero."
  }
  if (identical(plot_data$display$show_notes, FALSE)) basis_note <- ""
  overlay <- as.data.frame(plot_data$overlay %||% data.frame(), stringsAsFactors = FALSE)
  overlay_ok <- nrow(overlay) > 0L &&
    all(c("Theta", "Proportion", "Category") %in% names(overlay))
  shared_legend <- length(groups) > 1L || length(categories) > 5L

  draw_panel <- function(group, panel_title) {
    panel <- prob_df[as.character(prob_df$CurveGroup) == group, , drop = FALSE]
    plot(
      x = range(prob_df$Theta, finite = TRUE),
      y = c(0, 1),
      type = "n",
      xlab = "Theta / Logit",
      ylab = "Probability",
      main = panel_title
    )
    abline(v = 0, lty = 3, col = "gray75")
    for (category in categories) {
      sub <- panel[as.character(panel$Category) == category, , drop = FALSE]
      if (nrow(sub) == 0L) next
      lines(
        sub$Theta, sub$Probability,
        col = cols[category], lwd = 1.6, lty = line_types[category]
      )
    }
    if (overlay_ok) {
      max_n <- max(overlay$N, na.rm = TRUE)
      for (category in categories) {
        sub <- overlay[as.character(overlay$Category) == category, , drop = FALSE]
        if (nrow(sub) == 0L) next
        points(
          sub$Theta, sub$Proportion,
          pch = point_shapes[category], bg = cols[category], col = "gray20",
          cex = 0.72 + 0.5 * sqrt(sub$N / max(1, max_n))
        )
      }
    }
    if (!shared_legend) legend(
      "topright",
      legend = paste("Category", categories),
      col = cols[categories], lty = line_types[categories], lwd = 1.6,
      pch = if (overlay_ok) point_shapes[categories] else NA_integer_,
      pt.bg = cols[categories],
      bty = "n", cex = 0.72
    )
  }

  old_par <- graphics::par()[if (length(groups) > 1L) {
    c("mfrow", "cex", "mex", "mar", "oma")
  } else "mar"]
  on.exit(graphics::par(old_par), add = TRUE)
  legend_width <- if (shared_legend) {
    max(graphics::strwidth(paste("Category", categories), units = "inches", cex = 0.65)) + 0.5
  } else 0
  legend_margin <- legend_width / (graphics::par("csi") * graphics::par("mex"))
  note_width <- if (length(groups) > 1L) graphics::par("din")[1] else graphics::par("fin")[1]
  basis_lines <- strwrap(basis_note, width = max(30L, floor((note_width - legend_width - 0.8) * 16)))
  if (length(groups) > 1L) {
    n_cols <- ceiling(sqrt(length(groups)))
    n_rows <- ceiling(length(groups) / n_cols)
    graphics::par(
      mfrow = c(n_rows, n_cols),
      mar = c(4.8, 3.8, 3.2, 0.8),
      oma = c(if (length(basis_lines)) 1.4 + 0.7 * length(basis_lines) else 0.2,
              0.2, if (nzchar(display_title)) 3.1 else 0.2, max(0.2, legend_margin))
    )
    group_titles <- stats::setNames(truncate_axis_label(groups, width = 28L), groups)
    for (group in groups) {
      draw_panel(group, group_titles[group])
    }
    mtext(display_title, side = 3, outer = TRUE, line = 1.35, font = 2)
    if (length(basis_lines)) mtext(basis_lines, side = 1, outer = TRUE,
                    line = 0.4 + 0.7 * seq_along(basis_lines), cex = 0.68)
  } else {
    graphics::par(mar = c(if (length(basis_lines)) 4.7 + 0.7 * length(basis_lines) else 4.1,
                          4.1, if (nzchar(display_title)) 3.2 else 1,
                          if (shared_legend) legend_margin + 0.2 else 1.1))
    draw_panel(groups[1L], display_title)
    if (length(basis_lines)) mtext(basis_lines, side = 1,
                    line = 3.1 + 0.7 * seq_along(basis_lines), cex = 0.70)
  }
  if (shared_legend) legend(
    x = graphics::grconvertX(0.99, from = if (length(groups) > 1L) "ndc" else "nfc", to = "user"),
    y = graphics::grconvertY(0.86, from = if (length(groups) > 1L) "ndc" else "nfc", to = "user"),
    legend = paste("Category", categories), xjust = 1, yjust = 1,
    col = cols[categories], lty = line_types[categories], lwd = 1.6,
    pch = if (overlay_ok) point_shapes[categories] else NA_integer_,
    pt.bg = cols[categories],
    bty = "n", cex = 0.65, xpd = NA
  )
}

draw_person_plot <- function(person_tbl, bins, title = "Person measure distribution", palette = NULL) {
  pal <- resolve_palette(palette = palette, defaults = c(person_hist = "#2C7FB8"))
  hist(
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
  step_tbl <- calc_step_order(step_tbl)
  if (!nrow(step_tbl)) stop("No threshold estimates available.", call. = FALSE)
  estimates <- ifelse(is.finite(step_tbl$Estimate), step_tbl$Estimate, NA_real_)
  limits <- if (any(is.finite(estimates))) range(estimates, finite = TRUE) else c(-1, 1)
  labels <- if (isTRUE(all(step_tbl$StepFacet == "Common"))) as.character(step_tbl$Step) else
    paste(step_tbl$StepFacet, step_tbl$Step, sep = ": ")
  labels <- truncate_axis_label(labels, width = 24L)
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar <- old_mar
  angle <- label_angle * pi / 180
  label_height <- abs(sin(angle)) * graphics::strwidth(labels, units = "inches", cex = 0.84) +
    abs(cos(angle)) * graphics::strheight(labels, units = "inches", cex = 0.84)
  line_height <- graphics::par("csi") * graphics::par("mex")
  mar[1] <- max(mar[1], (max(label_height) + 0.24) / line_height)
  if (!length(title) || !nzchar(title)) mar[3] <- min(mar[3], 1)
  graphics::par(mar = mar)
  plot(
    x = seq_len(nrow(step_tbl)),
    y = estimates,
    ylim = limits,
    type = "p",
    pch = 16,
    xaxt = "n",
    xlab = "",
    ylab = "Estimate",
    main = title,
    col = pal["step_line"]
  )
  abline(h = pretty(limits, n = 5), col = pal["grid"], lty = 1)
  ends <- which(!is.na(step_tbl$Ordered))
  if (length(ends)) {
    segments(ends - 1L, estimates[ends - 1L], ends, estimates[ends], col = pal["step_line"])
  }
  if (anyNA(estimates)) {
    text(which(is.na(estimates)), mean(limits), "Unavailable", cex = 0.7)
  }
  draw_rotated_x_labels(
    at = seq_len(nrow(step_tbl)),
    labels = labels,
    srt = label_angle,
    cex = 0.84,
    line_offset = 0.12 / graphics::par("pin")[2]
  )
}

draw_facet_plot <- function(facet_tbl,
                            title,
                            palette = NULL,
                            show_ci = FALSE,
                            ci_level = 0.95) {
  pal <- resolve_palette(palette = palette, defaults = c(facet_bar = "gray70"))
  labels <- paste0(facet_tbl$Facet, ": ", facet_tbl$Level)
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  margin <- old_mar
  margin[2] <- max(margin[2], max(graphics::strwidth(labels, units = "inches",
    cex = graphics::par("cex.axis"), font = graphics::par("font.axis"))) /
    (graphics::par("csi") * graphics::par("mex")) + 1.2)
  graphics::par(mar = margin)
  limits <- range(c(0, facet_tbl$Estimate), finite = TRUE)
  if (isTRUE(show_ci) && "SE" %in% names(facet_tbl)) {
    z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
    limits <- range(c(limits, facet_tbl$Estimate - z * facet_tbl$SE,
      facet_tbl$Estimate + z * facet_tbl$SE), finite = TRUE)
  }
  mids <- barplot(
    height = facet_tbl$Estimate,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    main = title,
    xlab = "Estimate",
    xlim = limits,
    col = pal["facet_bar"],
    border = "white"
  )
  abline(v = 0, lty = 2, col = "gray40")
  if (isTRUE(show_ci) && "SE" %in% names(facet_tbl)) {
    z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
    ci_ok <- is.finite(facet_tbl$SE) & facet_tbl$SE > 0
    if (any(ci_ok)) {
      ci_lo <- facet_tbl$Estimate[ci_ok] - z * facet_tbl$SE[ci_ok]
      ci_hi <- facet_tbl$Estimate[ci_ok] + z * facet_tbl$SE[ci_ok]
      arrows(
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
#'   Additional diagnostic and reporting views are listed under "Single-fit
#'   views" in [mfrmr_plot_migration], with replacements for all former helpers.
#' @param facet Optional facet name for `type = "facet"`.
#' @param top_n Maximum number of non-person facet locations retained by the
#'   native Wright map for compact displays. Step transitions are always
#'   retained separately and do not consume this limit;
#'   any omitted facet locations are counted in the returned `retention` table
#'   and disclosed in the plot subtitle/note. Use `Inf` for a complete
#'   all-level final map. The FACETS-style payload always retains every fitted
#'   facet and step location.
#' @param highlight Optional named list of non-person facet level IDs, for
#'   example `list(Rater = c("R01", "R02"))`, for Wright maps only (including
#'   `view = "unified"` and the Wright component of a bundle). Selected labels
#'   receive a `> ` prefix in both renderers and saved ggplot conversion.
#'   Levels outside the usual `top_n` selection are added, so the displayed
#'   count may exceed that limit. Fitted locations and uncertainty do not
#'   change; the automatic scale uses all fitted locations before compaction.
#'   Selection is not a statistical flag. Inspect `locations$Highlighted`,
#'   `retention` and `notes` in the returned plot data. Unknown IDs or levels
#'   without finite locations are rejected. `NULL` selects no levels.
#' @param theta_range Numeric length-2 range for pathway, CCC, and
#'   category-surface plot data.
#' @param theta_points Number of theta grid points used for pathway, CCC, and
#'   category-surface plot data.
#' @param title Optional custom title.
#' @param palette Optional color overrides.
#' @param label_angle Rotation angle for x-axis labels where applicable.
#'   Native Wright maps reserve bottom space for full labels when abbreviating
#'   them would make different groups indistinguishable.
#' @param show_ci Whether to add approximate confidence intervals when
#'   available. `NULL` selects the display-specific default: `TRUE` for the
#'   primary native Wright map and `type = "fit_pathway"`, and `FALSE` for the
#'   FACETS-style ruler and other plot types. With `renderer = "facets"`,
#'   explicitly setting `show_ci = TRUE` creates a hybrid display: FACETS-style
#'   ruler grammar with mfrmr uncertainty intervals.
#' @param ci_level Confidence level used when `show_ci = TRUE`.
#' @param level Named alternative to `ci_level`, a single number strictly
#'   between 0 and 1. Supply only one of these names, even when equal.
#'   Both omitted retain 0.95. This controls the fitted-model plot intervals;
#'   attached interval results retain the level chosen when they were computed.
#' @param group Optional grouping for `type = "wright"` to overlay
#'   per-group person-density curves (DIF / DFF screening view).
#'   Either a column name (looked up first in `group_data` when
#'   supplied through `...`, then in `fit$prep$data`) or a vector
#'   aligned with `fit$facets$person`. Ignored for other `type`
#'   values and for `wright_style = "facets_style"`, whose person column is a
#'   single FACETS-style star frequency. To pass the source data alongside, use
#'   `plot(fit, type = "wright", group = "MyCol", group_data = <df>)`.
#' @param diagnostics Optional matching output from [diagnose_mfrm()].
#'   Recompute it after refitting; mismatched or outdated readiness records
#'   are rejected. When supplied,
#'   Wright-map standard errors and precision metadata reuse matching rows from
#'   `diagnostics$measures` without replacing fitted coordinates,
#'   while pathway plot data reuse `fit_measures`, `fit_status`, and
#'   `curve_fit_status` instead of recomputing diagnostics.
#'   Must match this fitted analysis; see [mfrmr_compatibility_layer] for
#'   migration of saved diagnostics.
#' @param include_fit_measures If `TRUE` (default), pathway plot data include
#'   tidy fit-measure and fit-status tables for custom R graphics. Set to
#'   `FALSE` when only the curve coordinates are needed.
#' @param fit_stat For `type = "fit_pathway"`, the horizontal fit statistic:
#'   `"Infit"` (default) or `"Outfit"`.
#' @param fit_scale For `type = "fit_pathway"`, use mean squares (`"mnsq"`)
#'   or standardized fit (`"zstd"`) on the horizontal axis.
#' @param zstd_method For a ZSTD fit pathway, `"facets"` (default) uses the
#'   FACETS/Wright-Masters fourth-moment df, matching [diagnose_mfrm()].
#'   `"engine"` retains the former package-native df convention. Neither
#'   establishes complete FACETS output equivalence. Refreshing companions
#'   from supplied diagnostics preserves their recorded WHEXACT choice.
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
#' @param show_title Logical; display the main plot title (default `FALSE`). Structural panel
#'   headings, axis labels, legends, and data labels remain visible.
#' @param show_notes Logical; display explanatory subtitles, footnotes, and
#'   review-only title markers. Notes remain available in the returned object,
#'   and R warnings are still issued. Default `FALSE`; use
#'   `plot_data(p, "notes")` for the descriptions and cautions outside the figure.
#' @param view Optional display variant for consolidated routes: for example,
#'   `"pairs"` or `"heatmap"` for agreement, `"unified"` for the former unified
#'   Wright helper, or `"funnel"` for shrinkage. See [mfrmr_plot_migration].
#' @param metric Quantity for a consolidated route, for example `"reliability"`,
#'   `"separation"` or `"strata"`, or `"exact"` / `"correlation"` for agreement.
#' @param comparison Another fitted model for paired Wright/CCC displays, or an
#'   ordered list of additional fits for `type = "trajectory"`. Plotting does
#'   not equate their scales or perform a model-selection test.
#' @param labels Labels for paired models or all trajectory waves, reference first.
#' @param ... View-specific controls listed in [mfrmr_plot_migration]. The
#'   consolidated routes reject irrelevant or duplicated arguments.
#'
#' @details
#' Start with `plot(fit)`: it draws a Wright map of person abilities, rater
#' severities, other facet locations, and category thresholds on one logit
#' scale. No plot options are needed for this first view.
#'
#' For adjusted-score comparisons, use `plot(fit, type = "fair_average",
#' show_ci = TRUE)`. `metric = "FairM"` uses mean reference measures;
#' `metric = "FairZ"` uses a zero reference, not a z-score. PCM averages
#' expectations over administered scales with the table's observation weights
#' and `totalscore` selection. Its conditional delta-method intervals use the
#' same retained observations as the table: finite direct anchors in ordinary
#' additive RSM/PCM JML are retained even with all-minimum/all-maximum responses,
#' including when `totalscore = FALSE`. The intervals use the
#' derivative of that same average. RSM/PCM intervals hold thresholds, other
#' effects and reference means fixed; they are diagnostic approximations,
#' not full calibration-uncertainty intervals. Inspect the returned data and
#' notes for the interval method and unavailable rows. See [fair_average_table()]
#' for reference definitions and [plot_fair_average] for view-specific controls.
#' For ordinary additive RSM/PCM JML, `xtreme = .3, xtreme_method = "facets"`
#' requests conditional total-score extrapolation of extreme measures and
#' their fair scores. Primary estimates and fixed anchors are retained.
#' Group-mean anchors are unsupported. The measure SE describes the corrected
#' total; fair-score intervals remain unavailable for positive `xtreme`.
#' The default `xtreme_method = "display"` retains the old average-score
#' display convention and leaves fair scores unchanged.
#'
#' This S3 plotting method provides the core fit-family visuals for
#' `mfrmr`. When `type` is omitted, it returns the Wright map alone as
#' an `mfrm_plot_data` object (the most useful single figure for a
#' first inspection). Pass `type = "bundle"` (or `"all"` / `"default"`)
#' to obtain the legacy three-plot `mfrm_plot_bundle` containing a
#' Wright map, pathway map, and category characteristic curves. The
#' compact native default records any omitted facet locations in
#' `data$retention` and annotates the subtitle and drawn figure; use
#' `top_n = Inf` to retain every fitted location in the final Wright map.
#' Every retained native location is labelled. `LabelX` / `LabelY` provide
#' initial text positions; base graphics further adjusts these positions for
#' the device's text dimensions, avoiding nearby labels and fitted points where
#' space permits. Leader lines connect displaced text to the unchanged point.
#' `label_points` retains fitted coordinates and initial text positions for
#' custom renderers. Step transitions share one vertical ladder and their
#' labels include the fitted threshold logit. When the fit records
#' boundary-separated facet levels and no
#' display range was supplied, the native and FACETS-style maps use the same
#' robust automatic range and place those levels at ruler ends. Exact fitted
#' values and intervals remain in `OriginalEstimate`, `CI_Lower`, and
#' `CI_Upper`; `DisplayEstimate`, `DisplayCI_Lower`, `DisplayCI_Upper`, and the
#' clipping-status columns describe only the rendered coordinates. Endpoint
#' triangles and the plot footer disclose any omitted or clipped interval.
#' The returned object always carries machine-readable metadata through
#' the `mfrm_plot_data` contract, even when the plot is drawn
#' immediately.
#' Set `show_title = FALSE` and `show_notes = FALSE` for a figure whose title
#' and explanation will be supplied by the surrounding document. The returned
#' `data$notes` table contains `Type` and `Text` columns for interpretation,
#' reference-profile conditions, and display/retention notes where applicable;
#' `print()` also prints these notes. `data$display` records the two flags.
#' Original titles, subtitles, coordinates, and readiness metadata are retained.
#' Device-dependent crowding warnings are issued during drawing rather than
#' stored in this draw-free notes table. `as_ggplot()` respects these flags and
#' retains the notes in `attr(plot, "mfrmr_notes")`.
#'
#' Fit-family Wright, pathway, and CCC displays use a shared, CUD-informed
#' palette. The first eight series use an Okabe-Ito foundation with bright
#' yellow omitted and darker pink/amber/cyan entries for light backgrounds.
#' Beyond eight series, colours come from the dark half of the viridis HCL
#' palette. Expected-score and category curves also vary line type; Wright
#' locations retain their distinct point shapes, and small data labels use
#' dark neutral text. Six line types cycle for larger series sets, so dense
#' plots still require labels, panels, or more space. This is not a guarantee
#' of perceptual separation for every viewer or device.
#' `preset = "monochrome"` applies to these series in both base and ggplot
#' renderers. Custom `palette` values are retained in `data$palette` and used
#' during ggplot conversion; colour overrides keep the non-colour encodings.
#' See \url{https://jfly.uni-koeln.de/color/} for Color Universal Design guidance.
#'
#' Every fit-derived payload also carries `data$fit_readiness`,
#' `data$interpretation_status`, and `data$interpretation_note`. Availability
#' and interpretability are separate: when a numerical, data, design, or
#' stability status requires review, the coordinates remain available for
#' diagnosis, but the call warns and marks the returned subtitle and drawn
#' title `REVIEW ONLY` by default. Hiding that annotation changes presentation
#' only; it does not change the fit's interpretation status or suppress warnings.
#' A prior-regularized extreme MML EAP remains in `fit$facets$person`. If its
#' source fit is blocked, however, it is omitted from the Wright-map scale so
#' that a finite but non-interpretable trace cannot collapse the display. The
#' exact excluded rows and reason remain in `data$person_exclusions`, and the
#' omission is stated in `data$retention_note`.
#'
#' `type = "unexpected"` uses the same default as
#' [unexpected_response_table()]: `abs_z_min = 3, rule = "residual"`.
#' The optional `rule = "either"` or `"both"` also uses `prob_max`.
#' Saved table plots retain the recorded thresholds. A fully evaluated screen
#' with no flags yields an empty display and retains its screening summary.
#'
#' `type = "wright"` shows persons, facet levels, and step thresholds on
#' a shared logit scale. Estimates are plotted as fitted, so the sign
#' convention follows the fit: higher Person values predict higher coded
#' scores in these positive-discrimination models. Whether that means more
#' comprehensible or less comprehensible speech depends on the rubric coding.
#' Under the default negative orientation, higher non-Person facet values
#' reduce expected coded scores; calling this greater severity/difficulty
#' presumes that higher scores represent better performance.
#' Facets listed in `fit_mfrm(positive_facets = ...)` are reversed
#' (higher values raise expected scores); state the active orientation in
#' figure captions when reporting.
#' The returned `notes`, also shown by `print()`, state the recorded facet
#' directions and the score-coding caveat. Missing orientation in an older or
#' imported result is not guessed. Category labels do not recode the data.
#' See [mfrmr_visual_diagnostics] for reverse-coded comprehensibility, fit,
#' conditional bias and Step explanations.
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
#' Column headings are retained even on narrow devices. If headings or person
#' frequency stars need more space, drawing warns; use a wider device or increase
#' `persons_per_star` for a more compact frequency column.
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
#' category response probabilities. Multiple curve groups are faceted rather
#' than overplotted by the native renderer, and category-specific legends use
#' the same colours across panels. Multiple groups or more than five categories use one legend
#' beside the plotting area; colour presets use distinct default colours
#' beyond eight categories. For `GPCM`, these curves retain the
#' estimated slope. With separate owners, each step-level/slope-level pair is
#' a labeled curve group. Such joint profiles have no single-facet Infit/Outfit
#' flag; inspect `fit_measures` or `type = "fit_pathway"` separately.
#' All curve families are reference-profile curves
#' with additive facet main effects and fitted interactions fixed at zero; the
#' native footer and ggplot subtitle disclose that conditioning.
#' Expected-score pathways use the same reference profile and expose the
#' same `curve_basis` table. These curves do not average over the observed
#' rater assignments or show a particular fitted interaction cell.
#' The draw-free `curve_basis`, `CurveBasis`, and `PredictorOffset` fields make
#' that conditioning explicit. `type = "ccc_surface"` or
#' `type = "category_surface"` returns 3D-ready category-probability surface
#' data for external rendering; it deliberately does not add a
#' plotly/rgl dependency or replace the 2D CCC/pathway reporting figures. The
#' returned object includes `category_support`, `interpretation_guide`, and
#' `reporting_policy` tables so retained zero-frequency categories and
#' manuscript-use boundaries remain visible to beginners. The remaining
#' types (`"facet"`, `"person"`, `"step"`, `"shrinkage"`) provide
#' compact location-specific displays.
#'
#' @section Facet locations and practical equivalence:
#' `plot(fit, type = "rater_severity", facet = "Rater")` orders locations
#' and shows available normal uncertainty intervals. Pass matching
#' `diagnostics` to reuse their precision metadata. The default has no
#' magnitude bands; `show_bands = TRUE` explicitly adds gray guides at
#' +/-0.5 and +/-1.0 logits. These distances from the fitted zero are not
#' equivalence bounds, evidence of homogeneity or calibration standards.
#' The returned `data` table uses numerical absolute-location ranges in
#' `Band`; a common range does not imply small pairwise differences.
#'
#' Use `p <- plot(fit, type = "rater_severity", diagnostics = diagnostics)`
#' and [plot_data()] to read `"data"` and `"notes"`. Intervals retain their
#' `CIUse`, `SE_Method` and `IntervalStatus`; individual whisker overlap is
#' not a contrast interval. [analyze_facet_equivalence()] separately examines
#' pairwise differences with their full covariance and a prespecified practical
#' bound in its supported MML scope. It does not use the profile's guides.
#' The zero depends on the model constraints/anchors, so read the signed
#' differences and their uncertainty instead of inferring quality from color.
#' Older saved profiles retain their original colors and labels when converted
#' with [as_ggplot()]; regenerate the figure from the saved fit and matching
#' diagnostics for the new display. No refitting is needed.
#'
#' @section Corrected JML plots:
#' A fit with an explicit `jml_correction_order` uses `type = "slopes"`
#' (default), `"locations"` or `"steps"`. Locations show one `facet`, defaulting
#' to the first declared facet. `style = "points"` shows estimates without
#' intervals; `"distribution"` shows their empirical cumulative distribution.
#' These displays do not classify rater quality. They support `title`, `caption`,
#' `show_title`, `show_notes`, `show_labels`, `sort`, `palette`, `preset`,
#' `text_scale`, `point_size`, and draw-free [as_ggplot()] conversion. Slopes
#' use a reference of one; location/step displays use zero. Ordinary plot types,
#' `show_ci = TRUE`, and unrelated plot controls are rejected for this estimator.
#' The figure retains point estimates if their RootSE is unavailable; it cannot
#' be drawn if no unambiguous point solution was obtained.
#'
#' @section Graphics layout:
#' Single-panel plots advance through a caller's `par(mfrow = ...)` or
#' `layout()` arrangement and restore the style and margins they change.
#' Native Wright maps, multi-group CCC plots, and faceted fit pathways create
#' their own page layouts; call these outside a custom `layout()` arrangement.
#' Use `renderer = "facets"` for a single-panel Wright map in a custom grid.
#' A 7 by 5 inch device is a useful starting size for standalone plots;
#' dense PCM labels benefit from a wider device.
#' Native Wright/pathway plots warn when labels cannot be separated within
#' the available space. Enlarge the device (for example, to 12 by 8 inches)
#' and inspect the result; this warning concerns readability, not model fit.
#' Abbreviated labels retain distinguishing text where possible, falling back
#' to full names if needed. Original names remain available in the plot data.
#' For non-Latin labels, select a graphics device and font that support the
#' characters before plotting; the package retains the caller's font family.
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
#'   data payload, plus a `notes` table and `display` settings. It also includes
#'   a one-row `scale_contract` table recording
#'   the fitted coordinate basis, population SD when applicable,
#'   discrimination basis, and GPCM MML identification convention.
#' @seealso [fit_mfrm()], [plot(x, type = "wright", view = "unified")][plot.mfrm_fit], [plot(x, type = "bubble")][plot.mfrm_fit],
#'   [mfrmr_visual_diagnostics]
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept shrinkage
#' @examples
#' library(mfrmr)
#'
#' # Saved synthetic MML/RSM analysis of example_operational.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "operational-mml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#'
#' # Run each plot command separately to inspect its figure
#' plot(fit) # Wright map: persons, facets, and category thresholds
#'
#' # Rater severity estimates (higher means stricter in this example)
#' plot(fit, type = "facet", facet = "Rater")
#'
#' # Probability of each score category across the ability scale
#' plot(fit, type = "ccc")
#'
#' # Optional: get plot data instead of drawing a figure
#' wright <- plot(fit, draw = FALSE)
#' head(wright$data$locations)
#' # Find a named judge without refitting or dropping the surrounding scale
#' judge <- subset(wright$data$locations, Group == "Rater")$Label[1]
#' focused <- plot(fit, highlight = list(Rater = judge), draw = FALSE)
#' subset(plot_data(focused, "locations"), Highlighted)
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
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
                          zstd_method = c("facets", "engine"),
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
                          show_title = FALSE,
                          show_notes = FALSE,
                          ..., level = NULL, view = NULL, metric = NULL,
                          comparison = NULL, labels = NULL, highlight = NULL) {
  if (mfrm_has_jml_adjustment(x)) {
    supplied <- names(as.list(match.call())[-1L])
    allowed <- c("x","type","facet","draw","title","show_title","show_notes","palette","preset",
      "show_ci","style","sort","caption","show_labels","text_scale","point_size")
    if (length(setdiff(supplied,allowed))) stop("Unsupported corrected-JML plot arguments: ",
      paste(setdiff(supplied,allowed),collapse=", "),call.=FALSE)
    if (!is.null(show_ci) && !identical(show_ci,FALSE))
      stop("Corrected-JML plots show point estimates; structural confidence intervals are not available.",call.=FALSE)
    if(missing(preset)) preset <- .mfrm_default_plot_preset()
    preset <- resolve_plot_preset(preset)$name
    dots <- list(...)
    if(is.null(dots$text_scale)) dots$text_scale <- switch(preset,compact=.9,publication=1.05,1)
    return(do.call(mfrm_jml_plot,c(list(x=x,type=type,facet=facet,draw=draw,title=title,
      show_title=show_title,show_notes=show_notes,palette=palette %||% if(preset=="monochrome") "mono" else "accessible"),dots)))
  }
  stop_if_product_slopes(x, "plot.mfrm_fit()")
  if (!is.null(diagnostics)) {
    mfrm_results_validate_diagnostics_identity(x, diagnostics, helper = "plot.mfrm_fit()")
  }
  selected <- if (is.null(type)) "wright" else .mfrm_plot_selector(type, "type")
  selected <- tolower(selected)
  if (!is.null(highlight) && !selected %in% c("wright", "bundle", "all", "default")) {
    stop("`highlight` is supported only for Wright maps (including the bundle).", call. = FALSE)
  }
  route <- .mfrm_fit_plot_route(selected, view)
  if (!is.na(route) || !is.null(comparison) || selected == "trajectory") {
    if (!inherits(x, "mfrm_fit") || inherits(x, "mfrm_imported_fit"))
      stop("This plot route requires a native, supported mfrm_fit object.", call. = FALSE)
    supplied <- names(as.list(match.call(expand.dots = FALSE)))[-1L]
    supplied <- setdiff(supplied, c("x", "type", "..."))
    args <- c(mget(supplied, envir = environment(), inherits = FALSE), list(...))
    args <- .mfrm_plot_rename(args, "level", "ci_level")
    return(.mfrm_dispatch_fit_plot(x, selected, args))
  }
  if (!is.null(view) || !is.null(metric) || !is.null(labels))
    stop("`view`, `metric` or `labels` is not supported for this plot type. See ?mfrmr_plot_migration.", call. = FALSE)
  if (!missing(level)) {
    if (!missing(ci_level)) stop("Supply only one of `level` and `ci_level`.", call. = FALSE)
    if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
        is.na(level) || !is.finite(level) || level <= 0 || level >= 1) {
      stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
    }
    ci_level <- level
  }
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(x, "mfrm_fit")) {
    stop("`x` must be an mfrm_fit object from fit_mfrm().")
  }
  if (inherits(x, "mfrm_imported_fit")) {
    .validate_imported_metric(x)
    if (!is.null(type) && !tolower(type[1L]) %in% "wright") {
      stop("Imported fits support point-only Wright maps. Inspect source fit statistics in the imported tables; use the source package for model-dependent displays.", call. = FALSE)
    }
    if (isTRUE(show_ci)) {
      stop("Source uncertainty estimates do not support a shared confidence-interval calculation. Inspect source-package intervals; imported Wright maps show point estimates only.", call. = FALSE)
    }
    show_ci <- FALSE
  }
  for (arg in c("show_title", "show_notes")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
    }
  }
  show_ci_auto <- (missing(show_ci) || is.null(show_ci)) &&
    !inherits(x, "mfrm_imported_fit")
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
  wright_colors <- .plot_series_colors(c("facet_level", "step_threshold"), style$name)
  style$accent_tertiary <- unname(wright_colors[1L])
  style$accent_secondary <- unname(wright_colors[2L])
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
    data <- data[!duplicated_named]
    data$scale_contract <- mfrm_fit_scale_contract(x)
    if (inherits(x, "mfrm_imported_fit")) {
      data$axis_label <- "Source ability scale"
      data$subtitle <- "Source person and item/design-cell estimates; point display only"
      data$reference_lines <- new_reference_lines("h", 0, "Source scale zero", "dashed", "reference")
      data$uncertainty_display <- "source_points_only"
    }
    out <- new_mfrm_plot_data(name, data)
    out <- .mfrm_attach_plot_readiness(out, fit_plot_readiness)
    out$data$display <- list(show_title = show_title, show_notes = show_notes)
    out$data$palette <- palette
    out$data$notes <- .mfrm_plot_notes(out$data)
    out$data$notes <- rbind(out$data$notes, .mfrm_score_direction_notes(x))
    if (name %in% c("wright_map", "pathway_map", "step",
                    "category_characteristic_curves", "category_characteristic_curves_overlay",
                    "category_probability_surface")) {
      out$data$notes <- rbind(out$data$notes,
        data.frame(Type = "step", Text = .mfrm_step_interpretation_note()))
    }
    if (name %in% c("fit_pathway", "pathway_map")) {
      out$data$notes <- rbind(out$data$notes, data.frame(Type = "fit", Text = c(
        "Infit/Outfit mean squares compare residual variation with model expectation (1): lower values indicate overfit, higher values underfit. Flagging uses the selected review limits, not every departure from 1.",
        "Similar judge severities or high agreement after calibration do not by themselves imply overfit. Review residual patterns, category use and rating independence; a fit flag is not a cause or an automatic exclusion rule.")))
    }
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
    legend_style <- style
    colors <- resolve_palette(palette, c(facet_level = style$accent_tertiary,
      step_threshold = style$accent_secondary, person_hist = style$fill_muted))
    legend_style$accent_tertiary <- unname(colors["facet_level"])
    legend_style$accent_secondary <- unname(colors["step_threshold"])
    legend_style$fill_muted <- unname(colors["person_hist"])
    .wright_map_legend(
      plot_core,
      style = legend_style,
      show_ci = show_ci,
      ci_level = ci_level
    )
  }

  se_tbl_ci <- if (!inherits(x, "mfrm_imported_fit") &&
                   (isTRUE(show_ci) || !is.null(diagnostics))) {
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
      highlight = highlight,
      se_tbl = se_tbl_ci,
      wright_style = wright_style,
      category_labels = category_labels,
      rows_per_logit = rows_per_logit,
      wright_range = wright_range,
      extreme_placement = extreme_placement,
      persons_per_star = persons_per_star
    )
    ccc_core <- build_ccc_data(
      x,
      theta_range = theta_range,
      theta_points = theta_points
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
          subtitle = .ccc_reference_subtitle("Expected score and dominant categories across theta"),
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
        ccc_core,
        list(
          title = title %||% "Category characteristic curves",
          subtitle = .ccc_reference_subtitle(),
          preset = style$name,
          legend = .ccc_plot_legend(ccc_core, style = style, palette = palette),
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
        palette = palette
      )
      draw_ccc(
        out$category_characteristic_curves$data,
        title = title,
        palette = palette
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
                      "Shrinkage direction", "Zero shrinkage target")
    legend_role <- c("location", "location", "arrow", "reference")
    legend_aesthetic <- c("point", "point", "arrow", "line")
    legend_value <- c(style$accent_primary, style$accent_tertiary,
                      style$neutral, style$neutral)
    if (isTRUE(show_ci)) {
      legend_label <- c(legend_label, sprintf("%g%% plug-in band",
                                              round(100 * ci_level)))
      legend_role <- c(legend_role, "interval")
      legend_aesthetic <- c(legend_aesthetic, "line")
      legend_value <- c(legend_value, style$accent_primary)
    }
    out <- as_plot_data("shrinkage", list(
      data = data_list$table,
      shrinkage_report = data_list$report,
      mode = data_list$mode,
      interpretation = shrinkage_interpretation(),
      show_ci = isTRUE(show_ci),
      ci_level = if (isTRUE(show_ci)) ci_level else NA_real_,
      title = title %||% "Empirical-Bayes shrinkage",
      subtitle = "Original and shrunk estimates; conditional adjustment only",
      preset = style$name,
      legend = new_plot_legend(
        label = legend_label,
        role = legend_role,
        aesthetic = legend_aesthetic,
        value = legend_value
      ),
      reference_lines = new_reference_lines(
        "v", 0, "Zero shrinkage target", "dashed", "reference"
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
      highlight = highlight,
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
        subtitle = .ccc_reference_subtitle("Expected score and dominant categories across theta"),
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
      draw_pathway_map(out$data, title = title, palette = palette)
    }
    return(invisible(out))
  }
  if (type == "ccc") {
    ccc_core <- build_ccc_data(
      x,
      theta_range = theta_range,
      theta_points = theta_points
    )
    out <- as_plot_data("category_characteristic_curves", c(
      ccc_core,
      list(
        title = title %||% "Category characteristic curves",
        subtitle = .ccc_reference_subtitle(),
        preset = style$name,
        legend = .ccc_plot_legend(ccc_core, style = style, palette = palette),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_ccc(out$data, title = title, palette = palette)
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
        subtitle = .ccc_reference_subtitle(sprintf(
          "Model curves vs observed proportions in %d theta bin(s)",
          length(unique(overlay$Bin))
        )),
        preset = style$name,
        legend = .ccc_plot_legend(
          ccc_payload,
          style = style,
          palette = palette,
          overlay = TRUE
        ),
        reference_lines = new_reference_lines("v", 0, "Centered theta reference",
                                               "dashed", "reference")
      )
    ))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      draw_ccc(out$data, title = title, palette = palette)
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
    se_join <- se_tbl_ci[, intersect(
      c(
        "Facet", "Level", "SE", "SE_Method", "PrecisionTier",
        "SupportsFormalInference", "SEUse", "CIBasis", "CIUse",
        "CIEligible", "CILabel", "Measure_Source"
      ),
      names(se_tbl_ci)
    ), drop = FALSE]
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
