.mfrm_docx_table <- function(doc, heading, table, digits) {
  if (NROW(table) == 0L) return(doc)
  doc <- officer::body_add_par(doc, heading, style = "heading 2")
  table <- .tidymfrm_round_df(table, digits)
  if (requireNamespace("flextable", quietly = TRUE)) {
    ft <- flextable::autofit(flextable::theme_booktabs(flextable::flextable(table)))
    return(flextable::body_add_flextable(doc, ft))
  }
  officer::body_add_table(doc, value = table, style = NULL)
}

.mfrm_rater_card <- function(wf, rater, criterion_facet, max_cases) {
  if (!is.character(rater) || length(rater) != 1L || is.na(rater) || !nzchar(rater))
    stop("rater must be one nonempty rater ID.", call. = FALSE)
  sp <- wf$spec
  rt <- wf$diagnostics$raters
  selected <- match(rater, rt$Level)
  if (is.na(selected)) stop("Run step_mfrm_diagnose() and supply a known rater ID.", call. = FALSE)
  rx <- rt[selected, , drop = FALSE]
  ready <- .mfrm_readiness(wf$fit)$ready
  rating_rows <- wf$diagnostics$observations
  obs <- rating_rows[as.character(rating_rows[[sp$rater_facet]]) == rater, , drop = FALSE]
  cmp <- as.character(rt$Level[is.finite(rt$Estimate)])
  severity_rank <- if (ready && isTRUE(wf$connectivity$rater_linkage$n_components == 1L) && rater %in% cmp)
    1 + sum(rt$Estimate[rt$Level %in% cmp] < rx$Estimate)
  severity <- data.frame(SeverityLogit = rx$Estimate, ModelSE = rx$ModelSE, RealSE = rx$RealSE,
                         RankFromLeastSevere = severity_rank %||% NA_real_,
                         RatersInRanking = if (length(severity_rank)) length(cmp) else NA_integer_)
  audit <- data.frame(RatingRows = nrow(obs), Persons = length(unique(obs[[sp$person]])),
                      ObservedAverage = rx$ObsAvg, FairM = rx$FairM,
                      Infit = rx$InfitMnSq, Outfit = rx$OutfitMnSq)
  notes <- c(paste("Rater:", rater, "|", rx$DiagnosticFlag),
             paste("Recommended review:", rx$RecommendedAction),
             "Higher severity implies lower expected scores when other fitted effects are fixed. The comparison panel uses the observed rating network.")

  correlation <- function(a, b, persons) {
    ok <- is.finite(a) & is.finite(b)
    if (ready && length(unique(persons[ok])) >= 3L &&
        length(unique(a[ok])) > 1L && length(unique(b[ok])) > 1L)
      stats::cor(a[ok], b[ok]) else NA_real_
  }
  possible <- setdiff(sp$facets, sp$rater_facet)
  if (!is.null(criterion_facet) && (length(criterion_facet) != 1L || !criterion_facet %in% possible))
    stop("criterion_facet must select a non-rater facet.", call. = FALSE)
  if (is.null(criterion_facet))
    criterion_facet <- if (length(possible) == 1L) possible else if ("Criterion" %in% possible) "Criterion"
  center_rows <- if (is.null(criterion_facet)) list(All = seq_len(nrow(obs))) else
    split(seq_len(nrow(obs)), as.character(obs[[criterion_facet]]))
  center <- do.call(rbind, lapply(names(center_rows), function(level) {
    z <- obs[center_rows[[level]], , drop = FALSE]
    usable <- is.finite(z$Residual) & is.finite(z$Expected)
    constant <- length(unique(z$Observed[usable])) < 2L
    value <- if (constant) NA_real_ else correlation(z$Residual, z$Expected, z[[sp$person]])
    data.frame(Criterion = level, ResidualExpectedCorrelation = value,
      Persons = length(unique(z[[sp$person]][usable])), RatingRows = sum(usable),
      Status = if (constant) "no_observed_score_variation" else
        if (is.finite(value)) "descriptive" else "unavailable")
  }))
  rownames(center) <- NULL
  halo <- NULL
  if (!is.null(criterion_facet)) {
    criterion <- as.character(obs[[criterion_facet]])
    contexts <- c(sp$person, setdiff(sp$facets, c(sp$rater_facet, criterion_facet)))
    if (!is.null(sp$occasion) && !sp$occasion %in% contexts) {
      obs[[sp$occasion]] <- wf$data[[sp$occasion]][obs$.row_id]
      contexts <- c(contexts, sp$occasion)
    }
    rows <- which(stats::complete.cases(obs[c(contexts, criterion_facet)]))
    key <- .mfrm_mvgt_group(lapply(obs[rows, contexts, drop = FALSE], factor))
    obs <- obs[rows, , drop = FALSE]
    criterion <- criterion[rows]
    if (length(unique(criterion)) >= 2L) {
      pairs <- utils::combn(unique(criterion), 2L)
      halo <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
        ia <- which(criterion == pairs[1L, i])
        ib <- which(criterion == pairs[2L, i])
        ka <- key[ia]; kb <- key[ib]
        shared <- intersect(ka, kb)
        ambiguous <- anyDuplicated(ka[ka %in% shared]) || anyDuplicated(kb[kb %in% shared])
        a <- ia[match(shared, ka)]; b <- ib[match(shared, kb)]
        ok <- is.finite(obs$StdResidual[a]) & is.finite(obs$StdResidual[b])
        paired_constant <- length(unique(obs$Observed[a][ok])) < 2L ||
          length(unique(obs$Observed[b][ok])) < 2L
        value <- if (!ambiguous && !paired_constant)
          correlation(obs$StdResidual[a], obs$StdResidual[b], obs[[sp$person]][a]) else NA_real_
        data.frame(CriterionA = pairs[1L, i], CriterionB = pairs[2L, i],
                   MatchedPersons = length(unique(obs[[sp$person]][a][ok])),
                   ResidualCorrelation = value,
                   Status = if (ambiguous) "ambiguous_repeated_context" else
                     if (!any(ok)) "no_matched_scores" else
                     if (paired_constant) "no_observed_score_variation" else
                     if (is.finite(value)) "descriptive" else "unavailable")
      }))
    }
  }
  notes <- c(notes,
    "Residual versus expected-score correlations are descriptive within each criterion. Negative values may reflect score compression, model misfit or conditioning on estimated person scores; they do not establish central tendency.",
    "Criterion-pair residual correlations are descriptive. Their null baseline depends on the fitted model, person scoring and rating design and can be negative. Neither a zero cutoff nor -1/(K-1) is calibrated here. Halo claims require a design-matched reference and review of shared performance features and model misfit.",
    if (is.null(criterion_facet)) "No criterion facet was selected; the central-tendency table pools scores and may confound criterion difficulty with compression.",
    if (is.null(halo)) "Matched criterion pairs are unavailable for a halo screen.")
  cases <- wf$diagnostics$unexpected
  cases <- cases[as.character(cases[[sp$rater_facet]]) == rater, , drop = FALSE]
  case_cols <- intersect(c(".row_id", possible, "Observed", "Expected", "StdResidual"), names(cases))
  cases <- utils::head(cases[, case_cols, drop = FALSE], max_cases)
  names(cases)[names(cases) == ".row_id"] <- "InputRow"
  list(notes = notes, tables = list("Severity and comparison panel" = severity,
       "Rater fit and observed rubric scores" = audit, "Central tendency screen" = center,
       "Criterion halo screen" = halo, "Largest unexpected ratings" = cases))
}

.mfrm_render_rater_card <- function(wf, rater, criterion_facet, max_cases, format, file, title, digits) {
  card <- .mfrm_rater_card(wf, rater, criterion_facet, max_cases)
  if (identical(format, "docx")) {
    if (!requireNamespace("officer", quietly = TRUE)) stop("Install officer to export Word.", call. = FALSE)
    if (is.null(file) || !nzchar(file)) stop("Supply file for Word output.", call. = FALSE)
    doc <- officer::body_add_par(officer::read_docx(), title, style = "heading 1")
    for (line in card$notes) doc <- officer::body_add_par(doc, line)
    for (heading in names(card$tables)) doc <- .mfrm_docx_table(doc, heading, card$tables[[heading]], digits)
    print(doc, target = file)
    return(invisible(file))
  }
  lines <- c(paste("#", title), "", paste0("- ", card$notes), "")
  for (heading in names(card$tables)) {
    table <- card$tables[[heading]]
    if (NROW(table)) lines <- c(lines, paste("##", heading), "", .tidymfrm_fmt_df_md(table, digits), "")
  }
  report <- paste(lines, collapse = "\n")
  if (identical(format, "html")) report <- .tidymfrm_md_to_html(report, title)
  if (!is.null(file)) writeLines(report, file, useBytes = TRUE)
  report
}
