# ==============================================================================
# Reproducible APA / JARS-Quant & Individual Rater Feedback Report Generator
# Optional `flextable`, `kableExtra`, and `officer` exporters in `Suggests`
# ==============================================================================

.tidymfrm_round_df <- function(df, digits = 3L) {
  d <- as.data.frame(df, stringsAsFactors = FALSE)
  is_num <- vapply(d, function(col) is.numeric(col) && !is.integer(col), logical(1))
  d[is_num] <- lapply(d[is_num], round, digits = digits)
  d
}

.tidymfrm_fmt_df_md <- function(df, digits = 3L) {
  if (is.null(df) || nrow(df) == 0L) return("No rows to display.")
  d <- as.data.frame(df, stringsAsFactors = FALSE)
  fmt_col <- function(v) {
    if (is.logical(v)) return(ifelse(is.na(v), "NA", ifelse(v, "Yes", "No")))
    if (is.integer(v)) return(ifelse(is.na(v), "NA", as.character(v)))
    if (is.numeric(v)) {
      return(ifelse(
        is.na(v), "NA",
        ifelse(!is.finite(v), as.character(v), sprintf(paste0("%.", digits, "f"), v))
      ))
    }
    as.character(ifelse(is.na(v), "NA", v))
  }
  escape_cell <- function(x) vapply(x, function(value) {
    value <- gsub("[\r\n]+", " ", value)
    chars <- strsplit(value, "", fixed = TRUE)[[1L]]
    paste(ifelse(chars %in% c("\\", "|"), paste0("\\", chars), chars), collapse = "")
  }, character(1), USE.NAMES = FALSE)
  cols <- lapply(d, function(x) escape_cell(fmt_col(x)))
  hdr <- paste0("| ", paste(escape_cell(names(d)), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |")
  rows <- vapply(seq_len(nrow(d)), function(i) {
    paste0("| ", paste(vapply(cols, `[`, character(1), i), collapse = " | "), " |")
  }, character(1))
  paste(c(hdr, sep, rows), collapse = "\n")
}

.tidymfrm_escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

.tidymfrm_md_table_cells <- function(line) {
  # An even run of backslashes leaves the following pipe as a delimiter.
  inside <- substr(line, 2L, nchar(line) - 1L)
  cells <- strsplit(inside, '(?<!\\\\)(?:\\\\\\\\)*\\K\\|', perl = TRUE)[[1L]]
  trimws(gsub('\\\\([\\\\|])', '\\1', cells, perl = TRUE))
}

.tidymfrm_md_to_html <- function(md_text, title) {
  lines <- strsplit(md_text, "\n", fixed = TRUE)[[1L]]
  out <- character(0)
  in_table <- FALSE
  in_list <- FALSE

  flush_list <- function() {
    if (in_list) {
      out <<- c(out, "</ul>")
      in_list <<- FALSE
    }
  }
  flush_table <- function() {
    if (in_table) {
      out <<- c(out, "</tbody></table>")
      in_table <<- FALSE
    }
  }

  inline_fmt <- function(s) {
    s <- .tidymfrm_escape_html(s)
    s <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", s)
    s <- gsub("`([^`]+)`", "<code>\\1</code>", s)
    s
  }

  i <- 1L
  n <- length(lines)
  while (i <= n) {
    ln <- lines[i]
    if (grepl("^\\|", ln)) {
      flush_list()
      cells <- .tidymfrm_md_table_cells(ln)
      if (!in_table) {
        out <- c(out, "<table><thead><tr>", paste0("<th>", vapply(cells, .tidymfrm_escape_html, character(1)), "</th>", collapse = ""), "</tr></thead><tbody>")
        in_table <- TRUE
        if (i + 1L <= n && grepl("^\\|\\s*---", lines[i + 1L])) i <- i + 1L
      } else {
        out <- c(out, "<tr>", paste0("<td>", vapply(cells, .tidymfrm_escape_html, character(1)), "</td>", collapse = ""), "</tr>")
      }
    } else {
      flush_table()
      if (grepl("^# ", ln)) {
        flush_list()
        out <- c(out, paste0("<h1>", inline_fmt(sub("^# ", "", ln)), "</h1>"))
      } else if (grepl("^## ", ln)) {
        flush_list()
        out <- c(out, paste0("<h2>", inline_fmt(sub("^## ", "", ln)), "</h2>"))
      } else if (grepl("^### ", ln)) {
        flush_list()
        out <- c(out, paste0("<h3>", inline_fmt(sub("^### ", "", ln)), "</h3>"))
      } else if (grepl("^- ", ln)) {
        if (!in_list) {
          out <- c(out, "<ul>")
          in_list <- TRUE
        }
        out <- c(out, paste0("<li>", inline_fmt(sub("^- ", "", ln)), "</li>"))
      } else if (grepl("^> ", ln)) {
        flush_list()
        out <- c(out, paste0("<blockquote>", inline_fmt(sub("^> ", "", ln)), "</blockquote>"))
      } else if (nzchar(trimws(ln))) {
        flush_list()
        out <- c(out, paste0("<p>", inline_fmt(ln), "</p>"))
      }
    }
    i <- i + 1L
  }
  flush_list()
  flush_table()

  css <- paste(
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;max-width:960px;margin:2rem auto;padding:0 1.2rem;color:#1e293b;line-height:1.55;background:#ffffff;}",
    "h1,h2,h3{color:#0f172a;margin-top:1.6rem;}",
    "h1{border-bottom:2px solid #cbd5e1;padding-bottom:0.4rem;}",
    "h2{border-bottom:1px solid #e2e8f0;padding-bottom:0.25rem;}",
    "table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:0.92rem;}",
    "th,td{border:1px solid #cbd5e1;padding:0.42rem 0.6rem;text-align:left;}",
    "th{background:#f8fafc;font-weight:600;}",
    "code{background:#f1f5f9;padding:0.1rem 0.3rem;border-radius:3px;font-size:0.88em;}",
    "blockquote{border-left:4px solid #0072B2;background:#f8fafc;margin:1rem 0;padding:0.6rem 1rem;color:#334155;}",
    sep = "\n"
  )

  paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<title>",
    .tidymfrm_escape_html(title),
    "</title>\n<style>\n", css, "\n</style>\n</head>\n<body>\n",
    paste(out, collapse = "\n"),
    "\n</body>\n</html>\n"
  )
}
























format_mfrm_table <- function(x,
                              component = c("raters", "facets", "persons", "steps", "slopes",
                                            "g_study", "d_study", "dff", "drift"),
                              engine = c("markdown", "kableExtra", "flextable"),
                              digits = 3L) {
  component <- match.arg(component)
  engine <- match.arg(engine)
  df <- .tidymfrm_round_df(tidy(x, component = component), digits = digits)
  if (identical(engine, "markdown")) {
    return(.tidymfrm_fmt_df_md(df, digits = digits))
  }
  if (identical(engine, "kableExtra")) {
    if (!requireNamespace("knitr", quietly = TRUE) || !requireNamespace("kableExtra", quietly = TRUE)) {
      stop("Install suggested packages 'knitr' and 'kableExtra' to use `engine = 'kableExtra'`.", call. = FALSE)
    }
    kb <- knitr::kable(df, format = "html", digits = digits)
    return(kableExtra::kable_styling(kb, full_width = FALSE))
  }
  if (!requireNamespace("flextable", quietly = TRUE)) {
    stop("Install suggested package 'flextable' to use `engine = 'flextable'`.", call. = FALSE)
  }
  ft <- flextable::flextable(df)
  flextable::autofit(flextable::theme_booktabs(ft))
}

.tidymfrm_plan_overview <- function(gt, digits) {
  if (is.null(gt)) return(list(notes = paste(
    "No G-Study result is attached. Check whether the assignments identify the",
    "requested variance components before comparing D-Study plans."
  ), table = NULL))

  ds <- gt$d_study
  budgets <- numeric()
  notes <- paste("Use G to compare relative standing and Phi for decisions",
                 "against a fixed standard on the observed rubric scale.")
  if (length(gt$nesting)) {
    nested_counts <- paste(names(gt$nesting), "per", gt$nesting)
    notes <- c(notes, paste0("Read nested counts as ", paste(nested_counts, collapse = ", "), "."))
  }
  for (metric in c("Phi", "G")) {
    target <- gt$targets[[metric]]
    n_available <- sum(is.finite(ds[[metric]]))
    chosen <- gt[[paste0("recommended_for_", tolower(metric))]]
    if (!n_available) {
      why <- paste(unique(ds[[paste0(metric, "Status")]]), collapse = ", ")
      note <- sprintf("All %d candidate plans have unavailable %s (%s). Attainment of the %.2f target cannot be assessed.",
                      nrow(ds), metric, why, target)
    } else if (!nrow(chosen)) {
      note <- sprintf("None of the %d evaluable plans reaches %s = %.2f (%d plans unavailable).",
                      n_available, metric, target, nrow(ds) - n_available)
    } else {
      allocation <- paste(paste(gt$facets, unlist(chosen[1L, gt$facets], use.names = FALSE)), collapse = " x ")
      budgets <- c(budgets, chosen$TotalRatingsPerPerson[1L])
      note <- sprintf(paste("For the %s target of %.2f, the lowest rating count among",
        "evaluable plans that reach the target is %.0f per person (%s, projected %s = %.*f)."),
        metric, target, chosen$TotalRatingsPerPerson[1L], allocation, metric, digits, chosen[[metric]][1L])
    }
    notes <- c(notes, note)
  }
  equal_count <- ds[ds$TotalRatingsPerPerson %in% budgets,
    c(gt$facets, "TotalRatingsPerPerson", "G", "Phi", "MeetsTargetG", "MeetsTargetPhi"), drop = FALSE]
  equal_count <- equal_count[order(equal_count$TotalRatingsPerPerson), , drop = FALSE]
  if (length(budgets)) notes <- c(notes, paste("These are point projections for the evaluated",
    "designs. Check staffing costs and the availability of comparable raters or criteria."))
  list(notes = notes, table = equal_count)
}














































render_mfrm_report <- function(x,
                               format = c("markdown", "html", "docx"),
                               file = NULL,
                               title = "Many-Facet Measurement, Generalizability, and Rater Audit Report",
                               include_rater_cards = TRUE,
                               max_unexpected = 10L,
                               digits = 3L,
                               rater = NULL,
                               criterion_facet = NULL,
                               reference_plan = NULL,
                               fair_level = 0.95) {
  if (!inherits(x, "tidymfrm")) {
    stop("`x` must be a `tidymfrm` object created by `mfrm_analysis()` or `mfrm_workflow()`.", call. = FALSE)
  }
  x <- .tidymfrm_refresh(x)
  format <- match.arg(format)
  if (length(max_unexpected) != 1L || !is.numeric(max_unexpected) || !is.finite(max_unexpected) ||
      max_unexpected < 0 || max_unexpected != floor(max_unexpected)) {
    stop("max_unexpected must be a nonnegative integer.", call. = FALSE)
  }
  if (!is.null(rater)) {
    if (!is.null(reference_plan)) {
      stop("Supply reference_plan to the assessment report. Use the rater card for the fitted rater's Fair(M).", call. = FALSE)
    }
    if (missing(title)) title <- "Individual rater feedback"
    return(.mfrm_render_rater_card(x, rater, criterion_facet, max_unexpected, format, file, title, digits))
  }
  conn <- x$connectivity
  sp <- x$spec
  rmin <- if (!is.null(x$fit)) x$fit$rating_min else sp$rating_min %||% min(as.numeric(x$data[[sp$score]]), na.rm = TRUE)
  rmax <- if (!is.null(x$fit)) x$fit$rating_max else sp$rating_max %||% max(as.numeric(x$data[[sp$score]]), na.rm = TRUE)
  fit_summary <- if (!is.null(x$fit)) summary(x$fit)
  readiness <- fit_summary$readiness
  model_label <- if (!is.null(fit_summary)) .tidymfrm_model_label(x$fit, fit_summary)
  calibration <- if (is.null(readiness)) {
    "No calibration result is attached. Run step_mfrm_measure() to estimate the measures."
  } else {
    reason <- if (length(readiness$reason_codes)) paste0(" (", paste(readiness$reason_codes, collapse = ", "), ")") else ""
    sprintf("Inference readiness is %s%s for the %s calibration. %s", readiness$status, reason, model_label,
      if (isTRUE(readiness$ready)) "Check model fit and rating design before interpreting intervals."
      else "Review the calibration before using interval or screening results.")
  }
  plan <- .tidymfrm_plan_overview(x$gtheory, digits)
  if (!is.null(x$gtheory)) {
    n_used <- x$gtheory$n_obs
    n_input <- x$gtheory$input_rows %||% nrow(x$data)
    plan$notes <- c(sprintf("The G-Study used %d of %d input rows. It omitted %d rows incomplete for its selected identifiers or scores. Source: %s.",
                            n_used, n_input, n_input - n_used, x$gtheory$input_source %||% "ratings"), plan$notes)
  } else if (!is.null(x$gtheory_reason)) {
    plan$notes <- paste("G-Study unavailable:", x$gtheory_reason)
  }
  if (!is.null(x$bootstrap_reason)) {
    calibration <- c(calibration, paste("Bootstrap review:", x$bootstrap_reason))
  }
  dg <- x$diagnostics
  review <- "No rater review is attached. Run step_mfrm_diagnose() after calibration."
  unexpected <- NULL
  if (!is.null(dg)) {
    flagged <- .tidymfrm_flag_count(dg$raters)
    unavailable <- sum(dg$raters$DiagnosticFlag == "Review required", na.rm = TRUE)
    review <- if (is.na(flagged)) sprintf("Rater screening is unavailable for %d raters.", unavailable) else
      sprintf("%d of %d evaluated raters triggered the screening rule (%d unavailable).",
              flagged, nrow(dg$raters) - unavailable, unavailable)
    review <- c(review, sprintf(paste("There are %d unexpected ratings with |standardized residual| >= 2.",
      "Review the score record and rubric before taking action."), nrow(dg$unexpected)))
    if (!is.null(dg$bias_reason)) review <- c(review, paste("Pairwise bias screening is unavailable.", dg$bias_reason))
    if (nrow(dg$unexpected) && max_unexpected > 0L) {
      fields <- intersect(c(".row_id", sp$person, sp$facets, "Observed", "Expected", "StdResidual", "ObsProb"), names(dg$unexpected))
      unexpected <- utils::head(dg$unexpected[fields], as.integer(max_unexpected))
      names(unexpected)[names(unexpected) == ".row_id"] <- "InputRow"
    }
  }
  design <- sprintf("%d complete ratings across %d persons. The rating graph has %d connected subset(s).",
                    conn$n_obs, conn$n_person, conn$n_subsets)
  if (nrow(x$data) > conn$n_obs) design <- paste(design, sprintf(
    "%d input rows were omitted from the connectivity check.", nrow(x$data) - conn$n_obs))
  anchors <- if (!is.null(x$fit)) {
    enc <- x$fit$encoding
    fixed <- sum(vapply(enc$anchor_fixed, function(a) sum(a == 1L), integer(1)))
    if (fixed || NROW(enc$group_anchors)) {
      sprintf("Calibration includes %d fixed elements and %d group anchor row(s).",
              fixed, NROW(enc$group_anchors))
    }
  }
  score_columns <- c("Level", "ObservedMin", "ObservedMax", "AtMinimum", "AtMaximum", "ScorePattern")
  slope_columns <- c("Facet", "Level", "Estimate", "ModelSE", "SandwichSE", "Formula")
  metadata <- .mfrm_report_metadata(x)
  report_tables <- .mfrm_report_tables(x, reference_plan, fair_level)
  notes <- c(metadata$notes,
    if (NROW(report_tables[["Separation, strata, and fixed-effect homogeneity"]]))
      "Homogeneity uses diagonal ModelSE variances. Read its nominal p-value as a screen, particularly for EAP or constrained measures.",
    if (NROW(report_tables[["Expected agreement and model-based kappa"]]))
      "Model-based agreement assumes independent ratings conditional on each fitted person's primary score and the fitted facets. Its exact kappa uses expected exact agreement on the same pairs.",
    if (NROW(report_tables[["Anchor displacement"]]))
      "Displacement is a one-step conditional change of a location or common group shift. Other measures stay fixed. ConditionalSE uses local information.",
    if (NROW(report_tables[["Linked benchmark comparison"]]))
      "Benchmark change is calculated on the stable reference anchored scale. Its standard error is conditional on the baseline calibration.",
    if (NROW(report_tables[["Common-plan Fair Average"]]))
      "Common-plan intervals transform normal ability bounds with calibration and ratings plan held fixed. EAP uses a conditional normal posterior approximation.",
    if (NROW(report_tables[["Effective weights at observed rating counts"]]))
      "Effective weights use the observed harmonic facet counts. A non-positive or invalid target covariance has no reported weights.",
    if (NROW(report_tables[["Rating network robustness"]]))
      "Rating network statistics describe observed assignments with one edge per person-rater pair.")
  details <- c(list("Calibration methods" = metadata$calibration), report_tables, list(
    "GMFRM Discrimination Parameters" = x$fit$slopes[, slope_columns, drop = FALSE],
    "Dynamic Markov Rater Drift Summary" = if (!is.null(x$drift)) as.data.frame(x$drift, component = "severity"),
    "Differential Facet Functioning (DFF) Group Screening" = dg$dff,
    "Observed Score Range and Endpoint Counts by Rater" = dg$raters[, score_columns, drop = FALSE],
    "Observed Rating Counts by Rater and Category" = dg$category_counts
  ))
  roster <- x$planned_coverage
  roster_note <- if (NROW(roster)) {
    sprintf("Planned roster: %d expected cells, %d not observed, %d observed outside the roster.",
            roster$ExpectedCells, roster$MissingExpectedCells, roster$UnexpectedObservedCells)
  }
  start_here <- c(design, roster_note, .tidymfrm_linkage_notes(x), calibration, anchors, plan$notes, review, notes)

  if (identical(format, "docx")) {
    if (!requireNamespace("officer", quietly = TRUE)) {
      stop("Install suggested package 'officer' to export `.docx` reports.", call. = FALSE)
    }
    if (is.null(file) || !nzchar(file)) {
      stop("Supply `file` (e.g., `tempfile(tmpdir = tempdir(), fileext = '.docx')`) when `format = 'docx'`.", call. = FALSE)
    }
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, title, style = "heading 1")
    doc <- officer::body_add_par(doc, "Start here", style = "heading 2")
    for (note in start_here) doc <- officer::body_add_par(doc, note, style = "Normal")
    add_table <- function(doc, heading, table) .mfrm_docx_table(doc, heading, table, digits)
    doc <- add_table(doc, "D-Study plans at the same rating count", plan$table)
    if (!is.null(x$fit)) {
      doc <- add_table(doc, "Facet Parameter Estimates",
                          x$fit$facets[, c("Facet", "Level", "ObsCount", "FairM", "Estimate", "ModelSE", "RealSE", "InfitMnSq")])
      doc <- add_table(doc, "Andrich Step Thresholds",
                          x$fit$steps[, c("ScaleLevel", "Step", "Estimate", "ModelSE")])
      doc <- add_table(doc, "Conditional Parametric Bootstrap Facet Uncertainty", x$bootstrap)
    }
    if (!is.null(x$gtheory)) {
      doc <- add_table(doc, "G-Study Variance Components", tidy(x$gtheory))
      doc <- add_table(doc, "D-Study Projections", x$gtheory$d_study)
    }
    if (!is.null(x$diagnostics)) {
      doc <- add_table(doc, "Rater Diagnostic Audit",
                          x$diagnostics$raters[, c("Level", "ObsCount", "FairM", "Estimate", "RealSE", "InfitMnSq", "DiagnosticFlag")])
    }
    for (heading in names(details)) doc <- add_table(doc, heading, details[[heading]])
    doc <- add_table(doc, "Unexpected ratings to review", unexpected)
    print(doc, target = file)
    return(invisible(file))
  }

  sec <- c(paste0("# ", title), "", "## Start here", "", paste0("- ", start_here), "")
  if (NROW(plan$table)) {
    sec <- c(sec, "### D-Study plans at the same rating count", "",
             .tidymfrm_fmt_df_md(plan$table, digits = digits), "")
  }

  # Section 1: Connectivity & Design
  sec <- c(
    sec,
    "## 1. Rating Design & Multi-Facet Graph Connectivity",
    sprintf("- **Total complete ratings**: %d across %d persons (`%s`).", conn$n_obs, conn$n_person, sp$person),
    sprintf("- **Non-person facets**: %s (designated rater facet: `%s`).",
            paste(paste0("`", sp$facets, "`"), collapse = ", "), sp$rater_facet),
    sprintf("- **Rubric category range**: %.0f to %.0f (%s).",
            rmin, rmax,
            if (isTRUE(conn$interior_categories_supported)) "all interior categories observed" else "WARNING: empty interior category detected"),
    sprintf("- **Bipartite / multi-facet connectivity**: **%s** (`n_subsets = %d`).",
            if (isTRUE(conn$connected)) "Multi-facet graph connected" else "DISCONNECTED SUBSETS",
            conn$n_subsets)
  )
  if (!isTRUE(conn$connected)) {
    sec <- c(
      sec, "",
      "> Caution: Disconnected rating subsets prevent direct cross-subset logit comparisons without common anchor raters, tasks, or benchmark performances (Linacre, 1989)."
    )
  }
  sec <- c(sec, "")

  # Section 2: Measurement Model
  if (!is.null(x$fit)) {
    f <- x$fit
    n_ext <- sum(f$person$ExtremeFlag != 0L)
    sec <- c(
      sec,
      "## 2. Many-Facet Measurement (Pillar 1)",
      sprintf("- **Model & Estimator**: `%s` estimated via `%s` (primary person measure: `%s`; converged = `%s` in %d iterations; `logLik = %.2f`).",
              model_label, f$method, f$score_method, f$converged, f$iterations, f$loglik),
      "- **Scale identification**: Non-person facets sum-centered to 0.00 logits unless anchored; facet signs follow MFRM convention (`Sign = -1` for severity/difficulty so higher values denote stricter raters or harder tasks).",
      sprintf("- **Extreme person scores**: %d/%d persons at minimum/maximum possible score.", n_ext, f$n_person)
    )
    if (is.finite(f$population$eap_reliability)) {
      sec <- c(
        sec,
        sprintf("- **Latent person distribution**: Population SD = %.3f logits; marginal latent EAP reliability = %.4f.",
                f$population$sigma, f$population$eap_reliability)
      )
    }
    sec <- c(sec, "", "### Non-Person Facet Parameter Summary", "")
    f_disp <- f$facets[, c("Facet", "Level", "ObsCount", "ObsAvg", "FairM", "Estimate", "ModelSE", "RealSE", "InfitMnSq", "OutfitMnSq"), drop = FALSE]
    sec <- c(sec, .tidymfrm_fmt_df_md(f_disp, digits = digits), "")

    sec <- c(sec, "### Andrich Step Thresholds", "")
    s_disp <- f$steps[, c("ScaleLevel", "Step", "FromCategory", "ToCategory", "Estimate", "ModelSE"), drop = FALSE]
    sec <- c(sec, .tidymfrm_fmt_df_md(s_disp, digits = digits), "")

    if (!is.null(f$explanatory) && nrow(f$explanatory) > 0L) {
      sec <- c(sec, "### Explanatory Facet Covariates", "", .tidymfrm_fmt_df_md(f$explanatory, digits = digits), "")
    }
    if (!is.null(x$method_comparison)) {
      sec <- c(sec, "### Multi-Estimator Sensitivity Check", "", .tidymfrm_fmt_df_md(x$method_comparison, digits = digits + 1L), "")
    }
    if (!is.null(x$bootstrap)) {
      sec <- c(sec, "### Conditional Parametric Bootstrap Facet Uncertainty", "",
               .tidymfrm_fmt_df_md(x$bootstrap, digits = digits), "")
    }
  }

  # Section 3: G-Study & D-Study
  if (!is.null(x$gtheory)) {
    gt <- x$gtheory
    gs <- gt$g_study_summary
    sec <- c(
      sec,
      "## 3. Generalizability (G-Study) & Decision (D-Study) Planning (Pillar 2)",
      sprintf("- **Observed-score generalizability**: Relative `G = %.4f` (`RelSEM = %.3f` rubric points); Absolute dependability `Phi = %.4f` (`AbsSEM = %.3f` rubric points).",
              gs$G, gs$RelSEM, gs$Phi, gs$AbsSEM),
      "> Methodological note: Observed-score `G` and `Phi` decompose raw rating variance via Henderson Method 1 (`urGENOVA`; Brennan, 2001) and include additive rater/task main-effect variance in `Phi`. They answer assessment-design allocation questions and are distinct from model-adjusted latent logit reliability (`EAP` reliability).",
      "",
      "### Variance Components",
      "",
      .tidymfrm_fmt_df_md(tidy(gt), digits = digits + 1L),
      "",
      "### D-Study Candidate Designs",
      "",
      .tidymfrm_fmt_df_md(gt$d_study, digits = digits),
      ""
    )
  }

  # Section 4: Rater Diagnostics, Fair(M), Real SE, Bias, DFF, & Rater Cards
  if (!is.null(x$diagnostics)) {
    dg <- x$diagnostics
    rt <- dg$raters
    sec <- c(
      sec,
      "## 4. Rater Diagnostics, Fair(M), Real SE, & Quality Assurance (Pillar 3)",
      sprintf(paste("A fit flag requires `Infit MnSq` outside `[%.2f, %.2f]`",
        "and either `|Infit ZSTD| >= %.1f` or unavailable ZSTD. Severity flags use a separate rule."),
              dg$mnsq_range[1L], dg$mnsq_range[2L], dg$zstd_limit),
      sprintf("- **Unavailable rater screens**: %d.", sum(rt$DiagnosticFlag == "Review required")),
      "",
      "### Rater Diagnostic Summary",
      "If you start with FACETS Table 7, read the measure (`Estimate`), `FairM`, SEs and Infit/Outfit here.",
      ""
    )
    rt_disp <- rt[, c("Level", "ObsCount", "ObsAvg", "FairM", "Estimate", "ModelSE", "RealSE",
                      "InfitMnSq", "InfitZSTD", "OutfitMnSq", "DiagnosticFlag", "InferenceReason"), drop = FALSE]
    sec <- c(sec, .tidymfrm_fmt_df_md(rt_disp, digits = digits), "")

    if (!is.null(dg$bias) && nrow(dg$bias$pairwise) > 0L) {
      pw_top <- utils::head(dg$bias$pairwise[order(-abs(dg$bias$pairwise$t)), , drop = FALSE], 10L)
      sec <- c(
        sec,
        sprintf("### Top Pairwise Bias Contrasts (%s)", dg$bias$test_label),
        "Like FACETS Table 14, these contrasts compare an element across two contexts. `pAdj` adjusts for multiple comparisons.",
        "",
        .tidymfrm_fmt_df_md(pw_top[, c("TargetFacet", "TargetLevel", "Context1", "Context2", "Contrast", "ContrastSE", "t", "df", "pAdj")], digits = digits),
        ""
      )
    }

    if (!is.null(dg$decisions) && nrow(dg$decisions) > 0L) {
      sec <- c(
        sec,
        "### Cut-Score Classification Sensitivity (`ModelSE` vs `RealSE`)",
        "These are conditional normal-approximation sensitivity counts, holding calibration fixed. Availability counts are the denominators; NA is not zero. CutLogit is NA when person-specific task mixtures imply different thresholds; CutLogitMin/Max show their range.",
        "",
        .tidymfrm_fmt_df_md(dg$decisions, digits = digits),
        ""
      )
    }

    if (NROW(unexpected)) {
      sec <- c(
        sec,
        sprintf("### Top Unexpected Ratings (`|StdResidual| >= 2.0`, up to %d rows)", max_unexpected),
        if ("InputRow" %in% names(unexpected)) "For a FACETS Table 4-style review, use `InputRow` to look up the record in the supplied data frame, counting from 1.",
        "",
        .tidymfrm_fmt_df_md(unexpected, digits = digits),
        ""
      )
    }

    sec <- c(
      sec,
      "> Validity & Operational Guidance: Statistical fit indices (`Infit`/`Outfit`), severity estimates, and DFF contrasts support formative rater training, rubric anchor clarification, and targeted second-scoring. They must not be used for automated rater removal without qualitative review of scoring rationales and construct representation (Eckes, 2015).",
      ""
    )

    if (isTRUE(include_rater_cards) && nrow(rt) > 0L) {
      sec <- c(sec, "## 5. Individual Rater Calibration Feedback Cards", "")
      for (i in seq_len(nrow(rt))) {
        r_id <- rt$Level[i]
        u_r <- dg$unexpected[as.character(dg$unexpected[[sp$rater_facet]]) == as.character(r_id), , drop = FALSE]
        sec <- c(
          sec,
          sprintf("### Rater `%s` -- %s", r_id, rt$DiagnosticFlag[i]),
          sprintf("- **Ratings scored**: %d (Observed average = %.2f; Model-adjusted `Fair(M)` = %.2f).",
                  as.integer(rt$ObsCount[i]), rt$ObsAvg[i], rt$FairM[i]),
          sprintf("- **Severity**: %.3f logits on the fitted scale (`ModelSE = %.3f`, `RealSE = %.3f`). Higher values indicate stricter ratings.",
                  rt$Estimate[i], rt$ModelSE[i], rt$RealSE[i]),
          sprintf("- **Consistency**: `Infit MnSq = %.3f` (`ZSTD = %.2f`), `Outfit MnSq = %.3f` (`ZSTD = %.2f`).",
                  rt$InfitMnSq[i], rt$InfitZSTD[i], rt$OutfitMnSq[i], rt$OutfitZSTD[i]),
          sprintf("- **Suggested focus**: %s (Unexpected ratings `|z| >= 2`: %d).",
                  rt$RecommendedAction[i], nrow(u_r)),
          ""
        )
      }
    }
  }

  if (any(lengths(details) > 0L)) {
    sec <- c(sec, "## Additional Diagnostic Tables", "")
    for (heading in names(details)) {
      table <- details[[heading]]
      if (NROW(table)) sec <- c(sec, paste0("### ", heading), "", .tidymfrm_fmt_df_md(table, digits = digits), "")
    }
  }
  md_out <- paste(sec, collapse = "\n")
  final_out <- if (format == "html") .tidymfrm_md_to_html(md_out, title = title) else md_out

  if (!is.null(file)) {
    con <- file(file, open = "wb", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(enc2utf8(final_out), con = con, useBytes = TRUE)
    attr(final_out, "file") <- file
    return(invisible(final_out))
  }
  final_out
}
