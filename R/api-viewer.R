mfrm_results_viewer_include_hint <- function(section) {
  hints <- c(
    diagnostics_summary = "Use include = \"standard\" or include = \"diagnostics\".",
    precision_review = "Use include = \"validation\" or include = \"precision\".",
    reporting_checklist = "Use include = \"publication\" or include = \"reporting\".",
    apa_outputs = "Use include = \"publication\" or include = \"apa\".",
    bias_screen = "Use include = \"bias\".",
    misfit_review = "Use include = \"misfit_review\".",
    response_time_review = "Use include = \"response_time\" with `response_time` metadata; for fitted objects also pass response_time_data."
  )
  hints[[section]] %||% "Use the corresponding include preset if needed."
}

mfrm_results_viewer_tab_status <- function(x,
                                           sections = character(0),
                                           plots = character(0),
                                           table_pattern = NULL) {
  status <- as.data.frame(x$status %||% data.frame(), stringsAsFactors = FALSE)
  plot_map <- as.data.frame(x$plot_map %||% data.frame(), stringsAsFactors = FALSE)
  tables <- x$tables %||% list()
  rows <- list()
  add <- function(source, name, status_value, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Source = as.character(source),
      Name = as.character(name),
      Status = as.character(status_value),
      Detail = as.character(detail %||% ""),
      stringsAsFactors = FALSE
    )
  }

  for (section in as.character(sections)) {
    row <- if (nrow(status) > 0L && all(c("Section", "Status") %in% names(status))) {
      status[status$Section %in% section, , drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(row) > 0L) {
      add("section", section, row$Status[1], row$Detail[1] %||% "")
    } else if (section %in% names(x$components %||% list())) {
      add("section", section, "ok", "Available in components.")
    } else {
      add(
        "section",
        section,
        "not_requested",
        paste0(
          "This section is not stored in this mfrm_results object. ",
          mfrm_results_viewer_include_hint(section)
        )
      )
    }
  }

  for (plot_type in as.character(plots)) {
    row <- if (nrow(plot_map) > 0L && all(c("Type", "Available") %in% names(plot_map))) {
      plot_map[plot_map$Type %in% plot_type, , drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(row) > 0L) {
      available <- isTRUE(row$Available[1])
      add(
        "plot",
        plot_type,
        if (available) "ok" else "not_available",
        row$Detail[1] %||% if (available) "Plot route is available." else "Plot route is not available."
      )
    } else {
      add("plot", plot_type, "not_available", "This plot route is not present in plot_map.")
    }
  }

  if (!is.null(table_pattern)) {
    table_names <- names(tables)
    matching <- table_names[grepl(table_pattern, table_names, ignore.case = TRUE)]
    add(
      "table",
      table_pattern,
      if (length(matching) > 0L) "ok" else "not_available",
      if (length(matching) > 0L) {
        paste0("Matching tables: ", paste(utils::head(matching, 8L), collapse = ", "), ".")
      } else {
        "No matching table was collected for this result object."
      }
    )
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  rank <- c(not_available = 1L, not_requested = 2L, review = 3L, ok = 4L)
  ord <- rank[out$Status]
  ord[is.na(ord)] <- 5L
  out[order(ord, out$Source, out$Name), , drop = FALSE]
}

mfrm_results_viewer_tab_statuses <- function(x) {
  list(
    qc = mfrm_results_viewer_tab_status(
      x,
      sections = c("diagnostics", "diagnostics_summary", "precision_review"),
      plots = "qc",
      table_pattern = "^(diagnostics_summary|precision_review|fit_measures_(status|summary|profile|underfit|overfit|mixed))"
    ),
    report = mfrm_results_viewer_tab_status(
      x,
      sections = c("reporting_checklist", "apa_outputs"),
      table_pattern = "^(reporting_checklist|apa_outputs)"
    ),
    bias = mfrm_results_viewer_tab_status(
      x,
      sections = "bias_screen",
      table_pattern = "(bias|dff)"
    ),
    misfit = mfrm_results_viewer_tab_status(
      x,
      sections = "misfit_review",
      plots = "pathway",
      table_pattern = "(unexpected|displacement|misfit_review|pathway)"
    ),
    response_time = mfrm_results_viewer_tab_status(
      x,
      sections = "response_time_review",
      plots = "response_time",
      table_pattern = "^response_time_review"
    )
  )
}

mfrm_results_viewer_response_time_boundary <- function() {
  data.frame(
    Area = c("Source", "Use", "Do not use as", "Next route"),
    Detail = c(
      "Response-time metadata supplied outside the fitted MFRM likelihood.",
      "Descriptive QC for rapid/slow timing patterns by event, person, facet, or score.",
      "A fitted speed parameter, joint speed-accuracy model, modified logit estimate, or automatic exclusion rule.",
      "Use response_time_review(), plot_response_time_review(..., draw = FALSE), and plot_data_components() for reusable timing QC tables."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_results_viewer_payload <- function(x, top_n = 100L) {
  if (!inherits(x, "mfrm_results")) {
    stop(
      "`x` must be an mfrm_results object. Call `mfrm_results()` first, ",
      "then pass that result to `launch_mfrmr_viewer()`.",
      call. = FALSE
    )
  }
  top_n <- max(1L, as.integer(top_n[1] %||% 100L))
  sx <- summary(x, top_n = top_n)
  tables <- x$tables %||% list()
  plot_map <- as.data.frame(x$plot_map %||% data.frame(), stringsAsFactors = FALSE)
  plot_choices <- character(0)
  if (nrow(plot_map) > 0L && all(c("Type", "Available") %in% names(plot_map))) {
    plot_choices <- unique(as.character(plot_map$Type[plot_map$Available %in% TRUE]))
  }
  qc_tables <- mfrm_results_viewer_table_subset(
    tables,
    "^(diagnostics_summary|precision_review|fit_measures_(status|summary|profile|underfit|overfit|mixed))"
  )
  report_tables <- mfrm_results_viewer_table_subset(
    tables,
    "^(reporting_checklist|apa_outputs)"
  )
  response_time_tables <- mfrm_results_viewer_table_subset(
    tables,
    "^response_time_review"
  )
  apa <- x$components$apa_outputs %||% NULL
  report_text <- if (inherits(apa, "mfrm_apa_outputs") && length(apa$report_text %||% character(0)) > 0L) {
    paste(as.character(apa$report_text), collapse = "\n\n")
  } else {
    paste(
      "APA-style report text is not available in this mfrm_results object.",
      "Create the object with include = \"publication\" or include = \"apa\",",
      "or call build_apa_outputs(fit, diagnostics) directly."
    )
  }
  report_notes <- if (inherits(apa, "mfrm_apa_outputs") && length(apa$table_figure_notes %||% character(0)) > 0L) {
    paste(as.character(apa$table_figure_notes), collapse = "\n\n")
  } else {
    "No APA table/figure notes are available in this mfrm_results object."
  }
  report_captions <- if (inherits(apa, "mfrm_apa_outputs") && length(apa$table_figure_captions %||% character(0)) > 0L) {
    paste(as.character(apa$table_figure_captions), collapse = "\n\n")
  } else {
    "No APA figure captions are available in this mfrm_results object."
  }
  bias_table <- mfrm_results_viewer_bias_table(x)
  unexpected_table <- mfrm_results_viewer_unexpected_table(x, top_n = top_n)
  list(
    summary = sx,
    tables = tables,
    table_names = names(tables),
    qc_tables = qc_tables,
    qc_table_names = names(qc_tables),
    report_tables = report_tables,
    report_table_names = names(report_tables),
    response_time_tables = response_time_tables,
    response_time_table_names = names(response_time_tables),
    response_time_boundary = mfrm_results_viewer_response_time_boundary(),
    report_text = report_text,
    report_notes = report_notes,
    report_captions = report_captions,
    bias_table = bias_table,
    bias_guidance = mfrm_results_viewer_bias_guidance(x),
    unexpected_table = unexpected_table,
    unexpected_choices = mfrm_results_viewer_unexpected_choices(unexpected_table),
    tab_status = mfrm_results_viewer_tab_statuses(x),
    plot_choices = plot_choices,
    replay_code = paste(as.character(sx$reproducible_code$Code %||% ""), collapse = "\n")
  )
}

mfrm_results_viewer_table_subset <- function(tables, pattern) {
  tables <- tables %||% list()
  if (length(tables) == 0L) {
    return(list())
  }
  keep <- grepl(pattern, names(tables), ignore.case = TRUE)
  tables[keep]
}

mfrm_results_viewer_table_or_empty <- function(x) {
  x <- as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(x) == 0L && ncol(x) == 0L) {
    return(data.frame(Message = "No rows available.", stringsAsFactors = FALSE))
  }
  x
}

mfrm_results_viewer_bias_table <- function(x) {
  direct <- x$components$bias_screen$table %||%
    x$components$bias_interaction %||%
    x$components$bias %||%
    x$diagnostics$bias %||%
    x$tables$bias_table %||%
    NULL
  if (is.null(direct)) {
    return(data.frame())
  }
  as.data.frame(direct, stringsAsFactors = FALSE)
}

mfrm_results_viewer_bias_guidance <- function(x) {
  has_bias <- nrow(mfrm_results_viewer_bias_table(x)) > 0L
  data.frame(
    Area = c("Facet-level bias screen", "Interaction bias screen", "Reporting boundary"),
    Status = c(
      if (has_bias) "available" else "not_available",
      "requires_explicit_facet_pair",
      "screening_not_final_fairness_decision"
    ),
    Route = c(
      "diagnose_mfrm(fit)$bias or estimate_bias(fit)",
      "estimate_bias(fit, facet1 = ..., facet2 = ...) -> bias_interaction_report()",
      "Use reporting_checklist() and substantive review before fairness claims."
    ),
    Detail = c(
      if (has_bias) "A bias table is available for this result object." else "No bias table is stored in this result object.",
      "The viewer does not choose a facet pair automatically.",
      "Bias outputs are conditional screening layers, not standalone validity decisions."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_results_viewer_unexpected_table <- function(x, top_n = 100L) {
  tbl <- x$components$unexpected$table %||%
    x$components$misfit_review$unexpected$table %||%
    x$components$unexpected %||%
    x$tables$unexpected_table %||%
    x$tables$misfit_review_unexpected_table %||%
    NULL
  if (is.null(tbl)) {
    return(data.frame())
  }
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (nrow(tbl) > top_n) {
    tbl <- utils::head(tbl, top_n)
  }
  tbl
}

mfrm_results_viewer_unexpected_choices <- function(tbl) {
  tbl <- as.data.frame(tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    return(character(0))
  }
  label_cols <- intersect(c("Row", "Person", "Rater", "Criterion", "Severity", "StdResidual"), names(tbl))
  labels <- vapply(seq_len(nrow(tbl)), function(i) {
    parts <- vapply(label_cols, function(col) paste0(col, "=", as.character(tbl[[col]][i])), character(1))
    paste(parts, collapse = " | ")
  }, character(1))
  stats::setNames(as.character(seq_len(nrow(tbl))), labels)
}

mfrm_results_shiny_app <- function(x, top_n = 100L) {
  payload <- mfrm_results_viewer_payload(x, top_n = top_n)
  has_tables <- length(payload$table_names) > 0L
  has_plots <- length(payload$plot_choices) > 0L
  has_qc_tables <- length(payload$qc_table_names) > 0L
  has_report_tables <- length(payload$report_table_names) > 0L
  has_unexpected <- nrow(payload$unexpected_table) > 0L
  has_response_time_tables <- length(payload$response_time_table_names) > 0L
  has_response_time_plot <- "response_time" %in% payload$plot_choices

  table_panel <- if (has_tables) {
    shiny::tagList(
      shiny::selectInput("table_name", "Table", choices = payload$table_names),
      shiny::downloadButton("download_table", "Download CSV"),
      shiny::tableOutput("selected_table")
    )
  } else {
    shiny::tags$p("No table outputs are available for this mfrm_results object.")
  }

  plot_panel <- if (has_plots) {
    shiny::tagList(
      shiny::selectInput("plot_type", "Plot", choices = payload$plot_choices),
      shiny::plotOutput("result_plot", height = "650px")
    )
  } else {
    shiny::tags$p("No plot routes are available for this mfrm_results object.")
  }

  qc_table_panel <- if (has_qc_tables) {
    shiny::tagList(
      shiny::selectInput("qc_table_name", "QC table", choices = payload$qc_table_names),
      shiny::tableOutput("qc_table")
    )
  } else {
    shiny::tags$p("No dedicated QC tables are available for this mfrm_results object.")
  }

  report_table_panel <- if (has_report_tables) {
    shiny::tagList(
      shiny::selectInput("report_table_name", "Report table", choices = payload$report_table_names),
      shiny::tableOutput("report_table")
    )
  } else {
    shiny::tags$p("No reporting or APA tables are available for this mfrm_results object.")
  }

  response_time_panel <- if (has_response_time_tables) {
    shiny::tagList(
      shiny::selectInput("response_time_table_name", "Response-time table",
                         choices = payload$response_time_table_names),
      shiny::tableOutput("response_time_table")
    )
  } else {
    shiny::tags$p("No response-time tables are available for this mfrm_results object.")
  }

  misfit_panel <- if (has_unexpected) {
    shiny::tagList(
      shiny::selectInput("unexpected_row", "Unexpected response", choices = payload$unexpected_choices),
      shiny::h4("Selected observation"),
      shiny::tableOutput("unexpected_selected"),
      shiny::h4("Unexpected-response table"),
      shiny::tableOutput("unexpected_table")
    )
  } else {
    shiny::tags$p("No unexpected-response rows are available for this mfrm_results object.")
  }

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML(
        paste(
          ".mfrmr-viewer .tab-content { padding-top: 16px; }",
          ".mfrmr-viewer table { font-size: 13px; }",
          ".mfrmr-viewer pre { white-space: pre-wrap; }"
        )
      ))
    ),
    shiny::div(
      class = "mfrmr-viewer",
      shiny::titlePanel("mfrmr results"),
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Overview",
          shiny::tableOutput("overview"),
          shiny::h4("Triage"),
          shiny::tableOutput("triage")
        ),
        shiny::tabPanel(
          "Status",
          shiny::tableOutput("status"),
          shiny::h4("Next actions"),
          shiny::tableOutput("next_actions")
        ),
        shiny::tabPanel(
          "QC",
          shiny::h4("Section status"),
          shiny::tableOutput("qc_status"),
          shiny::plotOutput("qc_plot", height = "650px"),
          shiny::h4("QC evidence"),
          qc_table_panel
        ),
        shiny::tabPanel(
          "Report",
          shiny::h4("Section status"),
          shiny::tableOutput("report_status"),
          shiny::h4("APA-style report text"),
          shiny::verbatimTextOutput("report_text"),
          shiny::h4("Table and figure notes"),
          shiny::verbatimTextOutput("report_notes"),
          shiny::h4("Figure captions"),
          shiny::verbatimTextOutput("report_captions"),
          shiny::h4("Reporting tables"),
          report_table_panel
        ),
        shiny::tabPanel(
          "Bias",
          shiny::h4("Section status"),
          shiny::tableOutput("bias_status"),
          shiny::h4("Bias guidance"),
          shiny::tableOutput("bias_guidance"),
          shiny::h4("Available bias table"),
          shiny::tableOutput("bias_table")
        ),
        shiny::tabPanel(
          "Pathway/Misfit",
          shiny::h4("Section status"),
          shiny::tableOutput("misfit_status"),
          shiny::plotOutput("pathway_plot", height = "650px"),
          misfit_panel
        ),
        shiny::tabPanel(
          "Response Time",
          shiny::h4("Section status"),
          shiny::tableOutput("response_time_status"),
          shiny::h4("Interpretation boundary"),
          shiny::tableOutput("response_time_boundary"),
          shiny::plotOutput("response_time_plot", height = "500px"),
          shiny::h4("Response-time tables"),
          response_time_panel
        ),
        shiny::tabPanel("Tables", table_panel),
        shiny::tabPanel("Plots", plot_panel),
        shiny::tabPanel(
          "Replay",
          shiny::verbatimTextOutput("replay_code")
        )
      )
    )
  )

  server <- function(input, output, session) {
    output$overview <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$summary$overview),
      rownames = FALSE
    )
    output$triage <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$summary$triage),
      rownames = FALSE
    )
    output$status <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$summary$status),
      rownames = FALSE
    )
    output$next_actions <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$summary$next_actions),
      rownames = FALSE
    )
    output$replay_code <- shiny::renderText(payload$replay_code)
    output$report_text <- shiny::renderText(payload$report_text)
    output$report_notes <- shiny::renderText(payload$report_notes)
    output$report_captions <- shiny::renderText(payload$report_captions)
    output$qc_status <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$tab_status$qc),
      rownames = FALSE
    )
    output$report_status <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$tab_status$report),
      rownames = FALSE
    )
    output$bias_status <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$tab_status$bias),
      rownames = FALSE
    )
    output$misfit_status <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$tab_status$misfit),
      rownames = FALSE
    )
    output$response_time_status <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$tab_status$response_time),
      rownames = FALSE
    )
    output$response_time_boundary <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$response_time_boundary),
      rownames = FALSE
    )
    output$bias_guidance <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$bias_guidance),
      rownames = FALSE
    )
    output$bias_table <- shiny::renderTable(
      mfrm_results_viewer_table_or_empty(payload$bias_table),
      rownames = FALSE
    )

    output$qc_plot <- shiny::renderPlot({
      if (!inherits(x$diagnostics, "mfrm_diagnostics")) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "QC plot requires diagnostics in the mfrm_results object.")
        return(invisible(NULL))
      }
      tryCatch(
        plot_qc_dashboard(x$fit, diagnostics = x$diagnostics, preset = "publication"),
        error = function(e) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, conditionMessage(e))
        }
      )
    })

    output$pathway_plot <- shiny::renderPlot({
      tryCatch(
        plot(x$fit, type = "pathway", preset = "publication"),
        error = function(e) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, conditionMessage(e))
        }
      )
    })

    output$response_time_plot <- shiny::renderPlot({
      if (!isTRUE(has_response_time_plot)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Response-time plot requires include = \"response_time\" with timing metadata.")
        return(invisible(NULL))
      }
      tryCatch(
        plot(x, type = "response_time", preset = "publication"),
        error = function(e) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, conditionMessage(e))
        }
      )
    })

    if (has_qc_tables) {
      output$qc_table <- shiny::renderTable({
        nm <- input$qc_table_name
        if (is.null(nm) || !nm %in% names(payload$qc_tables)) {
          return(data.frame(Message = "Select a QC table.", stringsAsFactors = FALSE))
        }
        mfrm_results_viewer_table_or_empty(payload$qc_tables[[nm]])
      }, rownames = FALSE)
    }

    if (has_report_tables) {
      output$report_table <- shiny::renderTable({
        nm <- input$report_table_name
        if (is.null(nm) || !nm %in% names(payload$report_tables)) {
          return(data.frame(Message = "Select a report table.", stringsAsFactors = FALSE))
        }
        mfrm_results_viewer_table_or_empty(payload$report_tables[[nm]])
      }, rownames = FALSE)
    }

    if (has_response_time_tables) {
      output$response_time_table <- shiny::renderTable({
        nm <- input$response_time_table_name
        if (is.null(nm) || !nm %in% names(payload$response_time_tables)) {
          return(data.frame(Message = "Select a response-time table.", stringsAsFactors = FALSE))
        }
        mfrm_results_viewer_table_or_empty(payload$response_time_tables[[nm]])
      }, rownames = FALSE)
    }

    if (has_unexpected) {
      output$unexpected_selected <- shiny::renderTable({
        idx <- suppressWarnings(as.integer(input$unexpected_row %||% NA_integer_))
        if (!is.finite(idx) || idx < 1L || idx > nrow(payload$unexpected_table)) {
          return(data.frame(Message = "Select an unexpected response.", stringsAsFactors = FALSE))
        }
        payload$unexpected_table[idx, , drop = FALSE]
      }, rownames = FALSE)
      output$unexpected_table <- shiny::renderTable(
        mfrm_results_viewer_table_or_empty(payload$unexpected_table),
        rownames = FALSE
      )
    }

    if (has_tables) {
      output$selected_table <- shiny::renderTable({
        nm <- input$table_name
        if (is.null(nm) || !nm %in% names(payload$tables)) {
          return(data.frame(Message = "Select a table.", stringsAsFactors = FALSE))
        }
        mfrm_results_viewer_table_or_empty(payload$tables[[nm]])
      }, rownames = FALSE)

      output$download_table <- shiny::downloadHandler(
        filename = function() {
          nm <- input$table_name %||% "mfrmr_table"
          paste0(gsub("[^A-Za-z0-9_.-]+", "_", nm), ".csv")
        },
        content = function(file) {
          nm <- input$table_name
          tbl <- if (!is.null(nm) && nm %in% names(payload$tables)) {
            payload$tables[[nm]]
          } else {
            data.frame()
          }
          utils::write.csv(as.data.frame(tbl, stringsAsFactors = FALSE), file, row.names = FALSE)
        }
      )
    }

    if (has_plots) {
      output$result_plot <- shiny::renderPlot({
        type <- input$plot_type
        if (is.null(type) || !type %in% payload$plot_choices) {
          return(invisible(NULL))
        }
        tryCatch(
          plot(x, type = type),
          error = function(e) {
            graphics::plot.new()
            graphics::text(0.5, 0.5, conditionMessage(e))
          }
        )
      })
    }
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Launch a local Shiny viewer for an mfrm_results object
#'
#' @description
#' `launch_mfrmr_viewer()` opens a local Shiny app for reading the object
#' returned by [mfrm_results()]. It is intentionally a viewer over an existing
#' comprehensive results object, not a new estimation interface.
#'
#' @param x An [mfrm_results()] object.
#' @param top_n Maximum number of table-index rows shown in the summary payload.
#' @param launch.browser Passed to `shiny::runApp()`.
#' @param port Optional port passed to `shiny::runApp()`. `NULL` lets Shiny
#'   choose its default.
#' @param host Host passed to `shiny::runApp()`. Defaults to the local host.
#' @param display.mode Passed to `shiny::runApp()`.
#' @param return_app Logical; if `TRUE`, return the Shiny app object without
#'   running it. This is useful for embedding or testing.
#' @param ... Additional arguments passed to `shiny::runApp()` when
#'   `return_app = FALSE`.
#'
#' @details
#' The viewer assumes that fitting, diagnostics, and section selection have
#' already happened through [mfrm_results()]. This keeps GUI exploration
#' separate from reproducible analysis setup: the Replay tab displays the
#' `mfrm_results()` scaffold stored in the result object.
#'
#' The app includes tabs for overview/triage, QC evidence, APA-style report
#' text when `include = "publication"` or `"apa"` was used, available bias
#' screens, pathway plotting, unexpected-response inspection, generic tables,
#' generic plot routes, and replay code. QC, Report, Bias, and Pathway/Misfit
#' tabs show local section-status tables so unavailable or not-requested
#' sections are visible where users look for them. Bias-interaction follow-up
#' still requires an explicit facet-pair decision outside the viewer.
#'
#' `shiny` is an optional dependency. Install it before using this viewer:
#' `install.packages("shiny")`.
#'
#' @return Invisibly returns the value from `shiny::runApp()`, or the Shiny app
#'   object when `return_app = TRUE`.
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 30
#' )
#' res <- mfrm_results(fit, include = c("fit", "diagnostics", "tables"))
#'
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   launch_mfrmr_viewer(res)
#' }
#' @seealso [mfrm_results()], [mfrm_results_interactive()]
#' @export
launch_mfrmr_viewer <- function(x,
                                top_n = 100L,
                                launch.browser = TRUE,
                                port = NULL,
                                host = getOption("shiny.host", "127.0.0.1"),
                                display.mode = c("auto", "normal", "showcase"),
                                return_app = FALSE,
                                ...) {
  if (!inherits(x, "mfrm_results")) {
    stop(
      "`x` must be an mfrm_results object. Call `mfrm_results()` first, ",
      "then pass that result to `launch_mfrmr_viewer()`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "`launch_mfrmr_viewer()` requires the optional package `shiny`. ",
      "Install it with `install.packages(\"shiny\")`.",
      call. = FALSE
    )
  }
  display.mode <- match.arg(display.mode)
  app <- mfrm_results_shiny_app(x, top_n = top_n)
  if (isTRUE(return_app)) {
    return(app)
  }

  args <- c(
    list(
      appDir = app,
      launch.browser = launch.browser,
      host = host,
      display.mode = display.mode
    ),
    list(...)
  )
  if (!is.null(port)) {
    args$port <- port
  }
  invisible(do.call(shiny::runApp, args))
}
