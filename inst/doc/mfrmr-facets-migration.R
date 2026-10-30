## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----facets-import-data-------------------------------------------------------
# library(mfrmr)
# score_directory <- system.file("extdata", "facets-example", package = "mfrmr")
# facets_result <- read_facets_results(score_directory,
#   facet_map = c("Person", "Rater", "Task"), source_version = "4.5.0")
# summary(facets_result)
# facets_result$table[, c("Facet", "ElementID", "Label", "Estimate", "SourceStatus")]

## ----facets-import-measures, fig.alt="Reported FACETS measures for six synthetic Persons. Downward and upward triangles mark extreme-score adjustments, a square marks a direct/group anchor flag whose subtype is not determined by the import, and circles mark estimated measures. No confidence intervals are shown."----
# plot(facets_result, type = "measure", facets = "Person")

## ----facets-import-availability-----------------------------------------------
# fit_points <- plot_data(facets_result, component = "plot_data",
#   type = "infit_outfit", facets = "Person")
# fit_points[, c("Label", "Infit", "Outfit", "Available", "Reason")]
# plot_data(facets_result, component = "notes")

## ----facets-import-display-labels---------------------------------------------
# display_names <- data.frame(Facet = "Person", ElementID = c("1", "2"),
#   DisplayLabel = c("受検者 甲", "受検者 乙"))
# labelled <- plot(facets_result, facets = "Person",
#   label_map = display_names, draw = FALSE)
# labelled$data$plot_data[, c("ElementID", "Label", "DisplayLabel")]

## ----facets-import-save-japanese, eval=FALSE----------------------------------
# png("facets-japanese.png", width = 7, height = 5, units = "in", res = 200,
#   type = "cairo", family = "Meiryo")
# plot(facets_result, facets = "Person", label_map = display_names)
# dev.off()
# 
# cairo_pdf("facets-japanese.pdf", width = 7, height = 5, family = "Meiryo")
# plot(facets_result, facets = "Person", label_map = display_names)
# dev.off()

## ----facets-import-graph------------------------------------------------------
# graph_path <- system.file("extdata", "facets-graph-example", "graph.txt",
#   package = "mfrmr")
# facets_curves <- read_facets_results(graph_file = graph_path,
#   source_version = "4.5.0")
# summary(facets_curves, component = "curves")

## ----facets-import-category-curves, fig.alt="Four reported category-probability curves for Scale 1 of a synthetic FACETS PCM analysis, on the original Graphfile measure axis. Categories 0 through 3 are distinguished by color and line type. Lines join exported points; there are no confidence intervals."----
# plot(facets_curves, type = "ccc", scales = 1)

## ----facets-import-graph-values-----------------------------------------------
# expected_points <- plot_data(facets_curves, "plot_data",
#   type = "expected", scales = 1)
# head(expected_points[, c("Scale", "Measure", "Value", "Available", "Reason")])
# head(facets_curves$graph$diagnostics)

## ----facets-mode--------------------------------------------------------------
# library(mfrmr)
# data("mfrmr_example_operational", package = "mfrmr")
# 
# run <- run_mfrm_facets(
#   data = mfrmr_example_operational,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   model = "RSM",
#   method = "JML"
# )
# 
# names(run)

## ----facets-mode-summary------------------------------------------------------
# jml_status <- summary(run$fit, profile = "fit", detail = "brief")
# jml_status$overview[, c(
#   "Model", "Method", "Converged", "InferenceReady",
#   "ConvergenceSeverity"
# )]
# jml_status$readiness
# head(run$fair_average)

## ----anchor-origin------------------------------------------------------------
# ratings <- load_mfrmr_data("example_core")
# rater_anchor <- data.frame(Facet = "Rater", Level = "R01", Anchor = 0.7)
# anchored_fit <- fit_mfrm(
#   ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#   method = "JML", anchors = rater_anchor, centering_scope = "anchor_origin"
# )
# anchored_fit$config$anchor_review$facet_summary

## ----group-anchor-origin------------------------------------------------------
# ratings <- load_mfrmr_data("example_core")
# group <- data.frame(
#   Facet = "Rater", Level = c("R01", "R02"), Group = "panel",
#   GroupValue = 0.6
# )
# group_fit <- fit_mfrm(
#   ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#   method = "JML", group_anchors = group, centering_scope = "anchor_origin",
#   anchor_policy = "error"
# )
# group_fit$config$anchor_review$facet_summary
# subset(group_fit$facets$others, Facet == "Rater")

