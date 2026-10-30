## ----include = FALSE--------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----plot-gallery, eval = TRUE, echo = FALSE, results = 'asis'--------------
figures <- mfrmr::mfrmr_output_guide("plots")
# Only presentation locations live here; questions and support come from the guide.
gallery <- data.frame(
  PlotName = c("wright_map", "fit_pathway", "category_characteristic_curves",
    "subset_connectivity", "cluster_dendrogram", "multivariate_d_study"),
  Chunk = c("wright", "fit-pathway", "category-curves", "linking",
    "feature-hierarchy", "d-study-planning"),
  Preview = c("Wright map", "Infit versus location", "Category probability curves",
    "Subset coverage", "Feature dendrogram", "D-study coefficient curves"))
rows <- figures[match(gallery$PlotName, figures$PlotName), ]
stopifnot(!anyNA(rows$Question))
support <- c(dedicated = "Base R and dedicated ggplot conversion",
  native = "Returns ggplot directly",
  generic = "Base R; ggplot offers a different, table-based view",
  unavailable = "Base R; custom graphics from the saved data")
cat('::: {.mfrmr-gallery}

')
for (i in seq_len(nrow(gallery))) {
  preview <- if (!is_cran_check) paste0("!", "[", gallery$Preview[i],
    " preview]", "(", knitr::fig_chunk(gallery$Chunk[i], "png", number = 1),
    "){loading=\"lazy\"} ") else ""
  cat(sprintf('::: {.mfrmr-gallery-card}

[%s**%s**](#figure-%s)

%s

:::

',
    preview, rows$Question[i], gallery$Chunk[i], support[rows$GGPlot[i]]))
}
cat(':::
')

## ----plot-guide, eval = TRUE------------------------------------------------
figures <- mfrmr::mfrmr_output_guide("plots")
knitr::kable(figures[, c("Question", "ResultFunction", "GGPlot")])

## ----population-interval-views, eval=FALSE----------------------------------
# summary(ci)                           # estimates, intervals and status
# p <- plot(ci, type = "contrasts")     # the saved comparison, without refitting
# plot_data(p, "contrast_weights")      # which conditional means were compared?
# plot_data(p, "population_coding")     # factor levels and coding, if applicable
# plot_data(p, "notes")                 # assumptions and experimental limits
# as_ggplot(p)                          # the same saved values for customization

## ----setup------------------------------------------------------------------
# library(mfrmr)
# 
# toy <- load_mfrmr_data("example_operational")
# 
# fit <- fit_mfrm(
#   toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM"
# )
# 
# diag <- diagnose_mfrm(fit, residual_pca = "none")
# checklist <- reporting_checklist(fit, diagnostics = diag)
# subset(
#   checklist$checklist,
#   Section == "Visual Displays",
#   c("Item", "Available", "NextAction")
# )

## ----wright, dev="png", fig.alt = "Wright map of person abilities, facet levels and category steps on the fitted logit scale, with requested uncertainty intervals. Higher facet locations indicate stricter scoring or greater difficulty."----
# wright <- plot(fit, type = "wright", preset = "publication", show_ci = TRUE)

## ----wright-values----------------------------------------------------------
# head(plot_data(wright, component = "locations"))

## ----wright-interpretation--------------------------------------------------
# knitr::kable(plot_data(wright, "notes"), row.names = FALSE)

## ----wright-highlight-------------------------------------------------------
# judge <- as.character(toy$Rater[1])
# focused <- plot(fit, type = "wright", highlight = list(Rater = judge),
#                 draw = FALSE)
# subset(plot_data(focused, "locations"), Highlighted)

## ----wright-facets-style, fig.alt = "Text-style Wright display with rubric labels Beginning, Developing, Secure and Advanced. Locations share a fitted logit scale; labels describe ordered categories rather than equal observed-score intervals."----
# plot(
#   fit,
#   type = "wright",
#   renderer = "facets",
#   category_labels = c(
#     `1` = "Beginning", `2` = "Developing", `3` = "Secure", `4` = "Advanced"
#   )
# )

## ----pathway, fig.alt = "Pathway display of fitted response probabilities across ability, separated by display layer. It shows model predictions rather than tests of rater bias."----
# plot(fit, type = "pathway", preset = "publication")

## ----fit-pathway, dev="png", fig.alt = "Infit mean squares against fitted locations for all facet levels and the 12 persons furthest from expected fit, with requested location intervals. Labels are limited to flagged facet levels; the accompanying table identifies the points. Mean-square bands are descriptive review thresholds, not automatic decisions."----
# fit_pathway <- plot(
#   fit,
#   type = "fit_pathway",
#   diagnostics = diag,
#   fit_stat = "Infit",
#   fit_scale = "mnsq",
#   include_person = TRUE,
#   top_n_person = 12,
#   person_labels = "none",
#   facet_labels = "flagged",
#   show_ci = TRUE,
#   preset = "publication"
# )

## ----fit-pathway-values-----------------------------------------------------
# fit_values <- plot_data(fit_pathway, component = "table")
# fit_values[, c("Facet", "Level", "FitValue", "Measure", "CIEligible")]

## ----category-curves, dev="png", fig.alt="Predicted category probabilities across ability for the fitted rating scale. Line patterns and colours identify categories; the curves use a reference profile with other additive facet effects fixed at zero."----
# category_curves <- plot(fit, type = "ccc", preset = "monochrome")

## ----category-values--------------------------------------------------------
# head(plot_data(category_curves, component = "probabilities"))
# # Optional customization without refitting:
# # as_ggplot(category_curves)

## ----unexpected, fig.alt = "Observed ratings screened by standardized residual size and fitted response probability. Highlighted observations merit contextual review rather than automatic deletion."----
# plot(fit, type = "unexpected", diagnostics = diag, abs_z_min = 1.5, prob_max = 0.4, view = "scatter", preset = "publication")

## ----displacement, fig.alt = "Level-specific displacement summaries for the fitted facets, including unanchored levels. The displayed shifts are diagnostic summaries, not evidence that a level changed over time."----
# plot(fit, type = "displacement", diagnostics = diag, anchored_only = FALSE, view = "lollipop", preset = "publication")

## ----strict-marginal, fig.alt = "Marginal fit summaries using posterior-averaged response probabilities. These diagnostics have a different probability basis from conditional residual screens."----
# fit_strict <- fit_mfrm(
#   toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7,
#   maxit = 40
# )
# 
# diag_strict <- diagnose_mfrm(
#   fit_strict,
#   residual_pca = "none",
#   diagnostic_mode = "both"
# )
# 
# strict_checklist <- reporting_checklist(fit_strict, diagnostics = diag_strict)
# subset(
#   strict_checklist$checklist,
#   Section == "Visual Displays" &
#     Item %in% c("QC / facet dashboard", "Strict marginal visuals"),
#   c("Item", "Available", "NextAction")
# )
# 
# plot(diag_strict, type = "marginal_fit", top_n = 12, preset = "publication")

## ----linking, dev="png", fig.alt = "Observation shares above a matrix of relative facet-level coverage across observed subsets. Matrix cells compare each subset with the largest observed count for that facet; they are not completion rates for planned assignments."----
# sc <- subset_connectivity_report(fit, diagnostics = diag)
# coverage <- plot(sc, type = "design_matrix", preset = "publication")

## ----coverage-values--------------------------------------------------------
# plot_data(coverage, component = "matrix")
# head(plot_data(coverage, component = "table"))
# # Only the matrix is converted; the observation-share panel is not included.
# # as_ggplot(coverage, component = "matrix")

## ----eval = FALSE-----------------------------------------------------------
# drift <- detect_anchor_drift(current_fit, baseline = baseline_anchors)
# plot(drift, type = "heatmap", preset = "publication")

## ----residual-pca, fig.alt = "Scree display of residual principal components when the requested decomposition is available. Residual components are diagnostic summaries, not automatically substantive dimensions."----
# diag_pca <- diagnose_mfrm(fit, residual_pca = "none")
# pca <- analyze_residual_pca(diag_pca, mode = "both", pca_max_factors = 4)
# summary(pca)
# if (nrow(pca$overall_table) > 0L) {
#   plot(pca, type = "residual_pca", mode = "overall", view = "scree", preset = "publication")
# }

## ----dimensionality-workflow, eval=FALSE------------------------------------
# diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
# pca <- analyze_residual_pca(diag, mode = "both", parallel = TRUE,
#                            parallel_reps = 200, seed = 9028)
# summary(pca)
# # Replace Criterion with the fitted facet relevant to the question.
# q3 <- q3_statistic(fit, diagnostics = diag, facet = "Criterion")
# q3$pairs
# q3$notes

## ----dimensionality-model-reference, eval=FALSE-----------------------------
# # fit must be an RSM/PCM MML fit with fixed standard-normal population,
# # fixed quadrature, additive facets, unit weights, and no anchors or shrinkage.
# reference <- analyze_residual_pca(
#   fit, mode = "both", parallel = TRUE,
#   parallel_method = "model_bootstrap", parallel_reps = 200, seed = 9028
# )
# reference$bootstrap_settings
# reference$bootstrap_trials       # All attempts, including failures and warnings
# reference$parallel_status        # Availability separately for each PCA scope
# plot(reference, type = "overall_parallel_scree")
# # If the reference is available:
# # as_ggplot(plot(reference, type = "overall_parallel_scree", draw = FALSE))

## ----dimensionality-criterion-reference, eval=FALSE-------------------------
# reference$parallel_status
# # For a prespecified criterion-level question, if that reference is available:
# plot(reference, type = "facet_parallel_scree", facet = "Criterion")

## ----residual-network-route, eval=FALSE-------------------------------------
# # Reuse a fitted model and its matching diagnostics.
# net <- mfrm_residual_network(fit, diagnostics = diag, facet = "Criterion",
#   missing = "complete", clusters = 2, min_abs = .2)
# net$settings
# net$pairs                    # include unavailable pairs and shared-person N
# net$person_usage             # complete-person inclusion, not rating-row N
# plot(net, type = "network")  # sign/strength; circular layout is not a distance
# plot(net, type = "matrix")   # all correlations, before the display threshold
# if (net$clustering$status == "available") plot(net, type = "dendrogram")
# report <- mfrm_report(net)
# export_summary_appendix(report, output_dir = "residual-review")
# saveRDS(net, "residual-review/network.rds")

## ----bias, fig.alt = "Rater-by-criterion interaction profiles from the fitted bias analysis. Differences identify combinations for review and do not establish a causal explanation."----
# bias_df <- load_mfrmr_data("example_bias")
# 
# fit_bias <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7
# )
# 
# diag_bias <- diagnose_mfrm(fit_bias, residual_pca = "none")
# bias <- estimate_bias(fit_bias, diag_bias, facet_a = "Rater", facet_b = "Criterion")
# 
# bias_figure <- plot(bias, type = "bias", view = "facet_profile", preset = "publication")

## ----custom-plot-data, fig.alt = "Monochrome Wright map of person, facet and step locations on the common logit scale. The accompanying code extracts the plotted locations for custom displays."----
# plot(fit, type = "wright", preset = "monochrome")
# 
# wright_payload <- plot(fit, type = "wright", draw = FALSE, preset = "publication")
# plot_data_components(wright_payload)
# 
# locations <- plot_data(wright_payload, component = "locations")
# head(locations)
# 
# pathway_long <- plot_data(
#   fit,
#   type = "pathway",
#   component = "pathway_long",
#   preset = "publication"
# )
# head(pathway_long[, c("Layer", "CurveGroup", "Theta", "Value")])

## ----custom-guidance--------------------------------------------------------
# names(wright_payload$data)
# wright_payload$data$reference_lines

## ----save-figure------------------------------------------------------------
# figure_file <- tempfile(fileext = ".rds")
# figure_record <- list(plot = wright,
#   rubric = c(`1` = "Beginning", `2` = "Developing", `3` = "Secure", `4` = "Advanced"))
# saveRDS(figure_record, figure_file)
# saved_record <- readRDS(figure_file)
# saved_figure <- saved_record$plot
# saved_record$rubric
# print(saved_figure)
# identical(plot_data(saved_figure, "locations"), plot_data(wright, "locations"))

## ----reopen-figure, eval = FALSE--------------------------------------------
# as_ggplot(saved_figure)

## ----session-preset---------------------------------------------------------
# local({
#   old <- options(mfrmr.plot_preset = "monochrome")
#   on.exit(options(old))
#   saved <- plot(fit, type = "ccc", draw = FALSE)
#   special <- plot(fit, type = "ccc", preset = "publication", draw = FALSE)
#   c(session_default = saved$data$preset, explicit_override = special$data$preset)
# })

## ----title-controls, fig.alt="Infit and Outfit mean squares for fitted facet levels, with the title omitted. Bubble sizes reflect rating counts; reference lines and the saved descriptive review band remain unchanged."----
# plain <- plot(fit, type = "bubble", diagnostics = diag, view = "infit_outfit", title = NULL, preset = "monochrome", draw = FALSE)
# as_ggplot(plain)
# # Reuse the saved values, screening references and interval settings.
# plain$data$reference_lines

## ----fair-score-data, eval=FALSE--------------------------------------------
# p_fair <- plot(fit, type = "fair_average", diagnostics = diag, facet = "Rater", metric = "FairZ", view = "measure", show_ci = TRUE, preset = "monochrome", show_title = FALSE, show_notes = FALSE, draw = FALSE)
# p_fair$data$notes
# p_fair$data$plot_data
# p_fair$data$excluded
# # Requires ggplot2:
# g_fair <- as_ggplot(p_fair)
# g_fair
# attr(g_fair, "mfrmr_notes")

## ----saved-agreement-view, fig.alt="Exact agreement for matched rater pairs with the fitted-model expected agreement overlaid. Review flags and the horizontal guide use the thresholds saved in the agreement result."----
# agreement <- interrater_agreement_table(fit, diagnostics = diag, rater_facet = "Rater")
# p_agreement <- plot(agreement)
# # Alternative saved view: plot(agreement, type = "agreement", metric = "correlation")
# # Optional ggplot conversion: as_ggplot(p_agreement)

## ----response-time-review, fig.alt = c("Distribution of synthetic response times with selected lower and upper quantile screens. These thresholds describe this example and do not diagnose careless responding.", "Person-level summaries of synthetic response times. Compare flagged patterns with their rating counts and context before drawing substantive conclusions.")----
# toy_rt <- toy
# toy_rt$ResponseTime <- 12 + (seq_len(nrow(toy_rt)) %% 7) +
#   as.numeric(toy_rt$Score)
# toy_rt$ResponseTime[1] <- 2
# toy_rt$ResponseTime[2] <- 38
# 
# rt <- response_time_review(
#   toy_rt,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   time = "ResponseTime",
#   rapid_quantile = 0.10,
#   slow_quantile = 0.90
# )
# 
# summary(rt)
# plot(rt, type = "distribution", preset = "publication")
# plot(rt, type = "person", preset = "publication")

## ----shrinkage-funnel, fig.alt = "Raw and empirically shrunken facet estimates with their displayed uncertainty. Movement toward a facet center reflects shrinkage, not a change in observed ratings."----
# fit_eb <- apply_empirical_bayes_shrinkage(fit)
# 
# shrink <- plot(fit_eb, type = "shrinkage", view = "funnel", show_ci = TRUE, ci_level = 0.95, preset = "publication", draw = FALSE)
# 
# head(shrink$data$table[, c(
#   "Facet", "Level", "RawEstimate", "RawCI_Lower", "RawCI_Upper",
#   "ShrunkEstimate", "ShrunkCI_Lower", "ShrunkCI_Upper",
#   "ShrinkageFactor"
# )])
# 
# plot(fit_eb, type = "shrinkage", view = "funnel", show_ci = TRUE, ci_level = 0.95, preset = "publication")

## ----feature-input----------------------------------------------------------
# backgrounds <- data.frame(
#   Rater = paste0("R", 1:8),
#   Years = c(1, 2, 4, 6, 8, 10, 12, 15),
#   Training = ordered(c("Basic", "Basic", "Advanced", "Basic",
#     "Advanced", "Advanced", "Specialist", "Specialist"),
#     levels = c("Basic", "Advanced", "Specialist")))
# features <- mfrm_features(backgrounds, "Rater", c("Years", "Training"))
# hierarchy <- mfrm_cluster_hierarchical(features, k = 3, linkage = "average")

## ----feature-hierarchy, dev="png", fig.alt="Average-linkage dendrogram for eight fictional raters, based on Gower dissimilarities in experience and ordered training level. Labelled leaves identify raters and boxes show the requested three groups; heights are dissimilarities, not significance levels."----
# hierarchy_plot <- plot(hierarchy, preset = "monochrome")

## ----hierarchy-values-------------------------------------------------------
# plot_data(hierarchy_plot, component = "table")
# tree <- plot_data(hierarchy_plot, component = "tree")

## ----hierarchy-ggplot, fig.width=8, fig.height=6, fig.alt="The same saved hierarchy of eight fictional raters rendered with ggplot. Leaf order, Gower merge heights and the three group boxes are unchanged. The caption explains that heights are not branch support."----
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   hierarchy_figure <- as_ggplot(hierarchy_plot, component = "tree") +
#     ggplot2::labs(title = NULL, subtitle = NULL)
#   print(hierarchy_figure)
# }

## ----d-study-input----------------------------------------------------------
# tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
# gstudy <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
# planned <- mfrm_d_study(gstudy,
#   data.frame(Tasks = c(3, 6, 12)), weights = c(V = 0.5, W = 0.5))

## ----d-study-planning, dev="png", fig.alt="G and Phi coefficients for an equally weighted V and W composite at three, six and twelve planned Tasks. Separate panels distinguish relative and absolute decisions; points identify requested designs, and lines guide the eye."----
# planning_plot <- plot(planned, type = "coefficients", preset = "monochrome")

## ----planning-values--------------------------------------------------------
# plot_data(planning_plot, component = "table")
# # Customizable conversion of this selected score/composite:
# # as_ggplot(planning_plot)
# # Error in score units is a different view:
# # plot(planned, type = "sem")

## ----eval = FALSE-----------------------------------------------------------
# summary(scores)
# plot(scores, type = "interval", preset = "publication")
# plot(scores, type = "precision", preset = "publication")
# plot(scores, type = "edge_mass", preset = "publication")

