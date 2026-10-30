## ----setup, include=FALSE-----------------------------------------------------
library(mfrmr)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 8,
  fig.height = 5, eval = isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false"))))

## ----observed-design----------------------------------------------------------
# ratings <- load_mfrmr_data("example_operational")
# if (requireNamespace("igraph", quietly = TRUE)) {
#   observed <- mfrm_network_analysis(data = ratings, person = "Person",
#     facets = "Rater", score = "Score")
#   observed$summary
#   head(observed$source_connectivity$row_usage)
# }

## ----planned-panels-----------------------------------------------------------
# assignment <- do.call(rbind, lapply(1:48, function(i) {
#   panel <- (i - 1L) %/% 16L
#   local <- (i - 1L) %% 4L
#   data.frame(Person = paste0("P", i),
#     Rater = paste0("R", panel * 4L + c(local + 1L, (local + 1L) %% 4L + 1L)))
# }))
# linked_assignment <- rbind(assignment, data.frame(Person = "P16", Rater = "R5"))
# if (requireNamespace("igraph", quietly = TRUE)) {
#   before <- mfrm_network_analysis(data = assignment, person = "Person",
#     facets = "Rater", basis = "planned")
#   after <- mfrm_network_analysis(data = linked_assignment, person = "Person",
#     facets = "Rater", basis = "planned")
#   rbind(Before = before$summary, Added_link = after$summary)
#   after$bridge_edges[, c("From", "To", "Weight")]
# }

## ----panel-figure, fig.height=4.8, fig.cap="Planned person-rater assignments. Columns represent three panels; circles are persons and squares are raters. The red line links P16 to R5. The third panel remains disconnected; the added single link is fragile."----
# if (requireNamespace("igraph", quietly = TRUE)) {
#   nodes <- after$source_connectivity$nodes
#   number <- as.integer(sub("^[PR]", "", nodes$Level))
#   is_person <- nodes$Facet == "Person"
#   block_size <- ifelse(is_person, 16L, 4L)
#   panel <- (number - 1L) %/% block_size
#   coords <- cbind(panel * 3 + ifelse(is_person, 0, 1),
#     1 - ((number - 1L) %% block_size) / (block_size - 1L))
#   rownames(coords) <- nodes$Node
#   edges <- after$source_connectivity$edges
#   extra <- edges$From == "Person:P16" & edges$To == "Rater:R5"
#   picture <- plot(after, type = "network", layout = coords,
#     graphical = list(vertex.shape = ifelse(is_person, "circle", "square"),
#       vertex.size = ifelse(is_person, 3.5, 7),
#       vertex.label = ifelse(!is_person | nodes$Level == "P16", nodes$Level, NA),
#       edge.color = ifelse(extra, "#B2182B", "#BBBBBB"),
#       edge.width = ifelse(extra, 2.5, .8),
#       vertex.label.dist = ifelse(is_person, 1.5, 0), asp = 0))
# }

## ----fit-for-score-networks---------------------------------------------------
# if (requireNamespace("igraph", quietly = TRUE)) {
#   fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#     method = "JML", maxit = 300)
#   diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
#   agreement <- rater_network_analysis(fit, diagnostics, mode = "agreement")
#   directional <- rater_network_analysis(fit, diagnostics,
#     mode = "severity_direction", score_diff_tolerance = 0, severity_continuity = .5)
# }

## ----agreement-direction, fig.width=8, fig.height=4, fig.cap="Observed agreement (left) and directional discrepancies (right), with the same rater positions. Edge widths represent the relevant proportions. In the right panel larger nodes have higher network severity indices, not greater statistical significance."----
# if (requireNamespace("igraph", quietly = TRUE)) {
#   raters <- sort(directional$node_metrics$Rater)
#   angle <- seq(0, 2*pi, length.out = length(raters) + 1L)[seq_along(raters)]
#   coords <- cbind(cos(angle), sin(angle)); rownames(coords) <- raters
#   old <- par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))
#   plot(agreement, type = "network", layout = coords,
#     graphical = list(vertex.size = 13))
#   sev <- directional$node_metrics$SeverityIndex
#   sizes <- rep(10, length(sev)); ok <- is.finite(sev)
#   if (any(ok)) sizes[ok] <- 10 + 8 * (sev[ok] - min(sev[ok])) / max(.01, diff(range(sev[ok])))
#   plot(directional, type = "network", layout = coords,
#     graphical = list(vertex.size = sizes, edge.curved = .2))
#   par(old)
#   comparison <- merge(directional$node_metrics[, c("Rater", "SeverityIndex")],
#     subset(diagnostics$measures, Facet == "Rater", select = c(Level, Estimate, SE)),
#     by.x = "Rater", by.y = "Level", all = TRUE)
#   comparison
# }

## ----signed-alpha-and-network-fit, fig.alt="Rater nodes joined by pairwise Krippendorff alpha edges. Width represents magnitude; negative edges, when present, are dashed. Unavailable pairs have no edge."----
# if (requireNamespace("igraph", quietly = TRUE)) {
#   signed <- rater_network_analysis(fit, diagnostics,
#     weight_metric = "KrippendorffAlpha", alpha_level = "interval")
#   signed$pair_metrics[, c("Rater1", "Rater2", "N", "KrippendorffAlpha", "AlphaStatus")]
#   plot(signed) # Negative agreement uses dashed edges; undefined pairs stay NA.
#   magnitude <- rater_network_analysis(fit, diagnostics,
#     mode = "severity_direction", weight_metric = "DirectionalMagnitude",
#     severity_continuity = 0)
#   magnitude$network_fit
# }

## ----halo-figure, fig.cap="Rater-by-criterion observed-score network. Nodes are grouped by rater (C: Content, L: Language, O: Organization); blue lines join criteria scored by the same rater and grey lines join different raters. Width shows correlation magnitude; dashed lines indicate negative correlations. Correlation does not identify a causal halo effect."----
# if (requireNamespace("igraph", quietly = TRUE)) {
#   halo <- rater_halo_network_analysis(fit, diagnostics, method = "spearman",
#     alpha = 1, min_pair_n = 5, positive_only = FALSE)
#   nodes <- halo$node_metrics
#   raters <- unique(nodes$Rater); criteria <- unique(nodes$Criterion)
#   center <- 2*pi*(match(nodes$Rater, raters) - 1)/length(raters)
#   corner <- 2*pi*(match(nodes$Criterion, criteria) - 1)/length(criteria)
#   coords <- cbind(3*cos(center) + .65*cos(corner), 3*sin(center) + .65*sin(corner))
#   rownames(coords) <- nodes$Node
#   within_rater <- halo$edge_metrics$EdgeType == "halo"
#   plot(halo, type = "network", layout = coords,
#     graphical = list(vertex.label = paste(nodes$Rater, substr(nodes$Criterion, 1, 1), sep = "-"),
#       vertex.size = 10, vertex.label.cex = .85,
#       edge.color = ifelse(within_rater, "#2166AC", "#D4D4D4"),
#       edge.width = .3 + 2.5 * abs(halo$edge_metrics$SignedWeight)))
#   halo$halo_summary_by_rater
# }

## ----halo-bootstrap, eval=FALSE-----------------------------------------------
# halo_ci <- rater_halo_network_analysis(fit, diagnostics, alpha = 1,
#   positive_only = FALSE, bootstrap = 999, seed = 47,
#   halo_reference = "different_criteria")
# halo_ci$halo_contrast
# halo_ci$bootstrap_draws # All attempts, including undefined fixed-pair draws.
# mfrm_report(halo_ci)

## ----ergm-example, eval=FALSE-------------------------------------------------
# # An explicit tie definition: pairwise interval alpha > 0.35.
# # Missing pairs are modeled only under an ignorable missingness assumption.
# model <- fit_rater_network(signed, ~ edges, threshold = 0.35,
#   missing = "model", seed = 48, gof = ~ degree,
#   gof_control = ergm::control.gof.ergm(nsim = 199, seed = 49))
# model$dyads # NA is unobserved; 0 is an observed non-tie.
# model$coefficients
# model$diagnostics
# plot(model, type = "gof")
# report <- mfrm_report(model)
# export_summary_appendix(report, output_dir = "network-report", prefix = "ergm")

## ----export-recipe, eval=FALSE------------------------------------------------
# write.csv(plot_data(picture, "nodes"), "design-nodes.csv", row.names = FALSE)
# write.csv(plot_data(picture, "source_edges"), "design-edges.csv", row.names = FALSE)
# saveRDS(list(network = after, plot = picture), "design-review.rds")
# # For a fitted clustering result:
# raw <- plot_data(groups, "observations", type = "profile", feature = "ExperienceYears")
# write.csv(raw, "feature-observations.csv", row.names = FALSE)

