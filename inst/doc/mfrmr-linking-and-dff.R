## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----anchor-data--------------------------------------------------------------
# library(mfrmr)
# ratings <- load_mfrmr_data("example_core")
# head(ratings[, c("Person", "Rater", "Criterion", "Score")])
# people <- unique(ratings$Person)
# baseline_data <- ratings[ratings$Person %in% people[1:24], , drop = FALSE]
# new_data <- ratings[ratings$Person %in% people[25:48], , drop = FALSE]

## ----anchor-reference---------------------------------------------------------
# baseline <- fit_mfrm(
#   baseline_data, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML", rating_min = 1, rating_max = 4,
#   category_policy = "preserve"
# )
# anchors <- make_anchor_table(baseline, facets = "Criterion")
# anchors

## ----anchor-review------------------------------------------------------------
# anchor_review <- review_mfrm_anchors(
#   new_data, "Person", c("Rater", "Criterion"), "Score", anchors = anchors,
#   rating_min = 1, rating_max = 4, category_policy = "preserve"
# )
# anchor_review$issue_counts
# anchor_review$anchors
# anchor_review$recommendations

## ----anchor-new-cohort--------------------------------------------------------
# linked <- anchor_to_baseline(
#   new_data, baseline, "Person", c("Rater", "Criterion"), "Score",
#   anchor_facets = "Criterion", anchor_policy = "error",
#   rating_min = 1, rating_max = 4, category_policy = "preserve"
# )
# summary(linked)
# linked$fit$facets$others
# linked$drift[, c("Facet", "Level", "Baseline", "New", "Drift")]

## ----migration-legacy---------------------------------------------------------
# migration_anchor <- data.frame(Facet = "Rater", Level = "R01", Anchor = 0.70)
# legacy_fit <- fit_mfrm(
#   ratings, "Person", c("Rater", "Criterion"), "Score",
#   model = "RSM", method = "MML", anchors = migration_anchor,
#   noncenter_facet = "Person", centering_scope = "remainder",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   quad_points = 31, maxit = 400, reltol = 1e-9
# )
# fit_file <- tempfile(fileext = ".rds")
# saveRDS(legacy_fit, fit_file)
# restored_fit <- readRDS(fit_file)
# unlink(fit_file)
# spec <- specifications_report(restored_fit)$data_spec
# spec[spec$Setting %in% c("NonCenteredFacet", "CenteringScope"), ]

## ----migration-whole-facet----------------------------------------------------
# updated_fit <- fit_mfrm(
#   ratings, "Person", c("Rater", "Criterion"), "Score",
#   model = "RSM", method = "MML", anchors = migration_anchor,
#   noncenter_facet = "Person", centering_scope = "facet",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   quad_points = 31, maxit = 400, reltol = 1e-9
# )
# summary(updated_fit, compute = "never", detail = "brief")
# old_raters <- subset(restored_fit$facets$others, Facet == "Rater")
# new_raters <- subset(updated_fit$facets$others, Facet == "Rater")
# data.frame(
#   Convention = c("remainder", "facet"),
#   FixedR01 = c(old_raters$Estimate[old_raters$Level == "R01"],
#                new_raters$Estimate[new_raters$Level == "R01"]),
#   OtherRaterSum = c(sum(old_raters$Estimate[old_raters$Level != "R01"]),
#                     sum(new_raters$Estimate[new_raters$Level != "R01"]))
# )

## ----setup--------------------------------------------------------------------
# library(mfrmr)
# 
# bias_df <- load_mfrmr_data("example_bias")
# 
# fit <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   rating_min = 1, rating_max = 4, category_policy = "preserve"
# )
# 
# diag <- diagnose_mfrm(fit, residual_pca = "none")

## ----connectivity, fig.alt = "Observed rating connections across design subsets. Disconnected blocks indicate where a common measurement scale needs additional linking evidence."----
# sc <- subset_connectivity_report(fit, diagnostics = diag)
# 
# sc$summary[, c("Subset", "Observations", "ObservationPercent")]
# plot(sc, type = "design_matrix", preset = "publication")

## ----anchors------------------------------------------------------------------
# anchors <- make_anchor_table(fit, facets = "Criterion")
# head(anchors)

## ----dff-residual, fig.alt = "Ordered criterion contrasts in mean observed-minus-expected scores, Group1 minus Group2. These descriptive differences are in score units, without confidence intervals or a formal DFF decision."----
# dff_resid <- analyze_dff(
#   fit,
#   diag,
#   facet = "Criterion",
#   group = "Group",
#   data = bias_df,
#   method = "residual"
# )
# 
# screen <- plot(dff_resid, type = "summary")
# comparisons <- plot_data(screen, "pairs")
# comparisons[, c("Level", "Group1", "Group2", "Contrast",
#                 "N_Group1", "N_Group2", "PlotStatus")]

## ----dff-response-support-----------------------------------------------------
# cases <- plot_data(screen, "cases")
# table(cases$Reason)
# cases[!cases$Included, ]
# plot_data(screen, "notes")

## ----dff-refit----------------------------------------------------------------
# dff_refit <- analyze_dff(
#   fit,
#   diag,
#   facet = "Criterion",
#   group = "Group",
#   data = bias_df,
#   method = "refit"
# )
# 
# dff_refit$summary
# dff_refit$dif_table[, c("Level", "Group1", "Group2", "Contrast",
#                        "ScaleLinkStatus", "Fixed1", "Fixed2")]
# dff_refit$dif_table[, c("Level", "LinkingAnchoredLevels1",
#                        "LinkingAnchoredLevels2", "LinkingThreshold")]

## ----school-comparison-review-------------------------------------------------
# school_ratings <- expand.grid(
#   Person = paste0("P", 1:8), Rater = paste0("R", 1:4),
#   Criterion = c("C1", "C2"), stringsAsFactors = FALSE
# )
# school_ratings$School <- ifelse(
#   school_ratings$Person %in% paste0("P", 1:4), "A", "B"
# )
# school_ratings$Score <- (as.integer(factor(school_ratings$Person)) +
#   as.integer(factor(school_ratings$Rater)) +
#   as.integer(factor(school_ratings$Criterion))) %% 3
# school_ratings <- subset(school_ratings, !(Rater == "R4" & School == "B"))
# 
# school_questions <- list(
#   facet = "Rater", group = "School", groups = c("A", "B"),
#   model = "RSM", method = "JML",
#   contrasts = rbind(
#     R1_vs_panel = c(R1 = .75, R2 = -.25, R3 = -.25, R4 = -.25),
#     R1_vs_R2 = c(R1 = 1, R2 = -1, R3 = 0, R4 = 0)
#   )
# )
# school_review <- describe_mfrm_data(
#   school_ratings, "Person", c("Rater", "Criterion", "School"), "Score",
#   rating_min = 0, rating_max = 2, category_policy = "preserve",
#   include_agreement = FALSE, comparisons = school_questions
# )
# school_review$comparison_review$table[, c(
#   "Comparison", "Group1", "Group2", "State", "N_Group1", "N_Group2"
# )]
# school_review$comparison_review$cells

## ----school-comparison-assumptions--------------------------------------------
# population_reviews <- lapply(
#   list(common_normal = ~1, school_means = ~School),
#   function(population) {
#     question <- school_questions
#     question$method <- "MML"
#     question$population_formula <- population
#     describe_mfrm_data(
#       school_ratings, "Person", c("Rater", "Criterion", "School"), "Score",
#       rating_min = 0, rating_max = 2, category_policy = "preserve",
#       include_agreement = FALSE, comparisons = question
#     )$comparison_review
#   }
# )
# lapply(population_reviews, function(review) review$table[, c(
#   "Comparison", "ConditionalState", "FreePersonState", "IdentificationBasis", "State"
# )])

## ----dff-comparison-review----------------------------------------------------
# rater_levels <- sort(unique(as.character(bias_df$Rater)))
# panel_weights <- matrix(-1 / length(rater_levels), 1, length(rater_levels),
#   dimnames = list("First rater versus panel", rater_levels))
# panel_weights[1, 1] <- panel_weights[1, 1] + 1
# prospective <- describe_mfrm_data(
#   bias_df, "Person", c("Rater", "Criterion", "Group"), "Score",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   include_agreement = FALSE,
#   comparisons = list(facet = "Rater", group = "Group", groups = c("A", "B"),
#     contrasts = panel_weights, model = "RSM", method = "MML")
# )
# rater_screen <- analyze_dff(
#   fit, facet = "Rater", group = "Group", data = bias_df,
#   comparison_review = prospective
# )
# summary(rater_screen)$comparison_review$table
# support <- plot(rater_screen, draw = FALSE)
# plot_data(support, "comparison_review")$weights
# report <- dif_report(rater_screen)
# report$comparison_review$table
# saveRDS(rater_screen, file.path(tempdir(), "rater-screen-with-review.rds"))

## ----joint-comparison-estimates-----------------------------------------------
# joint_questions <- prospective$comparison_review$specification
# joint_questions$population_formula <- ~Group
# pair_weights <- setNames(rep(0, length(rater_levels)), rater_levels)
# pair_weights[1:2] <- c(1, -1)
# joint_questions$contrasts <- rbind(
#   panel_weights, "First versus second rater" = pair_weights
# )
# joint_review <- describe_mfrm_data(
#   bias_df, "Person", c("Rater", "Criterion", "Group"), "Score",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   include_agreement = FALSE, comparisons = joint_questions
# )
# joint_fit <- do.call(fit_mfrm, c(
#   list(data = bias_df, quad_points = 63L),
#   joint_review$comparison_review$fit_arguments
# ))
# joint_results <- summary(joint_fit, comparison_review = joint_review)
# joint_comparisons <- joint_results$interaction_comparisons
# joint_comparisons[, c(
#   "Comparison", "Group1", "Group2", "Estimate", "EstimateStatus",
#   "State", "ConditionalState", "QuadratureState", "FormalInferenceEligible"
# )]
# attr(joint_comparisons, "fit_checks")$numerical
# attr(joint_comparisons, "comparison_review")$weights
# saveRDS(list(fit = joint_fit, review = joint_review, results = joint_results),
#   file.path(tempdir(), "joint-school-comparisons.rds"))

## ----partial-comparison-estimates---------------------------------------------
# incomplete_ratings <- bias_df[!(bias_df$Rater == tail(rater_levels, 1) &
#                                 bias_df$Group == "B"), ]
# partial_questions <- joint_questions
# partial_questions$method <- "JML"
# partial_questions$population_formula <- NULL
# partial_review <- describe_mfrm_data(
#   incomplete_ratings, "Person", c("Rater", "Criterion", "Group"), "Score",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   include_agreement = FALSE, comparisons = partial_questions
# )
# partial_fit <- do.call(fit_mfrm, c(
#   list(data = incomplete_ratings, fit_scope = "comparisons",
#        comparison_review = partial_review, interaction_policy = "silent"),
#   partial_review$comparison_review$fit_arguments
# ))
# partial_comparisons <- interaction_effect_table(partial_fit)
# partial_comparisons[, c(
#   "Comparison", "ConditionalState", "Estimate", "EstimateStatus",
#   "BoundaryState", "NumericalState", "EstimateReason"
# )]
# partial_fit$checks$numerical
# saveRDS(partial_fit, file.path(tempdir(), "partial-school-comparisons.rds"))

## ----partial-cml-comparisons--------------------------------------------------
# cml_data <- load_mfrmr_data("example_bias")
# cml_raters <- sort(unique(as.character(cml_data$Rater)))
# cml_ratings <- cml_data[!(cml_data$Rater == tail(cml_raters, 1) &
#                           cml_data$Group == "B"), ]
# cml_weights <- rbind(
#   panel = c(1, rep(0, length(cml_raters) - 1)) - 1 / length(cml_raters),
#   pair = c(1, -1, rep(0, length(cml_raters) - 2))
# )
# colnames(cml_weights) <- cml_raters
# cml_review <- describe_mfrm_data(
#   cml_ratings, "Person", c("Rater", "Criterion", "Group"), "Score",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   include_agreement = FALSE,
#   comparisons = list(facet = "Rater", group = "Group", groups = c("A", "B"),
#                      contrasts = cml_weights, model = "RSM", method = "CML")
# )
# cml_fit <- do.call(fit_mfrm, c(
#   list(data = cml_ratings, fit_scope = "comparisons",
#        comparison_review = cml_review, interaction_policy = "silent"),
#   cml_review$comparison_review$fit_arguments
# ))
# cml_table <- interaction_effect_table(cml_fit)
# cml_table[, c("Comparison", "ConditionalState", "Estimate", "SE", "Lower",
#               "Upper", "FormalInferenceEligible", "EstimateReason")]
# cml_fit$covariance
# saveRDS(cml_fit, file.path(tempdir(), "cml-school-comparisons.rds"))

## ----partial-cml-interval-display, fig.width=7, fig.height=4.2, fig.alt="CML rater-by-origin contrasts: the pair estimate has a pointwise 95 percent conditional-information interval; the unidentified panel row is labeled not estimated and has no point."----
# cml_intervals <- mfrm_facet_intervals(cml_fit)
# plot(cml_intervals)
# apa_table(cml_intervals, which = "intervals")$table
# saveRDS(cml_intervals, file.path(tempdir(), "cml-school-intervals.rds"))

## ----standardized-interaction-intervals, warning=FALSE------------------------
# ratings <- load_mfrmr_data("example_bias")
# raters <- sort(unique(as.character(ratings$Rater)))
# weights <- rbind(
#   panel = c(1, rep(0, length(raters) - 1)) - 1 / length(raters),
#   pair = c(1, -1, rep(0, length(raters) - 2))
# )
# colnames(weights) <- raters
# sd_review <- describe_mfrm_data(
#   ratings, "Person", c("Rater", "Criterion", "Group"), "Score",
#   rating_min = 1, rating_max = 4, category_policy = "preserve",
#   include_agreement = FALSE,
#   comparisons = list(facet = "Rater", group = "Group", groups = c("A", "B"),
#     contrasts = weights, model = "RSM", method = "MML", population_formula = ~1)
# )
# sd_fit <- do.call(fit_mfrm, c(
#   list(data = ratings, quad_points = 63, interaction_policy = "silent"),
#   sd_review$comparison_review$fit_arguments
# ))
# sd_intervals <- mfrm_facet_intervals(
#   sd_fit, comparison_review = sd_review, scale = "standardized"
# )
# sd_intervals$settings$population_sd
# summary(sd_intervals)[, c("Target", "NativeEstimate", "Estimate", "SE", "Lower", "Upper")]
# saveRDS(sd_intervals, file.path(tempdir(), "standardized-school-intervals.rds"))

## ----standardized-interaction-plot, fig.width=7, fig.height=4.2, fig.alt="Two rater-by-origin interaction contrasts with pointwise 95 percent model-information delta-method intervals in estimated population standard-deviation units; the zero line is a reference."----
# plot(sd_intervals)

## ----standardized-interaction-report------------------------------------------
# sd_results <- mfrm_results(
#   sd_fit, intervals = list(school = sd_intervals),
#   include = c("fit", "plots"), compute = "never"
# )
# sd_results$tables$facet_school_intervals[, c("Target", "Estimate", "ScaleReference")]
# saveRDS(sd_results, file.path(tempdir(), "standardized-school-results.rds"))

## ----dff-comparison-input, eval = !is_cran_check && requireNamespace("ordinal", quietly = TRUE)----
# library(mfrmr)
# set.seed(241010)
# group_order <- c("Reference", "Focal")
# person_group <- rep(group_order, each = 120)
# ability <- rnorm(240, mean = ifelse(person_group == "Focal", 0.5, -0.5))
# shared <- expand.grid(Person = seq_len(240), Item = paste0("I", 1:8),
#                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
# shared$Group <- person_group[shared$Person]
# difficulty <- setNames(seq(-1, 1, length.out = 8), paste0("I", 1:8))
# shared$Score <- rbinom(nrow(shared), 1,
#   plogis(ability[shared$Person] - difficulty[shared$Item]))
# stopifnot(!anyNA(shared), all(table(shared$Person, shared$Item) == 1L))
# shared$Total <- ave(shared$Score, shared$Person, FUN = sum)
# 
# shared_fit <- fit_mfrm(shared, person = "Person", facets = "Item",
#                        score = "Score", method = "JML", model = "RSM")
# residual <- analyze_dff(shared_fit, facet = "Item", group = "Group",
#   data = shared, focal = group_order[2], method = "residual")
# mh <- mantel_haenszel_table(shared, person = "Person", target = "Item",
#   score = "Score", group = "Group", strata = "Total", groups = group_order)
# ologit <- ordinal_logistic_table(shared, person = "Person", target = "Item",
#   score = "Score", group = "Group", matching = "Total",
#   groups = group_order, categories = 0:1)

## ----dff-comparison-audit, eval = !is_cran_check && requireNamespace("ordinal", quietly = TRUE)----
# case_audit <- rbind(
#   data.frame(Method = "Residual", Row = residual$cases$InputRow,
#     Included = residual$cases$Included, Reason = residual$cases$Reason),
#   data.frame(Method = "MH", mh$cases[c("Row", "Included", "Reason")]),
#   data.frame(Method = "Logistic", ologit$cases[c("Row", "Included", "Reason")])
# )
# case_audit$Target <- shared$Item[case_audit$Row]
# case_audit$Person <- shared$Person[case_audit$Row]
# case_audit$Group <- shared$Group[case_audit$Row]
# case_audit$Count <- 1L
# case_counts <- aggregate(Count ~ Target + Method + Group + Included + Reason,
#                          data = case_audit, FUN = sum)
# case_counts[case_counts$Target == "I2", ]

## ----dff-comparison-effects, eval = !is_cran_check && requireNamespace("ordinal", quietly = TRUE)----
# r <- residual$dif_table
# stopifnot(all(r$Group1 == group_order[1]), all(r$Group2 == group_order[2]))
# effects <- data.frame(Target = r$Level, ResidualDifference = r$Contrast,
#   ResidualStatus = ifelse(is.finite(r$Contrast), "descriptive", "unavailable"))
# m <- mh$table[c("Target", "OddsRatio", "Status")]
# names(m)[2:3] <- c("MH_ExactOR", "MH_Status")
# effects <- merge(effects, m, by = "Target", all = TRUE)
# mid <- subset(ologit$contrasts, Position == "OverlapMid",
#               select = c(Target, Matching, OddsRatio, Status))
# names(mid)[2:4] <- c("TotalAtMidpoint", "Logistic_MidOR", "Logistic_Status")
# effects <- merge(effects, mid, by = "Target", all = TRUE)
# effects[effects$Target == "I2", ]
# 
# # Preserve endpoint effects; a midpoint alone can conceal a crossing effect.
# subset(ologit$contrasts, Target == "I2")
# subset(mh$table, Target == "I2", select = c(Target, NReference, NFocal,
#   Lower, Upper, PValue, PAdjusted, Status, IntervalStatus))
# subset(ologit$table, Target == "I2", select = c(Target, Comparison,
#   DeltaNagelkerkeR2, PValue, PAdjusted, Status))

## ----dff-comparison-save, eval = FALSE----------------------------------------
# comparison <- list(data = shared, groups = group_order, categories = 0:1,
#   matching_rule = "Complete common eight-item total, including the target",
#   fit = shared_fit, residual = residual, mh = mh, ordinal = ologit,
#   effects = effects, case_audit = case_audit, case_counts = case_counts,
#   session = sessionInfo())
# saveRDS(comparison, "dff-comparison.rds")

## ----dff-follow-up------------------------------------------------------------
# dit <- dif_interaction_table(
#   fit,
#   diag,
#   facet = "Criterion",
#   group = "Group",
#   data = bias_df
# )
# 
# head(dit$table)
# 
# dr <- dif_report(dff_resid)
# cat(dr$narrative)

## ----exploratory-bias-report, eval = FALSE------------------------------------
# ratings <- load_mfrmr_data("example_bias")
# fit_add <- fit_mfrm(
#   ratings, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML", model = "RSM", quad_points = 31
# )
# screens <- estimate_all_bias(fit_add)
# report <- mfrm_report(screens, study_context = c(
#   question = "Which rater-criterion cells need substantive review?",
#   score_direction = "Higher coded scores indicate better performance"
# ))
# report$tables$overview
# head(report$tables$bias_cells)
# plot(screens)
# # Keep complete tables, including failed pairs and unavailable cells.
# export_summary_appendix(report, output_dir = "bias-review")

## ----model-estimated-interaction, eval = FALSE--------------------------------
# # Continue with ratings and fit_add from the preceding block.
# fit_interaction <- fit_mfrm(
#   ratings,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   facet_interactions = "Rater:Criterion",
#   quad_points = 31
# )
# 
# interaction_effect_table(fit_interaction)
# comparison <- compare_mfrm(Additive = fit_add, Interaction = fit_interaction,
#   labels = c("Additive", "Interaction"), nested = TRUE)
# comparison
# model_report <- mfrm_report(build_model_choice_review(
#   Additive = fit_add, Interaction = fit_interaction
# ))
# model_report$tables$comparison_table
# # The nested LRT remains in comparison$lrt; the review above reports ICs.

## ----drift-route, eval = FALSE------------------------------------------------
# declared_common_facets <- c("Criterion")
# 
# fit1 <- fit_mfrm(
#   wave1, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML"
# )
# 
# anchored <- anchor_to_baseline(
#   wave2,
#   fit1,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   anchor_facets = declared_common_facets
# )
# 
# fit2 <- fit_mfrm(
#   wave2, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML"
# )
# drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2),
#                            facets = declared_common_facets)
# plot(drift, type = "drift", preset = "publication")

## ----assignment-examples------------------------------------------------------
# library(mfrmr)
# 
# common <- build_mfrm_sim_spec(
#   n_person = 40, n_rater = 8, n_criterion = 3, raters_per_person = 2,
#   assignment = "sparse_linked",
#   sparse_controls = list(link_persons = 4, link_raters_per_person = 8,
#                          assignment_mode = "balanced")
# )
# rotating <- build_mfrm_sim_spec(
#   n_person = 40, n_rater = 8, n_criterion = 3, raters_per_person = 2,
#   assignment = "rotating"
# )
# random <- build_mfrm_sim_spec(
#   n_person = 40, n_rater = 8, n_criterion = 3, raters_per_person = 2,
#   assignment = "sparse_linked",
#   sparse_controls = list(link_persons = 0, assignment_mode = "random")
# )
# 
# # Use common, rotating, or random here to inspect the chosen assignment.
# ratings <- simulate_mfrm_data(sim_spec = common, seed = 27309)
# data_review <- describe_mfrm_data(
#   ratings, person = "Person", facets = c("Rater", "Criterion"),
#   score = "Score", include_agreement = FALSE
# )
# data_review$design_connectivity

