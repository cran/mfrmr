## ----include = FALSE--------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----declared-context, eval=FALSE-------------------------------------------
# context <- c(rating_unit = "One performance by judge by criterion rating",
#   score_direction = "Higher Comprehensibility means harder to understand",
#   rubric = "1 = very easy; 9 = very difficult", transformation = "None")
# report <- mfrm_report(results, study_context = context)
# report$tables$study_context

## ----gpcm-template, eval=FALSE----------------------------------------------
# template <- system.file("templates", "gpcm-report.Rmd", package = "mfrmr")
# file.copy(template, "gpcm-report.Rmd")
# rmarkdown::render("gpcm-report.Rmd")
# 
# # Then use your own data: Person, Rater, Criterion, Score; one rating per cell.
# rmarkdown::render("gpcm-report.Rmd", params = list(
#   ratings_csv = "ratings.csv", rating_min = 1, rating_max = 4,
#   output_dir = "my-gpcm-report-files"
# ))

## ----setup------------------------------------------------------------------
# # Load the package and synthetic ratings
# library(mfrmr)
# toy <- load_mfrmr_data("example_operational")
# head(toy)
# 
# # Estimate person ability, rater severity, and criterion difficulty
# fit <- fit_mfrm(
#   data = toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM"
# )
# 
# # Reuse these diagnostics in the tables and report below
# diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
# results <- summary(fit, diagnostics = diag)

## ----analysis-record--------------------------------------------------------
# results$overview[, c("Model", "Method", "N", "Persons", "Converged")]
# results$settings_overview[, c(
#   "QuadPoints", "RatingMin", "RatingMax", "RatingRangeSource", "NoncenterFacet"
# )]
# results$row_retention
# results$decision
# 
# # Record the iteration limit, integration grid, and numerical stopping details
# spec <- specifications_report(fit)
# spec$convergence_control
# results$overview[, c("OptimizerMethod", "ConvergenceStatus", "RequestedReltol", "EffectiveReltol")]

## ----rating-design----------------------------------------------------------
# data("mfrmr_example_operational_design", package = "mfrmr")
# data_review <- describe_mfrm_data(
#   data = toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   rating_min = 1,
#   rating_max = 4,
#   keep_original = TRUE,
#   expected_design = mfrmr_example_operational_design
# )
# data_summary <- summary(data_review)
# data_summary$missing
# data_summary$structural_missingness
# data_summary$design_connectivity
# 
# # Describe the scores and the amount of rating evidence
# data_summary$score_distribution
# data_summary$facet_overview
# data.frame(MeanScore = mean(toy$Score), SDScore = sd(toy$Score))
# table(table(toy$Person)) # Number of persons with each observed rating count
# as.data.frame(results$person_overview[, c(
#   "Persons", "DistributionN", "ReviewExcludedExtremeEAPs", "Mean", "SD", "MeanPosteriorSD"
# )])

## ----precision--------------------------------------------------------------
# prec <- precision_review_report(fit, diagnostics = diag)
# prec$profile[, c(
#   "PrecisionTier", "SupportsFormalInference", "PersonSEBasis", "NonPersonSEBasis"
# )]
# prec$checks
# 
# # These diagnostics address departures from the fitted response model
# diag_summary <- summary(diag)
# diag_summary$overall_fit
# diag_summary$diagnostic_basis[, c("DiagnosticPath", "Status", "ReportingUse")]

## ----residual-counts--------------------------------------------------------
# z <- diag$obs$StdResidual
# residual_counts <- data.frame(
#   AbsoluteZCutoff = c(2, 3),
#   Evaluated = sum(is.finite(z)),
#   Unavailable = sum(!is.finite(z)),
#   Flagged = c(sum(abs(z[is.finite(z)]) >= 2),
#               sum(abs(z[is.finite(z)]) >= 3))
# )
# residual_counts$Percent <- 100 * residual_counts$Flagged / residual_counts$Evaluated
# residual_counts

## ----marginal-reporting-----------------------------------------------------
# diag_summary$marginal_fit[, c(
#   "OverallRMSD", "StepGroupsFlagged", "FacetLevelsFlagged",
#   "PairwiseFlaggedLevelPairs", "ReportingUse"
# )]
# diag$marginal_fit$thresholds
# diag$marginal_fit$pairwise$thresholds
# subset(diag$marginal_fit$facet_level$summary_stats, Flagged,
#        c("Facet", "Level", "GroupCount", "RMSD", "MaxAbsStdResidual"))
# subset(diag$marginal_fit$pairwise$pair_stats, Flagged,
#        c("Facet", "Level1", "Level2", "LevelPairCount", "ExactAgreement", "ExpectedExactAgreement"))

## ----rater-estimates, fig.alt = "Person abilities, facet locations and shared category steps on a common logit scale. Facet whiskers show 95% normal intervals from the fitted MML model; they do not test pairwise rater differences."----
# measurements <- fit_measures_table(fit, diagnostics = diag, ci_level = 0.95)
# rater_table <- subset(
#   measurements$table, Facet == "Rater",
#   c("Level", "N", "Measure", "SE", "CI_Lower", "CI_Upper", "Infit", "Outfit")
# )
# rater_table <- rater_table[order(rater_table$Measure), ]
# 
# rater_apa <- apa_table(
#   rater_table,
#   digits = 3,
#   caption = "Rater severity and fit in the synthetic MML/RSM example",
#   note = paste(
#     "282 ratings of 48 persons by six raters on three criteria using a 1-4 rubric.",
#     "N counts rating rows. Measures and SEs are in logits; higher values mean stricter ratings.",
#     "Intervals are 95% normal approximations from observed-information SEs.",
#     "Rater effects are centered within facet. Infit and Outfit have reference value 1."
#   )
# )
# rater_apa
# 
# # Display the estimates on the common logit scale
# plot(fit, diagnostics = diag, show_ci = TRUE)

## ----reliability-and-agreement----------------------------------------------
# apa_table(
#   diag_summary$reliability[, c("Facet", "Levels", "Separation", "Strata", "Reliability", "RealReliability")],
#   digits = 3,
#   caption = "Separation and reliability by facet",
#   note = "High Person separation reliability indicates distinguishable abilities; high rater separation reliability indicates distinguishable severities, not agreement. Low rater values can also reflect imprecision."
# )
# 
# agreement <- interrater_agreement_table(fit, diagnostics = diag)
# agreement$summary[, c(
#   "OpportunityCount", "ExactAgreement", "ExpectedExactAgreement", "AdjacentAgreement"
# )]

## ----reliability-basis------------------------------------------------------
# diag$reliability[, c("Facet", "Levels", "ObservedVariance", "ModelErrorVariance", "Reliability")]

## ----categories, fig.alt = "Category-use and threshold summaries for the synthetic four-category rubric. Review ordering and sparse use together rather than treating a single display as validation of the scale."----
# scale <- rating_scale_table(fit, diagnostics = diag)
# scale$category_table[, c("Category", "Count", "AvgMeasure", "ExpectedMeasure", "Infit", "Outfit")]
# plot(scale)
# scale$notes

## ----category-owner---------------------------------------------------------
# scale$model_context
# scale$owner_category_table[, c(
#   "StepFacet", "CategoryLabel", "RawCount", "Count", "Percent", "UsageStatus"
# )]

## ----threshold-uncertainty--------------------------------------------------
# steps <- scale$threshold_table
# steps[, c("StepFacet", "Step", "Estimate", "SE", "CI_Lower", "CI_Upper")]
# steps[, c("StepFacet", "Step", "SE_Status", "CIEligible", "CIUse",
#           "UncertaintySourceStatus")]

## ----observed-and-adjusted--------------------------------------------------
# fair <- fair_average_table(fit, diagnostics = diag, reference = "zero")
# person_scores <- fair$raw_by_facet$Person
# person_scores <- person_scores[order(person_scores$Level), ]
# head(person_scores[, c(
#   "Level", "TotalCount", "ObservedAverage", "FairZ", "Measure", "ModelSE"
# )], 6)

## ----decision-plan----------------------------------------------------------
# decision_plan <- expand.grid(
#   Rater = sort(unique(as.character(toy$Rater)))[1:2],
#   Criterion = sort(unique(as.character(toy$Criterion))),
#   stringsAsFactors = FALSE
# )
# decisions <- classification_table(
#   fit, rating_plan = decision_plan,
#   cut_scores = c(Pass = 0), score_cuts = 9, operational_cuts = 0,
#   theta = c(-1, -0.25, 0, 0.25, 1),
#   theta_weights = c(0.10, 0.20, 0.30, 0.25, 0.15)
# )
# decisions$cuts
# decisions$table[, c("Theta", "Method", "Accuracy", "Consistency", "SEM", "Status", "SourceStatus")]
# decisions$marginal

## ----decision-errors--------------------------------------------------------
# error_rates <- do.call(rbind, lapply(
#   split(decisions$classification_matrix, decisions$classification_matrix$Method),
#   function(d) {
#     false_pass <- sum(d$Probability[d$TrueCategory == 1 & d$ObservedCategory == 2])
#     false_fail <- sum(d$Probability[d$TrueCategory == 2 & d$ObservedCategory == 1])
#     joint <- c(false_pass, false_fail)
#     true_group <- c(sum(d$Probability[d$TrueCategory == 1]),
#                     sum(d$Probability[d$TrueCategory == 2]))
#     decision_group <- c(sum(d$Probability[d$ObservedCategory == 2]),
#                         sum(d$Probability[d$ObservedCategory == 1]))
#     # A missing or zero conditioning mass gives an undefined rate, not zero.
#     given_true <- given_decision <- rep(NA_real_, 2)
#     ok <- is.finite(true_group) & true_group > 0
#     given_true[ok] <- joint[ok] / true_group[ok]
#     ok <- is.finite(decision_group) & decision_group > 0
#     given_decision[ok] <- joint[ok] / decision_group[ok]
#     data.frame(
#       Method = d$Method[1], Error = c("False pass", "False fail"),
#       JointProbability = joint, TrueGroupProbability = true_group,
#       RateGivenTrueGroup = given_true, DecisionGroupProbability = decision_group,
#       RateGivenDecisionGroup = given_decision,
#       Status = d$Status[1], SourceStatus = decisions$source$status
#     )
#   }
# ))
# rownames(error_rates) <- NULL
# error_rates

## ----checklist--------------------------------------------------------------
# chk <- reporting_checklist(fit, diagnostics = diag, include_references = TRUE)
# chk$section_summary
# subset(chk$checklist, !DraftReady, c("Section", "Item", "NextAction"))
# 
# res <- mfrm_results(fit, diagnostics = diag, include = "publication")
# report <- mfrm_report(res, style = "apa")
# report$first_screen[, c("Area", "Status", "MainIssue", "NextAction")]

## ----apa--------------------------------------------------------------------
# apa <- build_apa_outputs(
#   fit,
#   diagnostics = diag,
#   context = list(
#     assessment = "Synthetic writing ratings",
#     setting = "Illustration of the reporting workflow",
#     scale_desc = "1-4 rubric scale",
#     rater_facet = "Rater"
#   )
# )
# cat(apa$report_text)

## ----export-----------------------------------------------------------------
# # A new temporary folder makes this example repeatable.
# # For your study, choose a new permanent project folder instead.
# output_dir <- tempfile("mfrmr-report-")
# bundle <- export_mfrm_bundle(
#   fit,
#   diagnostics = diag,
#   data = toy,
#   output_dir = output_dir,
#   prefix = "analysis01",
#   include = c(
#     "core_tables", "checklist", "dashboard", "apa",
#     "summary_tables", "manifest", "script", "html"
#   ),
#   acknowledge_sensitive = TRUE
# )
# selected_files <- subset(
#   bundle$written_files,
#   Component %in% c("core_measures", "replay_script", "bundle_html")
# )
# data.frame(
#   Component = selected_files$Component,
#   File = basename(selected_files$Path)
# )

## ----manuscript-files-------------------------------------------------------
# # Keep full precision in CSV; the displayed table is rounded to three decimals
# write.csv(rater_table, file.path(output_dir, "Table1-raters.csv"), row.names = FALSE)
# writeLines(capture.output(print(rater_apa)), file.path(output_dir, "Table1-raters.txt"))
# write.csv(steps, file.path(output_dir, "Step-uncertainty.csv"), row.names = FALSE)
# write.csv(scale$owner_category_table, file.path(output_dir, "Category-owner-counts.csv"), row.names = FALSE)
# writeLines(scale$notes, file.path(output_dir, "Category-notes.txt"))
# saveRDS(list(classification = decisions, error_rates = error_rates), file.path(output_dir, "Decision-analysis.rds"))
# write.csv(error_rates, file.path(output_dir, "Decision-error-rates.csv"), row.names = FALSE)
# writeLines(c(decisions$notes,
#   "False pass conditions on true failure or declared pass; false fail on true pass or declared failure.",
#   "JointProbability uses the whole supplied population; missing/zero conditioning mass gives an undefined conditional rate.",
#   "Use Decision-analysis.rds for the rating plan, category map, cuts, population weights and source assessment."),
#   file.path(output_dir, "Decision-notes.txt"))
# write.csv(mfrmr_example_operational_design, file.path(output_dir, "Planned-ratings.csv"), row.names = FALSE)
# write.csv(data_summary$structural_missingness, file.path(output_dir, "Planned-missingness.csv"), row.names = FALSE)
# writeLines(apa$report_text, file.path(output_dir, "Methods-Results-draft.txt"))
# 
# # Draw to the file device, then close it even if drawing fails
# png(file.path(output_dir, "Figure1-Wright.png"), width = 1800, height = 1500, res = 200)
# tryCatch(
#   figure <- plot(fit, diagnostics = diag, show_ci = TRUE),
#   finally = dev.off()
# )
# write.csv(plot_data(figure, "notes"),
#   file.path(output_dir, "Figure1-Wright-notes.csv"), row.names = FALSE)
# writeLines(c(
#   "Figure 1. Person, facet, and step locations in the synthetic MML/RSM analysis.",
#   "282 ratings; 48 persons, six raters, three criteria; rubric 1-4.",
#   "The histogram shows person EAP estimates. All locations are in logits.",
#   "Higher rater/criterion values mean stricter scoring/greater difficulty; effects are centered within facet.",
#   "Facet whiskers are 95% normal intervals using MML observed-information SEs from the same diagnostics as Table 1.",
#   "Steps are shared adjacent-category transitions shown as point estimates; see Step-uncertainty.csv for their intervals.",
#   "The display does not test rater contrasts or establish unidimensionality."
# ), file.path(output_dir, "Figure1-Wright-caption.txt"))
# 
# # Filenames are shown here; output_dir retains the folder location.
# list.files(output_dir, pattern = "^(Table1|Step-uncertainty|Category|Decision|Planned|Figure1|Methods-Results)")

## ----citations--------------------------------------------------------------
# # Get the citation for the package version actually used, and for R
# citation("mfrmr")
# citation()
# 
# # Methodological reading prompts, not complete bibliography entries
# chk$references
# 
# # Retain the computing environment with the analysis record
# sessionInfo()

