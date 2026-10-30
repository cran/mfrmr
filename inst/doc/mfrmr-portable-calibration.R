## ----include = FALSE--------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  eval = !is_cran_check
)
options(width = 78)

## ----capabilities-----------------------------------------------------------
# library(mfrmr)
# 
# mfrm_calibration_capabilities()[, c(
#   "Model", "Estimator", "ScoringBasis", "PortableCalibration"
# )]

## ----create-calibration-----------------------------------------------------
# synthetic <- load_mfrmr_data("example_core")
# person_ids <- unique(as.character(synthetic$Person))
# training_ids <- person_ids[seq_len(18L)]
# training <- synthetic[
#   as.character(synthetic$Person) %in% training_ids,
#   ,
#   drop = FALSE
# ]
# 
# fit <- fit_mfrm(
#   training,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   model = "RSM",
#   method = "MML",
#   quad_points = 31
# )
# 
# quadrature_review <- mml_quadrature_sensitivity(
#   fit,
#   training,
#   quad_points = c(31, 61),
#   theta_points = 41
# )
# summary(quadrature_review)
# fit <- quadrature_review$fits$q61
# 
# draft <- extract_mfrm_calibration(
#   fit,
#   calibration_id = "portable-rsm-example",
#   source_fit_id = "synthetic-rsm-fit",
#   scoring_quad_points = 121,
#   quadrature_review = quadrature_review
# )
# review_mfrm_calibration(draft)
# 
# validated <- validate_mfrm_calibration(draft)
# frozen <- freeze_mfrm_calibration(validated)
# frozen
# 
# artifact_file <- tempfile(fileext = ".rds")
# save_mfrm_calibration(frozen, artifact_file)

## ----score-from-artifact, fig.alt = "Posterior EAP intervals for two newly scored Persons. A blue circle marks a conditionally ready score, an orange triangle marks a score requiring review, and a dashed vertical line marks the zero-logit scale origin."----
# new_rows <- synthetic[
#   as.character(synthetic$Person) %in% person_ids[19:20],
#   c("Person", "Rater", "Criterion", "Score"),
#   drop = FALSE
# ]
# new_ids <- stats::setNames(
#   c("NEW_PERSON_1", "NEW_PERSON_2"), person_ids[19:20]
# )
# new_rows$Person <- unname(new_ids[as.character(new_rows$Person)])
# rownames(new_rows) <- NULL
# 
# # Create one explicit review case for the tutorial figure. In operational
# # scoring, never alter observed responses to manufacture a disposition.
# new_rows$Score[new_rows$Person == "NEW_PERSON_2"] <- max(training$Score)
# 
# rm(fit, quadrature_review, draft, validated, frozen, training, synthetic)
# invisible(gc())
# 
# calibration <- load_mfrm_calibration(artifact_file)
# scores <- score_mfrm_calibration(
#   calibration, new_rows, adaptive_quad_points = c(31, 61),
#   cut_scores = c(Basic = -0.5, Advanced = 0.8)
# )
# summary(scores)
# 
# plot(
#   scores,
#   type = "interval",
#   preset = "publication"
# )
# 
# unlink(artifact_file)

## ----score-review-payload, eval = FALSE-------------------------------------
# # Valid response rows versus posterior SD.
# plot(scores, type = "precision", preset = "publication")
# 
# # Posterior mass on the outer quadrature nodes versus the recorded threshold.
# edge_review <- plot(scores, type = "edge_mass", draw = FALSE)
# edge_review$data$selection_summary
# edge_review$data$unplotted_dispositions

## ----conditional-decisions--------------------------------------------------
# decision_table <- apa_table(scores, which = "decisions")
# decision_table$table[, c("Person", "Type", "Label", "Lower", "Upper",
#                          "Probability", "Status", "Reason")]
# 
# saved_scores <- tempfile(fileext = ".rds")
# saveRDS(scores, saved_scores)
# restored_scores <- readRDS(saved_scores)
# # The appendix retains estimates, review dispositions, cuts and exact probabilities.
# appendix_dir <- tempfile()
# appendix <- export_summary_appendix(restored_scores, output_dir = appendix_dir,
#                                    include_html = FALSE)
# unlink(appendix_dir, recursive = TRUE)
# unlink(saved_scores)

## ----score-integration-review-----------------------------------------------
# summary(scores)$quadrature_overview
# # Keep full precision when inspecting individual differences.
# head(scores$quadrature_review[c(
#   "Person", "AdaptiveNodes", "EAPChange", "PosteriorSDChange",
#   "AdaptiveEAPChangeFromPrevious", "AdaptiveSDChangeFromPrevious", "Status"
# )])

## ----portable-gpcm-fit, eval=FALSE------------------------------------------
# library(mfrmr)
# calibration_responses <- read.csv(
#   "calibration-responses.csv", colClasses = "character", na.strings = "",
#   check.names = FALSE, fileEncoding = "UTF-8-BOM"
# )
# gpcm_fit <- fit_mfrm(
#   calibration_responses,
#   person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#   model = "GPCM", method = "MML",
#   step_facet = "Rater", slope_facet = "Criterion",
#   rating_min = 0, rating_max = 2, category_policy = "preserve",
#   mml_integration = "adaptive", quad_points = 31,
#   maxit = 400, reltol = 1e-10
# )
# summary(gpcm_fit, profile = "fit", detail = "brief")

## ----portable-gpcm, eval=FALSE----------------------------------------------
# draft <- extract_mfrm_calibration(
#   gpcm_fit, calibration_id = "rubric-calibration",
#   source_fit_id = "training-cohort", scoring_quad_points = 141
# )
# review_mfrm_calibration(draft)
# calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
# save_mfrm_calibration(calibration, "gpcm-calibration.rds")

## ----portable-gpcm-score, eval=FALSE----------------------------------------
# library(mfrmr)
# # A later session needs only this file and the new response rows.
# calibration <- load_mfrm_calibration("gpcm-calibration.rds")
# new_responses <- read.csv(
#   "new-responses.csv", colClasses = "character", na.strings = "",
#   check.names = FALSE, fileEncoding = "UTF-8-BOM"
# )
# scores <- score_mfrm_calibration(
#   calibration, new_responses, missing_response = "omit"
# )
# score_review <- summary(scores)
# score_review
# score_review$review
# score_review$settings$score_integration_review
# plot(scores, type = "interval", preset = "publication", top_n = Inf)
# 
# # Keep full-precision scores and every Person's disposition, including those
# # without scores. The complete object supports later review and plotting.
# write.csv(scores$estimates, "person-scores.csv", row.names = FALSE, na = "")
# write.csv(scores$person_dispositions, "person-review.csv", row.names = FALSE, na = "")
# saveRDS(scores, "scored-cohort.rds")

## ----portable-gpcm-prior, eval=FALSE----------------------------------------
# retained <- scores$settings$retained_prior_identity
# 
# # Sensitivity analysis on the same calibration scale, not a new population fit.
# alternative <- score_mfrm_calibration(
#   calibration, new_responses, missing_response = "omit",
#   scoring_prior = list(mean = retained$mean + 0.5, sd = retained$sd * 1.25)
# )
# summary(alternative)

## ----portable-jml, eval=FALSE-----------------------------------------------
# library(mfrmr)
# jml_data <- load_mfrmr_data("example_core")
# jml_data <- jml_data[jml_data$Person %in% unique(jml_data$Person)[1:14], ]
# jml_fit <- fit_mfrm(jml_data, "Person", c("Rater", "Criterion"), "Score",
#   model = "PCM", method = "JML", step_facet = "Criterion",
#   maxit = 500, reltol = 1e-10)
# jml_draft <- extract_mfrm_calibration(jml_fit, scoring_quad_points = 141)
# jml_calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(jml_draft))
# save_mfrm_calibration(jml_calibration, "jml-calibration.rds")
# 
# # In a new session, supply your later cohort's ratings instead of these example rows.
# jml_new <- jml_data[jml_data$Person == unique(jml_data$Person)[1], ]
# jml_new$Person <- "NEXT01"
# jml_saved <- load_mfrm_calibration("jml-calibration.rds")
# jml_scores <- score_mfrm_calibration(jml_saved, jml_new)
# summary(jml_scores)
# plot(jml_scores, main = "New-cohort EAP scores")
# 
# jml_alternative <- score_mfrm_calibration(jml_saved, jml_new,
#   scoring_prior = list(mean = 0.25, sd = 1.2))
# jml_scores$estimates[, c("Person", "Estimate", "SD", "PriorMean", "PriorSD")]
# jml_alternative$estimates[, c("Person", "Estimate", "SD", "PriorMean", "PriorSD")]

## ----portable-gpcm-jml, eval=FALSE------------------------------------------
# gpcm_jml_fit <- fit_mfrm(jml_data, "Person", c("Rater", "Criterion"), "Score",
#   model = "GPCM", method = "JML", step_facet = "Criterion", slope_facet = "Criterion",
#   maxit = 600, reltol = 1e-11)
# gpcm_jml_calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(
#   extract_mfrm_calibration(gpcm_jml_fit, scoring_quad_points = 141)))
# save_mfrm_calibration(gpcm_jml_calibration, "gpcm-jml-calibration.rds")
# gpcm_jml_scores <- score_mfrm_calibration(
#   load_mfrm_calibration("gpcm-jml-calibration.rds"), jml_new)
# summary(gpcm_jml_scores)
# plot(gpcm_jml_scores, main = "New-cohort GPCM JML EAP scores")

