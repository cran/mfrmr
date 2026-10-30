## ----setup, include=FALSE---------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(eval = !is_cran_check, collapse = TRUE, comment = "#>", dev.args = list(bg = "white"))
library(mfrmr)

## ----incomplete-example-----------------------------------------------------
# roster <- expand.grid(Condition = "Null", Replicate = 1:4,
#                      Target = c("R1", "R2"), stringsAsFactors = FALSE)
# roster$Affected <- FALSE
# results <- roster[c("Condition", "Replicate", "Target")]
# results$Flag <- c(FALSE, TRUE, FALSE, NA, FALSE, NA, FALSE, NA)
# performance <- mfrm_screening_performance(roster, results,
#   rule = "Infit or Outfit outside [0.5, 1.5]; both raters; no exclusion/refit")
# summary(performance)$by_family

## ----incomplete-plot, fig.width=8, fig.height=4.8, fig.alt="One known positive family event among three available trials out of four planned. The conditional rate is one third; unresolved-outcome bounds span one quarter to one half and are not confidence intervals."----
# plot(performance)

## ----generate-and-screen----------------------------------------------------
# spec <- build_mfrm_sim_spec(n_person = 60, n_rater = 4, n_criterion = 3,
#                            raters_per_person = 2, assignment = "rotating")
# planned <- expand.grid(Condition = "Rotating null", Replicate = 1:6,
#                        Target = sprintf("R%02d", 1:4),
#                        stringsAsFactors = FALSE)
# planned$Affected <- FALSE
# observed <- planned[c("Condition", "Replicate", "Target")]
# observed$Flag <- NA
# observed$FitReady <- FALSE
# observed$Message <- ""
# observed$Infit <- observed$Outfit <- NA_real_
# observed$InfitZSTD <- observed$OutfitZSTD <- NA_real_
# 
# for (i in 1:6) {
#   rows <- which(observed$Replicate == i)
#   messages <- character()
#   tryCatch(withCallingHandlers({
#     ratings <- simulate_mfrm_data(sim_spec = spec, seed = 700 + i)
#     fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#                     model = "RSM", quad_points = 31, maxit = 200)
#     observed$FitReady[rows] <- isTRUE(fit$summary$InferenceReady[1])
#     diagnostic <- diagnose_mfrm(fit, residual_pca = "none")
#     raters <- diagnostic$measures[diagnostic$measures$Facet == "Rater", ]
#     raters <- raters[match(observed$Target[rows], raters$Level), ]
#     observed[rows, c("Infit", "Outfit", "InfitZSTD", "OutfitZSTD")] <-
#       raters[c("Infit", "Outfit", "InfitZSTD", "OutfitZSTD")]
#     outside <- function(x) ifelse(is.finite(x), x < .5 | x > 1.5, NA)
#     observed$Flag[rows] <- outside(raters$Infit) | outside(raters$Outfit)
#   }, warning = function(w) {
#     messages <<- c(messages, conditionMessage(w))
#     invokeRestart("muffleWarning")
#   }), error = function(e) {
#     messages <<- c(messages, conditionMessage(e))
#   })
#   observed$Message[rows] <- paste(unique(messages), collapse = " | ")
# }
# 
# screen <- mfrm_screening_performance(planned, observed,
#   rule = paste("Infit OR Outfit outside [0.5, 1.5]; all four observed raters;",
#                "computable descriptive screens retained; no exclusion/refit"))
# screen$by_target
# screen$by_family[screen$by_family$Targets > 0, ]
# unique(observed[c("Replicate", "FitReady", "Message")])

## ----threshold-sensitivity--------------------------------------------------
# bands <- data.frame(
#   Profile = c("0.5-1.5", "0.7-1.3", "0.8-1.2", "0.4-1.2"),
#   Lower = c(.5, .7, .8, .4), Upper = c(1.5, 1.3, 1.2, 1.2)
# )
# sensitivity <- mfrm_screening_sensitivity(planned, observed, bands,
#   rule = "Ordinary RSM MML; EAP plug-in residuals; all raters; no exclusion/refit")
# subset(sensitivity$by_family, Targets > 0)

## ----directional-tiles, fig.width=9, fig.height=3.6, fig.alt="Separate underfit and overfit family false-flag rates for declared threshold profiles. Each cell shows the conditional percentage and known/planned replication count; these six illustrative trials do not qualify screening accuracy."----
# plot(sensitivity, direction = c("underfit", "overfit"))

## ----directional-curves, fig.width=7, fig.height=4.5, fig.alt="Directional false-flag rates across discrete threshold profiles, with pointwise Monte Carlo intervals. Different symbols and line types distinguish underfit and overfit."----
# plot(sensitivity, style = "curves", direction = c("underfit", "overfit"))

## ----separate-zstd, eval=FALSE----------------------------------------------
# review <- fit_measures_table(fit, flag_basis = "mnsq")
# combined_review <- fit_measures_table(fit, flag_basis = "mnsq_or_zstd")
# combined_study <- mfrm_screening_sensitivity(planned, observed, bands,
#   rule = "Same saved residuals; engine ZSTD; no exclusion/refit", zstd_cut = 2)

## ----stress-overfit, echo=FALSE, out.width="100%", fig.alt="Overfit family false-flag rates in six null conditions across five threshold bands. Mean-square-only rates rise sharply as the lower bound increases. Adding engine ZSTD yields 44 to 100 percent false flags under the 0.5 to 1.5 band. Cell counts are known outcomes out of 100 planned trials; inference-review fits are retained."----
# knitr::include_graphics("figures/screening-overfit-stress.png")

## ----stress-detection, echo=FALSE, out.width="100%", fig.alt="R6 detection curves show the tradeoff between thresholds and departures, with pointwise Monte Carlo intervals. Uniform-category contamination produces upper-side flags; modal-category replacement produces lower-side flags. A 0.4 to 1.2 band detects 65 percent of the former and 18 percent of the latter."----
# knitr::include_graphics("figures/screening-detection-stress.png")

