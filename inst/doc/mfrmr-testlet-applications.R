## ----setup, include=FALSE---------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(eval = !is_cran_check, collapse = TRUE, comment = "#>", cache = FALSE,
                     dev.args = list(bg = "white"))
library(mfrmr)

## ----data-------------------------------------------------------------------
# set.seed(92327001)
# ratings <- expand.grid(
#   Criterion = paste0("C", 1:5), Task = paste0("T", 1:3),
#   Person = sprintf("P%03d", 1:80), stringsAsFactors = FALSE
# )
# ability <- rnorm(80)
# local_effect <- matrix(rnorm(80 * 3, sd = sqrt(.8)), 80, 3)
# p <- match(ratings$Person, sprintf("P%03d", 1:80))
# t <- match(ratings$Task, paste0("T", 1:3))
# eta <- ability[p] + local_effect[cbind(p, t)]
# w <- exp(cbind(0, eta + .6, 2 * eta))
# probability <- w / rowSums(w)
# u <- runif(nrow(ratings))
# ratings$Score <- as.integer(u > probability[, 1]) +
#   as.integer(u > rowSums(probability[, 1:2]))

## ----calibration------------------------------------------------------------
# testlet_fit <- fit_mfrm_testlet(
#   ratings, "Person", "Score", "Task", score_levels = 0:2,
#   quad_points = 61
# )
# independent_fit <- fit_mfrm_testlet(
#   ratings, "Person", "Score", "Task", score_levels = 0:2,
#   testlet_variance = 0, quad_points = 121
# )
# for (fit in list(testlet_fit, independent_fit)) {
#   stopifnot(fit$checks$NumericalReady, fit$checks$InformationPositive)
# }
# data.frame(
#   Model = c("Testlet", "Independent"),
#   AbilityVariance = c(testlet_fit$calibration$person_variance,
#                      independent_fit$calibration$person_variance),
#   LocalVariance = c(testlet_fit$calibration$variance, 0)
# )

## ----load-plan-functions, ref.label='plan-functions', include=FALSE---------
# normal_rule <- function(n) {
#   jacobi <- matrix(0, n, n)
#   jacobi[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
#   e <- eigen(jacobi + t(jacobi), symmetric = TRUE)
#   list(nodes = e$values, weights = e$vectors[1, ]^2)
# }
# 
# pattern_mass <- function(fit, membership, patterns, order) {
#   rule <- normal_rule(order)
#   theta <- fit$calibration$person_sd * rule$nodes
#   gamma <- sqrt(fit$calibration$variance) * rule$nodes
#   answer <- numeric(nrow(patterns))
#   for (h in seq_along(theta)) {
#     logw <- outer(theta[h] + gamma, 0:2) -
#       matrix(c(0, cumsum(fit$calibration$steps)), order, 3, byrow = TRUE)
#     prob <- exp(logw - apply(logw, 1, max))
#     prob <- prob / rowSums(prob)
#     conditional <- rep(1, nrow(patterns))
#     for (block in split(seq_along(membership), membership)) {
#       block_prob <- matrix(1, order, nrow(patterns))
#       for (j in block) {
#         block_prob <- block_prob * prob[, patterns[, j] + 1L, drop = FALSE]
#       }
#       conditional <- conditional * as.vector(crossprod(rule$weights, block_prob))
#     }
#     answer <- answer + rule$weights[h] * conditional
#   }
#   answer
# }
# 
# evaluate_parallel_plan <- function(fit, membership) {
#   stopifnot(inherits(fit, "mfrm_testlet"),
#     length(fit$input$columns$facets) == 0L,
#     identical(as.numeric(fit$input$score_levels), 0:2 + 0),
#     length(membership) == 4L, !anyNA(membership))
#   patterns <- as.matrix(expand.grid(rep(list(0:2), 4)))
#   data <- data.frame(Person = rep(sprintf("Y%02d", 1:81), each = 4),
#     Task = rep(membership, 81), Score = as.vector(t(patterns)))
#   scores <- predict(fit, data, quad_points = 61)
#   stopifnot(identical(scores$table$Person, sprintf("Y%02d", 1:81)),
#             all(scores$table$Status == "available_conditional"))
#   low <- pattern_mass(fit, membership, patterns, 61)
#   probability <- pattern_mass(fit, membership, patterns, 121)
#   tab <- scores$table
#   v <- fit$calibration$person_variance
#   expected_variance <- sum(probability * tab$ConditionalSD^2)
#   stopifnot(abs(sum(probability) - 1) < 1e-8,
#     max(abs(low - probability)) < 1e-8,
#     abs(sum(probability * tab$Estimate)) < 1e-7,
#     abs(sum(probability * (tab$ConditionalSD^2 + tab$Estimate^2)) - v) < 1e-7)
#   list(scores = scores, probability = probability, patterns = patterns,
#     membership = membership, summary = data.frame(
#       RMSPosteriorSD = sqrt(expected_variance),
#       ModelImpliedEAPReliability = 1 - expected_variance / v))
# }

## ----plan-comparison--------------------------------------------------------
# plans <- list(
#   "1 task x 4 criteria" = rep("T1", 4),
#   "2 tasks x 2 criteria" = rep(c("T1", "T2"), each = 2),
#   "4 tasks x 1 criterion" = paste0("T", 1:4)
# )
# fits <- list(Testlet = testlet_fit, Independent = independent_fit)
# plan_results <- list()
# for (model in names(fits)) for (plan in names(plans)) {
#   key <- paste(model, plan, sep = ": ")
#   plan_results[[key]] <- evaluate_parallel_plan(fits[[model]], plans[[plan]])
#   plan_results[[key]]$summary$Model <- model
#   plan_results[[key]]$summary$Plan <- plan
# }
# plan_summary <- do.call(rbind, lapply(plan_results, `[[`, "summary"))
# rownames(plan_summary) <- NULL
# knitr::kable(plan_summary[c("Model", "Plan", "RMSPosteriorSD",
#                            "ModelImpliedEAPReliability")], digits = 3)

## ----unequal-blocks---------------------------------------------------------
# base <- data.frame(
#   Task = c(rep("T1", 5), rep("T2", 2)), Score = rep(1L, 7),
#   Criterion = c(paste0("C", 1:5), paste0("C", 1:2))
# )
# dense <- sparse <- base
# dense$Score[1] <- 2L
# sparse$Score[6] <- 2L
# profiles <- rbind(
#   transform(base, Person = "Baseline"),
#   transform(dense, Person = "5-criterion task +1"),
#   transform(sparse, Person = "2-criterion task +1")
# )
# unequal_scores <- lapply(fits, predict, newdata = profiles, quad_points = 61)
# changes <- do.call(rbind, lapply(names(unequal_scores), function(model) {
#   tab <- unequal_scores[[model]]$table
#   stopifnot(all(tab$Status == "available_conditional"))
#   data.frame(Model = model, Profile = tab$Person, Estimate = tab$Estimate,
#     ChangeFromBaseline = tab$Estimate - tab$Estimate[tab$Person == "Baseline"])
# }))
# knitr::kable(changes, digits = 3)

## ----unequal-plot, fig.width=8, fig.height=5, fig.alt='Conditional Person estimates for the baseline and two one-point changes under the testlet model. The full table above gives the values. Whiskers hold calibration fixed and are not task-weight or difference intervals.'----
# plot(unequal_scores$Testlet, sort = "estimate", palette = "mono",
#      title = "Where an additional point occurs matters",
#      reference = NULL)

## ----save-------------------------------------------------------------------
# archive <- list(fits = fits, plans = plans, plan_results = plan_results,
#                 profiles = profiles, unequal_scores = unequal_scores)
# path <- tempfile("testlet-planning-", fileext = ".rds")
# saveRDS(archive, path)
# restored <- readRDS(path)
# stopifnot(identical(restored$unequal_scores$Testlet$table,
#                     unequal_scores$Testlet$table))

## ----plan-functions, eval=FALSE---------------------------------------------
# normal_rule <- function(n) {
#   jacobi <- matrix(0, n, n)
#   jacobi[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
#   e <- eigen(jacobi + t(jacobi), symmetric = TRUE)
#   list(nodes = e$values, weights = e$vectors[1, ]^2)
# }
# 
# pattern_mass <- function(fit, membership, patterns, order) {
#   rule <- normal_rule(order)
#   theta <- fit$calibration$person_sd * rule$nodes
#   gamma <- sqrt(fit$calibration$variance) * rule$nodes
#   answer <- numeric(nrow(patterns))
#   for (h in seq_along(theta)) {
#     logw <- outer(theta[h] + gamma, 0:2) -
#       matrix(c(0, cumsum(fit$calibration$steps)), order, 3, byrow = TRUE)
#     prob <- exp(logw - apply(logw, 1, max))
#     prob <- prob / rowSums(prob)
#     conditional <- rep(1, nrow(patterns))
#     for (block in split(seq_along(membership), membership)) {
#       block_prob <- matrix(1, order, nrow(patterns))
#       for (j in block) {
#         block_prob <- block_prob * prob[, patterns[, j] + 1L, drop = FALSE]
#       }
#       conditional <- conditional * as.vector(crossprod(rule$weights, block_prob))
#     }
#     answer <- answer + rule$weights[h] * conditional
#   }
#   answer
# }
# 
# evaluate_parallel_plan <- function(fit, membership) {
#   stopifnot(inherits(fit, "mfrm_testlet"),
#     length(fit$input$columns$facets) == 0L,
#     identical(as.numeric(fit$input$score_levels), 0:2 + 0),
#     length(membership) == 4L, !anyNA(membership))
#   patterns <- as.matrix(expand.grid(rep(list(0:2), 4)))
#   data <- data.frame(Person = rep(sprintf("Y%02d", 1:81), each = 4),
#     Task = rep(membership, 81), Score = as.vector(t(patterns)))
#   scores <- predict(fit, data, quad_points = 61)
#   stopifnot(identical(scores$table$Person, sprintf("Y%02d", 1:81)),
#             all(scores$table$Status == "available_conditional"))
#   low <- pattern_mass(fit, membership, patterns, 61)
#   probability <- pattern_mass(fit, membership, patterns, 121)
#   tab <- scores$table
#   v <- fit$calibration$person_variance
#   expected_variance <- sum(probability * tab$ConditionalSD^2)
#   stopifnot(abs(sum(probability) - 1) < 1e-8,
#     max(abs(low - probability)) < 1e-8,
#     abs(sum(probability * tab$Estimate)) < 1e-7,
#     abs(sum(probability * (tab$ConditionalSD^2 + tab$Estimate^2)) - v) < 1e-7)
#   list(scores = scores, probability = probability, patterns = patterns,
#     membership = membership, summary = data.frame(
#       RMSPosteriorSD = sqrt(expected_variance),
#       ModelImpliedEAPReliability = 1 - expected_variance / v))
# }

