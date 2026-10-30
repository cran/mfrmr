## ----setup, include=FALSE-----------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
available <- !is_cran_check && requireNamespace("mice", quietly = TRUE) &&
  requireNamespace("cluster", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = !is_cran_check && available)

## ----attributes---------------------------------------------------------------
# library(mfrmr)
# set.seed(20260921)
# n <- 120L
# raters <- data.frame(
#   Rater = sprintf("R%03d", seq_len(n)),
#   ExperienceYears = sample(0:30, n, replace = TRUE),
#   AnnualRatings = round(exp(rnorm(n, log(250), 0.6))),
#   Specialty = factor(sample(c("Language", "Science", "Social science"),
#                             n, replace = TRUE)),
#   TrainingLevel = ordered(sample(c("Introductory", "Intermediate", "Advanced"),
#                                  n, replace = TRUE),
#                           levels = c("Introductory", "Intermediate", "Advanced")),
#   Certified = sample(c(FALSE, TRUE), n, replace = TRUE),
#   WorkshopHours = sample(0:80, n, replace = TRUE)
# )
# has_mentoring_role <- sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(.7, .3))
# raters$MentoringYears <- ifelse(has_mentoring_role,
#   pmin(raters$ExperienceYears, sample(0:10, n, replace = TRUE)), NA_real_)
# 
# # Introduce unrecorded values only after generating the fictional attributes.
# feature_names <- c("ExperienceYears", "AnnualRatings", "Specialty",
#                    "TrainingLevel", "Certified")
# missing_counts <- c(12L, 18L, 6L, 10L, 10L)
# for (j in seq_along(feature_names)) {
#   raters[[feature_names[j]]][sample.int(n, missing_counts[j])] <- NA
# }
# head(raters)

## ----missing-reasons----------------------------------------------------------
# all_features <- c(feature_names, "MentoringYears")
# initial <- mfrm_features(raters, "Rater", all_features)
# reasons <- initial$missing
# reasons$Reason <- ifelse(reasons$Feature == "MentoringYears",
#                         "Not applicable: no mentoring role", "Not recorded")
# review <- mfrm_features(raters, "Rater", all_features, reasons)
# summary(review)
# table(review$missing$Reason)

## ----selected-features--------------------------------------------------------
# selected_reasons <- subset(review$missing, Feature %in% feature_names)
# features <- mfrm_features(raters, "Rater", feature_names, selected_reasons)
# features$feature_summary
# table(features$row_summary$Complete)

## ----imputation, fig.width=8, fig.height=6, fig.alt="Imputation-chain means and standard deviations across iterations for experience years and annual rating workload. These traces help review mixing; they do not establish that the imputation model is appropriate."----
# eligible <- subset(features$missing, Reason == "Not recorded")
# imputation_data <- raters[c("Rater", feature_names, "WorkshopHours")]
# where <- matrix(FALSE, nrow(imputation_data), ncol(imputation_data),
#                 dimnames = dimnames(imputation_data))
# where[cbind(match(eligible$ID, imputation_data$Rater),
#             match(eligible$Feature, names(imputation_data)))] <- TRUE
# methods <- mice::make.method(imputation_data, where = where)
# methods[feature_names] <- c("pmm", "pmm", "polyreg", "polr", "logreg")
# methods["Rater"] <- ""
# predictors <- mice::make.predictorMatrix(imputation_data)
# predictors[, "Rater"] <- 0
# predictors["Rater", ] <- 0
# model <- mice::mice(imputation_data, m = 5, maxit = 5,
#   method = methods, predictorMatrix = predictors, where = where,
#   seed = 724, printFlag = FALSE)
# model$loggedEvents
# plot(model, c("ExperienceYears", "AnnualRatings"), layout = c(2, 2))

## ----grouping-----------------------------------------------------------------
# experience_weights <- setNames(rep(1, length(feature_names)), feature_names)
# experience_weights["ExperienceYears"] <- 3
# analyses <- list(
#   TwoGroups = mfrm_cluster_imputed(features, model, eligible, k = 2),
#   ThreeGroups = mfrm_cluster_imputed(features, model, eligible, k = 3),
#   ExperienceWeighted = mfrm_cluster_imputed(features, model, eligible, k = 3,
#                                           weights = experience_weights)
# )
# comparison <- mfrm_cluster_compare(analyses)
# comparison$analysis_summary
# summary(comparison)
# weight_effect <- subset(comparison$comparison_summary,
#   First == "ThreeGroups" & Second == "ExperienceWeighted")

## ----feature-selection--------------------------------------------------------
# background_names <- setdiff(feature_names, "AnnualRatings")
# background_features <- mfrm_features(raters, "Rater", background_names,
#   subset(selected_reasons, Feature %in% background_names))
# background_only <- mfrm_cluster_imputed(background_features, model,
#   subset(eligible, Feature %in% background_names), k = 3)
# feature_comparison <- mfrm_cluster_compare(list(
#   WithWorkload = analyses$ThreeGroups, WithoutWorkload = background_only))
# feature_comparison$analysis_summary
# feature_comparison$weights
# summary(feature_comparison)

## ----interpretation-----------------------------------------------------------
# one_completion <- analyses$ThreeGroups$analyses[[1]]
# one_completion$profiles$numeric
# subset(one_completion$profiles$categorical, Feature == "TrainingLevel")
# head(one_completion$membership)
# analyses$ThreeGroups$co_membership[1:5, 1:5]

## ----silhouette, fig.width=7, fig.height=5, fig.alt="Rater silhouette widths grouped by the descriptive partition in one completion. The dashed line is the overall mean; negative widths indicate profiles to inspect, not poor rater quality."----
# plot(one_completion)

## ----profiles, fig.width=7, fig.height=5, fig.alt=c("Experience in years summarized by descriptive rater group, with means and medians in the original units and no confidence intervals.", "Proportions of ordered training levels within each descriptive rater group. These profiles describe one imputed table, not effects of training.")----
# plot(one_completion, type = "profile", feature = "ExperienceYears")
# plot(one_completion, type = "profile", feature = "TrainingLevel")

## ----profile-ggplot, fig.width=7, fig.height=5, fig.alt="Experience in years summarized by the saved groups from one completed feature table. Circles mark means and triangles mark medians, offset vertically so equal values remain visible. Group labels give counts; no confidence intervals are shown."----
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   profile_figure <- as_ggplot(one_completion, type = "profile",
#     feature = "ExperienceYears", preset = "monochrome") +
#     ggplot2::labs(title = NULL)
#   print(profile_figure)
# }

## ----co-membership, fig.width=7, fig.height=6, fig.alt="Pairwise co-membership proportions across imputations for the first twenty rater IDs. The scale runs from zero to one; unavailable pairs are grey. These proportions are imputation sensitivity, not sampling stability."----
# plot(analyses$ThreeGroups, ids = raters$Rater[1:20])

## ----plot-values--------------------------------------------------------------
# view <- plot(analyses$ThreeGroups, ids = raters$Rater[1:5], draw = FALSE)
# plot_data(view)$matrix

## ----co-membership-ggplot, fig.width=7, fig.height=6, fig.alt="Co-membership fractions for the five selected rater IDs, in the saved order on both axes. The scale is fixed from zero to one; crosses mark unavailable cells, if any. Values use all supplied imputations and do not measure sampling stability."----
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   co_membership_figure <- as_ggplot(view, component = "matrix") +
#     ggplot2::labs(title = NULL)
#   print(co_membership_figure)
# }

## ----profile-observations-----------------------------------------------------
# values <- plot_data(one_completion, "observations", type = "profile",
#   feature = "ExperienceYears")
# head(values)

## ----raincloud, fig.width=7, fig.height=4.5-----------------------------------
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   displayed <- subset(values, Included & is.finite(Value))
#   displayed$Group <- factor(displayed$Cluster)
#   displayed$Position <- as.integer(displayed$Group)
#   clouds <- lapply(split(displayed, displayed$Group), function(d) {
#     if (nrow(d) < 5L || length(unique(d$Value)) < 2L) return(NULL)
#     density <- stats::density(d$Value, from = min(d$Value), to = max(d$Value))
#     center <- d$Position[1]
#     data.frame(X = c(center, center + .38*density$y/max(density$y), center),
#       Value = c(density$x[1], density$x, tail(density$x, 1)),
#       Group = d$Group[1])
#   })
#   clouds <- do.call(rbind, clouds)
#   raincloud <- ggplot2::ggplot(displayed,
#     ggplot2::aes(x = Position, y = Value, group = Group))
#   if (!is.null(clouds)) raincloud <- raincloud +
#     ggplot2::geom_polygon(data = clouds,
#       ggplot2::aes(x = X, y = Value, group = Group, fill = Group),
#       inherit.aes = FALSE, alpha = .35, colour = NA)
#   raincloud <- raincloud +
#     ggplot2::geom_boxplot(width = .08, outlier.shape = NA) +
#     ggplot2::geom_point(ggplot2::aes(x = Position - .15), alpha = .5, size = 1.2,
#       position = ggplot2::position_jitter(width = .045, height = 0, seed = 1)) +
#     ggplot2::scale_x_continuous(breaks = seq_along(levels(displayed$Group)),
#       labels = paste("Group", levels(displayed$Group))) +
#     ggplot2::labs(x = NULL, y = "Rating experience (years)") +
#     ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
#   print(raincloud)
# }

## ----save-profile-values, eval=FALSE------------------------------------------
# write.csv(values, "rater-experience-long.csv", row.names = FALSE)
# saveRDS(one_completion, "feature-groups.rds")
# ggplot2::ggsave("experience-raincloud.pdf", raincloud, width = 7, height = 4.5)

## ----hierarchical-comparison--------------------------------------------------
# average <- mfrm_cluster_imputed(features, model, eligible, k = 3,
#   method = "hierarchical", linkage = "average")
# complete <- mfrm_cluster_imputed(features, model, eligible, k = 3,
#   method = "hierarchical", linkage = "complete")
# method_comparison <- mfrm_cluster_compare(list(
#   PAM = analyses$ThreeGroups, Average = average, Complete = complete))
# method_comparison$analysis_summary
# summary(method_comparison)

## ----dendrogram, fig.width=7, fig.height=5, fig.alt="Average-linkage hierarchy from Gower dissimilarities in one completed rater-feature table. Boxes mark three groups; all 120 raters are included. Heights are dissimilarities, not significance or branch support."----
# plot(average$analyses[[1]])

## ----other-facet-features-----------------------------------------------------
# set.seed(7241)
# persons <- data.frame(
#   Person = sprintf("P%03d", 1:240),
#   LanguageStudyYears = sample(0:15, 240, replace = TRUE),
#   PreparationHours = round(exp(rnorm(240, log(20), 0.5))),
#   InstructionMode = factor(sample(c("Classroom", "Online", "Hybrid"),
#                                   240, replace = TRUE)))
# tasks <- data.frame(
#   Task = sprintf("T%02d", 1:12),
#   TimeLimitMinutes = sample(c(10, 20, 30), 12, replace = TRUE),
#   WordLimit = sample(c(100, 200, 400), 12, replace = TRUE),
#   ResponseFormat = factor(rep(c("Summary", "Argument", "Explanation"), 4)))
# person_groups <- mfrm_cluster_pam(mfrm_features(persons, "Person",
#   c("LanguageStudyYears", "PreparationHours", "InstructionMode")), k = 3)
# task_groups <- mfrm_cluster_pam(mfrm_features(tasks, "Task",
#   c("TimeLimitMinutes", "WordLimit", "ResponseFormat")), k = 3)

## ----sparse-rating-roster-----------------------------------------------------
# assignments <- do.call(rbind, lapply(persons$Person, function(id) {
#   do.call(rbind, lapply(sample(tasks$Task, 3), function(task) {
#     data.frame(Person = id, Rater = sample(raters$Rater, 2), Task = task)
#   }))
# }))
# ratings <- assignments
# ratings$Score <- sample(1:4, nrow(ratings), replace = TRUE)
# ratings$Score[sample(nrow(ratings), 72)] <- NA_integer_
# rating_review <- describe_mfrm_data(ratings, person = "Person",
#   facets = c("Rater", "Task"), score = "Score",
#   rating_min = 1, rating_max = 4, keep_original = TRUE,
#   include_agreement = FALSE, expected_design = assignments)
# rating_review$structural_missingness$summary
# rating_review$design_connectivity

## ----link-facet-groups--------------------------------------------------------
# rating_key <- c("Person", "Rater", "Task")
# stopifnot(!anyDuplicated(assignments[rating_key]),
#           !anyDuplicated(ratings[rating_key]))
# labeled <- merge(assignments, ratings, by = rating_key,
#                  all.x = TRUE, sort = FALSE)
# labeled$PersonGroup <- person_groups$membership$Cluster[
#   match(labeled$Person, person_groups$membership$ID)]
# labeled$RaterGroup <- one_completion$membership$Cluster[
#   match(labeled$Rater, one_completion$membership$ID)]
# labeled$TaskGroup <- task_groups$membership$Cluster[
#   match(labeled$Task, task_groups$membership$ID)]
# group_columns <- c("PersonGroup", "RaterGroup", "TaskGroup")
# stopifnot(!anyNA(labeled[group_columns]))
# for (column in group_columns) labeled[[column]] <- factor(labeled[[column]], 1:3)
# labeled$Assigned <- 1L
# labeled$Observed <- as.integer(!is.na(labeled$Score))
# planned_counts <- xtabs(Assigned ~ PersonGroup + RaterGroup + TaskGroup,
#                         data = labeled)
# observed_counts <- xtabs(Observed ~ PersonGroup + RaterGroup + TaskGroup,
#                          data = labeled)
# group_coverage <- as.data.frame(planned_counts, responseName = "Assigned")
# group_coverage$Observed <- as.vector(observed_counts)
# group_coverage$Coverage <- ifelse(group_coverage$Assigned > 0,
#   group_coverage$Observed / group_coverage$Assigned, NA_real_)
# group_coverage

## ----numeric-grouping---------------------------------------------------------
# numeric_names <- c("ExperienceYears", "AnnualRatings", "WorkshopHours")
# numeric_review <- mfrm_features(raters, "Rater", numeric_names)
# numeric_direct <- mfrm_cluster_imputed(numeric_review, model,
#   numeric_review$missing, k = 3, method = "kmeans", seed = 42)
# numeric_reduced <- mfrm_cluster_imputed(numeric_review, model,
#   numeric_review$missing, k = 3, method = "kmeans", components = 2, seed = 42)
# numeric_comparison <- mfrm_cluster_compare(list(
#   StandardizedFeatures = numeric_direct, TwoComponents = numeric_reduced))
# numeric_comparison$analysis_summary
# summary(numeric_comparison)

## ----numeric-pca-views, fig.width=8, fig.height=6, fig.alt=c("Explained-variance scree plot for numeric rater features in one completion. Filled points identify retained components and open points identify omitted components.", "Scores on the first two principal components, with descriptive k-means groups distinguished by colour and point shape. This projection is not a latent ability or rater-quality scale.", "First-component loadings for the transformed numeric rater features. Their signs are arbitrary and the coefficients are not original-unit correlations.", "Experience in years summarized by the k-means groups fitted in the retained principal-component space. The original-unit profile describes this completion only.")----
# first_numeric <- numeric_reduced$analyses[[1]]
# summary(first_numeric$pca)
# first_numeric$pca$loadings
# plot(first_numeric$pca)  # The full explained-variance table, including omitted PCs.
# plot(first_numeric$pca, type = "scores", groups = first_numeric, labels = FALSE)
# plot(first_numeric$pca, type = "loadings", components = 1)
# plot(first_numeric, type = "profile", feature = "ExperienceYears")

## ----numeric-pca-ggplot, fig.width=8, fig.height=6, fig.alt="Scores on the first two external-feature principal components for one completed feature dataset. K-means groups use grayscale colours and different point shapes. Headings and entity labels are omitted; the saved data still identify the completion's transformation and excluded entities."----
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   numeric_figure <- as_ggplot(first_numeric$pca, type = "scores",
#     groups = first_numeric, labels = FALSE, preset = "monochrome") +
#     ggplot2::labs(title = NULL, subtitle = NULL)
#   print(numeric_figure)
#   plot_data(numeric_figure)$components
# }

## ----feature-save-reopen------------------------------------------------------
# feature_file <- tempfile(fileext = ".rds")
# feature_view <- plot(first_numeric$pca, type = "scores", groups = first_numeric,
#   labels = FALSE, preset = "monochrome", draw = FALSE)
# saveRDS(list(comparison = numeric_comparison, view = feature_view), feature_file)
# saved_features <- readRDS(feature_file)
# summary(saved_features$comparison)
# saved_features$comparison$analysis_summary  # Each setting and completion.
# plot_data(saved_features$view)    # Inspect this first-completion view.
# feature_summary_file <- tempfile(fileext = ".csv")
# write.csv(summary(saved_features$comparison), feature_summary_file, row.names = FALSE)
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   saved_figure <- as_ggplot(saved_features$view)
# }

