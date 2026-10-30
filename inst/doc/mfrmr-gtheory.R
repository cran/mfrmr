## ----setup, include=FALSE-----------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
  fig.height = 4.5, dev.args = list(bg = "white"), eval = !is_cran_check)


## ----crossed-data-------------------------------------------------------------
library(mfrmr)
data("mfrmr_example_gtheory_crossed", package = "mfrmr")
ratings <- mfrmr_example_gtheory_crossed
head(ratings)


## ----crossed-gstudy-----------------------------------------------------------
g <- mfrm_g_study(ratings, scores = c("Content", "Organization"),
  person = "Person", rater = "Rater", task = "Task")
print(g, digits = 3)


## ----crossed-score-components, fig.height=5, fig.alt="Two panels show the seven variance components for Content and Organization on shared axes. They include Person, Rater, Task, three two-way interactions and the highest interaction combined with residual."----
components <- plot(g, score = c("Content", "Organization"), ncol = 2,
  preset = "standard")


## ----crossed-plans------------------------------------------------------------
plans <- data.frame(Raters = c(2, 3, 4), Tasks = c(6, 4, 3))
d <- mfrm_d_study(g, design_grid = plans,
  weights = c(Content = .5, Organization = .5))
s <- summary(d, intervals = TRUE, parm = c("G", "Phi"), assumption = "normal")
print(s, digits = 3)


## ----crossed-error-components, fig.height=7, fig.alt="Relative and absolute error-variance contributions for three twelve-rating plans. Person-by-rater error falls with more raters; Person-by-task error rises when fewer tasks are used. The highest interaction and residual contribution is unchanged because the product of counts remains twelve."----
errors <- plot(d, type = "error_components")
plot_data(errors, "table") # Divisors, variance and covariance contributions.
plot_data(errors, "totals") # Original totals, SEM, G/Phi and availability.


## ----crossed-panels, fig.height=4.5, fig.alt="Two panels show absolute dependability Phi for the same composite as task count increases. Each panel fixes two or four raters. Shaded regions are approximate 95 percent MLS intervals, not simultaneous confidence bands."----
d_grid <- mfrm_d_study(g,
  expand.grid(Raters = c(2, 4), Tasks = c(3, 6, 12)),
  weights = c(Content = .5, Organization = .5))
p <- plot(d_grid, metric = "Phi", x_var = "Tasks", panel_by = "Raters",
  ncol = 2, intervals = TRUE, assumption = "normal", preset = "standard")


## ----finer-planning-grid, eval=FALSE------------------------------------------
# d_fine <- mfrm_d_study(g,
#   expand.grid(Raters = c(2, 4), Tasks = 1:12),
#   weights = c(Content = .5, Organization = .5))
# plot(d_fine, metric = "Phi", x_var = "Tasks", panel_by = "Raters",
#   ncol = 2, preset = "standard")


## ----crossed-comparison-------------------------------------------------------
comparison <- mfrm_multivariate_d_compare(d, reference = 1,
  assumption = "normal")
print(comparison, digits = 3)


## ----save-and-reopen----------------------------------------------------------
analysis <- list(study = d, summary = s, comparison = comparison, plot = p,
  decisions = list(score = "Equal-weight Content/Organization mean",
    comparison = "Twelve ratings per person; task burden and costs differ",
    assumption = "Independent normal random effects"))
result_file <- tempfile(fileext = ".rds")
saveRDS(analysis, result_file)
reopened <- readRDS(result_file)
identical(reopened$summary, s)  # TRUE: the full summary was preserved.


## ----gtheory-report-route, eval=FALSE-----------------------------------------
# report <- mfrm_report(d, intervals = s$intervals, style = "apa")
# bundle <- build_summary_table_bundle(report)
# apa_table(report, which = "coefficients")
# export_summary_appendix(bundle, output_dir = "planning-report")
# saveRDS(report, "planning-report/report.rds")
# 
# # A complete editable manuscript starting point, using bundled data:
# template <- system.file("templates", "gtheory-report.Rmd", package = "mfrmr")
# file.copy(template, "gtheory-report.Rmd")
# rmarkdown::render("gtheory-report.Rmd")


## ----criterion-relations------------------------------------------------------
g$components$Person


## ----criterion-covariance, fig.height=4.5, fig.alt="One panel shows Content-Organization covariance estimates and conservative pointwise 95 percent intervals for all seven sources. Every interval includes zero; the Task and Rater intervals are much wider than the others. These are covariance units, not correlations or halo scores."----
relations <- plot(g, type = "covariance", score = c("Content", "Organization"),
  intervals = TRUE, method = "regression_region", assumption = "normal")


## ----criterion-covariance-table-----------------------------------------------
relation_intervals <- plot_data(relations, "intervals")[[1]]
relation_intervals[, c("Source", "Estimate", "Lower", "Upper", "Status")]


## ----save-criterion-evidence--------------------------------------------------
analysis$criterion_relations <- relations
saveRDS(analysis, result_file)
reopened <- readRDS(result_file)
identical(plot_data(reopened$criterion_relations, "intervals")[[1]],
  relation_intervals)  # TRUE: intervals and their qualifications were retained.


## ----task-data----------------------------------------------------------------
tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
g_tasks <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
d_tasks <- mfrm_d_study(g_tasks, data.frame(Tasks = c(6, 12)),
  weights = c(V = .5, W = .5))
change <- summary(d_tasks, intervals = TRUE, parm = "Phi", reference = 1,
  assumption = "normal")
print(change, digits = 3)


## ----nested-intervals---------------------------------------------------------
data("mfrmr_example_gtheory_nested", package = "mfrmr")
g_nested <- mfrm_g_study(mfrmr_example_gtheory_nested,
  c("Content", "Organization"), nesting = c(Rater = "Task"))
d_nested <- mfrm_d_study(g_nested,
  data.frame(Raters = 3, Tasks = c(6, 12)),
  weights = c(Content = .5, Organization = .5))
s_nested <- summary(d_nested, intervals = TRUE, parm = "G",
  assumption = "normal")
print(s_nested, digits = 3)


## ----crossed-interval-methods, fig.alt="One forest plot groups three future plans vertically. Adjacent rows compare selected interval methods for Phi on the same axis; color, shape and line type distinguish the methods."----
methods <- c("mls", "delta")
if (requireNamespace("lpSolve", quietly = TRUE))
  methods <- c(methods, "mean_square_region")
ci_methods <- confint(d, parm = "Phi", method = methods, assumption = "normal")
plot(ci_methods)
ci_methods[c("Scenario", "Method", "Estimate", "SE", "Lower", "Upper")]


## ----crossed-small-source-intervals, fig.alt="The same Phi method comparison now uses a small complete source. Wide or unavailable intervals remain visible on the shared axis; the plot does not measure coverage."----
small <- subset(ratings, Person %in% head(unique(Person), 10) &
  Rater %in% head(unique(Rater), 2) & Task %in% head(unique(Task), 4))
g_small <- mfrm_g_study(small, c("Content", "Organization"))
d_small <- mfrm_d_study(g_small, plans,
  weights = c(Content = .5, Organization = .5))
ci_small <- confint(d_small, parm = "Phi", method = methods, assumption = "normal")
plot(ci_small, preset = "monochrome")


## ----incomplete-residual------------------------------------------------------
observed <- ratings[-c(1, 8, 23), ]
sparse_g <- mfrm_g_study(observed,
  scores = c("Content", "Organization"), method = "minque0")
residual_ci <- confint(sparse_g, parm = "Residual", score = "Content",
  method = "residual_contrast", assumption = "normal")
residual_ci[, c("Estimate", "IntervalCenter", "Lower", "Upper", "Status")]
attr(residual_ci, "residual_projection")$df


## ----incomplete-relative-sem--------------------------------------------------
observed_tasks <- tasks[(tasks$Person + tasks$Task) %% 3 != 0, ]
task_g <- mfrm_g_study(observed_tasks, "V", rater = NULL,
  method = "minque0")
task_d <- mfrm_d_study(task_g, data.frame(Tasks = c(6, 12)))
sem_summary <- summary(task_d, intervals = TRUE, parm = "RelativeSEM",
  method = "residual_contrast", assumption = "normal")
print(sem_summary, digits = 3)
sem_change <- confint(task_d, parm = "RelativeSEM", reference = 1,
  method = "residual_contrast", assumption = "normal")
sem_change[, c("Tasks", "Estimate", "IntervalCenter", "Lower", "Upper")]


## ----incomplete-relative-g, fig.alt="Relative G increases with task count. The shaded interval uses a normal-model spectral F pivot from incomplete source ratings. Its wide lower range shows source uncertainty; it is not the interval for the difference between plans."----
g_change <- confint(task_d, parm = "G", reference = 1,
  method = "spectral_f", assumption = "normal")
g_change[, c("Tasks", "Estimate", "Lower", "Upper", "Status")]
plot(task_d, metric = "G", intervals = TRUE, method = "spectral_f",
  assumption = "normal")

