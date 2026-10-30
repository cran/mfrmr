test_that("Q3 counts unordered pairs and retains unavailable comparisons", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  levels <- sort(unique(as.character(diag$obs$Rater)))
  person <- as.integer(factor(diag$obs$Person))
  diag$obs$StdResidual <- ifelse(diag$obs$Rater == levels[1], person,
    ifelse(diag$obs$Rater == levels[2], -person,
      ifelse(diag$obs$Rater == levels[3], 0, NA_real_)))
  q <- q3_statistic(fit, diag)
  expect_match(q$notes, "Standardized residuals")
  expect_match(q$notes, "not raw-residual Yen Q3")
  expect_equal(nrow(q$pairs), choose(length(levels), 2))
  expect_equal(q$summary$AvailablePairs, 1)
  expect_equal(q$summary$UnavailablePairs, nrow(q$pairs) - 1)
  expect_equal(q$summary$YenFlagged, 1)
  expect_equal(q$pairs$Q3[is.finite(q$pairs$Q3)], -1)
  expect_true(all(is.na(q$pairs$YenFlag[is.na(q$pairs$Q3)])))
  expect_true(any(grepl("No residual variation", q$pairs$Interpretation)))
  expect_true(any(grepl("Insufficient shared persons", q$pairs$Interpretation)))
  heat <- plot_local_dependence_heatmap(fit, diag, draw = FALSE)
  expect_identical(q$notes, heat$data$caption)
  expect_equal(heat$data$matrix, t(heat$data$matrix))
  expect_true(is.na(heat$data$matrix[levels[4], levels[4]]))

  diag$obs$StdResidual[] <- NA_real_
  empty <- q3_statistic(fit, diag)
  expect_equal(empty$summary$AvailablePairs, 0)
  expect_true(is.na(empty$summary$YenFlagged))
  expect_true(is.na(empty$summary$MaxAbsQ3))
  expect_output(print(empty), "0 available")
  q$calculation_version <- NULL
  expect_error(print(q), "Recreate this Q3")

  diag$obs$Score[1] <- diag$obs$Score[1] + 1
  expect_error(q3_statistic(fit, diag), "fit/diagnostics mismatch")
  expect_error(plot_local_dependence_heatmap(fit, min_pairs = 3.5), "integer")
})

test_that("PCA preserves valid correlations and refuses fabricated structure", {
  mat <- cbind(A = 1:6, B = c(1, 3, 2, 6, 4, 5), C = c(3, 6, 1, 2, 5, 4))
  obs <- data.frame(Person = rep(seq_len(nrow(mat)), times = ncol(mat)),
                    Item = rep(colnames(mat), each = nrow(mat)),
                    StdResidual = as.vector(mat))
  diagnostics <- list(obs = obs, facet_names = "Item")
  p <- analyze_residual_pca(diagnostics, mode = "both", pca_max_factors = 1)
  expect_equal(p$overall$cor_matrix, cor(mat))
  expect_equal(p$overall_table$Eigenvalue, eigen(cor(mat))$values)
  expect_equal(p$by_facet$Item$cor_matrix, cor(mat))
  expect_equal(ncol(p$overall$pca$loadings), 1L)

  diagnostics$residual_pca_overall <- p$overall
  diagnostics$residual_pca_by_facet <- p$by_facet
  expanded <- analyze_residual_pca(diagnostics, pca_max_factors = 2)
  expect_equal(ncol(expanded$overall$pca$loadings), 2L)
  diagnostics$residual_pca_overall$calculation_version <- NULL
  diagnostics$residual_pca_overall$pca$values[] <- 99
  current <- analyze_residual_pca(diagnostics, mode = "overall", pca_max_factors = 1)
  expect_equal(current$overall_table$Eigenvalue, eigen(cor(mat))$values)
  p$calculation_version <- NULL
  expect_error(plot_residual_pca(p, draw = FALSE), "Recreate this residual PCA")
  expect_error(summary(p), "Recreate this residual PCA")

  diagnostics$obs$StdResidual[diagnostics$obs$Item == "C"] <- 0
  invalid <- analyze_residual_pca(diagnostics, mode = "both", parallel = TRUE, parallel_reps = 3)
  expect_equal(nrow(invalid$overall_table), 0L)
  expect_equal(nrow(invalid$by_facet_table), 0L)
  expect_match(invalid$errors$overall, "residual variation")
  expect_output(print(summary(invalid)), "residual variation")
  expect_equal(nrow(summary(invalid)$preview), 0L)
  expect_false(any(invalid$parallel_status$ParallelAvailable))
  expect_error(plot_residual_pca(invalid, draw = FALSE), "No eigenvalues")

  disconnected <- cbind(A = c(-1, 1, NA, NA), B = c(NA, NA, -1, 1))
  expect_match(residual_pca_correlation(disconnected)$error, "shared-person")
  indefinite <- cbind(A = c(-1, 1, -1, 1, NA, NA),
                     B = c(-1, 1, NA, NA, -1, 1),
                     C = c(NA, NA, -1, 1, 1, -1))
  bad <- residual_pca_correlation(indefinite)
  expect_match(bad$error, "positive semidefinite")
  expect_equal(min(eigen(bad$cor_matrix)$values), -1)
  expect_null(residual_pca_correlation(cbind(A = 1:4, B = 4:1))$error)
  ambiguous <- data.frame(Person = 1:3, A = c("a_b", "a", "c"),
                          B = c("c", "b_c", "d"), StdResidual = 1:3)
  expect_match(compute_pca_overall(ambiguous, c("A", "B"))$error, "ambiguous")
})

test_that("permutation comparison withholds a selected successful subset", {
  calls <- 0L
  local_mocked_bindings(residual_parallel_eigenvalues = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) c(1.5, 0.5) else numeric(0)
  })
  pa <- compute_residual_parallel_analysis(cbind(A = 1:4, B = c(1, 3, 2, 4)),
    observed_eigenvalues = c(1.5, 0.5), reps = 2, seed = 10)
  expect_equal(pa$successful_reps, 1L)
  expect_equal(nrow(pa$table), 0L)
  expect_match(pa$error, "1 of 2")
  expect_match(pa$error, "reference distribution")
})

test_that("PCA displays avoid automatic dimensionality decisions and internal codes", {
  obs <- expand.grid(Person = 1:8, Item = c("A", "B", "C"))
  obs$StdResidual <- c(1:8, c(1, 3, 5, 7, 2, 4, 6, 8), 8:1)
  p <- analyze_residual_pca(list(obs = obs, facet_names = "Item"),
                            parallel = TRUE, parallel_reps = 3, seed = 17)
  text <- paste(capture.output(print(summary(p))), collapse = "\n")
  expect_match(text, "fitted-model uncertainty")
  expect_false(grepl("InferenceTier|ReportingUse|screening_only|residual_permutation|mfrm_residual_pca", text))
  plot <- plot_residual_pca(p, draw = FALSE)
  expect_equal(nrow(plot$data$reference_lines), 1L)
  expect_false(any(grepl("second dim|Critical minimum", plot$data$legend$label)))
  old <- summary(p)
  old$calculation_version <- NULL
  expect_error(print(old), "Recreate this residual PCA")
})

test_that("reports retain unavailable PCA scopes without recomputing old results", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  level <- as.character(diag$obs$Rater[1])
  diag$obs$StdResidual[diag$obs$Rater == level] <- 0
  p <- analyze_residual_pca(diag, mode = "both")
  diag$residual_pca_overall <- p$overall
  diag$residual_pca_by_facet <- p$by_facet
  diag$residual_pca_mode <- "both"
  apa <- build_apa_outputs(fit, diag)
  expect_false(apa$contract$availability$has_pca_overall)
  text <- gsub("[[:space:]]+", " ", apa$report_text)
  expect_match(text, "Facet-specific residual PCA was unavailable: Rater")
  diag$residual_pca_overall$calculation_version <- NULL
  local_mocked_bindings(compute_pca_overall = function(...) stop("Unexpected PCA"),
                        compute_pca_by_facet = function(...) stop("Unexpected PCA"))
  result <- safe_residual_pca(diag)
  expect_match(result$errors$overall, "Recreate the stored residual PCA")
})
