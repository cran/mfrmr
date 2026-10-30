# Retained example_core fits and an independent marginal-information calculation.
# No optimization or data generation is needed to check these outputs.
local({
  fixtures <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))
  for (case in names(fixtures)) test_that(paste(case, "anchored intervals preserve joint uncertainty"), {
    x <- fixtures[[case]]; fit <- x$fit; J <- x$rater_jacobian
    info <- compute_mml_parameter_covariance(fit)
    scores <- mfrm_person_likelihood_scores(fit)
    expect_identical(info$status, "ok")
    expect_lt(max(abs(info$cov-x$covariance)), 1e-6)
    expect_lt(max(abs(unname(scores)-x$scores)), 1e-8)
    local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
      compute_mml_parameter_covariance = function(res, ...) {
        expect_identical(res, fit)
        info
      }, mfrm_person_likelihood_scores = function(res) {
        expect_identical(res, fit)
        scores
      }, .package = "mfrmr")
    S <- x$covariance %*% crossprod(x$scores) %*% x$covariance
    model <- mfrm_facet_intervals(fit, "Rater")
    sandwich <- mfrm_facet_intervals(fit, "Rater", method = "sandwich")
    expect_equal(unname(model$covariance), J %*% x$covariance %*% t(J), tolerance = 1e-6)
    expect_equal(unname(sandwich$covariance), J %*% S %*% t(J), tolerance = 1e-6)
    expect_lt(max(abs(rowSums(model$covariance))), 1e-12)
    fixed <- if (case == "RSM_group") c(TRUE, FALSE, FALSE, TRUE) else c(TRUE, FALSE, FALSE, FALSE)
    expect_identical(model$table$Status == "fixed", fixed)
    expect_true(all(model$table$SE[fixed] == 0 & is.na(model$table$Lower[fixed])))
    expect_true(all(model$table$SE[!fixed] > 0 & is.finite(model$table$Lower[!fixed])))
    C <- x$contrasts
    ci <- mfrm_facet_intervals(fit, "Rater", C[, 4:1], method = "sandwich", level = .9)
    expected <- C %*% J %*% S %*% t(J) %*% t(C)
    expect_equal(ci$covariance, expected, tolerance = 1e-6)
    expect_identical(ci$table$Status[3], "fixed")
    expect_identical(ci$table$SE[3], 0)
    expect_true(is.na(ci$table$Lower[3]))
    expect_identical(ci$table$Status[1:2], rep("available", 2))
    expect_equal(ci$table$Lower[1:2], ci$table$Estimate[1:2] - qnorm(.95)*ci$table$SE[1:2])
    if (case == "RSM_group") {
      expect_identical(ci$table$Status[4], "fixed")
      expect_equal(ci$table$Estimate[4], .1, tolerance = 1e-12)
    }
    clusters <- data.frame(Person = rev(rownames(scores)), Cluster = rep(1:24, each = 2))
    clustered <- mfrm_facet_intervals(fit, "Rater", C, method = "sandwich",
      clusters = clusters, adjust = TRUE)
    U <- rowsum(x$scores, clusters$Cluster[match(rownames(scores), clusters$Person)], reorder = FALSE)
    V <- x$covariance %*% crossprod(U) %*% x$covariance * 24/23
    expect_equal(clustered$covariance, C %*% J %*% V %*% t(J) %*% t(C), tolerance = 1e-6)
    expect_equal(clustered$settings$adjustment_factor, 24/23)
    clusters$Cluster <- rep(1:2, each = 24)
    refused <- mfrm_facet_intervals(fit, "Rater", C, method = "sandwich", clusters = clusters)
    expect_true(all(refused$table$Status[1:2] == "insufficient_cluster_rank"))
    expect_true(all(is.na(refused$table$Lower[1:2])))
    expect_identical(refused$table$Status[3], "fixed")
    path <- tempfile(fileext = ".rds")
    saveRDS(ci, path); reopened <- readRDS(path); unlink(path)
    expect_identical(reopened, ci)
    expect_identical(plot_data(plot(reopened, draw = FALSE))$table, ci$table)
    expect_identical(apa_table(reopened)$table$Status, ci$table$Status)
    results <- mfrm_results(fit, intervals = list(contrast = reopened), include = "fit", compute = "never")
    expect_identical(results$tables$facet_contrast_intervals$Status, ci$table$Status)
  })
})
