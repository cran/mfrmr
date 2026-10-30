# Consumer-contract tests: inject a verified review at the covariance boundary.
# This is not a claim that the example ratings have weak observed information.
local({
  ratings <- load_mfrmr_data("example_core")
  fits <- lapply(c("RSM", "PCM"), function(model) fit_mfrm(ratings,
    "Person", c("Rater", "Criterion"), "Score", model = model,
    step_facet = if (model == "PCM") "Criterion" else NULL))
  original <- compute_mml_parameter_covariance
  h <- diag(c(1e-10, 2))
  review <- mfrm_refine_mml_information(c(0, 0),
    function(p) sum(p * (h %*% p))/2, function(p) drop(h %*% p), h)$review
  reviewed <- function(fit) {
    info <- original(fit)
    info$solution_information$inverse_review <- review
    info$detail <- review$Detail
    info
  }

  test_that("fixed-facet cautions survive summaries, graphics and saved reports", {
    ordinary <- lapply(fits, mfrm_facet_intervals, facet = "Rater")
    local_mocked_bindings(compute_mml_parameter_covariance = reviewed)
    for (i in seq_along(fits)) {
      f <- fits[[i]]
      expect_warning(ci <- mfrm_facet_intervals(f, "Rater"), "Ill-conditioned")
      expect_equal(ci$table[names(ordinary[[i]]$table)], ordinary[[i]]$table)
      expect_equal(ci$parameter_covariance, ordinary[[i]]$parameter_covariance)
      expect_identical(ci$information_review, review)
      expect_identical(unique(summary(ci)$InferenceCaution), review$Detail)
      expect_output(print(ci), "Caution: Ill-conditioned")
      expect_warning(robust <- mfrm_facet_intervals(f, "Rater", method = "sandwich"),
        "Ill-conditioned")
      expect_identical(robust$cautions, review$Detail)
      res <- mfrm_results(f, intervals = ci, include = c("fit", "plots"), compute = "never")
      expect_identical(res$tables$facet_inference_information_review, review)
      expect_true(review$Detail %in% res$notes)
      expect_identical(apa_table(ci)$table$InferenceCaution, ci$table$InferenceCaution)
      report <- mfrm_report(res)
      expect_match(report$markdown, "Ill-conditioned")
      plot <- plot_data(ci)
      expect_identical(plot$cautions, review$Detail)
      expect_match(plot$subtitle, "Weak information")
      expect_null(plot_data(ci, subtitle = NULL)$subtitle)
      expect_identical(plot_data(ci, subtitle = NULL)$cautions, review$Detail)
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_match(paste(plot_data(as_ggplot(ci), "notes")$Text, collapse = "\n"), "Weak information")
      }
      file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
      saveRDS(res, file)
      expect_identical(readRDS(file)$facet_intervals$inference$cautions, review$Detail)
    }
  })

  test_that("equivalence cautions survive curated printing and extracted tables", {
    ordinary <- lapply(fits, analyze_facet_equivalence, facet = "Rater")
    local_mocked_bindings(compute_mml_parameter_covariance = reviewed)
    for (i in seq_along(fits)) {
      expect_warning(eq <- analyze_facet_equivalence(fits[[i]], facet = "Rater"), "Ill-conditioned")
      expect_equal(eq$pairwise[names(ordinary[[i]]$pairwise)], ordinary[[i]]$pairwise)
      expect_identical(eq$information_review, review)
      expect_true(review$Detail %in% summary(eq)$notes)
      expect_output(print(eq), "Ill-conditioned")
      expect_identical(unique(eq$pairwise$InferenceCaution), review$Detail)
      expect_match(plot(eq, draw = FALSE)$note, "Ill-conditioned")
      expect_identical(plot(eq, type = "rope", draw = FALSE)$information_review, review)
    }
  })

  test_that("failed information remains unavailable in fixed-facet consumers", {
    local_mocked_bindings(compute_mml_parameter_covariance = function(fit) {
      info <- original(fit); info$status <- "fallback"; info$cov <- NULL
      info$solution_information$inverse_review <- transform(review, Verified = FALSE)
      info
    })
    for (f in fits) {
      expect_warning(ci <- mfrm_facet_intervals(f, "Rater"), "Unregularized")
      expect_true(all(ci$table$Status == "information_unavailable"))
      expect_true(all(is.finite(ci$table$Estimate)))
      expect_true(all(is.na(ci$table$ModelSE)))
      expect_error(analyze_facet_equivalence(f, facet = "Rater"), "unregularized")
      expect_length(mfrm_mml_information_caution(compute_mml_parameter_covariance(f)), 0L)
    }
  })

  test_that("MI pooling preserves imputation identity and refuses failed information", {
    data <- ratings; data$Event <- paste0("E", seq_len(nrow(data)))
    completed <- list(data, data)
    data$Score[2] <- NA
    imputations <- mfrm_response_imputations(data, completed, "Person", c("Rater", "Criterion"),
      "Score", "Event", data$Event[2], 1:4,
      imputation_model = list(method = "consumer-contract fixture only"))
    for (f in fits) {
      # Reuse matching fitted completions; no extra optimizer runs are required.
      x <- structure(list(fits = list(f, f), imputations = imputations), class = "mfrm_imputed_fits")
      ordinary <- pool_mfrm_imputed(x, "Rater")
      local({
        calls <- 0L
        local_mocked_bindings(compute_mml_parameter_covariance = function(fit) {
          calls <<- calls + 1L
          if (calls == 2L) reviewed(fit) else original(fit)
        })
        expect_warning(out <- pool_mfrm_imputed(x, "Rater"), "Imputation\\(s\\) 2:")
        expect_equal(out$table[names(ordinary$table)], ordinary$table)
        expect_equal(out$total, ordinary$total)
        expect_length(out$information_review, 2)
        expect_null(out$information_review[[1]])
        expect_identical(out$information_review[[2]], review)
        expect_match(summary(out)$InferenceCaution[1], "Imputation\\(s\\) 2:")
        expect_output(print(out), "Imputation\\(s\\) 2:")
        expect_match(plot_data(out)$subtitle, "Weak information")
        expect_identical(plot_data(out)$cautions, out$cautions)
        expect_null(plot_data(out, title = NULL, subtitle = NULL)$subtitle)
        file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
        saveRDS(out, file); expect_identical(readRDS(file), out)
        # A failed completion must stop the entire pool, never omit a draw.
        local_mocked_bindings(compute_mml_parameter_covariance = function(fit) {
          info <- original(fit); info$status <- "fallback"; info$cov <- NULL; info
        })
        expect_error(pool_mfrm_imputed(x, "Rater"), "Imputation 1 requires")
      })
    }
  })

  test_that("diagnostic consumers retain details from the same information calculation", {
    for (f in fits) {
      info <- reviewed(f)
      expect_match(compute_mml_facet_model_se(f, info)$detail, "Ill-conditioned")
      expect_match(compute_mml_structural_parameter_se(f, covariance = info)$steps$SE_Detail[1],
        "Ill-conditioned")
    }
    local_mocked_bindings(compute_mml_parameter_covariance = reviewed)
    expect_identical(mfrmr_gqs_raw_information(fits[[1]])$covariance_detail, review$Detail)
  })

  test_that("GPCM diagnostic delta-method outputs retain the common caution", {
    f <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
      model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion")
    info <- reviewed(f)
    scorefile <- data.frame(Observed = f$prep$data$Score)
    ordinary <- add_gpcm_scorefile_delta_se(scorefile, f, covariance = original(f))
    scored <- add_gpcm_scorefile_delta_se(scorefile, f, covariance = info)
    expect_equal(scored$ExpectedScoreSE, ordinary$ExpectedScoreSE)
    expect_true(all(grepl("Ill-conditioned", scored$ScoreUncertaintyDetail)))
    raw <- suppressWarnings(fair_average_table(f))$raw_by_facet
    fair <- add_gpcm_fair_average_delta_se(raw, f, covariance = info)
    for (facet in c("Rater", "Criterion")) {
      expect_true(any(is.finite(fair[[facet]]$FairMSE)))
      expect_true(all(grepl("Ill-conditioned",
        fair[[facet]]$FairM_SE_Detail[is.finite(fair[[facet]]$FairMSE)])))
    }
  })

})
