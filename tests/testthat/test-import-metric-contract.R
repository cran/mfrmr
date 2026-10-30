test_that("TAM coordinates follow category logits rather than coefficient names", {
  source <- list(
    item = data.frame(item = c("Item_Cat1", "007")),
    xsi = data.frame(xsi = c(.2, .8, -.6, .4), se.xsi = rep(.2, 4)),
    A = array(0, c(2, 3, 4)), B = array(rep(0:2, each = 2), c(2, 3, 1)),
    se.B = array(0, c(2, 3, 1)),
    AXsi = rbind(c(0, -.2, -1), c(0, .6, .2))
  )
  source$A[1, 2, 1] <- source$A[1, 3, 1:2] <- -1
  source$A[2, 2, 3] <- source$A[2, 3, 3:4] <- -1
  extracted <- mfrmr:::.tam_extract_single(source, "Item")
  expect_identical(extracted$others$Level, c("Item_Cat1", "007"))
  expect_equal(extracted$others$Estimate, c(.5, -.1))
  expect_true(all(is.na(extracted$others$SE)))
  expect_equal(extracted$steps$Estimate, c(.2, .8, -.6, .4))

  # Equivalent item-location/step-deviation design: one marginal coefficient
  # now represents each item location, so its source SE can be retained.
  location_design <- source
  location_design$A[] <- 0
  location_design$A[1, 2, c(1, 3)] <- -1
  location_design$A[1, 3, 1] <- -2
  location_design$A[2, 2, c(2, 4)] <- -1
  location_design$A[2, 3, 2] <- -2
  location_design$xsi$xsi <- c(.5, -.1, -.3, -.5)
  alternate <- mfrmr:::.tam_extract_single(location_design, "Item")
  expect_identical(alternate$steps, extracted$steps)
  expect_equal(alternate$others$SE, c(.2, .2))

  scaled <- source
  scaled$B <- scaled$B * 2
  converted <- mfrmr:::.tam_extract_single(scaled, "Item")
  expect_equal(converted$others$Estimate, c(.25, -.05))
  for (i in 1:2) for (theta in c(-1, 0, 1)) {
    native <- exp(scaled$B[i, , 1] * theta + scaled$AXsi[i, ])
    steps <- converted$steps$Estimate[converted$steps$Level == source$item$item[i]]
    imported <- exp(c(0, cumsum(2 * (theta - steps))))
    expect_equal(imported / sum(imported), native / sum(native), tolerance = 1e-13)
  }
  location_design$se.B[] <- .1
  expect_true(all(is.na(mfrmr:::.tam_extract_single(location_design, "Item")$others$SE)))
  invalid <- source
  invalid$B[1, 3, 1] <- 3
  expect_error(mfrmr:::.tam_extract_single(invalid, "Item"), "constant adjacent-category slope")
})

test_that("TAM response flags use fitted support and retain Person alignment", {
  source <- list(
    pid = c("low", "high", "mixed", "missing", "partial", "unused", "masked"),
    item = data.frame(item = c("binary", "polytomous")),
    AXsi = rbind(c(0, -.2, NA), c(0, -.3, -.8)),
    resp = rbind(c(0, 0), c(1, 2), c(1, 0), c(NA, NA),
      c(NA, 2), c(1, 1), c(1, 0)))
  colnames(source$resp) <- source$item$item
  source$resp.ind <- !is.na(source$resp)
  source$resp.ind[7, 2] <- FALSE
  ids <- rev(source$pid)
  status <- mfrmr:::.tam_person_response_status(source, ids)
  expect_identical(status$Extreme, rev(c("low", "high", "none", NA,
    "high", "none", "high")))
  expect_identical(status$NObserved, rev(c(2L, 2L, 2L, 0L, 1L, 2L, 1L)))
  expect_equal(status$ResponseStatus[4], "no_observed_responses")
  # An observed ceiling of one cannot remove the fitted category two.
  source$resp[source$resp == 2 & !is.na(source$resp)] <- 1
  source$resp <- source$resp[, 2:1]
  source$resp.ind <- source$resp.ind[, 2:1]
  updated <- mfrmr:::.tam_person_response_status(source, source$pid)
  expect_equal(updated$Extreme[2], "none")
  expect_equal(updated$Extreme[5], "none")
  expect_equal(updated$Extreme[7], "high")
  source$resp[1, 1] <- 3
  expect_equal(mfrmr:::.tam_person_response_status(source, source$pid)$ResponseStatus[1], "unavailable")
  source$pid[1] <- source$pid[2]
  expect_true(all(is.na(mfrmr:::.tam_person_response_status(source, ids)$Extreme)))
  expect_true(all(is.na(mfrmr:::.tam_person_response_status(list(), ids)$Extreme)))
})

test_that("eRm locations and thresholds retain their difficulty meaning", {
  skip_if_not_installed("eRm")
  set.seed(47L)
  y <- matrix(sample(0:2, 120L, replace = TRUE), nrow = 30L)
  colnames(y) <- c("Item_Cat1", "007", "Item.3", "Item 4")
  for (model in c("PCM", "RSM", "RM")) {
    fit <- switch(model, PCM = eRm::PCM(y), RSM = eRm::RSM(y), RM = eRm::RM(y > 0))
    imported <- import_erm_fit(fit, model = if (model == "PCM") "PCM" else "RSM")
    expect_identical(imported$facets$others$Level, colnames(fit$X))
    expect_equal(nrow(imported$facets$others), ncol(fit$X))
    if (model != "RM") {
      native <- eRm::thresholds(fit)$threshtable[[1L]]
      expect_equal(unname(imported$facets$others$Estimate), unname(native[, "Location"]))
      expect_equal(imported$steps$Estimate, as.numeric(t(native[, -1L])))
      covariance <- fit$W %*% solve(fit$hessian) %*% t(fit$W)
      locations <- matrix(0, ncol(y), length(fit$betapar))
      locations[cbind(1:4, seq(2, 8, by = 2))] <- -.5
      expected_se <- sqrt(diag(locations %*% covariance %*% t(locations)))
      expect_equal(unname(imported$facets$others$SE), unname(expected_se))
      expect_error(import_erm_fit(fit, model = if (model == "PCM") "RSM" else "PCM"), "matching model label")
    } else {
      expect_equal(imported$facets$others$Estimate, unname(-fit$betapar))
      expect_equal(imported$facets$others$SE, unname(fit$se.beta))
    }
    expect_error(import_erm_fit(fit, model = "GPCM"), "not supported")
  }
})

test_that("mirt rating-scale offsets are included in absolute thresholds", {
  skip_if_not_installed("mirt")
  set.seed(47L)
  y <- matrix(sample(0:2, 120L, replace = TRUE), nrow = 30L)
  colnames(y) <- paste0("I", 1:4)
  rownames(y) <- paste0("Candidate", seq_len(nrow(y)))
  fit <- suppressWarnings(mirt::mirt(y, 1, itemtype = "rsm", verbose = FALSE))
  imported <- import_mirt_fit(fit, "RSM")
  expect_true(all(grepl("Absolute adjacent-category", imported$steps$Parameterization)))
  expect_identical(imported$facets$person$Person, paste0("P", seq_len(nrow(y))))
  expect_match(imported$source$person_identification, "cannot be recovered")
  native_scores <- mirt::fscores(fit, method = "EAP", full.scores.SE = TRUE, verbose = FALSE)
  expect_equal(imported$facets$person$Estimate, unname(native_scores[, 1]))
  expect_equal(imported$facets$person$SE, unname(native_scores[, 2]))
  theta <- matrix(c(-1, 0, 1), ncol = 1)
  for (i in 1:4) {
    native <- mirt::probtrace(mirt::extract.item(fit, i), theta)
    thresholds <- imported$steps$Estimate[imported$steps$Level == colnames(y)[i]]
    expected <- t(vapply(theta[, 1], function(value) {
      p <- exp(c(0, cumsum(value - thresholds)))
      p / sum(p)
    }, numeric(3)))
    expect_equal(unname(expected), unname(native), tolerance = 1e-12)
    expect_equal(imported$facets$others$Estimate[i], mean(thresholds))
  }
  multidimensional <- fit
  multidimensional@Model$nfact <- 2L
  expect_error(import_mirt_fit(multidimensional), "unidimensional")
  incompatible <- fit
  incompatible@Model$itemtype[] <- "graded"
  expect_error(import_mirt_fit(incompatible), "not interchangeable")
})

test_that("mixed mirt item types retain every item's thresholds and probabilities", {
  skip_if_not_installed("mirt")
  types <- c("Rasch", "2PL", "gpcm", "gpcmIRT", "gpcm", "Rasch")
  counts <- c(2L, 2L, 3L, 4L, 2L, 3L)
  y <- vapply(counts, function(k) rep(seq_len(k) - 1L, length.out = 60L), integer(60))
  colnames(y) <- paste0("Item", seq_along(types))
  # Fixed source parameters isolate the import mapping from estimation quality.
  pars <- mirt::mirt(y, 1, itemtype = types, pars = "values", verbose = FALSE)
  pars$value[pars$name == "a1"] <- c(1, 1.3, .8, 1.2, .9, 1)
  pars$value[pars$name == "c"] <- .35
  pars$est <- FALSE
  fit <- mirt::mirt(y, 1, itemtype = types, pars = pars, verbose = FALSE)
  imported <- import_mirt_fit(fit, model = "GPCM")
  theta <- matrix(c(-2, -.3, 0, 1.7), ncol = 1L)
  expect_equal(nrow(imported$steps), sum(counts - 1L))
  expect_identical(imported$source$item_types, types)
  for (i in seq_along(types)) {
    rows <- imported$steps[imported$steps$Level == colnames(y)[i], ]
    expect_identical(rows$Step, seq_len(counts[i] - 1L))
    expect_true(all(is.finite(rows$Estimate)))
    a <- imported$facets$others$Slope[i]
    expected <- t(vapply(theta[, 1L], function(t) {
      logits <- c(0, cumsum(a * (t - rows$Estimate)))
      p <- exp(logits - max(logits)); p / sum(p)
    }, numeric(counts[i])))
    native <- mirt::probtrace(mirt::extract.item(fit, i), theta)
    expect_equal(unname(expected), unname(native), tolerance = 1e-12)
    expect_equal(imported$facets$others$Estimate[i], mean(rows$Estimate))
  }
  expect_error(import_mirt_fit(fit, model = "PCM"), "Non-unit source item slopes")
  source <- mirt::coef(fit, simplify = TRUE, IRTpars = TRUE)$items
  actual <- source
  # An inapplicable NA must be ignored, but a missing actual threshold must fail.
  actual[3L, "b1"] <- NA_real_
  local_mocked_bindings(.mfrm_mirt_coef = function(...) list(items = actual), .package = "mfrmr")
  expect_error(import_mirt_fit(fit, model = "GPCM"), "finite source adjacent-category")
  actual <- source
  actual[4L, "a1"] <- NA_real_
  expect_error(import_mirt_fit(fit, model = "GPCM"), "finite positive source slope")
  actual <- source
  actual[4L, "c"] <- NA_real_
  expect_error(import_mirt_fit(fit, model = "GPCM"), "rating-scale offsets are unavailable")
})

test_that("imported summaries and plots do not invent native model assumptions", {
  skip_if_not_installed("TAM")
  set.seed(44L)
  y <- matrix(stats::rbinom(80L * 6L, 1, .5), nrow = 80L)
  colnames(y) <- paste0("I", seq_len(ncol(y)))
  y[1, ] <- 0
  y[2, ] <- 1
  y[1, 1] <- NA
  fit <- TAM::tam.mml(y, verbose = FALSE)
  imported <- import_tam_fit(fit, "RSM", compute_fit = TRUE)
  expect_equal(imported$diagnostics$n_obs, sum(!is.na(y)))
  expect_identical(imported$facets$person$Extreme[1:2], c("low", "high"))
  expect_identical(imported$facets$person$NObserved[1:2], c(5L, 6L))
  expect_identical(imported$facets$person$Estimate, fit$person$EAP)
  expect_identical(imported$facets$person$SE, fit$person$SD.EAP)
  masked <- fit
  masked$resp[is.na(masked$resp)] <- 0
  masked_import <- import_tam_fit(masked, "RSM", compute_fit = TRUE)
  expect_equal(masked_import$diagnostics$n_obs, sum(!is.na(y)))
  expect_identical(masked_import$facets$person, imported$facets$person)
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(imported, saved)
  expect_identical(readRDS(saved), imported)
  person_rel <- imported$diagnostics$reliability
  expect_true(is.na(person_rel$Reliability[person_rel$Facet == "Person"]))
  expect_match(person_rel$SummaryNote[person_rel$Facet == "Person"], "Posterior SD")
  expect_true(all(imported$facets$person$FitBasis == "Source WLE person fit; estimates remain EAP"))
  expect_true(all(imported$diagnostics$fit$FitBasis[imported$diagnostics$fit$Facet == "Person"] == "Source WLE person fit; estimates remain EAP"))
  overview <- summary(imported)
  expect_s3_class(overview, "summary.mfrm_imported_fit")
  expect_identical(overview$scale_contract$CoordinateBasis, "source_package_scale")
  expect_true(is.na(overview$scale_contract$PopulationSD))
  displayed <- paste(capture.output({print(imported); print(overview); print(summary(imported$diagnostics))}), collapse = "\n")
  expect_false(grepl("legacy_unknown|optimizer_polish|fixed N\\(0,1\\)|imported_tam_", displayed))
  expect_match(displayed, "source ability scale")
  expect_warning(map <- plot(imported, draw = FALSE), "Imported point estimates")
  expect_identical(map$data$scale_contract$CoordinateBasis, "source_package_scale")
  expect_identical(map$data$axis_label, "Source ability scale")
  expect_identical(map$data$reference_lines$label, "Source scale zero")
  expect_false(map$data$show_ci)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_identical(as_ggplot(map)$labels$y, "Source ability scale")
  }
  expect_error(plot(imported, show_ci = TRUE, draw = FALSE), "shared confidence-interval")
  expect_error(plot(imported, type = "ccc", draw = FALSE), "model-dependent")
  expect_error(diagnose_mfrm(imported), "Native response-level")
  expect_error(mfrm_results(imported, compute = "never"), "Comprehensive mfrmr reports require a native fit")
  old <- imported
  old$source$metric_version <- NULL
  expect_error(summary(old), "Re-import")
  expect_error(plot(old, draw = FALSE), "Re-import")
  old_diag <- imported$diagnostics
  old_diag$metric_version <- NULL
  expect_error(summary(old_diag), "Re-import")
})
