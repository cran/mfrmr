# Tests for DIF analysis module (Phase 2) and compare_mfrm enhancements (Phase 3)

# ---------- shared fixtures ----------
local_dif_fixtures <- function(env = parent.frame()) {
  toy <- load_mfrmr_data("study1")
  persons <- unique(toy$Person)
  half <- ceiling(length(persons) / 2)
  grp_map <- setNames(
    c(rep("A", half), rep("B", length(persons) - half)),
    persons
  )
  toy$Group <- grp_map[toy$Person]

  fit <- fit_mfrm(toy, person = "Person", facets = c("Rater", "Criterion"),
                  score = "Score", method = "JML")
  diag <- diagnose_mfrm(fit, residual_pca = "none")

  assign("toy",  toy,  envir = env)
  assign("fit",  fit,  envir = env)
  assign("diag", diag, envir = env)
}

# ================================================================
# Phase 2: DIF diagnostic module
# ================================================================

test_that("analyze_dff residual method returns expected structure", {
  local_dif_fixtures()

  dif <- analyze_dff(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")

  expect_s3_class(dif, "mfrm_dff")
  expect_s3_class(dif, "mfrm_dif")
  expect_true(is.data.frame(dif$dif_table))
  expect_true(is.data.frame(dif$cell_table))
  expect_true(nrow(dif$dif_table) > 0)
  expect_true(nrow(dif$cell_table) > 0)

  # Required columns in dif_table
  required_cols <- c("Level", "Group1", "Group2", "Contrast", "SE", "t",
                     "df", "p_value", "Classification", "ClassificationSystem",
                     "ETS", "Method", "SEBasis", "StatisticLabel",
                     "ProbabilityMetric", "DFBasis", "ContrastBasis",
                     "ReportingUse", "PrimaryReportingEligible")
  expect_true(all(required_cols %in% names(dif$dif_table)))

  # Method must be "residual"
  expect_true(all(dif$dif_table$Method == "residual"))
  expect_true(all(dif$dif_table$ClassificationSystem == "screening"))
  expect_true(all(dif$dif_table$StatisticLabel == "Welch screening t"))
  expect_true(all(dif$dif_table$DFBasis == "Welch-Satterthwaite approximation"))
  expect_true(all(is.na(dif$dif_table$ETS)))
  expect_true(all(dif$dif_table$ReportingUse == "screening_only"))
  expect_true(all(!dif$dif_table$PrimaryReportingEligible))
  expect_equal(dif$config$method, "residual")
})

test_that("analyze_dff alias is backward compatible with analyze_dif", {
  local_dif_fixtures()

  dff <- analyze_dff(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")
  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")

  expect_equal(dff$dif_table, dif$dif_table)
  expect_equal(dff$summary, dif$summary)
  expect_equal(dff$config$functioning_label, "DIF")
})

test_that("analyze_dif refit keeps JML contrasts descriptive even when linked", {
  local_dif_fixtures()

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "refit")

  expect_s3_class(dif, "mfrm_dif")
  expect_true(is.data.frame(dif$dif_table))
  expect_true(nrow(dif$dif_table) > 0)
  expect_true(all(c("ContrastComparable", "ScaleLinkStatus", "LinkingFacets",
                    "FormalInferenceEligible", "InferenceTier",
                    "ReportingUse", "PrimaryReportingEligible",
                    "BaselineConverged", "SubgroupConverged1", "SubgroupConverged2",
                    "BaselineMethod", "BaselinePrecisionTier",
                    "BaselineSupportsFormalInference") %in% names(dif$dif_table)))
  expect_true(all(c("LinkingStatus", "LinkingDetail", "ETS_Eligible",
                    "LinkComparable", "PrecisionTier", "SupportsFormalInference") %in% names(dif$group_fits[[1]])))
  expect_true("Converged" %in% names(dif$group_fits[[1]]))
  expect_true(all(dif$dif_table$ClassificationSystem == "descriptive"))
  expect_true(all(dif$dif_table$StatisticLabel == "linked descriptive contrast"))
  expect_true(all(dif$dif_table$SEBasis == "not reported without model-based subgroup precision"))
  expect_true(all(is.na(dif$dif_table$ETS)))
  expect_true(all(dif$dif_table$ContrastComparable))
  expect_true(all(!dif$dif_table$FormalInferenceEligible))
  expect_true(all(!dif$dif_table$PrimaryReportingEligible))
  expect_true(all(dif$dif_table$InferenceTier == "exploratory"))
  expect_true(all(dif$dif_table$ScaleLinkStatus == "linked"))
  expect_true(all(dif$dif_table$ReportingUse == "screening_only"))
  expect_equal(dif$config$method, "refit")
})

test_that("analyze_dif refit surfaces subgroup diagnostics failures", {
  local_dif_fixtures()

  testthat::local_mocked_bindings(
    diagnose_mfrm = function(...) stop("forced subgroup diagnostics failure"),
    .package = "mfrmr"
  )

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "refit")

  group_fit_tbl <- dif$group_fits[[1]]
  expect_true(all(c("DiagnosticsStatus", "DiagnosticsDetail") %in% names(group_fit_tbl)))
  expect_true(all(group_fit_tbl$DiagnosticsStatus == "failed"))
  expect_true(all(group_fit_tbl$PrecisionTier == "diagnostics_unavailable"))
  expect_true(all(grepl("forced subgroup diagnostics failure", group_fit_tbl$DiagnosticsDetail, fixed = TRUE)))
})

test_that("analyze_dif refit uses ETS only for linked model-based MML contrasts", {
  dat <- mfrmr:::sample_mfrm_data(seed = 42)
  persons <- unique(dat$Person)
  half <- ceiling(length(persons) / 2)
  grp_map <- setNames(
    c(rep("A", half), rep("B", length(persons) - half)),
    persons
  )
  dat$Group <- grp_map[dat$Person]

  fit_mml <- suppressWarnings(fit_mfrm(
    dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    quad_points = 7,
    maxit = 25
  ))
  diag_mml <- diagnose_mfrm(fit_mml, residual_pca = "none")
  dif_mml <- analyze_dif(
    fit_mml,
    diag_mml,
    facet = "Criterion",
    group = "Group",
    data = dat,
    method = "refit"
  )

  expect_true(any(dif_mml$dif_table$ClassificationSystem == "ETS"))
  ets_rows <- dif_mml$dif_table$ClassificationSystem == "ETS"
  expect_true(all(dif_mml$dif_table$ContrastComparable[ets_rows]))
  expect_true(all(dif_mml$dif_table$FormalInferenceEligible[ets_rows]))
  expect_true(all(dif_mml$dif_table$PrimaryReportingEligible[ets_rows]))
  expect_true(all(dif_mml$dif_table$InferenceTier[ets_rows] == "model_based"))
  expect_true(all(dif_mml$dif_table$ReportingUse[ets_rows] == "primary_reporting"))
  expect_true(all(stats::na.omit(dif_mml$dif_table$ETS) %in% c("A", "B", "C")))
})

test_that("analyze_dif refit demotes ETS when subgroup refits lack linking facets", {
  local_dif_fixtures()

  fit_one <- fit_mfrm(toy, person = "Person", facets = "Criterion",
                      score = "Score", method = "JML", maxit = 20)
  diag_one <- diagnose_mfrm(fit_one, residual_pca = "none")
  dif_one <- analyze_dif(fit_one, diag_one, facet = "Criterion", group = "Group",
                         data = toy, method = "refit")

  expect_true(all(dif_one$dif_table$ClassificationSystem == "descriptive"))
  expect_true(all(is.na(dif_one$dif_table$ETS)))
  expect_true(all(dif_one$dif_table$ScaleLinkStatus == "unlinked"))
  expect_true(all(!dif_one$dif_table$ContrastComparable))
  expect_match(dif_one$summary$Classification[nrow(dif_one$summary)], "insufficient linking", ignore.case = TRUE)
})

test_that("analyze_dif min_obs filter works", {
  local_dif_fixtures()

  # With a very high min_obs, all cells should be sparse
  dif_high <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                          data = toy, method = "residual", min_obs = 99999)
  expect_true(all(dif_high$cell_table$sparse))
})

test_that("analyze_dif p_adjust works for all methods", {
  local_dif_fixtures()

  for (m in c("holm", "fdr", "bonferroni", "none")) {
    dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                       data = toy, method = "residual", p_adjust = m)
    expect_true("p_adjusted" %in% names(dif$dif_table))
  }
})

test_that("residual method uses screening labels instead of ETS categories", {
  local_dif_fixtures()

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")
  expect_true(all(dif$summary$Classification %in% c("Screen positive", "Screen negative", "Unclassified")))
  expect_true(all(is.na(dif$dif_table$ETS)))
})

test_that("residual method uses Welch-Satterthwaite degrees of freedom", {
  local_dif_fixtures()

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")
  first_row <- dif$dif_table[1, , drop = FALSE]
  c1 <- subset(dif$cell_table, Level == first_row$Level & GroupValue == first_row$Group1)
  c2 <- subset(dif$cell_table, Level == first_row$Level & GroupValue == first_row$Group2)
  comp1 <- c1$Var_sum[1] / c1$N[1]^2
  comp2 <- c2$Var_sum[1] / c2$N[1]^2
  expected_df <- (comp1 + comp2)^2 / ((comp1^2) / (c1$N[1] - 1) + (comp2^2) / (c2$N[1] - 1))

  expect_equal(first_row$df, expected_df, tolerance = 1e-8)
})

test_that("refit method ETS classification is valid", {
  dat <- mfrmr:::sample_mfrm_data(seed = 99)
  dat$Group <- ifelse(dat$Person %in% unique(dat$Person)[1:30], "A", "B")
  fit_mml <- suppressWarnings(fit_mfrm(
    dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    quad_points = 7,
    maxit = 25
  ))
  diag_mml <- diagnose_mfrm(fit_mml, residual_pca = "none")
  dif <- analyze_dif(fit_mml, diag_mml, facet = "Criterion", group = "Group",
                     data = dat, method = "refit")
  expect_true(all(stats::na.omit(dif$dif_table$ETS) %in% c("A", "B", "C")))
})

test_that("dif_interaction_table returns expected structure", {
  local_dif_fixtures()

  int <- dif_interaction_table(fit, diag, facet = "Criterion", group = "Group",
                               data = toy, min_obs = 2)

  expect_s3_class(int, "mfrm_dif_interaction")
  expect_true(is.data.frame(int$table))
  expect_true(nrow(int$table) > 0)

  required_cols <- c("Level", "GroupValue", "N", "ObsScore", "ExpScore",
                     "ObsExpAvg")
  expect_true(all(required_cols %in% names(int$table)))
})

test_that("dif_interaction_table min_obs filter works", {
  local_dif_fixtures()

  int <- dif_interaction_table(fit, diag, facet = "Criterion", group = "Group",
                               data = toy, min_obs = 99999)
  expect_true(all(int$table$sparse))
})

test_that("plot_dif_heatmap returns matrix when draw = FALSE", {
  local_dif_fixtures()

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")

  for (m in c("obs_exp", "t", "contrast")) {
    mat <- plot_dif_heatmap(dif, metric = m, draw = FALSE)
    expect_true(is.matrix(mat))
  }
})

test_that("plot_dif_heatmap works with dif_interaction_table", {
  local_dif_fixtures()

  int <- dif_interaction_table(fit, diag, facet = "Criterion", group = "Group",
                               data = toy, min_obs = 2)

  mat <- plot_dif_heatmap(int, metric = "obs_exp", draw = FALSE)
  expect_true(is.matrix(mat))
})

test_that("dif_report produces interpretable output", {
  local_dif_fixtures()

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")

  rpt <- dif_report(dif)
  expect_s3_class(rpt, "mfrm_dif_report")
  expect_true(is.character(rpt$narrative))
  expect_true(nchar(rpt$narrative) > 0)
  expect_match(rpt$narrative, "screening", ignore.case = TRUE)
})

test_that("print and summary S3 methods work for DIF objects", {
  local_dif_fixtures()

  dif <- analyze_dif(fit, diag, facet = "Criterion", group = "Group",
                     data = toy, method = "residual")
  expect_output(print(dif))
  s <- summary(dif)
  expect_output(print(s))

  int <- dif_interaction_table(fit, diag, facet = "Criterion", group = "Group",
                               data = toy, min_obs = 2)
  expect_output(print(int))
  s2 <- summary(int)
  expect_output(print(s2))

  rpt <- dif_report(dif)
  expect_output(print(rpt))
  s3 <- summary(rpt)
  expect_output(print(s3))
})


# ================================================================
# Phase 3: compare_mfrm enhancements
# ================================================================

test_that("compare_mfrm reports comparable IC quantities on a common basis", {
  toy_small <- load_mfrmr_data("example_core")

  fit_rsm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_pcm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "PCM",
    step_facet = "Criterion",
    quad_points = 5,
    maxit = 15
  ))
  fit_rsm$summary$Converged[1] <- TRUE
  fit_pcm$summary$Converged[1] <- TRUE

  comp <- compare_mfrm(RSM = fit_rsm, PCM = fit_pcm)

  expect_s3_class(comp, "mfrm_comparison")
  tbl <- comp$table
  expect_true(all(tbl$ICComparable))
  expect_true("Delta_AIC" %in% names(tbl))
  expect_true("Delta_BIC" %in% names(tbl))
  expect_true("AkaikeWeight" %in% names(tbl))
  expect_true("BICWeight" %in% names(tbl))

  # Delta should have at least one zero (best model)
  expect_equal(min(tbl$Delta_AIC), 0)
  expect_equal(min(tbl$Delta_BIC), 0)

  # Weights should sum to 1
  expect_equal(sum(tbl$AkaikeWeight), 1, tolerance = 1e-10)
  expect_equal(sum(tbl$BICWeight), 1, tolerance = 1e-10)
  expect_true(isTRUE(comp$comparison_basis$ic_comparable))
  expect_true(isTRUE(comp$comparison_basis$same_data))
})

test_that("compare_mfrm evidence_ratios are reciprocal", {
  toy_small <- load_mfrmr_data("example_core")

  fit_rsm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_pcm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "PCM",
    step_facet = "Criterion",
    quad_points = 5,
    maxit = 15
  ))
  fit_rsm$summary$Converged[1] <- TRUE
  fit_pcm$summary$Converged[1] <- TRUE

  comp <- compare_mfrm(RSM = fit_rsm, PCM = fit_pcm)

  er <- comp$evidence_ratios
  expect_true(is.data.frame(er))
  expect_true(nrow(er) > 0)
  expect_true("EvidenceRatio" %in% names(er))
  expect_true(all(er$EvidenceRatio > 0))
})

test_that("compare_mfrm print and summary work with new fields", {
  toy_small <- load_mfrmr_data("example_core")

  fit_rsm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_pcm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "PCM",
    step_facet = "Criterion",
    quad_points = 5,
    maxit = 15
  ))
  fit_rsm$summary$Converged[1] <- TRUE
  fit_pcm$summary$Converged[1] <- TRUE

  comp <- compare_mfrm(RSM = fit_rsm, PCM = fit_pcm)

  expect_output(print(comp))
  s <- summary(comp)
  expect_output(print(s))
})

test_that("compare_mfrm suppresses IC ranking outside the formal MML path and requires nested = TRUE for LRT", {
  toy_small <- load_mfrmr_data("example_core")

  fit_jml <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 15
  ))
  fit_mml <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_jml$summary$Converged[1] <- TRUE
  fit_mml$summary$Converged[1] <- TRUE

  expect_warning(
    comp <- compare_mfrm(JML = fit_jml, MML = fit_mml),
    "different estimation methods"
  )

  expect_false(isTRUE(comp$comparison_basis$ic_comparable))
  expect_true(all(is.na(comp$table$Delta_AIC)))
  expect_true(all(is.na(comp$table$AkaikeWeight)))
  expect_null(comp$lrt)

  fit_pcm <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    model = "PCM",
    step_facet = "Criterion",
    maxit = 15
  ))
  fit_pcm$summary$Converged[1] <- TRUE

  expect_warning(
    comp_no_lrt <- compare_mfrm(RSM = fit_jml, PCM = fit_pcm),
    "limited to converged MML fits"
  )
  expect_null(comp_no_lrt$lrt)

  expect_warning(
    comp_lrt <- compare_mfrm(RSM = fit_jml, PCM = fit_pcm, nested = TRUE),
    "formal MML likelihood basis"
  )
  expect_true(isTRUE(comp_lrt$comparison_basis$nested_requested))
  expect_null(comp_lrt$lrt)
  expect_true(isTRUE(comp_lrt$comparison_basis$nesting_audit$eligible))
  expect_identical(as.character(comp_lrt$comparison_basis$nesting_audit$relation), "RSM_in_PCM")

  fit_rsm_2 <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 15
  ))
  fit_rsm_2$summary$Converged[1] <- TRUE

  expect_warning(
    comp_same <- compare_mfrm(RSM1 = fit_jml, RSM2 = fit_rsm_2, nested = TRUE),
    "formal MML likelihood basis"
  )
  expect_null(comp_same$lrt)
  expect_false(isTRUE(comp_same$comparison_basis$ic_comparable))
  expect_false(isTRUE(comp_same$comparison_basis$nesting_audit$eligible))
  expect_identical(as.character(comp_same$comparison_basis$nesting_audit$relation), "same_model")
})

test_that("compare_mfrm suppresses IC ranking for JML-only comparisons", {
  local_dif_fixtures()

  fit2 <- fit_mfrm(toy, person = "Person", facets = c("Rater", "Criterion"),
                   score = "Score", method = "JML",
                   model = "PCM", step_facet = "Criterion")
  fit$summary$Converged[1] <- TRUE
  fit2$summary$Converged[1] <- TRUE

  expect_warning(
    comp <- compare_mfrm(RSM = fit, PCM = fit2),
    "limited to converged MML fits"
  )

  expect_false(isTRUE(comp$comparison_basis$ic_comparable))
  expect_false(isTRUE(comp$comparison_basis$all_mml))
  expect_true(all(is.na(comp$table$Delta_AIC)))
  expect_true(all(is.na(comp$table$AkaikeWeight)))
  expect_null(comp$evidence_ratios)
})

test_that("compare_mfrm suppresses IC ranking when a fit is marked unconverged", {
  local_dif_fixtures()

  fit2 <- fit_mfrm(toy, person = "Person", facets = c("Rater", "Criterion"),
                   score = "Score", method = "JML",
                   model = "PCM", step_facet = "Criterion")
  fit2$summary$Converged[1] <- FALSE

  expect_warning(
    comp <- compare_mfrm(RSM = fit, PCM = fit2),
    "did not converge"
  )

  expect_false(isTRUE(comp$comparison_basis$ic_comparable))
  expect_false(isTRUE(comp$comparison_basis$all_converged))
  expect_true(all(is.na(comp$table$Delta_AIC)))
  expect_true(all(is.na(comp$table$AkaikeWeight)))
})

test_that("compare_mfrm requires the same prepared response data for IC ranking", {
  toy_small <- load_mfrmr_data("example_core")

  fit_a <- suppressWarnings(fit_mfrm(
    toy_small,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_a$summary$Converged[1] <- TRUE

  toy_perm <- toy_small[sample.int(nrow(toy_small)), , drop = FALSE]
  fit_perm <- suppressWarnings(fit_mfrm(
    toy_perm,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_perm$summary$Converged[1] <- TRUE

  comp_same <- compare_mfrm(RSM1 = fit_a, RSM2 = fit_perm)
  expect_true(isTRUE(comp_same$comparison_basis$same_data))
  expect_true(all(comp_same$table$ICComparable))

  toy_shuf <- toy_small
  toy_shuf$Score <- sample(toy_shuf$Score)
  fit_shuf <- suppressWarnings(fit_mfrm(
    toy_shuf,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 15
  ))
  fit_shuf$summary$Converged[1] <- TRUE

  expect_warning(
    comp_diff <- compare_mfrm(A = fit_a, B = fit_shuf),
    "same prepared response data"
  )
  expect_false(isTRUE(comp_diff$comparison_basis$same_data))
  expect_false(any(comp_diff$table$ICComparable))
  expect_true(all(is.na(comp_diff$table$Delta_AIC)))
  expect_true(all(is.na(comp_diff$table$AkaikeWeight)))
})
