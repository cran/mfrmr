population_interval_fixture <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      fixture <- readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))
      saved <<- list(fixture=fixture,
        native=suppressWarnings(mfrm_facet_intervals(fixture$fit, "Rater")))
    }
    saved
  }
})

test_that("unanchored native interval sources accept both centering conventions", {
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    mfrm_native_location_inference = function(...) stop("Reached numerical qualification"),
    .package = "mfrmr")
  for (name in c("rsm-zero-boundary", "pcm-regression", "fixed-gpcm")) {
    saved <- readRDS(test_path("fixtures", paste0("native-location-", name, ".rds")))
    legacy <- if (inherits(saved, "mfrm_fit")) saved else saved$fit
    expected <- mfrm_native_location_source(legacy)
    current <- legacy
    current$config$centering_scope <- "facet"
    for (facet in current$config$facet_names) {
      current$config$facet_specs[[facet]] <- build_facet_constraint(
        current$prep$levels[[facet]], centering_scope = "facet")
    }
    path <- tempfile(fileext = ".rds")
    saveRDS(current, path)
    current <- readRDS(path)
    unlink(path)
    expect_identical(mfrm_native_location_source(current), expected)
    expect_error(mfrm_facet_intervals(current, "Rater"), "Reached numerical qualification")
    if (current$config$model != "GPCM") {
      expect_error(mfrm_population_intervals(current), "Reached numerical qualification")
    }
    # Equivalent metadata must not admit different constraint equations.
    for (fit in list(legacy, current)) {
      anchored <- fit
      anchored$config$facet_specs$Rater$anchors[1] <- 0
      expect_error(mfrm_native_location_source(anchored), "without anchors or groups")
      grouped <- fit
      grouped$config$facet_specs$Rater$groups[1:2] <- "g"
      expect_error(mfrm_native_location_source(grouped), "without anchors or groups")
      uncentered <- fit
      uncentered$config$facet_specs$Rater$centered <- FALSE
      expect_error(mfrm_native_location_source(uncentered), "sum-to-zero")
    }
    changed <- current
    changed$config$facet_specs$Rater$remainder_sum <- .1
    expect_error(mfrm_native_location_source(changed), "sum-to-zero")
  }
})

test_that("population intervals reuse full information in the declared coordinates", {
  x <- population_interval_fixture()
  local_mocked_bindings(fit_mfrm=function(...) stop("No refitting"),
    mfrm_native_location_inference=function(...) stop("No repeated qualification"),
    compute_mml_parameter_covariance=function(...) stop("No repeated information"), .package="mfrmr")
  ci <- suppressWarnings(mfrm_population_intervals(x$native, level=.9))
  expect_identical(ci$fit, x$fixture$fit)
  expect_true(ci$settings$inference_reused)
  expect_true(all(ci$table$CIEligible))
  expected <- x$fixture$reference_covariance[11:13,11:13]
  expect_equal(vcov(ci), expected, tolerance=1e-5, ignore_attr=TRUE)
  expect_identical(rownames(vcov(ci)), c("beta:(Intercept)", "beta:X", "log_residual_variance"))
  tab <- summary(ci); critical <- qnorm(.95)
  expect_equal(tab$Lower[1:2],tab$Estimate[1:2]-critical*tab$WorkingSE[1:2])
  expect_equal(tab$Lower[3],exp(log(tab$Estimate[3])-critical*tab$WorkingSE[3]))
  expect_equal(tab$Upper[3],exp(log(tab$Estimate[3])+critical*tab$WorkingSE[3]))
  expect_equal(tab$SE[3],tab$Estimate[3]*tab$WorkingSE[3])
  expect_gt(abs(tab$Lower[3]-(tab$Estimate[3]-critical*tab$SE[3])),1e-3)
  expect_equal(tab$Estimate[3],x$fixture$fit$population$sigma2)
  expect_false(any(startsWith(ci$cautions, "Native location intervals:")))
  legacy <- x$native
  legacy$cautions <- sub("^Native location intervals:", "Experimental native location intervals:",
    legacy$cautions)
  expect_identical(suppressWarnings(mfrm_population_intervals(legacy, level = .9))$cautions,
    ci$cautions)
})

test_that("direct population requests use the same numerical qualification and preserve refusal", {
  x <- population_interval_fixture(); calls <- 0L; refused <- FALSE
  local_mocked_bindings(mfrm_native_location_inference=function(fit) {
    calls <<- calls+1L
    expect_identical(fit,x$fixture$fit)
    list(covariance=if(refused) NULL else x$native$parameter_covariance,
      checks=data.frame(Check="Injected numerical check",Passed=!refused,Detail="Finer information unavailable"),
      check=list(eligible=!refused,review="Finer information unavailable"))
  }, .package="mfrmr")
  ci <- suppressWarnings(mfrm_population_intervals(x$fixture$fit))
  expect_false(ci$settings$inference_reused)
  expect_true(all(ci$table$CIEligible))
  refused <- TRUE
  failed <- suppressWarnings(mfrm_population_intervals(x$fixture$fit))
  expect_equal(calls,2L)
  expect_equal(failed$table$Estimate,ci$table$Estimate)
  expect_true(all(!failed$table$CIEligible))
  expect_true(all(is.na(failed$table$SE) & is.na(failed$table$Lower) & is.na(failed$table$Upper)))
  expect_true(all(is.na(vcov(failed))))
  expect_match(failed$table$InferenceReview[1],"Finer information unavailable")
})

test_that("population covariance reuse requires a consistent qualified source", {
  x <- population_interval_fixture(); source <- x$native
  local_mocked_bindings(compute_mml_parameter_covariance=function(...) stop("No information"), .package="mfrmr")
  for (level in list(0,1,NA_real_,c(.9,.95),1i))
    expect_error(mfrm_population_intervals(source,level),"0 < level < 1")
  bad <- source; bad$source$engine <- "changed"
  expect_error(mfrm_population_intervals(bad),"must match")
  bad <- source; bad$numerical_checks <- NULL
  expect_error(mfrm_population_intervals(bad),"complete numerical qualification")
  bad <- source; bad$parameter_covariance[1,1] <- -1
  expect_error(mfrm_population_intervals(bad),"positive definite")
  bad <- source$fit; bad$config$model <- "GPCM"
  expect_error(mfrm_population_intervals(bad),"native RSM/PCM")
})

test_that("a better zero-variance candidate prevents population and facet intervals", {
  x <- readRDS(test_path("fixtures", "native-location-rsm-zero-boundary.rds"))
  local_mocked_bindings(fit_mfrm=function(...) stop("No refitting"), .package="mfrmr")
  # Complete pattern frequencies attain the saturated maximum at variance zero.
  # The historical local-only checks nevertheless admitted this small positive fit.
  expect_true(all(x$legacy_population$table$CIEligible))
  expect_gt(x$fit$population$sigma2, 0)
  ci <- suppressWarnings(mfrm_population_intervals(x$fit))
  facet <- suppressWarnings(mfrm_facet_intervals(x$fit, "Rater"))
  for (result in list(ci, facet)) {
    expect_identical(result$fit, x$fit)
    expect_true(all(!result$table$CIEligible))
    expect_true(all(is.finite(result$table$Estimate)))
    expect_true(all(is.na(result$table$SE) & is.na(result$table$Lower) & is.na(result$table$Upper)))
    expect_match(result$table$InferenceReview[1], "zero-residual-variance candidate")
    # The retained nuisance estimates differ slightly from the known zero optimum.
    expect_equal(result$numerical_checks$ZeroVarianceCandidateNLL, x$zero_candidate_nll, tolerance=1e-12)
    expect_gte(result$numerical_checks$ZeroVarianceCandidateNLL, x$zero_nll)
    expect_false(result$numerical_checks$BoundaryCandidatePassed)
    expect_gt(result$numerical_checks$BoundaryCandidateGain,
      result$numerical_checks$BoundaryCandidateTolerance)
  }
  # Reconstruct the saved facet-inference fields from the historical full
  # covariance and checks; an old reusable covariance must not bypass the guard.
  old <- x$legacy_population
  legacy <- old[c("fit", "source", "checks", "numerical_checks", "parameter_covariance",
    "information_review", "comparison_information_review", "cautions")]
  legacy$settings <- list(procedure="mml_native_location_model_v1")
  class(legacy) <- "mfrm_facet_intervals"
  local_mocked_bindings(mfrm_native_location_inference=function(...) stop("No repeated qualification"),
    compute_mml_parameter_covariance=function(...) stop("No differentiation"), .package="mfrmr")
  reused <- suppressWarnings(mfrm_population_intervals(legacy))
  expect_identical(reused$table, ci$table)
  # Reuse preserves the historical curvature diagnostics exactly. A fresh
  # Hessian on another platform need not have bitwise-identical roundoff.
  saved_checks <- names(old$numerical_checks)
  expect_identical(reused$numerical_checks[saved_checks], old$numerical_checks)
  boundary_checks <- setdiff(names(ci$numerical_checks), saved_checks)
  expect_identical(reused$numerical_checks[boundary_checks], ci$numerical_checks[boundary_checks])
  expect_true(all(is.na(vcov(reused))))
  expect_true(reused$settings$inference_reused)
})

test_that("an unavailable boundary candidate comparison cannot admit saved covariance", {
  x <- population_interval_fixture()
  local_mocked_bindings(make_mfrm_direct_evaluator=function(...) stop("Injected likelihood failure"),
    compute_mml_parameter_covariance=function(...) stop("No repeated information"), .package="mfrmr")
  ci <- suppressWarnings(mfrm_population_intervals(x$native))
  expect_true(all(!ci$table$CIEligible))
  expect_true(all(is.na(ci$table$Lower)))
  expect_match(ci$table$InferenceReview[1], "Injected likelihood failure")
})

test_that("the boundary candidate check does not mistake roundoff for likelihood improvement", {
  x <- readRDS(test_path("fixtures", "native-location-rsm-zero-boundary.rds"))
  old <- x$legacy_population
  inference <- list(covariance=old$parameter_covariance,checks=old$checks,
    numerical_checks=old$numerical_checks,check=list(eligible=TRUE))
  # Probe the likelihood limit only; this does not qualify information or
  # normal intervals at this altered, nearly degenerate parameter vector.
  near <- x$fit; near$opt$par[length(near$opt$par)] <- -32
  checked <- mfrm_native_location_boundary_candidate(near, inference)
  expect_true(checked$numerical_checks$BoundaryCandidatePassed)
  expect_lt(abs(checked$numerical_checks$BoundaryCandidateGain), 1e-10)
  expect_identical(checked$covariance, inference$covariance)
})

test_that("prespecified population contrasts use full nuisance-adjusted covariance", {
  x <- population_interval_fixture()
  local_mocked_bindings(fit_mfrm=function(...) stop("No refitting"),
    mfrm_native_location_inference=function(...) stop("No repeated qualification"),
    compute_mml_parameter_covariance=function(...) stop("No repeated information"), .package="mfrmr")
  # Reversed input columns; include a mean and dependent differences.
  L <- rbind(at_X1=c(1,1), X1_minus_Xminus1=c(2,0), reverse=c(-2,0))
  colnames(L) <- c("X","(Intercept)")
  ci <- suppressWarnings(mfrm_population_intervals(x$native, level=.9, contrasts=L))
  base <- suppressWarnings(mfrm_population_intervals(x$native, level=.9))
  ordered <- L[,c("(Intercept)","X")]
  expected <- ordered %*% x$fixture$reference_covariance[11:12,11:12] %*% t(ordered)
  expect_equal(ci$contrast_covariance, expected, tolerance=1e-5, ignore_attr=TRUE)
  expect_identical(ci$contrasts, ordered)
  expect_identical(vcov(ci), vcov(base))
  expect_identical(ci$fit, base$fit)
  tab <- subset(summary(ci), Type=="coefficient_contrast")
  expect_identical(tab$Target, paste0("contrast:",rownames(L)))
  expect_true(all(tab$CIEligible))
  expect_equal(tab$Estimate, drop(ordered %*% x$fixture$fit$population$coefficients), ignore_attr=TRUE)
  expect_equal(tab$SE, sqrt(diag(expected)), tolerance=1e-5, ignore_attr=TRUE)
  expect_equal(tab$Lower, tab$Estimate-qnorm(.95)*tab$SE)
  expect_equal(tab$Upper, tab$Estimate+qnorm(.95)*tab$SE)
  expect_equal(tab$Lower[2], -tab$Upper[3])
  expect_equal(tab$Upper[2], -tab$Lower[3])
  retained <- subset(ci$table, Type!="coefficient_contrast")
  rownames(retained) <- NULL
  expect_identical(retained, base$table)
  # Using only coefficient-wise SEs would lose this cross-covariance.
  expect_gt(abs(expected[1,1]-sum(diag(x$fixture$reference_covariance[11:12,11:12]))),1e-5)
})

test_that("population contrast input is checked before numerical inference", {
  x <- population_interval_fixture()
  local_mocked_bindings(mfrm_native_location_inference=function(...) stop("Inference called"), .package="mfrmr")
  L <- matrix(c(0,2),1,dimnames=list("difference",c("(Intercept)","X")))
  bad <- list(unname(L), as.data.frame(L), L+NA_real_, L*0, L+1i,
    L[FALSE,,drop=FALSE], L[,1,drop=FALSE], rbind(L,L))
  wrong <- L; colnames(wrong) <- c("(Intercept)","Y"); bad <- c(bad,list(wrong))
  for (label in c("", " ", NA_character_)) {
    wrong <- L; rownames(wrong) <- label; bad <- c(bad,list(wrong))
  }
  wrong <- L; colnames(wrong) <- c("X","X"); bad <- c(bad,list(wrong))
  for (item in bad) expect_error(mfrm_population_intervals(x$fixture$fit, contrasts=item),
    "all fitted coefficient column names")
})

test_that("population contrasts retain refusal and finite point estimates", {
  x <- population_interval_fixture(); bad <- x$native
  bad$checks$Passed[nrow(bad$checks)] <- FALSE
  bad$checks$Detail[nrow(bad$checks)] <- "Finer information unavailable"
  bad$parameter_covariance <- NULL
  L <- matrix(c(0,2),1,dimnames=list("difference",c("(Intercept)","X")))
  ci <- suppressWarnings(mfrm_population_intervals(bad, contrasts=L))
  tab <- subset(ci$table, Type=="coefficient_contrast")
  expect_equal(tab$Estimate, 2*x$fixture$fit$population$coefficients[["X"]])
  expect_false(tab$CIEligible)
  expect_identical(tab$Status,"numerical_review_failed")
  expect_match(tab$InferenceReview,"Finer information unavailable")
  expect_true(all(is.na(tab[c("SE","Lower","Upper")])) && all(is.na(ci$contrast_covariance)))
  # Finite input weights can still overflow/underflow the reported covariance.
  for (weight in c(1e308,1e-300)) {
    L[1,2] <- weight
    extreme <- suppressWarnings(mfrm_population_intervals(x$native,contrasts=L))
    row <- subset(extreme$table,Type=="coefficient_contrast")
    expect_false(row$CIEligible)
    expect_identical(row$NaturalScaleStatus,"unavailable")
    expect_true(all(is.na(row[c("SE","Lower","Upper")])))
    expect_true(all(extreme$table$CIEligible[extreme$table$Type!="coefficient_contrast"]))
  }
})

test_that("saved population contrasts flow into reports and full-precision exports", {
  x <- population_interval_fixture()
  L <- matrix(c(1,1.234567890123),1,dimnames=list("Target",c("(Intercept)","X")))
  ci <- suppressWarnings(mfrm_population_intervals(x$native, contrasts=L))
  local_mocked_bindings(fit_mfrm=function(...) stop("No refitting"),
    mfrm_population_intervals=function(...) stop("No repeated intervals"),
    compute_mml_parameter_covariance=function(...) stop("No information on replay"), .package="mfrmr")
  res <- mfrm_results(ci$fit,intervals=list(effects=ci),include="fit",compute="never")
  reopened <- unserialize(serialize(res,NULL))
  report <- mfrm_report(reopened)
  expect_match(report$markdown,"not individual ability differences or causal effects",fixed=TRUE)
  expect_identical(report$tables$population_effects_intervals$Estimate,ci$table$Estimate)
  expect_identical(report$tables$population_effects_interpretation_notes$Note,ci$cautions)
  expect_identical(report$tables$population_effects_contrast_weights[["beta:X"]],L[1,"X"])
  expect_identical(names(report$tables$population_effects_contrast_covariance),c("Target","contrast:Target"))
  expect_identical(reopened$population_intervals$effects$fit$population,ci$fit$population)
  expect_s3_class(apa_table(ci,which="contrast_weights"),"apa_table")
  path <- tempfile("population-contrasts-"); on.exit(unlink(path,recursive=TRUE),add=TRUE)
  export_mfrm_results(reopened,output_dir=path,include="tables",acknowledge_sensitive=TRUE)
  files <- list.files(path,recursive=TRUE,full.names=TRUE)
  weights <- files[grepl("population_effects_contrast_weights[.]csv$",files)]
  estimates <- files[grepl("population_effects_intervals[.]csv$",files)]
  notes <- files[grepl("population_effects_interpretation_notes[.]csv$",files)]
  expect_length(weights,1L); expect_length(estimates,1L)
  expect_length(notes,1L)
  expect_identical(read.csv(notes)$Note,ci$cautions)
  expect_equal(read.csv(weights,check.names=FALSE)[["beta:X"]],L[1,"X"],tolerance=1e-14)
  expect_equal(read.csv(estimates)$Estimate,ci$table$Estimate,tolerance=1e-14)
})

test_that("variance display overflow preserves the finite log-scale interval", {
  x <- population_interval_fixture(); V <- diag(length(x$fixture$fit$opt$par))
  V[nrow(V),ncol(V)] <- 1e6
  local_mocked_bindings(mfrm_native_location_inference=function(...) list(
    covariance=V, checks=data.frame(Check="Injected range check",Passed=TRUE,Detail="Numerical range fixture"),
    check=list(eligible=TRUE,review="Numerical range fixture")), .package="mfrmr")
  ci <- suppressWarnings(mfrm_population_intervals(x$fixture$fit))
  last <- nrow(ci$table)
  expect_true(ci$table$CIEligible[last])
  expect_identical(ci$table$NaturalScaleStatus[last],"outside_floating_point_range")
  expect_equal(ci$table$Lower[last],0)
  expect_equal(ci$table$Upper[last],Inf)
  expect_true(all(is.finite(unlist(ci$table[last,c("WorkingLower","WorkingUpper")]))))
  expect_equal(ci$table$WorkingUpper[last]-ci$table$WorkingLower[last],2*qnorm(.975)*1000)
})

test_that("population results retain targets and unavailable states through reporting", {
  x <- population_interval_fixture()
  local_mocked_bindings(fit_mfrm=function(...) stop("No refitting"),
    compute_mml_parameter_covariance=function(...) stop("No information on replay"), .package="mfrmr")
  ci <- suppressWarnings(mfrm_population_intervals(x$native))
  res <- mfrm_results(ci$fit,intervals=list(effects=ci,raters=x$native),include="fit",compute="never")
  expect_identical(res$tables$population_effects_intervals$Lower,ci$table$Lower)
  expect_identical(mfrm_report(res)$tables$population_effects_intervals$Lower,ci$table$Lower)
  expect_match(mfrm_report(res)$markdown,"Conditional-normal population uncertainty")
  expect_true(mfrm_results_uses_saved_replay(res))
  expect_s3_class(apa_table(ci),"apa_table")
  expect_identical(unserialize(serialize(res,NULL))$population_intervals$effects$table,ci$table)
  bad <- ci$fit; bad$population$sigma2 <- bad$population$sigma2+.1
  expect_error(mfrm_results(bad,intervals=ci,compute="never"),"Saved population")
  expect_error(mfrm_results(ci$fit,intervals=list(a=ci,A=x$native),compute="never"),"distinct")
  bad <- x$native; bad$checks$Passed[nrow(bad$checks)] <- FALSE
  bad$checks$Detail[nrow(bad$checks)] <- "Finer information unavailable"
  bad$parameter_covariance <- NULL
  failed <- suppressWarnings(mfrm_population_intervals(bad))
  refused <- mfrm_results(failed$fit,intervals=failed,include="fit",compute="never")
  expect_true(all(!refused$tables$population_inference_intervals$CIEligible))
})
