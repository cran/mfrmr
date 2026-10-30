gmfrm_interval_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      saved <- readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))
      problem <- mfrm_gmfrm_problem(saved$data, 2L, gauss_hermite_normal(31L))
      result <- mfrm_gmfrm_em(problem, start=saved$parameters, maxit=1L)
      # The retained optimum meets the current score tolerance without an EM step.
      fit <- mfrm_gmfrm_fit_result(problem, result)
      ci <- suppressWarnings(confint(fit, parm = "slopes"))
      cached <<- list(fit=fit, ci=ci, saved=saved)
    }
    cached
  }
})

test_that("two-family fits require an explicit component-slope target", {
  expect_error(confint(gmfrm_interval_fixture()$fit), "two-family GPCM.*parm = 'slopes'")
})

test_that("a small optimization residual gives a retained numerical caution", {
  x <- readRDS(test_path("fixtures","gmfrm-small-optimization-residual.rds"))
  problem <- mfrm_gmfrm_problem(x$data,2L,gauss_hermite_normal(31L))
  fit <- mfrm_gmfrm_fit_result(problem,mfrm_gmfrm_em(problem,start=x$parameters,maxit=1L))
  before <- fit
  expect_true(fit$summary$Converged)
  expect_warning(ci <- confint(fit, parm = "slopes"),"Small optimization residual")
  expect_identical(fit,before)
  check <- attr(ci,"numerical_checks")
  expect_gt(check$CurvatureScaledGradient,1e-4)
  expect_lt(check$CurvatureScaledGradient,.01)
  expect_true(check$OptimizationCaution)
  expect_true(all(attr(ci,"diagnostics")$CIEligible))
  expect_equal(attr(ci,"diagnostics")$Estimate,fit$slopes$Estimate,tolerance=1e-14)
  reference_se <- attr(x$reference_intervals,"diagnostics")$LogSE
  expect_lt(max(abs(log(ci)-log(x$reference_intervals))/reference_se),.001)
  expect_output(print(ci),"Small optimization residual")
  expect_match(paste(plot_data(as_ggplot(ci), "notes")$Text, collapse = "\n"),"Small optimization residual")
  res <- mfrm_results(fit,include="fit",compute="never",intervals=list(slopes=ci))
  expect_match(mfrm_report(res)$markdown,"Small optimization residual")
  expect_identical(res$tables$gpcm_slopes_numerical_checks,check)
  # A tiny mean score can hide a material correction near a weak slope.
  unstable <- mfrm_gmfrm_problem(x$unstable$data,2L,gauss_hermite_normal(31L))
  bad <- mfrm_gmfrm_fit_result(unstable,mfrm_gmfrm_em(unstable,
    start=x$unstable$parameters,maxit=1L))
  expect_true(bad$summary$Converged)
  unavailable <- confint(bad, parm = "slopes")
  expect_true(all(is.na(unavailable)))
  check_bad <- attr(unavailable,"numerical_checks")
  expect_lt(check_bad$MaximumMeanScore,1e-6)
  expect_gt(check_bad$CurvatureScaledGradient,.01)
  expect_match(attr(unavailable,"diagnostics")$InferenceReview[1],"at most 0.01")
})

test_that("two-family intervals use the inverse full observed marginal information", {
  x <- gmfrm_interval_fixture(); ci <- x$ci; fit <- x$fit
  expect_true(fit$summary$Converged)
  tab <- attr(ci,"diagnostics")
  expect_true(all(tab$CIEligible))
  # Construct the target map independently: first family sums to zero, second is free.
  jac <- matrix(0,6,13)
  jac[1:3,9:10] <- rbind(c(1,0),c(0,1),c(-1,-1))
  jac[4:6,11:13] <- diag(3)
  covariance <- jac %*% solve(x$saved$information) %*% t(jac)
  expect_equal(unname(attr(ci,"covariance")),unname(covariance),tolerance=1e-8)
  expect_gt(max(abs(covariance[1:3,4:6])),1e-5)
  expect_equal(colSums(covariance[1:3,]),rep(0,6),tolerance=1e-12)
  se <- sqrt(diag(covariance)); estimate <- as.numeric(jac %*% fit$opt$par)
  expect_equal(tab$LogSE,se,tolerance=1e-8)
  expect_equal(as.numeric(ci[,1]),exp(estimate-qnorm(.975)*se),tolerance=1e-8)
  expect_equal(as.numeric(ci[,2]),exp(estimate+qnorm(.975)*se),tolerance=1e-8)
  expect_true(all(is.na(tab$PValue) & is.na(tab$NullValue)))
  expect_identical(tab$SlopeOwner,rep(c("Task","Rater"),each=3))
  expect_identical(tab$ScaleReference,rep(c("geometric_mean_one","fixed_standard_normal"),each=3))
  expect_identical(attr(ci,"settings")$scale,"standardized")
  expect_true(all(attr(ci,"checks")$Passed))
  checks <- attr(ci,"numerical_checks")
  expect_identical(checks$ScoreRank,13L)
  expect_equal(checks$ComparisonPoints,61L)
  expect_lt(checks$QuadratureScoreShift,.01)
  expect_lt(checks$QuadratureCovarianceChange,.01)
  expect_warning(confint(fit, parm = "slopes",level=.9,simultaneous="bonferroni"),"Experimental two-family")
  expect_false(fit$summary$InferenceReady)
  expect_false(fit$summary$ICEligible)
  expect_true(all(!fit$slopes$CIEligible))
  for (args in list(list(scale="relative"),list(method="sandwich"),
      list(contrasts=matrix(1)),list(adjust=TRUE))) {
    expect_error(do.call(confint,c(list(fit, parm = "slopes"),args)),"Two-family intervals currently require")
  }
  expect_error(confint(fit, parm = "slopes",method="profile"),"identify one owner and level")
})

test_that("failed checks retain estimates, missing bounds and an actionable reason", {
  fit <- gmfrm_interval_fixture()$fit
  fail <- function(x, reason) {
    ci <- confint(x, parm = "slopes")
    expect_true(all(is.na(ci)))
    expect_equal(attr(ci,"diagnostics")$Estimate,fit$slopes$Estimate)
    expect_true(any(!attr(ci,"checks")$Passed))
    expect_match(attr(ci,"diagnostics")$InferenceReview[1],reason)
  }
  x <- fit; x$summary$Converged <- FALSE; fail(x,"per-Person")
  x <- fit; x$opt$value <- x$opt$value+1; fail(x,"likelihood must match")
  x <- fit; x$prep$data$Score[1] <- (x$prep$data$Score[1]+1) %% 3
  fail(x,"original fitted model")
  x <- fit; x$config$facet_specs[[1]]$center <- !isTRUE(x$config$facet_specs[[1]]$center)
  fail(x,"original fitted model")
  withr::local_options(mfrmr.max_information_bytes=1)
  fail(fit,"unregularized inverse")
})

test_that("a converged public fit can still fail the integration check", {
  d <- gmfrm_interval_fixture()$saved$data
  args <- list(data=d,person="Person",facets=c("Task","Rater"),score="Score",
    model="GPCM",method="MML",slope_facet=c("Task","Rater"),step_facet="Rater",
    noncenter_facet="Rater",gpcm_mml_identification="fixed_standard_normal",
    mml_engine="em",maxit=500L)
  coarse <- do.call(fit_mfrm,c(args,list(quad_points=5L)))
  expect_true(coarse$summary$Converged)
  ci <- confint(coarse, parm = "slopes")
  expect_true(all(is.na(ci)))
  checks <- attr(ci,"checks")
  expect_identical(checks$Check[!checks$Passed],"Quadrature sensitivity")
  expect_gt(attr(ci,"numerical_checks")$QuadratureScoreShift,1)
  fine <- do.call(fit_mfrm,c(args,list(quad_points=31L)))
  expect_warning(ci_fine <- confint(fine, parm = "slopes"),"Experimental two-family")
  expect_true(all(attr(ci_fine,"diagnostics")$CIEligible))
  expect_equal(as.numeric(ci_fine),as.numeric(gmfrm_interval_fixture()$ci),tolerance=1e-5)
  loose <- do.call(fit_mfrm,c(args,list(quad_points=31L,em_score_tol=1e-3)))
  expect_true(loose$summary$Converged)
  ci_loose <- confint(loose, parm = "slopes")
  expect_true(all(is.na(ci_loose)))
  checks <- attr(ci_loose,"checks")
  expect_identical(checks$Check[!checks$Passed],"Stationary local maximum")
  expect_gt(attr(ci_loose,"numerical_checks")$MaximumMeanScore,1e-6)
  expect_true(is.na(attr(ci_loose,"numerical_checks")$QuadratureScoreShift))
})

test_that("interval identities survive renamed facets and overlapping level labels", {
  x <- gmfrm_interval_fixture(); d <- x$saved$data
  levels(d$Task) <- c("00:shared","t2","t3")
  levels(d$Rater) <- c("00:shared","r2","r3")
  names(d) <- c("Candidate","観点 名","Judge-ID","Rating")
  problem <- mfrm_gmfrm_problem(d,2L,gauss_hermite_normal(31L),
    slope_facets=c("観点 名","Judge-ID"),person="Candidate",score="Rating")
  fit <- mfrm_gmfrm_fit_result(problem,mfrm_gmfrm_em(problem,start=x$saved$parameters,maxit=1L))
  expect_warning(ci <- confint(fit, parm = "slopes"),"Experimental two-family")
  tab <- attr(ci,"diagnostics")
  expect_true(all(tab$CIEligible))
  expect_identical(tab$SlopeOwner,rep(c("観点 名","Judge-ID"),each=3))
  expect_identical(tab$SlopeLevel[c(1,4)],rep("00:shared",2))
  expect_false(anyDuplicated(rownames(ci)) > 0)
  expect_equal(as.numeric(ci),as.numeric(x$ci),tolerance=1e-8)
})

test_that("saved two-family intervals keep owner identity and checks without recomputation", {
  x <- gmfrm_interval_fixture(); ci <- x$ci
  local_mocked_bindings(fit_mfrm=function(...) stop("unexpected fit"),
    compute_mml_parameter_covariance=function(...) stop("unexpected covariance"),.package="mfrmr")
  expect_output(print(ci),"Experimental Joint-information")
  expect_output(print(ci),"Task: geometric mean one; Rater: free slopes",fixed=TRUE)
  expect_match(paste(plot_data(as_ggplot(ci), "notes")$Text, collapse = "\n"),"Experimental")
  expect_match(paste(plot_data(as_ggplot(ci), "notes")$Text, collapse = "\n"),"different reference meanings",fixed=TRUE)
  expect_null(as_ggplot(ci,caption=NULL)$labels$caption)
  expect_identical(plot_data(ci,"table")$SlopeOwner,attr(ci,"diagnostics")$SlopeOwner)
  expect_identical(plot_data(as_ggplot(ci),"checks"),attr(ci,"checks"))
  expect_identical(plot_data(ci,"numerical_checks"),attr(ci,"numerical_checks"))
  expect_identical(apa_table(ci,which="checks")$table$Passed,attr(ci,"checks")$Passed)
  expect_identical(apa_table(ci)$table$SlopeLevel,attr(ci,"diagnostics")$SlopeLevel)
  res <- mfrm_results(x$fit,include=c("fit","plots"),compute="never",intervals=list(slopes=ci))
  expect_identical(res$tables$gpcm_slopes_checks,attr(ci,"checks"))
  report <- mfrm_report(res)
  expect_match(report$markdown,"6 experimental slope intervals are available",fixed=TRUE)
  expect_match(report$markdown,"Curve intervals remain unavailable",fixed=TRUE)
  expect_identical(report$tables$gpcm_slopes_numerical_checks,attr(ci,"numerical_checks"))
  folder <- tempfile("gmfrm-ci-"); withr::defer(unlink(folder,recursive=TRUE))
  exported <- export_mfrm_results(res,output_dir=folder,preset="starter",
    acknowledge_sensitive=TRUE,plot_width=1200,plot_height=900)
  expect_equal(nrow(exported$plot_errors),0L)
  expect_true(file.exists(file.path(folder,"mfrmr_results_plot_gpcm_slopes.png")))
  paths <- list.files(folder,recursive=TRUE,full.names=TRUE)
  saved <- readRDS(paths[grepl("\\.rds$",paths)][1])
  expect_identical(saved$gpcm_inference$slopes,ci)
  expect_true(any(grepl("gpcm_slopes_checks.*csv$",paths)))
  expect_identical(plot_data(as_ggplot(saved,type="gpcm_slopes"), "subtitle"),plot_data(as_ggplot(ci), "subtitle"))
})
