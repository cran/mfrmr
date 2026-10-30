jml_public_args <- function() {
  x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  list(data=x$data,person=x$person,facets=x$facets,score=x$score,
    model="GPCM",method="JML",step_facet=x$owner,
    rating_min=0,rating_max=x$rating_max,category_policy="preserve",
    jml_correction_order=2L)
}

test_that("corrected JML retains its estimator through public output and saved export", {
  skip_if_not_installed("nleqslv")
  args <- jml_public_args()
  # A shifted score origin must leave the corrected equations unchanged.
  args$data[[args$score]] <- args$data[[args$score]]+1
  args$rating_min <- 1; args$rating_max <- 4
  expect_warning(fit <- do.call(fit_mfrm,args),"duplicated Person")
  reference <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  expect_s3_class(fit,"mfrm_fit")
  expect_equal(fit$jml_adjustment$point$beta,reference$results[["2"]]$beta,tolerance=1e-7)
  expect_equal(fit$jml_adjustment$covariance$result$vcov,
    reference$results[["2"]]$covariance$vcov,tolerance=1e-8)
  expect_true(fit$summary$Converged)
  expect_false(fit$summary$InferenceReady || fit$summary$ICEligible)
  expect_true(all(is.na(fit$summary[c("LogLik","AIC","BIC")])) )
  expect_null(fit$opt$par)
  s <- summary(fit,include_person=TRUE)
  expect_s3_class(s,"summary.mfrm_fit")
  expect_output(print(fit),"Corrected JML")
  expect_match(s$estimation_note,"residual bias",ignore.case=TRUE)
  expect_match(s$tables$uncertainty$RepeatedSampling,"New independent Persons.*counts fixed")
  expect_match(s$tables$uncertainty$AbilityDistribution,"allowed to differ")
  expect_match(s$tables$uncertainty$SamePersonReassessment,"Not estimated")
  expect_output(print(s),"RootSE sampling: New independent Persons")
  expect_true(all(is.na(s$tables$persons$SE)))
  expect_true(all(!s$tables$slopes$CIEligible))
  estimates <- as.data.frame(fit)
  expect_true(all(estimates$CorrectionOrder==2L))
  expect_true(all(is.na(estimates$RootSE[estimates$Facet=="Person"])))
  expect_match(attr(estimates,"estimation_note"),"RootSE")
  expect_error(summary(fit,include_person=NA),"TRUE or FALSE")
  expect_error(summary(fit,digits=NA_real_),"digits")
  res <- mfrm_results(fit)
  expect_s3_class(res,"mfrm_results")
  expect_identical(res$tables$slopes,fit$slopes)
  expect_true(all(res$plot_map$Available))
  expect_output(print(summary(res)),"slopes")
  # Keep the live corrected fit/covariance and scope refusals on CRAN.
  # Detailed graphics, HTML and replay combinations run in the full suite.
  if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    for(type in c("slopes","locations","steps")) {
      p <- plot(res,type=type,draw=FALSE)
      expect_s3_class(p,"mfrm_plot_data")
      if(requireNamespace("ggplot2",quietly=TRUE)) {
        g <- as_ggplot(p)
        expect_s3_class(g,"ggplot")
        expect_s3_class(ggplot2::ggplot_build(g),"ggplot_built")
      }
    }
    distribution <- plot(fit,style="distribution",draw=FALSE)
    expect_s3_class(distribution,"mfrm_plot_data")
    expect_false(grepl("prior",distribution$data$caption,ignore.case=TRUE))
    expect_gt(distribution$data$display$limits[1],0)
    report <- mfrm_report(res)
    expect_s3_class(report,"mfrm_report")
    expect_match(report$markdown,"No confidence intervals")
    expect_match(report$markdown,"New independent Persons")
    expect_identical(report$tables$uncertainty,s$tables$uncertainty)
    expect_match(mfrm_report(res,output="html")$html,"RootSE")
    path <- tempfile("adjusted-export-")
    on.exit(unlink(path,recursive=TRUE),add=TRUE)
    exported <- export_mfrm_results(res,output_dir=path,preset="starter",acknowledge_sensitive=TRUE)
    expect_equal(nrow(exported$plot_errors),0L)
    index <- paste(readLines(file.path(path,"index.html")),collapse="\n")
    expect_match(index,"Corrected|adjusted")
    expect_false(grepl("required Wright|provisional fitted curves",index))
    saved <- exported$written_files$Path[exported$written_files$Component=="results_rds"]
    restored <- readRDS(saved)
    expect_identical(restored$fit$jml_adjustment,fit$jml_adjustment)
    expect_identical(restored$tables$uncertainty,s$tables$uncertainty)
    replay <- exported$written_files$Path[exported$written_files$Component=="replay_code"]
    withr::with_dir(path,sys.source(replay,envir=new.env(parent=globalenv())))
  } else {
    expect_s3_class(plot(res,type="slopes",draw=FALSE),"mfrm_plot_data")
  }
  # Ordinary inference/diagnostics must never be attached to these roots.
  expect_error(confint(fit),"corrected JML")
  expect_error(diagnose_mfrm(fit),"corrected JML")
  expect_error(fit_measures_table(fit,diagnostics=list()),"corrected JML")
  expect_error(compare_mfrm(fit,fit),"corrected JML")
  expect_error(compute_information(fit),"corrected JML")
  expect_error(mfrm_curve_intervals(fit,fit$prep$data),"corrected JML")
  expect_error(apa_table(fit),"corrected JML")
  expect_error(mfrm_results(fit,diagnostics=list()),"saved estimates only")
  expect_error(mfrm_report(res,style="rater"),"corrected JML")
  expect_error(plot(fit,show_ci=TRUE,draw=FALSE),"interval")
  expect_error(plot(fit,type="wright",draw=FALSE),"arg")
  expect_error(plot(mfrm_results(fit,include="fit"),draw=FALSE),"plots")
})

test_that("corrected JML refuses unsupported settings without changing defaults", {
  args <- jml_public_args()
  bad <- function(name,value) { a <- args; a[name] <- list(value); do.call(fit_mfrm,a) }
  expect_error(bad("model","PCM"),"requires model")
  expect_error(bad("method","MML"),"requires model")
  expect_error(bad("slope_facet",args$facets[1]),"same slope_facet")
  for(option in c("anchors","weight","reltol","optimizer","attach_diagnostics")) {
    value <- switch(option,reltol=1e-7,attach_diagnostics=TRUE,optimizer="BFGS","unsupported")
    expect_error(bad(option,value),"not supported by corrected JML")
  }
  expect_error(bad("rating_min",NULL),"Declare rating_min")
  for(change in c("missing","empty","space")) {
    d <- args$data
    d[[args$person]][1] <- switch(change,missing=NA_character_,empty="",space=" p1")
    expect_error(bad("data",d),"observed Person|IDs must")
  }
  raw <- args; raw$jml_correction_order <- NULL; raw$jml_correction_sampling <- "random_rosters"
  expect_error(do.call(fit_mfrm,raw),"requires an explicit")
  expect_null(formals(fit_mfrm)$jml_correction_order)
  expect_identical(formals(fit_mfrm)$method,quote(c("MML", "JML", "JMLE", "CML")))
})

test_that("unresolved corrected roots remain reportable without optimizer placeholders", {
  # Extended numerical/reference case; shorter live routes remain on CRAN.
  skip_on_cran()
  skip_if_not_installed("nleqslv")
  args <- jml_public_args(); args$maxit <- 1L
  expect_warning(fit <- do.call(fit_mfrm,args),"duplicated Person")
  expect_false(fit$jml_adjustment$point$available)
  expect_true(all(is.na(fit$slopes$Estimate)))
  expect_true(all(is.na(fit$facets$person$Estimate)))
  expect_true(length(fit$jml_adjustment$attempts)>=2L)
  expect_false(fit$summary$Converged)
  expect_false(any(mfrm_results(fit)$plot_map$Available))
  expect_match(mfrm_report(mfrm_results(fit))$markdown,"No supported local solution")
  expect_error(plot(fit,draw=FALSE),"No unambiguous")
})

test_that("covariance failure cannot erase corrected public point estimates", {
  # Extended numerical/reference case; shorter live routes remain on CRAN.
  skip_on_cran()
  skip_if_not_installed("nleqslv")
  local_mocked_bindings(mfrm_jml_adjustment_covariance=function(...)
    stop("Controlled covariance failure"),.package="mfrmr")
  args <- jml_public_args(); args$method <- "JMLE"
  expect_warning(fit <- do.call(fit_mfrm,args),"duplicated Person")
  expect_true(fit$jml_adjustment$point$available)
  expect_false(fit$jml_adjustment$covariance$available)
  expect_true(all(is.finite(fit$slopes$Estimate)))
  expect_true(all(is.na(fit$slopes$RootSE)))
  expect_match(fit$slopes$CovarianceReason,"Controlled covariance failure")
  res <- mfrm_results(fit,include=c("fit","plots","precision"))
  expect_true(all(res$plot_map$Available))
  expect_identical(res$status$Status[res$status$Section=="precision"],"not_available")
  expect_match(mfrm_report(res)$markdown,"Controlled covariance failure")
  expect_s3_class(plot(res,draw=FALSE),"mfrm_plot_data")
})
