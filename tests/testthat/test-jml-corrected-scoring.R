corrected_scoring_fixture <- local({
  saved <- NULL
  function() {
    skip_if_not_installed("nleqslv")
    if (is.null(saved)) {
      x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
      fit <- suppressWarnings(fit_mfrm(x$data,x$person,x$facets,x$score,
        model="GPCM",method="JML",step_facet=x$owner,rating_min=0,rating_max=3,
        category_policy="preserve",jml_correction_order=2L))
      rows <- x$data[1:6,]; rows[[x$person]] <- "001"
      low <- high <- sparse <- missing <- rows
      low[[x$person]] <- "low"; low[[x$score]] <- 0
      high[[x$person]] <- "high"; high[[x$score]] <- 3
      sparse <- sparse[1,]; sparse[[x$person]] <- "sparse"
      missing <- missing[1,]; missing[[x$person]] <- "missing"; missing[[x$score]] <- NA_real_
      rows <- rbind(rows,low,high,sparse,missing)
      rows$Event <- seq_len(nrow(rows))
      cal <- freeze_mfrm_calibration(validate_mfrm_calibration(
        extract_mfrm_calibration(fit,scoring_quad_points=101L)))
      saved <<- list(fit=fit,cal=cal,rows=rows,x=x)
    }
    saved
  }
})

test_that("corrected calibration scores agree with independent continuous posterior integration", {
  z <- corrected_scoring_fixture(); f <- z$fit; x <- z$x
  before <- serialize(f,NULL)
  expect_identical(z$cal$header$schema_version,5L)
  expect_length(mfrmr:::mfrmr_calibration_find_prohibited(z$cal),0)
  expect_equal(nrow(review_mfrm_calibration(z$cal)),0)
  b <- x$results[["2"]]$beta
  judge <- c(b[1:2],-sum(b[1:2])); criterion <- c(b[3],-b[3])
  steps <- rbind(c(b[4:5],-sum(b[4:5])),c(b[6:7],-sum(b[6:7])))
  slopes <- exp(c(b[8],-b[8]))
  for (prior in list(NULL,list(mean=.4,sd=1.3))) {
    native <- predict_mfrm_units(f,z$rows[!is.na(z$rows[[x$score]]),],
      scoring_quad_points=101L,scoring_prior=prior)
    expect_identical(native$fit_identity, prediction_fit_identity(f))
    scored <- score_mfrm_calibration(z$cal,z$rows,event_id="Event",
      missing_response="omit",scoring_prior=prior)
    expect_s3_class(summary(native),"summary.mfrm_unit_prediction")
    expect_equal(build_summary_table_bundle(native)$tables$estimates$CorrectionOrder,
      native$estimates$CorrectionOrder)
    expect_s3_class(summary(scored),"summary.mfrm_calibration_score")
    fields <- c("Estimate","SD","Lower","Upper")
    pos <- match(scored$estimates$Person,native$estimates$Person)
    expect_equal(as.matrix(scored$estimates[fields]),as.matrix(native$estimates[pos,fields]),ignore_attr=TRUE,tolerance=1e-10)
    pr <- prior %||% list(mean=0,sd=1)
    for (id in scored$estimates$Person) {
      d <- z$rows[z$rows[[x$person]]==id,,drop=FALSE]
      # Literal saved eight-coordinate response equation; no package kernel,
      # stored quadrature rule or interval routine supplies the reference.
      density <- function(theta) {
        ll <- rep(0,length(theta))
        for (i in seq_len(nrow(d))) {
          r <- as.integer(d[[x$facets[1]]][i]); c <- as.integer(d[[x$facets[2]]][i])
          lw <- slopes[c]*(outer(theta-judge[r]-criterion[c],0:3)-
            matrix(c(0,cumsum(steps[c,])),length(theta),4,byrow=TRUE))
          mx <- apply(lw,1,max)
          ll <- ll+lw[,d[[x$score]][i]+1]-mx-log(rowSums(exp(lw-mx)))
        }
        exp(ll+dnorm(theta,pr$mean,pr$sd,log=TRUE))
      }
      integrate_value <- function(fn,upper=Inf) integrate(fn,-Inf,upper,
        rel.tol=1e-10,abs.tol=1e-12,subdivisions=500)$value
      normalizer <- integrate_value(density)
      mean <- integrate_value(function(t) t*density(t))/normalizer
      sd <- sqrt(integrate_value(function(t) (t-mean)^2*density(t))/normalizer)
      ends <- vapply(c(.025,.975),function(p) uniroot(function(t)
        integrate_value(density,t)/normalizer-p,c(-15,15),tol=1e-9)$root,numeric(1))
      actual <- as.numeric(scored$estimates[scored$estimates$Person==id,fields])
      expect_equal(actual,c(mean,sd,ends),tolerance=1e-5)
    }
    expect_true(all(scored$estimates$ScoreIntegrationReady))
    expect_true(all(scored$estimates$CorrectionOrder==2L))
    expect_true(all(native$estimates$CalibrationMethod=="Corrected JML"))
    expect_identical(scored$person_dispositions$Disposition[
      scored$person_dispositions$Person=="missing"],"not_scored")
    expect_true(all(scored$person_dispositions$Disposition[
      scored$person_dispositions$Person %in% c("low","high","sparse")]=="scored_review"))
    expect_output(print(scored),"Experimental corrected JML")
    expect_output(print(z$cal),"corrected JML")
    p <- plot(scored,draw=FALSE)
    expect_s3_class(p,"mfrm_plot_data")
    expect_s3_class(ggplot2::ggplot_build(as_ggplot(p)),"ggplot_built")
  }
  pv <- sample_mfrm_plausible_values(f,z$rows[1:6,],scoring_quad_points=101L,n_draws=10,seed=7)
  expect_true(all(pv$values$CalibrationMethod=="Corrected JML"))
  expect_true(all(pv$values$CorrectionOrder==2L))
  expect_s3_class(summary(pv),"summary.mfrm_plausible_values")
  path <- tempfile(); dir.create(path); withr::defer(unlink(path,recursive=TRUE))
  export_summary_appendix(native,output_dir=path,include_html=FALSE)
  csv <- list.files(path,"estimates.*[.]csv$",recursive=TRUE,full.names=TRUE)
  expect_length(csv,1)
  expect_true(all(read.csv(csv)$CalibrationMethod=="Corrected JML"))
  expect_true(all(read.csv(csv)$CorrectionOrder==2L))
  broken <- native; broken$estimates$CorrectionOrder <- 4L
  expect_error(summary(broken),"numerical review records")
  broken <- scored; broken$estimates$CorrectionOrder <- 4L
  expect_error(summary(broken),"numerical-check records")
  expect_identical(serialize(f,NULL),before)
})

test_that("corrected scoring requires its own intact equation and rejects unsupported inputs", {
  z <- corrected_scoring_fixture(); f <- z$fit
  for (field in c("beta","table","data","order","point","weight","owner","anchor","levels")) {
    bad <- f
    if (field=="beta") bad$jml_adjustment$point$beta[1] <- bad$jml_adjustment$point$beta[1]+.1
    if (field=="table") bad$slopes$Estimate[1] <- bad$slopes$Estimate[1]+.1
    if (field=="data") bad$prep$data$score_k[1] <- (bad$prep$data$score_k[1]+1)%%4
    if (field=="order") bad$jml_adjustment$estimator$order <- 0L
    if (field=="point") bad$jml_adjustment$point$available <- FALSE
    if (field=="weight") bad$prep$data$Weight[1] <- 2
    if (field=="owner") bad$config$slope_facet <- bad$config$facet_names[1]
    if (field=="anchor") bad$config$theta_spec <- list(anchors=0)
    if (field=="levels") bad$prep$levels[[bad$config$facet_names[1]]] <- rev(bad$prep$levels[[bad$config$facet_names[1]]])
    expect_error(extract_mfrm_calibration(bad),"SOURCE_READINESS_INELIGIBLE")
  }
  shifted <- f
  shifted$config$rating_min <- shifted$prep$rating_min <- 5
  shifted$config$rating_max <- shifted$prep$rating_max <- 8
  shifted$prep$data$Score <- shifted$prep$data$Score+5
  shifted$prep$score_map$OriginalScore <- shifted$prep$score_map$OriginalScore+5
  shifted$prep$score_map$InternalScore <- shifted$prep$score_map$InternalScore+5
  moved <- freeze_mfrm_calibration(validate_mfrm_calibration(extract_mfrm_calibration(shifted,scoring_quad_points=101)))
  rows <- z$rows; rows[[z$x$score]] <- rows[[z$x$score]]+5
  expect_equal(score_mfrm_calibration(moved,rows,event_id="Event",missing_response="omit")$estimates$Estimate,
    score_mfrm_calibration(z$cal,z$rows,event_id="Event",missing_response="omit")$estimates$Estimate,tolerance=1e-12)
  failed <- f; failed$jml_adjustment$point$available <- FALSE
  expect_error(predict_mfrm_units(failed,z$rows,readiness_policy="review"),"unambiguous")
  no_cov <- f; no_cov$jml_adjustment$covariance$available <- FALSE
  c <- extract_mfrm_calibration(no_cov)
  expect_identical(c$parameters,z$cal$parameters)
  reordered <- f; reordered$steps <- f$steps[nrow(f$steps):1,]
  expect_identical(extract_mfrm_calibration(reordered)$parameters,z$cal$parameters)
  for (field in c("checks","version","order","coordinates")) {
    bad <- z$cal
    if (field=="checks") bad$eligibility$source_scoring_evidence$local_calibration_review$checks[1] <- .01
    if (field=="version") bad$header$schema_version <- 4L
    if (field=="order") bad$eligibility$source_scoring_evidence$local_calibration_review$correction_order <- 0L
    if (field=="coordinates") bad$parameters$coordinates$Value[1] <- 10
    expect_gt(nrow(review_mfrm_calibration(bad)),0)
    expect_error(score_mfrm_calibration(bad,z$rows,event_id="Event",missing_response="omit"),class="mfrm_calibration_error")
  }
  rows <- z$rows[1:6,]; rows$Weight <- 2
  expect_error(predict_mfrm_units(f,rows,weight="Weight"),"unit observation weights")
  expect_error(score_mfrm_calibration(z$cal,rows,weight="Weight",event_id="Event"),"SCORING_WEIGHT_UNSUPPORTED")
  rows <- z$rows[1:6,]; rows[[z$x$facets[1]]] <- "unknown"
  expect_error(score_mfrm_calibration(z$cal,rows,event_id="Event"),"LEVEL")
  expect_error(score_mfrm_calibration(z$cal,z$rows,missing_response="omit"),"EVENT_DUPLICATE")
  empty <- score_mfrm_calibration(z$cal,z$rows[is.na(z$rows[[z$x$score]]),],missing_response="omit")
  expect_equal(nrow(empty$estimates),0)
  expect_s3_class(summary(empty),"summary.mfrm_calibration_score")
  low <- freeze_mfrm_calibration(validate_mfrm_calibration(extract_mfrm_calibration(f,scoring_quad_points=5)))
  expect_error(score_mfrm_calibration(low,z$rows,event_id="Event",missing_response="omit"),"SCORING_INTEGRATION_FAILED")
})

test_that("corrected calibration replays in a fresh session without source fit or estimator", {
  skip_if_not_installed("callr")
  z <- corrected_scoring_fixture()
  path <- tempfile(fileext=".rds"); withr::defer(unlink(path))
  save_mfrm_calibration(z$cal,path)
  expected <- score_mfrm_calibration(z$cal,z$rows,event_id="Event",missing_response="omit")
  worker <- function(root,path,rows) {
    if(file.exists(file.path(root,"R","api-calibration.R"))) pkgload::load_all(root,quiet=TRUE,compile=FALSE) else library("mfrmr",lib.loc=dirname(root),character.only=TRUE)
    testthat::local_mocked_bindings(fit_mfrm=function(...) stop("unexpected refit"),
      mfrm_jml_scoring_components=function(...) stop("unexpected training check"),
      predict_mfrm_units=function(...) stop("unexpected native scoring"),.package="mfrmr")
    s <- mfrmr::score_mfrm_calibration(mfrmr::load_mfrm_calibration(path),rows,event_id="Event",missing_response="omit")
    list(estimates=s$estimates,settings=s$settings,summary=summary(s))
  }
  environment(worker) <- baseenv()
  out <- callr::r(worker,args=list(normalizePath(find.package("mfrmr")),path,z$rows),libpath=.libPaths())
  expect_identical(out$estimates,expected$estimates)
  expect_identical(out$settings,expected$settings)
  expect_identical(out$summary,summary(expected))
})
