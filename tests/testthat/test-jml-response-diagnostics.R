jml_response_fixture <- local({
  fit <- NULL
  function() {
    skip_if_not_installed("nleqslv")
    if (is.null(fit)) {
      x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
      expect_warning(fit <<- fit_mfrm(x$data,x$person,x$facets,x$score,
        model="GPCM",method="JML",step_facet=x$owner,rating_min=0,rating_max=3,
        category_policy="preserve",jml_correction_order=2L),"duplicated Person")
    }
    fit
  }
})

test_that("corrected conditional probabilities agree with the literal response equation", {
  f <- jml_response_fixture()
  x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  b <- x$results[["2"]]$beta
  judge <- c(b[1:2],-sum(b[1:2])); criterion <- c(b[3],-b[3])
  steps <- rbind(c(b[4:5],-sum(b[4:5])),c(b[6:7],-sum(b[6:7])))
  slopes <- exp(c(b[8],-b[8]))
  theta <- f$facets$person$Estimate[match(x$data[[x$person]],f$facets$person$Person)]
  r <- as.integer(x$data[[x$facets[1]]]); c <- as.integer(x$data[[x$facets[2]]])
  expected <- t(vapply(seq_len(nrow(x$data)),function(i) {
    if (!is.finite(theta[i])) return(if(theta[i]<0) c(1,0,0,0) else c(0,0,0,1))
    z <- slopes[c[i]]*((theta[i]-judge[r[i]]-criterion[c[i]])*(0:3)-c(0,cumsum(steps[c[i],])))
    w <- exp(z-max(z)); w/sum(w)
  },numeric(4)))
  d <- mfrm_response_diagnostics(f,group_by=x$facets)
  expect_equal(unname(d$probabilities),expected,tolerance=1e-8)
  expect_equal(rowSums(d$probabilities),rep(1,nrow(x$data)),ignore_attr=TRUE,tolerance=1e-14)
  expect_equal(d$rows$ExpectedScore,as.vector(expected %*% (0:3)),tolerance=1e-8)
  v <- rowSums(expected*(matrix(0:3,nrow(expected),4,byrow=TRUE)-as.vector(expected %*% (0:3)))^2)
  expect_equal(d$rows$PredictiveVariance,v,tolerance=1e-8)
  # Profiles were estimated from all the Person's fitted rows, including repeats.
  profile_score <- rowsum(slopes[c]*d$rows$Residual,x$data[[x$person]])
  expect_lt(max(abs(profile_score)),1e-7)
  expect_identical(names(d$source_data),names(x$data))
  expect_equal(nrow(d$rows),nrow(x$data))
  expect_identical(d$settings$probability_method,"corrected_jml_plugin")
  expect_false(d$settings$person_uncertainty || d$settings$calibration_uncertainty)
  expect_output(print(d),"Corrected JML conditional")
  expect_match(summary(d)$settings$integration,"No integration")
  selected <- c(19L,2L,1L)
  one <- mfrm_response_diagnostics(f,rows=selected)
  expect_equal(one$probabilities,d$probabilities[selected,,drop=FALSE])
  expect_identical(one$rows$InputRow,selected)
  expect_equal(mfrm_response_diagnostics(f,rows=1)$probabilities,d$probabilities[1,,drop=FALSE])
  shifted <- f
  shifted$config$rating_min <- 5; shifted$config$rating_max <- 8
  shifted$prep$data$Score <- shifted$prep$data$Score+5
  moved <- mfrm_response_diagnostics(shifted,group_by=x$facets)
  expect_equal(unname(moved$probabilities),unname(d$probabilities))
  expect_equal(moved$rows$ExpectedScore,d$rows$ExpectedScore+5)
  expect_equal(moved$rows$PredictiveVariance,d$rows$PredictiveVariance)
  expect_equal(moved$rows$Residual,d$rows$Residual)
  reordered <- f; reordered$steps <- f$steps[nrow(f$steps):1,]
  reordered$slopes <- f$slopes[nrow(f$slopes):1,]
  expect_equal(mfrm_response_diagnostics(reordered,group_by=x$facets)$probabilities,d$probabilities)
  no_cov <- f; no_cov$jml_adjustment$covariance$available <- FALSE
  expect_equal(mfrm_response_diagnostics(no_cov,group_by=x$facets)$probabilities,d$probabilities)
  expect_error(mfrm_response_diagnostics(f,rows=c(1,1)),"distinct")
  expect_error(mfrm_response_diagnostics(f,rows=Inf),"distinct")
  expect_error(mfrm_response_diagnostics(f,group_by=x$score),"identifier")
  expect_error(mfrm_response_diagnostics(f,quad_points=31),"omit quad_points")
  failed <- f; failed$jml_adjustment$point$available <- FALSE
  expect_error(mfrm_response_diagnostics(failed),"unambiguous")
})

test_that("extreme profiles retain probabilities and defined Infit without inventing Outfit", {
  f <- jml_response_fixture(); facet <- f$config$facet_names[1]
  d <- mfrm_response_diagnostics(f,group_by=facet)
  zero <- d$rows$Status=="zero_variance"
  expect_true(any(zero) && any(!zero))
  expect_true(all(d$rows$ProbabilityAvailable[zero]))
  expect_true(all(d$rows$PredictiveVariance[zero]==0))
  expect_true(all(d$rows$Residual[zero]==0))
  expect_true(all(is.na(d$rows$StandardizedResidual[zero])))
  expect_match(d$rows$Reason[zero][1],"limiting point mass")
  expect_true(all(d$measures$Status=="partially_available"))
  expect_true(all(is.finite(d$measures$Infit)))
  expect_true(all(is.na(d$measures$Outfit)))
  for (id in d$measures$Level) {
    j <- d$source_data[[facet]]==id
    expect_equal(d$measures$Infit[d$measures$Level==id],
      sum(d$rows$Residual[j]^2)/sum(d$rows$PredictiveVariance[j]))
  }
  finite <- mfrm_response_diagnostics(f,rows=which(!zero),group_by=facet)
  expect_true(all(is.finite(finite$measures$Outfit)))
  for(id in finite$measures$Level) {
    j <- finite$source_data[[facet]][finite$rows$InputRow]==id
    expect_equal(finite$measures$Outfit[finite$measures$Level==id],
      mean(finite$rows$Residual[j]^2/finite$rows$PredictiveVariance[j]))
  }
  extreme_only <- mfrm_response_diagnostics(f,rows=which(zero),group_by=facet)
  expect_true(all(is.na(extreme_only$measures$Infit)))
  expect_true(all(is.na(extreme_only$measures$Outfit)))
  extreme_results <- mfrm_results(f, response_diagnostics = extreme_only)
  overview <- extreme_results$tables$response_overview
  expect_equal(overview$ZeroVariance, sum(zero))
  expect_equal(overview$Unresolved, 0L)
  expect_equal(overview$Available, 0L)
  expect_identical(overview$Status, "review")
  expect_identical(extreme_results$status$Status[extreme_results$status$Section == "response_diagnostics"], "review")
  old <- extreme_results; old$tables$response_overview <- NULL
  old$status$Status[old$status$Section == "response_diagnostics"] <- "available"
  expect_identical(summary(old)$status$Status[summary(old)$status$Section == "response_diagnostics"], "review")
  expect_identical(mfrm_report(old)$tables$response_overview, overview)
  expect_false(any(c("Flag","ZSTD","p_value") %in% names(d$measures)))
  for(style in c("paired","scatter")) {
    p <- plot(d,style=style,draw=FALSE)
    expect_match(p$data$caption,"conditional probabilities")
    expect_match(p$data$caption,"Zero-variance ratings prevent Outfit")
    expect_false(grepl("posterior",p$data$caption))
    expect_equal(p$data$table,d$measures)
    grDevices::pdf(tempfile(fileext=".pdf"))
    expect_silent(plot(d,style=style)); grDevices::dev.off()
    if(requireNamespace("ggplot2",quietly=TRUE)) {
      g <- as_ggplot(p)
      expect_silent(ggplot2::ggplotGrob(g))
      if(style=="paired") expect_equal(nrow(ggplot2::ggplot_build(g)$data[[2]]),nrow(d$measures))
      hidden <- as_ggplot(plot(d,style=style,show_title=FALSE,show_notes=FALSE,draw=FALSE,palette="mono"))
      expect_null(hidden$labels$title); expect_null(hidden$labels$caption)
    }
  }
})

test_that("saved corrected response diagnostics retain source identity through output", {
  f <- jml_response_fixture(); d <- mfrm_response_diagnostics(f,group_by=f$config$facet_names)
  changed <- list(f,f,f,f)
  changed[[1]]$slopes$Estimate[1] <- 1.01*changed[[1]]$slopes$Estimate[1]
  changed[[2]]$facets$person$Estimate[1] <- changed[[2]]$facets$person$Estimate[1]+.1
  changed[[3]]$prep$data <- changed[[3]]$prep$data[nrow(f$prep$data):1,]
  changed[[4]]$prep$data$Score[1] <- 1
  for(bad in changed) expect_error(mfrm_results(bad,response_diagnostics=d),"exact source roster")
  wrong <- d; wrong$settings$probability_method <- "posterior"
  expect_error(mfrm_results(f,response_diagnostics=wrong),"matching calibration")
  # Collecting or displaying a saved object must not fit or calculate probabilities.
  local_mocked_bindings(mfrm_jml_response_diagnostics=function(...) stop("must not recompute"),
    category_prob_gpcm=function(...) stop("must not recompute"),.package="mfrmr")
  res <- mfrm_results(f,response_diagnostics=d,compute="never")
  expect_identical(res$tables$response_measures,d$measures)
  expect_identical(res$tables$response_residuals,d$rows)
  expect_identical(res$status$Status[res$status$Section == "response_diagnostics"], "review")
  expect_s3_class(plot(res,type="response_diagnostics",draw=FALSE),"mfrm_plot_data")
  report <- mfrm_report(res)
  expect_match(report$markdown,"Conditional fitted probabilities")
  expect_false(grepl("corrected_jml_plugin",report$markdown))
  expect_match(report$markdown,"Outfit is undefined")
  path <- tempfile("jml-response-export-"); on.exit(unlink(path,recursive=TRUE),add=TRUE)
  exported <- export_mfrm_results(res,output_dir=path,preset="starter",acknowledge_sensitive=TRUE)
  expect_equal(nrow(exported$plot_errors),0L)
  expect_true("plot_response_diagnostics" %in% exported$written_files$Component)
  restored <- readRDS(exported$written_files$Path[exported$written_files$Component=="results_rds"])
  expect_identical(restored$response_diagnostics,d)
  replay <- exported$written_files$Path[exported$written_files$Component=="replay_code"]
  withr::with_dir(path,sys.source(replay,envir=new.env(parent=globalenv())))
  expect_error(plot(mfrm_results(f),type="response_diagnostics",draw=FALSE),"Supply saved")
  expect_error(diagnose_mfrm(f),"corrected JML")
  expect_error(confint(f),"corrected JML")
  expect_error(mfrm_report(res,style="rater"),"corrected JML")
})
