decision_fixture <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      d <- load_mfrmr_data("example_core")
      d <- d[d$Person %in% unique(d$Person)[1:18],]
      fit <- fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",quad_points=61,maxit=300)
      pcm <- readRDS(test_path("fixtures","ordinary-score-integration.rds"))
      jml <- pcm$fits[["PCM-JML"]]
      calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(
        extract_mfrm_calibration(jml,scoring_quad_points=121)))
      saved <<- list(rsm=fit,data=d,pcm=pcm,calibration=calibration)
    }
    saved
  }
})

test_that("cut inputs reject ambiguous boundaries and unsupported routes", {
  for (cuts in list(numeric(),c(1,0),c(0,0),NA_real_,Inf,1i,"0",matrix(0),
      c(a=0,a=1),setNames(0,""),setNames(0,NA_character_)))
    expect_error(mfrm_validate_cut_scores(cuts),"cut_scores")
  cuts <- c(Basic=-.5,Advanced=.8)
  expect_identical(mfrm_validate_cut_scores(cuts),cuts)
  fit <- readRDS(test_path("fixtures","mfrm-conditional-scoring-gpcm.rds"))$fit
  expect_error(predict_mfrm_units(fit,data.frame(),cut_scores=cuts),"ordinary RSM/PCM")
  expect_error(score_mfrm_persons(structure(list(),class="mfrm_testlet"),cut_scores=cuts),"ordinary RSM/PCM")
})

test_that("continuous cut probabilities match an independent ordinal posterior", {
  cuts <- c(-40,-1,.25,4,4.000001,40)
  for (model in c("RSM","PCM")) for (n in c(4L,40L)) {
    owner <- rep(1:2,length.out=n)
    idx <- list(person=rep(1L,n),score_k=if(n==4L) 0:3 else rep(3L,n),
      weight=rep(1,n),step_idx=owner)
    config <- list(model=model,n_cat=4L,n_person=1L)
    params <- list(steps=c(-.8,0,.8),steps_mat=rbind(c(-1.2,.3,.9),c(-.4,-.1,.5)))
    base <- rep(c(-.8,.8),length.out=n)
    pop <- list(active=TRUE,design_matrix=matrix(1),coefficients=.6,sigma2=.81,person_lookup=1L)
    actual <- mfrmr_person_posterior_intervals(idx,config,params,base,"P",pop,cut_scores=cuts)
    baseline <- mfrmr_person_posterior_intervals(idx,config,params,base,"P",pop)
    expect_identical(actual[,1:2,drop=FALSE],baseline)
    steps <- if(model=="RSM") matrix(params$steps,n,3,byrow=TRUE) else params$steps_mat[owner,]
    cumulative <- t(apply(steps,1,function(x) c(0,cumsum(x))))
    logd <- function(theta) vapply(theta,function(value) {
      logits <- outer(value+base,0:3)-cumulative
      high <- apply(logits,1,max)
      sum(logits[cbind(seq_len(n),idx$score_k+1L)]-high-log(rowSums(exp(logits-high))))+
        dnorm(value,.6,.9,log=TRUE)
    },0)
    peak <- optimize(logd,c(-12,12),maximum=TRUE)$objective
    density <- function(x) exp(logd(x)-peak)
    integrate_ref <- function(a,b) integrate(density,a,b,rel.tol=1e-11,abs.tol=0)$value
    mass <- integrate_ref(-12,12)
    expected <- vapply(cuts,function(cut) if(cut<=-12) 1 else if(cut>=12) 0 else integrate_ref(cut,12)/mass,0)
    expect_equal(unname(actual[1,paste0("Above",seq_along(cuts))]),expected,tolerance=1e-9)
    expect_equal(unname(actual[1,"Band5"]),integrate_ref(4,4.000001)/mass,tolerance=1e-8)
    expect_gt(actual[1,"Band5"],0)
    expect_equal(sum(actual[1,grep("Band",colnames(actual))]),1,tolerance=1e-10)
    shifted <- pop; shifted$coefficients <- pop$coefficients+7
    translated <- mfrmr_person_posterior_intervals(idx,config,params,base-7,"P",shifted,cut_scores=cuts+7)
    expect_equal(translated[,-(1:2),drop=FALSE],actual[,-(1:2),drop=FALSE],tolerance=1e-9)
  }
})

test_that("known-model decision probabilities obey the prior predictive identity", {
  # Exact enumeration avoids adding a simulation study for a posterior integral.
  patterns <- expand.grid(a=0:1,b=0:1); cuts <- c(-.7,.6)
  weighted <- numeric(5); total <- 0; brier <- risk <- numeric(2)
  for (i in seq_len(nrow(patterns))) {
    y <- as.numeric(patterns[i,]); base <- c(-.4,.8)
    density <- function(theta) vapply(theta,function(t) {
      p <- plogis(t+base)
      prod(ifelse(y==1,p,1-p))*dnorm(t)
    },0)
    mass <- integrate(density,-Inf,Inf,rel.tol=1e-11,abs.tol=0)$value
    result <- mfrmr_person_posterior_intervals(list(person=c(1L,1L),score_k=y,step_idx=c(1L,1L)),
      list(model="RSM",n_cat=2L,n_person=1L),list(steps=0),base,"P",cut_scores=cuts)
    weighted <- weighted+mass*result[1,-(1:2)]; total <- total+mass
    p <- result[1,paste0("Above",1:2)]
    joint <- vapply(cuts,function(cut) integrate(density,cut,Inf,rel.tol=1e-11,abs.tol=0)$value,0)
    brier <- brier+joint*(1-p)^2+(mass-joint)*p^2
    risk <- risk+mass*p*(1-p)
  }
  expect_equal(total,1,tolerance=1e-10)
  expect_equal(unname(weighted),c(pnorm(cuts,lower.tail=FALSE),diff(pnorm(c(-Inf,cuts,Inf)))),tolerance=1e-9)
  expect_equal(brier,risk,tolerance=1e-10)
})

test_that("extreme cuts and narrow bands retain normal prior probabilities", {
  cuts <- c(-1e300,-1,0,1,5,5+1e-6,1e300)
  actual <- mfrm_cut_probabilities(function(a,b) integrate(dnorm,a,b,rel.tol=1e-10,abs.tol=0)$value,cuts,1)
  expect_equal(unname(actual[seq_along(cuts)]),pnorm(cuts,lower.tail=FALSE),tolerance=1e-12)
  expect_gt(actual["Band6"],0)
  expect_equal(unname(actual["Band6"]),pnorm(5,lower.tail=FALSE)-pnorm(5+1e-6,lower.tail=FALSE),tolerance=1e-8)
})

test_that("source RSM and native/portable PCM decisions preserve their scoring contracts", {
  x <- decision_fixture(); cuts <- c(Basic=-.513456789,Advanced=.823456789)
  persons <- unique(as.character(x$data$Person))[c(2,1)]
  source <- score_mfrm_persons(x$rsm,persons=persons,cut_scores=cuts)
  baseline <- score_mfrm_persons(x$rsm,persons=persons)
  expect_identical(source$table,baseline$table)
  expect_identical(summary(source)$decisions,source$decisions)
  rows <- x$data[x$data$Person %in% persons,]
  native <- predict_mfrm_units(x$rsm,rows,scoring_quad_points=121,cut_scores=cuts)
  a <- source$decisions; b <- native$decisions
  expect_equal(a$Probability,b$Probability[match(paste(a$Person,a$Type,a$Label),paste(b$Person,b$Type,b$Label))],tolerance=1e-8)

  # These saved fits expose a low-grid scoring failure; cuts must not bypass it.
  fit <- x$pcm$fits[["PCM-JML"]]; rows <- x$pcm$new_data
  expect_error(predict_mfrm_units(fit,rows,cut_scores=cuts),"integration did not pass")
  review <- predict_mfrm_units(fit,rows,cut_scores=cuts,readiness_policy="review")
  expect_true(any(review$decisions$Status=="review_only"))
  prior <- list(mean=.25,sd=1.2)
  rows <- rows[rows$Person %in% unique(rows$Person)[1:2],]
  rows$Event <- as.character(seq_len(nrow(rows)))
  p <- predict_mfrm_units(fit,rows,cut_scores=cuts,scoring_prior=prior,scoring_quad_points=121)
  portable <- score_mfrm_calibration(x$calibration,rows,cut_scores=cuts,scoring_prior=prior,event_id="Event")
  expect_equal(portable$decisions$Probability,p$decisions$Probability,tolerance=1e-10)
  expect_identical(summary(portable)$decisions,portable$decisions)
  expect_identical(summary(p)$decisions,p$decisions)
  expressions <- parse(text=build_mfrm_replay_script(fit,unit_prediction=p)$script)
  selected <- which(vapply(as.list(expressions),function(e)
    is.call(e) && identical(e[[1]],as.name("<-")) && identical(e[[2]],as.name("unit_prediction")),TRUE))
  expect_length(selected,1)
  args <- as.list(expressions[[selected]][[3]])[-1]
  expect_identical(eval(args$cut_scores),cuts)
  expect_identical(eval(args$scoring_prior),prior)
  expect_identical(eval(args$adaptive_quad_points),p$settings$adaptive_quad_points)
  expect_identical(apa_table(portable,which="decisions")$table,portable$decisions)
  expect_identical(apa_table(source,which="decisions")$table,source$decisions)
  expect_true(mfrm_decision_note() %in% summary(p)$notes)
  expect_true(mfrm_decision_note() %in% summary(portable)$notes)
  missing <- rows[1,]; missing$Person <- "missing-only"; missing$Score <- NA_real_
  empty <- score_mfrm_calibration(x$calibration,missing,missing_response="omit",cut_scores=cuts)
  expect_true(all(empty$decisions$Status=="unavailable"))
  expect_true(all(is.na(empty$decisions$Probability)))
  expect_identical(summary(empty)$decisions,empty$decisions)
})

test_that("portable appendices retain decisions and unscored Persons without rescoring", {
  x <- decision_fixture(); rows <- x$pcm$new_data
  rows <- rows[rows$Person %in% unique(rows$Person)[1:2],]
  rows$Event <- as.character(seq_len(nrow(rows)))
  missing <- rows[1,]; missing$Person <- "missing-only"; missing$Score <- NA_real_
  rows <- rbind(rows, missing)
  cuts <- c(Low=-.513456789,High=.823456789)
  scored <- score_mfrm_calibration(x$calibration, rows, event_id="Event",
    missing_response="omit", cut_scores=cuts)
  no_cuts <- score_mfrm_calibration(x$calibration, rows, event_id="Event",
    missing_response="omit")
  local_mocked_bindings(mfrmr_score_calibration=function(...) stop("No rescoring"),
    fit_mfrm=function(...) stop("No refitting"),.package="mfrmr")
  bundle <- build_summary_table_bundle(scored)
  expect_identical(bundle$tables$decisions, scored$decisions)
  expect_identical(bundle$table_index$Role[match(c("decisions","review"),bundle$table_index$Table)],
    c("conditional_decisions","score_dispositions"))
  expect_identical(apa_table(bundle,which="decisions")$table, scored$decisions)
  expect_identical(build_summary_table_bundle(summary(scored))$tables, bundle$tables)
  expect_true("missing-only" %in% bundle$tables$review$Person)
  expect_true(all(is.na(bundle$tables$decisions$Probability[
    bundle$tables$decisions$Person=="missing-only"])))
  expect_false("decisions" %in% names(build_summary_table_bundle(no_cuts)$tables))
  expect_true("decisions" %in% names(build_summary_table_bundle(scored,
    appendix_preset="recommended")$tables))
  folder <- tempfile(); dir.create(folder); withr::defer(unlink(folder,recursive=TRUE))
  export_summary_appendix(scored,output_dir=folder,include_html=FALSE)
  files <- list.files(folder,"decisions.*[.]csv$",recursive=TRUE,full.names=TRUE)
  expect_length(files,1)
  exported <- read.csv(files[1])
  expect_equal(exported$Probability,scored$decisions$Probability,tolerance=1e-14)
  expect_equal(exported$Lower,scored$decisions$Lower,tolerance=1e-14)
  broken <- summary(scored); broken$settings$cut_scores[1] <- -2
  expect_error(build_summary_table_bundle(broken),"Saved decision probabilities")
  expect_output(print(summary(scored)),"first 10; full table in")
})

test_that("ordinary missing-only and zero-variance states stay explicit", {
  x <- decision_fixture(); input <- mfrm_ordinary_response_input(x$rsm)
  input$assigned_data <- rbind(input$assigned_data,input$assigned_data[1,])
  last <- nrow(input$assigned_data)
  input$assigned_data[last,input$columns$person] <- "missing-only"
  input$assigned_data[last,input$columns$score] <- NA_real_
  input$mean <- .6; input$sd <- .9
  local_mocked_bindings(mfrm_ordinary_response_input=function(...) input,.package="mfrmr")
  cuts <- c(-1,.25,1e300)
  prior <- score_mfrm_persons(x$rsm,persons="missing-only",cut_scores=cuts)
  expect_true(all(prior$decisions$Status=="prior_only"))
  expect_equal(prior$decisions$Probability[1:3],pnorm(cuts,.6,.9,lower.tail=FALSE),tolerance=1e-12)
  expect_identical(summary(prior)$decisions,prior$decisions)
  level <- 1 - .Machine$double.eps / 2
  high <- score_mfrm_persons(x$rsm, persons="missing-only", level=level, cut_scores=cuts)
  expect_equal(high$table$Lower, .6-.9*sqrt(qchisq(level, 1)))
  expect_equal(high$table$Upper, .6+.9*sqrt(qchisq(level, 1)))
  expect_identical(high$decisions, prior$decisions)
  input$sd <- 0
  absent <- score_mfrm_persons(x$rsm,persons="missing-only",cut_scores=cuts)
  expect_true(all(absent$decisions$Status=="unavailable"))
  expect_true(all(is.na(absent$decisions$Probability)))
  expect_identical(summary(absent)$decisions,absent$decisions)
})

test_that("regression priors and saved decision exports retain exact cuts without rescoring", {
  fit <- readRDS(test_path("fixtures","native-location-pcm-regression.rds"))$fit
  people <- fit$population$person_table[1:2,]
  rows <- fit$prep$data[fit$prep$data$Person %in% people$Person,c("Person","Rater","Criterion","Score")]
  cuts <- c(Low=-.513456789,High=.823456789)
  p <- predict_mfrm_units(fit,rows,person_data=people,scoring_quad_points=61,cut_scores=cuts)
  expect_equal(p$decisions$PriorMean,p$estimates$PriorMean[match(p$decisions$Person,p$estimates$Person)])
  path <- tempfile(fileext=".rds"); withr::defer(unlink(path)); saveRDS(p,path)
  local_mocked_bindings(compute_person_posterior_summary=function(...) stop("No rescoring"),
    fit_mfrm=function(...) stop("No refitting"),.package="mfrmr")
  loaded <- readRDS(path)
  expect_identical(summary(loaded)$decisions,p$decisions)
  expect_identical(apa_table(loaded,which="decisions")$table,p$decisions)
  folder <- tempfile(); dir.create(folder); withr::defer(unlink(folder,recursive=TRUE))
  export_summary_appendix(loaded,output_dir=folder,include_html=FALSE)
  files <- list.files(folder,"decisions.*[.]csv$",recursive=TRUE,full.names=TRUE)
  expect_length(files,1)
  expect_equal(read.csv(files[1])$Probability,p$decisions$Probability,tolerance=1e-14)
  expect_equal(read.csv(files[1])$Lower,p$decisions$Lower,tolerance=1e-14)
  for (field in c("decisions","settings")) {
    broken <- p
    if(field=="decisions") broken$decisions$Probability[1] <- .99 else broken$settings$cut_scores[1] <- -2
    expect_error(summary(broken),"Saved decision probabilities")
    expect_error(apa_table(broken,which="decisions"),"Saved decision probabilities")
  }
})
