test_that("anchor origins reproduce native sparse RSM/PCM constraints", {
  input <- read.csv(test_path("fixtures", "tam-jml-0.2.4-response.csv"))
  reference <- read.csv(test_path("fixtures", "facets-anchor-origin.csv"))
  for (model in c("RSM", "PCM")) {
    d <- input[startsWith(input$Case, model), ]
    anchor <- data.frame(Facet="Person", Level="P03", Anchor=-.4)
    fit <- fit_mfrm(d,"Person",c("Rater","Task"),"Score",method="JML",
      model=model,step_facet=if(model=="PCM") "Task" else NULL,
      rating_min=0,rating_max=3,anchors=anchor,noncenter_facet="Rater",
      centering_scope="anchor_origin",maxit=500,reltol=1e-11,optimizer="BFGS")
    ref <- reference[reference$Model==model, ]
    pa <- mfrmr:::expand_params(fit$opt$par,mfrmr:::build_param_sizes(fit$config),fit$config)
    estimates <- c(pa$theta,pa$facets$Rater,pa$facets$Task)
    keys <- paste(rep(c("Person","Rater","Task"),c(length(pa$theta),
      length(pa$facets$Rater),length(pa$facets$Task))),names(estimates))
    idx <- match(keys,paste(ref$Facet,ref$Level))
    expect_false(anyNA(idx))
    expect_equal(length(unique(idx)),nrow(ref))
    expect_lt(max(abs(unname(estimates)-ref$Native[idx])), 1e-4)
    expect_equal(pa$theta["P03"],c(P03=-.4))
    expect_false(fit$config$theta_spec$centered)
    expect_false(fit$config$facet_specs$Rater$centered)
    expect_true(fit$config$facet_specs$Task$centered)
    expect_gt(abs(mean(pa$facets$Rater)),1)
    expect_equal(sum(pa$facets$Task),0,tolerance=1e-12)
    # Independent score equations, including retained missing assignments.
    eta <- pa$theta[d$Person]-pa$facets$Rater[as.character(d$Rater)]-
      pa$facets$Task[as.character(d$Task)]
    steps <- if(model=="PCM") pa$steps_mat[as.character(d$Task),,drop=FALSE] else
      matrix(rep(pa$steps,each=nrow(d)),nrow(d))
    logits <- cbind(0,eta-steps[,1],2*eta-rowSums(steps[,1:2,drop=FALSE]),3*eta-rowSums(steps))
    pr <- exp(logits-apply(logits,1,max));pr<-pr/rowSums(pr)
    expected <- drop(pr %*% (0:3))
    scores <- tapply(d$Score-expected,d$Person,sum)
    expect_lt(max(abs(scores[names(scores)!="P03"])),1e-4)
    expect_lt(max(abs(tapply(d$Score-expected,d$Rater,sum))),1e-4)
    before <- serialize(fit,NULL)
    dx <- diagnose_mfrm(fit,residual_pca="none",diagnostic_mode="legacy")
    fair <- fair_average_table(fit,dx)
    expect_no_error(summary(fair))
    sm <- summary(fit,compute="never",detail="brief")
    expect_identical(sm$settings_overview$CenteringScope,"anchor_origin")
    expect_output(print(sm),"[Dd]irect or group anchors set")
    review <- fit$config$anchor_review
    fs <- review$facet_summary
    expect_identical(fs$OriginBasis[fs$Facet=="Person"],"direct_anchors")
    expect_false(fs$Centered[fs$Facet=="Person"])
    expect_match(paste(review$recommendations,collapse=" "),"Direct or group anchors set")
    recipe <- build_mfrm_replay_script(fit,diagnostics=dx)
    expect_match(recipe$script,'centering_scope = "anchor_origin"',fixed=TRUE)
    expect_identical(serialize(fit,NULL),before)
  }
})

test_that("a Rater origin preserves the likelihood and allows nonzero fixed totals", {
  d <- load_mfrmr_data("example_core")
  args <- list(data=d,person="Person",facets=c("Rater","Criterion"),score="Score",
    method="JML",optimizer="BFGS",maxit=500,reltol=1e-11)
  raw <- do.call(fit_mfrm,args)
  args$anchors <- data.frame(Facet="Rater",Level="R01",Anchor=.7)
  args$centering_scope <- "anchor_origin"
  fit <- do.call(fit_mfrm,args)
  r <- subset(raw$facets$others,Facet=="Rater")
  a <- subset(fit$facets$others,Facet=="Rater")
  shift <- .7-r$Estimate[r$Level=="R01"]
  expect_equal(a$Estimate,r$Estimate+shift,tolerance=1e-4)
  expect_equal(fit$facets$person$Estimate,raw$facets$person$Estimate+shift,tolerance=1e-4)
  expect_equal(fit$opt$value,raw$opt$value,tolerance=1e-7)
  expect_equal(unname(a$Estimate[a$Level=="R01"]),.7,tolerance=1e-12)
  # Fully fixed nonzero totals are permitted because this facet is not centered.
  args$anchors <- data.frame(Facet="Rater",Level=a$Level,Anchor=a$Estimate)
  fixed <- do.call(fit_mfrm,args)
  expect_equal(subset(fixed$facets$others,Facet=="Rater")$Estimate,a$Estimate,tolerance=1e-12)
  expect_equal(fixed$opt$value,fit$opt$value,tolerance=1e-7)
})

test_that("anchor origins are explicit and unsupported combinations stop", {
  input <- read.csv(test_path("fixtures","tam-jml-0.2.4-response.csv"))
  d <- input[input$Case=="RSM_sparse", ]
  anchor <- data.frame(Facet="Person",Level="P03",Anchor=-.4)
  args <- list(data=d,person="Person",facets=c("Rater","Task"),score="Score",
    method="JML",anchors=anchor,centering_scope="anchor_origin")
  expect_error(do.call(fit_mfrm,args),"unanchored, active")
  args$noncenter_facet <- "Rater"
  for(method in c("MML","CML")) {
    args$method <- method
    expect_error(do.call(fit_mfrm,args),"ordinary additive RSM/PCM JML")
  }
  args$method <- "JML"; args$model <- "GPCM"
  expect_error(do.call(fit_mfrm,args),"ordinary additive RSM/PCM JML")
  args$model <- "RSM"; args$facet_interactions <- "Rater:Task"
  expect_error(do.call(fit_mfrm,args),"ordinary additive RSM/PCM JML")
  expect_error(review_mfrm_anchors(d,"Person",c("Rater","Task"),"Score",
    group_anchors=data.frame(Facet="Person",Level=c("P03","P04"),Group="g",GroupValue=0),
    noncenter_facet="Person",centering_scope="anchor_origin"),"unanchored, active")
  # Existing constraints still describe a different, intentionally restricted model.
  legacy <- review_mfrm_anchors(d,"Person",c("Rater","Task"),"Score",anchors=anchor,
    noncenter_facet="Rater",centering_scope="facet")
  expect_true(legacy$facet_summary$Centered[legacy$facet_summary$Facet=="Person"])
  expect_identical(eval(formals(fit_mfrm)$centering_scope)[1],"facet")
})

test_that("an extreme sole anchor cannot supply a false finite calibration", {
  d <- expand.grid(Person=c("A","B","C"),Rater=c("R1","R2"),Task=c("T1","T2"),
    stringsAsFactors=FALSE)
  d$Score <- c(0,0,1,0,1,0,0,1,0,0,0,1)
  # Increasing free Persons and both Raters preserves their predictors while
  # driving anchored all-zero A toward probability one: no finite joint MLE.
  expect_error(fit_mfrm(d,"Person",c("Rater","Task"),"Score",method="JML",
    anchors=data.frame(Facet="Person",Level="A",Anchor=0),noncenter_facet="Rater",
    centering_scope="anchor_origin",maxit=200),"No finite joint maximum")
})

test_that("same-boundary Person anchors stop before optimization or LP auditing", {
  local_mocked_bindings(
    run_mfrm_optimization = function(...) stop("Unexpected optimization"),
    audit_mfrm_jml_joint_recession = function(...) stop("Unexpected LP audit"),
    .package = "mfrmr")
  input <- read.csv(test_path("fixtures", "tam-jml-0.2.4-response.csv"))
  for (model in c("RSM", "PCM")) for (score in c(0L, 3L)) {
    d <- input[startsWith(input$Case, model), ]
    d$Score[d$Person == "P03"] <- score
    expect_error(fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
      method = "JML", model = model,
      step_facet = if (model == "PCM") "Task" else NULL,
      rating_min = 0, rating_max = 3,
      anchors = data.frame(Facet = "Person", Level = "P03", Anchor = -.4),
      noncenter_facet = "Rater", centering_scope = "anchor_origin"),
      "No finite joint maximum")
  }
})

test_that("the analytic anchor boundary uses retained support and allows either facet sign", {
  local_mocked_bindings(run_mfrm_optimization = function(...) stop("Unexpected optimization"),
    .package = "mfrmr")
  input <- read.csv(test_path("fixtures", "tam-jml-0.2.4-response.csv"))
  d <- input[input$Case == "RSM_sparse", ]
  d$Score <- d$Score + 2L
  d$Score[d$Person %in% c("P03", "P04")] <- 2L
  d$Weight <- .5
  # An omitted interior response must not conceal the same-boundary support.
  omitted <- d[which(d$Person == "P03")[1], ]
  omitted$Score <- 3L
  omitted$Weight <- 0
  d <- rbind(d, omitted)
  expect_error(suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
    method = "JML", weight = "Weight", positive_facets = "Rater",
    rating_min = 2, rating_max = 5,
    anchors = data.frame(Facet = "Person", Level = c("P03", "P04"), Anchor = c(-.4, .7)),
    noncenter_facet = "Rater", centering_scope = "anchor_origin")),
    "declared minimum.*No finite joint maximum")
})

test_that("additional anchor support is not mistaken for a same-boundary recession", {
  d <- expand.grid(Person = c("A", "B", "C", "D"), Rater = c("R1", "R2"),
    Task = c("T1", "T2"), stringsAsFactors = FALSE)
  d$Score <- c(0,0,0,1, 0,1,1,0, 0,0,1,1, 0,1,0,0)
  anchors <- data.frame(Facet = "Person", Level = c("A", "B"), Anchor = c(-.4, .7))
  for (opposite in c(FALSE, TRUE)) {
    ratings <- d
    if (opposite) ratings$Score[ratings$Person == "B"] <- 1L
    fit <- fit_mfrm(ratings, "Person", c("Rater", "Task"), "Score",
      method = "JML", rating_min = 0, rating_max = 1, anchors = anchors,
      noncenter_facet = "Rater", centering_scope = "anchor_origin",
      optimizer = "BFGS", maxit = 500, reltol = 1e-11)
    expect_s3_class(fit, "mfrm_fit")
    expect_true(all(is.finite(fit$facets$others$Estimate)))
    expect_equal(fit$facets$person$Estimate[match(anchors$Level,fit$facets$person$Person)],
      anchors$Anchor)
  }
})

test_that("supported extreme anchors retain score rows and PCM interval derivatives", {
  input <- read.csv(test_path("fixtures", "tam-jml-0.2.4-response.csv"))
  ref <- read.csv(test_path("fixtures", "facets-anchor-origin.csv"))
  for (model in c("RSM", "PCM")) for (endpoint in c(0L, 3L)) {
    d <- input[startsWith(input$Case, model), ]
    d$Score[d$Person == "P03"] <- endpoint
    anchor <- data.frame(Facet = "Person", Level = c("P03", "P24"),
      Anchor = c(-.4, ref$Native[ref$Model == model & ref$Facet == "Person" & ref$Level == "P24"]))
    fit <- fit_mfrm(d, "Person", c("Rater", "Task"), "Score", method = "JML",
      model = model, step_facet = if (model == "PCM") "Task" else NULL,
      rating_min = 0, rating_max = 3, anchors = anchor, noncenter_facet = "Rater",
      centering_scope = "anchor_origin", optimizer = "BFGS", maxit = 500, reltol = 1e-11)
    dx <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy")
    before <- serialize(list(fit, dx), NULL)
    for (convention in c("display", "facets")) for (correction in c(0, .3)) {
      result <- fair_average_table(fit, dx, totalscore = FALSE,
        xtreme = correction, xtreme_method = convention)
      for (facet in c("Person", "Rater", "Task")) {
        tab <- result$raw_by_facet[[facet]]
        expect_equal(tab$TotalCount, as.numeric(table(d[[facet]])[tab$Level]))
        expect_equal(tab$TotalScore, as.numeric(tapply(d$Score, d[[facet]], sum)[tab$Level]))
      }
      anchored <- subset(result$raw_by_facet$Person, Level == "P03")
      expect_equal(anchored$PrimaryMeasure, -.4)
      expect_equal(anchored$Measure, -.4)
      expect_equal(anchored$ExtremeAdjustment, 0)
      full <- fair_average_table(fit, dx, totalscore = TRUE,
        xtreme = correction, xtreme_method = convention)
      expect_equal(result$raw_by_facet$Rater$FairM, full$raw_by_facet$Rater$FairM)
    }
    # The conditional plot SE must use the same retained PCM scale mixture.
    if (model == "PCM") {
      pa <- mfrmr:::expand_params(fit$opt$par,mfrmr:::build_param_sizes(fit$config),fit$config)
      for (metric in c("FairM", "FairZ")) {
        shown <- plot(fit, type = "fair_average", diagnostics = dx, facet = "Rater",
          metric = metric, totalscore = FALSE, show_ci = TRUE, draw = FALSE)$data$data
        derivative <- vapply(shown$Level, function(level) {
          rows <- as.character(d$Rater) == level
          eta <- -pa$facets$Rater[level] + if (metric == "FairM")
            mean(pa$theta)-mean(pa$facets$Task) else 0
          variance <- vapply(as.character(d$Task[rows]), function(task) {
            k <- 0:3; logits <- k*eta-c(0,cumsum(pa$steps_mat[task,]))
            pr <- exp(logits-max(logits));pr<-pr/sum(pr)
            sum(pr*(k-sum(k*pr))^2)
          }, numeric(1))
          mean(variance)
        }, numeric(1))
        expect_equal(shown$CI_SE, unname(derivative)*shown$ModelSE, tolerance = 1e-8)
      }
    }
    expect_identical(serialize(list(fit, dx), NULL), before)
  }
})
