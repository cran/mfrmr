fair_bounded_case <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      d <- load_mfrmr_data("example_core")
      f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
        method="JML", rating_min=1, rating_max=4, category_policy="preserve",
        optimizer="BFGS", reltol=1e-12, maxit=500)
      p <- unique(d[c("Rater", "Criterion")]); p$Weight <- seq_len(nrow(p))
      saved <<- list(data=d, fit=f, plan=p)
    }
    saved
  }
})

test_that("bounded intervals use actual scoring responses and preserve JML points", {
  x <- fair_bounded_case(); f <- x$fit
  fair <- fair_average_table(f, reference_plan=x$plan)
  persons <- f$config$person_levels[1:2]
  set.seed(55); rng <- .Random.seed
  ci <- confint(fair, parm=persons, method="bounded")
  expect_identical(.Random.seed, rng)
  expect_true(all(ci$Available & ci$CIEligible))
  expect_true(all(is.na(ci$SE) & is.na(ci$Planned)))
  expect_equal(ci$Estimate, fair$reference_scores$Estimate[match(persons, fair$reference_scores$Person)])
  expect_true(all(attr(ci, "support")$CompletePlan))
  expect_identical(attr(ci,"settings")$score_support,
    list(range=c(lower=1,upper=4),provenance="declared_rubric",coding="identity"))
  # Independent direct cell lookup, not the production key builder.
  w <- x$plan$Weight/sum(x$plan$Weight)
  expected <- vapply(persons, function(person) sum(w * vapply(seq_len(nrow(x$plan)), function(i)
    mean(x$data$Score[x$data$Person == person & x$data$Rater == x$plan$Rater[i] &
      x$data$Criterion == x$plan$Criterion[i]]), 0)), 0)
  expect_equal(ci$DesignEstimate, unname(expected), tolerance=1e-12)
  radius <- diff(range(x$data$Score)) * sqrt(log(40)/2 * sum(w^2))
  expect_equal(ci$SamplingRadius, rep(radius, 2), tolerance=1e-12)
  expect_equal(ci$CI_Lower, pmax(min(x$data$Score), expected-radius), ignore_attr=TRUE)
  expect_equal(ci$CI_Upper, pmin(max(x$data$Score), expected+radius), ignore_attr=TRUE)
  L <- matrix(c(1,-1), 1, dimnames=list("difference",persons))
  difference <- confint(fair, parm=persons, method="bounded", contrasts=L)
  expect_equal(difference$SamplingRadius, sqrt(2)*radius, tolerance=1e-12)
  expect_equal(difference$DesignEstimate, unname(expected[1]-expected[2]))
  expect_match(plot(difference, draw=FALSE)$data$xlab, "contrast")
  expect_false(grepl("contrast", plot(ci, draw=FALSE)$data$xlab))
  expect_match(plot(ci, draw=FALSE)$data$notes$Text, "independent bounded")
  expect_output(print(ci), "bounded responses")
  path <- tempfile(fileext=".rds"); withr::defer(unlink(path))
  saveRDS(ci,path); expect_identical(readRDS(path),ci)
  expect_identical(summary(readRDS(path)),ci)
  expect_error(confint(fair, method="bounded", nsim=2), "omit")
  expect_error(confint(fair, method="bounded", anchor_uncertainty="fixed"), "omit")
  expect_error(confint(fair, method="bounded", gradient_tol=1), "omit")
})

test_that("extreme Persons and estimated anchors do not become degenerate generators", {
  x <- fair_bounded_case(); d <- x$data; persons <- x$fit$config$person_levels[1:3]
  d$Score[d$Person == persons[1]] <- min(d$Score)
  d$Score[d$Person == persons[2]] <- max(d$Score)
  f <- fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",method="JML",
    rating_min=1,rating_max=4,category_policy="preserve",
    optimizer="BFGS",reltol=1e-12,maxit=1000)
  fair <- fair_average_table(f, reference_plan=x$plan)
  ci <- confint(fair, parm=persons, method="bounded")
  expect_true(all(is.infinite(f$facets$person$PrimaryEstimate[1:2])))
  expect_true(all(ci$Available))
  expect_equal(ci$Estimate[1:2], range(d$Score))
  expect_gt(ci$CI_Upper[1],ci$CI_Lower[1])
  expect_gt(ci$CI_Upper[2],ci$CI_Lower[2])
  L <- matrix(c(-1,1,0),1,dimnames=list("high minus low",persons))
  difference <- confint(fair,parm=persons,contrasts=L,method="bounded")
  expect_gt(difference$CI_Lower,0)
  expect_equal(difference$CI_Upper,diff(range(d$Score)))
  a <- subset(f$facets$others, Facet == "Rater")[1,c("Facet","Level","Estimate")]
  names(a)[3] <- "Anchor"
  args <- f$config$replay_inputs; args$package_version <- NULL; args$data <- d; args$anchors <- a
  anchored <- do.call(fit_mfrm,args)
  bounded <- confint(fair_average_table(anchored,reference_plan=x$plan),parm=persons,method="bounded")
  expect_equal(bounded$CI_Lower,ci$CI_Lower,tolerance=0)
  expect_equal(bounded$CI_Upper,ci$CI_Upper,tolerance=0)
  expect_identical(attr(bounded,"settings")$anchor_uncertainty,"anchors_not_used_by_bound")
  # Existing default and delta refusals are not relaxed.
  expect_false(confint(fair,parm=persons[3],method="delta")$Available)
})

test_that("bounded inference requires declared support and a fixed score coding", {
  x <- fair_bounded_case()
  args <- x$fit$config$replay_inputs; args$package_version <- NULL; args$data <- x$data
  for (endpoint in list(c("rating_min","rating_max"), "rating_min", "rating_max")) {
    a <- args; a[endpoint] <- NULL
    f <- do.call(fit_mfrm,a)
    fair <- fair_average_table(f,reference_plan=x$plan)
    expect_error(confint(fair,method="bounded"),"declared rating_min and rating_max")
  }
  # Old objects cannot acquire a known-support guarantee from missing metadata.
  f <- x$fit; f$prep$rating_min_source <- NULL
  expect_error(mfrm_fair_bounded_support(f),"declared rating_min and rating_max")
  # An unchanged numeric range alone is insufficient after category recoding.
  a <- args; a$data$Score[a$data$Score == 2] <- 1
  a$category_policy <- "collapse"; a$keep_original <- FALSE
  expect_warning(f <- do.call(fit_mfrm,a),"recoded internally")
  expect_false(identical(f$prep$score_map$OriginalScore,f$prep$score_map$InternalScore))
  fair <- fair_average_table(f,reference_plan=x$plan)
  expect_error(confint(fair,method="bounded"),"unchanged full rubric coding")
  # The known upper endpoint remains 4 even when this sample never observes it.
  a <- args; a$data$Score[a$data$Score == 4] <- 3
  expect_warning(f <- do.call(fit_mfrm,a),"Category support is retained")
  expect_equal(max(f$prep$data$Score),3)
  ci <- confint(fair_average_table(f,reference_plan=x$plan),
    parm=f$config$person_levels[1:2],method="bounded")
  expect_true(all(ci$Available & ci$CIEligible))
  expect_equal(ci$SupportUpper,c(4,4))
  expect_equal(attr(ci,"settings")$score_support$range,c(lower=1,upper=4))
})

test_that("sample extrema are not valid known bounds in a finite response law", {
  # No simulation: six independent ratings, known rubric 0:3. On this event
  # both sample categories occur, but clipping to [0,1] excludes E[Y] = 1.31.
  prob <- c(.10,.69,.01,.20)
  truth <- sum((0:3)*prob)
  patterns <- as.matrix(expand.grid(rep(list(0:1),6)))
  event <- rowSums(patterns)>0 & rowSums(patterns)<6
  mass <- apply(patterns,1,function(y)prod(prob[y+1L]))
  failure <- sum(mass[event])
  expect_gt(truth,1)
  expect_equal(failure,sum(prob[1:2])^6-sum(prob[1:2]^6),tolerance=1e-14)
  expect_gt(failure,.05)
})

test_that("missing plan cells and duplicate observations retain the correct weights", {
  data <- data.frame(Person=c("p","p","q"), Rater=c("a:b","a:b","a"),
    Task=c("c","c","b:c"), Score=c(2,4,3), Weight=1)
  plan <- data.frame(Rater=c("a:b","a"),Task=c("c","b:c"),NormalizedWeight=c(.25,.75))
  design <- mfrm_fair_bounded_weights(data,plan,c("Rater","Task"),c("p","q"))
  expect_equal(as.matrix(design$weights),matrix(c(.125,.125,0,0,0,.75),3,2))
  expect_equal(design$support$UnobservedPlanWeight,c(.75,.25))
  expect_false(any(design$support$CompletePlan))
  f <- list(prep=list(data=data,rating_min=2,rating_max=4,
    rating_min_source="declared",rating_max_source="declared",
    score_map=data.frame(OriginalScore=2:4,InternalScore=2:4)),config=list(facet_names=c("Rater","Task")))
  local_mocked_bindings(prediction_fit_identity=function(...) "fixture",.package="mfrmr")
  L <- diag(2); dimnames(L) <- list(c("p","q"),c("p","q"))
  ci <- mfrm_fair_bounded_result(f,plan,c("p","q"),L,c(3,3),.95,FALSE)
  expect_true(all(is.na(ci$DesignEstimate)))
  expect_equal(ci$UnobservedLower,c(1.5,.5))
  expect_equal(ci$UnobservedUpper,c(3,1))
  # A plan entirely absent for p gives the full score range, not false precision.
  absent <- plan[2,,drop=FALSE]; absent$NormalizedWeight <- 1
  full <- mfrm_fair_bounded_result(f,absent,"p",matrix(1,1,1,dimnames=list("p","p")),3,.95,FALSE)
  expect_equal(c(full$CI_Lower,full$CI_Upper),c(2,4))
  expect_true(full$Available && full$CIEligible)
  expect_false(full$Informative)
  expect_identical(full$Status,"uninformative_bound")
  base <- mfrm_fair_bounded_result(f,plan,c("p","q"),L,c(3,3),.95,FALSE)
  for (multiplier in c(1e-200,1e200)) {
    scaled <- mfrm_fair_bounded_result(f,plan,c("p","q"),L*multiplier,c(3,3)*multiplier,.95,FALSE)
    expect_equal(scaled$SamplingRadius/multiplier,base$SamplingRadius,tolerance=1e-12)
    expect_equal(scaled$CI_Lower/multiplier,base$CI_Lower,tolerance=1e-12)
    expect_equal(scaled$CI_Upper/multiplier,base$CI_Upper,tolerance=1e-12)
  }
  f$prep$data$Weight[1] <- 2
  expect_error(mfrm_fair_bounded_result(f,plan,c("p","q"),L,c(3,3),.95,FALSE),"unit weights")
})

test_that("finite-pattern enumeration respects weighted and signed coverage bounds", {
  local_mocked_bindings(prediction_fit_identity=function(...) "enumeration",.package="mfrmr")
  patterns <- as.matrix(expand.grid(rep(list(0:1),6)))
  d <- data.frame(Person=rep(c("p","q"),each=3),Item=rep(c("a","a","b"),2),Score=0,Weight=1)
  L <- rbind(p=c(1,0),q=c(0,1),difference=c(1,-1)); colnames(L) <- c("p","q")
  probabilities <- list(rep(.5,6),c(.05,.05,.9,.8,.8,.2),c(.99,.99,.8,.1,.1,.01))
  for (missing in c(FALSE,TRUE)) {
    plan <- if (missing) data.frame(Item=c("a","b","c"),NormalizedWeight=c(.3,.5,.2)) else
      data.frame(Item=c("a","b"),NormalizedWeight=c(.4,.6))
    for (level in c(.5,.95)) {
      bounds <- lapply(seq_len(nrow(patterns)),function(i) {
        d$Score <- patterns[i,]
        f <- list(prep=list(data=d,rating_min=0,rating_max=1,
          rating_min_source="declared",rating_max_source="declared",
          score_map=data.frame(OriginalScore=0:1,InternalScore=0:1)),config=list(facet_names="Item"))
        mfrm_fair_bounded_result(f,plan,c("p","q"),L,rep(NA_real_,3),level,TRUE)
      })
      lower <- do.call(rbind,lapply(bounds,`[[`,"CI_Lower")); upper <- do.call(rbind,lapply(bounds,`[[`,"CI_Upper"))
      design <- mfrm_fair_bounded_weights(d,plan,"Item",c("p","q"))
      for (p in probabilities) for (unseen in list(c(0,1),c(.3,.7),c(1,0))) {
        mass <- apply(patterns,1,function(y)prod(ifelse(y==1,p,1-p)))
        truth <- drop(L %*% (as.numeric(Matrix::crossprod(design$weights,p)) +
          design$support$UnobservedPlanWeight*unseen))
        covered <- sweep(lower,2,truth,"<=") & sweep(upper,2,truth,">=")
        coverage <- colSums(covered*mass)
        expect_gte(min(coverage),level-1e-12)
        expect_equal(sum(mass),1,tolerance=1e-12)
      }
    }
  }
})
