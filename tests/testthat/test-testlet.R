testlet_fixture <- function() {
  y <- rbind(c(0,1,2,2,1,0), c(2,2,1,1,0,0), c(1,0,1,2,2,1),
    rep(0,6), rep(2,6), c(NA,NA,NA,1,2,0))
  data.frame(Person = rep(paste0("P",1:6),each=6),
    Rater = rep(rep(c("A","B"),each=3),6), Criterion = rep(1:3,12), Score = as.vector(t(y)))
}

test_that("selected testlet Persons retain the complete kernel coordinates", {
  d <- testlet_fixture()
  reference <- mfrmr:::mfrm_testlet_data(d, "Person", "Score", "Rater",
    c("Rater", "Criterion"), 0:2, "omit")
  empty <- d[1:6, ]; empty$Person <- "Unobserved"; empty$Score <- NA_real_
  input <- mfrmr:::mfrm_testlet_data(rbind(d, empty), "Person", "Score", "Rater",
    c("Rater", "Criterion"), 0:2, "omit", reference = reference)
  rule <- mfrmr:::gauss_hermite_normal(7L)
  for (variance in c(0, .7)) for (gradient in c(FALSE, TRUE)) {
    par <- c(.15 * sin(seq_len(ncol(input$X))), -.6, .6, variance)
    full <- mfrmr:::mfrm_testlet_kernel(input, par, c(-1, .2), rule, gradient)
    for (persons in list(6L, c(6L, 1L), c(6L, 1L, 6L), 7L,
      integer(), rev(seq_along(input$persons)))) {
      expected <- full
      expected$loglik <- full$loglik[, persons, drop = FALSE]
      if (gradient) {
        expected$score <- full$score[, persons, , drop = FALSE]
        expected$theta_score <- full$theta_score[, persons, drop = FALSE]
        expected$theta_curvature <- full$theta_curvature[, persons, drop = FALSE]
      }
      actual <- mfrmr:::mfrm_testlet_kernel(input, par, c(-1, .2), rule,
        gradient, persons)
      expect_identical(actual, expected)
    }
  }
})

test_that("testlet preparation preserves assignments and separates local from fixed roles", {
  d <- testlet_fixture()
  expect_error(fit_mfrm_testlet(d,"Person","Score","Rater",c("Rater","Criterion"),0:2), "Assigned scores")
  input <- mfrm_testlet_data(d,"Person","Score","Rater",c("Rater","Criterion"),0:2,"omit")
  expect_equal(input$omitted_rows, 31:33)
  expect_equal(nrow(input$assigned_data),36)
  expect_equal(nrow(input$data),33)
  expect_equal(nrow(input$blocks),11)
  expect_equal(colSums(input$basis$Criterion), c(C1=0,C2=0), tolerance=1e-12)
  bad <- d; bad$Other <- bad$Rater
  expect_error(mfrm_testlet_data(bad,"Person","Score","Rater",c("Rater","Other"),0:2,"omit"),"aliased")
  bad <- d; bad$Rater[31] <- NA
  expect_error(mfrm_testlet_data(bad,"Person","Score","Rater","Criterion",0:2,"omit"),"Identifiers")
  bad <- d; bad$Block <- seq_len(nrow(bad))
  expect_error(mfrm_testlet_data(bad,"Person","Score","Block","Criterion",0:2,"omit"),"repeated")
  expect_error(mfrm_testlet_data(d,"Person","Score","Rater","Criterion",c(0,2),"omit"),"consecutive")
  expect_error(fit_mfrm_testlet(d,"Person","Score","Rater",score_levels=0:2,testlet_variance=-1),"nonnegative")
  expect_error(fit_mfrm_testlet(d,"Person","Score","Rater",score_levels=0:2,quad_points=3),"quad_points")
})

test_that("three-block four-category likelihood and derivatives match an independent tensor integral", {
  d <- expand.grid(Criterion=c("c1","c2"),Station=c("s1","s2","s3"),Person=c("a:b","a"))
  d$Block <- d$Station; d$Score <- c(0,1,2,3,1,0,3,2,1,1,0,2)
  d$Score[2] <- NA
  x <- mfrm_testlet_data(d,"Person","Score","Block",c("Station","Criterion"),0:3,"omit")
  par <- c(as.vector(crossprod(x$basis$Station,c(-.2,.1,.1))),
    as.vector(crossprod(x$basis$Criterion,c(-.3,.3))),-.8,0,.6,.4)
  # Both independent tensor and nested integration use the same rule. Seven
  # nodes retain every derivative/variance-face check without a 13^4 grid.
  order <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) 13L else 7L
  jac <- matrix(0,order,order)
  jac[cbind(1:(order-1),2:order)] <- sqrt(1:(order-1))
  eig <- eigen(jac+t(jac),symmetric=TRUE); nodes <- eig$values; weights <- eig$vectors[1,]^2
  grid <- as.matrix(expand.grid(rep(list(seq_len(order)),4)))
  joint_weight <- apply(matrix(weights[grid],nrow(grid)),1,prod)
  direct <- function(p) {
    station <- as.vector(x$basis$Station %*% p[1:2])
    criterion <- as.vector(x$basis$Criterion %*% p[3])
    out <- 0
    for (id in unique(d$Person)) {
      likelihood <- rep(1,nrow(grid))
      for (i in which(d$Person==id & !is.na(d$Score))) {
        s <- match(d$Station[i],rownames(x$basis$Station)); c <- match(d$Criterion[i],rownames(x$basis$Criterion))
        ability_sd <- if (length(p) == 8L) sqrt(p[8]) else 1
        eta <- ability_sd*nodes[grid[,1]] + sqrt(p[7])*nodes[grid[,s+1]] - station[s] - criterion[c]
        w <- cbind(1,exp(eta-p[4]),exp(2*eta-sum(p[4:5])),exp(3*eta-sum(p[4:6])))
        likelihood <- likelihood*(w/rowSums(w))[,d$Score[i]+1]
      }
      out <- out+log(sum(joint_weight*likelihood))
    }
    out
  }
  a <- mfrm_testlet_evaluate(x,par,gauss_hermite_normal(order))
  expect_equal(a$loglik,direct(par),tolerance=1e-11)
  reference_gradient <- vapply(seq_along(par),function(i){ h<-rep(0,length(par));h[i]<-1e-5
    (direct(par+h)-direct(par-h))/2e-5 },numeric(1))
  expect_equal(a$gradient,reference_gradient,tolerance=1e-7)
  zero <- par; zero[length(zero)] <- 0
  f0 <- direct(zero); h <- 1e-4
  a0 <- mfrm_testlet_evaluate(x,zero,gauss_hermite_normal(order))
  one <- zero; two <- zero; one[length(one)] <- h; two[length(two)] <- h/2
  slope <- 2*(direct(two)-f0)/(h/2)-(direct(one)-f0)/h
  expect_equal(tail(a0$gradient,1),slope,tolerance=1e-6)
  # Independent tensor integration also checks the new ability-variance
  # derivative, including its one-sided score on either variance face.
  for (variances in list(c(.4,.49), c(.4,2.25), c(0,.7), c(.4,0), c(0,0))) {
    zpar <- c(par[1:6], variances)
    evaluated <- mfrm_testlet_evaluate(x,zpar,gauss_hermite_normal(order))
    expect_equal(evaluated$loglik,direct(zpar),tolerance=1e-11)
    reference <- vapply(seq_along(zpar),function(i) {
      h <- rep(0,length(zpar)); h[i] <- 1e-4
      if (i >= 7 && zpar[i] == 0) {
        (4*(direct(zpar+h/2)-direct(zpar))-(direct(zpar+h)-direct(zpar)))/h[i]
      } else (direct(zpar+h)-direct(zpar-h))/(2*h[i])
    },numeric(1))
    expect_equal(evaluated$gradient,reference,tolerance=2e-6)
  }
  # Reusing labels across Persons must not introduce shared effects.
  renamed <- d; renamed$Block <- paste(renamed$Person,renamed$Block,sep="|")
  z <- mfrm_testlet_data(renamed,"Person","Score","Block",c("Station","Criterion"),0:3,"omit")
  expect_equal(mfrm_testlet_evaluate(z,par,gauss_hermite_normal(order))$loglik,a$loglik,tolerance=1e-12)
  shuffled <- d[nrow(d):1,]
  z <- mfrm_testlet_data(shuffled,"Person","Score","Block",c("Station","Criterion"),0:3,"omit")
  expect_equal(mfrm_testlet_evaluate(z,par,gauss_hermite_normal(order))$loglik,a$loglik,tolerance=1e-12)
})

test_that("fitting and continuous scoring agree with retained independent reference values", {
  # Extended numerical/reference case; shorter live routes remain on CRAN.
  skip_on_cran()
  d <- testlet_fixture()
  fit <- fit_mfrm_testlet(d,"Person","Score","Rater",c("Rater","Criterion"),0:2,quad_points=61,missing="omit",person_sd=1)
  expect_true(fit$checks$NumericalReady); expect_true(fit$checks$InformationPositive)
  expect_equal(fit$calibration$variance,.98832920,tolerance=2e-5)
  expect_equal(fit$loglik,-30.4046955,tolerance=1e-7)
  expect_length(fit$runs,3)
  expect_equal(fit$covariance,t(fit$covariance),tolerance=1e-12)
  scores <- predict(fit)
  expect_true(all(scores$table$Status=="available_conditional"))
  expect_equal(as.matrix(scores$table[c("Estimate","ConditionalSD")]),fit$moments,tolerance=1e-6,ignore_attr=TRUE)
  expect_false(scores$settings$calibration_uncertainty)
  expect_true(all(scores$table$Lower<scores$table$Upper))
  # Continuous CDF mass at the returned bounds, using a finer local rule.
  row <- scores$table[1,]; rule <- gauss_hermite_normal(181)
  kernel <- function(t) exp(mfrm_testlet_kernel(fit$input,fit$parameters,t,rule,persons=1)$loglik[,1]+dnorm(t,log=TRUE))
  mass <- integrate(kernel,-Inf,Inf,rel.tol=1e-9)$value
  expect_equal(integrate(kernel,-Inf,row$Lower,rel.tol=1e-9)$value/mass,.025,tolerance=1e-7)
  expect_equal(integrate(kernel,row$Upper,Inf,rel.tol=1e-9)$value/mass,.025,tolerance=1e-7)
  new <- d[d$Person=="P2",]; new$Person <- "New"; new$Score <- NA_real_
  prior <- predict(fit,new,missing="omit")
  expect_identical(prior$table$Status,"prior_only")
  expect_equal(unlist(prior$table[c("Estimate","ConditionalSD","Lower","Upper")]),
    c(Estimate=0,ConditionalSD=1,Lower=qnorm(.025),Upper=qnorm(.975)))
  new$Rater <- "Unknown"
  expect_error(predict(fit,new,missing="omit"),"unknown")
  path <- tempfile(); saveRDS(fit,path)
  expect_equal(predict(readRDS(path),d[d$Person=="P1",]),predict(fit,d[d$Person=="P1",]))
  expect_equal(predict(fit,d[d$Person=="P1",])$table$Estimate,scores$table$Estimate[1],tolerance=1e-10)
  device <- dev.cur(); p <- plot(scores,draw=FALSE)
  expect_identical(dev.cur(),device); expect_identical(plot_data(p)$table,scores$table)
  p <- plot(fit,facet="Rater",draw=FALSE)
  expect_identical(plot_data(p)$table$Level,c("A","B"))
  if (requireNamespace("ggplot2", quietly=TRUE)) expect_s3_class(as_ggplot(p),"ggplot")
  bad <- fit; bad$checks$NumericalReady <- FALSE
  expect_error(predict(bad),"checks")
})

test_that("estimated boundaries, unresolved scoring and prior-only rows remain explicit", {
  d <- testlet_fixture(); permutations <- rbind(c(0,1,2),c(0,2,1),c(1,0,2),c(1,2,0),c(2,0,1),c(2,1,0))
  d$Score <- as.vector(t(cbind(permutations,permutations)))
  fit <- fit_mfrm_testlet(d,"Person","Score","Rater",c("Rater","Criterion"),0:2,quad_points=61,person_sd=1)
  expect_true(fit$checks$NumericalReady); expect_true(fit$checks$EstimatedVarianceBoundary)
  expect_equal(fit$calibration$variance,0); expect_lte(fit$checks$VarianceScore,0)
  expect_true(all(is.na(fit$calibration_table$SE)))
  scores <- predict(fit,d[d$Person=="P1",])
  expect_true(scores$settings$estimated_variance_boundary)
  expect_identical(scores$table$Status,"available_conditional")
  d$Score[d$Person=="P6"] <- NA
  local_mocked_bindings(mfrm_testlet_person_interval=function(...) stop("deliberately unresolved integration"))
  failed <- predict(fit,d,missing="omit")
  expect_equal(nrow(failed$table),6)
  expect_equal(failed$table$Status,c(rep("unavailable",5),"prior_only"))
  expect_true(all(is.na(failed$table$Lower[1:5])))
  expect_true(all(grepl("unresolved",failed$table$Reason[1:5])))
  grDevices::pdf(tempfile(fileext=".pdf")); old <- par("mar")
  expect_silent(plot(failed)); expect_equal(par("mar"),old); grDevices::dev.off()
  expect_identical(plot_data(plot(failed,draw=FALSE))$table$Person,failed$table$Person)
})

test_that("known-zero binary models without fixed facets score new memberships", {
  d <- expand.grid(Item=1:2,Block=c("a","b"),Person=c("p1","p2","p3"))
  d$Score <- rep(c(0,1),6)
  fit <- fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,
    testlet_variance=0,quad_points=61,person_sd=1)
  expect_true(fit$checks$NumericalReady); expect_true(fit$checks$InformationPositive)
  expect_false(fit$checks$EstimatedVarianceBoundary)
  expect_length(fit$runs,1); expect_equal(fit$calibration$steps,0,tolerance=1e-6)
  expect_true(is.finite(fit$calibration_table$SE))
  new <- d[1:4,]; new$Person <- "new"; new$Block <- "new_block"
  score <- predict(fit,new)
  expect_equal(score$table$Testlets,1)
  expect_identical(score$table$Status,"available_conditional")
  # At known zero dependence only the four Bernoulli likelihoods remain.
  density <- function(t) dnorm(t)*plogis(t)^2*plogis(-t)^2
  mass <- integrate(density,-Inf,Inf,rel.tol=1e-11)$value
  expect_equal(score$table$Estimate,0,tolerance=1e-7)
  expect_equal(integrate(density,-Inf,score$table$Lower,rel.tol=1e-11)$value/mass,.025,tolerance=1e-7)
  expect_error(plot(fit,draw=FALSE),"fixed facet")
})
