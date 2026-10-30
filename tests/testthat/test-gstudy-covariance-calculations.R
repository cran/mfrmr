test_that("covariance generalized pivots match the published Cholesky map", {
  S <- matrix(c(2,-.6,-.6,1),2)
  nu <- 9; M <- 120L
  withr::local_seed(816)
  observed_seed <- .Random.seed
  got <- .mfrm_mvgt_covariance_pivots(matrix(c(2,-.6,1),1),nu,M)
  assign('.Random.seed',observed_seed,envir=.GlobalEnv)
  q <- rchisq(M,nu); r <- rchisq(M,nu-1); z <- rnorm(M)
  C <- t(chol(nu*S))
  expected <- vapply(seq_len(M),function(i) {
    V <- matrix(c(sqrt(q[i]),z[i],0,sqrt(r[i])),2)
    G <- C %*% solve(V)
    tcrossprod(G)[1,2]
  },numeric(1))
  expect_equal(drop(got),expected,tolerance=1e-12)
  # The observed value of the generalized pivot recovers the true parameter.
  Sigma <- matrix(c(3,.4,.4,1),2)
  V <- solve(t(chol(Sigma)),C)
  G <- C %*% solve(V)
  expect_equal(tcrossprod(G),Sigma,tolerance=1e-12)
  expect_true(all(is.na(.mfrm_mvgt_covariance_pivots(matrix(c(1,.3,1),1),1,10))))
})

test_that("regression region agrees with independent lm and variance intervals", {
  x <- c(-2,1,3,0,-1,4)
  y <- c(1,-2,2,-3,2,1)
  nu <- length(x); level <- .95
  S <- crossprod(cbind(x,y))/nu
  expected_slope <- confint(lm(x~y-1),level=sqrt(level))
  tail <- (1-sqrt(level))/2
  expected_variance <- sum(y^2)/qchisq(c(1-tail,tail),nu)
  got <- .mfrm_mvgt_wishart_covariance_region(matrix(c(S[1,1],S[1,2],S[2,2]),1),nu,level)
  expect_equal(drop(got$bounds),range(outer(expected_slope,expected_variance)),tolerance=1e-12)
  swapped <- .mfrm_mvgt_wishart_covariance_region(matrix(c(S[2,2],S[1,2],S[1,1]),1),nu,level)
  expect_equal(swapped,got)
  for(scale in c(-100,.01,1e50)) {
    z <- matrix(c(S[1,1],scale*S[1,2],scale^2*S[2,2]),1)
    changed <- .mfrm_mvgt_wishart_covariance_region(z,nu,level)
    expect_equal(sort(drop(changed$bounds)/scale),drop(got$bounds),tolerance=1e-10)
  }
  one <- .mfrm_mvgt_wishart_covariance_region(matrix(c(4,2,1),1),1,level)
  expect_equal(drop(one$bounds),c(-2,2)/qchisq((1-level)/2,1))
  expect_match(one$method,"Cauchy")
})

test_that("covariance SE follows the independently solved Wishart moment identity", {
  nu <- c(2,5,20)
  S <- cbind(c(3,2,1),c(-.8,0,.6),c(2,5,1))
  A <- rbind(c(1,-1,0)/3,c(0,1,-1)/4,c(0,0,1))
  ans <- .mfrm_mvgt_covariance_bounds(S,A,nu,.95)
  # Solve the two second-moment equations instead of restating the code.
  v <- vapply(seq_along(nu),function(i) {
    n <- nu[i]
    moment_map <- rbind(c(1,2/n),c(1/n,1+1/n))
    estimated_products <- solve(moment_map,c(S[i,1]*S[i,3],S[i,2]^2))
    sum(estimated_products)/n
  },numeric(1))
  expect_equal(ans$se^2,drop(A^2%*%v),tolerance=1e-12)
  # The moment expression has the correct expectation at arbitrary covariances.
  for(n in c(2,3,20))for(cov in c(-.9,0,.9)) {
    a <- 2; b <- 1
    E <- c(a*b+2*cov^2/n,cov^2+(a*b+cov^2)/n)
    expect_equal((n*E[1]+(n-2)*E[2])/((n-1)*(n+2)),(a*b+cov^2)/n)
  }
  singular_df <- .mfrm_mvgt_covariance_bounds(S,A,c(1,5,20),.95)
  expect_true(is.na(singular_df$se[1]))
  expect_true(all(is.finite(singular_df$bounds)))
  expect_equal(singular_df$se[-1],ans$se[-1])
})

test_that("covariance SE remains representable when its square is not", {
  ms <- matrix(c(2,-.6,1),1)
  base <- .mfrm_mvgt_covariance_bounds(ms,matrix(1),9,.95)
  for (scale in c(1e-200,1e200)) {
    ci <- .mfrm_mvgt_covariance_bounds(ms*scale,matrix(1),9,.95)
    expect_equal(ci$se/scale,base$se,tolerance=1e-12)
    expect_equal(ci$bounds/scale,base$bounds,tolerance=1e-12)
    expect_identical(ci$se_status,base$se_status)
  }
})

test_that("component bounds preserve the shared signed linear map", {
  ms <- rbind(c(3,-.5,2),c(2,.3,4),c(1,-.6,1))
  A <- rbind(c(1,-1,0)/4,c(0,1,-1)/6,c(0,0,1))
  df <- c(9,5,45)
  ans <- .mfrm_mvgt_covariance_bounds(ms,A,df,.95)
  for(i in 1:3) {
    use <- which(A[i,]!=0); a <- A[i,use]
    z <- .mfrm_mvgt_wishart_covariance_region(ms[use,,drop=FALSE],df[use],.95^(1/length(use)))
    vertices <- expand.grid(rep(list(1:2),length(use)))
    values <- apply(vertices,1,function(indices)sum(a*z$bounds[cbind(seq_along(use),indices)]))
    expect_equal(drop(ans$bounds[i,]),range(values))
  }
  reversed <- ms; reversed[,2] <- -reversed[,2]
  neg <- .mfrm_mvgt_covariance_bounds(reversed,A,df,.95)
  expect_equal(neg$bounds,-ans$bounds[,2:1])
  expect_equal(neg$se,ans$se)
  expect_true(any(ans$bounds<0))
})

covariance_interval_fixture <- function() {
  tasks <- read.csv(system.file('extdata','mgenova-table12.csv',package='mfrmr',mustWork=TRUE))
  mfrm_g_study(tasks,c('V','W'),rater=NULL)
}

test_that("two-score intervals preserve targets, scaling, seeds and saved sources", {
  g <- covariance_interval_fixture(); original <- serialize(g,NULL)
  withr::local_seed(94)
  before <- .Random.seed
  for(method in c('gci','regression_region')) {
    ci <- confint(g,score=c('V','W'),assumption='normal',method=method)
    expect_identical(.Random.seed,before)
    expect_identical(attr(ci,'kind'),'Covariance')
    expect_equal(ci$Estimate,unname(vapply(g$components,function(m)m[1,2],numeric(1))))
    expect_identical(ci$Lower,ci$RawLower); expect_identical(ci$Upper,ci$RawUpper)
    reverse <- confint(g,score=c('W','V'),assumption='normal',method=method)
    for(n in c('Estimate','SE','Lower','Upper'))expect_equal(reverse[[n]],ci[[n]])
    saved <- g; saved$data <- NULL
    expect_identical(confint(saved,score=c('V','W'),assumption='normal',method=method),ci)
    for(scale in c(1e-60,1e60)) {
      changed <- g
      changed$components <- lapply(g$components,function(m)m*scale)
      changed$mean_products <- lapply(g$mean_products,function(m)m*scale)
      z <- confint(changed,score=c('V','W'),assumption='normal',method=method)
      for(n in c('Estimate','SE','Lower','Upper'))expect_equal(z[[n]]/scale,ci[[n]],tolerance=1e-10)
    }
    expect_equal(confint(g,parm='Residual',score=c('V','W'),assumption='normal',method=method)$Estimate,ci$Estimate[3])
    if(method=='gci') {
      band <- attr(ci,'monte_carlo')$endpoint_bands
      expect_true(all(band$LowerMCLower<=ci$Lower & ci$Lower<=band$LowerMCUpper))
      expect_true(all(band$UpperMCLower<=ci$Upper & ci$Upper<=band$UpperMCUpper))
    }
  }
  expect_identical(serialize(g,NULL),original)
})

test_that("df-one covariance SE and interval availability remain distinct", {
  tasks <- read.csv(system.file('extdata','mgenova-table12.csv',package='mfrmr',mustWork=TRUE))
  tasks <- tasks[tasks$Task %in% unique(tasks$Task)[1:2],]
  g <- mfrm_g_study(tasks,c('V','W'),rater=NULL)
  conservative <- confint(g,score=c('V','W'),assumption='normal',method='regression_region')
  generalized <- confint(g,score=c('V','W'),assumption='normal',method='gci')
  row <- conservative$Source=='Task'
  expect_true(is.na(conservative$SE[row]))
  expect_true(is.finite(conservative$Lower[row]) && is.finite(conservative$Upper[row]))
  expect_true(is.na(generalized$Lower[row]))
  expect_match(conservative$SEStatus[row],'df is one')
  expect_match(generalized$Status[row],'df one')
  expect_true(all(is.finite(generalized$Lower[!row])))
})

test_that("reciprocal score units cannot erase covariance information", {
  g <- covariance_interval_fixture()
  for (method in c("gci","regression_region")) {
    call <- function(x)confint(x,score=c("V","W"),assumption="normal",method=method)
    base <- call(g)
    for (exponent in c(-100,60,80,100)) {
      scale <- c(10^exponent,10^-exponent)
      changed <- g
      changed$components <- lapply(g$components,function(m)m*outer(scale,scale))
      changed$mean_products <- lapply(g$mean_products,function(m)m*outer(scale,scale))
      ci <- call(changed)
      for (column in c("Estimate","SE","Lower","Upper"))
        expect_equal(ci[[column]],base[[column]],tolerance=1e-11)
      expect_identical(ci$Status,base$Status)
      expect_identical(ci$SEStatus,base$SEStatus)
      expect_equal(attr(ci,"mean_products")$normalized,attr(base,"mean_products")$normalized,
        tolerance=1e-12)
      expect_equal(attr(ci,"mean_products")$score_scale/scale,
        attr(base,"mean_products")$score_scale,tolerance=1e-12,ignore_attr=TRUE)
      if (method=="gci")
        expect_equal(attr(ci,"monte_carlo")$endpoint_bands,
          attr(base,"monte_carlo")$endpoint_bands,tolerance=1e-11)
      # Scaling must not hide a materially inconsistent saved covariance.
      changed$components$Person[1,2] <- changed$components$Person[2,1] <- 1e4
      expect_error(call(changed),"disagree")
    }
  }
})

test_that("unresolved covariance intervals retain points and independent SEs", {
  g <- covariance_interval_fixture()
  for (method in c("regression_region","gci")) {
    call <- function(level)confint(g,score=c("V","W"),assumption="normal",
      method=method,level=level)
    base <- call(.95); ci <- call(1e-100)
    expect_identical(ci$Estimate,base$Estimate)
    expect_identical(ci$SE,base$SE)
    expect_true(all(is.na(ci$Lower) & is.na(ci$Upper)))
    expect_true(all(grepl("numeric resolution",ci$Status)))
  }
  unresolved <- confint(g,score=c("V","W"),assumption="normal",method="gci",
    draws=10001,level=1e-8)
  expect_true(all(is.na(unresolved$Lower) & is.na(unresolved$Upper)))
  expect_true(all(grepl("draw count",unresolved$Status)))
})

test_that("covariance input checks reject altered summaries and ignored settings", {
  g <- covariance_interval_fixture()
  call <- function(...)confint(g,score=c('V','W'),assumption='normal',...)
  expect_error(call(),'explicitly select')
  expect_error(call(method='regression_region',draws=5000),'only')
  expect_error(call(method='gci',draws=1000),'100 expected')
  expect_error(call(method='gci',seed=-1),'seed')
  expect_error(call(method=NA_character_),'method')
  expect_error(call(method='gci',parm=c('Task','Task')),'distinct')
  expect_error(confint(g,score=c('V','V'),method='gci',assumption='normal'),'distinct')
  expect_error(confint(g,score='V',method='gci',assumption='normal'),'Single-score')
  expect_error(confint(g,score='V',draws=8000,assumption='normal'),'Single-score')
  bad <- g; bad$mean_products$Person[1,2] <- bad$mean_products$Person[2,1] <- 1e4
  expect_error(confint(bad,score=c('V','W'),method='gci',assumption='normal'),'semidefinite')
  bad <- g; bad$mean_products$Residual[,] <- matrix(c(1e-100,2e-100,2e-100,1e-100),2)
  expect_error(confint(bad,score=c('V','W'),method='gci',assumption='normal'),'semidefinite')
  bad <- g; bad$components$Person[1,2] <- bad$components$Person[2,1] <- 1e4
  expect_error(confint(bad,score=c('V','W'),method='gci',assumption='normal'),'disagree')
  bad <- g; bad$design$method <- 'Multivariate MINQUE(0)'
  expect_error(confint(bad,score=c('V','W'),method='gci',assumption='normal'),'MINQUE')
})

test_that("covariance forest preserves signed intervals and saved conversion", {
  skip_if_not_installed('ggplot2')
  g <- covariance_interval_fixture()
  withr::local_pdf(NULL,width=7,height=5)
  for(method in c('gci','regression_region')) {
    ci <- confint(g,score=c('V','W'),assumption='normal',method=method)
    p <- plot(g,type='covariance',score=c('V','W'),intervals=TRUE,
      method=method,assumption='normal')
    expect_identical(p$data$kind,'Covariance')
    expect_equal(nrow(p$data$panels),1)
    expect_identical(p$data$intervals[[1]],ci)
    expect_identical(p$data$table$Lower,ci$Lower)
    expect_match(p$data$x_label,'Covariance')
    expect_false(any(grepl('negative population variances',p$data$notes$Text)))
    p <- unserialize(serialize(p,NULL))
    gg <- as_ggplot(p)
    expect_identical(plot_data(gg,'table'),p$data$table)
    expect_null(gg$labels$title);expect_null(gg$labels$caption)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
  expect_error(plot(g,type='covariance',score='V',intervals=TRUE,method='gci',assumption='normal'),'exactly two')
  expect_error(plot(g,method='gci'),'only when')
})

test_that("crossed and nested covariance outputs retain maps and score units", {
  fixture <- new.env()
  eval(parse(test_path('test-multivariate-gtheory.R'))[[1]],fixture)
  crossed <- mfrm_g_study(fixture$mvgt_fixture()$data,c('Content','Organization'))
  nested <- mfrm_g_study(mfrmr_example_gtheory_nested,c('Content','Organization'),
    facets=c(Occasion='Task',Assessor='Rater'),nesting=c(Assessor='Occasion'))
  for(g in list(crossed,nested))for(method in c('gci','regression_region')) {
    ci <- confint(g,score=c('Content','Organization'),method=method,assumption='normal')
    expect_identical(ci$Source,names(g$components))
    expect_equal(ci$Estimate,unname(vapply(g$components,function(m)m[1,2],numeric(1))))
    expect_true(all(is.finite(ci$Lower)))
    changed <- g
    rescale <- function(m) {
      w <- c(1,10)
      m*outer(w,w)
    }
    changed$components <- lapply(g$components,rescale)
    changed$mean_products <- lapply(g$mean_products,rescale)
    rescaled <- confint(changed,score=c('Content','Organization'),method=method,assumption='normal')
    for(n in c('Estimate','SE','Lower','Upper'))expect_equal(rescaled[[n]]/10,ci[[n]],tolerance=1e-10)
    expect_identical(attr(ci,'source_design')$nesting,g$design$nesting)
  }
})
