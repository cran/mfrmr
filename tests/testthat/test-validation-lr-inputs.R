lr_input_env <- new.env(parent=globalenv())
lr_path <- test_path('..','..','inst','validation','mfrm-wide-map-lr-inputs-20261002.R')
skip_if_not(file.exists(lr_path),'Repository-only D-LR generator')
sys.source(lr_path,lr_input_env)

lr_parent <- function(n=20L) {
  root <- normalizePath(test_path('..','..'))
  path <- file.path(root,'validation-results/mfrm-wide-map-r50-20261001/ab-inputs')
  skip_if_not(file.exists(file.path(path,'rng-registry.rds')),'Saved A/B inputs')
  r <- readRDS(file.path(path,'rng-registry.rds'))
  unit <- r$units[r$units$N==n & r$units$SD==1 & r$units$Truth=='S-PCM' & r$units$Replicate==1, ]
  unit$ABBundle <- file.path(path,unit$Bundle)
  unit$ABBundleMD5 <- unname(tools::md5sum(unit$ABBundle))
  unit$Bundle <- 'test.rds'
  seed <- parallel::nextRNGStream(r$next_unused_unit_stream)
  list(unit=unit,seeds=list(covariate=seed,
    future_covariate_first_unused=parallel::nextRNGSubStream(seed)))
}

test_that('D-LR covariate pairing preserves RNG, Persons and exact conditional targets', {
  withr::local_seed(32)
  before <- .Random.seed; kind <- RNGkind()
  p <- lr_parent(); x <- lr_input_env$dlr_generate(p$unit,p$seeds)
  expect_identical(.Random.seed,before); expect_identical(RNGkind(),kind)
  expect_identical(x,lr_input_env$dlr_generate(p$unit,p$seeds))
  expect_true(lr_input_env$dlr_check(x,p$unit,p$seeds))
  ab <- readRDS(p$unit$ABBundle)
  for (key in names(x$person_data)) {
    t <- x$population_truth[[key]]; X <- x$person_data[[key]]$X
    expect_identical(x$person_data[[key]]$Person,names(ab$theta))
    # Independent Gaussian quadratic-form identity verifies the conditional law.
    expect_equal(unname(ab$theta^2/p$unit$SD^2+x$z^2),
      unname(X^2+(ab$theta-t$coefficients['X']*X)^2/t$residual_variance),tolerance=1e-12)
    expect_equal(unname(t$residual_variance+t$coefficients['X']^2),t$marginal_variance)
  }
  expect_identical(x$person_data$rho0$X,x$z)
  expect_false(isTRUE(all.equal(mean(x$person_data[['rho0.5']]$X),0,tolerance=1e-8)))
  bad <- p$unit; bad$ABBundleMD5 <- 'changed'
  expect_error(lr_input_env$dlr_generate(bad,p$seeds))
  y <- x; y$person_data[['rho0.5']]$X[1] <- 0
  expect_error(lr_input_env$dlr_check(y,p$unit,p$seeds))
})

test_that('D-LR fitting receives only observed X and unchanged response rows', {
  p <- lr_parent(); x <- lr_input_env$dlr_generate(p$unit,p$seeds)
  out <- tempfile(); dir.create(out); withr::defer(unlink(out,recursive=TRUE))
  lr_input_env$dlr_save(x,file.path(out,p$unit$Bundle))
  row <- p$unit; row$View <- 'paired-L1'; row$CovariateView <- 'rho0.5'; row$Rho <- .5
  row$BundleMD5 <- unname(tools::md5sum(file.path(out,row$Bundle)))
  row$ControlConditionId <- 'AB:paired-L1:N20:SD1:S-PCM'
  row$ConditionId <- paste('DLR',row$CovariateView,row$ControlConditionId,sep=':')
  for(model in c('PCM','GPCM')) {
    a <- lr_input_env$dlr_fit_input(row,out,model)
    expect_identical(names(a$person_data),c('Person','X'))
    expect_identical(names(a$data),c('Person','Rater','Criterion','Score'))
    expect_identical(a$data,readRDS(row$ABBundle)$views[['paired-L1']][names(a$data)])
    expect_identical(a$person_data,x$person_data[['rho0.5']])
    expect_identical(all.vars(a$population_formula),'X')
    expect_null(a$quad_points)
  }
  a <- lr_input_env$dlr_fit_input(row,out,'PCM')
  row$View <- 'paired-L2'; row$ControlConditionId <- 'AB:paired-L2:N20:SD1:S-PCM'
  row$ConditionId <- paste('DLR',row$CovariateView,row$ControlConditionId,sep=':')
  b <- lr_input_env$dlr_fit_input(row,out,'PCM')
  expect_identical(a$person_data,b$person_data); expect_equal(nrow(b$data),2*nrow(a$data))
  row$Rho <- 0
  expect_error(lr_input_env$dlr_fit_input(row,out,'PCM'))
  expect_error(lr_input_env$dlr_fit_input(row,out,'RSM'))
  expect_identical(lr_input_env$dlr_models('S-RSM'),'RSM')
  expect_identical(lr_input_env$dlr_models('S-GPCM'),'GPCM')
})
