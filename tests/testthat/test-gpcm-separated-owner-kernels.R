test_that("separate slope and step indices preserve MML kernels and the PCM reduction", {
  # Independent kernel checks complement the fitted MML workflow tests.
  d <- expand.grid(Person=paste0('P',1:6),Rater=paste0('R',1:3),Criterion=paste0('C',1:2),stringsAsFactors=FALSE)
  d$Score <- (as.integer(sub('P','',d$Person))+2*as.integer(sub('R','',d$Rater))+as.integer(sub('C','',d$Criterion))) %% 3
  rows <- list()
  for (sparse in c(FALSE,TRUE)) for (step in c('Rater','Criterion')) {
    slope <- setdiff(c('Rater','Criterion'),step)
    data <- if(sparse) d[as.integer(sub('R','',d$Rater)) != 1+as.integer(sub('P','',d$Person))%%3,] else d
    prep <- mfrmr:::prepare_mfrm_data(data,'Person',c('Rater','Criterion'),'Score')
    build <- function(model) mfrmr:::build_estimation_config(prep,model,'MML',step,
      if(model=='GPCM') slope else NULL,weight_col=NULL,
      facet_signs=mfrmr:::build_facet_signs(prep$facet_names)$signs,positive_facets=character(),
      noncenter_facet='Person',dummy_facets=character(),anchor_df=NULL,group_anchor_df=NULL)
    cg <- build('GPCM'); cp <- build('PCM')
    idx <- mfrmr:::build_indices(prep,step,slope)
    p <- mfrmr:::build_initial_param_vector(cg$config,cg$sizes)+seq(-.15,.2,length.out=sum(unlist(cg$sizes)))
    params <- mfrmr:::expand_params(p,cg$sizes,cg$config)
    # Separate scalar category recursion, integrating each Person with adaptive R integration.
    base <- Reduce(`+`,lapply(prep$facet_names,function(f) cg$config$facet_signs[[f]]*params$facets[[f]][idx$facets[[f]]]))
    ref <- -sum(vapply(split(seq_along(idx$score_k),idx$person),function(observations) {
      integrand <- function(theta) vapply(theta,function(t) {
        lp <- vapply(observations,function(i) {
          adj <- params$slopes[idx$slope_idx[i]]*(t+base[i]-params$steps_mat[idx$step_idx[i],])
          logits <- c(0,cumsum(adj)); m <- max(logits)
          logits[idx$score_k[i]+1]-m-log(sum(exp(logits-m)))
        },numeric(1))
        exp(sum(lp))*dnorm(t)
      },numeric(1))
      ans <- integrate(integrand,-Inf,Inf,rel.tol=1e-10,stop.on.error=TRUE)
      log(ans$value)
    },numeric(1)))
    for(integration in c('fixed','adaptive')) {
      cg$config$estimation_control <- cp$config$estimation_control <- list(mml_integration=integration,quad_points=61L)
      evaluator <- function(c,indices) mfrmr:::make_mfrm_direct_evaluator('MML',
        mfrmr:::make_param_cache(c$sizes,c$config,indices,is_mml=TRUE),indices,c$config,c$sizes,mfrmr:::gauss_hermite_normal(61L))
      ev <- evaluator(cg,idx); pv <- evaluator(cp,mfrmr:::build_indices(prep,step,NULL))
      value_error <- abs(ev$value(p)-ref)
      numeric_grad <- vapply(seq_along(p),function(j) {
        hi <- lo <- p; hi[j] <- hi[j]+1e-5; lo[j] <- lo[j]-1e-5
        (ev$value(hi)-ev$value(lo))/2e-5
      },numeric(1))
      gradient_error <- max(abs(ev$gradient(p)-numeric_grad))
      pieces <- mfrmr:::split_params(p,cg$sizes); pieces$log_slopes[] <- 0
      unit <- unlist(pieces[names(cg$sizes)],use.names=FALSE)
      pcm <- unlist(pieces[names(cp$sizes)],use.names=FALSE)
      unit_error <- abs(ev$value(unit)-pv$value(pcm))
      gg <- mfrmr:::split_params(ev$gradient(unit),cg$sizes)
      unit_gradient_error <- max(abs(unlist(gg[names(cp$sizes)],use.names=FALSE)-pv$gradient(pcm)))
      rows[[length(rows)+1L]] <- data.frame(Sparse=sparse,StepOwner=step,SlopeOwner=slope,
        Integration=integration,Rows=nrow(data),SlopeParameters=cg$sizes$log_slopes,
        LikelihoodError=value_error,GradientError=gradient_error,PCMError=unit_error,
        PCMGradientError=unit_gradient_error)
    }
  }
  result <- do.call(rbind,rows)
  expect_lt(max(result$LikelihoodError), 1e-7)
  expect_lt(max(result$GradientError), 1e-6)
  expect_lt(max(result$PCMError), 1e-9)
  expect_lt(max(result$PCMGradientError), 1e-8)
})
