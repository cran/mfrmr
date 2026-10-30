.mfrm_mvgt_covariance_intervals <- function(object, score, level, method, draws, seed) {
  map <- .mfrm_mvgt_anova_summary(object)
  if (anyNA(score) || anyDuplicated(score) || !all(score %in% object$design$scores))
    stop("`score` must name two distinct original scores for covariance intervals.", call. = FALSE)
  marginal <- vapply(object$mean_products,function(m)diag(m[score,score]),numeric(2))
  if (any(marginal < 0))
    stop("Stored mean products are not positive semidefinite for the selected pair.", call. = FALSE)
  score_scale <- sqrt(apply(marginal,1L,max))
  score_scale[score_scale == 0] <- 1
  names(score_scale) <- score
  unit <- prod(score_scale)
  # A common variance scale can erase the smaller score's variance. Scale
  # each score separately; only covariance outputs use their product.
  ms <- t(vapply(object$mean_products,function(m)c(
    m[score[1],score[1]]/score_scale[1]/score_scale[1],
    m[score[1],score[2]]/score_scale[1]/score_scale[2],
    m[score[2],score[2]]/score_scale[2]/score_scale[2]),numeric(3)))
  estimate <- vapply(object$components,function(m)m[score[1],score[2]],numeric(1))
  tol <- 256*.Machine$double.eps
  if (any(!is.finite(ms)) ||
      any(abs(ms[,2]) > sqrt(ms[,1])*sqrt(ms[,3])*(1+tol)))
    stop("Stored mean products are not positive semidefinite for the selected pair.", call. = FALSE)
  if (any(abs(estimate/unit-drop(map$A%*%ms[,2])) > tol*pmax(1,rowSums(abs(map$A)))))
    stop("Stored covariance components disagree with their ANOVA mean products; recreate the G-study.", call. = FALSE)
  ans <- .mfrm_mvgt_covariance_bounds(ms,map$A,map$df,level,intervals=method!="gci")
  mc <- NULL
  if (method == "gci") {
    integer_scalar <- function(x)is.numeric(x) && !is.object(x) && !is.complex(x) &&
      is.null(dim(x)) && length(x)==1L && is.finite(x) && x==floor(x)
    tail <- (1-level)/2
    if (!integer_scalar(draws) || draws < 1000 || draws > .Machine$integer.max || draws*tail < 100)
      stop("`draws` must be an integer >= 1000 with at least 100 expected draws in each tail.",call.=FALSE)
    if (!integer_scalar(seed) || seed<0 || seed>.Machine$integer.max)
      stop("`seed` must be a nonnegative integer.",call.=FALSE)
    z <- with_preserved_rng_seed(seed,.mfrm_mvgt_covariance_pivots(ms,map$df,draws))
    ans$bounds[,] <- NA_real_
    bands <- matrix(NA_real_,nrow(map$A),4L)
    ans$status[] <- "Contributing mean-square df one or singular matrix"
    for(i in seq_len(nrow(map$A))) {
      use <- which(map$A[i,]!=0)
      if (any(!is.finite(z[,use,drop=FALSE]))) next
      target <- drop(z[,use,drop=FALSE]%*%map$A[i,use])
      if (any(!is.finite(target))) { ans$status[i] <- "Numerical generalized interval failure"; next }
      target <- sort(target)
      ans$bounds[i,] <- stats::quantile(target,c(tail,1-tail),type=1,names=FALSE)
      rank <- c(stats::qbinom(.0125,draws,tail),stats::qbinom(.9875,draws,tail)+1,
        stats::qbinom(.0125,draws,1-tail),stats::qbinom(.9875,draws,1-tail)+1)
      bands[i,] <- c(-Inf,target,Inf)[rank+1]
      ans$status[i] <- "Computed generalized interval"
      if (ans$bounds[i,1L] == ans$bounds[i,2L]) {
        ans$bounds[i,] <- NA_real_
        ans$status[i] <- if (tail == 1-tail)
          "Generalized quantile probabilities below numeric resolution" else
          "Generalized quantiles unresolved at this draw count"
      }
    }
    colnames(bands) <- c("LowerMCLower","LowerMCUpper","UpperMCLower","UpperMCUpper")
    mc <- list(draws=draws,seed=seed,rng_kind=RNGkind(),quantile_type=1L,
      band_level=.95,endpoint_bands=data.frame(Source=rownames(map$A),bands*unit))
  }
  bounds <- ans$bounds*unit; se <- ans$se*unit
  invalid <- apply(is.infinite(bounds) | (!is.na(ans$bounds) & ans$bounds != 0 & bounds == 0),1L,any)
  bounds[invalid,] <- NA_real_; ans$status[invalid] <- "Unrepresentable covariance endpoints"
  collapsed <- is.finite(bounds[,1L]) & is.finite(bounds[,2L]) & bounds[,1L] == bounds[,2L]
  bounds[collapsed,] <- NA_real_; ans$status[collapsed] <- "Covariance interval below numeric resolution"
  invalid_se <- !is.na(ans$se) & (!is.finite(se) | (ans$se > 0 & se == 0))
  ans$se_status[invalid_se] <- "Unrepresentable SE"
  se[invalid_se] <- NA_real_
  out <- data.frame(Source=rownames(map$A),Score=paste(score,collapse=" / "),
    Score1=score[1L],Score2=score[2L],Estimate=unname(estimate),SE=se,
    SEStatus=ans$se_status,Lower=bounds[,1L],Upper=bounds[,2L],
    RawLower=bounds[,1L],RawUpper=bounds[,2L],Level=level,
    Method=if(method=="gci")"Generalized covariance interval" else "Regression-region covariance enclosure",
    Status=ans$status,row.names=NULL)
  attr(out,"kind") <- "Covariance"; attr(out,"score_pair") <- score
  attr(out,"assumption") <- "normal"
  attr(out,"source_design") <- object$design[c("facets","counts","nesting","child_counts",
    "balanced","complete","method","score_convention")]
  attr(out,"source_rows") <- object$data_usage$counts
  attr(out,"mean_products") <- list(normalized=ms,scale=unit,score_scale=score_scale,
    df=map$df,component_map=map$A)
  attr(out,"diagnostics") <- if(method=="regression_region")ans$blocks else NULL
  attr(out,"monte_carlo") <- mc
  attr(out,"notes") <- c(
    "Signed covariance components for two prespecified scores under balanced independent normal effects. Units are products of the score units, not correlations. Negative covariance can be legitimate; neither estimates nor endpoints are clipped or repaired.",
    "SE is sampling uncertainty in a covariance component. Its square is unbiased when all contributing mean-square df exceed one. SEStatus and Status describe separate outputs; an unavailable SE does not remove an available interval.",
    "Each score is normalized separately before calculation. mean_products$score_scale retains those divisors; scale is their product for covariance outputs. Unrepresentable or collapsed intervals remain unavailable with an explicit Status, while estimates and calculable SEs are retained.",
    if(method=="gci") paste("The generalized interval uses Krishnamoorthy and Xia's Wishart-Cholesky pivot, then our independent-block ANOVA component map. It has no exact frequentist-coverage guarantee and requires df > 1 and nonsingular mean-product pairs for every contributing block.",
      "The monte_carlo attribute retains 95% joint numerical bands for the two endpoints of each component interval. These describe quantile error, not sampling uncertainty; more draws do not add source information.") else
      paste("The regression-region enclosure has at least nominal pointwise coverage under the positive-marginal-variance normal model. Independent t and chi-square pivots protect each mean-product block; one-df or singular blocks use marginal variance upper limits and Cauchy-Schwarz.",
        "It can be extremely wide with few facet levels. This includes conservatism from enclosing nuisance-parameter regions, not just limited source information. It is not a tight interval or a guarantee of useful precision."),
    "These are not simultaneous across components or score pairs, not distribution robust, and not intervals for an incomplete MINQUE source. Save the complete result with saveRDS() to retain assumptions and diagnostics.")
  out
}

# Signed covariance bounds from independent Wishart blocks.
.mfrm_mvgt_covariance_bounds <- function(ms, A, df, level, intervals = TRUE) {
  raw <- matrix(NA_real_, nrow(A), 2L)
  se <- rep(NA_real_, nrow(A))
  status <- rep("Unavailable mean-product region", nrow(A))
  se_status <- rep("Unavailable: contributing mean-square df is one", nrow(A))
  blocks <- vector("list", nrow(A))
  for (i in seq_len(nrow(A))) {
    use <- which(A[i, ] != 0)
    a <- A[i, use]; nu <- df[use]
    if (all(nu > 1)) {
      z <- ms[use, , drop = FALSE]
      # Sum squared SD contributions after scaling, without first squaring
      # very small/large covariances or multiplying marginal variances.
      product <- abs(a)*(sqrt(z[,1L])*sqrt(z[,3L]))
      covariance <- abs(a*z[,2L])
      scale <- max(product,covariance)
      se[i] <- if (scale == 0) 0 else sqrt(sum(
        (product/scale)^2*(nu/(nu-1)/(nu+2)) +
        (covariance/scale)^2*((nu-2)/(nu-1)/(nu+2))))*scale
      if (!is.finite(se[i])) se[i] <- NA_real_
      se_status[i] <- if (is.finite(se[i])) "Unbiased sampling-variance estimate" else "Numerical SE failure"
    }
    if (!intervals) next
    region <- .mfrm_mvgt_wishart_covariance_region(ms[use, , drop = FALSE], nu,
      exp(log(level)/length(use)))
    blocks[[i]] <- region
    lower <- region$bounds[, 1L]; upper <- region$bounds[, 2L]
    raw[i, ] <- c(sum(a*ifelse(a > 0, lower, upper)),
                  sum(a*ifelse(a > 0, upper, lower)))
    if (all(is.finite(raw[i, ]))) {
      if (raw[i, 1L] < raw[i, 2L]) status[i] <- "Conservative normal-model interval"
      else {
        raw[i, ] <- NA_real_
        status[i] <- "Covariance interval below numeric resolution"
      }
    } else {
      unavailable <- !apply(is.finite(region$bounds),1L,all)
      if (any(unavailable))
        status[i] <- paste("Unavailable:",paste(unique(region$method[unavailable]),collapse="; "))
    }
  }
  list(bounds=raw,se=se,status=status,se_status=se_status,blocks=blocks)
}

.mfrm_mvgt_wishart_covariance_region <- function(ms, df, level) {
  bound <- matrix(NA_real_, nrow(ms), 2L)
  method <- rep("Regression t and chi-square", nrow(ms))
  # T for the regression slope is independent of the predictor sum of squares.
  tail <- -expm1(log(level)/2)/2
  for (j in seq_len(nrow(ms))) {
    z <- ms[j, ]; nu <- df[j]
    product_sd <- sqrt(z[1L])*sqrt(z[3L])
    root_det <- sqrt(max(0, product_sd-abs(z[2L])))*sqrt(product_sd+abs(z[2L]))
    if (product_sd <= 0) { method[j] <- "Degenerate marginal mean squares"; next }
    if (nu == 1 || root_det == 0) {
      # Only upper marginal limits are needed: split the block error in two.
      t <- (1-level)/2
      upper <- product_sd*nu/stats::qchisq(t, nu)
      bound[j, ] <- c(-upper, upper)
      method[j] <- "Marginal variances and Cauchy-Schwarz"
    } else {
      width <- stats::qt(tail, nu-1, lower.tail=FALSE)*root_det/sqrt(nu-1)
      slope <- z[2L] + c(-width, width)
      variance_ratio <- nu/c(stats::qchisq(tail,nu,lower.tail=FALSE),stats::qchisq(tail,nu))
      bound[j, ] <- range(outer(slope, variance_ratio))
    }
    if (any(!is.finite(bound[j, ]))) {
      bound[j, ] <- NA_real_; method[j] <- "Unrepresentable covariance region"
    } else if (bound[j,1L] == bound[j,2L]) {
      bound[j, ] <- NA_real_; method[j] <- "Covariance region below numeric resolution"
    }
  }
  list(bounds=bound, method=method, level=level)
}

# Off-diagonal entry of (chol(nu*S)' V^-1)(chol(nu*S)' V^-1)'.
# Krishnamoorthy & Xia (2007), equation (4); our independent-block ANOVA map.
.mfrm_mvgt_covariance_pivots <- function(ms, df, draws) {
  out <- matrix(NA_real_, draws, nrow(ms))
  for (j in which(df > 1)) {
    z <- ms[j, ]; product_sd <- sqrt(z[1L])*sqrt(z[3L])
    root_det <- sqrt(max(0,product_sd-abs(z[2L])))*sqrt(product_sd+abs(z[2L]))
    if (product_sd <= 0 || root_det <= 0) next
    q <- stats::rchisq(draws,df[j])
    r <- stats::rchisq(draws,df[j]-1)
    out[,j] <- df[j]*(z[2L]-root_det*stats::rnorm(draws)/sqrt(r))/q
  }
  out
}
