# Expected-information score with extended nuisance mean-square expectations.
# The public scope is balanced two-crossed-facet G only.
.mfrm_mvgt_score_problem <- function(ms, df, h) {
  stopifnot(length(ms) == length(df), length(h) == length(ms),
    all(is.finite(c(ms, df, h))), all(ms > 0), all(df > 0), any(h != 0))
  active <- which(h != 0)
  scale <- max(ms[active]); hs <- max(abs(h[active]))
  m <- ms[active]/scale; d <- df[active]; a <- abs(h[active])/hs
  groups <- list(which(h[active] > 0), which(h[active] < 0))
  if (any(lengths(groups) == 0L)) return(NULL)
  free <- lapply(groups, function(g) head(g, -1L))
  dimensions <- lengths(free); offsets <- c(0L, cumsum(dimensions))
  probabilities <- function(z) {
    p <- numeric(length(active))
    for (j in 1:2) {
      v <- if (dimensions[j]) c(z[seq.int(offsets[j]+1L, offsets[j+1L])], 0) else 0
      v <- exp(v-max(v)); p[groups[[j]]] <- v/sum(v)
    }
    p
  }
  encode <- function(p) unlist(lapply(groups, function(g)
    if (length(g) > 1L) log(p[head(g, -1L)]/p[tail(g, 1L)]) else numeric()), use.names=FALSE)
  evaluate <- function(z, gradient=FALSE) {
    p <- probabilities(z)
    q <- d*m*a/p; Q <- sum(q); D <- sum(d); tau <- Q/D
    if (gradient) return(unlist(lapply(seq_along(groups), function(j) {
      g <- groups[[j]]
      (d-D*q/Q-p*(sum(d[g])-D*sum(q[g])/Q))[free[[j]]]
    }), use.names=FALSE))
    # Subtract the unrestricted saturated criterion term by term. Large
    # irrelevant log-scale constants otherwise trigger premature stopping.
    u <- log(tau)+log(p)-log(a)-log(m)
    sum(d*(u+expm1(-u)))
  }
  decode <- function(z) {
    p <- probabilities(z); tau <- sum(d*m*a/p)/sum(d)
    p*tau/a
  }
  starts <- lapply(list(m*a, rep(1,length(a)), d), function(x) {
    for (g in groups) x[g] <- x[g]/sum(x[g])
    encode(x)
  })
  list(active=active, scale=scale, h=h[active]/hs, m=m, df=d,
    groups=groups, evaluate=evaluate, decode=decode, starts=starts,
    dimension=sum(dimensions))
}

.mfrm_mvgt_score_fit <- function(ms, df, h) {
  problem <- .mfrm_mvgt_score_problem(ms, df, h)
  if (is.null(problem)) return(list(valid=FALSE, statistic=NA_real_,
    status='Null has no positive mean-square expectation'))
  p <- problem
  if (p$dimension == 0L) {
    fits <- list(list(par=numeric(), objective=p$evaluate(numeric()), convergence=0L))
  } else {
    fits <- lapply(p$starts, function(start) tryCatch(
      stats::nlminb(start, p$evaluate, gradient=function(z) p$evaluate(z, TRUE),
        control=list(iter.max=300L, eval.max=600L, rel.tol=1e-12, x.tol=1e-10)),
      error=function(e) list(par=start, objective=Inf, convergence=99L, message=conditionMessage(e))))
  }
  # In information-standardized coordinates the nuisance score is the
  # projection orthogonal to h; this check is scale and null-coordinate invariant.
  checks <- lapply(fits, function(f) {
    lambda <- p$decode(f$par)
    z <- sqrt(p$df/2)*(p$m/lambda-1)
    v <- sqrt(2/p$df)*lambda*p$h
    multiplier <- sum(z*v)/sum(v^2)
    residual <- max(abs(z-multiplier*v))
    null_error <- abs(sum(p$h*lambda))/sum(abs(p$h*lambda))
    T <- sum(z*v)^2/sum(v^2)
    # The score construction needs a verified nuisance root, not an
    # optimizer's status label. Preserve the label separately below.
    valid <- is.finite(f$objective) && all(is.finite(c(lambda,T,residual,null_error))) &&
      all(lambda > 0) && residual < 1e-5 && null_error < 1e-10
    list(valid=valid, lambda=lambda, multiplier=multiplier, residual=residual,
      null_error=null_error, statistic=T)
  })
  eligible <- which(vapply(checks, `[[`, logical(1), 'valid'))
  if (!length(eligible)) eligible <- seq_along(fits)
  index <- eligible[which.min(vapply(fits[eligible], `[[`, numeric(1), 'objective'))]
  check <- checks[[index]]
  lambda <- check$lambda; multiplier <- check$multiplier
  residual <- check$residual; null_error <- check$null_error
  T <- check$statistic; valid <- check$valid
  full <- ms; full[p$active] <- lambda*p$scale
  list(valid=valid, statistic=if (valid) T else NA_real_, candidate_statistic=T,
    status=if (valid) 'Verified extended nuisance root' else 'Unresolved nuisance root',
    lambda=full, multiplier=multiplier, nuisance_residual=residual,
    null_error=null_error, active=p$active, asymptotic_bound=1/sqrt(min(p$df)),
    deviance=sum(df*(log(full)+ms/full)), optimizer_status=fits[[index]]$convergence,
    attempts=fits, root_checks=checks)
}

.mfrm_mvgt_score_target <- function(ms, map, error, rho) {
  stopifnot(is.finite(rho), rho >= 0, rho <= 1)
  h <- (1-rho)*map$A[1L, ]-rho*drop(error %*% map$A)
  .mfrm_mvgt_score_fit(ms, map$df, h)
}

# Scan twice using a shared cache; preserve failures, empty sets and pieces.
.mfrm_mvgt_score_set <- function(ms,map,error,level=.95) {
  cutoff<-stats::qchisq(level,1)
  tolerance<-min(1e-5,1e-5*cutoff)
  resolution<-64*.Machine$double.eps*max(1,cutoff)
  root_tolerance<-min(1e-9,1e-9*sqrt(cutoff))
  grid_tolerance<-min(1e-6,1e-6*sqrt(cutoff))
  a<-map$A[1,];delta<-drop(error%*%map$A)
  switches<-a/(a+delta);switches<-switches[is.finite(switches)&switches>0&switches<1]
  theta<-drop(map$A%*%ms);point<-theta[1]/sum((c(1,rep(0,length(theta)-1))+error)*theta)
  if (!is.finite(cutoff) || cutoff <= 0 || tolerance < resolution) {
    failed<-list(valid=FALSE,bounds=c(NA_real_,NA_real_),pieces=matrix(numeric(),0,2),
      roots=numeric(),cutoff_error=NA_real_,status="Score cutoff below numeric resolution")
    return(list(coarse=failed,fine=failed,grid_agreement=FALSE,profiles=list(),
      point=point,theta=theta,critical=cutoff,switches=switches,
      cutoff_tolerance=tolerance,numeric_resolution=resolution,
      root_tolerance=root_tolerance,grid_tolerance=grid_tolerance))
  }
  cache<-new.env(parent=emptyenv())
  evaluate<-function(rho) {
    key<-sprintf('%.17g',rho)
    if(!exists(key,cache,inherits=FALSE)) {
      fit<-.mfrm_mvgt_score_target(ms,map,error,rho)
      if(identical(fit$status,'Null has no positive mean-square expectation')) {
        fit$valid<-TRUE;fit$statistic<-Inf;fit$structural_rejection<-TRUE
      } else fit$structural_rejection<-FALSE
      fit$target<-rho;assign(key,fit,cache)
    }
    get(key,cache,inherits=FALSE)
  }
  statistic<-function(rho) {
    fit<-evaluate(rho)
    if(!fit$valid)stop('Unresolved nuisance root')
    fit$statistic-cutoff
  }
  invert<-function(nodes)tryCatch({
    grid<-sort(unique(c(seq(0,1,length.out=nodes),switches,
      if(is.finite(point)&&point>=0&&point<=1)point)))
    value<-vapply(grid,statistic,numeric(1))
    crossing<-which(sign(value[-length(value)])*sign(value[-1])<0)
    roots<-vapply(crossing,function(i)uniroot(statistic,c(grid[i],grid[i+1]),tol=root_tolerance)$root,numeric(1))
    ends<-sort(unique(c(0,roots,grid[value==0],1)))
    segments<-cbind(head(ends,-1),tail(ends,-1))
    accepted<-vapply(rowMeans(segments),function(rho)statistic(rho)<=0,logical(1))
    pieces<-segments[accepted,,drop=FALSE]
    cutoff_error<-if(length(roots))max(abs(vapply(roots,statistic,numeric(1)))) else 0
    bounds<-if(nrow(pieces))range(pieces) else c(NA_real_,NA_real_)
    valid<-is.finite(cutoff_error)&&cutoff_error<tolerance
    list(valid=valid,bounds=bounds,pieces=pieces,roots=roots,cutoff_error=cutoff_error,
      status=if(!valid)'Cutoff check failed' else if(!nrow(pieces))'Empty scanned set' else 'Computed scanned set')
  },error=function(e)list(valid=FALSE,bounds=c(NA_real_,NA_real_),pieces=matrix(numeric(),0,2),
    roots=numeric(),cutoff_error=NA_real_,status=conditionMessage(e)))
  coarse<-invert(65);fine<-invert(129)
  agree<-coarse$valid&&fine$valid&&identical(dim(coarse$pieces),dim(fine$pieces))&&
    (!length(fine$pieces)||max(abs(coarse$pieces-fine$pieces))<grid_tolerance)
  list(coarse=coarse,fine=fine,grid_agreement=agree,profiles=as.list(cache),
    point=point,theta=theta,critical=cutoff,switches=switches,
    cutoff_tolerance=tolerance,numeric_resolution=resolution,
    root_tolerance=root_tolerance,grid_tolerance=grid_tolerance)
}

.mfrm_mvgt_score_table <- function(d, tab, p, level) {
  rows <- diagnostics <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    counts <- as.numeric(d$design_grid[i, ])
    error <- c(0, 0, 0, 1/counts[1L], 1/counts[2L], 0, 1/prod(counts))
    duplicate <- which(vapply(diagnostics[seq_len(i-1L)], function(x)
      identical(x$error_coefficients, error), logical(1)))
    ans <- if (length(duplicate)) diagnostics[[duplicate[1L]]] else
      if (any(!is.finite(p$ms) | p$ms <= 0)) list(
        fine = list(valid = FALSE, bounds = c(NA_real_, NA_real_),
          pieces = matrix(numeric(), 0, 2), status = "Positive mean squares required"),
        grid_agreement = FALSE, profiles = list()) else
          .mfrm_mvgt_score_set(p$ms, p$design, error, level)
    ans$error_coefficients <- error
    resolved <- isTRUE(ans$fine$valid) && isTRUE(ans$grid_agreement)
    empty <- resolved && nrow(ans$fine$pieces) == 0L
    row <- tab[i, c("Scenario", names(d$design_grid), "Kind", "Score"), drop = FALSE]
    row$Metric <- "G"; row$Estimate <- tab$G[i]; row$PointStatus <- tab$GStatus[i]
    row$Lower <- if (resolved) ans$fine$bounds[1L] else NA_real_
    row$Upper <- if (resolved) ans$fine$bounds[2L] else NA_real_
    row$CandidateLower <- ans$fine$bounds[1L]; row$CandidateUpper <- ans$fine$bounds[2L]
    row$NumericalComplete <- resolved; row$Empty <- empty
    row$Level <- level; row$Method <- "Extended score (approximate)"
    row$Status <- if (empty) "Empty numerical score set" else if (resolved)
      if (nrow(ans$fine$pieces) == 1L) "Computed score approximation" else
        "Envelope of disconnected numerical score sets" else
      if (isTRUE(ans$fine$valid)) "Unresolved score-grid comparison" else ans$fine$status
    row$Validation <- "Experimental score calibration; finite-sample coverage not guaranteed"
    rows[[i]] <- row; diagnostics[[i]] <- ans
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "weights") <- p$weights; attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- d$gstudy$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- d$gstudy$data_usage$counts
  attr(out, "mean_squares") <- list(normalized = p$ms, scale = p$restore,
    df = p$design$df, component_map = p$design$A)
  attr(out, "diagnostics") <- diagnostics
  attr(out, "notes") <- c(
    "Explicit experimental score intervals for G from a balanced Person-by-two-crossed-random-facet source under independent normal effects and positive mean-square expectations. Select one score or prespecified composite. Not a Phi, SEM, nested, incomplete-source or plan-difference method; it is never an automatic default.",
    "Estimate and PointStatus retain the original ANOVA D-study. Nuisance component values used by the score may be negative provided the restricted marginal covariance remains positive definite; these computational values do not replace the reported components. No component is individually clipped.",
    "The score uses expected efficient information after solving the null nuisance equations, then an approximate chi-square(1) cutoff. Verified root residuals are separate from optimizer status labels and do not prove global optimality or asymptotic consistency. More future ratings do not add source information.",
    "Limited normal-model evaluation with 360 independent sources per condition gave true-value inclusion of 93.1-97.8% (Monte Carlo SE 0.8-1.3 percentage points). These are truth-constraint checks, not full endpoint coverage or a guarantee of 95% coverage. The method is not distribution robust; few source Persons can produce wide intervals or poor approximation.",
    "The inversion compares finite grids and checks cutoff crossings, retaining accepted pieces and their envelope. It cannot exclude narrower undetected pieces. NumericalComplete records this numerical check, not validated coverage. Empty identifies a resolved empty scanned set, not zero uncertainty or an estimated G of zero.",
    "An unresolved chi-square cutoff leaves missing bounds with NumericalComplete = FALSE and Empty = FALSE, retaining the original point and diagnostic cutoff/endpoint tolerances. This is arithmetic precision control, not a minimum sample size or evidence of an empty confidence set.",
    "An empty set or failed numerical check leaves Lower/Upper unavailable. CandidateLower/Upper and diagnostics retain any exploratory bounds; failed candidates are not accepted intervals. Available means that nonempty finite bounds were returned. No automatic fallback is applied.",
    "Intervals are pointwise for a fixed score/composite and future plan; their overlap does not test a plan difference. Save the complete result with saveRDS(). Existing forests and curve ribbons use saved bounds and preserve missing points without refitting.")
  .mfrm_mvgt_interval_result(out, d, p, "score")
}
