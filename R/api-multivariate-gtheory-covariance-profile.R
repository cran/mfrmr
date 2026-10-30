# Joint normal REML for a signed covariance component; original estimates are retained.
.mfrm_mvgt_joint_groups <- function(g) {
  ids <- c(g$design$person, unname(g$design$facets))
  groups <- lapply(g$data[ids], factor)
  subsets <- if (length(ids) == 2L) list(1L, 2L, 1:2) else
    list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
  if (!is.null(g$design$nesting)) {
    labels <- c('Person', names(g$design$facets))
    child <- match(names(g$design$nesting), labels)
    parent <- match(unname(g$design$nesting), labels)
    subsets <- list(1L, parent, c(child, parent), c(1L, parent), 1:3)
  }
  setNames(lapply(subsets, function(s) .mfrm_mvgt_group(groups[s])),
    names(g$components))
}

.mfrm_mvgt_joint_sparse <- function(y, groups) {
  n <- nrow(y); k <- length(groups); df <- 2*n-2
  stopifnot(ncol(y) == 2L, all(is.finite(y)), length(unique(groups[[k]])) == n)
  d <- data.frame(response = as.vector(y), trait = factor(rep(c('A','B'),each=n)))
  terms <- paste0('group',seq_len(k-1L))
  for (j in seq_along(terms)) d[[terms[j]]] <- factor(rep(groups[[j]],2))
  formula <- stats::reformulate(c('0+trait', paste0('(0+trait|',terms,')')),
    response='response')
  parsed <- lme4::lFormula(formula,data=d,REML=TRUE)
  devfun <- do.call(lme4::mkLmerDevfun,parsed)
  env <- environment(devfun)
  order <- match(names(parsed$reTrms$cnms), terms)
  stopifnot(!anyNA(order), all(lengths(parsed$reTrms$cnms)==2L),
    length(parsed$reTrms$theta)==3L*(k-1L))
  function(factors, details=FALSE) {
    residual <- factors[[k]]
    stopifnot(length(factors)==k, all(diag(residual)>0))
    yw <- t(forwardsolve(residual,t(y)))
    theta <- unlist(lapply(order,function(j) {
      b <- forwardsolve(residual,factors[[j]])
      b[c(1,2,4)]
    }), use.names=FALSE)
    env$resp$setResp(as.vector(yw))
    profiled <- devfun(theta)
    sigma2 <- (env$pp$sqrL(1)+env$resp$wrss())/df
    logratio <- -log(sigma2)
    fixed <- profiled+df*(logratio+expm1(-logratio))
    correction <- (n-1)*2*sum(log(diag(residual)))
    if (details) return(list(objective=fixed+correction,
      whitened=fixed, correction=correction, sigma2=sigma2))
    fixed+correction
  }
}

.mfrm_mvgt_joint_coordinates <- function(k, target=NULL, tau=NULL) {
  sizes <- rep(3L,k); if (!is.null(target)) sizes[target] <- 2L
  ends <- cumsum(sizes); starts <- ends-sizes+1L
  decode <- function(par) lapply(seq_len(k),function(j) {
    z <- par[starts[j]:ends[j]]
    if (identical(j,target)) {
      if (tau==0) {
        if(j==k)z<-exp(z)
        return(.mfrm_mvgt_joint_factor(z[1],0,z[2]))
      }
      a <- exp(z[1]); c <- if(j==k)exp(z[2]) else z[2]
      return(.mfrm_mvgt_joint_factor(a,tau/a,c))
    }
    if(j==k)z[c(1,3)]<-exp(z[c(1,3)])
    .mfrm_mvgt_joint_factor(z[1],z[2],z[3])
  })
  encode <- function(factors) unlist(lapply(seq_len(k),function(j) {
    z <- factors[[j]][c(1,2,4)]
    if (identical(j,target)) {
      if(tau==0) {
        z <- c(z[1],sqrt(z[2]^2+z[3]^2))
        return(if(j==k)log(pmax(z,1e-8)) else z)
      }
      a <- max(z[1],sqrt(abs(tau)),1e-4)
      c <- sqrt(max(z[2]^2+z[3]^2-(tau/a)^2,.01))
      return(c(log(a),if(j==k)log(c) else c))
    }
    if(j==k)z[c(1,3)]<-log(pmax(z[c(1,3)],1e-8))
    z
  }),use.names=FALSE)
  list(lower=rep(-Inf,sum(sizes)),decode=decode,encode=encode)
}

.mfrm_mvgt_joint_cone_check <- function(factors, dense, target=NULL) {
  ans <- dense(factors,gradient=TRUE); k <- length(factors)
  cone <- lapply(seq_len(k),function(j) {
    sigma <- tcrossprod(factors[[j]]); g <- ans$gradient[[j]]
    lambda <- 0
    if(identical(j,target)) {
      h <- matrix(c(0,.5,.5,0),2); gs <- g%*%sigma; hs <- h%*%sigma
      lambda <- if(sum(hs*hs)>1e-16)-sum(gs*hs)/sum(hs*hs) else -2*g[1,2]
      g <- g+lambda*h
    }
    mineigen <- min(eigen(g,symmetric=TRUE,only.values=TRUE)$values)
    complementarity <- max(abs(g%*%sigma))
    # Positive residual is interior; random component matrices may be singular.
    residual_gradient <- if(j==k)max(abs(g)) else 0
    data.frame(Component=j,MinGradientEigenvalue=mineigen,
      Complementarity=complementarity,ResidualGradient=residual_gradient,Lambda=lambda,
      Rank=sum(eigen(sigma,symmetric=TRUE,only.values=TRUE)$values>1e-8),
      Valid=mineigen>=-.002 && complementarity<.002 && residual_gradient<.002)
  })
  list(objective=ans$objective,table=do.call(rbind,cone))
}

# Covariance derivatives using only PSD-preserving rank-one additions.
.mfrm_mvgt_joint_psd_factor <- function(m) {
  if(m[1,1]==0)return(.mfrm_mvgt_joint_factor(0,0,sqrt(m[2,2])))
  a<-sqrt(m[1,1]);b<-m[2,1]/a;remainder<-m[2,2]-b*b
  if(remainder < -128*.Machine$double.eps*max(1,abs(m)))stop('Invalid PSD factor')
  .mfrm_mvgt_joint_factor(a,b,sqrt(max(0,remainder)))
}
.mfrm_mvgt_joint_cone_gradient <- function(factors,objective) {
  value <- objective(factors);error <- consistency <- 0
  directions <- list(c(1,0),c(0,1),c(1,1)/sqrt(2),c(1,-1)/sqrt(2))
  gradients <- lapply(seq_along(factors),function(j) {
    m<-tcrossprod(factors[[j]])
    derivatives <- vapply(directions,function(u) {
      h<-1e-5*max(1,abs(m));z<-tcrossprod(u)
      changed <- function(delta) {
        f<-factors;f[[j]]<-.mfrm_mvgt_joint_psd_factor(m+delta*z);objective(f)
      }
      one<-changed(h);two<-changed(2*h)
      coarse<-(-3*value+4*one-two)/(2*h)
      best<-c(NA_real_,Inf)
      for(refinement in 0:4) {
        half<-changed(h/2);fine<-(-3*value+4*half-one)/h
        e<-abs(fine-coarse)/3
        if(is.finite(e) && e<best[2])best<-c((4*fine-coarse)/3,e)
        if(is.finite(e) && e<.0005)break
        h<-h/2;one<-half;coarse<-fine
      }
      error<<-max(error,best[2]);best[1]
    },numeric(1))
    consistency<<-max(consistency,abs(sum(derivatives[1:2])-sum(derivatives[3:4])))
    off<-(derivatives[3]-derivatives[4])/2
    matrix(c(derivatives[1],off,off,derivatives[2]),2)
  })
  list(objective=value,gradient=gradients,error=error,consistency=consistency,
    valid=all(is.finite(unlist(gradients))) && error<.0005 && consistency<.002)
}

.mfrm_mvgt_joint_factor <- function(a,b,c) matrix(c(a,b,0,c),2)

.mfrm_mvgt_joint_source <- function(g, score) {
  for (package in c("lme4", "minqa")) if (!requireNamespace(package, quietly = TRUE))
    stop("`method = 'reml_profile'` requires the optional ", package, " package.", call. = FALSE)
  ids <- c(g$design$person, unname(g$design$facets))
  if (!is.data.frame(g$data) || !all(c(ids, score) %in% names(g$data)) ||
      nrow(g$data) != g$data_usage$counts[["UsedRows"]] ||
      !all(vapply(g$data[score], is.numeric, logical(1))))
    stop("REML covariance profiles require the original used ratings; recreate the G-study.", call. = FALSE)
  y <- as.matrix(g$data[score])
  if (any(!is.finite(y))) stop("Invalid stored source ratings.", call. = FALSE)
  unit <- apply(abs(y), 2L, max); unit[unit == 0] <- 1
  y <- sweep(y, 2L, unit, "/"); y <- sweep(y, 2L, colMeans(y), "-")
  scale <- sqrt(colSums(y*y)/(nrow(y)-1))
  ans <- list(score = score, score_scale = unit*scale, valid = FALSE,
    status = "Degenerate selected score pair")
  if (any(scale <= 0 | !is.finite(scale))) return(ans)
  y <- sweep(y, 2L, scale, "/")
  singular_values <- svd(y, nu = 0L, nv = 0L)$d
  if (min(singular_values) <= 64*.Machine$double.eps*max(singular_values)) return(ans)
  groups <- .mfrm_mvgt_joint_groups(g)
  estimate <- lapply(g$components, function(m)
    sweep(sweep(m[score, score, drop = FALSE], 1L, ans$score_scale, "/"),
      2L, ans$score_scale, "/"))
  reconstructed <- .mfrm_mvgt_minque0(y, groups)$components
  if (any(!is.finite(unlist(estimate))) ||
      max(abs(unlist(estimate)-unlist(reconstructed))) > 1e-7*max(1, abs(unlist(estimate))))
    stop("Stored covariance components disagree with their used ratings; recreate the G-study.", call. = FALSE)
  initial <- lapply(estimate, function(m) {
    e <- eigen(m, symmetric = TRUE)
    t(chol(e$vectors %*% (pmax(e$values, .02)*t(e$vectors))))
  })
  ans$y <- y; ans$groups <- groups; ans$initial <- initial
  ans$objective <- .mfrm_mvgt_joint_sparse(y, groups)
  ans$fit <- .mfrm_mvgt_joint_fit(ans)
  ans$valid <- ans$fit$valid
  ans$status <- if (ans$valid) "Numerically checked joint REML source" else "Unresolved joint REML source optimum"
  ans
}

.mfrm_mvgt_joint_fit <- function(source, target = NULL, tau = NULL, start = NULL) {
  k <- length(source$groups)
  coordinates <- .mfrm_mvgt_joint_coordinates(k, target, tau)
  objective <- function(par) {
    f <- coordinates$decode(par)
    if (any(!is.finite(unlist(f))) || any(diag(f[[k]]) <= 0)) return(Inf)
    source$objective(f)
  }
  first <- if (is.null(start)) source$initial else start$factors
  starts <- lapply(list(first, rep(list(diag(sqrt(.2), 2)), k),
    rep(list(diag(sqrt(.7), 2)), k)), coordinates$encode)
  fit <- .mfrm_mvgt_profile_minimize(starts, objective, lower = coordinates$lower)
  if (!is.finite(fit$best$objective))
    return(list(fit = fit$best, attempts = fit$attempts, valid = FALSE,
      target = target, tau = tau, lr = NA_real_))
  inspect <- function(candidate) {
    factors <- coordinates$decode(candidate$par)
    checks <- tryCatch({
      derivative <- .mfrm_mvgt_joint_cone_gradient(factors, source$objective)
      cone <- if (derivative$valid) .mfrm_mvgt_joint_cone_check(factors, function(...) derivative, target)
        else list(table = data.frame(Component = seq_len(k), Valid = FALSE))
      gradient <- .mfrm_mvgt_profile_gradient(candidate$par, objective,
        lower = coordinates$lower, refinements = 4L)
      list(derivative = derivative, cone = cone, gradient = gradient, error = NULL)
    }, error = function(e) list(derivative = list(valid = FALSE),
      cone = list(table = data.frame(Component = seq_len(k), Valid = FALSE)),
      gradient = list(valid = FALSE), error = conditionMessage(e)))
    covariance <- vapply(factors, function(l) tcrossprod(l)[1L, 2L], numeric(1))
    target_error <- if (is.null(target)) 0 else abs(covariance[target]-tau)
    list(factors = factors, covariance = covariance, fit = candidate,
      gradient = checks$gradient, matrix_derivative = checks$derivative, cone = checks$cone,
      check_error = checks$error, target = target, tau = tau, target_error = target_error,
      valid = candidate$convergence == 0 && checks$gradient$valid && checks$derivative$valid &&
        all(checks$cone$table$Valid) && target_error < (if (is.null(target)) 1e-10 else 1e-10*max(1, abs(tau))))
  }
  original <- which(vapply(fit$attempts, identical, logical(1), fit$best))[1L]
  selected <- original; ans <- inspect(fit$best)
  checked <- list(); checked[[as.character(original)]] <- ans
  values <- vapply(fit$attempts, `[[`, numeric(1), "objective")
  minimum <- min(values[is.finite(values)])
  # Reuse the optimizer's existing deviance-equivalence tolerance, not a worse fit.
  if (!ans$valid) for (i in order(values)) {
    candidate <- fit$attempts[[i]]
    if (identical(i, original) || !is.finite(values[i]) ||
        values[i] > minimum+1e-7 || candidate$convergence != 0) next
    checked[[as.character(i)]] <- inspect(candidate)
    if (checked[[as.character(i)]]$valid) {
      ans <- checked[[as.character(i)]]; selected <- i; break
    }
  }
  ans$attempts <- fit$attempts
  ans$selection <- list(original = original, selected = selected,
    minimum_objective = minimum, objective_tolerance = 1e-7, checks = checked)
  ans
}

.mfrm_mvgt_joint_set <- function(source, target, level) {
  ans <- list(bounds = c(NA_real_, NA_real_), candidate_bounds = c(NA_real_, NA_real_),
    center = NA_real_, status = source$status, target = target)
  if (is.null(source$fit$covariance)) return(ans)
  ans$center <- source$fit$covariance[target]
  if (!source$valid) return(ans)
  sigma <- tcrossprod(source$fit$factors[[target]])
  anchor <- max(sqrt(prod(diag(sigma))), .05)
  cutoff <- .mfrm_mvgt_profile_cutoff(level, source$fit$fit$objective)
  ans$cutoff <- cutoff; ans$critical <- critical <- cutoff$critical
  if (!cutoff$valid) {
    ans$status <- "Likelihood-ratio cutoff below numeric resolution"
    return(ans)
  }
  cache <- list()
  at <- function(tau) {
    key <- sprintf("%.17g", tau)
    if (!is.null(cache[[key]])) return(cache[[key]])
    eligible <- Filter(function(x) isTRUE(x$valid) && !is.null(x$factors), cache)
    start <- if (length(eligible)) eligible[[which.min(vapply(eligible,
      function(x) abs(x$tau-tau), numeric(1)))]] else source$fit
    f <- .mfrm_mvgt_joint_fit(source, target, tau, start)
    f$lr <- f$fit$objective-source$fit$fit$objective
    f$valid <- f$valid && is.finite(f$lr) && f$lr >= -1e-6
    cache[[key]] <<- f
    f
  }
  status <- "Numerically checked REML covariance profile (approximate)"
  center <- at(ans$center)
  if (!center$valid || center$lr > min(1e-5, critical)) status <- "Unresolved constrained center"
  else for (side in 1:2) {
    previous <- ans$center
    for (step in 1:80) {
      next_target <- ans$center+c(-1, 1)[side]*anchor*sinh(step*.25)
      f <- at(next_target)
      if (!f$valid) { status <- "Unresolved required covariance profile"; break }
      if (f$lr > critical) {
        root <- tryCatch(stats::uniroot(function(tau) {
          z <- at(tau)
          if (!z$valid) stop("Unresolved root profile")
          z$lr-critical
        }, sort(c(previous, next_target)), tol = 1e-8, check.conv = TRUE)$root,
        error = function(e) NA_real_)
        ans$candidate_bounds[side] <- root
        if (is.finite(root) && at(root)$valid && abs(at(root)$lr-critical) < cutoff$tolerance)
          ans$bounds[side] <- root else status <- "Unresolved covariance endpoint"
        break
      }
      previous <- next_target
    }
    if (!is.finite(ans$bounds[side]) && grepl("^Numerically", status)) status <- "Unclosed covariance tail"
  }
  if (all(is.finite(ans$bounds)) && ans$bounds[1L] >= ans$bounds[2L]) {
    ans$bounds[] <- NA_real_
    status <- "Covariance profile interval below numeric resolution"
  }
  if (!all(is.finite(ans$bounds))) ans$bounds[] <- NA_real_
  ans$status <- status; ans$profiles <- cache; ans$critical <- critical; ans
}

.mfrm_mvgt_covariance_profile <- function(g, parm, score, level) {
  scores <- g$design$scores; sources <- names(g$components)
  if (!inherits(g, "mfrm_multivariate_gstudy") || !is.character(scores) ||
      !length(scores) || anyNA(scores) || anyDuplicated(scores) ||
      !is.list(g$components) || !length(sources) || anyNA(sources) || anyDuplicated(sources) ||
      !all(vapply(g$components, function(m) is.matrix(m) && is.numeric(m) &&
        !is.complex(m) && identical(dimnames(m), list(scores, scores)) &&
        all(is.finite(m)) && identical(m, t(m)), logical(1))))
    stop("The G-study covariance matrices or score identities are incomplete or altered.", call. = FALSE)
  if (!is.character(score) || !is.null(dim(score)) || length(score) != 2L ||
      anyNA(score) || anyDuplicated(score) || !all(score %in% scores))
    stop("`score` must name two distinct original scores for covariance intervals.", call. = FALSE)
  factors <- names(g$design$counts)
  if (!identical(factors[1L], "Person") || !.mfrm_mvgt_valid_facets(factors[-1L]) ||
      !identical(sources, .mfrm_mvgt_sources(factors, .mfrm_mvgt_nesting(g$design$nesting, factors[-1L]))))
    stop("The G-study component identities do not match its design; recreate the result.", call. = FALSE)
  selected <- .mfrm_mvgt_component_selection(parm, sources)
  source <- .mfrm_mvgt_joint_source(g, score)
  restore <- function(z) (z*source$score_scale[1L])*source$score_scale[2L]
  rows <- diagnostics <- vector("list", length(selected))
  for (j in seq_along(selected)) {
    i <- selected[j]; ans <- .mfrm_mvgt_joint_set(source, i, level)
    reported <- .mfrm_mvgt_profile_restore(ans, restore)
    values <- reported$values
    rows[[j]] <- data.frame(Source = sources[i], Score = paste(score, collapse = " / "),
      Score1 = score[1L], Score2 = score[2L], Estimate = g$components[[i]][score[1L], score[2L]],
      IntervalCenter = values[1L], SE = NA_real_, SEStatus = "Not supplied by profile inversion",
      Lower = values[2L], Upper = values[3L], RawLower = NA_real_, RawUpper = NA_real_,
      CandidateLower = values[4L], CandidateUpper = values[5L], Level = level,
      Method = "Joint REML covariance profile (approximate)", MethodCode = "reml_profile",
      Status = reported$status, Available = reported$available,
      Validation = "Experimental chi-square calibration; boundary and finite-sample coverage not guaranteed")
    diagnostics[[j]] <- ans
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL; names(diagnostics) <- sources[selected]
  source$objective <- NULL
  attr(out, "reml_source") <- source; attr(out, "diagnostics") <- diagnostics
  attr(out, "kind") <- "Covariance"; attr(out, "score_pair") <- score
  attr(out, "assumption") <- "normal"; attr(out, "source_data") <- g$data
  attr(out, "source_summary") <- g[c("design", "components", "mean_products", "degrees_of_freedom")]
  attr(out, "source_design") <- g$design[c("facets", "counts", "nesting", "child_counts",
    "balanced", "complete", "method", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "notes") <- c(
    "Explicit experimental signed covariance intervals for two prespecified scores under independent bivariate normal random effects and a common mean per score. Each component matrix is positive semidefinite; the residual matrix is positive definite. Source omissions must be noninformative conditional on the declared assignment.",
    "Estimate retains the original ANOVA/MINQUE covariance. IntervalCenter is a separate joint REML estimate of the same population covariance; the original point may lie outside its interval. Neither the original component matrices nor signed endpoints are clipped or repaired. Units are products of the two score units, not correlations.",
    "SE is unavailable: profile inversion does not compute a sampling SE for ANOVA, MINQUE or REML estimates. RawLower/RawUpper are inapplicable here and remain NA; CandidateLower/Upper hold numerical candidates separately.",
    "Fix the signed covariance, optimize all other variances and covariances in the joint restricted likelihood, and invert a chi-square(1) cutoff. Zero covariance is permitted exactly and is different from a zero variance or singular component matrix. Changing the parameterization does not remove statistical boundary problems.",
    "Available denotes numerical completion only. Chi-square calibration is approximate and can under-cover with few facet levels; more persons alone do not supply more independent raters or tasks. Finite-sample coverage, especially for incomplete sources and singular component matrices, is not established. The method is not selected automatically and does not replace the balanced GCI or conservative regression-region methods.",
    "Search from the REML center to the first checked cutoff crossing in each direction. The connected numerical interval does not enumerate remote accepted regions. Multiple starts, a second optimizer, derivative and covariance-cone checks do not prove global optima or exclude undetected crossings. Unresolved fits, evaluations, endpoints or tails leave accepted bounds unavailable; a finite search cap is not an infinite endpoint.",
    "These are pointwise normal-model intervals, not simultaneous or distribution robust. Raw used ratings are required for fitting. Save the complete result with saveRDS() to retain diagnostics and assumptions; saved plot conversion does not refit. Diagnostic covariance units are standardized: multiply by both reml_source$score_scale entries to restore score-product units.",
    "A likelihood-ratio cutoff below numerical resolution, or a nonzero interval lost when restoring score units, leaves bounds unavailable with a specific Status. Estimates and representable REML centers are retained; diagnostics$cutoff records the criterion tolerance and resolution.")
  out
}
