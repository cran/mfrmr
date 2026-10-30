#' Sampling covariance and SEs of G-study component estimates
#'
#' Quantify how estimated sources of score variation would vary across repeated
#' G-studies with the same assignment design. Supports balanced ANOVA and
#' identifiable MINQUE(0), including incomplete and supported nested sources.
#'
#' @param object A result from [mfrm_g_study()].
#' @param score One original score, or exactly two distinct original scores
#'   for the sampling covariance of their between-score component covariances.
#'   The default selects the first original score.
#' @param weights Optional named weights for one prespecified composite;
#'   name every original score. Do not also supply `score`.
#' @param assumption Required: `"normal"` asserts independent normal random
#'   effects, a common mean for each score and constant component covariances.
#'   Assignment and omissions must preserve that model when conditioning on
#'   the observed assignment. The assumption is not tested.
#' @param ... Reserved; extra arguments are rejected.
#'
#' @details
#' Rows and columns refer to **estimated components**, such as Person, Task
#' and Residual. This is neither the covariance matrix between individual
#' ratings nor the component matrices between score columns in `object`.
#' Use `attr(result, "summary")` to read the original estimates, sampling
#' variance estimates, SEs and availability. A component SE has squared score
#' units; the returned covariance matrix has fourth-power score units. For
#' a between-score covariance, the SE has the product of the two score units.
#' These SEs differ from the D-study's measurement SEM in original score units.
#'
#' The component estimates are unchanged. This method uses Gaussian product
#' moments to correct plug-in bias in their estimated sampling covariance.
#' It does not fit a likelihood, resample data or construct confidence intervals.
#' In particular, `Estimate +/- 1.96 * SE` is not a generally justified interval.
#' Use [confint.mfrm_multivariate_gstudy()] for its separately qualified
#' interval methods. A REML profile's center belongs to a different estimator;
#' the SE returned here describes the original ANOVA/MINQUE estimator.
#'
#' @section Calculation:
#' With mean-removing MINQUE quadratic forms \eqn{q_s=y^T A_s y},
#' \eqn{E(q_s)=\theta_s} and
#' \eqn{C_{st}=2\mathrm{tr}(A_s V A_t V)}, where
#' \eqn{V=\sum_u\theta_u K_u}. Let the linear operator P map
#' \eqn{\mathrm{vec}(\theta\theta^T)} to \eqn{\mathrm{vec}(C)}.
#' The normal product-moment identity gives
#' \eqn{E\{\mathrm{vec}(qq^T)\}=(I+P)\mathrm{vec}(\theta\theta^T)}.
#' Inverting this moment system gives an unbiased estimator of C for its
#' identifiable entries. This is a consequence of the cited moment identity,
#' not a claim that the reference proposes this package's G-theory procedure.
#' The squared SE, when available, is the corresponding diagonal estimate;
#' the square root is not itself an unbiased SE estimator.
#'
#' For two scores, write x and y for their variance-component vectors and c
#' for their covariance-component vector. On symmetric matrix coordinates,
#' \eqn{E\{\mathrm{sym}(\hat x\hat y^T)\}=\mathrm{sym}(xy^T)+P(cc^T)} and
#' \eqn{E(\hat c\hat c^T)=(I+P/2)(cc^T)+(P/2)\mathrm{sym}(xy^T)}.
#' Solve these coupled equations for the sampling covariance
#' \eqn{(P/2)\{\mathrm{sym}(xy^T)+cc^T\}}. Applying only the scalar
#' correction to a covariance component would be incorrect.
#'
#' Balanced ANOVA reduces to the existing independent-mean-square result
#' `A %*% diag(2 * MS^2 / (df + 2)) %*% t(A)` for a variance/fixed composite,
#' and the Wishart product-moment correction for a between-score covariance.
#' MINQUE uses the actual observed grouping kernels. Column-block calculations
#' avoid allocating observation-square matrices, but time still increases
#' quadratically with the number of observed ratings. No sampling simulation
#' is performed. Balanced MINQUE uses the same ANOVA limit.
#'
#' @section Unavailable or negative results:
#' A moment system can fail to identify some sampling-covariance entries even
#' when point components are identifiable. A contributing balanced stratum
#' with df one is an example for between-score covariance SEs. An SVD checks
#' each target entry's estimability; only unidentified entries become `NA`.
#' Its rank tolerance is `sqrt(.Machine$double.eps)` relative to the largest
#' singular value, with relative target reconstruction error below `1e-8`.
#'
#' An unbiased sampling-variance estimate may be negative. It is retained in
#' the returned matrix and `SamplingVariance`, while its `SE` is `NA` and
#' `SEStatus` explains why. The matrix can also be indefinite with positive
#' diagonal entries. Nothing is clipped or repaired to make it positive
#' semidefinite. A mathematically unbiased variance estimate does not ensure
#' a useful SE in a small sample or correct normal-reference coverage.
#'
#' @return A named matrix estimating the sampling covariance of component
#'   estimates. Attribute `summary` is a data frame with `Source`, `Estimate`,
#'   `SamplingVariance`, `SE`, and `SEStatus`. `diagnostics` retains the moment
#'   rank, singular values, entrywise estimability and matrix PSD check.
#'   `score`, `weights`, `assumption`, `source_design`, `source_rows` and `notes`
#'   retain interpretation. Save the full result using `saveRDS()`.
#'   MINQUE requires its retained grouping assignments, moment equations and
#'   quadratic products; balanced ANOVA uses saved summaries without raw data.
#' @references
#' Magnus, J. R. (1978). The moments of products of quadratic forms in normal
#' variables. *Statistica Neerlandica, 32*, 201--210. Lemma 6.2, p.209.
#' \doi{10.1111/j.1467-9574.1978.tb01399.x}.
#'
#' Rao, C. R. (1971). Estimation of variance and covariance components--MINQUE
#' theory. *Journal of Multivariate Analysis, 1*, 257--275.
#' \doi{10.1016/0047-259X(71)90001-7}.
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("vcov.mfrm_multivariate_gstudy", run.dontrun = TRUE)
#' \dontrun{
#' data("mfrmr_example_gtheory_crossed", package = "mfrmr")
#' ratings <- mfrmr_example_gtheory_crossed
#' # Complete balanced ratings: compute both targets directly.
#' g <- mfrm_g_study(ratings, c("Content", "Organization"),
#'                             method = "anova")
#' sampling <- vcov(g, score = "Content", assumption = "normal")
#' attr(sampling, "summary") # Component estimates and their sampling SEs.
#' attr(sampling, "notes")
#' # Uncertainty in the covariance between two criteria, by source.
#' cross_score <- vcov(g, score = c("Content", "Organization"),
#'                     assumption = "normal")
#' attr(cross_score, "summary")
#'
#' # Incomplete assignments also work, but their moment calculation is longer.
#' # These saved results use the same ratings with rows 1, 8 and 23 omitted.
#' # Omitted ratings are absent observations, never zero scores.
#' incomplete <- readRDS(system.file("examples", "gstudy-sampling.rds",
#'                                   package = "mfrmr", mustWork = TRUE))
#' incomplete$study$design
#' attr(incomplete$variance, "summary")
#' attr(incomplete$covariance, "summary")
#' attr(incomplete$covariance, "notes")
#' # To repeat the full incomplete-data calculation:
#' # source(system.file("examples", "gstudy-sampling.R", package = "mfrmr"))
#' # incomplete <- gstudy_sampling_example()
#' }
#' @export
vcov.mfrm_multivariate_gstudy <- function(object, score = NULL, weights = NULL,
                                         assumption, ...) {
  rlang::check_dots_empty()
  if (missing(assumption) || !identical(assumption, "normal"))
    stop("Sampling covariance requires `assumption = 'normal'`; the assumption is not tested.", call. = FALSE)
  g <- object; scores <- g$design$scores; sources <- names(g$components)
  if (!inherits(g, "mfrm_multivariate_gstudy") ||
      !is.character(scores) || !length(scores) || anyNA(scores) || anyDuplicated(scores) ||
      !.mfrm_mvgt_valid_facets(names(g$design$facets)) ||
      !identical(sources, .mfrm_mvgt_sources(c("Person", names(g$design$facets)), g$design$nesting)) ||
      !all(vapply(g$components, function(m) is.matrix(m) && is.numeric(m) &&
        !is.complex(m) && identical(dimnames(m), list(scores, scores)) &&
        all(is.finite(m)) && identical(m, t(m)), logical(1))))
    stop("The G-study components or identities are incomplete or altered; recreate the result.", call. = FALSE)
  pair <- is.character(score) && length(score) == 2L
  if (pair) {
    if (!is.null(weights) || anyNA(score) || anyDuplicated(score) || !all(score %in% scores))
      stop("Select two distinct original scores, without weights, for covariance-component SEs.", call. = FALSE)
  } else {
    selection <- .mfrm_mvgt_interval_score(scores, score, weights)
    score <- selection$score; weights <- selection$weights
  }
  k <- length(sources)
  if (identical(g$design$method, "Balanced multivariate ANOVA")) {
    design <- .mfrm_mvgt_anova_summary(g)
    for (m in g$mean_products) {
      unit <- max(abs(m))
      if (unit > 0 && min(eigen(m/unit, symmetric=TRUE, only.values=TRUE)$values) <
          -64*.Machine$double.eps*nrow(m))
        stop("Stored ANOVA mean products must be positive semidefinite; recreate the G-study.",call.=FALSE)
    }
    # Validate the stored point/mean-product correspondence for every score pair.
    mapped <- design$A %*% do.call(rbind, lapply(g$mean_products, as.vector))
    original <- do.call(rbind, lapply(g$components, as.vector))
    if (!isTRUE(all.equal(unname(mapped), unname(original), tolerance = 1e-10)))
      stop("Stored components disagree with their mean products; recreate the G-study.", call. = FALSE)
  } else if (identical(g$design$method, "Multivariate MINQUE(0)")) {
    groups <- .mfrm_mvgt_moment_groups(g)
    design <- if (isTRUE(g$design$balanced) && isTRUE(g$design$complete))
      .mfrm_mvgt_anova_interval_design(g$design) else NULL
  } else stop("Unsupported G-study estimator.", call. = FALSE)
  if (!is.null(design)) {
    A <- design$A; L <- solve(A)
    left <- vapply(seq_len(k), function(j) as.vector(tcrossprod(A[,j])), numeric(k*k))
    right <- t(vapply(seq_len(k), function(j) as.vector(tcrossprod(L[j,])), numeric(k*k)))
    P <- sweep(left, 2L, 2/design$df, "*") %*% right
  } else P <- .mfrm_mvgt_moment_polynomial(groups, g$estimation$kernel_gram)
  map <- .mfrm_mvgt_moment_map(P, pair)
  if (pair) {
    v1 <- vapply(g$components, function(m) m[score[1L], score[1L]], numeric(1))
    v2 <- vapply(g$components, function(m) m[score[2L], score[2L]], numeric(1))
    estimate <- vapply(g$components, function(m) m[score[1L], score[2L]], numeric(1))
    unit1 <- max(abs(v1)); unit2 <- max(abs(v2))
    if (unit1 == 0) unit1 <- 1
    if (unit2 == 0) unit2 <- 1
    unit <- sqrt(unit1)*sqrt(unit2)
    result <- .mfrm_mvgt_moment_correct(map, estimate/unit, v1/unit1, v2/unit2)
  } else {
    ws <- max(abs(weights)); w <- weights/ws
    active <- which(w != 0)
    ms <- max(abs(unlist(lapply(g$components, function(m) m[active,active,drop=FALSE]))))
    if (ms == 0) ms <- 1
    theta <- vapply(g$components, function(m)
      drop(crossprod(w[active], (m[active,active,drop=FALSE]/ms) %*% w[active])), numeric(1))
    unit <- ms*ws*ws; estimate <- theta*unit
    result <- .mfrm_mvgt_moment_correct(map, theta)
    if (!is.null(design)) {
      # The independent-MS form avoids cancellation at zero variance boundaries.
      means <- if (identical(g$design$method, "Balanced multivariate ANOVA"))
        vapply(g$mean_products, function(m)
          drop(crossprod(w[active], (m[active,active,drop=FALSE]/ms) %*% w[active])), numeric(1)) else
        drop(L %*% theta)
      result$covariance <- tcrossprod(sweep(A, 2L, means*sqrt(2/(design$df+2)), "*"))
      result$se <- sqrt(diag(result$covariance))
      result$estimable[] <- TRUE
    }
  }
  if (!is.finite(unit) || unit <= 0 || any(!is.finite(estimate)))
    stop("Sampling-covariance units exceed numeric range; rescale scores or weights.", call. = FALSE)
  out <- (result$covariance*unit)*unit
  if (any(is.infinite(out)) || any(!is.na(out) & result$covariance != 0 & out == 0))
    stop("Sampling covariance exceeds numeric range; rescale scores or weights.", call. = FALSE)
  dimnames(out) <- list(sources, sources)
  variance <- diag(out); se <- result$se*unit
  status <- ifelse(!diag(result$estimable), "Sampling variance not identifiable",
    ifelse(variance < 0, "Negative sampling-variance estimate",
      ifelse(variance == 0, "Zero sampling-variance estimate", "Available")))
  attr(out, "summary") <- data.frame(Source=sources, Estimate=unname(estimate),
    SamplingVariance=unname(variance), SE=unname(se), SEStatus=unname(status))
  eigenvalues <- if (all(is.finite(result$covariance)))
    eigen(result$covariance, symmetric=TRUE, only.values=TRUE)$values else NA_real_
  attr(out, "diagnostics") <- list(moment_rank=map$rank, moment_dimension=map$dimension,
    singular_values=map$singular_values, estimable=result$estimable,
    reconstruction_error=map$error,
    positive_semidefinite=if (all(is.finite(eigenvalues)))
      min(eigenvalues) >= -sqrt(.Machine$double.eps)*max(abs(eigenvalues)) else NA,
    normalized_covariance=result$covariance, component_scale=unit)
  attr(out, "score") <- score; attr(out, "weights") <- weights
  attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "notes") <- c(
    "Sampling covariance of original ANOVA/MINQUE component estimates, conditional on the declared assignment under independent normal random effects. This is not measurement SEM, covariance between individual ratings or the SE of a separate REML estimator.",
    "Gaussian product-moment correction gives unbiased identifiable sampling-covariance entries under the model. Taking square roots is not an unbiasedness guarantee for SEs and does not justify normal-reference confidence intervals.",
    "Negative sampling variances are retained, with unavailable SEs. Unidentified entries are NA. An indefinite matrix is not repaired; inspect summary and diagnostics before propagating it to another target.",
    "Scores and composite weights are prespecified; omission must preserve the model conditional on assignment. No robustness, informative-missingness correction or general confidence-interval coverage is claimed. Save the complete result with saveRDS().")
  out
}

.mfrm_mvgt_moment_groups <- function(g) {
  ids <- c(g$design$person, unname(g$design$facets)); data <- g$data
  if (anyDuplicated(ids) || !is.data.frame(data) || !all(ids %in% names(data)) || anyNA(data[ids]) ||
      nrow(data) != g$design$rows || anyDuplicated(data[ids]))
    stop("MINQUE sampling covariance needs the retained observed assignments; recreate the G-study.", call. = FALSE)
  factors <- lapply(data[ids], factor)
  actual_counts <- vapply(factors,nlevels,integer(1))
  subsets <- if (length(ids)==2L) list(1L,2L,1:2) else
    list(1L,2L,3L,c(1L,2L),c(1L,3L),c(2L,3L),1:3)
  if (!is.null(g$design$nesting)) {
    labels <- c("Person",names(g$design$facets))
    child <- match(names(g$design$nesting),labels); parent <- match(unname(g$design$nesting),labels)
    subsets <- list(1L,parent,c(child,parent),c(1L,parent),1:3)
    children <- vapply(split(factors[[child]],factors[[parent]]),function(x)length(unique(x)),integer(1))
    actual_counts[child] <- sum(children)
  }
  potential <- if (is.null(g$design$nesting)) prod(as.double(actual_counts)) else
    actual_counts[1L]*actual_counts[child]
  complete <- unname(nrow(data)==potential)
  balanced <- complete && (is.null(g$design$nesting) || length(unique(children))==1L)
  if (!identical(as.numeric(g$design$counts),as.numeric(actual_counts)) ||
      !identical(isTRUE(g$design$complete),complete) ||
      !identical(isTRUE(g$design$balanced),balanced) ||
      (!is.null(g$design$nesting) && !identical(as.numeric(g$design$child_counts),as.numeric(children))))
    stop("Stored design counts or completeness do not match the observed assignments.",call.=FALSE)
  groups <- lapply(subsets,function(s).mfrm_mvgt_group(factors[s]))
  k <- length(groups); N <- nrow(data)
  # Recompute assignment-only traces; neither scores nor point estimates are refit.
  size <- lapply(groups,function(g)as.double(tabulate(g)[g]))
  total <- vapply(size,sum,numeric(1));S <- matrix(0,k,k)
  for(i in seq_len(k))for(j in seq_len(i)) {
    counts <- tabulate(.mfrm_mvgt_group(list(groups[[i]],groups[[j]])))
    S[i,j] <- S[j,i] <- sum(as.double(counts)^2)-2*sum(size[[i]]*size[[j]])/N+total[i]*total[j]/N^2
  }
  gram <- g$estimation$kernel_gram; products <- g$estimation$quadratic_products
  if (!is.matrix(gram) || !identical(dim(gram),c(k,k)) || any(!is.finite(gram)) ||
      !isTRUE(all.equal(unname(gram),S,tolerance=1e-10)) ||
      !is.list(products) || length(products)!=k)
    stop("Stored MINQUE equations do not match the observed assignment; recreate the G-study.", call. = FALSE)
  rhs <- do.call(rbind,lapply(products,as.vector))
  theta <- do.call(rbind,lapply(g$components,as.vector))
  if (!isTRUE(all.equal(unname(gram%*%theta),unname(rhs),tolerance=1e-10)))
    stop("Stored MINQUE components disagree with their quadratic products; recreate the G-study.", call. = FALSE)
  groups
}

# Gaussian product-moment correction for component sampling covariance.
.mfrm_mvgt_moment_polynomial <- function(groups, gram, block=32L) {
  N <- length(groups[[1L]]); k <- length(groups)
  scale <- sqrt(diag(gram))
  inverse <- solve(gram/outer(scale,scale))/outer(scale,scale)
  apply_K <- function(X,g) rowsum(X,g)[g,,drop=FALSE]
  apply_A <- function(X) {
    X <- sweep(X,2,colMeans(X),'-')
    B <- vapply(groups,function(g) {
      x <- apply_K(X,g)
      as.vector(sweep(x,2,colMeans(x),'-'))
    },numeric(length(X)))
    B%*%t(inverse)
  }
  trace <- matrix(0,k*k,k*k)
  block <- max(1L,min(block,floor(1e6/(N*k*k))))
  for(first in seq.int(1L,N,by=block)) {
    columns <- seq.int(first,min(N,first+block-1L))
    E <- matrix(0,N,length(columns));E[cbind(columns,seq_along(columns))]<-1
    right <- vapply(groups,function(g)as.vector(apply_A(apply_K(E,g))),numeric(length(E)*k))
    A <- matrix(apply_A(E),nrow=N)
    left <- vapply(groups,function(g)as.vector(apply_K(A,g)),numeric(length(E)*k))
    trace <- trace+crossprod(matrix(right,nrow=length(E)),matrix(left,nrow=length(E)))
  }
  P <- 2*matrix(aperm(array(trace,rep(k,4)),c(1,3,2,4)),k*k,k*k)
  transposed <- as.vector(t(matrix(seq_len(k*k),k)))
  P <- (P+P[transposed,])/2
  P <- (P+P[,transposed])/2
  P
}

.mfrm_mvgt_moment_symmetric <- function(P) {
  k <- as.integer(sqrt(nrow(P)))
  ij <- which(upper.tri(matrix(0,k,k),diag=TRUE),arr.ind=TRUE)
  index <- ij[,1]+(ij[,2]-1L)*k
  E <- matrix(0,k*k,nrow(ij))
  for(j in seq_len(nrow(ij))) {
    x<-matrix(0,k,k);x[ij[j,1],ij[j,2]]<-x[ij[j,2],ij[j,1]]<-1
    E[,j]<-as.vector(x)
  }
  list(k=k,index=index,expand=E,L=P[index,,drop=FALSE]%*%E,
    diagonal=which(ij[,1]==ij[,2]))
}

.mfrm_mvgt_moment_map <- function(P,pair=FALSE) {
  s<-.mfrm_mvgt_moment_symmetric(P);L<-s$L;I<-diag(nrow(L))
  M<-if(pair)rbind(cbind(I,L),cbind(L/2,I+L/2)) else I+L
  target<-if(pair)cbind(L/2,L/2) else L
  singular<-svd(M);keep<-singular$d>sqrt(.Machine$double.eps)*max(singular$d)
  inverse<-sweep(singular$v[,keep,drop=FALSE],2,singular$d[keep],'/')%*%t(singular$u[,keep,drop=FALSE])
  weights<-target%*%inverse
  error<-sqrt(rowSums((weights%*%M-target)^2))/pmax(sqrt(rowSums(target^2)),.Machine$double.eps)
  list(symmetric=s,weights=weights,estimable=error<1e-8,error=error,rank=sum(keep),
    dimension=nrow(M),singular_values=singular$d,moment=M,target=target,pair=pair)
}

.mfrm_mvgt_moment_correct <- function(map,theta,variance1=NULL,variance2=NULL) {
  s<-map$symmetric
  obs<-as.vector(tcrossprod(theta))[s$index]
  if(map$pair) {
    B<-(outer(variance1,variance2)+outer(variance2,variance1))/2
    obs<-c(B[s$index],obs)
  }
  corrected<-drop(map$weights%*%obs)
  raw<-matrix(s$expand%*%corrected,s$k)
  available<-matrix(s$expand%*%as.numeric(map$estimable),s$k)>0
  value<-raw;value[!available]<-NA_real_
  list(covariance=value,raw_covariance=raw,estimable=available,
    se=ifelse(is.finite(diag(value))&diag(value)>=0,sqrt(pmax(diag(value),0)),NA_real_))
}
