# Remove all non-residual effects; no likelihood fit or n-by-n projector.
.mfrm_mvgt_residual_block <- function(g) {
  scores <- g$design$scores; ids <- c(g$design$person, unname(g$design$facets))
  if (!inherits(g, "mfrm_multivariate_gstudy") || !is.data.frame(g$data) ||
      !all(c(ids, scores) %in% names(g$data)) || anyNA(g$data[ids]) ||
      anyDuplicated(g$data[ids]) || nrow(g$data) != g$design$rows ||
      nrow(g$data) != g$data_usage$counts[["UsedRows"]] ||
      !all(vapply(g$data[scores], function(x) is.numeric(x) && !is.complex(x) &&
        all(is.finite(x)), logical(1))))
    stop("Residual contrasts require the original used ratings and assignments; recreate the G-study.", call. = FALSE)
  if (identical(g$design$method, "Balanced multivariate ANOVA")) {
    .mfrm_mvgt_anova_summary(g)
    counts <- vapply(g$data[ids], function(x) length(unique(x)), integer(1))
    if (!is.null(g$design$nesting)) {
      labels <- c("Person", names(g$design$facets))
      child <- match(names(g$design$nesting), labels)
      parent <- match(unname(g$design$nesting), labels)
      counts[child] <- nrow(unique(g$data[ids[c(child, parent)]]))
    }
    if (!isTRUE(all.equal(unname(counts), unname(g$design$counts), tolerance = 0)))
      stop("Stored design counts do not match the used ratings.", call. = FALSE)
    groups <- .mfrm_mvgt_joint_groups(g)
  } else if (identical(g$design$method, "Multivariate MINQUE(0)")) {
    groups <- .mfrm_mvgt_moment_groups(g)
  } else stop("Residual contrasts require a supported ANOVA or MINQUE(0) G-study.", call. = FALSE)
  names(groups) <- names(g$components)
  y <- as.matrix(g$data[scores]); unit <- apply(abs(y), 2L, max)
  unit[unit == 0] <- 1
  y <- sweep(y, 2L, unit, "/"); y <- sweep(y, 2L, colMeans(y), "-")
  # Check saved points against these ratings on a common scale, not by refitting.
  reconstructed <- .mfrm_mvgt_minque0(y, groups)$components
  stored <- lapply(g$components, function(m) sweep(sweep(m, 1L, unit, "/"), 2L, unit, "/"))
  if (!identical(names(stored), names(reconstructed)) ||
      any(!is.finite(unlist(stored))) ||
      max(abs(unlist(stored)-unlist(reconstructed))) > 1e-8*max(1, abs(unlist(stored))))
    stop("Stored components disagree with their used ratings; recreate the G-study.", call. = FALSE)
  n <- nrow(y)
  z <- do.call(cbind, lapply(groups[-length(groups)], function(a) {
    block <- matrix(0, n, max(a)); block[cbind(seq_len(n), a)] <- 1
    block
  }))
  decomposition <- qr(z, tol = 1e-10)
  singular <- svd(z, nu = 0L, nv = 0L)$d
  rank <- decomposition$rank; df <- n-rank
  ans <- list(df = df, rank = rank, rows = n, columns = ncol(z),
    singular_values = singular, score_scale = unit, valid = FALSE,
    status = "Numerically ambiguous nuisance rank")
  if (rank != sum(singular > max(singular)*1e-10)) return(ans)
  if (df == 0L) { ans$status <- "No residual contrasts in this assignment"; return(ans) }
  residual <- qr.resid(decomposition, y)
  ans$orthogonality_error <- max(abs(crossprod(z, residual)))/max(1, sqrt(sum(y*y)))
  if (ans$orthogonality_error > 1e-10) {
    ans$status <- "Residual projection failed numerical checks"; return(ans)
  }
  ans$mean_product <- crossprod(residual)/df
  ans$score_energy <- colSums(y*y)
  ans$relative_resolution <- (64*.Machine$double.eps*max(dim(z)))^2
  ans$valid <- TRUE; ans$status <- "Checked residual projection"
  ans
}

.mfrm_mvgt_residual_intervals <- function(g, parm, score, weights, level) {
  scores <- g$design$scores; sources <- names(g$components)
  selected <- .mfrm_mvgt_component_selection(parm, sources)
  pair <- is.character(score) && length(score) == 2L
  if (pair) {
    if (!is.null(weights) || anyNA(score) || anyDuplicated(score) || !all(score %in% scores))
      stop("Select two distinct original scores without weights.", call. = FALSE)
    kind <- "Covariance"; label <- paste(score, collapse = " / ")
    estimate <- vapply(g$components, function(m) m[score[1L], score[2L]], numeric(1))
  } else {
    selection <- .mfrm_mvgt_interval_score(scores, score, weights)
    weights <- selection$weights; kind <- selection$kind; label <- selection$score
    scale <- max(abs(weights)); w <- weights/scale
    estimate <- vapply(g$components, function(m) drop(crossprod(w, m %*% w))*scale*scale, numeric(1))
  }
  if (any(!is.finite(estimate))) stop("Component units exceed numeric range; rescale scores or weights.", call. = FALSE)
  block <- .mfrm_mvgt_residual_block(g)
  center <- lower <- upper <- rep(NA_real_, length(sources))
  status <- rep("Method applies only to Residual", length(sources)); j <- match("Residual", sources)
  status[j] <- block$status
  if (block$valid) {
    if (pair) {
      s <- block$mean_product[score, score, drop = FALSE]
      unit <- block$score_scale[score]
      center[j] <- s[1L, 2L]*unit[1L]*unit[2L]
      usable <- all(diag(s)*block$df > block$relative_resolution*block$score_energy[score])
      bounds <- c(NA_real_, NA_real_)
      status[j] <- "Residual information at numerical resolution"
      if (usable) {
        region <- .mfrm_mvgt_wishart_covariance_region(matrix(c(s[1,1], s[1,2], s[2,2]), 1L), block$df, level)
        bounds <- region$bounds[1L, ]*unit[1L]*unit[2L]
        block$region <- region
        status[j] <- if (all(is.finite(region$bounds))) "Conservative normal-model interval" else
          paste("Unavailable:",region$method[1L])
        if (all(is.finite(region$bounds)) &&
            (any(!is.finite(bounds)) || bounds[1L] == bounds[2L] ||
              any(bounds == 0 & region$bounds[1L, ] != 0))) {
          bounds[] <- NA_real_; status[j] <- "Unrepresentable residual interval"
        }
      }
    } else {
      a <- weights*block$score_scale; scale <- max(abs(a))
      if (!is.finite(scale) || scale == 0) stop("Composite units exceed numeric range; rescale scores or weights.", call. = FALSE)
      a <- a/scale
      ms <- drop(crossprod(a, block$mean_product %*% a))
      center[j] <- ms*scale*scale
      energy <- sum(abs(a)*sqrt(block$score_energy))^2
      usable <- ms*block$df > block$relative_resolution*energy
      bounds <- c(NA_real_, NA_real_)
      status[j] <- "Residual information at numerical resolution"
      if (usable) {
        interval <- .mfrm_mvgt_linear_variance_interval(ms, block$df, 1, level)
        bounds <- interval$bounds*scale*scale
        status[j] <- if (interval$status == "Available") "Exact normal-model interval" else interval$status
      }
    }
    if (usable && status[j] %in% c("Exact normal-model interval", "Conservative normal-model interval") &&
        (any(!is.finite(bounds)) || (!pair && (any(bounds == 0) || bounds[1L] == bounds[2L])))) {
      bounds[] <- NA_real_; status[j] <- "Unrepresentable residual interval"
    }
    if (!is.finite(center[j])) center[j] <- NA_real_
    lower[j] <- bounds[1L]; upper[j] <- bounds[2L]
  }
  out <- data.frame(Source = sources, Score = label, Estimate = unname(estimate),
    SE = NA_real_, SEStatus = "Not supplied for the original estimator",
    IntervalCenter = center, Lower = lower, Upper = upper,
    RawLower = lower, RawUpper = upper, Level = level,
    Method = if (pair) "Residual-contrast covariance enclosure" else "Residual-contrast chi-square",
    MethodCode = "residual_contrast", Status = status,
    Available = is.finite(lower) & is.finite(upper),
    Validation = if (pair) "Conservative pointwise normal-model coverage; may be very wide" else
      "Exact pointwise normal-model coverage for positive residual variance")
  if (pair) { out$Score1 <- score[1L]; out$Score2 <- score[2L] }
  out$Validation[!out$Available] <- "No interval supplied"
  out <- out[selected, , drop = FALSE]; rownames(out) <- NULL
  attr(out, "kind") <- kind; attr(out, "weights") <- weights
  if (pair) attr(out, "score_pair") <- score
  attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "nesting", "child_counts",
    "balanced", "complete", "method", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "residual_projection") <- block
  attr(out, "notes") <- c(
    "Residual-only intervals conditional on the observed assignment under independent homoskedastic normal residuals. All non-residual incidence columns, including interactions and the mean, are projected out. Omission must preserve this model.",
    "Estimate retains the original ANOVA/MINQUE point; IntervalCenter is the separate residual-contrast estimate. SE is not supplied for the original estimator; use vcov(g) for its separately qualified sampling SEs. No likelihood optimization or bootstrap is used.",
    if (pair) "The residual mean-product matrix has a Wishart law. Independent regression-t and chi-square pivots give a conservative covariance enclosure; one-df or singular pairs use marginal variance limits and Cauchy-Schwarz. It may be very wide, is not an exact equal-tail interval, and retains signed covariance." else
      "A chi-square pivot gives an exact pointwise interval for positive residual variance under the normal model. This is not distribution robust. Component estimates and limits have squared score units, not measurement-SEM units.",
    "If an interval loses its width or nonzero endpoints through numerical rounding, its limits remain unavailable with an explicit Status; the original estimate and projected center are retained.",
    "Other component intervals are unavailable under this method. No residual contrasts, numerically unresolved residual variation and projection failures remain explicit; no interval is not zero uncertainty. Residual combines the highest interaction and observation error, which cannot be separated without cell replication.",
    "Coverage statements are pointwise for prespecified scores/weights, not simultaneous or selection adjusted. Save the complete result and residual_projection diagnostics with saveRDS().")
  out
}
