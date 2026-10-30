.mfrmr_gt_geometry <- function(data, terms, reml) {
  groups <- lapply(data[terms], function(x) as.integer(factor(x)))
  groups$Residual <- seq_len(nrow(data))
  gram <- .mfrm_gt_kernel_gram(groups, centered = reml)
  scale <- sqrt(diag(gram))
  if (any(!is.finite(scale)) || any(scale <= 0))
    stop("The observed design cannot distinguish the requested variance components.", call. = FALSE)
  eig <- eigen(gram / outer(scale, scale), symmetric = TRUE)
  tolerance <- sqrt(.Machine$double.eps)
  keep <- eig$values > tolerance*max(eig$values)
  null <- eig$vectors[, !keep, drop = FALSE]
  rownames(null) <- names(groups)
  list(kernel_gram = gram, scale = scale, null_space = null,
    rank = sum(keep), components = length(groups), scaled_eigenvalues = eig$values,
    relative_tolerance = tolerance,
    basis = if (reml) "REML score contrasts" else "ML marginal covariance")
}

.mfrmr_gt_estimable <- function(weights, geometry) {
  if (is.null(geometry)) return(NA)
  if (!identical(names(weights), names(geometry$scale)))
    stop("G-study component identities do not match the saved design information; recreate the G-study.", call. = FALSE)
  if (!ncol(geometry$null_space)) return(TRUE)
  z <- weights / geometry$scale
  if (!any(z != 0)) return(TRUE)
  z <- z / max(abs(z))
  sqrt(sum(crossprod(geometry$null_space, z)^2)) <=
    geometry$relative_tolerance*sqrt(sum(z^2))
}

.mfrmr_gt_projection <- function(v, object, facets, counts, scaling, geometry = NULL) {
  sources <- c(object, facets, "Residual")
  v <- stats::setNames(as.numeric(v[sources]), sources)
  valid <- is.finite(v) & v >= 0
  v[!valid] <- NA_real_
  unit <- stats::setNames(rep(0, length(sources)), sources)
  u <- unit; u[object] <- 1
  e <- unit; e["Residual"] <- 1
  object_ok <- .mfrmr_gt_estimable(u, geometry)
  error_ok <- .mfrmr_gt_estimable(e, geometry)
  rows <- lapply(seq_len(nrow(counts)), function(i) {
    n <- as.numeric(counts[i, ])
    divisors <- switch(scaling[i], highest_order = n, single_condition = min(n), none = 1)
    divisor <- prod(divisors)
    # Divide successively: a product may overflow even when the projected
    # variance is representable (e.g. 1e308 / 1e200 / 1e200).
    residual <- v[["Residual"]]
    for (k in divisors) residual <- residual/k
    main <- v[facets]/n
    a <- unit; a[facets] <- 1/n
    a["Residual"] <- 1/divisor
    absolute_ok <- .mfrmr_gt_estimable(a, geometry)
    object_reason <- if (isFALSE(object_ok)) "object_variance_not_identified" else
      if (!is.finite(v[[object]])) "object_variance_missing_or_invalid" else ""
    relative_reason <- if (isFALSE(error_ok)) "relative_error_not_identified" else
      if (!is.finite(residual)) "residual_variance_missing_or_invalid" else
        if (residual == 0 && isTRUE(v[["Residual"]] > 0)) "relative_error_below_numerical_range" else ""
    absolute_reason <- if (isFALSE(absolute_ok)) "absolute_error_not_identified" else
      if (any(!is.finite(c(main, residual)))) "error_components_missing_or_invalid" else
        if (any(main == 0 & v[facets] > 0) || (residual == 0 && v[["Residual"]] > 0))
          "error_component_below_numerical_range" else ""
    ratio <- function(errors, reason) {
      why <- if (nzchar(object_reason)) object_reason else reason
      if (nzchar(why)) return(list(value = NA_real_, reason = why))
      scale <- max(c(v[[object]], errors))
      if (scale == 0) return(list(value = NA_real_, reason = "zero_total_variance"))
      terms <- c(v[[object]], errors)/scale
      value <- terms[1L]/sum(terms)
      if (value == 0 && v[[object]] > 0)
        return(list(value = NA_real_, reason = "coefficient_below_numerical_range"))
      list(value = value, reason = "")
    }
    g <- ratio(residual, relative_reason)
    phi <- ratio(c(main, residual), absolute_reason)
    absolute <- sum(c(main, residual))
    if (!nzchar(absolute_reason) && !is.finite(absolute))
      absolute_reason <- "absolute_error_exceeds_numerical_range"
    data.frame(ResidualDivisor = if (is.finite(divisor)) divisor else NA_real_,
      ResidualDivisorReason = if (is.finite(divisor)) "" else "divisor_exceeds_numerical_range",
      ObjectVariance = if (!nzchar(object_reason)) v[[object]] else NA_real_,
      ObjectVarianceReason = object_reason,
      RelativeErrorVariance = if (!nzchar(relative_reason)) residual else NA_real_,
      RelativeErrorReason = relative_reason,
      AbsoluteErrorVariance = if (!nzchar(absolute_reason)) absolute else NA_real_,
      AbsoluteErrorReason = absolute_reason,
      G = g$value, Phi = phi$value, GReason = g$reason, PhiReason = phi$reason,
      row.names = NULL)
  })
  do.call(rbind, rows)
}
