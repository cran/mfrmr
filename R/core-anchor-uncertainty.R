# Independent source calibration contributes through the constrained score root.
mfrmr_anchor_uncertainty_request <- function(fit, covariance, independent) {
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent))
    stop("`independent_anchors` must be TRUE or FALSE.", call. = FALSE)
  if (is.null(covariance)) {
    if (independent) stop("Supply `anchor_covariance` when declaring independent anchors.", call. = FALSE)
    return(NULL)
  }
  if (!independent)
    stop("Anchor-source propagation requires `independent_anchors = TRUE`: source calibration and receiving responses must be independent. Shared respondents or clusters require cross-sample covariance, which is unsupported.", call. = FALSE)
  cfg <- fit$config
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !cfg$model %in% c("RSM", "PCM") || !identical(cfg$method, "MML") ||
      is.null(cfg$population_spec$fixed_sd) ||
      !identical(cfg$population_spec$design_columns, "(Intercept)") ||
      identical(cfg$noncenter_facet, "Person") ||
      !identical(cfg$estimation_control$mml_engine_used, "direct") ||
      mfrmr_adaptive_integration(cfg) || length(cfg$interaction_specs) ||
      length(cfg$dummy_facets) || !identical(cfg$facet_shrinkage, "none"))
    stop("Anchor-source intervals require the direct fixed-grid RSM/PCM anchor-origin route, an intercept-only population with fixed SD, and no interactions or shrinkage.", call. = FALSE)
  mfrm_native_location_source(fit)
  anchors <- extract_anchor_tables(cfg)$anchors
  if (!nrow(anchors) || any(anchors$Facet != cfg$noncenter_facet))
    stop("All supplied anchors must be on the one noncentered facet for source-error propagation.", call. = FALSE)
  labels <- anchors$Level
  if (!is.matrix(covariance) || !is.numeric(covariance) || is.complex(covariance) ||
      !identical(dim(covariance), rep(length(labels), 2L)) || any(!is.finite(covariance)) ||
      is.null(rownames(covariance)) || is.null(colnames(covariance)) ||
      anyDuplicated(rownames(covariance)) || anyDuplicated(colnames(covariance)) ||
      !setequal(rownames(covariance), labels) || !setequal(colnames(covariance), labels))
    stop("`anchor_covariance` must be a finite square matrix with row and column names matching every supplied anchor level, in native ability units squared.", call. = FALSE)
  covariance <- covariance[labels, labels, drop = FALSE]
  tolerance <- 128 * .Machine$double.eps * max(abs(covariance))
  if (max(abs(covariance - t(covariance))) > tolerance)
    stop("`anchor_covariance` must be symmetric.", call. = FALSE)
  eig <- eigen(symmetrize_matrix(covariance), symmetric = TRUE)
  if (min(eig$values) < -128 * .Machine$double.eps * max(abs(eig$values)))
    stop("`anchor_covariance` must be positive semidefinite; anchor correlations cannot be omitted or replaced by an invalid covariance.", call. = FALSE)
  root <- sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), "*")
  used <- tcrossprod(root); dimnames(used) <- dimnames(covariance)
  list(anchors = as.data.frame(anchors), covariance = used,
    supplied_covariance = covariance, covariance_roundoff = max(abs(used-covariance)),
    independent = TRUE, source_identity = mfrm_native_location_identity(fit))
}

mfrmr_anchor_sensitivity <- function(fit, covariance, anchors) {
  cfg <- fit$config; sizes <- build_param_sizes(cfg)
  idx <- build_indices(fit$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
  p <- fit$opt$par; owner <- cfg$noncenter_facet
  q <- cfg$estimation_control$quad_points
  cross_information <- function(order, relative_step) {
    quad <- gauss_hermite_normal(order)
    vapply(seq_len(nrow(anchors)), function(j) {
      h <- relative_step * max(1, abs(anchors$Anchor[j]))
      plus <- minus <- cfg
      plus$facet_specs[[owner]]$anchors[anchors$Level[j]] <- anchors$Anchor[j] + h
      minus$facet_specs[[owner]]$anchors[anchors$Level[j]] <- anchors$Anchor[j] - h
      (mfrm_grad_mml(p, idx, plus, sizes, quad) -
        mfrm_grad_mml(p, idx, minus, sizes, quad)) / (2*h)
    }, numeric(length(p)))
  }
  coarse <- -covariance %*% cross_information(q, 1e-4)
  sensitivity <- -covariance %*% cross_information(q, 5e-5)
  fine <- -covariance %*% cross_information(2*q-1L, 5e-5)
  denominator <- max(1, norm(sensitivity, "F"))
  step_change <- norm(coarse-sensitivity, "F") / denominator
  grid_change <- norm(fine-sensitivity, "F") / denominator
  checks <- data.frame(Check = c("Anchor derivative step refinement", "Anchor derivative quadrature refinement"),
    Passed = c(is.finite(step_change) && step_change <= 1e-5,
      is.finite(grid_change) && grid_change <= 1e-3),
    Value = c(step_change, grid_change), Tolerance = c(1e-5, 1e-3))
  dimnames(sensitivity) <- list(mfrm_checkpoint_parameter_names(sizes), anchors$Level)
  list(sensitivity = sensitivity, checks = checks, eligible = all(checks$Passed))
}

mfrmr_anchor_uncertainty_validate_saved <- function(x) {
  a <- x$anchor_uncertainty
  if (is.null(a)) return(invisible(NULL))
  fingerprint <- a$fingerprint; a$fingerprint <- NULL
  if (!identical(a$source_identity, mfrm_native_location_identity(x$fit)) ||
      !identical(fingerprint, mfrm_checkpoint_fingerprint(a, canonical = TRUE)))
    stop("Saved anchor-source uncertainty must match its unchanged fit, covariance, sensitivity and qualification record.", call. = FALSE)
  invisible(NULL)
}

mfrmr_anchor_uncertainty_intervals <- function(out, request) {
  if (is.null(request)) return(out)
  fit <- out$fit; population <- inherits(out, "mfrm_population_intervals")
  V <- out$parameter_covariance
  a <- request
  a$conditional_parameter_covariance <- V
  a$fingerprint <- NULL
  if (is.null(a$checks) && !is.null(V)) {
    calculation <- tryCatch(mfrmr_anchor_sensitivity(fit, V, a$anchors),
      error = function(e) list(sensitivity = NULL, eligible = FALSE,
        checks = data.frame(Check = "Anchor derivative calculation", Passed = FALSE,
          Detail = conditionMessage(e))))
    a[names(calculation)] <- calculation
  }
  if (is.null(V)) {
    a$eligible <- FALSE
    a$checks <- data.frame(Check = "Conditional joint information", Passed = FALSE,
      Detail = "Receiving-fit information is unavailable; no source-error propagation was performed.")
  }
  tab <- out$table
  tab$ConditionalSE <- tab$SE; tab$ConditionalLower <- tab$Lower
  tab$ConditionalUpper <- tab$Upper; tab$ConditionalStatus <- tab$Status
  tab$AnchorSourceSE <- NA_real_
  sizes <- build_param_sizes(fit$config); slices <- build_param_slices(sizes)
  if (population) {
    J <- matrix(0, nrow(tab), length(fit$opt$par))
    beta <- slices$beta
    J[seq_along(beta), beta] <- diag(length(beta))
    if (!is.null(out$contrasts))
      J[length(beta)+seq_len(nrow(out$contrasts)), beta] <- out$contrasts
    direct <- matrix(0, nrow(tab), nrow(a$anchors))
  } else {
    J <- out$target_jacobian
    direct <- matrix(0, nrow(tab), nrow(a$anchors))
    if (identical(out$settings$facet, fit$config$noncenter_facet)) {
      direct <- out$contrasts[, a$anchors$Level, drop = FALSE]
      if (identical(out$settings$scale, "standardized_population_sd"))
        direct <- direct / out$settings$population_sd
    }
  }
  selected <- matrix(NA_real_, nrow(tab), nrow(tab))
  propagated <- NULL
  if (isTRUE(a$eligible) && !is.null(V)) {
    D <- J %*% a$sensitivity + direct
    added <- D %*% a$covariance %*% t(D)
    selected <- symmetrize_matrix(J %*% V %*% t(J) + added)
    propagated <- symmetrize_matrix(V + a$sensitivity %*% a$covariance %*% t(a$sensitivity))
    if (any(!is.finite(propagated))) {
      a$eligible <- FALSE
      a$checks <- rbind(a$checks[, c("Check", "Passed")],
        data.frame(Check = "Finite propagated covariance", Passed = FALSE))
      selected[,] <- NA_real_; propagated <- NULL
    } else {
      tab$AnchorSourceSE <- sqrt(pmax(diag(added), 0))
      invalid_target <- !is.finite(diag(selected)) | !is.finite(tab$AnchorSourceSE)
      selected[invalid_target, ] <- selected[, invalid_target] <- NA_real_
      tab$AnchorSourceSE[invalid_target] <- NA_real_
      a$target_sensitivity <- D
    }
  }
  se <- sqrt(pmax(diag(selected), 0))
  fixed <- tab$ConditionalStatus == "fixed" & is.finite(se) & se == 0
  if (population) fixed[tab$Type == "residual_variance"] <- TRUE
  available <- isTRUE(a$eligible) & is.finite(se) & se > 0 &
    tab$ConditionalStatus %in% c("available", "fixed")
  critical <- qnorm((1-out$settings$level)/2, lower.tail = FALSE)
  tab$SE <- ifelse(available, se, ifelse(fixed & !population, 0, NA_real_))
  tab$Lower <- ifelse(available, tab$Estimate-critical*se, NA_real_)
  tab$Upper <- ifelse(available, tab$Estimate+critical*se, NA_real_)
  available <- available & is.finite(tab$Lower) & is.finite(tab$Upper)
  tab[!available & !fixed, c("SE", "Lower", "Upper")] <- NA_real_
  tab$Status <- ifelse(fixed, "fixed", ifelse(available, "available", "anchor_propagation_unavailable"))
  tab$CIEligible <- available
  tab$CoverageQualified <- FALSE
  note <- paste("First-order intervals include supplied anchor-source covariance under independent source and receiving samples.",
    "Source-covariance validity, anchor invariance and fixed population SD are assumptions; drift, dependence and covariance-estimation uncertainty are not corrected. Coverage is not guaranteed.")
  tab$InferenceReview <- ifelse(available, note, ifelse(fixed,
    "Fixed under the stated model and supplied covariance; no inferential interval.",
    "Source-error propagation is unavailable; finite points and conditional intervals remain separately recorded."))
  out$table <- tab
  if (population) {
    out$conditional_covariance <- out$covariance
    out$covariance <- if (is.null(propagated)) matrix(NA_real_, length(beta), length(beta)) else
      propagated[beta, beta, drop = FALSE]
    dimnames(out$covariance) <- dimnames(out$conditional_covariance)
    if (!is.null(out$contrasts)) out$contrast_covariance <-
      out$contrasts %*% out$covariance %*% t(out$contrasts)
    coefficient <- tab$Type != "residual_variance"
    for (pair in list(c("WorkingSE", "SE"), c("WorkingLower", "Lower"), c("WorkingUpper", "Upper")))
      out$table[coefficient, pair[1]] <- tab[coefficient, pair[2]]
    out$table$NaturalScaleStatus[coefficient] <- ifelse(available[coefficient], "finite", "unavailable")
  } else {
    dimnames(selected) <- dimnames(out$covariance)
    out$covariance <- selected
  }
  out$parameter_covariance <- propagated
  out$settings$uncertainty <- "independent_anchor_delta"
  out$settings$anchor_independence <- TRUE
  omitted_note <- mfrm_gpcm_anchor_note(fit)
  out$cautions <- unique(c(out$cautions[!out$cautions %in% omitted_note], note))
  out$table$InferenceCaution <- paste(out$cautions, collapse = " ")
  a$fingerprint <- mfrm_checkpoint_fingerprint(a, canonical = TRUE)
  out$anchor_uncertainty <- a
  out
}

mfrmr_anchor_uncertainty_tables <- function(x) {
  if (is.null(x$anchor_uncertainty)) return(list())
  mfrmr_anchor_uncertainty_validate_saved(x)
  a <- x$anchor_uncertainty
  list(anchor_source_covariance = data.frame(Level = rownames(a$covariance),
    a$covariance, check.names = FALSE, row.names = NULL),
    anchor_source_values = a$anchors, anchor_sensitivity_checks = a$checks)
}
