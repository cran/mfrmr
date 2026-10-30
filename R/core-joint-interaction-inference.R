# Expand the fitted zero-margin interaction before transforming the joint inverse.
mfrm_joint_interaction_targets <- function(fit, interaction, contrasts) {
  spec <- fit$config$interaction_specs[[interaction]]
  cells <- as.data.frame(fit$interactions$effects)
  labels <- make.unique(paste0(cells$FacetA, "=", cells$FacetA_Level,
    " | ", cells$FacetB, "=", cells$FacetB_Level))
  cells$Target <- labels
  if (is.null(contrasts)) {
    contrasts <- diag(length(labels))
    dimnames(contrasts) <- list(labels, labels)
  }
  if (!is.matrix(contrasts) || !is.numeric(contrasts) || is.complex(contrasts) ||
      !nrow(contrasts) || ncol(contrasts) != length(labels) || any(!is.finite(contrasts)) ||
      is.null(colnames(contrasts)) || anyDuplicated(colnames(contrasts)) ||
      !setequal(colnames(contrasts), labels) || is.null(rownames(contrasts)) ||
      anyNA(rownames(contrasts)) || any(!nzchar(trimws(rownames(contrasts)))) ||
      anyDuplicated(rownames(contrasts)) || any(rowSums(abs(contrasts)) == 0))
    stop("Interaction contrasts require distinct target row names, every cell Target as a column name, and finite nonzero rows.", call. = FALSE)
  contrasts <- contrasts[, labels, drop = FALSE]
  G <- vapply(seq_len(spec$n_params), function(j) {
    unit <- numeric(spec$n_params); unit[j] <- 1
    as.vector(expand_two_way_interaction(unit, spec))
  }, numeric(spec$n_cells))
  J <- matrix(0, nrow(contrasts), length(fit$opt$par),
    dimnames = list(rownames(contrasts), mfrm_checkpoint_parameter_names(build_param_sizes(fit$config))))
  J[, build_param_slices(build_param_sizes(fit$config))$interactions] <- contrasts %*% G
  list(contrasts = contrasts, jacobian = J, cells = cells,
    fixed = rowSums(abs(J)) == 0, constant = rep(0, nrow(contrasts)))
}
