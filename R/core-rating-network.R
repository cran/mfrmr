.mfrm_network_robustness <- function(edges, audit, panels) {
  persons <- unique(edges$Person)
  raters <- unique(edges$Rater)
  n <- length(persons) + length(raters)
  i <- match(edges$Person, persons)
  j <- length(persons) + match(edges$Rater, raters)
  degree <- tabulate(c(i, j), n)
  L <- Matrix::Diagonal(n, degree) -
    Matrix::sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = c(n, n), symmetric = TRUE)
  components <- audit$summary$Components
  lambda <- if (components > 1L) {
    0
  } else if (n <= 128L) {
    sort(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)[2L]
  } else if (requireNamespace("RSpectra", quietly = TRUE)) {
    tryCatch({
      # RSpectra dispatches on general compressed sparse matrices.
      sparse <- methods::as(L, "generalMatrix")
      values <- RSpectra::eigs_sym(sparse, k = 2L, which = "SA", opts = list(retvec = FALSE))$values
      if (length(values) == 2L) max(values) else NA_real_
    }, error = function(e) NA_real_)
  } else {
    NA_real_
  }

  status <- if (components > 1L) "disconnected" else if (is.finite(lambda)) "available" else "sparse_eigenvalue_unavailable"
  overview <- data.frame(Persons = length(persons), Raters = length(raters),
                         Edges = nrow(edges), Components = components,
                         FiedlerValue = if (is.finite(lambda)) max(lambda, 0) else NA_real_,
                         EdgeWeight = "one per observed person-rater pair", Status = status)
  resistance <- NULL
  if (!is.null(panels)) {
    if (!is.character(panels) || is.null(names(panels)) || anyNA(panels) ||
        any(!nzchar(panels)) || anyNA(names(panels)) || anyDuplicated(names(panels)) ||
        length(unique(panels)) < 2L || any(!names(panels) %in% raters)) {
      stop("rater_panels must map unique observed rater IDs to at least two panel labels.", call. = FALSE)
    }
    members <- split(names(panels), panels)
    pair <- utils::combn(names(members), 2L)
    rcomp <- audit$membership[raters]
    pcomp <- rcomp[match(edges$Rater[match(persons, edges$Person)], raters)]
    component <- c(pcomp, rcomp)
    keep <- which(!duplicated(component))
    free <- setdiff(seq_len(n), keep)
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(L[free, free]), LDL = FALSE)
    resistance <- do.call(rbind, lapply(seq_len(ncol(pair)), function(k) {
      a <- members[[pair[1L, k]]]
      b <- members[[pair[2L, k]]]
      contrast <- numeric(n)
      contrast[length(persons) + match(a, raters)] <- 1 / length(a)
      contrast[length(persons) + match(b, raters)] <- -1 / length(b)
      balanced <- all(abs(rowsum(contrast, component)) < 1e-12)
      value <- if (balanced)
        max(0, sum(contrast[free] * Matrix::solve(factor, contrast[free]))) else Inf
      data.frame(PanelA = pair[1L, k], PanelB = pair[2L, k], RatersA = length(a),
                 RatersB = length(b), EffectiveResistance = value,
                 Status = if (balanced) "available" else "disconnected_contrast")
    }))
  }
  list(overview = overview, panel_resistance = resistance)
}
