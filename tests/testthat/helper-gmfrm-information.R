# Independent literal response equation and Louis observed information.
# Nodes/weights are supplied externally. No package response, parameter-expansion,
# score or Hessian function is called here; all location/step nuisance terms stay.
gmfrm_louis_reference <- function(spec, par, nodes, log_weights) {
  d <- spec$data; owners <- spec$slope_facets; maximum <- spec$max_score
  levels <- lapply(d[c(spec$person, owners)], function(x) sort(unique(as.character(x))))
  ids <- Map(function(x, lev) match(as.character(x), lev), d[names(levels)], levels)
  ni <- length(levels[[owners[1]]]); nr <- length(levels[[owners[2]]])
  centered <- rbind(diag(ni - 1L), rep(-1, ni - 1L))
  Z <- cbind(centered[ids[[owners[1]]], , drop = FALSE],
    diag(nr)[ids[[owners[2]]], , drop = FALSE])
  ng <- ncol(Z); nb <- ng + nr * (maximum - 1L); p <- nb + ng
  stopifnot(length(par) == p, nrow(nodes) == length(levels[[spec$person]]),
    identical(dim(nodes), dim(log_weights)))
  D <- lapply(0:maximum, function(k) {
    steps <- matrix(0, nrow(d), nr * (maximum - 1L))
    if (k > 0L && maximum > 1L) {
      contrast <- rbind(diag(maximum - 1L), rep(-1, maximum - 1L))
      for (r in seq_len(nr)) {
        cols <- (r - 1L) * (maximum - 1L) + seq_len(maximum - 1L)
        rows <- which(ids[[owners[2]]] == r)
        steps[rows, cols] <- rep(colSums(contrast[seq_len(k), , drop = FALSE]), each = length(rows))
      }
    }
    cbind(k * Z, steps)
  })
  b <- seq_len(nb); g <- nb + seq_len(ng)
  slope <- exp(drop(Z %*% par[g]))
  information <- matrix(0, p, p)
  scores <- matrix(0, nrow(nodes), p)
  log_marginal <- numeric(nrow(nodes))
  for (i in seq_len(nrow(nodes))) {
    rows <- which(ids[[spec$person]] == i)
    z <- Z[rows, , drop = FALSE]; a <- slope[rows]
    design <- lapply(D, function(x) x[rows, , drop = FALSE])
    offset <- vapply(design, function(x) drop(x %*% par[b]), numeric(length(rows)))
    at <- cbind(seq_along(rows), d[[spec$score]][rows] + 1L)
    complete <- lapply(seq_len(ncol(nodes)), function(q) {
      logits <- (matrix(nodes[i, q] * (0:maximum), length(rows), maximum + 1L,
        byrow = TRUE) - offset) * a
      high <- apply(logits, 1L, max)
      logp <- logits - high - log(rowSums(exp(logits - high)))
      probability <- exp(logp)
      residual <- -probability; residual[at] <- residual[at] + 1
      derivatives <- lapply(seq_along(design), function(k)
        cbind(-a * design[[k]], logits[, k] * z))
      mean <- Reduce(`+`, lapply(seq_along(design), function(k)
        probability[, k] * derivatives[[k]]))
      score <- colSums(Reduce(`+`, lapply(seq_along(design), function(k)
        residual[, k] * derivatives[[k]])))
      hessian <- Reduce(`+`, lapply(seq_along(design), function(k)
        crossprod(sqrt(probability[, k]) * derivatives[[k]]))) - crossprod(mean)
      cross <- crossprod(a * Reduce(`+`, lapply(seq_along(design), function(k)
        residual[, k] * design[[k]])), z)
      hessian[b, g] <- hessian[b, g] + cross
      hessian[g, b] <- hessian[g, b] + t(cross)
      hessian[g, g] <- hessian[g, g] - crossprod(rowSums(residual * logits) * z, z)
      list(logp = sum(logp[at]), score = score, information = hessian)
    })
    joint <- log_weights[i, ] + vapply(complete, `[[`, 0, "logp")
    high <- max(joint); log_marginal[i] <- high + log(sum(exp(joint - high)))
    posterior <- exp(joint - log_marginal[i])
    scores[i, ] <- Reduce(`+`, Map(function(x, w) w * x$score, complete, posterior))
    # I_obs = E(I_complete | y) - Var(U_complete | y), including cross blocks.
    information <- information + Reduce(`+`, Map(function(x, w)
      w * (x$information - tcrossprod(x$score)), complete, posterior)) +
      tcrossprod(scores[i, ])
  }
  list(information = information, scores = scores, log_marginal = log_marginal)
}
