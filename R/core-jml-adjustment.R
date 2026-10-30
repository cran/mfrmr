# Internal observed-data equation for shared-owner GPCM JML score adjustment.
# U_k = (I-P_beta)^k U_0, with MLE plug-in P, not posterior averaging.
# The exact owner-total recurrence generalizes the frozen research reference
# in inst/validation/jml-total-expectation-20260927.R. It does not establish
# removal of incidental-parameter bias, choose an order, or fit a public model.
mfrm_jml_adjustment_problem <- function(data, person, facets, score, owner,
    rating_max, max_states = 5000L) {
  names_used <- c(person, facets, score)
  if (!is.character(person) || length(person) != 1L ||
      !is.character(score) || length(score) != 1L ||
      !is.character(facets) || !length(facets) || anyNA(names_used) ||
      any(!nzchar(names_used)) || anyDuplicated(names_used) ||
      !is.character(owner) || length(owner) != 1L || is.na(owner) || !owner %in% facets ||
      !is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      !all(names_used %in% names(data))) stop("Declare distinct person, facet and score columns and one shared slope/step owner.")
  if (!is.numeric(rating_max) || length(rating_max) != 1L || !is.finite(rating_max) ||
      rating_max < 1 || rating_max != floor(rating_max) ||
      !is.numeric(max_states) || length(max_states) != 1L || !is.finite(max_states) ||
      max_states < 1 || max_states != floor(max_states)) stop("Invalid category range or owner-total state limit.")
  y <- data[[score]]
  if (!is.numeric(y) || is.complex(y) || anyNA(y) || any(!is.finite(y)) ||
      any(y < 0 | y > rating_max | y != floor(y)) || anyNA(data[names_used]))
    stop("Supply observed integer scores on the declared zero-based scale; do not impute unassigned rows.")
  ids <- lapply(data[c(person, facets)], as.character)
  if (any(vapply(ids, function(x) any(!nzchar(x) | x != trimws(x)), logical(1))))
    stop("Person and facet IDs must be nonempty and have no surrounding whitespace.")
  levels <- lapply(ids, function(x) sort(unique(x), method = "radix"))
  if (any(lengths(levels[facets]) < 2L)) stop("Each declared facet needs at least two observed levels.")
  ids <- Map(match, ids, levels)
  cell_ids <- as.data.frame(unname(ids[facets]), check.names = FALSE)
  names(cell_ids) <- facets
  cells <- unique(cell_ids)
  cells <- cells[do.call(order, unname(cells)), , drop = FALSE]
  key <- function(x) do.call(paste, c(unname(x), sep = ":")) # Integer indices, never user labels.
  cell <- match(key(cell_ids), key(cells))
  np <- length(levels[[person]]); nc <- nrow(cells); K <- as.integer(rating_max + 1L)
  if (!is.finite(K) || 8 * as.double(np) * nc * K > getOption("mfrmr.max_information_bytes", 256 * 1024^2))
    stop("Observed-count workspace exceeds mfrmr.max_information_bytes.")
  counts <- lapply(seq_len(nc), function(j) {
    selected <- which(cell == j)
    matrix(tabulate(ids[[person]][selected] + np * y[selected], nbins = np * K), np, K)
  })
  exposure <- vapply(counts, rowSums, numeric(np))
  exposure <- matrix(exposure, np, nc)
  roster_key <- apply(exposure, 1L, paste, collapse = ":")
  roster <- match(roster_key, unique(roster_key))
  owner_index <- cells[[owner]]; ng <- length(levels[[owner]])
  location <- do.call(cbind, lapply(facets, function(f)
    sum_zero_jacobian(length(levels[[f]]))[cells[[f]], , drop = FALSE]))
  steps <- sum_zero_jacobian(K - 1L)
  ns <- ng * (K - 2L)
  linear <- lapply(0:rating_max, function(k) {
    out <- matrix(0, nc, ns)
    if (k > 0L && K > 2L) for (g in seq_len(ng)) {
      cols <- (g - 1L) * (K - 2L) + seq_len(K - 2L)
      out[owner_index == g, cols] <- rep(-colSums(steps[seq_len(k), , drop = FALSE]),
        each = sum(owner_index == g))
    }
    cbind(-k * location, out)
  })
  slopes <- sum_zero_jacobian(ng)[owner_index, , drop = FALSE]
  parameter <- do.call(rbind, c(lapply(facets, function(f)
    data.frame(Type = "location", Facet = f, Level = head(levels[[f]], -1L), Step = NA_integer_)),
    if (K > 2L) list(data.frame(Type = "step", Facet = owner,
      Level = rep(levels[[owner]], each = K - 2L), Step = rep(seq_len(K - 2L), ng))),
    list(data.frame(Type = "log_slope", Facet = owner, Level = head(levels[[owner]], -1L), Step = NA_integer_))))
  rownames(parameter) <- NULL
  groups <- split(seq_len(np), roster)
  problems <- lapply(groups, function(rows) {
    used <- which(exposure[rows[1L], ] > 0)
    design <- list(linear = lapply(linear, function(x) x[used, , drop = FALSE]),
      slopes = slopes[used, , drop = FALSE], owner = owner_index[used], owners = ng, categories = 0:rating_max)
    mfrm_jml_adjustment_roster(design, exposure[rows[1L], used],
      lapply(counts[used], function(x) x[rows, , drop = FALSE]), max_states)
  })
  evaluate <- function(beta, order = 1L) {
    if (!is.numeric(beta) || is.complex(beta) || length(beta) != nrow(parameter) ||
        any(!is.finite(beta)) || !is.numeric(order) || length(order) != 1L ||
        !is.finite(order) || order < 0 || order != floor(order)) stop("Invalid structural coordinates or correction order.")
    U <- matrix(0, np, length(beta)); q <- theta <- numeric(np); residual <- 0
    for (g in seq_along(problems)) {
      z <- problems[[g]](beta, order)
      U[groups[[g]], ] <- z$value; q[groups[[g]]] <- z$q; theta[groups[[g]]] <- z$theta
      residual <- max(residual, z$root_residual)
    }
    list(value = U, q = q, theta = theta, root_residual = residual,
      order = order, equation = "finite_MLE_plugin_profile_score_adjustment")
  }
  list(parameters = parameter, persons = levels[[person]], levels = levels,
    cells = cells, exposure = exposure, roster = roster, counts = counts,
    total_states = vapply(problems, attr, integer(1), "total_states"),
    evaluate = evaluate, mean_score = function(beta, order = 1L) colMeans(evaluate(beta, order)$value),
    specification = list(model = "GPCM", owner = owner, person = person, facets = facets,
      score = score, rating_max = rating_max, identification = "centered_locations_steps_and_log_slopes",
      sampling_unit = "independent_Person", max_states = max_states))
}

mfrm_jml_adjustment_roster <- function(design, exposure, counts, max_states) {
  nc <- length(exposure); K <- length(design$categories); ng <- design$owners
  nl <- ncol(design$linear[[1L]]); na <- ncol(design$slopes); p <- nl + na
  maximum <- vapply(seq_len(ng), function(g) (K - 1L) * sum(exposure[design$owner == g]), numeric(1))
  if (sum(log(maximum + 1)) > log(max_states) + 1e-12)
    stop("Owner-total space exceeds max_states; no states were omitted.")
  totals <- as.matrix(expand.grid(lapply(maximum, function(n) 0:n)))
  G <- nrow(totals); extreme <- rowSums(totals) %in% c(0, sum(maximum)); ix <- which(!extreme)
  N <- nrow(counts[[1L]])
  observed <- matrix(0, N, ng)
  for (j in seq_len(nc)) observed[, design$owner[j]] <- observed[, design$owner[j]] + drop(counts[[j]] %*% design$categories)
  stride <- cumprod(c(1, head(maximum + 1, -1L)))
  group <- as.integer(1 + observed %*% stride)

  evaluate <- function(beta, order) {
    a <- exp(drop(design$slopes %*% tail(beta, na)))
    if (any(!is.finite(a) | a <= 0)) stop("Nonfinite component slope.")
    base <- vapply(design$linear, function(x) drop(x %*% head(beta, nl)), numeric(nc))
    base <- matrix(base, nc, K)
    ag <- rep(1, ng)
    for (j in seq_len(nc)) ag[design$owner[j]] <- a[j]
    logits <- function(theta, j) a[j] * (outer(theta, design$categories) +
      matrix(base[j, ], length(theta), K, byrow = TRUE))
    log_probability <- function(theta, j) {
      z <- logits(theta, j)
      z <- z - if (all(is.finite(z))) row_max_fast(z) else apply(z, 1L, max)
      z - log(rowSums(exp(z)))
    }
    # Positive forward convolution; no subtraction or probability pruning.
    table <- lapply(seq_len(ng), function(g) {
      probability <- 1; moment <- matrix(0, 1L, nc * K)
      for (j in which(design$owner == g)) for (repeat_id in seq_len(exposure[j])) {
        w <- as.numeric(exp(log_probability(0, j))); M <- length(probability)
        next_p <- numeric(M + K - 1L); next_m <- matrix(0, M + K - 1L, nc * K)
        for (k in seq_len(K)) {
          rows <- seq_len(M) + k - 1L; column <- K * (j - 1L) + k
          next_p[rows] <- next_p[rows] + w[k] * probability
          contribution <- moment; contribution[, column] <- contribution[, column] + probability
          next_m[rows, ] <- next_m[rows, , drop = FALSE] + w[k] * contribution
        }
        probability <- next_p; moment <- next_m
      }
      if (any(!is.finite(probability) | probability <= 0)) stop("Underflow in owner-total probabilities; no states were omitted.")
      list(probability = probability, counts = moment / probability)
    })
    expected_score <- function(theta) {
      value <- numeric(length(theta))
      for (j in seq_len(nc)) value <- value + exposure[j] * a[j] *
        drop(exp(log_probability(theta, j)) %*% design$categories)
      value
    }
    theta <- numeric(G); residual <- 0
    if (length(ix)) {
      target <- drop(totals[ix, , drop = FALSE] %*% ag); bound <- 1
      while (any(expected_score(-bound) >= target) || any(expected_score(bound) <= target)) {
        bound <- 2 * bound
        if (bound > 1e6) stop("Unresolved finite Person profile.")
      }
      lo <- rep(-bound, length(ix)); hi <- rep(bound, length(ix))
      for (iteration in seq_len(100L)) {
        mid <- (lo + hi) / 2; delta <- expected_score(mid) - target
        lo[delta < 0] <- mid[delta < 0]; hi[delta >= 0] <- mid[delta >= 0]
        if (max(hi - lo) < 1e-12) break
      }
      if (max(hi - lo) >= 1e-12) stop("Person profile failed bisection tolerance.")
      theta[ix] <- (lo + hi) / 2; residual <- max(abs(expected_score(theta[ix]) - target))
    }
    raw <- function(observed_counts, theta, boundary) {
      n <- length(theta); U <- matrix(0, n, p); q <- numeric(n)
      for (j in seq_len(nc)) {
        logp <- log_probability(theta, j); z <- logits(theta, j)
        r <- exposure[j] * exp(logp) - observed_counts[[j]]
        q <- q - rowSums(observed_counts[[j]] * logp)
        for (k in seq_len(K)) U[, seq_len(nl)] <- U[, seq_len(nl), drop = FALSE] +
          tcrossprod(r[, k], a[j] * design$linear[[k]][j, ])
        U[, nl + seq_len(na)] <- U[, nl + seq_len(na), drop = FALSE] +
          tcrossprod(rowSums(r * z), design$slopes[j, ])
      }
      U[boundary, ] <- 0; q[boundary] <- 0
      list(value = U, q = q)
    }
    expected_counts <- lapply(seq_len(nc), function(j)
      table[[design$owner[j]]]$counts[totals[, design$owner[j]] + 1L,
        K * (j - 1L) + seq_len(K), drop = FALSE])
    mean_raw <- raw(expected_counts, theta, extreme)$value
    actual <- raw(counts, theta[group], extreme[group])
    blocks <- if (length(ix)) lapply(seq.int(1L, length(ix), by = 16L),
      function(begin) ix[begin:min(length(ix), begin + 15L)]) else list()
    transition_probability <- function(ii) {
      probability <- matrix(1, G, length(ii))
      for (g in seq_len(ng)) {
        z <- outer(0:maximum[g], ag[g] * theta[ii]) + log(table[[g]]$probability)
        z <- sweep(z, 2L, apply(z, 2L, max)); mass <- exp(z)
        mass <- sweep(mass, 2L, colSums(mass), "/")
        probability <- probability * mass[totals[, g] + 1L, , drop = FALSE]
      }
      if (max(abs(colSums(probability) - 1)) > 1e-10) stop("Owner-total transition is not normalized.")
      probability
    }
    # Reuse only within this beta evaluation. Reserve three quarters of the
    # workspace allowance for counts, scores and temporary calculation buffers.
    probabilities <- if (order > 1L && 8 * as.double(G) * length(ix) <=
        getOption("mfrmr.max_information_bytes", 256 * 1024^2) / 4)
      lapply(blocks, transition_probability) else NULL
    apply_transition <- function(values) {
      answer <- matrix(0, G, p)
      for (block in seq_along(blocks)) {
        ii <- blocks[[block]]
        probability <- if (is.null(probabilities)) transition_probability(ii) else probabilities[[block]]
        answer[ii, ] <- crossprod(probability, values)
      }
      answer
    }
    adjustment <- matrix(0, G, p)
    for (k in seq_len(order)) adjustment <- adjustment + apply_transition(mean_raw - adjustment)
    actual$value <- actual$value - adjustment[group, , drop = FALSE]
    theta[rowSums(totals) == 0] <- -Inf
    theta[rowSums(totals) == sum(maximum)] <- Inf
    c(actual, list(theta = theta[group], root_residual = residual))
  }
  attr(evaluate, "total_states") <- as.integer(G)
  evaluate
}

# Full nonsymmetric equation derivative; conditional means never replace the
# actual Person contributions. Inference concerns the equation's root, which
# may still differ from the generating structural parameter.
# Within-roster centering describes independent new Persons from each roster's
# population, including ability composition. It does not condition on every
# Person's true ability. This is the empirical asymptotic meat, without an
# n_g/(n_g-1) correction; two Persons is only a computational minimum.
mfrm_jml_adjustment_covariance <- function(problem, beta, order,
    sampling = c("fixed_rosters", "random_rosters")) {
  sampling <- match.arg(sampling)
  rosters <- if (sampling == "fixed_rosters") split(seq_along(problem$roster), problem$roster)
  if (any(lengths(rosters) < 2L))
    stop("Fixed-roster covariance needs at least two Persons per observed roster.")
  fn <- function(b) problem$mean_score(b, order)
  coarse <- mfrmr_numeric_transformation_jacobian(fn, beta, relative_step = 1e-4)
  fine <- mfrmr_numeric_transformation_jacobian(fn, beta, relative_step = 5e-5)
  U <- problem$evaluate(beta, order)$value
  N <- nrow(U); p <- ncol(U); A <- fine$jacobian
  centered <- sweep(U, 2L, colMeans(U))
  for (rows in rosters) {
    centered[rows, ] <- sweep(U[rows, , drop = FALSE], 2L, colMeans(U[rows, , drop = FALSE]))
  }
  if (!coarse$valid || !fine$valid)
    stop("Covariance unavailable: the complete equation derivative could not be evaluated.")
  if (N <= p)
    stop("Covariance unavailable: the number of independent Persons must exceed the structural dimension.")
  if (max(abs(colMeans(U))) > 1e-7)
    stop("Covariance unavailable: the adjusted mean equation is not solved to the declared tolerance.")
  if (max(abs(coarse$jacobian - A)) / max(1, max(abs(A))) > 1e-5)
    stop("Covariance unavailable: the complete equation derivative changes materially with the numerical step.")
  if (min(svd(A, nu = 0, nv = 0)$d) <= 1e-6)
    stop("Covariance unavailable: the full equation Jacobian is singular or below the declared conditioning threshold.")
  if (qr(centered)$rank < p)
    stop("Covariance unavailable: the actual Person score contributions do not span the structural coordinates.")
  inverse <- solve(A)
  influence <- -centered %*% t(inverse)
  V <- crossprod(influence) / N^2
  list(vcov = V, jacobian = A, meat = crossprod(centered) / N,
    scores = U, influence = influence, order = order, sampling = sampling,
    target = "Local variation around the adjusted-equation root; residual structural bias is not removed.")
}
