# Literal two-family response equation and continuous integrals. This oracle
# uses no package response, parameter-expansion, mode or quadrature functions.
gmfrm_response_reference <- function(fit, rows) {
  spec <- fit$gmfrm$specification; data <- spec$data; owners <- spec$slope_facets
  person <- as.character(data[[spec$person]])
  stopifnot(length(unique(person[rows])) == 1L)
  indices <- which(person == person[rows[1L]])
  location <- Reduce(`+`, lapply(owners, function(f) {
    tab <- fit$facets$others[fit$facets$others$Facet == f, ]
    tab$Estimate[match(as.character(data[[f]][indices]), tab$Level)]
  }))
  slope <- Reduce(`*`, lapply(owners, function(f) {
    tab <- fit$slopes[fit$slopes$SlopeOwner == f, ]
    tab$Estimate[match(as.character(data[[f]][indices]), tab$SlopeFacet)]
  }))
  categories <- 0:spec$max_score
  cumulative <- t(vapply(as.character(data[[owners[2]]][indices]), function(level) {
    tab <- fit$steps[fit$steps$StepFacet == level, ]
    c(0, cumsum(tab$Estimate[order(tab$Step)]))
  }, numeric(length(categories))))
  observed <- cbind(seq_along(indices), data[[spec$score]][indices] + 1L)
  logp <- function(theta) {
    logits <- (outer(theta - location, categories) - cumulative) * slope
    high <- apply(logits, 1L, max)
    logits - high - log(rowSums(exp(logits - high)))
  }
  log_joint <- function(theta) sum(logp(theta)[observed]) + dnorm(theta, log = TRUE)
  score <- function(theta) sum(slope * (data[[spec$score]][indices] -
    drop(exp(logp(theta)) %*% categories))) - theta
  mode <- uniroot(score, c(-1, 1), extendInt = "downX", tol = 1e-12)$root
  probability <- exp(logp(mode)); mean <- drop(probability %*% categories)
  scale <- 1 / sqrt(1 + sum(slope^2 * rowSums(probability *
    (rep(categories, each = length(indices)) - mean)^2)))
  height <- log_joint(mode)
  density <- function(z) exp(log_joint(mode + scale * z) - height) * scale
  integrate_one <- function(fun) {
    answer <- integrate(Vectorize(fun), -Inf, Inf, rel.tol = 1e-10,
      abs.tol = 1e-12, subdivisions = 1000L)
    stopifnot(identical(answer$message, "OK"))
    answer$value
  }
  denominator <- integrate_one(density)
  probabilities <- t(vapply(rows, function(row) vapply(categories, function(k)
    integrate_one(function(z) density(z) * exp(logp(mode + scale * z)[match(row, indices), k + 1L])) /
      denominator, 0), numeric(length(categories))))
  expected <- drop(probabilities %*% categories)
  conditional_variance <- vapply(rows, function(row) integrate_one(function(z) {
    p <- exp(logp(mode + scale * z)[match(row, indices), ])
    mu <- sum(categories * p)
    density(z) * sum((categories - mu)^2 * p)
  }) / denominator, 0)
  eap <- mode + scale * integrate_one(function(z) z * density(z)) / denominator
  list(probabilities = probabilities, mean = expected,
    variance = rowSums(probabilities * (rep(categories, each = length(rows)) - expected)^2),
    conditional_variance = conditional_variance,
    plugin = exp(logp(eap)[match(rows, indices), , drop = FALSE]),
    mode = mode, local_sd = scale)
}

# Literal category logits and continuous scalar integration, independent of
# package indices, response kernels, adaptive nodes and interval inversion.
gmfrm_scoring_reference <- function(fit, data, level = .95) {
  owners <- fit$config$slope_facet
  location <- Reduce(`+`, lapply(owners, function(owner) {
    tab <- fit$facets$others[fit$facets$others$Facet == owner, ]
    tab$Estimate[match(as.character(data[[owner]]), tab$Level)]
  }))
  slope <- Reduce(`*`, lapply(owners, function(owner) {
    tab <- fit$slopes[fit$slopes$SlopeOwner == owner, ]
    tab$OptimizerEstimate[match(as.character(data[[owner]]), tab$SlopeFacet)]
  }))
  k <- 0:(fit$config$n_cat - 1L)
  cumulative <- t(vapply(as.character(data[[fit$config$step_facet]]), function(id) {
    tab <- fit$steps[fit$steps$StepFacet == id, ]
    c(0, cumsum(tab$Estimate[order(tab$Step)]))
  }, numeric(length(k))))
  observed <- cbind(seq_len(nrow(data)), data$Score + 1L)
  joint <- function(theta) {
    logits <- (outer(theta - location, k) - cumulative) * slope
    high <- apply(logits, 1L, max)
    sum(logits[observed] - high - log(rowSums(exp(logits - high)))) + dnorm(theta, log = TRUE)
  }
  mode <- optimize(joint, c(-40, 40), maximum = TRUE, tol = 1e-10)$maximum
  stopifnot(abs(mode) < 39)
  density <- Vectorize(function(z) exp(joint(mode + z) - joint(mode)))
  integral <- function(fun, lower = -Inf, upper = Inf) {
    integrate(fun, lower, upper, rel.tol = 1e-10, abs.tol = 1e-12,
      subdivisions = 1000L)$value
  }
  mass <- integral(density)
  center <- integral(function(z) z * density(z)) / mass
  variance <- integral(function(z) (z - center)^2 * density(z)) / mass
  left <- integral(density, upper = 0)
  quantile <- function(p) mode + uniroot(function(z)
    (if (z <= 0) integral(density, upper = z) else left + integral(density, 0, z)) /
      mass - p, c(-30, 30), tol = 1e-9)$root
  c(Estimate = mode + center, SD = sqrt(variance),
    Lower = quantile((1 - level) / 2), Upper = quantile((1 + level) / 2))
}

gmfrm_scoring_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      saved <- readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))
      p31 <- mfrm_gmfrm_problem(saved$data, 2L, gauss_hermite_normal(31L))
      f31 <- mfrm_gmfrm_fit_result(p31, mfrm_gmfrm_em(p31, start = saved$parameters, maxit = 1L))
      p61 <- mfrm_gmfrm_problem(saved$data, 2L, gauss_hermite_normal(61L))
      f61 <- mfrm_gmfrm_fit_result(p61, mfrm_gmfrm_em(p61, start = saved$parameters, maxit = 100L))
      d <- saved$data[saved$data$Person %in% unique(saved$data$Person)[1:3], ]
      d$Person <- paste0("NEW", match(d$Person, unique(d$Person)))
      d$Score[d$Person == "NEW1"] <- 0L
      d$Score[d$Person == "NEW2"] <- 2L
      d <- d[!(d$Person == "NEW3" & d$Task == levels(d$Task)[1]), ]
      d <- d[rev(seq_len(nrow(d))), ]
      cached <<- list(fit = f61, coarse = f31, data = d)
    }
    cached
  }
})
