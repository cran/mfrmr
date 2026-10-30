#' Principal components of numeric external features
#'
#' Express numeric attributes as principal components: new variables ordered by
#' how much variation they describe. For example, summarize rater experience,
#' training and workload together. Choose `components` when you want to retain
#' fewer variables; the default retains all nonzero components. PCA does not
#' estimate rater quality, ability or groups. Identifiers and omitted rows are
#' retained.
#' @param x For `mfrm_pca()`: A table reviewed with [mfrm_features()]. For [mfrm_cluster_kmeans()],
#'   a result from `mfrm_pca()` is also accepted.
#'
#' For `plot()`: A result from [mfrm_pca()].
#' @param components For `mfrm_pca()`: Number of leading principal components to retain. `NULL`
#'   retains the numerical rank, using `sqrt(.Machine$double.eps)` relative to
#'   the largest singular value. A requested count cannot exceed that rank.
#'
#' For `plot()`: For scores, two distinct retained component numbers.
#'   For loadings, one retained component number. Defaults to `c(1, 2)` for
#'   scores and `1` for loadings. Unavailable axes cause an error.
#' @param scale Divide centered numeric features by their sample standard
#'   deviations (`TRUE`, default). With `FALSE`, distances retain the original
#'   feature units. Numeric coding must have a meaningful Euclidean interpretation;
#'   factors, ordered factors, characters and logical features are refused.
#' @inheritParams mfrm_cluster
#' @return An `mfrm_pca` object with ID-aligned `scores` (omitted rows are `NA`),
#'   `loadings`, a `variance` table (variance, proportion, cumulative proportion
#'   and retained status), `transformation` (centers, divisors and weight factors),
#'   original `feature_data` and `settings`. `summary()` returns `variance`.
#' @details Uses [stats::prcomp()] on centered, optionally standardized features,
#'   multiplied by `sqrt(weights / max(weights))`. Thus weights apply to
#'   squared Euclidean distances. Weight names
#'   identify original features. Constant features, nonfinite transformations
#'   and numerically unrepresentable variance are refused. Missingness stops by
#'   default; explicit omission retains every excluded ID and its original reason.
#'
#'   PCA finds directions of feature variance, not clusters or latent MFRM
#'   abilities. All nonzero components preserve distances in the transformed
#'   feature space, within numerical precision. Retaining fewer components
#'   changes the clustering objective and can remove informative directions.
#'   High explained variance does not establish valid groups. Signs of loadings
#'   and scores are arbitrary; tied eigenvalues can also rotate their subspace.
#'   Loadings are eigenvector coefficients in the transformed feature space,
#'   not correlations in the original units.
#'
#'   Use the same fitted PCA to review scores, loadings and groups. Passing its
#'   result to [mfrm_cluster_kmeans()] uses the retained scores without further
#'   standardization or whitening. PCA is not a required step for k-means.
#'   Group comparison and imputation sensitivity retain their descriptive scope.
#' @seealso [mfrm_cluster_kmeans()], [mfrm_cluster_compare()], [plot.mfrm_pca()],
#'   [stats::prcomp()]
#' @examples raters <- data.frame(Rater = paste0("R", 1:8),
#'   ExperienceYears = c(1, 2, 3, 4, 10, 12, 13, 15),
#'   WorkshopHours = c(4, 12, 8, 20, 12, 28, 16, 32),
#'   AnnualRatings = c(80, 120, 90, 160, 280, 350, 300, 400))
#' features <- mfrm_features(raters, "Rater", names(raters)[-1])
#' pca <- mfrm_pca(features, components = 2)
#' summary(pca)
#' pca$loadings
#' groups <- mfrm_cluster_kmeans(pca, k = 2, seed = 17, silhouette = FALSE)
#' plot(pca, type = "scores", groups = groups)
#' plot(groups, type = "profile", feature = "ExperienceYears")
#' @export
#' @order 1
mfrm_pca <- function(x, components = NULL, scale = TRUE, weights = NULL,
                     missing = c("error", "omit")) {
  prepared <- prepare_numeric_features(x, scale, weights, match.arg(missing))
  z <- prepared$coordinates
  if (!is.null(components)) {
    check_feature_integer(components, "components", 1, min(nrow(z) - 1L, ncol(z)))
  }
  fit <- stats::prcomp(z, center = FALSE, scale. = FALSE, rank. = components)
  rank <- sum(fit$sdev > sqrt(.Machine$double.eps) * fit$sdev[1L])
  if (rank < 1L || any(!is.finite(fit$sdev)) ||
      !is.finite(sum(fit$sdev^2)) || sum(fit$sdev^2) <= 0) {
    stop("PCA variance is not numerically representable; rescale the features.", call. = FALSE)
  }
  if (is.null(components)) components <- rank
  check_feature_integer(components, "components", 1, rank)
  used <- seq_len(as.integer(components))
  scores <- matrix(NA_real_, nrow(prepared$x$data), length(used),
    dimnames = list(NULL, paste0("PC", used)))
  scores[prepared$keep, ] <- fit$x[, used, drop = FALSE]
  variance <- fit$sdev^2
  proportion <- variance / sum(variance)
  out <- list(scores = data.frame(ID = prepared$x$row_summary$ID, scores,
      check.names = FALSE), loadings = fit$rotation[, used, drop = FALSE],
    variance = data.frame(Component = paste0("PC", seq_along(variance)),
      Variance = variance, Proportion = proportion, Cumulative = cumsum(proportion),
      Retained = seq_along(variance) %in% used),
    transformation = prepared$transformation, feature_data = prepared$x,
    settings = list(components = as.integer(components), rank = rank,
      scale = scale, weights = prepared$weights, missing = match.arg(missing),
      included = sum(prepared$keep), excluded = sum(!prepared$keep)))
  class(out) <- "mfrm_pca"
  out
}

#' @rdname mfrm_pca
#' @param object An object returned by `mfrm_pca()`.
#' @param ... For `print()`: Reserved for method compatibility.
#'   For `plot()`: Reserved for future use; additional arguments are rejected.
#' @export
print.mfrm_pca <- function(x, ...) {
  cat("Principal components of numeric external features\n")
  cat(x$settings$included, "entities included;", x$settings$excluded, "excluded;",
    x$settings$components, "components retained\n")
  print(x$variance, row.names = FALSE)
  cat("Components describe feature variation, not ability, rater quality or validated groups.\n")
  invisible(x)
}

#' @rdname mfrm_pca
#' @export
summary.mfrm_pca <- function(object, ...) object$variance

#' K-means groups from numeric features or retained principal components
#'
#' Divide people, raters or tasks into k groups with similar numeric attributes.
#' The algorithm minimizes squared distances within groups, using the numeric
#' features or retained principal components supplied by you. The groups are
#' descriptive; they do not measure ability or rater quality.
#'
#' @inheritParams mfrm_pca
#' @inheritParams mfrm_cluster
#' @param nstart Number of random starts; default 25.
#' @param iter.max Maximum iterations per start; default 100.
#' @param seed Nonnegative integer seed. The default is 1. The caller's random
#'   number state is restored; the RNG kind is retained in `settings`.
#' @param silhouette Compute Euclidean silhouettes (`TRUE`, default) using the
#'   optional `cluster` package. This requires pairwise distances and is limited
#'   to 5,000 included entities. Explicit `FALSE` avoids that quadratic allocation;
#'   silhouettes then remain unavailable, not zero. Other memory/time limits
#'   still depend on the dimensions and numerical workload.
#' @return An `mfrm_clusters` object with memberships, original-unit profiles,
#'   cluster sizes and silhouettes, `centers` in the fitted space, per-cluster
#'   `withinss`, `totss`, `tot.withinss`, the input PCA when used, original
#'   feature data, transformation and settings. K-means has no medoids.
#'   Numeric profiles use original feature units and sample SDs: SD is `NA`
#'   for a one-entity group and zero for identical values in a larger group.
#'   A singleton's silhouette is zero by convention, not evidence of separation.
#' @details Uses [stats::kmeans()] with Hartigan-Wong updates and the stated
#'   random starts. Warnings, nonconvergence and nonfinite sums of squares stop
#'   the analysis; no failed run is silently presented as a usable result.
#'   Multiple starts reduce sensitivity to initialization without guaranteeing
#'   a global optimum. Review alternative seeds separately when needed. Labels
#'   are arbitrary and results can depend on input order, ties and RNG kind.
#'
#'   Numeric preprocessing matches [mfrm_pca()]. With a PCA input, its retained
#'   scores are used exactly: do not respecify `weights`, `scale` or `missing`.
#'   Full-rank PCA preserves the distance objective; truncation changes it.
#'   No component whitening is performed. The centers of a truncated PCA are
#'   coordinates in that subspace; `profiles` describe the original attributes.
#'   Categorical codes must not be converted to numeric solely to pass this
#'   interface. Use [mfrm_cluster_pam()] for Gower/PAM mixed-feature groups.
#'
#'   Compare fitted partitions using [mfrm_cluster_compare()]. Silhouettes
#'   calculated in different geometries are not directly comparable evidence
#'   of which feature set, scale or component count is correct. Grouping does
#'   not propagate measurement uncertainty or establish rater quality.
#' @section Reading clustering metrics:
#' `withinss` is the sum of squared Euclidean distances from members to their
#' cluster center; `tot.withinss` sums it across clusters. `totss` uses the
#' grand center instead. The descriptive ratio 1 - tot.withinss/totss is the
#' fraction of fitted-space scatter between groups, not predictive accuracy,
#' latent reliability or a test of the chosen k. Increasing k usually improves
#' this ratio even without meaningful groups. These sums have squared units
#' of the transformed coordinates (SD-scaled/weighted features or retained PCs),
#' not necessarily original feature units. Compare them only for the same rows
#' and geometry. For silhouette a(i), b(i), range, singleton conventions and
#' profile denominators, see the metric section of [mfrm_cluster_pam()].
#' @seealso [mfrm_pca()], [mfrm_cluster_pam()], [mfrm_cluster_imputed()],
#'   [mfrm_cluster_compare()], [plot.mfrm_clusters()], [stats::kmeans()]
#' @examples
#' raters <- data.frame(Rater = paste0("R", 1:8),
#'   ExperienceYears = c(1, 2, 3, 4, 10, 12, 13, 15),
#'   WorkshopHours = c(4, 12, 8, 20, 12, 28, 16, 32),
#'   AnnualRatings = c(80, 120, 90, 160, 280, 350, 300, 400))
#' features <- mfrm_features(raters, "Rater", names(raters)[-1])
#' direct <- mfrm_cluster_kmeans(features, 2, seed = 17, silhouette = FALSE)
#' reduced <- mfrm_cluster_kmeans(mfrm_pca(features, components = 2), 2,
#'   seed = 17, silhouette = FALSE)
#' summary(mfrm_cluster_compare(list(Features = direct, TwoPCs = reduced)))
#' direct$centers
#' direct$profiles$numeric
#' @export
mfrm_cluster_kmeans <- function(x, k, weights = NULL,
                                missing = c("error", "omit"), scale = TRUE,
                                nstart = 25, iter.max = 100, seed = 1,
                                silhouette = TRUE) {
  pca <- NULL
  if (inherits(x, "mfrm_pca")) {
    if (!is.null(weights) || !missing(scale) || !missing(missing)) {
      stop("Set weights, scale and missing when fitting mfrm_pca(), not on its k-means input.", call. = FALSE)
    }
    pca <- x
    prepared <- prepare_numeric_features(x$feature_data, x$settings$scale,
      x$settings$weights, x$settings$missing)
    if (!identical(x$scores$ID, prepared$x$row_summary$ID) ||
        ncol(x$scores) != x$settings$components + 1L ||
        !identical(unname(is.na(as.matrix(x$scores[-1]))),
          matrix(!prepared$keep, nrow(x$scores), x$settings$components))) {
      stop("PCA scores must retain their original IDs and omission pattern.", call. = FALSE)
    }
    z <- as.matrix(x$scores[prepared$keep, -1, drop = FALSE])
    missing <- x$settings$missing
  } else {
    missing <- match.arg(missing)
    prepared <- prepare_numeric_features(x, scale, weights, missing)
    z <- prepared$coordinates
  }
  if (any(!is.finite(z))) stop("K-means coordinates must be finite.", call. = FALSE)
  n <- nrow(z)
  check_feature_integer(k, "k", 2, n - 1L)
  check_feature_integer(nstart, "nstart", 1)
  check_feature_integer(iter.max, "iter.max", 1)
  check_feature_integer(seed, "seed", 0)
  if (k > nrow(unique(z))) stop("`k` exceeds the number of distinct fitted profiles.", call. = FALSE)
  if (!is.logical(silhouette) || length(silhouette) != 1L || is.na(silhouette)) {
    stop("`silhouette` must be TRUE or FALSE.", call. = FALSE)
  }
  if (silhouette && n > 5000L) {
    stop("Silhouettes are limited to 5,000 included entities; use silhouette = FALSE to avoid pairwise distances.", call. = FALSE)
  }
  if (silhouette && !requireNamespace("cluster", quietly = TRUE)) {
    stop("Install 'cluster' for silhouettes or use silhouette = FALSE.", call. = FALSE)
  }
  fit <- with_preserved_rng_seed(seed, withCallingHandlers(
    kmeans(z, centers = as.integer(k), nstart = as.integer(nstart),
      iter.max = as.integer(iter.max), algorithm = "Hartigan-Wong"),
    warning = function(w) stop("K-means did not complete cleanly: ",
      conditionMessage(w), call. = FALSE)))
  if ((!is.null(fit$ifault) && fit$ifault != 0L) ||
      any(!is.finite(c(fit$centers, fit$withinss, fit$totss, fit$tot.withinss)))) {
    stop("K-means did not return finite, converged results; review scales and iterations.", call. = FALSE)
  }
  x <- prepared$x
  keep <- prepared$keep
  member <- data.frame(ID = x$row_summary$ID, Cluster = NA_integer_, Silhouette = NA_real_)
  member$Cluster[keep] <- as.integer(fit$cluster)
  widths <- rep(NA_real_, n)
  if (silhouette) {
    distance <- stats::dist(z)
    if (any(!is.finite(distance))) stop("Euclidean distances overflow; rescale the features.", call. = FALSE)
    widths <- unname(cluster::silhouette(fit$cluster, distance)[, "sil_width"])
  }
  member$Silhouette[keep] <- widths
  out <- list(membership = member,
    cluster_summary = data.frame(Cluster = seq_len(k), N = as.integer(fit$size),
      MeanSilhouette = vapply(seq_len(k), function(g) mean(widths[fit$cluster == g]), numeric(1))),
    profiles = cluster_feature_profiles(x$data[keep, x$features, drop = FALSE], fit$cluster),
    medoids = NULL, centers = fit$centers, withinss = fit$withinss,
    totss = fit$totss, tot.withinss = fit$tot.withinss, pca = pca,
    transformation = prepared$transformation, feature_data = x,
    settings = list(k = as.integer(k), method = "k-means", distance = "Euclidean",
      space = if (is.null(pca)) "Numeric features" else "Principal components",
      components = if (is.null(pca)) NULL else pca$settings$components,
      weights = prepared$weights, scale = prepared$scale, missing = missing,
      nstart = as.integer(nstart), iter.max = as.integer(iter.max), seed = as.integer(seed),
      iterations = fit$iter, convergence = fit$ifault %||% 0L,
      rng_kind = RNGkind(), silhouette = silhouette, included = n, excluded = sum(!keep)))
  class(out) <- "mfrm_clusters"
  out
}

prepare_numeric_features <- function(x, scale, weights, missing) {
  if (!inherits(x, "mfrm_features")) stop("`x` must be prepared with mfrm_features().", call. = FALSE)
  x <- mfrm_features(x$data, x$id, x$features, x$missing)
  if (!all(x$feature_summary$Type == "Numeric")) {
    stop("Select numeric external features with meaningful Euclidean units; categorical and logical features are not accepted.", call. = FALSE)
  }
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  keep <- x$row_summary$Complete
  if (any(!keep) && missing == "error") stop("Selected features contain missing values; review reasons or use missing = 'omit'.", call. = FALSE)
  if (sum(keep) < 2L) stop("At least two complete entities are required.", call. = FALSE)
  values <- as.matrix(x$data[keep, x$features, drop = FALSE])
  if (any(vapply(as.data.frame(values), function(v) length(unique(v)) < 2L, logical(1)))) {
    stop("Every selected feature must vary among included entities; remove constant features.", call. = FALSE)
  }
  if (is.null(weights)) weights <- stats::setNames(rep(1, ncol(values)), x$features)
  if (!is.numeric(weights) || is.complex(weights) || !is.null(dim(weights)) ||
      length(weights) != ncol(values) || is.null(names(weights)) || anyDuplicated(names(weights)) ||
      !setequal(names(weights), x$features) || any(!is.finite(weights) | weights <= 0)) {
    stop("`weights` must be a named, positive finite numeric vector covering each selected feature once.", call. = FALSE)
  }
  weights <- weights[x$features]
  relative <- weights / max(weights)
  if (any(relative == 0)) stop("Relative feature weights underflow; revise the weight ratios.", call. = FALSE)
  center <- colMeans(values)
  z <- sweep(values, 2L, center, "-")
  divisor <- if (scale) sqrt(colSums(z^2) / (nrow(z) - 1L)) else rep(1, ncol(z))
  z <- sweep(sweep(z, 2L, divisor, "/"), 2L, sqrt(relative), "*")
  if (any(!is.finite(z)) || any(!is.finite(divisor) | divisor <= 0) ||
      any(!is.finite(colSums(z^2)) | colSums(z^2) <= 0)) {
    stop("Numeric transformation or variance is not representable; rescale the features.", call. = FALSE)
  }
  list(x = x, keep = keep, coordinates = z, weights = weights, scale = scale,
    transformation = data.frame(Feature = x$features, Center = unname(center),
      Divisor = unname(divisor), WeightFactor = unname(sqrt(relative))))
}

check_feature_integer <- function(x, name, lower, upper = .Machine$integer.max) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop("`", name, "` must be an integer from ", lower, " to ", upper, ".", call. = FALSE)
  }
}
