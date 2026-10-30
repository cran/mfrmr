#' Compare exploratory groups across settings and clustering methods
#'
#' Compare existing external-feature clustering results without refitting.
#' Review how group membership changes when the selected features, number of
#' groups, feature weights, or clustering method changes, including paired comparisons across
#' the same imputations.
#'
#' @param analyses A named list of at least two results from [mfrm_cluster_pam()]
#'   and/or [mfrm_cluster_hierarchical()] or [mfrm_cluster_kmeans()],
#'   or a named list of results from [mfrm_cluster_imputed()]. Do not mix the two
#'   result types. Names must be unique and nonblank. All results must use the
#'   same entity IDs and included entities. Selected features may differ;
#'   features with the same name must retain their original values and types.
#'   Row and feature order may differ. IDs must refer to the same entities,
#'   not, for example, persons and raters that happen to share ID labels.
#'   For imputed results, shared completed features must also match by
#'   imputation number. When feature selections differ, all results must retain
#'   the same `mids` object, including its data and settings; each retained
#'   completion is checked against it using the optional `mice` package.
#'   Reuse one fitted imputation model. Results saved in separate RDS files
#'   can be compared after reopening when their retained imputation data,
#'   settings and formula expressions are unchanged. Formula environments
#'   are not evaluated in this comparison; no completions are regenerated.
#' @return An `mfrm_cluster_comparison` object containing:
#'   \itemize{
#'   \item `analysis_summary`: method, distance, fitted space, scaling, retained
#'     component count, linkage (unavailable for PAM and k-means), group
#'     count, selected feature count (`Features`), included/excluded entity counts,
#'     smallest/largest group sizes, and mean silhouette for each analysis and
#'     imputation. `Imputation` is `NA` for ordinary clustering results.
#'   \item `weights`: selected features and their supplied weights for each
#'     analysis. Unselected features have no row; they are not zero-weight inputs.
#'   \item `comparisons`: all pairs of analyses, with one row per imputation.
#'     `Pairs` counts unordered pairs of included entities, excluding self-pairs.
#'     `SplitPairs` were together in `First` and apart in `Second`; `JoinedPairs`
#'     were apart in `First` and together in `Second`. `ChangedFraction` is
#'     their sum divided by `Pairs`. `AdjustedRand` is the adjusted Rand index.
#'   \item `comparison_summary`: number of compared partitions, mean/minimum/
#'     maximum changed fraction, and mean adjusted Rand index for each analysis
#'     pair. `summary()` returns this table.
#'   \item `analyses`: the original results, including omitted IDs and all
#'     imputation-specific partitions, profiles, and settings.
#'   }
#' @details Comparisons align entities by ID and do not compare numeric group
#'   labels directly. Zero changed fraction and an adjusted Rand index of one
#'   indicate identical partitions, allowing arbitrary renumbering of groups.
#'   The changed fraction can be small when many pairs are separated in both
#'   partitions. The adjusted Rand index corrects agreement against random
#'   partitions with fixed group sizes; it can be negative. Neither measure
#'   identifies the correct group count or establishes group validity.
#'
#'   Mean silhouette and group sizes help describe each partition. Silhouette
#'   values calculated with different features or weights use different distances;
#'   their maximum is not an automatic criterion for selecting features or weights. Review
#'   group profiles in the retained `analyses` alongside the comparison.
#'
#'   With multiple imputations, feature selection is applied to the same
#'   completed data on both sides, without refitting the imputation model.
#'   A feature omitted from clustering can remain an imputation predictor;
#'   this comparison does not evaluate its removal from the imputation model.
#'   Every imputation is retained. Summary means and ranges
#'   describe sensitivity to settings across those imputations; they are not
#'   Rubin-pooled estimates, confidence intervals, or sampling stability.
#'   Conflicting shared feature values or different inclusion masks cause an error rather than a
#'   silent intersection of entities or imputations. No preferred setting,
#'   consensus partition, or hypothesis test is returned.
#'
#'   For an executable rater-attribute example that pairs multiple imputations
#'   across settings, see `vignette("mfrmr-external-features", package = "mfrmr")`.
#' @references Hubert, L. and Arabie, P. (1985). Comparing partitions.
#'   Journal of Classification, 2, 193--218. \doi{10.1007/BF01908075}.
#' @seealso [mfrm_features()], [mfrm_cluster_pam()], [mfrm_cluster_imputed()],
#'   [mfrm_cluster_hierarchical()]
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # Fictional attributes; experience is measured in completed years.
#'   raters <- data.frame(Rater = paste0("R", 1:8),
#'     ExperienceYears = c(1, 2, 3, 4, 11, 12, 13, 14),
#'     Specialty = rep(c("Language", "Science"), 4))
#'   features <- mfrm_features(raters, "Rater", c("ExperienceYears", "Specialty"))
#'   fits <- list(
#'     TwoGroups = mfrm_cluster_pam(features, k = 2),
#'     ThreeGroups = mfrm_cluster_pam(features, k = 3),
#'     ExperienceWeighted = mfrm_cluster_pam(features, k = 2,
#'       weights = c(ExperienceYears = 3, Specialty = 1)))
#'   comparison <- mfrm_cluster_compare(fits)
#'   summary(comparison)
#'   comparison$analysis_summary
#'   comparison$comparisons
#'   # Hold the entities and group count fixed while changing the feature set.
#'   experience <- mfrm_features(raters, "Rater", "ExperienceYears")
#'   feature_comparison <- mfrm_cluster_compare(list(
#'     BothFeatures = fits$TwoGroups,
#'     ExperienceOnly = mfrm_cluster_pam(experience, k = 2)))
#'   feature_comparison$analysis_summary
#'   feature_comparison$weights
#' }
#' @export
mfrm_cluster_compare <- function(analyses) {
  if (!is.list(analyses) || length(analyses) < 2L || is.null(names(analyses)) ||
      anyNA(names(analyses)) || any(!nzchar(trimws(names(analyses)))) ||
      anyDuplicated(names(analyses))) {
    stop("`analyses` must be a list of at least two results with unique, nonblank names.", call. = FALSE)
  }
  imputed <- all(vapply(analyses, inherits, logical(1), "mfrm_imputed_clusters"))
  if (!imputed && !all(vapply(analyses, inherits, logical(1), "mfrm_clusters"))) {
    stop("Supply only mfrm_cluster_pam() / mfrm_cluster_hierarchical() / mfrm_cluster_kmeans() results, or only mfrm_cluster_imputed() results.", call. = FALSE)
  }
  partitions <- if (imputed) lapply(analyses, `[[`, "analyses") else lapply(analyses, list)
  m <- length(partitions[[1L]])
  if (m < 1L || any(lengths(partitions) != m) ||
      (imputed && any(vapply(analyses, function(a) !identical(length(a$analyses),
        a$settings$imputations), logical(1))))) {
    stop("All results must retain the same number of imputations.", call. = FALSE)
  }
  original <- analyses[[1L]]$feature_data
  ids <- original$row_summary$ID
  feature_sets <- lapply(analyses, function(a) a$feature_data$features)
  different_features <- any(!vapply(feature_sets, setequal, logical(1), feature_sets[[1L]]))
  if (imputed && different_features) {
    model_record <- function(a) {
      model <- a$imputation_model
      # Comparison uses retained completions, never formula evaluation environments.
      for (i in seq_along(model$formulas)) {
        if (inherits(model$formulas[[i]], "formula"))
          environment(model$formulas[[i]]) <- emptyenv()
      }
      model
    }
    reference_model <- model_record(analyses[[1L]])
    if (any(!vapply(analyses, function(a)
        identical(model_record(a), reference_model), logical(1)))) {
      stop("Different feature selections must reuse the same fitted mids object.", call. = FALSE)
    }
  }
  if (imputed && different_features && !requireNamespace("mice", quietly = TRUE)) {
    stop("Install the optional 'mice' package to compare imputed feature selections.", call. = FALSE)
  }
  # Reuse feature validation, aligning both original and completed tables by ID.
  feature_table <- function(x) {
    x <- mfrm_features(x$data, x$id, x$features, x$missing)
    if (!setequal(x$row_summary$ID, ids)) {
      stop("All results must use the same entity IDs.", call. = FALSE)
    }
    data <- x$data[match(ids, x$row_summary$ID), x$features, drop = FALSE]
    rownames(data) <- NULL
    data
  }
  # Accumulate every feature, so conflicts between later analyses are also checked.
  check_shared <- function(reference, data, message) {
    shared <- intersect(names(data), names(reference))
    if (length(shared) &&
        !isTRUE(all.equal(as.list(data[shared]), reference[shared], tolerance = 0))) {
      stop(message, call. = FALSE)
    }
    reference[names(data)] <- as.list(data)
    reference
  }
  original_data <- list()
  completed_data <- rep(list(list()), m)
  included <- !is.na(partitions[[1L]][[1L]]$membership$Cluster[
    match(ids, partitions[[1L]][[1L]]$membership$ID)])
  memberships <- summaries <- weight_tables <- vector("list", length(analyses))
  for (i in seq_along(analyses)) {
    reviewed <- analyses[[i]]$feature_data
    if (imputed && different_features) {
      original_data <- check_shared(original_data,
        feature_table(mfrm_features(analyses[[i]]$imputation_model$data,
          reviewed$id, reviewed$features)),
        "Shared features must use the same original feature values and types.")
    }
    original_data <- check_shared(original_data, feature_table(reviewed),
      "Shared features must use the same original feature values and types.")
    memberships[[i]] <- vector("list", m)
    summaries[[i]] <- vector("list", m)
    for (j in seq_len(m)) {
      a <- partitions[[i]][[j]]
      if (!inherits(a, "mfrm_clusters")) {
        stop("Every retained partition must be an external-feature clustering result.", call. = FALSE)
      }
      data <- feature_table(a$feature_data)
      if (!setequal(names(data), feature_sets[[i]])) {
        stop("Each retained partition must use its analysis's selected features.", call. = FALSE)
      }
      message <- paste0("Completed feature values and types must match at imputation ", j, ".")
      completed_data[[j]] <- check_shared(completed_data[[j]], data, message)
      if (imputed && different_features) {
        # Disjoint selections have no shared values to establish completion pairing.
        expected <- mice::complete(analyses[[i]]$imputation_model, action = j)
        for (name in reviewed$features) {
          if (is.logical(reviewed$data[[name]]) && is.numeric(expected[[name]]) &&
              all(is.na(expected[[name]]) | expected[[name]] %in% c(0, 1))) {
            expected[[name]] <- as.logical(expected[[name]])
          }
        }
        expected <- feature_table(mfrm_features(expected, reviewed$id, names(data)))
        if (!isTRUE(all.equal(data, expected, tolerance = 0))) stop(message, call. = FALSE)
      }
      member <- a$membership
      if (anyDuplicated(member$ID) || !setequal(member$ID, ids)) {
        stop("Membership IDs must match the original feature table.", call. = FALSE)
      }
      labels <- member$Cluster[match(ids, member$ID)]
      if (!identical(!is.na(labels), included)) {
        stop("All partitions must include the same entities; omitted IDs cannot be silently dropped.", call. = FALSE)
      }
      labels <- labels[included]
      if (!is.numeric(labels) || any(!is.finite(labels)) ||
          length(unique(labels)) != a$settings$k || a$settings$k < 2L ||
          a$settings$k >= length(labels)) {
        stop("Each partition must retain valid group memberships and its group count.", call. = FALSE)
      }
      memberships[[i]][[j]] <- labels
      sizes <- table(labels)
      summaries[[i]][[j]] <- data.frame(Analysis = names(analyses)[i],
        Method = a$settings$method, Linkage = a$settings$linkage %||% NA_character_,
        Imputation = if (imputed) j else NA_integer_, K = a$settings$k,
        Features = length(feature_sets[[i]]),
        Included = sum(included), Excluded = sum(!included),
        MinGroupSize = min(sizes), MaxGroupSize = max(sizes),
        MeanSilhouette = mean_available_silhouette(member$Silhouette),
        Distance = a$settings$distance, Space = a$settings$space %||% "Mixed features",
        Scaling = if (identical(a$settings$distance, "Gower")) "Range" else
          if (isTRUE(a$settings$scale)) "Sample SD" else "Original units",
        Components = a$settings$components %||% NA_integer_)
    }
    weights <- partitions[[i]][[1L]]$settings$weights
    weight_tables[[i]] <- data.frame(Analysis = names(analyses)[i],
      Feature = names(weights), Weight = unname(weights))
  }
  pairs <- utils::combn(seq_along(analyses), 2L)
  comparisons <- comparison_summary <- vector("list", ncol(pairs))
  for (p in seq_len(ncol(pairs))) {
    first <- pairs[1L, p]
    second <- pairs[2L, p]
    rows <- vector("list", m)
    for (j in seq_len(m)) {
      cross <- table(memberships[[first]][[j]], memberships[[second]][[j]])
      total <- choose(sum(cross), 2)
      together <- sum(choose(cross, 2))
      together_first <- sum(choose(rowSums(cross), 2))
      together_second <- sum(choose(colSums(cross), 2))
      expected <- together_first * together_second / total
      rows[[j]] <- data.frame(First = names(analyses)[first], Second = names(analyses)[second],
        Imputation = if (imputed) j else NA_integer_, Included = sum(included), Pairs = total,
        SplitPairs = together_first - together, JoinedPairs = together_second - together,
        ChangedFraction = (together_first + together_second - 2 * together) / total,
        AdjustedRand = (together - expected) / ((together_first + together_second) / 2 - expected))
    }
    comparisons[[p]] <- do.call(rbind, rows)
    values <- comparisons[[p]]
    comparison_summary[[p]] <- data.frame(First = names(analyses)[first], Second = names(analyses)[second],
      Partitions = m, Included = sum(included), Pairs = values$Pairs[1L],
      MeanChangedFraction = mean(values$ChangedFraction),
      MinChangedFraction = min(values$ChangedFraction), MaxChangedFraction = max(values$ChangedFraction),
      MeanAdjustedRand = mean(values$AdjustedRand))
  }
  out <- list(analysis_summary = do.call(rbind, lapply(summaries, function(x) do.call(rbind, x))),
    weights = do.call(rbind, weight_tables), comparisons = do.call(rbind, comparisons),
    comparison_summary = do.call(rbind, comparison_summary), analyses = analyses)
  class(out) <- "mfrm_cluster_comparison"
  out
}

#' @rdname mfrm_cluster_compare
#' @param x,object An object returned by [mfrm_cluster_compare()].
#' @param ... Reserved for method compatibility.
#' @export
print.mfrm_cluster_comparison <- function(x, ...) {
  cat("Exploratory grouping sensitivity to settings and clustering methods\n")
  print(x$comparison_summary, row.names = FALSE)
  cat("ChangedFraction counts pairs whose together/apart status changes; group numbers are arbitrary.\n")
  cat("Summaries are descriptive, not pooled inference, sampling stability, or automatic setting selection.\n")
  invisible(x)
}

#' @rdname mfrm_cluster_compare
#' @export
summary.mfrm_cluster_comparison <- function(object, ...) object$comparison_summary
