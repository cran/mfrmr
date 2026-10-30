#' Prepare external features for exploratory grouping
#'
#' Keep one row per Person, rater, or other entity, review feature availability,
#' and optionally attach user-supplied reasons for missing values.
#'
#' @param data A data frame with one row per entity. Repeated rating rows must
#'   be summarized or joined to an entity-level table before calling this function.
#' @param id Name of the unique, nonmissing identifier column. Character, factor,
#'   or finite numeric identifiers are stored as character strings without trimming.
#' @param features Explicit character vector of feature column names, excluding
#'   the identifier. Numeric, character, factor, ordered factor, and logical
#'   columns are supported. Character columns become nominal factors; logical
#'   columns represent symmetric binary features. Ordered factor levels retain
#'   their declared order. Dates, lists, and matrices must be converted explicitly.
#' @param missing_reasons Optional data frame with columns `ID`, `Feature`, and
#'   `Reason`, one row per annotated missing cell. IDs refer to `id`; features
#'   must be selected columns. Reasons for observed or unknown cells are refused.
#'   Unannotated missing cells are labelled "Not supplied". Reasons are supplied
#'   by the user, not inferred missing-data mechanisms.
#' @return An `mfrm_features` object containing `data`, `id`, `features`, a
#'   `feature_summary`, a `row_summary`, and `missing` (all missing cells and
#'   their reasons). No values are imputed or rows discarded.
#' @details Numeric `NA` and `NaN` are missing. Infinite numeric values and blank
#'   categorical labels are refused; replace missing markers with `NA` explicitly.
#'   Constant and entirely missing features remain available for review but
#'   cannot be used by [mfrm_cluster_pam()]. IDs and unselected columns do not enter
#'   distances. These functions are intended for external attributes such as
#'   training, experience, or specialization. They do not propagate uncertainty
#'   from estimated ability, severity, or fit statistics.
#'   See `vignette("mfrmr-external-features")` for a complete rater-attribute
#'   example including missingness review and multiple imputation.
#' @seealso [mfrm_cluster_pam()], [mfrm_cluster_imputed()], [mfrm_cluster_compare()],
#'   [mfrm_pca()], [mfrm_cluster_kmeans()]
#' @examples
#' # Fictional rater attributes; experience is measured in completed years.
#' raters <- data.frame(
#'   Rater = paste0("R", 1:6),
#'   ExperienceYears = c(1, NA, 3, 12, 13, 14),
#'   Specialty = c("Language", "Language", "Language", "Science", "Science", "Science")
#' )
#' reasons <- data.frame(ID = "R2", Feature = "ExperienceYears", Reason = "Not recorded")
#' features <- mfrm_features(raters, "Rater", c("ExperienceYears", "Specialty"),
#'                           missing_reasons = reasons)
#' summary(features)
#' features$missing
#' features$row_summary
#' @export
mfrm_features <- function(data, id, features, missing_reasons = NULL) {
  if (!is.data.frame(data) || nrow(data) == 0L || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data)))) {
    stop("`data` must be a nonempty data frame with unique, nonmissing column names.", call. = FALSE)
  }
  if (!is.character(id) || length(id) != 1L || is.na(id) || !id %in% names(data)) {
    stop("`id` must name one column in `data`.", call. = FALSE)
  }
  if (!is.character(features) || !length(features) || anyNA(features) ||
      anyDuplicated(features) || any(!features %in% names(data)) || id %in% features) {
    stop("`features` must name distinct columns in `data`, excluding `id`.", call. = FALSE)
  }
  ids <- data[[id]]
  if (!(is.character(ids) || is.factor(ids) || (is.numeric(ids) && !is.object(ids) && !is.complex(ids))) ||
      !is.null(dim(ids)) || anyNA(ids) ||
      (is.numeric(ids) && any(!is.finite(ids)))) {
    stop("Identifiers must be nonmissing character, factor, or finite numeric values.", call. = FALSE)
  }
  ids <- as.character(ids)
  if (any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("Identifiers must be nonblank and unique: supply one row per entity.", call. = FALSE)
  }
  selected <- as.data.frame(data[, c(id, features), drop = FALSE])
  selected[[id]] <- ids
  types <- character(length(features))
  for (j in seq_along(features)) {
    name <- features[j]
    value <- selected[[name]]
    if (!is.null(dim(value)) ||
        !(is.factor(value) || (!is.object(value) &&
          ((is.numeric(value) && !is.complex(value)) || is.character(value) || is.logical(value))))) {
      stop("Feature '", name, "' has an unsupported type; convert it explicitly.", call. = FALSE)
    }
    if (is.numeric(value) && any(is.infinite(value))) {
      stop("Feature '", name, "' contains infinite values.", call. = FALSE)
    }
    if ((is.character(value) || is.factor(value)) &&
        any(!is.na(value) & !nzchar(trimws(as.character(value))))) {
      stop("Feature '", name, "' contains blank labels; use NA for missing values.", call. = FALSE)
    }
    types[j] <- if (is.ordered(value)) "Ordinal" else if (is.factor(value) ||
      is.character(value)) "Nominal" else if (is.logical(value)) "Symmetric binary" else "Numeric"
    if (is.character(value)) selected[[name]] <- factor(value)
  }
  values <- selected[, features, drop = FALSE]
  absent <- is.na(values)
  cells <- which(absent, arr.ind = TRUE)
  missing <- data.frame(ID = ids[cells[, 1L]], Feature = features[cells[, 2L]],
                        stringsAsFactors = FALSE)
  if (!is.null(missing_reasons)) {
    if (!is.data.frame(missing_reasons) || anyDuplicated(names(missing_reasons)) ||
        !all(c("ID", "Feature", "Reason") %in% names(missing_reasons))) {
      stop("`missing_reasons` must be a data frame with ID, Feature, and Reason columns.", call. = FALSE)
    }
    reasons <- missing_reasons[, c("ID", "Feature", "Reason"), drop = FALSE]
    if (is.numeric(reasons$ID) && !is.object(reasons$ID) && !is.complex(reasons$ID) &&
        all(is.finite(reasons$ID))) reasons$ID <- as.character(reasons$ID)
    if (!all(vapply(reasons, function(v) (is.character(v) || is.factor(v)) &&
      is.null(dim(v)), logical(1)))) {
      stop("Missing-reason ID, Feature, and Reason columns must be character or factor.", call. = FALSE)
    }
    reasons[] <- lapply(reasons, as.character)
    if (anyNA(reasons) || any(!nzchar(trimws(reasons$Reason))) ||
        anyDuplicated(reasons[c("ID", "Feature")])) {
      stop("Missing reasons must be nonblank, nonmissing, and unique per ID and feature.", call. = FALSE)
    }
    row <- match(reasons$ID, ids)
    column <- match(reasons$Feature, features)
    if (anyNA(row) || anyNA(column) || !all(absent[cbind(row, column)])) {
      stop("Missing reasons must refer only to missing cells in the selected features.", call. = FALSE)
    }
    missing <- dplyr::left_join(missing, reasons, by = c("ID", "Feature"))
  } else {
    missing$Reason <- rep(NA_character_, nrow(missing))
  }
  missing$Reason[is.na(missing$Reason)] <- "Not supplied"
  out <- list(
    data = selected, id = id, features = features,
    feature_summary = data.frame(
      Feature = features, Type = types, Observed = colSums(!absent),
      Missing = colSums(absent),
      Distinct = vapply(values, function(v) length(unique(v[!is.na(v)])), integer(1)),
      row.names = NULL
    ),
    row_summary = data.frame(ID = ids, MissingFeatures = rowSums(absent),
                             Complete = rowSums(absent) == 0L),
    missing = missing
  )
  class(out) <- "mfrm_features"
  out
}

#' @rdname mfrm_features
#' @param x,object An object returned by the corresponding function.
#' @param ... Reserved for method compatibility.
#' @export
print.mfrm_features <- function(x, ...) {
  cat("External feature review\n")
  cat(nrow(x$data), "entities;", length(x$features), "features;",
      sum(x$row_summary$Complete), "complete rows\n")
  print(x$feature_summary, row.names = FALSE)
  cat("Missing values are retained. Reasons are user-supplied; no imputation is performed.\n")
  invisible(x)
}

#' @rdname mfrm_features
#' @export
summary.mfrm_features <- function(object, ...) object$feature_summary

#' Group mixed external features using Gower distances and PAM
#'
#' Partition an entity-level feature table using Gower dissimilarities and
#' partitioning around medoids (PAM) from the optional `cluster` package.
#' Use `mfrm_cluster_pam()` for this method. `mfrm_cluster()` is retained as
#' an identical compatibility alias; it does not select an algorithm for you.
#'
#' @param x An object returned by [mfrm_features()].
#' @param k Number of groups, an integer from 2 to one less than the number of
#'   included entities, and no greater than their number of distinct profiles.
#'   The number is chosen by the user, not optimized automatically.
#' @param weights Optional named, strictly positive finite numeric vector with
#'   one weight per selected feature. Defaults to equal weights. Names, not
#'   vector order, identify features. Remove a feature to exclude it.
#' @param missing Either `"error"` (default) or `"omit"`. The latter excludes
#'   incomplete entities explicitly, retaining their IDs, reasons, and missing
#'   group membership in the result. No values are imputed.
#' @return An `mfrm_clusters` object containing `membership` (ID, Cluster,
#'   Medoid, Silhouette), `cluster_summary`, numeric and categorical `profiles`,
#'   `medoids`, the reviewed `feature_data`, and `settings`. Silhouette values
#'   describe within-sample separation, not membership probabilities or
#'   resampling stability. `summary()` returns the cluster-size/silhouette table.
#'   Numeric profiles report the mean, median and sample SD in the original
#'   feature units. SD is `NA` for a one-entity group and zero for identical
#'   values in a group with at least two entities. A singleton's silhouette is
#'   zero by convention; it is not evidence of group separation.
#' @details Numeric differences are divided by the feature range among included
#'   entities; nominal features use match/mismatch and logical features use
#'   symmetric binary differences. Ordered factors use equally spaced integer
#'   codes in their declared level order, divided by the observed code range.
#'   Unused intermediate levels retain their spacing. The Gower dissimilarity
#'   is the weighted mean of these feature contributions, dividing by the sum
#'   of the selected weights. Weights apply to features, not entities.
#'   Numeric 0/1 features are treated as numeric, not asymmetric presence/absence.
#'   Every selected feature must vary among included entities. Gower scaling,
#'   feature types, weights, omission policy, and numeric ranges are retained.
#'   Numeric ranges and relative weights must be representable without overflow
#'   or underflow; rescale features or revise extreme weight ratios if refused.
#'
#'   PAM uses deterministic BUILD/SWAP initialization; tied distances may admit
#'   alternative partitions. Group numbers are arbitrary labels. These are
#'   exploratory groups, not latent classes, ability estimates, or assessments
#'   of rater quality. No inference, uncertainty propagation, new-entity
#'   classification, or resampling stability is provided. Omission may change
#'   both the sample and numeric ranges and does not correct missing-data bias.
#'
#'   Pairwise distances require quadratic memory. This interface is limited
#'   to 5,000 included entities. It does not silently sample larger inputs.
#'   This input limit does not guarantee low memory use or acceptable run time.
#'
#'   For feature selection, missingness review, and multiple-imputation examples,
#'   see `vignette("mfrmr-external-features", package = "mfrmr")`.
#' @section Reading clustering metrics:
#' For entity i, a(i) is its mean dissimilarity to other members of its own
#' group; b(i) is the smallest mean dissimilarity to any other group.
#' Silhouette = (b - a) / max(a, b), between -1 and 1. Values near 1 indicate
#' compact, separated placement in this fitted geometry; near zero indicates
#' a boundary; negative values indicate closer average proximity to another
#' group. A singleton receives 0 by convention; zero does not establish overlap.
#' A zero denominator also receives 0. Overall means weight entities equally;
#' cluster means use that cluster's included entities. Excluded IDs have NA.
#' A silhouette is neither a membership probability nor evidence of stability.
#'
#' Gower distance is unitless in `[0, 1]` here: numeric absolute differences
#' divided by included-sample ranges, nominal/binary mismatches, and ordinal
#' rank differences are combined using the specified feature weights.
#' The PAM medoid is an actual included entity minimizing within-group total
#' dissimilarity; it is not a feature-wise mean. Profile Mean/Median/SD retain
#' original units; categorical Proportion uses included group size as denominator.
#' A large group or a high silhouette is not a better-performing group.
#' Compare partitions under a fixed feature/distance definition; changing the
#' sample, weights or scaling changes what these numbers measure.
#' See [cluster::silhouette()] and Rousseeuw (1987),
#' \doi{10.1016/0377-0427(87)90125-7}, for the silhouette definition.
#' @seealso [mfrm_features()], [mfrm_cluster_imputed()], [mfrm_cluster_compare()],
#'   [mfrm_cluster_hierarchical()], [plot.mfrm_clusters()], [cluster::daisy()], [cluster::pam()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("mfrm_cluster", run.dontrun = TRUE)
#' \dontrun{
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # Fictional raters; R2 has unrecorded experience, not zero years.
#'   raters <- data.frame(Rater = paste0("R", 1:6),
#'     ExperienceYears = c(1, NA, 3, 12, 13, 14),
#'     Specialty = rep(c("Language", "Science"), each = 3))
#'   features <- mfrm_features(raters, "Rater", c("ExperienceYears", "Specialty"))
#'   # The default stops on missing features. Here omission is explicit.
#'   groups <- mfrm_cluster_pam(features, k = 2,
#'     weights = c(ExperienceYears = 2, Specialty = 1), missing = "omit")
#'   groups$membership  # R2 remains present with unavailable membership.
#'   summary(groups)
#'   groups$profiles
#'   groups$medoids
#' }
#' }
#' @name mfrm_cluster
#' @rdname mfrm_cluster
#' @export
mfrm_cluster_pam <- function(x, k, weights = NULL, missing = c("error", "omit")) {
  cluster_external_features(x, k, weights, match.arg(missing))
}

#' @rdname mfrm_cluster
#' @export
mfrm_cluster <- mfrm_cluster_pam

cluster_external_features <- function(x, k, weights, missing, linkage = NULL) {
  if (!inherits(x, "mfrm_features")) {
    stop("`x` must be prepared with mfrm_features().", call. = FALSE)
  }
  x <- mfrm_features(x$data, x$id, x$features, x$missing)
  missing <- match.arg(missing, c("error", "omit"))
  keep <- x$row_summary$Complete
  if (any(!keep) && missing == "error") {
    stop("Selected features contain missing values. Inspect `x$missing`; use missing = 'omit' to exclude incomplete entities explicitly.", call. = FALSE)
  }
  values <- x$data[keep, x$features, drop = FALSE]
  n <- nrow(values)
  if (!is.numeric(k) || is.complex(k) || length(k) != 1L || !is.finite(k) || k != floor(k) || k < 2 || k >= n) {
    stop("`k` must be an integer from 2 to one less than the number of included entities.", call. = FALSE)
  }
  # ponytail: quadratic distances; use a scalable mixed-data backend if larger cohorts are needed.
  if (n > 5000L) stop("Clustering is limited to 5,000 included entities because pairwise distances require quadratic memory.", call. = FALSE)
  if (any(vapply(values, function(v) length(unique(v)) < 2L, logical(1)))) {
    stop("Every selected feature must vary among included entities; remove constant features and retry.", call. = FALSE)
  }
  if (k > nrow(unique(values))) stop("`k` exceeds the number of distinct included feature profiles.", call. = FALSE)
  if (is.null(weights)) weights <- stats::setNames(rep(1, length(x$features)), x$features)
  if (!is.numeric(weights) || is.complex(weights) || !is.null(dim(weights)) || length(weights) != length(x$features) ||
      is.null(names(weights)) || anyDuplicated(names(weights)) ||
      !setequal(names(weights), x$features) || any(!is.finite(weights) | weights <= 0)) {
    stop("`weights` must be a named, positive finite numeric vector covering each selected feature once.", call. = FALSE)
  }
  weights <- weights[x$features]
  distance_weights <- weights / max(weights)
  if (any(distance_weights == 0)) {
    stop("Relative feature weights underflow to zero; revise the weight ratios.", call. = FALSE)
  }
  numeric_names <- x$features[vapply(values, is.numeric, logical(1))]
  ranges <- data.frame(Feature = numeric_names,
    Minimum = vapply(values[numeric_names], min, numeric(1)),
    Maximum = vapply(values[numeric_names], max, numeric(1)), row.names = NULL)
  if (any(!is.finite(ranges$Maximum - ranges$Minimum))) {
    stop("Numeric feature ranges overflow; rescale these features before clustering.", call. = FALSE)
  }
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Install the optional 'cluster' package for external-feature clustering.", call. = FALSE)
  }
  rownames(values) <- as.character(seq_len(n))
  binary <- which(vapply(values, is.logical, logical(1)))
  distance <- cluster::daisy(values, metric = "gower", weights = unname(distance_weights),
                             type = list(symm = binary), warnBin = FALSE)
  if (any(!is.finite(distance))) stop("The selected features do not produce finite distances; review their scales.", call. = FALSE)
  if (is.null(linkage)) {
    fit <- cluster::pam(distance, k = as.integer(k), diss = TRUE, variant = "original",
                      keep.diss = FALSE, keep.data = FALSE)
    labels <- fit$clustering
    silhouette <- numeric(n)
    silhouette[as.integer(rownames(fit$silinfo$widths))] <- fit$silinfo$widths[, "sil_width"]
  } else {
    tree <- hclust(distance, method = linkage)
    tree$labels <- x$row_summary$ID[keep]
    labels <- cutree(tree, k = as.integer(k))
    silhouette <- unname(cluster::silhouette(labels, distance)[, "sil_width"])
  }
  membership <- data.frame(ID = x$row_summary$ID, Cluster = NA_integer_,
                           Medoid = NA, Silhouette = NA_real_)
  membership$Cluster[keep] <- as.integer(labels)
  if (is.null(linkage)) {
    membership$Medoid[keep] <- seq_len(n) %in% fit$id.med
  } else {
    membership$Medoid <- NULL
  }
  membership$Silhouette[keep] <- silhouette
  sizes <- as.integer(tabulate(labels, nbins = k))
  cluster_summary <- data.frame(Cluster = seq_len(k), N = sizes,
    MeanSilhouette = as.numeric(tapply(silhouette, labels, mean)))
  profiles <- cluster_feature_profiles(values, labels)
  out <- list(membership = membership, cluster_summary = cluster_summary,
    profiles = profiles,
    medoids = if (is.null(linkage)) x$data[which(keep)[fit$id.med], c(x$id, x$features), drop = FALSE] else NULL,
    feature_data = x,
    settings = list(k = as.integer(k), method = "PAM", distance = "Gower",
      weights = weights, numeric_ranges = ranges, missing = missing,
      included = n, excluded = sum(!keep)))
  class(out) <- "mfrm_clusters"
  if (!is.null(linkage)) {
    out$tree <- tree
    out$settings$method <- "Hierarchical"
    out$settings$linkage <- linkage
    class(out) <- c("mfrm_hierarchical_clusters", "mfrm_clusters")
  }
  out
}

cluster_feature_profiles <- function(values, labels) {
  numeric_profiles <- categorical_profiles <- list()
  for (g in seq_len(max(labels))) for (name in names(values)) {
    value <- values[[name]][labels == g]
    if (is.numeric(value)) {
      spread <- stats::sd(value)
      if (length(value) > 1L && any(value != value[1L]) &&
          (!is.finite(spread) || spread < sqrt(.Machine$double.xmin))) {
        # SD can be representable even when its squared intermediate is not.
        shifted <- value - value[1L]
        unit <- max(abs(shifted))
        spread <- stats::sd(shifted / unit) * unit
      }
      numeric_profiles[[length(numeric_profiles) + 1L]] <- data.frame(
        Cluster = g, Feature = name, N = length(value), Mean = mean(value),
        Median = stats::median(value), SD = spread)
    } else {
      counts <- table(if (is.logical(value)) factor(value, levels = c(FALSE, TRUE)) else value)
      categorical_profiles[[length(categorical_profiles) + 1L]] <- data.frame(
        Cluster = g, Feature = name, Level = names(counts), N = as.integer(counts),
        Proportion = as.numeric(counts) / length(value))
    }
  }
  list(
      numeric = dplyr::bind_rows(data.frame(Cluster = integer(), Feature = character(),
        N = integer(), Mean = numeric(), Median = numeric(), SD = numeric()), numeric_profiles),
      categorical = dplyr::bind_rows(data.frame(Cluster = integer(), Feature = character(),
        Level = character(), N = integer(), Proportion = numeric()), categorical_profiles))
}

#' @rdname mfrm_cluster
#' @param object An object returned by [mfrm_cluster_pam()].
#' @param ... Reserved for method compatibility.
#' @export
print.mfrm_clusters <- function(x, ...) {
  cat(sprintf("Exploratory feature groups (%s distance, %s%s)\n", x$settings$distance, x$settings$method,
    if (!is.null(x$settings$linkage)) paste0(" / ", x$settings$linkage, " linkage") else ""))
  cat(x$settings$included, "entities included;", x$settings$excluded, "excluded for missing features\n")
  print(x$cluster_summary, row.names = FALSE)
  cat("Group labels are arbitrary. Silhouette describes separation in this sample; stability is not assessed.\n")
  cat("These groups do not establish ability levels or rater quality.\n")
  invisible(x)
}

#' @rdname mfrm_cluster
#' @export
summary.mfrm_clusters <- function(object, ...) object$cluster_summary

#' Compare exploratory groups across external-feature imputations
#'
#' Apply the same external-feature clustering analysis to each completed data set from `mice`,
#' retaining all analyses and the fraction of imputations in which each pair
#' of entities belongs to the same group.
#'
#' @param x An original feature table reviewed with [mfrm_features()].
#' @param imputed A `mids` object from the optional `mice` package, containing
#'   the identifier, selected features, and at least two imputations. Its original
#'   data must match `x` by ID. Auxiliary variables may be included in the
#'   imputation model without becoming clustering features.
#' @param impute A data frame with `ID` and `Feature` columns explicitly listing
#'   the missing cells to impute. For example, select appropriate rows from
#'   `x$missing` after reviewing their reasons. Extra columns are ignored.
#'   The selected-feature entries in `imputed$where` must match this selection;
#'   other missing cells must remain missing in every completed data set.
#'   If all selected features are already complete, supply the empty
#'   `x$missing` table. This retains one partition per imputation so the result
#'   can be compared with other feature selections from the same model.
#' @param k,weights As in [mfrm_cluster_pam()]. The same choices apply to every
#'   imputation.
#' @param missing Either `"error"` (default) or `"omit"`, applied to missing
#'   features remaining after imputation. Explicit omission retains excluded
#'   IDs with unavailable memberships and pairwise proportions.
#' @param method `"pam"` (default), `"hierarchical"`, or `"kmeans"`. The same method is used
#'   for every completion. Hierarchical analyses retain separate trees.
#' @param linkage For `method = "hierarchical"`, `"average"` (the default when
#'   `NULL`) or `"complete"`; see [mfrm_cluster_hierarchical()]. Must be `NULL`
#'   for PAM and k-means, which have no linkage.
#' @param components For k-means, `NULL` uses numeric features directly. An integer
#'   fits PCA separately in each completion and uses that many leading components.
#'   No whitening, pooled PCA basis or component-sign matching is performed.
#' @param scale,nstart,iter.max,seed,silhouette For `method = "kmeans"`, as in
#'   [mfrm_cluster_kmeans()]. These arguments must not be supplied for other methods.
#'   The same seed and settings are used for every completion; sensitivity still
#'   includes possible changes in local optima. Review alternative seeds separately.
#' @return An `mfrm_imputed_clusters` object containing `analyses` (one
#'   `mfrm_clusters` object per imputation), `co_membership` (a symmetric matrix
#'   indexed by ID), `analysis_summary`, original `feature_data`, `imputed_cells`
#'   with original reasons, the full `imputation_model`, and `settings`.
#'   Every available matrix entry uses all `settings$imputations` analyses as
#'   its denominator; pairs involving excluded entities are `NA`, including
#'   their diagonal entries. `summary()` returns per-imputation counts and mean
#'   silhouette widths. Numeric ranges are retained in each analysis.
#'   Hierarchical partitions inherit from `mfrm_clusters` and retain their
#'   trees. Method and linkage are retained in `settings`.
#' @details Fit and review the imputation model using [mice::mice()] before
#'   calling this function. Choose methods, predictors (including relevant
#'   auxiliary variables), iteration count, and number of imputations for the
#'   intended analysis. Exclude the identifier from imputation and prediction.
#'   Inspect model diagnostics, including `loggedEvents` and chain behavior;
#'   this adapter checks data consistency, not convergence or model adequacy.
#'   The original model and its diagnostics remain available in the result.
#'
#'   Missing reasons do not identify a statistical missing-data mechanism.
#'   Do not impute structurally undefined attributes such as "not applicable".
#'   Specify eligible cells through `where` when fitting `mice` and list the
#'   same cells in `impute`. An incomplete predictor that is not imputed can
#'   prevent imputation of other variables; configure predictors accordingly.
#'   Assumptions about nonresponse require substantive justification; this
#'   function does not correct bias automatically or impute rating responses.
#'   Nonresponse related to unobserved values needs separate sensitivity
#'   assumptions; observed-data checks cannot establish their adequacy.
#'
#'   IDs, observed feature values, feature types, and factor levels/order must
#'   be preserved. Numeric 0/1 completions of originally logical features are
#'   restored to logical values. Row order is restored by ID. Every selected
#'   cell must be completed in every imputation. Any invalid completion or failed clustering
#'   stops the comparison with its imputation number; no failures are discarded.
#'   Remaining missingness, and thus the included sample, is the same across
#'   imputations. Gower numeric ranges are recalculated in each completed sample;
#'   differences may reflect changes in both feature values and scaling.
#'   Numeric k-means likewise recalculates centers, standard deviations and, when
#'   requested, the PCA basis in each completion. Only partitions are compared;
#'   component scores/loadings and group labels are not averaged or Rubin-pooled.
#'   Hierarchical trees belong to individual completions; no pooled tree,
#'   consensus hierarchy, or branch-support estimate is returned.
#'
#'   Co-membership proportions are invariant to arbitrary group numbering.
#'   They describe sensitivity to the supplied imputations, conditional on the
#'   imputation model, features, weights, and group count. They are not posterior
#'   membership probabilities, sampling stability, or Rubin-pooled estimates.
#'   No consensus partition, confidence interval, or automatic group selection
#'   is produced. Both the full ID-indexed matrix and pairwise distances require
#'   quadratic memory, so this comparison is limited to 5,000 total entities.
#'   This input limit is not a memory or run-time guarantee. Retaining all
#'   completed analyses also increases memory use with the imputation count.
#'
#'   The short example below illustrates the interface using `mice::nhanes2`.
#'   For a complete example with fictional rater attributes, explicit missingness
#'   reasons, imputation diagnostics, and comparisons of group counts and weights,
#'   see `vignette("mfrmr-external-features", package = "mfrmr")`.
#' @seealso [mfrm_features()], [mfrm_cluster_pam()], [mfrm_cluster_compare()],
#'   [mfrm_cluster_hierarchical()], [plot.mfrm_clusters()],
#'   [mice::mice()], [mice::complete()]
#' @examples
#' # Optional integration example (loads additional modelling or rendering packages).
#' # example("mfrm_cluster_imputed", run.dontrun = TRUE)
#' \dontrun{
#' if (requireNamespace("mice", quietly = TRUE) &&
#'     requireNamespace("cluster", quietly = TRUE)) {
#'   attributes <- mice::nhanes2
#'   attributes$Person <- paste0("P", seq_len(nrow(attributes)))
#'   review <- mfrm_features(attributes, "Person", c("bmi", "chl"))
#'   # Here all missing selected attributes are assumed eligible after review.
#'   cells <- review$missing
#'   method <- mice::make.method(attributes)
#'   method["Person"] <- ""
#'   predictors <- mice::make.predictorMatrix(attributes)
#'   predictors[, "Person"] <- 0
#'   predictors["Person", ] <- 0
#'   # Small settings illustrate the API, not an adequacy recommendation.
#'   model <- mice::mice(attributes, m = 3, maxit = 2, method = method,
#'     predictorMatrix = predictors, seed = 42, printFlag = FALSE)
#'   result <- mfrm_cluster_imputed(review, model, cells, k = 2)
#'   summary(result)
#'   result$co_membership[1:4, 1:4]
#'   result$imputation_model$loggedEvents
#' }
#' }
#' @export
mfrm_cluster_imputed <- function(x, imputed, impute, k, weights = NULL,
                                 missing = c("error", "omit"),
                                 method = c("pam", "hierarchical", "kmeans"), linkage = NULL,
                                 components = NULL, scale = TRUE, nstart = 25,
                                 iter.max = 100, seed = 1, silhouette = TRUE) {
  method <- match.arg(method)
  if (method != "kmeans" && (!is.null(components) || !missing(scale) ||
      !missing(nstart) || !missing(iter.max) || !missing(seed) || !missing(silhouette))) {
    stop("Numeric preprocessing and initialization arguments require method = 'kmeans'.", call. = FALSE)
  }
  if (method == "hierarchical") {
    linkage <- match.arg(linkage, c("average", "complete"))
  } else if (!is.null(linkage)) {
    stop(if (method == "pam") "PAM has no linkage; `linkage` must be NULL for method = 'pam'." else
      "K-means has no linkage; `linkage` must be NULL for method = 'kmeans'.", call. = FALSE)
  }
  if (!inherits(x, "mfrm_features")) {
    stop("`x` must be prepared with mfrm_features().", call. = FALSE)
  }
  x <- mfrm_features(x$data, x$id, x$features, x$missing)
  missing <- match.arg(missing)
  if (nrow(x$data) > 5000L) {
    stop("Imputation comparisons are limited to 5,000 total entities because the ID-indexed matrix requires quadratic memory.", call. = FALSE)
  }
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Install the optional 'mice' package to use mfrm_cluster_imputed().", call. = FALSE)
  }
  if (!inherits(imputed, "mids") || !is.numeric(imputed$m) || is.complex(imputed$m) ||
      length(imputed$m) != 1L || !is.finite(imputed$m) ||
      imputed$m < 2 || imputed$m != floor(imputed$m)) {
    stop("`imputed` must be a mice mids object with at least two imputations.", call. = FALSE)
  }
  if (!is.data.frame(impute) || anyDuplicated(names(impute)) ||
      !all(c("ID", "Feature") %in% names(impute)) ||
      (!nrow(impute) && nrow(x$missing) > 0L)) {
    stop("`impute` must explicitly list missing cells in ID and Feature columns.", call. = FALSE)
  }
  cells <- as.data.frame(impute[, c("ID", "Feature"), drop = FALSE])
  if (!all(vapply(cells, function(v) is.atomic(v) && is.null(dim(v)), logical(1)))) {
    stop("Imputation cell IDs and features must be vectors of labels.", call. = FALSE)
  }
  cells[] <- lapply(cells, as.character)
  row <- match(cells$ID, x$row_summary$ID)
  column <- match(cells$Feature, x$features)
  absent <- is.na(x$data[x$features])
  if (anyNA(cells) || anyNA(row) || anyNA(column) ||
      anyDuplicated(cells) || !all(absent[cbind(row, column)])) {
    stop("`impute` must list unique missing cells in `x`, identified by ID and Feature.", call. = FALSE)
  }
  eligible <- absent & FALSE
  eligible[cbind(row, column)] <- TRUE
  remaining <- absent & !eligible
  if (any(remaining) && missing == "error") {
    stop("Some missing features are not selected for imputation. Review the selection or use missing = 'omit' explicitly.", call. = FALSE)
  }

  # Validate by ID, not row position; numeric storage mode may change on completion.
  align <- function(data, completion = FALSE) {
    # mice's logistic imputer returns numeric 0/1 for logical columns.
    if (completion) for (name in x$features) {
      value <- data[[name]]
      if (is.logical(x$data[[name]]) && is.numeric(value) && !is.object(value) &&
          !is.complex(value) && is.null(dim(value)) &&
          all(is.na(value) | value %in% c(0, 1))) data[[name]] <- as.logical(value)
    }
    reviewed <- mfrm_features(data, x$id, x$features)
    if (!setequal(reviewed$row_summary$ID, x$row_summary$ID)) {
      stop("Entity IDs must match the original feature table.", call. = FALSE)
    }
    reviewed$data[match(x$row_summary$ID, reviewed$row_summary$ID), , drop = FALSE]
  }
  check_values <- function(data, expected_missing) {
    if (!identical(unname(is.na(data[x$features])), unname(expected_missing))) {
      stop("Missing cells do not match the explicit imputation selection.", call. = FALSE)
    }
    for (name in x$features) {
      before <- x$data[[name]]
      after <- data[[name]]
      numeric_pair <- is.numeric(before) && is.numeric(after)
      if (!(numeric_pair || (identical(class(before), class(after)) &&
                            identical(levels(before), levels(after))))) {
        stop("Feature types and factor levels/order must match the original table.", call. = FALSE)
      }
      observed <- !is.na(before)
      if (!identical(unname(before[observed]), unname(after[observed])) &&
          !(numeric_pair && identical(as.numeric(before[observed]), as.numeric(after[observed])))) {
        stop("Observed feature values must remain unchanged.", call. = FALSE)
      }
    }
  }
  original <- align(imputed$data)
  check_values(original, absent)
  where <- imputed$where
  if (!is.matrix(where) || !is.logical(where) || anyNA(where) ||
      !identical(dim(where), dim(imputed$data)) ||
      !identical(colnames(where), names(imputed$data))) {
    stop("`imputed$where` must identify the cells requested from mice.", call. = FALSE)
  }
  source_order <- match(x$row_summary$ID, as.character(imputed$data[[x$id]]))
  if (any(where[, x$id]) ||
      !identical(unname(where[source_order, x$features, drop = FALSE]), unname(eligible))) {
    stop("`imputed$where` must match `impute` for selected features and must not impute IDs or observed feature values.", call. = FALSE)
  }
  reasons <- x$missing[remaining[cbind(match(x$missing$ID, x$row_summary$ID),
                                      match(x$missing$Feature, x$features))], , drop = FALSE]
  analyses <- vector("list", imputed$m)
  ids <- x$row_summary$ID
  co_membership <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in seq_len(imputed$m)) {
    analyses[[i]] <- tryCatch({
      completed <- align(mice::complete(imputed, action = i), completion = TRUE)
      check_values(completed, remaining)
      reviewed <- mfrm_features(completed, x$id, x$features, reasons)
      if (method == "pam") {
        mfrm_cluster_pam(reviewed, k = k, weights = weights, missing = missing)
      } else if (method == "hierarchical") {
        mfrm_cluster_hierarchical(reviewed, k = k, weights = weights,
          missing = missing, linkage = linkage)
      } else if (is.null(components)) {
        mfrm_cluster_kmeans(reviewed, k, weights, missing, scale,
          nstart, iter.max, seed, silhouette)
      } else {
        pca <- mfrm_pca(reviewed, components, scale, weights, missing)
        mfrm_cluster_kmeans(pca, k, nstart = nstart, iter.max = iter.max,
          seed = seed, silhouette = silhouette)
      }
    }, error = function(e) {
      stop("Imputation ", i, ": ", conditionMessage(e), call. = FALSE)
    })
    labels <- analyses[[i]]$membership$Cluster
    co_membership <- co_membership + outer(labels, labels, "==")
  }
  co_membership <- co_membership / imputed$m
  analysis_summary <- data.frame(
    Imputation = seq_len(imputed$m),
    Included = vapply(analyses, function(a) a$settings$included, integer(1)),
    Excluded = vapply(analyses, function(a) as.integer(a$settings$excluded), integer(1)),
    MeanSilhouette = vapply(analyses, function(a) mean_available_silhouette(a$membership$Silhouette), numeric(1)))
  out <- list(analyses = analyses, co_membership = co_membership,
    analysis_summary = analysis_summary, feature_data = x,
    imputed_cells = dplyr::left_join(cells, x$missing, by = c("ID", "Feature")),
    imputation_model = imputed,
    settings = list(imputations = as.integer(imputed$m), k = as.integer(k),
      method = analyses[[1L]]$settings$method, linkage = linkage,
      weights = analyses[[1L]]$settings$weights, missing = missing,
      components = components, scale = if (method == "kmeans") scale else NULL,
      nstart = if (method == "kmeans") nstart else NULL,
      iter.max = if (method == "kmeans") iter.max else NULL,
      seed = if (method == "kmeans") seed else NULL,
      silhouette = if (method == "kmeans") silhouette else TRUE,
      scaling = if (method == "kmeans")
        "Numeric transformation and any PCA refitted within each completed sample" else
        "Gower ranges recalculated within each completed sample"))
  class(out) <- "mfrm_imputed_clusters"
  out
}

#' @rdname mfrm_cluster_imputed
#' @param object An object returned by [mfrm_cluster_imputed()].
#' @param ... Reserved for method compatibility.
#' @export
print.mfrm_imputed_clusters <- function(x, ...) {
  cat("Exploratory groups across external-feature imputations\n")
  cat("Method:", x$settings$method %||% "PAM",
    if (!is.null(x$settings$linkage)) paste0(" / ", x$settings$linkage, " linkage") else "", "\n")
  cat(x$settings$imputations, "imputations;", nrow(x$imputed_cells), "selected missing cells\n")
  print(x$analysis_summary, row.names = FALSE)
  events <- x$imputation_model$loggedEvents
  cat(if (is.null(events)) 0L else nrow(events), "mice logged events; inspect the retained model and diagnostics.\n")
  cat("Co-membership describes sensitivity to these imputations, not membership probabilities or sampling stability.\n")
  invisible(x)
}

#' @rdname mfrm_cluster_imputed
#' @name summary.mfrm_imputed_clusters
#' @aliases summary.mfrm_imputed_clusters
#' @usage \method{summary}{mfrm_imputed_clusters}(object, ...)
#' @rawNamespace S3method(summary,mfrm_imputed_clusters,summary.mfrm_imputed_fits)
NULL

mean_available_silhouette <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}
