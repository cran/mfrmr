# ==============================================================================
# Import adapters: mirt / TAM / eRm fits to mfrmr-compatible bundles
# ==============================================================================
#
# Importers expose the measurement information that maps reliably to the
# `mfrm_fit` contract: item/step parameters, person scores, and basic fit
# statistics. Bias, DIF, anchor-review, and replay components are outside this
# adapter contract when the source package does not expose the underlying data.
#
# Imported bundles retain `mfrm_fit` inheritance for measurement displays,
# with a distinct first class and source-scale metadata for restricted dispatch.

#' Import an `mirt` fit to an mfrmr-compatible bundle
#'
#' Extracts item, step, and person parameters from a [mirt::mirt()]
#' fit and returns an `mfrm_imported_fit` object. The returned
#' object has the public slots `summary`, `facets$person`,
#' `facets$others`, `steps`, `config`, and `source` that the mfrmr
#' plot and table helpers expect. Only unidimensional Rasch and partial-credit
#' response models with positive slopes and ordinary category scores are
#' supported. Graded-response, guessing and multidimensional models are refused.
#' With `compute_fit = TRUE`, source Infit / Outfit statistics are attached.
#'
#' @param fit An object returned by [mirt::mirt()] (a
#'   `SingleGroupClass`).
#' @param model One of `"RSM"`, `"PCM"`, `"GPCM"`. The importer
#'   does not reconstruct all source constraints; pass the model that was
#'   estimated. Non-unit slopes require `"GPCM"`. A polytomous `"RSM"`
#'   import requires source item type `"rsm"`.
#' @param item_facet Name to assign to the item facet in the
#'   imported bundle (default `"Item"`).
#' @param compute_fit Logical. When `TRUE`, run [mirt::itemfit()]
#'   and [mirt::personfit()] to populate Infit / Outfit / OutfitZSTD
#'   columns on the returned facet tables, plus build a
#'   measurement-side diagnostics bundle. Person fit uses source EAP scores.
#'   Default `FALSE` extracts parameters without calculating fit statistics.
#'
#' @return An `mfrm_imported_fit` object. Slots:
#' \describe{
#'   \item{`summary`}{Model / method / N / LogLik / AIC / BIC.}
#'   \item{`facets$person`}{Person ID, Estimate, SE, Extreme, plus
#'     Infit / Outfit / OutfitZSTD / Zh when `compute_fit = TRUE`.}
#'   \item{`facets$others`}{Item-level estimates and slopes; with
#'     `compute_fit = TRUE`, also available Infit / Outfit statistics.}
#'   \item{`steps`}{Absolute adjacent-category thresholds on the source ability
#'     scale, labelled in `Parameterization`; these are not centered step
#'     deviations. Rating-scale offsets are included.}
#'   \item{`config`}{List with the declared `model` and facet names
#'     used for the import; downstream plot and table helpers consult
#'     this to dispatch correctly on the imported bundle.}
#'   \item{`diagnostics`}{`mfrm_diagnostics`-shape bundle when
#'     `compute_fit = TRUE`; `NULL` otherwise.}
#'   \item{`source`}{Imported-from metadata.}
#' }
#'
#' @section Source scale:
#' Item difficulty is the mean of its absolute adjacent-category thresholds.
#' Source identification and slopes are retained without rescaling. For mirt
#' `gpcmIRT` and `rsm`, the category offset is included as `b - c / a`.
#' Binary and polytomous items may be mixed; thresholds are read using each
#' item's response model and category count. Inapplicable padding columns in
#' the source coefficient table are not missing estimated thresholds.
#' Person estimates are EAP; the `SE` column contains conditional posterior
#' SDs, not sampling SEs. Person labels use retained source row names or
#' P-prefixed row positions. Original identifiers discarded by mirt cannot be
#' recovered. Imported summaries describe these conventions without assuming
#' a native mfrmr population distribution or slope normalization.
#'
#' @section Imported uncertainty:
#' Imported SEs retain the source package's interpretation. The measurement-side
#' diagnostics do not reconstruct the joint parameter covariance, so joint facet
#' chi-square statistics, degrees of freedom and p-values are unavailable.
#' Posterior SDs do not supply sampling SEs for separation reliability.
#' Other separation summaries require valid SEs for every finite estimate and
#' remain descriptive. Imported Wright maps show points only: source uncertainty
#' conventions do not establish one common confidence-interval calculation.
#' Re-import older saved bundles from the existing source-package fit to update
#' difficulties, thresholds and uncertainty labels. The mirt and TAM importers
#' accept `compute_fit = TRUE` when source fit statistics are needed; no model
#' re-estimation is required.
#'
#' @section Scope:
#' Use `summary()` for source-scale tables and `plot()` for a point-only Wright
#' map. Available source fit statistics remain in the facet and diagnostic
#' tables. Native model curves, comprehensive [mfrm_results()] reports,
#' response-level diagnostics, [run_qc_pipeline()], bias/DIF analysis, anchoring
#' and portable calibration are unavailable for imported bundles. This is a
#' one-way fitted-object import of the documented fields.
#' @seealso [import_tam_fit()], [import_erm_fit()]
#' @examples
#' # Optional integration example (loads additional modelling or rendering packages).
#' # example("import_mirt_fit", run.dontrun = TRUE)
#' \dontrun{
#' if (requireNamespace("mirt", quietly = TRUE)) {
#'   response_matrix <- matrix(sample(0:1, 120, replace = TRUE), nrow = 40)
#'   colnames(response_matrix) <- paste0("Item", seq_len(ncol(response_matrix)))
#'   fit <- mirt::mirt(response_matrix, 1, itemtype = "Rasch", verbose = FALSE)
#'   imported <- import_mirt_fit(fit, model = "RSM")
#'   imported$summary
#' }
#' }
#' @export
import_mirt_fit <- function(fit, model = c("RSM", "PCM", "GPCM"),
                             item_facet = "Item",
                             compute_fit = FALSE) {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("`import_mirt_fit()` requires the `mirt` package (suggested).",
         call. = FALSE)
  }
  if (!methods::is(fit, "SingleGroupClass")) {
    stop("`fit` must be an mirt SingleGroupClass result.", call. = FALSE)
  }
  model <- match.arg(model)

  if (!identical(as.integer(mirt::extract.mirt(fit, "nfact")), 1L)) {
    stop("Only unidimensional mirt fits can be imported; additional factor scores are not standard errors.",
         call. = FALSE)
  }
  item_types <- as.character(mirt::extract.mirt(fit, "itemtype"))
  if (!length(item_types) ||
      any(!item_types %in% c("Rasch", "2PL", "gpcm", "gpcmIRT", "rsm"))) {
    stop("Import supports Rasch and partial-credit item models only; graded, guessing and custom response models are not interchangeable with them.",
         call. = FALSE)
  }
  scoring_matrices <- mirt::extract.mirt(fit, "gpcm_mats")
  if (length(scoring_matrices) && any(vapply(scoring_matrices, function(x) {
    !is.null(x) && (ncol(as.matrix(x)) != 1L ||
      !isTRUE(all.equal(as.numeric(x), seq_len(nrow(as.matrix(x))) - 1,
                        check.attributes = FALSE)))
  }, logical(1L)))) {
    stop("Custom mirt category scoring cannot be imported as ordinary partial credit.", call. = FALSE)
  }

  items <- .mirt_extract_items(fit)
  if (model != "GPCM" && any(abs(items$Slope - 1) > 1e-8)) {
    stop("Non-unit source item slopes require model = 'GPCM'; importing does not rescale the source model.", call. = FALSE)
  }
  if (model == "RSM" && any(mirt::extract.mirt(fit, "K") > 2L) &&
      !all(item_types == "rsm")) {
    stop("A polytomous RSM import requires a source rating-scale model; use model = 'PCM' or 'GPCM' for item-specific thresholds.", call. = FALSE)
  }
  steps <- .mirt_extract_steps(items, item_facet = item_facet)
  persons <- .mirt_extract_persons(fit)

  facet_others <- data.frame(
    Facet = item_facet,
    Level = items$Level,
    Estimate = items$Difficulty,
    Slope = items$Slope,
    stringsAsFactors = FALSE
  )

  fit_attached <- list(person = NULL, item = NULL)
  if (isTRUE(compute_fit)) {
    fit_attached <- .mirt_compute_fit_stats(fit, persons)
    if (!is.null(fit_attached$item) && nrow(fit_attached$item) > 0L) {
      m <- match(facet_others$Level, fit_attached$item$Level)
      ok <- !is.na(m)
      attach_cols <- setdiff(names(fit_attached$item), "Level")
      for (col in attach_cols) {
        facet_others[[col]] <- NA_real_
        facet_others[[col]][ok] <- fit_attached$item[[col]][m[ok]]
      }
    }
    if (!is.null(fit_attached$person) && nrow(fit_attached$person) > 0L) {
      m <- match(persons$Person, fit_attached$person$Person)
      ok <- !is.na(m)
      attach_cols <- setdiff(names(fit_attached$person), "Person")
      for (col in attach_cols) {
        persons[[col]] <- NA_real_
        persons[[col]][ok] <- fit_attached$person[[col]][m[ok]]
      }
    }
  }

  summary_tbl <- data.frame(
    Model = model,
    Method = "MML",
    Source = "mirt",
    N = nrow(persons),
    Persons = nrow(persons),
    Facets = 1L,
    Categories = NA_integer_,
    LogLik = as.numeric(fit@Fit$logLik %||% NA_real_),
    AIC = as.numeric(fit@Fit$AIC %||% NA_real_),
    BIC = as.numeric(fit@Fit$BIC %||% NA_real_),
    Converged = isTRUE(fit@OptimInfo$converged %||% NA),
    ConvergenceStatus = if (isTRUE(fit@OptimInfo$converged %||% NA)) "ok" else "review",
    stringsAsFactors = FALSE
  )

  diagnostics <- NULL
  if (isTRUE(compute_fit)) {
    diagnostics <- .synthesize_imported_diagnostics(
      facet_others = facet_others,
      persons = persons,
      facet_names = item_facet,
      n_obs = sum(!is.na(mirt::extract.mirt(fit, "data"))),
      source = "mirt"
    )
  }

  out <- list(
    summary = summary_tbl,
    facets = list(person = persons, others = facet_others),
    steps = steps,
    diagnostics = diagnostics,
    config = list(model = model, method = "MML",
                  facet_names = item_facet,
                  source = "mirt"),
    source = list(package = "mirt",
                  package_version = as.character(utils::packageVersion("mirt")),
                  source_object_class = class(fit)[1],
                  metric_version = 1L,
                  metric = "Item difficulty and adjacent-category thresholds on the source ability scale; source identification and item slopes retained.",
                  person_scoring = "EAP with conditional posterior SD; calibration uncertainty excluded.",
                  person_identification = "Person labels use retained source row names or P-prefixed row positions; original identifiers discarded by mirt cannot be recovered.",
                  item_types = item_types,
                  compute_fit = isTRUE(compute_fit))
  )
  class(out) <- c("mfrm_imported_fit", "mfrm_fit", "list")
  out
}

# --- mirt internal extractors --------------------------------------------

.mirt_extract_items <- function(fit) {
  irt_pars <- tryCatch(
    .mfrm_mirt_coef(fit, simplify = TRUE, IRTpars = TRUE)$items,
    error = function(e) NULL
  )
  if (is.null(irt_pars) || !is.matrix(irt_pars)) {
    stop("Source IRT difficulties could not be extracted; raw intercepts cannot be substituted for thresholds.", call. = FALSE)
  }
  items_mat <- as.data.frame(irt_pars, stringsAsFactors = FALSE)
  items_mat$Level <- rownames(items_mat)
  rownames(items_mat) <- NULL
  # Keep the source discrimination without normalizing it to a native scale.
  slope_cols <- intersect(c("a", "a1"), names(items_mat))
  items_mat$Slope <- vapply(seq_len(nrow(items_mat)), function(i) {
    values <- as.numeric(unlist(items_mat[i, slope_cols, drop = FALSE], use.names = FALSE))
    values <- values[is.finite(values)]
    if (length(values) == 1L) values else NA_real_
  }, numeric(1L))
  if (any(!is.finite(items_mat$Slope) | items_mat$Slope <= 0)) {
    stop("Every imported item must have a finite positive source slope.", call. = FALSE)
  }
  counts <- as.integer(mirt::extract.mirt(fit, "K")) - 1L
  types <- as.character(mirt::extract.mirt(fit, "itemtype"))
  if (length(counts) != nrow(items_mat) || length(types) != nrow(items_mat) ||
      anyNA(counts) || any(counts < 1L)) {
    stop("Source item category counts cannot be aligned with their thresholds.", call. = FALSE)
  }
  # Mixed source tables use b for dichotomous items and b1, b2, ... for
  # partial-credit items; their other columns contain inapplicable NA padding.
  items_mat$Thresholds <- lapply(seq_len(nrow(items_mat)), function(i) {
    columns <- if (types[i] %in% c("Rasch", "2PL") && counts[i] == 1L) "b" else
      paste0("b", seq_len(counts[i]))
    if (!all(columns %in% names(items_mat)))
      stop("Source adjacent-category thresholds are unavailable.", call. = FALSE)
    values <- as.numeric(items_mat[i, columns, drop = TRUE])
    if (types[i] %in% c("gpcmIRT", "rsm")) {
      if (!"c" %in% names(items_mat) || !is.finite(items_mat$c[i]))
        stop("Source rating-scale offsets are unavailable.", call. = FALSE)
      values <- values - items_mat$c[i] / items_mat$Slope[i]
    }
    if (any(!is.finite(values)))
      stop("Every imported item must have finite source adjacent-category thresholds.", call. = FALSE)
    values
  })
  items_mat$Difficulty <- vapply(items_mat$Thresholds, mean, numeric(1L))
  if (any(!is.finite(items_mat$Difficulty))) {
    stop("Every imported item must have finite source adjacent-category thresholds.", call. = FALSE)
  }
  items_mat
}

.mirt_extract_steps <- function(items_mat, item_facet) {
  rows <- lapply(seq_len(nrow(items_mat)), function(i) {
    data.frame(
      StepFacet = item_facet,
      Level = items_mat$Level[i],
      Step = seq_along(items_mat$Thresholds[[i]]),
      Estimate = items_mat$Thresholds[[i]],
      Parameterization = "Absolute adjacent-category threshold on source ability scale",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  # Preserve the historical step-major ordering for homogeneous imports.
  out <- out[order(out$Step, match(out$Level, items_mat$Level)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.mirt_extract_persons <- function(fit) {
  fscores <- mirt::fscores(fit, method = "EAP", full.scores.SE = TRUE, verbose = FALSE)
  if (ncol(fscores) != 2L || !startsWith(colnames(fscores)[2L], "SE_")) {
    stop("The source EAP score and posterior-SD columns cannot be identified safely.", call. = FALSE)
  }
  ids <- rownames(mirt::extract.mirt(fit, "data"))
  if (is.null(ids) || identical(ids, as.character(seq_len(nrow(fscores))))) {
    ids <- paste0("P", seq_len(nrow(fscores)))
  }
  person_tbl <- data.frame(
    Person = ids,
    Estimate = as.numeric(fscores[, 1]),
    SE = if (ncol(fscores) >= 2L) as.numeric(fscores[, 2]) else NA_real_,
    Extreme = "none",
    EstimateBasis = "Source EAP",
    UncertaintyBasis = "Conditional posterior SD; calibration uncertainty excluded",
    stringsAsFactors = FALSE
  )
  person_tbl
}

.mirt_compute_fit_stats <- function(fit, persons) {
  item_fit <- tryCatch(
    suppressMessages(suppressWarnings(
      mirt::itemfit(fit, fit_stats = c("infit"),
                     na.rm = TRUE)
    )),
    error = function(e) NULL
  )
  person_fit <- tryCatch(
    suppressMessages(suppressWarnings(
      mirt::personfit(fit, method = "EAP")
    )),
    error = function(e) NULL
  )
  item_df <- if (!is.null(item_fit)) {
    df <- as.data.frame(item_fit, stringsAsFactors = FALSE)
    keep <- intersect(c("item", "infit", "outfit", "z.infit", "z.outfit",
                          "S_X2", "df.S_X2", "RMSEA.S_X2", "p.S_X2"),
                       names(df))
    df <- df[, keep, drop = FALSE]
    if ("item" %in% names(df)) names(df)[names(df) == "item"] <- "Level"
    if ("infit" %in% names(df)) names(df)[names(df) == "infit"] <- "Infit"
    if ("outfit" %in% names(df)) names(df)[names(df) == "outfit"] <- "Outfit"
    if ("z.infit" %in% names(df)) names(df)[names(df) == "z.infit"] <- "InfitZSTD"
    if ("z.outfit" %in% names(df)) names(df)[names(df) == "z.outfit"] <- "OutfitZSTD"
    df
  } else NULL
  person_df <- if (!is.null(person_fit)) {
    df <- as.data.frame(person_fit, stringsAsFactors = FALSE)
    df$Person <- persons$Person
    keep <- intersect(c("Person", "infit", "outfit", "z.infit", "z.outfit", "Zh"),
                       names(df))
    df <- df[, keep, drop = FALSE]
    if ("infit" %in% names(df)) names(df)[names(df) == "infit"] <- "Infit"
    if ("outfit" %in% names(df)) names(df)[names(df) == "outfit"] <- "Outfit"
    if ("z.infit" %in% names(df)) names(df)[names(df) == "z.infit"] <- "InfitZSTD"
    if ("z.outfit" %in% names(df)) names(df)[names(df) == "z.outfit"] <- "OutfitZSTD"
    df
  } else NULL
  list(item = item_df, person = person_df)
}

#' Import a `TAM` fit to an mfrmr-compatible bundle
#'
#' Extracts item / step / person parameters from a unidimensional
#' [TAM::tam.mml()] or [TAM::tam.mml.mfr()] fit. Difficulties and absolute
#' adjacent-category thresholds are derived from the source category logits,
#' rather than inferred from coefficient names. Source slopes and scale
#' identification are retained.
#'
#' @param fit An object returned by `TAM::tam.mml()` or
#'   `TAM::tam.mml.mfr()`.
#' @param model Declared response model: `"RSM"`, `"PCM"`, or `"GPCM"`.
#'   Non-unit source slopes require `"GPCM"`; import does not reconstruct
#'   all source constraints or certify equivalence to a native mfrmr model.
#' @param item_facet Name to assign to the item facet for the
#'   single-facet path. Ignored when the input is a multi-facet
#'   `tam.mml.mfr` fit, whose combined response conditions are labelled
#'   `"DesignCell"`.
#' @param compute_fit Logical. When `TRUE`, run [TAM::msq.itemfit()]
#'   and [TAM::tam.personfit()] to populate Infit / Outfit columns
#'   on the returned facet tables, plus build a measurement-side
#'   `mfrm_diagnostics` bundle. Default `FALSE`.
#'
#' @details
#' Each item difficulty is the mean of its absolute adjacent-category thresholds.
#' Positive constant adjacent-category slopes are required. For multi-facet
#' fits, each returned difficulty combines all facet effects for that response
#' condition. Separate facet coordinates are not reconstructed; the original
#' coefficient table remains in `source$native_parameters`.
#'
#' A transformed item SE is retained only when the location is a scalar multiple
#' of one source coefficient and its slope is fixed. Otherwise it is missing:
#' marginal coefficient SEs cannot replace the required joint covariance.
#' Persons retain source EAP and conditional posterior SD. Item fit is averaged
#' over source posteriors; `TAM::tam.personfit()` uses WLE scores. `FitBasis`
#' records this distinction without replacing the imported EAP estimates.
#' Person `Extreme` flags use observed responses and each item's fitted category
#' support, including unobserved categories. `NObserved` counts responses used
#' for this check. `ResponseStatus` distinguishes observed, entirely missing,
#' and unavailable response information. Entirely missing or unverifiable rows
#' have `Extreme = NA`; finite posterior scores do not imply nonextreme data.
#' Re-import older TAM objects to recover these flags; saved imported objects
#' without the source responses cannot reconstruct them.
#'
#' The public imported-fit surface is deliberately unidimensional and MML-only.
#' A `tam.jml` object is not silently relabelled as MML, and a TAM fit with
#' `ndim > 1` is rejected rather than flattened into one mfrmr scale. Keep
#' multidimensional TAM fits in a separate validation workflow.
#'
#' TAM-native AIC, BIC, and adjusted BIC are retained as explicitly named
#' `Native*` fields. Compatibility `AIC` and `BIC` columns still mirror the
#' native TAM values, but the imported object has
#' `ICStatus = "imported_native_descriptive"`, `ICEligible = FALSE`, and cannot
#' enter [compare_mfrm()] as a current mfrmr IC contract. In particular, TAM's
#' native `aBIC` is not relabelled as the package's Sclove `SABIC`. The source
#' metadata also retains the TAM version, dimension count, iterations, and
#' iteration ceiling used for the conservative imported convergence status.
#'
#' @inheritSection import_mirt_fit Imported uncertainty
#' @inheritSection import_mirt_fit Scope
#'
#' @return An `mfrm_imported_fit` object. Slots mirror
#'   [import_mirt_fit()], with explicit TAM-native IC provenance in `summary`
#'   and `source`.
#' @seealso [import_mirt_fit()], [import_erm_fit()]
#' @examples
#' if (requireNamespace("TAM", quietly = TRUE)) {
#'   response_matrix <- matrix(sample(0:3, 60, replace = TRUE), nrow = 20)
#'   colnames(response_matrix) <- paste0("Item", seq_len(ncol(response_matrix)))
#'   fit <- TAM::tam.mml(resp = response_matrix, irtmodel = "PCM")
#'   imported <- import_tam_fit(fit, model = "PCM")
#'   imported$summary
#' }
#' @export
import_tam_fit <- function(fit, model = c("RSM", "PCM", "GPCM"),
                            item_facet = "Item",
                            compute_fit = FALSE) {
  if (!requireNamespace("TAM", quietly = TRUE)) {
    stop("`import_tam_fit()` requires the `TAM` package (suggested).",
         call. = FALSE)
  }
  import_scope <- .tam_validate_import_scope(fit)
  if (!is.list(fit) || is.null(fit$xsi) || is.null(fit$person)) {
    stop("`fit` does not look like a TAM result (missing `xsi` / `person`).",
         call. = FALSE)
  }
  model <- match.arg(model)

  # `tam.mml.mfr()` returns an object with class `tam.mml` (no
  # dedicated subclass) but with an `xsi.facets` slot that identifies
  # the facet of every fitted parameter. Use that slot's presence as
  # the multi-facet detector.
  is_mfr <- !is.null(fit$xsi.facets) &&
            is.data.frame(fit$xsi.facets) &&
            nrow(fit$xsi.facets) > 0L &&
            "facet" %in% names(fit$xsi.facets) &&
            length(unique(as.character(fit$xsi.facets$facet))) >= 1L

  persons <- .tam_extract_persons(fit)

  if (is_mfr) {
    extracted <- .tam_extract_single(fit, item_facet = "DesignCell")
    facet_others <- extracted$others
    steps <- extracted$steps
    facet_names <- extracted$facet_names
  } else {
    extracted <- .tam_extract_single(fit, item_facet = item_facet)
    facet_others <- extracted$others
    steps <- extracted$steps
    facet_names <- item_facet
  }
  if (model != "GPCM" && any(abs(facet_others$Slope - 1) > 1e-8)) {
    stop("Non-unit TAM item slopes require model = 'GPCM'; importing does not rescale the source model.", call. = FALSE)
  }

  fit_attached <- list(person = NULL, item = NULL)
  if (isTRUE(compute_fit)) {
    fit_attached <- .tam_compute_fit_stats(fit, persons = persons)
    if (!is.null(fit_attached$item) && nrow(fit_attached$item) > 0L &&
        nrow(facet_others) > 0L) {
      m <- match(facet_others$Level, fit_attached$item$Level)
      ok <- !is.na(m)
      attach_cols <- setdiff(names(fit_attached$item), "Level")
      for (col in attach_cols) {
        facet_others[[col]] <- NA_real_
        facet_others[[col]][ok] <- fit_attached$item[[col]][m[ok]]
      }
    }
    if (!is.null(fit_attached$person) && nrow(fit_attached$person) > 0L) {
      m <- match(persons$Person, fit_attached$person$Person)
      ok <- !is.na(m)
      attach_cols <- setdiff(names(fit_attached$person), "Person")
      for (col in attach_cols) {
        persons[[col]] <- NA_real_
        persons[[col]][ok] <- fit_attached$person[[col]][m[ok]]
      }
    }
  }

  native_ic <- .tam_native_ic_record(fit, persons = persons)
  convergence <- .tam_import_convergence_state(fit)
  summary_tbl <- data.frame(
    Model = model,
    Method = "MML",
    Source = "TAM",
    N = nrow(persons),
    Persons = nrow(persons),
    Facets = length(unique(facet_others$Facet)),
    Categories = NA_integer_,
    LogLik = native_ic$NativeLogLik,
    AIC = native_ic$NativeAIC,
    BIC = native_ic$NativeBIC,
    ICContractVersion = "external_native_tam_v1",
    ICEligible = FALSE,
    ICSelectable = FALSE,
    ICStatus = "imported_native_descriptive",
    NativeDeviance = native_ic$NativeDeviance,
    NativeLogLik = native_ic$NativeLogLik,
    NativeNpar = native_ic$NativeNpar,
    NativeICSampleSize = native_ic$NativeICSampleSize,
    NativeAIC = native_ic$NativeAIC,
    NativeBIC = native_ic$NativeBIC,
    NativeABIC = native_ic$NativeABIC,
    NativeAICFormula = "tam_deviance_plus_2k",
    NativeBICFormula = "tam_deviance_plus_log_n_k",
    NativeABICFormula = "tam_deviance_plus_log_n_minus_2_over_24_k",
    NativeAICFormulaVerified = native_ic$NativeAICFormulaVerified,
    NativeBICFormulaVerified = native_ic$NativeBICFormulaVerified,
    NativeABICFormulaVerified = native_ic$NativeABICFormulaVerified,
    NativeLogLikDevianceConsistent = native_ic$NativeLogLikDevianceConsistent,
    NativePersonCountConsistent = native_ic$NativePersonCountConsistent,
    Converged = convergence$Converged,
    ConvergenceStatus = convergence$Status,
    Iterations = convergence$Iterations,
    IterationCeiling = convergence$IterationCeiling,
    stringsAsFactors = FALSE
  )

  diagnostics <- NULL
  if (isTRUE(compute_fit)) {
    diagnostics <- .synthesize_imported_diagnostics(
      facet_others = facet_others,
      persons = persons,
      facet_names = facet_names,
      n_obs = sum(persons$NObserved),
      source = "TAM"
    )
  }

  out <- list(
    summary = summary_tbl,
    facets = list(person = persons, others = facet_others),
    steps = steps,
    diagnostics = diagnostics,
    config = list(model = model, method = "MML",
                  facet_names = facet_names,
                  source = "TAM",
                  ndim = import_scope$Dimensions,
                  multi_facet = is_mfr),
    source = list(package = "TAM",
                  package_version = native_ic$TAMVersion,
                  source_object_class = class(fit)[1],
                  dimensions = import_scope$Dimensions,
                  metric_version = 1L,
                  metric = if (is_mfr) {
                    "Response-design cell difficulties include all source facet effects; separate facet coordinates are not reconstructed."
                  } else {
                    "Mean adjacent-category thresholds on the source ability scale; source identification and item slopes retained."
                  },
                  person_scoring = "Source EAP with conditional posterior SD; calibration uncertainty excluded.",
                  native_parameters = fit$xsi.facets %||% fit$xsi,
                  ic_contract = "external_native_tam_v1",
                  convergence = convergence,
                  multi_facet = is_mfr,
                  compute_fit = isTRUE(compute_fit))
  )
  class(out) <- c("mfrm_imported_fit", "mfrm_fit", "list")
  out
}

# --- TAM internal extractors ---------------------------------------------

.tam_validate_import_scope <- function(fit) {
  if (inherits(fit, "tam.jml")) {
    stop(
      "`import_tam_fit()` does not import `tam.jml` objects as MML. Keep TAM JML output in an explicitly JML-labelled external workflow.",
      call. = FALSE
    )
  }
  if (!inherits(fit, "tam.mml")) {
    stop("`fit` must inherit from `tam.mml`.", call. = FALSE)
  }
  dimensions_raw <- suppressWarnings(as.numeric(fit$ndim %||% NA_real_)[1])
  if (!is.finite(dimensions_raw) && !is.null(fit$beta)) {
    beta_dim <- dim(fit$beta)
    if (length(beta_dim) == 2L) dimensions_raw <- beta_dim[2]
  }
  if (!is.finite(dimensions_raw) || dimensions_raw < 1 ||
      abs(dimensions_raw - round(dimensions_raw)) >
        sqrt(.Machine$double.eps)) {
    stop("The TAM latent dimension could not be verified.", call. = FALSE)
  }
  dimensions <- as.integer(round(dimensions_raw))
  if (dimensions != 1L) {
    stop(
      "`import_tam_fit()` supports unidimensional TAM MML only; `ndim = ",
      dimensions,
      "` must remain in a separate multidimensional validation workflow and cannot be exposed as one mfrmr scale.",
      call. = FALSE
    )
  }
  list(Method = "MML", Dimensions = dimensions)
}

.tam_ic_scalar <- function(ic, names, fallback = NA_real_) {
  if (is.null(ic)) return(as.numeric(fallback)[1])
  for (name in names) {
    if (name %in% colnames(ic)) {
      return(suppressWarnings(as.numeric(ic[[name]][1])))
    }
  }
  as.numeric(fallback)[1]
}

.tam_native_ic_record <- function(fit,
                                  persons,
                                  tolerance = 1e-8) {
  ic <- tryCatch(
    as.data.frame(fit$ic, stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  deviance <- .tam_ic_scalar(ic, c("deviance", "Deviance"), fit$deviance)
  loglik <- .tam_ic_scalar(ic, c("loglike", "logLik"), -deviance / 2)
  npar <- .tam_ic_scalar(ic, c("np", "Npars"))
  native_n <- .tam_ic_scalar(ic, "n", nrow(persons))
  native_aic <- .tam_ic_scalar(ic, "AIC")
  native_bic <- .tam_ic_scalar(ic, "BIC")
  native_abic <- .tam_ic_scalar(ic, "aBIC")
  matches <- function(actual, expected) {
    is.finite(actual) && is.finite(expected) &&
      abs(actual - expected) <= tolerance * max(1, abs(expected))
  }
  expected_aic <- deviance + 2 * npar
  expected_bic <- deviance + log(native_n) * npar
  expected_abic <- if (is.finite(native_n) && native_n > 2) {
    deviance + log((native_n - 2) / 24) * npar
  } else {
    NA_real_
  }
  list(
    TAMVersion = as.character(utils::packageVersion("TAM")),
    NativeDeviance = deviance,
    NativeLogLik = loglik,
    NativeNpar = npar,
    NativeICSampleSize = native_n,
    NativeAIC = native_aic,
    NativeBIC = native_bic,
    NativeABIC = native_abic,
    NativeAICFormulaVerified = matches(native_aic, expected_aic),
    NativeBICFormulaVerified = matches(native_bic, expected_bic),
    NativeABICFormulaVerified = matches(native_abic, expected_abic),
    NativeLogLikDevianceConsistent = matches(deviance, -2 * loglik),
    NativePersonCountConsistent = is.finite(native_n) &&
      identical(as.integer(round(native_n)), as.integer(nrow(persons)))
  )
}

.tam_import_convergence_state <- function(fit) {
  iterations <- suppressWarnings(as.integer(fit$iter %||% NA_integer_)[1])
  max_iterations <- suppressWarnings(as.integer(
    fit$control$maxiter %||% NA_integer_
  )[1])
  if (is.finite(iterations) && is.finite(max_iterations)) {
    converged <- iterations < max_iterations
    return(list(
      Converged = converged,
      Status = if (converged) {
        "imported_tam_stopped_before_iteration_ceiling"
      } else {
        "imported_tam_iteration_ceiling_or_convergence_unverified"
      },
      Iterations = iterations,
      IterationCeiling = max_iterations
    ))
  }
  list(
    Converged = NA,
    Status = "imported_tam_convergence_unverified",
    Iterations = iterations,
    IterationCeiling = max_iterations
  )
}

.tam_extract_persons <- function(fit) {
  person_in <- as.data.frame(fit$person, stringsAsFactors = FALSE)
  person_id <- if ("pid" %in% names(person_in)) {
    as.character(person_in$pid)
  } else {
    paste0("P", seq_len(nrow(person_in)))
  }
  person_eap <- if ("EAP" %in% names(person_in)) person_in$EAP else NA_real_
  person_se <- if ("SD.EAP" %in% names(person_in)) person_in$SD.EAP else NA_real_
  data.frame(
    Person = person_id,
    Estimate = as.numeric(person_eap),
    SE = as.numeric(person_se),
    .tam_person_response_status(fit, person_id),
    EstimateBasis = "Source EAP",
    UncertaintyBasis = "Conditional posterior SD; calibration uncertainty excluded",
    stringsAsFactors = FALSE
  )
}

.tam_person_response_status <- function(fit, person_id) {
  out <- data.frame(NObserved = rep(NA_integer_, length(person_id)),
    Extreme = rep(NA_character_, length(person_id)),
    ResponseStatus = rep("unavailable", length(person_id)), stringsAsFactors = FALSE)
  if (is.null(fit$resp) || is.null(fit$AXsi)) return(out)
  response <- as.matrix(fit$resp)
  intercepts <- as.matrix(fit$AXsi)
  ids <- as.character(fit$pid)
  items <- as.character(fit$item$item)
  if (!(is.numeric(response) || is.logical(response)) || !is.numeric(intercepts) ||
      nrow(response) != length(ids) || anyNA(ids) || anyDuplicated(ids) ||
      anyNA(person_id) || anyDuplicated(person_id) ||
      !all(person_id %in% ids) || ncol(response) != length(items) ||
      nrow(intercepts) != length(items) || anyNA(items) || anyDuplicated(items) ||
      is.null(colnames(response)) || anyDuplicated(colnames(response)) ||
      !setequal(colnames(response), items)) return(out)
  support <- is.finite(intercepts)
  maximum <- rowSums(support) - 1L
  if (any(maximum < 1L) ||
      !all(support == outer(maximum, seq_len(ncol(intercepts)) - 1L, `>=`))) return(out)
  observed <- fit$resp.ind
  if (is.null(observed)) observed <- !is.na(response)
  if (!is.matrix(observed) || !identical(dim(observed), dim(response)) ||
      anyNA(observed) || !all(observed %in% c(0, 1))) return(out)
  row <- match(person_id, ids)
  column <- match(items, colnames(response))
  response <- response[row, column, drop = FALSE]
  observed <- observed[row, column, drop = FALSE] != 0
  ceiling <- matrix(maximum, nrow(response), ncol(response), byrow = TRUE)
  valid <- is.finite(response) & response >= 0 & response <= ceiling &
    response == floor(response)
  verified <- rowSums(observed & !valid) == 0L
  count <- rowSums(observed)
  out$NObserved[verified] <- as.integer(count[verified])
  out$ResponseStatus[verified & count == 0L] <- "no_observed_responses"
  used <- verified & count > 0L
  out$ResponseStatus[used] <- "observed"
  out$Extreme[used] <- "none"
  out$Extreme[used & rowSums(observed & (is.na(response) | response != 0)) == 0L] <- "low"
  out$Extreme[used & rowSums(observed & (is.na(response) | response != ceiling)) == 0L] <- "high"
  out
}

.tam_extract_single <- function(fit, item_facet) {
  intercepts <- as.matrix(fit$AXsi)
  loadings <- fit$B
  design <- fit$A
  item_ids <- as.character(fit$item$item)
  if (!is.numeric(intercepts) || !is.array(loadings) ||
      length(dim(loadings)) != 3L || dim(loadings)[3L] != 1L ||
      !identical(dim(intercepts), dim(loadings)[1:2]) ||
      !is.array(design) || length(dim(design)) != 3L ||
      !identical(dim(design)[1:2], dim(intercepts)) ||
      length(item_ids) != nrow(intercepts) || anyNA(item_ids) ||
      any(!nzchar(item_ids)) || anyDuplicated(item_ids) ||
      ncol(intercepts) < 2L || ncol(as.matrix(fit$xsi)) < 1L ||
      dim(design)[3L] != nrow(fit$xsi)) {
    stop("The TAM item labels, category intercepts, loadings and parameter design are not aligned; a common-scale import is unavailable.", call. = FALSE)
  }
  estimates <- standard_errors <- slopes <- rep(NA_real_, length(item_ids))
  step_rows <- vector("list", length(item_ids))
  for (i in seq_along(item_ids)) {
    available <- which(is.finite(intercepts[i, ]))
    if (length(available) < 2L ||
        !identical(available, seq_len(length(available)))) {
      stop("Each TAM item must retain consecutive categories starting at zero.", call. = FALSE)
    }
    count <- length(available) - 1L
    slope_steps <- diff(loadings[i, available, 1L])
    slope <- slope_steps[1L]
    if (!all(is.finite(slope_steps)) || slope <= 0 ||
        any(abs(slope_steps - slope) > 1e-8)) {
      stop("TAM import requires a positive constant adjacent-category slope per item; nominal or custom category scoring cannot be treated as partial credit.", call. = FALSE)
    }
    thresholds <- -diff(intercepts[i, available]) / slope
    estimates[i] <- mean(thresholds)
    slopes[i] <- slope
    # Without joint covariance, a transformed SE is retained only for a
    # single source coefficient and a fixed slope.
    weights <- -(design[i, count + 1L, ] - design[i, 1L, ]) / (count * slope)
    used <- which(abs(weights) > 1e-12)
    slope_se <- fit$se.B
    slope_fixed <- is.array(slope_se) && identical(dim(slope_se), dim(loadings)) &&
      all(is.finite(slope_se[i, available, 1L])) &&
      all(slope_se[i, available, 1L] == 0)
    if (length(used) == 1L && slope_fixed &&
        "se.xsi" %in% names(fit$xsi)) {
      se <- as.numeric(fit$xsi$se.xsi[used])
      if (is.finite(se) && se >= 0) standard_errors[i] <- abs(weights[used]) * se
    }
    step_rows[[i]] <- data.frame(
      StepFacet = item_facet, Level = item_ids[i], Step = seq_len(count),
      Estimate = thresholds,
      Parameterization = "Absolute adjacent-category threshold on source ability scale",
      stringsAsFactors = FALSE
    )
  }
  list(
    others = data.frame(Facet = item_facet, Level = item_ids,
      Estimate = estimates, SE = standard_errors, Slope = slopes,
      stringsAsFactors = FALSE),
    steps = do.call(rbind, step_rows), facet_names = item_facet
  )
}

.tam_compute_fit_stats <- function(fit, persons) {
  item_fit <- tryCatch(
    suppressMessages(suppressWarnings(TAM::msq.itemfit(fit))),
    error = function(e) NULL
  )
  person_fit <- tryCatch(
    suppressMessages(suppressWarnings(TAM::tam.personfit(fit))),
    error = function(e) NULL
  )
  item_df <- if (!is.null(item_fit)) {
    raw <- as.data.frame(item_fit$itemfit %||% item_fit, stringsAsFactors = FALSE)
    if ("parameter" %in% names(raw)) {
      raw$Level <- as.character(raw$parameter)
    } else if ("item" %in% names(raw)) {
      raw$Level <- as.character(raw$item)
    } else {
      raw$Level <- rownames(raw)
    }
    keep <- intersect(c("Level", "Outfit", "Outfit_t", "Infit", "Infit_t"),
                       names(raw))
    if (length(keep) == 0L) return(NULL)
    df <- raw[, keep, drop = FALSE]
    if ("Outfit_t" %in% names(df)) names(df)[names(df) == "Outfit_t"] <- "OutfitZSTD"
    if ("Infit_t" %in% names(df)) names(df)[names(df) == "Infit_t"] <- "InfitZSTD"
    df$FitBasis <- "Source posterior-averaged item fit"
    df
  } else NULL
  person_df <- if (!is.null(person_fit)) {
    raw <- as.data.frame(person_fit, stringsAsFactors = FALSE)
    raw$Person <- if ("pid" %in% names(raw)) as.character(raw$pid) else persons$Person[seq_len(nrow(raw))]
    keep <- intersect(c("Person", "outfitPerson", "infitPerson",
                          "outfitPerson_t", "infitPerson_t"),
                       names(raw))
    if (length(keep) == 0L) return(NULL)
    df <- raw[, keep, drop = FALSE]
    if ("outfitPerson" %in% names(df)) names(df)[names(df) == "outfitPerson"] <- "Outfit"
    if ("infitPerson" %in% names(df)) names(df)[names(df) == "infitPerson"] <- "Infit"
    if ("outfitPerson_t" %in% names(df)) names(df)[names(df) == "outfitPerson_t"] <- "OutfitZSTD"
    if ("infitPerson_t" %in% names(df)) names(df)[names(df) == "infitPerson_t"] <- "InfitZSTD"
    df$FitBasis <- "Source WLE person fit; estimates remain EAP"
    df
  } else NULL
  list(item = item_df, person = person_df)
}

#' Import an `eRm` fit to an mfrmr-compatible bundle
#'
#' Extracts item / person parameters from an [eRm::PCM()],
#' [eRm::RM()] or [eRm::RSM()] fit. Source cumulative easiness coefficients
#' are converted to absolute adjacent-category difficulties. One item location
#' is returned per item: the mean of its thresholds. Source identification is
#' retained. Current `eRm` person tables use `Person Parameter` and
#' `Std.Error`; historical `theta` / `thetapar` estimate labels are also
#' accepted. Unknown or internally misaligned person-table schemas stop with an
#' explicit error rather than silently recycling rows. Same caveats as
#' [import_mirt_fit()].
#'
#' @param fit An object returned by `eRm::PCM()`, `eRm::RM()`, or
#'   `eRm::RSM()`.
#' @param model Matching `"PCM"` or `"RSM"` label; either is accepted for a
#'   binary `RM`. `"GPCM"` and linear extensions are unsupported.
#' @param item_facet Name to assign to the item facet.
#'
#' @return An `mfrm_imported_fit` object.
#' @details
#' Item-location SEs use the corresponding scalar transformation of the source
#' cumulative coefficient SE. Person maximum-likelihood estimates and conditional
#' SEs retain source conventions, including labelled extreme-score extrapolations.
#' The original coefficient table is retained in `source$native_parameters`.
#' @inheritSection import_mirt_fit Imported uncertainty
#' @inheritSection import_mirt_fit Scope
#' @seealso [import_mirt_fit()], [import_tam_fit()]
#' @examples
#' if (requireNamespace("eRm", quietly = TRUE)) {
#'   response_matrix <- matrix(sample(0:3, 60, replace = TRUE), nrow = 20)
#'   colnames(response_matrix) <- paste0("Item", seq_len(ncol(response_matrix)))
#'   fit <- eRm::PCM(response_matrix)
#'   imported <- import_erm_fit(fit, model = "PCM")
#'   imported$summary
#' }
#' @export
import_erm_fit <- function(fit, model = c("RSM", "PCM", "GPCM"),
                            item_facet = "Item") {
  if (!requireNamespace("eRm", quietly = TRUE)) {
    stop("`import_erm_fit()` requires the `eRm` package (suggested).",
         call. = FALSE)
  }
  if (!is.list(fit) || is.null(fit$betapar)) {
    stop("`fit` does not look like an eRm result (missing `betapar`).",
         call. = FALSE)
  }
  model <- match.arg(model)
  if (!inherits(fit, "eRm") || !fit$model %in% c("RM", "PCM", "RSM") ||
      model == "GPCM" || (fit$model == "PCM" && model != "PCM") ||
      (fit$model == "RSM" && model != "RSM")) {
    stop("eRm import requires an RM, PCM or RSM fit and its matching model label; linear extensions and GPCM are not supported.", call. = FALSE)
  }
  beta <- fit$betapar
  counts <- apply(fit$X, 2L, max, na.rm = TRUE)
  if (any(!is.finite(counts) | counts < 1 | counts != floor(counts)) ||
      sum(counts) != length(beta) || is.null(colnames(fit$X))) {
    stop("The eRm cumulative category coefficients cannot be aligned to source items.", call. = FALSE)
  }
  last <- cumsum(counts)
  se_beta <- as.numeric(fit$se.beta %||% rep(NA_real_, length(beta)))
  if (length(se_beta) != length(beta)) {
    stop("The eRm coefficient and standard-error rows are not aligned.", call. = FALSE)
  }
  facet_others <- data.frame(
    Facet = item_facet,
    Level = colnames(fit$X),
    Estimate = -as.numeric(beta[last]) / counts,
    SE = se_beta[last] / counts,
    stringsAsFactors = FALSE
  )
  steps <- do.call(rbind, lapply(seq_along(counts), function(i) {
    positions <- seq.int(last[i] - counts[i] + 1L, last[i])
    data.frame(StepFacet = item_facet, Level = colnames(fit$X)[i],
      Step = seq_len(counts[i]), Estimate = -diff(c(0, beta[positions])),
      Parameterization = "Absolute adjacent-category threshold on source ability scale",
      stringsAsFactors = FALSE)
  }))
  pp <- tryCatch(eRm::person.parameter(fit), error = function(e) NULL)
  if (!is.null(pp)) {
    person_tbl <- .erm_extract_person_table(pp)
  } else {
    person_tbl <- data.frame(Person = character(0), Estimate = numeric(0),
                              SE = numeric(0), Extreme = character(0),
                              stringsAsFactors = FALSE)
  }
  summary_tbl <- data.frame(
    Model = model,
    Method = "CML",
    Source = "eRm",
    N = nrow(person_tbl),
    Persons = nrow(person_tbl),
    Facets = 1L,
    Categories = NA_integer_,
    LogLik = as.numeric(fit$loglik %||% NA_real_),
    AIC = NA_real_,
    BIC = NA_real_,
    Converged = if (length(fit$convergence) == 1L && !is.na(fit$convergence)) fit$convergence %in% c(1L, 2L) else NA,
    ConvergenceStatus = "imported",
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_tbl,
    facets = list(person = person_tbl, others = facet_others),
    steps = steps,
    config = list(model = model, method = "CML",
                  facet_names = item_facet,
                  source = "eRm"),
    source = list(package = "eRm",
                  package_version = as.character(utils::packageVersion("eRm")),
                  source_object_class = class(fit)[1],
                  metric_version = 1L,
                  metric = "Item difficulties and adjacent-category thresholds on the source ability scale; easiness signs and cumulative category coefficients converted.",
                  person_scoring = "Source maximum likelihood; extreme-response extrapolations retain their labels. Person SEs condition on the fitted items.",
                  native_parameters = data.frame(Parameter = names(beta),
                    Estimate = as.numeric(beta), SE = se_beta))
  )
  class(out) <- c("mfrm_imported_fit", "mfrm_fit", "list")
  out
}

.erm_normalize_schema_name <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(as.character(x)))
}

.erm_extract_person_table <- function(person_parameters) {
  theta_table <- person_parameters$theta.table %||% NULL
  if (!is.data.frame(theta_table) || nrow(theta_table) == 0L) {
    stop(
      paste0(
        "`eRm::person.parameter()` returned an unsupported person-parameter ",
        "schema: `theta.table` must be a non-empty data frame."
      ),
      call. = FALSE
    )
  }

  normalized_names <- .erm_normalize_schema_name(names(theta_table))
  estimate_columns <- which(normalized_names %in% c(
    "theta", "thetapar", "personparameter"
  ))
  if (length(estimate_columns) != 1L) {
    available <- paste(names(theta_table), collapse = ", ")
    if (!nzchar(available)) available <- "<none>"
    stop(
      paste0(
        "`eRm::person.parameter()` returned an unsupported `theta.table` ",
        "schema. Expected exactly one person-estimate column named `theta`, ",
        "`thetapar`, or `Person Parameter`; available columns: ", available,
        "."
      ),
      call. = FALSE
    )
  }

  estimate_raw <- theta_table[[estimate_columns]]
  estimates <- suppressWarnings(as.numeric(estimate_raw))
  if (length(estimates) != nrow(theta_table)) {
    stop(
      "The eRm person-estimate column could not be aligned to `theta.table` rows.",
      call. = FALSE
    )
  }

  id_columns <- which(normalized_names %in% c("person", "personid", "id"))
  if (length(id_columns) > 1L) {
    stop(
      "The eRm `theta.table` contains multiple possible person-ID columns.",
      call. = FALSE
    )
  }
  if (length(id_columns) == 1L) {
    person_ids <- as.character(theta_table[[id_columns]])
  } else {
    person_ids <- rownames(theta_table)
    default_rows <- identical(person_ids, as.character(seq_len(nrow(theta_table))))
    if (is.null(person_ids) || length(person_ids) != nrow(theta_table) ||
        default_rows || anyNA(person_ids) || any(!nzchar(person_ids))) {
      estimate_names <- names(estimate_raw)
      if (!is.null(estimate_names) && length(estimate_names) == nrow(theta_table) &&
          !anyNA(estimate_names) && all(nzchar(estimate_names))) {
        person_ids <- estimate_names
      } else {
        person_ids <- paste0("P", seq_len(nrow(theta_table)))
      }
    }
  }

  se_columns <- which(normalized_names %in% c(
    "se", "setheta", "stderr", "stderror", "standarderror"
  ))
  if (length(se_columns) > 1L) {
    stop(
      "The eRm `theta.table` contains multiple possible standard-error columns.",
      call. = FALSE
    )
  }
  if (length(se_columns) == 1L) {
    standard_errors <- suppressWarnings(as.numeric(theta_table[[se_columns]]))
  } else {
    standard_errors <- unlist(
      person_parameters$se.theta %||% numeric(0),
      use.names = TRUE
    )
    standard_errors <- suppressWarnings(as.numeric(standard_errors))
  }
  if (length(standard_errors) == 0L) {
    standard_errors <- rep(NA_real_, length(estimates))
  } else if (length(standard_errors) != length(estimates)) {
    stop(
      paste0(
        "`eRm::person.parameter()` returned ", length(estimates),
        " person estimates but ", length(standard_errors),
        " standard errors; the rows cannot be aligned safely."
      ),
      call. = FALSE
    )
  }

  data.frame(
    Person = person_ids,
    Estimate = estimates,
    SE = standard_errors,
    Extreme = ifelse(as.logical(theta_table$Interpolated %||% FALSE),
                     "source_extrapolated", "none"),
    Interpolated = as.logical(theta_table$Interpolated %||% FALSE),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Synthetic diagnostics bundle for imported fits
# ==============================================================================
#
# `.synthesize_imported_diagnostics()` builds an `mfrm_diagnostics`-shape
# object from the measure / fit columns extracted by
# `import_*_fit(compute_fit = TRUE)`. The bundle is intentionally
# minimal -- it carries the slots that downstream plot helpers
# (Wright map, QC dashboard, summary methods) actually look up,
# without re-running the residual / interaction / marginal-fit
# layers that would require the original observation table.
.synthesize_imported_diagnostics <- function(facet_others,
                                              persons,
                                              facet_names,
                                              n_obs,
                                              source = "imported") {
  facet_others <- as.data.frame(facet_others, stringsAsFactors = FALSE)
  persons <- as.data.frame(persons, stringsAsFactors = FALSE)

  # `measures`: the canonical Person + non-person facet measure
  # table that downstream helpers consume.
  measures_facets <- if (nrow(facet_others) > 0L) {
    df <- data.frame(
      Facet = as.character(facet_others$Facet),
      Level = as.character(facet_others$Level),
      Estimate = suppressWarnings(as.numeric(facet_others$Estimate)),
      SE = suppressWarnings(as.numeric(facet_others$SE %||% NA_real_)),
      ModelSE = suppressWarnings(as.numeric(facet_others$SE %||% NA_real_)),
      stringsAsFactors = FALSE
    )
    if ("Infit" %in% names(facet_others)) df$Infit <- facet_others$Infit
    if ("Outfit" %in% names(facet_others)) df$Outfit <- facet_others$Outfit
    if ("InfitZSTD" %in% names(facet_others)) df$InfitZSTD <- facet_others$InfitZSTD
    if ("OutfitZSTD" %in% names(facet_others)) df$OutfitZSTD <- facet_others$OutfitZSTD
    if ("FitBasis" %in% names(facet_others)) df$FitBasis <- facet_others$FitBasis
    df
  } else data.frame()
  measures_persons <- if (nrow(persons) > 0L) {
    df <- data.frame(
      Facet = "Person",
      Level = as.character(persons$Person),
      Estimate = suppressWarnings(as.numeric(persons$Estimate)),
      SE = suppressWarnings(as.numeric(persons$SE %||% NA_real_)),
      ModelSE = suppressWarnings(as.numeric(persons$SE %||% NA_real_)),
      stringsAsFactors = FALSE
    )
    if ("Infit" %in% names(persons)) df$Infit <- persons$Infit
    if ("Outfit" %in% names(persons)) df$Outfit <- persons$Outfit
    if ("InfitZSTD" %in% names(persons)) df$InfitZSTD <- persons$InfitZSTD
    if ("OutfitZSTD" %in% names(persons)) df$OutfitZSTD <- persons$OutfitZSTD
    if ("EstimateBasis" %in% names(persons)) df$EstimateBasis <- persons$EstimateBasis
    if ("UncertaintyBasis" %in% names(persons)) df$UncertaintyBasis <- persons$UncertaintyBasis
    if ("FitBasis" %in% names(persons)) df$FitBasis <- persons$FitBasis
    df
  } else data.frame()
  measures <- if (nrow(measures_facets) == 0L) measures_persons else
    if (nrow(measures_persons) == 0L) measures_facets else
      dplyr::bind_rows(measures_persons, measures_facets)

  fit_tbl <- if (nrow(measures) > 0L &&
                 all(c("Infit", "Outfit") %in% names(measures))) {
    measures[, c("Facet", "Level", "Infit", "Outfit",
                  intersect(c("InfitZSTD", "OutfitZSTD", "FitBasis"), names(measures))),
             drop = FALSE]
  } else data.frame()

  reliability <- .synthesize_reliability(measures)
  facets_chisq <- .synthesize_chisq(measures)

  out <- list(
    measures = measures,
    fit = fit_tbl,
    reliability = reliability,
    facets_chisq = facets_chisq,
    overall_fit = data.frame(
      Source = source,
      Method = "imported",
      Infit = NA_real_,
      Outfit = NA_real_,
      stringsAsFactors = FALSE
    ),
    obs = data.frame(),  # Original observation table is not available
    interrater = list(summary = data.frame(), pairs = data.frame()),
    diagnostic_basis = data.frame(
      DiagnosticPath = "imported",
      Status = "synthesised",
      Basis = paste0("Measurement-side bundle imported from `",
                      source, "`. The observation, residual, and ",
                      "interaction layers are not available; downstream ",
                      "helpers that depend on them should be re-run on ",
                      "an mfrmr-native fit."),
      stringsAsFactors = FALSE
    ),
    precision_profile = data.frame(
      Method = "imported",
      PrecisionTier = "imported",
      SupportsFormalInference = FALSE,
      RecommendedUse = "Imported estimates and SEs retain their source-package interpretation; no native uncertainty or joint facet test is established by import.",
      stringsAsFactors = FALSE
    ),
    precision_review = data.frame(),
    facet_names = facet_names,
    diagnostic_mode = "legacy",
    residual_pca_mode = "none",
    n_obs = as.integer(n_obs),
    metric_version = 1L,
    imported = TRUE,
    source = source
  )
  class(out) <- c("mfrm_imported_diagnostics", "mfrm_diagnostics", "list")
  out
}

.synthesize_reliability <- function(measures) {
  if (!is.data.frame(measures) || nrow(measures) == 0L ||
      !all(c("Facet", "Estimate", "SE") %in% names(measures))) {
    return(data.frame())
  }
  facets <- split(measures, measures$Facet)
  rows <- lapply(names(facets), function(fct) {
    stats <- summarize_precision_basis(facets[[fct]], "SE", "sample")
    posterior_sd <- "UncertaintyBasis" %in% names(facets[[fct]]) &&
      any(grepl("Conditional posterior SD", facets[[fct]]$UncertaintyBasis, fixed = TRUE))
    if (posterior_sd) {
      stats$Separation <- stats$Strata <- stats$Reliability <- NA_real_
      stats$SummaryNote <- "Posterior SD is conditional scoring uncertainty, not a sampling SE for separation reliability."
    }
    data.frame(
      Facet = fct, Levels = nrow(facets[[fct]]),
      EstimateAvailable = stats$EstimateAvailable, SEAvailable = stats$SEAvailable,
      ExcludedEstimates = stats$ExcludedEstimates, SummaryNote = stats$SummaryNote,
      Separation = stats$Separation, Strata = stats$Strata, Reliability = stats$Reliability,
      SupportsFormalInference = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.synthesize_chisq <- function(measures) {
  if (!is.data.frame(measures) || nrow(measures) == 0L ||
      !all(c("Facet", "Estimate", "SE") %in% names(measures))) {
    return(data.frame())
  }
  facets <- split(measures, measures$Facet)
  rows <- lapply(names(facets), function(fct) {
    rows_f <- facets[[fct]]
    est <- suppressWarnings(as.numeric(rows_f$Estimate))
    est <- est[is.finite(est)]
    data.frame(
      Facet = fct,
      Levels = nrow(rows_f),
      MeanMeasure = if (length(est)) mean(est) else NA_real_,
      SD = if (length(est) > 1L) stats::sd(est) else NA_real_,
      FixedChiSq = NA_real_,
      FixedDF = NA_real_,
      FixedProb = NA_real_,
      SupportsFormalInference = FALSE,
      TestNote = "A joint facet test is unavailable: imported marginal SEs do not supply the required covariance and estimation assumptions.",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @export
print.mfrm_imported_fit <- function(x, ...) {
  .validate_imported_metric(x)
  cat("Imported measurement results (limited functionality)\n")
  cat(sprintf("  Source: %s (%s)\n",
              x$source$package %||% "unknown",
              x$source$source_object_class %||% "unknown"))
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    ov <- x$summary[1, , drop = FALSE]
    cat(sprintf("  Model: %s | Method: %s | Persons: %s\n",
                ov$Model, ov$Method, ov$Persons))
  }
  print_wrapped_line(x$source$metric)
  print_wrapped_line(x$source$person_scoring)
  if (!is.null(x$source$person_identification)) print_wrapped_line(x$source$person_identification)
  cat("  Use summary(x) for source-scale tables or plot(x) for a point-only Wright map.\n")
  cat("  SEs retain source-package conventions; joint facet tests and native QC are unavailable.\n")
  cat("  Bias, DIF, anchoring, and reproducibility workflows require a native fit_mfrm() object.\n")
  invisible(x)
}

.validate_imported_metric <- function(x) {
  if (!identical(x$source$metric_version, 1L)) {
    stop("Re-import the existing source-package fit to update item difficulties, thresholds and uncertainty labels; no re-estimation is needed.", call. = FALSE)
  }
  invisible(x)
}

#' @rdname import_mirt_fit
#' @param object,x An imported measurement bundle.
#' @param digits Number of digits for displayed estimates.
#' @param ... Additional arguments (unused by imported summaries).
#' @export
summary.mfrm_imported_fit <- function(object, digits = 3L, ...) {
  .validate_imported_metric(object)
  digits <- mfrmr_calibration_score_digits(digits)
  display <- function(table) {
    columns <- intersect(c("Estimate", "SE", "Slope"), names(table))
    table[columns] <- lapply(table[columns], round, digits = digits)
    table
  }
  structure(list(
    overview = object$summary, facets = lapply(object$facets, display),
    steps = display(object$steps), source = object$source,
    scale_contract = mfrm_fit_scale_contract(object)
  ), class = c("summary.mfrm_imported_fit", "list"))
}

#' @rdname import_mirt_fit
#' @export
print.summary.mfrm_imported_fit <- function(x, ...) {
  cat("Imported measurement summary\n")
  cat("  Source: ", x$source$package, "; Persons: ", nrow(x$facets$person),
      "; item/design-cell measures: ", nrow(x$facets$others), "\n", sep = "")
  print_wrapped_line(x$source$metric)
  print_wrapped_line(x$source$person_scoring)
  if (!is.null(x$source$person_identification)) print_wrapped_line(x$source$person_identification)
  cat("\nItem/design-cell measures (first 10)\n")
  print(utils::head(x$facets$others, 10L), row.names = FALSE)
  print_wrapped_line("SEs retain source conventions; unavailable transformed SEs remain missing. Import does not establish joint facet tests or the model's suitability for inference.")
  cat("Person and absolute threshold tables are retained in $facets$person and $steps.\n")
  invisible(x)
}
