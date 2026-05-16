# ==============================================================================
# Import adapters: mirt / TAM / eRm fits to mfrmr-compatible bundles
# ==============================================================================
#
# The current import path is **partial**: each helper extracts the
# subset of measurement information that maps cleanly to the
# `mfrm_fit` contract (item / step parameters, person scores, basic
# fit statistics) and returns a thin object that the mfrmr plot and
# reporting helpers can consume. Full bundle import (bias / DIF /
# anchor review / replay) is deferred to a future release because the
# source packages do not always expose the underlying data.
#
# All importers refuse to claim the `mfrm_fit` class outright; they
# return an `mfrm_imported_fit` object that downstream helpers can
# detect and handle conservatively.

#' Import an `mirt` fit to an mfrmr-compatible bundle
#'
#' Extracts item, step, and person parameters from a [mirt::mirt()]
#' fit and returns an `mfrm_imported_fit` object. The returned
#' object has the public slots `summary`, `facets$person`,
#' `facets$others`, `steps`, `config`, and `source` that the mfrmr
#' plot and table helpers expect. With `compute_fit = TRUE` the
#' importer also runs [mirt::itemfit()] and [mirt::personfit()] so
#' Infit / Outfit columns are populated, and synthesises a
#' `mfrm_diagnostics`-shape `diagnostics` slot consumable by
#' downstream plot helpers (Wright map, QC dashboard, etc.).
#'
#' @param fit An object returned by [mirt::mirt()] (a
#'   `SingleGroupClass`).
#' @param model One of `"RSM"`, `"PCM"`, `"GPCM"`. The importer
#'   does not infer the model from the mirt object; pass the model
#'   that was estimated.
#' @param item_facet Name to assign to the item facet in the
#'   imported bundle (default `"Item"`).
#' @param compute_fit Logical. When `TRUE`, run [mirt::itemfit()]
#'   and [mirt::personfit()] to populate Infit / Outfit / OutfitZSTD
#'   columns on the returned facet tables, plus build a
#'   measurement-side `mfrm_diagnostics` bundle consumable by
#'   `summary()`, `plot.mfrm_fit()`, `plot_qc_dashboard()`, etc.
#'   Default `FALSE` keeps the importer fast (skeleton only).
#'
#' @return An `mfrm_imported_fit` object. Slots:
#' \describe{
#'   \item{`summary`}{Model / method / N / LogLik / AIC / BIC.}
#'   \item{`facets$person`}{Person ID, Estimate, SE, Extreme, plus
#'     Infit / Outfit / OutfitZSTD / Zh when `compute_fit = TRUE`.}
#'   \item{`facets$others`}{Item-level estimates and slopes; with
#'     `compute_fit = TRUE`, also Infit / Outfit / S_X2 / RMSEA / df
#'     from `mirt::itemfit()`.}
#'   \item{`steps`}{Per-item threshold parameters extracted from the
#'     IRT parameterisation (`b1`, ..., `b(K-1)`).}
#'   \item{`config`}{List with the resolved `model` and `item_facet`
#'     used for the import; downstream plot and table helpers consult
#'     this to dispatch correctly on the imported bundle.}
#'   \item{`diagnostics`}{`mfrm_diagnostics`-shape bundle when
#'     `compute_fit = TRUE`; `NULL` otherwise.}
#'   \item{`source`}{Imported-from metadata.}
#' }
#'
#' @section Scope:
#' Bundles bias / DIF / anchor / replay slots are explicitly not
#' populated; full bidirectional import / export is planned for a
#' future release.
#' @seealso [import_tam_fit()], [import_erm_fit()]
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

  items <- .mirt_extract_items(fit)
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
    fit_attached <- .mirt_compute_fit_stats(fit)
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
      n_obs = nrow(persons) * nrow(facet_others),
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
                  source_object_class = class(fit)[1],
                  compute_fit = isTRUE(compute_fit))
  )
  class(out) <- c("mfrm_imported_fit", "mfrm_fit", "list")
  out
}

# --- mirt internal extractors --------------------------------------------

.mirt_extract_items <- function(fit) {
  irt_pars <- tryCatch(
    mirt::coef(fit, simplify = TRUE, IRTpars = TRUE)$items,
    error = function(e) NULL
  )
  if (is.null(irt_pars)) {
    irt_pars <- mirt::coef(fit, simplify = TRUE, IRTpars = FALSE)$items
  }
  items_mat <- as.data.frame(irt_pars, stringsAsFactors = FALSE)
  items_mat$Level <- rownames(items_mat)
  rownames(items_mat) <- NULL
  # Slope (1PL / Rasch returns no slope column; populate as 1).
  slope_col <- intersect(c("a", "a1"), names(items_mat))[1]
  items_mat$Slope <- if (!is.na(slope_col)) {
    suppressWarnings(as.numeric(items_mat[[slope_col]]))
  } else {
    rep(1, nrow(items_mat))
  }
  # Item difficulty: the average of the b-thresholds for graded /
  # gpcm models, or the single `b` for binary Rasch.
  b_cols <- grep("^b[0-9]*$", names(items_mat), value = TRUE)
  if (length(b_cols) == 0L) {
    items_mat$Difficulty <- NA_real_
  } else {
    b_mat <- vapply(b_cols, function(col) {
      suppressWarnings(as.numeric(items_mat[[col]]))
    }, numeric(nrow(items_mat)))
    if (is.matrix(b_mat) && ncol(b_mat) > 1L) {
      items_mat$Difficulty <- rowMeans(b_mat, na.rm = TRUE)
    } else {
      items_mat$Difficulty <- as.numeric(b_mat)
    }
  }
  items_mat
}

.mirt_extract_steps <- function(items_mat, item_facet) {
  b_cols <- grep("^b[0-9]+$", names(items_mat), value = TRUE)
  b_cols <- b_cols[order(as.integer(sub("^b", "", b_cols)))]
  if (length(b_cols) == 0L) {
    # Single-threshold Rasch: emit one step per item.
    if ("b" %in% names(items_mat)) {
      return(data.frame(
        StepFacet = item_facet,
        Level = items_mat$Level,
        Step = 1L,
        Estimate = suppressWarnings(as.numeric(items_mat$b)),
        stringsAsFactors = FALSE
      ))
    }
    return(data.frame(
      StepFacet = character(0), Level = character(0),
      Step = integer(0), Estimate = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(seq_along(b_cols), function(k) {
    data.frame(
      StepFacet = item_facet,
      Level = items_mat$Level,
      Step = k,
      Estimate = suppressWarnings(as.numeric(items_mat[[b_cols[k]]])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.mirt_extract_persons <- function(fit) {
  fscores <- mirt::fscores(fit, full.scores.SE = TRUE, verbose = FALSE)
  person_tbl <- data.frame(
    Person = paste0("P", seq_len(nrow(fscores))),
    Estimate = as.numeric(fscores[, 1]),
    SE = if (ncol(fscores) >= 2L) as.numeric(fscores[, 2]) else NA_real_,
    Extreme = "none",
    stringsAsFactors = FALSE
  )
  person_tbl
}

.mirt_compute_fit_stats <- function(fit) {
  item_fit <- tryCatch(
    suppressMessages(suppressWarnings(
      mirt::itemfit(fit, fit_stats = c("infit"),
                     na.rm = TRUE)
    )),
    error = function(e) NULL
  )
  person_fit <- tryCatch(
    suppressMessages(suppressWarnings(
      mirt::personfit(fit, method = "MAP")
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
    df$Person <- paste0("P", seq_len(nrow(df)))
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
#' Extracts item / step / person parameters from a [TAM::tam.mml()],
#' [TAM::tam.jml()], or [TAM::tam.mml.mfr()] fit. The multi-facet
#' `tam.mml.mfr()` path is detected automatically and each
#' non-person facet is mapped onto a row of `fit$facets$others`
#' so downstream MFRM helpers (e.g. `plot_qc_dashboard()`) work
#' on the imported object.
#'
#' @param fit An object returned by `TAM::tam.mml()`,
#'   `TAM::tam.jml()`, or `TAM::tam.mml.mfr()`.
#' @param model Same as [import_mirt_fit()].
#' @param item_facet Name to assign to the item facet for the
#'   single-facet path. Ignored when the input is a multi-facet
#'   `tam.mml.mfr` fit (the original facet names are preserved).
#' @param compute_fit Logical. When `TRUE`, run [TAM::tam.fit()]
#'   and [TAM::tam.personfit()] to populate Infit / Outfit columns
#'   on the returned facet tables, plus build a measurement-side
#'   `mfrm_diagnostics` bundle. Default `FALSE`.
#'
#' @return An `mfrm_imported_fit` object. Slots mirror
#'   [import_mirt_fit()].
#' @seealso [import_mirt_fit()], [import_erm_fit()]
#' @export
import_tam_fit <- function(fit, model = c("RSM", "PCM", "GPCM"),
                            item_facet = "Item",
                            compute_fit = FALSE) {
  if (!requireNamespace("TAM", quietly = TRUE)) {
    stop("`import_tam_fit()` requires the `TAM` package (suggested).",
         call. = FALSE)
  }
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
    extracted <- .tam_extract_mfr(fit, persons = persons,
                                    fallback_facet = item_facet)
    facet_others <- extracted$others
    steps <- extracted$steps
    facet_names <- extracted$facet_names
  } else {
    extracted <- .tam_extract_single(fit, item_facet = item_facet)
    facet_others <- extracted$others
    steps <- extracted$steps
    facet_names <- item_facet
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

  summary_tbl <- data.frame(
    Model = model,
    Method = "MML",
    Source = "TAM",
    N = nrow(persons),
    Persons = nrow(persons),
    Facets = length(unique(facet_others$Facet)),
    Categories = NA_integer_,
    LogLik = as.numeric(fit$ic$loglike %||% NA_real_),
    AIC = as.numeric(fit$ic$AIC %||% NA_real_),
    BIC = as.numeric(fit$ic$BIC %||% NA_real_),
    Converged = TRUE,
    ConvergenceStatus = "imported",
    stringsAsFactors = FALSE
  )

  diagnostics <- NULL
  if (isTRUE(compute_fit)) {
    diagnostics <- .synthesize_imported_diagnostics(
      facet_others = facet_others,
      persons = persons,
      facet_names = facet_names,
      n_obs = nrow(persons) * nrow(facet_others),
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
                  multi_facet = is_mfr),
    source = list(package = "TAM",
                  source_object_class = class(fit)[1],
                  multi_facet = is_mfr,
                  compute_fit = isTRUE(compute_fit))
  )
  class(out) <- c("mfrm_imported_fit", "mfrm_fit", "list")
  out
}

# --- TAM internal extractors ---------------------------------------------

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
    Extreme = "none",
    stringsAsFactors = FALSE
  )
}

.tam_extract_single <- function(fit, item_facet) {
  xsi <- as.data.frame(fit$xsi, stringsAsFactors = FALSE)
  xsi$Level <- rownames(xsi)
  rownames(xsi) <- NULL
  facet_others <- data.frame(
    Facet = item_facet,
    Level = xsi$Level,
    Estimate = suppressWarnings(as.numeric(xsi$xsi)),
    SE = suppressWarnings(as.numeric(xsi$se.xsi %||% NA_real_)),
    stringsAsFactors = FALSE
  )
  # Step parameters: TAM stores per-item-step thresholds in $A and
  # $xsi together. We approximate by detecting `_Cat` suffix in the
  # parameter name (TAM convention) and parsing the step index.
  steps <- .tam_steps_from_xsi(xsi, step_facet = item_facet)
  list(others = facet_others, steps = steps,
       facet_names = item_facet)
}

.tam_extract_mfr <- function(fit, persons, fallback_facet = "Item") {
  # `tam.mml.mfr()` exposes a tidy `xsi.facets` table that tells us
  # which facet each parameter belongs to. Each row has columns
  # `parameter` (xsi name like "I1" or "raterR1"), `facet`
  # ("item" / "rater" / ...), `xsi`, `se.xsi`.
  xf <- as.data.frame(fit$xsi.facets, stringsAsFactors = FALSE)
  if (nrow(xf) == 0L) {
    return(list(others = data.frame(Facet = character(0),
                                     Level = character(0),
                                     Estimate = numeric(0),
                                     SE = numeric(0),
                                     stringsAsFactors = FALSE),
                steps = data.frame(StepFacet = character(0),
                                    Level = character(0),
                                    Step = integer(0),
                                    Estimate = numeric(0),
                                    stringsAsFactors = FALSE),
                facet_names = fallback_facet))
  }
  xf$facet <- as.character(xf$facet)
  xf$parameter <- as.character(xf$parameter)
  xf$xsi <- suppressWarnings(as.numeric(xf$xsi))
  xf$se <- suppressWarnings(as.numeric(xf$se.xsi %||% NA_real_))

  # Polytomous step parameters carry a `_Cat<k>` suffix on the
  # parameter name. Split into main effects vs steps so the
  # `mfrm_fit` shape (others + steps) is preserved. We use
  # `grepl()` for the boolean mask (`regexpr()` + `regmatches()`
  # would silently drop non-matching positions, leaving us with a
  # zero-length boolean for binary models).
  is_step <- grepl("_Cat[0-9]+$", xf$parameter)
  step_idx <- ifelse(
    is_step,
    suppressWarnings(as.integer(sub(".*_Cat([0-9]+)$", "\\1", xf$parameter))),
    NA_integer_
  )

  # Strip the facet prefix from the parameter name to recover the
  # human-readable level. TAM names parameters as
  # `<facet><level>` (no separator) for facets other than `item`,
  # and `<level>` for `item`.
  level_label <- xf$parameter
  for (idx in seq_len(nrow(xf))) {
    fct <- xf$facet[idx]
    if (!identical(fct, "item") && nzchar(fct) &&
        startsWith(level_label[idx], fct)) {
      level_label[idx] <- substr(level_label[idx],
                                  nchar(fct) + 1L,
                                  nchar(level_label[idx]))
    }
  }
  # Strip the `_Cat<k>` suffix from the step parameter labels so
  # they share the level name with their main-effect partner.
  level_no_step <- level_label
  level_no_step[is_step] <- sub("_Cat[0-9]+$", "",
                                  level_label[is_step])

  facet_others <- data.frame(
    Facet = xf$facet[!is_step],
    Level = level_no_step[!is_step],
    Estimate = xf$xsi[!is_step],
    SE = xf$se[!is_step],
    stringsAsFactors = FALSE
  )
  facet_others <- unique(facet_others)
  # Capitalise the facet label so it matches the package's typical
  # `Rater` / `Item` / `Criterion` casing.
  facet_others$Facet <- vapply(facet_others$Facet, function(s) {
    if (nchar(s) == 0L) return(s)
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }, character(1))

  steps_df <- if (any(is_step)) {
    df <- data.frame(
      StepFacet = xf$facet[is_step],
      Level = level_no_step[is_step],
      Step = step_idx[is_step],
      Estimate = xf$xsi[is_step],
      stringsAsFactors = FALSE
    )
    df$StepFacet <- vapply(df$StepFacet, function(s) {
      if (nchar(s) == 0L) return(s)
      paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
    }, character(1))
    df
  } else {
    data.frame(StepFacet = character(0), Level = character(0),
               Step = integer(0), Estimate = numeric(0),
               stringsAsFactors = FALSE)
  }

  facet_names <- unique(facet_others$Facet)
  if (length(facet_names) == 0L) facet_names <- fallback_facet
  list(others = facet_others, steps = steps_df,
       facet_names = facet_names)
}

.tam_steps_from_xsi <- function(xsi, step_facet) {
  is_step <- grepl("_Cat[0-9]+$", xsi$Level)
  if (!any(is_step)) {
    return(data.frame(StepFacet = character(0), Level = character(0),
                      Step = integer(0), Estimate = numeric(0),
                      stringsAsFactors = FALSE))
  }
  level <- sub("_Cat[0-9]+$", "", xsi$Level[is_step])
  step_idx <- suppressWarnings(as.integer(
    sub(".*_Cat([0-9]+)$", "\\1", xsi$Level[is_step])
  ))
  data.frame(
    StepFacet = step_facet,
    Level = level,
    Step = step_idx,
    Estimate = suppressWarnings(as.numeric(xsi$xsi[is_step])),
    stringsAsFactors = FALSE
  )
}

.tam_compute_fit_stats <- function(fit, persons) {
  item_fit <- tryCatch(
    suppressMessages(suppressWarnings(TAM::tam.fit(fit))),
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
    df
  } else NULL
  list(item = item_df, person = person_df)
}

#' Import an `eRm` fit to an mfrmr-compatible bundle
#'
#' Extracts item / person parameters from an [eRm::PCM()] /
#' [eRm::RM()] fit. Same caveats as [import_mirt_fit()].
#'
#' @param fit An object returned by `eRm::PCM()`, `eRm::RM()`, or
#'   `eRm::RSM()`.
#' @param model Same as [import_mirt_fit()].
#' @param item_facet Name to assign to the item facet.
#'
#' @return An `mfrm_imported_fit` object.
#' @seealso [import_mirt_fit()], [import_tam_fit()]
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
  beta <- fit$betapar
  facet_others <- data.frame(
    Facet = item_facet,
    Level = names(beta),
    Estimate = as.numeric(beta),
    SE = as.numeric(fit$se.beta %||% NA_real_),
    stringsAsFactors = FALSE
  )
  pp <- tryCatch(eRm::person.parameter(fit), error = function(e) NULL)
  if (!is.null(pp)) {
    theta <- as.numeric(pp$theta.table$theta)
    person_tbl <- data.frame(
      Person = paste0("P", seq_along(theta)),
      Estimate = theta,
      SE = as.numeric(pp$se.theta[[1]] %||% NA_real_),
      Extreme = "none",
      stringsAsFactors = FALSE
    )
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
    Converged = TRUE,
    ConvergenceStatus = "imported",
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_tbl,
    facets = list(person = person_tbl, others = facet_others),
    steps = data.frame(StepFacet = character(0), Step = integer(0),
                       Level = character(0), Estimate = numeric(0),
                       stringsAsFactors = FALSE),
    config = list(model = model, method = "CML",
                  facet_names = item_facet,
                  source = "eRm"),
    source = list(package = "eRm",
                  source_object_class = class(fit)[1])
  )
  class(out) <- c("mfrm_imported_fit", "mfrm_fit", "list")
  out
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
    df
  } else data.frame()
  measures <- if (nrow(measures_facets) == 0L) measures_persons else
    if (nrow(measures_persons) == 0L) measures_facets else
      dplyr::bind_rows(measures_persons, measures_facets)

  fit_tbl <- if (nrow(measures) > 0L &&
                 all(c("Infit", "Outfit") %in% names(measures))) {
    measures[, c("Facet", "Level", "Infit", "Outfit",
                  intersect(c("InfitZSTD", "OutfitZSTD"), names(measures))),
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
      stringsAsFactors = FALSE
    ),
    precision_review = data.frame(),
    facet_names = facet_names,
    diagnostic_mode = "legacy",
    residual_pca_mode = "none",
    n_obs = as.integer(n_obs),
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
    rows_f <- facets[[fct]]
    est <- suppressWarnings(as.numeric(rows_f$Estimate))
    se <- suppressWarnings(as.numeric(rows_f$SE))
    if (sum(is.finite(est)) < 2L) {
      return(data.frame(
        Facet = fct, Levels = sum(is.finite(est)),
        Separation = NA_real_, Strata = NA_real_, Reliability = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    obs_var <- stats::var(est, na.rm = TRUE)
    err_var <- if (any(is.finite(se))) {
      mean(se[is.finite(se)]^2, na.rm = TRUE)
    } else NA_real_
    true_var <- if (is.finite(err_var)) max(obs_var - err_var, 0) else NA_real_
    rmse <- if (is.finite(err_var)) sqrt(err_var) else NA_real_
    sep <- if (is.finite(rmse) && rmse > 0 && is.finite(true_var)) {
      sqrt(true_var) / rmse
    } else NA_real_
    rel <- if (is.finite(obs_var) && obs_var > 0 && is.finite(true_var)) {
      true_var / obs_var
    } else NA_real_
    strata <- if (is.finite(sep)) (4 * sep + 1) / 3 else NA_real_
    data.frame(
      Facet = fct, Levels = nrow(rows_f),
      Separation = sep, Strata = strata, Reliability = rel,
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
    se <- suppressWarnings(as.numeric(rows_f$SE))
    n_lev <- length(est)
    df <- max(0L, n_lev - 1L)
    chi <- if (df > 0L && any(is.finite(se) & se > 0)) {
      mean_est <- mean(est, na.rm = TRUE)
      sum(((est - mean_est) / se)^2, na.rm = TRUE)
    } else NA_real_
    pval <- if (is.finite(chi) && df > 0) {
      stats::pchisq(chi, df, lower.tail = FALSE)
    } else NA_real_
    data.frame(
      Facet = fct,
      Levels = n_lev,
      MeanMeasure = mean(est, na.rm = TRUE),
      SD = stats::sd(est, na.rm = TRUE),
      FixedChiSq = chi,
      FixedDF = df,
      FixedProb = pval,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @export
print.mfrm_imported_fit <- function(x, ...) {
  cat("mfrmr imported fit (measurement-side slots only)\n")
  cat(sprintf("  Source: %s (%s)\n",
              x$source$package %||% "unknown",
              x$source$source_object_class %||% "unknown"))
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    ov <- x$summary[1, , drop = FALSE]
    cat(sprintf("  Model: %s | Method: %s | Persons: %s\n",
                ov$Model, ov$Method, ov$Persons))
  }
  cat("  Use mfrmr plot helpers for measurement-side views;\n")
  cat("  bias / DIF / anchor / replay slots are not available.\n")
  invisible(x)
}
