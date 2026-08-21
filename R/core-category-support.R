# ==============================================================================
# Category and step-support audit
# ==============================================================================
#
# The current package fits one observed score scale at a time.  This audit uses
# an explicit internal scale key so that the one-scale reduction does not bake
# an anonymous global assumption into the later multiple-scale design.
#
# Under the current within-ladder sum-zero parameterization, a missing internal
# category creates an exact recession direction: increase the step immediately
# below that category and decrease the step immediately above it.  The missing
# category probability then tends to zero while every other category exponent
# is unchanged.  Boundary-category absence does not create this step-contrast
# direction and remains review evidence for the later element-boundary audit.
# Threshold anchors are not part of the 0.2.3 model.

mfrmr_category_contract_version <- function() {
  mfrmr_readiness_contract_version()
}

mfrmr_category_values_text <- function(x) {
  x <- as.character(x)
  if (length(x) == 0L) "" else paste(x, collapse = ";")
}

mfrmr_category_normalized_entropy <- function(counts) {
  counts <- as.numeric(counts)
  total <- sum(counts)
  if (!is.finite(total) || total <= 0 || length(counts) <= 1L) {
    return(NA_real_)
  }
  probabilities <- counts[counts > 0] / total
  -sum(probabilities * log(probabilities)) / log(length(counts))
}

mfrmr_category_counts <- function(score, declared) {
  as.integer(vapply(
    declared,
    function(category) sum(score == category, na.rm = TRUE),
    integer(1)
  ))
}

mfrmr_category_weighted_counts <- function(score, weight, declared) {
  as.numeric(vapply(
    declared,
    function(category) sum(weight[score == category], na.rm = TRUE),
    numeric(1)
  ))
}

mfrmr_category_zero_type <- function(category_index, n_categories,
                                      global_count, global_weight,
                                      within_scope_count, within_scope_weight,
                                      shared_scope) {
  if (is.finite(within_scope_weight) && within_scope_weight > 0) return("none")
  position <- if (category_index %in% c(1L, n_categories)) {
    "boundary"
  } else {
    "internal"
  }
  zero_class <- if (global_count <= 0) {
    "sampling_zero_global"
  } else if (global_weight <= 0) {
    "nonpositive_weight_global"
  } else if (within_scope_count <= 0) {
    if (isTRUE(shared_scope)) "sampling_zero_global" else "sampling_zero_step_scope"
  } else {
    if (isTRUE(shared_scope)) "nonpositive_weight_global" else "nonpositive_weight_step_scope"
  }
  paste(zero_class, position, sep = "_")
}

mfrmr_category_scope_specs <- function(prep, config) {
  n <- nrow(prep$data)
  if (identical(config$model, "RSM")) {
    return(list(list(
      label = "shared_rating_scale",
      level = "shared",
      row_mask = rep(TRUE, n),
      shared = TRUE
    )))
  }

  step_facet <- as.character(config$step_facet %||% "")
  if (!nzchar(step_facet) || !step_facet %in% names(prep$data)) {
    stop("Category support audit requires a valid PCM/GPCM `step_facet`.",
         call. = FALSE)
  }
  levels <- as.character(config$facet_levels[[step_facet]] %||% character(0))
  lapply(levels, function(level) {
    list(
      label = paste0(step_facet, "=", level),
      level = level,
      row_mask = as.character(prep$data[[step_facet]]) == level,
      shared = FALSE
    )
  })
}

mfrmr_category_local_support <- function(prep, declared) {
  if (length(prep$facet_names %||% character(0)) == 0L) {
    return(data.frame())
  }
  rows <- list()
  cursor <- 0L
  for (facet in prep$facet_names) {
    levels <- as.character(prep$levels[[facet]] %||% character(0))
    for (level in levels) {
      cursor <- cursor + 1L
      mask <- as.character(prep$data[[facet]]) == level
      counts <- mfrmr_category_counts(prep$data$Score[mask], declared)
      weighted_counts <- mfrmr_category_weighted_counts(
        prep$data$Score[mask], prep$data$Weight[mask], declared
      )
      missing <- declared[weighted_counts <= 0]
      observed_counts <- counts[weighted_counts > 0]
      weak <- length(missing) > 0L || any(observed_counts <= 1L)
      rows[[cursor]] <- data.frame(
        ScaleScope = "single_observed_scale",
        Facet = as.character(facet),
        Level = as.character(level),
        Observations = sum(mask),
        ObservedCategories = mfrmr_category_values_text(declared[counts > 0L]),
        MissingCategories = mfrmr_category_values_text(missing),
        ObservedPositiveWeightCategories = mfrmr_category_values_text(
          declared[weighted_counts > 0]
        ),
        BoundaryCategoryMissing = any(missing %in% declared[c(1L, length(declared))]),
        MinimumObservedCategoryCount = if (length(observed_counts)) {
          min(observed_counts)
        } else {
          0L
        },
        MinimumPositiveWeightedN = if (any(weighted_counts > 0)) {
          min(weighted_counts[weighted_counts > 0])
        } else {
          0
        },
        MaximumCategoryFraction = if (sum(counts) > 0L) {
          max(counts) / sum(counts)
        } else {
          NA_real_
        },
        MaximumWeightedCategoryFraction = if (sum(weighted_counts) > 0) {
          max(weighted_counts) / sum(weighted_counts)
        } else {
          NA_real_
        },
        NormalizedCategoryEntropy = mfrmr_category_normalized_entropy(counts),
        WeightedNormalizedCategoryEntropy = mfrmr_category_normalized_entropy(
          weighted_counts
        ),
        InformationState = if (weak) "weak_information" else "adequate",
        ReasonCode = if (weak) "weak_category_information" else "",
        stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  rownames(out) <- NULL
  out
}

audit_mfrm_category_support <- function(prep, config, sizes) {
  declared <- seq.int(as.integer(prep$rating_min), as.integer(prep$rating_max))
  n_categories <- length(declared)
  n_steps <- max(n_categories - 1L, 0L)
  free_per_scope <- sum_zero_param_count(n_steps)
  if (n_categories < 2L || n_steps < 1L) {
    stop("Category support audit requires at least two declared categories.",
         call. = FALSE)
  }

  score <- as.integer(prep$data$Score)
  weight <- as.numeric(prep$data$Weight)
  global_counts <- mfrmr_category_counts(score, declared)
  global_weighted <- mfrmr_category_weighted_counts(score, weight, declared)
  global_observed <- declared[global_weighted > 0]
  global_missing <- declared[global_weighted <= 0]
  scope_specs <- mfrmr_category_scope_specs(prep, config)
  step_slice <- build_param_slices(sizes)$steps %||% integer(0)

  support_rows <- vector("list", length(scope_specs))
  category_rows <- vector("list", length(scope_specs))
  transition_rows <- vector("list", length(scope_specs))
  contrast_rows <- vector("list", length(scope_specs))

  for (scope_index in seq_along(scope_specs)) {
    scope <- scope_specs[[scope_index]]
    mask <- as.logical(scope$row_mask)
    scope_score <- score[mask]
    scope_weight <- weight[mask]
    counts <- mfrmr_category_counts(scope_score, declared)
    weighted_counts <- mfrmr_category_weighted_counts(
      scope_score, scope_weight, declared
    )
    observed <- declared[weighted_counts > 0]
    missing <- declared[weighted_counts <= 0]
    zero_types <- vapply(seq_along(declared), function(category_index) {
      mfrmr_category_zero_type(
        category_index = category_index,
        n_categories = n_categories,
        global_count = global_counts[category_index],
        global_weight = global_weighted[category_index],
        within_scope_count = counts[category_index],
        within_scope_weight = weighted_counts[category_index],
        shared_scope = isTRUE(scope$shared)
      )
    }, character(1))
    internal_indices <- if (n_categories > 2L) {
      seq.int(2L, n_categories - 1L)
    } else {
      integer(0)
    }
    unsupported_category_indices <- internal_indices[
      weighted_counts[internal_indices] <= 0
    ]
    if (free_per_scope == 0L) {
      unsupported_category_indices <- integer(0)
    }

    transition_scope_rows <- vector("list", n_steps)
    weak_transition <- logical(n_steps)
    for (transition in seq_len(n_steps)) {
      lower_count <- sum(counts[seq_len(transition)])
      upper_count <- sum(counts[(transition + 1L):n_categories])
      lower_weight <- sum(weighted_counts[seq_len(transition)])
      upper_weight <- sum(weighted_counts[(transition + 1L):n_categories])
      supported <- lower_weight > 0 && upper_weight > 0

      adjacent_indices <- c(transition, transition + 1L)
      adjacent_zero <- any(weighted_counts[adjacent_indices] <= 0)
      singleton_side <- min(lower_count, upper_count) <= 1L
      weak_transition[transition] <- supported &&
        (adjacent_zero || singleton_side)
      unsupported_contrast <- any(
        unsupported_category_indices %in% adjacent_indices
      )

      is_free_component <- transition <= free_per_scope
      optimizer_index <- if (is_free_component) {
        offset <- (scope_index - 1L) * free_per_scope + transition
        if (offset <= length(step_slice)) as.integer(step_slice[offset]) else NA_integer_
      } else {
        NA_integer_
      }
      coordinate_prefix <- if (isTRUE(scope$shared)) {
        "shared_rating_scale"
      } else {
        paste0(config$step_facet, ":", scope$level)
      }
      status <- if (free_per_scope == 0L) {
        "not_estimated"
      } else if (unsupported_contrast) {
        "unsupported"
      } else if (weak_transition[transition]) {
        "weak_information"
      } else {
        "estimable"
      }
      transition_reason <- character(0)
      if (!isTRUE(scope$shared) && length(missing) > 0L) {
        transition_reason <- c(
          transition_reason, "step_scope_category_unobserved"
        )
      }
      if (unsupported_contrast) {
        transition_reason <- c(
          transition_reason,
          "step_scope_category_unobserved",
          "unsupported_step_coordinate"
        )
      } else if (weak_transition[transition]) {
        transition_reason <- c(transition_reason, "weak_category_information")
      }
      transition_scope_rows[[transition]] <- data.frame(
        ReadinessScope = "parameter",
        ScaleScope = "single_observed_scale",
        StepScope = scope$label,
        ParameterClass = "step",
        Coordinate = paste0(coordinate_prefix, ":Step", transition),
        Step = paste0("Step", transition),
        Transition = transition,
        LowerCategory = declared[transition],
        UpperCategory = declared[transition + 1L],
        RetainedForFit = TRUE,
        Fixed = FALSE,
        ConstraintRole = if (free_per_scope == 0L) {
          "derived_constraint_only"
        } else if (is_free_component) {
          "free_coordinate_component"
        } else {
          "sum_zero_reference_component"
        },
        OptimizerIndex = optimizer_index,
        LowerCount = lower_count,
        UpperCount = upper_count,
        LowerWeightedN = lower_weight,
        UpperWeightedN = upper_weight,
        CrossingSupport = supported,
        FreeCoordinateAffected = is_free_component && unsupported_contrast,
        ParameterStatus = status,
        ReasonCodes = paste(unique(transition_reason), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
    transition_scope <- do.call(rbind, transition_scope_rows)
    transition_rows[[scope_index]] <- transition_scope

    contrast_scope <- lapply(unsupported_category_indices, function(category_index) {
      lower_step <- category_index - 1L
      upper_step <- category_index
      affected_free <- intersect(
        c(lower_step, upper_step),
        seq_len(free_per_scope)
      )
      data.frame(
        ReadinessScope = "parameter",
        ScaleScope = "single_observed_scale",
        StepScope = scope$label,
        ParameterClass = "step_contrast",
        Category = declared[category_index],
        Direction = paste0("Step", lower_step, " - Step", upper_step),
        PlusStep = paste0("Step", lower_step),
        MinusStep = paste0("Step", upper_step),
        FreeCoordinates = mfrmr_category_values_text(
          paste0("FreeStep", affected_free)
        ),
        ParameterStatus = "unsupported",
        ReasonCodes = paste(
          "step_scope_category_unobserved",
          "unsupported_step_coordinate",
          sep = ";"
        ),
        stringsAsFactors = FALSE
      )
    })
    contrast_rows[[scope_index]] <- if (length(contrast_scope)) {
      do.call(rbind, contrast_scope)
    } else {
      NULL
    }

    has_unsupported <- length(unsupported_category_indices) > 0L
    observed_counts <- counts[weighted_counts > 0]
    scope_weak <- !has_unsupported && (
      any(weighted_counts <= 0) || any(observed_counts <= 1L) ||
        any(weak_transition)
    )
    information_state <- if (has_unsupported) {
      "unsupported_coordinate"
    } else if (scope_weak) {
      "weak_information"
    } else {
      "adequate"
    }
    reason_codes <- character(0)
    if (length(global_missing) > 0L) {
      reason_codes <- c(reason_codes, "declared_category_unobserved_global")
    }
    if (any(global_missing %in% declared[c(1L, n_categories)])) {
      reason_codes <- c(reason_codes, "declared_boundary_category_unobserved")
    }
    if (!isTRUE(scope$shared) && length(missing) > 0L) {
      reason_codes <- c(reason_codes, "step_scope_category_unobserved")
    }
    if (has_unsupported) {
      reason_codes <- c(
        reason_codes,
        "step_scope_category_unobserved",
        "unsupported_step_coordinate"
      )
    } else if (scope_weak) {
      reason_codes <- c(reason_codes, "weak_category_information")
    }

    unsupported_index <- which(
      vapply(seq_len(n_steps), function(transition) {
        any(unsupported_category_indices %in% c(transition, transition + 1L))
      }, logical(1))
    )
    unsupported_labels <- if (length(unsupported_index) > 0L) {
      paste0("Step", unsupported_index)
    } else {
      character(0)
    }
    support_rows[[scope_index]] <- data.frame(
      ScaleScope = "single_observed_scale",
      StepScope = scope$label,
      DeclaredCategories = mfrmr_category_values_text(declared),
      ObservedGlobal = mfrmr_category_values_text(global_observed),
      ObservedWithinScope = mfrmr_category_values_text(observed),
      RetainedForFit = mfrmr_category_values_text(declared),
      FreeStepCount = free_per_scope,
      FixedStepCount = 0L,
      DerivedStepCount = n_steps - free_per_scope,
      UnsupportedFreeStepCount = sum(
        seq_len(free_per_scope) %in% unsupported_index
      ),
      UnsupportedCategory = mfrmr_category_values_text(
        declared[unsupported_category_indices]
      ),
      UnsupportedStep = mfrmr_category_values_text(unsupported_labels),
      ZeroType = mfrmr_category_values_text(unique(zero_types[zero_types != "none"])),
      MinimumObservedCategoryCount = if (length(observed_counts)) {
        min(observed_counts)
      } else {
        0L
      },
      MaximumCategoryFraction = if (sum(counts) > 0L) {
        max(counts) / sum(counts)
      } else {
        NA_real_
      },
      MaximumWeightedCategoryFraction = if (sum(weighted_counts) > 0) {
        max(weighted_counts) / sum(weighted_counts)
      } else {
        NA_real_
      },
      NormalizedCategoryEntropy = mfrmr_category_normalized_entropy(counts),
      WeightedNormalizedCategoryEntropy = mfrmr_category_normalized_entropy(
        weighted_counts
      ),
      InformationState = information_state,
      ReasonCode = paste(unique(reason_codes), collapse = ";"),
      stringsAsFactors = FALSE
    )

    category_rows[[scope_index]] <- data.frame(
      ScaleScope = "single_observed_scale",
      StepScope = scope$label,
      Category = declared,
      BoundaryCategory = seq_along(declared) %in% c(1L, n_categories),
      GlobalCount = global_counts,
      GlobalWeightedN = global_weighted,
      WithinScopeCount = counts,
      WithinScopeWeightedN = weighted_counts,
      ObservedGlobalRaw = global_counts > 0L,
      ObservedWithinScopeRaw = counts > 0L,
      ObservedGlobal = global_weighted > 0,
      ObservedWithinScope = weighted_counts > 0,
      RetainedForFit = TRUE,
      ZeroType = zero_types,
      SingletonObserved = counts == 1L,
      stringsAsFactors = FALSE
    )
  }

  support_table <- do.call(rbind, support_rows)
  category_table <- do.call(rbind, category_rows)
  step_status <- do.call(rbind, transition_rows)
  unsupported_contrasts <- do.call(rbind, contrast_rows)
  if (is.null(unsupported_contrasts)) {
    unsupported_contrasts <- data.frame(
      ReadinessScope = character(0),
      ScaleScope = character(0),
      StepScope = character(0),
      ParameterClass = character(0),
      Category = integer(0),
      Direction = character(0),
      PlusStep = character(0),
      MinusStep = character(0),
      FreeCoordinates = character(0),
      ParameterStatus = character(0),
      ReasonCodes = character(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(support_table) <- NULL
  rownames(category_table) <- NULL
  rownames(step_status) <- NULL
  rownames(unsupported_contrasts) <- NULL

  local_support <- mfrmr_category_local_support(prep, declared)
  rsm_local_weak <- identical(config$model, "RSM") &&
    nrow(local_support) > 0L &&
    any(local_support$InformationState == "weak_information")
  scope_states <- as.character(support_table$InformationState)
  state <- if (any(scope_states == "unsupported_coordinate")) {
    "unsupported_coordinate"
  } else if (any(scope_states == "weak_information") || rsm_local_weak) {
    "weak_information"
  } else {
    "adequate"
  }
  audit_reasons <- unique(unlist(strsplit(
    support_table$ReasonCode[nzchar(support_table$ReasonCode)],
    ";", fixed = TRUE
  )))
  if (rsm_local_weak) {
    audit_reasons <- unique(c(audit_reasons, "weak_category_information"))
  }
  audit_reasons <- audit_reasons[nzchar(audit_reasons)]
  unsupported_coordinates <- step_status$Coordinate[
    step_status$ParameterStatus == "unsupported"
  ]

  readiness <- data.frame(
    ReadinessContractVersion = mfrmr_category_contract_version(),
    ReadinessScope = "fit",
    CategoryState = state,
    ReasonCodes = paste(audit_reasons, collapse = ";"),
    Complete = TRUE,
    StepScopes = nrow(support_table),
    UnsupportedStepCoordinates = length(unsupported_coordinates),
    UnsupportedCategoryContrasts = nrow(unsupported_contrasts),
    WeakStepScopes = sum(scope_states == "weak_information"),
    WeakLocalScopes = if (nrow(local_support) > 0L) {
      sum(local_support$InformationState == "weak_information")
    } else {
      0L
    },
    stringsAsFactors = FALSE
  )

  list(
    contract_version = mfrmr_category_contract_version(),
    model = as.character(config$model),
    method = as.character(config$method),
    scale_scope = "single_observed_scale",
    step_facet = as.character(config$step_facet %||% NA_character_),
    rating_range_source = as.character(prep$rating_range_source %||% NA_character_),
    score_map = as.data.frame(prep$score_map, stringsAsFactors = FALSE),
    declared_categories = declared,
    observed_global = global_observed,
    retained_categories = declared,
    readiness = readiness,
    support_table = support_table,
    category_table = category_table,
    step_status = step_status,
    unsupported_contrasts = unsupported_contrasts,
    local_support = local_support,
    unsupported_step_coordinates = as.character(unsupported_coordinates),
    complete = TRUE,
    exact_support_rule = paste(
      "A free step contrast is unsupported when an internal category has",
      "zero positive-weight observations within its fitted ladder scope;",
      "the recession direction increases the step below that category and",
      "decreases the step above it while leaving other category exponents unchanged."
    ),
    weak_information_rule = paste(
      "Empty boundary categories, categories absent only in a non-step local",
      "facet scope, singleton observed cells, and singleton transition sides",
      "are review evidence only;",
      "concentration thresholds remain pending pilot calibration."
    )
  )
}

mfrmr_category_condition <- function(message, audit,
                                      type = c("error", "warning")) {
  type <- match.arg(type)
  subclass <- if (identical(type, "error")) {
    "mfrmr_category_readiness_error"
  } else {
    "mfrmr_category_readiness_warning"
  }
  structure(
    list(
      message = as.character(message),
      call = NULL,
      readiness = audit$readiness,
      category_support = audit
    ),
    class = c(
      subclass,
      "mfrmr_readiness_condition",
      type,
      "condition"
    )
  )
}

mfrmr_stop_unsupported_category <- function(audit) {
  unsupported <- nrow(audit$unsupported_contrasts %||% data.frame())
  scopes <- sum(
    audit$support_table$InformationState == "unsupported_coordinate",
    na.rm = TRUE
  )
  message <- paste0(
    "The declared ", audit$model, " category structure contains ",
    unsupported, " unsupported internal-category step contrast(s) across ",
    scopes, " fitted step scope(s). Each affected internal category has zero ",
    "positive-weight observations, creating a recession direction in the ",
    "sum-zero step parameterization. Optimization was not run. Collect ",
    "observations in every affected internal category, revise the score ",
    "support or response model, or separate the fitted ladder; threshold ",
    "anchoring is not available ",
    "in mfrmr 0.2.3."
  )
  stop(mfrmr_category_condition(message, audit, type = "error"))
}

mfrmr_warn_category_support <- function(audit) {
  message <- paste(
    "Category support is retained but requires review: at least one fitted or",
    "local scope contains an empty or singleton category/transition cell.",
    "The fit may be inspected, but category-information strength has not been",
    "certified; inspect `fit$data_review$category_support` before inference."
  )
  warning(mfrmr_category_condition(message, audit, type = "warning"))
  invisible(audit)
}
