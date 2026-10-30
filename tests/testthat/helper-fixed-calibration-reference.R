# Shared declarations for current normalization and independent typed-anchor checks.

g2_mock_config <- function(model = c("RSM", "PCM"), reversed = FALSE) {
  model <- match.arg(model)
  facet_levels <- list(
    Rater = c("R1", "R2", "R3"),
    Criterion = c("C1", "C2", "C3")
  )
  config <- list(
    model = model,
    method = "MML",
    n_cat = 4L,
    facet_names = names(facet_levels),
    facet_levels = facet_levels,
    step_facet = if (identical(model, "PCM")) "Criterion" else NULL,
    interaction_specs = list(),
    facet_interactions = character(0),
    population_spec = list(active = FALSE),
    facet_signs = c(Rater = -1, Criterion = -1),
    score_map = data.frame(
      OriginalScore = if (reversed) 3:0 else 0:3,
      InternalScore = 0:3,
      stringsAsFactors = FALSE
    )
  )
  config$facet_specs <- lapply(facet_levels, function(levels) {
    mfrmr:::build_facet_constraint(levels)
  })
  config$step_specs <- mfrmr:::mfrmr_default_step_specs(config)
  config
}

g2_anchor <- function(type, parameter_class, owner = NA_character_,
                      level = NA_character_, step = NA_character_,
                      group = NA_character_, value = 0, id = "input::1",
                      order = 1L, coordinate_system = "expanded_logit") {
  data.frame(
    AnchorId = id,
    AnchorType = type,
    ParameterClass = parameter_class,
    OwnerFacet = owner,
    Level = level,
    Step = step,
    GroupId = group,
    Value = as.numeric(value),
    CoordinateSystem = coordinate_system,
    DeclarationOrder = as.integer(order),
    stringsAsFactors = FALSE
  )
}

g2_bind_anchors <- function(...) {
  rows <- list(...)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$AnchorId <- sprintf("input::%d", seq_len(nrow(out)))
  out$DeclarationOrder <- as.integer(seq_len(nrow(out)))
  out
}

g2_capture_error <- function(expr) {
  tryCatch(force(expr), mfrm_calibration_error = function(error) error)
}
