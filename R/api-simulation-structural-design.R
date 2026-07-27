# ==============================================================================
# Named-facet design-schema layer for structural simulation review
# ==============================================================================
#
# Internal helpers that build, normalize, and report the deterministic
# design schema used by the simulation engine (`api-simulation.R`). Split out
# so the design-schema structural review machinery has a dedicated home; the
# main public simulator entry
# (`simulate_mfrm_data`) and the design / diagnostic / signal
# evaluation entries (`evaluate_mfrm_*`) remain in `api-simulation.R`.
# Functions here are internal (no @export); they are called from
# `evaluate_mfrm_design()` and the surrounding orchestration helpers.

simulation_structural_design_contract <- function(sim_spec = NULL, facet_names = NULL) {
  structural_facet_table <- if (is.null(facet_names)) {
    simulation_structural_facet_table(sim_spec)
  } else {
    simulation_structural_facet_table(facet_names = facet_names)
  }
  structural_design_template <- if (is.null(facet_names)) {
    simulation_structural_design_template(sim_spec)
  } else {
    simulation_structural_design_template(facet_names = facet_names)
  }

  alias_source <- if (!is.null(sim_spec)) {
    sim_spec
  } else if (!is.null(facet_names)) {
    list(facet_names = simulation_validate_output_facet_names(facet_names))
  } else {
    NULL
  }

  facet_axes <- structural_facet_table |>
    dplyr::transmute(
      input_key = .data$structural_facet_key,
      facet = .data$facet,
      facet_kind = .data$facet_kind,
      axis_class = .data$design_axis_class,
      level_count = .data$level_count,
      canonical_design_variable = .data$current_planning_count_variable,
      public_design_alias = .data$current_planning_count_alias,
      current_planner_role_supported = .data$current_planner_role_supported,
      arbitrary_facet_design_candidate = .data$arbitrary_facet_design_candidate
    )

  assignment_axis <- list(
    current_planning_count_variable = "raters_per_person",
    current_planning_count_alias = unname(simulation_design_variable_aliases(alias_source)[["raters_per_person"]]),
    design_input_key = "assignment",
    axis_class = "assignment_count",
    depends_on_input_key = as.character(facet_axes$input_key[match("facet_level_count", facet_axes$axis_class)] %||% "rater"),
    depends_on_canonical_design_variable = as.character(
      facet_axes$canonical_design_variable[match("facet_level_count", facet_axes$axis_class)] %||% "n_rater"
    ),
    depends_on_public_design_alias = as.character(
      facet_axes$public_design_alias[match("facet_level_count", facet_axes$axis_class)] %||% "n_rater"
    )
  )

  schema <- list(
    schema_contract = "arbitrary_facet_design_schema",
    schema_stage = "deterministic_review",
    facet_axes = facet_axes,
    assignment_axis = assignment_axis,
    input_keys = c(as.character(facet_axes$input_key), assignment_axis$design_input_key),
    canonical_design_variables = c(
      as.character(facet_axes$canonical_design_variable),
      assignment_axis$current_planning_count_variable
    ),
    default_design = structural_design_template,
    note = paste(
      "Deterministic design-schema object for the named-facet structural review,",
      "bundling stable facet-count axes and the assignment axis in one",
      "machine-readable contract."
    )
  )
  schema$grid_semantics <- simulation_structural_design_grid_semantics(schema)
  schema
}

simulation_structural_design_grid_semantics <- function(design_schema) {
  design_schema <- simulation_coerce_structural_design_schema(design_schema)
  facet_axes <- design_schema$facet_axes
  assignment_axis <- design_schema$assignment_axis
  assignment_var <- as.character(
    assignment_axis$current_planning_count_variable %||% "raters_per_person"
  )
  assignment_key <- as.character(
    assignment_axis$design_input_key %||% "assignment"
  )
  canonical_columns <- c(
    "design_id",
    as.character(facet_axes$canonical_design_variable),
    assignment_var
  )
  public_columns <- c(
    "design_id",
    as.character(facet_axes$public_design_alias),
    as.character(assignment_axis$current_planning_count_alias)
  )
  branch_columns <- c(
    "design_id",
    as.character(facet_axes$input_key),
    assignment_key
  )

  list(
    semantics_contract = "arbitrary_facet_design_grid_semantics",
    semantics_stage = as.character(design_schema$schema_stage %||% "deterministic_review"),
    id_variable = "design_id",
    canonical_columns = canonical_columns,
    public_columns = public_columns,
    branch_columns = branch_columns,
    feasibility_rule = paste0(
      assignment_var,
      " <= ",
      as.character(assignment_axis$depends_on_canonical_design_variable %||% "n_rater")
    ),
    row_meaning = paste(
      "Each row is one deterministic structural-design condition formed by",
      "crossing facet-count axes and the assignment axis after applying the",
      "current feasibility rule."
    ),
    default_id_prefix = "F"
  )
}

simulation_coerce_structural_design_schema <- function(structural_design_schema) {
  required_fields <- c(
    "schema_contract",
    "schema_stage",
    "facet_axes",
    "assignment_axis",
    "input_keys",
    "canonical_design_variables",
    "default_design",
    "note"
  )
  attach_grid_semantics <- function(x) {
    if (is.list(x$grid_semantics) &&
        identical(x$grid_semantics$semantics_contract, "arbitrary_facet_design_grid_semantics")) {
      return(c(x[required_fields], list(grid_semantics = x$grid_semantics)))
    }

    facet_axes <- x$facet_axes
    assignment_axis <- x$assignment_axis
    assignment_var <- as.character(
      assignment_axis$current_planning_count_variable %||% "raters_per_person"
    )
    assignment_key <- as.character(
      assignment_axis$design_input_key %||% "assignment"
    )

    c(
      x[required_fields],
      list(
        grid_semantics = list(
          semantics_contract = "arbitrary_facet_design_grid_semantics",
          semantics_stage = as.character(x$schema_stage %||% "deterministic_review"),
          id_variable = "design_id",
          canonical_columns = c(
            "design_id",
            as.character(facet_axes$canonical_design_variable),
            assignment_var
          ),
          public_columns = c(
            "design_id",
            as.character(facet_axes$public_design_alias),
            as.character(assignment_axis$current_planning_count_alias)
          ),
          branch_columns = c(
            "design_id",
            as.character(facet_axes$input_key),
            assignment_key
          ),
          feasibility_rule = paste0(
            assignment_var,
            " <= ",
            as.character(assignment_axis$depends_on_canonical_design_variable %||% "n_rater")
          ),
          row_meaning = paste(
            "Each row is one deterministic structural-design condition formed by",
            "crossing facet-count axes and the assignment axis after applying the",
            "current feasibility rule."
          ),
          default_id_prefix = "F"
        )
      )
    )
  }
  if (is.list(structural_design_schema) &&
      all(required_fields %in% names(structural_design_schema)) &&
      is.data.frame(structural_design_schema$facet_axes) &&
      is.list(structural_design_schema$assignment_axis) &&
      is.list(structural_design_schema$default_design)) {
    return(attach_grid_semantics(structural_design_schema))
  }

  design_schema <- structural_design_schema$design_schema %||% NULL
  if (is.list(design_schema) &&
      all(required_fields %in% names(design_schema)) &&
      is.data.frame(design_schema$facet_axes) &&
      is.list(design_schema$assignment_axis) &&
      is.list(design_schema$default_design)) {
    return(attach_grid_semantics(design_schema))
  }

  if (!is.list(structural_design_schema) ||
      !is.data.frame(structural_design_schema$facet_table) ||
      !is.list(structural_design_schema$assignment_axis) ||
      !is.list(structural_design_schema$design_template)) {
    stop("`structural_design_schema` is not a valid deterministic arbitrary-facet design contract.",
         call. = FALSE)
  }

  facet_axes <- structural_design_schema$facet_table |>
    dplyr::transmute(
      input_key = .data$structural_facet_key,
      facet = .data$facet,
      facet_kind = .data$facet_kind,
      axis_class = .data$design_axis_class,
      level_count = .data$level_count,
      canonical_design_variable = .data$current_planning_count_variable,
      public_design_alias = .data$current_planning_count_alias,
      current_planner_role_supported = .data$current_planner_role_supported,
      arbitrary_facet_design_candidate = .data$arbitrary_facet_design_candidate
    )
  assignment_axis <- structural_design_schema$assignment_axis

  attach_grid_semantics(list(
    schema_contract = "arbitrary_facet_design_schema",
    schema_stage = as.character(structural_design_schema$planner_stage %||% "deterministic_review"),
    facet_axes = facet_axes,
    assignment_axis = assignment_axis,
    input_keys = c(as.character(facet_axes$input_key), assignment_axis$design_input_key),
    canonical_design_variables = c(
      as.character(facet_axes$canonical_design_variable),
      assignment_axis$current_planning_count_variable
    ),
    default_design = structural_design_schema$design_template,
    note = paste(
      "Compatibility coercion of a deterministic named-facet structural review",
      "contract into the bundled design-schema object."
    )
  ))
}

simulation_build_structural_design_grid <- function(design_schema,
                                                       design = NULL,
                                                       id_prefix = "F") {
  design_schema <- simulation_coerce_structural_design_schema(design_schema)
  parsed <- simulation_parse_structural_design(
    design = design,
    structural_design_schema = list(design_schema = design_schema),
    arg_name = "design"
  )

  defaults <- design_schema$default_design %||% list()
  default_facets <- defaults$facets %||% list()
  facet_axes <- design_schema$facet_axes
  assignment_axis <- design_schema$assignment_axis
  grid_semantics <- design_schema$grid_semantics

  values <- list()
  for (i in seq_len(nrow(facet_axes))) {
    canonical <- as.character(facet_axes$canonical_design_variable[i])
    input_key <- as.character(facet_axes$input_key[i])
    value <- parsed[[canonical]] %||% default_facets[[input_key]]
    values[[canonical]] <- simulation_validate_count_values(value, canonical, min_value = 2L)
  }

  assignment_var <- as.character(assignment_axis$current_planning_count_variable %||% "raters_per_person")
  assignment_key <- as.character(assignment_axis$design_input_key %||% "assignment")
  assignment_default <- defaults[[assignment_key]] %||%
    values[[as.character(assignment_axis$depends_on_canonical_design_variable %||% "n_rater")]]
  values[[assignment_var]] <- simulation_validate_count_values(
    parsed[[assignment_var]] %||% assignment_default,
    assignment_var,
    min_value = 1L
  )

  canonical_order <- setdiff(
    as.character(grid_semantics$canonical_columns %||% c("design_id", as.character(facet_axes$canonical_design_variable), assignment_var)),
    "design_id"
  )
  design_grid <- expand.grid(
    values[canonical_order],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  dependency_var <- as.character(
    assignment_axis$depends_on_canonical_design_variable %||% "n_rater"
  )
  design_grid <- design_grid[
    design_grid[[assignment_var]] <= design_grid[[dependency_var]],
    ,
    drop = FALSE
  ]
  if (nrow(design_grid) == 0L) {
    stop(
      "No valid structural-design rows remain after enforcing `",
      assignment_var,
      " <= ",
      dependency_var,
      "`.",
      call. = FALSE
    )
  }

  prefix <- as.character(grid_semantics$default_id_prefix %||% "F")
  design_grid$design_id <- sprintf("%s%02d", as.character(id_prefix[1] %||% prefix), seq_len(nrow(design_grid)))
  design_grid <- tibble::as_tibble(design_grid[, c("design_id", canonical_order), drop = FALSE])

  branch_grid <- tibble::tibble(design_id = design_grid$design_id)
  for (i in seq_len(nrow(facet_axes))) {
    branch_grid[[as.character(facet_axes$input_key[i])]] <-
      design_grid[[as.character(facet_axes$canonical_design_variable[i])]]
  }
  branch_grid[[assignment_key]] <- design_grid[[assignment_var]]

  aliases <- stats::setNames(
    c(as.character(facet_axes$public_design_alias), as.character(assignment_axis$current_planning_count_alias)),
    canonical_order
  )

  list(
    canonical = design_grid,
    public = simulation_append_design_alias_columns(design_grid, aliases),
    design = branch_grid,
    design_schema = design_schema,
    grid_semantics = grid_semantics
  )
}

simulation_structural_design_preview <- function(structural_design_schema,
                                             design = NULL,
                                             id_prefix = "F") {
  grid_contract <- simulation_structural_design_grid_contract(
    structural_design_schema = structural_design_schema,
    design = design,
    id_prefix = id_prefix
  )
  preview_fields <- c(
    "preview_available",
    "reason",
    "default_design",
    "canonical",
    "public",
    "design",
    "design_schema",
    "grid_semantics"
  )

  names(grid_contract)[names(grid_contract) == "contract"] <- "preview_contract"
  names(grid_contract)[names(grid_contract) == "stage"] <- "preview_stage"
  grid_contract[c("preview_contract", "preview_stage", preview_fields)]
}

simulation_structural_design_grid_contract <- function(structural_design_schema,
                                                   design = NULL,
                                                   id_prefix = "F") {
  design_schema <- simulation_coerce_structural_design_schema(structural_design_schema)
  default_design <- design_schema$default_design %||% list()
  grid_semantics <- design_schema$grid_semantics
  facet_defaults <- default_design$facets %||% list()
  assignment_key <- as.character(
    design_schema$assignment_axis$design_input_key %||% "assignment"
  )
  assignment_default <- default_design[[assignment_key]] %||% NULL

  facet_complete <- length(facet_defaults) == nrow(design_schema$facet_axes) &&
    all(vapply(
      facet_defaults,
      function(x) length(x) == 1L && is.numeric(x) && is.finite(x),
      logical(1)
    ))
  assignment_complete <- length(assignment_default) == 1L &&
    is.numeric(assignment_default) && is.finite(assignment_default)

  preview_design <- design %||% default_design
  if (!facet_complete || !assignment_complete) {
    return(list(
      contract = "arbitrary_facet_design_grid_contract",
      stage = as.character(design_schema$schema_stage %||% "deterministic_review"),
      planner_contract = as.character(
        structural_design_schema$planner_contract %||% "named_facet_structural_review"
      ),
      input_contract = as.character(
        structural_design_schema$input_contract %||% "design$facets(named counts)"
      ),
      preview_available = FALSE,
      reason = paste(
        "The deterministic named-facet structural review does not yet have",
        "complete default facet counts and assignment counts, so no",
        "deterministic design grid can be materialized."
      ),
      default_design = default_design,
      design_schema = design_schema,
      grid_semantics = grid_semantics,
      note = paste(
        "Deterministic structural-design grid contract is unavailable until the",
        "default nested design carries finite facet and assignment counts."
      )
    ))
  }

  preview_grid <- simulation_build_structural_design_grid(
    design_schema = design_schema,
    design = preview_design,
    id_prefix = id_prefix
  )

  list(
    contract = "arbitrary_facet_design_grid_contract",
    stage = as.character(design_schema$schema_stage %||% "deterministic_review"),
    planner_contract = as.character(
      structural_design_schema$planner_contract %||% "named_facet_structural_review"
    ),
    input_contract = as.character(
      structural_design_schema$input_contract %||% "design$facets(named counts)"
    ),
    preview_available = TRUE,
    reason = paste(
      "Deterministic design grid built from the future arbitrary-facet",
      "branch design schema and its current default design."
    ),
    default_design = default_design,
    canonical = preview_grid$canonical,
    public = preview_grid$public,
    design = preview_grid$design,
    design_schema = design_schema,
    grid_semantics = preview_grid$grid_semantics,
    note = paste(
      "Deterministic structural-design grid contract bundling the default branch",
      "grid, design schema, and grid semantics without activating arbitrary-",
      "facet planner logic."
    )
  )
}

simulation_coerce_structural_design_grid_contract <- function(x) {
  required_fields <- c(
    "contract",
    "stage",
    "planner_contract",
    "input_contract",
    "preview_available",
    "reason",
    "default_design",
    "design_schema",
    "grid_semantics",
    "note"
  )
  optional_grid_fields <- c("canonical", "public", "design")

  if (is.list(x) &&
      identical(x$contract %||% "", "arbitrary_facet_design_grid_contract") &&
      all(required_fields %in% names(x)) &&
      is.list(x$design_schema) &&
      is.list(x$grid_semantics)) {
    out <- x[unique(c(required_fields, optional_grid_fields[optional_grid_fields %in% names(x)]))]
    return(out)
  }

  nested <- x$grid_contract %||% x$structural_design_grid_contract %||% NULL
  if (is.list(nested)) {
    return(simulation_coerce_structural_design_grid_contract(nested))
  }

  branch <- x$structural_design_schema %||% x
  if (!is.list(branch)) {
    stop("`x` is not a valid deterministic structural-design grid contract.", call. = FALSE)
  }

  simulation_structural_design_grid_contract(branch)
}

simulation_build_structural_design_grid_from_contract <- function(grid_contract,
                                                                     design = NULL,
                                                                     id_prefix = NULL) {
  grid_contract <- simulation_coerce_structural_design_grid_contract(grid_contract)
  prefix <- as.character(id_prefix %||% grid_contract$grid_semantics$default_id_prefix %||% "F")

  if (is.null(design) &&
      isTRUE(grid_contract$preview_available) &&
      all(c("canonical", "public", "design") %in% names(grid_contract))) {
    return(list(
      canonical = grid_contract$canonical,
      public = grid_contract$public,
      design = grid_contract$design,
      design_schema = grid_contract$design_schema,
      grid_semantics = grid_contract$grid_semantics,
      grid_contract = grid_contract
    ))
  }

  built <- simulation_build_structural_design_grid(
    design_schema = grid_contract$design_schema,
    design = design %||% grid_contract$default_design,
    id_prefix = prefix
  )
  built$grid_contract <- grid_contract
  built
}

simulation_materialize_structural_design_grid <- function(x,
                                                      design = NULL,
                                                      id_prefix = NULL) {
  grid_contract <- NULL

  if (is.list(x) &&
      identical(x$contract %||% "", "arbitrary_facet_design_grid_contract")) {
    grid_contract <- x
  }

  if (is.null(grid_contract) && is.list(x)) {
    nested <- x$structural_design_grid_contract %||%
      x$grid_contract %||%
      x$structural_design_schema %||%
      x$planning_schema %||%
      x$settings$planning_schema %||%
      x$sim_spec$planning_schema %||%
      x$ademp$data_generating_mechanism$planning_schema %||%
      NULL
    if (is.list(nested)) {
      grid_contract <- simulation_coerce_structural_design_grid_contract(nested)
    }
  }

  if (is.null(grid_contract)) {
    schema <- tryCatch(
      simulation_object_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$structural_design_grid_contract)) {
      grid_contract <- schema$structural_design_grid_contract
    }
  }

  if (is.null(grid_contract)) {
    schema <- tryCatch(
      simulation_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$structural_design_grid_contract)) {
      grid_contract <- schema$structural_design_grid_contract
    }
  }

  if (is.null(grid_contract)) {
    stop(
      "`x` is not a valid planning object, planning schema, sim spec, or ",
      "structural-design grid contract.",
      call. = FALSE
    )
  }

  simulation_build_structural_design_grid_from_contract(
    grid_contract = grid_contract,
    design = design,
    id_prefix = id_prefix
  )
}

simulation_structural_design_grid_bundle <- function(x,
                                                 design = NULL,
                                                 id_prefix = NULL) {
  grid_contract <- simulation_coerce_structural_design_grid_contract(x)

  if (is.null(design) && !isTRUE(grid_contract$preview_available)) {
    return(list(
      bundle_contract = "arbitrary_facet_design_grid_bundle",
      bundle_stage = as.character(grid_contract$stage %||% "deterministic_review"),
      planner_contract = as.character(
        grid_contract$planner_contract %||% "named_facet_structural_review"
      ),
      input_contract = as.character(
        grid_contract$input_contract %||% "design$facets(named counts)"
      ),
      grid_available = FALSE,
      reason = as.character(
        grid_contract$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      default_design = grid_contract$default_design,
      design_schema = grid_contract$design_schema,
      grid_semantics = grid_contract$grid_semantics,
      grid_contract = grid_contract,
      note = paste(
        "Deterministic structural-design bundle carrying the design grid contract,",
        "schema, and semantics without a materialized grid because default",
        "counts are not yet available."
      )
    ))
  }

  built <- simulation_build_structural_design_grid_from_contract(
    grid_contract = grid_contract,
    design = design,
    id_prefix = id_prefix
  )

  list(
    bundle_contract = "arbitrary_facet_design_grid_bundle",
    bundle_stage = as.character(grid_contract$stage %||% "deterministic_review"),
    planner_contract = as.character(
      grid_contract$planner_contract %||% "named_facet_structural_review"
    ),
    input_contract = as.character(
      grid_contract$input_contract %||% "design$facets(named counts)"
    ),
    grid_available = TRUE,
    reason = paste(
      "Deterministic structural-design bundle materialized from the design grid",
      "contract and matching design schema."
    ),
    canonical = built$canonical,
    public = built$public,
    design = built$design,
    default_design = grid_contract$default_design,
    design_schema = built$design_schema,
    grid_semantics = built$grid_semantics,
    grid_contract = grid_contract,
    note = paste(
      "Deterministic structural-design bundle that materializes canonical, public,",
      "and branch-facing design grids from one authoritative branch-grid",
      "contract."
    )
  )
}

simulation_coerce_structural_design_grid_bundle <- function(x) {
  required_fields <- c(
    "bundle_contract",
    "bundle_stage",
    "planner_contract",
    "input_contract",
    "grid_available",
    "reason",
    "default_design",
    "design_schema",
    "grid_semantics",
    "grid_contract",
    "note"
  )
  optional_grid_fields <- c("canonical", "public", "design")

  if (is.list(x) &&
      identical(x$bundle_contract %||% "", "arbitrary_facet_design_grid_bundle") &&
      all(required_fields %in% names(x)) &&
      is.list(x$design_schema) &&
      is.list(x$grid_semantics) &&
      is.list(x$grid_contract)) {
    return(x[unique(c(required_fields, optional_grid_fields[optional_grid_fields %in% names(x)]))])
  }

  nested <- x$grid_bundle %||% x$structural_design_grid_bundle %||% NULL
  if (is.list(nested)) {
    return(simulation_coerce_structural_design_grid_bundle(nested))
  }

  simulation_structural_design_grid_bundle(x)
}

simulation_materialize_structural_design_grid_bundle <- function(x,
                                                             design = NULL,
                                                             id_prefix = NULL) {
  grid_bundle <- NULL

  if (is.list(x) &&
      identical(x$bundle_contract %||% "", "arbitrary_facet_design_grid_bundle")) {
    grid_bundle <- x
  }

  if (is.null(grid_bundle) && is.list(x)) {
    nested <- x$structural_design_grid_bundle %||%
      x$grid_bundle %||%
      x$structural_design_schema %||%
      x$planning_schema %||%
      x$settings$planning_schema %||%
      x$sim_spec$planning_schema %||%
      x$ademp$data_generating_mechanism$planning_schema %||%
      NULL
    if (is.list(nested)) {
      grid_bundle <- simulation_coerce_structural_design_grid_bundle(nested)
    }
  }

  if (is.null(grid_bundle)) {
    schema <- tryCatch(
      simulation_object_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$structural_design_grid_bundle)) {
      grid_bundle <- schema$structural_design_grid_bundle
    }
  }

  if (is.null(grid_bundle)) {
    schema <- tryCatch(
      simulation_planning_schema(x),
      error = function(e) NULL
    )
    if (is.list(schema) && is.list(schema$structural_design_grid_bundle)) {
      grid_bundle <- schema$structural_design_grid_bundle
    }
  }

  if (is.null(grid_bundle)) {
    stop(
      "`x` is not a valid planning object, planning schema, sim spec, or ",
      "structural-design grid bundle.",
      call. = FALSE
    )
  }

  if (is.null(design) && isTRUE(grid_bundle$grid_available)) {
    return(grid_bundle)
  }

  simulation_structural_design_grid_bundle(
    x = grid_bundle$grid_contract %||% grid_bundle,
    design = design,
    id_prefix = id_prefix
  )
}

simulation_structural_design_grid_view <- function(x,
                                               view = c("canonical", "public", "design"),
                                               design = NULL,
                                               id_prefix = NULL) {
  view <- match.arg(view)
  grid_bundle <- simulation_materialize_structural_design_grid_bundle(
    x = x,
    design = design,
    id_prefix = id_prefix
  )

  if (!isTRUE(grid_bundle$grid_available) || !is.data.frame(grid_bundle[[view]])) {
    stop(
      "The requested `", view, "` structural-design grid view is not currently ",
      "available from this deterministic bundle.",
      call. = FALSE
    )
  }

  grid_bundle[[view]]
}

simulation_structural_design_grid_context <- function(x,
                                                  design = NULL,
                                                  id_prefix = NULL) {
  grid_bundle <- simulation_materialize_structural_design_grid_bundle(
    x = x,
    design = design,
    id_prefix = id_prefix
  )

  design_schema <- grid_bundle$design_schema
  grid_semantics <- grid_bundle$grid_semantics
  facet_axes <- design_schema$facet_axes
  assignment_axis <- design_schema$assignment_axis

  axis_table <- dplyr::bind_rows(
    facet_axes |>
      dplyr::transmute(
        axis_source = "facet",
        input_key = .data$input_key,
        canonical_design_variable = .data$canonical_design_variable,
        public_design_alias = .data$public_design_alias,
        axis_class = .data$axis_class,
        facet = .data$facet
      ),
    tibble::tibble(
      axis_source = "assignment",
      input_key = as.character(assignment_axis$design_input_key %||% "assignment"),
      canonical_design_variable = as.character(
        assignment_axis$current_planning_count_variable %||% "raters_per_person"
      ),
      public_design_alias = as.character(
        assignment_axis$current_planning_count_alias %||% "raters_per_person"
      ),
      axis_class = as.character(assignment_axis$axis_class %||% "assignment_count"),
      facet = NA_character_
    )
  )

  if (!isTRUE(grid_bundle$grid_available) || !is.data.frame(grid_bundle$canonical)) {
    axis_table$values <- vector("list", nrow(axis_table))
    axis_table$n_values <- rep(NA_integer_, nrow(axis_table))
    axis_table$varying <- rep(NA, nrow(axis_table))
    axis_table$fixed_value <- rep(NA_integer_, nrow(axis_table))

    return(list(
      context_contract = "arbitrary_facet_design_grid_context",
      context_stage = as.character(grid_bundle$bundle_stage %||% "deterministic_review"),
      grid_available = FALSE,
      reason = as.character(grid_bundle$reason %||% "No structural-design grid is available."),
      id_variable = as.character(grid_semantics$id_variable %||% "design_id"),
      canonical_columns = as.character(grid_semantics$canonical_columns %||% character(0)),
      public_columns = as.character(grid_semantics$public_columns %||% character(0)),
      branch_columns = as.character(grid_semantics$branch_columns %||% character(0)),
      axis_table = axis_table,
      varying_canonical = character(0),
      fixed_canonical = character(0),
      varying_input_keys = character(0),
      fixed_input_keys = character(0),
      grid_bundle = grid_bundle,
      note = paste(
        "Deterministic structural-design context exposes axis metadata even when a",
        "materialized grid is not yet available."
      )
    ))
  }

  canonical_grid <- grid_bundle$canonical
  axis_values <- lapply(axis_table$canonical_design_variable, function(var) {
    sort(unique(as.integer(canonical_grid[[var]])))
  })
  n_values <- vapply(axis_values, length, integer(1))
  varying <- n_values > 1L
  fixed_value <- vapply(axis_values, function(vals) {
    if (length(vals) == 1L) vals[[1]] else NA_integer_
  }, integer(1))

  axis_table$values <- axis_values
  axis_table$n_values <- as.integer(n_values)
  axis_table$varying <- as.logical(varying)
  axis_table$fixed_value <- as.integer(fixed_value)

  list(
    context_contract = "arbitrary_facet_design_grid_context",
    context_stage = as.character(grid_bundle$bundle_stage %||% "deterministic_review"),
    grid_available = TRUE,
    reason = as.character(
      grid_bundle$reason %||% "Named-facet design grid context built from a materialized bundle."
    ),
    id_variable = as.character(grid_semantics$id_variable %||% "design_id"),
    canonical_columns = as.character(grid_semantics$canonical_columns %||% names(canonical_grid)),
    public_columns = as.character(grid_semantics$public_columns %||% names(grid_bundle$public %||% canonical_grid)),
    branch_columns = as.character(grid_semantics$branch_columns %||% names(grid_bundle$design %||% canonical_grid)),
    axis_table = axis_table,
    varying_canonical = as.character(axis_table$canonical_design_variable[axis_table$varying]),
    fixed_canonical = as.character(axis_table$canonical_design_variable[!axis_table$varying]),
    varying_input_keys = as.character(axis_table$input_key[axis_table$varying]),
    fixed_input_keys = as.character(axis_table$input_key[!axis_table$varying]),
    grid_bundle = grid_bundle,
    note = paste(
      "Deterministic structural-design context summarizing which branch axes vary",
      "within the currently materialized design grid."
    )
  )
}

simulation_structural_design_grid_summary <- function(x,
                                                  design = NULL,
                                                  id_prefix = NULL) {
  grid_context <- simulation_structural_design_grid_context(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  axis_table <- tibble::as_tibble(grid_context$axis_table)
  varying_rows <- !is.na(axis_table$varying) & as.logical(axis_table$varying)
  fixed_rows <- !is.na(axis_table$varying) & !as.logical(axis_table$varying)
  fixed_values <- stats::setNames(
    as.list(axis_table$fixed_value[fixed_rows]),
    as.character(axis_table$canonical_design_variable[fixed_rows])
  )

  if (!isTRUE(grid_context$grid_available)) {
    return(list(
      summary_contract = "arbitrary_facet_design_grid_summary",
      summary_stage = as.character(grid_context$context_stage %||% "deterministic_review"),
      planner_contract = as.character(
        grid_context$grid_bundle$planner_contract %||% "named_facet_structural_review"
      ),
      input_contract = as.character(
        grid_context$grid_bundle$input_contract %||% "design$facets(named counts)"
      ),
      grid_available = FALSE,
      reason = as.character(
        grid_context$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      n_designs = 0L,
      n_varying_axes = sum(varying_rows),
      n_fixed_axes = sum(fixed_rows),
      varying_canonical = as.character(grid_context$varying_canonical),
      fixed_canonical = as.character(grid_context$fixed_canonical),
      varying_input_keys = as.character(grid_context$varying_input_keys),
      fixed_input_keys = as.character(grid_context$fixed_input_keys),
      fixed_values = fixed_values,
      axis_table = axis_table,
      grid_context = grid_context,
      grid_bundle = grid_context$grid_bundle,
      note = paste(
        "Deterministic structural-design summary is unavailable because no",
        "materialized design grid exists yet, but axis metadata are still",
        "returned for design-review planning helpers."
      )
    ))
  }

  canonical <- tibble::as_tibble(grid_context$grid_bundle$canonical)
  public <- tibble::as_tibble(grid_context$grid_bundle$public)
  branch <- tibble::as_tibble(grid_context$grid_bundle$design)

  list(
    summary_contract = "arbitrary_facet_design_grid_summary",
    summary_stage = as.character(grid_context$context_stage %||% "deterministic_review"),
    planner_contract = as.character(
      grid_context$grid_bundle$planner_contract %||% "named_facet_structural_review"
    ),
    input_contract = as.character(
      grid_context$grid_bundle$input_contract %||% "design$facets(named counts)"
    ),
    grid_available = TRUE,
    reason = as.character(
      grid_context$reason %||%
        "Deterministic structural-design grid is materialized."
    ),
    n_designs = nrow(canonical),
    n_varying_axes = sum(varying_rows),
    n_fixed_axes = sum(fixed_rows),
    varying_canonical = as.character(grid_context$varying_canonical),
    fixed_canonical = as.character(grid_context$fixed_canonical),
    varying_input_keys = as.character(grid_context$varying_input_keys),
    fixed_input_keys = as.character(grid_context$fixed_input_keys),
    fixed_values = fixed_values,
    axis_table = axis_table,
    canonical = canonical,
    public = public,
    design = branch,
    grid_context = grid_context,
    grid_bundle = grid_context$grid_bundle,
    note = paste(
      "Deterministic structural-design summary bundling the materialized grid with",
      "axis-level varying/fixed metadata for later arbitrary-facet branch",
      "helpers."
    )
  )
}

simulation_structural_design_axis_lookup <- function(axis_table) {
  axis_table <- tibble::as_tibble(axis_table)
  lookup <- c(
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$canonical_design_variable)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$public_design_alias)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$input_key)
    )
  )
  lookup <- lookup[!is.na(names(lookup)) & nzchar(names(lookup))]
  lookup[!duplicated(names(lookup))]
}

simulation_resolve_structural_design_axis <- function(value,
                                                  axis_table,
                                                  arg_name = "axis",
                                                  allow_null = FALSE) {
  if (allow_null && is.null(value)) {
    return(NULL)
  }
  lookup <- simulation_structural_design_axis_lookup(axis_table)
  value <- as.character(value[1] %||% "")
  resolved <- unname(lookup[[value]])
  if (!is.null(resolved) && nzchar(resolved)) {
    return(resolved)
  }
  stop(
    "`", arg_name, "` must be one of: ",
    paste(sort(unique(names(lookup))), collapse = ", "),
    ".",
    call. = FALSE
  )
}

simulation_structural_design_axis_label <- function(canonical,
                                                axis_table,
                                                view = c("public", "canonical", "design")) {
  view <- match.arg(view)
  axis_table <- tibble::as_tibble(axis_table)
  row <- axis_table[axis_table$canonical_design_variable == canonical, , drop = FALSE]
  if (nrow(row) == 0L) {
    return(as.character(canonical))
  }
  if (identical(view, "canonical")) {
    return(as.character(canonical))
  }
  if (identical(view, "design")) {
    label <- as.character(row$input_key[[1]])
    if (nzchar(label)) return(label)
  }
  label <- as.character(row$public_design_alias[[1]])
  if (nzchar(label)) return(label)
  as.character(canonical)
}

simulation_structural_design_grid_recommendation <- function(x,
                                                         design = NULL,
                                                         prefer = NULL,
                                                         id_prefix = NULL) {
  grid_summary <- simulation_structural_design_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )

  axis_table <- tibble::as_tibble(grid_summary$axis_table)
  lookup <- c(
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$canonical_design_variable)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$public_design_alias)
    ),
    stats::setNames(
      as.character(axis_table$canonical_design_variable),
      as.character(axis_table$input_key)
    )
  )
  lookup <- lookup[!duplicated(names(lookup))]

  if (!isTRUE(grid_summary$grid_available)) {
    return(list(
      recommendation_contract = "arbitrary_facet_design_grid_recommendation",
      recommendation_stage = as.character(grid_summary$summary_stage %||% "deterministic_review"),
      planner_contract = as.character(
        grid_summary$planner_contract %||% "named_facet_structural_review"
      ),
      recommendation_available = FALSE,
      reason = as.character(
        grid_summary$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      selection_rule = paste(
        "Unavailable because no deterministic design grid exists yet.",
        "This helper does not fabricate facet counts."
      ),
      prefer = character(0),
      rank_order = character(0),
      grid_summary = grid_summary,
      note = paste(
        "Deterministic structural-design recommendation is unavailable until the",
        "underlying design grid is materialized."
      )
    ))
  }

  canonical_order <- setdiff(
    as.character(
      grid_summary$grid_bundle$grid_semantics$canonical_columns %||%
        names(grid_summary$canonical)
    ),
    "design_id"
  )
  default_prefer <- if (length(grid_summary$varying_canonical) > 0L) {
    intersect(canonical_order, grid_summary$varying_canonical)
  } else {
    canonical_order
  }

  if (is.null(prefer)) {
    resolved_prefer <- default_prefer
  } else {
    prefer <- unique(as.character(prefer))
    resolved_prefer <- unname(lookup[prefer])
    resolved_prefer <- unique(resolved_prefer[!is.na(resolved_prefer) & nzchar(resolved_prefer)])
    if (length(resolved_prefer) == 0L) {
      stop(
        "`prefer` must resolve to at least one structural-design axis. Valid names: ",
        paste(sort(unique(names(lookup))), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  rank_order <- unique(c(resolved_prefer, setdiff(canonical_order, resolved_prefer)))

  ranked <- dplyr::arrange(grid_summary$canonical, !!!rlang::syms(rank_order))
  recommended_canonical <- dplyr::slice_head(ranked, n = 1)
  recommended_id <- as.character(recommended_canonical$design_id[[1]])
  recommended_public <- dplyr::filter(grid_summary$public, .data$design_id == recommended_id)
  recommended_branch <- dplyr::filter(grid_summary$design, .data$design_id == recommended_id)

  list(
    recommendation_contract = "arbitrary_facet_design_grid_recommendation",
    recommendation_stage = as.character(grid_summary$summary_stage %||% "deterministic_review"),
    planner_contract = as.character(
      grid_summary$planner_contract %||% "named_facet_structural_review"
    ),
    recommendation_available = TRUE,
    reason = paste(
      "Deterministic structural-design baseline design selected by deterministic",
      "lexicographic ordering over the requested branch axes."
    ),
    selection_rule = paste(
      "Pick the lexicographically smallest feasible design row after sorting",
      "the deterministic canonical grid by `prefer`, then by the remaining",
      "canonical design axes. This is a deterministic baseline pick, not a",
      "performance-based recommendation."
    ),
    prefer = resolved_prefer,
    rank_order = rank_order,
    recommended_design_id = recommended_id,
    recommended_canonical = recommended_canonical,
    recommended_public = recommended_public,
    recommended_branch = recommended_branch,
    grid_summary = grid_summary,
    note = paste(
      "Deterministic structural-design baseline pick derived from the currently",
      "materialized design grid without activating arbitrary-facet planner",
      "logic or performance criteria."
    )
  )
}

simulation_structural_design_grid_table <- function(x,
                                                design = NULL,
                                                prefer = NULL,
                                                view = c("public", "canonical", "design"),
                                                id_prefix = NULL) {
  view <- match.arg(view)
  grid_summary <- simulation_structural_design_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  recommendation <- simulation_structural_design_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )

  if (!isTRUE(grid_summary$grid_available)) {
    return(list(
      table_contract = "arbitrary_facet_design_grid_table",
      table_stage = as.character(grid_summary$summary_stage %||% "deterministic_review"),
      planner_contract = as.character(
        grid_summary$planner_contract %||% "named_facet_structural_review"
      ),
      grid_available = FALSE,
      reason = as.character(
        grid_summary$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      view = view,
      view_label = view,
      table = tibble::tibble(),
      recommended_design_id = character(0),
      grid_summary = grid_summary,
      grid_recommendation = recommendation,
      note = paste(
        "Deterministic structural-design table is unavailable until the",
        "underlying design grid is materialized."
      )
    ))
  }

  table <- switch(
    view,
    canonical = tibble::as_tibble(grid_summary$canonical),
    public = tibble::as_tibble(grid_summary$public),
    design = tibble::as_tibble(grid_summary$design)
  )
  table$recommended <- table$design_id == recommendation$recommended_design_id

  list(
    table_contract = "arbitrary_facet_design_grid_table",
    table_stage = as.character(grid_summary$summary_stage %||% "deterministic_review"),
    planner_contract = as.character(
      grid_summary$planner_contract %||% "named_facet_structural_review"
    ),
    grid_available = TRUE,
    reason = as.character(
      grid_summary$reason %||%
        "Deterministic structural-design grid is materialized."
    ),
    view = view,
    view_label = view,
    table = table,
    recommended_design_id = recommendation$recommended_design_id,
    grid_summary = grid_summary,
    grid_recommendation = recommendation,
    note = paste(
      "Deterministic structural-design table exposing one grid view together with",
      "the deterministic baseline design flag."
    )
  )
}

simulation_structural_design_grid_plot_payload <- function(x,
                                                       design = NULL,
                                                       x_var = NULL,
                                                       group_var = NULL,
                                                       prefer = NULL,
                                                       view = c("public", "canonical", "design"),
                                                       id_prefix = NULL) {
  view <- match.arg(view)
  grid_summary <- simulation_structural_design_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  recommendation <- simulation_structural_design_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )
  axis_table <- tibble::as_tibble(grid_summary$axis_table)

  if (!isTRUE(grid_summary$grid_available)) {
    return(new_mfrm_plot_data(
      "structural_design_grid_schema",
      list(
        title = "Deterministic Named-Facet Design Grid",
        subtitle = as.character(
          grid_summary$reason %||%
            "No deterministic design grid is currently materialized."
        ),
        plot_available = FALSE,
        reason = as.character(
          grid_summary$reason %||%
            "No deterministic design grid is currently materialized."
        ),
        view = view,
        x_var = NULL,
        x_label = NULL,
        group_var = NULL,
        group_label = NULL,
        recommended_design_id = character(0),
        data = tibble::tibble(),
        axis_table = axis_table,
        grid_summary = grid_summary,
        grid_recommendation = recommendation,
        note = paste(
          "Draw-free plotting data are unavailable until the named-facet",
          "design grid can be materialized."
        )
      )
    ))
  }

  lookup <- simulation_structural_design_axis_lookup(axis_table)
  varying <- as.character(grid_summary$varying_canonical)
  canonical_default_x <- if (length(varying) > 0L) varying[[1]] else as.character(axis_table$canonical_design_variable[[1]])
  x_canonical <- if (is.null(x_var)) {
    canonical_default_x
  } else {
    simulation_resolve_structural_design_axis(x_var, axis_table, arg_name = "x_var")
  }

  group_default <- setdiff(varying, x_canonical)
  group_canonical <- if (is.null(group_var)) {
    if (length(group_default) > 0L) group_default[[1]] else NULL
  } else {
    simulation_resolve_structural_design_axis(group_var, axis_table, arg_name = "group_var", allow_null = TRUE)
  }
  if (!is.null(group_canonical) && identical(group_canonical, x_canonical)) {
    stop("`group_var` must differ from `x_var`.", call. = FALSE)
  }

  display_table <- switch(
    view,
    canonical = tibble::as_tibble(grid_summary$canonical),
    public = tibble::as_tibble(grid_summary$public),
    design = tibble::as_tibble(grid_summary$design)
  )
  canonical_axes <- tibble::as_tibble(grid_summary$canonical[, c("design_id", unique(c(x_canonical, group_canonical))), drop = FALSE])
  names(canonical_axes)[-1] <- paste0(names(canonical_axes)[-1], "__canonical")
  plot_data <- dplyr::left_join(display_table, canonical_axes, by = "design_id")
  plot_data$x_value <- plot_data[[paste0(x_canonical, "__canonical")]]
  if (!is.null(group_canonical)) {
    plot_data$group_value <- as.character(plot_data[[paste0(group_canonical, "__canonical")]])
  } else {
    plot_data$group_value <- "All designs"
  }
  plot_data$recommended <- plot_data$design_id == recommendation$recommended_design_id
  plot_data <- dplyr::arrange(plot_data, .data$x_value, .data$group_value, dplyr::desc(.data$recommended))

  x_label <- simulation_structural_design_axis_label(x_canonical, axis_table, view = view)
  group_label <- if (is.null(group_canonical)) "Design set" else simulation_structural_design_axis_label(group_canonical, axis_table, view = view)
  fixed_rows <- !is.na(axis_table$varying) & !as.logical(axis_table$varying)
  fixed_text <- if (any(fixed_rows)) {
    fixed_bits <- vapply(which(fixed_rows), function(i) {
      paste0(
        simulation_structural_design_axis_label(axis_table$canonical_design_variable[[i]], axis_table, view = view),
        "=",
        axis_table$fixed_value[[i]]
      )
    }, character(1))
    paste(fixed_bits, collapse = ", ")
  } else {
    "No fixed design axes"
  }

  legend <- if (!is.null(group_canonical)) {
    group_levels <- unique(as.character(plot_data$group_value))
    new_plot_legend(
      label = group_levels,
      role = rep("group", length(group_levels)),
      aesthetic = rep("line", length(group_levels)),
      value = group_levels
    )
  } else {
    new_plot_legend(
      label = "All designs",
      role = "group",
      aesthetic = "line",
      value = "All designs"
    )
  }

  new_mfrm_plot_data(
    "structural_design_grid_schema",
    list(
      title = "Deterministic Named-Facet Design Grid",
      subtitle = fixed_text,
      plot_available = TRUE,
      reason = as.character(
        grid_summary$reason %||%
          "Deterministic structural-design grid is materialized."
      ),
      view = view,
      x_var = x_canonical,
      x_label = x_label,
      group_var = group_canonical,
      group_label = group_label,
      recommended_design_id = recommendation$recommended_design_id,
      selection_rule = recommendation$selection_rule,
      data = plot_data,
      axis_table = axis_table,
      grid_summary = grid_summary,
      grid_recommendation = recommendation,
      legend = legend,
      note = paste(
        "Draw-free plotting data for the named-facet design review. They",
        "summarize the current grid and deterministic baseline selection",
        "without implying support for arbitrary-facet simulation."
      )
    )
  )
}

simulation_structural_design_cached_default <- function(x,
                                                    field,
                                                    design = NULL,
                                                    prefer = NULL,
                                                    view = NULL,
                                                    view_default = NULL,
                                                    mode = NULL,
                                                    mode_default = NULL,
                                                    surface = NULL,
                                                    surface_default = NULL,
                                                    table_component = NULL,
                                                    table_component_default = NULL,
                                                    component = NULL,
                                                    component_default = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  if (!is.list(x) || !is.list(x[[field]])) {
    return(NULL)
  }
  if (!is.null(design) || !is.null(prefer) || !is.null(x_var) ||
      !is.null(group_var) || !is.null(id_prefix)) {
    return(NULL)
  }
  if (!is.null(view_default) &&
      !identical(as.character(view[[1]] %||% view_default), as.character(view_default))) {
    return(NULL)
  }
  if (!is.null(mode_default) &&
      !identical(as.character(mode[[1]] %||% mode_default), as.character(mode_default))) {
    return(NULL)
  }
  if (!is.null(surface_default) &&
      !identical(as.character(surface[[1]] %||% surface_default), as.character(surface_default))) {
    return(NULL)
  }
  if (!is.null(table_component_default) &&
      !identical(
        as.character(table_component[[1]] %||% table_component_default),
        as.character(table_component_default)
      )) {
    return(NULL)
  }
  if (!is.null(component_default) &&
      !identical(as.character(component[[1]] %||% component_default), as.character(component_default))) {
    return(NULL)
  }
  x[[field]]
}

simulation_structural_design_report_bundle <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_bundle",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  grid_summary <- simulation_structural_design_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  recommendation <- simulation_structural_design_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )
  tables <- list(
    canonical = simulation_structural_design_grid_table(
      x = x,
      design = design,
      prefer = prefer,
      view = "canonical",
      id_prefix = id_prefix
    ),
    public = simulation_structural_design_grid_table(
      x = x,
      design = design,
      prefer = prefer,
      view = "public",
      id_prefix = id_prefix
    ),
    design = simulation_structural_design_grid_table(
      x = x,
      design = design,
      prefer = prefer,
      view = "design",
      id_prefix = id_prefix
    )
  )
  plots <- list(
    canonical = simulation_structural_design_grid_plot_payload(
      x = x,
      design = design,
      x_var = x_var,
      group_var = group_var,
      prefer = prefer,
      view = "canonical",
      id_prefix = id_prefix
    ),
    public = simulation_structural_design_grid_plot_payload(
      x = x,
      design = design,
      x_var = x_var,
      group_var = group_var,
      prefer = prefer,
      view = "public",
      id_prefix = id_prefix
    ),
    design = simulation_structural_design_grid_plot_payload(
      x = x,
      design = design,
      x_var = x_var,
      group_var = group_var,
      prefer = prefer,
      view = "design",
      id_prefix = id_prefix
    )
  )

  overview_table <- tibble::as_tibble(grid_summary$axis_table)

  if (!isTRUE(grid_summary$grid_available)) {
    return(list(
      report_contract = "arbitrary_facet_design_grid_report_bundle",
      report_stage = as.character(grid_summary$summary_stage %||% "deterministic_review"),
      planner_contract = as.character(
        grid_summary$planner_contract %||% "named_facet_structural_review"
      ),
      report_available = FALSE,
      reason = as.character(
        grid_summary$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      overview_table = overview_table,
      grid_summary = grid_summary,
      grid_recommendation = recommendation,
      tables = tables,
      plots = plots,
      note = paste(
        "Deterministic structural-design report bundle is unavailable until the",
        "underlying design grid is materialized. The contract still bundles",
        "summary metadata, design-review table contracts, and draw-free plot",
        "data views in one shared object."
      )
    ))
  }

  list(
    report_contract = "arbitrary_facet_design_grid_report_bundle",
    report_stage = as.character(grid_summary$summary_stage %||% "deterministic_review"),
    planner_contract = as.character(
      grid_summary$planner_contract %||% "named_facet_structural_review"
    ),
    report_available = TRUE,
    reason = as.character(
      grid_summary$reason %||%
        "Deterministic structural-design grid is materialized."
    ),
    recommended_design_id = as.character(recommendation$recommended_design_id %||% character(0)),
    overview_table = overview_table,
    grid_summary = grid_summary,
    grid_recommendation = recommendation,
    tables = tables,
    plots = plots,
    note = paste(
      "Deterministic structural-design report bundle combining the shared",
      "summary, deterministic baseline recommendation, canonical/public/",
      "branch table views, and draw-free plot data without implying an",
      "active arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_summary <- function(x,
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_summary",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  report_bundle <- simulation_structural_design_report_bundle(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  table_views <- names(report_bundle$tables %||% list())
  plot_views <- names(report_bundle$plots %||% list())
  table_available <- vapply(
    report_bundle$tables %||% list(),
    function(obj) isTRUE(obj$grid_available),
    logical(1)
  )
  plot_available <- vapply(
    report_bundle$plots %||% list(),
    function(obj) isTRUE(obj$data$plot_available),
    logical(1)
  )
  first_or_na <- function(value) {
    if (length(value) == 0L) {
      return(NA_character_)
    }
    as.character(value[[1]] %||% NA_character_)
  }

  table_index <- tibble::tibble(
    component_type = "table",
    view = table_views,
    available = as.logical(table_available),
    n_rows = vapply(
      report_bundle$tables %||% list(),
      function(obj) if (is.data.frame(obj$table)) nrow(obj$table) else 0L,
      integer(1)
    ),
    x_var = NA_character_,
    group_var = NA_character_,
    recommended_design_id = vapply(
      report_bundle$tables %||% list(),
      function(obj) first_or_na(obj$recommended_design_id),
      character(1)
    )
  )
  plot_index <- tibble::tibble(
    component_type = "plot",
    view = plot_views,
    available = as.logical(plot_available),
    n_rows = vapply(
      report_bundle$plots %||% list(),
      function(obj) if (is.data.frame(obj$data$data)) nrow(obj$data$data) else 0L,
      integer(1)
    ),
    x_var = vapply(
      report_bundle$plots %||% list(),
      function(obj) first_or_na(obj$data$x_var),
      character(1)
    ),
    group_var = vapply(
      report_bundle$plots %||% list(),
      function(obj) first_or_na(obj$data$group_var),
      character(1)
    ),
    recommended_design_id = vapply(
      report_bundle$plots %||% list(),
      function(obj) first_or_na(obj$data$recommended_design_id),
      character(1)
    )
  )
  component_index <- dplyr::bind_rows(
    tibble::tibble(
      component_type = "summary",
      view = NA_character_,
      available = isTRUE(report_bundle$grid_summary$grid_available),
      n_rows = as.integer(report_bundle$grid_summary$n_designs %||% 0L),
      x_var = NA_character_,
      group_var = NA_character_,
      recommended_design_id = NA_character_
    ),
    tibble::tibble(
      component_type = "recommendation",
      view = NA_character_,
      available = isTRUE(report_bundle$grid_recommendation$recommendation_available),
      n_rows = NA_integer_,
      x_var = NA_character_,
      group_var = NA_character_,
      recommended_design_id = first_or_na(report_bundle$grid_recommendation$recommended_design_id)
    ),
    table_index,
    plot_index
  )

  if (!isTRUE(report_bundle$report_available)) {
    return(list(
      report_summary_contract = "arbitrary_facet_design_grid_report_summary",
      report_summary_stage = as.character(report_bundle$report_stage %||% "deterministic_review"),
      planner_contract = as.character(
        report_bundle$planner_contract %||% "named_facet_structural_review"
      ),
      report_available = FALSE,
      reason = as.character(
        report_bundle$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      n_designs = 0L,
      available_table_views = as.character(table_views[table_available]),
      available_plot_views = as.character(plot_views[plot_available]),
      component_index = component_index,
      report_bundle = report_bundle,
      note = paste(
        "Deterministic structural-design report summary is unavailable until the",
        "underlying design grid is materialized. The component index still",
        "describes which design-review report surfaces would be populated."
      )
    ))
  }

  list(
    report_summary_contract = "arbitrary_facet_design_grid_report_summary",
    report_summary_stage = as.character(report_bundle$report_stage %||% "deterministic_review"),
    planner_contract = as.character(
      report_bundle$planner_contract %||% "named_facet_structural_review"
    ),
    report_available = TRUE,
    reason = as.character(
      report_bundle$reason %||%
        "Deterministic structural-design grid is materialized."
    ),
    n_designs = as.integer(report_bundle$grid_summary$n_designs %||% 0L),
    recommended_design_id = first_or_na(report_bundle$recommended_design_id),
    available_table_views = as.character(table_views[table_available]),
    available_plot_views = as.character(plot_views[plot_available]),
    component_index = component_index,
    report_bundle = report_bundle,
    note = paste(
      "Deterministic structural-design report summary exposing a compact component",
      "index over the current report bundle without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_overview_table <- function(x,
                                                           design = NULL,
                                                           prefer = NULL,
                                                           x_var = NULL,
                                                           group_var = NULL,
                                                           id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_overview_table",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  report_summary <- simulation_structural_design_report_summary(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  report_bundle <- report_summary$report_bundle
  overview_source <- tibble::as_tibble(report_bundle$overview_table %||% tibble::tibble())
  has_fixed_value <- "fixed_value" %in% names(overview_source)

  metrics_table <- tibble::tibble(
    report_available = isTRUE(report_summary$report_available),
    n_designs = as.integer(report_summary$n_designs %||% 0L),
    recommended_design_id = as.character(report_summary$recommended_design_id %||% NA_character_),
    n_table_views = length(report_summary$available_table_views %||% character(0)),
    n_plot_views = length(report_summary$available_plot_views %||% character(0)),
    available_table_views = paste(report_summary$available_table_views %||% character(0), collapse = ", "),
    available_plot_views = paste(report_summary$available_plot_views %||% character(0), collapse = ", ")
  )

  axis_overview_table <- if (nrow(overview_source) == 0L) {
    overview_source
  } else {
    overview_source |>
      dplyr::mutate(
        axis_state = dplyr::case_when(
          is.na(.data$varying) ~ "unavailable",
          .data$varying ~ "varying",
          TRUE ~ "fixed"
        ),
        value_summary = vapply(.data$values, function(vals) {
          if (length(vals) == 0L || all(is.na(vals))) {
            return(NA_character_)
          }
          paste(as.character(vals), collapse = ", ")
        }, character(1)),
        fixed_value = if (has_fixed_value) as.integer(.data$fixed_value) else NA_integer_
      ) |>
      dplyr::select(
        dplyr::all_of(c(
          "axis_source",
          "input_key",
          "canonical_design_variable",
          "public_design_alias",
          "axis_class",
          "facet",
          "axis_state",
          "n_values",
          "value_summary",
          "fixed_value"
        ))
      )
  }

  if (!isTRUE(report_summary$report_available)) {
    return(list(
      overview_contract = "arbitrary_facet_design_report_overview_table",
      overview_stage = as.character(report_summary$report_summary_stage %||% "deterministic_review"),
      planner_contract = as.character(
        report_summary$planner_contract %||% "named_facet_structural_review"
      ),
      report_available = FALSE,
      reason = as.character(
        report_summary$reason %||%
          "No deterministic structural-design grid is currently materialized."
      ),
      metrics_table = metrics_table,
      axis_overview_table = axis_overview_table,
      component_index = tibble::as_tibble(report_summary$component_index %||% tibble::tibble()),
      report_summary = report_summary,
      report_bundle = report_bundle,
      note = paste(
        "Deterministic structural-design overview table is unavailable until the",
        "underlying design grid is materialized. Metrics and axis metadata",
        "are still returned from the current report summary contract."
      )
    ))
  }

  list(
    overview_contract = "arbitrary_facet_design_report_overview_table",
    overview_stage = as.character(report_summary$report_summary_stage %||% "deterministic_review"),
    planner_contract = as.character(
      report_summary$planner_contract %||% "named_facet_structural_review"
    ),
    report_available = TRUE,
    reason = as.character(
      report_summary$reason %||%
        "Deterministic structural-design grid is materialized."
    ),
    metrics_table = metrics_table,
    axis_overview_table = axis_overview_table,
    component_index = tibble::as_tibble(report_summary$component_index %||% tibble::tibble()),
    report_summary = report_summary,
    report_bundle = report_bundle,
    note = paste(
      "Deterministic structural-design overview table exposing compact report",
      "metrics and axis-level state summaries without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_overview_view <- function(x,
                                                          component = c("metrics", "axes", "components"),
                                                          design = NULL,
                                                          prefer = NULL,
                                                          x_var = NULL,
                                                          group_var = NULL,
                                                          id_prefix = NULL) {
  component <- match.arg(component)
  overview <- simulation_structural_design_report_overview_table(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  table <- switch(
    component,
    metrics = tibble::as_tibble(overview$metrics_table %||% tibble::tibble()),
    axes = tibble::as_tibble(overview$axis_overview_table %||% tibble::tibble()),
    components = tibble::as_tibble(overview$component_index %||% tibble::tibble())
  )
  component_label <- switch(
    component,
    metrics = "report metrics",
    axes = "axis overview",
    components = "component index"
  )

  list(
    overview_view_contract = "arbitrary_facet_design_report_overview_view",
    overview_stage = as.character(overview$overview_stage %||% "deterministic_review"),
    planner_contract = as.character(
      overview$planner_contract %||% "named_facet_structural_review"
    ),
    component = component,
    component_label = component_label,
    report_available = isTRUE(overview$report_available),
    reason = as.character(
      overview$reason %||%
        "No deterministic structural-design report overview is currently available."
    ),
    table = table,
    report_overview = overview,
    note = paste(
      "Deterministic structural-design overview view exposing the selected",
      component_label,
      "table from the compact report-overview contract without implying an",
      "active arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_catalog <- function(x,
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_catalog",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  component_ids <- c("metrics", "axes", "components")
  views <- stats::setNames(vector("list", length(component_ids)), component_ids)

  for (component in component_ids) {
    views[[component]] <- simulation_structural_design_report_overview_view(
      x = x,
      component = component,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  }

  overview <- views[[1]]$report_overview %||% simulation_structural_design_report_overview_table(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  recommended_design_id <- as.character(
    overview$metrics_table$recommended_design_id[[1]] %||% NA_character_
  )
  surface_index <- tibble::tibble(
    component = component_ids,
    component_label = vapply(views, function(obj) {
      as.character(obj$component_label %||% NA_character_)
    }, character(1)),
    available = vapply(views, function(obj) {
      isTRUE(obj$report_available)
    }, logical(1)),
    n_rows = vapply(views, function(obj) {
      if (is.data.frame(obj$table)) nrow(obj$table) else 0L
    }, integer(1)),
    recommended_design_id = rep(recommended_design_id, length(component_ids))
  )

  list(
    catalog_contract = "arbitrary_facet_design_report_catalog",
    catalog_stage = as.character(overview$overview_stage %||% "deterministic_review"),
    planner_contract = as.character(
      overview$planner_contract %||% "named_facet_structural_review"
    ),
    report_available = isTRUE(overview$report_available),
    reason = as.character(
      overview$reason %||%
        "No deterministic structural-design report catalog is currently available."
    ),
    recommended_design_id = recommended_design_id,
    surface_index = surface_index,
    views = views,
    report_overview = overview,
    note = paste(
      "Deterministic structural-design report catalog enumerating the compact",
      "overview surfaces that design-review reporting code can request from",
      "the current overview contract without implying an active arbitrary-",
      "facet planner."
    )
  )
}

simulation_structural_design_report_digest <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_digest",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  report_catalog <- simulation_structural_design_report_catalog(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  report_overview <- report_catalog$report_overview %||% simulation_structural_design_report_overview_table(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  metrics_table <- tibble::as_tibble(report_overview$metrics_table %||% tibble::tibble())
  axis_table <- tibble::as_tibble(report_overview$axis_overview_table %||% tibble::tibble())
  surface_index <- tibble::as_tibble(report_catalog$surface_index %||% tibble::tibble())

  first_chr <- function(value, default = NA_character_) {
    if (length(value) == 0L) {
      return(as.character(default))
    }
    as.character(value[[1]] %||% default)
  }
  first_int <- function(value, default = 0L) {
    if (length(value) == 0L) {
      return(as.integer(default))
    }
    as.integer(value[[1]] %||% default)
  }

  available_surfaces <- if (nrow(surface_index) == 0L) {
    character(0)
  } else {
    as.character(surface_index$component[replace(as.logical(surface_index$available), is.na(surface_index$available), FALSE)])
  }
  varying_axes <- if (nrow(axis_table) == 0L) {
    character(0)
  } else {
    as.character(axis_table$canonical_design_variable[axis_table$axis_state %in% "varying"])
  }
  fixed_axes <- if (nrow(axis_table) == 0L) {
    character(0)
  } else {
    as.character(axis_table$canonical_design_variable[axis_table$axis_state %in% "fixed"])
  }

  digest_table <- tibble::tibble(
    report_available = isTRUE(report_catalog$report_available),
    n_designs = first_int(metrics_table$n_designs, 0L),
    recommended_design_id = first_chr(metrics_table$recommended_design_id),
    n_available_surfaces = length(available_surfaces),
    available_surfaces = paste(available_surfaces, collapse = ", "),
    varying_axes = paste(varying_axes, collapse = ", "),
    fixed_axes = paste(fixed_axes, collapse = ", ")
  )

  list(
    digest_contract = "arbitrary_facet_design_report_digest",
    digest_stage = as.character(report_catalog$catalog_stage %||% "deterministic_review"),
    planner_contract = as.character(
      report_catalog$planner_contract %||% "named_facet_structural_review"
    ),
    report_available = isTRUE(report_catalog$report_available),
    reason = as.character(
      report_catalog$reason %||%
        "No deterministic structural-design report digest is currently available."
    ),
    recommended_design_id = first_chr(metrics_table$recommended_design_id),
    available_surfaces = available_surfaces,
    varying_axes = varying_axes,
    fixed_axes = fixed_axes,
    digest_table = digest_table,
    report_catalog = report_catalog,
    report_overview = report_overview,
    note = paste(
      "Deterministic structural-design report digest exposing headline report",
      "availability, baseline design metadata, and compact surface/axis",
      "summaries from one shared contract without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_surface <- function(x,
                                                    surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  surface <- match.arg(surface)

  source_object <- switch(
    surface,
    digest = simulation_structural_design_report_digest(
      x = x,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    catalog = simulation_structural_design_report_catalog(
      x = x,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    metrics = simulation_structural_design_report_overview_view(
      x = x,
      component = "metrics",
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    axes = simulation_structural_design_report_overview_view(
      x = x,
      component = "axes",
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    components = simulation_structural_design_report_overview_view(
      x = x,
      component = "components",
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  )

  surface_label <- switch(
    surface,
    digest = "report digest",
    catalog = "report catalog",
    metrics = "report metrics",
    axes = "axis overview",
    components = "component index"
  )
  source_contract <- switch(
    surface,
    digest = as.character(source_object$digest_contract %||% "arbitrary_facet_design_report_digest"),
    catalog = as.character(source_object$catalog_contract %||% "arbitrary_facet_design_report_catalog"),
    metrics = as.character(source_object$overview_view_contract %||% "arbitrary_facet_design_report_overview_view"),
    axes = as.character(source_object$overview_view_contract %||% "arbitrary_facet_design_report_overview_view"),
    components = as.character(source_object$overview_view_contract %||% "arbitrary_facet_design_report_overview_view")
  )
  recommended_design_id <- switch(
    surface,
    digest = as.character(source_object$recommended_design_id %||% NA_character_),
    catalog = as.character(source_object$recommended_design_id %||% NA_character_),
    metrics = as.character(source_object$table$recommended_design_id[[1]] %||% NA_character_),
    axes = as.character(source_object$report_overview$metrics_table$recommended_design_id[[1]] %||% NA_character_),
    components = as.character(source_object$report_overview$metrics_table$recommended_design_id[[1]] %||% NA_character_)
  )
  table <- switch(
    surface,
    digest = tibble::as_tibble(source_object$digest_table %||% tibble::tibble()),
    catalog = {
      tbl <- tibble::as_tibble(source_object$surface_index %||% tibble::tibble())
      if ("component" %in% names(tbl) && !("surface" %in% names(tbl))) {
        names(tbl)[names(tbl) == "component"] <- "surface"
      }
      tbl
    },
    metrics = tibble::as_tibble(source_object$table %||% tibble::tibble()),
    axes = tibble::as_tibble(source_object$table %||% tibble::tibble()),
    components = tibble::as_tibble(source_object$table %||% tibble::tibble())
  )
  report_available <- switch(
    surface,
    digest = isTRUE(source_object$report_available),
    catalog = isTRUE(source_object$report_available),
    metrics = isTRUE(source_object$report_available),
    axes = isTRUE(source_object$report_available),
    components = isTRUE(source_object$report_available)
  )
  reason <- switch(
    surface,
    digest = as.character(source_object$reason %||% NA_character_),
    catalog = as.character(source_object$reason %||% NA_character_),
    metrics = as.character(source_object$reason %||% NA_character_),
    axes = as.character(source_object$reason %||% NA_character_),
    components = as.character(source_object$reason %||% NA_character_)
  )
  surface_stage <- switch(
    surface,
    digest = as.character(source_object$digest_stage %||% "deterministic_review"),
    catalog = as.character(source_object$catalog_stage %||% "deterministic_review"),
    metrics = as.character(source_object$overview_stage %||% "deterministic_review"),
    axes = as.character(source_object$overview_stage %||% "deterministic_review"),
    components = as.character(source_object$overview_stage %||% "deterministic_review")
  )
  planner_contract <- switch(
    surface,
    digest = as.character(source_object$planner_contract %||% "named_facet_structural_review"),
    catalog = as.character(source_object$planner_contract %||% "named_facet_structural_review"),
    metrics = as.character(source_object$planner_contract %||% "named_facet_structural_review"),
    axes = as.character(source_object$planner_contract %||% "named_facet_structural_review"),
    components = as.character(source_object$planner_contract %||% "named_facet_structural_review")
  )

  list(
    surface_contract = "arbitrary_facet_design_report_surface",
    surface_stage = surface_stage,
    planner_contract = planner_contract,
    surface = surface,
    surface_label = surface_label,
    source_contract = source_contract,
    report_available = report_available,
    reason = reason,
    recommended_design_id = recommended_design_id,
    table = table,
    source_object = source_object,
    note = paste(
      "Deterministic structural-design report surface exposing the selected",
      surface_label,
      "from one compact design-review operation without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_surface_registry <- function(x,
                                                             design = NULL,
                                                             prefer = NULL,
                                                             x_var = NULL,
                                                             group_var = NULL,
                                                             id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_surface_registry",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface_ids <- c("digest", "catalog", "metrics", "axes", "components")
  surfaces <- stats::setNames(vector("list", length(surface_ids)), surface_ids)

  for (surface in surface_ids) {
    surfaces[[surface]] <- simulation_structural_design_report_surface(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  }

  recommended_design_id <- as.character(
    surfaces$digest$recommended_design_id %||%
      surfaces$catalog$recommended_design_id %||%
      NA_character_
  )
  surface_index <- tibble::tibble(
    surface = surface_ids,
    surface_label = vapply(surfaces, function(obj) {
      as.character(obj$surface_label %||% NA_character_)
    }, character(1)),
    source_contract = vapply(surfaces, function(obj) {
      as.character(obj$source_contract %||% NA_character_)
    }, character(1)),
    available = vapply(surfaces, function(obj) {
      isTRUE(obj$report_available)
    }, logical(1)),
    n_rows = vapply(surfaces, function(obj) {
      if (is.data.frame(obj$table)) nrow(obj$table) else 0L
    }, integer(1)),
    recommended_design_id = rep(recommended_design_id, length(surface_ids))
  )

  list(
    registry_contract = "arbitrary_facet_design_report_surface_registry",
    registry_stage = as.character(surfaces[[1]]$surface_stage %||% "deterministic_review"),
    planner_contract = as.character(
      surfaces[[1]]$planner_contract %||% "named_facet_structural_review"
    ),
    report_available = isTRUE(surfaces[[1]]$report_available),
    reason = as.character(
      surfaces[[1]]$reason %||%
        "No deterministic structural-design report surfaces are currently available."
    ),
    recommended_design_id = recommended_design_id,
    surface_index = surface_index,
    surfaces = surfaces,
    note = paste(
      "Deterministic structural-design report surface registry exposing digest,",
      "catalog, and compact overview surfaces from one shared contract",
      "without implying an active arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_panel <- function(x,
                                                  surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                  design = NULL,
                                                  prefer = NULL,
                                                  x_var = NULL,
                                                  group_var = NULL,
                                                  id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_panel",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  registry <- simulation_structural_design_report_surface_registry(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  selected_surface <- registry$surfaces[[surface]] %||% simulation_structural_design_report_surface(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  digest_table <- tibble::as_tibble(
    registry$surfaces$digest$table %||% tibble::tibble()
  )

  list(
    panel_contract = "arbitrary_facet_design_report_panel",
    panel_stage = as.character(registry$registry_stage %||% "deterministic_review"),
    planner_contract = as.character(
      registry$planner_contract %||% "named_facet_structural_review"
    ),
    surface = surface,
    surface_label = as.character(selected_surface$surface_label %||% surface),
    report_available = isTRUE(selected_surface$report_available),
    reason = as.character(
      selected_surface$reason %||%
        "No deterministic structural-design report panel is currently available."
    ),
    recommended_design_id = as.character(
      selected_surface$recommended_design_id %||%
        registry$recommended_design_id %||%
        NA_character_
    ),
    digest_table = digest_table,
    surface_index = tibble::as_tibble(registry$surface_index %||% tibble::tibble()),
    selected_table = tibble::as_tibble(selected_surface$table %||% tibble::tibble()),
    selected_surface = selected_surface,
    registry = registry,
    note = paste(
      "Deterministic structural-design report panel exposing one selected compact",
      "surface together with headline digest metadata and the surface index",
      "from one shared contract without implying an active arbitrary-",
      "facet planner."
    )
  )
}

simulation_structural_design_report_operation <- function(x,
                                                      surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                      design = NULL,
                                                      prefer = NULL,
                                                      x_var = NULL,
                                                      group_var = NULL,
                                                      id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_operation",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  panel <- simulation_structural_design_report_panel(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  registry <- panel$registry %||% simulation_structural_design_report_surface_registry(
    x = x,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  digest_surface <- registry$surfaces$digest %||% simulation_structural_design_report_surface(
    x = x,
    surface = "digest",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  metrics_surface <- registry$surfaces$metrics %||% simulation_structural_design_report_surface(
    x = x,
    surface = "metrics",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  axes_surface <- registry$surfaces$axes %||% simulation_structural_design_report_surface(
    x = x,
    surface = "axes",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  components_surface <- registry$surfaces$components %||% simulation_structural_design_report_surface(
    x = x,
    surface = "components",
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  list(
    operation_contract = "arbitrary_facet_design_report_operation",
    operation_stage = as.character(
      panel$panel_stage %||% registry$registry_stage %||% "deterministic_review"
    ),
    planner_contract = as.character(
      panel$planner_contract %||%
        registry$planner_contract %||%
        "named_facet_structural_review"
    ),
    surface = surface,
    surface_label = as.character(panel$surface_label %||% surface),
    report_available = isTRUE(panel$report_available),
    reason = as.character(
      panel$reason %||%
        "No deterministic structural-design report operation is currently available."
    ),
    recommended_design_id = as.character(
      panel$recommended_design_id %||%
        registry$recommended_design_id %||%
        NA_character_
    ),
    digest_table = tibble::as_tibble(digest_surface$table %||% tibble::tibble()),
    surface_index = tibble::as_tibble(registry$surface_index %||% tibble::tibble()),
    metrics_table = tibble::as_tibble(metrics_surface$table %||% tibble::tibble()),
    axis_overview_table = tibble::as_tibble(axes_surface$table %||% tibble::tibble()),
    component_index = tibble::as_tibble(components_surface$table %||% tibble::tibble()),
    selected_table = tibble::as_tibble(panel$selected_table %||% tibble::tibble()),
    selected_surface = panel$selected_surface,
    report_panel = panel,
    report_surface_registry = registry,
    note = paste(
      "Deterministic structural-design report operation exposing digest, compact",
      "overview tables, the current surface registry, and one selected",
      "surface from one shared contract without implying an active",
      "arbitrary-facet planner."
    )
  )
}

simulation_structural_design_report_snapshot <- function(x,
                                                     surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                     design = NULL,
                                                     prefer = NULL,
                                                     x_var = NULL,
                                                     group_var = NULL,
                                                     id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_snapshot",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  operation <- simulation_structural_design_report_operation(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  first_chr <- function(value, default = NA_character_) {
    if (length(value) == 0L) {
      return(as.character(default))
    }
    as.character(value[[1]] %||% default)
  }

  digest_tbl <- tibble::as_tibble(operation$digest_table %||% tibble::tibble())
  surface_index <- tibble::as_tibble(operation$surface_index %||% tibble::tibble())
  selected_tbl <- tibble::as_tibble(operation$selected_table %||% tibble::tibble())

  available_surfaces <- if (nrow(surface_index) == 0L) {
    character(0)
  } else {
    as.character(surface_index$surface[replace(as.logical(surface_index$available), is.na(surface_index$available), FALSE)])
  }
  varying_axes <- if (!("varying_axes" %in% names(digest_tbl))) {
    character(0)
  } else {
    strsplit(first_chr(digest_tbl$varying_axes, ""), ", ", fixed = TRUE)[[1]]
  }
  varying_axes <- varying_axes[nzchar(varying_axes)]
  fixed_axes <- if (!("fixed_axes" %in% names(digest_tbl))) {
    character(0)
  } else {
    strsplit(first_chr(digest_tbl$fixed_axes, ""), ", ", fixed = TRUE)[[1]]
  }
  fixed_axes <- fixed_axes[nzchar(fixed_axes)]

  list(
    snapshot_contract = "arbitrary_facet_design_report_snapshot",
    snapshot_stage = as.character(operation$operation_stage %||% "deterministic_review"),
    planner_contract = as.character(
      operation$planner_contract %||% "named_facet_structural_review"
    ),
    surface = surface,
    surface_label = as.character(operation$surface_label %||% surface),
    report_available = isTRUE(operation$report_available),
    reason = as.character(
      operation$reason %||%
        "No deterministic structural-design report snapshot is currently available."
    ),
    recommended_design_id = as.character(
      operation$recommended_design_id %||% NA_character_
    ),
    n_designs = if ("n_designs" %in% names(digest_tbl)) as.integer(digest_tbl$n_designs[[1]] %||% 0L) else 0L,
    available_surfaces = available_surfaces,
    varying_axes = varying_axes,
    fixed_axes = fixed_axes,
    digest_table = digest_tbl,
    surface_index = surface_index,
    selected_table = selected_tbl,
    metrics_table = tibble::as_tibble(operation$metrics_table %||% tibble::tibble()),
    axis_overview_table = tibble::as_tibble(operation$axis_overview_table %||% tibble::tibble()),
    component_index = tibble::as_tibble(operation$component_index %||% tibble::tibble()),
    note = paste(
      "Deterministic structural-design report snapshot exposing compact headline",
      "metadata, compact overview tables, and one selected surface without",
      "carrying the deeper nested operation graph."
    )
  )
}

simulation_structural_design_report_brief <- function(x,
                                                  surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                  design = NULL,
                                                  prefer = NULL,
                                                  x_var = NULL,
                                                  group_var = NULL,
                                                  id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_brief",
    design = design,
    prefer = prefer,
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  surface <- match.arg(surface)
  snapshot <- simulation_structural_design_report_snapshot(
    x = x,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  headline_table <- tibble::tibble(
    report_available = isTRUE(snapshot$report_available),
    n_designs = as.integer(snapshot$n_designs %||% 0L),
    recommended_design_id = as.character(
      snapshot$recommended_design_id %||% NA_character_
    ),
    surface = as.character(snapshot$surface %||% surface),
    surface_label = as.character(snapshot$surface_label %||% surface),
    available_surfaces = paste(snapshot$available_surfaces %||% character(0), collapse = ", "),
    varying_axes = paste(snapshot$varying_axes %||% character(0), collapse = ", "),
    fixed_axes = paste(snapshot$fixed_axes %||% character(0), collapse = ", ")
  )

  list(
    brief_contract = "arbitrary_facet_design_report_brief",
    brief_stage = as.character(snapshot$snapshot_stage %||% "deterministic_review"),
    planner_contract = as.character(
      snapshot$planner_contract %||% "named_facet_structural_review"
    ),
    surface = surface,
    surface_label = as.character(snapshot$surface_label %||% surface),
    report_available = isTRUE(snapshot$report_available),
    reason = as.character(
      snapshot$reason %||%
        "No deterministic structural-design report brief is currently available."
    ),
    headline_table = headline_table,
    selected_table = tibble::as_tibble(snapshot$selected_table %||% tibble::tibble()),
    surface_index = tibble::as_tibble(snapshot$surface_index %||% tibble::tibble()),
    note = paste(
      "Deterministic structural-design report brief exposing one selected surface",
      "together with a headline table and surface index, without carrying",
      "the broader snapshot or nested operation graph."
    )
  )
}

simulation_structural_design_report_consume <- function(x,
                                                    mode = c("brief", "snapshot", "operation"),
                                                    surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                    design = NULL,
                                                    prefer = NULL,
                                                    x_var = NULL,
                                                    group_var = NULL,
                                                    id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "report_consumer",
    design = design,
    prefer = prefer,
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  mode <- match.arg(mode)
  surface <- match.arg(surface)

  payload <- switch(
    mode,
    brief = simulation_structural_design_report_brief(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    snapshot = simulation_structural_design_report_snapshot(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    ),
    operation = simulation_structural_design_report_operation(
      x = x,
      surface = surface,
      design = design,
      prefer = prefer,
      x_var = x_var,
      group_var = group_var,
      id_prefix = id_prefix
    )
  )

  payload_contract <- switch(
    mode,
    brief = as.character(payload$brief_contract %||% NA_character_),
    snapshot = as.character(payload$snapshot_contract %||% NA_character_),
    operation = as.character(payload$operation_contract %||% NA_character_)
  )
  stage <- switch(
    mode,
    brief = as.character(payload$brief_stage %||% "deterministic_review"),
    snapshot = as.character(payload$snapshot_stage %||% "deterministic_review"),
    operation = as.character(payload$operation_stage %||% "deterministic_review")
  )
  selected_table <- switch(
    mode,
    brief = tibble::as_tibble(payload$selected_table %||% tibble::tibble()),
    snapshot = tibble::as_tibble(payload$selected_table %||% tibble::tibble()),
    operation = tibble::as_tibble(payload$selected_table %||% tibble::tibble())
  )
  surface_index <- switch(
    mode,
    brief = tibble::as_tibble(payload$surface_index %||% tibble::tibble()),
    snapshot = tibble::as_tibble(payload$surface_index %||% tibble::tibble()),
    operation = tibble::as_tibble(payload$surface_index %||% tibble::tibble())
  )
  recommended_design_id <- switch(
    mode,
    brief = as.character(payload$headline_table$recommended_design_id[[1]] %||% NA_character_),
    snapshot = as.character(payload$recommended_design_id %||% NA_character_),
    operation = as.character(payload$recommended_design_id %||% NA_character_)
  )
  n_designs <- switch(
    mode,
    brief = as.integer(payload$headline_table$n_designs[[1]] %||% 0L),
    snapshot = as.integer(payload$n_designs %||% 0L),
    operation = if ("n_designs" %in% names(payload$digest_table)) {
      as.integer(payload$digest_table$n_designs[[1]] %||% 0L)
    } else {
      0L
    }
  )

  list(
    consumer_contract = "arbitrary_facet_design_report_consumer",
    consumer_stage = stage,
    planner_contract = as.character(
      payload$planner_contract %||% "named_facet_structural_review"
    ),
    mode = mode,
    payload_contract = payload_contract,
    surface = surface,
    surface_label = as.character(payload$surface_label %||% surface),
    report_available = isTRUE(payload$report_available),
    reason = as.character(
      payload$reason %||%
        "No deterministic structural-design report data are currently available."
    ),
    recommended_design_id = recommended_design_id,
    n_designs = n_designs,
    selected_table = selected_table,
    surface_index = surface_index,
    payload = payload,
    note = paste(
      "Deterministic structural-design report consumer dispatching to the selected",
      mode,
      "contract so design-review code can switch report-data weight without",
      "re-implementing the dispatch logic."
    )
  )
}

simulation_structural_design_report_mode_registry <- function(x) {
  cached <- if (is.list(x)) x$report_mode_registry %||% NULL else NULL
  if (is.list(cached) || (is.data.frame(cached) && nrow(cached) >= 0L)) {
    return(cached)
  }

  tibble::tibble(
    mode = c("brief", "snapshot", "operation"),
    payload_contract = c(
      "arbitrary_facet_design_report_brief",
      "arbitrary_facet_design_report_snapshot",
      "arbitrary_facet_design_report_operation"
    ),
    default_surface = "digest",
    carries_nested_graph = c(FALSE, FALSE, TRUE),
    note = c(
      "Selected surface plus headline table and surface index.",
      "Compact headline metadata plus one selected surface.",
      "Full compact report operation with selected surface and nested registry."
    )
  )
}

simulation_structural_design_compact_review <- function(x,
                                           design = NULL,
                                           prefer = NULL,
                                           view = c("public", "canonical", "design"),
                                           mode = c("brief", "snapshot", "operation"),
                                           surface = c("digest", "catalog", "metrics", "axes", "components"),
                                           x_var = NULL,
                                           group_var = NULL,
                                           id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "review",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  view <- match.arg(view)
  mode <- match.arg(mode)
  surface <- match.arg(surface)

  grid_summary <- simulation_structural_design_grid_summary(
    x = x,
    design = design,
    id_prefix = id_prefix
  )
  grid_recommendation <- simulation_structural_design_grid_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    id_prefix = id_prefix
  )
  grid_table <- simulation_structural_design_grid_table(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    id_prefix = id_prefix
  )
  plot_payload <- simulation_structural_design_grid_plot_payload(
    x = x,
    design = design,
    x_var = x_var,
    group_var = group_var,
    prefer = prefer,
    view = view,
    id_prefix = id_prefix
  )
  report_consumer <- simulation_structural_design_report_consume(
    x = x,
    mode = mode,
    surface = surface,
    design = design,
    prefer = prefer,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  review_available <- isTRUE(grid_summary$grid_available) && isTRUE(report_consumer$report_available)

  list(
    compact_review_contract = "named_facet_compact_review",
    review_stage = if (review_available) "review_ready" else "deterministic_review",
    planner_contract = as.character(
      grid_summary$planner_contract %||%
        report_consumer$planner_contract %||%
        "named_facet_structural_review"
    ),
    review_available = review_available,
    reason = as.character(
      if (review_available) {
        "Deterministic structural-design review is materialized from the current grid and report consumer contracts."
      } else {
        grid_summary$reason %||%
          report_consumer$reason %||%
          "No deterministic structural-design review is currently available."
      }
    ),
    view = view,
    mode = mode,
    surface = surface,
    recommended_design_id = as.character(
      report_consumer$recommended_design_id %||%
        grid_recommendation$recommended_design_id %||%
        NA_character_
    ),
    n_designs = as.integer(
      report_consumer$n_designs %||%
        grid_summary$n_designs %||%
        0L
    ),
    grid_table = tibble::as_tibble(grid_table$table %||% tibble::tibble()),
    plot_payload = plot_payload,
    report_consumer = report_consumer,
    grid_summary = grid_summary,
    grid_recommendation = grid_recommendation,
    note = paste(
      "Internal deterministic structural-design review bundling one materialized grid",
      "view, one draw-free plot-data object, and one mode-selected report consumer.",
      "This is a compact structural-review object, but it remains a",
      "deterministic review rather than a performance-based arbitrary-facet planner."
    )
  )
}

simulation_structural_design_compact_review_summary <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   view = c("public", "canonical", "design"),
                                                   mode = c("brief", "snapshot", "operation"),
                                                   surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "compact_review_summary",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  review <- simulation_structural_design_compact_review(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  headline_table <- tibble::tibble(
    review_available = isTRUE(review$review_available),
    n_designs = as.integer(review$n_designs %||% 0L),
    recommended_design_id = as.character(review$recommended_design_id %||% NA_character_),
    view = as.character(review$view %||% NA_character_),
    mode = as.character(review$mode %||% NA_character_),
    surface = as.character(review$surface %||% NA_character_),
    grid_rows = if (is.data.frame(review$grid_table)) nrow(review$grid_table) else 0L,
    plot_available = isTRUE(review$plot_payload$data$plot_available %||% FALSE),
    report_payload_contract = as.character(
      review$report_consumer$payload_contract %||% NA_character_
    )
  )

  list(
    compact_review_summary_contract = "arbitrary_facet_planning_compact_review_summary",
    review_stage = as.character(review$review_stage %||% "deterministic_review"),
    planner_contract = as.character(
      review$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(review$review_available),
    reason = as.character(
      review$reason %||%
        "No deterministic structural-design review summary is currently available."
    ),
    headline_table = headline_table,
    review = review,
    note = paste(
      "Compact summary for the deterministic structural-design review, exposing one",
      "headline table over the current grid and report layer without implying",
      "an active performance-based arbitrary-facet planner."
    )
  )
}

simulation_structural_design_compact_review_table <- function(x,
                                                 component = c("grid", "report", "surface_index"),
                                                 design = NULL,
                                                 prefer = NULL,
                                                 view = c("public", "canonical", "design"),
                                                 mode = c("brief", "snapshot", "operation"),
                                                 surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                 x_var = NULL,
                                                 group_var = NULL,
                                                 id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "compact_review_table",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    component = component,
    component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  component <- match.arg(component)
  review <- simulation_structural_design_compact_review(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  table <- switch(
    component,
    grid = tibble::as_tibble(review$grid_table %||% tibble::tibble()),
    report = tibble::as_tibble(review$report_consumer$selected_table %||% tibble::tibble()),
    surface_index = tibble::as_tibble(review$report_consumer$surface_index %||% tibble::tibble())
  )

  list(
    compact_review_table_contract = "arbitrary_facet_planning_compact_review_table",
    review_stage = as.character(review$review_stage %||% "deterministic_review"),
    planner_contract = as.character(
      review$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(review$review_available),
    reason = as.character(
      review$reason %||%
        "No deterministic structural-design review table is currently available."
    ),
    component = component,
    table = table,
    review = review,
    note = paste(
      "Selected table view from the deterministic structural-design review.",
      "This exposes either the grid table, the selected report table, or",
      "the compact surface index from the current review layer."
    )
  )
}

simulation_structural_design_compact_review_plot <- function(x,
                                                design = NULL,
                                                prefer = NULL,
                                                view = c("public", "canonical", "design"),
                                                mode = c("brief", "snapshot", "operation"),
                                                surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                x_var = NULL,
                                                group_var = NULL,
                                                id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "compact_review_plot",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  review <- simulation_structural_design_compact_review(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  list(
    compact_review_plot_contract = "arbitrary_facet_planning_compact_review_plot",
    review_stage = as.character(review$review_stage %||% "deterministic_review"),
    planner_contract = as.character(
      review$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(review$review_available),
    reason = as.character(
      review$reason %||%
        "No deterministic structural-design review plot is currently available."
    ),
    plot = review$plot_payload,
    review = review,
    note = paste(
      "Draw-free plot data from the deterministic structural-design review.",
      "This preserves the existing plot contract while routing access through",
      "the review layer."
    )
  )
}

simulation_structural_design_review <- function(x,
                                                   design = NULL,
                                                   prefer = NULL,
                                                   view = c("public", "canonical", "design"),
                                                   mode = c("brief", "snapshot", "operation"),
                                                   surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                   table_component = c("grid", "report", "surface_index"),
                                                   x_var = NULL,
                                                   group_var = NULL,
                                                   id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  if (inherits(x, "mfrm_structural_design_review") ||
      identical(as.character(x$review_contract %||% NA_character_), "arbitrary_facet_planning_design_review")) {
    return(x)
  }

  view <- match.arg(view)
  mode <- match.arg(mode)
  surface <- match.arg(surface)
  table_component <- match.arg(table_component)

  compact_review_summary <- simulation_structural_design_compact_review_summary(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  compact_review_table <- simulation_structural_design_compact_review_table(
    x = x,
    component = table_component,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  compact_review_plot <- simulation_structural_design_compact_review_plot(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  review_available <- isTRUE(compact_review_summary$review_available) &&
    isTRUE(compact_review_table$review_available) &&
    isTRUE(compact_review_plot$review_available)
  canonical_grid <- tibble::as_tibble(
    compact_review_summary$review$grid_summary$canonical %||% tibble::tibble()
  )

  structure(list(
    review_contract = "arbitrary_facet_planning_design_review",
    review_stage = if (review_available) "review_ready" else "deterministic_review",
    planner_contract = as.character(
      compact_review_summary$planner_contract %||%
        compact_review_table$planner_contract %||%
        compact_review_plot$planner_contract %||%
        "named_facet_structural_review"
    ),
    review_available = review_available,
    reason = as.character(
      if (review_available) {
        "Deterministic structural-design review is materialized from the review summary, table, and plot contracts."
      } else {
        compact_review_summary$reason %||%
          compact_review_table$reason %||%
          compact_review_plot$reason %||%
          "No deterministic structural-design review is currently available."
      }
    ),
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    recommended_design_id = as.character(
      compact_review_summary$headline_table$recommended_design_id[[1]] %||% NA_character_
    ),
    n_designs = as.integer(
      compact_review_summary$headline_table$n_designs[[1]] %||% 0L
    ),
    summary = compact_review_summary,
    table = compact_review_table,
    plot = compact_review_plot,
    canonical_grid = canonical_grid,
    note = paste(
      "Compact structural-design review object bundling review-level summary, one",
      "selected table component, and one draw-free plot-data object from the",
      "current named-facet review layer."
    )
  ), class = c("mfrm_structural_design_review", "list"))
}

simulation_structural_design_review_canonical_grid <- function(branch) {
  candidates <- list(
    branch$canonical_grid,
    branch$summary$review$grid_summary$canonical,
    branch$table$review$grid_summary$canonical,
    branch$plot$review$grid_summary$canonical
  )

  for (candidate in candidates) {
    if (is.data.frame(candidate) && nrow(candidate) > 0L) {
      return(tibble::as_tibble(candidate))
    }
  }

  tibble::tibble()
}

simulation_structural_design_review_metric_registry <- function() {
  tibble::tibble(
    metric = c(
      "total_observations",
      "observations_per_person",
      "observations_per_criterion",
      "expected_observations_per_rater",
      "assignment_fraction"
    ),
    basis_class = c(
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "balanced_expectation",
      "density_ratio"
    ),
    formula = c(
      "n_person * raters_per_person * n_criterion",
      "raters_per_person * n_criterion",
      "n_person * raters_per_person",
      "(n_person * raters_per_person * n_criterion) / n_rater",
      "raters_per_person / n_rater"
    ),
    interpretation = c(
      "Total number of scored observations implied by the current design row.",
      "Number of scored observations contributed by one person under the current design row.",
      "Number of scored observations contributed to one criterion across persons.",
      "Average expected rater load under balanced assignment, not a guaranteed realized count.",
      "Fraction of available raters assigned to each person."
    ),
    psychometric = FALSE
  )
}

simulation_structural_design_review_load_balance_registry <- function() {
  tibble::tibble(
    metric = c(
      "expected_observations_per_rater",
      "expected_person_assignments_per_rater",
      "observation_load_floor",
      "observation_load_ceiling",
      "observation_load_remainder",
      "perfect_integer_observation_balance"
    ),
    basis_class = c(
      "balanced_expectation",
      "balanced_expectation",
      "integer_balance_bound",
      "integer_balance_bound",
      "exact_identity",
      "exact_identity"
    ),
    formula = c(
      "(n_person * raters_per_person * n_criterion) / n_rater",
      "(n_person * raters_per_person) / n_rater",
      "floor((n_person * raters_per_person * n_criterion) / n_rater)",
      "ceiling((n_person * raters_per_person * n_criterion) / n_rater)",
      "(n_person * raters_per_person * n_criterion) %% n_rater",
      "as.integer(((n_person * raters_per_person * n_criterion) %% n_rater) == 0)"
    ),
    interpretation = c(
      "Average observation load per rater under balanced assignment, not a guaranteed realized count.",
      "Average person-assignment load per rater under balanced assignment, not a guaranteed realized count.",
      "Lower integer observation count per rater under the most even observation-level split implied by the current design row.",
      "Upper integer observation count per rater under the most even observation-level split implied by the current design row.",
      "Number of observations left over after equal integer splitting across raters.",
      "Indicator that the total observation count is divisible by the number of raters."
    ),
    psychometric = FALSE
  )
}

simulation_structural_design_review_coverage_registry <- function() {
  tibble::tibble(
    metric = c(
      "person_criterion_cells",
      "criterion_replications_per_person",
      "rater_pair_overlap_per_cell",
      "total_rater_pair_overlaps",
      "pair_coverage_fraction_per_cell",
      "redundant_scoring"
    ),
    basis_class = c(
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "exact_ratio",
      "exact_identity"
    ),
    formula = c(
      "n_person * n_criterion",
      "raters_per_person",
      "choose(raters_per_person, 2)",
      "n_person * n_criterion * choose(raters_per_person, 2)",
      "ifelse(choose(n_rater, 2) > 0, choose(raters_per_person, 2) / choose(n_rater, 2), 0)",
      "as.integer(raters_per_person > 1)"
    ),
    interpretation = c(
      "Number of person-by-criterion cells implied by the current design row.",
      "Number of ratings attached to each person-by-criterion cell.",
      "Number of distinct rater pairs overlapping within one scored cell.",
      "Total number of rater-pair overlaps implied across all person-by-criterion cells.",
      "Fraction of available rater pairs represented within one scored cell; set to 0 when fewer than two raters are available.",
      "Indicator that the design uses more than one rater per scored cell."
    ),
    psychometric = FALSE
  )
}

simulation_structural_design_review_profile <- function(x,
                                                           design = NULL,
                                                           prefer = NULL,
                                                           view = c("public", "canonical", "design"),
                                                           mode = c("brief", "snapshot", "operation"),
                                                           surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                           table_component = c("grid", "report", "surface_index"),
                                                           x_var = NULL,
                                                           group_var = NULL,
                                                           id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_profile",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  branch <- simulation_structural_design_review(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  canonical <- simulation_structural_design_review_canonical_grid(branch)

  if (!isTRUE(branch$review_available) || nrow(canonical) == 0L) {
    return(list(
      profile_contract = "arbitrary_facet_planning_design_review_profile",
      profile_stage = as.character(branch$review_stage %||% "deterministic_review"),
      planner_contract = as.character(
        branch$planner_contract %||% "named_facet_structural_review"
      ),
      review_available = FALSE,
      reason = as.character(
        branch$reason %||%
          "No deterministic structural-design profile is currently available."
      ),
      metric_registry = simulation_structural_design_review_metric_registry(),
      profile_summary_table = tibble::tibble(),
      profile_table = tibble::tibble(),
      design_review = branch,
      note = paste(
        "Deterministic design-profile metrics are unavailable until the",
        "structural-design review layer can be materialized."
      )
    ))
  }

  profile_table <- canonical |>
    dplyr::mutate(
      total_observations = as.numeric(.data$n_person) * as.numeric(.data$raters_per_person) * as.numeric(.data$n_criterion),
      observations_per_person = as.numeric(.data$raters_per_person) * as.numeric(.data$n_criterion),
      observations_per_criterion = as.numeric(.data$n_person) * as.numeric(.data$raters_per_person),
      expected_observations_per_rater = .data$total_observations / as.numeric(.data$n_rater),
      assignment_fraction = as.numeric(.data$raters_per_person) / as.numeric(.data$n_rater),
      recommended = .data$design_id == as.character(branch$recommended_design_id %||% NA_character_)
    )
  metric_registry <- simulation_structural_design_review_metric_registry()
  recommended_row <- profile_table[profile_table$recommended %in% TRUE, , drop = FALSE]
  profile_summary_table <- metric_registry |>
    dplyr::mutate(
      min = vapply(.data$metric, function(metric) {
        min(profile_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$metric, function(metric) {
        max(profile_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$metric, function(metric) {
        mean(profile_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$metric, function(metric) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[metric]][[1]])
      }, numeric(1))
    )

  list(
    profile_contract = "arbitrary_facet_planning_design_review_profile",
    profile_stage = as.character(branch$review_stage %||% "deterministic_review"),
    planner_contract = as.character(
      branch$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = TRUE,
    reason = as.character(
      branch$reason %||%
        "Deterministic structural-design profile is materialized."
    ),
    metric_registry = metric_registry,
    profile_summary_table = profile_summary_table,
    profile_table = profile_table,
    design_review = branch,
    note = paste(
      "Deterministic design-bookkeeping profile for the structural-design review",
      "review layer. These quantities summarize observation counts, expected rater",
      "load, and assignment density; they are not psychometric performance estimates."
    )
  )
}

simulation_structural_design_review_load_balance <- function(x,
                                                                design = NULL,
                                                                prefer = NULL,
                                                                view = c("public", "canonical", "design"),
                                                                mode = c("brief", "snapshot", "operation"),
                                                                surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                table_component = c("grid", "report", "surface_index"),
                                                                x_var = NULL,
                                                                group_var = NULL,
                                                                id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_load_balance",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  profile <- simulation_structural_design_review_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  profile_table <- tibble::as_tibble(profile$profile_table %||% tibble::tibble())
  metric_registry <- simulation_structural_design_review_load_balance_registry()

  if (!isTRUE(profile$review_available) || nrow(profile_table) == 0L) {
    return(list(
      diagnostics_contract = "arbitrary_facet_planning_design_review_load_balance",
      diagnostics_stage = as.character(profile$profile_stage %||% "deterministic_review"),
      planner_contract = as.character(
        profile$planner_contract %||% "named_facet_structural_review"
      ),
      review_available = FALSE,
      reason = as.character(
        profile$reason %||%
          "No deterministic structural-design load/balance diagnostics are currently available."
      ),
      metric_registry = metric_registry,
      diagnostic_summary_table = tibble::tibble(),
      diagnostic_table = tibble::tibble(),
      design_review_profile = profile,
      note = paste(
        "Deterministic load/balance diagnostics are unavailable until the",
        "structural-design review layer can be materialized."
      )
    ))
  }

  diagnostic_table <- profile_table |>
    dplyr::mutate(
      expected_person_assignments_per_rater = (as.numeric(.data$n_person) * as.numeric(.data$raters_per_person)) / as.numeric(.data$n_rater),
      observation_load_floor = floor(as.numeric(.data$total_observations) / as.numeric(.data$n_rater)),
      observation_load_ceiling = ceiling(as.numeric(.data$total_observations) / as.numeric(.data$n_rater)),
      observation_load_remainder = as.numeric(.data$total_observations) %% as.numeric(.data$n_rater),
      perfect_integer_observation_balance = as.integer(.data$observation_load_remainder == 0)
    )
  recommended_row <- diagnostic_table[diagnostic_table$recommended %in% TRUE, , drop = FALSE]
  diagnostic_summary_table <- metric_registry |>
    dplyr::mutate(
      min = vapply(.data$metric, function(metric) {
        min(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$metric, function(metric) {
        max(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$metric, function(metric) {
        mean(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$metric, function(metric) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[metric]][[1]])
      }, numeric(1))
    )

  list(
    diagnostics_contract = "arbitrary_facet_planning_design_review_load_balance",
    diagnostics_stage = as.character(profile$profile_stage %||% "deterministic_review"),
    planner_contract = as.character(
      profile$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = TRUE,
    reason = as.character(
      profile$reason %||%
        "Deterministic structural-design load/balance diagnostics are materialized."
    ),
    metric_registry = metric_registry,
    diagnostic_summary_table = diagnostic_summary_table,
    diagnostic_table = diagnostic_table,
    design_review_profile = profile,
    note = paste(
      "Deterministic load/balance diagnostics for the structural-design review",
      "review layer. Balanced-load quantities are labeled as expectations, while",
      "integer split diagnostics are exact combinatorial summaries of the",
      "current observation counts."
    )
  )
}

simulation_structural_design_review_overview <- function(x,
                                                            design = NULL,
                                                            prefer = NULL,
                                                            view = c("public", "canonical", "design"),
                                                            mode = c("brief", "snapshot", "operation"),
                                                            surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                            table_component = c("grid", "report", "surface_index"),
                                                            x_var = NULL,
                                                            group_var = NULL,
                                                            id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  branch <- simulation_structural_design_review(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  profile <- simulation_structural_design_review_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  headline_table <- tibble::tibble(
    review_available = isTRUE(branch$review_available),
    n_designs = as.integer(branch$n_designs %||% 0L),
    recommended_design_id = as.character(branch$recommended_design_id %||% NA_character_),
    view = as.character(branch$view %||% NA_character_),
    mode = as.character(branch$mode %||% NA_character_),
    surface = as.character(branch$surface %||% NA_character_),
    table_component = as.character(branch$table_component %||% NA_character_),
    selected_table_rows = if (is.data.frame(branch$table$table)) nrow(branch$table$table) else 0L,
    plot_available = isTRUE(branch$plot$plot$data$plot_available %||% FALSE),
    n_metrics = if (is.data.frame(profile$metric_registry)) nrow(profile$metric_registry) else 0L
  )

  list(
    overview_contract = "arbitrary_facet_planning_design_review_overview",
    overview_stage = as.character(branch$review_stage %||% "deterministic_review"),
    planner_contract = as.character(
      branch$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(branch$review_available) && isTRUE(profile$review_available),
    reason = as.character(
      branch$reason %||%
        profile$reason %||%
        "No structural-design overview is currently available."
    ),
    headline_table = headline_table,
    metric_registry = tibble::as_tibble(profile$metric_registry %||% tibble::tibble()),
    metric_summary_table = tibble::as_tibble(profile$profile_summary_table %||% tibble::tibble()),
    design_review = branch,
    design_review_profile = profile,
    note = paste(
      "Compact overview for the structural-design review layer, combining the",
      "branch headline with metric-basis-aware deterministic design summaries."
    )
  )
}

simulation_structural_design_review_load_balance_overview <- function(x,
                                                                         design = NULL,
                                                                         prefer = NULL,
                                                                         view = c("public", "canonical", "design"),
                                                                         mode = c("brief", "snapshot", "operation"),
                                                                         surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                         table_component = c("grid", "report", "surface_index"),
                                                                         x_var = NULL,
                                                                         group_var = NULL,
                                                                         id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_load_balance_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  diagnostics <- simulation_structural_design_review_load_balance(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  diagnostic_table <- tibble::as_tibble(diagnostics$diagnostic_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    review_available = isTRUE(diagnostics$review_available),
    n_designs = nrow(diagnostic_table),
    n_metrics = if (is.data.frame(diagnostics$metric_registry)) nrow(diagnostics$metric_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(diagnostic_table) &&
                                any(diagnostic_table$recommended %in% TRUE)) {
      as.character(diagnostic_table$design_id[which(diagnostic_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_perfect_integer_balance = if ("perfect_integer_observation_balance" %in% names(diagnostic_table)) {
      sum(diagnostic_table$perfect_integer_observation_balance %in% 1L, na.rm = TRUE)
    } else {
      0L
    },
    n_nondivisible_designs = if ("perfect_integer_observation_balance" %in% names(diagnostic_table)) {
      sum(diagnostic_table$perfect_integer_observation_balance %in% 0L, na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_design_review_load_balance_overview",
    overview_stage = as.character(diagnostics$diagnostics_stage %||% "deterministic_review"),
    planner_contract = as.character(
      diagnostics$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(diagnostics$review_available),
    reason = as.character(
      diagnostics$reason %||%
        "No structural-design load/balance overview is currently available."
    ),
    headline_table = headline_table,
    metric_registry = tibble::as_tibble(diagnostics$metric_registry %||% tibble::tibble()),
    diagnostic_summary_table = tibble::as_tibble(diagnostics$diagnostic_summary_table %||% tibble::tibble()),
    design_review_load_balance = diagnostics,
    note = paste(
      "Compact load/balance overview for the structural-design review layer,",
      "combining a one-row divisibility headline with basis-aware deterministic",
      "observation-load diagnostics."
    )
  )
}

simulation_structural_design_review_coverage <- function(x,
                                                            design = NULL,
                                                            prefer = NULL,
                                                            view = c("public", "canonical", "design"),
                                                            mode = c("brief", "snapshot", "operation"),
                                                            surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                            table_component = c("grid", "report", "surface_index"),
                                                            x_var = NULL,
                                                            group_var = NULL,
                                                            id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_coverage",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  profile <- simulation_structural_design_review_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  profile_table <- tibble::as_tibble(profile$profile_table %||% tibble::tibble())
  metric_registry <- simulation_structural_design_review_coverage_registry()

  if (!isTRUE(profile$review_available) || nrow(profile_table) == 0L) {
    return(list(
      diagnostics_contract = "arbitrary_facet_planning_design_review_coverage",
      diagnostics_stage = as.character(profile$profile_stage %||% "deterministic_review"),
      planner_contract = as.character(
        profile$planner_contract %||% "named_facet_structural_review"
      ),
      review_available = FALSE,
      reason = as.character(
        profile$reason %||%
          "No deterministic structural-design coverage/connectivity diagnostics are currently available."
      ),
      metric_registry = metric_registry,
      diagnostic_summary_table = tibble::tibble(),
      diagnostic_table = tibble::tibble(),
      design_review_profile = profile,
      note = paste(
        "Deterministic coverage/connectivity diagnostics are unavailable until",
        "the structural-design review layer can be materialized."
      )
    ))
  }

  diagnostic_table <- profile_table |>
    dplyr::mutate(
      person_criterion_cells = as.numeric(.data$n_person) * as.numeric(.data$n_criterion),
      criterion_replications_per_person = as.numeric(.data$raters_per_person),
      available_rater_pairs = choose(as.numeric(.data$n_rater), 2),
      rater_pair_overlap_per_cell = choose(as.numeric(.data$raters_per_person), 2),
      total_rater_pair_overlaps = as.numeric(.data$person_criterion_cells) * as.numeric(.data$rater_pair_overlap_per_cell),
      pair_coverage_fraction_per_cell = dplyr::if_else(
        .data$available_rater_pairs > 0,
        as.numeric(.data$rater_pair_overlap_per_cell) / as.numeric(.data$available_rater_pairs),
        0
      ),
      redundant_scoring = as.integer(as.numeric(.data$raters_per_person) > 1)
    ) |>
    dplyr::select(-"available_rater_pairs")

  recommended_row <- diagnostic_table[diagnostic_table$recommended %in% TRUE, , drop = FALSE]
  diagnostic_summary_table <- metric_registry |>
    dplyr::mutate(
      min = vapply(.data$metric, function(metric) {
        min(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$metric, function(metric) {
        max(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$metric, function(metric) {
        mean(diagnostic_table[[metric]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$metric, function(metric) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[metric]][[1]])
      }, numeric(1))
    )

  list(
    diagnostics_contract = "arbitrary_facet_planning_design_review_coverage",
    diagnostics_stage = as.character(profile$profile_stage %||% "deterministic_review"),
    planner_contract = as.character(
      profile$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = TRUE,
    reason = as.character(
      profile$reason %||%
        "Deterministic structural-design coverage/connectivity diagnostics are materialized."
    ),
    metric_registry = metric_registry,
    diagnostic_summary_table = diagnostic_summary_table,
    diagnostic_table = diagnostic_table,
    design_review_profile = profile,
    note = paste(
      "Deterministic coverage/connectivity diagnostics for the structural-design review",
      "review layer. These quantities summarize scored-cell counts and rater-pair",
      "overlap identities without implying psychometric performance."
    )
  )
}

simulation_structural_design_review_coverage_overview <- function(x,
                                                                     design = NULL,
                                                                     prefer = NULL,
                                                                     view = c("public", "canonical", "design"),
                                                                     mode = c("brief", "snapshot", "operation"),
                                                                     surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                     table_component = c("grid", "report", "surface_index"),
                                                                     x_var = NULL,
                                                                     group_var = NULL,
                                                                     id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_coverage_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  diagnostics <- simulation_structural_design_review_coverage(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  diagnostic_table <- tibble::as_tibble(diagnostics$diagnostic_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    review_available = isTRUE(diagnostics$review_available),
    n_designs = nrow(diagnostic_table),
    n_metrics = if (is.data.frame(diagnostics$metric_registry)) nrow(diagnostics$metric_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(diagnostic_table) &&
                                any(diagnostic_table$recommended %in% TRUE)) {
      as.character(diagnostic_table$design_id[which(diagnostic_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_redundant_designs = if ("redundant_scoring" %in% names(diagnostic_table)) {
      sum(diagnostic_table$redundant_scoring %in% 1L, na.rm = TRUE)
    } else {
      0L
    },
    n_single_rater_designs = if ("redundant_scoring" %in% names(diagnostic_table)) {
      sum(diagnostic_table$redundant_scoring %in% 0L, na.rm = TRUE)
    } else {
      0L
    },
    n_pair_connected_designs = if ("rater_pair_overlap_per_cell" %in% names(diagnostic_table)) {
      sum(diagnostic_table$rater_pair_overlap_per_cell > 0, na.rm = TRUE)
    } else {
      0L
    },
    n_zero_pair_overlap_designs = if ("rater_pair_overlap_per_cell" %in% names(diagnostic_table)) {
      sum(diagnostic_table$rater_pair_overlap_per_cell <= 0, na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_design_review_coverage_overview",
    overview_stage = as.character(diagnostics$diagnostics_stage %||% "deterministic_review"),
    planner_contract = as.character(
      diagnostics$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(diagnostics$review_available),
    reason = as.character(
      diagnostics$reason %||%
        "No structural-design coverage/connectivity overview is currently available."
    ),
    headline_table = headline_table,
    metric_registry = tibble::as_tibble(diagnostics$metric_registry %||% tibble::tibble()),
    diagnostic_summary_table = tibble::as_tibble(diagnostics$diagnostic_summary_table %||% tibble::tibble()),
    design_review_coverage = diagnostics,
    note = paste(
      "Compact coverage/connectivity overview for the structural-design review",
      "review layer, combining a one-row redundancy headline with basis-aware",
      "deterministic rater-overlap diagnostics."
    )
  )
}

simulation_structural_design_review_guardrail_registry <- function() {
  tibble::tibble(
    guardrail = c(
      "rater_linking_regime",
      "pair_coverage_regime",
      "integer_balance_regime",
      "redundancy_regime"
    ),
    basis_class = c(
      "exact_classification",
      "exact_classification",
      "exact_classification",
      "exact_classification"
    ),
    rule = c(
      "ifelse(raters_per_person <= 1, 'single_rater', ifelse(raters_per_person >= n_rater, 'fully_crossed', 'partial_overlap'))",
      "ifelse(pair_coverage_fraction_per_cell <= 0, 'no_pair_coverage', ifelse(pair_coverage_fraction_per_cell >= 1, 'full_pair_coverage', 'partial_pair_coverage'))",
      "ifelse(perfect_integer_observation_balance == 1, 'integer_balanced', 'integer_unbalanced')",
      "ifelse(redundant_scoring == 1, 'redundant', 'single_rater_only')"
    ),
    interpretation = c(
      "Exact cell-level linking regime implied by the current assignments-per-person relative to the available rater count.",
      "Exact pair-coverage regime implied by the current within-cell rater overlap fraction.",
      "Exact integer-balance regime induced by the current total observation count and rater count.",
      "Exact indicator of whether each scored cell uses more than one rater."
    ),
    psychometric = FALSE
  )
}

simulation_structural_design_review_guardrails <- function(x,
                                                              design = NULL,
                                                              prefer = NULL,
                                                              view = c("public", "canonical", "design"),
                                                              mode = c("brief", "snapshot", "operation"),
                                                              surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                              table_component = c("grid", "report", "surface_index"),
                                                              x_var = NULL,
                                                              group_var = NULL,
                                                              id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_guardrails",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  coverage <- simulation_structural_design_review_coverage(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  load_balance <- simulation_structural_design_review_load_balance(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  coverage_table <- tibble::as_tibble(coverage$diagnostic_table %||% tibble::tibble())
  load_balance_table <- tibble::as_tibble(load_balance$diagnostic_table %||% tibble::tibble())
  guardrail_registry <- simulation_structural_design_review_guardrail_registry()

  if (!isTRUE(coverage$review_available) || nrow(coverage_table) == 0L ||
      !isTRUE(load_balance$review_available) || nrow(load_balance_table) == 0L) {
    return(list(
      guardrail_contract = "arbitrary_facet_planning_design_review_guardrails",
      guardrail_stage = as.character(
        coverage$diagnostics_stage %||% load_balance$diagnostics_stage %||% "deterministic_review"
      ),
      planner_contract = as.character(
        coverage$planner_contract %||%
          load_balance$planner_contract %||%
          "named_facet_structural_review"
      ),
      review_available = FALSE,
      reason = as.character(
        coverage$reason %||%
          load_balance$reason %||%
          "No deterministic structural-design guardrail classifications are currently available."
      ),
      guardrail_registry = guardrail_registry,
      guardrail_summary_table = tibble::tibble(),
      guardrail_table = tibble::tibble(),
      design_review_coverage = coverage,
      design_review_load_balance = load_balance,
      note = paste(
        "Deterministic guardrail classifications are unavailable until the",
        "structural-design review layer can be materialized."
      )
    ))
  }

  join_cols <- c("design_id", "recommended")
  joined <- dplyr::left_join(
    coverage_table,
    load_balance_table[, c(join_cols, "perfect_integer_observation_balance"), drop = FALSE],
    by = join_cols
  )

  guardrail_table <- joined |>
    dplyr::mutate(
      rater_linking_regime = dplyr::case_when(
        as.numeric(.data$criterion_replications_per_person) <= 1 ~ "single_rater",
        as.numeric(.data$criterion_replications_per_person) >= as.numeric(.data$n_rater) ~ "fully_crossed",
        TRUE ~ "partial_overlap"
      ),
      pair_coverage_regime = dplyr::case_when(
        as.numeric(.data$pair_coverage_fraction_per_cell) <= 0 ~ "no_pair_coverage",
        as.numeric(.data$pair_coverage_fraction_per_cell) >= 1 ~ "full_pair_coverage",
        TRUE ~ "partial_pair_coverage"
      ),
      integer_balance_regime = dplyr::if_else(
        as.integer(.data$perfect_integer_observation_balance) == 1L,
        "integer_balanced",
        "integer_unbalanced"
      ),
      redundancy_regime = dplyr::if_else(
        as.integer(.data$redundant_scoring) == 1L,
        "redundant",
        "single_rater_only"
      )
    )

  guardrail_summary_table <- guardrail_registry |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_levels = dplyr::n_distinct(guardrail_table[[.data$guardrail]]),
      recommended_level = {
        row <- guardrail_table[guardrail_table$recommended %in% TRUE, , drop = FALSE]
        if (nrow(row) == 0L) NA_character_ else as.character(row[[.data$guardrail]][[1]])
      },
      observed_levels = paste(sort(unique(as.character(guardrail_table[[.data$guardrail]]))), collapse = ", ")
    ) |>
    dplyr::ungroup()

  list(
    guardrail_contract = "arbitrary_facet_planning_design_review_guardrails",
    guardrail_stage = as.character(
      coverage$diagnostics_stage %||% load_balance$diagnostics_stage %||% "deterministic_review"
    ),
    planner_contract = as.character(
      coverage$planner_contract %||%
        load_balance$planner_contract %||%
        "named_facet_structural_review"
    ),
    review_available = TRUE,
    reason = "Deterministic structural-design guardrail classifications are materialized.",
    guardrail_registry = guardrail_registry,
    guardrail_summary_table = guardrail_summary_table,
    guardrail_table = guardrail_table,
    design_review_coverage = coverage,
    design_review_load_balance = load_balance,
    note = paste(
      "Deterministic guardrail classifications for the structural-design review",
      "review layer. These are exact design-regime labels derived from overlap and",
      "integer-balance structure; they do not estimate psychometric performance."
    )
  )
}

simulation_structural_design_review_guardrail_overview <- function(x,
                                                                      design = NULL,
                                                                      prefer = NULL,
                                                                      view = c("public", "canonical", "design"),
                                                                      mode = c("brief", "snapshot", "operation"),
                                                                      surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                      table_component = c("grid", "report", "surface_index"),
                                                                      x_var = NULL,
                                                                      group_var = NULL,
                                                                      id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_guardrail_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  guardrails <- simulation_structural_design_review_guardrails(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  guardrail_table <- tibble::as_tibble(guardrails$guardrail_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    review_available = isTRUE(guardrails$review_available),
    n_designs = nrow(guardrail_table),
    n_guardrails = if (is.data.frame(guardrails$guardrail_registry)) nrow(guardrails$guardrail_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(guardrail_table) &&
                                any(guardrail_table$recommended %in% TRUE)) {
      as.character(guardrail_table$design_id[which(guardrail_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_single_rater_designs = if ("rater_linking_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$rater_linking_regime %in% "single_rater", na.rm = TRUE)
    } else {
      0L
    },
    n_partial_overlap_designs = if ("rater_linking_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$rater_linking_regime %in% "partial_overlap", na.rm = TRUE)
    } else {
      0L
    },
    n_fully_crossed_designs = if ("rater_linking_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$rater_linking_regime %in% "fully_crossed", na.rm = TRUE)
    } else {
      0L
    },
    n_integer_balanced_designs = if ("integer_balance_regime" %in% names(guardrail_table)) {
      sum(guardrail_table$integer_balance_regime %in% "integer_balanced", na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_design_review_guardrail_overview",
    overview_stage = as.character(guardrails$guardrail_stage %||% "deterministic_review"),
    planner_contract = as.character(
      guardrails$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(guardrails$review_available),
    reason = as.character(
      guardrails$reason %||%
        "No structural-design guardrail overview is currently available."
    ),
    headline_table = headline_table,
    guardrail_registry = tibble::as_tibble(guardrails$guardrail_registry %||% tibble::tibble()),
    guardrail_summary_table = tibble::as_tibble(guardrails$guardrail_summary_table %||% tibble::tibble()),
    design_review_guardrails = guardrails,
    note = paste(
      "Compact guardrail overview for the structural-design review layer,",
      "combining exact design-regime counts with a basis-aware guardrail index."
    )
  )
}

simulation_structural_design_review_readiness_registry <- function() {
  tibble::tibble(
    indicator = c(
      "supports_multi_rater_cells",
      "supports_pair_overlap",
      "supports_integer_balanced_load",
      "supports_full_pair_coverage"
    ),
    basis_class = c(
      "exact_indicator",
      "exact_indicator",
      "exact_indicator",
      "exact_indicator"
    ),
    rule = c(
      "as.integer(redundancy_regime == 'redundant')",
      "as.integer(pair_coverage_regime != 'no_pair_coverage')",
      "as.integer(integer_balance_regime == 'integer_balanced')",
      "as.integer(pair_coverage_regime == 'full_pair_coverage')"
    ),
    interpretation = c(
      "Exact indicator that each scored cell uses more than one rater.",
      "Exact indicator that the current design yields at least one within-cell rater pair overlap.",
      "Exact indicator that the total observation count is evenly divisible across raters.",
      "Exact indicator that every available rater pair is represented within each scored cell."
    ),
    psychometric = FALSE
  )
}

simulation_structural_design_review_readiness <- function(x,
                                                             design = NULL,
                                                             prefer = NULL,
                                                             view = c("public", "canonical", "design"),
                                                             mode = c("brief", "snapshot", "operation"),
                                                             surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                             table_component = c("grid", "report", "surface_index"),
                                                             x_var = NULL,
                                                             group_var = NULL,
                                                             id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_readiness",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  guardrails <- simulation_structural_design_review_guardrails(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  guardrail_table <- tibble::as_tibble(guardrails$guardrail_table %||% tibble::tibble())
  indicator_registry <- simulation_structural_design_review_readiness_registry()

  if (!isTRUE(guardrails$review_available) || nrow(guardrail_table) == 0L) {
    return(list(
      readiness_contract = "arbitrary_facet_planning_design_review_readiness",
      readiness_stage = as.character(guardrails$guardrail_stage %||% "deterministic_review"),
      planner_contract = as.character(
        guardrails$planner_contract %||% "named_facet_structural_review"
      ),
      review_available = FALSE,
      reason = as.character(
        guardrails$reason %||%
          "No deterministic structural-design structural readiness summary is currently available."
      ),
      indicator_registry = indicator_registry,
      readiness_summary_table = tibble::tibble(),
      readiness_table = tibble::tibble(),
      design_review_guardrails = guardrails,
      note = paste(
        "Deterministic structural readiness indicators are unavailable until",
        "the structural-design review layer can be materialized."
      )
    ))
  }

  readiness_table <- guardrail_table |>
    dplyr::mutate(
      supports_multi_rater_cells = as.integer(.data$redundancy_regime == "redundant"),
      supports_pair_overlap = as.integer(.data$pair_coverage_regime != "no_pair_coverage"),
      supports_integer_balanced_load = as.integer(.data$integer_balance_regime == "integer_balanced"),
      supports_full_pair_coverage = as.integer(.data$pair_coverage_regime == "full_pair_coverage"),
      structural_tier = dplyr::case_when(
        .data$supports_pair_overlap <= 0 ~ "single_rater_only",
        .data$supports_full_pair_coverage >= 1 & .data$supports_integer_balanced_load >= 1 ~ "full_overlap_balanced",
        .data$supports_full_pair_coverage >= 1 ~ "full_overlap_unbalanced",
        .data$supports_integer_balanced_load >= 1 ~ "partial_overlap_balanced",
        TRUE ~ "partial_overlap_unbalanced"
      )
    )

  recommended_row <- readiness_table[readiness_table$recommended %in% TRUE, , drop = FALSE]
  readiness_summary_table <- indicator_registry |>
    dplyr::mutate(
      min = vapply(.data$indicator, function(indicator) {
        min(readiness_table[[indicator]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(.data$indicator, function(indicator) {
        max(readiness_table[[indicator]], na.rm = TRUE)
      }, numeric(1)),
      mean = vapply(.data$indicator, function(indicator) {
        mean(readiness_table[[indicator]], na.rm = TRUE)
      }, numeric(1)),
      recommended_value = vapply(.data$indicator, function(indicator) {
        if (nrow(recommended_row) == 0L) return(NA_real_)
        as.numeric(recommended_row[[indicator]][[1]])
      }, numeric(1))
    )

  list(
    readiness_contract = "arbitrary_facet_planning_design_review_readiness",
    readiness_stage = as.character(guardrails$guardrail_stage %||% "deterministic_review"),
    planner_contract = as.character(
      guardrails$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = TRUE,
    reason = "Deterministic structural-design structural readiness summary is materialized.",
    indicator_registry = indicator_registry,
    readiness_summary_table = readiness_summary_table,
    readiness_table = readiness_table,
    design_review_guardrails = guardrails,
    note = paste(
      "Deterministic structural readiness indicators for the structural-design review",
      "review layer. These summarize exact overlap and balance preconditions only;",
      "they do not estimate psychometric performance."
    )
  )
}

simulation_structural_design_review_readiness_overview <- function(x,
                                                                      design = NULL,
                                                                      prefer = NULL,
                                                                      view = c("public", "canonical", "design"),
                                                                      mode = c("brief", "snapshot", "operation"),
                                                                      surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                      table_component = c("grid", "report", "surface_index"),
                                                                      x_var = NULL,
                                                                      group_var = NULL,
                                                                      id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_readiness_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  readiness <- simulation_structural_design_review_readiness(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  readiness_table <- tibble::as_tibble(readiness$readiness_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    review_available = isTRUE(readiness$review_available),
    n_designs = nrow(readiness_table),
    n_indicators = if (is.data.frame(readiness$indicator_registry)) nrow(readiness$indicator_registry) else 0L,
    recommended_design_id = if ("design_id" %in% names(readiness_table) &&
                                any(readiness_table$recommended %in% TRUE)) {
      as.character(readiness_table$design_id[which(readiness_table$recommended %in% TRUE)[1]])
    } else {
      NA_character_
    },
    n_single_rater_only_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% "single_rater_only", na.rm = TRUE)
    } else {
      0L
    },
    n_partial_overlap_balanced_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% "partial_overlap_balanced", na.rm = TRUE)
    } else {
      0L
    },
    n_partial_overlap_unbalanced_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% "partial_overlap_unbalanced", na.rm = TRUE)
    } else {
      0L
    },
    n_full_overlap_tiers = if ("structural_tier" %in% names(readiness_table)) {
      sum(readiness_table$structural_tier %in% c("full_overlap_balanced", "full_overlap_unbalanced"), na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_design_review_readiness_overview",
    overview_stage = as.character(readiness$readiness_stage %||% "deterministic_review"),
    planner_contract = as.character(
      readiness$planner_contract %||% "named_facet_structural_review"
    ),
    review_available = isTRUE(readiness$review_available),
    reason = as.character(
      readiness$reason %||%
        "No structural-design structural readiness overview is currently available."
    ),
    headline_table = headline_table,
    indicator_registry = tibble::as_tibble(readiness$indicator_registry %||% tibble::tibble()),
    readiness_summary_table = tibble::as_tibble(readiness$readiness_summary_table %||% tibble::tibble()),
    design_review_readiness = readiness,
    note = paste(
      "Compact structural readiness overview for the structural-design review",
      "review layer, combining exact overlap/balance tiers with a basis-aware",
      "indicator index."
    )
  )
}

simulation_structural_design_review_recommendation <- function(x,
                                                                  design = NULL,
                                                                  prefer = NULL,
                                                                  view = c("public", "canonical", "design"),
                                                                  mode = c("brief", "snapshot", "operation"),
                                                                  surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                  table_component = c("grid", "report", "surface_index"),
                                                                  x_var = NULL,
                                                                  group_var = NULL,
                                                                  id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_recommendation",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  readiness <- simulation_structural_design_review_readiness(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  profile <- simulation_structural_design_review_profile(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  readiness_table <- tibble::as_tibble(readiness$readiness_table %||% tibble::tibble())
  profile_table <- tibble::as_tibble(profile$profile_table %||% tibble::tibble())

  if (!isTRUE(readiness$review_available) || nrow(readiness_table) == 0L ||
      !isTRUE(profile$review_available) || nrow(profile_table) == 0L) {
    return(list(
      recommendation_contract = "arbitrary_facet_planning_design_review_recommendation",
      recommendation_stage = as.character(
        readiness$readiness_stage %||% profile$profile_stage %||% "deterministic_review"
      ),
      planner_contract = as.character(
        readiness$planner_contract %||%
          profile$planner_contract %||%
          "named_facet_structural_review"
      ),
      recommendation_available = FALSE,
      reason = as.character(
        readiness$reason %||%
          profile$reason %||%
          "No deterministic structural-design structural recommendation is currently available."
      ),
      ranking_rule = paste(
        "Unavailable because the structural-design readiness review and profile",
        "contracts are not both materialized."
      ),
      recommended_design_id = character(0),
      recommended_tier = character(0),
      recommendation_table = tibble::tibble(),
      design_review_readiness = readiness,
      design_review_profile = profile,
      note = paste(
        "Structural recommendation is unavailable until the structural-design review",
        "review layer can be materialized."
      )
    ))
  }

  tier_priority <- c(
    full_overlap_balanced = 1L,
    full_overlap_unbalanced = 2L,
    partial_overlap_balanced = 3L,
    partial_overlap_unbalanced = 4L,
    single_rater_only = 5L
  )

  recommendation_table <- readiness_table
  required_profile_cols <- c(
    "total_observations", "n_person", "n_rater", "n_criterion", "raters_per_person"
  )
  missing_profile_cols <- setdiff(required_profile_cols, names(recommendation_table))
  if (length(missing_profile_cols) > 0L) {
    recommendation_table <- dplyr::left_join(
      recommendation_table,
      profile_table[, c("design_id", missing_profile_cols), drop = FALSE],
      by = "design_id"
    )
  }

  recommendation_table <- recommendation_table |>
    dplyr::mutate(
      structural_priority = unname(tier_priority[as.character(.data$structural_tier)]),
      structural_priority = dplyr::coalesce(.data$structural_priority, length(tier_priority) + 1L)
    ) |>
    dplyr::arrange(
      .data$structural_priority,
      .data$total_observations,
      .data$n_person,
      .data$n_criterion,
      .data$raters_per_person,
      .data$n_rater,
      .data$design_id
    ) |>
    dplyr::mutate(recommended = dplyr::row_number() == 1L)

  recommended_row <- dplyr::slice_head(recommendation_table, n = 1)
  recommended_id <- as.character(recommended_row$design_id[[1]] %||% NA_character_)
  recommended_tier <- as.character(recommended_row$structural_tier[[1]] %||% NA_character_)

  list(
    recommendation_contract = "arbitrary_facet_planning_design_review_recommendation",
    recommendation_stage = as.character(
      readiness$readiness_stage %||% profile$profile_stage %||% "deterministic_review"
    ),
    planner_contract = as.character(
      readiness$planner_contract %||%
        profile$planner_contract %||%
        "named_facet_structural_review"
    ),
    recommendation_available = TRUE,
    reason = paste(
      "Structural recommendation is available from the structural-design review",
      "readiness and profile contracts."
    ),
    ranking_rule = paste(
      "Rank by structural tier in the order",
      "full_overlap_balanced > full_overlap_unbalanced > partial_overlap_balanced >",
      "partial_overlap_unbalanced > single_rater_only, then minimize total_observations",
      "and remaining canonical count variables."
    ),
    recommended_design_id = recommended_id,
    recommended_tier = recommended_tier,
    recommendation_table = recommendation_table,
    design_review_readiness = readiness,
    design_review_profile = profile,
    note = paste(
      "Conservative structural recommendation for the structural-design review",
      "review layer. This ranking uses exact overlap/balance tiers and total",
      "observation counts only; it is not a psychometric optimization."
    )
  )
}

simulation_structural_design_review_recommendation_overview <- function(x,
                                                                           design = NULL,
                                                                           prefer = NULL,
                                                                           view = c("public", "canonical", "design"),
                                                                           mode = c("brief", "snapshot", "operation"),
                                                                           surface = c("digest", "catalog", "metrics", "axes", "components"),
                                                                           table_component = c("grid", "report", "surface_index"),
                                                                           x_var = NULL,
                                                                           group_var = NULL,
                                                                           id_prefix = NULL) {
  cached <- simulation_structural_design_cached_default(
    x = x,
    field = "design_review_recommendation_overview",
    design = design,
    prefer = prefer,
    view = view,
    view_default = "public",
    mode = mode,
    mode_default = "brief",
    surface = surface,
    surface_default = "digest",
    table_component = table_component,
    table_component_default = "grid",
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )
  if (!is.null(cached)) {
    return(cached)
  }

  recommendation <- simulation_structural_design_review_recommendation(
    x = x,
    design = design,
    prefer = prefer,
    view = view,
    mode = mode,
    surface = surface,
    table_component = table_component,
    x_var = x_var,
    group_var = group_var,
    id_prefix = id_prefix
  )

  recommendation_table <- tibble::as_tibble(recommendation$recommendation_table %||% tibble::tibble())
  headline_table <- tibble::tibble(
    recommendation_available = isTRUE(recommendation$recommendation_available),
    n_designs = nrow(recommendation_table),
    recommended_design_id = as.character(recommendation$recommended_design_id %||% NA_character_),
    recommended_tier = as.character(recommendation$recommended_tier %||% NA_character_),
    n_top_tier_candidates = if ("structural_priority" %in% names(recommendation_table)) {
      min_priority <- min(recommendation_table$structural_priority, na.rm = TRUE)
      sum(recommendation_table$structural_priority == min_priority, na.rm = TRUE)
    } else {
      0L
    }
  )

  list(
    overview_contract = "arbitrary_facet_planning_design_review_recommendation_overview",
    overview_stage = as.character(recommendation$recommendation_stage %||% "deterministic_review"),
    planner_contract = as.character(
      recommendation$planner_contract %||% "named_facet_structural_review"
    ),
    recommendation_available = isTRUE(recommendation$recommendation_available),
    reason = as.character(
      recommendation$reason %||%
        "No structural-design structural recommendation overview is currently available."
    ),
    headline_table = headline_table,
    recommendation_table = recommendation_table,
    design_review_recommendation = recommendation,
    note = paste(
      "Compact structural recommendation overview for the structural-design review",
      "review layer, summarizing the conservative exact-tier ranking outcome."
    )
  )
}

structural_design_review_table_index <- function(active,
                                             overview,
                                             load_balance,
                                             coverage,
                                             guardrails,
                                             readiness,
                                             recommendation) {
  rows_or_zero <- function(tbl) {
    if (!is.data.frame(tbl)) return(0L)
    as.integer(nrow(tbl))
  }

  tibble::tibble(
    Table = c(
      "overview",
      "profile_summary",
      "load_balance_summary",
      "coverage_summary",
      "guardrail_summary",
      "readiness_summary",
      "recommendation_table"
    ),
    Rows = c(
      rows_or_zero(overview$headline_table),
      rows_or_zero(overview$metric_summary_table),
      rows_or_zero(load_balance$diagnostic_summary_table),
      rows_or_zero(coverage$diagnostic_summary_table),
      rows_or_zero(guardrails$guardrail_summary_table),
      rows_or_zero(readiness$readiness_summary_table),
      rows_or_zero(recommendation$recommendation_table)
    ),
    Role = c(
      "overview",
      "profile",
      "load_balance",
      "coverage",
      "guardrails",
      "readiness",
      "recommendation"
    ),
    Description = c(
      "Headline active-branch configuration and selected branch surface.",
      "Deterministic observation/load profile for the current design grid.",
      "Balanced-expectation and integer split summaries.",
      "Coverage/connectivity summaries implied by the current design grid.",
      "Exact overlap/balance regime summaries.",
      "Exact structural readiness indicators and tier summaries.",
      "Conservative deterministic recommendation ranking."
    ),
    stringsAsFactors = FALSE
  )
}

structural_design_appendix_selection_tables <- function(bundle,
                                                    label = "structural_design_review",
                                                    presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  presets <- unique(as.character(presets))
  original_bundles <- stats::setNames(list(bundle), label)

  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  selection_catalog <- bind_tables(lapply(presets, function(preset) {
    selected <- export_select_summary_table_bundles_for_appendix(
      original_bundles,
      preset = preset
    )
    export_summary_table_selection_catalog(
      original_bundles = original_bundles,
      selected_bundles = selected,
      preset = preset
    )
  }))

  list(
    selection_catalog = selection_catalog,
    selection_summary = bind_tables(lapply(presets, function(preset) {
      export_summary_table_selection_summary(selection_catalog, preset = preset)
    })),
    selection_role_summary = bind_tables(lapply(presets, function(preset) {
      export_summary_table_selection_role_summary(selection_catalog, preset = preset)
    })),
    selection_section_summary = bind_tables(lapply(presets, function(preset) {
      export_summary_table_selection_section_summary(selection_catalog, preset = preset)
    }))
  )
}

structural_design_selection_table_summary <- function(selection_catalog,
                                                  preset = NULL) {
  tbl <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    return(data.frame())
  }

  if (!is.null(preset)) {
    tbl <- tbl[as.character(tbl$Preset %||% "") %in% unique(as.character(preset)), , drop = FALSE]
  }
  tbl <- tbl[tbl$Selected %in% TRUE, , drop = FALSE]
  if (nrow(tbl) == 0L || !"Table" %in% names(tbl)) {
    return(data.frame())
  }

  compact_unique <- function(x, max_n = 4L) {
    summary_table_bundle_compact_labels(unique(as.character(x %||% character(0))), max_n = max_n)
  }
  first_non_missing <- function(x, default = NA_character_) {
    x <- as.character(x %||% character(0))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
      return(default)
    }
    x[[1]]
  }
  first_numeric <- function(x, default = NA_real_) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(default)
    }
    x[[1]]
  }

  split_tbl <- split(tbl, as.character(tbl$Table))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(table_nm) {
      part <- split_tbl[[table_nm]]
      data.frame(
        Table = as.character(table_nm),
        PresetsSelected = length(unique(as.character(part$Preset %||% ""))),
        Presets = compact_unique(part$Preset, max_n = 7L),
        Rows = as.integer(first_numeric(part$Rows, default = 0)),
        Role = first_non_missing(part$Role),
        AppendixSection = first_non_missing(part$AppendixSection),
        PreferredAppendixOrder = as.integer(first_numeric(part$PreferredAppendixOrder, default = 9999)),
        PlotReady = any(part$PlotReady %in% TRUE, na.rm = TRUE),
        ExportReady = any(part$ExportReady %in% TRUE, na.rm = TRUE),
        ApaTableReady = any(part$ApaTableReady %in% TRUE, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$PreferredAppendixOrder, out$Table, na.last = TRUE), , drop = FALSE]
}

structural_design_selection_table_preset_summary <- function(selection_catalog,
                                                         presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_table_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

structural_design_selection_handoff_table_summary <- function(selection_catalog,
                                                          presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_table_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

structural_design_selection_handoff_summary <- function(selection_catalog,
                                                    presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  tbl <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    return(data.frame())
  }

  presets <- unique(as.character(presets))
  tbl <- tbl[tbl$Selected %in% TRUE & as.character(tbl$Preset %||% "") %in% presets, , drop = FALSE]
  if (nrow(tbl) == 0L || !"AppendixSection" %in% names(tbl)) {
    return(data.frame())
  }

  compact_unique <- function(x, max_n = 4L) {
    summary_table_bundle_compact_labels(unique(as.character(x %||% character(0))), max_n = max_n)
  }

  split_tbl <- split(tbl, paste(as.character(tbl$Preset), as.character(tbl$AppendixSection), sep = "\r"))
  out <- do.call(
    rbind,
    lapply(split_tbl, function(part) {
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = as.character(part$Preset[[1]] %||% ""),
        AppendixSection = as.character(part$AppendixSection[[1]] %||% ""),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        RolesCovered = compact_unique(part$Role, max_n = 4L),
        KeyTables = compact_unique(part$Table, max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$AppendixSection), , drop = FALSE]
}

structural_design_selection_handoff_bundle_summary <- function(selection_catalog,
                                                           presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_bundle_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

structural_design_selection_handoff_preset_summary <- function(selection_catalog,
                                                           presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_preset_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

structural_design_selection_handoff_role_summary <- function(selection_catalog,
                                                         presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_role_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

structural_design_selection_handoff_role_section_summary <- function(selection_catalog,
                                                                 presets = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  bind_tables <- function(parts) {
    keep <- vapply(parts, function(df) is.data.frame(df) && nrow(df) > 0L, logical(1))
    parts <- parts[keep]
    if (length(parts) == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }

  presets <- unique(as.character(presets))
  bind_tables(lapply(presets, function(preset) {
    export_summary_table_selection_handoff_role_section_summary(
      selection_catalog = selection_catalog,
      preset = preset
    )
  }))
}

structural_design_review_plot_index_from_bundle <- function(plot_index) {
  idx <- as.data.frame(plot_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(idx) == 0L) {
    return(idx)
  }
  if (!all(c("Table", "DefaultPlotTypes") %in% names(idx))) {
    return(idx)
  }

  direct_routes <- c(
    structural_design_profile = "profile_metrics",
    structural_design_load_balance = "load_balance",
    structural_design_coverage = "coverage",
    structural_design_readiness = "readiness_tiers",
    structural_design_appendix_roles = "appendix_roles",
    structural_design_appendix_sections = "appendix_sections",
    structural_design_appendix_presets = "appendix_presets",
    structural_design_selection_table_presets = "selection_tables",
    structural_design_selection_handoff_presets = "selection_handoff_presets",
    structural_design_selection_handoff = "selection_handoff",
    structural_design_selection_handoff_bundles = "selection_handoff_bundles",
    structural_design_selection_handoff_roles = "selection_handoff_roles",
    structural_design_selection_handoff_role_sections = "selection_handoff_role_sections",
    structural_design_selection_tables = "selection_tables",
    structural_design_selection_summary = "selection_bundles",
    structural_design_selection_roles = "selection_roles",
    structural_design_selection_sections = "selection_sections"
  )

  idx$DefaultPlotTypes <- vapply(seq_len(nrow(idx)), function(i) {
    tbl <- as.character(idx$Table[i] %||% "")
    base_types <- strsplit(as.character(idx$DefaultPlotTypes[i] %||% ""), ",", fixed = TRUE)[[1]]
    base_types <- trimws(base_types)
    base_types <- base_types[nzchar(base_types)]
    extra <- unname(direct_routes[tbl])
    extra <- extra[!is.na(extra) & nzchar(extra)]
    paste(unique(c(extra, base_types)), collapse = ", ")
  }, character(1))

  idx
}

#' Summarize a deterministic named-facet design review
#'
#' @param object The structural design-review object stored in the
#'   normalized `planning_schema$structural_design_review` field.
#' @param digits Number of digits used in numeric summaries.
#' @param top_n Maximum number of recommendation rows to print in the preview.
#' @param detail Console detail: `"brief"` (default) prints the design overview,
#'   readiness, and recommendation; `"full"` prints all supporting tables.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is intentionally conservative. It aggregates deterministic
#' design quantities: observation bookkeeping, load/balance, coverage,
#' structural checks, and a conservative recommendation ranking.
#' It also exposes the same manuscript-facing table/appendix metadata used by
#' [build_summary_table_bundle()] so the design review can be inspected directly
#' without first routing through planning summaries. In addition to bundle-level
#' appendix presets and section counts, it includes export-like appendix
#' selection summaries by preset, reporting role, manuscript section,
#' bundle-aware handoff summaries, preset-specific table surface, and a
#' table-level handoff crosswalk, plus direct `role_summary` / `table_profile`
#' surfaces for table-shape review.
#' It does not report psychometric recovery or Monte Carlo performance.
#'
#' @return An object of class `summary.mfrm_structural_design_review`.
#' @seealso [summary.mfrm_design_evaluation()]
#' @noRd
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   design = list(person = 16, rater = 3, criterion = 2, assignment = 2),
#'   assignment = "rotating"
#' )
#' active <- spec$planning_schema$structural_design_review
#' summary(active)
#' }
#' @export
summary.mfrm_structural_design_review <- function(object,
                                                     digits = 3,
                                                     top_n = 8,
                                                     detail = c("brief", "full"),
                                                     ...) {
  active <- simulation_structural_design_review(object)
  overview <- simulation_structural_design_review_overview(active)
  load_balance <- simulation_structural_design_review_load_balance_overview(active)
  coverage <- simulation_structural_design_review_coverage_overview(active)
  guardrails <- simulation_structural_design_review_guardrail_overview(active)
  readiness <- simulation_structural_design_review_readiness_overview(active)
  recommendation <- simulation_structural_design_review_recommendation_overview(active)

  digits <- max(0L, as.integer(digits[1]))
  top_n <- max(1L, as.integer(top_n[1]))
  detail <- match.arg(detail)

  round_df <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  headline <- tibble::tibble(
    review_available = isTRUE(active$review_available),
    n_designs = as.integer(active$n_designs %||% 0L),
    recommended_design_id = as.character(active$recommended_design_id %||% NA_character_),
    view = as.character(active$view %||% NA_character_),
    mode = as.character(active$mode %||% NA_character_),
    surface = as.character(active$surface %||% NA_character_),
    table_component = as.character(active$table_component %||% NA_character_)
  )

  recommendation_table <- tibble::as_tibble(recommendation$recommendation_table %||% tibble::tibble())
  if (nrow(recommendation_table) > 0L) {
    recommendation_table <- utils::head(recommendation_table, n = top_n)
  }

  notes <- c(
    paste(
      "Design planning varies the documented person-count, non-person",
      "facet-count, and assignment axes."
    ),
    paste(
      "Readiness summarizes structural overlap and balance; it does not",
      "estimate parameter recovery or Monte Carlo performance."
    )
  )

  out <- list(
    overview = round_df(headline),
    table_index = structural_design_review_table_index(
      active = active,
      overview = overview,
      load_balance = load_balance,
      coverage = coverage,
      guardrails = guardrails,
      readiness = readiness,
      recommendation = recommendation
    ),
    profile_summary = round_df(tibble::as_tibble(overview$metric_summary_table %||% tibble::tibble())),
    load_balance_summary = round_df(tibble::as_tibble(load_balance$diagnostic_summary_table %||% tibble::tibble())),
    coverage_summary = round_df(tibble::as_tibble(coverage$diagnostic_summary_table %||% tibble::tibble())),
    guardrail_summary = round_df(tibble::as_tibble(guardrails$guardrail_summary_table %||% tibble::tibble())),
    readiness_summary = round_df(tibble::as_tibble(readiness$readiness_summary_table %||% tibble::tibble())),
    recommendation_table = round_df(recommendation_table),
    notes = notes,
    digits = digits,
    detail = detail
  )
  temp_core <- out
  class(temp_core) <- "summary.mfrm_structural_design_review"
  bundle_core <- build_summary_table_bundle(temp_core, include_empty = TRUE)
  selection_tables <- structural_design_appendix_selection_tables(bundle_core)
  out$selection_table_summary <- structural_design_selection_table_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_table_preset_summary <- structural_design_selection_table_preset_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_table_summary <- structural_design_selection_handoff_table_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_preset_summary <- structural_design_selection_handoff_preset_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_summary <- structural_design_selection_handoff_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_bundle_summary <- structural_design_selection_handoff_bundle_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_role_summary <- structural_design_selection_handoff_role_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_handoff_role_section_summary <- structural_design_selection_handoff_role_section_summary(
    selection_catalog = selection_tables$selection_catalog
  )
  out$selection_summary <- selection_tables$selection_summary
  out$selection_role_summary <- selection_tables$selection_role_summary
  out$selection_section_summary <- selection_tables$selection_section_summary
  out$selection_catalog <- selection_tables$selection_catalog

  temp_final <- out
  class(temp_final) <- "summary.mfrm_structural_design_review"
  bundle_final <- build_summary_table_bundle(temp_final, include_empty = TRUE)
  out$table_index <- as.data.frame(bundle_final$table_index %||% data.frame(), stringsAsFactors = FALSE)
  out$plot_index <- structural_design_review_plot_index_from_bundle(
    as.data.frame(bundle_final$plot_index %||% data.frame(), stringsAsFactors = FALSE)
  )
  out$table_catalog <- summary_table_bundle_catalog(bundle_final)
  out$table_profile <- summary_table_bundle_profile(bundle_final)
  if (nrow(out$table_profile) > 0L) {
    ord <- order(out$table_profile$Rows, out$table_profile$Cols, decreasing = TRUE, na.last = TRUE)
    out$table_profile <- out$table_profile[ord, , drop = FALSE]
    out$table_profile <- utils::head(out$table_profile, n = top_n)
  }
  out$role_summary <- data.frame()
  if (nrow(out$table_index) > 0L && "Role" %in% names(out$table_index)) {
    roles <- split(out$table_index, out$table_index$Role %||% "")
    out$role_summary <- do.call(
      rbind,
      lapply(names(roles), function(role_nm) {
        part <- roles[[role_nm]]
        data.frame(
          Role = as.character(role_nm),
          Tables = nrow(part),
          TotalRows = sum(suppressWarnings(as.numeric(part$Rows)), na.rm = TRUE),
          TotalCols = sum(suppressWarnings(as.numeric(part$Cols)), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
    out$role_summary <- out$role_summary[
      order(out$role_summary$Tables, out$role_summary$Role, decreasing = TRUE),
      ,
      drop = FALSE
    ]
    rownames(out$role_summary) <- NULL
  }
  out$appendix_presets <- summary_table_bundle_appendix_presets(out$table_catalog)
  out$appendix_role_summary <- summary_table_bundle_appendix_role_summary(out$table_catalog)
  out$appendix_section_summary <- summary_table_bundle_appendix_section_summary(out$table_catalog)
  out$reporting_map <- summary_table_bundle_reporting_map(bundle_final, out$table_catalog)
  out$overview$RecommendedAppendixTables <- sum(out$table_catalog$RecommendedAppendix %in% TRUE, na.rm = TRUE)
  out$overview$CompactAppendixTables <- sum(out$table_catalog$CompactAppendix %in% TRUE, na.rm = TRUE)
  out$overview$NumericTables <- sum(out$table_profile$NumericColumns > 0, na.rm = TRUE)
  out$overview$AnyNumericTable <- nrow(out$table_profile) > 0L &&
    any(out$table_profile$NumericColumns > 0, na.rm = TRUE)
  class(out) <- "summary.mfrm_structural_design_review"
  out
}

#' @export
print.summary.mfrm_structural_design_review <- function(x, ...) {
  digits <- max(0L, as.integer(x$digits %||% 3L))

  cat("mfrmr Structural Design Review\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (identical(as.character(x$detail %||% "brief"), "brief")) {
    if (!is.null(x$readiness_summary) && nrow(x$readiness_summary) > 0L) {
      cat("\nReadiness summary\n")
      print(round_numeric_df(as.data.frame(x$readiness_summary), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$recommendation_table) && nrow(x$recommendation_table) > 0L) {
      cat("\nRecommendation table\n")
      print(round_numeric_df(as.data.frame(x$recommendation_table), digits = digits), row.names = FALSE)
    }
    if (length(x$notes %||% character(0)) > 0L) {
      cat("\nNotes\n")
      for (line in x$notes) cat(" - ", line, "\n", sep = "")
    }
    cat("\nAdditional design tables remain in the structured summary; use `detail = \"full\"` to print them.\n")
    return(invisible(x))
  }
  if (!is.null(x$table_index) && nrow(x$table_index) > 0L) {
    cat("\nTable index\n")
    print(as.data.frame(x$table_index), row.names = FALSE)
  }
  if (!is.null(x$plot_index) && nrow(x$plot_index) > 0L) {
    cat("\nPlot index\n")
    print(as.data.frame(x$plot_index), row.names = FALSE)
  }
  if (!is.null(x$table_catalog) && nrow(x$table_catalog) > 0L) {
    cat("\nTable catalog\n")
    print(round_numeric_df(as.data.frame(x$table_catalog), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$role_summary) && nrow(x$role_summary) > 0L) {
    cat("\nRole summary\n")
    print(round_numeric_df(as.data.frame(x$role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$table_profile) && nrow(x$table_profile) > 0L) {
    cat("\nTable profile\n")
    print(round_numeric_df(as.data.frame(x$table_profile), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$profile_summary) && nrow(x$profile_summary) > 0L) {
    cat("\nProfile summary\n")
    print(round_numeric_df(as.data.frame(x$profile_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$load_balance_summary) && nrow(x$load_balance_summary) > 0L) {
    cat("\nLoad/balance summary\n")
    print(round_numeric_df(as.data.frame(x$load_balance_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$coverage_summary) && nrow(x$coverage_summary) > 0L) {
    cat("\nCoverage summary\n")
    print(round_numeric_df(as.data.frame(x$coverage_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$guardrail_summary) && nrow(x$guardrail_summary) > 0L) {
    cat("\nGuardrail summary\n")
    print(round_numeric_df(as.data.frame(x$guardrail_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$readiness_summary) && nrow(x$readiness_summary) > 0L) {
    cat("\nReadiness summary\n")
    print(round_numeric_df(as.data.frame(x$readiness_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$recommendation_table) && nrow(x$recommendation_table) > 0L) {
    cat("\nRecommendation table\n")
    print(round_numeric_df(as.data.frame(x$recommendation_table), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_presets) && nrow(x$appendix_presets) > 0L) {
    cat("\nAppendix presets\n")
    print(round_numeric_df(as.data.frame(x$appendix_presets), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_role_summary) && nrow(x$appendix_role_summary) > 0L) {
    cat("\nAppendix role summary\n")
    print(round_numeric_df(as.data.frame(x$appendix_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$appendix_section_summary) && nrow(x$appendix_section_summary) > 0L) {
    cat("\nAppendix section summary\n")
    print(round_numeric_df(as.data.frame(x$appendix_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_summary) && nrow(x$selection_summary) > 0L) {
    cat("\nSelection summary\n")
    print(round_numeric_df(as.data.frame(x$selection_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_table_summary) && nrow(x$selection_table_summary) > 0L) {
    cat("\nSelection table summary\n")
    print(round_numeric_df(as.data.frame(x$selection_table_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_table_preset_summary) && nrow(x$selection_table_preset_summary) > 0L) {
    cat("\nSelection table preset summary\n")
    print(round_numeric_df(as.data.frame(x$selection_table_preset_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_table_summary) && nrow(x$selection_handoff_table_summary) > 0L) {
    cat("\nSelection handoff table summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_table_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_preset_summary) && nrow(x$selection_handoff_preset_summary) > 0L) {
    cat("\nSelection handoff preset summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_preset_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_summary) && nrow(x$selection_handoff_summary) > 0L) {
    cat("\nSelection handoff summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_bundle_summary) && nrow(x$selection_handoff_bundle_summary) > 0L) {
    cat("\nSelection handoff bundle summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_bundle_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_role_summary) && nrow(x$selection_handoff_role_summary) > 0L) {
    cat("\nSelection handoff role summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_handoff_role_section_summary) && nrow(x$selection_handoff_role_section_summary) > 0L) {
    cat("\nSelection handoff role-section summary\n")
    print(round_numeric_df(as.data.frame(x$selection_handoff_role_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_role_summary) && nrow(x$selection_role_summary) > 0L) {
    cat("\nSelection role summary\n")
    print(round_numeric_df(as.data.frame(x$selection_role_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$selection_section_summary) && nrow(x$selection_section_summary) > 0L) {
    cat("\nSelection section summary\n")
    print(round_numeric_df(as.data.frame(x$selection_section_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$reporting_map) && nrow(x$reporting_map) > 0L) {
    cat("\nReporting map\n")
    print(as.data.frame(x$reporting_map), row.names = FALSE)
  }
  if (length(x$notes %||% character(0)) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

simulation_compact_structural_design_review_summary <- function(x,
                                                            digits = 3,
                                                            top_n = 6L) {
  planning_schema <- simulation_object_planning_schema(x)
  design_review <- planning_schema$structural_design_review %||% NULL
  if (!is.list(design_review)) {
    return(NULL)
  }

  out <- tryCatch(
    summary.mfrm_structural_design_review(
      design_review,
      digits = digits,
      top_n = top_n
    ),
    error = function(e) NULL
  )
  if (!inherits(out, "summary.mfrm_structural_design_review")) {
    return(NULL)
  }

  out
}

print_compact_structural_design_review_summary <- function(x,
                                                       digits = 3,
                                                       heading = "Structural design review") {
  if (!inherits(x, "summary.mfrm_structural_design_review")) {
    return(invisible(NULL))
  }

  cat("\n", heading, "\n", sep = "")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$readiness_summary) && nrow(x$readiness_summary) > 0L) {
    cat("\nReadiness summary\n")
    print(round_numeric_df(as.data.frame(x$readiness_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$recommendation_table) && nrow(x$recommendation_table) > 0L) {
    cat("\nRecommendation table\n")
    print(round_numeric_df(as.data.frame(x$recommendation_table), digits = digits), row.names = FALSE)
  }

  invisible(x)
}

#' Plot a deterministic named-facet design review
#'
#' @param x The structural design-review object stored in the
#'   normalized `planning_schema$structural_design_review` field.
#' @param y Unused placeholder for generic compatibility.
#' @param type Plot type: `"profile_metrics"` for recommended deterministic
#'   profile values by metric, `"load_balance"` for recommended load/balance
#'   values by metric, `"coverage"` for recommended coverage/connectivity values
#'   by metric, `"readiness_tiers"` for counts of structural tiers across the
#'   current design grid, `"table_rows"` / `"role_tables"` /
#'   `"appendix_roles"` for summary-table bundle QC,
#'   `"appendix_sections"` / `"appendix_presets"` for manuscript-facing
#'   appendix selection counts, `"selection_handoff_presets"` for preset-level
#'   appendix handoff counts, `"selection_tables"` for appendix-selected
#'   design-review tables ranked by row count within a preset,
#'   `"selection_handoff"` for section-aware plot-ready appendix handoff counts,
#'   `"selection_handoff_bundles"` for section-and-bundle plot-ready appendix
#'   handoff counts,
#'   `"selection_handoff_roles"` for role-aware plot-ready appendix handoff
#'   counts, `"selection_handoff_role_sections"` for role-by-section plot-ready
#'   appendix handoff counts,
#'   or `"selection_bundles"` / `"selection_roles"` / `"selection_sections"`
#'   for preset-filtered appendix selection summaries.
#' @param appendix_preset Appendix preset used for `selection_*` plot types.
#' @param selection_value For `selection_*` plot types, whether to plot exact
#'   counts (`"count"`) or the matching exact fraction (`"fraction"`) when that
#'   surface exposes one. `selection_tables` remains count-only because it
#'   represents table row counts rather than a normalized selection surface.
#' @param draw If `TRUE`, draw with base graphics; otherwise return plotting data.
#' @param main Optional title override.
#' @param palette Optional named color overrides.
#' @param label_angle Axis-label rotation angle.
#' @param ... Reserved for generic compatibility.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @noRd
#' @examples
#' \donttest{
#' spec <- build_mfrm_sim_spec(
#'   design = list(person = 16, rater = 3, criterion = 2, assignment = 2),
#'   assignment = "rotating"
#' )
#' active <- spec$planning_schema$structural_design_review
#' plot(active, type = "readiness_tiers", draw = FALSE)
#' }
#' @export
plot.mfrm_structural_design_review <- function(x,
                                                  y = NULL,
                                                  type = c("profile_metrics", "load_balance", "coverage", "readiness_tiers", "table_rows", "role_tables", "appendix_roles", "appendix_sections", "appendix_presets", "selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections"),
                                                  appendix_preset = c("recommended", "compact", "all", "methods", "results", "diagnostics", "reporting"),
                                                  selection_value = c("count", "fraction"),
                                                  draw = TRUE,
                                                  main = NULL,
                                                  palette = NULL,
                                                  label_angle = 45,
                                                  ...) {
  type <- match.arg(type)
  appendix_preset <- match.arg(
    tolower(as.character(appendix_preset[1])),
    choices = c("recommended", "compact", "all", "methods", "results", "diagnostics", "reporting")
  )
  selection_value <- match.arg(selection_value)
  active <- simulation_structural_design_review(x)

  if (type %in% c("table_rows", "role_tables", "appendix_roles", "appendix_sections", "appendix_presets")) {
    bundle <- build_summary_table_bundle(active)
    return(plot.mfrm_summary_table_bundle(
      bundle,
      type = type,
      selection_value = selection_value,
      main = main,
      palette = palette,
      label_angle = label_angle,
      draw = draw,
      ...
    ))
  }

  if (type %in% c("selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections")) {
    sx <- summary.mfrm_structural_design_review(active)
    if (type == "selection_handoff_presets") {
      tbl <- as.data.frame(sx$selection_handoff_preset_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Preset", "PlotReadyTables") %in% names(tbl))) {
        stop("No structural-design appendix handoff-preset summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Preset)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by preset for `", appendix_preset, "`")
      default_palette <- c(selection_handoff_presets = "#f4a261", grid = "#ececec")
      plot_name <- "selection_handoff_presets"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_tables") {
      tbl <- structural_design_selection_table_summary(
        selection_catalog = sx$selection_catalog,
        preset = appendix_preset
      )
      if (nrow(tbl) == 0L || !all(c("Table", "Rows") %in% names(tbl))) {
        stop("No structural-design appendix table selection is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Table)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Selected appendix tables for preset `", appendix_preset, "`")
      default_palette <- c(selection_tables = "#e76f51", grid = "#ececec")
      plot_name <- "selection_tables"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff") {
      tbl <- as.data.frame(sx$selection_handoff_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "PlotReadyTables") %in% names(tbl))) {
        stop("No structural-design appendix handoff summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$AppendixSection)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by section for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff = "#ff9f1c", grid = "#ececec")
      plot_name <- "selection_handoff"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff_bundles") {
      tbl <- as.data.frame(sx$selection_handoff_bundle_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Bundle", "PlotReadyTables") %in% names(tbl))) {
        stop("No structural-design appendix handoff-bundle summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Bundle))
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by section and bundle for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff_bundles = "#5c677d", grid = "#ececec")
      plot_name <- "selection_handoff_bundles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff_roles") {
      tbl <- as.data.frame(sx$selection_handoff_role_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Role", "PlotReadyTables") %in% names(tbl))) {
        stop("No structural-design appendix handoff-role summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Role)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by role for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff_roles = "#9c6644", grid = "#ececec")
      plot_name <- "selection_handoff_roles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_handoff_role_sections") {
      tbl <- as.data.frame(sx$selection_handoff_role_section_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Role", "PlotReadyTables") %in% names(tbl))) {
        stop("No structural-design appendix handoff role-section summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Role))
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Plot-ready appendix handoff by section and role for preset `", appendix_preset, "`")
      default_palette <- c(selection_handoff_role_sections = "#7f5539", grid = "#ececec")
      plot_name <- "selection_handoff_role_sections"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_bundles") {
      tbl <- as.data.frame(sx$selection_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Bundle", "TablesSelected") %in% names(tbl))) {
        stop("No structural-design appendix selection summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Bundle)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Appendix tables by bundle for preset `", appendix_preset, "`")
      default_palette <- c(selection_bundles = "#54a24b", grid = "#ececec")
      plot_name <- "selection_bundles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else if (type == "selection_roles") {
      tbl <- as.data.frame(sx$selection_role_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("Role", "Tables") %in% names(tbl))) {
        stop("No structural-design appendix role summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$Role)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Selected appendix roles for preset `", appendix_preset, "`")
      default_palette <- c(selection_roles = "#b279a2", grid = "#ececec")
      plot_name <- "selection_roles"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    } else {
      tbl <- as.data.frame(sx$selection_section_summary %||% data.frame(), stringsAsFactors = FALSE)
      tbl <- tbl[as.character(tbl$Preset %||% "") %in% appendix_preset, , drop = FALSE]
      if (nrow(tbl) == 0L || !all(c("AppendixSection", "Tables") %in% names(tbl))) {
        stop("No structural-design appendix section summary is available for preset `", appendix_preset, "`.", call. = FALSE)
      }
      labels <- as.character(tbl$AppendixSection)
      measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
      values <- measure$values
      subtitle <- paste0("Selected appendix sections for preset `", appendix_preset, "`")
      default_palette <- c(selection_sections = "#2a9d8f", grid = "#ececec")
      plot_name <- "selection_sections"
      legend_label <- measure$legend_label
      ylab <- measure$ylab
    }

    keep <- is.finite(values) & nzchar(labels)
    if (!any(keep)) {
      stop("Selected structural-design appendix plot has no finite values to display.", call. = FALSE)
    }
    labels <- labels[keep]
    values <- values[keep]
    tbl <- tbl[keep, , drop = FALSE]
    pal <- resolve_palette(palette = palette, defaults = default_palette)
    plot_title <- if (is.null(main)) paste("Structural design review:", gsub("_", " ", plot_name)) else as.character(main[1])

    if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = pal[plot_name],
        main = plot_title,
        ylab = ylab,
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }

    return(invisible(new_mfrm_plot_data(
      "structural_design_review",
      list(
        plot = plot_name,
        selection_value = measure$selection_value,
        appendix_preset = appendix_preset,
        table = tbl,
        title = plot_title,
        subtitle = subtitle,
        legend = new_plot_legend(legend_label, "structural_design", "bar", pal[plot_name]),
        reference_lines = new_reference_lines("h", 0, "Zero-table reference", "dashed", "reference")
      )
    )))
  }

  if (type == "profile_metrics") {
    tbl <- simulation_structural_design_review_profile(active)$profile_summary_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !all(c("metric", "recommended_value") %in% names(tbl))) {
      stop("No structural-design profile metrics are available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$metric)
    values <- suppressWarnings(as.numeric(tbl$recommended_value))
    subtitle <- "Recommended deterministic profile values"
    default_palette <- c(profile_metrics = "#3a7ca5", grid = "#ececec")
    plot_name <- "profile_metrics"
    legend_label <- "Recommended value"
  } else if (type == "load_balance") {
    tbl <- simulation_structural_design_review_load_balance(active)$diagnostic_summary_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !all(c("metric", "recommended_value") %in% names(tbl))) {
      stop("No structural-design load/balance metrics are available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$metric)
    values <- suppressWarnings(as.numeric(tbl$recommended_value))
    subtitle <- "Recommended deterministic load/balance values"
    default_palette <- c(load_balance = "#2a9d8f", grid = "#ececec")
    plot_name <- "load_balance"
    legend_label <- "Recommended value"
  } else if (type == "coverage") {
    tbl <- simulation_structural_design_review_coverage(active)$diagnostic_summary_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !all(c("metric", "recommended_value") %in% names(tbl))) {
      stop("No structural-design coverage metrics are available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$metric)
    values <- suppressWarnings(as.numeric(tbl$recommended_value))
    subtitle <- "Recommended deterministic coverage/connectivity values"
    default_palette <- c(coverage = "#f4a261", grid = "#ececec")
    plot_name <- "coverage"
    legend_label <- "Recommended value"
  } else {
    tbl <- simulation_structural_design_review_readiness(active)$readiness_table
    tbl <- tibble::as_tibble(tbl %||% tibble::tibble())
    if (nrow(tbl) == 0L || !"structural_tier" %in% names(tbl)) {
      stop("No structural-design readiness tiers are available for plotting.", call. = FALSE)
    }
    counts <- sort(table(as.character(tbl$structural_tier)), decreasing = TRUE)
    labels <- names(counts)
    values <- as.numeric(counts)
    subtitle <- "Structural tier counts across candidate designs"
    default_palette <- c(readiness_tiers = "#b279a2", grid = "#ececec")
    plot_name <- "readiness_tiers"
    legend_label <- "Designs"
  }

  keep <- is.finite(values) & nzchar(labels)
  if (!any(keep)) {
    stop("Selected structural-design plot has no finite values to display.", call. = FALSE)
  }
  labels <- labels[keep]
  values <- values[keep]
  pal <- resolve_palette(palette = palette, defaults = default_palette)
  color_name <- names(default_palette)[1]
  plot_title <- if (is.null(main)) "Structural design review" else as.character(main[1])

  if (isTRUE(draw)) {
    barplot_rot45(
      height = values,
      labels = labels,
      col = pal[color_name],
      main = plot_title,
      ylab = legend_label,
      label_angle = label_angle,
      mar_bottom = 9
    )
    graphics::abline(h = 0, col = pal["grid"], lty = 2)
  }

  invisible(new_mfrm_plot_data(
    "structural_design_review",
    list(
      plot = plot_name,
      label = labels,
      value = values,
      title = plot_title,
      subtitle = subtitle,
      recommended_design_id = as.character(active$recommended_design_id %||% NA_character_),
      legend = new_plot_legend(legend_label, "structural_design", "bar", pal[color_name]),
      reference_lines = new_reference_lines("h", 0, "Zero reference", "dashed", "reference")
    )
  ))
}

simulation_compact_structural_design_report_cache <- function(x) {
  if (!is.list(x) || is.data.frame(x) || inherits(x, "mfrm_plot_data")) {
    return(x)
  }

  nested_report_fields <- c(
    "report_bundle",
    "report_summary",
    "report_overview",
    "report_catalog",
    "report_surface_registry",
    "report_panel",
    "source_object",
    "registry",
    "selected_surface"
  )
  x[intersect(names(x), nested_report_fields)] <- NULL

  for (nm in names(x)) {
    if (is.list(x[[nm]]) && !is.data.frame(x[[nm]]) && !inherits(x[[nm]], "mfrm_plot_data")) {
      x[[nm]] <- simulation_compact_structural_design_report_cache(x[[nm]])
    }
  }

  x
}

simulation_compact_structural_design_schema_report_cache <- function(branch_schema) {
  report_fields <- c(
    "report_summary",
    "report_overview_table",
    "report_catalog",
    "report_digest",
    "report_surface_registry",
    "report_panel",
    "report_operation",
    "report_snapshot",
    "report_brief",
    "report_consumer"
  )
  for (nm in intersect(report_fields, names(branch_schema))) {
    branch_schema[[nm]] <- simulation_compact_structural_design_report_cache(branch_schema[[nm]])
  }
  branch_schema
}

simulation_compact_structural_design_schema_active_cache <- function(branch_schema) {
  if (is.list(branch_schema$design_review)) {
    if (is.list(branch_schema$design_review$summary)) {
      branch_schema$design_review$summary$review <- NULL
    }
    if (is.list(branch_schema$design_review$table)) {
      branch_schema$design_review$table$review <- NULL
    }
    if (is.list(branch_schema$design_review$plot)) {
      branch_schema$design_review$plot$review <- NULL
    }
  }

  drop_map <- list(
    compact_review_summary = "review",
    compact_review_table = "review",
    compact_review_plot = "review",
    design_review_profile = "design_review",
    design_review_load_balance = "design_review_profile",
    design_review_overview = c("design_review", "design_review_profile"),
    design_review_load_balance_overview = "design_review_load_balance",
    design_review_coverage = "design_review_profile",
    design_review_coverage_overview = "design_review_coverage",
    design_review_guardrails = c("design_review_coverage", "design_review_load_balance"),
    design_review_guardrail_overview = "design_review_guardrails",
    design_review_readiness = "design_review_guardrails",
    design_review_readiness_overview = "design_review_readiness",
    design_review_recommendation = c("design_review_readiness", "design_review_profile"),
    design_review_recommendation_overview = "design_review_recommendation"
  )

  for (nm in intersect(names(drop_map), names(branch_schema))) {
    branch_schema[[nm]][intersect(drop_map[[nm]], names(branch_schema[[nm]]))] <- NULL
  }

  branch_schema
}

simulation_structural_design_schema <- function(sim_spec = NULL, facet_names = NULL) {
  structural_facet_table <- if (is.null(facet_names)) {
    simulation_structural_facet_table(sim_spec)
  } else {
    simulation_structural_facet_table(facet_names = facet_names)
  }
  structural_design_template <- if (is.null(facet_names)) {
    simulation_structural_design_template(sim_spec)
  } else {
    simulation_structural_design_template(facet_names = facet_names)
  }
  design_schema <- if (is.null(facet_names)) {
    simulation_structural_design_contract(sim_spec)
  } else {
    simulation_structural_design_contract(facet_names = facet_names)
  }
  branch_schema <- list(
    design_review_contract = "structural_design_review",
    design_review_stage = "supported",
    planner_contract = "named_facet_structural_review",
    planner_stage = "deterministic_review",
    input_contract = "design$facets(named counts)",
    facet_table = structural_facet_table,
    design_template = structural_design_template,
    assignment_axis = design_schema$assignment_axis,
    design_schema = design_schema,
    grid_semantics = design_schema$grid_semantics,
    note = paste(
      "Normalized structural design schema bundling the named",
      "facet-count table, the matching `design$facets(named counts)` template,",
      "and preview-ready design-review metadata."
    )
  )
  branch_schema$grid_contract <- simulation_structural_design_grid_contract(branch_schema)
  branch_schema$preview <- simulation_structural_design_preview(branch_schema)
  branch_schema$grid_bundle <- simulation_structural_design_grid_bundle(branch_schema)
  branch_schema$grid_context <- simulation_structural_design_grid_context(branch_schema)
  branch_schema$report_bundle <- simulation_structural_design_report_bundle(branch_schema)
  branch_schema$report_summary <- simulation_structural_design_report_summary(branch_schema)
  branch_schema$report_overview_table <- simulation_structural_design_report_overview_table(branch_schema)
  branch_schema$report_catalog <- simulation_structural_design_report_catalog(branch_schema)
  branch_schema$report_digest <- simulation_structural_design_report_digest(branch_schema)
  branch_schema$report_surface_registry <- simulation_structural_design_report_surface_registry(branch_schema)
  branch_schema$report_panel <- simulation_structural_design_report_panel(branch_schema)
  branch_schema$report_operation <- simulation_structural_design_report_operation(branch_schema)
  branch_schema$report_snapshot <- simulation_structural_design_report_snapshot(branch_schema)
  branch_schema$report_brief <- simulation_structural_design_report_brief(branch_schema)
  branch_schema$report_mode_registry <- simulation_structural_design_report_mode_registry(branch_schema)
  branch_schema$report_consumer <- simulation_structural_design_report_consume(branch_schema)
  branch_schema$review <- simulation_structural_design_compact_review(branch_schema)
  branch_schema$compact_review_summary <- simulation_structural_design_compact_review_summary(branch_schema)
  branch_schema$compact_review_table <- simulation_structural_design_compact_review_table(branch_schema)
  branch_schema$compact_review_plot <- simulation_structural_design_compact_review_plot(branch_schema)
  branch_schema$design_review <- simulation_structural_design_review(branch_schema)
  branch_schema$design_review_profile <- simulation_structural_design_review_profile(branch_schema)
  branch_schema$design_review_load_balance <- simulation_structural_design_review_load_balance(branch_schema)
  branch_schema$design_review_overview <- simulation_structural_design_review_overview(branch_schema)
  branch_schema$design_review_load_balance_overview <- simulation_structural_design_review_load_balance_overview(branch_schema)
  branch_schema$design_review_coverage <- simulation_structural_design_review_coverage(branch_schema)
  branch_schema$design_review_coverage_overview <- simulation_structural_design_review_coverage_overview(branch_schema)
  branch_schema$design_review_guardrails <- simulation_structural_design_review_guardrails(branch_schema)
  branch_schema$design_review_guardrail_overview <- simulation_structural_design_review_guardrail_overview(branch_schema)
  branch_schema$design_review_readiness <- simulation_structural_design_review_readiness(branch_schema)
  branch_schema$design_review_readiness_overview <- simulation_structural_design_review_readiness_overview(branch_schema)
  branch_schema$design_review_recommendation <- simulation_structural_design_review_recommendation(branch_schema)
  branch_schema$design_review_recommendation_overview <- simulation_structural_design_review_recommendation_overview(branch_schema)
  branch_schema <- simulation_compact_structural_design_schema_report_cache(branch_schema)
  branch_schema <- simulation_compact_structural_design_schema_active_cache(branch_schema)
  branch_schema
}

simulation_planning_scope <- function(sim_spec = NULL, facet_names = NULL) {
  if (is.null(facet_names)) {
    facet_names <- simulation_spec_output_facet_names(sim_spec)
    aliases <- simulation_design_variable_aliases(sim_spec)
    facet_manifest <- simulation_facet_manifest(sim_spec)
  } else {
    facet_names <- simulation_validate_output_facet_names(facet_names)
    aliases <- simulation_design_variable_aliases(list(facet_names = facet_names))
    facet_manifest <- simulation_facet_manifest(facet_names = facet_names)
  }
  list(
    design_review_contract = "role_based_named_facets",
    design_review_stage = "supported",
    planner_contract = "role_based_two_non_person_facets",
    planner_stage = "role_based",
    planning_mode = "role_based",
    supports_arbitrary_facet_planning = FALSE,
    supports_arbitrary_facet_estimation = TRUE,
    supported_non_person_roles = c("rater", "criterion"),
    supported_non_person_facet_count = 2L,
    role_labels = unname(facet_names),
    design_variables = c("n_person", "n_rater", "n_criterion", "raters_per_person"),
    design_variable_aliases = aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")],
    facet_manifest = facet_manifest,
    structural_review_contract = "named_facet_structural_review",
    structural_review_stage = "deterministic_review",
    structural_design_input_contract = "design$facets(named counts)",
    note = paste0(
      "Planning helpers vary one person count and two named non-person facet roles (",
      facet_names[1], " and ", facet_names[2], "). ",
      "Estimation may contain additional facets, but planning and forecasting are limited to this role-based design."
    )
  )
}

simulation_planning_scope_note <- function(scope) {
  if (is.list(scope) && is.character(scope$note) && length(scope$note) > 0L) {
    note <- as.character(scope$note[1])
    if (nzchar(note)) {
      return(note)
    }
  }
  character(0)
}

simulation_planning_constraints <- function(sim_spec = NULL) {
  all_vars <- c("n_person", "n_rater", "n_criterion", "raters_per_person")
  if (is.null(sim_spec) || !inherits(sim_spec, "mfrm_sim_spec")) {
    return(list(
      mutable_design_variables = all_vars,
      locked_design_variables = character(0),
      lock_reasons = stats::setNames(character(0), character(0)),
      feasibility_rule = "`raters_per_person <= n_rater`",
      note = "Current scalar-argument planning paths allow `n_person`, `n_rater`, `n_criterion`, and `raters_per_person` to vary subject to `raters_per_person <= n_rater`."
    ))
  }

  mutable <- all_vars
  reasons <- character(0)
  assignment <- as.character(sim_spec$assignment %||% "rotating")

  if (identical(assignment, "resampled")) {
    reasons["n_rater"] <- paste0(
      "`assignment = \"resampled\"` reuses empirical person-level ",
      simulation_spec_output_facet_names(sim_spec)[1],
      " profiles."
    )
    reasons["raters_per_person"] <- reasons[["n_rater"]]
  }
  if (identical(assignment, "skeleton")) {
    skeleton_reason <- "`assignment = \"skeleton\"` reuses the observed person-by-facet response skeleton."
    reasons["n_rater"] <- skeleton_reason
    reasons["n_criterion"] <- skeleton_reason
    reasons["raters_per_person"] <- skeleton_reason
  }

  threshold_mode <- simulation_spec_threshold_mode(sim_spec)
  if (identical(threshold_mode, "step_facet_specific")) {
    role <- simulation_step_facet_role(sim_spec, step_facet = sim_spec$step_facet)
    lock_var <- if (identical(role, "criterion")) "n_criterion" else if (identical(role, "rater")) "n_rater" else NULL
    if (!is.null(lock_var)) {
      reasons[lock_var] <- paste0(
        "`sim_spec` contains step-facet-specific thresholds for `",
        sim_spec$step_facet,
        "`."
      )
    }
  }

  slope_mode <- simulation_spec_slope_mode(sim_spec)
  if (identical(slope_mode, "slope_facet_specific")) {
    role <- simulation_step_facet_role(sim_spec, step_facet = sim_spec$slope_facet)
    lock_var <- if (identical(role, "criterion")) "n_criterion" else if (identical(role, "rater")) "n_rater" else NULL
    if (!is.null(lock_var)) {
      reasons[lock_var] <- paste0(
        "The bounded `GPCM` specification stores slope values for `",
        sim_spec$slope_facet,
        "` levels."
      )
    }
  }

  reason_names <- names(reasons)
  if (is.null(reason_names)) {
    reason_names <- character(0)
  }
  locked <- intersect(all_vars, reason_names)
  mutable <- setdiff(all_vars, locked)
  note <- if (length(locked) == 0L) {
    "All current design variables remain mutable subject to `raters_per_person <= n_rater`."
  } else {
    paste0(
      "Current planning path allows changing ",
      paste(mutable, collapse = ", "),
      "; locked variables are ",
      paste(locked, collapse = ", "),
      "."
    )
  }

  list(
    mutable_design_variables = mutable,
    locked_design_variables = locked,
    lock_reasons = reasons[locked],
    feasibility_rule = "`raters_per_person <= n_rater`",
    note = note
  )
}

simulation_planning_constraints_note <- function(constraints) {
  if (is.list(constraints) && is.character(constraints$note) && length(constraints$note) > 0L) {
    note <- as.character(constraints$note[1])
    if (nzchar(note)) {
      return(note)
    }
  }
  character(0)
}

simulation_planning_schema <- function(sim_spec = NULL, facet_names = NULL) {
  descriptor <- if (is.null(facet_names)) {
    simulation_design_descriptor(sim_spec)
  } else {
    simulation_design_descriptor(list(facet_names = simulation_validate_output_facet_names(facet_names)))
  }
  scope <- if (is.null(facet_names)) {
    simulation_planning_scope(sim_spec)
  } else {
    simulation_planning_scope(facet_names = facet_names)
  }
  constraints <- simulation_planning_constraints(sim_spec)
  lock_lookup <- constraints$lock_reasons %||% stats::setNames(character(0), character(0))

  role_table <- descriptor |>
    dplyr::mutate(
      axis_class = c("person_count", "facet_level_count", "facet_level_count", "assignment_count"),
      depends_on_role = c(NA_character_, NA_character_, NA_character_, "rater"),
      mutable = .data$canonical %in% constraints$mutable_design_variables,
      locked = .data$canonical %in% constraints$locked_design_variables,
      lock_reason = unname(lock_lookup[.data$canonical])
    )
  role_table$lock_reason[is.na(role_table$lock_reason)] <- ""
  facet_manifest <- if (is.null(facet_names)) {
    simulation_facet_manifest(sim_spec)
  } else {
    simulation_facet_manifest(facet_names = facet_names)
  }
  structural_facet_table <- if (is.null(facet_names)) {
    simulation_structural_facet_table(sim_spec)
  } else {
    simulation_structural_facet_table(facet_names = facet_names)
  }
  structural_design_template <- if (is.null(facet_names)) {
    simulation_structural_design_template(sim_spec)
  } else {
    simulation_structural_design_template(facet_names = facet_names)
  }
  structural_design_schema <- if (is.null(facet_names)) {
    simulation_structural_design_schema(sim_spec)
  } else {
    simulation_structural_design_schema(facet_names = facet_names)
  }

  list(
    planner_contract = scope$planner_contract,
    planner_stage = scope$planner_stage,
    supports_arbitrary_facet_planning = scope$supports_arbitrary_facet_planning,
    supports_arbitrary_facet_estimation = scope$supports_arbitrary_facet_estimation,
    role_labels = scope$role_labels,
    design_variables = scope$design_variables,
    design_variable_aliases = scope$design_variable_aliases,
    facet_manifest = facet_manifest,
    structural_review_contract = scope$structural_review_contract,
    structural_review_stage = scope$structural_review_stage,
    structural_design_input_contract = scope$structural_design_input_contract,
    structural_facet_table = structural_facet_table,
    structural_design_template = structural_design_template,
    structural_design_schema = structural_design_schema,
    structural_design_preview = structural_design_schema$preview,
    structural_design_grid_semantics = structural_design_schema$grid_semantics,
    structural_design_grid_contract = structural_design_schema$grid_contract,
    structural_design_grid_bundle = structural_design_schema$grid_bundle,
    structural_design_grid_context = structural_design_schema$grid_context,
    structural_design_report_bundle = structural_design_schema$report_bundle,
    structural_design_report_summary = structural_design_schema$report_summary,
    structural_design_report_overview_table = structural_design_schema$report_overview_table,
    structural_design_report_catalog = structural_design_schema$report_catalog,
    structural_design_report_digest = structural_design_schema$report_digest,
    structural_design_report_surface_registry = structural_design_schema$report_surface_registry,
    structural_design_report_panel = structural_design_schema$report_panel,
    structural_design_report_operation = structural_design_schema$report_operation,
    structural_design_report_snapshot = structural_design_schema$report_snapshot,
    structural_design_report_brief = structural_design_schema$report_brief,
    structural_design_report_mode_registry = structural_design_schema$report_mode_registry,
    structural_design_report_consumer = structural_design_schema$report_consumer,
    structural_design_compact_review = structural_design_schema$review,
    structural_design_compact_review_summary = structural_design_schema$compact_review_summary,
    structural_design_compact_review_table = structural_design_schema$compact_review_table,
    structural_design_compact_review_plot = structural_design_schema$compact_review_plot,
    structural_design_review = structural_design_schema$design_review,
    structural_design_review_profile = structural_design_schema$design_review_profile,
    structural_design_review_load_balance = structural_design_schema$design_review_load_balance,
    structural_design_review_overview = structural_design_schema$design_review_overview,
    structural_design_review_load_balance_overview = structural_design_schema$design_review_load_balance_overview,
    structural_design_review_coverage = structural_design_schema$design_review_coverage,
    structural_design_review_coverage_overview = structural_design_schema$design_review_coverage_overview,
    structural_design_review_guardrails = structural_design_schema$design_review_guardrails,
    structural_design_review_guardrail_overview = structural_design_schema$design_review_guardrail_overview,
    structural_design_review_readiness = structural_design_schema$design_review_readiness,
    structural_design_review_readiness_overview = structural_design_schema$design_review_readiness_overview,
    structural_design_review_recommendation = structural_design_schema$design_review_recommendation,
    structural_design_review_recommendation_overview = structural_design_schema$design_review_recommendation_overview,
    role_table = role_table,
    feasibility_rules = constraints$feasibility_rule,
    note = paste(
      "Named-facet structural design metadata for person count, non-person",
      "facet counts, and assignments per person. These deterministic design",
      "summaries do not establish arbitrary-facet simulation support or",
      "parameter-recovery performance."
    )
  )
}

simulation_planning_schema_note <- function(schema) {
  if (is.list(schema) && is.character(schema$note) && length(schema$note) > 0L) {
    note <- as.character(schema$note[1])
    if (nzchar(note)) {
      return(note)
    }
  }
  character(0)
}

simulation_object_design_variable_aliases <- function(x) {
  aliases <- x$design_variable_aliases %||% x$settings$design_variable_aliases %||% NULL
  if (is.character(aliases) &&
      identical(sort(names(aliases)), sort(c("n_person", "n_rater", "n_criterion", "raters_per_person")))) {
    return(aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")])
  }
  sim_spec <- x$settings$sim_spec %||% NULL
  simulation_design_variable_aliases(sim_spec)
}

simulation_object_design_descriptor <- function(x) {
  descriptor <- x$design_descriptor %||% x$settings$design_descriptor %||% x$ademp$data_generating_mechanism$design_descriptor %||% NULL
  required_cols <- c("role", "canonical", "alias", "facet", "quantity", "description")
  if (is.data.frame(descriptor) && all(required_cols %in% names(descriptor))) {
    descriptor <- tibble::as_tibble(descriptor)
    return(descriptor[, required_cols])
  }
  sim_spec <- x$settings$sim_spec %||% NULL
  simulation_design_descriptor(sim_spec)
}

simulation_object_planning_scope <- function(x) {
  scope <- x$planning_scope %||% x$settings$planning_scope %||%
    x$ademp$data_generating_mechanism$planning_scope %||%
    x$sim_spec$planning_scope %||% NULL
  required_fields <- c(
    "planner_contract",
    "planner_stage",
    "supports_arbitrary_facet_planning",
    "supports_arbitrary_facet_estimation",
    "supported_non_person_roles",
    "supported_non_person_facet_count",
    "role_labels",
    "design_variables",
    "design_variable_aliases",
    "facet_manifest",
    "structural_review_contract",
    "structural_review_stage",
    "structural_design_input_contract",
    "note"
  )
  if (is.list(scope) &&
      all(required_fields %in% names(scope)) &&
      is.data.frame(scope$facet_manifest)) {
    return(scope[required_fields])
  }
  sim_spec <- x$settings$sim_spec %||% x$sim_spec %||% NULL
  simulation_planning_scope(sim_spec)
}

simulation_object_planning_constraints <- function(x) {
  constraints <- x$planning_constraints %||% x$settings$planning_constraints %||%
    x$ademp$data_generating_mechanism$planning_constraints %||%
    x$sim_spec$planning_constraints %||% NULL
  required_fields <- c(
    "mutable_design_variables",
    "locked_design_variables",
    "lock_reasons",
    "feasibility_rule",
    "note"
  )
  if (is.list(constraints) && all(required_fields %in% names(constraints))) {
    return(constraints[required_fields])
  }
  sim_spec <- x$settings$sim_spec %||% x$sim_spec %||% NULL
  simulation_planning_constraints(sim_spec)
}

simulation_object_planning_schema <- function(x) {
  schema <- x$planning_schema %||% x$settings$planning_schema %||%
    x$ademp$data_generating_mechanism$planning_schema %||%
    x$sim_spec$planning_schema %||% NULL
  required_fields <- c(
    "planner_contract",
    "planner_stage",
    "supports_arbitrary_facet_planning",
    "supports_arbitrary_facet_estimation",
    "role_labels",
    "design_variables",
    "design_variable_aliases",
    "facet_manifest",
    "structural_review_contract",
    "structural_review_stage",
    "structural_design_input_contract",
    "structural_facet_table",
    "structural_design_template",
    "structural_design_schema",
    "structural_design_preview",
    "structural_design_grid_semantics",
    "structural_design_grid_contract",
    "structural_design_grid_bundle",
    "structural_design_grid_context",
    "structural_design_report_bundle",
    "structural_design_report_summary",
    "structural_design_report_overview_table",
    "structural_design_report_catalog",
    "structural_design_report_digest",
    "structural_design_report_surface_registry",
    "structural_design_report_panel",
    "structural_design_report_operation",
    "structural_design_report_snapshot",
    "structural_design_report_brief",
    "structural_design_report_mode_registry",
    "structural_design_report_consumer",
    "structural_design_compact_review",
    "structural_design_compact_review_summary",
    "structural_design_compact_review_table",
    "structural_design_compact_review_plot",
    "structural_design_review",
    "structural_design_review_profile",
    "structural_design_review_load_balance",
    "structural_design_review_overview",
    "structural_design_review_load_balance_overview",
    "structural_design_review_coverage",
    "structural_design_review_coverage_overview",
    "structural_design_review_guardrails",
    "structural_design_review_guardrail_overview",
    "structural_design_review_readiness",
    "structural_design_review_readiness_overview",
    "structural_design_review_recommendation",
    "structural_design_review_recommendation_overview",
    "role_table",
    "feasibility_rules",
    "note"
  )
  if (is.list(schema) &&
      all(required_fields %in% names(schema)) &&
      is.data.frame(schema$role_table) &&
      is.data.frame(schema$facet_manifest) &&
      is.data.frame(schema$structural_facet_table) &&
      is.list(schema$structural_design_template) &&
      is.list(schema$structural_design_schema) &&
      is.list(schema$structural_design_preview) &&
      is.list(schema$structural_design_grid_semantics) &&
      is.list(schema$structural_design_grid_contract) &&
      is.list(schema$structural_design_grid_bundle) &&
      is.list(schema$structural_design_grid_context) &&
      is.list(schema$structural_design_report_bundle) &&
      is.list(schema$structural_design_report_summary) &&
      is.list(schema$structural_design_report_overview_table) &&
      is.list(schema$structural_design_report_catalog) &&
      is.list(schema$structural_design_report_digest) &&
      is.list(schema$structural_design_report_surface_registry) &&
      is.list(schema$structural_design_report_panel) &&
      is.list(schema$structural_design_report_operation) &&
      is.list(schema$structural_design_report_snapshot) &&
      is.list(schema$structural_design_report_brief) &&
      is.data.frame(schema$structural_design_report_mode_registry) &&
      is.list(schema$structural_design_report_consumer) &&
      is.list(schema$structural_design_compact_review) &&
      is.list(schema$structural_design_compact_review_summary) &&
      is.list(schema$structural_design_compact_review_table) &&
      is.list(schema$structural_design_compact_review_plot) &&
      is.list(schema$structural_design_review) &&
      is.list(schema$structural_design_review_profile) &&
      is.list(schema$structural_design_review_load_balance) &&
      is.list(schema$structural_design_review_overview) &&
      is.list(schema$structural_design_review_load_balance_overview) &&
      is.list(schema$structural_design_review_coverage) &&
      is.list(schema$structural_design_review_coverage_overview) &&
      is.list(schema$structural_design_review_guardrails) &&
      is.list(schema$structural_design_review_guardrail_overview) &&
      is.list(schema$structural_design_review_readiness) &&
      is.list(schema$structural_design_review_readiness_overview) &&
      is.list(schema$structural_design_review_recommendation) &&
      is.list(schema$structural_design_review_recommendation_overview)) {
    return(schema[required_fields])
  }
  sim_spec <- x$settings$sim_spec %||% x$sim_spec %||% NULL
  simulation_planning_schema(sim_spec)
}

simulation_design_variable_choices <- function(aliases, descriptor = NULL) {
  role_names <- if (is.data.frame(descriptor) && "role" %in% names(descriptor)) as.character(descriptor$role) else character(0)
  unique(c(names(aliases), unname(aliases), role_names))
}

simulation_design_variable_lookup <- function(aliases, descriptor = NULL) {
  canonical <- names(aliases)
  lookup <- c(stats::setNames(canonical, canonical), stats::setNames(canonical, unname(aliases)))
  if (is.data.frame(descriptor) && all(c("role", "canonical") %in% names(descriptor))) {
    role_lookup <- stats::setNames(as.character(descriptor$canonical), as.character(descriptor$role))
    lookup <- c(lookup, role_lookup)
  }
  lookup[!duplicated(names(lookup))]
}

simulation_resolve_design_variable <- function(value, aliases, arg_name, descriptor = NULL) {
  value <- as.character(value[1])
  lookup <- simulation_design_variable_lookup(aliases, descriptor = descriptor)
  resolved <- unname(lookup[[value]])
  if (!is.null(resolved) && nzchar(resolved)) {
    return(resolved)
  }
  stop(
    "`", arg_name, "` must be one of: ",
    paste(simulation_design_variable_choices(aliases, descriptor = descriptor), collapse = ", "),
    ".",
    call. = FALSE
  )
}

simulation_resolve_design_variable_vector <- function(values, aliases, descriptor = NULL) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    return(character(0))
  }
  lookup <- simulation_design_variable_lookup(aliases, descriptor = descriptor)
  resolved <- unname(lookup[values])
  unique(resolved[!is.na(resolved) & nzchar(resolved)])
}

simulation_design_variable_label <- function(value, aliases) {
  label <- unname(aliases[[value]])
  if (is.null(label) || is.na(label) || !nzchar(label)) {
    return(value)
  }
  label
}

simulation_parse_structural_facet_input <- function(facets,
                                                structural_facet_table,
                                                arg_name = "design$facets") {
  if (is.null(facets)) {
    return(list())
  }
  if (is.data.frame(facets)) {
    if (nrow(facets) != 1L) {
      stop("`", arg_name, "` must be a one-row data frame when supplied as a data frame.",
           call. = FALSE)
    }
    facets <- as.list(facets[1, , drop = FALSE])
  } else if (is.atomic(facets) && !is.list(facets)) {
    facets <- as.list(facets)
  } else if (!is.list(facets)) {
    stop(
      "`", arg_name, "` must be a named list, named vector, or one-row data frame keyed by future facet names.",
      call. = FALSE
    )
  }

  raw_names <- names(facets)
  valid_names <- as.character(structural_facet_table$structural_facet_key)
  if (is.null(raw_names) || any(!nzchar(raw_names))) {
    stop(
      "`", arg_name, "` must use names such as ",
      paste(valid_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  lookup <- stats::setNames(
    as.character(structural_facet_table$current_planning_count_variable),
    valid_names
  )
  resolved <- unname(lookup[raw_names])
  if (any(is.na(resolved) | !nzchar(resolved))) {
    bad <- unique(raw_names[is.na(resolved) | !nzchar(resolved)])
    stop(
      "`", arg_name, "` must use names such as ",
      paste(valid_names, collapse = ", "),
      ". Invalid names: ",
      paste(bad, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(resolved)) {
    dup <- unique(resolved[duplicated(resolved)])
    stop(
      "`", arg_name, "` supplies the same facet-count variable more than once: ",
      paste(dup, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  stats::setNames(unname(facets), resolved)
}

simulation_parse_structural_design <- function(design,
                                                  structural_design_schema,
                                                  arg_name = "design") {
  if (is.null(structural_design_schema) || !is.list(structural_design_schema)) {
    return(list())
  }
  design_schema <- simulation_coerce_structural_design_schema(structural_design_schema)

  parsed <- list()
  if ("facets" %in% names(design)) {
    parsed <- c(
      parsed,
      simulation_parse_structural_facet_input(
        facets = design[["facets"]],
        structural_facet_table = dplyr::transmute(
          design_schema$facet_axes,
          structural_facet_key = .data$input_key,
          current_planning_count_variable = .data$canonical_design_variable
        ),
        arg_name = paste0(arg_name, "$facets")
      )
    )
  }

  assignment_axis <- design_schema$assignment_axis
  assignment_key <- as.character(assignment_axis$design_input_key %||% "assignment")
  assignment_var <- as.character(
    assignment_axis$current_planning_count_variable %||% "raters_per_person"
  )
  if (assignment_key %in% names(design)) {
    parsed[[assignment_var]] <- design[[assignment_key]]
  }

  parsed
}

simulation_parse_design_input <- function(design,
                                          aliases,
                                          descriptor = NULL,
                                          structural_facet_table = NULL,
                                          structural_design_schema = NULL,
                                          arg_name = "design") {
  if (is.null(design)) {
    return(list())
  }

  if (is.data.frame(design)) {
    if (nrow(design) != 1L) {
      stop("`", arg_name, "` must be a one-row data frame when supplied as a data frame.",
           call. = FALSE)
    }
    design <- as.list(design[1, , drop = FALSE])
  } else if (is.atomic(design) && !is.list(design)) {
    design <- as.list(design)
  } else if (!is.list(design)) {
    stop(
      "`", arg_name, "` must be a named list, named vector, or one-row data frame. Valid names: ",
      paste(simulation_design_variable_choices(aliases, descriptor = descriptor), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  parsed_structural_design <- list()
  if ("facets" %in% names(design)) {
    if (!is.null(structural_design_schema) && is.list(structural_design_schema)) {
      parsed_structural_design <- simulation_parse_structural_design(
        design = design,
        structural_design_schema = structural_design_schema,
        arg_name = arg_name
      )
      assignment_key <- as.character(
        structural_design_schema$assignment_axis$design_input_key %||% "assignment"
      )
      design[["facets"]] <- NULL
      if (assignment_key %in% names(design)) {
        design[[assignment_key]] <- NULL
      }
    } else {
      if (is.null(structural_facet_table) || !is.data.frame(structural_facet_table)) {
        stop("`", arg_name, "$facets` is not available for this planning object.", call. = FALSE)
      }
      parsed_structural_design <- simulation_parse_structural_facet_input(
        facets = design[["facets"]],
        structural_facet_table = structural_facet_table,
        arg_name = paste0(arg_name, "$facets")
      )
      design[["facets"]] <- NULL
    }
  } else if (!is.null(structural_design_schema) && is.list(structural_design_schema)) {
    assignment_key <- as.character(
      structural_design_schema$assignment_axis$design_input_key %||% "assignment"
    )
    if (assignment_key %in% names(design)) {
      parsed_structural_design <- simulation_parse_structural_design(
        design = design,
        structural_design_schema = structural_design_schema,
        arg_name = arg_name
      )
      design[[assignment_key]] <- NULL
    }
  }

  raw_names <- names(design)
  if (length(design) == 0L) {
    raw_names <- character(0)
  } else if (is.null(raw_names) || any(!nzchar(raw_names))) {
    stop(
      "`", arg_name, "` must use names such as ",
      paste(simulation_design_variable_choices(aliases, descriptor = descriptor), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  parsed_top_level <- list()
  if (length(raw_names) > 0L) {
    canonical_names <- vapply(
      raw_names,
      simulation_resolve_design_variable,
      aliases = aliases,
      arg_name = arg_name,
      descriptor = descriptor,
      FUN.VALUE = character(1)
    )
    if (anyDuplicated(canonical_names)) {
      dup <- unique(canonical_names[duplicated(canonical_names)])
      stop(
        "`", arg_name, "` supplies the same design variable more than once: ",
        paste(dup, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    parsed_top_level <- stats::setNames(unname(design), canonical_names)
  }

  overlap <- intersect(names(parsed_top_level), names(parsed_structural_design))
  if (length(overlap) > 0L) {
    stop(
      "`", arg_name, "` supplies the same design variable through both top-level names and `",
      arg_name,
      "$facets`: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  c(parsed_top_level, parsed_structural_design)
}

simulation_resolve_design_counts <- function(sim_spec = NULL,
                                             n_person = NULL,
                                             n_rater = NULL,
                                             n_criterion = NULL,
                                             raters_per_person = NULL,
                                             design = NULL,
                                             defaults = NULL,
                                             design_arg = "design",
                                             explicit_scalar_names = NULL) {
  aliases <- simulation_design_variable_aliases(sim_spec)
  descriptor <- simulation_design_descriptor(sim_spec)
  parsed_design <- simulation_parse_design_input(
    design = design,
    aliases = aliases,
    descriptor = descriptor,
    structural_facet_table = simulation_structural_facet_table(sim_spec),
    structural_design_schema = simulation_structural_design_schema(sim_spec),
    arg_name = design_arg
  )

  explicit_scalars <- list(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person
  )
  if (is.null(explicit_scalar_names)) {
    explicit_scalars <- explicit_scalars[!vapply(explicit_scalars, is.null, logical(1))]
  } else {
    explicit_scalar_names <- intersect(names(explicit_scalars), explicit_scalar_names)
    explicit_scalars <- explicit_scalars[explicit_scalar_names]
    explicit_scalars <- explicit_scalars[!vapply(explicit_scalars, is.null, logical(1))]
  }
  overlap <- intersect(names(parsed_design), names(explicit_scalars))
  if (length(overlap) > 0L) {
    stop(
      "Do not supply the same design variable through both scalar arguments and `",
      design_arg,
      "`: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  values <- as.list(defaults %||% list())
  values[names(parsed_design)] <- parsed_design
  values[names(explicit_scalars)] <- explicit_scalars

  required <- c("n_person", "n_rater", "n_criterion")
  missing_required <- required[!required %in% names(values)]
  if (length(missing_required) > 0L) {
    stop(
      "Missing required design counts: ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!"raters_per_person" %in% names(values) || is.null(values$raters_per_person)) {
    values$raters_per_person <- values$n_rater
  }

  out <- tibble::tibble(
    n_person = simulation_validate_count(values$n_person, "n_person", min_value = 2L),
    n_rater = simulation_validate_count(values$n_rater, "n_rater", min_value = 2L),
    n_criterion = simulation_validate_count(values$n_criterion, "n_criterion", min_value = 2L),
    raters_per_person = simulation_validate_count(values$raters_per_person, "raters_per_person", min_value = 1L)
  )
  if (out$raters_per_person > out$n_rater) {
    stop("`raters_per_person` cannot exceed `n_rater`.", call. = FALSE)
  }
  out
}

simulation_validate_count_values <- function(x, arg_name, min_value = 1L) {
  values <- as.integer(x)
  if (length(values) < 1L || any(!is.finite(values)) || any(values < min_value)) {
    stop("`", arg_name, "` must contain integer values >= ", min_value, ".", call. = FALSE)
  }
  values
}

simulation_resolve_design_grid_values <- function(sim_spec = NULL,
                                                  n_person = NULL,
                                                  n_rater = NULL,
                                                  n_criterion = NULL,
                                                  raters_per_person = NULL,
                                                  design = NULL,
                                                  defaults = NULL,
                                                  design_arg = "design",
                                                  explicit_scalar_names = NULL) {
  aliases <- simulation_design_variable_aliases(sim_spec)
  descriptor <- simulation_design_descriptor(sim_spec)
  parsed_design <- simulation_parse_design_input(
    design = design,
    aliases = aliases,
    descriptor = descriptor,
    structural_facet_table = simulation_structural_facet_table(sim_spec),
    structural_design_schema = simulation_structural_design_schema(sim_spec),
    arg_name = design_arg
  )

  if (is.null(explicit_scalar_names)) {
    overlap_candidates <- c("n_person", "n_rater", "n_criterion", "raters_per_person")
  } else {
    overlap_candidates <- intersect(
      c("n_person", "n_rater", "n_criterion", "raters_per_person"),
      explicit_scalar_names
    )
  }
  overlap <- intersect(names(parsed_design), overlap_candidates)
  if (length(overlap) > 0L) {
    stop(
      "Do not supply the same design variable through both scalar arguments and `",
      design_arg,
      "`: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  values <- as.list(defaults %||% list(
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person
  ))
  values[names(parsed_design)] <- parsed_design

  required <- c("n_person", "n_rater", "n_criterion")
  missing_required <- required[!required %in% names(values)]
  if (length(missing_required) > 0L) {
    stop(
      "Missing required design values: ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!"raters_per_person" %in% names(values) || is.null(values$raters_per_person)) {
    values$raters_per_person <- values$n_rater
  }

  list(
    n_person = simulation_validate_count_values(values$n_person, "n_person", min_value = 2L),
    n_rater = simulation_validate_count_values(values$n_rater, "n_rater", min_value = 2L),
    n_criterion = simulation_validate_count_values(values$n_criterion, "n_criterion", min_value = 2L),
    raters_per_person = simulation_validate_count_values(values$raters_per_person, "raters_per_person", min_value = 1L)
  )
}

simulation_design_canonical_variables <- function(descriptor = NULL) {
  if (is.data.frame(descriptor) && "canonical" %in% names(descriptor)) {
    vals <- as.character(descriptor$canonical)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) > 0L) {
      return(vals)
    }
  }
  c("n_person", "n_rater", "n_criterion", "raters_per_person")
}

simulation_design_group_variables <- function(descriptor = NULL, include_design_id = TRUE) {
  vars <- simulation_design_canonical_variables(descriptor)
  if (isTRUE(include_design_id)) {
    c("design_id", vars)
  } else {
    vars
  }
}

simulation_append_design_alias_columns <- function(tbl, aliases) {
  tbl <- tibble::as_tibble(tbl)
  if (!is.character(aliases) || length(aliases) == 0L) {
    return(tbl)
  }
  for (canonical in names(aliases)) {
    alias <- as.character(aliases[[canonical]])
    if (!canonical %in% names(tbl)) next
    if (is.na(alias) || !nzchar(alias) || identical(alias, canonical) || alias %in% names(tbl)) next
    tbl[[alias]] <- tbl[[canonical]]
  }
  tbl
}

simulation_build_design_grid <- function(n_person,
                                         n_rater,
                                         n_criterion,
                                         raters_per_person,
                                         design = NULL,
                                         sim_spec = NULL,
                                         id_prefix = "D",
                                         explicit_scalar_names = NULL) {
  aliases <- simulation_design_variable_aliases(sim_spec)
  descriptor <- simulation_design_descriptor(sim_spec)
  canonical_order <- as.character(descriptor$canonical)
  design_values <- simulation_resolve_design_grid_values(
    sim_spec = sim_spec,
    n_person = n_person,
    n_rater = n_rater,
    n_criterion = n_criterion,
    raters_per_person = raters_per_person,
    design = design,
    defaults = list(
      n_person = n_person,
      n_rater = n_rater,
      n_criterion = n_criterion,
      raters_per_person = raters_per_person
    ),
    design_arg = "design",
    explicit_scalar_names = explicit_scalar_names
  )

  design_grid <- expand.grid(
    n_person = design_values$n_person,
    n_rater = design_values$n_rater,
    n_criterion = design_values$n_criterion,
    raters_per_person = design_values$raters_per_person,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design_grid <- design_grid[design_grid$raters_per_person <= design_grid$n_rater, , drop = FALSE]
  if (nrow(design_grid) == 0) {
    msg <- "No valid design rows remain after enforcing `raters_per_person <= n_rater`."
    public_rule <- paste0("`", unname(aliases[["raters_per_person"]]), " <= ", unname(aliases[["n_rater"]]), "`")
    if (!identical(public_rule, "`raters_per_person <= n_rater`")) {
      msg <- paste0(msg, " Public aliases: ", public_rule, ".")
    }
    stop(msg, call. = FALSE)
  }
  design_grid$design_id <- sprintf("%s%02d", as.character(id_prefix[1]), seq_len(nrow(design_grid)))
  design_grid <- tibble::as_tibble(design_grid[, c("design_id", canonical_order), drop = FALSE])

  list(
    canonical = design_grid,
    public = simulation_append_design_alias_columns(design_grid, aliases),
    aliases = aliases,
    descriptor = descriptor
  )
}
