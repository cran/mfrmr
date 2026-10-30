# Keep Person-level covariates and response assignments on a declared sampling law.
simulation_person_template_plan <- function(sim_spec) {
  population <- sim_spec$population %||% list()
  active <- isTRUE(population$active)
  mode <- match.arg(sim_spec$population_assignment %||% "independent",
    c("matched", "independent"))
  roster <- switch(sim_spec$assignment %||% "crossed",
    skeleton = sim_spec$design_skeleton, resampled = sim_spec$assignment_profiles, NULL)
  people <- if (active) population$covariate_template else NULL
  person_id <- population$person_id %||% "TemplatePerson"
  groups <- function(data, label) {
    if (is.null(data) || !"Group" %in% names(data)) return(NULL)
    value <- as.character(data$Group)
    if (anyNA(value) || any(!nzchar(value))) {
      stop(label, " must contain non-missing, non-empty Group labels.", call. = FALSE)
    }
    value
  }
  roster_groups <- groups(roster, "Assignment templates")
  population_groups <- groups(people, "Population templates")
  if (!is.null(roster_groups)) {
    mapping <- unique(data.frame(Person = roster$TemplatePerson, Group = roster_groups))
    if (anyDuplicated(mapping$Person)) stop("Assignment templates must assign one Group per TemplatePerson.", call. = FALSE)
  }
  supplied_groups <- sim_spec$group_levels
  if (!is.null(supplied_groups) &&
      any(!c(roster_groups, population_groups) %in% supplied_groups)) {
    stop("Template Group labels must be members of `group_levels`.", call. = FALSE)
  }
  reuse <- FALSE
  if (!is.null(roster) && "TemplatePersonReuse" %in% names(roster)) {
    reuse_values <- unique(simulation_sparse_design_active(roster$TemplatePersonReuse))
    if (length(reuse_values) != 1L) {
      stop("`TemplatePersonReuse` must be uniformly TRUE or FALSE across the skeleton.", call. = FALSE)
    }
    reuse <- isTRUE(reuse_values)
  }
  matched <- active && !is.null(roster) && identical(mode, "matched")
  templates <- if (matched) unique(as.character(roster$TemplatePerson)) else NULL
  if (matched) {
    pop_ids <- as.character(people[[person_id]])
    if (anyDuplicated(pop_ids) || !setequal(templates, pop_ids)) {
      stop("Matched population and assignment templates must cover the same unique `TemplatePerson` IDs. Use `population_assignment = \"independent\"` only for an intentionally independent design.", call. = FALSE)
    }
    if (!is.null(roster_groups) && !is.null(population_groups) &&
        any(roster_groups != population_groups[match(roster$TemplatePerson, pop_ids)])) {
      stop("Population and assignment Group labels must agree for each TemplatePerson.", call. = FALSE)
    }
    if (reuse && length(templates) != sim_spec$n_person) {
      stop("`TemplatePersonReuse = TRUE` requires the template count to match `n_person`.", call. = FALSE)
    }
  } else if (active && !is.null(roster_groups) && !is.null(population_groups)) {
    stop("Independent population and assignment templates cannot both define Group. Choose matched sampling or rename a distinct population covariate.", call. = FALSE)
  }
  list(mode = if (!active) "none" else if (matched) "matched" else "independent",
    templates = templates, reuse = reuse)
}
