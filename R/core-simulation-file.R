# File reuse is deliberately stricter than loading an arbitrary saved result.
simulation_file_controls <- function() {
  list(mfrmr.use_cpp11_backend = getOption("mfrmr.use_cpp11_backend", TRUE),
       mfrmr.max_information_bytes = getOption("mfrmr.max_information_bytes", 256 * 1024^2),
       mfrmr.misfit_lower = getOption("mfrmr.misfit_lower", 0.5),
       mfrmr.misfit_upper = getOption("mfrmr.misfit_upper", 1.5),
       contrasts = getOption("contrasts"), na.action = getOption("na.action"))
}

simulation_file_inputs <- function(fun, frame, call) {
  args <- setdiff(names(formals(fun)), c("file", "file_refit", "progress", "parallel"))
  list(arguments = mget(args, envir = frame, inherits = FALSE),
       supplied = sort(intersect(names(call)[-1],
         c("model", "n_person", "n_rater", "n_criterion", "raters_per_person"))))
}

simulation_file_engine <- function() {
  ns <- asNamespace("mfrmr")
  names <- sort(ls(ns, all.names = TRUE))
  code <- lapply(names, function(name) {
    obj <- get(name, envir = ns, inherits = FALSE)
    if (!is.function(obj) || !identical(environment(obj), ns)) return(NULL)
    list(name = name, formals = deparse(formals(obj)), body = deparse(body(obj)))
  })
  dll <- getLoadedDLLs()[["mfrmr"]]
  native <- if (is.null(dll)) NA_character_ else unname(tools::md5sum(dll[["path"]]))
  desc <- utils::packageDescription("mfrmr")
  packages <- trimws(gsub("\\s*\\(.*\\)", "", strsplit(desc$Imports, ",")[[1]]))
  packages <- sort(unique(c("mfrmr", packages, "RTMB", "nleqslv", "lpSolve")))
  versions <- vapply(packages, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  list(code = mfrm_checkpoint_fingerprint(Filter(Negate(is.null), code)),
       native = native, packages = versions, R = R.version.string,
       platform = R.version$platform, numerical_libraries = extSoftVersion(),
       controls = simulation_file_controls())
}

simulation_file_open <- function(file, file_refit, api, inputs, cells, seed) {
  file_refit <- match.arg(file_refit, c("never", "on_change", "always"))
  if (is.null(file)) {
    return(list(path = NULL, seeds = with_preserved_rng_seed(seed,
      sample.int(.Machine$integer.max, size = cells, replace = FALSE))))
  }
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be NULL or one non-empty path.", call. = FALSE)
  }
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) file <- paste0(file, ".rds")
  file <- path.expand(file)
  if (!dir.exists(dirname(file))) {
    stop("Simulation file directory does not exist: ", dirname(file), call. = FALSE)
  }
  file <- file.path(normalizePath(dirname(file), mustWork = TRUE), basename(file))
  lock <- paste0(file, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Simulation file is locked: ", lock,
      ". After a terminated process, remove this directory only after confirming no writer is active.",
      call. = FALSE)
  }
  opened <- FALSE
  on.exit(if (!opened) unlink(lock, recursive = TRUE), add = TRUE)
  identity <- list(schema = 1L, api = api,
    inputs = mfrm_checkpoint_fingerprint(inputs, canonical = TRUE),
    rng_kind = RNGkind(), engine = simulation_file_engine())
  saved <- if (file.exists(file) && file_refit != "always") {
    tryCatch(readRDS(file), error = function(e) stop(
      "Cannot read simulation file: ", conditionMessage(e),
      ". Use a new path or `file_refit = 'always'` to replace it.", call. = FALSE))
  } else NULL
  if (!is.null(saved)) {
    metadata <- if (is.list(saved)) saved$reproducibility else NULL
    if (!is.list(metadata) || !identical(metadata$identity, identity)) {
      if (file_refit == "never") {
        stop("Simulation file does not match the inputs, RNG kind or implementation/environment. ",
          "Use a new path or `file_refit = 'on_change'` to recompute.", call. = FALSE)
      }
      message("Simulation settings or implementation/environment changed; starting a new run.")
      saved <- NULL
    }
  }
  if (!is.null(saved)) {
    metadata <- saved$reproducibility
    if (!is.integer(metadata$seeds) || length(metadata$seeds) != cells ||
        anyNA(metadata$seeds) || any(metadata$seeds < 1L) ||
        !is.logical(metadata$complete) || length(metadata$complete) != 1L ||
        is.na(metadata$complete) ||
        (!metadata$complete && (!inherits(saved, "mfrm_simulation_checkpoint") ||
          !is.list(saved$cells) || length(saved$cells) != cells))) {
      stop("Invalid simulation checkpoint contents; use a new path or `file_refit = 'always'`.",
        call. = FALSE)
    }
    expected_class <- switch(api, simulate_mfrm_recovery = "mfrm_recovery_simulation",
      evaluate_mfrm_design = "mfrm_design_evaluation", NA_character_)
    if (metadata$complete && (is.na(expected_class) || !inherits(saved, expected_class))) {
      stop("Simulation file does not contain the expected result class.", call. = FALSE)
    }
    state <- list(path = file, lock = lock, seeds = metadata$seeds,
      metadata = metadata, cells = if (metadata$complete) vector("list", cells) else saved$cells,
      result = if (metadata$complete) saved else NULL)
    message(if (metadata$complete) "Loaded completed simulation: " else
      "Resuming saved simulation: ", file)
  } else {
    seeds <- with_preserved_rng_seed(seed,
      sample.int(.Machine$integer.max, size = cells, replace = FALSE))
    state <- list(path = file, lock = lock, seeds = seeds,
      metadata = list(identity = identity, inputs = inputs, seeds = seeds,
        complete = FALSE), cells = vector("list", cells), result = NULL)
    simulation_file_write(state)
  }
  opened <- TRUE
  state
}

simulation_file_write <- function(state) {
  if (!is.null(state$path)) {
    mfrm_atomic_save_checkpoint(structure(list(reproducibility = state$metadata,
      cells = state$cells), class = "mfrm_simulation_checkpoint"), state$path)
  }
  invisible(NULL)
}

simulation_file_finish <- function(state, result) {
  if (!is.null(state$path)) {
    state$metadata$complete <- TRUE
    result$reproducibility <- state$metadata
    mfrm_atomic_save_checkpoint(result, state$path)
  }
  result
}
