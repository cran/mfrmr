mfrm_write_export_zip <- function(path, files, overwrite = FALSE) {
  if (dir.exists(path)) stop("Output path is a directory: ", path, call. = FALSE)
  if (file.exists(path) && !overwrite)
    stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
  path <- file.path(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), basename(path))
  files <- unique(normalizePath(files, winslash = "/", mustWork = TRUE))
  if (path %in% files)
    stop("The ZIP destination must differ from every exported file.", call. = FALSE)
  if (anyDuplicated(basename(files)))
    stop("ZIP inputs must have distinct file names.", call. = FALSE)

  # Updating an existing ZIP retains obsolete members; build a fresh archive.
  pending <- tempfile(".mfrmr-export-", tmpdir = dirname(path), fileext = ".zip")
  on.exit(unlink(pending), add = TRUE)
  tryCatch({
    status <- zip(pending, files, flags = "-q9X", extras = "-j")
    if (length(status) != 1L || is.na(status) || status != 0)
      stop("ZIP command failed (exit status: ", paste(status, collapse = ", "), ").")
    if (!file.exists(pending)) stop("ZIP command did not create an archive.")
    entries <- utils::unzip(pending, list = TRUE)
    expected <- basename(files)
    if (!identical(sort(entries$Name), sort(expected)) ||
        !identical(as.numeric(entries$Length[match(expected, entries$Name)]),
                   as.numeric(file.info(files)$size)))
      stop("ZIP contents do not match the files written by this export.")
    if (!isTRUE(suppressWarnings(file.rename(pending, path))))
      stop("Could not move the completed ZIP archive to its destination.")
    list(path = path, note = NULL)
  }, error = function(e) list(path = NULL, note = conditionMessage(e)),
     warning = function(w) list(path = NULL, note = conditionMessage(w)))
}
