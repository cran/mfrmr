# Internal text-file checksum helper.
#
# Git and archive tools can represent text line endings differently across
# platforms. Normalize CRLF and lone CR line endings to LF before calculating
# the checksum so the result reflects text content rather than checkout policy.
mfrmr_normalized_text_md5 <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be one non-empty file path.", call. = FALSE)
  }
  info <- file.info(path)
  if (is.na(info$size[[1L]]) || isTRUE(info$isdir[[1L]])) {
    stop("`path` must identify an existing file.", call. = FALSE)
  }
  if (info$size[[1L]] > .Machine$integer.max) {
    stop("`path` is too large for the text checksum helper.", call. = FALSE)
  }

  bytes <- readBin(path, what = "raw", n = as.integer(info$size[[1L]]))
  values <- as.integer(bytes)
  if (length(values) > 0L) {
    cr_before_lf <- values == 13L & c(values[-1L] == 10L, FALSE)
    values <- values[!cr_before_lf]
    values[values == 13L] <- 10L
    bytes <- as.raw(values)
  }

  normalized_path <- tempfile("mfrmr-normalized-text-")
  on.exit(unlink(normalized_path), add = TRUE)
  con <- file(normalized_path, open = "wb")
  tryCatch(
    writeBin(bytes, con),
    finally = close(con)
  )
  unname(tools::md5sum(normalized_path))
}
