public_mfrm_method_label <- function(method, default = NA_character_) {
  default <- as.character(default)[1] %||% NA_character_
  if (is.null(method)) {
    return(default)
  }

  method_chr <- as.character(method)
  if (length(method_chr) == 0L) {
    return(default)
  }

  method_chr <- toupper(trimws(method_chr))
  missing_idx <- is.na(method_chr) | !nzchar(method_chr)
  method_chr[method_chr == "JMLE"] <- "JML"
  method_chr[missing_idx] <- default
  method_chr
}

resolve_public_mfrm_method <- function(summary_method = NULL,
                                       method_input = NULL,
                                       method_used = NULL,
                                       default = NA_character_) {
  candidates <- c(as.character(summary_method), as.character(method_input), as.character(method_used))
  candidates <- trimws(candidates)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 0L) {
    return(public_mfrm_method_label(default, default = default))
  }
  public_mfrm_method_label(candidates[1], default = default)
}
