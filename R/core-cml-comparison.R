mfrmr_cml_support_design <- function(problem) {
    p <- problem$parameters
    if (as.double(p) * p > 2e+06)
        return(NULL)
    rows <- list()
    for (g in problem$groups) {
        totals <- which(g$counts > 0) - 1L
        m <- length(g$features)
        K <- nrow(g$features[[1L]]) - 1L
        if (as.double(max(totals) + 1L) * p > 2e+06)
            return(NULL)
        previous <- matrix(0, 1L, p)
        for (j in seq_len(m)) {
            current <- matrix(NA_real_, min(j * K, max(totals)) + 1L, p)
            for (s in seq_len(nrow(current)) - 1L) {
                if (!any(totals >= s & totals <= s + (m - j) * K))
                  next
                for (k in seq.int(max(0L, s - nrow(previous) + 1L), min(K, s))) {
                  prefix <- previous[s - k + 1L, ]
                  if (anyNA(prefix))
                    next
                  candidate <- prefix + g$features[[j]][k + 1L, ]
                  if (anyNA(current[s + 1L, ]))
                    current[s + 1L, ] <- candidate
                  else {
                    delta <- candidate - current[s + 1L, ]
                    if (any(delta != 0)) {
                      if ((as.double(length(rows)) + 1) * p > 2e+06)
                        return(NULL)
                      rows[[length(rows) + 1L]] <- delta
                    }
                  }
                }
            }
            previous <- current
        }
    }
    if (length(rows))
        unique(do.call(rbind, rows))
    else matrix(0, 1L, p)
}
