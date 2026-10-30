# Keep optional backend failures testable without modifying another namespace.
.mfrm_boot_mer <- function(...) lme4::bootMer(...)
.mfrm_icc_lmer <- function(...) lme4::lmer(...)
.mfrm_mirt_coef <- function(...) mirt::coef(...)
.mfrm_plot_igraph <- function(...) igraph::plot.igraph(...)
