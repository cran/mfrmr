# Conditional fitted probabilities at the saved corrected calibration and
# reprofiled abilities. No posterior integration or new fitting is performed.
mfrm_jml_response_input <- function(fit) {
  if (!mfrm_has_jml_adjustment(fit) || !isTRUE(fit$jml_adjustment$point$available))
    stop("Corrected-JML response diagnostics require an unambiguous saved point solution.",call.=FALSE)
  cfg <- fit$config; d <- fit$prep$data; columns <- cfg$source_columns
  data <- as.data.frame(d[c("Person",cfg$facet_names,"Score")])
  names(data) <- c(columns$person,cfg$facet_names,columns$score)
  data[setdiff(names(data),columns$score)] <- lapply(data[setdiff(names(data),columns$score)],as.character)
  rownames(data) <- NULL
  theta <- fit$facets$person$Estimate[match(as.character(d$Person),fit$facets$person$Person)]
  offset <- rep(0,nrow(d))
  for (f in cfg$facet_names) {
    tab <- fit$facets$others[fit$facets$others$Facet==f,,drop=FALSE]
    if (anyDuplicated(tab$Level)) stop("Saved facet levels are duplicated.",call.=FALSE)
    offset <- offset+tab$Estimate[match(as.character(d[[f]]),tab$Level)]
  }
  owners <- as.character(fit$slopes$SlopeFacet); slopes <- fit$slopes$Estimate
  ns <- cfg$n_cat-1L
  steps <- t(vapply(owners,function(id) {
    tab <- fit$steps[fit$steps$StepFacet==id,,drop=FALSE]
    order <- match(paste0("Step_",seq_len(ns)),tab$Step)
    if (nrow(tab)!=ns || anyDuplicated(tab$Step) || anyNA(order))
      stop("Saved category steps do not match the declared ladder.",call.=FALSE)
    tab <- tab[order,,drop=FALSE]
    c(0,cumsum(tab$Estimate))
  },numeric(ns+1L)))
  owner <- match(as.character(d[[cfg$step_facet]]),owners)
  if (anyDuplicated(owners) || anyNA(owner) || anyNA(theta) || any(!is.finite(offset)) ||
      any(!is.finite(steps)) || any(!is.finite(slopes) | slopes<=0))
    stop("Saved corrected-JML locations, steps, slopes or Person profiles are unavailable.",call.=FALSE)
  list(data=data,columns=columns,score_levels=seq(cfg$rating_min,cfg$rating_max),
    eta=theta-offset,theta=theta,step_cum=steps,owner=owner,slopes=slopes)
}

mfrm_jml_response_diagnostics <- function(fit, rows, group_by) {
  input <- mfrm_jml_response_input(fit); data <- input$data
  rows <- rows %||% seq_len(nrow(data))
  if (!is.numeric(rows) || is.complex(rows) || !is.null(dim(rows)) || !length(rows) ||
      anyNA(rows) || any(rows!=floor(rows)) || any(rows<1 | rows>nrow(data)) || anyDuplicated(rows))
    stop("`rows` must contain distinct original input row numbers.",call.=FALSE)
  identifiers <- setdiff(names(data),input$columns$score)
  group_by <- group_by %||% c(input$columns$facets,input$columns$person)
  if (!is.character(group_by) || !length(group_by) || anyNA(group_by) ||
      anyDuplicated(group_by) || !all(group_by %in% identifiers))
    stop("Choose distinct fitted identifier columns for group_by.",call.=FALSE)
  categories <- input$score_levels; nc <- length(categories)
  p <- matrix(NA_real_,length(rows),nc,dimnames=list(as.character(rows),as.character(categories)))
  finite <- which(is.finite(input$theta[rows])); extreme <- which(is.infinite(input$theta[rows]))
  if (length(finite)) p[finite,] <- category_prob_gpcm(input$eta[rows[finite]],
    input$step_cum,input$owner[rows[finite]],input$slopes)
  if (length(extreme)) {
    p[extreme,] <- 0
    p[cbind(extreme,ifelse(input$theta[rows[extreme]]<0,1L,nc))] <- 1
  }
  probability_ok <- apply(p,1L,function(z) all(is.finite(z)) && all(z>=0) && abs(sum(z)-1)<1e-10)
  p[!probability_ok,] <- NA_real_
  mu <- as.vector(p %*% categories)
  variance <- rowSums(p*(matrix(categories,length(rows),nc,byrow=TRUE)-mu)^2)
  residual <- data[[input$columns$score]][rows]-mu
  standardized <- rep(NA_real_,length(rows))
  usable <- probability_ok & is.finite(variance) & variance>0
  standardized[usable] <- residual[usable]/sqrt(variance[usable])
  usable <- usable & is.finite(standardized^2) & is.finite(residual^2)
  standardized[!usable] <- NA_real_
  zero <- probability_ok & is.finite(variance) & variance==0
  table <- data.frame(InputRow=rows,Score=data[[input$columns$score]][rows],
    ExpectedScore=mu,PredictiveVariance=variance,Residual=residual,
    SquaredResidual=residual^2,StandardizedResidual=standardized,
    ProbabilityAvailable=probability_ok,
    Status=ifelse(usable,"available_conditional",ifelse(zero,"zero_variance","unavailable")),
    Reason=ifelse(usable,"",ifelse(zero,
      "Conditional variance is zero; standardized residuals are undefined. Probabilities and raw residuals are retained.",
      "Conditional probabilities or standardized residuals are numerically unavailable.")))
  table$Reason[extreme] <- paste("Extreme Person profile gives a limiting point mass;",
    "standardized residuals are undefined. Probabilities and raw residuals are retained.")
  structure(list(rows=table,probabilities=p,
    measures=mfrm_response_group_measures(table,data,rows,group_by),
    source=mfrm_response_source(fit),source_data=data,source_observed=data,
    settings=list(model="Corrected JML GPCM",group_by=group_by,
      probability_method="corrected_jml_plugin",correction_order=fit$jml_adjustment$estimator$order,
      calibration_uncertainty=FALSE,person_uncertainty=FALSE,
      target="Same-data conditional probabilities at corrected calibration and reprofiled Person estimates",
      integration="No integration; saved calibration and Person profiles held fixed",
      limitation=paste("Descriptive only; no established expectation-one reference, cutoffs, ZSTD, p-values or rater-quality classification.",
        "Residual bias may remain; calibration and Person estimation uncertainty are excluded.",
        "Outfit requires every selected standardized residual. Infit can retain zero-variance rows when all raw residuals are known and total variance is positive."),
      selection="Selected original observed rows only; Person profiles retain all their fitted observations")),
    class="mfrm_response_diagnostics")
}
