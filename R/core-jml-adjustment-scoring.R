# Fixed-calibration scoring has its own admission rule: the adjusted equation,
# not the unadjusted likelihood, must have a locally resolved root.
mfrm_jml_scoring_components <- function(fit) {
  z <- fit$jml_adjustment; cfg <- fit$config
  fail <- function() stop("Corrected-JML scoring requires an intact, unambiguous adjusted-equation solution; refit the source if its data or calibration have changed.", call.=FALSE)
  if (!mfrm_has_jml_adjustment(fit) || !identical(cfg$model,"GPCM") ||
      !identical(cfg$method,"JML") || !identical(cfg$slope_facet,cfg$step_facet) ||
      !isTRUE(z$point$available) || !z$point$status %in% c("consistent_roots","root_with_unresolved_starts") ||
      !identical(z$estimator$name,"finite_MLE_plugin_profile_score_adjustment") ||
      !identical(z$estimator$order_selection,"explicit") ||
      length(cfg$interaction_specs) || isTRUE(fit$population$active) ||
      any(is.finite(cfg$theta_spec$anchors)) || length(cfg$theta_spec$group_values) ||
      nrow(mfrmr_calibration_extract_anchors(fit))>0L ||
      is.null(fit$prep$data$Weight) || anyNA(fit$prep$data$Weight) || any(fit$prep$data$Weight != 1) ||
      !identical(z$policy$version,"explicit_order_v1")) fail()
  problem <- mfrm_jml_adjustment_problem(fit$prep$data,"Person",cfg$facet_names,
    "score_k",cfg$step_facet,cfg$n_cat-1L)
  if (!identical(problem$parameters,z$parameters) || !identical(problem$levels,z$levels) ||
      !identical(problem$specification,z$specification) ||
      !identical(cfg$facet_levels,problem$levels[cfg$facet_names]) ||
      !identical(fit$prep$levels[cfg$facet_names],problem$levels[cfg$facet_names]) ||
      !identical(as.numeric(fit$prep$data$Score-cfg$rating_min),as.numeric(fit$prep$data$score_k)) ||
      cfg$rating_min != fit$prep$rating_min || cfg$rating_max != fit$prep$rating_max) fail()
  order <- z$estimator$order
  if (length(order)!=1L || !is.numeric(order) || !is.finite(order) || order<1 || order!=floor(order)) fail()
  fn <- function(b) problem$mean_score(b,order)
  value <- fn(z$point$beta)
  derivative <- mfrmr_numeric_transformation_jacobian(fn,z$point$beta,relative_step=5e-5)
  if (!derivative$valid) fail()
  singular <- svd(derivative$jacobian,nu=0L,nv=0L)$d
  rank_ok <- length(singular)==length(value) &&
    min(singular)>max(singular)*length(value)*.Machine$double.eps
  if (!rank_ok) fail()
  checks <- c(equation_residual=max(abs(value)),
    remaining_step=max(abs(solve(derivative$jacobian,value))))
  if (any(!is.finite(checks)) || checks[1]>1e-7 || checks[2]>1e-5) fail()
  # Covariance availability does not control fixed-point scoring.
  z$covariance$available <- FALSE
  tab <- mfrm_jml_adjustment_tables(z)
  matching <- function(saved, expected, keys) {
    if (!is.data.frame(saved) || nrow(saved)!=nrow(expected)) return(FALSE)
    canonical <- function(x) {
      out <- as.data.frame(x[c(keys,"Estimate")],stringsAsFactors=FALSE)
      out <- out[do.call(base::order,out[keys]),,drop=FALSE]; rownames(out) <- NULL
      out
    }
    !anyDuplicated(saved[keys]) &&
      isTRUE(all.equal(canonical(saved),canonical(expected),tolerance=1e-12))
  }
  if (!matching(fit$facets$others,tab$facets$others,c("Facet","Level")) ||
      !matching(fit$steps,tab$steps,c("StepFacet","Step")) ||
      !matching(fit$slopes,tab$slopes,"SlopeFacet")) fail()
  config <- build_estimation_config(fit$prep,"GPCM","JML",cfg$step_facet,cfg$slope_facet,
    weight_col=NULL,facet_signs=build_facet_signs(cfg$facet_names)$signs,positive_facets=character(),
    noncenter_facet="Person",dummy_facets=character(),anchor_df=NULL,group_anchor_df=NULL)$config
  # Preserve estimator identity; no optimizer vector or ordinary MLE is made.
  config$jml_adjustment <- TRUE
  config$source_columns <- cfg$source_columns
  config$posterior_basis <- "post_hoc_standard_normal"
  config$estimation_control <- cfg$estimation_control
  params <- list(facets=setNames(lapply(cfg$facet_names,function(f)
      tab$facets$others$Estimate[tab$facets$others$Facet==f]),cfg$facet_names),
    steps_mat=matrix(tab$steps$Estimate,nrow=length(cfg$facet_levels[[cfg$step_facet]]),byrow=TRUE),
    slopes=tab$slopes$Estimate,log_slopes=log(tab$slopes$Estimate),interactions=list())
  local <- list(eligible=TRUE,basis="adjusted_equation_root_v1",
    review="The corrected equation and its full Jacobian support this local root; global uniqueness is not established.",
    caution="Experimental corrected JML: residual calibration bias may remain. EAP uses a separate normal reference prior; calibration uncertainty is excluded.",
    correction_order=as.integer(order),sampling=z$estimator$sampling,
    point_status=z$point$status,checks=checks)
  source <- list(policy_basis="corrected_jml_reference_eap_v1",status="conditional",
    local_calibration_review=local,
    source_audit_states=c(identification="local_equation_rank",boundary="not_evaluated",numerical="ready"),
    inference_ready=FALSE,ready=TRUE,reason_codes=character(),fit_readiness=data.frame(),
    readiness_contract_version=mfrmr_readiness_contract_version())
  if (!mfrm_jml_scoring_evidence_valid(source[1:5])) fail()
  list(config=config,params=params,source=source)
}

mfrm_jml_scoring_evidence_valid <- function(evidence) {
  tryCatch({
    z <- evidence$local_calibration_review; checks <- z$checks
    isTRUE(identical(names(evidence),c("policy_basis","status","local_calibration_review","source_audit_states","inference_ready")) &&
      identical(evidence$policy_basis,"corrected_jml_reference_eap_v1") &&
      identical(evidence$status,"conditional") && identical(evidence$inference_ready,FALSE) &&
      identical(evidence$source_audit_states,c(identification="local_equation_rank",boundary="not_evaluated",numerical="ready")) &&
      identical(names(z),c("eligible","basis","review","caution","correction_order","sampling","point_status","checks")) &&
      identical(z$eligible,TRUE) && identical(z$basis,"adjusted_equation_root_v1") &&
      is.character(z$review) && length(z$review)==1L && nzchar(z$review) &&
      is.character(z$caution) && length(z$caution)==1L && nzchar(z$caution) &&
      is.integer(z$correction_order) && length(z$correction_order)==1L && z$correction_order>0L &&
      length(z$sampling)==1L && z$sampling %in% c("fixed_rosters","random_rosters") &&
      length(z$point_status)==1L && z$point_status %in% c("consistent_roots","root_with_unresolved_starts") &&
      is.numeric(checks) && identical(names(checks),c("equation_residual","remaining_step")) &&
      all(is.finite(checks) & checks>=0 & checks<=c(1e-7,1e-5)))
  },error=function(e) FALSE)
}

mfrm_jml_scoring_note <- function(evidence) paste0(
  "Experimental corrected JML calibration (order ",evidence$local_calibration_review$correction_order,
  "). Scores are post-hoc EAP, not corrected Person maxima. Residual calibration bias may remain; calibration uncertainty is excluded.")
