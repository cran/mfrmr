#' Wright and Pathway maps for extended RSMs
#'
#' Plot a saved [mfrm_results()] object using `type = "wright"` or
#' `type = "fit_pathway"`. These routes support testlet and shared-rater models;
#' ordinary-model maps keep their existing definitions.
#'
#' @section Locations and their meaning:
#' The fitted equation is log(P(k)/P(k-1)) = ability + local effect -
#' sum(facet severities) - step(k). A fixed facet's displayed location is its
#' severity plus the unweighted mean of the fitted steps. It is the mean
#' adjacent-category boundary when all other facet severities and local effects
#' are zero. It is not the ability producing the middle expected score.
#' Shared raters use the fitted conditional rater mode plus that step mean;
#' replacement-rater marginal predictions are a different target. A testlet
#' effect is Person-local and has no separate global facet column.
#'
#' Person positions are saved conditional EAPs on the fitted mean-zero ability
#' scale. The separate category-boundary column shows fitted steps at zero
#' facet/local effects. These locations are conditional references: integrating
#' latent effects does not generally preserve these category crossings.
#' `SourceEstimate` and `StepCenter` retain each transformation in plot tables.
#' Do not add displayed facet locations to obtain a combined rating boundary:
#' this would count the mean step repeatedly. Use the fitted equation instead.
#' Whiskers, when requested, show **only Person conditional intervals**.
#' Composite facet/step uncertainty is not obtained by shifting an existing
#' facet interval; use the separate estimate plots for their original targets.
#'
#' @section Saved inputs:
#' Supply source-roster scores through `mfrm_results(fit, scores = scores)`
#' for either extension (testlet `predictions = scores` remains supported).
#' Calculate them explicitly using [score_mfrm_persons()]. Wright maps require
#' matching scores; older testlet scores without `scoring_data` need rescoring,
#' not refitting. The source calibration and complete rating-event multiset
#' must match. Shared-rater scores computed from a reduced roster cannot be
#' combined with rater modes from the full fitted roster.
#'
#' Pathways also require separately saved [mfrm_response_diagnostics()] output.
#' `fit_stat = "Infit"` (default) or `"Outfit"` chooses the horizontal axis;
#' the vertical axis uses the locations defined above. Selected-row counts and
#' full-roster conditioning are retained. A residual index for a subset of
#' ratings describes that subset, not the group's entire workload. Missing or
#' unresolved indices remain in tables but cannot be drawn. No expectation-one
#' line, acceptance band, ZSTD, bias test or rater-quality classification is used.
#' A pathway can display fixed facets/raters without Person scores.
#'
#' @section Display controls:
#' `facet = NULL` shows all available location panels in a Wright map and
#' all located diagnostic groups in a pathway; supply fitted column names to
#' select panels. `persons = NULL` retains all saved Person scores; a character
#' vector selects saved IDs without rescoring. `show_steps = TRUE` includes
#' the separate Wright category column; steps do not appear in pathways.
#' `show_intervals = TRUE` shows available Person conditional intervals.
#' `show_labels` defaults to `FALSE` for Wright maps and `TRUE` for pathways;
#' set it to `FALSE` for crowded pathways. `palette` is
#' `"accessible"` or `"mono"`; shapes also distinguish the location types
#' and prior-only results. A prior-only point is not a measured ability.
#' `title`/`caption` replace defaults; `show_title`/`show_notes` hide annotations
#' without removing interpretation metadata. `text_scale` and `point_size`
#' change sizes, and `draw = FALSE` returns saved plot data. [as_ggplot()]
#' preserves these choices and supplies alternative text. Plot tables retain
#' all selected unavailable rows; `locations` also retains undisplayed panels.
#' No fitting, scoring or integration is performed during plotting/export.
#'
#' @examples
#' example <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr"))
#' res <- mfrm_results(example$testlet$fit,
#'   scores = example$testlet$scores,
#'   response_diagnostics = example$testlet$diagnostics, compute = "never")
#' plot(res, type = "wright")
#' pathway <- plot(res, type = "fit_pathway", facet = "Rater", draw = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) as_ggplot(pathway)
#' @name mfrmr_model_maps
NULL

mfrm_model_locations <- function(fit, scores = NULL) {
  if (!mfrm_extended_fit(fit) || !isTRUE(fit$checks$NumericalReady) || !isTRUE(fit$checks$InformationPositive)) stop("Resolve numerical checks before drawing model locations.",call.=FALSE)
  steps <- fit$calibration$steps
  if (length(steps)!=length(fit$input$score_levels)-1L || any(!is.finite(steps))) stop("Complete finite step calibration is required for model locations.",call.=FALSE)
  center <- mean(steps)
  make <- function(facet,level,value,kind,shift=center) data.frame(Facet=facet,Level=as.character(level),
    SourceEstimate=value,StepCenter=shift,Estimate=value+shift,Kind=kind,
    Lower=NA_real_,Upper=NA_real_,Status=ifelse(is.finite(value),"available_location","unavailable"),
    Reason=ifelse(is.finite(value),"","Location unavailable"),stringsAsFactors=FALSE)
  fixed <- fit$calibration_table[fit$calibration_table$Parameter=="Fixed facet",,drop=FALSE]
  parts <- list()
  if(nrow(fixed)) parts[[length(parts)+1L]] <- make(fixed$Facet,fixed$Level,fixed$Estimate,"Fixed facet")
  if(inherits(fit,"mfrm_random_rater")) parts[[length(parts)+1L]] <- make(fit$input$columns$rater,
    fit$raters$Rater,fit$raters$Estimate,"Observed rater")
  if(!is.null(scores)) {
    mfrm_validate_person_scores(fit,scores)
    z <- scores$table
    p <- make(fit$input$columns$person,z$Person,z$Estimate,"Person EAP",0)
    p$Lower<-z$Lower; p$Upper<-z$Upper; p$Status<-z$Status; p$Reason<-z$Reason
    p$Kind[p$Status=="prior_only"] <- "Prior only"
    p$Estimate[p$Status=="unavailable"] <- NA_real_
    parts <- c(list(p),parts)
  }
  parts[[length(parts)+1L]] <- make(NA_character_,fit$input$score_levels[-1L],steps,"Category boundary",0)
  table <- do.call(rbind,parts);rownames(table)<-NULL
  step_panel <- tail(make.unique(c(unique(stats::na.omit(table$Facet)),"Category boundaries"),sep=" "),1L)
  table$Panel <- ifelse(is.na(table$Facet),step_panel,table$Facet)
  table
}

mfrm_model_map <- function(x, type, facet=NULL, persons=NULL, fit_stat=c("Infit","Outfit"),
    show_steps=TRUE, show_intervals=TRUE, draw=TRUE, palette=c("accessible","mono"),
    title=NULL, caption=NULL, show_title=FALSE, show_notes=FALSE, show_labels=type=="fit_pathway",
    text_scale=1, point_size=2.5, ...) {
  rlang::check_dots_empty();metric<-tolower(match.arg(fit_stat));palette<-match.arg(palette)
  for(key in c("show_steps","show_intervals","draw","show_title","show_notes","show_labels")) if(!is.logical(get(key)) || length(get(key))!=1L || is.na(get(key))) stop("Display switches must be TRUE or FALSE.",call.=FALSE)
  for(key in c("title","caption")) if(!is.null(get(key)) && (!is.character(get(key)) || length(get(key))!=1L || is.na(get(key)))) stop("Titles/captions must be NULL or one string.",call.=FALSE)
  for(key in c("text_scale","point_size")) if(!is.numeric(get(key)) || length(get(key))!=1L || !is.finite(get(key)) || get(key)<=0) stop("Text and point sizes must be positive finite numbers.",call.=FALSE)
  scores <- x$scores %||% if(inherits(x$fit,"mfrm_testlet")) x$predictions
  if(type=="wright") mfrm_validate_person_scores(x$fit,scores)
  if(type=="fit_pathway" && !is.null(scores) && !isTRUE(tryCatch(mfrm_validate_person_scores(x$fit,scores),error=function(e) FALSE))) scores<-NULL
  locations <- mfrm_model_locations(x$fit,scores)
  person_col <- x$fit$input$columns$person
  if(!is.null(persons) && (!is.character(persons) || !length(persons) || anyNA(persons) || anyDuplicated(persons) || !all(persons %in% locations$Level[locations$Facet %in% person_col]))) stop("Choose distinct saved Person IDs.",call.=FALSE)
  available_facets <- unique(stats::na.omit(locations$Facet))
  if(type=="fit_pathway") {
    if(is.null(x$diagnostics)) stop("Supply saved posterior response diagnostics for pathways.",call.=FALSE)
    mfrm_validate_response_diagnostics(x$fit,x$diagnostics)
    available_facets <- intersect(available_facets,x$diagnostics$settings$group_by)
  }
  facet <- facet %||% available_facets
  if(!is.character(facet) || !length(facet) || anyNA(facet) || anyDuplicated(facet) || !all(facet %in% available_facets)) stop("Choose available location/diagnostic facet columns.",call.=FALSE)
  tab <- locations[locations$Facet %in% facet | (type=="wright" & show_steps & is.na(locations$Facet)),,drop=FALSE]
  if(!is.null(persons)) tab <- tab[!(tab$Facet %in% person_col) | tab$Level %in% persons,,drop=FALSE]
  panels <- c(facet,if(type=="wright" && show_steps) unique(locations$Panel[is.na(locations$Facet)]))
  tab$Panel <- factor(tab$Panel,levels=panels)
  tab$X <- as.numeric(tab$Panel)
  if(type=="wright") {
    for(panel in panels) {
      ids <- which(tab$Panel==panel)
      if(length(ids)>1L) tab$X[ids] <- tab$X[ids] + seq(-.18,.18,length.out=length(ids))
    }
  } else {
    tab$X <- NA_real_; tab$Selected <- tab$Observed <- tab$Missing <- NA_integer_
    tab$DiagnosticStatus <- tab$DiagnosticReason <- NA_character_
    for(f in facet) {
      at <- which(tab$Facet==f); d <- x$diagnostics$measures[x$diagnostics$measures$Facet==f,,drop=FALSE]
      pos <- match(tab$Level[at],d$Level)
      tab$X[at] <- d[[if(metric=="infit") "Infit" else "Outfit"]][pos]
      for(field in c("Selected","Observed","Missing")) tab[[field]][at] <- d[[field]][pos]
      tab$DiagnosticStatus[at] <- d$Status[pos]; tab$DiagnosticReason[at] <- d$Reason[pos]
      tab$X[at][is.na(pos) | d$Status[pos]!="descriptive_only"] <- NA_real_
    }
  }
  tab$Included <- is.finite(tab$X) & is.finite(tab$Estimate)
  tab$DisplayReason <- ifelse(tab$Included,"",ifelse(!is.finite(tab$Estimate),"Location unavailable","No available selected-row residual index"))
  tab$IntervalShown <- show_intervals & tab$Included & tab$Kind %in% c("Person EAP","Prior only") & is.finite(tab$Lower) & is.finite(tab$Upper)
  tab$Shape <- c("Person EAP"=16,"Prior only"=1,"Fixed facet"=17,"Observed rater"=15,"Category boundary"=18)[tab$Kind]
  tab$Colour <- if(palette=="mono") "#222222" else c("Person EAP"="#0072B2","Prior only"="#0072B2","Fixed facet"="#9C4700","Observed rater"="#7B3294","Category boundary"="#444444")[tab$Kind]
  note <- paste("Facet/rater positions include mean step; other effects zero. Whiskers: Person conditional intervals only.",
    "Open circles: prior only.",
    if(type=="fit_pathway") "Residuals describe selected ratings; same-data posterior, no fit cutoffs.",
    sprintf("%d of %d selected locations displayed; prior-only points are not measured abilities.",sum(tab$Included),nrow(tab)))
  selection <- data.frame(SourcePersons=length(unique((x$fit$input$assigned_data %||% x$fit$input$data)[[person_col]])),
    SavedPersons=if(is.null(scores)) 0L else nrow(scores$table),SelectedPersons=sum(tab$Facet %in% person_col),
    DisplayedPersons=sum(tab$Facet %in% person_col & tab$Included))
  if(person_col %in% facet) note <- paste(note,sprintf("Persons: %d of %d source IDs selected.",selection$SelectedPersons,selection$SourcePersons))
  payload <- list(table=tab,locations=locations,person_selection=selection,models=x$model_family,
    title=title %||% paste(x$model_family,if(type=="wright") "Wright map" else "descriptive pathway"),
    caption=caption %||% note,notes=data.frame(Note=note),
    alt_text=paste(if(type=="wright") "Wright reference-location map." else "Descriptive residual versus reference-location pathway.",note),
    display=list(type=type,metric=metric,panels=panels,palette=palette,show_title=show_title,
      show_notes=show_notes,show_labels=show_labels,text_scale=text_scale,point_size=point_size))
  out<-new_mfrm_plot_data("extended_model_map",payload)
  if(draw) mfrm_draw_model_map(payload)
  invisible(out)
}

mfrm_draw_model_map <- function(payload) {
  opt<-payload$display;tab<-payload$table
  range_finite <- function(z,pad=.3) {z<-z[is.finite(z)];if(!length(z)) z<-c(-1,1);r<-range(z);r+c(-1,1)*max(pad,diff(r)*.15)}
  y<-range_finite(c(tab$Estimate[tab$Included],tab$Lower[tab$IntervalShown],tab$Upper[tab$IntervalShown]))
  panels<-if(opt$type=="wright") "all" else opt$panels
  nc<-min(3,ceiling(sqrt(length(panels))))
  old<-graphics::par(mfrow=c(ceiling(length(panels)/nc),nc),mar=c(5,5,3,1),oma=c(if(opt$show_notes) 5 else 0,0,if(opt$show_title) 3 else 0,0),cex=opt$text_scale)
  on.exit(graphics::par(old),add=TRUE)
  for(panel in panels) {
    d<-if(opt$type=="wright") tab else tab[tab$Panel==panel,,drop=FALSE]
    xr<-if(opt$type=="wright") c(.5,length(opt$panels)+.5) else range_finite(d$X[d$Included],.05)
    graphics::plot(NA_real_,xlim=xr,ylim=y,xaxt=if(opt$type=="wright") "n" else "s",
      xlab=if(opt$type=="wright") "" else paste("Descriptive",if(opt$metric=="infit") "Infit" else "Outfit"),
      ylab="Conditional reference location (logits)",main=if(opt$type=="wright") "" else panel,bty="n")
    if(opt$type=="wright") graphics::axis(1,at=seq_along(opt$panels),labels=opt$panels)
    ci<-d[d$IntervalShown,,drop=FALSE];used<-d[d$Included,,drop=FALSE]
    graphics::segments(ci$X,ci$Lower,ci$X,ci$Upper,col=ci$Colour)
    graphics::points(used$X,used$Estimate,pch=used$Shape,col=used$Colour,cex=opt$point_size/2.5)
    if(opt$show_labels) graphics::text(used$X,used$Estimate,used$Level,pos=3,cex=.7)
    if(!nrow(used)) graphics::text(mean(xr),mean(y),"No available locations")
  }
  if(opt$show_title) graphics::mtext(payload$title,side=3,outer=TRUE,line=1)
  if(opt$show_notes) graphics::mtext(paste(strwrap(payload$caption,105),collapse="\n"),side=1,outer=TRUE,line=1,cex=.7)
  invisible(payload)
}

.mfrmr_gg_model_map <- function(payload) {
  opt<-payload$display;tab<-payload$table;used<-tab[tab$Included,,drop=FALSE];ci<-tab[tab$IntervalShown,,drop=FALSE]
  kinds<-unique(tab$Kind);colours<-setNames(tab$Colour[match(kinds,tab$Kind)],kinds);shapes<-setNames(tab$Shape[match(kinds,tab$Kind)],kinds)
  p<-ggplot2::ggplot(tab,ggplot2::aes(.data$X,.data$Estimate))+
    ggplot2::geom_linerange(data=ci,ggplot2::aes(ymin=.data$Lower,ymax=.data$Upper,colour=.data$Kind),show.legend=FALSE)+
    ggplot2::geom_point(data=used,ggplot2::aes(colour=.data$Kind,shape=.data$Kind),size=opt$point_size)+
    ggplot2::scale_colour_manual(name=NULL,values=colours,drop=FALSE)+
    ggplot2::scale_shape_manual(name=NULL,values=shapes,drop=FALSE)+
    ggplot2::theme_minimal(base_size=11*opt$text_scale)+ggplot2::theme(legend.position="bottom",panel.spacing=ggplot2::unit(1.5,"lines"))
  if(opt$show_labels) p<-p+ggplot2::geom_text(data=used,ggplot2::aes(label=.data$Level),vjust=-1,size=3*opt$text_scale)
  if(opt$type=="wright") p<-p+ggplot2::scale_x_continuous(breaks=seq_along(opt$panels),labels=opt$panels,limits=c(.5,length(opt$panels)+.5)) else
    p<-p+ggplot2::facet_wrap(~Panel,ncol=min(3,ceiling(sqrt(length(opt$panels)))),drop=FALSE)
  if(opt$type=="wright" && !nrow(used)) p<-p+ggplot2::annotate("text",x=mean(seq_along(opt$panels)),y=0,label="No available locations")
  if(opt$type=="fit_pathway") {
    empty<-setdiff(opt$panels,as.character(unique(used$Panel)))
    if(length(empty)) p<-p+ggplot2::geom_text(data=data.frame(Panel=factor(empty,levels=opt$panels),
      X=if(nrow(used)) mean(range(used$X)) else 0,Estimate=if(nrow(used)) mean(range(used$Estimate)) else 0),
      label="No available locations")
  }
  p<-.mfrmr_gg_labs(p,payload,x=if(opt$type=="wright") NULL else paste("Descriptive",if(opt$metric=="infit") "Infit" else "Outfit"),y="Conditional reference location (logits)")+ggplot2::labs(alt=payload$alt_text)
  attr(p,"mfrmr_alt_text")<-payload$alt_text
  p
}
