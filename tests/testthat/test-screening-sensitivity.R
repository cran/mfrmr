screening_sensitivity_fixture <- function() {
  roster <- expand.grid(Condition='Null',Replicate=1:4,Target=c('A','B'),stringsAsFactors=FALSE)
  roster$Affected <- FALSE
  measures <- roster[c('Condition','Replicate','Target')]
  measures$Infit <- c(.4,.8,1.7,NA, .9,1.1,.9,NA)
  measures$Outfit <- c(.9,.9,1.1,.9, .9,1.1,.9,NA)
  measures$InfitZSTD <- c(-3,-3,4,NA,0,0,0,NA)
  measures$OutfitZSTD <- c(0,0,0,0,0,0,0,NA)
  list(roster=roster,measures=measures,
    bands=data.frame(Profile=c('Broad','Narrow'),Lower=c(.5,.85),Upper=c(1.5,1.2)))
}

test_that('directional screens keep planned families and missing outcomes', {
  f <- screening_sensitivity_fixture()
  x <- mfrm_screening_sensitivity(f$roster,f$measures,f$bands,rule='Prespecified demonstration')
  b <- subset(x$by_family,Profile=='Broad' & Metric=='Any unaffected target flagged')
  expect_equal(b$Positive,c(1,1,2))
  expect_equal(b$Available,c(3,3,3));expect_equal(b$Planned,c(4,4,4))
  expect_equal(b$AllTrialsLower,c(.25,.25,.5))
  expect_equal(b$AllTrialsUpper,c(.5,.5,.75))
  expect_equal(b$MCLower[3],unname(binom.test(2,3)$conf.int[1]))
  expect_equal(subset(x$by_family,Profile=='Narrow' & Direction=='overfit' & Targets>0)$Positive,2)
  expect_identical(x$measures,f$measures)
  expect_identical(summary(x)$by_family,x$by_family)
  dropped <- mfrm_screening_sensitivity(f$roster,f$measures[-8,],f$bands,'same rule')
  expect_identical(dropped$by_family,x$by_family)
  f$measures$Infit[4] <- 1.8
  partial <- mfrm_screening_sensitivity(f$roster,f$measures,f$bands,'same rule')
  expect_equal(subset(partial$by_family,Profile=='Broad' & Direction=='either' & Targets>0)$Available,4)
})

test_that('ZSTD is an explicit separate rule and individual indices are supported', {
  f <- screening_sensitivity_fixture()
  x <- mfrm_screening_sensitivity(f$roster,f$measures,f$bands,'combined',zstd_cut=2)
  expect_equal(subset(x$by_family,Profile=='Broad' & Direction=='overfit' & Targets>0)$Positive,2)
  y <- mfrm_screening_sensitivity(f$roster,f$measures[c('Condition','Replicate','Target','Outfit')],f$bands,'outfit',statistic='outfit')
  expect_equal(subset(y$by_family,Profile=='Broad' & Direction=='either' & Targets>0)$Available,3)
  expect_equal(subset(y$by_family,Profile=='Broad' & Direction=='either' & Targets>0)$Positive,0)
  # Mean-square equality is inside; ZSTD equality is an explicit inclusive flag.
  t <- data.frame(Infit=c(.5,1.5,NA),Outfit=c(1,1,NA),InfitZSTD=c(-2,2,NA),OutfitZSTD=c(0,0,NA))
  expect_identical(mfrm_fit_flags(t,.5,1.5)$either,c(FALSE,FALSE,NA))
  expect_identical(mfrm_fit_flags(t,.5,1.5,zstd_cut=2)$either,c(TRUE,TRUE,NA))
  bad<-f$measures;bad$Infit[1]<-Inf
  expect_true(is.na(subset(mfrm_screening_sensitivity(f$roster,bad,f$bands,'nonfinite')$outcomes,
    Profile=='Broad'&Direction=='either'&Target=='A'&Replicate=='1')$Flag))
  expect_error(mfrm_screening_sensitivity(f$roster,f$measures,f$bands,'bad',zstd_cut=Inf),'positive finite')
  bad<-f$bands;bad$Lower[1]<-2
  expect_error(mfrm_screening_sensitivity(f$roster,f$measures,bad,'bad'),'Lower < Upper')
  bad<-f$measures;bad$Infit[1]<--1
  expect_error(mfrm_screening_sensitivity(f$roster,bad,f$bands,'bad'),'nonnegative')
  expect_error(mfrm_screening_sensitivity(f$roster,rbind(f$measures,f$measures[1,]),f$bands,'bad'),'Duplicate')
  empty<-mfrm_screening_sensitivity(f$roster,f$measures[FALSE,],f$bands,'all failed')
  expect_true(all(is.na(empty$by_family$Rate)))
})

test_that('default fit table separates low MnSq from ZSTD-only evidence', {
  m <- data.frame(Facet='Rater',Level=c('zonly','low','partial','high','mixed'),
    Infit=c(.9,.4,1,1.6,.4),Outfit=c(.9,.9,NA,1.1,1.6),
    InfitZSTD=c(-3,-4,0,3,-4),OutfitZSTD=c(-3,0,NA,0,3))
  current <- fit_measures_table(list(measures=m),sort_by='level',threshold_profiles='all')
  t <- current$table[match(m$Level,current$table$Level),]
  expect_identical(t$FitStatus,c('within_band','overfit','not_available','underfit','mixed'))
  expect_identical(t$ZSTDOnly,c(TRUE,FALSE,NA,FALSE,FALSE))
  expect_false(t$ScreenComplete[3]);expect_true(is.na(t$Underfit[3]))
  expect_false(any(grepl('ZSTD low',current$overfit$ReviewReason)))
  old <- fit_measures_table(list(measures=m),flag_basis='mnsq_or_zstd',sort_by='level')
  expect_identical(subset(old$table,Level=='zonly')$FitStatus,'overfit')
  expect_true(any(grepl('ZSTD low',old$overfit$ReviewReason)))
  expect_true(all(current$profile_summary$FlagBasis=='mnsq'))
  expect_error(fit_measures_table(list(measures=m),flag_basis='unknown'),'arg')
  expect_identical(fit_measure_status_label(-.1,.5,1.5,2,'mnsq'),'not_available')
})

test_that('sensitivity views preserve values, missing cells and display controls', {
  f <- screening_sensitivity_fixture();x <- mfrm_screening_sensitivity(f$roster,f$measures,f$bands,'example')
  grDevices::pdf(tempfile(fileext='.pdf'),width=9,height=6);on.exit(grDevices::dev.off())
  before<-par(c('mar','oma','mfrow','cex'))
  for(style in c('tiles','curves')) {
    expect_silent(plot(x,style=style,direction=c('underfit','overfit'),palette='mono'))
    expect_equal(par(c('mar','oma','mfrow','cex')),before)
    p<-plot(x,style=style,direction='overfit',show_title=FALSE,show_notes=FALSE,draw=FALSE)
    expect_equal(p$data$table$Value,c(1/3,2/3))
    expect_match(p$data$alt_text,'known/planned')
    if(requireNamespace('ggplot2',quietly=TRUE)) {
      gg<-as_ggplot(p);expect_null(gg$labels$title);expect_null(gg$labels$caption)
      expect_silent(ggplot2::ggplotGrob(gg))
    }
  }
  empty<-mfrm_screening_sensitivity(f$roster,f$measures[FALSE,],f$bands,'all unavailable')
  expect_silent(plot(empty))
  if(requireNamespace('ggplot2',quietly=TRUE)) expect_silent(ggplot2::ggplotGrob(as_ggplot(plot(empty,draw=FALSE))))
  p<-plot(empty,quantity='unavailable',draw=FALSE);expect_equal(p$data$table$Value,c(1,1))
  tiny<-x;tiny$by_family$Rate[tiny$by_family$Rate>0]<-.0001
  expect_true(any(grepl('0.01%',plot(tiny,draw=FALSE)$data$table$Label,fixed=TRUE)))
  for(palette in c('accessible','mono')) expect_true(all(
    mfrm_screening_tile_colours(c(seq(0,1,length.out=101),NA),palette)$contrast>=4.5))
  file<-tempfile();saveRDS(x,file);expect_identical(readRDS(file),x)
})
