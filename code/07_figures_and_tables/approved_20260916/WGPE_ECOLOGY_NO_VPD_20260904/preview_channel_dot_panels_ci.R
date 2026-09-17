.libPaths(.libPaths())
base <- file.path(Sys.getenv('WGPE_PROJECT_ROOT', unset='.'),'WGPE_ECOLOGY_NO_VPD_20260904')
exprs <- parse(file.path(base, 'preview_channel_dot_panels.R'))
attach_ci <- function(dat, u, metrics, divisor=1) {
  key <- paste(dat$region, metrics[as.character(dat$channel)])
  uk <- paste(u$region,u$metric)
  stopifnot(!anyDuplicated(uk))
  j <- match(key,uk)
  stopifnot(!anyNA(j),max(abs(dat$value-u$estimate[j]/divisor)) < 1e-7)
  dat$lower95 <- u$lower95[j]/divisor
  dat$upper95 <- u$upper95[j]/divisor
  stopifnot(all(is.finite(dat$lower95)),all(dat$lower95<=dat$upper95))
  dat
}
for (ex in exprs) {
  lhs <- if(is.call(ex)&&identical(ex[[1]],as.name('<-'))) as.character(ex[[2]])[1] else ''
  if(lhs=='OUT') {
    OUT <- file.path(base,'CHANNEL_DOT_PREVIEW_CI')
    next
  }
  if(lhs=='p1') {
    u1 <- readRDS(file.path(ROOT,'output_TOP100_revision/derived_data/figure_source_data/Figure1_TOP100_adapter.rds'))$U_CORE
    u1 <- u1[u1$height_reference=='absolute',]
    a <- attach_ci(a,u1,c('IWV'='decomp_IWV','<z>'='decomp_z','Interaction'='decomp_interaction'))
    u2 <- readRDS(file.path(ROOT,'output_TOP100_revision/derived_data/figure_source_data/Figure2_TOP100_adapter.rds'))$enso_uncertainty
    u2 <- u2[u2$height_reference=='absolute' & u2$event_type=='ENSO_difference',]
    b <- attach_ci(b,u2,c('IWV'='IWV_channel','<z>'='z_channel','Interaction'='interaction_channel'),1000)
    old_plot <- make_plot
    make_plot <- function(dat,title,xlab) {
      p <- old_plot(dat,title,xlab)
      p$layers <- append(p$layers,list(geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.09,colour='black',linewidth=.4,show.legend=FALSE)),after=2)
      p
    }
  }
  eval(ex,envir=.GlobalEnv)
}
