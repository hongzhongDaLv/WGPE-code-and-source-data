ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D <- file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
outroot <- file.path(D,'ORIGINAL_CODE_NEW_RESULTS')
Sys.setenv(TOP100_CLOSURE_ROOT=outroot)
new_grid <- utils::read.csv(file.path(D,'GRID_RESULTS.csv'),check.names=FALSE)
new_grid <- new_grid[new_grid$model=='Without_VPD',]
new_profiles <- utils::read.csv(file.path(D,'original_Fig3_layout_no_VPD/LATITUDE_PROFILES_NO_VPD.csv'),check.names=FALSE)
spec <- data.frame(file=c('Figure3a_W-GPE_GPP.csv','Figure3b_z_GPP.csv','Figure3c_W-GPE_LAI.csv','Figure3d_z_LAI.csv'),target=c('GPP','GPP','LAI','LAI'),predictor=c('W-GPE','<z>','W-GPE','<z>'))
# Input adapter only. The archived R source is never written.
read.csv <- function(file,...) {
 if(basename(file)=='Fig3f_source.csv') {
  z<-utils::read.csv(file.path(D,'TEMPORAL_F_NO_VPD/F_CLASS_SUMMARY.csv'))
  z<-z[z$model=='Without_VPD',]
  return(data.frame(target=z$target,metric=paste0('lag3:',z$class),estimate=z$grid_pct))
 }
 k<-match(basename(file),spec$file)
 if(!is.na(k)) {
  z<-new_grid[new_grid$target==spec$target[k]&new_grid$predictor==spec$predictor[k],]
  z$estimate<-z$r
  return(z)
 }
 utils::read.csv(file,...)
}
readRDS <- function(file,...) {
 x<-base::readRDS(file,...)
 if(basename(file)=='Figure3_TOP100_adapter.rds') x$U_EPROF<-new_profiles
 x
}
src<-file.path(ROOT,'CANONICAL_FREEZE/USER_APPROVED_MAIN_FIGURES_20260717_FINAL/plot_scripts/build_Figure3_LAYOUT_v5.R')
for(ex in parse(src)) {
 lhs<-if(is.call(ex)&&identical(ex[[1]],as.name('<-')))as.character(ex[[2]])[1] else ''
 if(lhs=='gsum')break
 # Export only the final requested preview, avoiding old intermediate candidates.
 if(is.call(ex)&&as.character(ex[[1]])[1] %in% c('save_png','save_pdf'))next
 eval(ex)
 if(lhs=='colorbar') {
  cbsource<-file.path(ROOT,'CANONICAL_FREEZE/USER_APPROVED_MAIN_FIGURES_20260717_FINAL/plot_scripts/build_Figure2_LAYOUT_v5_TRUE_INCREMENT.R')
  cbenv<-new.env(parent=environment())
  for(ce in parse(cbsource))if(is.call(ce)&&identical(ce[[1]],as.name('<-'))&&identical(ce[[2]],as.name('colorbar')))eval(ce,envir=cbenv)
  colorbar<-function(lim,pal,title_expr,ticks,accuracy=.1){
   fn<-cbenv$colorbar;en<-new.env(parent=environment(fn));en$pal_blue<-pal;environment(fn)<-en
   fn(lim,title_expr,ticks,accuracy,extend=FALSE)
  }
 }
 if(lhs=='fig4e')fig4e<-fig4e+theme(axis.text.x=element_text(angle=45,hjust=1,size=PANEL_TEXT_SIZE))
}
# Execute the original f-panel and full composition expressions with the new f input.
active<-FALSE
for(ex in parse(src)) {
 lhs<-if(is.call(ex)&&identical(ex[[1]],as.name('<-')))as.character(ex[[2]])[1] else ''
 if(lhs=='temporal_class_colours_eco')eval(ex)
 if(lhs=='fsource')active<-TRUE
 if(!active)next
 if(is.call(ex)&&as.character(ex[[1]])[1] %in% c('save_png','save_pdf'))next
 eval(ex)
 if(lhs=='fig4_final')break
}
save_png(fig4_final,'Figure3_ORIGINAL_CODE_NO_VPD_FULL.png',215,205,600)
save_pdf(fig4_final,'Figure3_ORIGINAL_CODE_NO_VPD_FULL.pdf',215,205)
write.csv(hdf,file.path(OUT,'panel_e_source.csv'),row.names=FALSE)
stopifnot(abs(get_eco_sum('GPP','wgpe')$mean_r_area-.146420775501812)<1e-10)
cat('DONE original source + new a-f data. Stars retain original >=5% significant-area rule.\n')
