# Candidate only; reuse the aligned-map pipeline and exact panel allocations.
ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D <- file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
style_hook <- function(env,fig,lhs) {
 if(fig==1 && lhs=='panel12d' && inherits(env$panel12d,'ggplot') && length(env$panel12d$layers)==3) {
  p<-env$panel12d
  p$layers[[2]]<-geom_point(shape=21,size=2.1,stroke=.25,colour='#FFFFFF')
  p$layers[[3]]<-geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.12,linewidth=.4,colour='#000000')
  p$layers<-p$layers[c(1,3,2)]
  env$panel12d<-p
 }
 if(fig==1 && lhs=='panel12f' && inherits(env$panel12f,'ggplot') && length(env$panel12f$layers)==3) {
  p<-env$panel12f;p$layers[[2]]$geom_params$width<-.40
  p$layers[[2]]$aes_params$linewidth<-.08
  p$layers[[3]]$aes_params$linewidth<-.4;p$layers[[3]]$aes_params$colour<-'#000000'
  p$layers[[3]]$geom_params$width<-.15
  env$panel12f<-p
 }
 if(fig==1 && lhs=='f_bottom') {
  p<-ggdraw()
  labels<-c('Mass\nredistribution','Within-layer\nmean-height','Higher-order\ninteraction','Actual')
  for(k in 1:4)p<-p+draw_label(labels[k],x=c(.20,.42,.64,.85)[k],y=.90,
      angle=45,hjust=1,vjust=1,lineheight=.95,fontfamily=env$FONT,size=8.2)
  env$f_bottom<-p
 }
 if(fig==2 && lhs=='c' && inherits(env$c,'ggplot')) {
  p<-env$c;p$data$h<-p$data$h*.8
  p$layers[[3]]$aes_params$linewidth<-.4
  p$layers[[3]]$aes_params$colour<-'#000000'
  env$c<-p
 }
 if(fig==2 && lhs=='e' && inherits(env$e,'ggplot')) {
  p<-env$e
  p$layers[[2]]$position<-position_dodge(width=.34)
  p$layers[[2]]$geom_params$width<-.10
  p$layers[[2]]$aes_params$colour<-'#000000'
  p$layers[[2]]$aes_params$linewidth<-.4
  p$layers[[2]]$show.legend<-FALSE
  p$layers[[3]]$position<-position_dodge(width=.34)
  p$layers[[3]]$aes_params$size<-1.7
  p<-p+guides(colour=guide_legend(nrow=1,override.aes=list(shape=16,size=1.7)))
  env$e<-p
 }
}
inject <- function(x) {
 if(is.call(x)&&identical(x[[1]],as.name('if'))&&grepl('file.exists',paste(deparse(x[[2]]),collapse=' ')))return(quote(invisible(NULL)))
 if(identical(x,quote(eval(ex,envir=env))))return(quote({eval(ex,envir=env);style_hook(env,fig,lhs)}))
 if(is.call(x))for(i in seq_along(x))if(i>1)x[i]<-list(inject(x[[i]]))
 x
}
for(z in parse(file.path(D,'fix_fig12_profiles_v1.R'))) {
 lhs<-if(is.call(z)&&identical(z[[1]],as.name('<-')))as.character(z[[2]])[1] else ''
 if(lhs=='DEST') {DEST<-file.path(D,'FULL_FIGURE_STYLE_PREVIEW_v2');next}
 eval(inject(z))
}
