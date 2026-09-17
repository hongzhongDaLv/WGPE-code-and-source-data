# New candidate, preserving all earlier scripts and outputs.
ROOT<-Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D<-file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
.libPaths(.libPaths())
library(ggplot2);library(cowplot);library(grid)
for(z in parse(file.path(D,'preview_full_figures_refined_v2.R'))) {
 if(is.call(z)&&identical(z[[1]],as.name('for')))break
 eval(z)
}
previous_hook<-style_hook
# Fixed physical thickness in the drawing viewport, independent of data scales.
fixed_rect <- function(parent,vertical=FALSE) ggproto(NULL,parent,
 draw_panel=function(self,data,panel_params,coord,...) {
  g<-parent$draw_panel(data,panel_params,coord,...)
  adjust<-function(a) {
   if(inherits(a,'rect')) {
    if(vertical){a$x<-a$x+a$width/2;a$width<-unit(2.6,'mm');a$hjust<-.5}
    else {a$y<-a$y-a$height/2;a$height<-unit(2.6,'mm');a$vjust<-.5}
   }
   if(!is.null(a$children))for(j in seq_along(a$children))a$children[[j]]<-adjust(a$children[[j]])
   a
  }
  adjust(g)
 })
style_hook<-function(env,fig,lhs) {
 previous_hook(env,fig,lhs)
 key<-if(fig==1&&lhs=='left_e')'left_e' else if(fig==1&&lhs=='panel12g')'panel12g' else
      if(fig==2&&lhs=='c')'c' else if(fig==2&&lhs=='d')'d' else
      if(fig==1&&lhs=='panel12f')'panel12f' else ''
 if(nzchar(key)) {
  p<-env[[key]]
  if(inherits(p,'ggplot'))for(k in seq_along(p$layers)) {
   l<-p$layers[[k]]
   if(inherits(l$geom,'GeomCol')||inherits(l$geom,'GeomRect')) {
    l$geom<-fixed_rect(l$geom,vertical=(fig==1&&key=='panel12f'))
    l$aes_params$colour<-NA
   }
   if(inherits(l$geom,'GeomErrorbar')) {l$aes_params$colour<-'#000000';l$aes_params$linewidth<-.4}
   p$layers[[k]]<-l
  }
  if(inherits(p,'ggplot'))p<-p+theme(panel.grid.major=element_blank())
  env[[key]]<-p
 }
 # Let the axis layout engine align labels to bars instead of manual positions.
 if(fig==1&&lhs=='panel12f_core')env$panel12f_core<-env$panel12f+
    labs(title=NULL,x=NULL)+theme(axis.text.x=element_text(size=env$PANEL_TEXT_SIZE,
    angle=45,hjust=1,vjust=1,lineheight=.92),axis.ticks.x=element_blank(),plot.margin=margin(1,4,3,4))
 if(fig==1&&lhs=='f_bottom')env$f_bottom<-ggdraw()
 if(fig==1&&lhs=='panel12f'&&inherits(env$panel12f,'ggplot')&&length(env$panel12f$layers)==3&&
    exists('f_bottom',envir=env,inherits=FALSE)) {
   env$panel12f<-ggdraw()+draw_plot(env$f_title,x=0,y=.910,width=1,height=.090)+
      draw_plot(env$panel12f_core,x=0,y=0,width=1,height=.910)
 }
}
old_inject<-inject
inject3<-function(x) {
 if(is.call(x)&&identical(x[[1]],as.name('cairo_pdf')))
   return(quote(png(file.path(env$OUT,'measurement.png'),width=W,height=H,units='mm',res=600)))
 # Same panel width in millimetres; compensate outer margins separately.
 if(identical(x,quote(sg<-ggplotGrob(sp))))return(quote({
   sg<-ggplotGrob(sp)
   profile_panel_index<-sg$layout[sg$layout$name=='panel',]
   sg$widths[profile_panel_index$l]<-unit(11,'mm')
   profile_width_mm<-sum(convertWidth(sg$widths,'mm',valueOnly=TRUE))
 }))
 if(is.call(x)&&identical(x[[1]],as.name('draw_grob'))) {
   x$width<-quote(profile_width_mm/W)
 }
 if(is.call(x))for(i in seq_along(x))if(i>1)x[i]<-list(inject3(x[[i]]))
 x
}
for(z in parse(file.path(D,'fix_fig12_profiles_v1.R'))) {
 lhs<-if(is.call(z)&&identical(z[[1]],as.name('<-')))as.character(z[[2]])[1] else ''
 if(lhs=='DEST'){DEST<-file.path(D,'FULL_FIGURE_STYLE_PREVIEW_v3');next}
 eval(inject3(old_inject(z)))
}
