# Preserve plot allocations; render the three titles on one shared baseline.
ROOT<-Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D<-file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
for(ex1 in head(parse(file.path(D,'assemble_channel_dots_ci.R')),-1)) {
 eval(ex1)
}
before_title_hook<-style_hook
style_hook<-function(env,fig,lhs) {
 before_title_hook(env,fig,lhs)
 if(fig==2&&lhs=='bottom_row') {
   panels<-lapply(list(env$c,env$d,env$e),function(p)p+theme(plot.title=element_text(colour='transparent')))
   row<-plot_grid(plotlist=panels,nrow=1,rel_widths=c(1,1,1),align='h',axis='tb')
   row<-ggdraw(row)
   titles<-list(expression(bold(c)~'Wet/dry decomposition'),expression(bold(d)~'ENSO decomposition'),expression(bold(e)~'Information gain beyond IWV'))
   for(k in 1:3) row<-row+draw_label(titles[[k]],x=(k-1)/3+.009,y=.982,hjust=0,vjust=1,fontfamily=env$FONT,size=env$PANEL_TITLE_SIZE)
   env$bottom_row<-row
 }
}
only_fig2<-function(x) {
 if(is.call(x)&&identical(x[[1]],as.name('for'))&&identical(x[[2]],as.name('fig')))x[[3]]<-2L
 if(is.call(x))for(i in seq_along(x))if(i>1)x[i]<-list(only_fig2(x[[i]]))
 x
}
for(ex1 in parse(file.path(D,'fix_fig12_profiles_v1.R'))) {
 lhs1<-if(is.call(ex1)&&identical(ex1[[1]],as.name('<-')))as.character(ex1[[2]])[1] else ''
 if(lhs1=='DEST'){DEST<-file.path(D,'FULL_FIGURES_CHANNEL_DOTS_CI_ALIGNED');dir.create(DEST,recursive=TRUE,showWarnings=FALSE);next}
 eval(only_fig2(inject_width(inject3(old_inject(ex1)))))
}
