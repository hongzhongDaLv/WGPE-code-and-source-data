# Generate width variants from the same composed vector grobs, without rerunning
# scientific calculations for each option or changing any other layout element.
ROOT<-Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D<-file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
for(z in parse(file.path(D,'preview_full_figures_refined_v3.R'))) {
 if(is.call(z)&&identical(z[[1]],as.name('for'))&&grepl('fix_fig12',paste(deparse(z),collapse=' ')))break
 eval(z)
}
resize_bars<-function(g,mm) {
 if(inherits(g,'rect'))for(n in c('width','height')) {
  u<-g[[n]]
  if(length(u)==1L&&identical(unitType(u),'mm')&&abs(as.numeric(u)-2.6)<1e-8)g[[n]]<-unit(mm,'mm')
 }
 if(!is.null(g$children))for(i in seq_along(g$children))g$children[[i]]<-resize_bars(g$children[[i]],mm)
 if(!is.null(g$grobs))for(i in seq_along(g$grobs))g$grobs[[i]]<-resize_bars(g$grobs[[i]],mm)
 g
}
export_options<-function(base,fig,W,H) {
 g<-grid.force(as_grob(base))
 for(mm in c(2.6,3.2,3.8,4.4)) {
  out<-file.path(DEST,sprintf('width_%.1fmm',mm));dir.create(out,recursive=TRUE,showWarnings=FALSE)
  gg<-resize_bars(g,mm)
  ggsave(file.path(out,paste0('Figure',fig,'.png')),gg,width=W,height=H,units='mm',dpi=600,bg='white')
 }
}
inject_width<-function(x) {
 if(is.call(x)&&identical(x[[1]],as.name('ggsave'))) {
  if(grepl('latitude_aligned.png',paste(deparse(x),collapse=' ')))return(quote(export_options(base,fig,W,H)))
  return(quote(invisible(NULL)))
 }
 if(is.call(x))for(i in seq_along(x))if(i>1)x[i]<-list(inject_width(x[[i]]))
 x
}
for(z in parse(file.path(D,'fix_fig12_profiles_v1.R'))) {
 lhs<-if(is.call(z)&&identical(z[[1]],as.name('<-')))as.character(z[[2]])[1] else ''
 if(lhs=='DEST'){DEST<-file.path(D,'BAR_WIDTH_OPTIONS');dir.create(DEST,recursive=TRUE,showWarnings=FALSE);next}
 eval(inject_width(inject3(old_inject(z))))
}
# Contact sheet: run make_bar_width_contact_sheet.py after these exports.
