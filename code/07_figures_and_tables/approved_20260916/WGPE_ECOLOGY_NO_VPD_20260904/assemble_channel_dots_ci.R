# New full-figure candidate: replace only Fig.1g and Fig.2d.
ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D <- file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
for(ex0 in parse(file.path(D,'preview_bar_width_options.R'))) {
 if(is.call(ex0)&&identical(ex0[[1]],as.name('for'))&&grepl('fix_fig12_profiles',paste(deparse(ex0),collapse=' ')))break
 eval(ex0)
}
dot_env <- new.env(parent=globalenv())
for(ex0 in parse(file.path(D,'preview_channel_dot_panels.R'))) {
 lhs0<-if(is.call(ex0)&&identical(ex0[[1]],as.name('<-')))as.character(ex0[[2]])[1] else ''
 if(lhs0=='p1')break
 eval(ex0,envir=dot_env)
}
dot_a <- read.csv(file.path(D,'CHANNEL_DOT_PREVIEW_CI/Fig1g_values.csv'))
dot_b <- read.csv(file.path(D,'CHANNEL_DOT_PREVIEW_CI/Fig2d_values.csv'))
hook_before_dots <- style_hook
style_hook <- function(env,fig,lhs) {
 hook_before_dots(env,fig,lhs)
 if((fig==1&&lhs=='panel12g')||(fig==2&&lhs=='d')) {
   dat<-if(fig==1)dot_a else dot_b
   dat$channel<-factor(dat$channel,levels=c('IWV','<z>','Interaction'))
   title<-if(fig==1)expression(bold(g)~'Regional W-GPE composition') else expression(bold(d)~'ENSO decomposition')
   label<-if(fig==1)expression('Contribution (J m'^{-2}~'yr'^{-1}*')') else expression('Channel anomaly (10'^3~'J m'^{-2}*')')
   p<-dot_env$make_plot(dat,title,label)
   p$layers<-append(p$layers,list(geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.09,colour='black',linewidth=.4,show.legend=FALSE)),after=2)
   p<-p+theme(text=element_text(family=env$FONT),
      axis.text=element_text(size=env$PANEL_TEXT_SIZE),axis.title=element_text(size=env$PANEL_TEXT_SIZE),
      legend.text=element_text(size=env$PANEL_TEXT_SIZE),plot.title=element_text(size=env$PANEL_TITLE_SIZE),plot.title.position='plot',
      plot.margin=margin(2,4,3,4))
   env[[lhs]]<-p
 }
}
export_options <- function(base,fig,W,H) {
 g<-resize_bars(grid.force(as_grob(base)),3.2)
 ggsave(file.path(DEST,paste0('Figure',fig,'_channel_dots_CI.png')),g,width=W,height=H,units='mm',dpi=600,bg='white')
}
for(ex0 in parse(file.path(D,'fix_fig12_profiles_v1.R'))) {
 lhs0<-if(is.call(ex0)&&identical(ex0[[1]],as.name('<-')))as.character(ex0[[2]])[1] else ''
 if(lhs0=='DEST'){DEST<-file.path(D,'FULL_FIGURES_CHANNEL_DOTS_CI');dir.create(DEST,recursive=TRUE,showWarnings=FALSE);next}
 eval(inject_width(inject3(old_inject(ex0))))
}
