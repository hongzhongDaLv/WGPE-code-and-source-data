ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
DEST <- file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904/FIG12_PROFILE_ALIGNMENT')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
.libPaths(.libPaths())
library(ggplot2);library(cowplot);library(grid)
pos <- new.env(); profiles <- list(); idx <- 0L
drawDetails.record12 <- function(x,recording=TRUE) {
 a<-deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE)
 b<-deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
 pos[[x$key]]<-c(a$x,a$y,b$x,b$y)*25.4
}
build_top <- function(mp,side,sp) {
 idx <<- idx+1L; key<-as.character(idx);profiles[[key]]<<-sp
 mg<-ggplotGrob(mp);p<-mg$layout[mg$layout$name=='panel',]
 mg<-gtable::gtable_add_grob(mg,grob(key=paste0(key,'m'),cl='record12'),t=p$t,l=p$l,b=p$b,r=p$r)
 sg<-ggplotGrob(side);sg$grobs<-lapply(sg$grobs,function(x)nullGrob())
 sg<-gtable::gtable_add_grob(sg,grob(key=paste0(key,'s'),cl='record12'),t=1,l=1,b=length(sg$heights),r=length(sg$widths))
 plot_grid(mg,sg,nrow=1,rel_widths=c(.80,.20),align='h',axis='tb')
}
for(fig in 1:2) {
 if(file.exists(file.path(DEST,'Figure1/Figure1_latitude_aligned.png')) && fig==1)next
 idx<-0L;profiles<-list();pos<-new.env()
 env<-new.env(parent=globalenv())
 src<-file.path(ROOT,sprintf('WGPE_V23_APPROVED_WORD_BUILD/01_SCRIPTS/build_Figure%d_v23_labels.R',fig))
 expressions<-parse(src,keep.source=TRUE)
 for(ii in seq_along(expressions)) {
  ex<-expressions[[ii]]
  ln<-attr(expressions,'srcref')[[ii]][1]
  if(fig==1 && ln>=622 && ln<724 && !grepl('^wgpe_std <-',paste(deparse(ex),collapse=' ')))next
  hd<-if(is.call(ex))as.character(ex[[1]])[1] else ''
  lhs<-if(hd=='<-')as.character(ex[[2]])[1] else ''
  if(hd %in% c('save_png','save_pdf','write.csv','writeLines'))next
  if(fig==1 && lhs %in% c('precip','eco','event','granger','fdr_summary','fdr_grid','fdr_lag'))next
  if(lhs=='ROOT'){env$ROOT<-ROOT;next}
  if(lhs=='OUT'){env$OUT<-file.path(DEST,paste0('Figure',fig));dir.create(env$OUT,recursive=TRUE,showWarnings=FALSE);next}
  eval(ex,envir=env)
  if(fig==2 && lhs=='side_profile_cont') {
    env$side_profile_original<-env$side_profile_cont
    env$side_profile_cont<-eval(quote(function(...) {
      last_sp <<- side_profile_original(...)
      last_sp
    }),envir=env)
  }
  funname<-if(fig==1)'map_block' else 'map_block_cont'
  if(lhs==funname) {
   fn<-env[[funname]];bo<-as.list(body(fn))
   for(j in seq_along(bo))if(is.call(bo[[j]])&&identical(bo[[j]][[1]],as.name('<-'))) {
    target<-as.character(bo[[j]][[2]])[1]
    if(target=='top')bo[[j]]<-if(fig==1)quote(top<-build_top(mp,side,side_core)) else quote(top<-build_top(mp,side,last_sp))
   }
   body(fn)<-as.call(bo);env[[funname]]<-fn
  }
  if(lhs==if(fig==1)'merged12' else 'main')break
 }
 base<-env[[if(fig==1)'merged12' else 'main']]
 W<-if(fig==1)210 else 225;H<-if(fig==1)225 else 160
 # Preserve the original export inset: removing it would change the maps.
 base<-ggdraw()+draw_plot(base,x=.018,y=.018,width=.964,height=.964)
 cairo_pdf(file.path(env$OUT,'measurement.pdf'),width=W/25.4,height=H/25.4)
 print(base);grid.force()
 yy<-function(x)get('robinson_project',envir=env,inherits=TRUE)(rep(0,length(x)),x)$y*.52
 tr<-scales::trans_new('Robinson latitude',yy,approxfun(yy(seq(-90,90,.01)),seq(-90,90,.01),rule=2),domain=c(-90,90))
 rows<-list()
 for(key in names(profiles)) {
  m<-pos[[paste0(key,'m')]];s<-pos[[paste0(key,'s')]]
  if(is.null(m))next
  sp<-profiles[[key]]+scale_y_continuous(trans=tr,limits=c(-90,90),breaks=NULL,expand=expansion(mult=0,add=.012))+
      labs(y=NULL)+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
  sp$coordinates<-coord_cartesian(xlim=sp$coordinates$limits$x,expand=TRUE,clip='on')
  stopifnot(max(abs(ggplot_build(sp)$layout$panel_params[[1]]$y.range-c(-.532,.532)))<1e-9)
  sg<-ggplotGrob(sp);p<-sg$layout[sg$layout$name=='panel',];sg$heights[p$t]<-unit(m[4]-m[2],'mm')
  below<-sum(convertHeight(sg$heights[(p$b+1):length(sg$heights)],'mm',valueOnly=TRUE))
  total<-sum(convertHeight(sg$heights,'mm',valueOnly=TRUE))
  base<-ggdraw(base)+draw_grob(sg,x=s[1]/W,y=(m[2]-below)/H,width=(s[3]-s[1])/W,height=total/H)
  rows[[key]]<-data.frame(panel=key,map_x0=m[1],map_y0=m[2],map_x1=m[3],map_y1=m[4])
 }
 dev.off()
 ggsave(file.path(env$OUT,paste0('Figure',fig,'_latitude_aligned.png')),base,width=W,height=H,units='mm',dpi=600,bg='white',device='png')
 ggsave(file.path(env$OUT,paste0('Figure',fig,'_latitude_aligned.pdf')),base,width=W,height=H,units='mm',device=cairo_pdf,bg='white')
 write.csv(do.call(rbind,rows),file.path(env$OUT,'MAP_PANEL_POSITIONS_MM.csv'),row.names=FALSE)
}
