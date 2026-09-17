# Preserve v5 map viewports exactly; replace only companion profiles.
ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D <- file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
for(ex in parse(file.path(D,'layout_bottom_v5.R'))) {
  head <- if(is.call(ex)) as.character(ex[[1]])[1] else ''
  if(!head %in% c('ggsave','cat')) eval(ex)
}
OUT <- file.path(D,'LATITUDE_ALIGNMENT_v7')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
positions <- new.env(); companions <- list()
drawDetails.record_panel <- function(x,recording=TRUE) {
  a <- deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE)
  b <- deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
  positions[[x$key]] <- c(x0=a$x,y0=a$y,x1=b$x,y1=b$y)*25.4
}
lat_y <- function(x) robinson_project(rep(0,length(x)),x)$y*.52
lat_trans <- scales::trans_new('Robinson latitude',lat_y,
  approxfun(lat_y(seq(-90,90,.01)),seq(-90,90,.01),rule=2),domain=c(-90,90))
map_block <- function(field,grid,q,lim,pal,tag,title_expr,cbar_title_expr,cbar_ticks,
 fdr=FALSE,accuracy=.1,profile_accuracy=accuracy,profile_colour=COL$pos,subtitle_text=NULL) {
  b <- bundle_from_field(field,grid,q)
  mp <- map_plot(b,lim,pal,tag,title_expr,fdr,subtitle_text)
  sp <- map_profile_companion(b,lim,cbar_title_expr,profile_colour,profile_accuracy)
  # Keep the original side grob's layout dimensions for the v5 alignment pass.
  side <- plot_grid(NULL,sp,NULL,ncol=1,rel_heights=c(.12,.76,.12))
  sg <- ggplotGrob(side); sg$grobs <- lapply(sg$grobs,function(x)nullGrob())
  mg <- ggplotGrob(mp); pp <- mg$layout[mg$layout$name=='panel',]
  mg <- gtable::gtable_add_grob(mg,grob(key=tag,cl='record_panel'),t=pp$t,l=pp$l,b=pp$b,r=pp$r)
  top <- plot_grid(mg,NULL,sg,nrow=1,rel_widths=c(.780,.020,.200),align='h',axis='tb')
  sp <- sp+scale_y_continuous(trans=lat_trans,limits=c(-90,90),breaks=NULL,
       expand=expansion(mult=0,add=.012))+labs(y=NULL)+
       theme(axis.title.y=element_blank(),axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),axis.line.y=element_blank())
  sp$coordinates <- coord_cartesian(xlim=sp$coordinates$limits$x,expand=TRUE,clip='on')
  stopifnot(max(abs(ggplot_build(sp)$layout$panel_params[[1]]$y.range-c(-.532,.532)))<1e-10)
  companions[[tag]] <<- ggplotGrob(sp)
  cb <- plot_grid(NULL,colorbar(lim,pal,cbar_title_expr,cbar_ticks,accuracy),NULL,
                 nrow=1,rel_widths=c(.04,.72,.24))
  ggdraw()+draw_plot(top,x=0,y=.24,width=.965,height=.76)+
     draw_plot(cb,x=0,y=.01,width=.965,height=.16)
}
ECO_PROFILE_CALL <- 0L
fig4a <- eco_map('GPP','wgpe','a','W-GPE-GPP partial correlation')
fig4b <- eco_map('GPP','zbar','b','<z>-GPP partial correlation')
fig4c <- eco_map('LAI','wgpe','c','W-GPE-LAI partial correlation')
fig4d <- eco_map('LAI','zbar','d','<z>-LAI partial correlation')
fig4_maps_only <- plot_grid(fig4a,fig4b,fig4c,fig4d,ncol=2)
for(ex in parse(file.path(D,'layout_bottom_v5.R'))) {
 lhs <- if(is.call(ex)&&identical(ex[[1]],as.name('<-')))as.character(ex[[2]])[1] else ''
 if(lhs %in% c('core3','full3'))eval(ex)
}
# Measure actual map panels at the final output size, without moving them.
cairo_pdf(file.path(OUT,'measurement.pdf'),width=215/25.4,height=205/25.4)
print(full3); grid.force()
base <- full3
for(tag in c('a','b','c','d')) {
  pos <- positions[[tag]]; stopifnot(length(pos)==4)
  sg <- companions[[tag]]; pp <- sg$layout[sg$layout$name=='panel',]
  sg$heights[pp$t] <- unit(pos['y1']-pos['y0'],'mm')
  below <- sum(convertHeight(sg$heights[(pp$b+1):length(sg$heights)],'mm',valueOnly=TRUE))
  total <- sum(convertHeight(sg$heights,'mm',valueOnly=TRUE))
  x <- if(tag %in% c('a','c')) .965*.8/2 else .5+.965*.8/2
  base <- ggdraw(base)+draw_grob(sg,x=x,y=(pos['y0']-below)/205,width=.965*.2/2,height=total/205)
}
dev.off()
ggsave(file.path(OUT,'Figure3_NO_VPD_latitude_aligned_v7.png'),base,width=215,height=205,units='mm',dpi=600,bg='white')
ggsave(file.path(OUT,'Figure3_NO_VPD_latitude_aligned_v7.pdf'),base,width=215,height=205,units='mm',device=cairo_pdf,bg='white')
write.csv(do.call(rbind,lapply(c('a','b','c','d'),function(k)data.frame(panel=k,t(positions[[k]])))),
          file.path(OUT,'MAP_PANEL_POSITIONS_MM.csv'),row.names=FALSE)
