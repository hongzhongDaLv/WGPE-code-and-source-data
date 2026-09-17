# Isolated alternative to stacked columns; estimates and signs are unchanged.
.libPaths(.libPaths())
library(ggplot2);library(cowplot);library(grid)
ROOT<-Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D<-file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
OUT<-file.path(D,'CHANNEL_DOT_PREVIEW');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
g<-read.csv(file.path(D,'FULL_FIGURE_STYLE_PREVIEW_v3/Figure1/standard_latband_WGPE_trend_decomposition.csv'))
d<-readRDS(file.path(ROOT,'output_TOP100_revision/derived_data/figure_source_data/Figure2_TOP100_adapter.rds'))$enso_decomposition
regions<-c('Global','N high latitude','N temperate','Tropics','S temperate','S high latitude')
ys<-setNames(c(6.6,5,4,3,2,1),regions)
cols<-c('IWV'='#2C7BB6','<z>'='#D66A55','Interaction'='#8E2A62')
long<-function(df,names,divisor=1) {
 do.call(rbind,lapply(seq_along(names),function(i)data.frame(region=df$region,
  channel=factor(names(cols)[i],levels=names(cols)),value=df[[names[i]]]/divisor,
  y=unname(ys[df$region])+c(.21,0,-.21)[i])))
}
a<-long(g,c('iwv_channel_J_m2_yr','z_channel_J_m2_yr','interaction_channel_J_m2_yr'))
b<-long(d,c('IWV_channel_J_m2','z_channel_J_m2','interaction_channel_J_m2'))
b$value<-b$value/1000
make_plot<-function(dat,title,xlab) {
 ggplot(dat,aes(value,y,colour=channel))+
 geom_hline(yintercept=c(1,2,3,4,5,6.6),colour='#EFEFEF',linewidth=.18)+
 geom_vline(xintercept=0,colour='#777777',linewidth=.3)+
 geom_point(shape=16,size=1.8)+
 scale_colour_manual(values=cols,name=NULL)+
 scale_y_continuous(breaks=unname(ys),labels=regions,limits=c(.5,7.1),expand=c(0,0))+
 scale_x_continuous(expand=expansion(mult=c(.08,.20)))+
 labs(title=title,x=xlab,y=NULL)+theme_classic(base_size=9.6,base_family='sans')+
 theme(plot.title=element_text(size=9.6,hjust=0,margin=margin(b=5)),
       axis.text=element_text(size=9.6,colour='#202020'),axis.title=element_text(size=9.6),
       axis.line.y=element_blank(),axis.ticks.y=element_blank(),
       axis.line.x=element_line(linewidth=.25,colour='#777777'),
       axis.ticks=element_line(linewidth=.22),legend.position='top',legend.justification='left',
       legend.text=element_text(size=9.6),legend.key.width=unit(3,'mm'),legend.key.height=unit(2,'mm'),
       legend.spacing.x=unit(1,'mm'),legend.margin=margin(0,0,3,0),
       plot.margin=margin(3,4,3,2))+guides(colour=guide_legend(nrow=1,override.aes=list(size=1.8)))
}
p1<-make_plot(a,expression(bold(g)~'Regional W-GPE composition'),expression('Contribution (J m'^{-2}~'yr'^{-1}*')'))
p2<-make_plot(b,expression(bold(d)~'ENSO decomposition'),expression('Channel anomaly (10'^3~'J m'^{-2}*')'))
ggsave(file.path(OUT,'Fig1g_channel_dots.png'),p1,width=75,height=80,units='mm',dpi=600,bg='white')
ggsave(file.path(OUT,'Fig2d_channel_dots.png'),p2,width=75,height=80,units='mm',dpi=600,bg='white')
both<-plot_grid(p1,p2,nrow=1,align='h',axis='tb')
ggsave(file.path(OUT,'Channel_dot_comparison.png'),both,width=150,height=80,units='mm',dpi=600,bg='white')
write.csv(a,file.path(OUT,'Fig1g_values.csv'),row.names=FALSE)
write.csv(b,file.path(OUT,'Fig2d_values.csv'),row.names=FALSE)
