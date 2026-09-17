# Layout-only v3: keep numerical data and archived drawing source untouched.
ROOT<-Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
D<-file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')
for(ex in parse(file.path(D,'run_original_with_new_results_colorbar.R'))){
 lhs<-if(is.call(ex)&&identical(ex[[1]],as.name('<-')))as.character(ex[[2]])[1] else ''
 if(lhs=='outroot'){outroot<-file.path(D,'LAYOUT_BOTTOM_v5');next}
 if(is.call(ex)&&as.character(ex[[1]])[1] %in% c('save_png','save_pdf'))next
 eval(ex)
}
e3<-fig4e+labs(title=NULL)+
 guides(fill=guide_colorbar(direction='horizontal',title.position='left',label.position='bottom',
   barwidth=unit(57,'mm'),barheight=unit(3,'mm'),ticks.colour='#111111',frame.colour='#555555'))+
 theme(plot.margin=margin(0,3,0,0),legend.margin=margin(3,0,0,0),legend.box.margin=margin(0,0,0,0))
f3<-fig4f_body
f3$layers[[1]]$geom_params$width<-.30
f3<-f3+labs(title=NULL)+theme(plot.margin=margin(0,3,0,0),panel.grid.major.x=element_blank())
core3<-ggdraw()+
 draw_plot(fig4_maps_only,x=0,y=.36,width=1,height=.64)+
 draw_plot(e3,x=-.065,y=.015,width=.49,height=.305)+
 draw_plot(f3,x=.505,y=.13,width=.36,height=.20)+
 draw_plot(fig4f_legend,x=.515,y=.04,width=.36,height=.065)+
 draw_label(expression(bold(e)~'Ecological synthesis by latitude band'),x=.020,y=.343,hjust=0,vjust=1,fontfamily=FONT,size=PANEL_TITLE_SIZE)+
 draw_label(expression(bold(f)~'Lag-3 temporal precedence'),x=.505,y=.343,hjust=0,vjust=1,fontfamily=FONT,size=PANEL_TITLE_SIZE)
full3<-plot_grid(core3,NULL,ncol=1,rel_heights=c(.94,.06))
# Direct vector-object export: no outer 0.964 shrink and no raster resampling.
ggsave(file.path(OUT,'Figure3_NO_VPD_layout_v5.png'),full3,width=215,height=205,units='mm',dpi=600,bg='white')
ggsave(file.path(OUT,'Figure3_NO_VPD_layout_v5.pdf'),full3,width=215,height=205,units='mm',device=cairo_pdf,bg='white')
cat('DONE v5 Figure2 colorbar function; direct export without outer shrink\n')
