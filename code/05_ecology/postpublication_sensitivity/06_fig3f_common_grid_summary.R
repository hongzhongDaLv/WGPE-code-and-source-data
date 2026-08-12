ROOT<-"__WGPE_PROJECT_ROOT__"
CANONICAL<-file.path(ROOT,"CANONICAL_WGPE_TOP100_v2")
OUT<-file.path(CANONICAL,"08_REPRODUCTION_OUTPUT","Figure3f")
cls<-read.csv(file.path(OUT,"FIG3F_LAG3_CLASSIFICATION_GRID.csv"),check.names=FALSE)
eff<-read.csv(file.path(ROOT,"output_TOP100_revision","results","temporal_precedence","TOP100_effect_sizes","ECOLOGY_TEMPORAL_EFFECT_SIZE_GRID_TOP100.csv"),check.names=FALSE)
eff<-eff[eff$lag==3&eff$persistent&eff$fit_status==1,]
classes<-c("forward-only","reverse-only","bidirectional","none");rows<-list();k<-0L
for(tg in c("GPP","LAI"))for(me in c("F","HAC3","HAC12")){
 d<-cls[cls$target==tg&cls$method==me,];keep<-unique(eff$grid[eff$target==tg]);d<-d[d$grid_id%in%keep,];w<-cos(d$lat*pi/180)
 for(cl in classes){k<-k+1L;rows[[k]]<-data.frame(target=tg,method=me,class=cl,area_fraction_pct=100*sum(w[d$class==cl])/sum(w),n_common_grid=nrow(d))}
}
write.csv(do.call(rbind,rows),file.path(OUT,"FIG3F_LAG3_COMMON_CV_GRID_CLASSIFICATION.csv"),row.names=FALSE)

set.seed(20260722);B<-2000L;boot_rows<-list();k<-0L
for(tg in c("GPP","LAI")){
 d<-cls[cls$target==tg&cls$method=="F",];w<-cos(d$lat*pi/180);bid<-paste(floor(((d$lon+180)%%360)/10),floor((d$lat+90)/10),sep="_");ids<-unique(bid)
 for(cl in classes){est<-100*sum(w[d$class==cl])/sum(w);bv<-rep(NA_real_,B);for(b in seq_len(B)){ss<-sample(ids,length(ids),TRUE);take<-unlist(lapply(ss,function(x)which(bid==x)),use.names=FALSE);bv[b]<-100*sum(w[take][d$class[take]==cl])/sum(w[take])};k<-k+1L;boot_rows[[k]]<-data.frame(target=tg,class=cl,estimate=est,lower95=quantile(bv,.025),upper95=quantile(bv,.975),n_grid=nrow(d),method="ordinary nested F-test with time-preserving lag construction; 10-degree spatial-block bootstrap")}
}
write.csv(do.call(rbind,boot_rows),file.path(OUT,"FIG3F_SOURCE_FINAL.csv"),row.names=FALSE)
