options(stringsAsFactors=FALSE,warn=1)
ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
OUT <- Sys.getenv('WGPE_ECOLOGY_OUTPUT', unset=file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904'))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
setwd(ROOT)
sink(file.path(OUT,'RUN_LOG.txt'),split=TRUE)
library(ggplot2)
preprocess_monthly_cpp <- function(X,month,time){
 for(m in 1:12){ii<-which(month==m);X[,ii]<-X[,ii,drop=FALSE]-rowMeans(X[,ii,drop=FALSE],na.rm=TRUE)}
 for(i in seq_len(nrow(X))){ok<-is.finite(X[i,]);if(sum(ok)<10){X[i,]<-NA;next};tt<-time[ok]-mean(time[ok]);y<-X[i,ok];b<-sum(tt*(y-mean(y)))/sum(tt^2);y<-y-b*tt;y<-y-mean(y);ss<-sd(y);X[i,ok]<-if(is.finite(ss)&&ss>0)y/ss else NA_real_};X
}
partial_two_canonical_cpp <- function(Y,X1,X2,controls){
 out<-matrix(NA_real_,nrow(Y),10,dimnames=list(NULL,c('partial_r_X1','p_X1','neff_X1','partial_r_X2','p_X2','neff_X2','n_complete','rank_controls','condition_number_controls','fit_status')))
 for(i in seq_len(nrow(Y))){Z<-cbind(1,do.call(cbind,lapply(controls,function(x)x[i,])));ok<-is.finite(Y[i,])&is.finite(X1[i,])&is.finite(X2[i,])&rowSums(!is.finite(Z))==0;n<-sum(ok);out[i,7]<-n;out[i,10]<-0;if(n<36)next;Z<-Z[ok,,drop=FALSE];q<-qr(Z);if(q$rank<ncol(Z))next;res<-qr.resid(q,cbind(Y[i,ok],X1[i,ok],X2[i,ok]));out[i,8]<-q$rank;out[i,9]<-kappa(Z);out[i,10]<-1
  for(j in 1:2){a<-res[,1];b<-res[,j+1];r<-cor(a,b);aa<-cor(head(a,-1),tail(a,-1));bb<-cor(head(b,-1),tail(b,-1));if(!is.finite(aa))aa<-0;if(!is.finite(bb))bb<-0;ab<-max(-.99,min(.99,aa))*max(-.99,min(.99,bb));ne<-max(length(controls)+4,n*(1-ab)/(1+ab));df<-max(2,ne-length(controls)-2);p<-2*pt(-abs(r*sqrt(df/max(1e-12,1-r*r))),df);out[i,(j-1)*3+1:3]<-c(r,p,ne)}
 };out
}
SRC <- 'output_TOP100_revision/FINAL_COORDINATE_FIXED_CLOSURE'
controls <- c('precipitation','T2m','SSRD','RZSM','CO2','VPD')
band <- function(x)ifelse(x>66.5,'N high latitude',ifelse(x>23.5,'N temperate',ifelse(x>=-23.5,'Tropics',ifelse(x>=-66.5,'S temperate','S high latitude'))))
boot <- function(v,w,lon,lat){b<-paste(floor((lon%%360)/10),floor((lat+90)/10)); n<-tapply(v*w,b,sum); d<-tapply(w,b,sum); set.seed(20260904); quantile(replicate(1000,{s<-sample(seq_along(n),length(n),TRUE);sum(n[s])/sum(d[s])}),c(.025,.975),names=FALSE)}
grid <- list(); summaries <- list(); qa <- list()
locked <- read.csv(file.path(SRC,'results/ecology/ECOLOGY_PARTIAL_SUMMARY_CANONICAL.csv'))
for(target in c('GPP','LAI')){
 cat('\nLoading ',target,'\n');D<-readRDS(file.path(SRC,'derived_data/ecology',paste0(target,'_CANONICAL_INPUTS.rds')))
 vars<-c('y','wgpe_absolute','z_absolute',controls)
 for(v in vars){D[[v]]<-preprocess_monthly_cpp(D[[v]],as.integer(format(D$dates,'%m')),as.numeric(D$dates));cat('Preprocessed ',v,'\n')}
 # Freeze full-model complete cases before removing VPD.
 valid<-Reduce('&',lapply(vars,function(v)is.finite(D[[v]])))
 for(v in vars)D[[v]][!valid]<-NA_real_
 rr<-list(With_VPD=partial_two_canonical_cpp(D$y,D$wgpe_absolute,D$z_absolute,lapply(controls,function(v)D[[v]])),Without_VPD=partial_two_canonical_cpp(D$y,D$wgpe_absolute,D$z_absolute,lapply(head(controls,-1),function(v)D[[v]])))
 stopifnot(identical(rr[[1]][,'n_complete'],rr[[2]][,'n_complete']))
 meta<-D$meta;keep<-meta$target_valid_min36 & meta$persistent_vegetated_domain & rr[[1]][,'fit_status']==1 & rr[[2]][,'fit_status']==1
 for(pred in c('W-GPE','<z>')){
  j<-if(pred=='W-GPE')1 else 2
  for(model in names(rr)){
   a<-rr[[model]];g<-data.frame(target=target,predictor=pred,model=model,lon=meta$lon[keep],lat=meta$lat[keep],r=a[keep,paste0('partial_r_X',j)],p=a[keep,paste0('p_X',j)],n=a[keep,'n_complete'])
   g$q_BH<-p.adjust(g$p,'BH');g$q_BY<-p.adjust(g$p,'BY');g$region<-band(g$lat);grid[[length(grid)+1]]<-g
   for(reg in c('Global',unique(g$region))){z<-g[reg=='Global'|g$region==reg,];w<-cos(z$lat*pi/180);ci<-boot(z$r,w,z$lon,z$lat);summaries[[length(summaries)+1]]<-data.frame(target,predictor=pred,model,region=reg,n_grid=nrow(z),mean_r=weighted.mean(z$r,w),lower95=ci[1],upper95=ci[2],BH_positive_pct=100*sum(w[z$q_BH<=.05 & z$r>0])/sum(w),BH_negative_pct=100*sum(w[z$q_BH<=.05 & z$r<0])/sum(w),BY_positive_pct=100*sum(w[z$q_BY<=.05 & z$r>0])/sum(w))}
   if(model=='With_VPD'){old<-subset(locked,target==g$target[1]&predictor==pred&height_reference=='absolute'&model=='Model3'&domain=='persistent_vegetated_domain'&region=='Global');qa[[length(qa)+1]]<-data.frame(target,predictor=pred,recomputed=weighted.mean(g$r,cos(g$lat*pi/180)),locked=old$mean_partial_r)}
  }
 }
 rm(D);gc()
}
g<-do.call(rbind,grid);s<-do.call(rbind,summaries);q<-do.call(rbind,qa);q$difference<-q$recomputed-q$locked
write.csv(g,file.path(OUT,'GRID_RESULTS.csv'),row.names=FALSE);write.csv(s,file.path(OUT,'SUMMARY_RESULTS.csv'),row.names=FALSE);write.csv(q,file.path(OUT,'BASELINE_QA.csv'),row.names=FALSE)
stopifnot(all(abs(q$difference)<1e-8))
diffs<-list()
for(t in unique(g$target))for(pr in unique(g$predictor)){a<-subset(g,target==t&predictor==pr&model=='Without_VPD');b<-subset(g,target==t&predictor==pr&model=='With_VPD');stopifnot(identical(a$lon,b$lon),identical(a$lat,b$lat),identical(a$n,b$n));d<-a$r-b$r;w<-cos(a$lat*pi/180);ci<-boot(d,w,a$lon,a$lat);diffs[[length(diffs)+1]]<-data.frame(target=t,predictor=pr,delta_r_without_minus_with=weighted.mean(d,w),lower95=ci[1],upper95=ci[2])}
delta<-do.call(rbind,diffs);write.csv(delta,file.path(OUT,'PAIRED_DIFFERENCES.csv'),row.names=FALSE)
gp<-subset(g,model=='Without_VPD');gp$lon<-((gp$lon+180)%%360)-180
geom_raster <- function(...) geom_tile(...,width=1,height=1)
p<-ggplot(gp,aes(lon,lat))+geom_raster(aes(fill=r))+facet_grid(target~predictor)+coord_quickmap(xlim=c(-180,180),ylim=c(-60,85),expand=FALSE)+scale_fill_gradient2(low='#a6611a',mid='white',high='#018571',limits=c(-1,1),name='Partial r')+labs(title='Ecological associations without VPD adjustment',subtitle='Deseasonalized and detrended; controls: precipitation, T2m, SSRD, RZSM and CO2',x=NULL,y=NULL)+theme_bw(base_size=11)+theme(panel.grid=element_blank(),legend.position='bottom')
ggsave(file.path(OUT,'ECOLOGY_NO_VPD_MAPS.png'),p,width=12,height=7,dpi=220);ggsave(file.path(OUT,'ECOLOGY_NO_VPD_MAPS.pdf'),p,width=12,height=7)
sp<-subset(s,region=='Global');p2<-ggplot(sp,aes(model,mean_r,colour=model))+geom_point(size=3)+geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.12)+facet_grid(target~predictor,scales='free_y')+geom_hline(yintercept=0,linetype=2,colour='grey60')+labs(x=NULL,y='Area-weighted mean partial r',title='Same-sample comparison: with and without VPD')+theme_bw(base_size=12)+theme(legend.position='none')
ggsave(file.path(OUT,'VPD_PAIRED_COMPARISON.png'),p2,width=9,height=6,dpi=200)
if(requireNamespace('openxlsx',quietly=TRUE))openxlsx::write.xlsx(list(Global=subset(s,region=='Global'),All_regions=s,Paired_differences=delta,Baseline_QA=q,Grid_results=g),file.path(OUT,'ECOLOGY_NO_VPD_MASTER.xlsx'),overwrite=TRUE)
writeLines(c('# No-VPD ecological recalculation','', 'Status: baseline numerical reproduction PASS.','', 'Controls retained: precipitation, T2m, SSRD, RZSM, CO2. VPD alone removed.','Full-model complete-case months and persistent vegetation mask held fixed. Deseasonalization, linear detrending and canonical AR(1)/BH/BY algorithm retained.','Spatial intervals: paired 10-degree blocks, 1000 replicates.','',paste(capture.output(print(sp,row.names=FALSE)),collapse='\n'),'',paste(capture.output(print(delta,row.names=FALSE)),collapse='\n'),'','This run recalculates concurrent partial associations only. Temporal-precedence panels, OOF performance and circular-shift field tests have not been recalculated; old results for those analyses are not relabelled as no-VPD results.','The inherited significance algorithm estimates AR(1) on complete-case residual sequences; it does not explicitly preserve calendar gaps.','No manuscript, original figure or canonical result overwritten.'),file.path(OUT,'ECOLOGY_NO_VPD_MASTER_REPORT.md'))
cat('\nCOMPLETE\n');sink()
