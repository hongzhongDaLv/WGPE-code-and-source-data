# Purpose: Recompute wet/dry composites with an exact grid-month product
#          decomposition under the formal TOP100 definition.
# Inputs: TOP100 IWV/<z> monthly fields and ERA5 precipitation.
# Outputs: exact grid composites, regional summaries, closure QA and spatial CIs.
# Parameters: wet/dry are >90th/<10th percentiles of deseasoned, detrended
#             precipitation at each grid cell; component fields are linearly
#             trend-adjusted before calendar-month reference states are formed.
# Dependencies: base R.
# Overwrite policy: writes only exact_TOP100 outputs and never replaces legacy tables.

options(stringsAsFactors=FALSE,warn=1)
ROOT<-"__WGPE_PROJECT_ROOT__";source(file.path(ROOT,"scripts_TOP100_revision","18_coordinate_harmonization","coordinate_harmonization.R"));OUT<-Sys.getenv("TOP100_CLOSURE_WETDRY_OUT",unset=file.path(ROOT,"output_TOP100_revision","results","precip_wetdry","exact_TOP100"));dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
g<-9.80665
e<-new.env(parent=emptyenv());load(file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100","global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData"),envir=e)
lon<-as.numeric(e$lon);lat<-as.numeric(e$lat);dates<-as.Date(e$all_dates);nt<-length(dates);ncell<-length(lon)*length(lat);Imat<-matrix(as.numeric(e$iwv_st),nrow=ncell,ncol=nt);Zmat<-matrix(as.numeric(e$zbar_st),nrow=ncell,ncol=nt);rm(e);gc()
s<-readRDS(file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds"));p_aligned<-reorder_exact_grid(s$tp_native_m,s$lon,s$lat,lon,lat);Pmat<-matrix(as.numeric(p_aligned),nrow=ncell,ncol=nt);rm(s,p_aligned);gc()
cell<-expand.grid(lon=lon,lat=lat);cell$cell_id<-seq_len(ncell);cell$weight<-cos(cell$lat*pi/180)
month<-as.integer(format(dates,"%m"));td<-as.numeric(dates);tc<-td-mean(td)

adjust_and_anomaly<-function(m){
  an<-m
  for(mm in 1:12){ii<-which(month==mm);an[,ii]<-sweep(m[,ii,drop=FALSE],1,rowMeans(m[,ii,drop=FALSE],na.rm=TRUE),"-")}
  sl<-as.numeric(an%*%tc/sum(tc^2));adj<-m-outer(sl,tc)
  clim<-matrix(NA_real_,nrow(m),12)
  for(mm in 1:12)clim[,mm]<-rowMeans(adj[,month==mm,drop=FALSE],na.rm=TRUE)
  list(adjusted=adj,climatology=clim,anomaly=adj-clim[,month])
}

rows<-list();k<-0L;batch<-1000L
for(st in seq.int(1L,ncell,by=batch)){
  en<-min(ncell,st+batch-1L);id<-st:en;message("exact wet/dry cells ",st,"-",en)
  I<-Imat[id,,drop=FALSE]
  Z<-Zmat[id,,drop=FALSE]
  P<-Pmat[id,,drop=FALSE]
  ai<-adjust_and_anomaly(I);az<-adjust_and_anomaly(Z);ap<-adjust_and_anomaly(P)
  cov_clim<-matrix(NA_real_,nrow(I),12)
  prod_anom<-ai$anomaly*az$anomaly
  for(mm in 1:12)cov_clim[,mm]<-rowMeans(prod_anom[,month==mm,drop=FALSE],na.rm=TRUE)
  actual<-g*(ai$adjusted*az$adjusted-ai$climatology[,month]*az$climatology[,month]-cov_clim[,month])
  iwv_ch<-g*ai$anomaly*az$climatology[,month]
  z_ch<-g*az$anomaly*ai$climatology[,month]
  interaction<-g*(prod_anom-cov_clim[,month])
  q90<-apply(ap$anomaly,1,quantile,.9,na.rm=TRUE);q10<-apply(ap$anomaly,1,quantile,.1,na.rm=TRUE)
  for(ev in c("Wet","Dry")){
    mask<-if(ev=="Wet")ap$anomaly>q90 else ap$anomaly<q10
    mm<-function(x){x[!mask]<-NA_real_;rowMeans(x,na.rm=TRUE)}
    a<-mm(actual);ic<-mm(iwv_ch);zc<-mm(z_ch);it<-mm(interaction);pa<-mm(ap$anomaly)
    k<-k+1L;rows[[k]]<-data.frame(cell[id,],event=ev,actual=a,IWV_channel=ic,z_channel=zc,interaction_channel=it,closure_residual=a-ic-zc-it,precip_anomaly=pa,n_event=rowSums(mask,na.rm=TRUE),stringsAsFactors=FALSE)
  }
  rm(I,Z,P,ai,az,ap,cov_clim,prod_anom,actual,iwv_ch,z_ch,interaction);gc()
}
grid<-do.call(rbind,rows);write.csv(grid,file.path(OUT,"WET_DRY_EXACT_GRID_TOP100.csv"),row.names=FALSE)

band<-function(x)ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))));grid$band<-band(grid$lat);regions<-c("Global","S high latitude","S temperate","Tropics","N temperate","N high latitude")
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
block_boot<-function(d,cols,B=1000L){d$block<-paste(floor(((d$lon+180)%%360)/10),floor((d$lat+90)/10),sep="_");ids<-unique(d$block);den<-tapply(d$weight,d$block,sum)[ids];num<-lapply(cols,function(v)tapply(d$weight*d[[v]],d$block,sum)[ids]);set.seed(20260714+nrow(d));z<-matrix(NA_real_,B,length(cols));for(b in seq_len(B)){ii<-sample(seq_along(ids),length(ids),TRUE);z[b,]<-vapply(num,function(x)sum(x[ii])/sum(den[ii]),numeric(1))};z}
cols<-c("actual","IWV_channel","z_channel","interaction_channel","closure_residual","precip_anomaly");sumrows<-list();kk<-0L
for(ev in c("Wet","Dry"))for(reg in regions){d<-grid[grid$event==ev&(reg=="Global"|grid$band==reg),];boot<-block_boot(d,cols);for(j in seq_along(cols)){kk<-kk+1L;sumrows[[kk]]<-data.frame(event=ev,region=reg,metric=cols[j],estimate=wm(d[[cols[j]]],d$weight),lower95=quantile(boot[,j],.025,na.rm=TRUE),upper95=quantile(boot[,j],.975,na.rm=TRUE),n_grid=sum(is.finite(d[[cols[j]]])),mean_event_count=mean(d$n_event,na.rm=TRUE),uncertainty_method="10-degree spatial-block bootstrap, 1000 resamples")}}
summary<-do.call(rbind,sumrows);write.csv(summary,file.path(OUT,"WET_DRY_EXACT_SUMMARY_TOP100.csv"),row.names=FALSE)
qa<-data.frame(max_abs_grid_closure=max(abs(grid$closure_residual),na.rm=TRUE),global_wet_closure=summary$estimate[summary$event=="Wet"&summary$region=="Global"&summary$metric=="closure_residual"],global_dry_closure=summary$estimate[summary$event=="Dry"&summary$region=="Global"&summary$metric=="closure_residual"],wet_definition="> grid-cell 90th percentile",dry_definition="< grid-cell 10th percentile",precip_processing="calendar-month climatology removed and linear trend removed",component_processing="linear trend removed at factor level; interaction additionally centered by calendar-month IWV-z covariance",actual_definition="g*(IWV_adjusted*z_adjusted-calendar-month climatology of IWV_adjusted*z_adjusted)")
write.csv(qa,file.path(OUT,"WET_DRY_EXACT_CLOSURE_QA_TOP100.csv"),row.names=FALSE)
cat("Exact wet/dry TOP100 decomposition complete\n")
