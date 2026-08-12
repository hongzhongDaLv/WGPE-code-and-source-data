options(stringsAsFactors=FALSE, warn=1)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
Sys.setenv(OMP_NUM_THREADS="8")
ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- Sys.getenv("TOP100_CLOSURE_HYDRO_INCREMENT_OUT",file.path(ROOT,"output_TOP100_revision","results","incremental_value","hydroclimate_TOP100"))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
source(file.path(ROOT,"scripts_TOP100_revision","18_coordinate_harmonization","coordinate_harmonization.R"))
if(!requireNamespace("Rcpp",quietly=TRUE)||!requireNamespace("RcppArmadillo",quietly=TRUE)) stop("Rcpp/RcppArmadillo required")
Rcpp::sourceCpp(file.path(ROOT,"scripts_core_validation","core_validation_cv.cpp"),rebuild=TRUE,showOutput=FALSE)

core_file <- file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100","global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file <- file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
e <- new.env(parent=emptyenv()); load(core_file,envir=e)
lon <- as.numeric(e$lon); lat <- as.numeric(e$lat); dates <- as.Date(e$all_dates)
nlon <- length(lon); nlat <- length(lat); nt <- length(dates); ncell <- nlon*nlat
cell <- expand.grid(lon=lon,lat=lat); cell$cell_id <- seq_len(ncell); cell$weight <- cos(cell$lat*pi/180)
to_mat <- function(a) matrix(as.numeric(a),nrow=ncell,ncol=nt)
m_i <- to_mat(e$iwv_st); m_z <- to_mat(e$zbar_st); m_w <- to_mat(e$wgpe_st); rm(e);gc()
s <- readRDS(single_file)
m_p <- to_mat(reorder_exact_grid(s$tp_native_m,s$lon,s$lat,lon,lat)); rm(s);gc()
mon <- as.integer(format(dates,"%m")); yr <- as.integer(format(dates,"%Y")); tdec <- yr+(mon-.5)/12

prep <- function(m) {
  good <- rowSums(is.finite(m))==ncol(m)
  m[!good,] <- NA_real_
  for(mm in 1:12){id<-which(mon==mm);cl<-rowMeans(m[,id,drop=FALSE],na.rm=TRUE);m[,id]<-sweep(m[,id,drop=FALSE],1,cl,"-")}
  tc <- tdec-mean(tdec); mu<-rowMeans(m,na.rm=TRUE); yc<-sweep(m,1,mu,"-")
  sl <- as.numeric(yc%*%tc/sum(tc^2)); m <- m-outer(sl,tc)
  mu<-rowMeans(m,na.rm=TRUE); m<-sweep(m,1,mu,"-"); sdv<-sqrt(rowSums(m^2,na.rm=TRUE)/pmax(1,rowSums(is.finite(m))-1))
  ok<-is.finite(sdv)&sdv>0;m[ok,]<-m[ok,,drop=FALSE]/sdv[ok];m[!ok,]<-NA_real_;m
}
message("Preprocessing final monthly matrices")
m_i<-prep(m_i);m_z<-prep(m_z);m_w<-prep(m_w);m_p<-prep(m_p);m_int<-m_i*m_z;gc()
fold <- as.integer((yr-min(yr))%/%5)+1L
gaps <- c(0L,3L)

run_linear <- function(gap,batch=2500L){
  z<-list();kk<-0L
  for(st in seq.int(1L,ncell,by=batch)){en<-min(ncell,st+batch-1L);id<-st:en;kk<-kk+1L
    message("linear gap ",gap," batch ",kk,"/",ceiling(ncell/batch))
    rr<-cv_linear_nested_cpp(m_p[id,,drop=FALSE],m_i[id,,drop=FALSE],m_z[id,,drop=FALSE],m_int[id,,drop=FALSE],m_w[id,,drop=FALSE],list(),fold,gap,FALSE)
    z[[kk]]<-data.frame(cell[id,],gap_months=gap,rr,check.names=FALSE)
  }
  do.call(rbind,z)
}
lin_cache<-file.path(OUT,"PRECIP_MODEL_PERFORMANCE_grid.csv")
if(file.exists(lin_cache)){message("Using completed linear-model cache");lin<-read.csv(lin_cache,check.names=FALSE)}else{lin <- do.call(rbind,lapply(gaps,run_linear));write.csv(lin,lin_cache,row.names=FALSE)}

band <- function(x){ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))))}
cell$climate_band<-band(cell$lat)
regions<-c("Global","N high latitude","N temperate","Tropics","S temperate","S high latitude")
mask_reg<-function(d,r) if(r=="Global")rep(TRUE,nrow(d)) else d$climate_band==r
block_ci<-function(x,w,lo,la,B=1000L){ok<-is.finite(x)&is.finite(w);if(sum(ok)<5)return(c(NA,NA));x<-x[ok];w<-w[ok];bid<-paste(floor(((lo[ok]+180)%%360)/10),floor((la[ok]+90)/10),sep="_");bn<-tapply(x*w,bid,sum);bd<-tapply(w,bid,sum);ids<-names(bn);set.seed(20260713+length(x));v<-rep(NA_real_,B);for(b in seq_len(B)){s<-sample(ids,length(ids),replace=TRUE);v[b]<-sum(bn[s])/sum(bd[s])};quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)}
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}

model_perf<-list();incsum<-list();k1<-k2<-0L
for(gap in gaps){d<-lin[lin$gap_months==gap,];d$climate_band<-band(d$lat)
  for(m in 0:3)for(r in regions){ii<-mask_reg(d,r);w<-d$weight[ii]
    for(metric in c("cv_r2","rmse","mae","pred_cor")){v<-d[[paste0("M",m,"_",metric)]][ii];ci<-block_ci(v,w,d$lon[ii],d$lat[ii]);k1<-k1+1L;model_perf[[k1]]<-data.frame(module="precipitation",gap_months=gap,model=paste0("M",m),metric=metric,region=r,estimate=wm(v,w),lower95=ci[1],upper95=ci[2],median=median(v,na.rm=TRUE),n_grid=sum(is.finite(v)),n_month=nt)}
  }
  pairs<-list(z_increment=c(0,1),interaction_increment=c(1,2),total_vertical_position_increment=c(0,2))
  pcols<-c(z_increment="increment_z_p",interaction_increment="increment_interaction_p",total_vertical_position_increment="increment_total_p")
  for(nm in names(pairs)){a<-pairs[[nm]][1];b<-pairs[[nm]][2];dv<-d[[paste0("M",b,"_cv_r2")]]-d[[paste0("M",a,"_cv_r2")]];dr<-d[[paste0("M",b,"_rmse")]]-d[[paste0("M",a,"_rmse")]];dm<-d[[paste0("M",b,"_mae")]]-d[[paste0("M",a,"_mae")]];dc<-d[[paste0("M",b,"_pred_cor")]]-d[[paste0("M",a,"_pred_cor")]];pv<-d[[pcols[nm]]];qv<-rep(NA_real_,length(pv));ok<-is.finite(pv);qv[ok]<-p.adjust(pv[ok],"BH")
    if(gap==3){lin$delta_cv_r2[lin$gap_months==gap & match(lin$cell_id,d$cell_id)>0]<-NA_real_}
    for(r in regions){ii<-mask_reg(d,r);w<-d$weight[ii];ci<-block_ci(dv[ii],w,d$lon[ii],d$lat[ii]);den<-sum(w[is.finite(dv[ii])]);k2<-k2+1L;incsum[[k2]]<-data.frame(module="precipitation",gap_months=gap,increment=nm,region=r,delta_CV_R2_mean=wm(dv[ii],w),delta_CV_R2_median=median(dv[ii],na.rm=TRUE),lower95=ci[1],upper95=ci[2],delta_RMSE_mean=wm(dr[ii],w),delta_MAE_mean=wm(dm[ii],w),delta_prediction_correlation_mean=wm(dc[ii],w),positive_delta_CV_R2_area_pct=100*sum(w[is.finite(dv[ii])&dv[ii]>0])/den,significant_positive_delta_CV_R2_area_pct=100*sum(w[is.finite(qv[ii])&qv[ii]<=.05&dv[ii]>0])/den,improved_RMSE_area_pct=100*sum(w[is.finite(dr[ii])&dr[ii]<0])/den,n_grid=sum(is.finite(dv[ii])),n_month=nt)}
  }
}
model_perf<-do.call(rbind,model_perf);incsum<-do.call(rbind,incsum)
write.csv(model_perf,file.path(OUT,"PRECIP_MODEL_PERFORMANCE.csv"),row.names=FALSE)
write.csv(incsum,file.path(OUT,"PRECIP_INCREMENT_SUMMARY.csv"),row.names=FALSE)

coll<-lin[,c("cell_id","lon","lat","gap_months","predictor_cor_I_Z","vif_I_Z","condition_number","rank","coef_sd_z","coef_sd_interaction","fit_status")]
write.csv(coll,file.path(OUT,"COLLINEARITY_DIAGNOSTICS_precip.csv"),row.names=FALSE)
cvdiag<-data.frame(module="precipitation",gap_months=gaps,n_fold=max(fold),fold_years=5,n_grid=sapply(gaps,function(g)sum(lin$gap_months==g & lin$fit_status==1,na.rm=TRUE)),rank_or_fit_warning=sapply(gaps,function(g)sum(lin$gap_months==g & lin$fit_status==2,na.rm=TRUE)),notes="shared continuous test folds; training excludes the test block and specified adjacent-month gap")
write.csv(cvdiag,file.path(OUT,"CV_DIAGNOSTICS_precip.csv"),row.names=FALSE)

run_logistic_event<-function(event_name,Y,batch=1000L){all<-list();kk<-0L
  for(gap in gaps)for(st in seq.int(1L,ncell,by=batch)){en<-min(ncell,st+batch-1L);id<-st:en;kk<-kk+1L;message(event_name," gap ",gap," cells ",st,"-",en);rr<-cv_logistic_nested_cpp(Y[id,,drop=FALSE],m_i[id,,drop=FALSE],m_z[id,,drop=FALSE],m_int[id,,drop=FALSE],fold,gap);all[[kk]]<-data.frame(cell[id,],event=event_name,gap_months=gap,rr,check.names=FALSE)}
  do.call(rbind,all)
}
q90<-apply(m_p,1,quantile,.9,na.rm=TRUE);wet<-1*(m_p>q90);wet[!is.finite(m_p)]<-NA_real_;wetres<-run_logistic_event("wet",wet);rm(wet);gc()
q10<-apply(m_p,1,quantile,.1,na.rm=TRUE);dry<-1*(m_p<q10);dry[!is.finite(m_p)]<-NA_real_;dryres<-run_logistic_event("dry",dry);rm(dry);gc()
loggrid<-rbind(wetres,dryres);write.csv(loggrid,file.path(OUT,"WET_DRY_MODEL_PERFORMANCE_grid.csv"),row.names=FALSE)

summ_log<-function(dd,event){rows<-list();k<-0L;dd$climate_band<-band(dd$lat);pairs<-list(z_increment=c(0,1),interaction_increment=c(1,2),total_vertical_position_increment=c(0,2));for(gap in gaps){d<-dd[dd$gap_months==gap,];for(nm in names(pairs)){a<-pairs[[nm]][1];b<-pairs[[nm]][2];da<-d[[paste0("M",b,"_auc")]]-d[[paste0("M",a,"_auc")]];db<-d[[paste0("M",a,"_brier")]]-d[[paste0("M",b,"_brier")]];dl<-d[[paste0("M",a,"_logloss")]]-d[[paste0("M",b,"_logloss")]];for(r in regions){ii<-mask_reg(d,r);w<-d$weight[ii];ci<-block_ci(da[ii],w,d$lon[ii],d$lat[ii]);den<-sum(w[is.finite(da[ii])]);k<-k+1L;rows[[k]]<-data.frame(module=event,gap_months=gap,increment=nm,region=r,delta_AUC=wm(da[ii],w),lower95=ci[1],upper95=ci[2],delta_Brier_improvement=wm(db[ii],w),delta_logloss_improvement=wm(dl[ii],w),positive_delta_AUC_area_pct=100*sum(w[is.finite(da[ii])&da[ii]>0])/den,improved_Brier_area_pct=100*sum(w[is.finite(db[ii])&db[ii]>0])/den,n_grid=sum(is.finite(da[ii])),n_event=sum(d$n_event[ii],na.rm=TRUE),n_month=nt)}}};do.call(rbind,rows)}
wets<-summ_log(wetres,"wet");drys<-summ_log(dryres,"dry")
write.csv(wets,file.path(OUT,"WET_INCREMENT_SUMMARY.csv"),row.names=FALSE)
write.csv(drys,file.path(OUT,"DRY_INCREMENT_SUMMARY.csv"),row.names=FALSE)

maps<-merge(lin[lin$gap_months==3,c("cell_id","lon","lat","M0_cv_r2","M1_cv_r2","M2_cv_r2","M0_rmse","M1_rmse","M2_rmse")],wetres[wetres$gap_months==3,c("cell_id","M0_auc","M1_auc","M2_auc")],by="cell_id",all.x=TRUE,suffixes=c("","_wet"));names(maps)[names(maps)%in%c("M0_auc","M1_auc","M2_auc")]<-paste0(c("M0","M1","M2"),"_auc_wet");maps<-merge(maps,dryres[dryres$gap_months==3,c("cell_id","M0_auc","M1_auc","M2_auc")],by="cell_id",all.x=TRUE);names(maps)[names(maps)%in%c("M0_auc","M1_auc","M2_auc")]<-paste0(c("M0","M1","M2"),"_auc_dry");maps$delta_CV_R2_M2_M0<-maps$M2_cv_r2-maps$M0_cv_r2;maps$delta_RMSE_M2_M0<-maps$M2_rmse-maps$M0_rmse;maps$wet_delta_AUC_M2_M0<-maps$M2_auc_wet-maps$M0_auc_wet;maps$dry_delta_AUC_M2_M0<-maps$M2_auc_dry-maps$M0_auc_dry
write.csv(maps,file.path(OUT,"FIGURE_SOURCE_DATA_hydroclimate_maps.csv"),row.names=FALSE)
cat("Precipitation and wet/dry incremental analyses complete\n")
