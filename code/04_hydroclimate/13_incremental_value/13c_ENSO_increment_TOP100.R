# Purpose: Quantify TOP100 vertical-position incremental value for Niño3.4,
#          ENSO-winter precipitation, and El Niño/La Niña phase classification.
# Inputs: formal TOP100 monthly core, ERA5 precipitation, NOAA Niño3.4 and the
#         locked ONI-standard winter inventory.
# Outputs: grid-level model performance and climate-zone increment summaries.
# Parameters: M0=IWV; M1=IWV+<z>; M2=IWV+<z>+IWV:<z>; M3=W-GPE;
#             monthly deseasoning, linear detrending, within-grid scaling;
#             shared folds, gaps 0 and 3 months, 10-degree spatial bootstrap.
# Dependencies: Rcpp, RcppArmadillo, Rtools45, core_validation_cv.cpp.
# Overwrite policy: writes only TOP100-specific ENSO incremental outputs.

options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
Sys.setenv(OMP_NUM_THREADS="8")
ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- Sys.getenv("TOP100_CLOSURE_ENSO_INCREMENT_OUT",file.path(ROOT,"output_TOP100_revision","results","incremental_value","ENSO_TOP100"))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
source(file.path(ROOT,"scripts_TOP100_revision","18_coordinate_harmonization","coordinate_harmonization.R"))
Rcpp::sourceCpp(file.path(ROOT,"scripts_core_validation","core_validation_cv.cpp"),rebuild=TRUE,showOutput=FALSE)

core_file <- file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100",
  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file <- file.path(ROOT,"output_attribution_minimal",
  "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
nino_file <- file.path(ROOT,"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard",
  "11_ENSO_ONI_STANDARD_UPDATE","00_external_data","NOAA_PSL_Nino34_monthly_clean.csv")
event_file <- file.path(ROOT,"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard",
  "11_ENSO_ONI_STANDARD_UPDATE","01_event_definition","ONI_winter_realizations_1979_2024.csv")

e<-new.env(parent=emptyenv());load(core_file,envir=e)
lon<-as.numeric(e$lon);lat<-as.numeric(e$lat);dates<-as.Date(e$all_dates)
nlon<-length(lon);nlat<-length(lat);nt<-length(dates);ncell<-nlon*nlat
cell<-expand.grid(lon=lon,lat=lat);cell$cell_id<-seq_len(ncell);cell$weight<-cos(cell$lat*pi/180)
to_mat<-function(a)matrix(as.numeric(a),nrow=ncell,ncol=nt)
m_i<-to_mat(e$iwv_st);m_z<-to_mat(e$zbar_st);m_w<-to_mat(e$wgpe_st);rm(e);gc()
s<-readRDS(single_file);m_p<-to_mat(reorder_exact_grid(s$tp_native_m,s$lon,s$lat,lon,lat));rm(s);gc()

mon<-as.integer(format(dates,"%m"));yr<-as.integer(format(dates,"%Y"));tdec<-yr+(mon-.5)/12
prep<-function(m){
  good<-rowSums(is.finite(m))==ncol(m);m[!good,]<-NA_real_
  for(mm in 1:12){id<-which(mon==mm);cl<-rowMeans(m[,id,drop=FALSE],na.rm=TRUE);m[,id]<-sweep(m[,id,drop=FALSE],1,cl,"-")}
  tc<-tdec-mean(tdec);mu<-rowMeans(m,na.rm=TRUE);yc<-sweep(m,1,mu,"-");sl<-as.numeric(yc%*%tc/sum(tc^2));m<-m-outer(sl,tc)
  mu<-rowMeans(m,na.rm=TRUE);m<-sweep(m,1,mu,"-");sdv<-sqrt(rowSums(m^2,na.rm=TRUE)/pmax(1,rowSums(is.finite(m))-1));ok<-is.finite(sdv)&sdv>0
  m[ok,]<-m[ok,,drop=FALSE]/sdv[ok];m[!ok,]<-NA_real_;m
}
message("Preprocessing TOP100 ENSO incremental matrices")
m_i<-prep(m_i);m_z<-prep(m_z);m_w<-prep(m_w);m_p<-prep(m_p);m_int<-m_i*m_z;gc()

nino<-read.csv(nino_file);stopifnot(identical(as.character(as.Date(nino$date)),as.character(dates)))
y_nino<-as.numeric(nino$nino34_standardized)
real<-read.csv(event_file);real<-real[tolower(as.character(real$included_main))%in%c("true","t","1"),]
real$winter_start<-as.Date(real$winter_start);real$winter_end<-as.Date(real$winter_end)
idx_by_winter<-lapply(seq_len(nrow(real)),function(k)which(dates>=real$winter_start[k]&dates<=real$winter_end[k]))
stopifnot(all(lengths(idx_by_winter)==6L))
idx_winter<-unlist(idx_by_winter)
phase_vec<-unlist(lapply(seq_len(nrow(real)),function(k)rep(ifelse(real$event_type[k]=="ElNino",1,0),6)))
fold_episode<-unlist(lapply(seq_len(nrow(real)),function(k)rep(match(real$episode_id[k],unique(real$episode_id)),6)))
fold_full<-as.integer((yr-min(yr))%/%5)+1L
gaps<-c(0L,3L)

band<-function(x)ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))))
cell$climate_band<-band(cell$lat);regions<-c("Global","N high latitude","N temperate","Tropics","S temperate","S high latitude")
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
block_ci<-function(x,w,lo,la,B=1000L){ok<-is.finite(x)&is.finite(w);if(sum(ok)<5)return(c(NA,NA));x<-x[ok];w<-w[ok];bid<-paste(floor(((lo[ok]+180)%%360)/10),floor((la[ok]+90)/10),sep="_");bn<-tapply(x*w,bid,sum);bd<-tapply(w,bid,sum);ids<-names(bn);set.seed(20260714+length(x));v<-rep(NA_real_,B);for(b in seq_len(B)){ss<-sample(ids,length(ids),replace=TRUE);v[b]<-sum(bn[ss])/sum(bd[ss])};quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)}

run_linear<-function(target_name,Y,I,Z,INT,W,fold,batch=2500L){out<-list();k<-0L;for(gap in gaps)for(st in seq.int(1L,ncell,by=batch)){en<-min(ncell,st+batch-1L);id<-st:en;k<-k+1L;message(target_name," linear gap ",gap," cells ",st,"-",en);YY<-if(is.null(dim(Y)))matrix(Y,nrow=length(id),ncol=length(Y),byrow=TRUE)else Y[id,,drop=FALSE];rr<-cv_linear_nested_cpp(YY,I[id,,drop=FALSE],Z[id,,drop=FALSE],INT[id,,drop=FALSE],W[id,,drop=FALSE],list(),fold,gap,FALSE);out[[k]]<-data.frame(cell[id,],target=target_name,gap_months=gap,rr,check.names=FALSE)};do.call(rbind,out)}
run_logistic<-function(target_name,Y,I,Z,INT,fold,batch=1000L){out<-list();k<-0L;for(gap in gaps)for(st in seq.int(1L,ncell,by=batch)){en<-min(ncell,st+batch-1L);id<-st:en;k<-k+1L;message(target_name," logistic gap ",gap," cells ",st,"-",en);YY<-matrix(Y,nrow=length(id),ncol=length(Y),byrow=TRUE);rr<-cv_logistic_nested_cpp(YY,I[id,,drop=FALSE],Z[id,,drop=FALSE],INT[id,,drop=FALSE],fold,gap);out[[k]]<-data.frame(cell[id,],target=target_name,gap_months=gap,rr,check.names=FALSE)};do.call(rbind,out)}

f_n<-file.path(OUT,"NINO34_MODEL_PERFORMANCE_GRID_TOP100.csv")
f_w<-file.path(OUT,"ENSO_WINTER_PRECIP_MODEL_PERFORMANCE_GRID_TOP100.csv")
f_p<-file.path(OUT,"ENSO_PHASE_MODEL_PERFORMANCE_GRID_TOP100.csv")
if(file.exists(f_n))rn<-read.csv(f_n,check.names=FALSE)else{rn<-run_linear("Nino3.4",y_nino,m_i,m_z,m_int,m_w,fold_full);write.csv(rn,f_n,row.names=FALSE)}
if(file.exists(f_w))rw<-read.csv(f_w,check.names=FALSE)else{rw<-run_linear("ENSO_winter_precipitation",m_p[,idx_winter,drop=FALSE],m_i[,idx_winter,drop=FALSE],m_z[,idx_winter,drop=FALSE],m_int[,idx_winter,drop=FALSE],m_w[,idx_winter,drop=FALSE],fold_episode);write.csv(rw,f_w,row.names=FALSE)}
if(file.exists(f_p))rp<-read.csv(f_p,check.names=FALSE)else{rp<-run_logistic("ENSO_phase",phase_vec,m_i[,idx_winter,drop=FALSE],m_z[,idx_winter,drop=FALSE],m_int[,idx_winter,drop=FALSE],fold_episode);write.csv(rp,f_p,row.names=FALSE)}

summ_linear<-function(d){o<-list();k<-0L;pairs<-list(z_increment=c(0,1),interaction_increment=c(1,2),total_vertical_position_increment=c(0,2));for(gap in gaps){x<-d[d$gap_months==gap,];x$climate_band<-band(x$lat);for(nm in names(pairs)){a<-pairs[[nm]][1];b<-pairs[[nm]][2];dv<-x[[paste0("M",b,"_cv_r2")]]-x[[paste0("M",a,"_cv_r2")]];dr<-x[[paste0("M",b,"_rmse")]]-x[[paste0("M",a,"_rmse")]];for(reg in regions){ii<-if(reg=="Global")rep(TRUE,nrow(x))else x$climate_band==reg;ww<-x$weight[ii];ci<-block_ci(dv[ii],ww,x$lon[ii],x$lat[ii]);den<-sum(ww[is.finite(dv[ii])]);k<-k+1L;o[[k]]<-data.frame(target=x$target[1],gap_months=gap,increment=nm,region=reg,delta_CV_R2=wm(dv[ii],ww),lower95=ci[1],upper95=ci[2],delta_RMSE=wm(dr[ii],ww),positive_increment_area_pct=100*sum(ww[is.finite(dv[ii])&dv[ii]>0])/den,n_grid=sum(is.finite(dv[ii])))}}};do.call(rbind,o)}
summ_logistic<-function(d){o<-list();k<-0L;pairs<-list(z_increment=c(0,1),interaction_increment=c(1,2),total_vertical_position_increment=c(0,2));for(gap in gaps){x<-d[d$gap_months==gap,];x$climate_band<-band(x$lat);for(nm in names(pairs)){a<-pairs[[nm]][1];b<-pairs[[nm]][2];da<-x[[paste0("M",b,"_auc")]]-x[[paste0("M",a,"_auc")]];for(reg in regions){ii<-if(reg=="Global")rep(TRUE,nrow(x))else x$climate_band==reg;ww<-x$weight[ii];ci<-block_ci(da[ii],ww,x$lon[ii],x$lat[ii]);den<-sum(ww[is.finite(da[ii])]);k<-k+1L;o[[k]]<-data.frame(target=x$target[1],gap_months=gap,increment=nm,region=reg,delta_AUC=wm(da[ii],ww),lower95=ci[1],upper95=ci[2],positive_increment_area_pct=100*sum(ww[is.finite(da[ii])&da[ii]>0])/den,n_grid=sum(is.finite(da[ii])))}}};do.call(rbind,o)}

write.csv(summ_linear(rn),file.path(OUT,"NINO34_INCREMENT_SUMMARY_TOP100.csv"),row.names=FALSE)
write.csv(summ_linear(rw),file.path(OUT,"ENSO_WINTER_PRECIP_INCREMENT_SUMMARY_TOP100.csv"),row.names=FALSE)
write.csv(summ_logistic(rp),file.path(OUT,"ENSO_PHASE_INCREMENT_SUMMARY_TOP100.csv"),row.names=FALSE)
write.csv(data.frame(target=c("Nino3.4","ENSO_winter_precipitation","ENSO_phase"),period=c("1979-2024 monthly","locked ONI-standard Oct-Mar winters","locked ONI-standard Oct-Mar winters"),fold_definition=c("continuous 5-year folds","episode-cluster folds","episode-cluster folds"),n_month=c(nt,length(idx_winter),length(idx_winter)),n_winter=nrow(real),n_episode=length(unique(real$episode_id))),file.path(OUT,"ENSO_INCREMENT_METHOD_QA_TOP100.csv"),row.names=FALSE)
cat("ENSO incremental analyses complete\n")
