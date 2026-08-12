# Purpose: Add conditional effect sizes to TOP100 W-GPE-leading-GPP/LAI tests.
# Inputs: formal TOP100 W-GPE, GPP/LAI panels, precipitation, T2m, SSRD,
#         root-zone soil moisture, CO2, VPD and persistent-vegetation masks.
# Outputs: lag 1--4 partial-R2 and blocked-CV delta-R2/RMSE maps and summaries.
# Parameters: Model-3 controls; deseasoned/detrended standardized series;
#             continuous 5-year folds with 3-month gap.
# Dependencies: Rcpp, RcppArmadillo, Rtools45.
# Overwrite policy: TOP100-specific outputs only.

options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
ROOT<-"__WGPE_PROJECT_ROOT__";OUT<-file.path(ROOT,"output_TOP100_revision","results","temporal_precedence","TOP100_effect_sizes");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(ROOT,"scripts_core_validation","core_validation_cv.cpp"),rebuild=TRUE,showOutput=FALSE)
Rcpp::sourceCpp(file.path(ROOT,"scripts_TOP100_revision","08_temporal_precedence","temporal_precedence_effect_cv.cpp"),rebuild=TRUE,showOutput=FALSE)

make_target<-function(target){cache<-file.path(ROOT,"output_TOP100_revision","results","incremental_value","ecology_TOP100",paste0(target,"_target_matrix_cache.rds"));if(file.exists(cache))return(readRDS(cache));stop("Missing audited target cache: ",cache)}
TG<-make_target("GPP");TL<-make_target("LAI")
core<-new.env(parent=emptyenv());load(file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100","global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData"),envir=core)
lon<-as.numeric(core$lon);lat<-as.numeric(core$lat);dates<-as.Date(core$all_dates);ncell<-length(lon)*length(lat)
map_target<-function(T){ii<-vapply(round(T$meta$lon),function(x)which.min(abs(lon-x)),integer(1));jj<-vapply(round(T$meta$lat),function(x)which.min(abs(lat-x)),integer(1));list(row=ii+(jj-1L)*length(lon),time=match(as.character(T$dates),as.character(dates)))}
MG<-map_target(TG);ML<-map_target(TL);slice<-function(a,M)matrix(as.numeric(a),nrow=ncell,ncol=length(dates))[M$row,M$time,drop=FALSE]
DG<-list(target="GPP",y=TG$y,meta=TG$meta,dates=TG$dates,wgpe=slice(core$wgpe_st,MG));DL<-list(target="LAI",y=TL$y,meta=TL$meta,dates=TL$dates,wgpe=slice(core$wgpe_st,ML));rm(core,TG,TL);gc()
s<-readRDS(file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds"));for(v in c("tp_native_m","t2m_K","ssrd_native_J_m2")){nm<-c(tp_native_m="tp",t2m_K="t2m",ssrd_native_J_m2="ssrd")[[v]];DG[[nm]]<-slice(s[[v]],MG);DL[[nm]]<-slice(s[[v]],ML);s[[v]]<-NULL;gc()};rm(s);gc()
c0<-readRDS(file.path(ROOT,"output_attribution_minimal","FULL_RERUN_FINAL_absolute_relative_v3","01_inputs_and_processing","processed_controls_FINAL_v3.rds"));for(v in c("rzsm_arr","vpd_arr")){nm<-c(rzsm_arr="rzsm",vpd_arr="vpd")[[v]];DG[[nm]]<-slice(c0[[v]],MG);DL[[nm]]<-slice(c0[[v]],ML);c0[[v]]<-NULL;gc()};DG$co2<-matrix(c0$co2[MG$time],nrow=nrow(DG$y),ncol=ncol(DG$y),byrow=TRUE);DL$co2<-matrix(c0$co2[ML$time],nrow=nrow(DL$y),ncol=ncol(DL$y),byrow=TRUE);rm(c0);gc()

mask<-read.csv(file.path(ROOT,"output_attribution_minimal","FULL_RERUN_ECOLOGY_MASK_AND_PACKAGING_v5","01_ecology_domain_QA","ecology_domain_masks_v5.csv"),check.names=FALSE);mask$key<-paste(mask$lon,mask$lat,sep="_");asbool<-function(x)tolower(as.character(x))%in%c("true","t","1")
attach_domain<-function(D){key<-paste(round(D$meta$lon),round(D$meta$lat),sep="_");j<-match(key,mask$key);D$meta$persistent<-asbool(mask[[paste0("persistent_vegetated_domain_",D$target)]][j]);D}
DG<-attach_domain(DG);DL<-attach_domain(DL);rm(mask);gc()
prepD<-function(D){mo<-as.integer(format(D$dates,"%m"));td<-as.numeric(D$dates);for(nm in c("y","wgpe","tp","t2m","ssrd","rzsm","co2","vpd")){message(D$target," preprocessing ",nm);D[[nm]]<-preprocess_monthly_cpp(D[[nm]],mo,td);gc()};D}
DG<-prepD(DG);DL<-prepD(DL)

run_target<-function(D){fold<-as.integer((as.integer(format(D$dates,"%Y"))-min(as.integer(format(D$dates,"%Y"))))%/%5)+1L;C<-lapply(c("tp","t2m","ssrd","rzsm","co2","vpd"),function(v)D[[v]]);o<-list();for(lag in 1:4){message(D$target," effect lag ",lag);rr<-temporal_effect_cv_cpp(D$y,D$wgpe,C,fold,lag,3L);o[[lag]]<-data.frame(D$meta,target=D$target,predictor="W-GPE",lag=lag,rr,check.names=FALSE)};do.call(rbind,o)}
RG<-run_target(DG);rm(DG);gc();RL<-run_target(DL);rm(DL);gc();grid<-rbind(RG,RL);write.csv(grid,file.path(OUT,"ECOLOGY_TEMPORAL_EFFECT_SIZE_GRID_TOP100.csv"),row.names=FALSE)

band<-function(x)ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))));regions<-c("Global","S high latitude","S temperate","Tropics","N temperate","N high latitude")
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
block_ci<-function(x,w,lo,la,B=1000L){ok<-is.finite(x)&is.finite(w);if(sum(ok)<5)return(c(NA,NA));x<-x[ok];w<-w[ok];bid<-paste(floor(((lo[ok]+180)%%360)/10),floor((la[ok]+90)/10),sep="_");bn<-tapply(x*w,bid,sum);bd<-tapply(w,bid,sum);ids<-names(bn);set.seed(20260714+length(x));v<-rep(NA_real_,B);for(b in seq_len(B)){ss<-sample(ids,length(ids),TRUE);v[b]<-sum(bn[ss])/sum(bd[ss])};quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)}
grid$band<-band(grid$lat);metrics<-c("partial_R2","delta_cv_R2","delta_RMSE_improvement");rows<-list();z<-0L
for(target in c("GPP","LAI"))for(lag in 1:4)for(reg in regions)for(metric in metrics){d<-grid[grid$target==target&grid$lag==lag&grid$persistent&(reg=="Global"|grid$band==reg),];v<-d[[metric]];w<-cos(d$lat*pi/180);ci<-block_ci(v,w,d$lon,d$lat);z<-z+1L;rows[[z]]<-data.frame(target=target,predictor="W-GPE",lag=lag,region=reg,metric=metric,estimate=wm(v,w),lower95=ci[1],upper95=ci[2],positive_area_pct=100*sum(w[is.finite(v)&v>0])/sum(w[is.finite(v)]),n_grid=sum(is.finite(v)),controls="precipitation+T2m+SSRD+RZSM+CO2+VPD",domain="persistent_vegetated",uncertainty_method="10-degree spatial-block bootstrap, 1000 resamples")}
write.csv(do.call(rbind,rows),file.path(OUT,"ECOLOGY_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv"),row.names=FALSE)
cat("Ecology temporal effect sizes complete\n")
