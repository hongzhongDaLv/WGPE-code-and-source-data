# Purpose: Add effect sizes to TOP100 W-GPE-leading-precipitation lag tests.
# Inputs: formal TOP100 monthly W-GPE and ERA5 precipitation.
# Outputs: lag 1--4 grid maps and climate-zone partial-R2/blocked-CV summaries.
# Parameters: deseasoned and detrended standardized series; 5-year blocked folds;
#             3-month training gap; 10-degree spatial-block bootstrap.
# Dependencies: Rcpp, RcppArmadillo, Rtools45.
# Overwrite policy: TOP100-specific outputs only.

options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
ROOT<-"__WGPE_PROJECT_ROOT__"
OUT<-Sys.getenv("TOP100_CLOSURE_PRECIP_TEMP_EFFECT_OUT",file.path(ROOT,"output_TOP100_revision","results","temporal_precedence","TOP100_effect_sizes"))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
source(file.path(ROOT,"scripts_TOP100_revision","18_coordinate_harmonization","coordinate_harmonization.R"))
Rcpp::sourceCpp(file.path(ROOT,"scripts_core_validation","core_validation_cv.cpp"),rebuild=TRUE,showOutput=FALSE)
Rcpp::sourceCpp(file.path(ROOT,"scripts_TOP100_revision","08_temporal_precedence","temporal_precedence_effect_cv.cpp"),rebuild=TRUE,showOutput=FALSE)

e<-new.env(parent=emptyenv());load(file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100","global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData"),envir=e)
lon<-as.numeric(e$lon);lat<-as.numeric(e$lat);dates<-as.Date(e$all_dates);nt<-length(dates);ncell<-length(lon)*length(lat)
cell<-expand.grid(lon=lon,lat=lat);cell$cell_id<-seq_len(ncell);cell$weight<-cos(cell$lat*pi/180)
W<-matrix(as.numeric(e$wgpe_st),nrow=ncell,ncol=nt);rm(e);gc()
s<-readRDS(file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds"))
P<-matrix(as.numeric(reorder_exact_grid(s$tp_native_m,s$lon,s$lat,lon,lat)),nrow=ncell,ncol=nt);rm(s);gc()
month<-as.integer(format(dates,"%m"));td<-as.numeric(dates)
message("Preprocessing precipitation temporal-effect matrices")
W<-preprocess_monthly_cpp(W,month,td);P<-preprocess_monthly_cpp(P,month,td);gc()
fold<-as.integer((as.integer(format(dates,"%Y"))-1979)%/%5)+1L

all<-list();k<-0L;batch<-2500L
for(lag in 1:4)for(st in seq.int(1L,ncell,by=batch)){en<-min(ncell,st+batch-1L);id<-st:en;k<-k+1L;message("lag ",lag," cells ",st,"-",en);rr<-temporal_effect_cv_cpp(P[id,,drop=FALSE],W[id,,drop=FALSE],list(),fold,lag,3L);all[[k]]<-data.frame(cell[id,],target="precipitation",predictor="W-GPE",lag=lag,rr,check.names=FALSE)}
grid<-do.call(rbind,all);write.csv(grid,file.path(OUT,"PRECIP_TEMPORAL_EFFECT_SIZE_GRID_TOP100.csv"),row.names=FALSE)

band<-function(x)ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))))
regions<-c("Global","N high latitude","N temperate","Tropics","S temperate","S high latitude")
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
block_ci<-function(x,w,lo,la,B=1000L){ok<-is.finite(x)&is.finite(w);if(sum(ok)<5)return(c(NA,NA));x<-x[ok];w<-w[ok];bid<-paste(floor(((lo[ok]+180)%%360)/10),floor((la[ok]+90)/10),sep="_");bn<-tapply(x*w,bid,sum);bd<-tapply(w,bid,sum);ids<-names(bn);set.seed(20260714+length(x));v<-rep(NA_real_,B);for(b in seq_len(B)){ss<-sample(ids,length(ids),TRUE);v[b]<-sum(bn[ss])/sum(bd[ss])};quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)}
grid$band<-band(grid$lat);metrics<-c("partial_R2","delta_cv_R2","delta_RMSE_improvement");rows<-list();z<-0L
for(lag in 1:4)for(reg in regions)for(metric in metrics){d<-grid[grid$lag==lag&(reg=="Global"|grid$band==reg),];v<-d[[metric]];ci<-block_ci(v,d$weight,d$lon,d$lat);z<-z+1L;rows[[z]]<-data.frame(target="precipitation",predictor="W-GPE",lag=lag,region=reg,metric=metric,estimate=wm(v,d$weight),lower95=ci[1],upper95=ci[2],positive_area_pct=100*sum(d$weight[is.finite(v)&v>0])/sum(d$weight[is.finite(v)]),n_grid=sum(is.finite(v)),uncertainty_method="10-degree spatial-block bootstrap, 1000 resamples")}
write.csv(do.call(rbind,rows),file.path(OUT,"PRECIP_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv"),row.names=FALSE)
cat("Precipitation temporal effect sizes complete\n")
