# Purpose: Correct the lag-filtered spatial summaries from the completed
#          ecology temporal-effect grid cache without rerunning grid models.
# Inputs: ECOLOGY_TEMPORAL_EFFECT_SIZE_GRID_TOP100.csv.
# Outputs: corrected ECOLOGY_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv.
# Parameters: persistent vegetated domain; 10-degree block bootstrap, B=1000.
# Dependencies: base R.
# Overwrite policy: replaces only the erroneous TOP100 summary; the grid cache
#                   is read-only and all legacy files remain untouched.

options(stringsAsFactors=FALSE,warn=1)
ROOT<-"__WGPE_PROJECT_ROOT__";OUT<-file.path(ROOT,"output_TOP100_revision","results","temporal_precedence","TOP100_effect_sizes")
grid<-read.csv(file.path(OUT,"ECOLOGY_TEMPORAL_EFFECT_SIZE_GRID_TOP100.csv"),check.names=FALSE)
band<-function(x)ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))));regions<-c("Global","S high latitude","S temperate","Tropics","N temperate","N high latitude")
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
block_ci<-function(x,w,lo,la,B=1000L){ok<-is.finite(x)&is.finite(w);if(sum(ok)<5)return(c(NA,NA));x<-x[ok];w<-w[ok];bid<-paste(floor(((lo[ok]+180)%%360)/10),floor((la[ok]+90)/10),sep="_");bn<-tapply(x*w,bid,sum);bd<-tapply(w,bid,sum);ids<-names(bn);set.seed(20260714+length(x));v<-rep(NA_real_,B);for(b in seq_len(B)){ss<-sample(ids,length(ids),TRUE);v[b]<-sum(bn[ss])/sum(bd[ss])};quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)}
grid$band<-band(grid$lat);metrics<-c("partial_R2","delta_cv_R2","delta_RMSE_improvement");rows<-list();z<-0L
for(target in c("GPP","LAI"))for(lag in 1:4)for(reg in regions)for(metric in metrics){d<-grid[grid$target==target&grid$lag==lag&grid$persistent&(reg=="Global"|grid$band==reg),];v<-d[[metric]];w<-cos(d$lat*pi/180);ci<-block_ci(v,w,d$lon,d$lat);z<-z+1L;rows[[z]]<-data.frame(target=target,predictor="W-GPE",lag=lag,region=reg,metric=metric,estimate=wm(v,w),lower95=ci[1],upper95=ci[2],positive_area_pct=100*sum(w[is.finite(v)&v>0])/sum(w[is.finite(v)]),n_grid=sum(is.finite(v)),controls="precipitation+T2m+SSRD+RZSM+CO2+VPD",domain="persistent_vegetated",uncertainty_method="10-degree spatial-block bootstrap, 1000 resamples")}
write.csv(do.call(rbind,rows),file.path(OUT,"ECOLOGY_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv"),row.names=FALSE)
cat("Corrected lag-filtered ecology temporal summaries\n")
