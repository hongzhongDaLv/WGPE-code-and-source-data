options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
Sys.setenv(OMP_NUM_THREADS="8")
Sys.setenv(PATH=paste("__WGPE_RTOOLS__/usr/bin","__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin",Sys.getenv("PATH"),sep=";"))
ROOT<-"__WGPE_PROJECT_ROOT__";OUT<-file.path(ROOT,"output_TOP100_revision","results","incremental_value","ecology_TOP100");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(ROOT,"scripts_core_validation","core_validation_cv.cpp"),rebuild=TRUE,showOutput=FALSE)

make_target_cache<-function(target){
  cache<-file.path(OUT,paste0(target,"_target_matrix_cache.rds"));if(file.exists(cache))return(readRDS(cache))
  if(target=="GPP"){f<-file.path(ROOT,"output 1deg v3","step5_greening","step5_fixest.RData");obj<-"panel_gpp";d0<-as.Date("2000-03-01");nt<-298L}else{f<-file.path(ROOT,"output 1deg v3","step6_lai","step6_fixest.RData");obj<-"panel_lai";d0<-as.Date("2000-01-01");nt<-252L}
  e<-new.env(parent=emptyenv());load(f,envir=e);p<-e[[obj]];gr<-unique(p$grid);gi<-match(p$grid,gr);y<-matrix(NA_real_,length(gr),nt);ok<-is.finite(p$time)&p$time>=1&p$time<=nt;y[cbind(gi[ok],p$time[ok])]<-p$gpp[ok];first<-match(gr,p$grid);meta<-data.frame(grid=gr,lon=p$lon_val[first],lat=p$lat_val[first],stringsAsFactors=FALSE);rm(e,p);gc();z<-list(target=target,y=y,meta=meta,dates=seq(d0,by="month",length.out=nt));saveRDS(z,cache,compress="gzip");z
}
message("Building compact target matrices")
TG<-make_target_cache("GPP");TL<-make_target_cache("LAI")

core_file<-file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100","global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file<-file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
ctrl_file<-file.path(ROOT,"output_attribution_minimal","FULL_RERUN_FINAL_absolute_relative_v3","01_inputs_and_processing","processed_controls_FINAL_v3.rds")
e<-new.env(parent=emptyenv());load(core_file,envir=e);lon<-as.numeric(e$lon);lat<-as.numeric(e$lat);dates<-as.Date(e$all_dates);ncell<-length(lon)*length(lat)
map_target<-function(T,grid_lon,grid_lat){ii<-vapply(T$meta$lon,function(x)which.min(abs(((grid_lon-x+180)%%360)-180)),integer(1));jj<-vapply(T$meta$lat,function(x)which.min(abs(grid_lat-x)),integer(1));list(row=ii+(jj-1L)*length(grid_lon),time=match(as.character(T$dates),as.character(dates)))}
MG<-map_target(TG,lon,lat);ML<-map_target(TL,lon,lat)
slice<-function(a,M){matrix(as.numeric(a),nrow=ncell,ncol=length(dates))[M$row,M$time,drop=FALSE]}
DG<-list(y=TG$y,meta=TG$meta,dates=TG$dates);DL<-list(y=TL$y,meta=TL$meta,dates=TL$dates)
for(v in c("iwv_st","zbar_st","wgpe_st")){nm<-c(iwv_st="iwv",zbar_st="z",wgpe_st="wgpe")[[v]];message("Slicing ",nm);DG[[nm]]<-slice(e[[v]],MG);DL[[nm]]<-slice(e[[v]],ML);e[[v]]<-NULL;gc()};rm(e);gc()
s<-readRDS(single_file);MGx<-map_target(TG,as.numeric(s$lon),as.numeric(s$lat));MLx<-map_target(TL,as.numeric(s$lon),as.numeric(s$lat));for(v in c("tp_native_m","t2m_K","ssrd_native_J_m2")){nm<-c(tp_native_m="tp",t2m_K="t2m",ssrd_native_J_m2="ssrd")[[v]];message("Slicing ",nm);DG[[nm]]<-slice(s[[v]],MGx);DL[[nm]]<-slice(s[[v]],MLx);s[[v]]<-NULL;gc()};rm(s);gc()
c0<-readRDS(ctrl_file);for(v in c("rzsm_arr","vpd_arr")){nm<-c(rzsm_arr="rzsm",vpd_arr="vpd")[[v]];message("Slicing ",nm);DG[[nm]]<-slice(c0[[v]],MGx);DL[[nm]]<-slice(c0[[v]],MLx);c0[[v]]<-NULL;gc()};DG$co2<-matrix(c0$co2[MG$time],nrow=nrow(DG$y),ncol=ncol(DG$y),byrow=TRUE);DL$co2<-matrix(c0$co2[ML$time],nrow=nrow(DL$y),ncol=ncol(DL$y),byrow=TRUE);rm(c0);gc()

maskf<-file.path(ROOT,"output_attribution_minimal","FULL_RERUN_ECOLOGY_MASK_AND_PACKAGING_v5","01_ecology_domain_QA","ecology_domain_masks_v5.csv")
mask<-read.csv(maskf,check.names=FALSE);mask$key<-paste(mask$lon,mask$lat,sep="_");asbool<-function(x)tolower(as.character(x))%in%c("true","t","1")
attach_domain<-function(D,target){key<-paste(round(D$meta$lon),round(D$meta$lat),sep="_");j<-match(key,mask$key);D$meta$persistent<-asbool(mask[[paste0("persistent_vegetated_domain_",target)]][j]);D$meta$strict<-asbool(mask[[paste0("strict_stable_vegetated_domain_",target)]][j]);D$meta$biome_code<-mask$modal_IGBP_class[j];D$meta$biome_name<-mask$modal_IGBP_name[j];D}
DG<-attach_domain(DG,"GPP");DL<-attach_domain(DL,"LAI")

prepD<-function(D){td<-as.numeric(D$dates);mo<-as.integer(format(D$dates,"%m"));for(nm in c("y","iwv","z","wgpe","tp","t2m","ssrd","rzsm","co2","vpd")){message(D$target," preprocessing ",nm);D[[nm]]<-preprocess_monthly_cpp(D[[nm]],mo,td);gc()};D$interaction<-D$iwv*D$z;D}
DG$target<-"GPP";DL$target<-"LAI";DG<-prepD(DG);DL<-prepD(DL)

controls_by_model<-list(Control0=c("tp","t2m","ssrd"),Control1=c("tp","t2m","ssrd","rzsm"),Control2=c("tp","t2m","ssrd","rzsm","co2"),Control3=c("tp","t2m","ssrd","rzsm","co2","vpd"))
run_target<-function(D){fold<-as.integer((as.integer(format(D$dates,"%Y"))-min(as.integer(format(D$dates,"%Y"))))%/%5)+1L;all<-list();k<-0L;for(gap in c(0L,3L))for(cn in names(controls_by_model)){message(D$target," ",cn," gap ",gap);C<-lapply(controls_by_model[[cn]],function(v)D[[v]]);rr<-cv_linear_nested_cpp(D$y,D$iwv,D$z,D$interaction,D$wgpe,C,fold,gap,TRUE);k<-k+1L;all[[k]]<-data.frame(D$meta,target=D$target,control_set=cn,gap_months=gap,rr,check.names=FALSE)};do.call(rbind,all)}
RG<-run_target(DG);rm(DG);gc();RL<-run_target(DL);rm(DL);gc();grid<-rbind(RG,RL);write.csv(grid,file.path(OUT,"ECOLOGY_MODEL_PERFORMANCE_grid.csv"),row.names=FALSE)

band<-function(x)ifelse(x< -66.5,"S high latitude",ifelse(x< -23.5,"S temperate",ifelse(x<=23.5,"Tropics",ifelse(x<=66.5,"N temperate","N high latitude"))))
regions<-c("Global","S high latitude","S temperate","Tropics","N temperate","N high latitude")
wm<-function(x,w){ok<-is.finite(x)&is.finite(w);if(!any(ok))NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
block_ci<-function(x,w,lo,la,B=1000L){ok<-is.finite(x)&is.finite(w);if(sum(ok)<5)return(c(NA,NA));x<-x[ok];w<-w[ok];bid<-paste(floor(((lo[ok]+180)%%360)/10),floor((la[ok]+90)/10),sep="_");bn<-tapply(x*w,bid,sum);bd<-tapply(w,bid,sum);ids<-names(bn);set.seed(20260713+length(x));v<-rep(NA_real_,B);for(b in seq_len(B)){s<-sample(ids,length(ids),replace=TRUE);v[b]<-sum(bn[s])/sum(bd[s])};quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)}
summaries<-list();biomes<-list();kk<-bb<-0L
for(target in c("GPP","LAI"))for(control in names(controls_by_model))for(gap in c(0,3)){d<-grid[grid$target==target&grid$control_set==control&grid$gap_months==gap,];d$band<-band(d$lat);pairs<-list(z_increment=c(0,1),interaction_increment=c(1,2),total_vertical_position_increment=c(0,2));for(domain in c("persistent","strict")){dm<-as.logical(d[[domain]]);for(nm in names(pairs)){a<-pairs[[nm]][1];b<-pairs[[nm]][2];dv<-d[[paste0("M",b,"_cv_r2")]]-d[[paste0("M",a,"_cv_r2")]];dr<-d[[paste0("M",b,"_rmse")]]-d[[paste0("M",a,"_rmse")]];pv<-if(nm=="z_increment")d$increment_z_p else if(nm=="interaction_increment")d$increment_interaction_p else d$increment_total_p;qv<-rep(NA_real_,length(pv));idq<-dm&is.finite(pv);qv[idq]<-p.adjust(pv[idq],"BH");pr<-if(nm=="z_increment")d$partial_r_z_given_I_controls else if(nm=="interaction_increment")d$partial_r_int_given_I_z_controls else rep(NA_real_,nrow(d));pp<-if(nm=="z_increment")d$partial_p_z else if(nm=="interaction_increment")d$partial_p_interaction else rep(NA_real_,nrow(d));pq<-rep(NA_real_,length(pp));ip<-dm&is.finite(pp);pq[ip]<-p.adjust(pp[ip],"BH")
      for(r in regions){ir<-dm&(if(r=="Global")TRUE else d$band==r)&is.finite(dv);w<-cos(d$lat[ir]*pi/180);ci<-block_ci(dv[ir],w,d$lon[ir],d$lat[ir]);den<-sum(w);nmon<-if(any(ir))max(d$n_complete[ir],na.rm=TRUE)else NA_real_;kk<-kk+1L;summaries[[kk]]<-data.frame(target,control_set=control,domain=if(domain=="persistent")"persistent_vegetated_domain" else "strict_stable_vegetated_domain",gap_months=gap,increment=nm,region=r,delta_CV_R2=wm(dv[ir],w),lower95=ci[1],upper95=ci[2],delta_RMSE=wm(dr[ir],w),positive_increment_area_pct=if(den>0)100*sum(w[dv[ir]>0])/den else NA_real_,significant_positive_increment_area_pct=if(den>0)100*sum(w[is.finite(qv[ir])&qv[ir]<=.05&dv[ir]>0])/den else NA_real_,mean_partial_r=wm(pr[ir],w),partial_r_positive_significant_area_pct=if(den>0)100*sum(w[is.finite(pq[ir])&pq[ir]<=.05&pr[ir]>0])/den else NA_real_,n_grid=sum(ir),n_month=nmon,cv_scheme="continuous 5-year folds")}
      }
      if(gap==3&&domain=="persistent"){for(code in sort(unique(d$biome_code[dm]))){ib<-dm&d$biome_code==code&is.finite(dv);if(sum(ib)<100)next;w<-cos(d$lat[ib]*pi/180);ci<-block_ci(dv[ib],w,d$lon[ib],d$lat[ib]);bb<-bb+1L;biomes[[bb]]<-data.frame(target,control_set=control,increment=nm,biome_code=code,biome_name=d$biome_name[which(ib)[1]],delta_CV_R2=wm(dv[ib],w),lower95=ci[1],upper95=ci[2],positive_increment_area_pct=100*sum(w[dv[ib]>0])/sum(w),n_grid=sum(ib),status="eligible")}}
    }}
sm<-do.call(rbind,summaries);bio<-do.call(rbind,biomes)
write.csv(sm[sm$target=="GPP",],file.path(OUT,"GPP_INCREMENT_BY_CONTROL.csv"),row.names=FALSE)
write.csv(sm[sm$target=="LAI",],file.path(OUT,"LAI_INCREMENT_BY_CONTROL.csv"),row.names=FALSE)
write.csv(bio,file.path(OUT,"BIOME_INCREMENT.csv"),row.names=FALSE)
coll<-grid[,c("grid","lon","lat","target","control_set","gap_months","predictor_cor_I_Z","vif_I_Z","condition_number","rank","coef_sd_z","coef_sd_interaction","fit_status")];write.csv(coll,file.path(OUT,"COLLINEARITY_DIAGNOSTICS_ecology.csv"),row.names=FALSE)
cv<-aggregate(cbind(n_complete,n_oof,fit_status)~target+control_set+gap_months,grid,function(x)c(n=length(x),mean=mean(x,na.rm=TRUE),warning=sum(x==2,na.rm=TRUE)));write.csv(cv,file.path(OUT,"CV_DIAGNOSTICS_ecology.csv"),row.names=FALSE)
maps<-grid[grid$gap_months==3&grid$control_set%in%names(controls_by_model),c("grid","lon","lat","target","control_set","persistent","strict","M0_cv_r2","M1_cv_r2","M2_cv_r2","M0_rmse","M1_rmse","M2_rmse")];maps$delta_CV_R2_M1_M0<-maps$M1_cv_r2-maps$M0_cv_r2;maps$delta_CV_R2_M2_M1<-maps$M2_cv_r2-maps$M1_cv_r2;maps$delta_CV_R2_M2_M0<-maps$M2_cv_r2-maps$M0_cv_r2;write.csv(maps,file.path(OUT,"FIGURE_SOURCE_DATA_ecology_maps.csv"),row.names=FALSE)
cat("Ecology incremental analysis complete\n")
