options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
ROOT<-"__WGPE_PROJECT_ROOT__"
CANONICAL<-file.path(ROOT,"CANONICAL_WGPE_TOP100_v2")
OUT<-file.path(CANONICAL,"08_REPRODUCTION_OUTPUT","Figure3f")
V3<-file.path(ROOT,"output_attribution_minimal","FULL_RERUN_FINAL_absolute_relative_v3")
core_file<-file.path(ROOT,"output_TOP100_revision","derived_data","ERA5_TOP100","global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file<-file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
ctrl_file<-file.path(V3,"01_inputs_and_processing","processed_controls_FINAL_v3.rds")
gpp_file<-file.path(ROOT,"output 1deg v3","step5_greening","step5_fixest.RData")
lai_file<-file.path(ROOT,"output 1deg v3","step6_lai","step6_fixest.RData")

e<-new.env(parent=emptyenv());load(core_file,envir=e);lon<-as.numeric(e$lon);lat<-as.numeric(e$lat);dates<-as.Date(e$all_dates);wgpe<-e$wgpe_st;rm(e);gc()
single<-readRDS(single_file);ii_lon<-match((lon+360)%%360,(as.numeric(single$lon)+360)%%360);ii_lat<-match(lat,as.numeric(single$lat))
tp<-single$tp_native_m[ii_lon,ii_lat,,drop=FALSE];t2m<-single$t2m_K[ii_lon,ii_lat,,drop=FALSE];ssrd<-single$ssrd_native_J_m2[ii_lon,ii_lat,,drop=FALSE];rm(single);gc()
ctrl<-readRDS(ctrl_file);rzsm<-ctrl$rzsm_arr[ii_lon,ii_lat,,drop=FALSE];vpd<-ctrl$vpd_arr[ii_lon,ii_lat,,drop=FALSE];co2<-as.numeric(ctrl$co2);rm(ctrl);gc();nlon<-length(lon);nlat<-length(lat)

detrend_group<-function(v){ok<-is.finite(v);out<-rep(NA_real_,length(v));if(sum(ok)<10)return(out);tt<-seq_along(v);out[ok]<-residuals(lm(v[ok]~tt[ok]));out}
augment<-function(panel,panel_dates){panel$date<-panel_dates[panel$time];panel$month<-as.integer(format(panel$date,"%m"));ii<-vapply(panel$lon_val,function(v)which.min(abs(((lon-v+180)%%360)-180)),integer(1));jj<-vapply(panel$lat_val,function(v)which.min(abs(lat-v)),integer(1));kk<-match(as.character(panel$date),as.character(dates));good<-!is.na(ii)&!is.na(jj)&!is.na(kk);panel<-panel[good,];ii<-ii[good];jj<-jj[good];kk<-kk[good];lin<-as.numeric(ii)+(as.numeric(jj)-1)*nlon+(as.numeric(kk)-1)*nlon*nlat;panel$wgpe<-as.numeric(wgpe)[lin];panel$tp<-as.numeric(tp)[lin];panel$t2m<-as.numeric(t2m)[lin];panel$ssrd<-as.numeric(ssrd)[lin];panel$rzsm<-as.numeric(rzsm)[lin];panel$vpd<-as.numeric(vpd)[lin];panel$co2<-co2[kk];for(vn in c("gpp","wgpe","tp","t2m","ssrd","rzsm","vpd","co2")){panel[[vn]]<-panel[[vn]]-ave(panel[[vn]],panel$grid,panel$month,FUN=function(x)mean(x,na.rm=TRUE));panel[[vn]]<-ave(panel[[vn]],panel$grid,FUN=detrend_group)};panel}
lagmat<-function(x,p){n<-length(x);o<-matrix(NA_real_,n-p,p);for(k in 1:p)o[,k]<-x[(p+1-k):(n-k)];o}
hac_vcov<-function(X,u,L){n<-nrow(X);k<-ncol(X);S<-X*u;meat<-crossprod(S);for(l in 1:min(L,n-1)){w<-1-l/(L+1);G<-crossprod(S[(l+1):n,,drop=FALSE],S[1:(n-l),,drop=FALSE]);meat<-meat+w*(G+t(G))};bread<-tryCatch(solve(crossprod(X)),error=function(e)MASS::ginv(crossprod(X)));bread%*%meat%*%bread*n/(n-k)}
fit_one<-function(target,added,controls,p=3){y<-target[(p+1):length(target)];Xr<-cbind(1,lagmat(target,p));for(cc in controls)Xr<-cbind(Xr,lagmat(cc,p));Xa<-lagmat(added,p);Xf<-cbind(Xr,Xa);ok<-is.finite(y)&apply(is.finite(Xf),1,all);y<-y[ok];Xr<-Xr[ok,,drop=FALSE];Xf<-Xf[ok,,drop=FALSE];if(length(y)<ncol(Xf)+8)return(NULL);nr<-ncol(Xr);for(j in 2:ncol(Xf)){mu<-mean(Xf[,j]);sdv<-sd(Xf[,j]);if(!is.finite(sdv)||sdv<=0)return(NULL);Xf[,j]<-(Xf[,j]-mu)/sdv};Xr<-Xf[,seq_len(nr),drop=FALSE];fr<-lm.fit(Xr,y);ff<-lm.fit(Xf,y);if(fr$rank<ncol(Xr)||ff$rank<ncol(Xf))return(NULL);rssr<-sum(fr$residuals^2);rssf<-sum(ff$residuals^2);df1<-p;df2<-length(y)-ncol(Xf);Fv<-max(0,((rssr-rssf)/df1)/(rssf/df2));idx<-(ncol(Xf)-p+1):ncol(Xf);wald<-function(L){V<-hac_vcov(Xf,ff$residuals,L)[idx,idx,drop=FALSE];b<-ff$coefficients[idx];q<-as.numeric(t(b)%*%tryCatch(solve(V),error=function(e)MASS::ginv(V))%*%b);c(stat=q,p=pchisq(q,df=p,lower.tail=FALSE))};w3<-wald(3);w12<-wald(12);u<-ff$residuals;c(p_F=unname(pf(Fv,df1,df2,lower.tail=FALSE)),p_HAC3=unname(w3["p"]),p_HAC12=unname(w12["p"]),n_obs=length(y),residual_AR1=if(length(u)>2)cor(u[-1],u[-length(u)]) else NA,ljung_box_p=Box.test(u,lag=min(12,floor(length(u)/5)),type="Ljung-Box",fitdf=0)$p.value)}
run_panel<-function(panel,target_name){groups<-split(seq_len(nrow(panel)),panel$grid);rows<-vector("list",length(groups));z<-0L;controls<-c("tp","t2m","ssrd","rzsm","co2","vpd");for(gg in names(groups)){d<-panel[groups[[gg]],c("gpp","wgpe",controls,"lon_val","lat_val")];if(sum(complete.cases(d))<48)next;cc<-lapply(controls,function(v)d[[v]]);fw<-fit_one(d$gpp,d$wgpe,cc);rv<-fit_one(d$wgpe,d$gpp,cc);if(is.null(fw)||is.null(rv))next;z<-z+1L;rows[[z]]<-data.frame(grid_id=gg,lon=mean(d$lon_val),lat=mean(d$lat_val),target=target_name,direction=c("forward","reverse"),rbind(fw,rv),row.names=NULL);if(z%%1000==0)message(target_name," fitted ",z)};do.call(rbind,rows[seq_len(z)])}

panel_cache<-file.path(OUT,"FIG3F_LAG3_PREPROCESSED_PANELS.rds")
if(file.exists(panel_cache)){
  pp<-readRDS(panel_cache);pg<-pp$GPP;pl<-pp$LAI;rm(pp)
}else{
  load(gpp_file);gpp_dates<-seq(as.Date("2000-03-01"),by="month",length.out=298);pg<-augment(panel_gpp,gpp_dates);rm(panel_gpp);gc()
  load(lai_file);lai_dates<-seq(as.Date("2000-01-01"),by="month",length.out=252);pl<-augment(panel_lai,lai_dates);rm(panel_lai);gc()
  saveRDS(list(GPP=pg,LAI=pl),panel_cache,compress="gzip")
}
grid<-rbind(run_panel(pg,"GPP"),run_panel(pl,"LAI"));rm(pg,pl);gc()
for(tg in unique(grid$target))for(di in c("forward","reverse"))for(pn in c("p_F","p_HAC3","p_HAC12")){qn<-sub("p_","q_",pn);if(!qn%in%names(grid))grid[[qn]]<-NA_real_;ii<-grid$target==tg&grid$direction==di;grid[[qn]][ii]<-p.adjust(grid[[pn]][ii],"BH")}
classify<-function(method){f<-grid[grid$direction=="forward",];r<-grid[grid$direction=="reverse",];r<-r[match(paste(f$target,f$grid_id),paste(r$target,r$grid_id)),];q<-paste0("q_",method);data.frame(target=f$target,grid_id=f$grid_id,lon=f$lon,lat=f$lat,method=method,class=ifelse(f[[q]]<=.05&r[[q]]<=.05,"bidirectional",ifelse(f[[q]]<=.05,"forward-only",ifelse(r[[q]]<=.05,"reverse-only","none"))),residual_AR1_forward=f$residual_AR1,ljung_box_p_forward=f$ljung_box_p)}
cls<-do.call(rbind,lapply(c("F","HAC3","HAC12"),classify));write.csv(grid,file.path(OUT,"FIG3F_LAG3_F_HAC_GRID.csv"),row.names=FALSE);write.csv(cls,file.path(OUT,"FIG3F_LAG3_CLASSIFICATION_GRID.csv"),row.names=FALSE)
summ<-list();k<-0L;classes<-c("forward-only","reverse-only","bidirectional","none")
for(tg in unique(cls$target)){base<-cls[cls$target==tg&cls$method=="F",];for(me in c("F","HAC3","HAC12")){d<-cls[cls$target==tg&cls$method==me,];d<-d[match(base$grid_id,d$grid_id),];w<-cos(d$lat*pi/180);for(cl in classes){k<-k+1L;summ[[k]]<-data.frame(target=tg,method=me,class=cl,area_fraction_pct=100*sum(w[d$class==cl])/sum(w),n_grid=nrow(d))};if(me!="F"){k<-k+1L;summ[[k]]<-data.frame(target=tg,method=paste0("F_vs_",me),class="agreement",area_fraction_pct=100*sum(w[d$class==base$class])/sum(w),n_grid=nrow(d))}}}
summary<-do.call(rbind,summ);write.csv(summary,file.path(OUT,"FIG3F_LAG3_F_HAC_SUMMARY.csv"),row.names=FALSE)
diag<-aggregate(cbind(residual_AR1_forward,ljung_box_p_forward)~target+method,data=cls,FUN=function(x)mean(x,na.rm=TRUE));write.csv(diag,file.path(OUT,"FIG3F_LAG3_RESIDUAL_DIAGNOSTICS.csv"),row.names=FALSE)
