options(stringsAsFactors=FALSE)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
suppressPackageStartupMessages(library(ncdf4))
ROOT<-"__WGPE_PROJECT_ROOT__";OUT<-file.path(ROOT,"output_reviewer_P0_resolution")
tops<-c(300,200,100);vars<-c("iwv","zbar","wgpe");B<-1000L;L<-60L;set.seed(20260713)
dates<-seq(as.Date("1979-01-01"),as.Date("2024-12-01"),by="month")
x<-{y<-as.integer(format(dates,"%Y"));y0<-as.Date(sprintf("%04d-01-01",y));y1<-as.Date(sprintf("%04d-01-01",y+1));y+as.numeric(dates-y0)/as.numeric(y1-y0)}
slope<-function(y){xc<-x-mean(x);sum(xc*(y-mean(y)))/sum(xc^2)}
deseason<-function(y)y-ave(y,rep(1:12,length.out=length(y)),FUN=mean)
mbb<-function(y){b<-slope(y);f<-mean(y)-b*mean(x)+b*x;r<-y-f;r<-r-mean(r);starts<-1:(length(y)-L+1)
 replicate(B,{ss<-sample(starts,ceiling(length(y)/L),TRUE);rb<-unlist(lapply(ss,function(s)r[s:(s+L-1)]),use.names=FALSE)[1:length(y)];slope(f+rb)})}
series<-list();sid<-1
for(top in tops){nc<-nc_open(file.path(OUT,sprintf("ERA5_TOP%d_monthly.nc",top)));lat<-ncvar_get(nc,"latitude");wlat<-cos(lat*pi/180)
 for(v in vars){a<-ncvar_get(nc,v) # lon lat time
  for(lo in seq(-90,85,5)){hi<-lo+5;ii<-which(lat>=lo & if(hi<90)lat<hi else lat<=hi);ww<-matrix(rep(wlat[ii],each=360),nrow=360)
   y<-vapply(1:552,function(tt)weighted.mean(as.vector(a[,ii,tt]),as.vector(ww),na.rm=TRUE),numeric(1));yd<-deseason(y)
   series[[length(series)+1]]<-data.frame(top_hPa=top,variable=v,lat_start=lo,lat_end=hi,lat_mid=lo+2.5,date=dates,value_raw=y,value_deseasonalized=yd)
  }
 }
 nc_close(nc)
}
s<-do.call(rbind,series);write.csv(s,file.path(OUT,"TOP_BOUNDARY_LAT5_MONTHLY_SERIES.csv"),row.names=FALSE)
rows<-list();rid<-1
for(top in tops)for(v in vars)for(mid in sort(unique(s$lat_mid))){d<-s[s$top_hPa==top&s$variable==v&s$lat_mid==mid,];bt<-mbb(d$value_deseasonalized)
 rows[[length(rows)+1]]<-data.frame(result_id=sprintf("P0TLP%05d",rid),result_type="top_version",top_hPa=top,reference_top_hPa=NA,variable=v,lat_mid=mid,
  estimate=slope(d$value_deseasonalized),ci_low=unname(quantile(bt,.025)),ci_high=unname(quantile(bt,.975)),bootstrap_replicates=B,time_block_months=L);rid<-rid+1}
for(pa in list(c(300,100),c(200,100),c(300,200)))for(v in vars)for(mid in sort(unique(s$lat_mid))){a<-s[s$top_hPa==pa[1]&s$variable==v&s$lat_mid==mid,];b<-s[s$top_hPa==pa[2]&s$variable==v&s$lat_mid==mid,];stopifnot(all(a$date==b$date));bt<-mbb(a$value_deseasonalized-b$value_deseasonalized)
 rows[[length(rows)+1]]<-data.frame(result_id=sprintf("P0TLP%05d",rid),result_type="paired_difference",top_hPa=pa[1],reference_top_hPa=pa[2],variable=v,lat_mid=mid,
  estimate=slope(a$value_deseasonalized-b$value_deseasonalized),ci_low=unname(quantile(bt,.025)),ci_high=unname(quantile(bt,.975)),bootstrap_replicates=B,time_block_months=L);rid<-rid+1}
write.csv(do.call(rbind,rows),file.path(OUT,"TOP_BOUNDARY_LAT5_UNCERTAINTY.csv"),row.names=FALSE)
cat("done\n")
