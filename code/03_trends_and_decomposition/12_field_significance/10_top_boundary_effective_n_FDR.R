options(stringsAsFactors=FALSE)
.libPaths(c("__WGPE_R_LIBRARY__",.libPaths()))
suppressPackageStartupMessages(library(ncdf4))
ROOT<-"__WGPE_PROJECT_ROOT__";OUT<-file.path(ROOT,"output_reviewer_P0_resolution")
tops<-c(300,200,100);vars<-c("iwv","zbar","wgpe")
dates<-seq(as.Date("1979-01-01"),as.Date("2024-12-01"),by="month")
decimal_year<-function(d){y<-as.integer(format(d,"%Y"));y0<-as.Date(sprintf("%04d-01-01",y));y1<-as.Date(sprintf("%04d-01-01",y+1L));y+as.numeric(d-y0)/as.numeric(y1-y0)}

calc<-function(arr){
 d<-dim(arr);nlon<-d[1];nlat<-d[2];nt<-d[3];ncell<-nlon*nlat;dim(arr)<-c(ncell,nt)
 x<-decimal_year(dates);xc<-x-mean(x);sxx<-sum(xc^2);month<-as.integer(format(dates,"%m"))
 trend<-rho<-ne<-pe<-rep(NA_real_,ncell)
 for(start in seq(1,ncell,by=5000)){ii<-start:min(ncell,start+4999);y<-arr[ii,,drop=FALSE];a<-y
  for(m in 1:12){jj<-which(month==m);a[,jj]<-y[,jj,drop=FALSE]-rowMeans(y[,jj,drop=FALSE])}
  b<-drop(a%*%xc)/sxx;intercept<-rowMeans(a)-b*mean(x);res<-a-(intercept+outer(b,x));sse<-rowSums(res^2)
  u<-res[,1:(nt-1),drop=FALSE];v<-res[,2:nt,drop=FALSE];nl<-nt-1;su<-rowSums(u);sv<-rowSums(v)
  num<-rowSums(u*v)-su*sv/nl;den<-sqrt((rowSums(u^2)-su^2/nl)*(rowSums(v^2)-sv^2/nl));rr<-num/den;rr[!is.finite(rr)]<-0;rr<-pmax(-.99,pmin(.99,rr))
  nn<-pmin(nt,pmax(4,nt*(1-rr)/(1+rr)));df<-nn-2;se<-sqrt((sse/df)/sxx);pp<-2*pt(abs(b/se),df=df,lower.tail=FALSE)
  trend[ii]<-b;rho[ii]<-rr;ne[ii]<-nn;pe[ii]<-pp
 }
 q<-p.adjust(pe,method="BH");to_mat<-function(v)matrix(v,nrow=nlon,ncol=nlat)
 list(trend=to_mat(trend),residual_ar1=to_mat(rho),n_eff=to_mat(ne),p=to_mat(pe),q=to_mat(q))
}
res<-list();sumrows<-list();rid<-1
for(top in tops){nc<-nc_open(file.path(OUT,sprintf("ERA5_TOP%d_monthly.nc",top)))
 for(v in vars){a<-ncvar_get(nc,v);r<-calc(a);res[[paste0("TOP",top,"_",v)]]<-r
  lat<-ncvar_get(nc,"latitude");w<-matrix(rep(cos(lat*pi/180),each=360),nrow=360,ncol=181)
  pos<-r$q<.05&r$trend>0;neg<-r$q<.05&r$trend<0
  sumrows[[length(sumrows)+1]]<-data.frame(result_id=sprintf("P0TS%05d",rid),top_hPa=top,variable=v,n_grid=sum(is.finite(r$q)),
   positive_grid_pct=100*mean(pos),negative_grid_pct=100*mean(neg),positive_area_pct=100*sum(w[pos])/sum(w),negative_area_pct=100*sum(w[neg])/sum(w),
   median_n_eff=median(r$n_eff,na.rm=TRUE),method="residual-AR1 effective N; per-variable BH q<0.05");rid<-rid+1
 }
 nc_close(nc)
}
saveRDS(res,file.path(OUT,"TOP_BOUNDARY_EFFECTIVE_N_FDR.rds"),compress="xz")
write.csv(do.call(rbind,sumrows),file.path(OUT,"TOP_BOUNDARY_EFFECTIVE_N_FDR_SUMMARY.csv"),row.names=FALSE)

pairs<-list(c(300,100),c(200,100),c(300,200));pr<-list();rid<-1
for(pa in pairs)for(v in vars){a<-res[[paste0("TOP",pa[1],"_",v)]];b<-res[[paste0("TOP",pa[2],"_",v)]];sa<-a$q<.05;sb<-b$q<.05;both<-sa&sb;union<-sa|sb
 pr[[length(pr)+1]]<-data.frame(result_id=sprintf("P0TSG%05d",rid),comparison=sprintf("TOP%d-TOP%d",pa[1],pa[2]),variable=v,
  n_significant_both=sum(both),significant_both_sign_agreement_pct=if(sum(both))100*mean(sign(a$trend[both])==sign(b$trend[both])) else NA,
  n_significant_union=sum(union),significant_status_and_sign_agreement_pct=100*mean((sa==sb)&(!both|sign(a$trend)==sign(b$trend))),
  definition="both: denominator cells significant in both versions; status-and-sign: all cells, requires same significance status and same sign when significant");rid<-rid+1}
write.csv(do.call(rbind,pr),file.path(OUT,"TOP_BOUNDARY_SIGNIFICANT_SIGN_AGREEMENT.csv"),row.names=FALSE)
cat("done\n")
