options(stringsAsFactors = FALSE)
ROOT <- "__WGPE_PROJECT_ROOT__"
CANONICAL <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
OUT <- file.path(CANONICAL, "08_REPRODUCTION_OUTPUT", "Figure2e")
idx <- read.csv(file.path(OUT,"FIG2E_MONTHLY_OOF_INDEX.csv"),check.names=FALSE)
grid <- read.csv(file.path(OUT,"FIG2E_GRID_OOF_RESOLUTION.csv"),check.names=FALSE)
comparisons <- list(`M1-M0`=c(1,0),`M2-M0`=c(2,0),`M3-M0`=c(3,0))

fold_rows <- list(); fk<-0L; pooled_rows<-list(); pk<-0L
for(target in unique(idx$target)) {
  z <- idx[idx$target==target,]
  # sufficient statistics for pooled within-grid-centred OOF R2
  ss <- setNames(vector("list",4),paste0("M",0:3)); for(m in 0:3) ss[[m+1]]<-c(sse=0,sst=0,n=0,wobs=0)
  for(ii in seq_len(nrow(z))) {
    o<-readRDS(z$file[ii]); Y<-o$observed; w<-o$cell$weight; folds<-sort(unique(o$fold))
    for(ff in folds) {
      tt<-which(o$fold==ff)
      y<-Y[,tt,drop=FALSE]; ym<-rowMeans(y,na.rm=TRUE); sst<-rowSums((y-ym)^2,na.rm=TRUE)
      for(nm in names(comparisons)) {
        a<-comparisons[[nm]][1]; b<-comparisons[[nm]][2]
        pa<-o[[paste0("pred_M",a)]][,tt,drop=FALSE]; pb<-o[[paste0("pred_M",b)]][,tt,drop=FALSE]
        da<-(rowSums((y-pb)^2,na.rm=TRUE)-rowSums((y-pa)^2,na.rm=TRUE))/sst
        ok<-is.finite(da)&is.finite(w)&w>0
        fk<-fk+1L; fold_rows[[fk]]<-data.frame(target=target,comparison=nm,fold=ff,
          numerator=sum(da[ok]*w[ok]),denominator=sum(w[ok]),n_grid=sum(ok))
      }
    }
    ymfull<-rowMeans(Y,na.rm=TRUE); sstfull<-rowSums((Y-ymfull)^2,na.rm=TRUE)
    for(m in 0:3) {
      P<-o[[paste0("pred_M",m)]]; sse<-rowSums((Y-P)^2,na.rm=TRUE); ok<-is.finite(sse)&is.finite(sstfull)&is.finite(w)&w>0
      ss[[m+1]]["sse"]<-ss[[m+1]]["sse"]+sum(w[ok]*sse[ok]); ss[[m+1]]["sst"]<-ss[[m+1]]["sst"]+sum(w[ok]*sstfull[ok]); ss[[m+1]]["n"]<-ss[[m+1]]["n"]+sum(ok); ss[[m+1]]["wobs"]<-ss[[m+1]]["wobs"]+sum(w[ok])*ncol(Y)
    }
  }
  for(m in 0:3) { pk<-pk+1L; pooled_rows[[pk]]<-data.frame(target=target,model=paste0("M",m),pooled_space_time_oof_r2=1-ss[[m+1]]["sse"]/ss[[m+1]]["sst"],pooled_rmse=sqrt(ss[[m+1]]["sse"]/ss[[m+1]]["wobs"]),n_grid=ss[[m+1]]["n"]) }
}
fold_raw<-do.call(rbind,fold_rows)
fold_sum<-aggregate(cbind(numerator,denominator,n_grid)~target+comparison+fold,fold_raw,sum)
fold_sum$estimate<-fold_sum$numerator/fold_sum$denominator
write.csv(fold_sum,file.path(OUT,"FIG2E_TEMPORAL_FOLD_ESTIMATES.csv"),row.names=FALSE)
pooled<-do.call(rbind,pooled_rows)
for(target in unique(pooled$target)) for(nm in names(comparisons)) {
  a<-paste0("M",comparisons[[nm]][1]);b<-paste0("M",comparisons[[nm]][2]); ia<-pooled$target==target&pooled$model==a;ib<-pooled$target==target&pooled$model==b
  pooled$comparison[ia]<-nm; pooled$pooled_delta_r2[ia]<-pooled$pooled_space_time_oof_r2[ia]-pooled$pooled_space_time_oof_r2[ib]
  pooled$pooled_rmse_improvement[ia]<-pooled$pooled_rmse[ib]-pooled$pooled_rmse[ia]
}
write.csv(pooled,file.path(OUT,"FIG2E_POOLED_SPACE_TIME_OOF.csv"),row.names=FALSE)

# Separate paired temporal-fold and paired 10°x10° spatial-block bootstrap CIs.
set.seed(20260722); B<-2000L; out<-list();k<-0L
for(target in unique(grid$target)) for(nm in names(comparisons)) {
  fs<-fold_sum[fold_sum$target==target&fold_sum$comparison==nm,]
  tv<-rep(NA_real_,B); for(bi in seq_len(B)) tv[bi]<-mean(sample(fs$estimate,nrow(fs),replace=TRUE))
  d<-grid[grid$target==target,]; a<-comparisons[[nm]][1];bb<-comparisons[[nm]][2]
  dv<-d[[paste0("M",a,"_cv_r2")]]-d[[paste0("M",bb,"_cv_r2")]]; ok<-is.finite(dv)&is.finite(d$weight)&d$weight>0
  bid<-paste(floor(((d$lon[ok]+180)%%360)/10),floor((d$lat[ok]+90)/10),sep="_")
  bn<-tapply(dv[ok]*d$weight[ok],bid,sum); bd<-tapply(d$weight[ok],bid,sum); ids<-names(bn)
  sv<-rep(NA_real_,B); for(bi in seq_len(B)){q<-sample(ids,length(ids),replace=TRUE);sv[bi]<-sum(bn[q])/sum(bd[q])}
  k<-k+1L; out[[k]]<-data.frame(target=target,comparison=nm,
    temporal_estimate=mean(fs$estimate),temporal_lower95=quantile(tv,.025),temporal_upper95=quantile(tv,.975),
    spatial_estimate=sum(bn)/sum(bd),spatial_lower95=quantile(sv,.025),spatial_upper95=quantile(sv,.975),
    n_fold=nrow(fs),n_spatial_block=length(ids),B=B)
}
write.csv(do.call(rbind,out),file.path(OUT,"FIG2E_PAIRED_TEMPORAL_SPATIAL_BOOTSTRAP.csv"),row.names=FALSE)
