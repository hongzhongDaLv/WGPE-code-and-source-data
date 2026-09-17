options(stringsAsFactors=FALSE,warn=1)
ROOT<-Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
OUT<-file.path(Sys.getenv('WGPE_ECOLOGY_OUTPUT', unset=file.path(ROOT,'WGPE_ECOLOGY_NO_VPD_20260904')),'TEMPORAL_F_NO_VPD')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
ctl<-c('precipitation','T2m','SSRD','RZSM','CO2','VPD')
prep<-function(X,mo,tt){
 for(m in 1:12){j<-which(mo==m);X[,j]<-X[,j,drop=FALSE]-rowMeans(X[,j,drop=FALSE],na.rm=TRUE)}
 for(i in seq_len(nrow(X))){ok<-is.finite(X[i,]);if(sum(ok)<10){X[i,]<-NA;next};t<-tt[ok]-mean(tt[ok]);y<-X[i,ok];y<-y-sum(t*(y-mean(y)))/sum(t*t)*t;y<-y-mean(y);s<-sd(y);X[i,ok]<-if(is.finite(s)&&s>0)y/s else NA_real_};X
}
fit<-function(y,R,F){
 n<-length(y);if(n<48)return(c(n=n,p=NA,partial_R2=NA,rank=NA))
 ar<-lm.fit(R,y);af<-lm.fit(F,y)
 if(ar$rank<ncol(R)||af$rank<ncol(F))return(c(n=n,p=NA,partial_R2=NA,rank=0))
 sr<-sum(ar$residuals^2);sf<-sum(af$residuals^2);stat<-max(0,((sr-sf)/3)/(sf/(n-ncol(F))))
 c(n=n,p=pf(stat,3,n-ncol(F),lower.tail=FALSE),partial_R2=(sr-sf)/sr,rank=1)
}
all<-list()
for(target in c('GPP','LAI')){
 D<-readRDS(file.path(ROOT,'output_TOP100_revision/FINAL_COORDINATE_FIXED_CLOSURE/derived_data/ecology',paste0(target,'_CANONICAL_INPUTS.rds')))
 keep<-D$meta$persistent_vegetated_domain&D$meta$target_valid_min36
 meta<-D$meta[keep,];mo<-as.integer(format(D$dates,'%m'));tt<-as.numeric(D$dates)
 stopifnot(all(diff(as.integer(format(D$dates,'%Y'))*12+as.integer(format(D$dates,'%m')))==1))
 for(v in c('y','wgpe_absolute',ctl)){cat(target,'preprocess',v,'\n');D[[v]]<-prep(D[[v]][keep,,drop=FALSE],mo,tt)}
 out<-list();n<-nrow(meta);times<-4:length(D$dates)
 for(mod in c('With_VPD','Without_VPD')){
  ans<-matrix(NA_real_,n,8);colnames(ans)<-c('n_forward','p_forward','partial_R2_forward','rank_forward','n_reverse','p_reverse','partial_R2_reverse','rank_reverse')
  for(i in seq_len(n)){
   y<-D$y[i,];x<-D$wgpe_absolute[i,]
   yl<-sapply(1:3,function(l)y[times-l]);xl<-sapply(1:3,function(l)x[times-l])
   cc<-do.call(cbind,lapply(ctl,function(v)sapply(1:3,function(l)D[[v]][i,times-l])))
   # Six-control valid lag windows fixed before removal of VPD; true calendar lags.
   common<-rowSums(!is.finite(cbind(yl,xl,cc)))==0
   usec<-if(mod=='With_VPD')cc else cc[,1:15,drop=FALSE]
   fw<-common&is.finite(y[times]);rv<-common&is.finite(x[times])
   R<-cbind(1,yl,usec);F<-cbind(R,xl)
   a<-fit(y[times][fw],R[fw,,drop=FALSE],F[fw,,drop=FALSE])
   R<-cbind(1,xl,usec);F<-cbind(R,yl)
   b<-fit(x[times][rv],R[rv,,drop=FALSE],F[rv,,drop=FALSE])
   ans[i,]<-c(a,b)
   if(i%%3000==0)cat(target,mod,i,'/',n,'\n')
  }
  g<-cbind(meta[,c('lon','lat')],target=target,model=mod,as.data.frame(ans))
  ok<-is.finite(g$p_forward)&is.finite(g$p_reverse)
  g$q_forward<-g$q_reverse<-NA_real_;g$q_forward[ok]<-p.adjust(g$p_forward[ok],'BH');g$q_reverse[ok]<-p.adjust(g$p_reverse[ok],'BH')
  g$class<-NA_character_;g$class[ok]<-ifelse(g$q_forward[ok]<=.05,ifelse(g$q_reverse[ok]<=.05,'bidirectional','forward-only'),ifelse(g$q_reverse[ok]<=.05,'reverse-only','none'))
  out[[mod]]<-g
 }
 stopifnot(identical(out[[1]]$n_forward,out[[2]]$n_forward),identical(out[[1]]$n_reverse,out[[2]]$n_reverse))
 all[[target]]<-do.call(rbind,out);write.csv(all[[target]],file.path(OUT,paste0(target,'_GRID.csv')),row.names=FALSE)
 rm(D);gc()
}
grid<-do.call(rbind,all);lev<-c('forward-only','reverse-only','bidirectional','none')
sm<-do.call(rbind,lapply(split(grid,list(grid$target,grid$model),drop=TRUE),function(g){g<-g[!is.na(g$class),];w<-cos(g$lat*pi/180);data.frame(target=g$target[1],model=g$model[1],class=lev,n_tested=nrow(g),n_class=sapply(lev,function(k)sum(g$class==k)),grid_pct=sapply(lev,function(k)100*mean(g$class==k)),area_pct=sapply(lev,function(k)100*sum(w[g$class==k])/sum(w)))}))
write.csv(sm,file.path(OUT,'F_CLASS_SUMMARY.csv'),row.names=FALSE)
old<-read.csv(file.path(ROOT,'output_TOP100_revision/FINAL_COORDINATE_FIXED_CLOSURE/results/ecology/ECOLOGY_TEMPORAL_GRID_CANONICAL.csv'))
old<-old[old$predictor=='W-GPE'&old$height_reference=='absolute'&old$lag==3,]
qa<-do.call(rbind,lapply(c('GPP','LAI'),function(t){g<-grid[grid$target==t&grid$model=='With_VPD',];o<-old[old$target==t,];ii<-match(paste(g$lon,g$lat),paste(o$lon,o$lat));stopifnot(!anyNA(ii));o<-o[ii,];data.frame(target=t,max_p_forward_diff=max(abs(g$p_forward-o$p_forward),na.rm=TRUE),max_p_reverse_diff=max(abs(g$p_reverse-o$p_reverse),na.rm=TRUE),n_forward_mismatch=sum(g$n_forward!=o$n_forward),n_reverse_mismatch=sum(g$n_reverse!=o$n_reverse),class_mismatch=sum(g$class!=o$class_BH,na.rm=TRUE))}))
write.csv(qa,file.path(OUT,'BASELINE_QA.csv'),row.names=FALSE)
print(qa);print(sm)
stopifnot(all(qa$max_p_forward_diff<1e-7),all(qa$max_p_reverse_diff<1e-7),all(qa$class_mismatch==0),all(qa$n_forward_mismatch==0),all(qa$n_reverse_mismatch==0))
cat('PASS baseline and paired samples. Classification only; no new OOF computed.\n')
