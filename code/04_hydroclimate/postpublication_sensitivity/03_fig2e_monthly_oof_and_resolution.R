options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "8")
Sys.setenv(
  PATH = paste("__WGPE_RTOOLS__/usr/bin",
               "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin",
               Sys.getenv("PATH"), sep = .Platform$path.sep),
  BINPREF = "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin/",
  COMPILER_PATH = "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin",
  LIBRARY_PATH = "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/lib"
)
Sys.setenv(MAKEFLAGS = paste(
  "SHELL=__WGPE_RTOOLS__/usr/bin/sh.exe",
  "PATH=/usr/bin:/x86_64-w64-mingw32.static.posix/bin"
))

ROOT <- "__WGPE_PROJECT_ROOT__"
CANONICAL <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
OUT <- file.path(CANONICAL, "08_REPRODUCTION_OUTPUT", "Figure2e")
PRED <- file.path(OUT, "FIG2E_MONTHLY_OOF_RDS")
dir.create(PRED, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts_TOP100_revision", "18_coordinate_harmonization", "coordinate_harmonization.R"))
stopifnot(requireNamespace("Rcpp", quietly = TRUE), requireNamespace("RcppArmadillo", quietly = TRUE))
Rcpp::sourceCpp(file.path(CANONICAL, "07_SCRIPTS", "support",
                          "fig2e_oof_predictions_training_only.cpp"),
                rebuild = TRUE, showOutput = FALSE)

core_file <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                       "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file <- file.path(ROOT, "output_attribution_minimal", "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
nino_file <- file.path(ROOT, "WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL", "11_ENSO_ONI_STANDARD_UPDATE",
                       "00_external_data", "NOAA_PSL_Nino34_monthly_clean.csv")
e <- new.env(parent = emptyenv()); load(core_file, envir = e)
lon <- as.numeric(e$lon); lat <- as.numeric(e$lat); dates <- as.Date(e$all_dates)
nlon <- length(lon); nlat <- length(lat); nt <- length(dates); ncell <- nlon*nlat
cell <- expand.grid(lon = lon, lat = lat); cell$cell_id <- seq_len(ncell); cell$weight <- cos(cell$lat*pi/180)
to_mat <- function(a) matrix(as.numeric(a), nrow = ncell, ncol = nt)
I <- to_mat(e$iwv_st); Z <- to_mat(e$zbar_st); W <- to_mat(e$wgpe_st); rm(e); gc()
s <- readRDS(single_file)
P <- to_mat(reorder_exact_grid(s$tp_native_m, s$lon, s$lat, lon, lat)); rm(s); gc()

mon <- as.integer(format(dates, "%m")); yr <- as.integer(format(dates, "%Y")); tdec <- yr+(mon-.5)/12
prep <- function(m) {
  good <- rowSums(is.finite(m)) == ncol(m); m[!good,] <- NA_real_
  for (mm in 1:12) { id <- which(mon == mm); cl <- rowMeans(m[,id,drop=FALSE],na.rm=TRUE); m[,id] <- sweep(m[,id,drop=FALSE],1,cl,"-") }
  tc <- tdec-mean(tdec); mu <- rowMeans(m,na.rm=TRUE); sl <- as.numeric(sweep(m,1,mu,"-") %*% tc / sum(tc^2))
  m <- m-outer(sl,tc); mu <- rowMeans(m,na.rm=TRUE); m <- sweep(m,1,mu,"-")
  sdv <- sqrt(rowSums(m^2,na.rm=TRUE)/pmax(1,rowSums(is.finite(m))-1)); ok <- is.finite(sdv)&sdv>0
  m[ok,] <- m[ok,,drop=FALSE]/sdv[ok]; m[!ok,] <- NA_real_; m
}
message("Preprocessing common anomaly/detrend inputs")
I <- prep(I); Z <- prep(Z); W <- prep(W); P <- prep(P); Int <- I*Z; gc()
nino <- read.csv(nino_file); stopifnot(identical(as.character(as.Date(nino$date)), as.character(dates)))
nino_y <- as.numeric(scale(as.numeric(nino$nino34_standardized)))
fold <- as.integer((yr-min(yr)) %/% 5)+1L; gap <- 3L

index_rows <- list(); ik <- 0L
run_target <- function(name, Y, batch=2500L) {
  for (st in seq.int(1L,ncell,by=batch)) {
    en <- min(ncell,st+batch-1L); id <- st:en
    fn <- file.path(PRED, sprintf("%s_cells_%06d_%06d.rds", gsub(" ","_",name),st,en))
    if (!file.exists(fn)) {
      message(name," ",st,"-",en)
      YY <- if (is.null(dim(Y))) matrix(Y,nrow=length(id),ncol=nt,byrow=TRUE) else Y[id,,drop=FALSE]
      oo <- cv_oof_predictions_training_only_cpp(YY,I[id,,drop=FALSE],Z[id,,drop=FALSE],Int[id,,drop=FALSE],W[id,,drop=FALSE],fold,gap)
      obj <- c(list(cell=cell[id,],dates=dates,fold=fold,gap_months=gap,
                    preprocessing="calendar-month anomaly and detrend; fold-specific predictor/target standardization; fold-specific M3 W-GPE orthogonalization"),oo)
      saveRDS(obj,fn,compress="gzip")
      rm(oo,obj,YY); gc()
    }
    ik <<- ik+1L; index_rows[[ik]] <<- data.frame(target=name,file=fn,cell_start=st,cell_end=en,n_cell=length(id),n_time=nt)
  }
}
run_target("Precipitation",P)
run_target("Continuous Nino3.4",nino_y)
idx <- do.call(rbind,index_rows); write.csv(idx,file.path(OUT,"FIG2E_MONTHLY_OOF_INDEX.csv"),row.names=FALSE)

band <- function(x) ifelse(x < -66.5,"S high latitude",ifelse(x < -23.5,"S temperate",ifelse(x <= 23.5,"Tropics",ifelse(x <= 66.5,"N temperate","N high latitude"))))
wquant <- function(x,w,p) { ok <- is.finite(x)&is.finite(w)&w>0; x<-x[ok];w<-w[ok];o<-order(x);x<-x[o];w<-w[o]; x[pmax(1,findInterval(p*sum(w),cumsum(w))+1)] }
perf <- function(y,p) { ok<-is.finite(y)&is.finite(p); y<-y[ok];p<-p[ok]; c(r2=1-sum((y-p)^2)/sum((y-mean(y))^2),rmse=sqrt(mean((y-p)^2))) }

grid_rows <- list(); gk <- 0L; pool_acc <- list()
for (ii in seq_len(nrow(idx))) {
  o <- readRDS(idx$file[ii]); target <- idx$target[ii]
  for (j in seq_len(nrow(o$cell))) {
    yy <- o$observed[j,]; pp <- lapply(0:3,function(m)o[[paste0("pred_M",m)]][j,]); pm <- lapply(pp,function(x)perf(yy,x))
    gk <- gk+1L; grid_rows[[gk]] <- data.frame(target=target,o$cell[j,],climate_band=band(o$cell$lat[j]),
      M0_cv_r2=pm[[1]]["r2"],M1_cv_r2=pm[[2]]["r2"],M2_cv_r2=pm[[3]]["r2"],M3_cv_r2=pm[[4]]["r2"],
      M0_rmse=pm[[1]]["rmse"],M1_rmse=pm[[2]]["rmse"],M2_rmse=pm[[3]]["rmse"],M3_rmse=pm[[4]]["rmse"])
  }
}
grid <- do.call(rbind,grid_rows); write.csv(grid,file.path(OUT,"FIG2E_GRID_OOF_RESOLUTION.csv"),row.names=FALSE)

comparisons <- list(`M1-M0`=c(1,0),`M2-M0`=c(2,0),`M3-M0`=c(3,0))
summary_rows <- list(); sk<-0L
for (target in unique(grid$target)) for (reg in c("Global","N high latitude","N temperate","Tropics","S temperate","S high latitude")) {
  d<-grid[grid$target==target & (reg=="Global" | grid$climate_band==reg),]
  for (nm in names(comparisons)) { a<-comparisons[[nm]][1];b<-comparisons[[nm]][2]; dr<-d[[paste0("M",a,"_cv_r2")]]-d[[paste0("M",b,"_cv_r2")]]; drm<-d[[paste0("M",b,"_rmse")]]-d[[paste0("M",a,"_rmse")]]; w<-d$weight; ok<-is.finite(dr)&is.finite(w)&w>0
    sk<-sk+1L; summary_rows[[sk]]<-data.frame(target=target,comparison=nm,region=reg,n_grid=sum(ok),weighted_mean_delta_r2=sum(dr[ok]*w[ok])/sum(w[ok]),weighted_median_delta_r2=wquant(dr,w,.5),q25=wquant(dr,w,.25),q75=wquant(dr,w,.75),positive_area_pct=100*sum(w[ok&dr>0])/sum(w[ok]),negative_area_pct=100*sum(w[ok&dr<0])/sum(w[ok]),weighted_mean_rmse_improvement=sum(drm[ok]*w[ok])/sum(w[ok]),relative_rmse_reduction_pct=100*sum((drm/d[[paste0("M",b,"_rmse")]])[ok]*w[ok])/sum(w[ok]))
  }
}
summ<-do.call(rbind,summary_rows); write.csv(summ,file.path(OUT,"FIG2E_FINAL_RESOLUTION_SUMMARY.csv"),row.names=FALSE)

writeLines(c("# Figure 2e monthly OOF QA","","- OOF predictions are stored as batch RDS files indexed by `FIG2E_MONTHLY_OOF_INDEX.csv`.","- All M0–M3 models use identical complete cases, folds and a 3-month gap.","- Predictor and target standardization is estimated from the training fold only.","- M3 W-GPE orthogonalization against IWV is estimated from the training fold only.","- The orthogonalized M3 spans the same two-predictor linear space as IWV + W-GPE; it improves numerical attribution without changing the fitted estimand.","- Outputs are same-month out-of-sample associations, not operational forecasts and not causal estimates."),file.path(OUT,"FIG2E_MONTHLY_OOF_QA.md"))
