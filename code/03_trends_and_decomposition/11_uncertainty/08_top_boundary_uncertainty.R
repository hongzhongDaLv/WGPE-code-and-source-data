options(stringsAsFactors = FALSE)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages({library(ncdf4); library(sandwich)})

ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_reviewer_P0_resolution")
N_BOOT <- 1000L; BLOCK_MONTHS <- 60L; set.seed(20260713)

decimal_year <- function(d) {
  y <- as.integer(format(d, "%Y")); y0 <- as.Date(sprintf("%04d-01-01", y)); y1 <- as.Date(sprintf("%04d-01-01", y+1L))
  y + as.numeric(d-y0)/as.numeric(y1-y0)
}
deseason <- function(y) y - ave(y, rep(1:12, length.out=length(y)), FUN=mean)
slope <- function(y, x) {
  xc <- x-mean(x); sum(xc*(y-mean(y)))/sum(xc^2)
}
mbb_residual_slopes <- function(y, x, B=1000L, L=60L) {
  b <- slope(y,x); a0 <- mean(y)-b*mean(x); f <- a0+b*x; r <- y-f; r <- r-mean(r); n <- length(r); starts <- seq_len(n-L+1L)
  replicate(B, {
    ss <- sample(starts, ceiling(n/L), replace=TRUE)
    rb <- unlist(lapply(ss, function(s) r[s:(s+L-1L)]), use.names=FALSE)[seq_len(n)]
    slope(f+rb, x)
  })
}
ci2 <- function(v) unname(quantile(v, c(.025,.975), na.rm=TRUE, type=7))
fit_ci <- function(y,x) {
  fit <- lm(y~x); est <- unname(coef(fit)[2]); ols <- unname(confint(fit,"x",level=.95))
  V <- NeweyWest(fit, lag=NULL, prewhite=FALSE, adjust=TRUE); se <- sqrt(V[2,2]); nw <- est+c(-1,1)*qt(.975,df.residual(fit))*se
  mb <- ci2(mbb_residual_slopes(y,x))
  c(est=est,ols_lo=ols[1],ols_hi=ols[2],nw_lo=nw[1],nw_hi=nw[2],mbb_lo=mb[1],mbb_hi=mb[2])
}

z <- read.csv(file.path(OUT,"TOP_BOUNDARY_ZONE_MONTHLY_SERIES.csv"), check.names=FALSE)
names(z)[1] <- "date"
z$date <- as.Date(z$date); z$x <- decimal_year(z$date)

nc <- nc_open(file.path(OUT,"TOP_BOUNDARY_10DEG_BLOCK_MONTHLY.nc"))
a <- ncvar_get(nc,"block_area_mean") # lonblock, latblock, time, variable, top
bw <- as.vector(ncvar_get(nc,"block_area_weight"))
latmid <- ncvar_get(nc,"lat_block_mid")
topvals <- ncvar_get(nc,"top_hPa"); varvals <- ncvar_get(nc,"variable_code")
nc_close(nc)
blocklat <- rep(latmid, each=36)

zone_specs <- list(
  Global=c(-90,90), `S high latitude`=c(-90,-66.5), `S temperate`=c(-66.5,-23.5),
  Tropics=c(-23.5,23.5), `N temperate`=c(23.5,66.5), `N high latitude`=c(66.5,90)
)
zone_block <- function(region) {
  b <- zone_specs[[region]]; ok <- blocklat>=b[1]&blocklat<=b[2]
  if(region!="Global" && b[2]!=90) ok <- ok & blocklat<b[2]
  which(ok)
}
block_matrix <- function(top,var) {
  it <- match(top,topvals); iv <- match(var,varvals)
  t(matrix(a[,,,iv,it], nrow=36*18, ncol=552)) # time x block
}
block_deseason <- function(m) {
  out <- m
  for(mm in 1:12) {ii <- seq(mm,nrow(m),12); out[ii,] <- sweep(m[ii,,drop=FALSE],2,colMeans(m[ii,,drop=FALSE]),"-")}
  out
}
spatial_and_combined <- function(mat,x,idx,B=1000L) {
  mat <- mat[,idx,drop=FALSE]; ww <- bw[idx]
  dm <- block_deseason(mat)
  block_slopes <- apply(dm,2,slope,x=x)
  sp <- numeric(B); comb <- numeric(B); n <- length(idx)
  for(bb in seq_len(B)) {
    jj <- sample(seq_len(n),n,replace=TRUE); wj <- ww[jj]
    sp[bb] <- sum(block_slopes[jj]*wj)/sum(wj)
    ys <- rowSums(sweep(mat[,jj,drop=FALSE],2,wj,"*"))/sum(wj)
    ys <- deseason(ys)
    comb[bb] <- sample(mbb_residual_slopes(ys,x,B=1L),1)
  }
  c(spatial_lo=ci2(sp)[1],spatial_hi=ci2(sp)[2],combined_lo=ci2(comb)[1],combined_hi=ci2(comb)[2])
}

rows <- list(); rid <- 1L
for(top in c(300,200,100)) for(var in c("iwv","zbar","wgpe")) for(region in names(zone_specs)) {
  d <- z[z$top_hPa==top & z$variable==var & z$region==region,]
  fc <- fit_ci(d$area_mean_deseasonalized,d$x)
  sc <- spatial_and_combined(block_matrix(top,var),d$x,zone_block(region),B=N_BOOT)
  rows[[length(rows)+1L]] <- data.frame(
    result_id=sprintf("P0TU%05d",rid),result_type="top_version",top_hPa=top,reference_top_hPa=NA,
    variable=var,region=region,estimate=fc["est"],ols_ci_low=fc["ols_lo"],ols_ci_high=fc["ols_hi"],
    newey_west_ci_low=fc["nw_lo"],newey_west_ci_high=fc["nw_hi"],
    moving_block_5yr_ci_low=fc["mbb_lo"],moving_block_5yr_ci_high=fc["mbb_hi"],
    spatial_10deg_ci_low=sc["spatial_lo"],spatial_10deg_ci_high=sc["spatial_hi"],
    combined_ci_low=sc["combined_lo"],combined_ci_high=sc["combined_hi"],
    n_month=552,bootstrap_replicates=N_BOOT,time_block_months=BLOCK_MONTHS,
    notes="Newey-West automatic lag; residual moving-block bootstrap; fixed 10-degree spatial blocks; combined spatial resampling plus residual time-block bootstrap",
    check.names=FALSE
  ); rid <- rid+1L
}

# Paired top-boundary differences, preserving month and spatial block pairing.
for(pair in list(c(300,100),c(200,100),c(300,200))) for(var in c("iwv","zbar","wgpe")) for(region in names(zone_specs)) {
  ta<-pair[1]; tb<-pair[2]
  da<-z[z$top_hPa==ta&z$variable==var&z$region==region,]; db<-z[z$top_hPa==tb&z$variable==var&z$region==region,]
  stopifnot(all(da$date==db$date)); yd<-da$area_mean_deseasonalized-db$area_mean_deseasonalized
  fc<-fit_ci(yd,da$x)
  md<-block_matrix(ta,var)-block_matrix(tb,var); sc<-spatial_and_combined(md,da$x,zone_block(region),B=N_BOOT)
  rows[[length(rows)+1L]]<-data.frame(
    result_id=sprintf("P0TU%05d",rid),result_type="paired_difference",top_hPa=ta,reference_top_hPa=tb,
    variable=var,region=region,estimate=fc["est"],ols_ci_low=fc["ols_lo"],ols_ci_high=fc["ols_hi"],
    newey_west_ci_low=fc["nw_lo"],newey_west_ci_high=fc["nw_hi"],
    moving_block_5yr_ci_low=fc["mbb_lo"],moving_block_5yr_ci_high=fc["mbb_hi"],
    spatial_10deg_ci_low=sc["spatial_lo"],spatial_10deg_ci_high=sc["spatial_hi"],
    combined_ci_low=sc["combined_lo"],combined_ci_high=sc["combined_hi"],
    n_month=552,bootstrap_replicates=N_BOOT,time_block_months=BLOCK_MONTHS,
    notes="Paired by common month and common spatial-block resample; Newey-West automatic lag",
    check.names=FALSE
  ); rid<-rid+1L
}

res <- do.call(rbind,rows)
write.csv(res,file.path(OUT,"TOP_BOUNDARY_UNCERTAINTY.csv"),row.names=FALSE,fileEncoding="UTF-8")
cat("rows",nrow(res),"\n")
