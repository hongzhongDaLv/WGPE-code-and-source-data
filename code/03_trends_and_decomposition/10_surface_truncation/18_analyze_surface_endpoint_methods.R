options(stringsAsFactors = FALSE)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages(library(ncdf4))

ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_reviewer_P0_resolution")
G <- 9.80665
P0_FILE <- file.path(OUT, "ERA5_TOP100_monthly.nc")
PX_FILE <- file.path(OUT, "ERA5_TOP100_SURFACE_ENDPOINT_METHODS_monthly.nc")

to_lon_lat_time <- function(nc, vn) {
  x <- ncvar_get(nc, vn)
  dn <- vapply(nc$var[[vn]]$dim, function(d) d$name, character(1))
  target <- match(c("longitude", "latitude", "time"), dn)
  if (anyNA(target)) stop(sprintf("Unexpected dimensions for %s: %s", vn, paste(dn, collapse=",")))
  aperm(x, target)
}

deseason <- function(m, month) {
  out <- m
  for (mm in 1:12) {
    ii <- which(month == mm)
    out[, ii] <- out[, ii, drop=FALSE] - rowMeans(out[, ii, drop=FALSE], na.rm=TRUE)
  }
  out
}

trend_cells <- function(m, year, month) {
  a <- deseason(m, month)
  tc <- year - mean(year)
  den <- sum(tc^2)
  ok <- rowSums(is.finite(a)) == ncol(a)
  out <- rep(NA_real_, nrow(a))
  out[ok] <- as.vector(a[ok,,drop=FALSE] %*% tc / den)
  out
}

weighted_mean_safe <- function(x, w, mask = rep(TRUE, length(x))) {
  ok <- mask & is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

nc0 <- nc_open(P0_FILE)
ncx <- nc_open(PX_FILE)
lon <- nc0$dim$longitude$vals
lat <- nc0$dim$latitude$vals
time_vals <- nc0$dim$time$vals
time_units <- nc0$dim$time$units
origin <- sub("^[A-Za-z]+ since ", "", time_units)
dates <- as.POSIXct(origin, tz="UTC") + time_vals * ifelse(grepl("days", time_units), 86400, 1)
year <- as.integer(format(dates, "%Y")) + (as.integer(format(dates, "%m")) - 0.5) / 12
month <- as.integer(format(dates, "%m"))
cell <- expand.grid(lon = lon, lat = lat)
cell$w <- cos(cell$lat*pi/180)

band_levels <- c("Global","S high latitude","S temperate","Tropics","N temperate","N high latitude")
band_mask <- function(reg) {
  if (reg=="Global") rep(TRUE,nrow(cell)) else if (reg=="S high latitude") cell$lat < -66.5 else
    if (reg=="S temperate") cell$lat >= -66.5 & cell$lat < -23.5 else if (reg=="Tropics") cell$lat >= -23.5 & cell$lat <= 23.5 else
      if (reg=="N temperate") cell$lat > 23.5 & cell$lat <= 66.5 else cell$lat > 66.5
}

methods <- c(P0_current_crossing="P0", P1_surface_anchored="P1", P2_one_sided="P2")
variables <- c(iwv="kg m-2 yr-1", zbar="m yr-1", wgpe="J m-2 yr-1")
trends <- list(); globals <- list(); monthly <- list(); idx_result <- 0L

for (mn in names(methods)) for (vn in names(variables)) {
  idx_result <- idx_result + 1L
  arr <- if (mn == "P0_current_crossing") to_lon_lat_time(nc0, vn) else to_lon_lat_time(ncx, paste0(methods[[mn]], "_", vn))
  mat <- matrix(arr, nrow=length(lon)*length(lat), ncol=length(dates))
  tr <- trend_cells(mat, year, month)
  trends[[paste(mn,vn,sep="_")]] <- tr
  gm <- vapply(seq_along(dates), function(tt) weighted_mean_safe(mat[,tt], cell$w), numeric(1))
  monthly[[paste(mn,vn,sep="_")]] <- gm
  for (reg in band_levels) {
    mask <- band_mask(reg)
    globals[[length(globals)+1L]] <- data.frame(
      result_id=sprintf("P0SM%05d",idx_result*10+match(reg,band_levels)), method=mn, variable=vn,
      region=reg, period="1979-2024", n_grid=sum(mask & is.finite(tr)),
      trend=weighted_mean_safe(tr,cell$w,mask), unit=variables[[vn]])
  }
  rm(arr,mat); gc()
}
summary_tab <- do.call(rbind, globals)
write.csv(summary_tab, file.path(OUT,"SURFACE_METHOD_COMPARISON.csv"), row.names=FALSE)

monthly_tab <- data.frame(date=format(dates,"%Y-%m"))
for (nm in names(monthly)) monthly_tab[[nm]] <- monthly[[nm]]
write.csv(monthly_tab,file.path(OUT,"SURFACE_METHOD_GLOBAL_MONTHLY.csv"),row.names=FALSE)

cmp_rows <- list(); rid <- 0L
for (other in c("P1_surface_anchored","P2_one_sided")) for (vn in names(variables)) {
  rid <- rid+1L
  x <- trends[[paste("P0_current_crossing",vn,sep="_")]]
  y <- trends[[paste(other,vn,sep="_")]]
  ok <- is.finite(x)&is.finite(y)
  wx <- cell$w[ok]
  mx <- weighted_mean_safe(x,cell$w,ok); my <- weighted_mean_safe(y,cell$w,ok)
  xw <- x[ok]-weighted.mean(x[ok],wx); yw <- y[ok]-weighted.mean(y[ok],wx)
  rw <- sum(wx*xw*yw)/sqrt(sum(wx*xw^2)*sum(wx*yw^2))
  cmp_rows[[length(cmp_rows)+1L]] <- data.frame(
      result_id=sprintf("P0SC%05d",rid),comparison=paste("P0_current_crossing",other,sep="-vs-"),variable=vn,
    n_grid=sum(ok),P0_trend=mx,comparison_trend=my,difference_P0_minus_comparison=mx-my,
    relative_difference_pct=100*(mx-my)/my,pattern_r=cor(x[ok],y[ok]),pattern_r_coslat=rw,
    sign_agreement_pct=100*mean(sign(x[ok])==sign(y[ok])),unit=variables[[vn]])
}
cmp_tab <- do.call(rbind,cmp_rows)
write.csv(cmp_tab,file.path(OUT,"SURFACE_METHOD_PATTERN_AND_DIFFERENCE.csv"),row.names=FALSE)

# Reconstruct the exact 1-degree surface-height subset used by P1.
zg <- nc_open(file.path(ROOT,"geopotential.nc"))
zlon <- zg$dim$longitude$vals; zlat <- zg$dim$latitude$vals
ix <- vapply(lon,function(x)which.min(abs(zlon-x)),integer(1)); iy <- vapply(lat,function(x)which.min(abs(zlat-x)),integer(1))
zs <- ncvar_get(zg,"z",start=c(1,1,1),count=c(-1,-1,1))[ix,iy]/G
nc_close(zg)
cell$surface_height_m <- as.vector(zs)

terrain_masks <- list(
  all_valid=rep(TRUE,nrow(cell)),
  lowland_le_500m=cell$surface_height_m<=500,
  mid_500_1500m=cell$surface_height_m>500&cell$surface_height_m<=1500,
  high_1500_2500m=cell$surface_height_m>1500&cell$surface_height_m<=2500,
  very_high_gt_2500m=cell$surface_height_m>2500,
  Tibetan_Plateau=cell$lon>=73&cell$lon<=105&cell$lat>=25&cell$lat<=40&cell$surface_height_m>1500,
  Andes=cell$lon>=270&cell$lon<=300&cell$lat>=-55&cell$lat<=15&cell$surface_height_m>1000,
  Rocky_Mountains=cell$lon>=235&cell$lon<=255&cell$lat>=30&cell$lat<=60&cell$surface_height_m>1000,
  East_African_Highlands=cell$lon>=25&cell$lon<=45&cell$lat>=-12&cell$lat<=18&cell$surface_height_m>1000
)
terr_rows<-list(); rid<-0L
for (tm in names(terrain_masks)) for (mn in names(methods)) for (vn in names(variables)) {
  rid<-rid+1L; tr<-trends[[paste(mn,vn,sep="_")]]; mask<-terrain_masks[[tm]]
  terr_rows[[length(terr_rows)+1L]]<-data.frame(result_id=sprintf("P0HT%05d",rid),terrain=tm,method=mn,variable=vn,
    n_grid=sum(mask&is.finite(tr)),trend=weighted_mean_safe(tr,cell$w,mask),unit=variables[[vn]])
}
write.csv(do.call(rbind,terr_rows),file.path(OUT,"HIGH_TERRAIN_VALIDATION_P0_P1_P2.csv"),row.names=FALSE)

# Compact gridded trend fields for audit and future SI plotting.
lon_dim<-ncdim_def("longitude","degrees_east",lon)
lat_dim<-ncdim_def("latitude","degrees_north",lat)
defs<-list()
for(mn in names(methods))for(vn in names(variables))defs[[length(defs)+1L]]<-ncvar_def(paste0(methods[[mn]],"_",vn,"_trend"),variables[[vn]],list(lon_dim,lat_dim),-9999,prec="float")
fo<-nc_create(file.path(OUT,"SURFACE_METHOD_TREND_FIELDS.nc"),defs,force_v4=TRUE)
for(mn in names(methods))for(vn in names(variables))ncvar_put(fo,paste0(methods[[mn]],"_",vn,"_trend"),matrix(trends[[paste(mn,vn,sep="_")]],nrow=length(lon),ncol=length(lat)))
nc_close(fo)

audit <- data.frame(
  item=c("P0_below_ground_values","P1_below_ground_values","P2_below_ground_values","P1_surface_z","P1_surface_q","pressure_top","grid","period"),
  finding=c("used only to pressure-linearly interpolate crossing endpoint","not used","not used","surface geopotential / g","diagnosed from 2 m dewpoint and surface pressure","100 hPa","1 degree 360x181","1979-2024")
)
write.csv(audit,file.path(OUT,"SURFACE_TRUNCATION_AUDIT.csv"),row.names=FALSE)

nc_close(nc0); nc_close(ncx)
