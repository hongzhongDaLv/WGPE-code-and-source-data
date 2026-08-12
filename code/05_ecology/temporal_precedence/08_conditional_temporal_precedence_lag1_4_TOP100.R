# Purpose: Recompute conditional ecological temporal precedence for the formal
#          TOP100 W-GPE definition, in both directions and at lags 1--4.
# Inputs: formal TOP100 ERA5 core, processed Model-3 controls, GPP and LAI panels.
# Outputs: grid-level and area-weighted four-class lag summaries under
#          output_TOP100_revision/results/temporal_precedence/TOP100_conditional.
# Parameters: calendar-month anomaly + linear detrending; BH-FDR q <= 0.05;
#             controls = precipitation + T2m + SSRD + RZSM + CO2 + VPD.
# Dependencies: base R; sufficient memory for the monthly global panels.
# Overwrite policy: writes only TOP100-specific outputs and never modifies legacy files.

options(stringsAsFactors = FALSE, warn = 2)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))

ROOT <- "__WGPE_PROJECT_ROOT__"
V3 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_FINAL_absolute_relative_v3")
OUT <- file.path(ROOT, "output_TOP100_revision", "results", "temporal_precedence", "TOP100_conditional")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "logs", "conditional_temporal_precedence_TOP100.log"), split = TRUE)
on.exit(sink(), add = TRUE)

g <- 9.80665
core_file <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                       "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file <- file.path(ROOT, "output_attribution_minimal", "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
ctrl_file <- file.path(V3, "01_inputs_and_processing", "processed_controls_FINAL_v3.rds")
gpp_file <- file.path(ROOT, "output 1deg v3", "step5_greening", "step5_fixest.RData")
lai_file <- file.path(ROOT, "output 1deg v3", "step6_lai", "step6_fixest.RData")

e <- new.env(parent = emptyenv()); load(core_file, envir = e)
lon <- as.numeric(e$lon); lat <- as.numeric(e$lat); dates <- as.Date(e$all_dates)
wgpe_abs <- e$wgpe_st; z_abs <- e$zbar_st; iwv <- e$iwv_st
rm(e); gc()
single <- readRDS(single_file)
stopifnot(identical(as.character(dates), as.character(as.Date(single$dates))))
ii_lon <- match((lon + 360) %% 360, (as.numeric(single$lon) + 360) %% 360)
ii_lat <- match(lat, as.numeric(single$lat))
if (anyNA(ii_lon) || anyNA(ii_lat)) stop("Cannot align single-level inputs to TOP100 grid")
tp <- single$tp_native_m[ii_lon,ii_lat,,drop=FALSE]
t2m <- single$t2m_K[ii_lon,ii_lat,,drop=FALSE]
ssrd <- single$ssrd_native_J_m2[ii_lon,ii_lat,,drop=FALSE]
rm(single); gc()
ctrl <- readRDS(ctrl_file)
surface_height <- ctrl$surface_height[ii_lon,ii_lat,drop=FALSE]
rzsm <- ctrl$rzsm_arr[ii_lon,ii_lat,,drop=FALSE]
vpd <- ctrl$vpd_arr[ii_lon,ii_lat,,drop=FALSE]
co2 <- as.numeric(ctrl$co2)
rm(ctrl); gc()
z_rel <- sweep(z_abs, c(1,2), surface_height, "-")
wgpe_rel <- g * iwv * z_rel
nlon <- length(lon); nlat <- length(lat)

nearest_to <- function(values, reference) reference[vapply(values, function(v) which.min(abs(reference-v)), integer(1))]
detrend_group <- function(v) {
  ok <- is.finite(v); out <- rep(NA_real_, length(v))
  if (sum(ok) < 10L) return(out)
  tt <- seq_along(v); out[ok] <- residuals(lm(v[ok] ~ tt[ok])); out
}
augment <- function(panel, panel_dates, reference) {
  panel$date <- panel_dates[panel$time]
  panel$month <- as.integer(format(panel$date, "%m"))
  ii <- vapply(panel$lon_val, function(v) which.min(abs(((lon - v + 180) %% 360) - 180)), integer(1))
  jj <- vapply(panel$lat_val, function(v) which.min(abs(lat - v)), integer(1))
  kk <- match(as.character(panel$date), as.character(dates))
  good <- !is.na(ii) & !is.na(jj) & !is.na(kk)
  panel <- panel[good, ]; ii <- ii[good]; jj <- jj[good]; kk <- kk[good]
  lin <- as.numeric(ii) + (as.numeric(jj)-1) * as.numeric(nlon) +
    (as.numeric(kk)-1) * as.numeric(nlon) * as.numeric(nlat)
  w <- if (reference == "absolute") wgpe_abs else wgpe_rel
  panel$wgpe <- as.numeric(w)[lin]
  panel$tp <- as.numeric(tp)[lin]; panel$t2m <- as.numeric(t2m)[lin]
  panel$ssrd <- as.numeric(ssrd)[lin]; panel$rzsm <- as.numeric(rzsm)[lin]
  panel$vpd <- as.numeric(vpd)[lin]; panel$co2 <- co2[kk]
  vars <- c("gpp","wgpe","tp","t2m","ssrd","rzsm","vpd","co2")
  for (vn in vars) {
    panel[[vn]] <- panel[[vn]] - ave(panel[[vn]], panel$grid, panel$month,
      FUN = function(x) mean(x, na.rm = TRUE))
    panel[[vn]] <- ave(panel[[vn]], panel$grid, FUN = detrend_group)
  }
  panel
}

lag_matrix <- function(x, lag) {
  n <- length(x); out <- matrix(NA_real_, n-lag, lag)
  for (k in seq_len(lag)) out[,k] <- x[(lag+1-k):(n-k)]
  out
}
nested_p <- function(target, added, controls, lag) {
  y <- target[(lag+1):length(target)]
  xr <- cbind(1, lag_matrix(target, lag))
  for (cc in controls) xr <- cbind(xr, lag_matrix(cc, lag))
  xf <- cbind(xr, lag_matrix(added, lag))
  ok <- is.finite(y) & apply(is.finite(xf), 1, all)
  if (sum(ok) < ncol(xf)+5L) return(c(p=NA_real_, n=sum(ok), rank_ok=0))
  fr <- lm.fit(xr[ok,,drop=FALSE], y[ok]); ff <- lm.fit(xf[ok,,drop=FALSE], y[ok])
  if (fr$rank < ncol(xr) || ff$rank < ncol(xf)) return(c(p=NA_real_, n=sum(ok), rank_ok=0))
  rssr <- sum(fr$residuals^2); rssf <- sum(ff$residuals^2)
  df1 <- ncol(xf)-ncol(xr); df2 <- sum(ok)-ncol(xf)
  f <- max(0, ((rssr-rssf)/df1)/(rssf/df2))
  c(p=pf(f, df1, df2, lower.tail=FALSE), n=sum(ok), rank_ok=1)
}

run_fourclass <- function(panel, target_name, reference, lag) {
  controls <- c("tp","t2m","ssrd","rzsm","co2","vpd")
  groups <- split(seq_len(nrow(panel)), panel$grid)
  rows <- vector("list", length(groups)); z <- 0L
  for (gg in names(groups)) {
    d <- panel[groups[[gg]], c("gpp","wgpe",controls,"lat_val","lon_val"), drop=FALSE]
    ok <- complete.cases(d)
    if (sum(ok) < 48L) next
    d <- d[ok, ]
    cc <- lapply(controls, function(v) d[[v]])
    fw <- nested_p(d$gpp, d$wgpe, cc, lag)
    rv <- nested_p(d$wgpe, d$gpp, cc, lag)
    if (!is.finite(fw["p"]) || !is.finite(rv["p"])) next
    z <- z+1L
    rows[[z]] <- data.frame(grid_id=as.character(gg), lon=mean(d$lon_val), lat=mean(d$lat_val),
      target=target_name, predictor="W-GPE", height_reference=reference, lag=lag,
      p_forward=fw["p"], p_reverse=rv["p"], n_obs=min(fw["n"],rv["n"]),
      controls="precipitation+T2m+SSRD+RZSM+CO2+VPD", stringsAsFactors=FALSE)
  }
  out <- do.call(rbind, rows[seq_len(z)])
  out$q_forward <- p.adjust(out$p_forward, "BH")
  out$q_reverse <- p.adjust(out$p_reverse, "BH")
  out$class <- ifelse(out$q_forward <= .05 & out$q_reverse <= .05, "bidirectional",
    ifelse(out$q_forward <= .05, "forward-only", ifelse(out$q_reverse <= .05, "reverse-only", "none")))
  out
}

summarize_classes <- function(d) {
  classes <- c("forward-only","reverse-only","bidirectional","none")
  do.call(rbind, lapply(split(d, list(d$height_reference,d$target,d$lag), drop=TRUE), function(x) {
    w <- cos(x$lat*pi/180); den <- sum(w)
    data.frame(height_reference=x$height_reference[1], target=x$target[1], predictor="W-GPE", lag=x$lag[1],
      class=classes, area_fraction_pct=vapply(classes,function(cl)100*sum(w[x$class==cl])/den,numeric(1)),
      n_grid=nrow(x), controls=x$controls[1])
  }))
}

stability_check <- function(panel,target_name,reference) {
  vars <- c("wgpe","tp","t2m","ssrd","rzsm","co2","vpd")
  ok <- complete.cases(panel[,c("gpp",vars)])
  idx <- which(ok); if(length(idx)>200000L) idx <- idx[round(seq(1,length(idx),length.out=200000L))]
  d <- panel[idx,c("gpp",vars)]
  sx <- scale(d[,vars]); cm <- cor(sx)
  cond <- kappa(cbind(1,sx),exact=FALSE)
  vif <- vapply(seq_along(vars),function(j) {
    y<-sx[,j]; xx<-sx[,-j,drop=FALSE]; r2<-summary(lm(y~xx))$r.squared; 1/(1-r2)
  },numeric(1))
  valid_grid <- tapply(ok,panel$grid,function(z)sum(z)>=48L)
  cp <- which(upper.tri(cm),arr.ind=TRUE)
  cmet <- paste0("cor_",rownames(cm)[cp[,1]],"__",colnames(cm)[cp[,2]])
  data.frame(height_reference=reference,target=target_name,
    metric=c(paste0("VIF_",vars),cmet,"condition_number","CO2_residual_variance","VPD_residual_variance","valid_grid_fraction"),
    value=c(vif,cm[cp],cond,var(d$co2),var(d$vpd),mean(valid_grid)),model="Model3",controls="precipitation+T2m+SSRD+RZSM+CO2+VPD")
}

load(gpp_file); gpp_dates <- seq(as.Date("2000-03-01"), by="month", length.out=298)
load(lai_file); lai_dates <- seq(as.Date("2000-01-01"), by="month", length.out=252)
all_rows <- list(); stability_rows <- list()
for (ref in c("absolute","surface_relative")) {
  message("augmenting GPP ",ref); pg <- augment(panel_gpp,gpp_dates,ref)
  message("augmenting LAI ",ref); pl <- augment(panel_lai,lai_dates,ref)
  stability_rows[[length(stability_rows)+1L]] <- stability_check(pg,"GPP",ref)
  stability_rows[[length(stability_rows)+1L]] <- stability_check(pl,"LAI",ref)
  for (lag in 1:4) {
    message("four-class ",ref," lag ",lag," GPP")
    all_rows[[length(all_rows)+1L]] <- run_fourclass(pg,"GPP",ref,lag)
    message("four-class ",ref," lag ",lag," LAI")
    all_rows[[length(all_rows)+1L]] <- run_fourclass(pl,"LAI",ref,lag)
  }
  rm(pg,pl); gc()
}
grid <- do.call(rbind, all_rows); summary <- summarize_classes(grid)
write.csv(grid[grid$height_reference=="absolute",], file.path(OUT,"conditional_temporal_precedence_absolute_TOP100.csv"), row.names=FALSE)
write.csv(grid[grid$height_reference=="surface_relative",], file.path(OUT,"conditional_temporal_precedence_relative_TOP100.csv"), row.names=FALSE)
write.csv(summary, file.path(OUT,"conditional_temporal_precedence_summary_TOP100.csv"), row.names=FALSE)
stability <- do.call(rbind,stability_rows)
write.csv(stability,file.path(OUT,"diagnostics","ecology_model_stability_QA_TOP100.csv"),row.names=FALSE)
co2r <- stability[stability$metric=="CO2_residual_variance",]
writeLines(c("# CO2 residual variance QA TOP100","",paste0("- ",co2r$height_reference," / ",co2r$target,": ",format(co2r$value,scientific=TRUE)),"","CO2 was retained only where the complete Model-3 design matrix was full rank; grid models failing rank or observation-count checks were excluded rather than silently assigned coefficients."),file.path(OUT,"diagnostics","CO2_residual_variance_QA_TOP100.md"))
cat("DONE rows=",nrow(grid),"\n")
