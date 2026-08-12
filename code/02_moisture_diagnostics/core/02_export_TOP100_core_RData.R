# Purpose: Export the formal TOP100 monthly ERA5 core NetCDF to the object names
#          expected by the audited downstream R workflow.
# Inputs:  output_TOP100_revision/derived_data/ERA5_TOP100/ERA5_TOP100_monthly_core_metrics.nc
# Outputs: output_TOP100_revision/derived_data/ERA5_TOP100/global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData
# Parameters: g = 9.80665 m s-2; 1979-01 through 2024-12.
# Dependencies: R >= 4.5, ncdf4.
# Overwrite policy: the output is TOP100-specific and may be replaced only by rerunning this script.

options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Missing package: ncdf4")

ROOT <- "__WGPE_PROJECT_ROOT__"
INFILE <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                    "ERA5_TOP100_monthly_core_metrics.nc")
OUTFILE <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                     "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
if (!file.exists(INFILE)) stop("Missing formal TOP100 NetCDF: ", INFILE)

nc <- ncdf4::nc_open(INFILE)
on.exit(ncdf4::nc_close(nc), add = TRUE)

pick_var <- function(candidates) {
  hit <- candidates[candidates %in% names(nc$var)]
  if (!length(hit)) stop("None of these variables exists: ", paste(candidates, collapse = ", "))
  hit[[1]]
}

lon <- as.numeric(ncdf4::ncvar_get(nc, if ("lon" %in% names(nc$dim)) "lon" else "longitude"))
lat <- as.numeric(ncdf4::ncvar_get(nc, if ("lat" %in% names(nc$dim)) "lat" else "latitude"))
time_name <- if ("time" %in% names(nc$dim)) "time" else stop("No time dimension")
time_raw <- as.numeric(ncdf4::ncvar_get(nc, time_name))
time_units <- nc$dim[[time_name]]$units
origin <- sub(".*since[[:space:]]+", "", time_units)
if (grepl("days since", time_units, ignore.case = TRUE)) {
  all_dates <- as.Date(origin) + time_raw
} else if (grepl("hours since", time_units, ignore.case = TRUE)) {
  all_dates <- as.Date(as.POSIXct(origin, tz = "UTC") + time_raw * 3600)
} else {
  stop("Unsupported time units: ", time_units)
}

iwv_st <- ncdf4::ncvar_get(nc, pick_var(c("IWV", "iwv", "iwv_all", "iwv_st")))
zbar_st <- ncdf4::ncvar_get(nc, pick_var(c("zbar", "z_bar", "zbar_all", "zbar_st")))
wgpe_st <- ncdf4::ncvar_get(nc, pick_var(c("WGPE", "wgpe", "wgpe_all", "wgpe_st")))

expected <- c(length(lon), length(lat), length(all_dates))
for (nm in c("iwv_st", "zbar_st", "wgpe_st")) {
  d <- dim(get(nm))
  if (!identical(as.integer(d), as.integer(expected))) {
    stop(nm, " dimension mismatch: ", paste(d, collapse = "x"),
         " expected ", paste(expected, collapse = "x"))
  }
}

attr(iwv_st, "units") <- "kg m-2"
attr(zbar_st, "units") <- "m"
attr(wgpe_st, "units") <- "J m-2"
attr(wgpe_st, "definition") <- "g * IWV * <z>; local-surface truncated to 100 hPa"
save(wgpe_st, zbar_st, iwv_st, lon, lat, all_dates, file = OUTFILE, compress = "xz")
message("Wrote: ", OUTFILE)
