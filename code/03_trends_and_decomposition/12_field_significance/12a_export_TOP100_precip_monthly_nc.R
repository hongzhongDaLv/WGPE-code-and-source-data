# Purpose: Export the formal monthly precipitation input used by the TOP100
#          downstream workflow to a compact NetCDF for field-significance tests.
# Inputs: ERA5 single-level monthly RDS, 1979-2024.
# Outputs: output_TOP100_revision/derived_data/significance_fields/
#          ERA5_precipitation_monthly_1deg_1979_2024.nc
# Parameters: no transformation; units remain metres per month as stored.
# Dependencies: R >= 4.5, ncdf4.
# Overwrite policy: writes only the TOP100-revision derived-data namespace.

options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
stopifnot(requireNamespace("ncdf4", quietly = TRUE))

ROOT <- "__WGPE_PROJECT_ROOT__"
source(file.path(ROOT,"scripts_TOP100_revision","18_coordinate_harmonization","coordinate_harmonization.R"))
INFILE <- file.path(ROOT, "output_attribution_minimal",
                    "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
OUTDIR <- Sys.getenv("TOP100_CLOSURE_SIGNIFICANCE_DERIVED",file.path(ROOT, "output_TOP100_revision", "derived_data", "significance_fields"))
OUTFILE <- file.path(OUTDIR, "ERA5_precipitation_monthly_1deg_1979_2024.nc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(INFILE)
stopifnot(!is.null(x$tp_native_m), !is.null(x$dates))
cg<-readRDS(file.path(ROOT,"output_TOP100_revision","FINAL_COORDINATE_FIXED_CLOSURE","derived_data","TOP100_CANONICAL_GRID.rds"))
tp <- reorder_exact_grid(x$tp_native_m,x$lon,x$lat,cg$lon,cg$lat)
dates <- as.Date(x$dates)
rm(x); gc()

d <- dim(tp)
stopifnot(length(d) == 3L, d[3] == length(dates))
lon <- cg$lon
lat <- cg$lat
time <- as.numeric(dates - as.Date("1970-01-01"))

dlon <- ncdf4::ncdim_def("lon", "degrees_east", lon)
dlat <- ncdf4::ncdim_def("lat", "degrees_north", lat)
dtime <- ncdf4::ncdim_def("time", "days since 1970-01-01", time, unlim = FALSE)
v <- ncdf4::ncvar_def("precipitation", "m month-1", list(dlon, dlat, dtime),
                      missval = 9.96921e36, prec = "float",
                      longname = "ERA5 total precipitation monthly amount used by TOP100 workflow")
nc <- ncdf4::nc_create(OUTFILE, v, force_v4 = TRUE)
ncdf4::ncvar_put(nc, v, tp)
ncdf4::ncatt_put(nc, 0, "processing", "Exact periodic coordinate reorder from audited single-level RDS to TOP100_CANONICAL_GRID; no interpolation")
ncdf4::ncatt_put(nc, 0, "time_coverage_start", as.character(min(dates)))
ncdf4::ncatt_put(nc, 0, "time_coverage_end", as.character(max(dates)))
ncdf4::nc_close(nc)
message("Wrote: ", OUTFILE)
