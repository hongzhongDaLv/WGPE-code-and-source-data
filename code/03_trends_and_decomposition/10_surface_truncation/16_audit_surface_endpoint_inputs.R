options(stringsAsFactors = FALSE)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages(library(ncdf4))

root <- "__WGPE_PROJECT_ROOT__"
out <- file.path(root, "output_reviewer_P0_resolution")
files <- c(
  surface_geopotential = file.path(root, "geopotential.nc"),
  dewpoint_2m = file.path(root, "2m dewpoint temperature.nc"),
  surface_pressure = file.path(out, "ERA5_surface_pressure_aligned_0_359.nc")
)

decode_time <- function(nc) {
  candidate <- intersect(c("valid_time", "time"), names(nc$dim))
  if (!length(candidate)) candidate <- intersect(c("valid_time", "time"), names(nc$var))
  if (!length(candidate)) return(c(NA_character_, NA_character_, NA_character_))
  nm <- candidate[1]
  vals <- if (nm %in% names(nc$dim)) nc$dim[[nm]]$vals else ncvar_get(nc, nm)
  units <- ncatt_get(nc, nm, "units")$value
  origin <- sub("^[A-Za-z]+ since ", "", units)
  mult <- if (grepl("hours since", units)) 3600 else if (grepl("days since", units)) 86400 else 1
  tt <- as.POSIXct(origin, tz = "UTC") + vals * mult
  c(format(min(tt), "%Y-%m"), format(max(tt), "%Y-%m"), as.character(length(tt)))
}

rows <- list()
for (role in names(files)) {
  f <- files[[role]]
  nc <- nc_open(f)
  tr <- decode_time(nc)
  for (vn in names(nc$var)) {
    v <- nc$var[[vn]]
    dims <- paste(vapply(v$dim, function(x) paste0(x$name, "=", x$len), character(1)), collapse = ";")
    rows[[length(rows) + 1L]] <- data.frame(
      role = role,
      file = normalizePath(f, winslash = "/", mustWork = TRUE),
      variable = vn,
      units = as.character(ncatt_get(nc, vn, "units")$value),
      dimensions = dims,
      time_start = tr[1], time_end = tr[2], n_time = tr[3],
      fill_value = as.character(ncatt_get(nc, vn, "_FillValue")$value),
      scale_factor = as.character(ncatt_get(nc, vn, "scale_factor")$value),
      add_offset = as.character(ncatt_get(nc, vn, "add_offset")$value)
    )
  }
  nc_close(nc)
}

tab <- do.call(rbind, rows)
write.csv(tab, file.path(out, "SURFACE_ENDPOINT_INPUT_AUDIT.csv"), row.names = FALSE, na = "")

cat("# Surface endpoint input audit\n\n", file = file.path(out, "SURFACE_ENDPOINT_INPUT_AUDIT.md"))
cat("The files were opened with `ncdf4`; entries below are read from their actual dimensions and attributes.\n\n", file = file.path(out, "SURFACE_ENDPOINT_INPUT_AUDIT.md"), append = TRUE)
for (i in seq_len(nrow(tab))) {
  cat(sprintf("- %s / `%s`: %s; %s; time %s to %s (%s).\n", tab$role[i], tab$variable[i], tab$units[i], tab$dimensions[i], tab$time_start[i], tab$time_end[i], tab$n_time[i]), file = file.path(out, "SURFACE_ENDPOINT_INPUT_AUDIT.md"), append = TRUE)
}

