# Purpose: Recompute the locked ONI-standard ENSO composites, exact W-GPE
#          decomposition, uncertainty, and Nino3.4 regressions using TOP100.
# Inputs: formal TOP100 core; locked ONI event realizations; NOAA PSL Nino3.4;
#         precipitation and preprocessed surface height.
# Outputs: output_TOP100_revision/results/ENSO/ONI_standard_TOP100/*
# Parameters: October-March six-month windows; locked CPC ONI event list;
#             equal-winter means; episode-cluster x 10-degree block bootstrap.
# Dependencies: R >= 4.5, ncdf4; optional sandwich for HAC.
# Overwrite policy: only the TOP100 ENSO module directory is replaced.

options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "__WGPE_PROJECT_ROOT__"
SRC <- file.path(ROOT, "WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL",
                 "11_ENSO_ONI_STANDARD_UPDATE", "12_scripts",
                 "03_compute_ONIstandard_composites.R")
BASE <- file.path(ROOT, "output_TOP100_revision", "results", "ENSO",
                  "ONI_standard_TOP100")
CORE <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
EVENT_SRC <- file.path(ROOT, "WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard",
                       "11_ENSO_ONI_STANDARD_UPDATE", "01_event_definition",
                       "ONI_winter_realizations_1979_2024.csv")
NINO_SRC <- file.path(ROOT, "WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard",
                      "11_ENSO_ONI_STANDARD_UPDATE", "00_external_data",
                      "NOAA_PSL_Nino34_monthly_clean.csv")
PATCHED <- file.path(BASE, "12_scripts", "03_compute_ONIstandard_composites_TOP100_EXEC.R")
stopifnot(file.exists(SRC), file.exists(CORE), file.exists(EVENT_SRC), file.exists(NINO_SRC))
dir.create(dirname(PATCHED), recursive = TRUE, showWarnings = FALSE)
x <- readLines(SRC, warn = FALSE, encoding = "UTF-8")

set_line <- function(pattern, replacement, fixed = TRUE) {
  i <- grep(pattern, x, fixed = fixed)
  if (length(i) != 1L) stop("Expected one match for ", pattern, "; got ", length(i))
  x[i] <<- replacement
}
set_line('ROOT<-"__WGPE_PROJECT_ROOT__";BASE<-file.path(ROOT,"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard","11_ENSO_ONI_STANDARD_UPDATE")',
         sprintf('ROOT<-"%s";BASE<-"%s"', ROOT, gsub("\\\\", "/", BASE)))
set_line('core_file<-file.path(ROOT,"output_attribution_minimal","global_wgpe_1deg_surface_truncated_300hPa_FIXED_FULL.RData");single_file<-file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds");control_file<-file.path(ROOT,"output_attribution_minimal","FULL_RERUN_FINAL_absolute_relative_v3","01_inputs_and_processing","processed_controls_FINAL_v3.rds")',
         sprintf('core_file<-"%s";single_file<-file.path(ROOT,"output_attribution_minimal","ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds");control_file<-file.path(ROOT,"output_attribution_minimal","FULL_RERUN_FINAL_absolute_relative_v3","01_inputs_and_processing","processed_controls_FINAL_v3.rds")', gsub("\\\\", "/", CORE)))
set_line('real<-read.csv(file.path(BASE,"01_event_definition","ONI_winter_realizations_1979_2024.csv"));real<-real[as.logical(real$included_main),];real$winter_start<-as.Date(real$winter_start);real$winter_end<-as.Date(real$winter_end);nr<-nrow(real);ridx<-lapply(seq_len(nr),function(i)which(dates>=real$winter_start[i]&dates<=real$winter_end[i]));stopifnot(all(lengths(ridx)==6));phase<-real$event_type;phases<-c("ElNino","LaNina")',
         sprintf('real<-read.csv("%s");real<-real[as.logical(real$included_main),];real$winter_start<-as.Date(real$winter_start);real$winter_end<-as.Date(real$winter_end);nr<-nrow(real);ridx<-lapply(seq_len(nr),function(i)which(dates>=real$winter_start[i]&dates<=real$winter_end[i]));stopifnot(all(lengths(ridx)==6));phase<-real$event_type;phases<-c("ElNino","LaNina")', gsub("\\\\", "/", EVENT_SRC)))
set_line('nino<-read.csv(file.path(BASE,"00_external_data","NOAA_PSL_Nino34_monthly_clean.csv"));stopifnot(all(as.Date(nino$date)==dates));xni<-nino$nino34_anom_detrended',
         sprintf('nino<-read.csv("%s");stopifnot(all(as.Date(nino$date)==dates));xni<-nino$nino34_anom_detrended', gsub("\\\\", "/", NINO_SRC)))

writeLines(c('# AUTO-PATCHED TOP100 ONI-STANDARD EXECUTION COPY', x), PATCHED, useBytes = TRUE)
status <- system2("__WGPE_RSCRIPT__", PATCHED)
if (!identical(status, 0L)) stop("TOP100 ONI-standard module failed: ", status)
