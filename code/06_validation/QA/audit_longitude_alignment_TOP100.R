options(stringsAsFactors = FALSE, warn = 1)

ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_TOP100_revision", "QA",
                 "LONGITUDE_ALIGNMENT_RECHECK_TOP100")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

core_file <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                       "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
single_file <- file.path(ROOT, "output_attribution_minimal",
                         "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
e <- new.env(parent = emptyenv()); load(core_file, envir = e)
s <- readRDS(single_file)

coord <- rbind(
  data.frame(dataset = "ERA5 TOP100 core", lon_first = e$lon[1], lon_last = tail(e$lon, 1),
             lon_min = min(e$lon), lon_max = max(e$lon), n_lon = length(e$lon),
             lat_first = e$lat[1], lat_last = tail(e$lat, 1), n_lat = length(e$lat),
             convention = "0--359", source_file = core_file),
  data.frame(dataset = "ERA5 single-level bundle", lon_first = s$lon[1], lon_last = tail(s$lon, 1),
             lon_min = min(s$lon), lon_max = max(s$lon), n_lon = length(s$lon),
             lat_first = s$lat[1], lat_last = tail(s$lat, 1), n_lat = length(s$lat),
             convention = "-179--180", source_file = single_file)
)
write.csv(coord, file.path(OUT, "COORDINATE_INVENTORY_TOP100.csv"), row.names = FALSE)

idx <- match((as.numeric(s$lon) + 360) %% 360, (as.numeric(e$lon) + 360) %% 360)
mapping <- data.frame(reference_position = seq_along(s$lon), reference_lon = s$lon,
                      top100_position = idx, top100_lon = e$lon[idx],
                      same_raw_array_position = seq_along(s$lon) == idx)
write.csv(mapping, file.path(OUT, "LONGITUDE_ONE_TO_ONE_MAPPING_TOP100.csv"), row.names = FALSE)

old <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                          "Trend_pattern_correlations_FINAL_v3.csv"), check.names = FALSE)
old <- old[old$height_reference == "absolute", c("variable", "pearson_r")]
names(old)[2] <- "misaligned_estimate"
new <- read.csv(file.path(ROOT, "output_TOP100_revision", "derived_data", "figure_source_data",
                          "Figure1d_spatial_block_bootstrap_LONGITUDE_FIXED_TOP100.csv"),
                check.names = FALSE)
cmp <- merge(new[, c("variable", "estimate", "lower95", "upper95", "longitude_alignment")],
             old, by = "variable", all.x = TRUE)
names(cmp)[names(cmp) == "estimate"] <- "longitude_fixed_estimate"
cmp$difference_fixed_minus_misaligned <- cmp$longitude_fixed_estimate - cmp$misaligned_estimate
write.csv(cmp, file.path(OUT, "FIGURE1D_LONGITUDE_FIXED_COMPARISON.csv"), row.names = FALSE)

fig2 <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "precip_wetdry",
                           "Fig2a_deseasoned_only_TOP100",
                           "Fig2a_WGPE_precip_correlation_summary_TOP100.csv"), check.names = FALSE)
fig2 <- fig2[fig2$region == "Global", ]
write.csv(fig2, file.path(OUT, "FIGURE2A_LONGITUDE_FIXED_GLOBAL.csv"), row.names = FALSE)

status <- data.frame(
  module = c(
    "TOP100 core trends and maps", "TOP100 six native layers", "TOP100 W-GPE product decomposition",
    "JRA-55 comparison using explicit coordinate handling", "ONI W-GPE composite and exact core channels",
    "Fig1d IWV trend-pattern correlation", "Fig1d T2m/precipitation/SSRD trend-pattern correlations",
    "Fig2a W-GPE--precipitation correlation", "Wet/dry composites and decomposition",
    "Precipitation temporal precedence", "Hydroclimate incremental value",
    "Surface-relative-height branches", "Ecology partial correlations",
    "Ecology temporal precedence", "Ecology incremental value", "Ecology field significance"
  ),
  crosses_external_spatial_grid = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
  current_status = c("VALID","VALID","VALID","VALID","VALID","VALID",
                     "CORRECTED_TARGETED","CORRECTED_TARGETED","CORRECTED_RERUN","PENDING_RERUN",
                     "PENDING_RERUN","PENDING_RERUN","PENDING_RERUN","PENDING_RERUN",
                     "PENDING_RERUN","PENDING_RERUN"),
  main_text_use = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,rep(FALSE,7)),
  reason = c(
    "all variables share the TOP100 core grid", "all layer fields share the TOP100 core grid",
    "IWV, z and W-GPE share the TOP100 core grid", "dedicated JRA/ERA coordinate audit",
    "ONI is temporal; W-GPE and channel fields share the core grid", "IWV is in the same TOP100 core",
    "single-level arrays were explicitly reordered before spatial pairing",
    "TOP100 W-GPE was explicitly reordered one-to-one to the precipitation grid",
    "rerun completed after exact one-to-one longitude reconciliation",
    "formal script still pairs TOP100 and precipitation by raw array position",
    "formal script still pairs TOP100 and precipitation by raw array position",
    "surface-height/control grid convention must be reconciled before subtraction",
    "panel-to-core matching uses non-circular longitude and ambiguous half-degree nearest mapping",
    "same ecology coordinate defect plus external controls", "same ecology coordinate defect plus external controls",
    "derived from invalid ecology correlations"
  )
)
write.csv(status, file.path(OUT, "MODULE_VALIDITY_AFTER_LONGITUDE_AUDIT.csv"), row.names = FALSE)

trend <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                            "Trend_results_FINAL_v3.csv"), check.names = FALSE)
trend <- trend[trend$height_reference == "absolute" & trend$region == "Global", ]
dec <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                          "WGPE_decomposition_FINAL_v3.csv"), check.names = FALSE)
dec <- dec[dec$height_reference == "absolute" & dec$region == "Global", ]
wet <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "precip_wetdry",
                          "exact_TOP100", "WET_DRY_EXACT_SUMMARY_TOP100.csv"), check.names = FALSE)
wet <- wet[wet$region == "Global" & wet$metric %in% c("actual", "IWV_channel", "z_channel", "interaction_channel"), ]
enso <- read.csv(file.path(ROOT, "output_TOP100_revision", "derived_data", "figure_source_data",
                           "Figure2_ENSO_standard_band_decomposition_TOP100.csv"), check.names = FALSE)
enso <- enso[enso$region == "Global", ]
f2main <- fig2[fig2$method == "deseasonalized_only", ]
f2dt <- fig2[fig2$method == "deseasonalized_and_detrended", ]

headline <- rbind(
  data.frame(module="core", metric=c("zbar_trend","IWV_trend","WGPE_trend","WGPE_IWV_channel","WGPE_z_channel","WGPE_interaction"),
             value=c(trend$z_trend_m_yr,trend$iwv_trend_kg_m2_yr,trend$wgpe_trend_J_m2_yr,dec$IWV,dec$z,dec$interaction),
             unit=c("m yr^-1","kg m^-2 yr^-1",rep("J m^-2 yr^-1",4)), status="VALID"),
  data.frame(module="Fig1d", metric=paste0("trend_pattern_r_",cmp$variable), value=cmp$longitude_fixed_estimate,
             unit="Pearson r", status="CORRECTED_TARGETED"),
  data.frame(module="Fig2a", metric=c("WGPE_precip_r_deseasonalized_only","WGPE_precip_r_detrended_sensitivity"),
             value=c(f2main$mean_r_area,f2dt$mean_r_area),unit="Pearson r",status="CORRECTED_TARGETED"),
  data.frame(module="wet_dry",metric=paste(wet$event,wet$metric,sep="_"),value=wet$estimate,
             unit="J m^-2",status="CORRECTED_RERUN"),
  data.frame(module="ENSO",metric=c("ENSO_actual","ENSO_IWV_channel","ENSO_z_channel","ENSO_interaction"),
             value=c(enso$IWV_channel_J_m2+enso$z_channel_J_m2+enso$interaction_channel_J_m2,
                     enso$IWV_channel_J_m2,enso$z_channel_J_m2,enso$interaction_channel_J_m2),
             unit="J m^-2",status="VALID")
)
headline$source <- c(rep("results/core",6),rep("Figure1d_spatial_block_bootstrap_LONGITUDE_FIXED_TOP100.csv",nrow(cmp)),
                     rep("Fig2a_WGPE_precip_correlation_summary_TOP100.csv",2),
                     rep("WET_DRY_EXACT_SUMMARY_TOP100.csv",nrow(wet)),
                     rep("Figure2_ENSO_standard_band_decomposition_TOP100.csv",4))
write.csv(headline,file.path(OUT,"LATEST_VALID_HEADLINE_RESULTS_TOP100.csv"),row.names=FALSE)

fmt <- function(x, d = 6) formatC(as.numeric(x), digits = d, format = "f")
lines <- c(
  "# TOP100 经度对齐审计后的最新整合结果",
  "",
  "## 结论",
  "",
  "当前 `PASS WITH INTERPRETIVE CAVEATS` 不能作为全项目科学通过状态。原脚本把 0--359° TOP100 核心与 -179--180° 外部数组按位置直接配对，导致约 180° 经度错位。该错误不影响 TOP100 核心趋势、六层贡献及核心内部精确分解，但影响所有未经坐标匹配的跨数据集分析。",
  "",
  "## 仍然有效的核心结果",
  "",
  paste0("- Global <z> trend: ", fmt(trend$z_trend_m_yr, 6), " m yr^-1."),
  paste0("- Global IWV trend: ", fmt(trend$iwv_trend_kg_m2_yr, 6), " kg m^-2 yr^-1."),
  paste0("- Global W-GPE trend: ", fmt(trend$wgpe_trend_J_m2_yr, 6), " J m^-2 yr^-1."),
  paste0("- W-GPE channels: IWV ", fmt(dec$IWV, 6), "; <z> ", fmt(dec$z, 6),
         "; interaction ", fmt(dec$interaction, 6), " J m^-2 yr^-1."),
  "- Native six-layer vertical decomposition remains valid. The main-text display combines 300--200 and 200--100 hPa as 300--100 hPa; native intervals remain in SI/source data.",
  "",
  "## 已完成的定向修正",
  "",
  paste0("- Fig.1d corrected correlations: IWV ", fmt(cmp$longitude_fixed_estimate[cmp$variable == "IWV"], 4),
         ", T2m ", fmt(cmp$longitude_fixed_estimate[cmp$variable == "T2m"], 4),
         ", precipitation ", fmt(cmp$longitude_fixed_estimate[cmp$variable == "Precipitation"], 4),
         ", SSRD ", fmt(cmp$longitude_fixed_estimate[cmp$variable == "SSRD"], 4), "."),
  paste0("- Fig.2a main (calendar-month climatology removed; no detrending): global area-weighted r = ",
         fmt(f2main$mean_r_area, 6), "."),
  paste0("- Correctly aligned detrended sensitivity: r = ", fmt(f2dt$mean_r_area, 6), "."),
  "- Therefore r = 0.01124 is a longitude-misalignment artefact, not evidence that detrending or IWV control removes the W-GPE--precipitation association.",
  "",
  "## 暂停使用的结果",
  "",
  "Wet/dry exact decomposition has been rerun after coordinate alignment and closes numerically. 降水 temporal precedence、hydroclimate incremental value、surface-relative branches，以及全部 GPP/LAI partial-correlation、temporal-precedence、incremental-value 和 field-significance 结果仍须在统一坐标映射后重算。当前强生态相关不能进入正文。",
  "",
  "## 图形口径",
  "",
  "- Fig.1e main display: 1000--850, 850--700, 700--500, 500--300, 300--100 hPa.",
  "- Fig.1g and Fig.2d display order: Global, N high latitude, N temperate, Tropics, S temperate, S high latitude.",
  "- This audit does not alter the native six-layer scientific source table."
)
writeLines(lines, file.path(OUT, "LATEST_INTEGRATED_RESULTS_LONGITUDE_FIXED_CN.md"), useBytes = TRUE)

qa <- c(
  "# Longitude alignment QA",
  "",
  paste0("- TOP100 longitude convention: ", min(e$lon), " to ", max(e$lon), "."),
  paste0("- Single-level longitude convention: ", min(s$lon), " to ", max(s$lon), "."),
  paste0("- One-to-one longitude mapping: ", !anyNA(idx) && length(unique(idx)) == length(idx), "."),
  paste0("- Raw positions physically identical: ", sum(mapping$same_raw_array_position), "/", nrow(mapping), "."),
  "- Fig1d external correlations recomputed after coordinate reconciliation: YES.",
  "- Fig2a recomputed after coordinate reconciliation: YES.",
  "- Current all-module PASS state retained: NO; corrected modules are marked separately and unrecomputed affected modules remain PENDING_RERUN.",
  "- Core TOP100 trend/decomposition results invalidated: NO."
)
writeLines(qa, file.path(OUT, "LONGITUDE_ALIGNMENT_QA_TOP100.md"), useBytes = TRUE)

rm(e, s); gc()
message("Wrote longitude audit: ", OUT)
