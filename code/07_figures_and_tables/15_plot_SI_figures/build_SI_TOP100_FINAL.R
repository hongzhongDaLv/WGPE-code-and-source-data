options(stringsAsFactors = FALSE)

ROOT <- normalizePath("__WGPE_PROJECT_ROOT__", winslash = "/", mustWork = TRUE)
OUT <- file.path(ROOT, "output_TOP100_revision")
SI <- file.path(OUT, "10_SUPPLEMENTARY_INFORMATION")
DIRS <- file.path(SI, c("figures", "figure_scripts", "source_data", "captions", "tables", "QA"))
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

BLUE <- "#2468A2"; RED <- "#C44E52"; GOLD <- "#B2761B"
GREEN <- "#087567"; GREY <- "#5B6573"; LIGHT <- "#E9ECEF"; BLACK <- "#222222"
rd <- function(...) read.csv(file.path(...), check.names = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
ttl <- function(letter, text) title(main = paste0(letter, "  ", text), adj = 0, font.main = 2, cex.main = 1)
newfig <- function(n, w = 8.2, h = 6.4) {
  assign(".si_name", n, envir = .GlobalEnv); assign(".si_w", w, envir = .GlobalEnv); assign(".si_h", h, envir = .GlobalEnv)
  png(file.path(SI, "figures", paste0("Figure_", n, "_TOP100.png")),
      width = w, height = h, units = "in", res = 600, type = "cairo")
  par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.0), oma = c(0, 0, 0.4, 0),
      las = 1, cex = 0.78, family = "sans")
}
done <- function() {
  rp <- recordPlot(); invisible(dev.off())
  pdf(file.path(SI, "figures", paste0("Figure_", get(".si_name", envir = .GlobalEnv), "_TOP100.pdf")),
      width = get(".si_w", envir = .GlobalEnv), height = get(".si_h", envir = .GlobalEnv), useDingbats = FALSE)
  replayPlot(rp); invisible(dev.off())
}
errplot <- function(est, lo, hi, labs, col = BLUE, ylab = "Estimate", horiz = FALSE) {
  est <- num(est); lo <- num(lo); hi <- num(hi); x <- seq_along(est)
  lim <- range(c(lo, hi, 0), na.rm = TRUE)
  if (!horiz) {
    plot(x, est, ylim = lim, xaxt = "n", xlab = "", ylab = ylab, pch = 19, col = col)
    arrows(x, lo, x, hi, angle = 90, code = 3, length = 0.035, col = BLACK)
    axis(1, x, labs, las = 2, cex.axis = 0.65); abline(h = 0, col = "grey55", lwd = 0.8)
  } else {
    plot(est, x, xlim = lim, yaxt = "n", ylab = "", xlab = ylab, pch = 19, col = col)
    arrows(lo, x, hi, x, angle = 90, code = 3, length = 0.035, col = BLACK)
    axis(2, x, labs, las = 1, cex.axis = 0.65); abline(v = 0, col = "grey55", lwd = 0.8)
  }
}
stack4 <- function(d, main = "", ylab = "Area fraction (%)") {
  cls <- c("forward-only", "reverse-only", "bidirectional", "none")
  cols <- c(BLUE, GOLD, GREEN, "#F2F2F2")
  lags <- sort(unique(num(d$lag)))
  mat <- matrix(0, nrow = length(cls), ncol = length(lags), dimnames = list(cls, lags))
  for (j in seq_along(lags)) {
    z <- d[num(d$lag) == lags[j], ]
    ii <- match(as.character(z$class), cls)
    mat[ii[!is.na(ii)], j] <- num(z$area_fraction_pct[!is.na(ii)])
  }
  mat[is.na(mat)] <- 0
  barplot(mat, col = cols, border = NA, names.arg = lags,
          xlab = "Lag (months)", ylab = ylab, ylim = c(0, 100))
  legend("topright", c("Forward", "Reverse", "Both", "None"), fill = cols, bty = "n", cex = 0.62)
  title(main = main, adj = 0, font.main = 2)
}

P <- list(
  core_abs = file.path(OUT, "results/baseline_FULL_RERUN_TOP100/02_absolute_height_results/Trend_results_FINAL_v3.csv"),
  core_rel = file.path(OUT, "results/baseline_FULL_RERUN_TOP100/03_surface_relative_results/Trend_results_FINAL_v3.csv"),
  top_global = file.path(OUT, "results/top_boundary/TOP_BOUNDARY_GLOBAL_RESULTS.csv"),
  top_diff = file.path(OUT, "results/top_boundary/TOP_BOUNDARY_DIFFERENCES.csv"),
  top_domain = file.path(OUT, "results/top_boundary/TOP_BOUNDARY_DOMAIN_RECONCILIATION.csv"),
  layer = file.path(OUT, "results/top_boundary/TOP_BOUNDARY_LAYER_CONTRIBUTIONS.csv"),
  layer_lat = file.path(OUT, "results/top_boundary/TOP_BOUNDARY_LAYER_LAT5.csv"),
  precip = file.path(OUT, "results/precip_wetdry/Precipitation_correlation_summary_FINAL_v3.csv"),
  fieldsig = file.path(OUT, "results/correlation_significance/FIELD_SIGNIFICANCE_SUMMARY_TOP100.csv"),
  increm = file.path(OUT, "results/incremental_value/INCREMENTAL_VALUE_ROBUSTNESS_TOP100.csv"),
  wetdry = file.path(OUT, "results/precip_wetdry/exact_TOP100/WET_DRY_EXACT_SUMMARY_TOP100.csv"),
  enso_events = file.path(OUT, "results/ENSO/ONI_standard_TOP100/03_absolute_results/ENSO_ONIstandard_absolute_eventwise_values.csv"),
  enso_decomp = file.path(OUT, "results/ENSO/ONI_standard_TOP100/03_absolute_results/ENSO_ONIstandard_absolute_channel_decomposition.csv"),
  enso_ci = file.path(OUT, "results/ENSO/ONI_standard_TOP100/05_bootstrap/ENSO_ONIstandard_clustered_bootstrap_CI.csv"),
  precip_temp = file.path(OUT, "results/temporal_precedence/Precip_temporal_precedence_FINAL_v3.csv"),
  precip_temp_rel = file.path(OUT, "results/baseline_FULL_RERUN_TOP100/03_surface_relative_results/Precip_temporal_precedence_FINAL_v3.csv"),
  ecology = file.path(OUT, "results/ecology/Ecology_partial_correlations_summary_FINAL_v3.csv"),
  eco_ci = file.path(OUT, "derived_data/figure_source_data/Figure3_ecology_summary_uncertainty_TOP100.csv"),
  eco_prof = file.path(OUT, "derived_data/figure_source_data/Figure3_ecology_profile_uncertainty_TOP100.csv"),
  stable = file.path(OUT, "results/ecology/Stable_landcover_sensitivity_FINAL_v3.csv"),
  changed = file.path(OUT, "results/ecology/Changed_landcover_exclusion_FINAL_v3.csv"),
  biome = file.path(OUT, "results/ecology/Biome_summary_FINAL_v3.csv"),
  cond_abs = file.path(OUT, "results/temporal_precedence/TOP100_conditional/conditional_temporal_precedence_absolute_TOP100.csv"),
  cond_rel = file.path(OUT, "results/temporal_precedence/TOP100_conditional/conditional_temporal_precedence_relative_TOP100.csv"),
  cond_sum = file.path(OUT, "results/temporal_precedence/TOP100_conditional/conditional_temporal_precedence_summary_TOP100.csv"),
  vif = file.path(OUT, "results/incremental_value/ecology_TOP100/COLLINEARITY_DIAGNOSTICS_ecology.csv"),
  native = file.path(ROOT, "output_attribution_minimal/highres_z_vertical_profile_QC/native_interval_vertical_profile.csv"),
  agg4 = file.path(ROOT, "output_attribution_minimal/highres_z_vertical_profile_QC/aggregation_to_4layer_check.csv")
)
stopifnot(all(file.exists(unlist(P))))
invisible(file.copy(unlist(P), file.path(SI, "source_data", basename(unlist(P))), overwrite = TRUE))

abs_tr <- read.csv(P$core_abs); rel_tr <- read.csv(P$core_rel)
top <- read.csv(P$top_global); lay <- read.csv(P$layer); wet <- read.csv(P$wetdry)
eco <- read.csv(P$ecology); eco_ci <- read.csv(P$eco_ci); biome <- read.csv(P$biome)
pt <- read.csv(P$precip_temp); ptr <- read.csv(P$precip_temp_rel); cs <- read.csv(P$cond_sum)

# S1: coverage, surface truncation, coordinate definitions and workflow.
newfig("S1")
barplot(c(46, 25, 21), names.arg = c("ERA5", "GPP", "LAI"), horiz = TRUE,
        col = c(BLUE, GREEN, GOLD), border = NA, xlab = "Years represented"); ttl("a", "Data coverage")
plot(NA, xlim = c(0, 1), ylim = c(1050, 80), xaxt = "n", xlab = "", ylab = "Pressure (hPa)")
rect(.18, 1000, .82, 100, col = "#D7E8F5", border = BLUE); polygon(c(.1, .9, .9, .1), c(1050, 760, 1050, 1050), col = "#B7A27A", border = NA)
text(.5, 460, "q and z: valid lower boundary to 100 hPa"); text(.5, 930, "local-surface truncation"); ttl("b", "Pressure integration")
barplot(c(1, .72), names.arg = c("absolute z", "z - surface"), col = c(BLUE, GOLD), border = NA, ylab = "Schematic height"); ttl("c", "Height references")
plot.new(); ttl("d", "Processing workflow"); steps <- c("Monthly pressure-level fields", "Surface-truncated TOP100 column", "Calendar-month anomalies", "Grid-cell diagnostics", "cos(latitude) summaries")
for (i in seq_along(steps)) { y <- .92 - i * .16; rect(.1, y, .9, y + .1, col = LIGHT, border = GREY); text(.5, y + .05, steps[i]); if (i < length(steps)) arrows(.5, y, .5, y - .04, length = .05) }
done()

# S2: absolute versus surface-relative comparison.
newfig("S2")
regions <- c("Global", "S temperate", "Tropics", "N temperate")
for (i in seq_along(regions)) {
  a <- abs_tr[abs_tr$region == regions[i], ]; r <- rel_tr[rel_tr$region == regions[i], ]
  vals <- rbind(c(a$z_trend_m_yr, r$z_trend_m_yr), c(a$wgpe_trend_J_m2_yr, r$wgpe_trend_J_m2_yr))
  bp <- barplot(vals, beside = TRUE, col = c(BLUE, GOLD), border = NA,
                names.arg = c("<z>", "W-GPE"), ylab = "Scaled trend", main = "")
  text(bp, vals, labels = format(round(vals, 3), trim = TRUE), pos = 3, cex = .58)
  ttl(letters[i], regions[i]); if (i == 1) legend("topleft", c("absolute", "surface-relative"), fill = c(BLUE, GOLD), bty = "n", cex = .65)
}
done()

# S3: top-boundary, domain and uncertainty sensitivity.
newfig("S3")
for (i in 1:3) {
  v <- c("iwv", "zbar", "wgpe")[i]; z <- top[top$variable == v, ]; z <- z[order(num(z$top_hPa)), ]
  plot(num(z$top_hPa), num(z$global_trend_deseasonalized), type = "b", pch = 19, col = BLUE,
       xlab = "Top boundary (hPa)", ylab = z$trend_unit[1]); grid(col = "grey90"); ttl(letters[i], paste(v, "trend"))
}
dom <- read.csv(P$top_domain); z <- dom[dom$variable == "wgpe", ]
barplot(num(z$trend), names.arg = paste(z$top_hPa, sub("_.*", "", z$domain), sep = "\n"), col = c(BLUE, GOLD), border = NA, las = 2, ylab = "J m-2 yr-1"); ttl("d", "Averaging-domain distinction")
done()

# S4: six TOP100 layers plus native 1000-300 hPa QC profile.
newfig("S4")
g <- lay[lay$region == "Global", ]; g$layer <- factor(g$layer, levels = rev(c("1000-850", "850-700", "700-500", "500-300", "300-200", "200-100")))
barplot(num(g$total_first_order_contribution_m_yr), names.arg = g$layer, horiz = TRUE,
        col = ifelse(num(g$total_first_order_contribution_m_yr) >= 0, BLUE, GOLD), border = BLACK, xlab = "Contribution (m yr-1)"); abline(v = 0); ttl("a", "Six standard pressure layers")
nat <- read.csv(P$native); nat$density100 <- num(nat$contribution_to_z_trend_m_per_yr) * 100 / num(nat$layer_thickness_hPa)
plot(nat$density100, num(nat$p_mid), type = "b", pch = 16, cex = .45, col = BLUE,
                                ylim = c(1000, 300), xlab = "Density (m yr-1 per 100 hPa)", ylab = "Pressure (hPa)"); abline(v = 0); ttl("b", "Native intervals, 1000-300 hPa QC")
agg <- read.csv(P$agg4); plot(num(agg$reference_4layer_m_per_yr), num(agg$aggregated_highres_m_per_yr), pch = 19, col = BLUE,
                             xlab = "Canonical layer", ylab = "Aggregated native"); abline(0, 1); ttl("c", "Native-to-four-layer closure")
lat <- read.csv(P$layer_lat); q <- aggregate(num(lat$total_first_order_contribution_m_yr), list(lat$layer), mean, na.rm = TRUE)
barplot(q$x, names.arg = q$Group.1, las = 2, col = ifelse(q$x >= 0, BLUE, GOLD), border = NA, ylab = "Mean contribution"); abline(h = 0); ttl("d", "Climate-latitude aggregate")
done()

# S5: precipitation-association robustness and field significance.
newfig("S5")
ps <- read.csv(P$precip); q <- ps[ps$region == "Global", ]; barplot(num(q$mean_r_area), names.arg = q$height_reference, col = c(BLUE, GOLD), border = NA, ylab = "Global Pearson r"); abline(h = 0); ttl("a", "Height-reference sensitivity")
fs <- read.csv(P$fieldsig); z <- fs[fs$test == "common temporal circular shift", ]; errplot(z$observed, z$null_lower95, z$null_upper95, z$metric, BLUE, "Observed and null 95% interval"); ttl("b", "Field significance")
inc <- read.csv(P$increm); z <- inc[inc$target %in% c("precipitation", "wet_event", "dry_event"), ]; errplot(z$global_estimate, z$global_lower95, z$global_upper95, paste(z$target, z$increment, sep = "\n"), BLUE, "Incremental metric"); ttl("c", "Incremental-value robustness")
plot.new(); ttl("d", "Formal testing sequence"); text(.5, .75, "Deseasonalize and detrend"); text(.5, .52, "Effective sample size"); text(.5, .31, "BH-FDR by map"); text(.5, .12, "Common-shift field null")
done()

# S6: exact wet/dry decomposition with spatial-block CIs.
newfig("S6")
for (i in 1:2) {
  ev <- c("Wet", "Dry")[i]; z <- wet[wet$event == ev & wet$region == "Global", ]; z <- z[match(c("actual", "IWV_channel", "z_channel", "interaction_channel"), z$metric), ]
  errplot(z$estimate, z$lower95, z$upper95, c("Total", "IWV", "<z>", "Interaction"), c(GREY, BLUE, GOLD, GREEN), "J m-2"); ttl(letters[i], paste(ev, "global composite"))
}
for (i in 1:2) {
  ev <- c("Wet", "Dry")[i]; z <- wet[wet$event == ev & wet$metric == "actual", ]
  errplot(z$estimate, z$lower95, z$upper95, z$region, if (i == 1) BLUE else GOLD, "J m-2"); ttl(letters[i + 2], paste(ev, "regional totals"))
}
done()

# S7: locked ONI inventory, event variability, exact decomposition and CIs.
newfig("S7")
ee <- read.csv(P$enso_events); eg <- ee[ee$region == "Global", ]; tab <- table(eg$event_type); barplot(tab, col = c(RED, BLUE), border = NA, ylab = "Episodes"); ttl("a", "Locked ONI episode inventory")
boxplot(actual ~ event_type, data = eg, col = c("#E8856A", "#86B6D8"), border = GREY, ylab = "W-GPE anomaly (J m-2)"); abline(h = 0); ttl("b", "Eventwise global values")
ed <- read.csv(P$enso_decomp); z <- ed[ed$region == "Global", ]; barplot(num(z$estimate), names.arg = z$metric, col = c(GREY, BLUE, GOLD, GREEN)[seq_len(nrow(z))], border = NA, las = 2, ylab = "J m-2"); abline(h = 0); ttl("c", "Exact ENSO decomposition")
eci <- read.csv(P$enso_ci); z <- eci[eci$height_reference == "absolute" & eci$region == "Global", ]; errplot(z$estimate_boot, z$lower95, z$upper95, z$metric, BLUE, "J m-2"); ttl("d", "Clustered-episode 95% CIs")
done()

# S8: precipitation temporal-precedence lag sensitivity.
newfig("S8")
for (i in 1:4) {
  ref <- c("absolute", "absolute", "surface_relative", "surface_relative")[i]
  reg <- c("Global", "Tropics", "Global", "Tropics")[i]
  dd <- if (ref == "absolute") pt else ptr
  z <- dd[dd$height_reference == ref & dd$region == reg, ]
  long <- do.call(rbind, lapply(seq_len(nrow(z)), function(k) data.frame(lag = z$lag[k], class = c("forward-only", "reverse-only", "bidirectional", "none"), area_fraction_pct = num(z[k, c("forward_only_valid_area_pct", "reverse_only_valid_area_pct", "bidirectional_valid_area_pct", "none_valid_area_pct")]))))
  stack4(long, paste(letters[i], gsub("_", " ", ref), reg))
}
done()

# S9: ecology Model 0-3 and collinearity diagnostics.
newfig("S9")
for (i in 1:2) {
  tar <- c("GPP", "LAI")[i]; z <- eco[eco$height_reference == "absolute" & eco$region == "Global" & eco$target == tar, ]
  plot(1:4, num(z$mean_partial_r_area[z$predictor == "W-GPE"]), type = "b", pch = 19, col = BLUE, ylim = range(num(z$mean_partial_r_area)), xaxt = "n", xlab = "Control model", ylab = "Partial r")
  lines(1:4, num(z$mean_partial_r_area[z$predictor == "<z>"]), type = "b", pch = 19, col = GOLD); axis(1, 1:4, paste0("M", 0:3)); abline(h = 0, col = "grey60")
  legend("topright", c("W-GPE", "<z>"), col = c(BLUE, GOLD), lty = 1, pch = 19, bty = "n", cex = .65); ttl(letters[i], paste(tar, "Model 0-3"))
}
vd <- read.csv(P$vif); v <- aggregate(cbind(vif_I_Z = num(vd$vif_I_Z), condition_number = num(vd$condition_number)), list(vd$control_set), median, na.rm = TRUE)
barplot(v$vif_I_Z, names.arg = v$Group.1, col = BLUE, border = NA, las = 2, ylab = "Median VIF"); abline(h = 5, col = RED, lty = 2); ttl("c", "VIF by control set")
barplot(v$condition_number, names.arg = v$Group.1, col = GREEN, border = NA, las = 2, ylab = "Median condition number"); ttl("d", "Design-matrix conditioning")
done()

# S10: ecology-domain and land-cover sensitivity.
newfig("S10")
st <- read.csv(P$stable); ch <- read.csv(P$changed); allm <- eco[eco$height_reference == "absolute" & eco$model == "Model3" & eco$region == "Global", c("target", "predictor", "mean_partial_r_area", "n_grid")]
for (i in 1:2) {
  tar <- c("GPP", "LAI")[i]; a <- allm[allm$target == tar, ]; s <- st[st$height_reference == "absolute" & st$model == "Model3" & st$target == tar, ]; x <- ch[ch$height_reference == "absolute" & ch$model == "Model3" & ch$target == tar, ]
  vals <- c(a$mean_partial_r_area[a$predictor == "W-GPE"], s$mean_partial_r_area[s$predictor == "W-GPE"], x$mean_partial_r_area[x$predictor == "W-GPE"])
  barplot(vals, names.arg = c("all persistent", "stable only", "changed excluded"), col = c(BLUE, GREEN, GOLD), border = NA, las = 2, ylab = "Partial r"); abline(h = 0); ttl(letters[i], paste(tar, "W-GPE domain sensitivity"))
}
barplot(c(sum(num(st$n_grid), na.rm = TRUE), sum(num(ch$n_grid), na.rm = TRUE)), names.arg = c("stable", "changed excluded"), col = c(GREEN, GOLD), border = NA, ylab = "Summed grid count"); ttl("c", "Domain sample sizes")
plot.new(); ttl("d", "Land-cover rule"); text(.5, .72, "IGBP code treated as categorical"); text(.5, .48, "Water, urban, snow/ice and barren excluded"); text(.5, .24, "No linear interpolation of classes")
done()

# S11: biome summaries; nonvegetated classes excluded and n < 100 flagged.
newfig("S11", 9.0, 7.0)
nonveg <- grepl("water|urban|snow|ice|barren", biome$class_name, ignore.case = TRUE)
b <- biome[biome$height_reference == "absolute" & biome$model == "Model3" & !nonveg, ]
abbr <- c("Evergreen Needleleaf Forests"="ENF", "Evergreen Broadleaf Forests"="EBF", "Deciduous Needleleaf Forests"="DNF", "Deciduous Broadleaf Forests"="DBF", "Mixed Forests"="MF", "Closed Shrublands"="CSH", "Open Shrublands"="OSH", "Woody Savannas"="WSAV", "Savannas"="SAV", "Grasslands"="GRA", "Permanent Wetlands"="PWL", "Croplands"="CRO", "Cropland/Natural Vegetation Mosaics"="MOS")
comb <- unique(b[, c("target", "predictor")])
for (i in seq_len(min(4, nrow(comb)))) {
  z <- b[b$target == comb$target[i] & b$predictor == comb$predictor[i], ]; z <- z[order(num(z$mean_partial_r_area)), ]
  labs <- unname(abbr[z$class_name]); labs[is.na(labs)] <- z$class_name[is.na(labs)]; labs[num(z$n_grid) < 100] <- paste0(labs[num(z$n_grid) < 100], "*")
  barplot(num(z$mean_partial_r_area), names.arg = labs, horiz = TRUE, las = 1, col = ifelse(num(z$n_grid) < 100, "#CCCCCC", ifelse(num(z$mean_partial_r_area) >= 0, BLUE, GOLD)), border = NA, xlab = "Partial r", cex.names = .68); abline(v = 0); ttl(letters[i], paste(comb$target[i], comb$predictor[i], sep = " | "))
}
done()

# S12: conditional ecological temporal precedence, lag 1-4 and height reference.
newfig("S12")
for (i in 1:4) {
  ref <- c("absolute", "absolute", "surface_relative", "surface_relative")[i]
  d <- cs[cs$height_reference == ref, ]; tar <- c("GPP", "LAI", "GPP", "LAI")[i]
  z <- d[d$target == tar & d$predictor == "W-GPE", ]; stack4(z, paste(letters[i], if (i < 3) "absolute" else "surface-relative", tar))
}
done()

# Machine-readable SI index, panel-to-source mapping and table sources.
index <- data.frame(
  figure = paste0("Figure S", 1:12), panel = "a-d",
  scientific_question = c("coverage, truncation, height definitions and workflow", "absolute versus relative height", "top-boundary and averaging-domain sensitivity", "six-layer and native-interval vertical structure", "precipitation association robustness", "wet/dry exact decomposition", "ONI ENSO inventory and exact decomposition", "precipitation temporal-precedence lag sensitivity", "ecological control hierarchy", "ecological-domain and land-cover sensitivity", "biome-stratified ecological associations", "conditional ecological temporal precedence"),
  source_data = c(P$core_abs, paste(P$core_abs, P$core_rel, sep = "; "), paste(P$top_global, P$top_domain, sep = "; "), paste(P$layer, P$native, P$agg4, sep = "; "), paste(P$precip, P$fieldsig, P$increm, sep = "; "), P$wetdry, paste(P$enso_events, P$enso_decomp, P$enso_ci, sep = "; "), P$precip_temp, paste(P$ecology, P$vif, sep = "; "), paste(P$stable, P$changed, sep = "; "), P$biome, paste(P$cond_abs, P$cond_rel, sep = "; ")),
  source_script = file.path(ROOT, "scripts_TOP100_revision/15_plot_SI_figures/build_SI_TOP100_FINAL.R"),
  height_reference = c("absolute principal", "absolute and relative", rep("absolute principal", 8), "absolute", "absolute and relative"),
  domain = c("global", "global and climate bands", "common/all-valid sensitivity", "global and latitude bands", "global", "global and standard climate bands", "global and standard climate bands", "global and Tropics", "persistent vegetated", "persistent/stable vegetated", "vegetated biomes", "persistent vegetated"),
  period = c("1979-2024; ecology product periods", rep("1979-2024", 4), "1979-2024 events", "locked ONI episodes", "1979-2024", "GPP 2000-2024; LAI 2000-2020", "ecology product periods", "ecology product periods", "ecology product periods"),
  uncertainty_method = c("descriptive", "paired comparison", "time-block/bootstrap summaries", "first-order and aggregation closure", "effective-n, FDR and field null", "10-degree spatial-block bootstrap", "clustered episode bootstrap", "FDR four-class fractions", "spatial-block bootstrap where available", "domain sensitivity", "n_grid threshold; bootstrap source retained", "FDR four-class fractions"),
  main_or_optional = c(rep("main SI", 12)), status = "complete", notes = "TOP100 formal revision; no PDF"
)
write.csv(index, file.path(SI, "SI_FIGURE_INDEX.csv"), row.names = FALSE, fileEncoding = "UTF-8")
mapping <- do.call(rbind, lapply(seq_len(nrow(index)), function(i) data.frame(figure = index$figure[i], panel = strsplit(index$panel[i], "-")[[1]][1], source_data = index$source_data[i], source_script = index$source_script[i], notes = index$notes[i])))
write.csv(mapping, file.path(SI, "source_data", "SI_PANEL_TO_SOURCE_MAPPING.csv"), row.names = FALSE, fileEncoding = "UTF-8")

table_sources <- list(
  Table_S1_datasets = data.frame(product = c("ERA5", "JRA-55", "GPP", "LAI", "Land cover"), role = c("formal TOP100 core", "independent validation", "ecology target", "ecology target", "categorical domain/biome")),
  Table_S2_method_parameters = data.frame(parameter = c("top boundary", "lower boundary", "trend", "area weighting", "FDR"), value = c("100 hPa", "local valid boundary; max pressure level 1000 hPa", "calendar-month anomaly OLS", "cos(latitude)", "Benjamini-Hochberg by map/combination")),
  Table_S3_global_climate_trends = abs_tr,
  Table_S4_exact_decompositions = read.csv(file.path(OUT, "results/core/WGPE_decomposition_FINAL_v3.csv")),
  Table_S5_wet_dry = wet,
  Table_S6_ENSO = read.csv(P$enso_decomp),
  Table_S7_precip_temporal = pt,
  Table_S8_ecology_models = eco,
  Table_S9_ecology_domain_biome = biome,
  Table_S10_conditional_temporal = read.csv(P$cond_sum),
  Table_S11_absolute_relative = merge(abs_tr, rel_tr, by = "region", suffixes = c("_absolute", "_relative")),
  Table_S12_deprecated_final_audit = read.csv(file.path(OUT, "manifests/RESULT_REGISTRY_EXPORT.csv"))
)
for (nm in names(table_sources)) write.csv(table_sources[[nm]], file.path(SI, "tables", paste0(nm, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")

cap_cn <- c("# TOP100 补充图图注", "", "所有图均采用局地有效低层边界至 100 hPa 的绝对位势高度定义；相对地表高度仅作为敏感性。", "",
            paste0("- Figure S", 1:12, ". ", index$scientific_question, "。数据与不确定度口径见 SI_FIGURE_INDEX。"), "", "Figure S11 土地覆盖缩写：ENF 常绿针叶林；EBF 常绿阔叶林；DNF 落叶针叶林；DBF 落叶阔叶林；MF 混交林；CSH 闭合灌丛；OSH 开放灌丛；WSAV 木本稀树草原；SAV 稀树草原；GRA 草地；PWL 永久湿地；CRO 农田；MOS 农田—自然植被镶嵌。星号表示 n_grid < 100，样本不足。")
writeLines(cap_cn, file.path(SI, "captions", "SI_Figure_captions_CN.md"), useBytes = TRUE)
cap_en <- c("# TOP100 Supplementary figure captions", "", "All principal results use the absolute geopotential-height definition from the local valid lower boundary to 100 hPa; surface-relative height is a sensitivity.", "", paste0("- Figure S", 1:12, ". ", index$scientific_question, ". Sources and uncertainty are listed in SI_FIGURE_INDEX."), "", "In Figure S11, abbreviated IGBP biome labels are expanded in the Chinese caption and source table; an asterisk denotes n_grid < 100 (insufficient sample).")
writeLines(cap_en, file.path(SI, "captions", "SI_Figure_captions_EN.md"), useBytes = TRUE)

qa <- c("# SI FIGURE QA", "", paste0("Generated: ", Sys.time()), "", "- Figures S1-S12 generated: YES", "- Definition: local valid lower boundary to 100 hPa; absolute height principal: YES", "- Surface-relative results used only as sensitivity: YES", "- Wet/dry and ENSO use exact decompositions: YES", "- Temporal precedence uses the full four-class definition: YES", "- Ecology uses Model 0-3 and persistent/stable vegetated sensitivity: YES", "- Nonvegetated IGBP classes excluded from biome summaries: YES", "- n_grid < 100 marked insufficient_sample in Figure S11: YES", "- Native pressure-interval curve is explicitly QC for 1000-300 hPa and is not presented as native 20 hPa TOP100 resolution: YES", "- All output figures are 600 dpi PNG: YES", "- Vector PDF companions generated: YES", "- DOCX/manuscript/workbook modified: NO; this is a new TOP100 SI directory")
writeLines(qa, file.path(SI, "QA", "SI_FIGURE_QA.md"), useBytes = TRUE)

cat("Generated", length(list.files(file.path(SI, "figures"), pattern = "png$")), "TOP100 SI figures\n")
