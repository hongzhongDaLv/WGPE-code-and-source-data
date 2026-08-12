import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const ROOT = "__WGPE_PROJECT_ROOT__";
const OUT = path.join(ROOT, "output_TOP100_revision");
const PREVIEW = path.join(OUT, "QA", "workbook_previews");
await fs.mkdir(PREVIEW, { recursive: true });

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:F1").merge();
readme.getRange("A1").values = [["WGPE TOP100 — results and methods detail"]];
readme.getRange("A1:F1").format = {
  fill: "#174A7E", font: { bold: true, color: "#FFFFFF", size: 16 },
  verticalAlignment: "center"
};
readme.getRange("A1:F1").format.rowHeight = 30;
readme.getRange("A3:B11").values = [
  ["Formal definition", "Local valid lower boundary to 100 hPa; absolute geopotential-height coordinate"],
  ["ERA5 period", "1979-01 to 2024-12"],
  ["Primary grid", "1 degree global longitude-latitude"],
  ["Area weighting", "cos(latitude)"],
  ["Trend method", "Calendar-month anomalies; grid-cell OLS; area-weighted summaries"],
  ["Main result status", "Only main_approved records are eligible for manuscript text"],
  ["Registry row count", null],
  ["Registry ID uniqueness", "Verified in RESULT_ID_AUDIT.csv"],
  ["Generated", new Date()]
];
readme.getRange("B9").formulas = [["=COUNTA('Result Registry'!A2:A20000)"]];
readme.getRange("A3:A11").format = { fill: "#DCE6F1", font: { bold: true, color: "#1F1F1F" } };
readme.getRange("A3:B11").format.borders = { preset: "inside", style: "thin", color: "#D9E1E8" };
readme.getRange("A3:A11").format.columnWidth = 24;
readme.getRange("B3:B11").format.columnWidth = 78;
readme.getRange("B3:B11").format.wrapText = true;
readme.getRange("B11").format.numberFormat = "yyyy-mm-dd hh:mm";

const csvSheets = [
  ["Result Registry", "output_TOP100_revision/manifests/RESULT_REGISTRY_EXPORT.csv"],
  ["Writing Values", "output_TOP100_revision/derived_data/figure_source_data/TOP100_WRITING_SOURCE_VALUES.csv"],
  ["Core Trends", "output_TOP100_revision/results/core/Trend_results_FINAL_v3.csv"],
  ["WGPE Decomposition", "output_TOP100_revision/results/core/WGPE_decomposition_FINAL_v3.csv"],
  ["Vertical Layers", "output_TOP100_revision/results/top_boundary/TOP_BOUNDARY_LAYER_CONTRIBUTIONS.csv"],
  ["Top Sensitivity", "output_TOP100_revision/results/top_boundary/TOP_BOUNDARY_DIFFERENCES.csv"],
  ["Precip Association", "output_TOP100_revision/results/precip_wetdry/Precipitation_correlation_summary_FINAL_v3.csv"],
  ["Wet Dry", "output_TOP100_revision/results/precip_wetdry/exact_TOP100/WET_DRY_EXACT_SUMMARY_TOP100.csv"],
  ["ENSO Decomposition", "output_TOP100_revision/results/ENSO/ONI_standard_TOP100/03_absolute_results/ENSO_ONIstandard_absolute_channel_decomposition.csv"],
  ["JRA55 Validation", "output_TOP100_revision/results/JRA55/JRA55_TOP100_GLOBAL_TRENDS.csv"],
  ["JRA55 Patterns", "output_TOP100_revision/results/JRA55/ERA5_JRA55_TOP100_PATTERN_METRICS.csv"],
  ["Ecology Summary", "output_TOP100_revision/results/ecology/Ecology_partial_correlations_summary_FINAL_v3.csv"],
  ["Temporal Summary", "output_TOP100_revision/results/temporal_precedence/TOP100_conditional/conditional_temporal_precedence_summary_TOP100.csv"],
  ["Temporal Effects", "output_TOP100_revision/results/temporal_precedence/TOP100_effect_sizes/ECOLOGY_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv"],
  ["Incremental Value", "output_TOP100_revision/results/incremental_value/INCREMENTAL_VALUE_ROBUSTNESS_TOP100.csv"],
  ["Field Significance", "output_TOP100_revision/results/correlation_significance/FIELD_SIGNIFICANCE_SUMMARY_TOP100.csv"],
  ["Figure Source", "output_TOP100_revision/derived_data/figure_source_data/TOP100_WRITING_SOURCE_VALUES.csv"],
  ["FILE_INDEX", "output_TOP100_revision/manifests/FILE_INDEX.csv"],
  ["SCRIPT_INDEX", "output_TOP100_revision/manifests/SCRIPT_MANIFEST.csv"],
  ["FIGURE_INDEX", "output_TOP100_revision/manifests/FIGURE_MANIFEST.csv"],
  ["METHOD_FILE_INDEX", "output_TOP100_revision/manifests/METHOD_FILE_INDEX.csv"],
  ["DERIVED_DATA_INDEX", "output_TOP100_revision/manifests/DERIVED_DATA_INDEX.csv"],
  ["SI_FIGURE_INDEX", "output_TOP100_revision/10_SUPPLEMENTARY_INFORMATION/SI_FIGURE_INDEX.csv"],
  ["SI_PANEL_SOURCE", "output_TOP100_revision/10_SUPPLEMENTARY_INFORMATION/source_data/SI_PANEL_TO_SOURCE_MAPPING.csv"]
];

// Parse CSV locally and write values into newly created worksheets. The
// artifact-tool instance method `fromCSV()` hydrates a complete collaborative
// document and cannot be called repeatedly on a non-empty workbook.
function parseCSV(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        cell += '"';
        i++;
      } else if (ch === '"') {
        quoted = false;
      } else {
        cell += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ',') {
      row.push(cell);
      cell = "";
    } else if (ch === '\n') {
      row.push(cell);
      rows.push(row);
      row = [];
      cell = "";
    } else if (ch !== '\r') {
      cell += ch;
    }
  }
  if (cell.length || row.length) {
    row.push(cell);
    rows.push(row);
  }
  if (rows.length && rows[0].length) rows[0][0] = rows[0][0].replace(/^\uFEFF/, "");
  const numeric = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
  return rows.map((r, ri) => r.map((v) => {
    if (ri === 0) return v;
    const t = v.trim();
    if (t === "") return null;
    if (/^(TRUE|FALSE)$/i.test(t)) return /^TRUE$/i.test(t);
    if (numeric.test(t) && !/^0\d+/.test(t)) return Number(t);
    return v;
  }));
}

for (const [sheetName, rel] of csvSheets) {
  const full = path.join(ROOT, rel);
  const csv = await fs.readFile(full, "utf8");
  const rows = parseCSV(csv);
  const colCount = Math.max(...rows.map((r) => r.length));
  const padded = rows.map((r) => [...r, ...Array(colCount - r.length).fill(null)]);
  const sheet = workbook.worksheets.add(sheetName);
  sheet.getRangeByIndexes(0, 0, padded.length, colCount).values = padded;
}

const params = JSON.parse(await fs.readFile(path.join(OUT, "methods", "METHOD_PARAMETERS.json"), "utf8"));
const flatten = (obj, prefix = "") => {
  const rows = [];
  for (const [k, v] of Object.entries(obj)) {
    const key = prefix ? `${prefix}.${k}` : k;
    if (v && typeof v === "object" && !Array.isArray(v)) rows.push(...flatten(v, key));
    else rows.push([key, Array.isArray(v) ? v.join("; ") : v]);
  }
  return rows;
};
const paramSheet = workbook.worksheets.add("Method Parameters");
const paramRows = [["parameter", "value"], ...flatten(params)];
paramSheet.getRangeByIndexes(0, 0, paramRows.length, 2).values = paramRows;

const products = workbook.worksheets.add("Data Products");
const productRows = [
  ["product", "path", "role"],
  ["ERA5 TOP100 monthly core", "output_TOP100_revision/derived_data/ERA5_TOP100/global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData", "formal atmospheric core"],
  ["Top-boundary trend fields", "output_TOP100_revision/derived_data/trend_fields/TOP_BOUNDARY_GRIDCELL_FIELDS.nc", "TOP100/TOP200/TOP300 sensitivity"],
  ["Six-layer grid fields", "output_TOP100_revision/derived_data/decomposition_fields/TOP100_LAYER_GRIDCELL_FIELDS.nc", "vertical contribution diagnostic"],
  ["Figure 1", "output_TOP100_revision/figures/main/Figure1_TOP100_FINAL/Figure1_TOP100_FINAL.png", "main figure"],
  ["Figure 2", "output_TOP100_revision/figures/main/Figure2_TOP100_FINAL/Figure2_TOP100_FINAL.png", "main figure"],
  ["Figure 3", "output_TOP100_revision/figures/main/Figure3_TOP100_FINAL/Figure3_TOP100_FINAL.png", "main figure"],
  ["Methods", "output_TOP100_revision/methods/METHODS_FULL_TOP100_CN.md", "formal methods text"],
  ["Registry", "output_TOP100_revision/manifests/RESULT_REGISTRY_EXPORT.csv", "complete result registry"]
];
products.getRangeByIndexes(0, 0, productRows.length, 3).values = productRows;

const sheetInfo = await workbook.inspect({ kind: "sheet", include: "id,name", maxChars: 6000 });
await fs.writeFile(path.join(OUT, "QA", "workbook_sheet_inventory.ndjson"), sheetInfo.ndjson ?? String(sheetInfo), "utf8");

for (let i = 0; i < workbook.worksheets.items.length; i++) {
  const sheet = workbook.worksheets.getItemAt(i);
  sheet.showGridLines = false;
  const used = sheet.getUsedRange();
  if (!used) continue;
  const rowCount = used.rowCount;
  const colCount = used.columnCount;
  if (sheet.name !== "README") {
    sheet.getRangeByIndexes(0, 0, 1, colCount).format = {
      fill: "#174A7E", font: { bold: true, color: "#FFFFFF" },
      wrapText: true, verticalAlignment: "center"
    };
    sheet.getRangeByIndexes(0, 0, 1, colCount).format.rowHeight = 28;
    sheet.freezePanes.freezeRows(1);
    used.format.autofitColumns();
    const headers = sheet.getRangeByIndexes(0, 0, 1, colCount).values[0];
    for (let c = 0; c < colCount; c++) {
      const h = String(headers[c] ?? "").toLowerCase();
      let width = 16;
      if (/source|path|notes|method|uncertainty|controls|description/.test(h)) width = 38;
      else if (/result_id|registry_id|variable|predictor|target|region/.test(h)) width = 20;
      else if (/estimate|lower|upper|trend|fraction|percent|mean|r$/.test(h)) width = 15;
      sheet.getRangeByIndexes(0, c, rowCount, 1).format.columnWidth = width;
    }
    if (rowCount > 1) sheet.getRangeByIndexes(1, 0, rowCount - 1, colCount).format.rowHeight = 18;
  }
}

const registryCheck = await workbook.inspect({
  kind: "table", range: "'Result Registry'!A1:H12", include: "values,formulas",
  tableMaxRows: 12, tableMaxCols: 8, maxChars: 5000
});
await fs.writeFile(path.join(OUT, "QA", "workbook_registry_inspect.ndjson"), registryCheck.ndjson ?? String(registryCheck), "utf8");
const errors = await workbook.inspect({
  kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 }, summary: "final formula error scan"
});
await fs.writeFile(path.join(OUT, "QA", "workbook_formula_error_scan.ndjson"), errors.ndjson ?? String(errors), "utf8");

for (let i = 0; i < workbook.worksheets.items.length; i++) {
  const sheet = workbook.worksheets.getItemAt(i);
  const preview = await workbook.render({ sheetName: sheet.name, range: "A1:H25", scale: 0.8, format: "png" });
  const safe = sheet.name.replace(/[^A-Za-z0-9_-]/g, "_");
  await fs.writeFile(path.join(PREVIEW, `${String(i + 1).padStart(2, "0")}_${safe}.png`),
                     new Uint8Array(await preview.arrayBuffer()));
}

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(path.join(OUT, "ALL_RESULTS_AND_METHODS_DETAIL.xlsx"));
console.log(`Exported ${path.join(OUT, "ALL_RESULTS_AND_METHODS_DETAIL.xlsx")}`);
