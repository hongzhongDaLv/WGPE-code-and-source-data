import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const ROOT = "__WGPE_PROJECT_ROOT__";
const OUT = path.join(ROOT, "output_TOP100_revision");
const workbook = Workbook.create();

function parseCSV(text) {
  const rows = []; let row = [], cell = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') { cell += '"'; i++; }
      else if (ch === '"') quoted = false; else cell += ch;
    } else if (ch === '"') quoted = true;
    else if (ch === ',') { row.push(cell); cell = ""; }
    else if (ch === '\n') { row.push(cell); rows.push(row); row = []; cell = ""; }
    else if (ch !== '\r') cell += ch;
  }
  if (cell.length || row.length) { row.push(cell); rows.push(row); }
  if (rows.length) rows[0][0] = rows[0][0].replace(/^\uFEFF/, "");
  return rows;
}

async function addCSV(name, rel) {
  const rows = parseCSV(await fs.readFile(path.join(ROOT, rel), "utf8"));
  const ncol = Math.max(...rows.map(r => r.length));
  const padded = rows.map(r => [...r, ...Array(ncol-r.length).fill(null)]);
  const s = workbook.worksheets.add(name);
  s.getRangeByIndexes(0,0,padded.length,ncol).values = padded;
  s.showGridLines = false; s.freezePanes.freezeRows(1);
  s.getRangeByIndexes(0,0,1,ncol).format = {fill:"#174A7E",font:{bold:true,color:"#FFFFFF"},wrapText:true};
  const used=s.getUsedRange(); used.format.autofitColumns();
  for(let c=0;c<ncol;c++) s.getRangeByIndexes(0,c,padded.length,1).format.columnWidth=Math.min(42,Math.max(14,s.getRangeByIndexes(0,c,1,1).values[0][0]?.length||14));
}

const readme = workbook.worksheets.add("README"); readme.showGridLines=false;
readme.getRange("A1:F1").merge(); readme.getRange("A1").values=[["WGPE TOP100 data products and processing"]];
readme.getRange("A1:F1").format={fill:"#174A7E",font:{bold:true,color:"#FFFFFF",size:16}};
readme.getRange("A3:B10").values=[
  ["Formal definition","Local valid lower boundary to 100 hPa; absolute geopotential height"],
  ["Convergence sensitivity","TOP200"],["Truncation sensitivity","TOP300"],
  ["Primary period","ERA5 1979-2024"],["JRA-55 validation","1982-2023"],
  ["Grid","1 degree ERA5; 1.25 degree JRA-55 native and harmonized comparisons"],
  ["Missing values","Below-ground levels excluded; crossing intervals surface-truncated"],
  ["Traceability","Paths and hashes are recorded in FILE_INDEX and CHECKSUMS_SHA256"]];
readme.getRange("A3:A10").format={fill:"#DCE6F1",font:{bold:true}}; readme.getRange("A3:A10").format.columnWidth=24; readme.getRange("B3:B10").format.columnWidth=78; readme.getRange("B3:B10").format.wrapText=true;

await addCSV("Input Inventory","output_TOP100_revision/temp/INPUT_INVENTORY.csv");
await addCSV("Top Boundary Audit","output_TOP100_revision/results/top_boundary/TOP_BOUNDARY_INPUT_AUDIT.csv");
await addCSV("Surface Endpoint Audit","output_TOP100_revision/results/surface_truncation/SURFACE_ENDPOINT_INPUT_AUDIT.csv");
await addCSV("Derived Data Index","output_TOP100_revision/manifests/DERIVED_DATA_INDEX.csv");
await addCSV("File Index","output_TOP100_revision/manifests/FILE_INDEX.csv");

const params = JSON.parse(await fs.readFile(path.join(OUT,"methods","METHOD_PARAMETERS.json"),"utf8"));
const flat=[]; const walk=(o,p="")=>{for(const [k,v] of Object.entries(o)){const q=p?`${p}.${k}`:k;if(v&&typeof v==="object"&&!Array.isArray(v))walk(v,q);else flat.push([q,Array.isArray(v)?v.join("; "):String(v)]);}}; walk(params);
const ps=workbook.worksheets.add("Method Parameters"); ps.getRangeByIndexes(0,0,flat.length+1,2).values=[["parameter","value"],...flat]; ps.getRange("A1:B1").format={fill:"#174A7E",font:{bold:true,color:"#FFFFFF"}}; ps.getRange("A:A").format.columnWidth=36; ps.getRange("B:B").format.columnWidth=70;

const xlsx=await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(path.join(OUT,"methods","DATA_PRODUCTS_AND_PROCESSING.xlsx"));
console.log("Exported DATA_PRODUCTS_AND_PROCESSING.xlsx");

