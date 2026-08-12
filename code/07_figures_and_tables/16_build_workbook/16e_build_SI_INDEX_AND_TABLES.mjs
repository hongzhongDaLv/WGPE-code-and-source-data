import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const ROOT="__WGPE_PROJECT_ROOT__";
const SI=path.join(ROOT,"output_TOP100_revision","10_SUPPLEMENTARY_INFORMATION");
function parseCSV(text){const rows=[];let row=[],cell="",q=false;const clean=v=>v.replace(/[\x00-\x08\x0B\x0C\x0E-\x1F]/g,"");for(let i=0;i<text.length;i++){const c=text[i];if(q){if(c==='"'&&text[i+1]==='"'){cell+='"';i++;}else if(c==='"')q=false;else cell+=c;}else if(c==='"')q=true;else if(c===','){row.push(clean(cell));cell="";}else if(c==='\n'){row.push(clean(cell));rows.push(row);row=[];cell="";}else if(c!=='\r')cell+=c;}if(cell.length||row.length){row.push(clean(cell));rows.push(row);}if(rows.length)rows[0][0]=rows[0][0].replace(/^\uFEFF/,"");return rows;}
async function addCSV(wb,name,file){const rows=parseCSV(await fs.readFile(file,"utf8"));const n=Math.max(...rows.map(r=>r.length));const p=rows.map(r=>[...r,...Array(n-r.length).fill(null)]);const s=wb.worksheets.add(name.slice(0,31));s.getRangeByIndexes(0,0,p.length,n).values=p;s.showGridLines=false;s.freezePanes.freezeRows(1);s.getRangeByIndexes(0,0,1,n).format={fill:"#174A7E",font:{bold:true,color:"#FFFFFF"},wrapText:true};s.getUsedRange().format.autofitColumns();for(let c=0;c<n;c++)s.getRangeByIndexes(0,c,p.length,1).format.columnWidth=20;}
async function one(csv,out,sheet){const wb=Workbook.create();await addCSV(wb,sheet,csv);const x=await SpreadsheetFile.exportXlsx(wb);await x.save(out);}

await one(path.join(SI,"SI_FIGURE_INDEX.csv"),path.join(SI,"SI_FIGURE_INDEX.xlsx"),"SI Figure Index");
await one(path.join(SI,"source_data","SI_PANEL_TO_SOURCE_MAPPING.csv"),path.join(SI,"source_data","SI_PANEL_TO_SOURCE_MAPPING.xlsx"),"Panel Source Mapping");

const wb=Workbook.create();
const files=(await fs.readdir(path.join(SI,"tables"))).filter(x=>x.endsWith(".csv")).sort();
for(const f of files)await addCSV(wb,f.replace(/\.csv$/,""),path.join(SI,"tables",f));
const x=await SpreadsheetFile.exportXlsx(wb);await x.save(path.join(SI,"tables","SI_TABLES_TOP100.xlsx"));
console.log(`Exported SI index, mapping and ${files.length}-sheet SI table workbook`);
