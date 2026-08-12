param(
  [string]$Config = "config/config.yml",
  [switch]$Execute,
  [string]$PythonExe = "python"
)
$ErrorActionPreference = "Stop"
$Repo = (Resolve-Path -LiteralPath $PSScriptRoot).Path
$ConfigPath = if ([IO.Path]::IsPathRooted($Config)) {$Config} else {Join-Path $Repo $Config}
$Validator = Join-Path $Repo "code/support/release_tools/prepare_working_copy.py"
Write-Host "WGPE v1.0.0 full-analysis preflight"
Write-Host "Configuration: $ConfigPath"
& $PythonExe $Validator --config $ConfigPath
$preflight = $LASTEXITCODE
$stages = Import-Csv (Join-Path $Repo "docs/WORKFLOW_STAGES.csv")
Write-Host "Planned stages: $($stages.Count)"
$stages | ForEach-Object { Write-Host (" - {0}: {1}" -f $_.stage_id,$_.description) }
if (-not $Execute) {
  Write-Host "DRY RUN ONLY. No scientific analysis was executed."
  Write-Host "READY_FOR_EXECUTION=$($preflight -eq 0)"
  exit 0
}
if ($preflight -ne 0) { throw "External data/configuration incomplete; no formal output was generated." }
& $PythonExe $Validator --config $ConfigPath --prepare
if ($LASTEXITCODE -ne 0) { throw "Working-copy preparation failed." }
$cfg = @{}
Get-Content -LiteralPath $ConfigPath | ForEach-Object {
  $line = $_.Trim()
  if ($line -and -not $line.StartsWith('#') -and $line.Contains(':')) {
    $parts = $line.Split(':',2); $cfg[$parts[0].Trim()] = $parts[1].Trim().Trim('"').Trim("'")
  }
}
$WorkRoot = $cfg['project_root']
if (-not [IO.Path]::IsPathRooted($WorkRoot)) { $WorkRoot = Join-Path $Repo $WorkRoot }
$Rscript = if ($cfg['rscript_executable']) {$cfg['rscript_executable']} else {'Rscript'}
$Main = Join-Path $WorkRoot 'scripts_TOP100_revision/run_all_TOP100_revision.R'
if (-not (Test-Path -LiteralPath $Main -PathType Leaf)) { throw "Prepared main workflow not found: $Main" }
Write-Host "Starting the registered full workflow. This may require substantial time, memory and storage."
& $Rscript $Main
if ($LASTEXITCODE -ne 0) { throw "Full analysis failed; inspect the workflow log in the configured output root." }
Write-Host "Full analysis completed successfully."
