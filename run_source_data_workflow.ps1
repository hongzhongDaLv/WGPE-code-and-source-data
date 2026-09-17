param([string]$PythonExe = "python")
$ErrorActionPreference = "Stop"
$Repo = (Resolve-Path -LiteralPath $PSScriptRoot).Path
& $PythonExe (Join-Path $Repo "tests/smoke_test.py")
if ($LASTEXITCODE -ne 0) { throw "Source-data workflow smoke test failed." }
& $PythonExe (Join-Path $Repo "tests/test_september_update.py")
if ($LASTEXITCODE -ne 0) { throw "September ecology source-data verification failed." }
Write-Host "Source-data verification and test-panel generation completed."
