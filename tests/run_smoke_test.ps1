param([string]$PythonExe='python')
& $PythonExe (Join-Path $PSScriptRoot 'smoke_test.py')
exit $LASTEXITCODE
