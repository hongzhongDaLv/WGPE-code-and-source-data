$envs = [Environment]::GetEnvironmentVariables("Process")
$pathValue = $envs["Path"]
if ([string]::IsNullOrWhiteSpace($pathValue)) {
    $pathValue = $envs["PATH"]
}
[Environment]::SetEnvironmentVariable("PATH", $null, "Process")
[Environment]::SetEnvironmentVariable("Path", $pathValue, "Process")

$root = "__WGPE_PROJECT_ROOT__"
$script = Join-Path $root "CANONICAL_WGPE_TOP100_v2\07_SCRIPTS\08_wetdry_spatiotemporal_bootstrap.R"
$out = Join-Path $root "CANONICAL_WGPE_TOP100_v2\08_REPRODUCTION_OUTPUT\wetdry_spatiotemporal_uncertainty"
$stdout = Join-Path $out "run_stdout.log"
$stderr = Join-Path $out "run_stderr.log"

& "__WGPE_RSCRIPT__" $script 1> $stdout 2> $stderr
exit $LASTEXITCODE
