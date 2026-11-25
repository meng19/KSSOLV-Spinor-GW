# RUN_GW  Run a GW calculation from a JSON case file, from any directory.
#
#   .\run_gw.ps1 C:\path\to\case.json
#
# The repository may live anywhere: the launcher adds its own directory to
# the MATLAB path, so no absolute repository path is stored in any file.
# Relative paths inside the JSON are resolved against the JSON file.
#
# Environment variables:
#   GW_JSON          case file, used when no positional argument is given
#   GW_WORK_DIR      working directory (default: the JSON file's directory)
#   MATLAB_BIN       MATLAB executable (default: matlab found on PATH)
#   KSSOLV_ROOT      repository root (default: the directory of this script)
#   GW_MATLAB_ARGS   extra MATLAB options, e.g. "-singleCompThread"

[CmdletBinding()]
param(
    [Parameter(Position = 0)][string]$Json = $env:GW_JSON,
    [string]$RepoRoot = $env:KSSOLV_ROOT,
    [string]$MatlabBin = $env:MATLAB_BIN,
    [string]$WorkDir = $env:GW_WORK_DIR,
    [string[]]$MatlabArgs
)

$ErrorActionPreference = 'Stop'

if (-not $RepoRoot) { $RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path }
if (-not $Json) {
    throw "Usage: run_gw.ps1 <case.json>   (or set GW_JSON)"
}
if (-not (Test-Path -LiteralPath $Json -PathType Leaf)) {
    throw "JSON configuration not found: $Json"
}
$Json = (Resolve-Path -LiteralPath $Json).Path

if (-not $WorkDir) { $WorkDir = Split-Path -Parent $Json }
if (-not (Test-Path -LiteralPath $WorkDir)) {
    New-Item -ItemType Directory -Path $WorkDir | Out-Null
}

if (-not $MatlabBin) {
    $command = Get-Command matlab -ErrorAction SilentlyContinue
    if ($command) { $MatlabBin = $command.Source }
}
if (-not $MatlabBin) {
    throw "MATLAB executable not found; set MATLAB_BIN"
}

if (-not $MatlabArgs -and $env:GW_MATLAB_ARGS) {
    $MatlabArgs = $env:GW_MATLAB_ARGS -split '\s+' | Where-Object { $_ }
}

# Pass the paths through the environment: getenv keeps spaces intact, which a
# hand-built MATLAB string literal would not.
$env:KSSOLV_ROOT = $RepoRoot
$env:GW_JSON = $Json

Write-Host "GW_JSON:     $env:GW_JSON"
Write-Host "KSSOLV_ROOT: $env:KSSOLV_ROOT"
Write-Host "MATLAB_BIN:  $MatlabBin"
Write-Host "WORK DIR:    $WorkDir"

$statement = "addpath(getenv('KSSOLV_ROOT')); gw_run_json(getenv('GW_JSON'));"

Push-Location -LiteralPath $WorkDir
try {
    & $MatlabBin @MatlabArgs -batch $statement
    $exit_code = $LASTEXITCODE
}
finally {
    Pop-Location
}
exit $exit_code