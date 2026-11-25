#!/usr/bin/env bash
# RUN_GW  Run a GW calculation from a JSON case file, from any directory.
#
#   ./run_gw.sh /path/to/case.json
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
#   GW_MATLAB_ARGS   extra MATLAB options, e.g. "-singleCompThread -logfile run.log"
#   GW_ENV_FILE      settings file to source (default: ./gw_env.local)
#
# Example (SLURM):
#   export KSSOLV_ROOT=$HOME/opt/KSSOLV-Spinor-GW
#   export GW_JSON=$HOME/cases/si16.json
#   export MATLAB_BIN=/opt/matlab/R2024a/bin/matlab
#   srun "$KSSOLV_ROOT/run_gw.sh"

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Optional site settings (ignored by git), so one install can carry its paths
# without exporting them for every call.
env_file="${GW_ENV_FILE:-$script_dir/gw_env.local}"
if [[ -f "$env_file" ]]; then
    # shellcheck source=/dev/null
    source "$env_file"
fi

repo_root="${KSSOLV_ROOT:-$script_dir}"
json="${1:-${GW_JSON:-}}"

if [[ -z "$json" ]]; then
    echo "usage: $0 <case.json>   (or set GW_JSON)" >&2
    exit 2
fi
if [[ ! -f "$json" ]]; then
    echo "error: JSON configuration not found: $json" >&2
    exit 2
fi
# Make the JSON path absolute before changing the working directory.
json="$(cd -- "$(dirname -- "$json")" && pwd)/$(basename -- "$json")"

work_dir="${GW_WORK_DIR:-$(dirname -- "$json")}"
mkdir -p "$work_dir"

matlab_bin="${MATLAB_BIN:-matlab}"
if ! command -v "$matlab_bin" >/dev/null 2>&1; then
    echo "error: MATLAB executable not found: $matlab_bin (set MATLAB_BIN)" >&2
    exit 2
fi

# Pass the paths through the environment: getenv keeps spaces and quotes in
# file names intact, which a hand-built MATLAB string literal would not.
export KSSOLV_ROOT="$repo_root"
export GW_JSON="$json"

matlab_args=()
if [[ -n "${GW_MATLAB_ARGS:-}" ]]; then
    read -r -a matlab_args <<< "$GW_MATLAB_ARGS"
fi

echo "GW_JSON:     $GW_JSON"
echo "KSSOLV_ROOT: $KSSOLV_ROOT"
echo "MATLAB_BIN:  $matlab_bin"
echo "WORK DIR:    $work_dir"

cd -- "$work_dir"
# -batch must be the last MATLAB option; it also forwards the exit status.
exec "$matlab_bin" ${matlab_args[@]+"${matlab_args[@]}"} -batch \
    "addpath(getenv('KSSOLV_ROOT')); gw_run_json(getenv('GW_JSON'));"