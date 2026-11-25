#!/bin/bash
#SBATCH -J gw
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --gres=gpu:1
#SBATCH -o %j.log
#SBATCH -e %j.err
#SBATCH -t 24:00:00
#
# 自包含作业脚本：改好下面的"用户配置"，然后直接
#     sbatch slurm
# 资源需求用命令行覆盖，不必改文件，例如
#     sbatch -p long -c 1 --gres=gpu:0 slurm
#
# %j.log / %j.err 在提交目录（Slurm 20.02 及以后）或脚本目录（更早版本），
# 需要指定位置时用 sbatch -o / sbatch -e 覆盖。

set -euo pipefail

# ========================= 用户配置（只改这里） =========================
# 算例 JSON。相对路径相对本脚本所在目录；换体系只改这一行。
GW_JSON=$HOME/work/spinor_gw/speed/si16.json

# JSON 里 ${...} 引用所需的环境变量，在这里 export。
# 例如 JSON 中 "qe_path": "${QE_DATA_ROOT}/si16/Si16_FR/02-wfn/Si16_FR.save"
export QE_DATA_ROOT=$HOME/work/spinor_gw/speed

# 仓库根目录；留空表示本脚本所在目录。
KSSOLV_ROOT=

# MATLAB：MATLAB_BIN 优先（绝对路径），否则用 GW_MATLAB_MODULE 加载模块。
MATLAB_BIN=
GW_MATLAB_MODULE=apps/matlab/R2024a

# 额外 MATLAB 选项，例如 -singleCompThread。
GW_MATLAB_ARGS=

# MATLAB 工作目录；留空表示 JSON 所在目录。
GW_WORK_DIR=

# 完成记录目录；留空则不记录。
SLURM_LOG_DIR=$HOME/slurm_log
# =======================================================================

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="${KSSOLV_ROOT:-$script_dir}"

json="${1:-$GW_JSON}"
case "$json" in
    /*) ;;
    *) json="$script_dir/$json" ;;
esac
if [[ ! -f "$json" ]]; then
    echo "case file not found: $json" >&2
    echo "edit GW_JSON at the top of $0" >&2
    exit 2
fi

export GW_JSON="$json"
export MATLAB_BIN GW_MATLAB_MODULE GW_MATLAB_ARGS GW_WORK_DIR

if [[ -n "$GW_MATLAB_MODULE" ]] && command -v module >/dev/null 2>&1; then
    module purge || true
    module load "$GW_MATLAB_MODULE"
fi

# 让 squeue 显示具体算例名，而不是统一的 gw。
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    scontrol update JobID="$SLURM_JOB_ID" \
        JobName="gw:$(basename "$json" .json)" >/dev/null 2>&1 || true
fi

echo "SLURM_JOB_PARTITION=${SLURM_JOB_PARTITION:-}"
echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-}"
echo "SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE:-}"
echo "SUBMIT_DIR=$PWD"
echo "REPO_ROOT=$repo_root"
echo "GW_JSON=$json"
echo "Time is $(date)"

ulimit -s unlimited

# bash 而不是 exec：run_gw.sh 不需要可执行位，且下面的完成记录仍会执行。
status=0
bash "$repo_root/run_gw.sh" "$json" || status=$?

if [[ -n "$SLURM_LOG_DIR" ]]; then
    mkdir -p "$SLURM_LOG_DIR"
    echo "  $PWD  exit=$status  Finished at  $(date)" >> "$SLURM_LOG_DIR/$(date +%Y-%m-%d).log"
fi

exit "$status"