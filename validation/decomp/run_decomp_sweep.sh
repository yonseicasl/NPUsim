#!/usr/bin/env bash
# Execute the Section D decompression sweep described by runs.csv (see gen_decomp_sweep.py).
# Each run is one `model run-ir <accel.cfg> <network.cfg> <mapping.map>` launched inside its
# own output directory so the per-layer result files never collide. Re-runnable: an existing
# non-empty output directory is skipped unless FORCE=1.
#
#   bash validation/decomp/run_decomp_sweep.sh            # all runs
#   EXP=A bash validation/decomp/run_decomp_sweep.sh      # only experiment A
#   FORCE=1 bash validation/decomp/run_decomp_sweep.sh    # re-run everything
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)
model="$repo/models/model"
manifest="$here/runs.csv"
out_root="$here/out"
exp_filter="${EXP:-}"

[[ -x "$model" ]] || { echo "model not built: run ./npusim.sh build model" >&2; exit 1; }
[[ -f "$manifest" ]] || { echo "manifest missing: run gen_decomp_sweep.py" >&2; exit 1; }

export LD_LIBRARY_PATH="$repo/ext/DRAMsim3:$repo/ext/nebula/library:${LD_LIBRARY_PATH:-}"

n=0; ran=0; skipped=0
# Skip the CSV header, read the columns we need (run_id, experiment, paths).
while IFS=, read -r run_id experiment variant accel_cfg network_cfg mapping workload rest; do
    [[ "$run_id" == "run_id" ]] && continue
    [[ -z "$run_id" ]] && continue
    n=$((n+1))
    [[ -n "$exp_filter" && "$experiment" != "$exp_filter" ]] && continue
    dir="$out_root/$run_id"
    if [[ "${FORCE:-0}" != "1" && -n "$(ls -A "$dir"/*_layer_0.txt 2>/dev/null || true)" ]]; then
        skipped=$((skipped+1)); continue
    fi
    rm -rf "$dir"; mkdir -p "$dir"
    # The network config's dataset list paths are relative to models/ (datasets/gemm/...).
    # Symlink it so the relative path resolves while outputs stay isolated in $dir.
    ln -sfn "$repo/models/datasets" "$dir/datasets"
    ( cd "$dir" && "$model" run-ir "$accel_cfg" "$network_cfg" "$mapping" >run.log 2>&1 ) \
        && ran=$((ran+1)) \
        || echo "FAILED: $run_id (see $dir/run.log)" >&2
done < "$manifest"

echo "runs total=$n executed=$ran skipped(existing)=$skipped filter=${exp_filter:-none}"
echo "outputs under $out_root"
