#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [[ ! -x "$repo_dir/models/model" ]]; then
    echo "NPUsim model is not built; run ./npusim.sh build npusim first" >&2
    exit 1
fi

gemmini_points=(
    gemm_64x64x64
    gemm_128x128x128
    gemm_256x256x256
    gemm_16x512x512
    gemm_512x512x64
    gemm_512x512x512
)

for workload in "${gemmini_points[@]}"; do
    "$repo_dir/npusim.sh" run gemmini "$workload" ws >/dev/null
done
"$repo_dir/npusim.sh" run eyeriss alexnet silicon >/dev/null

python3 "$repo_dir/validation/check_timing.py" --check-baseline
