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

# P4-12: SCALE-Sim v2 cross-simulator gate. The third external reference we hold had never run
# automatically, so a model change could move it unnoticed. SCALE-Sim's own reports are
# committed under validation/phase1/out, so no SCALE-Sim installation is needed here.
"$repo_dir/npusim.sh" run scalesim alexnet matched >/dev/null
python3 "$repo_dir/validation/phase1/gate.py" --check-baseline
