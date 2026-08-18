#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [[ ! -x "$repo_dir/models/model" ]]; then
    echo "NPUsim model is not built; run ./npusim.sh build npusim first" >&2
    exit 1
fi

# P2-1/P2-2: a small 2-chip synthetic fixture that needs no dataset download, so
# CE6/CE4 multi-chip aggregation has an executable live regression.
"$repo_dir/npusim.sh" run gemmini_2chip gemm_64x64x64 ws >/dev/null

python3 "$repo_dir/validation/multi_chip/check.py"
