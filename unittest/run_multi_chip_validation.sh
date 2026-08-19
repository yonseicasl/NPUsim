#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [[ ! -x "$repo_dir/models/model" ]]; then
    echo "NPUsim model is not built; run ./npusim.sh build npusim first" >&2
    exit 1
fi

# P2-1/P2-2: two small 2-chip synthetic fixtures that need no dataset download, so
# CE6/CE4 multi-chip aggregation and the NoP source-read sharing contract have an
# executable live regression. The two runs are the same hardware and the same total work
# under mirrored mappings (B split vs K split), which is what makes the broadcast /
# partitioned diff attributable to the split dimension alone -- see validation/multi_chip/check.py.
"$repo_dir/npusim.sh" run gemmini_2chip gemm_64x64x64 ws >/dev/null
"$repo_dir/npusim.sh" run gemmini_2chip_ksplit gemm_64x64x64 ws >/dev/null
# L7: the same hardware and mapping with NoP multicast disabled, so the multicast
# link-sharing contract has an A/B that isolates it from every other cost.
"$repo_dir/npusim.sh" run gemmini_2chip_unicast gemm_64x64x64 ws >/dev/null

python3 "$repo_dir/validation/multi_chip/check.py"
