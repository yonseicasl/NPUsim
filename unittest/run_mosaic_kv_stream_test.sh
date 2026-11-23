#!/usr/bin/env bash
# MOSAIC DRAM stream helper checks (needs NPUsim built with MOSAIC_ROOT and that MOSAIC checkout):
#   1. two DRAM configurations in ONE process give the timings of two fresh processes (the helper's memo is
#      keyed by the configuration);
#   2. diagnostics for the attention consumer's K / V read streams, written to OUT (default: stdout only):
#      two fresh half streams vs one stream of both halves, and the long-stream rule against an exact 64 MiB
#      stream. ACCEPTANCE, fixed before the first run: |predicted - exact| / exact <= 1% at 64 MiB.
#      This gates a computational approximation; it says nothing about TPU accuracy.
# Usage: MOSAIC_ROOT=<MOSAIC checkout> bash unittest/run_mosaic_kv_stream_test.sh [OUT_DIR]
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${MOSAIC_ROOT:?set MOSAIC_ROOT to the MOSAIC checkout NPUsim was built against}"
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT
out=${1:-$test_dir/out}
mkdir -p "$out"
[[ ! -e $out/mosaic_kv_stream_diag.csv ]] || { echo "$out already holds a diagnostic; nothing was changed" >&2; exit 1; }

a=$MOSAIC_ROOT/configs/architecture/dram/hbm2e_tpu_v5e_refsb.yaml
b=$repo_dir/validation/tpu_v5e/hbm2e_v5e.yaml
libs="$repo_dir/library:$repo_dir/ext/nebula/library:$repo_dir/ext/DRAMsim3:$MOSAIC_ROOT/simulator/build/lib"
opencv_flags=()
if pkg-config --exists opencv4; then
    read -r -a opencv_flags <<<"$(pkg-config --libs opencv4)"
fi
"${CXX:-g++}" -std=c++11 -O1 -Wall -Wextra -Werror "$repo_dir/unittest/mosaic_kv_stream_test.cc" \
    -L"$repo_dir/library" -L"$repo_dir/ext/nebula/library" -L"$repo_dir/ext/DRAMsim3" \
    -L"$MOSAIC_ROOT/simulator/build/lib" \
    -lnpusim -lnebula -ldramsim3 -lramulator -lopenblas -lpthread -lz "${opencv_flags[@]}" -o "$test_dir/t"
run() { LD_LIBRARY_PATH="$libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$test_dir/t" "$@"; }

bytes=$((1 << 20))
fa=$(run one "$a" $bytes); fb=$(run one "$b" $bytes)
read -r pa pb pa2 <<<"$(run pair "$a" "$b" $bytes)"
echo "1 MiB read: fresh processes $fa / $fb ns; one process $pa / $pb / $pa2 ns (A, B, A again)"
[[ $fa != "$fb" ]] || { echo "the two configurations time alike: the check would be vacuous" >&2; exit 1; }
[[ $pa == "$fa" && $pb == "$fb" && $pa2 == "$fa" ]] || { echo "FAILED: a configuration reused another one's timing" >&2; exit 1; }

run diag "$a" > "$out/mosaic_kv_stream_diag.csv"
python3 - "$out/mosaic_kv_stream_diag.csv" <<'PY'
import sys
rows = [l.strip().split(",") for l in open(sys.argv[1])]
for r in rows:
    if r[0] == "unaligned":
        print(f"unaligned stream of {r[1]} B: {float(r[2]):.3f} ns")
    if r[0] == "two_streams":
        half, th, tt = int(r[1]), float(r[2]), float(r[3])
        print(f"half {half >> 20:3d} MiB: two fresh streams {2 * th:12.3f} ns, one stream of both {tt:12.3f} ns, "
              f"difference {2 * th - tt:+8.3f} ns ({100 * (2 * th / tt - 1):+.3f}%)")
    if r[0] == "extrapolation":
        t16, t32, pred, t64 = map(float, r[1:])
        err = 100 * (pred / t64 - 1)
        print(f"64 MiB: predicted from 16 / 32 MiB {pred:.3f} ns, exact {t64:.3f} ns, {err:+.4f}% (acceptance 1%)")
        assert t32 > t16 > 0 and abs(err) <= 1.0, "long-stream rule outside its acceptance"
PY
echo "MOSAIC KV stream checks passed"
