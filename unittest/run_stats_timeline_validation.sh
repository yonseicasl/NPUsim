#!/usr/bin/env bash
set -euo pipefail

# L10: asymmetric multi-chip aggregation regression. Needs the built NPUsim library because
# it drives the real stats_t scaling/merge/finalize chain rather than a pure helper.
repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [[ ! -f "$repo_dir/library/libnpusim.so" ]]; then
    echo "NPUsim library is not built; run ./npusim.sh build npusim first" >&2
    exit 1
fi

test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT

cxx=${CXX:-g++}
opencv_flags=()
if pkg-config --exists opencv4; then
    read -r -a opencv_flags <<<"$(pkg-config --libs opencv4)"
fi

"$cxx" -std=c++11 -O1 -Wall -Wextra -Werror \
    -I"$repo_dir/scheduler" -I"$repo_dir/utils" -I"$repo_dir/components" \
    -I"$repo_dir/ext/nebula/common" -I"$repo_dir/ext/nebula/models/layers" \
    -I"$repo_dir/ext/nebula/models/networks" -I"$repo_dir/ext/DRAMsim3/src" \
    "$repo_dir/unittest/stats_timeline_test.cc" \
    -L"$repo_dir/library" -L"$repo_dir/ext/nebula/library" -L"$repo_dir/ext/DRAMsim3" \
    -lnpusim -lnebula -ldramsim3 -lopenblas -lpthread -lz "${opencv_flags[@]}" \
    -o "$test_dir/stats_timeline_test"

LD_LIBRARY_PATH="$repo_dir/library:$repo_dir/ext/nebula/library:$repo_dir/ext/DRAMsim3${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
    "$test_dir/stats_timeline_test"
