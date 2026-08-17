#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

cxx=${CXX:-g++}
mapfile -d '' npusim_sources < <(
    find "$repo_dir/scheduler" "$repo_dir/utils" "$repo_dir/components" \
        -maxdepth 1 -type f -name '*.cc' -print0 | sort -z
)

opencv_flags=()
if pkg-config --exists opencv4; then
    read -r -a opencv_flags <<<"$(pkg-config --libs opencv4)"
fi

"$cxx" -std=c++11 -O1 -g -Wall -Wextra \
    -fsanitize=address,undefined -fno-omit-frame-pointer \
    -I"$repo_dir/scheduler" -I"$repo_dir/utils" -I"$repo_dir/components" \
    -I"$repo_dir/ext/nebula/common" -I"$repo_dir/ext/nebula/models/layers" \
    -I"$repo_dir/ext/nebula/models/networks" -I"$repo_dir/ext/DRAMsim3/src" \
    "$repo_dir/models/main.cc" "${npusim_sources[@]}" \
    -L"$repo_dir/ext/nebula/library" -L"$repo_dir/ext/DRAMsim3" \
    -lnebula -ldramsim3 -lopenblas -lpthread -lz "${opencv_flags[@]}" \
    -o "$build_dir/model-sanitized"

export LD_LIBRARY_PATH="$repo_dir/ext/nebula/library:$repo_dir/ext/DRAMsim3${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export NPUSIM_CONFIG_ROOT="$repo_dir/configs"
export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=0:halt_on_error=1}"
export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}"

# Nebula resolves dataset and weight paths relative to the process working directory.
# Symlinks keep all sanitizer-generated result files inside the temporary directory.
ln -s "$repo_dir/models/datasets" "$build_dir/datasets"
ln -s "$repo_dir/models/weights" "$build_dir/weights"

(
    cd "$build_dir"
    ./model-sanitized run gemmini gemm_64x64x64 ws >/dev/null
    ./model-sanitized run eyeriss alexnet silicon >/dev/null
)

echo "Full NPUsim ASan/UBSan smoke tests passed"
