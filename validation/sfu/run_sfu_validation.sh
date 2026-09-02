#!/usr/bin/env bash
# SFU validation driver (plan/plan_sfu.md test strategy).
#
# Builds validation/sfu/sfu_test.cc against the SFU component and runs:
#   * formula/boundary checks (bypass, lane tail, latency/II, chunking)
#   * event-count identities (elements = B x K x P x Q, chip-partition sum)
#   * energy linearity and per-operation profile independence
#   * softmax multi-pass hand calculation
#   * functional evaluators vs reference math
#   * fail-fast rejection of invalid [sfu] configs (subprocess exit codes)
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT

cxx=${CXX:-g++}
sanitizers=${NPUSIM_TEST_SANITIZERS:-}
# Phase-6 item 3: ext/nebula/common/activations.cc is linked in so the evaluators are
# compared against the ACTUAL Nebula implementations, not re-derived formulas.
cxxflags=(-std=c++11 -Wall -Wextra -Werror
          -I"$repo_dir/utils" -I"$repo_dir/components" -I"$repo_dir/ext/nebula/common")
if [[ $sanitizers == 1 ]]; then
    cxxflags+=(-fsanitize=address,undefined -fno-omit-frame-pointer)
fi

# Nebula is third-party style code (unused parameters in no-op gradients); compile it
# without -Werror so ITS warnings cannot mask a warning-free NPUsim build.
"$cxx" -std=c++11 -Wall -Wno-unused-parameter \
    ${sanitizers:+-fsanitize=address,undefined -fno-omit-frame-pointer} \
    -c "$repo_dir/ext/nebula/common/activations.cc" -o "$test_dir/nebula_activations.o"

"$cxx" "${cxxflags[@]}" \
    "$repo_dir/validation/sfu/sfu_test.cc" \
    "$repo_dir/components/sfu.cc" \
    "$repo_dir/utils/config.cc" \
    "$repo_dir/utils/utils.cc" \
    "$repo_dir/utils/datatype.cc" \
    "$repo_dir/utils/energy_units.cc" \
    "$test_dir/nebula_activations.o" \
    -o "$test_dir/sfu_test"

"$test_dir/sfu_test"

expect_failure() {
    if "$test_dir/sfu_test" "$1" >/dev/null 2>&1; then
        echo "Expected SFU fail-fast for mode $1" >&2
        exit 1
    fi
}

# Plan: 0 lane/II, unknown placement/fusion, queue depths beyond the serial contract and
# profiles on the linear bypass are config errors, never silent defaults; strict mode
# refuses an active operation with no declared profile.
expect_failure invalid-lanes
expect_failure invalid-units
expect_failure invalid-queue
expect_failure invalid-residency
expect_failure invalid-placement
expect_failure invalid-fusion
expect_failure invalid-zero-ii
expect_failure invalid-linear-profile
expect_failure invalid-approximation
expect_failure invalid-linear-approx
expect_failure invalid-supported-ops
expect_failure unsupported-op
expect_failure unsupported-softmax
expect_failure strict-precision-mismatch
expect_failure strict-missing-approximation
expect_failure strict-missing-profile

echo "NPUsim SFU validation checks passed"
