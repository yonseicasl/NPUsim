#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT

cxx=${CXX:-g++}
sanitizers=${NPUSIM_TEST_SANITIZERS:-}
cxxflags=(-std=c++11 -Wall -Wextra -Werror -I"$repo_dir/utils" -I"$repo_dir/scheduler")
if [[ $sanitizers == 1 ]]; then
    cxxflags+=(-fsanitize=address,undefined -fno-omit-frame-pointer)
fi

"$cxx" "${cxxflags[@]}" \
    "$repo_dir/unittest/validation_test.cc" \
    "$repo_dir/utils/config.cc" \
    "$repo_dir/utils/utils.cc" \
    "$repo_dir/utils/datatype.cc" \
    "$repo_dir/utils/interconnect_timing.cc" \
    "$repo_dir/utils/energy_units.cc" \
    "$repo_dir/scheduler/mapping_table.cc" \
    "$repo_dir/utils/pe_lane.cc" \
    -o "$test_dir/validation_test"

mapfile -d '' config_files < <(
    find "$repo_dir/configs" -type f \( -name '*.cfg' -o -name '*.map' \) -print0 | sort -z
)
if [[ ${#config_files[@]} -eq 0 ]]; then
    echo "No configuration files found" >&2
    exit 1
fi
"$test_dir/validation_test" "${config_files[@]}"

expect_failure() {
    if "$test_dir/validation_test" "$1" >/dev/null 2>&1; then
        echo "Expected validation failure for $1" >&2
        exit 1
    fi
}

expect_success() {
    if ! "$test_dir/validation_test" "$1" >/dev/null 2>&1; then
        echo "Expected validation success for $1" >&2
        "$test_dir/validation_test" "$1" >&2 || true
        exit 1
    fi
}

printf 'key=value\n[section]\nother=value\n' >"$test_dir/setting-before-section.cfg"
printf '[section\nkey=value\n' >"$test_dir/malformed-section.cfg"
printf '[section]\nkey=value\nkey=duplicate\n' >"$test_dir/duplicate-key.cfg"
# E8: an energy unit cost that cannot mean anything must be rejected, for any energy key --
# negative, non-numeric, non-finite, or declared with no value at all. The accelerator
# validator keys off an "/accelerators/" path component, so the fixtures live under one.
mkdir -p "$test_dir/accelerators"
printf '[accelerator]\nnum_chips=1\n[shared]\nread_energy=1.0:-2.0:3.0\n' >"$test_dir/accelerators/negative-energy.cfg"
printf '[accelerator]\nnum_chips=1\n[shared]\nwrite_energy=abc\n' >"$test_dir/accelerators/nonnumeric-energy.cfg"
printf '[accelerator]\nnum_chips=1\n[systolic_array]\ncomputation_energy=-0.5\n' >"$test_dir/accelerators/negative-compute-energy.cfg"
printf '[accelerator]\nnum_chips=1\n[shared]\nstatic_energy=1e400\n' >"$test_dir/accelerators/infinite-energy.cfg"
# RE5: schema-level rejection. A misspelled key, a short/long per-datatype vector and an empty
# middle field all used to pass value validation and silently leave a cost at zero.
printf '[accelerator]\nnum_chips=1\n[systolic_array]\ncomputaton_energy=0.5\n' >"$test_dir/accelerators/energy-key-typo.cfg"
printf '[accelerator]\nnum_chips=1\n[shared]\nread_energy=1.0:2.0\n' >"$test_dir/accelerators/short-energy-vector.cfg"
printf '[accelerator]\nnum_chips=1\n[shared]\nread_energy=1.0:2.0:3.0:4.0\n' >"$test_dir/accelerators/long-energy-vector.cfg"
printf '[accelerator]\nnum_chips=1\n[shared]\nread_energy=1.0::3.0\n' >"$test_dir/accelerators/empty-energy-field.cfg"
printf '[accelerator]\nnum_chips=1\n[systolic_array]\ncomputation_energy=1.0:2.0\n' >"$test_dir/accelerators/vector-for-scalar-energy.cfg"
# RE3: a setup energy with no setup execution has no event source.
printf '[accelerator]\nnum_chips=1\n[systolic_array]\nlayer_setup_energy=5.0\nlayer_setup_cycle=0\n' >"$test_dir/accelerators/setup-energy-no-event.cfg"
# ... and a well-formed zero cost must still be ACCEPTED: 0 is a legitimate modeled cost, so the
# guard must not turn "declared zero" into an error.
printf '[accelerator]\nnum_chips=1\ncompression_type=dense\n[systolic_array]\nmemory_type=separate\ninput_size=1\nweight_size=1\noutput_size=1\nnumber_of_macs=1\nmac_width=1\nbitwidth=8\nheight=1\nwidth=1\nnoc=bus\ncomputation_energy=0\nnoc_energy=0\n[shared]\nmemory_size=1\nbitwidth=8\nnoc=bus\nread_energy=0:0:0\n[multi_chip]\nmemory_type=shared\nmemory_size=1\nheight=1\nwidth=1\nbitwidth=8\nnoc=bus\n[dram]\nbitwidth=8\n' >"$test_dir/accelerators/zero-energy-ok.cfg"
printf '[layer]\nmac=1,1,1,1,1,1,1,1,1,1\n' >"$test_dir/too-few.map"
printf '[layer]\nmac=1,1,1,1,1,1,1,1,1,1,1,1\n' >"$test_dir/too-many.map"
printf '[layer]\nmac=1,1,1,1,1,1,1,1,1,invalid,1\n' >"$test_dir/non-numeric.map"
printf '[layer]\nmac=1,1,1,1,1,1,1,0,0,1,0\n' >"$test_dir/zero-stride.map"
printf '[layer]\nmac=4294967295,1,1,1,1,1,1,0,0,1,1\npe=2,1,1,1,1,1,1,0,0,1,1\n' >"$test_dir/overflow.map"

expect_failure "$test_dir/setting-before-section.cfg"
expect_failure "$test_dir/malformed-section.cfg"
expect_failure "$test_dir/duplicate-key.cfg"
expect_failure "$test_dir/accelerators/negative-energy.cfg"
expect_failure "$test_dir/accelerators/nonnumeric-energy.cfg"
expect_failure "$test_dir/accelerators/negative-compute-energy.cfg"
expect_failure "$test_dir/accelerators/infinite-energy.cfg"
expect_failure "$test_dir/accelerators/energy-key-typo.cfg"
expect_failure "$test_dir/accelerators/short-energy-vector.cfg"
expect_failure "$test_dir/accelerators/long-energy-vector.cfg"
expect_failure "$test_dir/accelerators/empty-energy-field.cfg"
expect_failure "$test_dir/accelerators/vector-for-scalar-energy.cfg"
expect_failure "$test_dir/accelerators/setup-energy-no-event.cfg"
# E8: a declared ZERO energy cost is a legitimate modeled value and must still be accepted --
# the guard rejects meaningless costs, not cheap ones.
expect_success "$test_dir/accelerators/zero-energy-ok.cfg"
expect_failure "$test_dir/too-few.map"
expect_failure "$test_dir/too-many.map"
expect_failure "$test_dir/non-numeric.map"
expect_failure "$test_dir/zero-stride.map"
expect_failure "$test_dir/overflow.map"

bash -n "$repo_dir/npusim.sh"
bash -n "$repo_dir/models/move.sh"
if "$repo_dir/npusim.sh" run '../invalid' alexnet energy >/dev/null 2>&1; then
    echo "Unsafe run identifier was accepted" >&2
    exit 1
fi

if "$repo_dir/models/move.sh" '../invalid' network mapping >/dev/null 2>&1; then
    echo "Unsafe result identifier was accepted" >&2
    exit 1
fi
if grep -R -E 'memset\([^;]*,[[:space:]]*1\.0' "$repo_dir/components" >/dev/null; then
    echo "Numeric initialization through memset reappeared" >&2
    exit 1
fi
# plan_sfu.md Phase-0: the removed hardcoded per-PE ReLU must never come back. The method
# is deleted from pe.h so a reintroduced CALL is a compile error; this guards against the
# method itself (declaration or unqualified/member call form) reappearing. Activation cost
# belongs to the SFU event contract and functional activation to Nebula.
if grep -R -n -E '(^|[^A-Za-z0-9_])activation\(\)' "$repo_dir/components" "$repo_dir/scheduler" \
        "$repo_dir/models" "$repo_dir/utils" --include='*.cc' --include='*.h' >/dev/null; then
    echo "Removed pe_t::activation() (declaration or call) reappeared" >&2
    exit 1
fi
if grep -q 'USE_INTEGER\|USE_FLOAT' "$repo_dir/npusim.sh"; then
    echo "Legacy build-option variable reappeared" >&2
    exit 1
fi
if grep -q 'schduler' "$repo_dir/library/Makefile"; then
    echo "Scheduler dependency typo reappeared" >&2
    exit 1
fi

python3 -m unittest discover -s "$repo_dir/unittest/python" -p 'test_*.py' -v

echo "NPUsim validation checks passed"
