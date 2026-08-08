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
    "$repo_dir/scheduler/mapping_table.cc" \
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

printf 'key=value\n[section]\nother=value\n' >"$test_dir/setting-before-section.cfg"
printf '[section\nkey=value\n' >"$test_dir/malformed-section.cfg"
printf '[section]\nkey=value\nkey=duplicate\n' >"$test_dir/duplicate-key.cfg"
printf '[layer]\nmac=1,1,1,1,1,1,1,1,1,1\n' >"$test_dir/too-few.map"
printf '[layer]\nmac=1,1,1,1,1,1,1,1,1,1,1,1\n' >"$test_dir/too-many.map"
printf '[layer]\nmac=1,1,1,1,1,1,1,1,1,invalid,1\n' >"$test_dir/non-numeric.map"
printf '[layer]\nmac=1,1,1,1,1,1,1,0,0,1,0\n' >"$test_dir/zero-stride.map"
printf '[layer]\nmac=4294967295,1,1,1,1,1,1,0,0,1,1\npe=2,1,1,1,1,1,1,0,0,1,1\n' >"$test_dir/overflow.map"

expect_failure "$test_dir/setting-before-section.cfg"
expect_failure "$test_dir/malformed-section.cfg"
expect_failure "$test_dir/duplicate-key.cfg"
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
if grep -q 'USE_INTEGER\|USE_FLOAT' "$repo_dir/npusim.sh"; then
    echo "Legacy build-option variable reappeared" >&2
    exit 1
fi
if grep -q 'schduler' "$repo_dir/library/Makefile"; then
    echo "Scheduler dependency typo reappeared" >&2
    exit 1
fi

echo "NPUsim validation checks passed"
