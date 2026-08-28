#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_binary="$(mktemp /tmp/npusim-workload-graph-test.XXXXXX)"
trap 'rm -f "$test_binary"' EXIT

python3 "$repo_root/unittest/python/test_pytorch_frontend.py" -v
python3 "$repo_root/unittest/python/test_pytorch_lowering.py" -v

g++ -std=c++11 -Wall -I"$repo_root/scheduler" \
    "$repo_root/scheduler/workload_graph.cc" \
    "$repo_root/scheduler/workload_graph_validation.cc" \
    "$repo_root/unittest/workload_graph_test.cc" \
    -o "$test_binary"
"$test_binary" "$repo_root/frontend/fixtures/linear_relu.exec.json"
"$test_binary" "$repo_root/frontend/fixtures/conv_relu.exec.json"
