#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 ACCELERATOR NETWORK MAPPING" >&2
    exit 1
fi

target=$1
network=$2
metric=$3
identifier_pattern='^[A-Za-z0-9._-]+$'
if [[ ! $target =~ $identifier_pattern || ! $network =~ $identifier_pattern || ! $metric =~ $identifier_pattern ]]; then
    echo "Error: result identifiers may contain only letters, numbers, '.', '_', and '-'." >&2
    exit 1
fi

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
model_dir="$repo_dir/models"
result_dir="$repo_dir/result/$target/$network/$metric"
mkdir -p "$result_dir"

for ((i = 0; i <= 198; i++)); do
    layer_result="$model_dir/${target}_${network}_layer_${i}.txt"
    if [[ -f $layer_result ]]; then
        mv "$layer_result" "$result_dir/layer_${i}.txt"
    fi
done

network_result="$model_dir/${target}_${network}.txt"
if [[ ! -f $network_result ]]; then
    echo "Error: simulation result not found: $network_result" >&2
    exit 1
fi
mv "$network_result" "$result_dir/network.txt"

dram_result="$model_dir/${target}_DRAM/dramsim3.txt"
if [[ -f $dram_result ]]; then
    mv "$dram_result" "$model_dir/${target}_DRAM/${target}_${network}-${metric}.txt"
fi
