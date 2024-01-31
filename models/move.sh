#!/bin/bash

if [[ "$#" -lt 2 ]]; then
    echo -e "Error"
    exit 0
fi

target=$1; shift
network=$1; shift
metric=$1; shift

cd ..

if [[ ! -d result ]]; then
    mkdir result
fi

if [[ ! -d result/$target ]]; then
    mkdir result/$target
fi

if [[ ! -d result/$target/$network ]]; then
    mkdir result/$target/$network
fi

if [[ ! -d result/$target/$network/$metric ]]; then
    mkdir result/$target/$network/$metric
fi

for i in {0..198}
do 
    if [[ -f models/${target}_${network}_layer_${i}.txt ]]; then
        mv models/${target}_${network}_layer_${i}.txt result/$target/$network/$metric/layer_${i}.txt
    fi
done

mv models/${target}_${network}.txt result/$target/$network/$metric/network.txt
if [[ -f models/${target}_DRAM/dramsim3.txt ]]; then
    mv models/${target}_DRAM/dramsim3.txt models/${target}_DRAM/${target}_${network}-${metric}.txt
fi

