#!/bin/bash

if [[ "$#" -lt 1 ]]; then 
    echo -e "Error"
    exit 0
fi


network=$1; shift

if [[ ! -d weights ]]; then
    mkdir weights
fi

if [[ ! -f weights/$network.wgh ]]; then
    echo -e "$network.wgh not exist"
    ../ext/nebula/weight_dropbox.sh $network large
    mv weights/input.wgh weights/$network.wgh
fi

if [[ ! -d datasets ]]; then
    mkdir datasets
fi

if [[ ! -d datasets/imagenet ]]; then 
    echo -e "Dataset imagenet not exist"
    ../ext/nebula/dataset_dropbox.sh imagenet large
    mv datasets/imagenet_large datasets/imagenet
fi

