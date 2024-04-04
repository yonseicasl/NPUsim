#!/bin/bash

if [[ "$#" -lt 1 ]]; then 
    echo -e "Error"
    exit 0
fi


network=$1; shift

if [[ ! -d weights ]]; then
    mkdir weights
fi

../ext/nebula/weight_dropbox.sh $network large
mv input.wgh weights/$network.wgh 


if [[ ! -d datasets ]]; then
    mkdir datasets
fi


../ext/nebula/dataset_dropbox.sh imagenet large
