import os
import sys
pwd_path = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
util_path = os.path.join(pwd_path, 'python_util')
sys.path.append(util_path)

import argparse
import re

import parsing
import network_builder
import data_loader
import execution

from torchvision import models
import torch 
import torch.nn as nn


def network_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', 
                        metavar='network.cfg', 
                        type=str,
                        help='Network configuration file', default='NONE',
                        )
    parser.add_argument('-d',
                        type=str,
                        help='Dataset type', default='cifar10',
                        ) 
    return parser

def init(m_network):
    # Initialize the name of DNN model
    t_network = m_network
    sub_list = ['networks/', '.cfg']
    for sub in sub_list:
        t_network = t_network.replace(sub,'')    

    # Initialize DNN data
    dataset_str = parsing.parsing_dataset(t_network)
    dataset = data_loader.load_data(dataset_str)

    DNN_layers, DNN_layers_name = network_builder.build_network(t_network)

    #for i in DNN_layers_name:
    #    print(i)

    return DNN_layers, DNN_layers_name

def forward(layer, x):
    layer_name = network_builder.extract_layer_name(layer) 
    print(layer_name)

if __name__ == '__main__':

    parser = network_parser()
    args = parser.parse_args()
    network = args.n
    dataset = args.d

    DNN_layers, DNN_layers_name = init(network)

    for i in DNN_layers:
        forward(i, 0)
