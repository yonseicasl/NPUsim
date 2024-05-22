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
    m_dataset = parsing.parsing_dataset(t_network)
    data_loader.load_data(m_dataset)

    print(t_network)
    DNN_model, DNN_layers_list = network_builder.build_network(t_network)
    return DNN_model

def forward(m_layer, x):
    
    execution.forward(m_layer, x) 

    return m_layer

if __name__ == '__main__':

    print('run Python main function')
    parser = network_parser()
    args = parser.parse_args()
    network = args.n
    dataset = args.d

    DNN_model = init(network)

