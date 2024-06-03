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

    DNN_Model, DNN_layers, DNN_layers_name = network_builder.build_network(t_network)

    return DNN_Model, DNN_layers, DNN_layers_name


def load_data(m_network):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    t_network = m_network
    sub_list = ['networks/', '.cfg']
    for sub in sub_list:
        t_network = t_network.replace(sub,'')    

    # Initialize DNN data
    dataset_str = parsing.parsing_dataset(t_network)
    dataset = data_loader.load_data(dataset_str)

    return dataset


def forward(m_network, m_dataset, m_iteration, m_index):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    m_network.to(device)

    m_network.eval()
    correct  = 0
    total = 0
    count = 0

    for data in m_dataset:
        if count == m_iteration :
            images, labels = data[0].to(device), data[1].to(device)
            outputs = m_network(images)
            _, predicted = torch.max(outputs, 1)
            correct += (predicted == labels).sum().item()
            total += 1
            break
        count += 1
    print (str(correct) + "/" + str(total))


def layerwise_forward(m_layer, m_input_data, m_index) :

    output = []
    return output
     

if __name__ == '__main__':

    parser = network_parser()
    args = parser.parse_args()
    network = args.n

    DNN_model, DNN_layers, DNN_layers_name = init(network)
    dataset = load_data(network)

    #output = forward(DNN_model, dataset, 0, 0)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    for data in dataset:
        image, label = data[0].to(device), data[1].to(device)
        for i in range(len(DNN_layers)):
            print(i)
            if i == 0 :
                output = layerwise_forward(DNN_layers[i], image, i)
            else :
                output = layerwise_forward(DNN_layers[i], output, i)

    #for i in range(len(DNN_layers)):
    #    if i == 0 :
    #        x = forward(DNN_model, image, i)
    #    else :
    #        x = forward(DNN_model, x, i)

