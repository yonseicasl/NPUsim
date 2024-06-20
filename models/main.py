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

from torchsummary import summary

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


def load_data(m_network, m_iteration):
    device = 'cpu'
    t_network = m_network
    sub_list = ['networks/', '.cfg']
    for sub in sub_list:
        t_network = t_network.replace(sub,'')    

    # Initialize DNN data
    dataset_str = parsing.parsing_dataset(t_network)
    dataset = data_loader.load_data(dataset_str)

    count = 0
    for data in dataset:
        if count == m_iteration:
            images, labels = data[0].to(device), data[1].to(device)
            break;
        count += 1
    return images, labels

def print_output(m_layers, m_index) :
    print(m_layers[m_index].out)

def forward(m_network, m_image, m_label, m_iteration):
    device = 'cpu'
    m_network.to(device)

    m_network.eval()
    correct  = 0
    total = 0
    count = 0

    outputs = m_network(m_image)
    layers = network_builder.get_layers(m_network)


    _, predicted = torch.max(outputs, 1)
    correct += (predicted == m_label).sum().item()
    total += 1
    print (str(correct) + "/" + str(total))

def flatten(m_input) :
    output = torch.flatten(m_input, 1)
    return output

def layerwise_forward(m_layer, m_layer_name, m_input, m_index) :

    output = m_layer[m_index](m_input)
    return output

def init_weight(m_layer, m_index) :
    return m_layer[m_index].weight

def print_result(m_input, m_label):
    correct = 0
    total = 0
    _, predicted = torch.max(m_input, 1)
    correct += (predicted == m_label).sum().item()
    total += 1
    print (str(correct) + "/" + str(total))
     

if __name__ == '__main__':

    parser = network_parser()
    args = parser.parse_args()
    network = args.n

    DNN_model, DNN_layers, DNN_layers_name = init(network)

    print(len(DNN_layers_name))
    #for i in range(len(DNN_layers)):
    #    if DNN_layers_name[i] == "Conv2d":
    #        print(DNN_layers_name[i] + " " + str(DNN_layers[i].in_channels) + " " + str(DNN_layers[i].out_channels))
    #    elif DNN_layers_name[i] == "Linear":
    #        print(DNN_layers_name[i] + " " + str(DNN_layers[i].in_features) + " " + str(DNN_layers[i].out_features))

    image, label = load_data(network, 0)

    forward(DNN_model, image, label, 0)

    outputs = []
    output = None
    total = 0
    correct = 0

    for i in range(len(DNN_layers)):
        if i == 0 :
            output = layerwise_forward(DNN_layers, DNN_layers_name, image, i)
        else :
            output = layerwise_forward(DNN_layers, DNN_layers_name, output, i)
    print_result(output, label)
