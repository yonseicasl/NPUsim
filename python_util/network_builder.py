# Package for neural network model
from torchvision import models
import torch
import torch.nn as nn
import torch.nn.functional as function

#from torchvision import models

# Package for dataset and transformation

import parsing

class network_configuration:
    # Initialize network configuration
    def __init__(self, section):

        # Initialize input data configuration
        self.batch_size = int(section.get_configuration_value("batch") if section.get_configuration_value("batch") != None else self.batch_size)
        self.input_height = int(section.get_configuration_value("height") if section.get_configuration_value("height") != None else self.input_height)
        self.input_width = int(section.get_configuration_value("width") if section.get_configuration_value("width") != None else self.input_width)
        self.input_channel = int(section.get_configuration_value("channels") if section.get_configuration_value("channels") != None else self.input_channel)

        self.input_size = self.batch_size * self.input_height * self.input_width * self.input_channel

        print(str(self.input_height) + ", " + str(self.input_width) + ", " + str(self.input_channel))

class layer_configuration:
    def __init__(self, network, prev_layer, section):
        self.name = section.name

        self.activation_type = section.get_configuration_value("activation") if section.get_configuration_value("activation") != None else "linear"
        self.batch_normalize = bool(section.get_configuration_value("batch_normalize") if section.get_configuration_value("batch_normalize") != None else "false")

    def activation(self):
        if self.activation_type == "elu":
            return nn.ELU()
        elif self.activation_type == "gelu":
            return nn.GELU()
        elif self.activation_type == "leaky":
            return nn.LeakyReLU()
        elif self.activation_type == "linear":
            return None
        elif self.activation_type == "mish":
            return nn.Mish()
        elif self.activation_type == "relu":
            return nn.ReLU()
        elif self.activation_type == "sigmoid":
            return nn.Sigmoid()
        elif self.activation_type == "silu":
            return nn.SiLU()
        elif self.activation_type == "softplus":
            return nn.SoftPlus()
        elif self.activation_type == "tanh":
            return nn.Tanh()

class convolutional_layer_configuration(layer_configuration):
    # Initialize layer configuration
    def __init__(self, network, prev_layer, section):
        super().__init__(network, prev_layer, section)

        self.name = section.name    # Layer name
        self.prev_layer = prev_layer
        
        self.activation = layer_configuration.activation(self)

        self.batch_size = network.batch_size # Batch size

        # Initialize input data configuration
        self.input_height = int(prev_layer.output_height if prev_layer else network.input_height)
        self.input_width = int(prev_layer.output_width if prev_layer else network.input_width)
        self.input_channel = int(prev_layer.output_channel if prev_layer else network.input_channel)
        
        # Initialize weight configuration
        self.weight_height = int(section.get_configuration_value("size") if section.get_configuration_value("size") != None else self.weight_height)
        self.weight_width = self.weight_height
        self.weight_height = int(section.get_configuration_value("filter_height") if section.get_configuration_value("filter_height") != None else self.weight_height)
        self.weight_width = int(section.get_configuration_value("filter_width") if section.get_configuration_value("filter_width") != None else self.weight_width)

        # Initialize stride of weight
        self.stride_height = int(section.get_configuration_value("stride") if section.get_configuration_value("stride") != None else self.stride_height)
        self.stride_width = self.stride_height
        self.stride_height = int(section.get_configuration_value("stride_height") if section.get_configuration_value("stride_height") != None else self.stride_height)
        self.stride_width = int(section.get_configuration_value("stride_width") if section.get_configuration_value("stride_width") != None else self.stride_width)

        # Initialize padding of input data 
        self.padding_height = int(self.weight_height/2)
        self.padding_width = int(self.weight_width/2)

        # Initialize output data configuration
        self.output_height = int((self.input_height + 2*self.padding_height - self.weight_height) / self.stride_height + 1)
        self.output_width = int((self.input_width + 2*self.padding_width - self.weight_width) / self.stride_width + 1)
        self.output_channel = int(section.get_configuration_value("filters") if section.get_configuration_value("filters") != None else self.output_channel)

        self.input_size = self.batch_size * self.input_height * self.input_width * self.input_channel
        self.weight_size = self.input_channel * self.weight_height * self.weight_width * self.output_channel
        self.output_size = self.batch_size * self.output_height * self.output_width * self.output_channel

        self.build_convolutional_layer()

        print(self.name + ": " + str(self.input_height) + "*" + str(self.input_width) + "*" + str(self.input_channel) + ", " + str(self.weight_height) + "*" + str(self.weight_width) + "*" + str(self.output_channel) + ", " + str(self.output_height) + "*" + str(self.output_width) + "*" + str(self.output_channel))

    def build_convolutional_layer(self):
        if self.batch_normalize:
            self.model = nn.Sequential(
                         nn.Conv2d(in_channels=self.input_channel, out_channels=self.output_channel, kernel_size=(self.weight_height, self.weight_width), stride=(self.stride_height, self.stride_width), padding=(self.padding_height, self.padding_width)),
                         nn.BatchNorm2d(self.output_channel),
                         self.activation
                         )
        else:
            self.model = nn.Sequential(
                         nn.Conv2d(in_channels=self.input_channel, out_channels=self.output_channel, kernel_size=(self.weight_height, self.weight_width), stride=(self.stride_height, self.stride_width), padding=(self.padding_height, self.padding_width)),
                         self.activation
                         )

    def forward(self, x):
        self.output = self.model(x)
        return self.output

class connected_layer_configuration(layer_configuration):
    def __init__(self, network, prev_layer, section):
        super().__init__(network, prev_layer, section)
        self.name = section.name
        self.prev_layer = prev_layer

        self.input_size = int(prev_layer.output_size if prev_layer else network.input_size)
        self.output_size = int(section.get_configuration_value("output") if prev_layer else self.output_size)
        self.weight_size = self.input_size * self.output_size

        self.build_connected_layer()

        print(self.name + ": " + str(self.input_size) + ", " + str(self.weight_size) + ", " + str(self.output_size))
    def build_connected_layer(self):
        if self.batch_normalize:
            self.model = nn.Sequential(
                         nn.Linear(in_features = self.input_size, out_features = self.output_size),
                         #activation()
                         )
        else:
            self.model == nn.Sequential(
                          nn.Linear(in_features = self.input_size, out_features = self.output_size),
                          #activation()
                          ) 
    def forward(self, x):
        self.output = self.model(x)
        return self.output

class pooling_layer_configuration(layer_configuration):
    def __init__(self, network, prev_layer, section):

        super().__init__(network, prev_layer, section)
        self.name = section.name 
        self.prev_layer = prev_layer

        self.batch_size = network.batch_size

        self.input_height = int(prev_layer.output_height if prev_layer else network.input_height)
        self.input_width = int(prev_layer.output_width if prev_layer else network.input_width)
        self.input_channel = int(prev_layer.output_channel if prev_layer else network.input_channel)

        self.weight_height = int(section.get_configuration_value("size") if section.get_configuration_value("size") != None else self.weight_height)
        self.weight_width = self.weight_height
        self.weight_height = int(section.get_configuration_value("weight_height") if section.get_configuration_value("weight_height") != None else self.weight_height)
        self.weight_width = int(section.get_configuration_value("weight_width") if section.get_configuration_value("weight_width") != None else self.weight_width)

        self.stride_height = int(section.get_configuration_value("stride") if section.get_configuration_value("stride") != None else self.stride_height)
        self.stride_width = self.stride_height
        self.stride_height = int(section.get_configuration_value("stride_height") if section.get_configuration_value("stride_height") != None else self.stride_height)
        self.stride_width = int(section.get_configuration_value("stride_width") if section.get_configuration_value("stride_width") != None else self.stride_width)

        self.padding_height = int(self.weight_height/2)
        self.padding_width = int(self.weight_width/2)
        self.padding_height = int(section.get_configuration_value("padding") if section.get_configuration_value("padding") != None else self.padding_height)
        self.padding_width = int(section.get_configuration_value("padding") if section.get_configuration_value("padding") != None else self.padding_width)
        self.padding_height = int(section.get_configuration_value("padding_height") if section.get_configuration_value("padding_height") != None else self.padding_height)
        self.padding_width = int(section.get_configuration_value("padding_width") if section.get_configuration_value("padding_width") != None else self.padding_width)
        
        self.output_height = int((self.input_height + 2*self.padding_height - self.weight_height) / self.stride_height + 1)
        self.output_width = int((self.input_width + 2*self.padding_width - self.weight_width) / self.stride_width +1)
        self.output_channel = self.input_channel

        self.input_size = self.batch_size * self.input_height * self.input_width * self.input_channel
        self.output_size = self.batch_size * self.output_height * self.output_width * self.output_channel

        self.build_pooling_layer()

        print(self.name + ": " + str(self.input_height) + "*" + str(self.input_width) + "*" + str(self.input_channel) + ", " + str(self.weight_height) + "*" + str(self.weight_width) + "*" + str(self.output_channel) + ", " + str(self.output_height) + "*" + str(self.output_width) + "*" + str(self.output_channel))
    
    def build_pooling_layer(self):
        if self.name == "maxpool":
            self.model = nn.MaxPool2d(kernel_size = (self.weight_height, self.weight_width), stride = (self.stride_height, self.stride_width), padding = (self.padding_height, self.padding_width))
        elif self.name == "avgpool":
            self.model = nn.AvgPool2d(kernel_size = (self.weight_height, self.weight_width), stride = (self.stride_height, self.stride_width), padding = (self.padding_height, self.padding_width))
        else:
            print("Error: Invalid layer name" + self.name)
            exit()
    def forward(self, x):
        self.output = self.model(x)
        return self.output

class shortcut_layer_configuration(layer_configuration):
    def __init__(self, network, prev_layer, section):
        super().__init__(network, prev_layer, section)
        self.name = section.name
        self.prev_layer = prev_layer

        self.batch_size = network.batch_size

        self.hops = int(section.get_configuration_value("hops") if section.get_configuration_value("hops") != None else 1)
        self.connection = self
        for i in range(self.hops):
            self.connection = self.connection.prev_layer


        self.output_height = prev_layer.output_height if prev_layer else network.input_height
        self.output_width = prev_layer.output_width if prev_layer else network.input_width
        self.output_channel = prev_layer.output_channel if prev_layer else network.input_channel
        self.output_size = self.batch_size * self.output_height * self.output_width *self.output_channel
        self.input_size = self.output_size

        print(self.name + ": " + str(self.output_height) + "*" + str(self.output_width) + "*" + str(self.output_channel))

    def forward(self, x):
        self.output = x + self.connection.output
        return self.output
        

def build_network(m_network):
    if m_network.lower() == "inceptionv3":
        DNN_model = models.inceptionv3(pretrained=False)
    elif m_network.lower() == "resnet50":
        DNN_model = models.resnet50(pretrained=False)
    elif m_network.lower() == "alexnet":
        DNN_model = models.alexnet(pretrained=True)
    else:
        DNN_model = build_network_from_scratch(m_network.lower())

    DNN_layers = get_layers(DNN_model)
    DNN_layers_name = []
    for i in DNN_layers:
        DNN_layers_name.append(extract_layer_name(i))
    #return DNN_model, DNN_layers_name
    return DNN_layers, DNN_layers_name


def build_network_from_scratch(network):
    sections = parsing.parsing_config(network)
    network
    layers = []
    for section in sections:
        # Parsing network configuration
        if section.name == 'net':
            network = network_configuration(section)
        elif section.name == 'data':
            continue
        elif section.name == 'convolutional':
            layer = convolutional_layer_configuration(network, layers[len(layers)-1] if len(layers) else None, section)
            layers.append(layer)
        elif section.name == 'connected':
            layer = connected_layer_configuration(network, layers[len(layers)-1] if len(layers) else None, section)
            layers.append(layer)
        elif section.name == 'avgpool' or section.name == 'maxpool':
            layer = pooling_layer_configuration(network, layers[len(layers)-1] if len(layers) else None, section)
            layers.append(layer)
        elif section.name == 'shortcut':
            layer = shortcut_layer_configuration(network, layers[len(layers)-1] if len(layers) else None, section)
            layers.append(layer)
        else:
            continue
    return layers

def get_layers(m_network):
    layers = list(m_network.children())
    layer_list = []

    if layers == []:
        return m_network
    else:
        for layer in layers:
            try:
                layer_list.extend(get_layers(layer))
            except TypeError:
                layer_list.append(get_layers(layer))
    return layer_list

def extract_layer_name(m_layer):
    layer_str = str(m_layer)
    index = layer_str.find('(')
    return layer_str[:index]
