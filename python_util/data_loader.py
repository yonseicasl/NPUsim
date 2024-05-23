import torchvision
from torch.utils.data import DataLoader, Dataset
from torchvision import datasets, transforms

import matplotlib.pyplot as plt
import numpy as np

import cv2


def load_data(m_dataset_str):

    if m_dataset_str.lower() == 'imagenet':
        mean = [0.485, 0.456, 0.406]
        std = [0.229, 0.224, 0.225]
    elif m_dataset_str.lower() == 'cifar10':
        mean = [0.5, 0.5, 0.5]
        std = [0.5, 0.5, 0.5]
    elif m_dataset_str.lower() == 'mnist':
        mean = [0.5]
        std = [1.0]

    transform = transforms.Compose([
        transforms.Resize(256),
        transforms.CenterCrop(227),
        transforms.ToTensor(),
        transforms.Normalize(mean, std)
    ])

    if m_dataset_str.lower() == "imagenet":
        print("ImageNet")
        test_set = datasets.ImageNet(root='./datasets', split='val',
                                     transform=transform)

    elif m_dataset_str.lower() == "cifar10":
        print("CIFAR10")
        test_set = datasets.CIFAR10(root='./datasets', train=False,
                                     download=True, transform=transform)
    elif m_dataset_str.lower() == "mnist":
        test_set = datasets.MNIST(root='./datasets', train=False,
                                   download=True, transform=transform)
    else :
        print('Not available ' + m_dataset_str)

    data_loader = DataLoader(test_set)
    #dataiter = iter(data_loader)
    #images, labels = next(dataiter)

    return data_loader

def load_data_from_scratch(m_dataset):
    return m_dataset
