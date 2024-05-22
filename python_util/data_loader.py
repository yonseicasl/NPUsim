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
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(mean, std)
    ])

    if m_dataset_str.lower() == "imagenet":
        m_dataset = datasets.ImageNet(root='./datasets', train=False,
                                      download=True, transform=transform)
    elif m_dataset_str.lower() == "cifar10":
        m_dataset = datasets.CIFAR10(root='./datasets', train=False,
                                     download=True, transform=transform)
    elif m_dataset_str.lower() == "mnist":
        m_dataset = datasets.MNIST(root='./datasets', train=False,
                                   download=True, transform=transform)
    else :
        print('Not available ' + m_dataset_str)
    data_loader = DataLoader(m_dataset)
    #dataiter = iter(data_loader)
    #images, labels = next(dataiter)
    #print(images)
    #img = torchvision.utils.make_grid(images)
    #img = img/2 + 0.5
    #npimg = img.numpy()
    #plt.imshow(np.transpose(npimg, (1,2,0)))
    #plt.show()

    #m_dataset = datasets.ImageNet('datasets/imagenet/test')
    #data_loader = DataLoader(m_dataset)
    #return data_loader


def load_data_from_scratch(m_dataset):
    return m_dataset
