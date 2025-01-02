# NPUsim: Full-Model, Cycle-Level, and Value-Aware Functional Simulator for DNN Accelerators
Developed by Bogil Kim, Chanho Park, and William J. Song\
Computer Architecture and Systems Lab, Yonsei University\
Current release: v1.0 (Feb. 2024)

## Table of Contents
1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Download](#download)
4. [Build](#build)
5. [Simulation](#simulation)
6. [Reference and Contact](#reference-and-contact)

## Introduction
Rapid advancements in deep learning techniques have sparked competitive races to develop specialized hardware for accelerating deep neural networks (DNNs). As the development of DNN accelerators becomes increasingly complex to incorporate advanced optimization features exploiting variable dataflows, bit precisions, data sparsity, etc., there is a strong demand for simulation frameworks for DNN accelerators to facilitate their design space explorations and rapid pre-silicon analyses. However, simulating DNN accelerators addresses different challenges from modeling conventional processes (e.g., CPUs, GPUs) in that accelerator executions are i) driven by dataflows rather than instructions, ii) associated with spatiotemporal data mapping in the hardware, and iii) influenced by value-based optimization controls, such as zero skipping.

Neural computations are typically expressed as nested loops in program code. Previous studies have employed loop-centric models, where hardware characteristics are implicitly embedded in the loops. This approach requires rearranging the entire loop sequence to modify execution flows or hardware features, which significantly restricts the tools’ usability and modeling coverage. The open-source NPUsim framework aims to overcome these modeling challenges and improve self-evaluate-only research practices in the accelerator domain.

_NPUsim_ is a modeling framework for full-model, cycle-level, and value-aware functional simulations of DNN accelerators. Unlike the prior tools that built loop-centric models with the implicit embedding of hardware characteristics in the loops, NPUsim implements an architecture-oriented simulation framework where hardware components are explicitly created outside the nested loops of neural computations. NPUsim employs a novel simulation methodology that i) provides comprehensive accelerator models, ii) offers configurable dataflows and mappings that cover a wide range of accelerator designs by effectively separating architecture models from the loop sequence of neural layers, iii) enables cycle-level timing simulations, and iv) supports in-model functional execution with actual neural data, allowing for integrating value-based optimizations into the framework.

## Prerequisites
NPUsim requires g++ to compile C++ code.
    
    $ sudo apt install build-essential

The simulation framework includes several extensions, such as Nebula, DRAMsim3, and PyTorch. For more information, please refer to https://github.com/yonsei-icsl/nebula, https://github.com/umd-memsys/DRAMsim3, and https://pytorch.org.

## Download
The latest release of the NPUsim framework is v1.0 (as of Feb. 2024). To download NPUsim v1.0, use the following command in the terminal.

    $ git clone --branch v1.0 https://github.com/yonsei-icsl/npusim

Alternatively, you can download the latest development version by cloning the Git repo without specifying a branch.

    $ git clone https://github.com/yonsei-icsl/npusim

## Build
NPUsim provides a script file named `npusim.sh` to facilitate the building and running of the simulator. To build all modules of NPUsim, enter the following commands in the main directory of NPUsim.

    $ cd npusim/
    $ ./npusim.sh build all

To build a specific NPUsim module, it can be executed as follows.

    $ ./npusim.sh build <module>

## Simulation
Running NPUsim takes a set of inputs, including accelerator specifications (`<accelerator>`), DNN network configurations (`<DNN model>`), and a scheduling table (`<scheduling table>`) that specifies dataflow and mapping methods for the target accelerator. The `npusim.sh` script is again used for running the simulation.

    $ ./npusim.sh run <accelerator> <DNN model> <scheduling table>

For example, the following command simulates the execution of AlexNet on Eyeriss with energy optimized scheduling scheme.

    $./npusim.sh run eyeriss alexnet energy

## Reference and Contact
To reference NPUsim, please use our ModSim workshop white paper.

    @misc{kim_modsim2021,
        author      = {B. Kim and C. Park and T. Lim and W. Song},
        title       = {{NPUsim: Full-System, Cycle-Accurate, Functional Simulations of Deep Neural Network Accelerators}},
        booktitle   = {US DOE Workshop on Modeling and Simulation of Systems and Applications (ModSim)}, 
        month       = {Oct.},
        year        = {2021},
    }
For troubleshooting, bug reports, or any questions regarding the NPUsim simulation framework, please contact Bogil Kim via email: bogilkim {\at} yonsei {\dot} ac {\dot} kr. Or, visit our lab webpage: https://casl.yonsei.ac.kr.
