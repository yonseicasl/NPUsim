# NPUsim: Full-System, Cycle-Accurate, Value-Aware Functional Simulations of DNN Accelerators
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
Rapid advances in deep learning techniques have sparked competitive races for developing dedicated hardware for deep neural network (DNN) acceleration. As the development of DNN accelerators becomes increasingly complicated to incorporate advanced optimization features exploiting variable dataflows, bit precision, data sparsity, etc. There are strong demands on simulation frameworks for DNN accelerators to facilitate their design space exploration and rapid pre-silicon analyses. Simulations for DNN accelerators address different modeling challenges compared to conventional processors (i.e., CPUs and GPUs). The execution of DNN on the accelerator is driven by dataflow and mapping instead of instructions. Hence, implementing variable dataflows and numerous data mapping options is crucial to cover a wide variety of accelerator designs and DNN models. In addition, neural accelerators heavily utilize value-dependent optimization features to enhance execution efficiency, such as zero-skipping, compression, quantization, etc. Enabling these advanced optimization techniques requires value-aware functional executions with cycle-accurate timing simulations.

_NPUsim_ is a simulation framework for DNN accelerator that features full-system, cycle-accurate, and value-aware functional simulations. The framework consists of three modules; input, main, and output modules. The input module processes system specifications to construct a target accelerator, including dataflow and mapping schemes. Interfacing with PyTorch, the input module also handles DNN workload configurations and neural datasets. The main module includes a simulation engine and execution schedulers. The simulation engine models an accelerator component hierarchy and performs the cycle-accurate simulations. Execution schedulers control data movements between accelerator components along the memory hierarchy to drive an intended execution flow. The output module displays simulation results at the end of the simulation.

## Prerequisites
NPUsim uses g++ to compile C++ codes. If the gcc compiler is not installed, type the following command in terminal.
    
    $ sudo apt install build-essential

The simulation framework contains several extensions including Nebula, DRAMsim3, and PyTorch. Please refer to https://github.com/yonsei-icsl/nebula, https://github.com/umd-memsys/DRAMsim3, and https://pytorch.org.

## Download
The latest release of NPUsim simulation framework is v1.0 (as of Feb. 2024). To obtain a copy of NPUsim v1.0, use the following command in a terminal 

    $ git clone --branch v1.0 https://github.com/yonsei-icsl/npusim

Or, you can download the latest development version, clone the git repository as is.

    $ git clone https://github.com/yonsei-icsl/npusim

## Build
NPUsim provides a script file name npusim.sh to facilitate the buid and run the simulator. To build the entire modules of NPUsim, type the following commnds in the main directory of NPUsim.

    $ cd npusim/
    $ ./npusim.sh build all

Or, you can specify NPUsim modules by typeing a command in the following format

    $ ./npusim.sh build <module>

## Simulation
After NPUsim framework is built, the simulator becomes ready to perform DNN simulation on an accelerator. The accelerator specification, DNN model, and scheduling table is required to dispatch simulation. The npusim.sh file facilitates the simulation. A simulation command follows the format shown below. The <accelerator> field represents the accelerator specification, the <DNN> model indicates the configuration of DNN, and <scheduling table> shows scheduling schemes including dataflows and data mappings.

    $ ./npusim.sh run <accelerator> <DNN model> <scheduling table>

For example, the following command simulates the execution of AlexNet on Eyeriss with energy optimized scheduling scheme.

    $./npusim.sh run eyeriss alexnet energy

## Reference and Contact
To reference NPUsim, please use our ModSim paper.

    @misc{kim_modsim2021,
        author      = {B. Kim and C. park and T. Lim and W. Song},
        title       = {{NPUsim: Full-System, Cycle-Accurate, Functional Simulations of Deep Neural Network Accelerators}},
        booktitle   = {Workshop on Modeling and Simulation of Systems and Applications}, 
        month       = {Oct.},
        year        = {2021},
        pages       = {1-2},
    }
For troubleshooting, bug reports, or any questions regarding NPUsim simualtion framework, please contact Bogil Kim via email: bogilkim {\at} yonsei {\dot} ac {\dot} kr. Or, visit our lab webpage: https://casl.yonsei.ac.kr