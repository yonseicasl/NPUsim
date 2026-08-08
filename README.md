# NPUsim: Cycle-Level and Value-Aware Simulator for DNN Accelerators
Developed by Bogil Kim, Chanho Park, and William J. Song\
Computer Architecture and Systems Lab, Yonsei University\
Current release: v1.0 (Feb. 2024)

## Table of Contents
1. [Introduction](#introduction)
2. [Current Implementation Status](#current-implementation-status)
3. [Prerequisites](#prerequisites)
4. [Download](#download)
5. [Build](#build)
6. [Simulation](#simulation)
7. [Validation](#validation)
8. [Reference and Contact](#reference-and-contact)

## Introduction
Rapid advancements in deep learning techniques have sparked competitive races to develop specialized hardware for accelerating deep neural networks (DNNs). As the development of DNN accelerators becomes increasingly complex to incorporate advanced optimization features exploiting variable dataflows, bit precisions, data sparsity, etc., there is a strong demand for simulation frameworks for DNN accelerators to facilitate their design space explorations and rapid pre-silicon analyses. However, simulating DNN accelerators addresses different challenges from modeling conventional processes (e.g., CPUs, GPUs) in that accelerator executions are i) driven by dataflows rather than instructions, ii) associated with spatiotemporal data mapping in the hardware, and iii) influenced by value-based optimization controls, such as zero skipping.

Neural computations are typically expressed as nested loops in program code. Previous studies have employed loop-centric models, where hardware characteristics are implicitly embedded in the loops. This approach requires rearranging the entire loop sequence to modify execution flows or hardware features, which significantly restricts the tools’ usability and modeling coverage. The open-source NPUsim framework aims to overcome these modeling challenges and improve self-evaluate-only research practices in the accelerator domain.

_NPUsim_ is an architecture-oriented modeling framework in which accelerator components are explicit and dataflow mappings are separated from neural-layer loop implementations. The current code provides cycle-level timing/energy simulation for convolutional and fully connected layers, plus an optional value-aware functional path. See the scope below before interpreting network totals.

## Current Implementation Status

- Accelerator timing/energy is modeled for convolutional and fully connected layers.
- Pooling, shortcut/residual, concat, activation, softmax, upsample, and excitation layers are not assigned accelerator timing or memory-traffic cost. The simulator emits a warning for every excluded layer and marks the network total as partial.
- Functional forwarding is optional (`FUNCTIONAL=1`) and does not make excluded accelerator costs part of the timing total.
- Multi-chip object counts and active mapping dimensions are validated against physical capacity before execution.
- `configs/` is the authoritative configuration root. Files formerly copied under `models/` are not used by `npusim.sh run`.
- Included accelerator and mapping files are configuration assets, not calibration evidence. Transformer mappings do not currently imply end-to-end Transformer execution support.
- There is no committed golden output or paper/hardware calibration dataset yet. Do not treat uncalibrated totals as validated reference results.

## Prerequisites
NPUsim requires g++ to compile C++ code.
    
    $ sudo apt install build-essential

The simulation framework uses Nebula and DRAMSim3, and optionally integrates with system-installed PyTorch/OpenCV. `npusim.sh` checks out pinned Nebula and DRAMSim3 revisions for reproducible builds; the revisions may be overridden with `NEBULA_COMMIT` and `DRAMSIM3_COMMIT`. PyTorch is not cloned or built by NPUsim.

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

## Runtime data formats

Data precision is selected at runtime in the `[accelerator]` section. `data_format` applies to all tensors; `input_format`, `weight_format`, and `output_format` override it. Supported storage/accounting formats are `fp32`, `fp16`, `bf16`, `int8`, `int4`, `int2`, `uint8`, `mxfp8`, and `mxfp4`. MXFP uses 32-element blocks with an 8-bit shared scale, and its scale metadata is included in PE local-buffer capacity accounting.

```ini
[accelerator]
data_format = bf16
weight_format = mxfp4
output_format = int8
accumulator_format = fp32
```

The current runtime implementation establishes format parsing and physical bit accounting. Format-aware packed storage, conversion, and MAC IP latency/energy remain follow-up work described in `plan/plan_datatype.md`.

The default build is a timing/energy build (`FUNCTIONAL=0`, `DRAMSIM3=0`). These options are defined near the top of `npusim.sh`; enable only the modes and dependencies required for the experiment.

## Simulation
Running NPUsim takes an accelerator name, DNN network name, and mapping name. The launcher resolves all three from the root `configs/` directory and fails before initialization if any selected file is missing.

Network totals cover convolutional and fully connected layers only. Non-supported layers still participate in optional Nebula functional forwarding, but their accelerator latency, energy, and memory traffic are excluded and reported as such.

    $ ./npusim.sh run <accelerator> <DNN model> <scheduling table>

For example, the following command simulates the execution of AlexNet on Eyeriss with energy optimized scheduling scheme.

    $ ./npusim.sh run eyeriss alexnet energy

## Validation

Run the dependency-free validation suite before submitting a change:

    $ ./tests/run_validation.sh

It parses all committed accelerator/network/mapping files, checks accelerator component cardinality and multi-chip dimensions, exercises strict parser failure cases, and checks build/run script regressions. GitHub Actions also runs this suite with AddressSanitizer and UndefinedBehaviorSanitizer.

This suite does not replace an end-to-end golden-result comparison, full simulator sanitizer run, functional comparison against PyTorch/Nebula, or latency/energy calibration. Those require the external dependencies and reference datasets and remain future validation work.

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
