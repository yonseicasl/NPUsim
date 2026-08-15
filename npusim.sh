#!/bin/bash

##### NPUsim directories #####
# NPUsim main directory
npusimdir=$PWD
# NPUsim scheduler directory 
schedulerdir=$npusimdir/scheduler
# NPUsim Utility directory
utildir=$npusimdir/utils
# Accelerator component directory
componentdir=$npusimdir/components
# Extension directory 
extdir=$npusimdir/ext
# Library directory
libdir=$npusimdir/library
# NPUsim library
lib=$libdir/libnpusim.so

##### External directory #####
# Nebula directory
nebuladir=$extdir/nebula
# Pytorch directory
pytorchdir=$extdir/pytorch
# DRAMsim directory
dramsim3dir=$extdir/DRAMsim3
# NPUsim model directory
modeldir=$npusimdir/models


##### NPUsim compile options #####
# C++ compiler 
cc=g++
# Linker
lc=$cc
# Makefile options
mopt="CC=$cc DIR=$npusimdir LIB=$lib"
# Compiler options
ccopt="-Wall -O3 -fPIC -I$schedulerdir -I$utildir -I$componentdir"
# Linker options
ldopt="-L$libdir"
libopt="-lnpusim"
# C++ version
stdc="-std=c++11"

# Options for Nebula and DRAMsim3
ccopt+=" -I$nebuladir/common -I$nebuladir/models/layers -I$nebuladir/models/networks -I$dramsim3dir/src"
ldopt+=" -L$nebuladir/library -L$dramsim3dir"
libopt+=" -lnebula -ldramsim3 -lopenblas -lpthread -lz `pkg-config --libs opencv4`"


##### NPUsim build options ###

# Print out the operations sequence
PRINT=0
# Debug.
DEBUG=0
# Functional simulation
FUNCTIONAL=0
# Using integer
USER_INTEGER=0
# Using float 
USER_FLOAT=0
# DRAMsim
DRAMSIM3=0
# PyTorch
Pytorch=0

##### Append Makefile options #####

# Message printing
if [[ $PRINT -eq 1 ]]; then
    ccopt+=" -DPRINT"
fi

# Debugging
if [[ $DEBUG -eq 1 ]]; then
    ccopt+=" -g"
fi

# Functional simulation
if [[ $FUNCTIONAL -eq 1 ]]; then
    ccopt+=" -DFUNCTIONAL"
fi

# Using integer
if [[ $USER_INTEGER -eq 1 ]]; then
    ccopt+=" -DUSER_INTEGER"
fi

#Using float
if [[ $USER_FLOAT -eq 1 ]]; then
    ccopt +=" -DUSER_FLOAT"
fi

# DRAMsim3 Connection
if [[ $DRAMSIM3 -eq 1 ]]; then
    ccopt+=" -DDRAMSIM3"
fi

# PyTorch connection
if [[ $Pytorch -eq 1 ]]; then
    #ccopt+=" -DPytorch `python-config --includes`"
    #libopt+=" -lm -ldl -lutil `python-config --libs`"
    ccopt+=" -I/usr/include/python3.8 -DPytorch"
    libopt+=" -lpython3.8 -lm -ldl -lutil `pkg-config --cflags --libs python3`"
fi

# Makefile MFLAG
mflag="$mopt LC=$lc"
# Makefile CCFLAG
ccflag="CCFLAG=\"$ccopt\""
# Makefile LDFLAG
ldflag="LDFLAG=\"$ldopt\""
# Makefile LIBFLAG
libflag="LIBFLAG=\"$libopt\""
# Makefile STD
std="STD=\"$stdc\""

##### NPUsim usage #####
function print_help {
    echo -e "Usage: $0 [build/clean/run] [accelerator target] [network target] [Mapping options]"
    exit 0
}

##### NPUsim build function #####
function build_model {
    #echo -e "NPUsim build $1"
    case "$1" in 
        all)
            # Make extension directory
            if [[ ! -d $extdir ]]; then
                mkdir $extdir
            fi

            # Build Nebula library
            # Install Nebula software framework from Github
            if [[ ! -d $nebuladir ]]; then
                echo -e "\n# Install Nebula framework from Github"
                cd $extdir; git clone --branch npusim --single-branch https://github.com/yonsei-icsl/nebula
            fi

            # Build Nebula library
            echo -e "\n# Build Nebula software framework"
            cd $nebuladir; ./nebula.sh build lib

            # Install PyTorch framework from Github
            if [[ ! -d $pytorchdir ]]; then
                echo -e "\n# Install PyTorch software framework"
                cd $extdir; git clone https://github.com/pytorch/pytorch.git
            fi

            # Build DRAMSim3
            # Install DRAMSim3 from Github
            if [[ ! -d $dramsim3dir ]]; then
                echo -e "\n# Install DRAMSim3"
                cd $extdir; git clone --branch master --single-branch https://github.com/umd-memsys/DRAMsim3.git
            fi
            # Build DRAMSim3
            echo -e "\n# Build DRAMsim3"
            cd $dramsim3dir; make libdramsim3.so

            # Build NPUsim library
            echo -e "\n# Build NPUsim library"
            cd $libdir; eval $mflag $ccflag $std make

            # Build NPUsim executable file
            echo -e "\n# Build NPUsim model"
            cd $modeldir; eval $mflag $ccflag $ldflag $libflag $std EXE='model' make
            ;;
        dramsim3)
            # Make extension directory
            if [[ ! -d $extdir ]]; then
                mkdir $extdir
            fi

            # Build DRAMSim3
            # Install DRAMSim3 from Github
            if [[ ! -d $dramsim3dir ]]; then
                echo -e "\n# Install DRAMSim3"
                cd $extdir; git clone --branch master --single-branch https://github.com/umd-memsys/DRAMsim3.git
            fi
            # Build DRAMSim3
            echo -e "\n# Build DRAMsim3"
            cd $dramsim3dir; make libdramsim3.so
            ;;
        ext)
            # Make extension directory
            if [[ ! -d $extdir ]]; then
                mkdir $extdir
            fi

            # Build Nebula library
            # Install Nebula software framework from Github
            if [[ ! -d $nebuladir ]]; then
                echo -e "\n# Install Nebula framework from Github"
                cd $extdir; git clone --branch npusim --single-branch https://github.com/yonsei-icsl/nebula
            fi

            # Build Nebula library
            echo -e "\n# Build Nebula software framework"
            cd $nebuladir; ./nebula.sh build lib

            # Install PyTorch framework from Github
            if [[ ! -d $pytorchdir ]]; then
                echo -e "\n# Install PyTorch software framework"
                cd $extdir; git clone https://github.com/pytorch/pytorch.git
            fi

            # Build DRAMSim3
            # Install DRAMSim3 from Github
            if [[ ! -d $dramsim3dir ]]; then
                echo -e "\n# Install DRAMSim3"
                cd $extdir; git clone --branch master --single-branch https://github.com/umd-memsys/DRAMsim3.git
            fi
            # Build DRAMSim3
            echo -e "\n# Build DRAMsim3"
            cd $dramsim3dir; make libdramsim3.so
            ;;
        lib)
            # Build NPUsim library
            echo -e "\n# Build NPUsim library"
            cd $libdir; eval $mflag $ccflag $std make
            ;;
        model)
            # Build NPUsim executable file
            echo -e "\n# Build NPUsim model"
            cd $modeldir; eval $mflag $ccflag $ldflag $libflag $std EXE='model' make
            ;;
        nebula)
            # Make extension directory
            if [[ ! -d $extdir ]]; then
                mkdir $extdir
            fi

            # Build Nebula library
            # Install Nebula software framework from Github
            if [[ ! -d $nebuladir ]]; then
                echo -e "\n# Install Nebula framework from Github"
                cd $extdir; git clone --branch npusim --single-branch https://github.com/yonsei-icsl/nebula
            fi

            # Build Nebula library
            echo -e "\n# Build Nebula software framework"
            cd $nebuladir; ./nebula.sh build lib
            ;;
        npusim)
            # Build NPUsim library and executable file
            # Build NPUsim library
            echo -e "\n# Build NPUsim library"
            cd $libdir; eval $mflag $ccflag $std make

            # Build NPUsim executable file
            echo -e "\n# Build NPUsim model"
            cd $modeldir; eval $mflag $ccflag $ldflag $libflag $std EXE='model' make
            ;;
        pytorch)
            # Make extension directory
            if [[ ! -d $extdir ]]; then
                mkdir $extdir
            fi
            # Install PyTorch framework from Github
            if [[ ! -d $pytorchdir ]]; then
                echo -e "\n# Install PyTorch software framework"
                cd $extdir; git clone https://github.com/pytorch/pytorch.git
            fi
            ;;
        *)
            # Print out error message
            echo -e "# Usage : ./npusim.sh build [options]"
            echo -e "\n List of [options]"
            echo -e "# all      : All"
            echo -e "# dramsim3 : DRAMsim3"
            echo -e "# ext      : Extensions libraries (DRAMSim3, Nebula)"
            echo -e "# lib      : NPUsim library"
            echo -e "# nebula   : Nebula software framework"
            echo -e "# npusim   : NPUsim library and executable file"
            echo -e "# model    : NPUsim model"
            exit 0;
            ;;
    esac
}

##### NPUsim clean function #####
function clean_model {
    case "$1" in 
        all)
            # Cleaning Nebula 
            echo -e "\n# Cleaning DNN framework"
            cd $nebuladir; ./nebula.sh clean lib;

            # Cleaning DRAMsim3 
            echo -e "\n# Cleaning DRAMsim3"
            cd $dramsim3dir; make clean;

            # Cleaning NPUsim library
            echo -e "\n# Cleaning NPUsim library"
            cd $libdir; eval $mflag make clean;
            
            # Cleaning NPUsim executable file
            echo -e "\n# Cleaning NPUsim model"
            cd $modeldir; eval EXE='model' make clean;
            ;;
        dramsim3)
            # Cleaning DRAMsim3
            echo -e "\n# Cleaning DRAMsim3"
            cd $dramsim3dir; make clean;
            ;;
        ext)
            # Cleaning Nebula
            echo -e "\n# Cleaning DNN software framework"
            cd $nebuladir; ./nebula.sh clean lib;

            # Cleaning DRAMsim3
            echo -e "\n# Cleaning DRAMsim3"
            cd $dramsim3dir; make clean;
            ;;
        lib)
            # Cleaning NPUsim library
            echo -e "\n# Cleaning NPUsim library"
            cd $libdir; eval $mflag make clean;
            ;;
        nebula)
            # Cleaning Nebula
            echo -e "\n# Cleaning DNN software framework"
            cd $nebuladir; ./nebula.sh clean lib;
            ;;
        npusim)
            # Cleaning NPUsim library and executable file
            # Cleaning NPUsim library
            echo -e "\n# Cleaning NPUsim library"
            cd $libdir; eval $mflag make clean;

            # Cleaning NPUsim executable file
            echo -e "\n# Cleaning NPUsim model"
            cd $modeldir; eval EXE='model' make clean;
            ;;
        model)
            # Cleaning NPUsim executable file
            echo -e "\n# Cleaning NPUsim model"
            cd $modeldir; eval EXE='model' make clean;
            ;;
        *)
            # Print out error message
            echo -e "# Usage : ./npusim.sh clean [options]"
            echo -e "\n List of [options]"
            echo -e "# all      : All"
            echo -e "# dramsim3 : DRAMSim3 Library"
            echo -e "# ext      : Extensions library (DramSim3, Nebula)" 
            echo -e "# lib      : NPUsim library"
            echo -e "# nebula   : Nebula software framework"
            echo -e "# npusim   : NPUsim library and executable file"
            echo -e "# model    : NPUsim model"
            exit 0;
            ;;
    esac
}

function run_npusim {
    echo -e "\n# Running $network on $target"
    cd $modeldir
    eval ./model run $target $network $mapping
}

# Print usage help.
if [[ $1 = '-h' || $1 = '-help' ]]; then 
    print_help
fi

action=$1; shift
target=$1; shift
network=$1; shift
mapping=$1; shift

case "$action" in 
    build) 
        # Action for build
        build_model $target
        ;;
    clean)
        # Action for clean
        clean_model $target
        ;;
    run)
        # Action for running
        run_npusim $target $network $mapping
        ;;
    help) 
        # Action for help
        print_help
        ;;
    *)
        echo -e "# Error : unknown action $action"
        echo -e "# Usage : ./npusim.sh [action]"
        echo -e "\nList of [action]"
        echo -e "# build : Compile"
        echo -e "# clean : Clean" 
        echo -e "# help  : NPUsim commands overview"
        echo -e "# run   : Execute simulation"
        exit 0
esac

