#!/bin/bash
# Commands to install HDF5 libraries for BioGPU

# Update package lists
sudo apt-get update

# Install HDF5 C/C++ development libraries
sudo apt-get install -y \
    libhdf5-dev \
    libhdf5-103-1t64 \
    libhdf5-cpp-103-1t64 \
    libhdf5-serial-dev \
    hdf5-tools

# The jsoncpp library (already installed on your system)
# sudo apt-get install -y libjsoncpp-dev libjsoncpp25

echo "Installation complete!"
echo ""
echo "HDF5 headers will be installed at:"
echo "  /usr/include/hdf5/serial/H5Cpp.h"
echo "  /usr/include/hdf5/serial/hdf5.h"
echo ""
echo "Libraries will be at:"
echo "  /usr/lib/x86_64-linux-gnu/hdf5/serial/"