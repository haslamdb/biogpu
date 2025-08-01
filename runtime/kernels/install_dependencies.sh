#!/bin/bash
# Script to install HDF5 and other dependencies for BioGPU pipelines

set -e

echo "========================================="
echo "Installing BioGPU Pipeline Dependencies"
echo "========================================="
echo ""
echo "This script will install:"
echo "- HDF5 C/C++ development libraries"
echo "- JSON C++ library (jsoncpp)"
echo ""
echo "You will need sudo privileges."
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Installation cancelled."
    exit 1
fi

echo ""
echo "Updating package lists..."
sudo apt-get update

echo ""
echo "Installing HDF5 libraries..."
sudo apt-get install -y \
    libhdf5-dev \
    libhdf5-103-1t64 \
    libhdf5-cpp-103-1t64 \
    libhdf5-serial-dev \
    hdf5-tools

echo ""
echo "Installing JSON C++ library..."
sudo apt-get install -y \
    libjsoncpp-dev \
    libjsoncpp25

echo ""
echo "Verifying installations..."
echo ""

# Check HDF5
echo -n "HDF5 C headers: "
if [ -f "/usr/include/hdf5/serial/hdf5.h" ]; then
    echo "✓ Found at /usr/include/hdf5/serial/hdf5.h"
else
    echo "✗ Not found"
fi

echo -n "HDF5 C++ headers: "
if [ -f "/usr/include/hdf5/serial/H5Cpp.h" ]; then
    echo "✓ Found at /usr/include/hdf5/serial/H5Cpp.h"
else
    echo "✗ Not found"
fi

# Check jsoncpp
echo -n "JSON C++ headers: "
if [ -f "/usr/include/jsoncpp/json/json.h" ]; then
    echo "✓ Found at /usr/include/jsoncpp/json/json.h"
else
    echo "✗ Not found"
fi

echo ""
echo "Setting up environment..."

# Create a helper script for compiler flags
cat > hdf5_flags.sh << 'EOF'
#!/bin/bash
# Helper script to get HDF5 compiler flags

# HDF5 flags
HDF5_CFLAGS="-I/usr/include/hdf5/serial"
HDF5_LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5"

# JSON flags
JSON_CFLAGS="-I/usr/include/jsoncpp"
JSON_LDFLAGS="-ljsoncpp"

echo "HDF5 compile flags: $HDF5_CFLAGS"
echo "HDF5 link flags: $HDF5_LDFLAGS"
echo "JSON compile flags: $JSON_CFLAGS"
echo "JSON link flags: $JSON_LDFLAGS"
EOF

chmod +x hdf5_flags.sh

echo ""
echo "Installation complete!"
echo ""
echo "To use HDF5 in your builds, add these flags:"
echo "  Compile: -I/usr/include/hdf5/serial"
echo "  Link: -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5"
echo ""
echo "For jsoncpp:"
echo "  Compile: -I/usr/include/jsoncpp"
echo "  Link: -ljsoncpp"
echo ""
echo "A helper script 'hdf5_flags.sh' has been created with these flags."
echo ""
echo "You can now build the resistance pipeline!"