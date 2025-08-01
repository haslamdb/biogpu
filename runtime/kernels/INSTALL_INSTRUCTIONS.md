# Installation Instructions for BioGPU Dependencies

## Required Dependencies

The BioGPU unified pipeline requires the following libraries to build successfully:

### 1. HDF5 C/C++ Libraries
HDF5 is used for efficient data storage and alignment output.

**To install:**
```bash
sudo apt-get update
sudo apt-get install -y \
    libhdf5-dev \
    libhdf5-103-1t64 \
    libhdf5-cpp-103-1t64 \
    libhdf5-serial-dev \
    hdf5-tools
```

### 2. JSON C++ Library (jsoncpp)
Already installed on your system at `/usr/include/jsoncpp/`

## Quick Installation

Run the provided installation script:
```bash
cd /home/david/projects/biogpu/runtime/kernels
sudo ./install_dependencies.sh
```

Or install manually:
```bash
# Update package lists
sudo apt-get update

# Install HDF5
sudo apt-get install -y libhdf5-dev libhdf5-cpp-103-1t64 libhdf5-serial-dev

# jsoncpp is already installed
```

## Building After Installation

Once dependencies are installed:

### 1. Build Resistance Pipeline
```bash
cd resistance
chmod +x build_resistance_with_deps.sh
./build_resistance_with_deps.sh
```

### 2. Test Unified Pipeline
```bash
cd ../build_complete
./test_integrated_pipeline.sh
```

## Verification

After installation, these files should exist:
- `/usr/include/hdf5/serial/H5Cpp.h` - HDF5 C++ headers
- `/usr/include/hdf5/serial/hdf5.h` - HDF5 C headers
- `/usr/include/jsoncpp/json/json.h` - JSON C++ headers ✓ (already present)

## Compiler Flags

When building with HDF5:
- **Include path**: `-I/usr/include/hdf5/serial`
- **Library path**: `-L/usr/lib/x86_64-linux-gnu/hdf5/serial`
- **Link flags**: `-lhdf5_cpp -lhdf5`

For jsoncpp:
- **Include path**: `-I/usr/include/jsoncpp`
- **Link flags**: `-ljsoncpp`

## Troubleshooting

If you encounter linking errors:
1. Check that all packages installed correctly
2. Run `ldconfig` to update library cache
3. Verify library paths with `pkg-config --libs hdf5`

## Next Steps

After installing dependencies:
1. Build the resistance pipeline
2. Run the integrated pipeline with both components
3. Test with real FASTQ data