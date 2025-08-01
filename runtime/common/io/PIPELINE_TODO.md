# BioGPU Pipeline Unification TODO

## Overview
This checklist tracks the steps needed to unify the resistance and genes pipelines into a single executable that will be placed at `/usr/local/bin/bio_gpu_pipeline` for the streaming service to execute.

## Prerequisites
- [ ] Ensure CUDA toolkit is installed
- [ ] Verify CMake version 3.18+ is available
- [ ] Check that reference databases exist:
  - [ ] `/data/microbial_ref_db`
  - [ ] `/data/quinolone_resistance_db`

## Phase 1: Create Shared Infrastructure (Day 1)

### 1.1 Directory Structure
- [ ] Create shared components directory:
  ```bash
  mkdir -p /home/david/projects/biogpu/runtime/kernels/shared
  ```

### 1.2 Move Duplicate Components
- [ ] Move `sample_csv_parser.{h,cpp}` to `shared/`
  ```bash
  mv runtime/kernels/resistance/sample_csv_parser.* runtime/kernels/shared/
  rm runtime/kernels/genes/sample_csv_parser.*  # Remove duplicate
  ```
- [ ] Move `bloom_filter.cu` to `shared/`
  ```bash
  mv runtime/kernels/resistance/bloom_filter.cu runtime/kernels/shared/
  rm runtime/kernels/genes/bloom_filter.cu  # Remove duplicate
  ```
- [ ] Move `fastq_reader.{h,cpp}` to `shared/`

### 1.3 Update Build Configurations
- [ ] Update `runtime/kernels/resistance/CMakeLists.txt` to reference shared components
- [ ] Update `runtime/kernels/genes/CMakeLists.txt` to reference shared components
- [ ] Create `runtime/kernels/shared/CMakeLists.txt` for shared library

## Phase 2: Create Unified Pipeline (Day 2-3)

### 2.1 Design Unified Entry Point
- [ ] Create `runtime/kernels/unified_amr_pipeline.cpp` with command-line parsing:
  ```cpp
  // Must accept these arguments (from biogpu_service.py):
  // --r1 <path>
  // --r2 <path>
  // --output-dir <path>
  // --reference-db <path>
  // --resistance-db <path>
  // --gpu-device <id>
  // --sample-id <string>
  // --progress-json
  // --threads <num>
  ```

### 2.2 Implement Progress Reporting
- [ ] Create JSON progress output format:
  ```json
  {
    "stage": "minimizer_screening",
    "progress": 25,
    "message": "Processing batch 1000/4000"
  }
  ```
- [ ] Add progress hooks to each pipeline stage

### 2.3 Refactor Existing Pipelines
- [ ] Extract resistance detection logic into `ResistanceDetector` class
- [ ] Extract AMR gene detection logic into `AMRGeneDetector` class
- [ ] Define shared interfaces:
  ```cpp
  class DetectorInterface {
    virtual void processBatch(const ReadBatch& batch) = 0;
    virtual void getResults(json& results) = 0;
  };
  ```

### 2.4 Implement Single-Pass FASTQ Processing
- [ ] Create unified FASTQ reader that feeds both detectors
- [ ] Implement shared k-mer/minimizer generation
- [ ] Add batch queuing system for GPU processing

## Phase 3: Fix Known Issues (Day 3-4)

### 3.1 GPU Memory Management
- [ ] Fix GPU memory errors in translated search (from genes/CLAUDE.md)
- [ ] Implement shared GPU memory pools
- [ ] Add proper CUDA error checking

### 3.2 Performance Optimizations
- [ ] Load protein database once at startup (not per batch)
- [ ] Implement bloom filter caching
- [ ] Optimize batch sizes for GPU utilization

## Phase 4: Build and Install (Day 4)

### 4.1 Create Build System
- [ ] Create master `CMakeLists.txt` for unified pipeline
- [ ] Add all necessary CUDA flags and optimizations
- [ ] Include all required libraries (zlib, pthread, etc.)

### 4.2 Build Executable
- [ ] Build the unified pipeline:
  ```bash
  cd runtime/kernels
  mkdir build && cd build
  cmake ..
  make -j8
  ```
- [ ] Test executable locally with sample data

### 4.3 Install to System
- [ ] Copy executable to system location:
  ```bash
  sudo cp bio_gpu_pipeline /usr/local/bin/
  sudo chmod +x /usr/local/bin/bio_gpu_pipeline
  ```
- [ ] Verify installation:
  ```bash
  /usr/local/bin/bio_gpu_pipeline --help
  ```

## Phase 5: Integration Testing (Day 5)

### 5.1 Test with Streaming Service
- [ ] Start infrastructure services:
  ```bash
  docker-compose up -d redis kafka zookeeper
  ```
- [ ] Run streaming service:
  ```bash
  source venv/bin/activate
  python3 run_service.py
  ```
- [ ] Submit test job via manifest
- [ ] Verify WebSocket updates are received
- [ ] Check output files are generated correctly

### 5.2 Validation Tests
- [ ] Compare results with individual pipelines
- [ ] Verify both resistance mutations and AMR genes are detected
- [ ] Test with various FASTQ file sizes
- [ ] Benchmark performance improvements

### 5.3 Error Handling
- [ ] Test with corrupted FASTQ files
- [ ] Test with missing reference databases
- [ ] Verify error messages are properly streamed
- [ ] Check retry logic works correctly

## Phase 6: Production Deployment (Day 6)

### 6.1 Documentation
- [ ] Update README with unified pipeline usage
- [ ] Document all command-line arguments
- [ ] Create example output format documentation

### 6.2 Performance Tuning
- [ ] Profile GPU utilization
- [ ] Optimize batch sizes based on GPU memory
- [ ] Set appropriate thread counts

### 6.3 Final Verification
- [ ] Run full test suite
- [ ] Process production-scale datasets
- [ ] Monitor resource usage
- [ ] Verify streaming updates work at scale

## Success Criteria

The unified pipeline is complete when:
1. ✓ Single executable handles both resistance and gene detection
2. ✓ Executable is installed at `/usr/local/bin/bio_gpu_pipeline`
3. ✓ Streaming service successfully executes pipeline via Kafka jobs
4. ✓ Real-time progress updates stream via WebSocket
5. ✓ Results match individual pipeline outputs
6. ✓ Performance is improved (single FASTQ read, shared GPU memory)
7. ✓ All tests pass

## Quick Commands Reference

```bash
# Build unified pipeline
cd /home/david/projects/biogpu/runtime/kernels
mkdir build && cd build
cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda
make -j8

# Install
sudo cp bio_gpu_pipeline /usr/local/bin/

# Test
/usr/local/bin/bio_gpu_pipeline \
  --r1 /data/test_R1.fastq.gz \
  --r2 /data/test_R2.fastq.gz \
  --output-dir /tmp/test_output \
  --reference-db /data/microbial_ref_db \
  --resistance-db /data/quinolone_resistance_db \
  --gpu-device 0 \
  --sample-id TEST001 \
  --progress-json \
  --threads 8

# Start streaming service
cd /home/david/projects/biogpu/runtime/common/io
source venv/bin/activate
python3 run_service.py
```

## Notes

- The pipeline must output JSON progress updates to stdout for the streaming service to capture
- All actual results should be written to files in the output directory
- The executable should return 0 on success, non-zero on failure
- GPU device selection must be respected for multi-GPU systems