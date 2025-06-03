# BioGPU: GPU-Accelerated Microbial Profiling and Antimicrobial Resistance Detection

## Overview

BioGPU is a domain-specific programming language and computational framework designed to accelerate microbial community profiling and antimicrobial resistance (AMR) detection using GPU computing. Initially focused on fluoroquinolone resistance detection through mutations in DNA gyrase (gyrA/gyrB) and topoisomerase IV (parC/parE) genes, BioGPU aims to provide real-time metagenomic analysis capabilities for clinical decision-making. 

## Project Goals

### Primary Objectives
1. **Rapid Microbial Community Profiling**: Map sequencing reads to comprehensive microbial genome databases to determine community structure and relative abundances in real-time
2. **Fluoroquinolone Resistance Detection**: Identify known resistance mutations in quinolone resistance-determining regions (QRDRs) and plasmid-mediated resistance genes
3. **Clinical Integration**: Generate actionable reports linking detected resistance to specific organisms with confidence scores, ultimately for informed treatment decisions

### Technical Features
- Native GPU acceleration for bioinformatics operations
- Domain-specific language that abstracts GPU complexity while maintaining performance
- Automatic optimization of sequence alignment algorithms for GPU architecture
- Real-time processing capabilities suitable for clinical environments

## Architecture

### Pipeline Structure
```
Raw Reads → Quality Control → Community Profiling → AMR Detection → Organism Attribution → Clinical Report
```

### Key Components

1. **BioGPU DSL**: High-level language for expressing bioinformatics pipelines
   - Intuitive syntax for pipeline definition
   - Built-in bioinformatics operations
   - Automatic GPU optimization

2. **Core Data Structures** (`biogpu_core.h`):
   - Compressed read storage for GPU efficiency
   - Resistance graph representation
   - GPU-optimized mutation indices
   - Community and resistance profile structures

3. **GPU Kernels**:
   - Sequence compression and indexing
   - Mutation scanning with Bloom filters
   - Community profiling algorithms
   - Resistance attribution logic

4. **VS Code Extension**: Syntax highlighting and code completion for BioGPU language

## Example Usage

```biogpu
pipeline FluoroquinoloneResistance {
    input: {
        patient_reads: fastq_file,
        patient_id: string,
        fq_concentration: float optional
    }
    
    output: {
        resistance_report: json,
        abundance_table: csv,
        clinical_summary: pdf
    }
    
    @gpu_kernel
    stage profile_microbiome {
        alignments = parallel_map(filtered_reads, microbe_genomes) {
            algorithm: minimap2_gpu,
            min_identity: 0.95
        }
        
        abundance = calculate_abundance(alignments) {
            method: "relative_abundance"
        }
    }
    
    @gpu_kernel
    stage detect_fq_resistance {
        mutations = scan_mutations(target_reads, fq_mutations) {
            positions: {
                gyrA: [83, 87],
                parC: [80, 84]
            }
        }
    }
    
    @gpu_kernel
    stage clinical_interpretation {
        resistance_calls = interpret_mutations(alignments) {
            thresholds: {
                high_confidence: "≥95% identity, ≥10x coverage",
                moderate_confidence: "≥90% identity, ≥5x coverage"
            }
        }
    }
}
```

## Clinical Applications

### 🧬 **Clinical Deployment**
- **Real-time FQ resistance screening**: Specialized pediatric infectious disease applications
- **Antimicrobial stewardship support**: Evidence-based therapy guidance

### 📊 **Validated Performance Metrics**
- **Speed**: Processing ~10 to 50M reads efficiently in batched workflow
- **Sensitivity**: High candidate detection rate 
- **Specificity**: Excellent alignment filtering 
- **Scalability**: Proven batch processing handles large clinical datasets

## Performance Targets

- ✅ Process >10 million reads in <2 minutes on single GPU
- ✅ >99% sensitivity for known resistance mutations
- ✅ >95% accuracy in organism-resistance attribution
- ✅ Support for real-time analysis during sequencing

## TODO List

### 1. Production Optimization and Clinical Integration

#### Immediate Workflow Optimizations
- [ ] **Expand FQDR sites in enhanced_kmer_builder.py**: Currenly only lists sites for E coli and E faecium
- [ ] **Batch size tuning**: Optimize batch size (currently 10K reads) for different GPU memory configurations


#### Enhanced Clinical Interpretation
- [ ] **Confidence scoring system**: Implement tiered confidence levels for resistance calls
       **Clinical output formats**: Add structured clinical reporting (FHIR, HL7)

### 2. Advanced Resistance Detection

#### Resistance Gene Quantification
- [ ] **Copy number estimation**: Quantify resistance gene abundance
- [ ] **Allele frequency analysis**: Support mixed populations and heteroresistance
- [ ] **Temporal tracking**: Monitor resistance evolution during treatment
- [ ] **Plasmid detection**: Identify mobile genetic elements carrying resistance

#### Expanded Resistance Classes
- [ ] **Beta-lactamases**: Carbapenem and ESBL resistance detection
- [ ] **Plasmid-mediated quinolone resistance (PMQR)**: qnr, aac(6')-Ib-cr genes
- [ ] **Vancomycin resistance genes**: vanA, vanB detection for enterococci
- [ ] **Aminoglycoside resistance**: aac, aph, ant gene families
- [ ] **Macrolide resistance**: erm, mef genes for respiratory pathogens

#### Possible future goals
- [ ] **Hospital LIS integration**: Develop interfaces for laboratory information systems
- [ ] **Real-time monitoring dashboard**: Clinical decision support interface
- [ ] **MIC prediction models**: Correlate mutations with quantitative resistance levels
- [ ] **Treatment recommendations**: Evidence-based antibiotic selection guidance
- [ ] **Resistance trend tracking**: Temporal analysis for hospital surveillance


### 3. Algorithm Enhancements

#### Microbiome profiling with K-mer Matching with Extension (Kraken2-style)
- [ ] **Three-phase k-mer matching algorithm**: Implement prefilter → match → extension workflow
  - **Phase 1: Bloom Filter Pre-screening**:
    - [ ] Space-efficient probabilistic filtering to eliminate non-matching k-mers
    - [ ] Multiple hash functions (typically 4-6) for low false positive rate
    - [ ] GPU implementation using bit arrays in shared memory
    - [ ] Target: <1% false positive rate with 4GB filter size
    ```cuda
    // Bloom filter structure
    struct BloomFilter {
        uint64_t* bit_array;
        size_t size_bits;
        int num_hash_functions;
        
        __device__ bool possibly_contains(uint64_t kmer) {
            for (int i = 0; i < num_hash_functions; i++) {
                uint64_t hash = murmur_hash(kmer, i);
                if (!test_bit(bit_array, hash % size_bits))
                    return false;
            }
            return true;  // Possibly in set (may be false positive)
        }
    };
    ```
  
  - **Phase 2: K-mer Matching with Minimizers**:
    - [ ] Extract minimizers from reads (smallest k-mer in window)
    - [ ] Match against pre-built genome k-mer index
    - [ ] Handle multiple genome matches (LCA resolution)
    - [ ] GPU kernel for parallel k-mer lookup
    ```cuda
    struct KmerHit {
        uint32_t read_position;
        uint32_t genome_id;
        uint32_t genome_position;
        float initial_score;
    };
    ```
  
  - **Phase 3: Extension and Verification**:
    - [ ] Extend matches by checking surrounding k-mers
    - [ ] Score based on consecutive matches and match density
    - [ ] Resolve ambiguous assignments using extension scores
    - [ ] Early termination for poor matches
    ```cuda
    __global__ void extend_kmer_matches(
        const char* reads,
        const KmerHit* initial_hits,
        const GenomeIndex* genome_index,
        ExtendedMatch* best_matches
    ) {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        KmerHit hit = initial_hits[tid];
        
        const int EXTENSION_LENGTH = 100;  // Check 100bp each direction
        const int K = 31;
        const int STEP = 5;  // Check every 5th k-mer
        
        float best_score = 0;
        int consecutive_matches = 0;
        
        // Extend in both directions
        for (int offset = -EXTENSION_LENGTH; offset <= EXTENSION_LENGTH; offset += STEP) {
            uint64_t read_kmer = extract_kmer(reads, hit.read_position + offset, K);
            uint64_t genome_kmer = extract_kmer(genome_index, hit.genome_id, 
                                              hit.genome_position + offset, K);
            
            if (read_kmer == genome_kmer) {
                consecutive_matches++;
                best_score += 1.0;
            } else {
                best_score -= 0.5;  // Penalize mismatches
                consecutive_matches = 0;
            }
            
            // Bonus for long consecutive matches
            if (consecutive_matches > 5) {
                best_score += 2.0;
            }
        }
        
        best_matches[tid].score = best_score;
        best_matches[tid].genome_id = hit.genome_id;
    }
    ```
  
  - **BioGPU DSL Implementation**:
    ```biogpu
    @gpu_kernel
    stage rapid_profile_with_extension {
        # Step 1: Bloom filter pre-screening
        potential_kmers = bloom_filter_check(reads) {
            filter_size: 4GB,
            hash_functions: 4,
            false_positive_rate: 0.01
        }
        
        # Step 2: K-mer matching with minimizers
        initial_matches = kmer_match(reads, genome_index) {
            k: 31,
            minimizer_window: 15,
            use_bloom_filter: potential_kmers
        }
        
        # Step 3: Extension and verification
        verified_matches = extend_matches(reads, initial_matches, genome_index) {
            extension_length: 100,
            scoring: {
                match: 1.0,
                mismatch: -0.5,
                consecutive_bonus: 2.0
            },
            min_extension_score: 10.0
        }
        
        # Step 4: Resolve ambiguous assignments
        final_assignments = resolve_taxonomy(verified_matches) {
            method: "weighted_lca",  # Or "best_hit"
            min_confidence: 0.8
        }
    }
    ```
  
  - **Performance Optimizations**:
    - [ ] Coalesced memory access by grouping similar k-mer hits
    - [ ] Shared memory caching of genome regions during extension
    - [ ] Warp-level primitives for score aggregation
    - [ ] Early termination when extension score drops below threshold
  
  - **Advantages over simple k-mer matching**:
    - Higher accuracy by verifying initial matches
    - Better strain-level resolution
    - Robust to sequencing errors
    - Can detect chimeric reads through poor extension scores

#### Translated Search Implementation (for FQ resistance testing)
- [ ] **Nucleotide-to-peptide alignment**: Implement 6-frame translation search
  - **Preprocessing Phase**:
    - [ ] SIMD optimization for codon lookup tables
    - [ ] Optional GPU translation kernel for very large datasets
    - [ ] K-mer pre-filtering to reduce search space
  - **GPU Kernel Design**:
    - [ ] Parallel alignment of translated peptides vs AMR protein database
    - [ ] Tiled loading of reference proteins into shared memory
    - [ ] Banded Smith-Waterman for efficient protein alignment
    - [ ] Early termination when alignment score drops below threshold
  - **Memory Optimization**:
    - [ ] Single large device buffer allocation with manual management
    - [ ] Pinned memory for efficient host-GPU transfers
    - [ ] Texture memory for BLOSUM62 substitution matrix
    - [ ] Coalesced memory access patterns for peptide sequences
  - **Specific Tools and Techniques**:
    - [ ] Build system: CMake + NVCC integration
    - [ ] Profiling: Nsight Compute for kernel optimization
    - [ ] Debugging: CUDA-GDB for kernel debugging, valgrind for host code
    - [ ] Performance targets: Process 1M translated frames in <30 seconds
- [ ] **Codon-aware alignment**: Optimize for protein-coding sequences
- [ ] **Stop codon handling**: Proper translation boundary detection
- [ ] **Genetic code variations**: Support alternative genetic codes for different organisms
- [ ] **Performance optimization**: GPU-accelerated translated search kernels

#### Advanced Sequence Analysis
- [ ] **GPU-accelerated Smith-Waterman**: Optimal local alignment on GPU
- [ ] **Profile HMM search**: Position-specific scoring matrices for gene families
- [ ] **Structural variant detection**: Large insertions/deletions in resistance genes
- [ ] **Phylogenetic analysis**: Strain typing and outbreak tracking

### 4. Database Development and Management

#### Resistance Database Expansion
- [ ] **CARD integration**: Comprehensive Antibiotic Resistance Database
- [ ] **Clinical isolate database**: Local hospital resistance patterns
- [ ] **Mutation effect prediction**: Functional impact scoring
- [ ] **Literature integration**: Automated PubMed resistance annotation

#### Genome Database Optimization
- [ ] **Incremental updates**: Live database updating system
- [ ] **Custom clinical genomes**: Integration of local isolate sequences
- [ ] **Taxonomic classification**: Species-level organism identification
- [ ] **Contamination detection**: Quality control for clinical samples

### 5. Clinical Validation and Quality Control

#### Validation Studies
- [ ] **Clinical isolate validation**: Compare with phenotypic AST results
- [ ] **Sensitivity/specificity studies**: ROC analysis for different resistance mechanisms
- [ ] **Inter-laboratory validation**: Multi-site clinical evaluation
- [ ] **Longitudinal studies**: Patient treatment outcome correlation

#### Quality Assurance
- [ ] **Control sample integration**: Positive/negative controls in workflows
- [ ] **Contamination monitoring**: Cross-sample contamination detection
- [ ] **Performance benchmarking**: Automated QC metrics and alerts
- [ ] **Proficiency testing**: External quality assessment participation

### 6. Infrastructure and Deployment

#### Scalability and Performance
- [ ] **Multi-GPU support**: Scale across multiple graphics cards
- [ ] **Cloud deployment**: AWS/Azure/GCP optimization
- [ ] **Container orchestration**: Kubernetes-based scaling
- [ ] **Load balancing**: Distribute clinical workloads efficiently

#### Security and Compliance
- [ ] **HIPAA compliance**: Patient data protection measures
- [ ] **Audit logging**: Complete analysis trail for regulatory compliance
- [ ] **Access control**: Role-based permissions for clinical users
- [ ] **Data encryption**: End-to-end security for sensitive medical data

### 7. Build Resistance Gene and Mutation Database
- [ ] ✅ **COMPLETED**: Design GPU-friendly mutation index structure from files in data/fq_resistance_db
- [ ] ✅ **COMPLETED**: Simplified data structure and index building
- [ ] ✅ **COMPLETED**: Working draft: `./build/fq_pipeline_gpu data/indices/fq_mutations/fq_mutation_index.h5`
- [ ] ✅ **COMPLETED**: Fixed k-mer matching phase
- [ ] ✅ **COMPLETED**: Mapping phase optimization (67% success rate achieved)
- [ ] Implement mutation confidence scoring system

### 8. Build Microbial Genome Database
- [ ] Create GPU-optimized k-mer index structure from downloaded genomes in data/genomes
- [ ] Implement genome database serialization format
- [ ] Add support for custom clinical isolate genomes
- [ ] Build incremental update mechanism

### 9. Adapt Current CPU Code to GPU
- [x] ✅ **COMPLETED**: Port sequence alignment algorithms to CUDA
- [x] ✅ **COMPLETED**: Implement GPU memory management for large datasets
- [x] ✅ **COMPLETED**: Create GPU kernel for parallel mutation scanning
- [x] ✅ **COMPLETED**: Optimize memory coalescing for sequence data
- [ ] Add CPU fallback for systems without GPU

### 10. Core Algorithm Implementation
- [x] ✅ **COMPLETED**: K-mer counting and indexing on GPU
- [x] ✅ **COMPLETED**: Minimizer-based sequence matching
- [x] ✅ **COMPLETED**: Bloom filter cascade for mutation filtering
- [ ] GPU-accelerated Smith-Waterman alignment
- [ ] Local assembly around resistance sites

### 11. Pipeline Development
- [x] ✅ **COMPLETED**: FASTQ parser with GPU streaming
- [x] ✅ **COMPLETED**: Quality control and filtering kernels
- [ ] Abundance calculation algorithms
- [x] ✅ **COMPLETED**: Resistance-organism attribution logic
- [ ] Clinical report generation

### 12. Language Infrastructure
- [x] ✅ **COMPLETED**: Implement lexer/parser for BioGPU syntax
- [x] ✅ **COMPLETED**: Create AST representation
- [x] ✅ **COMPLETED**: Build CUDA code generator
- [ ] Add LLVM-based optimization passes
- [ ] Implement type system for biological data

### 13. Testing and Validation
- [x] ✅ **COMPLETED**: Create synthetic test datasets with known resistance
- [x] ✅ **COMPLETED**: Validate with characterized clinical isolates
- [ ] Benchmark against existing tools (BWA, Kraken2)
- [ ] Compare with phenotypic AST results
- [ ] Performance profiling and optimization

### 14. Clinical Integration
- [ ] FHIR-compatible output format
- [ ] LIS system integration
- [ ] Real-time monitoring dashboard
- [ ] Automated alert system for critical resistance
- [ ] Clinical decision support rules

### 15. Documentation and Training
- [ ] API documentation
- [ ] Clinical user guide
- [ ] Performance tuning guide
- [ ] Example pipelines for common use cases
- [ ] Video tutorials for clinical users

### 16. Future Expansions
- [ ] Extend to other antibiotic classes (beta-lactams, aminoglycosides)
- [ ] Add long-read sequencing support
- [ ] Implement plasmid reconstruction
- [ ] Machine learning for novel resistance prediction
- [ ] Integration with structural biology for resistance mechanism analysis

### 17. GUI Development and User Interface

#### Web-Based Clinical Interface (Recommended)
- [ ] **Frontend Framework**: React or Vue.js for responsive clinical interface
  - [ ] Component library: Material-UI or Ant Design for professional look
  - [ ] Real-time updates: WebSocket connection for progress monitoring
  - [ ] Interactive visualizations: Plotly.js for results display
  - [ ] File management: Drag-and-drop upload with resumable transfers

- [ ] **Backend API**: FastAPI or Flask for Python integration
  - [ ] REST endpoints for job submission and status
  - [ ] WebSocket support for real-time progress streaming
  - [ ] Job queue management with Celery or RQ
  - [ ] Authentication: OAuth2/SAML for hospital SSO integration
  ```python
  # Example API structure
  @app.post("/api/v1/analysis")
  async def submit_analysis(
      files: List[UploadFile],
      config: AnalysisConfig,
      user: User = Depends(get_current_user)
  ):
      job_id = await biogpu_queue.submit(files, config, user)
      return {"job_id": job_id, "status": "queued"}
  
  @app.websocket("/ws/pipeline/{job_id}")
  async def pipeline_progress(websocket: WebSocket, job_id: str):
      await websocket.accept()
      async for update in monitor_job(job_id):
          await websocket.send_json(update)
  ```

#### Desktop Application Alternative
- [ ] **PyQt6 Implementation**: Native desktop app for offline use
  - [ ] Wizard-style interface for clinical users
  - [ ] Local file browser with batch processing support
  - [ ] Embedded results viewer with PDF export
  - [ ] System tray integration for background processing
  ```python
  class BiogpuWizard(QWizard):
      pages = [
          FileSelectionPage(),      # FASTQ file selection
          AnalysisConfigPage(),     # Resistance panels, settings
          GPUSelectionPage(),       # GPU device and memory
          ReviewPage(),             # Confirm settings
          ProgressPage(),           # Real-time monitoring
          ResultsPage()             # Interactive results
      ]
  ```

#### Streamlit Prototype (Quick Development)
- [ ] **Rapid prototyping interface**: For testing and development
  - [ ] Simple file upload interface
  - [ ] Checkbox-based analysis selection
  - [ ] Real-time progress tracking
  - [ ] Interactive results display
  ```python
  # Streamlit app structure
  st.set_page_config(page_title="BioGPU Clinical", layout="wide")
  
  # Sidebar: Configuration
  with st.sidebar:
      uploaded_files = st.file_uploader("FASTQ Files", accept_multiple_files=True)
      analysis_type = st.multiselect("Resistance Panels", 
          ["Fluoroquinolones", "Carbapenems", "Aminoglycosides"])
      
  # Main: Results display
  if st.button("Run Analysis"):
      progress = st.progress(0)
      status = st.empty()
      # Run pipeline with progress updates
  ```

#### Remote Data Streaming
- [ ] **S3/Cloud Storage Integration**:
  - [ ] Direct streaming from S3 buckets without local download
  - [ ] Support for presigned URLs for secure access
  - [ ] Chunked transfer with retry logic
  - [ ] Progress tracking for large files
  ```python
  async def stream_from_s3(bucket: str, key: str, pipeline: BiogpuPipeline):
      s3_client = aioboto3.client('s3')
      async with s3_client.get_object(Bucket=bucket, Key=key) as response:
          async for chunk in response['Body'].iter_chunks(chunk_size=1024*1024):
              await pipeline.feed_data(chunk)
  ```

- [ ] **HTTP/HTTPS Streaming**:
  - [ ] Support for remote FASTQ URLs
  - [ ] Resume capability for interrupted transfers
  - [ ] Bandwidth throttling for shared networks
  - [ ] Certificate validation for secure transfers

- [ ] **Real-time Sequencer Integration**:
  - [ ] Direct connection to Illumina/ONT sequencers
  - [ ] Live basecalling integration
  - [ ] Real-time quality metrics
  - [ ] Early stopping based on results

#### GUI Feature Implementation
- [ ] **Dynamic DSL Generation**:
  - [ ] Convert GUI selections to BioGPU DSL code
  - [ ] Validation of parameter combinations
  - [ ] Template system for common workflows
  - [ ] Custom DSL editor with syntax highlighting
  ```python
  def generate_dsl_from_gui(config: GUIConfig) -> str:
      template = DSLTemplate()
      template.add_input(config.input_files, config.remote_urls)
      template.add_analysis(config.selected_panels)
      template.add_gpu_settings(config.gpu_device, config.batch_size)
      return template.render()
  ```

- [ ] **Progress Monitoring System**:
  - [ ] Real-time read count updates
  - [ ] Per-stage progress tracking
  - [ ] ETA calculation based on throughput
  - [ ] Resource usage monitoring (GPU, memory)
  - [ ] Error handling with retry options

- [ ] **Results Visualization**:
  - [ ] Interactive microbiome composition charts (Plotly)
  - [ ] Resistance gene heatmaps with clustering
  - [ ] Phylogenetic trees for strain identification
  - [ ] Clinical interpretation summary
  - [ ] Export to PDF/HTML reports

#### Security and Compliance Features
- [ ] **HIPAA Compliance**:
  - [ ] End-to-end encryption for data transfer
  - [ ] Audit logging for all operations
  - [ ] Role-based access control (RBAC)
  - [ ] Automatic PHI de-identification option
  ```python
  class HIPAACompliantGUI:
      def __init__(self):
          self.encryption = AES256Encryption()
          self.audit_logger = ComplianceLogger()
          self.access_control = RBACManager()
      
      def process_patient_data(self, data, user):
          self.audit_logger.log(user, "data_access", data.metadata)
          encrypted = self.encryption.encrypt(data)
          return self.pipeline.process(encrypted)
  ```

- [ ] **Clinical Integration**:
  - [ ] HL7/FHIR message generation
  - [ ] LIS system connectors
  - [ ] Automated result routing
  - [ ] Clinical decision support alerts

#### Performance Optimizations
- [ ] **Efficient Data Transfer**:
  - [ ] Zero-copy streaming to GPU
  - [ ] Compressed transfer formats
  - [ ] Parallel upload/download streams
  - [ ] Smart caching for repeated analyses

- [ ] **Responsive UI Design**:
  - [ ] Non-blocking UI updates
  - [ ] Progressive result loading
  - [ ] Cancelable operations
  - [ ] Background job management

## Getting Started

### Quick Start (Production)
```bash
# Clone repository
git clone https://github.com/yourusername/biogpu.git
cd biogpu

# Build the pipeline
./scripts/build_integrated.sh

# Run with clinical data
./build_integrated/fq_pipeline_gpu /data/fq_resistance_index \
    clinical_R1.fastq.gz clinical_R2.fastq.gz \
    resistance_results.json
```

### Development Setup
```bash
# Set up development environment
conda create -n biogpu python=3.10
conda activate biogpu
pip install -r requirements.txt

# Install CUDA toolkit (if not already installed)
# Follow NVIDIA installation guide for your system

# Build and test
mkdir build && cd build
cmake ..
make -j8
make test
```

## Contributing

We welcome contributions! Please see our contributing guidelines and feel free to submit issues or pull requests.

## License

MIT License - see LICENSE.txt for details

## Contact

David Haslam - dbhaslam@interface-labs.com

---

*Accelerating precision medicine through GPU-powered bioinformatics*
