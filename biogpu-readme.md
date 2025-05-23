# BioGPU

A domain-specific programming language for GPU-accelerated DNA sequence mapping to microbial genome databases.

## Overview

BioGPU provides native GPU support for bioinformatics operations with automatic optimization for sequence alignment algorithms. It's designed to accelerate microbial genome analysis and metagenomics workflows by leveraging modern GPU architectures.

### Key Features

- **Native GPU acceleration** for sequence alignment algorithms
- **Domain-specific syntax** designed for bioinformatics workflows
- **Automatic optimization** of memory access patterns for biological data
- **Built-in support** for common formats (FASTA, FASTQ, SAM/BAM)
- **Parallel primitives** optimized for genomic operations

## Example

```biogpu
pipeline MapMicrobiome {
    input: fastq_file reads,
    reference: genome_db microbes,
    output: abundance_table results
    
    @gpu_kernel
    stage align {
        parallel_map(reads, microbes) {
            algorithm: bwa_mem,
            min_score: 30,
            gpu_config: auto
        }
    }
}
```

## Project Status

🚧 **Under Active Development** - This project is in early prototype phase.

### Current Progress
- [x] Project architecture design
- [x] Language syntax specification
- [ ] Basic lexer/parser implementation
- [ ] CUDA kernel prototypes
- [ ] Standard library development

## Installation

### Prerequisites
- CUDA Toolkit 12.0+
- LLVM 15+
- CMake 3.20+
- Python 3.8+ (for prototype)

### Building from Source

```bash
git clone https://github.com/dbhaslam/biogpu.git
cd biogpu
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Documentation

Comprehensive documentation is being developed. See the [docs/](docs/) directory for:
- [Language Specification](docs/language-spec.md)
- [Getting Started Guide](docs/getting-started.md)
- [GPU Optimization Guide](docs/gpu-optimization.md)
- [API Reference](docs/api-reference.md)

## Performance

Target performance improvements over CPU implementations:
- 20x speedup for exact matching
- 10x speedup for approximate matching
- Real-time analysis of microbial communities

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas where we especially need help:
- GPU optimization experts
- Bioinformaticians with alignment algorithm expertise
- Compiler engineers
- Documentation and tutorials
- Testing with real-world datasets

## Use Cases

- Rapid pathogen identification in clinical samples
- Real-time microbiome analysis
- Large-scale metagenomic studies
- Antimicrobial resistance gene detection
- Viral genome surveillance

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Contact

- **Project Lead**: [Your Name]
- **Email**: [your.email@example.com]
- **Issues**: [GitHub Issues](https://github.com/dbhaslam/biogpu/issues)

## Acknowledgments

This project builds upon concepts from:
- NVIDIA NVBIO
- BWA-MEM
- CUSHAW
- The bioinformatics community

---

*Accelerating genomics with the power of GPUs* 🧬⚡