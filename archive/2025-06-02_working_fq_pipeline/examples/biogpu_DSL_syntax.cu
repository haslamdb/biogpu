# BioGPU Domain-Specific Language Syntax
# This shows how users would write analysis pipelines

# =============================================================================
# BASIC LEVOFLOXACIN RESISTANCE DETECTION
# =============================================================================

pipeline LevofloxacinResistanceScreen {
    # Input/output declarations
    input: {
        patient_reads: fastq_file,
        patient_id: string,
        fq_concentration: float optional  # If available from LCMS
    }
    
    output: {
        resistance_report: json,
        abundance_table: csv,
        clinical_summary: pdf
    }
    
    # Reference databases
    references: {
        microbe_genomes: genome_database("refseq_bacteria_complete"),
        fq_mutations: mutation_database("fluoroquinolone_qrdr_v2"),
        gene_sequences: gene_database("resistance_genes_v3")
    }
    
    # Stage 1: Quality control and preprocessing
    @gpu_kernel
    stage quality_control {
        # Built-in QC operations
        filtered_reads = filter_reads(patient_reads) {
            min_quality: 20,
            min_length: 100,
            remove_adapters: true,
            remove_duplicates: true
        }
        
        # Generate QC metrics
        qc_stats = calculate_stats(filtered_reads) {
            metrics: [read_count, quality_distribution, length_distribution, gc_content]
        }
        
        emit: filtered_reads, qc_stats
    }
    
    # Stage 2: Microbial community profiling
    @gpu_kernel(memory: 16GB, threads: 1024)
    stage profile_microbiome {
        requires: filtered_reads from quality_control
        
        # Map reads to reference genomes
        alignments = parallel_map(filtered_reads, microbe_genomes) {
            algorithm: minimap2_gpu,
            preset: "sr",  # short read
            min_identity: 0.95,
            min_alignment_length: 100
        }
        
        # Calculate abundance
        abundance = calculate_abundance(alignments) {
            method: "relative_abundance",
            min_coverage: 0.1,
            min_reads: 10
        }
        
        emit: alignments, abundance
    }
    
    # Stage 3: Targeted resistance detection
    @gpu_kernel(priority: high)
    stage detect_fq_resistance {
        requires: alignments from profile_microbiome
        
        # Extract reads mapping to resistance genes
        target_reads = extract_gene_reads(alignments) {
            genes: ["gyrA", "gyrB", "parC", "parE"],
            include_flanking: 200  # bp on each side
        }
        
        # Scan for known mutations
        mutations = scan_mutations(target_reads, fq_mutations) {
            # Check specific codon positions
            positions: {
                gyrA: [83, 87],     # S83L, D87N common in E. coli
                gyrB: [426, 447],   
                parC: [80, 84],     # S80I, E84V
                parE: [416, 420]
            },
            
            # Allow 1 mismatch for sequencing errors
            max_mismatches: 1,
            min_quality: 30
        }
        
        # Also check for novel variants
        novel = detect_novel_variants(target_reads) {
            reference: gene_sequences,
            report_synonymous: false,
            min_frequency: 0.01  # Report if >1% of reads
        }
        
        emit: mutations, novel
    }
    
    # Stage 4: Link resistance to organisms
    @gpu_kernel
    stage attribute_resistance {
        requires: mutations from detect_fq_resistance,
                 abundance from profile_microbiome
        
        # Connect mutations to source organisms
        attributed = link_mutations_to_species(mutations, alignments) {
            min_confidence: 0.95,
            handle_ambiguous: "proportional"  # Split by abundance
        }
        
        # Calculate per-organism resistance scores
        organism_resistance = calculate_resistance_profile(attributed) {
            weight_by_abundance: true,
            aggregate_mutations: "maximum"  # Use highest resistance
        }
        
        emit: organism_resistance
    }
    
    # Stage 5: Generate clinical report
    stage create_report {
        requires: organism_resistance from attribute_resistance,
                 abundance from profile_microbiome,
                 qc_stats from quality_control
        
        # Build structured report
        report = clinical_report_builder() {
            sections: [
                qc_summary(qc_stats),
                microbiome_overview(abundance),
                resistance_findings(organism_resistance),
                treatment_guidance(organism_resistance)
            ]
        }
        
        # Add interpretation
        interpret = add_clinical_interpretation(report) {
            threshold_rules: {
                high_risk: "any organism >10% abundance with FQ resistance",
                moderate_risk: "total FQ resistant organisms >20%",
                low_risk: "no FQ resistance detected above 5%"
            }
        }
        
        write: report to resistance_report,
               abundance to abundance_table,
               interpret to clinical_summary
    }
}

# =============================================================================
# LONGITUDINAL RESISTANCE TRACKING
# =============================================================================

pipeline LongitudinalFQTracking {
    input: {
        patient_id: string,
        timepoints: array<{
            sample_id: string,
            reads: fastq_file,
            collection_date: date,
            fq_concentration: float,  # From LCMS
            days_on_therapy: int
        }>
    }
    
    # Process each timepoint
    @parallel
    foreach timepoint in timepoints {
        # Run standard resistance detection
        resistance = LevofloxacinResistanceScreen(
            patient_reads: timepoint.reads,
            patient_id: patient_id,
            fq_concentration: timepoint.fq_concentration
        )
        
        collect: resistance.organism_resistance as timepoint_results
    }
    
    # Analyze resistance evolution
    @gpu_kernel
    stage analyze_evolution {
        # Track mutation frequencies over time
        mutation_trajectories = track_allele_frequencies(timepoint_results) {
            mutations_of_interest: [
                "gyrA_S83L", "gyrA_D87N", 
                "parC_S80I", "parC_E84V"
            ],
            smoothing: "loess",
            confidence_intervals: true
        }
        
        # Correlate with drug levels
        selection_analysis = correlate_with_drug(
            trajectories: mutation_trajectories,
            drug_levels: timepoints.fq_concentration
        ) {
            model: "logistic_growth",
            estimate_selection_coefficient: true
        }
        
        # Identify emerging resistance
        emergence_events = detect_emergence(mutation_trajectories) {
            threshold: 0.01,  # 1% frequency
            min_fold_change: 10
        }
        
        emit: selection_analysis, emergence_events
    }
    
    # Create longitudinal report
    stage longitudinal_report {
        visualizations = create_plots() {
            resistance_timeline: stacked_area_plot(mutation_trajectories),
            selection_pressure: scatter_plot(fq_concentration vs resistance_fraction),
            emergence_dynamics: line_plot(mutation_frequencies over time)
        }
        
        summary = summarize_findings() {
            key_findings: [
                time_to_resistance_emergence,
                dominant_resistance_mechanism,
                correlation_with_drug_levels,
                clinical_implications
            ]
        }
    }
}

# =============================================================================
# BATCH PROCESSING FOR RESEARCH
# =============================================================================

pipeline FQResistanceStudy {
    input: {
        sample_manifest: csv_file,  # Columns: sample_id, reads_file, metadata
        output_dir: directory
    }
    
    # Load and validate samples
    stage load_samples {
        samples = parse_manifest(sample_manifest)
        validate_files(samples.reads_file)
    }
    
    # Process all samples in parallel
    @parallel(max_concurrent: 10)
    foreach sample in samples {
        result = LevofloxacinResistanceScreen(
            patient_reads: sample.reads_file,
            patient_id: sample.sample_id
        )
        
        # Store results
        save: result to "{output_dir}/{sample.sample_id}/"
    }
    
    # Aggregate results for analysis
    @gpu_kernel
    stage aggregate_results {
        # Combine all resistance profiles
        all_resistance = merge_resistance_profiles(results) {
            group_by: sample.metadata.patient_group
        }
        
        # Statistical analysis
        statistics = analyze_resistance_patterns(all_resistance) {
            tests: [
                prevalence_by_group,
                mutation_associations,
                diversity_metrics
            ],
            multiple_testing_correction: "fdr"
        }
        
        # Machine learning analysis
        ml_results = train_resistance_predictor(all_resistance) {
            features: [
                microbiome_composition,
                patient_metadata,
                prior_antibiotic_exposure
            ],
            model: "random_forest",
            cross_validation: 5
        }
    }
}

# =============================================================================
# CUSTOM EXTENSIONS
# =============================================================================

# Users can define custom mutation databases
mutation_database custom_mutations {
    format: "tsv",
    columns: {
        gene: string,
        position: int,
        wild_type: string,
        mutant: string,
        drug: string,
        mic_change: float,
        evidence: string
    },
    
    validation_rules: {
        position: "must be valid codon number",
        wild_type: "must be valid amino acid",
        mic_change: "must be positive number"
    }
}

# Custom algorithms can be plugged in
algorithm custom_aligner implements SequenceAligner {
    parameters: {
        match_score: 2,
        mismatch_penalty: -1,
        gap_open: -3,
        gap_extend: -1
    }
    
    # Implementation would be in CUDA/C++
    implementation: "custom_aligner.cu"
}

# Extend existing pipelines
extend LevofloxacinResistanceScreen {
    # Add plasmid detection
    stage detect_plasmids after profile_microbiome {
        plasmids = find_circular_contigs(filtered_reads) {
            min_length: 1000,
            min_coverage: 10
        }
        
        # Check for plasmid-mediated quinolone resistance
        pmqr = scan_genes(plasmids, ["qnrA", "qnrB", "qnrS", "aac(6')-Ib-cr"])
    }
}

# =============================================================================
# CONFIGURATION AND DEPLOYMENT
# =============================================================================

configuration production {
    gpu_settings: {
        device_id: 0,
        memory_limit: "32GB",
        kernel_timeout: 300  # seconds
    }
    
    error_handling: {
        on_gpu_memory_error: "split_batch_and_retry",
        on_kernel_timeout: "fall_back_to_cpu",
        max_retries: 3
    }
    
    performance: {
        batch_size: "auto",  # Automatically optimize
        prefetch_batches: 2,
        enable_profiling: false
    }
}

# Deploy as service
service resistance_api {
    pipeline: LevofloxacinResistanceScreen,
    configuration: production,
    
    endpoints: {
        "/analyze": {
            method: "POST",
            accepts: "fastq",
            returns: "json"
        },
        
        "/batch": {
            method: "POST", 
            accepts: "manifest",
            returns: "zip"
        }
    },
    
    authentication: "api_key",
    rate_limit: "10 requests per minute"
}