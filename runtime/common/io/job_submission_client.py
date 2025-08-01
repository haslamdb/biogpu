#!/usr/bin/env python3
"""
BioGPU Job Submission Client
Processes manifest.json files and submits jobs to the streaming service
"""

import json
import uuid
import logging
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass
from kafka import KafkaProducer
from kafka.errors import KafkaError
import argparse
import sys

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class SampleEntry:
    """Represents a sample with paired FASTQ files"""
    sample_id: str
    read1_path: str
    read2_path: str
    metadata: Optional[Dict] = None


class ManifestProcessor:
    """Processes manifest.json files and extracts sample information"""
    
    def __init__(self, manifest_path: str):
        self.manifest_path = Path(manifest_path)
        if not self.manifest_path.exists():
            raise FileNotFoundError(f"Manifest file not found: {manifest_path}")
        
        self.manifest_data = None
        self.samples = []
        
    def load_manifest(self) -> Dict:
        """Load and validate manifest.json"""
        try:
            with open(self.manifest_path, 'r') as f:
                self.manifest_data = json.load(f)
            
            logger.info(f"Loaded manifest with {len(self.manifest_data.get('samples', []))} samples")
            return self.manifest_data
            
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in manifest: {e}")
            raise
    
    def extract_samples(self) -> List[SampleEntry]:
        """Extract sample entries from manifest"""
        if not self.manifest_data:
            self.load_manifest()
        
        samples = []
        
        # Handle different manifest formats
        if 'samples' in self.manifest_data:
            # Format 1: {"samples": [{"sample_id": "...", "read1": "...", "read2": "..."}]}
            for sample in self.manifest_data['samples']:
                entry = SampleEntry(
                    sample_id=sample.get('sample_id', sample.get('id')),
                    read1_path=sample.get('read1', sample.get('r1', sample.get('forward'))),
                    read2_path=sample.get('read2', sample.get('r2', sample.get('reverse'))),
                    metadata=sample.get('metadata', {})
                )
                samples.append(entry)
                
        elif 'files' in self.manifest_data:
            # Format 2: {"files": {"sample1": {"R1": "...", "R2": "..."}, ...}}
            for sample_id, files in self.manifest_data['files'].items():
                entry = SampleEntry(
                    sample_id=sample_id,
                    read1_path=files.get('R1', files.get('r1')),
                    read2_path=files.get('R2', files.get('r2')),
                    metadata=files.get('metadata', {})
                )
                samples.append(entry)
        
        # Validate all samples have required fields
        for sample in samples:
            if not all([sample.sample_id, sample.read1_path, sample.read2_path]):
                logger.warning(f"Incomplete sample entry: {sample}")
                continue
            
            # Check if files exist
            r1_exists = Path(sample.read1_path).exists()
            r2_exists = Path(sample.read2_path).exists()
            
            if not r1_exists or not r2_exists:
                logger.warning(f"Missing files for sample {sample.sample_id}: "
                             f"R1={r1_exists}, R2={r2_exists}")
        
        self.samples = samples
        logger.info(f"Extracted {len(samples)} valid samples")
        return samples


class BioGPUJobSubmitter:
    """Submits jobs to BioGPU streaming service via Kafka"""
    
    def __init__(self, kafka_servers: List[str], output_base_path: str = "/data/output"):
        self.kafka_servers = kafka_servers
        self.output_base_path = Path(output_base_path)
        self.producer = None
        
    def connect(self):
        """Initialize Kafka producer"""
        try:
            self.producer = KafkaProducer(
                bootstrap_servers=self.kafka_servers,
                value_serializer=lambda v: json.dumps(v).encode('utf-8'),
                acks='all',  # Wait for all replicas
                retries=3,
                max_in_flight_requests_per_connection=1
            )
            logger.info(f"Connected to Kafka brokers: {self.kafka_servers}")
        except KafkaError as e:
            logger.error(f"Failed to connect to Kafka: {e}")
            raise
    
    def submit_job(self, sample: SampleEntry, batch_id: Optional[str] = None) -> str:
        """Submit a single job to the processing queue"""
        if not self.producer:
            self.connect()
        
        # Generate unique job ID
        job_id = f"job-{uuid.uuid4()}"
        
        # Create output directory for this sample
        output_path = self.output_base_path / (batch_id or "batch") / sample.sample_id
        
        # Build job data
        job_data = {
            'job_id': job_id,
            'sample_id': sample.sample_id,
            'path_to_r1': sample.read1_path,
            'path_to_r2': sample.read2_path,
            'output_path': str(output_path),
            'metadata': sample.metadata or {},
            'batch_id': batch_id,
            'submission_time': datetime.utcnow().isoformat()
        }
        
        try:
            # Send to Kafka
            future = self.producer.send('new_jobs', job_data)
            
            # Wait for confirmation
            record_metadata = future.get(timeout=10)
            
            logger.info(f"Submitted job {job_id} for sample {sample.sample_id} "
                       f"to partition {record_metadata.partition} "
                       f"at offset {record_metadata.offset}")
            
            return job_id
            
        except KafkaError as e:
            logger.error(f"Failed to submit job for sample {sample.sample_id}: {e}")
            raise
    
    def submit_batch(self, samples: List[SampleEntry]) -> Dict[str, str]:
        """Submit multiple samples as a batch"""
        batch_id = f"batch-{uuid.uuid4()}"
        job_mapping = {}
        
        logger.info(f"Submitting batch {batch_id} with {len(samples)} samples")
        
        for sample in samples:
            try:
                job_id = self.submit_job(sample, batch_id)
                job_mapping[sample.sample_id] = job_id
            except Exception as e:
                logger.error(f"Failed to submit sample {sample.sample_id}: {e}")
                continue
        
        # Flush all pending messages
        if self.producer:
            self.producer.flush()
        
        logger.info(f"Batch submission complete. Submitted {len(job_mapping)}/{len(samples)} jobs")
        return job_mapping
    
    def close(self):
        """Close Kafka producer"""
        if self.producer:
            self.producer.close()


def main():
    """Main entry point for job submission"""
    parser = argparse.ArgumentParser(description='Submit BioGPU jobs from manifest file')
    parser.add_argument(
        'manifest',
        help='Path to manifest.json file'
    )
    parser.add_argument(
        '--kafka-servers',
        default='localhost:9092',
        help='Comma-separated list of Kafka brokers (default: localhost:9092)'
    )
    parser.add_argument(
        '--output-path',
        default='/data/output',
        help='Base path for output files (default: /data/output)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Parse manifest without submitting jobs'
    )
    
    args = parser.parse_args()
    
    try:
        # Process manifest
        processor = ManifestProcessor(args.manifest)
        samples = processor.extract_samples()
        
        if not samples:
            logger.error("No valid samples found in manifest")
            sys.exit(1)
        
        if args.dry_run:
            logger.info("Dry run - would submit the following samples:")
            for sample in samples:
                logger.info(f"  {sample.sample_id}: {sample.read1_path} + {sample.read2_path}")
            return
        
        # Submit jobs
        kafka_servers = args.kafka_servers.split(',')
        submitter = BioGPUJobSubmitter(kafka_servers, args.output_path)
        
        try:
            job_mapping = submitter.submit_batch(samples)
            
            # Save job mapping for reference
            mapping_file = Path(args.manifest).with_suffix('.jobs.json')
            with open(mapping_file, 'w') as f:
                json.dump(job_mapping, f, indent=2)
            
            logger.info(f"Job mapping saved to {mapping_file}")
            
        finally:
            submitter.close()
            
    except Exception as e:
        logger.error(f"Job submission failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    from datetime import datetime
    main()