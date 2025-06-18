#!/usr/bin/env python3
"""
Improved download script with better rate limiting and progress tracking
"""

import os
import sys
import time
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez, SeqIO
import argparse
import traceback
from datetime import datetime

# Priority pathogens for pediatric infectious diseases
PRIORITY_PATHOGENS = {
    # UTI pathogens
    "Escherichia coli": {"taxid": "562", "priority": 1},
    "Klebsiella pneumoniae": {"taxid": "573", "priority": 1},
    "Proteus mirabilis": {"taxid": "584", "priority": 2},
    "Enterococcus faecalis": {"taxid": "1351", "priority": 2},
    "Pseudomonas aeruginosa": {"taxid": "287", "priority": 2},
    
    # Respiratory pathogens
    "Streptococcus pneumoniae": {"taxid": "1313", "priority": 1},
    "Haemophilus influenzae": {"taxid": "727", "priority": 1},
    "Moraxella catarrhalis": {"taxid": "480", "priority": 2},
    "Staphylococcus aureus": {"taxid": "1280", "priority": 1},
    
    # GI pathogens
    "Clostridioides difficile": {"taxid": "1496", "priority": 1},
    "Salmonella enterica": {"taxid": "28901", "priority": 1},
    "Campylobacter jejuni": {"taxid": "197", "priority": 2},
    "Shigella flexneri": {"taxid": "623", "priority": 2},
    
    # Bloodstream/sepsis
    "Streptococcus agalactiae": {"taxid": "1311", "priority": 1},
    "Listeria monocytogenes": {"taxid": "1639", "priority": 2},
    "Neisseria meningitidis": {"taxid": "487", "priority": 2},
    
    # Other common pathogens
    "Enterobacter cloacae": {"taxid": "550", "priority": 2},
    "Citrobacter freundii": {"taxid": "546", "priority": 3},
    "Serratia marcescens": {"taxid": "615", "priority": 3},
    "Acinetobacter baumannii": {"taxid": "470", "priority": 2},
}

class ProgressTracker:
    def __init__(self, output_dir):
        self.progress_file = os.path.join(output_dir, "download_progress.json")
        self.load_progress()
    
    def load_progress(self):
        if os.path.exists(self.progress_file):
            with open(self.progress_file, 'r') as f:
                self.progress = json.load(f)
        else:
            self.progress = {
                "downloaded_genomes": {},
                "failed_downloads": {},
                "species_status": {},
                "last_updated": None
            }
    
    def save_progress(self):
        self.progress["last_updated"] = datetime.now().isoformat()
        with open(self.progress_file, 'w') as f:
            json.dump(self.progress, f, indent=2)
    
    def mark_downloaded(self, species, genome_id, file_path):
        if species not in self.progress["downloaded_genomes"]:
            self.progress["downloaded_genomes"][species] = []
        self.progress["downloaded_genomes"][species].append({
            "genome_id": genome_id,
            "file": file_path,
            "timestamp": datetime.now().isoformat()
        })
        self.save_progress()
    
    def mark_failed(self, species, genome_id, error):
        if species not in self.progress["failed_downloads"]:
            self.progress["failed_downloads"][species] = []
        self.progress["failed_downloads"][species].append({
            "genome_id": genome_id,
            "error": str(error),
            "timestamp": datetime.now().isoformat()
        })
        self.save_progress()
    
    def is_downloaded(self, species, genome_id):
        if species in self.progress["downloaded_genomes"]:
            for genome in self.progress["downloaded_genomes"][species]:
                if genome["genome_id"] == genome_id:
                    return True
        return False
    
    def get_species_count(self, species):
        if species in self.progress["downloaded_genomes"]:
            return len(self.progress["downloaded_genomes"][species])
        return 0

def retry_with_backoff(func, max_retries=5, initial_delay=1):
    """Retry a function with exponential backoff"""
    delay = initial_delay
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            if "429" in str(e) or "Too Many Requests" in str(e):
                if attempt < max_retries - 1:
                    print(f"[RATE LIMIT] Waiting {delay}s before retry (attempt {attempt+1}/{max_retries})")
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                else:
                    raise
            else:
                raise
    return None

def download_genome_list(taxid, species_name, limit=10):
    """Get list of RefSeq genome assemblies for a given taxon"""
    print(f"\n[SEARCH] {species_name} (taxid {taxid})...")

    # def search_genomes():
    # search_term = f"txid{taxid}[Organism:exp] AND (complete genome[Title] OR chromosome[Title]) AND RefSeq[Filter]"
    # handle = Entrez.esearch(db="nuccore", term=search_term, retmax=limit)
    # record = Entrez.read(handle)
    # handle.close()
    # return record["IdList"]
    
    # Use a more specific search for complete genomes
    def search_genomes():
        search_term = f"txid{taxid}[Organism:exp] AND (complete genome[Title] AND RefSeq[Filter]"
        handle = Entrez.esearch(db="nuccore", term=search_term, retmax=limit)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    
    try:
        ids = retry_with_backoff(search_genomes)
        print(f"[RESULT] Found {len(ids)} genome records for {species_name}")
        
        if len(ids) == 0:
            print(f"[WARNING] No complete genomes found, trying broader search...")
            def search_broader():
                search_term = f"txid{taxid}[Organism:exp] AND RefSeq[Filter]"
                handle = Entrez.esearch(db="nuccore", term=search_term, retmax=limit)
                record = Entrez.read(handle)
                handle.close()
                return record["IdList"]
            
            ids = retry_with_backoff(search_broader)
            print(f"[RESULT] Broader search found {len(ids)} records")
        
        return ids[:limit]
        
    except Exception as e:
        print(f"[ERROR] Failed to search for {species_name}: {e}")
        return []

def download_genome_fasta(genome_id, output_dir, species_name, tracker):
    """Download a single genome in FASTA format"""
    safe_name = species_name.replace(" ", "_").replace("/", "_")
    output_file = os.path.join(output_dir, f"{safe_name}_{genome_id}.fasta")
    
    # Check if already downloaded
    if tracker.is_downloaded(species_name, genome_id) or os.path.exists(output_file):
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"[SKIP] Already exists: {safe_name}_{genome_id} ({file_size:,} bytes)")
            if not tracker.is_downloaded(species_name, genome_id):
                tracker.mark_downloaded(species_name, genome_id, output_file)
        return output_file
    
    def download():
        handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="fasta", retmode="text")
        content = handle.read()
        handle.close()
        
        with open(output_file, 'w') as f:
            f.write(content)
        
        return output_file
    
    try:
        print(f"[DOWNLOAD] {species_name} - {genome_id}")
        result = retry_with_backoff(download, max_retries=3, initial_delay=2)
        
        if result and os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"[SUCCESS] Downloaded {safe_name}_{genome_id} ({file_size:,} bytes)")
            tracker.mark_downloaded(species_name, genome_id, output_file)
            return output_file
        else:
            raise Exception("Download failed")
            
    except Exception as e:
        print(f"[ERROR] Failed {genome_id}: {str(e)}")
        tracker.mark_failed(species_name, genome_id, str(e))
        if os.path.exists(output_file):
            os.remove(output_file)
        return None

def extract_kmers_from_genome(fasta_file, k=31, stride=1):
    """Extract all k-mers from a genome FASTA file"""
    kmers = set()
    
    try:
        seq_count = 0
        total_length = 0
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_count += 1
            seq = str(record.seq).upper()
            total_length += len(seq)
            
            # Extract k-mers with stride
            for i in range(0, len(seq) - k + 1, stride):
                kmer = seq[i:i+k]
                
                # Skip k-mers with N or other ambiguous bases
                if all(base in 'ACGT' for base in kmer):
                    # Add canonical k-mer
                    rc = reverse_complement(kmer)
                    canonical = min(kmer, rc)
                    kmers.add(canonical)
        
        return kmers
        
    except Exception as e:
        print(f"[ERROR] Failed to process {fasta_file}: {e}")
        return set()

def reverse_complement(seq):
    """Get reverse complement of a sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def print_status(tracker):
    """Print current download status"""
    print("\n" + "="*70)
    print("DOWNLOAD STATUS")
    print("="*70)
    
    total_downloaded = 0
    total_failed = 0
    
    for species in sorted(PRIORITY_PATHOGENS.keys()):
        downloaded = tracker.get_species_count(species)
        failed = len(tracker.progress["failed_downloads"].get(species, []))
        total_downloaded += downloaded
        total_failed += failed
        
        status = f"{downloaded} downloaded"
        if failed > 0:
            status += f", {failed} failed"
        
        print(f"{species:<40} {status}")
    
    print(f"\nTOTAL: {total_downloaded} downloaded, {total_failed} failed")
    print("="*70 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Download microbial genomes with improved rate limiting")
    parser.add_argument("output_dir", help="Output directory for genomes")
    parser.add_argument("--genomes-per-species", type=int, default=5)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--email", required=True)
    parser.add_argument("--delay", type=float, default=1.0, help="Delay between downloads (seconds)")
    parser.add_argument("--status-only", action="store_true", help="Just show download status")
    args = parser.parse_args()
    
    Entrez.email = args.email
    
    # Create output directories
    genome_dir = os.path.join(args.output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # Initialize progress tracker
    tracker = ProgressTracker(args.output_dir)
    
    if args.status_only:
        print_status(tracker)
        return
    
    print(f"[CONFIG] Settings:")
    print(f"  - Output: {args.output_dir}")
    print(f"  - Genomes per species: {args.genomes_per_species}")
    print(f"  - K-mer size: {args.k}, Stride: {args.stride}")
    print(f"  - Download delay: {args.delay}s")
    
    # Show current status
    print_status(tracker)
    
    # Collect download tasks
    download_tasks = []
    
    print("[PHASE 1] Searching for genomes...")
    for species, info in PRIORITY_PATHOGENS.items():
        # Skip if already have enough genomes
        current_count = tracker.get_species_count(species)
        if current_count >= args.genomes_per_species:
            print(f"[SKIP] {species} - already have {current_count} genomes")
            continue
        
        taxid = info["taxid"]
        needed = args.genomes_per_species - current_count
        genome_ids = download_genome_list(taxid, species, limit=needed)
        
        for genome_id in genome_ids:
            if not tracker.is_downloaded(species, genome_id):
                download_tasks.append((genome_id, species, taxid))
        
        time.sleep(0.5)  # Rate limit searches
    
    print(f"\n[PHASE 2] Downloading {len(download_tasks)} genomes...")
    
    # Sequential downloads with rate limiting
    for i, (genome_id, species, taxid) in enumerate(download_tasks):
        print(f"\n[{i+1}/{len(download_tasks)}]", end=" ")
        download_genome_fasta(genome_id, genome_dir, species, tracker)
        time.sleep(args.delay)  # Rate limit
    
    # Show final status
    print_status(tracker)
    
    # Process k-mers from all downloaded genomes
    print("\n[PHASE 3] Extracting k-mers...")
    
    all_downloads = []
    for species, genomes in tracker.progress["downloaded_genomes"].items():
        taxid = PRIORITY_PATHOGENS[species]["taxid"]
        for genome in genomes:
            if os.path.exists(genome["file"]):
                all_downloads.append({
                    'species': species,
                    'taxid': taxid,
                    'genome_id': genome["genome_id"],
                    'file': genome["file"]
                })
    
    print(f"Processing {len(all_downloads)} genomes...")
    
    kmer_file = os.path.join(args.output_dir, "database_kmers.txt")
    species_summary = {}
    
    with open(kmer_file, 'w') as f:
        f.write("# K-mer database for BioGPU microbiome profiler\n")
        f.write(f"# K={args.k}, Stride={args.stride}\n")
        f.write("# Format: kmer<tab>taxon_id\n\n")
        
        for i, download in enumerate(all_downloads):
            print(f"[{i+1}/{len(all_downloads)}] {download['species']} - {os.path.basename(download['file'])}")
            
            kmers = extract_kmers_from_genome(download['file'], k=args.k, stride=args.stride)
            
            species = download['species']
            if species not in species_summary:
                species_summary[species] = {'taxid': download['taxid'], 'kmers': set(), 'genomes': 0}
            
            species_summary[species]['kmers'].update(kmers)
            species_summary[species]['genomes'] += 1
            
            # Write k-mers
            for kmer in kmers:
                f.write(f"{kmer}\t{download['taxid']}\n")
    
    # Create summary
    summary_file = os.path.join(args.output_dir, "database_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("K-mer Database Summary\n")
        f.write("="*70 + "\n\n")
        
        total_kmers = 0
        for species, data in sorted(species_summary.items()):
            num_kmers = len(data['kmers'])
            total_kmers += num_kmers
            f.write(f"{species:<40} (taxid {data['taxid']:>6}): {num_kmers:>10,} k-mers ({data['genomes']} genomes)\n")
        
        f.write(f"\nTotal unique k-mers: {total_kmers:,}\n")
    
    print(f"\n[COMPLETE] Database prepared!")
    print(f"  - K-mer file: {kmer_file}")
    print(f"  - Summary: {summary_file}")

if __name__ == "__main__":
    main()