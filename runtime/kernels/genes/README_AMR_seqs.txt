
# these peptide and nucleotide seqences come from AMRFinder: 
# https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa
# https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS.fa
# https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt
# note the AMRProt-mutation table, which we will probably move to FQ resistance workflow
-rw-rw-r-- 1 david david  2183397 Jun  5 14:58 ReferenceGeneCatalog.txt
-rw-rw-r-- 1 david david   179598 Jun  5 14:58 AMRProt-mutation.tsv
-rw-rw-r-- 1 david david  4484691 Jun  5 14:58 AMRProt.fa
-rw-rw-r-- 1 david david 10611888 Jun  5 14:58 AMR_CDS.fa

# these come from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047
-rw-rw-r-- 1 david david  9147398 Jul  1 15:17 amr_nucl.fasta
-rw-rw-r-- 1 david david  3412761 Jul  1 15:17 amr_prot.fasta


# we can download the AMRFinderPlus data using this script:
/home/david/Documents/Code/biogpu/scripts/download_amrfinderplus_db.sh

# then manually renamed AMRProt.fa to AMR_protein.fa and moved to /data directory