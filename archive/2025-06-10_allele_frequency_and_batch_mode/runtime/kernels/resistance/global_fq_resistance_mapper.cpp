// global_fq_resistance_mapper.cpp
// Implementation of the unified global FQ resistance mapping system

#include "global_fq_resistance_mapper.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cuda_runtime.h>

GlobalFQResistanceMapper::GlobalFQResistanceMapper() 
    : gpu_mutations(nullptr), num_gpu_mutations(0), gpu_data_ready(false) {
    initializeNormalizationMaps();
}

GlobalFQResistanceMapper::~GlobalFQResistanceMapper() {
    if (gpu_mutations) {
        cudaFree(gpu_mutations);
    }
}

void GlobalFQResistanceMapper::initializeNormalizationMaps() {
    // Species name normalization (handle both spaces and underscores)
    species_normalization["Acinetobacter baumannii"] = "Acinetobacter_baumannii";
    species_normalization["Acinetobacter_baumannii"] = "Acinetobacter_baumannii";
    species_normalization["Burkholderia cenocepacia"] = "Burkholderia_cenocepacia";
    species_normalization["Burkholderia_cenocepacia"] = "Burkholderia_cenocepacia";
    species_normalization["Burkholderia dolosa"] = "Burkholderia_dolosa";
    species_normalization["Burkholderia_dolosa"] = "Burkholderia_dolosa";
    species_normalization["Campylobacter coli"] = "Campylobacter_coli";
    species_normalization["Campylobacter_coli"] = "Campylobacter_coli";
    species_normalization["Campylobacter jejuni"] = "Campylobacter_jejuni";
    species_normalization["Campylobacter_jejuni"] = "Campylobacter_jejuni";
    species_normalization["Campylobacter lari"] = "Campylobacter_lari";
    species_normalization["Campylobacter_lari"] = "Campylobacter_lari";
    species_normalization["Campylobacter upsaliensis"] = "Campylobacter_upsaliensis";
    species_normalization["Campylobacter_upsaliensis"] = "Campylobacter_upsaliensis";
    species_normalization["Citrobacter farmeri"] = "Citrobacter_farmeri";
    species_normalization["Citrobacter_farmeri"] = "Citrobacter_farmeri";
    species_normalization["Citrobacter freundii"] = "Citrobacter_freundii";
    species_normalization["Citrobacter_freundii"] = "Citrobacter_freundii";
    species_normalization["Clostridioides difficile"] = "Clostridioides_difficile";
    species_normalization["Clostridioides_difficile"] = "Clostridioides_difficile";
    species_normalization["Enterobacter cloacae"] = "Enterobacter_cloacae";
    species_normalization["Enterobacter_cloacae"] = "Enterobacter_cloacae";
    species_normalization["Enterococcus faecalis"] = "Enterococcus_faecalis";
    species_normalization["Enterococcus_faecalis"] = "Enterococcus_faecalis";
    species_normalization["Enterococcus faecium"] = "Enterococcus_faecium";
    species_normalization["Enterococcus_faecium"] = "Enterococcus_faecium";
    species_normalization["Escherichia coli"] = "Escherichia_coli";
    species_normalization["Escherichia_coli"] = "Escherichia_coli";
    species_normalization["Escherichia fergusonii"] = "Escherichia_fergusonii";
    species_normalization["Escherichia_fergusonii"] = "Escherichia_fergusonii";
    species_normalization["Haemophilus influenzae"] = "Haemophilus_influenzae";
    species_normalization["Haemophilus_influenzae"] = "Haemophilus_influenzae";
    species_normalization["Klebsiella aerogenes"] = "Klebsiella_aerogenes";
    species_normalization["Klebsiella_aerogenes"] = "Klebsiella_aerogenes";
    species_normalization["Klebsiella michiganensis"] = "Klebsiella_michiganensis";
    species_normalization["Klebsiella_michiganensis"] = "Klebsiella_michiganensis";
    species_normalization["Klebsiella oxytoca"] = "Klebsiella_oxytoca";
    species_normalization["Klebsiella_oxytoca"] = "Klebsiella_oxytoca";
    species_normalization["Klebsiella pneumoniae"] = "Klebsiella_pneumoniae";
    species_normalization["Klebsiella_pneumoniae"] = "Klebsiella_pneumoniae";
    species_normalization["Neisseria gonorrhoeae"] = "Neisseria_gonorrhoeae";
    species_normalization["Neisseria_gonorrhoeae"] = "Neisseria_gonorrhoeae";
    species_normalization["Neisseria meningitidis"] = "Neisseria_meningitidis";
    species_normalization["Neisseria_meningitidis"] = "Neisseria_meningitidis";
    species_normalization["Pseudomonas aeruginosa"] = "Pseudomonas_aeruginosa";
    species_normalization["Pseudomonas_aeruginosa"] = "Pseudomonas_aeruginosa";
    species_normalization["Raoultella ornithinolytica"] = "Raoultella_ornithinolytica";
    species_normalization["Raoultella_ornithinolytica"] = "Raoultella_ornithinolytica";
    species_normalization["Salmonella enterica"] = "Salmonella_enterica";
    species_normalization["Salmonella_enterica"] = "Salmonella_enterica";
    species_normalization["Serratia marcescens"] = "Serratia_marcescens";
    species_normalization["Serratia_marcescens"] = "Serratia_marcescens";
    species_normalization["Shigella flexneri"] = "Shigella_flexneri";
    species_normalization["Shigella_flexneri"] = "Shigella_flexneri";
    species_normalization["Staphylococcus aureus"] = "Staphylococcus_aureus";
    species_normalization["Staphylococcus_aureus"] = "Staphylococcus_aureus";
    species_normalization["Staphylococcus pseudintermedius"] = "Staphylococcus_pseudintermedius";
    species_normalization["Staphylococcus_pseudintermedius"] = "Staphylococcus_pseudintermedius";
    species_normalization["Streptococcus pneumoniae"] = "Streptococcus_pneumoniae";
    species_normalization["Streptococcus_pneumoniae"] = "Streptococcus_pneumoniae";
    species_normalization["Vibrio cholerae"] = "Vibrio_cholerae";
    species_normalization["Vibrio_cholerae"] = "Vibrio_cholerae";
    species_normalization["Vibrio parahaemolyticus"] = "Vibrio_parahaemolyticus";
    species_normalization["Vibrio_parahaemolyticus"] = "Vibrio_parahaemolyticus";
    species_normalization["Vibrio vulnificus"] = "Vibrio_vulnificus";
    species_normalization["Vibrio_vulnificus"] = "Vibrio_vulnificus";
    
    // Gene name normalization (handle case variations)
    gene_normalization["gyrA"] = "gyrA";
    gene_normalization["GyrA"] = "gyrA";
    gene_normalization["GYRA"] = "gyrA";
    gene_normalization["gyrB"] = "gyrB";
    gene_normalization["GyrB"] = "gyrB";
    gene_normalization["GYRB"] = "gyrB";
    gene_normalization["parC"] = "parC";
    gene_normalization["ParC"] = "parC";
    gene_normalization["PARC"] = "parC";
    gene_normalization["parE"] = "parE";
    gene_normalization["ParE"] = "parE";
    gene_normalization["PARE"] = "parE";
    gene_normalization["grlA"] = "grlA";
    gene_normalization["GrlA"] = "grlA";
    gene_normalization["GRLA"] = "grlA";
    gene_normalization["emrR"] = "emrR";
    gene_normalization["EmrR"] = "emrR";
    gene_normalization["EMRR"] = "emrR";
    gene_normalization["nfxB"] = "nfxB";
    gene_normalization["NfxB"] = "nfxB";
    gene_normalization["NFXB"] = "nfxB";
    gene_normalization["qnrB19"] = "qnrB19";
    gene_normalization["QnrB19"] = "qnrB19";
    gene_normalization["QNRB19"] = "qnrB19";
}

std::string GlobalFQResistanceMapper::normalizeSpeciesName(const std::string& species) {
    auto it = species_normalization.find(species);
    if (it != species_normalization.end()) {
        return it->second;
    }
    // If not found, try simple replacement of spaces with underscores
    std::string normalized = species;
    std::replace(normalized.begin(), normalized.end(), ' ', '_');
    return normalized;
}

std::string GlobalFQResistanceMapper::normalizeGeneName(const std::string& gene) {
    auto it = gene_normalization.find(gene);
    if (it != gene_normalization.end()) {
        return it->second;
    }
    return gene;  // Return as-is if not found
}

bool GlobalFQResistanceMapper::loadFromCSV(const std::string& csv_path) {
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open FQ resistance CSV file: " << csv_path << std::endl;
        return false;
    }
    
    all_mutations.clear();
    std::string line;
    bool header_skipped = false;
    
    while (std::getline(file, line)) {
        // Skip header
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }
        
        // Parse CSV line
        std::stringstream ss(line);
        std::string species, gene, location_str, wt, mut;
        
        // Read fields (species,gene,location,wt,mut)
        if (!std::getline(ss, species, ',')) continue;
        if (!std::getline(ss, gene, ',')) continue;
        if (!std::getline(ss, location_str, ',')) continue;
        if (!std::getline(ss, wt, ',')) continue;
        if (!std::getline(ss, mut, ',')) continue;
        
        // Skip NA entries
        if (location_str == "NA" || wt == "NA" || mut == "NA") continue;
        
        // Normalize names
        species = normalizeSpeciesName(species);
        gene = normalizeGeneName(gene);
        
        // Parse position
        uint16_t position = 0;
        try {
            position = std::stoi(location_str);
        } catch (...) {
            std::cerr << "WARNING: Invalid position '" << location_str << "' in CSV" << std::endl;
            continue;
        }
        
        // Create mutation entry
        FQResistanceMutation mutation;
        mutation.species = species;
        mutation.gene = gene;
        mutation.position = position;
        mutation.wildtype_aa = wt[0];  // Take first character
        mutation.mutant_aa = mut[0];   // Take first character
        
        all_mutations.push_back(mutation);
    }
    
    file.close();
    
    // Build lookup tables
    buildLookupTables();
    
    std::cout << "Loaded " << all_mutations.size() << " FQ resistance mutations from CSV" << std::endl;
    return true;
}

void GlobalFQResistanceMapper::buildLookupTables() {
    mutations_by_position.clear();
    resistant_aas_by_position.clear();
    wildtype_by_position.clear();
    
    for (const auto& mut : all_mutations) {
        std::string pos_key = mut.getPositionKey();
        
        // Add to position lookup
        mutations_by_position[pos_key].push_back(mut);
        
        // Add resistant AA
        resistant_aas_by_position[pos_key].insert(mut.mutant_aa);
        
        // Set wildtype (should be consistent for same position)
        wildtype_by_position[pos_key] = mut.wildtype_aa;
    }
}

bool GlobalFQResistanceMapper::loadDatabaseMappings(const std::string& protein_db_path) {
    // Load species and gene mappings from protein database metadata
    std::string metadata_path = protein_db_path + "/metadata.json";
    std::ifstream metadata_file(metadata_path);
    
    if (!metadata_file.good()) {
        std::cerr << "WARNING: Cannot load protein database metadata from " << metadata_path << std::endl;
        // Try alternative paths
        metadata_path = protein_db_path + "/species_gene_map.json";
        metadata_file.open(metadata_path);
        if (!metadata_file.good()) {
            return false;
        }
    }
    
    // Simple JSON parsing (would use proper JSON library in production)
    std::string content((std::istreambuf_iterator<char>(metadata_file)),
                       std::istreambuf_iterator<char>());
    metadata_file.close();
    
    // Parse species_map
    size_t species_map_start = content.find("\"species_map\"");
    if (species_map_start != std::string::npos) {
        size_t map_start = content.find("{", species_map_start);
        size_t map_end = content.find("}", map_start);
        
        if (map_start != std::string::npos && map_end != std::string::npos) {
            std::string species_section = content.substr(map_start + 1, map_end - map_start - 1);
            
            // Parse each "id": "name" pair
            size_t pos = 0;
            while ((pos = species_section.find("\"", pos)) != std::string::npos) {
                size_t id_start = pos + 1;
                size_t id_end = species_section.find("\"", id_start);
                if (id_end == std::string::npos) break;
                
                std::string id_str = species_section.substr(id_start, id_end - id_start);
                
                size_t name_start = species_section.find("\"", id_end + 1) + 1;
                size_t name_end = species_section.find("\"", name_start);
                if (name_end == std::string::npos) break;
                
                std::string name = species_section.substr(name_start, name_end - name_start);
                
                try {
                    uint32_t id = std::stoi(id_str);
                    std::string normalized_name = normalizeSpeciesName(name);
                    species_to_id[normalized_name] = id;
                    id_to_species[id] = normalized_name;
                } catch (...) {
                    // Skip invalid entries
                }
                
                pos = name_end + 1;
            }
        }
    }
    
    // Parse gene_map
    size_t gene_map_start = content.find("\"gene_map\"");
    if (gene_map_start != std::string::npos) {
        size_t map_start = content.find("{", gene_map_start);
        size_t map_end = content.find("}", map_start);
        
        if (map_start != std::string::npos && map_end != std::string::npos) {
            std::string gene_section = content.substr(map_start + 1, map_end - map_start - 1);
            
            // Parse each "id": "name" pair
            size_t pos = 0;
            while ((pos = gene_section.find("\"", pos)) != std::string::npos) {
                size_t id_start = pos + 1;
                size_t id_end = gene_section.find("\"", id_start);
                if (id_end == std::string::npos) break;
                
                std::string id_str = gene_section.substr(id_start, id_end - id_start);
                
                size_t name_start = gene_section.find("\"", id_end + 1) + 1;
                size_t name_end = gene_section.find("\"", name_start);
                if (name_end == std::string::npos) break;
                
                std::string name = gene_section.substr(name_start, name_end - name_start);
                
                try {
                    uint32_t id = std::stoi(id_str);
                    std::string normalized_name = normalizeGeneName(name);
                    gene_to_id[normalized_name] = id;
                    id_to_gene[id] = normalized_name;
                } catch (...) {
                    // Skip invalid entries
                }
                
                pos = name_end + 1;
            }
        }
    }
    
    // Build combined species-gene mappings
    // This assumes your database uses a combined ID system
    // You may need to adjust this based on your actual ID scheme
    uint32_t combined_id = 0;
    for (const auto& species_pair : species_to_id) {
        for (const auto& gene_pair : gene_to_id) {
            std::string key = species_pair.first + "_" + gene_pair.first;
            species_gene_to_id[key] = combined_id;
            id_to_species_gene[combined_id] = {species_pair.first, gene_pair.first};
            combined_id++;
        }
    }
    
    std::cout << "Loaded database mappings: " << species_to_id.size() << " species, " 
              << gene_to_id.size() << " genes" << std::endl;
    
    return true;
}

bool GlobalFQResistanceMapper::isResistanceMutation(const std::string& species, 
                                                   const std::string& gene,
                                                   uint16_t position, 
                                                   char wildtype_aa, 
                                                   char observed_aa) const {
    std::string norm_species = const_cast<GlobalFQResistanceMapper*>(this)->normalizeSpeciesName(species);
    std::string norm_gene = const_cast<GlobalFQResistanceMapper*>(this)->normalizeGeneName(gene);
    std::string pos_key = norm_species + "_" + norm_gene + "_" + std::to_string(position);
    
    auto it = resistant_aas_by_position.find(pos_key);
    if (it != resistant_aas_by_position.end()) {
        // Check if the observed AA is in the resistant set
        return it->second.find(observed_aa) != it->second.end();
    }
    
    return false;
}

std::vector<char> GlobalFQResistanceMapper::getResistantAAs(const std::string& species,
                                                           const std::string& gene,
                                                           uint16_t position) const {
    std::string norm_species = const_cast<GlobalFQResistanceMapper*>(this)->normalizeSpeciesName(species);
    std::string norm_gene = const_cast<GlobalFQResistanceMapper*>(this)->normalizeGeneName(gene);
    std::string pos_key = norm_species + "_" + norm_gene + "_" + std::to_string(position);
    
    std::vector<char> result;
    auto it = resistant_aas_by_position.find(pos_key);
    if (it != resistant_aas_by_position.end()) {
        result.assign(it->second.begin(), it->second.end());
    }
    
    return result;
}

char GlobalFQResistanceMapper::getWildtypeAA(const std::string& species,
                                            const std::string& gene,
                                            uint16_t position) const {
    std::string norm_species = const_cast<GlobalFQResistanceMapper*>(this)->normalizeSpeciesName(species);
    std::string norm_gene = const_cast<GlobalFQResistanceMapper*>(this)->normalizeGeneName(gene);
    std::string pos_key = norm_species + "_" + norm_gene + "_" + std::to_string(position);
    
    auto it = wildtype_by_position.find(pos_key);
    if (it != wildtype_by_position.end()) {
        return it->second;
    }
    
    return 'X';  // Unknown
}

bool GlobalFQResistanceMapper::isResistancePosition(const std::string& species,
                                                   const std::string& gene,
                                                   uint16_t position) const {
    std::string norm_species = const_cast<GlobalFQResistanceMapper*>(this)->normalizeSpeciesName(species);
    std::string norm_gene = const_cast<GlobalFQResistanceMapper*>(this)->normalizeGeneName(gene);
    std::string pos_key = norm_species + "_" + norm_gene + "_" + std::to_string(position);
    
    return wildtype_by_position.find(pos_key) != wildtype_by_position.end();
}

std::vector<FQResistanceMutation> GlobalFQResistanceMapper::getMutationsForGene(
    const std::string& species, const std::string& gene) const {
    
    std::string norm_species = const_cast<GlobalFQResistanceMapper*>(this)->normalizeSpeciesName(species);
    std::string norm_gene = const_cast<GlobalFQResistanceMapper*>(this)->normalizeGeneName(gene);
    
    std::vector<FQResistanceMutation> result;
    for (const auto& mut : all_mutations) {
        if (mut.species == norm_species && mut.gene == norm_gene) {
            result.push_back(mut);
        }
    }
    
    return result;
}

bool GlobalFQResistanceMapper::isResistanceMutationByID(uint32_t species_gene_id,
                                                       uint16_t position,
                                                       char wildtype_aa,
                                                       char observed_aa) const {
    auto it = id_to_species_gene.find(species_gene_id);
    if (it != id_to_species_gene.end()) {
        return isResistanceMutation(it->second.first, it->second.second, 
                                   position, wildtype_aa, observed_aa);
    }
    return false;
}

std::vector<char> GlobalFQResistanceMapper::getResistantAAsByID(uint32_t species_gene_id,
                                                               uint16_t position) const {
    auto it = id_to_species_gene.find(species_gene_id);
    if (it != id_to_species_gene.end()) {
        return getResistantAAs(it->second.first, it->second.second, position);
    }
    return std::vector<char>();
}

char GlobalFQResistanceMapper::getWildtypeAAByID(uint32_t species_gene_id,
                                                uint16_t position) const {
    auto it = id_to_species_gene.find(species_gene_id);
    if (it != id_to_species_gene.end()) {
        return getWildtypeAA(it->second.first, it->second.second, position);
    }
    return 'X';
}

uint32_t GlobalFQResistanceMapper::getSpeciesGeneID(const std::string& species,
                                                   const std::string& gene) const {
    std::string norm_species = const_cast<GlobalFQResistanceMapper*>(this)->normalizeSpeciesName(species);
    std::string norm_gene = const_cast<GlobalFQResistanceMapper*>(this)->normalizeGeneName(gene);
    std::string key = norm_species + "_" + norm_gene;
    
    auto it = species_gene_to_id.find(key);
    if (it != species_gene_to_id.end()) {
        return it->second;
    }
    
    return UINT32_MAX;  // Invalid ID
}

std::pair<std::string, std::string> GlobalFQResistanceMapper::getSpeciesGeneFromID(
    uint32_t species_gene_id) const {
    
    auto it = id_to_species_gene.find(species_gene_id);
    if (it != id_to_species_gene.end()) {
        return it->second;
    }
    
    return {"", ""};
}

bool GlobalFQResistanceMapper::prepareGPUData() {
    if (gpu_data_ready) return true;
    
    // Group mutations by species-gene-position
    std::map<std::string, std::vector<char>> grouped_mutations;
    std::map<std::string, FQResistanceMutationGPU> gpu_mutation_map;
    
    for (const auto& mut : all_mutations) {
        uint32_t species_gene_id = getSpeciesGeneID(mut.species, mut.gene);
        if (species_gene_id == UINT32_MAX) continue;
        
        std::string key = std::to_string(species_gene_id) + "_" + std::to_string(mut.position);
        
        if (gpu_mutation_map.find(key) == gpu_mutation_map.end()) {
            FQResistanceMutationGPU gpu_mut;
            gpu_mut.species_gene_id = species_gene_id;
            gpu_mut.position = mut.position;
            gpu_mut.wildtype_aa = mut.wildtype_aa;
            gpu_mut.num_mutants = 0;
            gpu_mut.resistance_score = 1.0f;  // Default score
            gpu_mutation_map[key] = gpu_mut;
        }
        
        // Add mutant AA if not already present
        FQResistanceMutationGPU& gpu_mut = gpu_mutation_map[key];
        if (gpu_mut.num_mutants < 8) {
            bool found = false;
            for (int i = 0; i < gpu_mut.num_mutants; i++) {
                if (gpu_mut.mutant_aas[i] == mut.mutant_aa) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                gpu_mut.mutant_aas[gpu_mut.num_mutants++] = mut.mutant_aa;
            }
        }
    }
    
    // Convert map to vector
    std::vector<FQResistanceMutationGPU> gpu_mutations_vec;
    for (const auto& pair : gpu_mutation_map) {
        gpu_mutations_vec.push_back(pair.second);
    }
    
    num_gpu_mutations = gpu_mutations_vec.size();
    
    // Allocate GPU memory
    cudaError_t err = cudaMalloc(&gpu_mutations, num_gpu_mutations * sizeof(FQResistanceMutationGPU));
    if (err != cudaSuccess) {
        std::cerr << "ERROR: Failed to allocate GPU memory for FQ mutations: " 
                  << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    // Copy to GPU
    err = cudaMemcpy(gpu_mutations, gpu_mutations_vec.data(),
                     num_gpu_mutations * sizeof(FQResistanceMutationGPU),
                     cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "ERROR: Failed to copy FQ mutations to GPU: " 
                  << cudaGetErrorString(err) << std::endl;
        cudaFree(gpu_mutations);
        gpu_mutations = nullptr;
        return false;
    }
    
    gpu_data_ready = true;
    std::cout << "Prepared " << num_gpu_mutations << " FQ resistance mutations for GPU" << std::endl;
    
    return true;
}

void GlobalFQResistanceMapper::printSummary() const {
    std::cout << "\n=== FQ Resistance Mapper Summary ===" << std::endl;
    std::cout << "Total mutations loaded: " << all_mutations.size() << std::endl;
    std::cout << "Unique positions: " << wildtype_by_position.size() << std::endl;
    
    // Count by species
    std::map<std::string, int> species_counts;
    for (const auto& mut : all_mutations) {
        species_counts[mut.species]++;
    }
    
    std::cout << "\nMutations by species:" << std::endl;
    for (const auto& pair : species_counts) {
        std::cout << "  " << pair.first << ": " << pair.second << std::endl;
    }
    
    // Count by gene
    std::map<std::string, int> gene_counts;
    for (const auto& mut : all_mutations) {
        gene_counts[mut.gene]++;
    }
    
    std::cout << "\nMutations by gene:" << std::endl;
    for (const auto& pair : gene_counts) {
        std::cout << "  " << pair.first << ": " << pair.second << std::endl;
    }
}

void GlobalFQResistanceMapper::printMutationsForGene(const std::string& species,
                                                   const std::string& gene) const {
    auto mutations = getMutationsForGene(species, gene);
    
    std::cout << "\n=== FQ Resistance Mutations for " << species << " " << gene << " ===" << std::endl;
    
    if (mutations.empty()) {
        std::cout << "No mutations found." << std::endl;
        return;
    }
    
    // Sort by position
    std::sort(mutations.begin(), mutations.end(),
              [](const FQResistanceMutation& a, const FQResistanceMutation& b) {
                  return a.position < b.position;
              });
    
    for (const auto& mut : mutations) {
        std::cout << "  Position " << mut.position << ": " 
                  << mut.wildtype_aa << " -> " << mut.mutant_aa << std::endl;
    }
}

// C interface implementation
extern "C" {
    int init_global_fq_mapper(const char* csv_path, const char* protein_db_path) {
        GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
        
        if (!mapper.loadFromCSV(csv_path)) {
            return -1;
        }
        
        if (protein_db_path && !mapper.loadDatabaseMappings(protein_db_path)) {
            std::cerr << "WARNING: Failed to load database mappings, ID-based queries may not work" << std::endl;
        }
        
        return 0;
    }
    
    bool is_fq_resistance_mutation(const char* species, const char* gene,
                                  uint16_t position, char wildtype, char observed) {
        return GlobalFQResistanceMapper::getInstance().isResistanceMutation(
            species, gene, position, wildtype, observed);
    }
    
    char get_fq_wildtype_aa(const char* species, const char* gene, uint16_t position) {
        return GlobalFQResistanceMapper::getInstance().getWildtypeAA(species, gene, position);
    }
    
    bool is_fq_resistance_mutation_by_id(uint32_t species_gene_id, uint16_t position,
                                        char wildtype, char observed) {
        return GlobalFQResistanceMapper::getInstance().isResistanceMutationByID(
            species_gene_id, position, wildtype, observed);
    }
    
    void* get_fq_gpu_mutations() {
        GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
        mapper.prepareGPUData();
        return mapper.getGPUMutations();
    }
    
    uint32_t get_fq_gpu_mutation_count() {
        return GlobalFQResistanceMapper::getInstance().getNumGPUMutations();
    }
    
    void cleanup_global_fq_mapper() {
        // Singleton will clean up on program exit
    }
}