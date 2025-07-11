// classification_report_generator.h
// Comprehensive report generator for Kraken classification results
// Generates MetaPhlAn-style, species counts, and summary reports

#ifndef CLASSIFICATION_REPORT_GENERATOR_H
#define CLASSIFICATION_REPORT_GENERATOR_H

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <filesystem>

// Structure to hold classification result from Kraken output
struct KrakenResult {
    char status;           // 'C' = classified, 'U' = unclassified
    std::string read_id;   // Read identifier
    uint32_t taxon_id;     // Taxonomic ID (0 if unclassified)
    std::string length;    // Read length(s) - may be "150" or "150|150" for paired
    float confidence;      // Classification confidence score
    std::string votes;     // Vote counts - may be "45" or "45|50" for paired
    std::string kmers;     // K-mer counts - may be "50" or "50|55" for paired
    
    bool is_paired = false;
    int read1_length = 0, read2_length = 0;
    int read1_votes = 0, read2_votes = 0;
    int read1_kmers = 0, read2_kmers = 0;
};

// Taxonomic information for a taxon
struct TaxonInfo {
    uint32_t taxon_id;
    std::string name;
    uint32_t parent_id;
    std::string rank;
    int read_count = 0;
    float abundance = 0.0f;
    
    // For building hierarchical paths
    std::vector<std::string> lineage;
    std::string full_lineage;
};

// Classification summary statistics
struct ClassificationStats {
    int total_reads = 0;
    int classified_reads = 0;
    int unclassified_reads = 0;
    float classification_rate = 0.0f;
    
    // Per-rank statistics
    std::map<std::string, int> reads_per_rank;
    std::map<std::string, int> taxa_per_rank;
    
    // Diversity metrics
    int unique_species = 0;
    int unique_genera = 0;
    float simpson_diversity = 0.0f;
    float shannon_diversity = 0.0f;
    
    // Confidence statistics
    float avg_confidence = 0.0f;
    float median_confidence = 0.0f;
    
    // Paired-end specific (if applicable)
    int concordant_pairs = 0;
    float avg_pair_concordance = 0.0f;
};

class ClassificationReportGenerator {
private:
    // Data storage
    std::vector<KrakenResult> results;
    std::unordered_map<uint32_t, TaxonInfo> taxonomy;
    ClassificationStats stats;
    
    // Configuration
    std::string taxonomy_file;
    std::string output_directory;
    bool include_unclassified = false;
    float min_abundance_threshold = 0.01f;  // 0.01% minimum
    
    // Rank hierarchy (NCBI standard)
    const std::vector<std::string> rank_hierarchy = {
        "superkingdom", "kingdom", "phylum", "class", "order", 
        "family", "genus", "species", "strain"
    };
    
    // Rank prefixes for MetaPhlAn format
    const std::map<std::string, std::string> rank_prefixes = {
        {"superkingdom", "k__"}, {"kingdom", "k__"}, {"phylum", "p__"},
        {"class", "c__"}, {"order", "o__"}, {"family", "f__"},
        {"genus", "g__"}, {"species", "s__"}, {"strain", "t__"}
    };

public:
    ClassificationReportGenerator(const std::string& output_dir) 
        : output_directory(output_dir) {
        std::filesystem::create_directories(output_directory);
    }
    
    // Main interface methods
    bool load_kraken_results(const std::string& kraken_output_file);
    bool load_taxonomy(const std::string& taxonomy_file);
    
    void set_abundance_threshold(float threshold) { 
        min_abundance_threshold = threshold; 
    }
    void include_unclassified_reads(bool include) { 
        include_unclassified = include; 
    }
    
    // Report generation methods
    bool generate_all_reports(const std::string& sample_name = "sample");
    bool generate_mpa_report(const std::string& output_file);
    bool generate_species_counts(const std::string& output_file);
    bool generate_genus_counts(const std::string& output_file);
    bool generate_summary_stats(const std::string& output_file);
    bool generate_detailed_report(const std::string& output_file);
    
    // Analysis methods
    void compute_statistics();
    void build_taxonomic_lineages();
    ClassificationStats get_statistics() const { return stats; }
    
    // Utility methods
    void print_summary() const;
    std::vector<std::pair<std::string, int>> get_top_species(int n = 20) const;
    std::vector<std::pair<std::string, int>> get_top_genera(int n = 20) const;

private:
    // Internal helper methods
    bool parse_kraken_line(const std::string& line, KrakenResult& result);
    void parse_paired_field(const std::string& field, int& val1, int& val2);
    std::string build_lineage_string(uint32_t taxon_id);
    std::string get_rank_at_level(uint32_t taxon_id, const std::string& target_rank);
    float calculate_shannon_diversity();
    float calculate_simpson_diversity();
    std::string format_abundance(float abundance, bool as_percentage = true);
};

// ================================================================
// IMPLEMENTATION
// ================================================================

bool ClassificationReportGenerator::load_kraken_results(const std::string& kraken_output_file) {
    std::cout << "Loading Kraken results from: " << kraken_output_file << std::endl;
    
    std::ifstream file(kraken_output_file);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open Kraken output file: " << kraken_output_file << std::endl;
        return false;
    }
    
    results.clear();
    std::string line;
    int line_number = 0;
    
    while (std::getline(file, line)) {
        line_number++;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        KrakenResult result;
        if (parse_kraken_line(line, result)) {
            results.push_back(result);
        } else {
            std::cerr << "Warning: Could not parse line " << line_number << ": " << line << std::endl;
        }
    }
    
    file.close();
    
    std::cout << "Loaded " << results.size() << " classification results" << std::endl;
    return !results.empty();
}

bool ClassificationReportGenerator::parse_kraken_line(const std::string& line, KrakenResult& result) {
    std::istringstream iss(line);
    std::string field;
    std::vector<std::string> fields;
    
    // Split by tabs
    while (std::getline(iss, field, '\t')) {
        fields.push_back(field);
    }
    
    if (fields.size() < 5) return false;
    
    // Parse basic fields
    result.status = fields[0][0];
    result.read_id = fields[1];
    result.taxon_id = std::stoul(fields[2]);
    result.length = fields[3];
    result.confidence = std::stof(fields[4]);
    
    if (fields.size() >= 6) result.votes = fields[5];
    if (fields.size() >= 7) result.kmers = fields[6];
    
    // Check if this is paired-end data
    result.is_paired = (result.length.find('|') != std::string::npos);
    
    if (result.is_paired) {
        // Parse paired-end fields
        parse_paired_field(result.length, result.read1_length, result.read2_length);
        parse_paired_field(result.votes, result.read1_votes, result.read2_votes);
        parse_paired_field(result.kmers, result.read1_kmers, result.read2_kmers);
    } else {
        // Single-end data
        result.read1_length = std::stoi(result.length);
        result.read1_votes = result.votes.empty() ? 0 : std::stoi(result.votes);
        result.read1_kmers = result.kmers.empty() ? 0 : std::stoi(result.kmers);
    }
    
    return true;
}

void ClassificationReportGenerator::parse_paired_field(const std::string& field, int& val1, int& val2) {
    size_t pipe_pos = field.find('|');
    if (pipe_pos != std::string::npos) {
        val1 = std::stoi(field.substr(0, pipe_pos));
        val2 = std::stoi(field.substr(pipe_pos + 1));
    } else {
        val1 = std::stoi(field);
        val2 = 0;
    }
}

bool ClassificationReportGenerator::load_taxonomy(const std::string& taxonomy_file_path) {
    taxonomy_file = taxonomy_file_path;
    std::cout << "Loading taxonomy from: " << taxonomy_file << std::endl;
    
    std::ifstream file(taxonomy_file);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open taxonomy file: " << taxonomy_file << std::endl;
        return false;
    }
    
    taxonomy.clear();
    std::string line;
    bool first_line = true;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        // Skip header line
        if (first_line && line.find("taxon_id") != std::string::npos) {
            first_line = false;
            continue;
        }
        first_line = false;
        
        std::istringstream iss(line);
        std::string field;
        std::vector<std::string> fields;
        
        while (std::getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() >= 3) {
            TaxonInfo taxon;
            taxon.taxon_id = std::stoul(fields[0]);
            taxon.name = fields[1];
            taxon.parent_id = std::stoul(fields[2]);
            
            // Try to infer rank from name or use default
            std::string name_lower = taxon.name;
            std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
            
            if (name_lower.find("species") != std::string::npos) taxon.rank = "species";
            else if (name_lower.find("genus") != std::string::npos) taxon.rank = "genus";
            else if (name_lower.find("family") != std::string::npos) taxon.rank = "family";
            else if (name_lower.find("order") != std::string::npos) taxon.rank = "order";
            else if (name_lower.find("class") != std::string::npos) taxon.rank = "class";
            else if (name_lower.find("phylum") != std::string::npos) taxon.rank = "phylum";
            else taxon.rank = "no rank";
            
            taxonomy[taxon.taxon_id] = taxon;
        }
    }
    
    file.close();
    
    std::cout << "Loaded " << taxonomy.size() << " taxonomic entries" << std::endl;
    return !taxonomy.empty();
}

void ClassificationReportGenerator::compute_statistics() {
    std::cout << "Computing classification statistics..." << std::endl;
    
    stats = ClassificationStats();  // Reset
    
    std::map<uint32_t, int> taxon_counts;
    std::vector<float> confidence_scores;
    std::set<uint32_t> unique_species_ids, unique_genus_ids;
    
    // Basic counts
    for (const auto& result : results) {
        stats.total_reads++;
        
        if (result.status == 'C' && result.taxon_id > 0) {
            stats.classified_reads++;
            taxon_counts[result.taxon_id]++;
            confidence_scores.push_back(result.confidence);
            
            // Track species and genera
            std::string species = get_rank_at_level(result.taxon_id, "species");
            std::string genus = get_rank_at_level(result.taxon_id, "genus");
            
            if (!species.empty()) unique_species_ids.insert(result.taxon_id);
            if (!genus.empty()) unique_genus_ids.insert(result.taxon_id);
            
        } else {
            stats.unclassified_reads++;
        }
    }
    
    // Store read counts in taxonomy
    for (const auto& [taxon_id, count] : taxon_counts) {
        if (taxonomy.find(taxon_id) != taxonomy.end()) {
            taxonomy[taxon_id].read_count = count;
            taxonomy[taxon_id].abundance = (float)count / stats.classified_reads * 100.0f;
        }
    }
    
    stats.classification_rate = (float)stats.classified_reads / stats.total_reads * 100.0f;
    stats.unique_species = unique_species_ids.size();
    stats.unique_genera = unique_genus_ids.size();
    
    // Confidence statistics
    if (!confidence_scores.empty()) {
        std::sort(confidence_scores.begin(), confidence_scores.end());
        
        float sum = 0.0f;
        for (float conf : confidence_scores) sum += conf;
        stats.avg_confidence = sum / confidence_scores.size();
        
        size_t median_idx = confidence_scores.size() / 2;
        stats.median_confidence = confidence_scores[median_idx];
    }
    
    // Diversity metrics
    stats.shannon_diversity = calculate_shannon_diversity();
    stats.simpson_diversity = calculate_simpson_diversity();
    
    std::cout << "Statistics computed: " << stats.classified_reads << "/" << stats.total_reads 
              << " reads classified (" << std::fixed << std::setprecision(1) 
              << stats.classification_rate << "%)" << std::endl;
}

float ClassificationReportGenerator::calculate_shannon_diversity() {
    float shannon = 0.0f;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count > 0) {
            float p = (float)taxon.read_count / stats.classified_reads;
            shannon -= p * log2(p);
        }
    }
    
    return shannon;
}

float ClassificationReportGenerator::calculate_simpson_diversity() {
    float simpson = 0.0f;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count > 0) {
            float p = (float)taxon.read_count / stats.classified_reads;
            simpson += p * p;
        }
    }
    
    return 1.0f - simpson;  // Simpson's diversity index
}

std::string ClassificationReportGenerator::get_rank_at_level(uint32_t taxon_id, const std::string& target_rank) {
    uint32_t current_id = taxon_id;
    
    // Walk up the taxonomy tree
    for (int depth = 0; depth < 20; depth++) {  // Prevent infinite loops
        if (taxonomy.find(current_id) == taxonomy.end()) break;
        
        const TaxonInfo& taxon = taxonomy[current_id];
        if (taxon.rank == target_rank) {
            return taxon.name;
        }
        
        if (taxon.parent_id == 0 || taxon.parent_id == current_id) break;
        current_id = taxon.parent_id;
    }
    
    return "";  // Not found
}

bool ClassificationReportGenerator::generate_mpa_report(const std::string& output_file) {
    std::cout << "Generating MetaPhlAn-style report: " << output_file << std::endl;
    
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create MetaPhlAn report file: " << output_file << std::endl;
        return false;
    }
    
    // Header
    out << "#mpa_v30_CHOCOPhlAn_201901\n";
    out << "#clade_name\trelative_abundance\n";
    
    // Build hierarchical profile
    std::map<std::string, float> lineage_abundances;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count == 0) continue;
        if (taxon.abundance < min_abundance_threshold) continue;
        
        // Build full lineage for this taxon
        std::string lineage = build_lineage_string(taxon_id);
        if (!lineage.empty()) {
            lineage_abundances[lineage] = taxon.abundance;
        }
    }
    
    // Sort by abundance (descending)
    std::vector<std::pair<std::string, float>> sorted_lineages;
    for (const auto& [lineage, abundance] : lineage_abundances) {
        sorted_lineages.push_back({lineage, abundance});
    }
    
    std::sort(sorted_lineages.begin(), sorted_lineages.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Write results
    for (const auto& [lineage, abundance] : sorted_lineages) {
        out << lineage << "\t" << std::fixed << std::setprecision(5) << abundance << "\n";
    }
    
    out.close();
    std::cout << "MetaPhlAn report written with " << sorted_lineages.size() << " entries" << std::endl;
    return true;
}

std::string ClassificationReportGenerator::build_lineage_string(uint32_t taxon_id) {
    std::vector<std::string> lineage_parts;
    uint32_t current_id = taxon_id;
    
    // Walk up the taxonomy tree
    for (int depth = 0; depth < 10; depth++) {
        if (taxonomy.find(current_id) == taxonomy.end()) break;
        
        const TaxonInfo& taxon = taxonomy[current_id];
        
        // Add this level if it has a meaningful rank
        if (rank_prefixes.find(taxon.rank) != rank_prefixes.end()) {
            std::string prefix = rank_prefixes.at(taxon.rank);
            lineage_parts.insert(lineage_parts.begin(), prefix + taxon.name);
        }
        
        if (taxon.parent_id == 0 || taxon.parent_id == current_id) break;
        current_id = taxon.parent_id;
    }
    
    if (lineage_parts.empty()) return "";
    
    // Join with |
    std::string result = lineage_parts[0];
    for (size_t i = 1; i < lineage_parts.size(); i++) {
        result += "|" + lineage_parts[i];
    }
    
    return result;
}

bool ClassificationReportGenerator::generate_species_counts(const std::string& output_file) {
    std::cout << "Generating species count report: " << output_file << std::endl;
    
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create species report file: " << output_file << std::endl;
        return false;
    }
    
    // Header
    out << "Species\tTaxon_ID\tRead_Count\tRelative_Abundance\n";
    
    // Collect species-level classifications
    std::map<std::string, std::pair<uint32_t, int>> species_counts;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count == 0) continue;
        
        std::string species = get_rank_at_level(taxon_id, "species");
        if (species.empty()) continue;
        
        if (species_counts.find(species) == species_counts.end()) {
            species_counts[species] = {taxon_id, taxon.read_count};
        } else {
            species_counts[species].second += taxon.read_count;
        }
    }
    
    // Sort by read count
    std::vector<std::tuple<std::string, uint32_t, int, float>> sorted_species;
    for (const auto& [species, data] : species_counts) {
        float abundance = (float)data.second / stats.classified_reads * 100.0f;
        if (abundance >= min_abundance_threshold) {
            sorted_species.push_back({species, data.first, data.second, abundance});
        }
    }
    
    std::sort(sorted_species.begin(), sorted_species.end(),
              [](const auto& a, const auto& b) { return std::get<2>(a) > std::get<2>(b); });
    
    // Write results
    for (const auto& [species, taxon_id, count, abundance] : sorted_species) {
        out << species << "\t" << taxon_id << "\t" << count << "\t" 
            << std::fixed << std::setprecision(3) << abundance << "%\n";
    }
    
    out.close();
    std::cout << "Species report written with " << sorted_species.size() << " species" << std::endl;
    return true;
}

bool ClassificationReportGenerator::generate_genus_counts(const std::string& output_file) {
    std::cout << "Generating genus count report: " << output_file << std::endl;
    
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create genus report file: " << output_file << std::endl;
        return false;
    }
    
    // Header
    out << "Genus\tTaxon_ID\tRead_Count\tRelative_Abundance\n";
    
    // Collect genus-level classifications
    std::map<std::string, std::pair<uint32_t, int>> genus_counts;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count == 0) continue;
        
        std::string genus = get_rank_at_level(taxon_id, "genus");
        if (genus.empty()) continue;
        
        if (genus_counts.find(genus) == genus_counts.end()) {
            genus_counts[genus] = {taxon_id, taxon.read_count};
        } else {
            genus_counts[genus].second += taxon.read_count;
        }
    }
    
    // Sort by read count
    std::vector<std::tuple<std::string, uint32_t, int, float>> sorted_genera;
    for (const auto& [genus, data] : genus_counts) {
        float abundance = (float)data.second / stats.classified_reads * 100.0f;
        if (abundance >= min_abundance_threshold) {
            sorted_genera.push_back({genus, data.first, data.second, abundance});
        }
    }
    
    std::sort(sorted_genera.begin(), sorted_genera.end(),
              [](const auto& a, const auto& b) { return std::get<2>(a) > std::get<2>(b); });
    
    // Write results
    for (const auto& [genus, taxon_id, count, abundance] : sorted_genera) {
        out << genus << "\t" << taxon_id << "\t" << count << "\t" 
            << std::fixed << std::setprecision(3) << abundance << "%\n";
    }
    
    out.close();
    std::cout << "Genus report written with " << sorted_genera.size() << " genera" << std::endl;
    return true;
}

bool ClassificationReportGenerator::generate_summary_stats(const std::string& output_file) {
    std::cout << "Generating summary statistics: " << output_file << std::endl;
    
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create summary file: " << output_file << std::endl;
        return false;
    }
    
    out << "Classification Summary Report\n";
    out << "============================\n\n";
    
    // Basic statistics
    out << "Basic Statistics:\n";
    out << "  Total reads: " << stats.total_reads << "\n";
    out << "  Classified reads: " << stats.classified_reads << " (" 
        << std::fixed << std::setprecision(1) << stats.classification_rate << "%)\n";
    out << "  Unclassified reads: " << stats.unclassified_reads << " (" 
        << std::fixed << std::setprecision(1) << (100.0f - stats.classification_rate) << "%)\n\n";
    
    // Diversity metrics
    out << "Diversity Metrics:\n";
    out << "  Unique species detected: " << stats.unique_species << "\n";
    out << "  Unique genera detected: " << stats.unique_genera << "\n";
    out << "  Shannon diversity index: " << std::fixed << std::setprecision(3) << stats.shannon_diversity << "\n";
    out << "  Simpson diversity index: " << std::fixed << std::setprecision(3) << stats.simpson_diversity << "\n\n";
    
    // Confidence metrics
    out << "Classification Quality:\n";
    out << "  Average confidence: " << std::fixed << std::setprecision(3) << stats.avg_confidence << "\n";
    out << "  Median confidence: " << std::fixed << std::setprecision(3) << stats.median_confidence << "\n\n";
    
    // Top species
    out << "Top 10 Species by Read Count:\n";
    auto top_species = get_top_species(10);
    for (size_t i = 0; i < top_species.size(); i++) {
        float pct = (float)top_species[i].second / stats.classified_reads * 100.0f;
        out << "  " << (i+1) << ". " << top_species[i].first 
            << " (" << top_species[i].second << " reads, " 
            << std::fixed << std::setprecision(2) << pct << "%)\n";
    }
    
    out << "\nTop 10 Genera by Read Count:\n";
    auto top_genera = get_top_genera(10);
    for (size_t i = 0; i < top_genera.size(); i++) {
        float pct = (float)top_genera[i].second / stats.classified_reads * 100.0f;
        out << "  " << (i+1) << ". " << top_genera[i].first 
            << " (" << top_genera[i].second << " reads, " 
            << std::fixed << std::setprecision(2) << pct << "%)\n";
    }
    
    out.close();
    std::cout << "Summary statistics written" << std::endl;
    return true;
}

std::vector<std::pair<std::string, int>> ClassificationReportGenerator::get_top_species(int n) const {
    std::map<std::string, int> species_counts;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count == 0) continue;
        
        std::string species = const_cast<ClassificationReportGenerator*>(this)->get_rank_at_level(taxon_id, "species");
        if (!species.empty()) {
            species_counts[species] += taxon.read_count;
        }
    }
    
    std::vector<std::pair<std::string, int>> sorted_species;
    for (const auto& [species, count] : species_counts) {
        sorted_species.push_back({species, count});
    }
    
    std::sort(sorted_species.begin(), sorted_species.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    if (sorted_species.size() > n) {
        sorted_species.resize(n);
    }
    
    return sorted_species;
}

std::vector<std::pair<std::string, int>> ClassificationReportGenerator::get_top_genera(int n) const {
    std::map<std::string, int> genus_counts;
    
    for (const auto& [taxon_id, taxon] : taxonomy) {
        if (taxon.read_count == 0) continue;
        
        std::string genus = const_cast<ClassificationReportGenerator*>(this)->get_rank_at_level(taxon_id, "genus");
        if (!genus.empty()) {
            genus_counts[genus] += taxon.read_count;
        }
    }
    
    std::vector<std::pair<std::string, int>> sorted_genera;
    for (const auto& [genus, count] : genus_counts) {
        sorted_genera.push_back({genus, count});
    }
    
    std::sort(sorted_genera.begin(), sorted_genera.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    if (sorted_genera.size() > n) {
        sorted_genera.resize(n);
    }
    
    return sorted_genera;
}

bool ClassificationReportGenerator::generate_all_reports(const std::string& sample_name) {
    std::cout << "\nGenerating comprehensive reports for sample: " << sample_name << std::endl;
    
    // Compute statistics first
    compute_statistics();
    
    // Generate all report types
    bool success = true;
    
    success &= generate_mpa_report(output_directory + "/" + sample_name + "_profile.mpa");
    success &= generate_species_counts(output_directory + "/" + sample_name + "_species.tsv");
    success &= generate_genus_counts(output_directory + "/" + sample_name + "_genus.tsv");
    success &= generate_summary_stats(output_directory + "/" + sample_name + "_summary.txt");
    
    if (success) {
        std::cout << "\n✓ All reports generated successfully in: " << output_directory << std::endl;
        print_summary();
    } else {
        std::cerr << "\n❌ Some reports failed to generate" << std::endl;
    }
    
    return success;
}

void ClassificationReportGenerator::print_summary() const {
    std::cout << "\n=== CLASSIFICATION SUMMARY ===" << std::endl;
    std::cout << "Total reads: " << stats.total_reads << std::endl;
    std::cout << "Classified: " << stats.classified_reads 
              << " (" << std::fixed << std::setprecision(1) << stats.classification_rate << "%)" << std::endl;
    std::cout << "Unique species: " << stats.unique_species << std::endl;
    std::cout << "Unique genera: " << stats.unique_genera << std::endl;
    std::cout << "Shannon diversity: " << std::fixed << std::setprecision(3) << stats.shannon_diversity << std::endl;
    std::cout << "Average confidence: " << std::fixed << std::setprecision(3) << stats.avg_confidence << std::endl;
}

#endif // CLASSIFICATION_REPORT_GENERATOR_H