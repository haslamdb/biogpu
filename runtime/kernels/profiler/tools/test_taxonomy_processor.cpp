// test_taxonomy_processor.cpp
// Test program for the taxonomy processor module

#include "../enhanced_k2like/taxonomy/taxonomy_processor.h"
#include "../enhanced_k2like/gpu_kraken_types.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <cassert>

// Create test taxonomy data
void create_test_taxonomy_file(const std::string& filename) {
    std::ofstream out(filename);
    out << "taxon_id\tparent_id\tname\trank\n";
    out << "1\t0\troot\troot\n";
    out << "2\t1\tBacteria\tsuperkingdom\n";
    out << "3\t1\tArchaea\tsuperkingdom\n";
    out << "131567\t1\tcellular organisms\tno rank\n";
    out << "1224\t2\tProteobacteria\tphylum\n";
    out << "1236\t1224\tGammaproteobacteria\tclass\n";
    out << "91347\t1236\tEnterobacterales\torder\n";
    out << "543\t91347\tEnterobacteriaceae\tfamily\n";
    out << "561\t543\tEscherichia\tgenus\n";
    out << "562\t561\tEscherichia coli\tspecies\n";
    out << "511145\t562\tEscherichia coli str. K-12 substr. MG1655\tno rank\n";
    out << "1239\t2\tFirmicutes\tphylum\n";
    out << "1385\t1239\tBacilli\tclass\n";
    out << "1386\t1385\tBacillales\torder\n";
    out << "90964\t1386\tStaphylococcaceae\tfamily\n";
    out << "1279\t90964\tStaphylococcus\tgenus\n";
    out << "1280\t1279\tStaphylococcus aureus\tspecies\n";
    out.close();
}

// Test SimpleTaxonomyProcessor
void test_simple_taxonomy() {
    std::cout << "\n=== Testing SimpleTaxonomyProcessor ===" << std::endl;
    
    // Create test file
    std::string test_file = "test_taxonomy.tsv";
    create_test_taxonomy_file(test_file);
    
    SimpleTaxonomyProcessor simple_tax;
    
    // Load taxonomy
    if (!simple_tax.load_taxonomy_tsv(test_file)) {
        std::cerr << "Failed to load test taxonomy" << std::endl;
        return;
    }
    
    std::cout << "Loaded " << simple_tax.size() << " taxa" << std::endl;
    
    // Test basic lookups
    std::cout << "\nTesting basic lookups:" << std::endl;
    uint32_t ecoli = 562;
    std::cout << "E. coli (" << ecoli << "):" << std::endl;
    std::cout << "  Name: " << simple_tax.get_name(ecoli) << std::endl;
    std::cout << "  Rank: " << simple_tax.get_rank(ecoli) << std::endl;
    std::cout << "  Parent: " << simple_tax.get_parent(ecoli) << std::endl;
    
    // Test LCA computation
    std::cout << "\nTesting LCA computation:" << std::endl;
    uint32_t ecoli_k12 = 511145;
    uint32_t staph = 1280;
    
    uint32_t lca1 = simple_tax.compute_simple_lca(ecoli, ecoli_k12);
    std::cout << "LCA(E. coli, E. coli K-12) = " << lca1 
              << " (" << simple_tax.get_name(lca1) << ")" << std::endl;
    assert(lca1 == ecoli); // K-12 is a strain of E. coli
    
    uint32_t lca2 = simple_tax.compute_simple_lca(ecoli, staph);
    std::cout << "LCA(E. coli, S. aureus) = " << lca2 
              << " (" << simple_tax.get_name(lca2) << ")" << std::endl;
    assert(lca2 == 2); // Both are bacteria
    
    // Test multiple LCA
    std::vector<uint32_t> taxa = {ecoli, ecoli_k12, staph};
    uint32_t lca_multi = simple_tax.compute_lca_of_list(taxa);
    std::cout << "LCA of multiple taxa = " << lca_multi 
              << " (" << simple_tax.get_name(lca_multi) << ")" << std::endl;
    
    // Cleanup
    std::remove(test_file.c_str());
    std::cout << "✓ SimpleTaxonomyProcessor tests passed" << std::endl;
}

// Test PhylogeneticUtils functions
void test_phylogenetic_utils() {
    std::cout << "\n=== Testing PhylogeneticUtils ===" << std::endl;
    
    // Test species validation
    std::vector<uint32_t> valid_species = {562, 1280, 1313};
    std::vector<uint32_t> invalid_species = {0, 999999999};
    
    assert(PhylogeneticUtils::validate_species_list(valid_species) == true);
    assert(PhylogeneticUtils::validate_species_list(invalid_species) == false);
    assert(PhylogeneticUtils::validate_species_list({}) == false);
    
    std::cout << "✓ Species validation tests passed" << std::endl;
    
    // Test taxonomy line parsing
    uint32_t taxon_id, parent_id;
    std::string name, rank;
    std::string test_line = "562\t561\tEscherichia coli\tspecies";
    
    bool parsed = PhylogeneticUtils::parse_taxonomy_line(test_line, taxon_id, parent_id, name, rank);
    assert(parsed == true);
    assert(taxon_id == 562);
    assert(parent_id == 561);
    assert(name == "Escherichia coli");
    assert(rank == "species");
    
    std::cout << "✓ Taxonomy line parsing tests passed" << std::endl;
}

// Test LCAAlgorithms
void test_lca_algorithms() {
    std::cout << "\n=== Testing LCAAlgorithms ===" << std::endl;
    
    // Create test parent map
    std::unordered_map<uint32_t, uint32_t> parents = {
        {1, 0},      // root
        {2, 1},      // bacteria -> root
        {3, 1},      // archaea -> root
        {10, 2},     // proteobacteria -> bacteria
        {11, 2},     // firmicutes -> bacteria
        {20, 10},    // e. coli -> proteobacteria
        {21, 10},    // salmonella -> proteobacteria
        {22, 11},    // bacillus -> firmicutes
        {23, 11},    // clostridium -> firmicutes
        {30, 3}      // methanogen -> archaea
    };
    
    // Test pair LCA
    uint32_t lca1 = LCAAlgorithms::compute_lca_pair(20, 21, parents);
    assert(lca1 == 10); // Both are proteobacteria
    std::cout << "LCA(20, 21) = " << lca1 << " ✓" << std::endl;
    
    uint32_t lca2 = LCAAlgorithms::compute_lca_pair(20, 22, parents);
    assert(lca2 == 2); // Common ancestor is bacteria
    std::cout << "LCA(20, 22) = " << lca2 << " ✓" << std::endl;
    
    uint32_t lca3 = LCAAlgorithms::compute_lca_pair(20, 30, parents);
    assert(lca3 == 1); // Common ancestor is root
    std::cout << "LCA(20, 30) = " << lca3 << " ✓" << std::endl;
    
    // Test multiple LCA
    std::vector<uint32_t> taxa = {20, 21, 22};
    uint32_t lca_multi = LCAAlgorithms::compute_lca_multiple(taxa, parents);
    assert(lca_multi == 2); // All are bacteria
    std::cout << "LCA({20, 21, 22}) = " << lca_multi << " ✓" << std::endl;
    
    // Test path to root
    std::vector<uint32_t> path = LCAAlgorithms::get_path_to_root(20, parents);
    assert(path.size() == 4); // 20 -> 10 -> 2 -> 1
    assert(path[0] == 20);
    assert(path[1] == 10);
    assert(path[2] == 2);
    assert(path[3] == 1);
    std::cout << "Path from 20 to root: ";
    for (uint32_t t : path) std::cout << t << " ";
    std::cout << "✓" << std::endl;
    
    std::cout << "✓ All LCA algorithm tests passed" << std::endl;
}

// Test SpeciesTrackingData
void test_species_tracking() {
    std::cout << "\n=== Testing SpeciesTrackingData ===" << std::endl;
    
    SpeciesTrackingData tracker;
    
    // Add some species occurrences
    tracker.add_species_occurrence(562, 5);   // E. coli - 5 genomes
    tracker.add_species_occurrence(1280, 3);  // S. aureus - 3 genomes
    tracker.add_species_occurrence(562, 2);   // More E. coli
    
    assert(tracker.total_species() == 2);
    assert(tracker.total_genomes() == 10);
    assert(tracker.get_genome_count(562) == 7);
    assert(tracker.get_genome_count(1280) == 3);
    
    std::cout << "Species tracking:" << std::endl;
    std::cout << "  Total species: " << tracker.total_species() << std::endl;
    std::cout << "  Total genomes: " << tracker.total_genomes() << std::endl;
    std::cout << "  E. coli genomes: " << tracker.get_genome_count(562) << std::endl;
    
    // Test species list
    auto species_list = tracker.get_species_list();
    assert(species_list.size() == 2);
    
    std::cout << "✓ Species tracking tests passed" << std::endl;
}

// Test PhylogeneticLCACandidate
void test_phylogenetic_candidate() {
    std::cout << "\n=== Testing PhylogeneticLCACandidate ===" << std::endl;
    
    PhylogeneticLCACandidate candidate;
    candidate.minimizer_hash = 0x123456789ABCDEF0;
    candidate.lca_taxon = 10;
    candidate.contributing_species = {20, 21, 22};
    candidate.genome_counts_per_species = {5, 3, 2};
    candidate.genome_count = 10;
    candidate.phylogenetic_spread = 5;
    candidate.max_phylogenetic_distance = 3;
    candidate.uniqueness_score = 0.8f;
    
    // Test serialization size
    size_t size = candidate.serialized_size();
    std::cout << "Serialized size: " << size << " bytes" << std::endl;
    assert(size > 0);
    
    // Test serialization and deserialization
    std::vector<uint8_t> buffer(size);
    candidate.serialize(buffer.data());
    
    PhylogeneticLCACandidate candidate2;
    size_t bytes_read = candidate2.deserialize(buffer.data());
    assert(bytes_read == size);
    
    // Verify deserialized data
    assert(candidate2.minimizer_hash == candidate.minimizer_hash);
    assert(candidate2.lca_taxon == candidate.lca_taxon);
    assert(candidate2.contributing_species.size() == 3);
    assert(candidate2.genome_counts_per_species.size() == 3);
    assert(candidate2.phylogenetic_spread == candidate.phylogenetic_spread);
    
    std::cout << "✓ Serialization/deserialization tests passed" << std::endl;
}

// Performance test
void test_performance() {
    std::cout << "\n=== Performance Testing ===" << std::endl;
    
    // Create large parent map
    std::unordered_map<uint32_t, uint32_t> parents;
    parents[1] = 0; // root
    
    // Create a deep tree
    const int num_taxa = 100000;
    const int branch_factor = 10;
    
    std::cout << "Building test taxonomy with " << num_taxa << " taxa..." << std::endl;
    
    uint32_t current_id = 2;
    std::vector<uint32_t> current_level = {1};
    std::vector<uint32_t> next_level;
    
    while (current_id < num_taxa) {
        next_level.clear();
        for (uint32_t parent : current_level) {
            for (int i = 0; i < branch_factor && current_id < num_taxa; i++) {
                parents[current_id] = parent;
                next_level.push_back(current_id);
                current_id++;
            }
        }
        current_level = next_level;
    }
    
    std::cout << "Created taxonomy with " << parents.size() << " nodes" << std::endl;
    
    // Benchmark LCA computations
    const int num_queries = 10000;
    std::vector<std::pair<uint32_t, uint32_t>> test_pairs;
    
    // Generate random pairs
    std::srand(42);
    for (int i = 0; i < num_queries; i++) {
        uint32_t t1 = 2 + (std::rand() % (num_taxa - 2));
        uint32_t t2 = 2 + (std::rand() % (num_taxa - 2));
        test_pairs.push_back({t1, t2});
    }
    
    std::cout << "\nBenchmarking " << num_queries << " LCA computations..." << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (const auto& [t1, t2] : test_pairs) {
        LCAAlgorithms::compute_lca_pair(t1, t2, parents);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    double per_query = (double)duration.count() / num_queries;
    double throughput = (double)num_queries / duration.count() * 1000;
    
    std::cout << "Total time: " << duration.count() << " ms" << std::endl;
    std::cout << "Per query: " << per_query << " ms" << std::endl;
    std::cout << "Throughput: " << throughput << " LCAs/second" << std::endl;
}

int main() {
    std::cout << "=== Taxonomy Processor Test Suite ===" << std::endl;
    
    // Run all tests
    test_simple_taxonomy();
    test_phylogenetic_utils();
    test_lca_algorithms();
    test_species_tracking();
    test_phylogenetic_candidate();
    test_performance();
    
    std::cout << "\n=== All tests completed successfully ===" << std::endl;
    
    return 0;
}// test_taxonomy_processor.cpp
// Test program for the taxonomy processor module

#include "../enhanced_k2like/taxonomy/taxonomy_processor.h"
#include "../enhanced_k2like/gpu_kraken_types.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <cassert>

// Create test taxonomy data
void create_test_taxonomy_file(const std::string& filename) {
    std::ofstream out(filename);
    out << "taxon_id\tparent_id\tname\trank\n";
    out << "1\t0\troot\troot\n";
    out << "2\t1\tBacteria\tsuperkingdom\n";
    out << "3\t1\tArchaea\tsuperkingdom\n";
    out << "131567\t1\tcellular organisms\tno rank\n";
    out << "1224\t2\tProteobacteria\tphylum\n";
    out << "1236\t1224\tGammaproteobacteria\tclass\n";
    out << "91347\t1236\tEnterobacterales\torder\n";
    out << "543\t91347\tEnterobacteriaceae\tfamily\n";
    out << "561\t543\tEscherichia\tgenus\n";
    out << "562\t561\tEscherichia coli\tspecies\n";
    out << "511145\t562\tEscherichia coli str. K-12 substr. MG1655\tno rank\n";
    out << "1239\t2\tFirmicutes\tphylum\n";
    out << "1385\t1239\tBacilli\tclass\n";
    out << "1386\t1385\tBacillales\torder\n";
    out << "90964\t1386\tStaphylococcaceae\tfamily\n";
    out << "1279\t90964\tStaphylococcus\tgenus\n";
    out << "1280\t1279\tStaphylococcus aureus\tspecies\n";
    out.close();
}

// Test SimpleTaxonomyProcessor
void test_simple_taxonomy() {
    std::cout << "\n=== Testing SimpleTaxonomyProcessor ===" << std::endl;
    
    // Create test file
    std::string test_file = "test_taxonomy.tsv";
    create_test_taxonomy_file(test_file);
    
    SimpleTaxonomyProcessor simple_tax;
    
    // Load taxonomy
    if (!simple_tax.load_taxonomy_tsv(test_file)) {
        std::cerr << "Failed to load test taxonomy" << std::endl;
        return;
    }
    
    std::cout << "Loaded " << simple_tax.size() << " taxa" << std::endl;
    
    // Test basic lookups
    std::cout << "\nTesting basic lookups:" << std::endl;
    uint32_t ecoli = 562;
    std::cout << "E. coli (" << ecoli << "):" << std::endl;
    std::cout << "  Name: " << simple_tax.get_name(ecoli) << std::endl;
    std::cout << "  Rank: " << simple_tax.get_rank(ecoli) << std::endl;
    std::cout << "  Parent: " << simple_tax.get_parent(ecoli) << std::endl;
    
    // Test LCA computation
    std::cout << "\nTesting LCA computation:" << std::endl;
    uint32_t ecoli_k12 = 511145;
    uint32_t staph = 1280;
    
    uint32_t lca1 = simple_tax.compute_simple_lca(ecoli, ecoli_k12);
    std::cout << "LCA(E. coli, E. coli K-12) = " << lca1 
              << " (" << simple_tax.get_name(lca1) << ")" << std::endl;
    assert(lca1 == ecoli); // K-12 is a strain of E. coli
    
    uint32_t lca2 = simple_tax.compute_simple_lca(ecoli, staph);
    std::cout << "LCA(E. coli, S. aureus) = " << lca2 
              << " (" << simple_tax.get_name(lca2) << ")" << std::endl;
    assert(lca2 == 2); // Both are bacteria
    
    // Test multiple LCA
    std::vector<uint32_t> taxa = {ecoli, ecoli_k12, staph};
    uint32_t lca_multi = simple_tax.compute_lca_of_list(taxa);
    std::cout << "LCA of multiple taxa = " << lca_multi 
              << " (" << simple_tax.get_name(lca_multi) << ")" << std::endl;
    
    // Cleanup
    std::remove(test_file.c_str());
    std::cout << "✓ SimpleTaxonomyProcessor tests passed" << std::endl;
}

// Test PhylogeneticUtils functions
void test_phylogenetic_utils() {
    std::cout << "\n=== Testing PhylogeneticUtils ===" << std::endl;
    
    // Test species validation
    std::vector<uint32_t> valid_species = {562, 1280, 1313};
    std::vector<uint32_t> invalid_species = {0, 999999999};
    
    assert(PhylogeneticUtils::validate_species_list(valid_species) == true);
    assert(PhylogeneticUtils::validate_species_list(invalid_species) == false);
    assert(PhylogeneticUtils::validate_species_list({}) == false);
    
    std::cout << "✓ Species validation tests passed" << std::endl;
    
    // Test taxonomy line parsing
    uint32_t taxon_id, parent_id;
    std::string name, rank;
    std::string test_line = "562\t561\tEscherichia coli\tspecies";
    
    bool parsed = PhylogeneticUtils::parse_taxonomy_line(test_line, taxon_id, parent_id, name, rank);
    assert(parsed == true);
    assert(taxon_id == 562);
    assert(parent_id == 561);
    assert(name == "Escherichia coli");
    assert(rank == "species");
    
    std::cout << "✓ Taxonomy line parsing tests passed" << std::endl;
}

// Test LCAAlgorithms
void test_lca_algorithms() {
    std::cout << "\n=== Testing LCAAlgorithms ===" << std::endl;
    
    // Create test parent map
    std::unordered_map<uint32_t, uint32_t> parents = {
        {1, 0},      // root
        {2, 1},      // bacteria -> root
        {3, 1},      // archaea -> root
        {10, 2},     // proteobacteria -> bacteria
        {11, 2},     // firmicutes -> bacteria
        {20, 10},    // e. coli -> proteobacteria
        {21, 10},    // salmonella -> proteobacteria
        {22, 11},    // bacillus -> firmicutes
        {23, 11},    // clostridium -> firmicutes
        {30, 3}      // methanogen -> archaea
    };
    
    // Test pair LCA
    uint32_t lca1 = LCAAlgorithms::compute_lca_pair(20, 21, parents);
    assert(lca1 == 10); // Both are proteobacteria
    std::cout << "LCA(20, 21) = " << lca1 << " ✓" << std::endl;
    
    uint32_t lca2 = LCAAlgorithms::compute_lca_pair(20, 22, parents);
    assert(lca2 == 2); // Common ancestor is bacteria
    std::cout << "LCA(20, 22) = " << lca2 << " ✓" << std::endl;
    
    uint32_t lca3 = LCAAlgorithms::compute_lca_pair(20, 30, parents);
    assert(lca3 == 1); // Common ancestor is root
    std::cout << "LCA(20, 30) = " << lca3 << " ✓" << std::endl;
    
    // Test multiple LCA
    std::vector<uint32_t> taxa = {20, 21, 22};
    uint32_t lca_multi = LCAAlgorithms::compute_lca_multiple(taxa, parents);
    assert(lca_multi == 2); // All are bacteria
    std::cout << "LCA({20, 21, 22}) = " << lca_multi << " ✓" << std::endl;
    
    // Test path to root
    std::vector<uint32_t> path = LCAAlgorithms::get_path_to_root(20, parents);
    assert(path.size() == 4); // 20 -> 10 -> 2 -> 1
    assert(path[0] == 20);
    assert(path[1] == 10);
    assert(path[2] == 2);
    assert(path[3] == 1);
    std::cout << "Path from 20 to root: ";
    for (uint32_t t : path) std::cout << t << " ";
    std::cout << "✓" << std::endl;
    
    std::cout << "✓ All LCA algorithm tests passed" << std::endl;
}

// Test SpeciesTrackingData
void test_species_tracking() {
    std::cout << "\n=== Testing SpeciesTrackingData ===" << std::endl;
    
    SpeciesTrackingData tracker;
    
    // Add some species occurrences
    tracker.add_species_occurrence(562, 5);   // E. coli - 5 genomes
    tracker.add_species_occurrence(1280, 3);  // S. aureus - 3 genomes
    tracker.add_species_occurrence(562, 2);   // More E. coli
    
    assert(tracker.total_species() == 2);
    assert(tracker.total_genomes() == 10);
    assert(tracker.get_genome_count(562) == 7);
    assert(tracker.get_genome_count(1280) == 3);
    
    std::cout << "Species tracking:" << std::endl;
    std::cout << "  Total species: " << tracker.total_species() << std::endl;
    std::cout << "  Total genomes: " << tracker.total_genomes() << std::endl;
    std::cout << "  E. coli genomes: " << tracker.get_genome_count(562) << std::endl;
    
    // Test species list
    auto species_list = tracker.get_species_list();
    assert(species_list.size() == 2);
    
    std::cout << "✓ Species tracking tests passed" << std::endl;
}

// Test PhylogeneticLCACandidate
void test_phylogenetic_candidate() {
    std::cout << "\n=== Testing PhylogeneticLCACandidate ===" << std::endl;
    
    PhylogeneticLCACandidate candidate;
    candidate.minimizer_hash = 0x123456789ABCDEF0;
    candidate.lca_taxon = 10;
    candidate.contributing_species = {20, 21, 22};
    candidate.genome_counts_per_species = {5, 3, 2};
    candidate.genome_count = 10;
    candidate.phylogenetic_spread = 5;
    candidate.max_phylogenetic_distance = 3;
    candidate.uniqueness_score = 0.8f;
    
    // Test serialization size
    size_t size = candidate.serialized_size();
    std::cout << "Serialized size: " << size << " bytes" << std::endl;
    assert(size > 0);
    
    // Test serialization and deserialization
    std::vector<uint8_t> buffer(size);
    candidate.serialize(buffer.data());
    
    PhylogeneticLCACandidate candidate2;
    size_t bytes_read = candidate2.deserialize(buffer.data());
    assert(bytes_read == size);
    
    // Verify deserialized data
    assert(candidate2.minimizer_hash == candidate.minimizer_hash);
    assert(candidate2.lca_taxon == candidate.lca_taxon);
    assert(candidate2.contributing_species.size() == 3);
    assert(candidate2.genome_counts_per_species.size() == 3);
    assert(candidate2.phylogenetic_spread == candidate.phylogenetic_spread);
    
    std::cout << "✓ Serialization/deserialization tests passed" << std::endl;
}

// Performance test
void test_performance() {
    std::cout << "\n=== Performance Testing ===" << std::endl;
    
    // Create large parent map
    std::unordered_map<uint32_t, uint32_t> parents;
    parents[1] = 0; // root
    
    // Create a deep tree
    const int num_taxa = 100000;
    const int branch_factor = 10;
    
    std::cout << "Building test taxonomy with " << num_taxa << " taxa..." << std::endl;
    
    uint32_t current_id = 2;
    std::vector<uint32_t> current_level = {1};
    std::vector<uint32_t> next_level;
    
    while (current_id < num_taxa) {
        next_level.clear();
        for (uint32_t parent : current_level) {
            for (int i = 0; i < branch_factor && current_id < num_taxa; i++) {
                parents[current_id] = parent;
                next_level.push_back(current_id);
                current_id++;
            }
        }
        current_level = next_level;
    }
    
    std::cout << "Created taxonomy with " << parents.size() << " nodes" << std::endl;
    
    // Benchmark LCA computations
    const int num_queries = 10000;
    std::vector<std::pair<uint32_t, uint32_t>> test_pairs;
    
    // Generate random pairs
    std::srand(42);
    for (int i = 0; i < num_queries; i++) {
        uint32_t t1 = 2 + (std::rand() % (num_taxa - 2));
        uint32_t t2 = 2 + (std::rand() % (num_taxa - 2));
        test_pairs.push_back({t1, t2});
    }
    
    std::cout << "\nBenchmarking " << num_queries << " LCA computations..." << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (const auto& [t1, t2] : test_pairs) {
        LCAAlgorithms::compute_lca_pair(t1, t2, parents);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    double per_query = (double)duration.count() / num_queries;
    double throughput = (double)num_queries / duration.count() * 1000;
    
    std::cout << "Total time: " << duration.count() << " ms" << std::endl;
    std::cout << "Per query: " << per_query << " ms" << std::endl;
    std::cout << "Throughput: " << throughput << " LCAs/second" << std::endl;
}

int main() {
    std::cout << "=== Taxonomy Processor Test Suite ===" << std::endl;
    
    // Run all tests
    test_simple_taxonomy();
    test_phylogenetic_utils();
    test_lca_algorithms();
    test_species_tracking();
    test_phylogenetic_candidate();
    test_performance();
    
    std::cout << "\n=== All tests completed successfully ===" << std::endl;
    
    return 0;
}
