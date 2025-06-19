// kraken_database_loader.h
#ifndef KRAKEN_DATABASE_LOADER_H
#define KRAKEN_DATABASE_LOADER_H

#include <string>
#include <cstdint>
#include <memory>

// Kraken2 database structures
struct KrakenDBHeader {
    uint64_t version;
    uint64_t k_size;
    uint64_t minimizer_len;
    uint64_t taxonomy_nodes;
    uint64_t database_size;
};

// Memory-mapped Kraken2 database interface
class Kraken2Database {
private:
    class Impl;  // Forward declaration for PIMPL
    std::unique_ptr<Impl> pImpl;
    
public:
    Kraken2Database(const std::string& db_path);
    ~Kraken2Database();
    
    // Lookup single minimizer
    uint32_t lookup_minimizer(uint64_t minimizer) const;
    
    // Batch lookup for better cache efficiency
    void lookup_batch(const uint64_t* minimizers, uint32_t* taxon_ids, size_t count) const;
    
    // Find lowest common ancestor
    uint32_t find_lca(uint32_t taxon1, uint32_t taxon2) const;
    
    // Get database statistics
    size_t get_num_entries() const;
    uint32_t get_k_size() const;
    uint32_t get_minimizer_len() const;
};

#endif // KRAKEN_DATABASE_LOADER_H