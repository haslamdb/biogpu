// fm_index_builder.cpp
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
#include "wavefront_alignment_types.h" // For the FMIndexDevice struct definition

/**
 * @brief A helper class to build the components of an FM-Index.
 *
 * This class is designed for the offline, CPU-based construction of the FM-Index.
 * It takes the full concatenated sequence of the pangenome graph as input and
 * generates the Burrows-Wheeler Transform (BWT), Suffix Array (SA), C-table,
 * and Occurrence (Occ) table. These components are then written to disk
 * to be loaded by the GPU pipeline.
 */
class FMIndexBuilder {
public:
    /**
     * @brief Constructs the FM-Index from a given text.
     * @param text The concatenated DNA sequence from the pangenome graph.
     * IMPORTANT: The text must end with a special character ('$')
     * that is lexicographically smaller than any other character in the text.
     */
    void build(const std::string& text) {
        if (text.empty() || text.back() != '$') {
            throw std::runtime_error("Input text for FM-Index must be non-empty and end with '$'.");
        }

        std::cout << "Starting FM-Index construction for text of size " << text.size() << "..." << std::endl;
        
        // 1. Construct the Suffix Array.
        // This is the most computationally intensive step.
        // A real implementation would use a sophisticated algorithm like SA-IS.
        // For this stub, we use a simple (but slow) lexicographical sort of suffixes.
        std::cout << "  [1/4] Building Suffix Array (using simple sort)..." << std::endl;
        constructSuffixArray(text);
        std::cout << "      Done." << std::endl;

        // 2. Construct the Burrows-Wheeler Transform (BWT) from the Suffix Array.
        std::cout << "  [2/4] Building BWT..." << std::endl;
        constructBWT(text);
        std::cout << "      Done." << std::endl;

        // 3. Construct the C-table and Occurrence table.
        std::cout << "  [3/4] Building C-table and Occurrence table..." << std::endl;
        constructAuxiliaryTables(text);
        std::cout << "      Done." << std::endl;
        
        std::cout << "  [4/4] FM-Index construction complete." << std::endl;
    }

    /**
     * @brief Writes the constructed index components to binary files.
     * @param output_prefix The prefix for the output files (e.g., "pangenome_fm").
     */
    void writeToDisk(const std::string& output_prefix) const {
        std::cout << "Writing FM-Index to disk with prefix: " << output_prefix << std::endl;
        
        // Write BWT
        std::ofstream bwt_file(output_prefix + ".bwt", std::ios::binary);
        bwt_file.write(bwt_.c_str(), bwt_.size());
        
        // Write Suffix Array
        std::ofstream sa_file(output_prefix + ".sa", std::ios::binary);
        sa_file.write(reinterpret_cast<const char*>(suffix_array_.data()), suffix_array_.size() * sizeof(uint32_t));
        
        // Write C-table
        std::ofstream c_file(output_prefix + ".ctab", std::ios::binary);
        c_file.write(reinterpret_cast<const char*>(C_table_.data()), C_table_.size() * sizeof(uint32_t));
        
        // Write Occurrence table
        std::ofstream occ_file(output_prefix + ".occ", std::ios::binary);
        // (A real implementation would have a more complex structure for the Occ table)
        
        std::cout << "  - " << output_prefix << ".bwt" << std::endl;
        std::cout << "  - " << output_prefix << ".sa" << std::endl;
        std::cout << "  - " << output_prefix << ".ctab" << std::endl;
        std::cout << "  - " << output_prefix << ".occ" << std::endl;
    }

private:
    std::vector<uint32_t> suffix_array_;
    std::string bwt_;
    std::vector<uint32_t> C_table_;
    // The Occurrence table can be very large. A real implementation uses
    // checkpointing to reduce its size. This is a simplified stub.
    std::map<char, std::vector<uint32_t>> occ_table_;

    // Dummy implementation of Suffix Array construction.
    void constructSuffixArray(const std::string& text) {
        const size_t n = text.size();
        suffix_array_.resize(n);
        std::vector<const char*> suffixes(n);
        for(size_t i = 0; i < n; ++i) {
            suffixes[i] = &text[i];
            suffix_array_[i] = i;
        }

        // Slow sort. A real implementation MUST use a linear-time algorithm like SA-IS.
        std::sort(suffix_array_.begin(), suffix_array_.end(),
            [&](uint32_t a, uint32_t b) {
                return strcmp(&text[a], &text[b]) < 0;
            });
    }

    // Constructs the BWT string from the Suffix Array.
    void constructBWT(const std::string& text) {
        const size_t n = text.size();
        bwt_.resize(n);
        for (size_t i = 0; i < n; ++i) {
            if (suffix_array_[i] == 0) {
                bwt_[i] = text.back(); // The character preceding the first suffix is the last char of text
            } else {
                bwt_[i] = text[suffix_array_[i] - 1];
            }
        }
    }

    // Constructs the C-table and Occurrence table.
    void constructAuxiliaryTables(const std::string& text) {
        // Alphabet: $, A, C, G, T, N
        std::map<char, uint32_t> counts;
        for (char c : text) {
            counts[c]++;
        }

        C_table_.resize(256, 0); // Assuming ASCII
        uint32_t cumulative_count = 0;
        char alphabet[] = {'$', 'A', 'C', 'G', 'N', 'T'}; // Lexicographically sorted
        for (char c : alphabet) {
            if (counts.count(c)) {
                C_table_[static_cast<unsigned char>(c)] = cumulative_count;
                cumulative_count += counts[c];
            }
        }
        
        // Build Occurrence table (simplified version)
        occ_table_.clear();
        for (char c : alphabet) {
            occ_table_[c].resize(text.size() + 1, 0);
        }

        for (size_t i = 0; i < text.size(); ++i) {
            for (auto& pair : occ_table_) {
                pair.second[i+1] = pair.second[i];
            }
            occ_table_[bwt_[i]][i+1]++;
        }
    }
};

// Main function to drive the builder.
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <pangenome_sequence_file> <output_prefix>" << std::endl;
        std::cerr << "  pangenome_sequence_file: A file containing the concatenated DNA of all graph nodes." << std::endl;
        std::cerr << "  output_prefix: The base name for the output index files (e.g., 'pangenome_fm')." << std::endl;
        return 1;
    }

    std::string sequence_file = argv[1];
    std::string output_prefix = argv[2];

    // 1. Read the concatenated pangenome sequence from disk.
    std::cout << "Reading pangenome sequence from " << sequence_file << "..." << std::endl;
    std::ifstream ifs(sequence_file);
    if (!ifs) {
        std::cerr << "Error: Cannot open sequence file " << sequence_file << std::endl;
        return 1;
    }
    std::string pangenome_text((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    pangenome_text += '$'; // Append the special terminator character.

    // 2. Build the FM-Index.
    FMIndexBuilder builder;
    try {
        builder.build(pangenome_text);
    } catch (const std::exception& e) {
        std::cerr << "Error during FM-Index construction: " << e.what() << std::endl;
        return 1;
    }

    // 3. Write the index to disk.
    builder.writeToDisk(output_prefix);
    
    std::cout << "\nFM-Index successfully generated." << std::endl;

    return 0;
}
