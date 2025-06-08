// include/biogpu/algorithms/resistance_attribution.h
#ifndef BIOGPU_ALGORITHMS_RESISTANCE_ATTRIBUTION_H
#define BIOGPU_ALGORITHMS_RESISTANCE_ATTRIBUTION_H

#include <string>
#include <vector>
#include <memory>

namespace biogpu {
namespace algorithms {

enum class SequenceType {
    RESISTANCE_GENE,
    ORGANISM,
    PLASMID,
    CONTIG,
    UNKNOWN
};

struct SequenceFragment {
    std::string id;
    SequenceType type;
    float coverage;
    float identity;
    int length;
    int kmerCount;
};

struct AlignmentResult {
    std::string queryId;
    std::string targetId;
    float identity;
    float coverage;
    int alignmentLength;
};

struct OrganismCandidate {
    std::string organismId;
    float confidence;
};

struct ResistanceAttributionResult {
    std::string resistanceGeneId;
    std::vector<OrganismCandidate> candidateOrganisms;
};

struct PathNode {
    int nodeId;
    float weight;
};

struct PathScore {
    std::vector<std::string> path;
    float score;
    std::string targetOrganism;
};

/**
 * GPU-accelerated resistance gene attribution using assembly graphs
 */
class ResistanceAttribution {
public:
    ResistanceAttribution();
    ~ResistanceAttribution();
    
    /**
     * Build assembly graph from sequence fragments and alignments
     */
    void buildAssemblyGraph(
        const std::vector<SequenceFragment>& fragments,
        const std::vector<AlignmentResult>& alignments
    );
    
    /**
     * Attribute resistance genes to organisms
     */
    std::vector<ResistanceAttributionResult> attributeResistance(
        const std::vector<std::string>& resistanceGeneIds
    );
    
    /**
     * Score all paths from a gene to potential organisms
     */
    std::vector<PathScore> scoreAllPaths(
        const std::string& geneId,
        float minScore = 0.5f
    );
    
private:
    class ResistanceAttributionImpl;
    std::unique_ptr<ResistanceAttributionImpl> pImpl;
};

} // namespace algorithms
} // namespace biogpu

#endif // BIOGPU_ALGORITHMS_RESISTANCE_ATTRIBUTION_H