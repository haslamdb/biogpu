// stdlib/algorithms/resistance_attribution.cpp
#include "biogpu/algorithms/resistance_attribution.h"
#include "biogpu/types/bio_types.h"
#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <cmath>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

namespace biogpu {
namespace algorithms {

#ifdef USE_CUDA
// GPU kernel for parallel path scoring
__global__ void scorePathsKernel(
    const PathNode* paths,
    const float* coverageScores,
    const float* identityScores,
    float* pathScores,
    int numPaths,
    int maxPathLength
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numPaths) return;
    
    float totalScore = 0.0f;
    int pathStart = tid * maxPathLength;
    
    for (int i = 0; i < maxPathLength; ++i) {
        int nodeIdx = paths[pathStart + i].nodeId;
        if (nodeIdx < 0) break;  // End of path
        
        float coverage = coverageScores[nodeIdx];
        float identity = identityScores[nodeIdx];
        float weight = paths[pathStart + i].weight;
        
        // Weighted score combining coverage, identity, and k-mer weight
        totalScore += weight * coverage * identity;
    }
    
    pathScores[tid] = totalScore;
}

// GPU kernel for graph traversal using BFS-like approach
__global__ void findConnectedComponentsKernel(
    const int* adjacencyList,
    const int* adjacencyOffsets,
    int* componentLabels,
    bool* changed,
    int numNodes
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numNodes) return;
    
    int currentLabel = componentLabels[tid];
    int start = adjacencyOffsets[tid];
    int end = adjacencyOffsets[tid + 1];
    
    for (int i = start; i < end; ++i) {
        int neighbor = adjacencyList[i];
        int neighborLabel = componentLabels[neighbor];
        
        if (neighborLabel < currentLabel) {
            atomicMin(&componentLabels[tid], neighborLabel);
            *changed = true;
        }
    }
}
#endif // USE_CUDA

// CPU implementation of path scoring
void scorePathsCPU(
    const std::vector<PathNode>& paths,
    const std::vector<float>& coverageScores,
    const std::vector<float>& identityScores,
    std::vector<float>& pathScores,
    int numPaths,
    int maxPathLength
) {
    #pragma omp parallel for
    for (int tid = 0; tid < numPaths; ++tid) {
        float totalScore = 0.0f;
        int pathStart = tid * maxPathLength;
        
        for (int i = 0; i < maxPathLength; ++i) {
            int nodeIdx = paths[pathStart + i].nodeId;
            if (nodeIdx < 0) break;  // End of path
            
            float coverage = coverageScores[nodeIdx];
            float identity = identityScores[nodeIdx];
            float weight = paths[pathStart + i].weight;
            
            // Weighted score combining coverage, identity, and k-mer weight
            totalScore += weight * coverage * identity;
        }
        
        pathScores[tid] = totalScore;
    }
}

// CPU implementation of connected components
void findConnectedComponentsCPU(
    const std::vector<int>& adjacencyList,
    const std::vector<int>& adjacencyOffsets,
    std::vector<int>& componentLabels,
    int numNodes
) {
    bool changed = true;
    while (changed) {
        changed = false;
        
        #pragma omp parallel for
        for (int tid = 0; tid < numNodes; ++tid) {
            int currentLabel = componentLabels[tid];
            int start = adjacencyOffsets[tid];
            int end = adjacencyOffsets[tid + 1];
            
            for (int i = start; i < end; ++i) {
                int neighbor = adjacencyList[i];
                int neighborLabel = componentLabels[neighbor];
                
                if (neighborLabel < currentLabel) {
                    componentLabels[tid] = neighborLabel;
                    changed = true;
                }
            }
        }
    }
}

class ResistanceAttribution::ResistanceAttributionImpl {
private:
    struct GraphNode {
        std::string sequenceId;
        SequenceType type;
        float coverage;
        float identity;
        std::vector<int> neighbors;
        int componentId;
        float kmerWeight;
    };
    
    std::vector<GraphNode> nodes;
    std::unordered_map<std::string, int> sequenceToNode;
    
#ifdef USE_CUDA
    // GPU memory pointers
    float* d_coverageScores = nullptr;
    float* d_identityScores = nullptr;
    float* d_pathScores = nullptr;
    int* d_adjacencyList = nullptr;
    int* d_adjacencyOffsets = nullptr;
    int* d_componentLabels = nullptr;
    
    void allocateGPUMemory(size_t numNodes, size_t numEdges) {
        cudaMalloc(&d_coverageScores, numNodes * sizeof(float));
        cudaMalloc(&d_identityScores, numNodes * sizeof(float));
        cudaMalloc(&d_adjacencyList, numEdges * sizeof(int));
        cudaMalloc(&d_adjacencyOffsets, (numNodes + 1) * sizeof(int));
        cudaMalloc(&d_componentLabels, numNodes * sizeof(int));
    }
    
    void freeGPUMemory() {
        if (d_coverageScores) cudaFree(d_coverageScores);
        if (d_identityScores) cudaFree(d_identityScores);
        if (d_pathScores) cudaFree(d_pathScores);
        if (d_adjacencyList) cudaFree(d_adjacencyList);
        if (d_adjacencyOffsets) cudaFree(d_adjacencyOffsets);
        if (d_componentLabels) cudaFree(d_componentLabels);
    }
#endif
    
public:
    ResistanceAttributionImpl() = default;
    
    ~ResistanceAttributionImpl() {
#ifdef USE_CUDA
        freeGPUMemory();
#endif
    }
    
    void buildGraph(
        const std::vector<SequenceFragment>& fragments,
        const std::vector<AlignmentResult>& alignments
    ) {
        // Clear existing graph
        nodes.clear();
        sequenceToNode.clear();
        
        // Add nodes for each fragment
        for (const auto& fragment : fragments) {
            GraphNode node;
            node.sequenceId = fragment.id;
            node.type = fragment.type;
            node.coverage = fragment.coverage;
            node.identity = fragment.identity;
            node.componentId = nodes.size();
            node.kmerWeight = static_cast<float>(fragment.kmerCount) / fragment.length;
            
            sequenceToNode[fragment.id] = nodes.size();
            nodes.push_back(node);
        }
        
        // Add edges based on alignments
        for (const auto& alignment : alignments) {
            if (sequenceToNode.count(alignment.queryId) && 
                sequenceToNode.count(alignment.targetId)) {
                int queryIdx = sequenceToNode[alignment.queryId];
                int targetIdx = sequenceToNode[alignment.targetId];
                
                // Add bidirectional edges with quality threshold
                if (alignment.identity >= 0.95f && alignment.coverage >= 0.9f) {
                    nodes[queryIdx].neighbors.push_back(targetIdx);
                    nodes[targetIdx].neighbors.push_back(queryIdx);
                }
            }
        }
        
        // Find connected components
        findConnectedComponents();
    }
    
    void findConnectedComponents() {
        int numNodes = nodes.size();
        if (numNodes == 0) return;
        
        // Build adjacency list representation
        std::vector<int> adjacencyList;
        std::vector<int> adjacencyOffsets(numNodes + 1, 0);
        
        for (int i = 0; i < numNodes; ++i) {
            adjacencyOffsets[i] = adjacencyList.size();
            for (int neighbor : nodes[i].neighbors) {
                adjacencyList.push_back(neighbor);
            }
        }
        adjacencyOffsets[numNodes] = adjacencyList.size();
        
        // Initialize component labels
        std::vector<int> componentLabels(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            componentLabels[i] = i;
        }
        
#ifdef USE_CUDA
        // GPU implementation
        allocateGPUMemory(numNodes, adjacencyList.size());
        
        // Copy data to GPU
        cudaMemcpy(d_adjacencyList, adjacencyList.data(), 
                   adjacencyList.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_adjacencyOffsets, adjacencyOffsets.data(), 
                   adjacencyOffsets.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_componentLabels, componentLabels.data(), 
                   numNodes * sizeof(int), cudaMemcpyHostToDevice);
        
        // Run kernel
        bool* d_changed;
        cudaMalloc(&d_changed, sizeof(bool));
        
        bool h_changed = true;
        while (h_changed) {
            h_changed = false;
            cudaMemcpy(d_changed, &h_changed, sizeof(bool), cudaMemcpyHostToDevice);
            
            int blockSize = 256;
            int gridSize = (numNodes + blockSize - 1) / blockSize;
            
            findConnectedComponentsKernel<<<gridSize, blockSize>>>(
                d_adjacencyList, d_adjacencyOffsets, d_componentLabels, 
                d_changed, numNodes
            );
            
            cudaMemcpy(&h_changed, d_changed, sizeof(bool), cudaMemcpyDeviceToHost);
        }
        
        // Copy results back
        cudaMemcpy(componentLabels.data(), d_componentLabels, 
                   numNodes * sizeof(int), cudaMemcpyDeviceToHost);
        
        cudaFree(d_changed);
#else
        // CPU implementation
        findConnectedComponentsCPU(adjacencyList, adjacencyOffsets, componentLabels, numNodes);
#endif
        
        // Update node component IDs
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].componentId = componentLabels[i];
        }
    }
    
    std::vector<ResistanceAttributionResult> attributeResistance(
        const std::vector<std::string>& resistanceGeneIds
    ) {
        std::vector<ResistanceAttributionResult> results;
        
        for (const auto& geneId : resistanceGeneIds) {
            if (sequenceToNode.count(geneId) == 0) continue;
            
            int geneNodeIdx = sequenceToNode[geneId];
            const GraphNode& geneNode = nodes[geneNodeIdx];
            
            ResistanceAttributionResult attribution;
            attribution.resistanceGeneId = geneId;
            
            // Find all organisms in the same component
            std::unordered_map<std::string, float> organismScores;
            
            for (size_t i = 0; i < nodes.size(); ++i) {
                if (nodes[i].componentId == geneNode.componentId &&
                    nodes[i].type == SequenceType::ORGANISM) {
                    
                    // Calculate confidence based on graph distance and coverage
                    float pathScore = calculatePathScore(geneNodeIdx, i);
                    if (pathScore > 0.5f) {
                        OrganismCandidate candidate;
                        candidate.organismId = nodes[i].sequenceId;
                        candidate.confidence = pathScore;
                        attribution.candidateOrganisms.push_back(candidate);
                    }
                }
            }
            
            // Sort by confidence
            std::sort(attribution.candidateOrganisms.begin(), 
                     attribution.candidateOrganisms.end(),
                     [](const OrganismCandidate& a, const OrganismCandidate& b) {
                         return a.confidence > b.confidence;
                     });
            
            results.push_back(attribution);
        }
        
        return results;
    }
    
    std::vector<PathScore> scoreAllPaths(
        const std::string& geneId,
        float minScore
    ) {
        std::vector<PathScore> results;
        
        if (sequenceToNode.count(geneId) == 0) return results;
        
        int geneNodeIdx = sequenceToNode[geneId];
        
        // Find all paths to organisms using BFS
        for (size_t targetIdx = 0; targetIdx < nodes.size(); ++targetIdx) {
            if (nodes[targetIdx].type != SequenceType::ORGANISM) continue;
            
            std::vector<std::string> path = findPath(geneNodeIdx, targetIdx);
            if (!path.empty()) {
                float score = calculatePathScore(geneNodeIdx, targetIdx);
                if (score >= minScore) {
                    PathScore ps;
                    ps.path = path;
                    ps.score = score;
                    ps.targetOrganism = nodes[targetIdx].sequenceId;
                    results.push_back(ps);
                }
            }
        }
        
        // Sort by score
        std::sort(results.begin(), results.end(),
                 [](const PathScore& a, const PathScore& b) {
                     return a.score > b.score;
                 });
        
        return results;
    }
    
private:
    float calculatePathScore(int startIdx, int endIdx) {
        // Simple BFS to find shortest path
        if (nodes[startIdx].componentId != nodes[endIdx].componentId) {
            return 0.0f;
        }
        
        std::queue<std::pair<int, float>> q;
        std::unordered_map<int, float> visited;
        
        q.push({startIdx, 1.0f});
        visited[startIdx] = 1.0f;
        
        while (!q.empty()) {
            auto [currentIdx, currentScore] = q.front();
            q.pop();
            
            if (currentIdx == endIdx) {
                return currentScore;
            }
            
            for (int neighborIdx : nodes[currentIdx].neighbors) {
                float edgeScore = nodes[neighborIdx].coverage * 
                                 nodes[neighborIdx].identity * 
                                 nodes[neighborIdx].kmerWeight;
                float newScore = currentScore * edgeScore * 0.9f; // Distance penalty
                
                if (visited.count(neighborIdx) == 0 || visited[neighborIdx] < newScore) {
                    visited[neighborIdx] = newScore;
                    q.push({neighborIdx, newScore});
                }
            }
        }
        
        return 0.0f;
    }
    
    std::vector<std::string> findPath(int startIdx, int endIdx) {
        std::vector<std::string> path;
        
        if (nodes[startIdx].componentId != nodes[endIdx].componentId) {
            return path;
        }
        
        // BFS to find path
        std::queue<int> q;
        std::unordered_map<int, int> parent;
        
        q.push(startIdx);
        parent[startIdx] = -1;
        
        while (!q.empty()) {
            int current = q.front();
            q.pop();
            
            if (current == endIdx) {
                // Reconstruct path
                int idx = endIdx;
                while (idx != -1) {
                    path.push_back(nodes[idx].sequenceId);
                    idx = parent[idx];
                }
                std::reverse(path.begin(), path.end());
                return path;
            }
            
            for (int neighbor : nodes[current].neighbors) {
                if (parent.count(neighbor) == 0) {
                    parent[neighbor] = current;
                    q.push(neighbor);
                }
            }
        }
        
        return path;
    }
};

// Public API implementation
ResistanceAttribution::ResistanceAttribution() 
    : pImpl(std::make_unique<ResistanceAttributionImpl>()) {}

ResistanceAttribution::~ResistanceAttribution() = default;

void ResistanceAttribution::buildAssemblyGraph(
    const std::vector<SequenceFragment>& fragments,
    const std::vector<AlignmentResult>& alignments
) {
    pImpl->buildGraph(fragments, alignments);
}

std::vector<ResistanceAttributionResult> ResistanceAttribution::attributeResistance(
    const std::vector<std::string>& resistanceGeneIds
) {
    return pImpl->attributeResistance(resistanceGeneIds);
}

std::vector<PathScore> ResistanceAttribution::scoreAllPaths(
    const std::string& geneId,
    float minScore
) {
    return pImpl->scoreAllPaths(geneId, minScore);
}

} // namespace algorithms
} // namespace biogpu