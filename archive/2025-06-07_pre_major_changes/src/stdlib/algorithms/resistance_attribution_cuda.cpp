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

class ResistanceAttributionImpl {
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
    
public:
    ResistanceAttributionImpl() = default;
    
    ~ResistanceAttributionImpl() {
        freeGPUMemory();
    }
    
    void buildGraph(const std::vector<SequenceFragment>& fragments,
                   const std::vector<AlignmentResult>& alignments) {
        nodes.clear();
        sequenceToNode.clear();
        
        // Create nodes for each unique sequence fragment
        for (const auto& fragment : fragments) {
            int nodeId = nodes.size();
            GraphNode node;
            node.sequenceId = fragment.id;
            node.type = fragment.type;
            node.coverage = fragment.coverage;
            node.identity = fragment.identity;
            node.componentId = nodeId;  // Initially each node is its own component
            node.kmerWeight = fragment.kmerCount / static_cast<float>(fragment.length);
            
            nodes.push_back(node);
            sequenceToNode[fragment.id] = nodeId;
        }
        
        // Add edges based on alignments and overlap information
        for (const auto& alignment : alignments) {
            auto it1 = sequenceToNode.find(alignment.queryId);
            auto it2 = sequenceToNode.find(alignment.targetId);
            
            if (it1 != sequenceToNode.end() && it2 != sequenceToNode.end()) {
                int node1 = it1->second;
                int node2 = it2->second;
                
                // Add bidirectional edge if alignment quality is sufficient
                if (alignment.identity >= 0.9 && alignment.coverage >= 0.5) {
                    nodes[node1].neighbors.push_back(node2);
                    nodes[node2].neighbors.push_back(node1);
                }
            }
        }
    }
    
    std::vector<int> findConnectedComponents() {
        int numNodes = nodes.size();
        if (numNodes == 0) return {};
        
        // Prepare adjacency list for GPU
        std::vector<int> adjacencyList;
        std::vector<int> adjacencyOffsets;
        adjacencyOffsets.push_back(0);
        
        for (const auto& node : nodes) {
            adjacencyList.insert(adjacencyList.end(), 
                               node.neighbors.begin(), 
                               node.neighbors.end());
            adjacencyOffsets.push_back(adjacencyList.size());
        }
        
        // Allocate GPU memory
        allocateGPUMemory(numNodes, adjacencyList.size());
        
        // Copy data to GPU
        cudaMemcpy(d_adjacencyList, adjacencyList.data(), 
                  adjacencyList.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_adjacencyOffsets, adjacencyOffsets.data(), 
                  adjacencyOffsets.size() * sizeof(int), cudaMemcpyHostToDevice);
        
        // Initialize component labels
        std::vector<int> componentLabels(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            componentLabels[i] = i;
        }
        cudaMemcpy(d_componentLabels, componentLabels.data(), 
                  numNodes * sizeof(int), cudaMemcpyHostToDevice);
        
        // GPU-accelerated connected components using label propagation
        bool h_changed = true;
        bool* d_changed;
        cudaMalloc(&d_changed, sizeof(bool));
        
        int blockSize = 256;
        int gridSize = (numNodes + blockSize - 1) / blockSize;
        
        while (h_changed) {
            h_changed = false;
            cudaMemcpy(d_changed, &h_changed, sizeof(bool), cudaMemcpyHostToDevice);
            
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
        
        // Update node components
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].componentId = componentLabels[i];
        }
        
        return componentLabels;
    }
    
    std::vector<ResistanceAttribution> attributeResistanceGenes(
        const std::vector<std::string>& resistanceGeneIds) {
        
        std::vector<ResistanceAttribution> attributions;
        
        // First, find connected components
        findConnectedComponents();
        
        // Group nodes by component
        std::unordered_map<int, std::vector<int>> componentNodes;
        for (int i = 0; i < nodes.size(); ++i) {
            componentNodes[nodes[i].componentId].push_back(i);
        }
        
        // For each resistance gene, find the best organism attribution
        for (const auto& geneId : resistanceGeneIds) {
            auto it = sequenceToNode.find(geneId);
            if (it == sequenceToNode.end()) continue;
            
            int geneNode = it->second;
            int componentId = nodes[geneNode].componentId;
            
            ResistanceAttribution attr;
            attr.resistanceGeneId = geneId;
            
            // Find all organisms in the same component
            const auto& component = componentNodes[componentId];
            
            // Score each organism based on path quality
            std::vector<std::pair<std::string, float>> organismScores;
            
            for (int nodeId : component) {
                if (nodes[nodeId].type == SequenceType::ORGANISM) {
                    float score = scoreOrganismAttribution(geneNode, nodeId);
                    organismScores.push_back({nodes[nodeId].sequenceId, score});
                }
            }
            
            // Sort by score and take top candidates
            std::sort(organismScores.begin(), organismScores.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            for (const auto& [organismId, score] : organismScores) {
                if (score > 0.5) {  // Confidence threshold
                    attr.candidateOrganisms.push_back({organismId, score});
                }
            }
            
            attributions.push_back(attr);
        }
        
        return attributions;
    }
    
private:
    float scoreOrganismAttribution(int geneNode, int organismNode) {
        // Find shortest path and calculate path score
        std::vector<int> path = findShortestPath(geneNode, organismNode);
        if (path.empty()) return 0.0f;
        
        float pathScore = 0.0f;
        float totalWeight = 0.0f;
        
        for (int nodeId : path) {
            float weight = nodes[nodeId].kmerWeight;
            float coverage = nodes[nodeId].coverage;
            float identity = nodes[nodeId].identity;
            
            pathScore += weight * coverage * identity;
            totalWeight += weight;
        }
        
        return totalWeight > 0 ? pathScore / totalWeight : 0.0f;
    }
    
    std::vector<int> findShortestPath(int start, int end) {
        if (start == end) return {start};
        
        std::queue<int> queue;
        std::unordered_map<int, int> parent;
        std::unordered_map<int, bool> visited;
        
        queue.push(start);
        visited[start] = true;
        parent[start] = -1;
        
        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();
            
            if (current == end) {
                // Reconstruct path
                std::vector<int> path;
                int node = end;
                while (node != -1) {
                    path.push_back(node);
                    node = parent[node];
                }
                std::reverse(path.begin(), path.end());
                return path;
            }
            
            for (int neighbor : nodes[current].neighbors) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    parent[neighbor] = current;
                    queue.push(neighbor);
                }
            }
        }
        
        return {};  // No path found
    }
};

// Public interface implementation
ResistanceAttribution::ResistanceAttribution() 
    : pImpl(std::make_unique<ResistanceAttributionImpl>()) {}

ResistanceAttribution::~ResistanceAttribution() = default;

void ResistanceAttribution::buildAssemblyGraph(
    const std::vector<SequenceFragment>& fragments,
    const std::vector<AlignmentResult>& alignments) {
    pImpl->buildGraph(fragments, alignments);
}

std::vector<ResistanceAttribution> ResistanceAttribution::attributeResistance(
    const std::vector<std::string>& resistanceGeneIds) {
    return pImpl->attributeResistanceGenes(resistanceGeneIds);
}

std::vector<PathScore> ResistanceAttribution::scoreAllPaths(
    const std::string& geneId,
    float minScore) {
    // Implementation would score all paths from gene to organisms
    std::vector<PathScore> scores;
    // ... implementation details ...
    return scores;
}

} // namespace algorithms
} // namespace biogpu