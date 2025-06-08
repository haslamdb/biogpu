// resistance_graph.cpp
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <queue>
#include <set>
#include <cassert>
#include <fstream>
// #include <json/json.h>  // Would use nlohmann/json in practice
// For now, we'll implement without JSON support

// Node types in the resistance graph
enum class NodeType {
    ORGANISM,
    GENE,
    MUTATION,
    DRUG,
    RESISTANCE_MECHANISM
};

// Edge types representing relationships
enum class EdgeType {
    HAS_GENE,           // Organism -> Gene
    CONTAINS_MUTATION,  // Gene -> Mutation
    CONFERS_RESISTANCE, // Mutation -> Drug
    CO_OCCURS,         // Mutation -> Mutation
    MECHANISM,         // Mutation -> Resistance Mechanism
    AFFECTS            // Resistance Mechanism -> Drug
};

// Base node structure
struct Node {
    uint32_t id;
    NodeType type;
    std::string name;
    std::unordered_map<std::string, std::string> properties;
    
    virtual ~Node() = default;
};

// Specialized node types
struct OrganismNode : public Node {
    std::string taxonomy;
    std::string strain;
    float typical_abundance;  // In healthy microbiome
    
    OrganismNode() { type = NodeType::ORGANISM; }
};

struct GeneNode : public Node {
    std::string gene_family;
    std::string function;
    std::vector<std::pair<int, int>> conserved_regions;
    std::vector<std::pair<int, int>> qrdr_regions;  // Quinolone resistance-determining regions
    
    GeneNode() { type = NodeType::GENE; }
};

struct MutationNode : public Node {
    uint32_t gene_id;
    int position;
    char wild_type;
    char mutant;
    std::string codon_change;  // e.g., "TCT->TTT"
    std::string aa_change;     // e.g., "S83F"
    float fitness_cost;        // 0.0 = no cost, 1.0 = lethal
    
    MutationNode() { type = NodeType::MUTATION; }
};

struct DrugNode : public Node {
    std::string drug_class;
    std::string mechanism_of_action;
    float typical_mic;  // Typical MIC for susceptible strains
    
    DrugNode() { type = NodeType::DRUG; }
};

struct ResistanceMechanismNode : public Node {
    std::string mechanism_type;  // "target_modification", "efflux", "enzymatic"
    std::string description;
    
    ResistanceMechanismNode() { type = NodeType::RESISTANCE_MECHANISM; }
};

// Edge structure
struct Edge {
    uint32_t source_id;
    uint32_t target_id;
    EdgeType type;
    float weight;  // Strength/confidence of relationship
    std::unordered_map<std::string, std::string> properties;
};

// Main resistance graph class
class ResistanceGraph {
private:
    std::unordered_map<uint32_t, std::shared_ptr<Node>> nodes;
    std::vector<Edge> edges;
    std::unordered_map<uint32_t, std::vector<size_t>> adjacency_list;  // Node ID -> edge indices
    uint32_t next_node_id = 0;
    
    // Indices for fast lookup
    std::unordered_map<std::string, uint32_t> name_to_id;
    std::unordered_multimap<NodeType, uint32_t> type_to_ids;
    std::unordered_map<std::string, std::vector<uint32_t>> gene_to_mutations;

public:
    // Add nodes
    uint32_t add_organism(const std::string& name, const std::string& taxonomy) {
        auto node = std::make_shared<OrganismNode>();
        node->id = next_node_id++;
        node->name = name;
        node->taxonomy = taxonomy;
        
        nodes[node->id] = node;
        name_to_id[name] = node->id;
        type_to_ids.insert({NodeType::ORGANISM, node->id});
        
        return node->id;
    }
    
    uint32_t add_gene(const std::string& name, const std::string& function) {
        auto node = std::make_shared<GeneNode>();
        node->id = next_node_id++;
        node->name = name;
        node->function = function;
        
        nodes[node->id] = node;
        name_to_id[name] = node->id;
        type_to_ids.insert({NodeType::GENE, node->id});
        
        return node->id;
    }
    
    uint32_t add_mutation(const std::string& gene_name, int position, 
                         char wild_type, char mutant, const std::string& aa_change) {
        auto node = std::make_shared<MutationNode>();
        node->id = next_node_id++;
        node->name = gene_name + "_" + aa_change;
        node->position = position;
        node->wild_type = wild_type;
        node->mutant = mutant;
        node->aa_change = aa_change;
        
        // Link to gene
        if (name_to_id.count(gene_name)) {
            node->gene_id = name_to_id[gene_name];
            gene_to_mutations[gene_name].push_back(node->id);
        }
        
        nodes[node->id] = node;
        name_to_id[node->name] = node->id;
        type_to_ids.insert({NodeType::MUTATION, node->id});
        
        return node->id;
    }
    
    uint32_t add_drug(const std::string& name, const std::string& drug_class) {
        auto node = std::make_shared<DrugNode>();
        node->id = next_node_id++;
        node->name = name;
        node->drug_class = drug_class;
        
        nodes[node->id] = node;
        name_to_id[name] = node->id;
        type_to_ids.insert({NodeType::DRUG, node->id});
        
        return node->id;
    }
    
    // Add edges
    void add_edge(uint32_t source_id, uint32_t target_id, EdgeType type, float weight = 1.0) {
        Edge edge;
        edge.source_id = source_id;
        edge.target_id = target_id;
        edge.type = type;
        edge.weight = weight;
        
        size_t edge_idx = edges.size();
        edges.push_back(edge);
        adjacency_list[source_id].push_back(edge_idx);
    }
    
    // Query functions
    std::vector<uint32_t> find_resistance_path(uint32_t organism_id, const std::string& drug_name) {
        std::vector<uint32_t> path;
        if (!name_to_id.count(drug_name)) return path;
        
        uint32_t drug_id = name_to_id[drug_name];
        
        // BFS to find path from organism to drug through mutations
        std::queue<std::vector<uint32_t>> queue;
        std::set<uint32_t> visited;
        
        queue.push({organism_id});
        visited.insert(organism_id);
        
        while (!queue.empty()) {
            auto current_path = queue.front();
            queue.pop();
            
            uint32_t current_node = current_path.back();
            
            if (current_node == drug_id) {
                return current_path;
            }
            
            // Explore neighbors
            if (adjacency_list.count(current_node)) {
                for (size_t edge_idx : adjacency_list[current_node]) {
                    const Edge& edge = edges[edge_idx];
                    if (visited.find(edge.target_id) == visited.end()) {
                        visited.insert(edge.target_id);
                        auto new_path = current_path;
                        new_path.push_back(edge.target_id);
                        queue.push(new_path);
                    }
                }
            }
        }
        
        return path;  // Empty if no path found
    }
    
    // Find co-occurring mutations
    std::vector<uint32_t> find_co_occurring_mutations(uint32_t mutation_id) {
        std::vector<uint32_t> co_occurring;
        
        if (adjacency_list.count(mutation_id)) {
            for (size_t edge_idx : adjacency_list[mutation_id]) {
                const Edge& edge = edges[edge_idx];
                if (edge.type == EdgeType::CO_OCCURS) {
                    co_occurring.push_back(edge.target_id);
                }
            }
        }
        
        return co_occurring;
    }
    
    // Calculate resistance score for an organism
    float calculate_resistance_score(uint32_t organism_id, const std::string& drug_name) {
        float score = 0.0;
        
        // Find all mutations in this organism
        std::vector<uint32_t> organism_mutations;
        
        if (adjacency_list.count(organism_id)) {
            for (size_t edge_idx : adjacency_list[organism_id]) {
                const Edge& edge = edges[edge_idx];
                if (edge.type == EdgeType::HAS_GENE) {
                    uint32_t gene_id = edge.target_id;
                    
                    // Find mutations in this gene
                    if (adjacency_list.count(gene_id)) {
                        for (size_t mut_edge_idx : adjacency_list[gene_id]) {
                            const Edge& mut_edge = edges[mut_edge_idx];
                            if (mut_edge.type == EdgeType::CONTAINS_MUTATION) {
                                organism_mutations.push_back(mut_edge.target_id);
                            }
                        }
                    }
                }
            }
        }
        
        // Check which mutations confer resistance to the drug
        uint32_t drug_id = name_to_id[drug_name];
        for (uint32_t mut_id : organism_mutations) {
            if (adjacency_list.count(mut_id)) {
                for (size_t edge_idx : adjacency_list[mut_id]) {
                    const Edge& edge = edges[edge_idx];
                    if (edge.type == EdgeType::CONFERS_RESISTANCE && 
                        edge.target_id == drug_id) {
                        score += edge.weight;  // Weight represents resistance strength
                    }
                }
            }
        }
        
        return score;
    }
    
    // Build fluoroquinolone resistance subgraph
    void build_fluoroquinolone_resistance_graph() {
        // Add organisms
        uint32_t ecoli_id = add_organism("E.coli", "Enterobacteriaceae");
        uint32_t kpneumo_id = add_organism("K.pneumoniae", "Enterobacteriaceae");
        uint32_t paeru_id = add_organism("P.aeruginosa", "Pseudomonadaceae");
        
        // Add genes
        uint32_t gyrA_id = add_gene("gyrA", "DNA gyrase subunit A");
        uint32_t gyrB_id = add_gene("gyrB", "DNA gyrase subunit B");
        uint32_t parC_id = add_gene("parC", "Topoisomerase IV subunit A");
        uint32_t parE_id = add_gene("parE", "Topoisomerase IV subunit B");
        
        // Add common mutations
        uint32_t gyrA_S83L = add_mutation("gyrA", 83, 'S', 'L', "S83L");
        uint32_t gyrA_D87N = add_mutation("gyrA", 87, 'D', 'N', "D87N");
        uint32_t parC_S80I = add_mutation("parC", 80, 'S', 'I', "S80I");
        uint32_t parC_E84V = add_mutation("parC", 84, 'E', 'V', "E84V");
        
        // Add drugs
        uint32_t cipro_id = add_drug("ciprofloxacin", "fluoroquinolone");
        uint32_t levo_id = add_drug("levofloxacin", "fluoroquinolone");
        
        // Add relationships
        // Organisms have genes
        add_edge(ecoli_id, gyrA_id, EdgeType::HAS_GENE);
        add_edge(ecoli_id, parC_id, EdgeType::HAS_GENE);
        add_edge(kpneumo_id, gyrA_id, EdgeType::HAS_GENE);
        add_edge(kpneumo_id, parC_id, EdgeType::HAS_GENE);
        
        // Genes contain mutations
        add_edge(gyrA_id, gyrA_S83L, EdgeType::CONTAINS_MUTATION);
        add_edge(gyrA_id, gyrA_D87N, EdgeType::CONTAINS_MUTATION);
        add_edge(parC_id, parC_S80I, EdgeType::CONTAINS_MUTATION);
        add_edge(parC_id, parC_E84V, EdgeType::CONTAINS_MUTATION);
        
        // Mutations confer resistance
        add_edge(gyrA_S83L, cipro_id, EdgeType::CONFERS_RESISTANCE, 8.0);  // 8-fold MIC increase
        add_edge(gyrA_S83L, levo_id, EdgeType::CONFERS_RESISTANCE, 4.0);
        add_edge(gyrA_D87N, cipro_id, EdgeType::CONFERS_RESISTANCE, 4.0);
        add_edge(parC_S80I, cipro_id, EdgeType::CONFERS_RESISTANCE, 4.0);
        
        // Co-occurring mutations
        add_edge(gyrA_S83L, parC_S80I, EdgeType::CO_OCCURS, 0.7);  // Often found together
    }
    
    // Visualize graph (outputs DOT format)
    void export_dot(const std::string& filename) {
        std::ofstream out(filename);
        out << "digraph ResistanceGraph {\n";
        out << "  rankdir=LR;\n";
        
        // Define node styles
        out << "  node [shape=box];\n";
        
        // Write nodes
        for (const auto& [id, node] : nodes) {
            std::string color = "white";
            std::string shape = "box";
            
            switch(node->type) {
                case NodeType::ORGANISM: color = "lightblue"; shape = "ellipse"; break;
                case NodeType::GENE: color = "lightgreen"; break;
                case NodeType::MUTATION: color = "yellow"; shape = "diamond"; break;
                case NodeType::DRUG: color = "pink"; shape = "ellipse"; break;
                default: break;
            }
            
            out << "  " << id << " [label=\"" << node->name 
                << "\", style=filled, fillcolor=" << color 
                << ", shape=" << shape << "];\n";
        }
        
        // Write edges
        for (const Edge& edge : edges) {
            std::string style = "solid";
            std::string color = "black";
            
            switch(edge.type) {
                case EdgeType::HAS_GENE: color = "blue"; break;
                case EdgeType::CONTAINS_MUTATION: color = "green"; break;
                case EdgeType::CONFERS_RESISTANCE: color = "red"; style = "bold"; break;
                case EdgeType::CO_OCCURS: color = "purple"; style = "dashed"; break;
                default: break;
            }
            
            out << "  " << edge.source_id << " -> " << edge.target_id 
                << " [color=" << color << ", style=" << style;
            
            if (edge.weight != 1.0) {
                out << ", label=\"" << edge.weight << "\"";
            }
            
            out << "];\n";
        }
        
        out << "}\n";
        out.close();
    }
};

// Example usage
#ifdef TEST
int main() {
    ResistanceGraph graph;
    
    // Build a sample fluoroquinolone resistance graph
    graph.build_fluoroquinolone_resistance_graph();
    
    // Export for visualization
    graph.export_dot("resistance_graph.dot");
    std::cout << "Graph exported to resistance_graph.dot\n";
    std::cout << "Visualize with: dot -Tpng resistance_graph.dot -o resistance_graph.png\n\n";
    
    // Test queries
    uint32_t ecoli_id = 0;  // First organism added
    
    // Find resistance path
    auto path = graph.find_resistance_path(ecoli_id, "ciprofloxacin");
    std::cout << "Resistance path from E.coli to ciprofloxacin:\n";
    for (uint32_t node_id : path) {
        std::cout << "  -> Node " << node_id << "\n";
    }
    
    // Calculate resistance score
    float score = graph.calculate_resistance_score(ecoli_id, "ciprofloxacin");
    std::cout << "\nE.coli resistance score for ciprofloxacin: " << score << "\n";
    
    return 0;
}
#endif // TEST