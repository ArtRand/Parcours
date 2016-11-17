//
// Created by Arthur Rand on 10/24/16
//

#ifndef PARCOURS_HMM_GRAPH_H
#define PARCOURS_HMM_GRAPH_H

#include "stl_includes.h"
#include "symbol_string.h"
#include "common.h"

#include "vertex.h"

// intermediate graph structure
class HmmGraph {
public:
    // TODO move these to cpp, for consistency 
    HmmGraph();
    
    ~HmmGraph() {}; 
   
    HmmGraph(HmmGraph& other); 

    HmmGraph& operator = (HmmGraph& other);

    HmmGraph(HmmGraph&& other) = default;
    
    HmmGraph& operator = (HmmGraph&& other) = default;

    int64_t AddVertex(const std::string *seq);

    bool ContainsVertex(int64_t i) const;

    bool ContainsVertex(Vertex *v) const;

    void AddArc(int64_t, int64_t);
    
    void SetNextVertexId(int64_t i);

    int64_t K();

    Vertex *VertexGetter(int64_t i);  // caution: returns raw pointer to vertex

    unsigned long VertexOutDegree(int64_t i);

    unsigned long VertexInDegree(int64_t i);

    const std::string *VertexSequence(int64_t i);
    
    const std::set<int64_t>& OutNeighbors(int64_t i);

    const std::set<int64_t>& InNeighbors(int64_t i);

    const std::vector<int64_t>& Vertices() const;

    void TopologicalSort(bool test=false);

    bool TestSort();

    bool IsSorted();

    std::set<int64_t> Sources();

    std::set<int64_t> Sinks();

    std::vector<std::deque<int64_t>> AllPaths();
    
    std::vector<SymbolString> ExtractSequences(const VertexPaths& vertex_paths);
    //friend HmmGraph;

    //friend Vertex;

private:
    // containers
    std::unordered_map<int64_t, std::set<int64_t>> adjacentcy_list;
    std::unordered_map<int64_t, std::unique_ptr<Vertex>> vertex_map;  // for looking up vertices
    std::vector<int64_t> vertex_list;  // for ordering
    //std::set<std::string> labels;
    //std::set<std::string *> seqs;
    //std::vector<int64_t> starts;
    std::vector<std::deque<int64_t>> paths;
    std::unordered_map<int64_t, uint64_t> path_lookup;  // (pathID, paths.at(path))  
    std::unordered_map<int64_t, double> path_scores;    // (pathID, score)

    // counters and flags
    int64_t nVertices;
    int64_t nArcs;
    int64_t next_vertex_id;
    int64_t nPaths;
    bool sorted = false;
    bool initialized_paths = false;

    // internal functions
    void find_paths();
    void copy_graph(HmmGraph& orig, HmmGraph& other);
    void clear_graph();
};

#endif
