//
// Created by Arthur Rand on 10/24/16
//

#ifndef PARCOURS_HMM_GRAPH_H
#define PARCOURS_HMM_GRAPH_H

#include "stl_includes.h"
#include "vertex.h"
#include "common.h"

// intermediate graph structure
class HmmGraph {
public:
    HmmGraph(): nVertices(0), nArcs(0), next_vertex_id(0) {};
    
    ~HmmGraph() {}; 
   
    HmmGraph(HmmGraph& other) {
        copy_graph(*this, other);
    }

    HmmGraph& operator = (HmmGraph& other) {
        if (&other != this) {
            copy_graph(*this, other);
        }
        return *this;
    }

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
    
    //friend HmmGraph;

    friend Vertex;

private:
    // containers
    std::unordered_map<int64_t, std::set<int64_t>> adjacentcy_list;
    std::unordered_map<int64_t, std::unique_ptr<Vertex>> vertex_map;  // for looking up vertices
    std::vector<int64_t> vertex_list;  // for ordering
    //std::set<std::string> labels;
    //std::set<std::string *> seqs;
    std::vector<int64_t> starts;
    std::vector<std::deque<int64_t>> paths;

    // counters and flags
    int64_t nVertices;
    int64_t nArcs;
    int64_t next_vertex_id;
    bool sorted = false;
    bool initialized_paths = false;

    // internal functions
    void find_paths();
    void copy_graph(HmmGraph& orig, HmmGraph& other);
};

#endif
