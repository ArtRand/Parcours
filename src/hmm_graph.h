//
// Created by Arthur Rand on 10/24/16
//

#ifndef PARCOURS_HMM_GRAPH_H
#define PARCOURS_HMM_GRAPH_H

/* Vertex class.
 * Basic class for vertices (aka nodes) in a DAG. 
 */

#include "stl_includes.h"
#include "common.h"

class Vertex {
public:
    // constructors
    Vertex(): id(-1), base('\0') {};  // default construct with ID = -1
    Vertex(int64_t i);
    Vertex(int64_t i, char b);
    
    // destructor TODO test for leak 
    ~Vertex() {};

    void SetId(int64_t i);  // shouldn't ever need this

    int64_t Id() const;

    char Base() const;

    bool operator < (Vertex other) const;

private:
    int64_t id;
    char base;
};

// intermediate graph structure
class HmmGraph {
public:
    HmmGraph(): nVertices(0), nb_arcs(0), next_vertex_id(0) {st_uglyf("-->Made HmmGraph!!\n");};

    ~HmmGraph() {};

    int64_t AddVertex(char base);

    bool ContainsVertex(int64_t i);

    bool ContainsVertex(Vertex& v);

    void AddArc(int64_t, int64_t);

    int64_t K();

    Vertex *VertexGetter(int64_t i);

    unsigned long VertexOutDegree(int64_t i);

    unsigned long VertexInDegree(int64_t i);

    std::vector<int64_t>& Vertices();

    void TopologicalSort();

    bool TestSort();

    bool IsSorted();

private:
    // containers
    std::unordered_map<int64_t, std::set<int64_t>> adjacentcy_list;
    std::unordered_map<int64_t, Vertex*> vertex_map;  // for looking up vertices
    std::vector<int64_t> vertex_list;  // for ordering
    std::set<std::string> labels;
    std::set<std::string> seqs;
    std::vector<int64_t> starts;

    // counters and flags
    int64_t nVertices;
    int64_t nb_arcs;
    int64_t next_vertex_id;
    bool sorted = false;
};

#endif
