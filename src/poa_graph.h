//
// Created by Arthur Rand on 10/6/16.
//

#ifndef PARCOURS_POAGRAPH_H
#define PARCOURS_POAGRAPH_H

#include "stl_includes.h"
#include "poa_exceptions.h"
#include "vertex.h"


typedef struct _sequence {
    std::string seq;
    std::string label;
} Sequence;

class PoaGraph {
public:
    PoaGraph(): nVertices(0), nb_arcs(0), next_vertex_id(0) {};

    PoaGraph(const Sequence &seq);

    // todo
    ~PoaGraph() {}

    void AddBaseSequence(const Sequence &seq);

    int64_t AddVertex(char base);

    bool ContainsVertex(int64_t i);

    bool ContainsVertex(Vertex& v);

    void AddArc(int64_t startId, int64_t endId, std::string label);

    int64_t K();

    Vertex *VertexGetter(int64_t i);

    unsigned long VertexOutDegree(int i) {
        return vertex_map[i]->OutDegree();
    }

    unsigned long VertexInDegree(int i) {
        return vertex_map[i]->InDegree();
    }

    std::vector<int64_t>& Vertices();



    std::vector<int64_t>& Starts();

    std::unordered_map<int64_t, DirectedArc*>& OutNeighbors(int64_t);

    void TopologicalSort();

    bool TestSort();

    bool isSorted();
private:
    // containers
    std::unordered_map<int64_t, std::set<int>> adjacentcyList;  // map vertexId -> vertexId
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

#endif //PARCOURS_POAGRAPH_H
