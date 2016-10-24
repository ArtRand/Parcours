// hmm_graph.cpp

#include "hmm_graph.h"


Vertex::Vertex(int64_t i, char b) {
    id = i;
    base = b;
}

int64_t Vertex::Id() const {
    return id;
}

int64_t HmmGraph::AddVertex(char base) {
    // get the next vertex id, a monotonically increasing int
    int64_t vId = next_vertex_id;

    Vertex *v = new Vertex(vId, base);

    // check, todo remove after tested
    assert(v->Id() == vId);

    // updates to state of graph
    vertex_map[vId] = v;
    vertex_list.push_back(vId);
    nVertices += 1;
    next_vertex_id += 1;
    sorted = false;

    return vId;
}

void HmmGraph::AddArc(int64_t fromId, int64_t toId) {
    adjacentcy_list[fromId].insert(toId);
}

int64_t HmmGraph::K() { return nVertices; }


