//
// Created by Arthur Rand on 10/6/16.
//

#include <iostream>
#include "poa_graph.h"


PoaGraph::PoaGraph(const Sequence& seq) : nVertices(0), nb_arcs(0), next_vertex_id(0)  {
    AddBaseSequence(seq);
}

bool PoaGraph::ContainsVertex(int64_t i) {
    if (i > nVertices || i >= next_vertex_id) {
        return false;
    }

    std::vector<int64_t>::iterator it = std::find(vertex_list.begin(), vertex_list.end(), i);
    return it != vertex_list.end();
}

bool PoaGraph::ContainsVertex(Vertex& v) {
    return ContainsVertex(v.Id());
}

int64_t PoaGraph::AddVertex(char base) {
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

Vertex *PoaGraph::VertexGetter(int i) {
    if (!ContainsVertex(i)) {
        throw GraphException("VertexGetter: Could not find vertex");
    }
    return vertex_map[i];
}

void PoaGraph::AddArc(int64_t startId, int64_t endId, std::string label) {
    if (startId < 0 || endId < 0) {
        return;
    }

    if (!ContainsVertex(startId) || !ContainsVertex(endId)) {
        throw GraphException("AddArc: Vertex Contains error, don't have vertex");
    }

    bool added_in_arc = vertex_map[startId]->AddOutArc(endId, label);
    bool added_out_arc = vertex_map[endId]->AddInArc(startId, label);

    // if we added in or out arcs increment the number of arcs in the graph (the counter)

    if (added_in_arc || added_out_arc) {
        nb_arcs += 1;
    }

    sorted = false;
}

void PoaGraph::AddBaseSequence(const Sequence& seq) {
    int64_t firstId, prevId;
    firstId = -1;
    prevId = -1;
    bool need_sort = sorted;

    for (char c : seq.seq) {
        int64_t vId = AddVertex(c);
        // condition for first character, make it the first vertex-id
        if (firstId < 0) {
            firstId = vId;
        }
        // condition for adding all nodes after the first, connect them to previous vertex
        if (prevId >= 0) {
            AddArc(prevId, vId, seq.label);
        }
        prevId = vId;  // update
    }

    sorted = need_sort;

    // todo might need logic here
    seqs.insert(seq.seq);
    labels.insert(seq.label);
    starts.push_back(firstId);
}

// helper methods
int64_t PoaGraph::K() { return nVertices; }
std::vector<int64_t>& PoaGraph::Starts() { return starts; }
std::vector<int64_t>& PoaGraph::Vertices() { return vertex_list; }