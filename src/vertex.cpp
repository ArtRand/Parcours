//
// Created by Arthur Rand on 10/8/16.
//

#include "vertex.h"

Vertex::Vertex(int64_t i) {
    id = i;
};

Vertex::Vertex(int64_t i, char b) {
    id = i;
    base = b;
}

void Vertex::SetId(int64_t i) {
    id = i;
}

int64_t Vertex::Id() const {
    return id;
}

bool Vertex::operator < (Vertex other) const {
    return id == other.id;
}

// returns true if we added the arc, false if it already existed and we just updated it
bool Vertex::AddArc(std::unordered_map<int64_t, DirectedArc*>& arcs, int64_t neighborId, std::string label) {
    // check if this vertex already contains an arc to neighbor, in that case just update the arc label
    if (arcs.count(neighborId) > 0) {
        arcs[neighborId]->AddLabel(label);
        return false;
    }
    DirectedArc *a = new DirectedArc(id, neighborId, label);
    arcs[neighborId] = a;
    // todo remove after testing
    assert(arcs.count(neighborId) == 1);
    return true;
}

bool Vertex::AddInArc(int64_t nId, std::string label) {
    return AddArc(in_arcs, nId, label);
}

bool Vertex::AddOutArc(int64_t nId, std::string label) {
    return AddArc(out_arcs, nId, label);
}

unsigned long Vertex::InDegree() {
    return in_arcs.size();
}

unsigned long Vertex::OutDegree() {
    return out_arcs.size();
}
