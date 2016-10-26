//
// Created by Arthur Rand on 10/8/16.
//

#include "vertex.h"

Vertex::Vertex(int64_t i) {
    id = i;
    node_sequence = nullptr;
}

Vertex::Vertex(int64_t i, std::string *seq) {
    if (seq->length() <= 0) {
        throw ParcoursException("Initializing vertex with empty sequence\n");
    }
    id = i;
    node_sequence = seq;
}

/*
Vertex::Vertex(const Vertex& other) {
    if (&other != this) {
        in_neighbors = other.InNeighbors();
        out_neighbors = other.OutNeighbors();
        id = other.Id();
        node_sequence = other.Sequence();
    }
}

Vertex& Vertex::operator= (const Vertex& other) {
    in_neighbors = other.InNeighbors();
    out_neighbors = other.OutNeighbors();
    id = other.Id();
    node_sequence = other.Sequence();
    return *this;
}

Vertex::Vertex(Vertex&& other): in_neighbors(other.InNeighbors()), out_neighbors(other.OutNeighbors()),
                                id(other.Id()), node_sequence(other.Sequence()) {}

*/

bool Vertex::operator< (Vertex other) const { return id == other.id; }

void Vertex::AddInNeighbor(int64_t i) { in_neighbors.insert(i); }

void Vertex::AddOutNeighbor(int64_t i) { out_neighbors.insert(i); }

bool Vertex::IsInNeighbor(int64_t n) { return in_neighbors.count(n) > 0; }

bool Vertex::IsOutNeighbor(int64_t n) { return out_neighbors.count(n) > 0; }

const std::set<int64_t>& Vertex::InNeighbors() const { return in_neighbors; }

const std::set<int64_t>& Vertex::OutNeighbors() const { return out_neighbors; }

void Vertex::SetId(int64_t i) { id = i; }

int64_t Vertex::Id() const { return id; }

std::string *Vertex::Sequence() const { return node_sequence; }

unsigned long Vertex::InDegree() { return in_neighbors.size(); }

unsigned long Vertex::OutDegree() { return out_neighbors.size(); }

