//
// Created by Arthur Rand on 10/8/16.
//

#ifndef PARCOURS_VERTEX_H
#define PARCOURS_VERTEX_H

#include "stl_includes.h"
#include "parcours_exceptions.h"
/*
 * Vertex class.
 * Basic class for vertices (or nodes) in a DAG, intentionally left lightweight and extensible
 */
class Vertex {
public:
    // out_arcs and in_arcs: keys are adjacent verices by ID, values are the connecting 
    // directed arcs (for storing things like labels and weights)
    //std::unordered_map<int64_t, DirectedArc*> out_arcs;
    //std::unordered_map<int64_t, DirectedArc*> in_arcs;

    // aligned to, other vertex IDs that this node is aligned to, for comparing bases
    //std::vector<int64_t> aligned_to;

    // constructors
    Vertex(int64_t i);
    Vertex(int64_t i, const std::string *seq);
    Vertex(const Vertex& other) = default;
    Vertex& operator = (const Vertex& other) = default;
    Vertex(Vertex&& other) = default;
    Vertex& operator = (Vertex&& other) = default;
    
    //Vertex(const Vertex& other);
    //Vertex& operator= (const Vertex& other);
    //Vertex(Vertex&& other);
    ~Vertex() {};

    bool operator < (Vertex other) const;  // tests for equality in vertex ID, does not respect ordering

    void AddInNeighbor(int64_t i);  // add a vertex ID to in_neighbors

    void AddOutNeighbor(int64_t i);  // add a vertex ID to out_neighbors

    bool IsInNeighbor(int64_t n);  // returns true iff n is in in_neighbors

    bool IsOutNeighbor(int64_t n);  // returns true iff n is in out_neighbors

    const std::set<int64_t>& InNeighbors() const;

    const std::set<int64_t>& OutNeighbors() const;

    void SetId(int64_t i);  // shouldn't ever need this

    int64_t Id() const;

    const std::string *Sequence() const;

    // TODO not sure if I need this anymore
    //bool AddArc(std::unordered_map<int64_t, DirectedArc*>& arcs, int64_t neighborId, std::string label);
   
    //bool AddInArc(int64_t nId, std::string label);

    //bool AddOutArc(int64_t nId, std::string label);

    unsigned long InDegree();

    unsigned long OutDegree();

    //std::unordered_map<int64_t, DirectedArc*>& OutNeighbors();
    // TODO make InNeighbors,

private:
    // keep track of adjacent nodes for all paths calculation
    std::set<int64_t> in_neighbors;
    std::set<int64_t> out_neighbors;

    int64_t id;
    const std::string *node_sequence;
};


#endif //PARCOURS_VERTEX_H
