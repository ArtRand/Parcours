//
// Created by Arthur Rand on 10/8/16.
//

#ifndef PARCOURS_VERTEX_H
#define PARCOURS_VERTEX_H

#include "stl_includes.h"
#include "directed_arc.h"

/* Vertex class.
 * Basic class for vertices (or nodes) in a DAG. 
 */
class Vertex {
public:
    // out_arcs and in_arcs: keys are ajacent verices by ID, values are the connecting 
    // directed arcs (for storing things like labels and weights)
    std::unordered_map<int64_t, DirectedArc*> out_arcs;
    std::unordered_map<int64_t, DirectedArc*> in_arcs;

    // aligned to, other vertex IDs that this node is aligned to, for comparing bases
    std::vector<int64_t> aligned_to;

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
    
    // all Add[In/Out]Arc functions return true if the arc *DIDN'T* exist already, false if we just update the labels
    // TODO figure out how to hide this function
    bool AddArc(std::unordered_map<int64_t, DirectedArc*>& arcs, int64_t neighborId, std::string label); 
   
    bool AddInArc(int64_t nId, std::string label);

    bool AddOutArc(int64_t nId, std::string label);

    unsigned long InDegree();

    unsigned long OutDegree();

    std::unordered_map<int64_t, DirectedArc*>& OutNeighbors();
    // TODO make InNeighbors,

private:
    int64_t id;
    char base;
};


#endif //PARCOURS_VERTEX_H
