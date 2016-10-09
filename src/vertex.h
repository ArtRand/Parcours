//
// Created by Arthur Rand on 10/8/16.
//

#ifndef PARCOURS_VERTEX_H
#define PARCOURS_VERTEX_H

#include "stl_includes.h"
#include "directed_arc.h"

class Vertex {
public:
    std::unordered_map<int64_t, DirectedArc*> out_arcs;
    std::unordered_map<int64_t, DirectedArc*> in_arcs;
    std::vector<Vertex*> aligned_to;

    Vertex(): id(-1), base('\0') {};
    Vertex(int64_t i);

    Vertex(int64_t i, char b);

    ~Vertex() {};

    void SetId(int64_t i);

    int64_t Id() const;

    bool operator < (Vertex other) const;

    bool AddArc(std::unordered_map<int64_t, DirectedArc*>& arcs, int64_t neighborId, std::string label);

    bool AddInArc(int64_t nId, std::string label);

    bool AddOutArc(int64_t nId, std::string label);

    unsigned long InDegree();

    unsigned long OutDegree();

private:
    int64_t id;
    char base;
};


#endif //PARCOURS_VERTEX_H
