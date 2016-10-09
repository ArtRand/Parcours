//
// Created by Arthur Rand on 10/8/16.
//

#ifndef PARCOURS_DIRECTEDARC_H
#define PARCOURS_DIRECTEDARC_H

#include "stl_includes.h"


class DirectedArc {
    int64_t head_vertex_id;
    int64_t tail_vertex_id;
    std::set<std::string> labels;
public:
    DirectedArc(): head_vertex_id(-1), tail_vertex_id(-1) {};

    DirectedArc(int64_t t_id, int64_t h_id);

    DirectedArc(int64_t t_id, int64_t h_id, std::string label);

    ~DirectedArc() {};

    int64_t To() const;

    int64_t From() const;

    std::set<std::string> LabelSet() const;

    // returns true if the label was already in the set, otherwise it gets added and returns false
    bool AddLabel(std::string label);

    bool operator == (DirectedArc other);
};
#endif //PARCOURS_DIRECTEDARC_H
