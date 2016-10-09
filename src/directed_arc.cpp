//
// Created by Arthur Rand on 10/8/16.
//

#include "directed_arc.h"

DirectedArc::DirectedArc(int64_t t_id, int64_t h_id)  {
    head_vertex_id = h_id;
    tail_vertex_id = t_id;
}

DirectedArc::DirectedArc(int64_t t_id, int64_t h_id, std::string label) {
    head_vertex_id = h_id;
    tail_vertex_id = t_id;
    labels.insert(label);
}

int64_t DirectedArc::To() const {
    return head_vertex_id;
};

int64_t DirectedArc::From() const {
    return tail_vertex_id;
};

std::set<std::string> DirectedArc::LabelSet() const {
    return labels;
}

// returns true if the label was already in the set, otherwise it gets added and returns false
bool DirectedArc::AddLabel(std::string label) {
    if (labels.count(label)) {
        return true;
    } else {
        labels.insert(label);
        return false;
    }
}

bool DirectedArc::operator == (DirectedArc other) {
    return (head_vertex_id == other.To()) && (tail_vertex_id == other.From());
}
