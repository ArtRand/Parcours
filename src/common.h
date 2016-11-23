//
// Created by Arthur Rand on 10/13/16.
//

#ifndef PARCOURS_COMMON_H
#define PARCOURS_COMMON_H

#include "stl_includes.h"

typedef enum _hidden_state {
    match = 0,
    shortGapX = 1,
    shortGapY = 2,
    longGapX = 3,
    longGapY = 4,
} HiddenState;

typedef enum _set_types { 
    nucleotide = 4,
} SetType;

typedef enum _statemachine_type {
    fiveState = 5,
} StateMachineType;


typedef std::vector<std::pair<int64_t, int64_t>> AnchorPairs;
typedef std::tuple<double, int64_t, int64_t> AlignedPair;  // (prob, x, y)
typedef std::vector<AlignedPair> AlignedPairs;

// a type for storing aligned pairs to a graph structure (prob, x, (vertex, offset))
typedef std::tuple<double, int64_t, std::pair<int64_t, int64_t>> GraphAlignedPair; 
typedef std::vector<GraphAlignedPair> GraphAlignedPairs;

typedef std::vector<std::deque<int64_t>> VertexPaths;

void st_uglyf(const char *string, ...);

AnchorPairs EmptyAnchors();

#endif //PARCOURS_COMMON_H
