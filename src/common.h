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

typedef enum {
    a=0,
    c=1,
    g=2,
    t=3,
    n=4,
} Symbol;

typedef std::vector<std::pair<int64_t, int64_t>> AnchorPairs;
typedef std::vector<std::tuple<double, int64_t, int64_t>> AlignedPairs;

typedef std::vector<Symbol> SymbolString;

void st_uglyf(const char *string, ...);

#endif //PARCOURS_COMMON_H
