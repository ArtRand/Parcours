#ifndef PARCOURS_SEQUENCE_H
#define PARCOURS_SEQUENCE_H

#include "stl_includes.h"

class Sequence {
public:
    std::string seq;
    std::string label;

    Sequence() {};

    Sequence(std::string seq, std::string lab);

};

#endif
