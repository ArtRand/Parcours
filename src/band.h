// band.h

#ifndef PARCOURS_BAND_H
#define PARCOURS_BAND_H

#include "diagonal.h"
#include "stl_includes.h"

class Band {
public:
    Band(std::vector<std::pair<int64_t, int64_t>> anchors, int64_t lX, int64_t lY, int64_t expansion);

    ~Band() {};

    Band(Band& other) = default;

    Band& operator = (Band& other) = default;

    Diagonal Next();

    Diagonal Previous();

private:
    // used a vector instead of a list here bc. random access is faster
    std::vector<Diagonal> diagonals;  
    int64_t maxLxLy;  // lXalY
    int64_t index;
};

#endif
