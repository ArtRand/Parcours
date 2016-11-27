// 
// banded_dp_matrix.h
//

#ifndef PARCOURS_DIAGONAL_H
#define PARCOURS_DIAGONAL_H

#include "stl_includes.h"
#include "parcours_exceptions.h"

// Diagonal
class Diagonal {
public:
    Diagonal(int64_t xay, int64_t xmyL, int64_t xmyR);

    ~Diagonal() {};

    int64_t Xay();

    int64_t MinXmy();

    int64_t MaxXmy();

    int64_t Width();

    void Set(int64_t xay, int64_t xmyL, int64_t xmyR);

    std::string ToString();
    
    static int64_t XCoordinate(int64_t, int64_t);

    static int64_t YCoordinate(int64_t, int64_t);

    bool operator == (Diagonal& other) const;

private:
    int64_t xay;
    int64_t xmyL;
    int64_t xmyR;
};

#endif // PARCOURS_DIAGONAL_H
