#ifndef PARCOURS_DP_MATRIX_H
#define PARCOURS_DP_MATRIX_H

#include "dpDiagonal.h"

template<class T, size_t sn>
class DpMatrix {
public:
    DpMatrix(int64_t);

    DpDiagonal<T, sn>& DpDiagonalGetter(int64_t xay);

    //void Init(Band<T, sn> band):

    int64_t ActiveDiagonals();

    bool DiagonalCheck(int64_t xay);

    void CreateDpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR);

    void DeleteDpDiagonal(int64_t xay);
    
    int64_t DiagonalNumber();
private:
    std::vector<DpDiagonal<T, sn>> dpDiagonals;
    int64_t diagonal_number;
    int64_t active_diagonals;
    int64_t _state_number = sn;
};

#endif // PARCOURS_DP_MATRIX_H
