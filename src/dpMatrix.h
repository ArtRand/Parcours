#ifndef PARCOURS_DP_MATRIX_H
#define PARCOURS_DP_MATRIX_H

#include "dpDiagonal.h"
#include "band.h"

template<class T, size_t sn>
class DpMatrix {
public:
    DpMatrix<T, sn>(int64_t lx, int64_t ly, AnchorPairs& anchors, int64_t expansion);

    DpMatrix<T, sn>(int64_t lx, int64_t ly);

    DpDiagonal<T, sn> *DpDiagonalGetter(int64_t xay);

    int64_t ActiveDiagonals();

    // checks for illegal input xay, cannot be < 0, > diagonal number or >= len(dpDiagonals),
    // also retruns false if the diagonal at xay is inactive 
    bool DiagonalCheck(int64_t xay);

    void CreateDpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR);

    void CreateDpDiagonal(Diagonal d);

    // "Deletes" a diagonal in the matrix, this means that the cells are cleared
    // the diagonal is set to (-1, -1, -1) and the active flag is set to false, 
    // the size of the dpDiagonals is the same (in fact it should always be the same)
    void DeleteDpDiagonal(int64_t xay);
    
    int64_t DiagonalNumber();

    T TotalProbability(std::function<double(HiddenState s, bool re)> StateValueGetter, bool forward);
private:
    int64_t diagonal_number;
    int64_t active_diagonals;
    int64_t _state_number = sn;
    int64_t lX;
    int64_t lY;

    std::vector<DpDiagonal<T, sn>> dpDiagonals;
    Band<T, sn> band;

};

#endif // PARCOURS_DP_MATRIX_H
