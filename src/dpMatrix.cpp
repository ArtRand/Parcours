#include "dpMatrix.h"

template<class T, size_t sn>
DpMatrix<T, sn>::DpMatrix(int64_t lx, int64_t ly, 
                          AnchorPairs& anchors, 
                          int64_t expansion): diagonal_number(lx + ly), 
                                              active_diagonals(0), 
                                              lX(lx), lY(ly),
                                              band(anchors, lx, ly, expansion) {
    dpDiagonals.resize(lx + ly + 1);
    
    for (int64_t i = 0; i <= diagonal_number; i++) {
        CreateDpDiagonal(band.Next());
    }    
}

template<class T, size_t sn>
DpMatrix<T, sn>::DpMatrix(int64_t lx, int64_t ly): diagonal_number(lx + ly), 
                                                   active_diagonals(0), 
                                                   lX(lx), lY(ly), 
                                                   band(lx, ly) {
    dpDiagonals.resize(lx + ly + 1);
}

template<class T, size_t sn>
DpDiagonal<T, sn> *DpMatrix<T, sn>::DpDiagonalGetter(int64_t xay) {
    if (!DiagonalCheck(xay)) return nullptr;
    return &dpDiagonals.at(xay);
}

template<class T, size_t sn>
bool DpMatrix<T, sn>::DiagonalCheck(int64_t xay) {
    if (xay < 0 || 
        xay > diagonal_number || 
        xay >= static_cast<int>(dpDiagonals.size()) || 
        !dpDiagonals.at(xay).IsActive()) {
        return false;
    } 
    return true;
}

template<class T, size_t sn>
void DpMatrix<T, sn>::CreateDpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR) {
    if (xay < 0) throw ParcoursException(
            "[DpMatrix::CreateDpDiagonal] Invalid, xay < 0\n");
    if (xay > diagonal_number) throw ParcoursException(
            "[DpMatrix::CreateDpDiagonal] Invalid, xay > diagonal_number");
    if (xay > static_cast<int>(dpDiagonals.size())) throw ParcoursException(
            "[DpMatrix::CreateDpDiagonal] Invalid, xay > len(dpDiagonals)");
    dpDiagonals.at(xay).Activate(xay, xmyL, xmyR);
    active_diagonals++;
}

template<class T, size_t sn>
void DpMatrix<T, sn>::CreateDpDiagonal(Diagonal d) {
    CreateDpDiagonal(d.Xay(), d.MinXmy(), d.MaxXmy());
}

template<class T, size_t sn>
void DpMatrix<T, sn>::DeleteDpDiagonal(int64_t xay) {
    if (xay < 0) throw ParcoursException(
            "[DpMatrix::DeleteDpDiagonal] Invalid xay > 0\n");
    if (xay > diagonal_number) throw ParcoursException(
            "[DpMatrix::DeleteDpDiagonal] Invalid xay > diagonal_number");
    if (!dpDiagonals.at(xay).IsActive()) throw ParcoursException(
            "[DpMatrix::DeleteDiagonal] Cannot delete deactive diagonal\n");
    dpDiagonals.at(xay).Deactivate();
    active_diagonals--;
    if (active_diagonals < 0) throw ParcoursException(
            "[DpMatrix::DeleteDiagonal] active diagonals cannot become negative\n");
}

template<class T, size_t sn>
T DpMatrix<T, sn>::TotalProbability(std::function<double(HiddenState s, bool re)> StateValueGetter, 
                                    bool forward) {
    auto dot_prod = [this] (double *cell, std::function<double(HiddenState s, bool re)> fc) -> double {
        double total_prob = cell[0] + fc(match, false);
        for (int64_t s = 1; s < _state_number; s++) {
            total_prob = LogAdd(total_prob, cell[s] + fc(static_cast<HiddenState>(s), false));
        }
        return total_prob;
    };
    
    int64_t final_diag = forward ? lX + lY : 0;
    DpDiagonal<T, sn> *d = DpDiagonalGetter(final_diag);
    int64_t final_cell = forward ? lX - lY : 0;
    double *cell = d->CellGetter(final_cell);
    return dot_prod(cell, StateValueGetter);
}

template<class T, size_t sn>
int64_t DpMatrix<T, sn>::ActiveDiagonals() { return active_diagonals; }

template<class T, size_t sn>
int64_t DpMatrix<T, sn>::DiagonalNumber() { return diagonal_number; }

template class DpMatrix<double, 5>;
