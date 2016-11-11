#include "dpMatrix.h"

template<class T, size_t sn>
DpMatrix<T, sn>::DpMatrix(int64_t lx, int64_t ly): diagonal_number(lx + ly), 
                                                   active_diagonals(0), 
                                                   lX(lx), lY(ly) {
    dpDiagonals.reserve(lx + ly);
}

template<class T, size_t sn>
DpDiagonal<T, sn> *DpMatrix<T, sn>::DpDiagonalGetter(int64_t xay) {
    if (!DiagonalCheck(xay)) {
        return nullptr;
    }
    return &dpDiagonals.at(xay);
}

template<class T, size_t sn>
bool DpMatrix<T, sn>::DiagonalCheck(int64_t xay) {
    if (xay < 0 || xay > diagonal_number || xay >= static_cast<int>(dpDiagonals.size())) {
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
    DpDiagonal<T, sn> d(xay, xmyL, xmyR);
    dpDiagonals.insert(begin(dpDiagonals) + xay, d);
    active_diagonals++;
}

template<class T, size_t sn>
void DpMatrix<T, sn>::AddDiagonal(Diagonal d) {
    CreateDpDiagonal(d.Xay(), d.MinXmy(), d.MaxXmy());
}

template<class T, size_t sn>
void DpMatrix<T, sn>::DeleteDpDiagonal(int64_t xay) {
    if (xay < 0) throw ParcoursException(
            "[DpMatrix::DeleteDpDiagonal] Invalid xay > 0\n");
    if (xay > diagonal_number) throw ParcoursException(
            "[DpMatrix::DeleteDpDiagonal] Invalid xay > diagonal_number");
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
            total_prob = logAdd(total_prob, cell[s] + fc(static_cast<HiddenState>(s), false));
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
