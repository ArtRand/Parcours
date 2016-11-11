#include "dpMatrix.h"

template<class T, size_t sn>
DpMatrix<T, sn>::DpMatrix(int64_t dn): diagonal_number(dn), active_diagonals(0) {}

template<class T, size_t sn>
DpDiagonal<T, sn> *DpMatrix<T, sn>::DpDiagonalGetter(int64_t xay) {
    if (!DiagonalCheck(xay)) {
        //throw ParcoursException("[DpMatrix::DpDiagonalGetter] No Diagonal at %" PRIi64 "\n", xay);
        return nullptr;
    }
    return &dpDiagonals.at(xay);
}

template<class T, size_t sn>
bool DpMatrix<T, sn>::DiagonalCheck(int64_t xay) {
    if (xay < 0 || xay > diagonal_number || xay >= dpDiagonals.size()) {
        return false;
    } 
    return true;
}

template<class T, size_t sn>
void DpMatrix<T, sn>::CreateDpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR) {
    if (xay < 0) throw ParcoursException(
            "[DpMatrix::CreateDpDiagonal] Invalid xay > 0\n");
    if (xay > diagonal_number) throw ParcoursException(
            "[DpMatrix::CreateDpDiagonal] Invalid xay > diagonal_number");
    if (xay > dpDiagonals.size()) throw ParcoursException(
            "[DpMatrix::CreateDpDiagonal] Invalid xay > len(dpDiagonals)");
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
int64_t DpMatrix<T, sn>::ActiveDiagonals() { return active_diagonals; }

template<class T, size_t sn>
int64_t DpMatrix<T, sn>::DiagonalNumber() { return diagonal_number; }

template class DpMatrix<double, 5>;
