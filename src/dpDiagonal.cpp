#include "dpDiagonal.h"

template<class T, size_t sn>
DpDiagonal<T, sn>::DpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR): diagonal(xay, xmyL, xmyR) {
    assert(diagonal.Width() >= 0);
    cells.resize(sn * diagonal.Width(), LOG_ZERO);
}

template<class T, size_t sn> 
Diagonal DpDiagonal<T, sn>::DiagonalGetter() const {
    return diagonal;
}

template<class T, size_t sn>
std::vector<T>& DpDiagonal<T, sn>::Cells() {
    return cells;
}

template<class T, size_t sn>
bool DpDiagonal<T, sn>::operator == (DpDiagonal& other) const {
    Diagonal other_diag = other.DiagonalGetter();
    std::vector<T> c = other.Cells();
    return diagonal == other_diag && cells == c && _state_number == other.StateNumber();
}

template<class T, size_t sn>
T DpDiagonal<T, sn>::CellGetter(int64_t xmy, HiddenState s) {
    if (xmy < diagonal.MinXmy() || xmy > diagonal.MaxXmy()) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    if ((diagonal.Xay() + xmy) % 2 != 0) {
        throw ParcoursException("[DpDiagonal::CellGetter] Illegal request %" PRIi64 "\n", xmy);
    }
    
    //return cells.at((((xmy - diagonal.MinXmy()) / 2) * _state_number) + s);
    return cells.at(state_index(xmy, s));
}

template<class T, size_t sn>
void DpDiagonal<T, sn>::CellSetter(int64_t xmy, HiddenState s, T value) {
    cells.at(state_index(xmy, s)) = value;
}

template<class T, size_t sn>
bool DpDiagonal<T, sn>::CellCheck(int64_t xmy) {
    if (xmy < diagonal.MinXmy() || xmy > diagonal.MaxXmy()) return false;
    return true;
}

template<class T, size_t sn> 
//void DpDiagonal<T, sn>::InitValues(StateMachine5& sm, double (*StateValueGetter)(HiddenState state, bool ragged)) {
void DpDiagonal<T, sn>::InitValues(std::function<double(HiddenState s, bool re)> StateValueGetter) {
    for (int64_t i = diagonal.MinXmy(); i <= diagonal.MaxXmy(); i +=2) {
        assert(CellCheck(i));
        for (int64_t s = 0; s < _state_number; s++) {
            CellSetter(i, static_cast<HiddenState>(s), StateValueGetter(static_cast<HiddenState>(s), false));
        }
    }
}

template<class T, size_t sn>
int64_t DpDiagonal<T, sn>::StateNumber() { return _state_number; }

template<class T, size_t sn>
int64_t DpDiagonal<T, sn>::state_index(int64_t xmy, HiddenState s) {
    return (((xmy - diagonal.MinXmy()) / 2) * _state_number) + s;
}

template class DpDiagonal<double, 5>;
