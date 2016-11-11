#include "dpDiagonal.h"

template<class T, size_t sn>
DpDiagonal<T, sn>::DpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR): diagonal(xay, xmyL, xmyR) {
    assert(diagonal.Width() >= 0);
    cells.resize(sn * diagonal.Width(), LOG_ZERO);
    active = true;
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
T DpDiagonal<T, sn>::CellGetVal(int64_t xmy, HiddenState s) {
    if (xmy < diagonal.MinXmy() || xmy > diagonal.MaxXmy() || !active) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    if ((diagonal.Xay() + xmy) % 2 != 0) {
        throw ParcoursException("[DpDiagonal::CellGetVal] Illegal request %" PRIi64 "\n", xmy);
    }
    //return cells.at((((xmy - diagonal.MinXmy()) / 2) * _state_number) + s);
    return cells.at(state_index(xmy, s));
}

template<class T, size_t sn>
T *DpDiagonal<T, sn>::CellGetter(int64_t xmy) {
    if (xmy < diagonal.MinXmy() || xmy > diagonal.MaxXmy() || !active) return nullptr;
    if ((diagonal.Xay() + xmy) % 2 != 0) throw ParcoursException(
            "[DpDiagonal::CellGetVal] Illegal request %" PRIi64 "\n", xmy);
    // return 'match' index so that we're at the start of the cell
    return &cells.at(state_index(xmy, match));
}
template<class T, size_t sn>
void DpDiagonal<T, sn>::CellSetter(int64_t xmy, HiddenState s, T value) {
    if (!active) throw ParcoursException("[DpDiagonal::CellSetter] DpDiagonal not active\n");
    cells.at(state_index(xmy, s)) = value;
}

template<class T, size_t sn>
bool DpDiagonal<T, sn>::CellCheck(int64_t xmy) {
    if (xmy < diagonal.MinXmy() || xmy > diagonal.MaxXmy()) return false;
    return true;
}

template<class T, size_t sn> 
void DpDiagonal<T, sn>::InitValues(std::function<double(HiddenState s, bool re)> StateValueGetter) {
    for (int64_t i = diagonal.MinXmy(); i <= diagonal.MaxXmy(); i +=2) {
        assert(CellCheck(i));
        for (int64_t s = 0; s < _state_number; s++) {
            CellSetter(i, static_cast<HiddenState>(s), StateValueGetter(static_cast<HiddenState>(s), false));
        }
    }
    active = true;
}

template<class T, size_t sn>
T DpDiagonal<T, sn>::Dot(DpDiagonal& d2) {
    if (!active) throw ParcoursException("[DpDiagonal::Dot(DpDiagonal)] Diagonal not active\n");
    if (!d2.IsActive()) throw ParcoursException("[DpDiagonal::Dot(DpDiagonal)] other diagonal not active\n");
    double total_prob = LOG_ZERO;
    int64_t xmy = diagonal.MinXmy();
    while (xmy <= diagonal.MaxXmy()) {
        double p = CellGetVal(xmy, match) + d2.CellGetVal(xmy, match);
        for (int64_t s = 1; s < _state_number; s++) {
            p = logAdd(p, (CellGetVal(xmy, static_cast<HiddenState>(s)) 
                           + d2.CellGetVal(xmy, static_cast<HiddenState>(s))));
        }
        total_prob = logAdd(total_prob, p);
        xmy += 2;
    }
    return total_prob;
}
/*
template<class T, size_t sn>
T DpDiagonal<T, sn>::Dot(int64_t xmy, std::function<double(HiddenState s, bool re)> StateValueGetter) {
    if (!active) throw ParcoursException("[DpDiagonal::Dot(DpDiagonal)] Diagonal not active\n");
    double 
}
*/
template<class T, size_t sn>
int64_t DpDiagonal<T, sn>::StateNumber() { return _state_number; }

template<class T, size_t sn>
bool DpDiagonal<T, sn>::IsActive() { return active; }

template<class T, size_t sn>
void DpDiagonal<T, sn>::Deactivate() { 
    if (!active) throw ParcoursException("[DpDiagonal::Deactivate] Cannot deactivate deactive diagonal\n");
    cells.clear();
    active = false; 
}

template<class T, size_t sn>
void DpDiagonal<T, sn>::Activate() { 
    if (active) throw ParcoursException("[DpDiagonal::Activate] Cannot activate active diagonal\n");
    active = true; 
}

template<class T, size_t sn>
int64_t DpDiagonal<T, sn>::state_index(int64_t xmy, HiddenState s) {
    if (s >= _state_number) throw ParcoursException("[DpDiagonal::state_index] illegal Hidden state %i", s);
    return (((xmy - diagonal.MinXmy()) / 2) * _state_number) + s;
}

template class DpDiagonal<double, 5>;
