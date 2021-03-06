// DpDiagonal.h

#ifndef PARCOURS_DP_DIAGONAL_H
#define PARCOURS_DP_DIAGONAL_H

#include "diagonal.h"
#include "common.h"
#include "log_add.h"

template<class T, size_t sn>
class DpDiagonal {
public:
    DpDiagonal();

    DpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR);

    bool operator == (DpDiagonal& other) const;
    
    Diagonal DiagonalGetter() const;
    
    std::vector<T>& Cells();

    bool CellCheck(int64_t xmy);

    T CellGetVal(int64_t xmy, HiddenState s);

    T *CellGetter(int64_t xmy);

    void CellSetter(int64_t xmy, HiddenState s, T value);
    
    void InitValues(std::function<double(HiddenState s, bool re)> StateValueGetter, bool ragged_end=false);
    
    T Dot(DpDiagonal& d2);

    int64_t StateNumber();

    bool IsActive();

    void Deactivate();

    void Activate(int64_t xay, int64_t xmyL, int64_t xmyR);

private:
    int64_t stateIndex(int64_t xmy, HiddenState s);
    Diagonal diagonal;
    std::vector<T> cells;
    int64_t _state_number = sn;
    bool active;
};

#endif // PARCOURS_DP_DIAGONAL_H
