// DpDiagonal.h

#ifndef PARCOURS_DP_DIAGONAL_H
#define PARCOURS_DP_DIAGONAL_H

#include "diagonal.h"
#include "stateMachine.h"
#include "logAdd.h"

template<class T, size_t sn>
class DpDiagonal {
public:
    DpDiagonal(int64_t xay, int64_t xmyL, int64_t xmyR);

    //~DpDiagonal() = default;

    //DpDiagonal(DpDiagonal& other) = default;

    //DpDiagonal& operator = (DpDiagonal& other) = default;
    
    bool operator == (DpDiagonal& other) const;
    
    Diagonal DiagonalGetter() const;
    
    std::vector<T>& Cells();

    bool CellCheck(int64_t xmy);

    T CellGetter(int64_t xmy, HiddenState s);

    void CellSetter(int64_t xmy, HiddenState s, T value);
    
    //void InitValues(StateMachine5<set_size>& sm, double (*StateValueGetter)(HiddenState state, bool ragged));
    void InitValues(std::function<double(HiddenState s, bool re)> StateValueGetter);
    
    int64_t StateNumber();

private:
    int64_t state_index(int64_t xmy, HiddenState s);
    Diagonal diagonal;
    std::vector<T> cells;
    int64_t _state_number = sn;
};

#endif // PARCOURS_DP_DIAGONAL_H
