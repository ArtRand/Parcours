
#ifndef PARCOURS_STATEMACHINE_H
#define PARCOURS_STATEMACHINE_H

#include "stl_includes.h"
#include "logAdd.h"

typedef enum _hidden_state {
    match = 0,
    insertX = 1,
    insertY = 2,
    longGapX = 3,
    longGapY = 4,
} HiddenState;

template<size_t set_size, size_t state_number>
class StateMachine {
protected:
    //StateMachine<set_size, state_number>() {};
    StateMachine<set_size, state_number>();

    virtual double StartStateProb(HiddenState state, bool ragged_end) = 0;
    
    //virtual double EndStateProb(int64_t state, bool ragged_end) = 0;
    
    // TODO fill this in.. unless it's always overridden 
    //virtual void CellCalculate() = 0;
    
    std::array<double, set_size> match_probs;

    std::array<double, set_size> x_gap_probs;

    std::array<double, set_size> y_gap_probs;

    int64_t StateNumber = state_number;

};

template<size_t set_size>
class StateMachine5 : public StateMachine<set_size, 5> {
public:
    //StateMachine5<set_size>() {};
    StateMachine5<set_size>();

    double StartStateProb(HiddenState state, bool ragged_end);

    //double EndStateProb(HiddenState state, bool ragged_end);

    //void CellCalculate();

    //int64_t StateNumberGetter() { return StateMachine<set_size, 5>::StateNumber; };
    int64_t StateNumberGetter();
};

//template class StateMachine<4, 5>;
//template class StateMachine5<4>;

#endif
