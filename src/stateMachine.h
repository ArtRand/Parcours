
#ifndef PARCOURS_STATEMACHINE_H
#define PARCOURS_STATEMACHINE_H

#include "stl_includes.h"
#include "logAdd.h"

typedef enum _hidden_state {
    match = 0,
    shortGapX = 1,
    shortGapY = 2,
    longGapX = 3,
    longGapY = 4,
} HiddenState;

typedef enum _set_types { 
    nucleotide = 4,
} SetType;

template<size_t set_size, size_t state_number>
class StateMachine {
protected:
    StateMachine<set_size, state_number>();

    virtual double StartStateProb(HiddenState state, bool ragged_end) = 0;
    
    virtual double EndStateProb(HiddenState state, bool ragged_end) = 0;
    
    // TODO fill this in.. unless it's always overridden 
    //virtual void CellCalculate() = 0;

    virtual int64_t StateNumber() = 0;
    
    std::array<double, set_size> match_probs;

    std::array<double, set_size> x_gap_probs;

    std::array<double, set_size> y_gap_probs;

    int64_t _state_number = state_number;
};

template<size_t set_size>
class StateMachine5 : public StateMachine<set_size, 5> {
public:
    StateMachine5<set_size>();

    double StartStateProb(HiddenState state, bool ragged_end);
    
    double EndStateProb(HiddenState state, bool ragged_end);

    std::function<double(HiddenState state, bool ragged_end)> EndStateProbFcn();
    
    std::function<double(HiddenState state, bool ragged_end)> StartStateProbFcn();

    int64_t StateNumber();

    //void CellCalculate();
    
private:
    // transitions
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_SHORT_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_X; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_X; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_X; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_MATCH_FROM_SHORT_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_Y; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_Y; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_Y; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_Y; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_Y; //0.0073673675173412815f;
};


#endif
