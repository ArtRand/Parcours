#ifndef PARCOURS_STATEMACHINE_H
#define PARCOURS_STATEMACHINE_H

#include "stl_includes.h"
#include "common.h"
#include "dpMatrix.h"
#include "logAdd.h"

template<size_t set_size, size_t state_number>
class StateMachine {
protected:
    StateMachine<set_size, state_number>();

    virtual double StartStateProb(HiddenState state, bool ragged_end) = 0;
    
    virtual double EndStateProb(HiddenState state, bool ragged_end) = 0;
    
    virtual double GapXProb(Symbol cX) = 0;

    virtual double GapYProb(Symbol cY) = 0;

    virtual double MatchProb(Symbol cX, Symbol cY) = 0;

    // TODO fill this in.. unless it's always overridden 
    //virtual void CellCalculate() = 0;

    virtual const int64_t StateNumber() const = 0;

    //virtual constexpr StateMachineType Type() const = 0;
    
    std::array<double, set_size * set_size > match_probs;

    std::array<double, set_size> x_gap_probs;

    std::array<double, set_size> y_gap_probs;

    int64_t _state_number = state_number;
};

template<size_t set_size>
class StateMachine5 : public StateMachine<set_size, fiveState> {
public:
    StateMachine5<set_size>();

    double StartStateProb(HiddenState state, bool ragged_end);
    
    double EndStateProb(HiddenState state, bool ragged_end);

    virtual double GapXProb(Symbol cX);

    virtual double GapYProb(Symbol cY);

    virtual double MatchProb(Symbol cX, Symbol cY);

    std::function<double(HiddenState state, bool ragged_end)> EndStateProbFcn();
    
    std::function<double(HiddenState state, bool ragged_end)> StartStateProbFcn();

    const int64_t StateNumber() const;
    
    StateMachineType type;
    //constexpr StateMachineType Type() const;
    //void CellCalculate(bool lower, bool middle, bool upper, Symbol& cX, Symbol &cY, 
    //                   DpMatrix<double, fiveState>& mat, 
    //                   std::function<double(double, double&, double, double)> do_transition);
    
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
