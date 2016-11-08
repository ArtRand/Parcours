#include "stateMachine.h"


// StateMachine base class
template<size_t set_size, size_t state_number>
StateMachine<set_size, state_number>::StateMachine() {  }

// StateMachine5 
template<size_t set_size>
StateMachine5<set_size>::StateMachine5() { 
    TRANSITION_MATCH_CONTINUE =         -0.030064059121770816; //0.9703833696510062f
    TRANSITION_MATCH_FROM_SHORT_GAP_X = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    TRANSITION_MATCH_FROM_LONG_GAP_X =  -5.673280173170473; //1.0 - gapExtend = 0.00343657420938
    TRANSITION_GAP_SHORT_OPEN_X =       -4.34381910900448; //0.0129868352330243
    TRANSITION_GAP_SHORT_EXTEND_X =     -0.3388262689231553; //0.7126062401851738f;
    TRANSITION_GAP_SHORT_SWITCH_TO_X =  -4.910694825551255; //0.0073673675173412815f;
    TRANSITION_GAP_LONG_OPEN_X =        -6.30810595366929; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    TRANSITION_GAP_LONG_EXTEND_X =      -0.003442492794189331; //0.99656342579062f;
    TRANSITION_GAP_LONG_SWITCH_TO_X =   -6.30810595366929; //0.99656342579062f;
    // make it symmetric
    TRANSITION_MATCH_FROM_SHORT_GAP_Y = TRANSITION_MATCH_FROM_SHORT_GAP_X;
    TRANSITION_MATCH_FROM_LONG_GAP_Y =  TRANSITION_MATCH_FROM_LONG_GAP_X;
    TRANSITION_GAP_SHORT_OPEN_Y =       TRANSITION_GAP_SHORT_OPEN_X;
    TRANSITION_GAP_SHORT_EXTEND_Y =     TRANSITION_GAP_SHORT_EXTEND_X;
    TRANSITION_GAP_SHORT_SWITCH_TO_Y =  TRANSITION_GAP_SHORT_SWITCH_TO_X;
    TRANSITION_GAP_LONG_OPEN_Y =        TRANSITION_GAP_LONG_OPEN_X;
    TRANSITION_GAP_LONG_EXTEND_Y =      TRANSITION_GAP_LONG_EXTEND_X;
    TRANSITION_GAP_LONG_SWITCH_TO_Y =   TRANSITION_GAP_LONG_SWITCH_TO_X;
}

template<size_t set_size>
int64_t StateMachine5<set_size>::StateNumber() {
    return StateMachine<set_size, 5>::_state_number;
}

template<size_t set_size>
double StateMachine5<set_size>::StartStateProb(HiddenState state, bool ragged_end) {
    if (ragged_end) {
        return (state == longGapX || state == longGapY) ? 0 : LOG_ZERO;
    } else {
        return state == match ? 0 : LOG_ZERO;
    }
}

template<size_t set_size>
double StateMachine5<set_size>::EndStateProb(HiddenState state, bool ragged_end) {
    auto ragged_end_prob = [&] () -> double {
        switch (state) {
        case match:
            return TRANSITION_GAP_LONG_OPEN_X;
        case shortGapX:
            return TRANSITION_GAP_LONG_OPEN_X;
        case shortGapY:
            return TRANSITION_GAP_LONG_OPEN_Y;
        case longGapX:
            return TRANSITION_GAP_LONG_EXTEND_X;
        case longGapY:
            return TRANSITION_GAP_LONG_EXTEND_Y;
        default:
            return LOG_ZERO;
        } 
    };
 
    auto end_prob = [&] () -> double {
        switch (state) {
        case match:
            return TRANSITION_MATCH_CONTINUE;
        case shortGapX:
            return TRANSITION_MATCH_FROM_SHORT_GAP_X;
        case shortGapY:
            return TRANSITION_MATCH_FROM_SHORT_GAP_Y;
        case longGapX:
            return TRANSITION_MATCH_FROM_LONG_GAP_X;
        case longGapY:
            return TRANSITION_MATCH_FROM_LONG_GAP_Y;
        default:
            return LOG_ZERO;
        }
    };

    if (ragged_end) { 
        return ragged_end_prob();
    } else {
        return end_prob();
    }
}

template<size_t set_size>
std::function<double(HiddenState s, bool re)> StateMachine5<set_size>::EndStateProbFcn() {
    std::function<double(HiddenState, bool)> lambda = [&] (HiddenState s, bool re) -> double {
            return EndStateProb(s, re);
    };
    return lambda;
}

template<size_t set_size>
std::function<double(HiddenState s, bool re)> StateMachine5<set_size>::StartStateProbFcn() {
    std::function<double(HiddenState, bool)> lambda = [&] (HiddenState s, bool re) -> double {
            return StartStateProb(s, re);
    };
    return lambda;
}

template class StateMachine5<4>;
template class StateMachine<4, 5>;
