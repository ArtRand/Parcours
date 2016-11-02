#include "stateMachine.h"


// StateMachine base class
template<size_t set_size, size_t state_number>
StateMachine<set_size, state_number>::StateMachine() {  }

// StateMachine5 
template<size_t set_size>
StateMachine5<set_size>::StateMachine5() {  }

template<size_t set_size>
int64_t StateMachine5<set_size>::StateNumberGetter() {
    return StateMachine<set_size, 5>::StateNumber;
}

template<size_t set_size>
double StateMachine5<set_size>::StartStateProb(HiddenState state, bool ragged_end) {
    if (ragged_end) {
        return (state == longGapX || state == longGapY) ? 0 : LOG_ZERO;
    } else {
        return state == match ? 0 : LOG_ZERO;
    }
}

template class StateMachine5<4>;
template class StateMachine<4, 5>;
