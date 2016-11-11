#include "stateMachine.h"


// StateMachine base class
template<size_t set_size, size_t state_number>
StateMachine<set_size, state_number>::StateMachine() {  }

// StateMachine5 
template<size_t set_size>
StateMachine5<set_size>::StateMachine5() { 
    TRANSITION_MATCH_CONTINUE          =-0.030064059121770816; //0.9703833696510062f
    TRANSITION_MATCH_FROM_SHORT_GAP_X  =-1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    TRANSITION_MATCH_FROM_LONG_GAP_X   =-5.673280173170473; //1.0 - gapExtend = 0.00343657420938
    TRANSITION_GAP_SHORT_OPEN_X        =-4.34381910900448; //0.0129868352330243
    TRANSITION_GAP_SHORT_EXTEND_X      =-0.3388262689231553; //0.7126062401851738f;
    TRANSITION_GAP_SHORT_SWITCH_TO_X   =-4.910694825551255; //0.0073673675173412815f;
    TRANSITION_GAP_LONG_OPEN_X         =-6.30810595366929; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    TRANSITION_GAP_LONG_EXTEND_X       =-0.003442492794189331; //0.99656342579062f;
    TRANSITION_GAP_LONG_SWITCH_TO_X    =-6.30810595366929; //0.99656342579062f;
    // make it symmetric
    TRANSITION_MATCH_FROM_SHORT_GAP_Y  =TRANSITION_MATCH_FROM_SHORT_GAP_X;
    TRANSITION_MATCH_FROM_LONG_GAP_Y   =TRANSITION_MATCH_FROM_LONG_GAP_X;
    TRANSITION_GAP_SHORT_OPEN_Y        =TRANSITION_GAP_SHORT_OPEN_X;
    TRANSITION_GAP_SHORT_EXTEND_Y      =TRANSITION_GAP_SHORT_EXTEND_X;
    TRANSITION_GAP_SHORT_SWITCH_TO_Y   =TRANSITION_GAP_SHORT_SWITCH_TO_X;
    TRANSITION_GAP_LONG_OPEN_Y         =TRANSITION_GAP_LONG_OPEN_X;
    TRANSITION_GAP_LONG_EXTEND_Y       =TRANSITION_GAP_LONG_EXTEND_X;
    TRANSITION_GAP_LONG_SWITCH_TO_Y    =TRANSITION_GAP_LONG_SWITCH_TO_X;

    type = fiveState;
}

template<size_t set_size>
const int64_t StateMachine5<set_size>::StateNumber() const {
    //return StateMachine<set_size, 5>::_state_number;
    return 5;
}

template<size_t set_size>
void StateMachine5<set_size>::InitializeEmissions(EmissionsInitFunction<set_size> initFunc) {
    initFunc(StateMachine<set_size, fiveState>::match_probs, 
             StateMachine<set_size, fiveState>::x_gap_probs, 
             StateMachine<set_size, fiveState>::y_gap_probs);
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

template<size_t set_size>
double StateMachine5<set_size>::GapXProb(Symbol cX) {
    if (cX > set_size) throw ParcoursException("[StateMachine5::GapXProb] Illegal symbol %i", cX);
    if (cX == n) return -1.386294361;  // log(0.25)
    return StateMachine<set_size, fiveState>::x_gap_probs.at(cX);
}

template<size_t set_size>
double StateMachine5<set_size>::GapYProb(Symbol cY) {
    if (cY > set_size) throw ParcoursException("[StateMachine5::GapYProb] Illegal symbol %i", cY);
    if (cY == n) return -1.386294361;  // log(0.25)
    return StateMachine<set_size, fiveState>::y_gap_probs.at(cY);
}

template<size_t set_size>
double StateMachine5<set_size>::MatchProb(Symbol cX, Symbol cY) {
    if (cX > set_size) throw ParcoursException("[StateMachine5::MatchProb] Illegal cX symbol %i", cX);
    if (cY > set_size) throw ParcoursException("[StateMachine5::MatchProb] Illegal cY symbol %i", cY);
    if (cX == n || cY == n) return -2.772588722; // log(0.25**2)
    return StateMachine<set_size, fiveState>::match_probs.at(cX * set_size + cY);
}

template<size_t set_size>
void StateMachine5<set_size>::CellCalculate(double *current, double *lower, double *middle, double *upper,
                                            const Symbol& cX, const Symbol& cY,
                                            TransitionFunction do_transition) {
    if (lower != nullptr) {
        double eP = GapXProb(cX);
        do_transition(lower, current, match, shortGapX, eP, TRANSITION_GAP_SHORT_OPEN_X);
        do_transition(lower, current, shortGapX, shortGapX, eP, TRANSITION_GAP_SHORT_EXTEND_X);
        do_transition(lower, current, match, longGapX, eP, TRANSITION_GAP_LONG_OPEN_X);
        do_transition(lower, current, longGapX, longGapX, eP, TRANSITION_GAP_LONG_EXTEND_X);
    }

    if (middle != nullptr) {
        double eP = MatchProb(cX, cY);
        do_transition(middle, current, match, match, eP, TRANSITION_MATCH_CONTINUE);
        do_transition(middle, current, shortGapX, match, eP, TRANSITION_MATCH_FROM_SHORT_GAP_X);
        do_transition(middle, current, shortGapY, match, eP, TRANSITION_MATCH_FROM_SHORT_GAP_Y);
        do_transition(middle, current, longGapX, match, eP, TRANSITION_MATCH_FROM_LONG_GAP_X);
        do_transition(middle, current, longGapY, match, eP, TRANSITION_MATCH_FROM_LONG_GAP_Y);
    }

    if (upper != nullptr) {
        double eP = GapYProb(cY);
        do_transition(upper, current, match, shortGapY, eP, TRANSITION_GAP_SHORT_OPEN_Y);
        do_transition(upper, current, shortGapY, shortGapY, eP, TRANSITION_GAP_SHORT_EXTEND_Y);
        do_transition(upper, current, match, longGapY, eP, TRANSITION_GAP_LONG_OPEN_Y);
        do_transition(upper, current, longGapY, longGapY, eP, TRANSITION_GAP_LONG_EXTEND_Y);
    }
}

template<size_t set_size>
void StateMachine5<set_size>::DpDiagonalCalculation(DpDiagonal<double, fiveState> *curr, 
                                                    DpDiagonal<double, fiveState> *m1, 
                                                    DpDiagonal<double, fiveState> *m2, 
                                                    const SymbolString& cX, const SymbolString& cY,
                                                    TransitionFunction do_transition) {
    int64_t xmy = curr->DiagonalGetter().MinXmy();
    
    auto get_position = [] (const SymbolString& S, int64_t XaY, int64_t XmY, 
                              int64_t (*coord_func)(int64_t, int64_t)) -> int64_t {
        int64_t x = coord_func(XaY, XmY);
        if (x < 0 || x > S.size()) throw ParcoursException(
                "[StateMachine5::DpDiagonalCalculation] Illegal coordinate %" PRIi64 "\n", x);
        return x;
    };

    while (xmy <= curr->DiagonalGetter().MaxXmy()) {
        // get the indices of the things
        int64_t iX = get_position(cX, curr->DiagonalGetter().Xay(), xmy, 
                                  diagonal_XCoordinate) - 1;  // minus 1 to preserve legacy awesomeness
        int64_t iY = get_position(cY, curr->DiagonalGetter().Xay(), xmy, 
                                  diagonal_YCoordinate) - 1;  // minus 1 to preserve legacy awesomeness
        
        Symbol x = iX >= 0 ? cX.at(iX) : n;
        Symbol y = iY >= 0 ? cY.at(iY) : n;

        // do the calculation
        double *current = curr->CellGetter(xmy);
        double *lower   = m1 == nullptr ? nullptr : m1->CellGetter((xmy - 1));
        double *middle  = m2 == nullptr ? nullptr : m2->CellGetter(xmy);
        double *upper   = m1 == nullptr ? nullptr : m1->CellGetter((xmy + 1));
        CellCalculate(current, lower, middle, upper, x, y, do_transition);
        xmy += 2;
    }
}


template<size_t set_size>
void StateMachine5<set_size>::DpDiagonalCalculation(int64_t xay, DpMatrix<double, fiveState>& mat, 
                                                    const SymbolString& sX, 
                                                    const SymbolString& sY, 
                                                    TransitionFunction do_transition) {
    DpDiagonalCalculation(mat.DpDiagonalGetter(xay), 
                          mat.DpDiagonalGetter(xay - 1), 
                          mat.DpDiagonalGetter(xay - 2), 
                          sX, sY, 
                          do_transition);
}


void DoTransitionForward::operator () (double *from_cells, double *to_cells, 
                                       HiddenState from, HiddenState to, 
                                       double eP, double tP) {
    to_cells[to] = logAdd(to_cells[to], from_cells[from] + (eP + tP));
}

void DoTransitionBackward::operator () (double *from_cells, double *to_cells, 
                                        HiddenState from, HiddenState to, 
                                        double eP, double tP) {
    from_cells[from] = logAdd(from_cells[from], to_cells[to] + (eP + tP));
}

EmissionsInitFunction<nucleotide> SetNucleotideEmissionsToDefauts() {
    std::function<void(std::array<double, nucleotide * nucleotide>&, 
                       std::array<double, nucleotide>&,
                       std::array<double, nucleotide>&)> lambda = 
        [&] (std::array<double, nucleotide * nucleotide>& matchprobs,
        std::array<double, nucleotide>& xgapprobs,
        std::array<double, nucleotide>& ygapprobs) -> void {
            // Set Match probs to default values
            double EMISSION_MATCH=-2.1149196655034745; //log(0.12064298095701059);
            double EMISSION_TRANSVERSION=-4.5691014376830479; //log(0.010367271172731285);
            double EMISSION_TRANSITION=-3.9833860032220842; //log(0.01862247669752685);

            std::array<double, nucleotide * nucleotide> M = { 
                {EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION,
                EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION,
                EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION,
                EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH}
            };
            matchprobs = M;

            // Set Gap probs to default values
            double EMISSION_GAP = -1.6094379124341003; //log(0.2)
            std::array<double, nucleotide> G = { {EMISSION_GAP, EMISSION_GAP, EMISSION_GAP, EMISSION_GAP} };
            xgapprobs = G;
            ygapprobs = G;
        };  
    return lambda;
}

template class StateMachine5<nucleotide>;
template class StateMachine<nucleotide, fiveState>;

