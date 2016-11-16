#ifndef PARCOURS_STATEMACHINE_H
#define PARCOURS_STATEMACHINE_H

#include "stl_includes.h"
#include "common.h"
#include "dpMatrix.h"
#include "dpDiagonal.h"
#include "logAdd.h"
#include "symbol_string.h"

// TransitionFunction, is used to designate forward or backward calculations 
typedef std::function<void(double *, double *, HiddenState, HiddenState, double, double)> TransitionFunction;

// EmissionsInitFunction, a wrapper around functions used to initialize the 
// emissions/transitions matrices of the state machine
template<size_t set_size>
using EmissionsInitFunction = std::function<void(std::array<double, set_size * set_size>& matchprobs,
                                                 std::array<double, set_size>& xgapprobs,
                                                 std::array<double, set_size>& ygapprobs)>;

// StateMachine interface for all other StateMachine models. The interface contains a 
// lot of the documentation for the child classes, a StateMachine is an abstraction on an HMM
template<size_t set_size, size_t state_number>
class StateMachine {
protected:
    StateMachine<set_size, state_number>();

    virtual double StartStateProb(HiddenState state, bool ragged_end) = 0;
    
    virtual double EndStateProb(HiddenState state, bool ragged_end) = 0;
    
    // Gap/Match probability functions, return the probability of matching 
    // an element of a sequence to a gap or to another element. These 
    // functions can be overridden to use other datatypes, e.g. continuous 
    // data or amino acids.
    virtual double GapXProb(Symbol cX) = 0;

    virtual double GapYProb(Symbol cY) = 0;

    virtual double MatchProb(Symbol cX, Symbol cY) = 0;

    // TODO fill this in.. unless it's always overridden 
    //virtual void CellCalculate() = 0;

    // TODO also add DpDiagonalCalculation interfaces?

    virtual int64_t StateNumber() const = 0;

    virtual std::array<HiddenState, state_number> States() const = 0;

    // initializes emissions for matches and gaps to defaults based on the 
    // type of alignment
    virtual void InitializeEmissions(EmissionsInitFunction<set_size>) = 0;

    std::array<double, set_size * set_size > match_probs;

    std::array<double, set_size> x_gap_probs;

    std::array<double, set_size> y_gap_probs;

    std::array<HiddenState, state_number> _states;

    int64_t _state_number = state_number;

    int64_t _set_size = set_size;
};

template<size_t set_size>
class StateMachine5 : public StateMachine<set_size, fiveState> {
public:
    StateMachine5<set_size>();

    double StartStateProb(HiddenState state, bool ragged_end);
    
    double EndStateProb(HiddenState state, bool ragged_end);

    double GapXProb(Symbol cX);

    double GapYProb(Symbol cY);

    double MatchProb(Symbol cX, Symbol cY);

    // returns a function wrappter around EndStateProb/StartStateProb
    std::function<double(HiddenState state, bool ragged_end)> EndStateProbFcn();
    
    std::function<double(HiddenState state, bool ragged_end)> StartStateProbFcn();

    int64_t StateNumber() const;

    std::array<HiddenState, fiveState> States() const;

    void InitializeEmissions(EmissionsInitFunction<set_size> initFunc);

    // Basic dynamic programming function between cells of the DP matrix
    // depending on `do_transition` the function either performs a forward 
    // or backward calculation, this is taken advantage of in the 
    // following two functions.
    void CellCalculate(double *current, double *lower, double *middle, double *upper,
                       const Symbol& cX, const Symbol& cY,
                       TransitionFunction do_transition);

    // Performs dynamic programming on `diagonals`, in other words performs the
    // DP on an entire anti-diagonal of the DP matrix at a given `xay` (for the
    // overloaded version). The `do_transition` parameter is given to `CellCalculate`
    void DpDiagonalCalculation(DpDiagonal<double, fiveState> *curr, 
                               DpDiagonal<double, fiveState> *m1, 
                               DpDiagonal<double, fiveState> *m2, 
                               const SymbolString& cX, const SymbolString& cY,
                               TransitionFunction do_transition);
    
    void DpDiagonalCalculation(int64_t xay, DpMatrix<double, fiveState>& mat, 
                               const SymbolString& sX, const SymbolString& sY, 
                               TransitionFunction do_transition);

    StateMachineType type;
    
private:
    // transitions for 5-state HMM
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

struct DoTransitionForward {
    void operator () (double *from_cells, double *to_cells, 
                      HiddenState from, HiddenState to, 
                      double eP, double tP);
};

struct DoTransitionBackward {  
    void operator () (double *from_cells, double *to_cells, 
                      HiddenState from, HiddenState to, 
                      double eP, double tP);
};

typedef StateMachine5<nucleotide> FiveStateSymbolHmm;

// Sets the state machine (HMM) emissions tables (GAPS and MATCHES) 
// to sensible defaults. 
EmissionsInitFunction<nucleotide> SetNucleotideEmissionsToDefauts();

#endif
