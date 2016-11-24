#ifndef PARCOURS_PAIRWISE_ALIGNER_H
#define PARCOURS_PAIRWISE_ALIGNER_H

#include "dpMatrix.h"
#include "stateMachine.h"

#define PAIR_ALIGNMENT_PROB_1 10000000

typedef struct _alignment_parameters {
    double threshold;
    int64_t expansion;
    bool ignore_gaps;
} AlignmentParameters;

// SetType is the 'kind' of alignment (eg. nucleotide), sn is the state number (a 
// way to specify the kind of HMM) this is essentially equivalent to <set_size, sn>
template<class Hmm, size_t sn>
class PairwiseAlignment {
public:
    PairwiseAlignment<Hmm, sn>(Hmm& hmm, SymbolString&, SymbolString&, AnchorPairs&, AlignmentParameters, 
                               bool ragged_end=false);
    
    void Align();

    AlignedPairs& AlignedPairsGetter();
    
    // Gives the average posterior match probability per base of the two sequences, 
    // treating bases in indels as having 0 match probability, or optionally ignores
    // the gaps
    double Score(bool ignore_gaps);

    double TotalProbability();

    static void PosteriorMatchProbabilities(int64_t xay, double total_probability, 
                                            double threshold, HiddenState match_state,
                                            DpMatrix<double, sn>& forward_mat, 
                                            DpMatrix<double, sn>& backward_mat, 
                                            AlignedPairs& aligned_pairs);
private:
    Hmm& model;
    SymbolString& sX;
    SymbolString& sY;

    AlignmentParameters params;
    AlignedPairs aligned_pairs;
    double total_probability = LOG_ZERO;
    DpMatrix<double, sn> forward_matrix;
    DpMatrix<double, sn> backward_matrix;
    void forward_algorithm(Hmm&, SymbolString&, SymbolString&);
    void backward_algorithm(Hmm&, SymbolString&, SymbolString&);
    void posterior_match_probs();
    void calculate_total_prob(Hmm&);
    bool aligned = false;
    
};

//typedef PairwiseAlignment<FiveStateSymbolHmm, fiveState> FiveStateSymbolAlignment;

//template class PairwiseAlignment<FiveStateSymbolHmm, fiveState>;


#endif // PARCOURS_PAIRWISE_ALIGNER_H
