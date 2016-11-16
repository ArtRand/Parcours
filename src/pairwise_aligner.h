#ifndef PARCOURS_PAIRWISE_ALIGNER_H
#define PARCOURS_PAIRWISE_ALIGNER_H

#include "dpMatrix.h"
#include "stateMachine.h"

#define PAIR_ALIGNMENT_PROB_1 10000000

typedef struct _alignment_parameters {
    double threshold;
    int64_t expansion;
} AlignmentParameters;

// SetType is the 'kind' of alignment (eg. nucleotide), sn is the state number (a 
// way to specify the kind of HMM)
template<class Hmm, size_t sn>
class PairwiseAlignment {
public:
    PairwiseAlignment<Hmm, sn>(Hmm& hmm, SymbolString&, SymbolString&, AnchorPairs&, AlignmentParameters);
    
    void Align(Hmm& hmm, SymbolString& sX, SymbolString& sY);

    AlignedPairs& AlignedPairsGetter();
    
    double Score();

    double TotalProbability();

    static void PosteriorMatchProbabilities(int64_t xay, double total_probability, 
                                            double threshold, HiddenState match_state,
                                            DpMatrix<double, sn>& forward_mat, 
                                            DpMatrix<double, sn>& backward_mat, 
                                            AlignedPairs& aligned_pairs);
private:
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

// generates aligned pairs, which are tuples (prob, x, y) where prob is the posterior
// probability that x is aligned to y: p(x<>y | x, y) = p(x, y, x<>y) / p(x, y)
// populates `aligned_pairs`, does not check for `forward_mat` or `backward_mat` containing
// actual alignment probabilities. All probabilities are multiplied by `PAIR_ALIGNMENT_PROB_1`
//template<class T, size_t sn>
//void PosteriorMatchProbabilities(int64_t xay, double total_probability, 
//                                 double threshold, HiddenState match_state,
//                                 DpMatrix<T, sn>& forward_mat, DpMatrix<T, sn>& backward_mat, 
//                                 AlignedPairs& aligned_pairs);

#endif // PARCOURS_PAIRWISE_ALIGNER_H
