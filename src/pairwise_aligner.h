#ifndef PARCOURS_PAIRWISE_ALIGNER_H
#define PARCOURS_PAIRWISE_ALIGNER_H

#include "dpMatrix.h"

#define PAIR_ALIGNMENT_PROB_1 10000000

// generates aligned pairs, which are tuples (prob, x, y) where prob is the posterior
// probability that x is aligned to y: p(x<>y | x, y) = p(x, y, x<>y) / p(x, y)
// populates `aligned_pairs`, does not check for `forward_mat` or `backward_mat` containing
// actual alignment probabilities. All probabilities are multiplied by `PAIR_ALIGNMENT_PROB_1`
template<class T, size_t sn>
void PosteriorMatchProbabilities(int64_t xay, double total_probability, 
                                 double threshold, HiddenState match_state,
                                 DpMatrix<T, sn>& forward_mat, DpMatrix<T, sn>& backward_mat, 
                                 AlignedPairs& aligned_pairs);

#endif // PARCOURS_PAIRWISE_ALIGNER_H
