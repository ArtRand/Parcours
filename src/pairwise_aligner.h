#ifndef PARCOURS_PAIRWISE_ALIGNER_H
#define PARCOURS_PAIRWISE_ALIGNER_H

#include "dpMatrix.h"

#define PAIR_ALIGNMENT_PROB_1 10000000

template<class T, size_t sn>
void PosteriorMatchProbabilities(int64_t xay, double total_probability, 
                                 double threshold, HiddenState match_state,
                                 DpMatrix<T, sn>& forward_dp_mat, DpMatrix<T, sn>& backward_dp_mat, 
                                 AlignedPairs& aligned_pairs);

#endif // PARCOURS_PAIRWISE_ALIGNER_H
