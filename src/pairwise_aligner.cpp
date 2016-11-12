#include "pairwise_aligner.h"

template<class T, size_t sn>
void PosteriorMatchProbabilities(int64_t xay, double total_probability, 
                                 double threshold, HiddenState match_state,
                                 DpMatrix<T, sn>& forward_mat, DpMatrix<T, sn>& backward_mat, 
                                 AlignedPairs& aligned_pairs) {
    DpDiagonal<T, sn> *forward_diagonal= forward_mat.DpDiagonalGetter(xay);
    DpDiagonal<T, sn> *backward_diagonal= backward_mat.DpDiagonalGetter(xay);
    int64_t xmy = forward_diagonal->DiagonalGetter().MinXmy();
    while (xmy <= forward_diagonal->DiagonalGetter().MaxXmy()) {
        int64_t x = diagonal_XCoordinate(forward_diagonal->DiagonalGetter().Xay(), xmy);
        int64_t y = diagonal_YCoordinate(forward_diagonal->DiagonalGetter().Xay(), xmy);
        if (x > 0 && y > 0) {
            double *cell_forward = forward_diagonal->CellGetter(xmy);
            double *cell_backward = backward_diagonal->CellGetter(xmy);
            // probabilities are in log-space, then exponentiated 
            double posterior_probability = std::exp(
                (cell_forward[match_state] + cell_backward[match_state]) - 
                total_probability);
            if (posterior_probability >= threshold) {
                if (posterior_probability > 1.0) posterior_probability = 1.0;
                // for numerical stability when computing the total probability of an alignment later
                posterior_probability = std::floor(posterior_probability * PAIR_ALIGNMENT_PROB_1);
                aligned_pairs.emplace_back(posterior_probability, x - 1, y -1);
            } 
        }
        xmy += 2;
    }
}

template void PosteriorMatchProbabilities<double, fiveState>(
    int64_t xay, double total_probability, 
    double threshold, HiddenState match_state,
    DpMatrix<double, fiveState>& forward_mat, DpMatrix<double, fiveState>& backward_mat, 
    AlignedPairs& aligned_pairs);
