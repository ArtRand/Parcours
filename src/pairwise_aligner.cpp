#include "pairwise_aligner.h"

template<class Hmm, size_t sn>
PairwiseAlignment<Hmm, sn>::PairwiseAlignment(Hmm& hmm,SymbolString& sx, SymbolString& sy,
                                              AnchorPairs& anchors, AlignmentParameters p, bool ragged_end): 
                                                  model(hmm), sX(sx), sY(sy), params(p),
                                                  forward_matrix(sx.size(), sy.size(), anchors, p.expansion), 
                                                  backward_matrix(sx.size(), sy.size(), anchors, p.expansion) {
    // check input and initialization of Dp matrices                                                     
    if (static_cast<uint64_t>(forward_matrix.DiagonalNumber()) != (sX.size() + sY.size())) {
        throw ParcoursException("[PairwiseAlignment::PairwiseAlignment] "
            "Forward matrix does not have correct diagonal number");
    } 
    if (static_cast<uint64_t>(backward_matrix.DiagonalNumber()) != (sX.size() + sY.size())) {
        throw ParcoursException("[PairwiseAlignment::PairwiseAlignment] "
                "Backward matrix does not have correct diagonal number");
    }
    if (forward_matrix.DiagonalNumber() != backward_matrix.DiagonalNumber()) { 
        throw ParcoursException("[PairwiseAlignment::PairwiseAlignment] "
                "Dp matrices don't have the same number of diagonals\n");
    }
    // initialize 
    forward_matrix.DpDiagonalGetter(0)->InitValues(model.StartStateProbFcn(), ragged_end);
    backward_matrix.DpDiagonalGetter(backward_matrix.DiagonalNumber())->InitValues(model.EndStateProbFcn(), ragged_end);
}

template<class Hmm, size_t sn>
void PairwiseAlignment<Hmm, sn>::Align() {
    forwardAlgorithm(model, sX, sY);
    backwardAlgorithm(model, sX, sY);
    calculateTotalProbability(model);
    posteriorMatchProbs();
    aligned = true;
}   

template<class Hmm, size_t sn>
AlignedPairs& PairwiseAlignment<Hmm, sn>::AlignedPairsGetter() {
    if (!aligned) Align();
    return aligned_pairs;
}


template<class Hmm, size_t sn>
double PairwiseAlignment<Hmm, sn>::Score(bool ignore_gaps) {
    if (!aligned) Align();
    double total_score = [this] () -> double {
        auto accumulator = [&] (double i, AlignedPair& p) {
            return i + std::get<0>(p);
        };
        double s = std::accumulate(begin(aligned_pairs), end(aligned_pairs), 0.0, accumulator);
        return s;
    }();

    auto score_ignoring_gaps = [this, &total_score] () -> double {
        return 100.0 * total_score / ((double )aligned_pairs.size() * PAIR_ALIGNMENT_PROB_1);
    };

    auto score_with_gaps = [this, &total_score] () -> double {
        int64_t lX = sX.size();
        int64_t lY = sY.size();
        int64_t lBoth = lX + lY;
        return 100.0 * ((lBoth) == 0 ? 0 : (2.0 * total_score) / (lBoth * PAIR_ALIGNMENT_PROB_1)); 
    };

    return ignore_gaps ? score_ignoring_gaps() : score_with_gaps();
}

template<class Hmm, size_t sn>
void PairwiseAlignment<Hmm, sn>::PosteriorMatchProbabilities(int64_t xay, double total_probability, 
                                                           double threshold, HiddenState match_state,
                                                           DpMatrix<double, sn>& forward_mat, 
                                                           DpMatrix<double, sn>& backward_mat, 
                                                           AlignedPairs& aligned_pairs) {
    DpDiagonal<double, sn> *forward_diagonal= forward_mat.DpDiagonalGetter(xay);
    DpDiagonal<double, sn> *backward_diagonal= backward_mat.DpDiagonalGetter(xay);
    int64_t xmy = forward_diagonal->DiagonalGetter().MinXmy();
    while (xmy <= forward_diagonal->DiagonalGetter().MaxXmy()) {
        int64_t x = Diagonal::XCoordinate(forward_diagonal->DiagonalGetter().Xay(), xmy);
        int64_t y = Diagonal::YCoordinate(forward_diagonal->DiagonalGetter().Xay(), xmy);
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

template<class Hmm, size_t sn>
void PairwiseAlignment<Hmm, sn>::forwardAlgorithm(Hmm& hmm, SymbolString& sX, SymbolString& sY) {
    for (int64_t i = 1; i <= forward_matrix.DiagonalNumber(); i++) {
        hmm.DpDiagonalCalculation(i, forward_matrix, sX, sY, DoTransitionForward());
    }
}

template<class Hmm, size_t sn>
void PairwiseAlignment<Hmm, sn>::backwardAlgorithm(Hmm& hmm, SymbolString& sX, SymbolString& sY) {
    for (int64_t i = backward_matrix.DiagonalNumber(); i > 0; i--) {
        hmm.DpDiagonalCalculation(i, backward_matrix, sX, sY, DoTransitionBackward());
    }
}

template<class Hmm, size_t sn>
void PairwiseAlignment<Hmm, sn>::calculateTotalProbability(Hmm& hmm) {
    double total_fw_prob = forward_matrix.TotalProbability(hmm.EndStateProbFcn(), true);
    double total_bw_prob = backward_matrix.TotalProbability(hmm.StartStateProbFcn(), false);
    if (std::abs(total_fw_prob - total_bw_prob) > 0.01) throw ParcoursException(
            "[PairwiseAlignment::calculateTotalProbability] FW and BW total probs too different"
            "FW prob: %f BW prob: %f\n", total_fw_prob, total_bw_prob);
    total_probability = total_fw_prob;
}

template<class Hmm, size_t sn>
void PairwiseAlignment<Hmm, sn>::posteriorMatchProbs() {
    for (int64_t i = 1; i <= forward_matrix.DiagonalNumber(); i++) {
        PosteriorMatchProbabilities(i, total_probability, params.threshold, match, 
                                    forward_matrix, backward_matrix, aligned_pairs);
    }
}

template class PairwiseAlignment<FiveStateSymbolHmm, fiveState>;

//template void PosteriorMatchProbabilities<double, fiveState>(
//    int64_t xay, double total_probability, 
//    double threshold, HiddenState match_state,
//    DpMatrix<double, fiveState>& forward_mat, DpMatrix<double, fiveState>& backward_mat, 
//    AlignedPairs& aligned_pairs);
