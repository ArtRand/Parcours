// Pairwise aligner tests
//

#include "log_add.h"
#include "diagonal.h"
#include "band.h"
#include "common.h"
#include "dpMatrix.h"
#include "stateMachine.h"
#include "pairwise_aligner.h"
#include "test_helpers.h"
#include "symbol_string.h"


TEST_CASE("Diagonal Tests", "[PairwiseAligner]") {
    int64_t xL = 10, yL = 20, xU = 30, yU = 0; //Coordinates of the upper and lower
    Diagonal d(xL + yL, xL - yL, xU - yU);
    Diagonal d2(0, 0, 0);
    REQUIRE(d.Xay() == xL + yL);
    REQUIRE(d.MinXmy() == xL - yL);
    REQUIRE(d.MaxXmy() == xU - yU);
    REQUIRE(d.Width() == ((xU - yU - (xL - yL)) / 2 + 1));
    REQUIRE(Diagonal::XCoordinate(xL + yL, xL - yL) == xL);
    REQUIRE(Diagonal::YCoordinate(xL + yL, xL - yL) == yL);
    REQUIRE(d == d);
    REQUIRE(!(d == d2));
    REQUIRE_THROWS_AS(Diagonal d(10, 5, 5), ParcoursException);
    REQUIRE_THROWS_AS(Diagonal d(10, 6, 4), ParcoursException);
}

TEST_CASE("Band Tests", "[PairwiseAligner]") {
    std::vector<std::pair<int64_t, int64_t>> anchors;
    anchors.push_back(std::make_pair(1, 0));
    anchors.push_back(std::make_pair(2, 1));
    anchors.push_back(std::make_pair(3, 3));
    int64_t lX = 6, lY = 5, expansion = 2;
    Band<double, 5> b(anchors, lX, lY, expansion);
    
    // Forward pass
    Diagonal d0(0, 0, 0); REQUIRE(b.Next() == d0);
    Diagonal d1(1, -1, 1); REQUIRE(b.Next() == d1);
    Diagonal d2(2, -2, 2); REQUIRE(b.Next() == d2);
    Diagonal d3(3, -1, 3); REQUIRE(b.Next() == d3);
    Diagonal d4(4, -2, 4); REQUIRE(b.Next() == d4);
    Diagonal d5(5, -1, 3); REQUIRE(b.Next() == d5);
    Diagonal d6(6, -2, 4); REQUIRE(b.Next() == d6);
    Diagonal d7(7, -3, 3); REQUIRE(b.Next() == d7);
    Diagonal d8(8, -2, 2); REQUIRE(b.Next() == d8);
    Diagonal d9(9, -1, 3); REQUIRE(b.Next() == d9);
    Diagonal d10(10, -0, 2); REQUIRE(b.Next() == d10);
    Diagonal d11(11, 1, 1); REQUIRE(b.Next() == d11);
    Diagonal d12(11, 1, 1); REQUIRE(b.Next() == d12);

    // Go backward
    REQUIRE(b.Previous() == d11);
    REQUIRE(b.Previous() == d10);
    REQUIRE(b.Previous() == d9);
    REQUIRE(b.Previous() == d8);
    REQUIRE(b.Previous() == d7);
    
    // Go forward again
    REQUIRE(b.Next() == d7);
    REQUIRE(b.Next() == d8);
    REQUIRE(b.Next() == d9);
    REQUIRE(b.Next() == d10);
    
    //Now carry on back again
    REQUIRE(b.Previous() == d10);
    REQUIRE(b.Previous() == d9);
    REQUIRE(b.Previous() == d8);
    REQUIRE(b.Previous() == d7);
    REQUIRE(b.Previous() == d6);
    REQUIRE(b.Previous() == d5);
    REQUIRE(b.Previous() == d4);
    REQUIRE(b.Previous() == d3);
    REQUIRE(b.Previous() == d2);
    REQUIRE(b.Previous() == d1);
    REQUIRE(b.Previous() == d0);
    
    REQUIRE(b.Next() == d0);
    REQUIRE(b.Next() == d1);

    REQUIRE(b.Previous() == d1);
    REQUIRE(b.Previous() == d0);
}

TEST_CASE("LogAdd Tests", "[NumericTests]") {
    for (int64_t test = 0; test < 100000; test++) {
        double i = RandomDouble();
        double j = RandomDouble();
        double k = i + j;
        double l = exp(LogAdd(log(i), log(j)));
        REQUIRE(l < k + 0.001);
        REQUIRE(l > k - 0.001);
    }
}

TEST_CASE("Test Cell", "[DpTests]") {
    StateMachine5<nucleotide> sM5;
    REQUIRE(sM5.StateNumber() == fiveState);
    sM5.InitializeEmissions(SetNucleotideEmissionsToDefauts());
    double lowerF[sM5.StateNumber()], middleF[sM5.StateNumber()], upperF[sM5.StateNumber()], currentF[sM5.StateNumber()];
    double lowerB[sM5.StateNumber()], middleB[sM5.StateNumber()], upperB[sM5.StateNumber()], currentB[sM5.StateNumber()];
    
    for (int64_t i = 0; i < sM5.StateNumber(); i++) {
        middleF[i] = sM5.StartStateProb(static_cast<HiddenState>(i), false);
        middleB[i] = LOG_ZERO;
        lowerF[i] = LOG_ZERO;
        lowerB[i] = LOG_ZERO;
        upperF[i] = LOG_ZERO;
        upperB[i] = LOG_ZERO;
        currentF[i] = LOG_ZERO;
        currentB[i] = sM5.EndStateProb(static_cast<HiddenState>(i), false);
    }
    
    Symbol x = a;
    Symbol y = a;

    auto dot_prd = [&sM5] (double *cell, std::function<double(HiddenState s, bool re)> fc) -> double {
        double total_prob = cell[0] + fc(match, false);
        for (int64_t s = 1; s < sM5.StateNumber(); s++) {
            total_prob = LogAdd(total_prob, cell[s] + fc(static_cast<HiddenState>(s), false));
        }
        return total_prob;
    };

    // Do forward
    sM5.CellCalculate(lowerF, nullptr, nullptr, middleF, x, y, DoTransitionForward());
    sM5.CellCalculate(upperF, middleF, nullptr, nullptr, x, y, DoTransitionForward());
    sM5.CellCalculate(currentF, lowerF, middleF, upperF, x, y, DoTransitionForward());
    
    // Do backward
    sM5.CellCalculate(currentB, lowerB, middleB, upperB, x, y, DoTransitionBackward());
    sM5.CellCalculate(upperB, middleB, nullptr, nullptr, x, y, DoTransitionBackward());
    sM5.CellCalculate(lowerB, nullptr, nullptr, middleB, x, y, DoTransitionBackward());

    double total_prob_forward = dot_prd(currentF, sM5.EndStateProbFcn());
    double total_prob_backward = dot_prd(middleB, sM5.StartStateProbFcn());
    //st_uglyf("total forward prob : %f total backward prob : %f\n", total_prob_forward, total_prob_backward);
    REQUIRE(std::abs(total_prob_forward - total_prob_backward) < 0.00001);
}

TEST_CASE("Test DpDiagonal", "[DpTests]") {
    SECTION("DpDiagonal is initialized, copied, and has functional getters and setters") {
        StateMachine5<nucleotide> sM5;
        // test equality
        DpDiagonal<double, 5> d(3, -1, 1);
        DpDiagonal<double, 5> d2 = d;
        REQUIRE(d == d2);
        DpDiagonal<double, 5> d3 = d;  // copy constructor
        REQUIRE(d3 == d);
        DpDiagonal<double, 5> d4(0, 0, 0);
        REQUIRE(!(d4 == d));
        d4 = d;  // copy assignment
        REQUIRE(d4 == d);
    
        // test cell getter
        auto check_cell = [&] (int64_t xay) -> void {
            for (int64_t s = 0; s < sM5.StateNumber(); s++) {
                double c1 = d.CellGetVal(xay, static_cast<HiddenState>(s));
                REQUIRE(!std::isnan(c1));
                REQUIRE((d.CellGetter(xay) + s) != NULL);
                double c1_1 = *(d.CellGetter(xay) + s);
                REQUIRE(c1 == c1_1);
            }
        };

        int64_t w, x, y, z;
        w = -1;
        check_cell(w);
        x = 1;
        check_cell(x);
        y = 3;
        z = -3;
        REQUIRE(std::isnan(d.CellGetVal(y, match)));
        REQUIRE(std::isnan(d.CellGetVal(z, match)));
        REQUIRE(d.CellGetter(y) == NULL);
        REQUIRE(d.CellGetter(z) == NULL);
    
        // test set values
        double q = RandomDouble();  
        REQUIRE(d.CellCheck(-1));
        d.CellSetter(-1, match, q);
        REQUIRE(d.CellGetVal(-1, match) == q);
    }
    SECTION("DpDiagonal initialize values and dot product work as expected") {
        // test initialize values
        StateMachine5<nucleotide> sM5;
        DpDiagonal<double, fiveState> d(3, -1, 1);
        d.InitValues(sM5.EndStateProbFcn());;
        double total_prob = LOG_ZERO;
        for (int64_t s = 0; s < sM5.StateNumber(); s++) { // TODO change this to ints (enum)
            double c1 = d.CellGetVal(-1, static_cast<HiddenState>(s));
            double c2 = d.CellGetVal(1, static_cast<HiddenState>(s));
            REQUIRE(c1 == sM5.EndStateProb(static_cast<HiddenState>(s), false));
            REQUIRE(c2 == sM5.EndStateProb(static_cast<HiddenState>(s), false));
            total_prob = LogAdd(total_prob, 2 * c1);
            total_prob = LogAdd(total_prob, 2 * c2);
        }
        DpDiagonal<double, 5> d2 = d;
        double tp = d.Dot(d2);
        REQUIRE(std::abs(total_prob - tp) < 0.001);
    }
}

TEST_CASE("Test DpMatrix", "[DpTests]") {
    SECTION("DpMatrix is constructed with correct defaults") {
        int64_t x = RandomInt(10, 20);
        int64_t y = RandomInt(10, 20);
        DpMatrix<double, fiveState> mat(x, y);
        REQUIRE(mat.ActiveDiagonals() == 0);
        REQUIRE(mat.DiagonalNumber() == x + y);
    }
    SECTION("DpMatrix main tests") {
        int64_t lX = 3;
        int64_t lY = 2;
        DpMatrix<double, fiveState> mat(lX, lY);

        REQUIRE(mat.ActiveDiagonals() == 0);
        
        // check for 'fantom' diagonals
        for (int64_t i = -1; i <= lX + lY + 10; i++) {
            REQUIRE(mat.DpDiagonalGetter(i) == nullptr);
        }
        // make some diagonals in the dpMatrix, and check them, then make sure that
        // the number of active diagonals is correct.
        for (int64_t i = 0; i <= lX + lY; i++) {
            mat.CreateDpDiagonal(i, -i, i);
            DpDiagonal<double, fiveState> d(i, -i, i);
            REQUIRE(*mat.DpDiagonalGetter(i) == d);
            REQUIRE(mat.ActiveDiagonals() == i + 1);
        }
        // delete the diagonals 
        for (int64_t i = lX + lY; i >= 0; i--) {
            mat.DeleteDpDiagonal(i);
            REQUIRE_THROWS_AS(mat.DeleteDpDiagonal(i), ParcoursException);
            REQUIRE(mat.DpDiagonalGetter(i) == nullptr);
            REQUIRE(mat.ActiveDiagonals() == i);
        }
    }
}

TEST_CASE("Test DpDiagonalCalculations", "[DpTests][alignment][current]") {
    SymbolString sX = {{a, g, c, g}};
    SymbolString sY = {{a, g, t, t, c, g}};
    
    int64_t lX = sX.size();
    int64_t lY = sY.size();

    StateMachine5<nucleotide> sM5;
    sM5.InitializeEmissions(SetNucleotideEmissionsToDefauts());
    
    DpMatrix<double, fiveState> forward_mat(lX, lY);
    DpMatrix<double, fiveState> backward_mat(lX, lY);
    
    AnchorPairs anchors;
    anchors.emplace_back(0, 0);
    anchors.emplace_back(1, 1);
    anchors.emplace_back(3, 5);

    int64_t expansion = 2; 
    Band<double, fiveState> band(anchors, lX, lY, expansion);
    
    REQUIRE(forward_mat.DiagonalNumber() == backward_mat.DiagonalNumber());
    REQUIRE(forward_mat.DiagonalNumber() == lX + lY);

    for (int64_t i = 0; i <= forward_mat.DiagonalNumber(); i++) {
        Diagonal d = band.Next();
        forward_mat.CreateDpDiagonal(d);
        backward_mat.CreateDpDiagonal(d);
    }

    forward_mat.DpDiagonalGetter(0)->InitValues(sM5.StartStateProbFcn());
    backward_mat.DpDiagonalGetter(backward_mat.DiagonalNumber())->InitValues(sM5.EndStateProbFcn());

    // forward 
    for (int64_t i = 1; i <= forward_mat.DiagonalNumber(); i++) {
        sM5.DpDiagonalCalculation(i, forward_mat, sX, sY, DoTransitionForward());
    }
    // backward
    for (int64_t i = backward_mat.DiagonalNumber(); i > 0; i--) {
        sM5.DpDiagonalCalculation(i, backward_mat, sX, sY, DoTransitionBackward());
    }
    
    double total_forward_prob = forward_mat.TotalProbability(sM5.EndStateProbFcn(), true);
    double total_backward_prob = backward_mat.TotalProbability(sM5.StartStateProbFcn(), false);
    //st_uglyf("total forward p :%f total backward p: %f\n", total_forward_prob, total_backward_prob);
    REQUIRE(std::abs(total_forward_prob - total_backward_prob) < 0.001);
    
    // Do the total probability calculations manually, then check the TotalProbability method
    // in StateMachine5
    DpDiagonal<double, fiveState> *dpFw  = forward_mat.DpDiagonalGetter(forward_mat.DiagonalNumber());
    DpDiagonal<double, fiveState> *dpBw  = backward_mat.DpDiagonalGetter(backward_mat.DiagonalNumber());
    DpDiagonal<double, fiveState> *dpFw0 = forward_mat.DpDiagonalGetter(0);
    DpDiagonal<double, fiveState> *dpBw0 = backward_mat.DpDiagonalGetter(0);
    
    double test_total_prob_at_diagonal_number = dpFw->Dot(*dpBw);
    double test_total_prob_at_zero            = dpFw0->Dot(*dpBw0);
    REQUIRE(std::abs(test_total_prob_at_diagonal_number - test_total_prob_at_zero) < 0.01);;

    double total_prob_at_diagonal_number = sM5.TotalProbability(forward_mat, backward_mat, false);
    double total_prob_at_zero            = sM5.TotalProbability(forward_mat, backward_mat, true);
    REQUIRE(std::abs(total_prob_at_diagonal_number - total_prob_at_zero) < 0.01);;
    REQUIRE(std::abs(total_prob_at_diagonal_number - test_total_prob_at_diagonal_number) < 0.01);;

    AlignedPairs aligned_pairs;
    REQUIRE(aligned_pairs.size() == 0);
    double threshold = 0.2;
    for (int64_t i = 1; i <= forward_mat.DiagonalNumber(); i++) {
        PairwiseAlignment<FiveStateSymbolHmm, fiveState>::
            PosteriorMatchProbabilities(i, total_forward_prob, 
                                        threshold, match, 
                                        forward_mat, backward_mat, 
                                        aligned_pairs);
    }
    REQUIRE(aligned_pairs.size() == 4);
    CheckAlignedPairs(aligned_pairs, lX, lY);   
    // extracts just the aligned pair coordinates, leaves the posteriors
    auto pairs_no_probs = [&aligned_pairs] () -> AnchorPairs {
        AnchorPairs pairs;
        for (auto p : aligned_pairs) {
            pairs.emplace_back(std::get<1>(p), std::get<2>(p));
        }
        return pairs;
    }();
    // make sure that worked
    REQUIRE(pairs_no_probs.size() == 4);

    AnchorPairs correct_pairs;
    correct_pairs.emplace_back(0, 0);
    correct_pairs.emplace_back(1, 1);
    correct_pairs.emplace_back(2, 4);
    correct_pairs.emplace_back(3, 5);
    REQUIRE(correct_pairs == pairs_no_probs);
}

TEST_CASE("Test PairwiseAlignment", "[alignment]") {
    SECTION("PairwiseAlignment works on simple sequences", "[alignment]") {
        SymbolString sX = {{a, g, c, g}};
        SymbolString sY = {{a, g, t, t, c, g}}; 
        FiveStateSymbolHmm hmm;
        hmm.InitializeEmissions(SetNucleotideEmissionsToDefauts());
        
        AnchorPairs anchors;
        anchors.emplace_back(0, 0);
        anchors.emplace_back(1, 1);
        anchors.emplace_back(3, 5); 
        
        AlignmentParameters p;
        p.expansion = 2;
        p.threshold = 0.2;
        PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln(hmm, sX, sY, anchors, p);
        AlignedPairs aligned_pairs = aln.AlignedPairsGetter();
        REQUIRE(aligned_pairs.size() == 4);
        CheckAlignedPairs(aligned_pairs, sX.size(), sY.size());   
        // extracts just the aligned pair coordinates, leaves the posteriors
        auto pairs_no_probs = [&aligned_pairs] () -> AnchorPairs {
            AnchorPairs pairs;
            for (auto p : aligned_pairs) {
                pairs.emplace_back(std::get<1>(p), std::get<2>(p));
            }
            return pairs;
        }();
        // make sure that worked
        REQUIRE(pairs_no_probs.size() == 4);

        AnchorPairs correct_pairs;
        correct_pairs.emplace_back(0, 0);
        correct_pairs.emplace_back(1, 1);
        correct_pairs.emplace_back(2, 4);
        correct_pairs.emplace_back(3, 5);
        REQUIRE(correct_pairs == pairs_no_probs);    
    }

    SECTION("PairwiseAlignment produces empty AlignedPairs when sequences don't have any") {
        SymbolString sX = {{c, c, c, c}};
        SymbolString sY = {{t, t, t, t, t, t}}; 
        FiveStateSymbolHmm hmm;
        hmm.InitializeEmissions(SetNucleotideEmissionsToDefauts());
        AnchorPairs anchors = EmptyAnchors();
        AlignmentParameters p;
        p.expansion = 2;
        p.threshold = 0.8;
        PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln(hmm, sX, sY, anchors, p);
        AlignedPairs aligned_pairs = aln.AlignedPairsGetter();
        REQUIRE(aligned_pairs.size() == 0);
    }

    SECTION("PairwiseAlignment works on random related sequences", "[alignment]") {
        // generate some random (similar) sequences
        for (int64_t t = 0; t < 100; t++) {
            // TODO find solution to 0 length sequence
            std::string sX = RandomNucleotides(RandomInt(1, 100));
            std::string sY = EvolveSequence(sX);
            //st_uglyf("Seq X : %s\nSeq Y : %s\n", sX.c_str(), sY.c_str());

            AnchorPairs anchors = RandomAnchorPairs(sX.size(), sY.size());
    
            //st_uglyf("Got %" PRIi64 " random anchor pairs\n", anchors.size());
    
            SymbolString SsX = SymbolStringFromString(sX);
            SymbolString SsY = SymbolStringFromString(sY);
    
            AlignmentParameters p;
            p.expansion = 6;
            p.threshold = 0.1;
    
            FiveStateSymbolHmm hmm;
            hmm.InitializeEmissions(SetNucleotideEmissionsToDefauts());
        
            PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln(hmm, SsX, SsY, anchors, p);
            
            AlignedPairs aligned_pairs = aln.AlignedPairsGetter();

            CheckAlignedPairs(aligned_pairs, sX.size(), sY.size());
        }    
    }

    SECTION("PairwiseAlignment scoring makes sense", "[alignment]") {
        SymbolString sX = {{a, g, c, g}};
        SymbolString sY = {{a, g, t, t, c, g}}; 
        SymbolString sX2 = sX;
        FiveStateSymbolHmm hmm;
        hmm.InitializeEmissions(SetNucleotideEmissionsToDefauts());
        AlignmentParameters p;
        p.expansion = 2;
        p.threshold = 0.2;
        AnchorPairs anchors;
        PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln(hmm, sX, sY, anchors, p);
        PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln2(hmm, sX, sX2, anchors, p);
        double ungapped_score = aln.Score(true);
        double ungapped_score_same_seqs = aln2.Score(true);
        double gapped_score = aln.Score(false);
        double gapped_score_same_seqs = aln2.Score(false);
        REQUIRE(ungapped_score < ungapped_score_same_seqs);
        REQUIRE(gapped_score < gapped_score_same_seqs);
        REQUIRE(std::abs(gapped_score_same_seqs - ungapped_score_same_seqs) < 0.001);
    }
}

TEST_CASE("PairwiseAlignment works with ragged ends (local alignment)", "[current]") {
    SECTION("Local PairwiseAlignment works for simple sequences") {
        SymbolString sX =       {{a, g, c, g}};
        SymbolString sY = {{t, t, a, g, c, g, t, t}}; 
        
        FiveStateSymbolHmm hmm;
        hmm.InitializeEmissions(SetNucleotideEmissionsToDefauts());
        
        AnchorPairs anchors;
        anchors.emplace_back(0, 2);
        anchors.emplace_back(1, 3);
        anchors.emplace_back(3, 5); 
        
        AlignmentParameters p;
        p.expansion = 2;
        p.threshold = 0.2;
        PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln(hmm, sX, sY, anchors, p, true);
        AlignedPairs aligned_pairs = aln.AlignedPairsGetter();
        REQUIRE(aligned_pairs.size() == 4);
        CheckAlignedPairs(aligned_pairs, sX.size(), sY.size());   
        // extracts just the aligned pair coordinates, leaves the posteriors
        auto pairs_no_probs = [&aligned_pairs] () -> AnchorPairs {
            AnchorPairs pairs;
            for (auto p : aligned_pairs) {
                pairs.emplace_back(std::get<1>(p), std::get<2>(p));
            }
            return pairs;
        }();
        // make sure that worked
        REQUIRE(pairs_no_probs.size() == 4);

        AnchorPairs correct_pairs;
        correct_pairs.emplace_back(0, 2);
        correct_pairs.emplace_back(1, 3);
        correct_pairs.emplace_back(2, 4);
        correct_pairs.emplace_back(3, 5);
        REQUIRE(correct_pairs == pairs_no_probs);    
    }
}


