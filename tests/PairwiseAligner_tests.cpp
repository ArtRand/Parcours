// Pairwise aligner tests
//

#include "logAdd.h"
#include "catch.hpp"
#include "diagonal.h"
#include "band.h"
#include "common.h"
#include "dpMatrix.h"
#include "stateMachine.h"
#include "test_helpers.h"

TEST_CASE("Diagonal Tests", "[PairwiseAligner]") {
    int64_t xL = 10, yL = 20, xU = 30, yU = 0; //Coordinates of the upper and lower
    Diagonal d(xL + yL, xL - yL, xU - yU);
    Diagonal d2(0, 0, 0);
    REQUIRE(d.Xay() == xL + yL);
    REQUIRE(d.MinXmy() == xL - yL);
    REQUIRE(d.MaxXmy() == xU - yU);
    REQUIRE(d.Width() == ((xU - yU - (xL - yL)) / 2 + 1));
    REQUIRE(diagonal_XCoordinate(xL + yL, xL - yL) == xL);
    REQUIRE(diagonal_YCoordinate(xL + yL, xL - yL) == yL);
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
        double l = exp(logAdd(log(i), log(j)));
        REQUIRE(l < k + 0.001);
        REQUIRE(l > k - 0.001);
    }
}

TEST_CASE("Test Cell", "[DpTests]") {
    StateMachine5<nucleotide> sM5;
    REQUIRE(sM5.StateNumber() == 5);

    
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
            total_prob = logAdd(total_prob, 2 * c1);
            total_prob = logAdd(total_prob, 2 * c2);
        }
        DpDiagonal<double, 5> d2 = d;
        double tp = d.Dot(d2);
        REQUIRE(std::abs(total_prob - tp) < 0.001);
    }
}

TEST_CASE("Test DpMatrix", "[DpTests]") {
    SECTION("DpMatrix is constructed with correct defaults") {
        int64_t x = RandomInt(10, 20);
        DpMatrix<double, fiveState> mat(x);
        REQUIRE(mat.ActiveDiagonals() == 0);
        REQUIRE(mat.DiagonalNumber() == x);
    }
    SECTION("DpMatrix main tests") {
        int64_t lX = 3;
        int64_t lY = 2;

        DpMatrix<double, fiveState> mat(lX + lY);

        REQUIRE(mat.ActiveDiagonals() == 0);
        
        // check for 'fantom' diagonals
        for (int64_t i = -1; i <= lX + lY + 10; i++) {
            REQUIRE_THROWS_AS(mat.DpDiagonalGetter(i), ParcoursException);
        }

        // make some diagonals in the dpMatrix, and check them, then make sure that
        // the number of active diagonals is correct.
        for (int64_t i = 0; i <= lX + lY; i++) {
            mat.CreateDpDiagonal(i, -i, i);
            DpDiagonal<double, fiveState> d(i, -i, i);
            REQUIRE(mat.DpDiagonalGetter(i) == d);
            REQUIRE(mat.ActiveDiagonals() == i + 1);
        }

        // delete the diagonals
        for (int64_t i = lX + lY; i >= 0; i--) {
            mat.DeleteDpDiagonal(i);
            REQUIRE(!mat.DpDiagonalGetter(i).IsActive());
            REQUIRE(mat.ActiveDiagonals() == i);
        }
    }
}
