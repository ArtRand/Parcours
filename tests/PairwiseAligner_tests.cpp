// Pairwise aligner tests
//


#include "catch.hpp"
#include "diagonal.h"
#include "band.h"
#include "common.h"

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
    Band b(anchors, lX, lY, expansion);
    
    // Forward pass
    
    Diagonal d(0, 0, 0); REQUIRE(b.Next() == d);
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(2, -2, 2)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(3, -1, 3)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(4, -2, 4)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(5, -1, 3)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(6, -2, 4)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));
    //CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));

}

