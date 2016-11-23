//
// Created by Arthur Rand on 10/6/16.
//

#ifndef PARCOURS_TEST_HELPERS_H
#define PARCOURS_TEST_HELPERS_H

#include "pairwise_aligner.h"
#include "catch.hpp"


int64_t RandomInt(int64_t min, int64_t max);
std::string RandomString(int64_t length);
double RandomDouble();
std::string RandomNucleotides(int64_t length);
std::string EvolveSequence(const std::string& startSequence); 
AnchorPairs RandomAnchorPairs(int64_t lX, int64_t lY);
void CheckAlignedPairs(AlignedPairs pairs, int64_t lX, int64_t lY);

#endif //PARCOURS_TEST_HELPERS_H
