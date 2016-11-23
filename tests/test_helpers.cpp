//
// Created by Arthur Rand on 10/6/16.
//
#include <chrono>
#include <random>
#include "test_helpers.h"

int64_t RandomInt(int64_t min, int64_t max) {
    return rand() % (max - min + 1) + min;
}

double RandomDouble() {
    thread_local static std::mt19937 rg{std::random_device{}()};
    std::default_random_engine generator(rg());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(generator);
}

// credit: http://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
static std::string randomString(int64_t length, std::string alphabet) {
    thread_local static std::mt19937 rg{std::random_device{}()};
    thread_local static std::uniform_int_distribution<> pick(0, alphabet.size() - 1);

    std::string s;

    s.reserve(length);

    while(length--)
        s += alphabet[pick(rg)];

    return s;
}

std::string RandomString(int64_t length) {
    static const std::string alphanums =
            "0123456789"
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    return randomString(length, alphanums);
}

std::string RandomNucleotides(int64_t length) {
    //static const std::string alphanums = "ACGT";
    static const std::string alphanums =
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"
                    "ACGTACGTACGTACGTACGTACGTACGT"
                    "ACGTACGTACGTACGTACGTACGTACGT";
    return randomString(length, alphanums);
}

std::string EvolveSequence(const std::string& startSequence) {
    //Copy sequence
    std::string seq = startSequence;

    //Do substitutions
    int subs = 0, indels = 0;
    for (uint64_t i = 0; i < seq.size(); i++) {
        if (RandomDouble() > 0.9) {
            char random_base = RandomNucleotides(1).at(0);
            seq.at(i) = random_base;
            subs++;
        }
    }
    
    //Do indels
    while (RandomDouble() > 0.1) {
        std::string to_replace = RandomNucleotides(RandomInt(2, 4));
        std::string replacement = RandomNucleotides(RandomInt(0, 10));
        seq = [&seq, &to_replace, &replacement, &indels] () -> std::string {  
            size_t found = seq.find(to_replace);
            if (found != std::string::npos) {
                seq.replace(found, to_replace.size(), replacement);
                indels++;
            }
            return seq;
        }();
    }
    //st_uglyf("EvolveSequence: made %i subs and %i indels\n", subs, indels);
    return seq;
}

AnchorPairs RandomAnchorPairs(int64_t lX, int64_t lY) {
    int64_t x = -1, y = -1;
    AnchorPairs pairs;
    while (1) {
        x += RandomInt(1, 20);
        y += RandomInt(1, 20);
        if (x >= lX || y >= lY) break;
        assert(x >= 0 && x < lX);
        assert(y >= 0 && y < lY); 
        pairs.emplace_back(x, y);
    }
    return pairs;
}

void CheckAlignedPairs(AlignedPairs pairs, int64_t lX, int64_t lY) {
    auto compare_aligned_pair = [] (AlignedPair p1, AlignedPair p2) -> bool {
        return std::get<0>(p1) == std::get<0>(p2) && 
               std::get<1>(p1) == std::get<1>(p2) &&
               std::get<2>(p1) == std::get<2>(p2);
    };
    
    // quickly check this
    AlignedPair _p  = std::make_tuple(1.0, 1, 1);
    AlignedPair _p2 = std::make_tuple(1.0, 1, 1);
    AlignedPair _p3 = std::make_tuple(1.0, 0, 1);
    AlignedPair _p4 = std::make_tuple(0.5, 0, 1);

    REQUIRE(compare_aligned_pair(_p, _p2));
    REQUIRE(!(compare_aligned_pair(_p, _p3)));
    REQUIRE(!(compare_aligned_pair(_p3, _p4)));

    //st_uglyf("got %lu pairs to check...", pairs.size());

    // check the validity of the pairs
    for (AlignedPair p : pairs) {
        REQUIRE(std::tuple_size<decltype(p)>::value == 3);  // correct length, trivial
        double score = std::get<0>(p);
        int64_t x    = std::get<1>(p);
        int64_t y    = std::get<2>(p);
        REQUIRE(score > 0);  // scores must be positive 
        REQUIRE(score <= PAIR_ALIGNMENT_PROB_1);  // scores cannot be more than this multiplier
        REQUIRE((x >= 0 && x < lX));  // pairs cannot be negative or 'beyond' the length of the sequence
        REQUIRE((y >= 0 && y < lY));  // pairs cannot be negative or 'beyond' the length of the sequence
    }
    // check that all pairs are unique
    auto compare_fcn = [] (AlignedPair p1, AlignedPair p2) -> bool {
        if (std::get<1>(p1) == std::get<1>(p2)) {
            return std::get<2>(p1) < std::get<2>(p2);
        } else {
            return std::get<1>(p1) < std::get<1>(p2);
        }
    };
    
    std::sort(begin(pairs), end(pairs), compare_fcn);
    
    AlignedPair p = pairs.at(0);
    for (uint64_t i = 1; i < pairs.size(); i++) {
        REQUIRE(!(compare_aligned_pair(pairs.at(i), p)));
        p = pairs.at(i);
    }
    //st_uglyf("OK\n");
}
