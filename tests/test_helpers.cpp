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
        if (RandomDouble() > 0.8) {
            char random_base = RandomNucleotides(1).at(0);
            seq.at(i) = random_base;
            subs++;
        }
    }
    
    //Do indels
    while (RandomDouble() > 0.2) {
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
    st_uglyf("EvolveSequence: made %i subs and %i indels\n", subs, indels);
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
