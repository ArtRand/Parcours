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
    std::default_random_engine generator(std::time(0));
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
