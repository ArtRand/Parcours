//
// Created by Arthur Rand on 10/11/16.
//

#include "dp_matrix.h"
#include "parcours_exceptions.h"

void DpMatrix::dp_matrix_setup(int64_t r, int64_t c) {
    i = r;
    j = c;
    size = r * c;
}

DpMatrix::DpMatrix(int64_t r, int64_t c) {
    dp_matrix_setup(r, c);
    //i = r;
    //j = c;
    //size = r * c;
    matrix = std::vector<double>(i * j, 0.0);
}

DpMatrix::DpMatrix(int64_t r, int64_t c, double defVal) {
    dp_matrix_setup(r, c);
    //i = r;
    //j = c;
    //size = r * c;
    matrix = std::vector<double>(i * j, defVal);
}

DpMatrix::~DpMatrix() {}

static inline bool checkCoords(int64_t r, int64_t c, int64_t size) {
    return r * c < size;
}

double DpMatrix::Getter(int64_t r, int64_t c) {
    if (!checkCoords(r, c , size)) {
        throw ParcoursException("DpMatrix::Getter - out of bounds\n");
    }
    return matrix[r * j + c];
}

void DpMatrix::Setter(int64_t r, int64_t c, double val) {
    if (!checkCoords(r, c , size)) {
        throw ParcoursException("DpMatrix::Setter - out of bounds\n");
    }
    matrix[r * j + c] = val;
}

void DpMatrix::ToString() {
    for (int64_t row = 0; row < i; row++) {
        for (int64_t col = 0; col < j; col++) {
            std::cout << Getter(row, col) << " ";
        }
        std::cout << std::endl;
    }
}

int64_t DpMatrix::Cols() const { return j; }

int64_t DpMatrix::Rows() const { return i; }

int64_t DpMatrix::Size() const { return size; }
