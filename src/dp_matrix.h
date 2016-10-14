//
// Created by Arthur Rand on 10/11/16.
//

#ifndef PARCOURS_DP_MATRIX_H
#define PARCOURS_DP_MATRIX_H


#include "stl_includes.h"

// ------->i--|
// |
// |
// j
// |
// |
class DpMatrix {
public:
    DpMatrix(int64_t i, int64_t j);

    DpMatrix(int64_t i, int64_t j, double defVal);

    ~DpMatrix();

    double Getter(int64_t , int64_t);

    void Setter(int64_t , int64_t, double);

    void ToString();

    int64_t Cols() const;

    int64_t Rows() const;

    int64_t Size() const;

private:

    void dp_matrix_setup(int64_t r, int64_t c);

    int64_t i, j;
    int64_t size;
    std::vector<double> matrix;
};




#endif //PARCOURS_DP_MATRIX_H
