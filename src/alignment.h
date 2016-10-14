//
// Created by Arthur Rand on 10/9/16.
//

#ifndef PARCOURS_ALIGNMENT_H
#define PARCOURS_ALIGNMENT_H


#include "poa_graph.h"
#include "dp_matrix.h"

typedef enum _move {
    INS = 0,
    MATCH = 1,
    DEL = 2,
} Move;

class SimpleAlignment {
public:
    SimpleAlignment(Sequence *S, PoaGraph *G, int64_t (*SubFcn)(char, char));

    void AlignSequenceToGraph();

    void PrintScoresMatrix();

private:
    std::vector<int64_t> PredVertexIds(int64_t vId);

    int64_t (*SubstitutionFunction)(char i, char j);
    PoaGraph *graph;
    Sequence *sequence;
    std::unordered_map<int64_t, int64_t> IndexToId;
    std::unordered_map<int64_t, int64_t> IdToIndex;
    int64_t gapOpen;
    int64_t gapExtend;
    DpMatrix *scores;
    DpMatrix *bt_stringIdx;
    DpMatrix *bt_graphIdx;
};

class Op {
public:
    double score;
    int64_t graphPtr;
    int64_t seqPtr;
    uint8_t move;

    Op(double s, int64_t grPt, int64_t sqPt, uint8_t m) {
        score = s;
        graphPtr = grPt;
        seqPtr = sqPt;
        move = m;
    }

    ~Op() {};

    void ToString() {
        auto opString = [](int64_t o) -> std::string {
            switch (o) {
                case INS :
                    return "INSERT";
                case MATCH :
                    return "MATCH";
                case DEL :
                    return "DEL";
                default:
                    return "ERROR";
            }
        };
        std::cout << "[" << score << ", " << graphPtr << ", " << seqPtr << ", " << opString(move) << "]";
    }
};

int64_t BasicMatchFcn(char i, char j);

#endif //PARCOURS_ALIGNMENT_H
