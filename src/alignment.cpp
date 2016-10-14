//
// Created by Arthur Rand on 10/9/16.
//

#include "alignment.h"

static inline void indexMapMaker(PoaGraph *G, std::unordered_map<int64_t, int64_t> &IndexToId,
                                 std::unordered_map<int64_t, int64_t> &IdToIndex) {
    if (!G->isSorted()) {
        G->TopologicalSort();
    }
    assert(G->TestSort());

    IndexToId[-1] = -1;

    uint64_t i = 0;
    for (int64_t vid : G->Vertices()) {  // loop over the *sorted* vertex IDs
        IndexToId[i] = vid;
        IdToIndex[vid] = i;
        i++;
    }
};

static inline void initializeDpMatrixGlobal(DpMatrix *M, int64_t gapOpen, int64_t gapExtend) {
    int64_t cols = M->Cols();
    int64_t rows = M->Rows();
    
    auto globGapScoreInit = [gapOpen, gapExtend] (int64_t step) -> double {
        return gapOpen + (step - 1) * gapExtend;
    };

    for (int64_t i = 1; i < rows; i++) {
        M->Setter(i, 0, globGapScoreInit(i));
    }

    for (int64_t j = 1; j < cols; j++) {
        M->Setter(0, j, globGapScoreInit(j));
    }

}

SimpleAlignment::SimpleAlignment(Sequence *S, PoaGraph *G, int64_t (*SubFcn)(char, char)): gapOpen(-4), gapExtend(-2) {
    sequence = S;
    graph = G;
    SubstitutionFunction = SubFcn;

    // dp_matrix_setup vertexId-to-index and index-to-vertexId maps
    indexMapMaker(G, IndexToId, IdToIndex);

    // dp_matrix_setup dynamic programming data structures
    int64_t l2 = S->seq.length() + 1;
    int64_t l1 = G->K() + 1;

    scores = new DpMatrix(l1, l2);
    initializeDpMatrixGlobal(scores, gapOpen, gapExtend);
    bt_stringIdx = new DpMatrix(l1, l2);
    bt_graphIdx = new DpMatrix(l1, l2);
}

void SimpleAlignment::PrintScoresMatrix() {
    scores->ToString();
}

std::vector<int64_t> SimpleAlignment::PredVertexIds(int64_t vId) {
    // instantiate the vector to the size of the in arcs in the vertex we're at
    Vertex *v = graph->VertexGetter(vId);
    std::vector<int64_t> pred;

    // iterate over the in-neighbors vertexIDs
    for (auto kv : v->in_arcs) {  // kv = (vertexId, DirectedArc*)
        // add the dp matrix index this node is at to Pred
        pred.push_back(IdToIndex[kv.first]);
    }

    if (pred.size() == 0) {
        pred.push_back(-1);
    }

    return pred;
}

static inline bool compareOp(Op *o, Op *p) {
    return o->score < p->score;
}

void SimpleAlignment::AlignSequenceToGraph() {
    if (!graph->isSorted()) {
        graph->TopologicalSort();
        assert(graph->TestSort());
    }

    DpMatrix *insertCost = new DpMatrix(scores->Rows(), scores->Cols(), gapOpen);
    DpMatrix *deleteCost = new DpMatrix(scores->Rows(), scores->Cols(), gapOpen);

    int64_t i = 0;
    for (auto vId : graph->Vertices()) {
        Vertex *v = graph->VertexGetter(vId);
        char pBase = v->Base();

        int64_t j = 0;
        for (auto sBase : sequence->seq) {
            std::vector<Op*> candidates;  // array of move options

            // make the candidate operation for inserting a base from the sequence
            double insertScore = scores->Getter(i + 1, j) + insertCost->Getter(i + 1, j);
            Op *insOp = new Op(insertScore, i + 1, j, INS);
            candidates.push_back(insOp);

            for (int64_t pred : PredVertexIds(vId)) {
                // handle match
                double matchScore = scores->Getter(pred + 1, j) + SubstitutionFunction(sBase, pBase);
                Op *matchOp = new Op(matchScore, pred + 1, j, MATCH);

                // handle delete
                double deleteScore = scores->Getter(pred + 1, j + 1) + deleteCost->Getter(pred + 1, j + 1);
                Op *delOp = new Op(deleteScore, pred + 1, j + 1, DEL);

                candidates.insert(candidates.end(), {matchOp, delOp});
            }

            Op *maxOp = *std::max_element(candidates.begin(), candidates.end(), compareOp);

            scores->Setter(i + 1, j + 1, maxOp->score);
            bt_graphIdx->Setter(i + 1, j + 1, maxOp->graphPtr);
            bt_stringIdx->Setter(i + 1, j + 1, maxOp->seqPtr);

            // debug print out candidates
            /*
            std::cout << i << ", " << j;
            for (Op* op: candidates) {
                op->ToString();
                std::cout << " ";
            }
            std::cout << "max: "; maxOp->ToString(); std::cout << maxOp->score;
            std::cout << std::endl;
            */

            if (maxOp->move == INS) {
                insertCost->Setter(i + 1, j + 1, gapExtend);
            }
            if (maxOp->move == DEL) {
                deleteCost->Setter(i + 1, j + 1, gapExtend);
            }
            j++;
        }
        i++;
    }
    //scores->ToString(); std::cout << std::endl;
    //bt_graphIdx->ToString(); std::cout << std::endl;
    //bt_stringIdx->ToString();
    //deleteCost->ToString();

}

int64_t BasicMatchFcn(char i, char j) {
    return i == j ? 4 : -2;
}

