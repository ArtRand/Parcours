//
// Created by Arthur Rand on 10/9/16.
//

#include "alignment.h"
#include "common.h"

static inline void indexMapMaker(PoaGraph *G, std::unordered_map<int64_t, int64_t> &IndexToId,
                                 std::unordered_map<int64_t, int64_t> &IdToIndex) {
    if (!G->IsSorted()) {
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

SimpleAlignment::SimpleAlignment(Sequence *S, PoaGraph *G, int64_t (*SubFcn)(char, char)): gapOpen(-4),
                                                                                           gapExtend(-2),
                                                                                           is_aligned(false) {
    if (S->seq.size() <= 0) {
        throw GraphException("SimpleAlignment constructor, passed empty sequence string error");
    }
    sequence = S;
    graph = G;
    substitutionFunction = SubFcn;

    // dp_matrix_setup vertexId-to-index and index-to-vertexId maps
    indexMapMaker(G, indexToId, idToIndex);

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
        pred.push_back(idToIndex[kv.first]);
    }

    if (pred.size() == 0) {
        pred.push_back(-1);
    }

    return pred;
}

static inline bool compareOp(Op *o, Op *p) {
    // if the scores are != than just compare by score
    if (o->score != p->score) {
        return o->score < p->score;
    }
    // iff the scores are the same && one is a match, prefer a match
    if ((o->move == MATCH) || (p->move == MATCH)) {
        return o->move != MATCH;
    }
    // if the scores are the same and neither is a match, just return arbitrary choice
    return o->score < p->score;
}
/*  old version without "pick match" logic
static inline bool compareOp(Op *o, Op *p) {
    if (o->score == p->score) {
        st_uglyf("same!");
    }
    return o->score < p->score;
}
*/

void SimpleAlignment::AlignSequenceToGraph() {
    if (!graph->IsSorted()) {
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
                double matchScore = scores->Getter(pred + 1, j) + substitutionFunction(sBase, pBase);
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
    traceback_global();

    is_aligned = true;
}

void SimpleAlignment::traceback_global() {
    int64_t best_j = scores->Cols() - 1;  // minus one to make 0-indexed
    int64_t best_i = scores->Rows() - 1;  // see -^^

    // find the terminal indices
    std::vector<int64_t> terminals;

    //st_uglyf("before loop-> best_i: %lld, best_j: %lld\n", best_i, best_j);

    for (int64_t vId : graph->Vertices()) {
        if (graph->VertexGetter(vId)->OutDegree() == 0) {  // it's an "end" or "terminal" vertex
            terminals.push_back(vId);
        }
    }
    // there should be at least 1 terminal vertex
    if (terminals.size() < 1) {
        st_uglyf("Got %lld terminals\n", terminals.size());
        throw GraphException("SimpleAlignment::traceback_global - terminals less than 1 error\n");
    }
    best_i = terminals.at(0);
    double bestScore = scores->Getter(best_i, best_j);
    // check the other terminals for better scores
    for (uint64_t i = 1; i < terminals.size(); i++) {  // I could not find a more c++11'y way to do this !!?
        double s = scores->Getter(terminals.at(i), best_j);
        // TODO check if >= is better? longer alignments?
        if (s > bestScore) {
            bestScore = s;
            best_i = terminals.at(i);
        }
    }

    //st_uglyf("after loop-> best_i: %lld, best_j: %lld, bestScore: %f\n", best_i, best_j, bestScore);

    // TODO can probably remove this later
    if (matches.size() != 0) {
        throw GraphException("SimpleAlignment::traceback_global - Matches vector is not cleared");
    }
    if (strIdxs.size() != 0) {
        throw GraphException("SimpleAlignment::traceback_global - String Indices vector is not cleared");
    }

    //bt_graphIdx->ToString();
    //std::cout << "\n";
    //bt_stringIdx->ToString();
    //std::cout << "\n";

    while (!(best_i == 0 && best_j == 0)) {
        int64_t next_i = (int64_t )bt_graphIdx->Getter(best_i, best_j);
        int64_t next_j = (int64_t )bt_stringIdx->Getter(best_i, best_j);
        int64_t cur_str_idx = best_j - 1;
        int64_t cur_vertex_id = indexToId[best_i - 1];

        strIdxs.push_front((next_j != best_j ? cur_str_idx : -1));
        matches.push_front((next_i != best_i ? cur_vertex_id : -1));

        best_i = next_i;
        best_j = next_j;
    }
    /* 
    st_uglyf("strIdxs: ");
    for (auto i : strIdxs) {
        st_uglyf("%lld, ", i);
    }
    std::cout << "\n";
    st_uglyf("matches: ");
    for (auto i : matches) {
        st_uglyf("%lld, ", i);
    }
    std::cout << "\n";
    */
}

std::pair<std::string, std::string> SimpleAlignment::AlignmentStrings() {
    if (!is_aligned) {
        AlignSequenceToGraph();
        assert(is_aligned);
    }

    std::string seqStr = "";
    std::string grphStr = "";

    for (auto i : strIdxs) {
        //seqStr += (i >= 0 ? sequence->seq.at(i) : "-");
        if (i >= 0) {
            seqStr += sequence->seq.at(i);
        } else {
            seqStr += "-";
        }
    }
    for (auto j : matches) {
        if (j >= 0) {
            grphStr += graph->VertexGetter(j)->Base();
        } else {
            grphStr += "-";
        }
    }

    return std::make_pair(seqStr, grphStr);
}

static inline std::vector<int64_t> AlignedStringIdxs(const std::deque<int64_t>& string_indices) {
    std::vector<int64_t> validStringIdxs (string_indices.size());
    // iterator pointing to the start of validStringIdxs that contains all of the aligned string indexes
    auto it = std::copy_if(string_indices.begin(), string_indices.end(), validStringIdxs.begin(),
                           [] (int64_t i) { return i >= 0;});

    validStringIdxs.resize(std::distance(validStringIdxs.begin(), it));

    return validStringIdxs;
}

void SimpleAlignment::AddAlignmentToGraph() {
    int64_t firstId = -1;
    int64_t headId = -1;
    int64_t tailId = -1;

    std::vector<int64_t> aligned_string_idxs = AlignedStringIdxs(strIdxs);

    int64_t startSeqIdx = aligned_string_idxs.front();
    int64_t endSeqIdx = aligned_string_idxs.back();

    // first we deal with un-aligned sequence in the 'string' or read being aligned to the graph
    if (startSeqIdx > 0) {  // there is a left ragged-end
        std::string unaligned = sequence->seq.substr(0, startSeqIdx - 1);
        graph->AddBaseSequence(unaligned, sequence->label, false, firstId, headId);
    }
    if (endSeqIdx < (int64_t )sequence->seq.size()) {  // there is a right ragged end
        int64_t slice_start = endSeqIdx + 1;
        int64_t extend_length = sequence->seq.size() - slice_start;
        std::string unaligned = sequence->seq.substr(slice_start, extend_length);
        int64_t _;
        graph->AddBaseSequence(unaligned, sequence->label, false, tailId, _);
    }

    // check for matches and stringIdxs being a different length, shouldn't be the case...
    if (strIdxs.size() != matches.size()) {
        throw GraphException("string indices and matches are different sizes, can't have that");
    }

    for (uint64_t j = 0; j < matches.size(); j++) {
        int64_t string_idx = strIdxs.at(j);  // the index of the string we're looing at
        int64_t vertex_id = matches.at(j);   // the vertex id in the graph we're looking at

        if (string_idx == -1) {
            continue;
        }

        char base = sequence->seq.at(string_idx);

        int64_t newNodeId;

        // if there isn't an alignment to this vertex ??
        if (vertex_id == -1) {
            newNodeId = graph->AddVertex(base);
        } else if (graph->VertexGetter(vertex_id)->Base() == base) {  // if there is a match
            newNodeId = vertex_id;
        } else {
            int64_t foundMatch = -1;

            // search all of the vertices this vertex is aligned to for a match
            for (int64_t al : graph->VertexGetter(vertex_id)->aligned_to) {
                if (graph->VertexGetter(al)->Base() == base) {
                    foundMatch = al;
                }
            }

            // didn't find a match
            if (foundMatch == -1) {
                // add new vertex to the graph
                newNodeId = graph->AddVertex(base); 
                // align it to this vertex 
                graph->VertexGetter(newNodeId)->aligned_to.push_back(vertex_id);
                // and the vertices it's aligned to
                for (int64_t al: graph->VertexGetter(vertex_id)->aligned_to) {
                    graph->VertexGetter(newNodeId)->aligned_to.push_back(al);
                }
                // now update all the vertices it's aligned to with this new vertex id
                for (int64_t otherId : graph->VertexGetter(newNodeId)->aligned_to) {
                    graph->VertexGetter(otherId)->aligned_to.push_back(newNodeId);
                }
            } else {  // found a match
                newNodeId = foundMatch;
            }
        }
        // connect the new node, or just update the label
        graph->AddArc(headId, newNodeId, sequence->label);
        headId = newNodeId;  // update for the next loop
        if (firstId == -1) {  // keep track of the new start to update later
            firstId = headId;
        }
    }

    graph->AddArc(headId, tailId, sequence->label);

    // resort
    graph->TopologicalSort();
    assert(graph->TestSort());

    graph->AddSequence(sequence->seq);
    graph->AddLabel(sequence->label);
    graph->AddStart(firstId);
}

int64_t BasicMatchFcn(char i, char j) { return i == j ? 4 : -2; }

