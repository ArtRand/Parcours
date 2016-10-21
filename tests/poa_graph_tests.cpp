//
// Created by Arthur Rand on 10/6/16.
//
#define CATCH_CONFIG_MAIN

#include "alignment.h"
#include "dp_matrix.h"
#include "catch.hpp"
#include "test_helpers.h"
#include "common.h"

TEST_CASE("Test vertex and arcs") {
    for (int64_t i = 0; i < 100; i++) {
        Vertex *v = new Vertex();
  
        // check proper initialization
        REQUIRE(v->Id() == -1);

        v->SetId(i);

        // check setter
        REQUIRE(v->Id() == i);

        delete v;
    }

    DirectedArc *a;
    a = new DirectedArc();

    // check proper initialization
    REQUIRE(a->To() == -1);
    REQUIRE(a->From() == -1);
    REQUIRE(a->LabelSet().size() == 0);

    for (int64_t i = 1; i < 100; i++) {
        std::string l = RandomString(10);

        // check add label logic and container size
        REQUIRE(!a->AddLabel(l));
        REQUIRE(a->AddLabel(l));
        REQUIRE(a->LabelSet().size() == i);
    }

    delete a;

    for (int64_t i = 0; i < 100; i++) {
        Vertex *v = new Vertex(i);
        Vertex *w = new Vertex(i + 1);
        a = new DirectedArc(v->Id(), w->Id());
        std::string label = RandomString(10);
        v->AddOutArc(w->Id(), label);
        w->AddInArc(v->Id(), label);

        // check in and out degree
        REQUIRE(v->out_arcs.size() == 1);
        REQUIRE(v->OutDegree() == 1);
        REQUIRE(v->in_arcs.size() == 0);
        REQUIRE(v->InDegree() == 0);

        REQUIRE(w->out_arcs.size() == 0);
        REQUIRE(w->OutDegree() == 0);
        REQUIRE(w->in_arcs.size() == 1);
        REQUIRE(w->InDegree() == 1);

        delete v;
        delete w;
        delete a;
    }
}

TEST_CASE("PoaGraph Basic") {
    // trivial case
    PoaGraph *G = new PoaGraph();
    REQUIRE(G->K() == 0);
    delete G;

    // make test base sequence and label
    for (int64_t i = 0; i < 100; i++) {
        int64_t seqLen = RandomInt(100, 1000);
        std::string label = RandomString(10);
        std::string sequence = RandomNucleotides(seqLen);
        Sequence *S = new Sequence();
        S->seq = sequence;
        S->label = label;
        G = new PoaGraph(*S);

        // test we made the correct number of nodes
        REQUIRE(G->K() == seqLen);
        
         // test the base sequence
        for (int64_t i = 0; i < G->Vertices().size(); i++) {
            REQUIRE(G->VertexGetter(i)->Base() == sequence.at(i));
        }

        // check the start node
        REQUIRE(G->Starts().size() == 1);
        REQUIRE(G->VertexGetter(0)->InDegree() == 0);
        REQUIRE(G->VertexGetter(0)->OutDegree() == 1);

        // check the end node
        REQUIRE_THROWS_AS(G->VertexGetter(seqLen), GraphException);
        REQUIRE_THROWS_AS(G->VertexGetter(seqLen), GraphException);

        REQUIRE(G->VertexGetter(seqLen - 1)->InDegree() == 1);
        REQUIRE(G->VertexGetter(seqLen - 1)->OutDegree() == 0);

        for (int64_t j = 1; j < seqLen - 1; j++) {
            REQUIRE(G->VertexOutDegree(j) == 1);
            REQUIRE(G->VertexInDegree(j) == 1);
        }
        
        // cleanup 
        delete G;
        delete S;
    }
}

TEST_CASE("PoaGraph TestSort") {
    //  0    1    2    3    4    5    6
    // (A)->(C)->(G)->(A)->(G)->(T)->(A)
    //        \____(A)->(C)____/^
    //              7    8
    std::string label = RandomString(10);
    std::string label2 = RandomString(10);
    std::string sequence = "ACGAGTA";
    Sequence *S = new Sequence();
    S->seq = sequence;
    S->label = label;
    PoaGraph *G = new PoaGraph(*S);

    REQUIRE(G->K() == 7);

    int64_t v = G->AddVertex('A');
    int64_t w = G->AddVertex('C');

    REQUIRE(G->K() == 9);
    REQUIRE(v == 7);
    REQUIRE(w == 8);

    REQUIRE(G->VertexGetter(1)->Base() == 'C');
    REQUIRE(G->VertexGetter(5)->Base() == 'T');

    G->AddArc(1, v, label2);
    G->AddArc(v, w, label2);
    G->AddArc(w, 5, label2);

    REQUIRE(!G->TestSort());
    REQUIRE(!G->IsSorted());

    G->TopologicalSort();

    REQUIRE(G->IsSorted());
    REQUIRE(G->TestSort());

    delete G;
}

TEST_CASE("Test DpMatrix") {
    DpMatrix *m = new DpMatrix(10, 10);
    for (uint64_t u = 10; u < 20; u++) {
        REQUIRE_THROWS_AS(m->Getter(u, 10), GraphException);
        REQUIRE_THROWS_AS(m->Getter(10, u), GraphException);
        REQUIRE_THROWS_AS(m->Setter(u, 10, 0.0), GraphException);
        REQUIRE_THROWS_AS(m->Setter(10, u, 0.0), GraphException);
    }
    for (uint64_t u = 0; u < 50; u++) {
        int64_t x = RandomInt(0, 9);
        int64_t y = RandomInt(0, 9);
        double val = rand();
        REQUIRE_NOTHROW(m->Setter(x, y, val));
        REQUIRE(m->Getter(x, y) == val);
    }
    delete m;
}

TEST_CASE("Test Initialize Simple Sequence Alignment With Deletes") {
    std::string base_string = "ACAAATAG";
    std::string base_label = RandomString(5);
    Sequence *base = new Sequence;
    base->seq = base_string;
    base->label = base_label;

    PoaGraph *G = new PoaGraph(*base);
    
    REQUIRE(G->K() == 8);
    delete base;

    std::string s = "ACATAG";
    std::string l = RandomString(5);
    Sequence *S = new Sequence;
    S->seq = s;
    S->label = l;

    SimpleAlignment *A = new SimpleAlignment(S, G, BasicMatchFcn);
    A->AlignSequenceToGraph();

    // check alignment strings
    std::pair<std::string, std::string> results = A->AlignmentStrings();
    REQUIRE(results.first == "AC--ATAG");
    REQUIRE(results.second == "ACAAATA-");
    
    A->AddAlignmentToGraph();

    // manually check the graph
    REQUIRE(G->K() == 9);  // added one node to the graph
    REQUIRE(G->IsSorted());
    REQUIRE(G->TestSort());
    // walk though the graph and check topology

    auto checkGraph0 = [G, l, base_label]() {
        int64_t checking = 0;
        // 0
        REQUIRE(G->VertexGetter(checking)->Base() == 'A');
        REQUIRE(G->VertexInDegree(checking) == 0);
        REQUIRE(G->VertexOutDegree(checking) == 1);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 1
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(1) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().size() == 2);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().count(l) == 1);

        // 1
        checking = 1;
        REQUIRE(G->VertexGetter(checking)->Base() == 'C');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 2);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 2
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(2) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().size() == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().count(l) == 0);
        // check connection to 4
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(4) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().size() == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(base_label) == 0);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l) == 1);

        // 2
        checking = 2;
        REQUIRE(G->VertexGetter(checking)->Base() == 'A');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 1);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 3
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(3) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().size() == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().count(l) == 0);

        // 3
        checking = 3;
        REQUIRE(G->VertexGetter(checking)->Base() == 'A');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 1);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 4
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(4) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().size() == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l) == 0);

        // 4
        checking = 4;
        REQUIRE(G->VertexGetter(checking)->Base() == 'A');
        REQUIRE(G->VertexInDegree(checking) == 2);
        REQUIRE(G->VertexOutDegree(checking) == 1);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 5
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(5) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().size() == 2);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().count(l) == 1);

        // 5
        checking = 5;
        REQUIRE(G->VertexGetter(checking)->Base() == 'T');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 1);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 6
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(6) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().size() == 2);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().count(l) == 1);

        // 6
        checking = 6;
        REQUIRE(G->VertexGetter(checking)->Base() == 'A');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 2);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
        // check connection to 7
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(7) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().size() == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().count(base_label) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().count(l) == 0);
        // check connection to 8
        REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(8) == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[8]->LabelSet().size() == 1);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[8]->LabelSet().count(base_label) == 0);
        REQUIRE(G->VertexGetter(checking)->OutNeighbors()[8]->LabelSet().count(l) == 1);

        // 7
        checking = 7;
        REQUIRE(G->VertexGetter(checking)->Base() == 'G');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 0);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);

        // 8
        checking = 8;
        REQUIRE(G->VertexGetter(checking)->Base() == 'G');
        REQUIRE(G->VertexInDegree(checking) == 1);
        REQUIRE(G->VertexOutDegree(checking) == 0);
        REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    };

    checkGraph0();

    // delete this alignment and it's associated sequence
    delete A;
    delete S;

    // check graph is still correct, and no seg-faults
    checkGraph0();

    // align another sequence
    std::string s2 = "ACATATAG";
    std::string l2 = RandomString(6);
    S = new Sequence;
    S->seq = s2;
    S->label = l2;

    A = new SimpleAlignment(S, G, BasicMatchFcn);
    A->AlignSequenceToGraph();
    results = A->AlignmentStrings();
    REQUIRE(results.first == "ACATATAG");
    REQUIRE(results.second == "ACAAATAG");
    A->AddAlignmentToGraph();
    delete A;
    delete S;

    // should have added another aligned-to node
    REQUIRE(G->K() == 10);
    // check

    int64_t checking = 0;
    // 0
    REQUIRE(G->VertexGetter(checking)->Base() == 'A');
    REQUIRE(G->VertexInDegree(checking) == 0);
    REQUIRE(G->VertexOutDegree(checking) == 1);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    // check connection to 1
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(1) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().size() == 3);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().count(l) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[1]->LabelSet().count(l2) == 1);

    // 1
    checking = 1;
    REQUIRE(G->VertexGetter(checking)->Base() == 'C');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 2);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    // check connection to 2
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(2) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().size() == 2);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().count(l2) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[2]->LabelSet().count(l) == 0);
    // check connection to 4
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(4) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().size() == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(base_label) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l) == 1);

    // 2
    checking = 2;
    REQUIRE(G->VertexGetter(checking)->Base() == 'A');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 2);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    // check connection to 3
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(3) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().size() == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().count(l) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[3]->LabelSet().count(l2) == 0);
    // check connection to 9
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(9) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[9]->LabelSet().size() == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[9]->LabelSet().count(base_label) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[9]->LabelSet().count(l) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[9]->LabelSet().count(l2) == 1);

    // 3
    checking = 3;
    REQUIRE(G->VertexGetter(checking)->Base() == 'A');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 1);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 1);
    REQUIRE(std::find(G->VertexGetter(checking)->aligned_to.begin(),
                      G->VertexGetter(checking)->aligned_to.end(), 9) != G->VertexGetter(checking)->aligned_to.end());
    // check connection to 4
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(4) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().size() == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l2) == 0);

    // 4
    checking = 4;
    REQUIRE(G->VertexGetter(checking)->Base() == 'A');
    REQUIRE(G->VertexInDegree(checking) == 3);
    REQUIRE(G->VertexOutDegree(checking) == 1);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    // check connection to 5
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(5) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().size() == 3);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().count(l) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[5]->LabelSet().count(l2) == 1);

    // 5
    checking = 5;
    REQUIRE(G->VertexGetter(checking)->Base() == 'T');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 1);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    // check connection to 6
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(6) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().size() == 3);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().count(l) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[6]->LabelSet().count(l2) == 1);

    // 6
    checking = 6;
    REQUIRE(G->VertexGetter(checking)->Base() == 'A');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 2);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);
    // check connection to 7
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(7) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().size() == 2);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().count(base_label) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().count(l) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[7]->LabelSet().count(l2) == 1);
    // check connection to 8
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(8) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[8]->LabelSet().size() == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[8]->LabelSet().count(base_label) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[8]->LabelSet().count(l) == 1);

    // 7
    checking = 7;
    REQUIRE(G->VertexGetter(checking)->Base() == 'G');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 0);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);

    // 8
    checking = 8;
    REQUIRE(G->VertexGetter(checking)->Base() == 'G');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 0);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 0);

    // 9
    checking = 9;
    REQUIRE(G->VertexGetter(checking)->Base() == 'T');
    REQUIRE(G->VertexInDegree(checking) == 1);
    REQUIRE(G->VertexOutDegree(checking) == 1);
    REQUIRE(G->VertexGetter(checking)->aligned_to.size() == 1);
    REQUIRE(std::find(G->VertexGetter(checking)->aligned_to.begin(),
                      G->VertexGetter(checking)->aligned_to.end(), 3) != G->VertexGetter(checking)->aligned_to.end());
    // check connection to 4
    REQUIRE(G->VertexGetter(checking)->OutNeighbors().count(4) == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().size() == 1);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(base_label) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l) == 0);
    REQUIRE(G->VertexGetter(checking)->OutNeighbors()[4]->LabelSet().count(l2) == 1);

    delete G;
}


