//
// Created by Arthur Rand on 10/6/16.
//
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "test_helpers.h"

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
        delete G;
        delete S;
    }


}
