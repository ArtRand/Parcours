//
// Created by Arthur Rand on 10/24/16.
//

#define CATCH_CONFIG_MAIN

#include "Vertex.h"
#include "hmm_graph.h"
#include "catch.hpp"
#include "test_helpers.h"

TEST_CASE("Basic Object Tests", "[lib]") {
    SECTION("Vertex objects are created and destroyed correctly") {
        // test initializing with empty string throws exception
        std::string s = "";
        REQUIRE_THROWS_AS(Vertex *v = new Vertex(0, &s), ParcoursException);
        
        // test with random strings
        for (int64_t t = 0; t < 10; t++) {
            std::string s = RandomNucleotides(RandomInt(10, 20));
            Vertex *v = new Vertex(0, &s);
            REQUIRE(*v->Sequence() == s);
            delete v;
        }
    }
    
    SECTION("Vertices show correct adjacency") {
        std::string x = RandomNucleotides(RandomInt(10, 100));
        std::string y = RandomNucleotides(RandomInt(10, 100));
        std::string z = RandomNucleotides(RandomInt(10, 100));
        Vertex v = Vertex(0, &x);
        Vertex w = Vertex(1, &y);
        Vertex u = Vertex(2, &z);
        REQUIRE(v.OutDegree() == 0);
        REQUIRE(w.OutDegree() == 0);
        REQUIRE(u.OutDegree() == 0);

        REQUIRE(v.InDegree() == 0);
        REQUIRE(w.InDegree() == 0);
        REQUIRE(u.InDegree() == 0);
        
        // make u-->v connection
        v.AddInNeighbor(u.Id());
        u.AddOutNeighbor(v.Id());

        REQUIRE(v.InDegree() == 1);
        REQUIRE(u.OutDegree() == 1);

        REQUIRE(v.IsInNeighbor(u.Id()));
        REQUIRE(u.IsOutNeighbor(v.Id()));

        REQUIRE(!v.IsInNeighbor(w.Id()));
    }

    SECTION("HmmGraph objects can have vertices added to them") {
        std::string v = RandomNucleotides(RandomInt(5, 10));
        std::string w = RandomNucleotides(RandomInt(5, 10));
        std::string x = RandomNucleotides(RandomInt(5, 10));
        std::string y = RandomNucleotides(RandomInt(5, 10));
        std::string z = RandomNucleotides(RandomInt(5, 10));

        HmmGraph G = HmmGraph();

        int64_t n0 = G.AddVertex(&v);
        int64_t n1 = G.AddVertex(&w);
        int64_t n2 = G.AddVertex(&x);
        int64_t n3 = G.AddVertex(&y);
        int64_t n4 = G.AddVertex(&z);

        REQUIRE(G.K() == 5);
    }
    
    SECTION("HmmGraph objects are created, destroyed and sorted correctly") {
        std::string v = RandomNucleotides(RandomInt(5, 10));
        std::string w = RandomNucleotides(RandomInt(5, 10));
        std::string x = RandomNucleotides(RandomInt(5, 10));
        std::string y = RandomNucleotides(RandomInt(5, 10));
        std::string z = RandomNucleotides(RandomInt(5, 10));

        HmmGraph G = HmmGraph();

        int64_t n0 = G.AddVertex(&v);
        int64_t n1 = G.AddVertex(&w);
        int64_t n2 = G.AddVertex(&x);
        int64_t n3 = G.AddVertex(&y);
        int64_t n4 = G.AddVertex(&z);

        REQUIRE(!G.TestSort());
        
        // n0->n3->n2-v
        //   \         n1
        //     n4-----^
        G.AddArc(n0, n3);
        G.AddArc(n3, n2);
        G.AddArc(n2, n1);
        G.AddArc(n4, n1);
        G.AddArc(n0, n4);
        
        REQUIRE(!G.TestSort());

        REQUIRE_NOTHROW(G.TopologicalSort(true));
        REQUIRE(G.TestSort());
        REQUIRE(G.Vertices().size() == G.K());
        
        std::string a = RandomNucleotides(RandomInt(5, 10));

        // n0->n3->n2-v   (add a cycle)
        //   \  ^-----n1
        //     n4-----^
        G.AddArc(n1, n3);
        REQUIRE_THROWS_AS(G.TopologicalSort(true), ParcoursException);
    }
    
    SECTION("HmmGraph finds sources and sinks correctly") {
        std::string t = RandomNucleotides(RandomInt(5, 10));
        std::string v = RandomNucleotides(RandomInt(5, 10));
        std::string w = RandomNucleotides(RandomInt(5, 10));
        std::string x = RandomNucleotides(RandomInt(5, 10));
        std::string y = RandomNucleotides(RandomInt(5, 10));
        std::string z = RandomNucleotides(RandomInt(5, 10));

        HmmGraph G = HmmGraph();

        int64_t n0 = G.AddVertex(&v);
        int64_t n1 = G.AddVertex(&w);
        int64_t n2 = G.AddVertex(&x);
        int64_t n3 = G.AddVertex(&y);
        int64_t n4 = G.AddVertex(&z);
        int64_t n5 = G.AddVertex(&t);
        
        // n0\v     />n4
        //    n1->n2
        // n3-^      \>n5
        G.AddArc(n0, n1);
        G.AddArc(n3, n1);
        G.AddArc(n1, n2);
        G.AddArc(n2, n4);
        G.AddArc(n2, n5);

        std::set<int64_t> sources = G.Sources();
        std::set<int64_t> sinks = G.Sinks();

        REQUIRE(sources.count(n0) == 1);
        REQUIRE(sources.count(n3) == 1);
        REQUIRE(sinks.count(n4) == 1);
        REQUIRE(sinks.count(n5) == 1);
    }

    SECTION("HmmGraph finds all paths correctly") {
        std::string a = RandomNucleotides(RandomInt(5, 10));
        std::string b = RandomNucleotides(RandomInt(5, 10));
        std::string c = RandomNucleotides(RandomInt(5, 10));
        std::string d = RandomNucleotides(RandomInt(5, 10));
        std::string e = RandomNucleotides(RandomInt(5, 10));
        std::string f = RandomNucleotides(RandomInt(5, 10));
        std::string g = RandomNucleotides(RandomInt(5, 10));

        HmmGraph G = HmmGraph();

        int64_t n0 = G.AddVertex(&a);
        int64_t n1 = G.AddVertex(&b);
        int64_t n2 = G.AddVertex(&c);
        int64_t n3 = G.AddVertex(&d);
        int64_t n4 = G.AddVertex(&e);
        int64_t n5 = G.AddVertex(&f);
        int64_t n6 = G.AddVertex(&g);

        //   />n1    />n4
        // n0    x>n3    x>n6
        //   \>n2    \>n5
        G.AddArc(n0, n1);
        G.AddArc(n0, n2);
        G.AddArc(n1, n3);
        G.AddArc(n2, n3);
        G.AddArc(n3, n4);
        G.AddArc(n3, n5);
        G.AddArc(n4, n6);
        G.AddArc(n5, n6);

        G.TopologicalSort(false);
        REQUIRE(G.TestSort());

        REQUIRE(G.Sources().count(n0) == 1);
        REQUIRE(G.Sinks().count(n6) == 1);
        REQUIRE(G.Sources().size() == 1);
        REQUIRE(G.Sinks().size() == 1);
        
        std::vector<std::deque<int64_t>> paths = G.AllPaths();
        
        std::deque<int64_t> p0 {0, 1, 3, 4, 6};
        std::deque<int64_t> p1 {0, 1, 3, 5, 6};
        std::deque<int64_t> p2 {0, 2, 3, 4, 6};
        std::deque<int64_t> p3 {0, 2, 3, 5, 6};
        
        REQUIRE(paths.size() == 4);

        REQUIRE(std::find(paths.begin(), paths.end(), p0) != paths.end());
        REQUIRE(std::find(paths.begin(), paths.end(), p1) != paths.end());
        REQUIRE(std::find(paths.begin(), paths.end(), p2) != paths.end());
        REQUIRE(std::find(paths.begin(), paths.end(), p3) != paths.end());

        // check that rule-of-five functions work
        HmmGraph G2 = G;
        HmmGraph G3;
        G3 = G;

        auto compare_graphs = [] (HmmGraph& orig, HmmGraph& oth) {
            REQUIRE(orig.K() == oth.K());
            for (auto i : orig.Vertices()) {
                REQUIRE(orig.VertexGetter(i)->Sequence() == oth.VertexGetter(i)->Sequence());
                REQUIRE(orig.VertexGetter(i)->OutDegree() == oth.VertexGetter(i)->OutDegree());
                REQUIRE(orig.VertexGetter(i)->InDegree() == oth.VertexGetter(i)->InDegree());
                if (orig.OutNeighbors(i).size() > 0) {
                    REQUIRE(oth.OutNeighbors(i).size() > 0);
                    REQUIRE(oth.OutNeighbors(i) == orig.OutNeighbors(i));
                } else {
                    REQUIRE(oth.OutNeighbors(i).size() == 0);
                }
                if (orig.InNeighbors(i).size() > 0) {
                    REQUIRE(oth.InNeighbors(i).size() > 0);
                    REQUIRE(oth.InNeighbors(i) == orig.InNeighbors(i));
                } else {
                    REQUIRE(oth.InNeighbors(i).size() == 0);
                }
            }
            REQUIRE(orig.Sources() == oth.Sources());
            REQUIRE(orig.Sinks() == oth.Sinks());
            auto orig_paths = orig.AllPaths();
            auto oth_paths = oth.AllPaths();
            REQUIRE(orig_paths.size() == oth_paths.size());
            for (auto path : orig_paths) {
                REQUIRE(std::find(begin(oth_paths), end(oth_paths), path) != end(oth_paths));
            }
        };
        compare_graphs(G, G2);
        compare_graphs(G, G3);
        compare_graphs(G2, G3);
    }
}
