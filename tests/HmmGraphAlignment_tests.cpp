// HmmGraphAlignment tests
//

#include "catch.hpp"
#include "test_helpers.h"
#include "hmm_graph.h"
#include "pairwise_aligner.h"

TEST_CASE("Multipath Alignment test", "[alignment]") {
    SECTION("Graph extracts all path sequences correctly") {
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

        auto path_sequences_hash = G.PathSequences();
    
        auto path1_seq = a + b + d + e + g;
        auto path2_seq = a + b + d + f + g;
        auto path3_seq = a + c + d + e + g;
        auto path4_seq = a + c + d + f + g;

        SymbolString s1 = SymbolStringFromString(path1_seq);
        SymbolString s2 = SymbolStringFromString(path2_seq);
        SymbolString s3 = SymbolStringFromString(path3_seq);
        SymbolString s4 = SymbolStringFromString(path4_seq);
        
        auto path_sequences = [&] () -> std::vector<SymbolString> {
            std::vector<SymbolString> S;
            for (auto kv : path_sequences_hash) {
                S.push_back(kv.second);
            }
            return S;
        }();

        REQUIRE(std::find(begin(path_sequences), end(path_sequences), s1) != end(path_sequences));
        REQUIRE(std::find(begin(path_sequences), end(path_sequences), s2) != end(path_sequences));
        REQUIRE(std::find(begin(path_sequences), end(path_sequences), s3) != end(path_sequences));
        REQUIRE(std::find(begin(path_sequences), end(path_sequences), s4) != end(path_sequences));
    }
    SECTION("Sequences align to graph paths correctly when read matches exactly") {
        // setup the graph and get the path sequences
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

        std::unordered_map<int64_t, SymbolString> path_sequences = G.PathSequences();
        
        // Generate a "read" that matches one of the paths exactly
        SymbolString path1_read = [&] () {
            std::string read_string = (a + b + d + e + g);
            return SymbolStringFromString(read_string);
        }();        

        FiveStateSymbolHmm hmm;
        hmm.InitializeEmissions(SetNucleotideEmissionsToDefauts());

        std::unordered_map<int64_t, double> scores;
        
        AlignmentParameters p;
        p.expansion = 2;
        p.threshold = 0.1;
        
        AnchorPairs anchors = EmptyAnchors();

        for (auto& kv : path_sequences) {
            PairwiseAlignment<FiveStateSymbolHmm, fiveState> aln(hmm, path1_read, kv.second, anchors, p);
            scores[kv.first] = aln.Score(false);
        }

        auto max_id = [&] () -> int64_t {
            double max_score = LOG_ZERO;
            int mId = -1;
            for (auto kv : scores) {
                if (kv.second > max_score) {
                    max_score = kv.second;
                    mId = kv.first;
                } else {
                    continue;
                }
            }
            REQUIRE(mId >= 0);
            return mId;
        }();
        REQUIRE(max_id == 0);

        // test internal alignment routine
        G.Align<FiveStateSymbolHmm, fiveState>(path1_read, anchors, p, hmm, false);
        REQUIRE(G.PathScores(false) == scores);

        // test normalization 
        auto total_prob = [scores] () {
            double t = 0.0;
            for (auto kv : scores) {
                t += kv.second;
            }
            return t;
        }();

        for (auto& kv : scores) {
            kv.second /= total_prob;
        }

        REQUIRE(G.PathScores(true) == scores);
    }
    SECTION("Correct path becomes most likely after reads are aligned to it") {
        for (int64_t test = 0; test < 100; test++) {
            // setup the graph and get the path sequences
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

            FiveStateSymbolHmm hmm(SetNucleotideEmissionsToDefauts());

            AlignmentParameters p;
            p.expansion = 6;
            p.threshold = 0.2;
    
            AnchorPairs anchors = EmptyAnchors();

            int64_t choose_path = RandomInt(0, G.NumberOfPaths() - 1);

            // generate 30 reads from one of the paths
            SymbolString path_sequence = G.PathSequences()[choose_path];
            std::string path_sequence_string = StringFromSymbolString(path_sequence);
            REQUIRE(path_sequence_string.size() == path_sequence.size());
            std::vector<SymbolString> reads;
            for (int64_t i = 0; i < 10; i++) {
                SymbolString path_read = [&] () {
                    std::string r = EvolveSequence(path_sequence_string);
                    SymbolString s = SymbolStringFromString(r);
                    return s;
                }();
                reads.push_back(path_read);
            }
            
            // align them to the graph's paths
            G.Align<FiveStateSymbolHmm, fiveState>(reads, anchors, p, hmm, false);
            
            auto most_probable_path = [&] () -> int64_t {
                double max_score = LOG_ZERO;
                int mId = -1;
                for (auto kv : G.PathScores(true)) {
                    if (kv.second > max_score) {
                        max_score = kv.second;
                        mId = kv.first;
                    } else {
                        continue;
                    }
                }
                REQUIRE(mId >= 0);
                return mId;
            }();
            REQUIRE(most_probable_path == choose_path);
        }
    }
}
