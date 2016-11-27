// HmmGraphAlignment tests
//
#include "test_helpers.h"
#include "hmm_graph.h"
#include "pairwise_aligner.h"

#ifdef _OPENMP
#include "omp.h"
#endif

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
    SECTION("Multipath alignment gets correct aligned pairs in simple example") {
        std::string s1 ("AG");
        std::string s2 ("TT");
        std::string s3 ("CG");

        HmmGraph G = HmmGraph();

        int64_t vid0 = G.AddVertex(&s1);
        int64_t vid1 = G.AddVertex(&s2);
        int64_t vid2 = G.AddVertex(&s3);

        G.AddArc(vid0, vid1);
        G.AddArc(vid1, vid2);
        G.AddArc(vid0, vid2);

        std::string read ("AGCG");

        AnchorPairs anchors = EmptyAnchors();

        AlignmentParameters p;
        p.expansion = 2;
        p.threshold = 0.2;
        p.ignore_gaps = false;

        FiveStateSymbolHmm hmm(SetNucleotideEmissionsToDefauts());

        G.Align<FiveStateSymbolHmm, fiveState>(read, anchors, p, hmm);
        
        auto pairs = G.PathAlignedPairs();
        
        REQUIRE(pairs.size() == G.NumberOfPaths());

        for (auto kv : pairs) {
            REQUIRE(kv.second.size() == 4);
            // for checking pairs, delete when annoying to look at
            //st_uglyf("Path %lld aligned pairs\n", kv.first);
            //for (auto p : kv.second) {
            //    st_uglyf("x: %lld vertex: %lld offset: %lld p: %f\n", 
            //            std::get<1>(p), std::get<2>(p).first, std::get<2>(p).second, std::get<0>(p));
            //}
        }

        auto aln_pairs_no_probs = [&pairs] () -> std::set<std::tuple<int64_t, int64_t, int64_t>> {
            std::set<std::tuple<int64_t, int64_t, int64_t>> set_of_pairs;
            for (auto kv : pairs) {
                for (auto p : kv.second) {
                    int64_t x      = std::get<1>(p);
                    int64_t vid    = std::get<2>(p).first;
                    int64_t offset = std::get<2>(p).second;
                    auto t = std::make_tuple(x, vid, offset);
                    set_of_pairs.insert(t);
                }
            }
            return set_of_pairs;
        }();

        REQUIRE(aln_pairs_no_probs.size() == 4);

    }
    SECTION("Paths have empty AlignedPairs where there aren't any") {
        std::string s1 ("AAAATTTT");
        std::string s2 ("GG");
        std::string s3 ("AAAATTTT");

        HmmGraph G = HmmGraph();

        int64_t vid0 = G.AddVertex(&s1);
        int64_t vid1 = G.AddVertex(&s2);
        int64_t vid2 = G.AddVertex(&s3);

        G.AddArc(vid0, vid1);
        G.AddArc(vid1, vid2);
        G.AddArc(vid0, vid2);

        std::string read ("CCCCCCCCCCCCCCCC");

        AnchorPairs anchors = EmptyAnchors();

        AlignmentParameters p;
        p.expansion = 4;
        p.threshold = 0.9;
        p.ignore_gaps = false;

        FiveStateSymbolHmm hmm(SetNucleotideEmissionsToDefauts());

        G.Align<FiveStateSymbolHmm, fiveState>(read, anchors, p, hmm);
        
        auto pairs = G.PathAlignedPairs();
        
        REQUIRE(pairs.size() == G.NumberOfPaths());

        for (auto kv : pairs) {
            REQUIRE(kv.second.size() == 0);
            // for checking pairs, delete when annoying to look at
            //st_uglyf("Path %lld aligned pairs\n", kv.first);
            //for (auto p : kv.second) {
            //    st_uglyf("x: %lld vertex: %lld offset: %lld p: %f\n", 
            //            std::get<1>(p), std::get<2>(p).first, std::get<2>(p).second, std::get<0>(p));
            //}
        }
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
        p.ignore_gaps = false;
        
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
        std::string read_string = StringFromSymbolString(path1_read);
        G.Align<FiveStateSymbolHmm, fiveState>(read_string, p, hmm);
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
        REQUIRE(G.MaxScorePath() == max_id);
    }
    SECTION("Correct path becomes most likely after reads are aligned to it") {
        int64_t incorrect = 0;
#ifdef _OPENMP
        int64_t test_cases = 40;
#else
        int64_t test_cases = 5;
#endif
        #pragma omp parallel for num_threads(4) reduction(+:incorrect) 
        for (int64_t test = 0; test < test_cases; test++) {
            // setup the graph and get the path sequences
            std::string a = RandomNucleotides(RandomInt(10, 20));
            std::string b = RandomNucleotides(RandomInt(5, 10));
            std::string c = RandomNucleotides(RandomInt(5, 10));
            std::string d = RandomNucleotides(RandomInt(10, 20));
            std::string e = RandomNucleotides(RandomInt(5, 10));
            std::string f = RandomNucleotides(RandomInt(5, 10));
            std::string g = RandomNucleotides(RandomInt(10, 20));

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
            p.ignore_gaps = false;
    
            AnchorPairs anchors = EmptyAnchors();

            int64_t choose_path = RandomInt(0, G.NumberOfPaths() - 1);

            // generate reads from one of the paths (simulate ~15X coverage)
            SymbolString path_sequence = G.PathSequences()[choose_path];
            std::string path_sequence_string = StringFromSymbolString(path_sequence);
            REQUIRE(path_sequence_string.size() == path_sequence.size());
            std::vector<SymbolString> reads;
            for (int64_t i = 0; i < 20; i++) {
                SymbolString path_read = [&] () {
                    std::string r = EvolveSequence(path_sequence_string);
                    SymbolString s = SymbolStringFromString(r);
                    //SymbolString s = SymbolStringFromString(path_sequence_string);
                    return s;
                }();
                reads.push_back(path_read);
            }
            
            // align them to the graph's paths
            G.Align<FiveStateSymbolHmm, fiveState>(reads, anchors, p, hmm);
            
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
            if (most_probable_path != choose_path) {
                //st_uglyf("INCORRECT most probable path\ngot %lld expected %lld\n", most_probable_path, choose_path);
                //st_uglyf("exp: %s\nobs: %s\n", path_sequence_string.c_str(), 
                //        StringFromSymbolString(G.PathSequences()[most_probable_path]).c_str());
                //st_uglyf("exp: %f\nobs: %f\n", G.PathScores()[choose_path], G.PathScores()[most_probable_path]);
                //#pragma omp critical (incorrect) 
                //{
                incorrect++;
                //}
            }
        }
        double accuracy = 100.0 - (((double )incorrect / (double )test_cases) * 100);
        //st_uglyf("incorrect %lld of %lld, accuracy %f\n", incorrect, test_cases, accuracy);
        REQUIRE(accuracy >= 95.0);
    }
}
