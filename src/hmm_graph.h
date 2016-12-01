//
// 
//

#ifndef PARCOURS_HMM_GRAPH_H
#define PARCOURS_HMM_GRAPH_H

#include "stl_includes.h"
#include "symbol_string.h"
#include "common.h"
#include "pairwise_aligner.h"
#include "vertex.h"

// intermediate graph structure
class HmmGraph {
public:
    /*
     * Graph manipulation methods
     */
    HmmGraph();
    
    ~HmmGraph(); 
   
    HmmGraph(HmmGraph& other); 

    HmmGraph& operator = (HmmGraph& other);

    HmmGraph(HmmGraph&& other) = default;
    
    HmmGraph& operator = (HmmGraph&& other) = default;

    int64_t AddVertex(const std::string *seq);

    bool ContainsVertex(int64_t i) const;

    bool ContainsVertex(Vertex *v) const;

    void AddArc(int64_t, int64_t);
    
    void SetNextVertexId(int64_t i);

    int64_t K();

    Vertex *VertexGetter(int64_t i);  // caution: returns raw pointer to vertex

    unsigned long VertexOutDegree(int64_t i);

    unsigned long VertexInDegree(int64_t i);

    const std::string *VertexSequence(int64_t i);
    
    const std::set<int64_t>& OutNeighbors(int64_t i);

    const std::set<int64_t>& InNeighbors(int64_t i);

    const std::vector<int64_t>& Vertices() const;

    void TopologicalSort(bool test=false);

    bool TestSort();

    bool IsSorted();

    std::set<int64_t> Sources();

    std::set<int64_t> Sinks();
    
    std::unordered_map<int64_t, std::deque<int64_t>> PathMap();

    std::unordered_map<int64_t, SymbolString> PathSequences();

    // returns just the vertex_paths (the walks) from each source to sink 
    // removes the pathIDs
    std::vector<std::deque<int64_t>> AllPaths();

    int64_t NumberOfPaths();

    std::unordered_map<int64_t, double> PathScores(bool normalize=true);
    
    int64_t MaxScorePath();
    /*
     * Alignment Methods
     */
    // returns the map of pathID to aligned pairs
    std::unordered_map<int64_t, GraphAlignedPairs>& PathAlignedPairs();
    // Aligns a sequence (SymbolString) to the path_sequences in the graph
    // modifies: path_scores, path_aligned_pairs (iff get_pairs==True)
    // calls: PairwiseAlignment::Score() PairwiseAlignment::AlignedPairsGetter()
    template<class Hmm, size_t sn>
    void Align(SymbolString& S, AnchorPairs& anchors, AlignmentParameters& p, Hmm& hmm, 
               bool get_pairs=true, bool ragged_end=false);
    // with string
    template<class Hmm, size_t sn>
    void Align(std::string& S, AnchorPairs& anchors, AlignmentParameters& p, Hmm& hmm, 
               bool get_pairs=true, bool ragged_end=false);
    // if you don't want to use anchors, good for small alignments
    template<class Hmm, size_t sn>
    void Align(std::string& S, AlignmentParameters& p, Hmm& hmm, bool get_pairs=true, bool ragged_end=false);
    // calls the above function for each SymbolString in the vector, continually updates scores and 
    // adds to path_aligned_pairs
    template<class Hmm, size_t sn>
    void Align(std::vector<SymbolString>& vS, AnchorPairs& anchors, AlignmentParameters& p, Hmm& hmm, 
               bool get_pairs=true, bool ragged_end=false);

    // Call the above template methods with implemented models for convenience
    void AlignWithFiveStateSymbolHmm(std::string& S, AnchorPairs& anchors, AlignmentParameters& p, 
                                     bool get_pairs=true, bool ragged_end=false);

    void AlignWithFiveStateSymbolHmm(std::string& S, AlignmentParameters& p, 
                                     bool get_pairs=true, bool ragged_end=false);


private:
    // containers for graph representation
    std::unordered_map<int64_t, std::set<int64_t>> adjacentcy_list;
    std::unordered_map<int64_t, std::unique_ptr<Vertex>> vertex_map;  // for looking up vertices
    std::vector<int64_t> vertex_list;  // for ordering

    // containers for alignment and keeping track of sequences through the graph.
    // All are maps where there is a unique pathID the identifies the vertex_path 
    // (the sequences of vertices defining the path from source to sink).
    std::unordered_map<int64_t, std::deque<int64_t>> paths;  // (pathID, vertex_path)
    // the path sequence are the sequence generated by concatenating the vertex sequences 
    // in the path
    std::unordered_map<int64_t, SymbolString> path_sequences;  // (pathID, SymbolString)
    // to map aligned pairs back to vertices, we maintain a mapping of pathID to pairs
    // where each pair is (vertex, offset) and the position of the pair in the vector
    // is the position that pair corresponds to in the path sequence, the vectros must
    // be the same length as the sequences.
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, int64_t>>> sequence_to_vertex;
    // contains path scores as scored by pairwise aligner, contains the score for the entire 
    // path sequence, given all the reads that have been aligned
    std::unordered_map<int64_t, double> path_scores;  // (pathID, score)
    // contains the aligned pairs for each path, should accumulate as more reads are aligned 
    std::unordered_map<int64_t, GraphAlignedPairs> path_aligned_pairs;

    // counters and flags
    //
    int64_t nVertices;
    int64_t nArcs;
    int64_t next_vertex_id;
    int64_t nPaths;
    int64_t most_probable_path;
    bool sorted = false;
    bool initialized_paths = false;
    bool normalized_path_scores = false;

    // Internal Methods
    //
    // fills the `paths`, `path_sequences`, and `sequence_to_vertex` containers throws an exception 
    // iff doule check is true and results from the called methods don't check ok, will automatically
    // check if `initialized_paths` is already true
    // modifies: paths, path_sequences, sequence_to_vertex
    // calls: find_paths(), extract_sequences()
    void initializePaths(bool double_check=true);  
    void extractSequences();
    // test_sort is passed to TopologicalSort, if `false` the find_paths will not throw an exceptiion
    // modifies: paths, nPaths
    // calls:    TopologicalSort(), Sources(), Sinks(), Vertex->OutNeighbors(),
    void findPaths(bool test_sort=true);
    void normalizePathScores();
    GraphAlignedPairs translatePairwiseAlignedPairs(AlignedPairs& pairs, int64_t pathId);
    // TODO 
    // - Sort paths by score
    // - Flattened GraphAlignedPairs
    // - GraphAlignedPairs to CIGAR?
    void copyGraph(HmmGraph& orig, HmmGraph& other);
    void clearGraph();
};

#endif
