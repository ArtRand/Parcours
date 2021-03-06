// hmm_graph.cpp

#include "hmm_graph.h"

// Public Methods

HmmGraph::HmmGraph(): nVertices(0), nArcs(0), next_vertex_id(0), nPaths(-1), most_probable_path(-1) {  }

HmmGraph::~HmmGraph() {  }

HmmGraph::HmmGraph(HmmGraph& other) {
    if (&other != this) copyGraph(*this, other);
}

HmmGraph& HmmGraph::operator = (HmmGraph& other) {
    if (&other != this) {
        copyGraph(*this, other);
    }
    return *this;
}

int64_t HmmGraph::AddVertex(const std::string *seq) {
    // get the next vertex id, a monotonically increasing int
    int64_t vId = next_vertex_id;

    //Vertex *v = new Vertex(vId, seq);
    std::unique_ptr<Vertex> v(new Vertex(vId, seq));

    // check, TODO remove after tested
    assert(v->Id() == vId);

    // updates to state of graph
    //vertex_map[vId] = v;
    vertex_map[vId] = std::move(v);
    vertex_list.push_back(vId);
    nVertices += 1;
    next_vertex_id += 1;
    sorted = false;

    return vId;
}

bool HmmGraph::ContainsVertex(int64_t i) const {
    return vertex_map.count(i) != 0;
}

bool HmmGraph::ContainsVertex(Vertex *v) const {
    return ContainsVertex(v->Id());
}

void HmmGraph::AddArc(int64_t fromId, int64_t toId) {
    if (!ContainsVertex(fromId)) {
        throw ParcoursException("[HmmGraph::AddArc]: Missing from vertex %" PRIi64 "", fromId);
    }

    if (!ContainsVertex(toId)) {
        throw ParcoursException("[HmmGraph::AddArc]: Missing to vertex %" PRIi64 "", toId);
    }
    // update adjacentcy list
    adjacentcy_list[fromId].insert(toId);
    // update the to-vertex
    vertex_map[toId]->AddInNeighbor(fromId);
    // update the from-vertex
    vertex_map[fromId]->AddOutNeighbor(toId);
    // update counter and sort status
    nArcs += 1;
    sorted = false;
    return;
}

void HmmGraph::SetNextVertexId(int64_t i) {
    if (nVertices > 0) {
        throw ParcoursException("[HmmGraph::SetNextVertexId] Cannot set next vertex ID when "
                                "a graph has already been setup, may screw up ordering, "
                                " requested %" PRIi64 ", currently have %" PRIi64 "vertices", i, nVertices);
    }
    assert(next_vertex_id == 0);
    next_vertex_id = i;
}

int64_t HmmGraph::K() { return nVertices; }

Vertex *HmmGraph::VertexGetter(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::VertexGetter]: Graph doesn't contain vertex %" PRIi64 "", i);
    }
    
    return vertex_map[i].get();

}

unsigned long HmmGraph::VertexOutDegree(int64_t i) {
    //return VertexGetter(i)->OutDegree();
    return vertex_map[i]->OutDegree();
}

unsigned long HmmGraph::VertexInDegree(int64_t i) {
    //return VertexGetter(i)->OutDegree();
    return vertex_map[i]->OutDegree();
}

const std::string *HmmGraph::VertexSequence(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::VertexSequence]: Graph does not contain vertex %" PRIi64 "", i);
    }
    return vertex_map[i]->Sequence();
}

const std::set<int64_t>& HmmGraph::InNeighbors(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::InNeighbors]: Graph does not contain vertex %" PRIi64 "", i);
    }

    //return VertexGetter(i)->InNeighbors();
    return vertex_map[i]->InNeighbors();
}

const std::set<int64_t>& HmmGraph::OutNeighbors(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::OutNeighbors]: Graph does not contain vertex %" PRIi64 "", i);
    }
    
    return vertex_map[i]->OutNeighbors();
}

const std::vector<int64_t>& HmmGraph::Vertices() const {
    return vertex_list;
}

static inline void dfs(HmmGraph *G, int64_t s, 
                       std::unordered_map<int64_t, bool>& discovered, 
                       std::vector<int64_t>& finished, 
                       std::unordered_map<int64_t, bool>& onStack) {
    discovered[s] = true;
    onStack[s] = true;
    for (int64_t id : G->OutNeighbors(s)) {
        // default of map is false
        if (!discovered[id]) {
            dfs(G, id, discovered, finished, onStack);
        }
        if (onStack[id]) {
            throw ParcoursException("[HmmGraph::dfs(internal)]: Cycle found.");
        }
    }
    onStack[s] = false;
    finished.push_back(s);
}

// standard DFS topo sort
void HmmGraph::TopologicalSort(bool test) {
    // map of vertices we've visited
    std::unordered_map<int64_t, bool> discovered;
    // vector of vertices that have been completely explored, contains the reverse of the topo sort
    std::vector<int64_t> finished;
    // vertices in the recursion
    std::unordered_map<int64_t, bool> onStack;

    auto isFinished = [&finished](int64_t q) -> bool {
        std::vector<int64_t>::iterator it;
        it = std::find(finished.begin(), finished.end(), q);
        return it != finished.end();
    };

    for (int64_t i : vertex_list) {
        if (!isFinished(i)) {
            dfs(this, i, discovered, finished, onStack);
        }
    }

    std::reverse(finished.begin(), finished.end());

    assert(finished.size() == static_cast<uint64_t>(K()));

    for (uint64_t i = 0; i < vertex_list.size(); i++) {
        vertex_list[i] = finished[i];
    }

    sorted = true;

    if (test) {
        bool check = TestSort();
        if (!check) {
            throw ParcoursException("[HmmGraph::TopologicalSort]: Did not pass sorting test");
        }
    }

    return;
}

bool HmmGraph::TestSort() {
    // decide that an empty graph cannot be sorted
    if (K() == 0) return false; // an empty graph is not sorted
    // if there aren't any arcs, but there is more than one vertex, it
    // cannot be sorted
    if (adjacentcy_list.size() == 0 && K() > 1) return false;

    // container to keep track of vertices we've visited
    std::set<int64_t> seen;
    
    // loop over the vertices
    for (int64_t i : vertex_list) {
        for (int64_t id : InNeighbors(i)) {
            bool check = seen.count(id) > 0;
            if (!check) {
                return false;
            }
        }
        seen.insert(i);
    }
    return true;
}

bool HmmGraph::IsSorted() { return sorted; }

std::set<int64_t> HmmGraph::Sources() {
    std::set<int64_t> sources;
    auto is_source = [&] (int64_t id) {
        if (vertex_map[id]->InDegree() == 0) sources.insert(id);
    };
    std::for_each(begin(vertex_list), end(vertex_list), is_source);
    return sources;
}

std::set<int64_t> HmmGraph::Sinks() {
    std::set<int64_t> sinks;
    auto is_sink = [&] (int64_t id) {
        if (vertex_map[id]->OutDegree() == 0) sinks.insert(id);
    };
    std::for_each(begin(vertex_list), end(vertex_list), is_sink);
    return sinks;
}

std::unordered_map<int64_t, std::deque<int64_t>> HmmGraph::PathMap() { 
    if (!initialized_paths) initializePaths(true);
    return paths; 
}

std::unordered_map<int64_t, SymbolString> HmmGraph::PathSequences() {
    if (!initialized_paths) initializePaths(true);
    return path_sequences;
}

std::vector<std::deque<int64_t>> HmmGraph::AllPaths() {
    if (!initialized_paths) initializePaths(true);
    
    std::vector<std::deque<int64_t>> vertex_paths;

    for (auto kv : paths) {
        vertex_paths.push_back(kv.second);
    }

    return vertex_paths;
}

int64_t HmmGraph::NumberOfPaths() {
    if (!initialized_paths) initializePaths();
    return nPaths;
}

std::unordered_map<int64_t, double> HmmGraph::PathScores(bool normalize) { 
    if (!normalized_path_scores && normalize) normalizePathScores();
    return path_scores; 
}

int64_t HmmGraph::MaxScorePath() { 
    if (!normalized_path_scores) normalizePathScores();
    return most_probable_path; 
}

std::unordered_map<int64_t, GraphAlignedPairs>& HmmGraph::PathAlignedPairs() { 
    return path_aligned_pairs; 
}

template<class Hmm, size_t sn>
void HmmGraph::Align(SymbolString& S, AnchorPairs& anchors, AlignmentParameters& p,  Hmm& hmm, 
                     bool get_aligned_pairs, bool ragged_end) {
    if (!initialized_paths) initializePaths(true);
    
    if ((path_sequences.size() == 0) || (S.size() == 0)) return;

    for (auto& kv : path_sequences) {
        PairwiseAlignment<Hmm, sn> aln(hmm, S, kv.second, anchors, p, ragged_end);
        // perform the alignment and get the aligned pairs, if requested
        if (get_aligned_pairs) {
            AlignedPairs aln_pairs = aln.AlignedPairsGetter();
            // translate the aligned pairs to reflect the graph coordinates
            GraphAlignedPairs g_pairs = translatePairwiseAlignedPairs(aln_pairs, kv.first);
            assert(g_pairs.size() == aln_pairs.size());
            // add these aligned pairs to the rest for this path
            path_aligned_pairs[kv.first].insert(end(path_aligned_pairs[kv.first]), begin(g_pairs), end(g_pairs));
        }
        path_scores[kv.first] += aln.Score(p.ignore_gaps);
    }
}

template<class Hmm, size_t sn>
void HmmGraph::Align(std::string& sx, AnchorPairs& anchors, AlignmentParameters& p, Hmm& hmm, 
                     bool get_pairs, bool ragged_end) {
    auto S = SymbolStringFromString(sx);
    Align<Hmm, sn>(S, anchors, p, hmm, get_pairs, ragged_end);
}

template<class Hmm, size_t sn>
void HmmGraph::Align(std::string& sx, AlignmentParameters& p, Hmm& hmm, bool get_pairs, bool ragged_end) {
    auto S = SymbolStringFromString(sx);
    AnchorPairs anchors = EmptyAnchors();
    Align<Hmm, sn>(S, anchors, p, hmm, get_pairs, ragged_end);
}

template<class Hmm, size_t sn>
void HmmGraph::Align(std::vector<SymbolString>& vS, AnchorPairs& anchors, AlignmentParameters& p, Hmm& hmm, 
                     bool get_pairs, bool ragged_end) {
    for (auto s : vS) {
        Align<Hmm, sn>(s, anchors, p, hmm, get_pairs, ragged_end);
    }
}

template void HmmGraph::Align<FiveStateSymbolHmm, fiveState>(SymbolString& S, 
                                                             AnchorPairs& anchors, 
                                                             AlignmentParameters& p, 
                                                             FiveStateSymbolHmm& hmm, 
                                                             bool get_aligned_pairs, 
                                                             bool ragged_end);

template void HmmGraph::Align<FiveStateSymbolHmm, fiveState>(std::string& sx, 
                                                             AnchorPairs& anchors, 
                                                             AlignmentParameters& p, 
                                                             FiveStateSymbolHmm& hmm, 
                                                             bool get_aligned_pairs, 
                                                             bool ragged_end);

template void HmmGraph::Align<FiveStateSymbolHmm, fiveState>(std::string& sx, 
                                                             AlignmentParameters& p, 
                                                             FiveStateSymbolHmm& hmm, 
                                                             bool get_aligned_pairs, 
                                                             bool ragged_end);

template void HmmGraph::Align<FiveStateSymbolHmm, fiveState>(std::vector<SymbolString>& S, 
                                                             AnchorPairs& anchors, 
                                                             AlignmentParameters& p, 
                                                             FiveStateSymbolHmm& hmm, 
                                                             bool get_aligned_pairs, 
                                                             bool ragged_end);

void HmmGraph::AlignWithFiveStateSymbolHmm(std::string& S, AnchorPairs& anchors, AlignmentParameters& p, 
                                           bool get_pairs, bool ragged_end) {
    FiveStateSymbolHmm hmm(SetNucleotideEmissionsToDefauts());
    Align<FiveStateSymbolHmm, fiveState>(S, anchors, p, hmm, get_pairs, ragged_end);
}

void HmmGraph::AlignWithFiveStateSymbolHmm(std::string& S, AlignmentParameters& p, 
                                           bool get_pairs, bool ragged_end) {
    FiveStateSymbolHmm hmm(SetNucleotideEmissionsToDefauts());
    Align<FiveStateSymbolHmm, fiveState>(S, p, hmm, get_pairs, ragged_end);
}


/*
* Internal Methods
*/

void HmmGraph::initializePaths(bool double_check) {
    auto check = [this] () -> bool {
        assert(initialized_paths);
        if (nPaths < 0) return false;                                        // not initialized
        if (nPaths != static_cast<int>(paths.size())) return false;          // missing paths or bad accounting
        if (paths.size() != path_sequences.size()) return false;             // don't have correct # sequences
        for (auto& kv : path_sequences) {
            int64_t pId    = kv.first;
            SymbolString s = kv.second;
            if (s.size() != sequence_to_vertex[pId].size()) return false;    // don't have vertex mapping correct
            else continue;
        }
        return true;
    };

    if (initialized_paths) {
        bool ok = check();
        if (!ok) throw ParcoursException("[HmmGraph::initializePaths] check not OK\n");
        else return;
    }
    
    try {
        findPaths(true);
        extractSequences();
        initialized_paths = true;
    } catch (ParcoursException& e) {
        std::cerr << e.what() << "\n";
    }
    
    if (double_check) { 
        try {
            assert(initialized_paths);
            initializePaths(false);
        } catch (ParcoursException& e) {
            std::cerr << e.what() << "\n";
        }
    }
}

void HmmGraph::extractSequences() {
    if (paths.empty()) return;
    for (auto& kv : paths) {  // (pathID, vertex_path)
        SymbolString S = [&] () {
            SymbolString s;
            // loop over the vertices, and convert the sequences into SymbolStrings, 
            // keep track of the offset and vertexID to map the SymbolString sequence
            // back to the vertex_sequence base it came from 
            for (int64_t vId : kv.second) {
                int64_t offset = 0;
                for (char b : *(VertexSequence(vId))) {
                    s.push_back(CharToSymbol(b));
                    sequence_to_vertex[kv.first].emplace_back(vId, offset++);
                }
            }
            return s;
        }();
        path_sequences[kv.first] = S;
        if (S.size() != sequence_to_vertex[kv.first].size()) throw ParcoursException(""
                "[HmmGraph::extractSequences] sequence %" PRIi64 " length not equal to vertex_to_sequence "
                "pairs length, expected %lu == %lu\n", kv.first, S.size(), sequence_to_vertex[kv.first].size());
    }
    
    if (path_sequences.size() != paths.size()) throw ParcoursException("[HmmGraph::extractSequences]"
            " error collecting path sequences");

    return;
}

void HmmGraph::findPaths(bool test_sort) {
    if (!paths.empty()) {
        st_uglyf("[HmmGraph::findPaths] clearing non-empty paths");
        paths.clear();
        nPaths = -1;
    }

    // we're going to find all the paths from each source to each sink
    std::set<int64_t> sources = Sources();
    std::set<int64_t> sinks = Sinks();
    
    TopologicalSort(test_sort);
    
    // reverse the vertex list, we're going to go from sink to source
    std::vector<int64_t> reverse_node_list = vertex_list;
    std::reverse(reverse_node_list.begin(), reverse_node_list.end());
        
    // TODO delete all this crap once I finish the multi-source/sink unittest
    //st_uglyf("forward node list: ");
    //for (auto i : vertex_list) std::cout << i << ", ";
    //st_uglyf("\n");
    
    //st_uglyf("reversed node list: ");
    //for (auto i : reverse_node_list) std::cout << i << ", ";
    //st_uglyf("\n");

    for (int64_t sink : sinks) {
        for (int64_t source : sources) {
            //st_uglyf("Performing dp for sink %lld source %lld\n", sink, source);

            // make a map to hold the dynamic programming intermediates
            std::map<int64_t, std::vector<std::deque<int64_t>>> path_hash;
            
            // base case, for the sink
            std::deque<int64_t> d { sink };
            path_hash[sink].push_back(d);

            // walk back to the source
            for (uint64_t i = 1; i < reverse_node_list.size(); i++) {
                // get the vertex id we're at
                int64_t vId = reverse_node_list.at(i);
                //st_uglyf("checking %lld\n", vId);
                // loop over the out-neighbors to this vertex, get their
                // paths, add this vertex id to them and add those paths to
                // this vertex's paths
                for (int64_t out_neighbor_id : vertex_map[vId]->OutNeighbors()) {
                    for (auto p : path_hash[out_neighbor_id]) {
                        // copy the path from the neighbor, add this vertex to 
                        // the start of it, and update the map
                        std::deque<int64_t> new_path(p);
                        new_path.push_front(vId);
                        path_hash[vId].push_back(new_path);
                    }
                }
                if (vId == source) {
                    for (auto& p : path_hash[source]) {
                        if (nPaths < 0) nPaths = 0;
                        paths[nPaths] = p;
                        nPaths++;
                    }
                    //st_uglyf("finished: \n");
                    //for (auto p : paths) {
                    //    for (auto v : p) std::cout << v << ", ";
                    //    std::cout << std::endl;
                    //}
                }
            }
        }
    }
    //assert(nPaths == paths.size());
}

void HmmGraph::normalizePathScores() {
    if (path_scores.empty()) throw ParcoursException("[HmmGraph::normalized_path_scores] path scores is empty");
    if (normalized_path_scores) std::cerr << "[HmmGraph::normalized_path_scores] WARNING normalizing scores"
                                             " that have already been normalized?" << std::endl;
    auto total_path_prob = [this] () -> double {
        double total = 0.0;
        for (auto kv : path_scores) {
            total += kv.second;
        }
        return total;
    }();

    double max_prob_score = LOG_ZERO;

    for (auto& kv : path_scores) {
        kv.second /= total_path_prob;
        if (kv.second > max_prob_score) {
            max_prob_score = kv.second;
            most_probable_path = kv.first;
        }
    }

    assert(most_probable_path >= 0);

    normalized_path_scores = true;
}

GraphAlignedPairs HmmGraph::translatePairwiseAlignedPairs(AlignedPairs& pairs, int64_t pathId) {
    // loop over the aligned pairs and replace the y-corrdinate with the vertex offset pair 
    GraphAlignedPairs g_pairs;
    g_pairs.reserve(pairs.size());
    auto translate = [this, pathId] (AlignedPair& aln_pair) -> GraphAlignedPair {
        double p     = std::get<0>(aln_pair);
        int64_t x    = std::get<1>(aln_pair);
        int64_t y    = std::get<2>(aln_pair);
        auto v_coord = sequence_to_vertex[pathId].at(y);
        return std::make_tuple(p, x, v_coord);
    };

    std::transform(begin(pairs), end(pairs), std::back_inserter(g_pairs), translate);

    return g_pairs;
}


void HmmGraph::clearGraph() {
    if (vertex_list.size() > 0 || nVertices > 0 || nArcs > 0) {
        std::cerr << "[HmmGraph::clearGraph] WARNING: clearning an initialized graph\n";
    }
    vertex_list.clear();
    adjacentcy_list.clear();
    vertex_map.clear();
    paths.clear();
    nVertices = 0;
    nArcs = 0;
    nPaths = -1;
    next_vertex_id = 0;
    sorted = false;
    initialized_paths = false;
}

void HmmGraph::copyGraph(HmmGraph& newer, HmmGraph& other) {
    newer.clearGraph();
    // determine the range of the other graph
    int64_t min_vertex_id = *(std::min_element(begin(other.Vertices()), end(other.Vertices())));
    int64_t max_vertex_id = *(std::max_element(begin(other.Vertices()), end(other.Vertices())));
    // move the vertex counter (next_vertex_id) to the minimum 
    newer.SetNextVertexId(min_vertex_id);
    // now add each vertex to this graph (+1 here because we want to include the final vertex)
    for (int64_t i = min_vertex_id; i < max_vertex_id + 1; i++) {
        try {
            newer.AddVertex(other.VertexSequence(i));
        } catch (ParcoursException& e) {
            std::cerr << e.what() << "\n";
        }
    }
    // now add the arcs
    for (int64_t id : newer.Vertices()) {
        for (int64_t n : other.OutNeighbors(id)) {
            if (newer.VertexSequence(id) != other.VertexSequence(id)) {
                throw ParcoursException("[HmmGraph copy constructor]: vertex %" PRIi64 "'s sequence does not"
                                        "match expected", id);
            }
            newer.AddArc(id, n);
        }
    } 
}


