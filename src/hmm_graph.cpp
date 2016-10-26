// hmm_graph.cpp

#include "hmm_graph.h"

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

bool HmmGraph::ContainsVertex(int64_t i) {
    return vertex_map.count(i) != 0;
}

bool HmmGraph::ContainsVertex(Vertex *v) {
    return ContainsVertex(v->Id());
}

void HmmGraph::AddArc(int64_t fromId, int64_t toId) {
    if (!ContainsVertex(fromId)) {
        throw ParcoursException("[HmmGraph::AddArc]: Missing from vertex %lld", fromId);
    }

    if (!ContainsVertex(toId)) {
        throw ParcoursException("[HmmGraph::AddArc]: Missing to vertex %lld", toId);
    }
    // update adjacentcy list
    adjacentcy_list[fromId].insert(toId);
    // update the to-vertex
    vertex_map[toId]->AddInNeighbor(fromId);
    // update the from-vertex
    vertex_map[fromId]->AddOutNeighbor(toId);
    // update sort status
    nArcs += 1;
    sorted = false;
    return;
}

int64_t HmmGraph::K() { return nVertices; }

Vertex *HmmGraph::VertexGetter(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::VertexGetter]: Graph doesn't contain vertex %lld", i);
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

const std::set<int64_t>& HmmGraph::InNeighbors(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::InNeighbors]: Graph does not contain vertex %lld", i);
    }

    //return VertexGetter(i)->InNeighbors();
    return vertex_map[i]->InNeighbors();
}

const std::set<int64_t>& HmmGraph::OutNeighbors(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::OutNeighbors]: Graph does not contain vertex %lld", i);
    }
    
    return vertex_map[i]->OutNeighbors();
}

const std::vector<int64_t>& HmmGraph::Vertices() {
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
    if (K() == 0) {
        return false;
    }
    // if there aren't any arcs, cannot be sorted
    if (adjacentcy_list.size() == 0) {
        return false;
    }

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
    for (int64_t id : vertex_list) {
        if (vertex_map[id]->InDegree() == 0) {
            sources.insert(id);
        } 
    }
    return sources;
}

std::set<int64_t> HmmGraph::Sinks() {
    std::set<int64_t> sinks;
    for (int64_t id : vertex_list) {
        if (vertex_map[id]->OutDegree() == 0) {
            sinks.insert(id);
        }
    }
    return sinks;
}

std::vector<std::deque<int64_t>> HmmGraph::AllPaths() {
    if (!initialized_paths) {
        find_paths();
    }
    return paths;
}

void HmmGraph::find_paths() {
    // we're going to find all the paths from each source to each sink
    std::set<int64_t> sources = Sources();
    std::set<int64_t> sinks = Sinks();
    
    // preliminary checks
    TopologicalSort(true);
    
    // reverse the vertex list, we're going to go from sink to source
    std::vector<int64_t> reverse_node_list = vertex_list;
    std::reverse(reverse_node_list.begin(), reverse_node_list.end());

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
            for (int64_t i = 1; i < reverse_node_list.size(); i++) {
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
                    paths.insert(paths.end(), path_hash[source].begin(), path_hash[source].end());
                    //st_uglyf("finished: \n");
                    //for (auto p : paths) {
                    //    for (auto v : p) std::cout << v << ", ";
                    //    std::cout << std::endl;
                    //}
                }
            }
        }
    }
    initialized_paths = true;
}
