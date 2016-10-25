// hmm_graph.cpp

#include "hmm_graph.h"

int64_t HmmGraph::AddVertex(std::string *seq) {
    // get the next vertex id, a monotonically increasing int
    int64_t vId = next_vertex_id;

    Vertex *v = new Vertex(vId, seq);

    // check, TODO remove after tested
    assert(v->Id() == vId);

    // updates to state of graph
    vertex_map[vId] = v;
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
    sorted = false;
    return;
}

Vertex *HmmGraph::VertexGetter(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::VertexGetter]: Graph doesn't contain vertex %lld", i);
    }
    
    return vertex_map[i];

}

unsigned long HmmGraph::VertexOutDegree(int64_t i) {
    return VertexGetter(i)->OutDegree();
}

unsigned long HmmGraph::VertexInDegree(int64_t i) {
    return VertexGetter(i)->OutDegree();
}

const std::set<int64_t>& HmmGraph::InNeighbors(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::InNeighbors]: Graph does not contain vertex %lld", i);
    }

    return VertexGetter(i)->InNeighbors();
}

const std::set<int64_t>& HmmGraph::OutNeighbors(int64_t i) {
    if (!ContainsVertex(i)) {
        throw ParcoursException("[HmmGraph::OutNeighbors]: Graph does not contain vertex %lld", i);
    }
    
    return VertexGetter(i)->OutNeighbors();
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

int64_t HmmGraph::K() { return nVertices; }
