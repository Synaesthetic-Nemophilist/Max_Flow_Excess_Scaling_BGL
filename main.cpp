#include <iostream>
#include <math.h>

#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/grid_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/graph_utility.hpp>



using namespace boost;
using namespace std;


//Type & Struct Definitions -----------------------------------------------------------------

// Bundled properties for graph
struct bundleVertex {
    int id;
    int excess;
    int distance;
};

struct bundleEdge {
    int residual;
};


// Define graph type
typedef adjacency_list<vecS, vecS, bidirectionalS, bundleVertex, bundleEdge> Graph;

//Define graph traits type
typedef graph_traits<Graph> Traits;

//Define descriptor for vertices
typedef Traits::vertex_descriptor Vertex;
//Define descriptor for edges
typedef Traits::edge_descriptor Edge;




// Function Definitions -----
void preprocess(Graph &G, Vertex source, Vertex sink);

// Utility function for selecting node with smallest distance among nodes with large excess
Vertex selectNode(Graph &G, Vertex source, Vertex sink, float delta);

// implements push of flow or relabel of distance of current node
void pushRelabel(Graph &G, Vertex current, float delta);



int main() {

    // Graph instantiation
    Graph G;

    // Graph creation
    Vertex v1 = add_vertex(bundleVertex{1, 0, 2}, G);
    Vertex v2 = add_vertex(bundleVertex{2, 0, 1}, G);
    Vertex v3 = add_vertex(bundleVertex{3, 0, 1}, G);
    Vertex v4 = add_vertex(bundleVertex{4, 0, 0}, G);

    add_edge(v1, v2, bundleEdge{2}, G);
    add_edge(v1, v3, bundleEdge{4}, G);
    add_edge(v2, v3, bundleEdge{3}, G);
    add_edge(v2, v4, bundleEdge{1}, G);
    add_edge(v3, v4, bundleEdge{5}, G);


    // Excess Scaling Algorithm-----

    //Preprocessing ----
    preprocess(G, v1, v4);  // TODO: REMOVE HARDCODE

    // Find the max capacity to calculate delta = 2 ^ log(maxCapacity)
    Graph::edge_iterator edgeIt, edgeEnd;
    tie(edgeIt, edgeEnd) = edges(G);
    std::vector<int> capacities(num_edges(G));
    for (; edgeIt!= edgeEnd; ++edgeIt) {
        capacities.push_back(G[*edgeIt].residual);
    }
    int maxCap = *max_element(capacities.begin(), capacities.end());
    float delta = (float) pow(2, ceil(log2(maxCap)));

    while(delta >= 1) {
        Vertex currentNode = selectNode(G, v1, v4, delta); // select node with smallest distance among nodes with large excess
        while(currentNode != 0) {
            // push flow or relabel distance of current node
            pushRelabel(G, currentNode, delta);
            // select another node from updated graph
            currentNode = selectNode(G, v1, v4, delta);
        }
        delta /= 2;  // update delta
    }

    print_graph(G, get(&bundleVertex::id, G));

    Graph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    cout << endl << endl;
    for (; vertexIt != vertexEnd; ++vertexIt) {
        cout << "excess of node "<< G[*vertexIt].id << " : " << G[*vertexIt].excess << endl;
    }

    //TODO: MAX FLOW IS EXCESS VALUE AT SINK...handle that


    return 0;
}



// Function Declarations -----
void preprocess(Graph &G, Vertex source, Vertex sink) {

    // Store initial distances in this vector
    std::vector<int> distances(num_vertices(G));

    // Computation of exact distance labels via backward BFS from sink
    reverse_graph<Graph> R = make_reverse_graph(G);
    breadth_first_search(R,
                         sink,
                         visitor(make_bfs_visitor(record_distances(make_iterator_property_map(distances.begin(),
                                                                                              get(vertex_index, G)),
                                                                   on_tree_edge()))));
    Graph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++vertexIt) {
        G[*vertexIt].distance = distances[*vertexIt];
    }

    // Saturate source's adjacent edges
    Graph::out_edge_iterator outedgeIt, outedgeEnd;
    tie(outedgeIt, outedgeEnd) = out_edges(source, G);
    for(; outedgeIt != outedgeEnd; ++outedgeIt) {
        Vertex src = boost::source(*outedgeIt, G);
        Vertex target = boost::target(*outedgeIt, G);
        int saturation = G[*outedgeIt].residual;

        add_edge(target, src, bundleEdge{saturation}, G);  // add reverse of saturated edge
        remove_edge(*outedgeIt, G); // remove saturated edge

        G[target].excess += saturation;  // update target node's excess property
    }

    // Update source node's distance label to num_vertices
    G[source].distance = (int) num_vertices(G);
}

Vertex selectNode(Graph &G, Vertex source, Vertex sink, float delta) {

    double min = std::numeric_limits<double>::infinity();
    Vertex minNode = 0;

    Graph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if(G[*vertexIt].excess >= delta/2 && *vertexIt != source && *vertexIt != sink) {
            if(G[*vertexIt].excess < min) {
                min = G[*vertexIt].excess;
                minNode = *vertexIt;
            }
        }
    }
    return minNode;
}

void pushRelabel(Graph &G, Vertex current, float delta) {

    // stores candidate distances of adj nodes for raising current node's distance label to
    std::vector<int> adjDist;

    Graph::out_edge_iterator outedgeIt, outedgeEnd;
    tie(outedgeIt, outedgeEnd) = out_edges(current, G);
    for(; outedgeIt != outedgeEnd; ++outedgeIt) {
        Vertex target = boost::target(*outedgeIt, G);

        if(G[current].distance == G[target].distance + 1) {  // if network contains an admissable arc
            int pushFlow = 0;
            if(delta - G[target].excess > 0) {
                pushFlow = std::min({G[current].excess, G[*outedgeIt].residual, (int) delta - G[target].excess});  // last arg ensures that excess doesnt exceed delta
            } else {
                pushFlow = std::min(G[current].excess, G[*outedgeIt].residual);
            }

            // implement Push
            G[current].excess -= pushFlow;
            G[target].excess += pushFlow;
            G[*outedgeIt].residual -= pushFlow;
            if(G[*outedgeIt].residual == 0)  // saturation
                remove_edge(*outedgeIt, G); // remove saturated edge
            bool found;
            Edge reverse;
            boost::tie(reverse, found) = edge(target, current, G);
            if(found) {
                G[reverse].residual += pushFlow;
            } else {
                add_edge(target, current, bundleEdge{pushFlow}, G);  // add reverse edge
            }

            return;
        } else {
            if(G[*outedgeIt].residual > 0)
                adjDist.push_back(G[target].distance + 1);
        }
    }

    // implement Relabel
    int minDist = *min_element(adjDist.begin(), adjDist.end());
    G[current].distance = minDist;

    return;
}












