//
// Created by eaton on 20/04/18.
//

#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
using namespace std;

template <typename vertex>
class weighted_graph {

private:

    class node {

    private:

        int index;
        vertex data;
        bool visited;

    public:

        node(int, vertex, bool);
        node();
        ~node();
        int get_index();
        vertex get_data();
        bool is_visited();
        void set_visited(bool);

    };

    class edge {

    private:

        node * current_vertex;
        node * destination_vertex;
        int weight;

    public:

        edge(node *, node *, int);
        edge(int);
        edge();
        ~edge();
        void set_weight(int);
        int get_weight_or_has_edge();

    };

    //You will need to add some data members here
    //to actually represent the graph internally,
    //and keep track of whatever you need to.

    //The graph_iterator class provides an iterator
    //over the vertices of the graph.
    //This is one of the harder parts, so if you're
    //not too comfortable with C++ leave this for last.
    //If you are, there are many ways of doing this,
    //as long as it passes the tests, it's okay.
    class graph_iterator {

    private:

        //You may need data members here.

    public:
        graph_iterator(const weighted_graph &);
        graph_iterator(const weighted_graph &, size_t);
        ~graph_iterator();
        graph_iterator operator=(const graph_iterator&);
        bool operator==(const graph_iterator&) const;
        bool operator!=(const graph_iterator&) const;
        graph_iterator operator++();
        graph_iterator operator++(int);
        const vertex operator*();
        const vertex* operator->();
    };

    //The neighbour_iterator class provides an iterator
    //over the neighbours of a given vertex. This is
    //probably harder (conceptually) than the graph_iterator.
    //Unless you know how iterators work.
    class neighbour_iterator {

    private:

        //You may need data members here.

    public:
        neighbour_iterator(const neighbour_iterator&);
        neighbour_iterator(const weighted_graph &, const vertex&);
        neighbour_iterator(const weighted_graph &, const vertex&, size_t);
        ~neighbour_iterator();
        neighbour_iterator operator=(const neighbour_iterator& it);
        bool operator==(const neighbour_iterator&) const;
        bool operator!=(const neighbour_iterator&) const;
        neighbour_iterator operator++();
        neighbour_iterator operator++(int);
        const std::pair<vertex, int> operator*();
        const std::pair<const vertex, int>* operator->();
    };

    size_t size;//the number of nodes in the graph
    weighted_graph<vertex>::node * start;
    vector<vector<edge>> * adj_matrix ;
    vector<node> * nodes;

public:


    weighted_graph(); //A constructor for weighted_graph. It should start empty.
    ~weighted_graph(); //A destructor. Depending on how you do things, this may
    //not be necessary.

    bool are_adjacent(const vertex&, const vertex&) const; //Returns true if the two vertices are
    //adjacent, false otherwise.
    bool has_vertex(const vertex&) const; //Returns true if the passed in vertex is
    //a vertex of the graph, false otherwise.

    void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
    void add_edge(const vertex&, const vertex&, const int&); //Adds an edge between the two vertices
    //with the given weight (as an int).

    void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
    void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.
    void set_edge_weight(const vertex&, const vertex&, const int&); //Changes the edge weight between the two
    //vertices to the new weight (the int).

    int get_edge_weight(const vertex&, const vertex&) const; //Returns the weight on the edge between the two vertices.
    int degree(const vertex&) const; //Returns the degree of the vertex.
    int weighted_degree(const vertex&); //Returns the sum of the weights on all the edges incident to the vertex.
    int num_vertices() const; //Returns the total number of vertices in the graph.
    int num_edges() const; //Returns the total number of edges in the graph (just the count, not the weight).
    int total_weight(); //Returns the sum of all the edge weights in the graph.

    std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
    std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

    //todo:before doing iterators and traversals, test the tree first
    graph_iterator begin(); //Returns a graph_iterator pointing to the start of the vertex set.
    graph_iterator end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

    neighbour_iterator neighbours_begin(const vertex&); //Returns a neighbour_iterator pointing to the start
    //of the neighbour set for the given vertex.
    neighbour_iterator neighbours_end(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end
    //of the neighbour set for the given vertex.

    std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they
    //are visited in by a depth-first traversal starting at
    //the given vertex.
    std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they
    //are visisted in by a breadth-first traversal starting
    //at the given vertex.

    weighted_graph<vertex> mst(); //Returns a minimum spanning tree of the graph.

    weighted_graph<vertex>::node * find_node_by_data(const vertex&) const;
    vertex find_vertex_by_index(const int &) const;
    weighted_graph<vertex>::edge * find_edge_by_vertices(const vertex & vertextOne, const vertex & vertextTwo) const ;
};

//Define all your methods down here (or move them up into the header, but be careful you don't double up).
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.


template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g){}
template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g, size_t start_pos){}
template <typename vertex> weighted_graph<vertex>::graph_iterator::~graph_iterator(){}
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator=(const graph_iterator& it){ auto g = graph_iterator(weighted_graph<vertex>()); return g; }
template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator==(const graph_iterator& it) const { return false; }
template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator!=(const graph_iterator& it) const { return false; }
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(){ auto g = graph_iterator(weighted_graph<vertex>()); return g; }
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(int){ auto g = graph_iterator(weighted_graph<vertex>()); return g; }
template <typename vertex> const vertex weighted_graph<vertex>::graph_iterator::operator*(){ auto v = vertex(); return v; }
template <typename vertex> const vertex* weighted_graph<vertex>::graph_iterator::operator->(){ return nullptr; }

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u) {}
template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u, size_t start_pos) {}
template <typename vertex> weighted_graph<vertex>::neighbour_iterator::~neighbour_iterator() {}
template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator=(const neighbour_iterator& it) { auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }
template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator==(const neighbour_iterator& it) const { return false; }
template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator!=(const neighbour_iterator& it) const { return false; }
template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++() { auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }
template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++(int){ auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }
template <typename vertex> const std::pair<vertex, int> weighted_graph<vertex>::neighbour_iterator::operator*(){ auto p = std::pair<vertex,int>(); return p; }
template <typename vertex> const std::pair<const vertex, int>* weighted_graph<vertex>::neighbour_iterator::operator->(){ return nullptr; }

template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::begin() {
    return graph_iterator(weighted_graph<vertex>());
}
template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::end() {
    return graph_iterator(weighted_graph<vertex>());
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_begin(const vertex& u) {
    return neighbour_iterator(*this, vertex());
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_end(const vertex& u) {
    return neighbour_iterator(weighted_graph<vertex>(), vertex());
}

template <typename vertex> weighted_graph<vertex>::node::node(int index, vertex data, bool visited) {
    this->index = index;
    this->data = data;
    this->visited = visited;
}

template <typename vertex> weighted_graph<vertex>::node::node() {}

template <typename vertex> weighted_graph<vertex>::node::~node() {

}

template <typename vertex> vertex weighted_graph<vertex>::node::get_data() {
    return this->data;
}

template <typename vertex> int weighted_graph<vertex>::node::get_index()  {
    return this->index;
}

template <typename vertex> bool weighted_graph<vertex>::node::is_visited() {
    return visited;
}

template <typename vertex> void weighted_graph<vertex>::node::set_visited(bool visited) {
    this->visited = visited;
}

template <typename vertex> weighted_graph<vertex>::edge::edge(node * current_vertex, node * destination_vertex, int weight) {
    this->current_vertex = current_vertex;
    this->destination_vertex = destination_vertex;
    this->weight = weight;
}

template <typename vertex> weighted_graph<vertex>::edge::edge(int weight) {
    this->weight = weight;
}

template <typename vertex> weighted_graph<vertex>::edge::edge() {}

template <typename vertex> weighted_graph<vertex>::edge::~edge() {

}

template <typename vertex> void weighted_graph<vertex>::edge::set_weight(int weight) {
    this->weight = weight;
}

template <typename vertex> int weighted_graph<vertex>::edge::get_weight_or_has_edge() {
    return weight;
}

template <typename vertex> weighted_graph<vertex>::weighted_graph() {
    start = nullptr;
    size = 0;
    adj_matrix = nullptr;
    nodes = nullptr;
}

template <typename vertex> weighted_graph<vertex>::~weighted_graph() {}

template <typename vertex> bool weighted_graph<vertex>::has_vertex(const vertex& u) const {
    return find(nodes->begin(),nodes->end(),u) == nullptr;//todo:zhebian findhaoxianghuizhijiefanhuinend
}

template <typename vertex> bool weighted_graph<vertex>::are_adjacent(const vertex& u, const vertex& v) const {
    //vertex are data
    edge * edge1 =
            adj_matrix[find_node_by_data(u)->get_index()][find_node_by_data(v)->get_index()];
    return (bool) edge1->get_weight_or_has_edge();
}

template <typename vertex> void weighted_graph<vertex>::add_vertex(const vertex& v) {
    if(!has_vertex) {
        node * node1 = new node((int)nodes->size()+1, v, false);
        nodes->push_back(*node1);
        size++;
        if(nodes->size() == 0) {
            start = node1;
        }
        delete node1;
        if(nodes->size() != 0) {
            //this will initialize the matrix with added number vertex all to 0.
            /*TODO:to be tested whether if the objects are the same so that changing one value won't change them all */
            edge edge1;
            edge1.set_weight(0);
            vector<vector<edge>> tem_vector
                    (nodes->size()+1,vector<edge>(nodes->size()+1,edge1));
            vector<vector<edge>> * tem_vector_ptr = & tem_vector;
            //then copy the source vector to the new vector
            for(int i = 0; i < nodes->size(); i++) {
                for(int j = 0; j < nodes->size(); j++) {
                    tem_vector_ptr[i][j] = adj_matrix[i][j];
                }
            }
            adj_matrix = tem_vector_ptr;
        }
    }
}

template <typename vertex> void weighted_graph<vertex>::add_edge(const vertex& u, const vertex& v, const int& weight) {
    if(!are_adjacent(u, v) && has_vertex(u) && has_vertex(v)){

        //using edge to encapsulate a matrix the following
        edge * edge1 = new edge(u, v, weight);
        adj_matrix[find_node_by_data(u)->get_index()][find_node_by_data(v)->get_index()] = edge1;
        delete(edge1);

        edge * edge2 = new edge(v, u, weight);
        adj_matrix[find_node_by_data(v)->get_index()][find_node_by_data(u)->get_index()] = edge2;
        delete(edge2);
    }
}

template <typename vertex> void weighted_graph<vertex>::remove_vertex(const vertex& u) {
    if(has_vertex(u)) {
        int index = find_node_by_data(u)->get_index();

        // delete associated adj_matrix as a row and remove edge
        for(int i = 0; i < nodes->size(); i++) {
            for(int j = 0; j < nodes->size(); j++) {
                if(index-1 == j) { // locate a row
                    edge edge1 = adj_matrix[i][j];
                    if(edge1.get_weight_or_has_edge()) {
                        adj_matrix[i].erase(adj_matrix[i].begin() + j);


                    }

                }
            }
        }
        //delete associate adj_matrix in a column
        adj_matrix->erase(adj_matrix->begin()+index);

        //delete nodes
        nodes->erase(nodes->begin()+index);
    }
}


template <typename vertex> void weighted_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
    edge * edge1 = find_edge_by_vertices(u,v);
    edge1->set_weight(0);
}

template <typename vertex> void weighted_graph<vertex>::set_edge_weight(const vertex& u, const vertex& v, const int& weight) {
    edge * edge1 = find_edge_by_vertices(u,v);
    edge1->set_weight(weight);
}

template <typename vertex> int weighted_graph<vertex>::get_edge_weight(const vertex& u, const vertex& v) const {
    edge * edge1 = find_edge_by_vertices(u,v);
    return edge1->get_weight_or_has_edge();
}

template <typename vertex> int weighted_graph<vertex>::degree(const vertex& u) const {
    int sum = 0;
    int index = find_node_by_data(u)->get_index();
    for(int i = 0; i < adj_matrix->size(); i++) {
        if(index-1 == i) {
            edge edge1 = adj_matrix[i][0];
            if(edge1.get_weight_or_has_edge()){
                sum++;
            }
        }
    }
    return sum;
}

template <typename vertex> int weighted_graph<vertex>::weighted_degree(const vertex& u) {
    int sum = 0;
    int index = find_node_by_data(u)->get_index();
    for(int i = 0; i < adj_matrix->size(); i++) {
        for (int j = 0; j < adj_matrix[i].size(); j++) {
            if(index-1 == i || index-1 == j) {
                edge edge1 = adj_matrix[i][j];
                if(edge1.get_weight_or_has_edge()){
                    sum += edge1.get_weight_or_has_edge();
                }
            }
        }
    }
    return sum;
}

template <typename vertex> int weighted_graph<vertex>::num_vertices() const {
    return (int)nodes->size();
}

template <typename vertex> int weighted_graph<vertex>::num_edges() const {
    int sum = 0;
    for(int i = 0; i < adj_matrix->size(); i++) {
        for (int j = 0; j < adj_matrix[i].size(); j++) {
            edge edge1 = adj_matrix[i][j];
            if(edge1.get_weight_or_has_edge()){
                sum++;
            }
        }
    }
    return sum;
}

template <typename vertex> int weighted_graph<vertex>::total_weight() {
    int sum = 0;
    for(int i = 0; i < adj_matrix->size(); i++) {
        for (int j = 0; j < adj_matrix[i].size(); j++) {
            edge edge1 = adj_matrix[i][j];
            edge1.get_weight_or_has_edge();
        }
    }
    return sum/2;
}

template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_vertices() {
    vector<vertex> verticies(0);
    for(int i = 0; i < nodes->size(); i++) {
        node node1 = nodes->at(nodes->begin()+i);
        verticies.push_back(node1.get_data());
    }
    return verticies;
}

template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_neighbours(const vertex& u) {
    vector<vertex> vertices(0);
    int index = find_node_by_data(u)->get_index();
    for(int j = 0; j < adj_matrix->size(); j++) {
        edge edge1 = adj_matrix[index-1][j];
        if(edge1.get_weight_or_has_edge()){
            vertices.push_back(find_vertex_by_index(j));
        }
    }
    return vertices;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::depth_first(const vertex& start_vertex){
    return std::vector<vertex>();
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::breadth_first(const vertex& start_vertex){
    return std::vector<vertex>();
}

template <typename vertex>	weighted_graph<vertex> weighted_graph<vertex>::mst() {
    return weighted_graph<vertex>();
}

template <typename vertex> typename weighted_graph<vertex>::node *  weighted_graph<vertex>::find_node_by_data(const vertex& vertex) const {
    for(vector<node>::iterator it = nodes->begin(); it != nodes->end(); it++) {
        if (*it->get_data() == vertex) {
            return &(*it);
        }
    }
    return nullptr;
}

template <typename vertex> vertex weighted_graph<vertex>::find_vertex_by_index(const int & index) const {
    for(vector<node>::iterator it = nodes->begin(); it != nodes->end(); it++) {
        if(it->get_index() == index) {
            return *it->get_data();
        }
    }
}

template <typename vertex> typename weighted_graph<vertex>::edge * weighted_graph<vertex>::find_edge_by_vertices(const vertex & vertexOne, const vertex & vertexTwo) const {
    int index1 = find_node_by_data(vertexOne)->get_index();
    int index2 = find_node_by_data(vertexTwo)->get_index();
    // delete associated adj_matrix as a row and remove edge
    for(int i = 0; i < nodes->size(); i++) {
        for(int j = 0; j < nodes->size(); j++) {
            if(index1-1 == i && index2-1 == j) { // locate a row
                  edge edge1 = adj_matrix[i][j];
                  edge * edge1_ptr = &edge1;
                  return edge1_ptr;
            }
        }
    }
}

#endif
