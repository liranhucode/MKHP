#pragma once
#include "utility.h"

using index = int32_t;

class Hypernode
{
public:
    Hypernode (int i, double w) : id(i), weight(w) {} 
    int id;
    int part_id = 0;
    double weight;
    std::vector<index> edges;
};

class Hyperedge
{
public:
    Hyperedge (int i, double w) : id(id), weight(w) {} 
    int id;
    double weight;
    std::vector<index> nodes;
    bool is_cutted = false;

};

class Hypergraph
{
public:
    Hypergraph(/* args */) {};
    ~Hypergraph() {};


    void Parser(std::string &file_name);
    void Report()
    {
        std::cout << " Num nodes: " << num_nodes_ << std::endl;
        std::cout << " Num edges: " << num_edges_ << std::endl;
    }

    std::vector<Hyperedge> GetHedges() { return hedges_; }
    std::vector<Hypernode> GetHnodes() { return hnodes_; }
    int GetNumEdges() const { return num_edges_; }
    int GetNumNodes() const { return num_nodes_; }

private:
    int num_nodes_;
    std::vector<Hypernode> hnodes_;
    int num_edges_;
    std::vector<Hyperedge> hedges_;
};
