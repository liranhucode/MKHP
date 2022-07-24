#pragma once

#include "utility.h"

using index = int32_t;

class Hypernode
{
public:
    int id;
    int weight;
};

class Hyperedge
{
public:
    int id;
    int weight;
    std::vector<index> nodes;

};

class Hypergraph
{
public:
    Hypergraph(/* args */) {};
    ~Hypergraph() {};

private:
    int num_nodes_;
    std::vector<Hypernode> hnodes_;
    int num_edges_;
    std::vector<Hyperedge> hedges_;
};
