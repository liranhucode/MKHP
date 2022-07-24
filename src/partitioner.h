#pragma once
#include "hypergraph.h"
#include "coarse.h"
#include "refine.h"

class Partitioner
{

public:
    Partitioner() {}
    Partitioner(std::shared_ptr<Hypergraph> hgraph) : hgraph_(hgraph) {}
    ~Partitioner() {}

    void run();
    void report() {

        std::cout << "Partitioner cost: " << cost_ << std::endl;
    }

    double GetCost() const { return cost_; }
    void CalculateCost();

private:
    std::shared_ptr<Hypergraph> hgraph_;
    Option option_;
    Coarse coarse_;
    Refine refine_;
    double cost_;

};