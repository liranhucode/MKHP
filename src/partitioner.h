#pragma once
#include "hypergraph.h"
#include "coarse.h"
#include "refine.h"

class Partitioner
{

public:
    Partitioner(Hypergraph *hgraph) : hgraph_(hgraph) {}
    ~Partitioner() {}

    void run() {}
    void report() {}

private:
    Hypergraph *hgraph_;
    Option *option_;

    Coarse *coarse_;
    Refine *refine_;

};