#include "partitioner.h"


void Partitioner::run()
{
    int n_part = 0;

    for (int i = 0; i < n_part; ++i) {



    }
    CalculateCost();
}


void Partitioner::CalculateCost()
{
    cost_ = 0;
    for (const auto &edge : hgraph_->GetHedges())
    {
        if (edge.is_cutted)
        {
            cost_ += edge.weight;
        }
    }
}