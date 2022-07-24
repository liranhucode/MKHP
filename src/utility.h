#pragma once
#include <memory>
#include <vector>
#include "timer.h"


class Option
{
public:
    int n_parts;
    int ub_factor;
    int n_runs;
    int max_coarsen_to;
    float reduct_ratio;
    std::string outfile;
};
