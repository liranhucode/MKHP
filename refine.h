#pragma once
#ifndef REFINE_H
#define REFINE_H

#include "hgraph.h"
#include "initpart.h"



class s_refine {

public:
	void compute_2way_refine_params(s_hgraph * hgraph);

	void refine_2way(s_hgraph * original_hgraph, s_hgraph * hgraph, vector<float> tpwgts2);

	void project_2way_partition(s_hgraph * hgraph);

	void project_2way_partition_c(s_hgraph * hgraph);
};


#endif