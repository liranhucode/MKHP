#pragma once

#ifndef FM_H
#define FM_H
#include <vector>
#include <queue>

#include "hgraph.h"
#include "mystruct.h"
#include "rpq.h"


#include "hgraph.h"


/*计算所有的顶点*/
namespace fm {

	void fm_2way_refine(s_hgraph * hgraph, vector<float> tpwgts2);

	void update_adjcent_gain_and_bnd_info(s_hgraph *hgraph, vector<s_rpq> &queue, int vertex, vector<int> vstatus);

	void bnd_insert(priority_queue <gain_node> &queue, float gain, int vertex);

	int rpq_get_top(priority_queue <gain_node> &queue);

	void rpq_reset(vector<s_rpq >& queue);

	int compute_cutsize_and_cutsign(s_hgraph * hgraph);

	void recompute_fm_params(s_hgraph *hgraph);

	void modified_cut_sign(s_hgraph *hgraph);
}

#endif