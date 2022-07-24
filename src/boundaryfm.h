#pragma once
#ifndef BOUNDARYFM_H
#define BOUNDARYFM_H

#include <vector>
#include <queue>

#include "hgraph.h"
#include "mystruct.h"
#include "rpq.h"



/*只计算边界顶点*/
namespace boundaryfm {

	void boundary_fm_2way_refine(s_hgraph * hgraph, vector<float> tpwgts2);


	void update_adjcent_gain_and_bnd_info(s_hgraph *hgraph, vector<s_rpq> &queue, int vertex, vector<int> moved);


	void bnd_insert(priority_queue <gain_node> &queue, float gain, int vertex);

	int rpq_get_top(priority_queue <gain_node> &queue);



	void rpq_reset(vector<s_rpq >& queue);

	int compute_cutsize_and_cutsign(s_hgraph * hgraph);

	void recompute_fm_params(s_hgraph *hgraph);


	void modified_cut_sign(s_hgraph *hgraph);


	void make_where_to_move(vector<s_rpq> &queue, vector<float> tpwgts, vector<int> pwgts, int &from, int &to);
}


#endif