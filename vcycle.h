#pragma once
#ifndef  VCYCLE_H
#define  VCYCLE_H

#include "hgraph.h"

namespace v_cycle
{
	/*针对已经计算出的分区进行进一步优化*/
	void v_cycle(s_hgraph *hgraph, int npart, vector<float> tpwgts2);

	s_hgraph *restricted_coarsen_hgraph(s_hgraph *original_hgraph, int npart);


	/*有限制的匹配方法,用于v-cycle阶段*/
	void restricted_hyperedge_match(s_hgraph *hgraph, double maxvwgt);

	int compute_cut(s_hgraph *hgraph);

	void creat_coarser_hgraph(s_hgraph * hgraph, vector<int>coarser_where, int cnvtxs, vector<int> &hyperedge_contracted);

	void creat_rpq_for_hyperedge(s_hgraph * hgraph, vector<hyperedge>& queue);
}


#endif