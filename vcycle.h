#pragma once
#ifndef  VCYCLE_H
#define  VCYCLE_H

#include "hgraph.h"

namespace v_cycle
{
	/*����Ѿ�������ķ������н�һ���Ż�*/
	void v_cycle(s_hgraph *hgraph, int npart, vector<float> tpwgts2);

	s_hgraph *restricted_coarsen_hgraph(s_hgraph *original_hgraph, int npart);


	/*�����Ƶ�ƥ�䷽��,����v-cycle�׶�*/
	void restricted_hyperedge_match(s_hgraph *hgraph, double maxvwgt);

	int compute_cut(s_hgraph *hgraph);

	void creat_coarser_hgraph(s_hgraph * hgraph, vector<int>coarser_where, int cnvtxs, vector<int> &hyperedge_contracted);

	void creat_rpq_for_hyperedge(s_hgraph * hgraph, vector<hyperedge>& queue);
}


#endif