#pragma once

#pragma once
#ifndef PARTITION_H
#define PARTITION_H

//#include "coarsening.h"
#include "hgraph.h"
#include "coarsen.h"
#include "initpart.h"
#include "refine.h"
#include "multiobjective.h"

class s_rpart {

public:

	s_rpart();
	~s_rpart();


	/*多级递归划分顶层函数*/
	void mlevel_recursive_partition(s_hgraph * hgraph, int npart, int ubfactor, vector<int> &part, int &objval);

	/*多级二分法*/
	void  mlevel_bisection(s_hgraph * hgraph, vector<float> tpwgts2, int ubfactor, int npart);

	/*多级递归二分函数*/
	int mlevel_recursive_bisection(s_hgraph * hgraph, int npart, vector<float> tpwgts, vector<int> &part, int ubfactor, int fpart);


	/*多级递归划分顶层函数*/
	void mlevel_recursive_partition_api(int nvtxs, int nhedges, vector<int>eptr, vector<int> eind, \
		vector<float> vwgts, vector<float> hewgts, int npart, int ubfactor, \
		vector<int> &part, int &edgecut); 

	/*计算划分的外部度之和*/
	void compute_soed(s_hgraph * hgraph);



	/*设置超图*/
	s_hgraph *setup_hgraph(int nvtxs, int nhedges, vector<int>eptr, vector<int> eind, vector<float> vwgts, vector<float> hewgts);


	/*分配分区权值比例,默认均分*/
	void allocate_tpwgts(vector<float> &tpwgts, int npart);


	void split_hypergraph(s_hgraph *hgraph, s_hgraph * &lhgraph, s_hgraph *&rhgraph);

	void split_hypergraph_c(s_hgraph * hgraph, s_hgraph * &lhgraph, s_hgraph * &rhgraph);

	void reallocate_tpwgts(vector<float> &tpwgts, int npart, vector<float> &left_tpwgts, vector<float> &right_tpwgts);

public:
	s_hgraph * hgraph;
	vector<float> tpwgts;

	s_coarsen *coarsen;

	s_initpart *init_2way_part;

	s_refine *refinement;


	//objective function: min-cut
	int cutsize;

	//result array
	vector<int> part;


	int soed;

};


#endif