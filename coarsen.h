#pragma once
#ifndef COARSEN_H
#define COARSEN_H


#include <queue>

#include "hgraph.h"
#include "mystruct.h"

#define MATCHED 1
#define UNMATCHED 0


#define YES 1
#define NO 0


#define CONTRACTED 1
#define NO_CONTRACTED 0

#define DISCARD 888

class s_coarsen
{

public:


	//=====================超边匹配算法=====================//

	/*顶层粗化函数,输入较大的原始超图,返回粗化后较小的超图*/
	s_hgraph *coarsen_hgraph(s_hgraph *original_hgraph, int npart);

	/*粗化过程中采用的匹配方法*/

	/*以超边权值大小排序的超边匹配方法*/
	void hyperedge_match(s_hgraph *hgraph, double maxvwgt);

	/*遍历顶点, 从顶点所关联的超边中,寻找其可以直接进行contracted的超边*/
	void random_match(s_hgraph *hgraph, double maxvwgt);

	void creat_rpq_for_hyperedge(s_hgraph * hgraph, vector<hyperedge>& queue);

	void creat_coarser_hgraph(s_hgraph *&hgraph, const int &cnvtxs, vector<int> &hyperedge_contracted);


	//=====================顶点匹配算法=====================//

	s_hgraph *CoarsenHypergraph(s_hgraph *original_hgraph, int npart);

	void VertexMatching( s_hgraph* &hypergraph );

	int findCandicate(s_hgraph* &hypergraph, int &u, vector<int> &unmatched);

	void constructNewCoarserHypergraph(s_hgraph* &hypergraph, int cnvtxs);

private:
	int coarsen_to;
	double RateOfContraction;
};



#endif
