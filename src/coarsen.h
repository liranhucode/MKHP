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


	//=====================����ƥ���㷨=====================//

	/*����ֻ�����,����ϴ��ԭʼ��ͼ,���شֻ����С�ĳ�ͼ*/
	s_hgraph *coarsen_hgraph(s_hgraph *original_hgraph, int npart);

	/*�ֻ������в��õ�ƥ�䷽��*/

	/*�Գ���Ȩֵ��С����ĳ���ƥ�䷽��*/
	void hyperedge_match(s_hgraph *hgraph, double maxvwgt);

	/*��������, �Ӷ����������ĳ�����,Ѱ�������ֱ�ӽ���contracted�ĳ���*/
	void random_match(s_hgraph *hgraph, double maxvwgt);

	void creat_rpq_for_hyperedge(s_hgraph * hgraph, vector<hyperedge>& queue);

	void creat_coarser_hgraph(s_hgraph *&hgraph, const int &cnvtxs, vector<int> &hyperedge_contracted);


	//=====================����ƥ���㷨=====================//

	s_hgraph *CoarsenHypergraph(s_hgraph *original_hgraph, int npart);

	void VertexMatching( s_hgraph* &hypergraph );

	int findCandicate(s_hgraph* &hypergraph, int &u, vector<int> &unmatched);

	void constructNewCoarserHypergraph(s_hgraph* &hypergraph, int cnvtxs);

private:
	int coarsen_to;
	double RateOfContraction;
};



#endif
