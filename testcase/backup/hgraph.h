#pragma once

#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <unordered_set>
#include <set>
#include "myutil.h"

using namespace std;

class s_hgraph
{
public:
	s_hgraph();
	~s_hgraph();

	int fmt;
	// fmt repsent type of hypergraph 
	// 00 or none : HGraphFile for unweighted hypergraph
	// 01 : HGraphFile for weighted hyperedges
	// 10 : HGraphFile for weighted vertices
	// 11 : HGraphFile for weighted hyperedges an vertices

	// Reading the hypergraph file: XXXXX.hgr
	void read_hypergraph(string hypergraph_file);

	void compute_total_weight();

	void compute_hyperedge_degree();

	//���㶥���
	void compute_vertex_degree();


	//�����ٽ�����
	void compute_incidence_hyperedge();

	//�����ÿ��������ٽ�����
	void compute_adjcent_vertices();

	void print_hypergraph_info(string file, int nparts, int ubfactor);

	void compute_cutsign();

	float average_hyperedge_deg();

	float average_vertex_deg();

	float RatioOfHedgeAndVertex();

	/*�жϳ����Ƿ��ǵ�������ĳ���*/
	bool is_single_hyperedge(int & num);

	/*�жϳ����Ƿ��Ƕ����ĳ���*/
	bool is_dependent_hyperedge(int &num);

public:

						
	int	nvtxs, nhedges;		/*nvtxs = ��ͼ�Ķ�����, nhedges = ��ͼ�ĳ�����*/


	size_s tvwgts;			/*�����Ȩֵ֮��*/


	vector<int> eptr;
	vector<int> eind;
	// eptr , eind
	// Two arrays that are used to describe the hyperedges in the graph


	vector<size_s> vwgts;
	//An array of size nvtxs that stores the weight of the vertices



	vector<size_s> hewgts;
	//An array of size nhedges that stores the weight of the hyperedges

	vector<int> hyperedge_degree;
	//An array of edge size, size is represent the number of incidence node

	vector<int> vertex_degree;


	vector <int> cmap;
	// An array that represent the maps from orginal hypergraph to coarser hypergraph 

	
	vector<int> label;			/*��ʾ��ǰ��ͼ�Ķ�����ԭͼ�ϵ����*/
	
	//Initial partition
	vector<int> cutsign;


	//adjacency list of nvtxs
	vector<unordered_set <int>> adjcent_vertices;


	//An array that stores the incident edges of vertices 
	//vector<vector<int>> incidence_hyperedge;

	//��ŵĸö�����ٽ�����
	vector<vector<int> > incidence_hyperedge;


	
	//fm�Ż��׶��õ���
	int nbnd;
	std::unordered_set<int> bndlist;
	vector<int> bndptr;

	vector<float> fs;
	vector<float> te;      /*fs[i]����������������i�ĳ��ߵ���������ȫ����i����һ��������fs[i] + 1*/
						   /*te[i]����������������i�ĳ��ߣ����г��������ж�����iͬ��һ��������te[i] + 1*/



	//partitions result
	vector<int> where;
	//number of the partition 0, 1, ..
	// This is an array of size nvtxs that returns the computed partition

	int mincut;
	// edgecut
	// This is an integer that returns the number of hyperedges that are being cut by the partitioning algorithm


	vector<int> pwgts;				/*��������Ȩֵ֮��*/


									//The coarser hypergraph and finer hypergraph
	s_hgraph *coarser, *finer;


	///* K-way refinement parameters */
	//ckrinfo_t *ckrinfo;   /*!< The per-vertex cut-based refinement info */
};


#endif