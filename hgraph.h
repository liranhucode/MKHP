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

	//计算顶点度
	void compute_vertex_degree();


	//计算临近超边
	void compute_incidence_hyperedge();

	//计算出每个顶点的临近顶点
	void compute_adjcent_vertices();

	void print_hypergraph_info(string file, int nparts, int ubfactor);

	void compute_cutsign();

	float average_hyperedge_deg();

	float average_vertex_deg();

	float RatioOfHedgeAndVertex();

	/*判断超边是否是单个顶点的超边*/
	bool is_single_hyperedge(int & num);

	/*判断超边是否是独立的超边*/
	bool is_dependent_hyperedge(int &num);

public:

						
	int	nvtxs, nhedges;		/*nvtxs = 超图的顶点数, nhedges = 超图的超边数*/


	size_s tvwgts;			/*顶点的权值之和*/


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

	
	vector<int> label;			/*表示当前子图的顶点在原图上的序号*/
	
	//Initial partition
	vector<int> cutsign;


	//adjacency list of nvtxs
	vector<unordered_set <int>> adjcent_vertices;


	//An array that stores the incident edges of vertices 
	//vector<vector<int>> incidence_hyperedge;

	//存放的该顶点的临近超边
	vector<vector<int>> incidence_hyperedge;


	
	//fm优化阶段用到的
	int nbnd;
	unordered_set<int> bndlist;
	vector<int> bndptr;

	vector<float> fs;
	vector<float> te;      /*fs[i]————包含顶点i的超边的其他顶点全部与i不在一个分区，fs[i] + 1*/
						   /*te[i]————包含顶点i的超边，当有超边中所有顶点与i同在一个分区，te[i] + 1*/



	//partitions result
	vector<int> where;
	//number of the partition 0, 1, ..
	// This is an array of size nvtxs that returns the computed partition

	int mincut;
	// edgecut
	// This is an integer that returns the number of hyperedges that are being cut by the partitioning algorithm


	vector<int> pwgts;				/*分区顶点权值之和*/


									//The coarser hypergraph and finer hypergraph
	s_hgraph *coarser, *finer;


	///* K-way refinement parameters */
	//ckrinfo_t *ckrinfo;   /*!< The per-vertex cut-based refinement info */
};


#endif