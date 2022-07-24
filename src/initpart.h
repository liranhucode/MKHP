#pragma once
#ifndef INITPART_H
#define INITPART_H

#include <vector>
#include <queue>

#include "hgraph.h"
#include "mystruct.h"
#include "rpq.h"



class s_initpart
{
public:

	/*Ëæ»ú¶þ·Ö*/
	void random_bisection(s_hgraph *hgraph, int ubfactor, vector<float> tpwgts2);

	/*region growing algorithm*/
	void grow_bisection(s_hgraph * hgraph, int ubfactor, vector<float> tpwgts2);

	void compute_params_for_bfm(s_hgraph *hgraph);

	void compute_params_for_fm(s_hgraph *hgraph);

	int find_maxdeg_vertex(s_hgraph *hgraph);


};

#endif