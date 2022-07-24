#pragma once

#ifndef CHECKHGRAPH_H
#define CHECKHGRAPH_H

#include "hgraph.h"
#include "myutil.h"


/*检查超图是否符合规则*/
void check_hgraph(s_hgraph* hgraph);

/*检查超边数组中是否包含了 {1 - nvtxs}中所有的顶点*/
void check_vertex(s_hgraph *hgraph);


/*检查超边数组中是否包含了单一顶点的超边*/
void check_hyperedge(s_hgraph *hgraph);


void check_dependent_hyperedge(s_hgraph *hgraph);


int check_fs_and_te(s_hgraph *hgraph, int higain);


int compute_mincut(s_hgraph * hgraph);


void check_cutsign(s_hgraph * hgraph);

void check_fs_and_te(s_hgraph * hgraph);

#endif