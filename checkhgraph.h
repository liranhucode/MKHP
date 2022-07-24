#pragma once

#ifndef CHECKHGRAPH_H
#define CHECKHGRAPH_H

#include "hgraph.h"
#include "myutil.h"


/*��鳬ͼ�Ƿ���Ϲ���*/
void check_hgraph(s_hgraph* hgraph);

/*��鳬���������Ƿ������ {1 - nvtxs}�����еĶ���*/
void check_vertex(s_hgraph *hgraph);


/*��鳬���������Ƿ�����˵�һ����ĳ���*/
void check_hyperedge(s_hgraph *hgraph);


void check_dependent_hyperedge(s_hgraph *hgraph);


int check_fs_and_te(s_hgraph *hgraph, int higain);


int compute_mincut(s_hgraph * hgraph);


void check_cutsign(s_hgraph * hgraph);

void check_fs_and_te(s_hgraph * hgraph);

#endif