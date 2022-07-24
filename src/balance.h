#pragma once

#ifndef BALANCE_H
#define BALANCE_H


#include <vector>
#include <queue>
#include "hgraph.h"
#include "mystruct.h"
#include "rpq.h"


void balance_2way(s_hgraph *hgraph, vector<float> tpwgts2);

/*�����ٽ������������Ϣ,������Ψһ��һ��������*/
void balance_update_adjcent_info(s_hgraph *hgraph, s_rpq &queue, int vertex, vector<int> moved);


#endif