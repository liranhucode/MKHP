#pragma once

#ifndef BALANCE_H
#define BALANCE_H


#include <vector>
#include <queue>
#include "hgraph.h"
#include "mystruct.h"
#include "rpq.h"


void balance_2way(s_hgraph *hgraph, vector<float> tpwgts2);

/*更新临近顶点的增益信息,保存至唯一的一个队列中*/
void balance_update_adjcent_info(s_hgraph *hgraph, s_rpq &queue, int vertex, vector<int> moved);


#endif