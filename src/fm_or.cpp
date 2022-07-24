#include "fm_or.h"


void fm_or::FMAlgorithm(s_hgraph *hgraph, vector<float> tpwgts2)
{

	int n, i, j, k, gain, from, to, higain, mincutorder, nswaps;
	int avgvwgt, origdiff, mindiff, limit;
	int maxpartwgt[2];
	int minpartwgt[2];
	vector<float> tpwgts(2, 0);

	tpwgts[0] = hgraph->tvwgts * tpwgts2[0];
	tpwgts[1] = hgraph->tvwgts - tpwgts[0];

	maxpartwgt[0] = tpwgts[0] * 1.05;
	maxpartwgt[1] = tpwgts[1] * 1.05;
	minpartwgt[0] = tpwgts[0] * 0.95;
	minpartwgt[1] = tpwgts[1] * 0.95;

	avgvwgt = hgraph->tvwgts / hgraph->nvtxs;

	origdiff = abs(tpwgts[0] - hgraph->pwgts[0]);

	int newcut, initcut, mincut;
	int testcut;

	limit = (int)min(max(0.01*hgraph->nvtxs, 50.0), 100.0);

	map<int, Node*>     bList[2];      // bucket list of partition A(0) and B(1)

	vector<int> vstatus(hgraph->nvtxs, RPQ_NOTPRESENT);
	vector<int> moved(hgraph->nvtxs, -1);
	vector<int> swaps(hgraph->nvtxs, 0);

	//computing initial cutsize
	for (int pass = 0; pass < 20; pass++)
	{

		mincutorder = -1;
		mincut = newcut = initcut = hgraph->mincut;
		mindiff = abs(tpwgts[0] - hgraph->pwgts[0]);

		buildBList(bList, hgraph, 2, 100);

		for (k = 1; k <= hgraph->nvtxs; ++k) {
			vstatus[k - 1] = RPQ_PRESENT;			/*将顶点插入到队列中，并标记为present*/
		}

		for (nswaps = 0; nswaps < hgraph->nvtxs; nswaps++)
		{

			from = (tpwgts[0] - hgraph->pwgts[0] < tpwgts[1] - hgraph->pwgts[1] ? 0 : 1);

			to = (from + 1) % 2;


			if ((higain = queue[from].rpq_get_top()) == -1) {
				break;
			}

			newcut -= (hgraph->fs[higain - 1] - hgraph->te[higain - 1]);

			hgraph->pwgts[from] -= hgraph->vwgts[higain - 1];

			hgraph->pwgts[to] += hgraph->vwgts[higain - 1];

			/*如果newcut小于mincut且p0分区权值与目标权值相差不超过 原始的差距 + 平均顶点权值,则更新mincut*/
			if (((newcut < mincut) && (abs(tpwgts[0] - hgraph->pwgts[0]) <= origdiff + avgvwgt)) ||
				((newcut == mincut) && (abs(tpwgts[0] - hgraph->pwgts[0]) < mindiff)) ||
				((newcut < mincut) && (hgraph->pwgts[0] < maxpartwgt[0]) && (hgraph->pwgts[0] > minpartwgt[0])) ||
				((newcut == mincut) && (hgraph->pwgts[0] < maxpartwgt[0]) && (hgraph->pwgts[0] > minpartwgt[0]))
				) {
				mincut = newcut;

				mindiff = abs(tpwgts[0] - hgraph->pwgts[0]);

				mincutorder = nswaps;
			}
			else if (nswaps - mincutorder > limit) { /*连续15次都没有满足约束，则退出本次pass */

				newcut += (hgraph->fs[higain - 1] - hgraph->te[higain - 1]);

				hgraph->pwgts[from] += hgraph->vwgts[higain - 1];

				hgraph->pwgts[to] -= hgraph->vwgts[higain - 1];

				break;
			}


			hgraph->where[higain - 1] = to;

			moved[higain - 1] = nswaps;

			swaps[nswaps] = higain;

			vstatus[higain - 1] = RPQ_EXTRACTED;			/*表示顶点higain从队列中取出*/

			/*testcut = compute_mincut(hgraph);

			cout << "v = " << higain << " where " << from << " -> " << to << \
				" newcut = " << newcut << " acutalcut= " << testcut << endl;*/



			/*更新与移动的顶点higain相关顶点的fs和te,以及bnd情况*/
			update_adjcent_gain_and_bnd_info(hgraph, queue, higain, vstatus);

			check_fs_and_te(hgraph);

			check_cutsign(hgraph);
		}

		/****************************************************************
		* Roll back computations //moved存放的移动的顺序 ，nswap代表交换的顺序
		*****************************************************************/
		for (i = 0; i < nswaps; i++)
		{
			moved[swaps[i] - 1] = -1;  /* reset moved array */
		}


		for (nswaps--; nswaps > mincutorder; nswaps--) {
			higain = swaps[nswaps];	//取出最后一步移动的
			to = hgraph->where[higain - 1] = (hgraph->where[higain - 1] + 1) % 2;	//如果higain在0分区，则返回到1

			hgraph->pwgts[to] += hgraph->vwgts[higain - 1];
			hgraph->pwgts[(to + 1) % 2] -= hgraph->vwgts[higain - 1];

		}

		hgraph->mincut = mincut;


		recompute_fm_params(hgraph);


		if (mincutorder <= 0 || mincut == initcut) {
			/*testcut = compute_mincut(hgraph);
			cout << "testcut = " << testcut << " mincut = " << mincut << endl;*/
			break;
		}




	}



}


void fm_or::buildBList(s_hgraph *hgraph, int length, int vertex_degree)
{
	bList[0].clear();
	bList[1].clear();

	
	for (int i = -vertex_degree; i <= vertex_degree; ++i) {
		if (bList[0].count(i) == 0)
			bList[0][i] = new Node(-1);    // dummy node
		if (bList[1].count(i) == 0)
			bList[1][i] = new Node(-1);    // dummy node
	}

	for (size_t i = 0; i < hgraph->nvtxs; ++i) {
		insertCell(bList, hgraph, i);
	}

	return;
}

void fm_or::insertCell(map<int, Node*>  bList[], s_hgraph *hgraph, int i)
{
	int gain = hgraph->fs[i-1] - hgraph->te[i-1];

	bool part = hgraph->where[i-1];
	Node* node = new Node(i);

	node->setPrev(bList[part][gain]);
	node->setNext(bList[part][gain]->getNext());
	bList[part][gain]->setNext(node);
	if (node->getNext() != NULL)
		node->getNext()->setPrev(node);
	return;
}



Node* fm_or::findMaxGainCell(map<int, Node*>  bList[], bool part, int max_degree)
{
	int maxGain = max_degree;
	while (maxGain >= -max_degree && bList[part][maxGain]->getNext() == NULL) {
		--maxGain;
	}
	Node* maxGainCell = _cellArray[_bList[part][maxGain]->getNext()->getId()];

	return maxGainCell;
}


void fm_or::removeCell(Cell* c)
{
	Node* node = c->getNode();
	node->getPrev()->setNext(node->getNext());
	if (node->getNext() != NULL)
		node->getNext()->setPrev(node->getPrev());
	return;
}

// called when gain of a cell is updated
void fm_or::moveCell(Cell* c)
{
	removeCell(c);
	insertCell(c);
	return;
}