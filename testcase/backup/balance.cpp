#include "balance.h"
#include <algorithm>

#include "boundaryfm.h"

void balance_2way(s_hgraph * hgraph, vector<float> tpwgts2)
{
	int i, k;
	int from, to;
	int mindiff, mincut, initcut;
	int nswaps, higain;
	int maxpartwgt[2];
	int minpartwgt[2];
	int target_pwgts[2];

	target_pwgts[0] = hgraph->tvwgts * tpwgts2[0];
	target_pwgts[1] = hgraph->tvwgts - target_pwgts[0];

	maxpartwgt[0] = target_pwgts[0] * 1.05;
	minpartwgt[0] = target_pwgts[0] * 0.95;

	maxpartwgt[1] = target_pwgts[1] * 1.05;
	minpartwgt[1] = target_pwgts[1] * 0.95;

	mindiff = abs(target_pwgts[0] - hgraph->pwgts[0]);

	/*从较大权值的分区向较小权值的分区移动*/
	from = (hgraph->pwgts[0] < target_pwgts[0] ? 1 : 0);
	to = (from + 1) % 2;

	vector<int> initpwgts = hgraph->pwgts;
	vector<int> moved(hgraph->nvtxs, -1);
	s_rpq queue;

	initcut = mincut = hgraph->mincut;

	/*将边界顶点计算增益后加入优先队列中*/
	for (auto it = hgraph->bndlist.begin(); it != hgraph->bndlist.end(); it++) {

		k = *it;
		if (hgraph->where[k-1] == from && hgraph->vwgts[k -1] <= mindiff)
			queue.rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
	}

	for (nswaps = 0; nswaps < hgraph->nvtxs; ++nswaps){

		if ((higain = queue.rpq_get_top()) == -1) {
			break;
		}


		/*判断to分区+待移动顶点权值 是否大于 to分区 目标权值*/
		if ((hgraph->fs[higain - 1] - hgraph->te[higain - 1] < 0) ||
			hgraph->pwgts[to] + hgraph->vwgts[higain - 1] > maxpartwgt[to])
			break;


		mincut  -= (hgraph->fs[higain - 1] - hgraph->te[higain - 1]);


		hgraph->where[higain-1] = to;
		hgraph->pwgts[from] -= hgraph->vwgts[higain - 1];
		hgraph->pwgts[to] += hgraph->vwgts[higain - 1];
		moved[higain-1] = nswaps;

		balance_update_adjcent_info(hgraph, queue, higain, moved);

	}

	hgraph->mincut = mincut;


#if DEBUG 	
	cout << "balance: [ " << initpwgts[0] << " " << initpwgts[1] << " " << initcut << " ]"\
		<< " --> [ " << hgraph->pwgts[0] << " " << hgraph->pwgts[1] << " "<< mincut << " ]"<< endl;
#endif
}


void balance_update_adjcent_info(s_hgraph * hgraph, s_rpq & queue, int vertex, vector<int> moved)
{

	int i, ii, j, k, n, nn;
	int deg;
	int part;
	vector<int> incidence_hyperedge = hgraph->incidence_hyperedge[vertex - 1];
	int count[2];


	//更新vertex顶点临近顶点的fs和te,通过遍历与vertex顶点相关的超边
	//总共分四种情况：
	//1、 移动前被cut -> 移动后被cut	  -> 有可能会增加某些顶点的fs,但不改变所有的顶点的te	
	//2、 移动前被cut -> 移动后不被cut  -> 有可能会减小某些顶点的fs,同时会增大所有顶点的te
	//3、 移动前未被cut -> 移动后被cut  -> 会减小所有顶点的te,同时会增大某些顶点的fs
	//4、 移动前未被cut -> 移动后也未被cut	 -> 所有顶点的fs和te均未改变
	for (i = 0; i < incidence_hyperedge.size(); ++i) {

		n = incidence_hyperedge[i];
		deg = hgraph->eptr[n + 1] - hgraph->eptr[n];

		count[0] = 0;
		count[1] = 0;

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {

			k = hgraph->eind[j];				/*临近顶点*/
			count[hgraph->where[k - 1]]++;
		}

		//cout << " n = " << n << ": "<< count[0] << " " << count[0] << endl;

		/*该超边已被切割过*/
		if (hgraph->cutsign[n] == 1) {

			/*现在该超边没被切割了*/
			if (count[0] * count[1] == 0) {

				hgraph->cutsign[n] = 0;		/*此时超边已被切割*/

				if (deg == 2) {	/*原先deg=2的超边被切割了，那么之前该超边对所有顶点的fs贡献为1*/

					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						hgraph->fs[k - 1] -= hgraph->hewgts[n];
					}
				}
				else {
					/*原先deg > 2的超边被切割了，那么之前该超边对移动的顶点vertex的fs贡献为1,
					所以移动后需要,减小fs,切所有顶点te都需要增加 */

					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						if (k == vertex) {
							hgraph->fs[k - 1] -= hgraph->hewgts[n];
						}

					}

				}
			}
			else  /*现在该超边依然被切割了*/
			{

				if (count[0] == 1) {		/*被切成{1,x}型,在0分区的顶点fs++*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];

						if (hgraph->where[k - 1] == 0) {
							hgraph->fs[k - 1] += hgraph->hewgts[n];
						}
						else {
							if (count[1] == 2)			/*在1分区的顶点数为2,则移动前在1分区的顶点数为1,会产生fs,所以fs--*/
							{
								if (k != vertex && hgraph->fs[k - 1] > 0)
									hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}

					}
				}
				else if (count[1] == 1) {			/*被切成{x,1}型, 在1分区的顶点fs++*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];

						if (hgraph->where[k - 1] == 1) {
							hgraph->fs[k - 1] += hgraph->hewgts[n];
						}
						else {
							if (count[0] == 2)		/*在0分区的顶点数为2, 则移动前在0分区的顶点数为1,会产生fs,所以fs--*/
							{
								if (k != vertex && hgraph->fs[k - 1] > 0)
									hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}
					}

				}/*被切成{x,x}型, 在1分区的顶点fs--*/
				else {

					part = (hgraph->where[vertex - 1] + 1) % 2;

					if (count[1] == 2 && part == 0)	/*具体来说应该是{x,2}型,且移动顶点从0->1*/
					{
						for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
							k = hgraph->eind[j];
							if (hgraph->where[k - 1] == 1 && k != vertex) {
								hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}

					}
					else if (count[0] == 2 && part == 1)	/*具体来说应该是{2,x}型,且移动顶点从1->0*/
					{
						for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
							k = hgraph->eind[j];
							if (hgraph->where[k - 1] == 0 && k != vertex) {
								hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}
					}


				}
			}

		}
		else if (hgraph->cutsign[n] == 0) {		/*原先的超边没有被切割*/

			if (count[0] == 1 && count[1] == 1) {		/*被切割成{1,1}型超边*/
				hgraph->cutsign[n] = 1;

				for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
					k = hgraph->eind[j];
					hgraph->te[k - 1] -= hgraph->hewgts[n];
					hgraph->fs[k - 1] += hgraph->hewgts[n];
				}
			}
			else if (count[0] == 1 && count[1] > 1) {	/*被切割成{1,x}型超边*/
				hgraph->cutsign[n] = 1;

				for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
					k = hgraph->eind[j];
					hgraph->te[k - 1] -= hgraph->hewgts[n];

					if (hgraph->where[k - 1] == 0) {
						hgraph->fs[k - 1] += hgraph->hewgts[n];
					}
				}
			}
			else if (count[0] > 1 && count[1] == 1) {	/*被切割成{x,1}型超边*/
				hgraph->cutsign[n] = 1;

				for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
					k = hgraph->eind[j];
					hgraph->te[k - 1] -= hgraph->hewgts[n];

					if (hgraph->where[k - 1] == 1) {
						hgraph->fs[k - 1] += hgraph->hewgts[n];
					}
				}
			}
		}

	}

	//更新bnd数组和更新queue队列
	//bndptr[i] = 1 && moved[i] != -1   ->   不在增益队列中
	//bndptr[i] = 0  ->  不在增益队列中
	//bndptr[i] = 1 && moved[i] = 0     ->   在增益队列中
	for (i = 0; i < incidence_hyperedge.size(); ++i) {
		n = incidence_hyperedge[i];

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {

			k = hgraph->eind[j];
			/*判断顶点是边界顶点*/
			float sum = 0.0;
			for (ii = 0; ii < hgraph->incidence_hyperedge[k - 1].size(); ++ii) {
				nn = hgraph->incidence_hyperedge[k - 1][ii];
				sum += hgraph->hewgts[nn];
			}


			if (hgraph->te[k - 1] == sum) {

				/*顶点关联的所有超边都贡献te,则为边界顶点*/
				if (hgraph->bndptr[k - 1] == 1) {
					/*是否存在于增益队列中*/
					hgraph->bndlist.erase(k);
					hgraph->bndptr[k - 1] = 0;
					if (moved[k - 1] == -1)
						queue.rpq_delete(k);
				}
				else {

					/*之前不存在与增益队列中*/
					hgraph->bndlist.insert(k);
					hgraph->bndptr[k - 1] = 1;
					queue.rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
				}

			}
			else {
				/*顶点关联的超边不全部贡献te,则为边界顶点*/
				if (hgraph->bndptr[k - 1] == 1) {
					if (moved[k - 1] == -1)
						queue.rpq_update(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
				}
				else {
					/*之前不存在与增益队列中*/
					hgraph->bndlist.insert(k);
					hgraph->bndptr[k - 1] = 1;
					queue.rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
				}


			}


		}
	}


	hgraph->nbnd = hgraph->bndlist.size();

}
