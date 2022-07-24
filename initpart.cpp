
#include <algorithm>
#include <queue>
#include "initpart.h"
#include "rpq.h"
#include "boundaryfm.h"
#include "fm.h"
#include "checkhgraph.h"
#include "balance.h"

/*首先采用随机二分,然后通过booundary算法优化*/
void s_initpart::random_bisection(s_hgraph * hgraph, int ubfactor, vector<float> tpwgts2)
{
	vector<int> biwhere;

	int n = 0;
	int nvtxs = hgraph->nvtxs;
	int i = 0;
	int ii = 0;
	int inbfs = 0;
	float tpwgts;
	int init_cut;
	int bestcut = 0, refine_cut;


	vector<int> bestwhere(nvtxs, -1);

	/*分区最大权值上限*/
	tpwgts = hgraph->tvwgts * tpwgts2[0];

	//check_hgraph(hgraph);

	for (inbfs = 0; inbfs < 10; inbfs++) //20
	{
		vector<int> perm;

		vector<int> init_pwgts(2, 0);

		vector<int> init_where(nvtxs, 1);		//所有分区初始化为1

		hgraph->pwgts = init_pwgts;

		hgraph->where = init_where;

		/*生成一个随机数序列*/
		rand_permute(perm, nvtxs, inbfs+10);

		hgraph->pwgts[1] = hgraph->tvwgts;

		hgraph->pwgts[0] = 0;

		/*随机分配*/
		for (ii = 0; ii < nvtxs; ii++) {

			i = perm[ii] - 1;
			if (hgraph->pwgts[0] > tpwgts)
				break;
			hgraph->where[i] = 0;
			hgraph->pwgts[0] += hgraph->vwgts[i];
			hgraph->pwgts[1] -= hgraph->vwgts[i];
		}

		init_cut = boundaryfm::compute_cutsize_and_cutsign(hgraph);
		//init_cut = fm::compute_cutsize_and_cutsign(hgraph);

		hgraph->mincut = init_cut;

		compute_params_for_fm(hgraph);

		//balance_2way(hgraph, tpwgts2);

#if DEBUG
		int testcut = compute_mincut(hgraph);
#endif		
		//计算优化后产生的边割
		boundaryfm::boundary_fm_2way_refine(hgraph, tpwgts2);
		//fm::fm_2way_refine(hgraph, tpwgts2);

#if DEBUG
		printf("[%5d  %5d %5d %5d]\n", hgraph->pwgts[0], hgraph->pwgts[1], init_cut, hgraph->mincut);
#endif


		if (inbfs == 0 || bestcut > hgraph->mincut) {
			bestcut = hgraph->mincut;

			for (n = 0; n < nvtxs; ++n){
				bestwhere[n] = hgraph->where[n];
			}
			if (bestcut == 0)
				break;
		}


	}

	for (n = 0; n < nvtxs; ++n) {
		hgraph->where[n] = bestwhere[n];
	}



	hgraph->mincut = bestcut;

#if DEBUG
	//int testcut = compute_mincut(hgraph);
	//cout << "bestcut = " << bestcut /* << " " << testcut*/ << endl;
#endif
}




void s_initpart::grow_bisection(s_hgraph * hgraph, int ubfactor, vector<float> tpwgts2)
{

	vector<int> biwhere;

	int k = 0;
	int n = 0;
	int nvtxs = hgraph->nvtxs;
	int i = 0;
	int ii = 0;
	int inbfs = 0;
	float zeromaxpwgt;
	int init_cut;
	int bestcut = 0, refine_cut;


	vector<int> bestwhere(nvtxs, -1);
	//vector<int> bestpwgts(2, 0);

	zeromaxpwgt = hgraph->tvwgts * tpwgts2[0] * 1.05;


	for (inbfs = 0; inbfs < 20; inbfs++) //20
	{
		vector<int> perm;

		vector<int> init_pwgts(2, 0);
		vector<int> init_where(nvtxs, 1);		//所有分区初始化为1

		hgraph->pwgts = init_pwgts;
		hgraph->where = init_where;


		rand_permute(perm, nvtxs, inbfs + 10);
		hgraph->pwgts[1] = hgraph->tvwgts;
		hgraph->pwgts[0] = 0;

		/*区域增长分配*/
		int init_vertex;

		init_vertex = find_maxdeg_vertex(hgraph);

		for (i = 0; i < hgraph->incidence_hyperedge[init_vertex - 1].size(); ++i){

			n = hgraph->incidence_hyperedge[init_vertex-1][i];

			for (ii = hgraph->eptr[n]; ii < hgraph->eptr[n + 1]; ++ii){

				k = hgraph->eind[ii];

				hgraph->where[k-1] = 0;

				hgraph->pwgts[0] += hgraph->vwgts[k - 1];

				hgraph->pwgts[1] -= hgraph->vwgts[k - 1];
			}

			if (hgraph->pwgts[0] > zeromaxpwgt)
				break;
		}
		


		init_cut = fm::compute_cutsize_and_cutsign(hgraph);
		//init_cut = fm::compute_cutsize_and_cutsign(hgraph);

		hgraph->mincut = init_cut;

		compute_params_for_fm(hgraph);

		//计算优化后产生的边割
		fm::fm_2way_refine(hgraph, tpwgts2);
		//fm::fm_2way_refine(hgraph, tpwgts2);


#if DEBUG
		printf("[%5d  %5d %5d %5d]\n", hgraph->pwgts[0], hgraph->pwgts[1], init_cut, hgraph->mincut);
#endif

		if (inbfs == 0 || bestcut > hgraph->mincut) {
			bestcut = hgraph->mincut;

			for (n = 0; n < nvtxs; ++n) {
				bestwhere[n] = hgraph->where[n];
			}
			if (bestcut == 0)
				break;
		}


	}



	for (n = 0; n < nvtxs; ++n) {
		hgraph->where[n] = bestwhere[n];
	}


	hgraph->mincut = bestcut;


#if DEBUG
	int testcut = compute_mincut(hgraph);
	cout << "initial partition bestcut = " << bestcut  << " " << testcut << endl;
#endif

}






//根据where数组的信息,计算超图hgraph产生了多少边割, 并且标记被切割的超边
//bndlist
//bndptr
void s_initpart::compute_params_for_bfm(s_hgraph * hgraph)
{
	int n = 0;
	int i = 0;
	int k = 0;
	int count0, count1;

	vector<int> bndptr(hgraph->nvtxs, 0);
	vector<float> fs(hgraph->nvtxs, 0.0);
	vector<float> te(hgraph->nvtxs, 0.0);

	hgraph->bndlist.clear();


	hgraph->bndptr = bndptr;
	hgraph->fs = fs;
	hgraph->te = te;

	for (n = 0; n < hgraph->nhedges; ++n)
	{


		if (hgraph->cutsign[n] == 1){

			vector<int> count(2, 0);

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i){

				k = hgraph->eind[i];

				count[hgraph->where[k - 1]]++;

				hgraph->bndlist.insert(k);			//被cut的超边其中所有顶点就是边界顶点

				hgraph->bndptr[k - 1] = 1;

			}


			if (count[0] == 1){

				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					k = hgraph->eind[i];

					if( hgraph->where[k-1] == 0)
						hgraph->fs[k - 1] += hgraph->hewgts[n];

				}
			}
			
			if (count[1] == 1) {

				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					k = hgraph->eind[i];

					if (hgraph->where[k - 1] == 1)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}

			}



		}
		else {
			//没被切割的超边，所有顶点的te += 超边权值
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];

				hgraph->te[k - 1] += hgraph->hewgts[n];
			}
		}

	}

	hgraph->nbnd = hgraph->bndlist.size();

}

void s_initpart::compute_params_for_fm(s_hgraph * hgraph)
{

	int n = 0;
	int i = 0;
	int k = 0;
	int count0, count1;

	vector<float> fs(hgraph->nvtxs, 0.0);
	vector<float> te(hgraph->nvtxs, 0.0);

	hgraph->bndlist.clear();


	hgraph->fs = fs;
	hgraph->te = te;

	for (n = 0; n < hgraph->nhedges; ++n)
	{


		if (hgraph->cutsign[n] == 1) {

			vector<int> count(2, 0);

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];

				count[hgraph->where[k - 1]]++;
			}


			if (count[0] == 1) {

				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					k = hgraph->eind[i];

					if (hgraph->where[k - 1] == 0)
						hgraph->fs[k - 1] += hgraph->hewgts[n];

				}
			}

			if (count[1] == 1) {

				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					k = hgraph->eind[i];

					if (hgraph->where[k - 1] == 1)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}

			}



		}
		else {
			//没被切割的超边，所有顶点的te += 超边权值
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
				k = hgraph->eind[i];
				hgraph->te[k - 1] += hgraph->hewgts[n];
			}
		}

	}

	

}

int s_initpart::find_maxdeg_vertex(s_hgraph * hgraph)
{
	int maxdeg = 0;
	int vertex;

	for (int i = 0; i < hgraph->nvtxs; ++i){

		
		if (hgraph->vertex_degree[i] > maxdeg){

			maxdeg = hgraph->vertex_degree[i];

			vertex = i + 1;

		}

	}

	return vertex;
}


