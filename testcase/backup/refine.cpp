#include "refine.h"
#include "boundaryfm.h"
#include "balance.h"
#include "checkhgraph.h"
#include "fm.h"
using namespace std;

void s_refine::refine_2way(s_hgraph * original_hgraph, s_hgraph * hgraph, vector<float> tpwgts2)
{


	int oldcut = compute_mincut(hgraph);

	/*计算出超图的cutsign\ bndlist\ bndptr\te\fs\nbnd*/
	compute_2way_refine_params(hgraph);

	int level = 0;

	for (;;) {

		level++;

		balance_2way(hgraph, tpwgts2);

		oldcut = hgraph->mincut;	/*test*/

		//cout << endl;

		boundaryfm::boundary_fm_2way_refine(hgraph, tpwgts2);

		//fm::fm_2way_refine(hgraph, tpwgts2);


#if DEBUG
		cout << " [ " << hgraph->pwgts[0] << "  " << hgraph->pwgts[1] << " " << oldcut << " " << hgraph->mincut << " ] " << endl;
#endif

		if (hgraph == original_hgraph)
		{
			break;
		}
			


		hgraph = hgraph->finer;

		project_2way_partition(hgraph);		/*将hgraph->coarser级划分结果投射到hgraph中,并计算bnd\fs\te等相关内容*/

		//cout << endl;

	}

#if DEBUG
	cout << "The final cutsize = " << hgraph->mincut << endl;
#endif

}


/*计算出超图的cutsign\ bndlist\ bndptr\te\fs\nbnd*/
void s_refine::compute_2way_refine_params(s_hgraph * hgraph)
{

	int n = 0;
	int i = 0;
	int k = 0;

	hgraph->fs.clear();
	hgraph->te.clear();
	hgraph->cutsign.clear();
	hgraph->bndptr.clear();
	hgraph->nbnd = 0;


	vector<float> fs(hgraph->nvtxs, 0.0);
	vector<float> te(hgraph->nvtxs, 0.0);
	vector<int> cutsign(hgraph->nhedges, 0);
	vector<int> bndptr(hgraph->nvtxs, 0);

	hgraph->fs = fs;
	hgraph->te = te;
	hgraph->cutsign = cutsign;
	hgraph->bndptr = bndptr;

	for (n = 0; n < hgraph->nhedges; ++n)
	{

		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}

		std::sort(temp.begin(), temp.end());


		/*判断超边是否被cut*/
		if (temp.front() != temp.back()) {

			hgraph->cutsign[n] = 1;					/*如果超边被切了，则cutsign设置为1*/


			vector<int> count(2, 0);
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];
				count[hgraph->where[k - 1]]++;
				hgraph->bndptr[k - 1] = 1;
			}

			if (count[0] == 1) {				/*分区0存在单独超边,则fs增加*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 0)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}
			}
			if (count[1] == 1) {			/*分区1存在单独超边,则fs增加*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 1)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}

			}

		}
		else {
			/*没被切割的超边，所有顶点的te增加*/
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
				k = hgraph->eind[i];
				hgraph->te[k - 1] += hgraph->hewgts[n];
			}

		}
	}

}



/********************************************************
/*输入参数： 超图数据结构
/*返回值：	空
/*功能：		将输入超图的coarser级超图的划分结果投射回来,并计算分区权值
*********************************************************/
void s_refine::project_2way_partition_c(s_hgraph * hgraph)
{
	int i = 0;
	int k = 0;
	s_hgraph *coarser = hgraph->coarser;
	int oldcut = coarser->mincut;
	vector<int> pwgts(2, 0);
	vector<int> where(hgraph->nvtxs, 0);

	for (i = 0; i < hgraph->nvtxs; ++i)
	{
		k = hgraph->cmap[i];
		if (k == -1) 
		{
			where[i] = 0;
			continue;
		}
		where[i] = coarser->where[k - 1];
	}

	for (i = 0; i < hgraph->nvtxs; ++i)
	{
		pwgts[where[i]] += hgraph->vwgts[i];
	}

	hgraph->pwgts = pwgts;
	hgraph->where = where;
	hgraph->mincut = compute_mincut(hgraph);

	cout << "投射前： " << coarser->mincut << " 投射后： " << hgraph->mincut << endl;

	compute_2way_refine_params(hgraph);

	delete coarser;
	coarser = NULL;
	hgraph->coarser = NULL;
}

/*计算出hgraph->finer超图的where/pwgts/cutsign/fs/te/bndlist/bndptr/nbnd*/
void s_refine::project_2way_partition(s_hgraph * hgraph)
{
	int n, i, j, k;
	int newcut;
	int oldcut;
	newcut = 0;

	s_hgraph *coarser = hgraph->coarser;
	oldcut = coarser->mincut;

	vector<int> count(2, 0);
	vector<int> pwgts(2, 0);
	vector<int> bndptr(hgraph->nvtxs, 0);
	vector<float> fs(hgraph->nvtxs, 0.0);
	vector<float> te(hgraph->nvtxs, 0.0);
	vector<int> cutsign(hgraph->nhedges, 0);
	vector<int> where(hgraph->nvtxs, 0);


	hgraph->pwgts = pwgts;
	hgraph->bndptr = bndptr;
	hgraph->fs = fs;
	hgraph->te = te;
	hgraph->cutsign = cutsign;
	hgraph->where = where;

	/*where数组*/
	for (i = 0; i < hgraph->nvtxs; ++i)
	{
		k = hgraph->cmap[i];
		if (k == -1){
			hgraph->where[i] = 0;
			continue;
		}
		hgraph->where[i] = coarser->where[k - 1];
	}

	hgraph->mincut = compute_mincut(hgraph);

	/*分区权值*/
	for (i = 0; i < hgraph->nvtxs; ++i)
	{
		hgraph->pwgts[hgraph->where[i]] += hgraph->vwgts[i];
	}

	//cout << "投射前： " << coarser->mincut << " 投射后： " << hgraph->mincut << endl;

	//if (hgraph->mincut != coarser->mincut)
	//{
	//	int ctest = compute_mincut(coarser);

	//	cout << "投射前后不一样" << ctest << endl;
	//}




	/*计算fs和ts,同时计算超边cutsign,bnd数组*/
	for (n = 0; n < hgraph->nhedges; ++n)
	{
		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
		{
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}

		std::sort(temp.begin(), temp.end());

		/*判断该超边是否被cut*/
		if (temp.front() != temp.back()) 
		{	

			/*cutsign设置为1*/
			hgraph->cutsign[n] = 1;					
			count[0] = 0;
			count[1] = 0;

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];
				count[hgraph->where[k-1]]++;
				hgraph->bndlist.insert(k);			/*被切割的超边所包含的顶点即为边界顶点*/
				hgraph->bndptr[k - 1] = 1;

			}

			if (count[0] == 1) 
			{	/*该超边在分区0存在唯一的一个顶点,则该顶点对应的fs增加*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
				{
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 0)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}
			}

			if (count[1] == 1) 
			{	/*该超边在分区1存在唯一的一个顶点,则该顶点对应的fs增加*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
				{
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 1)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}

			}

		}
		else 
		{	/*没被切割的超边*/
			hgraph->cutsign[n] = 0;
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
			{	/*所有顶点的te增加*/
				k = hgraph->eind[i];
				hgraph->te[k - 1] += hgraph->hewgts[n];
			}

		}
	}

	hgraph->nbnd = hgraph->bndlist.size();







	delete coarser;
	coarser = NULL;
	hgraph->coarser = NULL;
}