#include "refine.h"
#include "boundaryfm.h"
#include "balance.h"
#include "checkhgraph.h"
#include "fm.h"
using namespace std;

void s_refine::refine_2way(s_hgraph * original_hgraph, s_hgraph * hgraph, vector<float> tpwgts2)
{


	int oldcut = compute_mincut(hgraph);

	/*�������ͼ��cutsign\ bndlist\ bndptr\te\fs\nbnd*/
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

		project_2way_partition(hgraph);		/*��hgraph->coarser�����ֽ��Ͷ�䵽hgraph��,������bnd\fs\te���������*/

		//cout << endl;

	}

#if DEBUG
	cout << "The final cutsize = " << hgraph->mincut << endl;
#endif

}


/*�������ͼ��cutsign\ bndlist\ bndptr\te\fs\nbnd*/
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


		/*�жϳ����Ƿ�cut*/
		if (temp.front() != temp.back()) {

			hgraph->cutsign[n] = 1;					/*������߱����ˣ���cutsign����Ϊ1*/


			vector<int> count(2, 0);
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];
				count[hgraph->where[k - 1]]++;
				hgraph->bndptr[k - 1] = 1;
			}

			if (count[0] == 1) {				/*����0���ڵ�������,��fs����*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 0)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}
			}
			if (count[1] == 1) {			/*����1���ڵ�������,��fs����*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 1)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}

			}

		}
		else {
			/*û���и�ĳ��ߣ����ж����te����*/
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
				k = hgraph->eind[i];
				hgraph->te[k - 1] += hgraph->hewgts[n];
			}

		}
	}

}



/********************************************************
/*��������� ��ͼ���ݽṹ
/*����ֵ��	��
/*���ܣ�		�����볬ͼ��coarser����ͼ�Ļ��ֽ��Ͷ�����,���������Ȩֵ
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

	cout << "Ͷ��ǰ�� " << coarser->mincut << " Ͷ��� " << hgraph->mincut << endl;

	compute_2way_refine_params(hgraph);

	delete coarser;
	coarser = NULL;
	hgraph->coarser = NULL;
}

/*�����hgraph->finer��ͼ��where/pwgts/cutsign/fs/te/bndlist/bndptr/nbnd*/
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

	/*where����*/
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

	/*����Ȩֵ*/
	for (i = 0; i < hgraph->nvtxs; ++i)
	{
		hgraph->pwgts[hgraph->where[i]] += hgraph->vwgts[i];
	}

	//cout << "Ͷ��ǰ�� " << coarser->mincut << " Ͷ��� " << hgraph->mincut << endl;

	//if (hgraph->mincut != coarser->mincut)
	//{
	//	int ctest = compute_mincut(coarser);

	//	cout << "Ͷ��ǰ��һ��" << ctest << endl;
	//}




	/*����fs��ts,ͬʱ���㳬��cutsign,bnd����*/
	for (n = 0; n < hgraph->nhedges; ++n)
	{
		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
		{
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}

		std::sort(temp.begin(), temp.end());

		/*�жϸó����Ƿ�cut*/
		if (temp.front() != temp.back()) 
		{	

			/*cutsign����Ϊ1*/
			hgraph->cutsign[n] = 1;					
			count[0] = 0;
			count[1] = 0;

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];
				count[hgraph->where[k-1]]++;
				hgraph->bndlist.insert(k);			/*���и�ĳ����������Ķ��㼴Ϊ�߽綥��*/
				hgraph->bndptr[k - 1] = 1;

			}

			if (count[0] == 1) 
			{	/*�ó����ڷ���0����Ψһ��һ������,��ö����Ӧ��fs����*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
				{
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 0)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}
			}

			if (count[1] == 1) 
			{	/*�ó����ڷ���1����Ψһ��һ������,��ö����Ӧ��fs����*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
				{
					k = hgraph->eind[i];
					if (hgraph->where[k - 1] == 1)
						hgraph->fs[k - 1] += hgraph->hewgts[n];
				}

			}

		}
		else 
		{	/*û���и�ĳ���*/
			hgraph->cutsign[n] = 0;
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) 
			{	/*���ж����te����*/
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