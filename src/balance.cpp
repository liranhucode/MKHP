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

	/*�ӽϴ�Ȩֵ�ķ������СȨֵ�ķ����ƶ�*/
	from = (hgraph->pwgts[0] < target_pwgts[0] ? 1 : 0);
	to = (from + 1) % 2;

	vector<int> initpwgts = hgraph->pwgts;
	vector<int> moved(hgraph->nvtxs, -1);
	s_rpq queue;

	initcut = mincut = hgraph->mincut;

	/*���߽綥����������������ȶ�����*/
	for (auto it = hgraph->bndlist.begin(); it != hgraph->bndlist.end(); it++) {

		k = *it;
		if (hgraph->where[k-1] == from && hgraph->vwgts[k -1] <= mindiff)
			queue.rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
	}

	for (nswaps = 0; nswaps < hgraph->nvtxs; ++nswaps){

		if ((higain = queue.rpq_get_top()) == -1) {
			break;
		}


		/*�ж�to����+���ƶ�����Ȩֵ �Ƿ���� to���� Ŀ��Ȩֵ*/
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


	//����vertex�����ٽ������fs��te,ͨ��������vertex������صĳ���
	//�ܹ������������
	//1�� �ƶ�ǰ��cut -> �ƶ���cut	  -> �п��ܻ�����ĳЩ�����fs,�����ı����еĶ����te	
	//2�� �ƶ�ǰ��cut -> �ƶ��󲻱�cut  -> �п��ܻ��СĳЩ�����fs,ͬʱ���������ж����te
	//3�� �ƶ�ǰδ��cut -> �ƶ���cut  -> ���С���ж����te,ͬʱ������ĳЩ�����fs
	//4�� �ƶ�ǰδ��cut -> �ƶ���Ҳδ��cut	 -> ���ж����fs��te��δ�ı�
	for (i = 0; i < incidence_hyperedge.size(); ++i) {

		n = incidence_hyperedge[i];
		deg = hgraph->eptr[n + 1] - hgraph->eptr[n];

		count[0] = 0;
		count[1] = 0;

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {

			k = hgraph->eind[j];				/*�ٽ�����*/
			count[hgraph->where[k - 1]]++;
		}

		//cout << " n = " << n << ": "<< count[0] << " " << count[0] << endl;

		/*�ó����ѱ��и��*/
		if (hgraph->cutsign[n] == 1) {

			/*���ڸó���û���и���*/
			if (count[0] * count[1] == 0) {

				hgraph->cutsign[n] = 0;		/*��ʱ�����ѱ��и�*/

				if (deg == 2) {	/*ԭ��deg=2�ĳ��߱��и��ˣ���ô֮ǰ�ó��߶����ж����fs����Ϊ1*/

					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						hgraph->fs[k - 1] -= hgraph->hewgts[n];
					}
				}
				else {
					/*ԭ��deg > 2�ĳ��߱��и��ˣ���ô֮ǰ�ó��߶��ƶ��Ķ���vertex��fs����Ϊ1,
					�����ƶ�����Ҫ,��Сfs,�����ж���te����Ҫ���� */

					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						if (k == vertex) {
							hgraph->fs[k - 1] -= hgraph->hewgts[n];
						}

					}

				}
			}
			else  /*���ڸó�����Ȼ���и���*/
			{

				if (count[0] == 1) {		/*���г�{1,x}��,��0�����Ķ���fs++*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];

						if (hgraph->where[k - 1] == 0) {
							hgraph->fs[k - 1] += hgraph->hewgts[n];
						}
						else {
							if (count[1] == 2)			/*��1�����Ķ�����Ϊ2,���ƶ�ǰ��1�����Ķ�����Ϊ1,�����fs,����fs--*/
							{
								if (k != vertex && hgraph->fs[k - 1] > 0)
									hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}

					}
				}
				else if (count[1] == 1) {			/*���г�{x,1}��, ��1�����Ķ���fs++*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];

						if (hgraph->where[k - 1] == 1) {
							hgraph->fs[k - 1] += hgraph->hewgts[n];
						}
						else {
							if (count[0] == 2)		/*��0�����Ķ�����Ϊ2, ���ƶ�ǰ��0�����Ķ�����Ϊ1,�����fs,����fs--*/
							{
								if (k != vertex && hgraph->fs[k - 1] > 0)
									hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}
					}

				}/*���г�{x,x}��, ��1�����Ķ���fs--*/
				else {

					part = (hgraph->where[vertex - 1] + 1) % 2;

					if (count[1] == 2 && part == 0)	/*������˵Ӧ����{x,2}��,���ƶ������0->1*/
					{
						for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
							k = hgraph->eind[j];
							if (hgraph->where[k - 1] == 1 && k != vertex) {
								hgraph->fs[k - 1] -= hgraph->hewgts[n];
							}
						}

					}
					else if (count[0] == 2 && part == 1)	/*������˵Ӧ����{2,x}��,���ƶ������1->0*/
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
		else if (hgraph->cutsign[n] == 0) {		/*ԭ�ȵĳ���û�б��и�*/

			if (count[0] == 1 && count[1] == 1) {		/*���и��{1,1}�ͳ���*/
				hgraph->cutsign[n] = 1;

				for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
					k = hgraph->eind[j];
					hgraph->te[k - 1] -= hgraph->hewgts[n];
					hgraph->fs[k - 1] += hgraph->hewgts[n];
				}
			}
			else if (count[0] == 1 && count[1] > 1) {	/*���и��{1,x}�ͳ���*/
				hgraph->cutsign[n] = 1;

				for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
					k = hgraph->eind[j];
					hgraph->te[k - 1] -= hgraph->hewgts[n];

					if (hgraph->where[k - 1] == 0) {
						hgraph->fs[k - 1] += hgraph->hewgts[n];
					}
				}
			}
			else if (count[0] > 1 && count[1] == 1) {	/*���и��{x,1}�ͳ���*/
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

	//����bnd����͸���queue����
	//bndptr[i] = 1 && moved[i] != -1   ->   �������������
	//bndptr[i] = 0  ->  �������������
	//bndptr[i] = 1 && moved[i] = 0     ->   �����������
	for (i = 0; i < incidence_hyperedge.size(); ++i) {
		n = incidence_hyperedge[i];

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {

			k = hgraph->eind[j];
			/*�ж϶����Ǳ߽綥��*/
			float sum = 0.0;
			for (ii = 0; ii < hgraph->incidence_hyperedge[k - 1].size(); ++ii) {
				nn = hgraph->incidence_hyperedge[k - 1][ii];
				sum += hgraph->hewgts[nn];
			}


			if (hgraph->te[k - 1] == sum) {

				/*������������г��߶�����te,��Ϊ�߽綥��*/
				if (hgraph->bndptr[k - 1] == 1) {
					/*�Ƿ���������������*/
					hgraph->bndlist.erase(k);
					hgraph->bndptr[k - 1] = 0;
					if (moved[k - 1] == -1)
						queue.rpq_delete(k);
				}
				else {

					/*֮ǰ�����������������*/
					hgraph->bndlist.insert(k);
					hgraph->bndptr[k - 1] = 1;
					queue.rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
				}

			}
			else {
				/*��������ĳ��߲�ȫ������te,��Ϊ�߽綥��*/
				if (hgraph->bndptr[k - 1] == 1) {
					if (moved[k - 1] == -1)
						queue.rpq_update(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
				}
				else {
					/*֮ǰ�����������������*/
					hgraph->bndlist.insert(k);
					hgraph->bndptr[k - 1] = 1;
					queue.rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
				}


			}


		}
	}


	hgraph->nbnd = hgraph->bndlist.size();

}
