#include "fm.h"
#include "checkhgraph.h"
#include "mystruct.h"

	/*������FM�㷨,���еĶ��������������������*/

	//#define VPQSTATUS_PRESENT      1       /* The vertex is in the queue */
	//#define VPQSTATUS_EXTRACTED    2       /* The vertex has been extracted from the queue */
#define VPQSTATUS_NOTPRESENT   3       /* The vertex is not present in the queue and*/
void fm::fm_2way_refine(s_hgraph * hgraph, vector<float> tpwgts2)
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

	vector<s_rpq> queue;

	vector<int> vstatus(hgraph->nvtxs, RPQ_NOTPRESENT);
	vector<int> moved(hgraph->nvtxs, -1);
	vector<int> swaps(hgraph->nvtxs, 0);

	//computing initial cutsize
	for (int pass = 0; pass < 20; pass++)
	{

		mincutorder = -1;
		mincut = newcut = initcut = hgraph->mincut;
		mindiff = abs(tpwgts[0] - hgraph->pwgts[0]);

		rpq_reset(queue);

		for (k = 1; k <= hgraph->nvtxs; ++k) {

			queue[hgraph->where[k - 1]].rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

			vstatus[k - 1] = RPQ_PRESENT;			/*��������뵽�����У������Ϊpresent*/

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

			/*���newcutС��mincut��p0����Ȩֵ��Ŀ��Ȩֵ������ ԭʼ�Ĳ�� + ƽ������Ȩֵ,�����mincut*/
			if (((newcut < mincut) && (abs(tpwgts[0] - hgraph->pwgts[0]) <= origdiff + avgvwgt)) ||
				((newcut == mincut) && (abs(tpwgts[0] - hgraph->pwgts[0]) < mindiff)) ||
				((newcut < mincut) && (hgraph->pwgts[0] < maxpartwgt[0]) && (hgraph->pwgts[0] > minpartwgt[0])) ||
				((newcut == mincut) && (hgraph->pwgts[0] < maxpartwgt[0]) && (hgraph->pwgts[0] > minpartwgt[0]))
				) {
				mincut = newcut;

				mindiff = abs(tpwgts[0] - hgraph->pwgts[0]);

				mincutorder = nswaps;
			}
			else if (nswaps - mincutorder > limit) { /*����15�ζ�û������Լ�������˳�����pass */

				newcut += (hgraph->fs[higain - 1] - hgraph->te[higain - 1]);

				hgraph->pwgts[from] += hgraph->vwgts[higain - 1];

				hgraph->pwgts[to] -= hgraph->vwgts[higain - 1];

				break;
			}


			hgraph->where[higain - 1] = to;

			moved[higain - 1] = nswaps;

			swaps[nswaps] = higain;

			vstatus[higain - 1] = RPQ_EXTRACTED;			/*��ʾ����higain�Ӷ�����ȡ��*/

			/*testcut = compute_mincut(hgraph);

			cout << "v = " << higain << " where " << from << " -> " << to << \
				" newcut = " << newcut << " acutalcut= " << testcut << endl;*/



			/*�������ƶ��Ķ���higain��ض����fs��te,�Լ�bnd���*/
			update_adjcent_gain_and_bnd_info(hgraph, queue, higain, vstatus);

			check_fs_and_te(hgraph);

			check_cutsign(hgraph);
		}

		/****************************************************************
		* Roll back computations //moved��ŵ��ƶ���˳�� ��nswap��������˳��
		*****************************************************************/
		for (i = 0; i < nswaps; i++)
		{
			moved[swaps[i] - 1] = -1;  /* reset moved array */
		}


		for (nswaps--; nswaps > mincutorder; nswaps--) {
			higain = swaps[nswaps];	//ȡ�����һ���ƶ���
			to = hgraph->where[higain - 1] = (hgraph->where[higain - 1] + 1) % 2;	//���higain��0�������򷵻ص�1

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

void fm::update_adjcent_gain_and_bnd_info(s_hgraph * hgraph, vector<s_rpq>& queue, int vertex,  vector<int> vstatus)
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

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) 
		{			/*����0������1�����������*/
			k = hgraph->eind[j];				
			count[hgraph->where[k - 1]]++;
		}

		//cout << " n = " << n << ": "<< count[0] << " " << count[0] << endl;


		if (hgraph->cutsign[n] == 1) 
		{	/*vertex��������ĳ����Ѿ���cut*/

			
			if (count[0] * count[1] == 0) 
			{		/*vertex�ƶ���,�ó���δ��cut*/

				hgraph->cutsign[n] = 0;		/*��ʱ�����ѱ��и�*/

				if (deg == 2)
				{	/*deg=2�ĳ��߱��и���*/ 
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) 
					{/*֮ǰ�ó��߶����ж����fs����Ϊ1,�ƶ���,���е�fs-1,�Ҵ�ʱ����δ��cut, ��te�Ĺ���Ϊ1*/
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						hgraph->fs[k - 1] -= hgraph->hewgts[n];
					}
				}
				else
				{	/*deg > 2�ĳ��߱��и���*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) 
					{/*��ô֮ǰ�ó��߶��ƶ��Ķ���vertex��fs����Ϊ1,�����ƶ�����Ҫ,��Сfs,�����ж���te����Ҫ���� */
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						if (k == vertex) {
							hgraph->fs[k - 1] -= hgraph->hewgts[n];
						}

					}

				}
			}
			else  
			{	/*vertex�ƶ���,�ó�����Ȼ��cut*/

				if (count[0] == 1) 
				{		/*���г�{1,x}��*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) 
					{
						/*��0�����Ķ���fs++*/
						k = hgraph->eind[j];

						if (hgraph->where[k - 1] == 0) 
						{
							hgraph->fs[k - 1] += hgraph->hewgts[n];
						}
						else 
						{	/*��1�����Ķ�����Ϊ2,���ƶ�ǰ��1�����Ķ�����Ϊ1,�����fs,����fs--*/
							if (count[1] == 2)			
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


			if (vstatus[k - 1] == RPQ_PRESENT) 
			{		/*������Ȼ�ڶ�����*/
				queue[hgraph->where[k - 1]].rpq_update(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);
			}


			//if (hgraph->te[k - 1] == sum) {

			//	if (vstatus[k - 1] == RPQ_PRESENT) {
			//		queue[hgraph->where[k - 1]].rpq_delete(k);
			//	}

			//}
			//else {
			//	/*��������ĳ��߲�ȫ������te,��Ϊ�߽綥��*/

			//	if (vstatus[k - 1] == RPQ_PRESENT) {

			//		queue[hgraph->where[k - 1]].rpq_update(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

			//	}
			//	else if (vstatus[k - 1] == RPQ_NOTPRESENT){/*֮ǰ�����������������*/

			//			queue[hgraph->where[k - 1]].rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

			//	}
			//}


		}
	}


}


/********************************************************
/*��������� ��ͼ���ݽṹ 
/*����ֵ��	��
/*���ܣ�		���¼��㳬ͼ������fs��te,���㳬�ߵ��и����
*********************************************************/
void fm::recompute_fm_params(s_hgraph *hgraph)
{
	int n = 0;
	int i = 0;
	int k = 0;
	int count0, count1;

	vector<int> bndptr(hgraph->nvtxs, 0);
	vector<float> fs(hgraph->nvtxs, 0.0);
	vector<float> te(hgraph->nvtxs, 0.0);

	hgraph->bndptr = bndptr;
	hgraph->fs = fs;
	hgraph->te = te;
	hgraph->bndlist.clear();

#if DEBUG

	int test_cut = compute_mincut(hgraph);

	if (test_cut != hgraph->mincut) {

		cout << "Error: The cutsize of hypergraph does'nt equal with compute cutsize [ " << \
			test_cut << " " << hgraph->mincut << " ]" << endl;

	}
#endif 


	/*���¼���cutsign*/
	for (n = 0; n < hgraph->nhedges; ++n) {

		hgraph->cutsign[n] = 0;

		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}


		std::sort(temp.begin(), temp.end());

		if (temp.front() != temp.back()) {

			hgraph->cutsign[n] = 1;					//������߱����ˣ���cutsign����Ϊ1

		}

	}

	/*���¼���fs��te*/
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

					if (hgraph->where[k - 1] == 0) {

						hgraph->fs[k - 1] += hgraph->hewgts[n];

					}

				}
			}

			if (count[1] == 1) {

				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					k = hgraph->eind[i];

					if (hgraph->where[k - 1] == 1) {

						hgraph->fs[k - 1] += hgraph->hewgts[n];

					}

				}

			}



		}
		else {
			//û���и�ĳ��ߣ����ж����te += ����Ȩֵ
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
				k = hgraph->eind[i];
				hgraph->te[k - 1] += hgraph->hewgts[n];
			}
		}

	}

}


void fm::bnd_insert(priority_queue <gain_node> &queue, float gain, int vertex)
{

	gain_node *n = new gain_node;

	n->gain = gain;
	n->number = vertex;

	queue.push(*n);
}




int fm::rpq_get_top(priority_queue <gain_node> &queue)
{
	int number = queue.top().number;
	queue.pop();

	if (queue.empty()) {
		return -1;
	}

	return number;
}


void fm::rpq_reset(vector<s_rpq >& queue)
{


	s_rpq tmp;
	s_rpq q0;
	s_rpq q1;

	while (queue.size() != 0) {
		queue.pop_back();
	}


	queue.push_back(q0);
	queue.push_back(q1);

}


//����where�������Ϣ,���㳬ͼhgraph�����˶��ٱ߸�, ���ұ�Ǳ��и�ĳ���
int fm::compute_cutsize_and_cutsign(s_hgraph * hgraph)
{
	int  i, n, k;
	int cut = 0;
	int nhedges = hgraph->nhedges;

	vector<int> cutsign(nhedges, 0);
	hgraph->cutsign = cutsign;



	for (n = 0; n < nhedges; ++n) {

		vector<int> tempwhere;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

			k = hgraph->eind[i];

			tempwhere.push_back(hgraph->where[k - 1]);

		}


		sort(tempwhere.begin(), tempwhere.end());

		if (tempwhere.front() != tempwhere.back()) {

			hgraph->cutsign[n] = 1;					//������߱����ˣ���cutsign����Ϊ1 

			cut += hgraph->hewgts[n];

		}

	}


	return cut;

}




void fm::modified_cut_sign(s_hgraph *hgraph)
{
	int n, i, k;
	for (n = 0; n < hgraph->nhedges; ++n)
	{
		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}


		std::sort(temp.begin(), temp.end());

		if (temp.front() != temp.back()) {

			if (hgraph->cutsign[n] == 0)
			{
				cout << "cutsing error: 1" << endl;
				hgraph->cutsign[n] == 1;
			}
			continue;
		}


		if (hgraph->cutsign[n] == 1)
		{
			cout << "cutsing error: 0" << endl;
			hgraph->cutsign[n] == 0;
		}


	}
}
