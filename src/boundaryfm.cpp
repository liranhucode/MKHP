#include "boundaryfm.h"
#include "checkhgraph.h"

void boundaryfm::boundary_fm_2way_refine(s_hgraph * hgraph, vector<float> tpwgts2)
{

	int n, i, j, k, gain, from, to, higain, mincutorder, nswaps;
	int avgvwgt, origdiff, mindiff, limit;
	int maxpartwgt[2];
	int minpartwgt[2];
	vector<float> tpwgts(2, 0);

	tpwgts[0] = hgraph->tvwgts * tpwgts2[0];
	tpwgts[1] = hgraph->tvwgts - tpwgts[0];

	/*���������ֵ����Сֵ*/
	maxpartwgt[0] = tpwgts[0] * 1.05;
	maxpartwgt[1] = tpwgts[1] * 1.05;
	minpartwgt[0] = tpwgts[0] * 0.95;
	minpartwgt[1] = tpwgts[1] * 0.95;

	avgvwgt = hgraph->tvwgts / hgraph->nvtxs;
	//avgvwgt = min( (hgraph->pwgts[0] + hgraph->pwgts[1]) / 20, 2 * (hgraph->pwgts[0] + hgraph->pwgts[1]) / hgraph->nvtxs );

	origdiff = abs(tpwgts[0] - hgraph->pwgts[0]);

	int newcut, initcut, mincut;
	int testcut;

	limit = (int)min(max(0.01*hgraph->nvtxs, 15.0), 100.0);

	vector<s_rpq> queue;


	vector<int> moved(hgraph->nvtxs, -1);
	vector<int> swaps(hgraph->nvtxs, 0);

	//computing initial cutsize
	for (int pass = 0; pass < 10; pass++)
	{

		mincutorder = -1;
		mincut = newcut = initcut = hgraph->mincut;
		mindiff = abs(tpwgts[0] - hgraph->pwgts[0]);

		rpq_reset(queue);

		/*���߽綥����������,���������С��������*/
		for (unordered_set<int> ::iterator it = hgraph->bndlist.begin(); it != hgraph->bndlist.end(); it++) {
			i = *it;
			queue[hgraph->where[i - 1]].rpq_insert(i, hgraph->fs[i - 1] - hgraph->te[i - 1]);
		}



		for (nswaps = 0; nswaps < hgraph->nvtxs; nswaps++)
		{

			/*ȷ����ʼ�ƶ�����*/
			from = (tpwgts[0] - hgraph->pwgts[0] < tpwgts[1] - hgraph->pwgts[1] ? 0 : 1);		
			to = (from + 1) % 2;


			if ((higain = queue[from].rpq_get_top()) == -1 ) {
				break;
			}

			//int gain = check_fs_and_te(hgraph, higain);

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
			else if (nswaps - mincutorder > limit) {
				/*����1limit�ζ�û������Լ�������˳�����pass */

				newcut += (hgraph->fs[higain - 1] - hgraph->te[higain - 1]);
				hgraph->pwgts[from] += hgraph->vwgts[higain - 1];
				hgraph->pwgts[to] -= hgraph->vwgts[higain - 1];
				break;
			}


			hgraph->where[higain - 1] = to;
			moved[higain - 1] = nswaps;
			swaps[nswaps] = higain;


			/*testcut = compute_mincut(hgraph);
			cout << "v = " << higain << " where " << from << " -> " << to << \
							" newcut = " << newcut << " acutalcut= "<<testcut  << endl;*/


			/*�ƶ�higain������ٽ����������*/
			update_adjcent_gain_and_bnd_info(hgraph, queue, higain, moved);

			//check_fs_and_te(hgraph);

			//check_cutsign(hgraph);
		}


		/****************************************************************
		* ���ؼ���,�˻���������ƶ�
		*****************************************************************/
		for (i = 0; i < nswaps; i++)		//moved��ŵ��ƶ���˳�� ��nswap��������˳��
			moved[swaps[i] - 1] = -1;  /* reset moved array */

		for (nswaps--; nswaps > mincutorder; nswaps--) {
			higain = swaps[nswaps];	
			to = hgraph->where[higain - 1] = (hgraph->where[higain - 1] + 1) % 2;	
			hgraph->pwgts[to] += hgraph->vwgts[higain - 1];
			hgraph->pwgts[(to + 1) % 2] -= hgraph->vwgts[higain - 1];

		}

		hgraph->mincut = mincut;
		recompute_fm_params(hgraph);

		if (mincutorder <= 0 || mincut == initcut) {
#if DEBUG
			testcut = compute_mincut(hgraph);
			cout << "testcut = " << testcut << " mincut = " << mincut << endl;
#endif
			break;
		}




	}
}






/*����vertex���㣬 �ٽ������fs��te,������*/
void boundaryfm::update_adjcent_gain_and_bnd_info(s_hgraph *hgraph, vector<s_rpq> &queue, int vertex, vector<int> moved)
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

				if (deg == 2) {	

					/*ԭ��deg=2�ĳ��߱��и���,��ô֮ǰ�ó��߶����ж����fs����Ϊ1,����û�в����и�,��fs--,��ʱte++*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						hgraph->fs[k - 1] -= hgraph->hewgts[n];
					}
				}
				else {
					/*ԭ��deg > 2�ĳ��߱��и��ˣ���ô֮ǰ�ó��߶��Ѿ��ƶ��Ķ���vertex��fs����Ϊ1,
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

	/*����bnd����͸���queue����*/
	for (i = 0; i < incidence_hyperedge.size(); ++i) 
	{	

		n = incidence_hyperedge[i];

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) 
		{
			k = hgraph->eind[j];

			/*�ж϶����Ǳ߽綥��*/
			float sum = 0.0;
			for (ii = 0; ii < hgraph->incidence_hyperedge[k - 1].size(); ++ii) 
			{
				nn = hgraph->incidence_hyperedge[k-1][ii];
				sum += hgraph->hewgts[nn];
			}

			//����������У�
			//			bndptr[i] = 1 && moved[i] = -1   
			//������������У�
			//			bndptr[i] = 0  
			//			bndptr[i] = 1 && moved[i] = 0     

			if (hgraph->bndptr[k - 1] == 1)
			{	/*֮ǰ�Ǳ߽綥��*/

				if (hgraph->te[k - 1] == sum)
				{	/*�������ڲ�����,�Ǳ߽綥��*/

					hgraph->bndlist.erase(k);
					hgraph->bndptr[k - 1] = 0;

					if (moved[k - 1] == NOT_MOVED)	/*֮ǰ���ڶ�����,��ɾ��*/
						queue[hgraph->where[k - 1]].rpq_delete(k);

				}
				else 
				{	/*������Ȼ�Ǳ߽綥�㶥��*/

					if (moved[k - 1] == NOT_MOVED)	/*δ�ƶ����Ķ���һ�����ڶ����У�����*/
						queue[hgraph->where[k - 1]].rpq_update(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

				}

			}
			else 
			{	/*֮ǰ���ڲ�����*/

				if (hgraph->te[k - 1] != sum)
				{	/*���ڱ���Ǳ߽綥��*/

					hgraph->bndlist.insert(k);
					hgraph->bndptr[k - 1] = 1;
					queue[hgraph->where[k - 1]].rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

				}
				else
				{	/*������Ȼ���ڲ�����*/

				}

			}


		}
	}


	hgraph->nbnd = hgraph->bndlist.size();

}



/********************************************************
/*��������� ��ͼ���ݽṹ
/*����ֵ��	��
/*���ܣ�		���¼��㳬ͼ������fs��te,���㳬�ߵ��и�����Լ�
/*			�߽綥��
*********************************************************/
void boundaryfm::recompute_fm_params(s_hgraph *hgraph)
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




	for (n = 0; n < hgraph->nhedges; ++n)
	{
		if (hgraph->cutsign[n] == 1) {

			vector<int> count(2, 0);

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];


				count[hgraph->where[k - 1]]++;

				hgraph->bndlist.insert(k);			//��cut�ĳ����������ж�����Ǳ߽綥��
				hgraph->bndptr[k - 1] = 1;

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
			//û���и�ĳ��ߣ����ж����te += ����Ȩֵ
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
				k = hgraph->eind[i];
				hgraph->te[k - 1] += hgraph->hewgts[n];
			}
		}

	}


	hgraph->nbnd = hgraph->bndlist.size();
}


void boundaryfm::bnd_insert(priority_queue <gain_node> &queue, float gain, int vertex)
{

	gain_node *n = new gain_node;

	n->gain = gain;
	n->number = vertex;

	queue.push(*n);
}




int boundaryfm::rpq_get_top(priority_queue <gain_node> &queue)
{
	int number = queue.top().number;
	queue.pop();

	if (queue.empty()) {
		return -1;
	}

	return number;
}





void boundaryfm::rpq_reset(vector<s_rpq >& queue)
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
int boundaryfm::compute_cutsize_and_cutsign(s_hgraph * hgraph)
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




void boundaryfm::modified_cut_sign(s_hgraph *hgraph)
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

/*���ݶ��е��������,���ݷ�����Ȩֵ,ȷ��from��to��ֵ*/
void boundaryfm::make_where_to_move(vector<s_rpq> &queue, vector<float> tpwgts, vector<int> pwgts, int &from, int &to)
{

	vector<float> max_vwgts_part(2,0.0);

	max_vwgts_part[0] = tpwgts[0] * 1.05;
	max_vwgts_part[1] = tpwgts[1] * 1.05;


	if (queue[0].rpq_top() > queue[1].rpq_top() ) {

		from = 0;
		to = 1;
	}
	else {
		from = 1;
		to = 0;
	}


	//from = (max_vwgts_part[0] - pwgts[0] < max_vwgts_part[1] - pwgts[1] ? 0 : 1);
	//to = (from + 1) % 2;


}

