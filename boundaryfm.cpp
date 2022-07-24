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

	/*分区的最大值和最小值*/
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

		/*将边界顶点插入队列中,按照增益大小进行排序*/
		for (unordered_set<int> ::iterator it = hgraph->bndlist.begin(); it != hgraph->bndlist.end(); it++) {
			i = *it;
			queue[hgraph->where[i - 1]].rpq_insert(i, hgraph->fs[i - 1] - hgraph->te[i - 1]);
		}



		for (nswaps = 0; nswaps < hgraph->nvtxs; nswaps++)
		{

			/*确定起始移动分区*/
			from = (tpwgts[0] - hgraph->pwgts[0] < tpwgts[1] - hgraph->pwgts[1] ? 0 : 1);		
			to = (from + 1) % 2;


			if ((higain = queue[from].rpq_get_top()) == -1 ) {
				break;
			}

			//int gain = check_fs_and_te(hgraph, higain);

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
			else if (nswaps - mincutorder > limit) {
				/*连续1limit次都没有满足约束，则退出本次pass */

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


			/*移动higain顶点的临近顶点的增益*/
			update_adjcent_gain_and_bnd_info(hgraph, queue, higain, moved);

			//check_fs_and_te(hgraph);

			//check_cutsign(hgraph);
		}


		/****************************************************************
		* 返回计算,退回无增益的移动
		*****************************************************************/
		for (i = 0; i < nswaps; i++)		//moved存放的移动的顺序 ，nswap代表交换的顺序
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






/*更新vertex顶点， 临近顶点的fs和te,并计算*/
void boundaryfm::update_adjcent_gain_and_bnd_info(s_hgraph *hgraph, vector<s_rpq> &queue, int vertex, vector<int> moved)
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

				if (deg == 2) {	

					/*原先deg=2的超边被切割了,那么之前该超边对所有顶点的fs贡献为1,现在没有产生切割,则fs--,此时te++*/
					for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) {
						k = hgraph->eind[j];
						hgraph->te[k - 1] += hgraph->hewgts[n];
						hgraph->fs[k - 1] -= hgraph->hewgts[n];
					}
				}
				else {
					/*原先deg > 2的超边被切割了，那么之前该超边对已经移动的顶点vertex的fs贡献为1,
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

	/*更新bnd数组和更新queue队列*/
	for (i = 0; i < incidence_hyperedge.size(); ++i) 
	{	

		n = incidence_hyperedge[i];

		for (j = hgraph->eptr[n]; j < hgraph->eptr[n + 1]; ++j) 
		{
			k = hgraph->eind[j];

			/*判断顶点是边界顶点*/
			float sum = 0.0;
			for (ii = 0; ii < hgraph->incidence_hyperedge[k - 1].size(); ++ii) 
			{
				nn = hgraph->incidence_hyperedge[k-1][ii];
				sum += hgraph->hewgts[nn];
			}

			//在增益队列中：
			//			bndptr[i] = 1 && moved[i] = -1   
			//不在增益队列中：
			//			bndptr[i] = 0  
			//			bndptr[i] = 1 && moved[i] = 0     

			if (hgraph->bndptr[k - 1] == 1)
			{	/*之前是边界顶点*/

				if (hgraph->te[k - 1] == sum)
				{	/*现在是内部顶点,非边界顶点*/

					hgraph->bndlist.erase(k);
					hgraph->bndptr[k - 1] = 0;

					if (moved[k - 1] == NOT_MOVED)	/*之前存在队列中,则删除*/
						queue[hgraph->where[k - 1]].rpq_delete(k);

				}
				else 
				{	/*现在仍然是边界顶点顶点*/

					if (moved[k - 1] == NOT_MOVED)	/*未移动过的顶点一定存在队列中，更新*/
						queue[hgraph->where[k - 1]].rpq_update(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

				}

			}
			else 
			{	/*之前是内部顶点*/

				if (hgraph->te[k - 1] != sum)
				{	/*现在变成是边界顶点*/

					hgraph->bndlist.insert(k);
					hgraph->bndptr[k - 1] = 1;
					queue[hgraph->where[k - 1]].rpq_insert(k, hgraph->fs[k - 1] - hgraph->te[k - 1]);

				}
				else
				{	/*现在仍然是内部顶点*/

				}

			}


		}
	}


	hgraph->nbnd = hgraph->bndlist.size();

}



/********************************************************
/*输入参数： 超图数据结构
/*返回值：	空
/*功能：		重新计算超图分区的fs和te,计算超边的切割情况以及
/*			边界顶点
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

	/*重新计算cutsign*/
	for (n = 0; n < hgraph->nhedges; ++n) {

		hgraph->cutsign[n] = 0;

		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}


		std::sort(temp.begin(), temp.end());

		if (temp.front() != temp.back()) {
			hgraph->cutsign[n] = 1;					//如果超边被切了，则cutsign设置为1 
		}

	}




	for (n = 0; n < hgraph->nhedges; ++n)
	{
		if (hgraph->cutsign[n] == 1) {

			vector<int> count(2, 0);

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				k = hgraph->eind[i];


				count[hgraph->where[k - 1]]++;

				hgraph->bndlist.insert(k);			//被cut的超边其中所有顶点就是边界顶点
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
			//没被切割的超边，所有顶点的te += 超边权值
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


//根据where数组的信息,计算超图hgraph产生了多少边割, 并且标记被切割的超边
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

			hgraph->cutsign[n] = 1;					//如果超边被切了，则cutsign设置为1 

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

/*根据队列的增益情况,根据分区的权值,确定from和to的值*/
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

