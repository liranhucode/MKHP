#include <unordered_set>
#include <iomanip>
#include <unordered_map>
#include "coarsen.h"
#include "mystruct.h"
#include "myutil.h"
#include "checkhgraph.h"


/********************************************************
/*输入参数： 超图数据结构, 分区个数
/*返回值：	粗化后的超图
/*功能：		将大的超图变成小的超图
*********************************************************/
s_hgraph * s_coarsen::coarsen_hgraph(s_hgraph *original_hgraph, int npart)
{
	//cout << endl;

	int level = 0;
	double maxvwgt = 0.0;
	s_hgraph *hgraph = original_hgraph;

	/*设置最粗超图顶点数*/
	coarsen_to = npart * 30;

	/*设置允许顶点匹配的最大权值*/
	maxvwgt = ( 5 * hgraph->tvwgts ) / coarsen_to;

	do {

		/*给cmap数组分配空间*/
		if (hgraph->cmap.size() == 0)
		{  
			vector<int> cmap(hgraph->nvtxs, -1);
			hgraph->cmap = cmap;
		}

		/*以超边为单位进行匹配*/
		hyperedge_match(hgraph, maxvwgt);

		level++;

		hgraph = hgraph->coarser;

#if DEBUG
		check_hgraph(hgraph);
#endif

	} while (hgraph->nvtxs > 200
		&&hgraph->nvtxs < 0.85 * hgraph->finer->nvtxs
		/*&&hgraph->nhedges < 0.85 * hgraph->finer->nhedges*/
		/*&&hgraph->nhedges < hgraph->nvtxs / 2*/
		);




#if 1
	s_hgraph *tmp = original_hgraph;
	for (int i = 0; i <= level; i++){
		cout <<" [ "<< setw(5) << tmp->nvtxs << "  " << setw(5) << tmp->nhedges << "  " << setw(5) << setiosflags(ios::fixed)<< setprecision(2) << tmp->average_vertex_deg() << \
			"  " << setw(5) <<  tmp->average_hyperedge_deg() << "  " << setw(5) <<  (float)tmp->average_vertex_deg()/ tmp->average_hyperedge_deg() << " ] " << endl;
		tmp = tmp->coarser;
	}

#endif

	cout << endl;
	return hgraph;
}



/********************************************************
/*输入参数： 超图数据结构 顶点最大的权值约束
/*返回值：	空
/*功能：		将超图中的顶点进行匹配
*********************************************************/
void s_coarsen::hyperedge_match(s_hgraph *hgraph, double maxvwgt)
{

	int k = 0;
	int i = 0; 
	int ii = 0; 
	int n = 0;
	int num = 0;
	int vertex = 0;
	int cnvtxs = 0;
	int vertex_weight_sum = 0;

	int already_match_to = 0;
	int istart, iend;
	int all_vertices_match = 1;
	int all_vertices_no_match = 0;
	int count = 0;
	int nvtxs = hgraph->nvtxs;
	int nhedges = hgraph->nhedges;


	float average_vertex_deg = hgraph->average_vertex_deg();

	/*所有顶点未匹配, 所有超边未收缩*/
	vector<int> hyperedge_contracted(nhedges, NO_CONTRACTED);
	vector<int> matched(nvtxs, UNMATCHED);		

	/*生成超边排序队列*/
	vector<hyperedge> queue;
	creat_rpq_for_hyperedge(hgraph, queue);				

 	for (n = 0; n < nhedges; ++n)
	{



		vertex_weight_sum = 0;
		all_vertices_no_match = YES;

		/*依次遍历队列中已排序的超边*/
		num = queue[n].number;
		
		if (hgraph->eptr[num + 1] - hgraph->eptr[num] > 20)
			continue;

		/*判断超边的类型*/
		if (hgraph->is_single_hyperedge(num))		
		{
			/*单个顶点的超边*/
			hyperedge_contracted[num] = DISCARD;
			vertex = hgraph->eind[hgraph->eptr[num]];
			//cout << "single edge: " << vertex << endl;

			/*如果该顶点已经被匹配过,则跳过*/
			if (matched[vertex - 1] == UNMATCHED) {			
				/*将独立超边中的顶点在cmap数组中标记成-1,表示该顶点不会copy到下一级超图*/
				hgraph->cmap[vertex - 1] = -1;
				matched[vertex - 1] = DISCARD;
			}

		}
		else if ( hgraph->is_dependent_hyperedge(num) ) 
		{		
			//cout << "独立超边： " << num << endl;
			/*独立超边, 所有顶点度为1的超边*/
			hyperedge_contracted[num] = DISCARD;	

			for ( i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i ) {
				/*将独立超边中的顶点在cmap数组中标记成-1,表示该顶点不会copy到下一级超图*/
				vertex = hgraph->eind[i];
				if (matched[vertex - 1] == UNMATCHED) {

					hgraph->cmap[vertex - 1] = -1;
					matched[vertex - 1] = DISCARD;
				}
			}
		}
		else {

			/*判断该超边是否未匹配过*/
			for ( i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i ) {
				vertex = hgraph->eind[i];
				if (matched[vertex-1] == MATCHED) {						
					all_vertices_no_match = NO;	
					break;
				}
			}

			if (all_vertices_no_match == YES) {		/*该超边中所有顶点均未匹配*/

				for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; i++) {
					vertex = hgraph->eind[i];
					vertex_weight_sum += hgraph->vwgts[vertex - 1];
				}


				if (vertex_weight_sum <= maxvwgt) {			
					/*满足顶点权值约束,该超边可以进行匹配*/
					cnvtxs++;
					hyperedge_contracted[num] = CONTRACTED;

					for (ii = hgraph->eptr[num]; ii < hgraph->eptr[num + 1]; ii++) {

						vertex = hgraph->eind[ii];
						hgraph->cmap[vertex - 1] = cnvtxs;		/*顶点序号从1开始*/
						matched[vertex - 1] = MATCHED;
					}
				}			
			}  
		}

	}


	/*对未进行匹配的超边进行重新审查,进一步寻找可以进行匹配的顶点*/
	for (n = 0; n < nhedges; ++n)
	{
		if (hyperedge_contracted[n] == NO_CONTRACTED)
		{
			count = 0;
			vertex_weight_sum = 0;

			/*计算未匹配的顶点权值之和*/
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				vertex = hgraph->eind[i];
				if (matched[vertex - 1] == UNMATCHED) {
					vertex_weight_sum += hgraph->vwgts[vertex - 1];
				}
			}

			if (vertex_weight_sum <= maxvwgt)
			{
				/*权值符合约束,将该超边中未匹配的顶点匹配到一起*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {

						if (count == 0)
						{
							cnvtxs++;			/*顶点只增加一次*/
							count++;
						}
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}
			else
			{
				/*权值不符合约束,将该超边中未匹配的顶点匹配到一起*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {
														/*无法进一步匹配的顶点直接copy到下一级*/
						cnvtxs++;
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}
			/*统计该超边中的可以进行匹配的顶点个数 */
			/*for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				vertex = hgraph->eind[i];
				if (matched[vertex - 1] == UNMATCHED) {

					vertex_weight_sum += hgraph->vwgts[vertex - 1];
					if (vertex_weight_sum <= maxvwgt) {	
						count++;				
						max_len = i;
					}
				}
			}

			if (count < 2) {

				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i)
				{
					vertex = hgraph->eind[i];

					if (matched[vertex - 1] == UNMATCHED) {
						cnvtxs++;
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}
			else {

				cnvtxs++;

				for (i = hgraph->eptr[n]; i <= max_len; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {

						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}*/
		}
		
			
	}


	//*检测cnvtxs顶点是否正确*/
	//vector<int> check_map(nvtxs, 0);
	//vector<int> check_map(nvtxs, 0);
	//int count_of_vertices = 0;
	//for (int i = 0; i < nvtxs; ++i)
	//{
	//	if(hgraph->cmap[i] != -1 )
	//		check_map[hgraph->cmap[i]]++;
	//}
	//for (int i = 0; i < nvtxs; ++i)
	//{
	//	if (check_map[i] != 0)
	//		count_of_vertices++;
	//}


	//test
	for (int i = 0; i < nvtxs; ++i)
	{
		if (matched[i] == UNMATCHED)
		{
			cout << " error: unmatched!" << endl;
		}

	}

	//if (count_of_vertices != cnvtxs)
	//	cout << "粗化顶点错误： " << count_of_vertices << " " << cnvtxs << endl;

	/*根据cnvtxs、cmap数组和hyperedge_contracted数组生成更粗一级超图*/
	creat_coarser_hgraph(hgraph, cnvtxs, hyperedge_contracted);


}


/*遍历顶点,寻找其关联的超边中可以直接进行contract的超边*/
void s_coarsen::random_match(s_hgraph * hgraph, double maxvwgt)
{

}




/********************************************************
/*输入参数： 超图数据结构 更粗一级超图的顶点数 超边收缩标志数组
/*返回值：	空
/*功能：		生成更粗一级超图,地址保存至hgraph->coarser中
*********************************************************/
void s_coarsen::creat_coarser_hgraph(s_hgraph *& hgraph, const int &cnvtxs, vector<int>& hyperedge_contracted)
{
	int n;
	int i, ii;
	int vertex;
	int nvtxs, nhedges;
	int cnhedges = 0;
	int base = 0;
	int result = 0;
	int l = 0;
	int find;

	/*new coarser hgraph*/
	s_hgraph *chgraph = new s_hgraph();


	nvtxs = hgraph->nvtxs;
	nhedges = hgraph->nhedges;

	/* store the new hgraph*/
	vector<int> ceptr;
	vector<int> ceind;


	vector<size_s> chewgts;
	vector<size_s> cvwgts(cnvtxs, 0);



	unordered_map<string,int> hyperedge_table;		/*键值：超边的string格式  值：超边权值*/
	unordered_map<string, int>::iterator hyperedge_table_it;
	string hyperedge_string;



	for (n = 0; n < nhedges; n++)
	{
		vector<int> temp;
		int cnt = 0;


		if (hyperedge_contracted[n] == NO_CONTRACTED) {

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; i++){

				vertex = hgraph->eind[i];

				if(hgraph->cmap[vertex - 1] != -1)
					temp.push_back(hgraph->cmap[vertex-1]);

			}


			/*删除重复的顶点*/
			sort(temp.begin(), temp.end());

			temp.erase(unique(temp.begin(), temp.end()), temp.end());

			if (temp.size() < 2){
				
				/*修正hyperedge_contracted数组*/
				hyperedge_contracted[n] = CONTRACTED;
			}
			else{

				/*将超点转换成字符串型便于存储, 如temp = { 0, 1, 2} 则 string型 = "0,1,2" */
				hyperedge_string = vector_to_string(temp);

				if ((hyperedge_table_it = hyperedge_table.find(hyperedge_string)) != hyperedge_table.end()) {

					hyperedge_table_it->second += hgraph->hewgts[n];		/*重复超边权值+1*/
				}
				else {
					/*如果没找到,则将原超边的权值传递下去*/
					hyperedge_table.insert(make_pair(hyperedge_string, hgraph->hewgts[n]));	
				}

			}

		}
	}


	map_convert_to_hyperedge(hyperedge_table, chgraph->eind, chgraph->eptr, chewgts, cnhedges);

	/*设置cvwgts数组*/
	for (i = 1; i <= nvtxs; i++) {

		if(hgraph->cmap[i - 1] == -1){
			continue;
		}
		cvwgts[hgraph->cmap[i - 1] - 1] += hgraph->vwgts[i - 1];
	}

	/*检查cvwgts数组*/
	//for (i = 0; i < cnvtxs; ++i)
	//{
	//	if (cvwgts[i] == 0)
	//		cout << "权值小于0" << endl;
	//}


	/*设置粗化的超图chgraph*/
	chgraph->vwgts = cvwgts;
	chgraph->hewgts = chewgts;
	chgraph->nvtxs = cnvtxs;
	chgraph->nhedges = cnhedges;
	hgraph->coarser = chgraph;
	chgraph->finer = hgraph;



	chgraph->compute_hyperedge_degree();		/*计算超边的度*/
	chgraph->compute_total_weight();		/*计算顶点权值之和*/
	chgraph->compute_vertex_degree();		/*计算顶点度*/
	chgraph->compute_incidence_hyperedge();	/*计算顶点关联超边*/
	//chgraph->compute_adjcent_vertices();
}




/*顶点匹配粗化算法*/
s_hgraph * s_coarsen::CoarsenHypergraph(s_hgraph * original_hgraph, int npart)
{

	int level = 0;
	double maxvwgt = 0.0;
	s_hgraph *hgraph = original_hgraph;

	coarsen_to = 100;
	

	do {

		VertexMatching( hgraph );

		level++;

		hgraph = hgraph->coarser;

#if DEBUG
		check_hgraph(hgraph);
#endif 

	} while (hgraph->nvtxs > 100);




#if 1
	s_hgraph *tmp = original_hgraph;
	for (int i = 0; i <= level; i++) {
		cout << " [ " << setw(5) << tmp->nvtxs << "  " << setw(5) << tmp->nhedges << "  " << setw(5) << setiosflags(ios::fixed) << setprecision(2) << tmp->average_vertex_deg() << \
			"  " << setw(5) << tmp->average_hyperedge_deg() << "  " << setw(5) << (float)tmp->average_vertex_deg() / tmp->average_hyperedge_deg() << " ] " << endl;
		tmp = tmp->coarser;
	}

#endif

	cout << endl;
	return hgraph;

	return nullptr;
}


int s_coarsen::findCandicate(s_hgraph* &hypergraph, int &u, vector<int> &unmatched)
{
	int cand = -1;
	int v = 0;
	int i = 0, j = 0;

	unordered_set <int> adj = hypergraph->adjcent_vertices[u-1];

	vector<int> incidence = hypergraph->incidence_hyperedge[u-1];

	float uweight = hypergraph->vwgts[u - 1];

	double conn = 0, maxconn = 0;

	float weight = 0;

	//遍历临近顶点
	for (auto it = adj.begin(); it != adj.end(); ++it)
	{

		v = *it;
		weight = 0;

		if (unmatched[v - 1] == UNMATCHED)		//该邻近顶点未匹配
		{
			vector<int> adj_incidence = hypergraph->incidence_hyperedge[v - 1];

			//计算vertex与邻近顶点*it公共超边的degree之和
			for (i = 0; i < incidence.size(); ++i){
				
				int n = incidence[i];
				vector<int>::iterator result = find(adj_incidence.begin(), adj_incidence.end(), n); //查找n
				//找到
				if (result != adj_incidence.end())
					weight += (float)1 / hypergraph->hyperedge_degree[n];
			}

			if (weight == 0) {
				cout << "Error!" << endl;
				exit(0);
			}

			//寻找conn最大的顶点
			conn = uweight * hypergraph->vwgts[*it - 1] * weight;

			if (conn > maxconn){
				maxconn = conn;
				cand = v;
			}
			
		}
		

	}

	return cand;


}

void s_coarsen::constructNewCoarserHypergraph(s_hgraph *& hypergraph, int cnvtxs)
{
	int n;
	int i, ii;
	int v;
	int nvtxs, nhedges;
	int cnhedges = 0;
	int base = 0;
	int result = 0;
	int l = 0;
	int find;

	/*new coarser hgraph*/
	s_hgraph *chgraph = new s_hgraph();


	nvtxs = hypergraph->nvtxs;
	nhedges = hypergraph->nhedges;

	/* store the new hgraph*/
	vector<int> ceptr;
	vector<int> ceind;


	vector<size_s> chewgts;
	vector<size_s> cvwgts(cnvtxs, 0);



	unordered_map<string, int> hyperedge_table;		/*键值：超边的string格式  值：超边权值*/
	unordered_map<string, int>::iterator hyperedge_table_it;
	string hyperedge_string;



	for (n = 0; n < nhedges; n++)
	{
		vector<int> temp;
		int cnt = 0;


		for (i = hypergraph->eptr[n]; i < hypergraph->eptr[n + 1]; ++i )
		{
			v = hypergraph->eind[i];
			temp.push_back(hypergraph->cmap[v - 1]);
		}


		/*删除重复的顶点*/
		sort(temp.begin(), temp.end());

		temp.erase(unique(temp.begin(), temp.end()), temp.end());

		if (temp.size() < 2) {
			continue;
		}
		else {

			/*将超点转换成字符串型便于存储, 如temp = { 0, 1, 2} 则 string型 = "0,1,2" */
			hyperedge_string = vector_to_string(temp);

			if ((hyperedge_table_it = hyperedge_table.find(hyperedge_string)) != hyperedge_table.end()) {

				hyperedge_table_it->second += hypergraph->hewgts[n];		/*重复超边权值+1*/
			}
			else {
				/*如果没找到,则将原超边的权值传递下去*/
				hyperedge_table.insert(make_pair(hyperedge_string, hypergraph->hewgts[n]));
			}

		}

	}


	map_convert_to_hyperedge(hyperedge_table, chgraph->eind, chgraph->eptr, chewgts, cnhedges);

	/*设置cvwgts数组*/
	for (i = 1; i <= nvtxs; i++) 
	{
		if (hypergraph->cmap[i - 1] == -1) {
			continue;
		}
		cvwgts[hypergraph->cmap[i - 1] - 1] += hypergraph->vwgts[i - 1];
	}



	/*设置粗化的超图chgraph*/
	chgraph->vwgts = cvwgts;
	chgraph->hewgts = chewgts;
	chgraph->nvtxs = cnvtxs;
	chgraph->nhedges = cnhedges;
	hypergraph->coarser = chgraph;
	chgraph->finer = hypergraph;



	chgraph->compute_hyperedge_degree();		/*计算超边的度*/
	chgraph->compute_total_weight();		/*计算顶点权值之和*/
	chgraph->compute_vertex_degree();		/*计算顶点度*/
	chgraph->compute_incidence_hyperedge();	/*计算顶点关联超边*/

	chgraph->compute_adjcent_vertices();
}

void s_coarsen::VertexMatching(s_hgraph* &hypergraph)
{

	int nvtxs = hypergraph->nvtxs;

	vector<int> numPermute;
	rand_permute( numPermute, nvtxs, 0);

	vector<int> unmatched(nvtxs, UNMATCHED);

	vector<int> cmap(nvtxs, -1);

	int nMatched = 0;
	int k = 0;
	int i = 0;
	int u, v;

	while( i < nvtxs &&  (float) nMatched/nvtxs < 1)
	{
		u = numPermute[i];

		if (unmatched[u-1] == UNMATCHED)
		{

			cmap[u - 1] = ++k;

			v = findCandicate(hypergraph, u, unmatched);

			if ( v != -1 ){		/*找到v*/
				cmap[v - 1] = k;
				nMatched += 2;
				unmatched[v - 1] = MATCHED;
			}

			unmatched[u - 1] = MATCHED;
		}

		i++;
	}

	while (i < nvtxs)
	{
		u = numPermute[i];

		if (unmatched[u -1] == UNMATCHED)
		{
			cmap[u - 1] = ++k;
			unmatched[u - 1] = MATCHED;
		}
		i++;
	}

	hypergraph->cmap = cmap;
	
	constructNewCoarserHypergraph(hypergraph, k);


}
















//=================================================================================//
//push weight of edge into priorty queue
//=================================================================================//
void s_coarsen::creat_rpq_for_hyperedge(s_hgraph * hgraph, vector<hyperedge>& queue)
{

	int n;
	int nhedges = hgraph->nhedges;
	int total_weight;
	vector<int> perm;

	rand_permute(perm, nhedges, 0);
	

	//rand_permute(perm, nhedges, 0);
	int k = 0;

	for (n = 0; n < nhedges; n++)
	{
		hyperedge edge;
		//total_weight = 0;

		/*for (int i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
		k = hgraph->eind[i];
		total_weight += hgraph->vwgts[k - 1];
		}*/

		edge.number = n;	/*顶点序号*/

		edge.weight = hgraph->hewgts[n];		/*超边权值*/

												//edge.total_vertex_weigth = total_weight;			/*超边中顶点权值之和*/

		edge.degree = hgraph->hyperedge_degree[n];		/*超边度*/

		//edge.degree = k++;											//edge.priority = perm[n];
		edge.priority = perm[n];

		queue.push_back(edge);
	}

	/*将超边排序*/
	sort(queue.begin(), queue.end(), cmp);
}

//=================================================================================//


