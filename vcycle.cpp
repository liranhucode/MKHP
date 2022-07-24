#include "vcycle.h"
#include "checkhgraph.h"
#include <iomanip>
#include "coarsen.h"
#include "refine.h"
#include "fm.h"

using namespace std;



/********************************************************
/*输入参数： 超图数据结构, 分区个数
/*返回值：	粗化后的超图
/*功能：		循环优化算法,在粗化阶段,只收缩那些只在一个分区的超边
/*			精化阶段采用FM算法
*********************************************************/
void v_cycle::v_cycle(s_hgraph * hgraph, int npart, vector<float> tpwgts2)
{

	s_refine *v_cycle_refine = new s_refine();
	s_hgraph *restricted_coarsest_hgraph = new s_hgraph();

	/*将当前图与下一级图的映射关系清空*/
	hgraph->cmap.clear();

#ifdef DEBUG
	check_cutsign(hgraph);
#endif

	/*循环优化V-cycle的粗化阶段,生成的超图中已经有分区*/
	restricted_coarsest_hgraph = restricted_coarsen_hgraph(hgraph, npart);

	/*计算出分区初始权值*/
	vector<int> pwgts(2, 0);

	for (int i = 0; i < restricted_coarsest_hgraph->nvtxs; ++i) 
	{
		pwgts[restricted_coarsest_hgraph->where[i]] += restricted_coarsest_hgraph->vwgts[i];
	}

	restricted_coarsest_hgraph->pwgts = pwgts;

	/*计算分区的cut*/
	restricted_coarsest_hgraph->mincut = compute_cut(restricted_coarsest_hgraph);

	/*循环优化V-cycle的精化阶段*/
	v_cycle_refine->refine_2way(hgraph, restricted_coarsest_hgraph, tpwgts2);
	//refine();

}

s_hgraph * v_cycle::restricted_coarsen_hgraph(s_hgraph * original_hgraph, int npart)
{

	//cout << endl;

	int level = 0;
	int i = 0;
	double maxvwgt = 0.0;
	s_hgraph *hgraph = original_hgraph;

	/*设置最粗超图顶点数*/
	int coarsen_to = 100;

	/*设置允许顶点匹配的最大权值*/
	maxvwgt = (2.5 * hgraph->tvwgts) / coarsen_to;

	do {

		/*给cmap数组分配空间*/
		if (hgraph->cmap.size() == 0)
		{
			vector<int> cmap(hgraph->nvtxs, -1);
			hgraph->cmap = cmap;
		}

		restricted_hyperedge_match(hgraph, maxvwgt);

		level++;

		hgraph = hgraph->coarser;

		check_cutsign(hgraph);

#if DEBUG
		check_hgraph(hgraph);
#endif

	} while (hgraph->nvtxs > coarsen_to
		&&hgraph->nvtxs < 0.85 * hgraph->finer->nvtxs
		/*&&hgraph->nhedges < 0.85 * hgraph->finer->nhedges*/
		/*&&hgraph->nhedges < hgraph->nvtxs / 2*/
		);



#if DEBUG
	s_hgraph *tmp = original_hgraph;

	for (int i = 0; i <= level; i++) {
		cout << " [ " << setw(5) << tmp->nvtxs << "  " << setw(5) << tmp->nhedges << "  " << setw(5) << setiosflags(ios::fixed) << setprecision(2) << tmp->average_vertex_deg() << \
			"  " << setw(5) << tmp->average_hyperedge_deg() << "  " << setw(5) << (float)tmp->average_vertex_deg() / tmp->average_hyperedge_deg() << " ] " << endl;
		tmp = tmp->coarser;
	}

#endif

	return hgraph;

}

void v_cycle::restricted_hyperedge_match(s_hgraph * hgraph, double maxvwgt)
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
	vector<int> coarser_where;

	/*生成超边排序队列*/
	vector<hyperedge> queue;
	creat_rpq_for_hyperedge(hgraph, queue);

	for (n = 0; n < nhedges; ++n)
	{	/*依次遍历队列中已排序的超边*/

		vertex_weight_sum = 0;
		all_vertices_no_match = YES;
		num = queue[n].number;


		if (hgraph->cutsign[num] == 1)
		{	/*已经被切割的超边,将直接copy到下一级超图*/
			if (hgraph->is_dependent_hyperedge(num))
			{	/*独立超边被切割*/
				hyperedge_contracted[num] = DISCARD;

				for ( i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i ) 
				{	/*将独立超边中的顶点在cmap数组中标记成-1,表示该顶点不会copy到下一级超图*/
					vertex = hgraph->eind[i];
					hgraph->cmap[vertex - 1] = -1;
					matched[vertex - 1] = DISCARD;
				}
			}
			else 
			{	/*非独立超边被切割,将直接copy到下一级超图*/
				for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i) 
				{
					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED)
					{
						cnvtxs++;
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
						coarser_where.push_back(hgraph->where[vertex - 1]);
					}

				}
			}
			
		}
		else {
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
			else if (hgraph->is_dependent_hyperedge(num))
			{
				//cout << "独立超边： " << num << endl;
				/*独立超边, 所有顶点度为1的超边*/
				hyperedge_contracted[num] = DISCARD;

				for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i) {
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
				for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i) {
					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == MATCHED) {
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
						coarser_where.push_back(hgraph->where[hgraph->eind[hgraph->eptr[num]] - 1]);

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
		

	}



	/*对未进行匹配的超边进行重新审查,进一步寻找可以进行匹配的顶点*/
	for (n = 0; n < nhedges; ++n)
	{

		/*已经被切割的超边,将直接copy到下一级超图*/
		if (hgraph->cutsign[n] == 1)
		{
			continue;
		}

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
							coarser_where.push_back(hgraph->where[vertex-1]);
							count++;
						}
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}
			else
			{
				/*权值不符合约束*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {
						/*无法进一步匹配的顶点直接copy到下一级*/
						cnvtxs++;
						coarser_where.push_back(hgraph->where[vertex - 1]);
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}

		}


	}




	/*根据cnvtxs、cmap数组和hyperedge_contracted数组生成更粗一级超图*/
	creat_coarser_hgraph(hgraph, coarser_where,  cnvtxs, hyperedge_contracted );



}

int v_cycle::compute_cut(s_hgraph * hgraph)
{


	int  i, n, k;
	int cut = 0;
	int nhedges = hgraph->nhedges;


	for (n = 0; n < nhedges; ++n)
	{
		vector<int> temp;

		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {
			k = hgraph->eind[i];
			temp.push_back(hgraph->where[k - 1]);
		}


		std::sort(temp.begin(), temp.end());

		if (temp.front() != temp.back()) {

			cut += hgraph->hewgts[n];

			if (hgraph->cutsign[n] == 0)
			{
				hgraph->cutsign[n] = 1;
			}

		}
		else
		{
			if (hgraph->cutsign[n] == 1)
				hgraph->cutsign[n] = 0;
		}
	}

	return cut;
}

void v_cycle::creat_rpq_for_hyperedge(s_hgraph * hgraph, vector<hyperedge>& queue)
{
	int n, k;
	int nhedges = hgraph->nhedges;
	int total_weight;
	vector<int> perm;

	//rand_permute(perm, nhedges, 0);

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

														//edge.priority = perm[n];

		queue.push_back(edge);
	}

	/*将超边排序*/
	sort(queue.begin(), queue.end(), cmp);

}


void v_cycle::creat_coarser_hgraph(s_hgraph * hgraph,  vector<int>coarser_where, int cnvtxs, vector<int> &hyperedge_contracted)
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



	unordered_map<string, int> hyperedge_table;		/*键值：超边的string格式  值：超边权值*/
	unordered_map<string, int>::iterator hyperedge_table_it;
	string hyperedge_string;



	for (n = 0; n < nhedges; n++)
	{
		vector<int> temp;
		int cnt = 0;


		if (hyperedge_contracted[n] == NO_CONTRACTED) {

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; i++) {

				vertex = hgraph->eind[i];

				if (hgraph->cmap[vertex - 1] != -1)
					temp.push_back(hgraph->cmap[vertex - 1]);

			}


			/*删除重复的顶点*/
			sort(temp.begin(), temp.end());

			temp.erase(unique(temp.begin(), temp.end()), temp.end());

			if (temp.size() < 2) {

				/*修正hyperedge_contracted数组*/
				hyperedge_contracted[n] = CONTRACTED;
			}
			else {

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

		if (hgraph->cmap[i - 1] == -1) {
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

	chgraph->where = coarser_where;

	chgraph->compute_cutsign();

	chgraph->compute_hyperedge_degree();		/*计算超边的度*/
	chgraph->compute_total_weight();		/*计算顶点权值之和*/
	chgraph->compute_vertex_degree();		/*计算顶点度*/
	chgraph->compute_incidence_hyperedge();	/*计算顶点关联超边*/



}
