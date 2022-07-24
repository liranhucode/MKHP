#include "vcycle.h"
#include "checkhgraph.h"
#include <iomanip>
#include "coarsen.h"
#include "refine.h"
#include "fm.h"

using namespace std;



/********************************************************
/*��������� ��ͼ���ݽṹ, ��������
/*����ֵ��	�ֻ���ĳ�ͼ
/*���ܣ�		ѭ���Ż��㷨,�ڴֻ��׶�,ֻ������Щֻ��һ�������ĳ���
/*			�����׶β���FM�㷨
*********************************************************/
void v_cycle::v_cycle(s_hgraph * hgraph, int npart, vector<float> tpwgts2)
{

	s_refine *v_cycle_refine = new s_refine();
	s_hgraph *restricted_coarsest_hgraph = new s_hgraph();

	/*����ǰͼ����һ��ͼ��ӳ���ϵ���*/
	hgraph->cmap.clear();

#ifdef DEBUG
	check_cutsign(hgraph);
#endif

	/*ѭ���Ż�V-cycle�Ĵֻ��׶�,���ɵĳ�ͼ���Ѿ��з���*/
	restricted_coarsest_hgraph = restricted_coarsen_hgraph(hgraph, npart);

	/*�����������ʼȨֵ*/
	vector<int> pwgts(2, 0);

	for (int i = 0; i < restricted_coarsest_hgraph->nvtxs; ++i) 
	{
		pwgts[restricted_coarsest_hgraph->where[i]] += restricted_coarsest_hgraph->vwgts[i];
	}

	restricted_coarsest_hgraph->pwgts = pwgts;

	/*���������cut*/
	restricted_coarsest_hgraph->mincut = compute_cut(restricted_coarsest_hgraph);

	/*ѭ���Ż�V-cycle�ľ����׶�*/
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

	/*������ֳ�ͼ������*/
	int coarsen_to = 100;

	/*����������ƥ������Ȩֵ*/
	maxvwgt = (2.5 * hgraph->tvwgts) / coarsen_to;

	do {

		/*��cmap�������ռ�*/
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

	/*���ж���δƥ��, ���г���δ����*/
	vector<int> hyperedge_contracted(nhedges, NO_CONTRACTED);
	vector<int> matched(nvtxs, UNMATCHED);
	vector<int> coarser_where;

	/*���ɳ����������*/
	vector<hyperedge> queue;
	creat_rpq_for_hyperedge(hgraph, queue);

	for (n = 0; n < nhedges; ++n)
	{	/*���α���������������ĳ���*/

		vertex_weight_sum = 0;
		all_vertices_no_match = YES;
		num = queue[n].number;


		if (hgraph->cutsign[num] == 1)
		{	/*�Ѿ����и�ĳ���,��ֱ��copy����һ����ͼ*/
			if (hgraph->is_dependent_hyperedge(num))
			{	/*�������߱��и�*/
				hyperedge_contracted[num] = DISCARD;

				for ( i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i ) 
				{	/*�����������еĶ�����cmap�����б�ǳ�-1,��ʾ�ö��㲻��copy����һ����ͼ*/
					vertex = hgraph->eind[i];
					hgraph->cmap[vertex - 1] = -1;
					matched[vertex - 1] = DISCARD;
				}
			}
			else 
			{	/*�Ƕ������߱��и�,��ֱ��copy����һ����ͼ*/
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
			/*�жϳ��ߵ�����*/
			if (hgraph->is_single_hyperedge(num))
			{
				/*��������ĳ���*/
				hyperedge_contracted[num] = DISCARD;
				vertex = hgraph->eind[hgraph->eptr[num]];
				//cout << "single edge: " << vertex << endl;

				/*����ö����Ѿ���ƥ���,������*/
				if (matched[vertex - 1] == UNMATCHED) {
					/*�����������еĶ�����cmap�����б�ǳ�-1,��ʾ�ö��㲻��copy����һ����ͼ*/
					hgraph->cmap[vertex - 1] = -1;
					matched[vertex - 1] = DISCARD;
				}

			}
			else if (hgraph->is_dependent_hyperedge(num))
			{
				//cout << "�������ߣ� " << num << endl;
				/*��������, ���ж����Ϊ1�ĳ���*/
				hyperedge_contracted[num] = DISCARD;

				for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i) {
					/*�����������еĶ�����cmap�����б�ǳ�-1,��ʾ�ö��㲻��copy����һ����ͼ*/
					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {

						hgraph->cmap[vertex - 1] = -1;
						matched[vertex - 1] = DISCARD;
					}
				}
			}
			else {

				/*�жϸó����Ƿ�δƥ���*/
				for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i) {
					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == MATCHED) {
						all_vertices_no_match = NO;
						break;
					}
				}

				if (all_vertices_no_match == YES) {		/*�ó��������ж����δƥ��*/

					for (i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; i++) {
						vertex = hgraph->eind[i];
						vertex_weight_sum += hgraph->vwgts[vertex - 1];
					}


					if (vertex_weight_sum <= maxvwgt) {
						/*���㶥��ȨֵԼ��,�ó��߿��Խ���ƥ��*/
						cnvtxs++;
						coarser_where.push_back(hgraph->where[hgraph->eind[hgraph->eptr[num]] - 1]);

						hyperedge_contracted[num] = CONTRACTED;

						for (ii = hgraph->eptr[num]; ii < hgraph->eptr[num + 1]; ii++) {

							vertex = hgraph->eind[ii];
							hgraph->cmap[vertex - 1] = cnvtxs;		/*������Ŵ�1��ʼ*/
							matched[vertex - 1] = MATCHED;
						}
					}
				}
			}
		}
		

	}



	/*��δ����ƥ��ĳ��߽����������,��һ��Ѱ�ҿ��Խ���ƥ��Ķ���*/
	for (n = 0; n < nhedges; ++n)
	{

		/*�Ѿ����и�ĳ���,��ֱ��copy����һ����ͼ*/
		if (hgraph->cutsign[n] == 1)
		{
			continue;
		}

		if (hyperedge_contracted[n] == NO_CONTRACTED)
		{
			count = 0;
			vertex_weight_sum = 0;

			/*����δƥ��Ķ���Ȩֵ֮��*/
			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

				vertex = hgraph->eind[i];
				if (matched[vertex - 1] == UNMATCHED) {
					vertex_weight_sum += hgraph->vwgts[vertex - 1];
				}
			}

			if (vertex_weight_sum <= maxvwgt)
			{
				/*Ȩֵ����Լ��,���ó�����δƥ��Ķ���ƥ�䵽һ��*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {

						if (count == 0)
						{
							cnvtxs++;			/*����ֻ����һ��*/
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
				/*Ȩֵ������Լ��*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {
						/*�޷���һ��ƥ��Ķ���ֱ��copy����һ��*/
						cnvtxs++;
						coarser_where.push_back(hgraph->where[vertex - 1]);
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}

		}


	}




	/*����cnvtxs��cmap�����hyperedge_contracted�������ɸ���һ����ͼ*/
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

		edge.number = n;	/*�������*/

		edge.weight = hgraph->hewgts[n];		/*����Ȩֵ*/

												//edge.total_vertex_weigth = total_weight;			/*�����ж���Ȩֵ֮��*/

		edge.degree = hgraph->hyperedge_degree[n];		/*���߶�*/

														//edge.priority = perm[n];

		queue.push_back(edge);
	}

	/*����������*/
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



	unordered_map<string, int> hyperedge_table;		/*��ֵ�����ߵ�string��ʽ  ֵ������Ȩֵ*/
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


			/*ɾ���ظ��Ķ���*/
			sort(temp.begin(), temp.end());

			temp.erase(unique(temp.begin(), temp.end()), temp.end());

			if (temp.size() < 2) {

				/*����hyperedge_contracted����*/
				hyperedge_contracted[n] = CONTRACTED;
			}
			else {

				/*������ת�����ַ����ͱ��ڴ洢, ��temp = { 0, 1, 2} �� string�� = "0,1,2" */
				hyperedge_string = vector_to_string(temp);

				if ((hyperedge_table_it = hyperedge_table.find(hyperedge_string)) != hyperedge_table.end()) {

					hyperedge_table_it->second += hgraph->hewgts[n];		/*�ظ�����Ȩֵ+1*/
				}
				else {
					/*���û�ҵ�,��ԭ���ߵ�Ȩֵ������ȥ*/
					hyperedge_table.insert(make_pair(hyperedge_string, hgraph->hewgts[n]));
				}

			}

		}
	}


	map_convert_to_hyperedge(hyperedge_table, chgraph->eind, chgraph->eptr, chewgts, cnhedges);

	/*����cvwgts����*/
	for (i = 1; i <= nvtxs; i++) {

		if (hgraph->cmap[i - 1] == -1) {
			continue;
		}
		cvwgts[hgraph->cmap[i - 1] - 1] += hgraph->vwgts[i - 1];
	}

	/*���cvwgts����*/
	//for (i = 0; i < cnvtxs; ++i)
	//{
	//	if (cvwgts[i] == 0)
	//		cout << "ȨֵС��0" << endl;
	//}


	/*���ôֻ��ĳ�ͼchgraph*/
	chgraph->vwgts = cvwgts;
	chgraph->hewgts = chewgts;
	chgraph->nvtxs = cnvtxs;
	chgraph->nhedges = cnhedges;
	hgraph->coarser = chgraph;
	chgraph->finer = hgraph;

	chgraph->where = coarser_where;

	chgraph->compute_cutsign();

	chgraph->compute_hyperedge_degree();		/*���㳬�ߵĶ�*/
	chgraph->compute_total_weight();		/*���㶥��Ȩֵ֮��*/
	chgraph->compute_vertex_degree();		/*���㶥���*/
	chgraph->compute_incidence_hyperedge();	/*���㶥���������*/



}
