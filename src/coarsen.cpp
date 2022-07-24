#include <unordered_set>
#include <iomanip>
#include <unordered_map>
#include "coarsen.h"
#include "mystruct.h"
#include "myutil.h"
#include "checkhgraph.h"


/********************************************************
/*��������� ��ͼ���ݽṹ, ��������
/*����ֵ��	�ֻ���ĳ�ͼ
/*���ܣ�		����ĳ�ͼ���С�ĳ�ͼ
*********************************************************/
s_hgraph * s_coarsen::coarsen_hgraph(s_hgraph *original_hgraph, int npart)
{
	//cout << endl;

	int level = 0;
	double maxvwgt = 0.0;
	s_hgraph *hgraph = original_hgraph;

	/*������ֳ�ͼ������*/
	coarsen_to = npart * 30;

	/*����������ƥ������Ȩֵ*/
	maxvwgt = ( 5 * hgraph->tvwgts ) / coarsen_to;

	do {

		/*��cmap�������ռ�*/
		if (hgraph->cmap.size() == 0)
		{  
			vector<int> cmap(hgraph->nvtxs, -1);
			hgraph->cmap = cmap;
		}

		/*�Գ���Ϊ��λ����ƥ��*/
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
/*��������� ��ͼ���ݽṹ ��������ȨֵԼ��
/*����ֵ��	��
/*���ܣ�		����ͼ�еĶ������ƥ��
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

	/*���ж���δƥ��, ���г���δ����*/
	vector<int> hyperedge_contracted(nhedges, NO_CONTRACTED);
	vector<int> matched(nvtxs, UNMATCHED);		

	/*���ɳ����������*/
	vector<hyperedge> queue;
	creat_rpq_for_hyperedge(hgraph, queue);				

 	for (n = 0; n < nhedges; ++n)
	{



		vertex_weight_sum = 0;
		all_vertices_no_match = YES;

		/*���α���������������ĳ���*/
		num = queue[n].number;
		
		if (hgraph->eptr[num + 1] - hgraph->eptr[num] > 20)
			continue;

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
		else if ( hgraph->is_dependent_hyperedge(num) ) 
		{		
			//cout << "�������ߣ� " << num << endl;
			/*��������, ���ж����Ϊ1�ĳ���*/
			hyperedge_contracted[num] = DISCARD;	

			for ( i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i ) {
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
			for ( i = hgraph->eptr[num]; i < hgraph->eptr[num + 1]; ++i ) {
				vertex = hgraph->eind[i];
				if (matched[vertex-1] == MATCHED) {						
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


	/*��δ����ƥ��ĳ��߽����������,��һ��Ѱ�ҿ��Խ���ƥ��Ķ���*/
	for (n = 0; n < nhedges; ++n)
	{
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
							count++;
						}
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}
			else
			{
				/*Ȩֵ������Լ��,���ó�����δƥ��Ķ���ƥ�䵽һ��*/
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i) {

					vertex = hgraph->eind[i];
					if (matched[vertex - 1] == UNMATCHED) {
														/*�޷���һ��ƥ��Ķ���ֱ��copy����һ��*/
						cnvtxs++;
						hgraph->cmap[vertex - 1] = cnvtxs;
						matched[vertex - 1] = MATCHED;
					}
				}
			}
			/*ͳ�Ƹó����еĿ��Խ���ƥ��Ķ������ */
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


	//*���cnvtxs�����Ƿ���ȷ*/
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
	//	cout << "�ֻ�������� " << count_of_vertices << " " << cnvtxs << endl;

	/*����cnvtxs��cmap�����hyperedge_contracted�������ɸ���һ����ͼ*/
	creat_coarser_hgraph(hgraph, cnvtxs, hyperedge_contracted);


}


/*��������,Ѱ��������ĳ����п���ֱ�ӽ���contract�ĳ���*/
void s_coarsen::random_match(s_hgraph * hgraph, double maxvwgt)
{

}




/********************************************************
/*��������� ��ͼ���ݽṹ ����һ����ͼ�Ķ����� ����������־����
/*����ֵ��	��
/*���ܣ�		���ɸ���һ����ͼ,��ַ������hgraph->coarser��
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



	unordered_map<string,int> hyperedge_table;		/*��ֵ�����ߵ�string��ʽ  ֵ������Ȩֵ*/
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


			/*ɾ���ظ��Ķ���*/
			sort(temp.begin(), temp.end());

			temp.erase(unique(temp.begin(), temp.end()), temp.end());

			if (temp.size() < 2){
				
				/*����hyperedge_contracted����*/
				hyperedge_contracted[n] = CONTRACTED;
			}
			else{

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

		if(hgraph->cmap[i - 1] == -1){
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



	chgraph->compute_hyperedge_degree();		/*���㳬�ߵĶ�*/
	chgraph->compute_total_weight();		/*���㶥��Ȩֵ֮��*/
	chgraph->compute_vertex_degree();		/*���㶥���*/
	chgraph->compute_incidence_hyperedge();	/*���㶥���������*/
	//chgraph->compute_adjcent_vertices();
}




/*����ƥ��ֻ��㷨*/
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

	//�����ٽ�����
	for (auto it = adj.begin(); it != adj.end(); ++it)
	{

		v = *it;
		weight = 0;

		if (unmatched[v - 1] == UNMATCHED)		//���ڽ�����δƥ��
		{
			vector<int> adj_incidence = hypergraph->incidence_hyperedge[v - 1];

			//����vertex���ڽ�����*it�������ߵ�degree֮��
			for (i = 0; i < incidence.size(); ++i){
				
				int n = incidence[i];
				vector<int>::iterator result = find(adj_incidence.begin(), adj_incidence.end(), n); //����n
				//�ҵ�
				if (result != adj_incidence.end())
					weight += (float)1 / hypergraph->hyperedge_degree[n];
			}

			if (weight == 0) {
				cout << "Error!" << endl;
				exit(0);
			}

			//Ѱ��conn���Ķ���
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



	unordered_map<string, int> hyperedge_table;		/*��ֵ�����ߵ�string��ʽ  ֵ������Ȩֵ*/
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


		/*ɾ���ظ��Ķ���*/
		sort(temp.begin(), temp.end());

		temp.erase(unique(temp.begin(), temp.end()), temp.end());

		if (temp.size() < 2) {
			continue;
		}
		else {

			/*������ת�����ַ����ͱ��ڴ洢, ��temp = { 0, 1, 2} �� string�� = "0,1,2" */
			hyperedge_string = vector_to_string(temp);

			if ((hyperedge_table_it = hyperedge_table.find(hyperedge_string)) != hyperedge_table.end()) {

				hyperedge_table_it->second += hypergraph->hewgts[n];		/*�ظ�����Ȩֵ+1*/
			}
			else {
				/*���û�ҵ�,��ԭ���ߵ�Ȩֵ������ȥ*/
				hyperedge_table.insert(make_pair(hyperedge_string, hypergraph->hewgts[n]));
			}

		}

	}


	map_convert_to_hyperedge(hyperedge_table, chgraph->eind, chgraph->eptr, chewgts, cnhedges);

	/*����cvwgts����*/
	for (i = 1; i <= nvtxs; i++) 
	{
		if (hypergraph->cmap[i - 1] == -1) {
			continue;
		}
		cvwgts[hypergraph->cmap[i - 1] - 1] += hypergraph->vwgts[i - 1];
	}



	/*���ôֻ��ĳ�ͼchgraph*/
	chgraph->vwgts = cvwgts;
	chgraph->hewgts = chewgts;
	chgraph->nvtxs = cnvtxs;
	chgraph->nhedges = cnhedges;
	hypergraph->coarser = chgraph;
	chgraph->finer = hypergraph;



	chgraph->compute_hyperedge_degree();		/*���㳬�ߵĶ�*/
	chgraph->compute_total_weight();		/*���㶥��Ȩֵ֮��*/
	chgraph->compute_vertex_degree();		/*���㶥���*/
	chgraph->compute_incidence_hyperedge();	/*���㶥���������*/

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

			if ( v != -1 ){		/*�ҵ�v*/
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

		edge.number = n;	/*�������*/

		edge.weight = hgraph->hewgts[n];		/*����Ȩֵ*/

												//edge.total_vertex_weigth = total_weight;			/*�����ж���Ȩֵ֮��*/

		edge.degree = hgraph->hyperedge_degree[n];		/*���߶�*/

		//edge.degree = k++;											//edge.priority = perm[n];
		edge.priority = perm[n];

		queue.push_back(edge);
	}

	/*����������*/
	sort(queue.begin(), queue.end(), cmp);
}

//=================================================================================//


