#include "rpart.h"
#include "boundaryfm.h"
#include "checkhgraph.h"
#include "multiobjective.h"
#include "vcycle.h"

s_rpart::s_rpart()
{
	hgraph = new s_hgraph;
	coarsen = new s_coarsen;
	init_2way_part = new s_initpart;
	refinement = new s_refine;

	cutsize = 0;
}


s_rpart::~s_rpart()
{
	delete hgraph;
}


/*�ݹ黮���㷨���㺯��,����Ϊ��ͼ,��������,���ز�ƽ������,��������*/
/*���Ϊfpart����,���滮�ֵ����ս��*/
void s_rpart::mlevel_recursive_partition(s_hgraph * hgraph, int npart, int ubfactor, vector<int> &part, int &objval)
{

	int i = 0;
	int bestcut = 0;
	int bestsoed = 0;
	int nruns = 0;
	vector <float> tpwgts;


	if (hgraph->nvtxs == 0 && hgraph->nhedges == 0) {
		printf("\t***Cannot bisect a graph with 0 vertices!\n"
			"\t***You are trying to partition a graph into too many parts!\n");
		exit(-1);
	}

	vector <int> bestwhere(hgraph->nvtxs, 0);
	allocate_tpwgts(tpwgts, npart);


	/*��ʼ��label�����part����*/
	for (i = 0; i < hgraph->nvtxs; ++i){
		part.push_back(0);
	}


	/*����10��ѡȡ�������ŵĽ��*/
	for ( nruns = 0; nruns < 1; nruns++)
	{
		 
		/*�ݹ���ַ�*/
		objval = mlevel_recursive_bisection(hgraph, npart, tpwgts, part, ubfactor, 0);

		hgraph->where = part;

		vector<int> pwgts(npart, 0);
		hgraph->pwgts.clear();

		for (i = 0; i < hgraph->nvtxs; i++) {
			pwgts[part[i]] += hgraph->vwgts[i];
		}

		/*objval�����ظ�������cut, ���Բ�������ʵ��mincut*/
		hgraph->mincut = compute_mincut(hgraph);			/*ͬʱ����cutsign����*/

		/*�����ⲿ�Ķ�֮��*/
		compute_soed(hgraph);

		//soed = obj::refine_soed(hgraph);


		if (nruns == 0 || hgraph->mincut < bestcut) {

			bestcut = hgraph->mincut;
			bestsoed = soed;

			for (i = 0; i < hgraph->nvtxs; ++i) {
				bestwhere[i] = hgraph->where[i];
			}

			if (bestcut == 0)
				break;
		}


		hgraph->bndlist.clear();
		hgraph->bndptr.clear();
		hgraph->fs.clear();
		hgraph->te.clear();
		hgraph->cutsign.clear();
		hgraph->where.clear();
		hgraph->pwgts.clear();

	}

	part = bestwhere;

	hgraph->mincut = bestcut;
	soed = bestsoed;

	vector<int> pwgts(npart, 0);
	hgraph->pwgts.clear();

	for (i = 0; i < hgraph->nvtxs; i++){
		pwgts[part[i]] += hgraph->vwgts[i];
	}

	for (i = 0; i < npart; i++){
		hgraph->pwgts.push_back(pwgts[i]);
	}


}



int s_rpart::mlevel_recursive_bisection(s_hgraph * hgraph, int npart, vector<float> tpwgts, vector<int> &part, int ubfactor, int fpart)
{
	int i;
	int objval;
	float wsum;

	vector<float> tpwgts2(2,0.0);

	s_hgraph *lhgraph = new s_hgraph();
	s_hgraph *rhgraph = new s_hgraph();


	vector<float>  left_tpwgts;
	vector<float>  right_tpwgts;

	/*ȷ�����ַ��ı���*/
	tpwgts2[0] = mysum((npart >> 1), tpwgts);
	tpwgts2[1] = 1 - tpwgts2[0];


	/*���ж༶���ַ�*/
	mlevel_bisection(hgraph, tpwgts2, ubfactor, npart);	  

	objval = hgraph->mincut;


	/*�����ֽ����Ӧ�����������*/
	if (hgraph->label.size() > 0) 
	{
		for (int i = 0; i < hgraph->nvtxs; i++) {
			part[hgraph->label[i] - 1] = hgraph->where[i] + fpart;
		}
	}
	else 
	{
		for (int i = 0; i < hgraph->nvtxs; i++) {
			part[i] = hgraph->where[i];
		}
	}

	/*����ͼ�ֳ���ͼ����ͼ*/
	if (npart > 2) {
		split_hypergraph(hgraph, lhgraph, rhgraph);			/*�ж��㱻����,��δ���*/
	}

	/*���·������Ȩֵ*/
	reallocate_tpwgts(tpwgts, npart, left_tpwgts, right_tpwgts);


	if (npart > 3) {

		objval += mlevel_recursive_bisection(lhgraph, npart >> 1, left_tpwgts, part, ubfactor, fpart);

		objval += mlevel_recursive_bisection(rhgraph, npart - (npart >> 1), right_tpwgts, part, ubfactor, fpart + (npart >> 1));

	}
	else if (npart == 3) {

		delete lhgraph;

		objval += mlevel_recursive_bisection(rhgraph, npart - (npart >> 1), right_tpwgts, part, ubfactor, fpart + (npart >> 1));

	}

	return objval;
}


/*�༶���ַ��Ķ��㺯��*/
void s_rpart::mlevel_bisection(s_hgraph * hgraph, vector<float> tpwgts2, int ubfactor, int npart)
{
	int testcut;

	int i = 0;

	s_hgraph *coarsest_hgraph = new s_hgraph();;

	/*�����ı���,Ԥ��*/
	if (tpwgts2.size() == 0) {
		cout << " Default tpwgts = {0.5, 0.5} " << endl;
		allocate_tpwgts(tpwgts2, 2);
	}



	for (i = 0; i < 1; i++)
	{

		/*�ֻ��׶�,��ø���һ����ͼ*/
		//coarsest_hgraph = coarsen->CoarsenHypergraph(hgraph, npart);
		coarsest_hgraph = coarsen->coarsen_hgraph(hgraph, npart);

		/*��ʼ����*/
		init_2way_part->random_bisection(coarsest_hgraph, ubfactor, tpwgts2);


		/*ϸ���׶�*/
		refinement->refine_2way(hgraph, coarsest_hgraph, tpwgts2);


		//cout << "befor v-cycle: " << hgraph->mincut << endl;

		/*��ѭ���Ż��㷨�Ż�������*/
		v_cycle::v_cycle(hgraph, 2, tpwgts2);

 
	}

	//cout << "bestcut = " << hgraph->mincut << endl;

	//check_cutsign(hgraph);
	//modified_cut_sign(hgraph);

	//testcut = compute_mincut(hgraph);

	/*����������֮��*/
	//span->compute_sum_of_span(hgraph);

	//delete coarsest_hgraph;

}


s_hgraph * s_rpart::setup_hgraph(int nvtxs, int nhedges, vector<int> eptr, vector<int> eind, vector<float> vwgts, vector<float> hewgts)
{
	hgraph = new s_hgraph;


	hgraph->nvtxs = nvtxs;

	hgraph->nhedges = nhedges;

	hgraph->eptr = eptr;

	hgraph->eind = eind;

	hgraph->vwgts = vwgts;

	hgraph->hewgts = hewgts;

	hgraph->compute_total_weight();		/*���㶥��Ȩֵ�ܺ�*/

	hgraph->compute_hyperedge_degree();	/*���㳬�߶�*/

	hgraph->compute_vertex_degree();	/*���㶥���*/

	hgraph->compute_incidence_hyperedge();		/*���㶥���������ĳ���*/

	hgraph->compute_adjcent_vertices();

	return hgraph;


}

void s_rpart::allocate_tpwgts(vector<float>& tpwgts, int npart)
{
	int i;
	float r = (float)1 / npart;

	for (i = 0; i < npart; ++i) {
		tpwgts.push_back(r);
	}


}

/*����where����Ľ��,����ͼ�ֳ���ͼ����ͼ��������*/
void s_rpart::split_hypergraph(s_hgraph * hgraph, s_hgraph * &lhgraph, s_hgraph * &rhgraph)
{
	int i, n, k;

	int lvtxs = 1;
	int rvtxs = 1;

	int lnhedges = 0;
	int rnhedges = 0;
	int left_ptr_base = 0;
	int right_ptr_base = 0;


	int left_cnt = 0;
	int right_cnt = 0;

	/*eptr�����׸�Ԫ��Ϊ 0*/
	lhgraph->eptr.push_back(0);
	rhgraph->eptr.push_back(0);


	vector<int> visit(hgraph->nvtxs, -1);
	vector<int> project(hgraph->nvtxs, -1);

	/*����ԭʼ��ͼ*/
	for (n = 0; n < hgraph->nhedges; n++)
	{
		//������0������1�����Ķ�����
		left_cnt = 0;		
		right_cnt = 0;

		for ( i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i )
		{	
			k = hgraph->eind[i];

			if (hgraph->where[k - 1] == 0) 
			{	// 0�����Ķ��㣬copy����ͼ��
				
				if (visit[k - 1] == -1) 
				{	// ����kδ������

					lhgraph->eind.push_back(lvtxs);
					lhgraph->vwgts.push_back(hgraph->vwgts[k - 1]);

					if (hgraph->label.size() == 0)
						lhgraph->label.push_back(k);
					else
						lhgraph->label.push_back(hgraph->label[k - 1]);

					visit[k - 1] = 1;
					project[k - 1] = lvtxs;
					lvtxs++;

				}
				else 
				{	// ����k�Ѿ�������	
					lhgraph->eind.push_back(project[k - 1]);
				}

				left_cnt++;
			}

			else 
			{
				// 1�����Ķ��㣬copy����ͼ
				if (visit[k - 1] == -1) 
				{

					rhgraph->eind.push_back(rvtxs);
					rhgraph->vwgts.push_back(hgraph->vwgts[k - 1]);

					if ( hgraph->label.size() == 0 )
						rhgraph->label.push_back(k);
					else
						rhgraph->label.push_back(hgraph->label[k - 1]);

					visit[k - 1] = 1;
					project[k - 1] = rvtxs;
					rvtxs++;

				}
				else 
				{
					rhgraph->eind.push_back(project[k - 1]);

				}
				right_cnt++;
			}

		}

		if (left_cnt > 0) 
		{
			left_ptr_base += left_cnt;
			lhgraph->eptr.push_back(left_ptr_base);

			if (hgraph->cutsign[n] == 1) 
				lhgraph->hewgts.push_back(hgraph->hewgts[n]);
			else
				lhgraph->hewgts.push_back(hgraph->hewgts[n]);

			lnhedges++;
		}

		if (right_cnt > 0) 
		{
			right_ptr_base += right_cnt;
			rhgraph->eptr.push_back(right_ptr_base);
			if (hgraph->cutsign[n] == 1)
				rhgraph->hewgts.push_back(hgraph->hewgts[n]);
			else
				rhgraph->hewgts.push_back(hgraph->hewgts[n]);

			rnhedges++;
		}
	}


	lhgraph->nvtxs = lvtxs - 1;		/*ʵ�ʶ�����*/
	rhgraph->nvtxs = rvtxs - 1;

	lhgraph->nhedges = lnhedges;
	rhgraph->nhedges = rnhedges;


	lhgraph->compute_total_weight();			/*���㶥��Ȩֵ�ܺ�*/
	lhgraph->compute_hyperedge_degree();		/*���㳬�߶�*/
	lhgraph->compute_vertex_degree();			/*���㶥���*/
	lhgraph->compute_incidence_hyperedge();		/*���㶥���������ĳ���*/
	lhgraph->compute_adjcent_vertices();

	rhgraph->compute_total_weight();			/*���㶥��Ȩֵ�ܺ�*/
	rhgraph->compute_hyperedge_degree();		/*���㳬�߶�*/
	rhgraph->compute_vertex_degree();			/*���㶥���*/
	rhgraph->compute_incidence_hyperedge();		/*���㶥���������ĳ���*/
	rhgraph->compute_adjcent_vertices();
	//nvtxs  nhedges
	lhgraph->nvtxs = lvtxs - 1;
	rhgraph->nvtxs = rvtxs - 1;
	lhgraph->nhedges = lnhedges;
	rhgraph->nhedges = rnhedges;

}

/*����where����Ľ��,����ͼ�ֳ���ͼ����ͼ��������*/
void s_rpart::split_hypergraph_c(s_hgraph * hgraph, s_hgraph * &lhgraph, s_hgraph * &rhgraph)
{
	int i, n, k;
	int lvtxs = 1;
	int rvtxs = 1;
	int lnhedges = 0;
	int rnhedges = 0;
	int left_ptr_base = 0;
	int right_ptr_base = 0;


	int left_cnt = 0;
	int right_cnt = 0;

	/*eptr�����׸�Ԫ��Ϊ 0*/
	lhgraph->eptr.push_back(0);
	rhgraph->eptr.push_back(0);


	vector<int> visit(hgraph->nvtxs, -1);
	vector<int> project(hgraph->nvtxs, -1);

	/*����ԭʼ��ͼ*/
	for (n = 0; n < hgraph->nhedges; n++)
	{
		left_cnt = 0;
		right_cnt = 0;


		if (hgraph->cutsign[n] == 0)
		{	/*δ����cut��, �ڲ�����, ֱ��copy����ͼ��*/


		}


		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i)
		{
			k = hgraph->eind[i];

			if (hgraph->where[k - 1] == 0) {
				//copy to left graph
				if (visit[k - 1] == -1) {
					lhgraph->eind.push_back(lvtxs);
					lhgraph->vwgts.push_back(hgraph->vwgts[k - 1]);
					left_cnt++;

					if (hgraph->label.size() == 0)
					{
						lhgraph->label.push_back(k);
					}
					else
					{
						lhgraph->label.push_back(hgraph->label[k - 1]);
					}

					visit[k - 1] = 1;
					project[k - 1] = lvtxs;
					lvtxs++;
				}
				else {
					lhgraph->eind.push_back(project[k - 1]);
					left_cnt++;
				}
			}

			else {
				//copy to right graph
				if (visit[k - 1] == -1) {
					rhgraph->eind.push_back(rvtxs);
					rhgraph->vwgts.push_back(hgraph->vwgts[k - 1]);
					right_cnt++;


					if (hgraph->label.size() == 0)
					{
						rhgraph->label.push_back(k);
					}
					else
					{
						rhgraph->label.push_back(hgraph->label[k - 1]);
					}

					visit[k - 1] = 1;
					project[k - 1] = rvtxs;
					rvtxs++;
				}
				else {
					rhgraph->eind.push_back(project[k - 1]);
					right_cnt++;
				}
			}

		}

		if (left_cnt > 0) {
			left_ptr_base += left_cnt;
			lhgraph->eptr.push_back(left_ptr_base);

			if (hgraph->cutsign[n] == 1)
			{
				lhgraph->hewgts.push_back(hgraph->hewgts[n]);
			}
			else
			{
				lhgraph->hewgts.push_back(hgraph->hewgts[n]);
			}

			lnhedges++;
		}

		if (right_cnt > 0) {

			right_ptr_base += right_cnt;
			rhgraph->eptr.push_back(right_ptr_base);

			if (hgraph->cutsign[n] == 1)
			{
				rhgraph->hewgts.push_back(hgraph->hewgts[n]);
			}
			else
			{
				rhgraph->hewgts.push_back(hgraph->hewgts[n]);
			}

			rnhedges++;
		}
	}


	lhgraph->nvtxs = lvtxs - 1;
	rhgraph->nvtxs = rvtxs - 1;

	lhgraph->nhedges = lnhedges;
	rhgraph->nhedges = rnhedges;


	lhgraph->compute_total_weight();			/*���㶥��Ȩֵ�ܺ�*/
	lhgraph->compute_hyperedge_degree();		/*���㳬�߶�*/
	lhgraph->compute_vertex_degree();			/*���㶥���*/
	lhgraph->compute_incidence_hyperedge();		/*���㶥���������ĳ���*/

	rhgraph->compute_total_weight();			/*���㶥��Ȩֵ�ܺ�*/
	rhgraph->compute_hyperedge_degree();		/*���㳬�߶�*/
	rhgraph->compute_vertex_degree();			/*���㶥���*/
	rhgraph->compute_incidence_hyperedge();		/*���㶥���������ĳ���*/

												//nvtxs  nhedges
	lhgraph->nvtxs = lvtxs - 1;
	rhgraph->nvtxs = rvtxs - 1;
	lhgraph->nhedges = lnhedges;
	rhgraph->nhedges = rnhedges;

}


/*���·��������Ȩֵ*/
void s_rpart::reallocate_tpwgts(vector<float>& tpwgts, int npart, vector<float>& left_tpwgts, vector<float>& right_tpwgts)
{

	float val;

	/*��ͼ��Ȩֵ����*/
	val = (float) 1 / (npart - (npart >> 1));

	for (int i = 0; i < npart - (npart >> 1); ++i) {

		right_tpwgts.push_back(val);

	}

	/*��npart <= 3, ����ͼ����Ҫ��������*/
	if (npart > 3){	

		/*��ͼ��Ȩֵ����*/
		val = (float) 1 / (npart >> 1);

		for (int i = 0; i < npart >> 1; ++i) {

			left_tpwgts.push_back(val);
		}
	}

}





/*�ݹ黮���㷨���㺯��,����Ϊ��ͼ,��������,���ز�ƽ������,��������*/
/*���Ϊfpart����,���滮�ֵ����ս��*/
void s_rpart::mlevel_recursive_partition_api(int nvtxs, int nhedges, vector<int>eptr, vector<int> eind, vector<float> vwgts, vector<float> hewgts,
	int npart, int ubfactor, vector<int> &part, int &edgecut)
{
	

	s_hgraph *hgraph;

	if (nvtxs == 0 && nhedges == 0) {
		printf("\t***Cannot bisect a graph with 0 vertices!\n"
			"\t***You are trying to partition a graph into too many parts!\n");
		exit(-1);
	}

	hgraph = setup_hgraph(nvtxs, nhedges, eptr, eind, vwgts, hewgts);

	mlevel_recursive_partition(hgraph, npart,  5, part, edgecut);

	cout << hgraph->mincut << endl;
	hgraph->where = part;

	//computing_cutsign(hgraph, part);

	hgraph->pwgts.clear();


	/**************************************/
	/*************printf info*************/
	/**************************************/
	vector<int> pwgts(npart, 0);
	for (int i = 0; i < nvtxs; i++)
	{
		pwgts[part[i]] += hgraph->vwgts[i];
	}

	for (int i = 0; i < npart; i++)
	{
		hgraph->pwgts.push_back(pwgts[i]);
	}

}



void s_rpart::compute_soed(s_hgraph * hgraph)
{
	int n = 0;
	int i = 0;
	int k = 0;

	soed = 0;

	for ( n = 0; n < hgraph->nhedges; ++n)
	{
		vector<int> v;


		for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i)
		{
			k = hgraph->eind[i];
			v.push_back(hgraph->where[k - 1]);
		}


		sort(v.begin(), v.end());
		v.erase(unique(v.begin(), v.end()), v.end());

		if (v.front() != v.back()) {

			soed += v.size();
		}

	}

}
