#include "hgraph.h"

#include <fstream>
#include <sstream> 
#include <vector>
#include <algorithm>
#include <iomanip>

s_hgraph::s_hgraph()
{
	tvwgts = 0;
	nvtxs = 0;
	nhedges = 0;
	nbnd = 0;
	//n_border = 0;
	mincut = 0;
	coarser = NULL;
	finer = NULL;
}

s_hgraph::~s_hgraph()
{


}


//=================================================================================//
void s_hgraph::compute_total_weight()
{

	int i = 1;
	int sum = 0;

	for (i = 1; i <= nvtxs; i++)
		sum += vwgts[i - 1];


	tvwgts = sum;
}
//=================================================================================//




//=================================================================================//

//=================================================================================//
void s_hgraph::compute_hyperedge_degree()
{
	int n = 0;
	int degree = 0;


	for (n = 0; n < nhedges; n++)
	{
		degree = eptr[n + 1] - eptr[n];
		hyperedge_degree.push_back(degree);
	}
}






void s_hgraph::compute_vertex_degree()
{

	//��ʼ��
	for (int i = 0; i < nvtxs; ++i){
		vertex_degree.push_back(0);
	}

	for (int n = 0; n < eind.size(); ++n){
		vertex_degree[eind[n] - 1]++;
	}

}


void s_hgraph::compute_incidence_hyperedge()
{
	int i = 0;
	int n = 0;
	int k = 0;

	//��ʼ��adjncy
	for ( i = 0; i < nvtxs; i++ ){
		vector<int> adjedge;
		incidence_hyperedge.push_back(adjedge);
	}
	
	/*��������*/
	for ( n = 0; n < nhedges; ++n ) {

		for (i = eptr[n]; i < eptr[n + 1]; ++i)
		{
			k = eind[i];		/*kΪ�����ж���*/
			incidence_hyperedge[k-1].push_back(n);
		}
	}




}
//=================================================================================//


/********************************************************
/*��������� ��ͼ���ݽṹ ����һ����ͼ�Ķ����� ����������־����
/*����ֵ��	��
/*���ܣ�		���ɸ���һ����ͼ,��ַ������hgraph->coarser��
*********************************************************/
void s_hgraph::compute_adjcent_vertices()
{
	int i = 0;
	int ii = 0;
	int iii = 0;
	int n = 0;
	int k = 0;


	for (i = 0; i < nvtxs; i++) {
		
		unordered_set<int> s;

		vector<int> hyperedge = incidence_hyperedge[i];

		for (ii = 0; ii < hyperedge.size(); ii++)
		{
			n = hyperedge[ii];

			for (iii = eptr[n]; iii < eptr[n + 1]; iii++)
			{
				k = eind[iii];
				if (k == i+1) {
					continue;
				}

				s.insert(k);
			}
		}

		adjcent_vertices.push_back(s);

	
	}
}

//=================================================================================//
void s_hgraph::print_hypergraph_info(string file, int nparts, int ubfactor)
{
	string name;
	name = file.substr(file.rfind("\\", file.length()) + 1, file.length());

	cout << "\nHyperGraph Information---------------------------------------------------------------" << endl;
	cout << "Name: " << name << "  #Vtxs: " << nvtxs << "  #Hedges: " << nhedges <<
		"  #nparts: " << nparts << "  UBfactor: " << setprecision(3) << (float)ubfactor/100 << endl;

}

/********************************************************
/*��������� ��Ա�������������
/*����ֵ��	��
/*���ܣ�		����where����,���㳬���Ƿ��и�,����ɱ�־λ����
*********************************************************/
void s_hgraph::compute_cutsign()
{
	int n = 0;
	int i = 0;
	int k = 0;
	vector<int> cutsign_temp(nhedges, 0);

	if (where.size() == 0)
	{
		cout << "No partition result!" << endl;
		exit(0);
	}

	for (n = 0; n < nhedges; ++n) {

		vector<int> tempwhere;

		for (i = eptr[n]; i < eptr[n + 1]; ++i) {

			k = eind[i];

			tempwhere.push_back(where[k - 1]);

		}

		sort(tempwhere.begin(), tempwhere.end());

		if (tempwhere.front() != tempwhere.back()) {
			cutsign_temp[n] = 1;					//������߱����ˣ���cutsign����Ϊ1 
		}

		cutsign = cutsign_temp;
	}



}
//=================================================================================//\



float s_hgraph::average_hyperedge_deg()
{

	float average;
	int total;

	average = 0.0;
	total = 0;


	for (int n = 0; n < nhedges; n++) {
		total += hyperedge_degree[n];
	}

	average = (float)total / nhedges;

	return average;
}



float s_hgraph::average_vertex_deg()
{

	float average;
	int total;

	average = 0.0;
	total = 0;


	for (int i = 0; i < nvtxs; i++) {
		total += vertex_degree[i];
	}

	average = (float)total / nvtxs;

	return average;

}


float s_hgraph::RatioOfHedgeAndVertex()
{
	float result;

	result = 0.0;


	result = (float) nhedges / nvtxs;

	return result;
}
//=================================================================================//


bool s_hgraph::is_dependent_hyperedge( int & num)
{
	int vertex;

	for (int i = eptr[num]; i < eptr[num + 1]; ++i) {

		vertex = eind[i];
		
		/*���һ��������������һ������Ķȴ���1,��ó��߷Ƕ�������*/
		if (vertex_degree[vertex - 1] > 1) {
			return false;
		}

	}

	return true;

}


bool s_hgraph::is_single_hyperedge(int & num)
{

	return eptr[num + 1] - eptr[num] < 2 ? true : false;

}