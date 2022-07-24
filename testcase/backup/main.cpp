
#include <iostream>
#include <stdio.h>
#include <vector>
#include <time.h>

#include "parameter.h"
#include "hgraph.h"
#include "rpart.h"
#include "myutil.h"

using namespace std;

s_params *para = new s_params();
s_hgraph *input_hgraph = new s_hgraph();
s_rpart *bisection = new s_rpart();

s_rpart *rpart = new s_rpart();


int main(int argc, char *argv[])
{
	cout << "*************************************************************************************" << endl;
	cout << "***MKHP v1.0(c) april 2018-, by huliran" << endl;
	cout << "***Build data: Friday april 13 09:12:19 2018" << endl;


/*
	int value = 0;
	vector<int> part;
	int hyperedgecut;
	int cutsize;

	/*���������в���*/
	para->parse_command(argc, argv);

	//��ȡ��ͼ
	input_hgraph->read_hypergraph(para->filename);


	/*��ӡ��ͼ��Ϣ*/
	input_hgraph->print_hypergraph_info(para->filename, para->nparts, para->ubfactor);
	
	//check_graph(hgraph);

	clock_t start = clock();

	vector<float> tpwgts;

	//recursive_partitioning
	//bisection->mlevel_bisection(input_hgraph, tpwgts, para->ubfactor);

	/*�ݹ���ַ�����*/


	rpart->mlevel_recursive_partition_api(input_hgraph->nvtxs, input_hgraph->nhedges, input_hgraph->eptr, input_hgraph->eind, input_hgraph->vwgts, input_hgraph->hewgts,
		para->nparts, para->ubfactor, part, cutsize);

	clock_t end = clock();

	//cout << "\n---------------------------------------------------------------------" << endl;
	//cout << "Summary for the " << para->nparts << "-way partition:" << endl;
	//cout << "\tHyperedge Cut:              " << input_hgraph->mincut << endl;
	//cout << "\tSum of external degrees:    " << rpart->soed << endl;

	//cout << "\tpartition size:" << endl;
	//for (int i = 0; i < para->nparts; i++) {
	//	cout << "\t\tpart " << i << ": " << input_hgraph->pwgts[i] << endl;
	//}
	//cout << "\n----------------------------------------------------------------------" << endl;


	cout << "\nTiming information-----------------------------------------------------" << endl;
	cout << "  Partitioning Time\t" << (double)(end - start) / CLOCKS_PER_SEC << "sec" << endl;

	cout << "MKHP is Done!" << endl;


	system("pause");
	*/
	return 0;
}

