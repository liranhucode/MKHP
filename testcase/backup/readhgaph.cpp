#include "hgraph.h"
#include <fstream>
#include <sstream> 
#include <vector>
#include <algorithm>
#include "myutil.h"
#include <time.h>

//=================================================================================//
//Reading hypergraph from  hypergraphFile
//=================================================================================//
void s_hgraph::read_hypergraph(string hypergraph_file)
{

	int nfields = 0;
	int i = 0;
	int length = 0;
	int nums = 0;
	int base = 0;
	int value = 0;

	string spec;
	string line;


	const char *d = " ";
	const char *p;


	std::ifstream _fin(hypergraph_file);
	if (!_fin) {
		std::cout << "Error: can not open this file!" << std::endl;
		exit(0);
	}

	getline(_fin, spec);

	p = spec.data();

	//Read: nvtxs nhedges fmt
	nfields = sscanf_s(p, "%d" "%d" "%d", &nhedges, &nvtxs, &fmt);

	if (nfields < 2) {
		std::cout << "The input file does not specify the number of vertices and edges." << std::endl;
		exit(1);
	}

	if (nvtxs <= 0 || nhedges <= 0) {
		std::cout << "The supplied nvtxs: " << nvtxs << " and" << nhedges << " must be postive." << std::endl;
		exit(1);
	}

	if (fmt > 111) {
		std::cout << "Cannot read this type of fime format: fmt = " << fmt << std::endl;
		exit(1);
	}


	eptr.push_back(0);

	if (fmt == 0)
	{
		//reading the array eptr and eind

		while (getline(_fin, line))
		{
			vector<string> res;

			nums = getNumber(line, res);

			if (nums < 2){
				cout << "Line ERROR! " << endl;
			}

			eptr.push_back(base + nums);
			base += nums;


			for (i = 0; i < (signed)res.size(); i++)
			{
				value = atoi(res[i].c_str());
				//cout << value;
				eind.push_back(value);
			}

		}


		/*初始默认顶点和超边的权值均为1*/
		vector<float> tmp1(nvtxs, 1);
		vector<float> tmp2(nhedges, 1);
		vwgts = tmp1;
		hewgts = tmp2;




		compute_total_weight();		/*计算顶点权值总和*/
		compute_hyperedge_degree();	/*计算超边度*/
		compute_vertex_degree();	/*就算顶点度*/

		clock_t start = clock();
		compute_incidence_hyperedge();		/*计算顶点所关联的超边*/
		compute_adjcent_vertices();
		clock_t end = clock();
		cout <<  (double)(end - start) / ((clock_t)1000) << "sec" << endl;
	}


}
//=================================================================================//