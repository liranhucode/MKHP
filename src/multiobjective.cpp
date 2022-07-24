#include "multiobjective.h"

#include <algorithm>
#include <vector>


using namespace std;



/********************************************************
/*输入参数： 超图数据结构
/*返回值：	优化后的soed
/*功能：		优化外部度之和
*********************************************************/
int obj::refine_soed(s_hgraph * hgraph)
{
	 
	int i = 0;
	int n = 0;
	int k = 0;

	for (n = 0; n < hgraph->nhedges; ++n)
	{
		if (hgraph->cutsign[n] == 1)
		{
			vector<int> v;

			for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i)
			{
				k = hgraph->eind[i];
				v.push_back(hgraph->where[k - 1]);
			}


			sort(v.begin(), v.end());
			v.erase(unique(v.begin(), v.end()), v.end());


			/*外部度大于2的超边应该尝试优化*/
			if (v.size() > 2)
			{
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i)
				{
					k = hgraph->eind[i];
					cout << k << "(" << hgraph->vertex_degree[k - 1] << "） ";
				}

				cout << endl;
			}

		}

	}


	return 0;
}
