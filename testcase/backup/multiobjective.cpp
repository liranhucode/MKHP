#include "multiobjective.h"

#include <algorithm>
#include <vector>


using namespace std;



/********************************************************
/*��������� ��ͼ���ݽṹ
/*����ֵ��	�Ż����soed
/*���ܣ�		�Ż��ⲿ��֮��
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


			/*�ⲿ�ȴ���2�ĳ���Ӧ�ó����Ż�*/
			if (v.size() > 2)
			{
				for (i = hgraph->eptr[n]; i < hgraph->eptr[n + 1]; ++i)
				{
					k = hgraph->eind[i];
					cout << k << "(" << hgraph->vertex_degree[k - 1] << "�� ";
				}

				cout << endl;
			}

		}

	}


	return 0;
}
