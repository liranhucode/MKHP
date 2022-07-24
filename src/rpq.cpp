#include "rpq.h"
#include <algorithm>

/*������в���Ԫ��*/
void s_rpq::rpq_insert(int n, float g)
{
	gain_node newnode;
	newnode.number = n;
	newnode.gain = g;

	//for (set<gain_node>::iterator it = myset.begin(); it != myset.end(); it++) {

	//	if (it->number == n)
	//	{
	//		return;
	//	}
	//}

	myset.insert(newnode);
}

/*���¶����е�Ԫ�أ��ҵ�it->number = n��λ�ã�Ȼ��ɾ����������²���*/
void s_rpq::rpq_update(int n, float g)
{
	gain_node newnode{ n,g };

	for (set<gain_node>::iterator it = myset.begin(); it != myset.end(); it++) {

		if (it->number == n)
		{
			myset.erase(it);
			myset.insert(newnode);
			break;
		}
	}

}



/*ɾ�������� number = k��Ԫ��*/
void s_rpq::rpq_delete(int k)
{
	set<gain_node>::iterator itor2;
	for (set<gain_node>::iterator it = myset.begin(); it != myset.end(); it++) {

		if (it->number == k) {
			itor2 = it;
			myset.erase(itor2);
			break;
		}

	}

}

int s_rpq::rpq_length()
{
	return myset.size();
}

/*��������ȡ�����ֵ��һ�������Ѿ��ź�����,���Ԫ�ش�����±�Ϊ0��λ��*/
int s_rpq::rpq_get_top()
{

	if (myset.size() == 0)
		return -1;

	set<gain_node>::iterator it = myset.end();
	it--;

	int res = it->number;

	myset.erase(it);


	return res;

}

float s_rpq::rpq_top()
{
	if (myset.size() == 0)
		return -1;

	set<gain_node>::iterator it = myset.end();
	it--;



	return 	it->gain;
}
