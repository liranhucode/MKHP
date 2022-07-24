#include "rpq.h"
#include <algorithm>

/*向队列中插入元素*/
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

/*更新队列中的元素，找到it->number = n的位置，然后删除，最后重新插入*/
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



/*删除队列中 number = k的元素*/
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

/*从数组中取出最大值，一般数组已经排好序了,最大元素存放在下边为0的位置*/
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
