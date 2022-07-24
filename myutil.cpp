

#include <fstream>
#include <sstream> 
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <time.h>
#include <sys/timeb.h>
#include "myutil.h"
#include "mystruct.h"
using namespace std;




//=================================================================================//
//count number of integers of a line
//=================================================================================//
int getInt(const char * buf)
{
	const char *loc = buf;
	int res = 0;
	atoi(buf);
	loc = strstr(buf, " ");
	while (loc != NULL)
	{
		atoi(loc + 1);
		res++;
		loc = strstr(loc + 1, " ");
	}
	return res;
}


int getNumber(string & str, vector<string> &res)
{

	string result;
	stringstream input(str);

	while (input >> result)
		res.push_back(result);

	return res.size();
}


//int数组转换成，string型，每个整型中间用1隔开
string vector_to_string(vector<int> v)
{
	string s;

	s += to_string(v[0]);

	for (int i = 1; i < v.size(); i++) {

		s += ' ' + to_string(v[i]);

	}

	return s;
}

//将string类型转换成vector类型
void StringToVector(string s, vector<int>& v)
{
	int temp;
	char flag = ' ';

	int position = 0;
	int i = 0;

	while ((position = s.find(flag, position)) != string::npos)
	{
		temp = atoi(s.substr(i, position).c_str());
		v.push_back(temp);
		i = position++;
	}

	temp = atoi(s.substr(i, s.size()).c_str());
	v.push_back(temp);


}



//将unordered_set中的vector<int> 变量转换到eind和eptr数组中
void map_convert_to_hyperedge(unordered_map<string,int> s, vector<int>& eind, vector<int>& eptr , vector<size_s> &hewgts, int &nhedges)
{
	int base = 0;

	eptr.push_back(base);

	for (unordered_map<string, int>::iterator it = s.begin(); it != s.end(); it++) {

		vector<int> v;

		StringToVector(it->first, v);
		eind.insert(eind.end(), v.begin(), v.end());

		base += v.size();
		eptr.push_back(base);
		hewgts.push_back(it->second);
	}

	nhedges = s.size();
}


bool cmp(const hyperedge &a,  const hyperedge &b)
{


	//x为降序
	//y为升序
	//number为升序

	if (a.weight != b.weight) {		/*weight为升序*/

		return a.weight > b.weight;

	}

	else if (a.degree != b.degree) {		/*超边的degree为降序*/

		return a.degree < b.degree;

	}

	//else if (a.total_vertex_weigth != b.total_vertex_weigth){		/*超边中顶点权值之和为降序*/

	//	return a.total_vertex_weigth < b.total_vertex_weigth;

	//}

	else{			/*超边的序号为升序*/

					//return a.priority > b.priority;
		return a.priority >  b.priority;

	}



}


//为perm数组生成一个从（0 - nvtxs-1）的随机序列
void rand_permute(vector<int> &perm, int nums, int s)
{
	struct timeb time_seed;


	//srand(s);		//test

	ftime(&time_seed);
	srand(time_seed.time * 1000 + time_seed.millitm);


	//srand(1);

	for (int i = 1; i <= nums; ++i){
		perm.push_back(i);
	}

	random_shuffle(perm.begin(), perm.end());


}


int random_number(int min, int max)
{
	srand(time(0));

	int n = rand() / (RAND_MAX + 1) * (max - min) + min;

	vector<int> perm;

	for (int i = min; i < max; ++i){
		perm.push_back(i);
	}

	random_shuffle(perm.begin(), perm.end());

	return perm[0];

}


float mysum(int num, vector<float> tpwgts)
{
	float sum = 0;
	int i;

	if (num > tpwgts.size()){
		printf(" Tpwgts Error!");
		exit(0);
	}

	for ( i = 0; i < num; i++) {
		sum += tpwgts[i];
	}

	return sum;
}