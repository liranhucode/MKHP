#pragma once

#ifndef MYUTIL_H
#define MYUTIL_H

#include <unordered_map>
#include <vector>
#include <string>
#include "mystruct.h"

#define size_s float
#define DEBUG 0

//#define random(a,b) (rand()%(b-a+1)+a)

using namespace std;


	int getInt(const char *buf);

	int getNumber(string &str, vector<string> &res);

	string vector_to_string(vector<int> v);

	bool cmp(const hyperedge &a, const hyperedge &b);


	void StringToVector(string s, vector<int> &v);


	void  map_convert_to_hyperedge(unordered_map<string, int> s, vector<int>& eind, vector<int>& eptr, vector<size_s> &hewgts, int &nhedges);

	void rand_permute(vector<int> &perm, int nums, int s);

	/* 产生一个随机数, {0-num} 之间 */
	int random_number(int min, int max);

	float mysum(int num, vector<float> tpwgts);

#endif // !MYUTIL_H
