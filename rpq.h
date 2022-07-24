#pragma once

#ifndef RPQ_H
#define RPQ_H


#include <vector>
#include <queue>
#include <map>
#include <set>

#include "mystruct.h"

using namespace std;


class s_rpq {


public:

	void rpq_insert(int n, float g);
	void rpq_update(int n, float g);
	void rpq_delete(int n);
	int rpq_length();
	int rpq_get_top();
	float rpq_top();

private:
	//map<int, float> mymap;
	set<gain_node> myset;
};


#endif