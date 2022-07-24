#pragma once
#pragma once
#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

using namespace std;
// ======================================================================= //
class s_params
{
public:
	s_params();
	~s_params();

	void parse_command(int argc, char *argv[]);


	string filename;		/*���볬ͼ�ļ�*/


	int  nparts;			/*��������*/


	string outfile;		


	int ubfactor;


	int nruns;


	int coarseTo;


	float reduct_ratio;

};

// ======================================================================= //


#endif