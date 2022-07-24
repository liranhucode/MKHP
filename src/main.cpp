
#include "partitioner.h"

int main(int argc, char *argv[])
{
	std::cout << "*************************************************************************************" << std::endl;
	std::cout << "***MKHP v1.0(c) april 2018-, by huliran" << std::endl;
	std::cout << "***Build data: Friday april 13 09:12:19 2018" << std::endl;

	std::string file = argv[1];
	std::cout << file << std::endl;

	Hypergraph hgraph;
	hgraph.Parser(file);
	hgraph.Report();

	Partitioner partioner(&hgraph); 
	partioner.run();
	partioner.report();

	return 0;
}

