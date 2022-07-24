
#include "partitioner.h"

int main(int argc, char *argv[])
{
	std::cout << "*************************************************************************************" << std::endl;
	std::cout << "***MKHP v1.0(c) april 2022-, by huliran" << std::endl;

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

