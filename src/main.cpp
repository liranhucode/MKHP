#include "partitioner.h"

int main(int argc, char *argv[])
{
	std::cout << "*************************************************************************************" << std::endl;
	std::cout << "***MKHP v1.0(c) april 2022-, by huliran" << std::endl;

	std::string file = argv[1];
	std::cout << file << std::endl;

	Timer timer("Parsing input hypergraph file");
	Hypergraph hgraph;
	hgraph.Parser(file);
	hgraph.Report();
	timer.Report("");

	timer.Restart("Hypergraph partitoner");
	Partitioner partioner(&hgraph); 
	partioner.run();
	partioner.report();
	timer.Report("");

	return 0;
}

