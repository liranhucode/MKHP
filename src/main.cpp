#include "partitioner.h"

int main(int argc, char *argv[])
{
	std::cout << "*************************************************************************************" << std::endl;
	std::cout << "***MKHP v1.0(c) april 2022-, by huliran" << std::endl;
	std::cout << "*************************************************************************************" << std::endl;

	std::string file = argv[1];
	std::cout << file << std::endl;

	Timer timer("Parsing input hypergraph file");
	//Hypergraph hgraph;
	std::shared_ptr<Hypergraph> hgraph = std::make_shared<Hypergraph>(); 
	hgraph->Parser(file);
	hgraph->Report();
	timer.Report("");
	timer.Restart("Hypergraph partitoner");
	Partitioner partioner(std::move(hgraph)); 
	//partioner.run();
	//partioner.report();
	timer.Report("");
	std::cout << "*************************************************************************************" << std::endl;

	return 0;
}

