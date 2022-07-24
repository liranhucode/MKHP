
#include "hypergraph.h"
#include <fstream>
#include <sstream>


void Hypergraph::Parser(std::string &file_name)
{
    std::ifstream fin;
    fin.open(file_name.c_str(), std::ios::in);
    if (!fin.is_open())
    {
        std::cerr << "cannot open the file " << file_name << std::endl;
        exit(1);
    }
    std::string line;

	std::getline(fin, line);
    std::stringstream word(line);
    word >> num_edges_ >> num_nodes_;

    for (int i = 0; i < num_nodes_; ++i)
    {
       Hypernode node(i, 1.0);
       hnodes_.push_back(std::move(node)); 
    }

    int id = 0;
    while (std::getline(fin, line)) 
    {
        if (line.empty())
        {
            continue;
        }

        std::stringstream word(line);
        Hyperedge hedge(id, 1.0);
        int buff = 0;
        while (word >> buff) {
            hedge.nodes.push_back(buff - 1);
            hnodes_[buff-1].edges.push_back(id);
        }
        hedges_.push_back(std::move(hedge));
        id++;
    }



}