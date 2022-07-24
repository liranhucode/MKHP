#pragma once

#ifndef FM_OR_H

#include <vector>
#include <queue>

#include "hgraph.h"
#include "mystruct.h"
#include "rpq.h"


#include "hgraph.h"

class Node
{
	
	friend class fm_or;

public:
	// Constructor and destructor
	Node(const int& id) :
		_id(id), _prev(NULL), _next(NULL) { }
	~Node() { }

	// Basic access methods
	int getId() const { return _id; }
	Node* getPrev() const { return _prev; }
	Node* getNext() const { return _next; }

	// Set functions
	void setId(const int& id) { _id = id; }
	void setPrev(Node* prev) { _prev = prev; }
	void setNext(Node* next) { _next = next; }

private:
	int         _id;    // id of the node (indicating the cell)
	Node*       _prev;  // pointer to the previous node
	Node*       _next;  // pointer to the next node
};


class fm_or {


public:

	void FMAlgorithm(s_hgraph * hgraph, vector<float> tpwgts2);

	void insertCell(map<int, Node*>  bList[], s_hgraph *hgraph, int i);

	void buildBList(s_hgraph *hgraph, int length, int vertex_degree);


private:

	map<int, Node*>  bList[2];

};

#endif