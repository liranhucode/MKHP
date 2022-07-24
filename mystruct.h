#pragma once
#ifndef MYSTRUCT_H_
#define MYSTRUCT_H_


/* Types of vertex statues in the priority queue */
#define RPQ_NOTPRESENT   1       /* The vertex is in the queue */
#define RPQ_EXTRACTED    2       /* The vertex has been extracted from the queue */
#define RPQ_PRESENT   3       /* The vertex is not present in the queue and has not been extracted before */

#define NOT_MOVED -1


//用于粗化阶段
struct hyperedge
{
	int number;					/*保存超边序号*/
	int weight;					/*保存超边权值*/
	int degree;					/*保存超边的度*/
	//int total_vertex_weigth;	/*保存超边中顶点的权值之和*/
	int priority;
};

//用于fm算法
struct gain_node
{
	int number;		//顶点序号
	float gain;		//增益



	gain_node(int a = 0, float b = 0) :
		number(a), gain(b) {}

	//gain大的优先
	friend bool operator<(gain_node a, gain_node b)
	{
		if (a.gain == b.gain)
			return a.number > b.number;

		return a.gain < b.gain;
	}

};


struct mycmp {
	bool operator()(gain_node a, gain_node b) {
		return a.gain > b.gain;
	}
};



#endif