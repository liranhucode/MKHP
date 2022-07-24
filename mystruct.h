#pragma once
#ifndef MYSTRUCT_H_
#define MYSTRUCT_H_


/* Types of vertex statues in the priority queue */
#define RPQ_NOTPRESENT   1       /* The vertex is in the queue */
#define RPQ_EXTRACTED    2       /* The vertex has been extracted from the queue */
#define RPQ_PRESENT   3       /* The vertex is not present in the queue and has not been extracted before */

#define NOT_MOVED -1


//���ڴֻ��׶�
struct hyperedge
{
	int number;					/*���泬�����*/
	int weight;					/*���泬��Ȩֵ*/
	int degree;					/*���泬�ߵĶ�*/
	//int total_vertex_weigth;	/*���泬���ж����Ȩֵ֮��*/
	int priority;
};

//����fm�㷨
struct gain_node
{
	int number;		//�������
	float gain;		//����



	gain_node(int a = 0, float b = 0) :
		number(a), gain(b) {}

	//gain�������
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