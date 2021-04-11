// Node.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "Node.h"
#include "string.h"

namespace CASENA
{
	Node::Node(){Initialize();}
	Node::Node(const Node& node){*this = node;}
	Node::~Node(){Reset();}
	Node& Node::operator=(const Node& node)
	{
		id = node.id;
		memcpy(voltages,node.voltages,HistoryCount*sizeof(double));
		return *this;
	}
	void Node::Reset(){Initialize();}
	void Node::ID(const unsigned int& value){id = value;}
	unsigned int Node::ID() const{return id;}
	double Node::Voltage() const{return voltages[0];}
	double Node::PreviousVoltage() const{return voltages[1];}
	void Node::Voltage(const double& value)
	{
		for(unsigned int i = HistoryCount - 1 ; i > 0 ; i--)
		{
			voltages[i] = voltages[i - 1];
		}
		voltages[0] = value;
	}
	void Node::Initialize()
	{
		id = 0;
		memset(voltages,0,HistoryCount*sizeof(double));
	}
}



