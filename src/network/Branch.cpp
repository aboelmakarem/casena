// Branch.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "Branch.h"

namespace CASENA
{
	Branch::Branch(){Initialize();}
	Branch::Branch(const Branch& branch){*this = branch;}
	Branch::~Branch(){Reset();}
	Branch& Branch::operator=(const Branch& branch)
	{
		id = branch.id;
		for(unsigned int i = 0 ; i < HistoryCount ; i++)
		{
			currents[i] = branch.currents[i];
		}
		start = branch.start;
		end = branch.end;
		return *this;
	}
	void Branch::Reset(){Initialize();}
	void Branch::ID(const unsigned int& value){id = value;}
	unsigned int Branch::ID() const{return id;}
	double Branch::Current() const{return currents[0];}
	double Branch::PreviousCurrent() const{return currents[1];}
	void Branch::Current(const double& value)
	{
		for(unsigned int i = HistoryCount ; i > 0 ; i--)
		{
			currents[i] = currents[i - 1];
		}
		currents[0] = value;
	}
	void Branch::StartNode(const Node* node){start = node;}
	const Node* Branch::StartNode() const{return start;}
	void Branch::EndNode(const Node* node){end = node;}
	const Node* Branch::EndNode() const{return end;}
	void Branch::Initialize()
	{
		id = 0;
		for(unsigned int i = 0 ; i < HistoryCount ; i++)
		{
			currents[i] = 0.0;
		}
		start = 0;
		end = 0;
	}
}



