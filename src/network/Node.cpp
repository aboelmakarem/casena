// Node.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "Node.h"
#include "Branch.h"

namespace CASENA
{
	Node::Node(){Initialize();}
	Node::Node(const Node& node) : Component(node){*this = node;}
	Node::~Node(){Reset();}
	Node& Node::operator=(const Node& node)
	{
		Component::operator=(node);
		for(unsigned int i = 0 ; i < HistoryCount ; i++)
		{
			voltages[i] = node.voltages[i];
		}
		branches = node.branches;
		return *this;
	}
	void Node::Reset()
	{
		Initialize();
		Component::Reset();
	}
	double Node::Voltage() const{return voltages[0];}
	double Node::PreviousVoltage() const{return voltages[1];}
	void Node::Voltage(const double& value)
	{
		for(unsigned int i = HistoryCount ; i > 0 ; i--)
		{
			voltages[i] = voltages[i - 1];
		}
		voltages[0] = value;
	}
	void Node::AddBranch(const Branch* branch)
	{
		if(branch == 0)		return;
		branches.PushBack(branch);
	}
	unsigned int Node::BranchCount() const{return branches.Size();}
	const EZ::List<const Branch*>* Node::Branches() const{return &branches;}
	void Node::Equation(EZ::Math::Matrix& f) const
	{
		// all node and branch IDs are zero based
		// no checks on topology consistency are made here, the network 
		// is assumed to be consistent (no loops, node belongs to its 
		// branches and vice versa, ...)
		// sum currents entering into (-) and leaving (+) this node
		const Branch* branch = 0;
		double current_sum = 0.0;
		double branch_direction = 0.0;
		for(EZ::ListItem<const Branch*>* item = branches.Start() ; item != 0 ; item = item->Next())
		{
			branch = item->Data();
			if(branch->StartNode() == this)		branch_direction = 1.0;
			else if(branch->EndNode() == this)	branch_direction = -1.0;
			else 								continue;
			current_sum += (branch_direction*branch->Current());
		}
		f(ID(),0,current_sum);
	}
	void Node::Gradients(EZ::Math::Matrix& A) const
	{
		// all node and branch IDs are zero based
		// no checks on topology consistency are made here, the network 
		// is assumed to be consistent (no loops, node belongs to its 
		// branches and vice versa, ...)
		const Branch* branch = 0;
		unsigned int branch_id = 0;
		unsigned int node_id = ID();
		for(EZ::ListItem<const Branch*>* item = branches.Start() ; item != 0 ; item = item->Next())
		{
			branch = item->Data();
			branch_id = branch->ID();
			if(branch->StartNode() == this)		A(node_id,branch_id,1.0);
			else if(branch->EndNode() == this)	A(node_id,branch_id,-1.0);
		}
	}
	void Node::Initialize()
	{
		for(unsigned int i = 0 ; i < HistoryCount ; i++)
		{
			voltages[i] = 0.0;
		}
		branches.Reset();
	}
}



