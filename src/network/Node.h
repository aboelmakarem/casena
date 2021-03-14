// Node.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef NODE_H_
#define NODE_H_

#define HistoryCount 3

#include "List.h"

namespace CASENA
{
	class Branch;
	class Node
	{
	public:
		Node();
		~Node();
		void Reset();
		void ID(const unsigned int& value);
		unsigned int ID() const;
		double Voltage() const;
		double PreviousVoltage() const;
		void Voltage(const double& value);
		void AddBranch(const Branch* branch);
		unsigned int BranchCount() const;
		const EZ::List<const Branch*>* Branches() const;
		double Equation(const double* variables) const;
		void Gradients(const double* variables,double* gradients) const;
		
	private:
		Node(const Node& node);
		Node& operator=(const Node& node);
		void Initialize();
		unsigned int id;
		double voltages[HistoryCount];
		EZ::List<const Branch*> branches;
	};
}

#endif

