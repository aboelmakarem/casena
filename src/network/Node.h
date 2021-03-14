// Node.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef NODE_H_
#define NODE_H_

#include "Component.h"
#include "List.h"

namespace CASENA
{
	class Branch;
	class Node : public Component
	{
	public:
		Node();
		~Node();
		void Reset();
		double Voltage() const;
		double PreviousVoltage() const;
		void Voltage(const double& value);
		void AddBranch(const Branch* branch);
		unsigned int BranchCount() const;
		const EZ::List<const Branch*>* Branches() const;
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		Node(const Node& node);
		Node& operator=(const Node& node);
		void Initialize();
		double voltages[HistoryCount];
		EZ::List<const Branch*> branches;
	};
}

#endif

