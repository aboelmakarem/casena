// Node.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef NODE_H_
#define NODE_H_

#include "Component.h"
#include "List.h"

namespace CASENA
{
	class Node : public Component
	{
	public:
		Node();
		Node(const Node& node);
		~Node();
		Node& operator=(const Node& node);
		void Reset();
		double Voltage() const;
		double PreviousVoltage() const;
		void Voltage(const double& value);
		void AddBranch(const unsigned int& branch_id);
		unsigned int BranchCount() const;
		const EZ::List<unsigned int>* BranchIDs() const;
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
		double voltages[HistoryCount];
		EZ::List<unsigned int> branch_ids;
	};
}

#endif

