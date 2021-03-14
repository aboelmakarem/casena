// Branch.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef BRANCH_H_
#define BRANCH_H_

#include "Node.h"

namespace CASENA
{
	class Branch
	{
	public:
		Branch();
		Branch(const Branch& branch);
		~Branch();
		Branch& operator=(const Branch& branch);
		void Reset();
		void ID(const unsigned int& value);
		unsigned int ID() const;
		double Current() const;
		double PreviousCurrent() const;
		void Current(const double& value);
		void StartNode(const Node* node);
		const Node* StartNode() const;
		void EndNode(const Node* node);
		const Node* EndNode() const;
		virtual double Equation(const double* variables) const = 0;
		virtual void Gradients(const double* variables,double* gradients) const = 0;
		
	private:
		void Initialize();
		unsigned int id;
		double currents[HistoryCount];
		const Node* start;
		const Node* end;
	};
}

#endif

