// Node.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef NODE_H_
#define NODE_H_

#include "Component.h"

namespace CASENA
{
	class Node
	{
	public:
		Node();
		Node(const Node& node);
		~Node();
		Node& operator=(const Node& node);
		void Reset();
		void ID(const unsigned int& value);
		unsigned int ID() const;
		double Voltage() const;
		double PreviousVoltage() const;
		void Voltage(const double& value);
		
	private:
		void Initialize();
		unsigned int id;
		double voltages[HistoryCount];
	};
}

#endif

