// Component.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/14/2021

#ifndef COMPONENT_H_
#define COMPONENT_H_

#define HistoryCount 3

#include "Matrix.h"

namespace CASENA
{
	class Component
	{
	public:
		Component();
		Component(const Component& component);
		virtual ~Component();
		virtual Component& operator=(const Component& component);
		void Reset();
		void ID(const unsigned int& value);
		unsigned int ID() const;
		virtual void Equation(EZ::Math::Matrix& f) const = 0;
		virtual void Gradients(EZ::Math::Matrix& A) const = 0;
		
	private:
		void Initialize();
		unsigned int id;
	};
}

#endif