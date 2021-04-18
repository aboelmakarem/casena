// Component.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/14/2021

#ifndef COMPONENT_H_
#define COMPONENT_H_

#define HistoryCount 3

#include "Matrix.h"
#include "String.h"

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
		virtual bool Read(const char* line);
		void ID(const unsigned int& value);
		unsigned int ID() const;
		void Name(const EZ::String& target_name);
		const EZ::String& Name() const;
		virtual unsigned int ClaimIDs(const unsigned int& start_id);
		virtual void Equation(EZ::Math::Matrix& f) const = 0;
		virtual void Gradients(EZ::Math::Matrix& A) const = 0;
		virtual void Update(const EZ::Math::Matrix& x,const unsigned int& id_offset) = 0;
		virtual void Print() const = 0;
		
	private:
		void Initialize();
		unsigned int id;
		EZ::String name;
	};
}

#endif
