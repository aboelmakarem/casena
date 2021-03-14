// Node.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/14/2021

#include "Component.h"

namespace CASENA
{
	Component::Component(){Initialize();}
	Component::Component(const Component& component){*this = component;}
	Component::~Component(){Reset();}
	Component& Component::operator=(const Component& component)
	{
		id = component.id;
		return *this;
	}
	void Component::Reset(){Initialize();}
	void Component::ID(const unsigned int& value){id = value;}
	unsigned int Component::ID() const{return id;}
	void Component::Initialize(){id = 0;}
}

